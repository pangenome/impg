//! 1aln format parsing
//!
//! Handles parsing of ONE alignment format (.1aln files) using the onecode library.

use crate::alignment_record::{AlignmentRecord, Strand};
use crate::seqidx::SequenceIndex;
use onecode::OneFile;
use std::cell::RefCell;
use std::collections::HashMap;
use std::convert::TryFrom;
use log::{debug, warn};

thread_local! {
    static ONE_ALN_FILE_CACHE: RefCell<HashMap<String, OneFile>> = RefCell::new(HashMap::new());
}

/// 1aln file parser with metadata and O(1) seeking support
pub struct OneAlnParser {
    file_path: String,
    trace_spacing: i64,
    metadata: OneAlnMetadata,
}

struct OneAlnMetadata {
    // Query genome (gdb1 - first reference, or embedded if self-alignment)
    query_seq_names: HashMap<i64, String>,
    query_seq_lengths: HashMap<i64, i64>,
    query_contig_offsets: HashMap<i64, (i64, i64)>,

    // Target genome (gdb2 - second reference, or embedded skeleton)
    target_seq_names: HashMap<i64, String>,
    target_seq_lengths: HashMap<i64, i64>,
    target_contig_offsets: HashMap<i64, (i64, i64)>,
}

/// Error type for 1aln parsing
#[derive(Debug)]
pub enum ParseErr {
    InvalidFormat(String),
}

impl std::fmt::Display for ParseErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseErr::InvalidFormat(msg) => write!(f, "{}", msg),
        }
    }
}

impl std::error::Error for ParseErr {}

impl OneAlnParser {
    /// Open a 1aln file and read its metadata
    ///
    /// # Arguments
    /// * `file_path` - Path to the .1aln file
    /// * `sequence_file_hints` - Optional list of sequence file paths to help locate GDB files in their directories
    pub fn new(file_path: String, sequence_file_hints: Option<&[String]>) -> Result<Self, ParseErr> {
        let mut file = OneFile::open_read(&file_path, None, None, 1)
            .map_err(|e| ParseErr::InvalidFormat(format!("Failed to open 1aln file: {}", e)))?;

        // Check if there are reference paths to external GDB files
        let references = file.get_references();

        // Initialize empty metadata structures
        let mut query_seq_names = HashMap::new();
        let mut query_seq_lengths = HashMap::new();
        let mut query_contig_offsets = HashMap::new();
        let mut target_seq_names = HashMap::new();
        let mut target_seq_lengths = HashMap::new();
        let mut target_contig_offsets = HashMap::new();

        // The embedded GDB skeleton (if present) is the target genome (gdb2)
        let embedded_names = file.get_all_sequence_names();
        let embedded_lengths = file.get_all_sequence_lengths();
        let embedded_offsets = file.get_all_contig_offsets();

        // Process references:
        // - First reference (count=1) is the QUERY genome (gdb1/A-read)
        // - Second reference (count=2) is the TARGET genome (gdb2/B-read)
        // - If there's an embedded skeleton, it's for the TARGET (gdb2/B-read)

        let mut has_external_query = false;
        let mut has_external_target = false;

        for (ref_idx, (ref_path, ref_count)) in references.iter().enumerate() {
            if ref_path.is_empty() {
                continue;
            }

            // Skip if this is not a genome reference (count > 2 might be other metadata)
            if *ref_count > 2 {
                continue;
            }

            let is_query = *ref_count == 1;    // First reference is query (A-read)
            let is_target = *ref_count == 2;   // Second reference is target (B-read)

            debug!(
                "Processing reference {}: {} (count: {}, type: {})",
                ref_idx + 1,
                ref_path,
                ref_count,
                if is_query { "query" } else if is_target { "target" } else { "unknown" }
            );

            // Helper function to strip fasta and AGC extensions
            let strip_seq_ext = |p: &str| -> String {
                let path_str = p.to_string();
                // Try to remove common sequence file extensions (FASTA and AGC)
                for ext in &[".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna", ".agc"] {
                    if path_str.ends_with(ext) {
                        return path_str[..path_str.len() - ext.len()].to_string();
                    }
                }
                path_str
            };

            // Get alignment file directory for relative path resolution
            let aln_dir = std::path::Path::new(&file_path)
                .parent()
                .unwrap_or_else(|| std::path::Path::new("."));

            // Try multiple strategies to find the GDB file
            let mut gdb_path: Option<String> = None;

            // Strategy 1: Check in directories of sequence files provided via command-line hints
            if gdb_path.is_none() {
                if let Some(seq_files) = sequence_file_hints {
                    // Extract the base name from ref_path (just the filename without directory)
                    let ref_base_name = std::path::Path::new(ref_path)
                        .file_name()
                        .and_then(|n| n.to_str())
                        .unwrap_or("");
                    let ref_base_stripped = strip_seq_ext(ref_base_name);

                    // Check each sequence file directory
                    for seq_file in seq_files {
                        let seq_path = std::path::Path::new(seq_file);
                        if let Some(seq_dir) = seq_path.parent() {
                            // Check if this sequence file matches the reference
                            let seq_base_name = seq_path.file_name().and_then(|n| n.to_str()).unwrap_or("");
                            let seq_base_stripped = strip_seq_ext(seq_base_name);

                            // If basenames match (with or without extensions), look for GDB here
                            if !ref_base_stripped.is_empty() && (
                                seq_base_stripped == ref_base_stripped ||
                                seq_base_name == ref_base_name ||
                                seq_base_name == ref_base_stripped
                            ) {
                                for ext in &[".1gdb", ".gdb"] {
                                    let with_ext = seq_dir.join(format!("{}{}", seq_base_stripped, ext));
                                    if with_ext.exists() {
                                        gdb_path = Some(with_ext.to_string_lossy().to_string());
                                        debug!("Found GDB via sequence file hint: {}", gdb_path.as_ref().unwrap());
                                        break;
                                    }
                                }
                                if gdb_path.is_some() {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // Strategy 2: Try as absolute path (as-is) - but only if it's a GDB file
            if gdb_path.is_none() && std::path::Path::new(ref_path).exists() && (ref_path.ends_with(".1gdb") || ref_path.ends_with(".gdb")) {
                gdb_path = Some(ref_path.clone());
            }

            // Strategy 3: Try adding .1gdb or .gdb extension to the original path
            if gdb_path.is_none() {
                for ext in &[".1gdb", ".gdb"] {
                    let with_ext = format!("{}{}", ref_path, ext);
                    if std::path::Path::new(&with_ext).exists() {
                        gdb_path = Some(with_ext);
                        break;
                    }
                }
            }

            // Strategy 4: Strip sequence file extension and try with GDB extensions (absolute path)
            if gdb_path.is_none() {
                let base_path = strip_seq_ext(ref_path);
                if base_path != *ref_path {
                    for ext in &[".1gdb", ".gdb"] {
                        let with_ext = format!("{}{}", base_path, ext);
                        if std::path::Path::new(&with_ext).exists() {
                            gdb_path = Some(with_ext);
                            break;
                        }
                    }
                }
            }

            // Strategy 5: Try in the same directory as the reference sequence file
            if gdb_path.is_none() {
                let ref_path_obj = std::path::Path::new(ref_path);
                if let Some(ref_dir) = ref_path_obj.parent() {
                    if let Some(file_name) = ref_path_obj.file_name().and_then(|n| n.to_str()) {
                        let base_name = strip_seq_ext(file_name);
                        for ext in &[".1gdb", ".gdb"] {
                            let with_ext = ref_dir.join(format!("{}{}", base_name, ext));
                            if with_ext.exists() {
                                gdb_path = Some(with_ext.to_string_lossy().to_string());
                                break;
                            }
                        }
                    }
                }
            }

            // Strategy 6: Try relative path with sequence file extension stripped and GDB extension added
            if gdb_path.is_none() {
                let base_path = strip_seq_ext(ref_path);
                if base_path != *ref_path {
                    for ext in &[".1gdb", ".gdb"] {
                        let relative_with_ext = aln_dir.join(format!("{}{}", base_path, ext));
                        if relative_with_ext.exists() {
                            gdb_path = Some(relative_with_ext.to_string_lossy().to_string());
                            break;
                        }
                    }
                }
            }

            // Strategy 7: Try relative path with .1gdb or .gdb extension
            if gdb_path.is_none() {
                for ext in &[".1gdb", ".gdb"] {
                    let relative_with_ext = aln_dir.join(format!("{}{}", ref_path, ext));
                    if relative_with_ext.exists() {
                        gdb_path = Some(relative_with_ext.to_string_lossy().to_string());
                        break;
                    }
                }
            }

            // Strategy 8: Try relative to alignment file directory (as-is) - but only if it's a GDB file
            if gdb_path.is_none() {
                let relative_path = aln_dir.join(ref_path);
                if relative_path.exists() && (ref_path.ends_with(".1gdb") || ref_path.ends_with(".gdb")) {
                    gdb_path = Some(relative_path.to_string_lossy().to_string());
                }
            }

            let gdb_path = if let Some(found_path) = gdb_path {
                found_path
            } else {
                return Err(ParseErr::InvalidFormat(format!(
                    "GDB file not found for '{}'. Run: FAtoGDB {}",
                    ref_path, ref_path
                )));
            };

            // Try to load the GDB metadata
            match OneFile::read_gdb_metadata(&gdb_path) {
                Ok((ref_names, ref_lengths, ref_offsets)) => {
                    if is_query {
                        query_seq_names = ref_names;
                        query_seq_lengths = ref_lengths;
                        query_contig_offsets = ref_offsets;
                        has_external_query = true;
                        debug!("Loaded query genome metadata from: {} ({} sequences)", gdb_path, query_seq_names.len());
                    } else if is_target {
                        target_seq_names = ref_names;
                        target_seq_lengths = ref_lengths;
                        target_contig_offsets = ref_offsets;
                        has_external_target = true;
                        debug!("Loaded target genome metadata from: {} ({} sequences)", gdb_path, target_seq_names.len());
                    }
                }
                Err(e) => {
                    return Err(ParseErr::InvalidFormat(format!(
                        "Failed to read GDB file '{}': {}. Run: FAtoGDB {}",
                        gdb_path, e, ref_path
                    )));
                }
            }
        }

        // If we didn't load external target, use embedded skeleton
        if !has_external_target && !embedded_names.is_empty() {
            target_seq_names = embedded_names;
            target_seq_lengths = embedded_lengths;
            target_contig_offsets = embedded_offsets;
            debug!("Using embedded skeleton for target genome ({} sequences)", target_seq_names.len());
        }

        // If this is a self-alignment (no external query), use target for query too
        if !has_external_query && !target_seq_names.is_empty() {
            query_seq_names = target_seq_names.clone();
            query_seq_lengths = target_seq_lengths.clone();
            query_contig_offsets = target_contig_offsets.clone();
            warn!("Self-alignment detected: using target genome for query");
        }

        if query_seq_names.is_empty() && target_seq_names.is_empty() {
            warn!("Warning: No sequence metadata found in file or external references");
        }

        let metadata = OneAlnMetadata {
            query_seq_names,
            query_seq_lengths,
            query_contig_offsets,
            target_seq_names,
            target_seq_lengths,
            target_contig_offsets,
        };

        // Get trace spacing
        let mut trace_spacing = 100; // default
        loop {
            match file.read_line() {
                't' => {
                    trace_spacing = file.int(0);
                    break;
                }
                'A' | '\0' => break,
                _ => {}
            }
        }

        Ok(OneAlnParser {
            file_path,
            trace_spacing,
            metadata,
        })
    }

    /// Parse all alignments from the 1aln file into AlignmentRecords
    pub fn parse_alignments(
        &self,
        seq_index: &mut SequenceIndex,
    ) -> Result<Vec<AlignmentRecord>, ParseErr> {
        let (query_contig_to_seq_id, target_contig_to_seq_id) = self.prepare_sequence_index(seq_index)?;

        let mut file = OneFile::open_read(&self.file_path, None, None, 1)
            .map_err(|e| ParseErr::InvalidFormat(format!("Failed to open 1aln file: {}", e)))?;
        let mut records = Vec::new();
        let mut alignment_index: u64 = 0;
        let mut current_line = file.read_line();

        loop {
            match current_line {
                '\0' => break,
                'A' => {
                    let (record, next_line) =
                        self.parse_single_alignment(&mut file, &query_contig_to_seq_id, &target_contig_to_seq_id, alignment_index)?;
                    records.push(record);
                    alignment_index += 1;
                    current_line = next_line;
                }
                _ => {
                    current_line = file.read_line();
                }
            };
        }

        Ok(records)
    }

    /// Parse a single alignment from current file position (after reading 'A' line)
    /// Returns (AlignmentRecord, found_next_A)
    fn parse_single_alignment(
        &self,
        file: &mut OneFile,
        query_contig_to_seq_id: &HashMap<i64, u32>,
        target_contig_to_seq_id: &HashMap<i64, u32>,
        alignment_index: u64,
    ) -> Result<(AlignmentRecord, char), ParseErr> {
        // Read alignment coordinates from current 'A' line
        let query_contig_id = file.int(0);
        let target_contig_id = file.int(3);

        // Get sequence names and lengths from metadata
        let query_name = self
            .metadata
            .query_seq_names
            .get(&query_contig_id)
            .cloned()
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Query sequence with contig ID {} not found in metadata",
                    query_contig_id
                ))
            })?;
        let target_name = self
            .metadata
            .target_seq_names
            .get(&target_contig_id)
            .cloned()
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Target sequence with contig ID {} not found in metadata",
                    target_contig_id
                ))
            })?;

        // Lookup already-registered scaffold IDs from separate namespaces
        let query_id = query_contig_to_seq_id
            .get(&query_contig_id)
            .copied()
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Sequence index ID for query {query_name} (contig {query_contig_id}) missing"
                ))
            })?;
        let target_id = target_contig_to_seq_id
            .get(&target_contig_id)
            .copied()
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Sequence index ID for target {target_name} (contig {target_contig_id}) missing"
                ))
            })?;

        let query_contig_start = file.int(1) as usize;
        let query_contig_end = file.int(2) as usize;
        let mut target_contig_start = file.int(4) as usize;
        let mut target_contig_end = file.int(5) as usize;

        let (query_contig_offset, _query_contig_len) = self
            .metadata
            .query_contig_offsets
            .get(&query_contig_id)
            .copied()
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Contig offset for query contig {} missing from metadata",
                    query_contig_id
                ))
            })?;

        let (target_contig_offset, target_contig_len) = self
            .metadata
            .target_contig_offsets
            .get(&target_contig_id)
            .copied()
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Contig offset for target contig {} missing from metadata",
                    target_contig_id
                ))
            })?;

        let mut strand = Strand::Forward;
        let mut num_tracepoints = 0;

        // Read associated lines
        let next_line = loop {
            let line_type = file.read_line();
            match line_type {
                'R' => strand = Strand::Reverse,
                'D' => {
                    // Differences (ignore for now)
                }
                'T' => {
                    // Tracepoints
                    if let Some(tp_vec) = file.int_list() {
                        num_tracepoints = tp_vec.len();
                    }
                }
                'X' => {
                    // Trace diffs (ignore for now)
                }
                'A' | 'a' | 'g' | 'S' | '^' | '\0' => break line_type,
                _ => {
                    // Skip other line types (D, X, etc.)
                }
            }
        };

        if strand == Strand::Reverse {
            let contig_len = target_contig_len as usize;
            let orig_start = target_contig_start;
            let orig_end = target_contig_end;
            target_contig_start = contig_len - orig_end;
            target_contig_end = contig_len - orig_start;
        }

        let query_scaffold_start = usize::try_from(query_contig_offset + query_contig_start as i64)
            .map_err(|_| {
                ParseErr::InvalidFormat(format!(
                    "Query scaffold start overflow for {query_name}: offset {query_contig_offset}, start {query_contig_start}"
                ))
            })?;
        let query_scaffold_end = usize::try_from(query_contig_offset + query_contig_end as i64)
            .map_err(|_| {
                ParseErr::InvalidFormat(format!(
                    "Query scaffold end overflow for {query_name}: offset {query_contig_offset}, end {query_contig_end}"
                ))
            })?;
        let target_scaffold_start =
            usize::try_from(target_contig_offset + target_contig_start as i64).map_err(|_| {
                ParseErr::InvalidFormat(format!(
                    "Target scaffold start overflow for {target_name}: offset {target_contig_offset}, start {target_contig_start}"
                ))
            })?;
        let target_scaffold_end =
            usize::try_from(target_contig_offset + target_contig_end as i64).map_err(|_| {
                ParseErr::InvalidFormat(format!(
                    "Target scaffold end overflow for {target_name}: offset {target_contig_offset}, end {target_contig_end}"
                ))
            })?;

        let mut record = AlignmentRecord {
            query_id, // Scaffold ID from SequenceIndex
            query_start: query_scaffold_start,
            query_end: query_scaffold_end,
            target_id, // Scaffold ID from SequenceIndex
            target_start: target_scaffold_start,
            target_end: target_scaffold_end,
            strand_and_data_offset: alignment_index, // Store alignment index for O(1) seeking
            data_bytes: num_tracepoints,             // Store number of tracepoints
        };
        record.set_strand(strand);

        Ok((record, next_line))
    }

    fn prepare_sequence_index(
        &self,
        seq_index: &mut SequenceIndex,
    ) -> Result<(HashMap<i64, u32>, HashMap<i64, u32>), ParseErr> {
        let mut query_contig_to_seq_id = HashMap::with_capacity(self.metadata.query_seq_names.len());
        let mut target_contig_to_seq_id = HashMap::with_capacity(self.metadata.target_seq_names.len());

        // Register query sequences
        for (&contig_id, name) in &self.metadata.query_seq_names {
            let length = self
                .metadata
                .query_seq_lengths
                .get(&contig_id)
                .copied()
                .ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Query sequence length for contig {contig_id} ({name}) missing from metadata"
                    ))
                })? as usize;
            let seq_id = seq_index.get_or_insert_id(name, Some(length));
            query_contig_to_seq_id.insert(contig_id, seq_id);
        }

        // Register target sequences
        for (&contig_id, name) in &self.metadata.target_seq_names {
            let length = self
                .metadata
                .target_seq_lengths
                .get(&contig_id)
                .copied()
                .ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Target sequence length for contig {contig_id} ({name}) missing from metadata"
                    ))
                })? as usize;
            let seq_id = seq_index.get_or_insert_id(name, Some(length));
            target_contig_to_seq_id.insert(contig_id, seq_id);
        }

        Ok((query_contig_to_seq_id, target_contig_to_seq_id))
    }

    /// Get trace spacing from the .1aln file
    pub fn get_trace_spacing(&self) -> u32 {
        self.trace_spacing as u32
    }

    /// Get query sequence names mapping (contig_id -> name)
    pub fn get_query_sequence_names(&self) -> &HashMap<i64, String> {
        &self.metadata.query_seq_names
    }

    /// Get target sequence names mapping (contig_id -> name)
    pub fn get_target_sequence_names(&self) -> &HashMap<i64, String> {
        &self.metadata.target_seq_names
    }

    /// Get query contig offset information (sbeg, clen) for all contigs
    pub fn get_query_contig_offsets(&self) -> &HashMap<i64, (i64, i64)> {
        &self.metadata.query_contig_offsets
    }

    /// Get target contig offset information (sbeg, clen) for all contigs
    pub fn get_target_contig_offsets(&self) -> &HashMap<i64, (i64, i64)> {
        &self.metadata.target_contig_offsets
    }

    /// Seek to a specific alignment using O(1) access
    pub fn seek_alignment(&self, alignment_index: u64) -> Result<OneAlnAlignment, ParseErr> {
        self.with_cached_file(|file| {
            // Use O(1) binary index to jump to alignment
            file.goto('A', (alignment_index + 1) as i64).map_err(|e| {
                ParseErr::InvalidFormat(format!(
                    "Failed to seek to alignment {}: {}.",
                    alignment_index, e
                ))
            })?;

            file.read_line(); // Read the 'A' line we jumped to

            // Read alignment coordinates
            let query_contig_id = file.int(0);
            let target_contig_id = file.int(3);

            let query_name = self
                .metadata
                .query_seq_names
                .get(&query_contig_id)
                .cloned()
                .ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Query sequence with contig ID {} not found in metadata",
                        query_contig_id
                    ))
                })?;
            let query_length = self
                .metadata
                .query_seq_lengths
                .get(&query_contig_id)
                .copied()
                .ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Query sequence length for contig ID {} not found in metadata",
                        query_contig_id
                    ))
                })?;

            let target_name = self
                .metadata
                .target_seq_names
                .get(&target_contig_id)
                .cloned()
                .ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Target sequence with contig ID {} not found in metadata",
                        target_contig_id
                    ))
                })?;
            let target_length = self
                .metadata
                .target_seq_lengths
                .get(&target_contig_id)
                .copied()
                .ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Target sequence length for contig ID {} not found in metadata",
                        target_contig_id
                    ))
                })?;

            let (query_contig_offset, _) = self
                .metadata
                .query_contig_offsets
                .get(&query_contig_id)
                .copied()
                .ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Contig offset for query contig {} missing from metadata",
                        query_contig_id
                    ))
                })?;

            let (target_contig_offset, target_contig_len) = self
                .metadata
                .target_contig_offsets
                .get(&target_contig_id)
                .copied()
                .ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Contig offset for target contig {} missing from metadata",
                        target_contig_id
                    ))
                })?;

            let alignment = OneAlnAlignment {
                query_name,
                query_length,
                query_contig_start: file.int(1),
                query_contig_end: file.int(2),
                query_contig_offset,
                target_name,
                target_length,
                target_contig_start: file.int(4),
                target_contig_end: file.int(5),
                target_contig_offset,
                target_contig_len,
                strand: '+',
                differences: 0,
                tracepoints: Vec::new(),
                trace_diffs: Vec::new(),
                trace_spacing: self.trace_spacing,
            };

            // Read associated lines to get tracepoints
            self.read_alignment_details(file, alignment)
        })
    }

    fn with_cached_file<R, F>(&self, f: F) -> Result<R, ParseErr>
    where
        F: FnOnce(&mut OneFile) -> Result<R, ParseErr>,
    {
        ONE_ALN_FILE_CACHE.with(|cache| -> Result<R, ParseErr> {
            let mut cache = cache.borrow_mut();
            if !cache.contains_key(&self.file_path) {
                let file = OneFile::open_read(&self.file_path, None, None, 1).map_err(|e| {
                    ParseErr::InvalidFormat(format!("Failed to open 1aln file: {}", e))
                })?;
                cache.insert(self.file_path.clone(), file);
            }
            let file = cache.get_mut(&self.file_path).expect("entry inserted above");
            f(file)
        })
    }

    fn read_alignment_details(
        &self,
        file: &mut OneFile,
        mut alignment: OneAlnAlignment,
    ) -> Result<OneAlnAlignment, ParseErr> {
        loop {
            match file.read_line() {
                'R' => alignment.strand = '-',
                'D' => alignment.differences = file.int(0),
                'T' => {
                    alignment.tracepoints = file.int_list().map(|v| v.to_vec()).unwrap_or_default()
                }
                'X' => {
                    alignment.trace_diffs = file.int_list().map(|v| v.to_vec()).unwrap_or_default()
                }
                'A' | 'a' | 'g' | '^' | '\0' => break,
                _ => {}
            }
        }

        // PAFtoALN can produce zero-length self-alignments with:
        // - query_start == query_end, target_start == target_end, 
        // - no differences, and one tracepoint.
        // This causes bugs in ALNtoPAF and needs fixing.
        if alignment.query_contig_start == alignment.query_contig_end && alignment.target_contig_start == alignment.target_contig_end
            && alignment.differences == 0
            && alignment.tracepoints.len() == 1 && alignment.trace_diffs.len() == 1
            // Check it's the same sequence
            && alignment.query_contig_end == alignment.target_contig_end && alignment.query_name == alignment.target_name
        {
            warn!(
                "Zero-length self-alignment at {}:{} (PAFtoALN artifact) - fixing coordinates",
                alignment.query_name,
                alignment.query_contig_start
            );

            alignment.query_contig_start = 0;
            alignment.target_contig_start = 0;
        }

        if alignment.differences == 0 {
            let query_len = alignment.query_contig_end - alignment.query_contig_start;
            let target_len = alignment.target_contig_end - alignment.target_contig_start;
            if query_len != target_len {
                return Err(ParseErr::InvalidFormat(format!(
                    "Zero-difference alignment but mismatched lengths: query {} (len {}), target {} (len {})",
                    alignment.query_name, query_len, alignment.target_name, target_len
                )));
            }
        }
        
        Ok(alignment)
    }
}

/// Represents a complete 1aln alignment with tracepoints for CIGAR reconstruction
#[derive(Debug)]
pub struct OneAlnAlignment {
    pub query_name: String,
    pub query_length: i64,
    pub query_contig_start: i64,
    pub query_contig_end: i64,
    pub query_contig_offset: i64,
    pub target_name: String,
    pub target_length: i64,
    pub target_contig_start: i64,
    pub target_contig_end: i64,
    pub target_contig_offset: i64,
    pub target_contig_len: i64,
    pub strand: char,
    pub differences: i64,
    pub tracepoints: Vec<i64>,
    pub trace_diffs: Vec<i64>,
    pub trace_spacing: i64,
}

impl OneAlnAlignment {
    pub fn query_scaffold_span(&self) -> Result<(i32, i32), ParseErr> {
        let start = self
            .query_contig_offset
            .checked_add(self.query_contig_start)
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Overflow computing query scaffold start: offset {}, start {}",
                    self.query_contig_offset, self.query_contig_start
                ))
            })?;
        let end = self
            .query_contig_offset
            .checked_add(self.query_contig_end)
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Overflow computing query scaffold end: offset {}, end {}",
                    self.query_contig_offset, self.query_contig_end
                ))
            })?;
        let start_i32 = i32::try_from(start).map_err(|_| {
            ParseErr::InvalidFormat(format!(
                "Query scaffold start outside i32 range: {}",
                start
            ))
        })?;
        let end_i32 = i32::try_from(end).map_err(|_| {
            ParseErr::InvalidFormat(format!(
                "Query scaffold end outside i32 range: {}",
                end
            ))
        })?;
        Ok((start_i32, end_i32))
    }

    pub fn target_scaffold_span(&self) -> Result<(i32, i32), ParseErr> {
        // Adjust target contig-relative coordinates based on strand
        let (contig_start, contig_end) = match self.strand {
            '+' => (self.target_contig_start, self.target_contig_end),
            '-' => {
                let start = self
                    .target_contig_len
                    .checked_sub(self.target_contig_end)
                    .ok_or_else(|| {
                        ParseErr::InvalidFormat(format!(
                            "Invalid reverse-strand alignment: contig length {} smaller than end {}",
                            self.target_contig_len, self.target_contig_end
                        ))
                    })?;
                let end = self
                    .target_contig_len
                    .checked_sub(self.target_contig_start)
                    .ok_or_else(|| {
                        ParseErr::InvalidFormat(format!(
                            "Invalid reverse-strand alignment: contig length {} smaller than start {}",
                            self.target_contig_len, self.target_contig_start
                        ))
                    })?;
                (start, end)
            }
            other => {
                return Err(ParseErr::InvalidFormat(format!(
                    "Unexpected strand character '{other}' in OneAlnAlignment"
                )))
            }
        };

        let start = self
            .target_contig_offset
            .checked_add(contig_start)
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Overflow computing target scaffold start: offset {}, start {}",
                    self.target_contig_offset, contig_start
                ))
            })?;
        let end = self
            .target_contig_offset
            .checked_add(contig_end)
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Overflow computing target scaffold end: offset {}, end {}",
                    self.target_contig_offset, contig_end
                ))
            })?;
        let start_i32 = i32::try_from(start).map_err(|_| {
            ParseErr::InvalidFormat(format!(
                "Target scaffold start outside i32 range: {}",
                start
            ))
        })?;
        let end_i32 = i32::try_from(end).map_err(|_| {
            ParseErr::InvalidFormat(format!(
                "Target scaffold end outside i32 range: {}",
                end
            ))
        })?;
        Ok((start_i32, end_i32))
    }
}
