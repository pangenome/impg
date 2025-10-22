//! 1aln format parsing
//!
//! Handles parsing of ONE alignment format (.1aln files) using the onecode library.

use crate::alignment_record::{AlignmentRecord, Strand};
use crate::seqidx::SequenceIndex;
use onecode::OneFile;
use std::cell::RefCell;
use std::collections::HashMap;
use std::convert::TryFrom;

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
    seq_names: HashMap<i64, String>,
    seq_lengths: HashMap<i64, i64>,
    contig_offsets: HashMap<i64, (i64, i64)>, // Maps contig ID to (sbeg, clen)
}

/// Error type for 1aln parsing
#[derive(Debug)]
pub enum ParseErr {
    InvalidFormat(String),
}

impl OneAlnParser {
    /// Open a 1aln file and read its metadata
    pub fn new(file_path: String) -> Result<Self, ParseErr> {
        let mut file = OneFile::open_read(&file_path, None, None, 1)
            .map_err(|e| ParseErr::InvalidFormat(format!("Failed to open 1aln file: {}", e)))?;

        let metadata = OneAlnMetadata {
            seq_names: file.get_all_sequence_names(), // contig id (0-indexed) to scaffold name
            seq_lengths: file.get_all_sequence_lengths(), // contig id (0-indexed) to scaffold length
            contig_offsets: file.get_all_contig_offsets(), // contig id (0-indexed) to (scaffold_offset, contig_length)
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
        let contig_to_seq_id = self.prepare_sequence_index(seq_index)?;

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
                        self.parse_single_alignment(&mut file, &contig_to_seq_id, alignment_index)?;
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
        contig_to_seq_id: &HashMap<i64, u32>,
        alignment_index: u64,
    ) -> Result<(AlignmentRecord, char), ParseErr> {
        // Read alignment coordinates from current 'A' line
        let query_contig_id = file.int(0);
        let target_contig_id = file.int(3);

        // Get sequence names and lengths from metadata
        let query_name = self
            .metadata
            .seq_names
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
            .seq_names
            .get(&target_contig_id)
            .cloned()
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Target sequence with contig ID {} not found in metadata",
                    target_contig_id
                ))
            })?;

        // Lookup already-registered scaffold IDs
        let query_id = contig_to_seq_id
            .get(&query_contig_id)
            .copied()
            .ok_or_else(|| {
                ParseErr::InvalidFormat(format!(
                    "Sequence index ID for query {query_name} (contig {query_contig_id}) missing"
                ))
            })?;
        let target_id = contig_to_seq_id
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
            .contig_offsets
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
            .contig_offsets
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
                    "Query scaffold start overflow: offset {query_contig_offset}, start {query_contig_start}"
                ))
            })?;
        let query_scaffold_end = usize::try_from(query_contig_offset + query_contig_end as i64)
            .map_err(|_| {
                ParseErr::InvalidFormat(format!(
                    "Query scaffold end overflow: offset {query_contig_offset}, end {query_contig_end}"
                ))
            })?;
        let target_scaffold_start =
            usize::try_from(target_contig_offset + target_contig_start as i64).map_err(|_| {
                ParseErr::InvalidFormat(format!(
                    "Target scaffold start overflow: offset {target_contig_offset}, start {target_contig_start}"
                ))
            })?;
        let target_scaffold_end =
            usize::try_from(target_contig_offset + target_contig_end as i64).map_err(|_| {
                ParseErr::InvalidFormat(format!(
                    "Target scaffold end overflow: offset {target_contig_offset}, end {target_contig_end}"
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
    ) -> Result<HashMap<i64, u32>, ParseErr> {
        let mut contig_to_seq_id = HashMap::with_capacity(self.metadata.seq_names.len());
        for (&contig_id, name) in &self.metadata.seq_names {
            let length = self
                .metadata
                .seq_lengths
                .get(&contig_id)
                .copied()
                .ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Sequence length for contig {contig_id} ({name}) missing from metadata"
                    ))
                })? as usize;
            let seq_id = seq_index.get_or_insert_id(name, Some(length));
            contig_to_seq_id.insert(contig_id, seq_id);
        }
        Ok(contig_to_seq_id)
    }

    /// Get trace spacing from the .1aln file
    pub fn get_trace_spacing(&self) -> u32 {
        self.trace_spacing as u32
    }

    /// Get sequence names mapping (contig_id -> name)
    pub fn get_sequence_names(&self) -> &HashMap<i64, String> {
        &self.metadata.seq_names
    }

    /// Get contig offset information (sbeg, clen) for all contigs
    pub fn get_contig_offsets(&self) -> &HashMap<i64, (i64, i64)> {
        &self.metadata.contig_offsets
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
                .seq_names
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
                .seq_lengths
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
                .seq_names
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
                .seq_lengths
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
                .contig_offsets
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
                .contig_offsets
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

        if alignment.strand == '-' {
            let orig_start = alignment.target_contig_start;
            let orig_end = alignment.target_contig_end;
            alignment.target_contig_start = alignment.target_contig_len - orig_end;
            alignment.target_contig_end = alignment.target_contig_len - orig_start;
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
