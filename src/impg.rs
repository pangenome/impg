use crate::alignment_record::{AlignmentRecord, Strand};
use crate::forest_map::ForestMap;
use crate::graph::reverse_complement;
use crate::onealn::OneAlnParser;
use crate::paf::ParseErr;
use crate::seqidx::SequenceIndex;
use crate::sequence_index::SequenceIndex as _; // The as _ syntax imports the trait so its methods are available, but doesn't bring the name into scope (avoiding the naming conflict)
use crate::sequence_index::UnifiedSequenceIndex;
use coitrees::{BasicCOITree, Interval, IntervalTree};
use lib_tracepoints::{tracepoints_to_cigar_fastga_with_aligner, DistanceMode};
use lib_wfa2::affine_wavefront::AffineWavefronts;
use log::debug;
use noodles::bgzf;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::cmp::{max, min};
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::sync::Arc;
use std::sync::RwLock;

thread_local! {
    static EDIT_ALIGNER: RefCell<Option<AffineWavefronts>> = const { RefCell::new(None) };
}

/// Parse a CIGAR string into a vector of CigarOp
// Note that the query_delta is negative for reverse strand alignments
#[derive(Debug, Clone, PartialEq)]
pub struct CigarOp {
    val: u32,
}

impl CigarOp {
    pub fn new(len: i32, op: char) -> Self {
        let val = match op {
            '=' => 0,
            'X' => 1,
            'I' => 2,
            'D' => 3,
            'M' => 4,
            _ => panic!("Invalid CIGAR operation: {op}"),
        };
        Self {
            val: (val << 29) | (len as u32),
        }
    }

    pub fn op(&self) -> char {
        // two most significant bits in the val tell us the op
        match self.val >> 29 {
            0 => '=',
            1 => 'X',
            2 => 'I',
            3 => 'D',
            4 => 'M',
            _ => panic!("Invalid CIGAR operation: {}", self.val >> 29),
        }
    }

    pub fn len(&self) -> i32 {
        (self.val & ((1 << 29) - 1)) as i32
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn target_delta(&self) -> i32 {
        match self.op() {
            '=' | 'X' | 'D' | 'M' => self.len(),
            'I' => 0,
            _ => panic!("Invalid CIGAR operation: {}", self.op()),
        }
    }

    pub fn query_delta(&self, strand: Strand) -> i32 {
        match self.op() {
            '=' | 'X' | 'I' | 'M' => {
                if strand == Strand::Forward {
                    self.len()
                } else {
                    -self.len()
                }
            }
            'D' => 0,
            _ => panic!("Invalid CIGAR operation: {}", self.op()),
        }
    }

    fn adjust_len(&mut self, length_delta: i32) {
        self.val = (self.val & (7 << 29)) | ((self.len() + length_delta) as u32);
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct QueryMetadata {
    query_id: u32,
    target_start: i32,
    target_end: i32,
    target_contig_offset: i64,
    target_contig_len: i64,
    query_start: i32,
    query_end: i32,
    query_contig_offset: i64,
    query_contig_len: i64,
    alignment_file_index: u16,
    strand_and_data_offset: u64, // Track strand and cigar/tracepoints offset
    data_bytes: usize,
}

impl QueryMetadata {
    // Constants for bit manipulation
    const STRAND_BIT: u64 = 0x8000000000000000; // Most significant bit for u64

    fn strand(&self) -> Strand {
        if (self.strand_and_data_offset & Self::STRAND_BIT) != 0 {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }

    fn data_offset(&self) -> u64 {
        self.strand_and_data_offset & !Self::STRAND_BIT
    }

    /// Determine file type based on extension
    fn get_file_type(alignment_file: &str) -> FileType {
        if alignment_file.ends_with(".1aln") {
            FileType::OneAln
        } else {
            FileType::Paf
        }
    }

    /// Get alignment data bytes (CIGAR or tracepoints)
    fn get_data_bytes(&self, alignment_files: &[String]) -> Result<Vec<u8>, String> {
        // Get the correct file and type
        let alignment_file = &alignment_files[self.alignment_file_index as usize];
        let file_type = Self::get_file_type(alignment_file);

        // Check if data is available
        if self.data_bytes == 0 {
            let expected_data = match file_type {
                FileType::Paf => "CIGAR strings ('cg:Z' tag)",
                FileType::OneAln => "tracepoints",
            };
            return Err(format!(
                "The alignment file '{}' does not contain {}.",
                alignment_file, expected_data
            ));
        }

        let data_buffer = match file_type {
            FileType::OneAln => {
                // For .1aln files, use OneAlnParser to seek to the alignment

                // The strand_and_data_offset contains the alignment_index (not a file offset)
                let alignment_index = self.data_offset(); // This is actually alignment_index for .1aln

                let parser = OneAlnParser::new(alignment_file.clone()).map_err(|e| {
                    format!(
                        "Failed to create OneAlnParser for '{}': {:?}",
                        alignment_file, e
                    )
                })?;

                let alignment = parser.seek_alignment(alignment_index).map_err(|e| {
                    format!(
                        "Failed to seek to alignment {} in '{}': {:?}",
                        alignment_index, alignment_file, e
                    )
                })?;

                // Convert tracepoints to string format for consistency
                let tracepoints_str = alignment
                    .trace_diffs
                    .iter()
                    .zip(alignment.tracepoints.iter())
                    .map(|(&diff, &tp)| format!("{diff},{tp}"))
                    .collect::<Vec<_>>()
                    .join(";");

                tracepoints_str.into_bytes()
            }
            FileType::Paf => {
                // For PAF files, use byte-level seeking to read CIGAR string

                // Allocate space for cigar
                let mut data_buffer = vec![0; self.data_bytes];

                if [".gz", ".bgz"].iter().any(|e| alignment_file.ends_with(e)) {
                    // For compressed files, use virtual position directly
                    let mut reader =
                        bgzf::io::Reader::new(File::open(alignment_file).map_err(|e| {
                            format!("Failed to open compressed file '{}': {}", alignment_file, e)
                        })?);
                    let virtual_position = bgzf::VirtualPosition::from(self.data_offset());
                    reader.seek(virtual_position).map_err(|e| {
                        format!(
                            "Failed to seek in compressed file '{}': {}",
                            alignment_file, e
                        )
                    })?;
                    reader.read_exact(&mut data_buffer).map_err(|e| {
                        format!(
                            "Failed to read data from compressed file '{}': {}",
                            alignment_file, e
                        )
                    })?;
                } else {
                    // For uncompressed files, use byte offset
                    let mut reader = File::open(alignment_file)
                        .map_err(|e| format!("Failed to open file '{}': {}", alignment_file, e))?;
                    reader
                        .seek(SeekFrom::Start(self.data_offset()))
                        .map_err(|e| {
                            format!("Failed to seek in file '{}': {}", alignment_file, e)
                        })?;
                    reader.read_exact(&mut data_buffer).map_err(|e| {
                        format!("Failed to read data from file '{}': {}", alignment_file, e)
                    })?;
                }

                data_buffer
            }
        };

        Ok(data_buffer)
    }

    fn get_cigar_ops_from_bytes(&self, data_bytes: Vec<u8>) -> Vec<CigarOp> {
        let cigar_str = std::str::from_utf8(&data_bytes).unwrap();
        match parse_cigar_to_delta(cigar_str) {
            Ok(ops) => ops,
            Err(e) => panic!(
                "Failed to parse CIGAR string '{cigar_str}' in QueryMetadata: {:?}",
                e
            ),
        }
    }

    fn get_tracepoints_from_bytes(&self, data_bytes: Vec<u8>) -> Vec<(usize, usize)> {
        let tracepoints_str = std::str::from_utf8(&data_bytes).unwrap();
        match parse_tracepoints(tracepoints_str) {
            Ok(tp) => tp,
            Err(e) => panic!(
                "Failed to parse tracepoints '{tracepoints_str}' in QueryMetadata: {:?}",
                e
            ),
        }
    }
}

pub type AdjustedInterval = (Interval<u32>, Vec<CigarOp>, Interval<u32>);
type TreeMap = FxHashMap<u32, Arc<BasicCOITree<QueryMetadata, u32>>>;

#[derive(Debug, Clone, Copy, PartialEq)]
enum FileType {
    Paf,
    OneAln,
}

#[derive(Serialize, Deserialize)]
struct SerializableInterval {
    first: i32,
    last: i32,
    metadata: QueryMetadata,
}

#[derive(Default, Clone)]
pub struct SortedRanges {
    pub ranges: Vec<(i32, i32)>,
    sequence_length: i32,
    min_distance: i32,
}

impl SortedRanges {
    pub fn new(sequence_length: i32, min_distance: i32) -> Self {
        Self {
            ranges: Vec::new(),
            sequence_length,
            min_distance,
        }
    }

    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, (i32, i32)> {
        self.ranges.iter()
    }

    pub fn insert(&mut self, new_range: (i32, i32)) -> Vec<(i32, i32)> {
        let (mut start, mut end) = if new_range.0 <= new_range.1 {
            (new_range.0, new_range.1)
        } else {
            (new_range.1, new_range.0)
        };

        // Find nearby ranges to potentially merge with
        let mut i = match self.ranges.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(pos) => pos,
            Err(pos) => pos,
        };

        // Check previous range
        if i > 0 && (start - self.ranges[i - 1].1).abs() < self.min_distance {
            start = self.ranges[i - 1].1;
            i -= 1;
        } else if start < self.min_distance {
            start = 0;
        }

        // Check next range
        if i < self.ranges.len() && (self.ranges[i].0 - end).abs() < self.min_distance {
            end = self.ranges[i].0;
        } else if end > (self.sequence_length - self.min_distance) {
            end = self.sequence_length;
        }

        // Return regions that don't overlap with existing ranges
        let mut non_overlapping = Vec::new();
        let mut current = start;

        // Find the first range that could overlap
        let mut i = match self.ranges.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(pos) => pos,
            Err(pos) => pos,
        };

        // Check previous range for overlap
        if i > 0 && self.ranges[i - 1].1 > start {
            i -= 1;
        }

        // Process all potentially overlapping ranges
        while i < self.ranges.len() && current < end {
            let (range_start, range_end) = self.ranges[i];
            if range_start > end {
                break;
            }
            if current < range_start {
                non_overlapping.push((current, range_start));
            }
            current = max(current, range_end);
            i += 1;
        }

        if current < end {
            non_overlapping.push((current, end));
        }

        match self.ranges.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(pos) | Err(pos) => {
                // Check if we can merge with the previous range
                if pos > 0 && self.ranges[pos - 1].1 >= start {
                    self.ranges[pos - 1].1 = max(self.ranges[pos - 1].1, end);
                    self.merge_forward_from(pos - 1);
                } else if pos < self.ranges.len() && end >= self.ranges[pos].0 {
                    self.ranges[pos].0 = min(start, self.ranges[pos].0);
                    self.ranges[pos].1 = max(end, self.ranges[pos].1);
                    self.merge_forward_from(pos);
                } else {
                    self.ranges.insert(pos, (start, end));
                }
            }
        }

        // let adjusted_new_range = if new_range.0 <= new_range.1 {
        //     (start, end)
        // } else {
        //     (end, start)
        // };

        non_overlapping
    }

    pub fn merge_forward_from(&mut self, start_idx: usize) {
        let mut write = start_idx;
        let mut read = start_idx + 1;
        while read < self.ranges.len() {
            if self.ranges[write].1 >= self.ranges[read].0 {
                self.ranges[write].1 = max(self.ranges[write].1, self.ranges[read].1);
            } else {
                write += 1;
                self.ranges.swap(write, read);
            }
            read += 1;
        }
        self.ranges.truncate(write + 1);
    }
}

/// Information about a contig within a scaffold for .1aln files
#[derive(Debug, Clone, Copy)]
struct ContigInfo {
    sbeg: i64, // Scaffold begin (position where this contig starts)
    clen: i64, // Contig length
}

pub struct Impg {
    pub trees: RwLock<TreeMap>,
    pub seq_index: SequenceIndex,
    alignment_files: Vec<String>, // List of all alignment files (PAF or 1aln)
    pub forest_map: ForestMap,    // Forest map for lazy loading
    index_file_path: String,      // Path to the index file for lazy loading
    scaffold_contigs: FxHashMap<u32, Vec<ContigInfo>>, // Aggregated contig info per scaffold across files
    trace_spacing_by_file: Vec<u32>,                   // Per-file trace spacing for .1aln files
}

/// Execute a closure with a thread-local edit distance mode aligner
fn with_edit_aligner<F, R>(f: F) -> R
where
    F: FnOnce(&mut AffineWavefronts) -> R,
{
    EDIT_ALIGNER.with(|aligner_cell| {
        let mut aligner_opt = aligner_cell.borrow_mut();
        if aligner_opt.is_none() {
            let distance_mode = DistanceMode::Edit {
                mismatch: 1,
                gap_opening: 1,
            };
            *aligner_opt = Some(distance_mode.create_aligner());
        }
        f(aligner_opt.as_mut().unwrap())
    })
}

impl Impg {
    fn collect_contig_metadata(
        alignment_files: &[String],
        seq_index: &SequenceIndex,
    ) -> (FxHashMap<u32, Vec<ContigInfo>>, Vec<u32>) {
        let mut scaffold_contigs: FxHashMap<u32, Vec<ContigInfo>> = FxHashMap::default();
        let mut trace_spacing_by_file: Vec<u32> = Vec::with_capacity(alignment_files.len());

        for alignment_file in alignment_files {
            let mut file_trace_spacing = 100u32;

            if alignment_file.ends_with(".1aln") {
                if let Ok(parser) = OneAlnParser::new(alignment_file.clone()) {
                    file_trace_spacing = parser.get_trace_spacing();
                    let seq_names = parser.get_sequence_names();
                    let contig_offsets_map = parser.get_contig_offsets();

                    for (&contig_id, &(sbeg, clen)) in contig_offsets_map {
                        let scaffold_name = seq_names.get(&contig_id).unwrap();
                        let scaffold_id = seq_index.get_id(scaffold_name).unwrap();
                        let contig_info = ContigInfo { sbeg, clen };
                        scaffold_contigs
                            .entry(scaffold_id)
                            .or_insert_with(Vec::new)
                            .push(contig_info);
                    }
                }
            }

            trace_spacing_by_file.push(file_trace_spacing);
        }

        for contigs in scaffold_contigs.values_mut() {
            contigs.sort_by_key(|c| c.sbeg);
            contigs.dedup_by(|a, b| a.sbeg == b.sbeg && a.clen == b.clen);
        }

        (scaffold_contigs, trace_spacing_by_file)
    }

    /// Get contig layout for a sequence, falling back to a single contig if no data is available
    fn sequence_contigs(&self, sequence_id: u32) -> Vec<ContigInfo> {
        if let Some(contigs) = self.scaffold_contigs.get(&sequence_id) {
            if !contigs.is_empty() {
                return contigs.clone();
            }
        }

        let len = self
            .seq_index
            .get_len_from_id(sequence_id)
            .unwrap_or_else(|| {
                panic!(
                    "Missing sequence length for scaffold id {sequence_id} while building default contig layout"
                )
            }) as i64;
        vec![ContigInfo { sbeg: 0, clen: len }]
    }

    /// Convert a scaffold-level range into contig-level ranges with corresponding offsets
    fn scaffold_range_to_contig_ranges(
        &self,
        sequence_id: u32,
        range: (i32, i32),
    ) -> Vec<(i32, i32, i64)> {
        let mut start = range.0 as i64;
        let mut end = range.1 as i64;
        if start > end {
            std::mem::swap(&mut start, &mut end);
        }

        let mut ranges = Vec::new();
        for contig in self.sequence_contigs(sequence_id) {
            let contig_start = contig.sbeg;
            let contig_end = contig.sbeg + contig.clen;
            let overlap_start = start.max(contig_start);
            let overlap_end = end.min(contig_end);

            if overlap_start < overlap_end {
                let contig_range_start = (overlap_start - contig_start) as i32;
                let contig_range_end = (overlap_end - contig_start) as i32;
                ranges.push((contig_range_start, contig_range_end, contig_start));
            }
        }

        ranges
    }

    #[inline]
    fn to_scaffold_coordinate(&self, coord: i32, contig_offset: i64) -> i32 {
        let value = coord as i64 + contig_offset;
        value.clamp(i32::MIN as i64, i32::MAX as i64) as i32
    }

    fn process_tracepoints_data(
        &self,
        data_buffer: Vec<u8>,
        metadata: &QueryMetadata,
        target_id: u32,
        sequence_index: &UnifiedSequenceIndex,
        _penalties: (u8, u8, u8, u8, u8, u8),
        _requested_range: (i32, i32),
    ) -> (i32, i32, i32, i32, Vec<CigarOp>) {
        // Get the tracepoints
        let tracepoints = metadata.get_tracepoints_from_bytes(data_buffer);

        // Get alignment file index for looking up trace spacing
        let file_index = metadata.alignment_file_index as usize;

        // Get trace spacing from the specific file
        let trace_spacing = *self.trace_spacing_by_file.get(file_index).unwrap() as usize;

        // Fetch only the relevant portions of the sequences from scaffold coordinates
        let query_name = self.seq_index.get_name(metadata.query_id).unwrap();
        let query_fetch_start =
            self.to_scaffold_coordinate(metadata.query_start, metadata.query_contig_offset);
        let query_fetch_end =
            self.to_scaffold_coordinate(metadata.query_end, metadata.query_contig_offset);
        let query_seq_result =
            sequence_index.fetch_sequence(query_name, query_fetch_start, query_fetch_end);

        // Fetch target sequence
        let target_name = self.seq_index.get_name(target_id).unwrap();
        let target_fetch_start =
            self.to_scaffold_coordinate(metadata.target_start, metadata.target_contig_offset);
        let target_fetch_end =
            self.to_scaffold_coordinate(metadata.target_end, metadata.target_contig_offset);
        let target_seq_result = if metadata.strand() == Strand::Forward {
            sequence_index.fetch_sequence(target_name, target_fetch_start, target_fetch_end)
        } else {
            // Fetch and reverse complement for reverse strand alignments
            match sequence_index.fetch_sequence(target_name, target_fetch_start, target_fetch_end) {
                Ok(seq) => Ok(reverse_complement(&seq)),
                Err(e) => Err(e),
            }
        };

        // Determine contig-relative starting position compatible with the aligner
        let contig_target_start = if metadata.strand() == Strand::Forward {
            metadata.target_start
        } else {
            (metadata
                .target_contig_len
                .saturating_sub(metadata.target_end as i64)) as i32
        };

        // Convert tracepoints to CIGAR if we successfully fetched both sequences
        match (query_seq_result, target_seq_result) {
            (Ok(query_seq), Ok(target_seq)) => {
                let cigar_str = with_edit_aligner(|aligner| {
                    tracepoints_to_cigar_fastga_with_aligner(
                        &tracepoints,
                        trace_spacing,
                        &query_seq,
                        &target_seq,
                        metadata.query_start as usize,
                        contig_target_start as usize, // Use contig-level coordinates
                        metadata.strand() == Strand::Reverse,
                        aligner,
                    )
                });

                (
                    metadata.target_start,
                    metadata.target_end,
                    metadata.query_start,
                    metadata.query_end,
                    match parse_cigar_to_delta(&cigar_str) {
                        Ok(ops) => ops,
                        Err(e) => panic!(
                            "Failed to parse CIGAR string '{cigar_str}' during tracepoint processing: {:?}",
                            e
                        ),
                    },
                )
            }
            (Err(e), _) => {
                panic!("Failed to fetch query sequence: {e}");
            }
            (_, Err(e)) => {
                panic!("Failed to fetch target sequence: {e}");
            }
        }
    }

    pub fn from_multi_alignment_records(
        records_by_file: &[(Vec<AlignmentRecord>, String)],
        seq_index: SequenceIndex,
    ) -> Result<Self, ParseErr> {
        // Extract just the alignment file paths
        let alignment_files: Vec<String> = records_by_file
            .par_iter()
            .map(|(_, alignment_file)| alignment_file.clone())
            .collect();

        let (scaffold_contigs, trace_spacing_by_file) =
            Self::collect_contig_metadata(&alignment_files, &seq_index);

        let intervals: FxHashMap<u32, Vec<Interval<QueryMetadata>>> = records_by_file
            .par_iter()
            .enumerate() // Add enumeration to get the position as index
            .flat_map(|(file_index, (records, _))| {
                records
                    .par_iter()
                    .filter_map(|record| {
                        let query_metadata = QueryMetadata {
                            query_id: record.query_id,
                            target_start: record.target_start as i32,
                            target_end: record.target_end as i32,
                            target_contig_offset: record.target_contig_offset,
                            target_contig_len: record.target_contig_len,
                            query_start: record.query_start as i32,
                            query_end: record.query_end as i32,
                            query_contig_offset: record.query_contig_offset,
                            query_contig_len: record.query_contig_len,
                            alignment_file_index: file_index as u16,
                            strand_and_data_offset: record.strand_and_data_offset, // Already includes strand bit
                            data_bytes: record.data_bytes,
                        };

                        Some((
                            record.target_id,
                            Interval {
                                first: record.target_start as i32,
                                last: record.target_end as i32,
                                metadata: query_metadata,
                            },
                        ))
                    })
                    .collect::<Vec<_>>()
            })
            .fold(
                FxHashMap::default,
                |mut acc: FxHashMap<u32, Vec<Interval<QueryMetadata>>>, (target_id, interval)| {
                    acc.entry(target_id).or_default().push(interval);
                    acc
                },
            )
            .reduce(FxHashMap::default, |mut acc, part| {
                for (key, value) in part {
                    acc.entry(key).or_default().extend(value);
                }
                acc
            });

        let trees: TreeMap = intervals
            .into_par_iter()
            .map(|(target_id, interval_nodes)| {
                (
                    target_id,
                    Arc::new(BasicCOITree::new(interval_nodes.as_slice())),
                )
            })
            .collect();

        // Populate forest map with placeholder offsets for in-memory trees
        // This allows stats and other operations to work before serialization
        let mut forest_map = ForestMap::new();
        for &target_id in trees.keys() {
            forest_map.add_entry(target_id, 0); // Offset 0 = in-memory tree
        }

        Ok(Self {
            trees: RwLock::new(trees),
            seq_index,
            alignment_files,
            forest_map,
            index_file_path: String::new(), // All trees are in memory, no need for index file path
            scaffold_contigs,
            trace_spacing_by_file,
        })
    }

    /// Serialize the IMPG index to a writer with embedded forest map
    pub fn serialize_with_forest_map<W: std::io::Write + std::io::Seek>(
        &self,
        mut writer: W,
    ) -> std::io::Result<()> {
        const MAGIC: &[u8] = b"IMPGIDX1";
        let mut forest_map = ForestMap::new();

        // Write magic bytes
        writer.write_all(MAGIC)?;

        // Write placeholder for forest map offset (will be updated later)
        writer.write_all(&[0u8; 8])?;

        let mut current_offset = 16u64; // Header size

        // Serialize sequence index
        let seq_index_data =
            bincode::serde::encode_to_vec(&self.seq_index, bincode::config::standard()).map_err(
                |e| std::io::Error::other(format!("Failed to encode sequence index: {e:?}")),
            )?;

        writer.write_all(&seq_index_data)?;
        current_offset += seq_index_data.len() as u64;

        // Serialize each tree and track its offset
        let trees = self.trees.read().unwrap();
        for (&target_id, tree) in trees.iter() {
            // Record the offset before serializing this tree
            forest_map.add_entry(target_id, current_offset);

            // Convert tree to serializable format
            let intervals: Vec<SerializableInterval> = tree
                .iter()
                .map(|interval| SerializableInterval {
                    first: interval.first,
                    last: interval.last,
                    metadata: interval.metadata.clone(),
                })
                .collect();

            // Serialize the target_id and intervals
            let tree_data =
                bincode::serde::encode_to_vec(&(target_id, intervals), bincode::config::standard())
                    .map_err(|e| {
                        std::io::Error::other(format!(
                            "Failed to encode tree for target {target_id}: {e:?}"
                        ))
                    })?;

            writer.write_all(&tree_data)?;
            current_offset += tree_data.len() as u64;
        }

        // Now write the forest map at the end
        let forest_map_offset = current_offset;
        let forest_map_data =
            bincode::serde::encode_to_vec(&forest_map, bincode::config::standard()).map_err(
                |e| std::io::Error::other(format!("Failed to encode forest map: {e:?}")),
            )?;

        writer.write_all(&forest_map_data)?;

        // Go back and update the forest map offset in the header
        writer.seek(SeekFrom::Start(8))?;
        writer.write_all(&forest_map_offset.to_le_bytes())?;

        Ok(())
    }

    /// Load a specific tree from disk using the forest map
    fn load_tree_from_disk(&self, target_id: u32) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        if let Some(tree_offset) = self.forest_map.get_tree_offset(target_id) {
            let mut file = File::open(&self.index_file_path).expect("Failed to open index file");
            file.seek(std::io::SeekFrom::Start(tree_offset))
                .expect("Failed to seek in index file");

            let mut reader = BufReader::new(file);
            let (loaded_target_id, intervals): (u32, Vec<SerializableInterval>) =
                bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                    .unwrap_or_else(|_| {
                        panic!("Failed to deserialize tree for target {target_id}")
                    });

            // Verify we loaded the correct tree
            if loaded_target_id != target_id {
                panic!("Tree mismatch: expected {target_id}, got {loaded_target_id}");
            }

            // Reconstruct the tree
            let tree = BasicCOITree::new(
                intervals
                    .into_iter()
                    .map(|interval| Interval {
                        first: interval.first,
                        last: interval.last,
                        metadata: interval.metadata, // Move instead of clone
                    })
                    .collect::<Vec<_>>()
                    .as_slice(),
            );

            let arc_tree = Arc::new(tree);

            // Cache the tree for future use
            self.trees
                .write()
                .unwrap()
                .insert(target_id, Arc::clone(&arc_tree));

            Some(arc_tree)
        } else {
            None // Tree not found in forest map
        }
    }

    /// Get a tree from memory or load it from disk if necessary
    pub fn get_or_load_tree(
        &self,
        target_id: u32,
    ) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        // First check if the tree is already in memory
        if let Some(tree) = self.trees.read().unwrap().get(&target_id) {
            // Get a clone of the Arc<tree> (incrementing the reference count, a cheap operation)
            // We clone to avoid holding the RwLock for the duration of the operation using the tree

            return Some(Arc::clone(tree));
        }

        // Not in memory - try to load from disk
        self.load_tree_from_disk(target_id)
    }

    /// Load IMPG index from the format with embedded forest map at the end
    pub fn load_from_file<R: std::io::Read + std::io::Seek>(
        mut reader: R,
        alignment_files: &[String],
        index_file_path: String,
    ) -> std::io::Result<Self> {
        const MAGIC: &[u8] = b"IMPGIDX1";

        // Read and verify magic bytes
        let mut magic_buf = [0u8; 8];
        reader.read_exact(&mut magic_buf)?;
        if magic_buf != MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Invalid magic bytes - not a valid IMPG index file",
            ));
        }

        // Read forest map offset
        let mut offset_buf = [0u8; 8];
        reader.read_exact(&mut offset_buf)?;
        let forest_map_offset = u64::from_le_bytes(offset_buf);

        // Read sequence index
        let seq_index: SequenceIndex =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load sequence index: {e:?}"),
                    )
                })?;

        // Seek to forest map and read it
        reader.seek(SeekFrom::Start(forest_map_offset))?;
        let forest_map: ForestMap =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load forest map: {e:?}"),
                    )
                })?;

        let (scaffold_contigs, trace_spacing_by_file) =
            Self::collect_contig_metadata(alignment_files, &seq_index);

        Ok(Self {
            trees: RwLock::new(FxHashMap::default()), // Start with empty trees - load on demand
            seq_index,
            alignment_files: alignment_files.to_vec(),
            forest_map,
            index_file_path,
            scaffold_contigs,
            trace_spacing_by_file,
        })
    }

    pub fn query(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        penalties: Option<(u8, u8, u8, u8, u8, u8)>,
    ) -> Vec<AdjustedInterval> {
        let mut results = Vec::new();
        // Add the input range to the results
        results.push((
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
            if store_cigar {
                vec![CigarOp::new(range_end - range_start, '=')]
            } else {
                Vec::new()
            },
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
        ));

        debug!(
            "Querying region: {}:{}-{}, len: {}",
            self.seq_index.get_name(target_id).unwrap(),
            range_start,
            range_end,
            range_end - range_start
        );

        // Use default penalties if not provided
        let (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) =
            penalties.unwrap_or((0, 4, 6, 2, 26, 1));

        // Get or load the tree - if None, no overlaps exist for this target
        if let Some(tree) = self.get_or_load_tree(target_id) {
            let contig_ranges =
                self.scaffold_range_to_contig_ranges(target_id, (range_start, range_end));

            for (contig_start, contig_end, contig_offset) in contig_ranges {
                tree.query(contig_start, contig_end, |interval| {
                    let metadata = &interval.metadata;

                    if metadata.target_contig_offset != contig_offset {
                        return;
                    }

                    let alignment_file =
                        &self.alignment_files[metadata.alignment_file_index as usize];

                    // Get data bytes with all validation checks
                    let data_buffer = metadata
                        .get_data_bytes(&self.alignment_files)
                        .unwrap_or_else(|e| {
                            panic!("{}", e);
                        });

                    let (
                        adj_target_start,
                        adj_target_end,
                        adj_query_start,
                        adj_query_end,
                        cigar_ops,
                    ) = if QueryMetadata::get_file_type(alignment_file) == FileType::OneAln {
                        self.process_tracepoints_data(
                            data_buffer,
                            metadata,
                            target_id,
                            sequence_index.unwrap(),
                            (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2),
                            (contig_start, contig_end),
                        )
                    } else {
                        // Handle regular CIGAR
                        (
                            metadata.target_start,
                            metadata.target_end,
                            metadata.query_start,
                            metadata.query_end,
                            metadata.get_cigar_ops_from_bytes(data_buffer),
                        )
                    };

                    // Project through alignment in contig coordinates
                    let result = project_target_range_through_alignment(
                        (contig_start, contig_end),
                        (
                            adj_target_start,
                            adj_target_end,
                            adj_query_start,
                            adj_query_end,
                            metadata.strand(),
                        ),
                        &cigar_ops,
                    );

                    if let Some((
                        adjusted_query_start,
                        adjusted_query_end,
                        adjusted_cigar,
                        adjusted_target_start,
                        adjusted_target_end,
                    )) = result
                    {
                        // Check gap-compressed identity if threshold is provided
                        if let Some(threshold) = min_gap_compressed_identity {
                            if calculate_gap_compressed_identity(&adjusted_cigar) < threshold {
                                return; // Skip this result
                            }
                        }

                        let adjusted_interval = (
                            Interval {
                                first: self.to_scaffold_coordinate(
                                    adjusted_query_start,
                                    metadata.query_contig_offset,
                                ),
                                last: self.to_scaffold_coordinate(
                                    adjusted_query_end,
                                    metadata.query_contig_offset,
                                ),
                                metadata: metadata.query_id,
                            },
                            if store_cigar {
                                adjusted_cigar
                            } else {
                                Vec::new()
                            },
                            Interval {
                                first: self.to_scaffold_coordinate(
                                    adjusted_target_start,
                                    metadata.target_contig_offset,
                                ),
                                last: self.to_scaffold_coordinate(
                                    adjusted_target_end,
                                    metadata.target_contig_offset,
                                ),
                                metadata: target_id,
                            },
                        );
                        results.push(adjusted_interval);
                    }
                });
            }
        }

        results
    }

    pub fn query_transitive_dfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        penalties: Option<(u8, u8, u8, u8, u8, u8)>,
    ) -> Vec<AdjustedInterval> {
        // Initialize visited ranges from masked regions if provided
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter().map(|(&k, v)| (k, (*v).clone())).collect()
        } else {
            (0..self.seq_index.len() as u32)
                .into_par_iter() // Use parallel iterator
                .map(|id| {
                    let len = self.seq_index.get_len_from_id(id).unwrap();
                    (id, SortedRanges::new(len as i32, 0))
                })
                .collect()
        };

        // Filter input range
        let filtered_input_range = visited_ranges
            .entry(target_id)
            .or_default()
            .insert((range_start, range_end));

        let mut results = Vec::new();
        let mut stack = Vec::new();

        for (filtered_start, filtered_end) in filtered_input_range {
            // Add the filtered input range(s) to the results
            results.push((
                Interval {
                    first: filtered_start,
                    last: filtered_end,
                    metadata: target_id,
                },
                if store_cigar {
                    vec![CigarOp::new(filtered_end - filtered_start, '=')]
                } else {
                    Vec::new()
                },
                Interval {
                    first: filtered_start,
                    last: filtered_end,
                    metadata: target_id,
                },
            ));

            // Add the filtered input range(s) to the stack
            if (filtered_start - filtered_end).abs() >= min_transitive_len {
                stack.push((target_id, filtered_start, filtered_end, 0u16));
            }
        }

        // Use default penalties if not provided
        let (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) =
            penalties.unwrap_or((0, 4, 6, 2, 26, 1));

        while let Some((
            current_target_id,
            current_target_start,
            current_target_end,
            current_depth,
        )) = stack.pop()
        {
            // Check if we've reached max depth
            if max_depth > 0 && current_depth >= max_depth {
                continue;
            }

            debug!(
                "Querying region: {}:{}-{}, len: {}",
                self.seq_index.get_name(current_target_id).unwrap(),
                current_target_start,
                current_target_end,
                current_target_end - current_target_start
            );

            let prec_num_results = results.len();

            // Get or load the tree - if None, no overlaps exist for this target
            if let Some(tree) = self.get_or_load_tree(current_target_id) {
                // Collect intervals first to process them in parallel
                let mut intervals: Vec<(QueryMetadata, (i32, i32))> = Vec::new();

                for (contig_start, contig_end, contig_offset) in self
                    .scaffold_range_to_contig_ranges(
                        current_target_id,
                        (current_target_start, current_target_end),
                    )
                {
                    tree.query(contig_start, contig_end, |interval| {
                        if interval.metadata.target_contig_offset == contig_offset {
                            intervals.push((interval.metadata.clone(), (contig_start, contig_end)));
                        }
                    });
                }

                // Process the intervals in parallel
                let processed_results: Vec<_> = intervals
                    .into_par_iter()
                    .filter_map(|(metadata, requested_range)| {
                        let alignment_file =
                            &self.alignment_files[metadata.alignment_file_index as usize];

                        // Get data bytes with all validation checks
                        let data_buffer = metadata
                            .get_data_bytes(&self.alignment_files)
                            .unwrap_or_else(|e| {
                                panic!("{}", e);
                            });

                        let (
                            adj_target_start,
                            adj_target_end,
                            adj_query_start,
                            adj_query_end,
                            cigar_ops,
                        ) = if QueryMetadata::get_file_type(alignment_file) == FileType::OneAln {
                            self.process_tracepoints_data(
                                data_buffer,
                                &metadata,
                                current_target_id,
                                sequence_index.unwrap(),
                                (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2),
                                requested_range,
                            )
                        } else {
                            // Handle regular CIGAR
                            (
                                metadata.target_start,
                                metadata.target_end,
                                metadata.query_start,
                                metadata.query_end,
                                metadata.get_cigar_ops_from_bytes(data_buffer),
                            )
                        };

                        let result = project_target_range_through_alignment(
                            requested_range,
                            (
                                adj_target_start,
                                adj_target_end,
                                adj_query_start,
                                adj_query_end,
                                metadata.strand(),
                            ),
                            &cigar_ops,
                        );

                        if let Some((
                            adjusted_query_start,
                            adjusted_query_end,
                            adjusted_cigar,
                            adjusted_target_start,
                            adjusted_target_end,
                        )) = result
                        {
                            // Check gap-compressed identity if threshold is provided
                            if let Some(threshold) = min_gap_compressed_identity {
                                if calculate_gap_compressed_identity(&adjusted_cigar) < threshold {
                                    return None; // Skip this result
                                }
                            }

                            let adjusted_query_start_scaff = self.to_scaffold_coordinate(
                                adjusted_query_start,
                                metadata.query_contig_offset,
                            );
                            let adjusted_query_end_scaff = self.to_scaffold_coordinate(
                                adjusted_query_end,
                                metadata.query_contig_offset,
                            );
                            let adjusted_target_start_scaff = self.to_scaffold_coordinate(
                                adjusted_target_start,
                                metadata.target_contig_offset,
                            );
                            let adjusted_target_end_scaff = self.to_scaffold_coordinate(
                                adjusted_target_end,
                                metadata.target_contig_offset,
                            );

                            Some((
                                Interval {
                                    first: adjusted_query_start_scaff,
                                    last: adjusted_query_end_scaff,
                                    metadata: metadata.query_id,
                                },
                                if store_cigar {
                                    adjusted_cigar
                                } else {
                                    Vec::new()
                                },
                                Interval {
                                    first: adjusted_target_start_scaff,
                                    last: adjusted_target_end_scaff,
                                    metadata: current_target_id,
                                },
                                metadata.query_id,
                                adjusted_query_start_scaff,
                                adjusted_query_end_scaff,
                            ))
                        } else {
                            None
                        }
                    })
                    .collect();

                // Process results sequentially to maintain deterministic behavior
                for (
                    query_interval,
                    cigar,
                    target_interval,
                    query_id,
                    adjusted_query_start,
                    adjusted_query_end,
                ) in processed_results
                {
                    results.push((query_interval, cigar, target_interval));

                    // Only add non-overlapping portions to the stack for further exploration
                    if query_id != current_target_id {
                        let ranges = visited_ranges.entry(query_id).or_default();

                        let mut should_add = true;

                        // Check if the range is too close to any existing ranges
                        if min_distance_between_ranges > 0 {
                            let (new_min, new_max) = if adjusted_query_start <= adjusted_query_end {
                                (adjusted_query_start, adjusted_query_end)
                            } else {
                                (adjusted_query_end, adjusted_query_start)
                            };

                            // Find insertion point in sorted ranges
                            let idx = match ranges
                                .ranges
                                .binary_search_by_key(&new_min, |&(start, _)| start)
                            {
                                Ok(i) => i,
                                Err(i) => i,
                            };

                            // Only need to check adjacent ranges due to sorting
                            if idx > 0 {
                                // Check previous range
                                let (_, prev_end) = ranges.ranges[idx - 1];
                                if (new_min - prev_end).abs() < min_distance_between_ranges {
                                    should_add = false;
                                }
                            }
                            if idx < ranges.ranges.len() {
                                // Check next range
                                let (next_start, _) = ranges.ranges[idx];
                                if (next_start - new_max).abs() < min_distance_between_ranges {
                                    should_add = false;
                                }
                            }
                        }

                        if should_add {
                            let new_ranges =
                                ranges.insert((adjusted_query_start, adjusted_query_end));

                            // Add non-overlapping portions to stack
                            for (new_start, new_end) in new_ranges {
                                if (new_end - new_start).abs() >= min_transitive_len {
                                    stack.push((query_id, new_start, new_end, current_depth + 1));
                                }
                            }
                        }
                    }
                }
            }

            debug!("Collected {} results", results.len() - prec_num_results);

            // Merge contiguous/overlapping ranges with same sequence_id
            let stack_size = stack.len();
            stack.par_sort_by_key(|(id, start, _, _)| (*id, *start));

            let mut write = 0;
            for read in 1..stack.len() {
                if stack[write].0 == stack[read].0 &&   // Same sequence_id 
                    stack[write].2 >= stack[read].1
                // Overlapping or contiguous
                {
                    // Merge by extending end
                    stack[write].2 = stack[write].2.max(stack[read].2);
                } else {
                    write += 1;
                    stack.swap(write, read);
                }
            }
            stack.truncate(write + 1);
            debug!("Merged stack size from {} to {}", stack_size, stack.len());
        }

        results
    }

    pub fn query_transitive_bfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        penalties: Option<(u8, u8, u8, u8, u8, u8)>,
    ) -> Vec<AdjustedInterval> {
        // Initialize visited ranges from masked regions if provided
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter().map(|(&k, v)| (k, (*v).clone())).collect()
        } else {
            (0..self.seq_index.len() as u32)
                .into_par_iter() // Use parallel iterator
                .map(|id| {
                    let len = self.seq_index.get_len_from_id(id).unwrap();
                    (id, SortedRanges::new(len as i32, 0))
                })
                .collect()
        };

        // Filter input range
        let filtered_input_range = visited_ranges
            .entry(target_id)
            .or_default()
            .insert((range_start, range_end));

        let mut results = Vec::new();

        // Add the filtered input range(s) to the results
        for (filtered_start, filtered_end) in &filtered_input_range {
            results.push((
                Interval {
                    first: *filtered_start,
                    last: *filtered_end,
                    metadata: target_id,
                },
                if store_cigar {
                    vec![CigarOp::new(filtered_end - filtered_start, '=')]
                } else {
                    Vec::new()
                },
                Interval {
                    first: *filtered_start,
                    last: *filtered_end,
                    metadata: target_id,
                },
            ));
        }

        // Initialize ranges for first depth
        let mut current_depth = 0;
        let mut current_ranges = Vec::new();

        for (filtered_start, filtered_end) in filtered_input_range {
            if (filtered_start - filtered_end).abs() >= min_transitive_len {
                current_ranges.push((target_id, filtered_start, filtered_end));
            }
        }

        // Use default penalties if not provided
        let (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) =
            penalties.unwrap_or((0, 4, 6, 2, 26, 1));

        // Process by depth until max_depth or no more ranges
        while !current_ranges.is_empty() && (max_depth == 0 || current_depth < max_depth) {
            debug!(
                "Processing depth {} with {} ranges",
                current_depth,
                current_ranges.len()
            );

            // Process current depth ranges in parallel
            let query_results: Vec<Vec<(u32, i32, i32, Vec<CigarOp>, i32, i32, u32)>> =
                current_ranges
                    .par_iter()
                    .map(
                        |(current_target_id, current_target_start, current_target_end)| {
                            let mut local_results = Vec::new();

                            // Get or load the tree - if None, no overlaps exist for this target
                            if let Some(tree) = self.get_or_load_tree(*current_target_id) {
                                for (contig_start, contig_end, contig_offset) in self
                                    .scaffold_range_to_contig_ranges(
                                        *current_target_id,
                                        (*current_target_start, *current_target_end),
                                    )
                                {
                                    tree.query(contig_start, contig_end, |interval| {
                                        let metadata = &interval.metadata;

                                        if metadata.target_contig_offset != contig_offset {
                                            return;
                                        }

                                        let alignment_file = &self.alignment_files
                                            [metadata.alignment_file_index as usize];

                                        // Get data bytes with all validation checks
                                        let data_buffer = metadata
                                            .get_data_bytes(&self.alignment_files)
                                            .unwrap_or_else(|e| {
                                                panic!("{}", e);
                                            });

                                        let (
                                            adj_target_start,
                                            adj_target_end,
                                            adj_query_start,
                                            adj_query_end,
                                            cigar_ops,
                                        ) = if QueryMetadata::get_file_type(alignment_file)
                                            == FileType::OneAln
                                        {
                                            self.process_tracepoints_data(
                                                data_buffer,
                                                metadata,
                                                *current_target_id,
                                                sequence_index.unwrap(),
                                                (
                                                    _match, mismatch, gap_open1, gap_ext1,
                                                    gap_open2, gap_ext2,
                                                ),
                                                (contig_start, contig_end),
                                            )
                                        } else {
                                            // Handle regular CIGAR
                                            (
                                                metadata.target_start,
                                                metadata.target_end,
                                                metadata.query_start,
                                                metadata.query_end,
                                                metadata.get_cigar_ops_from_bytes(data_buffer),
                                            )
                                        };

                                        let result = project_target_range_through_alignment(
                                            (contig_start, contig_end),
                                            (
                                                adj_target_start,
                                                adj_target_end,
                                                adj_query_start,
                                                adj_query_end,
                                                metadata.strand(),
                                            ),
                                            &cigar_ops,
                                        );

                                        if let Some((
                                            adjusted_query_start,
                                            adjusted_query_end,
                                            adjusted_cigar,
                                            adjusted_target_start,
                                            adjusted_target_end,
                                        )) = result
                                        {
                                            // Check gap-compressed identity if threshold is provided
                                            if let Some(threshold) = min_gap_compressed_identity {
                                                if calculate_gap_compressed_identity(
                                                    &adjusted_cigar,
                                                ) < threshold
                                                {
                                                    return; // Skip this result
                                                }
                                            }

                                            let adjusted_query_start_scaff = self
                                                .to_scaffold_coordinate(
                                                    adjusted_query_start,
                                                    metadata.query_contig_offset,
                                                );
                                            let adjusted_query_end_scaff = self
                                                .to_scaffold_coordinate(
                                                    adjusted_query_end,
                                                    metadata.query_contig_offset,
                                                );
                                            let adjusted_target_start_scaff = self
                                                .to_scaffold_coordinate(
                                                    adjusted_target_start,
                                                    metadata.target_contig_offset,
                                                );
                                            let adjusted_target_end_scaff = self
                                                .to_scaffold_coordinate(
                                                    adjusted_target_end,
                                                    metadata.target_contig_offset,
                                                );

                                            let cigar_vec = if store_cigar {
                                                adjusted_cigar
                                            } else {
                                                Vec::new()
                                            };

                                            local_results.push((
                                                metadata.query_id,
                                                adjusted_query_start_scaff,
                                                adjusted_query_end_scaff,
                                                cigar_vec,
                                                adjusted_target_start_scaff,
                                                adjusted_target_end_scaff,
                                                *current_target_id,
                                            ));
                                        }
                                    });
                                }
                            }

                            local_results
                        },
                    )
                    .collect();

            // Prepare for next depth
            let mut next_depth_ranges = Vec::new();

            // Process results sequentially to update visited_ranges and results
            for query_result in query_results {
                for (
                    query_id,
                    adjusted_query_start,
                    adjusted_query_end,
                    adjusted_cigar,
                    adjusted_target_start,
                    adjusted_target_end,
                    current_target_id,
                ) in query_result
                {
                    // Add to results
                    results.push((
                        Interval {
                            first: adjusted_query_start,
                            last: adjusted_query_end,
                            metadata: query_id,
                        },
                        adjusted_cigar,
                        Interval {
                            first: adjusted_target_start,
                            last: adjusted_target_end,
                            metadata: current_target_id,
                        },
                    ));

                    // Only consider for next depth if it's a different sequence
                    if query_id != current_target_id {
                        let ranges = visited_ranges.entry(query_id).or_default();

                        let mut should_add = true;

                        // Check proximity to existing ranges
                        if min_distance_between_ranges > 0 {
                            let (new_min, new_max) = if adjusted_query_start <= adjusted_query_end {
                                (adjusted_query_start, adjusted_query_end)
                            } else {
                                (adjusted_query_end, adjusted_query_start)
                            };

                            // Find insertion point in sorted ranges
                            let idx = match ranges
                                .ranges
                                .binary_search_by_key(&new_min, |&(start, _)| start)
                            {
                                Ok(i) => i,
                                Err(i) => i,
                            };

                            // Only need to check adjacent ranges due to sorting
                            if idx > 0 {
                                // Check previous range
                                let (_, prev_end) = ranges.ranges[idx - 1];
                                if (new_min - prev_end).abs() < min_distance_between_ranges {
                                    should_add = false;
                                }
                            }

                            if should_add && idx < ranges.ranges.len() {
                                // Check next range
                                let (next_start, _) = ranges.ranges[idx];
                                if (next_start - new_max).abs() < min_distance_between_ranges {
                                    should_add = false;
                                }
                            }
                        }

                        if should_add {
                            let new_ranges =
                                ranges.insert((adjusted_query_start, adjusted_query_end));

                            // Add non-overlapping portions to next depth
                            for (new_start, new_end) in new_ranges {
                                if (new_end - new_start).abs() >= min_transitive_len {
                                    next_depth_ranges.push((query_id, new_start, new_end));
                                }
                            }
                        }
                    }
                }
            }

            // Move to next depth
            current_depth += 1;

            // Prepare ranges for next depth
            if !next_depth_ranges.is_empty() {
                // Sort and merge contiguous/overlapping ranges
                next_depth_ranges.par_sort_by_key(|(id, start, _)| (*id, *start));

                let mut write = 0;
                for read in 1..next_depth_ranges.len() {
                    if next_depth_ranges[write].0 == next_depth_ranges[read].0 &&  // Same sequence_id 
                       next_depth_ranges[write].2 >= next_depth_ranges[read].1
                    // Overlapping or contiguous
                    {
                        // Merge by extending end
                        next_depth_ranges[write].2 =
                            next_depth_ranges[write].2.max(next_depth_ranges[read].2);
                    } else {
                        write += 1;
                        next_depth_ranges.swap(write, read);
                    }
                }
                next_depth_ranges.truncate(write + 1);

                debug!(
                    "Next depth will process {} ranges after merging",
                    next_depth_ranges.len()
                );
            }

            // Set up for next iteration
            current_ranges = next_depth_ranges;
        }

        results
    }
}

fn project_target_range_through_alignment(
    requested_target_range: (i32, i32),
    record: (i32, i32, i32, i32, Strand),
    cigar_ops: &[CigarOp],
) -> Option<(i32, i32, Vec<CigarOp>, i32, i32)> {
    let (target_start, target_end, query_start, query_end, strand) = record;

    let dir = if strand == Strand::Forward { 1 } else { -1 };
    let mut query_pos = if strand == Strand::Forward {
        query_start
    } else {
        query_end
    };
    let mut target_pos = target_start;

    // Track CIGAR slice bounds
    let mut first_op_idx = 0;
    let mut last_op_idx = 0;
    let mut found_overlap = false;

    let mut projected_query_start = -1;
    let mut projected_query_end = -1;
    let mut projected_target_start = -1;
    let mut projected_target_end = -1;
    let mut first_op_offset = 0;
    let mut last_op_remaining = 0;

    // Calculate the last valid target position
    let last_target_pos = min(target_end, requested_target_range.1);

    for (curr_op_idx, cigar_op) in cigar_ops.iter().enumerate() {
        // If the target position is past the end of the range, we can stop
        if target_pos > last_target_pos {
            break;
        }

        match (cigar_op.target_delta(), cigar_op.query_delta(strand)) {
            (0, query_delta) => {
                // Insertion in query (deletions in target)
                if target_pos >= requested_target_range.0 {
                    if !found_overlap {
                        projected_query_start = query_pos;
                        projected_target_start = target_pos;
                        first_op_idx = curr_op_idx;
                        found_overlap = true;
                    }
                    projected_query_end = query_pos + query_delta;
                    projected_target_end = target_pos;
                    last_op_idx = curr_op_idx + 1;
                }
                query_pos += query_delta;
            }
            (target_delta, 0) => {
                // Deletion in query (insertions in target)
                let overlap_start = target_pos.max(requested_target_range.0);
                let overlap_end = (target_pos + target_delta).min(last_target_pos);

                if overlap_start < overlap_end {
                    if !found_overlap {
                        projected_query_start = query_pos;
                        projected_target_start = overlap_start;
                        first_op_idx = curr_op_idx;
                        first_op_offset = overlap_start - target_pos;
                        found_overlap = true;
                    }
                    projected_query_end = query_pos;
                    projected_target_end = overlap_end;
                    last_op_idx = curr_op_idx + 1;
                    last_op_remaining = overlap_end - (target_pos + target_delta);
                }
                target_pos += target_delta;
            }
            (target_delta, query_delta) => {
                // Match or mismatch
                let overlap_start = target_pos.max(requested_target_range.0);
                let overlap_end = (target_pos + target_delta).min(requested_target_range.1);

                if overlap_start < overlap_end {
                    let overlap_length = overlap_end - overlap_start;
                    let query_overlap_start = query_pos + (overlap_start - target_pos) * dir;
                    let query_overlap_end = query_overlap_start + overlap_length * dir;

                    if !found_overlap {
                        projected_query_start = query_overlap_start;
                        projected_target_start = overlap_start;
                        first_op_idx = curr_op_idx;
                        first_op_offset = overlap_start - target_pos;
                        found_overlap = true;
                    }
                    projected_query_end = query_overlap_end;
                    projected_target_end = overlap_end;
                    last_op_idx = curr_op_idx + 1;
                    last_op_remaining = overlap_end - (target_pos + target_delta);
                }

                target_pos += target_delta;
                query_pos += query_delta;
            }
        }
    }

    // If we had at least one overlap, the variables were set
    // projected_query_start == projected_query_end in deletions in the query
    // projected_target_start == projected_target_end in insertions in the query
    if found_overlap
        && projected_query_start != projected_query_end
        && projected_target_start != projected_target_end
    {
        let mut projected_cigar_ops = cigar_ops[first_op_idx..last_op_idx].to_vec();
        // Adjust first operation length
        if first_op_offset > 0 {
            projected_cigar_ops[0].adjust_len(-first_op_offset);
        }
        // Adjust last operation length
        if last_op_remaining < 0 {
            projected_cigar_ops[last_op_idx - first_op_idx - 1].adjust_len(last_op_remaining);
        }

        Some((
            projected_query_start,
            projected_query_end,
            projected_cigar_ops,
            projected_target_start,
            projected_target_end,
        ))
    } else {
        None
    }
}

fn parse_cigar_to_delta(cigar: &str) -> Result<Vec<CigarOp>, ParseErr> {
    let mut ops = Vec::new();
    let mut len: i32 = 0;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            len = len * 10 + (c as i32 - '0' as i32);
        } else {
            let op = CigarOp::new(len, c);
            ops.push(op);
            len = 0;
        }
    }

    Ok(ops)
}

fn parse_tracepoints(tp_str: &str) -> Result<Vec<(usize, usize)>, ParseErr> {
    tp_str
        .split(';')
        .filter(|s| !s.trim().is_empty())
        .map(|entry| {
            let mut parts = entry
                .split(',')
                .filter_map(|s| s.trim().parse::<usize>().ok());
            match (parts.next(), parts.next()) {
                (Some(x), Some(y)) if parts.next().is_none() => Ok((x, y)),
                _ => Err(ParseErr::InvalidFormat(
                    "Invalid tracepoint format".to_string(),
                )),
            }
        })
        .collect()
}

fn calculate_gap_compressed_identity(cigar_ops: &[CigarOp]) -> f64 {
    let (matches, mismatches, insertions, deletions) =
        cigar_ops
            .iter()
            .fold((0i32, 0i32, 0i32, 0i32), |(m, mm, i, d), op| {
                let len = op.len();
                match op.op() {
                    'M' | '=' => (m + len, mm, i, d), // Assume 'M' represents matches
                    'X' => (m, mm + len, i, d),
                    'I' => (m, mm, i + 1, d),
                    'D' => (m, mm, i, d + 1),
                    _ => (m, mm, i, d),
                }
            });

    let total = matches + mismatches + insertions + deletions;
    if total == 0 {
        0.0
    } else {
        (matches as f64) / (total as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paf::parse_paf;
    use std::io::BufReader;

    #[test]
    fn test_project_target_range_through_alignment_forward() {
        let target_range = (100, 200);
        let record = (100, 200, 0, 100, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let result =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();

        assert_eq!(result, (0, 100, cigar_ops.clone(), 100, 200));
    }

    #[test]
    fn test_project_target_range_through_alignment_reverse() {
        let target_range = (100, 200);
        let record = (100, 200, 0, 100, Strand::Reverse);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let result =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();

        assert_eq!(result, (100, 0, cigar_ops.clone(), 100, 200));
    }

    #[test]
    fn test_project_target_range_through_alignment() {
        let cigar_ops = vec![
            CigarOp::new(10, '='), // 10, 60
            CigarOp::new(5, 'I'),  // 10, 65
            CigarOp::new(5, 'D'),  // 15, 65
            CigarOp::new(50, '='), // 65, 115
            CigarOp::new(50, 'I'), // 65, 165
            CigarOp::new(35, '='), // 100, 200
        ];
        let base = (0, 100, 50, 200, Strand::Forward);
        {
            let result =
                project_target_range_through_alignment((0, 100), base, &cigar_ops).unwrap();
            assert_eq!(result, (50, 200, cigar_ops.clone(), 0, 100));
        }
        {
            let result =
                project_target_range_through_alignment((50, 55), base, &cigar_ops).unwrap();
            assert_eq!(result, (100, 105, vec![CigarOp::new(5, '=')], 50, 55));
        }
        {
            let result =
                project_target_range_through_alignment((50, 64), base, &cigar_ops).unwrap();
            assert_eq!(result, (100, 114, vec![CigarOp::new(14, '=')], 50, 64));
        }
        // We no longer output empty target ranges
        // {
        //     let result = project_target_range_through_alignment((65, 65), base, &cigar_ops).unwrap();
        //     assert_eq!(result, (115, 165, vec![CigarOp::new(50, 'I')], 65, 65));
        // }
        {
            let result =
                project_target_range_through_alignment((50, 65), base, &cigar_ops).unwrap();
            let cigar_ops = vec![CigarOp::new(15, '='), CigarOp::new(50, 'I')];
            assert_eq!(result, (100, 165, cigar_ops, 50, 65));
        }
        {
            let result =
                project_target_range_through_alignment((50, 66), base, &cigar_ops).unwrap();
            let cigar_ops = vec![
                CigarOp::new(15, '='),
                CigarOp::new(50, 'I'),
                CigarOp::new(1, '='),
            ];
            assert_eq!(result, (100, 166, cigar_ops, 50, 66));
        }
        {
            let result =
                project_target_range_through_alignment((70, 95), base, &cigar_ops).unwrap();
            assert_eq!(result, (170, 195, vec![CigarOp::new(25, '=')], 70, 95));
        }
    }

    // 1. Simple Forward Projection
    #[test]
    fn test_forward_projection_simple() {
        let target_range = (100, 200);
        let record = (100, 200, 100, 200, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let (query_start, query_end, cigar, target_start, target_end) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!(
            (query_start, query_end, cigar, target_start, target_end),
            (100, 200, vec![CigarOp::new(100, '=')], 100, 200)
        );
    }

    // 2. Simple Reverse Projection
    #[test]
    fn test_reverse_projection_simple() {
        let target_range = (100, 200);
        let record = (100, 200, 100, 200, Strand::Reverse);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let (query_start, query_end, cigar, target_start, target_end) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!(
            (query_start, query_end, cigar, target_start, target_end),
            (200, 100, vec![CigarOp::new(100, '=')], 100, 200)
        ); // Adjust for reverse calculation
    }

    // 3. Forward Projection with Insertions
    #[test]
    fn test_forward_projection_with_insertions() {
        let target_range = (50, 150);
        let record = (50, 150, 50, 160, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(50, '='), // Match
            CigarOp::new(10, 'I'), // Insertion
            CigarOp::new(50, '='), // Match
        ];
        let (start, end, cigar, _, _) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!((start, end, cigar), (50, 160, cigar_ops));
    }

    // 4. Forward Projection with Deletions
    #[test]
    fn test_forward_projection_with_deletions() {
        let target_range = (50, 150);
        let record = (50, 150, 50, 140, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(50, '='), // Match
            CigarOp::new(10, 'D'), // Deletion
            CigarOp::new(40, '='), // Match
        ];
        let (start, end, cigar, _, _) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!((start, end, cigar), (50, 140, cigar_ops));
    }

    // 5. Reverse Projection with Mixed Operations
    #[test]
    fn test_reverse_projection_with_mixed_operations() {
        let target_range = (150, 250);
        let record = (100, 200, 200, 300, Strand::Reverse);
        let cigar_ops = vec![
            CigarOp::new(50, '='), // 150, 250
            CigarOp::new(10, 'D'), // 150, 260
            CigarOp::new(10, 'I'), // 150, 250
            CigarOp::new(40, '='), // 150, 250
        ];
        let (start, end, cigar, _, _) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        let cigar_ops = vec![
            CigarOp::new(10, 'D'), // 150, 260
            CigarOp::new(10, 'I'), // 150, 250
            CigarOp::new(40, '='), // 150, 250
        ];
        assert_eq!((start, end, cigar), (250, 200, cigar_ops));
    }

    // 6. Edge Case Projection
    #[test]
    fn test_edge_case_projection() {
        let target_range = (0, 10);
        let record = (0, 50, 0, 40, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(10, '='), // Match
            CigarOp::new(20, 'D'), // Deletion in target
            CigarOp::new(8, '='),  // Match
            CigarOp::new(1, 'X'),  // Match
            CigarOp::new(1, '='),  // Match
            CigarOp::new(10, 'I'), // Insertion in query
            CigarOp::new(10, '='), // Match
        ];
        let (query_start, query_end, cigar, target_start, target_end) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!(
            (query_start, query_end, cigar, target_start, target_end),
            (0, 10, vec![CigarOp::new(10, '=')], 0, 10)
        );
    }

    #[test]
    fn test_parse_cigar_to_delta_basic() {
        let cigar = "10=5I5D";
        let cigar_ops = vec![
            CigarOp::new(10, '='),
            CigarOp::new(5, 'I'),
            CigarOp::new(5, 'D'),
        ];
        let ops = parse_cigar_to_delta(cigar).unwrap();
        assert_eq!(ops, cigar_ops);
    }

    // #[test]
    // fn test_parse_cigar_to_delta_invalid() {
    //     let cigar = "10=5Q"; // Q is not a valid CIGAR operation
    //     assert!(parse_cigar_to_delta(cigar).is_err());
    // }

    #[test]
    fn test_parse_paf_valid() {
        let paf_data = b"seq1\t100\t10\t20\t+\tt1\t200\t30\t40\t10\t20\t255\tcg:Z:10M\n";
        let mut seq_index = SequenceIndex::new();
        let query_id = seq_index.get_or_insert_id("seq1", Some(100));
        let target_id = seq_index.get_or_insert_id("t1", Some(200));
        let reader = BufReader::new(&paf_data[..]);
        let expected_records = vec![
            AlignmentRecord {
                query_id,
                query_start: 10,
                query_end: 20,
                query_contig_offset: 0,
                query_contig_len: 100,
                target_id,
                target_start: 30,
                target_end: 40,
                target_contig_offset: 0,
                target_contig_len: 200,
                strand_and_data_offset: 45, // Forward strand
                data_bytes: 3,
            },
            // Add more test records as needed
        ];
        let records = parse_paf(reader, &mut seq_index).unwrap();
        assert_eq!(records, expected_records);
    }
}
