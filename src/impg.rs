use crate::alignment_record::{AlignmentRecord, Strand};
use crate::forest_map::ForestMap;
use crate::graph::reverse_complement;
use crate::impg_index::ImpgIndex;
use crate::onealn::{OneAlnAlignment, OneAlnParser, ParseErr as OneAlnParseErr};
use crate::paf::read_cigar_data;
use crate::seqidx::SequenceIndex;
use crate::sequence_index::SequenceIndex as _; // The as _ syntax imports the trait so its methods are available, but doesn't bring the name into scope (avoiding the naming conflict)
use crate::sequence_index::UnifiedSequenceIndex;
use crate::subset_filter::SubsetFilter;
use coitrees::{BasicCOITree, Interval, IntervalTree};
use tracepoints::tracepoints_to_cigar_fastga_with_aligner;
use lib_wfa2::affine_wavefront::{AffineWavefronts, Distance};
use log::{debug, warn};
use onecode::OneFile;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::cmp::{max, min};
use std::fs::File;
use std::io::{self, BufReader, Seek, SeekFrom};
use std::sync::{Arc, RwLock};

// use libc;
// fn log_memory_usage(label: &str) {
//     let mut usage = MaybeUninit::<libc::rusage>::uninit();
//     let result = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
//     if result == 0 {
//         let usage = unsafe { usage.assume_init() };
//         debug!(
//             "mem[{label}] max_rss_kb={} ixrss_kb={} idrss_kb={}",
//             usage.ru_maxrss, usage.ru_ixrss, usage.ru_idrss
//         );
//     } else {
//         debug!("mem[{label}] getrusage_failed code={}", result);
//     }
// }

thread_local! {
    static EDIT_ALIGNER: RefCell<Option<AffineWavefronts>> = const { RefCell::new(None) };
    static ONEALN_HANDLE: RefCell<Option<(usize, OneFile)>> = const { RefCell::new(None) };
    static TARGET_SEQ_CACHE: RefCell<Option<((u32, i32, i32, bool), Vec<u8>)>> = const { RefCell::new(None) };
}

/// Execute a closure with a thread-local edit distance mode aligner
fn with_edit_aligner<F, R>(f: F) -> R
where
    F: FnOnce(&mut AffineWavefronts) -> R,
{
    EDIT_ALIGNER.with(|aligner_cell| {
        let mut aligner_opt = aligner_cell.borrow_mut();
        if aligner_opt.is_none() {
            *aligner_opt = Some(Distance::Edit.create_aligner(None));
        }
        f(aligner_opt.as_mut().unwrap())
    })
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
    query_start: i32,
    query_end: i32,
    alignment_file_index: u32,
    strand_and_data_offset: u64, // Track strand and cigar/tracepoints offset
    data_bytes: usize,
}

impl QueryMetadata {
    // Constants for bit manipulation
    const STRAND_BIT: u64 = 0x8000000000000000; // Most significant bit for u64
    const REVERSED_BIT: u64 = 0x4000000000000000; // Second-most significant bit for reversed frame

    pub fn query_id(&self) -> u32 {
        self.query_id
    }

    pub fn query_start(&self) -> i32 {
        self.query_start
    }

    pub fn query_end(&self) -> i32 {
        self.query_end
    }

    fn strand(&self) -> Strand {
        if (self.strand_and_data_offset & Self::STRAND_BIT) != 0 {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }

    /// Returns true if this is a reversed index entry (query-indexed rather than target-indexed).
    /// Reversed entries project from query coords back to target coords.
    fn is_reversed(&self) -> bool {
        (self.strand_and_data_offset & Self::REVERSED_BIT) != 0
    }

    fn data_offset(&self) -> u64 {
        self.strand_and_data_offset & !(Self::STRAND_BIT | Self::REVERSED_BIT)
    }

    /// Determine file type based on extension
    fn get_file_type(alignment_file: &str) -> FileType {
        if alignment_file.ends_with(".1aln") {
            FileType::OneAln
        } else {
            FileType::Paf
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

/// Results from scanning tracepoints to find segments overlapping a target range
#[derive(Debug)]
struct SubsettingResult {
    /// Index of first overlapping tracepoint
    first_idx: usize,
    /// Index of last overlapping tracepoint (inclusive)
    last_idx: usize,
    /// Query position at start of first overlapping segment (before refinement)
    first_query_pos: i32,
    /// Query position at end of last overlapping segment (before refinement)
    last_query_pos: i32,
    /// Number of overlapping segments
    num_overlapping_segments: usize,
    /// Total matches accumulated across overlapping segments
    total_matches: f64,
    /// Total mismatches accumulated across overlapping segments
    total_mismatches: f64,
    /// Info about first overlapping segment: (query_pos, query_delta, seg_start, seg_end, abs_target_delta, trace_diffs)
    first_segment_info: (i32, i32, i32, i32, i32, i32),
    /// Info about last overlapping segment: (query_pos, query_delta, seg_start, seg_end, abs_target_delta, trace_diffs)
    last_segment_info: (i32, i32, i32, i32, i32, i32),
}

/// Results from scanning tracepoints to find segments overlapping a query range (inverse of SubsettingResult)
#[derive(Debug)]
struct InverseSubsettingResult {
    /// Index of first overlapping tracepoint
    first_idx: usize,
    /// Index of last overlapping tracepoint (inclusive)
    last_idx: usize,
    /// Target position at start of first overlapping segment (before refinement)
    first_target_pos: i32,
    /// Target position at end of last overlapping segment (before refinement)
    last_target_pos: i32,
    /// Number of overlapping segments
    num_overlapping_segments: usize,
    /// Total matches accumulated across overlapping segments
    total_matches: f64,
    /// Total mismatches accumulated across overlapping segments
    total_mismatches: f64,
    /// Info about first overlapping segment: (target_pos, target_delta, seg_query_start, seg_query_end, abs_query_delta, trace_diffs)
    first_segment_info: (i32, i32, i32, i32, i32, i32),
    /// Info about last overlapping segment: (target_pos, target_delta, seg_query_start, seg_query_end, abs_query_delta, trace_diffs)
    last_segment_info: (i32, i32, i32, i32, i32, i32),
}

pub struct Impg {
    pub trees: RwLock<TreeMap>,
    pub seq_index: SequenceIndex,
    alignment_files: Vec<String>,
    pub forest_map: ForestMap,
    index_file_path: String,
    pub sequence_files: Vec<String>,
    /// Cache for trace_spacing values per .1aln file (indexed by alignment_file_index)
    /// Uses Vec for O(1) direct indexing. None = not yet loaded or not a .1aln file.
    trace_spacing_cache: RwLock<Vec<Option<i64>>>,
}

impl Impg {
    /// Get cached trace_spacing for a .1aln file (lazy-loaded on first access)
    fn get_trace_spacing(&self, file_index: usize) -> Result<i64, String> {
        let file_path = &self.alignment_files[file_index];

        // Fast path: check if already cached (just array indexing - extremely fast!)
        {
            let cache = self.trace_spacing_cache.read().unwrap();
            if let Some(spacing) = cache[file_index] {
                return Ok(spacing);
            }
        }

        // Slow path: read from file and cache
        let spacing = OneAlnParser::read_trace_spacing(file_path)
            .map_err(|e| format!("Failed to read trace_spacing from '{}': {:?}", file_path, e))?;

        // Cache it
        let mut cache = self.trace_spacing_cache.write().unwrap();
        cache[file_index] = Some(spacing);

        Ok(spacing)
    }

    fn with_onealn_reader<F, R>(&self, file_index: usize, f: F) -> Result<R, String>
    where
        F: FnOnce(&mut OneFile) -> Result<R, OneAlnParseErr>,
    {
        let file_path = &self.alignment_files[file_index];
        ONEALN_HANDLE.with(|handle_cell| {
            let mut slot = handle_cell.borrow_mut();
            let needs_open = match slot.as_ref() {
                Some((idx, _)) => *idx != file_index,
                None => true,
            };

            if needs_open {
                *slot = None;
                let file = OneFile::open_read(file_path, None, None, 1)
                    .map_err(|e| format!("Failed to open 1aln file '{}': {}", file_path, e))?;
                *slot = Some((file_index, file));
            }

            let (_, handle) = slot
                .as_mut()
                .expect("ONEALN_HANDLE must contain a handle after opening");
            f(handle).map_err(|e| e.to_string())
        })
    }

    fn get_target_sequence_cached(
        &self,
        sequence_index: &UnifiedSequenceIndex,
        target_name: &str,
        target_id: u32,
        start: i32,
        end: i32,
        is_reverse: bool,
    ) -> Vec<u8> {
        let key = (target_id, start, end, is_reverse);
        TARGET_SEQ_CACHE.with(|cache_cell| {
            let mut slot = cache_cell.borrow_mut();
            if let Some((cached_key, cached_seq)) = slot.as_mut() {
                if *cached_key == key {
                    return cached_seq.clone();
                }
            }

            let seq = sequence_index
                .fetch_sequence(target_name, start, end)
                .unwrap_or_else(|e| panic!("Failed to fetch target sequence: {e}"));
            let seq = if is_reverse {
                reverse_complement(&seq)
            } else {
                seq
            };

            *slot = Some((key, seq.clone()));
            seq
        })
    }

    /// Get CIGAR operations for an alignment (handles both PAF and 1aln files)
    fn get_cigar_ops(
        &self,
        metadata: &QueryMetadata,
        target_id: u32,
        sequence_index: Option<&UnifiedSequenceIndex>,
    ) -> Vec<CigarOp> {
        let alignment_file = &self.alignment_files[metadata.alignment_file_index as usize];

        match QueryMetadata::get_file_type(alignment_file) {
            FileType::Paf => {
                // For PAF files, read CIGAR directly from file
                if metadata.data_bytes == 0 {
                    panic!(
                        "The alignment file '{}' does not contain CIGAR strings ('cg:Z' tag).",
                        alignment_file
                    );
                }

                let mut data_buffer = vec![0; metadata.data_bytes];
                read_cigar_data(alignment_file, metadata.data_offset(), &mut data_buffer)
                    .unwrap_or_else(|e| panic!("{}", e));

                // get_cigar_ops_from_bytes
                let cigar_str =
                    std::str::from_utf8(&data_buffer).expect("Failed to parse CIGAR data as UTF-8");
                parse_cigar_to_delta(cigar_str).unwrap_or_else(|e| {
                    panic!(
                        "Failed to parse CIGAR string '{}' in QueryMetadata: {:?}",
                        cigar_str, e
                    )
                })
            }
            FileType::OneAln => {
                // For 1aln files, convert tracepoints to CIGAR
                let alignment = self
                    .get_onealn_alignment(metadata)
                    .unwrap_or_else(|e| panic!("{}", e));

                self.process_tracepoints_data(
                    &alignment,
                    metadata,
                    target_id,
                    sequence_index.expect("Sequence index required for .1aln files"),
                )
            }
        }
    }

    fn get_onealn_alignment(&self, metadata: &QueryMetadata) -> Result<OneAlnAlignment, String> {
        let file_index = metadata.alignment_file_index as usize;
        let alignment_file = &self.alignment_files[file_index];

        if !alignment_file.ends_with(".1aln") {
            return Err(format!(
                "Alignment '{}' is not a .1aln file; tracepoints unavailable",
                alignment_file
            ));
        }

        // Get cached trace_spacing (O(1) array access after first load!)
        let trace_spacing = self.get_trace_spacing(file_index)?;

        let alignment_index = metadata.data_offset();
        self.with_onealn_reader(file_index, |file| {
            OneAlnParser::fetch_alignment_from_reader(file, trace_spacing, alignment_index)
        })
        .map_err(|e| {
            format!(
                "Failed to fetch alignment {} from '{}': {}",
                alignment_index, alignment_file, e
            )
        })
    }

    /// Scan tracepoints to find segments overlapping the requested target range.
    fn scan_overlapping_tracepoints(
        &self,
        alignment: &OneAlnAlignment,
        metadata: &QueryMetadata,
        range_start: i32,
        range_end: i32,
    ) -> Option<SubsettingResult> {
        let trace_spacing = alignment.trace_spacing as i32;
        let is_reverse = metadata.strand() == Strand::Reverse;

        // Calculate first boundary: FASTGA-style first segment length
        let query_start_in_contig = alignment.query_contig_start as i32;
        let first_boundary =
            ((query_start_in_contig / trace_spacing) + 1) * trace_spacing - query_start_in_contig;

        // Metadata stores query coords in FORWARD order for both strands
        let working_query_start = metadata.query_start;

        // Initialize scanning positions
        let mut target_pos = if is_reverse {
            metadata.target_end // Start from target end, scan backward
        } else {
            metadata.target_start
        };
        let mut query_pos = working_query_start; // Always start from working start
        let target_dir = if is_reverse { -1 } else { 1 };

        // Track overlapping segments
        let mut first_idx: Option<usize> = None;
        let mut last_idx: usize = 0;
        let mut first_query_pos: Option<i32> = None;
        let mut last_query_pos: Option<i32> = None;
        let mut first_segment_info: Option<(i32, i32, i32, i32, i32, i32)> = None;
        let mut last_segment_info: Option<(i32, i32, i32, i32, i32, i32)> = None;
        let mut num_overlapping_segments = 0;

        // Statistics for identity calculation
        let mut total_matches = 0.0;
        let mut total_mismatches = 0.0;

        // Scan tracepoints and collect those that overlap the requested range
        for (idx, &tracepoint) in alignment.tracepoints.iter().enumerate() {
            // Calculate deltas for this segment
            let query_delta = if idx == 0 {
                first_boundary
            } else {
                trace_spacing
            };
            let target_delta = (tracepoint as i32) * target_dir;
            let abs_target_delta = tracepoint.abs() as i32;

            // Segment boundaries (forward coordinate space)
            let seg_start = target_pos.min(target_pos + target_delta);
            let seg_end = target_pos.max(target_pos + target_delta);

            // Check overlap with requested range
            if seg_start < range_end && seg_end > range_start {
                num_overlapping_segments += 1;
                let num_diffs = alignment.trace_diffs[idx] as i32;
                let seg_info = (
                    query_pos,
                    query_delta,
                    seg_start,
                    seg_end,
                    abs_target_delta,
                    num_diffs,
                );

                // Track first and last overlapping segments
                if first_query_pos.is_none() {
                    first_idx = Some(idx);
                    first_query_pos = Some(query_pos);
                    first_segment_info = Some(seg_info);
                }
                last_idx = idx;
                last_query_pos = Some(query_pos + query_delta);
                last_segment_info = Some(seg_info);

                // Accumulate identity statistics
                let aligned_len = (query_delta.min(abs_target_delta) as f64).max(0.0);
                total_matches += (aligned_len - num_diffs as f64).max(0.0);
                total_mismatches += num_diffs as f64;
            }

            // Advance to next segment
            target_pos += target_delta;
            query_pos += query_delta;

            // Early exit when past requested range
            if (is_reverse && target_pos <= range_start) || (!is_reverse && target_pos >= range_end)
            {
                break;
            }
        }

        // Return None if no overlapping segments found
        first_idx?;

        Some(SubsettingResult {
            first_idx: first_idx.unwrap(),
            last_idx,
            first_query_pos: first_query_pos.unwrap(),
            last_query_pos: last_query_pos.unwrap(),
            num_overlapping_segments,
            total_matches,
            total_mismatches,
            first_segment_info: first_segment_info.unwrap(),
            last_segment_info: last_segment_info.unwrap(),
        })
    }

    /// Scan tracepoints to find segments overlapping the requested query range (inverse of scan_overlapping_tracepoints).
    /// Used for reversed index entries where we query by original query coordinates.
    fn scan_overlapping_tracepoints_by_query(
        &self,
        alignment: &OneAlnAlignment,
        metadata: &QueryMetadata,
        range_start: i32,
        range_end: i32,
    ) -> Option<InverseSubsettingResult> {
        let trace_spacing = alignment.trace_spacing as i32;
        let is_reverse = metadata.strand() == Strand::Reverse;

        // Calculate first boundary: FASTGA-style first segment length
        let query_start_in_contig = alignment.query_contig_start as i32;
        let first_boundary =
            ((query_start_in_contig / trace_spacing) + 1) * trace_spacing - query_start_in_contig;

        // Metadata stores query coords in FORWARD order for both strands
        let working_query_start = metadata.query_start;

        // Initialize scanning positions
        let mut target_pos = if is_reverse {
            metadata.target_end // Start from target end, scan backward
        } else {
            metadata.target_start
        };
        let mut query_pos = working_query_start; // Always start from working start
        let target_dir = if is_reverse { -1 } else { 1 };

        // Normalize query range for comparison (handle both orientations)
        let (query_range_low, query_range_high) = (range_start.min(range_end), range_start.max(range_end));

        // Track overlapping segments
        let mut first_idx: Option<usize> = None;
        let mut last_idx: usize = 0;
        let mut first_target_pos: Option<i32> = None;
        let mut last_target_pos: Option<i32> = None;
        let mut first_segment_info: Option<(i32, i32, i32, i32, i32, i32)> = None;
        let mut last_segment_info: Option<(i32, i32, i32, i32, i32, i32)> = None;
        let mut num_overlapping_segments = 0;

        // Statistics for identity calculation
        let mut total_matches = 0.0;
        let mut total_mismatches = 0.0;

        // Scan tracepoints and collect those that overlap the requested query range
        for (idx, &tracepoint) in alignment.tracepoints.iter().enumerate() {
            // Calculate deltas for this segment
            let query_delta = if idx == 0 {
                first_boundary
            } else {
                trace_spacing
            };
            let target_delta = (tracepoint as i32) * target_dir;
            let abs_target_delta = tracepoint.abs() as i32;

            // Query segment boundaries (in forward coordinate space)
            let seg_query_start = query_pos;
            let seg_query_end = query_pos + query_delta;

            // Check overlap with requested query range
            if seg_query_start < query_range_high && seg_query_end > query_range_low {
                num_overlapping_segments += 1;
                let num_diffs = alignment.trace_diffs[idx] as i32;

                // For inverse: track target positions, segment info stores target-centric data
                let seg_info = (
                    target_pos,
                    target_delta,
                    seg_query_start,
                    seg_query_end,
                    query_delta,
                    num_diffs,
                );

                // Track first and last overlapping segments
                if first_target_pos.is_none() {
                    first_idx = Some(idx);
                    first_target_pos = Some(target_pos);
                    first_segment_info = Some(seg_info);
                }
                last_idx = idx;
                last_target_pos = Some(target_pos + target_delta);
                last_segment_info = Some(seg_info);

                // Accumulate identity statistics
                let aligned_len = (query_delta.min(abs_target_delta) as f64).max(0.0);
                total_matches += (aligned_len - num_diffs as f64).max(0.0);
                total_mismatches += num_diffs as f64;
            }

            // Advance to next segment
            target_pos += target_delta;
            query_pos += query_delta;

            // Early exit when past requested query range
            if query_pos >= query_range_high {
                break;
            }
        }

        // Return None if no overlapping segments found
        first_idx?;

        Some(InverseSubsettingResult {
            first_idx: first_idx.unwrap(),
            last_idx,
            first_target_pos: first_target_pos.unwrap(),
            last_target_pos: last_target_pos.unwrap(),
            num_overlapping_segments,
            total_matches,
            total_mismatches,
            first_segment_info: first_segment_info.unwrap(),
            last_segment_info: last_segment_info.unwrap(),
        })
    }

    fn process_tracepoints_data(
        &self,
        alignment: &OneAlnAlignment,
        metadata: &QueryMetadata,
        target_id: u32,
        sequence_index: &UnifiedSequenceIndex,
    ) -> Vec<CigarOp> {
        let (query_start, query_end) = (metadata.query_start, metadata.query_end);
        let (target_start, target_end) = (metadata.target_start, metadata.target_end);

        // if there are no differences, we can shortcut to a perfect match CIGAR
        if alignment.differences == 0 {
            let match_len = query_end - query_start;
            return vec![CigarOp::new(match_len, '=')];
        }

        let tracepoints: Vec<(usize, usize)> = alignment
            .trace_diffs
            .iter()
            .zip(alignment.tracepoints.iter())
            .map(|(&diff, &tp)| (diff as usize, tp as usize))
            .collect();
        let trace_spacing = alignment.trace_spacing as usize;

        // Fetch query sequence (not cached)
        let query_name = self.seq_index.get_name(metadata.query_id).unwrap();
        let query_seq = sequence_index
            .fetch_sequence(query_name, query_start, query_end)
            .unwrap_or_else(|e| panic!("Failed to fetch query sequence: {e}"));

        let target_name = self.seq_index.get_name(target_id).unwrap();

        let is_reverse = metadata.strand() == Strand::Reverse;
        let target_seq = self.get_target_sequence_cached(
            sequence_index,
            target_name,
            target_id,
            target_start,
            target_end,
            is_reverse,
        );

        let cigar_str = with_edit_aligner(|aligner| {
            tracepoints_to_cigar_fastga_with_aligner(
                &tracepoints,
                trace_spacing.try_into().unwrap(),
                &query_seq,
                &target_seq,
                alignment.query_contig_start.try_into().unwrap(),
                alignment.target_contig_start.try_into().unwrap(),
                metadata.strand() == Strand::Reverse,
                aligner,
            )
        });

        match parse_cigar_to_delta(&cigar_str) {
            Ok(ops) => ops,
            Err(e) => panic!(
                "Failed to parse CIGAR string '{cigar_str}' during tracepoint processing: {:?}",
                e
            ),
        }
    }

    /// Reconstruct CIGAR only for the subset of tracepoints that overlap the requested range.
    fn process_subset_tracepoints(
        &self,
        alignment: &OneAlnAlignment,
        metadata: &QueryMetadata,
        target_id: u32,
        subset: &SubsettingResult,
        sequence_index: &UnifiedSequenceIndex,
    ) -> Vec<CigarOp> {
        let is_reverse = metadata.strand() == Strand::Reverse;

        // Extract subset tracepoints [first_idx..=last_idx]
        let subset_tracepoints: Vec<(usize, usize)> = alignment.trace_diffs
            [subset.first_idx..=subset.last_idx]
            .iter()
            .zip(alignment.tracepoints[subset.first_idx..=subset.last_idx].iter())
            .map(|(&diff, &tp)| (diff as usize, tp as usize))
            .collect();

        // Calculate subset sequence ranges (use RAW boundaries from subsetting scan, no refinement)
        let subset_query_start = subset.first_query_pos;

        // If last_idx is the actual last tracepoint of the alignment, the last
        // segment may be shorter than trace_spacing. Use metadata boundary.
        let subset_query_end = if subset.last_idx == alignment.tracepoints.len() - 1 {
            metadata.query_end // Last segment of alignment - use exact boundary
        } else {
            subset.last_query_pos // Middle segment - use calculated position
        };

        // For target: extract from first/last segment info (seg_start, seg_end already in forward space)
        let (_, _, first_seg_start, first_seg_end, _, _) = subset.first_segment_info;
        let (_, _, last_seg_start, last_seg_end, _, _) = subset.last_segment_info;
        let subset_target_start = first_seg_start
            .min(first_seg_end)
            .min(last_seg_start.min(last_seg_end));

        // If last_idx is the actual last tracepoint, use metadata target boundary
        let subset_target_end = if subset.last_idx == alignment.tracepoints.len() - 1 {
            metadata.target_end // Last segment of alignment - use exact boundary
        } else {
            first_seg_start
                .max(first_seg_end)
                .max(last_seg_start.max(last_seg_end)) // Middle segment - use calculated
        };

        // Fetch only subset sequences
        let query_name = self.seq_index.get_name(metadata.query_id).unwrap();
        let query_seq = sequence_index
            .fetch_sequence(query_name, subset_query_start, subset_query_end)
            .unwrap_or_else(|e| panic!("Failed to fetch subset query sequence: {e}"));

        let target_name = self.seq_index.get_name(target_id).unwrap();
        let target_seq = self.get_target_sequence_cached(
            sequence_index,
            target_name,
            target_id,
            subset_target_start,
            subset_target_end,
            is_reverse,
        );

        // Adjust contig offsets for the subset
        let adjusted_query_offset = alignment.query_contig_start as usize
            + (subset_query_start - metadata.query_start) as usize;
        let adjusted_target_offset = alignment.target_contig_start as usize
            + (subset_target_start - metadata.target_start) as usize;

        // Reconstruct CIGAR for subset only
        let trace_spacing = alignment.trace_spacing as usize;
        let cigar_str = with_edit_aligner(|aligner| {
            tracepoints_to_cigar_fastga_with_aligner(
                &subset_tracepoints,
                trace_spacing.try_into().unwrap(),
                &query_seq,
                &target_seq,
                adjusted_query_offset,
                adjusted_target_offset,
                is_reverse,
                aligner,
            )
        });

        match parse_cigar_to_delta(&cigar_str) {
            Ok(ops) => ops,
            Err(e) => panic!(
                "Failed to parse CIGAR string '{cigar_str}' during subset tracepoint processing: {:?}",
                e
            ),
        }
    }

    fn project_overlapping_interval(
        &self,
        metadata: &QueryMetadata,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        sequence_index: Option<&UnifiedSequenceIndex>,
        min_gap_compressed_identity: Option<f64>,
    ) -> Option<(Interval<u32>, Vec<CigarOp>, Interval<u32>)> {
        // Check if this is a .1aln file that supports tracepoint subsetting
        let file_index = metadata.alignment_file_index as usize;
        let alignment_file = &self.alignment_files[file_index];
        let is_onealn = alignment_file.ends_with(".1aln");

        // Try subsetting approach for .1aln files
        if let (true, Some(sequence_index)) = (is_onealn, sequence_index) {
            // Fetch alignment with tracepoints
            let alignment = self.get_onealn_alignment(metadata).unwrap_or_else(|e| {
                panic!("Subsetting failed: cannot fetch .1aln alignment: {}", e)
            });

            let is_reversed = metadata.is_reversed();

            // Scan tracepoints to find overlapping segments
            // For reversed entries, scan by query coords; for original, scan by target coords
            let (subset_first_idx, subset_last_idx, subset_first_pos, subset_last_pos,
                 num_overlapping_segments, total_matches, total_mismatches,
                 first_segment_info, last_segment_info) = if is_reversed {
                let subset = self.scan_overlapping_tracepoints_by_query(&alignment, metadata, range_start, range_end)
                    .unwrap_or_else(|| panic!(
                        "Subsetting failed: no overlapping segments found for query range {}-{} (target_id={}, file_index={})",
                        range_start, range_end, target_id, metadata.alignment_file_index
                    ));
                // For reversed: positions are target coords, segment info is target-centric
                (subset.first_idx, subset.last_idx, subset.first_target_pos, subset.last_target_pos,
                 subset.num_overlapping_segments, subset.total_matches, subset.total_mismatches,
                 subset.first_segment_info, subset.last_segment_info)
            } else {
                let subset = self.scan_overlapping_tracepoints(&alignment, metadata, range_start, range_end)
                    .unwrap_or_else(|| panic!(
                        "Subsetting failed: no overlapping segments found for target range {}-{} (target_id={}, file_index={})",
                        range_start, range_end, target_id, metadata.alignment_file_index
                    ));
                // For original: positions are query coords, segment info is query-centric
                (subset.first_idx, subset.last_idx, subset.first_query_pos, subset.last_query_pos,
                 subset.num_overlapping_segments, subset.total_matches, subset.total_mismatches,
                 subset.first_segment_info, subset.last_segment_info)
            };

            // Create a temporary SubsettingResult for process_subset_tracepoints
            // (it expects the original format regardless of reversed status)
            let subset_for_cigar = SubsettingResult {
                first_idx: subset_first_idx,
                last_idx: subset_last_idx,
                first_query_pos: if is_reversed { metadata.query_start } else { subset_first_pos },
                last_query_pos: if is_reversed { metadata.query_end } else { subset_last_pos },
                num_overlapping_segments,
                total_matches,
                total_mismatches,
                first_segment_info: if is_reversed {
                    // Convert from inverse format to original format for CIGAR processing
                    // Original: (query_pos, query_delta, seg_start, seg_end, abs_target_delta, trace_diffs)
                    // Inverse:  (target_pos, target_delta, seg_query_start, seg_query_end, query_delta, trace_diffs)
                    let (target_pos, target_delta, seg_query_start, _seg_query_end, query_delta, trace_diffs) = first_segment_info;
                    (seg_query_start, query_delta, target_pos.min(target_pos + target_delta), target_pos.max(target_pos + target_delta), target_delta.abs(), trace_diffs)
                } else {
                    first_segment_info
                },
                last_segment_info: if is_reversed {
                    let (target_pos, target_delta, seg_query_start, _seg_query_end, query_delta, trace_diffs) = last_segment_info;
                    (seg_query_start, query_delta, target_pos.min(target_pos + target_delta), target_pos.max(target_pos + target_delta), target_delta.abs(), trace_diffs)
                } else {
                    last_segment_info
                },
            };

            // Reconstruct CIGAR from subset tracepoints
            let cigar_ops = self.process_subset_tracepoints(
                &alignment,
                metadata,
                target_id,
                &subset_for_cigar,
                sequence_index,
            );

            // Check identity threshold if specified
            if let Some(threshold) = min_gap_compressed_identity {
                if calculate_gap_compressed_identity(&cigar_ops) < threshold {
                    return None;
                }
            }

            // Calculate alignment boundaries for CIGAR projection
            let (alignment_query_start, alignment_query_end, alignment_target_start, alignment_target_end) = if is_reversed {
                // For reversed: we scanned by query, got target positions
                let (_, _, first_seg_query_start, first_seg_query_end, _, _) = first_segment_info;
                let (_, _, last_seg_query_start, last_seg_query_end, _, _) = last_segment_info;
                let query_start = first_seg_query_start.min(first_seg_query_end).min(last_seg_query_start.min(last_seg_query_end));
                let query_end = if subset_last_idx == alignment.tracepoints.len() - 1 {
                    metadata.query_end
                } else {
                    first_seg_query_start.max(first_seg_query_end).max(last_seg_query_start.max(last_seg_query_end))
                };
                let target_start = subset_first_pos;
                let target_end = if subset_last_idx == alignment.tracepoints.len() - 1 {
                    metadata.target_end
                } else {
                    subset_last_pos
                };
                (query_start, query_end, target_start, target_end)
            } else {
                // For original: we scanned by target, got query positions
                let subset_query_start = subset_first_pos;
                let subset_query_end = if subset_last_idx == alignment.tracepoints.len() - 1 {
                    metadata.query_end
                } else {
                    subset_last_pos
                };
                let (_, _, first_seg_start, first_seg_end, _, _) = first_segment_info;
                let (_, _, last_seg_start, last_seg_end, _, _) = last_segment_info;
                let subset_target_start = first_seg_start.min(first_seg_end).min(last_seg_start.min(last_seg_end));
                let subset_target_end = if subset_last_idx == alignment.tracepoints.len() - 1 {
                    metadata.target_end
                } else {
                    first_seg_start.max(first_seg_end).max(last_seg_start.max(last_seg_end))
                };
                (subset_query_start, subset_query_end, subset_target_start, subset_target_end)
            };

            // Project the requested range through the CIGAR
            let projection = if is_reversed {
                project_query_range_through_alignment(
                    (range_start, range_end),
                    (alignment_target_start, alignment_target_end, alignment_query_start, alignment_query_end, metadata.strand()),
                    &cigar_ops,
                )
            } else {
                project_target_range_through_alignment(
                    (range_start, range_end),
                    (alignment_target_start, alignment_target_end, alignment_query_start, alignment_query_end, metadata.strand()),
                    &cigar_ops,
                )
            }.unwrap_or_else(|| panic!(
                "Subsetting failed: projection failed for CIGAR (range={}-{}, query={}-{}, target={}-{}, reversed={})",
                range_start, range_end, alignment_query_start, alignment_query_end, alignment_target_start, alignment_target_end, is_reversed
            ));

            let (
                adjusted_query_start,
                adjusted_query_end,
                adjusted_cigar,
                adjusted_target_start,
                adjusted_target_end,
            ) = projection;

            let query_interval = Interval {
                first: adjusted_query_start,
                last: adjusted_query_end,
                metadata: metadata.query_id,
            };

            let target_interval = Interval {
                first: adjusted_target_start,
                last: adjusted_target_end,
                metadata: target_id,
            };

            return Some((query_interval, adjusted_cigar, target_interval));
        }

        // Fallback: Full CIGAR reconstruction + projection (for non-.1aln files or if subsetting failed)
        let mut cigar_ops = self.get_cigar_ops(metadata, target_id, sequence_index);

        // Dispatch projection based on whether this is a reversed (query-indexed) entry
        let projection = if metadata.is_reversed() {
            // Inverse projection: query coords → target coords
            project_query_range_through_alignment(
                (range_start, range_end),
                (
                    metadata.target_start,
                    metadata.target_end,
                    metadata.query_start,
                    metadata.query_end,
                    metadata.strand(),
                ),
                &cigar_ops,
            )
        } else {
            // Forward projection: target coords → query coords (existing behavior)
            project_target_range_through_alignment(
                (range_start, range_end),
                (
                    metadata.target_start,
                    metadata.target_end,
                    metadata.query_start,
                    metadata.query_end,
                    metadata.strand(),
                ),
                &cigar_ops,
            )
        };

        if let Some((
            adjusted_query_start,
            adjusted_query_end,
            adjusted_cigar,
            adjusted_target_start,
            adjusted_target_end,
        )) = projection
        {
            if let Some(threshold) = min_gap_compressed_identity {
                if calculate_gap_compressed_identity(&adjusted_cigar) < threshold {
                    return None;
                }
            }

            let query_interval = Interval {
                first: adjusted_query_start,
                last: adjusted_query_end,
                metadata: metadata.query_id,
            };

            let target_interval = Interval {
                first: adjusted_target_start,
                last: adjusted_target_end,
                metadata: target_id,
            };

            cigar_ops = adjusted_cigar;

            Some((query_interval, cigar_ops, target_interval))
        } else {
            debug!(
                "Projection failed for alignment (target_id={target_id}, file_index={}, range={}-{})",
                metadata.alignment_file_index,
                range_start,
                range_end
            );
            None
        }
    }

    /// Fast approximate projection using tracepoints without CIGAR computation or sequence fetching.
    /// Returns approximate intervals with identity metrics calculated from tracepoint statistics.
    fn project_overlapping_interval_fast(
        &self,
        metadata: &QueryMetadata,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        min_gap_compressed_identity: Option<f64>,
    ) -> Option<(Interval<u32>, Vec<CigarOp>, Interval<u32>)> {
        let is_reversed = metadata.is_reversed();
        let is_reverse_strand = metadata.strand() == Strand::Reverse;

        // Fetch tracepoints without sequence I/O
        let alignment = match self.get_onealn_alignment(metadata) {
            Ok(aln) => aln,
            Err(e) => {
                warn!("Failed to fetch alignment for fast projection: {}", e);
                return None;
            }
        };

        if is_reversed {
            // INVERSE PATH: query coords → target coords
            // For reversed entries, range_start/range_end are in query coordinate space
            // and we need to find the corresponding target coordinates
            let subset = self.scan_overlapping_tracepoints_by_query(
                &alignment,
                metadata,
                range_start,
                range_end,
            )?;

            let working_target_start = metadata.target_start;
            let working_target_end = metadata.target_end;

            // Refine target coordinates using indel-aware heuristic (inverse of normal refinement)
            let refine_target_boundary = |target_pos: i32,
                                          target_delta: i32,
                                          segment_query_start: i32,
                                          overlap_pos: i32,
                                          abs_query_delta: i32,
                                          trace_diffs: i32,
                                          num_segs: usize|
             -> i32 {
                let aligned_len = target_delta.abs().min(abs_query_delta);
                let segment_identity = if aligned_len > 0 {
                    ((aligned_len - trace_diffs) as f64 / aligned_len as f64).max(0.0)
                } else {
                    0.0
                };

                if abs_query_delta == 0 {
                    panic!(
                        "Cannot refine target position with zero query delta: num_segs={}, identity={:.3}, trace_diffs={}",
                        num_segs, segment_identity, trace_diffs
                    );
                }

                let query_fraction =
                    (overlap_pos - segment_query_start) as f64 / abs_query_delta as f64;
                let indel_ratio = target_delta.abs() as f64 / abs_query_delta as f64;
                let target_advance = query_fraction * abs_query_delta as f64 * indel_ratio;
                let refined_pos = target_pos + target_advance.round() as i32;

                refined_pos.max(working_target_start).min(working_target_end)
            };

            // Refine first target position
            // segment_info: (target_pos, target_delta, seg_query_start, seg_query_end, abs_query_delta, trace_diffs)
            let (target_pos, target_delta, seg_query_start, _seg_query_end, abs_query_delta, trace_diffs) =
                subset.first_segment_info;
            let overlap_start = seg_query_start.max(range_start);
            let refined_first = refine_target_boundary(
                target_pos,
                target_delta,
                seg_query_start,
                overlap_start,
                abs_query_delta,
                trace_diffs,
                subset.num_overlapping_segments,
            );

            // Refine last target position
            let (target_pos, target_delta, seg_query_start, seg_query_end, abs_query_delta, trace_diffs) =
                subset.last_segment_info;
            let overlap_end = seg_query_end.min(range_end);
            let refined_last = refine_target_boundary(
                target_pos,
                target_delta,
                seg_query_start,
                overlap_end,
                abs_query_delta,
                trace_diffs,
                subset.num_overlapping_segments,
            );

            // Create approximate CIGAR from accumulated statistics
            let mut approx_cigar = Vec::new();
            if subset.total_matches > 0.0 {
                approx_cigar.push(CigarOp::new(subset.total_matches.round() as i32, '='));
            }
            if subset.total_mismatches > 0.0 {
                approx_cigar.push(CigarOp::new(subset.total_mismatches.round() as i32, 'X'));
            }

            // Check identity threshold if specified
            if let Some(threshold) = min_gap_compressed_identity {
                let identity = calculate_gap_compressed_identity(&approx_cigar);
                if identity < threshold {
                    return None;
                }
            }

            // For reversed entries with reverse strand: swap target coordinates
            let (target_start, target_end) = if is_reverse_strand {
                (refined_last, refined_first)
            } else {
                (refined_first, refined_last)
            };

            // Validate coordinates are non-negative
            if refined_first < 0 || refined_last < 0 {
                panic!(
                    "Inverse projection resulted in negative target coordinates: {}-{}",
                    refined_first, refined_last
                );
            }

            // For reversed entries: query_interval contains the projected target coords,
            // target_interval contains the input query range
            let query_interval = Interval {
                first: target_start,
                last: target_end,
                metadata: metadata.query_id, // This is the original target_id for reversed entries
            };

            let target_interval = Interval {
                first: range_start,
                last: range_end,
                metadata: target_id, // This is the original query_id for reversed entries
            };

            Some((query_interval, approx_cigar, target_interval))
        } else {
            // ORIGINAL PATH: target coords → query coords
            let subset =
                self.scan_overlapping_tracepoints(&alignment, metadata, range_start, range_end)?;

            let working_query_start = metadata.query_start;
            let working_query_end = metadata.query_end;

            // Refine query coordinates using indel-aware heuristic
            let refine_query_boundary = |query_pos: i32,
                                         query_delta: i32,
                                         segment_target_start: i32,
                                         overlap_pos: i32,
                                         abs_target_delta: i32,
                                         trace_diffs: i32,
                                         num_segs: usize|
             -> i32 {
                let aligned_len = query_delta.min(abs_target_delta);
                let segment_identity = if aligned_len > 0 {
                    ((aligned_len - trace_diffs) as f64 / aligned_len as f64).max(0.0)
                } else {
                    0.0
                };

                if abs_target_delta == 0 {
                    panic!(
                        "Cannot refine query position with zero target delta: num_segs={}, identity={:.3}, trace_diffs={}",
                        num_segs, segment_identity, trace_diffs
                    );
                }

                let target_fraction =
                    (overlap_pos - segment_target_start) as f64 / abs_target_delta as f64;
                let indel_ratio = query_delta as f64 / abs_target_delta as f64;
                let query_advance = target_fraction * abs_target_delta as f64 * indel_ratio;
                let refined_pos = query_pos + query_advance.round() as i32;

                refined_pos.max(working_query_start).min(working_query_end)
            };

            // Refine first and last query positions
            let (query_pos, query_delta, segment_target_start, _, abs_target_delta, trace_diffs) =
                subset.first_segment_info;
            let overlap_start = segment_target_start.max(range_start);
            let refined_first = refine_query_boundary(
                query_pos,
                query_delta,
                segment_target_start,
                overlap_start,
                abs_target_delta,
                trace_diffs,
                subset.num_overlapping_segments,
            );

            let (
                query_pos,
                query_delta,
                segment_target_start,
                segment_target_end,
                abs_target_delta,
                trace_diffs,
            ) = subset.last_segment_info;
            let overlap_end = segment_target_end.min(range_end);
            let refined_last = refine_query_boundary(
                query_pos,
                query_delta,
                segment_target_start,
                overlap_end,
                abs_target_delta,
                trace_diffs,
                subset.num_overlapping_segments,
            );

            // Create approximate CIGAR from accumulated statistics (for identity calculation only)
            let mut approx_cigar = Vec::new();
            if subset.total_matches > 0.0 {
                approx_cigar.push(CigarOp::new(subset.total_matches.round() as i32, '='));
            }
            if subset.total_mismatches > 0.0 {
                approx_cigar.push(CigarOp::new(subset.total_mismatches.round() as i32, 'X'));
            }

            // Check identity threshold if specified
            if let Some(threshold) = min_gap_compressed_identity {
                let identity = calculate_gap_compressed_identity(&approx_cigar);
                if identity < threshold {
                    return None;
                }
            }

            // Get refined query coordinates
            let working_query_start = refined_first;
            let working_query_end = refined_last;

            // For reverse alignments: swap projected query coordinates to match normal mode output
            let (query_start, query_end) = if is_reverse_strand {
                (working_query_end, working_query_start) // Swap: first=end, last=start
            } else {
                (working_query_start, working_query_end)
            };

            // Target coordinates: use exact requested range
            let target_start = range_start;
            let target_end = range_end;

            // Validate coordinates are non-negative
            if working_query_start < 0 || working_query_end < 0 {
                panic!(
                    "Projection resulted in negative query coordinates: {}-{}",
                    working_query_start, working_query_end
                );
            }

            let query_interval = Interval {
                first: query_start,
                last: query_end,
                metadata: metadata.query_id,
            };

            let target_interval = Interval {
                first: target_start,
                last: target_end,
                metadata: target_id,
            };

            Some((query_interval, approx_cigar, target_interval))
        }
    }

    pub fn from_multi_alignment_records(
        records_by_file: &[(Vec<AlignmentRecord>, String)],
        seq_index: SequenceIndex,
        sequence_files: Option<&[String]>,
    ) -> io::Result<Self> {
        // Extract just the alignment file paths
        let alignment_files: Vec<String> = records_by_file
            .par_iter()
            .map(|(_, alignment_file)| alignment_file.clone())
            .collect();

        let intervals: FxHashMap<u32, Vec<Interval<QueryMetadata>>> = records_by_file
            .par_iter()
            .enumerate() // Add enumeration to get the position as index
            .flat_map(|(file_index, (records, _))| {
                records
                    .par_iter()
                    .flat_map(|record| {
                        // Create TWO entries per alignment for bidirectional indexing:
                        // 1. Original entry: keyed by target, stores query in metadata
                        // 2. Reversed entry: keyed by query, stores target in metadata (swapped roles)

                        let mut entries = Vec::with_capacity(2);

                        // 1. ORIGINAL entry (existing behavior)
                        let original_metadata = QueryMetadata {
                            query_id: record.query_id,
                            target_start: record.target_start as i32,
                            target_end: record.target_end as i32,
                            query_start: record.query_start as i32,
                            query_end: record.query_end as i32,
                            alignment_file_index: file_index as u32,
                            strand_and_data_offset: record.strand_and_data_offset, // Already includes strand bit
                            data_bytes: record.data_bytes,
                        };
                        entries.push((
                            record.target_id,
                            Interval {
                                first: record.target_start as i32,
                                last: record.target_end as i32,
                                metadata: original_metadata,
                            },
                        ));

                        // 2. REVERSED entry (new for bidirectional)
                        // Keyed by query_id with query coords as interval bounds
                        // Metadata stores target_id as "query_id" (the sequence we project TO)
                        // REVERSED_BIT is set to indicate inverse projection is needed
                        let reversed_metadata = QueryMetadata {
                            query_id: record.target_id,  // Swap: original target becomes "query" in metadata
                            target_start: record.target_start as i32,  // Keep original coords
                            target_end: record.target_end as i32,
                            query_start: record.query_start as i32,
                            query_end: record.query_end as i32,
                            alignment_file_index: file_index as u32,
                            strand_and_data_offset: record.strand_and_data_offset | QueryMetadata::REVERSED_BIT,
                            data_bytes: record.data_bytes,
                        };
                        entries.push((
                            record.query_id,  // Key by original QUERY sequence
                            Interval {
                                first: record.query_start as i32,  // Use QUERY coords for interval
                                last: record.query_end as i32,
                                metadata: reversed_metadata,
                            },
                        ));

                        entries
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

        let num_files = alignment_files.len();
        Ok(Self {
            trees: RwLock::new(trees),
            seq_index,
            alignment_files,
            forest_map,
            index_file_path: String::new(),
            sequence_files: sequence_files.map(|s| s.to_vec()).unwrap_or_default(),
            trace_spacing_cache: RwLock::new(vec![None; num_files]),
        })
    }

    /// Serialize the IMPG index to a writer with embedded forest map
    pub fn serialize_with_forest_map<W: std::io::Write + std::io::Seek>(
        &self,
        mut writer: W,
    ) -> std::io::Result<()> {
        const MAGIC: &[u8] = b"IMPGIDX2"; // Version 2: bidirectional indexing
        let mut forest_map = ForestMap::new();

        // Write magic bytes
        writer.write_all(MAGIC)?;

        // Write placeholder for forest map offset (will be updated later)
        writer.write_all(&[0u8; 8])?;

        let mut current_offset = 16u64; // Header size

        // Serialize sequence index
        let seq_index_data =
            bincode::serde::encode_to_vec(&self.seq_index, bincode::config::standard()).map_err(
                |e| std::io::Error::other(format!("Failed to encode sequence index: {e}")),
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
                            "Failed to encode tree for target {target_id}: {e}"
                        ))
                    })?;

            writer.write_all(&tree_data)?;
            current_offset += tree_data.len() as u64;
        }

        // Now write the forest map at the end
        let forest_map_offset = current_offset;
        let forest_map_data =
            bincode::serde::encode_to_vec(&forest_map, bincode::config::standard())
                .map_err(|e| std::io::Error::other(format!("Failed to encode forest map: {e}")))?;

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
        sequence_files: Option<&[String]>,
    ) -> std::io::Result<Self> {
        const MAGIC: &[u8] = b"IMPGIDX2"; // Version 2: bidirectional indexing

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
                        format!("Failed to load sequence index: {e}"),
                    )
                })?;

        // Seek to forest map and read it
        reader.seek(SeekFrom::Start(forest_map_offset))?;
        let forest_map: ForestMap =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load forest map: {e}"),
                    )
                })?;

        let num_files = alignment_files.len();
        Ok(Self {
            trees: RwLock::new(FxHashMap::default()),
            seq_index,
            alignment_files: alignment_files.to_vec(),
            forest_map,
            index_file_path,
            sequence_files: sequence_files.map(|s| s.to_vec()).unwrap_or_default(),
            trace_spacing_cache: RwLock::new(vec![None; num_files]),
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
        approximate_mode: bool,
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
            "Querying region{}: {}:{}-{}, len: {}",
            if approximate_mode {
                " (approximate)"
            } else {
                ""
            },
            self.seq_index.get_name(target_id).unwrap(),
            range_start,
            range_end,
            range_end - range_start
        );

        // Get or load the tree - if None, no overlaps exist for this target
        if let Some(tree) = self.get_or_load_tree(target_id) {
            tree.query(range_start, range_end, |interval| {
                let metadata = &interval.metadata;
                let projection = if approximate_mode {
                    // Approximate mode: fast projection without sequence I/O
                    self.project_overlapping_interval_fast(
                        metadata,
                        target_id,
                        range_start,
                        range_end,
                        min_gap_compressed_identity,
                    )
                } else {
                    // Normal mode: full CIGAR computation with sequences
                    self.project_overlapping_interval(
                        metadata,
                        target_id,
                        range_start,
                        range_end,
                        sequence_index,
                        min_gap_compressed_identity,
                    )
                };

                if let Some((query_interval, cigar_ops, target_interval)) = projection {
                    let cigar_vec = if store_cigar { cigar_ops } else { Vec::new() };
                    results.push((query_interval, cigar_vec, target_interval));
                }
            });
        }

        results
    }

    pub fn populate_cigar_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        _min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cache: &mut FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) {
        if let Some(tree) = self.get_or_load_tree(target_id) {
            tree.query(range_start, range_end, |interval| {
                let metadata = &interval.metadata;
                let cache_key = (metadata.alignment_file_index, metadata.data_offset());
                cache
                    .entry(cache_key)
                    .or_insert_with(|| self.get_cigar_ops(metadata, target_id, sequence_index));
            });
        }
    }

    pub fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) -> Vec<AdjustedInterval> {
        let mut results = Vec::new();
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

        if let Some(tree) = self.get_or_load_tree(target_id) {
            tree.query(range_start, range_end, |interval| {
                let metadata = &interval.metadata;
                let cache_key = (metadata.alignment_file_index, metadata.data_offset());
                let cigar_ops = cigar_cache
                    .get(&cache_key)
                    .cloned()
                    .unwrap_or_else(|| self.get_cigar_ops(metadata, target_id, sequence_index));

                let result = project_target_range_through_alignment(
                    (range_start, range_end),
                    (
                        metadata.target_start,
                        metadata.target_end,
                        metadata.query_start,
                        metadata.query_end,
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
                    if let Some(threshold) = min_gap_compressed_identity {
                        if calculate_gap_compressed_identity(&adjusted_cigar) < threshold {
                            return;
                        }
                    }

                    results.push((
                        Interval {
                            first: adjusted_query_start,
                            last: adjusted_query_end,
                            metadata: metadata.query_id,
                        },
                        if store_cigar {
                            adjusted_cigar
                        } else {
                            Vec::new()
                        },
                        Interval {
                            first: adjusted_target_start,
                            last: adjusted_target_end,
                            metadata: target_id,
                        },
                    ));
                }
            });
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
        min_output_length: Option<i32>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&crate::subset_filter::SubsetFilter>,
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

            //let prec_num_results = results.len();

            // Get or load the tree - if None, no overlaps exist for this target
            if let Some(tree) = self.get_or_load_tree(current_target_id) {
                let mut intervals: Vec<(QueryMetadata, (i32, i32))> = Vec::new();
                tree.query(current_target_start, current_target_end, |interval| {
                    let overlap_start = current_target_start.max(interval.first);
                    let overlap_end = current_target_end.min(interval.last);
                    if overlap_start < overlap_end {
                        intervals.push((interval.metadata.clone(), (overlap_start, overlap_end)));
                    }
                });

                let processed_results: Vec<_> = intervals
                    .into_par_iter()
                    .filter_map(|(metadata, overlap_range)| {
                        let projection_result = if approximate_mode {
                            self.project_overlapping_interval_fast(
                                &metadata,
                                current_target_id,
                                overlap_range.0,
                                overlap_range.1,
                                min_gap_compressed_identity,
                            )
                        } else {
                            self.project_overlapping_interval(
                                &metadata,
                                current_target_id,
                                overlap_range.0,
                                overlap_range.1,
                                sequence_index,
                                min_gap_compressed_identity,
                            )
                        };
                        projection_result.and_then(|(query_interval, cigar_ops, target_interval)| {
                            let query_id = query_interval.metadata;

                            // Apply subset filter: keep if it's the target or matches the filter
                            let should_keep = if let Some(filter) = subset_filter {
                                query_id == target_id ||
                                self.seq_index
                                    .get_name(query_id)
                                    .is_some_and(|name| filter.matches(name))
                            } else {
                                true
                            };

                            if should_keep {
                                let adjusted_query_start = query_interval.first;
                                let adjusted_query_end = query_interval.last;
                                let cigar_vec = if store_cigar { cigar_ops } else { Vec::new() };

                                Some((
                                    query_interval,
                                    cigar_vec,
                                    target_interval,
                                    query_id,
                                    adjusted_query_start,
                                    adjusted_query_end,
                                ))
                            } else {
                                None
                            }
                        })
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
                    let length = (query_interval.last - query_interval.first).abs();

                    // Add to results only if it passes min_output_length filter
                    let should_add_to_output = if let Some(min_len) = min_output_length {
                        length >= min_len
                    } else {
                        true
                    };
                    if should_add_to_output {
                        results.push((query_interval, cigar.clone(), target_interval));
                    }

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

            //debug!("Collected {} results", results.len() - prec_num_results);

            // Merge contiguous/overlapping ranges with same sequence_id
            //let stack_size = stack.len();
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
            //debug!("Merged stack size from {} to {}", stack_size, stack.len());
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
        min_output_length: Option<i32>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&crate::subset_filter::SubsetFilter>,
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

        // Process by depth until max_depth or no more ranges
        while !current_ranges.is_empty() && (max_depth == 0 || current_depth < max_depth) {
            // debug!(
            //     "Processing depth {} with {} ranges",
            //     current_depth,
            //     current_ranges.len()
            // );

            // Process current depth ranges in parallel
            let query_results: Vec<Vec<(u32, i32, i32, Vec<CigarOp>, i32, i32, u32)>> =
                current_ranges
                    .par_iter()
                    .map(
                        |(current_target_id, current_target_start, current_target_end)| {
                            let mut local_results = Vec::new();

                            // Get or load the tree - if None, no overlaps exist for this target
                            if let Some(tree) = self.get_or_load_tree(*current_target_id) {
                                tree.query(
                                    *current_target_start,
                                    *current_target_end,
                                    |interval| {
                                        let metadata = &interval.metadata;
                                        let overlap_start =
                                            (*current_target_start).max(interval.first);
                                        let overlap_end = (*current_target_end).min(interval.last);
                                        if overlap_start >= overlap_end {
                                            return;
                                        }

                                        let projection_result = if approximate_mode {
                                            self.project_overlapping_interval_fast(
                                                metadata,
                                                *current_target_id,
                                                overlap_start,
                                                overlap_end,
                                                min_gap_compressed_identity,
                                            )
                                        } else {
                                            self.project_overlapping_interval(
                                                metadata,
                                                *current_target_id,
                                                overlap_start,
                                                overlap_end,
                                                sequence_index,
                                                min_gap_compressed_identity,
                                            )
                                        };

                                        if let Some((query_interval, cigar_ops, target_interval)) =
                                            projection_result
                                        {
                                            let query_id = query_interval.metadata;

                                            // Apply subset filter: keep if it's the target or matches the filter
                                            let should_keep = if let Some(filter) = subset_filter {
                                                query_id == target_id ||
                                                self.seq_index
                                                    .get_name(query_id)
                                                    .is_some_and(|name| filter.matches(name))
                                            } else {
                                                true
                                            };

                                            if should_keep {
                                                let cigar_vec =
                                                    if store_cigar { cigar_ops } else { Vec::new() };

                                                local_results.push((
                                                    query_id,
                                                    query_interval.first,
                                                    query_interval.last,
                                                    cigar_vec,
                                                    target_interval.first,
                                                    target_interval.last,
                                                    *current_target_id,
                                                ));
                                            }
                                        }
                                    },
                                );
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
                    let length = (adjusted_query_end - adjusted_query_start).abs();

                    // Add to results only if it passes min_output_length filter
                    let should_add_to_output = if let Some(min_len) = min_output_length {
                        length >= min_len
                    } else {
                        true
                    };
                    if should_add_to_output {
                        results.push((
                            Interval {
                                first: adjusted_query_start,
                                last: adjusted_query_end,
                                metadata: query_id,
                            },
                            adjusted_cigar.clone(),
                            Interval {
                                first: adjusted_target_start,
                                last: adjusted_target_end,
                                metadata: current_target_id,
                            },
                        ));
                    }

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

                // debug!(
                //     "Next depth will process {} ranges after merging",
                //     next_depth_ranges.len()
                // );
            }

            // Set up for next iteration
            current_ranges = next_depth_ranges;
        }

        results
    }
}

// Implement ImpgIndex trait for Impg
impl ImpgIndex for Impg {
    fn seq_index(&self) -> &SequenceIndex {
        &self.seq_index
    }

    fn query(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
    ) -> Vec<AdjustedInterval> {
        // Delegate to the existing method
        Impg::query(
            self,
            target_id,
            range_start,
            range_end,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            approximate_mode,
        )
    }

    fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) -> Vec<AdjustedInterval> {
        Impg::query_with_cache(
            self,
            target_id,
            range_start,
            range_end,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            cigar_cache,
        )
    }

    fn populate_cigar_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cache: &mut FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) {
        Impg::populate_cigar_cache(
            self,
            target_id,
            range_start,
            range_end,
            min_gap_compressed_identity,
            sequence_index,
            cache,
        )
    }

    fn query_transitive_dfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        min_output_length: Option<i32>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> Vec<AdjustedInterval> {
        Impg::query_transitive_dfs(
            self,
            target_id,
            range_start,
            range_end,
            masked_regions,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            min_output_length,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            approximate_mode,
            subset_filter,
        )
    }

    fn query_transitive_bfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        min_output_length: Option<i32>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> Vec<AdjustedInterval> {
        Impg::query_transitive_bfs(
            self,
            target_id,
            range_start,
            range_end,
            masked_regions,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            min_output_length,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            approximate_mode,
            subset_filter,
        )
    }

    fn get_or_load_tree(
        &self,
        target_id: u32,
    ) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        Impg::get_or_load_tree(self, target_id)
    }

    fn target_ids(&self) -> Vec<u32> {
        self.forest_map.entries.keys().copied().collect()
    }

    fn remove_cached_tree(&self, target_id: u32) {
        self.trees.write().unwrap().remove(&target_id);
    }

    fn num_targets(&self) -> usize {
        self.forest_map.entries.len()
    }

    fn sequence_files(&self) -> &[String] {
        &self.sequence_files
    }
}

fn project_target_range_through_alignment(
    requested_target_range: (i32, i32),
    record: (i32, i32, i32, i32, Strand),
    cigar_ops: &[CigarOp],
) -> Option<(i32, i32, Vec<CigarOp>, i32, i32)> {
    let (target_start, target_end, query_start, query_end, strand) = record;
    debug!(
        "Projecting target range {}-{} through alignment with target {}-{}, query {}-{}, strand {:?}",
        requested_target_range.0,
        requested_target_range.1,
        target_start,
        target_end,
        query_start,
        query_end,
        strand
    );

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

/// Inverse of project_target_range_through_alignment.
/// Given a range in query coordinates, project it back to target coordinates.
/// Used for reversed index entries where we query by original query coords.
///
/// Returns: (target_start, target_end, cigar_subset, query_start, query_end)
/// where target coords are the projected result and query coords are the overlapping input range.
fn project_query_range_through_alignment(
    requested_query_range: (i32, i32),
    record: (i32, i32, i32, i32, Strand),
    cigar_ops: &[CigarOp],
) -> Option<(i32, i32, Vec<CigarOp>, i32, i32)> {
    let (target_start, _target_end, query_start, query_end, strand) = record;
    debug!(
        "Inverse projecting query range {}-{} through alignment with query {}-{}, strand {:?}",
        requested_query_range.0,
        requested_query_range.1,
        query_start,
        query_end,
        strand
    );

    // For CIGAR traversal:
    // - Target always advances forward from target_start
    // - Query advances forward for Forward strand, backward for Reverse strand
    let dir = if strand == Strand::Forward { 1 } else { -1 };
    let mut query_pos = if strand == Strand::Forward {
        query_start
    } else {
        query_end
    };
    let mut target_pos = target_start; // Target always starts at target_start for CIGAR

    // Track CIGAR slice bounds
    let mut first_op_idx = 0;
    let mut last_op_idx = 0;
    let mut found_overlap = false;

    let mut projected_target_start = -1;
    let mut projected_target_end = -1;
    let mut projected_query_start = -1;
    let mut projected_query_end = -1;
    let mut first_op_offset = 0;
    let mut last_op_remaining = 0;

    // For query-based overlap checking, we need to normalize coordinates for reverse strand
    // In forward strand: query coords go low->high as we traverse CIGAR
    // In reverse strand: query coords go high->low as we traverse CIGAR
    let (query_range_low, query_range_high) = if strand == Strand::Forward {
        (requested_query_range.0, requested_query_range.1)
    } else {
        // For reverse strand, the requested range might be given as (high, low) or (low, high)
        // Normalize to (low, high) for comparison
        (
            requested_query_range.0.min(requested_query_range.1),
            requested_query_range.0.max(requested_query_range.1),
        )
    };

    for (curr_op_idx, cigar_op) in cigar_ops.iter().enumerate() {
        let target_delta = cigar_op.target_delta();
        let query_delta = cigar_op.query_delta(strand);

        match (target_delta, query_delta) {
            (0, query_delta) => {
                // Insertion in query (no target movement)
                // This contributes query bases but doesn't move target
                let (q_start, q_end) = if dir > 0 {
                    (query_pos, query_pos + query_delta)
                } else {
                    (query_pos + query_delta, query_pos) // query_delta is negative for reverse
                };

                // Check overlap with requested query range
                let overlap_start = q_start.max(query_range_low);
                let overlap_end = q_end.min(query_range_high);

                if overlap_start < overlap_end {
                    if !found_overlap {
                        projected_target_start = target_pos;
                        projected_query_start = overlap_start;
                        first_op_idx = curr_op_idx;
                        first_op_offset = if dir > 0 {
                            overlap_start - q_start
                        } else {
                            q_end - overlap_end
                        };
                        found_overlap = true;
                    }
                    projected_target_end = target_pos; // No target movement
                    projected_query_end = overlap_end;
                    last_op_idx = curr_op_idx + 1;
                    last_op_remaining = if dir > 0 {
                        overlap_end - q_end
                    } else {
                        q_start - overlap_start
                    };
                }
                query_pos += query_delta;
            }
            (target_delta, 0) => {
                // Deletion in query (target moves, query doesn't)
                // This doesn't contribute query bases, so skip if we haven't found overlap yet
                // But if we're in the overlap, accumulate target movement
                if found_overlap {
                    projected_target_end = target_pos + target_delta;
                    last_op_idx = curr_op_idx + 1;
                }
                target_pos += target_delta;
            }
            (target_delta, query_delta) => {
                // Match or mismatch - both move
                let (q_start, q_end) = if dir > 0 {
                    (query_pos, query_pos + query_delta)
                } else {
                    (query_pos + query_delta, query_pos)
                };

                let overlap_start = q_start.max(query_range_low);
                let overlap_end = q_end.min(query_range_high);

                if overlap_start < overlap_end {
                    // Calculate corresponding target positions
                    // For forward strand: query offset directly maps to target offset
                    // For reverse strand: query goes high→low while target goes low→high
                    //   so distance from q_end maps to distance from target_start
                    let (target_overlap_start, target_overlap_end, op_offset) = if dir > 0 {
                        let query_offset_start = overlap_start - q_start;
                        let query_offset_end = overlap_end - q_start;
                        (target_pos + query_offset_start, target_pos + query_offset_end, query_offset_start)
                    } else {
                        // For reverse: overlap_start is the LOW query coord, which corresponds
                        // to the HIGH target coord (further along in the traversal)
                        let dist_from_q_end_start = q_end - overlap_end;   // How far overlap_end is from q_end
                        let dist_from_q_end_end = q_end - overlap_start;   // How far overlap_start is from q_end
                        // For CIGAR offset, use target-based offset since CIGAR follows target
                        (target_pos + dist_from_q_end_start, target_pos + dist_from_q_end_end, dist_from_q_end_start)
                    };

                    if !found_overlap {
                        projected_target_start = target_overlap_start;
                        projected_query_start = overlap_start;
                        first_op_idx = curr_op_idx;
                        first_op_offset = op_offset;
                        found_overlap = true;
                    }
                    projected_target_end = target_overlap_end;
                    projected_query_end = overlap_end;
                    last_op_idx = curr_op_idx + 1;
                    last_op_remaining = overlap_end - q_end; // Will be negative or zero
                }

                target_pos += target_delta;
                query_pos += query_delta;
            }
        }
    }

    // If we had at least one overlap
    if found_overlap
        && projected_target_start != projected_target_end
        && projected_query_start != projected_query_end
    {
        let mut projected_cigar_ops = cigar_ops[first_op_idx..last_op_idx].to_vec();
        // Adjust first operation length
        if first_op_offset > 0 {
            projected_cigar_ops[0].adjust_len(-first_op_offset);
        }
        // Adjust last operation length
        if last_op_remaining < 0 && !projected_cigar_ops.is_empty() {
            projected_cigar_ops[last_op_idx - first_op_idx - 1].adjust_len(last_op_remaining);
        }

        Some((
            projected_target_start,
            projected_target_end,
            projected_cigar_ops,
            projected_query_start,
            projected_query_end,
        ))
    } else {
        None
    }
}

fn parse_cigar_to_delta(cigar: &str) -> Result<Vec<CigarOp>, OneAlnParseErr> {
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

    // Tests for inverse projection (project_query_range_through_alignment)
    #[test]
    fn test_inverse_projection_simple_forward() {
        // Simple 100bp match, forward strand
        let query_range = (0, 100);
        let record = (100, 200, 0, 100, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let result =
            project_query_range_through_alignment(query_range, record, &cigar_ops).unwrap();
        // Should project query 0-100 to target 100-200
        assert_eq!(result, (100, 200, cigar_ops.clone(), 0, 100));
    }

    #[test]
    fn test_inverse_projection_simple_reverse() {
        // Simple 100bp match, reverse strand
        let query_range = (0, 100);
        let record = (100, 200, 0, 100, Strand::Reverse);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let result =
            project_query_range_through_alignment(query_range, record, &cigar_ops).unwrap();
        // For reverse strand, query coords go from 100 down to 0 as we traverse
        // Querying 0-100 should still project to target 100-200
        assert_eq!(result, (100, 200, cigar_ops.clone(), 0, 100));
    }

    #[test]
    fn test_inverse_projection_partial_range() {
        // Partial range query
        let query_range = (10, 50);
        let record = (0, 100, 0, 100, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let result =
            project_query_range_through_alignment(query_range, record, &cigar_ops).unwrap();
        // Query 10-50 should project to target 10-50 (1:1 mapping with match)
        assert_eq!(result, (10, 50, vec![CigarOp::new(40, '=')], 10, 50));
    }

    #[test]
    fn test_inverse_projection_with_insertion() {
        // Alignment with insertion in query
        // Target: 0-50, Query: 0-60 (10bp insertion in query)
        let record = (0, 50, 0, 60, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(25, '='), // Match
            CigarOp::new(10, 'I'), // Insertion in query (query advances, target doesn't)
            CigarOp::new(25, '='), // Match
        ];

        // Query the insertion region (query 25-35)
        // This is pure insertion, so target range should be just the boundary
        let query_range = (25, 35);
        let result = project_query_range_through_alignment(query_range, record, &cigar_ops);
        // The insertion contributes 10 query bases but 0 target bases
        // So we should get target = 25-25 (empty target range... but we filter those out)
        // Actually with our logic, this returns None because target_start == target_end
        assert!(result.is_none());

        // Query a range that spans the insertion (query 20-40)
        let query_range = (20, 40);
        let result =
            project_query_range_through_alignment(query_range, record, &cigar_ops).unwrap();
        // Query 20-25 maps to target 20-25 (5bp)
        // Query 25-35 is insertion (0 target bases)
        // Query 35-40 maps to target 25-30 (5bp)
        // Total: target 20-30 (10bp)
        assert_eq!(result.0, 20); // target_start
        assert_eq!(result.1, 30); // target_end
    }

    #[test]
    fn test_inverse_projection_with_deletion() {
        // Alignment with deletion in query
        // Target: 0-60, Query: 0-50 (10bp deletion in query)
        let record = (0, 60, 0, 50, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(25, '='), // Match
            CigarOp::new(10, 'D'), // Deletion in query (target advances, query doesn't)
            CigarOp::new(25, '='), // Match
        ];

        // Query 20-30 spans the deletion point
        let query_range = (20, 30);
        let result =
            project_query_range_through_alignment(query_range, record, &cigar_ops).unwrap();
        // Query 20-25 maps to target 20-25 (5bp)
        // Deletion: target 25-35 (10bp) with no query contribution
        // Query 25-30 maps to target 35-40 (5bp)
        // Total: target 20-40 (20bp) because deletion is included
        assert_eq!(result.0, 20); // target_start
        assert_eq!(result.1, 40); // target_end (includes the deletion span)
    }

    #[test]
    fn test_inverse_projection_complex_cigar() {
        // Complex CIGAR: 10=5I5D50=
        // Target: 0-65 (10 + 5 + 50 = 65)
        // Query: 0-65 (10 + 5 + 50 = 65)
        let record = (0, 65, 0, 65, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(10, '='),
            CigarOp::new(5, 'I'),  // Query advances 5, target stays
            CigarOp::new(5, 'D'),  // Target advances 5, query stays
            CigarOp::new(50, '='),
        ];

        // Query the full range
        let query_range = (0, 65);
        let result =
            project_query_range_through_alignment(query_range, record, &cigar_ops).unwrap();
        assert_eq!(result.0, 0);  // target_start
        assert_eq!(result.1, 65); // target_end
    }

    // Round-trip property tests: forward projection followed by inverse should return original
    #[test]
    fn test_round_trip_simple_match_forward() {
        // Simple 100bp match, forward strand
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let record = (0, 100, 0, 100, Strand::Forward);

        // Forward: target 20-60 → query
        let target_range = (20, 60);
        let (q_start, q_end, _, t_start_out, t_end_out) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();

        // For a perfect match with same coords, query should equal target
        assert_eq!((q_start, q_end), (20, 60));
        assert_eq!((t_start_out, t_end_out), (20, 60));

        // Inverse: query 20-60 → target (should match)
        let (t_start_back, t_end_back, _, q_start_back, q_end_back) =
            project_query_range_through_alignment((q_start, q_end), record, &cigar_ops).unwrap();

        assert_eq!((t_start_back, t_end_back), (t_start_out, t_end_out));
        assert_eq!((q_start_back, q_end_back), (q_start, q_end));
    }

    #[test]
    fn test_round_trip_simple_match_reverse() {
        // Simple 100bp match, reverse strand
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let record = (0, 100, 0, 100, Strand::Reverse);

        // Forward: target 20-60 → query (in reverse coords)
        let target_range = (20, 60);
        let result = project_target_range_through_alignment(target_range, record, &cigar_ops);
        assert!(result.is_some());
        let (q_start, q_end, _, t_start_out, t_end_out) = result.unwrap();

        // Inverse: query range → target (should get back to original target)
        let inverse_result = project_query_range_through_alignment((q_start, q_end), record, &cigar_ops);
        assert!(inverse_result.is_some());
        let (t_start_back, t_end_back, _, _, _) = inverse_result.unwrap();

        assert_eq!((t_start_back, t_end_back), (t_start_out, t_end_out));
    }

    #[test]
    fn test_round_trip_with_offset() {
        // Match with offset: target 100-200 → query 0-100
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let record = (100, 200, 0, 100, Strand::Forward);

        // Forward: target 120-180 → query
        let target_range = (120, 180);
        let (q_start, q_end, _, t_start_out, t_end_out) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();

        assert_eq!((q_start, q_end), (20, 80)); // Offset by 100
        assert_eq!((t_start_out, t_end_out), (120, 180));

        // Inverse: query 20-80 → target 120-180
        let (t_start_back, t_end_back, _, _, _) =
            project_query_range_through_alignment((q_start, q_end), record, &cigar_ops).unwrap();

        assert_eq!((t_start_back, t_end_back), (t_start_out, t_end_out));
    }

    #[test]
    fn test_round_trip_with_insertion() {
        // CIGAR: 50=10I50= (query is 10bp longer due to insertion)
        // Target: 0-100, Query: 0-110
        let cigar_ops = vec![
            CigarOp::new(50, '='),
            CigarOp::new(10, 'I'), // 10bp insertion in query
            CigarOp::new(50, '='),
        ];
        let record = (0, 100, 0, 110, Strand::Forward);

        // Forward: target 40-60 (spans the insertion point at 50)
        let target_range = (40, 60);
        let (q_start, q_end, _, t_start_out, t_end_out) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();

        // Query should be 40-70 (10bp for first part + 10bp insertion + 10bp for second part = 30bp)
        assert_eq!(t_start_out, 40);
        assert_eq!(t_end_out, 60);

        // Inverse: query range → target
        let inverse = project_query_range_through_alignment((q_start, q_end), record, &cigar_ops);
        assert!(inverse.is_some());
        let (t_start_back, t_end_back, _, _, _) = inverse.unwrap();

        // Should get back to original target range
        assert_eq!((t_start_back, t_end_back), (t_start_out, t_end_out));
    }

    #[test]
    fn test_round_trip_with_deletion() {
        // CIGAR: 50=10D50= (target is 10bp longer due to deletion in query)
        // Target: 0-110, Query: 0-100
        let cigar_ops = vec![
            CigarOp::new(50, '='),
            CigarOp::new(10, 'D'), // 10bp deletion in query (exists in target only)
            CigarOp::new(50, '='),
        ];
        let record = (0, 110, 0, 100, Strand::Forward);

        // Forward: target 40-70 (spans the deletion at 50-60)
        let target_range = (40, 70);
        let (q_start, q_end, _, t_start_out, t_end_out) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();

        // Query should be 40-60 (10bp match + deletion skipped + 10bp match = 20bp)
        assert_eq!(t_start_out, 40);
        assert_eq!(t_end_out, 70);

        // Inverse: query range → target
        let inverse = project_query_range_through_alignment((q_start, q_end), record, &cigar_ops);
        assert!(inverse.is_some());
        let (t_start_back, t_end_back, _, _, _) = inverse.unwrap();

        // Should get back to original target range (including the deletion)
        assert_eq!((t_start_back, t_end_back), (t_start_out, t_end_out));
    }

    #[test]
    fn test_round_trip_complex_cigar() {
        // Complex CIGAR: 10=5I5D50=
        let cigar_ops = vec![
            CigarOp::new(10, '='),
            CigarOp::new(5, 'I'),  // Query advances 5, target stays
            CigarOp::new(5, 'D'),  // Target advances 5, query stays
            CigarOp::new(50, '='),
        ];
        let record = (0, 65, 0, 65, Strand::Forward);

        // Test multiple ranges
        for (t_start, t_end) in [(5, 20), (0, 65), (15, 40), (50, 60)] {
            let target_range = (t_start, t_end);
            let forward = project_target_range_through_alignment(target_range, record, &cigar_ops);

            if let Some((q_start, q_end, _, t_start_out, t_end_out)) = forward {
                let inverse = project_query_range_through_alignment((q_start, q_end), record, &cigar_ops);

                if let Some((t_start_back, t_end_back, _, _, _)) = inverse {
                    // The round-trip should preserve the target range (at least the overlapping portion)
                    assert_eq!(
                        (t_start_back, t_end_back),
                        (t_start_out, t_end_out),
                        "Round-trip failed for target range {:?}",
                        target_range
                    );
                }
            }
        }
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
                target_id,
                target_start: 30,
                target_end: 40,
                strand_and_data_offset: 45, // Forward strand
                data_bytes: 3,
            },
            // Add more test records as needed
        ];
        let records = parse_paf(reader, &mut seq_index).unwrap();
        assert_eq!(records, expected_records);
    }

    // Tests for QueryMetadata bit operations (bidirectional indexing support)
    #[test]
    fn test_query_metadata_reversed_bit_default_false() {
        let metadata = QueryMetadata {
            query_id: 1,
            target_start: 0,
            target_end: 100,
            query_start: 0,
            query_end: 100,
            alignment_file_index: 0,
            strand_and_data_offset: 42, // No special bits set
            data_bytes: 10,
        };
        assert!(!metadata.is_reversed());
        assert_eq!(metadata.strand(), Strand::Forward);
        assert_eq!(metadata.data_offset(), 42);
    }

    #[test]
    fn test_query_metadata_reversed_bit_set() {
        let metadata = QueryMetadata {
            query_id: 1,
            target_start: 0,
            target_end: 100,
            query_start: 0,
            query_end: 100,
            alignment_file_index: 0,
            strand_and_data_offset: 42 | QueryMetadata::REVERSED_BIT,
            data_bytes: 10,
        };
        assert!(metadata.is_reversed());
        assert_eq!(metadata.strand(), Strand::Forward);
        assert_eq!(metadata.data_offset(), 42);
    }

    #[test]
    fn test_query_metadata_strand_and_reversed_independent() {
        // Both strand and reversed bits set
        let metadata = QueryMetadata {
            query_id: 1,
            target_start: 0,
            target_end: 100,
            query_start: 0,
            query_end: 100,
            alignment_file_index: 0,
            strand_and_data_offset: 42 | QueryMetadata::STRAND_BIT | QueryMetadata::REVERSED_BIT,
            data_bytes: 10,
        };
        assert!(metadata.is_reversed());
        assert_eq!(metadata.strand(), Strand::Reverse);
        assert_eq!(metadata.data_offset(), 42);
    }

    #[test]
    fn test_query_metadata_data_offset_masks_both_bits() {
        // Test that data_offset correctly masks both high bits
        let base_offset: u64 = 0x1234567890ABCDEF & !(QueryMetadata::STRAND_BIT | QueryMetadata::REVERSED_BIT);

        // Only strand bit
        let metadata1 = QueryMetadata {
            query_id: 1,
            target_start: 0,
            target_end: 100,
            query_start: 0,
            query_end: 100,
            alignment_file_index: 0,
            strand_and_data_offset: base_offset | QueryMetadata::STRAND_BIT,
            data_bytes: 10,
        };
        assert_eq!(metadata1.data_offset(), base_offset);

        // Only reversed bit
        let metadata2 = QueryMetadata {
            query_id: 1,
            target_start: 0,
            target_end: 100,
            query_start: 0,
            query_end: 100,
            alignment_file_index: 0,
            strand_and_data_offset: base_offset | QueryMetadata::REVERSED_BIT,
            data_bytes: 10,
        };
        assert_eq!(metadata2.data_offset(), base_offset);

        // Both bits
        let metadata3 = QueryMetadata {
            query_id: 1,
            target_start: 0,
            target_end: 100,
            query_start: 0,
            query_end: 100,
            alignment_file_index: 0,
            strand_and_data_offset: base_offset | QueryMetadata::STRAND_BIT | QueryMetadata::REVERSED_BIT,
            data_bytes: 10,
        };
        assert_eq!(metadata3.data_offset(), base_offset);
    }
}
