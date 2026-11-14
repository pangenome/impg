use crate::alignment_record::{AlignmentRecord, Strand};
use crate::forest_map::ForestMap;
use crate::graph::reverse_complement;
use crate::impg_index::ImpgIndex;
use crate::onealn::{
    OneAlnAlignment, OneAlnParser, ParseErr as OneAlnParseErr, TracepointModeData,
};
use crate::paf::read_cigar_data;
use crate::seqidx::SequenceIndex;
use crate::sequence_index::SequenceIndex as _; // The as _ syntax imports the trait so its methods are available, but doesn't bring the name into scope (avoiding the naming conflict)
use crate::sequence_index::UnifiedSequenceIndex;
use crate::subset_filter::SubsetFilter;
use crate::tpa_parser::{ParseErr as TpaParseErr, TpaParser};
use coitrees::{BasicCOITree, GenericInterval, Interval, IntervalTree};
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
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, RwLock};
use tracepoints::{tracepoints_to_cigar_fastga_with_aligner, tracepoints_to_cigar_with_aligner};

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
    static GAP_AFFINE_ALIGNER: RefCell<Option<AffineWavefronts>> = const { RefCell::new(None) };
    static GAP_AFFINE2P_ALIGNER: RefCell<Option<AffineWavefronts>> = const { RefCell::new(None) };
    static ONEALN_HANDLE: RefCell<Option<(String, OneFile)>> = const { RefCell::new(None) };
    static PAF_HANDLE: RefCell<Option<(String, std::fs::File)>> = const { RefCell::new(None) };
    static PAF_CIGAR_BUF: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
    static TPA_HANDLE: RefCell<Option<(String, tpa::TpaReader)>> = const { RefCell::new(None) };
    /// Last target sequence fetched on this thread.
    ///
    /// The key must include the sequence name, not only its numeric id:
    /// `MultiImpg` sub-indices use dense per-file ids, so `(id, start, end,
    /// reverse)` alone collides across genomes and would serve another
    /// sample's DNA to CIGAR reconstruction.
    static TARGET_SEQ_CACHE: RefCell<Option<TargetSeqCacheEntry>> = const { RefCell::new(None) };
}

struct TargetSeqCacheEntry {
    name: String,
    target_id: u32,
    start: i32,
    end: i32,
    is_reverse: bool,
    seq: Vec<u8>,
}

impl TargetSeqCacheEntry {
    fn matches(&self, name: &str, target_id: u32, start: i32, end: i32, is_reverse: bool) -> bool {
        self.target_id == target_id
            && self.start == start
            && self.end == end
            && self.is_reverse == is_reverse
            && self.name == name
    }
}

/// Execute a closure with a thread-local aligner matched to the given distance metric
fn with_aligner<F, R>(distance: &Distance, f: F) -> R
where
    F: FnOnce(&mut AffineWavefronts) -> R,
{
    let (cell, dist) = match distance {
        d @ Distance::Edit => (&EDIT_ALIGNER, d),
        d @ Distance::GapAffine { .. } => (&GAP_AFFINE_ALIGNER, d),
        d @ Distance::GapAffine2p { .. } => (&GAP_AFFINE2P_ALIGNER, d),
    };
    cell.with(|aligner_cell| {
        let mut opt = aligner_cell.borrow_mut();
        if opt.is_none() {
            *opt = Some(dist.create_aligner(None, None));
        }
        f(opt.as_mut().unwrap())
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

/// Invert CIGAR ops for bidirectional interpretation: swap I<->D, keep =/X/M,
/// and reverse the array on the reverse strand.
fn invert_cigar_ops_in_place(ops: &mut [CigarOp], strand: Strand) {
    for op in ops.iter_mut() {
        let new_op = match op.op() {
            'I' => 'D',
            'D' => 'I',
            other => other,
        };
        *op = CigarOp::new(op.len(), new_op);
    }
    if strand == Strand::Reverse {
        ops.reverse();
    }
}

fn invert_cigar_ops(ops: &[CigarOp], strand: Strand) -> Vec<CigarOp> {
    let mut inverted = ops.to_vec();
    invert_cigar_ops_in_place(&mut inverted, strand);
    inverted
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
    const STRAND_BIT: u64 = 0x8000000000000000; // Most significant bit for u64 (bit 63)
    const REVERSED_BIT: u64 = 0x4000000000000000; // Second-most significant bit for reversed flag (bit 62)

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

    fn data_offset(&self) -> u64 {
        self.strand_and_data_offset & !(Self::STRAND_BIT | Self::REVERSED_BIT)
    }

    /// Key for `cigar_cache`: preserves `REVERSED_BIT` alongside `(alignment_file_index, data_offset)`
    /// because cached CIGARs are stored in the perspective (forward vs reversed) of the entry.
    fn cigar_cache_key(&self) -> (u32, u64) {
        (
            self.alignment_file_index,
            self.strand_and_data_offset & !Self::STRAND_BIT,
        )
    }

    fn is_reversed(&self) -> bool {
        (self.strand_and_data_offset & Self::REVERSED_BIT) != 0
    }

    fn set_reversed(&mut self) {
        self.strand_and_data_offset |= Self::REVERSED_BIT;
    }

    /// Determine file type based on extension
    fn get_file_type(alignment_file: &str) -> FileType {
        if alignment_file.ends_with(".1aln") {
            FileType::OneAln
        } else if alignment_file.ends_with(".tpa") {
            FileType::Tpa
        } else {
            FileType::Paf
        }
    }
}

pub type AdjustedInterval = (Interval<u32>, Vec<CigarOp>, Interval<u32>);

/// Reject inputs that would only fail later, deep inside a query or a worker.
///
/// Coordinates are stored as `i32`, so a sequence or a coordinate above
/// `i32::MAX` would wrap to a negative value without any message. A PAF record
/// without a `cg:Z` tag cannot be projected at all: the projection reads the
/// CIGAR to place the query interval inside each segment, and the old code
/// panicked inside a rayon worker when it reached such a record.
fn validate_index_inputs(
    records_by_file: &[(Vec<AlignmentRecord>, String)],
    seq_index: &SequenceIndex,
) -> io::Result<()> {
    for id in 0..seq_index.len() as u32 {
        if let Some(len) = seq_index.get_len_from_id(id) {
            if len > i32::MAX as usize {
                let name = seq_index.get_name(id).unwrap_or("<unknown>");
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "sequence '{name}' is {len} bp, above the maximum supported length of {} bp",
                        i32::MAX
                    ),
                ));
            }
        }
    }

    for (records, path) in records_by_file {
        for record in records {
            let max_coord = record
                .query_start
                .max(record.query_end)
                .max(record.target_start)
                .max(record.target_end);
            if max_coord > i32::MAX as usize {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "an alignment record in '{path}' has the coordinate {max_coord}, above the maximum supported value {} bp",
                        i32::MAX
                    ),
                ));
            }
        }

        if !matches!(QueryMetadata::get_file_type(path), FileType::Paf) {
            continue;
        }
        let missing = records.iter().filter(|r| r.data_bytes == 0).count();
        if missing > 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "the alignment file '{path}' has {missing} PAF record(s) out of {} without a CIGAR ('cg:Z' tag). impg projects intervals through the CIGAR, so these records cannot be indexed. Re-run the alignment step so that it writes CIGAR strings.",
                    records.len()
                ),
            ));
        }
    }

    Ok(())
}
type TreeMap = FxHashMap<u32, Arc<BasicCOITree<QueryMetadata, u32>>>;

struct CachedTree {
    tree: Arc<BasicCOITree<QueryMetadata, u32>>,
    bytes: usize,
    /// Logical timestamp of the last access, for LRU eviction. Atomic so a lookup
    /// only needs a read lock.
    last_used: AtomicU64,
}

/// Per-target interval trees held in memory, with optional LRU eviction.
///
/// A query only needs the trees of the targets it touches, but a multi-target run
/// used to keep every tree it ever loaded. Trees read from an index file can be
/// dropped and reloaded on demand, so the cache is bounded. Indexes built in memory
/// have no file to reload from (their forest-map offsets are placeholders), so they
/// are never evicted.
pub struct TreeCache {
    /// Key: cache namespace (one per sub-index) in the high 32 bits, target id below
    map: FxHashMap<u64, CachedTree>,
    capacity_bytes: Option<usize>,
    bytes: usize,
    clock: AtomicU64,
}

impl TreeCache {
    /// Cache for an in-memory index: holds everything, never evicts.
    fn unbounded(trees: TreeMap) -> Self {
        let mut cache = Self {
            map: FxHashMap::default(),
            capacity_bytes: None,
            bytes: 0,
            clock: AtomicU64::new(0),
        };
        for (target_id, tree) in trees {
            cache.insert(target_id as u64, tree);
        }
        cache
    }

    /// Cache for a file-backed index: evicts least-recently-used trees past `capacity_bytes`.
    pub(crate) fn bounded(capacity_bytes: usize) -> Self {
        Self {
            map: FxHashMap::default(),
            capacity_bytes: Some(capacity_bytes),
            bytes: 0,
            clock: AtomicU64::new(0),
        }
    }

    /// Approximate resident cost of a tree. A coitrees node holds the metadata plus
    /// subtree_last, first, last and two child indices, so about 24 bytes of overhead
    /// per record. This tracks the steady-state size of cached trees; the transient
    /// peak while a tree is being decoded and built is several times larger.
    fn estimate_bytes(tree: &BasicCOITree<QueryMetadata, u32>) -> usize {
        tree.len() * (std::mem::size_of::<QueryMetadata>() + 24)
    }

    fn get(&self, key: u64) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        self.map.get(&key).map(|entry| {
            entry
                .last_used
                .store(self.clock.fetch_add(1, Ordering::Relaxed), Ordering::Relaxed);
            Arc::clone(&entry.tree)
        })
    }

    fn insert(&mut self, key: u64, tree: Arc<BasicCOITree<QueryMetadata, u32>>) {
        let bytes = Self::estimate_bytes(&tree);
        let stamp = self.clock.fetch_add(1, Ordering::Relaxed);
        if let Some(previous) = self.map.insert(
            key,
            CachedTree {
                tree,
                bytes,
                last_used: AtomicU64::new(stamp),
            },
        ) {
            self.bytes -= previous.bytes;
        }
        self.bytes += bytes;
        self.evict_to_capacity(key);
    }

    /// Drop least-recently-used entries until the cache fits. Never evicts
    /// `keep`, the entry just inserted, so a single oversized tree still works.
    fn evict_to_capacity(&mut self, keep: u64) {
        let Some(capacity) = self.capacity_bytes else {
            return;
        };
        if self.bytes <= capacity {
            return;
        }
        // One scan and one sort per overflow, not one scan per victim
        let mut by_age: Vec<(u64, u64)> = self
            .map
            .iter()
            .filter(|(&key, _)| key != keep)
            .map(|(&key, entry)| (entry.last_used.load(Ordering::Relaxed), key))
            .collect();
        by_age.sort_unstable();
        for (_, key) in by_age {
            if self.bytes <= capacity {
                break;
            }
            if let Some(entry) = self.map.remove(&key) {
                self.bytes -= entry.bytes;
                debug!("tree cache: evicted {key:#x} ({} B), {} B cached", entry.bytes, self.bytes);
            }
        }
    }

    fn remove(&mut self, key: u64) {
        if let Some(entry) = self.map.remove(&key) {
            self.bytes -= entry.bytes;
        }
    }

    fn iter(&self) -> impl Iterator<Item = (&u64, &Arc<BasicCOITree<QueryMetadata, u32>>)> {
        self.map.iter().map(|(id, entry)| (id, &entry.tree))
    }
}

/// Capacity for the file-backed tree cache: `IMPG_TREE_CACHE_MB` if set, else a
/// quarter of the memory this process can use (physical RAM, or the cgroup limit
/// when it is lower), with a 512 MB floor. Computed once per process. In per-file
/// mode all sub-indices share one cache with this budget (`MultiImpg`).
pub(crate) fn default_tree_cache_bytes() -> usize {
    static BYTES: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *BYTES.get_or_init(|| {
        const FLOOR: usize = 512 * 1024 * 1024;
        if let Ok(value) = std::env::var("IMPG_TREE_CACHE_MB") {
            match value.trim().parse::<usize>() {
                Ok(mb) => return mb.saturating_mul(1024 * 1024),
                Err(e) => warn!("IMPG_TREE_CACHE_MB={value:?} is not a whole number of MB ({e}); using the default"),
            }
        }
        // SAFETY: sysconf only reads system configuration
        let (pages, page_size) =
            unsafe { (libc::sysconf(libc::_SC_PHYS_PAGES), libc::sysconf(libc::_SC_PAGESIZE)) };
        let physical = if pages > 0 && page_size > 0 { pages as usize * page_size as usize } else { 0 };
        // cgroup v2, then v1; "max" or a huge v1 value means no limit
        let cgroup = ["/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"]
            .iter()
            .find_map(|p| std::fs::read_to_string(p).ok()?.trim().parse::<usize>().ok());
        let usable = match cgroup {
            Some(limit) if physical == 0 || limit < physical => limit,
            _ => physical,
        };
        let bytes = (usable / 4).max(FLOOR);
        debug!("tree cache budget {} MB (usable memory {} MB)", bytes >> 20, usable >> 20);
        bytes
    })
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum FileType {
    Paf,
    OneAln,
    Tpa,
}

#[derive(Serialize, Deserialize)]
struct SerializableInterval {
    first: i32,
    last: i32,
    metadata: QueryMetadata,
}

/// Lets a decoded interval be fed straight to `BasicCOITree::new`, without first
/// rebuilding the whole slice as `Vec<Interval<QueryMetadata>>`.
impl GenericInterval<QueryMetadata> for SerializableInterval {
    fn first(&self) -> i32 {
        self.first
    }
    fn last(&self) -> i32 {
        self.last
    }
    fn metadata(&self) -> &QueryMetadata {
        &self.metadata
    }
}

#[derive(Default, Clone)]
pub struct SortedRanges {
    pub ranges: Vec<(i32, i32)>,
    sequence_length: i32,
    min_distance: i32,
}

impl SortedRanges {
    pub fn with_distance(sequence_length: i32, min_distance: i32) -> Self {
        Self {
            ranges: Vec::new(),
            sequence_length,
            min_distance,
        }
    }

    pub fn new() -> Self {
        Self {
            ranges: Vec::new(),
            sequence_length: 0,
            min_distance: 0,
        }
    }

    /// Insert a range, snapping its boundaries onto adjacent existing ranges
    /// when they lie within `min_distance` (used for expansion ranges, to absorb
    /// alignment boundary jitter while keeping novel territory).
    ///
    /// User-requested ranges must go through `insert`, which never snaps.
    pub fn insert_with_boundary_snap(
        &mut self,
        new_range: (i32, i32),
        min_distance: i32,
    ) -> Vec<(i32, i32)> {
        if min_distance <= 0 {
            return self.insert(new_range);
        }
        self.insert_inner(new_range, min_distance)
    }

    fn insert_with_bounds(
        &mut self,
        new_range: (i32, i32),
        min_distance: i32,
        clamp_seq_ends: bool,
    ) -> Vec<(i32, i32)> {
        let (mut start, mut end) = if new_range.0 <= new_range.1 {
            (new_range.0, new_range.1)
        } else {
            (new_range.1, new_range.0)
        };

        let i = match self.ranges.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(pos) => pos,
            Err(pos) => pos,
        };

        // Snap `start` back onto the previous range's end when the gap is small.
        // Only ever extend: never pull `start` forward over novel territory.
        if i > 0 {
            let prev_end = self.ranges[i - 1].1;
            if start > prev_end && start - prev_end < min_distance {
                start = prev_end;
            }
        } else if clamp_seq_ends && self.sequence_length > 0 && start < min_distance {
            start = 0;
        }

        // Snap `end` forward onto the next range's start when the gap is small.
        // The successor must be searched from `end`: `ranges[i]` can start inside
        // the range being inserted, and `end` may only ever grow. Using
        // `(ranges[i].0 - end).abs()` truncated `end` and silently dropped the
        // novel tail of the inserted range.
        if let Some(&(next_start, _)) = self.ranges[i..].iter().find(|(s, _)| *s >= end) {
            if next_start > end && next_start - end < min_distance {
                end = next_start;
            }
        } else if clamp_seq_ends
            && self.sequence_length > 0
            && end > (self.sequence_length - min_distance)
        {
            end = self.sequence_length;
        }

        // Collect regions that don't overlap with existing ranges
        let mut non_overlapping = Vec::new();
        let mut current = start;
        let mut idx = match self.ranges.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(pos) => pos,
            Err(pos) => pos,
        };
        if idx > 0 && self.ranges[idx - 1].1 > start {
            idx -= 1;
        }

        while idx < self.ranges.len() && current < end {
            let (range_start, range_end) = self.ranges[idx];
            if range_start > end {
                break;
            }
            if current < range_start {
                non_overlapping.push((current, range_start));
            }
            current = max(current, range_end);
            idx += 1;
        }
        if current < end {
            non_overlapping.push((current, end));
        }

        match self.ranges.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(pos) | Err(pos) => {
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

        non_overlapping
    }

    fn insert_inner(&mut self, new_range: (i32, i32), min_distance: i32) -> Vec<(i32, i32)> {
        self.insert_with_bounds(new_range, min_distance, false)
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
        self.insert_with_bounds(new_range, self.min_distance, true)
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
    /// Scan direction along metadata target axis (`+1` or `-1`)
    scan_dir: i32,
    /// Project direction along metadata query axis (`+1` or `-1`)
    project_dir: i32,
}

impl SubsettingResult {
    /// Compute normalized `(scan_start, scan_end, proj_lo, proj_hi)` bounds for the subset.
    /// When the subset reaches the last tracepoint (`tracepoints_len - 1`), clamp the trailing
    /// endpoint along `scan_dir` / `project_dir` to the exact metadata boundary.
    fn bounds(&self, tracepoints_len: usize, metadata: &QueryMetadata) -> (i32, i32, i32, i32) {
        let (_, _, first_seg_start, first_seg_end, _, _) = self.first_segment_info;
        let (_, _, last_seg_start, last_seg_end, _, _) = self.last_segment_info;
        let mut scan_start = first_seg_start
            .min(first_seg_end)
            .min(last_seg_start.min(last_seg_end));
        let mut scan_end = first_seg_start
            .max(first_seg_end)
            .max(last_seg_start.max(last_seg_end));

        let mut proj_end = self.last_query_pos;
        if self.last_idx == tracepoints_len - 1 {
            if self.scan_dir == 1 {
                scan_end = metadata.target_end;
            } else {
                scan_start = metadata.target_start;
            }
            proj_end = if self.project_dir == 1 {
                metadata.query_end
            } else {
                metadata.query_start
            };
        }

        let proj_lo = self.first_query_pos.min(proj_end);
        let proj_hi = self.first_query_pos.max(proj_end);
        (scan_start, scan_end, proj_lo, proj_hi)
    }
}

/// One alignment file as recorded in the index header.
///
/// The index stores a file index and a byte offset for every CIGAR, so the
/// order of the file list is part of the index. Without this record, passing
/// the same files in a different order reads CIGARs from the wrong file
/// without any error. `size` is `None` when the size cannot be read, for
/// example for a process-substitution path; such an entry cannot be compared.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct AlignmentFileEntry {
    pub path: String,
    pub size: Option<u64>,
}

impl AlignmentFileEntry {
    fn from_path(path: &str) -> Self {
        let size = std::fs::metadata(path).ok().map(|m| m.len());
        // Canonical, so the same file compares equal however it is spelled
        let path = std::fs::canonicalize(path)
            .map(|p| p.to_string_lossy().into_owned())
            .unwrap_or_else(|_| path.to_string());
        AlignmentFileEntry { path, size }
    }
}

/// Header block written by the current index format.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct IndexHeaderInfo {
    /// Alignment files in the order the index was built from.
    pub alignment_files: Vec<AlignmentFileEntry>,
    /// True when the index holds both directions of every record.
    pub bidirectional: bool,
}

/// Describe a whole recorded list for an error message.
fn describe_files_from(files: &[AlignmentFileEntry]) -> String {
    files
        .iter()
        .map(|f| match f.size {
            Some(size) => format!("{} ({} bytes)", f.path, size),
            None => f.path.clone(),
        })
        .collect::<Vec<_>>()
        .join(", ")
}

/// Refuse an index built from a different alignment file list.
///
/// A reordered list is the dangerous case: every file index still resolves, so
/// the query reads CIGARs from the wrong file and returns wrong intervals with
/// no message.
pub(crate) fn validate_alignment_file_list(
    recorded: &[AlignmentFileEntry],
    given: &[String],
    index_file_path: &str,
) -> io::Result<()> {
    if recorded.len() != given.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "The index '{index_file_path}' was built from {} alignment file(s) but {} were given. \
                 The index stores CIGAR offsets into those files, so the list and its order must match. \
                 Rebuild the index with `impg index` on this list, or pass the original one.",
                recorded.len(),
                given.len()
            ),
        ));
    }

    let mismatch = recorded
        .iter()
        .zip(given.iter())
        .enumerate()
        .find(|(_, (recorded_entry, path))| {
            let given_entry = AlignmentFileEntry::from_path(path);
            // Both paths are canonical, so `data/x.paf`, `./x.paf` from inside
            // `data/`, an absolute path and a symlink compare equal. When the
            // recorded path is gone (the files were moved), fall back to the
            // file name.
            let name = |p: &str| std::path::Path::new(p).file_name().map(|n| n.to_owned());
            let same_path = recorded_entry.path == given_entry.path
                || (!std::path::Path::new(&recorded_entry.path).exists()
                    && name(&recorded_entry.path) == name(&given_entry.path));
            match (recorded_entry.size, given_entry.size) {
                (Some(a), Some(b)) => a != b || !same_path,
                // A size that could not be read (a pipe, for example) cannot be compared
                _ => false,
            }
        });

    if let Some((position, (recorded_entry, path))) = mismatch {
        let given_entry = AlignmentFileEntry::from_path(path);
        let given_entries: Vec<AlignmentFileEntry> = given
            .iter()
            .map(|p| AlignmentFileEntry::from_path(p))
            .collect();
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "The index '{index_file_path}' was built from a different alignment file list. \
                 Entry {position} is '{}' in the index but '{}' was given. \
                 The index stores CIGAR offsets into those files, so a different list or a different \
                 order reads the wrong CIGAR. \
                 Recorded: [{}]. Given now: [{}]. \
                 Rebuild the index with `impg index` on this list, or pass the original one.",
                describe_files_from(std::slice::from_ref(recorded_entry)),
                describe_files_from(std::slice::from_ref(&given_entry)),
                describe_files_from(recorded),
                describe_files_from(&given_entries)
            ),
        ));
    }

    Ok(())
}

/// Identify one edge of the alignment graph, for deduplication.
///
/// Two entries with the same key project to the same interval pair.
fn edge_key(interval: &Interval<QueryMetadata>) -> (u32, i32, i32, i32, i32, bool) {
    let metadata = &interval.metadata;
    (
        metadata.query_id(),
        metadata.target_start,
        metadata.target_end,
        metadata.query_start,
        metadata.query_end,
        metadata.strand() == Strand::Reverse,
    )
}

fn with_cached_handle<H: 'static, OpenFn, UseFn, R>(
    cell: &'static std::thread::LocalKey<RefCell<Option<(String, H)>>>,
    file_path: &str,
    open_fn: OpenFn,
    use_fn: UseFn,
) -> Result<R, String>
where
    OpenFn: FnOnce(&str) -> Result<H, String>,
    UseFn: FnOnce(&mut H) -> Result<R, String>,
{
    cell.with(|handle_cell| {
        let mut slot = handle_cell.borrow_mut();
        let needs_open = match slot.as_ref() {
            Some((cached_path, _)) => cached_path != file_path,
            None => true,
        };

        if needs_open {
            *slot = None;
            let handle = open_fn(file_path)?;
            *slot = Some((file_path.to_string(), handle));
        }

        let (_, handle) = slot
            .as_mut()
            .expect("cached handle slot must contain a handle after opening");
        use_fn(handle)
    })
}

pub(crate) fn make_self_interval(
    target_id: u32,
    range_start: i32,
    range_end: i32,
    store_cigar: bool,
) -> AdjustedInterval {
    (
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
    )
}

pub(crate) trait QueueRangeItem {
    fn key_start(&self) -> (u32, i32);
    fn range_end(&self) -> i32;
    fn extend_end(&mut self, new_end: i32);
}

impl QueueRangeItem for (u32, i32, i32) {
    fn key_start(&self) -> (u32, i32) {
        (self.0, self.1)
    }
    fn range_end(&self) -> i32 {
        self.2
    }
    fn extend_end(&mut self, new_end: i32) {
        self.2 = new_end;
    }
}

impl QueueRangeItem for (u32, i32, i32, u16) {
    fn key_start(&self) -> (u32, i32) {
        (self.0, self.1)
    }
    fn range_end(&self) -> i32 {
        self.2
    }
    fn extend_end(&mut self, new_end: i32) {
        self.2 = new_end;
    }
}

/// Merge contiguous or overlapping ranges in a slice already sorted by `(id, start)`.
/// Returns the compacted prefix length to pass to `truncate`.
pub(crate) fn compact_sorted_ranges<T: QueueRangeItem>(items: &mut [T]) -> usize {
    if items.is_empty() {
        return 0;
    }
    let mut write = 0;
    for read in 1..items.len() {
        if items[write].key_start().0 == items[read].key_start().0
            && items[write].range_end() >= items[read].key_start().1
        {
            let merged_end = items[write].range_end().max(items[read].range_end());
            items[write].extend_end(merged_end);
        } else {
            write += 1;
            items.swap(write, read);
        }
    }
    write + 1
}

pub(crate) struct LoadedIndexHeader {
    pub(crate) seq_index: SequenceIndex,
    pub(crate) forest_map: ForestMap,
    pub(crate) built_bidirectional: bool,
}

/// Read magic bytes, optional V4 header (validating the alignment file list),
/// sequence index, and forest map from an `.impg` index stream.
pub(crate) fn read_index_header<R: std::io::Read + std::io::Seek>(
    reader: &mut R,
    alignment_files: &[String],
    index_file_path: &str,
    warn_on_legacy: bool,
) -> std::io::Result<LoadedIndexHeader> {
    const MAGIC_V1: &[u8] = b"IMPGIDX1";
    const MAGIC_V2: &[u8] = b"IMPGIDX2";
    const MAGIC_V4: &[u8] = b"IMPGIDX4";

    let mut magic_buf = [0u8; 8];
    reader.read_exact(&mut magic_buf)?;
    if magic_buf != MAGIC_V1 && magic_buf != MAGIC_V2 && magic_buf != MAGIC_V4 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Invalid magic bytes - not a valid IMPG index file '{}' (got {:?})",
                index_file_path,
                String::from_utf8_lossy(&magic_buf)
            ),
        ));
    }

    if warn_on_legacy && magic_buf == MAGIC_V1 {
        warn!("Loading IMPGIDX1 format (unidirectional). Rebuild without --unidirectional for better transitive query reach (bidirectional mode is default).");
    }

    let mut offset_buf = [0u8; 8];
    reader.read_exact(&mut offset_buf)?;
    let forest_map_offset = u64::from_le_bytes(offset_buf);

    let header_info: Option<IndexHeaderInfo> = if magic_buf == MAGIC_V4 {
        Some(
            bincode::serde::decode_from_std_read(&mut *reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load index header from '{index_file_path}': {e}"),
                    )
                })?,
        )
    } else {
        if warn_on_legacy {
            warn!(
                "The index '{}' does not record its alignment file list, so impg cannot check that the files given now are the ones it was built from. Rebuild the index to enable the check.",
                index_file_path
            );
        }
        None
    };

    if let Some(header) = &header_info {
        validate_alignment_file_list(&header.alignment_files, alignment_files, index_file_path)?;
        if warn_on_legacy {
            debug!(
                "Index '{}' was built {} from {} alignment file(s)",
                index_file_path,
                if header.bidirectional {
                    "bidirectionally"
                } else {
                    "unidirectionally"
                },
                header.alignment_files.len()
            );
        }
    }

    let seq_index: SequenceIndex =
        bincode::serde::decode_from_std_read(&mut *reader, bincode::config::standard()).map_err(
            |e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to load sequence index: {e}"),
                )
            },
        )?;

    reader.seek(SeekFrom::Start(forest_map_offset))?;
    let forest_map: ForestMap =
        bincode::serde::decode_from_std_read(&mut *reader, bincode::config::standard()).map_err(
            |e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to load forest map: {e}"),
                )
            },
        )?;

    Ok(LoadedIndexHeader {
        seq_index,
        forest_map,
        built_bidirectional: header_info
            .map(|h| h.bidirectional)
            .unwrap_or(magic_buf != MAGIC_V1),
    })
}

pub struct Impg {
    pub trees: Arc<RwLock<TreeCache>>,
    /// Namespace of this index inside `trees`; 0 unless the cache is shared
    cache_ns: u32,
    pub seq_index: SequenceIndex,
    alignment_files: Vec<String>,
    /// True when this index holds both directions of every record.
    built_bidirectional: bool,
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
        let spacing = if file_path.ends_with(".tpa") {
            TpaParser::read_trace_spacing(file_path)
                .map_err(|e| format!("Failed to read trace_spacing from '{}': {}", file_path, e))?
        } else {
            OneAlnParser::read_trace_spacing(file_path).map_err(|e| {
                format!("Failed to read trace_spacing from '{}': {:?}", file_path, e)
            })?
        };

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
        with_cached_handle(
            &ONEALN_HANDLE,
            file_path,
            |path| {
                OneFile::open_read(path, None, None, 1)
                    .map_err(|e| format!("Failed to open 1aln file '{}': {}", path, e))
            },
            |handle| f(handle).map_err(|e| e.to_string()),
        )
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
        TARGET_SEQ_CACHE.with(|cache_cell| {
            let mut slot = cache_cell.borrow_mut();
            if let Some(entry) = slot.as_ref() {
                if entry.matches(target_name, target_id, start, end, is_reverse) {
                    return entry.seq.clone();
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

            *slot = Some(TargetSeqCacheEntry {
                name: target_name.to_string(),
                target_id,
                start,
                end,
                is_reverse,
                seq: seq.clone(),
            });
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

        let mut ops = match QueryMetadata::get_file_type(alignment_file) {
            FileType::Paf => {
                // For PAF files, read CIGAR directly from file
                if metadata.data_bytes == 0 {
                    panic!(
                        "The alignment file '{}' does not contain CIGAR strings ('cg:Z' tag).",
                        alignment_file
                    );
                }

                // Read the CIGAR bytes into a per-thread scratch buffer reused across
                // projections, avoiding a malloc/free of up to tens of KB per call.
                PAF_CIGAR_BUF.with(|buf_cell| {
                    let mut data_buffer = buf_cell.borrow_mut();
                    data_buffer.resize(metadata.data_bytes, 0);
                    read_paf_cigar_data(alignment_file, metadata.data_offset(), &mut data_buffer)
                        .unwrap_or_else(|e| panic!("{}", e));

                    // get_cigar_ops_from_bytes
                    let cigar_str = std::str::from_utf8(&data_buffer)
                        .expect("Failed to parse CIGAR data as UTF-8");
                    parse_cigar_to_delta(cigar_str).unwrap_or_else(|e| {
                        panic!(
                            "Failed to parse CIGAR string '{}' in QueryMetadata: {:?}",
                            cigar_str, e
                        )
                    })
                })
            }
            FileType::OneAln | FileType::Tpa => {
                // For 1aln/TPA files, convert tracepoints to CIGAR
                let alignment = self
                    .get_tracepoint_alignment(metadata)
                    .unwrap_or_else(|e| panic!("{}", e));

                self.process_tracepoints_data(
                    &alignment,
                    metadata,
                    target_id,
                    sequence_index.expect("Sequence index required for tracepoint files"),
                )
            }
        };

        // Apply CIGAR inversion if this is a reversed entry (in place, no extra alloc)
        if metadata.is_reversed() {
            invert_cigar_ops_in_place(&mut ops, metadata.strand());
        }
        ops
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

    fn with_tpa_reader<F, R>(&self, file_index: usize, f: F) -> Result<R, String>
    where
        F: FnOnce(&mut tpa::TpaReader) -> Result<R, TpaParseErr>,
    {
        let file_path = &self.alignment_files[file_index];
        with_cached_handle(
            &TPA_HANDLE,
            file_path,
            |path| {
                let mut reader = tpa::TpaReader::new(path)
                    .map_err(|e| format!("Failed to open TPA file '{}': {}", path, e))?;
                reader
                    .load_string_table()
                    .map_err(|e| format!("Failed to load TPA string table '{}': {}", path, e))?;
                Ok(reader)
            },
            |handle| f(handle).map_err(|e| e.to_string()),
        )
    }

    fn get_tpa_alignment(&self, metadata: &QueryMetadata) -> Result<OneAlnAlignment, String> {
        let file_index = metadata.alignment_file_index as usize;
        let alignment_file = &self.alignment_files[file_index];
        let record_id = metadata.data_offset();

        self.with_tpa_reader(file_index, |reader| {
            TpaParser::fetch_alignment(reader, record_id)
        })
        .map_err(|e| {
            format!(
                "Failed to fetch TPA alignment {} from '{}': {}",
                record_id, alignment_file, e
            )
        })
    }

    /// Get tracepoint alignment from either .1aln or .tpa file
    fn get_tracepoint_alignment(
        &self,
        metadata: &QueryMetadata,
    ) -> Result<OneAlnAlignment, String> {
        let file = &self.alignment_files[metadata.alignment_file_index as usize];
        if file.ends_with(".tpa") {
            self.get_tpa_alignment(metadata)
        } else if file.ends_with(".1aln") {
            self.get_onealn_alignment(metadata)
        } else {
            Err(format!("No tracepoint support for: {}", file))
        }
    }

    /// Compute `(scan_dir, project_dir)` along metadata's `(target, query)` axes
    /// as tracepoint segments are traversed from index `0` to `len - 1`.
    ///
    /// On forward strand (`+`), both formats traverse original `(query, target)`
    /// ascending (`+1, +1`). On reverse strand (`-`), FastGA traverses original
    /// query ascending (`+1`) vs original target descending (`-1`), whereas
    /// Standard/Variable/Mixed TPA stores PAF CIGAR order: original query
    /// descending (`-1`) vs original target ascending (`+1`). Reversed entries
    /// swap which axis is scanned vs projected.
    fn tracepoint_directions(alignment: &OneAlnAlignment, metadata: &QueryMetadata) -> (i32, i32) {
        let is_reverse = metadata.strand() == Strand::Reverse;
        let is_standard = matches!(alignment.mode_data, TracepointModeData::Standard { .. });
        let orig_q_dir = if is_reverse && is_standard { -1 } else { 1 };
        let orig_t_dir = if is_reverse && !is_standard { -1 } else { 1 };
        if metadata.is_reversed() {
            (orig_q_dir, orig_t_dir)
        } else {
            (orig_t_dir, orig_q_dir)
        }
    }

    /// Scan tracepoints to find segments overlapping the requested target range.
    ///
    /// For forward entries: scans along the target axis, projects onto the query axis.
    /// For reversed entries: scans along the query axis (stored as target in metadata),
    /// projects onto the target axis (stored as query in metadata).
    /// Same tracepoint segments, different scan dimension
    fn scan_overlapping_tracepoints(
        &self,
        alignment: &OneAlnAlignment,
        metadata: &QueryMetadata,
        range_start: i32,
        range_end: i32,
        is_reversed_entry: bool,
    ) -> Option<SubsettingResult> {
        let (scan_dir, project_dir) = Self::tracepoint_directions(alignment, metadata);

        // Initialize scan position (along range axis) and project position (along output axis)
        let mut scan_pos = if scan_dir == 1 {
            metadata.target_start
        } else {
            metadata.target_end
        };
        let mut project_pos = if project_dir == 1 {
            metadata.query_start
        } else {
            metadata.query_end
        };

        // Track overlapping segments
        let mut first_idx: Option<usize> = None;
        let mut last_idx: usize = 0;
        let mut first_project_pos: Option<i32> = None;
        let mut last_project_pos: Option<i32> = None;
        let mut first_segment_info: Option<(i32, i32, i32, i32, i32, i32)> = None;
        let mut last_segment_info: Option<(i32, i32, i32, i32, i32, i32)> = None;
        let mut num_overlapping_segments = 0;

        // Statistics for identity calculation
        let mut total_matches = 0.0;
        let mut total_mismatches = 0.0;

        // Pre-compute FASTGA-specific values and remaining original query span for clamping
        let (trace_spacing, first_boundary) = match &alignment.mode_data {
            TracepointModeData::Fastga { trace_spacing, .. } => {
                let ts = *trace_spacing as i32;
                let query_start_in_contig = alignment.query_contig_start as i32;
                let fb = ((query_start_in_contig / ts) + 1) * ts - query_start_in_contig;
                (ts, fb)
            }
            TracepointModeData::Standard { .. } => (0, 0), // not used
        };
        let orig_q_total = if is_reversed_entry {
            metadata.target_end - metadata.target_start
        } else {
            metadata.query_end - metadata.query_start
        };
        let mut consumed_orig_q = 0;

        // Scan tracepoints and collect those that overlap the requested range
        for (idx, &tracepoint) in alignment.tracepoints.iter().enumerate() {
            // Query delta (unsigned magnitude along original query axis)
            let query_delta = match &alignment.mode_data {
                TracepointModeData::Standard { query_deltas, .. } => query_deltas[idx] as i32,
                TracepointModeData::Fastga { .. } => {
                    let raw = if idx == 0 {
                        first_boundary
                    } else {
                        trace_spacing
                    };
                    let clamped = raw.min(orig_q_total - consumed_orig_q);
                    consumed_orig_q += clamped;
                    clamped
                }
            };
            let abs_tracepoint = tracepoint.abs() as i32;

            // Compute signed scan and project deltas based on entry type and directions
            let (scan_delta, project_delta, abs_scan_delta) = if is_reversed_entry {
                (query_delta * scan_dir, abs_tracepoint * project_dir, query_delta)
            } else {
                (abs_tracepoint * scan_dir, query_delta * project_dir, abs_tracepoint)
            };

            // Segment boundaries (forward coordinate space along scan axis)
            let seg_start = scan_pos.min(scan_pos + scan_delta);
            let seg_end = scan_pos.max(scan_pos + scan_delta);

            // Check overlap with requested range
            if seg_start < range_end && seg_end > range_start {
                num_overlapping_segments += 1;
                let num_diffs = match &alignment.mode_data {
                    TracepointModeData::Fastga { diffs, .. } => {
                        if idx < diffs.len() {
                            diffs[idx] as i32
                        } else {
                            0
                        }
                    }
                    TracepointModeData::Standard { max_complexity, .. } => {
                        // Estimate: use max_complexity unless one delta is 0
                        if query_delta == 0 || abs_tracepoint == 0 {
                            query_delta.max(abs_tracepoint)
                        } else {
                            *max_complexity as i32
                        }
                    }
                };
                let seg_info = (
                    project_pos,
                    project_delta,
                    seg_start,
                    seg_end,
                    abs_scan_delta,
                    num_diffs,
                );

                // Track first and last overlapping segments
                if first_project_pos.is_none() {
                    first_idx = Some(idx);
                    first_project_pos = Some(project_pos);
                    first_segment_info = Some(seg_info);
                }
                last_idx = idx;
                last_project_pos = Some(project_pos + project_delta);
                last_segment_info = Some(seg_info);

                // Accumulate identity statistics
                let aligned_len = (project_delta.abs().min(abs_scan_delta) as f64).max(0.0);
                total_matches += (aligned_len - num_diffs as f64).max(0.0);
                total_mismatches += num_diffs as f64;
            }

            // Advance to next segment
            scan_pos += scan_delta;
            project_pos += project_delta;

            // Early exit when past requested range
            if (scan_dir == -1 && scan_pos <= range_start)
                || (scan_dir == 1 && scan_pos >= range_end)
            {
                break;
            }
        }

        // Return None if no overlapping segments found
        first_idx?;

        Some(SubsettingResult {
            first_idx: first_idx.unwrap(),
            last_idx,
            first_query_pos: first_project_pos.unwrap(),
            last_query_pos: last_project_pos.unwrap(),
            num_overlapping_segments,
            total_matches,
            total_mismatches,
            first_segment_info: first_segment_info.unwrap(),
            last_segment_info: last_segment_info.unwrap(),
            scan_dir,
            project_dir,
        })
    }

    /// Fetch `(a_seq, b_seq)` in the orientation expected by tracepoint WFA reconstruction.
    ///
    /// Always fetches original query (`a_seq`) and original target (`b_seq`) even when
    /// `metadata.is_reversed()` has swapped `query_id`/`target_id` in metadata.
    /// On reverse strand, FastGA expects `(fwd(orig_q), RC(orig_t))`, while Standard/
    /// Variable/Mixed TPA expects `(RC(orig_q), fwd(orig_t))` (PAF CIGAR order).
    #[allow(clippy::too_many_arguments)]
    fn fetch_tracepoint_sequences(
        &self,
        alignment: &OneAlnAlignment,
        metadata: &QueryMetadata,
        target_id: u32,
        orig_q_lo: i32,
        orig_q_hi: i32,
        orig_t_lo: i32,
        orig_t_hi: i32,
        sequence_index: &UnifiedSequenceIndex,
    ) -> (Vec<u8>, Vec<u8>) {
        let is_reverse = metadata.strand() == Strand::Reverse;
        let is_standard = matches!(alignment.mode_data, TracepointModeData::Standard { .. });
        let rc_orig_q = is_reverse && is_standard;
        let rc_orig_t = is_reverse && !is_standard;

        let (orig_q_id, orig_t_id) = if metadata.is_reversed() {
            (target_id, metadata.query_id)
        } else {
            (metadata.query_id, target_id)
        };

        let a_name = self.seq_index.get_name(orig_q_id).unwrap();
        let a_seq_raw = sequence_index
            .fetch_sequence(a_name, orig_q_lo, orig_q_hi)
            .unwrap_or_else(|e| panic!("Failed to fetch query sequence: {e}"));
        let a_seq = if rc_orig_q {
            reverse_complement(&a_seq_raw)
        } else {
            a_seq_raw
        };

        let b_name = self.seq_index.get_name(orig_t_id).unwrap();
        let b_seq = self.get_target_sequence_cached(
            sequence_index,
            b_name,
            orig_t_id,
            orig_t_lo,
            orig_t_hi,
            rc_orig_t,
        );

        (a_seq, b_seq)
    }

    /// Reconstruct a CIGAR vector in the original alignment perspective from a
    /// contiguous tracepoint slice `[start_idx..end_idx]` and aligned sequence slices.
    #[allow(clippy::too_many_arguments)]
    fn reconstruct_cigar_from_tracepoints(
        alignment: &OneAlnAlignment,
        start_idx: usize,
        end_idx: usize,
        a_seq: &[u8],
        b_seq: &[u8],
        a_contig_offset: usize,
        b_contig_offset: usize,
        is_reverse: bool,
    ) -> Vec<CigarOp> {
        let cigar_str = match &alignment.mode_data {
            TracepointModeData::Standard {
                query_deltas,
                complexity_metric,
                max_complexity,
                distance,
            } => {
                let tracepoints: Vec<(usize, usize)> = query_deltas[start_idx..end_idx]
                    .iter()
                    .zip(alignment.tracepoints[start_idx..end_idx].iter())
                    .map(|(&qd, &td)| (qd as usize, td as usize))
                    .collect();
                let metric = *complexity_metric;
                let max_val = Some(*max_complexity as u32);
                with_aligner(distance, |aligner| {
                    tracepoints_to_cigar_with_aligner(
                        &tracepoints,
                        a_seq,
                        b_seq,
                        0,
                        0,
                        metric,
                        aligner,
                        max_val,
                    )
                })
            }
            TracepointModeData::Fastga {
                diffs,
                trace_spacing,
            } => {
                let tracepoints: Vec<(usize, usize)> = diffs[start_idx..end_idx]
                    .iter()
                    .zip(alignment.tracepoints[start_idx..end_idx].iter())
                    .map(|(&diff, &tp)| (diff as usize, tp as usize))
                    .collect();
                let ts = *trace_spacing as usize;
                with_aligner(&Distance::Edit, |aligner| {
                    tracepoints_to_cigar_fastga_with_aligner(
                        &tracepoints,
                        ts.try_into().unwrap(),
                        a_seq,
                        b_seq,
                        a_contig_offset,
                        b_contig_offset,
                        is_reverse,
                        aligner,
                        true,
                    )
                })
            }
        };

        parse_cigar_to_delta(&cigar_str).unwrap_or_else(|e| {
            panic!("Failed to parse CIGAR string '{cigar_str}' during tracepoint processing: {e:?}")
        })
    }

    fn process_tracepoints_data(
        &self,
        alignment: &OneAlnAlignment,
        metadata: &QueryMetadata,
        target_id: u32,
        sequence_index: &UnifiedSequenceIndex,
    ) -> Vec<CigarOp> {
        let (orig_q_lo, orig_q_hi, orig_t_lo, orig_t_hi) = if metadata.is_reversed() {
            (
                metadata.target_start,
                metadata.target_end,
                metadata.query_start,
                metadata.query_end,
            )
        } else {
            (
                metadata.query_start,
                metadata.query_end,
                metadata.target_start,
                metadata.target_end,
            )
        };

        let (a_seq, b_seq) = self.fetch_tracepoint_sequences(
            alignment,
            metadata,
            target_id,
            orig_q_lo,
            orig_q_hi,
            orig_t_lo,
            orig_t_hi,
            sequence_index,
        );

        // Shortcut to a pure-match CIGAR only when both diff metadata and the
        // fetched sequences confirm an exact match of equal length (guards against
        // stale TPA metadata or boundary normalization offsets).
        if alignment.differences == 0 && a_seq.len() == b_seq.len() && a_seq.eq_ignore_ascii_case(&b_seq) {
            return vec![CigarOp::new(a_seq.len() as i32, '=')];
        }

        Self::reconstruct_cigar_from_tracepoints(
            alignment,
            0,
            alignment.tracepoints.len(),
            &a_seq,
            &b_seq,
            alignment.query_contig_start.try_into().unwrap(),
            alignment.target_contig_start.try_into().unwrap(),
            metadata.strand() == Strand::Reverse,
        )
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
        let (subset_scan_start, subset_scan_end, subset_proj_lo, subset_proj_hi) =
            subset.bounds(alignment.tracepoints.len(), metadata);

        let (orig_q_lo, orig_q_hi, orig_t_lo, orig_t_hi, orig_q_start, orig_t_start) =
            if metadata.is_reversed() {
                (
                    subset_scan_start,
                    subset_scan_end,
                    subset_proj_lo,
                    subset_proj_hi,
                    metadata.target_start,
                    metadata.query_start,
                )
            } else {
                (
                    subset_proj_lo,
                    subset_proj_hi,
                    subset_scan_start,
                    subset_scan_end,
                    metadata.query_start,
                    metadata.target_start,
                )
            };

        let (a_seq, b_seq) = self.fetch_tracepoint_sequences(
            alignment,
            metadata,
            target_id,
            orig_q_lo,
            orig_q_hi,
            orig_t_lo,
            orig_t_hi,
            sequence_index,
        );

        let qoff = alignment.query_contig_start as usize + (orig_q_lo - orig_q_start) as usize;
        let toff = alignment.target_contig_start as usize + (orig_t_lo - orig_t_start) as usize;

        Self::reconstruct_cigar_from_tracepoints(
            alignment,
            subset.first_idx,
            subset.last_idx + 1,
            &a_seq,
            &b_seq,
            qoff,
            toff,
            is_reverse,
        )
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
        // Check if this is a tracepoint file that supports subsetting
        let file_index = metadata.alignment_file_index as usize;
        let alignment_file = &self.alignment_files[file_index];
        let is_tracepoint_file = matches!(
            QueryMetadata::get_file_type(alignment_file),
            FileType::OneAln | FileType::Tpa
        );

        // Try subsetting approach for tracepoint files
        if let (true, Some(sequence_index)) = (is_tracepoint_file, sequence_index) {
            // Pre-filter: coitrees uses closed-interval overlap, but coordinates are half-open.
            // Skip alignments that merely touch at an endpoint (no actual overlap).
            if metadata.target_start >= range_end || metadata.target_end <= range_start {
                return None;
            }

            // Fetch alignment with tracepoints
            let alignment = match self.get_tracepoint_alignment(metadata) {
                Ok(a) => a,
                Err(e) => {
                    debug!("Cannot fetch tracepoint alignment for target_id={}, file_index={}: {}, skipping",
                        target_id, metadata.alignment_file_index, e);
                    return None;
                }
            };

            // Scan tracepoints to find overlapping segments
            let subset = match self.scan_overlapping_tracepoints(
                &alignment,
                metadata,
                range_start,
                range_end,
                metadata.is_reversed(),
            ) {
                Some(s) => s,
                None => {
                    debug!(
                        "No overlapping tracepoint segments for range {}-{} (target_id={}, file_index={}), skipping",
                        range_start, range_end, target_id, metadata.alignment_file_index
                    );
                    return None;
                }
            };

            // Reconstruct CIGAR from subset tracepoints
            // Note: When last_idx is the actual last tracepoint, process_subset_tracepoints
            // uses metadata boundaries to handle the potentially shorter last segment correctly
            let cigar_ops = self.process_subset_tracepoints(
                &alignment,
                metadata,
                target_id,
                &subset,
                sequence_index,
            );

            // For reversed entries, the CIGAR was reconstructed in the original alignment
            // perspective. Invert it to the metadata's (reversed) perspective so that
            // project_target_range_through_alignment works correctly.
            let cigar_ops = if metadata.is_reversed() {
                invert_cigar_ops(&cigar_ops, metadata.strand())
            } else {
                cigar_ops
            };

            // Check identity threshold if specified
            if let Some(threshold) = min_gap_compressed_identity {
                if calculate_gap_compressed_identity(&cigar_ops) < threshold {
                    return None;
                }
            }

            let (
                alignment_target_start,
                alignment_target_end,
                alignment_query_start,
                alignment_query_end,
            ) = subset.bounds(alignment.tracepoints.len(), metadata);

            // Project the requested range through the CIGAR
            let projection = match project_target_range_through_alignment(
                (range_start, range_end),
                (
                    alignment_target_start,
                    alignment_target_end,
                    alignment_query_start,
                    alignment_query_end,
                    metadata.strand(),
                ),
                &cigar_ops,
            ) {
                Some(p) => p,
                None => {
                    debug!(
                        "Subsetting projection failed for tracepoint alignment (range={}-{}, query={}-{}, target={}-{})",
                        range_start, range_end, alignment_query_start, alignment_query_end, alignment_target_start, alignment_target_end
                    );
                    return None;
                }
            };

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

        let projection = project_target_range_through_alignment(
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
        // Pre-filter: coitrees uses closed-interval overlap, but coordinates are half-open.
        // Skip alignments that merely touch at an endpoint (no actual overlap).
        if metadata.target_start >= range_end || metadata.target_end <= range_start {
            return None;
        }

        // Fetch tracepoints without sequence I/O
        let alignment = match self.get_tracepoint_alignment(metadata) {
            Ok(a) => a,
            Err(e) => {
                debug!("Cannot fetch tracepoint alignment for target_id={}, file_index={}: {}, skipping",
                    target_id, metadata.alignment_file_index, e);
                return None;
            }
        };

        // Scan tracepoints to find overlapping segments
        let subset = match self.scan_overlapping_tracepoints(
            &alignment,
            metadata,
            range_start,
            range_end,
            metadata.is_reversed(),
        ) {
            Some(s) => s,
            None => {
                debug!(
                    "No overlapping tracepoint segments for range {}-{} (target_id={}, file_index={}), skipping",
                    range_start, range_end, target_id, metadata.alignment_file_index
                );
                return None;
            }
        };

        let is_reverse = metadata.strand() == Strand::Reverse;
        let working_query_start = metadata.query_start;
        let working_query_end = metadata.query_end;

        // Store pre-refinement values for debugging
        let pre_refinement_query = (subset.first_query_pos, subset.last_query_pos);
        let pre_refinement_target = (range_start, range_end);

        // Refine query coordinates using indel-aware heuristic.
        // When scanning in -1 direction, query_pos corresponds to segment_target_end
        // (and first_segment hits range_end while last_segment hits range_start).
        let refine_boundary = |query_pos: i32,
                               query_delta: i32,
                               segment_target_start: i32,
                               segment_target_end: i32,
                               overlap_pos: i32,
                               abs_target_delta: i32,
                               trace_diffs: i32,
                               boundary_name: &str|
         -> i32 {
            let aligned_len = query_delta.abs().min(abs_target_delta);
            let segment_identity = if aligned_len > 0 {
                ((aligned_len - trace_diffs) as f64 / aligned_len as f64).max(0.0)
            } else {
                0.0
            };

            if abs_target_delta == 0 {
                // Pure insertion in query: zero target bases consumed.
                let refined_pos = if boundary_name == "First" {
                    query_pos
                } else {
                    query_pos + query_delta
                };
                debug!(
                    "{} segment refinement: pure query insertion, abs_target_delta=0, query_delta={}, trace_diffs={}",
                    boundary_name, query_delta, trace_diffs
                );
                let lo = working_query_start.min(working_query_end);
                let hi = working_query_start.max(working_query_end);
                return refined_pos.max(lo).min(hi);
            }

            let target_offset = if subset.scan_dir == 1 {
                overlap_pos - segment_target_start
            } else {
                segment_target_end - overlap_pos
            };
            let target_fraction = target_offset as f64 / abs_target_delta as f64;
            let indel_ratio = query_delta as f64 / abs_target_delta as f64;
            let query_advance = target_fraction * abs_target_delta as f64 * indel_ratio;
            let refined_pos = query_pos + query_advance.round() as i32;

            debug!(
                "{} segment refinement: identity={:.3}, indel_ratio={:.3}, trace_diffs={}, target_frac={:.3}",
                boundary_name, segment_identity, indel_ratio, trace_diffs, target_fraction
            );

            let lo = working_query_start.min(working_query_end);
            let hi = working_query_start.max(working_query_end);
            refined_pos.max(lo).min(hi)
        };

        // Refine first and last query positions
        let (
            query_pos,
            query_delta,
            segment_target_start,
            segment_target_end,
            abs_target_delta,
            trace_diffs,
        ) = subset.first_segment_info;
        let first_overlap = if subset.scan_dir == 1 {
            segment_target_start.max(range_start)
        } else {
            segment_target_end.min(range_end)
        };
        let refined_first = refine_boundary(
            query_pos,
            query_delta,
            segment_target_start,
            segment_target_end,
            first_overlap,
            abs_target_delta,
            trace_diffs,
            "First",
        );

        let (
            query_pos,
            query_delta,
            segment_target_start,
            segment_target_end,
            abs_target_delta,
            trace_diffs,
        ) = subset.last_segment_info;
        let last_overlap = if subset.scan_dir == 1 {
            segment_target_end.min(range_end)
        } else {
            segment_target_start.max(range_start)
        };
        let refined_last = refine_boundary(
            query_pos,
            query_delta,
            segment_target_start,
            segment_target_end,
            last_overlap,
            abs_target_delta,
            trace_diffs,
            "Last",
        );

        // Debug: show before/after refinement
        let post_refinement_query = (refined_first, refined_last);
        let post_refinement_target = (range_start, range_end);
        let overall_identity = if subset.total_matches + subset.total_mismatches > 0.0 {
            subset.total_matches / (subset.total_matches + subset.total_mismatches)
        } else {
            0.0
        };
        let requested_range_len = range_end - range_start;
        let strand_char = if is_reverse { '-' } else { '+' };

        debug!(
            "Refinement [{}]: req_len={}, query {}-{} (Δ{}) -> {}-{} (Δ{}), target {}-{} (Δ{}) -> {}-{} (Δ{}), segs={}, id={:.3}, diffs={:.0}",
            strand_char,
            requested_range_len,
            pre_refinement_query.0, pre_refinement_query.1,
            (pre_refinement_query.1 - pre_refinement_query.0).abs(),
            post_refinement_query.0, post_refinement_query.1,
            (post_refinement_query.1 - post_refinement_query.0).abs(),
            pre_refinement_target.0, pre_refinement_target.1,
            pre_refinement_target.1 - pre_refinement_target.0,
            post_refinement_target.0, post_refinement_target.1,
            post_refinement_target.1 - post_refinement_target.0,
            subset.num_overlapping_segments,
            overall_identity,
            subset.total_mismatches
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

        // Ensure reverse-strand query intervals are returned as (high, low).
        // When project_dir == 1, refined_first < refined_last so swap them;
        // when project_dir == -1, refined_first is already the higher coordinate.
        let (query_start, query_end) = if is_reverse && subset.project_dir == 1 {
            (working_query_end, working_query_start)
        } else {
            (working_query_start, working_query_end)
        };

        // The target interval is the part of the requested range that this
        // alignment actually covers. Copying the requested range here claimed
        // homology outside the alignment: with a reverse record over target
        // 100-900, querying 0-1000 returned target 0-1000 while the CIGAR path
        // returned 100-900.
        let target_start = range_start.max(metadata.target_start);
        let target_end = range_end.min(metadata.target_end);

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

    pub fn from_multi_alignment_records(
        records_by_file: &[(Vec<AlignmentRecord>, String)],
        seq_index: SequenceIndex,
        sequence_files: Option<&[String]>,
        bidirectional: bool,
    ) -> io::Result<Self> {
        validate_index_inputs(records_by_file, &seq_index)?;

        // Extract just the alignment file paths
        let alignment_files: Vec<String> = records_by_file
            .par_iter()
            .map(|(_, alignment_file)| alignment_file.clone())
            .collect();

        if bidirectional {
            debug!("Creating bidirectional alignment entries");
        } else {
            debug!("Creating unidirectional alignment entries");
        }

        const RECORD_CHUNK_SIZE: usize = 4096;
        let record_chunks: Vec<(u32, &[AlignmentRecord])> = records_by_file
            .iter()
            .enumerate()
            .flat_map(|(file_index, (records, _))| {
                let file_idx = file_index as u32;
                records
                    .chunks(RECORD_CHUNK_SIZE)
                    .map(move |chunk| (file_idx, chunk))
            })
            .collect();

        let intervals: FxHashMap<u32, Vec<Interval<QueryMetadata>>> = record_chunks
            .into_par_iter()
            .fold(
                FxHashMap::default,
                |mut acc: FxHashMap<u32, Vec<Interval<QueryMetadata>>>, (file_index, chunk)| {
                    for record in chunk {
                        // ENTRY 1: Original direction (query → target)
                        let forward_metadata = QueryMetadata {
                            query_id: record.query_id,
                            target_start: record.target_start as i32,
                            target_end: record.target_end as i32,
                            query_start: record.query_start as i32,
                            query_end: record.query_end as i32,
                            alignment_file_index: file_index,
                            strand_and_data_offset: record.strand_and_data_offset,
                            data_bytes: record.data_bytes,
                        };
                        acc.entry(record.target_id).or_default().push(Interval {
                            first: record.target_start as i32,
                            last: record.target_end as i32,
                            metadata: forward_metadata,
                        });

                        // ENTRY 2: Reversed direction (target → query)
                        // Skip self-alignments to avoid duplicates
                        if bidirectional && record.query_id != record.target_id {
                            let mut reversed_metadata = QueryMetadata {
                                query_id: record.target_id,              // SWAPPED
                                target_start: record.query_start as i32, // SWAPPED
                                target_end: record.query_end as i32,     // SWAPPED
                                query_start: record.target_start as i32, // SWAPPED
                                query_end: record.target_end as i32,     // SWAPPED
                                alignment_file_index: file_index,
                                strand_and_data_offset: record.strand_and_data_offset,
                                data_bytes: record.data_bytes,
                            };
                            reversed_metadata.set_reversed(); // Mark with REVERSED_BIT

                            acc.entry(record.query_id).or_default().push(Interval {
                                first: record.query_start as i32,
                                last: record.query_end as i32,
                                metadata: reversed_metadata,
                            });
                        }
                    }
                    acc
                },
            )
            .reduce(FxHashMap::default, |mut acc, part| {
                for (key, mut value) in part {
                    acc.entry(key)
                        .and_modify(|existing| existing.extend(value.drain(..)))
                        .or_insert(value);
                }
                acc
            });

        // A reciprocal pair (`A -> B` and `B -> A` for the same alignment)
        // reaches the same edge twice: once as the reversed entry of the first
        // record and once as the forward entry of the second. Both describe the
        // same projection, and both reached the output, so bedpe and paf
        // printed the row twice. The released set has 1,483 such pairs.
        //
        // The key ignores the file index, the byte offset and the reversed
        // flag: those say which record a projection reads its CIGAR from, not
        // which edge this is.
        let trees: TreeMap = intervals
            .into_par_iter()
            .map(|(target_id, mut interval_nodes)| {
                interval_nodes.sort_unstable_by_key(edge_key);
                interval_nodes.dedup_by(|a, b| edge_key(a) == edge_key(b));
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
            trees: Arc::new(RwLock::new(TreeCache::unbounded(trees))),
            cache_ns: 0,
            seq_index,
            alignment_files,
            built_bidirectional: bidirectional,
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
        const MAGIC: &[u8] = b"IMPGIDX4";
        let mut forest_map = ForestMap::new();

        // Write magic bytes
        writer.write_all(MAGIC)?;

        // Write placeholder for forest map offset (will be updated later)
        writer.write_all(&[0u8; 8])?;

        let mut current_offset = 16u64; // Magic plus forest map offset

        // Record which alignment files this index was built from, in order: the
        // byte offsets in the trees are only valid for that list.
        let header_info = IndexHeaderInfo {
            alignment_files: self
                .alignment_files
                .iter()
                .map(|path| AlignmentFileEntry::from_path(path))
                .collect(),
            bidirectional: self.built_bidirectional,
        };
        let header_data =
            bincode::serde::encode_to_vec(&header_info, bincode::config::standard()).map_err(
                |e| std::io::Error::other(format!("Failed to encode index header: {e}")),
            )?;
        writer.write_all(&header_data)?;
        current_offset += header_data.len() as u64;

        // Serialize sequence index
        let seq_index_data =
            bincode::serde::encode_to_vec(&self.seq_index, bincode::config::standard()).map_err(
                |e| std::io::Error::other(format!("Failed to encode sequence index: {e}")),
            )?;

        writer.write_all(&seq_index_data)?;
        current_offset += seq_index_data.len() as u64;

        // Serialize each tree and track its offset
        let trees = self.trees.read().unwrap();
        for (&key, tree) in trees.iter() {
            if (key >> 32) as u32 != self.cache_ns {
                continue;
            }
            let target_id = key as u32;
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
        // An in-memory index has no file to reload from. Callers can evict a
        // tree and ask for it again (`impg stats` does), so return nothing
        // instead of opening an empty path.
        if self.index_file_path.is_empty() {
            return None;
        }
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

            // Reconstruct the tree directly from the decoded records
            let tree = BasicCOITree::new(&intervals);

            let arc_tree = Arc::new(tree);

            // Cache the tree for future use
            self.trees
                .write()
                .unwrap()
                .insert(self.cache_key(target_id), Arc::clone(&arc_tree));

            Some(arc_tree)
        } else {
            None // Tree not found in forest map
        }
    }

    /// Get a tree from memory or load it from disk if necessary
    fn cache_key(&self, target_id: u32) -> u64 {
        ((self.cache_ns as u64) << 32) | target_id as u64
    }

    /// Use `cache` instead of this index's own tree cache, under namespace `ns`.
    /// `MultiImpg` gives every sub-index the same cache, so one budget bounds
    /// them all and the least recently used tree of any sub-index goes first.
    pub(crate) fn share_tree_cache(&mut self, cache: Arc<RwLock<TreeCache>>, ns: u32) {
        self.trees = cache;
        self.cache_ns = ns;
    }

    pub fn get_or_load_tree(
        &self,
        target_id: u32,
    ) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        // First check if the tree is already in memory.
        // The cache hands back a clone of the Arc (cheap refcount bump) so the RwLock
        // is not held while the tree is used, and so eviction can never invalidate a
        // tree a caller is still holding.
        if let Some(tree) = self.trees.read().unwrap().get(self.cache_key(target_id)) {
            return Some(tree);
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
        let header = read_index_header(&mut reader, alignment_files, &index_file_path, true)?;
        let num_files = alignment_files.len();
        Ok(Self {
            trees: Arc::new(RwLock::new(TreeCache::bounded(default_tree_cache_bytes()))),
            cache_ns: 0,
            seq_index: header.seq_index,
            alignment_files: alignment_files.to_vec(),
            built_bidirectional: header.built_bidirectional,
            forest_map: header.forest_map,
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
        let mut results = vec![make_self_interval(
            target_id,
            range_start,
            range_end,
            store_cigar,
        )];

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
                let cache_key = metadata.cigar_cache_key();
                cache
                    .entry(cache_key)
                    .or_insert_with(|| self.get_cigar_ops(metadata, target_id, sequence_index));
            });
        }
    }

    /// Query a range, reusing CIGARs from `cigar_cache` instead of reading them
    /// again.
    ///
    /// With `approximate_mode`, the projection interpolates the query interval
    /// from the segment boundaries and never reads a CIGAR or a sequence, so
    /// `cigar_cache` and `sequence_index` are not used.
    #[allow(clippy::too_many_arguments)]
    pub fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
        approximate_mode: bool,
    ) -> Vec<AdjustedInterval> {
        let mut results = vec![make_self_interval(
            target_id,
            range_start,
            range_end,
            store_cigar,
        )];

        if let Some(tree) = self.get_or_load_tree(target_id) {
            tree.query(range_start, range_end, |interval| {
                let metadata = &interval.metadata;

                if approximate_mode {
                    // Fast path: no CIGAR and no sequence I/O.
                    if let Some((query_interval, cigar_ops, target_interval)) = self
                        .project_overlapping_interval_fast(
                            metadata,
                            target_id,
                            range_start,
                            range_end,
                            min_gap_compressed_identity,
                        )
                    {
                        let cigar_vec = if store_cigar { cigar_ops } else { Vec::new() };
                        results.push((query_interval, cigar_vec, target_interval));
                    }
                    return;
                }

                let cache_key = metadata.cigar_cache_key();
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

    /// Get or lazily create a sequence's visited-range set during a transitive walk.
    /// A new entry gets the sequence's real length when unmasked (matching the old
    /// eager pre-population) and length 0 otherwise (matching the old `or_default()`).
    #[inline]
    fn visited_entry<'a>(
        &self,
        visited_ranges: &'a mut FxHashMap<u32, SortedRanges>,
        id: u32,
        masked_none: bool,
    ) -> &'a mut SortedRanges {
        visited_ranges.entry(id).or_insert_with(|| {
            let len = if masked_none {
                self.seq_index.get_len_from_id(id).unwrap() as i32
            } else {
                0
            };
            SortedRanges::with_distance(len, 0)
        })
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
        // Seed from masked regions if given; other sequences are created lazily on
        // first touch (visited_entry) instead of pre-populating one SortedRanges per
        // sequence in the index on every query.
        let masked_none = masked_regions.is_none();
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter().map(|(&k, v)| (k, (*v).clone())).collect()
        } else {
            FxHashMap::default()
        };

        // Filter input range
        let filtered_input_range = self
            .visited_entry(&mut visited_ranges, target_id, masked_none)
            .insert((range_start, range_end));

        let mut results = Vec::new();
        let mut stack = Vec::new();

        for (filtered_start, filtered_end) in filtered_input_range {
            // Add the filtered input range(s) to the results
            results.push(make_self_interval(
                target_id,
                filtered_start,
                filtered_end,
                store_cigar,
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
                        projection_result.and_then(
                            |(query_interval, cigar_ops, target_interval)| {
                                let query_id = query_interval.metadata;

                                // Apply subset filter: keep if it's the target or matches the filter
                                let should_keep = if let Some(filter) = subset_filter {
                                    query_id == target_id
                                        || self
                                            .seq_index
                                            .get_name(query_id)
                                            .is_some_and(|name| filter.matches(name))
                                } else {
                                    true
                                };

                                if should_keep {
                                    let adjusted_query_start = query_interval.first;
                                    let adjusted_query_end = query_interval.last;
                                    let cigar_vec =
                                        if store_cigar { cigar_ops } else { Vec::new() };

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
                            },
                        )
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
                        results.push((query_interval, cigar, target_interval));
                    }

                    // Only add non-overlapping portions to the stack for further exploration
                    if query_id != current_target_id {
                        let ranges = self.visited_entry(&mut visited_ranges, query_id, masked_none);

                        let new_ranges = ranges.insert_with_boundary_snap(
                            (adjusted_query_start, adjusted_query_end),
                            min_distance_between_ranges,
                        );

                        // Add non-overlapping portions to stack
                        for (new_start, new_end) in new_ranges {
                            if (new_end - new_start).abs() >= min_transitive_len {
                                stack.push((query_id, new_start, new_end, current_depth + 1));
                            }
                        }
                    }
                }
            }

            //debug!("Collected {} results", results.len() - prec_num_results);

            // Merge contiguous/overlapping ranges with same sequence_id
            stack.par_sort_by_key(QueueRangeItem::key_start);
            let new_len = compact_sorted_ranges(&mut stack);
            stack.truncate(new_len);
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
        // Seed from masked regions if given; other sequences are created lazily on
        // first touch (visited_entry) instead of pre-populating one SortedRanges per
        // sequence in the index on every query.
        let masked_none = masked_regions.is_none();
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter().map(|(&k, v)| (k, (*v).clone())).collect()
        } else {
            FxHashMap::default()
        };

        // Filter input range
        let filtered_input_range = self
            .visited_entry(&mut visited_ranges, target_id, masked_none)
            .insert((range_start, range_end));

        let mut results = Vec::new();

        // Add the filtered input range(s) to the results
        for &(filtered_start, filtered_end) in &filtered_input_range {
            results.push(make_self_interval(
                target_id,
                filtered_start,
                filtered_end,
                store_cigar,
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
                                                query_id == target_id
                                                    || self
                                                        .seq_index
                                                        .get_name(query_id)
                                                        .is_some_and(|name| filter.matches(name))
                                            } else {
                                                true
                                            };

                                            if should_keep {
                                                let cigar_vec = if store_cigar {
                                                    cigar_ops
                                                } else {
                                                    Vec::new()
                                                };

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
                            adjusted_cigar,
                            Interval {
                                first: adjusted_target_start,
                                last: adjusted_target_end,
                                metadata: current_target_id,
                            },
                        ));
                    }

                    // Only consider for next depth if it's a different sequence
                    if query_id != current_target_id {
                        let ranges = self.visited_entry(&mut visited_ranges, query_id, masked_none);

                        let new_ranges = ranges.insert_with_boundary_snap(
                            (adjusted_query_start, adjusted_query_end),
                            min_distance_between_ranges,
                        );

                        // Add non-overlapping portions to next depth
                        for (new_start, new_end) in new_ranges {
                            if (new_end - new_start).abs() >= min_transitive_len {
                                next_depth_ranges.push((query_id, new_start, new_end));
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
                next_depth_ranges.par_sort_by_key(QueueRangeItem::key_start);
                let new_len = compact_sorted_ranges(&mut next_depth_ranges);
                next_depth_ranges.truncate(new_len);
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
    ) -> io::Result<Vec<AdjustedInterval>> {
        // Delegate to the existing method
        Ok(Impg::query(
            self,
            target_id,
            range_start,
            range_end,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            approximate_mode,
        ))
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
        approximate_mode: bool,
    ) -> io::Result<Vec<AdjustedInterval>> {
        Ok(Impg::query_with_cache(
            self,
            target_id,
            range_start,
            range_end,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            cigar_cache,
            approximate_mode,
        ))
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
    ) -> io::Result<Vec<AdjustedInterval>> {
        Ok(Impg::query_transitive_dfs(
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
        ))
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
    ) -> io::Result<Vec<AdjustedInterval>> {
        Ok(Impg::query_transitive_bfs(
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
        ))
    }

    fn get_or_load_tree(&self, target_id: u32) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        Impg::get_or_load_tree(self, target_id)
    }

    fn target_ids(&self) -> Vec<u32> {
        self.forest_map.entries.keys().copied().collect()
    }

    fn remove_cached_tree(&self, target_id: u32) {
        self.trees.write().unwrap().remove(self.cache_key(target_id));
    }

    fn num_targets(&self) -> usize {
        self.forest_map.entries.len()
    }

    fn sequence_files(&self) -> &[String] {
        &self.sequence_files
    }

    fn can_reload_tree(&self, target_id: u32) -> bool {
        // Offset 0 marks a tree that only exists in memory (see
        // `from_multi_alignment_records`), so there is nothing to read back.
        !self.index_file_path.is_empty()
            && self
                .forest_map
                .get_tree_offset(target_id)
                .is_some_and(|offset| offset > 0)
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

/// Read a PAF CIGAR span into `buffer`. Uncompressed PAF uses a per-thread cached
/// handle and `read_exact_at` (no open()+seek per projection); BGZF keeps the
/// virtual-position path in `read_cigar_data`.
fn read_paf_cigar_data(alignment_file: &str, offset: u64, buffer: &mut [u8]) -> Result<(), String> {
    let is_compressed = alignment_file.ends_with(".gz") || alignment_file.ends_with(".bgz");
    if is_compressed {
        return read_cigar_data(alignment_file, offset, buffer);
    }

    use std::os::unix::fs::FileExt;
    with_cached_handle(
        &PAF_HANDLE,
        alignment_file,
        |path| {
            std::fs::File::open(path)
                .map_err(|e| format!("Failed to open file '{}': {}", path, e))
        },
        |file| {
            file.read_exact_at(buffer, offset).map_err(|e| {
                format!(
                    "Failed to read {} CIGAR bytes at offset {} in '{}': {}",
                    buffer.len(),
                    offset,
                    alignment_file,
                    e
                )
            })
        },
    )
}

fn parse_cigar_to_delta(cigar: &str) -> Result<Vec<CigarOp>, OneAlnParseErr> {
    // CIGAR is pure ASCII; iterate bytes to skip UTF-8 decoding of every char.
    const MAX_CIGAR_OP_LEN: i32 = (1 << 29) - 1;
    let mut ops = Vec::new();
    let mut len: i32 = 0;
    let mut has_len = false;

    for &c in cigar.as_bytes() {
        if c.is_ascii_digit() {
            has_len = true;
            let next_len = len
                .checked_mul(10)
                .and_then(|v| v.checked_add((c - b'0') as i32))
                .ok_or_else(|| {
                    OneAlnParseErr::InvalidFormat(format!(
                        "CIGAR operation length overflows i32 in '{}'",
                        cigar
                    ))
                })?;
            if next_len > MAX_CIGAR_OP_LEN {
                return Err(OneAlnParseErr::InvalidFormat(format!(
                    "CIGAR operation length {} exceeds maximum {} in '{}'",
                    next_len, MAX_CIGAR_OP_LEN, cigar
                )));
            }
            len = next_len;
        } else {
            let op = c as char;
            if !matches!(op, '=' | 'X' | 'I' | 'D' | 'M') {
                return Err(OneAlnParseErr::InvalidFormat(format!(
                    "Invalid CIGAR operation '{}' in '{}'",
                    op, cigar
                )));
            }
            if !has_len {
                return Err(OneAlnParseErr::InvalidFormat(format!(
                    "Missing CIGAR operation length before '{}' in '{}'",
                    op, cigar
                )));
            }
            // Zero-length operation: no-op
            if len != 0 {
                ops.push(CigarOp::new(len, op));
            }
            len = 0;
            has_len = false;
        }
    }

    if has_len {
        return Err(OneAlnParseErr::InvalidFormat(format!(
            "Trailing CIGAR digits without operation in '{}'",
            cigar
        )));
    }

    Ok(ops)
}

fn calculate_gap_compressed_identity(cigar_ops: &[CigarOp]) -> f64 {
    let (matches, mismatches, insertions, deletions) =
        cigar_ops
            .iter()
            .fold((0i64, 0i64, 0i64, 0i64), |(m, mm, i, d), op| {
                let len = op.len() as i64;
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
mod tree_cache_tests {
    use super::*;

    /// Build a tree with `n` one-base intervals, so its accounted size is
    /// predictable: n * (size_of::<QueryMetadata>() + 24).
    fn tree_of(n: usize) -> Arc<BasicCOITree<QueryMetadata, u32>> {
        let intervals: Vec<Interval<QueryMetadata>> = (0..n)
            .map(|i| Interval {
                first: i as i32,
                last: i as i32,
                metadata: QueryMetadata {
                    query_id: i as u32,
                    target_start: 0,
                    target_end: 1,
                    query_start: 0,
                    query_end: 1,
                    alignment_file_index: 0,
                    strand_and_data_offset: 0,
                    data_bytes: 0,
                },
            })
            .collect();
        Arc::new(BasicCOITree::new(intervals.as_slice()))
    }

    fn bytes_for(n: usize) -> usize {
        n * (std::mem::size_of::<QueryMetadata>() + 24)
    }

    #[test]
    fn unbounded_cache_never_evicts() {
        let mut trees = TreeMap::default();
        for id in 0..10u32 {
            trees.insert(id, tree_of(50));
        }
        let cache = TreeCache::unbounded(trees);
        for id in 0..10u32 {
            assert!(cache.get(id as u64).is_some(), "target {id} was evicted");
        }
    }

    #[test]
    fn bounded_cache_evicts_least_recently_used() {
        // Room for two trees of 100 intervals, not three.
        let mut cache = TreeCache::bounded(bytes_for(250));
        cache.insert(1, tree_of(100));
        cache.insert(2, tree_of(100));
        // Touch 1 so that 2 becomes the least recently used.
        assert!(cache.get(1).is_some());
        cache.insert(3, tree_of(100));

        assert!(cache.get(1).is_some(), "recently used tree was evicted");
        assert!(cache.get(3).is_some(), "just-inserted tree was evicted");
        assert!(cache.get(2).is_none(), "LRU tree should have been evicted");
    }

    #[test]
    fn oversized_tree_is_still_served() {
        // A single tree larger than the whole capacity must not evict itself,
        // otherwise a big chromosome could never be loaded.
        let mut cache = TreeCache::bounded(bytes_for(10));
        cache.insert(7, tree_of(1000));
        assert!(cache.get(7).is_some());
    }

    #[test]
    fn reinserting_same_target_does_not_double_count() {
        let mut cache = TreeCache::bounded(bytes_for(250));
        cache.insert(1, tree_of(100));
        cache.insert(1, tree_of(100));
        cache.insert(2, tree_of(100));
        // Only two trees' worth of bytes are live, so nothing should be evicted.
        assert!(cache.get(1).is_some());
        assert!(cache.get(2).is_some());
    }

    #[test]
    fn remove_frees_capacity() {
        let mut cache = TreeCache::bounded(bytes_for(250));
        cache.insert(1, tree_of(100));
        cache.insert(2, tree_of(100));
        cache.remove(1);
        cache.insert(3, tree_of(100));
        // Removing 1 made room, so 2 must survive.
        assert!(cache.get(2).is_some(), "eviction happened despite free room");
        assert!(cache.get(3).is_some());
        assert!(cache.get(1).is_none());
    }

    #[test]
    fn shared_cache_evicts_across_namespaces() {
        // Two sub-indices share one budget: the LRU tree goes first, whoever owns it
        let key = |ns: u64, target: u64| (ns << 32) | target;
        let mut cache = TreeCache::bounded(bytes_for(250));
        cache.insert(key(0, 1), tree_of(100));
        cache.insert(key(1, 1), tree_of(100));
        cache.insert(key(2, 1), tree_of(100));
        assert!(cache.get(key(0, 1)).is_none(), "the oldest tree, of another sub-index, must go");
        assert!(cache.get(key(1, 1)).is_some());
        assert!(cache.get(key(2, 1)).is_some());
        assert!(cache.bytes <= bytes_for(250));
    }

    #[test]
    fn iter_yields_every_cached_tree() {
        let mut cache = TreeCache::bounded(bytes_for(10_000));
        for id in 0..5u64 {
            cache.insert(id, tree_of(10));
        }
        let mut ids: Vec<u64> = cache.iter().map(|(&id, _)| id).collect();
        ids.sort_unstable();
        assert_eq!(ids, vec![0, 1, 2, 3, 4]);
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

    #[test]
    fn test_zero_length_cigar_op_is_ignored() {
        let ops = parse_cigar_to_delta("10=0I5=").unwrap();
        // No merge: adjacent runs stay separate
        assert_eq!(ops, vec![CigarOp::new(10, '='), CigarOp::new(5, '=')]);
        assert_eq!(ops, parse_cigar_to_delta("10=5=").unwrap());
        assert_eq!(parse_cigar_to_delta("0D10=0X5=0I").unwrap(), ops);
        assert!(parse_cigar_to_delta("0=0I").unwrap().is_empty());
        assert!(parse_cigar_to_delta("10=0").is_err()); // trailing digits
        assert!(parse_cigar_to_delta("10=0Q").is_err()); // invalid op

        let plain = parse_cigar_to_delta("10=5=").unwrap();
        assert_eq!(
            calculate_gap_compressed_identity(&ops),
            calculate_gap_compressed_identity(&plain)
        );
        assert_eq!(calculate_gap_compressed_identity(&ops), 1.0);
        let record = (0, 15, 0, 15, Strand::Forward);
        for range in [(0, 15), (3, 12), (9, 11)] {
            assert_eq!(
                project_target_range_through_alignment(range, record, &ops),
                project_target_range_through_alignment(range, record, &plain),
                "{range:?}"
            );
        }
    }

    #[test]
    fn test_parse_cigar_to_delta_invalid() {
        assert!(parse_cigar_to_delta("10=5Q").is_err()); // Q is not a valid CIGAR operation
        assert!(parse_cigar_to_delta("=5M").is_err()); // missing length
        assert!(parse_cigar_to_delta("10=5").is_err()); // trailing digits
        assert!(parse_cigar_to_delta("600000000=").is_err()); // > 29-bit max length
        assert!(parse_cigar_to_delta("9999999999=").is_err()); // i32 overflow
    }

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

    #[test]
    fn test_invert_cigar_forward_strand() {
        let ops = vec![
            CigarOp::new(10, '='),
            CigarOp::new(5, 'I'),
            CigarOp::new(3, 'D'),
            CigarOp::new(7, 'X'),
        ];
        let inverted = invert_cigar_ops(&ops, Strand::Forward);

        // For forward strand: only I↔D swap, order unchanged
        assert_eq!(inverted.len(), 4);
        assert_eq!(inverted[0].op(), '=');
        assert_eq!(inverted[0].len(), 10);
        assert_eq!(inverted[1].op(), 'D'); // I became D
        assert_eq!(inverted[1].len(), 5);
        assert_eq!(inverted[2].op(), 'I'); // D became I
        assert_eq!(inverted[2].len(), 3);
        assert_eq!(inverted[3].op(), 'X');
        assert_eq!(inverted[3].len(), 7);
    }

    #[test]
    fn test_invert_cigar_reverse_strand() {
        let ops = vec![
            CigarOp::new(10, '='),
            CigarOp::new(5, 'I'),
            CigarOp::new(3, 'D'),
        ];
        let inverted = invert_cigar_ops(&ops, Strand::Reverse);

        // For reverse strand: I↔D swap AND array reversal
        assert_eq!(inverted.len(), 3);
        assert_eq!(inverted[0].op(), 'I'); // Was D at end, now I at start
        assert_eq!(inverted[0].len(), 3);
        assert_eq!(inverted[1].op(), 'D'); // Was I in middle, now D
        assert_eq!(inverted[1].len(), 5);
        assert_eq!(inverted[2].op(), '='); // Was = at start, now at end
        assert_eq!(inverted[2].len(), 10);
    }

    #[test]
    fn test_invert_cigar_empty() {
        let ops: Vec<CigarOp> = vec![];
        let inverted_forward = invert_cigar_ops(&ops, Strand::Forward);
        let inverted_reverse = invert_cigar_ops(&ops, Strand::Reverse);

        assert_eq!(inverted_forward.len(), 0);
        assert_eq!(inverted_reverse.len(), 0);
    }

    #[test]
    fn test_invert_cigar_matches_only() {
        let ops = vec![CigarOp::new(100, '='), CigarOp::new(50, 'X')];
        let inverted_forward = invert_cigar_ops(&ops, Strand::Forward);
        let inverted_reverse = invert_cigar_ops(&ops, Strand::Reverse);

        // Forward: no change (no I or D)
        assert_eq!(inverted_forward[0].op(), '=');
        assert_eq!(inverted_forward[1].op(), 'X');

        // Reverse: just reversed order
        assert_eq!(inverted_reverse[0].op(), 'X');
        assert_eq!(inverted_reverse[1].op(), '=');
    }

    #[test]
    fn test_index_rejects_paf_record_without_cigar() {
        // A PAF record without cg:Z cannot be projected. The index must say so
        // instead of panicking inside a worker at query time.
        let paf = "q\t100\t0\t100\t+\tt\t100\t0\t100\t60\t100\t255\n".to_string();
        let mut seq_index = SequenceIndex::new();
        let records = parse_paf(BufReader::new(paf.as_bytes()), &mut seq_index).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].data_bytes, 0, "fixture must have no CIGAR");

        let records_by_file = vec![(records, "test.paf".to_string())];
        let err = match Impg::from_multi_alignment_records(&records_by_file, seq_index, None, true) {
            Ok(_) => panic!("a PAF record without cg:Z must be rejected at index time"),
            Err(e) => e,
        };
        let msg = err.to_string();
        assert!(
            msg.contains("without a CIGAR"),
            "unexpected error message: {msg}"
        );
        assert!(msg.contains("test.paf"), "error must name the file: {msg}");
    }

    #[test]
    fn test_index_accepts_paf_record_with_cigar() {
        let paf = "q\t100\t0\t100\t+\tt\t100\t0\t100\t60\t100\t255\tcg:Z:100=\n".to_string();
        let mut seq_index = SequenceIndex::new();
        let records = parse_paf(BufReader::new(paf.as_bytes()), &mut seq_index).unwrap();
        let records_by_file = vec![(records, "test.paf".to_string())];
        Impg::from_multi_alignment_records(&records_by_file, seq_index, None, true)
            .expect("a PAF record with a CIGAR must index");
    }

    #[test]
    fn test_index_rejects_sequence_longer_than_i32() {
        // Coordinates are i32, so a longer sequence would wrap silently.
        let mut seq_index = SequenceIndex::new();
        let long = i32::MAX as usize + 1;
        seq_index.get_or_insert_id("chrBig", Some(long));
        seq_index.get_or_insert_id("chrSmall", Some(1000));

        let err = match Impg::from_multi_alignment_records(&[], seq_index, None, true) {
            Ok(_) => panic!("a sequence above i32::MAX bp must be rejected"),
            Err(e) => e,
        };
        let msg = err.to_string();
        assert!(msg.contains("chrBig"), "error must name the sequence: {msg}");
        assert!(
            msg.contains(&long.to_string()),
            "error must give the length: {msg}"
        );
    }

    #[test]
    fn test_index_rejects_coordinate_above_i32() {
        let mut seq_index = SequenceIndex::new();
        seq_index.get_or_insert_id("q", Some(1000));
        seq_index.get_or_insert_id("t", Some(1000));

        let record = AlignmentRecord {
            query_id: 0,
            query_start: 0,
            query_end: i32::MAX as usize + 10,
            target_id: 1,
            target_start: 0,
            target_end: 100,
            strand_and_data_offset: 0,
            data_bytes: 10,
        };
        let records_by_file = vec![(vec![record], "test.paf".to_string())];
        let err = match Impg::from_multi_alignment_records(&records_by_file, seq_index, None, true) {
            Ok(_) => panic!("a coordinate above i32::MAX must be rejected"),
            Err(e) => e,
        };
        assert!(
            err.to_string().contains("above the maximum supported value"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_in_memory_index_does_not_reload_trees_from_disk() {
        // `impg stats` walks every tree twice and evicts as it goes. An
        // in-memory index has no file to load a tree back from, so the second
        // walk must not try to open an empty path.
        let paf = "q\t100\t0\t100\t+\tt\t100\t0\t100\t60\t100\t255\tcg:Z:100=\n".to_string();
        let mut seq_index = SequenceIndex::new();
        let records = parse_paf(BufReader::new(paf.as_bytes()), &mut seq_index).unwrap();
        let records_by_file = vec![(records, "test.paf".to_string())];
        let impg =
            Impg::from_multi_alignment_records(&records_by_file, seq_index, None, true).unwrap();

        let target_id = *impg.target_ids().first().expect("one target sequence");
        assert!(
            !impg.can_reload_tree(target_id),
            "an in-memory index must report that its trees cannot be reloaded"
        );
        assert!(impg.get_or_load_tree(target_id).is_some());
        impg.remove_cached_tree(target_id);
        assert!(
            impg.get_or_load_tree(target_id).is_none(),
            "a tree evicted from an in-memory index must come back as None, not as a panic"
        );
    }

    #[test]
    fn in_memory_index_serves_two_walks_when_trees_are_not_evicted() {
        // The behaviour `impg stats` needs: when the index cannot reload a
        // tree, both walks must still see every interval.
        let paf = "q\t1000\t0\t500\t+\tt\t1000\t0\t500\t500\t500\t60\tcg:Z:500=\n\
                   t\t1000\t100\t600\t+\tr\t1000\t100\t600\t500\t500\t60\tcg:Z:500=\n"
            .to_string();
        let mut seq_index = SequenceIndex::new();
        let records = parse_paf(BufReader::new(paf.as_bytes()), &mut seq_index).unwrap();
        let records_by_file = vec![(records, "test.paf".to_string())];
        let impg =
            Impg::from_multi_alignment_records(&records_by_file, seq_index, None, true).unwrap();

        let mut first_walk = 0usize;
        let mut second_walk = 0usize;
        for id in impg.target_ids() {
            if let Some(tree) = impg.get_or_load_tree(id) {
                first_walk += tree.len();
                if impg.can_reload_tree(id) {
                    impg.remove_cached_tree(id);
                }
            }
        }
        for id in impg.target_ids() {
            if let Some(tree) = impg.get_or_load_tree(id) {
                second_walk += tree.iter().count();
                if impg.can_reload_tree(id) {
                    impg.remove_cached_tree(id);
                }
            }
        }

        assert_eq!(first_walk, 4, "two PAF records indexed bidirectionally");
        assert_eq!(
            second_walk, first_walk,
            "the second walk must see the same intervals as the first"
        );
    }

    #[test]
    fn index_header_records_the_alignment_file_list() {
        let temp = tempfile::TempDir::new().unwrap();
        let a = temp.path().join("a.paf");
        let b = temp.path().join("b.paf");
        std::fs::write(&a, "x".repeat(100)).unwrap();
        std::fs::write(&b, "y".repeat(200)).unwrap();
        let names = vec![a.to_string_lossy().into_owned(), b.to_string_lossy().into_owned()];

        let header = IndexHeaderInfo {
            alignment_files: names
                .iter()
                .map(|p| AlignmentFileEntry::from_path(p))
                .collect(),
            bidirectional: true,
        };
        assert_eq!(header.alignment_files[0].size, Some(100));
        assert_eq!(header.alignment_files[1].size, Some(200));

        // The same list in the same order is accepted, even when a path is
        // spelled differently.
        validate_alignment_file_list(&header.alignment_files, &names, "idx").unwrap();
        let respelled = vec![
            format!("{}/./a.paf", temp.path().display()),
            temp.path().join("..").join(temp.path().file_name().unwrap()).join("b.paf").to_string_lossy().into_owned(),
        ];
        validate_alignment_file_list(&header.alignment_files, &respelled, "idx")
            .expect("a different spelling of the same files must be accepted");

        // A reordered list must be refused: every file index still resolves, so
        // the query would read a CIGAR offset from the wrong file.
        let mut reordered = names.clone();
        reordered.reverse();
        let err = validate_alignment_file_list(&header.alignment_files, &reordered, "idx")
            .expect_err("a reordered file list must be refused");
        let msg = err.to_string();
        assert!(msg.contains("different alignment file list"), "{msg}");
        assert!(msg.contains("Entry 0"), "{msg}");

        // A different list length must be refused too.
        let err = validate_alignment_file_list(&header.alignment_files, &names[..1], "idx")
            .expect_err("a shorter file list must be refused");
        assert!(err.to_string().contains("was built from 2"), "{err}");

        // A file of a different size in the same position must be refused.
        std::fs::write(&a, "z".repeat(101)).unwrap();
        assert!(
            validate_alignment_file_list(&header.alignment_files, &names, "idx").is_err(),
            "a file that changed size must be refused"
        );
    }

    #[test]
    fn serialized_index_carries_the_file_list_and_refuses_a_mismatch() {
        let temp = tempfile::TempDir::new().unwrap();
        let paf = temp.path().join("a.paf");
        std::fs::write(
            &paf,
            "q\t100\t0\t100\t+\tt\t100\t0\t100\t60\t100\t255\tcg:Z:100=\n",
        )
        .unwrap();
        let files = vec![paf.to_string_lossy().into_owned()];

        let mut seq_index = SequenceIndex::new();
        let content = std::fs::read_to_string(&paf).unwrap();
        let records = parse_paf(BufReader::new(content.as_bytes()), &mut seq_index).unwrap();
        let impg =
            Impg::from_multi_alignment_records(&[(records, files[0].clone())], seq_index, None, true)
                .unwrap();

        let mut buffer = std::io::Cursor::new(Vec::new());
        impg.serialize_with_forest_map(&mut buffer).unwrap();
        let bytes = buffer.into_inner();

        // The same list loads.
        let loaded = Impg::load_from_file(
            std::io::Cursor::new(bytes.clone()),
            &files,
            "idx".to_string(),
            None,
        )
        .expect("the same file list must load");
        assert!(loaded.built_bidirectional);

        // A different file list is refused rather than reading wrong offsets.
        let mut other = files.clone();
        other.push(temp.path().join("b.paf").to_string_lossy().into_owned());
        std::fs::write(&other[1], "q\n").unwrap();
        let err = match Impg::load_from_file(
            std::io::Cursor::new(bytes),
            &other,
            "idx".to_string(),
            None,
        ) {
            Ok(_) => panic!("a different file list must be refused"),
            Err(e) => e,
        };
        assert!(
            err.to_string().contains("different alignment file list")
                || err.to_string().contains("was built from 1"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn reciprocal_records_collapse_to_one_edge() {
        // The same alignment written in both directions describes one edge. It
        // used to be indexed twice, so every projection appeared twice.
        let paf = "A\t1000\t0\t500\t+\tB\t1000\t0\t500\t500\t500\t60\tcg:Z:500=\n\
                   B\t1000\t0\t500\t+\tA\t1000\t0\t500\t500\t500\t60\tcg:Z:500=\n"
            .to_string();
        let mut seq_index = SequenceIndex::new();
        let records = parse_paf(BufReader::new(paf.as_bytes()), &mut seq_index).unwrap();
        let impg =
            Impg::from_multi_alignment_records(&[(records, "r.paf".to_string())], seq_index, None, true)
                .unwrap();

        let mut total = 0usize;
        for id in impg.target_ids() {
            if let Some(tree) = impg.get_or_load_tree(id) {
                total += tree.len();
            }
        }
        assert_eq!(
            total, 2,
            "one edge between A and B must be indexed once per direction, not twice"
        );

        // A single record still indexes both directions.
        let paf = "A\t1000\t0\t500\t+\tB\t1000\t0\t500\t500\t500\t60\tcg:Z:500=\n".to_string();
        let mut seq_index = SequenceIndex::new();
        let records = parse_paf(BufReader::new(paf.as_bytes()), &mut seq_index).unwrap();
        let impg =
            Impg::from_multi_alignment_records(&[(records, "r.paf".to_string())], seq_index, None, true)
                .unwrap();
        let mut total = 0usize;
        for id in impg.target_ids() {
            if let Some(tree) = impg.get_or_load_tree(id) {
                total += tree.len();
            }
        }
        assert_eq!(total, 2, "one record must keep both directions");
    }
    #[test]
    fn sorted_ranges_snap_keeps_novel_tail_past_an_inner_range() {
        // Regression: snapping `end` onto `ranges[i].0` truncated the inserted
        // range and silently dropped `110..140`.
        let mut ranges = SortedRanges::with_distance(1000, 50);
        ranges.insert((105, 110));

        let novel = ranges.insert_with_boundary_snap((10, 140), 50);

        assert_eq!(
            novel,
            vec![(10, 105), (110, 140)],
            "the tail of the inserted range must survive"
        );
        assert_eq!(ranges.iter().copied().collect::<Vec<_>>(), vec![(10, 140)]);
    }

    #[test]
    fn sorted_ranges_snap_never_shrinks_the_inserted_range() {
        let mut ranges = SortedRanges::with_distance(1000, 10);
        ranges.insert((105, 110));

        let novel = ranges.insert_with_boundary_snap((50, 200), 10);

        assert_eq!(novel, vec![(50, 105), (110, 200)]);
        assert_eq!(ranges.iter().copied().collect::<Vec<_>>(), vec![(50, 200)]);
    }

    #[test]
    fn sorted_ranges_snap_still_closes_small_gaps() {
        let mut ranges = SortedRanges::with_distance(1000, 20);
        ranges.insert((100, 200));

        // 5 bp gap on the left: `start` snaps onto the previous end, so the
        // jitter is absorbed but the novel tail stays novel.
        let novel = ranges.insert_with_boundary_snap((205, 300), 20);

        assert_eq!(novel, vec![(200, 300)]);
        assert_eq!(ranges.iter().copied().collect::<Vec<_>>(), vec![(100, 300)]);
    }

    #[test]
    fn sorted_ranges_plain_insert_ignores_snapping() {
        // `insert` is used for user-requested ranges: `SortedRanges::new()`
        // carries `min_distance == 0`, so nothing is ever absorbed.
        let mut ranges = SortedRanges::new();
        ranges.insert((0, 5));

        let novel = ranges.insert((10, 40));

        assert_eq!(novel, vec![(10, 40)]);
    }

    fn test_dna_seq(len: usize, seed: u64) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                bases[((state >> 33) % 4) as usize]
            })
            .collect()
    }

    fn test_revcomp_seq(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => *other,
            })
            .collect()
    }

    #[test]
    fn test_calculate_gap_compressed_identity_large_cigar_no_overflow() {
        // Total bases > i32::MAX (2.5e9) must not overflow with overflow-checks=true.
        let ops = vec![
            CigarOp::new(500_000_000, '='),
            CigarOp::new(500_000_000, '='),
            CigarOp::new(500_000_000, '='),
            CigarOp::new(500_000_000, '='),
            CigarOp::new(500_000_000, 'X'),
        ];
        let id = calculate_gap_compressed_identity(&ops);
        assert!((id - 0.8).abs() < 1e-12, "unexpected identity: {id}");
    }

    #[test]
    fn test_populate_cigar_cache_distinguishes_forward_and_reversed() {
        let temp = tempfile::TempDir::new().unwrap();
        let paf_path = temp.path().join("asym.paf");
        // Asymmetric deletion (20=5D20=): indexed bidirectionally, seqB's tree holds
        // the forward entry (REVERSED_BIT=0 -> 20=5D20=) and seqA's tree holds the
        // reversed entry (REVERSED_BIT=1 -> 20=5I20=) at the exact same file index
        // and CIGAR byte offset.
        std::fs::write(
            &paf_path,
            "seqA\t200\t0\t40\t+\tseqB\t200\t100\t145\t40\t45\t60\tcg:Z:20=5D20=\n",
        )
        .unwrap();

        let paf_str = paf_path.to_string_lossy().into_owned();
        let mut seq_index = SequenceIndex::new();
        let file = std::fs::File::open(&paf_path).unwrap();
        let records = parse_paf(BufReader::new(file), &mut seq_index).unwrap();
        let impg =
            Impg::from_multi_alignment_records(&[(records, paf_str)], seq_index, None, true)
                .unwrap();

        let seq_a_id = impg.seq_index.get_id("seqA").unwrap();
        let seq_b_id = impg.seq_index.get_id("seqB").unwrap();
        let mut cache = FxHashMap::default();
        impg.populate_cigar_cache(seq_b_id, 100, 145, None, None, &mut cache);
        impg.populate_cigar_cache(seq_a_id, 0, 40, None, None, &mut cache);
        assert_eq!(
            cache.len(),
            2,
            "forward and reversed entries sharing an alignment offset must have distinct cache keys"
        );

        let to_comparable = |res: Vec<AdjustedInterval>| {
            res.into_iter()
                .map(|(q, ops, t)| (q.metadata, q.first, q.last, ops, t.metadata, t.first, t.last))
                .collect::<Vec<_>>()
        };

        for (tid, start, end) in [(seq_b_id, 105, 140), (seq_a_id, 5, 35)] {
            let uncached =
                to_comparable(impg.query(tid, start, end, true, None, None, false));
            let cached = to_comparable(impg.query_with_cache(
                tid, start, end, true, None, None, &cache, false,
            ));
            assert_eq!(cached.len(), 2);
            assert_eq!(cached, uncached, "cached projection diverged for ({tid}, {start}..{end})");
        }
    }

    #[test]
    fn test_tpa_process_tracepoints_stale_differences_and_reversed_entry() {
        let temp = tempfile::TempDir::new().unwrap();
        let fa_path = temp.path().join("seqs.fa");
        // Two sequences of equal length with 1 mismatch in the middle
        let seq_q = test_dna_seq(120, 7);
        let mut seq_t = seq_q.clone();
        seq_t[60] = if seq_t[60] == b'A' { b'C' } else { b'A' };
        std::fs::write(
            &fa_path,
            format!(
                ">q\n{}\n>t\n{}\n",
                String::from_utf8_lossy(&seq_q),
                String::from_utf8_lossy(&seq_t)
            ),
        )
        .unwrap();

        let fa_str = fa_path.to_string_lossy().into_owned();
        let seq_index_uni = UnifiedSequenceIndex::from_files(&[fa_str.clone()]).unwrap();

        let mut seq_index = SequenceIndex::new();
        let q_id = seq_index.get_or_insert_id("q", Some(120));
        let t_id = seq_index.get_or_insert_id("t", Some(120));
        let rec = AlignmentRecord {
            query_id: q_id,
            query_start: 0,
            query_end: 120,
            target_id: t_id,
            target_start: 0,
            target_end: 120,
            strand_and_data_offset: 0,
            data_bytes: 0,
        };
        let impg = Impg::from_multi_alignment_records(
            &[(vec![rec], "dummy.1aln".to_string())],
            seq_index,
            Some(&[fa_str]),
            true,
        )
        .unwrap();

        // Simulate stale metadata claiming differences == 0 even though base 60 differs.
        let aln = crate::onealn::OneAlnAlignment {
            query_name: "q".to_string(),
            query_length: 120,
            query_contig_start: 0,
            query_contig_end: 120,
            query_contig_offset: 0,
            target_name: "t".to_string(),
            target_length: 120,
            target_contig_start: 0,
            target_contig_end: 120,
            target_contig_offset: 0,
            target_contig_len: 120,
            strand: '+',
            differences: 0,
            tracepoints: vec![60, 60],
            mode_data: crate::onealn::TracepointModeData::Fastga {
                diffs: vec![0, 0],
                trace_spacing: 60,
            },
        };

        // Reversed tree entry (`t` was indexed under `q`'s tree with REVERSED_BIT).
        let rev_meta = QueryMetadata {
            query_id: t_id,
            target_start: 0,
            target_end: 120,
            query_start: 0,
            query_end: 120,
            alignment_file_index: 0,
            strand_and_data_offset: QueryMetadata::REVERSED_BIT,
            data_bytes: 0,
        };

        let ops = impg.process_tracepoints_data(&aln, &rev_meta, q_id, &seq_index_uni);
        assert!(
            ops.iter().any(|op| op.op() == 'X'),
            "stale differences == 0 must not blindly emit matches when bases differ: {ops:?}"
        );
    }

    #[test]
    fn test_tpa_reverse_strand_asymmetric_indels_matches_paf_all_modes_and_approx() {
        let temp = tempfile::TempDir::new().unwrap();
        let dir = temp.path();

        // Build `t` (len 1200) and `q_rev` aligning on '-' strand to t[50..1130]
        // with asymmetric deletion (30D), mismatch (1X), and insertion (20I):
        //   PAF CIGAR (fwd(t[50..1130]) vs revcomp(q_rev)):
        //   200=30D170=1X179=20I500=
        let t = test_dna_seq(1200, 99);
        let ins = test_dna_seq(20, 123);
        let mut aligned_q_fwd = Vec::with_capacity(1070);
        aligned_q_fwd.extend_from_slice(&t[50..250]); // 200=
        // 30D skips t[250..280]
        aligned_q_fwd.extend_from_slice(&t[280..450]); // 170=
        aligned_q_fwd.push(if t[450] == b'A' { b'C' } else { b'A' }); // 1X at t[450]
        aligned_q_fwd.extend_from_slice(&t[451..630]); // 179=
        aligned_q_fwd.extend_from_slice(&ins); // 20I
        aligned_q_fwd.extend_from_slice(&t[630..1130]); // 500=
        assert_eq!(aligned_q_fwd.len(), 1070);

        let q_rev = test_revcomp_seq(&aligned_q_fwd);

        let fa_path = dir.join("seqs.fa");
        std::fs::write(
            &fa_path,
            format!(
                ">q_rev\n{}\n>t\n{}\n",
                String::from_utf8_lossy(&q_rev),
                String::from_utf8_lossy(&t),
            ),
        )
        .unwrap();
        let fa_str = fa_path.to_string_lossy().into_owned();
        let seq_index_uni = UnifiedSequenceIndex::from_files(&[fa_str.clone()]).unwrap();

        let paf_path = dir.join("aln.paf");
        std::fs::write(
            &paf_path,
            "q_rev\t1070\t0\t1070\t-\tt\t1200\t50\t1130\t1049\t1100\t60\tcg:Z:200=30D170=1X179=20I500=\n",
        )
        .unwrap();
        let paf_str = paf_path.to_string_lossy().into_owned();

        let mut paf_seq_idx = SequenceIndex::new();
        let paf_records =
            parse_paf(BufReader::new(std::fs::File::open(&paf_path).unwrap()), &mut paf_seq_idx)
                .unwrap();
        let paf_impg =
            Impg::from_multi_alignment_records(&[(paf_records, paf_str)], paf_seq_idx, None, true)
                .unwrap();

        let normalize_results =
            |impg: &Impg, res: Vec<AdjustedInterval>| -> Vec<(String, i32, i32, i32, i32)> {
                let mut rows: Vec<_> = res
                    .into_iter()
                    .map(|(q_iv, _, t_iv)| {
                        (
                            impg.seq_index.get_name(q_iv.metadata).unwrap().to_string(),
                            q_iv.first,
                            q_iv.last,
                            t_iv.first,
                            t_iv.last,
                        )
                    })
                    .collect();
                rows.sort();
                rows
            };

        use tpa::TracepointType as T;
        for tp_type in [T::Standard, T::Variable, T::Mixed, T::Fastga, T::FastgaNoDiff] {
            let tpa_path = dir.join(format!("{tp_type:?}.tpa"));
            let tpa_str = tpa_path.to_string_lossy().into_owned();
            tpa::paf_to_tpa(
                paf_path.to_str().unwrap(),
                tpa_path.to_str().unwrap(),
                tpa::CompressionConfig {
                    tp_type,
                    max_complexity: 50,
                    ..tpa::CompressionConfig::default()
                },
            )
            .unwrap();

            let parser = crate::tpa_parser::TpaParser::new(tpa_str.clone()).unwrap();
            let mut tpa_seq_idx = SequenceIndex::new();
            let tpa_records = parser.parse_alignments(&mut tpa_seq_idx, 1).unwrap();
            let tpa_impg = Impg::from_multi_alignment_records(
                &[(tpa_records, tpa_str)],
                tpa_seq_idx,
                Some(std::slice::from_ref(&fa_str)),
                true,
            )
            .unwrap();

            // Exact sub-range queries spanning start, middle, and trailing segments
            // on both t (forward entry) and q_rev (REVERSED_BIT entry).
            for (seq_name, start, end) in [
                ("t", 50, 1130),
                ("t", 100, 220),
                ("t", 320, 600),
                ("t", 700, 1100),
                ("t", 1000, 1130),
                ("q_rev", 0, 1070),
                ("q_rev", 50, 450),
                ("q_rev", 620, 820),
                ("q_rev", 950, 1070),
            ] {
                let paf_tid = paf_impg.seq_index.get_id(seq_name).unwrap();
                let tpa_tid = tpa_impg.seq_index.get_id(seq_name).unwrap();

                let expected = normalize_results(
                    &paf_impg,
                    paf_impg.query(paf_tid, start, end, false, None, None, false),
                );
                let actual = normalize_results(
                    &tpa_impg,
                    tpa_impg.query(tpa_tid, start, end, false, None, Some(&seq_index_uni), false),
                );
                assert_eq!(
                    actual, expected,
                    "{tp_type:?} projection diverged for {seq_name}:{start}-{end}"
                );

                // Approximate mode must preserve self-interval forward orientation
                // (`qs < qe`), reverse-strand projected orientation (`qs > qe`), and
                // exact full-span endpoint bounds.
                let approx = normalize_results(
                    &tpa_impg,
                    tpa_impg.query(tpa_tid, start, end, false, None, None, true),
                );
                assert_eq!(approx.len(), 2, "{tp_type:?} approx count for {seq_name}:{start}-{end}");
                for (q_name, qs, qe, ts, te) in &approx {
                    if q_name == seq_name {
                        assert!(qs < qe, "{tp_type:?} approx self query inverted: {qs}..{qe}");
                    } else {
                        assert!(
                            qs > qe,
                            "{tp_type:?} approx reverse projection expected qs > qe, got {qs}..{qe}"
                        );
                    }
                    assert!(ts < te, "{tp_type:?} approx target inverted: {ts}..{te}");
                }
                if (start, end) == (50, 1130) || (start, end) == (0, 1070) {
                    assert_eq!(approx, expected, "{tp_type:?} approx full-span endpoints mismatch");
                }
            }
        }
    }
}
