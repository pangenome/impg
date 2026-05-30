//! Dump a loaded `SyngIndex` to GFA (1.0 with P lines by default,
//! 1.1 with W lines if `--gfa-version 1.1`).
//!
//! Output: one S per syncmer (node id 1..N, canonical DNA from the
//! KmerHash) plus one S per inter-syncmer gap (id N+1, N+2, …). When
//! the bp offset between two consecutive syncmers in a path exceeds the
//! syncmer length, a gap segment of `offset - syncmer_length` bp is
//! inserted between them. Gap bases come from `--sequence-files` if
//! provided, otherwise the gap is filled with `N`s and a warning is
//! emitted. Gap segments are interned by exact sequence plus local
//! signed-syncmer context, so identical terminal/inter-syncmer DNA shared
//! by the same local graph context is one node, without collapsing unrelated
//! repeated sequence elsewhere. High-frequency syncmers are filtered out of
//! the raw syng walk before materialization, and the original sequence is
//! used to bridge retained anchors. The CLI default is blunt mode, which
//! materializes exact source-spelling 0M paths directly; raw mode is the
//! explicit syng-native overlap graph.

use std::io::{self, BufWriter, Write};
use std::time::Instant;

use log::{debug, info, warn};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use crate::syng::{Anchor, GbwtPathStart, HomologousIntervalWithAnchors, SyngIndex};
use crate::syng_ffi;

/// Default fraction of query-local syncmer nodes filtered out of the raw syng
/// topology. This mirrors the minimizer-style seed filter: only the worst
/// 0.05% of used nodes are treated as repetitive glue.
pub const DEFAULT_GFA_MASK_TOP_FRACTION: f64 = 0.0005;
/// Default minimum run of consecutive shared syncmers required before a local
/// syncmer node is allowed to act as graph glue. This mirrors the default
/// bounded-walk seed length used during syng query discovery.
pub const DEFAULT_GFA_MIN_SHARED_RUN: usize = crate::syng::DEFAULT_WALK_SEED_ANCHORS;
/// Clone rare repeated-copy contexts for otherwise dominant local syncmers.
/// This catches single-syncmer loops where a mostly single-copy node appears
/// a second time in a selected path and glues two paralogous contexts together.
pub const DEFAULT_GFA_LOCAL_REPEAT_MAX_MINOR: u32 = 2;
pub const DEFAULT_GFA_LOCAL_REPEAT_MIN_DOMINANT_FRACTION: f64 = 0.80;
/// Default minimum N-run length that splits syng GFA paths when N cutting is
/// enabled. A value of 1 means any ambiguous base breaks the emitted path.
pub const DEFAULT_GFA_CUT_N_MIN_RUN: usize = 1;

/// Frequency-aware syncmer node sharing policy for syng GFA materialization.
///
/// High-frequency syncmers are removed from the raw syng topology before
/// bluntification and the retained anchors are bridged by sequence. Rare
/// repeated local contexts are still emitted as per-occurrence clones so
/// ordinary variation is not erased while obvious single-copy glue is split.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SyngGfaFrequencyMask {
    pub drop_top_fraction: f64,
    pub max_occurrences: Option<u32>,
    pub min_shared_run: usize,
    pub min_sequence_span_bp: usize,
    pub local_repeat_max_minor: u32,
    pub local_repeat_min_dominant_fraction: f64,
    pub cut_n_gaps: bool,
    pub cut_n_min_run: usize,
}

impl SyngGfaFrequencyMask {
    pub fn disabled() -> Self {
        Self {
            drop_top_fraction: 0.0,
            max_occurrences: None,
            min_shared_run: 1,
            min_sequence_span_bp: 0,
            local_repeat_max_minor: 0,
            local_repeat_min_dominant_fraction: DEFAULT_GFA_LOCAL_REPEAT_MIN_DOMINANT_FRACTION,
            cut_n_gaps: false,
            cut_n_min_run: DEFAULT_GFA_CUT_N_MIN_RUN,
        }
    }

    pub fn local_default() -> Self {
        Self {
            drop_top_fraction: DEFAULT_GFA_MASK_TOP_FRACTION,
            max_occurrences: None,
            min_shared_run: DEFAULT_GFA_MIN_SHARED_RUN,
            min_sequence_span_bp: 0,
            local_repeat_max_minor: DEFAULT_GFA_LOCAL_REPEAT_MAX_MINOR,
            local_repeat_min_dominant_fraction: DEFAULT_GFA_LOCAL_REPEAT_MIN_DOMINANT_FRACTION,
            cut_n_gaps: false,
            cut_n_min_run: DEFAULT_GFA_CUT_N_MIN_RUN,
        }
    }

    pub fn enabled(self) -> bool {
        self.drop_top_fraction > 0.0
            || self.max_occurrences.is_some()
            || self.min_shared_run > 1
            || self.min_sequence_span_bp > 0
            || self.local_repeat_max_minor > 0
            || self.cut_n_gaps
    }

    fn frequency_filter_enabled(self) -> bool {
        self.drop_top_fraction > 0.0 || self.max_occurrences.is_some()
    }

    fn shared_run_filter_enabled(self) -> bool {
        self.min_shared_run > 1
    }

    fn sequence_context_filter_enabled(self) -> bool {
        self.min_sequence_span_bp > 0
    }

    fn shared_context_filter_enabled(self) -> bool {
        self.shared_run_filter_enabled() || self.sequence_context_filter_enabled()
    }
}

impl Default for SyngGfaFrequencyMask {
    fn default() -> Self {
        Self::disabled()
    }
}

/// A source path interval to render as one GFA path.
pub struct SyngGfaPathRange {
    pub path_idx: usize,
    pub name: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
}

/// Syng GFA graph shape.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SyngGfaMode {
    /// Native syng overlap graph, preserving nonzero link overlaps.
    Raw,
    /// Blunt graph produced by pangenome/bluntg; all link/path overlaps become 0M.
    Blunt,
}

impl SyngGfaMode {
    pub fn parse(s: &str) -> io::Result<Self> {
        match s.trim().replace('_', "-").to_ascii_lowercase().as_str() {
            "raw" | "syng:raw" | "syng-raw" | "syng-native:raw" => Ok(Self::Raw),
            "blunt" | "bluntg" | "syng" | "syng:blunt" | "syng:bluntg" | "syng-blunt"
            | "syng-bluntg" | "syng-native" | "syng-native:blunt" | "syng-native:bluntg" => {
                Ok(Self::Blunt)
            }
            other => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("unknown syng GFA mode '{other}' (expected raw or blunt)"),
            )),
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Raw => "raw",
            Self::Blunt => "blunt",
        }
    }
}

/// Which GFA spec version to emit.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GfaVersion {
    /// Strict GFA v1.0 — H/S/L/P only.
    V1_0,
    /// GFA v1.1 — H/S/L/W (paths emitted as walks with PanSN metadata).
    V1_1,
}

impl GfaVersion {
    pub fn parse(s: &str) -> io::Result<Self> {
        match s {
            "1.0" | "v1.0" | "1" => Ok(GfaVersion::V1_0),
            "1.1" | "v1.1" => Ok(GfaVersion::V1_1),
            other => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("unknown GFA version '{other}' (expected 1.0 or 1.1)"),
            )),
        }
    }

    fn header_tag(self) -> &'static str {
        match self {
            GfaVersion::V1_0 => "1.0",
            GfaVersion::V1_1 => "1.1",
        }
    }
}

/// Parsed PanSN-style path name.
struct PansnParts<'a> {
    sample: &'a str,
    hap: &'a str,
    contig: &'a str,
    start: u64,
    end: u64,
}

/// Try to parse `sample#hap#contig:start-end` or `sample#contig:start-end`.
///
/// Returns `None` if the name doesn't match. The hap index defaults to "0"
/// when only two `#`-fields are present (the chm13/grch38 convention used by
/// HPRC PanSN).
fn parse_pansn(name: &str) -> Option<PansnParts<'_>> {
    let hash_parts: Vec<&str> = name.split('#').collect();
    let (sample, hap, locus) = match hash_parts.as_slice() {
        [sample, locus] => (*sample, "0", *locus),
        [sample, hap, locus] => (*sample, *hap, *locus),
        _ => return None,
    };
    let (contig, range) = locus.rsplit_once(':')?;
    let (start_s, end_s) = range.split_once('-')?;
    let start: u64 = start_s.parse().ok()?;
    let end: u64 = end_s.parse().ok()?;
    Some(PansnParts {
        sample,
        hap,
        contig,
        start,
        end,
    })
}

/// One step in a P/W walk: either a syncmer (signed by strand) or an
/// inserted gap segment (always on the forward strand).
#[derive(Clone)]
enum PathStep {
    Syncmer(i32),
    Segment { id: String, sign: char },
    Gap(String),
}

/// One inserted gap segment ready to be written as an S line.
struct GapSegment {
    id: String,
    seq: Vec<u8>,
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
enum GapContext {
    Prefix {
        right: i32,
    },
    Between {
        left: i32,
        right: i32,
    },
    Suffix {
        left: i32,
    },
    PathOnly {
        path_idx: usize,
        start: u64,
        end: u64,
        strand: char,
    },
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct GapKey {
    context: GapContext,
    seq: Vec<u8>,
}

struct GapInterner {
    n_nodes: usize,
    next_gap_id: u64,
    keys: FxHashMap<GapKey, String>,
    segments: Vec<GapSegment>,
    total_segment_bp: u64,
}

impl GapInterner {
    fn new(n_nodes: usize) -> Self {
        Self {
            n_nodes,
            next_gap_id: 0,
            keys: FxHashMap::default(),
            segments: Vec::new(),
            total_segment_bp: 0,
        }
    }

    fn next_segment_id(&mut self) -> String {
        let id = format!("{}", self.n_nodes as u64 + 1 + self.next_gap_id);
        self.next_gap_id += 1;
        id
    }

    fn intern(&mut self, context: GapContext, seq: Vec<u8>) -> String {
        let key = GapKey { context, seq };
        if let Some(id) = self.keys.get(&key) {
            return id.clone();
        }

        let id = self.next_segment_id();
        self.total_segment_bp += key.seq.len() as u64;
        self.segments.push(GapSegment {
            id: id.clone(),
            seq: key.seq.clone(),
        });
        self.keys.insert(key, id.clone());
        id
    }

    fn len(&self) -> usize {
        self.segments.len()
    }

    fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }

    fn total_segment_bp(&self) -> u64 {
        self.total_segment_bp
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
enum ExactSegmentKey {
    Gap(GapKey),
    SyncmerSlice {
        signed: i32,
        left_trim: u32,
        right_trim: u32,
    },
    PrivateSyncmerSlice {
        path_name: String,
        occurrence: u32,
        signed: i32,
        left_trim: u32,
        right_trim: u32,
    },
}

struct ExactSegmentInterner {
    n_nodes: usize,
    next_extra_id: u64,
    keys: FxHashMap<ExactSegmentKey, String>,
    segments: Vec<GapSegment>,
    gap_segments: usize,
    total_gap_segment_bp: u64,
}

impl ExactSegmentInterner {
    fn new(n_nodes: usize) -> Self {
        Self {
            n_nodes,
            next_extra_id: 0,
            keys: FxHashMap::default(),
            segments: Vec::new(),
            gap_segments: 0,
            total_gap_segment_bp: 0,
        }
    }

    fn next_segment_id(&mut self) -> String {
        let id = format!("{}", self.n_nodes as u64 + 1 + self.next_extra_id);
        self.next_extra_id += 1;
        id
    }

    fn intern(&mut self, key: ExactSegmentKey, seq: Vec<u8>, is_gap: bool) -> String {
        if let Some(id) = self.keys.get(&key) {
            return id.clone();
        }

        let id = self.next_segment_id();
        if is_gap {
            self.gap_segments += 1;
            self.total_gap_segment_bp += seq.len() as u64;
        }
        self.segments.push(GapSegment {
            id: id.clone(),
            seq,
        });
        self.keys.insert(key, id.clone());
        id
    }

    fn intern_gap(&mut self, context: GapContext, seq: Vec<u8>) -> String {
        let key = ExactSegmentKey::Gap(GapKey {
            context,
            seq: seq.clone(),
        });
        self.intern(key, seq, true)
    }

    fn intern_syncmer_slice(
        &mut self,
        path_name: &str,
        syncmer: RawSyncmerOcc,
        left_trim: u32,
        right_trim: u32,
        seq: Vec<u8>,
        private: bool,
    ) -> String {
        let key = if private {
            ExactSegmentKey::PrivateSyncmerSlice {
                path_name: path_name.to_string(),
                occurrence: syncmer.occ,
                signed: syncmer.signed,
                left_trim,
                right_trim,
            }
        } else {
            ExactSegmentKey::SyncmerSlice {
                signed: syncmer.signed,
                left_trim,
                right_trim,
            }
        };
        self.intern(key, seq, false)
    }

    fn total_extra_segments(&self) -> usize {
        self.segments.len()
    }
}

/// Walked path plus any inserted gap segments.
struct WalkedPath {
    name: String,
    steps: Vec<PathStep>,
}

fn push_walked_path_segments(
    walked_paths: &mut Vec<WalkedPath>,
    name: String,
    segments: Vec<Vec<PathStep>>,
    had_break: bool,
) {
    let non_empty: Vec<Vec<PathStep>> = segments
        .into_iter()
        .filter(|steps| !steps.is_empty())
        .collect();
    if non_empty.is_empty() {
        return;
    }
    if !had_break && non_empty.len() == 1 {
        walked_paths.push(WalkedPath {
            name,
            steps: non_empty.into_iter().next().expect("non-empty segment"),
        });
        return;
    }
    for (idx, steps) in non_empty.into_iter().enumerate() {
        walked_paths.push(WalkedPath {
            name: format!("{name}|part{}", idx + 1),
            steps,
        });
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct RawSyncmerOcc {
    signed: i32,
    occ: u32,
    pos: u64,
    left_clip: u32,
    right_clip: u32,
}

impl RawSyncmerOcc {
    fn abs(self) -> u32 {
        self.signed.unsigned_abs()
    }
}

#[derive(Clone)]
enum RawEndpoint {
    Syncmer(RawSyncmerOcc),
    Gap(usize),
}

#[derive(Clone)]
enum RawPathStep {
    Syncmer(RawSyncmerOcc),
    Gap(usize),
    Break,
}

struct RawGap {
    context: GapContext,
    seq: Vec<u8>,
}

#[derive(Default)]
struct PushedGap {
    pieces: Vec<usize>,
    starts_with_cut: bool,
    ends_with_cut: bool,
    cut_runs: usize,
    cut_bp: u64,
}

impl PushedGap {
    fn first_endpoint(&self) -> Option<RawEndpoint> {
        self.pieces.first().copied().map(raw_gap_endpoint)
    }

    fn last_endpoint(&self) -> Option<RawEndpoint> {
        self.pieces.last().copied().map(raw_gap_endpoint)
    }
}

struct PathWork {
    name: String,
    gaps: Vec<RawGap>,
    steps: Vec<RawPathStep>,
    edges: Vec<(RawEndpoint, RawEndpoint, u32)>,
    gaps_filled_with_ns: usize,
    gap_occurrences: usize,
    gap_occurrence_bp: u64,
    n_cut_runs: usize,
    n_cut_bp: u64,
    skipped: bool,
}

fn push_syncmer_step(
    steps: &mut Vec<RawPathStep>,
    next_occ: &mut u32,
    signed: i32,
    pos: u64,
    left_clip: u32,
    right_clip: u32,
) -> RawSyncmerOcc {
    let syncmer = RawSyncmerOcc {
        signed,
        occ: *next_occ,
        pos,
        left_clip,
        right_clip,
    };
    *next_occ += 1;
    steps.push(RawPathStep::Syncmer(syncmer));
    syncmer
}

fn raw_syncmer_endpoint(syncmer: RawSyncmerOcc) -> RawEndpoint {
    RawEndpoint::Syncmer(syncmer)
}

fn raw_gap_endpoint(idx: usize) -> RawEndpoint {
    RawEndpoint::Gap(idx)
}

fn is_n_base(base: u8) -> bool {
    matches!(base, b'N' | b'n')
}

fn push_raw_gap_piece(
    gaps: &mut Vec<RawGap>,
    steps: &mut Vec<RawPathStep>,
    context: GapContext,
    seq: Vec<u8>,
) -> usize {
    let gap_idx = gaps.len();
    gaps.push(RawGap { context, seq });
    steps.push(RawPathStep::Gap(gap_idx));
    gap_idx
}

fn push_gap_with_n_cutting(
    gaps: &mut Vec<RawGap>,
    steps: &mut Vec<RawPathStep>,
    context: GapContext,
    seq: Vec<u8>,
    policy: SyngGfaFrequencyMask,
) -> PushedGap {
    if seq.is_empty() {
        return PushedGap::default();
    }

    if !policy.cut_n_gaps {
        return PushedGap {
            pieces: vec![push_raw_gap_piece(gaps, steps, context, seq)],
            ..PushedGap::default()
        };
    }

    let min_run = policy.cut_n_min_run.max(1);
    let mut pushed = PushedGap::default();
    let mut chunk_start = 0usize;
    let mut idx = 0usize;
    while idx < seq.len() {
        if !is_n_base(seq[idx]) {
            idx += 1;
            continue;
        }
        let run_start = idx;
        while idx < seq.len() && is_n_base(seq[idx]) {
            idx += 1;
        }
        let run_end = idx;
        let run_len = run_end - run_start;
        if run_len < min_run {
            continue;
        }

        if run_start > chunk_start {
            let piece = push_raw_gap_piece(
                gaps,
                steps,
                context.clone(),
                seq[chunk_start..run_start].to_vec(),
            );
            pushed.pieces.push(piece);
        }
        if run_start == 0 {
            pushed.starts_with_cut = true;
        }
        if run_end == seq.len() {
            pushed.ends_with_cut = true;
        }
        pushed.cut_runs += 1;
        pushed.cut_bp += run_len as u64;
        steps.push(RawPathStep::Break);
        chunk_start = run_end;
    }

    if chunk_start < seq.len() {
        let piece = push_raw_gap_piece(gaps, steps, context, seq[chunk_start..].to_vec());
        pushed.pieces.push(piece);
    }

    pushed
}

fn raw_endpoint_to_gfa(
    endpoint: &RawEndpoint,
    gap_ids: &[String],
    cloned_syncmer_ids: &FxHashMap<u32, String>,
) -> (String, char) {
    match endpoint {
        RawEndpoint::Syncmer(syncmer) => {
            let (id, sign) = signed_to_gfa(syncmer.signed);
            if let Some(clone_id) = cloned_syncmer_ids.get(&syncmer.occ) {
                (clone_id.clone(), sign)
            } else {
                (id.to_string(), sign)
            }
        }
        RawEndpoint::Gap(idx) => (gap_ids[*idx].clone(), '+'),
    }
}

fn path_step_endpoint(step: &PathStep) -> (String, char) {
    match step {
        PathStep::Syncmer(signed) => {
            let (id, sign) = signed_to_gfa(*signed);
            (id.to_string(), sign)
        }
        PathStep::Segment { id, sign } => (id.clone(), *sign),
        PathStep::Gap(id) => (id.clone(), '+'),
    }
}

fn add_zero_edges_for_steps(
    steps: &[PathStep],
    edges: &mut FxHashSet<(String, char, String, char, u32)>,
) {
    for pair in steps.windows(2) {
        let (from_id, from_sign) = path_step_endpoint(&pair[0]);
        let (to_id, to_sign) = path_step_endpoint(&pair[1]);
        edges.insert((from_id, from_sign, to_id, to_sign, 0));
    }
}

fn syncmer_path_step_exact(
    index: &SyngIndex,
    interner: &mut ExactSegmentInterner,
    used_syncmers: &mut FxHashSet<i32>,
    path_name: &str,
    syncmer: RawSyncmerOcc,
    incoming_overlap: u32,
    private: bool,
) -> Option<PathStep> {
    let syncmer_len = index.syncmer_length_bp() as u32;
    let left_trim = syncmer.left_clip.max(incoming_overlap).min(syncmer_len);
    let right_trim = syncmer.right_clip.min(syncmer_len);
    if left_trim.saturating_add(right_trim) >= syncmer_len {
        return None;
    }
    if left_trim == 0 && right_trim == 0 && !private {
        let (id, _) = signed_to_gfa(syncmer.signed);
        used_syncmers.insert(id);
        return Some(PathStep::Syncmer(syncmer.signed));
    }

    let mut seq = index.syncmer_seq(syncmer.signed);
    seq.make_ascii_uppercase();
    let end = seq.len().saturating_sub(right_trim as usize);
    let seq = seq[left_trim as usize..end].to_vec();
    let id = interner.intern_syncmer_slice(path_name, syncmer, left_trim, right_trim, seq, private);
    Some(PathStep::Segment { id, sign: '+' })
}

fn collect_path_work(
    index: &SyngIndex,
    name: &str,
    start_info: Option<&GbwtPathStart>,
    path_len_bp: u64,
    syncmer_len_u64: u64,
    gap_fill: Option<&UnifiedSequenceIndex>,
) -> io::Result<PathWork> {
    let Some(start_info) = start_info else {
        return Ok(PathWork {
            name: name.to_string(),
            gaps: Vec::new(),
            steps: Vec::new(),
            edges: Vec::new(),
            gaps_filled_with_ns: 0,
            gap_occurrences: 0,
            gap_occurrence_bp: 0,
            n_cut_runs: 0,
            n_cut_bp: 0,
            skipped: true,
        });
    };

    let nodes = index.walk_forward_path(start_info);
    let mut gaps: Vec<RawGap> = Vec::new();
    let mut steps: Vec<RawPathStep> = Vec::with_capacity(nodes.len());
    let mut edges: Vec<(RawEndpoint, RawEndpoint, u32)> = Vec::with_capacity(nodes.len() * 2);
    let mut gaps_filled_with_ns = 0usize;
    let mut gap_occurrences = 0usize;
    let mut gap_occurrence_bp = 0u64;
    let mut next_syncmer_occ = 0u32;
    let mut last_syncmer_occ: Option<RawSyncmerOcc> = None;

    let mut push_gap = |context: GapContext, seq: Vec<u8>, steps: &mut Vec<RawPathStep>| -> usize {
        gap_occurrences += 1;
        gap_occurrence_bp += seq.len() as u64;
        let gap_idx = gaps.len();
        gaps.push(RawGap { context, seq });
        steps.push(RawPathStep::Gap(gap_idx));
        gap_idx
    };

    // Prefix: bases [0, first_syncmer_pos) before the first syncmer.
    if let Some(&(first, first_pos)) = nodes.first() {
        let prefix_gap_idx = if first_pos > 0 {
            let gap_seq = if let Some(seq_idx) = gap_fill {
                fetch_gap_dna(seq_idx, name, 0, first_pos)?
            } else {
                gaps_filled_with_ns += 1;
                vec![b'N'; first_pos as usize]
            };
            Some(push_gap(
                GapContext::Prefix { right: first },
                gap_seq,
                &mut steps,
            ))
        } else {
            None
        };
        let first_occ =
            push_syncmer_step(&mut steps, &mut next_syncmer_occ, first, first_pos, 0, 0);
        if let Some(gap_idx) = prefix_gap_idx {
            edges.push((
                raw_gap_endpoint(gap_idx),
                raw_syncmer_endpoint(first_occ),
                0,
            ));
        }
        last_syncmer_occ = Some(first_occ);
    }

    for w in nodes.windows(2) {
        let (a, pos_a) = w[0];
        let (b, pos_b) = w[1];
        let a_occ = last_syncmer_occ.expect("path window without a previous syncmer occurrence");
        let off_u64: u64 = pos_b.saturating_sub(pos_a);
        if off_u64 <= syncmer_len_u64 {
            // Direct adjacency: overlap = syncmer_len - offset.
            let overlap = (syncmer_len_u64 - off_u64) as u32;
            let b_occ = push_syncmer_step(&mut steps, &mut next_syncmer_occ, b, pos_b, 0, 0);
            edges.push((
                raw_syncmer_endpoint(a_occ),
                raw_syncmer_endpoint(b_occ),
                overlap,
            ));
            last_syncmer_occ = Some(b_occ);
        } else {
            // Inter-syncmer gap of (offset - syncmer_len) bp.
            let gap_len = off_u64 - syncmer_len_u64;
            let gap_seq = if let Some(seq_idx) = gap_fill {
                let gap_start = pos_a + syncmer_len_u64;
                let gap_end = pos_b;
                fetch_gap_dna(seq_idx, name, gap_start, gap_end)?
            } else {
                gaps_filled_with_ns += 1;
                vec![b'N'; gap_len as usize]
            };
            let gap_idx = push_gap(
                GapContext::Between { left: a, right: b },
                gap_seq,
                &mut steps,
            );
            let b_occ = push_syncmer_step(&mut steps, &mut next_syncmer_occ, b, pos_b, 0, 0);
            edges.push((raw_syncmer_endpoint(a_occ), raw_gap_endpoint(gap_idx), 0));
            edges.push((raw_gap_endpoint(gap_idx), raw_syncmer_endpoint(b_occ), 0));
            last_syncmer_occ = Some(b_occ);
        }
    }

    // Suffix: bases [last_pos + syncmer_len, path_len) after the last syncmer.
    if let Some(&(last, last_pos)) = nodes.last() {
        let suffix_start = last_pos.saturating_add(syncmer_len_u64);
        if path_len_bp > suffix_start {
            let suffix_end = path_len_bp;
            let gap_seq = if let Some(seq_idx) = gap_fill {
                fetch_gap_dna(seq_idx, name, suffix_start, suffix_end)?
            } else {
                gaps_filled_with_ns += 1;
                vec![b'N'; (suffix_end - suffix_start) as usize]
            };
            let gap_idx = push_gap(GapContext::Suffix { left: last }, gap_seq, &mut steps);
            if let Some(last_occ) = last_syncmer_occ {
                edges.push((raw_syncmer_endpoint(last_occ), raw_gap_endpoint(gap_idx), 0));
            }
        }
    }

    Ok(PathWork {
        name: name.to_string(),
        gaps,
        steps,
        edges,
        gaps_filled_with_ns,
        gap_occurrences,
        gap_occurrence_bp,
        n_cut_runs: 0,
        n_cut_bp: 0,
        skipped: false,
    })
}

fn fetch_oriented_gap_dna(
    gap_fill: Option<&UnifiedSequenceIndex>,
    seq_name: &str,
    gap_start: u64,
    gap_end: u64,
    reverse: bool,
) -> io::Result<Vec<u8>> {
    let mut seq = if let Some(seq_idx) = gap_fill {
        fetch_gap_dna(seq_idx, seq_name, gap_start, gap_end)?
    } else {
        vec![b'N'; gap_end.saturating_sub(gap_start) as usize]
    };
    if reverse {
        seq = crate::graph::reverse_complement(&seq);
    }
    Ok(seq)
}

fn oriented_syncmer_clips(
    pos: u64,
    range_start: u64,
    range_end: u64,
    syncmer_len: u64,
    reverse: bool,
) -> (u32, u32) {
    let syncmer_end = pos.saturating_add(syncmer_len);
    let low_clip = range_start.saturating_sub(pos).min(syncmer_len);
    let high_clip = syncmer_end.saturating_sub(range_end).min(syncmer_len);
    let (left, right) = if reverse {
        (high_clip, low_clip)
    } else {
        (low_clip, high_clip)
    };
    (left as u32, right as u32)
}

fn collect_range_work(
    index: &SyngIndex,
    range: &SyngGfaPathRange,
    syncmer_len_u64: u64,
    gap_fill: Option<&UnifiedSequenceIndex>,
    filtered_syncmers: Option<&FxHashSet<u32>>,
    frequency_mask: SyngGfaFrequencyMask,
) -> io::Result<PathWork> {
    let source_name = index
        .name_map
        .path_to_name
        .get(range.path_idx)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("path index {} is outside the syng name map", range.path_idx),
            )
        })?;
    let reverse = range.strand == '-';
    let mut fwd_nodes = index.walk_path_range(range.path_idx, range.start, range.end)?;
    if let Some(filtered) = filtered_syncmers {
        fwd_nodes.retain(|(node, _)| !filtered.contains(&node.unsigned_abs()));
    }
    let oriented_nodes: Vec<(i32, u64)> = if reverse {
        fwd_nodes
            .iter()
            .rev()
            .map(|(node, pos)| (-*node, *pos))
            .collect()
    } else {
        fwd_nodes
    };

    let mut gaps: Vec<RawGap> = Vec::new();
    let mut steps: Vec<RawPathStep> = Vec::with_capacity(oriented_nodes.len());
    let mut edges: Vec<(RawEndpoint, RawEndpoint, u32)> =
        Vec::with_capacity(oriented_nodes.len() * 2);
    let mut gap_occurrences = 0usize;
    let mut gap_occurrence_bp = 0u64;
    let mut gaps_filled_with_ns = 0usize;
    let mut n_cut_runs = 0usize;
    let mut n_cut_bp = 0u64;
    let mut next_syncmer_occ = 0u32;
    let mut last_syncmer_occ: Option<RawSyncmerOcc> = None;

    let mut push_gap =
        |context: GapContext, seq: Vec<u8>, steps: &mut Vec<RawPathStep>| -> PushedGap {
            gap_occurrences += 1;
            gap_occurrence_bp += seq.len() as u64;
            if gap_fill.is_none() {
                gaps_filled_with_ns += 1;
            }
            let pushed = push_gap_with_n_cutting(&mut gaps, steps, context, seq, frequency_mask);
            n_cut_runs += pushed.cut_runs;
            n_cut_bp += pushed.cut_bp;
            pushed
        };

    if oriented_nodes.is_empty() {
        if range.end > range.start {
            let seq =
                fetch_oriented_gap_dna(gap_fill, source_name, range.start, range.end, reverse)?;
            push_gap(
                GapContext::PathOnly {
                    path_idx: range.path_idx,
                    start: range.start,
                    end: range.end,
                    strand: range.strand,
                },
                seq,
                &mut steps,
            );
        }
        return Ok(PathWork {
            name: range.name.clone(),
            gaps,
            steps,
            edges,
            gaps_filled_with_ns,
            gap_occurrences,
            gap_occurrence_bp,
            n_cut_runs,
            n_cut_bp,
            skipped: false,
        });
    }

    if !reverse {
        if let Some(&(first, first_pos)) = oriented_nodes.first() {
            let prefix_gap_idx = if first_pos > range.start {
                let gap_seq =
                    fetch_oriented_gap_dna(gap_fill, source_name, range.start, first_pos, false)?;
                Some(push_gap(
                    GapContext::Prefix { right: first },
                    gap_seq,
                    &mut steps,
                ))
            } else {
                None
            };
            let (left_clip, right_clip) =
                oriented_syncmer_clips(first_pos, range.start, range.end, syncmer_len_u64, false);
            let first_occ = push_syncmer_step(
                &mut steps,
                &mut next_syncmer_occ,
                first,
                first_pos,
                left_clip,
                right_clip,
            );
            if let Some(gap) = prefix_gap_idx {
                if !gap.ends_with_cut {
                    if let Some(endpoint) = gap.last_endpoint() {
                        edges.push((endpoint, raw_syncmer_endpoint(first_occ), 0));
                    }
                }
            }
            last_syncmer_occ = Some(first_occ);
        }

        for w in oriented_nodes.windows(2) {
            let (a, pos_a) = w[0];
            let (b, pos_b) = w[1];
            let a_occ =
                last_syncmer_occ.expect("range window without a previous syncmer occurrence");
            let off_u64 = pos_b.saturating_sub(pos_a);
            if off_u64 <= syncmer_len_u64 {
                let (left_clip, right_clip) =
                    oriented_syncmer_clips(pos_b, range.start, range.end, syncmer_len_u64, false);
                let b_occ = push_syncmer_step(
                    &mut steps,
                    &mut next_syncmer_occ,
                    b,
                    pos_b,
                    left_clip,
                    right_clip,
                );
                edges.push((
                    raw_syncmer_endpoint(a_occ),
                    raw_syncmer_endpoint(b_occ),
                    (syncmer_len_u64 - off_u64) as u32,
                ));
                last_syncmer_occ = Some(b_occ);
            } else {
                let gap_seq = fetch_oriented_gap_dna(
                    gap_fill,
                    source_name,
                    pos_a + syncmer_len_u64,
                    pos_b,
                    false,
                )?;
                let gap = push_gap(
                    GapContext::Between { left: a, right: b },
                    gap_seq,
                    &mut steps,
                );
                let (left_clip, right_clip) =
                    oriented_syncmer_clips(pos_b, range.start, range.end, syncmer_len_u64, false);
                let b_occ = push_syncmer_step(
                    &mut steps,
                    &mut next_syncmer_occ,
                    b,
                    pos_b,
                    left_clip,
                    right_clip,
                );
                if !gap.starts_with_cut {
                    if let Some(endpoint) = gap.first_endpoint() {
                        edges.push((raw_syncmer_endpoint(a_occ), endpoint, 0));
                    }
                }
                if !gap.ends_with_cut {
                    if let Some(endpoint) = gap.last_endpoint() {
                        edges.push((endpoint, raw_syncmer_endpoint(b_occ), 0));
                    }
                }
                last_syncmer_occ = Some(b_occ);
            }
        }

        if let Some(&(last, last_pos)) = oriented_nodes.last() {
            let suffix_start = last_pos.saturating_add(syncmer_len_u64);
            if range.end > suffix_start {
                let gap_seq =
                    fetch_oriented_gap_dna(gap_fill, source_name, suffix_start, range.end, false)?;
                let gap = push_gap(GapContext::Suffix { left: last }, gap_seq, &mut steps);
                if let Some(last_occ) = last_syncmer_occ {
                    if !gap.starts_with_cut {
                        if let Some(endpoint) = gap.first_endpoint() {
                            edges.push((raw_syncmer_endpoint(last_occ), endpoint, 0));
                        }
                    }
                }
            }
        }
    } else {
        if let Some(&(first_oriented, last_fwd_pos)) = oriented_nodes.first() {
            let prefix_fwd_start = last_fwd_pos.saturating_add(syncmer_len_u64);
            let prefix_gap_idx = if range.end > prefix_fwd_start {
                let gap_seq = fetch_oriented_gap_dna(
                    gap_fill,
                    source_name,
                    prefix_fwd_start,
                    range.end,
                    true,
                )?;
                Some(push_gap(
                    GapContext::Prefix {
                        right: first_oriented,
                    },
                    gap_seq,
                    &mut steps,
                ))
            } else {
                None
            };
            let (left_clip, right_clip) =
                oriented_syncmer_clips(last_fwd_pos, range.start, range.end, syncmer_len_u64, true);
            let first_occ = push_syncmer_step(
                &mut steps,
                &mut next_syncmer_occ,
                first_oriented,
                last_fwd_pos,
                left_clip,
                right_clip,
            );
            if let Some(gap) = prefix_gap_idx {
                if !gap.ends_with_cut {
                    if let Some(endpoint) = gap.last_endpoint() {
                        edges.push((endpoint, raw_syncmer_endpoint(first_occ), 0));
                    }
                }
            }
            last_syncmer_occ = Some(first_occ);
        }

        for w in oriented_nodes.windows(2) {
            let (a_oriented, a_fwd_pos) = w[0];
            let (b_oriented, b_fwd_pos) = w[1];
            let a_occ =
                last_syncmer_occ.expect("range window without a previous syncmer occurrence");
            let off_u64 = a_fwd_pos.saturating_sub(b_fwd_pos);
            if off_u64 <= syncmer_len_u64 {
                let (left_clip, right_clip) = oriented_syncmer_clips(
                    b_fwd_pos,
                    range.start,
                    range.end,
                    syncmer_len_u64,
                    true,
                );
                let b_occ = push_syncmer_step(
                    &mut steps,
                    &mut next_syncmer_occ,
                    b_oriented,
                    b_fwd_pos,
                    left_clip,
                    right_clip,
                );
                edges.push((
                    raw_syncmer_endpoint(a_occ),
                    raw_syncmer_endpoint(b_occ),
                    (syncmer_len_u64 - off_u64) as u32,
                ));
                last_syncmer_occ = Some(b_occ);
            } else {
                let gap_seq = fetch_oriented_gap_dna(
                    gap_fill,
                    source_name,
                    b_fwd_pos + syncmer_len_u64,
                    a_fwd_pos,
                    true,
                )?;
                let gap = push_gap(
                    GapContext::Between {
                        left: a_oriented,
                        right: b_oriented,
                    },
                    gap_seq,
                    &mut steps,
                );
                let (left_clip, right_clip) = oriented_syncmer_clips(
                    b_fwd_pos,
                    range.start,
                    range.end,
                    syncmer_len_u64,
                    true,
                );
                let b_occ = push_syncmer_step(
                    &mut steps,
                    &mut next_syncmer_occ,
                    b_oriented,
                    b_fwd_pos,
                    left_clip,
                    right_clip,
                );
                if !gap.starts_with_cut {
                    if let Some(endpoint) = gap.first_endpoint() {
                        edges.push((raw_syncmer_endpoint(a_occ), endpoint, 0));
                    }
                }
                if !gap.ends_with_cut {
                    if let Some(endpoint) = gap.last_endpoint() {
                        edges.push((endpoint, raw_syncmer_endpoint(b_occ), 0));
                    }
                }
                last_syncmer_occ = Some(b_occ);
            }
        }

        if let Some(&(last_oriented, first_fwd_pos)) = oriented_nodes.last() {
            if first_fwd_pos > range.start {
                let gap_seq = fetch_oriented_gap_dna(
                    gap_fill,
                    source_name,
                    range.start,
                    first_fwd_pos,
                    true,
                )?;
                let gap = push_gap(
                    GapContext::Suffix {
                        left: last_oriented,
                    },
                    gap_seq,
                    &mut steps,
                );
                if let Some(last_occ) = last_syncmer_occ {
                    if !gap.starts_with_cut {
                        if let Some(endpoint) = gap.first_endpoint() {
                            edges.push((raw_syncmer_endpoint(last_occ), endpoint, 0));
                        }
                    }
                }
            }
        }
    }

    Ok(PathWork {
        name: range.name.clone(),
        gaps,
        steps,
        edges,
        gaps_filled_with_ns,
        gap_occurrences,
        gap_occurrence_bp,
        n_cut_runs,
        n_cut_bp,
        skipped: false,
    })
}

#[cfg(test)]
fn collect_range_syncmer_nodes(
    index: &SyngIndex,
    range: &SyngGfaPathRange,
) -> io::Result<FxHashSet<u32>> {
    let mut nodes = FxHashSet::default();
    for (node, _) in index.walk_path_range(range.path_idx, range.start, range.end)? {
        nodes.insert(node.unsigned_abs());
    }
    Ok(nodes)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct LocalSyncmerWalkStep {
    node: u32,
    pos: u64,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct LocalOccurrenceKey {
    path_idx: usize,
    occ: u32,
}

#[derive(Clone, Copy, Debug)]
struct LocalSyncmerOccurrence {
    key: LocalOccurrenceKey,
    node: u32,
    pos: u64,
}

#[derive(Clone, Copy, Debug)]
struct PairAnchorOccurrence {
    query_occ: u32,
    target_occ: u32,
    node: u32,
    query_pos: u64,
    target_pos: u64,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct PairAnchorKey {
    query_pos: u64,
    target_pos: u64,
    node: u32,
}

fn range_oriented_syncmer_pos(
    pos: u64,
    range: &SyngGfaPathRange,
    syncmer_len: u64,
    reverse: bool,
) -> u64 {
    if reverse {
        range.end.saturating_sub(pos.saturating_add(syncmer_len))
    } else {
        pos.saturating_sub(range.start)
    }
}

fn collect_range_syncmer_walk(
    index: &SyngIndex,
    range: &SyngGfaPathRange,
    syncmer_len: u64,
) -> io::Result<Vec<LocalSyncmerWalkStep>> {
    let mut walk = index.walk_path_range(range.path_idx, range.start, range.end)?;
    let reverse = range.strand == '-';
    if reverse {
        walk.reverse();
    }
    Ok(walk
        .into_iter()
        .map(|(node, pos)| LocalSyncmerWalkStep {
            node: node.unsigned_abs(),
            pos: range_oriented_syncmer_pos(pos, range, syncmer_len, reverse),
        })
        .collect())
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct ExactSequenceRunKey {
    nodes: Vec<u32>,
    deltas: Vec<u64>,
}

#[derive(Default)]
struct SharedContextFilterResult {
    weak_occurrences: Vec<FxHashSet<u32>>,
    shared_occurrences: usize,
    scaffold_supported_occurrences: usize,
    sequence_supported_occurrences: usize,
}

fn step_delta(a: LocalSyncmerWalkStep, b: LocalSyncmerWalkStep) -> u64 {
    a.pos.abs_diff(b.pos)
}

fn canonical_exact_sequence_run_key(window: &[LocalSyncmerWalkStep]) -> ExactSequenceRunKey {
    let nodes: Vec<u32> = window.iter().map(|step| step.node).collect();
    let deltas: Vec<u64> = window.windows(2).map(|w| step_delta(w[0], w[1])).collect();
    let mut rev_nodes = nodes.clone();
    rev_nodes.reverse();
    let mut rev_deltas = deltas.clone();
    rev_deltas.reverse();

    let forward = ExactSequenceRunKey { nodes, deltas };
    let reverse = ExactSequenceRunKey {
        nodes: rev_nodes,
        deltas: rev_deltas,
    };
    forward.min(reverse)
}

fn shared_nodes_from_local_walks(walks: &[Vec<LocalSyncmerWalkStep>]) -> FxHashSet<u32> {
    let mut counts: FxHashMap<u32, u32> = FxHashMap::default();
    for walk in walks {
        for step in walk {
            *counts.entry(step.node).or_default() += 1;
        }
    }
    counts
        .into_iter()
        .filter_map(|(node, count)| (count > 1).then_some(node))
        .collect()
}

fn shared_occurrence_keys(
    walks: &[Vec<LocalSyncmerWalkStep>],
    shared_nodes: &FxHashSet<u32>,
) -> FxHashSet<LocalOccurrenceKey> {
    let mut keys = FxHashSet::default();
    for (path_idx, walk) in walks.iter().enumerate() {
        for (occ, step) in walk.iter().enumerate() {
            if shared_nodes.contains(&step.node) {
                keys.insert(LocalOccurrenceKey {
                    path_idx,
                    occ: occ as u32,
                });
            }
        }
    }
    keys
}

fn local_occurrences_by_node(
    walks: &[Vec<LocalSyncmerWalkStep>],
    shared_nodes: &FxHashSet<u32>,
) -> FxHashMap<u32, Vec<LocalSyncmerOccurrence>> {
    let mut by_node: FxHashMap<u32, Vec<LocalSyncmerOccurrence>> = FxHashMap::default();
    for (path_idx, walk) in walks.iter().enumerate() {
        for (occ, step) in walk.iter().enumerate() {
            if !shared_nodes.contains(&step.node) {
                continue;
            }
            by_node
                .entry(step.node)
                .or_default()
                .push(LocalSyncmerOccurrence {
                    key: LocalOccurrenceKey {
                        path_idx,
                        occ: occ as u32,
                    },
                    node: step.node,
                    pos: step.pos,
                });
        }
    }
    by_node
}

fn scaffold_supported_occurrences(
    walks: &[Vec<LocalSyncmerWalkStep>],
    syncmer_len: u64,
    min_shared_run: usize,
    shared_nodes: &FxHashSet<u32>,
) -> FxHashSet<LocalOccurrenceKey> {
    if min_shared_run <= 1 || walks.len() < 2 || shared_nodes.is_empty() {
        return FxHashSet::default();
    }

    let mut pair_anchors: FxHashMap<(usize, usize), Vec<PairAnchorOccurrence>> =
        FxHashMap::default();
    for (_node, occurrences) in local_occurrences_by_node(walks, shared_nodes) {
        for left_idx in 0..occurrences.len() {
            for right_idx in (left_idx + 1)..occurrences.len() {
                let left = occurrences[left_idx];
                let right = occurrences[right_idx];
                if left.key.path_idx == right.key.path_idx {
                    continue;
                }
                let (query, target) = if left.key.path_idx < right.key.path_idx {
                    (left, right)
                } else {
                    (right, left)
                };
                pair_anchors
                    .entry((query.key.path_idx, target.key.path_idx))
                    .or_default()
                    .push(PairAnchorOccurrence {
                        query_occ: query.key.occ,
                        target_occ: target.key.occ,
                        node: query.node,
                        query_pos: query.pos,
                        target_pos: target.pos,
                    });
            }
        }
    }

    let min_scaffold_length = (min_shared_run as u64).saturating_mul(syncmer_len.max(1));
    let mut supported = FxHashSet::default();
    for ((query_idx, target_idx), anchors) in pair_anchors {
        supported.extend(retained_pair_scaffold_occurrences(
            query_idx,
            target_idx,
            &anchors,
            syncmer_len,
            min_scaffold_length,
        ));
    }
    supported
}

fn retained_pair_scaffold_occurrences(
    query_idx: usize,
    target_idx: usize,
    anchors: &[PairAnchorOccurrence],
    syncmer_len: u64,
    min_scaffold_length: u64,
) -> FxHashSet<LocalOccurrenceKey> {
    if anchors.is_empty() {
        return FxHashSet::default();
    }

    let mut occurrence_by_anchor: FxHashMap<PairAnchorKey, Vec<(u32, u32)>> = FxHashMap::default();
    for anchor in anchors {
        occurrence_by_anchor
            .entry(PairAnchorKey {
                query_pos: anchor.query_pos,
                target_pos: anchor.target_pos,
                node: anchor.node,
            })
            .or_default()
            .push((anchor.query_occ, anchor.target_occ));
    }

    let mut chain_anchors: Vec<Anchor> = occurrence_by_anchor
        .keys()
        .map(|key| Anchor {
            query_pos: key.query_pos,
            target_pos: key.target_pos,
            node_id: key.node,
        })
        .collect();
    chain_anchors.sort_by(|a, b| {
        a.query_pos
            .cmp(&b.query_pos)
            .then(a.target_pos.cmp(&b.target_pos))
            .then(a.node_id.cmp(&b.node_id))
    });

    let start = chain_anchors
        .iter()
        .map(|anchor| anchor.target_pos)
        .min()
        .unwrap_or(0);
    let end = chain_anchors
        .iter()
        .map(|anchor| anchor.target_pos.saturating_add(syncmer_len))
        .max()
        .unwrap_or(start);
    let hits = vec![HomologousIntervalWithAnchors {
        genome: format!("__syng2gfa_target_{target_idx}"),
        start,
        end,
        strand: '+',
        anchors: chain_anchors,
    }];
    let retained = crate::syng_transitive::chain_anchors_with_sweepga_scaffold_mass(
        hits,
        syncmer_len,
        crate::syng_transitive::DEFAULT_EXTEND_BUDGET_BP,
        min_scaffold_length,
    );

    let mut supported = FxHashSet::default();
    for chain in retained {
        for anchor in chain.anchors {
            let key = PairAnchorKey {
                query_pos: anchor.query_pos,
                target_pos: anchor.target_pos,
                node: anchor.node_id,
            };
            if let Some(occurrences) = occurrence_by_anchor.get(&key) {
                for &(query_occ, target_occ) in occurrences {
                    supported.insert(LocalOccurrenceKey {
                        path_idx: query_idx,
                        occ: query_occ,
                    });
                    supported.insert(LocalOccurrenceKey {
                        path_idx: target_idx,
                        occ: target_occ,
                    });
                }
            }
        }
    }
    supported
}

fn sequence_supported_occurrences(
    walks: &[Vec<LocalSyncmerWalkStep>],
    syncmer_len: u64,
    min_sequence_span_bp: usize,
    shared_nodes: &FxHashSet<u32>,
) -> FxHashSet<LocalOccurrenceKey> {
    if min_sequence_span_bp == 0 || walks.is_empty() || shared_nodes.is_empty() {
        return FxHashSet::default();
    }
    if min_sequence_span_bp as u64 <= syncmer_len {
        return shared_occurrence_keys(walks, shared_nodes);
    }

    let min_span = min_sequence_span_bp as u64;
    let mut run_counts: FxHashMap<ExactSequenceRunKey, u32> = FxHashMap::default();
    let mut candidate_runs: Vec<(ExactSequenceRunKey, Vec<LocalOccurrenceKey>)> = Vec::new();
    for (path_idx, walk) in walks.iter().enumerate() {
        for start in 0..walk.len() {
            let mut end = start;
            while end + 1 < walk.len() {
                let delta = step_delta(walk[end], walk[end + 1]);
                if delta > syncmer_len {
                    break;
                }
                end += 1;
                let span = step_delta(walk[start], walk[end]) + syncmer_len;
                if span >= min_span {
                    let window = &walk[start..=end];
                    let key = canonical_exact_sequence_run_key(window);
                    let occurrences = window
                        .iter()
                        .enumerate()
                        .filter_map(|(offset, step)| {
                            shared_nodes
                                .contains(&step.node)
                                .then_some(LocalOccurrenceKey {
                                    path_idx,
                                    occ: (start + offset) as u32,
                                })
                        })
                        .collect();
                    *run_counts.entry(key.clone()).or_default() += 1;
                    candidate_runs.push((key, occurrences));
                    break;
                }
            }
        }
    }

    let mut supported = FxHashSet::default();
    for (key, occurrences) in candidate_runs {
        if run_counts.get(&key).copied().unwrap_or(0) > 1 {
            supported.extend(occurrences);
        }
    }
    supported
}

fn shared_context_filtered_occurrences(
    walks: &[Vec<LocalSyncmerWalkStep>],
    syncmer_len: u64,
    min_shared_run: usize,
    min_sequence_span_bp: usize,
) -> SharedContextFilterResult {
    if walks.is_empty() || (min_shared_run <= 1 && min_sequence_span_bp == 0) {
        return SharedContextFilterResult {
            weak_occurrences: vec![FxHashSet::default(); walks.len()],
            ..SharedContextFilterResult::default()
        };
    }

    let shared_nodes = shared_nodes_from_local_walks(walks);
    if shared_nodes.is_empty() {
        return SharedContextFilterResult {
            weak_occurrences: vec![FxHashSet::default(); walks.len()],
            ..SharedContextFilterResult::default()
        };
    }

    let shared_occurrences = shared_occurrence_keys(walks, &shared_nodes);
    let scaffold_supported =
        scaffold_supported_occurrences(walks, syncmer_len, min_shared_run, &shared_nodes);
    let sequence_supported =
        sequence_supported_occurrences(walks, syncmer_len, min_sequence_span_bp, &shared_nodes);
    let supported: FxHashSet<LocalOccurrenceKey> = scaffold_supported
        .union(&sequence_supported)
        .copied()
        .collect();
    let mut weak_occurrences = vec![FxHashSet::default(); walks.len()];
    for key in shared_occurrences.difference(&supported).copied() {
        if let Some(path_set) = weak_occurrences.get_mut(key.path_idx) {
            path_set.insert(key.occ);
        }
    }

    SharedContextFilterResult {
        weak_occurrences,
        shared_occurrences: shared_occurrences.len(),
        scaffold_supported_occurrences: scaffold_supported.len(),
        sequence_supported_occurrences: sequence_supported.len(),
    }
}

fn local_walks_after_frequency_mask(
    walks: &[Vec<LocalSyncmerWalkStep>],
    masked_syncmers: &FxHashSet<u32>,
) -> Vec<Vec<LocalSyncmerWalkStep>> {
    if masked_syncmers.is_empty() {
        return walks.to_vec();
    }
    walks
        .iter()
        .map(|walk| {
            walk.iter()
                .copied()
                .filter(|step| !masked_syncmers.contains(&step.node))
                .collect()
        })
        .collect()
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
enum RawContextSide {
    Boundary,
    Gap,
    Syncmer(i32),
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct RawSyncmerContext {
    left: RawContextSide,
    right: RawContextSide,
}

#[derive(Clone, Copy, Debug)]
struct RawSyncmerContextOcc {
    syncmer: RawSyncmerOcc,
    context: RawSyncmerContext,
}

fn raw_step_context_side(step: Option<&RawPathStep>) -> RawContextSide {
    match step {
        Some(RawPathStep::Syncmer(syncmer)) => RawContextSide::Syncmer(syncmer.signed),
        Some(RawPathStep::Gap(_)) => RawContextSide::Gap,
        Some(RawPathStep::Break) | None => RawContextSide::Boundary,
    }
}

fn raw_work_syncmer_context_occurrences(work: &PathWork) -> Vec<RawSyncmerContextOcc> {
    let mut occurrences = Vec::new();
    for (idx, step) in work.steps.iter().enumerate() {
        if let RawPathStep::Syncmer(syncmer) = step {
            occurrences.push(RawSyncmerContextOcc {
                syncmer: *syncmer,
                context: RawSyncmerContext {
                    left: raw_step_context_side(idx.checked_sub(1).and_then(|i| work.steps.get(i))),
                    right: raw_step_context_side(work.steps.get(idx + 1)),
                },
            });
        }
    }
    occurrences
}

fn dominant_local_repeat_contexts(
    path_work: &[PathWork],
    max_minor: u32,
    min_dominant_fraction: f64,
) -> FxHashMap<u32, RawSyncmerContext> {
    if max_minor == 0 || path_work.is_empty() {
        return FxHashMap::default();
    }

    let mut counts: FxHashMap<u32, FxHashMap<RawSyncmerContext, u32>> = FxHashMap::default();
    for work in path_work {
        for occ in raw_work_syncmer_context_occurrences(work) {
            *counts
                .entry(occ.syncmer.abs())
                .or_default()
                .entry(occ.context)
                .or_default() += 1;
        }
    }

    let mut dominant = FxHashMap::default();
    for (node, context_counts) in counts {
        if context_counts.len() <= 1 {
            continue;
        }
        let total: u32 = context_counts.values().sum();
        let mut ranked: Vec<(RawSyncmerContext, u32)> = context_counts.into_iter().collect();
        ranked.sort_by(|(_, a), (_, b)| b.cmp(a));
        let (context, max_count) = ranked[0];
        if ranked.get(1).is_some_and(|(_, count)| *count == max_count) {
            continue;
        }
        let minor = total.saturating_sub(max_count);
        let dominant_fraction = max_count as f64 / total as f64;
        if minor > 0 && minor <= max_minor && dominant_fraction >= min_dominant_fraction {
            dominant.insert(node, context);
        }
    }
    dominant
}

fn local_repeat_clone_occurrences(
    work: &PathWork,
    dominant_contexts: &FxHashMap<u32, RawSyncmerContext>,
) -> FxHashSet<u32> {
    if dominant_contexts.is_empty() {
        return FxHashSet::default();
    }

    let occurrences = raw_work_syncmer_context_occurrences(work);
    let mut path_counts: FxHashMap<u32, u32> = FxHashMap::default();
    for occ in &occurrences {
        *path_counts.entry(occ.syncmer.abs()).or_default() += 1;
    }

    let mut clone_occurrences = FxHashSet::default();
    for occ in occurrences {
        let node = occ.syncmer.abs();
        if path_counts.get(&node).copied().unwrap_or(0) <= 1 {
            continue;
        }
        if dominant_contexts
            .get(&node)
            .is_some_and(|dominant| *dominant != occ.context)
        {
            clone_occurrences.insert(occ.syncmer.occ);
        }
    }
    clone_occurrences
}

fn split_shared_context_clone_occurrences(
    work: &PathWork,
    split_occurrences: &FxHashSet<u32>,
) -> FxHashSet<u32> {
    if split_occurrences.is_empty() {
        return FxHashSet::default();
    }
    raw_work_syncmer_context_occurrences(work)
        .into_iter()
        .filter_map(|occ| {
            split_occurrences
                .contains(&occ.syncmer.occ)
                .then_some(occ.syncmer.occ)
        })
        .collect()
}

fn top_frequency_nodes(occurrence_counts: &[(u32, u32)], drop_top_fraction: f64) -> FxHashSet<u32> {
    if drop_top_fraction <= 0.0 || occurrence_counts.is_empty() {
        return FxHashSet::default();
    }
    let drop_count = ((occurrence_counts.len() as f64) * drop_top_fraction).floor() as usize;
    if drop_count == 0 {
        return FxHashSet::default();
    }
    let mut ranked = occurrence_counts.to_vec();
    ranked.sort_by(|(node_a, count_a), (node_b, count_b)| {
        count_b.cmp(count_a).then(node_a.cmp(node_b))
    });
    ranked
        .into_iter()
        .take(drop_count.min(occurrence_counts.len()))
        .map(|(node, _)| node)
        .collect()
}

fn frequency_masked_syncmers(
    index: &SyngIndex,
    used_nodes: &FxHashSet<u32>,
    mask: SyngGfaFrequencyMask,
) -> io::Result<FxHashSet<u32>> {
    if (mask.drop_top_fraction <= 0.0 && mask.max_occurrences.is_none()) || used_nodes.is_empty() {
        return Ok(FxHashSet::default());
    }

    let count_start = Instant::now();
    let mut occurrence_counts = Vec::with_capacity(used_nodes.len());
    let mut max_occurrences = 0u32;
    let mut nodes: Vec<u32> = used_nodes.iter().copied().collect();
    nodes.sort_unstable();
    for node in nodes {
        let count = index.syncmer_node_occurrence_count(node)?;
        max_occurrences = max_occurrences.max(count);
        occurrence_counts.push((node, count));
    }

    let mut masked = top_frequency_nodes(&occurrence_counts, mask.drop_top_fraction);
    if let Some(max_allowed) = mask.max_occurrences {
        for (node, count) in &occurrence_counts {
            if *count > max_allowed {
                masked.insert(*node);
            }
        }
    }
    if !masked.is_empty() {
        let masked_max = occurrence_counts
            .iter()
            .filter(|(node, _)| masked.contains(node))
            .map(|(_, count)| *count)
            .max()
            .unwrap_or(0);
        let masked_min = occurrence_counts
            .iter()
            .filter(|(node, _)| masked.contains(node))
            .map(|(_, count)| *count)
            .min()
            .unwrap_or(0);
        info!(
            "[syng2gfa] frequency-filtered {} / {} local syncmer node(s) before raw GFA materialization (top_fraction={}, max_occurrences={:?}, local max={}, filtered min/max={}/{}) in {:.3}s",
            masked.len(),
            occurrence_counts.len(),
            mask.drop_top_fraction,
            mask.max_occurrences,
            max_occurrences,
            masked_min,
            masked_max,
            count_start.elapsed().as_secs_f64()
        );
    }
    Ok(masked)
}

fn push_exact_walked_path_segments(
    walked_paths: &mut Vec<WalkedPath>,
    edges: &mut FxHashSet<(String, char, String, char, u32)>,
    name: String,
    segments: Vec<Vec<PathStep>>,
    had_break: bool,
) {
    let non_empty: Vec<Vec<PathStep>> = segments
        .into_iter()
        .filter(|steps| !steps.is_empty())
        .collect();
    if non_empty.is_empty() {
        return;
    }
    if !had_break && non_empty.len() == 1 {
        let steps = non_empty.into_iter().next().expect("non-empty segment");
        add_zero_edges_for_steps(&steps, edges);
        walked_paths.push(WalkedPath { name, steps });
        return;
    }
    for (idx, steps) in non_empty.into_iter().enumerate() {
        add_zero_edges_for_steps(&steps, edges);
        walked_paths.push(WalkedPath {
            name: format!("{name}|part{}", idx + 1),
            steps,
        });
    }
}

fn write_exact_blunt_path_work_gfa<W: Write>(
    index: &SyngIndex,
    path_work: Vec<PathWork>,
    writer: &mut W,
    version: GfaVersion,
    split_occurrences: &[FxHashSet<u32>],
    frequency_mask: SyngGfaFrequencyMask,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    unsafe { syng_ffi::impg_syng_suppress_debug() };
    writeln!(writer, "H\tVN:Z:{}", version.header_tag())?;

    let reduce_start = Instant::now();
    let mut walked_paths: Vec<WalkedPath> = Vec::with_capacity(path_work.len());
    let mut edges: FxHashSet<(String, char, String, char, u32)> = FxHashSet::default();
    let mut interner = ExactSegmentInterner::new(index.num_syncmer_nodes());
    let mut used_syncmers: FxHashSet<i32> = FxHashSet::default();
    let dominant_repeat_contexts = dominant_local_repeat_contexts(
        &path_work,
        frequency_mask.local_repeat_max_minor,
        frequency_mask.local_repeat_min_dominant_fraction,
    );
    if !dominant_repeat_contexts.is_empty() {
        info!(
            "[syng2gfa] detected {} near-single-copy syncmer node(s) with rare repeated local contexts (max_minor={}, min_dominant_fraction={})",
            dominant_repeat_contexts.len(),
            frequency_mask.local_repeat_max_minor,
            frequency_mask.local_repeat_min_dominant_fraction
        );
    }

    let mut n_skipped = 0usize;
    let mut gaps_filled_with_ns = 0usize;
    let mut gap_occurrences = 0usize;
    let mut gap_occurrence_bp = 0u64;
    let mut n_cut_runs = 0usize;
    let mut n_cut_bp = 0u64;
    let mut n_local_repeat_clones = 0usize;
    let mut n_scaffold_context_clones = 0usize;

    for (work_idx, work) in path_work.into_iter().enumerate() {
        if work.skipped {
            warn!(
                "[syng2gfa] path '{}' has no GBWT start info — skipping",
                work.name
            );
            n_skipped += 1;
            continue;
        }

        gaps_filled_with_ns += work.gaps_filled_with_ns;
        gap_occurrences += work.gap_occurrences;
        gap_occurrence_bp += work.gap_occurrence_bp;
        n_cut_runs += work.n_cut_runs;
        n_cut_bp += work.n_cut_bp;

        let split_clone_occurrences = split_occurrences
            .get(work_idx)
            .map(|occurrences| split_shared_context_clone_occurrences(&work, occurrences))
            .unwrap_or_default();
        let local_repeat_clone_occurrences =
            local_repeat_clone_occurrences(&work, &dominant_repeat_contexts);
        let mut clone_occurrences = split_clone_occurrences.clone();
        clone_occurrences.extend(local_repeat_clone_occurrences.iter().copied());
        n_scaffold_context_clones += split_clone_occurrences.len();
        n_local_repeat_clones += local_repeat_clone_occurrences.len();

        let mut incoming_overlap_by_occ: FxHashMap<u32, u32> = FxHashMap::default();
        for (_, to, overlap) in &work.edges {
            if let RawEndpoint::Syncmer(syncmer) = to {
                let entry = incoming_overlap_by_occ.entry(syncmer.occ).or_default();
                *entry = (*entry).max(*overlap);
            }
        }

        let gap_ids: Vec<String> = work
            .gaps
            .into_iter()
            .map(|gap| interner.intern_gap(gap.context, gap.seq))
            .collect();

        let mut segments = vec![Vec::new()];
        let mut had_break = false;
        for step in work.steps {
            match step {
                RawPathStep::Syncmer(syncmer) => {
                    let incoming_overlap = incoming_overlap_by_occ
                        .get(&syncmer.occ)
                        .copied()
                        .unwrap_or(0);
                    let private = clone_occurrences.contains(&syncmer.occ);
                    if let Some(path_step) = syncmer_path_step_exact(
                        index,
                        &mut interner,
                        &mut used_syncmers,
                        &work.name,
                        syncmer,
                        incoming_overlap,
                        private,
                    ) {
                        segments.last_mut().expect("path segment").push(path_step);
                    }
                }
                RawPathStep::Gap(idx) => {
                    segments
                        .last_mut()
                        .expect("path segment")
                        .push(PathStep::Gap(gap_ids[idx].clone()));
                }
                RawPathStep::Break => {
                    had_break = true;
                    segments.push(Vec::new());
                }
            }
        }

        push_exact_walked_path_segments(
            &mut walked_paths,
            &mut edges,
            work.name,
            segments,
            had_break,
        );
    }

    info!(
        "[syng2gfa] reduced exact blunt path walks into {} full untrimmed syncmer anchor segment(s), {} interned exact sequence segment(s) ({} gap, {} trimmed syncmer-slice; {} local-repeat clones, {} scaffold-context clones), {} edge(s) in {:.3}s",
        used_syncmers.len(),
        interner.total_extra_segments(),
        interner.gap_segments,
        interner.total_extra_segments().saturating_sub(interner.gap_segments),
        n_local_repeat_clones,
        n_scaffold_context_clones,
        edges.len(),
        reduce_start.elapsed().as_secs_f64()
    );

    if gaps_filled_with_ns > 0 {
        warn!(
            "[syng2gfa] filled {} gap occurrence(s) with 'N' ({} bp over paths; {} unique segment bp). \
                 Concatenating these paths will NOT reconstruct the original genome. \
                 Pass --sequence-files <FASTA> to splice in real DNA.",
            gaps_filled_with_ns,
            gap_occurrence_bp,
            interner.total_gap_segment_bp
        );
    } else if interner.gap_segments > 0 {
        info!(
            "[syng2gfa] spliced {} unique real-DNA gap segment(s) from {} gap occurrence(s), {} unique bp ({} bp over paths)",
            interner.gap_segments,
            gap_occurrences,
            interner.total_gap_segment_bp,
            gap_occurrence_bp
        );
    }
    if n_cut_runs > 0 {
        info!(
            "[syng2gfa] cut {} N-run(s) ({} bp) out of gap DNA and split emitted paths at those breaks",
            n_cut_runs,
            n_cut_bp
        );
    }

    let write_start = Instant::now();
    let mut used_nodes: Vec<i32> = used_syncmers.into_iter().collect();
    used_nodes.sort_unstable();
    for node_id in &used_nodes {
        let mut seq = index.syncmer_seq(*node_id);
        seq.make_ascii_uppercase();
        write_segment(writer, &node_id.to_string(), &seq)?;
    }
    for segment in &interner.segments {
        write_segment(writer, &segment.id, &segment.seq)?;
    }

    let mut edge_vec: Vec<(String, char, String, char, u32)> = edges.into_iter().collect();
    edge_vec.sort_unstable();
    for (a_id, a_sign, b_id, b_sign, overlap) in &edge_vec {
        writeln!(writer, "L\t{a_id}\t{a_sign}\t{b_id}\t{b_sign}\t{overlap}M")?;
    }

    for walked in &walked_paths {
        match version {
            GfaVersion::V1_0 => write_p_line(writer, &walked.name, &walked.steps)?,
            GfaVersion::V1_1 => {
                if let Some(parts) = parse_pansn(&walked.name) {
                    write_w_line(writer, &parts, &walked.steps)?;
                } else {
                    warn!(
                        "[syng2gfa] path '{}' is not PanSN-style; falling back to P line under GFA 1.1",
                        walked.name
                    );
                    write_p_line(writer, &walked.name, &walked.steps)?;
                }
            }
        }
    }

    let n_segments = used_nodes.len() + interner.total_extra_segments();
    let n_links = edge_vec.len();
    let n_paths_emitted = walked_paths.len();
    info!(
        "[syng2gfa] wrote exact blunt GFA: {} S, {} L, {} path line(s) in {:.3}s",
        n_segments,
        n_links,
        n_paths_emitted,
        write_start.elapsed().as_secs_f64()
    );

    Ok((
        n_segments,
        n_links,
        n_paths_emitted,
        n_skipped,
        interner.gap_segments,
        interner.total_gap_segment_bp,
    ))
}

/// Walk every forward path in `index`, gap-fill via FASTA/AGC if given,
/// and write a GFA to `writer`.
///
/// `gap_fill` selects how mid-path gaps are filled (FASTA vs. NNN).
/// Returns `(n_segments, n_links, n_paths_emitted, n_paths_skipped,
/// n_gap_segments, unique_gap_segment_bp)`.
pub fn write_gfa<W: Write>(
    index: &SyngIndex,
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    let n_nodes = index.num_syncmer_nodes();
    let syncmer_len = index.syncmer_length_bp();
    info!(
        "[syng2gfa] {} syncmer segments (length {} bp), GFA {}, gap fill: {}",
        n_nodes,
        syncmer_len,
        version.header_tag(),
        if gap_fill.is_some() {
            "FASTA/AGC"
        } else {
            "Ns"
        }
    );

    // Silence syng's C-side debug printfs (syngBWTpathStartOld et al.)
    unsafe { syng_ffi::impg_syng_suppress_debug() };

    writeln!(writer, "H\tVN:Z:{}", version.header_tag())?;

    // Walk every forward path and collect:
    //   - per-path step list (syncmers + gap segments) for P/W emission
    //   - unique L edges with overlaps
    //   - in-order list of gap segments to emit later as S lines
    let mut walked_paths: Vec<WalkedPath> = Vec::with_capacity(index.name_map.path_to_name.len());
    // Edge tuple: (from_id, from_sign, to_id, to_sign, overlap_bp).
    // `String` because endpoints are integer IDs (syncmers 1..N, gaps
    // N+1, N+2, …) and we serialize them as text.
    let mut edges: FxHashSet<(String, char, String, char, u32)> = FxHashSet::default();
    let mut gap_interner = GapInterner::new(n_nodes);
    let mut n_skipped = 0usize;
    let mut gaps_filled_with_ns: usize = 0;
    let mut gap_occurrences: usize = 0;
    let mut gap_occurrence_bp: u64 = 0;

    let syncmer_len_u64 = syncmer_len as u64;
    let collect_start = Instant::now();
    let path_work: Vec<PathWork> = (0..index.name_map.path_to_name.len())
        .into_par_iter()
        .map(|path_idx| {
            let name = &index.name_map.path_to_name[path_idx];
            let start_info = index
                .name_map
                .path_starts
                .get(path_idx)
                .and_then(|o| o.as_ref());
            let path_len_bp = index
                .name_map
                .path_to_length
                .get(path_idx)
                .copied()
                .unwrap_or(0);
            collect_path_work(
                index,
                name,
                start_info,
                path_len_bp,
                syncmer_len_u64,
                gap_fill,
            )
        })
        .collect::<io::Result<Vec<_>>>()?;
    info!(
        "[syng2gfa] collected {} path walks in {:.3}s using {} rayon threads",
        path_work.len(),
        collect_start.elapsed().as_secs_f64(),
        rayon::current_num_threads()
    );

    let reduce_start = Instant::now();
    for work in path_work {
        if work.skipped {
            warn!(
                "[syng2gfa] path '{}' has no GBWT start info — skipping (rebuild the syng index to include path starts)",
                work.name
            );
            n_skipped += 1;
            continue;
        }

        gaps_filled_with_ns += work.gaps_filled_with_ns;
        gap_occurrences += work.gap_occurrences;
        gap_occurrence_bp += work.gap_occurrence_bp;

        let gap_ids: Vec<String> = work
            .gaps
            .into_iter()
            .map(|gap| gap_interner.intern(gap.context, gap.seq))
            .collect();
        let cloned_syncmer_ids = FxHashMap::default();

        for (from, to, overlap) in work.edges {
            let (from_id, from_sign) = raw_endpoint_to_gfa(&from, &gap_ids, &cloned_syncmer_ids);
            let (to_id, to_sign) = raw_endpoint_to_gfa(&to, &gap_ids, &cloned_syncmer_ids);
            edges.insert((from_id, from_sign, to_id, to_sign, overlap));
        }

        let mut segments = vec![Vec::new()];
        let mut had_break = false;
        for step in work.steps {
            match step {
                RawPathStep::Syncmer(syncmer) => {
                    segments
                        .last_mut()
                        .expect("path segment")
                        .push(PathStep::Syncmer(syncmer.signed));
                }
                RawPathStep::Gap(idx) => {
                    segments
                        .last_mut()
                        .expect("path segment")
                        .push(PathStep::Gap(gap_ids[idx].clone()));
                }
                RawPathStep::Break => {
                    had_break = true;
                    segments.push(Vec::new());
                }
            }
        }
        push_walked_path_segments(&mut walked_paths, work.name, segments, had_break);
    }
    info!(
        "[syng2gfa] reduced path walks into {} unique gap segment(s) and {} edge(s) in {:.3}s",
        gap_interner.len(),
        edges.len(),
        reduce_start.elapsed().as_secs_f64()
    );

    if gaps_filled_with_ns > 0 {
        warn!(
            "[syng2gfa] filled {} gap occurrence(s) with 'N' ({} bp over paths; {} unique segment bp). \
                 Concatenating these paths will NOT reconstruct the original genome. \
                 Pass --sequence-files <FASTA> to splice in real DNA.",
            gaps_filled_with_ns,
            gap_occurrence_bp,
            gap_interner.total_segment_bp()
        );
    } else if !gap_interner.is_empty() {
        info!(
            "[syng2gfa] spliced {} unique real-DNA gap segment(s) from {} gap occurrence(s), {} unique bp ({} bp over paths)",
            gap_interner.len(),
            gap_occurrences,
            gap_interner.total_segment_bp(),
            gap_occurrence_bp
        );
    }

    // S lines: syncmer nodes 1..=N first, then gap segments in walk order.
    // Syncmer DNA from the C lib is lowercase; uppercase for consistency with
    // FASTA-fetched gap sequences (already uppercased upstream).
    let write_start = Instant::now();
    for node_id in 1..=n_nodes as i32 {
        let mut seq = index.syncmer_seq(node_id);
        seq.make_ascii_uppercase();
        write_segment(writer, &node_id.to_string(), &seq)?;
    }
    for gap in &gap_interner.segments {
        write_segment(writer, &gap.id, &gap.seq)?;
    }
    let n_segments = n_nodes + gap_interner.len();
    let n_gap_segments = gap_interner.len();

    debug!("[syng2gfa] collected {} unique edges", edges.len());

    // L lines, sorted for deterministic output.
    let mut edge_vec: Vec<(String, char, String, char, u32)> = edges.into_iter().collect();
    edge_vec.sort_unstable();
    for (a_id, a_sign, b_id, b_sign, overlap) in &edge_vec {
        writeln!(writer, "L\t{a_id}\t{a_sign}\t{b_id}\t{b_sign}\t{overlap}M")?;
    }
    let n_links = edge_vec.len();

    // P/W lines.
    let n_paths_emitted = walked_paths.len();
    for walked in &walked_paths {
        match version {
            GfaVersion::V1_0 => write_p_line(writer, &walked.name, &walked.steps)?,
            GfaVersion::V1_1 => {
                if let Some(parts) = parse_pansn(&walked.name) {
                    write_w_line(writer, &parts, &walked.steps)?;
                } else {
                    warn!(
                        "[syng2gfa] path '{}' is not PanSN-style; falling back to P line under GFA 1.1",
                        walked.name
                    );
                    write_p_line(writer, &walked.name, &walked.steps)?;
                }
            }
        }
    }

    info!(
        "[syng2gfa] wrote {} S, {} L, {} path line(s) in {:.3}s",
        n_segments,
        n_links,
        n_paths_emitted,
        write_start.elapsed().as_secs_f64()
    );

    Ok((
        n_segments,
        n_links,
        n_paths_emitted,
        n_skipped,
        n_gap_segments,
        gap_interner.total_segment_bp(),
    ))
}

/// Write every forward path as an exact, zero-overlap GFA.
fn write_gfa_blunt_exact<W: Write>(
    index: &SyngIndex,
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    if version != GfaVersion::V1_0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--gfa-mode blunt currently requires --gfa-version 1.0",
        ));
    }

    let syncmer_len = index.syncmer_length_bp();
    info!(
        "[syng2gfa] rendering exact blunt GFA for {} syncmer segments (length {} bp), gap fill: {}",
        index.num_syncmer_nodes(),
        syncmer_len,
        if gap_fill.is_some() {
            "FASTA/AGC"
        } else {
            "Ns"
        }
    );

    let syncmer_len_u64 = syncmer_len as u64;
    let collect_start = Instant::now();
    let path_work: Vec<PathWork> = (0..index.name_map.path_to_name.len())
        .into_par_iter()
        .map(|path_idx| {
            let name = &index.name_map.path_to_name[path_idx];
            let start_info = index
                .name_map
                .path_starts
                .get(path_idx)
                .and_then(|o| o.as_ref());
            let path_len_bp = index
                .name_map
                .path_to_length
                .get(path_idx)
                .copied()
                .unwrap_or(0);
            collect_path_work(
                index,
                name,
                start_info,
                path_len_bp,
                syncmer_len_u64,
                gap_fill,
            )
        })
        .collect::<io::Result<Vec<_>>>()?;
    info!(
        "[syng2gfa] collected {} path walks for exact blunt materialization in {:.3}s using {} rayon threads",
        path_work.len(),
        collect_start.elapsed().as_secs_f64(),
        rayon::current_num_threads()
    );

    write_exact_blunt_path_work_gfa(
        index,
        path_work,
        writer,
        version,
        &[],
        SyngGfaFrequencyMask::disabled(),
    )
}

fn write_range_gfa_blunt_exact<W: Write>(
    index: &SyngIndex,
    ranges: &[SyngGfaPathRange],
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
    frequency_mask: SyngGfaFrequencyMask,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    if version != GfaVersion::V1_0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--gfa-mode blunt currently requires --gfa-version 1.0",
        ));
    }

    let syncmer_len = index.syncmer_length_bp();
    info!(
        "[syng2gfa] rendering {} selected path range(s) as exact blunt GFA, gap fill: {}",
        ranges.len(),
        if gap_fill.is_some() {
            "FASTA/AGC"
        } else {
            "Ns"
        }
    );

    let syncmer_len_u64 = syncmer_len as u64;
    let need_local_walks =
        frequency_mask.frequency_filter_enabled() || frequency_mask.shared_context_filter_enabled();
    let local_walks = if need_local_walks {
        let mask_collect_start = Instant::now();
        let walks: Vec<Vec<LocalSyncmerWalkStep>> = ranges
            .par_iter()
            .map(|range| collect_range_syncmer_walk(index, range, syncmer_len_u64))
            .collect::<io::Result<Vec<_>>>()?;
        info!(
            "[syng2gfa] collected {} local syncmer walk(s) for exact blunt masking in {:.3}s using {} rayon threads",
            walks.len(),
            mask_collect_start.elapsed().as_secs_f64(),
            rayon::current_num_threads()
        );
        Some(walks)
    } else {
        None
    };
    let all_local_syncmers = local_walks.as_ref().map(|walks| {
        walks
            .iter()
            .flatten()
            .map(|step| step.node)
            .collect::<FxHashSet<_>>()
    });
    let frequency_masked_syncmers = if frequency_mask.frequency_filter_enabled() {
        frequency_masked_syncmers(
            index,
            all_local_syncmers.as_ref().expect("local syncmers"),
            frequency_mask,
        )?
    } else {
        FxHashSet::default()
    };
    let masked_syncmers = frequency_masked_syncmers.clone();
    let mut split_occurrences = vec![FxHashSet::default(); ranges.len()];
    if frequency_mask.shared_context_filter_enabled() {
        let context_start = Instant::now();
        let context_walks = local_walks_after_frequency_mask(
            local_walks.as_ref().expect("local syncmer walks"),
            &frequency_masked_syncmers,
        );
        let context_filter = shared_context_filtered_occurrences(
            &context_walks,
            syncmer_len_u64,
            frequency_mask.min_shared_run,
            frequency_mask.min_sequence_span_bp,
        );
        let weak_occurrences = context_filter
            .weak_occurrences
            .iter()
            .map(FxHashSet::len)
            .sum::<usize>();
        split_occurrences = context_filter.weak_occurrences;
        info!(
            "[syng2gfa] scaffold-context evaluated {} local shared syncmer occurrence(s) before exact blunt materialization (weak={}, split={}, min_run={}, scaffold-supported={}, sequence_k={}, sequence-supported={}) in {:.3}s",
            context_filter.shared_occurrences,
            weak_occurrences,
            weak_occurrences,
            frequency_mask.min_shared_run,
            context_filter.scaffold_supported_occurrences,
            frequency_mask.min_sequence_span_bp,
            context_filter.sequence_supported_occurrences,
            context_start.elapsed().as_secs_f64()
        );
        if weak_occurrences > 0 {
            info!(
                "[syng2gfa] scaffold-context split {} local syncmer occurrence(s) into private exact sequence segments",
                weak_occurrences
            );
        }
    }

    let collect_start = Instant::now();
    let path_work: Vec<PathWork> = ranges
        .par_iter()
        .map(|range| {
            collect_range_work(
                index,
                range,
                syncmer_len_u64,
                gap_fill,
                (!masked_syncmers.is_empty()).then_some(&masked_syncmers),
                frequency_mask,
            )
        })
        .collect::<io::Result<Vec<_>>>()?;
    info!(
        "[syng2gfa] collected {} selected path walks for exact blunt materialization in {:.3}s using {} rayon threads",
        path_work.len(),
        collect_start.elapsed().as_secs_f64(),
        rayon::current_num_threads()
    );

    write_exact_blunt_path_work_gfa(
        index,
        path_work,
        writer,
        version,
        &split_occurrences,
        frequency_mask,
    )
}

/// Write syng GFA in either explicit native-overlap form or exact blunt form.
pub fn write_gfa_with_mode<W: Write>(
    index: &SyngIndex,
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
    mode: SyngGfaMode,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    match mode {
        SyngGfaMode::Raw => write_gfa(index, writer, version, gap_fill),
        SyngGfaMode::Blunt => write_gfa_blunt_exact(index, writer, version, gap_fill),
    }
}

pub fn write_range_gfa_with_mode<W: Write>(
    index: &SyngIndex,
    ranges: &[SyngGfaPathRange],
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
    mode: SyngGfaMode,
    frequency_mask: SyngGfaFrequencyMask,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    match mode {
        SyngGfaMode::Raw => {
            write_range_gfa(index, ranges, writer, version, gap_fill, frequency_mask)
        }
        SyngGfaMode::Blunt => {
            write_range_gfa_blunt_exact(index, ranges, writer, version, gap_fill, frequency_mask)
        }
    }
}

pub fn write_range_gfa<W: Write>(
    index: &SyngIndex,
    ranges: &[SyngGfaPathRange],
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
    frequency_mask: SyngGfaFrequencyMask,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    write_range_gfa_internal(index, ranges, writer, version, gap_fill, frequency_mask)
}

fn write_range_gfa_internal<W: Write>(
    index: &SyngIndex,
    ranges: &[SyngGfaPathRange],
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
    frequency_mask: SyngGfaFrequencyMask,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    let syncmer_len = index.syncmer_length_bp();
    info!(
        "[syng2gfa] rendering {} selected path range(s), GFA {}, gap fill: {}",
        ranges.len(),
        version.header_tag(),
        if gap_fill.is_some() {
            "FASTA/AGC"
        } else {
            "Ns"
        }
    );

    unsafe { syng_ffi::impg_syng_suppress_debug() };
    writeln!(writer, "H\tVN:Z:{}", version.header_tag())?;

    let syncmer_len_u64 = syncmer_len as u64;
    let need_local_walks =
        frequency_mask.frequency_filter_enabled() || frequency_mask.shared_context_filter_enabled();
    let local_walks = if need_local_walks {
        let mask_collect_start = Instant::now();
        let walks: Vec<Vec<LocalSyncmerWalkStep>> = ranges
            .par_iter()
            .map(|range| collect_range_syncmer_walk(index, range, syncmer_len_u64))
            .collect::<io::Result<Vec<_>>>()?;
        info!(
            "[syng2gfa] collected {} local syncmer walk(s) for raw-layer masking in {:.3}s using {} rayon threads",
            walks.len(),
            mask_collect_start.elapsed().as_secs_f64(),
            rayon::current_num_threads()
        );
        Some(walks)
    } else {
        None
    };
    let all_local_syncmers = local_walks.as_ref().map(|walks| {
        walks
            .iter()
            .flatten()
            .map(|step| step.node)
            .collect::<FxHashSet<_>>()
    });
    let frequency_masked_syncmers = if frequency_mask.frequency_filter_enabled() {
        frequency_masked_syncmers(
            index,
            all_local_syncmers.as_ref().expect("local syncmers"),
            frequency_mask,
        )?
    } else {
        FxHashSet::default()
    };
    let masked_syncmers = frequency_masked_syncmers.clone();
    let mut split_occurrences = vec![FxHashSet::default(); ranges.len()];
    let mut shared_context_split = 0usize;
    if frequency_mask.shared_context_filter_enabled() {
        let context_start = Instant::now();
        let context_walks = local_walks_after_frequency_mask(
            local_walks.as_ref().expect("local syncmer walks"),
            &frequency_masked_syncmers,
        );
        let context_filter = shared_context_filtered_occurrences(
            &context_walks,
            syncmer_len_u64,
            frequency_mask.min_shared_run,
            frequency_mask.min_sequence_span_bp,
        );
        shared_context_split = context_filter
            .weak_occurrences
            .iter()
            .map(FxHashSet::len)
            .sum::<usize>();
        split_occurrences = context_filter.weak_occurrences;
        info!(
            "[syng2gfa] scaffold-context evaluated {} local shared syncmer occurrence(s) before raw GFA materialization (weak={}, split={}, min_run={}, scaffold-supported={}, sequence_k={}, sequence-supported={}) in {:.3}s",
            context_filter.shared_occurrences,
            shared_context_split,
            shared_context_split,
            frequency_mask.min_shared_run,
            context_filter.scaffold_supported_occurrences,
            frequency_mask.min_sequence_span_bp,
            context_filter.sequence_supported_occurrences,
            context_start.elapsed().as_secs_f64()
        );
        if shared_context_split > 0 {
            info!(
                "[syng2gfa] scaffold-context split {} local syncmer occurrence(s) into private per-occurrence topology before raw GFA materialization (shared occurrences={}, min_run={}, scaffold-supported={}, sequence_k={}, sequence-supported={}) in {:.3}s",
                shared_context_split,
                context_filter.shared_occurrences,
                frequency_mask.min_shared_run,
                context_filter.scaffold_supported_occurrences,
                frequency_mask.min_sequence_span_bp,
                context_filter.sequence_supported_occurrences,
                context_start.elapsed().as_secs_f64()
            );
        }
    };

    let collect_start = Instant::now();
    let path_work: Vec<PathWork> = ranges
        .par_iter()
        .map(|range| {
            collect_range_work(
                index,
                range,
                syncmer_len_u64,
                gap_fill,
                (!masked_syncmers.is_empty()).then_some(&masked_syncmers),
                frequency_mask,
            )
        })
        .collect::<io::Result<Vec<_>>>()?;
    info!(
        "[syng2gfa] collected {} selected path walks in {:.3}s using {} rayon threads",
        path_work.len(),
        collect_start.elapsed().as_secs_f64(),
        rayon::current_num_threads()
    );

    let reduce_start = Instant::now();
    let mut walked_paths: Vec<WalkedPath> = Vec::with_capacity(path_work.len());
    let mut edges: FxHashSet<(String, char, String, char, u32)> = FxHashSet::default();
    let mut gap_interner = GapInterner::new(index.num_syncmer_nodes());
    let mut used_syncmers: FxHashSet<i32> = FxHashSet::default();
    let dominant_repeat_contexts = dominant_local_repeat_contexts(
        &path_work,
        frequency_mask.local_repeat_max_minor,
        frequency_mask.local_repeat_min_dominant_fraction,
    );
    if !dominant_repeat_contexts.is_empty() {
        info!(
            "[syng2gfa] detected {} near-single-copy syncmer node(s) with rare repeated local contexts (max_minor={}, min_dominant_fraction={})",
            dominant_repeat_contexts.len(),
            frequency_mask.local_repeat_max_minor,
            frequency_mask.local_repeat_min_dominant_fraction
        );
    }
    let mut cloned_syncmers: Vec<(String, u32)> = Vec::new();
    let mut n_local_repeat_clones = 0usize;
    let mut n_scaffold_context_clones = 0usize;
    let mut n_skipped = 0usize;
    let mut gaps_filled_with_ns = 0usize;
    let mut gap_occurrences = 0usize;
    let mut gap_occurrence_bp = 0u64;
    let mut n_cut_runs = 0usize;
    let mut n_cut_bp = 0u64;

    for (work_idx, work) in path_work.into_iter().enumerate() {
        if work.skipped {
            warn!(
                "[syng2gfa] path '{}' has no GBWT start info — skipping",
                work.name
            );
            n_skipped += 1;
            continue;
        }

        gaps_filled_with_ns += work.gaps_filled_with_ns;
        gap_occurrences += work.gap_occurrences;
        gap_occurrence_bp += work.gap_occurrence_bp;
        n_cut_runs += work.n_cut_runs;
        n_cut_bp += work.n_cut_bp;

        let work_occurrences = raw_work_syncmer_context_occurrences(&work);
        let split_clone_occurrences = split_occurrences
            .get(work_idx)
            .map(|occurrences| split_shared_context_clone_occurrences(&work, occurrences))
            .unwrap_or_default();
        let local_repeat_clone_occurrences =
            local_repeat_clone_occurrences(&work, &dominant_repeat_contexts);
        let mut clone_occurrences = split_clone_occurrences.clone();
        clone_occurrences.extend(local_repeat_clone_occurrences.iter().copied());
        let mut cloned_syncmer_ids: FxHashMap<u32, String> = FxHashMap::default();
        if !clone_occurrences.is_empty() {
            let mut local_occurrences: Vec<_> = work_occurrences
                .iter()
                .filter(|occ| clone_occurrences.contains(&occ.syncmer.occ))
                .collect();
            local_occurrences.sort_by_key(|occ| (occ.syncmer.abs(), occ.syncmer.occ));
            for occ in local_occurrences {
                let node = occ.syncmer.abs();
                let clone_id = gap_interner.next_segment_id();
                cloned_syncmer_ids.insert(occ.syncmer.occ, clone_id.clone());
                cloned_syncmers.push((clone_id, node));
                if split_clone_occurrences.contains(&occ.syncmer.occ) {
                    n_scaffold_context_clones += 1;
                } else {
                    n_local_repeat_clones += 1;
                }
            }
        }

        let gap_ids: Vec<String> = work
            .gaps
            .into_iter()
            .map(|gap| gap_interner.intern(gap.context, gap.seq))
            .collect();

        for (from, to, overlap) in work.edges {
            if let RawEndpoint::Syncmer(node) = &from {
                let abs = node.abs();
                if !cloned_syncmer_ids.contains_key(&node.occ) {
                    used_syncmers.insert(abs as i32);
                }
            }
            if let RawEndpoint::Syncmer(node) = &to {
                let abs = node.abs();
                if !cloned_syncmer_ids.contains_key(&node.occ) {
                    used_syncmers.insert(abs as i32);
                }
            }
            let (from_id, from_sign) = raw_endpoint_to_gfa(&from, &gap_ids, &cloned_syncmer_ids);
            let (to_id, to_sign) = raw_endpoint_to_gfa(&to, &gap_ids, &cloned_syncmer_ids);
            edges.insert((from_id, from_sign, to_id, to_sign, overlap));
        }

        let mut segments = vec![Vec::new()];
        let mut had_break = false;
        for step in work.steps {
            match step {
                RawPathStep::Syncmer(syncmer) => {
                    let (id, sign) = signed_to_gfa(syncmer.signed);
                    if let Some(clone_id) = cloned_syncmer_ids.get(&syncmer.occ) {
                        segments
                            .last_mut()
                            .expect("path segment")
                            .push(PathStep::Segment {
                                id: clone_id.clone(),
                                sign,
                            });
                    } else {
                        used_syncmers.insert(id);
                        segments
                            .last_mut()
                            .expect("path segment")
                            .push(PathStep::Syncmer(syncmer.signed));
                    }
                }
                RawPathStep::Gap(idx) => {
                    segments
                        .last_mut()
                        .expect("path segment")
                        .push(PathStep::Gap(gap_ids[idx].clone()));
                }
                RawPathStep::Break => {
                    had_break = true;
                    segments.push(Vec::new());
                }
            }
        }
        push_walked_path_segments(&mut walked_paths, work.name, segments, had_break);
    }
    info!(
        "[syng2gfa] reduced selected path walks into {} shared syncmer node(s), {} cloned syncmer segment(s) ({} local-repeat, {} scaffold-context), {} frequency syncmer node(s) removed and {} scaffold-context occurrence(s) split before raw GFA, {} unique gap segment(s), {} edge(s) in {:.3}s",
        used_syncmers.len(),
        cloned_syncmers.len(),
        n_local_repeat_clones,
        n_scaffold_context_clones,
        frequency_masked_syncmers.len(),
        shared_context_split,
        gap_interner.len(),
        edges.len(),
        reduce_start.elapsed().as_secs_f64()
    );

    if gaps_filled_with_ns > 0 {
        warn!(
            "[syng2gfa] filled {} gap occurrence(s) with 'N' ({} bp over paths; {} unique segment bp). \
                 Concatenating these paths will NOT reconstruct the original genome. \
                 Pass --sequence-files <FASTA> to splice in real DNA.",
            gaps_filled_with_ns,
            gap_occurrence_bp,
            gap_interner.total_segment_bp()
        );
    } else if !gap_interner.is_empty() {
        info!(
            "[syng2gfa] spliced {} unique real-DNA gap segment(s) from {} gap occurrence(s), {} unique bp ({} bp over paths)",
            gap_interner.len(),
            gap_occurrences,
            gap_interner.total_segment_bp(),
            gap_occurrence_bp
        );
    }
    if n_cut_runs > 0 {
        info!(
            "[syng2gfa] cut {} N-run(s) ({} bp) out of gap DNA and split emitted paths at those breaks",
            n_cut_runs,
            n_cut_bp
        );
    }

    let write_start = Instant::now();
    let mut used_nodes: Vec<i32> = used_syncmers.into_iter().collect();
    used_nodes.sort_unstable();
    let n_syncmer_segments = used_nodes.len();
    for node_id in used_nodes {
        let mut seq = index.syncmer_seq(node_id);
        seq.make_ascii_uppercase();
        write_segment(writer, &node_id.to_string(), &seq)?;
    }
    for (clone_id, node_id) in &cloned_syncmers {
        let mut seq = index.syncmer_seq(*node_id as i32);
        seq.make_ascii_uppercase();
        write_segment(writer, clone_id, &seq)?;
    }
    for gap in &gap_interner.segments {
        write_segment(writer, &gap.id, &gap.seq)?;
    }
    let n_segments = n_syncmer_segments + cloned_syncmers.len() + gap_interner.len();
    let n_gap_segments = gap_interner.len();

    let mut edge_vec: Vec<(String, char, String, char, u32)> = edges.into_iter().collect();
    edge_vec.sort_unstable();
    for (a_id, a_sign, b_id, b_sign, overlap) in &edge_vec {
        writeln!(writer, "L\t{a_id}\t{a_sign}\t{b_id}\t{b_sign}\t{overlap}M")?;
    }
    let n_links = edge_vec.len();

    let n_paths_emitted = walked_paths.len();
    for walked in &walked_paths {
        match version {
            GfaVersion::V1_0 => write_p_line(writer, &walked.name, &walked.steps)?,
            GfaVersion::V1_1 => {
                if let Some(parts) = parse_pansn(&walked.name) {
                    write_w_line(writer, &parts, &walked.steps)?;
                } else {
                    warn!(
                        "[syng2gfa] path '{}' is not PanSN-style; falling back to P line under GFA 1.1",
                        walked.name
                    );
                    write_p_line(writer, &walked.name, &walked.steps)?;
                }
            }
        }
    }
    info!(
        "[syng2gfa] wrote selected range GFA: {} S, {} L, {} path line(s) in {:.3}s",
        n_segments,
        n_links,
        n_paths_emitted,
        write_start.elapsed().as_secs_f64()
    );

    Ok((
        n_segments,
        n_links,
        n_paths_emitted,
        n_skipped,
        n_gap_segments,
        gap_interner.total_segment_bp(),
    ))
}

#[cfg(test)]
fn bluntify_gfa_bytes(raw: &[u8], version: GfaVersion) -> io::Result<Vec<u8>> {
    if version != GfaVersion::V1_0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--gfa-mode blunt currently requires --gfa-version 1.0 because bluntg preserves P lines, not W lines",
        ));
    }

    let gfa = bluntg::parse_gfa(raw);
    let blunted = bluntg::bluntify_auto(gfa);
    let mut out = Vec::new();
    bluntg::write_gfa(&mut out, &blunted)?;
    Ok(out)
}

#[cfg(test)]
fn bluntify_gfa_bytes_with_stable_clones(
    raw: &[u8],
    version: GfaVersion,
    clone_blunt_sequences: &FxHashMap<String, Vec<u8>>,
) -> io::Result<Vec<u8>> {
    let blunted = bluntify_gfa_bytes(raw, version)?;
    if clone_blunt_sequences.is_empty() {
        return Ok(blunted);
    }
    patch_blunted_clone_segment_sequences(&blunted, clone_blunt_sequences)
}

#[cfg(test)]
fn patch_blunted_clone_segment_sequences(
    blunted: &[u8],
    clone_blunt_sequences: &FxHashMap<String, Vec<u8>>,
) -> io::Result<Vec<u8>> {
    let text = std::str::from_utf8(blunted).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("bluntg emitted non-UTF8 GFA: {err}"),
        )
    })?;
    let mut out = Vec::with_capacity(blunted.len());
    let mut replaced = 0usize;

    for line in text.split_inclusive('\n') {
        let has_newline = line.ends_with('\n');
        let body = if has_newline {
            &line[..line.len() - 1]
        } else {
            line
        };
        let mut fields = body.splitn(4, '\t');
        if fields.next() == Some("S") {
            let id = fields.next().unwrap_or("");
            let _old_seq = fields.next();
            if let Some(seq) = clone_blunt_sequences.get(id) {
                out.extend_from_slice(b"S\t");
                out.extend_from_slice(id.as_bytes());
                out.push(b'\t');
                out.extend_from_slice(seq);
                if let Some(rest) = fields.next() {
                    out.push(b'\t');
                    out.extend_from_slice(rest.as_bytes());
                }
                if has_newline {
                    out.push(b'\n');
                }
                replaced += 1;
                continue;
            }
        }
        out.extend_from_slice(line.as_bytes());
    }

    if replaced != clone_blunt_sequences.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "bluntg output omitted {} private syncmer clone segment(s)",
                clone_blunt_sequences.len() - replaced
            ),
        ));
    }
    Ok(out)
}

/// Pull the DNA between two syncmers from a sequence index.
fn fetch_gap_dna(
    seq_idx: &UnifiedSequenceIndex,
    seq_name: &str,
    gap_start: u64,
    gap_end: u64,
) -> io::Result<Vec<u8>> {
    if gap_end <= gap_start {
        return Ok(Vec::new());
    }
    if gap_start > i32::MAX as u64 || gap_end > i32::MAX as u64 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "gap coordinates exceed i32 range on '{}': [{}, {})",
                seq_name, gap_start, gap_end
            ),
        ));
    }
    seq_idx.fetch_sequence(seq_name, gap_start as i32, gap_end as i32)
}

fn signed_to_gfa(signed: i32) -> (i32, char) {
    if signed >= 0 {
        (signed, '+')
    } else {
        (-signed, '-')
    }
}

fn write_segment<W: Write>(writer: &mut W, id: &str, seq: &[u8]) -> io::Result<()> {
    writer.write_all(b"S\t")?;
    writer.write_all(id.as_bytes())?;
    writer.write_all(b"\t")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_p_line<W: Write>(writer: &mut W, name: &str, steps: &[PathStep]) -> io::Result<()> {
    writer.write_all(b"P\t")?;
    writer.write_all(name.as_bytes())?;
    writer.write_all(b"\t")?;
    for (i, step) in steps.iter().enumerate() {
        if i > 0 {
            writer.write_all(b",")?;
        }
        match step {
            PathStep::Syncmer(signed) => {
                let (id, sign) = signed_to_gfa(*signed);
                write!(writer, "{id}{sign}")?;
            }
            PathStep::Segment { id, sign } => {
                write!(writer, "{id}{sign}")?;
            }
            PathStep::Gap(id) => {
                write!(writer, "{id}+")?;
            }
        }
    }
    writer.write_all(b"\t*\n")?;
    Ok(())
}

fn write_w_line<W: Write>(
    writer: &mut W,
    parts: &PansnParts<'_>,
    steps: &[PathStep],
) -> io::Result<()> {
    write!(
        writer,
        "W\t{}\t{}\t{}\t{}\t{}\t",
        parts.sample, parts.hap, parts.contig, parts.start, parts.end
    )?;
    for step in steps {
        match step {
            PathStep::Syncmer(signed) => {
                let (id, sign) = signed_to_gfa(*signed);
                let arrow = if sign == '+' { '>' } else { '<' };
                write!(writer, "{arrow}{id}")?;
            }
            PathStep::Segment { id, sign } => {
                let arrow = if *sign == '+' { '>' } else { '<' };
                write!(writer, "{arrow}{id}")?;
            }
            PathStep::Gap(id) => {
                write!(writer, ">{id}")?;
            }
        }
    }
    writer.write_all(b"\n")?;
    Ok(())
}

/// Convenience: load a syng index from `prefix`, optionally build a
/// sequence index from `sequence_files`, and write the GFA to `out_path`
/// (`-` for stdout).
pub fn run(
    prefix: &str,
    out_path: &str,
    version: GfaVersion,
    sequence_files: &[String],
    mode: SyngGfaMode,
) -> io::Result<()> {
    info!("[syng2gfa] loading syng index from prefix '{}'", prefix);
    // SyncmerParams in load() is overridden by the .syng.meta sidecar, so any
    // placeholder works.
    let index = SyngIndex::load(prefix, crate::syng::SyncmerParams::default())?;

    let seq_idx: Option<UnifiedSequenceIndex> = if sequence_files.is_empty() {
        None
    } else {
        info!(
            "[syng2gfa] building sequence index for {} file(s) (gap fill)",
            sequence_files.len()
        );
        Some(UnifiedSequenceIndex::from_files(sequence_files)?)
    };

    let report = |s: usize,
                  l: usize,
                  p: usize,
                  skipped: usize,
                  gaps: usize,
                  gap_bp: u64,
                  where_: &str| {
        info!(
                "[syng2gfa] wrote {} S ({} syncmer + {} gap, {} gap bp), {} L, {} P/W lines ({} paths skipped) to {}",
                s,
                s - gaps,
                gaps,
                gap_bp,
                l,
                p,
                skipped,
                where_,
            );
    };

    if out_path == "-" {
        let stdout = io::stdout();
        let mut w = BufWriter::new(stdout.lock());
        let (s, l, p, skipped, gaps, gap_bp) =
            write_gfa_with_mode(&index, &mut w, version, seq_idx.as_ref(), mode)?;
        w.flush()?;
        report(
            s,
            l,
            p,
            skipped,
            gaps,
            gap_bp,
            &format!("stdout ({})", mode.label()),
        );
    } else {
        let file = std::fs::File::create(out_path)?;
        let mut w = BufWriter::new(file);
        let (s, l, p, skipped, gaps, gap_bp) =
            write_gfa_with_mode(&index, &mut w, version, seq_idx.as_ref(), mode)?;
        w.flush()?;
        report(
            s,
            l,
            p,
            skipped,
            gaps,
            gap_bp,
            &format!("'{}' ({})", out_path, mode.label()),
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blunt_path_map(gfa: &[u8]) -> FxHashMap<String, String> {
        let text = std::str::from_utf8(gfa).expect("GFA should be UTF-8");
        crate::resolution::path_sequences(text)
            .expect("blunt GFA path sequences")
            .into_iter()
            .collect()
    }

    fn assert_stable_private_clone_fixture(
        baseline_raw: &[u8],
        split_raw: &[u8],
        clone_id: &str,
        clone_blunt_seq: &[u8],
        mismatching_path: &str,
    ) {
        let baseline = bluntify_gfa_bytes(baseline_raw, GfaVersion::V1_0).unwrap();
        let unpatched_split = bluntify_gfa_bytes(split_raw, GfaVersion::V1_0).unwrap();
        let baseline_paths = blunt_path_map(&baseline);
        let unpatched_split_paths = blunt_path_map(&unpatched_split);
        assert_ne!(
            unpatched_split_paths.get(mismatching_path),
            baseline_paths.get(mismatching_path),
            "the reduced fixture should reproduce the pre-fix private-copy spelling drift"
        );

        let mut clone_sequences = FxHashMap::default();
        clone_sequences.insert(clone_id.to_string(), clone_blunt_seq.to_vec());
        let patched =
            bluntify_gfa_bytes_with_stable_clones(split_raw, GfaVersion::V1_0, &clone_sequences)
                .unwrap();
        let patched_paths = blunt_path_map(&patched);
        assert_eq!(
            patched_paths, baseline_paths,
            "private clone should preserve every post-blunt path spelling from the shared-node baseline"
        );
    }

    fn context_test_work(name: &str, nodes: &[i32]) -> PathWork {
        let steps: Vec<RawPathStep> = nodes
            .iter()
            .enumerate()
            .map(|(occ, signed)| {
                RawPathStep::Syncmer(RawSyncmerOcc {
                    signed: *signed,
                    occ: occ as u32,
                    pos: occ as u64,
                    left_clip: 0,
                    right_clip: 0,
                })
            })
            .collect();
        PathWork {
            name: name.to_string(),
            gaps: Vec::new(),
            steps,
            edges: Vec::new(),
            gaps_filled_with_ns: 0,
            gap_occurrences: 0,
            gap_occurrence_bp: 0,
            n_cut_runs: 0,
            n_cut_bp: 0,
            skipped: false,
        }
    }

    #[test]
    fn test_syng_gfa_mode_parse_defaults_syng_to_blunt() {
        assert_eq!(SyngGfaMode::parse("syng").unwrap(), SyngGfaMode::Blunt);
        assert_eq!(
            SyngGfaMode::parse("syng:blunt").unwrap(),
            SyngGfaMode::Blunt
        );
        assert_eq!(SyngGfaMode::parse("bluntg").unwrap(), SyngGfaMode::Blunt);
        assert_eq!(SyngGfaMode::parse("raw").unwrap(), SyngGfaMode::Raw);
        assert_eq!(SyngGfaMode::parse("syng:raw").unwrap(), SyngGfaMode::Raw);
    }

    #[test]
    fn test_bluntify_gfa_bytes_handles_variable_overlaps() {
        let raw = b"\
H\tVN:Z:1.0\n\
S\t1\tACGTAC\n\
S\t2\tTACGGA\n\
S\t3\tGGAACC\n\
L\t1\t+\t2\t+\t3M\n\
L\t2\t+\t3\t+\t1M\n\
P\tp\t1+,2+,3+\t3M,1M\n";
        let out = bluntify_gfa_bytes(raw, GfaVersion::V1_0).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("L\t1\t+\t2\t+\t0M"), "{out}");
        assert!(out.contains("L\t2\t+\t3\t+\t0M"), "{out}");
        assert!(out.contains("P\tp\t1+,2+,3+\t0M,0M"), "{out}");
        assert!(!out.contains("\t3M"), "{out}");
        assert!(!out.contains("\t1M"), "{out}");
    }

    #[test]
    fn test_private_clone_blunt_spelling_uses_shared_trim_coordinates() {
        assert_stable_private_clone_fixture(
            b"\
H\tVN:Z:1.0\n\
S\t1\tAA\n\
S\t2\tCCGGTT\n\
S\t3\tTT\n\
S\t4\tGGGGCC\n\
L\t1\t+\t2\t+\t0M\n\
L\t2\t+\t3\t+\t0M\n\
L\t4\t+\t2\t+\t4M\n\
P\ttarget\t1+,2+,3+\t*\n\
P\tcontext\t4+,2+\t*\n",
            b"\
H\tVN:Z:1.0\n\
S\t1\tAA\n\
S\t2\tCCGGTT\n\
S\t3\tTT\n\
S\t4\tGGGGCC\n\
S\t5\tCCGGTT\n\
L\t1\t+\t5\t+\t0M\n\
L\t5\t+\t3\t+\t0M\n\
L\t4\t+\t2\t+\t4M\n\
P\ttarget\t1+,5+,3+\t*\n\
P\tcontext\t4+,2+\t*\n",
            "5",
            b"GGTT",
            "target",
        );
    }

    #[test]
    fn test_sequence_k_191_chr1_reduced_fixture_preserves_blunt_path_spellings() {
        assert_stable_private_clone_fixture(
            b"\
H\tVN:Z:1.0\n\
S\t1\tAA\n\
S\t2\tCCGGTT\n\
S\t3\tTT\n\
S\t4\tGGGGCC\n\
L\t1\t+\t2\t+\t0M\n\
L\t2\t+\t3\t+\t0M\n\
L\t4\t+\t2\t+\t4M\n\
P\tchr1-sanity#0#chr1:0-10\t1+,2+,3+\t*\n\
P\tchr1-sanity#0#chr1-context:0-8\t4+,2+\t*\n",
            b"\
H\tVN:Z:1.0\n\
S\t1\tAA\n\
S\t2\tCCGGTT\n\
S\t3\tTT\n\
S\t4\tGGGGCC\n\
S\t5\tCCGGTT\n\
L\t1\t+\t5\t+\t0M\n\
L\t5\t+\t3\t+\t0M\n\
L\t4\t+\t2\t+\t4M\n\
P\tchr1-sanity#0#chr1:0-10\t1+,5+,3+\t*\n\
P\tchr1-sanity#0#chr1-context:0-8\t4+,2+\t*\n",
            "5",
            b"GGTT",
            "chr1-sanity#0#chr1:0-10",
        );
    }

    #[test]
    fn test_sequence_k_311_c4_reduced_fixture_preserves_blunt_path_spellings() {
        assert_stable_private_clone_fixture(
            b"\
H\tVN:Z:1.0\n\
S\t1\tAA\n\
S\t2\tCCGGTT\n\
S\t3\tTT\n\
S\t4\tTTGGGG\n\
L\t1\t+\t2\t+\t0M\n\
L\t2\t+\t3\t+\t0M\n\
L\t2\t+\t4\t+\t4M\n\
P\tC4#0#chr6:31744284-31744294\t1+,2+,3+\t*\n\
P\tC4#0#chr6-context:31744284-31744292\t2+,4+\t*\n",
            b"\
H\tVN:Z:1.0\n\
S\t1\tAA\n\
S\t2\tCCGGTT\n\
S\t3\tTT\n\
S\t4\tTTGGGG\n\
S\t5\tCCGGTT\n\
L\t1\t+\t5\t+\t0M\n\
L\t5\t+\t3\t+\t0M\n\
L\t2\t+\t4\t+\t4M\n\
P\tC4#0#chr6:31744284-31744294\t1+,5+,3+\t*\n\
P\tC4#0#chr6-context:31744284-31744292\t2+,4+\t*\n",
            "5",
            b"CCGG",
            "C4#0#chr6:31744284-31744294",
        );
    }

    #[test]
    fn test_cut_ns_splits_gap_segments_and_paths() {
        let mut gaps = Vec::new();
        let mut steps = Vec::new();
        let policy = SyngGfaFrequencyMask {
            cut_n_gaps: true,
            cut_n_min_run: 2,
            ..SyngGfaFrequencyMask::disabled()
        };
        let pushed = push_gap_with_n_cutting(
            &mut gaps,
            &mut steps,
            GapContext::PathOnly {
                path_idx: 0,
                start: 0,
                end: 10,
                strand: '+',
            },
            b"AAANNCCNNN".to_vec(),
            policy,
        );
        assert_eq!(pushed.pieces, vec![0, 1]);
        assert!(!pushed.starts_with_cut);
        assert!(pushed.ends_with_cut);
        assert_eq!(pushed.cut_runs, 2);
        assert_eq!(pushed.cut_bp, 5);
        assert_eq!(gaps[0].seq, b"AAA");
        assert_eq!(gaps[1].seq, b"CC");
        assert!(matches!(steps[0], RawPathStep::Gap(0)));
        assert!(matches!(steps[1], RawPathStep::Break));
        assert!(matches!(steps[2], RawPathStep::Gap(1)));
        assert!(matches!(steps[3], RawPathStep::Break));
    }

    #[test]
    fn test_walked_path_segments_are_named_when_split() {
        let mut walked_paths = Vec::new();
        push_walked_path_segments(
            &mut walked_paths,
            "sample#0#chr1:0-100".to_string(),
            vec![
                vec![PathStep::Syncmer(1)],
                Vec::new(),
                vec![PathStep::Gap("g1".to_string()), PathStep::Syncmer(-2)],
            ],
            true,
        );
        assert_eq!(walked_paths.len(), 2);
        assert_eq!(walked_paths[0].name, "sample#0#chr1:0-100|part1");
        assert_eq!(walked_paths[1].name, "sample#0#chr1:0-100|part2");
        assert_eq!(walked_paths[0].steps.len(), 1);
        assert_eq!(walked_paths[1].steps.len(), 2);
    }

    #[test]
    fn test_range_gfa_frequency_filter_drops_repetitive_syncmers_before_raw_gfa() {
        let params = crate::syng::SyncmerParams {
            k: 8,
            w: 55,
            seed: 7,
        };
        let seq = b"ACGT".repeat(300);
        let index = crate::syng::SyngIndex::build(
            params,
            vec![
                ("sample#0#chr1".to_string(), seq.clone()),
                ("sample#1#chr1".to_string(), seq.clone()),
            ]
            .into_iter(),
        );
        let ranges = vec![
            SyngGfaPathRange {
                path_idx: 0,
                name: "sample#0#chr1:0-1200".to_string(),
                start: 0,
                end: 1200,
                strand: '+',
            },
            SyngGfaPathRange {
                path_idx: 1,
                name: "sample#1#chr1:0-1200".to_string(),
                start: 0,
                end: 1200,
                strand: '+',
            },
        ];
        let mask = SyngGfaFrequencyMask {
            drop_top_fraction: 0.5,
            max_occurrences: None,
            min_shared_run: 1,
            min_sequence_span_bp: 0,
            local_repeat_max_minor: 0,
            local_repeat_min_dominant_fraction: DEFAULT_GFA_LOCAL_REPEAT_MIN_DOMINANT_FRACTION,
            cut_n_gaps: false,
            cut_n_min_run: DEFAULT_GFA_CUT_N_MIN_RUN,
        };
        let all_local_syncmers = ranges.iter().fold(FxHashSet::default(), |mut acc, range| {
            acc.extend(collect_range_syncmer_nodes(&index, range).unwrap());
            acc
        });
        let masked = frequency_masked_syncmers(&index, &all_local_syncmers, mask).unwrap();
        assert!(
            !masked.is_empty(),
            "synthetic repeat should produce masked syncmers"
        );

        let mut out = Vec::new();
        write_range_gfa_with_mode(
            &index,
            &ranges,
            &mut out,
            GfaVersion::V1_0,
            None,
            SyngGfaMode::Raw,
            mask,
        )
        .unwrap();
        let gfa = String::from_utf8(out).unwrap();
        assert!(
            !gfa.contains("\nS\thf"),
            "high-frequency syncmers should be filtered, not cloned into raw GFA:\n{gfa}"
        );
        let segment_ids: FxHashSet<String> = gfa
            .lines()
            .filter_map(|line| {
                let mut fields = line.split('\t');
                (fields.next() == Some("S")).then(|| fields.next().unwrap_or("").to_string())
            })
            .collect();
        let path_step_ids: FxHashSet<String> = gfa
            .lines()
            .filter_map(|line| {
                let mut fields = line.split('\t');
                (fields.next() == Some("P")).then(|| fields.nth(1).unwrap_or(""))
            })
            .flat_map(|steps| steps.split(','))
            .map(|step| step.trim_end_matches(|c| c == '+' || c == '-').to_string())
            .collect();
        for node in masked {
            let id = node.to_string();
            assert!(
                !segment_ids.contains(&id),
                "filtered syncmer {id} should not be emitted as an S line:\n{gfa}"
            );
            assert!(
                !path_step_ids.contains(&id),
                "filtered syncmer {id} should not appear in P-line walks:\n{gfa}"
            );
        }
    }

    fn local_walk(items: &[(u32, u64)]) -> Vec<LocalSyncmerWalkStep> {
        items
            .iter()
            .map(|&(node, pos)| LocalSyncmerWalkStep { node, pos })
            .collect()
    }

    #[test]
    fn test_scaffold_context_retains_chain_members_and_splits_off_chain_occurrences() {
        let walks = vec![
            local_walk(&[(40, 0), (41, 100), (42, 200), (42, 5_000)]),
            local_walk(&[(40, 0), (41, 100), (42, 200)]),
            local_walk(&[(42, 7_000)]),
        ];

        let filtered = shared_context_filtered_occurrences(&walks, 31, 3, 0);
        assert_eq!(filtered.shared_occurrences, 8);
        assert_eq!(
            filtered.scaffold_supported_occurrences, 6,
            "the three-anchor scaffold chain should retain both path copies"
        );

        for occ in [0, 1, 2] {
            assert!(
                !filtered.weak_occurrences[0].contains(&occ),
                "path 0 occurrence {occ} is part of the retained scaffold chain"
            );
            assert!(
                !filtered.weak_occurrences[1].contains(&occ),
                "path 1 occurrence {occ} is part of the retained scaffold chain"
            );
        }
        assert!(
            filtered.weak_occurrences[0].contains(&3),
            "the extra path-0 copy of node 42 is off-chain and must be private"
        );
        assert!(
            filtered.weak_occurrences[2].contains(&0),
            "the singleton path-2 copy of node 42 is off-chain and must be private"
        );
    }

    #[test]
    fn test_scaffold_context_keeps_indel_shifted_colinear_runs() {
        let walks = vec![
            local_walk(&[(1, 0), (2, 100), (3, 220), (4, 320)]),
            local_walk(&[(1, 0), (2, 100), (3, 260), (4, 360)]),
        ];

        let filtered = shared_context_filtered_occurrences(&walks, 31, 4, 0);
        assert_eq!(
            filtered.scaffold_supported_occurrences, 8,
            "SweepGA scaffold membership should retain all four shifted anchors on both paths"
        );
        assert!(
            filtered.weak_occurrences.iter().all(|set| set.is_empty()),
            "a colinear run with an indel-sized diagonal shift must not be split"
        );
    }

    #[test]
    fn test_sequence_k_filter_uses_exact_span_or_scaffold_support() {
        let walks = vec![
            local_walk(&[
                (1, 0),
                (2, 32),
                (3, 96),
                (10, 5_200),
                (11, 5_232),
                (12, 5_264),
            ]),
            local_walk(&[
                (4, 0),
                (2, 5_000),
                (5, 5_064),
                (10, 5_200),
                (11, 5_232),
                (12, 5_264),
            ]),
        ];

        let sequence_supported = shared_context_filtered_occurrences(&walks, 63, 1, 127);
        assert!(
            sequence_supported.weak_occurrences[0].contains(&1)
                && sequence_supported.weak_occurrences[1].contains(&1),
            "isolated shared node 2 should be filtered"
        );
        for occ in [3, 4, 5] {
            assert!(
                !sequence_supported.weak_occurrences[0].contains(&occ)
                    && !sequence_supported.weak_occurrences[1].contains(&occ),
                "occurrence {occ} is part of a repeated exact 127 bp syncmer span"
            );
        }
        assert_eq!(sequence_supported.sequence_supported_occurrences, 6);

        let too_short = shared_context_filtered_occurrences(&walks, 63, 1, 191);
        for occ in [3, 4, 5] {
            assert!(
                too_short.weak_occurrences[0].contains(&occ)
                    && too_short.weak_occurrences[1].contains(&occ),
                "127 bp evidence should not satisfy sequence-k=191"
            );
        }

        let rescued_by_scaffold = shared_context_filtered_occurrences(&walks, 63, 3, 191);
        for occ in [3, 4, 5] {
            assert!(
                !rescued_by_scaffold.weak_occurrences[0].contains(&occ)
                    && !rescued_by_scaffold.weak_occurrences[1].contains(&occ),
                "min-run=3 should rescue the consecutive syncmer scaffold"
            );
        }
        assert!(
            rescued_by_scaffold.weak_occurrences[0].contains(&1)
                && rescued_by_scaffold.weak_occurrences[1].contains(&1),
            "isolated shared node 2 is supported by neither sequence-k nor min-run"
        );
    }

    #[test]
    fn test_scaffold_context_split_clones_only_weak_occurrences() {
        let work = context_test_work("p", &[1, 2, -2, 3, 2]);
        let split_occurrences = [2].into_iter().collect::<FxHashSet<_>>();
        let clone_occurrences = split_shared_context_clone_occurrences(&work, &split_occurrences);
        assert_eq!(clone_occurrences.len(), 1);
        assert!(
            clone_occurrences.contains(&2),
            "only the weak off-chain occurrence should be split"
        );
        assert!(
            !clone_occurrences.contains(&1) && !clone_occurrences.contains(&4),
            "other occurrences of the same node can remain shared when scaffold-supported"
        );
    }

    #[test]
    fn test_local_repeat_context_clones_only_repeated_minor_context() {
        let works = vec![
            context_test_work("p0", &[1, 2, 3]),
            context_test_work("p1", &[1, 2, 3]),
            context_test_work("p2", &[1, 2, 3, 4, 2, 5]),
            context_test_work("p3", &[1, 2, 6]),
        ];
        let dominant = dominant_local_repeat_contexts(&works, 2, 0.60);
        assert!(dominant.contains_key(&2));

        let clones = local_repeat_clone_occurrences(&works[2], &dominant);
        assert!(
            clones.contains(&4),
            "the second copy of node 2 in p2 should be split"
        );
        assert!(
            !clones.contains(&1),
            "the dominant first copy of node 2 in p2 should stay shared"
        );

        let rare_single_copy = local_repeat_clone_occurrences(&works[3], &dominant);
        assert!(
            rare_single_copy.is_empty(),
            "single-copy rare contexts should remain ordinary variation, not cloned"
        );
    }

    #[test]
    fn test_gap_interner_reuses_only_same_sequence_and_context() {
        let mut gaps = GapInterner::new(10);

        let first = gaps.intern(GapContext::Prefix { right: 7 }, b"ACGT".to_vec());
        let same = gaps.intern(GapContext::Prefix { right: 7 }, b"ACGT".to_vec());
        assert_eq!(first, same);
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps.total_segment_bp(), 4);

        let different_context = gaps.intern(GapContext::Prefix { right: 8 }, b"ACGT".to_vec());
        assert_ne!(first, different_context);
        assert_eq!(gaps.len(), 2);

        let different_sequence = gaps.intern(GapContext::Prefix { right: 7 }, b"ACGA".to_vec());
        assert_ne!(first, different_sequence);
        assert_eq!(gaps.len(), 3);
    }

    #[test]
    fn test_extra_segment_ids_are_numeric_and_unique() {
        let mut gaps = GapInterner::new(10);

        let clone_id = gaps.next_segment_id();
        let first_gap = gaps.intern(GapContext::Prefix { right: 7 }, b"ACGT".to_vec());
        let same_gap = gaps.intern(GapContext::Prefix { right: 7 }, b"ACGT".to_vec());
        let second_gap = gaps.intern(GapContext::Suffix { left: 7 }, b"ACGT".to_vec());

        assert_eq!(clone_id, "11");
        assert_eq!(first_gap, "12");
        assert_eq!(same_gap, first_gap);
        assert_eq!(second_gap, "13");
        for id in [&clone_id, &first_gap, &second_gap] {
            assert!(id.bytes().all(|b| b.is_ascii_digit()));
        }
    }
}
