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
//! repeated sequence elsewhere. High-frequency syncmer occurrences are
//! private-split unless they are part of a configured supported run. The CLI
//! default is blunt mode, which materializes exact source-spelling 0M paths
//! directly; raw mode is the explicit syng-native overlap graph.

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
/// Default minimum high-frequency run that rescues otherwise private
/// high-frequency syncmer occurrences. High-frequency singleton reuse is
/// private by default; long collinear runs remain shared.
pub const DEFAULT_GFA_HIGH_FREQ_MIN_RUN: usize = 10;
/// Default exact sequence span for rescuing high-frequency syncmer occurrences.
/// A 1 kb default matches the C4 high-frequency sweep setting that first
/// materially reduced private-copy explosion while keeping isolated singleton
/// reuse private.
pub const DEFAULT_GFA_HIGH_FREQ_MIN_SEQUENCE_SPAN_BP: usize = 1_000;
/// Clone rare repeated-copy contexts for otherwise dominant local syncmers.
/// This catches single-syncmer loops where a mostly single-copy node appears
/// a second time in a selected path and glues two paralogous contexts together.
pub const DEFAULT_GFA_LOCAL_REPEAT_MAX_MINOR: u32 = 2;
pub const DEFAULT_GFA_LOCAL_REPEAT_MIN_DOMINANT_FRACTION: f64 = 0.80;
/// Local syncmers that are repeatedly reused within carrying paths and spread
/// over scaffold-scale path distance are not allowed to act as scaffold glue.
/// They are split per occurrence instead of removed, preserving path spelling
/// while avoiding dispersed repeat all-pairs glue.
pub const DEFAULT_GFA_SCAFFOLD_GLUE_MIN_OCC_PER_PATH_RATIO: f64 = 2.0;
pub const DEFAULT_GFA_SCAFFOLD_GLUE_MIN_OCCURRENCES: u32 = 64;
pub const DEFAULT_GFA_SCAFFOLD_GLUE_MIN_DISPERSION_BP: u64 =
    crate::syng_transitive::DEFAULT_EXTEND_BUDGET_BP;
/// Default minimum N-run length that splits syng GFA paths when N cutting is
/// enabled. A value of 1 means any ambiguous base breaks the emitted path.
pub const DEFAULT_GFA_CUT_N_MIN_RUN: usize = 1;

/// Frequency-aware syncmer node sharing policy for syng GFA materialization.
///
/// High-frequency syncmers are selected by node-level frequency, but the
/// default policy is occurrence-level: unsupported occurrences are emitted as
/// private clones while occurrences in long supported runs stay shared. A
/// compatibility flag keeps the historical node-removal behavior available for
/// debugging and comparisons. Rare repeated local contexts are still emitted as
/// per-occurrence clones so ordinary variation is not erased while obvious
/// single-copy glue is split.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SyngGfaFrequencyMask {
    pub drop_top_fraction: f64,
    pub max_occurrences: Option<u32>,
    pub high_freq_min_run: usize,
    pub high_freq_min_sequence_span_bp: usize,
    pub run_aware_frequency_mask: bool,
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
            high_freq_min_run: 0,
            high_freq_min_sequence_span_bp: 0,
            run_aware_frequency_mask: true,
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
            high_freq_min_run: DEFAULT_GFA_HIGH_FREQ_MIN_RUN,
            high_freq_min_sequence_span_bp: DEFAULT_GFA_HIGH_FREQ_MIN_SEQUENCE_SPAN_BP,
            run_aware_frequency_mask: true,
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

    fn run_aware_frequency_filter_enabled(self) -> bool {
        self.frequency_filter_enabled() && self.run_aware_frequency_mask
    }

    fn shared_run_filter_enabled(self) -> bool {
        self.min_shared_run > 1
    }

    fn sequence_context_filter_enabled(self) -> bool {
        self.min_sequence_span_bp > 0
    }

    fn shared_context_filter_enabled(self) -> bool {
        self.shared_run_filter_enabled()
            || self.sequence_context_filter_enabled()
            || self.high_freq_min_run > 0
            || self.high_freq_min_sequence_span_bp > 0
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
struct PairAnchorOccurrence {
    query_occ: u32,
    target_occ: u32,
    node: u32,
    query_pos: u64,
    target_pos: u64,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct PairAnchorOccurrenceKey {
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

#[derive(Clone, Copy, Debug)]
struct SharedWalkOccurrence {
    occ: u32,
    node: u32,
    pos: u64,
}

#[derive(Clone, Copy, Debug)]
struct ScaffoldRunSeed {
    path_idx: usize,
    start: usize,
}

const SCAFFOLD_CONTEXT_GAP_BP: u64 = crate::syng_transitive::DEFAULT_EXTEND_BUDGET_BP;
const MAX_SCAFFOLD_SIGNATURE_SEED_PAIRS: usize = if cfg!(test) { 32 } else { 1_000_000 };
const MAX_SCAFFOLD_TOTAL_CANDIDATE_ANCHORS: usize = if cfg!(test) { 10_000 } else { 250_000 };

#[derive(Clone, Copy, Debug)]
struct LocalSyncmerNodeSpectrum {
    node: u32,
    total_occurrences: u32,
    carrier_paths: u32,
    max_occurrences_per_path: u32,
    max_path_span_bp: u64,
}

impl LocalSyncmerNodeSpectrum {
    fn occ_per_path_ratio(self) -> f64 {
        self.total_occurrences as f64 / self.carrier_paths.max(1) as f64
    }
}

#[derive(Default)]
struct LocalSyncmerNodeAccumulator {
    total_occurrences: u32,
    carrier_paths: u32,
    max_occurrences_per_path: u32,
    max_path_span_bp: u64,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct LocalSyncmerCandidateSummary {
    nodes: usize,
    occurrences: usize,
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
    scaffold_candidate_anchors: usize,
    scaffold_dense_signatures: usize,
    path_copy_filtered_nodes: usize,
    path_copy_filtered_occurrences: usize,
    path_copy_supported_occurrences: usize,
    path_copy_run_supported_occurrences: usize,
    path_copy_sequence_supported_occurrences: usize,
    path_copy_private_split_occurrences: usize,
    path_copy_scaffold_candidate_anchors: usize,
    path_copy_scaffold_dense_signatures: usize,
}

#[derive(Default)]
struct HighFrequencyOccurrencePolicyResult {
    split_occurrences: Vec<FxHashSet<u32>>,
    supported_occurrence_keys: FxHashSet<LocalOccurrenceKey>,
    selected_nodes: usize,
    total_occurrences: usize,
    supported_occurrences: usize,
    private_split_occurrences: usize,
    run_supported_occurrences: usize,
    sequence_supported_occurrences: usize,
    scaffold_candidate_anchors: usize,
    scaffold_dense_signatures: usize,
}

struct FilteredLocalWalks {
    walks: Vec<Vec<LocalSyncmerWalkStep>>,
    original_occurrences: Vec<Vec<u32>>,
}

struct RangeMaskEvaluation {
    filtered_syncmers: FxHashSet<u32>,
    split_occurrences: Vec<FxHashSet<u32>>,
    shared_context_split: usize,
    high_frequency_split: usize,
}

impl RangeMaskEvaluation {
    fn empty(path_count: usize) -> Self {
        Self {
            filtered_syncmers: FxHashSet::default(),
            split_occurrences: vec![FxHashSet::default(); path_count],
            shared_context_split: 0,
            high_frequency_split: 0,
        }
    }
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

fn local_syncmer_node_spectrum(
    walks: &[Vec<LocalSyncmerWalkStep>],
    shared_nodes: &FxHashSet<u32>,
) -> Vec<LocalSyncmerNodeSpectrum> {
    if walks.is_empty() || shared_nodes.is_empty() {
        return Vec::new();
    }

    let mut accumulators: FxHashMap<u32, LocalSyncmerNodeAccumulator> = FxHashMap::default();
    for walk in walks {
        let mut path_occurrences: FxHashMap<u32, (u32, u64, u64)> = FxHashMap::default();
        for step in walk {
            if shared_nodes.contains(&step.node) {
                let entry = path_occurrences
                    .entry(step.node)
                    .or_insert((0, step.pos, step.pos));
                entry.0 = entry.0.saturating_add(1);
                entry.1 = entry.1.min(step.pos);
                entry.2 = entry.2.max(step.pos);
            }
        }
        for (node, (count, min_pos, max_pos)) in path_occurrences {
            let acc = accumulators.entry(node).or_default();
            acc.total_occurrences = acc.total_occurrences.saturating_add(count);
            acc.carrier_paths = acc.carrier_paths.saturating_add(1);
            acc.max_occurrences_per_path = acc.max_occurrences_per_path.max(count);
            acc.max_path_span_bp = acc.max_path_span_bp.max(max_pos.saturating_sub(min_pos));
        }
    }

    let mut spectrum = accumulators
        .into_iter()
        .map(|(node, acc)| LocalSyncmerNodeSpectrum {
            node,
            total_occurrences: acc.total_occurrences,
            carrier_paths: acc.carrier_paths,
            max_occurrences_per_path: acc.max_occurrences_per_path,
            max_path_span_bp: acc.max_path_span_bp,
        })
        .collect::<Vec<_>>();
    spectrum.sort_by_key(|stat| stat.node);
    spectrum
}

fn is_dispersed_scaffold_glue_candidate(
    stat: LocalSyncmerNodeSpectrum,
    min_occ_per_path_ratio: f64,
) -> bool {
    stat.total_occurrences >= DEFAULT_GFA_SCAFFOLD_GLUE_MIN_OCCURRENCES
        && stat.max_occurrences_per_path >= 2
        && stat.occ_per_path_ratio() >= min_occ_per_path_ratio
        && stat.max_path_span_bp >= DEFAULT_GFA_SCAFFOLD_GLUE_MIN_DISPERSION_BP
}

fn dispersed_scaffold_glue_summary(
    spectrum: &[LocalSyncmerNodeSpectrum],
    min_occ_per_path_ratio: f64,
) -> LocalSyncmerCandidateSummary {
    spectrum.iter().fold(
        LocalSyncmerCandidateSummary::default(),
        |mut summary, stat| {
            if is_dispersed_scaffold_glue_candidate(*stat, min_occ_per_path_ratio) {
                summary.nodes += 1;
                summary.occurrences += stat.total_occurrences as usize;
            }
            summary
        },
    )
}

fn fixed_selected_path_factor_summary(
    spectrum: &[LocalSyncmerNodeSpectrum],
    selected_paths: usize,
    factor: f64,
) -> LocalSyncmerCandidateSummary {
    if selected_paths == 0 {
        return LocalSyncmerCandidateSummary::default();
    }
    let max_occurrences = ((selected_paths as f64) * factor).ceil() as u32;
    spectrum.iter().fold(
        LocalSyncmerCandidateSummary::default(),
        |mut summary, stat| {
            if stat.total_occurrences > max_occurrences {
                summary.nodes += 1;
                summary.occurrences += stat.total_occurrences as usize;
            }
            summary
        },
    )
}

fn spectrum_scaffold_glue_nodes(spectrum: &[LocalSyncmerNodeSpectrum]) -> FxHashSet<u32> {
    spectrum
        .iter()
        .filter_map(|stat| {
            is_dispersed_scaffold_glue_candidate(
                *stat,
                DEFAULT_GFA_SCAFFOLD_GLUE_MIN_OCC_PER_PATH_RATIO,
            )
            .then_some(stat.node)
        })
        .collect()
}

fn format_spectrum_bins<F>(
    spectrum: &[LocalSyncmerNodeSpectrum],
    labels: &[&str],
    classify: F,
) -> String
where
    F: Fn(LocalSyncmerNodeSpectrum) -> usize,
{
    let mut bins = vec![LocalSyncmerCandidateSummary::default(); labels.len()];
    for stat in spectrum {
        let idx = classify(*stat).min(labels.len().saturating_sub(1));
        bins[idx].nodes += 1;
        bins[idx].occurrences += stat.total_occurrences as usize;
    }
    labels
        .iter()
        .zip(bins)
        .map(|(label, summary)| format!("{label}={}/{}", summary.nodes, summary.occurrences))
        .collect::<Vec<_>>()
        .join(",")
}

fn format_local_spectrum_ratio_bins(spectrum: &[LocalSyncmerNodeSpectrum]) -> String {
    format_spectrum_bins(
        spectrum,
        &["<=1.0", "<=1.25", "<=2.0", "<=4.0", "<=8.0", ">8.0"],
        |stat| {
            let ratio = stat.occ_per_path_ratio();
            if ratio <= 1.0 {
                0
            } else if ratio <= 1.25 {
                1
            } else if ratio <= 2.0 {
                2
            } else if ratio <= 4.0 {
                3
            } else if ratio <= 8.0 {
                4
            } else {
                5
            }
        },
    )
}

fn format_local_spectrum_max_copy_bins(spectrum: &[LocalSyncmerNodeSpectrum]) -> String {
    format_spectrum_bins(spectrum, &["1", "2", "3-4", "5-8", ">8"], |stat| match stat
        .max_occurrences_per_path
    {
        0 | 1 => 0,
        2 => 1,
        3 | 4 => 2,
        5..=8 => 3,
        _ => 4,
    })
}

fn format_local_spectrum_dispersion_bins(spectrum: &[LocalSyncmerNodeSpectrum]) -> String {
    format_spectrum_bins(
        spectrum,
        &["0", "<1kb", "1-10kb", "10-100kb", ">100kb"],
        |stat| match stat.max_path_span_bp {
            0 => 0,
            1..=999 => 1,
            1_000..=9_999 => 2,
            10_000..=99_999 => 3,
            _ => 4,
        },
    )
}

fn format_candidate_summary(summary: LocalSyncmerCandidateSummary) -> String {
    format!("{}/{}", summary.nodes, summary.occurrences)
}

fn log_local_syncmer_spectrum(spectrum: &[LocalSyncmerNodeSpectrum], selected_paths: usize) {
    if spectrum.is_empty() {
        return;
    }

    let total_occurrences: usize = spectrum
        .iter()
        .map(|stat| stat.total_occurrences as usize)
        .sum();
    let max_occurrences = spectrum
        .iter()
        .map(|stat| stat.total_occurrences)
        .max()
        .unwrap_or(0);
    let max_carrier_paths = spectrum
        .iter()
        .map(|stat| stat.carrier_paths)
        .max()
        .unwrap_or(0);
    let max_occurrences_per_path = spectrum
        .iter()
        .map(|stat| stat.max_occurrences_per_path)
        .max()
        .unwrap_or(0);
    let max_path_span_bp = spectrum
        .iter()
        .map(|stat| stat.max_path_span_bp)
        .max()
        .unwrap_or(0);
    info!(
        "[syng2gfa] local syng shared-node spectrum: nodes={}, occurrences={}, selected_paths={}, max_occurrences={}, max_carrier_paths={}, max_occ_per_path={}, max_path_span_bp={}, occ/path_bins(nodes/occ)=[{}], max-copy_bins(nodes/occ)=[{}], dispersion_bins(nodes/occ)=[{}]",
        spectrum.len(),
        total_occurrences,
        selected_paths,
        max_occurrences,
        max_carrier_paths,
        max_occurrences_per_path,
        max_path_span_bp,
        format_local_spectrum_ratio_bins(spectrum),
        format_local_spectrum_max_copy_bins(spectrum),
        format_local_spectrum_dispersion_bins(spectrum)
    );
    info!(
        "[syng2gfa] local syng scaffold-glue threshold candidates (nodes/occurrences): selected-path-factor>1.25={}, dispersed-ratio>=1.25={}, dispersed-ratio>=2.0={}, dispersed-ratio>=3.0={}, dispersed-ratio>=4.0={} (min_occurrences={}, min_path_span_bp={})",
        format_candidate_summary(fixed_selected_path_factor_summary(
            spectrum,
            selected_paths,
            1.25,
        )),
        format_candidate_summary(dispersed_scaffold_glue_summary(spectrum, 1.25)),
        format_candidate_summary(dispersed_scaffold_glue_summary(spectrum, 2.0)),
        format_candidate_summary(dispersed_scaffold_glue_summary(spectrum, 3.0)),
        format_candidate_summary(dispersed_scaffold_glue_summary(spectrum, 4.0)),
        DEFAULT_GFA_SCAFFOLD_GLUE_MIN_OCCURRENCES,
        DEFAULT_GFA_SCAFFOLD_GLUE_MIN_DISPERSION_BP,
    );
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

fn shared_walk_occurrences(
    walks: &[Vec<LocalSyncmerWalkStep>],
    shared_nodes: &FxHashSet<u32>,
) -> Vec<Vec<SharedWalkOccurrence>> {
    walks
        .iter()
        .map(|walk| {
            walk.iter()
                .enumerate()
                .filter_map(|(occ, step)| {
                    shared_nodes
                        .contains(&step.node)
                        .then_some(SharedWalkOccurrence {
                            occ: occ as u32,
                            node: step.node,
                            pos: step.pos,
                        })
                })
                .collect()
        })
        .collect()
}

fn scaffold_run_seeds_by_key(
    shared_walks: &[Vec<SharedWalkOccurrence>],
    min_shared_run: usize,
    max_anchor_step_bp: u64,
) -> FxHashMap<Vec<u32>, Vec<ScaffoldRunSeed>> {
    let mut seeds_by_key: FxHashMap<Vec<u32>, Vec<ScaffoldRunSeed>> = FxHashMap::default();
    if min_shared_run == 0 {
        return seeds_by_key;
    }

    for (path_idx, walk) in shared_walks.iter().enumerate() {
        if walk.len() < min_shared_run {
            continue;
        }
        for start in 0..=walk.len() - min_shared_run {
            let window = &walk[start..start + min_shared_run];
            if window
                .windows(2)
                .any(|pair| pair[0].pos.abs_diff(pair[1].pos) > max_anchor_step_bp)
            {
                continue;
            }
            let key = window.iter().map(|occ| occ.node).collect::<Vec<_>>();
            seeds_by_key
                .entry(key)
                .or_default()
                .push(ScaffoldRunSeed { path_idx, start });
        }
    }

    seeds_by_key
}

fn seed_pair_count(seeds_by_path: &[(usize, Vec<usize>)]) -> usize {
    let mut total = 0usize;
    let mut previous = 0usize;
    for (_path_idx, starts) in seeds_by_path {
        total = total.saturating_add(previous.saturating_mul(starts.len()));
        previous = previous.saturating_add(starts.len());
    }
    total
}

fn support_seed_occurrences(
    supported: &mut FxHashSet<LocalOccurrenceKey>,
    shared_walks: &[Vec<SharedWalkOccurrence>],
    seeds_by_path: &[(usize, Vec<usize>)],
    min_shared_run: usize,
) {
    for (path_idx, starts) in seeds_by_path {
        let Some(walk) = shared_walks.get(*path_idx) else {
            continue;
        };
        for &start in starts {
            for offset in 0..min_shared_run {
                if let Some(occ) = walk.get(start + offset) {
                    supported.insert(LocalOccurrenceKey {
                        path_idx: *path_idx,
                        occ: occ.occ,
                    });
                }
            }
        }
    }
}

fn scaffold_supported_occurrences(
    walks: &[Vec<LocalSyncmerWalkStep>],
    syncmer_len: u64,
    min_shared_run: usize,
    shared_nodes: &FxHashSet<u32>,
) -> (FxHashSet<LocalOccurrenceKey>, usize, usize) {
    if min_shared_run <= 1 || walks.len() < 2 || shared_nodes.is_empty() {
        return (FxHashSet::default(), 0, 0);
    }

    let shared_walks = shared_walk_occurrences(walks, shared_nodes);
    let max_anchor_step_bp = SCAFFOLD_CONTEXT_GAP_BP.saturating_add(syncmer_len.max(1));
    let seeds_by_key = scaffold_run_seeds_by_key(&shared_walks, min_shared_run, max_anchor_step_bp);
    if seeds_by_key.is_empty() {
        return (FxHashSet::default(), 0, 0);
    }

    let mut direct_supported = FxHashSet::default();
    let mut pair_anchor_keys: FxHashMap<(usize, usize), FxHashSet<PairAnchorOccurrenceKey>> =
        FxHashMap::default();
    let mut dense_signatures = 0usize;
    let mut projected_candidate_anchors = 0usize;

    for (_run_key, seeds) in seeds_by_key {
        let mut grouped: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
        for seed in seeds {
            grouped.entry(seed.path_idx).or_default().push(seed.start);
        }
        if grouped.len() < 2 {
            continue;
        }

        let mut seeds_by_path = grouped.into_iter().collect::<Vec<_>>();
        seeds_by_path.sort_unstable_by_key(|(path_idx, _starts)| *path_idx);
        support_seed_occurrences(
            &mut direct_supported,
            &shared_walks,
            &seeds_by_path,
            min_shared_run,
        );

        let seed_pairs = seed_pair_count(&seeds_by_path);
        let signature_candidate_anchors = seed_pairs.saturating_mul(min_shared_run);
        let exceeds_signature_budget = seed_pairs > MAX_SCAFFOLD_SIGNATURE_SEED_PAIRS;
        let exceeds_total_budget = projected_candidate_anchors
            .saturating_add(signature_candidate_anchors)
            > MAX_SCAFFOLD_TOTAL_CANDIDATE_ANCHORS;
        if exceeds_signature_budget || exceeds_total_budget {
            dense_signatures += 1;
            continue;
        }
        projected_candidate_anchors =
            projected_candidate_anchors.saturating_add(signature_candidate_anchors);

        for left_idx in 0..seeds_by_path.len() {
            for right_idx in (left_idx + 1)..seeds_by_path.len() {
                let (query_idx, query_starts) = &seeds_by_path[left_idx];
                let (target_idx, target_starts) = &seeds_by_path[right_idx];
                let query_walk = &shared_walks[*query_idx];
                let target_walk = &shared_walks[*target_idx];
                let pair_set = pair_anchor_keys
                    .entry((*query_idx, *target_idx))
                    .or_default();
                for &query_start in query_starts {
                    for &target_start in target_starts {
                        for offset in 0..min_shared_run {
                            let query = query_walk[query_start + offset];
                            let target = target_walk[target_start + offset];
                            debug_assert_eq!(query.node, target.node);
                            pair_set.insert(PairAnchorOccurrenceKey {
                                query_occ: query.occ,
                                target_occ: target.occ,
                                node: query.node,
                                query_pos: query.pos,
                                target_pos: target.pos,
                            });
                        }
                    }
                }
            }
        }
    }

    let min_scaffold_length = (min_shared_run as u64).saturating_mul(syncmer_len.max(1));
    let mut supported = direct_supported;
    let mut candidate_anchors = 0usize;
    for ((query_idx, target_idx), anchor_keys) in pair_anchor_keys {
        candidate_anchors = candidate_anchors.saturating_add(anchor_keys.len());
        let anchors = anchor_keys
            .into_iter()
            .map(|key| PairAnchorOccurrence {
                query_occ: key.query_occ,
                target_occ: key.target_occ,
                node: key.node,
                query_pos: key.query_pos,
                target_pos: key.target_pos,
            })
            .collect::<Vec<_>>();
        supported.extend(retained_pair_scaffold_occurrences(
            query_idx,
            target_idx,
            &anchors,
            syncmer_len,
            min_scaffold_length,
        ));
    }
    (supported, candidate_anchors, dense_signatures)
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
    high_freq_min_run: usize,
    high_freq_min_sequence_span_bp: usize,
) -> SharedContextFilterResult {
    let generic_context_enabled = min_shared_run > 1 || min_sequence_span_bp > 0;
    let high_frequency_context_enabled =
        high_freq_min_run > 0 || high_freq_min_sequence_span_bp > 0;
    if walks.is_empty() || (!generic_context_enabled && !high_frequency_context_enabled) {
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

    let spectrum = local_syncmer_node_spectrum(walks, &shared_nodes);
    log_local_syncmer_spectrum(&spectrum, walks.len());
    let scaffold_glue_nodes = spectrum_scaffold_glue_nodes(&spectrum);
    let scaffold_nodes: FxHashSet<u32> = shared_nodes
        .difference(&scaffold_glue_nodes)
        .copied()
        .collect();
    let shared_occurrences = shared_occurrence_keys(walks, &shared_nodes);
    let path_copy_filtered_occurrences = shared_occurrence_keys(walks, &scaffold_glue_nodes).len();
    if !scaffold_glue_nodes.is_empty() {
        info!(
            "[syng2gfa] spectrum-selected {} dispersed high-copy scaffold-glue node(s) covering {} local occurrence(s) for unified run/span masking (min_occ_per_path_ratio={}, min_occurrences={}, min_path_span_bp={})",
            scaffold_glue_nodes.len(),
            path_copy_filtered_occurrences,
            DEFAULT_GFA_SCAFFOLD_GLUE_MIN_OCC_PER_PATH_RATIO,
            DEFAULT_GFA_SCAFFOLD_GLUE_MIN_OCCURRENCES,
            DEFAULT_GFA_SCAFFOLD_GLUE_MIN_DISPERSION_BP
        );
    }
    let (scaffold_supported, scaffold_candidate_anchors, scaffold_dense_signatures) =
        scaffold_supported_occurrences(walks, syncmer_len, min_shared_run, &scaffold_nodes);
    let sequence_supported =
        sequence_supported_occurrences(walks, syncmer_len, min_sequence_span_bp, &scaffold_nodes);
    let spectrum_glue_policy = high_frequency_occurrence_policy(
        walks,
        syncmer_len,
        &scaffold_glue_nodes,
        high_freq_min_run,
        high_freq_min_sequence_span_bp,
    );
    let supported: FxHashSet<LocalOccurrenceKey> = scaffold_supported
        .union(&sequence_supported)
        .copied()
        .chain(
            spectrum_glue_policy
                .supported_occurrence_keys
                .iter()
                .copied(),
        )
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
        scaffold_candidate_anchors,
        scaffold_dense_signatures,
        path_copy_filtered_nodes: scaffold_glue_nodes.len(),
        path_copy_filtered_occurrences,
        path_copy_supported_occurrences: spectrum_glue_policy.supported_occurrences,
        path_copy_run_supported_occurrences: spectrum_glue_policy.run_supported_occurrences,
        path_copy_sequence_supported_occurrences: spectrum_glue_policy
            .sequence_supported_occurrences,
        path_copy_private_split_occurrences: spectrum_glue_policy.private_split_occurrences,
        path_copy_scaffold_candidate_anchors: spectrum_glue_policy.scaffold_candidate_anchors,
        path_copy_scaffold_dense_signatures: spectrum_glue_policy.scaffold_dense_signatures,
    }
}

fn high_frequency_occurrence_policy(
    walks: &[Vec<LocalSyncmerWalkStep>],
    syncmer_len: u64,
    high_frequency_nodes: &FxHashSet<u32>,
    min_run: usize,
    min_sequence_span_bp: usize,
) -> HighFrequencyOccurrencePolicyResult {
    let mut result = HighFrequencyOccurrencePolicyResult {
        split_occurrences: vec![FxHashSet::default(); walks.len()],
        selected_nodes: high_frequency_nodes.len(),
        ..HighFrequencyOccurrencePolicyResult::default()
    };
    if walks.is_empty() || high_frequency_nodes.is_empty() {
        return result;
    }

    let high_frequency_occurrences = shared_occurrence_keys(walks, high_frequency_nodes);
    result.total_occurrences = high_frequency_occurrences.len();
    if high_frequency_occurrences.is_empty() {
        return result;
    }

    let (run_supported, scaffold_candidate_anchors, scaffold_dense_signatures) = match min_run {
        0 => (FxHashSet::default(), 0, 0),
        1 => (high_frequency_occurrences.clone(), 0, 0),
        _ => scaffold_supported_occurrences(walks, syncmer_len, min_run, high_frequency_nodes),
    };
    let sequence_supported = sequence_supported_occurrences(
        walks,
        syncmer_len,
        min_sequence_span_bp,
        high_frequency_nodes,
    );
    let supported: FxHashSet<LocalOccurrenceKey> =
        run_supported.union(&sequence_supported).copied().collect();

    result.run_supported_occurrences = run_supported.len();
    result.sequence_supported_occurrences = sequence_supported.len();
    result.scaffold_candidate_anchors = scaffold_candidate_anchors;
    result.scaffold_dense_signatures = scaffold_dense_signatures;
    result.supported_occurrences = high_frequency_occurrences.intersection(&supported).count();
    result.supported_occurrence_keys = high_frequency_occurrences
        .intersection(&supported)
        .copied()
        .collect();

    for key in high_frequency_occurrences.difference(&supported).copied() {
        if let Some(path_set) = result.split_occurrences.get_mut(key.path_idx) {
            path_set.insert(key.occ);
            result.private_split_occurrences += 1;
        }
    }
    result
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

fn local_walks_without_nodes(
    walks: &[Vec<LocalSyncmerWalkStep>],
    excluded_nodes: &FxHashSet<u32>,
) -> FilteredLocalWalks {
    if excluded_nodes.is_empty() {
        return FilteredLocalWalks {
            walks: walks.to_vec(),
            original_occurrences: walks
                .iter()
                .map(|walk| (0..walk.len() as u32).collect())
                .collect(),
        };
    }

    let mut filtered_walks = Vec::with_capacity(walks.len());
    let mut original_occurrences = Vec::with_capacity(walks.len());
    for walk in walks {
        let mut filtered = Vec::new();
        let mut mapping = Vec::new();
        for (original_occ, step) in walk.iter().copied().enumerate() {
            if excluded_nodes.contains(&step.node) {
                continue;
            }
            filtered.push(step);
            mapping.push(original_occ as u32);
        }
        filtered_walks.push(filtered);
        original_occurrences.push(mapping);
    }
    FilteredLocalWalks {
        walks: filtered_walks,
        original_occurrences,
    }
}

fn remap_occurrences_to_original(
    occurrences: &[FxHashSet<u32>],
    original_occurrences: &[Vec<u32>],
) -> Vec<FxHashSet<u32>> {
    occurrences
        .iter()
        .enumerate()
        .map(|(path_idx, path_occurrences)| {
            let Some(mapping) = original_occurrences.get(path_idx) else {
                return FxHashSet::default();
            };
            path_occurrences
                .iter()
                .filter_map(|occ| mapping.get(*occ as usize).copied())
                .collect()
        })
        .collect()
}

fn merge_split_occurrences(destination: &mut [FxHashSet<u32>], source: &[FxHashSet<u32>]) -> usize {
    let mut inserted = 0usize;
    for (dest, src) in destination.iter_mut().zip(source) {
        for occ in src {
            if dest.insert(*occ) {
                inserted += 1;
            }
        }
    }
    inserted
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
            "[syng2gfa] selected {} / {} high-frequency local syncmer node(s) for frequency masking (top_fraction={}, max_occurrences={:?}, local max={}, selected min/max={}/{}) in {:.3}s",
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

fn evaluate_range_masking(
    index: &SyngIndex,
    local_walks: &[Vec<LocalSyncmerWalkStep>],
    syncmer_len: u64,
    frequency_mask: SyngGfaFrequencyMask,
    path_count: usize,
    materialization_label: &str,
) -> io::Result<RangeMaskEvaluation> {
    let mut evaluation = RangeMaskEvaluation::empty(path_count);
    if local_walks.is_empty() {
        return Ok(evaluation);
    }

    let all_local_syncmers = local_walks
        .iter()
        .flatten()
        .map(|step| step.node)
        .collect::<FxHashSet<_>>();
    let selected_frequency_syncmers = if frequency_mask.frequency_filter_enabled() {
        frequency_masked_syncmers(index, &all_local_syncmers, frequency_mask)?
    } else {
        FxHashSet::default()
    };
    let run_aware_frequency = frequency_mask.run_aware_frequency_filter_enabled();

    if run_aware_frequency {
        let high_freq_start = Instant::now();
        let high_frequency = high_frequency_occurrence_policy(
            local_walks,
            syncmer_len,
            &selected_frequency_syncmers,
            frequency_mask.high_freq_min_run,
            frequency_mask.high_freq_min_sequence_span_bp,
        );
        evaluation.high_frequency_split = high_frequency.private_split_occurrences;
        merge_split_occurrences(
            &mut evaluation.split_occurrences,
            &high_frequency.split_occurrences,
        );
        info!(
            "[syng2gfa] high-frequency occurrence mask selected {} node(s) covering {} local occurrence(s); rescued/supported {} occurrence(s) (run-supported={}, sequence-supported={}, scaffold-candidates={}, dense-signatures={}), privately split {} occurrence(s), settings: freq-run={}, freq-span={}bp before {} in {:.3}s",
            high_frequency.selected_nodes,
            high_frequency.total_occurrences,
            high_frequency.supported_occurrences,
            high_frequency.run_supported_occurrences,
            high_frequency.sequence_supported_occurrences,
            high_frequency.scaffold_candidate_anchors,
            high_frequency.scaffold_dense_signatures,
            high_frequency.private_split_occurrences,
            frequency_mask.high_freq_min_run,
            frequency_mask.high_freq_min_sequence_span_bp,
            materialization_label,
            high_freq_start.elapsed().as_secs_f64()
        );
    } else {
        evaluation.filtered_syncmers = selected_frequency_syncmers.clone();
        if frequency_mask.frequency_filter_enabled() {
            info!(
                "[syng2gfa] high-frequency occurrence mask disabled; using legacy node-level removal for {} selected high-frequency syncmer node(s) before {}",
                evaluation.filtered_syncmers.len(),
                materialization_label
            );
        }
    }

    if frequency_mask.shared_context_filter_enabled() {
        let context_start = Instant::now();
        let context_filter = if run_aware_frequency && !selected_frequency_syncmers.is_empty() {
            let filtered = local_walks_without_nodes(local_walks, &selected_frequency_syncmers);
            let mut context_filter = shared_context_filtered_occurrences(
                &filtered.walks,
                syncmer_len,
                frequency_mask.min_shared_run,
                frequency_mask.min_sequence_span_bp,
                frequency_mask.high_freq_min_run,
                frequency_mask.high_freq_min_sequence_span_bp,
            );
            context_filter.weak_occurrences = remap_occurrences_to_original(
                &context_filter.weak_occurrences,
                &filtered.original_occurrences,
            );
            context_filter
        } else {
            let context_walks =
                local_walks_after_frequency_mask(local_walks, &evaluation.filtered_syncmers);
            shared_context_filtered_occurrences(
                &context_walks,
                syncmer_len,
                frequency_mask.min_shared_run,
                frequency_mask.min_sequence_span_bp,
                frequency_mask.high_freq_min_run,
                frequency_mask.high_freq_min_sequence_span_bp,
            )
        };
        evaluation.shared_context_split = context_filter
            .weak_occurrences
            .iter()
            .map(FxHashSet::len)
            .sum::<usize>();
        let inserted = merge_split_occurrences(
            &mut evaluation.split_occurrences,
            &context_filter.weak_occurrences,
        );
        info!(
            "[syng2gfa] scaffold-context evaluated {} local shared syncmer occurrence(s) before {} (weak={}, split={}, min_run={}, spectrum-glue-filtered-nodes={}, spectrum-glue-filtered-occurrences={}, scaffold-supported={}, scaffold-candidates={}, dense-signatures={}, sequence_k={}, sequence-supported={}, spectrum-glue-rescued={}, spectrum-glue-run-supported={}, spectrum-glue-sequence-supported={}, spectrum-glue-private-split={}, spectrum-glue-scaffold-candidates={}, spectrum-glue-dense-signatures={}, freq-run={}, freq-span={}) in {:.3}s",
            context_filter.shared_occurrences,
            materialization_label,
            evaluation.shared_context_split,
            inserted,
            frequency_mask.min_shared_run,
            context_filter.path_copy_filtered_nodes,
            context_filter.path_copy_filtered_occurrences,
            context_filter.scaffold_supported_occurrences,
            context_filter.scaffold_candidate_anchors,
            context_filter.scaffold_dense_signatures,
            frequency_mask.min_sequence_span_bp,
            context_filter.sequence_supported_occurrences,
            context_filter.path_copy_supported_occurrences,
            context_filter.path_copy_run_supported_occurrences,
            context_filter.path_copy_sequence_supported_occurrences,
            context_filter.path_copy_private_split_occurrences,
            context_filter.path_copy_scaffold_candidate_anchors,
            context_filter.path_copy_scaffold_dense_signatures,
            frequency_mask.high_freq_min_run,
            frequency_mask.high_freq_min_sequence_span_bp,
            context_start.elapsed().as_secs_f64()
        );
        if inserted > 0 {
            info!(
                "[syng2gfa] scaffold-context split {} local syncmer occurrence(s) into private topology before {}",
                inserted,
                materialization_label
            );
        }
    }

    Ok(evaluation)
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
    let mask_policy_split = split_occurrences.iter().map(FxHashSet::len).sum::<usize>();

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
        "[syng2gfa] reduced exact blunt path walks into {} full untrimmed syncmer anchor segment(s), {} interned exact sequence segment(s) ({} gap, {} trimmed syncmer-slice; {} local-repeat clones, {} mask-policy clones from {} split occurrence(s)), {} edge(s) in {:.3}s",
        used_syncmers.len(),
        interner.total_extra_segments(),
        interner.gap_segments,
        interner.total_extra_segments().saturating_sub(interner.gap_segments),
        n_local_repeat_clones,
        n_scaffold_context_clones,
        mask_policy_split,
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
    let mask_evaluation = if let Some(walks) = local_walks.as_ref() {
        evaluate_range_masking(
            index,
            walks,
            syncmer_len_u64,
            frequency_mask,
            ranges.len(),
            "exact blunt materialization",
        )?
    } else {
        RangeMaskEvaluation::empty(ranges.len())
    };
    let masked_syncmers = mask_evaluation.filtered_syncmers;
    let split_occurrences = mask_evaluation.split_occurrences;

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
    let mask_evaluation = if let Some(walks) = local_walks.as_ref() {
        evaluate_range_masking(
            index,
            walks,
            syncmer_len_u64,
            frequency_mask,
            ranges.len(),
            "raw GFA materialization",
        )?
    } else {
        RangeMaskEvaluation::empty(ranges.len())
    };
    let high_frequency_split = mask_evaluation.high_frequency_split;
    let shared_context_split = mask_evaluation.shared_context_split;
    let masked_syncmers = mask_evaluation.filtered_syncmers;
    let split_occurrences = mask_evaluation.split_occurrences;
    let mask_policy_split = split_occurrences.iter().map(FxHashSet::len).sum::<usize>();

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
        "[syng2gfa] reduced selected path walks into {} shared syncmer node(s), {} cloned syncmer segment(s) ({} local-repeat, {} mask-policy), {} legacy frequency syncmer node(s) removed, {} high-frequency occurrence(s) split, and {} scaffold-context occurrence(s) split before raw GFA ({} total mask-policy split occurrence(s)), {} unique gap segment(s), {} edge(s) in {:.3}s",
        used_syncmers.len(),
        cloned_syncmers.len(),
        n_local_repeat_clones,
        n_scaffold_context_clones,
        masked_syncmers.len(),
        high_frequency_split,
        shared_context_split,
        mask_policy_split,
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
    fn test_range_gfa_legacy_frequency_filter_drops_repetitive_syncmers_before_raw_gfa() {
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
            high_freq_min_run: 0,
            high_freq_min_sequence_span_bp: 0,
            run_aware_frequency_mask: false,
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
            "legacy high-frequency syncmers should be filtered, not cloned into raw GFA:\n{gfa}"
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
    fn test_high_frequency_occurrence_policy_splits_isolated_reuse() {
        let walks = vec![local_walk(&[(99, 0)]), local_walk(&[(99, 10_000)])];
        let high_frequency_nodes = [99].into_iter().collect::<FxHashSet<_>>();

        let filtered = high_frequency_occurrence_policy(&walks, 31, &high_frequency_nodes, 3, 0);

        assert_eq!(filtered.selected_nodes, 1);
        assert_eq!(filtered.total_occurrences, 2);
        assert_eq!(filtered.supported_occurrences, 0);
        assert_eq!(filtered.private_split_occurrences, 2);
        assert!(filtered.split_occurrences[0].contains(&0));
        assert!(filtered.split_occurrences[1].contains(&0));
    }

    #[test]
    fn test_high_frequency_occurrence_policy_keeps_configured_run_shared() {
        let run = [(10, 0), (11, 100), (12, 200), (13, 300)];
        let walks = vec![local_walk(&run), local_walk(&run)];
        let high_frequency_nodes = [10, 11, 12, 13].into_iter().collect::<FxHashSet<_>>();

        let filtered = high_frequency_occurrence_policy(&walks, 31, &high_frequency_nodes, 4, 0);

        assert_eq!(filtered.total_occurrences, 8);
        assert_eq!(filtered.supported_occurrences, 8);
        assert_eq!(filtered.run_supported_occurrences, 8);
        assert_eq!(filtered.private_split_occurrences, 0);
        assert!(filtered.split_occurrences.iter().all(FxHashSet::is_empty));
    }

    #[test]
    fn test_high_frequency_occurrence_policy_keeps_configured_exact_span_shared() {
        let run = [(20, 0), (21, 32), (22, 64)];
        let walks = vec![local_walk(&run), local_walk(&run)];
        let high_frequency_nodes = [20, 21, 22].into_iter().collect::<FxHashSet<_>>();

        let filtered = high_frequency_occurrence_policy(&walks, 63, &high_frequency_nodes, 0, 127);

        assert_eq!(filtered.total_occurrences, 6);
        assert_eq!(filtered.supported_occurrences, 6);
        assert_eq!(filtered.sequence_supported_occurrences, 6);
        assert_eq!(filtered.private_split_occurrences, 0);
        assert!(filtered.split_occurrences.iter().all(FxHashSet::is_empty));
    }

    #[test]
    fn test_scaffold_context_retains_chain_members_and_splits_off_chain_occurrences() {
        let walks = vec![
            local_walk(&[(40, 0), (41, 100), (42, 200), (42, 5_000)]),
            local_walk(&[(40, 0), (41, 100), (42, 200)]),
            local_walk(&[(42, 7_000)]),
        ];

        let filtered = shared_context_filtered_occurrences(&walks, 31, 3, 0, 3, 0);
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
    fn test_scaffold_context_prunes_high_copy_singletons_without_global_blessing() {
        let mut walks = vec![
            local_walk(&[(10, 0), (11, 100), (99, 200), (12, 300), (13, 400)]),
            local_walk(&[(10, 0), (11, 100), (99, 200), (12, 300), (13, 400)]),
        ];
        for copy in 0..64 {
            walks.push(local_walk(&[(99, 10_000 + copy * 2_000)]));
        }

        let filtered = shared_context_filtered_occurrences(&walks, 31, 5, 0, 5, 0);
        assert_eq!(filtered.shared_occurrences, 74);
        assert_eq!(
            filtered.scaffold_supported_occurrences, 10,
            "only the two five-anchor run copies should be scaffold-supported"
        );
        assert!(
            filtered.scaffold_candidate_anchors <= 5,
            "high-copy singleton copies must not materialize all pairwise anchors"
        );

        for occ in 0..5 {
            assert!(
                !filtered.weak_occurrences[0].contains(&occ)
                    && !filtered.weak_occurrences[1].contains(&occ),
                "occurrence {occ} belongs to the retained scaffold run"
            );
        }
        for path_idx in 2..walks.len() {
            assert!(
                filtered.weak_occurrences[path_idx].contains(&0),
                "isolated high-copy path {path_idx} should remain private"
            );
        }
    }

    #[test]
    fn test_scaffold_context_splits_dispersed_high_copy_glue_from_spectrum() {
        let repeated = (0..16).map(|copy| (99, copy * 2_000)).collect::<Vec<_>>();
        let walks = (0..4).map(|_| local_walk(&repeated)).collect::<Vec<_>>();

        let filtered = shared_context_filtered_occurrences(&walks, 31, 3, 0, 3, 0);
        assert_eq!(filtered.shared_occurrences, 64);
        assert_eq!(filtered.path_copy_filtered_nodes, 1);
        assert_eq!(
            filtered.path_copy_filtered_occurrences, 64,
            "node 99 has a high occ/path ratio and scaffold-scale positional dispersion"
        );
        assert_eq!(filtered.scaffold_supported_occurrences, 0);
        assert_eq!(filtered.path_copy_supported_occurrences, 0);
        assert_eq!(filtered.path_copy_private_split_occurrences, 64);
        assert_eq!(
            filtered
                .weak_occurrences
                .iter()
                .map(FxHashSet::len)
                .sum::<usize>(),
            64
        );
    }

    #[test]
    fn test_scaffold_context_rescues_spectrum_glue_with_high_frequency_run_policy() {
        let repeated = (0..16).map(|copy| (99, copy * 100)).collect::<Vec<_>>();
        let walks = (0..4).map(|_| local_walk(&repeated)).collect::<Vec<_>>();

        let filtered = shared_context_filtered_occurrences(&walks, 31, 3, 0, 10, 0);

        assert_eq!(filtered.shared_occurrences, 64);
        assert_eq!(filtered.path_copy_filtered_nodes, 1);
        assert_eq!(filtered.path_copy_filtered_occurrences, 64);
        assert_eq!(
            filtered.path_copy_run_supported_occurrences, 64,
            "spectrum-selected nodes must use the same freq-run rescue as explicit high-frequency nodes"
        );
        assert_eq!(filtered.path_copy_supported_occurrences, 64);
        assert_eq!(filtered.path_copy_private_split_occurrences, 0);
        assert!(
            filtered.weak_occurrences.iter().all(FxHashSet::is_empty),
            "a spectrum-selected node in a supported run should remain shared"
        );
    }

    #[test]
    fn test_scaffold_context_rescues_spectrum_glue_with_high_frequency_span_policy() {
        let repeated = (0..20).map(|copy| (99, copy * 63)).collect::<Vec<_>>();
        let walks = (0..4).map(|_| local_walk(&repeated)).collect::<Vec<_>>();

        let filtered = shared_context_filtered_occurrences(&walks, 63, 1, 0, 0, 1_000);

        assert_eq!(filtered.shared_occurrences, 80);
        assert_eq!(filtered.path_copy_filtered_nodes, 1);
        assert_eq!(filtered.path_copy_filtered_occurrences, 80);
        assert_eq!(
            filtered.path_copy_sequence_supported_occurrences, 80,
            "spectrum-selected nodes must use the same freq-span rescue as explicit high-frequency nodes"
        );
        assert_eq!(filtered.path_copy_supported_occurrences, 80);
        assert_eq!(filtered.path_copy_private_split_occurrences, 0);
        assert!(
            filtered.weak_occurrences.iter().all(FxHashSet::is_empty),
            "a spectrum-selected node in a supported exact span should remain shared"
        );
    }

    #[test]
    fn test_scaffold_context_preserves_compact_repeated_copy_structure() {
        let walks = vec![
            local_walk(&[(99, 0), (99, 100), (99, 200)]),
            local_walk(&[(99, 0), (99, 100), (99, 200)]),
            local_walk(&[(99, 0), (99, 100), (99, 200)]),
            local_walk(&[(99, 0), (99, 100), (99, 200)]),
        ];

        let filtered = shared_context_filtered_occurrences(&walks, 31, 3, 0, 3, 0);
        assert_eq!(filtered.shared_occurrences, 12);
        assert_eq!(
            filtered.path_copy_filtered_occurrences, 0,
            "compact local repeat copies are not selected by the dispersed-glue spectrum rule"
        );
        assert_eq!(
            filtered.scaffold_supported_occurrences, 12,
            "bounded occurrence-level scaffold support can preserve compact repeated structure"
        );
        assert!(
            filtered.weak_occurrences.iter().all(FxHashSet::is_empty),
            "compact repeated structure should stay shared for later graph resolution"
        );
    }

    #[test]
    fn test_scaffold_context_dense_repeated_runs_stay_occurrence_level() {
        let repeated = [(1, 0), (2, 100), (3, 200), (4, 300), (5, 400)];
        let walks = (0..10).map(|_| local_walk(&repeated)).collect::<Vec<_>>();

        let filtered = shared_context_filtered_occurrences(&walks, 31, 5, 0, 5, 0);
        assert!(
            filtered.scaffold_dense_signatures > 0,
            "test-sized dense signatures should exercise the bounded fallback"
        );
        assert_eq!(filtered.path_copy_filtered_occurrences, 0);
        assert_eq!(filtered.scaffold_supported_occurrences, 50);
        assert!(
            filtered.weak_occurrences.iter().all(FxHashSet::is_empty),
            "dense exact runs are still supported by occurrence, not by global node id"
        );
    }

    #[test]
    fn test_scaffold_context_keeps_indel_shifted_colinear_runs() {
        let walks = vec![
            local_walk(&[(1, 0), (2, 100), (3, 220), (4, 320)]),
            local_walk(&[(1, 0), (2, 100), (3, 260), (4, 360)]),
        ];

        let filtered = shared_context_filtered_occurrences(&walks, 31, 4, 0, 4, 0);
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

        let sequence_supported = shared_context_filtered_occurrences(&walks, 63, 1, 127, 1, 127);
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

        let too_short = shared_context_filtered_occurrences(&walks, 63, 1, 191, 1, 191);
        for occ in [3, 4, 5] {
            assert!(
                too_short.weak_occurrences[0].contains(&occ)
                    && too_short.weak_occurrences[1].contains(&occ),
                "127 bp evidence should not satisfy sequence-k=191"
            );
        }

        let rescued_by_scaffold = shared_context_filtered_occurrences(&walks, 63, 3, 191, 3, 191);
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
