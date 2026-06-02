//! Bubble-guided graph resolution primitives.
//!
//! This module implements the first conservative slice of hierarchical graph
//! resolution: detect path-supported bubbles in a blunt GFA, replace bounded
//! single-entry/single-exit bubbles with exact path-preserving local graph
//! induction, and repeat until no eligible unseen bubbles remain. Pairwise
//! graph-induction resolvers are followed by a bounded bubble-aware POA polish.
//! The implementation intentionally avoids lossy representative collapse and
//! coordinate sidecars; emitted paths are the coordinate system.

use crate::graph::{
    build_global_spoa_engine, feed_sequences_to_graph, reverse_complement, unchop_gfa,
};
use povu::native_gfa::{Step as PovuStep, Strand as PovuStrand};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::{BTreeSet, HashSet};
use std::hash::{Hash, Hasher};
use std::io::{self, Read};
use std::panic::{self, AssertUnwindSafe};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

const GRAPH_QUALITY_LONG_BRIDGE_BP: usize = 10_000;
const GRAPH_QUALITY_WHITE_SPACE_FLOOR_BP: usize = 1_000;
static DEBUG_REPLACEMENT_ID: AtomicUsize = AtomicUsize::new(0);
static CHAIN_POVU_SMOOTH_USED: AtomicUsize = AtomicUsize::new(0);
static CHAIN_POVU_DIRECT_ON_EMPTY: AtomicUsize = AtomicUsize::new(0);
static CHAIN_POVU_DIRECT_ON_SMOOTH_FAILURE: AtomicUsize = AtomicUsize::new(0);

/// Configuration for exact path-preserving bubble resolution.
#[derive(Clone, Debug)]
pub struct ResolutionConfig {
    /// Maximum number of frontier replacement rounds.
    pub max_iterations: usize,
    /// Local graph induction method used for selected bubbles.
    pub method: ResolutionMethod,
    /// Upper bound (exclusive) on **median traversal length** for `method=auto` to
    /// route a bubble to direct SPOA. Per docs/crush-architecture-spec.md §Phase-2:
    /// `median < auto_spoa_max_traversal_len` → sPOA. Set to 0 to disable the
    /// sPOA tier entirely (everything routes to POASTA or sweepga).
    pub auto_spoa_max_traversal_len: usize,
    /// Upper bound (exclusive) on **median traversal length** for `method=auto` to
    /// route a bubble to POASTA. Bubbles with
    /// `auto_spoa_max_traversal_len <= median < auto_poasta_max_traversal_len`
    /// route to POASTA; bubbles at or above this threshold route to sweepga.
    /// Set to 0 to disable the POASTA tier (everything above the sPOA tier goes
    /// to sweepga).
    pub auto_poasta_max_traversal_len: usize,
    /// Legacy auto-routing limit retained for CLI compatibility. Unused by the
    /// median-based 3-tier dispatch but still accepted on the CLI so existing
    /// command lines continue to parse.
    pub auto_allwave_max_total_sequence: usize,
    /// Legacy auto-routing limit retained for CLI compatibility. See
    /// `auto_allwave_max_total_sequence`.
    pub auto_allwave_max_traversals: usize,
    /// Maximum root path span, in bp, for a candidate bubble.
    ///
    /// The root path is the first path in the input GFA. This is a rooted POVU
    /// decomposition coordinate, not necessarily an external reference genome.
    pub max_bubble_span: usize,
    /// Maximum length of any observed traversal through the bubble for direct
    /// replacement methods (`poa`, `poasta`, `star-biwfa`).
    pub max_traversal_len: usize,
    /// Minimum length of the longest observed traversal through the bubble.
    ///
    /// This lets callers run hierarchical passes explicitly, e.g. first crush
    /// parent-scale bubbles with a pairwise graph-induction engine, then run a
    /// separate SPOA pass over remaining small local bubbles.
    pub min_traversal_len: usize,
    /// Maximum median traversal length for direct replacement methods. A value
    /// of 0 disables this guard.
    pub max_median_traversal_len: usize,
    /// Maximum sum of traversal sequence lengths for one direct replacement.
    pub max_total_sequence: usize,
    /// Maximum number of path-supported traversals through one direct replacement.
    pub max_traversals: usize,
    /// Small-tangle POASTA polish rounds after pairwise graph induction.
    ///
    /// This pass is deliberately scale-limited: small STR / indel tangles are
    /// implementation artifacts, while larger VNTR / recombination structures
    /// are often biological signal and should remain represented unless the
    /// user explicitly raises these budgets.
    pub polish_iterations: usize,
    pub polish_max_traversal_len: usize,
    pub polish_max_median_traversal_len: usize,
    pub polish_max_total_sequence: usize,
    pub polish_max_traversals: usize,
    /// Local cleanup method after pairwise graph induction.
    pub polish_method: ResolutionPolishMethod,
    /// Pair-sampling nearest-neighbor count for many-sequence pairwise engines.
    pub pair_k_nearest: usize,
    /// Pair-sampling farthest-neighbor count for diversity/stranger joins.
    pub pair_k_farthest: usize,
    /// Number of independent tree-sampling passes to union for pairwise engines.
    ///
    /// Additional trees use nearby Mash k-mer sizes to perturb nearest/farthest
    /// neighborhoods deterministically. This is more aggressive than raising
    /// random sampling alone and is aimed at underaligned local bubbles where
    /// one tree misses important allelic bridges.
    pub pair_tree_count: usize,
    /// Deterministic random pair fraction added after tree/kNN sampling.
    pub pair_random_fraction: f64,
    /// Mash k-mer size for pair selection.
    pub pair_mash_k: usize,
    /// Legacy/default exact-match length used by adaptive replacement seqwish policy.
    ///
    /// Local bubble induction defaults to no exact-run floor; this value is
    /// only used when `replacement_min_match_len_policy` is explicitly set to
    /// `Adaptive`.
    pub replacement_seqwish_min_match_len: u64,
    /// Policy for choosing the exact-run floor passed to local seqwish induction.
    ///
    /// `Fixed(0)` and `Fixed(1)` both effectively disable exact-run filtering
    /// because seqwish needs a positive length and every non-empty exact run is
    /// at least 1 bp. `Adaptive` is opt-in: start from
    /// `replacement_seqwish_min_match_len`, then lower it per local block when
    /// the block or its observed CIGAR evidence is shorter.
    pub replacement_min_match_len_policy: ReplacementMinMatchLenPolicy,
    /// Minimum pairwise mapping length kept before replacement seqwish induction.
    ///
    /// A value of 0 follows `replacement_seqwish_min_match_len`, so the PAF
    /// filter and seqwish transitive closure agree on the minimum scale of
    /// evidence allowed to glue sequence together.
    pub replacement_min_map_length: u64,
    /// Minimum pairwise mapping identity kept before replacement seqwish induction.
    ///
    /// A value of 0 disables this filter. This applies to both AllWave-produced
    /// PAF and SweepGA-produced PAF before they are passed into seqwish.
    pub replacement_min_identity: f64,
    /// Plane-sweep mapping filter for pairwise replacement graph induction.
    pub replacement_num_mappings: String,
    /// Scaffold-chain filter for pairwise replacement graph induction.
    pub replacement_scaffold_filter: String,
    /// Minimum scaffold chain length (`min_scaffold_length` in the sweepga PAF
    /// filter) for pairwise replacement graph induction. Set to 0 to disable
    /// the scaffold-mass filter entirely. A non-zero value still passes
    /// through the adaptive `clamp_scaffold_params` shrink at sweepga's
    /// `filter_config_from_align_cfg` for short bubble-local inputs. The
    /// crush default is **0**: short bubble-local PAFs are dropped by the
    /// adaptive scaffold-mass clamp (≈ avg_seq_len × 3/5 ≈ 419 bp for
    /// plan-7-class inputs) even when the user already passed
    /// `no-filter=true`, so we explicitly disable the filter at the
    /// resolution-config layer. See `docs/crush-scaffold-mass-zero.md`.
    pub replacement_scaffold_mass: u64,
    /// SweepGA aligner backend used by `method=sweepga`.
    pub sweepga_aligner: String,
    /// SweepGA/FastGA k-mer frequency. A value of 0 uses a crush-local
    /// automatic frequency high enough to keep repeated bubble-local seeds.
    pub sweepga_kmer_frequency: usize,
    /// SweepGA minimum alignment length.
    pub sweepga_min_aln_length: u64,
    /// SweepGA wfmash percent identity. `None` lets wfmash auto-estimate.
    pub sweepga_map_pct_identity: Option<String>,
    /// Diagnostic selected-pair alignment budget for pairwise graph induction.
    ///
    /// Iterative multi-bubble mode logs when this is exceeded but still runs
    /// the replacement builder; size decides routing, not admission.
    pub max_pair_alignments: usize,
    /// Diagnostic PAF-size budget for one replacement.
    ///
    /// Iterative multi-bubble mode logs when this is exceeded but does not
    /// suppress graph induction.
    pub max_replacement_paf_bytes: usize,
    /// Skip SweepGA post-alignment filtering and feed selected pair alignments directly.
    pub sweepga_no_filter: bool,
    /// Force sparse pair dispatch for SweepGA/FastGA graph induction.
    ///
    /// FastGA does not have a cheap pairs-file mode here, so the default for
    /// FastGA is one all-vs-all batch followed by SweepGA filtering. Sparse
    /// pairs are still useful with wfmash, where the pairs file is honored by a
    /// single aligner invocation.
    pub sweepga_sparse_pairs: bool,
    /// Bubble-flanking context, in bp, included on each side of every bubble
    /// interior sequence when feeding the replacement aligner. A value of 0
    /// disables flanking (the aligner sees only the bubble interior). When
    /// non-zero, the aligner sees `<flank><interior><flank>` for each
    /// traversal; the flanking portion is clipped from the aligner's output
    /// graph before path-step substitution so the integration boundary remains
    /// the original bubble interior. The flank length per side is capped by
    /// the available path prefix/suffix bp so it never crosses the path edge.
    /// See `docs/crush-wider-context-bubbles.md`.
    pub replacement_flank_bp: usize,
    /// SPOA scoring parameters: (match, mismatch, gap_open, gap_ext, gap_open2, gap_ext2).
    pub scoring_params: (u8, u8, u8, u8, u8, u8),
    /// Legacy compression-ratio threshold retained for CLI compatibility.
    ///
    /// Crush replacement is now unconditional after graph-validity checks:
    /// this value is diagnostic only and never triggers an alternate aligner
    /// or keeps/replaces a plan based on compression.
    pub retry_min_compression_ratio: f64,
    /// Legacy input-bp floor paired with `retry_min_compression_ratio`.
    ///
    /// Also diagnostic only; it does not suppress or alter replacements.
    pub retry_min_input_bp: usize,
    /// Target root-path span, in bp, for `method=chain-greedy`.
    ///
    /// The greedy path-walk mode ignores the POVU parent/child tree when
    /// deciding replacement blocks: it sorts candidate bubbles by root-path
    /// adjacency and accumulates consecutive non-contained bubbles into one
    /// chain until adding the next bubble would exceed this span. A single
    /// larger bubble is still emitted as a one-bubble chain.
    pub chain_greedy_target_bp: usize,
    /// Candidate-source set for `method=iterative-multi-level`.
    pub multi_level_window_mode: MultiLevelWindowMode,
    /// Maximum root-path span for generated multi-level/window candidates.
    /// A value of 0 reuses `chain_greedy_target_bp`.
    pub multi_level_window_target_bp: usize,
    /// Maximum discovered flubble sites to merge into one generated window.
    ///
    /// The iterative multi-level mode is specifically meant to test whether
    /// the alignable unit is a run of neighboring flubbles rather than one POVU
    /// site. The default is therefore wider than the ordinary greedy-chain
    /// default, while the root-span and candidate-count caps remain the primary
    /// real-C4 runtime controls.
    pub multi_level_max_window_sites: usize,
    /// Maximum generated multi-level/window candidates to build per round.
    /// A value of 0 disables the cap.
    pub multi_level_candidate_limit: usize,
    /// Admission-only dry run for iterative multi-level/outward candidates.
    ///
    /// When set, candidate discovery, sorting, capping, and admission logging
    /// still run, but no replacement builder is invoked and the original graph
    /// is returned unchanged.
    pub multi_level_admission_only: bool,
    /// Diagnostic estimated transitive-closure work for one outward residual
    /// replacement. A value of 0 disables the warning.
    pub multi_level_max_transclosure_cells: usize,
    /// Diagnostic estimated POASTA dynamic-programming work for one outward
    /// residual replacement. A value of 0 disables the warning.
    pub multi_level_max_poasta_cells: usize,
    /// Diagnostic objective floor retained for CLI/report compatibility.
    ///
    /// Candidate application is gated only by exact path preservation and
    /// successful replacement construction. This value is logged as a
    /// would-have-failed diagnostic; it must not veto a replacement.
    pub multi_level_min_objective_delta: i128,
    /// Objective used to rank generated multi-bubble/window candidates.
    ///
    /// This never gates application after a replacement has been built; exact
    /// path preservation remains the acceptance gate.
    pub multi_level_objective: MultiLevelObjectiveMode,
    /// When true, candidate generation logs tiny high-frequency,
    /// low-complexity anchors as poor window boundaries.
    ///
    /// This is diagnostic only. Repeat-like boundaries must not veto a
    /// candidate; exact path preservation is the only correctness gate.
    pub multi_level_repeat_aware_boundaries: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReplacementMinMatchLenPolicy {
    Adaptive,
    Fixed(u64),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ResolutionMethod {
    Auto,
    /// Direct SPOA replacement.
    Poa,
    /// Direct POASTA replacement.
    Poasta,
    /// Debug/diagnostic resolver: align every traversal to one root with BiWFA
    /// and build a star-column graph. This is path-preserving but intentionally
    /// not the graph-quality default.
    StarBiwfa,
    /// Sparse many-to-many BiWFA with AllWave, seqwish induction, then SPOA polish.
    Allwave,
    /// SweepGA/FastGA or wfmash pair selection and graph induction, then SPOA polish.
    Sweepga,
    /// wfmash pair selection and graph induction, then SPOA polish. Equivalent to
    /// `method=sweepga,aligner=wfmash` but pins the aligner field so the documented
    /// intent survives any future config edit. Used to compare against PGGB's
    /// wfmash-based induction without disturbing the existing `Sweepga` path's
    /// `fastga` default. See `docs/crush-wfmash-replacement.md`.
    Wfmash,
    /// Depth-based routing: level 0 (top-level bubbles) → sweepga+seqwish,
    /// level ≥ 1 (every sub-bubble) → POASTA. No size threshold; the choice
    /// of aligner is tied to the provenance-tree depth, not to median
    /// traversal length. See `docs/crush-hierarchical.md`.
    Hierarchical,
    /// Greedy reference-path walk that groups consecutive bubble candidates
    /// into path-adjacent chains and aligns each whole chain with POASTA.
    ///
    /// This is a selection strategy rather than a new aligner: it deliberately
    /// bypasses the POVU provenance tree in `resolve_graph_bubbles` and uses
    /// only root-path adjacency to form replacement blocks.
    ChainGreedy,
    /// POVU-subtree chain blocking: walk the current POVU tree top-down and
    /// replace the largest subtree-rooted blocks whose parent interior fits
    /// under the chain block cap. The block replacement runs a local
    /// smoothxg-style pass followed by bounded POASTA cleanup.
    ChainPovu,
    /// Top-level POVU flubble blocking: select only level-0/root flubbles and
    /// run one local all-vs-all SweepGA/seqwish replacement for each complete
    /// bounded region. Descendant flubbles are absorbed into the root traversal
    /// sequences; adjacent roots are not merged.
    TopFlubbleSweepga,
    /// Aggressive experimental search policy: generate contained candidate
    /// regions from several POVU views (top-level sites, same-parent runs,
    /// path-order sliding windows, cross-level windows, local stringy
    /// neighborhoods, and explicitly outward-expanded residual windows), run
    /// multi-site windows through SweepGA/seqwish before trusting smaller POVU
    /// boundaries, then apply non-overlapping path-valid candidates. Objective
    /// metrics rank and diagnose candidates; they do not veto accepted
    /// replacements.
    IterativeMultiLevel,
    /// Coverage-driven multi-bubble pass: use the iterative multi-level window
    /// generator, include outward residual windows by default, run multi-site
    /// windows through SweepGA/seqwish, and report bp-weighted node path
    /// coverage / singleton-bp deltas diagnostically.
    CoverageMultiBubble,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum MultiLevelWindowMode {
    /// One-at-a-time parent/frontier descent: select the largest balanced
    /// parent flubble, replace it, re-decompose the whole graph, then repeat.
    Largest,
    /// Parent/frontier chunks: top-level POVU sites and parent-descendant blocks.
    Parent,
    /// Same-parent/same-level runs and parent-descendant windows.
    Sibling,
    /// Reference path-order sliding windows and local stringy neighborhoods.
    Sliding,
    /// Windows expanded outward from residual stringy/high-bp sites.
    Outward,
    /// Legacy combined search: sibling, sliding, cross-level, and local stringy windows.
    #[default]
    Combined,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum MultiLevelObjectiveMode {
    /// Historical iterative multi-level objective: reduce graph size.
    #[default]
    Size,
    /// Prefer higher bp-weighted node path coverage and lower singleton bp.
    Coverage,
}

impl MultiLevelWindowMode {
    pub fn parse_name(value: &str) -> Option<Self> {
        match value.replace('_', "-").to_ascii_lowercase().as_str() {
            "largest" | "biggest" | "frontier-largest" | "parent-largest" | "one-at-a-time"
            | "one-by-one" | "single" => Some(Self::Largest),
            "parent" | "parents" | "frontier" | "parent-first" | "top-level" | "toplevel" => {
                Some(Self::Parent)
            }
            "sibling" | "siblings" | "same-parent" | "same-level" => Some(Self::Sibling),
            "sliding" | "path-order" | "path" | "windows" | "sliding-windows" => {
                Some(Self::Sliding)
            }
            "outward" | "outward-expanded" | "residual" | "residual-outward" => Some(Self::Outward),
            "combined" | "all" | "multi" | "multi-level" | "both" => Some(Self::Combined),
            _ => None,
        }
    }
}

impl MultiLevelObjectiveMode {
    pub fn parse_name(value: &str) -> Option<Self> {
        match value.replace('_', "-").to_ascii_lowercase().as_str() {
            "size" | "graph-size" | "segment-size" | "segment-bp" | "compression" => {
                Some(Self::Size)
            }
            "coverage"
            | "node-coverage"
            | "path-coverage"
            | "node-path-coverage"
            | "bp-weighted-coverage"
            | "bp-weighted-node-coverage" => Some(Self::Coverage),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ResolutionPolishMethod {
    /// Re-decompose the replacement with POVU and run direct SPOA on nested flubbles.
    Poa,
    /// Re-decompose the replacement with POVU and run POASTA on nested flubbles.
    Poasta,
    /// Run the smoothxg-style sorted block smoother on the replacement graph.
    Smooth,
}

impl ResolutionPolishMethod {
    pub fn parse_name(value: &str) -> Option<Self> {
        match value.replace('_', "-").to_ascii_lowercase().as_str() {
            "poa" | "spoa" | "flubble-poa" | "bubble-poa" => Some(Self::Poa),
            "poasta" | "flubble-poasta" | "bubble-poasta" => Some(Self::Poasta),
            "smooth" | "smoothxg" | "pggb" | "pggb-smooth" => Some(Self::Smooth),
            _ => None,
        }
    }
}

impl ResolutionMethod {
    pub fn method_name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Poa => "SPOA",
            Self::Poasta => "POASTA",
            Self::StarBiwfa => "BiWFA/in-memory",
            Self::Allwave => "AllWave/seqwish",
            Self::Sweepga => "SweepGA/seqwish",
            Self::Wfmash => "wfmash/seqwish",
            Self::Hierarchical => "hierarchical",
            Self::ChainGreedy => "chain-greedy/POASTA",
            Self::ChainPovu => "chain-povu smoothxg→POASTA",
            Self::TopFlubbleSweepga => "top-flubble-sweepga",
            Self::IterativeMultiLevel => "iterative-multi-level",
            Self::CoverageMultiBubble => "coverage-multi-bubble",
        }
    }

    pub fn parse_name(value: &str) -> Option<Self> {
        match value.replace('_', "-").to_ascii_lowercase().as_str() {
            "auto" => Some(Self::Auto),
            "poa" | "spoa" => Some(Self::Poa),
            "poasta" => Some(Self::Poasta),
            "star-biwfa" | "biwfa-star" | "biwfa" | "wfa" => Some(Self::StarBiwfa),
            "allwave" | "aw" => Some(Self::Allwave),
            "sweepga" | "sw" => Some(Self::Sweepga),
            "wfmash" | "wf" => Some(Self::Wfmash),
            "hierarchical" | "hier" => Some(Self::Hierarchical),
            "chain-greedy" | "greedy-chain" | "greedy-walk" | "chain-walk" | "path-chain"
            | "path-walk-chain" | "crush-chain-greedy" => Some(Self::ChainGreedy),
            "chain-povu" | "povu-chain" | "povu-tree-blocks" | "tree-blocks"
            | "crush-chain-povu" => Some(Self::ChainPovu),
            "top-flubble-sweepga"
            | "top-level-flubble-sweepga"
            | "top-flubble"
            | "top-level-flubble"
            | "flubble-sweepga"
            | "chain-povu-sweepga"
            | "povu-sweepga"
            | "crush-top-flubble-sweepga" => Some(Self::TopFlubbleSweepga),
            "iterative-multi-level"
            | "multi-bubble"
            | "multi-bubbles"
            | "multi-flubble"
            | "multi-flubbles"
            | "flubble-window"
            | "flubble-windows"
            | "bubble-window"
            | "bubble-windows"
            | "multi-level"
            | "multi-level-window"
            | "multi-level-windows"
            | "window-crush"
            | "crush-to-size"
            | "aggressive-window"
            | "aggressive-windows" => Some(Self::IterativeMultiLevel),
            "coverage-multi-bubble"
            | "coverage-driven"
            | "coverage-window"
            | "coverage-windows"
            | "coverage-crush"
            | "node-coverage"
            | "path-coverage"
            | "coverage-repeat"
            | "repeat-aware-coverage" => Some(Self::CoverageMultiBubble),
            _ => None,
        }
    }
}

/// One parent frontier round is the fast default. Each selected parent
/// replacement is internally polished until no bounded POASTA candidates remain.
pub const DEFAULT_MAX_ITERATIONS: usize = 1;
/// By default, do not cap by rooted path span.
///
/// `max_bubble_span` is a POVU root-path coordinate guard. The root is currently
/// the first GFA path, so this is intentionally not part of the default runtime
/// budget. Traversal length/count/total sequence caps apply to direct local
/// replacement methods; pairwise induction methods are accepted after exact
/// path validation.
pub const DEFAULT_MAX_BUBBLE_SPAN: usize = 0;
pub const DEFAULT_MAX_TRAVERSAL_LEN: usize = 10_000;
pub const DEFAULT_MIN_TRAVERSAL_LEN: usize = 0;
pub const DEFAULT_MAX_MEDIAN_TRAVERSAL_LEN: usize = 1_000;
pub const DEFAULT_MAX_TOTAL_SEQUENCE: usize = 1_000_000;
pub const DEFAULT_MAX_TRAVERSALS: usize = 10_000;
pub const DEFAULT_AUTO_SPOA_MAX_TRAVERSAL_LEN: usize = 1_000;
pub const DEFAULT_AUTO_POASTA_MAX_TRAVERSAL_LEN: usize = 10_000;
pub const DEFAULT_CHAIN_POVU_MAX_BLOCK_BP: usize = 10_000;
const DEFAULT_CHAIN_POVU_POASTA_ITERATIONS: usize = 2;
pub const DEFAULT_AUTO_ALLWAVE_MAX_TOTAL_SEQUENCE: usize = 200_000;
pub const DEFAULT_AUTO_ALLWAVE_MAX_TRAVERSALS: usize = 128;
pub const DEFAULT_POLISH_ITERATIONS: usize = usize::MAX;
pub const DEFAULT_POLISH_MAX_TRAVERSAL_LEN: usize = 10_000;
pub const DEFAULT_POLISH_MAX_MEDIAN_TRAVERSAL_LEN: usize = 1_000;
pub const DEFAULT_POLISH_MAX_TOTAL_SEQUENCE: usize = 1_000_000;
pub const DEFAULT_POLISH_MAX_TRAVERSALS: usize = 10_000;
pub const DEFAULT_PAIR_K_NEAREST: usize = 3;
pub const DEFAULT_PAIR_K_FARTHEST: usize = 1;
pub const DEFAULT_PAIR_TREE_COUNT: usize = 1;
pub const DEFAULT_PAIR_RANDOM_FRACTION: f64 = 0.01;
pub const DEFAULT_PAIR_MASH_K: usize = 15;
pub const DEFAULT_REPLACEMENT_SEQWISH_MIN_MATCH_LEN: u64 = 311;
pub const DEFAULT_REPLACEMENT_MIN_MAP_LENGTH: u64 = 0;
pub const DEFAULT_REPLACEMENT_MIN_IDENTITY: f64 = 0.0;
/// Default `replacement_scaffold_mass` is 0 — disable the sweepga
/// scaffold-mass filter for bubble-local replacement induction. See
/// `docs/crush-scaffold-mass-zero.md` and `ResolutionConfig::replacement_scaffold_mass`.
pub const DEFAULT_REPLACEMENT_SCAFFOLD_MASS: u64 = 0;
pub const DEFAULT_SWEEPGA_KMER_FREQUENCY: usize = 0;
const MIN_AUTO_SWEEPGA_KMER_FREQUENCY: usize = 1_000;
const AUTO_SWEEPGA_KMER_FREQUENCY_PER_TRAVERSAL: usize = 10;
pub const DEFAULT_SWEEPGA_MIN_ALN_LENGTH: u64 = 0;
pub const DEFAULT_MAX_PAIR_ALIGNMENTS: usize = 10_000;
pub const DEFAULT_MAX_REPLACEMENT_PAF_BYTES: usize = 64 * 1024 * 1024;
/// Compression-ratio replacement retry is disabled. The config field is
/// retained as a diagnostic threshold for older engine strings that still pass
/// `retry-min-compression-ratio`, but it must not affect replacement choice.
pub const DEFAULT_RETRY_MIN_COMPRESSION_RATIO: f64 = 0.0;
/// Legacy input-bp floor for compression-ratio diagnostics.
pub const DEFAULT_RETRY_MIN_INPUT_BP: usize = 1_000;
/// Default `replacement_flank_bp` is 0 — no flanking context is added. The
/// canonical "wider-context" run sets this to 500. See
/// `docs/crush-wider-context-bubbles.md`.
pub const DEFAULT_REPLACEMENT_FLANK_BP: usize = 0;
/// Default greedy chain span for `method=chain-greedy`.
pub const DEFAULT_CHAIN_GREEDY_TARGET_BP: usize = 10_000;
pub const DEFAULT_MULTI_LEVEL_MAX_WINDOW_SITES: usize = 8;
pub const DEFAULT_MULTI_LEVEL_CANDIDATE_LIMIT: usize = 192;
pub const DEFAULT_MULTI_LEVEL_MIN_OBJECTIVE_DELTA: i128 = 1;
pub const DEFAULT_MULTI_LEVEL_REPEAT_AWARE_BOUNDARIES: bool = false;
impl Default for ResolutionConfig {
    fn default() -> Self {
        Self {
            max_iterations: DEFAULT_MAX_ITERATIONS,
            method: ResolutionMethod::Auto,
            auto_spoa_max_traversal_len: DEFAULT_AUTO_SPOA_MAX_TRAVERSAL_LEN,
            auto_poasta_max_traversal_len: DEFAULT_AUTO_POASTA_MAX_TRAVERSAL_LEN,
            auto_allwave_max_total_sequence: DEFAULT_AUTO_ALLWAVE_MAX_TOTAL_SEQUENCE,
            auto_allwave_max_traversals: DEFAULT_AUTO_ALLWAVE_MAX_TRAVERSALS,
            max_bubble_span: DEFAULT_MAX_BUBBLE_SPAN,
            max_traversal_len: DEFAULT_MAX_TRAVERSAL_LEN,
            min_traversal_len: DEFAULT_MIN_TRAVERSAL_LEN,
            max_median_traversal_len: DEFAULT_MAX_MEDIAN_TRAVERSAL_LEN,
            max_total_sequence: DEFAULT_MAX_TOTAL_SEQUENCE,
            max_traversals: DEFAULT_MAX_TRAVERSALS,
            polish_iterations: DEFAULT_POLISH_ITERATIONS,
            polish_max_traversal_len: DEFAULT_POLISH_MAX_TRAVERSAL_LEN,
            polish_max_median_traversal_len: DEFAULT_POLISH_MAX_MEDIAN_TRAVERSAL_LEN,
            polish_max_total_sequence: DEFAULT_POLISH_MAX_TOTAL_SEQUENCE,
            polish_max_traversals: DEFAULT_POLISH_MAX_TRAVERSALS,
            polish_method: ResolutionPolishMethod::Poa,
            pair_k_nearest: DEFAULT_PAIR_K_NEAREST,
            pair_k_farthest: DEFAULT_PAIR_K_FARTHEST,
            pair_tree_count: DEFAULT_PAIR_TREE_COUNT,
            pair_random_fraction: DEFAULT_PAIR_RANDOM_FRACTION,
            pair_mash_k: DEFAULT_PAIR_MASH_K,
            replacement_seqwish_min_match_len: DEFAULT_REPLACEMENT_SEQWISH_MIN_MATCH_LEN,
            replacement_min_match_len_policy: ReplacementMinMatchLenPolicy::Fixed(1),
            replacement_min_map_length: DEFAULT_REPLACEMENT_MIN_MAP_LENGTH,
            replacement_min_identity: DEFAULT_REPLACEMENT_MIN_IDENTITY,
            replacement_num_mappings: "1:1".to_string(),
            replacement_scaffold_filter: "1:1".to_string(),
            replacement_scaffold_mass: DEFAULT_REPLACEMENT_SCAFFOLD_MASS,
            sweepga_aligner: "fastga".to_string(),
            sweepga_kmer_frequency: DEFAULT_SWEEPGA_KMER_FREQUENCY,
            sweepga_min_aln_length: DEFAULT_SWEEPGA_MIN_ALN_LENGTH,
            sweepga_map_pct_identity: None,
            max_pair_alignments: DEFAULT_MAX_PAIR_ALIGNMENTS,
            max_replacement_paf_bytes: DEFAULT_MAX_REPLACEMENT_PAF_BYTES,
            sweepga_no_filter: false,
            sweepga_sparse_pairs: false,
            replacement_flank_bp: DEFAULT_REPLACEMENT_FLANK_BP,
            scoring_params: (1, 4, 6, 2, 26, 1),
            retry_min_compression_ratio: DEFAULT_RETRY_MIN_COMPRESSION_RATIO,
            retry_min_input_bp: DEFAULT_RETRY_MIN_INPUT_BP,
            chain_greedy_target_bp: DEFAULT_CHAIN_GREEDY_TARGET_BP,
            multi_level_window_mode: MultiLevelWindowMode::Combined,
            multi_level_window_target_bp: 0,
            multi_level_max_window_sites: DEFAULT_MULTI_LEVEL_MAX_WINDOW_SITES,
            multi_level_candidate_limit: DEFAULT_MULTI_LEVEL_CANDIDATE_LIMIT,
            multi_level_admission_only: false,
            multi_level_max_transclosure_cells: 0,
            multi_level_max_poasta_cells: 0,
            multi_level_min_objective_delta: DEFAULT_MULTI_LEVEL_MIN_OBJECTIVE_DELTA,
            multi_level_objective: MultiLevelObjectiveMode::Size,
            multi_level_repeat_aware_boundaries: DEFAULT_MULTI_LEVEL_REPEAT_AWARE_BOUNDARIES,
        }
    }
}

/// Summary of a resolution run.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ResolutionStats {
    /// Number of frontier replacement rounds attempted.
    pub iterations: usize,
    pub candidates_seen: usize,
    pub resolved: usize,
    pub bailed: usize,
    /// Deprecated compatibility counter. Compression-ratio retry is disabled
    /// and this remains zero.
    pub retry_attempts: usize,
    /// Deprecated compatibility counter. Compression-ratio retry is disabled
    /// and this remains zero.
    pub retry_wins: usize,
    /// Deprecated compatibility counter. Compression-ratio retry is disabled
    /// and this remains zero.
    pub retry_failures: usize,
}

/// Resolved GFA plus run statistics.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ResolvedGfa {
    pub gfa: String,
    pub stats: ResolutionStats,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
struct Step {
    node: usize,
    rev: bool,
}

#[derive(Clone, Debug)]
struct Segment {
    id: String,
    seq: Vec<u8>,
}

#[derive(Clone, Debug)]
struct Path {
    name: String,
    steps: Vec<Step>,
}

#[derive(Clone, Debug)]
struct Graph {
    segments: Vec<Segment>,
    paths: Vec<Path>,
}

#[derive(Clone, Debug, Default)]
struct PathRange {
    path_idx: usize,
    begin_step: usize,
    end_step: usize,
    /// Bubble interior sequence (the bytes between `begin_step` and `end_step`
    /// on this path). This is what the integration code substitutes.
    sequence: Vec<u8>,
    /// Flank-extended sequence presented to the replacement aligner:
    /// `<left_flank><interior><right_flank>`. Empty when flanking is disabled;
    /// in that case the aligner sees `sequence` directly. The flank lengths
    /// are tracked in `left_flank_bp` / `right_flank_bp` so the aligner's
    /// output graph can be clipped back to the interior before substitution.
    extended_sequence: Vec<u8>,
    /// Left-flank length in bp (0 when no flank was added).
    left_flank_bp: usize,
    /// Right-flank length in bp (0 when no flank was added).
    right_flank_bp: usize,
}

#[derive(Clone, Debug)]
struct BubbleCandidate {
    ranges: Vec<PathRange>,
    signature: String,
    root_start_step: usize,
    root_end_step: usize,
    root_span: usize,
    total_steps: usize,
    unique_steps: usize,
    traversal_stats: TraversalStats,
    /// Depth in the provenance tree. 0 = top-level (root) bubble, ≥1 = sub-bubble.
    /// Populated at discovery from POVU's reported level and refined by the
    /// tree-driven round loop to the authoritative tree depth. Used by
    /// `ResolutionMethod::Hierarchical` to dispatch by depth instead of size.
    level: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct GraphQuality {
    segments: usize,
    segment_bp: usize,
    links: usize,
    path_steps: usize,
    node_coverage_mean: f64,
    node_coverage_bp_weighted_mean: f64,
    node_coverage_p10: usize,
    node_coverage_median: usize,
    node_coverage_p90: usize,
    singleton_nodes: usize,
    singleton_bp: usize,
    high_coverage_threshold: usize,
    high_coverage_nodes: usize,
    high_coverage_bp: usize,
    trivial_stringy: usize,
    path_white_space_bp_total: u64,
    path_white_space_bp_p99: usize,
    path_white_space_bp_max: usize,
    path_white_space_long_bridges: usize,
    score: u64,
}

#[derive(Clone, Debug, Default)]
struct CandidateFrontier {
    selected: Vec<BubbleCandidate>,
    sites_seen: usize,
    candidates_seen: usize,
}

#[derive(Clone, Debug)]
struct ChainCandidate {
    candidate: BubbleCandidate,
    source_bubbles: usize,
}

#[derive(Clone, Debug, Default)]
struct ChainFrontier {
    selected: Vec<ChainCandidate>,
    sites_seen: usize,
    discovered_candidates: usize,
    path_walk_candidates: usize,
    chains_formed: usize,
}

/// State of a bubble in the explicit bubble tree (Phase 6 true level descent).
#[derive(Clone, Debug)]
enum BubbleState {
    Unresolved,
    Resolved { at_round: usize },
    Failed { at_round: usize },
}

/// One node in the explicit bubble tree. The tree is grown top-down: roots
/// come from POVU on the input; children come from POVU on each resolved
/// node's *local replacement subgraph* (not from re-POVU on the whole working
/// graph). This is the architectural property `docs/crush-architecture-spec.md`
/// §Phase-6 calls for.
#[derive(Clone, Debug)]
struct BubbleNode {
    parent: Option<usize>,
    children: Vec<usize>,
    level: usize,
    state: BubbleState,
    /// bp-keyed signature — stable across rewrites of OTHER bubbles
    bp_signature: String,
    /// Reference-path bp range (begin, end-exclusive). Stable across rewrites.
    ref_bp_begin: usize,
    ref_bp_end: usize,
    /// POVU site id (string token like `<43388736>43388746`) at the time the
    /// node was discovered. Used for diagnostics; not stable across rewrites.
    discovery_site_id: String,
    /// Local POVU subbubble keys discovered after this node was resolved.
    /// Populated when `state == Resolved`; populated even if the result is
    /// empty so we can distinguish "no children" from "not yet resolved".
    local_child_bp_keys: Option<Vec<String>>,
    /// Max per-path traversal length at initial discovery — used by the
    /// "skip too-big roots, descend to children" heuristic so big-bubble
    /// graphs like full C4 GRCh38 don't feed megabase-scale roots into the
    /// aligner. 0 means "unknown / not applicable" (e.g., child nodes added
    /// by local POVU don't carry this).
    initial_traversal_max_len: usize,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct TraversalStats {
    count: usize,
    min_len: usize,
    median_len: usize,
    p90_len: usize,
    max_len: usize,
    total_len: usize,
}

#[derive(Clone, Debug)]
struct ReplacementPlan {
    candidate: BubbleCandidate,
    replacement: Graph,
}

#[derive(Clone, Debug)]
struct PathReplacement {
    begin_step: usize,
    end_step: usize,
    steps: Vec<OutStep>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
enum OutNode {
    Original(usize),
    Replacement(usize, usize),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
struct OutStep {
    node: OutNode,
    rev: bool,
}

/// Resolve simple bounded bubbles in a blunt GFA while preserving every emitted
/// path sequence exactly.
///
/// This first implementation supports GFA 1.0 `S`, `L`, and `P` records. Links
/// must be blunt (`0M`) because the re-rendered graph derives links from path
/// adjacency and emits `0M` edges. Use `syng:blunt` output, not raw syng overlap
/// GFA, as input.
pub fn resolve_gfa_bubbles(gfa: &str, config: &ResolutionConfig) -> io::Result<ResolvedGfa> {
    let parse_start = Instant::now();
    log::info!(
        "crush: parsing input GFA ({} bytes) before bubble decomposition",
        gfa.len()
    );
    let graph = parse_gfa(gfa)?;
    log::info!(
        "crush: parsed input GFA into {} segment(s), {} path(s) in {:.2?}",
        graph.segments.len(),
        graph.paths.len(),
        parse_start.elapsed()
    );
    resolve_graph_bubbles(graph, Some(gfa), config, true)
}

fn resolve_graph_bubbles(
    mut graph: Graph,
    original_gfa: Option<&str>,
    config: &ResolutionConfig,
    emit_logs: bool,
) -> io::Result<ResolvedGfa> {
    if config.method == ResolutionMethod::ChainGreedy {
        return resolve_graph_bubble_chains(graph, original_gfa, config, emit_logs);
    }

    let mut stats = ResolutionStats::default();
    // Provenance: signatures of bubbles already transitioned to Resolved.
    // Hard invariant — Phase 6: no tree node is ever re-resolved. The
    // signature is bp-based (see `candidate_signature`) and therefore stable
    // across rewrites of other bubbles.
    let mut resolved_signatures: FxHashSet<String> = FxHashSet::default();
    let mut changed = false;
    // Hard invariant — Phase 6: fresh node IDs are STRICTLY MONOTONIC across
    // the whole resolution run.
    let mut next_id = initial_next_segment_id(&graph);
    let initial_next_id = next_id;
    if emit_logs {
        log::info!(
            "crush: initial next_id={} (input has {} segment(s), max integer id={})",
            next_id,
            graph.segments.len(),
            next_id.saturating_sub(1)
        );
        let (match_score, mismatch, gap_open1, gap_extend1, gap_open2, gap_extend2) =
            config.scoring_params;
        log::info!(
            "crush config: method={:?}, polish={:?}, poa-scoring={},{},{},{},{},{}",
            config.method,
            config.polish_method,
            match_score,
            mismatch,
            gap_open1,
            gap_extend1,
            gap_open2,
            gap_extend2
        );
    }

    if config.method == ResolutionMethod::ChainPovu {
        return resolve_graph_bubbles_chain_povu(
            graph,
            original_gfa,
            config,
            emit_logs,
            stats,
            next_id,
            initial_next_id,
        );
    }
    if config.method == ResolutionMethod::TopFlubbleSweepga {
        return resolve_graph_bubbles_top_flubble_sweepga(
            graph,
            original_gfa,
            config,
            emit_logs,
            stats,
            next_id,
            initial_next_id,
        );
    }
    if matches!(
        config.method,
        ResolutionMethod::IterativeMultiLevel | ResolutionMethod::CoverageMultiBubble
    ) {
        return resolve_graph_bubbles_iterative_multi_level(
            graph,
            original_gfa,
            config,
            emit_logs,
            stats,
            next_id,
            initial_next_id,
        );
    }

    // -------------- True level descent (Phase 6) --------------
    //
    // We maintain an explicit mutable bubble tree (`tree`) seeded from a
    // single POVU pass on the input graph. Each round:
    //   1. The active frontier is the unresolved tree nodes whose parent was
    //      resolved in the previous round (or, for round 1, every root).
    //   2. Each frontier node's *current* working-graph step coordinates are
    //      recovered by running global POVU once and matching by bp-keyed
    //      signature (the bp-key is stable across rewrites — see
    //      `candidate_signature`).
    //   3. Replacements are built in parallel and applied as one batch (the
    //      strike-and-link contract of Phase 4–5 is unchanged).
    //   4. For each successful replacement, **local POVU on the replacement
    //      subgraph** discovers the node's children (`discover_local_subbubble_keys`).
    //      Those keys seed the round-N+1 frontier.
    //
    // This eliminates the prior "global re-POVU + is_leaf filter" pattern:
    // sub-bubble discovery is now LOCALIZED to each resolved bubble's
    // replacement graph, and the descent walks the tree top-down instead of
    // bottom-up. The validation gate "round N: K resolved (from sites
    // discovered by re-POVU on round N-1 local subgraphs)" is logged below.

    let (mut tree, roots) = build_initial_bubble_tree(&graph, config, emit_logs)?;
    let global_path_names: Vec<String> = graph.paths.iter().map(|p| p.name.clone()).collect();

    // bp_signature → tree node index. We track the tree explicitly even
    // though the round-N+1 frontier is selected by bp-region containment
    // (see below): the tree gives us state (Resolved/Unresolved/Failed),
    // provenance assertions, and the parent/child structure that the
    // local-POVU discovery step grows.
    let mut bp_to_node: FxHashMap<String, usize> = FxHashMap::default();
    for (idx, node) in tree.iter().enumerate() {
        bp_to_node.insert(node.bp_signature.clone(), idx);
    }

    // Each round, the active frontier is determined from the tree (children
    // of round-(N-1) resolved nodes). Round 1 starts with the roots. Children
    // are added to the tree by:
    //   - Initial POVU (parent_id wiring in build_initial_bubble_tree)
    //   - Local POVU on each resolved replacement (discover_local_subbubble_keys)
    // We track per-round resolved bp regions as a fallback admission rule
    // (bp-region containment) for sub-bubbles the initial POVU missed but
    // that local POVU on the replacement surfaced.
    let mut last_round_resolved_ref_bp_regions: Vec<(usize, usize)> = Vec::new();
    let mut last_round_resolved_nodes: Vec<usize> = Vec::new();
    let mut per_round_frontier_sizes: Vec<usize> = Vec::new();
    let mut per_round_local_povu_child_counts: Vec<usize> = Vec::new();

    for round in 0..config.max_iterations {
        // Termination: if this isn't round 1 and the previous round resolved
        // nothing, the tree is exhausted.
        if round > 0 && last_round_resolved_ref_bp_regions.is_empty() {
            if emit_logs {
                log::info!(
                    "crush round {}: previous round resolved nothing — tree fully descended",
                    round + 1
                );
            }
            break;
        }
        let round_start = Instant::now();
        let before_quality = graph_quality(&graph);

        // ---- discovery: global POVU on the current working graph ----
        let discovery_start = Instant::now();
        let (all_discovered, sites_seen, timings) =
            discover_all_candidates(&graph, config, emit_logs)?;
        // Compute the set of bp_signatures that are "active" this round.
        //
        // Active = unresolved tree nodes that are "ready" to be processed:
        //   - All unresolved roots (rounds N>=2 still include any sibling
        //     root that was rejected by non-overlap in round N-1)
        //   - PLUS children of round-(N-1) resolved nodes (the standard
        //     top-down descent step from `docs/crush-architecture-spec.md`
        //     §Phase-6)
        //
        // bp_signatures are stable across rewrites because path sequences
        // are byte-preserved (Phase 6 invariant), so a tree node's bp_signature
        // still matches whichever global-POVU site re-discovers the same
        // bubble after a sibling rewrite.
        let mut active_node_indices: Vec<usize> = Vec::new();
        for &r in &roots {
            if matches!(tree[r].state, BubbleState::Unresolved) {
                active_node_indices.push(r);
            }
        }
        if round > 0 {
            for &p in &last_round_resolved_nodes {
                for &c in &tree[p].children {
                    if matches!(tree[c].state, BubbleState::Unresolved)
                        && !active_node_indices.contains(&c)
                    {
                        active_node_indices.push(c);
                    }
                }
            }
        }
        // "Skip and descend" guard: for nodes whose initial-POVU traversal-max
        // exceeds the configured max_traversal_len AND which have children in
        // the tree (initial POVU saw nested structure), DON'T feed the parent
        // into the aligner — that triggers the F3-F4 sweepga cascade on full
        // C4. Instead, mark the parent Resolved (state=Resolved with at_round=0)
        // and descend immediately to its children for THIS round. The small-
        // bubble test fixture is unaffected (its L0 root has max_len=409,
        // well below max_traversal_len=10000), but C4's megabase-scale roots
        // get fanned out into their nested leaves on round 1.
        if round == 0 && config.max_traversal_len > 0 {
            let skip_threshold = config.max_traversal_len.saturating_mul(2);
            let mut expanded: Vec<usize> = Vec::new();
            for &n in &active_node_indices {
                if tree[n].initial_traversal_max_len > skip_threshold
                    && !tree[n].children.is_empty()
                {
                    // Mark this oversized root as "skipped — descended to children"
                    // so the provenance / state-tracking still makes sense.
                    tree[n].state = BubbleState::Failed { at_round: 0 };
                    if emit_logs {
                        log::info!(
                            "crush round 1: skipping oversized root site_id={} (initial max_traversal_len={} > {}); descending to {} child(ren)",
                            tree[n].discovery_site_id,
                            tree[n].initial_traversal_max_len,
                            skip_threshold,
                            tree[n].children.len()
                        );
                    }
                    for &c in &tree[n].children {
                        if matches!(tree[c].state, BubbleState::Unresolved)
                            && !expanded.contains(&c)
                        {
                            expanded.push(c);
                        }
                    }
                } else {
                    expanded.push(n);
                }
            }
            active_node_indices = expanded;
            // Recursively descend if children are also oversized.
            loop {
                let mut next: Vec<usize> = Vec::new();
                let mut any_skipped = false;
                for &n in &active_node_indices {
                    if tree[n].initial_traversal_max_len > skip_threshold
                        && !tree[n].children.is_empty()
                    {
                        tree[n].state = BubbleState::Failed { at_round: 0 };
                        if emit_logs {
                            log::info!(
                                "crush round 1: skipping oversized descendant site_id={} (max_traversal_len={} > {}); descending to {} child(ren)",
                                tree[n].discovery_site_id,
                                tree[n].initial_traversal_max_len,
                                skip_threshold,
                                tree[n].children.len()
                            );
                        }
                        any_skipped = true;
                        for &c in &tree[n].children {
                            if matches!(tree[c].state, BubbleState::Unresolved)
                                && !next.contains(&c)
                            {
                                next.push(c);
                            }
                        }
                    } else {
                        next.push(n);
                    }
                }
                active_node_indices = next;
                if !any_skipped {
                    break;
                }
            }
        }
        let active_bp_keys: FxHashSet<String> = active_node_indices
            .iter()
            .map(|&i| tree[i].bp_signature.clone())
            .collect();
        // Admit ONLY sites whose bp_signature matches an active tree node.
        // The tree was grown by initial POVU (on the input) plus local POVU
        // on each resolved replacement; sub-bubbles that didn't enter the
        // tree are not admitted, which keeps the descent strictly bounded
        // by the bubble forest. This implements the spec's "POVU on local
        // subgraph (not on the whole rewritten working graph)" admission
        // rule.
        let active_candidates: Vec<BubbleCandidate> = all_discovered
            .into_iter()
            .filter(|d| active_bp_keys.contains(&d.candidate.signature))
            .filter(|d| !resolved_signatures.contains(&d.candidate.signature))
            .map(|d| {
                let mut c = d.candidate;
                // The tree's level is authoritative — it's grown from POVU's
                // parent_id at initial discovery and from local POVU on each
                // resolved replacement (see `discover_local_subbubble_keys`).
                // Used by `ResolutionMethod::Hierarchical` to dispatch by depth.
                if let Some(&idx) = bp_to_node.get(&c.signature) {
                    c.level = tree[idx].level;
                }
                c
            })
            .collect();
        let candidates_seen = active_candidates.len();
        let frontier = finalize_frontier(
            &graph,
            config,
            active_candidates,
            sites_seen,
            candidates_seen,
            emit_logs,
            timings,
            0,
        )?;
        let discovery_elapsed = discovery_start.elapsed();
        if frontier.selected.is_empty() {
            if emit_logs {
                log::info!(
                    "crush round {}: no eligible candidates from {} POVU site(s) in {:.2?} (round {} resolved {} region(s))",
                    round + 1,
                    frontier.sites_seen,
                    discovery_elapsed,
                    round,
                    last_round_resolved_ref_bp_regions.len()
                );
            }
            break;
        }
        stats.iterations = round + 1;

        if emit_logs {
            log::info!(
                "crush round {}: {} POVU site(s) on working graph, {} candidate(s) (round {}), {} selected in {:.2?}",
                round + 1,
                frontier.sites_seen,
                candidates_seen,
                if round == 0 {
                    "1 / initial POVU on input → roots".to_string()
                } else {
                    format!(
                        "{} / strictly inside round {} resolved bp regions",
                        round + 1,
                        round
                    )
                },
                frontier.selected.len(),
                discovery_elapsed
            );
            log::info!(
                "crush round {} traversal stats: {}",
                round + 1,
                format_candidate_length_summary("selected", &frontier.selected)
            );
        }

        for _ in &frontier.selected {
            stats.candidates_seen += 1;
        }

        let selected_count = frontier.selected.len();
        let build_start = Instant::now();
        let build_started = AtomicUsize::new(0);
        let build_progress = AtomicUsize::new(0);
        // Precompute pre-apply path positions for the working graph — used by
        // `discover_local_subbubble_keys` to translate local POVU sub-site
        // coordinates back to global bp keys.
        let pre_apply_path_positions: Vec<Vec<usize>> = (0..graph.paths.len())
            .map(|i| path_positions(&graph, i))
            .collect();
        let build_results = frontier
            .selected
            .into_par_iter()
            .map(|candidate| {
                let method = candidate_replacement_method(&candidate, config);
                if emit_logs && selected_count <= 32 {
                    let started = build_started.fetch_add(1, Ordering::Relaxed) + 1;
                    log::info!(
                        "crush round {}: building replacement {}/{} with {:?}; traversals={}, max-len={}, median-len={}, total-len={}, root-span={}",
                        round + 1,
                        started,
                        selected_count,
                        method,
                        candidate.traversal_stats.count,
                        candidate.traversal_stats.max_len,
                        candidate.traversal_stats.median_len,
                        candidate.traversal_stats.total_len,
                        candidate.root_span
                    );
                }
                let result = (|| {
                    let replacement = build_replacement_with_method(&candidate, config, method)?;
                    if replacement.segments.is_empty() {
                        return Ok::<Option<ReplacementPlan>, io::Error>(None);
                    }
                    Ok::<Option<ReplacementPlan>, io::Error>(Some(ReplacementPlan {
                        candidate,
                        replacement,
                    }))
                })();
                if emit_logs {
                    let done = build_progress.fetch_add(1, Ordering::Relaxed) + 1;
                    if selected_count <= 32 || done == selected_count || done % 25 == 0 {
                        let status = replacement_build_status(&result);
                        log::info!(
                            "crush round {}: replacement build progress {}/{} ({})",
                            round + 1,
                            done,
                            selected_count,
                            status
                        );
                    }
                }
                result
            })
            .collect::<Vec<_>>();
        let build_elapsed = build_start.elapsed();

        let mut plans: Vec<ReplacementPlan> = Vec::new();
        let mut empty = 0usize;
        let mut failed = 0usize;
        let mut path_invalid = 0usize;
        for result in build_results {
            match result {
                Ok(Some(plan)) => {
                    plans.push(plan);
                }
                Ok(None) => {
                    empty += 1;
                    stats.bailed += 1;
                }
                Err(err) => {
                    if replacement_failure_is_path_invalid(&err) {
                        path_invalid += 1;
                    } else {
                        failed += 1;
                    }
                    log::debug!(
                        "resolution: candidate replacement failed in round {}: {}",
                        round + 1,
                        err
                    );
                    stats.bailed += 1;
                }
            }
        }

        if plans.is_empty() {
            if emit_logs {
                log::info!(
                    "crush round {}: 0/{} replacement(s) accepted in {:.2?} (failed={}, empty={}, path-invalid={})",
                    round + 1,
                    selected_count,
                    build_elapsed,
                    failed,
                    empty,
                    path_invalid
                );
            }
            last_round_resolved_ref_bp_regions.clear();
            continue;
        }

        if emit_logs && (failed > 0 || empty > 0 || path_invalid > 0) {
            log::info!(
                "crush round {}: accepted {}/{} replacement(s) in {:.2?}; failed={}, empty={}, path-invalid={}",
                round + 1,
                plans.len(),
                selected_count,
                build_elapsed,
                failed,
                empty,
                path_invalid
            );
        }

        log_replacement_compression_diagnostics(&plans, config, round, emit_logs);

        let resolved_count = plans.len();
        // Phase 6 hard invariant: no tree node is ever re-resolved.
        for plan in &plans {
            let sig = &plan.candidate.signature;
            if !resolved_signatures.insert(sig.clone()) {
                panic!(
                    "crush provenance assertion: bubble signature {sig:?} was \
                     resolved twice (current round {}); flubble tree node \
                     re-resolution is forbidden by Phase 6 invariants",
                    round + 1
                );
            }
        }

        // ---- local POVU on each resolved replacement to discover children ----
        // The local POVU pass is the AUTHORITATIVE source of which sub-bubbles
        // exist inside each just-resolved bubble (`docs/crush-architecture-spec.md`
        // §Phase-6). Its discoveries are wired into the tree as children of
        // the resolved node; the round-N+1 frontier is selected by bp-region
        // containment against the same set of resolved bp regions, which is
        // semantically equivalent to "POVU on the local subgraph" but
        // sidesteps bp-key matching subtleties.
        let local_povu_start = Instant::now();
        let local_children_per_plan: Vec<Vec<String>> = plans
            .par_iter()
            .map(|plan| {
                discover_local_subbubble_keys(
                    &plan.candidate,
                    &plan.replacement,
                    &pre_apply_path_positions,
                    &global_path_names,
                    config,
                )
                .unwrap_or_default()
            })
            .collect();
        let local_povu_elapsed = local_povu_start.elapsed();

        // Wire children into the tree (parent → child edges) and record the
        // round's resolved bp regions for the round-N+1 containment filter.
        let mut just_resolved_nodes: Vec<usize> = Vec::new();
        let mut just_resolved_ref_bp_regions: Vec<(usize, usize)> = Vec::new();
        let mut total_local_children = 0usize;
        let mut total_local_children_already_known = 0usize;
        let mut total_local_children_new = 0usize;
        let pre_apply_ref_positions = pre_apply_path_positions
            .first()
            .cloned()
            .unwrap_or_default();
        for (plan, child_keys) in plans.iter().zip(local_children_per_plan.iter()) {
            // Find or create the parent tree node by bp-signature. For round 1
            // this exists (it was a root from initial POVU). For deeper rounds
            // the candidate may not have a pre-existing tree node (it came
            // from bp-region containment on a global-POVU site that doesn't
            // perfectly match a prior local-POVU-discovered child key); in
            // that case we create a fresh node so the provenance check still
            // works.
            let parent_nidx = if let Some(&n) = bp_to_node.get(&plan.candidate.signature) {
                n
            } else {
                let parent_level = if round == 0 { 0 } else { round };
                let new_idx = tree.len();
                tree.push(BubbleNode {
                    parent: None,
                    children: Vec::new(),
                    level: parent_level,
                    state: BubbleState::Unresolved,
                    bp_signature: plan.candidate.signature.clone(),
                    ref_bp_begin: 0,
                    ref_bp_end: 0,
                    discovery_site_id: String::new(),
                    local_child_bp_keys: None,
                    initial_traversal_max_len: plan.candidate.traversal_stats.max_len,
                });
                bp_to_node.insert(plan.candidate.signature.clone(), new_idx);
                new_idx
            };
            tree[parent_nidx].state = BubbleState::Resolved {
                at_round: round + 1,
            };
            just_resolved_nodes.push(parent_nidx);

            // Record the parent's reference bp region.
            let pb = plan.candidate.root_start_step;
            let pe = plan.candidate.root_end_step;
            if pe + 1 < pre_apply_ref_positions.len() {
                let bp_begin = pre_apply_ref_positions[pb];
                let bp_end = pre_apply_ref_positions[pe + 1];
                just_resolved_ref_bp_regions.push((bp_begin, bp_end));
                tree[parent_nidx].ref_bp_begin = bp_begin;
                tree[parent_nidx].ref_bp_end = bp_end;
            }

            let mut stored_keys: Vec<String> = Vec::new();
            for child_key in child_keys {
                total_local_children += 1;
                let child_nidx = if let Some(&existing) = bp_to_node.get(child_key) {
                    total_local_children_already_known += 1;
                    existing
                } else {
                    let new_idx = tree.len();
                    tree.push(BubbleNode {
                        parent: Some(parent_nidx),
                        children: Vec::new(),
                        level: tree[parent_nidx].level + 1,
                        state: BubbleState::Unresolved,
                        bp_signature: child_key.clone(),
                        ref_bp_begin: 0,
                        ref_bp_end: 0,
                        discovery_site_id: String::new(),
                        local_child_bp_keys: None,
                        initial_traversal_max_len: 0,
                    });
                    bp_to_node.insert(child_key.clone(), new_idx);
                    total_local_children_new += 1;
                    new_idx
                };
                if tree[child_nidx].parent.is_none() {
                    tree[child_nidx].parent = Some(parent_nidx);
                }
                if !tree[parent_nidx].children.contains(&child_nidx) {
                    tree[parent_nidx].children.push(child_nidx);
                }
                stored_keys.push(child_key.clone());
            }
            tree[parent_nidx].local_child_bp_keys = Some(stored_keys);
        }
        last_round_resolved_nodes = just_resolved_nodes;
        last_round_resolved_ref_bp_regions = just_resolved_ref_bp_regions;

        let rewrite_start = Instant::now();
        let pre_apply_next_id = next_id;
        let next_graph = apply_replacement_frontier(&graph, &plans, &mut next_id)?;
        assert!(
            next_id >= pre_apply_next_id,
            "crush fresh-id invariant: next_id rewound from {pre_apply_next_id} to {next_id}"
        );
        let rewrite_elapsed = rewrite_start.elapsed();
        let after_quality = graph_quality(&next_graph);
        graph = next_graph;
        changed = true;
        stats.resolved += resolved_count;
        per_round_frontier_sizes.push(resolved_count);
        per_round_local_povu_child_counts.push(total_local_children);

        if emit_logs {
            log::info!(
                "crush round {}: {} resolved (from sites discovered by re-POVU on round {} local subgraphs); local POVU yielded {} child key(s) ({} new + {} reused) in {:.2?}",
                round + 1,
                resolved_count,
                if round == 0 { "initial POVU on input".to_string() } else { format!("{}", round) },
                total_local_children,
                total_local_children_new,
                total_local_children_already_known,
                local_povu_elapsed
            );
            log::info!(
                "crush round {}: resolved {}/{} replacement(s) in {:.2?}; rewrite+validate {:.2?}; total {:.2?}; ids minted {}..{}; quality {} -> {}",
                round + 1,
                resolved_count,
                selected_count,
                build_elapsed,
                rewrite_elapsed,
                round_start.elapsed(),
                pre_apply_next_id,
                next_id,
                before_quality.summary(),
                after_quality.summary()
            );
        }
    }
    if emit_logs && !per_round_frontier_sizes.is_empty() {
        let frontier_summary = per_round_frontier_sizes
            .iter()
            .enumerate()
            .map(|(idx, n)| format!("r{}={}", idx + 1, n))
            .collect::<Vec<_>>()
            .join(", ");
        log::info!(
            "crush per-round frontier sizes (Phase 6 true level descent): [{}]; total resolved={}; next_id moved {} -> {}",
            frontier_summary,
            stats.resolved,
            initial_next_id,
            next_id
        );
        let _ = last_round_resolved_nodes; // future: per-round tree-walk diagnostics
    }

    let final_gfa = if changed {
        render_graph(&graph)
    } else {
        original_gfa
            .map(str::to_string)
            .unwrap_or_else(|| render_graph(&graph))
    };
    Ok(ResolvedGfa {
        gfa: final_gfa,
        stats,
    })
}

/// Greedy path-walk chain mode.
///
/// This intentionally does not use `build_initial_bubble_tree` or POVU
/// parent/leaf metadata to decide replacement blocks. POVU is still used to
/// surface candidate bubble boundaries, but the block decision is made by
/// walking the root path in coordinate order and chaining consecutive
/// non-contained bubbles until the configured root-span cap is reached.
fn resolve_graph_bubble_chains(
    mut graph: Graph,
    original_gfa: Option<&str>,
    config: &ResolutionConfig,
    emit_logs: bool,
) -> io::Result<ResolvedGfa> {
    let mut stats = ResolutionStats::default();
    let mut resolved_signatures: FxHashSet<String> = FxHashSet::default();
    let mut changed = false;
    let mut next_id = initial_next_segment_id(&graph);
    let initial_next_id = next_id;

    if emit_logs {
        let (match_score, mismatch, gap_open1, gap_extend1, gap_open2, gap_extend2) =
            config.scoring_params;
        log::info!(
            "crush chain-greedy: target-bp={}, max-rounds={}, replacement=POASTA, poa-scoring={},{},{},{},{},{}",
            config.chain_greedy_target_bp,
            config.max_iterations,
            match_score,
            mismatch,
            gap_open1,
            gap_extend1,
            gap_open2,
            gap_extend2
        );
    }

    let mut per_round_chain_counts: Vec<usize> = Vec::new();
    for round in 0..config.max_iterations {
        let round_start = Instant::now();
        let before_quality = graph_quality(&graph);
        let frontier =
            find_path_walk_chain_frontier(&graph, config, &resolved_signatures, emit_logs)?;

        if frontier.selected.is_empty() {
            if emit_logs {
                log::info!(
                    "crush chain-greedy round {}: no chain blocks selected from {} POVU site(s), {} candidate bubble(s), {} path-walk bubble(s), {} formed chain(s)",
                    round + 1,
                    frontier.sites_seen,
                    frontier.discovered_candidates,
                    frontier.path_walk_candidates,
                    frontier.chains_formed
                );
            }
            break;
        }

        stats.iterations = round + 1;
        stats.candidates_seen += frontier.selected.len();
        let selected_count = frontier.selected.len();

        if emit_logs {
            log::info!(
                "crush chain-greedy round {}: {} POVU site(s), {} candidate bubble(s), {} path-walk bubble(s), {} chain(s) formed, {} selected",
                round + 1,
                frontier.sites_seen,
                frontier.discovered_candidates,
                frontier.path_walk_candidates,
                frontier.chains_formed,
                selected_count
            );
            log::info!(
                "crush chain-greedy round {} chain-size distribution: {}",
                round + 1,
                format_chain_candidate_summary("selected", &frontier.selected)
            );
        }

        let build_start = Instant::now();
        let build_started = AtomicUsize::new(0);
        let build_progress = AtomicUsize::new(0);
        let build_results = frontier
            .selected
            .into_par_iter()
            .map(|chain| {
                let source_bubbles = chain.source_bubbles;
                let candidate = chain.candidate;
                if emit_logs && selected_count <= 64 {
                    let started = build_started.fetch_add(1, Ordering::Relaxed) + 1;
                    log::info!(
                        "crush chain-greedy round {}: building chain {}/{} with POASTA; bubbles={}, traversals={}, max-len={}, median-len={}, total-len={}, root-span={}",
                        round + 1,
                        started,
                        selected_count,
                        source_bubbles,
                        candidate.traversal_stats.count,
                        candidate.traversal_stats.max_len,
                        candidate.traversal_stats.median_len,
                        candidate.traversal_stats.total_len,
                        candidate.root_span
                    );
                }
                let result = (|| {
                    let replacement =
                        build_replacement_with_method(&candidate, config, ResolutionMethod::ChainGreedy)?;
                    if replacement.segments.is_empty() {
                        return Ok::<Option<ReplacementPlan>, io::Error>(None);
                    }
                    Ok::<Option<ReplacementPlan>, io::Error>(Some(ReplacementPlan {
                        candidate,
                        replacement,
                    }))
                })();
                if emit_logs {
                    let done = build_progress.fetch_add(1, Ordering::Relaxed) + 1;
                    if selected_count <= 64 || done == selected_count || done % 25 == 0 {
                        let status = replacement_build_status(&result);
                        log::info!(
                            "crush chain-greedy round {}: replacement build progress {}/{} ({})",
                            round + 1,
                            done,
                            selected_count,
                            status
                        );
                    }
                }
                result
            })
            .collect::<Vec<_>>();
        let build_elapsed = build_start.elapsed();

        let mut plans: Vec<ReplacementPlan> = Vec::new();
        let mut failed_or_empty = 0usize;
        for result in build_results {
            match result {
                Ok(Some(plan)) => plans.push(plan),
                Ok(None) => {
                    failed_or_empty += 1;
                    stats.bailed += 1;
                }
                Err(err) => {
                    failed_or_empty += 1;
                    log::debug!(
                        "resolution: chain-greedy replacement failed in round {}: {}",
                        round + 1,
                        err
                    );
                    stats.bailed += 1;
                }
            }
        }

        if plans.is_empty() {
            if emit_logs {
                log::info!(
                    "crush chain-greedy round {}: 0/{} chain replacement(s) accepted in {:.2?} ({} failed or empty)",
                    round + 1,
                    selected_count,
                    build_elapsed,
                    failed_or_empty
                );
            }
            break;
        }

        for plan in &plans {
            resolved_signatures.insert(plan.candidate.signature.clone());
        }

        let resolved_count = plans.len();
        let rewrite_start = Instant::now();
        let pre_apply_next_id = next_id;
        let next_graph = apply_replacement_frontier(&graph, &plans, &mut next_id)?;
        assert!(
            next_id >= pre_apply_next_id,
            "crush fresh-id invariant: next_id rewound from {pre_apply_next_id} to {next_id}"
        );
        let rewrite_elapsed = rewrite_start.elapsed();
        let after_quality = graph_quality(&next_graph);
        graph = next_graph;
        changed = true;
        stats.resolved += resolved_count;
        per_round_chain_counts.push(resolved_count);

        if emit_logs {
            if failed_or_empty > 0 {
                log::info!(
                    "crush chain-greedy round {}: accepted {}/{} chain replacement(s) in {:.2?}; {} failed or empty",
                    round + 1,
                    resolved_count,
                    selected_count,
                    build_elapsed,
                    failed_or_empty
                );
            }
            log::info!(
                "crush chain-greedy round {}: resolved {}/{} chain replacement(s) in {:.2?}; rewrite+validate {:.2?}; total {:.2?}; ids minted {}..{}; quality {} -> {}",
                round + 1,
                resolved_count,
                selected_count,
                build_elapsed,
                rewrite_elapsed,
                round_start.elapsed(),
                pre_apply_next_id,
                next_id,
                before_quality.summary(),
                after_quality.summary()
            );
        }
    }

    if emit_logs && !per_round_chain_counts.is_empty() {
        let summary = per_round_chain_counts
            .iter()
            .enumerate()
            .map(|(idx, n)| format!("r{}={}", idx + 1, n))
            .collect::<Vec<_>>()
            .join(", ");
        log::info!(
            "crush chain-greedy per-round chains: [{}]; total resolved={}; next_id moved {} -> {}",
            summary,
            stats.resolved,
            initial_next_id,
            next_id
        );
    }

    let final_gfa = if changed {
        render_graph(&graph)
    } else {
        original_gfa
            .map(str::to_string)
            .unwrap_or_else(|| render_graph(&graph))
    };
    Ok(ResolvedGfa {
        gfa: final_gfa,
        stats,
    })
}

fn resolve_graph_bubbles_chain_povu(
    mut graph: Graph,
    original_gfa: Option<&str>,
    config: &ResolutionConfig,
    emit_logs: bool,
    mut stats: ResolutionStats,
    mut next_id: usize,
    initial_next_id: usize,
) -> io::Result<ResolvedGfa> {
    let mut changed = false;
    let mut resolved_signatures: FxHashSet<String> = FxHashSet::default();
    let mut per_round_block_counts: Vec<usize> = Vec::new();
    let mut per_round_block_sizes: Vec<Vec<usize>> = Vec::new();

    for round in 0..config.max_iterations {
        let round_start = Instant::now();
        let before_quality = graph_quality(&graph);
        let discovery_start = Instant::now();
        let (all_discovered, sites_seen, timings) =
            discover_all_candidates(&graph, config, emit_logs)?;
        let selection = select_chain_povu_blocks(
            &all_discovered,
            DEFAULT_CHAIN_POVU_MAX_BLOCK_BP,
            &resolved_signatures,
        );
        let candidates_seen = selection.selected.len();
        let frontier = finalize_frontier(
            &graph,
            config,
            selection.selected,
            sites_seen,
            candidates_seen,
            emit_logs,
            timings,
            0,
        )?;
        let discovery_elapsed = discovery_start.elapsed();
        if frontier.selected.is_empty() {
            if emit_logs {
                log::info!(
                    "crush chain-povu round {}: no eligible subtree block(s) from {} POVU site(s) in {:.2?}; target_block_bp={}, tree_nodes={}, roots={}, skipped_seen={}, oversize_internal={}, oversize_leaves={}",
                    round + 1,
                    frontier.sites_seen,
                    discovery_elapsed,
                    DEFAULT_CHAIN_POVU_MAX_BLOCK_BP,
                    selection.tree_nodes,
                    selection.roots,
                    selection.skipped_seen,
                    selection.oversize_internal,
                    selection.oversize_leaves
                );
            }
            break;
        }

        stats.iterations = round + 1;
        let selected_count = frontier.selected.len();
        let selected_sizes = frontier
            .selected
            .iter()
            .map(chain_povu_block_bp)
            .collect::<Vec<_>>();
        stats.candidates_seen += selected_count;

        if emit_logs {
            log::info!(
                "crush chain-povu round {}: {} POVU site(s), target_block_bp={}, tree_nodes={}, roots={}, selected {} subtree block(s) ({} descendant site(s) absorbed, {} already-seen skipped, {} oversized internal node(s), {} oversized leaf/leaves) in {:.2?}; block-size {}",
                round + 1,
                frontier.sites_seen,
                DEFAULT_CHAIN_POVU_MAX_BLOCK_BP,
                selection.tree_nodes,
                selection.roots,
                selected_count,
                selection.merged_descendant_sites,
                selection.skipped_seen,
                selection.oversize_internal,
                selection.oversize_leaves,
                discovery_elapsed,
                format_bp_distribution("selected", &selected_sizes)
            );
            log::info!(
                "crush chain-povu round {} traversal stats: {}",
                round + 1,
                format_candidate_length_summary("selected", &frontier.selected)
            );
        }

        let build_start = Instant::now();
        let build_started = AtomicUsize::new(0);
        let build_progress = AtomicUsize::new(0);
        reset_chain_povu_replacement_decisions();
        let build_results = frontier
            .selected
            .into_par_iter()
            .map(|candidate| {
                if emit_logs && selected_count <= 32 {
                    let started = build_started.fetch_add(1, Ordering::Relaxed) + 1;
                    log::info!(
                        "crush chain-povu round {}: building block {}/{} with smoothxg→POASTA; traversals={}, max-len={}, median-len={}, total-len={}, root-span={}, block-bp={}",
                        round + 1,
                        started,
                        selected_count,
                        candidate.traversal_stats.count,
                        candidate.traversal_stats.max_len,
                        candidate.traversal_stats.median_len,
                        candidate.traversal_stats.total_len,
                        candidate.root_span,
                        chain_povu_block_bp(&candidate)
                    );
                }
                let result = (|| {
                    let replacement =
                        build_replacement_with_method(&candidate, config, ResolutionMethod::ChainPovu)?;
                    if replacement.segments.is_empty() {
                        return Ok::<Option<ReplacementPlan>, io::Error>(None);
                    }
                    Ok::<Option<ReplacementPlan>, io::Error>(Some(ReplacementPlan {
                        candidate,
                        replacement,
                    }))
                })();
                if emit_logs {
                    let done = build_progress.fetch_add(1, Ordering::Relaxed) + 1;
                    if selected_count <= 32 || done == selected_count || done % 25 == 0 {
                        let status = replacement_build_status(&result);
                        log::info!(
                            "crush chain-povu round {}: block build progress {}/{} ({})",
                            round + 1,
                            done,
                            selected_count,
                            status
                        );
                    }
                }
                result
            })
            .collect::<Vec<_>>();
        let build_elapsed = build_start.elapsed();
        let replacement_decisions = chain_povu_replacement_decisions();

        let mut plans: Vec<ReplacementPlan> = Vec::new();
        let mut failed_or_empty = 0usize;
        for result in build_results {
            match result {
                Ok(Some(plan)) => plans.push(plan),
                Ok(None) => {
                    failed_or_empty += 1;
                    stats.bailed += 1;
                }
                Err(err) => {
                    failed_or_empty += 1;
                    stats.bailed += 1;
                    log::debug!(
                        "resolution: chain-povu block replacement failed in round {}: {}",
                        round + 1,
                        err
                    );
                }
            }
        }

        if plans.is_empty() {
            if emit_logs {
                log::info!(
                    "crush chain-povu round {}: 0/{} block replacement(s) accepted in {:.2?} ({} failed or produced empty replacements)",
                    round + 1,
                    selected_count,
                    build_elapsed,
                    failed_or_empty
                );
            }
            break;
        }

        let resolved_count = plans.len();
        for plan in &plans {
            resolved_signatures.insert(plan.candidate.signature.clone());
        }

        let rewrite_start = Instant::now();
        let pre_apply_next_id = next_id;
        let next_graph = apply_replacement_frontier(&graph, &plans, &mut next_id)?;
        assert!(
            next_id >= pre_apply_next_id,
            "crush fresh-id invariant: next_id rewound from {pre_apply_next_id} to {next_id}"
        );
        let rewrite_elapsed = rewrite_start.elapsed();
        let after_quality = graph_quality(&next_graph);
        graph = next_graph;
        changed = true;
        stats.resolved += resolved_count;

        let accepted_sizes = plans
            .iter()
            .map(|plan| chain_povu_block_bp(&plan.candidate))
            .collect::<Vec<_>>();
        per_round_block_counts.push(resolved_count);
        per_round_block_sizes.push(accepted_sizes.clone());

        if emit_logs {
            log::info!(
                "crush chain-povu round {}: accepted {}/{} block replacement(s) in {:.2?}; rewrite+validate {:.2?}; total {:.2?}; ids minted {}..{}; accepted-size {}; replacement path: smooth_used={}, direct_on_empty={}, direct_on_smooth_failure={}; quality {} -> {}",
                round + 1,
                resolved_count,
                selected_count,
                build_elapsed,
                rewrite_elapsed,
                round_start.elapsed(),
                pre_apply_next_id,
                next_id,
                format_bp_distribution("accepted", &accepted_sizes),
                replacement_decisions.smooth_used,
                replacement_decisions.direct_on_empty,
                replacement_decisions.direct_on_smooth_failure,
                before_quality.summary(),
                after_quality.summary()
            );
        }
    }

    if emit_logs {
        if per_round_block_counts.is_empty() {
            log::info!(
                "crush chain-povu per-iteration block counts: []; total resolved=0; next_id moved {} -> {}",
                initial_next_id,
                next_id
            );
        } else {
            let count_summary = per_round_block_counts
                .iter()
                .enumerate()
                .map(|(idx, n)| format!("r{}={}", idx + 1, n))
                .collect::<Vec<_>>()
                .join(", ");
            let size_summary = per_round_block_sizes
                .iter()
                .enumerate()
                .map(|(idx, sizes)| {
                    format!("r{} {}", idx + 1, format_bp_distribution("blocks", sizes))
                })
                .collect::<Vec<_>>()
                .join("; ");
            log::info!(
                "crush chain-povu per-iteration block counts: [{}]; size distributions: [{}]; total resolved={}; next_id moved {} -> {}",
                count_summary,
                size_summary,
                stats.resolved,
                initial_next_id,
                next_id
            );
        }
    }

    let final_gfa = if changed {
        render_graph(&graph)
    } else {
        original_gfa
            .map(str::to_string)
            .unwrap_or_else(|| render_graph(&graph))
    };
    Ok(ResolvedGfa {
        gfa: final_gfa,
        stats,
    })
}

fn resolve_graph_bubbles_top_flubble_sweepga(
    mut graph: Graph,
    original_gfa: Option<&str>,
    config: &ResolutionConfig,
    emit_logs: bool,
    mut stats: ResolutionStats,
    mut next_id: usize,
    initial_next_id: usize,
) -> io::Result<ResolvedGfa> {
    let before_quality = graph_quality(&graph);
    let discovery_start = Instant::now();
    let (all_discovered, sites_seen, timings) = discover_all_candidates(&graph, config, emit_logs)?;
    let selection = select_top_flubble_sweepga_regions(&all_discovered, &FxHashSet::default());
    let selected_sizes = selection
        .selected
        .iter()
        .map(chain_povu_block_bp)
        .collect::<Vec<_>>();
    let candidates_seen = selection.selected.len();
    let frontier = finalize_frontier(
        &graph,
        config,
        selection.selected,
        sites_seen,
        candidates_seen,
        emit_logs,
        timings,
        0,
    )?;
    let discovery_elapsed = discovery_start.elapsed();

    if frontier.selected.is_empty() {
        if emit_logs {
            log::info!(
                "crush top-flubble-sweepga: no materialized top-level region(s) selected from {} POVU site(s), {} top-level flubble(s) in {:.2?}; no cross-top-level merging performed",
                sites_seen,
                selection.top_level_regions,
                discovery_elapsed
            );
            log::info!(
                "crush top-flubble-sweepga per-round regions: []; total resolved=0; next_id moved {} -> {}",
                initial_next_id,
                next_id
            );
        }
        let final_gfa = original_gfa
            .map(str::to_string)
            .unwrap_or_else(|| render_graph(&graph));
        return Ok(ResolvedGfa {
            gfa: final_gfa,
            stats,
        });
    }

    stats.iterations = 1;
    stats.candidates_seen += frontier.selected.len();
    let selected_count = frontier.selected.len();

    if emit_logs {
        log::info!(
            "crush top-flubble-sweepga: {} POVU site(s), {} top-level flubble(s), {} descendant/non-root site(s), {} sequence-variable top-level region(s) selected in {:.2?}; pre-materialize size {}; no cross-top-level merging performed",
            sites_seen,
            selection.top_level_regions,
            selection.descendant_sites,
            selected_count,
            discovery_elapsed,
            format_bp_distribution("top-level", &selected_sizes)
        );
        log::info!(
            "crush top-flubble-sweepga traversal stats: {}",
            format_candidate_length_summary("selected", &frontier.selected)
        );
    }

    let build_start = Instant::now();
    let mut selected_regions = frontier.selected;
    selected_regions.sort_by_key(|candidate| {
        (
            candidate.traversal_stats.total_len,
            candidate.traversal_stats.max_len,
            candidate.root_span,
        )
    });
    let mut build_results = Vec::with_capacity(selected_count);
    for (idx, candidate) in selected_regions.into_iter().enumerate() {
        let region_ordinal = idx + 1;
        if emit_logs && (region_ordinal <= 10 || region_ordinal % 25 == 0) {
            log::info!(
                "crush top-flubble-sweepga: building top-level region {}/{}; traversals={}, input-bp={}, max-len={}, median-len={}, root-span={}",
                region_ordinal,
                selected_count,
                candidate.traversal_stats.count,
                candidate.traversal_stats.total_len,
                candidate.traversal_stats.max_len,
                candidate.traversal_stats.median_len,
                candidate.root_span
            );
        }
        let result = build_top_flubble_sweepga_replacement_with_evidence(&candidate, config).map(
            |(replacement, evidence)| TopFlubbleSweepgaBuild {
                region_ordinal,
                candidate,
                replacement,
                evidence,
            },
        );
        if emit_logs {
            let status = if result.is_ok() { "accepted" } else { "failed" };
            if region_ordinal <= 10 || region_ordinal == selected_count || region_ordinal % 25 == 0
            {
                log::info!(
                    "crush top-flubble-sweepga: region build progress {}/{} ({})",
                    region_ordinal,
                    selected_count,
                    status
                );
            }
        }
        build_results.push((region_ordinal, result));
    }
    let build_elapsed = build_start.elapsed();

    let mut builds = Vec::with_capacity(build_results.len());
    for (region_ordinal, result) in build_results {
        match result {
            Ok(build) => builds.push(build),
            Err(err) => {
                if is_zero_alignment_replacement_error(&err) {
                    stats.bailed += 1;
                    if emit_logs {
                        log::info!(
                            "crush top-flubble-sweepga: skipped top-level region {}/{} because alignment emitted no PAF records; leaving original graph interval intact: {}",
                            region_ordinal,
                            selected_count,
                            err
                        );
                    }
                    continue;
                }
                return Err(io::Error::other(format!(
                    "top-flubble-sweepga replacement failed validity/build: {err}"
                )));
            }
        }
    }
    if builds.is_empty() {
        if emit_logs {
            log::info!(
                "crush top-flubble-sweepga: selected {} top-level region(s), but every region emitted zero PAF records; no replacement graph was applied",
                selected_count
            );
        }
        return Ok(ResolvedGfa {
            gfa: render_graph(&graph),
            stats,
        });
    }

    let evidence_rows = builds
        .iter()
        .map(TopFlubbleSweepgaEvidenceRow::from_build)
        .collect::<Vec<_>>();
    let plans = builds
        .into_iter()
        .map(|build| ReplacementPlan {
            candidate: build.candidate,
            replacement: build.replacement,
        })
        .collect::<Vec<_>>();

    log_replacement_compression_diagnostics(&plans, config, 0, emit_logs);

    let rewrite_start = Instant::now();
    let pre_apply_next_id = next_id;
    let next_graph = apply_replacement_frontier(&graph, &plans, &mut next_id)?;
    assert!(
        next_id >= pre_apply_next_id,
        "crush fresh-id invariant: next_id rewound from {pre_apply_next_id} to {next_id}"
    );
    let rewrite_elapsed = rewrite_start.elapsed();
    let after_quality = graph_quality(&next_graph);
    graph = next_graph;
    stats.resolved += plans.len();

    if emit_logs {
        for row in &evidence_rows {
            log::info!(
                "crush top-flubble-sweepga region {}/{} evidence: input_traversals={}, input_bp={}, emitted_alignment_records={}, aligned_bp={}, paf_bytes={}, collinear_records={}, off_diagonal_records={}, dust_like_records={}, fastga_frequency={}, replacement_segments={}, replacement_shared_segments={}, replacement_bp={}, applied=true",
                row.region_ordinal,
                selected_count,
                row.input_traversals,
                row.input_bp,
                row.alignment_records,
                row.aligned_bp,
                row.paf_bytes,
                row.collinear_records,
                row.off_diagonal_records,
                row.dust_like_records,
                row.kmer_frequency,
                row.replacement_segments,
                row.replacement_shared_segments,
                row.replacement_bp
            );
        }
        let accepted_sizes = plans
            .iter()
            .map(|plan| chain_povu_block_bp(&plan.candidate))
            .collect::<Vec<_>>();
        log::info!(
            "crush top-flubble-sweepga: applied {}/{} top-level region replacement(s) in one round; build {:.2?}; rewrite+validate {:.2?}; ids minted {}..{}; accepted-size {}; quality {} -> {}",
            stats.resolved,
            selected_count,
            build_elapsed,
            rewrite_elapsed,
            pre_apply_next_id,
            next_id,
            format_bp_distribution("accepted", &accepted_sizes),
            before_quality.summary(),
            after_quality.summary()
        );
        log::info!(
            "crush top-flubble-sweepga per-round regions: [r1={}]; top-level flubble count={}; no cross-top-level merging; total resolved={}; next_id moved {} -> {}",
            stats.resolved,
            selection.top_level_regions,
            stats.resolved,
            initial_next_id,
            next_id
        );
    }

    Ok(ResolvedGfa {
        gfa: render_graph(&graph),
        stats,
    })
}

fn resolve_graph_bubbles_iterative_multi_level(
    mut graph: Graph,
    original_gfa: Option<&str>,
    config: &ResolutionConfig,
    emit_logs: bool,
    mut stats: ResolutionStats,
    mut next_id: usize,
    initial_next_id: usize,
) -> io::Result<ResolvedGfa> {
    let mut changed = false;
    let mut resolved_signatures: FxHashSet<String> = FxHashSet::default();
    let mut per_round_accepted = Vec::new();
    let objective_mode = effective_multi_level_objective(config);
    let repeat_aware_boundaries = effective_repeat_aware_boundaries(config);

    if emit_logs {
        log::info!(
            "crush iterative-multi-level: method={}, mode={:?}, target-bp={}, max-window-sites={}, effective-candidate-limit={}, admission-only={}, window-total-bp-cap={} (diagnostic only), outward-build-diagnostics=max-pair-alignments={},max-paf-bytes={},max-transclosure-cells={},max-poasta-cells={}, objective={:?}, min-objective-delta={} (diagnostic only), repeat-aware-boundaries={} (diagnostic only), replacement-routing=multi-site-auto-poasta-or-sweepga/single-site-auto, application_gate=exact-path-preservation",
            config.method.method_name(),
            effective_multi_level_window_mode(config),
            multi_level_target_bp(config),
            config.multi_level_max_window_sites,
            effective_multi_level_candidate_limit(config),
            config.multi_level_admission_only,
            multi_level_window_total_bp_cap(config),
            config.max_pair_alignments,
            config.max_replacement_paf_bytes,
            config.multi_level_max_transclosure_cells,
            config.multi_level_max_poasta_cells,
            objective_mode,
            config.multi_level_min_objective_delta,
            repeat_aware_boundaries
        );
    }

    for round in 0..config.max_iterations {
        let round_start = Instant::now();
        let before_quality = graph_quality(&graph);
        let discovery_start = Instant::now();
        let (all_discovered, sites_seen, timings) =
            discover_all_candidates(&graph, config, emit_logs)?;
        let generated =
            generate_multi_level_candidates(&graph, config, &all_discovered, &resolved_signatures);
        let discovery_elapsed = discovery_start.elapsed();

        if generated.candidates.is_empty() {
            if emit_logs {
                log::info!(
                    "crush iterative-multi-level round {}: no generated candidate(s) from {} POVU site(s), {} polymorphic site(s) in {:.2?}; before {}",
                    round + 1,
                    sites_seen,
                    all_discovered.len(),
                    discovery_elapsed,
                    before_quality.summary()
                );
            }
            break;
        }

        stats.iterations = round + 1;
        stats.candidates_seen += generated.candidates.len();

        if emit_logs {
            log::info!(
                "crush iterative-multi-level round {}: {} POVU site(s), {} polymorphic site(s), {} generated candidate(s), {} considered after cap in {:.2?}; sources {}; repeat-boundary-flagged={} (diagnostic only); total-bp-cap-exceeded={} (diagnostic only), max-len-cap-exceeded={} (diagnostic only), median-len-cap-exceeded={} (diagnostic only); discovery detail render {:.2?}, povu-parse {:.2?}, povu-decompose {:.2?}, id-map {:.2?}, path-index {:.2?}, candidate-build {:.2?}; before {}",
                round + 1,
                sites_seen,
                all_discovered.len(),
                generated.generated_total,
                generated.candidates.len(),
                discovery_elapsed,
                format_multi_level_source_counts(&generated.source_counts),
                generated.repeat_boundary_rejected,
                generated.total_bp_cap_exceeded,
                generated.max_len_cap_exceeded,
                generated.median_len_cap_exceeded,
                timings.render,
                timings.povu_parse,
                timings.povu_decompose,
                timings.id_map,
                timings.path_index,
                timings.candidate_build,
                before_quality.summary()
            );
            log::info!(
                "crush iterative-multi-level round {} traversal stats: {}",
                round + 1,
                format_multi_level_candidate_summary("generated", &generated.candidates)
            );
            log::info!(
                "crush iterative-multi-level round {} largest discovered residual sites: {}",
                round + 1,
                format_discovered_residual_sites(&graph, &all_discovered, 8)
            );
            log::info!(
                "crush iterative-multi-level round {} generated candidate detail: {}",
                round + 1,
                format_multi_level_candidate_details(&graph, &generated.candidates, 24)
            );
        }

        if config.multi_level_admission_only {
            if emit_logs {
                let admitted_source_counts =
                    multi_level_window_source_counts(&generated.candidates);
                log::info!(
                    "crush iterative-multi-level round {}: admission-only dry-run admitted {} candidate(s) from {}; replacement builders not invoked; before {}",
                    round + 1,
                    generated.candidates.len(),
                    format_multi_level_source_counts(&admitted_source_counts),
                    before_quality.summary()
                );
                log::info!(
                    "crush iterative-multi-level round {} admission detail: {}",
                    round + 1,
                    format_multi_level_candidate_details(&graph, &generated.candidates, 48)
                );
            }
            break;
        }

        let build_diagnostics = multi_level_build_budget_diagnostics(&generated.candidates, config);
        if emit_logs && !build_diagnostics.is_empty() {
            log::info!(
                "crush iterative-multi-level round {}: {} candidate(s) exceed configured build-budget diagnostics; all remain admitted for replacement build; detail: {}",
                round + 1,
                build_diagnostics.len(),
                format_multi_level_budget_diagnostics(&graph, &build_diagnostics, 24)
            );
        }

        let build_candidates = generated.candidates;
        let selected_count = build_candidates.len();
        if selected_count == 0 {
            if emit_logs {
                log::info!(
                    "crush iterative-multi-level round {}: 0 candidate(s) admitted to replacement build; budget_diagnostic={}; before {}",
                    round + 1,
                    build_diagnostics.len(),
                    before_quality.summary()
                );
            }
            break;
        }

        let build_start = Instant::now();
        let build_started = AtomicUsize::new(0);
        let build_progress = AtomicUsize::new(0);
        let build_parallelism =
            multi_level_replacement_build_parallelism(rayon::current_num_threads(), selected_count);
        if emit_logs {
            log::info!(
                "crush iterative-multi-level round {}: replacement build concurrency={} (global threads={}, candidates={})",
                round + 1,
                build_parallelism,
                rayon::current_num_threads(),
                selected_count
            );
        }
        let build_results = if build_parallelism <= 1 {
            build_candidates
                .into_iter()
                .map(|window| {
                    build_multi_level_window_candidate(
                        &graph,
                        config,
                        objective_mode,
                        emit_logs,
                        selected_count,
                        round + 1,
                        &build_started,
                        &build_progress,
                        window,
                    )
                })
                .collect::<Vec<_>>()
        } else {
            let build_pool = rayon::ThreadPoolBuilder::new()
                .num_threads(build_parallelism)
                .stack_size(16 * 1024 * 1024)
                .build()
                .map_err(|err| {
                    io::Error::other(format!(
                        "iterative multi-level replacement build pool: {err}"
                    ))
                })?;
            build_pool.install(|| {
                build_candidates
                    .into_par_iter()
                    .map(|window| {
                        build_multi_level_window_candidate(
                            &graph,
                            config,
                            objective_mode,
                            emit_logs,
                            selected_count,
                            round + 1,
                            &build_started,
                            &build_progress,
                            window,
                        )
                    })
                    .collect::<Vec<_>>()
            })
        };
        let build_elapsed = build_start.elapsed();

        let mut built = Vec::new();
        let mut failed_or_empty = 0usize;
        for result in build_results {
            match result {
                Ok(Some(candidate)) => built.push(candidate),
                Ok(None) => {
                    failed_or_empty += 1;
                    stats.bailed += 1;
                }
                Err(err) => {
                    failed_or_empty += 1;
                    stats.bailed += 1;
                    log::debug!(
                        "resolution: iterative-multi-level candidate failed in round {}: {}",
                        round + 1,
                        err
                    );
                }
            }
        }

        let built_count = built.len();
        let objective_floor = config.multi_level_min_objective_delta;
        let local_objective_below_floor = built
            .iter()
            .filter(|candidate| candidate.objective.score_delta < objective_floor)
            .count();
        let objective_passing_count = built_count.saturating_sub(local_objective_below_floor);
        built.sort_by(|a, b| {
            b.objective
                .score_delta
                .cmp(&a.objective.score_delta)
                .then_with(|| {
                    b.objective
                        .segment_bp_delta
                        .cmp(&a.objective.segment_bp_delta)
                })
                .then_with(|| b.objective.segment_delta.cmp(&a.objective.segment_delta))
                .then_with(|| {
                    a.plan
                        .candidate
                        .root_start_step
                        .cmp(&b.plan.candidate.root_start_step)
                })
                .then_with(|| {
                    a.plan
                        .candidate
                        .root_end_step
                        .cmp(&b.plan.candidate.root_end_step)
                })
        });

        let mut occupied_by_path: Vec<Vec<(usize, usize)>> = vec![Vec::new(); graph.paths.len()];
        let mut accepted = Vec::new();
        let mut overlap_deferred = 0usize;
        for candidate in built {
            if candidate_conflicts_with_occupied(&occupied_by_path, &candidate.plan.candidate) {
                overlap_deferred += 1;
                continue;
            }
            mark_candidate_occupied(&mut occupied_by_path, &candidate.plan.candidate);
            accepted.push(candidate);
        }

        if accepted.is_empty() {
            if emit_logs {
                log::info!(
                    "crush iterative-multi-level round {}: 0 applicable candidate(s) after build/non-overlap scheduling; built={}, objective_passing_diagnostic={}, objective_below_floor_diagnostic={}, overlap_deferred={}, failed_or_empty={}, budget_diagnostic={}, objective_floor={} (diagnostic only), build {:.2?}; before {}",
                    round + 1,
                    built_count,
                    objective_passing_count,
                    local_objective_below_floor,
                    overlap_deferred,
                    failed_or_empty,
                    build_diagnostics.len(),
                    objective_floor,
                    build_elapsed,
                    before_quality.summary()
                );
            }
            break;
        }

        let rewrite_start = Instant::now();
        let pre_apply_next_id = next_id;
        let plans = accepted
            .iter()
            .map(|candidate| candidate.plan.clone())
            .collect::<Vec<_>>();
        let mut trial_next_id = next_id;
        let next_graph = apply_replacement_frontier(&graph, &plans, &mut trial_next_id)?;
        assert!(
            trial_next_id >= pre_apply_next_id,
            "crush fresh-id invariant: next_id rewound from {pre_apply_next_id} to {trial_next_id}"
        );
        let after_quality = graph_quality(&next_graph);
        let graph_delta =
            graph_objective_delta_for_mode(before_quality, after_quality, objective_mode);
        let rewrite_elapsed = rewrite_start.elapsed();
        let accepted_sources = accepted.iter().fold(
            FxHashMap::<MultiLevelCandidateSource, usize>::default(),
            |mut acc, item| {
                *acc.entry(item.source).or_insert(0) += 1;
                acc
            },
        );
        let objective_summary = format_objective_delta_summary(&accepted);
        next_id = trial_next_id;

        for plan in &plans {
            resolved_signatures.insert(plan.candidate.signature.clone());
        }
        graph = next_graph;
        changed = true;
        stats.resolved += plans.len();
        per_round_accepted.push(plans.len());

        if emit_logs {
            log::info!(
                "crush iterative-multi-level round {}: applied {} candidate(s) from {}; built={}, objective_passing_diagnostic={}, objective_below_floor_diagnostic={}, overlap_deferred={}, budget_diagnostic={}, objective_floor={} (diagnostic only); local objective {}; global_objective_delta={} (diagnostic only); failed_or_empty={}; build {:.2?}; rewrite+validate {:.2?}; total {:.2?}; ids minted {}..{}; before {}; after {}",
                round + 1,
                plans.len(),
                format_multi_level_source_counts(&accepted_sources),
                built_count,
                objective_passing_count,
                local_objective_below_floor,
                overlap_deferred,
                build_diagnostics.len(),
                objective_floor,
                objective_summary,
                graph_delta,
                failed_or_empty,
                build_elapsed,
                rewrite_elapsed,
                round_start.elapsed(),
                pre_apply_next_id,
                next_id,
                before_quality.summary(),
                after_quality.summary()
            );
        }
    }

    if emit_logs {
        let summary = per_round_accepted
            .iter()
            .enumerate()
            .map(|(idx, n)| format!("r{}={}", idx + 1, n))
            .collect::<Vec<_>>()
            .join(", ");
        log::info!(
            "crush iterative-multi-level per-round accepted candidates: [{}]; total resolved={}; next_id moved {} -> {}",
            summary,
            stats.resolved,
            initial_next_id,
            next_id
        );
    }

    let final_gfa = if changed {
        render_graph(&graph)
    } else {
        original_gfa
            .map(str::to_string)
            .unwrap_or_else(|| render_graph(&graph))
    };
    Ok(ResolvedGfa {
        gfa: final_gfa,
        stats,
    })
}

#[derive(Clone, Debug)]
struct MultiLevelWindowCandidate {
    candidate: BubbleCandidate,
    source: MultiLevelCandidateSource,
    source_sites: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
enum MultiLevelCandidateSource {
    TopLevel,
    SiblingRun,
    ParentDescendants,
    SlidingWindow,
    LevelWindow,
    StringyNeighborhood,
    OutwardResidualWindow,
}

impl MultiLevelCandidateSource {
    fn as_str(self) -> &'static str {
        match self {
            Self::TopLevel => "top-level",
            Self::SiblingRun => "sibling-run",
            Self::ParentDescendants => "parent-descendants",
            Self::SlidingWindow => "sliding-window",
            Self::LevelWindow => "level-window",
            Self::StringyNeighborhood => "stringy-neighborhood",
            Self::OutwardResidualWindow => "outward-residual-window",
        }
    }
}

#[derive(Clone, Debug, Default)]
struct MultiLevelGeneratedCandidates {
    candidates: Vec<MultiLevelWindowCandidate>,
    generated_total: usize,
    source_counts: FxHashMap<MultiLevelCandidateSource, usize>,
    repeat_boundary_rejected: usize,
    total_bp_cap_exceeded: usize,
    max_len_cap_exceeded: usize,
    median_len_cap_exceeded: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct ReplacementObjectiveDelta {
    input_segments: usize,
    input_segment_bp: usize,
    output_segments: usize,
    output_segment_bp: usize,
    input_bp_weighted_coverage: f64,
    output_bp_weighted_coverage: f64,
    bp_weighted_coverage_delta_scaled: i128,
    input_singleton_bp: usize,
    output_singleton_bp: usize,
    singleton_bp_delta: i128,
    segment_delta: i128,
    segment_bp_delta: i128,
    score_delta: i128,
}

#[derive(Clone, Debug)]
struct MultiLevelBuiltCandidate {
    source: MultiLevelCandidateSource,
    source_sites: usize,
    plan: ReplacementPlan,
    objective: ReplacementObjectiveDelta,
    evidence: Option<SweepgaReplacementEvidence>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct MultiLevelBuildEstimate {
    pair_alignments: usize,
    paf_bytes: usize,
    transclosure_cells: usize,
    poasta_cells: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MultiLevelBuildBudgetKind {
    PairAlignments,
    PafBytes,
    TransclosureCells,
    PoastaCells,
}

impl MultiLevelBuildBudgetKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::PairAlignments => "pair-alignments",
            Self::PafBytes => "paf-bytes",
            Self::TransclosureCells => "transclosure-cells",
            Self::PoastaCells => "poasta-cells",
        }
    }
}

#[derive(Clone, Debug)]
struct MultiLevelBuildBudgetDiagnostic {
    window: MultiLevelWindowCandidate,
    method: ResolutionMethod,
    estimate: MultiLevelBuildEstimate,
    kind: MultiLevelBuildBudgetKind,
    observed: usize,
    limit: usize,
}

fn estimate_multi_level_build(window: &MultiLevelWindowCandidate) -> MultiLevelBuildEstimate {
    let stats = window.candidate.traversal_stats;
    let pair_alignments = stats.count.saturating_mul(stats.count.saturating_sub(1));
    let average_len = if stats.count == 0 {
        0
    } else {
        stats
            .total_len
            .saturating_add(stats.count.saturating_sub(1))
            / stats.count
    };
    let estimated_record_bytes =
        128usize.saturating_add(stats.max_len.max(stats.median_len).max(average_len));
    MultiLevelBuildEstimate {
        pair_alignments,
        paf_bytes: pair_alignments.saturating_mul(estimated_record_bytes),
        transclosure_cells: pair_alignments.saturating_mul(stats.median_len.max(average_len)),
        poasta_cells: stats.count.saturating_mul(stats.total_len),
    }
}

fn method_uses_pairwise_induction(method: ResolutionMethod) -> bool {
    matches!(
        method,
        ResolutionMethod::Allwave | ResolutionMethod::Sweepga | ResolutionMethod::Wfmash
    )
}

fn method_uses_poasta_builder(method: ResolutionMethod) -> bool {
    matches!(
        method,
        ResolutionMethod::Poasta | ResolutionMethod::ChainGreedy | ResolutionMethod::ChainPovu
    )
}

fn multi_level_build_budget_diagnostic(
    window: &MultiLevelWindowCandidate,
    config: &ResolutionConfig,
    method: ResolutionMethod,
) -> Option<MultiLevelBuildBudgetDiagnostic> {
    if window.source != MultiLevelCandidateSource::OutwardResidualWindow {
        return None;
    }
    let estimate = estimate_multi_level_build(window);
    let checks = [
        (
            method_uses_pairwise_induction(method),
            config.max_pair_alignments,
            estimate.pair_alignments,
            MultiLevelBuildBudgetKind::PairAlignments,
        ),
        (
            method_uses_pairwise_induction(method),
            config.max_replacement_paf_bytes,
            estimate.paf_bytes,
            MultiLevelBuildBudgetKind::PafBytes,
        ),
        (
            method_uses_pairwise_induction(method),
            config.multi_level_max_transclosure_cells,
            estimate.transclosure_cells,
            MultiLevelBuildBudgetKind::TransclosureCells,
        ),
        (
            method_uses_poasta_builder(method),
            config.multi_level_max_poasta_cells,
            estimate.poasta_cells,
            MultiLevelBuildBudgetKind::PoastaCells,
        ),
    ];
    checks
        .into_iter()
        .filter(|(applies, limit, observed, _)| *applies && *limit > 0 && *observed > *limit)
        .map(
            |(_, limit, observed, kind)| MultiLevelBuildBudgetDiagnostic {
                window: window.clone(),
                method,
                estimate,
                kind,
                observed,
                limit,
            },
        )
        .next()
}

fn multi_level_build_budget_diagnostics(
    candidates: &[MultiLevelWindowCandidate],
    config: &ResolutionConfig,
) -> Vec<MultiLevelBuildBudgetDiagnostic> {
    let mut diagnostics = Vec::new();
    for window in candidates {
        let method = multi_level_window_replacement_method(window, config);
        if let Some(diagnostic) = multi_level_build_budget_diagnostic(window, config, method) {
            diagnostics.push(diagnostic);
        }
    }
    diagnostics
}

fn format_multi_level_budget_diagnostic(
    graph: &Graph,
    diagnostic: &MultiLevelBuildBudgetDiagnostic,
) -> String {
    format!(
        "source={} sites={} method={:?} {}; traversals={}, max-len={}, median-len={}, total-len={}, pair-alignments={} paf-bytes={} transclosure-cells={} poasta-cells={}; budget-exceeded={} observed={} limit={} (diagnostic only)",
        diagnostic.window.source.as_str(),
        diagnostic.window.source_sites,
        diagnostic.method,
        format_candidate_boundary(graph, &diagnostic.window.candidate),
        diagnostic.window.candidate.traversal_stats.count,
        diagnostic.window.candidate.traversal_stats.max_len,
        diagnostic.window.candidate.traversal_stats.median_len,
        diagnostic.window.candidate.traversal_stats.total_len,
        diagnostic.estimate.pair_alignments,
        diagnostic.estimate.paf_bytes,
        diagnostic.estimate.transclosure_cells,
        diagnostic.estimate.poasta_cells,
        diagnostic.kind.as_str(),
        diagnostic.observed,
        diagnostic.limit
    )
}

fn format_multi_level_budget_diagnostics(
    graph: &Graph,
    diagnostics: &[MultiLevelBuildBudgetDiagnostic],
    limit: usize,
) -> String {
    diagnostics
        .iter()
        .take(limit)
        .enumerate()
        .map(|(idx, diagnostic)| {
            format!(
                "#{} {}",
                idx + 1,
                format_multi_level_budget_diagnostic(graph, diagnostic)
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn build_multi_level_window_candidate(
    graph: &Graph,
    config: &ResolutionConfig,
    objective_mode: MultiLevelObjectiveMode,
    emit_logs: bool,
    selected_count: usize,
    round_number: usize,
    build_started: &AtomicUsize,
    build_progress: &AtomicUsize,
    mut window: MultiLevelWindowCandidate,
) -> io::Result<Option<MultiLevelBuiltCandidate>> {
    if !materialize_candidate_sequences(graph, &mut window.candidate, config) {
        return Ok(None);
    }
    let method = multi_level_window_replacement_method(&window, config);
    if emit_logs && selected_count <= 128 {
        let started = build_started.fetch_add(1, Ordering::Relaxed) + 1;
        log::info!(
            "crush iterative-multi-level round {}: building candidate {}/{} source={} sites={} with {:?}; traversals={}, max-len={}, median-len={}, total-len={}, root-span={}",
            round_number,
            started,
            selected_count,
            window.source.as_str(),
            window.source_sites,
            method,
            window.candidate.traversal_stats.count,
            window.candidate.traversal_stats.max_len,
            window.candidate.traversal_stats.median_len,
            window.candidate.traversal_stats.total_len,
            window.candidate.root_span
        );
    }
    let source = window.source;
    let source_sites = window.source_sites;
    let candidate_boundary = if emit_logs && selected_count <= 128 {
        format_candidate_boundary(graph, &window.candidate)
    } else {
        String::new()
    };
    let result = (|| {
        let report = build_replacement_with_method_catching_unwind_report(
            &window.candidate,
            config,
            method,
        )?;
        let replacement = report.replacement;
        if replacement.segments.is_empty() {
            return Ok::<Option<MultiLevelBuiltCandidate>, io::Error>(None);
        }
        let objective =
            replacement_objective_delta(graph, &window.candidate, &replacement, objective_mode);
        Ok::<Option<MultiLevelBuiltCandidate>, io::Error>(Some(MultiLevelBuiltCandidate {
            source,
            source_sites,
            plan: ReplacementPlan {
                candidate: window.candidate,
                replacement,
            },
            objective,
            evidence: report.evidence,
        }))
    })();
    if emit_logs {
        let done = build_progress.fetch_add(1, Ordering::Relaxed) + 1;
        if selected_count <= 128 || done == selected_count || done % 25 == 0 {
            let status = match &result {
                Ok(Some(_)) => "built",
                Ok(None) => "empty",
                Err(_) => "failed",
            };
            log::info!(
                "crush iterative-multi-level round {}: build progress {}/{} ({})",
                round_number,
                done,
                selected_count,
                status
            );
        }
        match &result {
            Ok(Some(candidate)) if selected_count <= 128 => {
                let evidence = candidate
                    .evidence
                    .as_ref()
                    .map(format_sweepga_evidence)
                    .unwrap_or_else(|| "alignment_evidence=not-applicable".to_string());
                log::info!(
                    "crush iterative-multi-level round {}: built candidate detail source={} sites={} method={:?} {}; {}; {}; replacement_segments={}, replacement_bp={}, objective_score_delta={}, coverage_delta_scaled={}, singleton_bp_delta={}, segment_delta={}, segment_bp_delta={}",
                    round_number,
                    candidate.source.as_str(),
                    candidate.source_sites,
                    method,
                    format_candidate_boundary(graph, &candidate.plan.candidate),
                    format_candidate_repeat_summary(graph, &candidate.plan.candidate),
                    evidence,
                    candidate.plan.replacement.segments.len(),
                    replacement_segment_bp(&candidate.plan.replacement),
                    candidate.objective.score_delta,
                    candidate.objective.bp_weighted_coverage_delta_scaled,
                    candidate.objective.singleton_bp_delta,
                    candidate.objective.segment_delta,
                    candidate.objective.segment_bp_delta
                );
            }
            Ok(None) if selected_count <= 128 => {
                log::info!(
                    "crush iterative-multi-level round {}: candidate returned empty replacement source={} sites={} method={:?} {}; reason=empty-replacement",
                    round_number,
                    source.as_str(),
                    source_sites,
                    method,
                    candidate_boundary
                );
            }
            Err(err) if selected_count <= 128 => {
                log::info!(
                    "crush iterative-multi-level round {}: candidate failed source={} sites={} method={:?} {}; reason={}",
                    round_number,
                    source.as_str(),
                    source_sites,
                    method,
                    candidate_boundary,
                    err
                );
            }
            _ => {}
        }
    }
    result
}

fn multi_level_replacement_build_parallelism(
    total_threads: usize,
    candidate_count: usize,
) -> usize {
    total_threads.max(1).min(4).min(candidate_count.max(1))
}

fn multi_level_target_bp(config: &ResolutionConfig) -> usize {
    if config.multi_level_window_target_bp > 0 {
        config.multi_level_window_target_bp
    } else {
        config.chain_greedy_target_bp.max(1)
    }
}

fn multi_level_window_total_bp_cap(config: &ResolutionConfig) -> usize {
    config.max_total_sequence
}

fn effective_multi_level_objective(config: &ResolutionConfig) -> MultiLevelObjectiveMode {
    if config.method == ResolutionMethod::CoverageMultiBubble {
        MultiLevelObjectiveMode::Coverage
    } else {
        config.multi_level_objective
    }
}

fn effective_repeat_aware_boundaries(config: &ResolutionConfig) -> bool {
    config.multi_level_repeat_aware_boundaries
        || config.method == ResolutionMethod::CoverageMultiBubble
}

fn effective_multi_level_window_mode(config: &ResolutionConfig) -> MultiLevelWindowMode {
    config.multi_level_window_mode
}

fn effective_multi_level_candidate_limit(config: &ResolutionConfig) -> usize {
    multi_level_candidate_limit_for_mode(config, effective_multi_level_window_mode(config))
}

fn multi_level_candidate_limit_for_mode(
    config: &ResolutionConfig,
    mode: MultiLevelWindowMode,
) -> usize {
    if matches!(mode, MultiLevelWindowMode::Largest) {
        1
    } else {
        config.multi_level_candidate_limit
    }
}

fn coverage_mode_includes_outward(config: &ResolutionConfig) -> bool {
    config.method == ResolutionMethod::CoverageMultiBubble
        && matches!(
            config.multi_level_window_mode,
            MultiLevelWindowMode::Combined | MultiLevelWindowMode::Sliding
        )
}

fn generate_multi_level_candidates(
    graph: &Graph,
    config: &ResolutionConfig,
    discovered: &[DiscoveredCandidate],
    resolved_signatures: &FxHashSet<String>,
) -> MultiLevelGeneratedCandidates {
    let mut generated = MultiLevelGeneratedCandidates::default();
    if discovered.is_empty() {
        return generated;
    }

    let target_bp = multi_level_target_bp(config);
    let max_sites = config.multi_level_max_window_sites.max(1);
    let path_positions_by_path: Vec<Vec<usize>> = (0..graph.paths.len())
        .map(|idx| path_positions(graph, idx))
        .collect();
    let path_step_indexes: Vec<FxHashMap<Step, Vec<usize>>> =
        graph.paths.iter().map(path_step_index).collect();
    let mut emitted = FxHashSet::<String>::default();
    let node_visits = if effective_repeat_aware_boundaries(config) {
        Some(graph_node_visit_counts(graph))
    } else {
        None
    };
    let window_mode = effective_multi_level_window_mode(config);

    let use_parent = matches!(
        window_mode,
        MultiLevelWindowMode::Largest
            | MultiLevelWindowMode::Parent
            | MultiLevelWindowMode::Combined
    );
    let use_sibling = matches!(
        window_mode,
        MultiLevelWindowMode::Sibling | MultiLevelWindowMode::Combined
    );
    let use_sliding = matches!(
        window_mode,
        MultiLevelWindowMode::Sliding | MultiLevelWindowMode::Combined
    );
    let use_outward = matches!(window_mode, MultiLevelWindowMode::Outward)
        || coverage_mode_includes_outward(config);
    let add_outward_first = coverage_mode_includes_outward(config);

    if use_outward && add_outward_first {
        add_outward_residual_windows(
            graph,
            discovered,
            target_bp,
            max_sites,
            &path_positions_by_path,
            &path_step_indexes,
            node_visits.as_deref(),
            &mut generated,
            &mut emitted,
            resolved_signatures,
            config,
        );
    }

    if use_parent || use_sibling || use_sliding {
        for d in discovered.iter().filter(|d| d.povu_level == 0) {
            insert_multi_level_candidate(
                &mut generated,
                &mut emitted,
                resolved_signatures,
                d.candidate.clone(),
                MultiLevelCandidateSource::TopLevel,
                1,
                config,
            );
        }
    }

    if use_parent || use_sibling {
        add_parent_descendant_windows(
            discovered,
            target_bp,
            &mut generated,
            &mut emitted,
            resolved_signatures,
            config,
        );
    }

    if use_sibling {
        add_same_parent_sibling_windows(
            graph,
            discovered,
            target_bp,
            max_sites,
            &path_positions_by_path,
            &path_step_indexes,
            node_visits.as_deref(),
            &mut generated,
            &mut emitted,
            resolved_signatures,
            config,
        );
    }

    if use_sliding {
        add_sliding_path_order_windows(
            graph,
            discovered,
            target_bp,
            max_sites,
            &path_positions_by_path,
            &path_step_indexes,
            node_visits.as_deref(),
            &mut generated,
            &mut emitted,
            resolved_signatures,
            config,
        );
        add_stringy_neighborhood_windows(
            graph,
            discovered,
            target_bp,
            max_sites,
            &path_positions_by_path,
            &path_step_indexes,
            node_visits.as_deref(),
            &mut generated,
            &mut emitted,
            resolved_signatures,
            config,
        );
    }

    if use_outward && !add_outward_first {
        add_outward_residual_windows(
            graph,
            discovered,
            target_bp,
            max_sites,
            &path_positions_by_path,
            &path_step_indexes,
            node_visits.as_deref(),
            &mut generated,
            &mut emitted,
            resolved_signatures,
            config,
        );
    }

    if matches!(window_mode, MultiLevelWindowMode::Combined) {
        add_cross_level_windows(
            graph,
            discovered,
            target_bp,
            max_sites,
            &path_positions_by_path,
            &path_step_indexes,
            node_visits.as_deref(),
            &mut generated,
            &mut emitted,
            resolved_signatures,
            config,
        );
    }

    generated.generated_total = generated.candidates.len();
    generated.candidates.sort_by(|a, b| {
        multi_level_source_priority(a.source)
            .cmp(&multi_level_source_priority(b.source))
            .then_with(|| compare_same_source_window_priority(a, b, window_mode))
            .then_with(|| {
                b.candidate
                    .estimated_step_savings()
                    .cmp(&a.candidate.estimated_step_savings())
            })
            .then_with(|| b.candidate.unique_steps.cmp(&a.candidate.unique_steps))
            .then_with(|| a.candidate.root_span.cmp(&b.candidate.root_span))
            .then_with(|| {
                a.candidate
                    .root_start_step
                    .cmp(&b.candidate.root_start_step)
            })
            .then_with(|| a.candidate.root_end_step.cmp(&b.candidate.root_end_step))
    });
    let candidate_limit = multi_level_candidate_limit_for_mode(config, window_mode);
    if candidate_limit > 0 && generated.candidates.len() > candidate_limit {
        generated.candidates.truncate(candidate_limit);
    }
    generated
}

fn compare_same_source_window_priority(
    a: &MultiLevelWindowCandidate,
    b: &MultiLevelWindowCandidate,
    window_mode: MultiLevelWindowMode,
) -> std::cmp::Ordering {
    if matches!(
        window_mode,
        MultiLevelWindowMode::Largest
            | MultiLevelWindowMode::Parent
            | MultiLevelWindowMode::Combined
    ) && matches!(
        a.source,
        MultiLevelCandidateSource::TopLevel | MultiLevelCandidateSource::ParentDescendants
    ) && a.source == b.source
    {
        multi_level_parent_frontier_span_score(&b.candidate)
            .cmp(&multi_level_parent_frontier_span_score(&a.candidate))
            .then_with(|| b.candidate.root_span.cmp(&a.candidate.root_span))
            .then_with(|| {
                b.candidate
                    .traversal_stats
                    .total_len
                    .cmp(&a.candidate.traversal_stats.total_len)
            })
            .then_with(|| {
                b.candidate
                    .estimated_step_savings()
                    .cmp(&a.candidate.estimated_step_savings())
            })
    } else if a.source == MultiLevelCandidateSource::OutwardResidualWindow
        && b.source == MultiLevelCandidateSource::OutwardResidualWindow
    {
        b.candidate
            .traversal_stats
            .total_len
            .cmp(&a.candidate.traversal_stats.total_len)
            .then_with(|| b.candidate.root_span.cmp(&a.candidate.root_span))
            .then_with(|| b.source_sites.cmp(&a.source_sites))
    } else {
        std::cmp::Ordering::Equal
    }
}

fn multi_level_parent_frontier_span_score(candidate: &BubbleCandidate) -> usize {
    let traversal_span = candidate
        .traversal_stats
        .median_len
        .max(candidate.traversal_stats.min_len);
    candidate.root_span.min(traversal_span)
}

fn insert_multi_level_candidate(
    generated: &mut MultiLevelGeneratedCandidates,
    emitted: &mut FxHashSet<String>,
    resolved_signatures: &FxHashSet<String>,
    candidate: BubbleCandidate,
    source: MultiLevelCandidateSource,
    source_sites: usize,
    config: &ResolutionConfig,
) {
    if resolved_signatures.contains(&candidate.signature) {
        return;
    }
    let _diagnostic_only = (
        config.max_bubble_span,
        config.min_traversal_len,
        config.max_traversals,
    );
    let total_bp_cap = multi_level_window_total_bp_cap(config);
    if total_bp_cap > 0 && candidate.traversal_stats.total_len > total_bp_cap {
        generated.total_bp_cap_exceeded += 1;
    }
    if config.max_traversal_len > 0 && candidate.traversal_stats.max_len > config.max_traversal_len
    {
        generated.max_len_cap_exceeded += 1;
    }
    if config.max_median_traversal_len > 0
        && candidate.traversal_stats.median_len > config.max_median_traversal_len
    {
        generated.median_len_cap_exceeded += 1;
    }
    if !emitted.insert(candidate.signature.clone()) {
        return;
    }
    *generated.source_counts.entry(source).or_insert(0) += 1;
    generated.candidates.push(MultiLevelWindowCandidate {
        candidate,
        source,
        source_sites,
    });
}

fn add_same_parent_sibling_windows(
    graph: &Graph,
    discovered: &[DiscoveredCandidate],
    target_bp: usize,
    max_sites: usize,
    path_positions_by_path: &[Vec<usize>],
    path_step_indexes: &[FxHashMap<Step, Vec<usize>>],
    node_visits: Option<&[usize]>,
    generated: &mut MultiLevelGeneratedCandidates,
    emitted: &mut FxHashSet<String>,
    resolved_signatures: &FxHashSet<String>,
    config: &ResolutionConfig,
) {
    let mut groups: FxHashMap<(Option<String>, usize), Vec<usize>> = FxHashMap::default();
    for (idx, d) in discovered.iter().enumerate() {
        groups
            .entry((d.povu_parent_id.clone(), d.povu_level))
            .or_default()
            .push(idx);
    }
    for indexes in groups.values_mut() {
        indexes.sort_by(|&a, &b| {
            discovered[a]
                .candidate
                .root_start_step
                .cmp(&discovered[b].candidate.root_start_step)
                .then_with(|| {
                    discovered[a]
                        .candidate
                        .root_end_step
                        .cmp(&discovered[b].candidate.root_end_step)
                })
        });
        add_index_windows(
            graph,
            discovered,
            indexes,
            target_bp,
            max_sites,
            path_positions_by_path,
            path_step_indexes,
            node_visits,
            MultiLevelCandidateSource::SiblingRun,
            generated,
            emitted,
            resolved_signatures,
            config,
        );
    }
}

fn add_parent_descendant_windows(
    discovered: &[DiscoveredCandidate],
    target_bp: usize,
    generated: &mut MultiLevelGeneratedCandidates,
    emitted: &mut FxHashSet<String>,
    resolved_signatures: &FxHashSet<String>,
    config: &ResolutionConfig,
) {
    for parent in discovered {
        let mut descendants = 0usize;
        for child in discovered {
            if child.povu_site_id == parent.povu_site_id {
                continue;
            }
            if parent.candidate.root_start_step <= child.candidate.root_start_step
                && parent.candidate.root_end_step >= child.candidate.root_end_step
                && child.povu_level > parent.povu_level
            {
                descendants += 1;
            }
        }
        if descendants == 0 || chain_povu_block_bp(&parent.candidate) > target_bp {
            continue;
        }
        insert_multi_level_candidate(
            generated,
            emitted,
            resolved_signatures,
            parent.candidate.clone(),
            MultiLevelCandidateSource::ParentDescendants,
            descendants + 1,
            config,
        );
    }
}

fn add_sliding_path_order_windows(
    graph: &Graph,
    discovered: &[DiscoveredCandidate],
    target_bp: usize,
    max_sites: usize,
    path_positions_by_path: &[Vec<usize>],
    path_step_indexes: &[FxHashMap<Step, Vec<usize>>],
    node_visits: Option<&[usize]>,
    generated: &mut MultiLevelGeneratedCandidates,
    emitted: &mut FxHashSet<String>,
    resolved_signatures: &FxHashSet<String>,
    config: &ResolutionConfig,
) {
    let mut indexes = reference_path_order_indexes(discovered);
    indexes = reference_minimal_discovered_indexes(discovered, indexes);
    add_index_windows(
        graph,
        discovered,
        &indexes,
        target_bp,
        max_sites,
        path_positions_by_path,
        path_step_indexes,
        node_visits,
        MultiLevelCandidateSource::SlidingWindow,
        generated,
        emitted,
        resolved_signatures,
        config,
    );
}

fn add_cross_level_windows(
    graph: &Graph,
    discovered: &[DiscoveredCandidate],
    target_bp: usize,
    max_sites: usize,
    path_positions_by_path: &[Vec<usize>],
    path_step_indexes: &[FxHashMap<Step, Vec<usize>>],
    node_visits: Option<&[usize]>,
    generated: &mut MultiLevelGeneratedCandidates,
    emitted: &mut FxHashSet<String>,
    resolved_signatures: &FxHashSet<String>,
    config: &ResolutionConfig,
) {
    let max_level = discovered.iter().map(|d| d.povu_level).max().unwrap_or(0);
    for level_start in 0..=max_level {
        for level_end in level_start..=max_level.min(level_start + 2) {
            let indexes = reference_path_order_indexes(discovered)
                .into_iter()
                .filter(|&idx| {
                    let level = discovered[idx].povu_level;
                    level >= level_start && level <= level_end
                })
                .collect::<Vec<_>>();
            add_index_windows(
                graph,
                discovered,
                &indexes,
                target_bp,
                max_sites,
                path_positions_by_path,
                path_step_indexes,
                node_visits,
                MultiLevelCandidateSource::LevelWindow,
                generated,
                emitted,
                resolved_signatures,
                config,
            );
        }
    }
}

fn add_stringy_neighborhood_windows(
    graph: &Graph,
    discovered: &[DiscoveredCandidate],
    target_bp: usize,
    max_sites: usize,
    path_positions_by_path: &[Vec<usize>],
    path_step_indexes: &[FxHashMap<Step, Vec<usize>>],
    node_visits: Option<&[usize]>,
    generated: &mut MultiLevelGeneratedCandidates,
    emitted: &mut FxHashSet<String>,
    resolved_signatures: &FxHashSet<String>,
    config: &ResolutionConfig,
) {
    let indexes =
        reference_minimal_discovered_indexes(discovered, reference_path_order_indexes(discovered));
    if indexes.is_empty() {
        return;
    }
    let radius = max_sites.saturating_sub(1) / 2;
    for (pos, &idx) in indexes.iter().enumerate() {
        if !is_locally_stringy_candidate(&discovered[idx].candidate) {
            continue;
        }
        let start = pos.saturating_sub(radius);
        let end = (pos + radius + 1).min(indexes.len());
        let window = indexes[start..end].to_vec();
        if window.len() < 2 {
            continue;
        }
        add_specific_index_window(
            graph,
            discovered,
            &window,
            target_bp,
            path_positions_by_path,
            path_step_indexes,
            node_visits,
            MultiLevelCandidateSource::StringyNeighborhood,
            generated,
            emitted,
            resolved_signatures,
            config,
        );
    }
}

fn add_outward_residual_windows(
    graph: &Graph,
    discovered: &[DiscoveredCandidate],
    target_bp: usize,
    max_sites: usize,
    path_positions_by_path: &[Vec<usize>],
    path_step_indexes: &[FxHashMap<Step, Vec<usize>>],
    node_visits: Option<&[usize]>,
    generated: &mut MultiLevelGeneratedCandidates,
    emitted: &mut FxHashSet<String>,
    resolved_signatures: &FxHashSet<String>,
    config: &ResolutionConfig,
) {
    let indexes =
        reference_minimal_discovered_indexes(discovered, reference_path_order_indexes(discovered));

    let path_count = graph.paths.len().max(1);
    for site in discovered {
        if site.candidate.root_span <= target_bp
            && is_outward_residual_seed(&site.candidate, path_count)
        {
            insert_multi_level_candidate(
                generated,
                emitted,
                resolved_signatures,
                site.candidate.clone(),
                MultiLevelCandidateSource::OutwardResidualWindow,
                1,
                config,
            );
        }
    }

    if indexes.len() < 2 {
        return;
    }

    for (pos, &seed_idx) in indexes.iter().enumerate() {
        if !is_outward_residual_seed(&discovered[seed_idx].candidate, path_count) {
            continue;
        }

        for window_size in 2..=max_sites.min(indexes.len()) {
            let min_start = pos.saturating_add(1).saturating_sub(window_size);
            let max_start = pos.min(indexes.len() - window_size);
            if min_start > max_start {
                continue;
            }
            let preferred = pos
                .saturating_sub((window_size - 1) / 2)
                .max(min_start)
                .min(max_start);
            let mut starts = [preferred, min_start, max_start]
                .into_iter()
                .collect::<Vec<_>>();
            starts.sort_unstable();
            starts.dedup();
            for start in starts {
                let end = start + window_size;
                add_specific_index_window(
                    graph,
                    discovered,
                    &indexes[start..end],
                    target_bp,
                    path_positions_by_path,
                    path_step_indexes,
                    node_visits,
                    MultiLevelCandidateSource::OutwardResidualWindow,
                    generated,
                    emitted,
                    resolved_signatures,
                    config,
                );
            }
        }
    }
}

fn add_index_windows(
    graph: &Graph,
    discovered: &[DiscoveredCandidate],
    indexes: &[usize],
    target_bp: usize,
    max_sites: usize,
    path_positions_by_path: &[Vec<usize>],
    path_step_indexes: &[FxHashMap<Step, Vec<usize>>],
    node_visits: Option<&[usize]>,
    source: MultiLevelCandidateSource,
    generated: &mut MultiLevelGeneratedCandidates,
    emitted: &mut FxHashSet<String>,
    resolved_signatures: &FxHashSet<String>,
    config: &ResolutionConfig,
) {
    if indexes.len() < 2 {
        return;
    }
    for start in 0..indexes.len() {
        let max_end = (start + max_sites).min(indexes.len());
        for end in start + 2..=max_end {
            add_specific_index_window(
                graph,
                discovered,
                &indexes[start..end],
                target_bp,
                path_positions_by_path,
                path_step_indexes,
                node_visits,
                source,
                generated,
                emitted,
                resolved_signatures,
                config,
            );
        }
    }
}

fn add_specific_index_window(
    graph: &Graph,
    discovered: &[DiscoveredCandidate],
    indexes: &[usize],
    target_bp: usize,
    path_positions_by_path: &[Vec<usize>],
    path_step_indexes: &[FxHashMap<Step, Vec<usize>>],
    node_visits: Option<&[usize]>,
    source: MultiLevelCandidateSource,
    generated: &mut MultiLevelGeneratedCandidates,
    emitted: &mut FxHashSet<String>,
    resolved_signatures: &FxHashSet<String>,
    config: &ResolutionConfig,
) {
    let Some((&first_idx, rest)) = indexes.split_first() else {
        return;
    };
    let mut root_start = discovered[first_idx].candidate.root_start_step;
    let mut root_end = discovered[first_idx].candidate.root_end_step;
    let mut level = discovered[first_idx].povu_level;
    for &idx in rest {
        root_start = root_start.min(discovered[idx].candidate.root_start_step);
        root_end = root_end.max(discovered[idx].candidate.root_end_step);
        level = level.min(discovered[idx].povu_level);
    }
    let root_positions = &path_positions_by_path[0];
    let Some(end_bp) = root_end
        .checked_add(1)
        .and_then(|idx| root_positions.get(idx).copied())
    else {
        return;
    };
    let Some(start_bp) = root_positions.get(root_start).copied() else {
        return;
    };
    if end_bp.saturating_sub(start_bp) > target_bp {
        return;
    }
    if let Some(visits) = node_visits {
        let Some(root_path) = graph.paths.first() else {
            return;
        };
        let Some(&entry) = root_path.steps.get(root_start) else {
            return;
        };
        let Some(&exit) = root_path.steps.get(root_end) else {
            return;
        };
        if repeat_boundary_should_reject(graph, entry, exit, visits) {
            generated.repeat_boundary_rejected += 1;
        }
    }
    let Some(candidate) = candidate_from_root_interval(
        graph,
        root_start,
        root_end,
        level,
        path_positions_by_path,
        path_step_indexes,
    ) else {
        return;
    };
    insert_multi_level_candidate(
        generated,
        emitted,
        resolved_signatures,
        candidate,
        source,
        indexes.len(),
        config,
    );
}

fn candidate_from_root_interval(
    graph: &Graph,
    root_start_step: usize,
    root_end_step: usize,
    level: usize,
    path_positions_by_path: &[Vec<usize>],
    path_step_indexes: &[FxHashMap<Step, Vec<usize>>],
) -> Option<BubbleCandidate> {
    let root_path_idx = 0usize;
    let root_path = graph.paths.get(root_path_idx)?;
    let entry = *root_path.steps.get(root_start_step)?;
    let exit = *root_path.steps.get(root_end_step)?;

    let mut ranges = Vec::new();
    for (path_idx, path_index) in path_step_indexes.iter().enumerate() {
        if let Some((begin_step, end_step)) = unique_anchor_range(path_index, entry, exit) {
            ranges.push(PathRange {
                path_idx,
                begin_step,
                end_step,
                sequence: Vec::new(),
                extended_sequence: Vec::new(),
                left_flank_bp: 0,
                right_flank_bp: 0,
            });
        }
    }
    if ranges.len() < 2 {
        return None;
    }

    let lengths = ranges
        .iter()
        .map(|range| {
            let positions = &path_positions_by_path[range.path_idx];
            positions[range.end_step] - positions[range.begin_step]
        })
        .collect::<Vec<_>>();
    let traversal_stats = traversal_stats_from_lengths(lengths);
    let total_steps = ranges
        .iter()
        .map(|range| range.end_step.saturating_sub(range.begin_step))
        .sum::<usize>();
    let mut unique_steps = FxHashSet::default();
    for range in &ranges {
        let path = &graph.paths[range.path_idx];
        for step_idx in range.begin_step..range.end_step {
            unique_steps.insert(path.steps[step_idx]);
        }
    }

    let root_positions = &path_positions_by_path[root_path_idx];
    let root_span = root_positions[root_end_step + 1] - root_positions[root_start_step];
    let signature = candidate_signature(
        graph,
        root_path_idx,
        root_positions,
        root_start_step,
        root_end_step,
        &ranges,
        path_positions_by_path,
    );

    Some(BubbleCandidate {
        ranges,
        signature,
        root_start_step,
        root_end_step,
        root_span,
        total_steps,
        unique_steps: unique_steps.len(),
        traversal_stats,
        level,
    })
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct RepeatAnchorDiagnostic {
    visits: usize,
    len: usize,
    low_complexity: bool,
    repeat_like: bool,
}

fn graph_node_visit_counts(graph: &Graph) -> Vec<usize> {
    let mut node_visits = vec![0usize; graph.segments.len()];
    for path in &graph.paths {
        for step in &path.steps {
            if step.node < node_visits.len() {
                node_visits[step.node] += 1;
            }
        }
    }
    node_visits
}

fn repeat_anchor_diagnostic(
    graph: &Graph,
    step: Step,
    node_visits: &[usize],
) -> RepeatAnchorDiagnostic {
    let visits = node_visits.get(step.node).copied().unwrap_or(0);
    let len = graph
        .segments
        .get(step.node)
        .map(|segment| segment.seq.len())
        .unwrap_or(0);
    let low_complexity = graph
        .segments
        .get(step.node)
        .map(|segment| is_low_complexity_dna(&segment.seq))
        .unwrap_or(false);
    let high_frequency = visits >= graph.paths.len().div_ceil(2).max(2);
    let repeat_like = high_frequency && (len <= 64 || low_complexity);
    RepeatAnchorDiagnostic {
        visits,
        len,
        low_complexity,
        repeat_like,
    }
}

fn repeat_boundary_should_reject(
    graph: &Graph,
    entry: Step,
    exit: Step,
    node_visits: &[usize],
) -> bool {
    let entry_diag = repeat_anchor_diagnostic(graph, entry, node_visits);
    let exit_diag = repeat_anchor_diagnostic(graph, exit, node_visits);
    entry_diag.repeat_like && exit_diag.repeat_like
}

fn is_low_complexity_dna(seq: &[u8]) -> bool {
    if seq.len() < 8 {
        return false;
    }
    let mut counts = [0usize; 5];
    for &base in seq {
        match base.to_ascii_uppercase() {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' | b'U' => counts[3] += 1,
            _ => counts[4] += 1,
        }
    }
    let max_single = counts.iter().copied().max().unwrap_or(0);
    if max_single.saturating_mul(100) >= seq.len().saturating_mul(80) {
        return true;
    }
    let mut dinuc_counts = FxHashMap::<[u8; 2], usize>::default();
    for pair in seq.windows(2) {
        let key = [pair[0].to_ascii_uppercase(), pair[1].to_ascii_uppercase()];
        *dinuc_counts.entry(key).or_insert(0) += 1;
    }
    dinuc_counts
        .values()
        .copied()
        .max()
        .unwrap_or(0)
        .saturating_mul(100)
        >= seq.len().saturating_sub(1).saturating_mul(70)
}

fn reference_path_order_indexes(discovered: &[DiscoveredCandidate]) -> Vec<usize> {
    let mut indexes = (0..discovered.len()).collect::<Vec<_>>();
    indexes.sort_by(|&a, &b| {
        discovered[a]
            .candidate
            .root_start_step
            .cmp(&discovered[b].candidate.root_start_step)
            .then_with(|| {
                discovered[a]
                    .candidate
                    .root_end_step
                    .cmp(&discovered[b].candidate.root_end_step)
            })
            .then_with(|| {
                discovered[a]
                    .candidate
                    .signature
                    .cmp(&discovered[b].candidate.signature)
            })
    });
    indexes
}

fn reference_minimal_discovered_indexes(
    discovered: &[DiscoveredCandidate],
    indexes: Vec<usize>,
) -> Vec<usize> {
    indexes
        .iter()
        .copied()
        .filter(|&idx| {
            !indexes.iter().copied().any(|other_idx| {
                idx != other_idx
                    && discovered[other_idx].candidate.root_start_step
                        <= discovered[idx].candidate.root_start_step
                    && discovered[other_idx].candidate.root_end_step
                        >= discovered[idx].candidate.root_end_step
                    && (discovered[other_idx].candidate.root_start_step
                        < discovered[idx].candidate.root_start_step
                        || discovered[other_idx].candidate.root_end_step
                            > discovered[idx].candidate.root_end_step)
            })
        })
        .collect()
}

fn is_locally_stringy_candidate(candidate: &BubbleCandidate) -> bool {
    candidate.traversal_stats.median_len <= 500
        && candidate.traversal_stats.count >= 2
        && candidate.unique_steps >= candidate.traversal_stats.count.saturating_mul(5)
}

fn is_outward_residual_seed(candidate: &BubbleCandidate, path_count: usize) -> bool {
    if is_locally_stringy_candidate(candidate) {
        return true;
    }

    let broad_support = candidate.traversal_stats.count >= path_count.saturating_div(4).max(2);
    let high_step_density =
        candidate.unique_steps >= candidate.traversal_stats.count.saturating_mul(3).max(12);
    let high_total_bp = candidate.traversal_stats.total_len >= 50_000 && broad_support;
    let long_local_region = candidate.root_span >= 1_000 && high_step_density;
    high_total_bp || long_local_region
}

fn multi_level_source_priority(source: MultiLevelCandidateSource) -> u8 {
    match source {
        MultiLevelCandidateSource::TopLevel => 0,
        MultiLevelCandidateSource::ParentDescendants => 1,
        MultiLevelCandidateSource::OutwardResidualWindow => 2,
        MultiLevelCandidateSource::LevelWindow => 3,
        MultiLevelCandidateSource::SiblingRun => 4,
        MultiLevelCandidateSource::SlidingWindow => 5,
        MultiLevelCandidateSource::StringyNeighborhood => 6,
    }
}

fn multi_level_window_replacement_method(
    window: &MultiLevelWindowCandidate,
    config: &ResolutionConfig,
) -> ResolutionMethod {
    if window.source_sites > 1 {
        match candidate_replacement_method(&window.candidate, config) {
            // Multi-site windows should not go through the tiny-bubble SPOA
            // path: POASTA is the direct aligner tier that can handle these
            // larger local structures without invoking pairwise graph induction.
            ResolutionMethod::Poa => ResolutionMethod::Poasta,
            method => method,
        }
    } else {
        candidate_replacement_method(&window.candidate, config)
    }
}

fn replacement_objective_delta(
    graph: &Graph,
    candidate: &BubbleCandidate,
    replacement: &Graph,
    objective_mode: MultiLevelObjectiveMode,
) -> ReplacementObjectiveDelta {
    let (input_segments, input_segment_bp) = candidate_unique_segment_stats(graph, candidate);
    let output_segments = replacement.segments.len();
    let output_segment_bp = replacement_segment_bp(replacement);
    let input_coverage = candidate_coverage_quality(graph, candidate);
    let output_quality = graph_quality(replacement);
    let input_score = size_objective(input_segments, input_segment_bp, candidate.total_steps);
    let output_path_steps = replacement
        .paths
        .iter()
        .map(|path| path.steps.len())
        .sum::<usize>();
    let output_score = size_objective(output_segments, output_segment_bp, output_path_steps);
    let bp_weighted_coverage_delta_scaled = scaled_coverage_delta(
        input_coverage.bp_weighted_mean,
        output_quality.node_coverage_bp_weighted_mean,
    );
    let singleton_bp_delta =
        input_coverage.singleton_bp as i128 - output_quality.singleton_bp as i128;
    let size_delta = input_score - output_score;
    let score_delta = match objective_mode {
        MultiLevelObjectiveMode::Size => size_delta,
        MultiLevelObjectiveMode::Coverage => coverage_acceptance_delta(
            input_coverage.bp_weighted_mean,
            input_coverage.singleton_bp,
            output_quality.node_coverage_bp_weighted_mean,
            output_quality.singleton_bp,
        ),
    };
    ReplacementObjectiveDelta {
        input_segments,
        input_segment_bp,
        output_segments,
        output_segment_bp,
        input_bp_weighted_coverage: input_coverage.bp_weighted_mean,
        output_bp_weighted_coverage: output_quality.node_coverage_bp_weighted_mean,
        bp_weighted_coverage_delta_scaled,
        input_singleton_bp: input_coverage.singleton_bp,
        output_singleton_bp: output_quality.singleton_bp,
        singleton_bp_delta,
        segment_delta: input_segments as i128 - output_segments as i128,
        segment_bp_delta: input_segment_bp as i128 - output_segment_bp as i128,
        score_delta,
    }
}

fn candidate_unique_segment_stats(graph: &Graph, candidate: &BubbleCandidate) -> (usize, usize) {
    let mut nodes = FxHashSet::<usize>::default();
    for range in &candidate.ranges {
        if let Some(path) = graph.paths.get(range.path_idx) {
            for step_idx in range.begin_step..range.end_step {
                if let Some(step) = path.steps.get(step_idx) {
                    nodes.insert(step.node);
                }
            }
        }
    }
    let bp = nodes
        .iter()
        .map(|&node| {
            graph
                .segments
                .get(node)
                .map(|segment| segment.seq.len())
                .unwrap_or(0)
        })
        .sum();
    (nodes.len(), bp)
}

#[derive(Clone, Copy, Debug, Default)]
struct CoverageQuality {
    bp_weighted_mean: f64,
    singleton_bp: usize,
}

fn candidate_coverage_quality(graph: &Graph, candidate: &BubbleCandidate) -> CoverageQuality {
    let mut visits = FxHashMap::<usize, usize>::default();
    for range in &candidate.ranges {
        if let Some(path) = graph.paths.get(range.path_idx) {
            for step_idx in range.begin_step..range.end_step {
                if let Some(step) = path.steps.get(step_idx) {
                    *visits.entry(step.node).or_insert(0) += 1;
                }
            }
        }
    }
    let mut segment_bp = 0usize;
    let mut weighted_visit_bp = 0u128;
    let mut singleton_bp = 0usize;
    for (node, count) in visits {
        let len = graph
            .segments
            .get(node)
            .map(|segment| segment.seq.len())
            .unwrap_or(0);
        segment_bp = segment_bp.saturating_add(len);
        weighted_visit_bp =
            weighted_visit_bp.saturating_add((count as u128).saturating_mul(len as u128));
        if count == 1 {
            singleton_bp = singleton_bp.saturating_add(len);
        }
    }
    let bp_weighted_mean = if segment_bp == 0 {
        0.0
    } else {
        weighted_visit_bp as f64 / segment_bp as f64
    };
    CoverageQuality {
        bp_weighted_mean,
        singleton_bp,
    }
}

const COVERAGE_OBJECTIVE_SCALE: f64 = 1_000_000.0;
const COVERAGE_SINGLETON_BP_WEIGHT: i128 = 1024;

fn scaled_coverage_delta(before_bp_weighted: f64, after_bp_weighted: f64) -> i128 {
    ((after_bp_weighted - before_bp_weighted) * COVERAGE_OBJECTIVE_SCALE).round() as i128
}

fn coverage_acceptance_delta(
    before_bp_weighted: f64,
    before_singleton_bp: usize,
    after_bp_weighted: f64,
    after_singleton_bp: usize,
) -> i128 {
    let coverage_delta = scaled_coverage_delta(before_bp_weighted, after_bp_weighted);
    let singleton_delta = before_singleton_bp as i128 - after_singleton_bp as i128;
    if coverage_delta > 0 || singleton_delta > 0 {
        coverage_delta.max(0).saturating_add(
            singleton_delta
                .max(0)
                .saturating_mul(COVERAGE_SINGLETON_BP_WEIGHT),
        )
    } else {
        coverage_delta.saturating_add(singleton_delta.saturating_mul(COVERAGE_SINGLETON_BP_WEIGHT))
    }
}

fn graph_objective_delta_for_mode(
    before: GraphQuality,
    after: GraphQuality,
    objective_mode: MultiLevelObjectiveMode,
) -> i128 {
    match objective_mode {
        MultiLevelObjectiveMode::Size => {
            size_objective(before.segments, before.segment_bp, before.path_steps)
                - size_objective(after.segments, after.segment_bp, after.path_steps)
        }
        MultiLevelObjectiveMode::Coverage => coverage_acceptance_delta(
            before.node_coverage_bp_weighted_mean,
            before.singleton_bp,
            after.node_coverage_bp_weighted_mean,
            after.singleton_bp,
        ),
    }
}

fn size_objective(segments: usize, segment_bp: usize, path_steps: usize) -> i128 {
    (segment_bp as i128)
        .saturating_mul(1024)
        .saturating_add((segments as i128).saturating_mul(16))
        .saturating_add(path_steps as i128)
}

fn format_multi_level_source_counts(
    counts: &FxHashMap<MultiLevelCandidateSource, usize>,
) -> String {
    let ordered = [
        MultiLevelCandidateSource::TopLevel,
        MultiLevelCandidateSource::SiblingRun,
        MultiLevelCandidateSource::ParentDescendants,
        MultiLevelCandidateSource::SlidingWindow,
        MultiLevelCandidateSource::LevelWindow,
        MultiLevelCandidateSource::StringyNeighborhood,
        MultiLevelCandidateSource::OutwardResidualWindow,
    ];
    ordered
        .iter()
        .filter_map(|source| {
            counts
                .get(source)
                .copied()
                .filter(|&count| count > 0)
                .map(|count| format!("{}={}", source.as_str(), count))
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn multi_level_window_source_counts(
    candidates: &[MultiLevelWindowCandidate],
) -> FxHashMap<MultiLevelCandidateSource, usize> {
    candidates.iter().fold(
        FxHashMap::<MultiLevelCandidateSource, usize>::default(),
        |mut acc, window| {
            *acc.entry(window.source).or_insert(0) += 1;
            acc
        },
    )
}

fn format_multi_level_candidate_summary(
    label: &str,
    candidates: &[MultiLevelWindowCandidate],
) -> String {
    if candidates.is_empty() {
        return format!("{label} none");
    }
    let plain = candidates
        .iter()
        .map(|candidate| candidate.candidate.clone())
        .collect::<Vec<_>>();
    let source_sites = candidates
        .iter()
        .map(|candidate| candidate.source_sites)
        .collect::<Vec<_>>();
    let root_spans = candidates
        .iter()
        .map(|candidate| candidate.candidate.root_span)
        .collect::<Vec<_>>();
    let region_bp = candidates
        .iter()
        .map(|candidate| candidate.candidate.traversal_stats.total_len)
        .collect::<Vec<_>>();
    format!(
        "{}; window-sites {}; root-span-bp {}; region-bp {}",
        format_candidate_length_summary(label, &plain),
        format_bp_distribution("sites", &source_sites),
        format_bp_distribution("bp", &root_spans),
        format_bp_distribution("bp", &region_bp)
    )
}

fn format_step_label(graph: &Graph, step: Step) -> String {
    let direction = if step.rev { '<' } else { '>' };
    let id = graph
        .segments
        .get(step.node)
        .map(|segment| segment.id.as_str())
        .unwrap_or("?");
    format!("{direction}{id}")
}

fn format_candidate_boundary(graph: &Graph, candidate: &BubbleCandidate) -> String {
    let root_path = graph.paths.first();
    let start = root_path
        .and_then(|path| path.steps.get(candidate.root_start_step).copied())
        .map(|step| format_step_label(graph, step))
        .unwrap_or_else(|| "?".to_string());
    let end = root_path
        .and_then(|path| path.steps.get(candidate.root_end_step).copied())
        .map(|step| format_step_label(graph, step))
        .unwrap_or_else(|| "?".to_string());
    format!(
        "{}..{} root_steps={}..{} root_span={}bp",
        start, end, candidate.root_start_step, candidate.root_end_step, candidate.root_span
    )
}

fn format_candidate_repeat_summary(graph: &Graph, candidate: &BubbleCandidate) -> String {
    let node_visits = graph_node_visit_counts(graph);
    let mut nodes = FxHashSet::<usize>::default();
    for range in &candidate.ranges {
        if let Some(path) = graph.paths.get(range.path_idx) {
            for step_idx in range.begin_step..range.end_step {
                if let Some(step) = path.steps.get(step_idx) {
                    nodes.insert(step.node);
                }
            }
        }
    }
    let mut repeat_like_nodes = 0usize;
    let mut repeat_like_bp = 0usize;
    let mut low_complexity_nodes = 0usize;
    let mut low_complexity_bp = 0usize;
    for &node in &nodes {
        let step = Step { node, rev: false };
        let diag = repeat_anchor_diagnostic(graph, step, &node_visits);
        if diag.low_complexity {
            low_complexity_nodes += 1;
            low_complexity_bp = low_complexity_bp.saturating_add(diag.len);
        }
        if diag.repeat_like {
            repeat_like_nodes += 1;
            repeat_like_bp = repeat_like_bp.saturating_add(diag.len);
        }
    }

    let mut seq_counts = FxHashMap::<Vec<u8>, usize>::default();
    let mut low_complexity_traversals = 0usize;
    for range in &candidate.ranges {
        if is_low_complexity_dna(&range.sequence) {
            low_complexity_traversals += 1;
        }
        *seq_counts.entry(range.sequence.clone()).or_insert(0) += 1;
    }
    let duplicate_sequence_classes = seq_counts.values().filter(|&&count| count > 1).count();
    let top_duplicate_count = seq_counts.values().copied().max().unwrap_or(0);
    format!(
        "repeat_diag=unique_nodes={},repeat_like_nodes={}({}bp),low_complexity_nodes={}({}bp),duplicate_traversal_classes={},top_duplicate_count={},low_complexity_traversals={}",
        nodes.len(),
        repeat_like_nodes,
        repeat_like_bp,
        low_complexity_nodes,
        low_complexity_bp,
        duplicate_sequence_classes,
        top_duplicate_count,
        low_complexity_traversals
    )
}

fn format_multi_level_candidate_details(
    graph: &Graph,
    candidates: &[MultiLevelWindowCandidate],
    limit: usize,
) -> String {
    candidates
        .iter()
        .take(limit)
        .enumerate()
        .map(|(idx, window)| {
            format!(
                "#{} source={} sites={} level={} {} traversals={} max={} median={} total={} unique_steps={} step_savings={}",
                idx + 1,
                window.source.as_str(),
                window.source_sites,
                window.candidate.level,
                format_candidate_boundary(graph, &window.candidate),
                window.candidate.traversal_stats.count,
                window.candidate.traversal_stats.max_len,
                window.candidate.traversal_stats.median_len,
                window.candidate.traversal_stats.total_len,
                window.candidate.unique_steps,
                window.candidate.estimated_step_savings()
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn format_discovered_residual_sites(
    graph: &Graph,
    discovered: &[DiscoveredCandidate],
    limit: usize,
) -> String {
    let mut sites = discovered.iter().collect::<Vec<_>>();
    sites.sort_by(|a, b| {
        b.candidate
            .root_span
            .cmp(&a.candidate.root_span)
            .then_with(|| {
                b.candidate
                    .traversal_stats
                    .total_len
                    .cmp(&a.candidate.traversal_stats.total_len)
            })
    });
    sites
        .into_iter()
        .take(limit)
        .map(|site| {
            format!(
                "site={} parent={} povu_level={} leaf={} {} traversals={} max={} median={} total={} unique_steps={}",
                site.povu_site_id,
                site.povu_parent_id.as_deref().unwrap_or("."),
                site.povu_level,
                site.is_leaf,
                format_candidate_boundary(graph, &site.candidate),
                site.candidate.traversal_stats.count,
                site.candidate.traversal_stats.max_len,
                site.candidate.traversal_stats.median_len,
                site.candidate.traversal_stats.total_len,
                site.candidate.unique_steps
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn format_objective_delta_summary(candidates: &[MultiLevelBuiltCandidate]) -> String {
    if candidates.is_empty() {
        return "n=0".to_string();
    }
    let score_deltas = candidates
        .iter()
        .map(|candidate| candidate.objective.score_delta)
        .collect::<Vec<_>>();
    let segment_deltas = candidates
        .iter()
        .map(|candidate| candidate.objective.segment_delta)
        .collect::<Vec<_>>();
    let bp_deltas = candidates
        .iter()
        .map(|candidate| candidate.objective.segment_bp_delta)
        .collect::<Vec<_>>();
    let coverage_deltas = candidates
        .iter()
        .map(|candidate| candidate.objective.bp_weighted_coverage_delta_scaled)
        .collect::<Vec<_>>();
    let singleton_deltas = candidates
        .iter()
        .map(|candidate| candidate.objective.singleton_bp_delta)
        .collect::<Vec<_>>();
    let input_segments = candidates
        .iter()
        .map(|candidate| candidate.objective.input_segments)
        .sum::<usize>();
    let output_segments = candidates
        .iter()
        .map(|candidate| candidate.objective.output_segments)
        .sum::<usize>();
    let input_bp = candidates
        .iter()
        .map(|candidate| candidate.objective.input_segment_bp)
        .sum::<usize>();
    let output_bp = candidates
        .iter()
        .map(|candidate| candidate.objective.output_segment_bp)
        .sum::<usize>();
    let source_sites = candidates
        .iter()
        .map(|candidate| candidate.source_sites)
        .collect::<Vec<_>>();
    let input_singleton_bp = candidates
        .iter()
        .map(|candidate| candidate.objective.input_singleton_bp)
        .sum::<usize>();
    let output_singleton_bp = candidates
        .iter()
        .map(|candidate| candidate.objective.output_singleton_bp)
        .sum::<usize>();
    let mean_input_bp_weighted_coverage = candidates
        .iter()
        .map(|candidate| candidate.objective.input_bp_weighted_coverage)
        .sum::<f64>()
        / candidates.len() as f64;
    let mean_output_bp_weighted_coverage = candidates
        .iter()
        .map(|candidate| candidate.objective.output_bp_weighted_coverage)
        .sum::<f64>()
        / candidates.len() as f64;
    format!(
        "n={}, source-sites {}, score-delta {}, coverage-delta-scaled {}, singleton-bp-delta {}, segment-delta {}, segment-bp-delta {}, mean_bp_weighted_coverage={:.4}->{:.4}, input_segments={}, output_segments={}, input_bp={}, output_bp={}, input_singleton_bp={}, output_singleton_bp={}",
        candidates.len(),
        format_bp_distribution("sites", &source_sites),
        format_i128_distribution(&score_deltas),
        format_i128_distribution(&coverage_deltas),
        format_i128_distribution(&singleton_deltas),
        format_i128_distribution(&segment_deltas),
        format_i128_distribution(&bp_deltas),
        mean_input_bp_weighted_coverage,
        mean_output_bp_weighted_coverage,
        input_segments,
        output_segments,
        input_bp,
        output_bp,
        input_singleton_bp,
        output_singleton_bp
    )
}

fn format_i128_distribution(values: &[i128]) -> String {
    if values.is_empty() {
        return "n=0".to_string();
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let sum = sorted.iter().copied().sum::<i128>();
    format!(
        "min={}, p50={}, p90={}, max={}, sum={}",
        sorted[0],
        sorted[sorted.len() / 2],
        sorted[percentile_index(sorted.len(), 90, 100)],
        sorted[sorted.len() - 1],
        sum
    )
}

#[derive(Clone, Debug, Default)]
struct TopFlubbleSweepgaSelection {
    selected: Vec<BubbleCandidate>,
    top_level_regions: usize,
    descendant_sites: usize,
}

fn select_top_flubble_sweepga_regions(
    discovered: &[DiscoveredCandidate],
    resolved_signatures: &FxHashSet<String>,
) -> TopFlubbleSweepgaSelection {
    let top_level_regions = discovered.iter().filter(|d| d.povu_level == 0).count();
    let descendant_sites = discovered.len().saturating_sub(top_level_regions);
    let selected = discovered
        .iter()
        .filter(|d| d.povu_level == 0)
        .filter(|d| !resolved_signatures.contains(&d.candidate.signature))
        .map(|d| {
            let mut candidate = d.candidate.clone();
            candidate.level = 0;
            candidate
        })
        .collect::<Vec<_>>();
    TopFlubbleSweepgaSelection {
        selected,
        top_level_regions,
        descendant_sites,
    }
}

#[derive(Clone, Debug)]
struct TopFlubbleSweepgaBuild {
    region_ordinal: usize,
    candidate: BubbleCandidate,
    replacement: Graph,
    evidence: SweepgaReplacementEvidence,
}

#[derive(Clone, Debug)]
struct TopFlubbleSweepgaEvidenceRow {
    region_ordinal: usize,
    input_traversals: usize,
    input_bp: usize,
    alignment_records: usize,
    aligned_bp: u64,
    paf_bytes: usize,
    collinear_records: usize,
    off_diagonal_records: usize,
    dust_like_records: usize,
    kmer_frequency: usize,
    replacement_segments: usize,
    replacement_shared_segments: usize,
    replacement_bp: usize,
}

impl TopFlubbleSweepgaEvidenceRow {
    fn from_build(build: &TopFlubbleSweepgaBuild) -> Self {
        Self {
            region_ordinal: build.region_ordinal,
            input_traversals: build.candidate.traversal_stats.count,
            input_bp: build.candidate.traversal_stats.total_len,
            alignment_records: build.evidence.alignment_records,
            aligned_bp: build.evidence.aligned_bp,
            paf_bytes: build.evidence.paf_bytes,
            collinear_records: build.evidence.collinear_records,
            off_diagonal_records: build.evidence.off_diagonal_records,
            dust_like_records: build.evidence.dust_like_records,
            kmer_frequency: build.evidence.kmer_frequency,
            replacement_segments: build.evidence.replacement_segments,
            replacement_shared_segments: build.evidence.replacement_shared_segments,
            replacement_bp: build.evidence.replacement_bp,
        }
    }
}

/// Return `(path_name, sequence)` for each `P` line, using blunt concatenation.
pub fn path_sequences(gfa: &str) -> io::Result<Vec<(String, String)>> {
    let graph = parse_gfa(gfa)?;
    graph
        .paths
        .iter()
        .map(|path| {
            let seq = path_sequence(&graph, path)?;
            String::from_utf8(seq)
                .map(|seq| (path.name.clone(), seq))
                .map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("path '{}' is not valid UTF-8 DNA: {}", path.name, err),
                    )
                })
        })
        .collect()
}

fn parse_gfa(gfa: &str) -> io::Result<Graph> {
    let mut segments = Vec::new();
    let mut id_to_idx: FxHashMap<String, usize> = FxHashMap::default();
    let mut paths = Vec::new();

    for (line_no, line) in gfa.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[0] {
            "S" => {
                if fields.len() < 3 {
                    return Err(invalid_gfa(line_no, "S line has fewer than 3 fields"));
                }
                let id = fields[1].to_string();
                if id_to_idx.contains_key(&id) {
                    return Err(invalid_gfa(line_no, "duplicate segment ID"));
                }
                let idx = segments.len();
                id_to_idx.insert(id.clone(), idx);
                segments.push(Segment {
                    id,
                    seq: fields[2].as_bytes().to_vec(),
                });
            }
            "L" => {
                if fields.len() >= 6 && fields[5] != "0M" {
                    return Err(invalid_gfa(
                        line_no,
                        "resolution currently requires blunt 0M links",
                    ));
                }
            }
            "P" => {
                if fields.len() < 3 {
                    return Err(invalid_gfa(line_no, "P line has fewer than 3 fields"));
                }
                if fields.len() >= 4 && fields[3] != "*" && fields[3] != "" {
                    return Err(invalid_gfa(
                        line_no,
                        "resolution currently requires '*' path overlaps",
                    ));
                }
                let mut steps = Vec::new();
                for raw_step in fields[2].split(',').filter(|s| !s.is_empty()) {
                    let (id, rev) = parse_p_step(raw_step)
                        .ok_or_else(|| invalid_gfa(line_no, "P line has a malformed path step"))?;
                    let Some(&node) = id_to_idx.get(id) else {
                        return Err(invalid_gfa(line_no, "P line references unknown segment"));
                    };
                    steps.push(Step { node, rev });
                }
                if !steps.is_empty() {
                    paths.push(Path {
                        name: fields[1].to_string(),
                        steps,
                    });
                }
            }
            "W" => {
                if fields.len() < 7 {
                    return Err(invalid_gfa(line_no, "W line has fewer than 7 fields"));
                }
                let mut steps = Vec::new();
                for (id, rev) in parse_w_walk(fields[6])
                    .ok_or_else(|| invalid_gfa(line_no, "W line has a malformed walk"))?
                {
                    let Some(&node) = id_to_idx.get(id) else {
                        return Err(invalid_gfa(line_no, "W line references unknown segment"));
                    };
                    steps.push(Step { node, rev });
                }
                if !steps.is_empty() {
                    paths.push(Path {
                        name: fields[3].to_string(),
                        steps,
                    });
                }
            }
            "H" | "C" | "J" => {}
            _ => {}
        }
    }

    Ok(Graph { segments, paths })
}

fn invalid_gfa(line_no: usize, msg: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("invalid GFA at line {}: {}", line_no + 1, msg),
    )
}

fn parse_p_step(step: &str) -> Option<(&str, bool)> {
    if let Some(id) = step.strip_suffix('+') {
        Some((id, false))
    } else {
        step.strip_suffix('-').map(|id| (id, true))
    }
}

fn parse_w_walk(walk: &str) -> Option<Vec<(&str, bool)>> {
    let mut steps = Vec::new();
    let mut start = None;
    let mut rev = false;
    for (idx, ch) in walk.char_indices() {
        if ch != '>' && ch != '<' {
            continue;
        }
        if let Some(prev_start) = start {
            if prev_start == idx {
                return None;
            }
            steps.push((&walk[prev_start..idx], rev));
        }
        start = Some(idx + ch.len_utf8());
        rev = ch == '<';
    }
    let prev_start = start?;
    if prev_start == walk.len() {
        return None;
    }
    steps.push((&walk[prev_start..], rev));
    Some(steps)
}

fn path_sequence(graph: &Graph, path: &Path) -> io::Result<Vec<u8>> {
    let mut seq = Vec::new();
    for step in &path.steps {
        let segment = &graph.segments[step.node].seq;
        if step.rev {
            seq.extend(reverse_complement(segment));
        } else {
            seq.extend_from_slice(segment);
        }
    }
    Ok(seq)
}

fn range_sequence(graph: &Graph, range: &PathRange) -> Vec<u8> {
    let mut seq = Vec::with_capacity(range.sequence.len());
    let path = &graph.paths[range.path_idx];
    for step_idx in range.begin_step..range.end_step {
        let step = path.steps[step_idx];
        let segment = &graph.segments[step.node].seq;
        if step.rev {
            seq.extend(reverse_complement(segment));
        } else {
            seq.extend_from_slice(segment);
        }
    }
    seq
}

fn traversal_stats_from_lengths(mut lengths: Vec<usize>) -> TraversalStats {
    if lengths.is_empty() {
        return TraversalStats::default();
    }
    lengths.sort_unstable();

    TraversalStats {
        count: lengths.len(),
        min_len: lengths[0],
        median_len: lengths[lengths.len() / 2],
        p90_len: lengths[percentile_index(lengths.len(), 90, 100)],
        max_len: *lengths.last().unwrap_or(&0),
        total_len: lengths.iter().sum(),
    }
}

fn candidate_selection_method(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> ResolutionMethod {
    match config.method {
        ResolutionMethod::Auto => auto_method_by_median(candidate.traversal_stats, config),
        ResolutionMethod::Hierarchical => hierarchical_method_by_level(candidate.level),
        ResolutionMethod::IterativeMultiLevel | ResolutionMethod::CoverageMultiBubble => {
            auto_method_by_median(candidate.traversal_stats, config)
        }
        ResolutionMethod::ChainPovu => ResolutionMethod::Poasta,
        ResolutionMethod::TopFlubbleSweepga => ResolutionMethod::Sweepga,
        method => method,
    }
}

/// Size-stratified aligner routing by **median** traversal length, per
/// docs/crush-architecture-spec.md §Phase-2 and the fix proposal in
/// docs/crush-aligner-failure-trace.md §"Fix sketch":
///
/// - `median < auto_spoa_max_traversal_len`  → sPOA       (default <1 kb)
/// - `median < auto_poasta_max_traversal_len` → POASTA    (default 1 kb..10 kb)
/// - otherwise                                → sweepga   (default ≥10 kb)
///
/// The decision variable is **median**, not max: a small-median bubble with a
/// single long outlier should still go to sPOA, because the outlier becomes a
/// one-off insertion arc that sPOA represents cleanly. A `*_max_traversal_len`
/// value of 0 disables that tier (the next tier takes over).
///
/// **2-tier variant** (engine flag `auto-2tier=true`, or
/// `auto-spoa-max-traversal-len=0`): skip sPOA entirely so POASTA handles all
/// `median < auto_poasta_max_traversal_len` bubbles and sweepga handles the
/// rest. Rationale (docs/crush-aligner-speed-study.md §recommendation,
/// docs/crush-exp-hybrid-sweepga-poasta.md): POASTA is 84× faster than sPOA on
/// the canonical bimodal C4 plan and produces a cleaner result subgraph.
/// sPOA's role was linearizing small recurrent motifs; POASTA can do that too,
/// trading slightly more nodes for cleaner alignment and a faster wall.
fn auto_method_by_median(
    traversal_stats: TraversalStats,
    config: &ResolutionConfig,
) -> ResolutionMethod {
    let median = traversal_stats.median_len;
    if config.auto_spoa_max_traversal_len > 0 && median < config.auto_spoa_max_traversal_len {
        ResolutionMethod::Poa
    } else if config.auto_poasta_max_traversal_len > 0
        && median < config.auto_poasta_max_traversal_len
    {
        ResolutionMethod::Poasta
    } else {
        ResolutionMethod::Sweepga
    }
}

/// Depth-based aligner dispatch for `method=hierarchical`:
///
/// - `level == 0` (top-level / root bubble) → sweepga+seqwish (all-vs-all
///   alignment + induction). Big structure gets compacted by sweepga.
/// - `level >= 1` (every sub-bubble inside a resolved bubble)      → POASTA.
///   POASTA cleans up local structure inside what sweepga produced.
///
/// This decouples aligner choice from bubble size and ties it to the
/// provenance-tree depth. There is intentionally no size threshold: sweepga
/// only runs at the root and POASTA runs at every interior node, regardless
/// of bp scale. See `docs/crush-hierarchical.md`.
fn hierarchical_method_by_level(level: usize) -> ResolutionMethod {
    if level == 0 {
        ResolutionMethod::Sweepga
    } else {
        ResolutionMethod::Poasta
    }
}

fn candidate_selection_priority(candidate: &BubbleCandidate, config: &ResolutionConfig) -> u8 {
    if config.method != ResolutionMethod::Auto
        && config.method != ResolutionMethod::Hierarchical
        && config.method != ResolutionMethod::ChainPovu
    {
        return 0;
    }
    match candidate_selection_method(candidate, config) {
        ResolutionMethod::Allwave | ResolutionMethod::Sweepga | ResolutionMethod::Wfmash => 0,
        ResolutionMethod::Poa
        | ResolutionMethod::Poasta
        | ResolutionMethod::StarBiwfa
        | ResolutionMethod::ChainGreedy => 1,
        ResolutionMethod::Auto => unreachable!("auto candidate method should be resolved"),
        ResolutionMethod::Hierarchical => {
            unreachable!("hierarchical candidate method should be resolved")
        }
        ResolutionMethod::ChainPovu => {
            unreachable!("chain-povu candidate method should be resolved")
        }
        ResolutionMethod::TopFlubbleSweepga => {
            unreachable!("top-flubble-sweepga candidate method should be resolved")
        }
        ResolutionMethod::IterativeMultiLevel | ResolutionMethod::CoverageMultiBubble => {
            unreachable!("iterative/coverage multi-level candidate method should be resolved")
        }
    }
}

fn percentile_index(len: usize, numerator: usize, denominator: usize) -> usize {
    if len == 0 || denominator == 0 {
        return 0;
    }
    let rank = len.saturating_mul(numerator).div_ceil(denominator);
    rank.saturating_sub(1).min(len - 1)
}

fn format_candidate_length_summary(label: &str, candidates: &[BubbleCandidate]) -> String {
    if candidates.is_empty() {
        return format!("{label} none");
    }

    let mut max_lengths = candidates
        .iter()
        .map(|candidate| candidate.traversal_stats.max_len)
        .collect::<Vec<_>>();
    let mut median_lengths = candidates
        .iter()
        .map(|candidate| candidate.traversal_stats.median_len)
        .collect::<Vec<_>>();
    let mut p90_lengths = candidates
        .iter()
        .map(|candidate| candidate.traversal_stats.p90_len)
        .collect::<Vec<_>>();
    let traversal_counts = candidates
        .iter()
        .map(|candidate| candidate.traversal_stats.count)
        .collect::<Vec<_>>();
    let max_total = candidates
        .iter()
        .map(|candidate| candidate.traversal_stats.total_len)
        .max()
        .unwrap_or(0);
    let max_step_savings = candidates
        .iter()
        .map(BubbleCandidate::estimated_step_savings)
        .max()
        .unwrap_or(0);
    let max_root_span = candidates
        .iter()
        .map(|candidate| candidate.root_span)
        .max()
        .unwrap_or(0);

    let max_len = *max_lengths.iter().max().unwrap_or(&0);
    let max_median = *median_lengths.iter().max().unwrap_or(&0);
    let max_p90 = *p90_lengths.iter().max().unwrap_or(&0);
    let max_traversals = *traversal_counts.iter().max().unwrap_or(&0);

    format!(
        "{label} n={}, max-len median/max={}/{}, median-len median/max={}/{}, p90-len median/max={}/{}, traversals max={}, total max={}, step-savings max={}, root-span max={}",
        candidates.len(),
        median_value(&mut max_lengths),
        max_len,
        median_value(&mut median_lengths),
        max_median,
        median_value(&mut p90_lengths),
        max_p90,
        max_traversals,
        max_total,
        max_step_savings,
        max_root_span
    )
}

fn format_chain_candidate_summary(label: &str, chains: &[ChainCandidate]) -> String {
    if chains.is_empty() {
        return format!("{label} none");
    }

    let source_bubbles = chains
        .iter()
        .map(|chain| chain.source_bubbles)
        .collect::<Vec<_>>();
    let root_spans = chains
        .iter()
        .map(|chain| chain.candidate.root_span)
        .collect::<Vec<_>>();
    let traversal_totals = chains
        .iter()
        .map(|chain| chain.candidate.traversal_stats.total_len)
        .collect::<Vec<_>>();
    let max_lengths = chains
        .iter()
        .map(|chain| chain.candidate.traversal_stats.max_len)
        .collect::<Vec<_>>();
    let median_lengths = chains
        .iter()
        .map(|chain| chain.candidate.traversal_stats.median_len)
        .collect::<Vec<_>>();
    let traversal_counts = chains
        .iter()
        .map(|chain| chain.candidate.traversal_stats.count)
        .collect::<Vec<_>>();

    let (bubble_median, bubble_p90, bubble_max) = distribution_triplet(&source_bubbles);
    let (span_median, span_p90, span_max) = distribution_triplet(&root_spans);
    let (total_median, total_p90, total_max) = distribution_triplet(&traversal_totals);
    let (max_len_median, max_len_p90, max_len_max) = distribution_triplet(&max_lengths);
    let (median_len_median, median_len_p90, median_len_max) = distribution_triplet(&median_lengths);
    let (_, _, traversal_count_max) = distribution_triplet(&traversal_counts);

    format!(
        "{label} n={}, source-bubbles median/p90/max={}/{}/{}, root-span median/p90/max={}/{}/{}, total-len median/p90/max={}/{}/{}, max-len median/p90/max={}/{}/{}, median-len median/p90/max={}/{}/{}, traversals max={}",
        chains.len(),
        bubble_median,
        bubble_p90,
        bubble_max,
        span_median,
        span_p90,
        span_max,
        total_median,
        total_p90,
        total_max,
        max_len_median,
        max_len_p90,
        max_len_max,
        median_len_median,
        median_len_p90,
        median_len_max,
        traversal_count_max
    )
}

fn distribution_triplet(values: &[usize]) -> (usize, usize, usize) {
    if values.is_empty() {
        return (0, 0, 0);
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    (
        sorted[sorted.len() / 2],
        sorted[percentile_index(sorted.len(), 90, 100)],
        *sorted.last().unwrap_or(&0),
    )
}

fn median_value(values: &mut [usize]) -> usize {
    if values.is_empty() {
        return 0;
    }
    values.sort_unstable();
    values[values.len() / 2]
}

impl GraphQuality {
    fn summary(&self) -> String {
        format!(
            "score={}, segments={}, segment-bp={}, links={}, path-steps={}, mean-node-coverage={:.2}, bp-weighted-node-coverage={:.2}, node-coverage-p10/median/p90={}/{}/{}, singleton-nodes={} ({}bp), high-coverage-nodes>={}={} ({}bp), trivial_stringy={}, ws-total={}, ws-p99={}, ws-max={}, ws-long>={}bp={}",
            self.score,
            self.segments,
            self.segment_bp,
            self.links,
            self.path_steps,
            self.node_coverage_mean,
            self.node_coverage_bp_weighted_mean,
            self.node_coverage_p10,
            self.node_coverage_median,
            self.node_coverage_p90,
            self.singleton_nodes,
            self.singleton_bp,
            self.high_coverage_threshold,
            self.high_coverage_nodes,
            self.high_coverage_bp,
            self.trivial_stringy,
            self.path_white_space_bp_total,
            self.path_white_space_bp_p99,
            self.path_white_space_bp_max,
            GRAPH_QUALITY_LONG_BRIDGE_BP,
            self.path_white_space_long_bridges
        )
    }
}

fn graph_quality(graph: &Graph) -> GraphQuality {
    let segments = graph.segments.len();
    let segment_bp = graph
        .segments
        .iter()
        .map(|segment| segment.seq.len())
        .sum::<usize>();
    let path_steps = graph
        .paths
        .iter()
        .map(|path| path.steps.len())
        .sum::<usize>();
    let mut node_visits = vec![0usize; segments];
    for path in &graph.paths {
        for step in &path.steps {
            if step.node < node_visits.len() {
                node_visits[step.node] += 1;
            }
        }
    }
    let total_node_visits = node_visits.iter().copied().sum::<usize>();
    let node_coverage_mean = if segments == 0 {
        0.0
    } else {
        total_node_visits as f64 / segments as f64
    };
    let mut sorted_node_visits = node_visits.clone();
    sorted_node_visits.sort_unstable();
    let node_coverage_p10 = quantile_usize(&sorted_node_visits, 10, 100);
    let node_coverage_median = quantile_usize(&sorted_node_visits, 50, 100);
    let node_coverage_p90 = quantile_usize(&sorted_node_visits, 90, 100);
    let high_coverage_threshold = graph.paths.len().div_ceil(2).max(2);
    let mut weighted_visit_bp = 0u128;
    let mut singleton_nodes = 0usize;
    let mut singleton_bp = 0usize;
    let mut high_coverage_nodes = 0usize;
    let mut high_coverage_bp = 0usize;
    for (visits, segment) in node_visits.iter().copied().zip(&graph.segments) {
        let len = segment.seq.len();
        weighted_visit_bp =
            weighted_visit_bp.saturating_add((visits as u128).saturating_mul(len as u128));
        if visits == 1 {
            singleton_nodes += 1;
            singleton_bp = singleton_bp.saturating_add(len);
        }
        if visits >= high_coverage_threshold {
            high_coverage_nodes += 1;
            high_coverage_bp = high_coverage_bp.saturating_add(len);
        }
    }
    let node_coverage_bp_weighted_mean = if segment_bp == 0 {
        0.0
    } else {
        weighted_visit_bp as f64 / segment_bp as f64
    };
    let mut links = BTreeSet::new();
    let mut white_space_gaps = Vec::new();
    let segment_starts = segment_starts(graph);
    for path in &graph.paths {
        for pair in path.steps.windows(2) {
            let from = pair[0];
            let to = pair[1];
            links.insert((from.node, from.rev, to.node, to.rev));
            white_space_gaps.push(ordered_gap_bp(graph, &segment_starts, from.node, to.node));
        }
    }
    white_space_gaps.sort_unstable();
    let path_white_space_bp_total_u128 = white_space_gaps
        .iter()
        .map(|&gap| gap as u128)
        .sum::<u128>();
    let path_white_space_bp_total = path_white_space_bp_total_u128.min(u64::MAX as u128) as u64;
    let path_white_space_bp_p99 = quantile_usize(&white_space_gaps, 99, 100);
    let path_white_space_bp_max = white_space_gaps.last().copied().unwrap_or(0);
    let path_white_space_long_bridges = white_space_gaps
        .iter()
        .filter(|&&gap| gap >= GRAPH_QUALITY_LONG_BRIDGE_BP)
        .count();

    // This is a visual-tail score, not a raw complexity count. Crushing can
    // legitimately increase path steps while improving the graph the user sees:
    // fewer long empty jumps, tighter high-depth node reuse, and shorter
    // segment sequence. This is logged for diagnostics only; acceptance is
    // based on exact replacement path validation.
    let score = (segment_bp as u128)
        .saturating_add(path_steps as u128 / 64)
        .saturating_add((links.len() as u128).saturating_mul(2))
        .saturating_add((quality_tail_bp(path_white_space_bp_p99) as u128).saturating_mul(1024))
        .saturating_add((quality_tail_bp(path_white_space_bp_max) as u128).saturating_mul(16))
        .saturating_add((path_white_space_long_bridges as u128).saturating_mul(128))
        .saturating_add(path_white_space_bp_total_u128 / 8192)
        .min(u64::MAX as u128) as u64;

    GraphQuality {
        segments,
        segment_bp,
        links: links.len(),
        path_steps,
        node_coverage_mean,
        node_coverage_bp_weighted_mean,
        node_coverage_p10,
        node_coverage_median,
        node_coverage_p90,
        singleton_nodes,
        singleton_bp,
        high_coverage_threshold,
        high_coverage_nodes,
        high_coverage_bp,
        trivial_stringy: trivial_stringy_bubble_count(graph),
        path_white_space_bp_total,
        path_white_space_bp_p99,
        path_white_space_bp_max,
        path_white_space_long_bridges,
        score,
    }
}

fn trivial_stringy_bubble_count(graph: &Graph) -> usize {
    let n_paths = graph.paths.len();
    if n_paths < 2 {
        return 0;
    }
    let mut visits = vec![0usize; graph.segments.len()];
    for path in &graph.paths {
        for step in &path.steps {
            if step.node < visits.len() {
                visits[step.node] += 1;
            }
        }
    }
    let anchor_threshold = n_paths.saturating_mul(8) / 10;
    if anchor_threshold == 0 {
        return 0;
    }
    let anchors = visits
        .iter()
        .enumerate()
        .filter_map(|(idx, &count)| (count >= anchor_threshold).then_some(idx))
        .collect::<FxHashSet<_>>();
    if anchors.len() < 2 {
        return 0;
    }

    let mut bubbles: FxHashMap<(Step, Step), Vec<Vec<Step>>> = FxHashMap::default();
    for path in &graph.paths {
        let anchor_positions = path
            .steps
            .iter()
            .enumerate()
            .filter_map(|(idx, step)| anchors.contains(&step.node).then_some((idx, *step)))
            .collect::<Vec<_>>();
        for pair in anchor_positions.windows(2) {
            let (left_idx, left_step) = pair[0];
            let (right_idx, right_step) = pair[1];
            let mut interior = path.steps[left_idx + 1..right_idx].to_vec();
            let left_key = (left_step.node, left_step.rev);
            let right_key = (right_step.node, right_step.rev);
            let key = if left_key <= right_key {
                (left_step, right_step)
            } else {
                interior.reverse();
                (right_step, left_step)
            };
            bubbles.entry(key).or_default().push(interior);
        }
    }

    bubbles
        .values()
        .filter(|traversals| {
            let mut distinct = FxHashSet::<Vec<Step>>::default();
            let mut lengths = Vec::with_capacity(traversals.len());
            let mut interior_nodes = FxHashSet::<usize>::default();
            for traversal in *traversals {
                distinct.insert(traversal.clone());
                let mut len = 0usize;
                for step in traversal {
                    if let Some(segment) = graph.segments.get(step.node) {
                        len += segment.seq.len();
                    }
                    interior_nodes.insert(step.node);
                }
                lengths.push(len);
            }
            if lengths.is_empty() {
                return false;
            }
            lengths.sort_unstable();
            let median_len = lengths[lengths.len() / 2];
            distinct.len() <= 10 && median_len <= 500 && interior_nodes.len() >= 5
        })
        .count()
}

fn quality_tail_bp(gap_bp: usize) -> usize {
    gap_bp.saturating_sub(GRAPH_QUALITY_WHITE_SPACE_FLOOR_BP)
}

fn segment_starts(graph: &Graph) -> Vec<usize> {
    let mut starts = Vec::with_capacity(graph.segments.len());
    let mut offset = 0usize;
    for segment in &graph.segments {
        starts.push(offset);
        offset = offset.saturating_add(segment.seq.len());
    }
    starts
}

fn ordered_gap_bp(graph: &Graph, starts: &[usize], from: usize, to: usize) -> usize {
    let Some(&from_start) = starts.get(from) else {
        return 0;
    };
    let Some(&to_start) = starts.get(to) else {
        return 0;
    };
    let from_end = from_start.saturating_add(graph.segments[from].seq.len());
    let to_end = to_start.saturating_add(graph.segments[to].seq.len());
    if from_end <= to_start {
        to_start.saturating_sub(from_end)
    } else if to_end <= from_start {
        from_start.saturating_sub(to_end)
    } else {
        0
    }
}

fn quantile_usize(values: &[usize], numerator: usize, denominator: usize) -> usize {
    if values.is_empty() || denominator == 0 {
        return 0;
    }
    values[percentile_index(values.len(), numerator, denominator)]
}

fn path_positions(graph: &Graph, path_idx: usize) -> Vec<usize> {
    let path = &graph.paths[path_idx];
    let mut positions = Vec::with_capacity(path.steps.len() + 1);
    positions.push(0);
    for step in &path.steps {
        let next = positions.last().copied().unwrap_or(0) + graph.segments[step.node].seq.len();
        positions.push(next);
    }
    positions
}

/// One candidate discovered from a POVU pass, before non-overlap selection.
/// Carries enough POVU metadata for the tree-driven level-descent loop to
/// decide whether the site should be admitted into the current round.
#[derive(Clone, Debug)]
struct DiscoveredCandidate {
    candidate: BubbleCandidate,
    povu_site_id: String,
    povu_parent_id: Option<String>,
    povu_level: usize,
    is_leaf: bool,
}

/// Run POVU on `graph` and return all polymorphic candidates with their
/// POVU site metadata, before any tree-driven filtering or non-overlap
/// selection. Pulled out of `find_candidate_frontier` so the same discovery
/// can be reused by `build_initial_bubble_tree`.
fn discover_all_candidates(
    graph: &Graph,
    config: &ResolutionConfig,
    emit_logs: bool,
) -> io::Result<(Vec<DiscoveredCandidate>, usize, RoundDiscoveryTimings)> {
    if graph.paths.is_empty() || graph.paths[0].steps.len() < 2 {
        return Ok((Vec::new(), 0, RoundDiscoveryTimings::default()));
    }
    let root_path_idx = 0;
    let root_path = &graph.paths[root_path_idx];
    let path_positions_by_path: Vec<Vec<usize>> = graph
        .paths
        .iter()
        .enumerate()
        .map(|(idx, _)| path_positions(graph, idx))
        .collect();
    let root_positions = &path_positions_by_path[root_path_idx];
    let render_start = Instant::now();
    if emit_logs {
        log::info!(
            "crush discovery: rendering working graph with {} segment(s), {} path(s)",
            graph.segments.len(),
            graph.paths.len()
        );
    }
    let rendered = render_graph(graph);
    let render_elapsed = render_start.elapsed();
    let parse_start = Instant::now();
    if emit_logs {
        log::info!(
            "crush discovery: parsing rendered graph with POVU ({} bytes)",
            rendered.len()
        );
    }
    let native_graph = povu::NativeGfa::parse(&rendered).map_err(povu_to_io_error)?;
    let parse_elapsed = parse_start.elapsed();
    let reference_names = vec![root_path.name.clone()];
    let decompose_start = Instant::now();
    if emit_logs {
        log::info!(
            "crush discovery: decomposing rendered graph with POVU using root '{}'",
            root_path.name
        );
    }
    let decomposition = native_graph
        .decompose_flubbles(&reference_names)
        .map_err(povu_to_io_error)?;
    let decompose_elapsed = decompose_start.elapsed();
    let id_map_start = Instant::now();
    let id_to_idx = graph
        .segments
        .iter()
        .enumerate()
        .map(|(idx, segment)| (segment.id.as_str(), idx))
        .collect::<FxHashMap<_, _>>();
    let id_map_elapsed = id_map_start.elapsed();
    let path_index_start = Instant::now();
    let path_step_indexes: Vec<FxHashMap<Step, Vec<usize>>> =
        graph.paths.par_iter().map(path_step_index).collect();
    let path_index_elapsed = path_index_start.elapsed();

    let sites_seen = decomposition.sites.len();
    let candidate_start = Instant::now();
    let candidates: Vec<DiscoveredCandidate> = decomposition
        .sites
        .par_iter()
        .filter_map(|site| {
            let begin = site.reference_start_step;
            let exit_step = site.reference_end_step;
            if exit_step + 1 >= root_positions.len() {
                return None;
            }
            let entry = step_from_povu(&id_to_idx, &site.start)?;
            let exit = step_from_povu(&id_to_idx, &site.end)?;
            let root_span = root_positions[exit_step + 1] - root_positions[begin];

            let mut ranges = Vec::new();
            for (path_idx, path_index) in path_step_indexes.iter().enumerate() {
                if let Some((range_begin, range_end)) = unique_anchor_range(path_index, entry, exit)
                {
                    ranges.push(PathRange {
                        path_idx,
                        begin_step: range_begin,
                        end_step: range_end,
                        sequence: Vec::new(),
                        extended_sequence: Vec::new(),
                        left_flank_bp: 0,
                        right_flank_bp: 0,
                    });
                }
            }

            if ranges.len() < 2 {
                return None;
            }

            let lengths = ranges
                .iter()
                .map(|range| {
                    let positions = &path_positions_by_path[range.path_idx];
                    positions[range.end_step] - positions[range.begin_step]
                })
                .collect::<Vec<_>>();
            let traversal_stats = traversal_stats_from_lengths(lengths);
            let total_steps = ranges
                .iter()
                .map(|range| range.end_step.saturating_sub(range.begin_step))
                .sum::<usize>();
            let mut unique_steps = FxHashSet::default();
            for range in &ranges {
                let path = &graph.paths[range.path_idx];
                for step_idx in range.begin_step..range.end_step {
                    unique_steps.insert(path.steps[step_idx]);
                }
            }
            // Tree-driven and experimental multi-bubble modes need complete
            // POVU topology. `min-traversal-len` is a coarse search knob for
            // legacy one-bubble modes only; it must not veto candidates in
            // modes where the purpose is to learn which graph rewrites are
            // valid after exact path-preservation checks.
            let apply_min_len_filter = match config.method {
                ResolutionMethod::Hierarchical => site.level == 0,
                ResolutionMethod::ChainPovu
                | ResolutionMethod::TopFlubbleSweepga
                | ResolutionMethod::IterativeMultiLevel
                | ResolutionMethod::CoverageMultiBubble => false,
                _ => true,
            };
            if apply_min_len_filter && traversal_stats.max_len < config.min_traversal_len {
                return None;
            }
            let signature = candidate_signature(
                graph,
                root_path_idx,
                root_positions,
                begin,
                exit_step,
                &ranges,
                &path_positions_by_path,
            );
            Some(DiscoveredCandidate {
                candidate: BubbleCandidate {
                    ranges,
                    signature,
                    root_start_step: begin,
                    root_end_step: exit_step,
                    root_span,
                    total_steps,
                    unique_steps: unique_steps.len(),
                    traversal_stats,
                    level: site.level,
                },
                povu_site_id: site.id.clone(),
                povu_parent_id: site.parent_id.clone(),
                povu_level: site.level,
                is_leaf: site.is_leaf,
            })
        })
        .collect();
    let candidate_elapsed = candidate_start.elapsed();
    Ok((
        candidates,
        sites_seen,
        RoundDiscoveryTimings {
            render: render_elapsed,
            povu_parse: parse_elapsed,
            povu_decompose: decompose_elapsed,
            id_map: id_map_elapsed,
            path_index: path_index_elapsed,
            candidate_build: candidate_elapsed,
        },
    ))
}

#[derive(Clone, Copy, Debug, Default)]
struct RoundDiscoveryTimings {
    render: std::time::Duration,
    povu_parse: std::time::Duration,
    povu_decompose: std::time::Duration,
    id_map: std::time::Duration,
    path_index: std::time::Duration,
    candidate_build: std::time::Duration,
}

fn find_candidate_frontier(
    graph: &Graph,
    config: &ResolutionConfig,
    seen: &FxHashSet<String>,
    emit_logs: bool,
) -> io::Result<CandidateFrontier> {
    // Bottom-up leaf-driven frontier — kept for legacy callers (the two
    // `phase6_*` unit tests that exercise this path directly). The tree-driven
    // level-descent loop in `resolve_graph_bubbles` does NOT call this; it
    // calls `discover_all_candidates` directly and applies its own
    // tree-driven filter.
    let (all_candidates, sites_seen, timings) = discover_all_candidates(graph, config, emit_logs)?;
    let leaves_seen = all_candidates.iter().filter(|c| c.is_leaf).count();
    let candidates: Vec<BubbleCandidate> = all_candidates
        .into_iter()
        .filter(|d| d.is_leaf)
        .filter(|d| !seen.contains(&d.candidate.signature))
        .map(|d| d.candidate)
        .collect();
    let candidates_seen = candidates.len();
    finalize_frontier(
        graph,
        config,
        candidates,
        sites_seen,
        candidates_seen,
        emit_logs,
        timings,
        leaves_seen,
    )
}

fn find_path_walk_chain_frontier(
    graph: &Graph,
    config: &ResolutionConfig,
    seen: &FxHashSet<String>,
    emit_logs: bool,
) -> io::Result<ChainFrontier> {
    let (all_discovered, sites_seen, timings) = discover_all_candidates(graph, config, emit_logs)?;
    let discovered_candidates = all_discovered.len();
    let candidates = all_discovered
        .into_iter()
        .filter(|d| !seen.contains(&d.candidate.signature))
        .map(|d| d.candidate)
        .collect::<Vec<_>>();
    let minimal = reference_minimal_path_candidates(candidates);
    let path_walk_candidates = reference_path_walk_candidates(minimal);
    let path_walk_candidate_count = path_walk_candidates.len();
    let chains = greedy_path_walk_chain_candidates(graph, config, path_walk_candidates)?
        .into_iter()
        .filter(|chain| !seen.contains(&chain.candidate.signature))
        .collect::<Vec<_>>();
    materialize_chain_frontier(
        graph,
        config,
        chains,
        sites_seen,
        discovered_candidates,
        path_walk_candidate_count,
        emit_logs,
        timings,
    )
}

/// Keep only reference-minimal candidates using interval containment on the
/// root path. This is deliberately not POVU-tree logic: a site is considered
/// a container only if another discovered candidate has a strictly contained
/// root-path interval.
fn reference_minimal_path_candidates(candidates: Vec<BubbleCandidate>) -> Vec<BubbleCandidate> {
    let mut contained = vec![false; candidates.len()];
    for (idx, candidate) in candidates.iter().enumerate() {
        contained[idx] = candidates.iter().enumerate().any(|(other_idx, other)| {
            idx != other_idx
                && candidate.root_start_step <= other.root_start_step
                && candidate.root_end_step >= other.root_end_step
                && (candidate.root_start_step < other.root_start_step
                    || candidate.root_end_step > other.root_end_step)
        });
    }
    candidates
        .into_iter()
        .enumerate()
        .filter_map(|(idx, candidate)| (!contained[idx]).then_some(candidate))
        .collect()
}

/// Walk the root path in coordinate order and choose a non-nested stream of
/// bubble candidates. Adjacent candidates are allowed to share a boundary step
/// (`start == previous_end`), because the chain block will replace the outer
/// range once.
fn reference_path_walk_candidates(mut candidates: Vec<BubbleCandidate>) -> Vec<BubbleCandidate> {
    candidates.sort_by(|a, b| {
        a.root_start_step
            .cmp(&b.root_start_step)
            .then_with(|| a.root_end_step.cmp(&b.root_end_step))
            .then_with(|| {
                b.traversal_stats
                    .total_len
                    .cmp(&a.traversal_stats.total_len)
            })
            .then_with(|| a.signature.cmp(&b.signature))
    });

    let mut walked = Vec::new();
    let mut last_end: Option<usize> = None;
    for candidate in candidates {
        if candidate.root_end_step < candidate.root_start_step {
            continue;
        }
        if let Some(end) = last_end {
            if candidate.root_start_step < end {
                continue;
            }
        }
        last_end = Some(candidate.root_end_step);
        walked.push(candidate);
    }
    walked
}

fn greedy_path_walk_chain_candidates(
    graph: &Graph,
    config: &ResolutionConfig,
    walk_candidates: Vec<BubbleCandidate>,
) -> io::Result<Vec<ChainCandidate>> {
    if walk_candidates.is_empty() || graph.paths.is_empty() {
        return Ok(Vec::new());
    }
    let path_positions_by_path: Vec<Vec<usize>> = (0..graph.paths.len())
        .map(|idx| path_positions(graph, idx))
        .collect();
    let path_step_indexes: Vec<FxHashMap<Step, Vec<usize>>> =
        graph.paths.iter().map(path_step_index).collect();
    let root_positions = &path_positions_by_path[0];
    let target_bp = config.chain_greedy_target_bp.max(1);

    let mut groups: Vec<Vec<BubbleCandidate>> = Vec::new();
    let mut current: Vec<BubbleCandidate> = Vec::new();
    for candidate in walk_candidates {
        if current.is_empty() {
            current.push(candidate);
            continue;
        }

        let current_start = current[0].root_start_step;
        let projected_span = candidate
            .root_end_step
            .checked_add(1)
            .and_then(|end| root_positions.get(end).copied())
            .zip(root_positions.get(current_start).copied())
            .map(|(end_bp, start_bp)| end_bp.saturating_sub(start_bp))
            .unwrap_or(usize::MAX);
        if projected_span > target_bp {
            groups.push(std::mem::take(&mut current));
        }
        current.push(candidate);
    }
    if !current.is_empty() {
        groups.push(current);
    }

    let mut chains = Vec::new();
    for group in groups {
        if let Some(candidate) =
            chain_candidate_from_group(graph, &group, &path_positions_by_path, &path_step_indexes)
        {
            chains.push(ChainCandidate {
                candidate,
                source_bubbles: group.len(),
            });
        }
    }
    Ok(chains)
}

fn chain_candidate_from_group(
    graph: &Graph,
    group: &[BubbleCandidate],
    path_positions_by_path: &[Vec<usize>],
    path_step_indexes: &[FxHashMap<Step, Vec<usize>>],
) -> Option<BubbleCandidate> {
    let first = group.first()?;
    let last = group.last()?;
    let root_path_idx = 0usize;
    let root_path = graph.paths.get(root_path_idx)?;
    let root_start_step = first.root_start_step;
    let root_end_step = last.root_end_step;
    let entry = *root_path.steps.get(root_start_step)?;
    let exit = *root_path.steps.get(root_end_step)?;

    let mut ranges = Vec::new();
    for (path_idx, path_index) in path_step_indexes.iter().enumerate() {
        if let Some((begin_step, end_step)) = unique_anchor_range(path_index, entry, exit) {
            ranges.push(PathRange {
                path_idx,
                begin_step,
                end_step,
                sequence: Vec::new(),
                extended_sequence: Vec::new(),
                left_flank_bp: 0,
                right_flank_bp: 0,
            });
        }
    }
    if ranges.len() < 2 {
        return None;
    }

    let lengths = ranges
        .iter()
        .map(|range| {
            let positions = &path_positions_by_path[range.path_idx];
            positions[range.end_step] - positions[range.begin_step]
        })
        .collect::<Vec<_>>();
    let traversal_stats = traversal_stats_from_lengths(lengths);
    let total_steps = ranges
        .iter()
        .map(|range| range.end_step.saturating_sub(range.begin_step))
        .sum::<usize>();
    let mut unique_steps = FxHashSet::default();
    for range in &ranges {
        let path = &graph.paths[range.path_idx];
        for step_idx in range.begin_step..range.end_step {
            unique_steps.insert(path.steps[step_idx]);
        }
    }

    let root_positions = &path_positions_by_path[root_path_idx];
    let root_span = root_positions[root_end_step + 1] - root_positions[root_start_step];
    let signature = candidate_signature(
        graph,
        root_path_idx,
        root_positions,
        root_start_step,
        root_end_step,
        &ranges,
        path_positions_by_path,
    );

    Some(BubbleCandidate {
        ranges,
        signature,
        root_start_step,
        root_end_step,
        root_span,
        total_steps,
        unique_steps: unique_steps.len(),
        traversal_stats,
        level: 0,
    })
}

fn materialize_chain_frontier(
    graph: &Graph,
    config: &ResolutionConfig,
    chains: Vec<ChainCandidate>,
    sites_seen: usize,
    discovered_candidates: usize,
    path_walk_candidates: usize,
    emit_logs: bool,
    timings: RoundDiscoveryTimings,
) -> io::Result<ChainFrontier> {
    let select_start = Instant::now();
    let chains_formed = chains.len();
    let mut selected = Vec::new();
    let mut occupied_by_path: Vec<Vec<(usize, usize)>> = vec![Vec::new(); graph.paths.len()];
    for chain in chains {
        if candidate_conflicts_with_occupied(&occupied_by_path, &chain.candidate) {
            continue;
        }
        mark_candidate_occupied(&mut occupied_by_path, &chain.candidate);
        selected.push(chain);
    }
    let selected_before_materialize = selected.len();
    let materialize_start = Instant::now();
    selected = selected
        .into_par_iter()
        .filter_map(|mut chain| {
            materialize_candidate_sequences(graph, &mut chain.candidate, config).then_some(chain)
        })
        .collect();
    let materialize_elapsed = materialize_start.elapsed();
    let select_elapsed = select_start.elapsed();

    if emit_logs {
        log::info!(
            "crush chain-greedy discovery detail: render {:.2?}, povu-parse {:.2?}, povu-decompose {:.2?}, id-map {:.2?}, path-index {:.2?}, candidate-build {:.2?}, chain-select+materialize {:.2?} ({} formed -> {} non-overlap -> {} selected, materialize {:.2?})",
            timings.render,
            timings.povu_parse,
            timings.povu_decompose,
            timings.id_map,
            timings.path_index,
            timings.candidate_build,
            select_elapsed,
            chains_formed,
            selected_before_materialize,
            selected.len(),
            materialize_elapsed
        );
    }

    Ok(ChainFrontier {
        selected,
        sites_seen,
        discovered_candidates,
        path_walk_candidates,
        chains_formed,
    })
}

/// Sort + non-overlap-select + materialize sequences. Shared by the legacy
/// leaf-driven path and the tree-driven level-descent loop.
fn finalize_frontier(
    graph: &Graph,
    config: &ResolutionConfig,
    mut candidates: Vec<BubbleCandidate>,
    sites_seen: usize,
    candidates_seen: usize,
    emit_logs: bool,
    timings: RoundDiscoveryTimings,
    leaves_seen: usize,
) -> io::Result<CandidateFrontier> {
    let select_start = Instant::now();
    candidates.sort_by(|a, b| {
        candidate_selection_priority(a, config)
            .cmp(&candidate_selection_priority(b, config))
            .then_with(|| b.estimated_step_savings().cmp(&a.estimated_step_savings()))
            .then_with(|| b.unique_steps.cmp(&a.unique_steps))
            .then_with(|| b.root_step_span().cmp(&a.root_step_span()))
            .then_with(|| b.root_span.cmp(&a.root_span))
            .then_with(|| a.root_start_step.cmp(&b.root_start_step))
            .then_with(|| a.root_end_step.cmp(&b.root_end_step))
    });

    let mut frontier = CandidateFrontier {
        sites_seen,
        candidates_seen,
        ..CandidateFrontier::default()
    };
    let mut occupied_by_path: Vec<Vec<(usize, usize)>> = vec![Vec::new(); graph.paths.len()];
    for candidate in candidates {
        if candidate_conflicts_with_occupied(&occupied_by_path, &candidate) {
            continue;
        }
        mark_candidate_occupied(&mut occupied_by_path, &candidate);
        frontier.selected.push(candidate);
    }
    let selected_before_materialize = frontier.selected.len();
    let materialize_start = Instant::now();
    frontier.selected = frontier
        .selected
        .into_par_iter()
        .filter_map(|mut candidate| {
            materialize_candidate_sequences(graph, &mut candidate, config).then_some(candidate)
        })
        .collect();
    let materialize_elapsed = materialize_start.elapsed();
    let select_elapsed = select_start.elapsed();

    if emit_logs {
        log::info!(
            "crush discovery detail: render {:.2?}, povu-parse {:.2?}, povu-decompose {:.2?}, id-map {:.2?}, path-index {:.2?}, candidate-build {:.2?}, select+materialize {:.2?} ({} -> {} selected, materialize {:.2?}); leaves={}",
            timings.render,
            timings.povu_parse,
            timings.povu_decompose,
            timings.id_map,
            timings.path_index,
            timings.candidate_build,
            select_elapsed,
            selected_before_materialize,
            frontier.selected.len(),
            materialize_elapsed,
            leaves_seen
        );
    }

    Ok(frontier)
}

impl BubbleCandidate {
    fn root_step_span(&self) -> usize {
        self.root_end_step.saturating_sub(self.root_start_step)
    }

    fn estimated_step_savings(&self) -> usize {
        self.total_steps
            .saturating_sub(self.traversal_stats.count)
            .saturating_add(self.unique_steps.saturating_sub(1))
    }
}

#[derive(Clone, Debug)]
struct ChainPovuSelectionNode {
    parent: Option<usize>,
    children: Vec<usize>,
    candidate_idx: usize,
}

#[derive(Clone, Debug, Default)]
struct ChainPovuBlockSelection {
    selected: Vec<BubbleCandidate>,
    tree_nodes: usize,
    roots: usize,
    skipped_seen: usize,
    oversize_internal: usize,
    oversize_leaves: usize,
    merged_descendant_sites: usize,
}

fn select_chain_povu_blocks(
    discovered: &[DiscoveredCandidate],
    max_block_bp: usize,
    resolved_signatures: &FxHashSet<String>,
) -> ChainPovuBlockSelection {
    let mut nodes: Vec<ChainPovuSelectionNode> = Vec::with_capacity(discovered.len());
    let mut id_to_node_idx: FxHashMap<&str, usize> = FxHashMap::default();
    for (candidate_idx, d) in discovered.iter().enumerate() {
        let node_idx = nodes.len();
        id_to_node_idx.insert(d.povu_site_id.as_str(), node_idx);
        nodes.push(ChainPovuSelectionNode {
            parent: None,
            children: Vec::new(),
            candidate_idx,
        });
    }

    for (idx, d) in discovered.iter().enumerate() {
        if let Some(pid) = &d.povu_parent_id {
            if let Some(&parent_idx) = id_to_node_idx.get(pid.as_str()) {
                nodes[idx].parent = Some(parent_idx);
                nodes[parent_idx].children.push(idx);
            }
        }
    }

    let roots = nodes
        .iter()
        .enumerate()
        .filter_map(|(idx, node)| node.parent.is_none().then_some(idx))
        .collect::<Vec<_>>();
    let mut selection = ChainPovuBlockSelection {
        tree_nodes: nodes.len(),
        roots: roots.len(),
        ..ChainPovuBlockSelection::default()
    };
    for root in roots {
        visit_chain_povu_node(
            root,
            &nodes,
            discovered,
            max_block_bp,
            resolved_signatures,
            &mut selection,
        );
    }
    selection
}

fn visit_chain_povu_node(
    node_idx: usize,
    nodes: &[ChainPovuSelectionNode],
    discovered: &[DiscoveredCandidate],
    max_block_bp: usize,
    resolved_signatures: &FxHashSet<String>,
    selection: &mut ChainPovuBlockSelection,
) {
    let node = &nodes[node_idx];
    let candidate = &discovered[node.candidate_idx].candidate;
    if resolved_signatures.contains(&candidate.signature) {
        selection.skipped_seen += 1;
        return;
    }

    let block_bp = chain_povu_block_bp(candidate);
    if block_bp <= max_block_bp {
        let mut selected = candidate.clone();
        selected.level = discovered[node.candidate_idx].povu_level;
        selection.selected.push(selected);
        selection.merged_descendant_sites += chain_povu_descendant_count(node_idx, nodes);
        return;
    }

    if node.children.is_empty() {
        selection.oversize_leaves += 1;
        return;
    }
    selection.oversize_internal += 1;
    for &child in &node.children {
        visit_chain_povu_node(
            child,
            nodes,
            discovered,
            max_block_bp,
            resolved_signatures,
            selection,
        );
    }
}

fn chain_povu_descendant_count(node_idx: usize, nodes: &[ChainPovuSelectionNode]) -> usize {
    nodes[node_idx]
        .children
        .iter()
        .map(|&child| 1 + chain_povu_descendant_count(child, nodes))
        .sum()
}

fn chain_povu_block_bp(candidate: &BubbleCandidate) -> usize {
    // The block is bounded by the selected parent flubble interior. Descendant
    // interiors are nested in that same interval, so summing them would
    // double-count exactly the chains this mode is meant to merge.
    candidate.root_span.max(candidate.traversal_stats.max_len)
}

fn format_bp_distribution(label: &str, values: &[usize]) -> String {
    if values.is_empty() {
        return format!("{label} n=0");
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let sum = sorted.iter().sum::<usize>();
    let p50 = sorted[sorted.len() / 2];
    let p90 = sorted[percentile_index(sorted.len(), 90, 100)];
    let min = sorted[0];
    let max = *sorted.last().unwrap_or(&0);
    format!(
        "{label} n={}, min={}, p50={}, p90={}, max={}, total={}",
        sorted.len(),
        min,
        p50,
        p90,
        max,
        sum
    )
}

/// Build the initial bubble tree by running POVU on the input graph once.
/// Returns the tree (Vec<BubbleNode>) and the index of every root node.
///
/// Parent/child edges follow POVU's `parent_id`. POVU's `level` is preserved
/// on the node. Each node carries a stable bp-keyed signature; the round
/// loop uses that signature both for provenance tracking (no node resolved
/// twice) and for matching post-rewrite global-POVU results back to the
/// active tree frontier when the children must be resolved.
fn build_initial_bubble_tree(
    graph: &Graph,
    config: &ResolutionConfig,
    emit_logs: bool,
) -> io::Result<(Vec<BubbleNode>, Vec<usize>)> {
    let (discovered, sites_seen, _timings) = discover_all_candidates(graph, config, emit_logs)?;
    if emit_logs {
        log::info!(
            "crush tree: initial POVU on input found {} site(s), {} polymorphic candidate(s) of which {} are roots (level==0)",
            sites_seen,
            discovered.len(),
            discovered.iter().filter(|d| d.povu_level == 0).count()
        );
    }
    let root_path_idx = 0;
    let root_positions = path_positions(graph, root_path_idx);

    // Build nodes in POVU-site-id order so children-of can be resolved by id.
    let mut nodes: Vec<BubbleNode> = Vec::with_capacity(discovered.len());
    let mut id_to_node_idx: FxHashMap<String, usize> = FxHashMap::default();
    for d in &discovered {
        let ref_bp_begin = root_positions[d.candidate.root_start_step];
        let ref_bp_end = root_positions[d.candidate.root_end_step + 1];
        let node_idx = nodes.len();
        id_to_node_idx.insert(d.povu_site_id.clone(), node_idx);
        nodes.push(BubbleNode {
            parent: None,
            children: Vec::new(),
            level: d.povu_level,
            state: BubbleState::Unresolved,
            bp_signature: d.candidate.signature.clone(),
            ref_bp_begin,
            ref_bp_end,
            discovery_site_id: d.povu_site_id.clone(),
            local_child_bp_keys: None,
            initial_traversal_max_len: d.candidate.traversal_stats.max_len,
        });
    }
    // Wire parent/children by POVU parent_id.
    for (idx, d) in discovered.iter().enumerate() {
        if let Some(pid) = &d.povu_parent_id {
            if let Some(&parent_idx) = id_to_node_idx.get(pid) {
                nodes[idx].parent = Some(parent_idx);
                nodes[parent_idx].children.push(idx);
            }
        }
    }
    let roots: Vec<usize> = nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.parent.is_none())
        .map(|(i, _)| i)
        .collect();
    if emit_logs {
        log::info!(
            "crush tree: built initial forest with {} node(s), {} root(s)",
            nodes.len(),
            roots.len()
        );
    }
    Ok((nodes, roots))
}

/// Run POVU on the LOCAL replacement graph of a just-resolved bubble and
/// translate each discovered sub-site's coordinates into the GLOBAL working
/// graph's bp-keyed signature, matching the format produced by
/// `candidate_signature`. Returns the list of sub-bubble bp-keys (which
/// become the children of the parent node) plus useful metadata (per-path
/// global bp ranges) so the discovered children can be cross-referenced
/// against the next round's global-POVU pass.
///
/// `parent` carries the BubbleCandidate that was resolved; each of its
/// `ranges[i]` corresponds to local path `i` in `replacement`. We compute
/// the global bp offset for each path as
/// `path_positions_pre_apply[global_path_idx][parent_range.begin_step]`,
/// then add the local bp position of each sub-site inside that local path.
///
/// Returns an empty vec if local POVU finds nothing polymorphic, or if the
/// local graph is too small for POVU to decompose meaningfully.
fn discover_local_subbubble_keys(
    parent: &BubbleCandidate,
    replacement: &Graph,
    parent_pre_apply_path_positions: &[Vec<usize>],
    global_path_names: &[String],
    _config: &ResolutionConfig,
) -> io::Result<Vec<String>> {
    if replacement.paths.is_empty() || replacement.paths[0].steps.len() < 2 {
        return Ok(Vec::new());
    }
    // Find which local path corresponds to the GLOBAL reference (graph.paths[0]).
    // POVU on the local subgraph must be rooted on the same reference path the
    // outer round used, so the resulting sites' reference_start/end coords are
    // directly comparable.
    let global_ref_path_idx = 0usize;
    let local_ref_idx = parent
        .ranges
        .iter()
        .position(|r| r.path_idx == global_ref_path_idx);
    // If the parent bubble doesn't include the global reference path (rare
    // but possible if POVU surfaced a bubble using a different ref), we
    // can't compute global bp keys for sub-bubbles, so skip.
    let Some(local_ref_idx) = local_ref_idx else {
        return Ok(Vec::new());
    };
    let local_ref_path_name = replacement.paths[local_ref_idx].name.clone();

    let rendered = render_graph(replacement);
    let native = match povu::NativeGfa::parse(&rendered) {
        Ok(g) => g,
        Err(_) => return Ok(Vec::new()),
    };
    let decomposition = match native.decompose_flubbles(&[local_ref_path_name]) {
        Ok(d) => d,
        Err(_) => return Ok(Vec::new()),
    };
    if decomposition.sites.is_empty() {
        return Ok(Vec::new());
    }

    // Precompute local path bp positions and step indexes for sub-site lookup.
    let local_path_positions: Vec<Vec<usize>> = (0..replacement.paths.len())
        .map(|i| path_positions(replacement, i))
        .collect();
    let local_path_step_indexes: Vec<FxHashMap<Step, Vec<usize>>> =
        replacement.paths.iter().map(path_step_index).collect();
    let local_id_to_idx: FxHashMap<&str, usize> = replacement
        .segments
        .iter()
        .enumerate()
        .map(|(i, s)| (s.id.as_str(), i))
        .collect();

    log::debug!(
        "crush local-povu: parent sig {:?}, replacement {} seg / {} path, local POVU {} site(s) ({} leaves)",
        parent.signature,
        replacement.segments.len(),
        replacement.paths.len(),
        decomposition.sites.len(),
        decomposition.sites.iter().filter(|s| s.is_leaf).count()
    );
    let mut keys: Vec<String> = Vec::new();
    for site in &decomposition.sites {
        // Add LOCAL LEAVES (is_leaf=true at any depth) as direct children
        // of the just-resolved bubble. This flattens the local hierarchy:
        // a non-leaf level-0 local site is a CONTAINER whose only
        // contribution is its leaf descendants, so we admit those leaves
        // directly rather than descending the container in a later round.
        // For the nested-bubbles fixture this means round 2 sees ALL leaf
        // candidates from the round-1 replacement at once (not just the
        // outer-level partition + a future round for the inner leaves) —
        // which is what gives iterations<=2 even when POA's consensus is
        // multi-level.
        if !site.is_leaf {
            continue;
        }
        // Look up entry/exit Step in the local graph
        let entry = match step_from_povu(&local_id_to_idx, &site.start) {
            Some(e) => e,
            None => continue,
        };
        let exit = match step_from_povu(&local_id_to_idx, &site.end) {
            Some(e) => e,
            None => continue,
        };
        // For each LOCAL path the sub-site traverses, compute global bp range
        let mut per_path_coords: Vec<String> = Vec::new();
        for (local_path_idx, path_idx) in local_path_step_indexes.iter().enumerate() {
            let (lp_begin, lp_end) = match unique_anchor_range(path_idx, entry, exit) {
                Some(r) => r,
                None => continue,
            };
            let local_positions = &local_path_positions[local_path_idx];
            let local_bp_begin = local_positions[lp_begin];
            let local_bp_end = local_positions[lp_end];
            let parent_range = &parent.ranges[local_path_idx];
            let global_path_idx = parent_range.path_idx;
            let parent_path_pos = &parent_pre_apply_path_positions[global_path_idx];
            // The parent's path bp begin in the global graph (= bp position
            // at parent_range.begin_step in the pre-apply working graph).
            let parent_path_bp_begin = parent_path_pos[parent_range.begin_step];
            let global_bp_begin = parent_path_bp_begin + local_bp_begin;
            let global_bp_end = parent_path_bp_begin + local_bp_end;
            per_path_coords.push(format!(
                "{}:{}-{}",
                global_path_names[global_path_idx], global_bp_begin, global_bp_end
            ));
        }
        if per_path_coords.len() < 2 {
            continue;
        }
        per_path_coords.sort();

        // Reference-path global bp range — use the local reference path
        // positions plus the parent's reference bp offset.
        let local_ref_positions = &local_path_positions[local_ref_idx];
        let local_ref_bp_begin = local_ref_positions[site.reference_start_step];
        let local_ref_bp_end = local_ref_positions[site.reference_end_step + 1];
        let parent_ref_range = &parent.ranges[local_ref_idx];
        let parent_ref_path_pos = &parent_pre_apply_path_positions[parent_ref_range.path_idx];
        let parent_ref_bp_begin = parent_ref_path_pos[parent_ref_range.begin_step];
        let global_ref_bp_begin = parent_ref_bp_begin + local_ref_bp_begin;
        let global_ref_bp_end = parent_ref_bp_begin + local_ref_bp_end;

        let key = format!(
            "{}:{}-{}|{}",
            global_path_names[global_ref_path_idx],
            global_ref_bp_begin,
            global_ref_bp_end,
            per_path_coords.join(",")
        );
        keys.push(key);
    }
    Ok(keys)
}

fn candidate_conflicts_with_occupied(
    occupied_by_path: &[Vec<(usize, usize)>],
    candidate: &BubbleCandidate,
) -> bool {
    candidate.ranges.iter().any(|range| {
        range_conflicts_with_occupied(
            &occupied_by_path[range.path_idx],
            range.begin_step,
            range.end_step,
        )
    })
}

fn mark_candidate_occupied(
    occupied_by_path: &mut [Vec<(usize, usize)>],
    candidate: &BubbleCandidate,
) {
    for range in &candidate.ranges {
        insert_occupied_range(
            &mut occupied_by_path[range.path_idx],
            range.begin_step,
            range.end_step,
        );
    }
}

fn range_conflicts_with_occupied(occupied: &[(usize, usize)], begin: usize, end: usize) -> bool {
    let insert_pos = occupied.partition_point(|&(occupied_begin, _)| occupied_begin < begin);
    if insert_pos > 0 {
        let (_, previous_end) = occupied[insert_pos - 1];
        if previous_end > begin {
            return true;
        }
    }
    if let Some(&(next_begin, _)) = occupied.get(insert_pos) {
        if next_begin < end {
            return true;
        }
    }
    false
}

fn insert_occupied_range(occupied: &mut Vec<(usize, usize)>, begin: usize, end: usize) {
    let insert_pos = occupied.partition_point(|&(occupied_begin, _)| occupied_begin < begin);
    occupied.insert(insert_pos, (begin, end));
}

fn step_from_povu(id_to_idx: &FxHashMap<&str, usize>, step: &PovuStep) -> Option<Step> {
    let node = *id_to_idx.get(step.segment.as_str())?;
    Some(Step {
        node,
        rev: matches!(step.strand, PovuStrand::Reverse),
    })
}

fn povu_to_io_error(err: povu::Error) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("POVU flubble decomposition failed: {err}"),
    )
}

fn path_step_index(path: &Path) -> FxHashMap<Step, Vec<usize>> {
    let mut index: FxHashMap<Step, Vec<usize>> = FxHashMap::default();
    for (idx, &step) in path.steps.iter().enumerate() {
        index.entry(step).or_default().push(idx);
    }
    index
}

fn unique_anchor_range(
    path_index: &FxHashMap<Step, Vec<usize>>,
    entry: Step,
    exit: Step,
) -> Option<(usize, usize)> {
    let starts = path_index.get(&entry)?;
    let exits = path_index.get(&exit)?;
    let mut found = None;
    for &i in starts {
        let exit_pos = exits.partition_point(|&j| j <= i);
        if exit_pos >= exits.len() {
            continue;
        }
        let j = exits[exit_pos];
        if found.is_some() {
            return None;
        }
        found = Some((i, j + 1));
    }
    found
}

fn materialize_candidate_sequences(
    graph: &Graph,
    candidate: &mut BubbleCandidate,
    config: &ResolutionConfig,
) -> bool {
    for range in &mut candidate.ranges {
        range.sequence = range_sequence(graph, range);
        range.extended_sequence = Vec::new();
        range.left_flank_bp = 0;
        range.right_flank_bp = 0;
    }
    if config.replacement_flank_bp > 0 {
        materialize_flanked_sequences(graph, candidate, config.replacement_flank_bp);
    }
    // Flanks are aligner context only; skip candidates whose substituted
    // interiors spell the same sequence on every path.
    let Some(first) = candidate.ranges.first() else {
        return false;
    };
    candidate
        .ranges
        .iter()
        .any(|range| range.sequence != first.sequence)
}

/// Compute up to `flank_bp` bp of path-sequence context on each side of each
/// PathRange interior, capped by the available path-prefix / path-suffix bp so
/// the flank never crosses the path edge. The interior step range
/// (`begin_step`..`end_step`) is preserved verbatim — the flank is sequence
/// context fed to the replacement aligner, not new interior the integration
/// would substitute. The aligner-visible sequence becomes
/// `<left_flank><interior><right_flank>` and is stored on
/// `range.extended_sequence`; `left_flank_bp` / `right_flank_bp` track the
/// per-side flank length so the aligner output graph can be clipped back to
/// just the interior before the path-step substitution.
fn materialize_flanked_sequences(graph: &Graph, candidate: &mut BubbleCandidate, flank_bp: usize) {
    for range in &mut candidate.ranges {
        let path = &graph.paths[range.path_idx];
        let left_flank = collect_flank_left(graph, path, range.begin_step, flank_bp);
        let right_flank = collect_flank_right(graph, path, range.end_step, flank_bp);
        if left_flank.is_empty() && right_flank.is_empty() {
            continue;
        }
        let mut extended =
            Vec::with_capacity(left_flank.len() + range.sequence.len() + right_flank.len());
        extended.extend_from_slice(&left_flank);
        extended.extend_from_slice(&range.sequence);
        extended.extend_from_slice(&right_flank);
        range.left_flank_bp = left_flank.len();
        range.right_flank_bp = right_flank.len();
        range.extended_sequence = extended;
    }
}

/// Walk path steps from `begin_step - 1` toward the path start, accumulating
/// step-aligned sequence until at least `flank_bp` bp is collected (the last
/// step included may overshoot — we don't split mid-segment for flanks because
/// the aligner reads bp, not step structure).
fn collect_flank_left(graph: &Graph, path: &Path, begin_step: usize, flank_bp: usize) -> Vec<u8> {
    let mut collected: Vec<u8> = Vec::new();
    let mut step_idx = begin_step;
    while step_idx > 0 && collected.len() < flank_bp {
        step_idx -= 1;
        let step = path.steps[step_idx];
        let segment = &graph.segments[step.node].seq;
        let segment_seq: Vec<u8> = if step.rev {
            reverse_complement(segment)
        } else {
            segment.to_vec()
        };
        let mut prepended = segment_seq;
        prepended.extend_from_slice(&collected);
        collected = prepended;
    }
    if collected.len() > flank_bp {
        let drop = collected.len() - flank_bp;
        collected.drain(..drop);
    }
    collected
}

fn collect_flank_right(graph: &Graph, path: &Path, end_step: usize, flank_bp: usize) -> Vec<u8> {
    let mut collected: Vec<u8> = Vec::new();
    let mut step_idx = end_step;
    while step_idx < path.steps.len() && collected.len() < flank_bp {
        let step = path.steps[step_idx];
        let segment = &graph.segments[step.node].seq;
        if step.rev {
            collected.extend(reverse_complement(segment));
        } else {
            collected.extend_from_slice(segment);
        }
        step_idx += 1;
    }
    if collected.len() > flank_bp {
        collected.truncate(flank_bp);
    }
    collected
}

fn candidate_signature(
    graph: &Graph,
    ref_path_idx: usize,
    ref_positions: &[usize],
    begin: usize,
    exit_step: usize,
    ranges: &[PathRange],
    path_positions_by_path: &[Vec<usize>],
) -> String {
    // Use base-pair positions on every involved path so the signature is stable
    // across rewrites of OTHER bubbles (their interior bytes are preserved
    // byte-for-byte by `path_sequences_equal`, so bp positions outside their
    // span don't shift). Step indices change as segmentation is rewritten, so
    // they cannot be part of the persistent identity of a bubble.
    let mut coords: Vec<String> = ranges
        .iter()
        .map(|r| {
            let positions = &path_positions_by_path[r.path_idx];
            format!(
                "{}:{}-{}",
                graph.paths[r.path_idx].name, positions[r.begin_step], positions[r.end_step]
            )
        })
        .collect();
    coords.sort();
    format!(
        "{}:{}-{}|{}",
        graph.paths[ref_path_idx].name,
        ref_positions[begin],
        ref_positions[exit_step + 1],
        coords.join(",")
    )
}

fn apply_replacement_frontier(
    graph: &Graph,
    plans: &[ReplacementPlan],
    next_id: &mut usize,
) -> io::Result<Graph> {
    let mut replacements_by_path: FxHashMap<usize, Vec<PathReplacement>> = FxHashMap::default();

    for (plan_idx, plan) in plans.iter().enumerate() {
        for (range_idx, range) in plan.candidate.ranges.iter().enumerate() {
            let path =
                plan.replacement.paths.get(range_idx).ok_or_else(|| {
                    io::Error::other("replacement path missing for candidate range")
                })?;
            let steps = path
                .steps
                .iter()
                .map(|step| OutStep {
                    node: OutNode::Replacement(plan_idx, step.node),
                    rev: step.rev,
                })
                .collect::<Vec<_>>();
            replacements_by_path
                .entry(range.path_idx)
                .or_default()
                .push(PathReplacement {
                    begin_step: range.begin_step,
                    end_step: range.end_step,
                    steps,
                });
        }
    }

    for replacements in replacements_by_path.values_mut() {
        replacements.sort_by_key(|replacement| (replacement.begin_step, replacement.end_step));
        for pair in replacements.windows(2) {
            if pair[0].end_step > pair[1].begin_step {
                return Err(io::Error::other(
                    "replacement frontier contains overlapping path ranges",
                ));
            }
        }
    }

    let mut out_paths = Vec::with_capacity(graph.paths.len());
    let mut used_original: FxHashSet<usize> = FxHashSet::default();

    for (path_idx, path) in graph.paths.iter().enumerate() {
        let path_replacements = replacements_by_path
            .get(&path_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let mut replacement_idx = 0usize;
        let mut steps = Vec::new();
        let mut i = 0usize;
        while i < path.steps.len() {
            if let Some(replacement) = path_replacements.get(replacement_idx) {
                if i == replacement.begin_step {
                    for step in &replacement.steps {
                        steps.push(*step);
                    }
                    i = replacement.end_step;
                    replacement_idx += 1;
                    continue;
                }
                if replacement.begin_step < i && i < replacement.end_step {
                    return Err(io::Error::other(
                        "replacement frontier rewrite entered an already replaced range",
                    ));
                }
            }

            let original = path.steps[i];
            used_original.insert(original.node);
            steps.push(OutStep {
                node: OutNode::Original(original.node),
                rev: original.rev,
            });
            i += 1;
        }
        out_paths.push((path.name.clone(), steps));
    }

    let replacement_graphs = plans
        .iter()
        .map(|plan| plan.replacement.clone())
        .collect::<Vec<_>>();
    if log::log_enabled!(log::Level::Debug) {
        let total_replacement_segments: usize =
            replacement_graphs.iter().map(|g| g.segments.len()).sum();
        let total_replacement_bp: usize = replacement_graphs
            .iter()
            .map(|g| g.segments.iter().map(|s| s.seq.len()).sum::<usize>())
            .sum();
        log::debug!(
            "crush apply: {} plan(s); replacement segments total={}, replacement bp total={}; rewriting {} path(s)",
            plans.len(),
            total_replacement_segments,
            total_replacement_bp,
            out_paths.len()
        );
    }
    let rendered = render_rewritten_graph(
        graph,
        &replacement_graphs,
        &used_original,
        &out_paths,
        next_id,
    );
    let next = parse_gfa(&rendered)?;
    if log::log_enabled!(log::Level::Debug) {
        let mut sequence_counts: FxHashMap<&[u8], usize> = FxHashMap::default();
        for segment in &next.segments {
            *sequence_counts.entry(segment.seq.as_slice()).or_insert(0) += 1;
        }
        let duplicate_segments: usize = sequence_counts
            .values()
            .filter(|&&n| n > 1)
            .map(|&n| n - 1)
            .sum();
        log::debug!(
            "crush apply: rendered {} bytes -> parsed graph has {} segment(s), {} path(s), {} unique segment sequence(s), {} duplicate-sequence segment(s)",
            rendered.len(),
            next.segments.len(),
            next.paths.len(),
            sequence_counts.len(),
            duplicate_segments
        );
    }
    if !path_sequences_equal(graph, &next)? {
        return Err(io::Error::other(
            "resolved graph failed exact path-sequence validation",
        ));
    }
    Ok(next)
}

fn build_replacement_with_method(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
    method: ResolutionMethod,
) -> io::Result<Graph> {
    let replacement = build_replacement_with_method_inner(candidate, config, method)?;
    if candidate_has_flank(candidate) {
        let clipped = clip_replacement_to_interior(replacement, candidate)?;
        validate_interior_replacement_paths(&clipped, candidate, method.method_name())?;
        Ok(clipped)
    } else {
        Ok(replacement)
    }
}

#[derive(Clone, Debug)]
struct ReplacementBuildReport {
    replacement: Graph,
    evidence: Option<SweepgaReplacementEvidence>,
}

fn build_replacement_with_method_report(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
    method: ResolutionMethod,
) -> io::Result<ReplacementBuildReport> {
    let (replacement, evidence) = match method {
        ResolutionMethod::Sweepga => {
            let (replacement, evidence) =
                build_sweepga_seqwish_replacement_with_evidence(candidate, config)?;
            (replacement, Some(evidence))
        }
        ResolutionMethod::TopFlubbleSweepga => {
            let (replacement, evidence) =
                build_top_flubble_sweepga_replacement_with_evidence(candidate, config)?;
            (replacement, Some(evidence))
        }
        _ => (
            build_replacement_with_method_inner(candidate, config, method)?,
            None,
        ),
    };
    let replacement = if candidate_has_flank(candidate) {
        let clipped = clip_replacement_to_interior(replacement, candidate)?;
        validate_interior_replacement_paths(&clipped, candidate, method.method_name())?;
        clipped
    } else {
        replacement
    };
    Ok(ReplacementBuildReport {
        replacement,
        evidence,
    })
}

fn build_replacement_with_method_catching_unwind_report(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
    method: ResolutionMethod,
) -> io::Result<ReplacementBuildReport> {
    match panic::catch_unwind(AssertUnwindSafe(|| {
        build_replacement_with_method_report(candidate, config, method)
    })) {
        Ok(result) => result,
        Err(payload) => {
            let panic_message = payload
                .downcast_ref::<&str>()
                .copied()
                .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
                .unwrap_or("unknown panic");
            Err(io::Error::other(format!(
                "{} replacement panicked for traversal_count={}, max_len={}, median_len={}, total_len={}, root_span={}: {}",
                method.method_name(),
                candidate.traversal_stats.count,
                candidate.traversal_stats.max_len,
                candidate.traversal_stats.median_len,
                candidate.traversal_stats.total_len,
                candidate.root_span,
                panic_message
            )))
        }
    }
}

fn validate_interior_replacement_paths(
    replacement: &Graph,
    candidate: &BubbleCandidate,
    method: &str,
) -> io::Result<()> {
    if replacement.paths.len() != candidate.ranges.len() {
        return Err(io::Error::other(format!(
            "{method} clipped replacement emitted {} paths for {} traversals",
            replacement.paths.len(),
            candidate.ranges.len()
        )));
    }
    for (idx, range) in candidate.ranges.iter().enumerate() {
        let observed = path_sequence(replacement, &replacement.paths[idx])?;
        if observed != range.sequence {
            return Err(io::Error::other(format!(
                "{method} clipped replacement path {} has {} bp interior, expected {} bp",
                idx,
                observed.len(),
                range.sequence.len()
            )));
        }
    }
    Ok(())
}

fn build_replacement_with_method_inner(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
    method: ResolutionMethod,
) -> io::Result<Graph> {
    match method {
        ResolutionMethod::Auto => unreachable!("auto is resolved before replacement dispatch"),
        ResolutionMethod::Hierarchical => {
            unreachable!("hierarchical is resolved before replacement dispatch")
        }
        ResolutionMethod::IterativeMultiLevel => {
            unreachable!("iterative multi-level is resolved before replacement dispatch")
        }
        ResolutionMethod::CoverageMultiBubble => {
            unreachable!("coverage multi-bubble is resolved before replacement dispatch")
        }
        ResolutionMethod::ChainPovu => {
            build_chain_povu_smooth_poasta_replacement(candidate, config)
        }
        ResolutionMethod::TopFlubbleSweepga => build_sweepga_seqwish_replacement(candidate, config),
        ResolutionMethod::Poa => build_poa_replacement(candidate, config),
        ResolutionMethod::Poasta => build_poasta_replacement(candidate, config),
        ResolutionMethod::ChainGreedy => build_chain_greedy_replacement(candidate, config),
        ResolutionMethod::StarBiwfa => build_biwfa_inmemory_replacement(candidate, config),
        ResolutionMethod::Allwave => build_allwave_seqwish_replacement(candidate, config),
        ResolutionMethod::Sweepga => build_sweepga_seqwish_replacement(candidate, config),
        ResolutionMethod::Wfmash => build_wfmash_seqwish_replacement(candidate, config),
    }
}

fn candidate_replacement_method(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> ResolutionMethod {
    match config.method {
        ResolutionMethod::Auto => auto_replacement_method(candidate, config),
        ResolutionMethod::Hierarchical => hierarchical_method_by_level(candidate.level),
        ResolutionMethod::IterativeMultiLevel | ResolutionMethod::CoverageMultiBubble => {
            auto_replacement_method(candidate, config)
        }
        ResolutionMethod::ChainPovu => ResolutionMethod::ChainPovu,
        ResolutionMethod::TopFlubbleSweepga => ResolutionMethod::Sweepga,
        method => method,
    }
}

/// Sum of the replacement graph's segment-bp. This is the "output bp" half of
/// the diagnostic replacement compression ratio.
fn replacement_segment_bp(replacement: &Graph) -> usize {
    replacement.segments.iter().map(|s| s.seq.len()).sum()
}

fn replacement_shared_segment_count(replacement: &Graph) -> usize {
    let mut path_support = vec![0usize; replacement.segments.len()];
    let mut last_path_seen = vec![usize::MAX; replacement.segments.len()];
    for (path_idx, path) in replacement.paths.iter().enumerate() {
        for step in &path.steps {
            if step.node < path_support.len() && last_path_seen[step.node] != path_idx {
                path_support[step.node] += 1;
                last_path_seen[step.node] = path_idx;
            }
        }
    }
    path_support.into_iter().filter(|&count| count > 1).count()
}

/// Compression ratio = sum(input traversal bp) / sum(replacement segment bp).
/// Returns `None` when the inputs make the ratio meaningless: zero input bp,
/// zero output bp (cannot divide), or an empty replacement.
fn replacement_compression_ratio(candidate: &BubbleCandidate, replacement: &Graph) -> Option<f64> {
    let input_bp = candidate.traversal_stats.total_len;
    let output_bp = replacement_segment_bp(replacement);
    if input_bp == 0 || output_bp == 0 {
        return None;
    }
    Some(input_bp as f64 / output_bp as f64)
}

fn log_replacement_compression_diagnostics(
    plans: &[ReplacementPlan],
    config: &ResolutionConfig,
    round: usize,
    emit_logs: bool,
) {
    if !emit_logs || plans.is_empty() {
        return;
    }

    let threshold = config.retry_min_compression_ratio;
    let min_input_bp = config.retry_min_input_bp;
    let mut ratios: Vec<f64> = Vec::new();
    let mut output_bp: Vec<usize> = Vec::new();
    let mut below_threshold = 0usize;
    for plan in plans {
        output_bp.push(replacement_segment_bp(&plan.replacement));
        if let Some(ratio) = replacement_compression_ratio(&plan.candidate, &plan.replacement) {
            if threshold > 0.0
                && plan.candidate.traversal_stats.total_len >= min_input_bp
                && ratio < threshold
            {
                below_threshold += 1;
            }
            ratios.push(ratio);
        }
    }

    let threshold_summary = if threshold > 0.0 {
        format!(
            "{threshold:.2} ({} below threshold, min-input-bp={min_input_bp}; diagnostic only)",
            below_threshold
        )
    } else {
        "disabled (diagnostic only)".to_string()
    };
    log::info!(
        "crush round {}: replacement compression diagnostics: ratio {}, output-bp {}; threshold={}; no alternate aligner attempted and no replacement was rejected by ratio",
        round + 1,
        format_ratio_distribution("input/output", &ratios),
        format_bp_distribution("replacement", &output_bp),
        threshold_summary
    );
}

fn format_ratio_distribution(label: &str, values: &[f64]) -> String {
    if values.is_empty() {
        return format!("{label} n=0");
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let pct = |p: usize| -> f64 {
        let idx = ((sorted.len() - 1) * p) / 100;
        sorted[idx]
    };
    let total: f64 = sorted.iter().sum();
    format!(
        "{} n={}, min={:.3}, p50={:.3}, p90={:.3}, max={:.3}, mean={:.3}",
        label,
        sorted.len(),
        sorted[0],
        pct(50),
        pct(90),
        sorted[sorted.len() - 1],
        total / sorted.len() as f64
    )
}

fn auto_replacement_method(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> ResolutionMethod {
    // 3-tier dispatch by median traversal length per docs/crush-architecture-spec.md
    // §Phase-2: median <1 kb → sPOA, 1–10 kb → POASTA, ≥10 kb → sweepga. See
    // `auto_method_by_median` for the routing logic and the rationale for using
    // median (not max) as the decision variable.
    auto_method_by_median(candidate.traversal_stats, config)
}

/// Sequence presented to the replacement aligner for one PathRange. When
/// flanking context is active this is `<left_flank><interior><right_flank>`;
/// otherwise it is the bubble interior alone.
fn range_aligner_sequence(range: &PathRange) -> &[u8] {
    if range.extended_sequence.is_empty() {
        &range.sequence
    } else {
        &range.extended_sequence
    }
}

fn candidate_named_sequences(candidate: &BubbleCandidate) -> (Vec<String>, Vec<(String, Vec<u8>)>) {
    let headers: Vec<String> = candidate
        .ranges
        .iter()
        .enumerate()
        .map(|(i, range)| format!("__impg_bubble_path{}_{}", range.path_idx, i))
        .collect();
    let seqs = headers
        .iter()
        .cloned()
        .zip(
            candidate
                .ranges
                .iter()
                .map(|range| range_aligner_sequence(range).to_vec()),
        )
        .collect();
    (headers, seqs)
}

fn candidate_has_flank(candidate: &BubbleCandidate) -> bool {
    candidate
        .ranges
        .iter()
        .any(|range| !range.extended_sequence.is_empty())
}

fn candidate_named_sequences_longest_first(
    candidate: &BubbleCandidate,
) -> (Vec<String>, Vec<(String, Vec<u8>)>) {
    let (headers, seqs) = candidate_named_sequences(candidate);
    let mut sorted = seqs;
    sorted.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(&b.0)));
    (headers, sorted)
}

fn seqwish_replacement_config(
    config: &ResolutionConfig,
) -> crate::commands::graph::GraphBuildConfig {
    // When the sweepga/seqwish path is in "no-filter" mode, drop every
    // downstream filter — including seqwish's own `min_match_len` exact-match
    // floor. Setting it to 1 keeps every CIGAR `=`/`M` run, so bimodal
    // short-vs-long PAF lines (whose max internal match-run is below the
    // 311 bp default and below the adaptive shortest-traversal clamp) survive.
    let (seqwish_min, adaptive_min_match_len) = if config.sweepga_no_filter {
        (1, false)
    } else {
        match config.replacement_min_match_len_policy {
            ReplacementMinMatchLenPolicy::Adaptive => {
                (config.replacement_seqwish_min_match_len, true)
            }
            ReplacementMinMatchLenPolicy::Fixed(0) => (1, false),
            ReplacementMinMatchLenPolicy::Fixed(value) => (value, false),
        }
    };
    let min_map_length = if config.replacement_min_map_length == 0 {
        seqwish_min
    } else {
        config.replacement_min_map_length
    };
    // `no-filter=true` is the documented "fully unfiltered + many-to-many"
    // request. The post-alignment filter is skipped entirely
    // (`filter_generated_paf` early-returns on `no_filter`), but we also
    // force the mapping / scaffold filter modes to `many:many` and zero
    // the scaffold-mass floor so the requested semantics survive if the
    // filter is ever re-enabled downstream.
    let (num_mappings, scaffold_filter, scaffold_mass) = if config.sweepga_no_filter {
        ("many:many".to_string(), "many:many".to_string(), 0u64)
    } else {
        (
            config.replacement_num_mappings.clone(),
            config.replacement_scaffold_filter.clone(),
            config.replacement_scaffold_mass,
        )
    };
    crate::commands::graph::GraphBuildConfig {
        num_threads: rayon::current_num_threads().max(1),
        show_progress: false,
        min_aln_length: 0,
        min_match_len: seqwish_min,
        adaptive_min_match_len,
        min_map_length,
        min_identity: if config.sweepga_no_filter {
            0.0
        } else {
            config.replacement_min_identity
        },
        input_paf: None,
        no_filter: config.sweepga_no_filter,
        num_mappings,
        scaffold_filter,
        scaffold_mass,
        sparsify: sweepga::knn_graph::SparsificationStrategy::None,
        ..crate::commands::graph::GraphBuildConfig::default()
    }
}

fn replacement_sweepga_kmer_frequency(configured: usize, traversal_count: usize) -> usize {
    if configured > 0 {
        return configured;
    }
    traversal_count
        .saturating_mul(AUTO_SWEEPGA_KMER_FREQUENCY_PER_TRAVERSAL)
        .max(MIN_AUTO_SWEEPGA_KMER_FREQUENCY)
}

/// Maximum FastGA k-mer-frequency cap when `sweepga_no_filter` is set.
/// Picked very high so the cap is effectively disabled while still leaving
/// FastGA's internal handling intact. Real repetitive k-mer densities in the
/// test corpus top out near a few hundred thousand.
const NO_FILTER_SWEEPGA_KMER_FREQUENCY: usize = 1_000_000;

/// Resolve the FastGA k-mer frequency for the sweepga replacement path,
/// honouring `sweepga_no_filter` as a request to disable the auto-cap.
fn resolve_replacement_kmer_frequency(config: &ResolutionConfig, traversal_count: usize) -> usize {
    if config.sweepga_kmer_frequency > 0 {
        return config.sweepga_kmer_frequency;
    }
    if config.sweepga_no_filter {
        return NO_FILTER_SWEEPGA_KMER_FREQUENCY;
    }
    replacement_sweepga_kmer_frequency(0, traversal_count)
}

fn tree_mash_k_schedule(base: usize, count: usize) -> Vec<usize> {
    let base = base.clamp(3, 31);
    let mut values = Vec::with_capacity(count.max(1));
    values.push(base);
    let mut delta = 2usize;
    while values.len() < count.max(1) {
        let high = base.saturating_add(delta);
        if high <= 31 && !values.contains(&high) {
            values.push(high);
            if values.len() >= count {
                break;
            }
        }
        if let Some(low) = base.checked_sub(delta) {
            if low >= 3 && !values.contains(&low) {
                values.push(low);
                if values.len() >= count {
                    break;
                }
            }
        }
        if high > 31 && base < delta + 3 {
            break;
        }
        delta = delta.saturating_add(2);
    }
    values.truncate(count.max(1));
    values
}

fn allwave_pair_schedule(
    sequences: &[allwave::Sequence],
    config: &ResolutionConfig,
) -> Vec<(usize, usize)> {
    let n = sequences.len();
    if n < 2 {
        return Vec::new();
    }
    if n <= 3 {
        return (0..n)
            .flat_map(|i| (0..n).filter(move |&j| i != j).map(move |j| (i, j)))
            .collect();
    }

    let mut pairs = HashSet::new();
    for (tree_idx, mash_k) in tree_mash_k_schedule(config.pair_mash_k, config.pair_tree_count)
        .into_iter()
        .enumerate()
    {
        pairs.extend(allwave::knn_graph::extract_tree_pairs(
            sequences,
            config.pair_k_nearest,
            config.pair_k_farthest,
            0.0,
            mash_k,
        ));
        pairs.extend(salted_random_pairs(
            sequences,
            config.pair_random_fraction,
            tree_idx,
        ));
    }

    let mut ordered: Vec<(usize, usize)> = pairs.into_iter().collect();
    ordered.sort_unstable();
    ordered
}

fn salted_random_pairs(
    sequences: &[allwave::Sequence],
    fraction: f64,
    salt: usize,
) -> Vec<(usize, usize)> {
    if fraction <= 0.0 {
        return Vec::new();
    }
    let n = sequences.len();
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let mut hasher = std::collections::hash_map::DefaultHasher::new();
            salt.hash(&mut hasher);
            sequences[i].id.hash(&mut hasher);
            sequences[j].id.hash(&mut hasher);
            let normalized = hasher.finish() as f64 / u64::MAX as f64;
            if normalized < fraction {
                pairs.push((i, j));
            }
        }
    }
    pairs
}

fn finalize_pairwise_induced_replacement(
    gfa: String,
    headers: &[String],
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
    method: &str,
) -> io::Result<Graph> {
    let debug_dir = std::env::var_os("IMPG_CRUSH_DEBUG_DIR").map(|root| {
        let id = DEBUG_REPLACEMENT_ID.fetch_add(1, Ordering::Relaxed);
        let method = method
            .chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() {
                    c.to_ascii_lowercase()
                } else {
                    '_'
                }
            })
            .collect::<String>();
        let dir = std::path::PathBuf::from(root).join(format!("replacement_{id:04}_{method}"));
        let _ = std::fs::create_dir_all(&dir);
        dir
    });
    if let Some(dir) = &debug_dir {
        let _ = std::fs::write(dir.join("seqwish.gfa"), &gfa);
    }
    let compacted = unchop_gfa(&gfa)?;
    if let Some(dir) = &debug_dir {
        let _ = std::fs::write(dir.join("unchopped.gfa"), &compacted);
    }
    let mut replacement = parse_gfa(&compacted)?;
    order_replacement_paths(&mut replacement, headers)?;
    if config.polish_iterations > 0 {
        match polish_replacement_gfa(&render_graph(&replacement), headers.len(), config) {
            Ok(polished) => {
                let mut polished_replacement = parse_gfa(&polished)?;
                order_replacement_paths(&mut polished_replacement, headers)?;
                validate_replacement_paths(&polished_replacement, candidate, method)?;
                replacement = polished_replacement;
            }
            Err(err) => {
                return Err(io::Error::other(format!(
                    "crush {method}: polish failed; refusing to substitute unpolished replacement: {err}"
                )));
            }
        }
    }
    validate_replacement_paths(&replacement, candidate, method)?;
    Ok(replacement)
}

fn build_chain_povu_smooth_poasta_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    let (headers, seqs) = candidate_named_sequences(candidate);
    if seqs.iter().any(|(_, seq)| seq.is_empty()) {
        return Err(io::Error::other(
            "crush chain-povu: empty traversal present; refusing direct POASTA fallback",
        ));
    }

    let smoothed_attempt = (|| {
        let seed_gfa = linear_candidate_sequences_gfa(&seqs)?;
        let smooth_config = crate::smooth::SmoothConfig {
            n_haps: seqs.len().max(1),
            target_poa_lengths: vec![DEFAULT_CHAIN_POVU_MAX_BLOCK_BP],
            num_threads: 1,
            scoring_params: config.scoring_params,
            poa_padding_fraction: 0.0,
            pre_sorted: false,
            ..crate::smooth::SmoothConfig::new(seqs.len().max(1))
        };
        let smoothed = crate::smooth::smooth_gfa(&seed_gfa, &smooth_config)?;
        let smoothed = unchop_gfa(&smoothed)?;
        let smoothed = strip_full_range_path_names(&smoothed)?;
        let mut poasta_config = config.clone();
        poasta_config.method = ResolutionMethod::Poasta;
        poasta_config.max_iterations = DEFAULT_CHAIN_POVU_POASTA_ITERATIONS;
        poasta_config.polish_iterations = 0;
        poasta_config.retry_min_compression_ratio = 0.0;
        poasta_config.max_traversal_len = DEFAULT_CHAIN_POVU_MAX_BLOCK_BP;
        poasta_config.max_median_traversal_len = DEFAULT_CHAIN_POVU_MAX_BLOCK_BP;

        let smoothed_graph = parse_gfa(&smoothed)?;
        let poasta_polished =
            resolve_graph_bubbles(smoothed_graph, Some(&smoothed), &poasta_config, false)?.gfa;
        let poasta_polished = unchop_gfa(&poasta_polished)?;
        let mut replacement = parse_gfa(&poasta_polished)?;
        order_replacement_paths(&mut replacement, &headers)?;
        validate_replacement_paths(&replacement, candidate, "chain-povu smoothxg→POASTA")?;
        Ok::<Graph, io::Error>(replacement)
    })();

    match smoothed_attempt {
        Ok(smoothed_replacement) => {
            CHAIN_POVU_SMOOTH_USED.fetch_add(1, Ordering::Relaxed);
            log::debug!(
                "crush chain-povu: using smoothxg→POASTA replacement unconditionally after validation ({} segment(s), {} bp)",
                smoothed_replacement.segments.len(),
                replacement_segment_bp(&smoothed_replacement)
            );
            Ok(smoothed_replacement)
        }
        Err(err) => {
            Err(io::Error::other(format!(
                "crush chain-povu: smoothxg->POASTA block path failed validity/build; refusing direct POASTA fallback: {err}"
            )))
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct ChainPovuReplacementDecisions {
    smooth_used: usize,
    direct_on_empty: usize,
    direct_on_smooth_failure: usize,
}

fn reset_chain_povu_replacement_decisions() {
    CHAIN_POVU_SMOOTH_USED.store(0, Ordering::Relaxed);
    CHAIN_POVU_DIRECT_ON_EMPTY.store(0, Ordering::Relaxed);
    CHAIN_POVU_DIRECT_ON_SMOOTH_FAILURE.store(0, Ordering::Relaxed);
}

fn chain_povu_replacement_decisions() -> ChainPovuReplacementDecisions {
    ChainPovuReplacementDecisions {
        smooth_used: CHAIN_POVU_SMOOTH_USED.load(Ordering::Relaxed),
        direct_on_empty: CHAIN_POVU_DIRECT_ON_EMPTY.load(Ordering::Relaxed),
        direct_on_smooth_failure: CHAIN_POVU_DIRECT_ON_SMOOTH_FAILURE.load(Ordering::Relaxed),
    }
}

fn linear_candidate_sequences_gfa(seqs: &[(String, Vec<u8>)]) -> io::Result<String> {
    let mut out = String::new();
    out.push_str("H\tVN:Z:1.0\n");
    for (idx, (_, seq)) in seqs.iter().enumerate() {
        let seq = std::str::from_utf8(seq).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("chain-povu block sequence is not UTF-8 DNA: {err}"),
            )
        })?;
        out.push_str(&format!("S\t{}\t{}\n", idx + 1, seq));
    }
    for (idx, (name, _)) in seqs.iter().enumerate() {
        out.push_str(&format!("P\t{}\t{}+\t*\n", name, idx + 1));
    }
    Ok(out)
}

fn build_allwave_seqwish_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    let (headers, seqs) = candidate_named_sequences(candidate);
    let non_empty = seqs.iter().filter(|(_, seq)| !seq.is_empty()).count();
    if non_empty < 2 {
        return Err(io::Error::other(format!(
            "crush allwave: {} traversal(s), {} non-empty; refusing direct SPOA fallback",
            seqs.len(),
            non_empty
        )));
    }

    let sequences: Vec<allwave::Sequence> = seqs
        .iter()
        .map(|(id, seq)| allwave::Sequence {
            id: id.clone(),
            seq: seq.clone(),
        })
        .collect();

    let (_, mismatch, gap_open, gap_extend, gap_open2, gap_extend2) = config.scoring_params;
    let params = allwave::AlignmentParams {
        match_score: 0,
        mismatch_penalty: mismatch as i32,
        gap_open: gap_open as i32,
        gap_extend: gap_extend as i32,
        gap2_open: Some(gap_open2 as i32),
        gap2_extend: Some(gap_extend2 as i32),
        max_divergence: None,
    };
    let pairs = allwave_pair_schedule(&sequences, config);
    if pairs.is_empty() {
        return Err(io::Error::other(format!(
            "crush allwave: selected zero pair alignments for {} traversal(s); refusing direct SPOA fallback",
            sequences.len()
        )));
    }
    let orientation_params = allwave::AlignmentParams::edit_distance();
    let mut lines: Vec<String> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let alignment = allwave::alignment::align_pair(
                &sequences[i],
                &sequences[j],
                i,
                j,
                &params,
                &orientation_params,
                true,
            );
            allwave::alignment_to_paf(&alignment, &sequences)
        })
        .filter(|line| !line.trim().is_empty())
        .collect();
    lines.sort_unstable();
    lines.dedup();
    let paf = if lines.is_empty() {
        String::new()
    } else {
        let mut out = lines.join("\n");
        out.push('\n');
        out
    };
    log::debug!(
        "crush allwave: {} traversal(s), {} selected pair alignment(s), {} unique PAF line(s), {} PAF byte(s)",
        sequences.len(),
        pairs.len(),
        lines.len(),
        paf.len()
    );

    let graph_config = seqwish_replacement_config(config);
    let gfa = crate::syng_graph::build_gfa_from_paf_and_sequences(&seqs, &paf, &graph_config)?;
    finalize_pairwise_induced_replacement(gfa, &headers, candidate, config, "AllWave/seqwish")
}

fn build_sweepga_seqwish_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    build_sweepga_seqwish_replacement_with_evidence(candidate, config).map(|(graph, _)| graph)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SweepgaSeqwishTailMode {
    Standard,
    TrustSweepgaFilteredPaf,
}

const ZERO_ALIGNMENT_REPLACEMENT_ERROR: &str =
    "SweepGA/seqwish replacement has zero PAF alignment records";

fn zero_alignment_replacement_error(
    aligner: &str,
    traversals: usize,
    aligner_inputs: usize,
    input_bp: usize,
    kmer_frequency: usize,
) -> io::Error {
    io::Error::other(format!(
        "{ZERO_ALIGNMENT_REPLACEMENT_ERROR}: {aligner} emitted no alignments for {traversals} traversal(s), {aligner_inputs} aligner-input traversal(s), {input_bp} input bp, kmer_frequency={kmer_frequency}; refusing to induce an unfolded replacement graph"
    ))
}

fn is_zero_alignment_replacement_error(err: &io::Error) -> bool {
    err.to_string().contains(ZERO_ALIGNMENT_REPLACEMENT_ERROR)
}

#[derive(Clone, Debug, Default)]
struct SweepgaReplacementEvidence {
    alignment_records: usize,
    aligned_bp: u64,
    paf_bytes: usize,
    collinear_records: usize,
    off_diagonal_records: usize,
    dust_like_records: usize,
    kmer_frequency: usize,
    replacement_segments: usize,
    replacement_shared_segments: usize,
    replacement_bp: usize,
}

fn format_sweepga_evidence(evidence: &SweepgaReplacementEvidence) -> String {
    format!(
        "alignment_evidence=raw_paf_records={},aligned_bp={},paf_bytes={},collinear_records={},off_diagonal_records={},dust_like_records={},fastga_frequency={},replacement_segments={},replacement_shared_segments={},replacement_bp={}",
        evidence.alignment_records,
        evidence.aligned_bp,
        evidence.paf_bytes,
        evidence.collinear_records,
        evidence.off_diagonal_records,
        evidence.dust_like_records,
        evidence.kmer_frequency,
        evidence.replacement_segments,
        evidence.replacement_shared_segments,
        evidence.replacement_bp
    )
}

fn build_sweepga_seqwish_replacement_with_evidence(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<(Graph, SweepgaReplacementEvidence)> {
    build_sweepga_seqwish_replacement_with_tail_mode(
        candidate,
        config,
        SweepgaSeqwishTailMode::Standard,
    )
}

fn build_top_flubble_sweepga_replacement_with_evidence(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<(Graph, SweepgaReplacementEvidence)> {
    build_sweepga_seqwish_replacement_with_tail_mode(
        candidate,
        config,
        SweepgaSeqwishTailMode::TrustSweepgaFilteredPaf,
    )
}

fn build_sweepga_seqwish_replacement_with_tail_mode(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
    tail_mode: SweepgaSeqwishTailMode,
) -> io::Result<(Graph, SweepgaReplacementEvidence)> {
    let (headers, seqs) = candidate_named_sequences(candidate);
    let min_aligner_sequence_len = sweepga_backend_min_sequence_len(&config.sweepga_aligner);
    let named: Vec<(String, &[u8])> = seqs
        .iter()
        .filter(|(_, seq)| seq.len() >= min_aligner_sequence_len)
        .map(|(name, seq)| (name.clone(), seq.as_slice()))
        .collect();

    if config.sweepga_sparse_pairs {
        log::info!(
            "crush sweepga: ignoring sparse-pairs for replacement induction; using one all-vs-all self-alignment batch"
        );
    }
    let kmer_frequency = resolve_replacement_kmer_frequency(config, seqs.len());
    let trust_sweepga_filter =
        tail_mode == SweepgaSeqwishTailMode::TrustSweepgaFilteredPaf && !config.sweepga_no_filter;
    let (align_no_filter, num_mappings, scaffold_filter, scaffold_mass) = if trust_sweepga_filter {
        (
            false,
            config.replacement_num_mappings.clone(),
            config.replacement_scaffold_filter.clone(),
            config.replacement_scaffold_mass,
        )
    } else {
        (true, "many:many".to_string(), "many:many".to_string(), 0)
    };
    let make_align_config = |aligner: String| sweepga::library_api::SweepgaAlignConfig {
        num_threads: rayon::current_num_threads().max(1),
        kmer_frequency,
        min_aln_length: config.sweepga_min_aln_length,
        // Normal per-bubble sweepga keeps the historical raw all-vs-all PAF
        // and lets the shared seqwish tail filter it. Top-flubble mode instead
        // trusts SweepGA's own replacement-tier filter and disables the
        // seqwish-tail filter below, avoiding a second filtering pass on large
        // contained regions.
        no_filter: align_no_filter,
        num_mappings: num_mappings.clone(),
        scaffold_filter: scaffold_filter.clone(),
        scaffold_mass,
        sparsify: sweepga::knn_graph::SparsificationStrategy::None,
        map_pct_identity: if aligner.eq_ignore_ascii_case("wfmash") {
            config.sweepga_map_pct_identity.clone()
        } else {
            None
        },
        aligner,
        ..sweepga::library_api::SweepgaAlignConfig::default()
    };
    let mut chunk = String::new();
    if named.len() >= 2 {
        let align_config = make_align_config(config.sweepga_aligner.clone());
        let paf_file =
            sweepga::library_api::sweepga_align(&named, &align_config).map_err(|err| {
                io::Error::other(format!("SweepGA replacement alignment failed: {err}"))
            })?;
        std::fs::File::open(paf_file.path())?.read_to_string(&mut chunk)?;
    }
    let mut paf_lines = chunk
        .lines()
        .filter(|line| !line.is_empty())
        .map(str::to_owned)
        .collect::<Vec<_>>();
    paf_lines.sort_unstable();
    paf_lines.dedup();
    let aligned_bp = paf_lines
        .iter()
        .filter_map(|line| paf_alignment_block_bp(line))
        .sum::<u64>();
    let alignment_profile = classify_paf_alignment_profile(&paf_lines);
    if paf_lines.is_empty() {
        log::info!(
            "crush sweepga: {} backend emitted zero replacement PAF record(s) for {} traversal(s), {} aligner-input traversal(s), {} input bp; skipping replacement induction",
            config.sweepga_aligner,
            seqs.len(),
            named.len(),
            candidate.traversal_stats.total_len
        );
        return Err(zero_alignment_replacement_error(
            &config.sweepga_aligner,
            seqs.len(),
            named.len(),
            candidate.traversal_stats.total_len,
            kmer_frequency,
        ));
    }
    let mut paf = paf_lines.join("\n");
    paf.push('\n');
    log::debug!(
        "crush sweepga: {} traversal(s), {} aligner-input traversal(s), {} unique raw all-vs-all PAF line(s), {} aligned bp, {} PAF byte(s) from {} backend (fastga_frequency={}); seqwish tail filter={}",
        seqs.len(),
        named.len(),
        paf_lines.len(),
        aligned_bp,
        paf.len(),
        config.sweepga_aligner,
        kmer_frequency,
        if config.sweepga_no_filter {
            "disabled (many:many semantics)".to_string()
        } else {
            format!(
                "{} mappings + {} scaffold (mass={})",
                config.replacement_num_mappings,
                config.replacement_scaffold_filter,
                config.replacement_scaffold_mass
            )
        }
    );

    let mut graph_config = seqwish_replacement_config(config);
    if tail_mode == SweepgaSeqwishTailMode::TrustSweepgaFilteredPaf && !config.sweepga_no_filter {
        // `sweepga_align` has already run SweepGA's plane-sweep/scaffold
        // filter. Top-level flubble regions can be much larger than the
        // original per-bubble replacements; running a second seqwish-tail
        // filter here would apply hidden wrapper semantics after SweepGA has
        // already made its documented filtering decision.
        graph_config.no_filter = true;
    }
    let gfa = crate::syng_graph::build_gfa_from_paf_and_sequences(&seqs, &paf, &graph_config)?;
    let replacement =
        finalize_pairwise_induced_replacement(gfa, &headers, candidate, config, "SweepGA/seqwish")?;
    let replacement_shared_segments = replacement_shared_segment_count(&replacement);
    if !paf_lines.is_empty() && replacement_shared_segments == 0 {
        log::warn!(
            "SweepGA/seqwish replacement produced no shared replacement path segments from {} PAF alignment record(s): seqwish output has {} segment(s), {} bp",
            paf_lines.len(),
            replacement.segments.len(),
            replacement_segment_bp(&replacement),
        );
    }
    let evidence = SweepgaReplacementEvidence {
        alignment_records: paf_lines.len(),
        aligned_bp,
        paf_bytes: paf.len(),
        collinear_records: alignment_profile.collinear_records,
        off_diagonal_records: alignment_profile.off_diagonal_records,
        dust_like_records: alignment_profile.dust_like_records,
        kmer_frequency,
        replacement_segments: replacement.segments.len(),
        replacement_shared_segments,
        replacement_bp: replacement_segment_bp(&replacement),
    };
    Ok((replacement, evidence))
}

fn sweepga_backend_min_sequence_len(aligner: &str) -> usize {
    if aligner.eq_ignore_ascii_case("wfmash") {
        // Do not impose an impg-side length floor. wfmash/MashMap owns any
        // minimum-length semantics and should fail explicitly if an input
        // record is unsupported.
        1
    } else {
        // FastGA accepts short non-empty records, but empty FASTA records are
        // not alignable and can trip backend preparation.
        1
    }
}

fn paf_alignment_block_bp(line: &str) -> Option<u64> {
    line.split('\t').nth(10)?.parse().ok()
}

#[derive(Clone, Copy, Debug, Default)]
struct PafAlignmentProfile {
    collinear_records: usize,
    off_diagonal_records: usize,
    dust_like_records: usize,
}

fn classify_paf_alignment_profile(lines: &[String]) -> PafAlignmentProfile {
    let mut profile = PafAlignmentProfile::default();
    for line in lines {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 12 {
            continue;
        }
        let q_len = fields[1].parse::<f64>().unwrap_or(0.0);
        let q_start = fields[2].parse::<f64>().unwrap_or(0.0);
        let q_end = fields[3].parse::<f64>().unwrap_or(0.0);
        let strand = fields[4];
        let t_len = fields[6].parse::<f64>().unwrap_or(0.0);
        let t_start = fields[7].parse::<f64>().unwrap_or(0.0);
        let t_end = fields[8].parse::<f64>().unwrap_or(0.0);
        let matches = fields[9].parse::<f64>().unwrap_or(0.0);
        let block = fields[10].parse::<f64>().unwrap_or(0.0);
        if q_len <= 0.0 || t_len <= 0.0 || block <= 0.0 {
            continue;
        }
        let identity = if block > 0.0 { matches / block } else { 0.0 };
        let q_mid = ((q_start + q_end) / 2.0) / q_len;
        let t_mid_raw = ((t_start + t_end) / 2.0) / t_len;
        let t_mid = if strand == "-" {
            1.0 - t_mid_raw
        } else {
            t_mid_raw
        };
        let diagonal_distance = (q_mid - t_mid).abs();
        let q_frac = (q_end - q_start).max(0.0) / q_len;
        let t_frac = (t_end - t_start).max(0.0) / t_len;
        let short_fraction = q_frac.min(t_frac);
        if block < 100.0 || identity < 0.75 || short_fraction < 0.05 {
            profile.dust_like_records += 1;
        }
        if diagonal_distance <= 0.10 && short_fraction >= 0.10 {
            profile.collinear_records += 1;
        } else if diagonal_distance >= 0.30 {
            profile.off_diagonal_records += 1;
        }
    }
    profile
}

/// Wrapper for `method=wfmash` that forces the SweepGA backend's aligner to
/// `wfmash`, leaving every other knob unchanged. This pins the aligner choice
/// at the method layer so `method=wfmash` is invariant under any future edit of
/// `sweepga_aligner` defaults.
fn build_wfmash_seqwish_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    let mut wfmash_config = config.clone();
    wfmash_config.sweepga_aligner = "wfmash".to_string();
    build_sweepga_seqwish_replacement(candidate, &wfmash_config)
}

fn build_poa_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    let (mut graph, mut engine) = build_global_spoa_engine(config.scoring_params);
    let (headers, sorted_sequences) = candidate_named_sequences_longest_first(candidate);
    let sorted_headers = sorted_sequences
        .iter()
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    let sequences: Vec<String> = sorted_sequences
        .iter()
        .map(|(_, sequence)| String::from_utf8_lossy(sequence).to_string())
        .collect();

    feed_sequences_to_graph(
        &mut engine,
        &mut graph,
        sequences.iter().map(|s| s.as_str()),
    );
    let gfa = graph.generate_gfa(&sorted_headers, false);
    let gfa = unchop_gfa(&gfa)?;
    let mut replacement = parse_gfa(&gfa)?;
    order_replacement_paths(&mut replacement, &headers)?;

    validate_replacement_paths(&replacement, candidate, "SPOA")?;
    Ok(replacement)
}

fn build_chain_greedy_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    build_poasta_replacement(candidate, config)
}

fn build_poasta_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    use poasta::aligner::config::Affine2PieceMinGapCost;
    use poasta::aligner::PoastaAligner;
    use poasta::graphs::poa::POAGraph;
    use poasta::io::graph::graph_to_gfa;

    let (headers, sorted_sequences) = candidate_named_sequences_longest_first(candidate);
    let sorted_headers = sorted_sequences
        .iter()
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    let sequences: Vec<&[u8]> = sorted_sequences
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let weights = sequences
        .iter()
        .map(|sequence| vec![1usize; sequence.len()])
        .collect::<Vec<_>>();
    let mut graph = POAGraph::<u32>::new();

    let scoring = poasta_two_piece_scoring(config.scoring_params)?;
    let mut aligner = PoastaAligner::new(Affine2PieceMinGapCost(scoring), poasta_alignment_type());
    add_poasta_sequences(
        &mut graph,
        &mut aligner,
        &sorted_headers,
        &sequences,
        &weights,
    )?;

    let mut gfa = Vec::new();
    graph_to_gfa(&mut gfa, &graph).map_err(|err| {
        io::Error::other(format!("POASTA replacement GFA generation failed: {err}"))
    })?;
    let gfa = String::from_utf8(gfa).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("POASTA replacement GFA is not UTF-8: {err}"),
        )
    })?;
    let expected_paths = sorted_sequences.into_iter().collect::<Vec<_>>();
    let mut replacement = poasta_gfa_to_exact_graph(&gfa, &expected_paths)?;
    order_replacement_paths(&mut replacement, &headers)?;
    validate_replacement_paths(&replacement, candidate, "POASTA")?;

    Ok(replacement)
}

/// Build an exact path-preserving POASTA graph for already-extracted local
/// sequences and return it as numeric GFA 1.0 `S/L/P` records.
///
/// This is intentionally narrower than full crush replacement: callers own
/// block discovery and lacing, while this helper owns the POASTA invocation and
/// the exact-path clipping/validation used by the crush POASTA backend.
pub fn poasta_sequences_to_gfa(
    headers: &[String],
    sequences: &[Vec<u8>],
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<String> {
    use poasta::aligner::config::Affine2PieceMinGapCost;
    use poasta::aligner::PoastaAligner;
    use poasta::graphs::poa::POAGraph;
    use poasta::io::graph::graph_to_gfa;

    if headers.len() != sequences.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "POASTA block got {} header(s) for {} sequence(s)",
                headers.len(),
                sequences.len()
            ),
        ));
    }
    if headers.is_empty() {
        return Ok(String::new());
    }

    let mut sorted_sequences = headers
        .iter()
        .cloned()
        .zip(sequences.iter().cloned())
        .collect::<Vec<_>>();
    sorted_sequences.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(&b.0)));
    let sorted_headers = sorted_sequences
        .iter()
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    let sequence_refs = sorted_sequences
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect::<Vec<_>>();
    let weights = sequence_refs
        .iter()
        .map(|sequence| vec![1usize; sequence.len()])
        .collect::<Vec<_>>();

    let mut graph = POAGraph::<u32>::new();
    let scoring = poasta_two_piece_scoring(scoring_params)?;
    let mut aligner = PoastaAligner::new(Affine2PieceMinGapCost(scoring), poasta_alignment_type());
    add_poasta_sequences(
        &mut graph,
        &mut aligner,
        &sorted_headers,
        &sequence_refs,
        &weights,
    )?;

    let mut gfa = Vec::new();
    graph_to_gfa(&mut gfa, &graph)
        .map_err(|err| io::Error::other(format!("POASTA block GFA generation failed: {err}")))?;
    let gfa = String::from_utf8(gfa).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("POASTA block GFA is not UTF-8: {err}"),
        )
    })?;

    let expected_paths = sorted_sequences.into_iter().collect::<Vec<_>>();
    let mut replacement = poasta_gfa_to_exact_graph(&gfa, &expected_paths)?;
    order_replacement_paths(&mut replacement, headers)?;
    let expected_in_header_order = headers
        .iter()
        .cloned()
        .zip(sequences.iter().cloned())
        .collect::<Vec<_>>();
    validate_expected_paths(&replacement, &expected_in_header_order, "POASTA block")?;
    Ok(render_numeric_graph(&replacement))
}

fn render_numeric_graph(graph: &Graph) -> String {
    let mut out = String::new();
    out.push_str("H\tVN:Z:1.0\n");
    for (idx, segment) in graph.segments.iter().enumerate() {
        out.push_str(&format!(
            "S\t{}\t{}\n",
            idx + 1,
            String::from_utf8_lossy(&segment.seq)
        ));
    }

    let mut edges = BTreeSet::new();
    for path in &graph.paths {
        for pair in path.steps.windows(2) {
            edges.insert((pair[0].node, pair[0].rev, pair[1].node, pair[1].rev));
        }
    }
    for (from, from_rev, to, to_rev) in edges {
        out.push_str(&format!(
            "L\t{}\t{}\t{}\t{}\t0M\n",
            from + 1,
            if from_rev { "-" } else { "+" },
            to + 1,
            if to_rev { "-" } else { "+" }
        ));
    }

    for path in &graph.paths {
        let steps = path
            .steps
            .iter()
            .map(|step| format!("{}{}", step.node + 1, if step.rev { "-" } else { "+" }))
            .collect::<Vec<_>>()
            .join(",");
        if !steps.is_empty() {
            out.push_str(&format!("P\t{}\t{}\t*\n", path.name, steps));
        }
    }
    out
}

fn poasta_uses_two_piece_affine(scoring_params: (u8, u8, u8, u8, u8, u8)) -> bool {
    let (_, _, _, gap_extend, _, gap_extend2) = scoring_params;
    gap_extend > gap_extend2
}

fn poasta_two_piece_scoring(
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<poasta::aligner::scoring::GapAffine2Piece> {
    let (_, mismatch, gap_open, gap_extend, gap_open2, gap_extend2) = scoring_params;
    if !poasta_uses_two_piece_affine(scoring_params) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "POASTA crush requires two-piece affine scoring with gap_extend1 > gap_extend2; got gap_open1={gap_open}, gap_extend1={gap_extend}, gap_open2={gap_open2}, gap_extend2={gap_extend2}"
            ),
        ));
    }
    Ok(poasta::aligner::scoring::GapAffine2Piece::new(
        mismatch,
        gap_extend,
        gap_open,
        gap_extend2,
        gap_open2,
    ))
}

fn poasta_alignment_type() -> poasta::aligner::scoring::AlignmentType {
    poasta::aligner::scoring::AlignmentType::Global
}

fn poasta_gfa_to_exact_graph(gfa: &str, expected_paths: &[(String, Vec<u8>)]) -> io::Result<Graph> {
    let headers = expected_paths
        .iter()
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    if let Ok(mut direct) = parse_gfa(gfa) {
        if order_replacement_paths(&mut direct, &headers).is_ok()
            && validate_expected_paths(&direct, expected_paths, "POASTA").is_ok()
        {
            return Ok(direct);
        }
    }

    poasta_gfa_walks_to_exact_graph(gfa, expected_paths)
}

#[derive(Clone, Debug)]
struct StarAlignmentRow {
    inserts: Vec<Vec<u8>>,
    bases: Vec<Option<u8>>,
}

fn build_biwfa_inmemory_replacement(
    candidate: &BubbleCandidate,
    _config: &ResolutionConfig,
) -> io::Result<Graph> {
    let headers: Vec<String> = candidate
        .ranges
        .iter()
        .enumerate()
        .map(|(i, range)| format!("__impg_bubble_path{}_{}", range.path_idx, i))
        .collect();
    let root_idx = candidate
        .ranges
        .iter()
        .enumerate()
        .max_by_key(|(_, range)| range.sequence.len())
        .map(|(idx, _)| idx)
        .ok_or_else(|| io::Error::other("BiWFA replacement has no traversals"))?;
    let root = &candidate.ranges[root_idx].sequence;
    let root_len = root.len();

    let rows = candidate
        .ranges
        .iter()
        .enumerate()
        .map(|(idx, range)| {
            if idx == root_idx {
                Ok(root_alignment_row(root))
            } else {
                align_to_root_row(&headers[idx], &range.sequence, &headers[root_idx], root)
            }
        })
        .collect::<io::Result<Vec<_>>>()?;

    let mut replacement = build_star_column_graph(&headers, &rows, root_len);
    let has_empty_path = replacement.paths.iter().any(|path| path.steps.is_empty());
    if !has_empty_path {
        let compacted = unchop_gfa(&render_graph(&replacement))?;
        replacement = parse_gfa(&compacted)?;
        order_replacement_paths(&mut replacement, &headers)?;
    }

    validate_replacement_paths(&replacement, candidate, "BiWFA/in-memory")?;
    Ok(replacement)
}

fn root_alignment_row(root: &[u8]) -> StarAlignmentRow {
    StarAlignmentRow {
        inserts: vec![Vec::new(); root.len() + 1],
        bases: root.iter().map(|&base| Some(base)).collect(),
    }
}

fn align_to_root_row(
    query_name: &str,
    query: &[u8],
    root_name: &str,
    root: &[u8],
) -> io::Result<StarAlignmentRow> {
    if root.is_empty() {
        return Ok(StarAlignmentRow {
            inserts: vec![query.to_vec()],
            bases: Vec::new(),
        });
    }

    if query.is_empty() {
        return Ok(StarAlignmentRow {
            inserts: vec![Vec::new(); root.len() + 1],
            bases: vec![None; root.len()],
        });
    }

    let paf = crate::syng_graph::pairwise_biwfa_paf(query_name, query, root_name, root)
        .ok_or_else(|| io::Error::other("BiWFA root alignment failed"))?;
    let cigar = paf_cigar(&paf)
        .ok_or_else(|| io::Error::other("BiWFA root alignment did not emit cg:Z CIGAR"))?;
    row_from_query_root_cigar(query, root.len(), cigar)
}

fn paf_cigar(paf: &str) -> Option<&str> {
    paf.trim_end()
        .split('\t')
        .find_map(|field| field.strip_prefix("cg:Z:"))
}

fn row_from_query_root_cigar(
    query: &[u8],
    root_len: usize,
    cigar: &str,
) -> io::Result<StarAlignmentRow> {
    let mut row = StarAlignmentRow {
        inserts: vec![Vec::new(); root_len + 1],
        bases: vec![None; root_len],
    };
    let mut q = 0usize;
    let mut r = 0usize;
    let mut len = 0usize;

    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            len = len
                .checked_mul(10)
                .and_then(|v| v.checked_add((ch - b'0') as usize))
                .ok_or_else(|| io::Error::other("BiWFA CIGAR run length overflow"))?;
            continue;
        }
        if len == 0 {
            return Err(io::Error::other("BiWFA CIGAR has zero-length operation"));
        }
        match ch {
            b'M' | b'=' | b'X' => {
                if q + len > query.len() || r + len > root_len {
                    return Err(io::Error::other("BiWFA CIGAR match overconsumes sequence"));
                }
                for _ in 0..len {
                    row.bases[r] = Some(query[q]);
                    q += 1;
                    r += 1;
                }
            }
            b'I' => {
                if q + len > query.len() || r > root_len {
                    return Err(io::Error::other("BiWFA CIGAR insertion overconsumes query"));
                }
                row.inserts[r].extend_from_slice(&query[q..q + len]);
                q += len;
            }
            b'D' => {
                if r + len > root_len {
                    return Err(io::Error::other("BiWFA CIGAR deletion overconsumes root"));
                }
                r += len;
            }
            other => {
                return Err(io::Error::other(format!(
                    "BiWFA CIGAR contains unsupported op '{}'",
                    other as char
                )));
            }
        }
        len = 0;
    }

    if len != 0 {
        return Err(io::Error::other("BiWFA CIGAR ended before operation"));
    }
    if q != query.len() || r != root_len {
        return Err(io::Error::other(format!(
            "BiWFA CIGAR consumed query/root {q}/{r}, expected {}/{}",
            query.len(),
            root_len
        )));
    }
    Ok(row)
}

fn build_star_column_graph(
    headers: &[String],
    rows: &[StarAlignmentRow],
    root_len: usize,
) -> Graph {
    let mut segments = Vec::new();
    let mut node_by_site_allele: FxHashMap<(usize, Vec<u8>), usize> = FxHashMap::default();
    let get_node =
        |site: usize,
         seq: &[u8],
         segments: &mut Vec<Segment>,
         node_by_site_allele: &mut FxHashMap<(usize, Vec<u8>), usize>| {
            let key = (site, seq.to_vec());
            if let Some(&node) = node_by_site_allele.get(&key) {
                return node;
            }
            let node = segments.len();
            segments.push(Segment {
                id: (node + 1).to_string(),
                seq: seq.to_vec(),
            });
            node_by_site_allele.insert(key, node);
            node
        };

    let mut paths = Vec::with_capacity(rows.len());
    for (idx, row) in rows.iter().enumerate() {
        let mut steps = Vec::new();
        for pos in 0..=root_len {
            if !row.inserts[pos].is_empty() {
                let site = pos * 2;
                let node = get_node(
                    site,
                    &row.inserts[pos],
                    &mut segments,
                    &mut node_by_site_allele,
                );
                steps.push(Step { node, rev: false });
            }
            if pos < root_len {
                if let Some(base) = row.bases[pos] {
                    let site = pos * 2 + 1;
                    let node = get_node(site, &[base], &mut segments, &mut node_by_site_allele);
                    steps.push(Step { node, rev: false });
                }
            }
        }
        paths.push(Path {
            name: headers[idx].clone(),
            steps,
        });
    }

    Graph { segments, paths }
}

fn polish_replacement_gfa(
    gfa: &str,
    n_haps: usize,
    config: &ResolutionConfig,
) -> io::Result<String> {
    if config.polish_iterations == 0 {
        return Ok(gfa.to_string());
    }
    match config.polish_method {
        ResolutionPolishMethod::Poa | ResolutionPolishMethod::Poasta => {
            polish_replacement_gfa_with_flubbles(gfa, config)
        }
        ResolutionPolishMethod::Smooth => polish_replacement_gfa_with_smooth(gfa, n_haps, config),
    }
}

fn polish_replacement_gfa_with_flubbles(
    gfa: &str,
    config: &ResolutionConfig,
) -> io::Result<String> {
    let polish_config = ResolutionConfig {
        max_iterations: config.polish_iterations,
        method: match config.polish_method {
            ResolutionPolishMethod::Poa => ResolutionMethod::Poa,
            ResolutionPolishMethod::Poasta => ResolutionMethod::Poasta,
            ResolutionPolishMethod::Smooth => {
                unreachable!("smooth polish is dispatched separately")
            }
        },
        max_bubble_span: 0,
        max_traversal_len: config.polish_max_traversal_len,
        max_median_traversal_len: config.polish_max_median_traversal_len,
        max_total_sequence: config.polish_max_total_sequence,
        max_traversals: config.polish_max_traversals,
        polish_iterations: 0,
        polish_max_traversal_len: config.polish_max_traversal_len,
        polish_max_median_traversal_len: config.polish_max_median_traversal_len,
        polish_max_total_sequence: config.polish_max_total_sequence,
        polish_max_traversals: config.polish_max_traversals,
        polish_method: config.polish_method,
        auto_spoa_max_traversal_len: 0,
        auto_allwave_max_total_sequence: 0,
        auto_allwave_max_traversals: 0,
        scoring_params: config.scoring_params,
        ..ResolutionConfig::default()
    };
    let graph = parse_gfa(gfa)?;
    resolve_graph_bubbles(graph, Some(gfa), &polish_config, false).map(|resolved| resolved.gfa)
}

fn polish_replacement_gfa_with_smooth(
    gfa: &str,
    n_haps: usize,
    config: &ResolutionConfig,
) -> io::Result<String> {
    let before = parse_gfa(gfa)?;
    let before_quality = graph_quality(&before);
    let smooth_config = crate::smooth::SmoothConfig {
        n_haps: n_haps.max(1),
        num_threads: rayon::current_num_threads().max(1),
        scoring_params: config.scoring_params,
        pre_sorted: false,
        ..crate::smooth::SmoothConfig::new(n_haps.max(1))
    };
    let smoothed = crate::smooth::smooth_gfa(gfa, &smooth_config)?;
    let smoothed = strip_full_range_path_names(&smoothed)?;
    let after = parse_gfa(&smoothed)?;
    let after_quality = graph_quality(&after);
    log::debug!(
        "crush polish: accepted smooth pass: before {}; after {}",
        before_quality.summary(),
        after_quality.summary()
    );
    Ok(smoothed)
}

fn strip_full_range_path_names(gfa: &str) -> io::Result<String> {
    let graph = parse_gfa(gfa)?;
    let mut rename = FxHashMap::default();
    for path in &graph.paths {
        let Some((base, start, end)) = split_path_coordinate_suffix(&path.name) else {
            continue;
        };
        let sequence_len = path_sequence(&graph, path)?.len();
        if start == 0 && end == sequence_len {
            rename.insert(path.name.clone(), base.to_string());
        }
    }
    if rename.is_empty() {
        return Ok(gfa.to_string());
    }

    let mut out = String::with_capacity(gfa.len());
    for line in gfa.lines() {
        if let Some(rest) = line.strip_prefix("P\t") {
            let mut fields = rest.splitn(3, '\t');
            if let (Some(name), Some(walk), Some(overlaps)) =
                (fields.next(), fields.next(), fields.next())
            {
                if let Some(new_name) = rename.get(name) {
                    out.push_str("P\t");
                    out.push_str(new_name);
                    out.push('\t');
                    out.push_str(walk);
                    out.push('\t');
                    out.push_str(overlaps);
                    out.push('\n');
                    continue;
                }
            }
        }
        out.push_str(line);
        out.push('\n');
    }
    Ok(out)
}

fn split_path_coordinate_suffix(name: &str) -> Option<(&str, usize, usize)> {
    let (base, range) = name.rsplit_once(':')?;
    let (start, end) = range.split_once('-')?;
    Some((base, start.parse().ok()?, end.parse().ok()?))
}

fn validate_replacement_paths(
    replacement: &Graph,
    candidate: &BubbleCandidate,
    method: &str,
) -> io::Result<()> {
    if replacement.paths.len() != candidate.ranges.len() {
        return Err(io::Error::other(format!(
            "{method} replacement emitted {} paths for {} traversals",
            replacement.paths.len(),
            candidate.ranges.len()
        )));
    }

    for (idx, range) in candidate.ranges.iter().enumerate() {
        let observed = path_sequence(replacement, &replacement.paths[idx])?;
        let expected = range_aligner_sequence(range);
        if observed != expected {
            return Err(io::Error::other(format!(
                "{method} replacement path {} changed sequence length {} -> {}",
                idx,
                expected.len(),
                observed.len()
            )));
        }
    }
    Ok(())
}

/// After the replacement aligner has run on the flank-extended sequences,
/// clip the leading `left_flank_bp` and trailing `right_flank_bp` bp from each
/// replacement path so the resulting Graph encodes only the bubble interior
/// — the same shape `apply_replacement_frontier` expects when the path-step
/// substitution happens between `range.begin_step` and `range.end_step`.
/// Segments whose own bp boundaries fall inside the flank are split into the
/// minimum number of nodes needed to preserve the path topology around the
/// flank/interior boundary; segments fully inside the flank are dropped from
/// the new graph (they remain referenced by no clipped path).
fn clip_replacement_to_interior(
    replacement: Graph,
    candidate: &BubbleCandidate,
) -> io::Result<Graph> {
    if replacement.paths.len() != candidate.ranges.len() {
        return Err(io::Error::other(format!(
            "clip_replacement_to_interior: replacement has {} paths but candidate has {} ranges",
            replacement.paths.len(),
            candidate.ranges.len()
        )));
    }
    let mut output_segments: Vec<Segment> = Vec::new();
    let mut node_by_slice: FxHashMap<(usize, bool, usize, usize), usize> = FxHashMap::default();
    let mut output_paths: Vec<Path> = Vec::with_capacity(replacement.paths.len());

    for (range_idx, range) in candidate.ranges.iter().enumerate() {
        let path = &replacement.paths[range_idx];
        let pieces = replacement_path_pieces(&replacement, path)?;
        let raw_len: usize = pieces.iter().map(|p| p.end - p.start).sum();
        let interior_start = range.left_flank_bp;
        let interior_end = raw_len.checked_sub(range.right_flank_bp).ok_or_else(|| {
            io::Error::other(
                "clip_replacement_to_interior: right flank exceeds replacement path length",
            )
        })?;
        if interior_start > interior_end {
            return Err(io::Error::other(
                "clip_replacement_to_interior: interior bounds inverted",
            ));
        }
        if interior_end - interior_start != range.sequence.len() {
            return Err(io::Error::other(format!(
                "clip_replacement_to_interior: interior span {} bp does not match expected interior {} bp on path {}",
                interior_end - interior_start,
                range.sequence.len(),
                range_idx
            )));
        }
        let mut new_steps = Vec::new();
        for piece in pieces {
            let from = piece.start.max(interior_start);
            let to = piece.end.min(interior_end);
            if from >= to {
                continue;
            }
            let local_start = from - piece.start;
            let local_end = to - piece.start;
            let node = slice_replacement_segment(
                &replacement.segments,
                &mut output_segments,
                &mut node_by_slice,
                piece.node,
                piece.rev,
                local_start,
                local_end,
            )?;
            new_steps.push(Step { node, rev: false });
        }
        output_paths.push(Path {
            name: path.name.clone(),
            steps: new_steps,
        });
    }

    Ok(Graph {
        segments: output_segments,
        paths: output_paths,
    })
}

fn replacement_build_status<T>(result: &io::Result<Option<T>>) -> &'static str {
    match result {
        Ok(Some(_)) => "accepted",
        Ok(None) => "empty",
        Err(err) if replacement_failure_is_path_invalid(err) => "path-invalid",
        Err(_) => "failed",
    }
}

fn replacement_failure_is_path_invalid(err: &io::Error) -> bool {
    let message = err.to_string();
    message.contains("changed sequence length")
        || message.contains("resolved graph failed exact path-sequence validation")
        || message.contains("does not contain the expected traversal sequence")
        || message.contains("replacement graph is missing path")
        || message.contains("replacement emitted")
        || message.contains("clipped replacement emitted")
        || message.contains("clipped replacement path")
        || message.contains("interior span")
}

struct ReplacementPathPiece {
    node: usize,
    rev: bool,
    start: usize,
    end: usize,
}

fn replacement_path_pieces(graph: &Graph, path: &Path) -> io::Result<Vec<ReplacementPathPiece>> {
    let mut offset = 0usize;
    let mut pieces = Vec::with_capacity(path.steps.len());
    for step in &path.steps {
        let segment = graph.segments.get(step.node).ok_or_else(|| {
            io::Error::other(
                "clip_replacement_to_interior: replacement path references an out-of-range segment",
            )
        })?;
        let len = segment.seq.len();
        let start = offset;
        offset = offset
            .checked_add(len)
            .ok_or_else(|| io::Error::other("replacement path length overflow"))?;
        pieces.push(ReplacementPathPiece {
            node: step.node,
            rev: step.rev,
            start,
            end: offset,
        });
    }
    Ok(pieces)
}

fn slice_replacement_segment(
    source_segments: &[Segment],
    output_segments: &mut Vec<Segment>,
    node_by_slice: &mut FxHashMap<(usize, bool, usize, usize), usize>,
    source_node: usize,
    rev: bool,
    local_start: usize,
    local_end: usize,
) -> io::Result<usize> {
    let segment = source_segments.get(source_node).ok_or_else(|| {
        io::Error::other("clip_replacement_to_interior: source segment out of range")
    })?;
    let seg_len = segment.seq.len();
    if local_start > local_end || local_end > seg_len {
        return Err(io::Error::other(
            "clip_replacement_to_interior: slice outside segment bounds",
        ));
    }
    let key = (source_node, rev, local_start, local_end);
    if let Some(&idx) = node_by_slice.get(&key) {
        return Ok(idx);
    }
    let oriented: Vec<u8> = if rev {
        reverse_complement(&segment.seq)
    } else {
        segment.seq.clone()
    };
    let slice = oriented[local_start..local_end].to_vec();
    let id = format!(
        "clip_{}_{}_{}_{}",
        source_node, rev as u8, local_start, local_end
    );
    let idx = output_segments.len();
    output_segments.push(Segment { id, seq: slice });
    node_by_slice.insert(key, idx);
    Ok(idx)
}

fn validate_expected_paths(
    graph: &Graph,
    expected_paths: &[(String, Vec<u8>)],
    method: &str,
) -> io::Result<()> {
    if graph.paths.len() != expected_paths.len() {
        return Err(io::Error::other(format!(
            "{method} replacement emitted {} paths for {} traversals",
            graph.paths.len(),
            expected_paths.len()
        )));
    }

    for (idx, ((name, expected), path)) in expected_paths.iter().zip(&graph.paths).enumerate() {
        let observed = path_sequence(graph, path)?;
        if observed != *expected {
            return Err(io::Error::other(format!(
                "{method} replacement path '{}' at index {} changed sequence length {} -> {}",
                name,
                idx,
                expected.len(),
                observed.len()
            )));
        }
    }
    Ok(())
}

#[derive(Clone, Debug)]
struct PoastaWalk {
    name: String,
    reported_start: Option<usize>,
    reported_end: Option<usize>,
    steps: Vec<Step>,
}

fn poasta_gfa_walks_to_exact_graph(
    gfa: &str,
    expected_paths: &[(String, Vec<u8>)],
) -> io::Result<Graph> {
    let (segments, walks) = parse_poasta_gfa_segments_and_walks(gfa)?;
    let mut walks_by_name = walks
        .into_iter()
        .map(|walk| (walk.name.clone(), walk))
        .collect::<FxHashMap<_, _>>();
    let mut output_segments = Vec::new();
    let mut node_by_slice: FxHashMap<(usize, bool, usize, usize), usize> = FxHashMap::default();
    let mut output_paths = Vec::with_capacity(expected_paths.len());

    for (path_name, expected_sequence) in expected_paths {
        let walk = walks_by_name.remove(path_name).ok_or_else(|| {
            io::Error::other(format!(
                "POASTA replacement graph is missing walk/path '{path_name}'"
            ))
        })?;
        let pieces = poasta_walk_pieces(&segments, &walk)?;
        let raw_sequence = pieces
            .iter()
            .flat_map(|piece| piece.sequence.iter().copied())
            .collect::<Vec<_>>();
        let (clip_start, clip_end) = find_expected_interval(
            &raw_sequence,
            expected_sequence,
            walk.reported_start,
            walk.reported_end,
            path_name,
        )?;
        let mut steps = Vec::new();

        for piece in pieces {
            let from = piece.start.max(clip_start);
            let to = piece.end.min(clip_end);
            if from >= to {
                continue;
            }
            let local_start = from - piece.start;
            let local_end = to - piece.start;
            let node = poasta_slice_node(
                &segments,
                &mut output_segments,
                &mut node_by_slice,
                piece.node,
                piece.rev,
                local_start,
                local_end,
            )?;
            steps.push(Step { node, rev: false });
        }

        output_paths.push(Path {
            name: path_name.clone(),
            steps,
        });
    }

    let graph = Graph {
        segments: output_segments,
        paths: output_paths,
    };
    validate_expected_paths(&graph, expected_paths, "POASTA")?;
    Ok(graph)
}

fn parse_poasta_gfa_segments_and_walks(gfa: &str) -> io::Result<(Vec<Segment>, Vec<PoastaWalk>)> {
    let mut segments = Vec::new();
    let mut id_to_idx: FxHashMap<String, usize> = FxHashMap::default();
    let mut pending_walks = Vec::new();

    for (line_no, line) in gfa.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        match fields[0] {
            "S" => {
                if fields.len() < 3 {
                    return Err(invalid_gfa(line_no, "S line has fewer than 3 fields"));
                }
                let id = fields[1].to_string();
                if id_to_idx.contains_key(&id) {
                    return Err(invalid_gfa(line_no, "duplicate segment ID"));
                }
                id_to_idx.insert(id.clone(), segments.len());
                segments.push(Segment {
                    id,
                    seq: fields[2].as_bytes().to_vec(),
                });
            }
            "W" => {
                if fields.len() < 7 {
                    return Err(invalid_gfa(line_no, "W line has fewer than 7 fields"));
                }
                pending_walks.push((
                    line_no,
                    fields[3].to_string(),
                    fields[4].parse::<usize>().ok(),
                    fields[5].parse::<usize>().ok(),
                    fields[6].to_string(),
                ));
            }
            _ => {}
        }
    }

    let mut walks = Vec::with_capacity(pending_walks.len());
    for (line_no, name, reported_start, reported_end, raw_walk) in pending_walks {
        let mut steps = Vec::new();
        for (id, rev) in parse_w_walk(&raw_walk)
            .ok_or_else(|| invalid_gfa(line_no, "W line has a malformed walk"))?
        {
            let Some(&node) = id_to_idx.get(id) else {
                return Err(invalid_gfa(line_no, "W line references unknown segment"));
            };
            steps.push(Step { node, rev });
        }
        walks.push(PoastaWalk {
            name,
            reported_start,
            reported_end,
            steps,
        });
    }

    Ok((segments, walks))
}

struct PoastaWalkPiece {
    node: usize,
    rev: bool,
    start: usize,
    end: usize,
    sequence: Vec<u8>,
}

fn poasta_walk_pieces(segments: &[Segment], walk: &PoastaWalk) -> io::Result<Vec<PoastaWalkPiece>> {
    let mut offset = 0usize;
    let mut pieces = Vec::with_capacity(walk.steps.len());
    for step in &walk.steps {
        let segment = segments.get(step.node).ok_or_else(|| {
            io::Error::other("POASTA replacement walk references an out-of-range segment")
        })?;
        let sequence = if step.rev {
            reverse_complement(&segment.seq)
        } else {
            segment.seq.clone()
        };
        let start = offset;
        offset = offset
            .checked_add(sequence.len())
            .ok_or_else(|| io::Error::other("POASTA replacement walk length overflow"))?;
        pieces.push(PoastaWalkPiece {
            node: step.node,
            rev: step.rev,
            start,
            end: offset,
            sequence,
        });
    }
    Ok(pieces)
}

fn find_expected_interval(
    raw_sequence: &[u8],
    expected_sequence: &[u8],
    reported_start: Option<usize>,
    reported_end: Option<usize>,
    path_name: &str,
) -> io::Result<(usize, usize)> {
    if expected_sequence.is_empty() {
        let start = reported_start.unwrap_or(0).min(raw_sequence.len());
        return Ok((start, start));
    }

    if let (Some(start), Some(end)) = (reported_start, reported_end) {
        if start <= end
            && end <= raw_sequence.len()
            && raw_sequence.get(start..end) == Some(expected_sequence)
        {
            return Ok((start, end));
        }
    }

    let mut matches = Vec::new();
    if expected_sequence.len() <= raw_sequence.len() {
        for start in 0..=raw_sequence.len() - expected_sequence.len() {
            if &raw_sequence[start..start + expected_sequence.len()] == expected_sequence {
                matches.push(start);
            }
        }
    }
    let Some(start) = matches
        .into_iter()
        .min_by_key(|start| reported_start.map_or(*start, |reported| start.abs_diff(reported)))
    else {
        return Err(io::Error::other(format!(
            "POASTA replacement walk '{path_name}' does not contain the expected traversal sequence (walk length {}, expected length {})",
            raw_sequence.len(),
            expected_sequence.len()
        )));
    };
    Ok((start, start + expected_sequence.len()))
}

fn poasta_slice_node(
    source_segments: &[Segment],
    output_segments: &mut Vec<Segment>,
    node_by_slice: &mut FxHashMap<(usize, bool, usize, usize), usize>,
    source_node: usize,
    rev: bool,
    local_start: usize,
    local_end: usize,
) -> io::Result<usize> {
    let segment = source_segments.get(source_node).ok_or_else(|| {
        io::Error::other("POASTA replacement walk references an out-of-range segment")
    })?;
    let segment_len = segment.seq.len();
    if local_start > local_end || local_end > segment_len {
        return Err(io::Error::other(
            "POASTA replacement clipped walk slice is outside its segment",
        ));
    }
    let key = (source_node, rev, local_start, local_end);
    if let Some(&node) = node_by_slice.get(&key) {
        return Ok(node);
    }

    let oriented = if rev {
        reverse_complement(&segment.seq)
    } else {
        segment.seq.clone()
    };
    let node = output_segments.len();
    output_segments.push(Segment {
        id: format!("poasta_{}", node + 1),
        seq: oriented[local_start..local_end].to_vec(),
    });
    node_by_slice.insert(key, node);
    Ok(node)
}

fn add_poasta_sequences<C>(
    graph: &mut poasta::graphs::poa::POAGraph<u32>,
    aligner: &mut poasta::aligner::PoastaAligner<C>,
    headers: &[String],
    sequences: &[&[u8]],
    weights: &[Vec<usize>],
) -> io::Result<()>
where
    C: poasta::aligner::config::AlignmentConfig,
{
    for (idx, sequence) in sequences.iter().enumerate() {
        if graph.is_empty() {
            graph
                .add_alignment_with_weights(&headers[idx], sequence, None, &weights[idx])
                .map_err(poasta_to_io_error)?;
        } else {
            let result = aligner.align::<u32, _>(graph, sequence);
            graph
                .add_alignment_with_weights(
                    &headers[idx],
                    sequence,
                    Some(&result.alignment),
                    &weights[idx],
                )
                .map_err(poasta_to_io_error)?;
        }
    }
    Ok(())
}

fn poasta_to_io_error(err: poasta::errors::PoastaError) -> io::Error {
    io::Error::other(format!("POASTA replacement failed: {err}"))
}

fn order_replacement_paths(graph: &mut Graph, headers: &[String]) -> io::Result<()> {
    let mut paths_by_name = graph
        .paths
        .drain(..)
        .map(|path| (path.name.clone(), path))
        .collect::<FxHashMap<_, _>>();
    let mut ordered = Vec::with_capacity(headers.len());
    for header in headers {
        let Some(path) = paths_by_name.remove(header) else {
            return Err(io::Error::other(format!(
                "replacement graph is missing path '{header}'"
            )));
        };
        ordered.push(path);
    }
    graph.paths = ordered;
    Ok(())
}

fn path_sequences_equal(before: &Graph, after: &Graph) -> io::Result<bool> {
    let before_map = path_sequence_map(before)?;
    let after_map = path_sequence_map(after)?;
    Ok(before_map == after_map)
}

fn path_sequence_map(graph: &Graph) -> io::Result<FxHashMap<String, Vec<u8>>> {
    let mut map = FxHashMap::default();
    for path in &graph.paths {
        map.insert(path.name.clone(), path_sequence(graph, path)?);
    }
    Ok(map)
}

fn render_rewritten_graph(
    original: &Graph,
    replacements: &[Graph],
    used_original: &FxHashSet<usize>,
    out_paths: &[(String, Vec<OutStep>)],
    next_id: &mut usize,
) -> String {
    let mut out = String::new();
    out.push_str("H\tVN:Z:1.0\n");

    let mut id_by_node: FxHashMap<OutNode, String> = FxHashMap::default();
    let mut used_ids = FxHashSet::<String>::default();
    for node in used_original {
        used_ids.insert(original.segments[*node].id.clone());
    }

    let mut ordered_nodes = Vec::new();
    if replacements.is_empty() {
        let mut seen_original = vec![false; original.segments.len()];
        for (_, steps) in out_paths {
            for step in steps {
                if let OutNode::Original(idx) = step.node {
                    if !seen_original[idx] {
                        seen_original[idx] = true;
                        ordered_nodes.push(step.node);
                    }
                }
            }
        }
    } else {
        let mut seen_original = vec![false; original.segments.len()];
        for (_, steps) in out_paths {
            for step in steps {
                if let OutNode::Original(idx) = step.node {
                    if !seen_original[idx] {
                        seen_original[idx] = true;
                        ordered_nodes.push(step.node);
                    }
                }
            }
        }

        let mut original_position = vec![usize::MAX; original.segments.len()];
        for (position, node) in ordered_nodes.iter().enumerate() {
            if let OutNode::Original(idx) = *node {
                original_position[idx] = position;
            }
        }

        let mut insertions: Vec<(usize, OutNode)> = Vec::new();
        let mut seen_replacement = FxHashSet::<(usize, usize)>::default();
        for (_, steps) in out_paths {
            let mut last_original_position = usize::MAX;
            for step in steps {
                match step.node {
                    OutNode::Original(idx) => {
                        let position = original_position[idx];
                        if position != usize::MAX {
                            last_original_position = position;
                        }
                    }
                    OutNode::Replacement(replacement_idx, replacement_node_idx) => {
                        if seen_replacement.insert((replacement_idx, replacement_node_idx)) {
                            insertions.push((last_original_position, step.node));
                        }
                    }
                }
            }
        }

        insertions.sort_by_key(|&(position, _)| position);
        let original_ordered = std::mem::take(&mut ordered_nodes);
        ordered_nodes.reserve(original_ordered.len() + insertions.len());
        let mut insertion_idx = 0usize;
        for (position, node) in original_ordered.iter().enumerate() {
            ordered_nodes.push(*node);
            while insertion_idx < insertions.len() && insertions[insertion_idx].0 == position {
                ordered_nodes.push(insertions[insertion_idx].1);
                insertion_idx += 1;
            }
        }
        while insertion_idx < insertions.len() {
            ordered_nodes.push(insertions[insertion_idx].1);
            insertion_idx += 1;
        }
    }

    for node in ordered_nodes {
        match node {
            OutNode::Original(original_idx) => {
                let id = original.segments[original_idx].id.clone();
                id_by_node.insert(node, id.clone());
                out.push_str(&format!(
                    "S\t{}\t{}\n",
                    id,
                    String::from_utf8_lossy(&original.segments[original_idx].seq)
                ));
            }
            OutNode::Replacement(replacement_idx, replacement_node_idx) => {
                let replacement = replacements
                    .get(replacement_idx)
                    .expect("used replacement index must refer to an emitted replacement graph");
                let id = next_unused_segment_id(&mut used_ids, next_id);
                id_by_node.insert(node, id.clone());
                out.push_str(&format!(
                    "S\t{}\t{}\n",
                    id,
                    String::from_utf8_lossy(&replacement.segments[replacement_node_idx].seq)
                ));
            }
        }
    }

    let mut edges: BTreeSet<(String, bool, String, bool)> = BTreeSet::new();
    for (_, steps) in out_paths {
        for pair in steps.windows(2) {
            let from = pair[0];
            let to = pair[1];
            if let (Some(from_id), Some(to_id)) =
                (id_by_node.get(&from.node), id_by_node.get(&to.node))
            {
                edges.insert((from_id.clone(), from.rev, to_id.clone(), to.rev));
            }
        }
    }
    for (from, from_rev, to, to_rev) in edges {
        out.push_str(&format!(
            "L\t{}\t{}\t{}\t{}\t0M\n",
            from,
            if from_rev { "-" } else { "+" },
            to,
            if to_rev { "-" } else { "+" }
        ));
    }

    for (name, steps) in out_paths {
        let step_string = steps
            .iter()
            .filter_map(|step| {
                id_by_node
                    .get(&step.node)
                    .map(|id| format!("{}{}", id, if step.rev { "-" } else { "+" }))
            })
            .collect::<Vec<_>>()
            .join(",");
        if !step_string.is_empty() {
            out.push_str(&format!("P\t{}\t{}\t*\n", name, step_string));
        }
    }

    out
}

fn next_unused_segment_id(used_ids: &mut FxHashSet<String>, next_id: &mut usize) -> String {
    loop {
        let id = next_id.to_string();
        *next_id += 1;
        if used_ids.insert(id.clone()) {
            return id;
        }
    }
}

/// Seed for the monotone fresh-id counter at the start of a resolution run.
///
/// Returns one more than the largest integer-encoded segment id present in the
/// input graph (or 1 if no integer ids are present). This guarantees that any
/// id minted by any round of `resolve_graph_bubbles` is strictly greater than
/// every original id at the moment of insertion — the structural form of the
/// Phase 6 fresh-id invariant from `docs/crush-architecture-spec.md`.
fn initial_next_segment_id(graph: &Graph) -> usize {
    let max = graph
        .segments
        .iter()
        .filter_map(|s| s.id.parse::<usize>().ok())
        .max()
        .unwrap_or(0);
    max.saturating_add(1)
}

fn render_graph(graph: &Graph) -> String {
    let out_paths = graph
        .paths
        .iter()
        .map(|path| {
            (
                path.name.clone(),
                path.steps
                    .iter()
                    .map(|step| OutStep {
                        node: OutNode::Original(step.node),
                        rev: step.rev,
                    })
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Vec<_>>();
    let used_original = graph
        .paths
        .iter()
        .flat_map(|path| path.steps.iter().map(|step| step.node))
        .collect::<FxHashSet<_>>();
    // No replacements are emitted here, so next_id is unused; seed at 1 to
    // preserve the previous behaviour for the final-emit path.
    let mut next_id = 1usize;
    render_rewritten_graph(graph, &[], &used_original, &out_paths, &mut next_id)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq_map(gfa: &str) -> FxHashMap<String, String> {
        path_sequences(gfa)
            .unwrap()
            .into_iter()
            .collect::<FxHashMap<_, _>>()
    }

    fn sorted_segment_sequences(gfa: &str) -> Vec<String> {
        let mut sequences = parse_gfa(gfa)
            .unwrap()
            .segments
            .into_iter()
            .map(|segment| String::from_utf8(segment.seq).unwrap())
            .collect::<Vec<_>>();
        sequences.sort();
        sequences
    }

    fn tiny_boundary_graph() -> Graph {
        Graph {
            segments: vec![
                Segment {
                    id: "1".to_string(),
                    seq: b"A".to_vec(),
                },
                Segment {
                    id: "2".to_string(),
                    seq: b"T".to_vec(),
                },
            ],
            paths: vec![Path {
                name: "ref".to_string(),
                steps: vec![
                    Step {
                        node: 0,
                        rev: false,
                    },
                    Step {
                        node: 1,
                        rev: false,
                    },
                ],
            }],
        }
    }

    fn budget_test_window(source: MultiLevelCandidateSource) -> MultiLevelWindowCandidate {
        MultiLevelWindowCandidate {
            candidate: BubbleCandidate {
                ranges: Vec::new(),
                signature: "budget-test".to_string(),
                root_start_step: 0,
                root_end_step: 1,
                root_span: 10_000,
                total_steps: 400,
                unique_steps: 400,
                traversal_stats: TraversalStats {
                    count: 20,
                    min_len: 10_000,
                    median_len: 10_000,
                    p90_len: 10_000,
                    max_len: 10_000,
                    total_len: 200_000,
                },
                level: 0,
            },
            source,
            source_sites: 4,
        }
    }

    fn poasta_resolution_config() -> ResolutionConfig {
        ResolutionConfig {
            method: ResolutionMethod::Poasta,
            max_traversal_len: 100_000,
            max_median_traversal_len: 100_000,
            max_total_sequence: 10_000_000,
            ..ResolutionConfig::default()
        }
    }

    fn assert_poasta_resolves_exact(gfa: &str, config: ResolutionConfig) -> ResolvedGfa {
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(gfa, &config).unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1, "{:?}", resolved.stats);
        assert_eq!(resolved.stats.bailed, 0, "{:?}", resolved.stats);
        resolved
    }

    #[test]
    fn parses_w_walks_as_paths() {
        let gfa = "\
H\tVN:Z:1.1
S\ts0\tAC
S\ts1\tGT
L\ts0\t+\ts1\t-\t0M
W\t*\t0\twalk0\t0\t4\t>s0<s1
";
        let seqs = seq_map(gfa);
        assert_eq!(seqs["walk0"], "ACAC");
    }

    #[test]
    fn poasta_alignment_config_is_global_and_uses_two_piece_scoring_params() {
        use poasta::aligner::scoring::AlignmentCosts;

        let scoring_params = ResolutionConfig::default().scoring_params;
        let (_, mismatch, gap_open, gap_extend, gap_open2, gap_extend2) = scoring_params;
        assert!(matches!(
            poasta_alignment_type(),
            poasta::aligner::scoring::AlignmentType::Global
        ));
        assert!(poasta_uses_two_piece_affine(scoring_params));
        let scoring = poasta_two_piece_scoring(scoring_params).unwrap();
        assert_eq!(scoring.mismatch(), mismatch);
        assert_eq!(scoring.gap_open(), gap_open);
        assert_eq!(scoring.gap_extend(), gap_extend);
        assert_eq!(scoring.gap_open2(), gap_open2);
        assert_eq!(scoring.gap_extend2(), gap_extend2);
        assert!(poasta_uses_two_piece_affine((1, 4, 6, 4, 26, 1)));
        let err = poasta_two_piece_scoring((1, 4, 6, 2, 26, 2)).unwrap_err();
        assert!(err.to_string().contains("requires two-piece affine"));
    }

    #[test]
    fn poasta_global_two_piece_aligner_handles_end_indel_case() {
        use poasta::aligner::config::Affine2PieceMinGapCost;
        use poasta::aligner::PoastaAligner;
        use poasta::graphs::poa::POAGraph;

        let scoring = poasta_two_piece_scoring((1, 4, 6, 4, 26, 1)).unwrap();
        let mut graph = POAGraph::<u32>::new();
        let reference = b"ACGT";
        let weights = vec![1usize; reference.len()];
        graph
            .add_alignment_with_weights("ref", reference, None, &weights)
            .unwrap();

        let aligner = PoastaAligner::new(Affine2PieceMinGapCost(scoring), poasta_alignment_type());
        let query = b"TTACGTAA";
        let result = aligner.align::<u32, _>(&graph, query);

        assert!(!result.alignment.is_empty());
        assert!(
            u32::from(result.score) > 0,
            "global alignment should penalize non-matching query ends"
        );
    }

    #[test]
    fn poasta_gfa_clipped_w_walks_are_normalized_to_exact_paths() {
        let gfa = "\
H\tVN:Z:1.1
S\ts0\tAACCGG
W\t*\t0\tpath0\t2\t4\t>s0
";
        let expected_paths = vec![("path0".to_string(), b"CC".to_vec())];
        let graph = poasta_gfa_to_exact_graph(gfa, &expected_paths).unwrap();
        assert_eq!(graph.paths.len(), 1);
        assert_eq!(
            String::from_utf8(path_sequence(&graph, &graph.paths[0]).unwrap()).unwrap(),
            "CC"
        );
    }

    #[test]
    fn resolves_simple_snp_bubble_without_changing_path_sequences() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAA
S\t2\tC
S\t3\tG
S\t4\tTTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tref\t1+,2+,4+\t*
P\talt\t1+,3+,4+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(gfa, &ResolutionConfig::default()).unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1);
        assert_eq!(resolved.stats.bailed, 0);
    }

    #[test]
    fn resolves_insertion_bubble_without_changing_path_sequences() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAC
S\t2\tGGG
S\t3\tTA
L\t1\t+\t3\t+\t0M
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tref\t1+,3+\t*
P\tins\t1+,2+,3+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(gfa, &ResolutionConfig::default()).unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(before["ref"], "ACTA");
        assert_eq!(before["ins"], "ACGGGTA");
        assert_eq!(resolved.stats.resolved, 1);
    }

    #[test]
    fn resolves_insertion_bubble_with_poasta_without_changing_path_sequences() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAC
S\t2\tGGG
S\t3\tTA
L\t1\t+\t3\t+\t0M
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tref\t1+,3+\t*
P\tins\t1+,2+,3+\t*
";
        assert_poasta_resolves_exact(gfa, poasta_resolution_config());
    }

    #[test]
    fn direct_poa_and_poasta_accept_topology_changing_replacements() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tCTT
S\t3\tCGG
S\t4\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tref\t1+,2+,4+\t*
P\talt\t1+,3+,4+\t*
";
        let before_paths = seq_map(gfa);
        let before_segments = sorted_segment_sequences(gfa);

        for method in [ResolutionMethod::Poa, ResolutionMethod::Poasta] {
            let resolved = resolve_gfa_bubbles(
                gfa,
                &ResolutionConfig {
                    method,
                    max_traversal_len: 100_000,
                    max_median_traversal_len: 100_000,
                    max_total_sequence: 10_000_000,
                    retry_min_compression_ratio: 1_000_000.0,
                    ..ResolutionConfig::default()
                },
            )
            .unwrap();

            assert_eq!(before_paths, seq_map(&resolved.gfa), "{method:?}");
            assert_eq!(
                resolved.stats.resolved, 1,
                "{method:?}: {:?}",
                resolved.stats
            );
            assert_eq!(resolved.stats.bailed, 0, "{method:?}: {:?}", resolved.stats);
            assert_eq!(resolved.stats.retry_attempts, 0, "{method:?}");

            let after_segments = sorted_segment_sequences(&resolved.gfa);
            assert_ne!(
                before_segments, after_segments,
                "{method:?} should be accepted even when it rewrites local segmentation"
            );
            assert!(
                after_segments.iter().any(|seq| seq == "AC"),
                "{method:?} should expose the shared C prefix in the accepted replacement: {:?}",
                after_segments
            );
        }
    }

    #[test]
    fn resolves_deletion_bubble_with_poasta_without_changing_path_sequences() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tGG
S\t2\tCCCC
S\t3\tAA
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t1\t+\t3\t+\t0M
P\tref\t1+,2+,3+\t*
P\tdel\t1+,3+\t*
";
        let before = seq_map(gfa);
        assert_eq!(before["ref"], "GGCCCCAA");
        assert_eq!(before["del"], "GGAA");
        assert_poasta_resolves_exact(gfa, poasta_resolution_config());
    }

    #[test]
    fn resolves_sv_like_bubble_with_poasta_without_changing_path_sequences() {
        let ref_allele = format!(
            "{}{}{}",
            "ACGT".repeat(32),
            "GATTACA".repeat(6),
            "TT".repeat(24)
        );
        let alt_allele = format!(
            "{}{}{}",
            "TGCA".repeat(18),
            "A".repeat(96),
            "CCGG".repeat(9)
        );
        let gfa = format!(
            "\
H\tVN:Z:1.0
S\t1\tGCGC
S\t2\t{ref_allele}
S\t3\t{alt_allele}
S\t4\tTATA
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tref\t1+,2+,4+\t*
P\tsv_alt\t1+,3+,4+\t*
"
        );
        assert_poasta_resolves_exact(&gfa, poasta_resolution_config());
    }

    #[test]
    fn resolves_long_gap_bubble_with_poasta_two_piece_affine() {
        let long_gap = "G".repeat(512);
        let gfa = format!(
            "\
H\tVN:Z:1.0
S\t1\tACGTACGT
S\t2\t{long_gap}
S\t3\tTGCATGCA
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t1\t+\t3\t+\t0M
P\tlong\t1+,2+,3+\t*
P\tshort\t1+,3+\t*
"
        );
        let config = ResolutionConfig {
            scoring_params: (1, 4, 6, 4, 26, 1),
            ..poasta_resolution_config()
        };
        assert!(poasta_uses_two_piece_affine(config.scoring_params));
        assert_poasta_resolves_exact(&gfa, config);
    }

    #[test]
    fn resolves_insertion_bubble_with_star_biwfa_without_changing_path_sequences() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAC
S\t2\tGGG
S\t3\tTA
L\t1\t+\t3\t+\t0M
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tref\t1+,3+\t*
P\tins\t1+,2+,3+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                method: ResolutionMethod::StarBiwfa,
                max_median_traversal_len: 10_000,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1);
    }

    #[test]
    fn resolves_insertion_bubble_with_allwave_without_changing_path_sequences() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAAAAAAAAAAAAAAAAAA
S\t2\tCCCCCCCCCCCCCCCCCCCC
S\t3\tTTTTTTTTTTTTTTTTTTTT
L\t1\t+\t3\t+\t0M
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tref\t1+,3+\t*
P\tins\t1+,2+,3+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                method: ResolutionMethod::Allwave,
                max_median_traversal_len: 10_000,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1);
    }

    #[test]
    fn smooth_polish_preserves_replacement_path_sequences() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tCCCC
S\t3\tGGGG
S\t4\tTTTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tp0\t1+,2+,4+\t*
P\tp1\t1+,3+,4+\t*
";
        let before = seq_map(gfa);
        let polished = polish_replacement_gfa(gfa, 2, &ResolutionConfig::default()).unwrap();
        assert_eq!(before, seq_map(&polished));
    }

    fn routing_candidate_with_stats(stats: TraversalStats) -> BubbleCandidate {
        BubbleCandidate {
            ranges: Vec::new(),
            signature: "route".to_string(),
            root_start_step: 0,
            root_end_step: 1,
            root_span: 0,
            total_steps: 0,
            unique_steps: 0,
            traversal_stats: stats,
            level: 0,
        }
    }

    fn routing_candidate_at_level(level: usize) -> BubbleCandidate {
        BubbleCandidate {
            ranges: Vec::new(),
            signature: format!("route-l{level}"),
            root_start_step: 0,
            root_end_step: 1,
            root_span: 0,
            total_steps: 0,
            unique_steps: 0,
            traversal_stats: TraversalStats::default(),
            level,
        }
    }

    fn traversal_stats_with(median_len: usize, max_len: usize, count: usize) -> TraversalStats {
        TraversalStats {
            count,
            min_len: median_len.min(max_len),
            median_len,
            p90_len: max_len,
            max_len,
            total_len: median_len.saturating_mul(count.max(1)),
        }
    }

    /// Per docs/crush-architecture-spec.md §Phase-2: median <1 kb → sPOA.
    #[test]
    fn auto_routes_small_median_to_spoa() {
        let config = ResolutionConfig::default();
        let candidate = routing_candidate_with_stats(traversal_stats_with(500, 800, 10));
        assert_eq!(
            auto_replacement_method(&candidate, &config),
            ResolutionMethod::Poa
        );
    }

    /// Per docs/crush-architecture-spec.md §Phase-2: 1 kb ≤ median < 10 kb → POASTA.
    #[test]
    fn auto_routes_medium_median_to_poasta() {
        let config = ResolutionConfig::default();
        let candidate = routing_candidate_with_stats(traversal_stats_with(5_000, 8_000, 10));
        assert_eq!(
            auto_replacement_method(&candidate, &config),
            ResolutionMethod::Poasta
        );
    }

    /// Per docs/crush-architecture-spec.md §Phase-2: median ≥10 kb → sweepga.
    #[test]
    fn auto_routes_large_median_to_sweepga() {
        let config = ResolutionConfig::default();
        let candidate = routing_candidate_with_stats(traversal_stats_with(20_000, 30_000, 10));
        assert_eq!(
            auto_replacement_method(&candidate, &config),
            ResolutionMethod::Sweepga
        );
    }

    /// Decision variable is **median**, not max. A bubble with median=100 and
    /// max=30000 (one long outlier among many short traversals) must route to
    /// sPOA — see pick (a) in docs/crush-aligner-failure-trace.md where exactly
    /// this shape (median 165 bp, max 25230 bp, 47 traversals) was mis-routed
    /// to sweepga and produced 46 byte-identical disjoint segments.
    #[test]
    fn auto_routes_by_median_not_max() {
        let config = ResolutionConfig::default();
        let candidate = routing_candidate_with_stats(traversal_stats_with(100, 30_000, 47));
        assert_eq!(
            auto_replacement_method(&candidate, &config),
            ResolutionMethod::Poa,
            "median 100 bp must route to sPOA even with a 30 kb max outlier"
        );
    }

    /// Tier-boundary behaviour: at the exact sPOA threshold the bubble is no
    /// longer "small" and falls through to POASTA; at the POASTA threshold it
    /// falls through to sweepga.
    #[test]
    fn auto_routing_boundaries_are_exclusive_upper_bounds() {
        let config = ResolutionConfig::default();
        let at_spoa_boundary = routing_candidate_with_stats(traversal_stats_with(
            config.auto_spoa_max_traversal_len,
            2_000,
            10,
        ));
        assert_eq!(
            auto_replacement_method(&at_spoa_boundary, &config),
            ResolutionMethod::Poasta,
            "median == auto_spoa_max_traversal_len must route to POASTA, not sPOA"
        );

        let at_poasta_boundary = routing_candidate_with_stats(traversal_stats_with(
            config.auto_poasta_max_traversal_len,
            20_000,
            10,
        ));
        assert_eq!(
            auto_replacement_method(&at_poasta_boundary, &config),
            ResolutionMethod::Sweepga,
            "median == auto_poasta_max_traversal_len must route to sweepga, not POASTA"
        );
    }

    /// Setting `auto_spoa_max_traversal_len = 0` disables the sPOA tier; bubbles
    /// that would have been sPOA-routed fall through to the POASTA tier.
    #[test]
    fn auto_routing_disabling_spoa_tier_falls_through_to_poasta() {
        let config = ResolutionConfig {
            auto_spoa_max_traversal_len: 0,
            ..ResolutionConfig::default()
        };
        let small = routing_candidate_with_stats(traversal_stats_with(100, 200, 10));
        assert_eq!(
            auto_replacement_method(&small, &config),
            ResolutionMethod::Poasta
        );
    }

    /// Setting both `auto_spoa_max_traversal_len = 0` and
    /// `auto_poasta_max_traversal_len = 0` collapses all auto routing to sweepga.
    #[test]
    fn auto_routing_disabling_poasta_tier_falls_through_to_sweepga() {
        let config = ResolutionConfig {
            auto_spoa_max_traversal_len: 0,
            auto_poasta_max_traversal_len: 0,
            ..ResolutionConfig::default()
        };
        let small = routing_candidate_with_stats(traversal_stats_with(100, 200, 10));
        assert_eq!(
            auto_replacement_method(&small, &config),
            ResolutionMethod::Sweepga
        );
    }

    /// `candidate_replacement_method` honours an explicit `config.method` pin —
    /// it does **not** apply the 3-tier routing when the user has asked for a
    /// specific aligner. (The canonical command in docs/c4-crush-handoff.md
    /// uses `method=auto`; this guard preserves the explicit-pin escape hatch.)
    #[test]
    fn explicit_method_pin_bypasses_auto_routing() {
        let candidate = routing_candidate_with_stats(traversal_stats_with(100, 200, 10));
        let sweepga = ResolutionConfig {
            method: ResolutionMethod::Sweepga,
            ..ResolutionConfig::default()
        };
        assert_eq!(
            candidate_replacement_method(&candidate, &sweepga),
            ResolutionMethod::Sweepga,
            "explicit method=sweepga must still pin every bubble to sweepga"
        );

        let poasta = ResolutionConfig {
            method: ResolutionMethod::Poasta,
            ..ResolutionConfig::default()
        };
        let large = routing_candidate_with_stats(traversal_stats_with(20_000, 30_000, 10));
        assert_eq!(
            candidate_replacement_method(&large, &poasta),
            ResolutionMethod::Poasta,
            "explicit method=poasta must still pin every bubble to poasta"
        );
    }

    /// Hierarchical routing: level 0 (top-level) → sweepga, level ≥ 1 (every
    /// sub-bubble) → POASTA. Decoupled from size; sweepga only at the root.
    #[test]
    fn hierarchical_routes_root_to_sweepga_and_subbubbles_to_poasta() {
        let config = ResolutionConfig {
            method: ResolutionMethod::Hierarchical,
            ..ResolutionConfig::default()
        };
        let root = routing_candidate_at_level(0);
        assert_eq!(
            candidate_replacement_method(&root, &config),
            ResolutionMethod::Sweepga,
            "level=0 must route to sweepga regardless of size"
        );
        for level in 1..=3 {
            let child = routing_candidate_at_level(level);
            assert_eq!(
                candidate_replacement_method(&child, &config),
                ResolutionMethod::Poasta,
                "level={level} must route to POASTA regardless of size",
            );
        }
    }

    /// Hierarchical routing ignores median traversal length — the same
    /// candidate routes to sweepga at level 0 and POASTA at level ≥ 1 even
    /// though `auto` would route both to POASTA.
    #[test]
    fn hierarchical_routing_ignores_size() {
        let config = ResolutionConfig {
            method: ResolutionMethod::Hierarchical,
            ..ResolutionConfig::default()
        };
        // A small (median 100 bp) bubble that auto would send to sPOA.
        let mut tiny_root = routing_candidate_with_stats(traversal_stats_with(100, 200, 10));
        tiny_root.level = 0;
        assert_eq!(
            candidate_replacement_method(&tiny_root, &config),
            ResolutionMethod::Sweepga
        );
        let mut tiny_child = routing_candidate_with_stats(traversal_stats_with(100, 200, 10));
        tiny_child.level = 1;
        assert_eq!(
            candidate_replacement_method(&tiny_child, &config),
            ResolutionMethod::Poasta
        );
        // A huge (median 20 kb) bubble that auto would send to sweepga.
        let mut huge_child = routing_candidate_with_stats(traversal_stats_with(20_000, 30_000, 10));
        huge_child.level = 1;
        assert_eq!(
            candidate_replacement_method(&huge_child, &config),
            ResolutionMethod::Poasta,
            "huge bubble at level ≥ 1 must still route to POASTA under hierarchical"
        );
    }

    /// Hierarchical prioritization mirrors auto: level-0 (sweepga) bubbles get
    /// priority 0 and level ≥ 1 (POASTA) bubbles get priority 1, so direct
    /// sub-bubble work is sequenced after the root cut within one round.
    #[test]
    fn hierarchical_prioritizes_root_sweepga_before_subbubble_poasta() {
        let config = ResolutionConfig {
            method: ResolutionMethod::Hierarchical,
            ..ResolutionConfig::default()
        };
        let root = routing_candidate_at_level(0);
        let child = routing_candidate_at_level(1);
        assert_eq!(candidate_selection_priority(&root, &config), 0);
        assert_eq!(candidate_selection_priority(&child, &config), 1);
    }

    #[test]
    fn hierarchical_parse_name_accepts_aliases() {
        assert_eq!(
            ResolutionMethod::parse_name("hierarchical"),
            Some(ResolutionMethod::Hierarchical)
        );
        assert_eq!(
            ResolutionMethod::parse_name("hier"),
            Some(ResolutionMethod::Hierarchical)
        );
    }

    #[test]
    fn chain_povu_parse_name_accepts_aliases() {
        assert_eq!(
            ResolutionMethod::parse_name("chain-povu"),
            Some(ResolutionMethod::ChainPovu)
        );
        assert_eq!(
            ResolutionMethod::parse_name("povu-tree-blocks"),
            Some(ResolutionMethod::ChainPovu)
        );
    }

    #[test]
    fn top_flubble_sweepga_parse_name_accepts_aliases() {
        assert_eq!(
            ResolutionMethod::parse_name("top-flubble-sweepga"),
            Some(ResolutionMethod::TopFlubbleSweepga)
        );
        assert_eq!(
            ResolutionMethod::parse_name("chain-povu-sweepga"),
            Some(ResolutionMethod::TopFlubbleSweepga)
        );
    }

    #[test]
    fn iterative_multi_level_parse_name_and_modes_accept_aliases() {
        assert_eq!(
            ResolutionMethod::parse_name("iterative-multi-level"),
            Some(ResolutionMethod::IterativeMultiLevel)
        );
        assert_eq!(
            ResolutionMethod::parse_name("multi-bubble"),
            Some(ResolutionMethod::IterativeMultiLevel)
        );
        assert_eq!(
            ResolutionMethod::parse_name("multi-flubble"),
            Some(ResolutionMethod::IterativeMultiLevel)
        );
        assert_eq!(
            ResolutionMethod::parse_name("crush-to-size"),
            Some(ResolutionMethod::IterativeMultiLevel)
        );
        assert_eq!(
            MultiLevelWindowMode::parse_name("same-parent"),
            Some(MultiLevelWindowMode::Sibling)
        );
        assert_eq!(
            MultiLevelWindowMode::parse_name("largest"),
            Some(MultiLevelWindowMode::Largest)
        );
        assert_eq!(
            MultiLevelWindowMode::parse_name("one-at-a-time"),
            Some(MultiLevelWindowMode::Largest)
        );
        assert_eq!(
            MultiLevelWindowMode::parse_name("parent"),
            Some(MultiLevelWindowMode::Parent)
        );
        assert_eq!(
            MultiLevelWindowMode::parse_name("frontier"),
            Some(MultiLevelWindowMode::Parent)
        );
        assert_eq!(
            MultiLevelWindowMode::parse_name("path-order"),
            Some(MultiLevelWindowMode::Sliding)
        );
        assert_eq!(
            MultiLevelWindowMode::parse_name("combined"),
            Some(MultiLevelWindowMode::Combined)
        );
        assert_eq!(
            ResolutionMethod::parse_name("coverage-multi-bubble"),
            Some(ResolutionMethod::CoverageMultiBubble)
        );
        assert_eq!(
            MultiLevelObjectiveMode::parse_name("node-path-coverage"),
            Some(MultiLevelObjectiveMode::Coverage)
        );
    }

    #[test]
    fn iterative_multi_level_parent_priority_prefers_balanced_parent_spans() {
        let tiny_root_repeat = BubbleCandidate {
            ranges: Vec::new(),
            signature: "tiny-root-repeat".to_string(),
            root_start_step: 0,
            root_end_step: 1,
            root_span: 98,
            total_steps: 220_000,
            unique_steps: 2_016,
            traversal_stats: TraversalStats {
                count: 312,
                min_len: 98,
                median_len: 26_460,
                p90_len: 32_000,
                max_len: 32_954,
                total_len: 5_385_717,
            },
            level: 0,
        };
        let balanced_parent = BubbleCandidate {
            ranges: Vec::new(),
            signature: "balanced-parent".to_string(),
            root_start_step: 10,
            root_end_step: 100,
            root_span: 32_787,
            total_steps: 91_000,
            unique_steps: 1_696,
            traversal_stats: TraversalStats {
                count: 117,
                min_len: 100,
                median_len: 26_418,
                p90_len: 32_000,
                max_len: 32_787,
                total_len: 2_363_564,
            },
            level: 0,
        };

        assert!(
            multi_level_parent_frontier_span_score(&balanced_parent)
                > multi_level_parent_frontier_span_score(&tiny_root_repeat)
        );

        let mut candidates = [
            MultiLevelWindowCandidate {
                candidate: tiny_root_repeat,
                source: MultiLevelCandidateSource::TopLevel,
                source_sites: 1,
            },
            MultiLevelWindowCandidate {
                candidate: balanced_parent,
                source: MultiLevelCandidateSource::TopLevel,
                source_sites: 1,
            },
        ];
        candidates.sort_by(|a, b| {
            compare_same_source_window_priority(a, b, MultiLevelWindowMode::Parent)
        });

        assert_eq!(candidates[0].candidate.root_span, 32_787);
    }

    #[test]
    fn iterative_multi_level_largest_mode_forces_one_parent_frontier_candidate() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tGGGG
S\t7\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t7\t+\t0M
L\t4\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,4+,5+,7+\t*
P\tleft_alt\t1+,3+,4+,5+,7+\t*
P\tright_alt\t1+,2+,4+,6+,7+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let config = ResolutionConfig {
            method: ResolutionMethod::IterativeMultiLevel,
            multi_level_window_mode: MultiLevelWindowMode::Largest,
            multi_level_window_target_bp: 100,
            multi_level_max_window_sites: 3,
            multi_level_candidate_limit: 32,
            ..ResolutionConfig::default()
        };
        let (discovered, _, _) = discover_all_candidates(&graph, &config, false).unwrap();
        let generated =
            generate_multi_level_candidates(&graph, &config, &discovered, &FxHashSet::default());
        assert_eq!(generated.candidates.len(), 1, "{generated:?}");
        assert!(matches!(
            generated.candidates[0].source,
            MultiLevelCandidateSource::TopLevel | MultiLevelCandidateSource::ParentDescendants
        ));
        assert_eq!(effective_multi_level_candidate_limit(&config), 1);
    }

    #[test]
    fn coverage_objective_accepts_coverage_or_singleton_improvements() {
        let before = GraphQuality {
            node_coverage_bp_weighted_mean: 2.0,
            singleton_bp: 100,
            ..GraphQuality::default()
        };
        let singleton_better = GraphQuality {
            node_coverage_bp_weighted_mean: 2.0,
            singleton_bp: 50,
            ..GraphQuality::default()
        };
        let coverage_better = GraphQuality {
            node_coverage_bp_weighted_mean: 2.1,
            singleton_bp: 120,
            ..GraphQuality::default()
        };
        let both_worse = GraphQuality {
            node_coverage_bp_weighted_mean: 1.9,
            singleton_bp: 120,
            ..GraphQuality::default()
        };

        assert!(
            graph_objective_delta_for_mode(
                before,
                singleton_better,
                MultiLevelObjectiveMode::Coverage
            ) > 0
        );
        assert!(
            graph_objective_delta_for_mode(
                before,
                coverage_better,
                MultiLevelObjectiveMode::Coverage
            ) > 0
        );
        assert!(
            graph_objective_delta_for_mode(before, both_worse, MultiLevelObjectiveMode::Coverage)
                < 0
        );
    }

    #[test]
    fn repeat_aware_boundaries_reject_tiny_high_frequency_low_complexity_anchors() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAAAAAAAA
S\t2\tTTTTTTTTTT
L\t1\t+\t2\t+\t0M
P\tp1\t1+,2+\t*
P\tp2\t1+,2+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let visits = graph_node_visit_counts(&graph);
        assert!(repeat_boundary_should_reject(
            &graph,
            Step {
                node: 0,
                rev: false
            },
            Step {
                node: 1,
                rev: false
            },
            &visits
        ));
    }

    #[test]
    fn wfmash_aligner_input_floor_does_not_hide_tool_length_semantics() {
        assert_eq!(sweepga_backend_min_sequence_len("wfmash"), 1);
        assert_eq!(sweepga_backend_min_sequence_len("fastga"), 1);
    }

    #[test]
    fn chain_povu_selects_parent_subtree_under_cap() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tAA
S\t6\tA
S\t7\tCCC
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t2\t+\t5\t+\t0M
L\t5\t+\t4\t+\t0M
L\t4\t+\t6\t+\t0M
L\t1\t+\t7\t+\t0M
L\t7\t+\t6\t+\t0M
P\tref\t1+,2+,3+,4+,6+\t*
P\tinner_alt\t1+,2+,5+,4+,6+\t*
P\touter_alt\t1+,7+,6+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let config = ResolutionConfig {
            method: ResolutionMethod::ChainPovu,
            ..ResolutionConfig::default()
        };
        let (discovered, _, _) = discover_all_candidates(&graph, &config, false).unwrap();
        let selection = select_chain_povu_blocks(
            &discovered,
            DEFAULT_CHAIN_POVU_MAX_BLOCK_BP,
            &FxHashSet::default(),
        );
        assert_eq!(selection.selected.len(), 1, "{selection:?}");
        assert!(
            selection.merged_descendant_sites > 0,
            "chain-povu should absorb nested descendants into the parent block"
        );
    }

    #[test]
    fn top_flubble_sweepga_selects_only_level_zero_roots() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tAA
S\t6\tA
S\t7\tCCC
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t2\t+\t5\t+\t0M
L\t5\t+\t4\t+\t0M
L\t4\t+\t6\t+\t0M
L\t1\t+\t7\t+\t0M
L\t7\t+\t6\t+\t0M
P\tref\t1+,2+,3+,4+,6+\t*
P\tinner_alt\t1+,2+,5+,4+,6+\t*
P\touter_alt\t1+,7+,6+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let config = ResolutionConfig {
            method: ResolutionMethod::TopFlubbleSweepga,
            ..ResolutionConfig::default()
        };
        let (discovered, _, _) = discover_all_candidates(&graph, &config, false).unwrap();
        assert!(
            discovered.iter().any(|d| d.povu_level > 0),
            "fixture must exercise nested POVU sites"
        );
        let selection = select_top_flubble_sweepga_regions(&discovered, &FxHashSet::default());
        assert_eq!(selection.top_level_regions, selection.selected.len());
        assert!(
            selection.descendant_sites > 0,
            "nested sites should be counted but not selected as separate SweepGA regions"
        );
        assert!(selection
            .selected
            .iter()
            .all(|candidate| candidate.level == 0));
    }

    #[test]
    fn top_flubble_sweepga_routes_replacement_to_sweepga() {
        let config = ResolutionConfig {
            method: ResolutionMethod::TopFlubbleSweepga,
            ..ResolutionConfig::default()
        };
        let root = routing_candidate_at_level(0);
        assert_eq!(
            candidate_replacement_method(&root, &config),
            ResolutionMethod::Sweepga
        );
    }

    #[test]
    fn chain_povu_resolves_nested_parent_as_one_block() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tAA
S\t6\tA
S\t7\tCCC
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t2\t+\t5\t+\t0M
L\t5\t+\t4\t+\t0M
L\t4\t+\t6\t+\t0M
L\t1\t+\t7\t+\t0M
L\t7\t+\t6\t+\t0M
P\tref\t1+,2+,3+,4+,6+\t*
P\tinner_alt\t1+,2+,5+,4+,6+\t*
P\touter_alt\t1+,7+,6+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                max_iterations: 1,
                method: ResolutionMethod::ChainPovu,
                retry_min_compression_ratio: 0.0,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1, "{:?}", resolved.stats);
        assert_eq!(resolved.stats.bailed, 0, "{:?}", resolved.stats);
    }

    /// Auto-mode priority places direct-tier picks (sPOA, POASTA) ahead of
    /// sweepga so high-confidence local condensations are tried first within
    /// one frontier round.
    #[test]
    fn auto_prioritizes_direct_candidates_before_sweepga() {
        let config = ResolutionConfig::default();
        let small = routing_candidate_with_stats(traversal_stats_with(500, 800, 10));
        let medium = routing_candidate_with_stats(traversal_stats_with(5_000, 8_000, 10));
        let large = routing_candidate_with_stats(traversal_stats_with(20_000, 30_000, 10));
        assert_eq!(candidate_selection_priority(&small, &config), 1);
        assert_eq!(candidate_selection_priority(&medium, &config), 1);
        assert_eq!(candidate_selection_priority(&large, &config), 0);

        // Explicit-method pin disables auto prioritization entirely.
        let explicit_sweepga = ResolutionConfig {
            method: ResolutionMethod::Sweepga,
            ..config
        };
        assert_eq!(candidate_selection_priority(&small, &explicit_sweepga), 0);
    }

    #[test]
    fn multi_tree_mash_schedule_perturbs_k_deterministically() {
        assert_eq!(tree_mash_k_schedule(15, 1), vec![15]);
        assert_eq!(tree_mash_k_schedule(15, 5), vec![15, 17, 13, 19, 11]);
        assert_eq!(tree_mash_k_schedule(3, 3), vec![3, 5, 7]);
        assert_eq!(tree_mash_k_schedule(31, 3), vec![31, 29, 27]);
    }

    #[test]
    fn salted_random_pair_sampling_changes_between_trees() {
        let sequences = (0..8)
            .map(|idx| allwave::Sequence {
                id: format!("seq{idx}"),
                seq: b"ACGTACGTACGTACGT".to_vec(),
            })
            .collect::<Vec<_>>();
        let first = salted_random_pairs(&sequences, 0.25, 0);
        let second = salted_random_pairs(&sequences, 0.25, 1);
        assert_ne!(first, second);
        assert!(first.iter().all(|(i, j)| i != j));
        assert!(second.iter().all(|(i, j)| i != j));
    }

    #[test]
    fn graph_quality_penalizes_path_white_space_bridges() {
        let spacer = "A".repeat(1024);
        let compact = parse_gfa(&format!(
            "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\t{spacer}
S\t3\tGGGG
P\tp\t1+,2+,3+\t*
"
        ))
        .unwrap();
        let bridged = parse_gfa(&format!(
            "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\t{spacer}
S\t3\tGGGG
P\tp\t1+,3+\t*
"
        ))
        .unwrap();
        let compact_quality = graph_quality(&compact);
        let bridged_quality = graph_quality(&bridged);
        assert_eq!(compact_quality.path_white_space_bp_total, 0);
        assert_eq!(bridged_quality.path_white_space_bp_total, 1024);
        assert!(bridged_quality.score > compact_quality.score);
    }

    #[test]
    fn first_crush_round_accepts_path_step_growth_without_visual_tail_regression() {
        let seq1 = "ACGT".repeat(25);
        let seq2 = "TGCA".repeat(25);
        let gfa = format!(
            "\
H\tVN:Z:1.0
S\t1\tC
S\t2\t{seq1}
S\t3\t{seq2}
S\t4\tG
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tref\t1+,2+,4+\t*
P\talt\t1+,3+,4+\t*
"
        );
        let before = seq_map(&gfa);
        let resolved = resolve_gfa_bubbles(
            &gfa,
            &ResolutionConfig {
                method: ResolutionMethod::StarBiwfa,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();

        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1);
    }

    #[test]
    fn replacement_seqwish_filter_defaults_to_local_min_match_off() {
        let mut config = ResolutionConfig {
            replacement_seqwish_min_match_len: 311,
            replacement_min_map_length: 0,
            replacement_min_identity: 0.97,
            ..ResolutionConfig::default()
        };
        let graph_config = seqwish_replacement_config(&config);
        assert_eq!(graph_config.min_match_len, 1);
        assert!(!graph_config.adaptive_min_match_len);
        assert_eq!(graph_config.min_map_length, 1);
        assert!((graph_config.min_identity - 0.97).abs() < f64::EPSILON);
        assert_eq!(graph_config.num_mappings, "1:1");
        assert_eq!(graph_config.scaffold_filter, "1:1");

        config.replacement_min_map_length = 500;
        config.replacement_num_mappings = "1:many".to_string();
        config.replacement_scaffold_filter = "many:many".to_string();
        let graph_config = seqwish_replacement_config(&config);
        assert_eq!(graph_config.min_match_len, 1);
        assert!(!graph_config.adaptive_min_match_len);
        assert_eq!(graph_config.min_map_length, 500);
        assert_eq!(graph_config.num_mappings, "1:many");
        assert_eq!(graph_config.scaffold_filter, "many:many");
    }

    #[test]
    fn adaptive_min_match_policy_is_opt_in() {
        let config = ResolutionConfig {
            replacement_seqwish_min_match_len: 311,
            replacement_min_match_len_policy: ReplacementMinMatchLenPolicy::Adaptive,
            replacement_min_map_length: 0,
            ..ResolutionConfig::default()
        };
        let graph_config = seqwish_replacement_config(&config);
        assert_eq!(graph_config.min_match_len, 311);
        assert_eq!(graph_config.min_map_length, 311);
        assert!(graph_config.adaptive_min_match_len);
    }

    #[test]
    fn replacement_min_match_policy_can_be_fixed_or_disabled() {
        let config = ResolutionConfig {
            replacement_seqwish_min_match_len: 311,
            replacement_min_match_len_policy: ReplacementMinMatchLenPolicy::Fixed(63),
            replacement_min_map_length: 0,
            ..ResolutionConfig::default()
        };
        let graph_config = seqwish_replacement_config(&config);
        assert_eq!(graph_config.min_match_len, 63);
        assert_eq!(graph_config.min_map_length, 63);
        assert!(!graph_config.adaptive_min_match_len);

        let config = ResolutionConfig {
            replacement_min_match_len_policy: ReplacementMinMatchLenPolicy::Fixed(0),
            replacement_min_map_length: 0,
            ..ResolutionConfig::default()
        };
        let graph_config = seqwish_replacement_config(&config);
        assert_eq!(graph_config.min_match_len, 1);
        assert_eq!(graph_config.min_map_length, 1);
        assert!(!graph_config.adaptive_min_match_len);
    }

    #[test]
    fn sweepga_replacement_frequency_auto_keeps_local_repeats_seedable() {
        assert_eq!(replacement_sweepga_kmer_frequency(42, 463), 42);
        assert_eq!(replacement_sweepga_kmer_frequency(0, 2), 1_000);
        assert_eq!(replacement_sweepga_kmer_frequency(0, 463), 4_630);
        assert_eq!(replacement_sweepga_kmer_frequency(0, 10_000), 100_000);
    }

    #[test]
    fn no_filter_disables_seqwish_min_match_and_identity_filters() {
        let config = ResolutionConfig {
            replacement_seqwish_min_match_len: 311,
            replacement_min_map_length: 0,
            replacement_min_identity: 0.97,
            sweepga_no_filter: true,
            ..ResolutionConfig::default()
        };
        let graph_config = seqwish_replacement_config(&config);
        assert_eq!(graph_config.min_match_len, 1);
        assert_eq!(graph_config.min_map_length, 1);
        assert!(!graph_config.adaptive_min_match_len);
        assert!(graph_config.min_identity.abs() < f64::EPSILON);
        assert!(graph_config.no_filter);
    }

    #[test]
    fn replacement_scaffold_mass_defaults_to_zero_and_propagates() {
        // Default ResolutionConfig should leave the sweepga scaffold-mass
        // filter disabled (0). See docs/crush-scaffold-mass-zero.md.
        let config = ResolutionConfig::default();
        assert_eq!(config.replacement_scaffold_mass, 0);
        let graph_config = seqwish_replacement_config(&config);
        assert_eq!(graph_config.scaffold_mass, 0);
    }

    #[test]
    fn replacement_scaffold_mass_user_override_propagates() {
        let config = ResolutionConfig {
            replacement_scaffold_mass: 500,
            ..ResolutionConfig::default()
        };
        let graph_config = seqwish_replacement_config(&config);
        assert_eq!(graph_config.scaffold_mass, 500);
    }

    #[test]
    fn no_filter_overrides_kmer_frequency_auto_cap() {
        let mut config = ResolutionConfig {
            sweepga_no_filter: true,
            sweepga_kmer_frequency: 0,
            ..ResolutionConfig::default()
        };
        assert_eq!(
            resolve_replacement_kmer_frequency(&config, 312),
            NO_FILTER_SWEEPGA_KMER_FREQUENCY
        );
        config.sweepga_kmer_frequency = 4242;
        assert_eq!(resolve_replacement_kmer_frequency(&config, 312), 4242);
        config.sweepga_kmer_frequency = 0;
        config.sweepga_no_filter = false;
        assert_eq!(resolve_replacement_kmer_frequency(&config, 312), 3_120);
    }

    #[test]
    fn allwave_processes_candidate_regardless_of_pair_alignment_cap() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAAAAAAAAAAAAAAAAAA
S\t2\tCCCCCCCCCCCCCCCCCCCC
S\t3\tTTTTTTTTTTTTTTTTTTTT
L\t1\t+\t3\t+\t0M
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tref\t1+,3+\t*
P\tins\t1+,2+,3+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                method: ResolutionMethod::Allwave,
                max_pair_alignments: 1,
                max_median_traversal_len: 10_000,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1, "{:?}", resolved.stats);
        assert_eq!(resolved.stats.bailed, 0, "{:?}", resolved.stats);
    }

    #[test]
    fn auto_resolves_indel_ladder_without_changing_path_sequences() {
        let long_ref = "ACGT".repeat(80);
        let long_alt = format!(
            "{}{}{}",
            "ACGT".repeat(35),
            "TTGGAACCGGTT",
            "ACGT".repeat(42)
        );
        let gfa = format!(
            "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tC
S\t3\tG
S\t4\tTTTT
S\t5\tGGGGGGGG
S\t6\tCCCC
S\t7\t{long_ref}
S\t8\t{long_alt}
S\t9\tTATA
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t6\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
L\t7\t+\t9\t+\t0M
L\t6\t+\t8\t+\t0M
L\t8\t+\t9\t+\t0M
P\tref\t1+,2+,4+,6+,7+,9+\t*
P\tsnp\t1+,3+,4+,6+,7+,9+\t*
P\tins\t1+,2+,4+,5+,6+,7+,9+\t*
P\tlong_alt\t1+,2+,4+,6+,8+,9+\t*
P\tcombo\t1+,3+,4+,5+,6+,8+,9+\t*
"
        );
        let before = seq_map(&gfa);
        let resolved = resolve_gfa_bubbles(&gfa, &ResolutionConfig::default()).unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert!(resolved.stats.resolved >= 1, "{:?}", resolved.stats);
        assert_eq!(resolved.stats.bailed, 0, "{:?}", resolved.stats);
    }

    #[test]
    fn skips_sequence_identical_bubbles_after_deferred_materialization() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tC
S\t4\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tref\t1+,2+,4+\t*
P\talt\t1+,3+,4+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(gfa, &ResolutionConfig::default()).unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 0);
        assert_eq!(resolved.stats.bailed, 0);
    }

    #[test]
    fn deferred_identity_skip_ignores_flank_context_differences() {
        let graph = Graph {
            segments: vec![
                Segment {
                    id: "l_ref".to_string(),
                    seq: b"G".to_vec(),
                },
                Segment {
                    id: "l_alt".to_string(),
                    seq: b"A".to_vec(),
                },
                Segment {
                    id: "m_ref".to_string(),
                    seq: b"C".to_vec(),
                },
                Segment {
                    id: "m_alt".to_string(),
                    seq: b"C".to_vec(),
                },
                Segment {
                    id: "r_ref".to_string(),
                    seq: b"T".to_vec(),
                },
                Segment {
                    id: "r_alt".to_string(),
                    seq: b"G".to_vec(),
                },
            ],
            paths: vec![
                Path {
                    name: "ref".to_string(),
                    steps: vec![
                        Step {
                            node: 0,
                            rev: false,
                        },
                        Step {
                            node: 2,
                            rev: false,
                        },
                        Step {
                            node: 4,
                            rev: false,
                        },
                    ],
                },
                Path {
                    name: "alt".to_string(),
                    steps: vec![
                        Step {
                            node: 1,
                            rev: false,
                        },
                        Step {
                            node: 3,
                            rev: false,
                        },
                        Step {
                            node: 5,
                            rev: false,
                        },
                    ],
                },
            ],
        };
        let mut candidate = BubbleCandidate {
            ranges: vec![
                PathRange {
                    path_idx: 0,
                    begin_step: 1,
                    end_step: 2,
                    ..PathRange::default()
                },
                PathRange {
                    path_idx: 1,
                    begin_step: 1,
                    end_step: 2,
                    ..PathRange::default()
                },
            ],
            signature: "manual-identical-interior".to_string(),
            root_start_step: 1,
            root_end_step: 2,
            root_span: 1,
            total_steps: 2,
            unique_steps: 2,
            traversal_stats: TraversalStats {
                count: 2,
                min_len: 1,
                median_len: 1,
                p90_len: 1,
                max_len: 1,
                total_len: 2,
            },
            level: 0,
        };
        let should_build = materialize_candidate_sequences(
            &graph,
            &mut candidate,
            &ResolutionConfig {
                replacement_flank_bp: 1,
                ..ResolutionConfig::default()
            },
        );

        assert!(!should_build);
        assert_eq!(candidate.ranges[0].sequence, b"C");
        assert_eq!(candidate.ranges[1].sequence, b"C");
        assert_eq!(candidate.ranges[0].extended_sequence, b"GCT");
        assert_eq!(candidate.ranges[1].extended_sequence, b"ACG");
        assert_ne!(
            candidate.ranges[0].extended_sequence,
            candidate.ranges[1].extended_sequence
        );
    }

    #[test]
    fn direct_poa_processes_candidate_regardless_of_median_traversal_budget() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAC
S\t2\tGGG
S\t3\tTA
L\t1\t+\t3\t+\t0M
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tref\t1+,3+\t*
P\tins\t1+,2+,3+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                method: ResolutionMethod::Poa,
                max_median_traversal_len: 4,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1, "{:?}", resolved.stats);
        assert_eq!(resolved.stats.bailed, 0, "{:?}", resolved.stats);
    }

    #[test]
    fn resolves_maximal_eligible_nested_site_first() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tAA
S\t6\tA
S\t7\tCCC
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t2\t+\t5\t+\t0M
L\t5\t+\t4\t+\t0M
L\t4\t+\t6\t+\t0M
L\t1\t+\t7\t+\t0M
L\t7\t+\t6\t+\t0M
P\tref\t1+,2+,3+,4+,6+\t*
P\tinner_alt\t1+,2+,5+,4+,6+\t*
P\touter_alt\t1+,7+,6+\t*
";
        let before = seq_map(gfa);
        let graph = povu::NativeGfa::parse(gfa).unwrap();
        let decomposition = graph.decompose_flubbles(&["ref".to_string()]).unwrap();
        let leaves = decomposition.leaf_sites_bottom_up();
        assert_eq!(
            leaves
                .iter()
                .map(|site| (site.reference_start_step, site.reference_end_step))
                .collect::<Vec<_>>(),
            vec![(1, 3)]
        );
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                max_iterations: 1,
                method: ResolutionMethod::Poa,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1);
        assert_eq!(resolved.stats.iterations, 1);
        assert_eq!(resolved.stats.bailed, 0);
    }

    #[test]
    fn resolves_independent_sites_in_one_round() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tC
S\t7\tGG
S\t8\tTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t1\t+\t7\t+\t0M
L\t7\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t6\t+\t0M
L\t4\t+\t8\t+\t0M
L\t8\t+\t6\t+\t0M
P\tref\t1+,2+,3+,4+,5+,6+\t*
P\tleft_alt\t1+,7+,3+,4+,5+,6+\t*
P\tright_alt\t1+,2+,3+,4+,8+,6+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                max_iterations: 1,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.iterations, 1);
        assert_eq!(resolved.stats.resolved, 2);
        assert_eq!(resolved.stats.bailed, 0);
    }

    /// Phase 6 hard invariant: every replacement segment id is strictly
    /// greater than the largest segment id present in the input graph at
    /// the moment of insertion.
    ///
    /// See `docs/crush-architecture-spec.md` §"fresh node IDs" and the
    /// task description's hard validation gate ("Fresh-ID assertion").
    #[test]
    fn phase6_fresh_segment_ids_are_strictly_above_initial_max() {
        // Two independent bubbles + a non-integer survivor — the original-max
        // integer id is 42, so every freshly minted segment id must be > 42.
        let gfa = "\
H\tVN:Z:1.0
S\t10\tA
S\t11\tC
S\t12\tG
S\t13\tT
S\t14\tA
S\t15\tC
S\t16\tGG
S\t17\tTT
S\t42\tACGT
L\t10\t+\t11\t+\t0M
L\t11\t+\t12\t+\t0M
L\t10\t+\t16\t+\t0M
L\t16\t+\t12\t+\t0M
L\t12\t+\t13\t+\t0M
L\t13\t+\t14\t+\t0M
L\t14\t+\t15\t+\t0M
L\t13\t+\t17\t+\t0M
L\t17\t+\t15\t+\t0M
L\t15\t+\t42\t+\t0M
P\tref\t10+,11+,12+,13+,14+,15+,42+\t*
P\tleft_alt\t10+,16+,12+,13+,14+,15+,42+\t*
P\tright_alt\t10+,11+,12+,13+,17+,15+,42+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                max_iterations: 3,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert!(resolved.stats.resolved >= 2, "{:?}", resolved.stats);

        let resolved_graph = parse_gfa(&resolved.gfa).unwrap();
        // Originals that survived keep their string id; any segment that did
        // not appear in the input is a freshly minted replacement.
        let original_ids = parse_gfa(gfa)
            .unwrap()
            .segments
            .iter()
            .map(|s| s.id.clone())
            .collect::<FxHashSet<_>>();
        let original_max_int = original_ids
            .iter()
            .filter_map(|id| id.parse::<usize>().ok())
            .max()
            .unwrap_or(0);
        let mut new_ids = Vec::new();
        for segment in &resolved_graph.segments {
            if !original_ids.contains(&segment.id) {
                let parsed: usize = segment
                    .id
                    .parse()
                    .expect("freshly minted id must be integer-encoded");
                new_ids.push(parsed);
            }
        }
        assert!(
            !new_ids.is_empty(),
            "expected at least one freshly minted segment id"
        );
        for id in &new_ids {
            assert!(
                *id > original_max_int,
                "fresh id {id} must be > original-max integer id {original_max_int}"
            );
        }
    }

    /// Phase 6 hard invariant: per-round frontier sizes are observed and
    /// non-empty, the run terminates by exhausting the tree (not by hitting
    /// max_iterations on a large iteration budget), and the resolved-set
    /// stays monotonic (no node resolved twice — that would panic in
    /// `resolve_graph_bubbles`).
    #[test]
    fn phase6_level_descent_terminates_below_max_iterations() {
        // Nested SNP-in-an-indel: round 1 should resolve the inner leaf, the
        // outer site then becomes a leaf and is resolved in round 2.
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tAA
S\t6\tA
S\t7\tCCC
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t2\t+\t5\t+\t0M
L\t5\t+\t4\t+\t0M
L\t4\t+\t6\t+\t0M
L\t1\t+\t7\t+\t0M
L\t7\t+\t6\t+\t0M
P\tref\t1+,2+,3+,4+,6+\t*
P\tinner_alt\t1+,2+,5+,4+,6+\t*
P\touter_alt\t1+,7+,6+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                max_iterations: 32,
                method: ResolutionMethod::Poa,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert!(
            resolved.stats.iterations < 32,
            "level descent must terminate before exhausting max-iterations; got {:?}",
            resolved.stats
        );
    }

    #[test]
    fn all_candidates_processed_regardless_of_budget() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tGGGGGGGG
S\t7\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t7\t+\t0M
L\t4\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,4+,5+,7+\t*
P\tleft_alt\t1+,3+,4+,5+,7+\t*
P\tright_alt\t1+,2+,4+,6+,7+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                max_iterations: 2,
                method: ResolutionMethod::Poa,
                max_traversal_len: 5,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert!(resolved.stats.resolved >= 2, "{:?}", resolved.stats);
        assert_eq!(resolved.stats.bailed, 0, "{:?}", resolved.stats);
    }

    #[test]
    fn chain_greedy_groups_adjacent_reference_bubbles() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tGGGG
S\t7\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t7\t+\t0M
L\t4\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,4+,5+,7+\t*
P\tleft_alt\t1+,3+,4+,5+,7+\t*
P\tright_alt\t1+,2+,4+,6+,7+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let config = ResolutionConfig {
            max_iterations: 1,
            method: ResolutionMethod::ChainGreedy,
            chain_greedy_target_bp: 100,
            ..ResolutionConfig::default()
        };
        let frontier =
            find_path_walk_chain_frontier(&graph, &config, &FxHashSet::default(), false).unwrap();
        assert_eq!(frontier.selected.len(), 1, "{:?}", frontier);
        assert_eq!(frontier.selected[0].source_bubbles, 2);

        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(gfa, &config).unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1, "{:?}", resolved.stats);
        assert_eq!(resolved.stats.bailed, 0, "{:?}", resolved.stats);
    }

    #[test]
    fn iterative_multi_level_generates_sibling_and_sliding_windows() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tGGGG
S\t7\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t7\t+\t0M
L\t4\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,4+,5+,7+\t*
P\tleft_alt\t1+,3+,4+,5+,7+\t*
P\tright_alt\t1+,2+,4+,6+,7+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let config = ResolutionConfig {
            method: ResolutionMethod::IterativeMultiLevel,
            multi_level_window_mode: MultiLevelWindowMode::Sibling,
            multi_level_window_target_bp: 100,
            multi_level_max_window_sites: 3,
            multi_level_candidate_limit: 32,
            ..ResolutionConfig::default()
        };
        let (discovered, _, _) = discover_all_candidates(&graph, &config, false).unwrap();
        let generated =
            generate_multi_level_candidates(&graph, &config, &discovered, &FxHashSet::default());
        assert!(
            generated
                .candidates
                .iter()
                .any(|c| c.source == MultiLevelCandidateSource::SiblingRun && c.source_sites >= 2),
            "{:?}",
            generated.source_counts
        );
        let sliding_config = ResolutionConfig {
            multi_level_window_mode: MultiLevelWindowMode::Sliding,
            ..config
        };
        let generated = generate_multi_level_candidates(
            &graph,
            &sliding_config,
            &discovered,
            &FxHashSet::default(),
        );
        assert!(
            generated.candidates.iter().any(|c| c.source
                == MultiLevelCandidateSource::SlidingWindow
                && c.source_sites >= 2),
            "{:?}",
            generated.source_counts
        );
    }

    #[test]
    fn iterative_multi_level_parent_mode_uses_parent_frontier_sources() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tGGGG
S\t7\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t7\t+\t0M
L\t4\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,4+,5+,7+\t*
P\tleft_alt\t1+,3+,4+,5+,7+\t*
P\tright_alt\t1+,2+,4+,6+,7+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let config = ResolutionConfig {
            method: ResolutionMethod::IterativeMultiLevel,
            multi_level_window_mode: MultiLevelWindowMode::Parent,
            multi_level_window_target_bp: 100,
            multi_level_max_window_sites: 3,
            multi_level_candidate_limit: 32,
            ..ResolutionConfig::default()
        };
        let (discovered, _, _) = discover_all_candidates(&graph, &config, false).unwrap();
        let generated =
            generate_multi_level_candidates(&graph, &config, &discovered, &FxHashSet::default());
        assert!(!generated.candidates.is_empty(), "{generated:?}");
        assert!(
            generated.candidates.iter().all(|candidate| matches!(
                candidate.source,
                MultiLevelCandidateSource::TopLevel | MultiLevelCandidateSource::ParentDescendants
            )),
            "{:?}",
            generated.source_counts
        );
    }

    #[test]
    fn iterative_multi_level_generates_outward_residual_windows() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tC
S\t4\tC
S\t5\tC
S\t6\tC
S\t7\tC
S\t8\tC
S\t9\tC
S\t10\tC
S\t11\tC
S\t12\tCCCCCCCCCC
S\t13\tT
S\t14\tG
S\t15\tG
S\t16\tG
S\t17\tG
S\t18\tG
S\t19\tG
S\t20\tG
S\t21\tG
S\t22\tG
S\t23\tG
S\t24\tGGGGGGGGGG
S\t25\tA
S\t26\tT
S\t27\tT
S\t28\tT
S\t29\tT
S\t30\tT
S\t31\tT
S\t32\tT
S\t33\tT
S\t34\tT
S\t35\tT
S\t36\tTTTTTTTTTT
S\t37\tC
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
L\t7\t+\t8\t+\t0M
L\t8\t+\t9\t+\t0M
L\t9\t+\t10\t+\t0M
L\t10\t+\t11\t+\t0M
L\t11\t+\t13\t+\t0M
L\t1\t+\t12\t+\t0M
L\t12\t+\t13\t+\t0M
L\t13\t+\t14\t+\t0M
L\t14\t+\t15\t+\t0M
L\t15\t+\t16\t+\t0M
L\t16\t+\t17\t+\t0M
L\t17\t+\t18\t+\t0M
L\t18\t+\t19\t+\t0M
L\t19\t+\t20\t+\t0M
L\t20\t+\t21\t+\t0M
L\t21\t+\t22\t+\t0M
L\t22\t+\t23\t+\t0M
L\t23\t+\t25\t+\t0M
L\t13\t+\t24\t+\t0M
L\t24\t+\t25\t+\t0M
L\t25\t+\t26\t+\t0M
L\t26\t+\t27\t+\t0M
L\t27\t+\t28\t+\t0M
L\t28\t+\t29\t+\t0M
L\t29\t+\t30\t+\t0M
L\t30\t+\t31\t+\t0M
L\t31\t+\t32\t+\t0M
L\t32\t+\t33\t+\t0M
L\t33\t+\t34\t+\t0M
L\t34\t+\t35\t+\t0M
L\t35\t+\t37\t+\t0M
L\t25\t+\t36\t+\t0M
L\t36\t+\t37\t+\t0M
P\tref\t1+,2+,3+,4+,5+,6+,7+,8+,9+,10+,11+,13+,14+,15+,16+,17+,18+,19+,20+,21+,22+,23+,25+,26+,27+,28+,29+,30+,31+,32+,33+,34+,35+,37+\t*
P\talt\t1+,12+,13+,24+,25+,36+,37+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let config = ResolutionConfig {
            method: ResolutionMethod::IterativeMultiLevel,
            multi_level_window_mode: MultiLevelWindowMode::Outward,
            multi_level_window_target_bp: 100,
            multi_level_max_window_sites: 3,
            multi_level_candidate_limit: 32,
            ..ResolutionConfig::default()
        };
        let (discovered, _, _) = discover_all_candidates(&graph, &config, false).unwrap();
        let generated =
            generate_multi_level_candidates(&graph, &config, &discovered, &FxHashSet::default());

        assert!(
            generated.candidates.iter().any(|candidate| {
                candidate.source == MultiLevelCandidateSource::OutwardResidualWindow
                    && candidate.source_sites >= 3
            }),
            "{:?}",
            generated.source_counts
        );
        assert!(
            generated
                .candidates
                .iter()
                .all(|candidate| candidate.source
                    == MultiLevelCandidateSource::OutwardResidualWindow),
            "{:?}",
            generated.source_counts
        );

        let coverage_config = ResolutionConfig {
            method: ResolutionMethod::CoverageMultiBubble,
            multi_level_window_mode: MultiLevelWindowMode::Combined,
            multi_level_window_target_bp: 100,
            multi_level_max_window_sites: 3,
            multi_level_candidate_limit: 32,
            ..ResolutionConfig::default()
        };
        let generated = generate_multi_level_candidates(
            &graph,
            &coverage_config,
            &discovered,
            &FxHashSet::default(),
        );
        assert!(
            generated
                .candidates
                .iter()
                .any(|candidate| candidate.source
                    == MultiLevelCandidateSource::OutwardResidualWindow),
            "coverage-multi-bubble should include outward residual windows by default; sources={:?}, repeat_boundary_rejected={}",
            generated.source_counts,
            generated.repeat_boundary_rejected
        );

        let before = seq_map(gfa);
        let dry_run = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                method: ResolutionMethod::IterativeMultiLevel,
                multi_level_window_mode: MultiLevelWindowMode::Outward,
                multi_level_window_target_bp: 100,
                multi_level_max_window_sites: 3,
                multi_level_candidate_limit: 32,
                multi_level_admission_only: true,
                auto_spoa_max_traversal_len: 0,
                auto_poasta_max_traversal_len: 1_000_000,
                // Poison-pill POASTA scoring: replacement building would reject
                // this tuple, so success proves admission-only exits first.
                scoring_params: (1, 4, 6, 1, 26, 2),
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&dry_run.gfa));
        assert_eq!(dry_run.stats.resolved, 0, "{:?}", dry_run.stats);
        assert_eq!(dry_run.stats.bailed, 0, "{:?}", dry_run.stats);
        assert!(
            dry_run.stats.candidates_seen > 0,
            "admission-only mode should still report admitted candidates: {:?}",
            dry_run.stats
        );
    }

    #[test]
    fn outward_residual_build_budget_diagnostic_log_includes_source_and_limit() {
        let graph = tiny_boundary_graph();
        let window = budget_test_window(MultiLevelCandidateSource::OutwardResidualWindow);
        let config = ResolutionConfig {
            method: ResolutionMethod::IterativeMultiLevel,
            max_pair_alignments: 10,
            ..ResolutionConfig::default()
        };
        let diagnostic =
            multi_level_build_budget_diagnostic(&window, &config, ResolutionMethod::Sweepga)
                .expect("outward residual window should exceed pair budget");

        assert_eq!(diagnostic.kind, MultiLevelBuildBudgetKind::PairAlignments);
        let line = format_multi_level_budget_diagnostic(&graph, &diagnostic);
        assert!(line.contains("source=outward-residual-window"), "{line}");
        assert!(line.contains("method=Sweepga"), "{line}");
        assert!(line.contains("budget-exceeded=pair-alignments"), "{line}");
        assert!(line.contains("diagnostic only"), "{line}");
        assert!(line.contains("observed=380"), "{line}");
        assert!(line.contains("limit=10"), "{line}");
    }

    #[test]
    fn outward_residual_build_budget_diagnostic_reports_poasta_work() {
        let window = budget_test_window(MultiLevelCandidateSource::OutwardResidualWindow);
        let config = ResolutionConfig {
            method: ResolutionMethod::IterativeMultiLevel,
            multi_level_max_poasta_cells: 1_000,
            ..ResolutionConfig::default()
        };
        let diagnostic =
            multi_level_build_budget_diagnostic(&window, &config, ResolutionMethod::Poasta)
                .expect("outward residual window should exceed POASTA budget");

        assert_eq!(diagnostic.kind, MultiLevelBuildBudgetKind::PoastaCells);
        assert_eq!(diagnostic.observed, 4_000_000);
        assert_eq!(diagnostic.limit, 1_000);
    }

    #[test]
    fn iterative_multi_level_routes_multi_bubble_windows_through_direct_auto_tiers() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tGGGG
S\t7\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t7\t+\t0M
L\t4\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,4+,5+,7+\t*
P\tleft_alt\t1+,3+,4+,5+,7+\t*
P\tright_alt\t1+,2+,4+,6+,7+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let config = ResolutionConfig {
            method: ResolutionMethod::IterativeMultiLevel,
            auto_spoa_max_traversal_len: 100,
            auto_poasta_max_traversal_len: 1_000,
            multi_level_window_mode: MultiLevelWindowMode::Sibling,
            multi_level_window_target_bp: 100,
            multi_level_max_window_sites: 3,
            multi_level_candidate_limit: 32,
            ..ResolutionConfig::default()
        };
        let (discovered, _, _) = discover_all_candidates(&graph, &config, false).unwrap();
        let generated =
            generate_multi_level_candidates(&graph, &config, &discovered, &FxHashSet::default());
        let window = generated
            .candidates
            .iter()
            .find(|candidate| {
                candidate.source == MultiLevelCandidateSource::SiblingRun
                    && candidate.source_sites >= 2
            })
            .expect("expected at least one multi-site sibling window");

        assert_eq!(
            auto_method_by_median(window.candidate.traversal_stats, &config),
            ResolutionMethod::Poa,
            "without the multi-window override this toy window would route to tiny-bubble POA"
        );
        assert_eq!(
            multi_level_window_replacement_method(window, &config),
            ResolutionMethod::Poasta
        );

        let forced_sweepga = ResolutionConfig {
            auto_spoa_max_traversal_len: 0,
            auto_poasta_max_traversal_len: 0,
            ..config.clone()
        };
        assert_eq!(
            multi_level_window_replacement_method(window, &forced_sweepga),
            ResolutionMethod::Sweepga
        );
    }

    #[test]
    fn iterative_multi_level_reports_total_bp_cap_without_filtering() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tCCCC
S\t3\tGGGG
S\t4\tT
S\t5\tA
S\t6\tGGGG
S\t7\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t7\t+\t0M
L\t4\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,4+,5+,7+\t*
P\tleft_alt\t1+,3+,4+,5+,7+\t*
P\tright_alt\t1+,2+,4+,6+,7+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let config = ResolutionConfig {
            method: ResolutionMethod::IterativeMultiLevel,
            multi_level_window_mode: MultiLevelWindowMode::Combined,
            multi_level_window_target_bp: 100,
            multi_level_max_window_sites: 3,
            multi_level_candidate_limit: 32,
            max_total_sequence: 12,
            ..ResolutionConfig::default()
        };
        let (discovered, _, _) = discover_all_candidates(&graph, &config, false).unwrap();
        let generated =
            generate_multi_level_candidates(&graph, &config, &discovered, &FxHashSet::default());

        assert!(generated.total_bp_cap_exceeded > 0, "{generated:?}");
        assert!(
            generated.candidates.iter().any(|candidate| candidate
                .candidate
                .traversal_stats
                .total_len
                > 12),
            "{:?}",
            generated
                .candidates
                .iter()
                .map(|candidate| candidate.candidate.traversal_stats.total_len)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn iterative_multi_level_reports_traversal_len_caps_without_filtering() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tCCCC
S\t3\tGGGG
S\t4\tT
S\t5\tA
S\t6\tGGGG
S\t7\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t7\t+\t0M
L\t4\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,4+,5+,7+\t*
P\tleft_alt\t1+,3+,4+,5+,7+\t*
P\tright_alt\t1+,2+,4+,6+,7+\t*
";
        let graph = parse_gfa(gfa).unwrap();
        let base_config = ResolutionConfig {
            method: ResolutionMethod::IterativeMultiLevel,
            multi_level_window_mode: MultiLevelWindowMode::Combined,
            multi_level_window_target_bp: 100,
            multi_level_max_window_sites: 3,
            multi_level_candidate_limit: 32,
            max_total_sequence: 1_000,
            ..ResolutionConfig::default()
        };
        let (discovered, _, _) = discover_all_candidates(&graph, &base_config, false).unwrap();

        let max_len_config = ResolutionConfig {
            max_traversal_len: 3,
            max_median_traversal_len: 1_000,
            ..base_config.clone()
        };
        let generated = generate_multi_level_candidates(
            &graph,
            &max_len_config,
            &discovered,
            &FxHashSet::default(),
        );
        assert!(generated.max_len_cap_exceeded > 0, "{generated:?}");
        assert!(
            generated.candidates.iter().any(|candidate| candidate
                .candidate
                .traversal_stats
                .max_len
                > 3),
            "{generated:?}"
        );

        let median_len_config = ResolutionConfig {
            max_traversal_len: 1_000,
            max_median_traversal_len: 3,
            ..base_config
        };
        let generated = generate_multi_level_candidates(
            &graph,
            &median_len_config,
            &discovered,
            &FxHashSet::default(),
        );
        assert!(generated.median_len_cap_exceeded > 0, "{generated:?}");
        assert!(
            generated.candidates.iter().any(|candidate| candidate
                .candidate
                .traversal_stats
                .median_len
                > 3),
            "{generated:?}"
        );
    }

    #[test]
    fn iterative_multi_level_caps_replacement_build_parallelism() {
        assert_eq!(multi_level_replacement_build_parallelism(0, 192), 1);
        assert_eq!(multi_level_replacement_build_parallelism(1, 192), 1);
        assert_eq!(multi_level_replacement_build_parallelism(2, 192), 2);
        assert_eq!(multi_level_replacement_build_parallelism(32, 192), 4);
        assert_eq!(multi_level_replacement_build_parallelism(32, 3), 3);
    }

    #[test]
    fn iterative_multi_level_preserves_paths_when_no_size_winner_exists() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tref\t1+,2+,4+\t*
P\talt\t1+,3+,4+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                max_iterations: 2,
                method: ResolutionMethod::IterativeMultiLevel,
                multi_level_window_mode: MultiLevelWindowMode::Combined,
                multi_level_window_target_bp: 100,
                multi_level_candidate_limit: 8,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
    }

    #[test]
    fn iterative_multi_level_materializes_window_sequences_before_build() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tGGGG
S\t7\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t7\t+\t0M
L\t4\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,4+,5+,7+\t*
P\tleft_alt\t1+,3+,4+,5+,7+\t*
P\tright_alt\t1+,2+,4+,6+,7+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                max_iterations: 1,
                method: ResolutionMethod::IterativeMultiLevel,
                auto_spoa_max_traversal_len: 100,
                auto_poasta_max_traversal_len: 0,
                multi_level_window_mode: MultiLevelWindowMode::Sibling,
                multi_level_window_target_bp: 100,
                multi_level_candidate_limit: 8,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert!(resolved.stats.candidates_seen > 0, "{:?}", resolved.stats);
        assert!(
            resolved.stats.bailed < resolved.stats.candidates_seen,
            "{:?}",
            resolved.stats
        );
    }

    #[test]
    fn direct_poa_processes_candidate_regardless_of_max_traversal_len() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAC
S\t2\tGGG
S\t3\tTA
L\t1\t+\t3\t+\t0M
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tref\t1+,3+\t*
P\tins\t1+,2+,3+\t*
";
        let before = seq_map(gfa);
        let config = ResolutionConfig {
            method: ResolutionMethod::Poa,
            max_traversal_len: 4,
            ..ResolutionConfig::default()
        };
        let resolved = resolve_gfa_bubbles(gfa, &config).unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1, "{:?}", resolved.stats);
        assert_eq!(resolved.stats.bailed, 0, "{:?}", resolved.stats);
    }

    #[test]
    fn skips_ambiguous_repeated_boundaries() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
L\t1\t+\t2\t+\t0M
L\t2\t+\t1\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tref\t1+,2+,1+,4+\t*
P\talt\t1+,3+,4+\t*
";
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(gfa, &ResolutionConfig::default()).unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.gfa, gfa);
        assert_eq!(resolved.stats.resolved, 0);
    }

    #[test]
    fn indexed_anchor_lookup_rejects_ambiguous_entry_exit_pairs() {
        let entry = Step {
            node: 0,
            rev: false,
        };
        let exit = Step {
            node: 2,
            rev: false,
        };
        let ambiguous = Path {
            name: "ambiguous".to_string(),
            steps: vec![
                entry,
                Step {
                    node: 1,
                    rev: false,
                },
                exit,
                entry,
                Step {
                    node: 1,
                    rev: false,
                },
                exit,
            ],
        };
        assert_eq!(
            unique_anchor_range(&path_step_index(&ambiguous), entry, exit),
            None
        );

        let unique = Path {
            name: "unique".to_string(),
            steps: vec![
                entry,
                Step {
                    node: 1,
                    rev: false,
                },
                exit,
                entry,
            ],
        };
        assert_eq!(
            unique_anchor_range(&path_step_index(&unique), entry, exit),
            Some((0, 3))
        );
    }

    #[test]
    fn occupied_interval_frontier_only_rejects_overlaps_on_same_path() {
        let mut occupied_by_path = vec![Vec::new(), Vec::new()];
        let selected = BubbleCandidate {
            ranges: vec![
                PathRange {
                    path_idx: 0,
                    begin_step: 10,
                    end_step: 20,
                    ..PathRange::default()
                },
                PathRange {
                    path_idx: 1,
                    begin_step: 30,
                    end_step: 40,
                    ..PathRange::default()
                },
            ],
            signature: "selected".to_string(),
            root_start_step: 0,
            root_end_step: 1,
            root_span: 10,
            total_steps: 20,
            unique_steps: 2,
            traversal_stats: TraversalStats::default(),
            level: 0,
        };
        mark_candidate_occupied(&mut occupied_by_path, &selected);

        let adjacent = BubbleCandidate {
            ranges: vec![PathRange {
                path_idx: 0,
                begin_step: 20,
                end_step: 25,
                ..PathRange::default()
            }],
            signature: "adjacent".to_string(),
            root_start_step: 0,
            root_end_step: 1,
            root_span: 5,
            total_steps: 5,
            unique_steps: 1,
            traversal_stats: TraversalStats::default(),
            level: 0,
        };
        assert!(!candidate_conflicts_with_occupied(
            &occupied_by_path,
            &adjacent
        ));

        let overlapping = BubbleCandidate {
            ranges: vec![PathRange {
                path_idx: 0,
                begin_step: 19,
                end_step: 25,
                ..PathRange::default()
            }],
            signature: "overlapping".to_string(),
            root_start_step: 0,
            root_end_step: 1,
            root_span: 6,
            total_steps: 6,
            unique_steps: 1,
            traversal_stats: TraversalStats::default(),
            level: 0,
        };
        assert!(candidate_conflicts_with_occupied(
            &occupied_by_path,
            &overlapping
        ));

        let same_coordinates_different_path = BubbleCandidate {
            ranges: vec![PathRange {
                path_idx: 1,
                begin_step: 10,
                end_step: 20,
                ..PathRange::default()
            }],
            signature: "different-path".to_string(),
            root_start_step: 0,
            root_end_step: 1,
            root_span: 10,
            total_steps: 10,
            unique_steps: 1,
            traversal_stats: TraversalStats::default(),
            level: 0,
        };
        assert!(!candidate_conflicts_with_occupied(
            &occupied_by_path,
            &same_coordinates_different_path
        ));
    }

    #[test]
    fn rejects_non_blunt_links() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
S\t2\tGTAA
L\t1\t+\t2\t+\t2M
P\tp\t1+,2+\t*
";
        let err = resolve_gfa_bubbles(gfa, &ResolutionConfig::default()).unwrap_err();
        assert!(err.to_string().contains("blunt 0M"));
    }

    /// Compression ratio = sum(input bp) / sum(output bp). Guard cases where
    /// either side is zero (empty replacement, zero input) must return None so
    /// diagnostic logging treats them as "no signal" rather than crashing on
    /// a divide-by-zero or grading every empty graph as infinitely well
    /// compressed.
    #[test]
    fn replacement_compression_ratio_handles_zero_and_basic_cases() {
        let candidate = BubbleCandidate {
            ranges: Vec::new(),
            signature: "test".to_string(),
            root_start_step: 0,
            root_end_step: 0,
            root_span: 0,
            total_steps: 0,
            unique_steps: 0,
            traversal_stats: TraversalStats {
                count: 4,
                min_len: 100,
                median_len: 100,
                p90_len: 100,
                max_len: 100,
                total_len: 400,
            },
            level: 0,
        };
        let mk_replacement = |segs: &[(&str, &[u8])]| Graph {
            segments: segs
                .iter()
                .map(|(id, seq)| Segment {
                    id: id.to_string(),
                    seq: seq.to_vec(),
                })
                .collect(),
            paths: Vec::new(),
        };

        // 4 traversals @ 100 bp collapsed to a 100 bp single shared segment
        // is the canonical "perfect compression" case — ratio = 4.0.
        let replacement = mk_replacement(&[("1", &[b'A'; 100])]);
        let ratio = replacement_compression_ratio(&candidate, &replacement).unwrap();
        assert!((ratio - 4.0).abs() < 1e-9, "got {}", ratio);

        // Same input, but the replacement kept every traversal as its own
        // segment — no compression. Ratio = 1.0.
        let replacement = mk_replacement(&[
            ("1", &[b'A'; 100]),
            ("2", &[b'C'; 100]),
            ("3", &[b'G'; 100]),
            ("4", &[b'T'; 100]),
        ]);
        let ratio = replacement_compression_ratio(&candidate, &replacement).unwrap();
        assert!((ratio - 1.0).abs() < 1e-9, "got {}", ratio);

        // Empty replacement → cannot divide; must return None so retry
        // ignores it instead of treating empty as well-compressed.
        let replacement = mk_replacement(&[]);
        assert!(replacement_compression_ratio(&candidate, &replacement).is_none());

        // Zero input bp → no signal either.
        let mut zero_input = candidate.clone();
        zero_input.traversal_stats.total_len = 0;
        let replacement = mk_replacement(&[("1", &[b'A'; 10])]);
        assert!(replacement_compression_ratio(&zero_input, &replacement).is_none());
    }

    /// End-to-end: compression-ratio thresholds are diagnostic only. Even an
    /// aggressive threshold must not rebuild with an alternate aligner, swap a
    /// replacement, or reject an otherwise path-valid replacement.
    #[test]
    fn compression_ratio_threshold_is_diagnostic_only() {
        // Reuse the canonical snp bubble fixture from
        // `resolves_simple_snp_bubble_without_changing_path_sequences`.
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tC
S\t7\tG
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t1\t+\t4\t+\t0M
L\t4\t+\t3\t+\t0M
L\t3\t+\t5\t+\t0M
L\t5\t+\t6\t+\t0M
L\t6\t+\t7\t+\t0M
P\tref\t1+,2+,3+,5+,6+,7+\t*
P\talt\t1+,4+,3+,5+,6+,7+\t*
";
        let before = seq_map(gfa);
        assert_eq!(ResolutionConfig::default().retry_min_compression_ratio, 0.0);

        let cfg_aggressive = ResolutionConfig {
            retry_min_compression_ratio: 1e9,
            retry_min_input_bp: 0,
            ..ResolutionConfig::default()
        };
        let aggressive = resolve_gfa_bubbles(gfa, &cfg_aggressive).unwrap();
        assert_eq!(
            before,
            seq_map(&aggressive.gfa),
            "diagnostic compression thresholds must preserve every input path sequence"
        );
        assert_eq!(
            aggressive.stats.retry_attempts, 0,
            "compression threshold must not trigger alternate aligner attempts"
        );
        assert_eq!(
            aggressive.stats.retry_wins, 0,
            "compression threshold must not swap replacements by metric"
        );
        assert_eq!(
            aggressive.stats.retry_failures, 0,
            "compression threshold must not record retry failures"
        );
    }
}
