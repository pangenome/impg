//! Bubble-guided graph resolution primitives.
//!
//! This module implements the first conservative slice of hierarchical graph
//! resolution: detect path-supported bubbles in a blunt GFA, replace bounded
//! single-entry/single-exit bubbles with exact path-preserving local graph
//! induction, and repeat until no eligible unseen bubbles remain. The default
//! resolver uses AllWave/seqwish followed by a small SPOA polish. It
//! intentionally avoids lossy representative collapse and coordinate sidecars;
//! emitted paths are the coordinate system.

use crate::graph::{
    build_global_spoa_engine, feed_sequences_to_graph, reverse_complement, unchop_gfa,
};
use povu::native_gfa::{Step as PovuStep, Strand as PovuStrand};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::{BTreeSet, HashSet};
use std::hash::{Hash, Hasher};
use std::io::{self, Read};
use std::time::Instant;

/// Configuration for exact path-preserving bubble resolution.
#[derive(Clone, Debug)]
pub struct ResolutionConfig {
    /// Maximum number of frontier replacement rounds.
    pub max_iterations: usize,
    /// Local graph induction method used for selected bubbles.
    pub method: ResolutionMethod,
    /// Maximum traversal length for `method=auto` to use direct SPOA.
    ///
    /// Larger eligible single-entry/single-exit bubbles go through AllWave /
    /// seqwish induction first, then the bounded SPOA polish pass. Set to 0 to
    /// force the AllWave-first path for every auto-selected bubble.
    pub auto_spoa_max_traversal_len: usize,
    /// Maximum total traversal sequence for `method=auto` to route to
    /// AllWave/seqwish. Dense many-haplotype local bubbles above this limit
    /// stay on direct SPOA, because seqwish transclosure can dominate runtime
    /// in repetitive loci. Set to 0 to disable this guard.
    pub auto_allwave_max_total_sequence: usize,
    /// Maximum number of traversals for `method=auto` to route to
    /// AllWave/seqwish. Dense high-copy bubbles above this limit stay on
    /// direct SPOA, which is more predictable for hundreds of short traversals.
    /// Set to 0 to disable this guard.
    pub auto_allwave_max_traversals: usize,
    /// Maximum root path span, in bp, for a candidate bubble.
    ///
    /// The root path is the first path in the input GFA. This is a rooted POVU
    /// decomposition coordinate, not necessarily an external reference genome.
    pub max_bubble_span: usize,
    /// Maximum length of any observed traversal through the bubble.
    pub max_traversal_len: usize,
    /// Minimum length of the longest observed traversal through the bubble.
    ///
    /// This lets callers run hierarchical passes explicitly, e.g. first crush
    /// parent-scale bubbles with a pairwise graph-induction engine, then run a
    /// separate SPOA pass over remaining small local bubbles.
    pub min_traversal_len: usize,
    /// Maximum median traversal length. A value of 0 disables this guard.
    pub max_median_traversal_len: usize,
    /// Maximum sum of traversal sequence lengths for one replacement.
    pub max_total_sequence: usize,
    /// Maximum number of path-supported traversals through one replacement.
    pub max_traversals: usize,
    /// One-pass small-tangle polish rounds after BiWFA in-memory induction.
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
    /// Minimum exact-match length used by seqwish when inducing replacement graphs.
    ///
    /// This is deliberately separate from the outer graph builder defaults:
    /// bubble-local induction in repetitive loci should ignore short
    /// off-diagonal matches while still allowing the pairwise alignment engine
    /// to align full traversal sequences end-to-end.
    pub replacement_seqwish_min_match_len: u64,
    /// SweepGA aligner backend used by `method=sweepga`.
    pub sweepga_aligner: String,
    /// SweepGA/FastGA k-mer frequency.
    pub sweepga_kmer_frequency: usize,
    /// SweepGA minimum alignment length.
    pub sweepga_min_aln_length: u64,
    /// SweepGA wfmash percent identity. `None` lets wfmash auto-estimate.
    pub sweepga_map_pct_identity: Option<String>,
    /// Maximum selected pair alignments for one pairwise graph induction.
    ///
    /// The selected-pair graph is only a guide for seqwish transitive closure.
    /// Past this point, repetitive parent bubbles can produce dense enough
    /// alignments that graph induction is dominated by closure rather than
    /// useful local condensation. Set to 0 to disable this guard.
    pub max_pair_alignments: usize,
    /// Maximum PAF bytes handed to seqwish for one replacement.
    ///
    /// This catches parent-scale candidates whose pair count looks acceptable
    /// but whose alignments cover too many repetitive positions. Set to 0 to
    /// disable this guard.
    pub max_replacement_paf_bytes: usize,
    /// Maximum allowed round-level graph-complexity score growth.
    ///
    /// A round whose score grows by more than this fraction is rolled back and
    /// stops the run. The score combines segment sequence, path steps, links,
    /// and path white-space bridges in current segment order.
    pub max_round_score_growth: f64,
    /// Required fractional graph-complexity score improvement per round.
    ///
    /// If positive, rounds with less improvement are rolled back and stop the
    /// run. This is useful for exploratory multi-round crushing where parent
    /// bubbles become expensive after the main condensation has happened.
    pub min_round_score_improvement: f64,
    /// Disable round-level quality acceptance checks.
    pub disable_round_quality_check: bool,
    /// Skip SweepGA post-alignment filtering and feed selected pair alignments directly.
    pub sweepga_no_filter: bool,
    /// Force sparse pair dispatch for SweepGA/FastGA graph induction.
    ///
    /// FastGA does not have a cheap pairs-file mode here, so the default for
    /// FastGA is one all-vs-all batch followed by SweepGA filtering. Sparse
    /// pairs are still useful with wfmash, where the pairs file is honored by a
    /// single aligner invocation.
    pub sweepga_sparse_pairs: bool,
    /// SPOA scoring parameters: (match, mismatch, gap_open, gap_ext, gap_open2, gap_ext2).
    pub scoring_params: (u8, u8, u8, u8, u8, u8),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ResolutionMethod {
    Auto,
    /// Direct SPOA replacement.
    Poa,
    /// POASTA replacement, falling back to SPOA if graph export fails validation.
    Poasta,
    /// Debug/diagnostic resolver: align every traversal to one root with BiWFA
    /// and build a star-column graph. This is path-preserving but intentionally
    /// not the graph-quality default.
    StarBiwfa,
    /// Sparse many-to-many BiWFA with AllWave, seqwish induction, then SPOA polish.
    Allwave,
    /// SweepGA/FastGA or wfmash pair selection and graph induction, then SPOA polish.
    Sweepga,
}

impl ResolutionMethod {
    pub fn parse_name(value: &str) -> Option<Self> {
        match value.replace('_', "-").to_ascii_lowercase().as_str() {
            "auto" => Some(Self::Auto),
            "poa" | "spoa" => Some(Self::Poa),
            "poasta" => Some(Self::Poasta),
            "star-biwfa" | "biwfa-star" | "biwfa" | "wfa" => Some(Self::StarBiwfa),
            "allwave" | "aw" => Some(Self::Allwave),
            "sweepga" | "sw" => Some(Self::Sweepga),
            _ => None,
        }
    }
}

/// One frontier round is the fast default. More rounds are explicit polishing.
pub const DEFAULT_MAX_ITERATIONS: usize = 1;
/// By default, do not cap by rooted path span.
///
/// `max_bubble_span` is a POVU root-path coordinate guard. The root is currently
/// the first GFA path, so this is intentionally not part of the default runtime
/// budget. Use traversal length/count/total sequence for the real work cap.
pub const DEFAULT_MAX_BUBBLE_SPAN: usize = 0;
pub const DEFAULT_MAX_TRAVERSAL_LEN: usize = 10_000;
pub const DEFAULT_MIN_TRAVERSAL_LEN: usize = 0;
pub const DEFAULT_MAX_MEDIAN_TRAVERSAL_LEN: usize = 1_000;
pub const DEFAULT_MAX_TOTAL_SEQUENCE: usize = 1_000_000;
pub const DEFAULT_MAX_TRAVERSALS: usize = 10_000;
pub const DEFAULT_AUTO_SPOA_MAX_TRAVERSAL_LEN: usize = 2_000;
pub const DEFAULT_AUTO_ALLWAVE_MAX_TOTAL_SEQUENCE: usize = 200_000;
pub const DEFAULT_AUTO_ALLWAVE_MAX_TRAVERSALS: usize = 128;
pub const DEFAULT_POLISH_ITERATIONS: usize = 1;
pub const DEFAULT_POLISH_MAX_TRAVERSAL_LEN: usize = 2_000;
pub const DEFAULT_POLISH_MAX_MEDIAN_TRAVERSAL_LEN: usize = 256;
pub const DEFAULT_POLISH_MAX_TOTAL_SEQUENCE: usize = 1_000_000;
pub const DEFAULT_POLISH_MAX_TRAVERSALS: usize = 10_000;
pub const DEFAULT_PAIR_K_NEAREST: usize = 3;
pub const DEFAULT_PAIR_K_FARTHEST: usize = 1;
pub const DEFAULT_PAIR_TREE_COUNT: usize = 1;
pub const DEFAULT_PAIR_RANDOM_FRACTION: f64 = 0.01;
pub const DEFAULT_PAIR_MASH_K: usize = 15;
pub const DEFAULT_REPLACEMENT_SEQWISH_MIN_MATCH_LEN: u64 = 79;
pub const DEFAULT_SWEEPGA_KMER_FREQUENCY: usize = 10;
pub const DEFAULT_SWEEPGA_MIN_ALN_LENGTH: u64 = 0;
pub const DEFAULT_MAX_PAIR_ALIGNMENTS: usize = 10_000;
pub const DEFAULT_MAX_REPLACEMENT_PAF_BYTES: usize = 64 * 1024 * 1024;
pub const DEFAULT_MAX_ROUND_SCORE_GROWTH: f64 = 0.02;
pub const DEFAULT_MIN_ROUND_SCORE_IMPROVEMENT: f64 = 0.0;

impl Default for ResolutionConfig {
    fn default() -> Self {
        Self {
            max_iterations: DEFAULT_MAX_ITERATIONS,
            method: ResolutionMethod::Auto,
            auto_spoa_max_traversal_len: DEFAULT_AUTO_SPOA_MAX_TRAVERSAL_LEN,
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
            pair_k_nearest: DEFAULT_PAIR_K_NEAREST,
            pair_k_farthest: DEFAULT_PAIR_K_FARTHEST,
            pair_tree_count: DEFAULT_PAIR_TREE_COUNT,
            pair_random_fraction: DEFAULT_PAIR_RANDOM_FRACTION,
            pair_mash_k: DEFAULT_PAIR_MASH_K,
            replacement_seqwish_min_match_len: DEFAULT_REPLACEMENT_SEQWISH_MIN_MATCH_LEN,
            sweepga_aligner: "fastga".to_string(),
            sweepga_kmer_frequency: DEFAULT_SWEEPGA_KMER_FREQUENCY,
            sweepga_min_aln_length: DEFAULT_SWEEPGA_MIN_ALN_LENGTH,
            sweepga_map_pct_identity: None,
            max_pair_alignments: DEFAULT_MAX_PAIR_ALIGNMENTS,
            max_replacement_paf_bytes: DEFAULT_MAX_REPLACEMENT_PAF_BYTES,
            max_round_score_growth: DEFAULT_MAX_ROUND_SCORE_GROWTH,
            min_round_score_improvement: DEFAULT_MIN_ROUND_SCORE_IMPROVEMENT,
            disable_round_quality_check: false,
            sweepga_no_filter: false,
            sweepga_sparse_pairs: false,
            scoring_params: (1, 4, 6, 2, 26, 1),
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

#[derive(Clone, Debug)]
struct PathRange {
    path_idx: usize,
    begin_step: usize,
    end_step: usize,
    sequence: Vec<u8>,
}

#[derive(Clone, Debug)]
struct BubbleCandidate {
    ranges: Vec<PathRange>,
    signature: String,
    within_budget: bool,
    root_start_step: usize,
    root_end_step: usize,
    root_span: usize,
    total_steps: usize,
    unique_steps: usize,
    traversal_stats: TraversalStats,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct GraphQuality {
    segments: usize,
    segment_bp: usize,
    links: usize,
    path_steps: usize,
    path_white_space_bp_total: u64,
    path_white_space_bp_p99: usize,
    path_white_space_bp_max: usize,
    score: u64,
}

#[derive(Clone, Debug, Default)]
struct CandidateFrontier {
    selected: Vec<BubbleCandidate>,
    bailed: Vec<BubbleCandidate>,
    sites_seen: usize,
    candidates_seen: usize,
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

fn resolve_gfa_bubbles_quiet(gfa: &str, config: &ResolutionConfig) -> io::Result<ResolvedGfa> {
    let graph = parse_gfa(gfa)?;
    resolve_graph_bubbles(graph, Some(gfa), config, false)
}

fn resolve_graph_bubbles(
    mut graph: Graph,
    original_gfa: Option<&str>,
    config: &ResolutionConfig,
    emit_logs: bool,
) -> io::Result<ResolvedGfa> {
    let mut stats = ResolutionStats::default();
    let mut seen: FxHashSet<String> = FxHashSet::default();
    let mut changed = false;

    for round in 0..config.max_iterations {
        let round_start = Instant::now();
        let before_quality = graph_quality(&graph);
        let discovery_start = Instant::now();
        let frontier = find_candidate_frontier(&graph, config, &seen, emit_logs)?;
        let discovery_elapsed = discovery_start.elapsed();
        if frontier.selected.is_empty() && frontier.bailed.is_empty() {
            if emit_logs {
                log::info!(
                    "crush round {}: no eligible candidates from {} POVU site(s) in {:.2?}",
                    round + 1,
                    frontier.sites_seen,
                    discovery_elapsed
                );
            }
            break;
        }
        stats.iterations = round + 1;

        if emit_logs {
            log::info!(
                "crush round {}: {} POVU site(s), {} unseen polymorphic candidate(s), {} selected, {} over budget in {:.2?}",
                round + 1,
                frontier.sites_seen,
                frontier.candidates_seen,
                frontier.selected.len(),
                frontier.bailed.len(),
                discovery_elapsed
            );
            log::info!(
                "crush round {} traversal stats: {}; {}",
                round + 1,
                format_candidate_length_summary("selected", &frontier.selected),
                format_candidate_length_summary("over-budget", &frontier.bailed)
            );
        }

        for _candidate in frontier.bailed {
            stats.candidates_seen += 1;
            stats.bailed += 1;
        }

        if frontier.selected.is_empty() {
            if emit_logs {
                log::info!(
                    "crush round {}: stopping because remaining candidate(s) are over budget",
                    round + 1
                );
            }
            break;
        }

        for candidate in &frontier.selected {
            seen.insert(candidate.signature.clone());
            stats.candidates_seen += 1;
        }

        let selected_count = frontier.selected.len();
        let build_start = Instant::now();
        let build_results = frontier
            .selected
            .into_par_iter()
            .map(|candidate| {
                build_replacement(&candidate, config).map(|replacement| ReplacementPlan {
                    candidate,
                    replacement,
                })
            })
            .collect::<Vec<_>>();
        let build_elapsed = build_start.elapsed();

        let mut plans = Vec::new();
        for result in build_results {
            match result {
                Ok(plan) if !plan.replacement.segments.is_empty() => plans.push(plan),
                Ok(_) => stats.bailed += 1,
                Err(err) => {
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
                    "crush round {}: 0/{} replacement(s) built in {:.2?}",
                    round + 1,
                    selected_count,
                    build_elapsed
                );
            }
            continue;
        }

        let resolved_count = plans.len();
        let rewrite_start = Instant::now();
        let next_graph = apply_replacement_frontier(&graph, &plans)?;
        let rewrite_elapsed = rewrite_start.elapsed();
        let after_quality = graph_quality(&next_graph);
        if !config.disable_round_quality_check {
            let decision = round_quality_decision(before_quality, after_quality, config);
            if let RoundQualityDecision::Reject { reason } = decision {
                stats.bailed += resolved_count;
                if emit_logs {
                    log::info!(
                        "crush round {}: rejecting {} replacement(s): {}; before {}; after {}",
                        round + 1,
                        resolved_count,
                        reason,
                        before_quality.summary(),
                        after_quality.summary()
                    );
                }
                break;
            }
        }
        graph = next_graph;
        changed = true;
        stats.resolved += resolved_count;
        if emit_logs {
            log::info!(
                "crush round {}: resolved {}/{} replacement(s) in {:.2?}; rewrite+validate {:.2?}; total {:.2?}; quality {} -> {}",
                round + 1,
                resolved_count,
                selected_count,
                build_elapsed,
                rewrite_elapsed,
                round_start.elapsed(),
                before_quality.summary(),
                after_quality.summary()
            );
        }
    }

    Ok(ResolvedGfa {
        gfa: if changed {
            render_graph(&graph)
        } else {
            original_gfa
                .map(str::to_string)
                .unwrap_or_else(|| render_graph(&graph))
        },
        stats,
    })
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

fn median_value(values: &mut [usize]) -> usize {
    if values.is_empty() {
        return 0;
    }
    values.sort_unstable();
    values[values.len() / 2]
}

enum RoundQualityDecision {
    Accept,
    Reject { reason: String },
}

impl GraphQuality {
    fn summary(&self) -> String {
        format!(
            "score={}, segments={}, segment-bp={}, links={}, path-steps={}, ws-total={}, ws-p99={}, ws-max={}",
            self.score,
            self.segments,
            self.segment_bp,
            self.links,
            self.path_steps,
            self.path_white_space_bp_total,
            self.path_white_space_bp_p99,
            self.path_white_space_bp_max
        )
    }
}

fn round_quality_decision(
    before: GraphQuality,
    after: GraphQuality,
    config: &ResolutionConfig,
) -> RoundQualityDecision {
    if before.score == 0 {
        return RoundQualityDecision::Accept;
    }

    let before_score = before.score as f64;
    let after_score = after.score as f64;
    let growth = (after_score - before_score) / before_score;
    if growth > config.max_round_score_growth {
        return RoundQualityDecision::Reject {
            reason: format!(
                "complexity score grew by {:.4}, above allowed {:.4}",
                growth, config.max_round_score_growth
            ),
        };
    }

    if config.min_round_score_improvement > 0.0 {
        let improvement = (before_score - after_score) / before_score;
        if improvement < config.min_round_score_improvement {
            return RoundQualityDecision::Reject {
                reason: format!(
                    "complexity score improved by {:.4}, below required {:.4}",
                    improvement, config.min_round_score_improvement
                ),
            };
        }
    }

    RoundQualityDecision::Accept
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

    // This intentionally approximates rendered graph cost while still pushing
    // down long path white-space bridges. The divisor keeps bridge cost large
    // enough to matter without overwhelming path-step/sequence condensation.
    let score = (segment_bp as u128)
        .saturating_add(path_steps as u128)
        .saturating_add((links.len() as u128).saturating_mul(4))
        .saturating_add(path_white_space_bp_total_u128 / 64)
        .min(u64::MAX as u128) as u64;

    GraphQuality {
        segments,
        segment_bp,
        links: links.len(),
        path_steps,
        path_white_space_bp_total,
        path_white_space_bp_p99,
        path_white_space_bp_max,
        score,
    }
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

fn find_candidate_frontier(
    graph: &Graph,
    config: &ResolutionConfig,
    seen: &FxHashSet<String>,
    emit_logs: bool,
) -> io::Result<CandidateFrontier> {
    if graph.paths.is_empty() || graph.paths[0].steps.len() < 2 {
        return Ok(CandidateFrontier::default());
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
    let mut candidates = decomposition
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
            if traversal_stats.max_len < config.min_traversal_len {
                return None;
            }
            let within_budget = (config.max_bubble_span == 0
                || root_span <= config.max_bubble_span)
                && ranges.len() <= config.max_traversals
                && traversal_stats.max_len <= config.max_traversal_len
                && (config.max_median_traversal_len == 0
                    || traversal_stats.median_len <= config.max_median_traversal_len)
                && traversal_stats.total_len <= config.max_total_sequence;
            let signature = candidate_signature(
                graph,
                root_path_idx,
                root_positions,
                begin,
                exit_step,
                &ranges,
            );
            if seen.contains(&signature) {
                return None;
            }

            Some(BubbleCandidate {
                ranges,
                signature,
                within_budget,
                root_start_step: begin,
                root_end_step: exit_step,
                root_span,
                total_steps,
                unique_steps: unique_steps.len(),
                traversal_stats,
            })
        })
        .collect::<Vec<_>>();
    let candidate_elapsed = candidate_start.elapsed();
    let candidates_seen = candidates.len();

    let select_start = Instant::now();
    candidates.sort_by(|a, b| {
        b.estimated_step_savings()
            .cmp(&a.estimated_step_savings())
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
        if !candidate.within_budget {
            frontier.bailed.push(candidate);
            continue;
        }
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
            materialize_candidate_sequences(graph, &mut candidate).then_some(candidate)
        })
        .collect();
    let materialize_elapsed = materialize_start.elapsed();
    let select_elapsed = select_start.elapsed();

    if emit_logs {
        log::info!(
            "crush discovery detail: render {:.2?}, povu-parse {:.2?}, povu-decompose {:.2?}, id-map {:.2?}, path-index {:.2?}, candidate-build {:.2?}, select+materialize {:.2?} ({} -> {} selected, materialize {:.2?})",
            render_elapsed,
            parse_elapsed,
            decompose_elapsed,
            id_map_elapsed,
            path_index_elapsed,
            candidate_elapsed,
            select_elapsed,
            selected_before_materialize,
            frontier.selected.len(),
            materialize_elapsed
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

fn materialize_candidate_sequences(graph: &Graph, candidate: &mut BubbleCandidate) -> bool {
    for range in &mut candidate.ranges {
        range.sequence = range_sequence(graph, range);
    }
    let Some(first) = candidate.ranges.first() else {
        return false;
    };
    candidate
        .ranges
        .iter()
        .any(|range| range.sequence != first.sequence)
}

fn candidate_signature(
    graph: &Graph,
    ref_path_idx: usize,
    ref_positions: &[usize],
    begin: usize,
    exit_step: usize,
    ranges: &[PathRange],
) -> String {
    let mut coords: Vec<String> = ranges
        .iter()
        .map(|r| format!("{}:{}-{}", r.path_idx, r.begin_step, r.end_step))
        .collect();
    coords.sort();
    format!(
        "{}:{}-{}:{}",
        graph.paths[ref_path_idx].name,
        ref_positions[begin],
        ref_positions[exit_step + 1],
        coords.join("|")
    )
}

fn apply_replacement_frontier(graph: &Graph, plans: &[ReplacementPlan]) -> io::Result<Graph> {
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
    let mut used_replacement: FxHashSet<(usize, usize)> = FxHashSet::default();

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
                        if let OutNode::Replacement(plan_idx, node_idx) = step.node {
                            used_replacement.insert((plan_idx, node_idx));
                        }
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
    let rendered = render_rewritten_graph(
        graph,
        &replacement_graphs,
        &used_original,
        &used_replacement,
        &out_paths,
    );
    let next = parse_gfa(&rendered)?;
    if !path_sequences_equal(graph, &next)? {
        return Err(io::Error::other(
            "resolved graph failed exact path-sequence validation",
        ));
    }
    Ok(next)
}

fn build_replacement(candidate: &BubbleCandidate, config: &ResolutionConfig) -> io::Result<Graph> {
    let method = match config.method {
        ResolutionMethod::Auto => auto_replacement_method(candidate, config),
        method => method,
    };
    match method {
        ResolutionMethod::Auto => unreachable!("auto is resolved before replacement dispatch"),
        ResolutionMethod::Poa => build_poa_replacement(candidate, config),
        ResolutionMethod::Poasta => build_poasta_replacement(candidate, config).or_else(|err| {
            log::debug!("POASTA replacement failed; falling back to SPOA: {err}");
            build_poa_replacement(candidate, config)
        }),
        ResolutionMethod::StarBiwfa => build_biwfa_inmemory_replacement(candidate, config),
        ResolutionMethod::Allwave => build_allwave_seqwish_replacement(candidate, config),
        ResolutionMethod::Sweepga => build_sweepga_seqwish_replacement(candidate, config),
    }
}

fn auto_replacement_method(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> ResolutionMethod {
    // Keep very small local tangles on the deterministic SPOA path. Above this
    // cutoff, use the intended many-to-many AllWave/seqwish path and then run
    // the bounded SPOA polish in `finalize_pairwise_induced_replacement`.
    if config.auto_spoa_max_traversal_len > 0
        && candidate.traversal_stats.max_len <= config.auto_spoa_max_traversal_len
    {
        ResolutionMethod::Poa
    } else if config.auto_allwave_max_total_sequence > 0
        && candidate.traversal_stats.total_len > config.auto_allwave_max_total_sequence
    {
        ResolutionMethod::Poa
    } else if config.auto_allwave_max_traversals > 0
        && candidate.traversal_stats.count > config.auto_allwave_max_traversals
    {
        ResolutionMethod::Poa
    } else {
        ResolutionMethod::Allwave
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
        .zip(candidate.ranges.iter().map(|range| range.sequence.clone()))
        .collect();
    (headers, seqs)
}

fn seqwish_replacement_config(
    config: &ResolutionConfig,
) -> crate::commands::graph::GraphBuildConfig {
    crate::commands::graph::GraphBuildConfig {
        num_threads: rayon::current_num_threads().max(1),
        show_progress: false,
        min_aln_length: 0,
        min_match_len: config.replacement_seqwish_min_match_len,
        input_paf: None,
        no_filter: config.sweepga_no_filter,
        num_mappings: "1:1".to_string(),
        scaffold_filter: "1:1".to_string(),
        sparsify: sweepga::knn_graph::SparsificationStrategy::None,
        ..crate::commands::graph::GraphBuildConfig::default()
    }
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
    let compacted = unchop_gfa(&gfa)?;
    let mut replacement = parse_gfa(&compacted)?;
    order_replacement_paths(&mut replacement, headers)?;
    if config.polish_iterations > 0 {
        let polished = polish_replacement_gfa(&render_graph(&replacement), config)?;
        replacement = parse_gfa(&polished)?;
        order_replacement_paths(&mut replacement, headers)?;
    }
    validate_replacement_paths(&replacement, candidate, method)?;
    Ok(replacement)
}

fn build_allwave_seqwish_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    let (headers, seqs) = candidate_named_sequences(candidate);
    let non_empty = seqs.iter().filter(|(_, seq)| !seq.is_empty()).count();
    if non_empty < 2 {
        log::debug!(
            "crush allwave: {} traversal(s), {} non-empty; falling back to SPOA",
            seqs.len(),
            non_empty
        );
        return build_poa_replacement(candidate, config);
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
        log::debug!(
            "crush allwave: selected zero pair alignments for {} traversal(s); falling back to SPOA",
            sequences.len()
        );
        return build_poa_replacement(candidate, config);
    }
    if config.max_pair_alignments > 0 && pairs.len() > config.max_pair_alignments {
        return Err(io::Error::other(format!(
            "AllWave replacement selected {} pair alignments, above max-pair-alignments {}",
            pairs.len(),
            config.max_pair_alignments
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
    if config.max_replacement_paf_bytes > 0 && paf.len() > config.max_replacement_paf_bytes {
        return Err(io::Error::other(format!(
            "AllWave replacement produced {} PAF bytes, above max-replacement-paf-bytes {}",
            paf.len(),
            config.max_replacement_paf_bytes
        )));
    }
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
    let (headers, seqs) = candidate_named_sequences(candidate);
    let named: Vec<(String, &[u8])> = seqs
        .iter()
        .map(|(name, seq)| (name.clone(), seq.as_slice()))
        .collect();
    let mut paf_lines = Vec::new();
    let sparse_pairs = config.sweepga_sparse_pairs || config.sweepga_aligner == "wfmash";
    let mash_schedule = if sparse_pairs {
        tree_mash_k_schedule(config.pair_mash_k, config.pair_tree_count)
    } else {
        vec![config.pair_mash_k]
    };
    for mash_k in mash_schedule {
        let align_config = sweepga::library_api::SweepgaAlignConfig {
            num_threads: rayon::current_num_threads().max(1),
            kmer_frequency: config.sweepga_kmer_frequency,
            min_aln_length: config.sweepga_min_aln_length,
            no_filter: true,
            sparsify: if sparse_pairs {
                sweepga::knn_graph::SparsificationStrategy::TreeSampling(
                    config.pair_k_nearest,
                    config.pair_k_farthest,
                    config.pair_random_fraction,
                )
            } else {
                sweepga::knn_graph::SparsificationStrategy::None
            },
            mash_params: sweepga::knn_graph::MashParams {
                kmer_size: mash_k,
                ..sweepga::knn_graph::MashParams::default()
            },
            aligner: config.sweepga_aligner.clone(),
            map_pct_identity: config.sweepga_map_pct_identity.clone(),
            ..sweepga::library_api::SweepgaAlignConfig::default()
        };
        let paf_file =
            sweepga::library_api::sweepga_align(&named, &align_config).map_err(|err| {
                io::Error::other(format!("SweepGA replacement alignment failed: {err}"))
            })?;
        let mut chunk = String::new();
        std::fs::File::open(paf_file.path())?.read_to_string(&mut chunk)?;
        paf_lines.extend(
            chunk
                .lines()
                .filter(|line| !line.is_empty())
                .map(str::to_owned),
        );
    }
    paf_lines.sort_unstable();
    paf_lines.dedup();
    let paf = if paf_lines.is_empty() {
        String::new()
    } else {
        let mut out = paf_lines.join("\n");
        out.push('\n');
        out
    };
    if config.max_pair_alignments > 0 && paf_lines.len() > config.max_pair_alignments {
        return Err(io::Error::other(format!(
            "SweepGA replacement produced {} pair alignments, above max-pair-alignments {}",
            paf_lines.len(),
            config.max_pair_alignments
        )));
    }
    if config.max_replacement_paf_bytes > 0 && paf.len() > config.max_replacement_paf_bytes {
        return Err(io::Error::other(format!(
            "SweepGA replacement produced {} PAF bytes, above max-replacement-paf-bytes {}",
            paf.len(),
            config.max_replacement_paf_bytes
        )));
    }
    log::debug!(
        "crush sweepga: {} traversal(s), {} unique PAF line(s), {} PAF byte(s) from {} backend ({})",
        seqs.len(),
        paf_lines.len(),
        paf.len(),
        config.sweepga_aligner,
        if sparse_pairs { "sparse pairs" } else { "all pairs" }
    );

    let graph_config = seqwish_replacement_config(config);
    let gfa = crate::syng_graph::build_gfa_from_paf_and_sequences(&seqs, &paf, &graph_config)?;
    finalize_pairwise_induced_replacement(gfa, &headers, candidate, config, "SweepGA/seqwish")
}

fn build_poa_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    let (mut graph, mut engine) = build_global_spoa_engine(config.scoring_params);
    let headers: Vec<String> = candidate
        .ranges
        .iter()
        .enumerate()
        .map(|(i, range)| format!("__impg_bubble_path{}_{}", range.path_idx, i))
        .collect();
    let sequences: Vec<String> = candidate
        .ranges
        .iter()
        .map(|range| String::from_utf8_lossy(&range.sequence).to_string())
        .collect();

    feed_sequences_to_graph(
        &mut engine,
        &mut graph,
        sequences.iter().map(|s| s.as_str()),
    );
    let gfa = graph.generate_gfa(&headers, false);
    let gfa = unchop_gfa(&gfa)?;
    let replacement = parse_gfa(&gfa)?;

    if replacement.paths.len() != candidate.ranges.len() {
        return Err(io::Error::other(format!(
            "SPOA replacement emitted {} paths for {} traversals",
            replacement.paths.len(),
            candidate.ranges.len()
        )));
    }

    for (idx, range) in candidate.ranges.iter().enumerate() {
        let observed = path_sequence(&replacement, &replacement.paths[idx])?;
        if observed != range.sequence {
            return Err(io::Error::other(format!(
                "SPOA replacement path {} changed sequence length {} -> {}",
                idx,
                range.sequence.len(),
                observed.len()
            )));
        }
    }

    Ok(replacement)
}

fn build_poasta_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    use poasta::aligner::config::{Affine2PieceMinGapCost, AffineMinGapCost};
    use poasta::aligner::scoring::{AlignmentType, GapAffine, GapAffine2Piece};
    use poasta::aligner::PoastaAligner;
    use poasta::graphs::poa::POAGraph;
    use poasta::io::graph::graph_to_gfa;

    let headers: Vec<String> = candidate
        .ranges
        .iter()
        .enumerate()
        .map(|(i, range)| format!("__impg_bubble_path{}_{}", range.path_idx, i))
        .collect();
    let sequences: Vec<&[u8]> = candidate
        .ranges
        .iter()
        .map(|range| range.sequence.as_slice())
        .collect();
    let weights = sequences
        .iter()
        .map(|sequence| vec![1usize; sequence.len()])
        .collect::<Vec<_>>();
    let (_, mismatch, gap_open, gap_extend, gap_open2, gap_extend2) = config.scoring_params;
    let mut graph = POAGraph::<u32>::new();

    if gap_extend > gap_extend2 {
        let scoring = GapAffine2Piece::new(mismatch, gap_extend, gap_open, gap_extend2, gap_open2);
        let mut aligner =
            PoastaAligner::new(Affine2PieceMinGapCost(scoring), AlignmentType::Global);
        add_poasta_sequences(&mut graph, &mut aligner, &headers, &sequences, &weights)?;
    } else {
        let scoring = GapAffine::new(mismatch, gap_extend, gap_open);
        let mut aligner = PoastaAligner::new(AffineMinGapCost(scoring), AlignmentType::Global);
        add_poasta_sequences(&mut graph, &mut aligner, &headers, &sequences, &weights)?;
    }

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
    let mut replacement = parse_gfa(&gfa)?;
    order_replacement_paths(&mut replacement, &headers)?;

    if replacement.paths.len() != candidate.ranges.len() {
        return Err(io::Error::other(format!(
            "POASTA replacement emitted {} paths for {} traversals",
            replacement.paths.len(),
            candidate.ranges.len()
        )));
    }

    for (idx, range) in candidate.ranges.iter().enumerate() {
        let observed = path_sequence(&replacement, &replacement.paths[idx])?;
        if observed != range.sequence {
            return Err(io::Error::other(format!(
                "POASTA replacement path {} changed sequence length {} -> {}",
                idx,
                range.sequence.len(),
                observed.len()
            )));
        }
    }

    Ok(replacement)
}

#[derive(Clone, Debug)]
struct StarAlignmentRow {
    inserts: Vec<Vec<u8>>,
    bases: Vec<Option<u8>>,
}

fn build_biwfa_inmemory_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
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

    if !has_empty_path {
        let polished = polish_replacement_gfa(&render_graph(&replacement), config)?;
        replacement = parse_gfa(&polished)?;
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

fn polish_replacement_gfa(gfa: &str, config: &ResolutionConfig) -> io::Result<String> {
    if config.polish_iterations == 0 {
        return Ok(gfa.to_string());
    }
    let polish_config = ResolutionConfig {
        max_iterations: config.polish_iterations,
        method: ResolutionMethod::Poa,
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
        auto_spoa_max_traversal_len: 0,
        auto_allwave_max_total_sequence: 0,
        auto_allwave_max_traversals: 0,
        scoring_params: config.scoring_params,
        ..ResolutionConfig::default()
    };
    resolve_gfa_bubbles_quiet(gfa, &polish_config).map(|resolved| resolved.gfa)
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
        if observed != range.sequence {
            return Err(io::Error::other(format!(
                "{method} replacement path {} changed sequence length {} -> {}",
                idx,
                range.sequence.len(),
                observed.len()
            )));
        }
    }
    Ok(())
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
    used_replacement: &FxHashSet<(usize, usize)>,
    out_paths: &[(String, Vec<OutStep>)],
) -> String {
    let mut out = String::new();
    out.push_str("H\tVN:Z:1.0\n");

    let mut id_by_node: FxHashMap<OutNode, String> = FxHashMap::default();
    let mut used_ids = FxHashSet::<String>::default();
    for node in used_original {
        used_ids.insert(original.segments[*node].id.clone());
    }
    let mut next_id = 1usize;

    let mut original_nodes: Vec<usize> = used_original.iter().copied().collect();
    original_nodes.sort_unstable();
    for node in original_nodes {
        let id = original.segments[node].id.clone();
        id_by_node.insert(OutNode::Original(node), id.clone());
        out.push_str(&format!(
            "S\t{}\t{}\n",
            id,
            String::from_utf8_lossy(&original.segments[node].seq)
        ));
    }

    let mut replacement_nodes: Vec<(usize, usize)> = used_replacement.iter().copied().collect();
    replacement_nodes.sort_unstable();
    for (replacement_idx, node) in replacement_nodes {
        let replacement = replacements
            .get(replacement_idx)
            .expect("used replacement index must refer to an emitted replacement graph");
        let id = next_unused_segment_id(&mut used_ids, &mut next_id);
        id_by_node.insert(OutNode::Replacement(replacement_idx, node), id.clone());
        out.push_str(&format!(
            "S\t{}\t{}\n",
            id,
            String::from_utf8_lossy(&replacement.segments[node].seq)
        ));
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
    render_rewritten_graph(
        graph,
        &[],
        &used_original,
        &FxHashSet::default(),
        &out_paths,
    )
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
        let before = seq_map(gfa);
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                method: ResolutionMethod::Poasta,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1);
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
                disable_round_quality_check: true,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 1);
    }

    #[test]
    fn auto_routes_small_bubbles_to_spoa_and_larger_bubbles_to_allwave() {
        let mut candidate = BubbleCandidate {
            ranges: Vec::new(),
            signature: "route".to_string(),
            within_budget: true,
            root_start_step: 0,
            root_end_step: 1,
            root_span: 0,
            total_steps: 0,
            unique_steps: 0,
            traversal_stats: TraversalStats {
                count: 2,
                max_len: 128,
                ..TraversalStats::default()
            },
        };
        let config = ResolutionConfig::default();
        assert_eq!(
            auto_replacement_method(&candidate, &config),
            ResolutionMethod::Poa
        );

        candidate.traversal_stats.max_len = config.auto_spoa_max_traversal_len + 1;
        candidate.traversal_stats.total_len = config.auto_allwave_max_total_sequence;
        assert_eq!(
            auto_replacement_method(&candidate, &config),
            ResolutionMethod::Allwave
        );

        candidate.traversal_stats.count = config.auto_allwave_max_traversals + 1;
        assert_eq!(
            auto_replacement_method(&candidate, &config),
            ResolutionMethod::Poa
        );

        candidate.traversal_stats.count = 2;
        candidate.traversal_stats.total_len = config.auto_allwave_max_total_sequence + 1;
        assert_eq!(
            auto_replacement_method(&candidate, &config),
            ResolutionMethod::Poa
        );

        let config = ResolutionConfig {
            auto_spoa_max_traversal_len: 0,
            auto_allwave_max_total_sequence: 0,
            auto_allwave_max_traversals: 0,
            ..ResolutionConfig::default()
        };
        candidate.traversal_stats.max_len = 1;
        assert_eq!(
            auto_replacement_method(&candidate, &config),
            ResolutionMethod::Allwave
        );
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
    fn round_quality_rejects_large_complexity_growth() {
        let before = GraphQuality {
            score: 1_000,
            ..GraphQuality::default()
        };
        let after = GraphQuality {
            score: 1_100,
            ..GraphQuality::default()
        };
        let config = ResolutionConfig {
            max_round_score_growth: 0.02,
            ..ResolutionConfig::default()
        };
        assert!(matches!(
            round_quality_decision(before, after, &config),
            RoundQualityDecision::Reject { .. }
        ));
    }

    #[test]
    fn first_crush_round_rejects_quality_regression() {
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
                max_round_score_growth: 0.0,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();

        assert_eq!(resolved.stats.resolved, 0);
        assert_eq!(resolved.stats.bailed, 1);
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(gfa, resolved.gfa);
    }

    #[test]
    fn allwave_pair_alignment_cap_bails_before_seqwish() {
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
        assert_eq!(resolved.stats.resolved, 0);
        assert_eq!(resolved.stats.bailed, 1);
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
    fn bails_out_when_candidate_exceeds_median_traversal_budget() {
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
                max_median_traversal_len: 4,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.resolved, 0);
        assert!(resolved.stats.bailed >= 1);
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
                disable_round_quality_check: true,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.stats.iterations, 1);
        assert_eq!(resolved.stats.resolved, 2);
        assert_eq!(resolved.stats.bailed, 0);
    }

    #[test]
    fn over_budget_candidates_are_not_marked_seen() {
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
        assert_eq!(resolved.stats.resolved, 1, "{:?}", resolved.stats);
        assert_eq!(resolved.stats.iterations, 2, "{:?}", resolved.stats);
        assert!(
            resolved.stats.bailed >= 2,
            "over-budget candidate should still be visible in the second round: {:?}",
            resolved.stats
        );
    }

    #[test]
    fn bails_out_when_candidate_exceeds_budget() {
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
            max_traversal_len: 4,
            ..ResolutionConfig::default()
        };
        let resolved = resolve_gfa_bubbles(gfa, &config).unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert_eq!(resolved.gfa, gfa);
        assert_eq!(resolved.stats.resolved, 0);
        assert!(resolved.stats.bailed >= 1);
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
                    sequence: Vec::new(),
                },
                PathRange {
                    path_idx: 1,
                    begin_step: 30,
                    end_step: 40,
                    sequence: Vec::new(),
                },
            ],
            signature: "selected".to_string(),
            within_budget: true,
            root_start_step: 0,
            root_end_step: 1,
            root_span: 10,
            total_steps: 20,
            unique_steps: 2,
            traversal_stats: TraversalStats::default(),
        };
        mark_candidate_occupied(&mut occupied_by_path, &selected);

        let adjacent = BubbleCandidate {
            ranges: vec![PathRange {
                path_idx: 0,
                begin_step: 20,
                end_step: 25,
                sequence: Vec::new(),
            }],
            signature: "adjacent".to_string(),
            within_budget: true,
            root_start_step: 0,
            root_end_step: 1,
            root_span: 5,
            total_steps: 5,
            unique_steps: 1,
            traversal_stats: TraversalStats::default(),
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
                sequence: Vec::new(),
            }],
            signature: "overlapping".to_string(),
            within_budget: true,
            root_start_step: 0,
            root_end_step: 1,
            root_span: 6,
            total_steps: 6,
            unique_steps: 1,
            traversal_stats: TraversalStats::default(),
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
                sequence: Vec::new(),
            }],
            signature: "different-path".to_string(),
            within_budget: true,
            root_start_step: 0,
            root_end_step: 1,
            root_span: 10,
            total_steps: 10,
            unique_steps: 1,
            traversal_stats: TraversalStats::default(),
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
}
