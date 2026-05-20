//! Bubble-guided graph resolution primitives.
//!
//! This module implements the first conservative slice of hierarchical graph
//! resolution: detect simple path-supported bubbles in a blunt GFA, replace a
//! bounded bubble with an exact path-preserving POA graph, and repeat until no
//! eligible unseen bubbles remain. It intentionally avoids lossy representative
//! collapse and coordinate sidecars; emitted paths are the coordinate system.

use crate::graph::{build_spoa_engine, feed_sequences_to_graph, reverse_complement, unchop_gfa};
use povu::native_gfa::{Step as PovuStep, Strand as PovuStrand};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::BTreeSet;
use std::io;
use std::time::Instant;

/// Configuration for exact path-preserving bubble resolution.
#[derive(Clone, Debug)]
pub struct ResolutionConfig {
    /// Maximum number of frontier replacement rounds.
    pub max_iterations: usize,
    /// Local graph induction method used for selected bubbles.
    pub method: ResolutionMethod,
    /// Maximum root path span, in bp, for a candidate bubble.
    ///
    /// The root path is the first path in the input GFA. This is a rooted POVU
    /// decomposition coordinate, not necessarily an external reference genome.
    pub max_bubble_span: usize,
    /// Maximum length of any observed traversal through the bubble.
    pub max_traversal_len: usize,
    /// Maximum median traversal length. A value of 0 disables this guard.
    pub max_median_traversal_len: usize,
    /// Maximum sum of traversal sequence lengths for one replacement.
    pub max_total_sequence: usize,
    /// Maximum number of path-supported traversals through one replacement.
    pub max_traversals: usize,
    /// SPOA scoring parameters: (match, mismatch, gap_open, gap_ext, gap_open2, gap_ext2).
    pub scoring_params: (u8, u8, u8, u8, u8, u8),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ResolutionMethod {
    Auto,
    Poa,
    Poasta,
}

impl ResolutionMethod {
    pub fn parse_name(value: &str) -> Option<Self> {
        match value.replace('_', "-").to_ascii_lowercase().as_str() {
            "auto" => Some(Self::Auto),
            "poa" | "spoa" => Some(Self::Poa),
            "poasta" => Some(Self::Poasta),
            _ => None,
        }
    }
}

pub const DEFAULT_MAX_ITERATIONS: usize = 64;
/// By default, do not cap by rooted path span.
///
/// `max_bubble_span` is a POVU root-path coordinate guard. The root is currently
/// the first GFA path, so this is intentionally not part of the default runtime
/// budget. Use traversal length/count/total sequence for the real work cap.
pub const DEFAULT_MAX_BUBBLE_SPAN: usize = 0;
pub const DEFAULT_MAX_TRAVERSAL_LEN: usize = 10_000;
pub const DEFAULT_MAX_MEDIAN_TRAVERSAL_LEN: usize = 1_000;
pub const DEFAULT_MAX_TOTAL_SEQUENCE: usize = 1_000_000;
pub const DEFAULT_MAX_TRAVERSALS: usize = 10_000;

impl Default for ResolutionConfig {
    fn default() -> Self {
        Self {
            max_iterations: DEFAULT_MAX_ITERATIONS,
            method: ResolutionMethod::Auto,
            max_bubble_span: DEFAULT_MAX_BUBBLE_SPAN,
            max_traversal_len: DEFAULT_MAX_TRAVERSAL_LEN,
            max_median_traversal_len: DEFAULT_MAX_MEDIAN_TRAVERSAL_LEN,
            max_total_sequence: DEFAULT_MAX_TOTAL_SEQUENCE,
            max_traversals: DEFAULT_MAX_TRAVERSALS,
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
    traversal_stats: TraversalStats,
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
    let mut graph = parse_gfa(gfa)?;
    log::info!(
        "crush: parsed input GFA into {} segment(s), {} path(s) in {:.2?}",
        graph.segments.len(),
        graph.paths.len(),
        parse_start.elapsed()
    );
    let mut stats = ResolutionStats::default();
    let mut seen: FxHashSet<String> = FxHashSet::default();
    let mut changed = false;

    for round in 0..config.max_iterations {
        let round_start = Instant::now();
        let discovery_start = Instant::now();
        let frontier = find_candidate_frontier(&graph, config, &seen)?;
        let discovery_elapsed = discovery_start.elapsed();
        if frontier.selected.is_empty() && frontier.bailed.is_empty() {
            log::info!(
                "crush round {}: no eligible candidates from {} POVU site(s) in {:.2?}",
                round + 1,
                frontier.sites_seen,
                discovery_elapsed
            );
            break;
        }
        stats.iterations = round + 1;

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

        for candidate in frontier.bailed {
            seen.insert(candidate.signature);
            stats.candidates_seen += 1;
            stats.bailed += 1;
        }

        if frontier.selected.is_empty() {
            continue;
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
            log::info!(
                "crush round {}: 0/{} replacement(s) built in {:.2?}",
                round + 1,
                selected_count,
                build_elapsed
            );
            continue;
        }

        let resolved_count = plans.len();
        let rewrite_start = Instant::now();
        graph = apply_replacement_frontier(&graph, &plans)?;
        let rewrite_elapsed = rewrite_start.elapsed();
        changed = true;
        stats.resolved += resolved_count;
        log::info!(
            "crush round {}: resolved {}/{} replacement(s) in {:.2?}; rewrite+validate {:.2?}; total {:.2?}",
            round + 1,
            resolved_count,
            selected_count,
            build_elapsed,
            rewrite_elapsed,
            round_start.elapsed()
        );
    }

    Ok(ResolvedGfa {
        gfa: if changed {
            render_graph(&graph)
        } else {
            gfa.to_string()
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

fn traversal_stats(ranges: &[PathRange]) -> TraversalStats {
    if ranges.is_empty() {
        return TraversalStats::default();
    }

    let mut lengths = ranges
        .iter()
        .map(|range| range.sequence.len())
        .collect::<Vec<_>>();
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
        "{label} n={}, max-len median/max={}/{}, median-len median/max={}/{}, p90-len median/max={}/{}, traversals max={}, total max={}, root-span max={}",
        candidates.len(),
        median_value(&mut max_lengths),
        max_len,
        median_value(&mut median_lengths),
        max_median,
        median_value(&mut p90_lengths),
        max_p90,
        max_traversals,
        max_total,
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
) -> io::Result<CandidateFrontier> {
    if graph.paths.is_empty() || graph.paths[0].steps.len() < 2 {
        return Ok(CandidateFrontier::default());
    }
    let root_path_idx = 0;
    let root_path = &graph.paths[root_path_idx];
    let root_positions = path_positions(graph, root_path_idx);
    let render_start = Instant::now();
    log::info!(
        "crush discovery: rendering working graph with {} segment(s), {} path(s)",
        graph.segments.len(),
        graph.paths.len()
    );
    let rendered = render_graph(graph);
    let render_elapsed = render_start.elapsed();
    let parse_start = Instant::now();
    log::info!(
        "crush discovery: parsing rendered graph with POVU ({} bytes)",
        rendered.len()
    );
    let native_graph = povu::NativeGfa::parse(&rendered).map_err(povu_to_io_error)?;
    let parse_elapsed = parse_start.elapsed();
    let reference_names = vec![root_path.name.clone()];
    let decompose_start = Instant::now();
    log::info!(
        "crush discovery: decomposing rendered graph with POVU using root '{}'",
        root_path.name
    );
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
            for path_idx in 0..graph.paths.len() {
                if let Some((range_begin, range_end)) =
                    unique_anchor_range(&graph.paths[path_idx], entry, exit)
                {
                    let mut range = PathRange {
                        path_idx,
                        begin_step: range_begin,
                        end_step: range_end,
                        sequence: Vec::new(),
                    };
                    range.sequence = range_sequence(graph, &range);
                    ranges.push(range);
                }
            }

            if ranges.len() < 2 {
                return None;
            }

            let distinct_sequences: BTreeSet<Vec<u8>> =
                ranges.iter().map(|r| r.sequence.clone()).collect();
            if distinct_sequences.len() < 2 {
                return None;
            }

            let signature = candidate_signature(
                graph,
                root_path_idx,
                &root_positions,
                begin,
                exit_step,
                &ranges,
            );
            if seen.contains(&signature) {
                return None;
            }

            let traversal_stats = traversal_stats(&ranges);
            let within_budget = (config.max_bubble_span == 0
                || root_span <= config.max_bubble_span)
                && ranges.len() <= config.max_traversals
                && traversal_stats.max_len <= config.max_traversal_len
                && (config.max_median_traversal_len == 0
                    || traversal_stats.median_len <= config.max_median_traversal_len)
                && traversal_stats.total_len <= config.max_total_sequence;

            Some(BubbleCandidate {
                ranges,
                signature,
                within_budget,
                root_start_step: begin,
                root_end_step: exit_step,
                root_span,
                traversal_stats,
            })
        })
        .collect::<Vec<_>>();
    let candidate_elapsed = candidate_start.elapsed();
    let candidates_seen = candidates.len();

    let select_start = Instant::now();
    candidates.sort_by(|a, b| {
        b.root_step_span()
            .cmp(&a.root_step_span())
            .then_with(|| b.root_span.cmp(&a.root_span))
            .then_with(|| a.root_start_step.cmp(&b.root_start_step))
            .then_with(|| a.root_end_step.cmp(&b.root_end_step))
    });

    let mut frontier = CandidateFrontier {
        sites_seen,
        candidates_seen,
        ..CandidateFrontier::default()
    };
    for candidate in candidates {
        if !candidate.within_budget {
            frontier.bailed.push(candidate);
            continue;
        }
        if frontier
            .selected
            .iter()
            .any(|selected| candidates_conflict(selected, &candidate))
        {
            continue;
        }
        frontier.selected.push(candidate);
    }
    let select_elapsed = select_start.elapsed();

    log::info!(
        "crush discovery detail: render {:.2?}, povu-parse {:.2?}, povu-decompose {:.2?}, id-map {:.2?}, candidate-build {:.2?}, select {:.2?}",
        render_elapsed,
        parse_elapsed,
        decompose_elapsed,
        id_map_elapsed,
        candidate_elapsed,
        select_elapsed
    );

    Ok(frontier)
}

impl BubbleCandidate {
    fn root_step_span(&self) -> usize {
        self.root_end_step.saturating_sub(self.root_start_step)
    }
}

fn candidates_conflict(a: &BubbleCandidate, b: &BubbleCandidate) -> bool {
    for a_range in &a.ranges {
        for b_range in &b.ranges {
            if a_range.path_idx == b_range.path_idx
                && ranges_overlap(
                    a_range.begin_step,
                    a_range.end_step,
                    b_range.begin_step,
                    b_range.end_step,
                )
            {
                return true;
            }
        }
    }
    false
}

fn ranges_overlap(a_begin: usize, a_end: usize, b_begin: usize, b_end: usize) -> bool {
    a_begin < b_end && b_begin < a_end
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

fn unique_anchor_range(path: &Path, entry: Step, exit: Step) -> Option<(usize, usize)> {
    let mut found = None;
    for (i, &step) in path.steps.iter().enumerate() {
        if step != entry {
            continue;
        }
        for j in i + 1..path.steps.len() {
            if path.steps[j] == exit {
                if found.is_some() {
                    return None;
                }
                found = Some((i, j + 1));
                break;
            }
        }
    }
    found
}

fn candidate_signature(
    graph: &Graph,
    ref_path_idx: usize,
    ref_positions: &[usize],
    begin: usize,
    exit_step: usize,
    ranges: &[PathRange],
) -> String {
    let mut seqs: Vec<String> = ranges
        .iter()
        .map(|r| String::from_utf8_lossy(&r.sequence).to_string())
        .collect();
    seqs.sort();
    format!(
        "{}:{}-{}:{}",
        graph.paths[ref_path_idx].name,
        ref_positions[begin],
        ref_positions[exit_step + 1],
        seqs.join("|")
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
        ResolutionMethod::Auto => ResolutionMethod::Poasta,
        method => method,
    };
    match method {
        ResolutionMethod::Auto => unreachable!("auto is resolved before replacement dispatch"),
        ResolutionMethod::Poa => build_poa_replacement(candidate, config),
        ResolutionMethod::Poasta => build_poasta_replacement(candidate, config).or_else(|err| {
            log::debug!("POASTA replacement failed; falling back to SPOA: {err}");
            build_poa_replacement(candidate, config)
        }),
    }
}

fn build_poa_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    let (mut graph, mut engine) = build_spoa_engine(config.scoring_params);
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
                max_iterations: 4,
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
