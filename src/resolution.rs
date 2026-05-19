//! Bubble-guided graph resolution primitives.
//!
//! This module implements the first conservative slice of hierarchical graph
//! resolution: detect simple path-supported bubbles in a blunt GFA, replace a
//! bounded bubble with an exact path-preserving POA graph, and repeat until no
//! eligible unseen bubbles remain. It intentionally avoids lossy representative
//! collapse and coordinate sidecars; emitted paths are the coordinate system.

use crate::graph::{build_spoa_engine, feed_sequences_to_graph, reverse_complement, unchop_gfa};
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::{BTreeMap, BTreeSet};
use std::io;

/// Configuration for exact path-preserving bubble resolution.
#[derive(Clone, Debug)]
pub struct ResolutionConfig {
    /// Maximum number of replacement iterations.
    pub max_iterations: usize,
    /// Maximum reference path span, in bp, for a candidate bubble.
    pub max_bubble_span: usize,
    /// Maximum length of any observed traversal through the bubble.
    pub max_traversal_len: usize,
    /// Maximum sum of traversal sequence lengths for one replacement.
    pub max_total_sequence: usize,
    /// Maximum number of path-supported traversals through one replacement.
    pub max_traversals: usize,
    /// SPOA scoring parameters: (match, mismatch, gap_open, gap_ext, gap_open2, gap_ext2).
    pub scoring_params: (u8, u8, u8, u8, u8, u8),
}

impl Default for ResolutionConfig {
    fn default() -> Self {
        Self {
            max_iterations: 8,
            max_bubble_span: 10_000,
            max_traversal_len: 10_000,
            max_total_sequence: 1_000_000,
            max_traversals: 128,
            scoring_params: (1, 4, 6, 2, 26, 1),
        }
    }
}

/// Summary of a resolution run.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ResolutionStats {
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
}

#[derive(Clone, Debug)]
struct SiteDraft {
    ref_start: usize,
    ref_end: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
enum OutNode {
    Original(usize),
    Replacement(usize),
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
    let mut graph = parse_gfa(gfa)?;
    let mut stats = ResolutionStats::default();
    let mut seen: FxHashSet<String> = FxHashSet::default();
    let mut changed = false;

    for iteration in 0..config.max_iterations {
        let Some(candidate) = find_next_candidate(&graph, config, &seen) else {
            break;
        };
        stats.iterations = iteration + 1;
        seen.insert(candidate.signature.clone());
        stats.candidates_seen += 1;

        if !candidate.within_budget {
            stats.bailed += 1;
            continue;
        }

        match replace_candidate(&graph, &candidate, config) {
            Ok(Some(next_graph)) => {
                graph = next_graph;
                changed = true;
                stats.resolved += 1;
            }
            Ok(None) => {
                stats.bailed += 1;
            }
            Err(err) => {
                log::debug!(
                    "resolution: candidate replacement failed at iteration {}: {}",
                    iteration + 1,
                    err
                );
                stats.bailed += 1;
            }
        }
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
            "H" | "W" | "C" | "J" => {}
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

fn find_next_candidate(
    graph: &Graph,
    config: &ResolutionConfig,
    seen: &FxHashSet<String>,
) -> Option<BubbleCandidate> {
    if graph.paths.is_empty() || graph.paths[0].steps.len() < 2 {
        return None;
    }
    let ref_path_idx = 0;
    let ref_path = &graph.paths[ref_path_idx];
    let ref_positions = path_positions(graph, ref_path_idx);

    for site in leaf_site_order(&candidate_sites(graph, ref_path_idx)) {
        let begin = site.ref_start;
        let exit_step = site.ref_end;
        let entry = ref_path.steps[begin];
        let exit = ref_path.steps[exit_step];
        let ref_span = ref_positions[exit_step + 1] - ref_positions[begin];

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
            continue;
        }

        let distinct_sequences: BTreeSet<Vec<u8>> =
            ranges.iter().map(|r| r.sequence.clone()).collect();
        if distinct_sequences.len() < 2 {
            continue;
        }

        let signature = candidate_signature(graph, ref_path_idx, begin, exit_step, &ranges);
        if seen.contains(&signature) {
            continue;
        }

        let max_traversal_len = ranges.iter().map(|r| r.sequence.len()).max().unwrap_or(0);
        let total_sequence = ranges.iter().map(|r| r.sequence.len()).sum::<usize>();
        let within_budget = ref_span <= config.max_bubble_span
            && ranges.len() <= config.max_traversals
            && max_traversal_len <= config.max_traversal_len
            && total_sequence <= config.max_total_sequence;

        return Some(BubbleCandidate {
            ranges,
            signature,
            within_budget,
        });
    }

    None
}

/// Discover candidate sites using the same reference-path collinear-match
/// model as POVU's native VCF extractor: every adjacent pair of shared,
/// collinear steps between the reference path and another path defines a
/// candidate if the internal traversals differ.
fn candidate_sites(graph: &Graph, ref_path_idx: usize) -> BTreeMap<(usize, usize), SiteDraft> {
    let reference = &graph.paths[ref_path_idx];
    let mut candidates = BTreeMap::<(usize, usize), SiteDraft>::new();

    for (path_idx, path) in graph.paths.iter().enumerate() {
        if path_idx == ref_path_idx {
            continue;
        }
        let matches = collinear_matches(&reference.steps, &path.steps);
        for window in matches.windows(2) {
            let (ref_start, path_start) = window[0];
            let (ref_end, path_end) = window[1];
            if ref_end <= ref_start || path_end <= path_start {
                continue;
            }
            let ref_internal = &reference.steps[ref_start + 1..ref_end];
            let alt_internal = &path.steps[path_start + 1..path_end];
            if ref_internal == alt_internal {
                continue;
            }
            candidates
                .entry((ref_start, ref_end))
                .or_insert(SiteDraft { ref_start, ref_end });
        }
    }

    candidates
}

fn collinear_matches(reference: &[Step], path: &[Step]) -> Vec<(usize, usize)> {
    let rows = reference.len();
    let cols = path.len();
    let mut dp = vec![vec![0usize; cols + 1]; rows + 1];
    for i in (0..rows).rev() {
        for j in (0..cols).rev() {
            dp[i][j] = if reference[i] == path[j] {
                1 + dp[i + 1][j + 1]
            } else {
                dp[i + 1][j].max(dp[i][j + 1])
            };
        }
    }

    let mut matches = Vec::new();
    let mut i = 0usize;
    let mut j = 0usize;
    while i < rows && j < cols {
        if reference[i] == path[j] {
            matches.push((i, j));
            i += 1;
            j += 1;
        } else if dp[i + 1][j] >= dp[i][j + 1] {
            i += 1;
        } else {
            j += 1;
        }
    }
    matches
}

fn leaf_site_order(candidates: &BTreeMap<(usize, usize), SiteDraft>) -> Vec<SiteDraft> {
    let mut leaves = candidates
        .values()
        .filter(|site| {
            !candidates.values().any(|other| {
                site.ref_start < other.ref_start && other.ref_end < site.ref_end
            })
        })
        .cloned()
        .collect::<Vec<_>>();
    leaves.sort_by_key(|site| (site.ref_end - site.ref_start, site.ref_start, site.ref_end));
    leaves
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
    begin: usize,
    exit_step: usize,
    ranges: &[PathRange],
) -> String {
    let positions = path_positions(graph, ref_path_idx);
    let mut seqs: Vec<String> = ranges
        .iter()
        .map(|r| String::from_utf8_lossy(&r.sequence).to_string())
        .collect();
    seqs.sort();
    format!(
        "{}:{}-{}:{}",
        graph.paths[ref_path_idx].name,
        positions[begin],
        positions[exit_step + 1],
        seqs.join("|")
    )
}

fn replace_candidate(
    graph: &Graph,
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Option<Graph>> {
    let replacement = build_poa_replacement(candidate, config)?;
    if replacement.segments.is_empty() {
        return Ok(None);
    }

    let mut replacement_by_path: FxHashMap<usize, Vec<OutStep>> = FxHashMap::default();
    for (range_idx, range) in candidate.ranges.iter().enumerate() {
        let Some(path) = replacement.paths.get(range_idx) else {
            return Ok(None);
        };
        let steps = path
            .steps
            .iter()
            .map(|step| OutStep {
                node: OutNode::Replacement(step.node),
                rev: step.rev,
            })
            .collect::<Vec<_>>();
        replacement_by_path.insert(range.path_idx, steps);
    }

    let range_by_path: FxHashMap<usize, (usize, usize)> = candidate
        .ranges
        .iter()
        .map(|r| (r.path_idx, (r.begin_step, r.end_step)))
        .collect();

    let mut out_paths = Vec::with_capacity(graph.paths.len());
    let mut used_original: FxHashSet<usize> = FxHashSet::default();
    let mut used_replacement: FxHashSet<usize> = FxHashSet::default();

    for (path_idx, path) in graph.paths.iter().enumerate() {
        let mut steps = Vec::new();
        let mut i = 0usize;
        while i < path.steps.len() {
            if let Some(&(begin, end)) = range_by_path.get(&path_idx) {
                if i == begin {
                    let replacement_steps =
                        replacement_by_path.get(&path_idx).ok_or_else(|| {
                            io::Error::other("replacement path missing for candidate range")
                        })?;
                    for step in replacement_steps {
                        if let OutNode::Replacement(idx) = step.node {
                            used_replacement.insert(idx);
                        }
                        steps.push(*step);
                    }
                    i = end;
                    continue;
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

    let rendered = render_rewritten_graph(
        graph,
        &replacement,
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
    Ok(Some(next))
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
    replacement: &Graph,
    used_original: &FxHashSet<usize>,
    used_replacement: &FxHashSet<usize>,
    out_paths: &[(String, Vec<OutStep>)],
) -> String {
    let mut out = String::new();
    out.push_str("H\tVN:Z:1.0\n");

    let mut id_by_node: FxHashMap<OutNode, String> = FxHashMap::default();
    let mut next_id = 1usize;

    let mut original_nodes: Vec<usize> = used_original.iter().copied().collect();
    original_nodes.sort_unstable();
    for node in original_nodes {
        let id = next_id.to_string();
        next_id += 1;
        id_by_node.insert(OutNode::Original(node), id.clone());
        out.push_str(&format!(
            "S\t{}\t{}\n",
            id,
            String::from_utf8_lossy(&original.segments[node].seq)
        ));
    }

    let mut replacement_nodes: Vec<usize> = used_replacement.iter().copied().collect();
    replacement_nodes.sort_unstable();
    for node in replacement_nodes {
        let id = next_id.to_string();
        next_id += 1;
        id_by_node.insert(OutNode::Replacement(node), id.clone());
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
        &Graph {
            segments: Vec::new(),
            paths: Vec::new(),
        },
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
    fn resolves_nested_leaf_sites_bottom_up() {
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
        let graph = parse_gfa(gfa).unwrap();
        let sites = candidate_sites(&graph, 0);
        let leaves = leaf_site_order(&sites);
        assert_eq!(
            leaves
                .iter()
                .map(|site| (site.ref_start, site.ref_end))
                .collect::<Vec<_>>(),
            vec![(1, 3)]
        );
        let resolved = resolve_gfa_bubbles(
            gfa,
            &ResolutionConfig {
                max_iterations: 4,
                ..ResolutionConfig::default()
            },
        )
        .unwrap();
        assert_eq!(before, seq_map(&resolved.gfa));
        assert!(resolved.stats.resolved >= 2);
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
