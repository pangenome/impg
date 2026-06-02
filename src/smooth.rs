//! Graph smoothing module implementing smoothxg-style block decomposition and POA smoothing.
//!
//! Takes a raw GFA (e.g., from seqwish), decomposes it into collinear blocks,
//! runs SPOA partial-order alignment on each block, and reassembles the smoothed
//! graph via lacing.

use crate::graph::{
    build_spoa_engine, feed_sequences_to_graph, reverse_complement, sort_gfa, unchop_gfa,
};
use bitvec::{bitvec, prelude::BitVec};
use log::info;
use rustc_hash::{FxHashMap, FxHashSet};
use std::{collections::BTreeSet, io, time::Instant};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for graph smoothing.
pub struct SmoothConfig {
    /// Number of haplotypes (used for per-pass max_block_weight = target * n_haps).
    pub n_haps: usize,
    /// Target POA length(s) in bp — one value per smoothing pass (default: [700, 1100],
    /// matching pggb's default `-G 700,1100`).  Each pass re-decomposes the graph with
    /// that target length and feeds its output into the next pass.
    pub target_poa_lengths: Vec<usize>,
    /// Chop nodes to this maximum length (default: 100).
    pub max_node_length: usize,
    /// Flanking sequence fraction for POA padding (default: 0.001).
    pub poa_padding_fraction: f64,
    /// Depth threshold for minimum ~311bp padding (default: 100).
    pub max_block_depth_for_padding_more: usize,
    /// SPOA scoring parameters: (match, mismatch, gap_open, gap_ext, gap_open2, gap_ext2).
    pub scoring_params: (u8, u8, u8, u8, u8, u8),
    /// Number of threads for parallel block processing.
    pub num_threads: usize,
    /// Optional temp directory for intermediate files.
    pub temp_dir: Option<String>,
    /// Skip the initial unchop+sort step (use when input is already unchopped and sorted).
    pub pre_sorted: bool,
    /// Source used to place smoothing blocks.
    pub block_source: SmoothBlockSource,
    /// Optional POVU reference path hints used by flubble-guided block placement.
    pub flubble_reference_names: Vec<String>,
}

/// Block placement strategy for smoothing.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SmoothBlockSource {
    /// smoothxg-style path-overlap decomposition.
    PathOverlap,
    /// POVU flubble extents (`level >= 1`) plus passthrough gap blocks.
    Flubble,
    /// POVU bubble neighborhoods merged greedily along the reference path and
    /// realigned with POASTA as local blocks.
    NeighborMergePoasta,
}

impl SmoothConfig {
    /// Create a new SmoothConfig with defaults derived from n_haps.
    /// Defaults match pggb's smoothxg invocation: two passes at 700 bp then 1100 bp
    /// (`-G 700,1100`), asm20 POA params, pad_max_depth=100.
    pub fn new(n_haps: usize) -> Self {
        SmoothConfig {
            n_haps,
            target_poa_lengths: vec![700, 1100],
            max_node_length: 100,
            poa_padding_fraction: 0.001,
            max_block_depth_for_padding_more: 100,
            scoring_params: (1, 4, 6, 2, 26, 1),
            num_threads: 4,
            temp_dir: None,
            pre_sorted: false,
            block_source: SmoothBlockSource::PathOverlap,
            flubble_reference_names: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// Internal graph representation
// ---------------------------------------------------------------------------

struct SmoothGraph {
    /// Node sequences, indexed by position in sorted GFA order.
    nodes: Vec<Vec<u8>>,
    /// Edges: set of (from_idx, from_rev, to_idx, to_rev).
    edges: FxHashSet<(usize, bool, usize, bool)>,
    /// Paths through the graph.
    paths: Vec<PathData>,
}

struct PathData {
    name: String,
    /// Steps: (node_idx, is_reverse).
    steps: Vec<(usize, bool)>,
}

// ---------------------------------------------------------------------------
// Block types
// ---------------------------------------------------------------------------

/// A sequence entry prepared for SPOA with padding metadata.
struct SeqEntry {
    sequence: String,
    core_sequence: String,
    path_name: String,
    left_pad: usize,
    right_pad: usize,
}

#[derive(Clone)]
struct PathRange {
    path_idx: usize,
    begin_step: usize,
    end_step: usize, // exclusive
    length: usize,   // total base pairs
}

#[derive(Clone)]
struct Block {
    path_ranges: Vec<PathRange>,
    source: BlockSource,
}

#[derive(Clone)]
enum BlockSource {
    PathOverlap,
    Flubble {
        site_id: String,
        level: usize,
        reference_start_step: usize,
        reference_end_step: usize,
        is_leaf: bool,
    },
    NeighborMerge {
        site_count: usize,
        reference_start_step: usize,
        reference_end_step: usize,
        reference_span_bp: usize,
    },
    Gap,
}

impl BlockSource {
    fn is_gap(&self) -> bool {
        matches!(self, Self::Gap)
    }

    fn uses_poasta(&self) -> bool {
        matches!(self, Self::NeighborMerge { .. })
    }
}

// ---------------------------------------------------------------------------
// Union-Find for topological split
// ---------------------------------------------------------------------------

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    fn union(&mut self, x: usize, y: usize) {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx != ry {
            if self.rank[rx] < self.rank[ry] {
                self.parent[rx] = ry;
            } else if self.rank[rx] > self.rank[ry] {
                self.parent[ry] = rx;
            } else {
                self.parent[ry] = rx;
                self.rank[rx] += 1;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Smooth a GFA graph using smoothxg-style block decomposition and POA.
///
/// Runs one pass per entry in `config.target_poa_lengths` (default: `[700, 1100]`,
/// matching pggb's `-G 700,1100`).  Each pass feeds its output into the next pass,
/// progressively resolving longer variation.
///
/// Per-pass pipeline: unchop+sort → chop → block decompose → break long blocks →
/// per-block SPOA → lace.
pub fn smooth_gfa(gfa_content: &str, config: &SmoothConfig) -> io::Result<String> {
    let n_passes = config.target_poa_lengths.len();
    if n_passes == 0 {
        return Ok(gfa_content.to_string());
    }

    let mut current = gfa_content.to_string();
    for (pass_idx, &target_poa_length) in config.target_poa_lengths.iter().enumerate() {
        let is_first = pass_idx == 0;
        let pass_start = Instant::now();
        let before_stats = gfa_basic_stats(&current);
        info!(
            "[smooth] Pass {}/{} (target_poa_length={}bp)",
            pass_idx + 1,
            n_passes,
            target_poa_length
        );
        let result = smooth_gfa_pass(
            &current,
            target_poa_length,
            is_first && config.pre_sorted,
            config,
        )?;
        current = result.gfa;
        if let Some(metrics) = result.neighbor_metrics {
            let after_stats = gfa_basic_stats(&current);
            info!(
                "[neighbor-merge-poasta] iteration {}/{}: block-size-cap={}bp, blocks-merged={} ({} site(s) in merged blocks), neighbor-blocks={}, gap-blocks={}, selected-sites={}, skipped-overlap={}, skipped-over-cap={}, segs-before={}, segs-after={}, segment-bp-before={}, segment-bp-after={}, wall={:.3}s",
                pass_idx + 1,
                n_passes,
                metrics.block_size_cap,
                metrics.blocks_merged,
                metrics.sites_in_merged_blocks,
                metrics.neighbor_blocks,
                metrics.gap_blocks,
                metrics.selected_sites,
                metrics.overlap_skipped_sites,
                metrics.over_cap_skipped_sites,
                before_stats.segments,
                after_stats.segments,
                before_stats.segment_bp,
                after_stats.segment_bp,
                pass_start.elapsed().as_secs_f64(),
            );
        }
    }
    Ok(current)
}

#[derive(Clone, Copy, Debug, Default)]
struct GfaBasicStats {
    segments: usize,
    segment_bp: usize,
}

fn gfa_basic_stats(gfa: &str) -> GfaBasicStats {
    let mut stats = GfaBasicStats::default();
    for line in gfa.lines() {
        let mut fields = line.split('\t');
        if fields.next() != Some("S") {
            continue;
        }
        fields.next();
        if let Some(seq) = fields.next() {
            stats.segments += 1;
            stats.segment_bp += seq.len();
        }
    }
    stats
}

#[derive(Clone, Copy, Debug, Default)]
struct NeighborMergeMetrics {
    block_size_cap: usize,
    selected_sites: usize,
    overlap_skipped_sites: usize,
    over_cap_skipped_sites: usize,
    blocks_merged: usize,
    sites_in_merged_blocks: usize,
    neighbor_blocks: usize,
    gap_blocks: usize,
}

struct SmoothPassResult {
    gfa: String,
    neighbor_metrics: Option<NeighborMergeMetrics>,
}

/// Run a single smoothing pass on a (not yet sorted) GFA string.
///
/// `pre_sorted` — skip the unchop+sort step (true only for the first pass when the
/// caller guarantees the input is already sorted).
fn smooth_gfa_pass(
    gfa_content: &str,
    target_poa_length: usize,
    pre_sorted: bool,
    config: &SmoothConfig,
) -> io::Result<SmoothPassResult> {
    let max_block_weight = target_poa_length * config.n_haps.max(1);
    let max_poa_length = 2 * target_poa_length;

    // Step 1: Unchop then sort input GFA (gives a 1D layout for block decomposition)
    let sorted_gfa = if pre_sorted {
        gfa_content.to_string()
    } else {
        sort_gfa(&unchop_gfa(gfa_content)?, config.num_threads)?
    };

    // Step 2: Parse sorted GFA
    let mut graph = parse_gfa(&sorted_gfa);
    let num_nodes_before_chop = graph.nodes.len();
    let num_paths = graph.paths.len();
    info!(
        "[smooth] Parsed sorted GFA: {} nodes, {} paths",
        num_nodes_before_chop, num_paths
    );

    if graph.nodes.is_empty() || graph.paths.is_empty() {
        return Ok(SmoothPassResult {
            gfa: gfa_content.to_string(),
            neighbor_metrics: None,
        });
    }

    // Step 3: Chop nodes to max_node_length
    chop_graph(&mut graph, config.max_node_length);
    info!(
        "[smooth] Chopped: {} → {} nodes (max {}bp)",
        num_nodes_before_chop,
        graph.nodes.len(),
        config.max_node_length
    );

    // Step 4: Build inverse index (node → path steps)
    let node_step_index = build_node_step_index(&graph);

    // Step 5: Compute cumulative path positions
    let path_positions = compute_path_positions(&graph);

    // Step 6: Block decomposition
    let mut neighbor_metrics = None;
    let mut blocks = match config.block_source {
        SmoothBlockSource::PathOverlap => {
            let blocks = smoothable_blocks(
                &graph,
                &node_step_index,
                max_block_weight,
                target_poa_length,
            );
            info!(
                "[smooth] Decomposed into {} path-overlap blocks",
                blocks.len()
            );
            blocks
        }
        SmoothBlockSource::Flubble => {
            let plan = flubble_guided_blocks(&graph, &path_positions, config)?;
            info!(
                "[smooth] POVU flubble-guided placement: {} site(s), {} level>=1 candidate(s), {} selected flubble-block(s), {} overlap-skipped, {} gap block(s)",
                plan.povu_sites,
                plan.level_ge_1_sites,
                plan.selected_flubble_blocks,
                plan.overlap_skipped_sites,
                plan.gap_blocks
            );
            plan.blocks
        }
        SmoothBlockSource::NeighborMergePoasta => {
            let plan = neighbor_merge_blocks(&graph, &path_positions, config, target_poa_length)?;
            info!(
                "[neighbor-merge-poasta] placement: {} POVU site(s), {} candidate site(s), {} selected site(s), {} selected group(s), {} neighbor block(s), {} merged block(s), {} site(s) inside merged blocks, {} gap block(s), {} overlap-skipped, {} over-cap-skipped",
                plan.povu_sites,
                plan.candidate_sites,
                plan.selected_sites,
                plan.selected_groups,
                plan.neighbor_blocks,
                plan.blocks_merged,
                plan.sites_in_merged_blocks,
                plan.gap_blocks,
                plan.overlap_skipped_sites,
                plan.over_cap_skipped_sites,
            );
            neighbor_metrics = Some(NeighborMergeMetrics {
                block_size_cap: target_poa_length,
                selected_sites: plan.selected_sites,
                overlap_skipped_sites: plan.overlap_skipped_sites,
                over_cap_skipped_sites: plan.over_cap_skipped_sites,
                blocks_merged: plan.blocks_merged,
                sites_in_merged_blocks: plan.sites_in_merged_blocks,
                neighbor_blocks: plan.neighbor_blocks,
                gap_blocks: plan.gap_blocks,
            });
            plan.blocks
        }
    };

    // Step 7: Break long blocks
    break_blocks(&mut blocks, &graph, max_poa_length);
    info!(
        "[smooth] After breaking: {} blocks, {} total path ranges",
        blocks.len(),
        blocks.iter().map(|b| b.path_ranges.len()).sum::<usize>()
    );
    if config.block_source == SmoothBlockSource::Flubble {
        let processing_flubble_blocks = blocks
            .iter()
            .filter(|block| matches!(block.source, BlockSource::Flubble { .. }))
            .count();
        let processing_gap_blocks = blocks
            .iter()
            .filter(|block| matches!(block.source, BlockSource::Gap))
            .count();
        info!(
            "[smooth] Processing {} flubble block(s) and {} gap block(s) after size splitting",
            processing_flubble_blocks, processing_gap_blocks,
        );
    }

    // Step 8: Per-block SPOA/POASTA smoothing (parallelized).
    //
    // Gap blocks are passthrough by construction: they are the coverage
    // complement between selected smoothable blocks. Smoothable blocks must
    // return an induced graph; substituting the original subgraph after an
    // aligner failure would silently change the requested smoothing semantics.
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.num_threads)
        .build()
        .map_err(|e| io::Error::other(format!("failed to build thread pool: {}", e)))?;
    let n_smoothed = AtomicUsize::new(0);
    let n_gap_passthrough = AtomicUsize::new(0);
    let block_results: Vec<io::Result<Option<String>>> = pool.install(|| {
        blocks
            .par_iter()
            .enumerate()
            .map(|(i, block)| -> io::Result<Option<String>> {
                let input_bp = block_input_bp(block);
                if block.source.is_gap() {
                    let s = passthrough_block_gfa(block, &graph, &path_positions)?;
                    if s.is_empty() {
                        return Ok(None);
                    }
                    n_gap_passthrough.fetch_add(1, Ordering::Relaxed);
                    log_flubble_block_compression(i, block, input_bp, gfa_segment_bp(&s), true);
                    return Ok(Some(s));
                }

                let smooth_result = if block.source.uses_poasta() {
                    poasta_block(block, &graph, config, &path_positions)
                } else {
                    smooth_block(block, &graph, config, &path_positions)
                };

                match smooth_result {
                    Ok(s) if !s.is_empty() => {
                        n_smoothed.fetch_add(1, Ordering::Relaxed);
                        log_flubble_block_compression(
                            i,
                            block,
                            input_bp,
                            gfa_segment_bp(&s),
                            false,
                        );
                        Ok(Some(s))
                    }
                    Ok(_) => Err(io::Error::other(smooth_empty_block_error(
                        i, block, input_bp,
                    ))),
                    Err(e) => Err(io::Error::other(format!(
                        "smooth block {i} failed; refusing passthrough fallback: {e}"
                    ))),
                }
            })
            .collect()
    });
    let mut block_gfas = Vec::with_capacity(block_results.len());
    for result in block_results {
        if let Some(gfa) = result? {
            block_gfas.push(gfa);
        }
    }

    info!(
        "[smooth] Smoothed {} blocks ({} gap passthrough)",
        n_smoothed.load(Ordering::Relaxed),
        n_gap_passthrough.load(Ordering::Relaxed),
    );

    if block_gfas.is_empty() {
        return Ok(SmoothPassResult {
            gfa: String::from("H\tVN:Z:1.0\n"),
            neighbor_metrics,
        });
    }

    // Step 9: Lace block GFAs together
    // Per-block gfaffix already normalized each block; the caller's normalize_and_sort
    // (gfaffix + sort) handles any remaining redundancies introduced at lacing boundaries.
    let laced = crate::commands::lace::lace_subgraphs(&block_gfas, config.temp_dir.as_deref())?;
    info!("[smooth] Lacing complete");

    Ok(SmoothPassResult {
        gfa: laced,
        neighbor_metrics,
    })
}

// ---------------------------------------------------------------------------
// GFA parsing
// ---------------------------------------------------------------------------

fn parse_gfa(gfa: &str) -> SmoothGraph {
    let mut nodes: Vec<Vec<u8>> = Vec::new();
    let mut node_id_to_idx: FxHashMap<u64, usize> = FxHashMap::default();
    let mut edges: FxHashSet<(usize, bool, usize, bool)> = FxHashSet::default();
    let mut paths: Vec<PathData> = Vec::new();

    for line in gfa.lines() {
        if line.is_empty() {
            continue;
        }
        let first_byte = line.as_bytes()[0];
        match first_byte {
            b'S' => {
                let fields: Vec<&str> = line.splitn(4, '\t').collect();
                if fields.len() >= 3 {
                    if let Ok(node_id) = fields[1].parse::<u64>() {
                        let idx = nodes.len();
                        nodes.push(fields[2].as_bytes().to_vec());
                        node_id_to_idx.insert(node_id, idx);
                    }
                }
            }
            b'L' => {
                let fields: Vec<&str> = line.splitn(7, '\t').collect();
                if fields.len() >= 6 {
                    if let (Ok(from_id), Ok(to_id)) =
                        (fields[1].parse::<u64>(), fields[3].parse::<u64>())
                    {
                        if let (Some(&from_idx), Some(&to_idx)) =
                            (node_id_to_idx.get(&from_id), node_id_to_idx.get(&to_id))
                        {
                            let from_rev = fields[2] == "-";
                            let to_rev = fields[4] == "-";
                            edges.insert((from_idx, from_rev, to_idx, to_rev));
                        }
                    }
                }
            }
            b'P' => {
                let fields: Vec<&str> = line.splitn(4, '\t').collect();
                if fields.len() >= 3 {
                    let name = fields[1].to_string();
                    let steps: Vec<(usize, bool)> = fields[2]
                        .split(',')
                        .filter_map(|step| {
                            let (id_str, is_rev) = if let Some(stripped) = step.strip_suffix('+') {
                                (stripped, false)
                            } else if let Some(stripped) = step.strip_suffix('-') {
                                (stripped, true)
                            } else {
                                (step, false)
                            };
                            let node_id: u64 = id_str.parse().ok()?;
                            let idx = *node_id_to_idx.get(&node_id)?;
                            Some((idx, is_rev))
                        })
                        .collect();
                    if !steps.is_empty() {
                        paths.push(PathData { name, steps });
                    }
                }
            }
            _ => {}
        }
    }

    SmoothGraph {
        nodes,
        edges,
        paths,
    }
}

// ---------------------------------------------------------------------------
// Node chopping
// ---------------------------------------------------------------------------

fn chop_graph(graph: &mut SmoothGraph, max_node_length: usize) {
    if max_node_length == 0 {
        return;
    }

    // Check if any node needs chopping
    let needs_chop = graph.nodes.iter().any(|n| n.len() > max_node_length);
    if !needs_chop {
        return;
    }

    let mut new_nodes: Vec<Vec<u8>> = Vec::new();
    // replacement_map[old_idx] = vec of new indices (in forward order)
    let mut replacement_map: Vec<Vec<usize>> = Vec::with_capacity(graph.nodes.len());

    for node_seq in graph.nodes.iter() {
        if node_seq.len() <= max_node_length {
            let new_idx = new_nodes.len();
            new_nodes.push(node_seq.clone());
            replacement_map.push(vec![new_idx]);
        } else {
            let pieces: Vec<usize> = node_seq
                .chunks(max_node_length)
                .map(|chunk| {
                    let idx = new_nodes.len();
                    new_nodes.push(chunk.to_vec());
                    idx
                })
                .collect();
            replacement_map.push(pieces);
        }
    }

    // Update paths
    for path in graph.paths.iter_mut() {
        let mut new_steps = Vec::new();
        for &(node_idx, is_reverse) in &path.steps {
            let pieces = &replacement_map[node_idx];
            if is_reverse {
                for &piece_idx in pieces.iter().rev() {
                    new_steps.push((piece_idx, true));
                }
            } else {
                for &piece_idx in pieces.iter() {
                    new_steps.push((piece_idx, false));
                }
            }
        }
        path.steps = new_steps;
    }

    // Update edges
    let mut new_edges: FxHashSet<(usize, bool, usize, bool)> = FxHashSet::default();
    for &(from_idx, from_rev, to_idx, to_rev) in &graph.edges {
        let from_pieces = &replacement_map[from_idx];
        let to_pieces = &replacement_map[to_idx];

        // Source piece: rightmost in forward, leftmost in reverse
        let from_piece = if from_rev {
            from_pieces[0]
        } else {
            *from_pieces.last().unwrap()
        };

        // Target piece: leftmost in forward, rightmost in reverse
        let to_piece = if to_rev {
            *to_pieces.last().unwrap()
        } else {
            to_pieces[0]
        };

        new_edges.insert((from_piece, from_rev, to_piece, to_rev));
    }

    // Add chain edges for chopped nodes
    for pieces in &replacement_map {
        for i in 0..pieces.len().saturating_sub(1) {
            new_edges.insert((pieces[i], false, pieces[i + 1], false));
        }
    }

    graph.nodes = new_nodes;
    graph.edges = new_edges;
}

// ---------------------------------------------------------------------------
// Index building
// ---------------------------------------------------------------------------

/// Build inverse index: node_idx → [(path_idx, step_idx)].
fn build_node_step_index(graph: &SmoothGraph) -> Vec<Vec<(usize, usize)>> {
    let mut index: Vec<Vec<(usize, usize)>> = vec![Vec::new(); graph.nodes.len()];
    for (path_idx, path) in graph.paths.iter().enumerate() {
        for (step_idx, &(node_idx, _)) in path.steps.iter().enumerate() {
            index[node_idx].push((path_idx, step_idx));
        }
    }
    index
}

/// Compute cumulative path positions: positions[path_idx][step_idx] = bp offset.
/// positions[path_idx] has length steps.len() + 1 (last entry = total path length).
fn compute_path_positions(graph: &SmoothGraph) -> Vec<Vec<usize>> {
    graph
        .paths
        .iter()
        .map(|path| {
            let mut positions = Vec::with_capacity(path.steps.len() + 1);
            positions.push(0);
            for &(node_idx, _) in &path.steps {
                let prev = *positions.last().unwrap();
                positions.push(prev + graph.nodes[node_idx].len());
            }
            positions
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Block decomposition (smoothxg smoothable_blocks)
// ---------------------------------------------------------------------------

fn smoothable_blocks(
    graph: &SmoothGraph,
    node_step_index: &[Vec<(usize, usize)>],
    max_block_weight: usize,
    target_poa_length: usize,
) -> Vec<Block> {
    let mut seen_steps: Vec<BitVec> = graph
        .paths
        .iter()
        .map(|p| bitvec![0; p.steps.len()])
        .collect();

    let mut blocks: Vec<Block> = Vec::new();

    // Current block accumulation
    let mut block_node_indices: Vec<usize> = Vec::new();
    let mut total_path_length: usize = 0;
    // path_idx → (total_bp, step_count) for the current block
    let mut path_coverage: FxHashMap<usize, (usize, usize)> = FxHashMap::default();
    let mut block_handles_count: usize = 0;

    for node_idx in 0..graph.nodes.len() {
        let node_len = graph.nodes[node_idx].len();

        // Compute sequence_to_add: sum of node_len for each unseen step
        let mut sequence_to_add: usize = 0;
        for &(path_idx, step_idx) in &node_step_index[node_idx] {
            if !seen_steps[path_idx][step_idx] {
                sequence_to_add += node_len;
            }
        }

        // Estimate max path length from path_coverage
        let max_path_length = if block_handles_count > 0 {
            path_coverage
                .values()
                .map(|&(total_bp, count)| {
                    if count > 0 {
                        total_bp * block_handles_count / count
                    } else {
                        0
                    }
                })
                .max()
                .unwrap_or(0)
        } else {
            0
        };

        // Check if current block should be finalized
        let should_finalize = !block_node_indices.is_empty()
            && (total_path_length + sequence_to_add > max_block_weight
                || max_path_length > target_poa_length);

        if should_finalize {
            let new_blocks =
                finalize_block(&block_node_indices, graph, node_step_index, &mut seen_steps);
            blocks.extend(new_blocks);

            block_node_indices.clear();
            total_path_length = 0;
            path_coverage.clear();
            block_handles_count = 0;
        }

        // Add node to current block
        block_node_indices.push(node_idx);
        total_path_length += sequence_to_add;
        block_handles_count += 1;

        // Update path_coverage
        for &(path_idx, step_idx) in &node_step_index[node_idx] {
            if !seen_steps[path_idx][step_idx] {
                let entry = path_coverage.entry(path_idx).or_insert((0, 0));
                entry.0 += node_len;
                entry.1 += 1;
            }
        }
    }

    // Finalize last block
    if !block_node_indices.is_empty() {
        let new_blocks =
            finalize_block(&block_node_indices, graph, node_step_index, &mut seen_steps);
        blocks.extend(new_blocks);
    }

    blocks
}

fn finalize_block(
    block_node_indices: &[usize],
    graph: &SmoothGraph,
    node_step_index: &[Vec<(usize, usize)>],
    seen_steps: &mut [BitVec],
) -> Vec<Block> {
    // 1. Collect all unseen step traversals on block handles
    let mut traversals: Vec<(usize, usize)> = Vec::new();
    for &node_idx in block_node_indices {
        for &(path_idx, step_idx) in &node_step_index[node_idx] {
            if !seen_steps[path_idx][step_idx] {
                traversals.push((path_idx, step_idx));
            }
        }
    }

    if traversals.is_empty() {
        return Vec::new();
    }

    // 2. Sort by (path_idx, step_idx)
    traversals.sort();

    // 3. Build path ranges, breaking when:
    //    - path changes
    //    - step_idx is not consecutive (gap; equivalent to max_path_jump=0)
    let mut path_ranges: Vec<PathRange> = Vec::new();
    let mut range_path: usize = traversals[0].0;
    let mut range_begin: usize = traversals[0].1;
    let mut prev_step: usize = traversals[0].1;

    for &(path_idx, step_idx) in traversals.iter().skip(1) {
        if path_idx != range_path || step_idx != prev_step + 1 {
            // Finalize current range
            let length = compute_range_length(graph, range_path, range_begin, prev_step + 1);
            if length > 0 {
                path_ranges.push(PathRange {
                    path_idx: range_path,
                    begin_step: range_begin,
                    end_step: prev_step + 1,
                    length,
                });
            }
            range_path = path_idx;
            range_begin = step_idx;
        }
        prev_step = step_idx;
    }

    // Finalize last range
    let length = compute_range_length(graph, range_path, range_begin, prev_step + 1);
    if length > 0 {
        path_ranges.push(PathRange {
            path_idx: range_path,
            begin_step: range_begin,
            end_step: prev_step + 1,
            length,
        });
    }

    if path_ranges.is_empty() {
        return Vec::new();
    }

    // 4. Mark kept steps as seen
    for range in &path_ranges {
        for step_idx in range.begin_step..range.end_step {
            seen_steps[range.path_idx].set(step_idx, true);
        }
    }

    // 5. Sort by length (longest first, for better SPOA quality)
    path_ranges.sort_by(|a, b| b.length.cmp(&a.length));

    // 6. Topological split: split disconnected components into separate blocks
    topological_split(path_ranges, graph)
}

/// Compute total bp length of a path range.
fn compute_range_length(
    graph: &SmoothGraph,
    path_idx: usize,
    begin_step: usize,
    end_step: usize,
) -> usize {
    let path = &graph.paths[path_idx];
    let mut length = 0;
    for step_idx in begin_step..end_step {
        let (node_idx, _) = path.steps[step_idx];
        length += graph.nodes[node_idx].len();
    }
    length
}

/// Split path_ranges into separate blocks for disconnected graph components.
fn topological_split(path_ranges: Vec<PathRange>, graph: &SmoothGraph) -> Vec<Block> {
    topological_split_with_source(path_ranges, graph, BlockSource::PathOverlap)
}

/// Split path_ranges into separate blocks for disconnected graph components.
fn topological_split_with_source(
    path_ranges: Vec<PathRange>,
    graph: &SmoothGraph,
    source: BlockSource,
) -> Vec<Block> {
    if path_ranges.len() <= 1 {
        return vec![Block {
            path_ranges,
            source,
        }];
    }

    // Collect all nodes referenced in path_ranges
    let mut node_set: FxHashSet<usize> = FxHashSet::default();
    for range in &path_ranges {
        for step_idx in range.begin_step..range.end_step {
            let (node_idx, _) = graph.paths[range.path_idx].steps[step_idx];
            node_set.insert(node_idx);
        }
    }

    let nodes: Vec<usize> = node_set.into_iter().collect();
    if nodes.len() <= 1 {
        return vec![Block {
            path_ranges,
            source,
        }];
    }

    let node_to_local: FxHashMap<usize, usize> =
        nodes.iter().enumerate().map(|(i, &n)| (n, i)).collect();

    let mut uf = UnionFind::new(nodes.len());

    // Union consecutive nodes in each path_range
    for range in &path_ranges {
        let mut prev_local: Option<usize> = None;
        for step_idx in range.begin_step..range.end_step {
            let (node_idx, _) = graph.paths[range.path_idx].steps[step_idx];
            let local = node_to_local[&node_idx];
            if let Some(prev) = prev_local {
                uf.union(prev, local);
            }
            prev_local = Some(local);
        }
    }

    // Group path_ranges by component
    let mut component_ranges: FxHashMap<usize, Vec<PathRange>> = FxHashMap::default();
    for range in path_ranges {
        let first_step = range.begin_step;
        let (node_idx, _) = graph.paths[range.path_idx].steps[first_step];
        let local = node_to_local[&node_idx];
        let root = uf.find(local);
        component_ranges.entry(root).or_default().push(range);
    }

    if component_ranges.len() == 1 {
        return component_ranges
            .into_values()
            .map(|path_ranges| Block {
                path_ranges,
                source: source.clone(),
            })
            .collect();
    }

    // Return one block per component
    component_ranges
        .into_values()
        .filter(|ranges| !ranges.is_empty())
        .map(|path_ranges| Block {
            path_ranges,
            source: source.clone(),
        })
        .collect()
}

struct FlubbleGuidedBlockPlan {
    blocks: Vec<Block>,
    povu_sites: usize,
    level_ge_1_sites: usize,
    selected_flubble_blocks: usize,
    overlap_skipped_sites: usize,
    gap_blocks: usize,
}

#[derive(Clone)]
struct FlubbleBlockCandidate {
    source: BlockSource,
    path_ranges: Vec<PathRange>,
    level: usize,
    reference_start_step: usize,
    reference_end_step: usize,
    reference_span_bp: usize,
}

struct NeighborMergeBlockPlan {
    blocks: Vec<Block>,
    povu_sites: usize,
    candidate_sites: usize,
    selected_sites: usize,
    selected_groups: usize,
    neighbor_blocks: usize,
    blocks_merged: usize,
    sites_in_merged_blocks: usize,
    overlap_skipped_sites: usize,
    over_cap_skipped_sites: usize,
    gap_blocks: usize,
}

#[derive(Clone)]
struct NeighborSiteCandidate {
    entry: (usize, bool),
    exit: (usize, bool),
    level: usize,
    reference_start_step: usize,
    reference_end_step: usize,
    reference_span_bp: usize,
}

struct NeighborSiteGroup {
    sites: Vec<NeighborSiteCandidate>,
    reference_start_step: usize,
    reference_end_step: usize,
    reference_span_bp: usize,
}

fn flubble_guided_blocks(
    graph: &SmoothGraph,
    path_positions: &[Vec<usize>],
    config: &SmoothConfig,
) -> io::Result<FlubbleGuidedBlockPlan> {
    if graph.paths.is_empty() {
        return Ok(FlubbleGuidedBlockPlan {
            blocks: Vec::new(),
            povu_sites: 0,
            level_ge_1_sites: 0,
            selected_flubble_blocks: 0,
            overlap_skipped_sites: 0,
            gap_blocks: 0,
        });
    }

    let rendered = render_smooth_graph_gfa(graph)?;
    let native = povu::NativeGfa::parse(&rendered).map_err(povu_to_io_error)?;
    let reference_names = if config.flubble_reference_names.is_empty() {
        vec![graph.paths[0].name.clone()]
    } else {
        config.flubble_reference_names.clone()
    };
    let decomposition = match native.decompose_flubbles(&reference_names) {
        Ok(decomposition) => decomposition,
        Err(err) if !config.flubble_reference_names.is_empty() => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "POVU flubble reference hint(s) {:?} did not match graph paths: {}",
                    config.flubble_reference_names, err
                ),
            ));
        }
        Err(err) => return Err(povu_to_io_error(err)),
    };
    let root_positions = path_positions
        .get(decomposition.reference_path_index)
        .ok_or_else(|| {
            io::Error::other("POVU reference path index is outside smooth graph paths")
        })?;

    let id_to_idx = (0..graph.nodes.len())
        .map(|idx| ((idx + 1).to_string(), idx))
        .collect::<FxHashMap<_, _>>();
    let path_step_indexes = graph
        .paths
        .iter()
        .map(smooth_path_step_index)
        .collect::<Vec<_>>();

    let mut candidates = Vec::new();
    let mut level_ge_1_sites = 0usize;
    for site in &decomposition.sites {
        if site.level < 1 {
            continue;
        }
        level_ge_1_sites += 1;
        let Some(entry) = smooth_step_from_povu(&id_to_idx, &site.start) else {
            continue;
        };
        let Some(exit) = smooth_step_from_povu(&id_to_idx, &site.end) else {
            continue;
        };
        let mut path_ranges = Vec::new();
        for (path_idx, path_index) in path_step_indexes.iter().enumerate() {
            if let Some((begin_step, end_step)) =
                unique_smooth_anchor_range(path_index, entry, exit)
            {
                let length =
                    path_positions[path_idx][end_step] - path_positions[path_idx][begin_step];
                if length > 0 {
                    path_ranges.push(PathRange {
                        path_idx,
                        begin_step,
                        end_step,
                        length,
                    });
                }
            }
        }
        if path_ranges.len() < 2 {
            continue;
        }
        path_ranges.sort_by(|a, b| b.length.cmp(&a.length));
        let reference_span_bp = if site.reference_end_step + 1 < root_positions.len()
            && site.reference_start_step < root_positions.len()
        {
            root_positions[site.reference_end_step + 1] - root_positions[site.reference_start_step]
        } else {
            0
        };
        candidates.push(FlubbleBlockCandidate {
            source: BlockSource::Flubble {
                site_id: site.id.clone(),
                level: site.level,
                reference_start_step: site.reference_start_step,
                reference_end_step: site.reference_end_step,
                is_leaf: site.is_leaf,
            },
            path_ranges,
            level: site.level,
            reference_start_step: site.reference_start_step,
            reference_end_step: site.reference_end_step,
            reference_span_bp,
        });
    }

    // Prefer the most local natural boundaries first. If a parent and child
    // flubble overlap, this keeps the deeper child and leaves the parent span to
    // gap passthrough around it, avoiding duplicate lacing ranges.
    candidates.sort_by_key(|candidate| {
        (
            std::cmp::Reverse(candidate.level),
            candidate.reference_span_bp,
            candidate.reference_start_step,
            candidate.reference_end_step,
        )
    });

    let mut claimed_steps: Vec<BitVec> = graph
        .paths
        .iter()
        .map(|path| bitvec![0; path.steps.len()])
        .collect();
    let mut blocks = Vec::new();
    let mut selected_flubble_blocks = 0usize;
    let mut overlap_skipped_sites = 0usize;

    for candidate in candidates {
        if ranges_overlap_claimed(&candidate.path_ranges, &claimed_steps) {
            overlap_skipped_sites += 1;
            continue;
        }
        mark_claimed_ranges(&candidate.path_ranges, &mut claimed_steps);
        let candidate_blocks =
            topological_split_with_source(candidate.path_ranges, graph, candidate.source);
        selected_flubble_blocks += candidate_blocks.len();
        blocks.extend(candidate_blocks);
    }

    let gap_ranges = unclaimed_path_ranges(graph, &claimed_steps);
    let mut gap_blocks = topological_split_with_source(gap_ranges, graph, BlockSource::Gap);
    let gap_block_count = gap_blocks.len();
    blocks.append(&mut gap_blocks);

    Ok(FlubbleGuidedBlockPlan {
        blocks,
        povu_sites: decomposition.sites.len(),
        level_ge_1_sites,
        selected_flubble_blocks,
        overlap_skipped_sites,
        gap_blocks: gap_block_count,
    })
}

fn neighbor_merge_blocks(
    graph: &SmoothGraph,
    path_positions: &[Vec<usize>],
    config: &SmoothConfig,
    max_block_bp: usize,
) -> io::Result<NeighborMergeBlockPlan> {
    if graph.paths.is_empty() || max_block_bp == 0 {
        return Ok(NeighborMergeBlockPlan {
            blocks: Vec::new(),
            povu_sites: 0,
            candidate_sites: 0,
            selected_sites: 0,
            selected_groups: 0,
            neighbor_blocks: 0,
            blocks_merged: 0,
            sites_in_merged_blocks: 0,
            overlap_skipped_sites: 0,
            over_cap_skipped_sites: 0,
            gap_blocks: 0,
        });
    }

    let rendered = render_smooth_graph_gfa(graph)?;
    let native = povu::NativeGfa::parse(&rendered).map_err(povu_to_io_error)?;
    let reference_names = if config.flubble_reference_names.is_empty() {
        vec![graph.paths[0].name.clone()]
    } else {
        config.flubble_reference_names.clone()
    };
    let decomposition = match native.decompose_flubbles(&reference_names) {
        Ok(decomposition) => decomposition,
        Err(err) if !config.flubble_reference_names.is_empty() => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "POVU flubble reference hint(s) {:?} did not match graph paths: {}",
                    config.flubble_reference_names, err
                ),
            ));
        }
        Err(err) => return Err(povu_to_io_error(err)),
    };
    let root_positions = path_positions
        .get(decomposition.reference_path_index)
        .ok_or_else(|| {
            io::Error::other("POVU reference path index is outside smooth graph paths")
        })?;

    let id_to_idx = (0..graph.nodes.len())
        .map(|idx| ((idx + 1).to_string(), idx))
        .collect::<FxHashMap<_, _>>();
    let path_step_indexes = graph
        .paths
        .iter()
        .map(smooth_path_step_index)
        .collect::<Vec<_>>();

    let mut candidates = Vec::new();
    for site in &decomposition.sites {
        let Some(entry) = smooth_step_from_povu(&id_to_idx, &site.start) else {
            continue;
        };
        let Some(exit) = smooth_step_from_povu(&id_to_idx, &site.end) else {
            continue;
        };
        let reference_end_idx = site.reference_end_step.saturating_add(1);
        let reference_span_bp = if reference_end_idx < root_positions.len()
            && site.reference_start_step < root_positions.len()
        {
            root_positions[reference_end_idx] - root_positions[site.reference_start_step]
        } else {
            continue;
        };

        let mut path_ranges = Vec::new();
        for (path_idx, path_index) in path_step_indexes.iter().enumerate() {
            if let Some((begin_step, end_step)) =
                unique_smooth_anchor_range(path_index, entry, exit)
            {
                let length =
                    path_positions[path_idx][end_step] - path_positions[path_idx][begin_step];
                if length > 0 {
                    path_ranges.push(PathRange {
                        path_idx,
                        begin_step,
                        end_step,
                        length,
                    });
                }
            }
        }
        if path_ranges.len() < 2 {
            continue;
        }
        path_ranges.sort_by(|a, b| b.length.cmp(&a.length));
        candidates.push(NeighborSiteCandidate {
            entry,
            exit,
            level: site.level,
            reference_start_step: site.reference_start_step,
            reference_end_step: site.reference_end_step,
            reference_span_bp,
        });
    }
    let candidate_sites = candidates.len();

    // Pick the smallest/deepest non-overlapping sites first, independent of
    // parent ID. That keeps natural local bubbles while still allowing sites
    // from unrelated tree positions to become neighbors in the later path-order
    // grouping step.
    candidates.sort_by_key(|candidate| {
        (
            candidate.reference_span_bp,
            std::cmp::Reverse(candidate.level),
            candidate.reference_start_step,
            candidate.reference_end_step,
        )
    });

    let mut occupied_reference: Vec<(usize, usize)> = Vec::new();
    let mut selected = Vec::new();
    let mut overlap_skipped_sites = 0usize;
    let mut over_cap_skipped_sites = 0usize;
    for candidate in candidates {
        if candidate.reference_span_bp > max_block_bp {
            over_cap_skipped_sites += 1;
            continue;
        }
        let end_exclusive = candidate
            .reference_end_step
            .max(candidate.reference_start_step.saturating_add(1));
        if step_interval_conflicts(
            &occupied_reference,
            candidate.reference_start_step,
            end_exclusive,
        ) {
            overlap_skipped_sites += 1;
            continue;
        }
        insert_step_interval(
            &mut occupied_reference,
            candidate.reference_start_step,
            end_exclusive,
        );
        selected.push(candidate);
    }

    selected.sort_by_key(|candidate| {
        (
            candidate.reference_start_step,
            candidate.reference_end_step,
            candidate.reference_span_bp,
        )
    });
    let selected_sites = selected.len();
    let groups = greedy_neighbor_site_groups(selected, root_positions, max_block_bp);
    let selected_groups = groups.len();
    let blocks_merged = groups.iter().filter(|group| group.sites.len() > 1).count();
    let sites_in_merged_blocks = groups
        .iter()
        .filter(|group| group.sites.len() > 1)
        .map(|group| group.sites.len())
        .sum();

    let mut claimed_steps: Vec<BitVec> = graph
        .paths
        .iter()
        .map(|path| bitvec![0; path.steps.len()])
        .collect();
    let mut blocks = Vec::new();
    let mut final_overlap_skipped_sites = 0usize;
    for group in groups {
        let Some(first) = group.sites.first() else {
            continue;
        };
        let Some(last) = group.sites.last() else {
            continue;
        };
        let mut path_ranges = Vec::new();
        for (path_idx, path_index) in path_step_indexes.iter().enumerate() {
            if let Some((begin_step, end_step)) =
                unique_smooth_anchor_range(path_index, first.entry, last.exit)
            {
                let length =
                    path_positions[path_idx][end_step] - path_positions[path_idx][begin_step];
                if length > 0 {
                    path_ranges.push(PathRange {
                        path_idx,
                        begin_step,
                        end_step,
                        length,
                    });
                }
            }
        }
        if path_ranges.len() < 2 || ranges_overlap_claimed(&path_ranges, &claimed_steps) {
            final_overlap_skipped_sites += group.sites.len();
            continue;
        }
        path_ranges.sort_by(|a, b| b.length.cmp(&a.length));
        mark_claimed_ranges(&path_ranges, &mut claimed_steps);
        let source = BlockSource::NeighborMerge {
            site_count: group.sites.len(),
            reference_start_step: group.reference_start_step,
            reference_end_step: group.reference_end_step,
            reference_span_bp: group.reference_span_bp,
        };
        blocks.extend(topological_split_with_source(path_ranges, graph, source));
    }
    overlap_skipped_sites += final_overlap_skipped_sites;
    let neighbor_blocks = blocks.len();

    let gap_ranges = unclaimed_path_ranges(graph, &claimed_steps);
    let mut gap_blocks = topological_split_with_source(gap_ranges, graph, BlockSource::Gap);
    let gap_block_count = gap_blocks.len();
    blocks.append(&mut gap_blocks);

    Ok(NeighborMergeBlockPlan {
        blocks,
        povu_sites: decomposition.sites.len(),
        candidate_sites,
        selected_sites,
        selected_groups,
        neighbor_blocks,
        blocks_merged,
        sites_in_merged_blocks,
        overlap_skipped_sites,
        over_cap_skipped_sites,
        gap_blocks: gap_block_count,
    })
}

fn greedy_neighbor_site_groups(
    selected: Vec<NeighborSiteCandidate>,
    root_positions: &[usize],
    max_block_bp: usize,
) -> Vec<NeighborSiteGroup> {
    let mut groups = Vec::new();
    let mut current: Vec<NeighborSiteCandidate> = Vec::new();
    let mut current_start_step = 0usize;
    let mut current_end_step = 0usize;
    let mut current_span_bp = 0usize;

    for candidate in selected {
        if current.is_empty() {
            current_start_step = candidate.reference_start_step;
            current_end_step = candidate.reference_end_step;
            current_span_bp = candidate.reference_span_bp;
            current.push(candidate);
            continue;
        }

        let proposed_end_step = current_end_step.max(candidate.reference_end_step);
        let proposed_span_bp = if proposed_end_step + 1 < root_positions.len()
            && current_start_step < root_positions.len()
        {
            root_positions[proposed_end_step + 1] - root_positions[current_start_step]
        } else {
            usize::MAX
        };
        let path_adjacent_or_downstream = candidate.reference_start_step >= current_end_step;
        if path_adjacent_or_downstream && proposed_span_bp <= max_block_bp {
            current_end_step = proposed_end_step;
            current_span_bp = proposed_span_bp;
            current.push(candidate);
            continue;
        }

        groups.push(NeighborSiteGroup {
            sites: std::mem::take(&mut current),
            reference_start_step: current_start_step,
            reference_end_step: current_end_step,
            reference_span_bp: current_span_bp,
        });
        current_start_step = candidate.reference_start_step;
        current_end_step = candidate.reference_end_step;
        current_span_bp = candidate.reference_span_bp;
        current.push(candidate);
    }

    if !current.is_empty() {
        groups.push(NeighborSiteGroup {
            sites: current,
            reference_start_step: current_start_step,
            reference_end_step: current_end_step,
            reference_span_bp: current_span_bp,
        });
    }
    groups
}

fn step_interval_conflicts(occupied: &[(usize, usize)], begin: usize, end: usize) -> bool {
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

fn insert_step_interval(occupied: &mut Vec<(usize, usize)>, begin: usize, end: usize) {
    let insert_pos = occupied.partition_point(|&(occupied_begin, _)| occupied_begin < begin);
    occupied.insert(insert_pos, (begin, end));
}

fn render_smooth_graph_gfa(graph: &SmoothGraph) -> io::Result<String> {
    let mut out = String::new();
    out.push_str("H\tVN:Z:1.0\n");
    for (idx, seq) in graph.nodes.iter().enumerate() {
        let seq = std::str::from_utf8(seq)
            .map_err(|err| io::Error::other(format!("non-UTF8 node sequence: {}", err)))?;
        out.push_str(&format!("S\t{}\t{}\n", idx + 1, seq));
    }

    let mut edges = BTreeSet::new();
    for &(from_idx, from_rev, to_idx, to_rev) in &graph.edges {
        edges.insert((from_idx, from_rev, to_idx, to_rev));
    }
    for path in &graph.paths {
        for pair in path.steps.windows(2) {
            edges.insert((pair[0].0, pair[0].1, pair[1].0, pair[1].1));
        }
    }
    for (from_idx, from_rev, to_idx, to_rev) in edges {
        out.push_str(&format!(
            "L\t{}\t{}\t{}\t{}\t0M\n",
            from_idx + 1,
            if from_rev { "-" } else { "+" },
            to_idx + 1,
            if to_rev { "-" } else { "+" },
        ));
    }

    for path in &graph.paths {
        let steps = path
            .steps
            .iter()
            .map(|&(node_idx, is_reverse)| {
                format!("{}{}", node_idx + 1, if is_reverse { "-" } else { "+" })
            })
            .collect::<Vec<_>>()
            .join(",");
        if !steps.is_empty() {
            out.push_str(&format!("P\t{}\t{}\t*\n", path.name, steps));
        }
    }
    Ok(out)
}

fn povu_to_io_error(err: povu::Error) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("POVU flubble decomposition failed: {err}"),
    )
}

fn smooth_path_step_index(path: &PathData) -> FxHashMap<(usize, bool), Vec<usize>> {
    let mut index: FxHashMap<(usize, bool), Vec<usize>> = FxHashMap::default();
    for (idx, &step) in path.steps.iter().enumerate() {
        index.entry(step).or_default().push(idx);
    }
    index
}

fn smooth_step_from_povu(
    id_to_idx: &FxHashMap<String, usize>,
    step: &povu::native_gfa::Step,
) -> Option<(usize, bool)> {
    let node_idx = *id_to_idx.get(&step.segment)?;
    Some((
        node_idx,
        matches!(step.strand, povu::native_gfa::Strand::Reverse),
    ))
}

fn unique_smooth_anchor_range(
    path_index: &FxHashMap<(usize, bool), Vec<usize>>,
    entry: (usize, bool),
    exit: (usize, bool),
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

fn ranges_overlap_claimed(path_ranges: &[PathRange], claimed_steps: &[BitVec]) -> bool {
    path_ranges.iter().any(|range| {
        (range.begin_step..range.end_step).any(|step_idx| claimed_steps[range.path_idx][step_idx])
    })
}

fn mark_claimed_ranges(path_ranges: &[PathRange], claimed_steps: &mut [BitVec]) {
    for range in path_ranges {
        for step_idx in range.begin_step..range.end_step {
            claimed_steps[range.path_idx].set(step_idx, true);
        }
    }
}

fn unclaimed_path_ranges(graph: &SmoothGraph, claimed_steps: &[BitVec]) -> Vec<PathRange> {
    let mut ranges = Vec::new();
    for (path_idx, path) in graph.paths.iter().enumerate() {
        let mut begin = None;
        for step_idx in 0..path.steps.len() {
            if !claimed_steps[path_idx][step_idx] {
                begin.get_or_insert(step_idx);
                continue;
            }
            if let Some(start) = begin.take() {
                let length = compute_range_length(graph, path_idx, start, step_idx);
                if length > 0 {
                    ranges.push(PathRange {
                        path_idx,
                        begin_step: start,
                        end_step: step_idx,
                        length,
                    });
                }
            }
        }
        if let Some(start) = begin.take() {
            let length = compute_range_length(graph, path_idx, start, path.steps.len());
            if length > 0 {
                ranges.push(PathRange {
                    path_idx,
                    begin_step: start,
                    end_step: path.steps.len(),
                    length,
                });
            }
        }
    }
    ranges.sort_by(|a, b| b.length.cmp(&a.length));
    ranges
}

// ---------------------------------------------------------------------------
// Block breaking (smoothxg break_blocks)
// ---------------------------------------------------------------------------

fn break_blocks(blocks: &mut Vec<Block>, graph: &SmoothGraph, max_poa_length: usize) {
    let mut new_blocks: Vec<Block> = Vec::with_capacity(blocks.len());

    for block in blocks.drain(..) {
        if block.source.is_gap() {
            new_blocks.push(block);
            continue;
        }

        // Check if any path_range exceeds max_poa_length
        let needs_break = block.path_ranges.iter().any(|r| r.length > max_poa_length);

        if !needs_break || block.path_ranges.len() <= 1 {
            new_blocks.push(block);
            continue;
        }

        // Try repeat detection on long sequences
        let mut cut_length = max_poa_length;
        for range in &block.path_ranges {
            if range.length < 2000 {
                // min_copy_length * 2
                continue;
            }
            let seq = extract_path_range_sequence(graph, range);
            if let Some(repeat_len) = detect_repeat(&seq, 1000, 20000, 5.0, 50) {
                cut_length = (repeat_len / 2).max(1);
                break;
            }
        }

        // Cut path ranges at cut_length intervals
        let mut new_ranges: Vec<PathRange> = Vec::new();
        for range in &block.path_ranges {
            if range.length <= cut_length {
                new_ranges.push(range.clone());
            } else {
                let mut current_start = range.begin_step;
                let mut current_length = 0;

                for step_idx in range.begin_step..range.end_step {
                    let (node_idx, _) = graph.paths[range.path_idx].steps[step_idx];
                    current_length += graph.nodes[node_idx].len();

                    if current_length >= cut_length && step_idx + 1 < range.end_step {
                        new_ranges.push(PathRange {
                            path_idx: range.path_idx,
                            begin_step: current_start,
                            end_step: step_idx + 1,
                            length: current_length,
                        });
                        current_start = step_idx + 1;
                        current_length = 0;
                    }
                }

                // Last piece
                if current_start < range.end_step && current_length > 0 {
                    new_ranges.push(PathRange {
                        path_idx: range.path_idx,
                        begin_step: current_start,
                        end_step: range.end_step,
                        length: current_length,
                    });
                }
            }
        }

        // Re-sort by length (longest first)
        new_ranges.sort_by(|a, b| b.length.cmp(&a.length));
        new_blocks.push(Block {
            path_ranges: new_ranges,
            source: block.source,
        });
    }

    *blocks = new_blocks;
}

/// Extract the nucleotide sequence for a path range.
fn extract_path_range_sequence(graph: &SmoothGraph, range: &PathRange) -> Vec<u8> {
    let path = &graph.paths[range.path_idx];
    let mut seq = Vec::with_capacity(range.length);
    for step_idx in range.begin_step..range.end_step {
        let (node_idx, is_reverse) = path.steps[step_idx];
        let node_seq = &graph.nodes[node_idx];
        if is_reverse {
            seq.extend(reverse_complement(node_seq));
        } else {
            seq.extend_from_slice(node_seq);
        }
    }
    seq
}

// ---------------------------------------------------------------------------
// Autocorrelation-based repeat detection (port of sautocorr)
// ---------------------------------------------------------------------------

/// Detect repeat period in a sequence via autocorrelation.
/// Returns the repeat period (lag) if a significant repeat is found, None otherwise.
fn detect_repeat(
    sequence: &[u8],
    min_copy_length: usize,
    max_copy_length: usize,
    min_z: f64,
    stride: usize,
) -> Option<usize> {
    if sequence.len() < 2 * min_copy_length {
        return None;
    }

    // Convert to numeric values
    let values: Vec<f64> = sequence
        .iter()
        .map(|&b| match b {
            b'A' | b'a' => 0.0,
            b'C' | b'c' => 1.0,
            b'G' | b'g' => 2.0,
            b'T' | b't' => 3.0,
            _ => 2.0,
        })
        .collect();

    let n = values.len();

    // Compute mean and variance
    let mean = values.iter().sum::<f64>() / n as f64;
    let variance = values.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / n as f64;

    if variance < 1e-10 {
        return None;
    }

    let min_lag = min_copy_length;
    let max_lag = max_copy_length.min(n / 2);

    if min_lag >= max_lag {
        return None;
    }

    // Compute autocorrelation at each lag with stride sampling
    let stride = stride.max(1);
    let mut autocorrs: Vec<f64> = Vec::with_capacity(max_lag - min_lag);
    for lag in min_lag..max_lag {
        let mut sum = 0.0;
        let mut count = 0usize;
        let mut i = 0;
        while i + lag < n {
            sum += (values[i] - mean) * (values[i + lag] - mean);
            count += 1;
            i += stride;
        }
        if count > 0 {
            autocorrs.push(sum / (count as f64 * variance));
        } else {
            autocorrs.push(0.0);
        }
    }

    if autocorrs.is_empty() {
        return None;
    }

    // Normalize to z-scores
    let acf_mean = autocorrs.iter().sum::<f64>() / autocorrs.len() as f64;
    let acf_var = autocorrs
        .iter()
        .map(|&x| (x - acf_mean) * (x - acf_mean))
        .sum::<f64>()
        / autocorrs.len() as f64;
    let acf_std = acf_var.sqrt();

    if acf_std < 1e-10 {
        return None;
    }

    // Find peak z-score above threshold
    let mut best_lag = 0;
    let mut best_z = 0.0;
    let mut found = false;

    for (i, &acf) in autocorrs.iter().enumerate() {
        let z = (acf - acf_mean) / acf_std;
        if z > min_z {
            if z > best_z {
                best_z = z;
                best_lag = min_lag + i;
            }
            found = true;
        } else if found {
            break;
        }
    }

    if found {
        Some(best_lag)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Per-block SPOA smoothing
// ---------------------------------------------------------------------------

/// Emit a block's path-ranges as a sub-GFA using the **un-smoothed** input
/// sequence content. Called when [`smooth_block`] returns an empty/Err result
/// so the block's bp coverage stays in the laced graph rather than vanishing.
///
/// Node IDs are local to the emitted sub-GFA (1..=N over the unique input
/// nodes touched). The lacer assigns globally fresh IDs anyway, so local
/// numbering doesn't collide. Edges include only step-to-step transitions
/// taken inside the block's path ranges — the parent graph's full edge set
/// would carry unrelated traversals through these nodes.
fn passthrough_block_gfa(
    block: &Block,
    graph: &SmoothGraph,
    path_positions: &[Vec<usize>],
) -> io::Result<String> {
    if block.path_ranges.is_empty() {
        return Ok(String::new());
    }

    let mut node_remap: FxHashMap<usize, u64> = FxHashMap::default();
    let mut next_id: u64 = 1;
    let mut out = String::new();

    for range in &block.path_ranges {
        let path = &graph.paths[range.path_idx];
        for step_idx in range.begin_step..range.end_step {
            let (node_idx, _) = path.steps[step_idx];
            if let std::collections::hash_map::Entry::Vacant(e) = node_remap.entry(node_idx) {
                e.insert(next_id);
                let seq = std::str::from_utf8(&graph.nodes[node_idx])
                    .map_err(|err| io::Error::other(format!("non-UTF8 node sequence: {}", err)))?;
                out.push_str(&format!("S\t{}\t{}\n", next_id, seq));
                next_id += 1;
            }
        }
    }

    let mut edges: FxHashSet<(u64, bool, u64, bool)> = FxHashSet::default();
    for range in &block.path_ranges {
        let path = &graph.paths[range.path_idx];
        for window_start in range.begin_step..range.end_step.saturating_sub(1) {
            let (from_idx, from_rev) = path.steps[window_start];
            let (to_idx, to_rev) = path.steps[window_start + 1];
            if let (Some(&from_id), Some(&to_id)) =
                (node_remap.get(&from_idx), node_remap.get(&to_idx))
            {
                edges.insert((from_id, from_rev, to_id, to_rev));
            }
        }
    }
    for (from_id, from_rev, to_id, to_rev) in &edges {
        out.push_str(&format!(
            "L\t{}\t{}\t{}\t{}\t0M\n",
            from_id,
            if *from_rev { "-" } else { "+" },
            to_id,
            if *to_rev { "-" } else { "+" },
        ));
    }

    for range in &block.path_ranges {
        let path = &graph.paths[range.path_idx];
        let (path_key, path_start) = parse_path_coords(&path.name);
        let range_start_bp = path_start + path_positions[range.path_idx][range.begin_step];
        let range_end_bp = path_start + path_positions[range.path_idx][range.end_step];
        let mut steps_str = String::new();
        for step_idx in range.begin_step..range.end_step {
            let (node_idx, is_rev) = path.steps[step_idx];
            let nid = node_remap[&node_idx];
            if !steps_str.is_empty() {
                steps_str.push(',');
            }
            steps_str.push_str(&format!("{}{}", nid, if is_rev { "-" } else { "+" }));
        }
        out.push_str(&format!(
            "P\t{}:{}-{}\t{}\t*\n",
            path_key, range_start_bp, range_end_bp, steps_str,
        ));
    }

    Ok(out)
}

fn block_input_bp(block: &Block) -> usize {
    block.path_ranges.iter().map(|range| range.length).sum()
}

fn gfa_segment_bp(gfa: &str) -> usize {
    gfa.lines()
        .filter_map(|line| {
            let mut fields = line.split('\t');
            if fields.next()? != "S" {
                return None;
            }
            Some(fields.nth(1)?.len())
        })
        .sum()
}

fn log_flubble_block_compression(
    block_idx: usize,
    block: &Block,
    input_bp: usize,
    output_bp: usize,
    passthrough: bool,
) {
    if let BlockSource::Flubble {
        site_id,
        level,
        reference_start_step,
        reference_end_step,
        is_leaf,
    } = &block.source
    {
        let ratio = if input_bp > 0 {
            output_bp as f64 / input_bp as f64
        } else {
            0.0
        };
        info!(
            "[smooth] flubble-block {} site={} level={} leaf={} ref_steps={}-{} ranges={} input_bp={} output_bp={} ratio={:.4} mode={}",
            block_idx,
            site_id,
            level,
            is_leaf,
            reference_start_step,
            reference_end_step,
            block.path_ranges.len(),
            input_bp,
            output_bp,
            ratio,
            if passthrough { "passthrough" } else { "poa" },
        );
    } else if let BlockSource::NeighborMerge {
        site_count,
        reference_start_step,
        reference_end_step,
        reference_span_bp,
    } = &block.source
    {
        let ratio = if input_bp > 0 {
            output_bp as f64 / input_bp as f64
        } else {
            0.0
        };
        info!(
            "[neighbor-merge-poasta] block {} sites={} ref_steps={}-{} ref_span_bp={} ranges={} input_bp={} output_bp={} ratio={:.4} mode={}",
            block_idx,
            site_count,
            reference_start_step,
            reference_end_step,
            reference_span_bp,
            block.path_ranges.len(),
            input_bp,
            output_bp,
            ratio,
            if passthrough { "passthrough" } else { "poasta" },
        );
    }
}

fn smooth_empty_block_error(block_idx: usize, block: &Block, input_bp: usize) -> String {
    let source = match &block.source {
        BlockSource::PathOverlap => "path-overlap",
        BlockSource::Flubble { .. } => "flubble",
        BlockSource::NeighborMerge { .. } => "neighbor-merge",
        BlockSource::Gap => "gap",
    };
    let max_range_bp = block
        .path_ranges
        .iter()
        .map(|range| range.length)
        .max()
        .unwrap_or(0);
    format!(
        "smooth block {block_idx} ({source}, ranges={}, input_bp={}, max_range_bp={}) produced no path-preserving GFA after padded/no-padding POA attempts; refusing graph-loss output",
        block.path_ranges.len(),
        input_bp,
        max_range_bp
    )
}

/// Realign a neighbor-merged local block with POASTA and emit lacer-compatible
/// GFA. The surrounding smooth pass owns block placement, path coordinate
/// naming, gap passthrough, smoothable-block error handling, and final lacing.
fn poasta_block(
    block: &Block,
    graph: &SmoothGraph,
    config: &SmoothConfig,
    path_positions: &[Vec<usize>],
) -> io::Result<String> {
    if block.path_ranges.is_empty() {
        return Ok(String::new());
    }

    let mut headers = Vec::with_capacity(block.path_ranges.len());
    let mut sequences = Vec::with_capacity(block.path_ranges.len());
    for range in &block.path_ranges {
        let path = &graph.paths[range.path_idx];
        let (path_key, path_start) = parse_path_coords(&path.name);
        let range_start_bp = path_start + path_positions[range.path_idx][range.begin_step];
        let range_end_bp = path_start + path_positions[range.path_idx][range.end_step];
        headers.push(format!("{}:{}-{}", path_key, range_start_bp, range_end_bp));
        sequences.push(extract_path_range_sequence(graph, range));
    }

    crate::resolution::poasta_sequences_to_gfa(&headers, &sequences, config.scoring_params)
}

fn smooth_entries_without_padding(
    entries: &[SeqEntry],
    config: &SmoothConfig,
) -> io::Result<String> {
    if entries.is_empty() {
        return Ok(String::new());
    }

    let mut ordered_entries: Vec<&SeqEntry> = entries.iter().collect();
    ordered_entries.sort_by(|a, b| b.core_sequence.len().cmp(&a.core_sequence.len()));

    let (mut spoa_graph, mut spoa_engine) = build_spoa_engine(config.scoring_params);
    feed_sequences_to_graph(
        &mut spoa_engine,
        &mut spoa_graph,
        ordered_entries
            .iter()
            .map(|entry| entry.core_sequence.as_str()),
    );
    let headers: Vec<String> = ordered_entries
        .iter()
        .map(|entry| entry.path_name.clone())
        .collect();
    let gfa = spoa_graph.generate_gfa(&headers, false);
    unchop_gfa(&gfa)
}

/// Smooth a single block: extract sequences → SPOA → GFA.
fn smooth_block(
    block: &Block,
    graph: &SmoothGraph,
    config: &SmoothConfig,
    path_positions: &[Vec<usize>],
) -> io::Result<String> {
    if block.path_ranges.is_empty() {
        return Ok(String::new());
    }

    // 1. Compute padding
    let avg_seq_len = block
        .path_ranges
        .iter()
        .map(|r| r.length as f64)
        .sum::<f64>()
        / block.path_ranges.len() as f64;

    let poa_padding: usize = if config.poa_padding_fraction > 0.0 {
        let base = if block.path_ranges.len() <= config.max_block_depth_for_padding_more {
            311
        } else {
            0
        };
        base.max((avg_seq_len * config.poa_padding_fraction) as usize)
    } else {
        0
    };

    // 2. Extract sequences with padding
    let mut entries: Vec<SeqEntry> = Vec::with_capacity(block.path_ranges.len());

    for range in &block.path_ranges {
        let path = &graph.paths[range.path_idx];
        let (path_key, path_start) = parse_path_coords(&path.name);

        // Extract core sequence
        let mut core_seq: Vec<u8> = Vec::with_capacity(range.length);
        for step_idx in range.begin_step..range.end_step {
            let (node_idx, is_reverse) = path.steps[step_idx];
            let node_seq = &graph.nodes[node_idx];
            if is_reverse {
                core_seq.extend(reverse_complement(node_seq));
            } else {
                core_seq.extend_from_slice(node_seq);
            }
        }

        // Extract left padding (walk backwards from begin_step)
        let mut left_pad_seq = Vec::new();
        if poa_padding > 0 {
            let mut remaining = poa_padding;
            let mut si = range.begin_step;
            while si > 0 && remaining > 0 {
                si -= 1;
                let (node_idx, is_reverse) = path.steps[si];
                let node_seq = &graph.nodes[node_idx];
                let take = node_seq.len().min(remaining);
                if is_reverse {
                    let rc = reverse_complement(node_seq);
                    // Take the LAST `take` bases (we're walking backwards)
                    left_pad_seq.extend_from_slice(&rc[rc.len() - take..]);
                } else {
                    left_pad_seq.extend_from_slice(&node_seq[node_seq.len() - take..]);
                }
                remaining -= take;
            }
            left_pad_seq.reverse();
        }

        // Extract right padding (walk forwards from end_step)
        let mut right_pad_seq = Vec::new();
        if poa_padding > 0 {
            let mut remaining = poa_padding;
            let mut si = range.end_step;
            while si < path.steps.len() && remaining > 0 {
                let (node_idx, is_reverse) = path.steps[si];
                let node_seq = &graph.nodes[node_idx];
                let take = node_seq.len().min(remaining);
                if is_reverse {
                    let rc = reverse_complement(node_seq);
                    right_pad_seq.extend_from_slice(&rc[..take]);
                } else {
                    right_pad_seq.extend_from_slice(&node_seq[..take]);
                }
                remaining -= take;
                si += 1;
            }
        }

        let actual_left_pad = left_pad_seq.len();
        let actual_right_pad = right_pad_seq.len();

        // Combine: left_pad + core + right_pad
        let mut full_seq = Vec::with_capacity(actual_left_pad + core_seq.len() + actual_right_pad);
        full_seq.extend_from_slice(&left_pad_seq);
        full_seq.extend_from_slice(&core_seq);
        full_seq.extend_from_slice(&right_pad_seq);

        // Compute path coordinates for naming
        let range_start_bp = path_start + path_positions[range.path_idx][range.begin_step];
        let range_end_bp = path_start + path_positions[range.path_idx][range.end_step];
        let path_name = format!("{}:{}-{}", path_key, range_start_bp, range_end_bp);

        entries.push(SeqEntry {
            sequence: String::from_utf8_lossy(&full_seq).to_string(),
            core_sequence: String::from_utf8_lossy(&core_seq).to_string(),
            path_name,
            left_pad: actual_left_pad,
            right_pad: actual_right_pad,
        });
    }

    if entries.is_empty() {
        return Ok(String::new());
    }

    // Sort by sequence length descending (SPOA benefits from longest-first)
    entries.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));

    let any_padding = entries.iter().any(|e| e.left_pad > 0 || e.right_pad > 0);
    if !any_padding {
        // No padding — generate GFA directly
        return smooth_entries_without_padding(&entries, config);
    }

    // 3. Run SPOA on padded sequences to find the core alignment window.
    let (mut spoa_graph, mut spoa_engine) = build_spoa_engine(config.scoring_params);
    feed_sequences_to_graph(
        &mut spoa_engine,
        &mut spoa_graph,
        entries.iter().map(|e| e.sequence.as_str()),
    );

    // 4. Handle padding trimming

    // With padding: use MSA-based trimming approach
    let msa = spoa_graph.generate_msa();

    // Find core column range (where all sequences have consumed their padding)
    let (core_start, core_end) = find_core_column_range(&msa, &entries);

    if core_start >= core_end {
        return smooth_entries_without_padding(&entries, config);
    }

    // Extract core subsequences and re-run SPOA
    let nseqs = msa.len().min(entries.len());
    let msa_bytes: Vec<&[u8]> = msa.iter().map(|s| s.as_bytes()).collect();

    let mut core_sequences: Vec<String> = Vec::with_capacity(nseqs);
    let mut n_empty_core = 0usize;
    for row in &msa_bytes[..nseqs] {
        let core_seq: String = row[core_start..core_end]
            .iter()
            .filter(|&&b| b != b'-')
            .map(|&b| b as char)
            .collect();
        if core_seq.is_empty() {
            n_empty_core += 1;
        }
        core_sequences.push(core_seq);
    }

    // Any entry with an empty core would be fed an empty string to SPOA, which
    // silently emits no P-line for that entry and loses that input path-range's
    // bp coverage. Even ONE empty core is a correctness problem, not just a
    // quality one: the laced graph's path-bp coverage contract is "every input
    // range survives somewhere". Retry the block as unpadded POA rather than
    // emit partial smoothed output or substitute the original subgraph.
    if n_empty_core > 0 {
        return smooth_entries_without_padding(&entries, config);
    }

    // Stronger check: `find_core_column_range`'s MAX-of-left-pad-end /
    // MIN-of-right-pad-start can collapse the core to a tiny intersection
    // when sequences have very different padding/length ratios — even when
    // no single core is empty. SPOA then emits P-lines whose name claims
    // the full input bp range (set above from path_positions) but whose
    // step walk covers only the tiny core. The lacer trusts the name, the
    // next pass parses the GFA and sees a shorter path: observed on C4
    // `impg graph`, block #160 named 60,802 bp but walked 1,962 bp (97%
    // silent shrinkage). The contract — name's bp range == walked content
    // bp — must hold. If the sum of trimmed core bp is far less than the
    // sum of input range bp, route the block through passthrough.
    let total_core_bp: usize = core_sequences.iter().map(|s| s.len()).sum();
    let total_input_bp: usize = entries
        .iter()
        .take(nseqs)
        .map(|e| e.sequence.len() - e.left_pad - e.right_pad)
        .sum();
    // Strict equality. Any bp loss from MSA-trim ends up as a P-line whose name
    // claims input bp coords but whose step walk is shorter. The lacer trusts
    // the name and so does the next pass — the bp shortfall becomes a real
    // coordinate error in the output graph. Use no-padding POA for this block
    // instead; that preserves the path range contract without passthrough.
    if total_core_bp != total_input_bp {
        return smooth_entries_without_padding(&entries, config);
    }

    // Build fresh SPOA graph from core sequences
    let (mut graph2, mut engine2) = build_spoa_engine(config.scoring_params);
    feed_sequences_to_graph(
        &mut engine2,
        &mut graph2,
        core_sequences.iter().map(|s| s.as_str()),
    );

    let headers: Vec<String> = entries
        .iter()
        .take(nseqs)
        .map(|e| e.path_name.clone())
        .collect();
    let gfa = graph2.generate_gfa(&headers, false);
    unchop_gfa(&gfa)
}

/// Find the MSA column range corresponding to the core (non-padding) region.
#[allow(clippy::needless_range_loop)]
fn find_core_column_range(msa: &[String], entries: &[SeqEntry]) -> (usize, usize) {
    let ncols = msa.first().map(|s| s.len()).unwrap_or(0);
    if ncols == 0 {
        return (0, 0);
    }

    let msa_bytes: Vec<&[u8]> = msa.iter().map(|s| s.as_bytes()).collect();
    let nseqs = msa_bytes.len().min(entries.len());

    // Find core_start: max over all sequences of the column where left padding ends
    let mut core_start = 0;
    for i in 0..nseqs {
        let left_pad = entries[i].left_pad;
        if left_pad == 0 {
            continue;
        }
        let mut non_gap_count = 0;
        for (col, &b) in msa_bytes[i].iter().take(ncols).enumerate() {
            if b != b'-' {
                non_gap_count += 1;
                if non_gap_count == left_pad {
                    core_start = core_start.max(col + 1);
                    break;
                }
            }
        }
    }

    // Find core_end: min over all sequences of the column where right padding begins
    let mut core_end = ncols;
    for i in 0..nseqs {
        let right_pad = entries[i].right_pad;
        if right_pad == 0 {
            continue;
        }
        let total_bases: usize = msa_bytes[i].iter().filter(|&&b| b != b'-').count();
        let core_base_count = total_bases.saturating_sub(right_pad);
        if core_base_count == 0 {
            core_end = 0;
            break;
        }
        let mut non_gap_count = 0;
        for (col, &b) in msa_bytes[i].iter().take(ncols).enumerate() {
            if b != b'-' {
                non_gap_count += 1;
                if non_gap_count == core_base_count {
                    core_end = core_end.min(col + 1);
                    break;
                }
            }
        }
    }

    if core_start >= core_end {
        let mid = ncols / 2;
        return (mid, mid);
    }

    (core_start, core_end)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Parse path name into (key, start_offset).
/// If name has ":start-end" suffix, returns (key, start).
/// Otherwise, returns (name, 0).
fn parse_path_coords(name: &str) -> (String, usize) {
    if let Some(last_colon) = name.rfind(':') {
        let (key, range_str) = name.split_at(last_colon);
        let range_str = &range_str[1..];
        if let Some((start_str, _end_str)) = range_str.split_once('-') {
            if let Ok(start) = start_str.parse::<usize>() {
                return (key.to_string(), start);
            }
        }
    }
    (name.to_string(), 0)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_path_coords() {
        assert_eq!(parse_path_coords("chr1:100-200"), ("chr1".to_string(), 100));
        assert_eq!(
            parse_path_coords("sample#1#chr1:0-5000"),
            ("sample#1#chr1".to_string(), 0)
        );
        assert_eq!(parse_path_coords("chr1"), ("chr1".to_string(), 0));
        assert_eq!(
            parse_path_coords("s:123:456-789"),
            ("s:123".to_string(), 456)
        );
    }

    #[test]
    fn test_parse_gfa() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
S\t2\tTG
S\t3\tAATT
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tseq1:0-10\t1+,2+,3+\t*
P\tseq2:0-8\t1+,3+\t*
";
        let graph = parse_gfa(gfa);
        assert_eq!(graph.nodes.len(), 3);
        assert_eq!(graph.paths.len(), 2);
        assert_eq!(graph.paths[0].steps.len(), 3);
        assert_eq!(graph.paths[1].steps.len(), 2);
        assert_eq!(graph.edges.len(), 2);
    }

    #[test]
    fn test_chop_graph() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGTACGT
S\t2\tTG
P\tseq:0-10\t1+,2+\t*
";
        let mut graph = parse_gfa(gfa);
        assert_eq!(graph.nodes.len(), 2);

        // Chop to max 4bp
        chop_graph(&mut graph, 4);

        // Node 1 (8bp) should be split into 2 pieces of 4bp each
        // Node 2 (2bp) stays as-is
        assert_eq!(graph.nodes.len(), 3);
        assert_eq!(graph.nodes[0], b"ACGT");
        assert_eq!(graph.nodes[1], b"ACGT");
        assert_eq!(graph.nodes[2], b"TG");

        // Path should now have 3 steps
        assert_eq!(graph.paths[0].steps.len(), 3);
        assert_eq!(graph.paths[0].steps[0], (0, false));
        assert_eq!(graph.paths[0].steps[1], (1, false));
        assert_eq!(graph.paths[0].steps[2], (2, false));
    }

    #[test]
    fn test_compute_path_positions() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
S\t2\tTG
S\t3\tAATT
P\tseq:0-10\t1+,2+,3+\t*
";
        let graph = parse_gfa(gfa);
        let positions = compute_path_positions(&graph);
        assert_eq!(positions.len(), 1);
        assert_eq!(positions[0], vec![0, 4, 6, 10]);
    }

    #[test]
    fn test_union_find() {
        let mut uf = UnionFind::new(5);
        uf.union(0, 1);
        uf.union(2, 3);
        assert_eq!(uf.find(0), uf.find(1));
        assert_ne!(uf.find(0), uf.find(2));
        uf.union(1, 3);
        assert_eq!(uf.find(0), uf.find(3));
    }

    #[test]
    fn test_detect_repeat_no_repeat() {
        // Random sequence — should not detect repeats
        let seq = b"ACGTACGATCGATCGTAGCTAGCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        // Too short for min_copy_length=1000
        assert!(detect_repeat(seq, 1000, 20000, 5.0, 50).is_none());
    }

    #[test]
    fn test_extract_path_range_sequence() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
S\t2\tTG
S\t3\tAATT
P\tseq:0-10\t1+,2+,3+\t*
";
        let graph = parse_gfa(gfa);
        let range = PathRange {
            path_idx: 0,
            begin_step: 1,
            end_step: 3,
            length: 6,
        };
        let seq = extract_path_range_sequence(&graph, &range);
        assert_eq!(seq, b"TGAATT");
    }

    #[test]
    fn test_smoothable_blocks_simple() {
        // Small graph that should form a single block
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
S\t2\tTG
P\tseq1:0-6\t1+,2+\t*
P\tseq2:0-6\t1+,2+\t*
";
        let graph = parse_gfa(gfa);
        let index = build_node_step_index(&graph);
        let config = SmoothConfig::new(2);
        let target_poa_length = config.target_poa_lengths[0];
        let max_block_weight = target_poa_length * config.n_haps.max(1);
        let blocks = smoothable_blocks(&graph, &index, max_block_weight, target_poa_length);

        // Should produce 1 block (total weight is small)
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].path_ranges.len(), 2);
    }

    #[test]
    fn test_passthrough_block_gfa_preserves_path_coverage() {
        // Regression for the smooth fragmentation bug: when smooth_block
        // returns Ok("") (degenerate-padding or many-empty-cores case), the
        // dispatch in smooth_gfa_pass routes the block through
        // passthrough_block_gfa so its path-ranges survive to the lacer.
        //
        // Verifies: (1) one P-line per input path-range, (2) names use the
        // path's bp coords, (3) emitted node IDs are local (no clash with
        // other blocks during lacing), (4) edges link consecutive steps.
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
S\t2\tTGCA
S\t3\tGGGG
P\tchr1:100-108\t1+,2+\t*
P\tchr2:200-212\t1+,2+,3+\t*
";
        let graph = parse_gfa(gfa);
        let path_positions = compute_path_positions(&graph);
        let block = Block {
            path_ranges: vec![
                PathRange {
                    path_idx: 0,
                    begin_step: 0,
                    end_step: 2,
                    length: 8,
                },
                PathRange {
                    path_idx: 1,
                    begin_step: 0,
                    end_step: 3,
                    length: 12,
                },
            ],
            source: BlockSource::PathOverlap,
        };
        let out =
            passthrough_block_gfa(&block, &graph, &path_positions).expect("passthrough failed");
        let p_lines: Vec<&str> = out.lines().filter(|l| l.starts_with("P\t")).collect();
        assert_eq!(p_lines.len(), 2, "one P-line per path-range");
        // Names carry full input bp range, so the lacer sees consecutive ranges
        // and can join across blocks.
        let names: Vec<&str> = p_lines
            .iter()
            .map(|l| l.split('\t').nth(1).unwrap())
            .collect();
        assert!(
            names.iter().any(|n| *n == "chr1:100-108"),
            "got {:?}",
            names
        );
        assert!(
            names.iter().any(|n| *n == "chr2:200-212"),
            "got {:?}",
            names
        );
        // 3 distinct nodes touched → 3 S-lines with local IDs 1..=3.
        let s_lines: Vec<&str> = out.lines().filter(|l| l.starts_with("S\t")).collect();
        assert_eq!(s_lines.len(), 3);
        // At least two L-lines: edges 1→2 (both paths) and 2→3 (chr2 only).
        let l_lines: Vec<&str> = out.lines().filter(|l| l.starts_with("L\t")).collect();
        assert!(l_lines.len() >= 2, "edges: {:?}", l_lines);
    }

    #[test]
    fn test_smooth_block_retries_without_padding_when_trim_would_empty_core() {
        // The lower-locus PGGB failure produced an empty smoothable block
        // because MSA-based padding trim would have dropped a path's entire
        // core. That is not an aligner failure: the block can still be
        // validly POA-smoothed if retried without flank padding.
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tAAAA
S\t3\tAAAA
L\t1\t+\t2\t+\t0M
P\tp1:0-8\t1+,2+\t*
P\tp2:0-4\t3+\t*
";
        let graph = parse_gfa(gfa);
        let path_positions = compute_path_positions(&graph);
        let block = Block {
            path_ranges: vec![
                PathRange {
                    path_idx: 0,
                    begin_step: 1,
                    end_step: 2,
                    length: 4,
                },
                PathRange {
                    path_idx: 1,
                    begin_step: 0,
                    end_step: 1,
                    length: 4,
                },
            ],
            source: BlockSource::PathOverlap,
        };
        let config = SmoothConfig {
            poa_padding_fraction: 1.0,
            num_threads: 1,
            ..SmoothConfig::new(2)
        };

        let out = smooth_block(&block, &graph, &config, &path_positions)
            .expect("block should retry without padding instead of failing empty");
        assert!(!out.is_empty(), "no-padding retry should emit a GFA");
        let smoothed = parse_gfa(&out);
        let paths: FxHashMap<_, _> = smoothed
            .paths
            .iter()
            .map(|path| {
                let seq = path
                    .steps
                    .iter()
                    .flat_map(|&(node_idx, is_reverse)| {
                        if is_reverse {
                            reverse_complement(&smoothed.nodes[node_idx])
                        } else {
                            smoothed.nodes[node_idx].clone()
                        }
                    })
                    .collect::<Vec<_>>();
                (path.name.as_str(), String::from_utf8(seq).unwrap())
            })
            .collect();
        assert_eq!(paths.get("p1:4-8").map(String::as_str), Some("AAAA"));
        assert_eq!(paths.get("p2:0-4").map(String::as_str), Some("AAAA"));
    }

    #[test]
    fn test_flubble_guided_blocks_use_level_ge_1_sites_and_cover_paths() {
        let gfa = include_str!("../tests/test_data/crush/nested_bubbles_real.gfa");
        let graph = parse_gfa(gfa);
        let positions = compute_path_positions(&graph);
        let config = SmoothConfig {
            block_source: SmoothBlockSource::Flubble,
            flubble_reference_names: vec!["CHM13#0#chr6".to_string()],
            ..SmoothConfig::new(graph.paths.len())
        };

        let plan = flubble_guided_blocks(&graph, &positions, &config).unwrap();
        assert!(plan.povu_sites > 0);
        assert!(
            plan.level_ge_1_sites > 0,
            "fixture must exercise nested flubbles"
        );
        assert!(
            plan.selected_flubble_blocks > 0,
            "expected at least one level>=1 flubble block"
        );
        assert!(plan.gap_blocks > 0, "outside-flubble gaps must be retained");

        let mut covered: Vec<BitVec> = graph
            .paths
            .iter()
            .map(|path| bitvec![0; path.steps.len()])
            .collect();
        for block in &plan.blocks {
            for range in &block.path_ranges {
                for step_idx in range.begin_step..range.end_step {
                    assert!(
                        !covered[range.path_idx][step_idx],
                        "step covered by more than one block"
                    );
                    covered[range.path_idx].set(step_idx, true);
                }
            }
        }
        for (path_idx, path) in graph.paths.iter().enumerate() {
            assert_eq!(
                covered[path_idx].count_ones(),
                path.steps.len(),
                "path {} coverage",
                path.name
            );
        }
    }

    #[test]
    fn test_flubble_guided_blocks_errors_when_reference_hint_misses() {
        let gfa = include_str!("../tests/test_data/crush/nested_bubbles_real.gfa");
        let graph = parse_gfa(gfa);
        let positions = compute_path_positions(&graph);
        let config = SmoothConfig {
            block_source: SmoothBlockSource::Flubble,
            flubble_reference_names: vec!["GRCh38#0#chr6".to_string()],
            ..SmoothConfig::new(graph.paths.len())
        };

        let err = match flubble_guided_blocks(&graph, &positions, &config) {
            Ok(_) => panic!("explicit reference hint miss should fail instead of falling back"),
            Err(err) => err,
        };
        assert!(
            err.to_string().contains("GRCh38#0#chr6"),
            "explicit reference hint miss should fail instead of falling back: {err}"
        );
    }

    #[test]
    fn test_flubble_guided_smooth_preserves_path_sequences_on_nested_fixture() {
        let gfa = include_str!("../tests/test_data/crush/nested_bubbles_real.gfa");
        let before = keyed_path_sequences(gfa);
        let config = SmoothConfig {
            target_poa_lengths: vec![700],
            block_source: SmoothBlockSource::Flubble,
            flubble_reference_names: vec!["CHM13#0#chr6".to_string()],
            pre_sorted: true,
            num_threads: 1,
            ..SmoothConfig::new(5)
        };

        let smoothed = smooth_gfa(gfa, &config).unwrap();
        let after = keyed_path_sequences(&smoothed);
        assert_eq!(before, after);
    }

    #[test]
    fn test_neighbor_site_groups_merge_adjacent_sites_up_to_cap() {
        fn site(start: usize, end: usize, span: usize) -> NeighborSiteCandidate {
            NeighborSiteCandidate {
                entry: (start, false),
                exit: (end, false),
                level: 0,
                reference_start_step: start,
                reference_end_step: end,
                reference_span_bp: span,
            }
        }

        let root_positions = vec![0, 100, 200, 300, 400, 500, 600];
        let groups = greedy_neighbor_site_groups(
            vec![site(0, 1, 200), site(1, 2, 200), site(4, 5, 200)],
            &root_positions,
            350,
        );

        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].sites.len(), 2);
        assert_eq!(groups[0].reference_span_bp, 300);
        assert_eq!(groups[1].sites.len(), 1);
    }

    #[test]
    fn test_neighbor_merge_blocks_select_povu_sites_on_nested_fixture() {
        let gfa = include_str!("../tests/test_data/crush/nested_bubbles_real.gfa");
        let graph = parse_gfa(gfa);
        let positions = compute_path_positions(&graph);
        let config = SmoothConfig {
            block_source: SmoothBlockSource::NeighborMergePoasta,
            flubble_reference_names: vec!["CHM13#0#chr6".to_string()],
            ..SmoothConfig::new(graph.paths.len())
        };

        let plan = neighbor_merge_blocks(&graph, &positions, &config, 10_000).unwrap();
        assert!(plan.povu_sites > 0);
        assert!(plan.candidate_sites > 0);
        assert!(
            plan.selected_sites > 0,
            "expected at least one <=10kb neighbor-merge candidate"
        );
        assert!(
            plan.neighbor_blocks > 0,
            "selected candidates should materialize into local blocks"
        );
        assert!(
            plan.gap_blocks > 0,
            "passthrough gaps should retain coverage"
        );
    }

    #[test]
    fn test_neighbor_merge_poasta_preserves_path_sequences_on_nested_fixture() {
        let gfa = include_str!("../tests/test_data/crush/nested_bubbles_real.gfa");
        let before = keyed_path_sequences(gfa);
        let config = SmoothConfig {
            target_poa_lengths: vec![10_000],
            block_source: SmoothBlockSource::NeighborMergePoasta,
            flubble_reference_names: vec!["CHM13#0#chr6".to_string()],
            pre_sorted: true,
            num_threads: 1,
            ..SmoothConfig::new(5)
        };

        let smoothed = smooth_gfa(gfa, &config).unwrap();
        let after = keyed_path_sequences(&smoothed);
        assert_eq!(before, after);
    }

    fn keyed_path_sequences(gfa: &str) -> FxHashMap<String, String> {
        let graph = parse_gfa(gfa);
        graph
            .paths
            .iter()
            .map(|path| {
                let (key, _) = parse_path_coords(&path.name);
                let seq = path
                    .steps
                    .iter()
                    .flat_map(|&(node_idx, is_reverse)| {
                        if is_reverse {
                            reverse_complement(&graph.nodes[node_idx])
                        } else {
                            graph.nodes[node_idx].clone()
                        }
                    })
                    .collect::<Vec<_>>();
                (key, String::from_utf8(seq).unwrap())
            })
            .collect()
    }
}
