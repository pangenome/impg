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
use std::io;

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

struct Block {
    path_ranges: Vec<PathRange>,
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
        info!(
            "[smooth] Pass {}/{} (target_poa_length={}bp)",
            pass_idx + 1,
            n_passes,
            target_poa_length
        );
        current = smooth_gfa_pass(
            &current,
            target_poa_length,
            is_first && config.pre_sorted,
            config,
        )?;
    }
    Ok(current)
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
) -> io::Result<String> {
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
        return Ok(gfa_content.to_string());
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
    let mut blocks = smoothable_blocks(
        &graph,
        &node_step_index,
        max_block_weight,
        target_poa_length,
    );
    info!("[smooth] Decomposed into {} blocks", blocks.len());

    // Step 7: Break long blocks
    break_blocks(&mut blocks, &graph, max_poa_length);
    info!(
        "[smooth] After breaking: {} blocks, {} total path ranges",
        blocks.len(),
        blocks.iter().map(|b| b.path_ranges.len()).sum::<usize>()
    );

    // Step 8: Per-block SPOA smoothing (parallelized)
    use rayon::prelude::*;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.num_threads)
        .build()
        .map_err(|e| io::Error::other(format!("failed to build thread pool: {}", e)))?;
    let block_gfas: Vec<String> = pool.install(|| {
        blocks
            .par_iter()
            .enumerate()
            .filter_map(
                |(i, block)| match smooth_block(block, &graph, config, &path_positions) {
                    Ok(gfa) if !gfa.is_empty() => Some(gfa),
                    Ok(_) => None,
                    Err(e) => {
                        log::warn!("[smooth] Block {} failed: {}", i, e);
                        None
                    }
                },
            )
            .collect()
    });

    info!("[smooth] Smoothed {} blocks successfully", block_gfas.len());

    if block_gfas.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }

    // Step 9: Lace block GFAs together
    // Per-block gfaffix already normalized each block; the caller's normalize_and_sort
    // (gfaffix + sort) handles any remaining redundancies introduced at lacing boundaries.
    let laced = crate::commands::lace::lace_subgraphs(&block_gfas, config.temp_dir.as_deref())?;
    info!("[smooth] Lacing complete");

    Ok(laced)
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
    if path_ranges.len() <= 1 {
        return vec![Block { path_ranges }];
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
        return vec![Block { path_ranges }];
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
            .map(|path_ranges| Block { path_ranges })
            .collect();
    }

    // Return one block per component
    component_ranges
        .into_values()
        .filter(|ranges| !ranges.is_empty())
        .map(|path_ranges| Block { path_ranges })
        .collect()
}

// ---------------------------------------------------------------------------
// Block breaking (smoothxg break_blocks)
// ---------------------------------------------------------------------------

fn break_blocks(blocks: &mut Vec<Block>, graph: &SmoothGraph, max_poa_length: usize) {
    let mut new_blocks: Vec<Block> = Vec::with_capacity(blocks.len());

    for block in blocks.drain(..) {
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

    // 3. Run SPOA on all sequences
    let (mut spoa_graph, mut spoa_engine) = build_spoa_engine(config.scoring_params);
    feed_sequences_to_graph(
        &mut spoa_engine,
        &mut spoa_graph,
        entries.iter().map(|e| e.sequence.as_str()),
    );

    // 4. Handle padding trimming
    let any_padding = entries.iter().any(|e| e.left_pad > 0 || e.right_pad > 0);

    if !any_padding {
        // No padding — generate GFA directly
        let headers: Vec<String> = entries.iter().map(|e| e.path_name.clone()).collect();
        let gfa = spoa_graph.generate_gfa(&headers, false);
        return unchop_gfa(&gfa);
    }

    // With padding: use MSA-based trimming approach
    let msa = spoa_graph.generate_msa();

    // Find core column range (where all sequences have consumed their padding)
    let (core_start, core_end) = find_core_column_range(&msa, &entries);

    if core_start >= core_end {
        // Degenerate case — return empty
        return Ok(String::new());
    }

    // Extract core subsequences and re-run SPOA
    let nseqs = msa.len().min(entries.len());
    let msa_bytes: Vec<&[u8]> = msa.iter().map(|s| s.as_bytes()).collect();

    let mut core_sequences: Vec<String> = Vec::with_capacity(nseqs);
    for i in 0..nseqs {
        let core_seq: String = msa_bytes[i][core_start..core_end]
            .iter()
            .filter(|&&b| b != b'-')
            .map(|&b| b as char)
            .collect();
        core_sequences.push(core_seq);
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
        for col in 0..ncols {
            if msa_bytes[i][col] != b'-' {
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
        for col in 0..ncols {
            if msa_bytes[i][col] != b'-' {
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
}
