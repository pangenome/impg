pub mod poa;

use crate::commands::align::{SparsificationStrategy, SweepgaAlignConfig};
use crate::graph::{
    prepare_sequences, sort_gfa, SequenceMetadata,
};
use crate::impg_index::ImpgIndex;
use crate::realize::poa::padded_poa_from_sequences;
use crate::sequence_index::UnifiedSequenceIndex;
use coitrees::Interval;
use log::{info, warn};
use rayon::prelude::*;
use std::io::{self, BufRead, BufReader};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Configuration for the realize engine.
pub struct RealizeConfig {
    /// Regions at or below this size (bp) go directly to POA.
    /// Above this, sweepga + recursive decomposition is used.
    /// Default: 1000
    pub poa_threshold: usize,

    /// Target chunk size (bp) when partitioning large regions.
    /// Default: 5000
    pub chunk_size: usize,

    /// Boundary padding (bp) added on each side of POA regions.
    /// Default: 100
    pub padding: usize,

    /// Maximum recursion depth.
    /// Default: 10
    pub max_depth: usize,

    /// Number of threads for parallel operations.
    pub num_threads: usize,

    /// SPOA scoring parameters: (match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2)
    /// Default: (5, 4, 6, 2, 24, 1)
    pub scoring_params: (u8, u8, u8, u8, u8, u8),

    /// Optional temp directory for intermediate files.
    pub temp_dir: Option<String>,

    /// Whether to sort the final GFA using gfasort.
    pub sort_output: bool,
}

impl Default for RealizeConfig {
    fn default() -> Self {
        RealizeConfig {
            poa_threshold: 1000,
            chunk_size: 5000,
            padding: 100,
            max_depth: 10,
            num_threads: 4,
            scoring_params: (5, 4, 6, 2, 24, 1),
            temp_dir: None,
            sort_output: true,
        }
    }
}

/// Result of realizing a region.
pub struct RealizeResult {
    /// The GFA graph as a string.
    pub gfa: String,
    /// Statistics about the realization.
    pub stats: RealizeStats,
}

/// Statistics collected during realize.
pub struct RealizeStats {
    /// Total sequences in the region.
    pub num_sequences: usize,
    /// Maximum recursion depth reached.
    pub max_depth_reached: usize,
    /// Number of POA leaf calls.
    pub poa_calls: usize,
    /// Number of sweepga alignment calls.
    pub sweepga_calls: usize,
    /// Total time in milliseconds.
    pub total_ms: u64,
}

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

/// A chunk of subsequences covering one window of the anchor sequence.
struct Chunk {
    /// Subsequences for this chunk (one per input sequence, possibly fewer if
    /// some sequences don't span this region).
    sequences: Vec<(String, SequenceMetadata)>,
    /// Position in the anchor coordinate system.
    #[allow(dead_code)]
    anchor_start: usize,
    /// Exclusive end position in the anchor coordinate system.
    #[allow(dead_code)]
    anchor_end: usize,
}

/// A parsed PAF record with the fields we need for partitioning.
struct PafRecord {
    query_name: String,
    #[allow(dead_code)]
    query_len: usize,
    query_start: usize,
    query_end: usize,
    target_name: String,
    #[allow(dead_code)]
    target_len: usize,
    target_start: usize,
    target_end: usize,
    #[allow(dead_code)]
    strand: char,
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Realize a set of homologous intervals into a GFA variation graph.
///
/// This is the main entry point for the realize engine. It extracts sequences
/// from the intervals, then delegates to the recursive orchestrator.
pub fn realize(
    impg: &impl ImpgIndex,
    intervals: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    config: &RealizeConfig,
) -> io::Result<RealizeResult> {
    let start = Instant::now();

    // Extract sequences from intervals (parallel, strand-normalized, sorted longest-first).
    let sequences = prepare_sequences(impg, intervals, sequence_index)?;

    if sequences.is_empty() {
        return Ok(RealizeResult {
            gfa: String::from("H\tVN:Z:1.0\n"),
            stats: RealizeStats {
                num_sequences: 0,
                max_depth_reached: 0,
                poa_calls: 0,
                sweepga_calls: 0,
                total_ms: start.elapsed().as_millis() as u64,
            },
        });
    }

    let num_sequences = sequences.len();

    let poa_calls = AtomicUsize::new(0);
    let sweepga_calls = AtomicUsize::new(0);
    let max_depth_reached = AtomicUsize::new(0);

    let gfa = realize_recursive(&sequences, config, 0, &poa_calls, &sweepga_calls, &max_depth_reached)?;

    // Optionally sort the final GFA.
    let gfa = if config.sort_output {
        sort_gfa(&gfa, config.num_threads)?
    } else {
        gfa
    };

    Ok(RealizeResult {
        gfa,
        stats: RealizeStats {
            num_sequences,
            max_depth_reached: max_depth_reached.load(Ordering::Relaxed),
            poa_calls: poa_calls.load(Ordering::Relaxed),
            sweepga_calls: sweepga_calls.load(Ordering::Relaxed),
            total_ms: start.elapsed().as_millis() as u64,
        },
    })
}

/// Realize directly from pre-extracted sequences (for use by other modules that
/// already have sequences in memory).
pub fn realize_from_sequences(
    sequences: &[(String, SequenceMetadata)],
    config: &RealizeConfig,
) -> io::Result<RealizeResult> {
    let start = Instant::now();

    if sequences.is_empty() {
        return Ok(RealizeResult {
            gfa: String::from("H\tVN:Z:1.0\n"),
            stats: RealizeStats {
                num_sequences: 0,
                max_depth_reached: 0,
                poa_calls: 0,
                sweepga_calls: 0,
                total_ms: start.elapsed().as_millis() as u64,
            },
        });
    }

    let num_sequences = sequences.len();
    let poa_calls = AtomicUsize::new(0);
    let sweepga_calls = AtomicUsize::new(0);
    let max_depth_reached = AtomicUsize::new(0);

    let gfa = realize_recursive(sequences, config, 0, &poa_calls, &sweepga_calls, &max_depth_reached)?;

    let gfa = if config.sort_output {
        sort_gfa(&gfa, config.num_threads)?
    } else {
        gfa
    };

    Ok(RealizeResult {
        gfa,
        stats: RealizeStats {
            num_sequences,
            max_depth_reached: max_depth_reached.load(Ordering::Relaxed),
            poa_calls: poa_calls.load(Ordering::Relaxed),
            sweepga_calls: sweepga_calls.load(Ordering::Relaxed),
            total_ms: start.elapsed().as_millis() as u64,
        },
    })
}

// ---------------------------------------------------------------------------
// Recursive orchestrator
// ---------------------------------------------------------------------------

/// Recursive realize: decide between POA base case and sweepga + partition.
fn realize_recursive(
    sequences: &[(String, SequenceMetadata)],
    config: &RealizeConfig,
    depth: usize,
    poa_calls: &AtomicUsize,
    sweepga_calls: &AtomicUsize,
    max_depth_reached: &AtomicUsize,
) -> io::Result<String> {
    // Track max depth.
    loop {
        let current = max_depth_reached.load(Ordering::Relaxed);
        if depth <= current {
            break;
        }
        if max_depth_reached
            .compare_exchange_weak(current, depth, Ordering::Relaxed, Ordering::Relaxed)
            .is_ok()
        {
            break;
        }
    }

    if sequences.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }

    // Single sequence → trivial GFA.
    if sequences.len() == 1 {
        poa_calls.fetch_add(1, Ordering::Relaxed);
        return make_trivial_gfa(&sequences[0]);
    }

    // Determine the longest sequence span (the anchor).
    let total_span = sequences.iter().map(|(s, _)| s.len()).max().unwrap_or(0);

    // Base case: small region or max depth → POA.
    if total_span <= config.poa_threshold || depth >= config.max_depth {
        if depth >= config.max_depth && total_span > config.poa_threshold {
            warn!(
                "realize: hit max depth {} with region of {}bp (poa_threshold={}), forcing POA",
                depth, total_span, config.poa_threshold
            );
        }
        poa_calls.fetch_add(1, Ordering::Relaxed);
        let result = padded_poa_from_sequences(sequences, config.scoring_params, config.padding)?;
        return Ok(result.gfa);
    }

    // Recursive case: sweepga align → partition → recurse → lace.
    sweepga_calls.fetch_add(1, Ordering::Relaxed);

    // Always use all-vs-all alignment in the realize engine.
    // FastGA is efficient when run all-vs-all on a combined FASTA (indexes once).
    // Sparsification is for the top-level `impg align` CLI with hundreds of genomes,
    // not for the small sequence sets inside recursive chunks.
    let align_config = build_sweepga_config(config);
    let named_seqs: Vec<(String, &[u8])> = sequences
        .iter()
        .map(|(seq, meta)| (meta.alignment_name(), seq.as_bytes()))
        .collect();

    let paf_temp = crate::commands::align::sweepga_align(&named_seqs, &align_config)?;

    // Parse PAF records for partitioning.
    let paf_records = parse_paf_file(paf_temp.path())?;

    // If no alignments were produced, fall back to POA.
    if paf_records.is_empty() {
        warn!(
            "realize: sweepga produced no alignments at depth {} for {}bp region with {} sequences, falling back to POA",
            depth, total_span, sequences.len()
        );
        poa_calls.fetch_add(1, Ordering::Relaxed);
        let result = padded_poa_from_sequences(sequences, config.scoring_params, config.padding)?;
        return Ok(result.gfa);
    }

    // Partition into chunks along the anchor (longest) sequence.
    let chunks = partition_into_chunks(sequences, &paf_records, config.chunk_size, config.padding);

    // If partitioning produced a single chunk that is no smaller, use POA directly.
    if chunks.len() <= 1 {
        poa_calls.fetch_add(1, Ordering::Relaxed);
        let result = padded_poa_from_sequences(sequences, config.scoring_params, config.padding)?;
        return Ok(result.gfa);
    }

    info!(
        "realize: depth={}, span={}bp, {} seqs → {} chunks",
        depth,
        total_span,
        sequences.len(),
        chunks.len()
    );

    // Recurse on each chunk in parallel.
    let sub_gfa_results: Vec<Option<io::Result<String>>> = chunks
        .par_iter()
        .map(|chunk| {
            if chunk.sequences.is_empty() {
                return None;
            }
            Some(realize_recursive(
                &chunk.sequences,
                config,
                depth + 1,
                poa_calls,
                sweepga_calls,
                max_depth_reached,
            ))
        })
        .collect();

    // Collect results, propagating any errors.
    let mut sub_gfas: Vec<String> = Vec::with_capacity(sub_gfa_results.len());
    for result in sub_gfa_results {
        if let Some(gfa_result) = result {
            sub_gfas.push(gfa_result?);
        }
    }

    if sub_gfas.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }

    if sub_gfas.len() == 1 {
        return Ok(sub_gfas.into_iter().next().unwrap());
    }

    // Lace sub-graphs together.
    lace_subgraphs(&sub_gfas, config.num_threads)
}

// ---------------------------------------------------------------------------
// Sweepga alignment helpers
// ---------------------------------------------------------------------------

/// Build a SweepgaAlignConfig from RealizeConfig.
///
/// Uses `no_filter: true` because the realize engine only needs raw alignment
/// coordinates for window projection — aggressive PAF filtering (overlap,
/// scaffolding, 1:1 mapping) would remove valid cross-sequence alignments that
/// are needed for partitioning. Self-alignments are still filtered out during
/// PAF parsing since they share query_name == target_name.
fn build_sweepga_config(config: &RealizeConfig) -> SweepgaAlignConfig {
    SweepgaAlignConfig {
        num_threads: config.num_threads,
        kmer_frequency: 10, // Reasonable default for small sets
        min_alignment_length: 100,
        no_filter: true,
        num_mappings: "1:1".to_string(),
        scaffold_jump: 50_000,
        scaffold_mass: 10_000,
        scaffold_filter: "1:1".to_string(),
        overlap: 0.95,
        min_identity: 0.0,
        scaffold_dist: 0,
        min_mapping_length: 0,
        temp_dir: config.temp_dir.clone(),
        sparsification: SparsificationStrategy::None, // Always all-vs-all for realize
    }
}

/// Parse a PAF file into PafRecord structs, filtering out self-alignments.
fn parse_paf_file(path: &std::path::Path) -> io::Result<Vec<PafRecord>> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let query_name = fields[0].to_string();
        let target_name = fields[5].to_string();

        // Skip self-alignments (same sequence aligned to itself).
        if query_name == target_name {
            continue;
        }

        let query_len: usize = fields[1].parse().unwrap_or(0);
        let query_start: usize = fields[2].parse().unwrap_or(0);
        let query_end: usize = fields[3].parse().unwrap_or(0);
        let strand = if fields[4] == "-" { '-' } else { '+' };
        let target_len: usize = fields[6].parse().unwrap_or(0);
        let target_start: usize = fields[7].parse().unwrap_or(0);
        let target_end: usize = fields[8].parse().unwrap_or(0);

        records.push(PafRecord {
            query_name,
            query_len,
            query_start,
            query_end,
            target_name,
            target_len,
            target_start,
            target_end,
            strand,
        });
    }

    Ok(records)
}

// ---------------------------------------------------------------------------
// Chunk partitioning
// ---------------------------------------------------------------------------

/// Partition sequences into chunks along the anchor (longest) sequence.
///
/// Uses the PAF alignments to project windows from the anchor coordinate space
/// onto each other sequence. Adjacent chunks overlap by `padding` bp.
fn partition_into_chunks(
    sequences: &[(String, SequenceMetadata)],
    paf_records: &[PafRecord],
    chunk_size: usize,
    padding: usize,
) -> Vec<Chunk> {
    if sequences.is_empty() {
        return Vec::new();
    }

    // The anchor is the longest sequence (sequences are pre-sorted longest-first).
    let (anchor_seq, anchor_meta) = &sequences[0];
    let anchor_len = anchor_seq.len();
    let anchor_name = anchor_meta.alignment_name();

    // If the anchor is smaller than one chunk, return a single chunk with all sequences.
    if anchor_len <= chunk_size + padding {
        return vec![Chunk {
            sequences: sequences.to_vec(),
            anchor_start: 0,
            anchor_end: anchor_len,
        }];
    }

    // Build name-to-index map for the input sequences (using alignment names to match PAF).
    let seq_name_map: std::collections::HashMap<String, usize> = sequences
        .iter()
        .enumerate()
        .map(|(i, (_, meta))| (meta.alignment_name(), i))
        .collect();

    // For each non-anchor sequence, collect alignment projections from the anchor.
    // projection[i] = Vec<(anchor_start, anchor_end, seq_start, seq_end, is_reverse)>
    let mut projections: Vec<Vec<(usize, usize, usize, usize, bool)>> =
        vec![Vec::new(); sequences.len()];

    for rec in paf_records {
        let is_reverse = rec.strand == '-';
        // We want alignments where the anchor is either query or target.
        if rec.query_name == anchor_name {
            if let Some(&idx) = seq_name_map.get(&rec.target_name) {
                projections[idx].push((
                    rec.query_start,
                    rec.query_end,
                    rec.target_start,
                    rec.target_end,
                    is_reverse,
                ));
            }
        } else if rec.target_name == anchor_name {
            if let Some(&idx) = seq_name_map.get(&rec.query_name) {
                projections[idx].push((
                    rec.target_start,
                    rec.target_end,
                    rec.query_start,
                    rec.query_end,
                    is_reverse,
                ));
            }
        }
    }

    // Sort projections by anchor start for binary search.
    for proj in &mut projections {
        proj.sort_by_key(|&(a_start, _, _, _, _)| a_start);
    }

    // Create windows along the anchor with overlap.
    let mut windows: Vec<(usize, usize)> = Vec::new();
    let mut pos = 0;
    while pos < anchor_len {
        let win_start = if pos > padding { pos - padding } else { 0 };
        let win_end = (pos + chunk_size + padding).min(anchor_len);
        windows.push((win_start, win_end));
        pos += chunk_size;
    }

    // For each window, extract the subsequences.
    let mut chunks: Vec<Chunk> = Vec::with_capacity(windows.len());

    for (win_start, win_end) in &windows {
        let mut chunk_seqs: Vec<(String, SequenceMetadata)> = Vec::new();

        // Always include the anchor's slice for this window.
        // Preserve the original metadata so that path_name() is consistent
        // across all chunks — lacing merges paths by name.
        let anchor_slice = &anchor_seq[*win_start..*win_end];
        chunk_seqs.push((anchor_slice.to_string(), anchor_meta.clone()));

        // For each non-anchor sequence, project the window through alignments.
        // If a sequence has no projection from the anchor (no direct alignment),
        // include it fully in every chunk as a fallback.
        for (seq_idx, (seq_str, seq_meta)) in sequences.iter().enumerate().skip(1) {
            if projections[seq_idx].is_empty() {
                // No anchor alignment — include full sequence in every chunk.
                chunk_seqs.push((seq_str.clone(), seq_meta.clone()));
            } else if let Some((seq_start, seq_end)) =
                project_window_to_sequence(&projections[seq_idx], *win_start, *win_end)
            {
                let clamped_start = seq_start.min(seq_str.len());
                let clamped_end = seq_end.min(seq_str.len());
                if clamped_end > clamped_start {
                    let sub_seq = &seq_str[clamped_start..clamped_end];
                    // Preserve the original metadata so path_name() is
                    // consistent across chunks for lacing.
                    chunk_seqs.push((sub_seq.to_string(), seq_meta.clone()));
                }
            }
        }

        // Sort chunk sequences by length descending (expected by POA).
        chunk_seqs.sort_by(|a, b| b.0.len().cmp(&a.0.len()));

        chunks.push(Chunk {
            sequences: chunk_seqs,
            anchor_start: *win_start,
            anchor_end: *win_end,
        });
    }

    chunks
}

/// Project an anchor window [win_start, win_end) onto a target sequence
/// using the alignment projections.
///
/// Returns the (start, end) range in the target sequence coordinates, or None
/// if no alignment covers this window.
///
/// For reverse-strand alignments, the mapping is inverted: moving forward in
/// anchor space maps to moving backward in target sequence space.
fn project_window_to_sequence(
    projections: &[(usize, usize, usize, usize, bool)],
    win_start: usize,
    win_end: usize,
) -> Option<(usize, usize)> {
    let mut best_start = usize::MAX;
    let mut best_end = 0usize;

    for &(a_start, a_end, s_start, s_end, is_reverse) in projections {
        // Check if this alignment overlaps the window.
        if a_start >= win_end || a_end <= win_start {
            continue;
        }

        // Linear interpolation: map [win_start, win_end) through this alignment block.
        let a_span = a_end.saturating_sub(a_start).max(1);
        let s_span = s_end.saturating_sub(s_start);

        let overlap_start = win_start.max(a_start);
        let overlap_end = win_end.min(a_end);

        let frac_start = (overlap_start - a_start) as f64 / a_span as f64;
        let frac_end = (overlap_end - a_start) as f64 / a_span as f64;

        let (projected_start, projected_end) = if is_reverse {
            // Reverse strand: moving forward in anchor space → backward in seq space.
            // frac_start (closer to a_start) maps to s_end, frac_end maps toward s_start.
            let p_end = s_end - (frac_start * s_span as f64) as usize;
            let p_start = s_end - (frac_end * s_span as f64) as usize;
            (p_start, p_end)
        } else {
            // Forward strand: straightforward linear interpolation.
            let p_start = s_start + (frac_start * s_span as f64) as usize;
            let p_end = s_start + (frac_end * s_span as f64) as usize;
            (p_start, p_end)
        };

        best_start = best_start.min(projected_start);
        best_end = best_end.max(projected_end);
    }

    if best_start < best_end {
        Some((best_start, best_end))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Lacing sub-graphs
// ---------------------------------------------------------------------------

/// Lace multiple GFA sub-graphs into a single GFA.
///
/// This is a simplified version that:
/// 1. Remaps node IDs to a unified namespace.
/// 2. Concatenates segments and links.
/// 3. Joins paths that share the same name across chunks by linking their
///    last/first nodes.
fn lace_subgraphs(sub_gfas: &[String], _num_threads: usize) -> io::Result<String> {
    if sub_gfas.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }

    if sub_gfas.len() == 1 {
        return Ok(sub_gfas[0].clone());
    }

    let mut next_id: u64 = 1;
    let mut all_segments: Vec<String> = Vec::new();
    let mut all_links: Vec<String> = Vec::new();

    // For each sub-GFA: map of path_name -> ordered list of (segment_id, orientation) tuples.
    // We collect these across all sub-GFAs to build the final merged paths.
    let mut path_segments_per_gfa: Vec<std::collections::HashMap<String, Vec<String>>> =
        Vec::with_capacity(sub_gfas.len());

    for gfa_str in sub_gfas {
        let mut id_remap: std::collections::HashMap<String, u64> = std::collections::HashMap::new();
        let mut local_path_segments: std::collections::HashMap<String, Vec<String>> =
            std::collections::HashMap::new();

        for line in gfa_str.lines() {
            if line.starts_with("S\t") {
                let fields: Vec<&str> = line.splitn(4, '\t').collect();
                if fields.len() >= 3 {
                    let old_id = fields[1].to_string();
                    let new_id = next_id;
                    next_id += 1;
                    id_remap.insert(old_id, new_id);
                    // Re-emit segment with new ID.
                    if fields.len() == 3 {
                        all_segments.push(format!("S\t{}\t{}", new_id, fields[2]));
                    } else {
                        all_segments.push(format!("S\t{}\t{}\t{}", new_id, fields[2], fields[3]));
                    }
                }
            } else if line.starts_with("L\t") {
                let fields: Vec<&str> = line.splitn(7, '\t').collect();
                if fields.len() >= 6 {
                    let from_id = id_remap
                        .get(fields[1])
                        .copied()
                        .unwrap_or_else(|| fields[1].parse().unwrap_or(0));
                    let to_id = id_remap
                        .get(fields[3])
                        .copied()
                        .unwrap_or_else(|| fields[3].parse().unwrap_or(0));
                    if fields.len() == 6 {
                        all_links.push(format!(
                            "L\t{}\t{}\t{}\t{}\t{}",
                            from_id, fields[2], to_id, fields[4], fields[5]
                        ));
                    } else {
                        all_links.push(format!(
                            "L\t{}\t{}\t{}\t{}\t{}\t{}",
                            from_id, fields[2], to_id, fields[4], fields[5], fields[6]
                        ));
                    }
                }
            } else if line.starts_with("P\t") {
                let fields: Vec<&str> = line.splitn(4, '\t').collect();
                if fields.len() >= 3 {
                    let path_name = fields[1].to_string();
                    let steps: Vec<String> = fields[2]
                        .split(',')
                        .map(|step| {
                            // Remap node ID in step (e.g. "5+" → "12+").
                            let (node_str, orient) =
                                if let Some(stripped) = step.strip_suffix('+') {
                                    (stripped, "+")
                                } else if let Some(stripped) = step.strip_suffix('-') {
                                    (stripped, "-")
                                } else {
                                    (step, "+")
                                };
                            let new_node = id_remap
                                .get(node_str)
                                .copied()
                                .unwrap_or_else(|| node_str.parse().unwrap_or(0));
                            format!("{}{}", new_node, orient)
                        })
                        .collect();
                    local_path_segments
                        .entry(path_name)
                        .or_default()
                        .extend(steps);
                }
            }
        }

        path_segments_per_gfa.push(local_path_segments);
    }

    // Merge paths: for each path name, concatenate segments from all sub-GFAs
    // and add linking edges between the last node of chunk N and first node of chunk N+1.
    let mut merged_paths: std::collections::HashMap<String, Vec<String>> =
        std::collections::HashMap::new();

    // Collect all path names in order of first appearance.
    let mut path_order: Vec<String> = Vec::new();
    {
        let mut seen = std::collections::HashSet::new();
        for gfa_paths in &path_segments_per_gfa {
            for name in gfa_paths.keys() {
                if seen.insert(name.clone()) {
                    path_order.push(name.clone());
                }
            }
        }
    }

    for path_name in &path_order {
        let mut full_path: Vec<String> = Vec::new();

        for gfa_paths in &path_segments_per_gfa {
            if let Some(steps) = gfa_paths.get(path_name) {
                if !steps.is_empty() {
                    // Add a linking edge between the last step of the accumulated path
                    // and the first step of this chunk.
                    if let (Some(last_step), Some(first_step)) = (full_path.last(), steps.first())
                    {
                        let (last_node, last_orient) = parse_step(last_step);
                        let (first_node, first_orient) = parse_step(first_step);
                        all_links.push(format!(
                            "L\t{}\t{}\t{}\t{}\t0M",
                            last_node, last_orient, first_node, first_orient
                        ));
                    }
                    full_path.extend(steps.iter().cloned());
                }
            }
        }

        if !full_path.is_empty() {
            merged_paths.insert(path_name.clone(), full_path);
        }
    }

    // Assemble final GFA.
    let mut output = String::new();
    output.push_str("H\tVN:Z:1.0\n");

    for seg in &all_segments {
        output.push_str(seg);
        output.push('\n');
    }

    // Deduplicate links while preserving insertion order for determinism.
    let mut seen_links = std::collections::HashSet::new();
    for link in &all_links {
        if seen_links.insert(link) {
            output.push_str(link);
            output.push('\n');
        }
    }

    for path_name in &path_order {
        if let Some(steps) = merged_paths.get(path_name) {
            output.push_str(&format!("P\t{}\t{}\t*\n", path_name, steps.join(",")));
        }
    }

    Ok(output)
}

/// Parse a GFA path step like "12+" into (node_id_str, orientation_str).
fn parse_step(step: &str) -> (&str, &str) {
    if let Some(stripped) = step.strip_suffix('+') {
        (stripped, "+")
    } else if let Some(stripped) = step.strip_suffix('-') {
        (stripped, "-")
    } else {
        (step, "+")
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Create a trivial GFA for a single sequence: one segment, one path.
fn make_trivial_gfa(seq_data: &(String, SequenceMetadata)) -> io::Result<String> {
    let (seq, meta) = seq_data;
    Ok(format!(
        "H\tVN:Z:1.0\nS\t1\t{}\nP\t{}\t1+\t*\n",
        seq, meta.path_name()
    ))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_seq(name: &str, seq: &str, start: i32, total: usize) -> (String, SequenceMetadata) {
        (
            seq.to_string(),
            SequenceMetadata {
                name: name.to_string(),
                start,
                size: seq.len() as i32,
                strand: '+',
                total_length: total,
            },
        )
    }

    #[test]
    fn test_make_trivial_gfa() {
        let seq = make_seq("chr1", "ACGTACGT", 100, 1000);
        let gfa = make_trivial_gfa(&seq).unwrap();
        assert!(gfa.contains("H\tVN:Z:1.0"));
        assert!(gfa.contains("S\t1\tACGTACGT"));
        assert!(gfa.contains("P\tchr1:100-108\t1+\t*"));
    }

    #[test]
    fn test_realize_from_sequences_empty() {
        let config = RealizeConfig::default();
        let result = realize_from_sequences(&[], &config).unwrap();
        assert!(result.gfa.contains("H\tVN:Z:1.0"));
        assert_eq!(result.stats.num_sequences, 0);
    }

    #[test]
    fn test_realize_from_sequences_single() {
        let mut config = RealizeConfig::default();
        config.sort_output = false;
        let seqs = vec![make_seq("s1", "ACGTACGTACGT", 0, 100)];
        let result = realize_from_sequences(&seqs, &config).unwrap();
        assert!(result.gfa.contains("S\t"));
        assert!(result.gfa.contains("P\t"));
        assert_eq!(result.stats.num_sequences, 1);
        assert_eq!(result.stats.poa_calls, 1);
    }

    #[test]
    fn test_realize_from_sequences_small_region_uses_poa() {
        // Two short sequences → should go directly to POA (below poa_threshold).
        let mut config = RealizeConfig::default();
        config.poa_threshold = 1000;
        config.padding = 0; // No padding for tiny test sequences.
        config.sort_output = false;

        let seqs = vec![
            make_seq("s1", "ACGTACGTACGTACGTACGT", 0, 100),
            make_seq("s2", "ACGTACGTACGTACGTACGT", 10, 100),
        ];
        let result = realize_from_sequences(&seqs, &config).unwrap();
        assert!(result.gfa.contains("S\t"));
        assert!(result.gfa.contains("P\t"));
        assert_eq!(result.stats.num_sequences, 2);
        assert_eq!(result.stats.poa_calls, 1);
        assert_eq!(result.stats.sweepga_calls, 0);
    }

    #[test]
    fn test_project_window_to_sequence_no_overlap() {
        let projections = vec![(100, 200, 50, 150, false)];
        // Window completely before the alignment.
        assert_eq!(project_window_to_sequence(&projections, 0, 50), None);
        // Window completely after the alignment.
        assert_eq!(project_window_to_sequence(&projections, 300, 400), None);
    }

    #[test]
    fn test_project_window_to_sequence_full_overlap() {
        // Alignment: anchor [100,200) → seq [50,150)
        let projections = vec![(100, 200, 50, 150, false)];
        // Window exactly covers the alignment.
        let result = project_window_to_sequence(&projections, 100, 200);
        assert_eq!(result, Some((50, 150)));
    }

    #[test]
    fn test_project_window_to_sequence_partial_overlap() {
        // Alignment: anchor [100,200) → seq [0,100)
        let projections = vec![(100, 200, 0, 100, false)];
        // Window covers first half: anchor [100,150) → seq [0,50)
        let result = project_window_to_sequence(&projections, 100, 150);
        assert_eq!(result, Some((0, 50)));
    }

    #[test]
    fn test_project_window_reverse_strand_full_overlap() {
        // Reverse-strand alignment: anchor [100,200) → seq [50,150) but inverted.
        // Moving forward in anchor → backward in seq.
        // So anchor 100 → seq 150, anchor 200 → seq 50.
        let projections = vec![(100, 200, 50, 150, true)];
        let result = project_window_to_sequence(&projections, 100, 200);
        assert_eq!(result, Some((50, 150)));
    }

    #[test]
    fn test_project_window_reverse_strand_first_half() {
        // Reverse-strand: anchor [100,200) → seq [0,100) inverted.
        // anchor 100 → seq 100, anchor 200 → seq 0.
        // Window [100,150): frac_start=0.0, frac_end=0.5
        // projected: p_end = 100 - 0.0*100 = 100, p_start = 100 - 0.5*100 = 50
        let projections = vec![(100, 200, 0, 100, true)];
        let result = project_window_to_sequence(&projections, 100, 150);
        assert_eq!(result, Some((50, 100)));
    }

    #[test]
    fn test_project_window_reverse_strand_second_half() {
        // Reverse-strand: anchor [100,200) → seq [0,100) inverted.
        // Window [150,200): frac_start=0.5, frac_end=1.0
        // projected: p_end = 100 - 0.5*100 = 50, p_start = 100 - 1.0*100 = 0
        let projections = vec![(100, 200, 0, 100, true)];
        let result = project_window_to_sequence(&projections, 150, 200);
        assert_eq!(result, Some((0, 50)));
    }

    #[test]
    fn test_project_window_mixed_strands() {
        // Two alignments: one forward, one reverse, both overlapping the window.
        // Forward: anchor [0,100) → seq [200,300)
        // Reverse: anchor [50,150) → seq [400,500)
        // Window [40,110):
        //   Forward overlap [40,100): frac 0.4..1.0 → seq [240,300)
        //   Reverse overlap [50,110): frac 0.0..0.6 → seq end=500-0=500, start=500-60=440 → [440,500)
        let projections = vec![
            (0, 100, 200, 300, false),
            (50, 150, 400, 500, true),
        ];
        let result = project_window_to_sequence(&projections, 40, 110);
        // best_start = min(240, 440) = 240, best_end = max(300, 500) = 500
        assert_eq!(result, Some((240, 500)));
    }

    #[test]
    fn test_partition_single_chunk() {
        // Sequences shorter than chunk_size → single chunk.
        let seqs = vec![
            make_seq("s1", "ACGTACGT", 0, 100),
            make_seq("s2", "ACGTACGT", 0, 100),
        ];
        let chunks = partition_into_chunks(&seqs, &[], 10000, 100);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].sequences.len(), 2);
    }

    #[test]
    fn test_lace_subgraphs_single() {
        let gfa = "H\tVN:Z:1.0\nS\t1\tACGT\nP\ts1:0-4\t1+\t*\n".to_string();
        let result = lace_subgraphs(&[gfa.clone()], 1).unwrap();
        assert_eq!(result, gfa);
    }

    #[test]
    fn test_lace_subgraphs_two() {
        let gfa1 = "H\tVN:Z:1.0\nS\t1\tACGT\nP\ts1:0-8\t1+\t*\n".to_string();
        let gfa2 = "H\tVN:Z:1.0\nS\t1\tTGCA\nP\ts1:0-8\t1+\t*\n".to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], 1).unwrap();

        // Should have header, two segments (remapped to 1 and 2), a link, and one merged path.
        assert!(result.contains("H\tVN:Z:1.0"));
        assert!(result.contains("S\t1\tACGT"));
        assert!(result.contains("S\t2\tTGCA"));
        assert!(result.contains("L\t1\t+\t2\t+\t0M"));
        assert!(result.contains("P\ts1:0-8\t1+,2+\t*"));
    }

    #[test]
    fn test_lace_subgraphs_empty() {
        let result = lace_subgraphs(&[], 1).unwrap();
        assert!(result.contains("H\tVN:Z:1.0"));
    }

    #[test]
    fn test_lace_subgraphs_reverse_orientation() {
        // GFA1 has a path step with reverse orientation; GFA2 has forward.
        let gfa1 = "H\tVN:Z:1.0\nS\t1\tACGT\nP\tp1\t1-\t*\n".to_string();
        let gfa2 = "H\tVN:Z:1.0\nS\t1\tTGCA\nP\tp1\t1+\t*\n".to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], 1).unwrap();

        assert!(result.contains("S\t1\tACGT"));
        assert!(result.contains("S\t2\tTGCA"));
        // Linking edge should preserve orientations: 1- → 2+
        assert!(result.contains("L\t1\t-\t2\t+\t0M"));
        assert!(result.contains("P\tp1\t1-,2+\t*"));
    }

    #[test]
    fn test_lace_subgraphs_multi_segment_paths() {
        // Each sub-GFA has a path with multiple segments.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tACGT\n\
                     S\t2\tTTTT\n\
                     L\t1\t+\t2\t+\t0M\n\
                     P\tp1\t1+,2+\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tGGGG\n\
                     S\t2\tCCCC\n\
                     L\t1\t+\t2\t-\t0M\n\
                     P\tp1\t1+,2-\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], 1).unwrap();

        // Segments remapped: gfa1 1→1, 2→2; gfa2 1→3, 2→4
        assert!(result.contains("S\t1\tACGT"));
        assert!(result.contains("S\t2\tTTTT"));
        assert!(result.contains("S\t3\tGGGG"));
        assert!(result.contains("S\t4\tCCCC"));
        // Internal links remapped
        assert!(result.contains("L\t1\t+\t2\t+\t0M"));
        assert!(result.contains("L\t3\t+\t4\t-\t0M"));
        // Linking edge between last of gfa1 (2+) and first of gfa2 (3+)
        assert!(result.contains("L\t2\t+\t3\t+\t0M"));
        // Merged path
        assert!(result.contains("P\tp1\t1+,2+,3+,4-\t*"));
    }

    #[test]
    fn test_lace_subgraphs_multiple_paths() {
        // Two different path names across two sub-GFAs.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tAAAA\n\
                     S\t2\tCCCC\n\
                     P\tpathA\t1+\t*\n\
                     P\tpathB\t2+\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tGGGG\n\
                     S\t2\tTTTT\n\
                     P\tpathA\t1+\t*\n\
                     P\tpathB\t2-\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], 1).unwrap();

        // Segments: gfa1 1→1, 2→2; gfa2 1→3, 2→4
        assert!(result.contains("S\t1\tAAAA"));
        assert!(result.contains("S\t2\tCCCC"));
        assert!(result.contains("S\t3\tGGGG"));
        assert!(result.contains("S\t4\tTTTT"));
        // pathA: 1+ from gfa1, 3+ from gfa2, linked 1+ → 3+
        assert!(result.contains("P\tpathA\t1+,3+\t*"));
        assert!(result.contains("L\t1\t+\t3\t+\t0M"));
        // pathB: 2+ from gfa1, 4- from gfa2, linked 2+ → 4-
        assert!(result.contains("P\tpathB\t2+,4-\t*"));
        assert!(result.contains("L\t2\t+\t4\t-\t0M"));
    }

    #[test]
    fn test_lace_subgraphs_three_chunks() {
        // Lace three sub-GFAs to verify chaining works across >2 chunks.
        let gfa1 = "H\tVN:Z:1.0\nS\t1\tAA\nP\tp\t1+\t*\n".to_string();
        let gfa2 = "H\tVN:Z:1.0\nS\t1\tCC\nP\tp\t1-\t*\n".to_string();
        let gfa3 = "H\tVN:Z:1.0\nS\t1\tGG\nP\tp\t1+\t*\n".to_string();
        let result = lace_subgraphs(&[gfa1, gfa2, gfa3], 1).unwrap();

        assert!(result.contains("S\t1\tAA"));
        assert!(result.contains("S\t2\tCC"));
        assert!(result.contains("S\t3\tGG"));
        // Linking edges: 1+ → 2-, 2- → 3+
        assert!(result.contains("L\t1\t+\t2\t-\t0M"));
        assert!(result.contains("L\t2\t-\t3\t+\t0M"));
        assert!(result.contains("P\tp\t1+,2-,3+\t*"));
    }

    #[test]
    fn test_lace_subgraphs_mixed_orientation_multi_segment() {
        // Sub-GFAs with multi-segment paths in mixed orientations across 3 chunks.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tAA\n\
                     S\t2\tCC\n\
                     L\t1\t+\t2\t-\t0M\n\
                     P\tp\t1+,2-\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tGG\n\
                     P\tp\t1-\t*\n"
            .to_string();
        let gfa3 = "H\tVN:Z:1.0\n\
                     S\t1\tTT\n\
                     S\t2\tAA\n\
                     L\t1\t-\t2\t+\t0M\n\
                     P\tp\t1-,2+\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2, gfa3], 1).unwrap();

        // Segments: gfa1 {1→1,2→2}, gfa2 {1→3}, gfa3 {1→4,2→5}
        assert!(result.contains("S\t1\tAA"));
        assert!(result.contains("S\t2\tCC"));
        assert!(result.contains("S\t3\tGG"));
        assert!(result.contains("S\t4\tTT"));
        assert!(result.contains("S\t5\tAA"));
        // Internal links remapped
        assert!(result.contains("L\t1\t+\t2\t-\t0M"));
        assert!(result.contains("L\t4\t-\t5\t+\t0M"));
        // Linking edges: chunk1 last=2- → chunk2 first=3-, chunk2 last=3- → chunk3 first=4-
        assert!(result.contains("L\t2\t-\t3\t-\t0M"));
        assert!(result.contains("L\t3\t-\t4\t-\t0M"));
        // Full merged path
        assert!(result.contains("P\tp\t1+,2-,3-,4-,5+\t*"));
    }

    #[test]
    fn test_lace_subgraphs_path_present_in_subset_of_chunks() {
        // pathA appears in both chunks, pathB only in the first chunk.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tACGT\n\
                     P\tpathA\t1+\t*\n\
                     P\tpathB\t1-\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tTGCA\n\
                     P\tpathA\t1+\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], 1).unwrap();

        // pathA should be merged with linking edge
        assert!(result.contains("P\tpathA\t1+,2+\t*"));
        assert!(result.contains("L\t1\t+\t2\t+\t0M"));
        // pathB should appear with just its single step (no linking edge needed)
        assert!(result.contains("P\tpathB\t1-\t*"));
    }

    #[test]
    fn test_lace_subgraphs_link_deduplication() {
        // Both sub-GFAs have an internal link that, after remapping, would be
        // different (1+→2+ vs 3+→4+). But a linking edge 2+→3+ is unique.
        // Also test that truly duplicate links are deduplicated.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tAA\n\
                     S\t2\tCC\n\
                     L\t1\t+\t2\t+\t0M\n\
                     L\t1\t+\t2\t+\t0M\n\
                     P\tp\t1+,2+\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tGG\n\
                     P\tp\t1+\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], 1).unwrap();

        // The duplicate L 1+→2+ should appear only once.
        let link_count = result.matches("L\t1\t+\t2\t+\t0M").count();
        assert_eq!(link_count, 1, "duplicate links should be deduplicated");
        // Linking edge 2+ → 3+
        assert!(result.contains("L\t2\t+\t3\t+\t0M"));
    }

    #[test]
    fn test_lace_subgraphs_segment_tags_preserved() {
        // Segments may have optional tag fields after the sequence.
        let gfa1 = "H\tVN:Z:1.0\nS\t1\tACGT\tLN:i:4\nP\tp\t1+\t*\n".to_string();
        let gfa2 = "H\tVN:Z:1.0\nS\t1\tTGCA\nP\tp\t1+\t*\n".to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], 1).unwrap();

        // Tag should be preserved on the first segment.
        assert!(result.contains("S\t1\tACGT\tLN:i:4"));
        assert!(result.contains("S\t2\tTGCA"));
    }

    // -----------------------------------------------------------------------
    // parse_paf_file tests
    // -----------------------------------------------------------------------

    /// Helper: write PAF content to a temp file and return its path handle.
    fn write_paf_temp(content: &str) -> tempfile::NamedTempFile {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn test_parse_paf_file_basic() {
        // Standard 12-column PAF line (mandatory fields).
        let paf = "qry1\t1000\t100\t500\t+\ttgt1\t2000\t200\t600\t350\t400\t60\n";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.query_name, "qry1");
        assert_eq!(r.query_len, 1000);
        assert_eq!(r.query_start, 100);
        assert_eq!(r.query_end, 500);
        assert_eq!(r.strand, '+');
        assert_eq!(r.target_name, "tgt1");
        assert_eq!(r.target_len, 2000);
        assert_eq!(r.target_start, 200);
        assert_eq!(r.target_end, 600);
    }

    #[test]
    fn test_parse_paf_file_reverse_strand() {
        let paf = "q\t100\t10\t90\t-\tt\t200\t20\t180\t70\t80\t60\n";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].strand, '-');
    }

    #[test]
    fn test_parse_paf_file_multiple_records() {
        let paf = "\
q1\t1000\t0\t500\t+\tt1\t2000\t0\t500\t400\t500\t60
q2\t800\t100\t700\t-\tt2\t1500\t200\t1200\t500\t600\t60
q3\t500\t50\t450\t+\tt3\t600\t60\t460\t350\t400\t60
";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].query_name, "q1");
        assert_eq!(records[1].query_name, "q2");
        assert_eq!(records[1].strand, '-');
        assert_eq!(records[2].query_name, "q3");
        assert_eq!(records[2].target_start, 60);
        assert_eq!(records[2].target_end, 460);
    }

    #[test]
    fn test_parse_paf_file_extra_columns() {
        // PAF files often have optional tag columns beyond the 12 mandatory fields.
        let paf = "q\t1000\t10\t90\t+\tt\t2000\t20\t180\t70\t80\t60\tcg:Z:80M\tNM:i:5\n";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_name, "q");
        assert_eq!(records[0].target_name, "t");
    }

    #[test]
    fn test_parse_paf_file_skips_short_lines() {
        // Lines with fewer than 12 tab-delimited fields should be skipped.
        let paf = "\
q1\t1000\t0\t500\t+\tt1\t2000\t0\t500\t400\t500\t60
too\tfew\tfields
q2\t800\t100\t700\t+\tt2\t1500\t200\t1200\t500\t600\t60
";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].query_name, "q1");
        assert_eq!(records[1].query_name, "q2");
    }

    #[test]
    fn test_parse_paf_file_empty_file() {
        let f = write_paf_temp("");
        let records = parse_paf_file(f.path()).unwrap();
        assert!(records.is_empty());
    }

    #[test]
    fn test_parse_paf_file_blank_lines() {
        // Blank lines should be skipped (fewer than 12 fields).
        let paf = "\n\nq\t100\t0\t50\t+\tt\t200\t0\t50\t40\t50\t60\n\n";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn test_parse_paf_file_nonexistent_path() {
        let result = parse_paf_file(std::path::Path::new("/nonexistent/path.paf"));
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_paf_file_non_numeric_fields_default_to_zero() {
        // If numeric fields can't be parsed, unwrap_or(0) kicks in.
        let paf = "q\tNaN\tabc\txyz\t+\tt\tNaN\tabc\txyz\t0\t0\t0\n";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_len, 0);
        assert_eq!(records[0].query_start, 0);
        assert_eq!(records[0].query_end, 0);
        assert_eq!(records[0].target_len, 0);
        assert_eq!(records[0].target_start, 0);
        assert_eq!(records[0].target_end, 0);
    }

    #[test]
    fn test_parse_paf_file_strand_defaults_to_plus() {
        // Any strand value other than "-" should be treated as '+'.
        let paf = "q\t100\t0\t50\t*\tt\t200\t0\t50\t40\t50\t60\n";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        assert_eq!(records[0].strand, '+');
    }

    #[test]
    fn test_parse_paf_file_realize_style_names() {
        // The realize engine generates names like "chr1:100-200(+)".
        let paf = "chr1:100-200(+)\t100\t0\t100\t+\tchr2:300-400(+)\t100\t0\t100\t90\t100\t60\n";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_name, "chr1:100-200(+)");
        assert_eq!(records[0].target_name, "chr2:300-400(+)");
    }

    #[test]
    fn test_parse_paf_file_skips_self_alignments() {
        // Self-alignments (query_name == target_name) should be filtered out.
        let paf = "\
q1\t1000\t0\t1000\t+\tq1\t1000\t0\t1000\t1000\t1000\t255
q1\t1000\t0\t1000\t+\tq2\t1000\t0\t1000\t990\t1000\t255
q2\t1000\t0\t1000\t+\tq2\t1000\t0\t1000\t1000\t1000\t255
";
        let f = write_paf_temp(paf);
        let records = parse_paf_file(f.path()).unwrap();
        // Only the cross-alignment (q1→q2) should survive.
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].query_name, "q1");
        assert_eq!(records[0].target_name, "q2");
    }

    // -----------------------------------------------------------------------
    // Recursive sweepga → partition → lace integration tests
    // -----------------------------------------------------------------------

    /// Generate a pseudorandom DNA sequence of given length using a simple LCG.
    /// The seed ensures reproducibility across test runs.
    fn make_random_dna(len: usize, seed: u64) -> String {
        const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
        let mut state = seed;
        (0..len)
            .map(|_| {
                // LCG: state = (a * state + c) mod m
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                BASES[((state >> 33) % 4) as usize]
            })
            .collect()
    }

    /// Create a variant of a base sequence by introducing a SNP every `interval` bp.
    fn mutate_every_n(base: &str, interval: usize, replacement: char) -> String {
        base.chars()
            .enumerate()
            .map(|(i, c)| if i > 0 && i % interval == 0 { replacement } else { c })
            .collect()
    }

    #[test]
    fn test_recursive_path_triggers_sweepga() {
        // Sequences larger than poa_threshold trigger the sweepga path.
        // Use 5000bp pseudorandom sequences — FastGA needs sufficient length for
        // k-mer seeding to produce alignments that survive filtering.
        let base = make_random_dna(5000, 42);
        let variant = mutate_every_n(&base, 500, 'T');

        let seqs = vec![
            make_seq("s1", &base, 0, 10000),
            make_seq("s2", &variant, 0, 10000),
        ];

        let mut config = RealizeConfig::default();
        config.poa_threshold = 1000;
        config.chunk_size = 2000;
        config.padding = 0;
        config.sort_output = false;
        config.max_depth = 5;


        let result = realize_from_sequences(&seqs, &config).unwrap();

        // Must produce valid GFA with segments and paths.
        assert!(result.gfa.contains("H\tVN:Z:1.0"), "missing GFA header");
        assert!(
            result.gfa.lines().any(|l| l.starts_with("S\t")),
            "no segments in GFA"
        );
        assert!(
            result.gfa.lines().any(|l| l.starts_with("P\t")),
            "no paths in GFA"
        );

        // The sweepga path should have been taken.
        assert!(
            result.stats.sweepga_calls >= 1,
            "expected sweepga_calls >= 1, got {}",
            result.stats.sweepga_calls
        );

        // POA should also have been called (as leaf base cases).
        assert!(
            result.stats.poa_calls >= 1,
            "expected poa_calls >= 1, got {}",
            result.stats.poa_calls
        );

        // Should have recursed (depth > 0 reached).
        assert!(
            result.stats.max_depth_reached >= 1,
            "expected recursion depth >= 1, got {}",
            result.stats.max_depth_reached
        );
    }

    #[test]
    fn test_recursive_path_multi_sequence() {
        // Three similar 5kb sequences → sweepga + partition + lace.
        let base = make_random_dna(5000, 42);
        let v1 = mutate_every_n(&base, 500, 'A');
        let v2 = mutate_every_n(&base, 700, 'G');

        let seqs = vec![
            make_seq("s1", &base, 0, 10000),
            make_seq("s2", &v1, 0, 10000),
            make_seq("s3", &v2, 0, 10000),
        ];

        let mut config = RealizeConfig::default();
        config.poa_threshold = 1000;
        config.chunk_size = 2000;
        config.padding = 0;
        config.sort_output = false;
        config.max_depth = 5;


        let result = realize_from_sequences(&seqs, &config).unwrap();

        assert_eq!(result.stats.num_sequences, 3);
        assert!(result.stats.sweepga_calls >= 1);
        assert!(result.stats.poa_calls >= 1);

        // All three paths should be present in the output GFA.
        let path_lines: Vec<&str> = result.gfa.lines().filter(|l| l.starts_with("P\t")).collect();
        assert!(
            path_lines.len() >= 3,
            "expected at least 3 paths, got {}: {:?}",
            path_lines.len(),
            path_lines
        );
    }

    #[test]
    fn test_recursive_max_depth_forces_poa() {
        // With max_depth=1, depth 0 uses sweepga, chunks at depth 1 must use POA
        // even though they exceed poa_threshold.
        let base = make_random_dna(5000, 42);
        let variant = mutate_every_n(&base, 500, 'T');

        let seqs = vec![
            make_seq("s1", &base, 0, 10000),
            make_seq("s2", &variant, 0, 10000),
        ];

        let mut config = RealizeConfig::default();
        config.poa_threshold = 1000;
        config.chunk_size = 2000;   // Chunks ~2000bp > poa_threshold.
        config.padding = 0;
        config.sort_output = false;
        config.max_depth = 1;       // Force POA at depth 1.


        let result = realize_from_sequences(&seqs, &config).unwrap();

        assert!(result.gfa.contains("H\tVN:Z:1.0"));
        assert!(
            result.gfa.lines().any(|l| l.starts_with("S\t")),
            "no segments"
        );

        // Depth 0: sweepga should be called.
        assert!(
            result.stats.sweepga_calls >= 1,
            "expected sweepga at depth 0"
        );
        // Depth 1: forced to POA (max_depth reached).
        assert!(
            result.stats.poa_calls >= 1,
            "expected POA calls from max_depth forcing"
        );
        // max_depth_reached should be at least 1 (chunks recurse to depth 1).
        assert!(
            result.stats.max_depth_reached >= 1,
            "expected max_depth_reached >= 1, got {}",
            result.stats.max_depth_reached
        );
    }

    #[test]
    fn test_recursive_stats_depth_tracking() {
        // Verify that depth tracking works correctly in the recursive case.
        let base = make_random_dna(5000, 42);
        let variant = mutate_every_n(&base, 500, 'T');

        let seqs = vec![
            make_seq("s1", &base, 0, 10000),
            make_seq("s2", &variant, 0, 10000),
        ];

        let mut config = RealizeConfig::default();
        config.poa_threshold = 1000;
        config.chunk_size = 2000;
        config.padding = 0;
        config.sort_output = false;
        config.max_depth = 10;


        let result = realize_from_sequences(&seqs, &config).unwrap();

        // Should have recursed at least once (depth > 0).
        assert!(
            result.stats.max_depth_reached >= 1,
            "expected recursion depth >= 1, got {}",
            result.stats.max_depth_reached
        );

        // Both sweepga and POA should have been called.
        assert!(result.stats.sweepga_calls >= 1);
        assert!(result.stats.poa_calls >= 1);
    }

    #[test]
    fn test_recursive_gfa_paths_cover_sequences() {
        // Verify that all input sequences appear as paths in the output GFA
        // and every path step references a valid segment.
        let base = make_random_dna(5000, 42);
        let v1 = mutate_every_n(&base, 500, 'A');

        let seqs = vec![
            make_seq("s1", &base, 0, 10000),
            make_seq("s2", &v1, 50, 10000),
        ];

        let mut config = RealizeConfig::default();
        config.poa_threshold = 1000;
        config.chunk_size = 2000;
        config.padding = 0;
        config.sort_output = false;
        config.max_depth = 5;


        let result = realize_from_sequences(&seqs, &config).unwrap();

        // Both sequences should have corresponding paths.
        let path_lines: Vec<&str> = result.gfa.lines().filter(|l| l.starts_with("P\t")).collect();
        assert!(
            path_lines.len() >= 2,
            "expected at least 2 paths, got {}: {:?}",
            path_lines.len(),
            path_lines
        );

        // Every segment referenced by a path should exist.
        let segment_ids: std::collections::HashSet<&str> = result
            .gfa
            .lines()
            .filter(|l| l.starts_with("S\t"))
            .filter_map(|l| l.split('\t').nth(1))
            .collect();

        for path_line in &path_lines {
            let steps_str = path_line.split('\t').nth(2).unwrap_or("");
            for step in steps_str.split(',') {
                let node = step.trim_end_matches(['+', '-']);
                assert!(
                    segment_ids.contains(node),
                    "path references non-existent segment '{}' in: {}",
                    node,
                    path_line
                );
            }
        }
    }

    #[test]
    fn test_recursive_links_reference_valid_segments() {
        // Verify that all links reference segments that exist.
        let base = make_random_dna(5000, 42);
        let variant = mutate_every_n(&base, 500, 'G');

        let seqs = vec![
            make_seq("s1", &base, 0, 10000),
            make_seq("s2", &variant, 0, 10000),
        ];

        let mut config = RealizeConfig::default();
        config.poa_threshold = 1000;
        config.chunk_size = 2000;
        config.padding = 0;
        config.sort_output = false;
        config.max_depth = 5;


        let result = realize_from_sequences(&seqs, &config).unwrap();

        let segment_ids: std::collections::HashSet<&str> = result
            .gfa
            .lines()
            .filter(|l| l.starts_with("S\t"))
            .filter_map(|l| l.split('\t').nth(1))
            .collect();

        for line in result.gfa.lines().filter(|l| l.starts_with("L\t")) {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 5 {
                assert!(
                    segment_ids.contains(fields[1]),
                    "link from-node '{}' not in segments: {}",
                    fields[1],
                    line
                );
                assert!(
                    segment_ids.contains(fields[3]),
                    "link to-node '{}' not in segments: {}",
                    fields[3],
                    line
                );
            }
        }
    }

    #[test]
    fn test_sweepga_empty_alignment_falls_back_to_poa() {
        // When sequences are too small/repetitive for sweepga to produce
        // usable alignments, the engine should gracefully fall back to POA.
        let base = make_random_dna(500, 42);
        let variant = mutate_every_n(&base, 50, 'T');

        let seqs = vec![
            make_seq("s1", &base, 0, 1000),
            make_seq("s2", &variant, 0, 1000),
        ];

        let mut config = RealizeConfig::default();
        config.poa_threshold = 100;  // Below sequence length → tries sweepga.
        config.chunk_size = 200;
        config.padding = 0;
        config.sort_output = false;
        config.max_depth = 1; // One sweepga level, then POA — avoids fragile inner recursion on CI.


        let result = realize_from_sequences(&seqs, &config).unwrap();

        // Should produce valid GFA regardless of whether sweepga succeeded.
        assert!(result.gfa.contains("H\tVN:Z:1.0"));
        assert!(result.gfa.lines().any(|l| l.starts_with("S\t")));
        assert!(result.gfa.lines().any(|l| l.starts_with("P\t")));
        assert!(result.stats.poa_calls >= 1);
    }

    #[test]
    fn test_partition_multiple_chunks_with_projections() {
        // Test partition_into_chunks with synthetic PAF records.
        // Anchor is 1000bp, chunk_size=300, padding=50.
        let anchor = make_random_dna(1000, 42);
        let other = make_random_dna(900, 99);

        let seqs = vec![
            make_seq("anchor", &anchor, 0, 2000),
            make_seq("other", &other, 0, 2000),
        ];

        // PAF: anchor aligns to other over [0,900) → [0,900).
        let paf_records = vec![PafRecord {
            query_name: "anchor:0-1000(+)".to_string(),
            query_len: 1000,
            query_start: 0,
            query_end: 900,
            target_name: "other:0-900(+)".to_string(),
            target_len: 900,
            target_start: 0,
            target_end: 900,
            strand: '+',
        }];

        let chunks = partition_into_chunks(&seqs, &paf_records, 300, 50);

        // 1000bp / 300bp chunks = ~4 windows.
        assert!(
            chunks.len() >= 3,
            "expected at least 3 chunks, got {}",
            chunks.len()
        );

        // Each chunk should contain the anchor slice.
        for chunk in &chunks {
            assert!(
                !chunk.sequences.is_empty(),
                "chunk should have at least the anchor"
            );
        }

        // At least some chunks should contain the other sequence (projected).
        let chunks_with_other: usize = chunks
            .iter()
            .filter(|c| c.sequences.len() > 1)
            .count();
        assert!(
            chunks_with_other >= 2,
            "expected multiple chunks to contain the projected sequence, got {}",
            chunks_with_other
        );
    }

    // -----------------------------------------------------------------------
    // P1: Multi-chunk partition and projection edge cases
    // -----------------------------------------------------------------------

    // --- project_window_to_sequence edge cases ---

    #[test]
    fn test_project_window_empty_projections() {
        assert_eq!(project_window_to_sequence(&[], 0, 100), None);
        assert_eq!(project_window_to_sequence(&[], 500, 1000), None);
    }

    #[test]
    fn test_project_window_gap_between_disjoint_blocks() {
        // Window falls entirely in a gap between two alignment blocks → None.
        let projections = vec![
            (0, 100, 0, 100, false),
            (200, 300, 200, 300, false),
        ];
        assert_eq!(project_window_to_sequence(&projections, 100, 200), None);
    }

    #[test]
    fn test_project_window_spans_disjoint_blocks() {
        // Window overlaps both disjoint blocks.
        let projections = vec![
            (0, 100, 0, 100, false),
            (200, 300, 200, 300, false),
        ];
        let result = project_window_to_sequence(&projections, 90, 210);
        // Block 1: [90,100) → seq [90,100)
        // Block 2: [200,210) → seq [200,210)
        assert_eq!(result, Some((90, 210)));
    }

    #[test]
    fn test_project_window_compression_scaling() {
        // 2:1 compression: anchor [0,100) → seq [0,50).
        // Window [25,75): frac 0.25..0.75 → seq [12,37)
        let projections = vec![(0, 100, 0, 50, false)];
        assert_eq!(project_window_to_sequence(&projections, 25, 75), Some((12, 37)));
    }

    #[test]
    fn test_project_window_expansion_scaling() {
        // 1:4 expansion: anchor [0,50) → seq [0,200).
        // Window [0,25): frac 0.0..0.5 → seq [0,100)
        let projections = vec![(0, 50, 0, 200, false)];
        assert_eq!(project_window_to_sequence(&projections, 0, 25), Some((0, 100)));
    }

    #[test]
    fn test_project_window_reverse_compression() {
        // Reverse strand with 2:1 compression.
        // anchor [0,100) → seq [0,50), reversed.
        // Window [0,50): frac_start=0.0, frac_end=0.5
        // p_end = 50 - 0*50 = 50, p_start = 50 - 25 = 25
        let projections = vec![(0, 100, 0, 50, true)];
        assert_eq!(project_window_to_sequence(&projections, 0, 50), Some((25, 50)));
    }

    #[test]
    fn test_project_window_touching_adjacent_blocks() {
        // Two blocks touch at position 100 with no gap.
        let projections = vec![
            (0, 100, 0, 100, false),
            (100, 200, 100, 200, false),
        ];
        assert_eq!(project_window_to_sequence(&projections, 80, 120), Some((80, 120)));
    }

    #[test]
    fn test_project_window_single_base() {
        let projections = vec![(0, 100, 0, 100, false)];
        assert_eq!(project_window_to_sequence(&projections, 50, 51), Some((50, 51)));
    }

    #[test]
    fn test_project_window_at_start_boundary() {
        let projections = vec![(100, 200, 50, 150, false)];
        assert_eq!(project_window_to_sequence(&projections, 100, 120), Some((50, 70)));
    }

    #[test]
    fn test_project_window_at_end_boundary() {
        let projections = vec![(100, 200, 50, 150, false)];
        assert_eq!(project_window_to_sequence(&projections, 180, 200), Some((130, 150)));
    }

    // --- partition_into_chunks edge cases ---

    #[test]
    fn test_partition_empty_input() {
        assert!(partition_into_chunks(&[], &[], 1000, 100).is_empty());
    }

    #[test]
    fn test_partition_multi_chunk_projected_slices() {
        // Anchor 300bp, chunk_size=100, padding=10 → 3 windows.
        // Verify anchor ranges and that s2 appears in each chunk.
        let anchor = "A".repeat(300);
        let seq2 = "C".repeat(280);
        let seqs = vec![
            make_seq("anchor", &anchor, 0, 300),
            make_seq("s2", &seq2, 0, 280),
        ];

        let anchor_name = seqs[0].1.alignment_name();
        let s2_name = seqs[1].1.alignment_name();
        let paf_records = vec![PafRecord {
            query_name: anchor_name,
            query_len: 300,
            query_start: 0,
            query_end: 300,
            target_name: s2_name,
            target_len: 280,
            target_start: 0,
            target_end: 280,
            strand: '+',
        }];

        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 10);
        assert_eq!(chunks.len(), 3);

        for chunk in &chunks {
            assert_eq!(chunk.sequences.len(), 2, "each chunk: anchor + s2");
        }

        // Windows: pos=0→[0,110), pos=100→[90,210), pos=200→[190,300)
        assert_eq!((chunks[0].anchor_start, chunks[0].anchor_end), (0, 110));
        assert_eq!((chunks[1].anchor_start, chunks[1].anchor_end), (90, 210));
        assert_eq!((chunks[2].anchor_start, chunks[2].anchor_end), (190, 300));
    }

    #[test]
    fn test_partition_unaligned_seq_included_as_fallback() {
        // s3 has no alignment to anchor → included fully in every chunk
        // as a fallback so it doesn't get dropped.
        let anchor = "A".repeat(300);
        let seqs = vec![
            make_seq("anchor", &anchor, 0, 300),
            make_seq("s2", &"C".repeat(280), 0, 280),
            make_seq("s3", &"G".repeat(100), 0, 100),
        ];

        let anchor_name = seqs[0].1.alignment_name();
        let s2_name = seqs[1].1.alignment_name();
        let paf_records = vec![PafRecord {
            query_name: anchor_name,
            query_len: 300,
            query_start: 0,
            query_end: 300,
            target_name: s2_name,
            target_len: 280,
            target_start: 0,
            target_end: 280,
            strand: '+',
        }];

        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 10);
        for (i, chunk) in chunks.iter().enumerate() {
            assert!(chunk.sequences.iter().any(|(_, m)| m.name == "s3"),
                "chunk {} should have s3 (fallback inclusion)", i);
        }
    }

    #[test]
    fn test_partition_partial_alignment_drops_from_late_chunks() {
        // s2 aligns only to anchor [0,150). Chunks past that range exclude s2.
        let anchor = "A".repeat(300);
        let seqs = vec![
            make_seq("anchor", &anchor, 0, 300),
            make_seq("s2", &"C".repeat(150), 0, 150),
        ];

        let anchor_name = seqs[0].1.alignment_name();
        let s2_name = seqs[1].1.alignment_name();
        let paf_records = vec![PafRecord {
            query_name: anchor_name,
            query_len: 300,
            query_start: 0,
            query_end: 150,
            target_name: s2_name,
            target_len: 150,
            target_start: 0,
            target_end: 150,
            strand: '+',
        }];

        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 10);
        assert!(chunks.len() >= 3);

        assert!(chunks[0].sequences.iter().any(|(_, m)| m.name == "s2"),
            "first chunk should include s2");
        // Last chunk [190,300) doesn't overlap [0,150).
        assert!(!chunks.last().unwrap().sequences.iter().any(|(_, m)| m.name == "s2"),
            "last chunk should not include s2");
    }

    #[test]
    fn test_partition_reverse_strand() {
        // s2 aligns in reverse across the full anchor.
        let anchor = "A".repeat(300);
        let seqs = vec![
            make_seq("anchor", &anchor, 0, 300),
            make_seq("s2", &"C".repeat(300), 0, 300),
        ];

        let anchor_name = seqs[0].1.alignment_name();
        let s2_name = seqs[1].1.alignment_name();
        let paf_records = vec![PafRecord {
            query_name: anchor_name,
            query_len: 300,
            query_start: 0,
            query_end: 300,
            target_name: s2_name,
            target_len: 300,
            target_start: 0,
            target_end: 300,
            strand: '-',
        }];

        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 10);
        assert!(chunks.len() >= 3);
        for (i, chunk) in chunks.iter().enumerate() {
            assert!(chunk.sequences.iter().any(|(_, m)| m.name == "s2"),
                "chunk {} should have s2 (reverse strand full coverage)", i);
        }
    }

    #[test]
    fn test_partition_anchor_as_paf_target() {
        // Anchor appears as target (not query) in PAF → still projects correctly.
        let anchor = "A".repeat(300);
        let seqs = vec![
            make_seq("anchor", &anchor, 0, 300),
            make_seq("s2", &"C".repeat(300), 0, 300),
        ];

        let anchor_name = seqs[0].1.alignment_name();
        let s2_name = seqs[1].1.alignment_name();
        let paf_records = vec![PafRecord {
            query_name: s2_name,
            query_len: 300,
            query_start: 0,
            query_end: 300,
            target_name: anchor_name,
            target_len: 300,
            target_start: 0,
            target_end: 300,
            strand: '+',
        }];

        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 10);
        assert!(chunks.len() >= 3);
        for (i, chunk) in chunks.iter().enumerate() {
            assert!(chunk.sequences.iter().any(|(_, m)| m.name == "s2"),
                "chunk {} should have s2 (anchor as PAF target)", i);
        }
    }

    #[test]
    fn test_partition_fragmented_alignment_gap() {
        // Two alignment blocks with a gap: anchor [0,100)→s2 [0,100), anchor [200,300)→s2 [100,200).
        // Middle chunk (anchor [100,200)) has no alignment → no s2.
        let anchor = "A".repeat(300);
        let seqs = vec![
            make_seq("anchor", &anchor, 0, 300),
            make_seq("s2", &"C".repeat(200), 0, 200),
        ];

        let anchor_name = seqs[0].1.alignment_name();
        let s2_name = seqs[1].1.alignment_name();
        let paf_records = vec![
            PafRecord {
                query_name: anchor_name.clone(), query_len: 300,
                query_start: 0, query_end: 100,
                target_name: s2_name.clone(), target_len: 200,
                target_start: 0, target_end: 100, strand: '+',
            },
            PafRecord {
                query_name: anchor_name, query_len: 300,
                query_start: 200, query_end: 300,
                target_name: s2_name, target_len: 200,
                target_start: 100, target_end: 200, strand: '+',
            },
        ];

        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 0);
        assert_eq!(chunks.len(), 3);

        assert!(chunks[0].sequences.iter().any(|(_, m)| m.name == "s2"),
            "first chunk should have s2");
        assert!(!chunks[1].sequences.iter().any(|(_, m)| m.name == "s2"),
            "middle chunk should NOT have s2 (alignment gap)");
        assert!(chunks[2].sequences.iter().any(|(_, m)| m.name == "s2"),
            "last chunk should have s2");
    }

    #[test]
    fn test_partition_length_sort_with_expansion() {
        // Projected s2 is longer than anchor slice → s2 should sort first.
        let anchor = "A".repeat(300);
        let seqs = vec![
            make_seq("anchor", &anchor, 0, 300),
            make_seq("s2", &"C".repeat(400), 0, 400),
        ];

        let anchor_name = seqs[0].1.alignment_name();
        let s2_name = seqs[1].1.alignment_name();
        let paf_records = vec![PafRecord {
            query_name: anchor_name, query_len: 300,
            query_start: 0, query_end: 300,
            target_name: s2_name, target_len: 400,
            target_start: 0, target_end: 400, strand: '+',
        }];

        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 10);
        for (i, chunk) in chunks.iter().enumerate() {
            let lens: Vec<usize> = chunk.sequences.iter().map(|(s, _)| s.len()).collect();
            for w in lens.windows(2) {
                assert!(w[0] >= w[1], "chunk {} not sorted desc: {:?}", i, lens);
            }
        }
    }

    #[test]
    fn test_partition_metadata_preserved_across_chunks() {
        // Verify that chunk metadata preserves original coordinates for
        // consistent path naming across chunks (needed for lacing).
        let anchor = "A".repeat(300);
        let seqs = vec![
            make_seq("anchor", &anchor, 100, 1000),
            make_seq("s2", &"C".repeat(300), 50, 500),
        ];

        let anchor_name = seqs[0].1.alignment_name();
        let s2_name = seqs[1].1.alignment_name();
        let paf_records = vec![PafRecord {
            query_name: anchor_name, query_len: 300,
            query_start: 0, query_end: 300,
            target_name: s2_name, target_len: 300,
            target_start: 0, target_end: 300, strand: '+',
        }];

        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 0);
        assert!(chunks.len() >= 3);

        // All chunks should preserve the original metadata so path_name()
        // is consistent for lacing.
        for (i, chunk) in chunks.iter().enumerate() {
            let anchor_meta = &chunk.sequences.iter().find(|(_, m)| m.name == "anchor").unwrap().1;
            assert_eq!(anchor_meta.start, 100,
                "chunk {} anchor should preserve original start", i);
            assert_eq!(anchor_meta.size, 300,
                "chunk {} anchor should preserve original size", i);
        }
    }

    #[test]
    fn test_partition_no_paf_includes_unaligned_as_fallback() {
        // Large anchor, no PAF records → each chunk has the anchor plus
        // unaligned sequences as fallback (to avoid dropping them).
        let seqs = vec![
            make_seq("anchor", &"A".repeat(300), 0, 300),
            make_seq("s2", &"C".repeat(100), 0, 100),
        ];

        let chunks = partition_into_chunks(&seqs, &[], 100, 10);
        assert!(chunks.len() >= 3);
        for (i, chunk) in chunks.iter().enumerate() {
            assert_eq!(chunk.sequences.len(), 2,
                "chunk {} should have anchor + s2 (fallback)", i);
        }
    }

    #[test]
    fn test_partition_five_sequences_all_aligned() {
        // All 5 sequences aligned → every chunk has 5 sequences.
        let anchor = "A".repeat(500);
        let seqs: Vec<_> = std::iter::once(make_seq("anchor", &anchor, 0, 500))
            .chain((1..=4).map(|i| make_seq(&format!("s{}", i), &"C".repeat(500), 0, 500)))
            .collect();

        let anchor_name = seqs[0].1.alignment_name();
        let paf_records: Vec<PafRecord> = (1..=4).map(|i| PafRecord {
            query_name: anchor_name.clone(), query_len: 500,
            query_start: 0, query_end: 500,
            target_name: seqs[i].1.alignment_name(), target_len: 500,
            target_start: 0, target_end: 500, strand: '+',
        }).collect();

        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 10);
        for (i, chunk) in chunks.iter().enumerate() {
            assert_eq!(chunk.sequences.len(), 5, "chunk {} should have all 5 seqs", i);
        }
    }

    #[test]
    fn test_partition_projection_clamped_to_seq_length() {
        // PAF says s2 maps to [0,300) but s2 is only 200bp. Should not panic.
        let seqs = vec![
            make_seq("anchor", &"A".repeat(300), 0, 300),
            make_seq("s2", &"C".repeat(200), 0, 200),
        ];

        let anchor_name = seqs[0].1.alignment_name();
        let s2_name = seqs[1].1.alignment_name();
        let paf_records = vec![PafRecord {
            query_name: anchor_name, query_len: 300,
            query_start: 0, query_end: 300,
            target_name: s2_name, target_len: 300,
            target_start: 0, target_end: 300, strand: '+',
        }];

        // Should not panic.
        let chunks = partition_into_chunks(&seqs, &paf_records, 100, 0);
        for chunk in &chunks {
            for (seq_str, _) in &chunk.sequences {
                assert!(!seq_str.is_empty());
            }
        }
    }

    #[test]
    fn test_partition_window_padding_clamped_at_boundaries() {
        // First window starts at 0 (not negative), last window ends at anchor_len.
        let seqs = vec![make_seq("anchor", &"A".repeat(250), 0, 250)];
        let chunks = partition_into_chunks(&seqs, &[], 100, 20);

        assert_eq!(chunks[0].anchor_start, 0, "first window should start at 0");
        assert_eq!(chunks.last().unwrap().anchor_end, 250, "last window should end at anchor_len");
    }

    #[test]
    fn test_partition_boundary_exactly_threshold() {
        // anchor_len == chunk_size + padding → single chunk (early return).
        let seqs = vec![
            make_seq("anchor", &"A".repeat(110), 0, 110),
            make_seq("s2", &"C".repeat(100), 0, 100),
        ];
        let chunks = partition_into_chunks(&seqs, &[], 100, 10);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].sequences.len(), 2);
    }

    #[test]
    fn test_partition_boundary_one_over_threshold() {
        // anchor_len = chunk_size + padding + 1 → multi-chunk.
        let seqs = vec![make_seq("anchor", &"A".repeat(111), 0, 111)];
        let chunks = partition_into_chunks(&seqs, &[], 100, 10);
        assert!(chunks.len() >= 2, "should produce multiple chunks");
    }

    #[test]
    fn test_project_window_zero_length_alignment_span() {
        // a_start == a_end → zero-length alignment block.
        // The code uses .max(1) on a_span to avoid division by zero.
        // With s_start == s_end, projection produces a zero-length range → None.
        let projections = vec![(50, 50, 100, 100, false)];
        assert_eq!(project_window_to_sequence(&projections, 0, 100), None);
    }

    #[test]
    fn test_project_window_zero_length_at_boundary() {
        // Zero-length alignment at the exact window boundary.
        let projections = vec![(100, 100, 50, 50, false)];
        // a_start(100) >= win_end(100) → skip. Returns None.
        assert_eq!(project_window_to_sequence(&projections, 0, 100), None);
    }

    #[test]
    fn test_partition_small_chunk_size_one() {
        // chunk_size = 1 should create many chunks without panicking.
        let anchor = "A".repeat(50);
        let other = "C".repeat(50);

        let seqs = vec![
            make_seq("anchor", &anchor, 0, 100),
            make_seq("other", &other, 0, 100),
        ];

        let paf_records = vec![PafRecord {
            query_name: seqs[0].1.alignment_name(),
            query_len: 50,
            query_start: 0,
            query_end: 50,
            target_name: seqs[1].1.alignment_name(),
            target_len: 50,
            target_start: 0,
            target_end: 50,
            strand: '+',
        }];

        // chunk_size=1, padding=0 → should create ~50 chunks.
        let chunks = partition_into_chunks(&seqs, &paf_records, 1, 0);
        assert!(chunks.len() >= 40, "expected many chunks with chunk_size=1, got {}", chunks.len());

        // No chunk should be empty (anchor always included).
        for chunk in &chunks {
            assert!(!chunk.sequences.is_empty());
        }
    }

    #[test]
    fn test_partition_chunk_size_larger_than_anchor() {
        // chunk_size >> anchor_len → single chunk (early return).
        let seqs = vec![
            make_seq("anchor", &"A".repeat(100), 0, 200),
            make_seq("s2", &"C".repeat(80), 0, 200),
        ];
        let chunks = partition_into_chunks(&seqs, &[], 10000, 0);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].sequences.len(), 2);
    }

    #[test]
    fn test_partition_large_padding_relative_to_chunk() {
        // padding = chunk_size → heavy overlap between adjacent chunks.
        let anchor = "A".repeat(500);
        let seqs = vec![make_seq("anchor", &anchor, 0, 500)];

        let chunks = partition_into_chunks(&seqs, &[], 100, 100);
        assert!(chunks.len() >= 2);

        // Every pair of adjacent chunks should overlap significantly.
        for pair in chunks.windows(2) {
            let overlap = pair[0].anchor_end.saturating_sub(pair[1].anchor_start);
            assert!(overlap >= 100, "expected >= 100bp overlap, got {}", overlap);
        }
    }

    #[test]
    fn test_project_window_many_small_blocks() {
        // Many small alignment blocks covering the window.
        let projections: Vec<_> = (0..10)
            .map(|i| {
                let start = i * 20;
                (start, start + 20, start, start + 20, false)
            })
            .collect();

        // Window covering all blocks [0, 200).
        let result = project_window_to_sequence(&projections, 0, 200);
        assert_eq!(result, Some((0, 200)));

        // Window covering middle blocks [40, 120).
        let result = project_window_to_sequence(&projections, 40, 120);
        assert_eq!(result, Some((40, 120)));
    }

    #[test]
    fn test_project_window_overlapping_reverse_forward() {
        // A forward and reverse block cover the same anchor region.
        // Forward: anchor [0,100) → seq [0,100)
        // Reverse: anchor [0,100) → seq [200,300)
        // Window [0,100):
        //   Forward: [0,100)
        //   Reverse: p_end = 300 - 0 = 300, p_start = 300 - 100 = 200 → [200,300)
        //   Union: [0, 300)
        let projections = vec![
            (0, 100, 0, 100, false),
            (0, 100, 200, 300, true),
        ];
        let result = project_window_to_sequence(&projections, 0, 100);
        assert_eq!(result, Some((0, 300)));
    }
}
