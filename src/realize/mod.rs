pub mod poa;

use crate::commands::align::{SparsificationStrategy, SweepgaAlignConfig};
use crate::graph::{
    prepare_sequences, sort_gfa, SeqwishConfig, SequenceMetadata,
};
use crate::impg_index::ImpgIndex;
use crate::realize::poa::padded_poa_from_sequences;
use crate::sequence_index::UnifiedSequenceIndex;
use coitrees::Interval;
use log::{info, warn};
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

    /// Sparsification strategy per recursion depth level.
    /// Index 0 = depth 0 (coarsest), higher indices = deeper levels.
    /// If depth exceeds the length, the last entry is reused.
    /// Default: [Connectivity(0.99)]
    pub sparsification: Vec<SparsificationStrategy>,

    /// SPOA scoring parameters: (match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2)
    /// Default: (5, 4, 6, 2, 24, 1)
    pub scoring_params: (u8, u8, u8, u8, u8, u8),

    /// SeqwishConfig for sweepga+seqwish alignment steps.
    pub seqwish_config: SeqwishConfig,

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
            sparsification: vec![SparsificationStrategy::Connectivity(0.99)],
            scoring_params: (5, 4, 6, 2, 24, 1),
            seqwish_config: SeqwishConfig::default(),
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

    let _strategy = config
        .sparsification
        .get(depth)
        .or(config.sparsification.last())
        .cloned()
        .unwrap_or(SparsificationStrategy::Connectivity(0.99));

    // Run sweepga all-vs-all alignment.
    let align_config = build_sweepga_config(config);
    let named_seqs: Vec<(String, &[u8])> = sequences
        .iter()
        .map(|(seq, meta)| {
            let name = format!(
                "{}:{}-{}({})",
                meta.name,
                meta.start,
                meta.start + meta.size,
                meta.strand
            );
            (name, seq.as_bytes())
        })
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

    // Recurse on each chunk.
    let mut sub_gfas: Vec<String> = Vec::with_capacity(chunks.len());
    for chunk in &chunks {
        if chunk.sequences.is_empty() {
            continue;
        }
        let sub_gfa = realize_recursive(
            &chunk.sequences,
            config,
            depth + 1,
            poa_calls,
            sweepga_calls,
            max_depth_reached,
        )?;
        sub_gfas.push(sub_gfa);
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
fn build_sweepga_config(config: &RealizeConfig) -> SweepgaAlignConfig {
    SweepgaAlignConfig {
        num_threads: config.num_threads,
        kmer_frequency: 10, // Reasonable default for small sets
        min_alignment_length: 100,
        no_filter: false,
        num_mappings: "1:1".to_string(),
        scaffold_jump: 50_000,
        scaffold_mass: 10_000,
        scaffold_filter: "1:1".to_string(),
        overlap: 0.95,
        min_identity: 0.0,
        scaffold_dist: 0,
        min_mapping_length: 0,
        temp_dir: config.temp_dir.clone(),
    }
}

/// Parse a PAF file into PafRecord structs.
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
        let query_len: usize = fields[1].parse().unwrap_or(0);
        let query_start: usize = fields[2].parse().unwrap_or(0);
        let query_end: usize = fields[3].parse().unwrap_or(0);
        let strand = if fields[4] == "-" { '-' } else { '+' };
        let target_name = fields[5].to_string();
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
    let anchor_name = format!(
        "{}:{}-{}({})",
        anchor_meta.name,
        anchor_meta.start,
        anchor_meta.start + anchor_meta.size,
        anchor_meta.strand
    );

    // If the anchor is smaller than one chunk, return a single chunk with all sequences.
    if anchor_len <= chunk_size + padding {
        return vec![Chunk {
            sequences: sequences.to_vec(),
            anchor_start: 0,
            anchor_end: anchor_len,
        }];
    }

    // Build name-to-index map for the input sequences.
    let seq_name_map: std::collections::HashMap<String, usize> = sequences
        .iter()
        .enumerate()
        .map(|(i, (_, meta))| {
            let name = format!(
                "{}:{}-{}({})",
                meta.name, meta.start, meta.start + meta.size, meta.strand
            );
            (name, i)
        })
        .collect();

    // For each non-anchor sequence, collect alignment projections from the anchor.
    // projection[i] = Vec<(anchor_start, anchor_end, seq_start, seq_end)>
    let mut projections: Vec<Vec<(usize, usize, usize, usize)>> =
        vec![Vec::new(); sequences.len()];

    for rec in paf_records {
        // We want alignments where the anchor is either query or target.
        if rec.query_name == anchor_name {
            if let Some(&idx) = seq_name_map.get(&rec.target_name) {
                projections[idx].push((
                    rec.query_start,
                    rec.query_end,
                    rec.target_start,
                    rec.target_end,
                ));
            }
        } else if rec.target_name == anchor_name {
            if let Some(&idx) = seq_name_map.get(&rec.query_name) {
                projections[idx].push((
                    rec.target_start,
                    rec.target_end,
                    rec.query_start,
                    rec.query_end,
                ));
            }
        }
    }

    // Sort projections by anchor start for binary search.
    for proj in &mut projections {
        proj.sort_by_key(|&(a_start, _, _, _)| a_start);
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
        let anchor_slice = &anchor_seq[*win_start..*win_end];
        let mut chunk_anchor_meta = anchor_meta.clone();
        chunk_anchor_meta.start = anchor_meta.start + *win_start as i32;
        chunk_anchor_meta.size = (*win_end - *win_start) as i32;
        chunk_seqs.push((anchor_slice.to_string(), chunk_anchor_meta));

        // For each non-anchor sequence, project the window through alignments.
        for (seq_idx, (seq_str, seq_meta)) in sequences.iter().enumerate().skip(1) {
            if let Some((seq_start, seq_end)) =
                project_window_to_sequence(&projections[seq_idx], *win_start, *win_end)
            {
                let clamped_start = seq_start.min(seq_str.len());
                let clamped_end = seq_end.min(seq_str.len());
                if clamped_end > clamped_start {
                    let sub_seq = &seq_str[clamped_start..clamped_end];
                    let mut sub_meta = seq_meta.clone();
                    sub_meta.start = seq_meta.start + clamped_start as i32;
                    sub_meta.size = (clamped_end - clamped_start) as i32;
                    chunk_seqs.push((sub_seq.to_string(), sub_meta));
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
fn project_window_to_sequence(
    projections: &[(usize, usize, usize, usize)],
    win_start: usize,
    win_end: usize,
) -> Option<(usize, usize)> {
    let mut best_start = usize::MAX;
    let mut best_end = 0usize;

    for &(a_start, a_end, s_start, s_end) in projections {
        // Check if this alignment overlaps the window.
        if a_start >= win_end || a_end <= win_start {
            continue;
        }

        // Linear interpolation: map [win_start, win_end) through this alignment block.
        let a_span = a_end.saturating_sub(a_start).max(1);
        let s_span = s_end.saturating_sub(s_start);

        let overlap_start = win_start.max(a_start);
        let overlap_end = win_end.min(a_end);

        // Project the overlap region.
        let frac_start = (overlap_start - a_start) as f64 / a_span as f64;
        let frac_end = (overlap_end - a_start) as f64 / a_span as f64;

        let projected_start = s_start + (frac_start * s_span as f64) as usize;
        let projected_end = s_start + (frac_end * s_span as f64) as usize;

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

    // Deduplicate links.
    let link_set: std::collections::HashSet<&String> = all_links.iter().collect();
    for link in &link_set {
        output.push_str(link);
        output.push('\n');
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
    let path_name = format!(
        "{}:{}-{}",
        meta.name,
        if meta.strand == '+' {
            meta.start
        } else {
            (meta.total_length as i32) - meta.start - meta.size
        },
        if meta.strand == '+' {
            meta.start + meta.size
        } else {
            (meta.total_length as i32) - meta.start
        }
    );
    Ok(format!(
        "H\tVN:Z:1.0\nS\t1\t{}\nP\t{}\t1+\t*\n",
        seq, path_name
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
        let projections = vec![(100, 200, 50, 150)];
        // Window completely before the alignment.
        assert_eq!(project_window_to_sequence(&projections, 0, 50), None);
        // Window completely after the alignment.
        assert_eq!(project_window_to_sequence(&projections, 300, 400), None);
    }

    #[test]
    fn test_project_window_to_sequence_full_overlap() {
        // Alignment: anchor [100,200) → seq [50,150)
        let projections = vec![(100, 200, 50, 150)];
        // Window exactly covers the alignment.
        let result = project_window_to_sequence(&projections, 100, 200);
        assert_eq!(result, Some((50, 150)));
    }

    #[test]
    fn test_project_window_to_sequence_partial_overlap() {
        // Alignment: anchor [100,200) → seq [0,100)
        let projections = vec![(100, 200, 0, 100)];
        // Window covers first half: anchor [100,150) → seq [0,50)
        let result = project_window_to_sequence(&projections, 100, 150);
        assert_eq!(result, Some((0, 50)));
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
}
