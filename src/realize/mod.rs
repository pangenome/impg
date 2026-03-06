pub mod poa;

use crate::commands::align::{SparsificationStrategy, SweepgaAlignConfig};
use crate::commands::lace::{
    CompactGraph, RangeInfo, link_contiguous_ranges, mark_nodes_for_removal, sort_and_filter_ranges,
    split_path_name, trim_range_overlaps,
};
use crate::graph::{
    prepare_sequences, sort_gfa, SequenceMetadata,
};
use crate::impg_index::ImpgIndex;
use crate::realize::poa::padded_poa_from_sequences;
use crate::sequence_index::UnifiedSequenceIndex;
use coitrees::Interval;
use handlegraph::handle::{Handle, NodeId};
use log::{info, warn};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::io::{self, BufRead, BufReader};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
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

    /// When the number of sequences exceeds this threshold, use seqwish
    /// instead of POA for base-case graph construction. POA (SPOA) is
    /// O(N×L) in memory and becomes impractical for hundreds of sequences.
    /// Default: 500
    pub seqwish_threshold: usize,

    /// Optional directory for saving intermediate debug files (PAFs, sub-GFAs, FASTAs).
    /// When set, each recursive step saves its inputs and outputs.
    pub debug_dir: Option<String>,

    /// Aligner backend: "wfmash" or "fastga"
    pub aligner: String,
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
            seqwish_threshold: 500,
            debug_dir: None,
            aligner: "fastga".to_string(),
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
    /// Number of seqwish leaf calls (used when sequence count exceeds seqwish_threshold).
    pub seqwish_calls: usize,
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
    let sequences = prepare_sequences(impg, intervals, sequence_index)?;
    realize_from_sequences(&sequences, config)
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
                seqwish_calls: 0,
                total_ms: start.elapsed().as_millis() as u64,
            },
        });
    }

    let num_sequences = sequences.len();
    let poa_calls = AtomicUsize::new(0);
    let sweepga_calls = AtomicUsize::new(0);
    let seqwish_calls = AtomicUsize::new(0);
    let max_depth_reached = AtomicUsize::new(0);

    let gfa = realize_recursive(sequences, config, 0, &poa_calls, &sweepga_calls, &seqwish_calls, &max_depth_reached)?;

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
            seqwish_calls: seqwish_calls.load(Ordering::Relaxed),
            total_ms: start.elapsed().as_millis() as u64,
        },
    })
}

// ---------------------------------------------------------------------------
// Recursive orchestrator
// ---------------------------------------------------------------------------

/// Base case: use seqwish for large sequence sets, POA for small ones.
fn realize_base_case(
    sequences: &[(String, SequenceMetadata)],
    config: &RealizeConfig,
    depth: usize,
    poa_calls: &AtomicUsize,
    seqwish_calls: &AtomicUsize,
) -> io::Result<String> {
    if sequences.len() > config.seqwish_threshold {
        info!(
            "realize: depth={}, {} seqs > seqwish_threshold={}, using seqwish instead of POA",
            depth, sequences.len(), config.seqwish_threshold
        );
        seqwish_calls.fetch_add(1, Ordering::Relaxed);
        let seqwish_config = crate::graph::SeqwishConfig {
            num_threads: config.num_threads,
            temp_dir: config.temp_dir.clone(),
            no_filter: false, // Must filter: raw PAF from many sequences is too large for seqwish
            // Use 1:1 filtering: with hundreds of closely-related sequences from the
            // same locus, keeping only the best alignment per query-target pair is
            // sufficient for transitive closure and drastically reduces PAF size.
            num_mappings: "1:1".to_string(),
            scaffold_filter: "1:1".to_string(),
            ..crate::graph::SeqwishConfig::default()
        };
        crate::graph::generate_gfa_seqwish_from_sequences(sequences, &seqwish_config)
    } else {
        poa_calls.fetch_add(1, Ordering::Relaxed);
        // Pass padding=0: the lacing step handles overlap trimming between adjacent
        // chunks, so we must NOT also trim padding in the POA. Doing both would shrink
        // the actual sequence while keeping the original coordinate range in path names,
        // creating gaps that the lacing can't bridge.
        let result = padded_poa_from_sequences(sequences, config.scoring_params, 0)?;
        Ok(result.gfa)
    }
}

/// Recursive realize: decide between POA base case and sweepga + partition.
fn realize_recursive(
    sequences: &[(String, SequenceMetadata)],
    config: &RealizeConfig,
    depth: usize,
    poa_calls: &AtomicUsize,
    sweepga_calls: &AtomicUsize,
    seqwish_calls: &AtomicUsize,
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

    // Too many sequences for POA-based recursion → go straight to seqwish.
    // This check MUST be before the recursive case: sweepga + par_iter on 1000+
    // sequences consumes O(N×L) memory across parallel chunks before any base
    // case is reached.  Seqwish handles large sequence counts efficiently.
    if sequences.len() > config.seqwish_threshold {
        info!(
            "realize: depth={}, {} seqs > seqwish_threshold={}, bypassing recursion for seqwish",
            depth, sequences.len(), config.seqwish_threshold
        );
        return realize_base_case(sequences, config, depth, poa_calls, seqwish_calls);
    }

    // Base case: small region or max depth → POA or seqwish.
    if total_span <= config.poa_threshold || depth >= config.max_depth {
        if depth >= config.max_depth && total_span > config.poa_threshold {
            warn!(
                "realize: hit max depth {} with region of {}bp (poa_threshold={}), forcing base case",
                depth, total_span, config.poa_threshold
            );
        }
        return realize_base_case(sequences, config, depth, poa_calls, seqwish_calls);
    }

    // Recursive case: sweepga align → partition → recurse → lace.
    //
    // Sort longest-first: partition_into_chunks assumes sequences[0] is the
    // anchor (longest).  SPOA also benefits from this order.  The sort is
    // intentionally placed *after* the seqwish-threshold check so that the
    // seqwish fallback receives sequences in their original (unsorted) order,
    // matching the output of `graph --engine seqwish`.
    let mut sequences = sequences.to_vec();
    sequences.sort_by(|a, b| b.0.len().cmp(&a.0.len()));
    let sequences = sequences; // rebind as immutable

    sweepga_calls.fetch_add(1, Ordering::Relaxed);

    // Always use all-vs-all alignment in the realize engine.
    // FastGA is efficient when run all-vs-all on a combined FASTA (indexes once).
    // Sparsification is for the top-level `impg align` CLI with hundreds of genomes,
    // not for the small sequence sets inside recursive chunks.
    let align_config = build_sweepga_config(config, sequences.len());
    let named_seqs: Vec<(String, &[u8])> = sequences
        .iter()
        .map(|(seq, meta)| (meta.alignment_name(), seq.as_bytes()))
        .collect();

    let paf_temp = crate::commands::align::sweepga_align(&named_seqs, &align_config)?;

    // Debug: save input FASTA and PAF
    if let Some(ref debug_dir) = config.debug_dir {
        let prefix = format!("{}/depth{}", debug_dir, depth);
        // Save input sequences as FASTA
        let fasta_path = format!("{}_input.fa", prefix);
        if let Ok(mut f) = std::fs::File::create(&fasta_path) {
            use std::io::Write;
            for (seq, meta) in sequences.iter() {
                let _ = writeln!(f, ">{}", meta.alignment_name());
                let _ = writeln!(f, "{}", seq);
            }
        }
        // Copy PAF
        let paf_path = format!("{}_sweepga.paf", prefix);
        let _ = std::fs::copy(paf_temp.path(), &paf_path);
        info!("realize: debug saved {}_input.fa ({} seqs) and {}_sweepga.paf", prefix, sequences.len(), prefix);
    }

    // Parse PAF records for partitioning.
    let paf_records = parse_paf_file(paf_temp.path())?;

    // If no alignments were produced, fall back to POA (or seqwish for large sets).
    if paf_records.is_empty() {
        warn!(
            "realize: sweepga produced no alignments at depth {} for {}bp region with {} sequences, falling back to base case",
            depth, total_span, sequences.len()
        );
        return realize_base_case(&sequences, config, depth, poa_calls, seqwish_calls);
    }

    // Partition into chunks along the anchor (longest) sequence.
    let chunks = partition_into_chunks(&sequences, &paf_records, config.chunk_size, config.padding);

    // If partitioning produced a single chunk that is no smaller, use base case directly.
    if chunks.len() <= 1 {
        return realize_base_case(&sequences, config, depth, poa_calls, seqwish_calls);
    }

    info!(
        "realize: depth={}, span={}bp, {} seqs → {} chunks",
        depth,
        total_span,
        sequences.len(),
        chunks.len()
    );

    // Debug: save chunk info
    if let Some(ref debug_dir) = config.debug_dir {
        let chunks_path = format!("{}/depth{}_chunks.tsv", debug_dir, depth);
        if let Ok(mut f) = std::fs::File::create(&chunks_path) {
            use std::io::Write;
            let _ = writeln!(f, "chunk\tnum_seqs\tanchor_start\tanchor_end\tmax_len\tmin_len");
            for (ci, chunk) in chunks.iter().enumerate() {
                let max_len = chunk.sequences.iter().map(|(s, _)| s.len()).max().unwrap_or(0);
                let min_len = chunk.sequences.iter().map(|(s, _)| s.len()).min().unwrap_or(0);
                let _ = writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}",
                    ci, chunk.sequences.len(), chunk.anchor_start, chunk.anchor_end, max_len, min_len);
            }
        }
    }

    // Recurse on each chunk in parallel.
    let sub_gfa_results: Vec<Option<io::Result<String>>> = chunks
        .par_iter()
        .enumerate()
        .map(|(ci, chunk)| {
            if chunk.sequences.is_empty() {
                return None;
            }
            Some(realize_recursive(
                &chunk.sequences,
                config,
                depth + 1,
                poa_calls,
                sweepga_calls,
                seqwish_calls,
                max_depth_reached,
            ).map(|gfa| {
                // Debug: save each sub-GFA
                if let Some(ref debug_dir) = config.debug_dir {
                    let gfa_path = format!("{}/depth{}_chunk{}.gfa", debug_dir, depth, ci);
                    let _ = std::fs::write(&gfa_path, &gfa);
                }
                gfa
            }))
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
    let laced = lace_subgraphs(&sub_gfas, config.temp_dir.as_deref())?;

    // Debug: save laced output
    if let Some(ref debug_dir) = config.debug_dir {
        let laced_path = format!("{}/depth{}_laced.gfa", debug_dir, depth);
        let _ = std::fs::write(&laced_path, &laced);
    }

    Ok(laced)
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
fn build_sweepga_config(config: &RealizeConfig, num_sequences: usize) -> SweepgaAlignConfig {
    // Scale k-mer frequency with sequence count.  With N nearly-identical
    // sequences from the same locus, most k-mers appear ~N times.  A fixed
    // frequency of 10 would discard almost all seeds when N >> 10, causing
    // FastGA to produce too few alignments for effective partitioning.
    // Using N * 10 mirrors what build_graph() does (num_genomes * multiplier).
    let kmer_frequency = (num_sequences * 10).max(10);
    SweepgaAlignConfig {
        num_threads: config.num_threads,
        kmer_frequency,
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
        aligner: config.aligner.clone(),
        map_pct_identity: Some("90".to_string()), // Override wfmash ANI auto-estimation
        sparsify: None,
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

        // Include the anchor's slice for this window with per-chunk metadata
        // so that path names encode per-chunk coordinates for overlap-aware lacing.
        let anchor_slice = &anchor_seq[*win_start..*win_end];
        let chunk_anchor_meta = SequenceMetadata {
            name: anchor_meta.name.clone(),
            start: anchor_meta.start + *win_start as i32,
            size: (*win_end - *win_start) as i32,
            strand: anchor_meta.strand,
            total_length: anchor_meta.total_length,
        };
        chunk_seqs.push((anchor_slice.to_string(), chunk_anchor_meta));

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
                    // Per-chunk metadata with adjusted coordinates for overlap-aware lacing.
                    let chunk_meta = SequenceMetadata {
                        name: seq_meta.name.clone(),
                        start: seq_meta.start + clamped_start as i32,
                        size: (clamped_end - clamped_start) as i32,
                        strand: seq_meta.strand,
                        total_length: seq_meta.total_length,
                    };
                    chunk_seqs.push((sub_seq.to_string(), chunk_meta));
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

/// Lace multiple GFA sub-graphs into a single GFA with proper overlap trimming.
///
/// Uses the lace pipeline from `commands::lace` to:
/// 1. Parse all sub-GFA strings into a unified CompactGraph with translated node IDs.
/// 2. Sort, deduplicate, and filter path ranges.
/// 3. Trim overlapping boundary regions between adjacent chunks.
/// 4. Link contiguous ranges.
/// 5. Remove unused nodes and compact IDs.
/// 6. Emit the final GFA as a string.
fn lace_subgraphs(sub_gfas: &[String], temp_dir: Option<&str>) -> io::Result<String> {
    if sub_gfas.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }

    if sub_gfas.len() == 1 {
        return Ok(sub_gfas[0].clone());
    }

    // 1. Build CompactGraph + path_key_ranges from in-memory GFA strings
    let mut graph = CompactGraph::new(temp_dir)?;
    let mut path_key_ranges: FxHashMap<String, Vec<RangeInfo>> = FxHashMap::default();

    for gfa_str in sub_gfas {
        let mut id_translation: FxHashMap<NodeId, NodeId> = FxHashMap::default();
        let mut temp_edges: Vec<(u64, bool, u64, bool)> = Vec::new();

        for line in gfa_str.lines() {
            if line.is_empty() || line.starts_with('#') || line.starts_with('H') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            match fields[0] {
                "S" => {
                    if fields.len() < 3 {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Invalid S line in sub-GFA: {line}"),
                        ));
                    }
                    let node_id: u64 = fields[1].parse().map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Invalid node ID '{}' in S line", fields[1]),
                        )
                    })?;
                    let sequence = fields[2].as_bytes();
                    let new_node_id = graph.add_node(sequence)?;
                    id_translation
                        .insert(NodeId::from(node_id), NodeId::from(new_node_id));
                }
                "L" => {
                    if fields.len() < 6 {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Invalid L line in sub-GFA: {line}"),
                        ));
                    }
                    let from_id: u64 = fields[1].parse().map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Invalid from node ID '{}' in L line", fields[1]),
                        )
                    })?;
                    let from_rev = fields[2] == "-";
                    let to_id: u64 = fields[3].parse().map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Invalid to node ID '{}' in L line", fields[3]),
                        )
                    })?;
                    let to_rev = fields[4] == "-";
                    temp_edges.push((from_id, from_rev, to_id, to_rev));
                }
                "P" => {
                    if fields.len() < 3 {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Invalid P line in sub-GFA: {line}"),
                        ));
                    }
                    let path_name = fields[1];
                    let nodes_str = fields[2];

                    if let Some((path_key, start, end)) = split_path_name(path_name) {
                        let mut translated_steps = Vec::new();
                        for step_str in nodes_str.split(',') {
                            if step_str.is_empty() {
                                continue;
                            }
                            let (node_str, orient) =
                                if let Some(stripped) = step_str.strip_suffix('+') {
                                    (stripped, false)
                                } else if let Some(stripped) = step_str.strip_suffix('-') {
                                    (stripped, true)
                                } else {
                                    return Err(io::Error::new(
                                        io::ErrorKind::InvalidData,
                                        format!(
                                            "Invalid step format '{step_str}' in path {path_name}"
                                        ),
                                    ));
                                };

                            let node_id: u64 = node_str.parse().map_err(|_| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!(
                                        "Invalid node ID '{node_str}' in path {path_name}"
                                    ),
                                )
                            })?;

                            if let Some(&translated_id) =
                                id_translation.get(&NodeId::from(node_id))
                            {
                                translated_steps.push(Handle::pack(translated_id, orient));
                            } else {
                                return Err(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!(
                                        "Node {node_id} in path {path_name} not found in translation map"
                                    ),
                                ));
                            }
                        }
                        if !translated_steps.is_empty() {
                            path_key_ranges
                                .entry(path_key)
                                .or_default()
                                .push(RangeInfo {
                                    start,
                                    end,
                                    steps: translated_steps,
                                });
                        }
                    }
                }
                _ => {}
            }
        }

        // Translate and add edges
        for (from_id, from_rev, to_id, to_rev) in temp_edges {
            if let (Some(&translated_from), Some(&translated_to)) = (
                id_translation.get(&NodeId::from(from_id)),
                id_translation.get(&NodeId::from(to_id)),
            ) {
                graph.add_edge(
                    translated_from.into(),
                    from_rev,
                    translated_to.into(),
                    to_rev,
                );
            }
        }
    }

    // 2. Sort and filter ranges per path key
    for ranges in path_key_ranges.values_mut() {
        sort_and_filter_ranges(ranges);
    }

    // 3. Trim overlaps and link contiguous ranges
    let graph_mutex = Arc::new(Mutex::new(graph));

    for ranges in path_key_ranges.values_mut() {
        trim_range_overlaps(ranges, &graph_mutex);
        link_contiguous_ranges(ranges, &graph_mutex);
    }

    // Unwrap graph from mutex
    let graph = Arc::try_unwrap(graph_mutex)
        .ok()
        .expect("Failed to unwrap graph mutex")
        .into_inner()
        .unwrap();

    // 4. Mark unused nodes for removal
    let nodes_to_remove = mark_nodes_for_removal(graph.node_count, &path_key_ranges);

    // 5. Build ID compaction mapping
    let max_id = graph.node_count as usize;
    let mut id_mapping: Vec<u64> = vec![0; max_id + 1];
    let mut new_id: u64 = 1;

    // 6. Emit GFA as String
    let mut output = String::new();
    output.push_str("H\tVN:Z:1.0\n");

    // Write S lines with compacted IDs (skip removed nodes)
    for node_id in 1..=graph.node_count {
        if !nodes_to_remove[node_id as usize] {
            id_mapping[node_id as usize] = new_id;
            let sequence =
                graph.get_sequence(Handle::pack(NodeId::from(node_id), false))?;
            let sequence_str =
                String::from_utf8(sequence).expect("Node sequence contains invalid UTF-8");
            output.push_str(&format!("S\t{new_id}\t{sequence_str}\n"));
            new_id += 1;
        }
    }

    // Write L lines with mapped IDs (skip edges to removed nodes)
    for edge in &graph.edges {
        let from_id = edge.source_id() as usize;
        let to_id = edge.target_id() as usize;

        if !nodes_to_remove[from_id] && !nodes_to_remove[to_id] {
            let from_mapped = id_mapping[from_id];
            let to_mapped = id_mapping[to_id];
            let from_orient = if edge.source_rev() { "-" } else { "+" };
            let to_orient = if edge.target_rev() { "-" } else { "+" };
            output.push_str(&format!(
                "L\t{from_mapped}\t{from_orient}\t{to_mapped}\t{to_orient}\t0M\n"
            ));
        }
    }

    // Write P lines: merge contiguous ranges per path key
    let mut path_keys: Vec<&String> = path_key_ranges.keys().collect();
    path_keys.sort();

    for path_key in path_keys {
        let ranges = &path_key_ranges[path_key];
        if ranges.is_empty() {
            continue;
        }

        let mut path_elements: Vec<String> = Vec::new();
        let mut path_start = ranges[0].start;
        let mut path_end = ranges[0].start;

        let mut i = 0;
        while i < ranges.len() {
            // Add steps from this range
            for handle in &ranges[i].steps {
                let node_id = id_mapping[u64::from(handle.id()) as usize];
                let orient = if handle.is_reverse() { "-" } else { "+" };
                path_elements.push(format!("{node_id}{orient}"));
            }
            let mut last_range_end = ranges[i].end;
            i += 1;

            // Merge contiguous ranges
            while i < ranges.len() && ranges[i - 1].is_contiguous_with(&ranges[i]) {
                for handle in &ranges[i].steps {
                    let node_id = id_mapping[u64::from(handle.id()) as usize];
                    let orient = if handle.is_reverse() { "-" } else { "+" };
                    path_elements.push(format!("{node_id}{orient}"));
                }
                last_range_end = ranges[i].end;
                i += 1;
            }
            path_end = last_range_end;

            // If there's a gap before the next block, emit current path and start new one
            if i < ranges.len() {
                if !path_elements.is_empty() {
                    let path_name = format!("{path_key}:{path_start}-{path_end}");
                    output.push_str(&format!("P\t{}\t{}\t*\n", path_name, path_elements.join(",")));
                    path_elements.clear();
                }
                path_start = ranges[i].start;
            }
        }

        // Write remaining path
        if !path_elements.is_empty() {
            let path_name = format!("{path_key}:{path_start}-{path_end}");
            output.push_str(&format!("P\t{}\t{}\t*\n", path_name, path_elements.join(",")));
        }
    }

    Ok(output)
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
        let result = lace_subgraphs(&[gfa.clone()], None).unwrap();
        assert_eq!(result, gfa);
    }

    #[test]
    fn test_lace_subgraphs_empty() {
        let result = lace_subgraphs(&[], None).unwrap();
        assert!(result.contains("H\tVN:Z:1.0"));
    }

    #[test]
    fn test_lace_subgraphs_contiguous_no_overlap() {
        // Two chunks with contiguous ranges: s:0-10 and s:10-20.
        // Should be linked but no trimming needed.
        let gfa1 = "H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nP\ts:0-10\t1+\t*\n".to_string();
        let gfa2 = "H\tVN:Z:1.0\nS\t1\tTGCATGCATG\nP\ts:10-20\t1+\t*\n".to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], None).unwrap();

        assert!(result.contains("H\tVN:Z:1.0"));
        assert!(result.contains("S\t1\tACGTACGTAC"));
        assert!(result.contains("S\t2\tTGCATGCATG"));
        // Contiguous ranges should be linked
        assert!(result.contains("L\t1\t+\t2\t+\t0M"));
        // Merged into a single path spanning 0-20
        assert!(result.contains("P\ts:0-20\t1+,2+\t*"));
    }

    #[test]
    fn test_lace_subgraphs_with_overlap() {
        // Two chunks with overlapping ranges.
        // Chunk 1: s:0-12 with nodes covering 12bp.
        // Chunk 2: s:8-20 with nodes covering 12bp.
        // The overlap region [8,12) (4bp) should be trimmed from chunk 2.
        //
        // Chunk 1: segment 1 (8bp "AAAAAAAA") + segment 2 (4bp "CCCC")
        //   path s:0-12 = 1+,2+
        // Chunk 2: segment 1 (4bp "CCCC") + segment 2 (8bp "GGGGGGGG")
        //   path s:8-20 = 1+,2+
        //
        // After trimming: chunk 2's first segment (4bp at position 8-12) is fully
        // in the overlap [8,12), so it's removed. Chunk 2 keeps only segment 2.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tAAAAAAAA\n\
                     S\t2\tCCCC\n\
                     L\t1\t+\t2\t+\t0M\n\
                     P\ts:0-12\t1+,2+\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tCCCC\n\
                     S\t2\tGGGGGGGG\n\
                     L\t1\t+\t2\t+\t0M\n\
                     P\ts:8-20\t1+,2+\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], None).unwrap();

        assert!(result.contains("H\tVN:Z:1.0"));
        // Chunk 1 segments preserved
        assert!(result.contains("S\t1\tAAAAAAAA"));
        assert!(result.contains("S\t2\tCCCC"));
        // Chunk 2's first segment (CCCC at 8-12) is removed by overlap trimming.
        // Chunk 2's second segment (GGGGGGGG at 12-20) is kept.
        // The trimmed node (node 3, original chunk2's node 1) is removed by ID compaction.
        // So the GGGGGGGG segment gets compacted ID 3.
        assert!(result.contains("S\t3\tGGGGGGGG"));
        // Path should span 0-20 with trimmed overlap
        assert!(
            result.contains("P\ts:0-20\t1+,2+,3+\t*"),
            "expected merged path, got:\n{result}"
        );
    }

    #[test]
    fn test_lace_subgraphs_overlap_with_node_split() {
        // Two chunks where the overlap boundary falls in the middle of a node.
        // Chunk 1: s:0-10 with single node (10bp "ACGTACGTAC")
        //   path s:0-10 = 1+
        // Chunk 2: s:8-18 with single node (10bp "ACGGGGGGGG")
        //   path s:8-18 = 1+
        //
        // Overlap is [8,10) (2bp). In chunk 2, the overlap falls within the
        // single node at positions [0,2) relative to the node. The node should
        // be split: the overlap part (first 2bp) is removed, the remaining 8bp
        // become a new node.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tACGTACGTAC\n\
                     P\ts:0-10\t1+\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tACGGGGGGGG\n\
                     P\ts:8-18\t1+\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], None).unwrap();

        assert!(result.contains("H\tVN:Z:1.0"));
        // Chunk 1's node
        assert!(result.contains("S\t1\tACGTACGTAC"));
        // Chunk 2's original node (10bp) is split: overlap [8,10) removes first 2bp.
        // Right part (8bp) becomes a new node: "GGGGGGGG"
        assert!(
            result.contains("S\t2\tGGGGGGGG"),
            "expected split node with right part, got:\n{result}"
        );
        // Merged path
        assert!(
            result.contains("P\ts:0-18\t1+,2+\t*"),
            "expected merged path, got:\n{result}"
        );
    }

    #[test]
    fn test_lace_subgraphs_multiple_paths() {
        // Two different path keys across two contiguous chunks.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tAAAA\n\
                     S\t2\tCCCC\n\
                     P\tpathA:0-4\t1+\t*\n\
                     P\tpathB:0-4\t2+\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tGGGG\n\
                     S\t2\tTTTT\n\
                     P\tpathA:4-8\t1+\t*\n\
                     P\tpathB:4-8\t2-\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], None).unwrap();

        // All segments present
        assert!(result.contains("S\t1\tAAAA"));
        assert!(result.contains("S\t2\tCCCC"));
        assert!(result.contains("S\t3\tGGGG"));
        assert!(result.contains("S\t4\tTTTT"));
        // pathA: contiguous [0-4) + [4-8) → merged to 0-8
        assert!(
            result.contains("P\tpathA:0-8\t1+,3+\t*"),
            "expected pathA merged, got:\n{result}"
        );
        assert!(result.contains("L\t1\t+\t3\t+\t0M"));
        // pathB: contiguous [0-4) + [4-8) → merged to 0-8
        assert!(
            result.contains("P\tpathB:0-8\t2+,4-\t*"),
            "expected pathB merged, got:\n{result}"
        );
        assert!(result.contains("L\t2\t+\t4\t-\t0M"));
    }

    #[test]
    fn test_lace_subgraphs_path_in_subset_of_chunks() {
        // pathA appears in both chunks, pathB only in the first chunk.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tACGT\n\
                     P\tpathA:0-4\t1+\t*\n\
                     P\tpathB:0-4\t1-\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tTGCA\n\
                     P\tpathA:4-8\t1+\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], None).unwrap();

        // pathA should be merged (contiguous 0-4 + 4-8)
        assert!(
            result.contains("P\tpathA:0-8\t1+,2+\t*"),
            "expected pathA merged, got:\n{result}"
        );
        assert!(result.contains("L\t1\t+\t2\t+\t0M"));
        // pathB should appear with just its single step
        assert!(
            result.contains("P\tpathB:0-4\t1-\t*"),
            "expected pathB unchanged, got:\n{result}"
        );
    }

    #[test]
    fn test_lace_subgraphs_link_deduplication() {
        // Duplicate links within a sub-GFA should be deduplicated via FxHashSet.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tAA\n\
                     S\t2\tCC\n\
                     L\t1\t+\t2\t+\t0M\n\
                     L\t1\t+\t2\t+\t0M\n\
                     P\ts:0-4\t1+,2+\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tGG\n\
                     P\ts:4-6\t1+\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], None).unwrap();

        // The duplicate L 1+→2+ should appear only once (FxHashSet dedup).
        let link_count = result.matches("L\t1\t+\t2\t+\t0M").count();
        assert_eq!(link_count, 1, "duplicate links should be deduplicated");
        // Contiguous linking edge 2+ → 3+
        assert!(result.contains("L\t2\t+\t3\t+\t0M"));
    }

    #[test]
    fn test_lace_subgraphs_id_compaction() {
        // Verify that unused nodes from overlap trimming are removed and IDs are compacted.
        // Chunk 1: s:0-8, two segments of 4bp each
        // Chunk 2: s:4-12, two segments of 4bp each
        // Overlap [4,8): chunk 2's first segment is fully in the overlap, removed.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tAAAA\n\
                     S\t2\tBBBB\n\
                     L\t1\t+\t2\t+\t0M\n\
                     P\ts:0-8\t1+,2+\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tBBBB\n\
                     S\t2\tCCCC\n\
                     L\t1\t+\t2\t+\t0M\n\
                     P\ts:4-12\t1+,2+\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], None).unwrap();

        // Total raw nodes: 4 (two per chunk). After trimming, chunk 2's node 1
        // (which maps to global node 3) is removed. So 3 nodes remain.
        // Compacted IDs: node 1 → 1, node 2 → 2, node 4 → 3 (node 3 removed).
        let seg_count = result.lines().filter(|l| l.starts_with("S\t")).count();
        assert_eq!(
            seg_count, 3,
            "expected 3 segments after compaction, got {seg_count}:\n{result}"
        );

        // IDs should be contiguous 1,2,3
        assert!(result.contains("S\t1\t"));
        assert!(result.contains("S\t2\t"));
        assert!(result.contains("S\t3\t"));
        // No node 4 should exist
        assert!(
            !result.contains("S\t4\t"),
            "node 4 should be compacted away"
        );
    }

    #[test]
    fn test_lace_subgraphs_three_contiguous_chunks() {
        // Lace three contiguous sub-GFAs to verify chaining works across >2 chunks.
        let gfa1 = "H\tVN:Z:1.0\nS\t1\tAA\nP\ts:0-2\t1+\t*\n".to_string();
        let gfa2 = "H\tVN:Z:1.0\nS\t1\tCC\nP\ts:2-4\t1-\t*\n".to_string();
        let gfa3 = "H\tVN:Z:1.0\nS\t1\tGG\nP\ts:4-6\t1+\t*\n".to_string();
        let result = lace_subgraphs(&[gfa1, gfa2, gfa3], None).unwrap();

        assert!(result.contains("S\t1\tAA"));
        assert!(result.contains("S\t2\tCC"));
        assert!(result.contains("S\t3\tGG"));
        // Linking edges between contiguous ranges
        assert!(result.contains("L\t1\t+\t2\t-\t0M"));
        assert!(result.contains("L\t2\t-\t3\t+\t0M"));
        assert!(
            result.contains("P\ts:0-6\t1+,2-,3+\t*"),
            "expected merged path, got:\n{result}"
        );
    }

    #[test]
    fn test_lace_subgraphs_multi_segment_contiguous() {
        // Each sub-GFA has a path with multiple segments, contiguous ranges.
        let gfa1 = "H\tVN:Z:1.0\n\
                     S\t1\tACGT\n\
                     S\t2\tTTTT\n\
                     L\t1\t+\t2\t+\t0M\n\
                     P\ts:0-8\t1+,2+\t*\n"
            .to_string();
        let gfa2 = "H\tVN:Z:1.0\n\
                     S\t1\tGGGG\n\
                     S\t2\tCCCC\n\
                     L\t1\t+\t2\t-\t0M\n\
                     P\ts:8-16\t1+,2-\t*\n"
            .to_string();
        let result = lace_subgraphs(&[gfa1, gfa2], None).unwrap();

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
        assert!(
            result.contains("P\ts:0-16\t1+,2+,3+,4-\t*"),
            "expected merged path, got:\n{result}"
        );
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
    #[cfg_attr(target_os = "macos", ignore)] // FastGA crashes on macOS
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
        config.max_depth = 1;


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
    #[cfg_attr(target_os = "macos", ignore)] // FastGA crashes on macOS
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
        config.max_depth = 1;


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
    #[cfg_attr(target_os = "macos", ignore)] // FastGA crashes on macOS
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
    #[cfg_attr(target_os = "macos", ignore)] // FastGA crashes on macOS
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
        config.max_depth = 1;


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
    #[cfg_attr(target_os = "macos", ignore)] // FastGA crashes on macOS
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
        config.max_depth = 1;


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
    #[cfg_attr(target_os = "macos", ignore)] // FastGA crashes on macOS
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
        config.max_depth = 1;


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
    #[cfg_attr(target_os = "macos", ignore)] // FastGA crashes on macOS
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
    fn test_partition_metadata_per_chunk_coordinates() {
        // Verify that chunk metadata has per-chunk coordinates reflecting
        // actual slice positions (needed for overlap-aware lacing).
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

        // Each chunk's anchor metadata should have adjusted start/size
        // reflecting the actual slice coordinates within that chunk.
        for (i, chunk) in chunks.iter().enumerate() {
            let (seq, anchor_meta) = chunk.sequences.iter().find(|(_, m)| m.name == "anchor").unwrap();
            // start should be original_start + win_start
            assert!(anchor_meta.start >= 100,
                "chunk {} anchor start should be >= original start 100, got {}", i, anchor_meta.start);
            // size should match the actual sequence length
            assert_eq!(anchor_meta.size as usize, seq.len(),
                "chunk {} anchor size should match sequence length", i);
            // name and total_length should be preserved
            assert_eq!(anchor_meta.name, "anchor",
                "chunk {} anchor name should be preserved", i);
            assert_eq!(anchor_meta.total_length, 1000,
                "chunk {} anchor total_length should be preserved", i);
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
