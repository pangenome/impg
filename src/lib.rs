// lib.rs
#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]
pub mod agc_index;
pub mod syng_ffi;
pub mod syng;
pub mod alignment_record;
pub mod commands;
pub mod faidx;
pub mod forest_map;
pub mod graph;
pub mod impg;
pub mod impg_index;
pub mod multi_impg;
pub mod onealn;
pub mod paf;
pub mod seqidx;
pub mod sequence_index;
pub mod smooth;
pub mod subset_filter;
pub mod tpa_parser;

/// GFA engine selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum GfaEngine {
    /// Alignment + seqwish graph induction + smoothing + gfaffix normalization
    Pggb,
    /// Alignment + seqwish graph induction + gfaffix normalization (no smoothing)
    Seqwish,
    /// Single-pass partial order alignment (POA)
    Poa,
}

/// Resolved engine configuration passed to subcommand functions.
pub struct EngineOpts {
    pub engine: GfaEngine,
    pub num_threads: usize,
    pub no_filter: bool,
    /// Optional directory to save intermediate debug files (PAFs, FASTAs, etc.)
    pub debug_dir: Option<String>,
    /// Unified sparsification strategy.
    pub sparsify: sweepga::knn_graph::SparsificationStrategy,
    /// Mash distance parameters for sparsification sketching.
    pub mash_params: sweepga::knn_graph::MashParams,
    // Alignment filtering parameters (shared with graph/align commands)
    pub aligner: String,
    pub num_mappings: String,
    pub scaffold_jump: u64,
    pub scaffold_mass: u64,
    pub scaffold_filter: String,
    pub overlap: f64,
    pub min_identity: f64,
    pub scaffold_dist: u64,
    pub min_map_length: u64,
    pub min_aln_length: u64,
    pub frequency_multiplier: usize,
    // Seqwish graph induction parameters
    pub repeat_max: u64,
    pub min_repeat_dist: u64,
    pub min_match_len: u64,
    pub sparse_factor: f32,
    pub transclose_batch: u64,
    pub disk_backed: bool,
    /// Optional temp directory for intermediate files
    pub temp_dir: Option<String>,
    /// Batch genome alignment to limit resource usage per batch (e.g. "2G", "500M").
    pub batch_bytes: Option<String>,
    // Smoothxg-style smoothing parameters (pggb engine)
    /// Target POA length(s) per pass — one value per smoothing pass (default: [700, 1100]).
    pub target_poa_lengths: Vec<usize>,
    pub max_node_length: usize,
    pub poa_padding_fraction: f64,
    /// When set, activates the partitioned GFA pipeline: build per-partition GFA,
    /// lace together, then run gfaffix once at the end.
    pub partition_size: Option<usize>,
}

/// Minimal ImpgIndex wrapper around a SequenceIndex for syng query path.
/// Only `seq_index()` is functional; query methods are not called.
struct SeqIndexWrapper {
    seq_index: seqidx::SequenceIndex,
}

impl impg_index::ImpgIndex for SeqIndexWrapper {
    fn seq_index(&self) -> &seqidx::SequenceIndex {
        &self.seq_index
    }

    fn query(
        &self, _: u32, _: i32, _: i32, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>, _: bool,
    ) -> Vec<impg::AdjustedInterval> {
        unimplemented!("not used in syng path")
    }

    fn query_with_cache(
        &self, _: u32, _: i32, _: i32, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>,
        _: &rustc_hash::FxHashMap<(u32, u64), Vec<impg::CigarOp>>,
    ) -> Vec<impg::AdjustedInterval> {
        unimplemented!("not used in syng path")
    }

    fn populate_cigar_cache(
        &self, _: u32, _: i32, _: i32, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>,
        _: &mut rustc_hash::FxHashMap<(u32, u64), Vec<impg::CigarOp>>,
    ) {
        unimplemented!("not used in syng path")
    }

    fn query_transitive_dfs(
        &self, _: u32, _: i32, _: i32,
        _: Option<&rustc_hash::FxHashMap<u32, impg::SortedRanges>>,
        _: u16, _: i32, _: i32, _: Option<i32>, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>, _: bool,
        _: Option<&subset_filter::SubsetFilter>,
    ) -> Vec<impg::AdjustedInterval> {
        unimplemented!("not used in syng path")
    }

    fn query_transitive_bfs(
        &self, _: u32, _: i32, _: i32,
        _: Option<&rustc_hash::FxHashMap<u32, impg::SortedRanges>>,
        _: u16, _: i32, _: i32, _: Option<i32>, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>, _: bool,
        _: Option<&subset_filter::SubsetFilter>,
    ) -> Vec<impg::AdjustedInterval> {
        unimplemented!("not used in syng path")
    }

    fn get_or_load_tree(&self, _: u32) -> Option<std::sync::Arc<coitrees::BasicCOITree<impg::QueryMetadata, u32>>> {
        unimplemented!("not used in syng path")
    }

    fn target_ids(&self) -> Vec<u32> {
        unimplemented!("not used in syng path")
    }

    fn remove_cached_tree(&self, _: u32) {
        unimplemented!("not used in syng path")
    }

    fn sequence_files(&self) -> &[String] {
        &[]
    }
}

/// Dispatch GFA generation using a SequenceIndex directly (for syng queries).
///
/// This avoids needing a full ImpgIndex when the intervals come from syng
/// rather than from alignment-based queries.
pub fn dispatch_gfa_engine_with_seq_index(
    seq_index: &seqidx::SequenceIndex,
    query_intervals: &[coitrees::Interval<u32>],
    sequence_index: &sequence_index::UnifiedSequenceIndex,
    scoring_params: Option<(u8, u8, u8, u8, u8, u8)>,
    engine_opts: &EngineOpts,
) -> std::io::Result<String> {
    let wrapper = SeqIndexWrapper {
        seq_index: seq_index.clone(),
    };
    dispatch_gfa_engine(&wrapper, query_intervals, sequence_index, scoring_params, engine_opts)
}

/// Build a `SeqwishConfig` from `EngineOpts`. Extracted so
/// `partitioned_gfa_pipeline` can build the config ONCE and reuse it for every
/// partition, instead of reconstructing (and re-cloning ~14 String fields) per
/// partition inside `dispatch_gfa_engine`.
fn build_seqwish_config(engine_opts: &EngineOpts) -> graph::SeqwishConfig {
    graph::SeqwishConfig {
        num_threads: engine_opts.num_threads,
        no_filter: engine_opts.no_filter,
        aligner: engine_opts.aligner.clone(),
        num_mappings: engine_opts.num_mappings.clone(),
        scaffold_jump: engine_opts.scaffold_jump,
        scaffold_mass: engine_opts.scaffold_mass,
        scaffold_filter: engine_opts.scaffold_filter.clone(),
        overlap: engine_opts.overlap,
        min_identity: engine_opts.min_identity,
        scaffold_dist: engine_opts.scaffold_dist,
        min_map_length: engine_opts.min_map_length,
        min_aln_length: engine_opts.min_aln_length,
        frequency_multiplier: engine_opts.frequency_multiplier,
        debug_dir: engine_opts.debug_dir.clone(),
        sparsify: engine_opts.sparsify.clone(),
        mash_params: engine_opts.mash_params.clone(),
        repeat_max: engine_opts.repeat_max,
        min_repeat_dist: engine_opts.min_repeat_dist,
        min_match_len: engine_opts.min_match_len,
        sparse_factor: engine_opts.sparse_factor,
        transclose_batch: engine_opts.transclose_batch,
        use_in_memory: !engine_opts.disk_backed,
        temp_dir: engine_opts.temp_dir.clone(),
        batch_bytes: engine_opts.batch_bytes.clone(),
        ..graph::SeqwishConfig::default()
    }
}

/// Dispatch GFA generation to the selected engine.
///
/// Shared by `query -o gfa` and `partition -o gfa` so the engine match
/// logic lives in one place.
pub fn dispatch_gfa_engine(
    impg: &impl impg_index::ImpgIndex,
    query_intervals: &[coitrees::Interval<u32>],
    sequence_index: &sequence_index::UnifiedSequenceIndex,
    scoring_params: Option<(u8, u8, u8, u8, u8, u8)>,
    engine_opts: &EngineOpts,
) -> std::io::Result<String> {
    // Create debug dir once — needed by both seqwish and pggb pipelines.
    if let Some(ref dir) = engine_opts.debug_dir {
        std::fs::create_dir_all(dir).map_err(|e| {
            std::io::Error::other(format!("Failed to create debug dir '{}': {}", dir, e))
        })?;
    }
    let seqwish_config = build_seqwish_config(engine_opts);
    dispatch_gfa_engine_with_config(
        impg,
        query_intervals,
        sequence_index,
        scoring_params,
        engine_opts,
        &seqwish_config,
    )
}

/// Internal variant of `dispatch_gfa_engine` that takes a pre-built
/// `SeqwishConfig`. Used by `partitioned_gfa_pipeline` to avoid rebuilding
/// the config on every partition.
fn dispatch_gfa_engine_with_config(
    impg: &impl impg_index::ImpgIndex,
    query_intervals: &[coitrees::Interval<u32>],
    sequence_index: &sequence_index::UnifiedSequenceIndex,
    scoring_params: Option<(u8, u8, u8, u8, u8, u8)>,
    engine_opts: &EngineOpts,
    seqwish_config: &graph::SeqwishConfig,
) -> std::io::Result<String> {
    let skip_normalize = engine_opts.partition_size.is_some();

    match engine_opts.engine {
        GfaEngine::Poa => {
            let params = scoring_params.ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "POA scoring parameters required for poa engine",
                )
            })?;
            let gfa =
                graph::generate_gfa_from_intervals(impg, query_intervals, sequence_index, params)?;
            if skip_normalize {
                Ok(gfa)
            } else {
                graph::normalize_and_sort(gfa, engine_opts.num_threads)
            }
        }
        GfaEngine::Seqwish => {
            let gfa = graph::generate_gfa_seqwish_from_intervals(
                impg,
                query_intervals,
                sequence_index,
                seqwish_config,
            )?;
            if skip_normalize {
                Ok(gfa)
            } else {
                graph::normalize_and_sort(gfa, engine_opts.num_threads)
            }
        }
        GfaEngine::Pggb => {
            // Step 1: seqwish graph induction (shared config, same as seqwish engine)
            let raw_gfa = graph::generate_gfa_seqwish_from_intervals(
                impg,
                query_intervals,
                sequence_index,
                seqwish_config,
            )?;

            // Sort for 1D layout (required by smoothxg block decomposition)
            let raw_gfa = graph::sort_gfa(&raw_gfa, engine_opts.num_threads)?;

            // Step 2: smooth
            let n_haps = query_intervals.len().max(1);
            let smooth_config = smooth::SmoothConfig {
                num_threads: engine_opts.num_threads,
                target_poa_lengths: engine_opts.target_poa_lengths.clone(),
                max_node_length: engine_opts.max_node_length,
                poa_padding_fraction: engine_opts.poa_padding_fraction,
                pre_sorted: true,
                ..smooth::SmoothConfig::new(n_haps)
            };
            let smoothed = smooth::smooth_gfa(&raw_gfa, &smooth_config)?;

            // Step 3: gfaffix + sort (skipped when partitioned — done once after lacing)
            if skip_normalize {
                Ok(smoothed)
            } else {
                graph::normalize_and_sort(smoothed, engine_opts.num_threads)
            }
        }
    }
}

/// Run the partitioned GFA pipeline: build per-partition GFA (sequentially,
/// each using all threads), lace them together, then run a single final
/// gfaffix normalization + sort.
pub fn partitioned_gfa_pipeline(
    partitions: &[Vec<coitrees::Interval<u32>>],
    impg: &impl impg_index::ImpgIndex,
    sequence_index: &sequence_index::UnifiedSequenceIndex,
    scoring_params: Option<(u8, u8, u8, u8, u8, u8)>,
    engine_opts: &EngineOpts,
) -> std::io::Result<String> {
    use log::info;
    use std::time::Instant;

    if partitions.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }

    // Create debug dir once here (not inside each dispatch_gfa_engine call).
    if let Some(ref dir) = engine_opts.debug_dir {
        std::fs::create_dir_all(dir).map_err(|e| {
            std::io::Error::other(format!("Failed to create debug dir '{}': {}", dir, e))
        })?;
    }

    // Perf 3: build SeqwishConfig ONCE, reuse across all partitions.
    let seqwish_config = build_seqwish_config(engine_opts);

    // Perf 5: compute per-partition bp in a single pass, derive total from the sum.
    let per_partition_bp: Vec<u64> = partitions
        .iter()
        .map(|ivs| {
            ivs.iter()
                .map(|iv| (iv.last - iv.first).unsigned_abs() as u64)
                .sum()
        })
        .collect();
    let total_bp: u64 = per_partition_bp.iter().sum();

    // 1. Generate per-partition GFAs sequentially (each partition uses all threads)
    let total_partitions = partitions.len();
    let mut sub_gfas: Vec<String> = Vec::with_capacity(total_partitions);
    let mut total_partitioned_bp: u64 = 0;
    let pipeline_start = Instant::now();
    for (idx, intervals) in partitions.iter().enumerate() {
        let num_regions = intervals.len();
        let current_partition_length = per_partition_bp[idx];
        let current_percentage = if total_bp > 0 {
            (current_partition_length as f64 / total_bp as f64) * 100.0
        } else {
            0.0
        };
        let total_percentage = if total_bp > 0 {
            (total_partitioned_bp as f64 / total_bp as f64) * 100.0
        } else {
            0.0
        };
        let current_percentage_str = if current_percentage < 0.0001 {
            format!("{current_percentage:.4e}%")
        } else {
            format!("{current_percentage:.4}%")
        };
        let total_percentage_str = if total_percentage < 0.0001 {
            format!("{total_percentage:.4e}%")
        } else {
            format!("{total_percentage:.4}%")
        };
        info!(
            "[partitioned] Building GFA for partition {}/{} with {} intervals: {} bp this partition ({}), {} bp total ({})",
            idx + 1,
            total_partitions,
            num_regions,
            current_partition_length,
            current_percentage_str,
            total_partitioned_bp,
            total_percentage_str,
        );
        // Perf 7: per-partition elapsed time.
        let part_start = Instant::now();
        let gfa = dispatch_gfa_engine_with_config(
            impg,
            intervals,
            sequence_index,
            scoring_params,
            engine_opts,
            &seqwish_config,
        )?;
        let part_elapsed = part_start.elapsed();
        let pipeline_elapsed = pipeline_start.elapsed();
        info!(
            "[partitioned] Completed partition {}/{} in {:.1}s (pipeline elapsed: {:.1}s)",
            idx + 1,
            total_partitions,
            part_elapsed.as_secs_f64(),
            pipeline_elapsed.as_secs_f64(),
        );
        sub_gfas.push(gfa);
        total_partitioned_bp += current_partition_length;
    }

    // 2. Lace all partition GFAs together
    info!(
        "[partitioned] Lacing {} partition GFAs",
        sub_gfas.len()
    );
    let laced = commands::lace::lace_subgraphs(&sub_gfas, None)?;
    // Free all sub-partition GFAs before gfaffix/sort
    drop(sub_gfas);

    // 3. Single final gfaffix normalization + sort
    info!("[partitioned] Running final gfaffix normalization");
    graph::normalize_and_sort(laced, engine_opts.num_threads)
}
