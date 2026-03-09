// lib.rs
pub mod agc_index;
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
pub mod realize;
pub mod seqidx;
pub mod sequence_index;
pub mod smooth;
pub mod subset_filter;
pub mod tpa_parser;

/// GFA engine selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum GfaEngine {
    /// Recursive partitioning + POA + lacing
    Recursive,
    /// Seqwish graph induction via transitive closure
    Seqwish,
    /// Flat single-pass partial order alignment
    Poa,
    /// Seqwish + smoothxg-style smoothing + gfaffix normalization
    Pggb,
}

/// Resolved engine configuration passed to subcommand functions.
pub struct EngineOpts {
    pub engine: GfaEngine,
    pub recursive_config: Option<realize::RealizeConfig>,
    pub num_threads: usize,
    pub no_filter: bool,
    /// Optional directory to save intermediate debug files (PAFs, FASTAs, etc.)
    pub debug_dir: Option<String>,
    /// Wfmash mapping sparsification: "auto" or a float string like "0.1".
    pub sparsify: Option<String>,
    // Seqwish graph induction parameters
    pub repeat_max: u64,
    pub min_repeat_dist: u64,
    pub min_match_len: u64,
    pub sparse_factor: f32,
    pub transclose_batch: u64,
    pub disk_backed: bool,
    // Smoothxg-style smoothing parameters (pggb engine)
    /// Target POA length(s) per pass — one value per smoothing pass (default: [700, 1100]).
    pub target_poa_lengths: Vec<usize>,
    pub max_node_length: usize,
    pub poa_padding_fraction: f64,
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

    // Shared seqwish config — identical for the seqwish and pggb engines; built
    // once here so neither arm repeats the field list.
    let seqwish_config = graph::SeqwishConfig {
        num_threads: engine_opts.num_threads,
        no_filter: engine_opts.no_filter,
        debug_dir: engine_opts.debug_dir.clone(),
        sparsify: engine_opts.sparsify.clone(),
        repeat_max: engine_opts.repeat_max,
        min_repeat_dist: engine_opts.min_repeat_dist,
        min_match_len: engine_opts.min_match_len,
        sparse_factor: engine_opts.sparse_factor,
        transclose_batch: engine_opts.transclose_batch,
        use_in_memory: !engine_opts.disk_backed,
        ..graph::SeqwishConfig::default()
    };

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
            graph::normalize_and_sort(gfa, engine_opts.num_threads)
        }
        GfaEngine::Seqwish => {
            let gfa = graph::generate_gfa_seqwish_from_intervals(
                impg,
                query_intervals,
                sequence_index,
                &seqwish_config,
            )?;
            graph::normalize_and_sort(gfa, engine_opts.num_threads)
        }
        GfaEngine::Recursive => {
            let config = engine_opts.recursive_config.as_ref().ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "RecursiveOpts config required for recursive engine",
                )
            })?;
            let result = realize::realize(impg, query_intervals, sequence_index, config)?;
            log::info!(
                "Recursive engine stats: {} sequences, max_depth={}, poa_calls={}, sweepga_calls={}, seqwish_calls={}, {}ms",
                result.stats.num_sequences,
                result.stats.max_depth_reached,
                result.stats.poa_calls,
                result.stats.sweepga_calls,
                result.stats.seqwish_calls,
                result.stats.total_ms,
            );
            graph::normalize_and_sort(result.gfa, engine_opts.num_threads)
        }
        GfaEngine::Pggb => {
            // Step 1: seqwish graph induction (shared config, same as seqwish engine)
            let raw_gfa = graph::generate_gfa_seqwish_from_intervals(
                impg,
                query_intervals,
                sequence_index,
                &seqwish_config,
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

            // Step 3: gfaffix + sort
            graph::normalize_and_sort(smoothed, engine_opts.num_threads)
        }
    }
}
