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
pub mod smooth;
pub mod tpa_parser;
pub mod seqidx;
pub mod sequence_index;
pub mod subset_filter;

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
    match engine_opts.engine {
        GfaEngine::Poa => {
            let params = scoring_params.ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "POA scoring parameters required for poa engine",
                )
            })?;
            graph::generate_gfa_from_intervals(
                impg,
                query_intervals,
                sequence_index,
                params,
                engine_opts.num_threads,
            )
        }
        GfaEngine::Seqwish => {
            if let Some(ref dir) = engine_opts.debug_dir {
                std::fs::create_dir_all(dir).map_err(|e| {
                    std::io::Error::other(format!("Failed to create debug dir '{}': {}", dir, e))
                })?;
            }
            let config = graph::SeqwishConfig {
                num_threads: engine_opts.num_threads,
                no_filter: engine_opts.no_filter,
                debug_dir: engine_opts.debug_dir.clone(),
                sparsify: engine_opts.sparsify.clone(),
                ..graph::SeqwishConfig::default()
            };
            graph::generate_gfa_seqwish_from_intervals(
                impg,
                query_intervals,
                sequence_index,
                &config,
            )
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
            Ok(result.gfa)
        }
        GfaEngine::Pggb => {
            // Step 1: run seqwish pipeline
            let seqwish_config = graph::SeqwishConfig {
                num_threads: engine_opts.num_threads,
                no_filter: engine_opts.no_filter,
                debug_dir: engine_opts.debug_dir.clone(),
                sparsify: engine_opts.sparsify.clone(),
                ..graph::SeqwishConfig::default()
            };
            let raw_gfa = graph::generate_gfa_seqwish_from_intervals(
                impg,
                query_intervals,
                sequence_index,
                &seqwish_config,
            )?;

            // Step 2: smooth
            let n_haps = query_intervals.len().max(1);
            let smooth_config = smooth::SmoothConfig {
                num_threads: engine_opts.num_threads,
                pre_sorted: true,
                ..smooth::SmoothConfig::new(n_haps)
            };
            let smoothed = smooth::smooth_gfa(&raw_gfa, &smooth_config)?;

            // Step 3: gfaffix normalization (optional, graceful fallback)
            let normalized = graph::run_gfaffix(&smoothed, engine_opts.num_threads)
                .unwrap_or_else(|e| {
                    log::debug!("[pggb] gfaffix skipped: {}", e);
                    smoothed
                });

            // Step 4: final sort
            graph::sort_gfa(&normalized, engine_opts.num_threads)
        }
    }
}
