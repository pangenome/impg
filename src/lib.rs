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
    pub pipeline: commands::graph::GraphBuildConfig,
    // Smoothxg-style smoothing parameters (pggb engine)
    /// Target POA length(s) per pass — one value per smoothing pass (default: [700, 1100]).
    pub target_poa_lengths: Vec<usize>,
    pub max_node_length: usize,
    pub poa_padding_fraction: f64,
    /// When set, activates the partitioned GFA pipeline: build per-partition GFA,
    /// lace together, then run gfaffix once at the end.
    pub partition_size: Option<usize>,
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
    if let Some(ref dir) = engine_opts.pipeline.debug_dir {
        std::fs::create_dir_all(dir).map_err(|e| {
            std::io::Error::other(format!("Failed to create debug dir '{}': {}", dir, e))
        })?;
    }
    dispatch_gfa_engine_inner(
        impg,
        query_intervals,
        sequence_index,
        scoring_params,
        engine_opts,
    )
}

/// Dispatch that assumes `debug_dir` (if any) has already been created.
fn dispatch_gfa_engine_inner(
    impg: &impl impg_index::ImpgIndex,
    query_intervals: &[coitrees::Interval<u32>],
    sequence_index: &sequence_index::UnifiedSequenceIndex,
    scoring_params: Option<(u8, u8, u8, u8, u8, u8)>,
    engine_opts: &EngineOpts,
) -> std::io::Result<String> {
    let skip_normalize = engine_opts.partition_size.is_some();
    let num_threads = engine_opts.pipeline.num_threads;

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
                graph::normalize_and_sort(gfa, num_threads)
            }
        }
        GfaEngine::Seqwish => {
            let gfa = graph::generate_gfa_seqwish_from_intervals(
                impg,
                query_intervals,
                sequence_index,
                &engine_opts.pipeline,
            )?;
            if skip_normalize {
                Ok(gfa)
            } else {
                graph::normalize_and_sort(gfa, num_threads)
            }
        }
        GfaEngine::Pggb => {
            // Step 1: seqwish graph induction (shared pipeline config)
            let raw_gfa = graph::generate_gfa_seqwish_from_intervals(
                impg,
                query_intervals,
                sequence_index,
                &engine_opts.pipeline,
            )?;

            // Normalize seqwish line order before sort_gfa (seqwish emits
            // L-edges in a thread-dependent order — same content, unstable
            // ordering — and ygs_sort is parse-order-sensitive). Must
            // match `run_graph_build_pggb` for the `query -o gfa` and
            // `graph --sequence-files` paths to produce the same layout.
            let raw_gfa = {
                let mut lines: Vec<&str> = raw_gfa.lines().collect();
                lines.sort_unstable();
                lines.join("\n") + "\n"
            };

            // Sort for 1D layout (required by smoothxg block decomposition)
            let raw_gfa = graph::sort_gfa(&raw_gfa, num_threads)?;

            // Step 2: smooth. `n_haps` is the number of unique haplotypes
            // in the input (distinct `SAMPLE#HAPLOTYPE` prefixes), NOT
            // the number of query intervals. `max_block_weight` scales
            // with this, and using the interval count (typically many
            // multiples larger) inflates block size → different graph
            // structure. Must match the `graph` subcommand, which counts
            // haplotypes from the input FASTAs via the same helper.
            let n_haps = sweepga::pansn::count_pansn_keys(
                query_intervals
                    .iter()
                    .filter_map(|iv| impg.seq_index().get_name(iv.metadata)),
                sweepga::pansn::PanSnLevel::Haplotype,
            );
            let smooth_config = smooth::SmoothConfig {
                num_threads,
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
                graph::normalize_and_sort(smoothed, num_threads)
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
    if let Some(ref dir) = engine_opts.pipeline.debug_dir {
        std::fs::create_dir_all(dir).map_err(|e| {
            std::io::Error::other(format!("Failed to create debug dir '{}': {}", dir, e))
        })?;
    }

    // Compute per-partition bp in a single pass, derive total from the sum.
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
        let gfa = dispatch_gfa_engine_inner(
            impg,
            intervals,
            sequence_index,
            scoring_params,
            engine_opts,
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
    graph::normalize_and_sort(laced, engine_opts.pipeline.num_threads)
}
