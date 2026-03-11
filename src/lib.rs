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
        ..graph::SeqwishConfig::default()
    };

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
                &seqwish_config,
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
    partitions: &[(usize, Vec<coitrees::Interval<u32>>)],
    impg: &impl impg_index::ImpgIndex,
    sequence_index: &sequence_index::UnifiedSequenceIndex,
    scoring_params: Option<(u8, u8, u8, u8, u8, u8)>,
    engine_opts: &EngineOpts,
) -> std::io::Result<String> {
    use log::info;

    if partitions.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }

    // 1. Generate per-partition GFAs sequentially (each partition uses all threads)
    let mut sub_gfas: Vec<String> = Vec::with_capacity(partitions.len());
    for (partition_num, intervals) in partitions {
        info!(
            "[partitioned] Building GFA for partition {} ({} intervals)",
            partition_num,
            intervals.len()
        );
        let gfa = dispatch_gfa_engine(impg, intervals, sequence_index, scoring_params, engine_opts)?;
        sub_gfas.push(gfa);
    }

    // 2. Lace all partition GFAs together
    info!(
        "[partitioned] Lacing {} partition GFAs",
        sub_gfas.len()
    );
    let laced = commands::lace::lace_subgraphs(&sub_gfas, None)?;

    // 3. Single final gfaffix normalization + sort
    info!("[partitioned] Running final gfaffix normalization");
    graph::normalize_and_sort(laced, engine_opts.num_threads)
}
