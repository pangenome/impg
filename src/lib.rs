// lib.rs
#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]
pub mod agc_index;
pub mod syng_ffi;
pub mod syng;
pub mod syng_parallel;
pub mod alignment_record;
pub mod commands;
pub mod faidx;
pub mod forest_map;
pub mod graph;
pub mod impg;
pub mod impg_index;
pub mod multi_impg;
pub mod onealn;
pub mod pack;
pub mod paf;
pub mod projection;
pub mod seqidx;
pub mod sequence_index;
pub mod smooth;
pub mod subset_filter;
pub mod syng_graph;
pub mod syng_graph_norm;
pub mod syng_transitive;
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
    /// Syng-native: anchor-seeded BiWFA between allwave-sparsified pairs,
    /// strand-groomed to partition reference frame, indel-left-aligned,
    /// then fed to seqwish. Fast path for many-haplotype low-divergence
    /// partitions — skips running an external aligner (wfmash/fastGA) by
    /// reusing the syng anchor chains the partition step already computed.
    ///
    /// Pipeline stages (to be filled in):
    ///   1. per-partition distance matrix from syng shared-anchor counts
    ///   2. pair selection via `sweepga::knn_graph::extract_tree_pairs_from_matrix`
    ///   3. strand groom: reverse-complement members to match partition-ref frame
    ///   4. anchor-seeded BiWFA: match blocks from syncmers, fill inter-anchor gaps
    ///   5. left-align indels in partition-ref frame
    ///   6. emit PAF → seqwish → GFA
    SyngNative,
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

/// Minimal ImpgIndex wrapper around a SequenceIndex for syng query path.
/// Only `seq_index()` is functional; query methods are not called.
struct SeqIndexWrapper {
    seq_index: seqidx::SequenceIndex,
}

impl impg_index::ImpgIndex for SeqIndexWrapper {
    // Bare-minimum wrapper for `dispatch_gfa_engine_with_seq_index`:
    // every method returns a safe empty value so this can sit inside
    // a heterogeneous `Vec<Arc<dyn ImpgIndex + Send + Sync>>` without
    // risking a panic. Only `seq_index()` is called in practice.
    fn seq_index(&self) -> &seqidx::SequenceIndex {
        &self.seq_index
    }

    fn query(
        &self, _: u32, _: i32, _: i32, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>, _: bool,
    ) -> Vec<impg::AdjustedInterval> {
        Vec::new()
    }

    fn query_with_cache(
        &self, _: u32, _: i32, _: i32, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>,
        _: &rustc_hash::FxHashMap<(u32, u64), Vec<impg::CigarOp>>,
    ) -> Vec<impg::AdjustedInterval> {
        Vec::new()
    }

    fn populate_cigar_cache(
        &self, _: u32, _: i32, _: i32, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>,
        _: &mut rustc_hash::FxHashMap<(u32, u64), Vec<impg::CigarOp>>,
    ) {
    }

    fn query_transitive_dfs(
        &self, _: u32, _: i32, _: i32,
        _: Option<&rustc_hash::FxHashMap<u32, impg::SortedRanges>>,
        _: u16, _: i32, _: i32, _: Option<i32>, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>, _: bool,
        _: Option<&subset_filter::SubsetFilter>,
    ) -> Vec<impg::AdjustedInterval> {
        Vec::new()
    }

    fn query_transitive_bfs(
        &self, _: u32, _: i32, _: i32,
        _: Option<&rustc_hash::FxHashMap<u32, impg::SortedRanges>>,
        _: u16, _: i32, _: i32, _: Option<i32>, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>, _: bool,
        _: Option<&subset_filter::SubsetFilter>,
    ) -> Vec<impg::AdjustedInterval> {
        Vec::new()
    }

    fn get_or_load_tree(&self, _: u32) -> Option<std::sync::Arc<coitrees::BasicCOITree<impg::QueryMetadata, u32>>> {
        None
    }

    fn target_ids(&self) -> Vec<u32> {
        (0..self.seq_index.len() as u32).collect()
    }

    fn remove_cached_tree(&self, _: u32) {
    }

    fn sequence_files(&self) -> &[String] {
        &[]
    }
}

/// Full ImpgIndex adapter around a SyngIndex.
///
/// Translates `query_transitive_dfs`/`bfs` calls into `SyngIndex::query_region`
/// so that `partition_alignments` (and any other code taking `&impl ImpgIndex`)
/// can work transparently with syng-based queries.
pub struct SyngImpgWrapper {
    syng_index: syng::SyngIndex,
    seq_index: seqidx::SequenceIndex,
    syng_padding: u64,
    /// Optional chain filter: only chains with at least this many anchors
    /// survive. `None` means no chain-based filtering (raw query_region).
    /// Populated via [`Self::with_chain_filter`].
    min_chain_anchors: Option<usize>,
    min_chain_fraction: f64,
}

impl SyngImpgWrapper {
    pub fn new(syng_index: syng::SyngIndex, seq_index: seqidx::SequenceIndex, syng_padding: u64) -> Self {
        Self {
            syng_index,
            seq_index,
            syng_padding,
            min_chain_anchors: None,
            min_chain_fraction: 0.0,
        }
    }

    /// Enable chain-based filtering in `query_via_syng`. When set, raw syncmer
    /// matches are first plane-sweep-chained (same algorithm as
    /// `query_transitive_ext`), then only chains meeting both thresholds
    /// survive. Weaker/paralog-noise chains are dropped before partition's
    /// union-find sees them.
    pub fn with_chain_filter(mut self, min_chain_anchors: usize, min_chain_fraction: f64) -> Self {
        self.min_chain_anchors = Some(min_chain_anchors.max(1));
        self.min_chain_fraction = min_chain_fraction.clamp(0.0, 1.0);
        self
    }

    /// Borrow the inner `SyngIndex` — for callers that need full
    /// `HomologousInterval` results (strand, genome name) rather than just
    /// coitrees intervals.
    pub fn syng_index(&self) -> &syng::SyngIndex {
        &self.syng_index
    }

    pub fn syng_padding(&self) -> u64 {
        self.syng_padding
    }

    /// Convert syng query_region results into AdjustedInterval format.
    ///
    /// If chain filtering is enabled (via [`Self::with_chain_filter`]),
    /// raw syncmer matches are plane-sweep-chained and filtered by anchor
    /// count / query-coverage fraction before being emitted — otherwise the
    /// raw per-syncmer intervals flow through unchanged.
    fn query_via_syng(&self, target_id: u32, range_start: i32, range_end: i32) -> Vec<impg::AdjustedInterval> {
        let name = match self.seq_index.get_name(target_id) {
            Some(n) => n,
            None => return vec![],
        };
        match self.min_chain_anchors {
            None => self.query_via_syng_raw(name, range_start, range_end),
            Some(min_anchors) => self.query_via_syng_filtered(
                name,
                range_start,
                range_end,
                min_anchors,
                self.min_chain_fraction,
            ),
        }
    }

    fn query_via_syng_raw(
        &self,
        name: &str,
        range_start: i32,
        range_end: i32,
    ) -> Vec<impg::AdjustedInterval> {
        let intervals = match self.syng_index.query_region(
            name,
            range_start as u64,
            range_end as u64,
            self.syng_padding,
        ) {
            Ok(ivs) => ivs,
            Err(e) => {
                log::warn!("syng query_region failed for {}:{}-{}: {}", name, range_start, range_end, e);
                return vec![];
            }
        };
        intervals
            .into_iter()
            .filter_map(|iv| {
                let id = self.seq_index.get_id(&iv.genome)?;
                let interval = coitrees::Interval::new(iv.start as i32, iv.end as i32, id);
                Some((interval, vec![], interval))
            })
            .collect()
    }

    fn query_via_syng_filtered(
        &self,
        name: &str,
        range_start: i32,
        range_end: i32,
        min_anchors: usize,
        min_fraction: f64,
    ) -> Vec<impg::AdjustedInterval> {
        let hits = match self.syng_index.query_region_with_anchors_ext(
            name,
            range_start as u64,
            range_end as u64,
            self.syng_padding,
            0,
        ) {
            Ok(h) => h,
            Err(e) => {
                log::warn!("syng query_region_with_anchors failed for {}:{}-{}: {}", name, range_start, range_end, e);
                return vec![];
            }
        };
        let syncmer_len = (self.syng_index.params.w + self.syng_index.params.k) as u64;
        // Reuse the same extend_budget default as `query_transitive_ext`.
        let chained = syng_transitive::chain_anchors(
            hits,
            syng_transitive::DEFAULT_EXTEND_BUDGET_BP,
            syncmer_len,
        );
        let query_range_len = (range_end - range_start).max(0) as u64;
        let min_extent_bp = (query_range_len as f64 * min_fraction.max(0.0)) as u64;
        chained
            .into_iter()
            .filter(|c| {
                if c.anchors.len() < min_anchors {
                    return false;
                }
                if min_extent_bp == 0 {
                    return true;
                }
                let mut qmin = u64::MAX;
                let mut qmax = 0u64;
                for a in &c.anchors {
                    qmin = qmin.min(a.query_pos);
                    qmax = qmax.max(a.query_pos);
                }
                qmax.saturating_sub(qmin).saturating_add(syncmer_len) >= min_extent_bp
            })
            .filter_map(|c| {
                let tid = self.seq_index.get_id(&c.genome)?;
                // Match raw-path AdjustedInterval shape: both fields are the
                // target (homolog) interval. partition uses the first element
                // as "the other genome's region that belongs in this
                // partition" — which is the homolog, not the query region.
                let t_iv = coitrees::Interval::new(c.start as i32, c.end as i32, tid);
                Some((t_iv, Vec::<impg::CigarOp>::new(), t_iv))
            })
            .collect()
    }
}

impl impg_index::ImpgIndex for SyngImpgWrapper {
    fn seq_index(&self) -> &seqidx::SequenceIndex {
        &self.seq_index
    }

    fn query(
        &self, target_id: u32, range_start: i32, range_end: i32, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>, _: bool,
    ) -> Vec<impg::AdjustedInterval> {
        self.query_via_syng(target_id, range_start, range_end)
    }

    // CIGAR-cache methods are PAF-native concepts; syng has no cached
    // alignment storage, so they no-op. Reshaped from `unimplemented!()`
    // panics so `SyngImpgWrapper` can safely live inside a heterogeneous
    // `Vec<Arc<dyn ImpgIndex + Send + Sync>>` alongside alignment backends.
    fn query_with_cache(
        &self, target_id: u32, range_start: i32, range_end: i32,
        _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>,
        _: &rustc_hash::FxHashMap<(u32, u64), Vec<impg::CigarOp>>,
    ) -> Vec<impg::AdjustedInterval> {
        self.query_via_syng(target_id, range_start, range_end)
    }

    fn populate_cigar_cache(
        &self, _: u32, _: i32, _: i32, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>,
        _: &mut rustc_hash::FxHashMap<(u32, u64), Vec<impg::CigarOp>>,
    ) {
    }

    fn query_transitive_dfs(
        &self, target_id: u32, range_start: i32, range_end: i32,
        _: Option<&rustc_hash::FxHashMap<u32, impg::SortedRanges>>,
        _: u16, _: i32, _: i32, _: Option<i32>, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>, _: bool,
        _: Option<&subset_filter::SubsetFilter>,
    ) -> Vec<impg::AdjustedInterval> {
        self.query_via_syng(target_id, range_start, range_end)
    }

    fn query_transitive_bfs(
        &self, target_id: u32, range_start: i32, range_end: i32,
        _: Option<&rustc_hash::FxHashMap<u32, impg::SortedRanges>>,
        _: u16, _: i32, _: i32, _: Option<i32>, _: bool, _: Option<f64>,
        _: Option<&sequence_index::UnifiedSequenceIndex>, _: bool,
        _: Option<&subset_filter::SubsetFilter>,
    ) -> Vec<impg::AdjustedInterval> {
        self.query_via_syng(target_id, range_start, range_end)
    }

    fn get_or_load_tree(&self, _: u32) -> Option<std::sync::Arc<coitrees::BasicCOITree<impg::QueryMetadata, u32>>> {
        // Syng has no alignment interval tree. Returning None is correct —
        // callers that need this (partition / similarity / stats on PAF)
        // should dispatch via `query()` instead.
        None
    }

    fn target_ids(&self) -> Vec<u32> {
        // Every syng-indexed path is a potential target. Map through this
        // wrapper's sequence index so IDs match the caller's namespace.
        self.syng_index
            .name_map
            .path_to_name
            .iter()
            .filter_map(|name| self.seq_index.get_id(name))
            .collect()
    }

    fn remove_cached_tree(&self, _: u32) {
    }

    fn sequence_files(&self) -> &[String] {
        &[]
    }

    fn syng_index_ref(&self) -> Option<&syng::SyngIndex> {
        Some(&self.syng_index)
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
        GfaEngine::SyngNative => {
            // v2.3: if we have access to the underlying syng index,
            // re-query from the first interval as seed and use
            // anchor-seeded gap-only BiWFA for pairs that share anchors.
            // Otherwise fall back to the v1 full-pair path.
            let sequences = graph::prepare_sequences(impg, query_intervals, sequence_index)?;
            if sequences.is_empty() {
                return Ok(String::from("H\tVN:Z:1.0\n"));
            }
            let seq_pairs: Vec<(String, Vec<u8>)> = sequences
                .iter()
                .map(|(seq, meta)| (meta.path_name().to_string(), seq.as_bytes().to_vec()))
                .collect();

            if let Some(syng_idx) = impg.syng_index_ref() {
                // Build Member metadata alongside sequences for the
                // anchor-seeded driver.
                let members: Vec<(String, Vec<u8>, syng_graph::Member)> = sequences
                    .iter()
                    .map(|(seq, meta)| {
                        let (fwd_start, fwd_end) = if meta.strand == '+' {
                            (meta.start as u64, (meta.start + meta.size) as u64)
                        } else {
                            (
                                (meta.total_length as i32 - meta.start - meta.size) as u64,
                                (meta.total_length as i32 - meta.start) as u64,
                            )
                        };
                        (
                            meta.path_name().to_string(),
                            seq.as_bytes().to_vec(),
                            syng_graph::Member {
                                chrom: meta.name.clone(),
                                fwd_start,
                                fwd_end,
                                strand: meta.strand,
                            },
                        )
                    })
                    .collect();
                // Pick first interval as seed. (The partition's seed
                // isn't tracked explicitly; using the first element
                // matches the typical greedy-expansion order.)
                let seed = &members[0].2;
                let paf = syng_graph::build_paf_anchor_seeded(
                    &members,
                    &seed.chrom,
                    seed.fwd_start,
                    seed.fwd_end,
                    syng_idx,
                    /* syng_padding */ 120,
                    /* max_depth */ 1,
                    /* k_near */ 3,
                    /* k_far */ 1,
                    /* random_fraction */ 0.01,
                );
                let t_induce = std::time::Instant::now();
                let gfa = syng_graph::build_gfa_from_paf_and_sequences(
                    &seq_pairs,
                    &paf,
                    &engine_opts.pipeline,
                )?;
                log::info!(
                    "syng_graph: seqwish induction {:.2}s ({} bytes PAF, {} sequences)",
                    t_induce.elapsed().as_secs_f64(),
                    paf.len(),
                    seq_pairs.len(),
                );
                if skip_normalize {
                    Ok(gfa)
                } else {
                    graph::normalize_and_sort(gfa, num_threads)
                }
            } else {
                let gfa = syng_graph::build_gfa_syng_native_from_sequences(
                    &seq_pairs,
                    &engine_opts.pipeline,
                )?;
                if skip_normalize {
                    Ok(gfa)
                } else {
                    graph::normalize_and_sort(gfa, num_threads)
                }
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

    // Compute per-partition bp once for deterministic output ordering below.
    let per_partition_bp: Vec<u64> = partitions
        .iter()
        .map(|ivs| {
            ivs.iter()
                .map(|iv| (iv.last - iv.first).unsigned_abs() as u64)
                .sum()
        })
        .collect();

    // 1. Generate per-partition GFAs with two-pool parallelism.
    //
    // Problem: naive nested `par_iter` thrashes. If the outer loop is a
    // `par_iter` on the global pool AND each inner step is also a
    // `par_iter` on the global pool, every outer task spawns inner
    // tasks that compete with the other outer tasks for the same 16
    // threads. On yeast235 this took wall from ~3s/partition solo to
    // >24 min/partition with zero completions over 24 min.
    //
    // Fix (well-documented rayon pattern; see Piotr Kołaczkowski /
    // Daniel Imfeld posts, and our own `smooth.rs:249` +
    // `align.rs:1194`): two *separate* thread pools. Outer pool drives
    // `par_iter` over partitions. Inside each partition task we
    // `inner_pool.install(|| ...)` so every inner `par_iter` routes to
    // the inner pool. No shared queue → no contention.
    //
    // Inner pool gets ALL cores; outer pool is tiny (4 threads) since
    // each outer thread mostly blocks on `install` waiting for its
    // partition's inner par_iter to drain on the shared inner pool.
    // Multiple outer threads submitting `install` concurrently is fine:
    // rayon steals pair-tasks from all of them into inner's global queue,
    // keeping all inner threads busy across partitions.
    use rayon::prelude::*;
    use std::sync::Arc;
    let total_threads = engine_opts.pipeline.num_threads.max(1);
    let outer_threads = 4.min(total_threads);
    let inner_threads = total_threads;
    // 16 MB stack: default rayon stack (~2 MB on Linux) is too small for
    // our nested seqwish + BiWFA call graph (observed stack overflow on
    // yeast235). Matches the cost of ~16 threads × 16 MB = 256 MB
    // reserved, which is negligible at our scale.
    let outer_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(outer_threads)
        .stack_size(16 * 1024 * 1024)
        .build()
        .map_err(|e| std::io::Error::other(format!("outer pool: {}", e)))?;
    let inner_pool = Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(inner_threads)
            .stack_size(16 * 1024 * 1024)
            .build()
            .map_err(|e| std::io::Error::other(format!("inner pool: {}", e)))?,
    );
    info!(
        "[partitioned] Two-pool scheduling: {} outer × {} inner = {} total threads",
        outer_threads, inner_threads, outer_threads * inner_threads,
    );

    let total_partitions = partitions.len();
    let completed = std::sync::atomic::AtomicUsize::new(0);
    let pipeline_start = Instant::now();

    let sub_gfas: Vec<String> = outer_pool.install(|| -> std::io::Result<Vec<String>> {
        partitions
            .par_iter()
            .enumerate()
            .map(|(idx, intervals)| -> std::io::Result<String> {
                let num_regions = intervals.len();
                let current_partition_length = per_partition_bp[idx];
                let part_start = Instant::now();
                let inner_pool = Arc::clone(&inner_pool);
                let gfa = inner_pool.install(|| {
                    dispatch_gfa_engine_inner(
                        impg,
                        intervals,
                        sequence_index,
                        scoring_params,
                        engine_opts,
                    )
                })?;
                let part_elapsed = part_start.elapsed();
                let done = completed.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                let pipeline_elapsed = pipeline_start.elapsed();
                info!(
                    "[partitioned] Completed partition {}/{} ({} intervals, {} bp) in {:.1}s (pipeline elapsed: {:.1}s)",
                    done, total_partitions, num_regions, current_partition_length,
                    part_elapsed.as_secs_f64(),
                    pipeline_elapsed.as_secs_f64(),
                );
                Ok(gfa)
            })
            .collect()
    })?;
    let total_partitioned_bp: u64 = per_partition_bp.iter().sum();
    let _ = total_partitioned_bp;

    // 2. Lace all partition GFAs together
    info!("[partitioned] Lacing {} partition GFAs", sub_gfas.len());
    let laced = commands::lace::lace_subgraphs(&sub_gfas, None)?;
    // Free all sub-partition GFAs before gfaffix/sort
    drop(sub_gfas);

    // 3. Single final gfaffix normalization + sort
    info!("[partitioned] Running final gfaffix normalization");
    graph::normalize_and_sort(laced, engine_opts.pipeline.num_threads)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::impg_index::ImpgIndex;

    // Syng C library has non-thread-safe global state; serialize all FFI tests.
    static SYNG_LOCK: std::sync::LazyLock<std::sync::Mutex<()>> =
        std::sync::LazyLock::new(|| std::sync::Mutex::new(()));
    fn lock_syng() -> std::sync::MutexGuard<'static, ()> {
        SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner())
    }

    fn make_test_sequence(len: usize, seed: u8) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut seq = Vec::with_capacity(len);
        let mut state: u32 = seed as u32;
        for _ in 0..len {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            seq.push(bases[((state >> 16) % 4) as usize]);
        }
        seq
    }

    /// Build a SyngIndex with shared-backbone sequences for testing.
    fn build_test_syng() -> (syng::SyngIndex, seqidx::SequenceIndex) {
        let params = syng::SyncmerParams { k: 8, w: 55, seed: 7 };
        let shared_len = 500;
        let total_len = 1000;
        let backbone = make_test_sequence(shared_len, 42);

        let mut seq_a = backbone.clone();
        seq_a.extend_from_slice(&make_test_sequence(total_len - shared_len, 1));

        let mut seq_b = backbone.clone();
        seq_b.extend_from_slice(&make_test_sequence(total_len - shared_len, 2));

        let seq_c = make_test_sequence(total_len, 99);

        let sequences = vec![
            ("genome_a".to_string(), seq_a),
            ("genome_b".to_string(), seq_b),
            ("genome_c".to_string(), seq_c),
        ];

        let index = syng::SyngIndex::build(params, sequences.into_iter());
        let seq_index = index.build_seq_index();
        (index, seq_index)
    }

    #[test]
    fn test_syng_impg_wrapper_seq_index() {
        let _guard = lock_syng();
        let (syng_index, seq_index) = build_test_syng();
        let wrapper = SyngImpgWrapper::new(syng_index, seq_index, 0);

        // All three genomes should be in the seq_index
        assert_eq!(wrapper.seq_index().len(), 3);
        assert!(wrapper.seq_index().get_id("genome_a").is_some());
        assert!(wrapper.seq_index().get_id("genome_b").is_some());
        assert!(wrapper.seq_index().get_id("genome_c").is_some());
    }

    #[test]
    fn test_syng_impg_wrapper_query_transitive_dfs() {
        let _guard = lock_syng();
        let (syng_index, seq_index) = build_test_syng();
        let wrapper = SyngImpgWrapper::new(syng_index, seq_index, 0);

        let target_id = wrapper.seq_index().get_id("genome_a").unwrap();

        // Query the shared backbone region via the ImpgIndex trait
        let results = wrapper.query_transitive_dfs(
            target_id, 0, 500, None, 1, 0, 0, None, false, None, None, false, None,
        );

        // Should find at least genome_a (self) and genome_b (shared backbone)
        let genomes_found: Vec<&str> = results
            .iter()
            .filter_map(|(iv, _, _)| wrapper.seq_index().get_name(iv.metadata))
            .collect();
        assert!(
            genomes_found.contains(&"genome_a"),
            "should find self-hit; found: {:?}",
            genomes_found
        );
        assert!(
            genomes_found.contains(&"genome_b"),
            "should find genome_b (shared backbone); found: {:?}",
            genomes_found
        );

        // All intervals should have valid coordinates
        for (iv, cigar, target_iv) in &results {
            assert!(iv.last > iv.first, "interval end > start");
            assert!(cigar.is_empty(), "syng path produces no CIGAR ops");
            assert_eq!(iv.first, target_iv.first, "query and target intervals should match");
            assert_eq!(iv.last, target_iv.last);
            assert_eq!(iv.metadata, target_iv.metadata);
        }
    }

    #[test]
    fn test_syng_impg_wrapper_query_bfs_same_as_dfs() {
        let _guard = lock_syng();
        let (syng_index, seq_index) = build_test_syng();
        let wrapper = SyngImpgWrapper::new(syng_index, seq_index, 0);

        let target_id = wrapper.seq_index().get_id("genome_a").unwrap();

        let dfs = wrapper.query_transitive_dfs(
            target_id, 0, 500, None, 1, 0, 0, None, false, None, None, false, None,
        );
        let bfs = wrapper.query_transitive_bfs(
            target_id, 0, 500, None, 1, 0, 0, None, false, None, None, false, None,
        );

        // Both should return the same results since they both call query_region
        assert_eq!(dfs.len(), bfs.len());
        for (d, b) in dfs.iter().zip(bfs.iter()) {
            assert_eq!(d.0.first, b.0.first);
            assert_eq!(d.0.last, b.0.last);
            assert_eq!(d.0.metadata, b.0.metadata);
        }
    }

    #[test]
    fn test_syng_impg_wrapper_unknown_target() {
        let _guard = lock_syng();
        let (syng_index, seq_index) = build_test_syng();
        let wrapper = SyngImpgWrapper::new(syng_index, seq_index, 0);

        // Query with an invalid target_id should return empty (not panic)
        let results = wrapper.query_transitive_dfs(
            999, 0, 500, None, 1, 0, 0, None, false, None, None, false, None,
        );
        assert!(results.is_empty());
    }

    #[test]
    fn test_syng_impg_wrapper_with_padding() {
        let _guard = lock_syng();
        let (syng_index, seq_index) = build_test_syng();
        let no_pad = SyngImpgWrapper::new(syng_index, seq_index, 0);

        // Re-build for the padded wrapper (SyngIndex is consumed by first wrapper)
        let (syng_index2, seq_index2) = build_test_syng();
        let with_pad = SyngImpgWrapper::new(syng_index2, seq_index2, 120);

        let target_a = no_pad.seq_index().get_id("genome_a").unwrap();

        let results_no_pad = no_pad.query_transitive_dfs(
            target_a, 0, 500, None, 1, 0, 0, None, false, None, None, false, None,
        );
        let results_with_pad = with_pad.query_transitive_dfs(
            target_a, 0, 500, None, 1, 0, 0, None, false, None, None, false, None,
        );

        // With padding, intervals should be at least as wide
        for (no_p, with_p) in results_no_pad.iter().zip(results_with_pad.iter()) {
            if no_p.0.metadata == with_p.0.metadata {
                assert!(
                    (with_p.0.last - with_p.0.first) >= (no_p.0.last - no_p.0.first),
                    "padded intervals should be >= unpadded"
                );
            }
        }
    }
}
