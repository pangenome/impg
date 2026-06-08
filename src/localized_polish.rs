//! Iterative localized graph polishing over SYNG-collected local sequences.
//!
//! This layer composes the explicit local seed graph builder, diagnostic dirty
//! region detector, and localized resolver. Graph-quality metrics are logged
//! and reported for diagnosis only; exact path-name and path-spelling
//! preservation is the hard correctness gate.

use crate::graph_badness::{
    self, DirtyChunk, DirtyRegionOptions, DirtyRegionReport, GraphBadnessMetrics,
};
use crate::local_seed::{
    induce_seed_graph, validate_seed_path_sequences, LocalSeedGraph, LocalSeedInductionConfig,
    LocalSeedReport, LocalSequenceSet,
};
use crate::resolution::{
    resolve_gfa_dirty_chunks, LocalizedResolutionStatus, LocalizedResolverConfig,
    LocalizedResolverReport,
};
use crate::syng;
use std::fmt::Write as _;
use std::io;
use std::path::Path;
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct LocalizedPolishConfig {
    pub max_iterations: usize,
    pub max_chunks_per_iteration: usize,
    pub max_total_chunks: usize,
    pub max_total_chunk_bp: Option<usize>,
    pub max_runtime_secs: Option<f64>,
    pub dirty_options: DirtyRegionOptions,
    pub resolver: LocalizedResolverConfig,
    pub debug_dir: Option<String>,
}

impl Default for LocalizedPolishConfig {
    fn default() -> Self {
        Self {
            max_iterations: 3,
            max_chunks_per_iteration: 1,
            max_total_chunks: 16,
            max_total_chunk_bp: Some(500_000),
            max_runtime_secs: None,
            dirty_options: DirtyRegionOptions::default(),
            resolver: LocalizedResolverConfig::default(),
            debug_dir: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct LocalizedPolishResult {
    pub gfa: String,
    pub report: LocalizedPolishReport,
}

#[derive(Clone, Debug)]
pub struct LocalizedPolishReport {
    pub seed: LocalSeedReport,
    pub final_status: String,
    pub final_path_validation: String,
    pub hard_rejection_gate: String,
    pub metrics_are_diagnostic_only: bool,
    pub iterations: Vec<LocalizedPolishIterationReport>,
    pub total_chunks_attempted: usize,
    pub total_chunks_applied: usize,
    pub total_elapsed_secs: f64,
}

#[derive(Clone, Debug)]
pub struct LocalizedPolishIterationReport {
    pub iteration: usize,
    pub status: String,
    pub dirty_sites: usize,
    pub candidate_chunks: usize,
    pub selected_chunk_ids: Vec<String>,
    pub deferred_chunks: usize,
    pub selected_chunk_bp: usize,
    pub metrics_before: GraphBadnessMetrics,
    pub metrics_after: Option<GraphBadnessMetrics>,
    pub resolver_reports: Vec<LocalizedResolverReport>,
    pub elapsed_secs: f64,
}

/// Build a local seed graph from a SYNG-collected sequence set, then iteratively
/// polish detector-selected dirty regions until clean or budget-limited.
pub fn polish_local_sequences(
    sequence_set: &LocalSequenceSet,
    seed_config: &LocalSeedInductionConfig,
    syng_index: Option<&syng::SyngIndex>,
    config: &LocalizedPolishConfig,
) -> io::Result<LocalizedPolishResult> {
    validate_config(config)?;
    let seed = induce_seed_graph(sequence_set, seed_config, syng_index)?;
    polish_seed_graph(sequence_set, seed, config)
}

/// Run the iterative loop on an already-built seed graph.
///
/// This is primarily useful for tests and reproducibility. The production
/// `syng-local:localized` path calls [`polish_local_sequences`] so seed
/// induction is still part of the runnable output path.
pub fn polish_seed_graph(
    sequence_set: &LocalSequenceSet,
    seed: LocalSeedGraph,
    config: &LocalizedPolishConfig,
) -> io::Result<LocalizedPolishResult> {
    validate_config(config)?;
    validate_seed_path_sequences(sequence_set, &seed.gfa).map_err(|err| {
        path_corruption_error(format!(
            "localized polish seed path validation failed: {err}"
        ))
    })?;

    let run_start = Instant::now();
    let mut gfa = seed.gfa;
    let mut iterations = Vec::new();
    let mut final_status = String::from("iteration-budget-exhausted");
    let mut total_chunks_attempted = 0usize;
    let mut total_chunks_applied = 0usize;
    let mut total_chunk_bp = 0usize;

    for iteration in 0..config.max_iterations {
        if runtime_exhausted(run_start, config.max_runtime_secs) {
            final_status = String::from("runtime-budget-exhausted");
            break;
        }

        let iteration_start = Instant::now();
        let source = format!("localized-polish.iter{iteration}");
        let dirty_report = graph_badness::analyze_gfa(source, &gfa, &config.dirty_options)?;
        log_dirty_summary(iteration, &dirty_report);
        write_dirty_debug(config, iteration, &dirty_report);

        if dirty_report.candidate_chunks.is_empty() {
            final_status = String::from("converged-clean");
            iterations.push(LocalizedPolishIterationReport {
                iteration,
                status: final_status.clone(),
                dirty_sites: dirty_report.dirty_sites.len(),
                candidate_chunks: dirty_report.candidate_chunks.len(),
                selected_chunk_ids: Vec::new(),
                deferred_chunks: 0,
                selected_chunk_bp: 0,
                metrics_before: dirty_report.metrics,
                metrics_after: None,
                resolver_reports: Vec::new(),
                elapsed_secs: iteration_start.elapsed().as_secs_f64(),
            });
            break;
        }

        let selection = select_chunks(
            &dirty_report.candidate_chunks,
            config,
            total_chunks_attempted,
            total_chunk_bp,
        );

        if selection.chunks.is_empty() {
            final_status = String::from("chunk-budget-exhausted");
            iterations.push(LocalizedPolishIterationReport {
                iteration,
                status: final_status.clone(),
                dirty_sites: dirty_report.dirty_sites.len(),
                candidate_chunks: dirty_report.candidate_chunks.len(),
                selected_chunk_ids: Vec::new(),
                deferred_chunks: dirty_report.candidate_chunks.len(),
                selected_chunk_bp: 0,
                metrics_before: dirty_report.metrics,
                metrics_after: None,
                resolver_reports: Vec::new(),
                elapsed_secs: iteration_start.elapsed().as_secs_f64(),
            });
            break;
        }

        total_chunks_attempted += selection.chunks.len();
        total_chunk_bp = total_chunk_bp.saturating_add(selection.selected_bp);

        let selected_ids = selection
            .chunks
            .iter()
            .map(|chunk| chunk.id.clone())
            .collect::<Vec<_>>();
        log::info!(
            "localized polish: iteration={} selecting {} / {} dirty chunk(s), selected_bp={}, deferred={}, hard_gate=exact_path_corruption, metrics_gate=false",
            iteration,
            selection.chunks.len(),
            dirty_report.candidate_chunks.len(),
            selection.selected_bp,
            selection.deferred,
        );

        let resolved = resolve_gfa_dirty_chunks(&gfa, &selection.chunks, &config.resolver)?;
        let path_invalid_reports = resolved
            .reports
            .iter()
            .filter(|report| report.status == LocalizedResolutionStatus::PathInvalid)
            .map(|report| report.chunk_id.clone())
            .collect::<Vec<_>>();
        if !path_invalid_reports.is_empty() {
            return Err(path_corruption_error(format!(
                "localized polish exact path corruption in iteration {iteration}: chunks={path_invalid_reports:?}"
            )));
        }

        validate_seed_path_sequences(sequence_set, &resolved.gfa).map_err(|err| {
            path_corruption_error(format!(
                "localized polish exact path corruption after iteration {iteration}: {err}"
            ))
        })?;

        let applied = resolved
            .reports
            .iter()
            .filter(|report| report.status == LocalizedResolutionStatus::Applied)
            .count();
        total_chunks_applied += applied;
        let changed = resolved.gfa != gfa;
        let after_source = format!("localized-polish.iter{iteration}.after");
        let after_report =
            graph_badness::analyze_gfa(after_source, &resolved.gfa, &config.dirty_options)?;
        log_after_summary(iteration, applied, changed, &after_report);

        gfa = resolved.gfa;
        let status = if applied == 0 || !changed {
            final_status = String::from("stalled-no-change");
            final_status.clone()
        } else if total_chunks_attempted >= config.max_total_chunks {
            final_status = String::from("chunk-budget-exhausted");
            final_status.clone()
        } else {
            String::from("applied")
        };

        iterations.push(LocalizedPolishIterationReport {
            iteration,
            status: status.clone(),
            dirty_sites: dirty_report.dirty_sites.len(),
            candidate_chunks: dirty_report.candidate_chunks.len(),
            selected_chunk_ids: selected_ids,
            deferred_chunks: selection.deferred,
            selected_chunk_bp: selection.selected_bp,
            metrics_before: dirty_report.metrics,
            metrics_after: Some(after_report.metrics),
            resolver_reports: resolved.reports,
            elapsed_secs: iteration_start.elapsed().as_secs_f64(),
        });

        if status != "applied" {
            break;
        }
    }

    validate_seed_path_sequences(sequence_set, &gfa).map_err(|err| {
        path_corruption_error(format!(
            "localized polish final exact path validation failed: {err}"
        ))
    })?;

    let report = LocalizedPolishReport {
        seed: seed.report,
        final_status,
        final_path_validation: String::from("pass"),
        hard_rejection_gate: String::from("exact_path_corruption"),
        metrics_are_diagnostic_only: true,
        iterations,
        total_chunks_attempted,
        total_chunks_applied,
        total_elapsed_secs: run_start.elapsed().as_secs_f64(),
    };
    log::info!(
        "localized polish: status={} attempted={} applied={} iterations={} path_validation={} elapsed={:.3}s",
        report.final_status,
        report.total_chunks_attempted,
        report.total_chunks_applied,
        report.iterations.len(),
        report.final_path_validation,
        report.total_elapsed_secs,
    );
    write_summary_debug(config, &report);

    Ok(LocalizedPolishResult { gfa, report })
}

fn validate_config(config: &LocalizedPolishConfig) -> io::Result<()> {
    if config.max_iterations == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "localized polish max iterations must be > 0",
        ));
    }
    if config.max_chunks_per_iteration == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "localized polish max chunks per iteration must be > 0",
        ));
    }
    if config.max_total_chunks == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "localized polish max total chunks must be > 0",
        ));
    }
    if let Some(max_runtime_secs) = config.max_runtime_secs {
        if max_runtime_secs <= 0.0 || !max_runtime_secs.is_finite() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "localized polish max runtime seconds must be positive and finite",
            ));
        }
    }
    Ok(())
}

#[derive(Clone, Debug, Default)]
struct ChunkSelection {
    chunks: Vec<DirtyChunk>,
    selected_bp: usize,
    deferred: usize,
}

fn select_chunks(
    chunks: &[DirtyChunk],
    config: &LocalizedPolishConfig,
    total_chunks_attempted: usize,
    total_chunk_bp: usize,
) -> ChunkSelection {
    let mut selection = ChunkSelection::default();
    for chunk in chunks {
        if selection.chunks.len() >= config.max_chunks_per_iteration
            || total_chunks_attempted + selection.chunks.len() >= config.max_total_chunks
        {
            selection.deferred += 1;
            continue;
        }

        let chunk_bp = chunk_bp(chunk);
        if let Some(max_total_bp) = config.max_total_chunk_bp {
            if total_chunk_bp
                .saturating_add(selection.selected_bp)
                .saturating_add(chunk_bp)
                > max_total_bp
            {
                selection.deferred += 1;
                continue;
            }
        }

        selection.selected_bp = selection.selected_bp.saturating_add(chunk_bp);
        selection.chunks.push(chunk.clone());
    }
    selection
}

fn chunk_bp(chunk: &DirtyChunk) -> usize {
    chunk
        .path_start_bp
        .zip(chunk.path_end_bp)
        .map(|(start, end)| end.saturating_sub(start))
        .or_else(|| {
            chunk
                .graph_start_bp
                .zip(chunk.graph_end_bp)
                .map(|(start, end)| end.saturating_sub(start))
        })
        .unwrap_or(0)
}

fn runtime_exhausted(run_start: Instant, max_runtime_secs: Option<f64>) -> bool {
    max_runtime_secs.is_some_and(|limit| run_start.elapsed().as_secs_f64() >= limit)
}

fn path_corruption_error(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

fn log_dirty_summary(iteration: usize, report: &DirtyRegionReport) {
    let m = &report.metrics;
    log::info!(
        "localized polish: iteration={} dirty_sites={} candidate_chunks={} metrics_diagnostic_only={} segments={} links={} paths={} path_steps={} compression={} singleton_bp={} low_depth_bp={} white_space_bp_total={} white_space_bp_max={} small_loop_sites={}",
        iteration,
        report.dirty_sites.len(),
        report.candidate_chunks.len(),
        report.metrics_are_diagnostic_only,
        m.segment_count,
        m.link_count,
        m.path_count,
        m.total_path_steps,
        format_option_f64(m.path_replay_compression_ratio),
        m.singleton_bp,
        m.low_depth_bp,
        m.white_space_proxy_bp_total,
        m.white_space_proxy_bp_max,
        m.small_loop_sites,
    );
}

fn log_after_summary(
    iteration: usize,
    applied: usize,
    changed: bool,
    after_report: &DirtyRegionReport,
) {
    let m = &after_report.metrics;
    log::info!(
        "localized polish: iteration={} resolver_applied={} changed={} after_candidate_chunks={} after_segments={} after_links={} after_path_steps={} after_compression={} after_singleton_bp={} after_low_depth_bp={} after_white_space_bp_total={}",
        iteration,
        applied,
        changed,
        after_report.candidate_chunks.len(),
        m.segment_count,
        m.link_count,
        m.total_path_steps,
        format_option_f64(m.path_replay_compression_ratio),
        m.singleton_bp,
        m.low_depth_bp,
        m.white_space_proxy_bp_total,
    );
}

fn format_option_f64(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.6}"))
        .unwrap_or_else(|| "n/a".to_string())
}

fn write_dirty_debug(
    config: &LocalizedPolishConfig,
    iteration: usize,
    dirty_report: &DirtyRegionReport,
) {
    let Some(debug_dir) = config.debug_dir.as_deref() else {
        return;
    };
    if let Err(err) = std::fs::create_dir_all(debug_dir) {
        log::warn!("localized polish: failed to create debug dir {debug_dir}: {err}");
        return;
    }
    for format in ["json", "tsv"] {
        let text = match graph_badness::format_dirty_report(dirty_report, format) {
            Ok(text) => text,
            Err(err) => {
                log::warn!("localized polish: failed to format dirty report as {format}: {err}");
                continue;
            }
        };
        let path = Path::new(debug_dir).join(format!(
            "localized_polish_iter_{iteration:03}_dirty_regions.{format}"
        ));
        if let Err(err) = std::fs::write(&path, text) {
            log::warn!(
                "localized polish: failed to write dirty report {}: {err}",
                path.display()
            );
        }
    }
}

fn write_summary_debug(config: &LocalizedPolishConfig, report: &LocalizedPolishReport) {
    let Some(debug_dir) = config.debug_dir.as_deref() else {
        return;
    };
    if let Err(err) = std::fs::create_dir_all(debug_dir) {
        log::warn!("localized polish: failed to create debug dir {debug_dir}: {err}");
        return;
    }
    let path = Path::new(debug_dir).join("localized_polish_summary.tsv");
    if let Err(err) = std::fs::write(&path, summary_tsv(report)) {
        log::warn!(
            "localized polish: failed to write summary {}: {err}",
            path.display()
        );
    }
}

fn summary_tsv(report: &LocalizedPolishReport) -> String {
    let mut out = String::from(
        "record_type\titeration\tstatus\tseed_source\tpath_validation\tdirty_sites\tcandidate_chunks\tselected_chunks\tdeferred_chunks\tselected_chunk_bp\tapplied\tsegments_before\tsegments_after\tpath_steps_before\tpath_steps_after\twhite_space_bp_before\twhite_space_bp_after\tsingleton_bp_before\tsingleton_bp_after\telapsed_secs\n",
    );
    let run_row = vec![
        "run".to_string(),
        String::new(),
        report.final_status.clone(),
        report.seed.seed_source.clone(),
        report.final_path_validation.clone(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        report.total_chunks_applied.to_string(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        format!("{:.6}", report.total_elapsed_secs),
    ];
    out.push_str(&run_row.join("\t"));
    out.push('\n');
    for iteration in &report.iterations {
        let after = iteration.metrics_after.as_ref();
        let applied = iteration
            .resolver_reports
            .iter()
            .filter(|resolver| resolver.status == LocalizedResolutionStatus::Applied)
            .count();
        let _ = writeln!(
            out,
            "iteration\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
            iteration.iteration,
            iteration.status,
            report.seed.seed_source,
            report.final_path_validation,
            iteration.dirty_sites,
            iteration.candidate_chunks,
            iteration.selected_chunk_ids.join(","),
            iteration.deferred_chunks,
            iteration.selected_chunk_bp,
            applied,
            iteration.metrics_before.segment_count,
            after.map(|m| m.segment_count).map(|v| v.to_string()).unwrap_or_default(),
            iteration.metrics_before.total_path_steps,
            after.map(|m| m.total_path_steps).map(|v| v.to_string()).unwrap_or_default(),
            iteration.metrics_before.white_space_proxy_bp_total,
            after
                .map(|m| m.white_space_proxy_bp_total)
                .map(|v| v.to_string())
                .unwrap_or_default(),
            iteration.metrics_before.singleton_bp,
            after.map(|m| m.singleton_bp).map(|v| v.to_string()).unwrap_or_default(),
            iteration.elapsed_secs,
        );
    }
    out
}
