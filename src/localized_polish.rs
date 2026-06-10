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
use std::collections::{BTreeMap, BTreeSet};
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
        let mut resolver = LocalizedResolverConfig::default();
        // Localized pairwise polishing should not silently inherit the global
        // crush 1:1 replacement filter. C4 dirty blocks are many-haplotype
        // regions; if a user wants a strict replacement filter they can still
        // request it with replacement-num-mappings/replacement-scaffold-filter.
        resolver.resolution.replacement_num_mappings = "many:many".to_string();
        resolver.resolution.replacement_scaffold_filter = "many:many".to_string();
        resolver.resolution.replacement_scaffold_mass = 0;
        Self {
            max_iterations: 3,
            max_chunks_per_iteration: 1,
            max_total_chunks: 16,
            max_total_chunk_bp: Some(500_000),
            max_runtime_secs: None,
            dirty_options: DirtyRegionOptions::default(),
            resolver,
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
    let mut polished_provenance = DirtyChunkProvenance::default();

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

        let selection = select_chunks_with_provenance(
            &dirty_report.candidate_chunks,
            config,
            total_chunks_attempted,
            total_chunk_bp,
            &polished_provenance,
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
        for chunk in &selection.chunks {
            if resolved.reports.iter().any(|report| {
                report.chunk_id == chunk.id && report.status == LocalizedResolutionStatus::Applied
            }) {
                polished_provenance.remember(chunk);
            }
        }
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

#[derive(Clone, Debug)]
struct CompoundDirtyChunk {
    first_index: usize,
    source_count: usize,
    chunk: DirtyChunk,
}

#[derive(Clone, Debug, Default)]
struct DirtyChunkProvenance {
    nodes: BTreeSet<String>,
    path_intervals: BTreeMap<String, Vec<(usize, usize)>>,
}

#[cfg(test)]
fn select_chunks(
    chunks: &[DirtyChunk],
    config: &LocalizedPolishConfig,
    total_chunks_attempted: usize,
    total_chunk_bp: usize,
) -> ChunkSelection {
    select_chunks_with_provenance(
        chunks,
        config,
        total_chunks_attempted,
        total_chunk_bp,
        &DirtyChunkProvenance::default(),
    )
}

fn select_chunks_with_provenance(
    chunks: &[DirtyChunk],
    config: &LocalizedPolishConfig,
    total_chunks_attempted: usize,
    total_chunk_bp: usize,
    polished_provenance: &DirtyChunkProvenance,
) -> ChunkSelection {
    let mut selection = ChunkSelection::default();
    let mut selected_provenance = polished_provenance.clone();
    let mut edge_compounds = Vec::new();

    for compound in compound_dirty_chunks(chunks, config) {
        if chunk_core_has_internal_path_anchors(&compound.chunk) {
            consider_compound_for_selection(
                compound,
                config,
                total_chunks_attempted,
                total_chunk_bp,
                &mut selection,
                &mut selected_provenance,
            );
        } else {
            edge_compounds.push(compound);
        }
    }

    for compound in edge_compounds {
        consider_compound_for_selection(
            compound,
            config,
            total_chunks_attempted,
            total_chunk_bp,
            &mut selection,
            &mut selected_provenance,
        );
    }
    selection
}

fn consider_compound_for_selection(
    compound: CompoundDirtyChunk,
    config: &LocalizedPolishConfig,
    total_chunks_attempted: usize,
    total_chunk_bp: usize,
    selection: &mut ChunkSelection,
    selected_provenance: &mut DirtyChunkProvenance,
) {
    let chunk = &compound.chunk;
    if selected_provenance.covers(chunk) {
        selection.deferred = selection.deferred.saturating_add(compound.source_count);
        return;
    }

    if selection.chunks.len() >= config.max_chunks_per_iteration
        || total_chunks_attempted + selection.chunks.len() >= config.max_total_chunks
    {
        selection.deferred = selection.deferred.saturating_add(compound.source_count);
        return;
    }

    let chunk_bp = chunk_bp(chunk);
    if let Some(max_total_bp) = config.max_total_chunk_bp {
        if total_chunk_bp
            .saturating_add(selection.selected_bp)
            .saturating_add(chunk_bp)
            > max_total_bp
        {
            selection.deferred = selection.deferred.saturating_add(compound.source_count);
            return;
        }
    }

    selection.selected_bp = selection.selected_bp.saturating_add(chunk_bp);
    selection.chunks.push(chunk.clone());
    selected_provenance.remember(chunk);
}

fn chunk_core_has_internal_path_anchors(chunk: &DirtyChunk) -> bool {
    let Some(path_len) = chunk.path_length_bp else {
        return true;
    };
    let Some((start, end)) = chunk_core_span(chunk) else {
        return true;
    };
    start > 0 && end < path_len
}

fn compound_dirty_chunks(
    chunks: &[DirtyChunk],
    config: &LocalizedPolishConfig,
) -> Vec<CompoundDirtyChunk> {
    let mut path_groups: BTreeMap<String, Vec<(usize, DirtyChunk)>> = BTreeMap::new();
    let mut compounds = Vec::new();

    for (index, chunk) in chunks.iter().cloned().enumerate() {
        if let Some(path_name) = chunk.path_name.clone() {
            if chunk_path_span(&chunk).is_some() {
                path_groups
                    .entry(path_name)
                    .or_default()
                    .push((index, chunk));
                continue;
            }
        }
        compounds.push(CompoundDirtyChunk {
            first_index: index,
            source_count: 1,
            chunk,
        });
    }

    let merge_distance_bp = config
        .dirty_options
        .merge_distance_bp
        .max(config.dirty_options.flank_bp);
    for (_path_name, mut group) in path_groups {
        group.sort_by(|left, right| {
            chunk_path_span(&left.1)
                .map(|(start, end)| (start, end, left.0))
                .cmp(&chunk_path_span(&right.1).map(|(start, end)| (start, end, right.0)))
        });

        let mut acc: Vec<(usize, DirtyChunk)> = Vec::new();
        let mut acc_end = 0usize;
        for item in group {
            let Some((next_start, next_end)) = chunk_path_span(&item.1) else {
                if !acc.is_empty() {
                    push_compound_dirty_chunk(&mut compounds, std::mem::take(&mut acc));
                }
                compounds.push(CompoundDirtyChunk {
                    first_index: item.0,
                    source_count: 1,
                    chunk: item.1,
                });
                continue;
            };

            if acc.is_empty() {
                acc_end = next_end;
                acc.push(item);
                continue;
            }

            if next_start <= acc_end.saturating_add(merge_distance_bp) {
                acc_end = acc_end.max(next_end);
                acc.push(item);
            } else {
                push_compound_dirty_chunk(&mut compounds, std::mem::take(&mut acc));
                acc_end = next_end;
                acc.push(item);
            }
        }
        if !acc.is_empty() {
            push_compound_dirty_chunk(&mut compounds, acc);
        }
    }

    compounds.sort_by_key(|compound| compound.first_index);
    compounds
}

fn push_compound_dirty_chunk(
    compounds: &mut Vec<CompoundDirtyChunk>,
    chunks: Vec<(usize, DirtyChunk)>,
) {
    if chunks.is_empty() {
        return;
    }
    let first_index = chunks.iter().map(|(index, _)| *index).min().unwrap_or(0);
    let source_count = chunks.len();
    let chunk = if source_count == 1 {
        chunks.into_iter().next().map(|(_, chunk)| chunk).unwrap()
    } else {
        merge_dirty_chunks(&chunks)
    };
    compounds.push(CompoundDirtyChunk {
        first_index,
        source_count,
        chunk,
    });
}

fn merge_dirty_chunks(chunks: &[(usize, DirtyChunk)]) -> DirtyChunk {
    let first = &chunks[0].1;
    let last = &chunks[chunks.len() - 1].1;
    let mut site_ids = Vec::new();
    let mut kinds = Vec::new();
    let mut nodes = Vec::new();
    let mut seen_site_ids = BTreeSet::new();
    let mut seen_kinds = BTreeSet::new();
    let mut seen_nodes = BTreeSet::new();

    for (_, chunk) in chunks {
        extend_unique(&mut site_ids, &mut seen_site_ids, &chunk.site_ids);
        extend_unique(&mut kinds, &mut seen_kinds, &chunk.kinds);
        extend_unique(&mut nodes, &mut seen_nodes, &chunk.nodes);
    }

    DirtyChunk {
        id: format!("{}..{}", first.id, last.id),
        path_name: first.path_name.clone(),
        path_start_bp: option_min(chunks.iter().map(|(_, chunk)| chunk.path_start_bp)),
        path_end_bp: option_max(chunks.iter().map(|(_, chunk)| chunk.path_end_bp)),
        core_path_start_bp: option_min(chunks.iter().map(|(_, chunk)| chunk.core_path_start_bp)),
        core_path_end_bp: option_max(chunks.iter().map(|(_, chunk)| chunk.core_path_end_bp)),
        path_length_bp: first
            .path_length_bp
            .or_else(|| chunks.iter().find_map(|(_, chunk)| chunk.path_length_bp)),
        graph_start_bp: option_min(chunks.iter().map(|(_, chunk)| chunk.graph_start_bp)),
        graph_end_bp: option_max(chunks.iter().map(|(_, chunk)| chunk.graph_end_bp)),
        start_order: option_min(chunks.iter().map(|(_, chunk)| chunk.start_order)),
        end_order: option_max(chunks.iter().map(|(_, chunk)| chunk.end_order)),
        site_count: chunks
            .iter()
            .map(|(_, chunk)| chunk.site_count)
            .sum::<usize>(),
        site_ids,
        kinds,
        nodes,
        max_severity: chunks
            .iter()
            .map(|(_, chunk)| chunk.max_severity)
            .fold(first.max_severity, f64::max),
        flank_bp: chunks
            .iter()
            .map(|(_, chunk)| chunk.flank_bp)
            .max()
            .unwrap_or(first.flank_bp),
        capped_by_budget: chunks.iter().any(|(_, chunk)| chunk.capped_by_budget),
    }
}

fn extend_unique(out: &mut Vec<String>, seen: &mut BTreeSet<String>, values: &[String]) {
    for value in values {
        if seen.insert(value.clone()) {
            out.push(value.clone());
        }
    }
}

fn option_min(values: impl Iterator<Item = Option<usize>>) -> Option<usize> {
    values.flatten().min()
}

fn option_max(values: impl Iterator<Item = Option<usize>>) -> Option<usize> {
    values.flatten().max()
}

fn chunk_path_span(chunk: &DirtyChunk) -> Option<(usize, usize)> {
    chunk
        .path_start_bp
        .zip(chunk.path_end_bp)
        .or_else(|| chunk.core_path_start_bp.zip(chunk.core_path_end_bp))
}

fn chunk_core_span(chunk: &DirtyChunk) -> Option<(usize, usize)> {
    chunk
        .core_path_start_bp
        .zip(chunk.core_path_end_bp)
        .or_else(|| chunk.path_start_bp.zip(chunk.path_end_bp))
}

impl DirtyChunkProvenance {
    fn covers(&self, chunk: &DirtyChunk) -> bool {
        self.covers_nodes(chunk) || self.covers_path_interval(chunk)
    }

    fn remember(&mut self, chunk: &DirtyChunk) {
        for node in &chunk.nodes {
            self.nodes.insert(node.clone());
        }
        if let (Some(path_name), Some((start, end))) =
            (chunk.path_name.as_ref(), chunk_core_span(chunk))
        {
            if end > start {
                self.path_intervals
                    .entry(path_name.clone())
                    .or_default()
                    .push((start, end));
            }
        }
    }

    fn covers_nodes(&self, chunk: &DirtyChunk) -> bool {
        if chunk.nodes.is_empty() {
            return false;
        }
        let unique_nodes = chunk.nodes.iter().collect::<BTreeSet<_>>();
        let overlap = unique_nodes
            .iter()
            .filter(|node| self.nodes.contains(node.as_str()))
            .count();
        overlap == unique_nodes.len() || (overlap >= 2 && overlap * 2 >= unique_nodes.len())
    }

    fn covers_path_interval(&self, chunk: &DirtyChunk) -> bool {
        let (Some(path_name), Some((start, end))) =
            (chunk.path_name.as_ref(), chunk_core_span(chunk))
        else {
            return false;
        };
        if end <= start {
            return false;
        }
        self.path_intervals.get(path_name).is_some_and(|intervals| {
            intervals.iter().any(|&(seen_start, seen_end)| {
                interval_covers_half(start, end, seen_start, seen_end)
            })
        })
    }
}

fn interval_covers_half(start: usize, end: usize, seen_start: usize, seen_end: usize) -> bool {
    let overlap_start = start.max(seen_start);
    let overlap_end = end.min(seen_end);
    if overlap_end <= overlap_start {
        return false;
    }
    overlap_end.saturating_sub(overlap_start) * 2 >= end.saturating_sub(start)
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
    let path = Path::new(debug_dir).join("localized_polish_alignment_yield.tsv");
    if let Err(err) = std::fs::write(&path, alignment_yield_tsv(report)) {
        log::warn!(
            "localized polish: failed to write alignment-yield summary {}: {err}",
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

fn alignment_yield_tsv(report: &LocalizedPolishReport) -> String {
    let mut out = String::new();
    out.push_str(
        "record_type\titeration\tchunk_id\tstatus\treason\tmethod\tpath_count\tcandidate_sequence_count\tcandidate_total_bp\tstep_range\tpath_range_bp\tpair_count_expected\traw_paf_records\tfiltered_paf_records\taligned_pair_fraction\tmedian_alignment_len\tmax_alignment_len\tmin_match_len\tmin_map_length\tnum_mappings\tscaffold_filter\tscaffold_mass\tno_filter\tseqwish_input_paf\treplacement_segments\treplacement_bp\treplacement_shared_segments\tseqwish_segments\tseqwish_segment_bp\tseqwish_links\tseqwish_paths\n",
    );
    for iteration in &report.iterations {
        for resolver in &iteration.resolver_reports {
            let evidence = resolver.alignment.as_ref();
            let _ = writeln!(
                out,
                "alignment\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                iteration.iteration,
                tsv_field(&resolver.chunk_id),
                resolver.status.as_str(),
                tsv_field(&resolver.reason),
                resolver
                    .method
                    .map(|method| method.method_name())
                    .unwrap_or("n/a"),
                resolver.input_paths,
                evidence.map(|e| e.candidate_sequence_count).unwrap_or(resolver.input_paths),
                evidence.map(|e| e.candidate_total_bp).unwrap_or(resolver.input_bp),
                evidence
                    .and_then(|e| e.step_range.as_ref())
                    .map(|s| tsv_field(s))
                    .unwrap_or_default(),
                evidence
                    .and_then(|e| e.path_range_bp.as_ref())
                    .map(|s| tsv_field(s))
                    .unwrap_or_default(),
                evidence.map(|e| e.pair_count_expected).unwrap_or(0),
                evidence.map(|e| e.raw_paf_records).unwrap_or(0),
                evidence.map(|e| e.filtered_paf_records).unwrap_or(0),
                evidence
                    .and_then(|e| e.aligned_pair_fraction)
                    .map(|v| format!("{v:.6}"))
                    .unwrap_or_default(),
                evidence
                    .and_then(|e| e.median_alignment_len)
                    .map(|v| v.to_string())
                    .unwrap_or_default(),
                evidence.map(|e| e.max_alignment_len).unwrap_or(0),
                evidence.map(|e| e.min_match_len).unwrap_or(0),
                evidence.map(|e| e.min_map_length).unwrap_or(0),
                evidence
                    .map(|e| tsv_field(&e.num_mappings))
                    .unwrap_or_default(),
                evidence
                    .map(|e| tsv_field(&e.scaffold_filter))
                    .unwrap_or_default(),
                evidence.map(|e| e.scaffold_mass).unwrap_or(0),
                evidence.map(|e| e.no_filter).unwrap_or(false),
                evidence
                    .and_then(|e| e.seqwish_input_paf.as_ref())
                    .map(|s| tsv_field(s))
                    .unwrap_or_default(),
                resolver.replacement_segments.unwrap_or(0),
                resolver.replacement_bp.unwrap_or(0),
                evidence.map(|e| e.replacement_shared_segments).unwrap_or(0),
                evidence.map(|e| e.seqwish_segments).unwrap_or(0),
                evidence.map(|e| e.seqwish_segment_bp).unwrap_or(0),
                evidence.map(|e| e.seqwish_links).unwrap_or(0),
                evidence.map(|e| e.seqwish_paths).unwrap_or(0),
            );
        }
    }
    out
}

fn tsv_field(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_dirty_chunk(
        id: &str,
        path_start: usize,
        path_end: usize,
        core_start: usize,
        core_end: usize,
        nodes: &[&str],
    ) -> DirtyChunk {
        test_dirty_chunk_on_path(id, "ref", path_start, path_end, core_start, core_end, nodes)
    }

    fn test_dirty_chunk_on_path(
        id: &str,
        path_name: &str,
        path_start: usize,
        path_end: usize,
        core_start: usize,
        core_end: usize,
        nodes: &[&str],
    ) -> DirtyChunk {
        DirtyChunk {
            id: id.to_string(),
            path_name: Some(path_name.to_string()),
            path_start_bp: Some(path_start),
            path_end_bp: Some(path_end),
            core_path_start_bp: Some(core_start),
            core_path_end_bp: Some(core_end),
            path_length_bp: Some(1_000),
            graph_start_bp: Some(path_start),
            graph_end_bp: Some(path_end),
            start_order: Some(path_start),
            end_order: Some(path_end),
            site_count: 1,
            site_ids: vec![format!("{id}_site")],
            kinds: vec!["synthetic-scar".to_string()],
            nodes: nodes.iter().map(|node| (*node).to_string()).collect(),
            max_severity: 1.0,
            flank_bp: 4,
            capped_by_budget: false,
        }
    }

    #[test]
    fn localized_polish_compounds_overlapping_dirty_windows_before_resolution() {
        let mut config = LocalizedPolishConfig::default();
        config.max_chunks_per_iteration = 8;
        config.max_total_chunks = 8;
        config.max_total_chunk_bp = Some(1_000);
        config.dirty_options.flank_bp = 4;

        let chunks = vec![
            test_dirty_chunk("chunk_0001", 10, 20, 14, 16, &["n2"]),
            test_dirty_chunk("chunk_0002", 18, 28, 22, 24, &["n3"]),
        ];

        let selection = select_chunks(&chunks, &config, 0, 0);
        assert_eq!(
            selection.chunks.len(),
            1,
            "overlapping flanked windows from one scar should be resolved as a single anchored compound window"
        );
        let compound = &selection.chunks[0];
        assert_eq!(compound.path_start_bp, Some(10));
        assert_eq!(compound.path_end_bp, Some(28));
        assert_eq!(compound.core_path_start_bp, Some(14));
        assert_eq!(compound.core_path_end_bp, Some(24));
        assert_eq!(compound.site_count, 2);
        assert_eq!(selection.selected_bp, 18);
    }

    #[test]
    fn localized_polish_defers_duplicate_dirty_interiors_by_provenance() {
        let mut config = LocalizedPolishConfig::default();
        config.max_chunks_per_iteration = 8;
        config.max_total_chunks = 8;
        config.max_total_chunk_bp = Some(1_000);

        let chunks = vec![
            test_dirty_chunk_on_path("chunk_0001", "GRCh38", 10, 20, 14, 16, &["n2", "n3"]),
            test_dirty_chunk_on_path("chunk_0002", "HG00097", 40, 50, 44, 46, &["n2", "n3"]),
        ];

        let selection = select_chunks(&chunks, &config, 0, 0);
        assert_eq!(
            selection.chunks.len(),
            1,
            "duplicate interiors on another path should not be reprocessed in the same selection"
        );
        assert_eq!(selection.chunks[0].id, "chunk_0001");
        assert_eq!(selection.deferred, 1);
    }

    #[test]
    fn localized_polish_prefers_anchorable_interior_chunks() {
        let mut config = LocalizedPolishConfig::default();
        config.max_chunks_per_iteration = 1;
        config.max_total_chunks = 1;
        config.max_total_chunk_bp = Some(1_000);
        config.dirty_options.merge_distance_bp = 0;
        config.dirty_options.flank_bp = 0;

        let chunks = vec![
            test_dirty_chunk("left_edge", 0, 10, 0, 10, &["n1"]),
            test_dirty_chunk("right_edge", 990, 1_000, 990, 1_000, &["n9"]),
            test_dirty_chunk("interior", 100, 140, 110, 130, &["n5"]),
        ];

        let selection = select_chunks(&chunks, &config, 0, 0);
        assert_eq!(
            selection.chunks.len(),
            1,
            "budget should be spent on a chunk that can have both flank anchors"
        );
        assert_eq!(selection.chunks[0].id, "interior");
        assert_eq!(selection.deferred, 2);
    }

    #[test]
    fn localized_polish_defaults_do_not_hide_one_to_one_pairwise_filter() {
        let config = LocalizedPolishConfig::default();
        assert_eq!(
            config.resolver.resolution.replacement_num_mappings,
            "many:many"
        );
        assert_eq!(
            config.resolver.resolution.replacement_scaffold_filter,
            "many:many"
        );
        assert_eq!(config.resolver.resolution.replacement_scaffold_mass, 0);
    }
}
