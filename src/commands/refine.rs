use crate::impg::{AdjustedInterval, Impg};
use crate::subset_filter::SubsetFilter;
use clap::ValueEnum;
use log::{debug, info, warn};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Ordering;
use std::io;

/// Configuration parameters for the refinement routine.
/// Mirrors CLI flags and constrains how aggressively flanks can be explored while
/// searching for loci that remain well supported at both boundaries.
pub struct RefineConfig<'a> {
    pub span_bp: i32,
    /// Maximum per-side expansion; <=1 interpreted as fraction of the locus, >1 as absolute bp.
    pub max_extension: f64,
    /// Aggregation mode used when counting boundary support.
    pub support_mode: SupportMode,
    pub extension_step: i32,
    pub merge_distance: i32,
    pub min_identity: Option<f64>,
    pub use_transitive_bfs: bool,
    pub use_transitive_dfs: bool,
    pub max_transitive_depth: u16,
    pub min_transitive_len: i32,
    pub min_distance_between_ranges: i32,
    pub subset_filter: Option<&'a SubsetFilter>,
}

/// Summary for each refined interval produced by [`run_refine`].
pub struct RefineRecord {
    pub chrom: String,
    pub refined_start: i32,
    pub refined_end: i32,
    pub original_start: i32,
    pub original_end: i32,
    pub label: String,
    pub applied_left_extension: i32,
    pub applied_right_extension: i32,
    pub support_count: usize,
    pub original_support_count: usize,
    pub support_entities: Vec<SupportEntity>,
}

#[derive(Clone, Debug)]
pub struct SupportEntity {
    pub sequence: String,
    pub start: i32,
    pub end: i32,
}

/// How to aggregate PanSN identifiers when counting support.
#[derive(Clone, Copy, Debug)]
pub enum SupportMode {
    Sequence,
    Sample,
    Haplotype,
}

/// CLI-level PanSN aggregation mode (sample/haplotype).
#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum RefineSupportArg {
    Sample,
    Haplotype,
}

fn support_mode_label(mode: SupportMode) -> &'static str {
    match mode {
        SupportMode::Sequence => "sequences",
        SupportMode::Sample => "samples",
        SupportMode::Haplotype => "haplotypes",
    }
}

impl From<RefineSupportArg> for SupportMode {
    fn from(value: RefineSupportArg) -> Self {
        match value {
            RefineSupportArg::Sample => SupportMode::Sample,
            RefineSupportArg::Haplotype => SupportMode::Haplotype,
        }
    }
}

struct SampleInterval {
    query_start: i32,
    query_end: i32,
    target_start: i32,
    target_end: i32,
}

/// Candidate solution capturing a concrete left/right expansion and its support.
#[derive(Clone)]
struct CandidateResult {
    start: i32,
    end: i32,
    left_extension: i32,
    right_extension: i32,
    support_count: usize,
    support_entities: Vec<SupportEntity>,
}

/// Run the refinement procedure on the provided ranges
pub fn run_refine(
    impg: &Impg,
    ranges: &[(String, (i32, i32), String)],
    config: RefineConfig<'_>,
) -> io::Result<Vec<RefineRecord>> {
    info!(
        "Refining {} range(s) with support aggregated by {} (span_bp={}, max_extension={}, extension_step={})",
        ranges.len(),
        support_mode_label(config.support_mode),
        config.span_bp,
        config.max_extension,
        config.extension_step
    );

    let intermediate: Vec<Result<(usize, RefineRecord), io::Error>> = ranges
        .par_iter()
        .enumerate()
        .map(|(idx, (chrom, (orig_start, orig_end), label))| {
            refine_single_range(impg, chrom, *orig_start, *orig_end, label, &config)
                .map(|record| (idx, record))
        })
        .collect();

    let mut ordered = Vec::with_capacity(intermediate.len());
    for entry in intermediate {
        let (idx, record) = entry?;
        ordered.push((idx, record));
    }

    ordered.sort_unstable_by_key(|(idx, _)| *idx);

    Ok(ordered.into_iter().map(|(_, record)| record).collect())
}

fn refine_single_range(
    impg: &Impg,
    chrom: &str,
    orig_start: i32,
    orig_end: i32,
    label: &str,
    config: &RefineConfig<'_>,
) -> io::Result<RefineRecord> {
    if orig_end <= orig_start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Invalid range {}:{}-{} (end must be greater than start)",
                chrom, orig_start, orig_end
            ),
        ));
    }

    let Some(target_id) = impg.seq_index.get_id(chrom) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Target sequence '{chrom}' not found in index"),
        ));
    };
    let Some(seq_len) = impg.seq_index.get_len_from_id(target_id) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Length for sequence '{chrom}' missing from index"),
        ));
    };
    let seq_len = seq_len as i32;

    let locus_len = (orig_end - orig_start).max(0);
    let max_extension_bp = if config.max_extension <= 1.0 {
        ((locus_len as f64) * config.max_extension)
            .ceil()
            .clamp(0.0, i32::MAX as f64) as i32
    } else {
        config.max_extension.ceil().clamp(0.0, i32::MAX as f64) as i32
    };
    let max_extension_bp = max_extension_bp.max(0);

    // Build the grid of candidate flank sizes based on dynamic constraints.
    let flanks = build_flanks(max_extension_bp, config.extension_step);
    debug!(
        "Evaluating {} flank steps up to {} bp per side (max_extension={}) for region {}:{}-{}",
        flanks.len(),
        max_extension_bp,
        if config.max_extension <= 1.0 {
            format!("{:.2}x locus", config.max_extension)
        } else {
            format!("{:.0} bp", config.max_extension)
        },
        chrom,
        orig_start,
        orig_end
    );

    let evaluate = |left: i32, right: i32| -> Option<CandidateResult> {
        evaluate_candidate(
            impg, target_id, chrom, orig_start, orig_end, seq_len, left, right, config,
        )
    };

    let mut best_candidate: Option<CandidateResult> = None;
    let baseline_candidate = evaluate(0, 0);
    let original_support_count = baseline_candidate
        .as_ref()
        .map(|candidate| candidate.support_count)
        .unwrap_or(0);

    if let Some(candidate) = baseline_candidate {
        best_candidate = update_best_candidate(best_candidate, candidate);
    }

    // First, search for the minimal left extension that still satisfies the span requirement on the left boundary.
    if let Some(candidate) = flanks
        .par_iter()
        .filter_map(|&left| evaluate(left, 0))
        .reduce_with(better_candidate)
    {
        best_candidate = update_best_candidate(best_candidate, candidate);
    }

    let left_fixed = best_candidate
        .as_ref()
        .map(|c| c.left_extension)
        .unwrap_or(0);

    // Next, keep the chosen left flank fixed and look for the best right expansion.
    if let Some(candidate) = flanks
        .par_iter()
        .filter_map(|&right| evaluate(left_fixed, right))
        .reduce_with(better_candidate)
    {
        best_candidate = update_best_candidate(best_candidate, candidate);
    }

    let right_fixed = best_candidate
        .as_ref()
        .map(|c| c.right_extension)
        .unwrap_or(0);

    // Finally, re-optimise the left flank while holding the right flank steady to capture asymmetric trade-offs.
    if let Some(candidate) = flanks
        .par_iter()
        .filter_map(|&left| evaluate(left, right_fixed))
        .reduce_with(better_candidate)
    {
        best_candidate = update_best_candidate(best_candidate, candidate);
    }

    let Some(best_candidate) = best_candidate else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "No valid flank sizes evaluated for region {}:{}-{}",
                chrom, orig_start, orig_end
            ),
        ));
    };

    debug!(
        "Selected flanks left={}bp right={}bp for region {}:{}-{} (supporting {}: {} -> {})",
        best_candidate.left_extension,
        best_candidate.right_extension,
        chrom,
        best_candidate.start,
        best_candidate.end,
        support_mode_label(config.support_mode),
        original_support_count,
        best_candidate.support_count
    );

    Ok(RefineRecord {
        chrom: chrom.to_string(),
        refined_start: best_candidate.start,
        refined_end: best_candidate.end,
        original_start: orig_start,
        original_end: orig_end,
        label: label.to_string(),
        applied_left_extension: best_candidate.left_extension,
        applied_right_extension: best_candidate.right_extension,
        support_count: best_candidate.support_count,
        original_support_count,
        support_entities: best_candidate.support_entities,
    })
}

/// Evaluate a single `(left_flank, right_flank)` combination and return its support summary.
fn evaluate_candidate(
    impg: &Impg,
    target_id: u32,
    chrom: &str,
    orig_start: i32,
    orig_end: i32,
    seq_len: i32,
    left_flank: i32,
    right_flank: i32,
    config: &RefineConfig<'_>,
) -> Option<CandidateResult> {
    let tentative_start = orig_start.saturating_sub(left_flank);
    let tentative_end = orig_end.saturating_add(right_flank);

    let start = tentative_start.max(0);
    let end = tentative_end.min(seq_len);

    if end <= start {
        warn!(
            "Skipping non-positive range {}:{}-{} after applying flanks L{} R{}",
            chrom, start, end, left_flank, right_flank
        );
        return None;
    }

    let mut overlaps = Vec::new();
    query_overlaps(impg, target_id, (start, end), config, &mut overlaps);

    apply_subset_filter(
        impg,
        target_id,
        &mut overlaps,
        chrom,
        (start, end),
        config.subset_filter,
    );

    let stats = compute_supporting_stats(
        impg,
        config.support_mode,
        target_id,
        &overlaps,
        start,
        end,
        config.span_bp,
        config.merge_distance,
    );

    let left_extension = orig_start.saturating_sub(start);
    let right_extension = end.saturating_sub(orig_end);

    let candidate = CandidateResult {
        start,
        end,
        left_extension,
        right_extension,
        support_count: stats.aggregated,
        support_entities: stats.survivors,
    };

    Some(candidate)
}

fn query_overlaps(
    impg: &Impg,
    target_id: u32,
    range: (i32, i32),
    config: &RefineConfig<'_>,
    buffer: &mut Vec<AdjustedInterval>,
) {
    buffer.clear();

    let mut results = if config.use_transitive_bfs {
        impg.query_transitive_bfs(
            target_id,
            range.0,
            range.1,
            None,
            config.max_transitive_depth,
            config.min_transitive_len,
            config.min_distance_between_ranges,
            false,
            config.min_identity,
            None,
        )
    } else if config.use_transitive_dfs {
        impg.query_transitive_dfs(
            target_id,
            range.0,
            range.1,
            None,
            config.max_transitive_depth,
            config.min_transitive_len,
            config.min_distance_between_ranges,
            false,
            config.min_identity,
            None,
        )
    } else {
        impg.query(
            target_id,
            range.0,
            range.1,
            false,
            config.min_identity,
            None,
        )
    };

    buffer.append(&mut results);
}

fn apply_subset_filter(
    impg: &Impg,
    target_id: u32,
    overlaps: &mut Vec<AdjustedInterval>,
    chrom: &str,
    range: (i32, i32),
    subset_filter: Option<&SubsetFilter>,
) {
    if let Some(filter) = subset_filter {
        let before = overlaps.len();
        overlaps.retain(|(query_interval, _, _)| {
            if query_interval.metadata == target_id {
                return true;
            }
            impg.seq_index
                .get_name(query_interval.metadata)
                .map(|name| filter.matches(name))
                .unwrap_or(false)
        });

        let filtered_out = before.saturating_sub(overlaps.len());
        if filtered_out > 0 {
            debug!(
                "Filtered out {} sequences outside subset for region {}:{}-{}",
                filtered_out, chrom, range.0, range.1
            );
        }

        if overlaps.len() <= 1 {
            warn!(
                "Subset filtering left no comparison sequences for region {}:{}-{}",
                chrom, range.0, range.1
            );
        }
    }
}

/// Return whichever candidate ranks higher according to `compare_candidates`.
fn better_candidate(a: CandidateResult, b: CandidateResult) -> CandidateResult {
    if compare_candidates(&a, &b) == Ordering::Greater {
        a
    } else {
        b
    }
}

/// Maintain the best-so-far candidate when combining sequential passes.
fn update_best_candidate(
    best: Option<CandidateResult>,
    candidate: CandidateResult,
) -> Option<CandidateResult> {
    match best {
        Some(current) => {
            if compare_candidates(&candidate, &current) == Ordering::Greater {
                Some(candidate)
            } else {
                Some(current)
            }
        }
        None => Some(candidate),
    }
}

fn compare_candidates(a: &CandidateResult, b: &CandidateResult) -> Ordering {
    a.support_count
        .cmp(&b.support_count)
        .then_with(|| {
            let a_total = a.left_extension + a.right_extension;
            let b_total = b.left_extension + b.right_extension;
            b_total.cmp(&a_total)
        })
        .then_with(|| {
            let a_max = a.left_extension.max(a.right_extension);
            let b_max = b.left_extension.max(b.right_extension);
            b_max.cmp(&a_max)
        })
        .then_with(|| {
            let a_len = a.end - a.start;
            let b_len = b.end - b.start;
            b_len.cmp(&a_len)
        })
}

struct SupportStats {
    aggregated: usize,
    survivors: Vec<SupportEntity>,
}

fn compute_supporting_stats(
    impg: &Impg,
    mode: SupportMode,
    target_id: u32,
    overlaps: &[AdjustedInterval],
    region_start: i32,
    region_end: i32,
    span_bp: i32,
    merge_distance: i32,
) -> SupportStats {
    let (aggregated, survivors) = compute_support_sets(
        impg,
        mode,
        target_id,
        overlaps,
        region_start,
        region_end,
        span_bp,
        merge_distance,
    );

    SupportStats {
        aggregated: aggregated.len(),
        survivors,
    }
}

fn compute_support_sets(
    impg: &Impg,
    mode: SupportMode,
    target_id: u32,
    overlaps: &[AdjustedInterval],
    region_start: i32,
    region_end: i32,
    span_bp: i32,
    merge_distance: i32,
) -> (FxHashSet<String>, Vec<SupportEntity>) {
    let mut aggregated = FxHashSet::default();
    let mut sequence_ranges: FxHashMap<String, (i32, i32)> = FxHashMap::default();

    if overlaps.len() <= 1 {
        return (aggregated, Vec::new());
    }

    let mut per_sample: FxHashMap<u32, Vec<SampleInterval>> = FxHashMap::default();

    for (query_interval, _, target_interval) in overlaps.iter() {
        if query_interval.metadata == target_id {
            continue;
        }

        let q_start = query_interval.first.min(query_interval.last);
        let q_end = query_interval.first.max(query_interval.last);
        let t_start = target_interval.first.min(target_interval.last);
        let t_end = target_interval.first.max(target_interval.last);

        per_sample
            .entry(query_interval.metadata)
            .or_default()
            .push(SampleInterval {
                query_start: q_start,
                query_end: q_end,
                target_start: t_start,
                target_end: t_end,
            });
    }

    let effective_span = (region_end - region_start).max(0).min(span_bp.max(0));
    let left_threshold = region_start + effective_span;
    let right_threshold = region_end - effective_span;

    for (sample_id, intervals) in per_sample {
        let merged = merge_intervals(intervals, merge_distance);

        let mut query_range: Option<(i32, i32)> = None;
        for interval in &merged {
            if covers_boundaries(
                interval.target_start,
                interval.target_end,
                region_start,
                region_end,
                left_threshold,
                right_threshold,
            ) {
                let q_start = interval.query_start.min(interval.query_end);
                let q_end = interval.query_start.max(interval.query_end);
                query_range = Some(match query_range {
                    Some((s, e)) => (s.min(q_start), e.max(q_end)),
                    None => (q_start, q_end),
                });
            }
        }

        if let (Some((q_start, q_end)), Some(name)) =
            (query_range, impg.seq_index.get_name(sample_id))
        {
            let seq_name = name.to_string();
            let entry = sequence_ranges
                .entry(seq_name.clone())
                .or_insert((q_start, q_end));
            entry.0 = entry.0.min(q_start);
            entry.1 = entry.1.max(q_end);
            if let Some(key) = pansn_key(name, mode) {
                aggregated.insert(key);
            }
        }
    }

    let mut survivors: Vec<SupportEntity> = sequence_ranges
        .into_iter()
        .map(|(sequence, (start, end))| SupportEntity {
            sequence,
            start,
            end,
        })
        .collect();
    survivors.sort_unstable_by(|a, b| {
        a.sequence
            .cmp(&b.sequence)
            .then_with(|| a.start.cmp(&b.start))
    });

    (aggregated, survivors)
}

fn covers_boundaries(
    interval_start: i32,
    interval_end: i32,
    region_start: i32,
    region_end: i32,
    left_threshold: i32,
    right_threshold: i32,
) -> bool {
    interval_start <= region_start
        && interval_end >= region_end
        && interval_end >= left_threshold
        && interval_start <= right_threshold
}

fn merge_intervals(mut intervals: Vec<SampleInterval>, merge_distance: i32) -> Vec<SampleInterval> {
    if intervals.is_empty() {
        return intervals;
    }

    if merge_distance < 0 {
        return intervals;
    }

    intervals.sort_by(|a, b| {
        a.query_start
            .cmp(&b.query_start)
            .then_with(|| a.query_end.cmp(&b.query_end))
    });

    let mut merged = Vec::with_capacity(intervals.len());
    let mut iter = intervals.into_iter();
    let mut current = iter.next().unwrap();

    for next in iter {
        if should_merge(&current, &next, merge_distance) {
            current.query_start = current.query_start.min(next.query_start);
            current.query_end = current.query_end.max(next.query_end);
            current.target_start = current.target_start.min(next.target_start);
            current.target_end = current.target_end.max(next.target_end);
        } else {
            merged.push(current);
            current = next;
        }
    }

    merged.push(current);
    merged
}

fn should_merge(a: &SampleInterval, b: &SampleInterval, merge_distance: i32) -> bool {
    if merge_distance < 0 {
        return false;
    }
    let distance = merge_distance as u32;
    let query_adjacent = a
        .query_end
        .abs_diff(b.query_start)
        .min(a.query_start.abs_diff(b.query_end))
        <= distance;
    let target_adjacent = a
        .target_end
        .abs_diff(b.target_start)
        .min(a.target_start.abs_diff(b.target_end))
        <= distance;
    query_adjacent || target_adjacent
}

fn build_flanks(max_extension: i32, step: i32) -> Vec<i32> {
    let mut flanks = Vec::new();
    let mut current = 0;

    if max_extension == 0 {
        flanks.push(0);
        return flanks;
    }

    while current <= max_extension {
        flanks.push(current);
        if max_extension - current < step {
            break;
        }
        current = current.saturating_add(step);
    }

    if flanks.last().copied() != Some(max_extension) {
        flanks.push(max_extension);
    }

    flanks.sort_unstable();
    flanks.dedup();
    flanks
}
fn pansn_key(name: &str, mode: SupportMode) -> Option<String> {
    match mode {
        SupportMode::Sequence => Some(name.to_string()),
        SupportMode::Sample => {
            let base = name.split(':').next().unwrap_or(name);
            let sample = base.split('#').next().unwrap_or(base).trim();
            if sample.is_empty() {
                None
            } else {
                Some(sample.to_string())
            }
        }
        SupportMode::Haplotype => {
            let base = name.split(':').next().unwrap_or(name);
            let mut parts = base.split('#');
            let sample = parts.next().unwrap_or("").trim();
            if sample.is_empty() {
                return None;
            }
            match parts.next() {
                Some(hap) if !hap.is_empty() => Some(format!("{sample}#{hap}")),
                _ => Some(sample.to_string()),
            }
        }
    }
}
