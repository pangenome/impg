use crate::impg::{AdjustedInterval, Impg};
use crate::subset_filter::SubsetFilter;
use log::{debug, info, warn};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::io;

/// Configuration parameters for the refinement routine
pub struct RefineConfig<'a> {
    pub span_bp: i32,
    pub max_extension: i32,
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

/// Summary for each refined interval
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
}

struct SampleInterval {
    query_start: i32,
    query_end: i32,
    target_start: i32,
    target_end: i32,
}

struct CandidateResult {
    start: i32,
    end: i32,
    left_extension: i32,
    right_extension: i32,
    support_count: usize,
}

/// Run the refinement procedure on the provided ranges
pub fn run_refine(
    impg: &Impg,
    ranges: &[(String, (i32, i32), String)],
    config: RefineConfig<'_>,
) -> io::Result<Vec<RefineRecord>> {
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

    let flanks = build_flanks(config.max_extension, config.extension_step);
    debug!(
        "Evaluating {}x{} flank sizes for region {}:{}-{}",
        flanks.len(),
        flanks.len(),
        chrom,
        orig_start,
        orig_end
    );

    let best_candidate = flanks
        .par_iter()
        .flat_map_iter(|&left_flank| {
            flanks
                .iter()
                .map(move |&right_flank| (left_flank, right_flank))
        })
        .filter_map(|(left_flank, right_flank)| {
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

            let support_count = compute_supporting_samples(
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
                support_count,
            };

            debug!(
                "Region {}:{}-{} flank L{} R{} -> {} supporting samples",
                chrom, start, end, left_extension, right_extension, candidate.support_count
            );

            Some(candidate)
        })
        .reduce_with(|best, candidate| {
            if compare_candidates(&candidate, &best) == Ordering::Greater {
                candidate
            } else {
                best
            }
        });

    let Some(best_candidate) = best_candidate else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "No valid flank sizes evaluated for region {}:{}-{}",
                chrom, orig_start, orig_end
            ),
        ));
    };

    info!(
        "Selected flanks left={}bp right={}bp for region {}:{}-{} ({} supporting samples)",
        best_candidate.left_extension,
        best_candidate.right_extension,
        chrom,
        best_candidate.start,
        best_candidate.end,
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
    })
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
        )
    } else {
        impg.query(target_id, range.0, range.1, false, config.min_identity)
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

fn compute_supporting_samples(
    target_id: u32,
    overlaps: &[AdjustedInterval],
    region_start: i32,
    region_end: i32,
    span_bp: i32,
    merge_distance: i32,
) -> usize {
    if overlaps.len() <= 1 {
        return 0;
    }

    let mut per_sample: FxHashMap<u32, Vec<SampleInterval>> = FxHashMap::default();

    for (query_interval, _, target_interval) in overlaps.iter() {
        if query_interval.metadata == target_id {
            continue; // Skip the target region itself
        }

        // Store absolute coordinates regardless of strand
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

    let mut support_ids = Vec::new();

    for (sample_id, intervals) in per_sample {
        let merged = merge_intervals(intervals, merge_distance);

        if merged.iter().any(|interval| {
            covers_boundaries(
                interval.target_start,
                interval.target_end,
                region_start,
                region_end,
                left_threshold,
                right_threshold,
            )
        }) {
            support_ids.push(sample_id);
        }
    }

    support_ids.sort_unstable();
    support_ids.dedup();
    support_ids.len()
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
