use crate::impg::{AdjustedInterval, Impg};
use crate::sequence_index::UnifiedSequenceIndex;
use log::{debug, info};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::{self, BufWriter, Write};

/// Configuration for the depth command
pub struct DepthConfig {
    pub transitive: bool,
    pub transitive_dfs: bool,
    pub max_depth: u16,
    pub min_transitive_len: i32,
    pub min_distance_between_ranges: i32,
    pub merge_adjacent: bool,
    pub approximate_mode: bool,
}

/// Tracks processed regions for each haplotype to avoid duplicate output
#[derive(Default)]
pub struct ProcessedRegions {
    /// haplotype_name -> sorted list of (start, end) intervals
    regions: FxHashMap<String, Vec<(i64, i64)>>,
}

impl ProcessedRegions {
    pub fn new() -> Self {
        Self {
            regions: FxHashMap::default(),
        }
    }

    /// Check if a region overlaps with any processed region for this haplotype
    /// Returns the unprocessed portion, or None if fully processed
    pub fn get_unprocessed(&self, haplotype: &str, start: i64, end: i64) -> Option<(i64, i64)> {
        if let Some(intervals) = self.regions.get(haplotype) {
            let mut current_start = start;

            for &(proc_start, proc_end) in intervals {
                if proc_end <= current_start {
                    continue; // Processed region is before our region
                }
                if proc_start >= end {
                    break; // Processed region is after our region
                }

                // There's overlap
                if proc_start <= current_start {
                    // Our start is inside a processed region
                    current_start = proc_end;
                    if current_start >= end {
                        return None; // Fully processed
                    }
                } else {
                    // There's an unprocessed gap before this processed region
                    return Some((current_start, proc_start.min(end)));
                }
            }

            if current_start < end {
                Some((current_start, end))
            } else {
                None
            }
        } else {
            Some((start, end))
        }
    }

    /// Mark a region as processed for a haplotype
    pub fn mark_processed(&mut self, haplotype: &str, start: i64, end: i64) {
        let intervals = self.regions.entry(haplotype.to_string()).or_default();

        // Insert and merge overlapping intervals
        let mut new_start = start;
        let mut new_end = end;
        let mut merged_indices = Vec::new();

        for (i, &(s, e)) in intervals.iter().enumerate() {
            if e < new_start {
                continue; // Before our region
            }
            if s > new_end {
                break; // After our region
            }
            // Overlapping or adjacent - merge
            new_start = new_start.min(s);
            new_end = new_end.max(e);
            merged_indices.push(i);
        }

        // Remove merged intervals (in reverse order to preserve indices)
        for i in merged_indices.into_iter().rev() {
            intervals.remove(i);
        }

        // Insert the new merged interval at the correct position
        let pos = intervals
            .iter()
            .position(|&(s, _)| s > new_start)
            .unwrap_or(intervals.len());
        intervals.insert(pos, (new_start, new_end));
    }

    /// Check if a region is completely processed
    pub fn is_fully_processed(&self, haplotype: &str, start: i64, end: i64) -> bool {
        self.get_unprocessed(haplotype, start, end).is_none()
    }
}

/// Represents a depth window - a region with a specific coverage depth
#[derive(Debug, Clone)]
pub struct DepthWindow {
    pub window_id: usize,
    pub depth: usize,
    pub ref_name: String,
    pub ref_start: i64,
    pub ref_end: i64,
    /// haplotype_name -> (start, end) position on the reference
    pub haplotype_positions: FxHashMap<String, (i64, i64)>,
}

/// Event for sweep-line algorithm to compute depth
#[derive(Debug, Clone, PartialEq, Eq)]
struct DepthEvent {
    position: i64,
    is_start: bool,
    haplotype: String,
}

impl Ord for DepthEvent {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| other.is_start.cmp(&self.is_start)) // Starts before ends at same position
            .then_with(|| self.haplotype.cmp(&other.haplotype))
    }
}

impl PartialOrd for DepthEvent {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Extract sample from PanSN format sequence name (sample#haplotype#chr -> sample)
pub fn extract_sample(seq_name: &str, separator: &str) -> String {
    let parts: Vec<&str> = seq_name.split(separator).collect();
    if !parts.is_empty() {
        parts[0].to_string()
    } else {
        seq_name.to_string()
    }
}
/// Merge adjacent windows with the same depth and same haplotype set
fn merge_adjacent_windows(mut windows: Vec<DepthWindow>) -> Vec<DepthWindow> {
    if windows.len() <= 1 {
        return windows;
    }

    let mut merged: Vec<DepthWindow> = Vec::new();
    let mut current = windows.remove(0);

    for window in windows {
        // Check if windows are adjacent and have same haplotype set
        let same_haplotypes = current.haplotype_positions.len() == window.haplotype_positions.len()
            && current
                .haplotype_positions
                .keys()
                .all(|k| window.haplotype_positions.contains_key(k));

        if current.ref_end == window.ref_start
            && current.depth == window.depth
            && same_haplotypes
            && current.ref_name == window.ref_name
        {
            // Merge windows
            current.ref_end = window.ref_end;
            // Update haplotype positions to span both windows
            for (hap, (_, new_end)) in window.haplotype_positions {
                if let Some(pos) = current.haplotype_positions.get_mut(&hap) {
                    pos.1 = new_end;
                }
            }
        } else {
            merged.push(current);
            current = window;
        }
    }
    merged.push(current);

    // Renumber window IDs
    for (i, w) in merged.iter_mut().enumerate() {
        w.window_id = i;
    }

    merged
}

/// Query overlaps for a single reference (used in parallel phase)
fn query_overlaps_for_ref(
    impg: &Impg,
    ref_id: u32,
    ref_length: i64,
    config: &DepthConfig,
    sequence_index: Option<&UnifiedSequenceIndex>,
) -> Vec<AdjustedInterval> {
    if config.transitive {
        if config.transitive_dfs {
            impg.query_transitive_dfs(
                ref_id,
                0,
                ref_length as i32,
                None,
                config.max_depth,
                config.min_transitive_len,
                config.min_distance_between_ranges,
                None,
                false,
                None,
                sequence_index,
                config.approximate_mode,
            )
        } else {
            impg.query_transitive_bfs(
                ref_id,
                0,
                ref_length as i32,
                None,
                config.max_depth,
                config.min_transitive_len,
                config.min_distance_between_ranges,
                None,
                false,
                None,
                sequence_index,
                config.approximate_mode,
            )
        }
    } else {
        impg.query(
            ref_id,
            0,
            ref_length as i32,
            false,
            None,
            sequence_index,
            config.approximate_mode,
        )
    }
}

/// Process overlaps into depth windows (used in sequential phase with deduplication)
fn process_overlaps_to_windows(
    impg: &Impg,
    _ref_id: u32,
    ref_name: &str,
    ref_length: i64,
    overlaps: &[AdjustedInterval],
    config: &DepthConfig,
    processed_regions: &mut ProcessedRegions,
    separator: &str,
) -> Vec<DepthWindow> {
    let ref_sample = extract_sample(ref_name, separator);

    // Check if this reference region is already processed
    if processed_regions.is_fully_processed(&ref_sample, 0, ref_length) {
        debug!("Reference {} is fully processed, skipping", ref_name);
        return Vec::new();
    }

    if overlaps.is_empty() {
        debug!("No overlaps found for reference {}", ref_name);
        processed_regions.mark_processed(&ref_sample, 0, ref_length);
        return vec![DepthWindow {
            window_id: 0,
            depth: 1,
            ref_name: ref_name.to_string(),
            ref_start: 0,
            ref_end: ref_length,
            haplotype_positions: {
                let mut map = FxHashMap::default();
                map.insert(ref_sample.clone(), (0, ref_length));
                map
            },
        }];
    }

    // Group overlaps by sample, collecting ALL intervals (union approach)
    let mut sample_intervals: FxHashMap<String, Vec<(i64, i64)>> = FxHashMap::default();

    for overlap in overlaps {
        let query_interval = &overlap.0;
        let target_interval = &overlap.2;

        let query_name = impg.seq_index.get_name(query_interval.metadata).unwrap();
        let sample = extract_sample(query_name, separator);

        let target_start = target_interval.first.min(target_interval.last) as i64;
        let target_end = target_interval.first.max(target_interval.last) as i64;

        sample_intervals
            .entry(sample)
            .or_default()
            .push((target_start, target_end));
    }

    // Add self (reference) as a sample
    sample_intervals
        .entry(ref_sample.clone())
        .or_default()
        .push((0, ref_length));

    // Merge overlapping intervals for each sample and filter out processed regions
    let mut filtered_intervals: FxHashMap<String, Vec<(i64, i64)>> = FxHashMap::default();
    for (sample, mut intervals) in sample_intervals {
        // Sort and merge overlapping intervals
        intervals.sort_by_key(|(s, _)| *s);
        let mut merged: Vec<(i64, i64)> = Vec::new();
        for (start, end) in intervals {
            if let Some(last) = merged.last_mut() {
                if start <= last.1 {
                    // Overlapping or adjacent, extend
                    last.1 = last.1.max(end);
                } else {
                    merged.push((start, end));
                }
            } else {
                merged.push((start, end));
            }
        }

        // Filter out already processed regions
        let mut unprocessed: Vec<(i64, i64)> = Vec::new();
        for (start, end) in merged {
            // Get all unprocessed portions of this interval
            let mut current = start;
            while current < end {
                if let Some((unproc_start, unproc_end)) =
                    processed_regions.get_unprocessed(&sample, current, end)
                {
                    unprocessed.push((unproc_start, unproc_end));
                    current = unproc_end;
                } else {
                    break;
                }
            }
        }

        if !unprocessed.is_empty() {
            filtered_intervals.insert(sample, unprocessed);
        }
    }

    if filtered_intervals.is_empty() {
        debug!("All intervals for {} are already processed", ref_name);
        return Vec::new();
    }

    // Create events for sweep-line algorithm (one pair per interval)
    let mut events: Vec<DepthEvent> = Vec::new();
    for (sample, intervals) in &filtered_intervals {
        for (start, end) in intervals {
            events.push(DepthEvent {
                position: *start,
                is_start: true,
                haplotype: sample.clone(),
            });
            events.push(DepthEvent {
                position: *end,
                is_start: false,
                haplotype: sample.clone(),
            });
        }
    }
    events.sort();

    // Sweep-line to compute depth windows
    // Use counter to track multiple overlapping intervals from same sample
    let mut windows: Vec<DepthWindow> = Vec::new();
    let mut active_sample_counts: FxHashMap<String, usize> = FxHashMap::default();
    let mut prev_pos: Option<i64> = None;
    let mut window_id = 0;

    for event in events {
        if let Some(prev) = prev_pos {
            if event.position > prev {
                // Collect samples that are currently active
                let active_samples: Vec<&String> = active_sample_counts
                    .iter()
                    .filter(|(_, &count)| count > 0)
                    .map(|(s, _)| s)
                    .collect();

                if !active_samples.is_empty() {
                    let mut sample_positions: FxHashMap<String, (i64, i64)> = FxHashMap::default();
                    for sample in &active_samples {
                        // For this window, use the window boundaries as positions
                        sample_positions.insert((*sample).clone(), (prev, event.position));
                    }

                    windows.push(DepthWindow {
                        window_id,
                        depth: sample_positions.len(),
                        ref_name: ref_name.to_string(),
                        ref_start: prev,
                        ref_end: event.position,
                        haplotype_positions: sample_positions,
                    });
                    window_id += 1;
                }
            }
        }

        // Update active counts
        if event.is_start {
            *active_sample_counts.entry(event.haplotype).or_insert(0) += 1;
        } else {
            if let Some(count) = active_sample_counts.get_mut(&event.haplotype) {
                *count = count.saturating_sub(1);
            }
        }

        prev_pos = Some(event.position);
    }

    // Mark all processed regions
    for (sample, intervals) in &filtered_intervals {
        for (start, end) in intervals {
            processed_regions.mark_processed(sample, *start, *end);
        }
    }

    // Optionally merge adjacent windows
    if config.merge_adjacent {
        windows = merge_adjacent_windows(windows);
    }

    windows
}

/// Main function to compute depth across all sequences
pub fn compute_depth(
    impg: &Impg,
    config: &DepthConfig,
    ref_order: Option<Vec<String>>,
    sequence_index: Option<&UnifiedSequenceIndex>,
    separator: &str,
    output_prefix: Option<&str>,
) -> io::Result<()> {
    info!("Computing depth coverage");

    // Determine processing order
    let mut seq_order: Vec<(u32, String, i64)> = Vec::new();

    // First add user-specified sequences
    if let Some(order) = ref_order {
        for name in order {
            if let Some(id) = impg.seq_index.get_id(&name) {
                let len = impg.seq_index.get_len_from_id(id).unwrap_or(0) as i64;
                seq_order.push((id, name, len));
            } else {
                debug!("Sequence '{}' from ref-order not found in index", name);
            }
        }
    }

    // Add remaining sequences in alphabetical order
    let mut all_seqs: Vec<(u32, String, i64)> = (0..impg.seq_index.len() as u32)
        .filter_map(|id| {
            let name = impg.seq_index.get_name(id)?;
            let len = impg.seq_index.get_len_from_id(id)? as i64;
            Some((id, name.to_string(), len))
        })
        .collect();
    all_seqs.sort_by(|a, b| a.1.cmp(&b.1));

    let specified_ids: FxHashSet<u32> = seq_order.iter().map(|(id, _, _)| *id).collect();
    for (id, name, len) in all_seqs {
        if !specified_ids.contains(&id) {
            seq_order.push((id, name, len));
        }
    }

    info!("Processing {} sequences as references", seq_order.len());

    // Collect all unique samples for header
    let mut all_samples: FxHashSet<String> = FxHashSet::default();
    for (_, name, _) in &seq_order {
        all_samples.insert(extract_sample(name, separator));
    }
    let mut sample_list: Vec<String> = all_samples.into_iter().collect();
    sample_list.sort();

    // Phase 1: Parallel query all overlaps
    info!("Phase 1: Querying overlaps in parallel...");
    let all_overlaps: Vec<Vec<AdjustedInterval>> = seq_order
        .par_iter()
        .map(|(ref_id, _ref_name, ref_len)| {
            query_overlaps_for_ref(impg, *ref_id, *ref_len, config, sequence_index)
        })
        .collect();

    // Phase 2: Sequential processing with deduplication
    info!("Phase 2: Processing depth windows with deduplication...");
    let mut processed_regions = ProcessedRegions::new();
    let mut all_windows: Vec<DepthWindow> = Vec::new();
    let mut global_window_id = 0;

    for (i, (ref_id, ref_name, ref_len)) in seq_order.iter().enumerate() {
        let mut windows = process_overlaps_to_windows(
            impg,
            *ref_id,
            ref_name,
            *ref_len,
            &all_overlaps[i],
            config,
            &mut processed_regions,
            separator,
        );

        // Update global window IDs
        for w in &mut windows {
            w.window_id = global_window_id;
            global_window_id += 1;
        }

        all_windows.extend(windows);
    }

    info!("Generated {} depth windows", all_windows.len());

    // Output results
    output_depth_tsv(&all_windows, &sample_list, output_prefix)?;

    Ok(())
}

/// Output depth windows as TSV
fn output_depth_tsv(
    windows: &[DepthWindow],
    sample_list: &[String],
    output_prefix: Option<&str>,
) -> io::Result<()> {
    let writer: Box<dyn Write> = if let Some(prefix) = output_prefix {
        let path = format!("{}.depth.tsv", prefix);
        Box::new(BufWriter::new(std::fs::File::create(&path)?))
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };
    let mut writer = BufWriter::new(writer);

    // Write header
    write!(writer, "window_id\tdepth")?;
    for sample in sample_list {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // Write data rows
    for window in windows {
        write!(writer, "{}\t{}", window.window_id, window.depth)?;

        for sample in sample_list {
            if let Some((start, end)) = window.haplotype_positions.get(sample) {
                write!(writer, "\t{}:{}-{}", window.ref_name, start, end)?;
            } else {
                write!(writer, "\tNA")?;
            }
        }
        writeln!(writer)?;
    }

    writer.flush()?;
    Ok(())
}

/// Load reference order from file
pub fn load_ref_order_file(path: &str) -> io::Result<Vec<String>> {
    let content = std::fs::read_to_string(path)?;
    Ok(content
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.trim().starts_with('#'))
        .map(|line| line.trim().to_string())
        .collect())
}
