use crate::impg::{AdjustedInterval, CigarOp, Impg};
use crate::sequence_index::UnifiedSequenceIndex;
use bitvec::prelude::*;
use coitrees::IntervalTree;
use log::{debug, info};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::{self, BufRead, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};

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

// ============================================================================
// Sample filtering for depth computation
// ============================================================================

/// Filter for samples to include in depth calculation
#[derive(Debug, Clone)]
pub struct SampleFilter {
    /// Set of sample names to include (empty = include all)
    samples: FxHashSet<String>,
}

impl SampleFilter {
    /// Create a new filter with no restrictions (all samples included)
    pub fn new() -> Self {
        Self {
            samples: FxHashSet::default(),
        }
    }

    /// Create from a list of sample names
    pub fn from_samples(samples: Vec<String>) -> Self {
        Self {
            samples: samples.into_iter().collect(),
        }
    }

    /// Load from a file (one sample name per line)
    pub fn from_file(path: &str) -> io::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let samples: FxHashSet<String> = content
            .lines()
            .map(|line| line.trim())
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(|s| s.to_string())
            .collect();
        Ok(Self { samples })
    }

    /// Check if a sample should be included
    pub fn includes(&self, sample: &str) -> bool {
        self.samples.is_empty() || self.samples.contains(sample)
    }

    /// Check if the filter is active (has restrictions)
    pub fn is_active(&self) -> bool {
        !self.samples.is_empty()
    }

    /// Get the number of samples in the filter
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Check if the filter is empty (no restrictions)
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Get the set of sample names
    pub fn get_samples(&self) -> &FxHashSet<String> {
        &self.samples
    }
}

impl Default for SampleFilter {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Depth statistics for stats mode
// ============================================================================

/// Statistics for depth distribution across all sequences
#[derive(Debug, Clone, Default)]
pub struct DepthStats {
    /// Total bases analyzed
    pub total_bases: i64,
    /// Distribution: depth -> total bases at that depth
    pub depth_distribution: FxHashMap<usize, i64>,
    /// Per-depth intervals for output files: depth -> Vec<(seq_name, start, end)>
    pub depth_intervals: FxHashMap<usize, Vec<(String, i64, i64)>>,
}

impl DepthStats {
    /// Create a new empty stats container
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an interval with a specific depth
    pub fn add_interval(&mut self, seq_name: &str, start: i64, end: i64, depth: usize) {
        let length = end - start;
        if length <= 0 {
            return;
        }
        self.total_bases += length;
        *self.depth_distribution.entry(depth).or_insert(0) += length;
        self.depth_intervals
            .entry(depth)
            .or_default()
            .push((seq_name.to_string(), start, end));
    }

    /// Merge another DepthStats into this one
    pub fn merge(&mut self, other: &DepthStats) {
        self.total_bases += other.total_bases;
        for (&depth, &bases) in &other.depth_distribution {
            *self.depth_distribution.entry(depth).or_insert(0) += bases;
        }
        for (depth, intervals) in &other.depth_intervals {
            self.depth_intervals
                .entry(*depth)
                .or_default()
                .extend(intervals.iter().cloned());
        }
    }

    /// Get sorted list of depths with their base counts
    pub fn get_sorted_distribution(&self) -> Vec<(usize, i64)> {
        let mut dist: Vec<_> = self
            .depth_distribution
            .iter()
            .map(|(&d, &c)| (d, c))
            .collect();
        dist.sort_by_key(|&(d, _)| d);
        dist
    }

    /// Get maximum depth
    pub fn max_depth(&self) -> usize {
        self.depth_distribution.keys().copied().max().unwrap_or(0)
    }

    /// Write summary to a writer
    pub fn write_summary<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "Total bases: {}", self.total_bases)?;
        writeln!(writer, "Depth distribution:")?;
        for (depth, bases) in self.get_sorted_distribution() {
            let pct = if self.total_bases > 0 {
                100.0 * bases as f64 / self.total_bases as f64
            } else {
                0.0
            };
            writeln!(writer, "  depth={}: {} bp ({:.2}%)", depth, bases, pct)?;
        }
        Ok(())
    }

    /// Write per-depth BED files
    pub fn write_depth_bed_files(&self, prefix: &str) -> io::Result<()> {
        for (&depth, intervals) in &self.depth_intervals {
            let path = format!("{}.depth{}.bed", prefix, depth);
            let file = std::fs::File::create(&path)?;
            let mut writer = BufWriter::new(file);
            for (seq_name, start, end) in intervals {
                writeln!(writer, "{}\t{}\t{}", seq_name, start, end)?;
            }
            writer.flush()?;
            info!("Wrote {} intervals to {}", intervals.len(), path);
        }
        Ok(())
    }
}

// ============================================================================
// Combined depth output with sample lists
// ============================================================================

/// Interval with depth and list of covering samples
#[derive(Debug, Clone)]
pub struct DepthIntervalWithSamples {
    /// Sequence name
    pub seq_name: String,
    /// Start position (0-based)
    pub start: i64,
    /// End position (exclusive)
    pub end: i64,
    /// Depth (number of unique samples)
    pub depth: usize,
    /// Sorted list of sample names covering this interval
    pub samples: Vec<String>,
}

impl DepthIntervalWithSamples {
    pub fn new(seq_name: String, start: i64, end: i64, depth: usize, samples: Vec<String>) -> Self {
        Self {
            seq_name,
            start,
            end,
            depth,
            samples,
        }
    }
}

/// Statistics with sample tracking for combined output
#[derive(Debug, Clone, Default)]
pub struct DepthStatsWithSamples {
    /// Total bases analyzed
    pub total_bases: i64,
    /// Distribution: depth -> total bases at that depth
    pub depth_distribution: FxHashMap<usize, i64>,
    /// All intervals with samples (unsorted, to be sorted at output time)
    pub intervals: Vec<DepthIntervalWithSamples>,
}

impl DepthStatsWithSamples {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an interval with depth and sample list
    pub fn add_interval(&mut self, seq_name: &str, start: i64, end: i64, depth: usize, samples: Vec<String>) {
        let length = end - start;
        if length <= 0 {
            return;
        }
        self.total_bases += length;
        *self.depth_distribution.entry(depth).or_insert(0) += length;
        self.intervals.push(DepthIntervalWithSamples::new(
            seq_name.to_string(),
            start,
            end,
            depth,
            samples,
        ));
    }

    /// Merge another DepthStatsWithSamples into this one
    pub fn merge(&mut self, other: DepthStatsWithSamples) {
        self.total_bases += other.total_bases;
        for (&depth, &bases) in &other.depth_distribution {
            *self.depth_distribution.entry(depth).or_insert(0) += bases;
        }
        self.intervals.extend(other.intervals);
    }

    /// Get sorted list of depths with their base counts
    pub fn get_sorted_distribution(&self) -> Vec<(usize, i64)> {
        let mut dist: Vec<_> = self
            .depth_distribution
            .iter()
            .map(|(&d, &c)| (d, c))
            .collect();
        dist.sort_by_key(|&(d, _)| d);
        dist
    }

    /// Get maximum depth
    pub fn max_depth(&self) -> usize {
        self.depth_distribution.keys().copied().max().unwrap_or(0)
    }

    /// Write summary to a writer
    pub fn write_summary<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "Total bases: {}", self.total_bases)?;
        writeln!(writer, "Depth distribution:")?;
        for (depth, bases) in self.get_sorted_distribution() {
            let pct = if self.total_bases > 0 {
                100.0 * bases as f64 / self.total_bases as f64
            } else {
                0.0
            };
            writeln!(writer, "  depth={}: {} bp ({:.2}%)", depth, bases, pct)?;
        }
        Ok(())
    }

    /// Write combined output file sorted by chromosome and position
    /// If merge_tolerance > 0, adjacent intervals with depth difference within tolerance are merged
    /// tolerance is a fraction (e.g., 0.05 = 5%)
    pub fn write_combined_output(&mut self, prefix: &str, merge_tolerance: f64) -> io::Result<()> {
        let path = format!("{}.combined.bed", prefix);
        let file = std::fs::File::create(&path)?;
        let mut writer = BufWriter::new(file);

        // Sort intervals by sequence name, then by start position
        self.intervals.sort_by(|a, b| {
            a.seq_name.cmp(&b.seq_name)
                .then_with(|| a.start.cmp(&b.start))
                .then_with(|| a.end.cmp(&b.end))
        });

        // Write header
        writeln!(writer, "#seq_name\tstart\tend\tdepth\tsamples")?;

        if merge_tolerance > 0.0 && !self.intervals.is_empty() {
            // Merge adjacent intervals within tolerance
            let merged = self.merge_intervals_with_tolerance(merge_tolerance);
            info!(
                "Merged {} intervals into {} intervals (tolerance: {:.1}%)",
                self.intervals.len(),
                merged.len(),
                merge_tolerance * 100.0
            );

            for interval in &merged {
                let samples_str = interval.samples.join(",");
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    interval.seq_name, interval.start, interval.end, interval.depth, samples_str
                )?;
            }

            writer.flush()?;
            info!(
                "Wrote {} intervals to {} (combined output with samples)",
                merged.len(),
                path
            );
        } else {
            // No merging, write as-is
            for interval in &self.intervals {
                let samples_str = interval.samples.join(",");
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    interval.seq_name, interval.start, interval.end, interval.depth, samples_str
                )?;
            }

            writer.flush()?;
            info!(
                "Wrote {} intervals to {} (combined output with samples)",
                self.intervals.len(),
                path
            );
        }

        Ok(())
    }

    /// Merge adjacent intervals where depth difference is within tolerance
    /// tolerance: fraction (e.g., 0.05 means 5%)
    /// When merging: depth = max, samples = union
    /// Tracks full min-max range of merged group: (max - min) / max <= tolerance
    fn merge_intervals_with_tolerance(&self, tolerance: f64) -> Vec<DepthIntervalWithSamples> {
        if self.intervals.is_empty() {
            return Vec::new();
        }

        let mut merged: Vec<DepthIntervalWithSamples> = Vec::new();

        // Track current merge group with min and max depth
        let mut current_seq = self.intervals[0].seq_name.clone();
        let mut current_start = self.intervals[0].start;
        let mut current_end = self.intervals[0].end;
        let mut current_min_depth = self.intervals[0].depth;
        let mut current_max_depth = self.intervals[0].depth;
        let mut current_samples: FxHashSet<String> = self.intervals[0].samples.iter().cloned().collect();

        for interval in self.intervals.iter().skip(1) {
            // Check if this interval can be merged with current group
            // The new range after merging would be [min(current_min, new), max(current_max, new)]
            // Check if this range is within tolerance
            let can_merge = interval.seq_name == current_seq
                && interval.start == current_end  // Adjacent
                && Self::within_tolerance_range(current_min_depth, current_max_depth, interval.depth, tolerance);

            if can_merge {
                // Extend current group
                current_end = interval.end;
                current_min_depth = current_min_depth.min(interval.depth);
                current_max_depth = current_max_depth.max(interval.depth);
                current_samples.extend(interval.samples.iter().cloned());
            } else {
                // Emit current group and start new one
                let mut samples_vec: Vec<String> = current_samples.into_iter().collect();
                samples_vec.sort();
                merged.push(DepthIntervalWithSamples::new(
                    current_seq,
                    current_start,
                    current_end,
                    current_max_depth,  // Use max depth for output
                    samples_vec,
                ));

                // Start new group
                current_seq = interval.seq_name.clone();
                current_start = interval.start;
                current_end = interval.end;
                current_min_depth = interval.depth;
                current_max_depth = interval.depth;
                current_samples = interval.samples.iter().cloned().collect();
            }
        }

        // Emit final group
        let mut samples_vec: Vec<String> = current_samples.into_iter().collect();
        samples_vec.sort();
        merged.push(DepthIntervalWithSamples::new(
            current_seq,
            current_start,
            current_end,
            current_max_depth,
            samples_vec,
        ));

        merged
    }

    /// Check if adding a new depth to an existing range would keep it within tolerance
    /// Returns true if (new_max - new_min) / new_max <= tolerance
    /// where new_min = min(current_min, new_depth) and new_max = max(current_max, new_depth)
    fn within_tolerance_range(current_min: usize, current_max: usize, new_depth: usize, tolerance: f64) -> bool {
        let new_min = current_min.min(new_depth);
        let new_max = current_max.max(new_depth);

        if new_max == 0 {
            return true;  // All zeros, can merge
        }

        let diff = (new_max - new_min) as f64 / new_max as f64;
        diff <= tolerance
    }
}

/// Compute depth using sweep-line algorithm, tracking which samples cover each window
fn compute_sweep_line_depth_with_samples(
    _anchor_seq: &str,
    seq_len: i64,
    sample_intervals: &[(i64, i64, String)],
) -> Vec<(i64, i64, usize, Vec<String>)> {
    if sample_intervals.is_empty() {
        return vec![(0, seq_len, 0, Vec::new())];
    }

    // Create events: (position, is_start, sample_name)
    let mut events: Vec<(i64, bool, &str)> = Vec::with_capacity(sample_intervals.len() * 2);
    for (start, end, sample) in sample_intervals {
        events.push((*start, true, sample.as_str()));
        events.push((*end, false, sample.as_str()));
    }

    // Sort by position, with end events before start events at same position
    events.sort_by(|a, b| {
        a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1))
    });

    let mut results: Vec<(i64, i64, usize, Vec<String>)> = Vec::new();
    let mut active_samples: FxHashMap<&str, usize> = FxHashMap::default();
    let mut prev_pos: Option<i64> = None;
    let mut prev_sample_set: Vec<String> = Vec::new();

    for (pos, is_start, sample) in events {
        if let Some(prev) = prev_pos {
            if pos > prev {
                // Emit window from prev to pos with current depth and samples
                let depth = active_samples.len();
                if depth > 0 {
                    results.push((prev, pos, depth, prev_sample_set.clone()));
                }
            }
        }

        // Update active samples
        if is_start {
            *active_samples.entry(sample).or_insert(0) += 1;
        } else if let Some(count) = active_samples.get_mut(sample) {
            *count -= 1;
            if *count == 0 {
                active_samples.remove(sample);
            }
        }

        // Update sample set for next window
        prev_sample_set = active_samples.keys().map(|s| s.to_string()).collect();
        prev_sample_set.sort();  // Keep sorted for consistent output
        prev_pos = Some(pos);
    }

    // Merge adjacent windows with same depth AND same sample set
    if results.len() <= 1 {
        return results;
    }

    let mut merged: Vec<(i64, i64, usize, Vec<String>)> = Vec::with_capacity(results.len());
    let mut current = results[0].clone();

    for (start, end, depth, samples) in results.into_iter().skip(1) {
        if current.1 == start && current.2 == depth && current.3 == samples {
            // Extend current window
            current.1 = end;
        } else {
            merged.push(current);
            current = (start, end, depth, samples);
        }
    }
    merged.push(current);

    merged
}

// ============================================================================
// Region query result for region query mode
// ============================================================================

/// Result of a region query with per-sample position tracking
#[derive(Debug, Clone)]
pub struct RegionDepthResult {
    /// Reference sequence name
    pub ref_seq: String,
    /// Reference start position
    pub ref_start: i64,
    /// Reference end position
    pub ref_end: i64,
    /// Depth (number of samples covering this region)
    pub depth: usize,
    /// Per-sample positions: sample -> Vec<(seq_name, start, end)>
    /// Multiple alignments per sample are possible
    pub sample_positions: FxHashMap<String, Vec<(String, i64, i64)>>,
}

impl RegionDepthResult {
    /// Create a new result
    pub fn new(ref_seq: String, ref_start: i64, ref_end: i64) -> Self {
        Self {
            ref_seq,
            ref_start,
            ref_end,
            depth: 0,
            sample_positions: FxHashMap::default(),
        }
    }

    /// Add a sample position
    pub fn add_sample_position(&mut self, sample: &str, seq_name: &str, start: i64, end: i64) {
        self.sample_positions
            .entry(sample.to_string())
            .or_default()
            .push((seq_name.to_string(), start, end));
    }

    /// Update depth based on sample_positions
    pub fn update_depth(&mut self) {
        self.depth = self.sample_positions.len();
    }

    /// Format sample positions as semicolon-separated strings
    pub fn format_sample_positions(&self, sample: &str) -> String {
        match self.sample_positions.get(sample) {
            Some(positions) if !positions.is_empty() => positions
                .iter()
                .map(|(seq, start, end)| format!("{}:{}-{}", seq, start, end))
                .collect::<Vec<_>>()
                .join(";"),
            _ => "NA".to_string(),
        }
    }
}

/// Write region depth results in tabular format
pub fn write_region_depth_output<W: Write>(
    writer: &mut W,
    results: &[RegionDepthResult],
    sample_order: &[String],
) -> io::Result<()> {
    // Write header
    write!(writer, "#ref_seq\tref_start\tref_end\tdepth")?;
    for sample in sample_order {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // Write data rows
    for result in results {
        write!(
            writer,
            "{}\t{}\t{}\t{}",
            result.ref_seq, result.ref_start, result.ref_end, result.depth
        )?;
        for sample in sample_order {
            write!(writer, "\t{}", result.format_sample_positions(sample))?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

// ============================================================================
// New data structures for sample-based depth computation
// ============================================================================

/// Stores sequence lengths and sample-to-sequence mappings
#[derive(Debug, Clone)]
pub struct SequenceLengths {
    /// seq_name -> length
    pub lengths: FxHashMap<String, i64>,
    /// sample_name -> list of seq_names belonging to this sample
    pub sample_to_seqs: FxHashMap<String, Vec<String>>,
    /// seq_name -> sample_name (reverse lookup)
    pub seq_to_sample: FxHashMap<String, String>,
}

impl SequenceLengths {
    /// Build from impg sequence index
    pub fn from_impg(impg: &Impg, separator: &str) -> Self {
        let mut lengths = FxHashMap::default();
        let mut sample_to_seqs: FxHashMap<String, Vec<String>> = FxHashMap::default();
        let mut seq_to_sample = FxHashMap::default();

        for id in 0..impg.seq_index.len() as u32 {
            if let (Some(name), Some(len)) = (
                impg.seq_index.get_name(id),
                impg.seq_index.get_len_from_id(id),
            ) {
                let sample = extract_sample(name, separator);
                lengths.insert(name.to_string(), len as i64);
                sample_to_seqs
                    .entry(sample.clone())
                    .or_default()
                    .push(name.to_string());
                seq_to_sample.insert(name.to_string(), sample);
            }
        }

        // Sort sequences within each sample for deterministic ordering
        for seqs in sample_to_seqs.values_mut() {
            seqs.sort();
        }

        SequenceLengths {
            lengths,
            sample_to_seqs,
            seq_to_sample,
        }
    }

    /// Build from FAI file list
    pub fn from_fai_list(fai_list_path: &str, separator: &str) -> io::Result<Self> {
        let mut lengths = FxHashMap::default();
        let mut sample_to_seqs: FxHashMap<String, Vec<String>> = FxHashMap::default();
        let mut seq_to_sample = FxHashMap::default();

        let content = std::fs::read_to_string(fai_list_path).map_err(|e| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Failed to read FAI list file '{}': {}", fai_list_path, e),
            )
        })?;

        for line in content.lines() {
            let fai_path = line.trim();
            if fai_path.is_empty() || fai_path.starts_with('#') {
                continue;
            }

            let fai_content = std::fs::read_to_string(fai_path).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("Failed to read FAI file '{}': {}", fai_path, e),
                )
            })?;

            for fai_line in fai_content.lines() {
                let parts: Vec<&str> = fai_line.split('\t').collect();
                if parts.len() >= 2 {
                    let seq_name = parts[0].to_string();
                    let seq_len: i64 = parts[1].parse().unwrap_or(0);

                    let sample = extract_sample(&seq_name, separator);
                    lengths.insert(seq_name.clone(), seq_len);
                    sample_to_seqs
                        .entry(sample.clone())
                        .or_default()
                        .push(seq_name.clone());
                    seq_to_sample.insert(seq_name, sample);
                }
            }
        }

        // Sort sequences within each sample for deterministic ordering
        for seqs in sample_to_seqs.values_mut() {
            seqs.sort();
        }

        Ok(SequenceLengths {
            lengths,
            sample_to_seqs,
            seq_to_sample,
        })
    }

    /// Get all samples
    pub fn get_samples(&self) -> Vec<String> {
        let mut samples: Vec<_> = self.sample_to_seqs.keys().cloned().collect();
        samples.sort();
        samples
    }

    /// Get sequences for a sample
    pub fn get_seqs_of_sample(&self, sample: &str) -> &[String] {
        self.sample_to_seqs
            .get(sample)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Get length of a sequence
    pub fn get_length(&self, seq_name: &str) -> Option<i64> {
        self.lengths.get(seq_name).copied()
    }
}

/// High-performance interval set with efficient intersection and subtraction
#[derive(Debug, Clone, Default)]
pub struct IntervalSet {
    /// Sorted, non-overlapping intervals [(start, end), ...]
    intervals: Vec<(i64, i64)>,
    /// Total covered length (for quick empty check)
    total_length: i64,
}

impl IntervalSet {
    /// Create a new empty interval set
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
            total_length: 0,
        }
    }

    /// Create interval set with a single interval
    pub fn new_single(start: i64, end: i64) -> Self {
        if start >= end {
            return Self::new();
        }
        Self {
            intervals: vec![(start, end)],
            total_length: end - start,
        }
    }

    /// Add an interval to this set (merging with existing if overlapping)
    pub fn add(&mut self, start: i64, end: i64) {
        if start >= end {
            return;
        }

        if self.intervals.is_empty() {
            self.intervals.push((start, end));
            self.total_length = end - start;
            return;
        }

        // Find insertion point and merge overlapping intervals
        let mut new_start = start;
        let mut new_end = end;
        let mut to_remove = Vec::new();

        for (i, &(iv_start, iv_end)) in self.intervals.iter().enumerate() {
            // Check if intervals overlap or are adjacent
            if iv_end >= new_start && iv_start <= new_end {
                // Merge
                new_start = new_start.min(iv_start);
                new_end = new_end.max(iv_end);
                to_remove.push(i);
            }
        }

        // Calculate new total length
        let added_length = new_end - new_start;
        let removed_length: i64 = to_remove
            .iter()
            .map(|&i| self.intervals[i].1 - self.intervals[i].0)
            .sum();

        // Remove merged intervals (in reverse order to preserve indices)
        for &i in to_remove.iter().rev() {
            self.intervals.remove(i);
        }

        // Insert new merged interval
        let insert_pos = self
            .intervals
            .iter()
            .position(|&(s, _)| s > new_start)
            .unwrap_or(self.intervals.len());
        self.intervals.insert(insert_pos, (new_start, new_end));

        self.total_length = self.total_length - removed_length + added_length;
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.total_length == 0
    }

    /// Get total covered length
    pub fn total_length(&self) -> i64 {
        self.total_length
    }

    /// Get number of intervals
    pub fn interval_count(&self) -> usize {
        self.intervals.len()
    }

    /// Get all intervals
    pub fn intervals(&self) -> &[(i64, i64)] {
        &self.intervals
    }

    /// Intersect with a query interval, return intersecting portions
    pub fn intersect(&self, query_start: i64, query_end: i64) -> Vec<(i64, i64)> {
        if query_start >= query_end || self.intervals.is_empty() {
            return Vec::new();
        }

        let mut result = Vec::new();

        // Binary search for the first interval that might overlap
        let start_idx = self
            .intervals
            .partition_point(|&(_, end)| end <= query_start);

        for &(iv_start, iv_end) in &self.intervals[start_idx..] {
            if iv_start >= query_end {
                break;
            }

            let intersect_start = iv_start.max(query_start);
            let intersect_end = iv_end.min(query_end);

            if intersect_start < intersect_end {
                result.push((intersect_start, intersect_end));
            }
        }

        result
    }

    /// Subtract a single interval from this set
    pub fn subtract(&mut self, sub_start: i64, sub_end: i64) {
        if sub_start >= sub_end || self.intervals.is_empty() {
            return;
        }

        let mut new_intervals = Vec::new();
        let mut removed_length: i64 = 0;

        for &(iv_start, iv_end) in &self.intervals {
            if iv_end <= sub_start || iv_start >= sub_end {
                // No overlap, keep as is
                new_intervals.push((iv_start, iv_end));
            } else {
                // There's overlap
                if iv_start < sub_start {
                    // Keep left portion
                    new_intervals.push((iv_start, sub_start));
                }
                if iv_end > sub_end {
                    // Keep right portion
                    new_intervals.push((sub_end, iv_end));
                }
                // Calculate removed length
                let overlap_start = iv_start.max(sub_start);
                let overlap_end = iv_end.min(sub_end);
                removed_length += overlap_end - overlap_start;
            }
        }

        self.intervals = new_intervals;
        self.total_length -= removed_length;
    }

    /// Batch subtract multiple intervals (more efficient than individual subtracts)
    pub fn subtract_batch(&mut self, to_remove: &[(i64, i64)]) {
        if to_remove.is_empty() || self.intervals.is_empty() {
            return;
        }

        // Merge overlapping intervals to remove
        let merged_remove = merge_intervals(to_remove);

        let mut new_intervals = Vec::new();
        let mut remove_idx = 0;
        let mut removed_length: i64 = 0;

        for &(iv_start, iv_end) in &self.intervals {
            let mut current_start = iv_start;

            while remove_idx < merged_remove.len() && current_start < iv_end {
                let (rm_start, rm_end) = merged_remove[remove_idx];

                if rm_end <= current_start {
                    // Remove interval is before current position
                    remove_idx += 1;
                    continue;
                }

                if rm_start >= iv_end {
                    // Remove interval is after current interval
                    break;
                }

                // There's overlap
                if rm_start > current_start {
                    // Keep [current_start, rm_start)
                    new_intervals.push((current_start, rm_start));
                }

                // Calculate removed portion
                let overlap_start = current_start.max(rm_start);
                let overlap_end = iv_end.min(rm_end);
                if overlap_end > overlap_start {
                    removed_length += overlap_end - overlap_start;
                }

                current_start = rm_end;

                if rm_end <= iv_end {
                    remove_idx += 1;
                }
            }

            if current_start < iv_end {
                new_intervals.push((current_start, iv_end));
            }
        }

        self.intervals = new_intervals;
        self.total_length -= removed_length;
    }
}

/// Merge overlapping/adjacent intervals
fn merge_intervals(intervals: &[(i64, i64)]) -> Vec<(i64, i64)> {
    if intervals.is_empty() {
        return Vec::new();
    }

    let mut sorted: Vec<_> = intervals.to_vec();
    sorted.sort_by_key(|&(s, _)| s);

    let mut merged = Vec::new();
    let mut current = sorted[0];

    for &(start, end) in &sorted[1..] {
        if start <= current.1 {
            // Overlapping or adjacent
            current.1 = current.1.max(end);
        } else {
            merged.push(current);
            current = (start, end);
        }
    }
    merged.push(current);

    merged
}

/// Pool: tracks unprocessed regions for all samples
#[derive(Debug)]
pub struct Pool {
    /// sample_name -> seq_name -> interval set of unprocessed regions
    regions: FxHashMap<String, FxHashMap<String, IntervalSet>>,
}

impl Pool {
    /// Initialize pool with full-length intervals for all sequences
    pub fn new(seq_lengths: &SequenceLengths) -> Self {
        let mut regions = FxHashMap::default();

        for (sample, seqs) in &seq_lengths.sample_to_seqs {
            let mut sample_regions = FxHashMap::default();
            for seq in seqs {
                if let Some(&len) = seq_lengths.lengths.get(seq) {
                    sample_regions.insert(seq.clone(), IntervalSet::new_single(0, len));
                }
            }
            regions.insert(sample.clone(), sample_regions);
        }

        Pool { regions }
    }

    /// Get remaining intervals for a specific sequence
    pub fn get_intervals(&self, sample: &str, seq: &str) -> Vec<(i64, i64)> {
        self.regions
            .get(sample)
            .and_then(|s| s.get(seq))
            .map(|is| is.intervals().to_vec())
            .unwrap_or_default()
    }

    /// Get all remaining intervals for a sample: Vec<(seq_name, start, end)>
    pub fn get_sample_intervals(&self, sample: &str) -> Vec<(String, i64, i64)> {
        let mut result = Vec::new();
        if let Some(seqs) = self.regions.get(sample) {
            for (seq, interval_set) in seqs {
                for &(start, end) in interval_set.intervals() {
                    result.push((seq.clone(), start, end));
                }
            }
        }
        result
    }

    /// Count total intervals for a sample
    pub fn count_intervals(&self, sample: &str) -> usize {
        self.regions
            .get(sample)
            .map(|seqs| seqs.values().map(|is| is.interval_count()).sum())
            .unwrap_or(0)
    }

    /// Check if a sample has no remaining intervals
    pub fn is_sample_empty(&self, sample: &str) -> bool {
        self.regions
            .get(sample)
            .map(|seqs| seqs.values().all(|is| is.is_empty()))
            .unwrap_or(true)
    }

    /// Check if entire pool is empty
    pub fn is_empty(&self) -> bool {
        self.regions
            .values()
            .all(|seqs| seqs.values().all(|is| is.is_empty()))
    }

    /// Subtract a single region
    pub fn subtract(&mut self, sample: &str, seq: &str, start: i64, end: i64) {
        if let Some(seqs) = self.regions.get_mut(sample) {
            if let Some(interval_set) = seqs.get_mut(seq) {
                interval_set.subtract(start, end);
            }
        }
    }

    /// Check if a region has any overlap with the pool for a sample/seq
    /// Returns the intersection intervals if any
    pub fn get_overlap(&self, sample: &str, seq: &str, start: i64, end: i64) -> Vec<(i64, i64)> {
        if let Some(seqs) = self.regions.get(sample) {
            if let Some(interval_set) = seqs.get(seq) {
                return interval_set.intersect(start, end);
            }
        }
        Vec::new()
    }

    /// Check if a region has any overlap with the pool
    pub fn has_overlap(&self, sample: &str, seq: &str, start: i64, end: i64) -> bool {
        !self.get_overlap(sample, seq, start, end).is_empty()
    }

    /// Batch subtract multiple regions: Vec<(sample, seq, start, end)>
    pub fn subtract_batch(&mut self, regions_to_remove: &[(String, String, i64, i64)]) {
        // Group by (sample, seq) for efficiency
        let mut grouped: FxHashMap<(&str, &str), Vec<(i64, i64)>> = FxHashMap::default();

        for (sample, seq, start, end) in regions_to_remove {
            grouped
                .entry((sample.as_str(), seq.as_str()))
                .or_default()
                .push((*start, *end));
        }

        // Batch update each (sample, seq)
        for ((sample, seq), intervals) in grouped {
            if let Some(seqs) = self.regions.get_mut(sample) {
                if let Some(interval_set) = seqs.get_mut(seq) {
                    interval_set.subtract_batch(&intervals);
                }
            }
        }
    }

    /// Get total remaining length across all samples
    pub fn total_remaining_length(&self) -> i64 {
        self.regions
            .values()
            .flat_map(|seqs| seqs.values())
            .map(|is| is.total_length())
            .sum()
    }
}

/// Reverse index: maps query sequences to their alignments
/// This enables efficient lookup of alignments where a sequence is the query
#[derive(Debug, Default)]
pub struct ReverseAlignmentIndex {
    /// seq_name -> Vec<(target_seq_name, target_start, target_end, query_start, query_end, is_reverse)>
    by_query_seq: FxHashMap<String, Vec<ReverseAlignmentEntry>>,
}

#[derive(Debug, Clone)]
pub struct ReverseAlignmentEntry {
    pub target_seq: String,
    pub target_start: i64,
    pub target_end: i64,
    pub query_start: i64,
    pub query_end: i64,
    pub is_reverse: bool,
}

impl ReverseAlignmentIndex {
    /// Build reverse index from impg
    pub fn build(impg: &Impg) -> Self {
        let mut by_query_seq: FxHashMap<String, Vec<ReverseAlignmentEntry>> = FxHashMap::default();

        // Iterate through all target sequences' interval trees
        let target_ids: Vec<u32> = impg.forest_map.entries.keys().copied().collect();

        for target_id in target_ids {
            let target_name = match impg.seq_index.get_name(target_id) {
                Some(name) => name.to_string(),
                None => continue,
            };

            // Get the actual tree
            let tree = match impg.get_or_load_tree(target_id) {
                Some(t) => t,
                None => continue,
            };

            // Query all intervals in this tree
            let target_len = impg.seq_index.get_len_from_id(target_id).unwrap_or(0);
            tree.query(0, target_len as i32, |interval| {
                let query_id = interval.metadata.query_id();
                if let Some(query_name) = impg.seq_index.get_name(query_id) {
                    let q_start = interval.metadata.query_start();
                    let q_end = interval.metadata.query_end();
                    let is_reverse = q_start > q_end;
                    let query_start = q_start.min(q_end) as i64;
                    let query_end = q_start.max(q_end) as i64;

                    by_query_seq
                        .entry(query_name.to_string())
                        .or_default()
                        .push(ReverseAlignmentEntry {
                            target_seq: target_name.clone(),
                            target_start: interval.first as i64,
                            target_end: interval.last as i64,
                            query_start,
                            query_end,
                            is_reverse,
                        });
                }
            });
        }

        ReverseAlignmentIndex { by_query_seq }
    }

    /// Query alignments where the given sequence is a query, filtering by range
    pub fn query_range(
        &self,
        query_seq: &str,
        range_start: i64,
        range_end: i64,
    ) -> Vec<&ReverseAlignmentEntry> {
        self.by_query_seq
            .get(query_seq)
            .map(|entries| {
                entries
                    .iter()
                    .filter(|e| {
                        // Check if alignment overlaps with query range
                        e.query_start < range_end && e.query_end > range_start
                    })
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get all alignments where the given sequence is a query
    pub fn get_all(&self, query_seq: &str) -> Option<&Vec<ReverseAlignmentEntry>> {
        self.by_query_seq.get(query_seq)
    }

    /// Count total alignments indexed
    pub fn total_alignments(&self) -> usize {
        self.by_query_seq.values().map(|v| v.len()).sum()
    }
}

/// Strategy for querying alignments
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueryStrategy {
    /// First round: collect all alignments for the sample from index
    FromIndex,
    /// Subsequent rounds: query from remaining pool intervals
    FromPool,
}

/// Result from processing a single sequence's alignments
#[derive(Debug)]
pub struct SeqProcessResult {
    /// Generated depth windows
    pub windows: Vec<DepthWindow>,
    /// Regions to remove from pool: (sample, seq, start, end)
    pub regions_to_remove: Vec<(String, String, i64, i64)>,
}

/// Result from processing a sample
#[derive(Debug)]
pub struct SampleProcessResult {
    /// All windows from this sample
    pub windows: Vec<DepthWindow>,
    /// All regions to remove from pool
    pub regions_to_remove: Vec<(String, String, i64, i64)>,
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

/// Tracks processed regions in query space for each (sample, query_seq) pair
#[derive(Default)]
pub struct QueryProcessedRegions {
    /// (sample, query_seq_name) -> sorted list of (start, end) intervals
    regions: FxHashMap<(String, String), Vec<(i64, i64)>>,
}

impl QueryProcessedRegions {
    pub fn new() -> Self {
        Self {
            regions: FxHashMap::default(),
        }
    }

    /// Get ALL unprocessed portions of a query region
    /// Returns a vector of (start, end) intervals that haven't been processed yet
    pub fn get_all_unprocessed(
        &self,
        sample: &str,
        query_seq: &str,
        start: i64,
        end: i64,
    ) -> Vec<(i64, i64)> {
        let key = (sample.to_string(), query_seq.to_string());
        if let Some(intervals) = self.regions.get(&key) {
            let mut result = Vec::new();
            let mut current_start = start;

            for &(proc_start, proc_end) in intervals {
                if proc_end <= current_start {
                    continue;
                }
                if proc_start >= end {
                    break;
                }

                if proc_start > current_start {
                    // There's an unprocessed gap before this processed region
                    result.push((current_start, proc_start.min(end)));
                }
                // Move past this processed region
                current_start = proc_end;
                if current_start >= end {
                    return result;
                }
            }

            // Don't forget the tail portion after all processed regions
            if current_start < end {
                result.push((current_start, end));
            }
            result
        } else {
            vec![(start, end)]
        }
    }

    /// Mark a query region as processed
    pub fn mark_processed(&mut self, sample: &str, query_seq: &str, start: i64, end: i64) {
        let key = (sample.to_string(), query_seq.to_string());
        let intervals = self.regions.entry(key).or_default();

        let mut new_start = start;
        let mut new_end = end;
        let mut merged_indices = Vec::new();

        for (i, &(s, e)) in intervals.iter().enumerate() {
            if e < new_start {
                continue;
            }
            if s > new_end {
                break;
            }
            new_start = new_start.min(s);
            new_end = new_end.max(e);
            merged_indices.push(i);
        }

        for i in merged_indices.into_iter().rev() {
            intervals.remove(i);
        }

        let pos = intervals
            .iter()
            .position(|&(s, _)| s > new_start)
            .unwrap_or(intervals.len());
        intervals.insert(pos, (new_start, new_end));
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
    /// sample_name -> (query_seq_name, query_start, query_end) position on the query sequence
    pub haplotype_positions: FxHashMap<String, (String, i64, i64)>,
}

/// Information about a sample's alignment (for tracking query coordinates)
#[derive(Debug, Clone)]
struct SampleAlignmentInfo {
    sample: String,
    query_name: String,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    /// True if the alignment is on the reverse strand (query.first > query.last)
    is_reverse: bool,
    /// CIGAR operations for precise coordinate mapping
    cigar: Vec<CigarOp>,
}

/// Map a target coordinate range to query coordinates using linear interpolation.
/// This provides smooth coordinate transitions between adjacent windows.
///
/// Note: CIGAR-precise mapping was considered but produces more discontinuities
/// because it reveals actual alignment structure differences. Linear interpolation
/// smooths these over at the cost of some coordinate accuracy.
#[allow(unused_variables)]
fn map_target_to_query_linear(
    cigar: &[CigarOp],
    aln_target_start: i64,
    aln_target_end: i64,
    aln_query_start: i64,
    aln_query_end: i64,
    window_target_start: i64,
    window_target_end: i64,
    is_reverse: bool,
) -> (i64, i64) {
    let target_len = aln_target_end - aln_target_start;
    let query_len = aln_query_end - aln_query_start;

    if target_len == 0 {
        return (aln_query_start, aln_query_end);
    }

    let window_target_start_clamped = window_target_start.max(aln_target_start);
    let window_target_end_clamped = window_target_end.min(aln_target_end);

    let start_offset = window_target_start_clamped - aln_target_start;
    let end_offset = window_target_end_clamped - aln_target_start;

    let ratio = query_len as f64 / target_len as f64;

    if is_reverse {
        let q_end = aln_query_end - (start_offset as f64 * ratio) as i64;
        let q_start = aln_query_end - (end_offset as f64 * ratio) as i64;
        (q_start.min(q_end), q_start.max(q_end))
    } else {
        let q_start = aln_query_start + (start_offset as f64 * ratio) as i64;
        let q_end = aln_query_start + (end_offset as f64 * ratio) as i64;
        (q_start.min(q_end), q_start.max(q_end))
    }
}

/// Event for sweep-line algorithm to compute depth
#[derive(Debug, Clone, PartialEq, Eq)]
struct DepthEvent {
    position: i64,
    is_start: bool,
    haplotype: String,
    /// Index into the alignment info array (for looking up query coordinates)
    alignment_idx: usize,
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
            // Format: (query_name, query_start, query_end)
            for (hap, (query_name, _, new_end)) in window.haplotype_positions {
                if let Some(pos) = current.haplotype_positions.get_mut(&hap) {
                    // Only extend if same query sequence
                    if pos.0 == query_name {
                        pos.2 = new_end;
                    }
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
                true, // Store CIGAR for precise coordinate mapping
                None,
                sequence_index,
                config.approximate_mode,
                None, // subset_filter
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
                true, // Store CIGAR for precise coordinate mapping
                None,
                sequence_index,
                config.approximate_mode,
                None, // subset_filter
            )
        }
    } else {
        impg.query(
            ref_id,
            0,
            ref_length as i32,
            true, // Store CIGAR for precise coordinate mapping
            None,
            sequence_index,
            config.approximate_mode,
        )
    }
}

/// Process overlaps into depth windows (used in sequential phase)
/// Note: Deduplication is handled in output_depth_tsv using QueryProcessedRegions
/// reverse_alignments: Vec of (ref_start, ref_end, target_id) for alignments where ref is query
fn process_overlaps_to_windows(
    impg: &Impg,
    _ref_id: u32,
    ref_name: &str,
    ref_length: i64,
    overlaps: &[AdjustedInterval],
    reverse_alignments: &[(i32, i32, u32)],
    config: &DepthConfig,
    separator: &str,
) -> Vec<DepthWindow> {
    let ref_sample = extract_sample(ref_name, separator);

    if overlaps.is_empty() && reverse_alignments.is_empty() {
        debug!("No overlaps found for reference {}", ref_name);
        return vec![DepthWindow {
            window_id: 0,
            depth: 1,
            ref_name: ref_name.to_string(),
            ref_start: 0,
            ref_end: ref_length,
            haplotype_positions: {
                let mut map = FxHashMap::default();
                map.insert(ref_sample.clone(), (ref_name.to_string(), 0, ref_length));
                map
            },
        }];
    }

    // Collect all alignment info (including query coordinates and CIGAR)
    let mut all_alignments: Vec<SampleAlignmentInfo> = Vec::new();

    // Forward direction: ref is target, other samples are queries
    for overlap in overlaps {
        let query_interval = &overlap.0;
        let cigar = &overlap.1;
        let target_interval = &overlap.2;

        let query_name = impg.seq_index.get_name(query_interval.metadata).unwrap();
        let sample = extract_sample(query_name, separator);

        // Check if reverse strand (query.first > query.last)
        let is_reverse = query_interval.first > query_interval.last;
        let query_start = query_interval.first.min(query_interval.last) as i64;
        let query_end = query_interval.first.max(query_interval.last) as i64;
        let target_start = target_interval.first.min(target_interval.last) as i64;
        let target_end = target_interval.first.max(target_interval.last) as i64;

        all_alignments.push(SampleAlignmentInfo {
            sample,
            query_name: query_name.to_string(),
            query_start,
            query_end,
            target_start,
            target_end,
            is_reverse,
            cigar: cigar.clone(),
        });
    }

    // Reverse direction: ref is query, other samples are targets
    // In these alignments, our ref was aligned TO other sequences
    // So the "covering sample" is the original target sequence
    for &(ref_start, ref_end, target_id) in reverse_alignments {
        // target_id is the sequence that our ref aligned to
        // That sequence's sample is the "covering sample"
        if let Some(target_name) = impg.seq_index.get_name(target_id) {
            let sample = extract_sample(target_name, separator);

            // Skip if same sample (shouldn't contribute to depth)
            if sample == ref_sample {
                continue;
            }

            // For reverse alignments:
            // - target_start/target_end on ref = ref_start/ref_end (the query coords)
            // - query_start/query_end = position on the covering sample's sequence (we use ref coords for simplicity)
            // Since we're computing depth on ref, the "target" coords are the ref coords
            all_alignments.push(SampleAlignmentInfo {
                sample,
                query_name: target_name.to_string(),
                query_start: ref_start as i64, // Approximation: use ref coords
                query_end: ref_end as i64,
                target_start: ref_start as i64,
                target_end: ref_end as i64,
                is_reverse: false,
                cigar: Vec::new(), // No CIGAR for reverse (linear mapping)
            });
        }
    }

    // Add self (reference) as a sample with empty CIGAR (identity mapping)
    all_alignments.push(SampleAlignmentInfo {
        sample: ref_sample.clone(),
        query_name: ref_name.to_string(),
        query_start: 0,
        query_end: ref_length,
        target_start: 0,
        target_end: ref_length,
        is_reverse: false,
        cigar: Vec::new(), // Empty CIGAR for self-alignment (identity)
    });

    // Create events for sweep-line algorithm
    // We need to track alignment indices to look up query info later
    let mut events: Vec<DepthEvent> = Vec::new();
    for (idx, aln) in all_alignments.iter().enumerate() {
        events.push(DepthEvent {
            position: aln.target_start,
            is_start: true,
            haplotype: aln.sample.clone(),
            alignment_idx: idx,
        });
        events.push(DepthEvent {
            position: aln.target_end,
            is_start: false,
            haplotype: aln.sample.clone(),
            alignment_idx: idx,
        });
    }
    events.sort();

    // Sweep-line to compute depth windows
    // Track active alignments for each sample (to get query info)
    let mut windows: Vec<DepthWindow> = Vec::new();
    let mut active_alignments: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    let mut prev_pos: Option<i64> = None;
    let mut window_id = 0;

    for event in events {
        if let Some(prev) = prev_pos {
            if event.position > prev {
                // Collect samples that are currently active (have at least one alignment)
                let active_samples: Vec<&String> = active_alignments
                    .iter()
                    .filter(|(_, alns)| !alns.is_empty())
                    .map(|(s, _)| s)
                    .collect();

                if !active_samples.is_empty() {
                    let mut sample_positions: FxHashMap<String, (String, i64, i64)> =
                        FxHashMap::default();
                    for sample in &active_samples {
                        // Get the best alignment for this sample (pick the one with largest overlap)
                        if let Some(aln_indices) = active_alignments.get(*sample) {
                            if let Some(&best_idx) = aln_indices.iter().max_by_key(|&&idx| {
                                let aln = &all_alignments[idx];
                                // Calculate overlap with current window
                                let overlap_start = prev.max(aln.target_start);
                                let overlap_end = event.position.min(aln.target_end);
                                overlap_end - overlap_start
                            }) {
                                let aln = &all_alignments[best_idx];

                                // Calculate query coordinates corresponding to this window
                                // Window on target: [prev, event.position]
                                let window_target_start = prev.max(aln.target_start);
                                let window_target_end = event.position.min(aln.target_end);

                                // Use linear interpolation for coordinate mapping
                                let (window_query_start, window_query_end) =
                                    map_target_to_query_linear(
                                        &aln.cigar,
                                        aln.target_start,
                                        aln.target_end,
                                        aln.query_start,
                                        aln.query_end,
                                        window_target_start,
                                        window_target_end,
                                        aln.is_reverse,
                                    );

                                sample_positions.insert(
                                    (*sample).clone(),
                                    (aln.query_name.clone(), window_query_start, window_query_end),
                                );
                            }
                        }
                    }

                    if !sample_positions.is_empty() {
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
        }

        // Update active alignments
        if event.is_start {
            active_alignments
                .entry(event.haplotype.clone())
                .or_default()
                .push(event.alignment_idx);
        } else if let Some(alns) = active_alignments.get_mut(&event.haplotype) {
            alns.retain(|&idx| idx != event.alignment_idx);
        }

        prev_pos = Some(event.position);
    }

    // Optionally merge adjacent windows
    if config.merge_adjacent {
        windows = merge_adjacent_windows(windows);
    }

    windows
}

/// Main function to compute depth across all sequences
/// Uses streaming architecture to minimize memory usage:
/// - Process one reference at a time
/// - Query overlaps → Process windows → Deduplicate → Output immediately
/// - Release memory before processing next reference
pub fn compute_depth(
    impg: &Impg,
    config: &DepthConfig,
    ref_order: Option<Vec<String>>,
    sequence_index: Option<&UnifiedSequenceIndex>,
    separator: &str,
    output_prefix: Option<&str>,
) -> io::Result<()> {
    info!("Computing depth coverage (streaming mode)");

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

    // Build seq_lengths map for uncovered regions
    let seq_lengths: FxHashMap<String, i64> = seq_order
        .iter()
        .map(|(_, name, len)| (name.clone(), *len))
        .collect();

    // Build string interning maps
    let sample_to_idx: FxHashMap<&str, u16> = sample_list
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i as u16))
        .collect();

    let seq_names: Vec<&String> = seq_lengths.keys().collect();
    let seq_to_idx: FxHashMap<&str, u32> = seq_names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i as u32))
        .collect();

    let idx_to_seq: Vec<&str> = {
        let mut v: Vec<(&str, u32)> = seq_to_idx.iter().map(|(&s, &i)| (s, i)).collect();
        v.sort_by_key(|(_, i)| *i);
        v.into_iter().map(|(s, _)| s).collect()
    };

    info!(
        "String interning: {} samples, {} sequences",
        sample_to_idx.len(),
        seq_to_idx.len()
    );

    // Open output file
    let writer: Box<dyn Write> = if let Some(prefix) = output_prefix {
        let path = format!("{}.depth.tsv", prefix);
        Box::new(BufWriter::new(std::fs::File::create(&path)?))
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };
    let mut writer = BufWriter::new(writer);

    // Write header
    write!(writer, "window_id\tdepth")?;
    for sample in &sample_list {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // Coverage trackers for each (sample_idx, seq_idx) pair
    // This tracks what we've already output to avoid duplicates
    let mut coverage_trackers: FxHashMap<(u16, u32), CoverageTracker> = FxHashMap::default();

    // Process in batches for parallelism while limiting memory
    // Batch size: process N references at a time
    const BATCH_SIZE: usize = 8;
    let mut global_window_id: usize = 0;
    let total_refs = seq_order.len();

    info!("Processing {} references in batches of {}", total_refs, BATCH_SIZE);

    for (batch_idx, batch) in seq_order.chunks(BATCH_SIZE).enumerate() {
        let batch_start = batch_idx * BATCH_SIZE;
        debug!(
            "Processing batch {}: refs {}-{} of {}",
            batch_idx,
            batch_start,
            batch_start + batch.len() - 1,
            total_refs
        );

        // Phase 1: Query overlaps for this batch in parallel
        let batch_overlaps: Vec<Vec<AdjustedInterval>> = batch
            .par_iter()
            .map(|(ref_id, _ref_name, ref_len)| {
                query_overlaps_for_ref(impg, *ref_id, *ref_len, config, sequence_index)
            })
            .collect();

        // Phase 2: Query reverse alignments in parallel (where ref is query in PAF)
        let batch_reverse_alignments: Vec<Vec<(i32, i32, u32)>> = batch
            .par_iter()
            .map(|(ref_id, _ref_name, _ref_len)| {
                impg.query_reverse_for_depth(*ref_id)
            })
            .collect();

        // Phase 3: Process overlaps into depth windows in parallel
        let batch_windows: Vec<Vec<DepthWindow>> = batch
            .par_iter()
            .zip(batch_overlaps.par_iter())
            .zip(batch_reverse_alignments.par_iter())
            .map(|(((ref_id, ref_name, ref_len), overlaps), reverse_alignments)| {
                process_overlaps_to_windows(
                    impg,
                    *ref_id,
                    ref_name,
                    *ref_len,
                    overlaps,
                    reverse_alignments,
                    config,
                    separator,
                )
            })
            .collect();

        // Phase 4: Deduplicate and output (sequential to maintain order)
        for mut windows in batch_windows {
            // Assign global window IDs
            for w in &mut windows {
                w.window_id = global_window_id;
                global_window_id += 1;
            }

            // Process each window and output deduplicated intervals
            for window in &windows {
                for (sample, (seq_name, start, end)) in &window.haplotype_positions {
                    if let (Some(&sample_idx), Some(&seq_idx)) = (
                        sample_to_idx.get(sample.as_str()),
                        seq_to_idx.get(seq_name.as_str()),
                    ) {
                        let key = (sample_idx, seq_idx);
                        let tracker = coverage_trackers.entry(key).or_default();

                        // Get uncovered portions
                        let uncovered = tracker.add_interval(*start, *end);

                        // Output each uncovered portion
                        for (unc_start, unc_end) in uncovered {
                            write!(writer, "{}\t{}", window.window_id, window.depth)?;
                            for (i, _) in sample_list.iter().enumerate() {
                                if i as u16 == sample_idx {
                                    write!(
                                        writer,
                                        "\t{}:{}-{}",
                                        idx_to_seq[seq_idx as usize], unc_start, unc_end
                                    )?;
                                } else {
                                    write!(writer, "\tNA")?;
                                }
                            }
                            writeln!(writer)?;
                        }
                    }
                }
            }
            // Windows for this reference are dropped here, releasing memory
        }
        // batch_overlaps and batch_windows are dropped here

        if (batch_idx + 1) % 10 == 0 || batch_idx == seq_order.chunks(BATCH_SIZE).count() - 1 {
            info!(
                "Progress: {}/{} references processed",
                (batch_idx + 1) * BATCH_SIZE.min(total_refs - batch_idx * BATCH_SIZE),
                total_refs
            );
        }
    }

    info!("Generated {} depth windows", global_window_id);

    // Phase 4: Find and output uncovered regions
    info!("Finding uncovered regions...");
    let mut uncovered_count = 0;

    for (seq_name, &seq_len) in &seq_lengths {
        let sample = extract_sample(seq_name, separator);
        if let (Some(&sample_idx), Some(&seq_idx)) = (
            sample_to_idx.get(sample.as_str()),
            seq_to_idx.get(seq_name.as_str()),
        ) {
            let key = (sample_idx, seq_idx);
            let mut current_pos: i64 = 0;

            if let Some(tracker) = coverage_trackers.get(&key) {
                for &(iv_start, iv_end) in tracker.get_intervals() {
                    if iv_start > current_pos {
                        // Output uncovered region
                        write!(writer, "uncovered\t1")?;
                        for (i, _) in sample_list.iter().enumerate() {
                            if i as u16 == sample_idx {
                                write!(
                                    writer,
                                    "\t{}:{}-{}",
                                    idx_to_seq[seq_idx as usize], current_pos, iv_start
                                )?;
                            } else {
                                write!(writer, "\tNA")?;
                            }
                        }
                        writeln!(writer)?;
                        uncovered_count += 1;
                    }
                    current_pos = current_pos.max(iv_end);
                }
            }

            // Check for uncovered region at the end
            if current_pos < seq_len {
                write!(writer, "uncovered\t1")?;
                for (i, _) in sample_list.iter().enumerate() {
                    if i as u16 == sample_idx {
                        write!(
                            writer,
                            "\t{}:{}-{}",
                            idx_to_seq[seq_idx as usize], current_pos, seq_len
                        )?;
                    } else {
                        write!(writer, "\tNA")?;
                    }
                }
                writeln!(writer)?;
                uncovered_count += 1;
            }
        }
    }

    if uncovered_count > 0 {
        info!("Found {} uncovered regions", uncovered_count);
    }

    writer.flush()?;
    info!("Output complete");

    Ok(())
}

/// Coverage tracker for a single (sample, seq) pair
/// Maintains a sorted list of non-overlapping covered intervals
#[derive(Debug, Default)]
struct CoverageTracker {
    /// Sorted, non-overlapping intervals: Vec<(start, end)>
    intervals: Vec<(i64, i64)>,
}

impl CoverageTracker {
    /// Add a new interval and return the uncovered portions
    /// Returns: Vec<(start, end)> of newly covered regions
    ///
    /// This implementation uses a simple linear scan to be robust against edge cases.
    /// It guarantees that returned intervals are non-overlapping with existing coverage.
    fn add_interval(&mut self, start: i64, end: i64) -> Vec<(i64, i64)> {
        if start >= end {
            return Vec::new();
        }

        let mut uncovered = Vec::new();
        let mut current = start;

        // Simple linear scan through all intervals
        for &(iv_start, iv_end) in &self.intervals {
            // Skip intervals that end before our current position
            if iv_end <= current {
                continue;
            }
            // Stop if interval starts at or after our end
            if iv_start >= end {
                break;
            }

            // There's some overlap or gap
            if current < iv_start {
                // Gap before this interval - this portion is uncovered
                let gap_end = iv_start.min(end);
                uncovered.push((current, gap_end));
            }

            // Move current position past this interval
            current = current.max(iv_end);

            // Early exit if we've passed our end
            if current >= end {
                break;
            }
        }

        // Final gap after all intervals
        if current < end {
            uncovered.push((current, end));
        }

        // Merge the new interval into our list if there are uncovered portions
        if !uncovered.is_empty() {
            self.merge_interval(start, end);
        }

        uncovered
    }

    /// Merge a new interval into the sorted list
    /// Maintains the invariant that intervals are sorted and non-overlapping
    fn merge_interval(&mut self, start: i64, end: i64) {
        if self.intervals.is_empty() {
            self.intervals.push((start, end));
            return;
        }

        // Find all intervals that overlap or are adjacent to [start, end]
        let mut merge_start = start;
        let mut merge_end = end;
        let mut first_overlap: Option<usize> = None;
        let mut last_overlap: Option<usize> = None;

        for (i, &(iv_start, iv_end)) in self.intervals.iter().enumerate() {
            // Check if overlaps or adjacent (iv_end >= start && iv_start <= end)
            if iv_end >= start && iv_start <= end {
                if first_overlap.is_none() {
                    first_overlap = Some(i);
                }
                last_overlap = Some(i);
                merge_start = merge_start.min(iv_start);
                merge_end = merge_end.max(iv_end);
            }
        }

        match (first_overlap, last_overlap) {
            (Some(first), Some(last)) => {
                // Remove overlapping intervals and insert merged one
                self.intervals.drain(first..=last);
                self.intervals.insert(first, (merge_start, merge_end));
            }
            _ => {
                // No overlap, insert at correct position
                let pos = self
                    .intervals
                    .iter()
                    .position(|&(s, _)| s > start)
                    .unwrap_or(self.intervals.len());
                self.intervals.insert(pos, (start, end));
            }
        }
    }

    /// Get all intervals (for final output)
    fn get_intervals(&self) -> &[(i64, i64)] {
        &self.intervals
    }
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

// ============================================================================
// New sample-based depth computation with hybrid strategy
// ============================================================================

/// Determine the optimal sample processing order
/// Samples with more alignments are processed first to maximize pool coverage early
#[allow(dead_code)]
fn determine_sample_order(
    seq_lengths: &SequenceLengths,
    impg: &Impg,
    separator: &str,
) -> Vec<String> {
    let mut sample_alignment_counts: FxHashMap<String, usize> = FxHashMap::default();

    // Count alignments per sample (as target)
    let target_ids: Vec<u32> = impg.forest_map.entries.keys().copied().collect();

    for target_id in target_ids {
        if let Some(target_name) = impg.seq_index.get_name(target_id) {
            let sample = extract_sample(target_name, separator);

            // Get the actual tree and count alignments
            if let Some(tree) = impg.get_or_load_tree(target_id) {
                let target_len = impg.seq_index.get_len_from_id(target_id).unwrap_or(0);
                let mut count = 0;
                tree.query(0, target_len as i32, |_| count += 1);
                *sample_alignment_counts.entry(sample).or_default() += count;
            }
        }
    }

    // Sort samples by alignment count (descending) for better early coverage
    let mut samples: Vec<_> = seq_lengths.get_samples();
    samples.sort_by(|a, b| {
        let count_a = sample_alignment_counts.get(a).copied().unwrap_or(0);
        let count_b = sample_alignment_counts.get(b).copied().unwrap_or(0);
        count_b.cmp(&count_a).then_with(|| a.cmp(b))
    });

    samples
}

/// Choose query strategy based on pool state and alignment density
fn choose_strategy(
    _sample: &str,
    _pool: &Pool,
    _impg: &Impg,
    _seq_lengths: &SequenceLengths,
    _is_first_round: bool,
) -> QueryStrategy {
    // Always use FromIndex to ensure ALL alignments are collected
    // This is necessary because FromPool would miss alignments for regions
    // that were removed from the pool when earlier samples were processed,
    // resulting in incomplete windows (missing samples)
    // The pool filtering in window generation prevents duplicate output
    QueryStrategy::FromIndex
}

/// Collect alignments for a sequence using FromIndex strategy
fn collect_alignments_from_index(
    impg: &Impg,
    anchor_seq: &str,
    anchor_sample: &str,
    reverse_index: &ReverseAlignmentIndex,
    config: &DepthConfig,
    sequence_index: Option<&UnifiedSequenceIndex>,
    separator: &str,
) -> Vec<SampleAlignmentInfo> {
    let mut alignments = Vec::new();

    let seq_id = match impg.seq_index.get_id(anchor_seq) {
        Some(id) => id,
        None => return alignments,
    };
    let seq_len = impg.seq_index.get_len_from_id(seq_id).unwrap_or(0) as i64;

    // 1. Query where anchor_seq is TARGET (other samples aligned to this)
    let overlaps = if config.transitive {
        if config.transitive_dfs {
            impg.query_transitive_dfs(
                seq_id,
                0,
                seq_len as i32,
                None,
                config.max_depth,
                config.min_transitive_len,
                config.min_distance_between_ranges,
                None,
                false, // Don't need CIGAR for now
                None,
                sequence_index,
                config.approximate_mode,
                None,
            )
        } else {
            impg.query_transitive_bfs(
                seq_id,
                0,
                seq_len as i32,
                None,
                config.max_depth,
                config.min_transitive_len,
                config.min_distance_between_ranges,
                None,
                false,
                None,
                sequence_index,
                config.approximate_mode,
                None,
            )
        }
    } else {
        impg.query(seq_id, 0, seq_len as i32, false, None, sequence_index, config.approximate_mode)
    };

    for overlap in overlaps {
        let query_interval = &overlap.0;
        let target_interval = &overlap.2;

        let query_name = match impg.seq_index.get_name(query_interval.metadata) {
            Some(name) => name,
            None => continue,
        };
        let sample = extract_sample(query_name, separator);

        let is_reverse = query_interval.first > query_interval.last;
        let query_start = query_interval.first.min(query_interval.last) as i64;
        let query_end = query_interval.first.max(query_interval.last) as i64;
        let target_start = target_interval.first.min(target_interval.last) as i64;
        let target_end = target_interval.first.max(target_interval.last) as i64;

        alignments.push(SampleAlignmentInfo {
            sample,
            query_name: query_name.to_string(),
            query_start,
            query_end,
            target_start,
            target_end,
            is_reverse,
            cigar: Vec::new(),
        });
    }

    // 2. Query where anchor_seq is QUERY (this sample aligned to others)
    if let Some(reverse_entries) = reverse_index.get_all(anchor_seq) {
        for entry in reverse_entries {
            let target_sample = extract_sample(&entry.target_seq, separator);

            // Skip self-alignments
            if target_sample == anchor_sample {
                continue;
            }

            alignments.push(SampleAlignmentInfo {
                sample: target_sample,
                query_name: entry.target_seq.clone(),
                query_start: entry.target_start,
                query_end: entry.target_end,
                target_start: entry.query_start,
                target_end: entry.query_end,
                is_reverse: entry.is_reverse,
                cigar: Vec::new(),
            });
        }
    }

    alignments
}

/// Collect alignments for pool intervals using FromPool strategy
fn collect_alignments_from_pool(
    impg: &Impg,
    anchor_seq: &str,
    anchor_sample: &str,
    pool_intervals: &[(i64, i64)],
    reverse_index: &ReverseAlignmentIndex,
    config: &DepthConfig,
    sequence_index: Option<&UnifiedSequenceIndex>,
    separator: &str,
) -> Vec<SampleAlignmentInfo> {
    let mut alignments = Vec::new();
    let mut seen_alignments: FxHashSet<(String, i64, i64, i64, i64)> = FxHashSet::default();

    let seq_id = match impg.seq_index.get_id(anchor_seq) {
        Some(id) => id,
        None => return alignments,
    };

    // Merge nearby pool intervals to reduce query count
    let merged_intervals = merge_nearby_intervals(pool_intervals, 1000);

    for (start, end) in merged_intervals {
        // 1. Query where anchor_seq is TARGET
        let overlaps = if config.transitive {
            if config.transitive_dfs {
                impg.query_transitive_dfs(
                    seq_id,
                    start as i32,
                    end as i32,
                    None,
                    config.max_depth,
                    config.min_transitive_len,
                    config.min_distance_between_ranges,
                    None,
                    false,
                    None,
                    sequence_index,
                    config.approximate_mode,
                    None,
                )
            } else {
                impg.query_transitive_bfs(
                    seq_id,
                    start as i32,
                    end as i32,
                    None,
                    config.max_depth,
                    config.min_transitive_len,
                    config.min_distance_between_ranges,
                    None,
                    false,
                    None,
                    sequence_index,
                    config.approximate_mode,
                    None,
                )
            }
        } else {
            impg.query(seq_id, start as i32, end as i32, false, None, sequence_index, config.approximate_mode)
        };

        for overlap in overlaps {
            let query_interval = &overlap.0;
            let target_interval = &overlap.2;

            let query_name = match impg.seq_index.get_name(query_interval.metadata) {
                Some(name) => name,
                None => continue,
            };
            let sample = extract_sample(query_name, separator);

            let is_reverse = query_interval.first > query_interval.last;
            let query_start = query_interval.first.min(query_interval.last) as i64;
            let query_end = query_interval.first.max(query_interval.last) as i64;
            let target_start = target_interval.first.min(target_interval.last) as i64;
            let target_end = target_interval.first.max(target_interval.last) as i64;

            // Deduplicate
            let key = (query_name.to_string(), query_start, query_end, target_start, target_end);
            if seen_alignments.insert(key) {
                alignments.push(SampleAlignmentInfo {
                    sample,
                    query_name: query_name.to_string(),
                    query_start,
                    query_end,
                    target_start,
                    target_end,
                    is_reverse,
                    cigar: Vec::new(),
                });
            }
        }

        // 2. Query reverse alignments for this range
        let reverse_entries = reverse_index.query_range(anchor_seq, start, end);
        for entry in reverse_entries {
            let target_sample = extract_sample(&entry.target_seq, separator);

            if target_sample == anchor_sample {
                continue;
            }

            let key = (
                entry.target_seq.clone(),
                entry.target_start,
                entry.target_end,
                entry.query_start,
                entry.query_end,
            );
            if seen_alignments.insert(key) {
                alignments.push(SampleAlignmentInfo {
                    sample: target_sample,
                    query_name: entry.target_seq.clone(),
                    query_start: entry.target_start,
                    query_end: entry.target_end,
                    target_start: entry.query_start,
                    target_end: entry.query_end,
                    is_reverse: entry.is_reverse,
                    cigar: Vec::new(),
                });
            }
        }
    }

    alignments
}

/// Merge nearby intervals (gap <= max_gap)
fn merge_nearby_intervals(intervals: &[(i64, i64)], max_gap: i64) -> Vec<(i64, i64)> {
    if intervals.is_empty() {
        return Vec::new();
    }

    let mut sorted: Vec<_> = intervals.to_vec();
    sorted.sort_by_key(|&(s, _)| s);

    let mut merged = Vec::new();
    let mut current = sorted[0];

    for &(start, end) in &sorted[1..] {
        if start <= current.1 + max_gap {
            current.1 = current.1.max(end);
        } else {
            merged.push(current);
            current = (start, end);
        }
    }
    merged.push(current);

    merged
}

/// Process alignments for a single anchor sequence, generating depth windows
fn process_seq_to_windows(
    anchor_sample: &str,
    anchor_seq: &str,
    anchor_len: i64,
    alignments: &[SampleAlignmentInfo],
    pool_intervals: &[(i64, i64)],
    config: &DepthConfig,
    _separator: &str,
) -> SeqProcessResult {
    let mut windows = Vec::new();
    let mut regions_to_remove = Vec::new();

    if alignments.is_empty() && pool_intervals.is_empty() {
        return SeqProcessResult {
            windows,
            regions_to_remove,
        };
    }

    // Build events for sweep-line
    let mut events: Vec<DepthEvent> = Vec::new();
    let mut all_alignments: Vec<SampleAlignmentInfo> = Vec::new();

    // Add anchor sample itself
    all_alignments.push(SampleAlignmentInfo {
        sample: anchor_sample.to_string(),
        query_name: anchor_seq.to_string(),
        query_start: 0,
        query_end: anchor_len,
        target_start: 0,
        target_end: anchor_len,
        is_reverse: false,
        cigar: Vec::new(),
    });

    // Add other alignments
    all_alignments.extend(alignments.iter().cloned());

    // Create events
    for (idx, aln) in all_alignments.iter().enumerate() {
        events.push(DepthEvent {
            position: aln.target_start,
            is_start: true,
            haplotype: aln.sample.clone(),
            alignment_idx: idx,
        });
        events.push(DepthEvent {
            position: aln.target_end,
            is_start: false,
            haplotype: aln.sample.clone(),
            alignment_idx: idx,
        });
    }
    events.sort();

    // Sweep-line to compute raw windows
    let mut raw_windows: Vec<DepthWindow> = Vec::new();
    let mut active_alignments: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    let mut prev_pos: Option<i64> = None;
    let mut window_id = 0;

    for event in events {
        if let Some(prev) = prev_pos {
            if event.position > prev {
                let active_samples: Vec<&String> = active_alignments
                    .iter()
                    .filter(|(_, alns)| !alns.is_empty())
                    .map(|(s, _)| s)
                    .collect();

                if !active_samples.is_empty() {
                    let mut sample_positions: FxHashMap<String, (String, i64, i64)> =
                        FxHashMap::default();

                    for sample in &active_samples {
                        if let Some(aln_indices) = active_alignments.get(*sample) {
                            if let Some(&best_idx) = aln_indices.iter().max_by_key(|&&idx| {
                                let aln = &all_alignments[idx];
                                let overlap_start = prev.max(aln.target_start);
                                let overlap_end = event.position.min(aln.target_end);
                                overlap_end - overlap_start
                            }) {
                                let aln = &all_alignments[best_idx];
                                let (window_query_start, window_query_end) =
                                    map_target_to_query_linear(
                                        &aln.cigar,
                                        aln.target_start,
                                        aln.target_end,
                                        aln.query_start,
                                        aln.query_end,
                                        prev,
                                        event.position,
                                        aln.is_reverse,
                                    );

                                sample_positions.insert(
                                    (*sample).clone(),
                                    (aln.query_name.clone(), window_query_start, window_query_end),
                                );
                            }
                        }
                    }

                    if !sample_positions.is_empty() {
                        raw_windows.push(DepthWindow {
                            window_id,
                            depth: sample_positions.len(),
                            ref_name: anchor_seq.to_string(),
                            ref_start: prev,
                            ref_end: event.position,
                            haplotype_positions: sample_positions,
                        });
                        window_id += 1;
                    }
                }
            }
        }

        if event.is_start {
            active_alignments
                .entry(event.haplotype.clone())
                .or_default()
                .push(event.alignment_idx);
        } else if let Some(alns) = active_alignments.get_mut(&event.haplotype) {
            alns.retain(|&idx| idx != event.alignment_idx);
        }

        prev_pos = Some(event.position);
    }

    // Optionally merge adjacent windows
    if config.merge_adjacent {
        raw_windows = merge_adjacent_windows(raw_windows);
    }

    // Filter windows by pool intervals and collect regions to remove
    let pool_interval_set = {
        let mut is = IntervalSet::new();
        for &(start, end) in pool_intervals {
            // Simple merge for filtering
            is.intervals.push((start, end));
        }
        is.intervals.sort_by_key(|&(s, _)| s);
        is
    };

    for window in raw_windows {
        // Check if window overlaps with any pool interval
        let intersections = pool_interval_set.intersect(window.ref_start, window.ref_end);

        for (int_start, int_end) in intersections {
            // Create clipped window
            let mut clipped_window = window.clone();
            clipped_window.ref_start = int_start;
            clipped_window.ref_end = int_end;

            // Update sample positions proportionally
            let window_len = window.ref_end - window.ref_start;
            if window_len > 0 {
                let ratio_start =
                    (int_start - window.ref_start) as f64 / window_len as f64;
                let ratio_end = (int_end - window.ref_start) as f64 / window_len as f64;

                for (sample, (seq, orig_start, orig_end)) in &window.haplotype_positions {
                    let orig_len = orig_end - orig_start;
                    let new_start = orig_start + (orig_len as f64 * ratio_start) as i64;
                    let new_end = orig_start + (orig_len as f64 * ratio_end) as i64;

                    clipped_window
                        .haplotype_positions
                        .insert(sample.clone(), (seq.clone(), new_start, new_end));
                }
            }

            // Record regions to remove from pool (for all samples in window)
            for (sample, (seq, start, end)) in &clipped_window.haplotype_positions {
                regions_to_remove.push((sample.clone(), seq.clone(), *start, *end));
            }

            windows.push(clipped_window);
        }
    }

    SeqProcessResult {
        windows,
        regions_to_remove,
    }
}

/// Main function: compute depth using sample-based iteration with hybrid strategy
///
/// When `ref_sample` is specified (targeted mode), only computes depth for the specified
/// sample's sequences as reference. When `ref_sample` is None (global mode), computes
/// depth for all sequences.
pub fn compute_depth_by_sample(
    impg: &Impg,
    config: &DepthConfig,
    sequence_index: Option<&UnifiedSequenceIndex>,
    separator: &str,
    output_prefix: Option<&str>,
    fai_list: Option<&str>,
    ref_sample: Option<&str>,
) -> io::Result<()> {
    if let Some(ref_name) = ref_sample {
        info!("Computing depth coverage (targeted mode for sample: {})", ref_name);
    } else {
        info!("Computing depth coverage (global mode)");
    }

    // Phase 0: Preprocessing
    info!("Phase 0: Building indices...");

    // 0.1 Build sequence lengths from impg (for alignment processing)
    let seq_lengths = SequenceLengths::from_impg(impg, separator);
    info!(
        "  Found {} samples, {} sequences in alignments",
        seq_lengths.sample_to_seqs.len(),
        seq_lengths.lengths.len()
    );

    // 0.2 Build reverse alignment index
    info!("  Building reverse alignment index...");
    let reverse_index = ReverseAlignmentIndex::build(impg);
    info!(
        "  Reverse index contains {} alignments",
        reverse_index.total_alignments()
    );

    // Clear tree cache after building reverse index - the index now contains all needed data
    impg.clear_tree_cache();
    info!("  Tree cache cleared after building reverse index");

    // 0.3 Initialize pool - use FAI files if provided for full coverage
    let (pool_seq_lengths, mut pool) = if let Some(fai_path) = fai_list {
        info!("  Loading sequence lengths from FAI files...");
        let fai_seq_lengths = SequenceLengths::from_fai_list(fai_path, separator)?;
        info!(
            "  FAI: {} samples, {} sequences",
            fai_seq_lengths.sample_to_seqs.len(),
            fai_seq_lengths.lengths.len()
        );
        let p = Pool::new(&fai_seq_lengths);
        (fai_seq_lengths, p)
    } else {
        let p = Pool::new(&seq_lengths);
        (seq_lengths.clone(), p)
    };
    let initial_pool_length = pool.total_remaining_length();
    info!("  Pool initialized with {} bp total", initial_pool_length);

    // 0.4 Determine sample processing order (from FAI if provided, otherwise from impg)
    let sample_order = pool_seq_lengths.get_samples();
    info!("  Sample processing order: {:?}", sample_order);

    // 0.4.1 Validate ref_sample if specified (targeted mode)
    if let Some(ref_name) = ref_sample {
        if !sample_order.contains(&ref_name.to_string()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Reference sample '{}' not found. Available samples: {:?}",
                    ref_name, sample_order
                ),
            ));
        }
        info!("  Targeted mode: only processing sequences from sample '{}'", ref_name);
    }

    // 0.5 Prepare output
    let writer: Box<dyn Write> = if let Some(prefix) = output_prefix {
        let path = format!("{}.depth.tsv", prefix);
        Box::new(BufWriter::new(std::fs::File::create(&path)?))
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };
    let mut writer = BufWriter::new(writer);

    // Write header
    write!(writer, "window_id\tdepth")?;
    for sample in &sample_order {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    let mut global_window_id: usize = 0;
    let total_samples = sample_order.len();

    // Phase 1: Process samples
    info!("Phase 1: Processing {} samples...", total_samples);

    for (sample_idx, anchor_sample) in sample_order.iter().enumerate() {
        // In targeted mode, only process the specified reference sample
        if let Some(ref_name) = ref_sample {
            if anchor_sample != ref_name {
                continue;
            }
        }

        // Check if sample has remaining intervals
        if pool.is_sample_empty(anchor_sample) {
            debug!("Sample {} has no remaining intervals, skipping", anchor_sample);
            continue;
        }

        let is_first_round = sample_idx == 0;
        let strategy = choose_strategy(anchor_sample, &pool, impg, &seq_lengths, is_first_round);

        debug!(
            "Processing sample {}/{}: {} (strategy: {:?})",
            sample_idx + 1,
            total_samples,
            anchor_sample,
            strategy
        );

        // Get sequences for this sample
        let seqs = seq_lengths.get_seqs_of_sample(anchor_sample);

        // Process sequences in parallel
        let seq_results: Vec<SeqProcessResult> = seqs
            .par_iter()
            .filter_map(|seq_name| {
                let pool_intervals = pool.get_intervals(anchor_sample, seq_name);
                if pool_intervals.is_empty() {
                    return None;
                }

                let seq_len = seq_lengths.get_length(seq_name).unwrap_or(0);

                // Collect alignments based on strategy
                let alignments = match strategy {
                    QueryStrategy::FromIndex => collect_alignments_from_index(
                        impg,
                        seq_name,
                        anchor_sample,
                        &reverse_index,
                        config,
                        sequence_index,
                        separator,
                    ),
                    QueryStrategy::FromPool => collect_alignments_from_pool(
                        impg,
                        seq_name,
                        anchor_sample,
                        &pool_intervals,
                        &reverse_index,
                        config,
                        sequence_index,
                        separator,
                    ),
                };

                // Process to windows
                let result = process_seq_to_windows(
                    anchor_sample,
                    seq_name,
                    seq_len,
                    &alignments,
                    &pool_intervals,
                    config,
                    separator,
                );

                Some(result)
            })
            .collect();

        // Collect all windows from parallel processing
        let mut all_windows: Vec<DepthWindow> = Vec::new();
        for result in seq_results {
            all_windows.extend(result.windows);
        }

        // Track regions output in this batch to prevent intra-batch overlaps
        // Key: (sample, seq) -> IntervalSet of output regions
        let mut batch_output: FxHashMap<(String, String), IntervalSet> = FxHashMap::default();
        let mut regions_to_remove: Vec<(String, String, i64, i64)> = Vec::new();

        for mut window in all_windows {
            // Check each sample's region against:
            // 1. The global pool (for inter-sample filtering)
            // 2. The batch_output (for intra-batch overlap prevention)
            let mut valid_positions: FxHashMap<String, (String, i64, i64)> =
                FxHashMap::default();

            for (sample, (seq, start, end)) in &window.haplotype_positions {
                // First check global pool
                let pool_overlaps = pool.get_overlap(sample, seq, *start, *end);
                if pool_overlaps.is_empty() {
                    continue;
                }

                // Get the largest overlap with pool
                let (pool_start, pool_end) = pool_overlaps
                    .iter()
                    .max_by_key(|(s, e)| e - s)
                    .copied()
                    .unwrap();

                // Then check if this region overlaps with already-output regions in this batch
                let key = (sample.clone(), seq.clone());
                let batch_overlaps = batch_output
                    .get(&key)
                    .map(|is| is.intersect(pool_start, pool_end))
                    .unwrap_or_default();

                if batch_overlaps.is_empty() {
                    // No overlap with batch output, this region is valid
                    valid_positions.insert(sample.clone(), (seq.clone(), pool_start, pool_end));
                } else {
                    // There's overlap - try to find a non-overlapping portion
                    // Subtract batch_overlaps from [pool_start, pool_end]
                    let mut remaining = vec![(pool_start, pool_end)];
                    for (ov_start, ov_end) in batch_overlaps {
                        let mut new_remaining = Vec::new();
                        for (r_start, r_end) in remaining {
                            if r_end <= ov_start || r_start >= ov_end {
                                // No overlap
                                new_remaining.push((r_start, r_end));
                            } else {
                                // Split around overlap
                                if r_start < ov_start {
                                    new_remaining.push((r_start, ov_start));
                                }
                                if r_end > ov_end {
                                    new_remaining.push((ov_end, r_end));
                                }
                            }
                        }
                        remaining = new_remaining;
                    }

                    // Use the largest remaining interval if any
                    if let Some(&(best_start, best_end)) =
                        remaining.iter().max_by_key(|(s, e)| e - s)
                    {
                        if best_end > best_start {
                            valid_positions
                                .insert(sample.clone(), (seq.clone(), best_start, best_end));
                        }
                    }
                }
            }

            // Skip window if no valid samples remain
            if valid_positions.is_empty() {
                continue;
            }

            // Update window with only valid positions
            window.haplotype_positions = valid_positions;
            window.depth = window.haplotype_positions.len();
            window.window_id = global_window_id;
            global_window_id += 1;

            // Write window
            write!(writer, "{}\t{}", window.window_id, window.depth)?;
            for sample in &sample_order {
                if let Some((seq, start, end)) = window.haplotype_positions.get(sample) {
                    write!(writer, "\t{}:{}-{}", seq, start, end)?;
                } else {
                    write!(writer, "\tNA")?;
                }
            }
            writeln!(writer)?;

            // Record output regions for both batch tracking and pool update
            for (sample, (seq, start, end)) in &window.haplotype_positions {
                let key = (sample.clone(), seq.clone());
                batch_output
                    .entry(key)
                    .or_default()
                    .add(*start, *end);
                regions_to_remove.push((sample.clone(), seq.clone(), *start, *end));
            }
        }

        // Batch update pool at end
        pool.subtract_batch(&regions_to_remove);

        // Progress reporting
        if (sample_idx + 1) % 10 == 0 || sample_idx + 1 == total_samples {
            let remaining = pool.total_remaining_length();
            let processed_pct =
                100.0 * (1.0 - remaining as f64 / initial_pool_length as f64);
            info!(
                "  Progress: {}/{} samples, {:.1}% of genome covered",
                sample_idx + 1,
                total_samples,
                processed_pct
            );
        }

        // Early termination if pool is empty
        if pool.is_empty() {
            info!(
                "  All regions covered after {} samples, stopping early",
                sample_idx + 1
            );
            break;
        }
    }

    // Phase 2: Fill remaining pool intervals with depth=1 (uncovered regions)
    let remaining_length = pool.total_remaining_length();
    if remaining_length > 0 {
        let uncovered_pct = 100.0 * remaining_length as f64 / initial_pool_length as f64;
        info!(
            "Phase 2: Filling uncovered regions ({} bp, {:.2}%) with depth=1...",
            remaining_length, uncovered_pct
        );

        let mut gap_windows = 0u64;
        let mut gap_bp = 0i64;

        // Output remaining pool intervals as depth=1
        for sample in &sample_order {
            let seqs = pool_seq_lengths.get_seqs_of_sample(sample);
            for seq_name in seqs {
                let remaining_intervals = pool.get_intervals(sample, seq_name);

                for (gap_start, gap_end) in remaining_intervals {
                    if gap_end <= gap_start {
                        continue;
                    }

                    // Create depth=1 window with only this sample
                    write!(writer, "{}\t1", global_window_id)?;
                    for s in &sample_order {
                        if s == sample {
                            write!(writer, "\t{}:{}-{}", seq_name, gap_start, gap_end)?;
                        } else {
                            write!(writer, "\tNA")?;
                        }
                    }
                    writeln!(writer)?;

                    global_window_id += 1;
                    gap_windows += 1;
                    gap_bp += gap_end - gap_start;
                }
            }
        }

        info!(
            "  Generated {} gap windows covering {} bp",
            gap_windows, gap_bp
        );
    }

    writer.flush()?;
    info!("Output complete: {} windows generated", global_window_id);

    Ok(())
}

// ============================================================================
// Memory-efficient implementation using compressed bitmap coverage tracking
// ============================================================================

/// Resolution for bitmap coverage tracking (bp per bit)
/// Lower value = more accurate but more memory
/// 1bp: 3TB genome = 375GB (too much)
/// 10bp: 3TB genome = 37.5GB (acceptable for most systems)
/// For TB-scale data, recommend 10-100bp resolution
const COVERAGE_RESOLUTION: i64 = 1;

/// Compressed coverage tracker using BitVec
/// Memory usage: ~375MB for 3TB genome @ 100bp resolution
#[derive(Debug)]
pub struct CompressedCoverageTracker {
    /// (sample_idx, seq_idx) -> BitVec where each bit = COVERAGE_RESOLUTION bp
    coverage: FxHashMap<(u16, u32), BitVec>,
    /// Sequence lengths for bounds checking
    seq_lengths: Vec<i64>,
    /// Sample index mapping
    sample_to_idx: FxHashMap<String, u16>,
    /// Sequence index mapping
    seq_to_idx: FxHashMap<String, u32>,
    /// Reverse mappings for output (reserved for future use)
    #[allow(dead_code)]
    idx_to_sample: Vec<String>,
    #[allow(dead_code)]
    idx_to_seq: Vec<String>,
}

impl CompressedCoverageTracker {
    /// Create a new coverage tracker from sequence lengths
    pub fn new(seq_lengths_map: &SequenceLengths) -> Self {
        let mut sample_to_idx = FxHashMap::default();
        let mut seq_to_idx = FxHashMap::default();
        let mut idx_to_sample = Vec::new();
        let mut idx_to_seq = Vec::new();
        let mut seq_lengths = Vec::new();

        // Build sample index
        let mut samples: Vec<_> = seq_lengths_map.sample_to_seqs.keys().collect();
        samples.sort();
        for (i, sample) in samples.iter().enumerate() {
            sample_to_idx.insert((*sample).clone(), i as u16);
            idx_to_sample.push((*sample).clone());
        }

        // Build sequence index
        let mut seqs: Vec<_> = seq_lengths_map.lengths.keys().collect();
        seqs.sort();
        for (i, seq) in seqs.iter().enumerate() {
            seq_to_idx.insert((*seq).clone(), i as u32);
            idx_to_seq.push((*seq).clone());
            seq_lengths.push(*seq_lengths_map.lengths.get(*seq).unwrap_or(&0));
        }

        // Initialize coverage bitmaps (all zeros = uncovered)
        let mut coverage = FxHashMap::default();
        for (seq_name, &length) in &seq_lengths_map.lengths {
            let sample = seq_lengths_map.seq_to_sample.get(seq_name).unwrap();
            if let (Some(&sample_idx), Some(&seq_idx)) = (
                sample_to_idx.get(sample),
                seq_to_idx.get(seq_name.as_str()),
            ) {
                let num_bits = ((length + COVERAGE_RESOLUTION - 1) / COVERAGE_RESOLUTION) as usize;
                coverage.insert((sample_idx, seq_idx), bitvec![0; num_bits]);
            }
        }

        CompressedCoverageTracker {
            coverage,
            seq_lengths,
            sample_to_idx,
            seq_to_idx,
            idx_to_sample,
            idx_to_seq,
        }
    }

    /// Convert position to bit index
    #[inline]
    fn pos_to_bit(pos: i64) -> usize {
        (pos / COVERAGE_RESOLUTION) as usize
    }

    /// Convert bit index to position range
    #[inline]
    fn bit_to_pos_range(bit_idx: usize) -> (i64, i64) {
        let start = bit_idx as i64 * COVERAGE_RESOLUTION;
        (start, start + COVERAGE_RESOLUTION)
    }

    /// Check if any portion of [start, end) is uncovered
    pub fn has_uncovered(&self, sample: &str, seq: &str, start: i64, end: i64) -> bool {
        let sample_idx = match self.sample_to_idx.get(sample) {
            Some(&idx) => idx,
            None => return false,
        };
        let seq_idx = match self.seq_to_idx.get(seq) {
            Some(&idx) => idx,
            None => return false,
        };

        if let Some(bits) = self.coverage.get(&(sample_idx, seq_idx)) {
            let start_bit = Self::pos_to_bit(start);
            let end_bit = Self::pos_to_bit(end.saturating_sub(1)).min(bits.len().saturating_sub(1));

            for bit_idx in start_bit..=end_bit {
                if bit_idx < bits.len() && !bits[bit_idx] {
                    return true;
                }
            }
        }
        false
    }

    /// Get uncovered intervals within [start, end)
    /// Returns merged intervals of uncovered regions
    pub fn get_uncovered_intervals(
        &self,
        sample: &str,
        seq: &str,
        start: i64,
        end: i64,
    ) -> Vec<(i64, i64)> {
        let sample_idx = match self.sample_to_idx.get(sample) {
            Some(&idx) => idx,
            None => return vec![(start, end)],
        };
        let seq_idx = match self.seq_to_idx.get(seq) {
            Some(&idx) => idx,
            None => return vec![(start, end)],
        };

        let bits = match self.coverage.get(&(sample_idx, seq_idx)) {
            Some(b) => b,
            None => return vec![(start, end)],
        };

        let seq_len = self.seq_lengths.get(seq_idx as usize).copied().unwrap_or(end);
        let clamped_end = end.min(seq_len);

        let start_bit = Self::pos_to_bit(start);
        let end_bit = Self::pos_to_bit(clamped_end.saturating_sub(1));

        let mut result = Vec::new();
        let mut current_start: Option<i64> = None;

        for bit_idx in start_bit..=end_bit {
            if bit_idx >= bits.len() {
                break;
            }

            let (bit_start, bit_end) = Self::bit_to_pos_range(bit_idx);
            let interval_start = bit_start.max(start);
            let _interval_end = bit_end.min(clamped_end);

            if !bits[bit_idx] {
                // Uncovered
                if current_start.is_none() {
                    current_start = Some(interval_start);
                }
            } else {
                // Covered - close current interval if open
                if let Some(s) = current_start {
                    result.push((s, interval_start));
                    current_start = None;
                }
            }
        }

        // Close final interval
        if let Some(s) = current_start {
            result.push((s, clamped_end));
        }

        result
    }

    /// Mark region [start, end) as covered
    pub fn mark_covered(&mut self, sample: &str, seq: &str, start: i64, end: i64) {
        let sample_idx = match self.sample_to_idx.get(sample) {
            Some(&idx) => idx,
            None => return,
        };
        let seq_idx = match self.seq_to_idx.get(seq) {
            Some(&idx) => idx,
            None => return,
        };

        if let Some(bits) = self.coverage.get_mut(&(sample_idx, seq_idx)) {
            let start_bit = Self::pos_to_bit(start);
            let end_bit = Self::pos_to_bit(end.saturating_sub(1));

            for bit_idx in start_bit..=end_bit {
                if bit_idx < bits.len() {
                    bits.set(bit_idx, true);
                }
            }
        }
    }

    /// Get all uncovered intervals for a (sample, seq) pair
    pub fn get_all_uncovered(&self, sample: &str, seq: &str) -> Vec<(i64, i64)> {
        let sample_idx = match self.sample_to_idx.get(sample) {
            Some(&idx) => idx,
            None => return Vec::new(),
        };
        let seq_idx = match self.seq_to_idx.get(seq) {
            Some(&idx) => idx,
            None => return Vec::new(),
        };

        let bits = match self.coverage.get(&(sample_idx, seq_idx)) {
            Some(b) => b,
            None => return Vec::new(),
        };

        let seq_len = self.seq_lengths.get(seq_idx as usize).copied().unwrap_or(0);

        let mut result = Vec::new();
        let mut current_start: Option<i64> = None;

        for (bit_idx, bit) in bits.iter().enumerate() {
            let (bit_start, bit_end) = Self::bit_to_pos_range(bit_idx);
            let _interval_end = bit_end.min(seq_len);

            if !*bit {
                if current_start.is_none() {
                    current_start = Some(bit_start);
                }
            } else if let Some(s) = current_start {
                result.push((s, bit_start));
                current_start = None;
            }
        }

        if let Some(s) = current_start {
            result.push((s, seq_len));
        }

        result
    }

    /// Get total uncovered length
    pub fn total_uncovered_length(&self) -> i64 {
        let mut total = 0i64;

        for ((_sample_idx, seq_idx), bits) in &self.coverage {
            let seq_len = self.seq_lengths.get(*seq_idx as usize).copied().unwrap_or(0);

            for (bit_idx, bit) in bits.iter().enumerate() {
                if !*bit {
                    let (bit_start, bit_end) = Self::bit_to_pos_range(bit_idx);
                    let interval_len = bit_end.min(seq_len) - bit_start;
                    if interval_len > 0 {
                        total += interval_len;
                    }
                }
            }
        }

        total
    }

    /// Get memory usage estimate in bytes
    pub fn memory_usage(&self) -> usize {
        let mut total = 0;
        for (_, bits) in &self.coverage {
            total += bits.len() / 8 + 1;
        }
        total
    }
}

/// Memory-efficient depth computation using compressed bitmap
///
/// Key improvements over compute_depth_by_sample:
/// 1. Uses BitVec instead of IntervalSet (memory: ~50GB -> ~500MB)
/// 2. No pre-built reverse index (memory: ~30GB -> 0)
/// 3. Streaming output (no batch accumulation)
///
/// When `ref_sample` is specified (targeted mode), only computes depth for the specified
/// sample's sequences as reference. When `ref_sample` is None (global mode), computes
/// depth for all sequences.
pub fn compute_depth_by_sample_v2(
    impg: &Impg,
    config: &DepthConfig,
    sequence_index: Option<&UnifiedSequenceIndex>,
    separator: &str,
    output_prefix: Option<&str>,
    fai_list: Option<&str>,
    ref_sample: Option<&str>,
) -> io::Result<()> {
    if let Some(ref_name) = ref_sample {
        info!("Computing depth coverage (memory-efficient mode v2, targeted for sample: {})", ref_name);
    } else {
        info!("Computing depth coverage (memory-efficient mode v2, global mode)");
    }

    // Phase 0: Build sequence lengths
    info!("Phase 0: Building indices...");

    let seq_lengths = SequenceLengths::from_impg(impg, separator);
    info!(
        "  Found {} samples, {} sequences in alignments",
        seq_lengths.sample_to_seqs.len(),
        seq_lengths.lengths.len()
    );

    // Use FAI files for pool if provided
    let pool_seq_lengths = if let Some(fai_path) = fai_list {
        info!("  Loading sequence lengths from FAI files...");
        let fai_seq_lengths = SequenceLengths::from_fai_list(fai_path, separator)?;
        info!(
            "  FAI: {} samples, {} sequences",
            fai_seq_lengths.sample_to_seqs.len(),
            fai_seq_lengths.lengths.len()
        );
        fai_seq_lengths
    } else {
        seq_lengths.clone()
    };

    // Initialize compressed coverage tracker
    let mut coverage = CompressedCoverageTracker::new(&pool_seq_lengths);
    info!(
        "  Coverage tracker initialized: {} bytes",
        coverage.memory_usage()
    );

    let initial_uncovered = coverage.total_uncovered_length();
    info!("  Total genome size: {} bp", initial_uncovered);

    // Sample processing order
    let sample_order = pool_seq_lengths.get_samples();
    info!("  Sample order: {:?}", sample_order);

    // Validate ref_sample if specified (targeted mode)
    if let Some(ref_name) = ref_sample {
        if !sample_order.contains(&ref_name.to_string()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Reference sample '{}' not found. Available samples: {:?}",
                    ref_name, sample_order
                ),
            ));
        }
        info!("  Targeted mode: only processing sequences from sample '{}'", ref_name);
    }

    // Build lightweight reverse index only for global mode
    // For --ref mode, skip this expensive step to save memory
    let query_to_targets = if ref_sample.is_none() {
        info!("  Building lightweight reverse index (global mode)...");
        let map = impg.build_query_to_targets_map();
        info!(
            "  Reverse index built: {} query sequences mapped",
            map.len()
        );
        Some(map)
    } else {
        info!("  Skipping reverse index (targeted mode - using forward queries only)");
        None
    };

    // Prepare output
    let writer: Box<dyn Write> = if let Some(prefix) = output_prefix {
        let path = format!("{}.depth.tsv", prefix);
        Box::new(BufWriter::new(std::fs::File::create(&path)?))
    } else {
        Box::new(BufWriter::new(std::io::stdout()))
    };
    let mut writer = BufWriter::new(writer);

    // Write header
    write!(writer, "window_id\tdepth")?;
    for sample in &sample_order {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    let mut global_window_id: usize = 0;
    let total_samples = sample_order.len();

    // Phase 1: Process each sample as anchor
    info!("Phase 1: Processing {} samples...", total_samples);

    for (sample_idx, anchor_sample) in sample_order.iter().enumerate() {
        // In targeted mode, only process the specified reference sample
        if let Some(ref_name) = ref_sample {
            if anchor_sample != ref_name {
                continue;
            }
        }

        let seqs = seq_lengths.get_seqs_of_sample(anchor_sample);

        debug!(
            "Processing sample {}/{}: {} ({} sequences)",
            sample_idx + 1,
            total_samples,
            anchor_sample,
            seqs.len()
        );

        // Process each sequence of this sample
        for anchor_seq in seqs {
            let seq_id = match impg.seq_index.get_id(anchor_seq) {
                Some(id) => id,
                None => continue,
            };
            let seq_len = impg.seq_index.get_len_from_id(seq_id).unwrap_or(0) as i64;

            // Check if there are any uncovered regions for this sequence
            let uncovered = coverage.get_all_uncovered(anchor_sample, anchor_seq);
            if uncovered.is_empty() {
                continue;
            }

            // Collect forward alignments (anchor_seq as target)
            let overlaps = if config.transitive {
                if config.transitive_dfs {
                    impg.query_transitive_dfs(
                        seq_id, 0, seq_len as i32, None,
                        config.max_depth, config.min_transitive_len,
                        config.min_distance_between_ranges, None, false, None,
                        sequence_index, config.approximate_mode, None,
                    )
                } else {
                    impg.query_transitive_bfs(
                        seq_id, 0, seq_len as i32, None,
                        config.max_depth, config.min_transitive_len,
                        config.min_distance_between_ranges, None, false, None,
                        sequence_index, config.approximate_mode, None,
                    )
                }
            } else {
                impg.query(seq_id, 0, seq_len as i32, false, None, sequence_index, config.approximate_mode)
            };

            // Collect reverse alignments (anchor_seq as query) using pre-built index
            // Skip in targeted mode (--ref) to save memory - forward queries capture most coverage
            let reverse_alignments = if let Some(ref map) = query_to_targets {
                impg.query_reverse_for_depth_with_map(seq_id, map)
            } else {
                Vec::new() // Skip reverse queries in targeted mode
            };

            // Build alignment info list
            let mut all_alignments: Vec<SampleAlignmentInfo> = Vec::new();

            // Add self
            all_alignments.push(SampleAlignmentInfo {
                sample: anchor_sample.clone(),
                query_name: anchor_seq.clone(),
                query_start: 0,
                query_end: seq_len,
                target_start: 0,
                target_end: seq_len,
                is_reverse: false,
                cigar: Vec::new(),
            });

            // Forward alignments
            for overlap in &overlaps {
                let query_interval = &overlap.0;
                let target_interval = &overlap.2;

                let query_name = match impg.seq_index.get_name(query_interval.metadata) {
                    Some(name) => name,
                    None => continue,
                };
                let sample = extract_sample(query_name, separator);

                let is_reverse = query_interval.first > query_interval.last;
                let query_start = query_interval.first.min(query_interval.last) as i64;
                let query_end = query_interval.first.max(query_interval.last) as i64;
                let target_start = target_interval.first.min(target_interval.last) as i64;
                let target_end = target_interval.first.max(target_interval.last) as i64;

                all_alignments.push(SampleAlignmentInfo {
                    sample,
                    query_name: query_name.to_string(),
                    query_start,
                    query_end,
                    target_start,
                    target_end,
                    is_reverse,
                    cigar: Vec::new(),
                });
            }

            // Reverse alignments
            for &(ref_start, ref_end, target_id) in &reverse_alignments {
                if let Some(target_name) = impg.seq_index.get_name(target_id) {
                    let sample = extract_sample(target_name, separator);
                    if sample == *anchor_sample {
                        continue;
                    }

                    all_alignments.push(SampleAlignmentInfo {
                        sample,
                        query_name: target_name.to_string(),
                        query_start: ref_start as i64,
                        query_end: ref_end as i64,
                        target_start: ref_start as i64,
                        target_end: ref_end as i64,
                        is_reverse: false,
                        cigar: Vec::new(),
                    });
                }
            }

            // Build sweep-line events
            let mut events: Vec<DepthEvent> = Vec::new();
            for (idx, aln) in all_alignments.iter().enumerate() {
                events.push(DepthEvent {
                    position: aln.target_start,
                    is_start: true,
                    haplotype: aln.sample.clone(),
                    alignment_idx: idx,
                });
                events.push(DepthEvent {
                    position: aln.target_end,
                    is_start: false,
                    haplotype: aln.sample.clone(),
                    alignment_idx: idx,
                });
            }
            events.sort();

            // Sweep-line to generate windows and output immediately
            let mut active_alignments: FxHashMap<String, Vec<usize>> = FxHashMap::default();
            let mut prev_pos: Option<i64> = None;

            for event in events {
                if let Some(prev) = prev_pos {
                    if event.position > prev {
                        // Check if this window overlaps any uncovered region
                        let window_start = prev;
                        let window_end = event.position;

                        // Get active samples
                        let active_samples: Vec<&String> = active_alignments
                            .iter()
                            .filter(|(_, alns)| !alns.is_empty())
                            .map(|(s, _)| s)
                            .collect();

                        if !active_samples.is_empty() {
                            // Build sample positions for this window
                            let mut sample_positions: FxHashMap<String, (String, i64, i64)> =
                                FxHashMap::default();

                            for sample in &active_samples {
                                if let Some(aln_indices) = active_alignments.get(*sample) {
                                    if let Some(&best_idx) = aln_indices.iter().max_by_key(|&&idx| {
                                        let aln = &all_alignments[idx];
                                        let overlap_start = prev.max(aln.target_start);
                                        let overlap_end = event.position.min(aln.target_end);
                                        overlap_end - overlap_start
                                    }) {
                                        let aln = &all_alignments[best_idx];
                                        let (q_start, q_end) = map_target_to_query_linear(
                                            &aln.cigar,
                                            aln.target_start, aln.target_end,
                                            aln.query_start, aln.query_end,
                                            window_start, window_end,
                                            aln.is_reverse,
                                        );
                                        sample_positions.insert(
                                            (*sample).clone(),
                                            (aln.query_name.clone(), q_start, q_end),
                                        );
                                    }
                                }
                            }

                            // Filter: only output samples with uncovered regions
                            let mut valid_positions: FxHashMap<String, (String, i64, i64)> =
                                FxHashMap::default();

                            for (sample, (seq, start, end)) in &sample_positions {
                                // Check if this sample has uncovered region
                                let uncovered_parts =
                                    coverage.get_uncovered_intervals(sample, seq, *start, *end);

                                if !uncovered_parts.is_empty() {
                                    // Use the largest uncovered portion
                                    if let Some(&(unc_start, unc_end)) = uncovered_parts
                                        .iter()
                                        .max_by_key(|(s, e)| e - s)
                                    {
                                        valid_positions.insert(
                                            sample.clone(),
                                            (seq.clone(), unc_start, unc_end),
                                        );
                                    }
                                }
                            }

                            // Output if any valid samples
                            if !valid_positions.is_empty() {
                                let depth = valid_positions.len();

                                write!(writer, "{}\t{}", global_window_id, depth)?;
                                for sample in &sample_order {
                                    if let Some((seq, start, end)) = valid_positions.get(sample) {
                                        write!(writer, "\t{}:{}-{}", seq, start, end)?;
                                    } else {
                                        write!(writer, "\tNA")?;
                                    }
                                }
                                writeln!(writer)?;
                                global_window_id += 1;

                                // Mark as covered
                                for (sample, (seq, start, end)) in &valid_positions {
                                    coverage.mark_covered(sample, seq, *start, *end);
                                }
                            }
                        }
                    }
                }

                // Update active alignments
                if event.is_start {
                    active_alignments
                        .entry(event.haplotype.clone())
                        .or_default()
                        .push(event.alignment_idx);
                } else if let Some(alns) = active_alignments.get_mut(&event.haplotype) {
                    alns.retain(|&idx| idx != event.alignment_idx);
                }

                prev_pos = Some(event.position);
            }
        }

        // Progress reporting (every 10 samples or at the end)
        if (sample_idx + 1) % 10 == 0 || sample_idx + 1 == total_samples {
            let remaining = coverage.total_uncovered_length();
            let covered_pct = 100.0 * (1.0 - remaining as f64 / initial_uncovered as f64);
            info!(
                "  Progress: {}/{} samples, {:.1}% covered, {} windows, {} cached trees",
                sample_idx + 1,
                total_samples,
                covered_pct,
                global_window_id,
                impg.cached_tree_count()
            );
        }

        // Clear tree cache after processing each sample to control memory usage
        // Trees will be reloaded from disk as needed for subsequent samples
        impg.clear_tree_cache();
    }

    // Phase 2: Output uncovered regions as depth=1
    let final_uncovered = coverage.total_uncovered_length();
    if final_uncovered > 0 {
        let uncovered_pct = 100.0 * final_uncovered as f64 / initial_uncovered as f64;
        info!(
            "Phase 2: Outputting uncovered regions ({} bp, {:.2}%) as depth=1...",
            final_uncovered, uncovered_pct
        );

        let mut gap_count = 0u64;

        for sample in &sample_order {
            let seqs = pool_seq_lengths.get_seqs_of_sample(sample);
            for seq_name in seqs {
                let uncovered = coverage.get_all_uncovered(sample, seq_name);
                for (gap_start, gap_end) in uncovered {
                    if gap_end <= gap_start {
                        continue;
                    }

                    write!(writer, "{}\t1", global_window_id)?;
                    for s in &sample_order {
                        if s == sample {
                            write!(writer, "\t{}:{}-{}", seq_name, gap_start, gap_end)?;
                        } else {
                            write!(writer, "\tNA")?;
                        }
                    }
                    writeln!(writer)?;

                    global_window_id += 1;
                    gap_count += 1;
                }
            }
        }

        info!("  Generated {} gap windows", gap_count);
    }

    writer.flush()?;
    info!(
        "Output complete: {} windows, memory tracker: {} bytes",
        global_window_id,
        coverage.memory_usage()
    );

    Ok(())
}

// ============================================================================
// Streaming depth computation - memory efficient for large datasets
// ============================================================================

/// Coverage entry storing both reference and query coordinates
#[derive(Debug, Clone)]
struct CoverageEntry {
    sample: String,
    query_seq: String,
    ref_start: i64,
    ref_end: i64,
    query_start: i64,
    query_end: i64,
    is_reverse: bool, // true if query is on reverse strand relative to ref
}

/// Streaming depth computation that iterates through index without loading all trees.
/// Memory usage is O(ref_genome_size × num_samples) instead of O(total_alignments).
///
/// Key insight: We load one tree at a time, extract coverage info, then release it.
/// Uses parallel processing for scanning trees.
pub fn compute_depth_streaming(
    impg: &Impg,
    ref_sample: &str,
    separator: &str,
    output_prefix: Option<&str>,
    merge_adjacent: bool,
    fai_list: Option<&str>,
) -> io::Result<()> {
    info!(
        "Computing depth (streaming mode) for reference sample: {}",
        ref_sample
    );

    // Phase 1: Identify reference sequences and collect all samples
    info!("Phase 1: Identifying reference sequences...");

    let seq_lengths = SequenceLengths::from_impg(impg, separator);
    let all_samples = seq_lengths.get_samples();

    // Get reference sequence IDs
    let ref_seq_ids: Vec<(u32, String, i64)> = (0..impg.seq_index.len() as u32)
        .filter_map(|id| {
            let name = impg.seq_index.get_name(id)?;
            let sample = extract_sample(name, separator);
            if sample == ref_sample {
                let len = impg.seq_index.get_len_from_id(id)? as i64;
                Some((id, name.to_string(), len))
            } else {
                None
            }
        })
        .collect();

    if ref_seq_ids.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("No sequences found for reference sample '{}'", ref_sample),
        ));
    }

    // If FAI list provided, use it to get complete sequence lengths
    let ref_seq_lengths: FxHashMap<String, i64> = if let Some(fai_path) = fai_list {
        info!("  Loading sequence lengths from FAI files...");
        let fai_lengths = load_fai_lengths_for_sample(fai_path, ref_sample, separator)?;
        info!("  Loaded {} sequences from FAI", fai_lengths.len());
        fai_lengths
    } else {
        ref_seq_ids.iter().map(|(_, name, len)| (name.clone(), *len)).collect()
    };

    info!(
        "  Found {} reference sequences, {} total samples",
        ref_seq_ids.len(),
        all_samples.len()
    );

    // Phase 2: Parallel scan of index to accumulate coverage
    info!("Phase 2: Scanning index for coverage (parallel)...");

    // Create a set of reference sequence IDs for quick lookup
    let ref_id_set: FxHashSet<u32> = ref_seq_ids.iter().map(|(id, _, _)| *id).collect();
    let ref_id_to_name: FxHashMap<u32, String> = ref_seq_ids.iter().map(|(id, name, _)| (*id, name.clone())).collect();

    // Get all target IDs from forest map
    let all_target_ids: Vec<u32> = impg.get_all_target_ids();
    let total_trees = all_target_ids.len();
    info!(
        "  Processing {} trees in parallel with {} threads...",
        total_trees,
        rayon::current_num_threads()
    );

    // Progress tracking for Phase 2
    let processed_trees = AtomicUsize::new(0);

    // Parallel processing: each target_id is processed independently
    let parallel_results: Vec<(FxHashMap<String, Vec<CoverageEntry>>, u64)> = all_target_ids
        .par_iter()
        .map(|&target_id| {
            let mut local_coverage: FxHashMap<String, Vec<CoverageEntry>> = FxHashMap::default();
            let mut local_alignments = 0u64;

            // Progress logging every 1000 trees
            let count = processed_trees.fetch_add(1, Ordering::Relaxed);
            if count % 1000 == 0 || count < 10 {
                debug!(
                    "    Progress: {}/{} trees (thread {:?})",
                    count,
                    total_trees,
                    std::thread::current().id()
                );
            }

            {
                let target_name = match impg.seq_index.get_name(target_id) {
                    Some(n) => n.to_string(),
                    None => return (local_coverage, local_alignments),
                };
                let target_sample = extract_sample(&target_name, separator);
                let target_is_ref = ref_id_set.contains(&target_id);

                // Load tree without caching (avoids write lock contention in parallel)
                let tree = match impg.load_tree_no_cache(target_id) {
                    Some(t) => t,
                    None => return (local_coverage, local_alignments),
                };

                // Iterate through all alignments in this tree
                for interval in tree.iter() {
                    let query_id = interval.metadata.query_id();
                    let query_name = match impg.seq_index.get_name(query_id) {
                        Some(n) => n.to_string(),
                        None => continue,
                    };
                    let query_sample = extract_sample(&query_name, separator);

                    let target_start = interval.first as i64;
                    let target_end = interval.last as i64;
                    let is_reverse = interval.metadata.query_start() > interval.metadata.query_end();
                    let query_start = interval.metadata.query_start().min(interval.metadata.query_end()) as i64;
                    let query_end = interval.metadata.query_start().max(interval.metadata.query_end()) as i64;

                    // Case 1: target is reference -> query sample covers target region
                    if target_is_ref {
                        local_coverage
                            .entry(target_name.clone())
                            .or_default()
                            .push(CoverageEntry {
                                sample: query_sample.clone(),
                                query_seq: query_name.clone(),
                                ref_start: target_start,
                                ref_end: target_end,
                                query_start,
                                query_end,
                                is_reverse,
                            });
                        local_alignments += 1;
                    }

                    // Case 2: query is reference -> target sample covers query region
                    if ref_id_set.contains(&query_id) {
                        if let Some(ref_name) = ref_id_to_name.get(&query_id) {
                            local_coverage
                                .entry(ref_name.clone())
                                .or_default()
                                .push(CoverageEntry {
                                    sample: target_sample.clone(),
                                    query_seq: target_name.clone(),
                                    ref_start: query_start,
                                    ref_end: query_end,
                                    query_start: target_start,
                                    query_end: target_end,
                                    is_reverse,
                                });
                            local_alignments += 1;
                        }
                    }
                }
            }

            (local_coverage, local_alignments)
        })
        .collect();

    // Note: No need to clear tree cache since load_tree_no_cache doesn't cache

    // Merge results from all threads
    let mut coverage: FxHashMap<String, Vec<CoverageEntry>> = FxHashMap::default();
    for (_, name, _) in &ref_seq_ids {
        coverage.insert(name.clone(), Vec::new());
    }

    let mut alignments_found = 0u64;
    for (local_coverage, local_alignments) in parallel_results {
        alignments_found += local_alignments;
        for (ref_seq, entries) in local_coverage {
            coverage
                .entry(ref_seq)
                .or_default()
                .extend(entries);
        }
    }

    info!(
        "  Processed {} trees, found {} relevant alignments",
        total_trees, alignments_found
    );

    // Add self-coverage for reference sample
    for (ref_seq, seq_len) in &ref_seq_lengths {
        if let Some(entries) = coverage.get_mut(ref_seq) {
            entries.push(CoverageEntry {
                sample: ref_sample.to_string(),
                query_seq: ref_seq.clone(),
                ref_start: 0,
                ref_end: *seq_len,
                query_start: 0,
                query_end: *seq_len,
                is_reverse: false, // self-coverage is always forward
            });
        }
    }

    // Phase 3: Generate depth output using sweep-line (parallel with chunking)
    info!("Phase 3: Generating depth output (parallel)...");

    // Sort samples for consistent output, with reference sample first
    let mut sample_order: Vec<String> = all_samples;
    sample_order.sort_by(|a, b| {
        // Reference sample comes first
        if a == ref_sample {
            std::cmp::Ordering::Less
        } else if b == ref_sample {
            std::cmp::Ordering::Greater
        } else {
            a.cmp(b)
        }
    });

    // Sort reference sequences for consistent output
    let mut ref_seq_list: Vec<(String, i64)> = ref_seq_lengths.into_iter().collect();
    ref_seq_list.sort_by(|a, b| a.0.cmp(&b.0));

    // Chunk size for parallel processing (10MB)
    const CHUNK_SIZE: i64 = 10_000_000;

    // Build list of all chunks across all reference sequences
    // Each chunk is (ref_seq_idx, chunk_idx, chunk_start, chunk_end)
    let mut all_chunks: Vec<(usize, usize, i64, i64)> = Vec::new();
    for (ref_idx, (_, seq_len)) in ref_seq_list.iter().enumerate() {
        let num_chunks = ((*seq_len + CHUNK_SIZE - 1) / CHUNK_SIZE) as usize;
        for chunk_idx in 0..num_chunks {
            let chunk_start = (chunk_idx as i64) * CHUNK_SIZE;
            let chunk_end = ((chunk_idx as i64 + 1) * CHUNK_SIZE).min(*seq_len);
            all_chunks.push((ref_idx, chunk_idx, chunk_start, chunk_end));
        }
    }

    info!(
        "  Processing {} chunks across {} sequences with {} threads",
        all_chunks.len(),
        ref_seq_list.len(),
        rayon::current_num_threads()
    );

    // Create temp directory for parallel output files
    let temp_dir = tempfile::tempdir()?;
    let temp_dir_path = temp_dir.path().to_path_buf();

    // Parallel processing: each chunk writes to its own temp file
    // Returns (ref_seq_idx, chunk_idx, temp_file_path, window_count)
    let processed_chunks = AtomicUsize::new(0);
    let total_chunks = all_chunks.len();

    let chunk_results: Vec<(usize, usize, std::path::PathBuf, usize)> = all_chunks
        .par_iter()
        .map(|&(ref_idx, chunk_idx, chunk_start, chunk_end)| {
            let temp_path = temp_dir_path.join(format!("depth_{:06}_{:06}.tmp", ref_idx, chunk_idx));
            let mut window_count = 0usize;

            let (ref_seq, seq_len) = &ref_seq_list[ref_idx];

            // Progress logging every 100 chunks
            let count = processed_chunks.fetch_add(1, Ordering::Relaxed);
            if count % 100 == 0 || count < 10 {
                debug!(
                    "    Progress: {}/{} chunks (thread {:?})",
                    count,
                    total_chunks,
                    std::thread::current().id()
                );
            }
            let entries = match coverage.get(ref_seq) {
                Some(e) => e,
                None => return (ref_idx, chunk_idx, temp_path, 0),
            };

            if entries.is_empty() {
                return (ref_idx, chunk_idx, temp_path, 0);
            }

            // Find entries that overlap this chunk
            let chunk_entries: Vec<(usize, &CoverageEntry)> = entries
                .iter()
                .enumerate()
                .filter(|(_, e)| e.ref_start < chunk_end && e.ref_end > chunk_start)
                .collect();

            if chunk_entries.is_empty() {
                // No coverage in this chunk - output self-coverage only if this is part of ref
                let file = match std::fs::File::create(&temp_path) {
                    Ok(f) => f,
                    Err(_) => return (ref_idx, chunk_idx, temp_path, 0),
                };
                let mut temp_writer = std::io::BufWriter::new(file);

                // Output uncovered region with depth=1 (self only)
                let _ = write!(temp_writer, "1");
                for s in &sample_order {
                    if s == ref_sample {
                        let _ = write!(temp_writer, "\t{}:{}-{}", ref_seq, chunk_start, chunk_end);
                    } else {
                        let _ = write!(temp_writer, "\tNA");
                    }
                }
                let _ = writeln!(temp_writer);
                let _ = temp_writer.flush();
                return (ref_idx, chunk_idx, temp_path, 1);
            }

            // Create temp file writer
            let file = match std::fs::File::create(&temp_path) {
                Ok(f) => f,
                Err(_) => return (ref_idx, chunk_idx, temp_path, 0),
            };
            let mut temp_writer = std::io::BufWriter::new(file);

            // Build sweep-line events for this chunk
            // Include events at chunk boundaries for entries that span the boundary
            let mut events: Vec<(i64, bool, usize)> = Vec::new();
            for &(orig_idx, entry) in &chunk_entries {
                // Event for entry start (or chunk start if entry started earlier)
                let effective_start = entry.ref_start.max(chunk_start);
                // Event for entry end (or chunk end if entry extends beyond)
                let effective_end = entry.ref_end.min(chunk_end);

                if effective_start < effective_end {
                    events.push((effective_start, true, orig_idx));
                    events.push((effective_end, false, orig_idx));
                }
            }

            // Sort events: by position, then ends before starts at same position
            events.sort_by(|a, b| {
                a.0.cmp(&b.0)
                    .then_with(|| a.1.cmp(&b.1))
            });

            // Sweep-line processing
            let mut active_entries: FxHashMap<String, usize> = FxHashMap::default();
            let mut prev_pos: Option<i64> = None;
            let mut prev_active: Option<FxHashMap<String, usize>> = None;

            for (pos, is_start, entry_idx) in events {
                if let Some(prev) = prev_pos {
                    if pos > prev && !active_entries.is_empty() {
                        let should_output = if merge_adjacent {
                            match &prev_active {
                                Some(pa) => {
                                    pa.len() != active_entries.len()
                                    || !pa.keys().all(|k| active_entries.contains_key(k))
                                }
                                None => true,
                            }
                        } else {
                            true
                        };

                        if should_output {
                            let depth = active_entries.len();
                            let _ = write!(temp_writer, "{}", depth);

                            for s in &sample_order {
                                if let Some(&eidx) = active_entries.get(s) {
                                    let entry = &entries[eidx];
                                    let ref_len = entry.ref_end - entry.ref_start;
                                    let query_len = entry.query_end - entry.query_start;

                                    let window_ref_start = prev.max(entry.ref_start);
                                    let window_ref_end = pos.min(entry.ref_end);

                                    let ratio_start = if ref_len > 0 {
                                        (window_ref_start - entry.ref_start) as f64 / ref_len as f64
                                    } else {
                                        0.0
                                    };
                                    let ratio_end = if ref_len > 0 {
                                        (window_ref_end - entry.ref_start) as f64 / ref_len as f64
                                    } else {
                                        1.0
                                    };

                                    let (wqs, wqe) = if entry.is_reverse {
                                        // Reverse strand: moving forward in ref means moving backward in query
                                        let s = entry.query_end - (query_len as f64 * ratio_start) as i64;
                                        let e = entry.query_end - (query_len as f64 * ratio_end) as i64;
                                        (s, e)
                                    } else {
                                        // Forward strand: simple linear mapping
                                        let s = entry.query_start + (query_len as f64 * ratio_start) as i64;
                                        let e = entry.query_start + (query_len as f64 * ratio_end) as i64;
                                        (s, e)
                                    };
                                    // Always ensure output has start < end
                                    let (window_query_start, window_query_end) = (wqs.min(wqe), wqs.max(wqe));

                                    let _ = write!(temp_writer, "\t{}:{}-{}", entry.query_seq, window_query_start, window_query_end);
                                } else {
                                    let _ = write!(temp_writer, "\tNA");
                                }
                            }
                            let _ = writeln!(temp_writer);
                            window_count += 1;

                            if merge_adjacent {
                                prev_active = Some(active_entries.clone());
                            }
                        }
                    }
                }

                // Update active entries
                let entry = &entries[entry_idx];
                if is_start {
                    active_entries
                        .entry(entry.sample.clone())
                        .and_modify(|existing_idx| {
                            let existing = &entries[*existing_idx];
                            let existing_len = existing.ref_end - existing.ref_start;
                            let new_len = entry.ref_end - entry.ref_start;
                            if new_len > existing_len {
                                *existing_idx = entry_idx;
                            }
                        })
                        .or_insert(entry_idx);
                } else if let Some(&active_idx) = active_entries.get(&entry.sample) {
                    if active_idx == entry_idx {
                        active_entries.remove(&entry.sample);
                    }
                }

                prev_pos = Some(pos);
            }

            // Handle uncovered tail region within this chunk
            let last_covered_pos = prev_pos.unwrap_or(chunk_start);
            if last_covered_pos < chunk_end {
                // Check if this is the last chunk for this sequence
                let is_last_chunk = chunk_end >= *seq_len;
                if is_last_chunk || last_covered_pos < chunk_end {
                    let _ = write!(temp_writer, "1");
                    for s in &sample_order {
                        if s == ref_sample {
                            let _ = write!(temp_writer, "\t{}:{}-{}", ref_seq, last_covered_pos, chunk_end);
                        } else {
                            let _ = write!(temp_writer, "\tNA");
                        }
                    }
                    let _ = writeln!(temp_writer);
                    window_count += 1;
                }
            }

            let _ = temp_writer.flush();
            (ref_idx, chunk_idx, temp_path, window_count)
        })
        .collect();

    // Sort results by (ref_idx, chunk_idx) to ensure correct ordering
    let mut sorted_results = chunk_results;
    sorted_results.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    // Convert to the format expected by merge step
    let temp_files: Vec<(std::path::PathBuf, usize)> = sorted_results
        .into_iter()
        .map(|(_, _, path, count)| (path, count))
        .collect();

    // Merge temp files into final output with global window IDs
    // Also handle deduplication at chunk boundaries
    info!("  Merging {} temp files...", temp_files.len());

    let writer: Box<dyn Write> = if let Some(prefix) = output_prefix {
        let path = format!("{}.depth.tsv", prefix);
        Box::new(std::io::BufWriter::new(std::fs::File::create(&path)?))
    } else {
        Box::new(std::io::BufWriter::new(std::io::stdout()))
    };
    let mut writer = std::io::BufWriter::new(writer);

    // Write header
    write!(writer, "window_id\tdepth")?;
    for sample in &sample_order {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // Helper function to extract coverage signature (depth + which samples have coverage)
    // Returns (depth, Vec<bool>) where bool indicates if sample has coverage (not NA)
    fn get_coverage_signature(line: &str) -> (String, Vec<bool>) {
        let parts: Vec<&str> = line.split('\t').collect();
        let depth = parts.first().map(|s| s.to_string()).unwrap_or_default();
        let has_coverage: Vec<bool> = parts.iter().skip(1).map(|s| *s != "NA").collect();
        (depth, has_coverage)
    }

    let mut global_window_id: usize = 0;
    let mut prev_line: Option<String> = None;
    let mut prev_signature: Option<(String, Vec<bool>)> = None;
    let mut at_chunk_boundary = false;
    let mut merged_count: usize = 0;

    for (temp_path, window_count) in &temp_files {
        if *window_count == 0 {
            at_chunk_boundary = true;
            continue;
        }

        let file = std::fs::File::open(temp_path)?;
        let reader = std::io::BufReader::new(file);
        let mut first_line_of_chunk = true;

        for line in reader.lines() {
            let line = line?;
            let current_signature = get_coverage_signature(&line);

            // Check for duplicate at chunk boundary
            let is_duplicate = if at_chunk_boundary && first_line_of_chunk {
                if let Some(ref prev_sig) = prev_signature {
                    // Same depth and same coverage pattern -> likely duplicate at boundary
                    prev_sig.0 == current_signature.0 && prev_sig.1 == current_signature.1
                } else {
                    false
                }
            } else {
                false
            };

            if is_duplicate {
                // Skip this line as it's a duplicate at chunk boundary
                merged_count += 1;
            } else {
                // Output the previous line if exists
                if let Some(ref prev) = prev_line {
                    writeln!(writer, "{}\t{}", global_window_id, prev)?;
                    global_window_id += 1;
                }
                prev_line = Some(line);
                prev_signature = Some(current_signature);
            }

            first_line_of_chunk = false;
        }

        at_chunk_boundary = true;
    }

    // Output the last line
    if let Some(ref prev) = prev_line {
        writeln!(writer, "{}\t{}", global_window_id, prev)?;
        global_window_id += 1;
    }

    writer.flush()?;

    // Temp directory and files are automatically cleaned up when temp_dir goes out of scope
    drop(temp_dir);

    if merged_count > 0 {
        info!(
            "  Merged {} duplicate windows at chunk boundaries",
            merged_count
        );
    }
    info!("Output complete: {} windows", global_window_id);

    Ok(())
}

/// Global streaming depth computation - processes all samples as reference one by one.
/// Each sample's sequences are processed once, ensuring complete coverage without overlap.
///
/// Algorithm:
/// 1. Get all unique samples from the index
/// 2. For each sample (in sorted order), compute coverage for its sequences
/// 3. Output combined results with reference sample first if specified
pub fn compute_depth_global_streaming(
    impg: &Impg,
    separator: &str,
    output_prefix: Option<&str>,
    merge_adjacent: bool,
    fai_list: Option<&str>,
    ref_sample_first: Option<&str>,  // Optional: put this sample's columns first
) -> io::Result<()> {
    info!("Computing depth (global streaming mode)...");

    // Phase 1: Get all samples and their sequences
    info!("Phase 1: Collecting all samples and sequences...");

    let seq_lengths = SequenceLengths::from_impg(impg, separator);
    let mut all_samples: Vec<String> = seq_lengths.get_samples();

    // Sort samples: ref_sample_first comes first if specified, then alphabetical
    all_samples.sort_by(|a, b| {
        if let Some(ref_first) = ref_sample_first {
            if a == ref_first {
                return std::cmp::Ordering::Less;
            } else if b == ref_first {
                return std::cmp::Ordering::Greater;
            }
        }
        a.cmp(b)
    });

    info!("  Found {} samples to process", all_samples.len());

    // Build mapping of sample -> sequences
    let mut sample_seq_ids: FxHashMap<String, Vec<(u32, String, i64)>> = FxHashMap::default();
    for id in 0..impg.seq_index.len() as u32 {
        if let Some(name) = impg.seq_index.get_name(id) {
            let sample = extract_sample(name, separator);
            if let Some(len) = impg.seq_index.get_len_from_id(id) {
                sample_seq_ids
                    .entry(sample.to_string())
                    .or_default()
                    .push((id, name.to_string(), len as i64));
            }
        }
    }

    // Load FAI lengths if provided
    let fai_lengths: Option<FxHashMap<String, i64>> = if let Some(fai_path) = fai_list {
        info!("  Loading sequence lengths from FAI files...");
        let fai_seq_lengths = SequenceLengths::from_fai_list(fai_path, separator)?;
        info!(
            "  FAI: {} samples, {} sequences",
            fai_seq_lengths.sample_to_seqs.len(),
            fai_seq_lengths.lengths.len()
        );
        Some(fai_seq_lengths.lengths)
    } else {
        None
    };

    // Phase 2: Process each sample's sequences
    info!("Phase 2: Scanning index for coverage...");

    // Get all target IDs once
    let all_target_ids: Vec<u32> = impg.get_all_target_ids();
    let total_trees = all_target_ids.len();
    info!(
        "  {} trees to scan, {} threads available",
        total_trees,
        rayon::current_num_threads()
    );

    // Build a global mapping: for each sequence ID, record all overlapping alignments
    // This is done once and used for all samples
    let processed_trees = AtomicUsize::new(0);

    // Collect all coverage entries in parallel
    // Key: (ref_seq_name) -> Vec<CoverageEntry>
    let parallel_results: Vec<FxHashMap<String, Vec<CoverageEntry>>> = all_target_ids
        .par_iter()
        .map(|&target_id| {
            let mut local_coverage: FxHashMap<String, Vec<CoverageEntry>> = FxHashMap::default();

            let count = processed_trees.fetch_add(1, Ordering::Relaxed);
            if count % 1000 == 0 || count < 10 {
                debug!(
                    "    Progress: {}/{} trees (thread {:?})",
                    count,
                    total_trees,
                    std::thread::current().id()
                );
            }

            let target_name = match impg.seq_index.get_name(target_id) {
                Some(n) => n.to_string(),
                None => return local_coverage,
            };
            let target_sample = extract_sample(&target_name, separator);

            // Load tree
            let tree = match impg.load_tree_no_cache(target_id) {
                Some(t) => t,
                None => return local_coverage,
            };

            // For each alignment in this tree
            for interval in tree.iter() {
                let query_id = interval.metadata.query_id();
                let query_name = match impg.seq_index.get_name(query_id) {
                    Some(n) => n.to_string(),
                    None => continue,
                };
                let query_sample = extract_sample(&query_name, separator);

                // Skip same-sample alignments (self-coverage added separately)
                if target_sample == query_sample {
                    continue;
                }

                let target_start = interval.first as i64;
                let target_end = interval.last as i64;
                let is_reverse = interval.metadata.query_start() > interval.metadata.query_end();
                let query_start = interval.metadata.query_start().min(interval.metadata.query_end()) as i64;
                let query_end = interval.metadata.query_start().max(interval.metadata.query_end()) as i64;

                // Record coverage bidirectionally for complete coverage
                // 1. On target's sequence: query sample provides coverage
                local_coverage
                    .entry(target_name.clone())
                    .or_default()
                    .push(CoverageEntry {
                        sample: query_sample.to_string(),
                        query_seq: query_name.clone(),
                        ref_start: target_start,
                        ref_end: target_end,
                        query_start,
                        query_end,
                        is_reverse,
                    });

                // 2. On query's sequence: target sample provides coverage (reverse direction)
                // This ensures the window can be output from the alphabetically-first sample's perspective
                local_coverage
                    .entry(query_name)  // Move instead of clone
                    .or_default()
                    .push(CoverageEntry {
                        sample: target_sample.to_string(),
                        query_seq: target_name.clone(),
                        ref_start: query_start,
                        ref_end: query_end,
                        query_start: target_start,
                        query_end: target_end,
                        is_reverse, // same relative orientation
                    });
            }

            local_coverage
        })
        .collect();

    // Merge parallel results
    info!("  Merging coverage data...");
    let mut coverage: FxHashMap<String, Vec<CoverageEntry>> = FxHashMap::default();
    for local in parallel_results {
        for (seq_name, entries) in local {
            coverage
                .entry(seq_name)
                .or_default()
                .extend(entries);
        }
    }

    // Add self-coverage for each sequence
    for (sample, seqs) in &sample_seq_ids {
        for (_, seq_name, seq_len) in seqs {
            let len = if let Some(ref fai) = fai_lengths {
                *fai.get(seq_name).unwrap_or(seq_len)
            } else {
                *seq_len
            };

            coverage
                .entry(seq_name.clone())
                .or_default()
                .push(CoverageEntry {
                    sample: sample.clone(),
                    query_seq: seq_name.clone(),
                    ref_start: 0,
                    ref_end: len,
                    query_start: 0,
                    query_end: len,
                    is_reverse: false, // self-coverage is always forward
                });
        }
    }

    info!(
        "  Coverage data collected for {} sequences",
        coverage.len()
    );

    // Phase 3: Generate depth output
    info!("Phase 3: Generating depth output (parallel with chunking)...");

    // Chunk size for parallel processing
    const CHUNK_SIZE: i64 = 10_000_000;

    // Build list of all sequences to process, grouped by sample
    // Process in sample order (ref_sample_first if specified)
    let mut all_seq_list: Vec<(String, String, i64)> = Vec::new();  // (sample, seq_name, seq_len)
    for sample in &all_samples {
        if let Some(seqs) = sample_seq_ids.get(sample) {
            for (_, seq_name, seq_len) in seqs {
                let len = if let Some(ref fai) = fai_lengths {
                    *fai.get(seq_name).unwrap_or(seq_len)
                } else {
                    *seq_len
                };
                all_seq_list.push((sample.clone(), seq_name.clone(), len));
            }
        }
    }

    // Sort sequences: by sample order, then by sequence name
    all_seq_list.sort_by(|a, b| {
        let sample_order_a = all_samples.iter().position(|s| s == &a.0).unwrap_or(usize::MAX);
        let sample_order_b = all_samples.iter().position(|s| s == &b.0).unwrap_or(usize::MAX);
        sample_order_a.cmp(&sample_order_b).then_with(|| a.1.cmp(&b.1))
    });

    // Build chunks for all sequences
    let mut all_chunks: Vec<(usize, usize, String, i64, i64)> = Vec::new();  // (seq_idx, chunk_idx, seq_name, chunk_start, chunk_end)
    for (seq_idx, (_, seq_name, seq_len)) in all_seq_list.iter().enumerate() {
        let num_chunks = ((*seq_len + CHUNK_SIZE - 1) / CHUNK_SIZE) as usize;
        for chunk_idx in 0..num_chunks.max(1) {
            let chunk_start = (chunk_idx as i64) * CHUNK_SIZE;
            let chunk_end = ((chunk_idx as i64 + 1) * CHUNK_SIZE).min(*seq_len);
            all_chunks.push((seq_idx, chunk_idx, seq_name.clone(), chunk_start, chunk_end));
        }
    }

    info!(
        "  Processing {} chunks across {} sequences with {} threads",
        all_chunks.len(),
        all_seq_list.len(),
        rayon::current_num_threads()
    );

    // Create temp directory
    let temp_dir = tempfile::tempdir()?;
    let temp_dir_path = temp_dir.path().to_path_buf();

    // Process chunks in parallel
    let processed_chunks = AtomicUsize::new(0);
    let total_chunks = all_chunks.len();

    let chunk_results: Vec<(usize, usize, std::path::PathBuf, usize)> = all_chunks
        .par_iter()
        .map(|(seq_idx, chunk_idx, seq_name, chunk_start, chunk_end)| {
            let temp_path = temp_dir_path.join(format!("depth_{:06}_{:06}.tmp", seq_idx, chunk_idx));
            let mut window_count = 0usize;

            let count = processed_chunks.fetch_add(1, Ordering::Relaxed);
            if count % 100 == 0 || count < 10 {
                debug!(
                    "    Progress: {}/{} chunks (thread {:?})",
                    count,
                    total_chunks,
                    std::thread::current().id()
                );
            }

            // Extract the reference sample from the sequence name
            let ref_sample = extract_sample(seq_name, separator);

            let entries = match coverage.get(seq_name) {
                Some(e) => e,
                None => return (*seq_idx, *chunk_idx, temp_path, 0),
            };

            if entries.is_empty() {
                return (*seq_idx, *chunk_idx, temp_path, 0);
            }

            // Find entries that overlap this chunk
            let chunk_entries: Vec<(usize, &CoverageEntry)> = entries
                .iter()
                .enumerate()
                .filter(|(_, e)| e.ref_start < *chunk_end && e.ref_end > *chunk_start)
                .collect();

            if chunk_entries.is_empty() {
                return (*seq_idx, *chunk_idx, temp_path, 0);
            }

            // Create temp file
            let file = match std::fs::File::create(&temp_path) {
                Ok(f) => f,
                Err(_) => return (*seq_idx, *chunk_idx, temp_path, 0),
            };
            let mut temp_writer = std::io::BufWriter::new(file);

            // Build sweep-line events
            let mut events: Vec<(i64, bool, usize)> = Vec::new();
            for &(orig_idx, entry) in &chunk_entries {
                let effective_start = entry.ref_start.max(*chunk_start);
                let effective_end = entry.ref_end.min(*chunk_end);
                if effective_start < effective_end {
                    events.push((effective_start, true, orig_idx));
                    events.push((effective_end, false, orig_idx));
                }
            }

            events.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

            // Sweep-line processing
            let mut active_entries: FxHashMap<String, usize> = FxHashMap::default();
            let mut prev_pos: Option<i64> = None;
            let mut prev_active: Option<FxHashMap<String, usize>> = None;

            for (pos, is_start, entry_idx) in events {
                if let Some(prev) = prev_pos {
                    if pos > prev && !active_entries.is_empty() {
                        let should_output = if merge_adjacent {
                            match &prev_active {
                                Some(pa) => {
                                    pa.len() != active_entries.len()
                                        || !pa.keys().all(|k| active_entries.contains_key(k))
                                }
                                None => true,
                            }
                        } else {
                            true
                        };

                        if should_output {
                            // Deduplication: only output if ref_sample is the first (alphabetically)
                            // among all samples covering this window. This ensures each homologous
                            // region is output exactly once.
                            // Note: ref_sample is not in active_entries (it's implicit), so we
                            // compare it against the first covering sample.
                            let first_covering = active_entries.keys().min().map(|s| s.as_str());
                            let is_first = match first_covering {
                                Some(first) => ref_sample.as_str() <= first,
                                None => true, // Only ref_sample covers this region
                            };
                            if !is_first {
                                prev_pos = Some(pos);
                                continue;
                            }

                            let depth = active_entries.len();
                            let _ = write!(temp_writer, "{}", depth);

                            for s in &all_samples {
                                // For the reference sample, use window coordinates directly
                                // This ensures contiguous coverage without overlaps
                                if s.as_str() == ref_sample {
                                    let _ = write!(
                                        temp_writer,
                                        "\t{}:{}-{}",
                                        seq_name, prev, pos
                                    );
                                } else if let Some(&eidx) = active_entries.get(s) {
                                    let entry = &entries[eidx];
                                    let ref_len = entry.ref_end - entry.ref_start;
                                    let query_len = entry.query_end - entry.query_start;

                                    let window_ref_start = prev.max(entry.ref_start);
                                    let window_ref_end = pos.min(entry.ref_end);

                                    let ratio_start = if ref_len > 0 {
                                        (window_ref_start - entry.ref_start) as f64 / ref_len as f64
                                    } else {
                                        0.0
                                    };
                                    let ratio_end = if ref_len > 0 {
                                        (window_ref_end - entry.ref_start) as f64 / ref_len as f64
                                    } else {
                                        1.0
                                    };

                                    let (wqs, wqe) = if entry.is_reverse {
                                        // Reverse strand: moving forward in ref means moving backward in query
                                        let s = entry.query_end - (query_len as f64 * ratio_start) as i64;
                                        let e = entry.query_end - (query_len as f64 * ratio_end) as i64;
                                        (s, e)
                                    } else {
                                        // Forward strand: simple linear mapping
                                        let s = entry.query_start + (query_len as f64 * ratio_start) as i64;
                                        let e = entry.query_start + (query_len as f64 * ratio_end) as i64;
                                        (s, e)
                                    };
                                    // Always ensure output has start < end
                                    let (window_query_start, window_query_end) = (wqs.min(wqe), wqs.max(wqe));

                                    let _ = write!(
                                        temp_writer,
                                        "\t{}:{}-{}",
                                        entry.query_seq, window_query_start, window_query_end
                                    );
                                } else {
                                    let _ = write!(temp_writer, "\tNA");
                                }
                            }
                            let _ = writeln!(temp_writer);
                            window_count += 1;

                            if merge_adjacent {
                                prev_active = Some(active_entries.clone());
                            }
                        }
                    }
                }

                // Update active entries
                let entry = &entries[entry_idx];
                if is_start {
                    active_entries
                        .entry(entry.sample.clone())
                        .and_modify(|existing_idx| {
                            let existing = &entries[*existing_idx];
                            let existing_len = existing.ref_end - existing.ref_start;
                            let new_len = entry.ref_end - entry.ref_start;
                            if new_len > existing_len {
                                *existing_idx = entry_idx;
                            }
                        })
                        .or_insert(entry_idx);
                } else if let Some(&active_idx) = active_entries.get(&entry.sample) {
                    if active_idx == entry_idx {
                        active_entries.remove(&entry.sample);
                    }
                }

                prev_pos = Some(pos);
            }

            let _ = temp_writer.flush();
            (*seq_idx, *chunk_idx, temp_path, window_count)
        })
        .collect();

    // Sort results
    let mut sorted_results = chunk_results;
    sorted_results.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    let temp_files: Vec<(std::path::PathBuf, usize)> = sorted_results
        .into_iter()
        .map(|(_, _, path, count)| (path, count))
        .collect();

    // Merge and output
    info!("  Merging {} temp files...", temp_files.len());

    let writer: Box<dyn Write> = if let Some(prefix) = output_prefix {
        let path = format!("{}.depth.tsv", prefix);
        Box::new(std::io::BufWriter::new(std::fs::File::create(&path)?))
    } else {
        Box::new(std::io::BufWriter::new(std::io::stdout()))
    };
    let mut writer = std::io::BufWriter::new(writer);

    // Write header
    write!(writer, "window_id\tdepth")?;
    for sample in &all_samples {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;

    // Merge with deduplication
    fn get_coverage_signature(line: &str) -> (String, Vec<bool>) {
        let parts: Vec<&str> = line.split('\t').collect();
        let depth = parts.first().map(|s| s.to_string()).unwrap_or_default();
        let has_coverage: Vec<bool> = parts.iter().skip(1).map(|s| *s != "NA").collect();
        (depth, has_coverage)
    }

    let mut global_window_id: usize = 0;
    let mut prev_line: Option<String> = None;
    let mut prev_signature: Option<(String, Vec<bool>)> = None;
    let mut at_chunk_boundary = false;
    let mut merged_count: usize = 0;

    for (temp_path, window_count) in &temp_files {
        if *window_count == 0 {
            at_chunk_boundary = true;
            continue;
        }

        let file = std::fs::File::open(temp_path)?;
        let reader = std::io::BufReader::new(file);
        let mut first_line_of_chunk = true;

        for line in reader.lines() {
            let line = line?;
            let current_signature = get_coverage_signature(&line);

            let is_duplicate = if at_chunk_boundary && first_line_of_chunk {
                if let Some(ref prev_sig) = prev_signature {
                    prev_sig.0 == current_signature.0 && prev_sig.1 == current_signature.1
                } else {
                    false
                }
            } else {
                false
            };

            if is_duplicate {
                merged_count += 1;
            } else {
                if let Some(ref prev) = prev_line {
                    writeln!(writer, "{}\t{}", global_window_id, prev)?;
                    global_window_id += 1;
                }
                prev_line = Some(line);
                prev_signature = Some(current_signature);
            }

            first_line_of_chunk = false;
        }

        at_chunk_boundary = true;
    }

    if let Some(ref prev) = prev_line {
        writeln!(writer, "{}\t{}", global_window_id, prev)?;
        global_window_id += 1;
    }

    writer.flush()?;
    drop(temp_dir);

    if merged_count > 0 {
        info!(
            "  Merged {} duplicate windows at chunk boundaries",
            merged_count
        );
    }
    info!("Output complete: {} windows", global_window_id);

    Ok(())
}

/// Load sequence lengths from FAI files for a specific sample
fn load_fai_lengths_for_sample(
    fai_list_path: &str,
    ref_sample: &str,
    separator: &str,
) -> io::Result<FxHashMap<String, i64>> {
    let mut lengths: FxHashMap<String, i64> = FxHashMap::default();

    let content = std::fs::read_to_string(fai_list_path).map_err(|e| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!("Failed to read FAI list file '{}': {}", fai_list_path, e),
        )
    })?;

    for line in content.lines() {
        let fai_path = line.trim();
        if fai_path.is_empty() || fai_path.starts_with('#') {
            continue;
        }

        let fai_content = std::fs::read_to_string(fai_path).map_err(|e| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Failed to read FAI file '{}': {}", fai_path, e),
            )
        })?;

        for fai_line in fai_content.lines() {
            let parts: Vec<&str> = fai_line.split('\t').collect();
            if parts.len() >= 2 {
                let seq_name = parts[0];
                let seq_len: i64 = parts[1].parse().unwrap_or(0);
                let sample = extract_sample(seq_name, separator);

                if sample == ref_sample {
                    lengths.insert(seq_name.to_string(), seq_len);
                }
            }
        }
    }

    Ok(lengths)
}

// ============================================================================
// New depth functions: stats mode and region query mode
// ============================================================================

/// Alignment info for sweep-line processing (tracks all alignments per sample)
#[derive(Debug, Clone)]
struct AlignmentInfoMulti {
    sample: String,
    query_name: String,
    query_start: i64,
    query_end: i64,
    target_start: i64,
    target_end: i64,
    is_reverse: bool,
}

/// Event for sweep-line algorithm (used in new depth functions)
#[derive(Debug, Clone, PartialEq, Eq)]
struct DepthEventMulti {
    position: i64,
    is_start: bool,
    sample: String,
    alignment_idx: usize,
}

impl Ord for DepthEventMulti {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| other.is_start.cmp(&self.is_start)) // Starts before ends at same position
            .then_with(|| self.sample.cmp(&other.sample))
    }
}

impl PartialOrd for DepthEventMulti {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Compute global depth statistics across all sequences
///
/// Algorithm:
/// 1. For each sample in the filter (or all samples if no filter):
///    a. For each sequence belonging to this sample:
///       - Find all alignments FROM other samples TO this sequence (forward)
///       - Find all alignments FROM this sequence TO other samples (reverse)
///    b. Use sweep-line algorithm to compute depth intervals
/// 2. Aggregate statistics across all sequences
/// 3. Output summary and per-depth BED files
pub fn compute_depth_stats(
    impg: &Impg,
    config: &DepthConfig,
    separator: &str,
    output_prefix: &str,
    sample_filter: Option<&SampleFilter>,
    fai_list: Option<&str>,
    ref_sample: Option<&str>,
) -> io::Result<DepthStats> {
    info!("Computing depth statistics (stats mode) - parallelized");

    // Build sequence lengths
    let seq_lengths = if let Some(fai_path) = fai_list {
        info!("Loading sequence lengths from FAI files...");
        SequenceLengths::from_fai_list(fai_path, separator)?
    } else {
        SequenceLengths::from_impg(impg, separator)
    };

    info!(
        "Found {} samples, {} sequences",
        seq_lengths.sample_to_seqs.len(),
        seq_lengths.lengths.len()
    );

    // Get samples to process
    // If ref_sample is specified, only process that sample's sequences (single reference mode)
    // Otherwise, apply sample_filter if present (global mode)
    let samples_to_process: Vec<String> = if let Some(ref_name) = ref_sample {
        info!("Single reference mode: only processing sequences from '{}'", ref_name);
        if seq_lengths.sample_to_seqs.contains_key(ref_name) {
            vec![ref_name.to_string()]
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Reference sample '{}' not found in alignment data", ref_name),
            ));
        }
    } else {
        seq_lengths
            .get_samples()
            .into_iter()
            .filter(|s| sample_filter.map_or(true, |f| f.includes(s)))
            .collect()
    };

    if ref_sample.is_some() {
        info!(
            "Single reference mode: processing {} sample",
            samples_to_process.len()
        );
    } else if sample_filter.is_some() {
        info!(
            "Sample filter active: processing {} of {} samples",
            samples_to_process.len(),
            seq_lengths.sample_to_seqs.len()
        );
    }

    // Flatten all sequences to process into a single list for parallel iteration
    let all_seqs: Vec<(String, String, i64)> = samples_to_process
        .iter()
        .flat_map(|sample| {
            seq_lengths
                .get_seqs_of_sample(sample)
                .iter()
                .filter_map(|seq| {
                    seq_lengths.get_length(seq).map(|len| (sample.clone(), seq.clone(), len))
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let total_seqs = all_seqs.len();
    info!("Processing {} sequences in parallel", total_seqs);

    let processed_counter = AtomicUsize::new(0);

    // Process sequences in parallel and collect partial stats
    let partial_stats: Vec<DepthStats> = all_seqs
        .par_iter()
        .map(|(anchor_sample, anchor_seq, seq_len)| {
            let mut local_stats = DepthStats::new();

            let count = processed_counter.fetch_add(1, Ordering::Relaxed) + 1;
            if count % 50 == 0 || count == total_seqs {
                debug!("Progress: {}/{} sequences processed", count, total_seqs);
            }

            let seq_id = match impg.seq_index.get_id(anchor_seq) {
                Some(id) => id,
                None => {
                    // Sequence not in alignment index, depth = 1 (self only)
                    local_stats.add_interval(anchor_seq, 0, *seq_len, 1);
                    return local_stats;
                }
            };

            // Collect samples covering this sequence (using lightweight representation)
            let mut sample_intervals: Vec<(i64, i64, String)> = Vec::new();

            // Forward direction: where anchor_seq is TARGET
            let overlaps = if config.transitive {
                if config.transitive_dfs {
                    impg.query_transitive_dfs(
                        seq_id,
                        0,
                        *seq_len as i32,
                        None,
                        config.max_depth,
                        config.min_transitive_len,
                        config.min_distance_between_ranges,
                        None,
                        false,
                        None,
                        None,
                        config.approximate_mode,
                        None,
                    )
                } else {
                    impg.query_transitive_bfs(
                        seq_id,
                        0,
                        *seq_len as i32,
                        None,
                        config.max_depth,
                        config.min_transitive_len,
                        config.min_distance_between_ranges,
                        None,
                        false,
                        None,
                        None,
                        config.approximate_mode,
                        None,
                    )
                }
            } else {
                impg.query(seq_id, 0, *seq_len as i32, false, None, None, config.approximate_mode)
            };

            for overlap in &overlaps {
                let query_interval = &overlap.0;
                let target_interval = &overlap.2;

                let query_name = match impg.seq_index.get_name(query_interval.metadata) {
                    Some(name) => name,
                    None => continue,
                };
                let sample = extract_sample(query_name, separator);

                // Apply sample filter
                if let Some(filter) = sample_filter {
                    if !filter.includes(&sample) {
                        continue;
                    }
                }

                let target_start = target_interval.first.min(target_interval.last) as i64;
                let target_end = target_interval.first.max(target_interval.last) as i64;

                sample_intervals.push((target_start, target_end, sample));
            }

            // Reverse direction: where anchor_seq is QUERY
            let reverse_alignments = impg.query_reverse_for_depth(seq_id);
            for (ref_start, ref_end, target_id) in reverse_alignments {
                if let Some(target_name) = impg.seq_index.get_name(target_id) {
                    let sample = extract_sample(target_name, separator);

                    // Skip self
                    if sample == *anchor_sample {
                        continue;
                    }

                    // Apply sample filter
                    if let Some(filter) = sample_filter {
                        if !filter.includes(&sample) {
                            continue;
                        }
                    }

                    sample_intervals.push((ref_start as i64, ref_end as i64, sample));
                }
            }

            // Add self (anchor sample)
            sample_intervals.push((0, *seq_len, anchor_sample.clone()));

            // Fast sweep-line to compute depth (only counts unique samples, no position tracking)
            let depth_windows = compute_sweep_line_depth_fast(anchor_seq, *seq_len, &sample_intervals);

            for (start, end, depth) in depth_windows {
                local_stats.add_interval(anchor_seq, start, end, depth);
            }

            local_stats
        })
        .collect();

    // Merge all partial stats
    info!("Merging results from {} sequences...", partial_stats.len());
    let mut stats = DepthStats::new();
    for partial in partial_stats {
        stats.merge(&partial);
    }

    // Write outputs
    info!("Writing depth statistics...");

    // Write summary
    let summary_path = format!("{}.summary.txt", output_prefix);
    let summary_file = std::fs::File::create(&summary_path)?;
    let mut summary_writer = BufWriter::new(summary_file);
    stats.write_summary(&mut summary_writer)?;
    summary_writer.flush()?;
    info!("Wrote summary to {}", summary_path);

    // Also print to stdout
    stats.write_summary(&mut std::io::stdout())?;

    // Write per-depth BED files
    stats.write_depth_bed_files(output_prefix)?;

    info!(
        "Stats complete: {} total bases, max depth = {}",
        stats.total_bases,
        stats.max_depth()
    );

    Ok(stats)
}

/// Compute depth statistics with sample tracking for combined output
/// This is similar to compute_depth_stats but tracks which samples cover each interval
/// merge_tolerance: fraction (0.0-1.0) for merging adjacent intervals with similar depth
///                  e.g., 0.05 = 5% tolerance, merge if (max-min)/max <= 0.05
pub fn compute_depth_stats_with_samples(
    impg: &Impg,
    config: &DepthConfig,
    separator: &str,
    output_prefix: &str,
    sample_filter: Option<&SampleFilter>,
    fai_list: Option<&str>,
    ref_sample: Option<&str>,
    merge_tolerance: f64,
) -> io::Result<DepthStatsWithSamples> {
    info!("Computing depth statistics with sample tracking (combined output mode)");

    // Build sequence lengths
    let seq_lengths = if let Some(fai_path) = fai_list {
        info!("Loading sequence lengths from FAI files...");
        SequenceLengths::from_fai_list(fai_path, separator)?
    } else {
        SequenceLengths::from_impg(impg, separator)
    };

    info!(
        "Found {} samples, {} sequences",
        seq_lengths.sample_to_seqs.len(),
        seq_lengths.lengths.len()
    );

    // Get samples to process
    let samples_to_process: Vec<String> = if let Some(ref_name) = ref_sample {
        info!("Single reference mode: only processing sequences from '{}'", ref_name);
        if seq_lengths.sample_to_seqs.contains_key(ref_name) {
            vec![ref_name.to_string()]
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Reference sample '{}' not found in alignment data", ref_name),
            ));
        }
    } else {
        seq_lengths
            .get_samples()
            .into_iter()
            .filter(|s| sample_filter.map_or(true, |f| f.includes(s)))
            .collect()
    };

    if ref_sample.is_some() {
        info!(
            "Single reference mode: processing {} sample",
            samples_to_process.len()
        );
    } else if sample_filter.is_some() {
        info!(
            "Sample filter active: processing {} of {} samples",
            samples_to_process.len(),
            seq_lengths.sample_to_seqs.len()
        );
    }

    // Flatten all sequences to process
    let all_seqs: Vec<(String, String, i64)> = samples_to_process
        .iter()
        .flat_map(|sample| {
            seq_lengths
                .get_seqs_of_sample(sample)
                .iter()
                .filter_map(|seq| {
                    seq_lengths.get_length(seq).map(|len| (sample.clone(), seq.clone(), len))
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let total_seqs = all_seqs.len();
    info!("Processing {} sequences in parallel (with sample tracking)", total_seqs);

    let processed_counter = AtomicUsize::new(0);

    // Process sequences in parallel and collect partial stats
    let partial_stats: Vec<DepthStatsWithSamples> = all_seqs
        .par_iter()
        .map(|(anchor_sample, anchor_seq, seq_len)| {
            let mut local_stats = DepthStatsWithSamples::new();

            let count = processed_counter.fetch_add(1, Ordering::Relaxed) + 1;
            if count % 50 == 0 || count == total_seqs {
                debug!("Progress: {}/{} sequences processed", count, total_seqs);
            }

            let seq_id = match impg.seq_index.get_id(anchor_seq) {
                Some(id) => id,
                None => {
                    // Sequence not in alignment index, depth = 1 (self only)
                    local_stats.add_interval(anchor_seq, 0, *seq_len, 1, vec![anchor_sample.clone()]);
                    return local_stats;
                }
            };

            // Collect samples covering this sequence
            let mut sample_intervals: Vec<(i64, i64, String)> = Vec::new();

            // Forward direction: where anchor_seq is TARGET
            let overlaps = if config.transitive {
                if config.transitive_dfs {
                    impg.query_transitive_dfs(
                        seq_id,
                        0,
                        *seq_len as i32,
                        None,
                        config.max_depth,
                        config.min_transitive_len,
                        config.min_distance_between_ranges,
                        None,
                        false,
                        None,
                        None,
                        config.approximate_mode,
                        None,
                    )
                } else {
                    impg.query_transitive_bfs(
                        seq_id,
                        0,
                        *seq_len as i32,
                        None,
                        config.max_depth,
                        config.min_transitive_len,
                        config.min_distance_between_ranges,
                        None,
                        false,
                        None,
                        None,
                        config.approximate_mode,
                        None,
                    )
                }
            } else {
                impg.query(seq_id, 0, *seq_len as i32, false, None, None, config.approximate_mode)
            };

            for overlap in &overlaps {
                let query_interval = &overlap.0;
                let target_interval = &overlap.2;

                let query_name = match impg.seq_index.get_name(query_interval.metadata) {
                    Some(name) => name,
                    None => continue,
                };
                let sample = extract_sample(query_name, separator);

                // Apply sample filter
                if let Some(filter) = sample_filter {
                    if !filter.includes(&sample) {
                        continue;
                    }
                }

                let target_start = target_interval.first.min(target_interval.last) as i64;
                let target_end = target_interval.first.max(target_interval.last) as i64;

                sample_intervals.push((target_start, target_end, sample));
            }

            // Reverse direction: where anchor_seq is QUERY
            let reverse_alignments = impg.query_reverse_for_depth(seq_id);
            for (ref_start, ref_end, target_id) in reverse_alignments {
                if let Some(target_name) = impg.seq_index.get_name(target_id) {
                    let sample = extract_sample(target_name, separator);

                    // Skip self
                    if sample == *anchor_sample {
                        continue;
                    }

                    // Apply sample filter
                    if let Some(filter) = sample_filter {
                        if !filter.includes(&sample) {
                            continue;
                        }
                    }

                    sample_intervals.push((ref_start as i64, ref_end as i64, sample));
                }
            }

            // Add self (anchor sample)
            sample_intervals.push((0, *seq_len, anchor_sample.clone()));

            // Sweep-line with sample tracking
            let depth_windows = compute_sweep_line_depth_with_samples(anchor_seq, *seq_len, &sample_intervals);

            for (start, end, depth, samples) in depth_windows {
                local_stats.add_interval(anchor_seq, start, end, depth, samples);
            }

            local_stats
        })
        .collect();

    // Merge all partial stats
    info!("Merging results from {} sequences...", partial_stats.len());
    let mut stats = DepthStatsWithSamples::new();
    for partial in partial_stats {
        stats.merge(partial);
    }

    // Write outputs
    info!("Writing depth statistics with sample lists...");
    if merge_tolerance > 0.0 {
        info!("Merge tolerance: {:.1}%", merge_tolerance * 100.0);
    }

    // Write summary
    let summary_path = format!("{}.summary.txt", output_prefix);
    let summary_file = std::fs::File::create(&summary_path)?;
    let mut summary_writer = BufWriter::new(summary_file);
    stats.write_summary(&mut summary_writer)?;
    summary_writer.flush()?;
    info!("Wrote summary to {}", summary_path);

    // Also print to stdout
    stats.write_summary(&mut std::io::stdout())?;

    // Write combined output file (with optional merging)
    stats.write_combined_output(output_prefix, merge_tolerance)?;

    info!(
        "Stats complete: {} total bases, {} intervals, max depth = {}",
        stats.total_bases,
        stats.intervals.len(),
        stats.max_depth()
    );

    Ok(stats)
}

/// Fast sweep-line algorithm for stats mode (only computes depth, no sample position tracking)
/// Returns Vec<(start, end, depth)>
fn compute_sweep_line_depth_fast(
    _anchor_seq: &str,
    seq_len: i64,
    sample_intervals: &[(i64, i64, String)],
) -> Vec<(i64, i64, usize)> {
    if sample_intervals.is_empty() {
        return vec![(0, seq_len, 0)];
    }

    // Create events: (position, is_start, sample_name)
    let mut events: Vec<(i64, bool, &str)> = Vec::with_capacity(sample_intervals.len() * 2);
    for (start, end, sample) in sample_intervals {
        events.push((*start, true, sample.as_str()));
        events.push((*end, false, sample.as_str()));
    }

    // Sort by position, with end events before start events at same position
    events.sort_by(|a, b| {
        a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1))
    });

    let mut results: Vec<(i64, i64, usize)> = Vec::new();
    let mut active_samples: FxHashMap<&str, usize> = FxHashMap::default();
    let mut prev_pos: Option<i64> = None;

    for (pos, is_start, sample) in events {
        if let Some(prev) = prev_pos {
            if pos > prev {
                // Emit window from prev to pos with current depth
                let depth = active_samples.len();
                if depth > 0 {
                    results.push((prev, pos, depth));
                }
            }
        }

        // Update active samples
        if is_start {
            *active_samples.entry(sample).or_insert(0) += 1;
        } else if let Some(count) = active_samples.get_mut(sample) {
            *count -= 1;
            if *count == 0 {
                active_samples.remove(sample);
            }
        }

        prev_pos = Some(pos);
    }

    // Merge adjacent windows with same depth
    if results.len() <= 1 {
        return results;
    }

    let mut merged: Vec<(i64, i64, usize)> = Vec::with_capacity(results.len());
    let mut current = results[0];

    for &(start, end, depth) in &results[1..] {
        if current.1 == start && current.2 == depth {
            // Merge
            current.1 = end;
        } else {
            merged.push(current);
            current = (start, end, depth);
        }
    }
    merged.push(current);

    merged
}

/// Sweep-line algorithm to compute depth windows (tracking all alignments per sample)
fn compute_sweep_line_depth_multi(
    anchor_seq: &str,
    _seq_len: i64,
    alignments: &[AlignmentInfoMulti],
    config: &DepthConfig,
) -> Vec<RegionDepthResult> {
    if alignments.is_empty() {
        return Vec::new();
    }

    // Create events
    let mut events: Vec<DepthEventMulti> = Vec::new();
    for (idx, aln) in alignments.iter().enumerate() {
        events.push(DepthEventMulti {
            position: aln.target_start,
            is_start: true,
            sample: aln.sample.clone(),
            alignment_idx: idx,
        });
        events.push(DepthEventMulti {
            position: aln.target_end,
            is_start: false,
            sample: aln.sample.clone(),
            alignment_idx: idx,
        });
    }
    events.sort();

    // Sweep-line
    let mut results: Vec<RegionDepthResult> = Vec::new();
    let mut active_alignments: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    let mut prev_pos: Option<i64> = None;

    for event in events {
        if let Some(prev) = prev_pos {
            if event.position > prev {
                // Create window from prev to event.position
                let mut result = RegionDepthResult::new(
                    anchor_seq.to_string(),
                    prev,
                    event.position,
                );

                // Collect all active samples and their alignments
                for (sample, aln_indices) in &active_alignments {
                    if aln_indices.is_empty() {
                        continue;
                    }
                    for &idx in aln_indices {
                        let aln = &alignments[idx];
                        // Map target coords to query coords (linear interpolation)
                        let (q_start, q_end) = map_coords_linear(
                            aln.target_start,
                            aln.target_end,
                            aln.query_start,
                            aln.query_end,
                            prev,
                            event.position,
                            aln.is_reverse,
                        );
                        result.add_sample_position(sample, &aln.query_name, q_start, q_end);
                    }
                }

                result.update_depth();
                if result.depth > 0 {
                    results.push(result);
                }
            }
        }

        // Update active alignments
        if event.is_start {
            active_alignments
                .entry(event.sample.clone())
                .or_default()
                .push(event.alignment_idx);
        } else {
            if let Some(alns) = active_alignments.get_mut(&event.sample) {
                alns.retain(|&idx| idx != event.alignment_idx);
            }
        }

        prev_pos = Some(event.position);
    }

    // Merge adjacent windows with same depth if configured
    if config.merge_adjacent {
        results = merge_adjacent_results(results);
    }

    results
}

/// Linear interpolation for coordinate mapping
fn map_coords_linear(
    aln_target_start: i64,
    aln_target_end: i64,
    aln_query_start: i64,
    aln_query_end: i64,
    window_target_start: i64,
    window_target_end: i64,
    is_reverse: bool,
) -> (i64, i64) {
    let target_len = aln_target_end - aln_target_start;
    let query_len = aln_query_end - aln_query_start;

    if target_len == 0 {
        return (aln_query_start, aln_query_end);
    }

    let start_offset = window_target_start.max(aln_target_start) - aln_target_start;
    let end_offset = window_target_end.min(aln_target_end) - aln_target_start;

    let ratio = query_len as f64 / target_len as f64;

    if is_reverse {
        let q_end = aln_query_end - (start_offset as f64 * ratio) as i64;
        let q_start = aln_query_end - (end_offset as f64 * ratio) as i64;
        (q_start.min(q_end), q_start.max(q_end))
    } else {
        let q_start = aln_query_start + (start_offset as f64 * ratio) as i64;
        let q_end = aln_query_start + (end_offset as f64 * ratio) as i64;
        (q_start.min(q_end), q_start.max(q_end))
    }
}

/// Merge adjacent results with same depth
fn merge_adjacent_results(mut results: Vec<RegionDepthResult>) -> Vec<RegionDepthResult> {
    if results.len() <= 1 {
        return results;
    }

    let mut merged: Vec<RegionDepthResult> = Vec::new();
    let mut current = results.remove(0);

    for next in results {
        // Check if adjacent and same depth
        if current.ref_end == next.ref_start
            && current.depth == next.depth
            && current.ref_seq == next.ref_seq
        {
            // Merge
            current.ref_end = next.ref_end;
            // Merge sample positions
            for (sample, positions) in next.sample_positions {
                current
                    .sample_positions
                    .entry(sample)
                    .or_default()
                    .extend(positions);
            }
        } else {
            merged.push(current);
            current = next;
        }
    }
    merged.push(current);

    merged
}

/// Query depth for a specific region with per-sample position tracking
///
/// Algorithm:
/// 1. Parse the target region (seq_name:start-end)
/// 2. Find all alignments overlapping this region
/// 3. For each overlapping alignment, track the sample and its position
/// 4. If sample has multiple alignments, track all of them
/// 5. Output in tabular format with per-sample columns
pub fn query_region_depth(
    impg: &Impg,
    config: &DepthConfig,
    target_seq: &str,
    target_start: i32,
    target_end: i32,
    separator: &str,
    sample_filter: Option<&SampleFilter>,
    sequence_index: Option<&UnifiedSequenceIndex>,
) -> io::Result<Vec<RegionDepthResult>> {
    debug!(
        "Querying region depth: {}:{}-{}",
        target_seq, target_start, target_end
    );

    let target_id = impg.seq_index.get_id(target_seq).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!("Sequence '{}' not found in index", target_seq),
        )
    })?;

    let target_sample = extract_sample(target_seq, separator);

    // Collect all alignments for this region
    let mut alignments: Vec<AlignmentInfoMulti> = Vec::new();

    // Forward direction: where target_seq is TARGET
    let overlaps = if config.transitive {
        if config.transitive_dfs {
            impg.query_transitive_dfs(
                target_id,
                target_start,
                target_end,
                None,
                config.max_depth,
                config.min_transitive_len,
                config.min_distance_between_ranges,
                None,
                false,
                None,
                sequence_index,
                config.approximate_mode,
                None,
            )
        } else {
            impg.query_transitive_bfs(
                target_id,
                target_start,
                target_end,
                None,
                config.max_depth,
                config.min_transitive_len,
                config.min_distance_between_ranges,
                None,
                false,
                None,
                sequence_index,
                config.approximate_mode,
                None,
            )
        }
    } else {
        impg.query(
            target_id,
            target_start,
            target_end,
            false,
            None,
            sequence_index,
            config.approximate_mode,
        )
    };

    for overlap in &overlaps {
        let query_interval = &overlap.0;
        let target_interval = &overlap.2;

        let query_name = match impg.seq_index.get_name(query_interval.metadata) {
            Some(name) => name,
            None => continue,
        };
        let sample = extract_sample(query_name, separator);

        // Apply sample filter
        if let Some(filter) = sample_filter {
            if !filter.includes(&sample) {
                continue;
            }
        }

        let is_reverse = query_interval.first > query_interval.last;
        let query_start = query_interval.first.min(query_interval.last) as i64;
        let query_end = query_interval.first.max(query_interval.last) as i64;
        let t_start = target_interval.first.min(target_interval.last) as i64;
        let t_end = target_interval.first.max(target_interval.last) as i64;

        alignments.push(AlignmentInfoMulti {
            sample,
            query_name: query_name.to_string(),
            query_start,
            query_end,
            target_start: t_start,
            target_end: t_end,
            is_reverse,
        });
    }

    // Reverse direction: where target_seq is QUERY
    let reverse_alignments = impg.query_reverse_for_depth(target_id);
    for (ref_start, ref_end, other_id) in reverse_alignments {
        // Check if overlaps with query region
        if ref_end <= target_start || ref_start >= target_end {
            continue;
        }

        if let Some(other_name) = impg.seq_index.get_name(other_id) {
            let sample = extract_sample(other_name, separator);

            // Skip self
            if sample == target_sample {
                continue;
            }

            // Apply sample filter
            if let Some(filter) = sample_filter {
                if !filter.includes(&sample) {
                    continue;
                }
            }

            alignments.push(AlignmentInfoMulti {
                sample,
                query_name: other_name.to_string(),
                query_start: ref_start as i64,
                query_end: ref_end as i64,
                target_start: ref_start as i64,
                target_end: ref_end as i64,
                is_reverse: false,
            });
        }
    }

    // Add self (target sample) if in filter
    if sample_filter.map_or(true, |f| f.includes(&target_sample)) {
        alignments.push(AlignmentInfoMulti {
            sample: target_sample,
            query_name: target_seq.to_string(),
            query_start: target_start as i64,
            query_end: target_end as i64,
            target_start: target_start as i64,
            target_end: target_end as i64,
            is_reverse: false,
        });
    }

    // Compute depth windows using sweep-line
    let seq_len = impg.seq_index.get_len_from_id(target_id).unwrap_or(0) as i64;
    let results = compute_sweep_line_depth_multi(target_seq, seq_len, &alignments, config);

    // Filter results to only include the query region
    let filtered_results: Vec<RegionDepthResult> = results
        .into_iter()
        .filter(|r| r.ref_end > target_start as i64 && r.ref_start < target_end as i64)
        .map(|mut r| {
            // Clip to query region
            r.ref_start = r.ref_start.max(target_start as i64);
            r.ref_end = r.ref_end.min(target_end as i64);
            r
        })
        .collect();

    Ok(filtered_results)
}

/// Parse target range string in format "seq_name:start-end"
pub fn parse_target_range_depth(target_range: &str) -> io::Result<(String, i32, i32)> {
    // Handle format: seq_name:start-end
    // Note: seq_name may contain colons (e.g., sample#hap#chr)
    let parts: Vec<&str> = target_range.rsplitn(2, ':').collect();
    if parts.len() != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Invalid target range format '{}'. Expected format: seq_name:start-end",
                target_range
            ),
        ));
    }

    let range_str = parts[0];
    let seq_name = parts[1].to_string();

    let range_parts: Vec<&str> = range_str.split('-').collect();
    if range_parts.len() != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Invalid range format '{}'. Expected format: start-end",
                range_str
            ),
        ));
    }

    let start: i32 = range_parts[0].parse().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Invalid start position: {}", range_parts[0]),
        )
    })?;

    let end: i32 = range_parts[1].parse().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Invalid end position: {}", range_parts[1]),
        )
    })?;

    if start >= end {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Start ({}) must be less than end ({})", start, end),
        ));
    }

    Ok((seq_name, start, end))
}

/// Parse BED file for region queries
pub fn parse_bed_file_depth(bed_path: &str) -> io::Result<Vec<(String, i32, i32)>> {
    let content = std::fs::read_to_string(bed_path)?;
    let mut regions = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }

        let seq_name = parts[0].to_string();
        let start: i32 = parts[1].parse().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid start position in BED file: {}", parts[1]),
            )
        })?;
        let end: i32 = parts[2].parse().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid end position in BED file: {}", parts[2]),
            )
        })?;

        regions.push((seq_name, start, end));
    }

    Ok(regions)
}
