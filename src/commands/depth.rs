use crate::impg::{CigarOp, SortedRanges};
use crate::impg_index::{ImpgIndex, RawAlignmentInterval};
use indicatif::{ProgressBar, ProgressStyle};
use crate::sequence_index::UnifiedSequenceIndex;
use bitvec::prelude::*;
use log::{debug, info};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::{self, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// Configuration for the depth command
pub struct DepthConfig {
    pub transitive: bool,
    pub transitive_dfs: bool,
    pub max_depth: u16,
    pub min_transitive_len: i32,
    pub min_distance_between_ranges: i32,
    pub merge_adjacent: bool,
    /// When true, use CIGAR-precise BFS for transitive depth (--use-BFS).
    /// When false (default), use raw-interval BFS with linear interpolation.
    pub use_cigar_bfs: bool,
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
        writeln!(writer, "#seq_name\tlength\tdepth\tsamples")?;

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
                    "{}\t{}\t{}\t{}",
                    interval.seq_name, interval.end - interval.start, interval.depth, samples_str
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
                    "{}\t{}\t{}\t{}",
                    interval.seq_name, interval.end - interval.start, interval.depth, samples_str
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
    pub fn from_impg(impg: &impl ImpgIndex, separator: &str) -> Self {
        let mut lengths = FxHashMap::default();
        let mut sample_to_seqs: FxHashMap<String, Vec<String>> = FxHashMap::default();
        let mut seq_to_sample = FxHashMap::default();

        for id in 0..impg.seq_index().len() as u32 {
            if let (Some(name), Some(len)) = (
                impg.seq_index().get_name(id),
                impg.seq_index().get_len_from_id(id),
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

// ============================================================================
// ID-based compact data structures for memory-efficient depth computation
// ============================================================================

/// Sample index: maps sample names to compact IDs (u16) and vice versa
/// Max 65535 samples, which is sufficient for most pangenome analyses
#[derive(Debug, Clone)]
pub struct SampleIndex {
    name_to_id: FxHashMap<String, u16>,
    id_to_name: Vec<String>,  // Use Vec for O(1) lookup
}

impl SampleIndex {
    /// Create a new empty sample index
    pub fn new() -> Self {
        Self {
            name_to_id: FxHashMap::default(),
            id_to_name: Vec::new(),
        }
    }

    /// Build from SequenceLengths
    pub fn from_seq_lengths(seq_lengths: &SequenceLengths) -> Self {
        let mut samples: Vec<_> = seq_lengths.sample_to_seqs.keys().cloned().collect();
        samples.sort(); // Deterministic ordering

        let mut name_to_id = FxHashMap::default();
        let mut id_to_name = Vec::with_capacity(samples.len());

        for (i, sample) in samples.into_iter().enumerate() {
            name_to_id.insert(sample.clone(), i as u16);
            id_to_name.push(sample);
        }

        Self { name_to_id, id_to_name }
    }

    /// Get sample ID from name
    #[inline]
    pub fn get_id(&self, name: &str) -> Option<u16> {
        self.name_to_id.get(name).copied()
    }

    /// Get sample name from ID
    #[inline]
    pub fn get_name(&self, id: u16) -> Option<&str> {
        self.id_to_name.get(id as usize).map(|s| s.as_str())
    }

    /// Get all sample names in ID order
    pub fn names(&self) -> &[String] {
        &self.id_to_name
    }

    /// Get number of samples
    pub fn len(&self) -> usize {
        self.id_to_name.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.id_to_name.is_empty()
    }
}

impl Default for SampleIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// Compact sequence lengths using numeric IDs
/// Memory usage: O(num_sequences) instead of O(num_sequences × avg_name_length)
#[derive(Debug, Clone)]
pub struct CompactSequenceLengths {
    /// Sequence lengths indexed by seq_id (from impg's SequenceIndex)
    lengths: Vec<i64>,
    /// sample_id -> list of seq_ids belonging to this sample
    sample_to_seqs: Vec<Vec<u32>>,
    /// seq_id -> sample_id (reverse lookup)
    seq_to_sample: Vec<u16>,
    /// Sample index for name<->ID conversion
    sample_index: SampleIndex,
}

impl CompactSequenceLengths {
    /// Build from impg sequence index and separator
    pub fn from_impg(impg: &impl ImpgIndex, separator: &str) -> Self {
        // First pass: collect all samples
        let mut sample_names: FxHashSet<String> = FxHashSet::default();
        for id in 0..impg.seq_index().len() as u32 {
            if let Some(name) = impg.seq_index().get_name(id) {
                sample_names.insert(extract_sample(name, separator));
            }
        }

        // Build sample index
        let mut sorted_samples: Vec<_> = sample_names.into_iter().collect();
        sorted_samples.sort();
        let mut name_to_id = FxHashMap::default();
        let mut id_to_name = Vec::with_capacity(sorted_samples.len());
        for (i, sample) in sorted_samples.into_iter().enumerate() {
            name_to_id.insert(sample.clone(), i as u16);
            id_to_name.push(sample);
        }
        let sample_index = SampleIndex { name_to_id, id_to_name };

        // Second pass: build mappings
        let num_seqs = impg.seq_index().len();
        let mut lengths = vec![0i64; num_seqs];
        let mut sample_to_seqs: Vec<Vec<u32>> = vec![Vec::new(); sample_index.len()];
        let mut seq_to_sample = vec![0u16; num_seqs];

        for seq_id in 0..num_seqs as u32 {
            if let (Some(name), Some(len)) = (
                impg.seq_index().get_name(seq_id),
                impg.seq_index().get_len_from_id(seq_id),
            ) {
                let sample_name = extract_sample(name, separator);
                if let Some(sample_id) = sample_index.get_id(&sample_name) {
                    lengths[seq_id as usize] = len as i64;
                    sample_to_seqs[sample_id as usize].push(seq_id);
                    seq_to_sample[seq_id as usize] = sample_id;
                }
            }
        }

        // Sort sequences within each sample for deterministic ordering
        for seqs in &mut sample_to_seqs {
            seqs.sort();
        }

        CompactSequenceLengths {
            lengths,
            sample_to_seqs,
            seq_to_sample,
            sample_index,
        }
    }

    /// Get sequence length by ID
    #[inline]
    pub fn get_length(&self, seq_id: u32) -> i64 {
        self.lengths.get(seq_id as usize).copied().unwrap_or(0)
    }

    /// Get sample ID for a sequence
    #[inline]
    pub fn get_sample_id(&self, seq_id: u32) -> u16 {
        self.seq_to_sample.get(seq_id as usize).copied().unwrap_or(0)
    }

    /// Get sequences for a sample
    #[inline]
    pub fn get_seqs_of_sample(&self, sample_id: u16) -> &[u32] {
        self.sample_to_seqs.get(sample_id as usize).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Get sample index
    #[inline]
    pub fn sample_index(&self) -> &SampleIndex {
        &self.sample_index
    }

    /// Get number of samples
    pub fn num_samples(&self) -> usize {
        self.sample_index.len()
    }

    /// Get all sample IDs in order
    pub fn sample_ids(&self) -> impl Iterator<Item = u16> {
        0..self.sample_index.len() as u16
    }
}

/// Compact alignment info using numeric IDs
/// Memory usage: ~60% less than SampleAlignmentInfo for large datasets
#[derive(Debug, Clone)]
pub struct CompactAlignmentInfo {
    /// Sample ID (from SampleIndex)
    pub sample_id: u16,
    /// Query sequence ID (from impg's SequenceIndex)
    pub query_id: u32,
    /// Query coordinates
    pub query_start: i64,
    pub query_end: i64,
    /// Target coordinates (on anchor sequence)
    pub target_start: i64,
    pub target_end: i64,
    /// True if alignment is on reverse strand
    pub is_reverse: bool,
}

impl CompactAlignmentInfo {
    /// Create a new compact alignment info
    pub fn new(
        sample_id: u16,
        query_id: u32,
        query_start: i64,
        query_end: i64,
        target_start: i64,
        target_end: i64,
        is_reverse: bool,
    ) -> Self {
        Self {
            sample_id,
            query_id,
            query_start,
            query_end,
            target_start,
            target_end,
            is_reverse,
        }
    }
}

/// Bitmap for tracking active samples in sweep-line algorithm
/// Memory efficient: uses 1 bit per sample instead of hash map entries
/// Supports counting multiple overlaps per sample
#[derive(Debug, Clone)]
pub struct SampleBitmap {
    /// Bit set for sample presence (sample_id -> present)
    bits: BitVec,
    /// Count of active alignments per sample (for overlapping alignments from same sample)
    counts: Vec<u32>,
    /// Number of unique active samples (cached for O(1) depth query)
    active_count: usize,
}

impl SampleBitmap {
    /// Create a new bitmap for the given number of samples
    pub fn new(num_samples: usize) -> Self {
        Self {
            bits: bitvec![0; num_samples],
            counts: vec![0; num_samples],
            active_count: 0,
        }
    }

    /// Add an alignment for a sample (increment count)
    #[inline]
    pub fn add(&mut self, sample_id: u16) {
        let idx = sample_id as usize;
        if idx < self.counts.len() {
            if self.counts[idx] == 0 {
                self.bits.set(idx, true);
                self.active_count += 1;
            }
            self.counts[idx] += 1;
        }
    }

    /// Remove an alignment for a sample (decrement count)
    #[inline]
    pub fn remove(&mut self, sample_id: u16) {
        let idx = sample_id as usize;
        if idx < self.counts.len() && self.counts[idx] > 0 {
            self.counts[idx] -= 1;
            if self.counts[idx] == 0 {
                self.bits.set(idx, false);
                self.active_count -= 1;
            }
        }
    }

    /// Check if a sample is active
    #[inline]
    pub fn is_active(&self, sample_id: u16) -> bool {
        let idx = sample_id as usize;
        idx < self.bits.len() && self.bits[idx]
    }

    /// Get the count of active alignments for a sample
    #[inline]
    pub fn count(&self, sample_id: u16) -> u32 {
        self.counts.get(sample_id as usize).copied().unwrap_or(0)
    }

    /// Get the number of unique active samples (depth)
    #[inline]
    pub fn depth(&self) -> usize {
        self.active_count
    }

    /// Iterate over active sample IDs
    pub fn active_samples(&self) -> impl Iterator<Item = u16> + '_ {
        self.bits.iter_ones().map(|idx| idx as u16)
    }

    /// Clear all active samples
    pub fn clear(&mut self) {
        self.bits.fill(false);
        self.counts.fill(0);
        self.active_count = 0;
    }
}

impl Default for SampleBitmap {
    fn default() -> Self {
        Self::new(0)
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
    /// Uses binary search for O(log n) lookup instead of O(n) linear scan
    pub fn add(&mut self, start: i64, end: i64) {
        if start >= end {
            return;
        }

        if self.intervals.is_empty() {
            self.intervals.push((start, end));
            self.total_length = end - start;
            return;
        }

        // Binary search: find the first interval whose end >= start (could overlap)
        // partition_point returns the first index where predicate is false
        let first_overlap_idx = self.intervals.partition_point(|&(_, iv_end)| iv_end < start);

        // Find the last overlapping interval: first interval whose start > end
        // Only search from first_overlap_idx onward
        let last_overlap_end = first_overlap_idx
            + self.intervals[first_overlap_idx..]
                .partition_point(|&(iv_start, _)| iv_start <= end);

        if first_overlap_idx == last_overlap_end {
            // No overlapping intervals - just insert at correct position
            self.intervals.insert(first_overlap_idx, (start, end));
            self.total_length += end - start;
        } else {
            // Merge with overlapping intervals [first_overlap_idx, last_overlap_end)
            let merged_start = start.min(self.intervals[first_overlap_idx].0);
            let merged_end = end.max(self.intervals[last_overlap_end - 1].1);

            // Calculate removed length
            let removed_length: i64 = self.intervals[first_overlap_idx..last_overlap_end]
                .iter()
                .map(|&(s, e)| e - s)
                .sum();

            // Replace overlapping range with merged interval using drain + insert
            self.intervals.drain(first_overlap_idx..last_overlap_end);
            self.intervals.insert(first_overlap_idx, (merged_start, merged_end));

            // Update total length
            let added_length = merged_end - merged_start;
            self.total_length = self.total_length - removed_length + added_length;
        }
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
    /// Uses binary search for O(log n) lookup of affected range
    pub fn subtract(&mut self, sub_start: i64, sub_end: i64) {
        if sub_start >= sub_end || self.intervals.is_empty() {
            return;
        }

        // Binary search: find the first interval whose end > sub_start (could overlap)
        let first_overlap_idx = self.intervals.partition_point(|&(_, iv_end)| iv_end <= sub_start);

        // Find the last overlapping interval: first interval whose start >= sub_end
        let last_overlap_end = first_overlap_idx
            + self.intervals[first_overlap_idx..]
                .partition_point(|&(iv_start, _)| iv_start < sub_end);

        if first_overlap_idx == last_overlap_end {
            // No overlapping intervals - nothing to subtract
            return;
        }

        // Process only the affected range [first_overlap_idx, last_overlap_end)
        let mut replacement = Vec::new();
        let mut removed_length: i64 = 0;

        for &(iv_start, iv_end) in &self.intervals[first_overlap_idx..last_overlap_end] {
            // Keep left portion if it exists
            if iv_start < sub_start {
                replacement.push((iv_start, sub_start));
            }
            // Keep right portion if it exists
            if iv_end > sub_end {
                replacement.push((sub_end, iv_end));
            }
            // Calculate removed length
            let overlap_start = iv_start.max(sub_start);
            let overlap_end = iv_end.min(sub_end);
            removed_length += overlap_end - overlap_start;
        }

        // Replace affected range with replacement intervals
        self.intervals.splice(first_overlap_idx..last_overlap_end, replacement);
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
        let q_end = aln_query_end - (start_offset as f64 * ratio).round() as i64;
        let q_start = aln_query_end - (end_offset as f64 * ratio).round() as i64;
        (q_start.min(q_end), q_start.max(q_end))
    } else {
        let q_start = aln_query_start + (start_offset as f64 * ratio).round() as i64;
        let q_end = aln_query_start + (end_offset as f64 * ratio).round() as i64;
        (q_start.min(q_end), q_start.max(q_end))
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
// ============================================================================
// Compact depth computation using ID-based data structures
// ============================================================================

/// Compact depth event for sweep-line algorithm using numeric IDs
#[derive(Debug, Clone, PartialEq, Eq)]
struct CompactDepthEvent {
    position: i64,
    is_start: bool,
    sample_id: u16,
    /// Index into the alignment info array
    alignment_idx: usize,
}

impl Ord for CompactDepthEvent {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| other.is_start.cmp(&self.is_start)) // Starts before ends at same position
            .then_with(|| self.sample_id.cmp(&other.sample_id))
    }
}

impl PartialOrd for CompactDepthEvent {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Check if an alignment is a same-sample transitive alignment (intra-genome duplication).
/// Returns true if the alignment should be filtered out.
#[inline]
fn is_self_alignment(query_sample_id: u16, anchor_sample_id: u16, query_id: u32, target_seq_id: u32) -> bool {
    query_sample_id == anchor_sample_id && query_id != target_seq_id
}


// ============================================================================
// Windowed depth computation - fixed window size, sample by sample
// ============================================================================

/// Default window size for windowed depth computation (50kb)
const DEFAULT_WINDOW_SIZE: i64 = 50_000;

// ============================================================================
// Windowed depth structures - sparse storage, streaming output
// ============================================================================

/// Sample position entry (sparse representation)
/// (sample_id, query_seq_id, query_start, query_end)
type SamplePosition = (u16, u32, i64, i64);

/// Depth interval with sparse sample storage
/// Memory: 24 bytes base + 16 bytes per actual sample (vs 16 bytes per ALL samples)
#[derive(Debug, Clone)]
struct SparseDepthInterval {
    start: i64,
    end: i64,
    /// Sparse sample positions - only stores samples that are present
    /// Sorted by sample_id for binary search lookup
    samples: Vec<SamplePosition>,
}

impl SparseDepthInterval {
    #[inline]
    fn depth(&self) -> usize {
        self.samples.len()
    }

    /// Get position for a specific sample (binary search)
    #[inline]
    fn get_sample(&self, sample_id: u16) -> Option<(u32, i64, i64)> {
        self.samples
            .binary_search_by_key(&sample_id, |s| s.0)
            .ok()
            .map(|idx| (self.samples[idx].1, self.samples[idx].2, self.samples[idx].3))
    }
}

/// Split intervals at fixed window boundaries.
/// Each interval crossing a window edge is cut into pieces,
/// with sample query coordinates proportionally adjusted.
fn split_intervals_by_window(intervals: Vec<SparseDepthInterval>, window_size: i64) -> Vec<SparseDepthInterval> {
    let mut result: Vec<SparseDepthInterval> = Vec::new();
    for interval in intervals {
        let mut pos = interval.start;
        while pos < interval.end {
            // Next window boundary after pos
            let next_boundary = ((pos / window_size) + 1) * window_size;
            let chunk_end = interval.end.min(next_boundary);

            if pos >= chunk_end {
                break;
            }

            let interval_len = interval.end - interval.start;
            // Proportionally adjust sample query coordinates for this chunk
            let samples: Vec<SamplePosition> = interval.samples.iter().map(|&(sid, qid, qs, qe)| {
                if interval_len <= 0 {
                    return (sid, qid, qs, qe);
                }
                let q_len = qe - qs;
                let frac_start = (pos - interval.start) as f64 / interval_len as f64;
                let frac_end = (chunk_end - interval.start) as f64 / interval_len as f64;
                let new_qs = qs + (frac_start * q_len as f64) as i64;
                let new_qe = qs + (frac_end * q_len as f64) as i64;
                (sid, qid, new_qs, new_qe)
            }).collect();

            result.push(SparseDepthInterval {
                start: pos,
                end: chunk_end,
                samples,
            });
            pos = chunk_end;
        }
    }
    result
}

/// Merge adjacent intervals with similar depth
/// tolerance: fraction (0.0-1.0), e.g., 0.05 = merge if depth diff <= 5%
fn merge_sparse_intervals(intervals: Vec<SparseDepthInterval>, tolerance: f64) -> Vec<SparseDepthInterval> {
    if intervals.is_empty() || tolerance <= 0.0 {
        return intervals;
    }

    let mut merged: Vec<SparseDepthInterval> = Vec::new();
    let mut iter = intervals.into_iter();

    let mut current = iter.next().unwrap();
    let mut current_min_depth = current.depth();
    let mut current_max_depth = current.depth();

    for interval in iter {
        // Check if can merge: adjacent and within tolerance
        let new_depth = interval.depth();
        let new_min = current_min_depth.min(new_depth);
        let new_max = current_max_depth.max(new_depth);

        let within_tolerance = if new_max == 0 {
            true
        } else {
            (new_max - new_min) as f64 / new_max as f64 <= tolerance
        };

        if interval.start == current.end && within_tolerance {
            // Merge: extend end, union samples, update depth range
            current.end = interval.end;
            current_min_depth = new_min;
            current_max_depth = new_max;

            // Union samples (merge by sample_id, extend coordinates)
            let mut sample_map: FxHashMap<u16, (u32, i64, i64)> = FxHashMap::default();
            for (sid, qid, qs, qe) in current.samples.drain(..) {
                sample_map.insert(sid, (qid, qs, qe));
            }
            for (sid, qid, qs, qe) in interval.samples {
                sample_map.entry(sid)
                    .and_modify(|e| {
                        // Extend query range
                        e.1 = e.1.min(qs);
                        e.2 = e.2.max(qe);
                    })
                    .or_insert((qid, qs, qe));
            }
            current.samples = sample_map.into_iter()
                .map(|(sid, (qid, qs, qe))| (sid, qid, qs, qe))
                .collect();
            current.samples.sort_by_key(|s| s.0);
        } else {
            // Emit current, start new
            merged.push(current);
            current = interval;
            current_min_depth = current.depth();
            current_max_depth = current.depth();
        }
    }

    merged.push(current);
    merged
}

/// Thread-safe processed region tracker using per-sequence locks.
/// Each sequence has its own Mutex, so concurrent access to different
/// sequences has zero contention. This enables fully parallel depth
/// computation across all sequences without batch synchronization barriers.
#[derive(Debug)]
struct ConcurrentProcessedTracker {
    /// Per-sequence IntervalSet behind a Mutex for thread-safe access
    processed: Vec<Mutex<IntervalSet>>,
}

impl ConcurrentProcessedTracker {
    /// Create a new tracker for the given number of sequences
    fn new(num_sequences: usize) -> Self {
        Self {
            processed: (0..num_sequences)
                .map(|_| Mutex::new(IntervalSet::new()))
                .collect(),
        }
    }

    /// Get unprocessed regions for a sequence within [start, end)
    fn get_unprocessed(&self, seq_id: u32, start: i64, end: i64) -> Vec<(i64, i64)> {
        let lock = self.processed[seq_id as usize].lock().unwrap();
        if lock.is_empty() {
            return vec![(start, end)];
        }
        if lock.total_length() >= (end - start) {
            return Vec::new();
        }
        // Subtract processed intervals from [start, end)
        let mut unprocessed = IntervalSet::new_single(start, end);
        for &(s, e) in lock.intervals() {
            if s >= end {
                break;
            }
            if e <= start {
                continue;
            }
            unprocessed.subtract(s, e);
        }
        unprocessed.intervals().to_vec()
    }

    /// Mark a region as processed
    fn mark_processed(&self, seq_id: u32, start: i64, end: i64) {
        let mut lock = self.processed[seq_id as usize].lock().unwrap();
        lock.add(start, end);
    }

    /// Mark a batch of regions as processed: Vec<(seq_id, start, end)>
    fn mark_processed_batch(&self, regions: &[(u32, i64, i64)]) {
        for &(seq_id, start, end) in regions {
            self.mark_processed(seq_id, start, end);
        }
    }
}

// ============================================================================
// Global depth computation - connected-component traversal
// ============================================================================

/// Pre-scan: compute alignment degree for each included sequence.
/// Degree = number of unique OTHER samples with direct alignments to this sequence.
/// Used to automatically identify hub sequences for Phase 1 when --ref is not specified.
fn compute_alignment_degrees(
    impg: &(impl ImpgIndex + Sync),
    compact_lengths: &CompactSequenceLengths,
    seq_included: &[bool],
) -> Vec<u16> {
    (0..seq_included.len() as u32)
        .into_par_iter()
        .map(|seq_id| {
            if !seq_included.get(seq_id as usize).copied().unwrap_or(false) {
                return 0u16;
            }
            let raw_alns = impg.query_raw_intervals(seq_id);
            if raw_alns.is_empty() {
                return 0;
            }
            let self_sample = compact_lengths.get_sample_id(seq_id);
            let mut samples = FxHashSet::default();
            for aln in &raw_alns {
                if !seq_included.get(aln.query_id as usize).copied().unwrap_or(false) {
                    continue;
                }
                let sample_id = compact_lengths.get_sample_id(aln.query_id);
                if sample_id != self_sample {
                    samples.insert(sample_id);
                }
            }
            samples.len().min(u16::MAX as usize) as u16
        })
        .collect()
}

/// Build sequence processing order: sorted by degree descending, then length descending.
/// If ref_sample is specified, that sample's sequences are moved to the front.
fn build_sequence_order(
    impg: &impl ImpgIndex,
    compact_lengths: &CompactSequenceLengths,
    ref_sample_id: Option<u16>,
    seq_included: &[bool],
    degrees: &[u16],
) -> Vec<u32> {
    let num_seqs = impg.seq_index().len();
    let mut seq_order: Vec<(u32, i64, bool, u16)> = (0..num_seqs as u32)
        .filter_map(|seq_id| {
            if !seq_included.get(seq_id as usize).copied().unwrap_or(false) {
                return None;
            }
            let len = compact_lengths.get_length(seq_id);
            if len <= 0 {
                return None;
            }
            let is_ref = ref_sample_id
                .map(|ref_id| compact_lengths.get_sample_id(seq_id) == ref_id)
                .unwrap_or(false);
            let degree = degrees.get(seq_id as usize).copied().unwrap_or(0);
            Some((seq_id, len, is_ref, degree))
        })
        .collect();

    // Sort: ref sample first, then by degree descending, then by length descending
    seq_order.sort_by(|a, b| {
        b.2.cmp(&a.2) // ref sequences first
            .then_with(|| b.3.cmp(&a.3)) // then by degree descending
            .then_with(|| b.1.cmp(&a.1)) // then by length descending
    });

    seq_order.into_iter().map(|(seq_id, _, _, _)| seq_id).collect()
}

/// Result of processing a single anchor region
struct AnchorRegionResult {
    /// Depth intervals for output
    intervals: Vec<SparseDepthInterval>,
    /// All discovered regions to mark as processed: (seq_id, start, end)
    discovered_regions: Vec<(u32, i64, i64)>,
    /// The anchor seq_id
    anchor_seq_id: u32,
    /// The anchor sample_id
    anchor_sample_id: u16,
}

/// A single hit from the depth-specific raw-interval BFS.
/// Stores query and target coordinates plus orientation.
struct DepthBfsHit {
    query_id: u32,
    query_start: i32,
    query_end: i32,
    target_id: u32,
    target_start: i32,
    target_end: i32,
    is_reverse: bool,
}

/// Depth-specific raw-interval BFS using `query_raw_overlapping()` at each hop.
///
/// Unlike the CIGAR-based BFS in `impg.rs`, this uses raw alignment extents with
/// linear interpolation to clip query coordinates. This is faster (no CIGAR walk)
/// and avoids chunk boundary gaps from CIGAR-precise projection, making it better
/// suited for depth computation where base-level precision is not needed.
///
/// Returns all discovered hits (excluding self-referential alignments).
fn depth_transitive_bfs(
    impg: &impl ImpgIndex,
    target_id: u32,
    range_start: i32,
    range_end: i32,
    max_depth: u16,
    min_transitive_len: i32,
    min_distance_between_ranges: i32,
    _use_dfs: bool,
) -> Vec<DepthBfsHit> {
    use std::collections::VecDeque;

    let mut results: Vec<DepthBfsHit> = Vec::new();

    // Lazily allocated visited ranges per sequence (not pre-allocated for all sequences)
    let mut visited_ranges: FxHashMap<u32, SortedRanges> = FxHashMap::default();

    // Initialize visited for the starting target
    let target_len = impg.seq_index().get_len_from_id(target_id).unwrap_or(0) as i32;
    let target_sorted = visited_ranges
        .entry(target_id)
        .or_insert_with(|| SortedRanges::new(target_len, 0));
    let filtered_input_range = target_sorted.insert((range_start, range_end));

    // BFS queue: (seq_id, start, end, current_depth)
    let mut queue: VecDeque<(u32, i32, i32, u16)> = VecDeque::new();
    for (s, e) in filtered_input_range {
        if (e - s).abs() >= min_transitive_len {
            queue.push_back((target_id, s, e, 0));
        }
    }

    while let Some((current_target_id, current_start, current_end, current_depth)) = queue.pop_front() {
        if max_depth > 0 && current_depth >= max_depth {
            continue;
        }

        let raw_alns = impg.query_raw_overlapping(current_target_id, current_start, current_end);

        // Collect next-depth ranges to sort and merge before adding to queue
        let mut next_ranges: Vec<(u32, i32, i32)> = Vec::new();

        for aln in &raw_alns {
            // Skip self-referential (same sequence)
            if aln.query_id == current_target_id {
                continue;
            }

            let aln_target_start = aln.target_start;
            let aln_target_end = aln.target_end;

            // Clip target to current range
            let clipped_target_start = aln_target_start.max(current_start);
            let clipped_target_end = aln_target_end.min(current_end);
            if clipped_target_start >= clipped_target_end {
                continue;
            }

            // Linear interpolation to compute query coordinates
            let target_len_aln = aln_target_end - aln_target_start;
            let query_len_aln = aln.query_end - aln.query_start;
            let ratio = if target_len_aln > 0 {
                query_len_aln as f64 / target_len_aln as f64
            } else {
                1.0
            };

            let (clipped_query_start, clipped_query_end) = if aln.is_reverse {
                let cqe = aln.query_end - ((clipped_target_start - aln_target_start) as f64 * ratio) as i32;
                let cqs = aln.query_end - ((clipped_target_end - aln_target_start) as f64 * ratio) as i32;
                (cqs.min(cqe), cqs.max(cqe))
            } else {
                let cqs = aln.query_start + ((clipped_target_start - aln_target_start) as f64 * ratio) as i32;
                let cqe = aln.query_start + ((clipped_target_end - aln_target_start) as f64 * ratio) as i32;
                (cqs.min(cqe), cqs.max(cqe))
            };

            // Record this hit (clipped coordinates for depth calculation)
            results.push(DepthBfsHit {
                query_id: aln.query_id,
                query_start: clipped_query_start,
                query_end: clipped_query_end,
                target_id: current_target_id,
                target_start: clipped_target_start,
                target_end: clipped_target_end,
                is_reverse: aln.is_reverse,
            });

            // For BFS exploration: use FULL alignment query extent (not clipped).
            // Linear interpolation can underestimate query range when indels are present,
            // causing hop 2+ to miss alignments at range boundaries. Using the full
            // extent guarantees we explore at least as widely as CIGAR-precise BFS.
            // visited_ranges prevents re-exploration, bounding the cost.
            let explore_start = aln.query_start.min(aln.query_end);
            let explore_end = aln.query_start.max(aln.query_end);

            let query_seq_len = impg.seq_index().get_len_from_id(aln.query_id).unwrap_or(0) as i32;
            let ranges = visited_ranges
                .entry(aln.query_id)
                .or_insert_with(|| SortedRanges::new(query_seq_len, 0));

            // Check proximity to existing ranges
            let mut should_add = true;
            if min_distance_between_ranges > 0 {
                let (new_min, new_max) = (explore_start, explore_end);
                let idx = match ranges.ranges.binary_search_by_key(&new_min, |&(start, _)| start) {
                    Ok(i) => i,
                    Err(i) => i,
                };

                if idx > 0 {
                    let (_, prev_end) = ranges.ranges[idx - 1];
                    if (new_min - prev_end).abs() < min_distance_between_ranges {
                        should_add = false;
                    }
                }
                if should_add && idx < ranges.ranges.len() {
                    let (next_start, _) = ranges.ranges[idx];
                    if (next_start - new_max).abs() < min_distance_between_ranges {
                        should_add = false;
                    }
                }
            }

            if should_add {
                let new_ranges = ranges.insert((explore_start, explore_end));
                for (new_start, new_end) in new_ranges {
                    if (new_end - new_start).abs() >= min_transitive_len {
                        next_ranges.push((aln.query_id, new_start, new_end));
                    }
                }
            }
        }

        // Sort and merge contiguous ranges before adding to queue
        if !next_ranges.is_empty() {
            next_ranges.sort_by_key(|(id, start, _)| (*id, *start));

            let mut write = 0;
            for read in 1..next_ranges.len() {
                if next_ranges[write].0 == next_ranges[read].0
                    && next_ranges[write].2 >= next_ranges[read].1
                {
                    next_ranges[write].2 = next_ranges[write].2.max(next_ranges[read].2);
                } else {
                    write += 1;
                    next_ranges.swap(write, read);
                }
            }
            next_ranges.truncate(write + 1);

            for (seq_id, start, end) in next_ranges {
                queue.push_back((seq_id, start, end, current_depth + 1));
            }
        }
    }

    results
}

/// Process a single anchor region using raw-interval BFS for transitive depth.
///
/// This is the default transitive path. Uses `depth_transitive_bfs()` which
/// operates on raw alignment extents with linear interpolation, avoiding
/// CIGAR walk overhead and chunk boundary gaps.
fn process_anchor_region_transitive_raw(
    impg: &impl ImpgIndex,
    config: &DepthConfig,
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
    seq_included: &[bool],
    min_seq_length: i64,
) -> AnchorRegionResult {
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::new();
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    let hits = depth_transitive_bfs(
        impg,
        anchor_seq_id,
        region_start as i32,
        region_end as i32,
        config.max_depth,
        config.min_transitive_len,
        config.min_distance_between_ranges,
        config.transitive_dfs,
    );

    let mut alignments: Vec<CompactAlignmentInfo> = Vec::new();

    // Self alignment (anchor covers itself)
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id, anchor_seq_id,
        region_start, region_end, region_start, region_end, false,
    ));

    // Pass 1: Build anchor coverage map from hop 0 results (target on anchor)
    let mut seq_anchor_coverage: FxHashMap<u32, (i64, i64, i64, i64)> = FxHashMap::default();
    for hit in &hits {
        if hit.target_id == anchor_seq_id {
            let q_start = hit.query_start.min(hit.query_end) as i64;
            let q_end = hit.query_start.max(hit.query_end) as i64;
            let t_start = hit.target_start.min(hit.target_end) as i64;
            let t_end = hit.target_start.max(hit.target_end) as i64;
            let entry = seq_anchor_coverage.entry(hit.query_id).or_insert((q_start, q_end, t_start, t_end));
            entry.0 = entry.0.min(q_start);
            entry.1 = entry.1.max(q_end);
            entry.2 = entry.2.min(t_start);
            entry.3 = entry.3.max(t_end);
        }
    }

    // Pass 2: Process all hits for depth
    for hit in &hits {
        // Filter by sequence inclusion (min_seq_length)
        if min_seq_length > 0 && !seq_included.get(hit.query_id as usize).copied().unwrap_or(false) {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(hit.query_id);

        // Skip self-alignment (same sample)
        if query_sample_id == anchor_sample_id {
            continue;
        }

        let query_start = hit.query_start.min(hit.query_end) as i64;
        let query_end = hit.query_start.max(hit.query_end) as i64;

        discovered_regions.push((hit.query_id, query_start, query_end));

        // Determine anchor coordinates
        let (a_start, a_end) = if hit.target_id == anchor_seq_id {
            let t_start = hit.target_start.min(hit.target_end) as i64;
            let t_end = hit.target_start.max(hit.target_end) as i64;
            (t_start.max(region_start), t_end.min(region_end))
        } else {
            let t_start = hit.target_start.min(hit.target_end) as i64;
            let t_end = hit.target_start.max(hit.target_end) as i64;
            if let Some(&(seq_start, seq_end, anc_start, anc_end)) =
                seq_anchor_coverage.get(&hit.target_id)
            {
                let seq_len = seq_end - seq_start;
                let anc_len = anc_end - anc_start;
                if seq_len > 0 && anc_len > 0 {
                    let frac_s = (t_start - seq_start) as f64 / seq_len as f64;
                    let frac_e = (t_end - seq_start) as f64 / seq_len as f64;
                    let proj_s = anc_start + (frac_s * anc_len as f64) as i64;
                    let proj_e = anc_start + (frac_e * anc_len as f64) as i64;
                    (proj_s.min(proj_e).max(region_start), proj_s.max(proj_e).min(region_end))
                } else {
                    (region_start, region_end)
                }
            } else {
                (region_start, region_end)
            }
        };

        if a_start >= a_end { continue; }

        alignments.push(CompactAlignmentInfo::new(
            query_sample_id, hit.query_id,
            query_start, query_end,
            a_start, a_end,
            hit.is_reverse,
        ));
    }

    // Sweep-line to compute depth intervals
    let seq_intervals = sweep_line_depth(
        &alignments, num_samples, region_start, region_end,
    );

    AnchorRegionResult {
        intervals: seq_intervals,
        discovered_regions,
        anchor_seq_id,
        anchor_sample_id,
    }
}

/// Process a single anchor region: query alignments, compute depth via sweep-line.
/// Returns depth intervals and all discovered regions (for marking as processed).
///
/// Query strategy:
/// - Non-transitive (default): direct query — O(log n), bounded results.
///   impg's bidirectional index means 1-hop already captures both alignment directions.
/// - Transitive (-x): routes to raw BFS (default) or CIGAR BFS (--use-BFS).
///   Caller should pass bounded-size regions (see TRANSITIVE_CHUNK_SIZE).
fn process_anchor_region(
    impg: &impl ImpgIndex,
    config: &DepthConfig,
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
    seq_included: &[bool],
    min_seq_length: i64,
) -> AnchorRegionResult {
    let is_transitive = config.transitive || config.transitive_dfs;

    // Transitive mode: route to raw BFS (default) or CIGAR BFS (--use-BFS)
    if is_transitive {
        if config.use_cigar_bfs {
            return process_anchor_region_transitive_cigar(
                impg, config, compact_lengths, num_samples,
                anchor_seq_id, anchor_sample_id, region_start, region_end,
                seq_included, min_seq_length,
            );
        } else {
            return process_anchor_region_transitive_raw(
                impg, config, compact_lengths, num_samples,
                anchor_seq_id, anchor_sample_id, region_start, region_end,
                seq_included, min_seq_length,
            );
        }
    }

    // Non-transitive mode: direct 1-hop query + sweep-line
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::new();
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    let overlaps = impg.query(
        anchor_seq_id, region_start as i32, region_end as i32,
        false, None, None, false,
    );

    let mut alignments: Vec<CompactAlignmentInfo> = Vec::new();

    // Add self (anchor sample covers the range)
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id, anchor_seq_id,
        region_start, region_end, region_start, region_end, false,
    ));

    for overlap in &overlaps {
        let query_interval = &overlap.0;
        let target_interval = &overlap.2;

        let query_id = query_interval.metadata;

        // Filter by sequence inclusion (min_seq_length)
        if min_seq_length > 0 && !seq_included.get(query_id as usize).copied().unwrap_or(false) {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(query_id);

        if is_self_alignment(query_sample_id, anchor_sample_id, query_id, anchor_seq_id) {
            continue;
        }

        let is_reverse = query_interval.first > query_interval.last;
        let query_start = query_interval.first.min(query_interval.last) as i64;
        let query_end = query_interval.first.max(query_interval.last) as i64;
        let target_start = target_interval.first.min(target_interval.last) as i64;
        let target_end = target_interval.first.max(target_interval.last) as i64;

        // Clip target to range boundaries
        let clipped_target_start = target_start.max(region_start);
        let clipped_target_end = target_end.min(region_end);

        if clipped_target_start >= clipped_target_end {
            continue;
        }

        // Proportionally adjust query coordinates for clipped target
        let target_len = target_end - target_start;
        let query_len = query_end - query_start;
        let ratio = if target_len > 0 { query_len as f64 / target_len as f64 } else { 1.0 };

        let (cq_start, cq_end) = if is_reverse {
            let cqe = query_end - ((clipped_target_start - target_start) as f64 * ratio) as i64;
            let cqs = query_end - ((clipped_target_end - target_start) as f64 * ratio) as i64;
            (cqs.min(cqe), cqs.max(cqe))
        } else {
            let cqs = query_start + ((clipped_target_start - target_start) as f64 * ratio) as i64;
            let cqe = query_start + ((clipped_target_end - target_start) as f64 * ratio) as i64;
            (cqs.min(cqe), cqs.max(cqe))
        };

        alignments.push(CompactAlignmentInfo::new(
            query_sample_id, query_id,
            cq_start, cq_end,
            clipped_target_start, clipped_target_end,
            is_reverse,
        ));

        discovered_regions.push((query_id, query_start, query_end));
    }

    // Sweep-line to compute depth intervals
    let seq_intervals = sweep_line_depth(
        &alignments, num_samples, region_start, region_end,
    );

    AnchorRegionResult {
        intervals: seq_intervals,
        discovered_regions,
        anchor_seq_id,
        anchor_sample_id,
    }
}

/// Process a single anchor region using pre-scanned raw alignment data.
/// No CIGAR projection — uses raw coordinates from interval tree metadata.
/// This is the fast path for non-transitive depth computation.
fn process_anchor_region_raw(
    raw_intervals: &[RawAlignmentInterval],
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
) -> AnchorRegionResult {
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::new();
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    let mut alignments: Vec<CompactAlignmentInfo> = Vec::new();

    // Self alignment (anchor covers itself)
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id, anchor_seq_id,
        region_start, region_end, region_start, region_end, false,
    ));

    // Binary search: only scan intervals where target_start < region_end
    // (alignment_table is sorted by target_start after pre-scan)
    let end_idx = raw_intervals.partition_point(|aln| (aln.target_start as i64) < region_end);
    for aln in &raw_intervals[..end_idx] {
        let target_start = aln.target_start as i64;
        let target_end = aln.target_end as i64;

        // Skip intervals that end before region starts
        if target_end <= region_start {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(aln.query_id);

        // Filter self-alignments (same sample, different sequence)
        if is_self_alignment(query_sample_id, anchor_sample_id, aln.query_id, anchor_seq_id) {
            continue;
        }

        let query_start = aln.query_start as i64;
        let query_end = aln.query_end as i64;

        // Clip target to region
        let clipped_target_start = target_start.max(region_start);
        let clipped_target_end = target_end.min(region_end);

        // Proportional query coordinate clipping (linear interpolation)
        let target_len = target_end - target_start;
        let query_len = query_end - query_start;
        let ratio = if target_len > 0 { query_len as f64 / target_len as f64 } else { 1.0 };

        let (clipped_query_start, clipped_query_end) = if aln.is_reverse {
            let cqe = query_end - ((clipped_target_start - target_start) as f64 * ratio) as i64;
            let cqs = query_end - ((clipped_target_end - target_start) as f64 * ratio) as i64;
            (cqs, cqe)
        } else {
            let cqs = query_start + ((clipped_target_start - target_start) as f64 * ratio) as i64;
            let cqe = query_start + ((clipped_target_end - target_start) as f64 * ratio) as i64;
            (cqs, cqe)
        };

        alignments.push(CompactAlignmentInfo::new(
            query_sample_id, aln.query_id,
            clipped_query_start.min(clipped_query_end),
            clipped_query_start.max(clipped_query_end),
            clipped_target_start, clipped_target_end,
            aln.is_reverse,
        ));

        // Record discovered query region
        discovered_regions.push((aln.query_id, query_start, query_end));
    }

    // Sweep-line to compute depth intervals
    let seq_intervals = sweep_line_depth(
        &alignments, num_samples, region_start, region_end,
    );

    AnchorRegionResult {
        intervals: seq_intervals,
        discovered_regions,
        anchor_seq_id,
        anchor_sample_id,
    }
}

/// Process a single anchor region using query_transitive_bfs/dfs for discovery,
/// then map results to anchor coordinates for the depth sweep-line.
///
/// Uses the same BFS as the query command (guaranteeing identical sample discovery).
/// Anchor coordinate mapping:
/// - Hop 0 results (target on anchor): precise CIGAR-projected coordinates
/// - Hop 1+ results (target on intermediate): projected back via the intermediate's
///   known anchor coverage from hop 0 results (linear interpolation)
/// - Deeper hops with unknown intermediate: full anchor region as fallback
fn process_anchor_region_transitive_cigar(
    impg: &impl ImpgIndex,
    config: &DepthConfig,
    compact_lengths: &CompactSequenceLengths,
    num_samples: usize,
    anchor_seq_id: u32,
    anchor_sample_id: u16,
    region_start: i64,
    region_end: i64,
    seq_included: &[bool],
    min_seq_length: i64,
) -> AnchorRegionResult {
    let mut discovered_regions: Vec<(u32, i64, i64)> = Vec::new();
    discovered_regions.push((anchor_seq_id, region_start, region_end));

    // Use the same CIGAR-precise BFS/DFS as the query command (--use-BFS path)
    let overlaps = if config.transitive_dfs {
        impg.query_transitive_dfs(
            anchor_seq_id, region_start as i32, region_end as i32, None,
            config.max_depth, config.min_transitive_len,
            config.min_distance_between_ranges, None, false, None,
            None, false, None,
        )
    } else {
        impg.query_transitive_bfs(
            anchor_seq_id, region_start as i32, region_end as i32, None,
            config.max_depth, config.min_transitive_len,
            config.min_distance_between_ranges, None, false, None,
            None, false, None,
        )
    };

    let mut alignments: Vec<CompactAlignmentInfo> = Vec::new();

    // Self alignment (anchor covers itself)
    alignments.push(CompactAlignmentInfo::new(
        anchor_sample_id, anchor_seq_id,
        region_start, region_end, region_start, region_end, false,
    ));

    // Pass 1: Build anchor coverage map from hop 0 results (target on anchor)
    let mut seq_anchor_coverage: FxHashMap<u32, (i64, i64, i64, i64)> = FxHashMap::default();
    for overlap in &overlaps {
        let query_interval = &overlap.0;
        let target_interval = &overlap.2;
        if target_interval.metadata == anchor_seq_id {
            let query_id = query_interval.metadata;
            let q_start = query_interval.first.min(query_interval.last) as i64;
            let q_end = query_interval.first.max(query_interval.last) as i64;
            let t_start = target_interval.first.min(target_interval.last) as i64;
            let t_end = target_interval.first.max(target_interval.last) as i64;
            let entry = seq_anchor_coverage.entry(query_id).or_insert((q_start, q_end, t_start, t_end));
            entry.0 = entry.0.min(q_start);
            entry.1 = entry.1.max(q_end);
            entry.2 = entry.2.min(t_start);
            entry.3 = entry.3.max(t_end);
        }
    }

    // Pass 2: Process all results for depth
    for overlap in &overlaps {
        let query_interval = &overlap.0;
        let target_interval = &overlap.2;

        let query_id = query_interval.metadata;

        // Filter by sequence inclusion (min_seq_length)
        if min_seq_length > 0 && !seq_included.get(query_id as usize).copied().unwrap_or(false) {
            continue;
        }

        let query_sample_id = compact_lengths.get_sample_id(query_id);

        // Skip self-alignment (same sample)
        if query_sample_id == anchor_sample_id {
            continue;
        }

        let is_reverse = query_interval.first > query_interval.last;
        let query_start = query_interval.first.min(query_interval.last) as i64;
        let query_end = query_interval.first.max(query_interval.last) as i64;

        discovered_regions.push((query_id, query_start, query_end));

        // Determine anchor coordinates
        let (a_start, a_end) = if target_interval.metadata == anchor_seq_id {
            let t_start = target_interval.first.min(target_interval.last) as i64;
            let t_end = target_interval.first.max(target_interval.last) as i64;
            (t_start.max(region_start), t_end.min(region_end))
        } else {
            let t_start = target_interval.first.min(target_interval.last) as i64;
            let t_end = target_interval.first.max(target_interval.last) as i64;
            if let Some(&(seq_start, seq_end, anc_start, anc_end)) =
                seq_anchor_coverage.get(&target_interval.metadata)
            {
                let seq_len = seq_end - seq_start;
                let anc_len = anc_end - anc_start;
                if seq_len > 0 && anc_len > 0 {
                    let frac_s = (t_start - seq_start) as f64 / seq_len as f64;
                    let frac_e = (t_end - seq_start) as f64 / seq_len as f64;
                    let proj_s = anc_start + (frac_s * anc_len as f64) as i64;
                    let proj_e = anc_start + (frac_e * anc_len as f64) as i64;
                    (proj_s.min(proj_e).max(region_start), proj_s.max(proj_e).min(region_end))
                } else {
                    (region_start, region_end)
                }
            } else {
                (region_start, region_end)
            }
        };

        if a_start >= a_end { continue; }

        alignments.push(CompactAlignmentInfo::new(
            query_sample_id, query_id,
            query_start, query_end,
            a_start, a_end,
            is_reverse,
        ));
    }

    // Augment discovered_regions with raw alignment extents for direct (hop 0) overlaps.
    // The BFS loop above records CIGAR-projected sub-ranges (line 1899) which may leave
    // gaps at chunk boundaries due to indels. Raw extents ensure the full alignment
    // coverage is marked as processed, preventing Phase 2 from re-processing these
    // regions and producing duplicate output (e.g., CHM13 appearing in Phase 2 rows).
    let raw_hop0 = impg.query_raw_overlapping(
        anchor_seq_id, region_start as i32, region_end as i32,
    );
    for aln in &raw_hop0 {
        if min_seq_length > 0 && !seq_included.get(aln.query_id as usize).copied().unwrap_or(false) {
            continue;
        }
        let query_sample_id = compact_lengths.get_sample_id(aln.query_id);
        if is_self_alignment(query_sample_id, anchor_sample_id, aln.query_id, anchor_seq_id) {
            continue;
        }
        discovered_regions.push((aln.query_id, aln.query_start as i64, aln.query_end as i64));
    }

    // Sweep-line to compute depth intervals
    let seq_intervals = sweep_line_depth(
        &alignments, num_samples, region_start, region_end,
    );

    AnchorRegionResult {
        intervals: seq_intervals,
        discovered_regions,
        anchor_seq_id,
        anchor_sample_id,
    }
}

/// Sweep-line algorithm: given alignments, produce depth intervals with sample tracking.
/// Callers must include the anchor sample as an alignment covering [region_start, region_end]
/// so that depth naturally counts all samples (including the anchor itself).
fn sweep_line_depth(
    alignments: &[CompactAlignmentInfo],
    num_samples: usize,
    region_start: i64,
    region_end: i64,
) -> Vec<SparseDepthInterval> {
    // Build sweep-line events
    let mut events: Vec<CompactDepthEvent> = Vec::with_capacity(alignments.len() * 2);
    for (idx, aln) in alignments.iter().enumerate() {
        events.push(CompactDepthEvent {
            position: aln.target_start,
            is_start: true,
            sample_id: aln.sample_id,
            alignment_idx: idx,
        });
        events.push(CompactDepthEvent {
            position: aln.target_end,
            is_start: false,
            sample_id: aln.sample_id,
            alignment_idx: idx,
        });
    }
    events.sort();

    // Sweep-line to compute depth intervals with SPARSE storage
    let mut seq_intervals: Vec<SparseDepthInterval> = Vec::new();
    let mut active_bitmap = SampleBitmap::new(num_samples);
    let mut active_alns: Vec<Vec<usize>> = vec![Vec::new(); num_samples];
    let mut prev_pos: Option<i64> = None;

    for event in events {
        if let Some(prev) = prev_pos {
            if event.position > prev && active_bitmap.depth() > 0 {
                let interval_start = prev;
                let interval_end = event.position;

                // Clip to anchor region
                let clipped_start = interval_start.max(region_start);
                let clipped_end = interval_end.min(region_end);

                if clipped_start < clipped_end {
                    // Build SPARSE sample positions
                    let mut samples: Vec<SamplePosition> = Vec::new();

                    for sample_id in active_bitmap.active_samples() {
                        let alns = &active_alns[sample_id as usize];
                        if let Some(&best_idx) = alns.iter().max_by_key(|&&idx| {
                            let aln = &alignments[idx];
                            let overlap_start = clipped_start.max(aln.target_start);
                            let overlap_end = clipped_end.min(aln.target_end);
                            overlap_end - overlap_start
                        }) {
                            let aln = &alignments[best_idx];
                            let (q_start, q_end) = map_target_to_query_linear(
                                &[],
                                aln.target_start, aln.target_end,
                                aln.query_start, aln.query_end,
                                clipped_start, clipped_end,
                                aln.is_reverse,
                            );
                            samples.push((sample_id, aln.query_id, q_start, q_end));
                        }
                    }

                    if !samples.is_empty() {
                        samples.sort_by_key(|s| s.0);
                        seq_intervals.push(SparseDepthInterval {
                            start: clipped_start,
                            end: clipped_end,
                            samples,
                        });
                    }
                }
            }
        }

        // Update active samples
        if event.is_start {
            active_bitmap.add(event.sample_id);
            active_alns[event.sample_id as usize].push(event.alignment_idx);
        } else {
            active_bitmap.remove(event.sample_id);
            active_alns[event.sample_id as usize].retain(|&idx| idx != event.alignment_idx);
        }

        prev_pos = Some(event.position);
    }

    seq_intervals
}

/// Maximum region size for a single transitive BFS query (5 MB).
/// Prevents BFS from visiting the entire pangenome when starting from a large chromosome.
/// Non-transitive queries don't need this limit (they're O(log n) per lookup).
const TRANSITIVE_CHUNK_SIZE: i64 = 5_000_000;

/// Reference-free global depth computation.
///
/// Algorithm:
/// 1. Sort all sequences by length descending (--ref sample gets priority)
/// 2. For each sequence, find unprocessed regions
/// 3. Query alignments for each unprocessed region (direct or transitive)
/// 4. Compute depth via sweep-line, mark all discovered regions as processed
/// 5. Every base in the pangenome is assigned to exactly one output row
///
/// Query modes:
/// - Default (non-transitive): direct 1-hop queries. Fast, bounded memory.
///   impg's bidirectional index captures both alignment directions in one lookup.
/// - Transitive (-x): BFS with 5MB region chunking to bound memory.
///
/// Key properties:
/// - Every base processed exactly once (no double-counting)
/// - Hub-first ordering: high-connectivity sequences anchor first (auto-detected or via --ref)
/// - --ref: ref-anchored mode, guarantees ref sample's coordinate system for covered regions
/// - --ref-only: ref-only mode, output filtered to ref sample's anchored regions only
pub fn compute_depth_global(
    impg: &impl ImpgIndex,
    config: &DepthConfig,
    separator: &str,
    output_prefix: Option<&str>,
    fai_list: Option<&str>,
    window_size: Option<i64>,
    merge_tolerance: f64,
    ref_sample: Option<&str>,
    ref_only: bool,
    min_seq_length: i64,
    stats_mode: bool,
    stats_combined: bool,
) -> io::Result<()> {
    let is_transitive = config.transitive || config.transitive_dfs;
    let user_window_size = window_size;
    let window_size = window_size.unwrap_or(DEFAULT_WINDOW_SIZE);
    if is_transitive {
        info!(
            "Computing depth (global mode, transitive BFS, {}MB query chunks)",
            TRANSITIVE_CHUNK_SIZE / 1_000_000
        );
    } else {
        info!("Computing depth (global mode, pre-scan + sweep-line)");
    }
    if merge_tolerance > 0.0 {
        info!("Merge tolerance: {:.1}%", merge_tolerance * 100.0);
    }

    // Build compact data structures
    let compact_lengths = CompactSequenceLengths::from_impg(impg, separator);
    let sample_index = compact_lengths.sample_index();
    let num_samples = sample_index.len();
    let num_sequences = impg.seq_index().len();

    debug!(
        "Found {} samples, {} sequences",
        num_samples, num_sequences
    );

    // Build sequence inclusion filter (for --min-seq-length)
    let seq_included: Vec<bool> = (0..num_sequences as u32)
        .map(|id| compact_lengths.get_length(id) >= min_seq_length)
        .collect();
    if min_seq_length > 0 {
        let included = seq_included.iter().filter(|&&v| v).count();
        let excluded = num_sequences - included;
        info!(
            "Sequence length filter (>= {} bp): {} included, {} excluded",
            min_seq_length, included, excluded
        );
    }

    // Resolve ref sample ID
    let ref_sample_id: Option<u16> = if let Some(ref_name) = ref_sample {
        let id = sample_index.get_id(ref_name).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Reference sample '{}' not found in alignment index", ref_name),
            )
        })?;
        if ref_only {
            info!("Ref-only mode: '{}' (output filtered to this sample's anchored regions)", ref_name);
        } else {
            info!("Ref-anchored mode: '{}' (Phase 1 anchor, coordinate system priority)", ref_name);
        }
        Some(id)
    } else {
        if ref_only {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--ref-only requires --ref to be specified",
            ));
        }
        None
    };

    // Load FAI sequences if provided
    let fai_seq_lengths: Option<SequenceLengths> = if let Some(fai_path) = fai_list {
        let fai = SequenceLengths::from_fai_list(fai_path, separator)?;
        debug!("FAI: {} samples, {} sequences", fai.sample_to_seqs.len(), fai.lengths.len());
        Some(fai)
    } else {
        None
    };

    // Pre-scan: compute alignment degrees for all included sequences.
    // Used for automatic hub detection (when --ref is not specified) and
    // for sorting sequences by connectivity (hubs first) in all modes.
    //
    // Disable tree caching during pre-scan: with per-file indexing, 64 threads
    // simultaneously loading all sub-indices and caching every interval tree
    // causes 200+ GB RSS. Trees are only needed transiently here.
    impg.set_tree_cache_enabled(false);
    info!("Pre-scanning alignment degrees...");
    let degrees = compute_alignment_degrees(impg, &compact_lengths, &seq_included);
    impg.clear_sub_index_cache();
    let max_degree = degrees.iter().copied().max().unwrap_or(0);
    let included_count = seq_included.iter().filter(|&&v| v).count();
    info!(
        "Degree scan complete: {} sequences, max degree = {}",
        included_count, max_degree
    );

    // Build sequence processing order (degree descending, length descending, ref sample first)
    let sequence_order = build_sequence_order(impg, &compact_lengths, ref_sample_id, &seq_included, &degrees);
    debug!("Sequence processing order: {} sequences", sequence_order.len());

    // Prepare output: TSV writer for normal mode, stats accumulators for --stats add-on
    let writer: Option<Mutex<BufWriter<Box<dyn Write + Send>>>> = if !stats_mode {
        let inner_writer: Box<dyn Write + Send> = if let Some(prefix) = output_prefix {
            let path = format!("{}.depth.tsv", prefix);
            Box::new(BufWriter::new(std::fs::File::create(&path)?))
        } else {
            Box::new(BufWriter::new(std::io::stdout()))
        };
        let mut w = BufWriter::new(inner_writer);

        // Write header: #id length depth Sample1 Sample2 ...
        write!(w, "#id\tlength\tdepth")?;
        for sample_id in 0..num_samples as u16 {
            if let Some(name) = sample_index.get_name(sample_id) {
                write!(w, "\t{}", name)?;
            }
        }
        writeln!(w)?;
        Some(Mutex::new(w))
    } else {
        None
    };

    // Stats accumulators (only allocated in --stats mode)
    let stats_accumulator: Option<Mutex<DepthStats>> = if stats_mode && !stats_combined {
        Some(Mutex::new(DepthStats::new()))
    } else {
        None
    };
    let stats_combined_acc: Option<Mutex<DepthStatsWithSamples>> = if stats_mode && stats_combined {
        Some(Mutex::new(DepthStatsWithSamples::new()))
    } else {
        None
    };

    let total_sequences = sequence_order.len();

    // Enable tree caching: trees loaded by get_or_load_tree are cached in sub-indices.
    // With jemalloc, this no longer causes RSS fragmentation. Trees are reused across
    // BFS nodes within the same sequence (avoiding repeated disk I/O), then freed
    // when clear_sub_index_cache drops the sub-index Arc (and all its cached trees).
    impg.set_tree_cache_enabled(true);

    let tracker = ConcurrentProcessedTracker::new(num_sequences);

    // =========================================================================
    // Two-phase parallel processing with hub-first guarantee.
    //
    // Phase 1 ensures high-connectivity (hub) sequences are processed first,
    // so their coordinate systems become the anchors for all homologous regions.
    //
    // Hub selection:
    //   --ref specified: ref sample's sequences are the hubs.
    //   --ref not specified: auto-detect hubs via alignment degree pre-scan.
    //     Sequences with degree >= max_degree/2 are classified as hubs.
    //     For star topology (hub degree ~800, leaf degree ~1), this cleanly
    //     separates the central sequence from leaves.
    //
    // Phase 1 parallelism:
    //   Transitive (-x): chunk-level (5MB chunks) for full thread utilization.
    //   Non-transitive: sequence-level (fast per-sequence, acceptable overhead).
    //
    // Phase 2: remaining sequences, most regions already claimed by Phase 1.
    // =========================================================================

    // Partition sequences into Phase 1 (hubs) and Phase 2 (rest)
    let (phase1_seqs, phase2_seqs): (Vec<u32>, Vec<u32>) = if let Some(ref_id) = ref_sample_id {
        // --ref specified: ref sequences go to Phase 1
        sequence_order.iter().partition(|&&seq_id| {
            compact_lengths.get_sample_id(seq_id) == ref_id
        })
    } else {
        // No --ref: auto-detect hubs by alignment degree
        let hub_threshold = (max_degree + 1) / 2; // ceiling of max_degree/2
        if hub_threshold > 1 {
            let (p1, p2): (Vec<u32>, Vec<u32>) = sequence_order.iter().partition(|&&seq_id| {
                degrees.get(seq_id as usize).copied().unwrap_or(0) >= hub_threshold
            });
            if !p1.is_empty() {
                info!(
                    "Auto-detected {} hub sequences (degree >= {}) for Phase 1",
                    p1.len(), hub_threshold
                );
            }
            (p1, p2)
        } else {
            // No clear hubs (max_degree <= 1), skip Phase 1
            (Vec::new(), sequence_order.clone())
        }
    };

    let processed_count = AtomicUsize::new(0);
    let row_counter = AtomicUsize::new(0);
    let intervals_counter = AtomicUsize::new(0);

    // Helper closure: format and write results for a batch of AnchorRegionResults.
    // In stats mode: accumulates depth distribution + intervals into stats accumulators.
    // In normal mode: formats TSV rows into buf for the caller to write.
    let write_results = |region_results: Vec<AnchorRegionResult>,
                         buf: &mut Vec<u8>| -> io::Result<()> {
        // Thread-local accumulators (merged into global at end of batch)
        let mut local_stats: Option<DepthStats> = if stats_accumulator.is_some() {
            Some(DepthStats::new())
        } else {
            None
        };
        let mut local_combined: Option<DepthStatsWithSamples> = if stats_combined_acc.is_some() {
            Some(DepthStatsWithSamples::new())
        } else {
            None
        };

        for result in region_results {
            let seq_name = impg.seq_index().get_name(result.anchor_seq_id).unwrap_or("?");
            let anchor_sample_id = result.anchor_sample_id;

            let should_output = if ref_only {
                ref_sample_id == Some(anchor_sample_id)
            } else {
                true
            };

            if !should_output {
                continue;
            }

            let windowed_intervals = if user_window_size.is_some() {
                split_intervals_by_window(result.intervals, window_size)
            } else {
                result.intervals
            };

            let final_intervals = if merge_tolerance > 0.0 {
                merge_sparse_intervals(windowed_intervals, merge_tolerance)
            } else {
                windowed_intervals
            };

            for interval in &final_intervals {
                if stats_mode {
                    // Stats add-on: accumulate depth distribution and intervals
                    let depth = interval.depth();
                    if let Some(ref mut ls) = local_stats {
                        ls.add_interval(seq_name, interval.start, interval.end, depth);
                    }
                    if let Some(ref mut lc) = local_combined {
                        let sample_names: Vec<String> = interval.samples.iter()
                            .filter_map(|&(sid, _, _, _)| {
                                sample_index.get_name(sid).map(|s| s.to_string())
                            })
                            .collect();
                        lc.add_interval(seq_name, interval.start, interval.end, depth, sample_names);
                    }
                } else {
                    // Normal mode: write TSV row to buf
                    let rid = row_counter.fetch_add(1, Ordering::Relaxed) + 1;

                    write!(buf, "{}\t{}\t{}",
                           rid, interval.end - interval.start, interval.depth())?;

                    for sid in 0..num_samples as u16 {
                        if sid == anchor_sample_id {
                            write!(buf, "\t{}:{}-{}", seq_name, interval.start, interval.end)?;
                        } else if let Some((query_id, start, end)) = interval.get_sample(sid) {
                            let query_name = impg.seq_index().get_name(query_id).unwrap_or("?");
                            write!(buf, "\t{}:{}-{}", query_name, start, end)?;
                        } else {
                            write!(buf, "\tNA")?;
                        }
                    }
                    writeln!(buf)?;

                    intervals_counter.fetch_add(1, Ordering::Relaxed);
                }
            }
        }

        // Merge thread-local stats into global accumulators (one lock per batch)
        if let Some(ls) = local_stats {
            if let Some(ref acc) = stats_accumulator {
                acc.lock().unwrap().merge(&ls);
            }
        }
        if let Some(lc) = local_combined {
            if let Some(ref acc) = stats_combined_acc {
                acc.lock().unwrap().merge(lc);
            }
        }

        Ok(())
    };

    // =========================================================================
    // Phase 1: Process hub sequences first (guaranteed to complete before Phase 2)
    // =========================================================================
    if !phase1_seqs.is_empty() {
        if is_transitive {
            // Transitive mode: chunk-level parallelism for full thread utilization.
            // Split hub sequences into 5MB chunks → ~1200 tasks for 128 threads.
            let phase1_chunks: Vec<(u32, u16, i64, i64)> = phase1_seqs.iter().flat_map(|&seq_id| {
                let seq_len = compact_lengths.get_length(seq_id);
                let sample_id = compact_lengths.get_sample_id(seq_id);
                let mut chunks = Vec::new();
                let mut pos = 0i64;
                while pos < seq_len {
                    let chunk_end = (pos + TRANSITIVE_CHUNK_SIZE).min(seq_len);
                    chunks.push((seq_id, sample_id, pos, chunk_end));
                    pos = chunk_end;
                }
                chunks
            }).collect();

            let num_chunks = phase1_chunks.len();
            info!(
                "Phase 1: {} hub sequences -> {} chunks ({}MB each), chunk-level parallelism",
                phase1_seqs.len(), num_chunks, TRANSITIVE_CHUNK_SIZE / 1_000_000
            );

            let pb_phase1 = ProgressBar::new(num_chunks as u64);
            pb_phase1.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} chunks ({eta}) | Phase 1: hub sequences")
                    .unwrap()
                    .progress_chars("#>-")
            );

            let phase1_count = AtomicUsize::new(0);

            phase1_chunks.par_iter().try_for_each(|&(seq_id, sample_id, chunk_start, chunk_end)| -> io::Result<()> {
                let result = process_anchor_region(
                    impg,
                    config,
                    &compact_lengths, num_samples,
                    seq_id, sample_id, chunk_start, chunk_end,
                    &seq_included, min_seq_length,
                );

                tracker.mark_processed_batch(&result.discovered_regions);

                let mut buf: Vec<u8> = Vec::new();
                write_results(vec![result], &mut buf)?;
                if !buf.is_empty() {
                    if let Some(ref w) = writer {
                        w.lock().unwrap().write_all(&buf)?;
                    }
                }

                impg.clear_sub_index_cache();

                let count = phase1_count.fetch_add(1, Ordering::Relaxed) + 1;
                pb_phase1.set_position(count as u64);

                Ok(())
            })?;

            pb_phase1.finish_and_clear();
        } else {
            // Non-transitive mode: sequence-level parallelism.
            // Hub sequences are few but each is fast (direct query + sweep-line).
            info!(
                "Phase 1: {} hub sequences, sequence-level parallelism",
                phase1_seqs.len()
            );

            let pb_phase1 = ProgressBar::new(phase1_seqs.len() as u64);
            pb_phase1.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} seqs ({eta}) | Phase 1: hub sequences")
                    .unwrap()
                    .progress_chars("#>-")
            );

            let phase1_count = AtomicUsize::new(0);

            phase1_seqs.par_iter().try_for_each(|&seq_id| -> io::Result<()> {
                let seq_len = compact_lengths.get_length(seq_id);
                if seq_len <= 0 {
                    let count = phase1_count.fetch_add(1, Ordering::Relaxed) + 1;
                    pb_phase1.set_position(count as u64);
                    return Ok(());
                }

                let sample_id = compact_lengths.get_sample_id(seq_id);

                let mut raw_alns = impg.query_raw_intervals(seq_id);
                if min_seq_length > 0 {
                    raw_alns.retain(|aln| {
                        seq_included.get(aln.query_id as usize).copied().unwrap_or(false)
                    });
                    raw_alns.shrink_to_fit();
                }
                raw_alns.sort_unstable_by_key(|aln| aln.target_start);

                let result = process_anchor_region_raw(
                    &raw_alns, &compact_lengths, num_samples,
                    seq_id, sample_id, 0, seq_len,
                );
                tracker.mark_processed_batch(&result.discovered_regions);

                let mut buf: Vec<u8> = Vec::new();
                write_results(vec![result], &mut buf)?;
                if !buf.is_empty() {
                    if let Some(ref w) = writer {
                        w.lock().unwrap().write_all(&buf)?;
                    }
                }

                let count = phase1_count.fetch_add(1, Ordering::Relaxed) + 1;
                pb_phase1.set_position(count as u64);

                Ok(())
            })?;

            pb_phase1.finish_and_clear();
        }

        // Mark all Phase 1 sequences as fully processed in the tracker
        for &seq_id in &phase1_seqs {
            let seq_len = compact_lengths.get_length(seq_id);
            tracker.mark_processed(seq_id, 0, seq_len);
        }

        processed_count.store(phase1_seqs.len(), Ordering::Relaxed);
        info!(
            "Phase 1 complete: {} hub sequences processed",
            phase1_seqs.len()
        );
    }

    // =========================================================================
    // Phase 2: Process remaining sequences (sequence-level parallelism)
    // =========================================================================
    let phase2_label = if phase1_seqs.is_empty() { "" } else { "Phase 2: " };
    let pb_depth = ProgressBar::new(total_sequences as u64);
    pb_depth.set_style(
        ProgressStyle::default_bar()
            .template(&format!(
                "{{spinner:.green}} [{{elapsed_precise}}] [{{wide_bar:.cyan/blue}}] {{pos}}/{{len}} seqs ({{eta}}) | {}remaining sequences",
                phase2_label
            ))
            .unwrap()
            .progress_chars("#>-")
    );
    pb_depth.set_position(processed_count.load(Ordering::Relaxed) as u64);

    phase2_seqs.par_iter().try_for_each(|&seq_id| -> io::Result<()> {
        let seq_len = compact_lengths.get_length(seq_id);
        if seq_len <= 0 {
            let count = processed_count.fetch_add(1, Ordering::Relaxed) + 1;
            pb_depth.set_position(count as u64);
            return Ok(());
        }

        let unprocessed = tracker.get_unprocessed(seq_id, 0, seq_len);
        if unprocessed.is_empty() {
            let count = processed_count.fetch_add(1, Ordering::Relaxed) + 1;
            pb_depth.set_position(count as u64);
            return Ok(());
        }

        let sample_id = compact_lengths.get_sample_id(seq_id);

        // Non-transitive: load raw intervals once per sequence (not per gap)
        let raw_alns_cached = if !is_transitive {
            let mut raw_alns = impg.query_raw_intervals(seq_id);
            if min_seq_length > 0 {
                raw_alns.retain(|aln| {
                    seq_included.get(aln.query_id as usize).copied().unwrap_or(false)
                });
                raw_alns.shrink_to_fit();
            }
            raw_alns.sort_unstable_by_key(|aln| aln.target_start);
            Some(raw_alns)
        } else {
            None
        };

        for (region_start, region_end) in unprocessed {
            let region_results = if is_transitive {
                let gap_len = region_end - region_start;

                // Fast path: gaps smaller than min_transitive_len can't trigger
                // transitive hops, so BFS setup (visited_ranges allocation for all
                // sequences, sub-index loading) is wasted. Use raw intervals instead.
                if gap_len < config.min_transitive_len as i64 {
                    let mut raw_alns = impg.query_raw_overlapping(
                        seq_id, region_start as i32, region_end as i32,
                    );
                    if min_seq_length > 0 {
                        raw_alns.retain(|aln| {
                            seq_included.get(aln.query_id as usize).copied().unwrap_or(false)
                        });
                    }
                    raw_alns.sort_unstable_by_key(|aln| aln.target_start);
                    let result = process_anchor_region_raw(
                        &raw_alns, &compact_lengths, num_samples,
                        seq_id, sample_id, region_start, region_end,
                    );
                    tracker.mark_processed_batch(&result.discovered_regions);
                    vec![result]
                } else if gap_len > TRANSITIVE_CHUNK_SIZE {
                    let mut results = Vec::new();
                    let mut pos = region_start;
                    while pos < region_end {
                        let chunk_end = (pos + TRANSITIVE_CHUNK_SIZE).min(region_end);
                        let result = process_anchor_region(
                            impg,
                            config,
                            &compact_lengths, num_samples,
                            seq_id, sample_id, pos, chunk_end,
                            &seq_included, min_seq_length,
                        );
                        tracker.mark_processed_batch(&result.discovered_regions);
                        results.push(result);
                        pos = chunk_end;
                    }
                    results
                } else {
                    let result = process_anchor_region(
                        impg,
                        config,
                        &compact_lengths, num_samples,
                        seq_id, sample_id, region_start, region_end,
                        &seq_included, min_seq_length,
                    );
                    tracker.mark_processed_batch(&result.discovered_regions);
                    vec![result]
                }
            } else {
                let result = process_anchor_region_raw(
                    raw_alns_cached.as_ref().unwrap(), &compact_lengths, num_samples,
                    seq_id, sample_id, region_start, region_end,
                );
                tracker.mark_processed_batch(&result.discovered_regions);
                vec![result]
            };

            // Format output into thread-local buffer, then write atomically
            let mut buf: Vec<u8> = Vec::new();
            write_results(region_results, &mut buf)?;
            if !buf.is_empty() {
                if let Some(ref w) = writer {
                    w.lock().unwrap().write_all(&buf)?;
                }
            }
        }

        // Phase 2: do NOT clear sub-index cache after each sequence.
        // Unlike Phase 1 (which processes large hub regions and can accumulate many trees),
        // Phase 2 sequences are mostly small/partially-processed. Keeping sub-indices loaded
        // across sequences avoids repeated disk deserialization — the total memory is bounded
        // by num_alignment_files × sub_index_size regardless of how many sequences we process.
        let count = processed_count.fetch_add(1, Ordering::Relaxed) + 1;
        pb_depth.set_position(count as u64);

        Ok(())
    })?;

    // Clear sub-index cache once after all Phase 2 processing
    if is_transitive {
        impg.clear_sub_index_cache();
    }

    pb_depth.finish_and_clear();

    // Handle FAI sequences not in alignment index (depth=1 for unaligned sequences)
    if let Some(ref fai) = fai_seq_lengths {
        let indexed_seqs: FxHashSet<&str> = (0..impg.seq_index().len() as u32)
            .filter_map(|id| impg.seq_index().get_name(id))
            .collect();

        for (seq_name, &seq_len) in &fai.lengths {
            if indexed_seqs.contains(seq_name.as_str()) {
                continue;
            }

            let sample_name = extract_sample(seq_name, separator);

            // ref_only filter for FAI sequences
            if ref_only {
                if let Some(ref_name) = ref_sample {
                    if sample_name != ref_name {
                        continue;
                    }
                }
            }

            if stats_mode {
                // Add depth=1 to stats accumulators
                if let Some(ref acc) = stats_accumulator {
                    acc.lock().unwrap().add_interval(seq_name, 0, seq_len, 1);
                }
                if let Some(ref acc) = stats_combined_acc {
                    acc.lock().unwrap().add_interval(
                        seq_name, 0, seq_len, 1,
                        vec![sample_name.to_string()],
                    );
                }
            } else if let Some(ref w) = writer {
                let rid = row_counter.fetch_add(1, Ordering::Relaxed) + 1;
                let mut w = w.lock().unwrap();

                // Output as depth=1 for entire sequence
                write!(w, "{}\t{}\t1", rid, seq_len)?;
                for sample_id in 0..num_samples as u16 {
                    let current_sample = sample_index.get_name(sample_id).unwrap_or("");
                    if current_sample == sample_name {
                        write!(w, "\t{}:0-{}", seq_name, seq_len)?;
                    } else {
                        write!(w, "\tNA")?;
                    }
                }
                writeln!(w)?;
            }
        }
    }

    // Finalize output
    if stats_mode {
        let prefix = output_prefix.unwrap_or("depth_stats");

        if let Some(acc) = stats_accumulator {
            let stats = acc.into_inner().unwrap();

            // Write summary
            let summary_path = format!("{}.summary.txt", prefix);
            let mut summary_writer = BufWriter::new(std::fs::File::create(&summary_path)?);
            stats.write_summary(&mut summary_writer)?;
            summary_writer.flush()?;
            info!("Wrote summary to {}", summary_path);

            // Print summary to stdout
            stats.write_summary(&mut std::io::stdout())?;

            // Write per-depth BED files
            stats.write_depth_bed_files(prefix)?;

            info!(
                "Stats complete: {} total bases, max depth = {}",
                stats.total_bases,
                stats.max_depth()
            );
        }

        if let Some(acc) = stats_combined_acc {
            let mut stats = acc.into_inner().unwrap();

            // Write summary
            let summary_path = format!("{}.summary.txt", prefix);
            let mut summary_writer = BufWriter::new(std::fs::File::create(&summary_path)?);
            stats.write_summary(&mut summary_writer)?;
            summary_writer.flush()?;
            info!("Wrote summary to {}", summary_path);

            // Print summary to stdout
            stats.write_summary(&mut std::io::stdout())?;

            // Write combined output file (with optional merging)
            stats.write_combined_output(prefix, merge_tolerance)?;

            info!(
                "Stats complete: {} total bases, {} intervals, max depth = {}",
                stats.total_bases,
                stats.intervals.len(),
                stats.max_depth()
            );
        }
    } else {
        // Normal mode: flush TSV writer
        if let Some(w) = writer {
            let mut w = w.into_inner().unwrap();
            w.flush()?;
        }

        let total_intervals = row_counter.load(Ordering::Relaxed);
        info!(
            "Depth complete: {} sequences, {} intervals output",
            total_sequences, total_intervals
        );
    }

    Ok(())
}

// ============================================================================
// Region query mode and sweep-line helpers
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
        let q_end = aln_query_end - (start_offset as f64 * ratio).round() as i64;
        let q_start = aln_query_end - (end_offset as f64 * ratio).round() as i64;
        (q_start.min(q_end), q_start.max(q_end))
    } else {
        let q_start = aln_query_start + (start_offset as f64 * ratio).round() as i64;
        let q_end = aln_query_start + (end_offset as f64 * ratio).round() as i64;
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
    impg: &impl ImpgIndex,
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

    let target_id = impg.seq_index().get_id(target_seq).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!("Sequence '{}' not found in index", target_seq),
        )
    })?;

    let target_sample = extract_sample(target_seq, separator);

    // Collect all alignments for this region
    let mut alignments: Vec<AlignmentInfoMulti> = Vec::new();

    let is_transitive = config.transitive || config.transitive_dfs;

    // Forward direction: where target_seq is TARGET
    if is_transitive && config.use_cigar_bfs {
        // CIGAR-precise BFS (--use-BFS): identical to query command's BFS
        let overlaps = if config.transitive_dfs {
            impg.query_transitive_dfs(
                target_id, target_start, target_end, None,
                config.max_depth, config.min_transitive_len,
                config.min_distance_between_ranges, None, false, None,
                sequence_index, false, None,
            )
        } else {
            impg.query_transitive_bfs(
                target_id, target_start, target_end, None,
                config.max_depth, config.min_transitive_len,
                config.min_distance_between_ranges, None, false, None,
                sequence_index, false, None,
            )
        };

        let region_start = target_start as i64;
        let region_end = target_end as i64;

        // Pass 1: Build a map of each sequence's anchor coverage from hop 0 results
        let mut seq_anchor_coverage: FxHashMap<u32, (i64, i64, i64, i64)> = FxHashMap::default();
        for overlap in &overlaps {
            let query_interval = &overlap.0;
            let target_interval = &overlap.2;
            if target_interval.metadata == target_id {
                let query_id = query_interval.metadata;
                let q_start = query_interval.first.min(query_interval.last) as i64;
                let q_end = query_interval.first.max(query_interval.last) as i64;
                let t_start = target_interval.first.min(target_interval.last) as i64;
                let t_end = target_interval.first.max(target_interval.last) as i64;
                let entry = seq_anchor_coverage.entry(query_id).or_insert((q_start, q_end, t_start, t_end));
                entry.0 = entry.0.min(q_start);
                entry.1 = entry.1.max(q_end);
                entry.2 = entry.2.min(t_start);
                entry.3 = entry.3.max(t_end);
            }
        }

        // Pass 2: Process all results for depth
        for overlap in &overlaps {
            let query_interval = &overlap.0;
            let target_interval = &overlap.2;

            let query_id = query_interval.metadata;
            let query_name = match impg.seq_index().get_name(query_id) {
                Some(name) => name,
                None => continue,
            };
            let sample = extract_sample(query_name, separator);

            if sample == target_sample {
                continue;
            }

            if let Some(filter) = sample_filter {
                if !filter.includes(&sample) {
                    continue;
                }
            }

            let is_reverse = query_interval.first > query_interval.last;
            let q_start = query_interval.first.min(query_interval.last) as i64;
            let q_end = query_interval.first.max(query_interval.last) as i64;

            let (a_start, a_end) = if target_interval.metadata == target_id {
                let t_start = target_interval.first.min(target_interval.last) as i64;
                let t_end = target_interval.first.max(target_interval.last) as i64;
                (t_start.max(region_start), t_end.min(region_end))
            } else {
                let t_start = target_interval.first.min(target_interval.last) as i64;
                let t_end = target_interval.first.max(target_interval.last) as i64;
                if let Some(&(seq_start, seq_end, anc_start, anc_end)) =
                    seq_anchor_coverage.get(&target_interval.metadata)
                {
                    let seq_len = seq_end - seq_start;
                    let anc_len = anc_end - anc_start;
                    if seq_len > 0 && anc_len > 0 {
                        let frac_s = (t_start - seq_start) as f64 / seq_len as f64;
                        let frac_e = (t_end - seq_start) as f64 / seq_len as f64;
                        let proj_s = anc_start + (frac_s * anc_len as f64) as i64;
                        let proj_e = anc_start + (frac_e * anc_len as f64) as i64;
                        (proj_s.min(proj_e).max(region_start), proj_s.max(proj_e).min(region_end))
                    } else {
                        (region_start, region_end)
                    }
                } else {
                    (region_start, region_end)
                }
            };

            if a_start >= a_end { continue; }

            alignments.push(AlignmentInfoMulti {
                sample,
                query_name: query_name.to_string(),
                query_start: q_start,
                query_end: q_end,
                target_start: a_start,
                target_end: a_end,
                is_reverse,
            });
        }
    } else if is_transitive {
        // Default transitive: raw-interval BFS with linear interpolation
        let hits = depth_transitive_bfs(
            impg,
            target_id,
            target_start,
            target_end,
            config.max_depth,
            config.min_transitive_len,
            config.min_distance_between_ranges,
            config.transitive_dfs,
        );

        let region_start = target_start as i64;
        let region_end = target_end as i64;

        // Pass 1: Build anchor coverage map from hop 0 results
        let mut seq_anchor_coverage: FxHashMap<u32, (i64, i64, i64, i64)> = FxHashMap::default();
        for hit in &hits {
            if hit.target_id == target_id {
                let q_start = hit.query_start.min(hit.query_end) as i64;
                let q_end = hit.query_start.max(hit.query_end) as i64;
                let t_start = hit.target_start.min(hit.target_end) as i64;
                let t_end = hit.target_start.max(hit.target_end) as i64;
                let entry = seq_anchor_coverage.entry(hit.query_id).or_insert((q_start, q_end, t_start, t_end));
                entry.0 = entry.0.min(q_start);
                entry.1 = entry.1.max(q_end);
                entry.2 = entry.2.min(t_start);
                entry.3 = entry.3.max(t_end);
            }
        }

        // Pass 2: Process all hits for depth
        for hit in &hits {
            let query_name = match impg.seq_index().get_name(hit.query_id) {
                Some(name) => name,
                None => continue,
            };
            let sample = extract_sample(query_name, separator);

            if sample == target_sample {
                continue;
            }

            if let Some(filter) = sample_filter {
                if !filter.includes(&sample) {
                    continue;
                }
            }

            let q_start = hit.query_start.min(hit.query_end) as i64;
            let q_end = hit.query_start.max(hit.query_end) as i64;

            let (a_start, a_end) = if hit.target_id == target_id {
                let t_start = hit.target_start.min(hit.target_end) as i64;
                let t_end = hit.target_start.max(hit.target_end) as i64;
                (t_start.max(region_start), t_end.min(region_end))
            } else {
                let t_start = hit.target_start.min(hit.target_end) as i64;
                let t_end = hit.target_start.max(hit.target_end) as i64;
                if let Some(&(seq_start, seq_end, anc_start, anc_end)) =
                    seq_anchor_coverage.get(&hit.target_id)
                {
                    let seq_len = seq_end - seq_start;
                    let anc_len = anc_end - anc_start;
                    if seq_len > 0 && anc_len > 0 {
                        let frac_s = (t_start - seq_start) as f64 / seq_len as f64;
                        let frac_e = (t_end - seq_start) as f64 / seq_len as f64;
                        let proj_s = anc_start + (frac_s * anc_len as f64) as i64;
                        let proj_e = anc_start + (frac_e * anc_len as f64) as i64;
                        (proj_s.min(proj_e).max(region_start), proj_s.max(proj_e).min(region_end))
                    } else {
                        (region_start, region_end)
                    }
                } else {
                    (region_start, region_end)
                }
            };

            if a_start >= a_end { continue; }

            alignments.push(AlignmentInfoMulti {
                sample,
                query_name: query_name.to_string(),
                query_start: q_start,
                query_end: q_end,
                target_start: a_start,
                target_end: a_end,
                is_reverse: hit.is_reverse,
            });
        }
    } else if config.use_cigar_bfs {
        // Non-transitive with --use-BFS: CIGAR-precise query
        let overlaps = impg.query(
            target_id, target_start, target_end,
            false, None, sequence_index, false,
        );

        for overlap in &overlaps {
            let query_interval = &overlap.0;
            let target_interval = &overlap.2;

            let query_name = match impg.seq_index().get_name(query_interval.metadata) {
                Some(name) => name,
                None => continue,
            };
            let sample = extract_sample(query_name, separator);

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
    } else {
        // Default non-transitive: raw intervals + linear interpolation
        let mut raw_alns = impg.query_raw_overlapping(target_id, target_start, target_end);
        raw_alns.sort_unstable_by_key(|a| a.target_start);

        let region_start = target_start as i64;
        let region_end = target_end as i64;

        for aln in &raw_alns {
            let query_name = match impg.seq_index().get_name(aln.query_id) {
                Some(name) => name,
                None => continue,
            };
            let sample = extract_sample(query_name, separator);

            if let Some(filter) = sample_filter {
                if !filter.includes(&sample) {
                    continue;
                }
            }

            let aln_target_start = aln.target_start as i64;
            let aln_target_end = aln.target_end as i64;
            let aln_query_start = aln.query_start as i64;
            let aln_query_end = aln.query_end as i64;

            // Clip target to region
            let clipped_target_start = aln_target_start.max(region_start);
            let clipped_target_end = aln_target_end.min(region_end);

            // Proportional query coordinate clipping (linear interpolation)
            let target_len = aln_target_end - aln_target_start;
            let query_len = aln_query_end - aln_query_start;
            let ratio = if target_len > 0 { query_len as f64 / target_len as f64 } else { 1.0 };

            let (clipped_query_start, clipped_query_end) = if aln.is_reverse {
                let cqe = aln_query_end - ((clipped_target_start - aln_target_start) as f64 * ratio) as i64;
                let cqs = aln_query_end - ((clipped_target_end - aln_target_start) as f64 * ratio) as i64;
                (cqs, cqe)
            } else {
                let cqs = aln_query_start + ((clipped_target_start - aln_target_start) as f64 * ratio) as i64;
                let cqe = aln_query_start + ((clipped_target_end - aln_target_start) as f64 * ratio) as i64;
                (cqs, cqe)
            };

            alignments.push(AlignmentInfoMulti {
                sample,
                query_name: query_name.to_string(),
                query_start: clipped_query_start.min(clipped_query_end),
                query_end: clipped_query_start.max(clipped_query_end),
                target_start: clipped_target_start,
                target_end: clipped_target_end,
                is_reverse: aln.is_reverse,
            });
        }
    }

    // Reverse direction: where target_seq is QUERY
    let reverse_alignments = impg.query_reverse_for_depth(target_id);
    for (ref_start, ref_end, other_id) in reverse_alignments {
        // Check if overlaps with query region
        if ref_end <= target_start || ref_start >= target_end {
            continue;
        }

        if let Some(other_name) = impg.seq_index().get_name(other_id) {
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
    let seq_len = impg.seq_index().get_len_from_id(target_id).unwrap_or(0) as i64;
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
