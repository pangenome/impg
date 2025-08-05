use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use bitvec::{bitvec, prelude::BitVec};
use gzp::{
    deflate::{Bgzf, Gzip}, // Both Gzip and Bgzf are in deflate module
    par::compress::{ParCompress, ParCompressBuilder},
    Compression,
};
use handlegraph::handle::{Handle, NodeId};
use log::{debug, error, info, warn};
use niffler::compression::Format;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use tempfile::NamedTempFile;
use zstd::stream::Encoder as ZstdEncoder;

// use std::process::Command;
// #[cfg(not(debug_assertions))]
// fn log_memory_usage(stage: &str) {
//     let output = Command::new("ps")
//         .args(&["-o", "rss=", "-p", &std::process::id().to_string()])
//         .output()
//         .expect("Failed to execute ps command");
//     let memory_kb = String::from_utf8_lossy(&output.stdout)
//         .trim()
//         .parse::<u64>()
//         .unwrap_or(0);
//     let memory_mb = memory_kb as f64 / 1024.0;
//     info!("Memory usage at {}: {:.2} MB", stage, memory_mb);
// }

fn get_compression_format(compress_arg: &str, output_path: &str) -> Format {
    match compress_arg.to_lowercase().as_str() {
        "none" => Format::No,
        "gzip" | "gz" => Format::Gzip,
        "bgzip" | "bgz" => Format::Bzip,
        "zstd" | "zst" => Format::Zstd,
        "auto" => {
            // Auto-detect based on file extension
            if output_path.ends_with(".gz") {
                Format::Gzip
            } else if output_path.ends_with(".bgz") {
                Format::Bzip
            } else if output_path.ends_with(".zst") {
                Format::Zstd
            } else {
                Format::No
            }
        }
        _ => {
            warn!(
                "Unsupported compression format '{compress_arg}', using none"
            );
            Format::No
        }
    }
}

// Compact edge representation using bit-packed orientations
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct CompactEdge {
    // Use top 2 bits for orientations, rest for node IDs
    from: u64, // bit 63: orientation, bits 0-62: node ID
    to: u64,   // bit 63: orientation, bits 0-62: node ID
}

impl CompactEdge {
    const ORIENT_MASK: u64 = 1u64 << 63;
    const ID_MASK: u64 = !Self::ORIENT_MASK;

    fn new(from_id: u64, from_rev: bool, to_id: u64, to_rev: bool) -> Self {
        let from = from_id | (if from_rev { Self::ORIENT_MASK } else { 0 });
        let to = to_id | (if to_rev { Self::ORIENT_MASK } else { 0 });
        CompactEdge { from, to }
    }

    fn from_id(&self) -> u64 {
        self.from & Self::ID_MASK
    }
    fn from_rev(&self) -> bool {
        (self.from & Self::ORIENT_MASK) != 0
    }
    fn to_id(&self) -> u64 {
        self.to & Self::ID_MASK
    }
    fn to_rev(&self) -> bool {
        (self.to & Self::ORIENT_MASK) != 0
    }
}

// Sequence storage with memory mapping
struct SequenceStore {
    sequences_file: File,
    offsets: Vec<(u64, u32)>, // (offset, length) for each node
}

impl SequenceStore {
    fn new(temp_dir: Option<&str>) -> io::Result<Self> {
        let temp_file = if let Some(dir) = temp_dir {
            NamedTempFile::new_in(dir)?
        } else {
            NamedTempFile::new()?
        };

        let sequences_file = temp_file.into_file();

        Ok(SequenceStore {
            sequences_file,
            offsets: Vec::new(),
        })
    }

    fn add_sequence(&mut self, seq: &[u8]) -> io::Result<usize> {
        let offset = self.sequences_file.seek(SeekFrom::End(0))?;
        self.sequences_file.write_all(seq)?;

        let idx = self.offsets.len();
        self.offsets.push((offset, seq.len() as u32));
        Ok(idx)
    }

    fn get_sequence(&self, idx: usize) -> io::Result<Vec<u8>> {
        if idx >= self.offsets.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid sequence index",
            ));
        }

        let (offset, length) = self.offsets[idx];
        let mut buffer = vec![0u8; length as usize];

        // Use pread for thread-safe reading without modifying file cursor
        use std::os::unix::fs::FileExt;
        self.sequences_file.read_exact_at(&mut buffer, offset)?;

        Ok(buffer)
    }

    // Get sequence length without reading the actual sequence data
    fn get_sequence_length(&self, idx: usize) -> io::Result<usize> {
        if idx < self.offsets.len() {
            let (_, length) = self.offsets[idx];
            Ok(length as usize)
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Sequence index {} out of bounds", idx),
            ))
        }
    }
}

// Simplified graph structure
struct CompactGraph {
    node_count: u64,
    edges: FxHashSet<CompactEdge>,
    sequence_store: SequenceStore,
}

impl CompactGraph {
    fn new(temp_dir: Option<&str>) -> io::Result<Self> {
        Ok(CompactGraph {
            node_count: 0,
            edges: FxHashSet::default(),
            sequence_store: SequenceStore::new(temp_dir)?,
        })
    }

    fn add_node(&mut self, seq: &[u8]) -> io::Result<u64> {
        self.sequence_store.add_sequence(seq)?;
        self.node_count += 1;
        Ok(self.node_count) // Return the actual node ID (1-based as per handlegraph convention)
    }

    fn add_edge(&mut self, from_id: u64, from_rev: bool, to_id: u64, to_rev: bool) {
        self.edges
            .insert(CompactEdge::new(from_id, from_rev, to_id, to_rev));
    }

    fn has_edge(&self, from_id: u64, from_rev: bool, to_id: u64, to_rev: bool) -> bool {
        // Check for the edge in both directions since handlegraph edges are bidirectional
        let edge1 = CompactEdge::new(from_id, from_rev, to_id, to_rev);
        let edge2 = CompactEdge::new(to_id, !to_rev, from_id, !from_rev);
        self.edges.contains(&edge1) || self.edges.contains(&edge2)
    }

    fn get_sequence(&self, handle: Handle) -> io::Result<Vec<u8>> {
        let seq = self
            .sequence_store
            .get_sequence((u64::from(handle.id()) - 1) as usize)?;
        Ok(if handle.is_reverse() {
            crate::graph::reverse_complement(&seq)
        } else {
            seq
        })
    }
}

#[derive(Debug, Clone)]
struct RangeInfo {
    start: usize,
    end: usize,
    //gfa_id: usize,      // GFA file ID this range belongs to
    steps: Vec<Handle>, // Path steps for this range
}

impl RangeInfo {
    /// Returns true if this range is immediately followed by another range
    /// with no gap between them
    #[inline]
    fn is_contiguous_with(&self, other: &Self) -> bool {
        self.end == other.start
    }

    /// Returns true if this range overlaps with another range
    /// Two ranges overlap if one starts before the other ends
    #[inline]
    fn overlaps_with(&self, other: &Self) -> bool {
        self.start < other.end && other.start < self.end
    }
}

/// Lace pangenome graphs together
pub fn run_gfa_lace(
    gfa_files: Option<Vec<String>>,
    gfa_list: Option<String>,
    output: &str,
    compress: &str,
    fill_gaps: u8,
    skip_validation: bool,
    temp_dir: Option<String>,
    sequence_index: Option<&UnifiedSequenceIndex>,
    verbose: u8,
) -> io::Result<()> {
    // Determine compression format
    let compression_format = get_compression_format(compress, output);

    // Resolve GFA files
    let gfa_files = resolve_gfa_files(gfa_files, gfa_list)?;
    if gfa_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No GFA files specified",
        ));
    }

    // Check if sequence index is available for gap filling
    if sequence_index.is_none() && fill_gaps > 0 {
        warn!("Gap filling requested but no sequence files provided");
    }

    // Validate temp directory
    if let Some(ref temp_path) = temp_dir {
        if !Path::new(temp_path).exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("Temporary directory '{temp_path}' does not exist"),
            ));
        }
    }

    // log_memory_usage("start");

    // Create a single combined graph without paths and a map of path key to ranges
    info!("Collecting metadata from {} GFA files", gfa_files.len());
    let (combined_graph, mut path_key_ranges) = read_gfa_files(&gfa_files, temp_dir.as_deref(), skip_validation)?;

    // log_memory_usage("after_reading_files");

    // Sort and deduplicate path ranges in parallel
    info!(
        "Sorting and deduplicating {} path ranges from {} path keys",
        path_key_ranges
            .values()
            .map(|ranges| ranges.len())
            .sum::<usize>(),
        path_key_ranges.len()
    );
    path_key_ranges
        .par_iter_mut()
        .for_each(|(_path_key, ranges)| {
            sort_and_filter_ranges(ranges);
        });

    // Process different path keys in parallel with minimal locking
    info!(
        "Trimming overlaps and linking contiguous ranges for {} path ranges from {} path keys",
        path_key_ranges.len(),
        path_key_ranges.len()
    );

    // Wrap graph in Arc<Mutex> for thread-safe access
    let graph_mutex = Arc::new(Mutex::new(combined_graph));

    path_key_ranges
        .par_iter_mut()
        .for_each(|(_path_key, ranges)| {
            trim_range_overlaps(ranges, &graph_mutex);
            link_contiguous_ranges(ranges, &graph_mutex);
        });

    // Unwrap the graph from Arc<Mutex>
    let mut combined_graph = Arc::try_unwrap(graph_mutex)
        .ok()
        .expect("Failed to unwrap graph mutex")
        .into_inner()
        .unwrap();

    info!(
        "Created {} nodes and {} edges",
        combined_graph.node_count,
        combined_graph.edges.len()
    );

    // log_memory_usage("before_writing");

    // Write the combined graph to output
    match write_graph_to_gfa(
        &mut combined_graph,
        &path_key_ranges,
        output,
        compression_format,
        fill_gaps,
        sequence_index,
        verbose > 1,
    ) {
        Ok(_) => info!(
            "Successfully wrote the combined graph to {output} ({compression_format:?} format)"
        ),
        Err(e) => error!("Error writing the GFA file: {e}"),
    }

    // log_memory_usage("end");

    Ok(())
}

/// Resolve GFA files from either --gfa-files or --gfa-list
fn resolve_gfa_files(
    gfa_files: Option<Vec<String>>,
    gfa_list: Option<String>,
) -> io::Result<Vec<String>> {
    match (gfa_files, gfa_list) {
        (Some(files), None) => {
            // Validate all files exist
            for file in &files {
                if !Path::new(file).exists() {
                    return Err(io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("GFA file '{file}' not found"),
                    ));
                }
            }
            Ok(files)
        }
        (None, Some(list_file)) => {
            let file = File::open(&list_file)?;
            let reader = BufReader::new(file);
            let mut files = Vec::new();

            for line in reader.lines() {
                let line = line?;
                let trimmed = line.trim();
                if !trimmed.is_empty() && !trimmed.starts_with('#') {
                    if !Path::new(trimmed).exists() {
                        return Err(io::Error::new(
                            io::ErrorKind::NotFound,
                            format!("GFA file '{trimmed}' not found"),
                        ));
                    }
                    files.push(trimmed.to_string());
                }
            }

            if files.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("No valid GFA files found in list file: {list_file}"),
                ));
            }

            Ok(files)
        }
        (None, None) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Either --gfa-files or --gfa-list must be provided",
        )),
        (Some(_), Some(_)) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Cannot specify both --gfa-files and --gfa-list",
        )),
    }
}

fn read_gfa_files(
    gfa_list: &[String],
    temp_dir: Option<&str>,
    skip_validation: bool,
) -> io::Result<(CompactGraph, FxHashMap<String, Vec<RangeInfo>>)> {
    let combined_graph = Arc::new(Mutex::new(CompactGraph::new(temp_dir)?));
    let path_key_ranges: Arc<Mutex<FxHashMap<String, Vec<RangeInfo>>>> =
        Arc::new(Mutex::new(FxHashMap::default()));

    let num_path_ranges = AtomicUsize::new(0);
    let num_path_range_steps = AtomicUsize::new(0);

    gfa_list
        .par_iter()
        .enumerate()
        .for_each(|(gfa_id, gfa_path)| {
            if let Ok(reader) = get_gfa_reader(gfa_path) {
                let mut id_translation: FxHashMap<NodeId, NodeId> = FxHashMap::default();
                let mut temp_edges: Vec<CompactEdge> = Vec::new();

                let mut node_count = 0usize;
                let mut edge_count = 0usize;

                for line in reader.lines().map_while(Result::ok) {
                    let line = line.trim();
                    if line.is_empty() || line.starts_with('#') {
                        continue;
                    }

                    let fields: Vec<&str> = line.split('\t').collect();
                    match fields[0] {
                        "S" => {
                            if fields.len() < 3 {
                                error!("Invalid S line: {line}");
                                std::process::exit(1);
                            }
                            let node_id: u64 = match fields[1].parse() {
                                Ok(id) => id,
                                Err(_) => {
                                    error!("Invalid node ID {} in line {}", fields[1], line);
                                    std::process::exit(1);
                                }
                            };
                            let sequence = fields[2].as_bytes();

                            let new_node_id = {
                                let mut graph = combined_graph.lock().unwrap();
                                graph.add_node(sequence).unwrap()
                            };
                            id_translation.insert(NodeId::from(node_id), NodeId::from(new_node_id));
                            node_count += 1;
                        }
                        "L" => {
                            if fields.len() < 6 {
                                error!("Invalid L line: {line}");
                                std::process::exit(1);
                            }

                            let from_id: u64 = match fields[1].parse() {
                                Ok(id) => id,
                                Err(_) => {
                                    error!("Invalid from node ID {} in line {}", fields[1], line);
                                    std::process::exit(1);
                                }
                            };
                            let from_rev = fields[2] == "-";
                            let to_id: u64 = match fields[3].parse() {
                                Ok(id) => id,
                                Err(_) => {
                                    error!("Invalid to node ID {} in line {}", fields[3], line);
                                    std::process::exit(1);
                                }
                            };
                            let to_rev = fields[4] == "-";
                            temp_edges.push(CompactEdge::new(from_id, from_rev, to_id, to_rev));
                            edge_count += 1;
                        }
                        "P" => {
                            if fields.len() < 3 {
                                error!("Invalid P line: {line}");
                                std::process::exit(1);
                            }
                            let path_name = fields[1];
                            let nodes_str = fields[2];

                            if let Some((sample_hap_name, start, end)) = split_path_name(path_name)
                            {
                                let mut translated_steps = Vec::new();
                                for step_str in nodes_str.split(',') {
                                    if step_str.is_empty() {
                                        error!("Empty step in path {path_name} in line {line}");
                                        std::process::exit(1);
                                    }
                                    let (node_str, orient) = if let Some(stripped) =
                                        step_str.strip_suffix('+')
                                    {
                                        (stripped, false)
                                    } else if let Some(stripped) = step_str.strip_suffix('-') {
                                        (stripped, true)
                                    } else {
                                        error!("Invalid step format {step_str} in line {line}");
                                        std::process::exit(1);
                                    };

                                    let node_id: u64 = match node_str.parse() {
                                        Ok(id) => id,
                                        Err(_) => {
                                            error!(
                                                "Invalid node ID in path {path_name} in line {line}"
                                            );
                                            std::process::exit(1);
                                        }
                                    };

                                    if let Some(&translated_id) =
                                        id_translation.get(&NodeId::from(node_id))
                                    {
                                        translated_steps.push(Handle::pack(translated_id, orient));
                                    } else {
                                        error!(
                                            "Node {node_id} in path {path_name} not found in translation map"
                                        );
                                        std::process::exit(1);
                                    }
                                }
                                if !translated_steps.is_empty() {
                                    num_path_ranges.fetch_add(1, Ordering::Relaxed);
                                    num_path_range_steps
                                        .fetch_add(translated_steps.len(), Ordering::Relaxed);

                                    let mut map = path_key_ranges.lock().unwrap();
                                    map.entry(sample_hap_name.to_string()).or_default().push(
                                        RangeInfo {
                                            start,
                                            end,
                                            steps: translated_steps,
                                        },
                                    );
                                }
                            }
                        }
                        _ => {}
                    }
                }

                for edge in temp_edges {
                    if let (Some(&from_id), Some(&to_id)) = (
                        id_translation.get(&NodeId::from(edge.from_id())),
                        id_translation.get(&NodeId::from(edge.to_id())),
                    ) {
                        let mut graph = combined_graph.lock().unwrap();
                        graph.add_edge(
                            from_id.into(),
                            edge.from_rev(),
                            to_id.into(),
                            edge.to_rev(),
                        );
                    }
                }

                debug!(
                    "GFA file {gfa_id} ({gfa_path}) processed: {node_count} nodes, {edge_count} edges"
                );
            } else {
                error!("Failed to open GFA file '{gfa_path}'");
                std::process::exit(1);
            }
        });

    let path_map = Arc::try_unwrap(path_key_ranges)
        .expect("More than one Arc pointer to path map")
        .into_inner()
        .unwrap();

    // Unwrap the graph for lock-free validation
    let graph = Arc::try_unwrap(combined_graph)
        .ok()
        .expect("More than one Arc pointer to graph")
        .into_inner()
        .unwrap();

    if !skip_validation {
        // Pre-compute all unique node sequence lengths in parallel
        info!("Caching sequence lengths for validation");
        let unique_node_ids: FxHashSet<u64> = path_map
            .par_iter()
            .flat_map(|(_, ranges)| ranges.par_iter())
            .flat_map(|range| range.steps.par_iter())
            .map(|handle| u64::from(handle.id()))
            .collect();
        
        let node_length_cache: FxHashMap<u64, usize> = unique_node_ids
            .par_iter()
            .map(|&node_id| {
                let index = (node_id - 1) as usize; // Convert to 0-based index
                match graph.sequence_store.get_sequence_length(index) {
                    Ok(length) => (node_id, length),
                    Err(_) => (node_id, 0), // Handle error case
                }
            })
            .collect();
        
        // Validate all path ranges in parallel using cached lengths
        info!("Validating path range lengths");
        let validation_errors: Vec<String> = path_map
            .par_iter()
            .flat_map(|(path_key, ranges)| {
                ranges
                    .par_iter()
                    .filter_map(|range| {
                        let expected_length = range.end - range.start;
                        
                        // Calculate actual length using cached sequence lengths
                        let actual_length: usize = range.steps.iter()
                            .map(|handle| node_length_cache.get(&u64::from(handle.id())).copied().unwrap_or(0))
                            .sum();
                        
                        if expected_length != actual_length {
                            Some(format!(
                                "Path range length mismatch for '{}:{}-{}': \
                                expected length {} but sum of step lengths is {}",
                                path_key, range.start, range.end, expected_length, actual_length
                            ))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<String>>()
            })
            .collect();

        // Check if there were any validation errors
        if !validation_errors.is_empty() {
            for error in validation_errors {
                error!("{error}");
            }
            std::process::exit(1);
        }
    } else {
        info!("Skipping path range length validation (--skip-validation flag set)");
    }

    info!(
        "Collected {} nodes, {} edges, {} path keys, {} path ranges and {} path steps",
        graph.node_count,
        graph.edges.len(),
        path_map.len(),
        num_path_ranges.load(Ordering::Relaxed),
        num_path_range_steps.load(Ordering::Relaxed)
    );

    Ok((graph, path_map))
}

fn get_gfa_reader(gfa_path: &str) -> io::Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(gfa_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to open file '{gfa_path}': {e}"),
        )
    })?;

    let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
        io::Error::other(format!("Failed to open reader for '{gfa_path}': {e}"))
    })?;

    Ok(Box::new(BufReader::new(reader)))
}

fn split_path_name(path_name: &str) -> Option<(String, usize, usize)> {
    // Find the last ':' to split the range from the key
    if let Some(last_colon) = path_name.rfind(':') {
        let (key, range_str) = path_name.split_at(last_colon);
        // Skip the ':' character
        let range_str = &range_str[1..];

        // Find the '-' in the range portion
        if let Some((start_str, end_str)) = range_str.split_once('-') {
            if let (Ok(start), Ok(end)) = (start_str.parse(), end_str.parse()) {
                return Some((key.to_string(), start, end));
            }
        }
    }
    None
}

fn sort_and_filter_ranges(ranges: &mut Vec<RangeInfo>) {
    // Sort ranges by start position
    ranges.sort_by_key(|r| (r.start, r.end));

    debug!("  Removing redundant ranges");

    // Remove ranges that are contained within other ranges
    let mut write_idx = 0;
    for read_idx in 1..ranges.len() {
        let (prev_start, prev_end) = (ranges[write_idx].start, ranges[write_idx].end);
        let (curr_start, curr_end) = (ranges[read_idx].start, ranges[read_idx].end);

        if curr_start == prev_start && curr_end == prev_end {
            // Skip duplicate range
            // Current range is a duplicate of the previous range, skip it
            debug!(
                "    Duplicate range detected: Range [start={curr_start}, end={curr_end}] is identical to previous range and will be removed."
            );

            continue;
        } else if curr_start >= prev_start && curr_end <= prev_end {
            // Skip range that is fully contained within previous range
            debug!(
                "    Contained range detected: Range [start={curr_start}, end={curr_end}] is fully contained within previous range [start={prev_start}, end={prev_end}] and will be removed."
            );

            continue;
        } else if prev_start >= curr_start && prev_end <= curr_end {
            // Previous range is fully contained within current range
            debug!(
                "    Containing range detected: Previous range [start={prev_start}, end={prev_end}] is fully contained within current range [start={curr_start}, end={curr_end}] and will be removed."
            );

            ranges.swap(write_idx, read_idx);
        } else if curr_start < prev_end {
            // Handle overlapping ranges - check both previous and next ranges
            let mut should_skip = false;

            if read_idx < ranges.len() - 1 {
                let next_start = ranges[read_idx + 1].start;

                // Check if current range is significantly overlapped by both neighbors
                if curr_start > prev_start && next_start < curr_end {
                    let overlap_with_prev = prev_end - curr_start;
                    let overlap_with_next = curr_end - next_start;
                    let range_length = curr_end - curr_start;

                    // Skip if the range is mostly covered by its neighbors
                    if overlap_with_prev + overlap_with_next > range_length {
                        should_skip = true;
                    }
                }
            }

            if !should_skip {
                debug!(
                    "    Overlapping range detected: Range [start={curr_start}, end={curr_end}] overlaps with previous range [start={prev_start}, end={prev_end}] and will be kept."
                );
                write_idx += 1;
                if write_idx != read_idx {
                    ranges.swap(write_idx, read_idx);
                }
            }
        } else {
            // No overlap - keep both ranges
            write_idx += 1;
            if write_idx != read_idx {
                ranges.swap(write_idx, read_idx);
            }
        }
    }
    ranges.truncate(write_idx + 1);

    // if debug {
    //     debug!("  Path key '{}' without redundancy", path_key);
    //     for range in ranges.iter() {
    //         //debug!("    Range: start={}, end={}, num.steps={}, gfa_id={}", range.start, range.end, range.steps.len(), range.gfa_id);
    //         debug!("    Range: start={}, end={}, num.steps={}", range.start, range.end, range.steps.len());
    //     }
    // }
}

fn trim_range_overlaps(ranges: &mut [RangeInfo], graph_mutex: &Arc<Mutex<CompactGraph>>) {
    debug!("  Trimming overlapping ranges");

    for i in 1..ranges.len() {
        let (left, right) = ranges.split_at_mut(i);
        let r1 = &mut left[left.len() - 1];
        let r2 = &mut right[0];

        if r1.overlaps_with(r2) {
            // Calculate the overlap region - use max/min to get precise overlap bounds
            let overlap_start = std::cmp::max(r1.start, r2.start);
            let overlap_end = std::cmp::min(r1.end, r2.end);

            debug!(
                "    Overlap detected: Range1 [start={}, end={}], Range2 [start={}, end={}], Overlap [start={}, end={}, size={}]",
                r1.start, r1.end, r2.start, r2.end, overlap_start, overlap_end, overlap_end - overlap_start
            );

            // Adjust r2 to remove the overlap
            let mut steps_to_remove = Vec::new();
            let mut step_to_split: Option<usize> = None;
            let mut cumulative_pos = r2.start;

            // First pass: identify steps to remove/split (read-only, brief lock)
            {
                let graph = graph_mutex.lock().unwrap();
                for (idx, &step_handle) in r2.steps.iter().enumerate() {
                    let step_start = cumulative_pos;
                    let node_seq = graph.get_sequence(step_handle).unwrap();
                    let node_length = node_seq.len();
                    cumulative_pos += node_length;
                    let step_end = cumulative_pos;

                    if step_end <= overlap_start {
                        continue;
                    } else if step_start >= overlap_end {
                        break;
                    } else if step_start >= overlap_start && step_end <= overlap_end {
                        steps_to_remove.push(idx);
                    } else {
                        if step_to_split.is_some() {
                            error!("Error: More than one step is partially overlapping, which is not allowed.");
                            std::process::exit(1);
                        }
                        step_to_split = Some(idx);
                    }
                }
            } // Lock released here

            // Initialize new vectors to store updated steps
            let mut new_steps = Vec::with_capacity(r2.steps.len() / 2);
            let mut range_new_start = None;
            let mut current_pos = None;

            // Reset cumulative position for second pass
            cumulative_pos = r2.start;

            // Second pass: Iterate over the original steps using incrementally computed positions
            for (idx, &step_handle) in r2.steps.iter().enumerate() {
                let step_start = cumulative_pos;

                // Get node sequence with lock
                let node_seq = {
                    let graph = graph_mutex.lock().unwrap();
                    graph.get_sequence(step_handle).unwrap()
                }; // Lock released here

                let node_length = node_seq.len();
                cumulative_pos += node_length;
                let step_end = cumulative_pos;

                if steps_to_remove.contains(&idx) {
                    // Skip steps to remove
                    continue;
                } else if step_to_split == Some(idx) {
                    // Split node for the single partially overlapping step
                    let overlap_within_step_start = std::cmp::max(step_start, overlap_start);
                    let overlap_within_step_end = std::cmp::min(step_end, overlap_end);

                    // Calculate offsets relative to the node sequence
                    let node_len = node_seq.len();
                    let overlap_start_offset =
                        (overlap_within_step_start - step_start).min(node_len);
                    let overlap_end_offset = (overlap_within_step_end - step_start).min(node_len);

                    debug!("      Splitting step {} [start={}, end={}, len={}] to remove overlap at [start={}, end={}], Overlap offsets: start={}, end={}",
                        idx, step_start, step_end, step_end - step_start, overlap_within_step_start, overlap_within_step_end, overlap_start_offset, overlap_end_offset);

                    if step_start < overlap_start {
                        debug!(
                            "      Adding left part of step [start={step_start}, end={overlap_within_step_start}]"
                        );
                        assert!(overlap_start_offset > 0);

                        // Keep left part
                        let new_seq = node_seq[0..overlap_start_offset].to_vec();

                        // Add node with lock
                        let node_id = {
                            let mut graph = graph_mutex.lock().unwrap();
                            graph.add_node(&new_seq).unwrap()
                        }; // Lock released here

                        let new_node = Handle::pack(NodeId::from(node_id), false);

                        new_steps.push(new_node);
                        if range_new_start.is_none() {
                            range_new_start = Some(step_start);
                            current_pos = Some(step_start);
                        }
                    } else if step_end > overlap_end {
                        debug!(
                            "      Adding right part of step [start={overlap_within_step_end}, end={step_end}]"
                        );
                        assert!(overlap_end_offset < node_len);

                        // Keep right part
                        let new_seq = node_seq[overlap_end_offset..].to_vec();

                        // Add node with lock
                        let node_id = {
                            let mut graph = graph_mutex.lock().unwrap();
                            graph.add_node(&new_seq).unwrap()
                        }; // Lock released here

                        let new_node = Handle::pack(NodeId::from(node_id), false);

                        new_steps.push(new_node);
                        if range_new_start.is_none() {
                            range_new_start = Some(overlap_end);
                            current_pos = Some(overlap_end);
                        }
                        current_pos = Some(current_pos.unwrap() + new_seq.len());
                    }
                } else {
                    // Keep steps that are not to be removed or split
                    new_steps.push(step_handle);
                    if range_new_start.is_none() {
                        range_new_start = Some(step_start);
                        current_pos = Some(step_start);
                    }
                    current_pos = Some(current_pos.unwrap() + node_length);
                }
            }

            // Update r2 with the new steps
            r2.steps = new_steps;

            // Update edges for the modified steps
            for idx in 0..r2.steps.len() {
                if idx > 0 {
                    let prev_step = r2.steps[idx - 1];
                    let curr_step = r2.steps[idx];

                    // Check and add edge with lock
                    let mut graph = graph_mutex.lock().unwrap();
                    if !graph.has_edge(
                        prev_step.id().into(),
                        prev_step.is_reverse(),
                        curr_step.id().into(),
                        curr_step.is_reverse(),
                    ) {
                        debug!(
                            "      Creating edge between steps: {} -> {}",
                            prev_step.id(),
                            curr_step.id()
                        );
                        graph.add_edge(
                            prev_step.id().into(),
                            prev_step.is_reverse(),
                            curr_step.id().into(),
                            curr_step.is_reverse(),
                        );
                    }
                    // Lock released at end of scope
                }
            }

            // Update r2.start and r2.end based on the new step positions
            if !r2.steps.is_empty() {
                r2.start = range_new_start.unwrap();
                r2.end = current_pos.unwrap();
            } else {
                // If no steps remain, set start and end to overlap_end to effectively remove this range
                r2.start = overlap_end;
                r2.end = overlap_end;
            }

            debug!(
                "      Updated overlaps: Range2 [start={}, end={}]",
                r2.start, r2.end
            );
        }
    }
}

fn link_contiguous_ranges(ranges: &[RangeInfo], graph_mutex: &Arc<Mutex<CompactGraph>>) {
    debug!("  Linking contiguous ranges");

    for i in 1..ranges.len() {
        let r1 = &ranges[i - 1];
        let r2 = &ranges[i];

        // Check if ranges are contiguous
        if r1.is_contiguous_with(r2) {
            // Get last handle from previous range and first handle from current range
            if let (Some(&last_handle), Some(&first_handle)) = (r1.steps.last(), r2.steps.first()) {
                // Lock only for checking and adding edge
                let mut graph = graph_mutex.lock().unwrap();
                if !graph.has_edge(
                    last_handle.id().into(),
                    last_handle.is_reverse(),
                    first_handle.id().into(),
                    first_handle.is_reverse(),
                ) {
                    debug!(
                        "    Creating edge between contiguous ranges at position {}: {} -> {}",
                        r1.end,
                        last_handle.id(),
                        first_handle.id()
                    );
                    graph.add_edge(
                        last_handle.id().into(),
                        last_handle.is_reverse(),
                        first_handle.id().into(),
                        first_handle.is_reverse(),
                    );
                }
                // Lock released at end of scope
            }
        }
    }
}

fn mark_nodes_for_removal(
    node_count: u64,
    path_key_ranges: &FxHashMap<String, Vec<RangeInfo>>,
) -> BitVec {
    // Create a bitvector with all nodes initially marked for removal
    let mut nodes_to_remove = bitvec![1; node_count as usize + 1]; // +1 to account for 0-indexing

    // Mark nodes used in path ranges as not to be removed (set bit to 0)
    for ranges in path_key_ranges.values() {
        for range in ranges {
            for handle in &range.steps {
                nodes_to_remove.set(u64::from(handle.id()) as usize, false);
            }
        }
    }

    nodes_to_remove
}

fn write_graph_to_gfa(
    combined_graph: &mut CompactGraph,
    path_key_ranges: &FxHashMap<String, Vec<RangeInfo>>,
    output_path: &str,
    compression_format: Format,
    fill_gaps: u8,
    sequence_index: Option<&UnifiedSequenceIndex>,
    debug: bool,
) -> std::io::Result<()> {
    info!("Marking unused nodes");
    let nodes_to_remove = mark_nodes_for_removal(combined_graph.node_count, path_key_ranges);
    debug!("Marked {} unused nodes", nodes_to_remove.count_ones() - 1);

    // Create the output file
    let output_file = File::create(output_path)?;

    // Create writer based on compression format
    let writer: Box<dyn Write> = match compression_format {
        Format::Gzip => {
            // Use parallel gzip compression
            let parz: ParCompress<Gzip> = ParCompressBuilder::new()
                .num_threads(rayon::current_num_threads())
                .map_err(|e| std::io::Error::other(format!("Failed to set threads: {e:?}")))?
                .compression_level(Compression::new(6))
                .from_writer(output_file);
            Box::new(parz)
        }
        Format::Bzip => {
            // Use parallel BGZF compression
            let parz: ParCompress<Bgzf> = ParCompressBuilder::new()
                .num_threads(rayon::current_num_threads())
                .map_err(|e| std::io::Error::other(format!("Failed to set threads: {e:?}")))?
                .compression_level(Compression::new(6))
                .from_writer(output_file);
            Box::new(parz)
        }
        Format::Zstd => {
            // Use multi-threaded zstd compression
            let mut encoder = ZstdEncoder::new(output_file, 6)?;
            encoder.multithread(rayon::current_num_threads() as u32)?;
            Box::new(encoder)
        }
        Format::No => {
            // No compression
            Box::new(output_file)
        }
        _ => {
            // Fallback to niffler for other formats
            niffler::get_writer(
                Box::new(output_file),
                compression_format,
                niffler::compression::Level::Six,
            )
            .map_err(std::io::Error::other)?
        }
    };

    let mut file = BufWriter::new(writer);

    // Write GFA version
    writeln!(file, "H\tVN:Z:1.0")?;

    // Write nodes by excluding marked ones and create the id_mapping
    info!("Writing used nodes by compacting their IDs");
    let max_id = combined_graph.node_count as usize;
    let mut id_mapping = vec![0; max_id + 1];
    let mut new_id = 1; // Start from 1

    for node_id in 1..=combined_graph.node_count {
        if !nodes_to_remove[node_id as usize] {
            id_mapping[node_id as usize] = new_id;

            let sequence =
                combined_graph.get_sequence(Handle::pack(NodeId::from(node_id), false))?;
            let sequence_str =
                String::from_utf8(sequence).expect("Node sequence contains invalid UTF-8");
            writeln!(file, "S\t{new_id}\t{sequence_str}")?;

            new_id += 1;
        }
    }

    info!("Writing edges connecting used nodes");
    for edge in &combined_graph.edges {
        let from_id = edge.from_id() as usize;
        let to_id = edge.to_id() as usize;

        if !nodes_to_remove[from_id] && !nodes_to_remove[to_id] {
            let from_mapped = id_mapping[from_id];
            let to_mapped = id_mapping[to_id];
            let from_orient = if edge.from_rev() { "-" } else { "+" };
            let to_orient = if edge.to_rev() { "-" } else { "+" };
            writeln!(
                file,
                "L\t{from_mapped}\t{from_orient}\t{to_mapped}\t{to_orient}\t0M"
            )?;
        }
    }

    // Write paths by processing ranges directly
    info!("Writing paths by merging contiguous path ranges");
    let mut path_key_vec: Vec<_> = path_key_ranges.keys().collect();
    path_key_vec.par_sort_unstable(); // Sort path keys for consistent output (for path keys, the order of equal elements doesn't matter since they're unique)

    let mut start_gaps = 0;
    let mut middle_gaps = 0;
    let mut end_gaps = 0;

    // Check if a valid sequence index is provided for end gap filling
    if fill_gaps == 2 && sequence_index.is_none() {
        warn!("Cannot fill end gaps without sequence files; trailing gaps will be skipped");
    }

    for path_key in path_key_vec {
        let ranges = &path_key_ranges[path_key];
        //if ranges.is_empty() { continue }

        if debug {
            debug!("Processing Path key '{path_key}'");

            let mut current_start = ranges[0].start;
            let mut current_end = ranges[0].end;

            for i in 1..ranges.len() {
                if ranges[i - 1].is_contiguous_with(&ranges[i]) {
                    // Extend current merged range
                    current_end = ranges[i].end;
                } else {
                    // Print current merged range
                    debug!(
                        "  Merged range: start={current_start}, end={current_end}"
                    );

                    if !ranges[i - 1].overlaps_with(&ranges[i]) {
                        // Calculate and print gap
                        let gap = ranges[i].start - current_end;
                        debug!("    Gap to next range: {gap} positions");
                    } else {
                        // Calculate and print overlap (IT SHOULD NOT HAPPEN)
                        let overlap = current_end - ranges[i].start;
                        debug!("    Overlap with next range: {overlap} positions");
                    }

                    // Start new merged range
                    current_start = ranges[i].start;
                    current_end = ranges[i].end;
                }
            }

            // Print final merged range
            debug!(
                "  Final merged range: start={current_start}, end={current_end}"
            );
        }

        let mut path_elements: Vec<String> = Vec::new();
        let mut path_start = ranges[0].start;
        let mut path_end = ranges[0].start;

        // Handle initial gap if it exists and gap filling is enabled
        if fill_gaps == 2 && ranges[0].start > 0 {
            start_gaps += 1;
            path_elements.push(create_gap_node(
                &mut file,
                (0, ranges[0].start),
                path_key,
                sequence_index,
                None, // No previous element for initial gap
                ranges[0].steps.first(),
                &id_mapping,
                &mut new_id,
            )?);
            path_start = 0;
        }

        // Process subsequent contiguous ranges or add gap nodes
        let mut i = 0;
        while i < ranges.len() {
            // Merge a block of contiguous ranges
            add_range_steps_to_path(&ranges[i], &id_mapping, &mut path_elements);
            let mut last_range_end = ranges[i].end;
            i += 1;

            while i < ranges.len() && ranges[i - 1].is_contiguous_with(&ranges[i]) {
                add_range_steps_to_path(&ranges[i], &id_mapping, &mut path_elements);
                last_range_end = ranges[i].end;
                i += 1;
            }
            path_end = last_range_end;

            // We're now at the end of a contiguous block
            if i < ranges.len() {
                // there is a gap before the next block
                let next_start = ranges[i].start;

                if fill_gaps > 0 {
                    // Bridge the gap
                    middle_gaps += 1;

                    path_elements.push(create_gap_node(
                        &mut file,
                        (last_range_end, next_start),
                        path_key,
                        sequence_index,
                        path_elements.last(),
                        ranges[i].steps.first(),
                        &id_mapping,
                        &mut new_id,
                    )?);
                } else {
                    // Finish current partial path and start a new one
                    if !path_elements.is_empty() {
                        let path_name = format!("{path_key}:{path_start}-{path_end}");
                        writeln!(file, "P\t{}\t{}\t*", path_name, path_elements.join(","))?;
                        path_elements.clear();
                    }
                    path_start = next_start;
                }
            }
        }

        // Handle final gap if it exists and gap filling is enabled
        if fill_gaps == 2 && sequence_index.is_some() {
            // Try to get the sequence length for this path_key
            match sequence_index.unwrap().get_sequence_length(path_key) {
                Ok(total_len) => {
                    if path_end < total_len {
                        end_gaps += 1;

                        path_elements.push(create_gap_node(
                            &mut file,
                            (path_end, total_len),
                            path_key,
                            sequence_index,
                            path_elements.last(),
                            None, // No next node for final gap
                            &id_mapping,
                            &mut new_id,
                        )?);
                        path_end = total_len;
                    } else if path_end > total_len {
                        warn!(
                            "Path '{path_key}' extends beyond sequence length ({path_end} > {total_len})"
                        );
                    }
                }
                Err(e) => {
                    warn!(
                        "Failed to get sequence length for path '{path_key}': {e}"
                    );
                }
            }
        }

        // Write the remaining path
        if !path_elements.is_empty() {
            let is_full_path = path_start == 0
                && sequence_index
                    .map(|idx| {
                        idx.get_sequence_length(path_key)
                            .map(|len| len == path_end)
                            .unwrap_or(false)
                    })
                    .unwrap_or(true); // Assume path_end matches full length if no sequence index is available

            let path_name = if is_full_path {
                path_key.to_string()
            } else {
                format!("{path_key}:{path_start}-{path_end}")
            };
            writeln!(file, "P\t{}\t{}\t*", path_name, path_elements.join(","))?;
        }
    }

    if fill_gaps == 2 {
        info!(
            "Filled {} gaps: {} start gaps, {} middle gaps, {} end gaps",
            start_gaps + middle_gaps + end_gaps,
            start_gaps,
            middle_gaps,
            end_gaps
        );
    } else if fill_gaps == 1 {
        info!("Filled {middle_gaps} middle gaps");
    }

    file.flush()?;
    Ok(())
}

fn create_gap_node<W: Write>(
    file: &mut BufWriter<W>,
    gap_range: (usize, usize),
    path_key: &str,
    sequence_index: Option<&UnifiedSequenceIndex>,
    last_element: Option<&String>,
    next_handle: Option<&Handle>,
    id_mapping: &[usize],
    new_id: &mut usize,
) -> io::Result<String> {
    let (gap_start, gap_end) = gap_range;
    let gap_size = gap_end - gap_start;

    // Get gap sequence either from sequence index or create string of N's
    let gap_sequence = if let Some(index) = sequence_index {
        match index.fetch_sequence(path_key, gap_start as i32, gap_end as i32) {
            Ok(seq) => String::from_utf8(seq).unwrap_or_else(|e| {
                error!("Failed to convert sequence to UTF-8: {e}");
                "N".repeat(gap_size)
            }),
            Err(e) => {
                error!("Failed to fetch sequence for '{path_key}': {e}");
                "N".repeat(gap_size)
            }
        }
    } else {
        "N".repeat(gap_size)
    };

    // Write gap node
    writeln!(file, "S\t{new_id}\t{gap_sequence}")?;

    // Add edge from previous node if it exists
    if let Some(last_element) = last_element {
        let last_id = last_element[..last_element.len() - 1]
            .parse::<usize>()
            .unwrap();
        let last_orient = &last_element[last_element.len() - 1..];
        writeln!(file, "L\t{last_id}\t{last_orient}\t{new_id}\t+\t0M")?;
    }

    // Add edge to next node if it exists
    if let Some(handle) = next_handle {
        let next_id = id_mapping[u64::from(handle.id()) as usize];
        let next_orient = if handle.is_reverse() { "-" } else { "+" };
        writeln!(file, "L\t{}\t+\t{}\t{}\t0M", *new_id, next_id, next_orient)?;
    }

    let path_element = format!("{new_id}+");
    *new_id += 1;

    Ok(path_element)
}

fn add_range_steps_to_path(
    range: &RangeInfo,
    id_mapping: &[usize],
    path_elements: &mut Vec<String>,
) {
    for handle in &range.steps {
        let node_id = id_mapping[u64::from(handle.id()) as usize];
        let orient = if handle.is_reverse() { "-" } else { "+" };
        path_elements.push(format!("{node_id}{orient}"));
    }
}


/// Parse VCF CHROM field to extract base contig and range
fn parse_vcf_chrom(chrom: &str) -> Option<(String, u64, u64)> {
    if let Some(colon_pos) = chrom.rfind(':') {
        let (base_contig, range_str) = chrom.split_at(colon_pos);
        let range_str = &range_str[1..]; // Skip the ':'
        
        if let Some(dash_pos) = range_str.find('-') {
            let (start_str, end_str) = range_str.split_at(dash_pos);
            let end_str = &end_str[1..]; // Skip the '-'
            
            if let (Ok(start), Ok(end)) = (start_str.parse::<u64>(), end_str.parse::<u64>()) {
                return Some((base_contig.to_string(), start, end));
            }
        }
    }
    None
}

/// Get chromosome sort key for human-friendly ordering
fn chr_sort_key(base_contig: &str) -> (u8, u64, String) {
    // Extract the last segment after '#', e.g., "chr19" from "CHM13#0#chr19"
    let chr_label = base_contig.split('#').last().unwrap_or(base_contig);
    
    if chr_label.starts_with("chr") {
        let suffix = &chr_label[3..];
        
        // Numbered chromosome
        if let Ok(num) = suffix.parse::<u64>() {
            if (1..=22).contains(&num) {
                return (0, num, String::new());
            }
        }
        
        // X chromosome
        if suffix == "X" {
            return (0, 23, String::new());
        }
        
        // Y chromosome
        if suffix == "Y" {
            return (0, 24, String::new());
        }
        
        // Mitochondrial: accept "M" or "MT"
        if suffix == "M" || suffix == "MT" {
            return (0, 25, String::new());
        }
    }
    
    // Fallback: anything not matching above goes after, sorted alphabetically
    (1, 0, chr_label.to_string())
}

/// Lace VCF files together
pub fn run_vcf_lace(
    vcf_files: Option<Vec<String>>,
    vcf_list: Option<String>,
    output: &str,
    compress: &str,
    threads: usize,
    verbose: u8,
    reference_index: Option<&UnifiedSequenceIndex>,
) -> io::Result<()> {
    // Determine compression format
    let compression_format = get_compression_format(compress, output);
    
    // Resolve VCF files
    let vcf_files = resolve_vcf_files(vcf_files, vcf_list)?;
    if vcf_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No VCF files specified",
        ));
    }
    
    info!("[1/2] Found {} VCF files to merge", vcf_files.len());
    
    // Process files in parallel and collect results as they complete
    let mut all_samples = FxHashSet::default();
    let mut contig_max_end = FxHashMap::default();
    let mut file_order_info = Vec::new();
    
    info!("[1/2] Scanning files with {} threads...", threads);
    
    // Use channels to process results as they become available
    use std::sync::mpsc;
    let (tx, rx) = mpsc::channel();
    
    // Process files in parallel using rayon::scope
    let result = rayon::scope(|s| -> io::Result<()> {
        s.spawn(|_| {
            vcf_files
                .par_iter()
                .enumerate()
                .for_each_with(tx, |tx, (_idx, path)| {
                    let result = process_vcf_file(path, _idx);
                    let _ = tx.send(result);
                });
        });
        
        // Process results as they become available
        let mut processed_count = 0;
        for result in rx {
            match result {
                Ok((path, local_samples, local_contigs, local_order)) => {
                    processed_count += 1;
                    
                    info!("[1/2] Processed {}/{}: {}", 
                          processed_count, vcf_files.len(), 
                          Path::new(&path).file_name().unwrap().to_string_lossy());
                    
                    // Update global sample set
                    all_samples.extend(local_samples);
                    
                    // Update global contig extents
                    for (base_contig, end_pos) in local_contigs {
                        let current_max = contig_max_end.entry(base_contig).or_insert(0);
                        if end_pos > *current_max {
                            *current_max = end_pos;
                        }
                    }
                    
                    // Record ordering info
                    file_order_info.push((path, local_order));
                    
                    // Break when all files are processed
                    if processed_count == vcf_files.len() {
                        break;
                    }
                }
                Err(e) => {
                    return Err(e);
                }
            }
        }
        Ok(())
    });
    
    // Handle any errors from the scope
    result?;
    
    let merged_samples: Vec<String> = {
        let mut samples = all_samples.into_iter().collect::<Vec<_>>();
        samples.sort();
        samples
    };
    
    info!("[1/2] Collected {} unique samples", merged_samples.len());
    info!("[1/2] Identified {} contigs", contig_max_end.len());
    
    // Sort files by ordering key
    file_order_info.sort_by(|a, b| a.1.cmp(&b.1));
    let sorted_paths: Vec<String> = file_order_info.into_iter().map(|(path, _)| path).collect();
    
    // Second pass: merge and write VCF
    write_merged_vcf(
        &sorted_paths,
        &merged_samples,
        &contig_max_end,
        output,
        compression_format,
        verbose,
        reference_index,
    )?;
    
    info!("[2/2] Done: merged {} files to {}", vcf_files.len(), output);
    Ok(())
}

/// Resolve VCF files from either --vcf-files or --vcf-list
fn resolve_vcf_files(
    vcf_files: Option<Vec<String>>,
    vcf_list: Option<String>,
) -> io::Result<Vec<String>> {
    match (vcf_files, vcf_list) {
        (Some(files), None) => {
            // Validate all files exist
            for file in &files {
                if !Path::new(file).exists() {
                    return Err(io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("VCF file '{}' not found", file),
                    ));
                }
            }
            Ok(files)
        }
        (None, Some(list_file)) => {
            let file = File::open(&list_file)?;
            let reader = BufReader::new(file);
            let mut files = Vec::new();
            
            for line in reader.lines() {
                let line = line?;
                let trimmed = line.trim();
                if !trimmed.is_empty() && !trimmed.starts_with('#') {
                    if !Path::new(trimmed).exists() {
                        return Err(io::Error::new(
                            io::ErrorKind::NotFound,
                            format!("VCF file '{}' not found", trimmed),
                        ));
                    }
                    files.push(trimmed.to_string());
                }
            }
            
            if files.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("No valid VCF files found in list file: {}", list_file),
                ));
            }
            
            Ok(files)
        }
        (None, None) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Either --vcf-files or --vcf-list must be provided",
        )),
        (Some(_), Some(_)) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Cannot specify both --vcf-files and --vcf-list",
        )),
    }
}

/// Process a single VCF file to extract metadata
fn process_vcf_file(
    path: &str,
    _file_idx: usize,
) -> io::Result<(String, Vec<String>, FxHashMap<String, u64>, (u8, u64, String, u64))> {
    let mut local_samples = Vec::new();
    let mut local_contigs = FxHashMap::default();
    let mut local_order: Option<(u8, u64, String, u64)> = None;
    
    let reader = get_vcf_reader(path)?;
    
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        
        if line.is_empty() {
            continue;
        }
        
        if line.starts_with("##") {
            continue;
        } else if line.starts_with("#CHROM") {
            // Parse header to get sample names
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() > 9 {
                local_samples = parts[9..].iter().map(|s| s.to_string()).collect();
            }
        } else {
            // Data line
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 2 {
                continue;
            }
            
            let chrom = parts[0];
            let pos_str = parts[1];
            
            if let Some((base_contig, start, end)) = parse_vcf_chrom(chrom) {
                // Parse position to ensure it's valid, but we don't use the value
                if let Ok(_pos) = pos_str.parse::<u64>() {
                    // Update contig max position using the end value from CHROM field
                    let current_max = local_contigs.entry(base_contig.clone()).or_insert(0);
                    if end > *current_max {
                        *current_max = end;
                    }
                    
                    // Update file ordering
                    let sort_key = chr_sort_key(&base_contig);
                    let file_key = (sort_key.0, sort_key.1, sort_key.2, start);
                    
                    if local_order.is_none() || file_key < *local_order.as_ref().unwrap() {
                        local_order = Some(file_key);
                    }
                }
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Unexpected CHROM format in {}: {}", path, chrom),
                ));
            }
        }
    }
    
    let order = local_order.unwrap_or((2, 0, String::new(), 0));
    Ok((path.to_string(), local_samples, local_contigs, order))
}

/// Get VCF reader that handles compression transparently
fn get_vcf_reader(path: &str) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    
    // Use niffler to automatically detect and handle compression
    let (reader, _format) = niffler::get_reader(Box::new(file))
        .map_err(|e| io::Error::other(format!("Failed to open reader for '{}': {}", path, e)))?;
    
    Ok(Box::new(BufReader::new(reader)))
}

/// Write merged VCF file
fn write_merged_vcf(
    sorted_paths: &[String],
    merged_samples: &[String],
    contig_max_end: &FxHashMap<String, u64>,
    output_path: &str,
    compression_format: Format,
    verbose: u8,
    reference_index: Option<&UnifiedSequenceIndex>,
) -> io::Result<()> {
    // Create output writer with compression
    let output_file = File::create(output_path)?;
    let writer: Box<dyn Write> = match compression_format {
        Format::Gzip => {
            let parz: ParCompress<Gzip> = ParCompressBuilder::new()
                .num_threads(rayon::current_num_threads())
                .map_err(|e| std::io::Error::other(format!("Failed to set threads: {:?}", e)))?
                .compression_level(Compression::new(6))
                .from_writer(output_file);
            Box::new(parz)
        }
        Format::No => Box::new(output_file),
        _ => {
            // Use niffler for other formats
            niffler::get_writer(
                Box::new(output_file),
                compression_format,
                niffler::compression::Level::Six,
            )
            .map_err(std::io::Error::other)?
        }
    };
    
    let mut file = BufWriter::new(writer);
    
    // Write VCF header
    writeln!(file, "##fileformat=VCFv4.2")?;
    
    // Read and write meta-lines from the first file
    if let Some(first_path) = sorted_paths.first() {
        let reader = get_vcf_reader(first_path)?;
        for line in reader.lines() {
            let line = line?;
            if line.starts_with("##") {
                // Skip fileformat line (already written) and contig lines (will be regenerated)
                if line.starts_with("##fileformat") || line.starts_with("##contig") {
                    continue;
                }
                writeln!(file, "{}", line)?;
            } else {
                break; // Stop at header line
            }
        }
    }
    
    // Write new contig lines sorted by chromosome order with validation
    let mut sorted_contigs: Vec<_> = contig_max_end.iter().collect();
    sorted_contigs.sort_by(|a, b| chr_sort_key(a.0).cmp(&chr_sort_key(b.0)));
    
    for (base_contig, &estimated_length) in sorted_contigs {
        let final_length = if let Some(ref_index) = reference_index {
            match ref_index.get_sequence_length(base_contig) {
                Ok(ref_length) => {
                    if estimated_length > ref_length as u64 {
                        warn!(
                            "Contig '{}': Estimated length {} exceeds reference length {}, using reference length",
                            base_contig, estimated_length, ref_length
                        );
                        ref_length as u64
                    } else if estimated_length < ref_length as u64 {
                        warn!(
                            "Contig '{}': Estimated length {} is shorter than reference length {}, using reference length",
                            base_contig, estimated_length, ref_length
                        );
                        ref_length as u64
                    } else {
                        estimated_length
                    }
                }
                Err(_) => {
                    warn!(
                        "Contig '{}' not found in reference, using estimated length {}",
                        base_contig, estimated_length
                    );
                    estimated_length
                }
            }
        } else {
            estimated_length
        };
        
        writeln!(file, "##contig=<ID={},length={}>", base_contig, final_length)?;
    }
    
    // Write header line
    let header_cols = vec!["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        .into_iter()
        .chain(merged_samples.iter().map(|s| s.as_str()))
        .collect::<Vec<_>>();
    writeln!(file, "{}", header_cols.join("\t"))?;
    
    info!("[2/2] Writing merged VCF records");
    
    // Process each file in sorted order
    for (idx, path) in sorted_paths.iter().enumerate() {
        if verbose > 0 {
            info!("[2/2] Merging {}/{}: {}", 
                  idx + 1, sorted_paths.len(), 
                  Path::new(path).file_name().unwrap().to_string_lossy());
        }
        
        merge_vcf_file_records(&mut file, path, merged_samples)?;
    }
    
    file.flush()?;
    Ok(())
}

/// Merge records from a single VCF file
fn merge_vcf_file_records<W: Write>(
    output_file: &mut BufWriter<W>,
    path: &str,
    merged_samples: &[String],
) -> io::Result<()> {
    let reader = get_vcf_reader(path)?;
    let mut this_samples = Vec::new();
    let mut in_header = true;
    
    // Cache for missing genotype strings per FORMAT field
    let mut format_cache: FxHashMap<String, String> = FxHashMap::default();
    
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        
        if line.is_empty() {
            continue;
        }
        
        if in_header {
            if line.starts_with("#CHROM") {
                // Parse sample names from this file
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() > 9 {
                    this_samples = parts[9..].iter().map(|s| s.to_string()).collect();
                }
                in_header = false;
            }
            continue;
        }
        
        // Process data record
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 {
            error!("Malformed VCF record in {}: {}", path, line);
            std::process::exit(1);
        }
        
        let chrom = parts[0];
        let pos_str = parts[1];
        let id = parts[2];
        let ref_allele = parts[3];
        let alt_alleles = parts[4];
        let qual = parts[5];
        let filter = parts[6];
        let info = parts[7];
        let format = parts[8];
        let sample_cols = &parts[9..];
        
        // Parse and adjust coordinates
        if let Some((base_contig, start, _end)) = parse_vcf_chrom(chrom) {
            if let Ok(old_pos) = pos_str.parse::<u64>() {
                let new_pos = start + old_pos;
                
                // Get cached missing genotype string or compute it
                let missing_str = format_cache.entry(format.to_string()).or_insert_with(|| {
                    let format_keys: Vec<&str> = format.split(':').collect();
                    let missing_fields: Vec<&str> = format_keys
                        .iter()
                        .map(|&key| if key == "GT" { "./." } else { "." })
                        .collect();
                    missing_fields.join(":")
                });
                
                // Map this file's samples to their genotype strings
                let sample_to_gt: FxHashMap<String, String> = this_samples
                    .iter()
                    .zip(sample_cols.iter())
                    .map(|(sample, gt)| (sample.clone(), gt.to_string()))
                    .collect();
                
                // Build output genotypes in merged sample order
                let out_genotypes: Vec<String> = merged_samples
                    .iter()
                    .map(|sample| {
                        sample_to_gt.get(sample)
                            .cloned()
                            .unwrap_or_else(|| missing_str.clone())
                    })
                    .collect();
                
                // Write output record
                let out_fields = vec![
                    base_contig,
                    new_pos.to_string(),
                    id.to_string(),
                    ref_allele.to_string(),
                    alt_alleles.to_string(),
                    qual.to_string(),
                    filter.to_string(),
                    info.to_string(),
                    format.to_string(),
                ]
                .into_iter()
                .chain(out_genotypes.into_iter())
                .collect::<Vec<_>>();
                
                writeln!(output_file, "{}", out_fields.join("\t"))?;
            } else {
                error!("Cannot parse POS in {}: {}", path, pos_str);
                std::process::exit(1);
            }
        } else {
            error!("Unexpected CHROM format in {}: {}", path, chrom);
            std::process::exit(1);
        }
    }
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create a simple RangeInfo for testing
    fn create_range_info(start: usize, end: usize) -> RangeInfo {
        RangeInfo {
            start,
            end,
            //gfa_id,
            steps: vec![], // Empty steps for testing
        }
    }

    #[test]
    fn test_range_containment_removal() {
        // Test cases
        let test_cases = vec![
            // Test case 1: Basic containment
            (
                vec![(10, 50), (20, 30)], // Input ranges (start, end)
                vec![(10, 50)],           // Expected result
                "Basic containment",
            ),
            // Test case 2: No containment
            (
                vec![(10, 20), (30, 40)],
                vec![(10, 20), (30, 40)],
                "No containment",
            ),
            // Test case 3: Multiple contained ranges
            (
                vec![(10, 100), (20, 30), (40, 50), (60, 70)],
                vec![(10, 100)],
                "Multiple contained ranges",
            ),
            // Test case 4: Identical ranges
            (vec![(10, 20), (10, 20)], vec![(10, 20)], "Identical ranges"),
            // Test case 5: Nested containment
            (
                vec![(10, 100), (20, 80), (30, 40)],
                vec![(10, 100)],
                "Nested containment",
            ),
            // Test case 6: Partial overlap (should keep both)
            (
                vec![(10, 30), (20, 40)],
                vec![(10, 30), (20, 40)],
                "Partial overlap",
            ),
            // Test case 7: Edge cases - touching ranges
            (
                vec![(10, 20), (20, 30)],
                vec![(10, 20), (20, 30)],
                "Touching ranges",
            ),
            // Test case 8: Overlapping ranges from same GFA
            (
                vec![(0, 11742), (9714, 13000), (11000, 19000)],
                vec![(0, 11742), (11000, 19000)],
                "Overlapping ranges from same GFA",
            ),
            // Test case 9: Overlapping ranges with different GFA IDs
            (
                vec![(0, 11742), (9714, 13000), (11000, 19000)],
                vec![(0, 11742), (11000, 19000)],
                "Overlapping ranges",
            ),
            // Test case 10: Overlapping ranges with different GFA IDs 2
            (
                vec![(0, 10), (8, 20), (15, 30)],
                vec![(0, 10), (8, 20), (15, 30)],
                "Overlapping ranges",
            ),
            // Test case 11: Overlapping ranges with different GFA IDs 3
            (
                vec![(8000, 11000), (9694, 12313), (10908, 13908)],
                vec![(8000, 11000), (10908, 13908)],
                "Overlapping ranges",
            ),
        ];

        // Run each test case
        for (case_index, (input_ranges, expected_ranges, case_name)) in
            test_cases.iter().enumerate()
        {
            println!("Running test case {}: {}", case_index + 1, case_name);

            // Create input ranges
            let mut ranges: Vec<RangeInfo> = input_ranges
                .iter()
                .map(|(start, end)| create_range_info(*start, *end))
                .collect();

            // Sort ranges by start position
            ranges.sort_by_key(|r| (r.start, r.end));

            let mut write_idx = 0;
            for read_idx in 1..ranges.len() {
                let (prev_start, prev_end) = (ranges[write_idx].start, ranges[write_idx].end);
                let (curr_start, curr_end) = (ranges[read_idx].start, ranges[read_idx].end);

                if curr_start == prev_start && curr_end == prev_end {
                    // Skip duplicate range
                    continue;
                } else if curr_start >= prev_start && curr_end <= prev_end {
                    // Skip range that is fully contained within previous range
                    continue;
                } else if prev_start >= curr_start && prev_end <= curr_end {
                    // Previous range is fully contained within current range
                    ranges.swap(write_idx, read_idx);
                } else if curr_start < prev_end {
                    // Handle overlapping ranges - check both previous and next ranges
                    let mut should_skip = false;

                    if read_idx < ranges.len() - 1 {
                        let next_start = ranges[read_idx + 1].start;

                        // Check if current range is significantly overlapped by both neighbors
                        if curr_start > prev_start && next_start < curr_end {
                            let overlap_with_prev = prev_end - curr_start;
                            let overlap_with_next = curr_end - next_start;
                            let range_length = curr_end - curr_start;

                            // Skip if the range is mostly covered by its neighbors
                            if overlap_with_prev + overlap_with_next > range_length {
                                should_skip = true;
                            }
                        }
                    }

                    if !should_skip {
                        write_idx += 1;
                        if write_idx != read_idx {
                            ranges.swap(write_idx, read_idx);
                        }
                    }
                } else {
                    // No overlap - keep both ranges
                    write_idx += 1;
                    if write_idx != read_idx {
                        ranges.swap(write_idx, read_idx);
                    }
                }
            }
            ranges.truncate(write_idx + 1);

            // Create expected ranges
            let expected: Vec<RangeInfo> = expected_ranges
                .iter()
                .map(|(start, end)| create_range_info(*start, *end))
                .collect();

            // Compare results
            assert_eq!(
                ranges.len(),
                expected.len(),
                "Test case '{case_name}': Wrong number of ranges after containment removal"
            );

            for (i, (result, expected)) in ranges.iter().zip(expected.iter()).enumerate() {
                assert_eq!(
                    (result.start, result.end),
                    (expected.start, expected.end),
                    "Test case '{case_name}': Mismatch at position {i}"
                );
            }

            println!("Test case {} passed: {}", case_index + 1, case_name);
        }
    }
}
