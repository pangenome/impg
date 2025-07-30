use crate::sequence_index::UnifiedSequenceIndex;
use handlegraph::handle::{Handle, NodeId};
use log::{debug, error, info, warn};
use niffler::compression::Format;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
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
                "Unsupported compression format '{}', using none",
                compress_arg
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

    fn get_sequence(&mut self, idx: usize) -> io::Result<Vec<u8>> {
        if idx >= self.offsets.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid sequence index",
            ));
        }

        let (offset, length) = self.offsets[idx];
        let mut buffer = vec![0u8; length as usize];

        self.sequences_file.seek(SeekFrom::Start(offset))?;
        self.sequences_file.read_exact(&mut buffer)?;

        Ok(buffer)
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

    fn get_sequence(&mut self, handle: Handle) -> io::Result<Vec<u8>> {
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
    //gfa_id: usize,          // GFA file ID this range belongs to
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
pub fn run_lace(
    gfa_files: Option<Vec<String>>,
    gfa_list: Option<String>,
    output: &str,
    compress: &str,
    fill_gaps: u8,
    temp_dir: Option<String>,
    sequence_index: Option<&UnifiedSequenceIndex>,
    verbose: u8,
) -> io::Result<()> {
    info!("Running lace command");

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
    if let Some(seq_index) = sequence_index {
        match seq_index {
            UnifiedSequenceIndex::Fasta(fasta_index) => {
                info!(
                    "Using FASTA sequence index with {} files for gap filling",
                    fasta_index.fasta_paths.len()
                );
            }
            #[cfg(feature = "agc")]
            UnifiedSequenceIndex::Agc(agc_index) => {
                info!(
                    "Using AGC sequence index with {} files for gap filling",
                    agc_index.agc_paths.len()
                );
            }
        }
    } else if fill_gaps > 0 {
        warn!("Gap filling requested but no sequence files provided");
    }

    // Validate temp directory
    if let Some(ref temp_path) = temp_dir {
        if !Path::new(temp_path).exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("Temporary directory '{}' does not exist", temp_path),
            ));
        }
    }

    // log_memory_usage("start");

    // Create a single combined graph without paths and a map of path key to ranges
    info!("Collecting metadata from {} GFA files", gfa_files.len());
    let (combined_graph, mut path_key_ranges) =
        read_gfa_files(&gfa_files, temp_dir.as_deref())?;

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
            "Successfully wrote the combined graph to {} ({:?} format)",
            output, compression_format
        ),
        Err(e) => error!("Error writing the GFA file: {}", e),
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
                        format!("GFA file '{}' not found", file),
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
                            format!("GFA file '{}' not found", trimmed),
                        ));
                    }
                    files.push(trimmed.to_string());
                }
            }

            if files.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("No valid GFA files found in list file: {}", list_file),
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
                                error!("Invalid S line: {}", line);
                                continue;
                            }
                            let node_id: u64 = match fields[1].parse() {
                                Ok(id) => id,
                                Err(_) => {
                                    error!("Invalid node ID {} in line {}", fields[1], line);
                                    continue;
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
                                error!("Invalid L line: {}", line);
                                continue;
                            }

                            let from_id: u64 = match fields[1].parse() {
                                Ok(id) => id,
                                Err(_) => {
                                    error!("Invalid from node ID {} in line {}", fields[1], line);
                                    continue;
                                }
                            };
                            let from_rev = fields[2] == "-";
                            let to_id: u64 = match fields[3].parse() {
                                Ok(id) => id,
                                Err(_) => {
                                    error!("Invalid to node ID {} in line {}", fields[3], line);
                                    continue;
                                }
                            };
                            let to_rev = fields[4] == "-";
                            temp_edges.push(CompactEdge::new(from_id, from_rev, to_id, to_rev));
                            edge_count += 1;
                        }
                        "P" => {
                            if fields.len() < 3 {
                                error!("Invalid P line: {}", line);
                                continue;
                            }
                            let path_name = fields[1];
                            let nodes_str = fields[2];

                            if let Some((sample_hap_name, start, end)) = split_path_name(path_name)
                            {
                                let mut translated_steps = Vec::new();
                                for step_str in nodes_str.split(',') {
                                    if step_str.is_empty() {
                                        error!("Empty step in path {} in line {}", path_name, line);
                                        continue;
                                    }
                                    let (node_str, orient) = if let Some(stripped) =
                                        step_str.strip_suffix('+')
                                    {
                                        (stripped, false)
                                    } else if let Some(stripped) = step_str.strip_suffix('-') {
                                        (stripped, true)
                                    } else {
                                        error!("Invalid step format {} in line {}", step_str, line);
                                        continue;
                                    };

                                    let node_id: u64 = match node_str.parse() {
                                        Ok(id) => id,
                                        Err(_) => {
                                            error!(
                                                "Invalid node ID in path {} in line {}",
                                                path_name, line
                                            );
                                            continue;
                                        }
                                    };

                                    if let Some(&translated_id) =
                                        id_translation.get(&NodeId::from(node_id))
                                    {
                                        translated_steps.push(Handle::pack(translated_id, orient));
                                    } else {
                                        error!(
                                            "Node {} in path {} not found in translation map",
                                            node_id, path_name
                                        );
                                        continue;
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
                    "GFA file {} ({}) processed: {} nodes, {} edges",
                    gfa_id, gfa_path, node_count, edge_count
                );
            } else {
                error!("Failed to open GFA file '{}'", gfa_path);
            }
        });

    let path_map = Arc::try_unwrap(path_key_ranges)
        .expect("More than one Arc pointer to path map")
        .into_inner()
        .unwrap();

    let graph = Arc::try_unwrap(combined_graph)
        .ok()
        .expect("More than one Arc pointer to graph")
        .into_inner()
        .unwrap();

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
            format!("Failed to open file '{}': {}", gfa_path, e),
        )
    })?;

    let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
        io::Error::other(format!("Failed to open reader for '{}': {}", gfa_path, e))
    })?;

    Ok(Box::new(BufReader::new(reader)))
}

fn split_path_name(path_name: &str) -> Option<(String, usize, usize)> {
    if let Some(last_colon) = path_name.rfind(':') {
        let (key, range_str) = path_name.split_at(last_colon);
        let range_str = &range_str[1..];

        if let Some((start_str, end_str)) = range_str.split_once('-') {
            if let (Ok(start), Ok(end)) = (start_str.parse(), end_str.parse()) {
                return Some((key.to_string(), start, end));
            }
        }
    }
    None
}

fn sort_and_filter_ranges(ranges: &mut Vec<RangeInfo>) {
    ranges.sort_by_key(|r| (r.start, r.end));

    debug!("  Removing redundant ranges");

    let mut write_idx = 0;
    for read_idx in 1..ranges.len() {
        let (prev_start, prev_end) = (ranges[write_idx].start, ranges[write_idx].end);
        let (curr_start, curr_end) = (ranges[read_idx].start, ranges[read_idx].end);

        if curr_start == prev_start && curr_end == prev_end {
            continue;
        } else if curr_start >= prev_start && curr_end <= prev_end {
            debug!(
                "    Contained range [{}, {}) within [{}, {})",
                curr_start, curr_end, prev_start, prev_end
            );
            continue;
        } else {
            write_idx += 1;
            if write_idx != read_idx {
                ranges[write_idx] = ranges[read_idx].clone();
            }
        }
    }

    ranges.truncate(write_idx + 1);
}

fn trim_range_overlaps(ranges: &mut Vec<RangeInfo>, graph_mutex: &Arc<Mutex<CompactGraph>>) {
    let mut i = 0;
    while i < ranges.len() {
        let mut j = i + 1;
        while j < ranges.len() {
            if ranges[i].overlaps_with(&ranges[j]) {
                let overlap_start = ranges[i].start.max(ranges[j].start);
                let overlap_end = ranges[i].end.min(ranges[j].end);
                let overlap_size = overlap_end - overlap_start;

                debug!(
                    "    Trimming overlap [{}, {}) of size {} between ranges [{}, {}) and [{}, {})",
                    overlap_start,
                    overlap_end,
                    overlap_size,
                    ranges[i].start,
                    ranges[i].end,
                    ranges[j].start,
                    ranges[j].end
                );

                // Use split_at_mut to avoid double borrow
                if i < j {
                    let (left, right) = ranges.split_at_mut(j);
                    split_overlapping_ranges(&mut left[i], &mut right[0], graph_mutex);
                }
            }
            j += 1;
        }
        i += 1;
    }
}

fn split_overlapping_ranges(
    range1: &mut RangeInfo,
    range2: &mut RangeInfo,
    graph_mutex: &Arc<Mutex<CompactGraph>>,
) {
    // Simplified overlap handling - just truncate the later range
    if range1.start < range2.start {
        if range1.end > range2.start {
            range1.end = range2.start;
            // Adjust steps accordingly (simplified)
            while let Some(last_step) = range1.steps.pop() {
                let seq_len = {
                    let mut graph = graph_mutex.lock().unwrap();
                    graph.get_sequence(last_step).unwrap_or_default().len()
                };
                if range1.end + seq_len <= range2.start {
                    range1.steps.push(last_step);
                    break;
                }
            }
        }
    }
}

fn link_contiguous_ranges(ranges: &mut Vec<RangeInfo>, graph_mutex: &Arc<Mutex<CompactGraph>>) {
    for i in 0..ranges.len().saturating_sub(1) {
        if ranges[i].is_contiguous_with(&ranges[i + 1]) {
            if let (Some(&last_handle), Some(&first_handle)) =
                (ranges[i].steps.last(), ranges[i + 1].steps.first())
            {
                let mut graph = graph_mutex.lock().unwrap();
                if !graph.has_edge(
                    last_handle.id().into(),
                    last_handle.is_reverse(),
                    first_handle.id().into(),
                    first_handle.is_reverse(),
                ) {
                    graph.add_edge(
                        last_handle.id().into(),
                        last_handle.is_reverse(),
                        first_handle.id().into(),
                        first_handle.is_reverse(),
                    );
                    debug!(
                        "    Added link between nodes {} and {}",
                        last_handle.id(),
                        first_handle.id()
                    );
                }
            }
        }
    }
}

fn write_graph_to_gfa(
    combined_graph: &mut CompactGraph,
    path_key_ranges: &FxHashMap<String, Vec<RangeInfo>>,
    output_path: &str,
    compression_format: Format,
    fill_gaps: u8,
    sequence_index: Option<&UnifiedSequenceIndex>,
    _debug: bool,
) -> std::io::Result<()> {
    // Create the output file
    let output_file = File::create(output_path)?;

    // Create writer based on compression format
    let writer: Box<dyn Write> = match compression_format {
        Format::Gzip => {
            use flate2::write::GzEncoder;
            use flate2::Compression;
            Box::new(GzEncoder::new(output_file, Compression::default()))
        }
        Format::Zstd => {
            let encoder = ZstdEncoder::new(output_file, 6)?;
            Box::new(encoder)
        }
        Format::No => {
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

    // Write all nodes
    for node_id in 1..=combined_graph.node_count {
        let sequence = combined_graph.get_sequence(Handle::pack(NodeId::from(node_id), false))?;
        let sequence_str = String::from_utf8(sequence).expect("Node sequence contains invalid UTF-8");
        writeln!(file, "S\t{}\t{}", node_id, sequence_str)?;
    }

    // Write all edges
    for edge in &combined_graph.edges {
        let from_orient = if edge.from_rev() { "-" } else { "+" };
        let to_orient = if edge.to_rev() { "-" } else { "+" };
        writeln!(
            file,
            "L\t{}\t{}\t{}\t{}\t0M",
            edge.from_id(), from_orient, edge.to_id(), to_orient
        )?;
    }

    // Write paths by processing ranges directly
    for (path_key, ranges) in path_key_ranges {
        for range in ranges {
            if !range.steps.is_empty() {
                let path_name = format!("{}:{}-{}", path_key, range.start, range.end);
                let steps_str = range
                    .steps
                    .iter()
                    .map(|handle| {
                        format!(
                            "{}{}",
                            handle.id(),
                            if handle.is_reverse() { "-" } else { "+" }
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(",");

                if fill_gaps > 0 && sequence_index.is_some() {
                    fill_gaps_in_path(
                        &mut file,
                        &path_name,
                        &steps_str,
                        range,
                        sequence_index.unwrap(),
                        fill_gaps,
                        combined_graph,
                    )?;
                } else {
                    writeln!(file, "P\t{}\t{}\t*", path_name, steps_str)?;
                }
            }
        }
    }

    file.flush()?;
    Ok(())
}

fn fill_gaps_in_path(
    writer: &mut BufWriter<Box<dyn Write>>,
    path_name: &str,
    steps_str: &str,
    _range: &RangeInfo,
    _sequence_index: &UnifiedSequenceIndex,
    fill_gaps: u8,
    _graph: &mut CompactGraph,
) -> io::Result<()> {
    // Simplified gap filling - just write the path as-is for now
    // In a full implementation, this would fetch sequences and fill gaps
    writeln!(writer, "P\t{}\t{}\t*", path_name, steps_str)?;
    
    if fill_gaps > 0 {
        // Gap filling logic would go here
        debug!("Gap filling requested but not fully implemented yet");
    }
    
    Ok(())
}
