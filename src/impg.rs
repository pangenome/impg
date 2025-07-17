use crate::forest_map::ForestMap;
use crate::graph::reverse_complement;
use crate::paf::{ParseErr, PartialPafRecord, Strand};
use crate::seqidx::SequenceIndex;
use crate::sequence_index::SequenceIndex as _; // The as _ syntax imports the trait so its methods are available, but doesn't bring the name into scope (avoiding the naming conflict)
use crate::sequence_index::UnifiedSequenceIndex;
use coitrees::{BasicCOITree, Interval, IntervalTree};
use lib_tracepoints::variable_tracepoints_to_cigar_with_aligner;
use lib_wfa2::affine_wavefront::AffineWavefronts;
use log::debug;
use noodles::bgzf;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::cmp::{max, min};
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::sync::Arc;
use std::sync::RwLock;
use std::cell::RefCell;

thread_local! {
    static ALIGNERS: RefCell<Vec<AffineWavefronts>> = RefCell::new(Vec::new());
}

/// Parse a CIGAR string into a vector of CigarOp
// Note that the query_delta is negative for reverse strand alignments
#[derive(Debug, Clone, PartialEq)]
pub struct CigarOp {
    val: u32,
}

impl CigarOp {
    pub fn new(len: i32, op: char) -> Self {
        let val = match op {
            '=' => 0,
            'X' => 1,
            'I' => 2,
            'D' => 3,
            'M' => 4,
            _ => panic!("Invalid CIGAR operation: {}", op),
        };
        Self {
            val: (val << 29) | (len as u32),
        }
    }

    pub fn op(&self) -> char {
        // two most significant bits in the val tell us the op
        match self.val >> 29 {
            0 => '=',
            1 => 'X',
            2 => 'I',
            3 => 'D',
            4 => 'M',
            _ => panic!("Invalid CIGAR operation: {}", self.val >> 29),
        }
    }

    pub fn len(&self) -> i32 {
        (self.val & ((1 << 29) - 1)) as i32
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn target_delta(&self) -> i32 {
        match self.op() {
            '=' | 'X' | 'D' | 'M' => self.len(),
            'I' => 0,
            _ => panic!("Invalid CIGAR operation: {}", self.op()),
        }
    }

    pub fn query_delta(&self, strand: Strand) -> i32 {
        match self.op() {
            '=' | 'X' | 'I' | 'M' => {
                if strand == Strand::Forward {
                    self.len()
                } else {
                    -self.len()
                }
            }
            'D' => 0,
            _ => panic!("Invalid CIGAR operation: {}", self.op()),
        }
    }

    fn adjust_len(&mut self, length_delta: i32) {
        self.val = (self.val & (7 << 29)) | ((self.len() + length_delta) as u32);
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct QueryMetadata {
    query_id: u32,
    target_start: i32,
    target_end: i32,
    query_start: i32,
    query_end: i32,
    paf_file_index: u16,
    strand_and_data_offset: u64, // Track strand and cigar/tracepoints offset
    data_bytes: usize,
}

impl QueryMetadata {
    // Constants for bit manipulation
    const STRAND_BIT: u64 = 0x8000000000000000; // Most significant bit for u64

    fn strand(&self) -> Strand {
        if (self.strand_and_data_offset & Self::STRAND_BIT) != 0 {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }

    fn data_offset(&self) -> u64 {
        self.strand_and_data_offset & !Self::STRAND_BIT
    }

    fn get_data_bytes(
        &self,
        paf_files: &[String],
        paf_gzi_indices: &[Option<bgzf::gzi::Index>],
    ) -> Vec<u8> {
        // Allocate space for cigar
        let mut data_buffer = vec![0; self.data_bytes];

        // Get the correct PAF file
        let paf_file_index = self.paf_file_index as usize;
        let paf_file = &paf_files[paf_file_index];

        // Get reader and seek start of cigar str
        if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
            // Get the GZI index for the PAF file
            let paf_gzi_index = paf_gzi_indices.get(paf_file_index).and_then(Option::as_ref);

            let mut reader = bgzf::io::Reader::new(File::open(paf_file).unwrap());
            reader
                .seek_by_uncompressed_position(paf_gzi_index.unwrap(), self.data_offset())
                .unwrap();
            reader.read_exact(&mut data_buffer).unwrap();
        } else {
            let mut reader = File::open(paf_file).unwrap();
            reader.seek(SeekFrom::Start(self.data_offset())).unwrap();
            reader.read_exact(&mut data_buffer).unwrap();
        };

        data_buffer
    }

    fn get_cigar_ops_from_bytes(&self, data_bytes: Vec<u8>) -> Vec<CigarOp> {
        let cigar_str = std::str::from_utf8(&data_bytes).unwrap();
        parse_cigar_to_delta(cigar_str).ok().unwrap_or_default()
    }

    fn get_tracepoints_from_bytes(&self, data_bytes: Vec<u8>) -> Vec<(usize, Option<usize>)> {
        let tracepoints_str = std::str::from_utf8(&data_bytes).unwrap();
        parse_variable_tracepoints(tracepoints_str)
            .ok()
            .unwrap_or_default()
    }

    /// Fast check to determine if data contains tracepoints or CIGAR
    /// Tracepoints: "123,456;789,012" or "123;456;789" or "123" or "123,456"
    /// CIGAR: "123M45I67D"
    fn is_tracepoints_data(data_bytes: &[u8]) -> bool {
        // Find the first non-digit character after skipping initial digits
        let mut i = 0;
        // Skip leading digits
        while i < data_bytes.len() && data_bytes[i].is_ascii_digit() {
            i += 1;
        }

        // Check the first non-digit character
        if i < data_bytes.len() {
            let c = data_bytes[i];
            // Tracepoints use comma or semicolon as separators
            // CIGAR uses operation characters (M, I, D, =, X)
            c == b',' || c == b';'
        } else {
            // If it's all digits, it's a tracepoint (single number)
            // CIGAR strings always have operation characters
            true
        }
    }
}

pub type AdjustedInterval = (Interval<u32>, Vec<CigarOp>, Interval<u32>);
type TreeMap = FxHashMap<u32, Arc<BasicCOITree<QueryMetadata, u32>>>;

#[derive(Serialize, Deserialize)]
struct SerializableInterval {
    first: i32,
    last: i32,
    metadata: QueryMetadata,
}

#[derive(Default, Clone)]
pub struct SortedRanges {
    pub ranges: Vec<(i32, i32)>,
    sequence_length: i32,
    min_distance: i32,
}

impl SortedRanges {
    pub fn new(sequence_length: i32, min_distance: i32) -> Self {
        Self {
            ranges: Vec::new(),
            sequence_length,
            min_distance,
        }
    }

    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, (i32, i32)> {
        self.ranges.iter()
    }

    pub fn insert(&mut self, new_range: (i32, i32)) -> Vec<(i32, i32)> {
        let (mut start, mut end) = if new_range.0 <= new_range.1 {
            (new_range.0, new_range.1)
        } else {
            (new_range.1, new_range.0)
        };

        // Find nearby ranges to potentially merge with
        let mut i = match self.ranges.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(pos) => pos,
            Err(pos) => pos,
        };

        // Check previous range
        if i > 0 && (start - self.ranges[i - 1].1).abs() < self.min_distance {
            start = self.ranges[i - 1].1;
            i -= 1;
        } else if start < self.min_distance {
            start = 0;
        }

        // Check next range
        if i < self.ranges.len() && (self.ranges[i].0 - end).abs() < self.min_distance {
            end = self.ranges[i].0;
        } else if end > (self.sequence_length - self.min_distance) {
            end = self.sequence_length;
        }

        // Return regions that don't overlap with existing ranges
        let mut non_overlapping = Vec::new();
        let mut current = start;

        // Find the first range that could overlap
        let mut i = match self.ranges.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(pos) => pos,
            Err(pos) => pos,
        };

        // Check previous range for overlap
        if i > 0 && self.ranges[i - 1].1 > start {
            i -= 1;
        }

        // Process all potentially overlapping ranges
        while i < self.ranges.len() && current < end {
            let (range_start, range_end) = self.ranges[i];
            if range_start > end {
                break;
            }
            if current < range_start {
                non_overlapping.push((current, range_start));
            }
            current = max(current, range_end);
            i += 1;
        }

        if current < end {
            non_overlapping.push((current, end));
        }

        match self.ranges.binary_search_by_key(&start, |&(s, _)| s) {
            Ok(pos) | Err(pos) => {
                // Check if we can merge with the previous range
                if pos > 0 && self.ranges[pos - 1].1 >= start {
                    self.ranges[pos - 1].1 = max(self.ranges[pos - 1].1, end);
                    self.merge_forward_from(pos - 1);
                } else if pos < self.ranges.len() && end >= self.ranges[pos].0 {
                    self.ranges[pos].0 = min(start, self.ranges[pos].0);
                    self.ranges[pos].1 = max(end, self.ranges[pos].1);
                    self.merge_forward_from(pos);
                } else {
                    self.ranges.insert(pos, (start, end));
                }
            }
        }

        // let adjusted_new_range = if new_range.0 <= new_range.1 {
        //     (start, end)
        // } else {
        //     (end, start)
        // };

        non_overlapping
    }

    pub fn merge_forward_from(&mut self, start_idx: usize) {
        let mut write = start_idx;
        let mut read = start_idx + 1;
        while read < self.ranges.len() {
            if self.ranges[write].1 >= self.ranges[read].0 {
                self.ranges[write].1 = max(self.ranges[write].1, self.ranges[read].1);
            } else {
                write += 1;
                self.ranges.swap(write, read);
            }
            read += 1;
        }
        self.ranges.truncate(write + 1);
    }
}

/// Execute a closure with a thread-local aligner
fn with_thread_aligner<F, R>(penalties: (u8, u8, u8, u8, u8, u8), f: F) -> R
where
    F: FnOnce(&mut AffineWavefronts) -> R,
{
    ALIGNERS.with(|aligners| {
        let mut aligners = aligners.borrow_mut();
        if aligners.is_empty() {
            let (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = penalties;
            aligners.push(AffineWavefronts::with_penalties_affine2p(
                0, // match score
                mismatch.into(),
                gap_open1.into(),
                gap_ext1.into(),
                gap_open2.into(),
                gap_ext2.into(),
            ));
        }
        f(&mut aligners[0])
    })
}

pub struct Impg {
    pub trees: RwLock<TreeMap>,
    pub seq_index: SequenceIndex,
    paf_files: Vec<String>,                         // List of all PAF files
    paf_gzi_indices: Vec<Option<bgzf::gzi::Index>>, // Corresponding GZI indices
    pub forest_map: ForestMap,                      // Forest map for lazy loading
    index_file_path: String,                        // Path to the index file for lazy loading
}

impl Impg {
    fn process_tracepoints_data(
        &self,
        data_buffer: Vec<u8>,
        metadata: &QueryMetadata,
        target_id: u32,
        sequence_index: Option<&UnifiedSequenceIndex>,
        penalties: (u8, u8, u8, u8, u8, u8),
        requested_range: (i32, i32),
    ) -> (i32, i32, i32, i32, Vec<CigarOp>) {
        if let Some(sequence_index) = sequence_index {
            // Get the tracepoints
            let tracepoints = metadata.get_tracepoints_from_bytes(data_buffer);

            // Scan tracepoints to find the relevant subset for the requested range
            let mut query_pos = if metadata.strand() == Strand::Forward {
                metadata.query_start
            } else {
                metadata.query_end
            };
            let mut target_pos = metadata.target_start;
            let mut start_idx = 0;
            let mut end_idx = tracepoints.len();
            let mut found_start = false;
            let mut subset_query_start = if metadata.strand() == Strand::Forward {
                metadata.query_start
            } else {
                metadata.query_end
            };
            let mut subset_target_start = metadata.target_start;
            let mut subset_query_end = if metadata.strand() == Strand::Forward {
                metadata.query_end
            } else {
                metadata.query_start
            };
            let mut subset_target_end = metadata.target_end;

            // Find the relevant subset of tracepoints
            for (i, &(a_len, b_len_opt)) in tracepoints.iter().enumerate() {
                debug!("Tracepoint {}: {:?}, {:?}", i, a_len, b_len_opt);
                debug!(
                    "Current query_pos: {}, target_pos: {}",
                    query_pos, target_pos
                );

                let query_delta = if metadata.strand() == Strand::Forward {
                    a_len as i32
                } else {
                    -(a_len as i32)
                };
                let target_delta = b_len_opt.unwrap_or(a_len) as i32;
                
                let next_query_pos = query_pos + query_delta;
                let next_target_pos = target_pos + target_delta;
                
                // Check if this tracepoint overlaps with the requested range
                let overlaps = target_pos < requested_range.1 && next_target_pos >= requested_range.0;
                
                if overlaps && !found_start {
                    debug!("-------------->Found overlapping tracepoint: {:?}", tracepoints[i]);
                    
                    // Found the start of the relevant region
                    start_idx = i;
                    subset_query_start = query_pos;
                    subset_target_start = target_pos;
                    found_start = true;
                }
                
                if found_start {
                    // Update the end of the relevant region as we go
                    end_idx = i + 1;
                    subset_query_end = next_query_pos;
                    subset_target_end = next_target_pos;
                    
                    // Only break if the current tracepoint no longer overlaps
                    if target_pos >= requested_range.1 {
                        debug!("-------------->Found end of relevant tracepoint region: {:?}", tracepoints[i]);
                        break;
                    }
                }
                
                query_pos = next_query_pos;
                target_pos = next_target_pos;
            }

            debug!("Subsets query: {}-{}, target: {}-{}",
                subset_query_start, subset_query_end,
                subset_target_start, subset_target_end
            );

            // Extract the relevant subset of tracepoints
            let subset_tracepoints = &tracepoints[start_idx..end_idx];

            debug!("Tracepoints from {} to {} ({:?})",
                tracepoints.len(), subset_tracepoints.len(), subset_tracepoints
            );
            debug!(
                "From {}/{} -> {}/{} to {}/{} -> {}/{}",
                metadata.query_start,
                metadata.target_start,
                metadata.query_end,
                metadata.target_end,
                subset_query_start,
                subset_target_start,
                subset_query_end,
                subset_target_end
            );

            // Fetch only the relevant portions of the sequences
            let query_name = self.seq_index.get_name(metadata.query_id).unwrap();
            let query_seq_result = if metadata.strand() == Strand::Forward {
                sequence_index.fetch_sequence(query_name, subset_query_start, subset_query_end)
            } else {
                // For reverse strand, coordinates are inverted, so we need min/max
                let (fetch_start, fetch_end) = if subset_query_start < subset_query_end {
                    (subset_query_start, subset_query_end)
                } else {
                    (subset_query_end, subset_query_start)
                };
                // Fetch and reverse complement
                match sequence_index.fetch_sequence(query_name, fetch_start, fetch_end) {
                    Ok(seq) => Ok(reverse_complement(&seq)),
                    Err(e) => Err(e),
                }
            };

            // Fetch target sequence
            let target_name = self.seq_index.get_name(target_id).unwrap();
            let target_seq_result = sequence_index.fetch_sequence(target_name, subset_target_start, subset_target_end);

            // Convert tracepoints to CIGAR if we successfully fetched both sequences
            match (query_seq_result, target_seq_result) {
                (Ok(query_seq), Ok(target_seq)) => {
                    // Convert tracepoints to CIGAR operations using thread-local aligner
                    let cigar_str = with_thread_aligner(penalties, |aligner| {
                        variable_tracepoints_to_cigar_with_aligner(
                            subset_tracepoints,
                            &query_seq,
                            &target_seq,
                            0,
                            0,
                            aligner,
                        )
                    });

                    (
                        subset_target_start,
                        subset_target_end,
                        // For reverse strand, coordinates should be in forward orientation (start < end)
                        if metadata.strand() == Strand::Forward {
                            subset_query_start
                        } else {
                            subset_query_end
                        },
                        if metadata.strand() == Strand::Forward {
                            subset_query_end
                        } else {
                            subset_query_start
                        },
                        // Parse the CIGAR string to CigarOp vector
                        parse_cigar_to_delta(&cigar_str).ok().unwrap_or_default(),
                    )
                }
                (Err(e), _) => {
                    panic!("Failed to fetch query sequence: {}", e);
                }
                (_, Err(e)) => {
                    panic!("Failed to fetch target sequence: {}", e);
                }
            }
        } else {
            #[cfg(feature = "agc")]
            let file_types = "FASTA/AGC";
            #[cfg(not(feature = "agc"))]
            let file_types = "FASTA";
            panic!(
                "Sequence data ({}) is required for tracepoints conversion. Use --sequence-files or --sequence-list",
                file_types
            )
        }
    }

    pub fn from_multi_paf_records(
        records_by_file: &[(Vec<PartialPafRecord>, String)],
        seq_index: SequenceIndex,
    ) -> Result<Self, ParseErr> {
        // Use par_iter to process the files in parallel and collect both pieces of information
        let (paf_files, paf_gzi_indices): (Vec<String>, Vec<Option<bgzf::gzi::Index>>) =
            records_by_file
                .par_iter()
                .map(|(_, paf_file)| {
                    let paf_gzi_index = if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
                        let paf_gzi_file = paf_file.to_owned() + ".gzi";
                        Some(
                            bgzf::gzi::fs::read(paf_gzi_file.clone())
                                .unwrap_or_else(|_| panic!("Could not open {}", paf_gzi_file)),
                        )
                    } else {
                        None
                    };

                    // Return both values as a tuple
                    (paf_file.clone(), paf_gzi_index)
                })
                .unzip(); // Separate the tuples into two vectors

        let intervals: FxHashMap<u32, Vec<Interval<QueryMetadata>>> = records_by_file
            .par_iter()
            .enumerate() // Add enumeration to get the position as index
            .flat_map(|(file_index, (records, _))| {
                records
                    .par_iter()
                    .filter_map(|record| {
                        let query_metadata = QueryMetadata {
                            query_id: record.query_id,
                            target_start: record.target_start as i32,
                            target_end: record.target_end as i32,
                            query_start: record.query_start as i32,
                            query_end: record.query_end as i32,
                            paf_file_index: file_index as u16,
                            strand_and_data_offset: record.strand_and_data_offset, // Already includes strand bit
                            data_bytes: record.data_bytes,
                        };

                        Some((
                            record.target_id,
                            Interval {
                                first: record.target_start as i32,
                                last: record.target_end as i32,
                                metadata: query_metadata,
                            },
                        ))
                    })
                    .collect::<Vec<_>>()
            })
            .fold(
                FxHashMap::default,
                |mut acc: FxHashMap<u32, Vec<Interval<QueryMetadata>>>, (target_id, interval)| {
                    acc.entry(target_id).or_default().push(interval);
                    acc
                },
            )
            .reduce(FxHashMap::default, |mut acc, part| {
                for (key, value) in part {
                    acc.entry(key).or_default().extend(value);
                }
                acc
            });

        let trees: TreeMap = intervals
            .into_par_iter()
            .map(|(target_id, interval_nodes)| {
                (
                    target_id,
                    Arc::new(BasicCOITree::new(interval_nodes.as_slice())),
                )
            })
            .collect();

        Ok(Self {
            trees: RwLock::new(trees),
            seq_index,
            paf_files,
            paf_gzi_indices,
            forest_map: ForestMap::new(), // All trees are in memory, no need for forest map
            index_file_path: String::new(), // All trees are in memory, no need for index file path
        })
    }

    /// Serialize the IMPG index to a writer with embedded forest map
    pub fn serialize_with_forest_map<W: std::io::Write + std::io::Seek>(
        &self,
        mut writer: W,
    ) -> std::io::Result<()> {
        const MAGIC: &[u8] = b"IMPGIDX1";
        let mut forest_map = ForestMap::new();

        // Write magic bytes
        writer.write_all(MAGIC)?;

        // Write placeholder for forest map offset (will be updated later)
        writer.write_all(&[0u8; 8])?;

        let mut current_offset = 16u64; // Header size

        // Serialize sequence index
        let seq_index_data =
            bincode::serde::encode_to_vec(&self.seq_index, bincode::config::standard()).map_err(
                |e| std::io::Error::other(format!("Failed to encode sequence index: {:?}", e)),
            )?;

        writer.write_all(&seq_index_data)?;
        current_offset += seq_index_data.len() as u64;

        // Serialize each tree and track its offset
        let trees = self.trees.read().unwrap();
        for (&target_id, tree) in trees.iter() {
            // Record the offset before serializing this tree
            forest_map.add_entry(target_id, current_offset);

            // Convert tree to serializable format
            let intervals: Vec<SerializableInterval> = tree
                .iter()
                .map(|interval| SerializableInterval {
                    first: interval.first,
                    last: interval.last,
                    metadata: interval.metadata.clone(),
                })
                .collect();

            // Serialize the target_id and intervals
            let tree_data =
                bincode::serde::encode_to_vec(&(target_id, intervals), bincode::config::standard())
                    .map_err(|e| {
                        std::io::Error::other(format!(
                            "Failed to encode tree for target {}: {:?}",
                            target_id, e
                        ))
                    })?;

            writer.write_all(&tree_data)?;
            current_offset += tree_data.len() as u64;
        }

        // Now write the forest map at the end
        let forest_map_offset = current_offset;
        let forest_map_data =
            bincode::serde::encode_to_vec(&forest_map, bincode::config::standard()).map_err(
                |e| std::io::Error::other(format!("Failed to encode forest map: {:?}", e)),
            )?;

        writer.write_all(&forest_map_data)?;

        // Go back and update the forest map offset in the header
        writer.seek(SeekFrom::Start(8))?;
        writer.write_all(&forest_map_offset.to_le_bytes())?;

        Ok(())
    }

    /// Load a specific tree from disk using the forest map
    fn load_tree_from_disk(&self, target_id: u32) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        if let Some(tree_offset) = self.forest_map.get_tree_offset(target_id) {
            let mut file = File::open(&self.index_file_path).expect("Failed to open index file");
            file.seek(std::io::SeekFrom::Start(tree_offset))
                .expect("Failed to seek in index file");

            let mut reader = BufReader::new(file);
            let (loaded_target_id, intervals): (u32, Vec<SerializableInterval>) =
                bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                    .unwrap_or_else(|_| {
                        panic!("Failed to deserialize tree for target {}", target_id)
                    });

            // Verify we loaded the correct tree
            if loaded_target_id != target_id {
                panic!(
                    "Tree mismatch: expected {}, got {}",
                    target_id, loaded_target_id
                );
            }

            // Reconstruct the tree
            let tree = BasicCOITree::new(
                intervals
                    .into_iter()
                    .map(|interval| Interval {
                        first: interval.first,
                        last: interval.last,
                        metadata: interval.metadata, // Move instead of clone
                    })
                    .collect::<Vec<_>>()
                    .as_slice(),
            );

            let arc_tree = Arc::new(tree);

            // Cache the tree for future use
            self.trees
                .write()
                .unwrap()
                .insert(target_id, Arc::clone(&arc_tree));

            Some(arc_tree)
        } else {
            None // Tree not found in forest map
        }
    }

    /// Get a tree from memory or load it from disk if necessary
    pub fn get_or_load_tree(
        &self,
        target_id: u32,
    ) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        // First check if the tree is already in memory
        if let Some(tree) = self.trees.read().unwrap().get(&target_id) {
            // Get a clone of the Arc<tree> (incrementing the reference count, a cheap operation)
            // We clone to avoid holding the RwLock for the duration of the operation using the tree

            return Some(Arc::clone(tree));
        }

        // Not in memory - try to load from disk
        self.load_tree_from_disk(target_id)
    }

    /// Load IMPG index from the format with embedded forest map at the end
    pub fn load_from_file<R: std::io::Read + std::io::Seek>(
        mut reader: R,
        paf_files: &[String],
        index_file_path: String,
    ) -> std::io::Result<Self> {
        const MAGIC: &[u8] = b"IMPGIDX1";

        // Read and verify magic bytes
        let mut magic_buf = [0u8; 8];
        reader.read_exact(&mut magic_buf)?;
        if magic_buf != MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Invalid magic bytes - not a valid IMPG index file",
            ));
        }

        // Read forest map offset
        let mut offset_buf = [0u8; 8];
        reader.read_exact(&mut offset_buf)?;
        let forest_map_offset = u64::from_le_bytes(offset_buf);

        // Read sequence index
        let seq_index: SequenceIndex =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load sequence index: {:?}", e),
                    )
                })?;

        // Seek to forest map and read it
        reader.seek(SeekFrom::Start(forest_map_offset))?;
        let forest_map: ForestMap =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load forest map: {:?}", e),
                    )
                })?;

        // Determine PAF GZI indices
        let paf_gzi_indices = paf_files
            .iter()
            .map(|paf_file| {
                if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
                    let paf_gzi_file = format!("{}.gzi", paf_file);
                    Some(
                        bgzf::gzi::fs::read(paf_gzi_file.clone())
                            .unwrap_or_else(|_| panic!("Could not open {}", paf_gzi_file)),
                    )
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        Ok(Self {
            trees: RwLock::new(FxHashMap::default()), // Start with empty trees - load on demand
            seq_index,
            paf_files: paf_files.to_vec(),
            paf_gzi_indices,
            forest_map,
            index_file_path,
        })
    }

    pub fn query(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        penalties: Option<(u8, u8, u8, u8, u8, u8)>,
    ) -> Vec<AdjustedInterval> {
        let mut results = Vec::new();
        // Add the input range to the results
        results.push((
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
            if store_cigar {
                vec![CigarOp::new(range_end - range_start, '=')]
            } else {
                Vec::new()
            },
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
        ));

        debug!(
            "Querying region: {}:{}-{}, len: {}",
            self.seq_index.get_name(target_id).unwrap(),
            range_start,
            range_end,
            range_end - range_start
        );

        // Use default penalties if not provided
        let (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) =
            penalties.unwrap_or((0, 4, 6, 2, 26, 1));

        // Get or load the tree - if None, no overlaps exist for this target
        if let Some(tree) = self.get_or_load_tree(target_id) {
            tree.query(range_start, range_end, |interval| {
                let metadata = &interval.metadata;
                let data_buffer = metadata.get_data_bytes(&self.paf_files, &self.paf_gzi_indices);

                let (adj_target_start, adj_target_end, adj_query_start, adj_query_end, cigar_ops) =
                    if QueryMetadata::is_tracepoints_data(&data_buffer) {
                        self.process_tracepoints_data(
                            data_buffer,
                            metadata,
                            target_id,
                            sequence_index,
                            (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2),
                            (range_start, range_end),
                        )
                    } else {
                        // Handle regular CIGAR
                        (
                            metadata.target_start,
                            metadata.target_end,
                            metadata.query_start,
                            metadata.query_end,
                            metadata.get_cigar_ops_from_bytes(data_buffer),
                        )
                    };

                // Project through alignment
                let result = project_target_range_through_alignment(
                    (range_start, range_end),
                    (
                        adj_target_start,
                        adj_target_end,
                        adj_query_start,
                        adj_query_end,
                        metadata.strand(),
                    ),
                    &cigar_ops,
                );

                if let Some((
                    adjusted_query_start,
                    adjusted_query_end,
                    adjusted_cigar,
                    adjusted_target_start,
                    adjusted_target_end,
                )) = result
                {
                    // Check gap-compressed identity if threshold is provided
                    if let Some(threshold) = min_gap_compressed_identity {
                        if calculate_gap_compressed_identity(&adjusted_cigar) < threshold {
                            return; // Skip this result
                        }
                    }

                    let adjusted_interval = (
                        Interval {
                            first: adjusted_query_start,
                            last: adjusted_query_end,
                            metadata: metadata.query_id,
                        },
                        if store_cigar {
                            adjusted_cigar
                        } else {
                            Vec::new()
                        },
                        Interval {
                            first: adjusted_target_start,
                            last: adjusted_target_end,
                            metadata: target_id,
                        },
                    );
                    results.push(adjusted_interval);
                }
            });
        }

        results
    }

    pub fn query_transitive_dfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        penalties: Option<(u8, u8, u8, u8, u8, u8)>,
    ) -> Vec<AdjustedInterval> {
        // Initialize visited ranges from masked regions if provided
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter().map(|(&k, v)| (k, (*v).clone())).collect()
        } else {
            (0..self.seq_index.len() as u32)
                .into_par_iter() // Use parallel iterator
                .map(|id| {
                    let len = self.seq_index.get_len_from_id(id).unwrap();
                    (id, SortedRanges::new(len as i32, 0))
                })
                .collect()
        };

        // Filter input range
        let filtered_input_range = visited_ranges
            .entry(target_id)
            .or_default()
            .insert((range_start, range_end));

        let mut results = Vec::new();
        let mut stack = Vec::new();

        for (filtered_start, filtered_end) in filtered_input_range {
            // Add the filtered input range(s) to the results
            results.push((
                Interval {
                    first: filtered_start,
                    last: filtered_end,
                    metadata: target_id,
                },
                if store_cigar {
                    vec![CigarOp::new(filtered_end - filtered_start, '=')]
                } else {
                    Vec::new()
                },
                Interval {
                    first: filtered_start,
                    last: filtered_end,
                    metadata: target_id,
                },
            ));

            // Add the filtered input range(s) to the stack
            if (filtered_start - filtered_end).abs() >= min_transitive_len {
                stack.push((target_id, filtered_start, filtered_end, 0u16));
            }
        }

        while let Some((
            current_target_id,
            current_target_start,
            current_target_end,
            current_depth,
        )) = stack.pop()
        {
            // Check if we've reached max depth
            if max_depth > 0 && current_depth >= max_depth {
                continue;
            }

            debug!(
                "Querying region: {}:{}-{}, len: {}",
                self.seq_index.get_name(current_target_id).unwrap(),
                current_target_start,
                current_target_end,
                current_target_end - current_target_start
            );

            let prec_num_results = results.len();

            // Use default penalties if not provided
            let (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) =
                penalties.unwrap_or((0, 4, 6, 2, 26, 1));

            // Get or load the tree - if None, no overlaps exist for this target
            if let Some(tree) = self.get_or_load_tree(current_target_id) {
                // Collect intervals first to process them in parallel
                let mut intervals = Vec::new();

                tree.query(current_target_start, current_target_end, |interval| {
                    intervals.push(interval.metadata.clone());
                });

                // Process the intervals in parallel
                let processed_results: Vec<_> = intervals
                    .into_par_iter()
                    .filter_map(|metadata| {
                        let data_buffer =
                            metadata.get_data_bytes(&self.paf_files, &self.paf_gzi_indices);

                        let (
                            adj_target_start,
                            adj_target_end,
                            adj_query_start,
                            adj_query_end,
                            cigar_ops,
                        ) = if QueryMetadata::is_tracepoints_data(&data_buffer) {
                            self.process_tracepoints_data(
                                data_buffer,
                                &metadata,
                                current_target_id,
                                sequence_index,
                                (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2),
                                (current_target_start, current_target_end),
                            )
                        } else {
                            // Handle regular CIGAR
                            (
                                metadata.target_start,
                                metadata.target_end,
                                metadata.query_start,
                                metadata.query_end,
                                metadata.get_cigar_ops_from_bytes(data_buffer),
                            )
                        };

                        let result = project_target_range_through_alignment(
                            (current_target_start, current_target_end),
                            (
                                adj_target_start,
                                adj_target_end,
                                adj_query_start,
                                adj_query_end,
                                metadata.strand(),
                            ),
                            &cigar_ops,
                        );

                        if let Some((
                            adjusted_query_start,
                            adjusted_query_end,
                            adjusted_cigar,
                            adjusted_target_start,
                            adjusted_target_end,
                        )) = result
                        {
                            // Check gap-compressed identity if threshold is provided
                            if let Some(threshold) = min_gap_compressed_identity {
                                if calculate_gap_compressed_identity(&adjusted_cigar) < threshold {
                                    return None; // Skip this result
                                }
                            }

                            Some((
                                Interval {
                                    first: adjusted_query_start,
                                    last: adjusted_query_end,
                                    metadata: metadata.query_id,
                                },
                                if store_cigar {
                                    adjusted_cigar
                                } else {
                                    Vec::new()
                                },
                                Interval {
                                    first: adjusted_target_start,
                                    last: adjusted_target_end,
                                    metadata: current_target_id,
                                },
                                metadata.query_id,
                                adjusted_query_start,
                                adjusted_query_end,
                            ))
                        } else {
                            None
                        }
                    })
                    .collect();

                // Process results sequentially to maintain deterministic behavior
                for (
                    query_interval,
                    cigar,
                    target_interval,
                    query_id,
                    adjusted_query_start,
                    adjusted_query_end,
                ) in processed_results
                {
                    results.push((query_interval, cigar, target_interval));

                    // Only add non-overlapping portions to the stack for further exploration
                    if query_id != current_target_id {
                        let ranges = visited_ranges.entry(query_id).or_default();

                        let mut should_add = true;

                        // Check if the range is too close to any existing ranges
                        if min_distance_between_ranges > 0 {
                            let (new_min, new_max) = if adjusted_query_start <= adjusted_query_end {
                                (adjusted_query_start, adjusted_query_end)
                            } else {
                                (adjusted_query_end, adjusted_query_start)
                            };

                            // Find insertion point in sorted ranges
                            let idx = match ranges
                                .ranges
                                .binary_search_by_key(&new_min, |&(start, _)| start)
                            {
                                Ok(i) => i,
                                Err(i) => i,
                            };

                            // Only need to check adjacent ranges due to sorting
                            if idx > 0 {
                                // Check previous range
                                let (_, prev_end) = ranges.ranges[idx - 1];
                                if (new_min - prev_end).abs() < min_distance_between_ranges {
                                    should_add = false;
                                }
                            }
                            if idx < ranges.ranges.len() {
                                // Check next range
                                let (next_start, _) = ranges.ranges[idx];
                                if (next_start - new_max).abs() < min_distance_between_ranges {
                                    should_add = false;
                                }
                            }
                        }

                        if should_add {
                            let new_ranges =
                                ranges.insert((adjusted_query_start, adjusted_query_end));

                            // Add non-overlapping portions to stack
                            for (new_start, new_end) in new_ranges {
                                if (new_end - new_start).abs() >= min_transitive_len {
                                    stack.push((query_id, new_start, new_end, current_depth + 1));
                                }
                            }
                        }
                    }
                }
            }

            debug!("Collected {} results", results.len() - prec_num_results);

            // Merge contiguous/overlapping ranges with same sequence_id
            let stack_size = stack.len();
            stack.par_sort_by_key(|(id, start, _, _)| (*id, *start));

            let mut write = 0;
            for read in 1..stack.len() {
                if stack[write].0 == stack[read].0 &&   // Same sequence_id 
                    stack[write].2 >= stack[read].1
                // Overlapping or contiguous
                {
                    // Merge by extending end
                    stack[write].2 = stack[write].2.max(stack[read].2);
                } else {
                    write += 1;
                    stack.swap(write, read);
                }
            }
            stack.truncate(write + 1);
            debug!("Merged stack size from {} to {}", stack_size, stack.len());
        }

        results
    }

    pub fn query_transitive_bfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        penalties: Option<(u8, u8, u8, u8, u8, u8)>,
    ) -> Vec<AdjustedInterval> {
        // Initialize visited ranges from masked regions if provided
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter().map(|(&k, v)| (k, (*v).clone())).collect()
        } else {
            (0..self.seq_index.len() as u32)
                .into_par_iter() // Use parallel iterator
                .map(|id| {
                    let len = self.seq_index.get_len_from_id(id).unwrap();
                    (id, SortedRanges::new(len as i32, 0))
                })
                .collect()
        };

        // Filter input range
        let filtered_input_range = visited_ranges
            .entry(target_id)
            .or_default()
            .insert((range_start, range_end));

        let mut results = Vec::new();

        // Add the filtered input range(s) to the results
        for (filtered_start, filtered_end) in &filtered_input_range {
            results.push((
                Interval {
                    first: *filtered_start,
                    last: *filtered_end,
                    metadata: target_id,
                },
                if store_cigar {
                    vec![CigarOp::new(filtered_end - filtered_start, '=')]
                } else {
                    Vec::new()
                },
                Interval {
                    first: *filtered_start,
                    last: *filtered_end,
                    metadata: target_id,
                },
            ));
        }

        // Initialize ranges for first depth
        let mut current_depth = 0;
        let mut current_ranges = Vec::new();

        for (filtered_start, filtered_end) in filtered_input_range {
            if (filtered_start - filtered_end).abs() >= min_transitive_len {
                current_ranges.push((target_id, filtered_start, filtered_end));
            }
        }

        // Use default penalties if not provided
        let (_match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) =
            penalties.unwrap_or((0, 4, 6, 2, 26, 1));

        // Process by depth until max_depth or no more ranges
        while !current_ranges.is_empty() && (max_depth == 0 || current_depth < max_depth) {
            debug!(
                "Processing depth {} with {} ranges",
                current_depth,
                current_ranges.len()
            );

            // Process current depth ranges in parallel
            let query_results: Vec<Vec<(u32, i32, i32, Vec<CigarOp>, i32, i32, u32)>> =
                current_ranges
                    .par_iter()
                    .map(
                        |(current_target_id, current_target_start, current_target_end)| {
                            let mut local_results = Vec::new();

                            // Get or load the tree - if None, no overlaps exist for this target
                            if let Some(tree) = self.get_or_load_tree(*current_target_id) {
                                tree.query(
                                    *current_target_start,
                                    *current_target_end,
                                    |interval| {
                                        let metadata = &interval.metadata;
                                        let data_buffer = metadata
                                            .get_data_bytes(&self.paf_files, &self.paf_gzi_indices);

                                        let (
                                            adj_target_start,
                                            adj_target_end,
                                            adj_query_start,
                                            adj_query_end,
                                            cigar_ops,
                                        ) = if QueryMetadata::is_tracepoints_data(&data_buffer) {
                                            self.process_tracepoints_data(
                                                data_buffer,
                                                metadata,
                                                *current_target_id,
                                                sequence_index,
                                                (
                                                    _match, mismatch, gap_open1, gap_ext1,
                                                    gap_open2, gap_ext2,
                                                ),
                                                (*current_target_start, *current_target_end),
                                            )
                                        } else {
                                            // Handle regular CIGAR
                                            (
                                                metadata.target_start,
                                                metadata.target_end,
                                                metadata.query_start,
                                                metadata.query_end,
                                                metadata.get_cigar_ops_from_bytes(data_buffer),
                                            )
                                        };

                                        let result = project_target_range_through_alignment(
                                            (*current_target_start, *current_target_end),
                                            (
                                                adj_target_start,
                                                adj_target_end,
                                                adj_query_start,
                                                adj_query_end,
                                                metadata.strand(),
                                            ),
                                            &cigar_ops,
                                        );

                                        if let Some((
                                            adjusted_query_start,
                                            adjusted_query_end,
                                            adjusted_cigar,
                                            adjusted_target_start,
                                            adjusted_target_end,
                                        )) = result
                                        {
                                            // Check gap-compressed identity if threshold is provided
                                            if let Some(threshold) = min_gap_compressed_identity {
                                                if calculate_gap_compressed_identity(
                                                    &adjusted_cigar,
                                                ) < threshold
                                                {
                                                    return; // Skip this result
                                                }
                                            }

                                            local_results.push((
                                                metadata.query_id,
                                                adjusted_query_start,
                                                adjusted_query_end,
                                                if store_cigar {
                                                    adjusted_cigar
                                                } else {
                                                    Vec::new()
                                                },
                                                adjusted_target_start,
                                                adjusted_target_end,
                                                *current_target_id,
                                            ));
                                        }
                                    },
                                );
                            }

                            local_results
                        },
                    )
                    .collect();

            // Prepare for next depth
            let mut next_depth_ranges = Vec::new();

            // Process results sequentially to update visited_ranges and results
            for query_result in query_results {
                for (
                    query_id,
                    adjusted_query_start,
                    adjusted_query_end,
                    adjusted_cigar,
                    adjusted_target_start,
                    adjusted_target_end,
                    current_target_id,
                ) in query_result
                {
                    // Add to results
                    results.push((
                        Interval {
                            first: adjusted_query_start,
                            last: adjusted_query_end,
                            metadata: query_id,
                        },
                        adjusted_cigar,
                        Interval {
                            first: adjusted_target_start,
                            last: adjusted_target_end,
                            metadata: current_target_id,
                        },
                    ));

                    // Only consider for next depth if it's a different sequence
                    if query_id != current_target_id {
                        let ranges = visited_ranges.entry(query_id).or_default();

                        let mut should_add = true;

                        // Check proximity to existing ranges
                        if min_distance_between_ranges > 0 {
                            let (new_min, new_max) = if adjusted_query_start <= adjusted_query_end {
                                (adjusted_query_start, adjusted_query_end)
                            } else {
                                (adjusted_query_end, adjusted_query_start)
                            };

                            // Find insertion point in sorted ranges
                            let idx = match ranges
                                .ranges
                                .binary_search_by_key(&new_min, |&(start, _)| start)
                            {
                                Ok(i) => i,
                                Err(i) => i,
                            };

                            // Only need to check adjacent ranges due to sorting
                            if idx > 0 {
                                // Check previous range
                                let (_, prev_end) = ranges.ranges[idx - 1];
                                if (new_min - prev_end).abs() < min_distance_between_ranges {
                                    should_add = false;
                                }
                            }

                            if should_add && idx < ranges.ranges.len() {
                                // Check next range
                                let (next_start, _) = ranges.ranges[idx];
                                if (next_start - new_max).abs() < min_distance_between_ranges {
                                    should_add = false;
                                }
                            }
                        }

                        if should_add {
                            let new_ranges =
                                ranges.insert((adjusted_query_start, adjusted_query_end));

                            // Add non-overlapping portions to next depth
                            for (new_start, new_end) in new_ranges {
                                if (new_end - new_start).abs() >= min_transitive_len {
                                    next_depth_ranges.push((query_id, new_start, new_end));
                                }
                            }
                        }
                    }
                }
            }

            // Move to next depth
            current_depth += 1;

            // Prepare ranges for next depth
            if !next_depth_ranges.is_empty() {
                // Sort and merge contiguous/overlapping ranges
                next_depth_ranges.par_sort_by_key(|(id, start, _)| (*id, *start));

                let mut write = 0;
                for read in 1..next_depth_ranges.len() {
                    if next_depth_ranges[write].0 == next_depth_ranges[read].0 &&  // Same sequence_id 
                       next_depth_ranges[write].2 >= next_depth_ranges[read].1
                    // Overlapping or contiguous
                    {
                        // Merge by extending end
                        next_depth_ranges[write].2 =
                            next_depth_ranges[write].2.max(next_depth_ranges[read].2);
                    } else {
                        write += 1;
                        next_depth_ranges.swap(write, read);
                    }
                }
                next_depth_ranges.truncate(write + 1);

                debug!(
                    "Next depth will process {} ranges after merging",
                    next_depth_ranges.len()
                );
            }

            // Set up for next iteration
            current_ranges = next_depth_ranges;
        }

        results
    }
}

fn project_target_range_through_alignment(
    requested_target_range: (i32, i32),
    record: (i32, i32, i32, i32, Strand),
    cigar_ops: &[CigarOp],
) -> Option<(i32, i32, Vec<CigarOp>, i32, i32)> {
    let (target_start, target_end, query_start, query_end, strand) = record;

    let dir = if strand == Strand::Forward { 1 } else { -1 };
    let mut query_pos = if strand == Strand::Forward {
        query_start
    } else {
        query_end
    };
    let mut target_pos = target_start;

    // Track CIGAR slice bounds
    let mut first_op_idx = 0;
    let mut last_op_idx = 0;
    let mut found_overlap = false;

    let mut projected_query_start = -1;
    let mut projected_query_end = -1;
    let mut projected_target_start = -1;
    let mut projected_target_end = -1;
    let mut first_op_offset = 0;
    let mut last_op_remaining = 0;

    // Calculate the last valid target position
    let last_target_pos = min(target_end, requested_target_range.1);

    for (curr_op_idx, cigar_op) in cigar_ops.iter().enumerate() {
        // If the target position is past the end of the range, we can stop
        if target_pos > last_target_pos {
            break;
        }

        match (cigar_op.target_delta(), cigar_op.query_delta(strand)) {
            (0, query_delta) => {
                // Insertion in query (deletions in target)
                if target_pos >= requested_target_range.0 {
                    if !found_overlap {
                        projected_query_start = query_pos;
                        projected_target_start = target_pos;
                        first_op_idx = curr_op_idx;
                        found_overlap = true;
                    }
                    projected_query_end = query_pos + query_delta;
                    projected_target_end = target_pos;
                    last_op_idx = curr_op_idx + 1;
                }
                query_pos += query_delta;
            }
            (target_delta, 0) => {
                // Deletion in query (insertions in target)
                let overlap_start = target_pos.max(requested_target_range.0);
                let overlap_end = (target_pos + target_delta).min(last_target_pos);

                if overlap_start < overlap_end {
                    if !found_overlap {
                        projected_query_start = query_pos;
                        projected_target_start = overlap_start;
                        first_op_idx = curr_op_idx;
                        first_op_offset = overlap_start - target_pos;
                        found_overlap = true;
                    }
                    projected_query_end = query_pos;
                    projected_target_end = overlap_end;
                    last_op_idx = curr_op_idx + 1;
                    last_op_remaining = overlap_end - (target_pos + target_delta);
                }
                target_pos += target_delta;
            }
            (target_delta, query_delta) => {
                // Match or mismatch
                let overlap_start = target_pos.max(requested_target_range.0);
                let overlap_end = (target_pos + target_delta).min(requested_target_range.1);

                if overlap_start < overlap_end {
                    let overlap_length = overlap_end - overlap_start;
                    let query_overlap_start = query_pos + (overlap_start - target_pos) * dir;
                    let query_overlap_end = query_overlap_start + overlap_length * dir;

                    if !found_overlap {
                        projected_query_start = query_overlap_start;
                        projected_target_start = overlap_start;
                        first_op_idx = curr_op_idx;
                        first_op_offset = overlap_start - target_pos;
                        found_overlap = true;
                    }
                    projected_query_end = query_overlap_end;
                    projected_target_end = overlap_end;
                    last_op_idx = curr_op_idx + 1;
                    last_op_remaining = overlap_end - (target_pos + target_delta);
                }

                target_pos += target_delta;
                query_pos += query_delta;
            }
        }
    }

    // If we had at least one overlap, the variables were set
    // projected_query_start == projected_query_end in deletions in the query
    // projected_target_start == projected_target_end in insertions in the query
    if found_overlap
        && projected_query_start != projected_query_end
        && projected_target_start != projected_target_end
    {
        let mut projected_cigar_ops = cigar_ops[first_op_idx..last_op_idx].to_vec();
        // Adjust first operation length
        if first_op_offset > 0 {
            projected_cigar_ops[0].adjust_len(-first_op_offset);
        }
        // Adjust last operation length
        if last_op_remaining < 0 {
            projected_cigar_ops[last_op_idx - first_op_idx - 1].adjust_len(last_op_remaining);
        }

        Some((
            projected_query_start,
            projected_query_end,
            projected_cigar_ops,
            projected_target_start,
            projected_target_end,
        ))
    } else {
        None
    }
}

fn parse_cigar_to_delta(cigar: &str) -> Result<Vec<CigarOp>, ParseErr> {
    let mut ops = Vec::new();
    let mut len: i32 = 0;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            len = len * 10 + (c as i32 - '0' as i32);
        } else {
            let op = CigarOp::new(len, c);
            ops.push(op);
            len = 0;
        }
    }

    Ok(ops)
}

fn parse_variable_tracepoints(tp_str: &str) -> Result<Vec<(usize, Option<usize>)>, ParseErr> {
    Ok(tp_str
        .split(';')
        .filter_map(|s| {
            if s.contains(',') {
                // This has both coordinates
                let parts: Vec<&str> = s.split(',').collect();
                Some((parts[0].parse().unwrap(), Some(parts[1].parse().unwrap())))
            } else {
                // This has only first coordinate
                match s.parse() {
                    Ok(a) => Some((a, None)),
                    Err(_) => None,
                }
            }
        })
        .collect())
}

fn calculate_gap_compressed_identity(cigar_ops: &[CigarOp]) -> f64 {
    let (matches, mismatches, insertions, deletions) =
        cigar_ops
            .iter()
            .fold((0i32, 0i32, 0i32, 0i32), |(m, mm, i, d), op| {
                let len = op.len();
                match op.op() {
                    'M' | '=' => (m + len, mm, i, d), // Assume 'M' represents matches
                    'X' => (m, mm + len, i, d),
                    'I' => (m, mm, i + 1, d),
                    'D' => (m, mm, i, d + 1),
                    _ => (m, mm, i, d),
                }
            });

    let total = matches + mismatches + insertions + deletions;
    if total == 0 {
        0.0
    } else {
        (matches as f64) / (total as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paf::parse_paf;
    use std::io::BufReader;

    #[test]
    fn test_project_target_range_through_alignment_forward() {
        let target_range = (100, 200);
        let record = (100, 200, 0, 100, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let result =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();

        assert_eq!(result, (0, 100, cigar_ops.clone(), 100, 200));
    }

    #[test]
    fn test_project_target_range_through_alignment_reverse() {
        let target_range = (100, 200);
        let record = (100, 200, 0, 100, Strand::Reverse);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let result =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();

        assert_eq!(result, (100, 0, cigar_ops.clone(), 100, 200));
    }

    #[test]
    fn test_project_target_range_through_alignment() {
        let cigar_ops = vec![
            CigarOp::new(10, '='), // 10, 60
            CigarOp::new(5, 'I'),  // 10, 65
            CigarOp::new(5, 'D'),  // 15, 65
            CigarOp::new(50, '='), // 65, 115
            CigarOp::new(50, 'I'), // 65, 165
            CigarOp::new(35, '='), // 100, 200
        ];
        let base = (0, 100, 50, 200, Strand::Forward);
        {
            let result =
                project_target_range_through_alignment((0, 100), base, &cigar_ops).unwrap();
            assert_eq!(result, (50, 200, cigar_ops.clone(), 0, 100));
        }
        {
            let result =
                project_target_range_through_alignment((50, 55), base, &cigar_ops).unwrap();
            assert_eq!(result, (100, 105, vec![CigarOp::new(5, '=')], 50, 55));
        }
        {
            let result =
                project_target_range_through_alignment((50, 64), base, &cigar_ops).unwrap();
            assert_eq!(result, (100, 114, vec![CigarOp::new(14, '=')], 50, 64));
        }
        // We no longer output empty target ranges
        // {
        //     let result = project_target_range_through_alignment((65, 65), base, &cigar_ops).unwrap();
        //     assert_eq!(result, (115, 165, vec![CigarOp::new(50, 'I')], 65, 65));
        // }
        {
            let result =
                project_target_range_through_alignment((50, 65), base, &cigar_ops).unwrap();
            let cigar_ops = vec![CigarOp::new(15, '='), CigarOp::new(50, 'I')];
            assert_eq!(result, (100, 165, cigar_ops, 50, 65));
        }
        {
            let result =
                project_target_range_through_alignment((50, 66), base, &cigar_ops).unwrap();
            let cigar_ops = vec![
                CigarOp::new(15, '='),
                CigarOp::new(50, 'I'),
                CigarOp::new(1, '='),
            ];
            assert_eq!(result, (100, 166, cigar_ops, 50, 66));
        }
        {
            let result =
                project_target_range_through_alignment((70, 95), base, &cigar_ops).unwrap();
            assert_eq!(result, (170, 195, vec![CigarOp::new(25, '=')], 70, 95));
        }
    }

    // 1. Simple Forward Projection
    #[test]
    fn test_forward_projection_simple() {
        let target_range = (100, 200);
        let record = (100, 200, 100, 200, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let (query_start, query_end, cigar, target_start, target_end) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!(
            (query_start, query_end, cigar, target_start, target_end),
            (100, 200, vec![CigarOp::new(100, '=')], 100, 200)
        );
    }

    // 2. Simple Reverse Projection
    #[test]
    fn test_reverse_projection_simple() {
        let target_range = (100, 200);
        let record = (100, 200, 100, 200, Strand::Reverse);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let (query_start, query_end, cigar, target_start, target_end) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!(
            (query_start, query_end, cigar, target_start, target_end),
            (200, 100, vec![CigarOp::new(100, '=')], 100, 200)
        ); // Adjust for reverse calculation
    }

    // 3. Forward Projection with Insertions
    #[test]
    fn test_forward_projection_with_insertions() {
        let target_range = (50, 150);
        let record = (50, 150, 50, 160, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(50, '='), // Match
            CigarOp::new(10, 'I'), // Insertion
            CigarOp::new(50, '='), // Match
        ];
        let (start, end, cigar, _, _) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!((start, end, cigar), (50, 160, cigar_ops));
    }

    // 4. Forward Projection with Deletions
    #[test]
    fn test_forward_projection_with_deletions() {
        let target_range = (50, 150);
        let record = (50, 150, 50, 140, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(50, '='), // Match
            CigarOp::new(10, 'D'), // Deletion
            CigarOp::new(40, '='), // Match
        ];
        let (start, end, cigar, _, _) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!((start, end, cigar), (50, 140, cigar_ops));
    }

    // 5. Reverse Projection with Mixed Operations
    #[test]
    fn test_reverse_projection_with_mixed_operations() {
        let target_range = (150, 250);
        let record = (100, 200, 200, 300, Strand::Reverse);
        let cigar_ops = vec![
            CigarOp::new(50, '='), // 150, 250
            CigarOp::new(10, 'D'), // 150, 260
            CigarOp::new(10, 'I'), // 150, 250
            CigarOp::new(40, '='), // 150, 250
        ];
        let (start, end, cigar, _, _) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        let cigar_ops = vec![
            CigarOp::new(10, 'D'), // 150, 260
            CigarOp::new(10, 'I'), // 150, 250
            CigarOp::new(40, '='), // 150, 250
        ];
        assert_eq!((start, end, cigar), (250, 200, cigar_ops));
    }

    // 6. Edge Case Projection
    #[test]
    fn test_edge_case_projection() {
        let target_range = (0, 10);
        let record = (0, 50, 0, 40, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(10, '='), // Match
            CigarOp::new(20, 'D'), // Deletion in target
            CigarOp::new(8, '='),  // Match
            CigarOp::new(1, 'X'),  // Match
            CigarOp::new(1, '='),  // Match
            CigarOp::new(10, 'I'), // Insertion in query
            CigarOp::new(10, '='), // Match
        ];
        let (query_start, query_end, cigar, target_start, target_end) =
            project_target_range_through_alignment(target_range, record, &cigar_ops).unwrap();
        assert_eq!(
            (query_start, query_end, cigar, target_start, target_end),
            (0, 10, vec![CigarOp::new(10, '=')], 0, 10)
        );
    }

    #[test]
    fn test_parse_cigar_to_delta_basic() {
        let cigar = "10=5I5D";
        let cigar_ops = vec![
            CigarOp::new(10, '='),
            CigarOp::new(5, 'I'),
            CigarOp::new(5, 'D'),
        ];
        let ops = parse_cigar_to_delta(cigar).unwrap();
        assert_eq!(ops, cigar_ops);
    }

    // #[test]
    // fn test_parse_cigar_to_delta_invalid() {
    //     let cigar = "10=5Q"; // Q is not a valid CIGAR operation
    //     assert!(parse_cigar_to_delta(cigar).is_err());
    // }

    #[test]
    fn test_parse_paf_valid() {
        let paf_data = b"seq1\t100\t10\t20\t+\tt1\t200\t30\t40\t10\t20\t255\tcg:Z:10M\n";
        let mut seq_index = SequenceIndex::new();
        let query_id = seq_index.get_or_insert_id("seq1", Some(100));
        let target_id = seq_index.get_or_insert_id("t1", Some(200));
        let reader = BufReader::new(&paf_data[..]);
        let expected_records = vec![
            PartialPafRecord {
                query_id,
                query_start: 10,
                query_end: 20,
                target_id,
                target_start: 30,
                target_end: 40,
                strand_and_data_offset: 45, // Forward strand
                data_bytes: 3,
            },
            // Add more test records as needed
        ];
        let records = parse_paf(reader, &mut seq_index).unwrap();
        assert_eq!(records, expected_records);
    }

    #[test]
    fn test_is_tracepoints_data() {
        // Test CIGAR strings (should return false)
        assert!(!QueryMetadata::is_tracepoints_data(b"10M"));
        assert!(!QueryMetadata::is_tracepoints_data(b"123M45I67D"));
        assert!(!QueryMetadata::is_tracepoints_data(b"5=3X2I"));
        assert!(!QueryMetadata::is_tracepoints_data(b"100="));
        assert!(!QueryMetadata::is_tracepoints_data(b"50I"));
        assert!(!QueryMetadata::is_tracepoints_data(b"25D"));
        assert!(!QueryMetadata::is_tracepoints_data(b"10X"));

        // Test tracepoints with comma (should return true)
        assert!(QueryMetadata::is_tracepoints_data(b"123,456"));
        assert!(QueryMetadata::is_tracepoints_data(b"123,456;789,012"));
        assert!(QueryMetadata::is_tracepoints_data(b"0,10;20,30;40,50"));

        // Test tracepoints with semicolon only (should return true)
        assert!(QueryMetadata::is_tracepoints_data(b"123;456"));
        assert!(QueryMetadata::is_tracepoints_data(b"123;456;789"));
        assert!(QueryMetadata::is_tracepoints_data(b"0;10;20;30"));

        // Test mixed tracepoints (with some having only first coordinate)
        assert!(QueryMetadata::is_tracepoints_data(b"123,456;789;012,345"));

        // Test single number tracepoints (should return true)
        assert!(QueryMetadata::is_tracepoints_data(b"123"));
        assert!(QueryMetadata::is_tracepoints_data(b"0"));
        assert!(QueryMetadata::is_tracepoints_data(b"999"));

        // Edge cases (should return false)
        assert!(!QueryMetadata::is_tracepoints_data(b"abc"));
        assert!(!QueryMetadata::is_tracepoints_data(b"M123"));
        assert!(!QueryMetadata::is_tracepoints_data(b"=50"));
    }
}
