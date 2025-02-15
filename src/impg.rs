use rustc_hash::FxHashMap;
use coitrees::{BasicCOITree, Interval, IntervalTree};
use crate::paf::{PafRecord, ParseErr, Strand};
use crate::seqidx::SequenceIndex;
use serde::{Serialize, Deserialize};
use std::io::{Read, SeekFrom, Seek};
use std::fs::File;
use rayon::prelude::*;
use noodles::bgzf;
use regex::Regex;
use std::cmp::{min,max};
use log::debug;

/// Parse a CIGAR string into a vector of CigarOp
// Note that the query_delta is negative for reverse strand alignments
#[derive(Clone, Debug)]
#[derive(PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CigarOp {
    pub val: u32,
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
        Self { val: (val << 29) | (len as u32) }
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
            '=' | 'X' | 'I' | 'M' => if strand == Strand::Forward { self.len() } else { -self.len() },
            'D' => 0,
            _ => panic!("Invalid CIGAR operation: {}", self.op()),
        }
    }

    fn adjust_len(&mut self, length_delta: i32) {
        self.val = (self.val & (7 << 29)) | ((self.len() + length_delta) as u32);
    }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct QueryMetadata {
    query_id: u32,
    target_start: i32,
    target_end: i32,
    query_start: i32,
    query_end: i32,
    strand: Strand,
    cigar_offset: u64,
    cigar_bytes: usize,
}

impl QueryMetadata {
    fn get_cigar_ops(&self, paf_file: &String, paf_gzi_index: Option<&bgzf::gzi::Index>) -> Vec<CigarOp> {
        // Allocate space for cigar
        let mut cigar_buffer = vec![0; self.cigar_bytes];

        // Get reader and seek start of cigar str
        if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
            let mut reader = bgzf::Reader::new(File::open(paf_file).unwrap());
            reader.seek_by_uncompressed_position(paf_gzi_index.unwrap(), self.cigar_offset).unwrap();
            reader.read_exact(&mut cigar_buffer).unwrap();
        } else {
            let mut reader = File::open(paf_file).unwrap();
            reader.seek(SeekFrom::Start(self.cigar_offset)).unwrap();
            reader.read_exact(&mut cigar_buffer).unwrap();
        };

        let cigar_str: &str = std::str::from_utf8(&cigar_buffer).unwrap();
        parse_cigar_to_delta(cigar_str).ok().unwrap_or_default()
    }
}

pub type AdjustedInterval = (Interval<u32>, Vec<CigarOp>, Interval<u32>);
type TreeMap = FxHashMap<u32, BasicCOITree<QueryMetadata, u32>>;
pub type SerializableImpg = (FxHashMap<u32, Vec<SerializableInterval>>, SequenceIndex);

#[derive(Clone, Serialize, Deserialize)]
pub struct SerializableInterval {
    first: i32,
    last: i32,
    metadata: QueryMetadata,
}

#[derive(Debug, Default, Clone)]
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
        if i > 0 && (start - self.ranges[i-1].1).abs() < self.min_distance {
            start = self.ranges[i-1].1;
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

#[derive(Clone)]
pub struct Impg {
    pub trees: TreeMap,
    pub seq_index: SequenceIndex,
    pub paf_file: String,
    pub paf_gzi_index: Option<bgzf::gzi::Index>,
}

impl Impg {
    pub fn from_paf_records(records: &[PafRecord], paf_file: &str) -> Result<Self, ParseErr> {
        let paf_gzi_index: Option<bgzf::gzi::Index> = if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
            let paf_gzi_file = paf_file.to_owned() + ".gzi";
            Some(bgzf::gzi::read(paf_gzi_file.clone()).unwrap_or_else(|_| panic!("Could not open {}", paf_gzi_file)))
        } else {
            None
        };

        let mut seq_index = SequenceIndex::new();
        for record in records {
            seq_index.get_or_insert_id(&record.query_name, Some(record.query_length));
            seq_index.get_or_insert_id(&record.target_name, Some(record.target_length));
        }
        
        let intervals: FxHashMap<u32, Vec<Interval<QueryMetadata>>> = records.par_iter()
            .filter_map(|record| {
                let query_id = seq_index.get_id(&record.query_name).expect("Query name not found in index");
                let target_id = seq_index.get_id(&record.target_name).expect("Target name not found in index");

                let query_metadata = QueryMetadata {
                    query_id,
                    target_start: record.target_start as i32,
                    target_end: record.target_end as i32,
                    query_start: record.query_start as i32,
                    query_end: record.query_end as i32,
                    strand: record.strand,
                    cigar_offset: record.cigar_offset,
                    cigar_bytes: record.cigar_bytes
                };

                Some((target_id, Interval {
                    first: record.target_start as i32,
                    last: record.target_end as i32,
                    metadata: query_metadata,
                }))
            })  // Use fold and reduce to achieve grouping
            .fold(FxHashMap::default, |mut acc: FxHashMap<u32, Vec<Interval<QueryMetadata>>>, (target_id, interval)| {
                acc.entry(target_id).or_default().push(interval);
                acc
            })
            .reduce(FxHashMap::default, |mut acc, part| {
                for (key, value) in part {
                    acc.entry(key).or_default().extend(value);
                }
                acc
            });

        let trees: TreeMap = intervals.into_iter().map(|(target_id, interval_nodes)| {
            (target_id, BasicCOITree::new(interval_nodes.as_slice()))
        }).collect();

        Ok(Self { trees, seq_index, paf_file: paf_file.to_string(), paf_gzi_index })
    }

    pub fn to_serializable(&self) -> SerializableImpg {
        let serializable_trees = self.trees.iter().map(|(target_id, tree)| {
            let intervals = tree.iter().map(|interval| SerializableInterval {
                first: interval.first,
                last: interval.last,
                metadata: interval.metadata.clone(),
            }).collect();
            (*target_id, intervals)
        }).collect();
        (serializable_trees, self.seq_index.clone())
    }

    pub fn from_paf_and_serializable(paf_file: &str, serializable: SerializableImpg) -> Self {
        let (serializable_trees, seq_index) = serializable;
        let paf_gzi_index: Option<bgzf::gzi::Index> = if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
            let paf_gzi_file = paf_file.to_owned() + ".gzi";
            Some(bgzf::gzi::read(paf_gzi_file.clone()).unwrap_or_else(|_| panic!("Could not open {}", paf_gzi_file)))
        } else {
            None
        };
        let trees = serializable_trees.into_iter().map(|(target_id, intervals)| {
            let tree = BasicCOITree::new(intervals.iter().map(|interval| Interval {
                first: interval.first,
                last: interval.last,
                metadata: interval.metadata.clone(),
            }).collect::<Vec<_>>().as_slice());
            (target_id, tree)
        }).collect();
        Self { trees, seq_index, paf_file: paf_file.to_string(), paf_gzi_index }
    }

    pub fn query(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32
    ) -> Vec<AdjustedInterval> {
        let mut results = Vec::new();
        // Add the input range to the results
        results.push((
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
            vec![CigarOp::new(range_end - range_start, '=')],
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id
            }
        ));

        debug!("Querying region: {}:{}-{}, len: {}", self.seq_index.get_name(target_id).unwrap(), range_start, range_end, range_end - range_start);

        if let Some(tree) = self.trees.get(&target_id) {
            tree.query(range_start, range_end, |interval| {
                let metadata = &interval.metadata;
                let result = project_target_range_through_alignment(
                    (range_start, range_end),
                    (metadata.target_start, metadata.target_end, metadata.query_start, metadata.query_end, metadata.strand),
                    &metadata.get_cigar_ops(&self.paf_file, self.paf_gzi_index.as_ref()),
                    true
                );
                if let Some((adjusted_query_start, adjusted_query_end, adjusted_cigar, adjusted_target_start, adjusted_target_end)) = result {
                    let adjusted_interval = (
                        Interval {
                            first: adjusted_query_start,
                            last: adjusted_query_end,
                            metadata: metadata.query_id
                        },
                        adjusted_cigar,
                        Interval {
                            first: adjusted_target_start,
                            last: adjusted_target_end,
                            metadata: target_id
                        }
                    );
                    results.push(adjusted_interval);
                }
            });
        }
        
        debug!("Collected {} results", results.len());

        results
    }

    pub fn query_transitive(
        &self, 
        target_id: u32, 
        range_start: i32, 
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_region_size: i32,
        min_distance_between_ranges: i32,
        store_cigar: bool
    ) -> Vec<AdjustedInterval> {
        let mut results = Vec::new();
        // Add the input range to the results
        results.push((
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
            if store_cigar { vec![CigarOp::new(range_end - range_start, '=')] } else { Vec::new() },
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id
            }
        ));
        // Initialize stack with first query
        let mut stack = vec![(target_id, range_start, range_end, 0u16)];
        // Initialize visited ranges from masked regions if provided
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter()
            .map(|(&k, v)| (k, (*v).clone()))
            .collect()
        } else {
            (0..self.seq_index.len() as u32)
            .map(|id| {
                let len = self.seq_index.get_len_from_id(id).unwrap();
                (id, SortedRanges::new(len as i32, min_distance_between_ranges))
            })
            .collect()
        };
        // Initialize first visited range for target_id if not already present
        visited_ranges.entry(target_id)
            .or_default()
            .insert((range_start, range_end));

        while let Some((current_target_id, current_target_start, current_target_end, current_depth)) = stack.pop() {
            // Check if we've reached max depth
            if max_depth > 0 && current_depth >= max_depth {
                continue;
            }

            debug!("Querying region: {}:{}-{}, len: {}", self.seq_index.get_name(current_target_id).unwrap(), current_target_start, current_target_end, current_target_end - current_target_start);

            let prec_num_results = results.len();

            if let Some(tree) = self.trees.get(&current_target_id) {
                tree.query(current_target_start, current_target_end, |interval| {
                    let metadata = &interval.metadata;
                    let result = project_target_range_through_alignment(
                        (current_target_start, current_target_end),
                        (metadata.target_start, metadata.target_end, metadata.query_start, metadata.query_end, metadata.strand),
                        &metadata.get_cigar_ops(&self.paf_file, self.paf_gzi_index.as_ref()),
                        store_cigar
                    );
                    if let Some((adjusted_query_start, adjusted_query_end, adjusted_cigar, adjusted_target_start, adjusted_target_end)) = result {
                        let adjusted_interval = (
                            Interval {
                                first: adjusted_query_start,
                                last: adjusted_query_end,
                                metadata: metadata.query_id
                            },
                            adjusted_cigar,
                            Interval {
                                first: adjusted_target_start,
                                last: adjusted_target_end,
                                metadata: current_target_id
                            }
                        );
                        results.push(adjusted_interval);

                        // Only add non-overlapping portions to the stack for further exploration
                        if metadata.query_id != current_target_id {
                            let ranges = visited_ranges.entry(metadata.query_id)
                                .or_default();

                            let mut should_add = true;
                
                            // Check if the range is too close to any existing ranges
                            if min_distance_between_ranges > 0 {
                                let (new_min, new_max) = if adjusted_query_start <= adjusted_query_end {
                                    (adjusted_query_start, adjusted_query_end)
                                } else {
                                    (adjusted_query_end, adjusted_query_start)
                                };

                                // Find insertion point in sorted ranges
                                let idx = match ranges.ranges.binary_search_by_key(&new_min, |&(start, _)| start) {
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
                                let new_ranges = ranges.insert((adjusted_query_start, adjusted_query_end));

                                // Add non-overlapping portions to stack
                                for (new_start, new_end) in new_ranges {
                                    if (new_end - new_start).abs() >= min_transitive_region_size {
                                        stack.push((metadata.query_id, new_start, new_end, current_depth + 1));
                                    }
                                }
                            }
                        }
                    }
                });

                // Merge contiguous/overlapping ranges with same sequence_id
                let stack_size = stack.len();
                stack.sort_by_key(|(id, start, _, _)| (*id, *start));
                let mut write = 0;
                for read in 1..stack.len() {
                    if stack[write].0 == stack[read].0 &&   // Same sequence_id 
                        stack[write].2 >= stack[read].1 {   // Overlapping or contiguous
                        stack[write].2 = stack[write].2.max(stack[read].2);
                    } else {
                        write += 1;
                        stack.swap(write, read);
                    }
                }
                stack.truncate(write + 1);
                debug!("Merged stack size from {} to {}", stack_size, stack.len());
            }
            
            debug!("Collected {} results", results.len() - prec_num_results);
        }

        results
    }

}

fn project_target_range_through_alignment(
    requested_target_range: (i32, i32),
    record: (i32, i32, i32, i32, Strand),
    cigar_ops: &[CigarOp],
    store_cigar: bool
) -> Option<(i32, i32, Vec<CigarOp>, i32, i32)> {
    let (target_start, target_end, query_start, query_end, strand) = record;

    let dir = if strand == Strand::Forward { 1 } else { -1 };
    let mut query_pos = if strand == Strand::Forward { query_start } else { query_end };
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
            (0, query_delta) => { // Insertion in query (deletions in target)
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
            },
            (target_delta, 0) => { // Deletion in query (insertions in target)
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
            },
            (target_delta, query_delta) => { // Match or mismatch
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
            },
        }
    }

    // If we had at least one overlap, the variables were set
    // projected_query_start == projected_query_end in deletions in the query
    // projected_target_start == projected_target_end in insertions in the query
    (found_overlap && projected_query_start != projected_query_end && projected_target_start != projected_target_end).then(|| {
        let projected_cigar_ops = if store_cigar {
            let mut projected_cigar_ops = cigar_ops[first_op_idx..last_op_idx].to_vec();
            
            // Adjust first operation length
            if first_op_offset > 0 {
                projected_cigar_ops[0].adjust_len(-first_op_offset);
            }

            // Adjust last operation length
            if last_op_remaining < 0 {
                projected_cigar_ops[last_op_idx - first_op_idx - 1].adjust_len(last_op_remaining);
            }

            projected_cigar_ops
        } else {
            Vec::new()
        };

        (
            projected_query_start,
            projected_query_end,
            projected_cigar_ops,
            projected_target_start,
            projected_target_end,
        )
    })
}

fn parse_cigar_to_delta(cigar: &str) -> Result<Vec<CigarOp>, ParseErr> {
    let mut ops = Vec::new();
    let mut len: i32 = 0;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            len = len*10 + (c as i32 - '0' as i32);
        } else {
            let op = CigarOp::new(len, c);
            ops.push(op);
            len = 0;
        }
    }

    Ok(ops)
}

fn is_valid_cigar(cigar: &[CigarOp]) -> Result<(), String> {
    let cigar_str: String = cigar.iter().map(|op| format!("{}{}", op.len(), op.op())).collect();

    let re = Regex::new(r"^(\d+[MX=ID])+$").unwrap();
    if !re.is_match(&cigar_str) {
        return Err("Invalid format: non-standard or not-yet-supported operations, or formatting errors detected.".to_string());
    }

    let mut last_type = None;
    for op in cigar {
        let op_type = op.op();
        if let Some(last) = last_type {
            if "ISHP".contains(last) && "ISHP".contains(op_type) {
                return Err(format!("Consecutive non-reference-consuming operations detected: {} followed by {}.", last, op_type));
            }
        }
        last_type = Some(op_type);
    }

    Ok(())
}

fn parse_cigar(cigar: &[CigarOp]) -> (i32, i32) {
    let (query_length, target_length) = cigar.iter().fold((0, 0), |(query_len, target_len), op| {
        let len = op.len();
        match op.op() {
            'M' | 'X' | '=' | 'E' => (query_len + len, target_len + len),
            'I' | 'S' => (query_len + len, target_len),
            'D' | 'N' => (query_len, target_len + len),
            _ => (query_len, target_len),
        }
    });
    (query_length, target_length)
}

pub fn check_intervals(impg: &Impg, results: &Vec<AdjustedInterval>) -> Vec<(String, String)> {
    let mut invalid = Vec::new();

    for (overlap_query, cigar, overlap_target) in results {
        let query_name = impg.seq_index.get_name(overlap_query.metadata).unwrap();
        let query_len = impg.seq_index.get_len_from_id(overlap_query.metadata).unwrap();
        let target_name = impg.seq_index.get_name(overlap_target.metadata).unwrap();
        let target_len = impg.seq_index.get_len_from_id(overlap_target.metadata).unwrap();

        let (query_start, query_end) = (overlap_query.first, overlap_query.last);
        let (target_start, target_end) = (overlap_target.first, overlap_target.last);

        let full_cigar: String = cigar.iter().map(|op| format!("{}{}", op.len(), op.op())).collect();
        let first_chunk_cigar = if full_cigar.len() > 20 {
            format!("{}...", &full_cigar[..20])
        } else {
            full_cigar
        };

        let (calc_query_len, calc_target_len) = parse_cigar(cigar);

        let mut error_details = Vec::new();
        if calc_query_len != (query_end - query_start).abs() {
            error_details.push(format!("Query length mismatch: expected {} from the query range [{}-{}), got {} from the CIGAR string", (query_end - query_start).abs(), query_start, query_end, calc_query_len));
        }
        if calc_target_len != (target_end - target_start).abs() {
            error_details.push(format!("Target length mismatch: expected {} from the target range [{}-{}), got {} from the CIGAR string", (target_end - target_start).abs(), target_start, target_end, calc_target_len));
        }

        match is_valid_cigar(cigar) {
            Ok(()) => {
                if !error_details.is_empty() {
                    let error_reason = error_details.join("; ");
                    invalid.push((format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_name, query_len, query_start, query_end, if query_start <= query_end { '+' } else { '-' }, target_name, target_len, target_start, target_end, first_chunk_cigar), error_reason));
                }
            }
            Err(error_msg) => {
                let error_reason = if error_details.is_empty() {
                    error_msg
                } else {
                    format!("{}; {}", error_msg, error_details.join("; "))
                };
                invalid.push((format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", query_name, query_len, query_start, query_end, if query_start <= query_end { '+' } else { '-' }, target_name, target_len, target_start, target_end, first_chunk_cigar), error_reason));
            }
        }
    }

    invalid
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;
    use crate::paf::parse_paf;

    #[test]
    fn test_project_target_range_through_alignment_forward() {
        let target_range = (100, 200);
        let record = (100, 200, 0, 100, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let result = project_target_range_through_alignment(target_range, record, &cigar_ops, true).unwrap();

        assert_eq!(result, (0, 100, cigar_ops.clone(), 100, 200));
    }

    #[test]
    fn test_project_target_range_through_alignment_reverse() {
        let target_range = (100, 200);
        let record = (100, 200, 0, 100, Strand::Reverse);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let result = project_target_range_through_alignment(target_range, record, &cigar_ops, true).unwrap();

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
            let result = project_target_range_through_alignment((0, 100), base, &cigar_ops, true).unwrap();
            assert_eq!(result, (50, 200, cigar_ops.clone(), 0, 100));
        }
        {
            let result = project_target_range_through_alignment((50, 55), base, &cigar_ops, true).unwrap();
            assert_eq!(result, (100, 105, vec![CigarOp::new(5, '=')], 50, 55));
        }
        {
            let result = project_target_range_through_alignment((50, 64), base, &cigar_ops, true).unwrap();
            assert_eq!(result, (100, 114, vec![CigarOp::new(14, '=')], 50, 64));
        }
        // We no longer output empty target ranges
        // {
        //     let result = project_target_range_through_alignment((65, 65), base, &cigar_ops).unwrap();
        //     assert_eq!(result, (115, 165, vec![CigarOp::new(50, 'I')], 65, 65));
        // }
        {
            let result = project_target_range_through_alignment((50, 65), base, &cigar_ops, true).unwrap();
            let cigar_ops = vec![
                CigarOp::new(15, '='),
                CigarOp::new(50, 'I')
            ];
            assert_eq!(result, (100, 165, cigar_ops, 50, 65));
        }
        {
            let result = project_target_range_through_alignment((50, 66), base, &cigar_ops, true).unwrap();
            let cigar_ops = vec![
                CigarOp::new(15, '='),
                CigarOp::new(50, 'I'),
                CigarOp::new(1, '=')
            ];
            assert_eq!(result, (100, 166, cigar_ops, 50, 66));
        }
        {
            let result = project_target_range_through_alignment((70, 95), base, &cigar_ops, true).unwrap();
            assert_eq!(result, (170, 195, vec![CigarOp::new(25, '=')], 70, 95));
        }
    }

    // 1. Simple Forward Projection
    #[test]
    fn test_forward_projection_simple() {
        let target_range = (100, 200);
        let record = (100, 200, 100, 200, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let (query_start, query_end, cigar, target_start, target_end) = project_target_range_through_alignment(target_range, record, &cigar_ops, true).unwrap();
        assert_eq!((query_start, query_end, cigar, target_start, target_end), (100, 200, vec![CigarOp::new(100, '=')], 100, 200));
    }

    // 2. Simple Reverse Projection
    #[test]
    fn test_reverse_projection_simple() {
        let target_range = (100, 200);
        let record = (100, 200, 100, 200, Strand::Reverse);
        let cigar_ops = vec![CigarOp::new(100, '=')];
        let (query_start, query_end, cigar, target_start, target_end) = project_target_range_through_alignment(target_range, record, &cigar_ops, true).unwrap();
        assert_eq!((query_start, query_end, cigar, target_start, target_end), (200, 100, vec![CigarOp::new(100, '=')], 100, 200)); // Adjust for reverse calculation
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
        let (start, end, cigar, _, _) = project_target_range_through_alignment(target_range, record, &cigar_ops, true).unwrap();
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
        let (start, end, cigar, _, _) = project_target_range_through_alignment(target_range, record, &cigar_ops, true).unwrap();
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
        let (start, end, cigar, _, _) = project_target_range_through_alignment(target_range, record, &cigar_ops, true).unwrap();
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
        let (query_start, query_end, cigar, target_start, target_end) = project_target_range_through_alignment(target_range, record, &cigar_ops, true).unwrap();
        assert_eq!((query_start, query_end, cigar, target_start, target_end), (0, 10, vec![CigarOp::new(10, '=')], 0, 10));
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
        let reader = BufReader::new(&paf_data[..]);
        let expected_records = vec![
            PafRecord {
                query_name: "seq1".to_string(),
                query_length: 100,
                query_start: 10,
                query_end: 20,
                target_name: "t1".to_string(),
                target_length: 200,
                target_start: 30,
                target_end: 40,
                cigar_offset: 45,
                cigar_bytes: 3,
                strand: Strand::Forward,
            },
            // Add more test records as needed
        ];
        let records = parse_paf(reader).unwrap();
        assert_eq!(records, expected_records);
    }

}
