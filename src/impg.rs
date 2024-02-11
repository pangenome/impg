use coitrees::{Interval, COITree, IntervalTree};
use std::collections::{HashMap, HashSet};
use crate::paf::{PafRecord, ParseErr};

#[derive(Clone, Debug)]
struct CigarOp {
    target_delta: i32,
    query_delta: i32,
}

#[derive(Clone, Debug)]
struct QueryMetadata {
    query_id: u32,
    cigar_ops: Vec<CigarOp>,
    mapping_target_start: i32,
    mapping_target_end: i32,
    mapping_query_start: i32,
    mapping_query_end: i32,
}

type QueryInterval = Interval<QueryMetadata>;
type TreeMap = HashMap<u32, COITree<QueryInterval, u32>>;

struct SequenceIndex {
    name_to_id: HashMap<String, u32>,
    next_id: u32,
}

impl SequenceIndex {
    fn new() -> Self {
        SequenceIndex {
            name_to_id: HashMap::new(),
            next_id: 0,
        }
    }

    fn get_or_insert_id(&mut self, name: &str) -> u32 {
        *self.name_to_id.entry(name.to_owned()).or_insert_with(|| {
            let id = self.next_id;
            self.next_id += 1;
            id
        })
    }

    fn get_id(&self, name: &str) -> Option<u32> {
        self.name_to_id.get(name).copied()
    }
}

pub struct Impg {
    trees: TreeMap,
}

impl Impg {
    pub fn from_paf_records(records: &[PafRecord], seq_index: &SequenceIndex) -> Result<Self, ParseErr> {
        let mut intervals: HashMap<u32, Vec<Interval<QueryMetadata>>> = HashMap::new();

        for record in records {
            let cigar_ops = record.cigar.as_ref().map(|x: &std::string::String| parse_cigar_to_delta(x)).transpose()?.unwrap_or_else(Vec::new);
            let query_id = seq_index.get_id(&record.query_name).expect("Query name not found in index");
            let target_id = seq_index.get_id(&record.target_name).expect("Target name not found in index");

            let query_metadata = QueryMetadata {
                query_id,
                cigar_ops,
                mapping_target_start: record.target_start as i32,
                mapping_target_end: record.target_end as i32,
                mapping_query_start: record.query_start as i32,
                mapping_query_end: record.query_end as i32,
            };

            intervals.entry(target_id).or_default().push(Interval {
                first: record.target_start as i32,
                last: record.target_end as i32,
                metadata: query_metadata,
            });
        }

        let trees: TreeMap = intervals.into_iter().map(|(target_id, interval_nodes)| (target_id, COITree::new(interval_nodes))).collect();
        Ok(Self { trees })
    }

    pub fn query(&self, target_id: u32, query_start: i32, query_end: i32) -> Vec<QueryInterval> {
        let mut results = Vec::new();
        if let Some(tree) = self.trees.get(&target_id) {
            tree.query(query_start, query_end, |interval| {
                let (adjusted_start, adjusted_end) = project_target_range_through_alignment((interval.first, interval.last), (interval.metadata.mapping_target_start, interval.metadata.mapping_target_end, interval.metadata.mapping_query_start, interval.metadata.mapping_query_end), &interval.metadata.cigar_ops);
                let adjusted_interval = QueryInterval {
                    first: adjusted_start,
                    last: adjusted_end,
                    metadata: interval.metadata.clone(),
                };
                results.push(adjusted_interval);
            });
        }
        results
    }

    pub fn query_transitive(&self, target_id: u32, query_start: i32, query_end: i32) -> Vec<QueryInterval> {
        let mut results = Vec::new();
        let mut stack = vec![(target_id, query_start, query_end)];
        let mut visited = HashSet::new();

        while let Some((current_target, current_start, current_end)) = stack.pop() {
            if !visited.insert((current_target, current_start, current_end)) {
                continue; // Skip if this query has already been processed
            }

            if let Some(tree) = self.trees.get(&current_target) {
                tree.query(current_start, current_end, |interval| {
                    let (adjusted_start, adjusted_end) = project_target_range_through_alignment((interval.first, interval.last), (interval.metadata.mapping_target_start, interval.metadata.mapping_target_end, interval.metadata.mapping_query_start, interval.metadata.mapping_query_end), &interval.metadata.cigar_ops);

                    let adjusted_interval = QueryInterval {
                        first: adjusted_start,
                        last: adjusted_end,
                        metadata: interval.metadata.clone(),
                    };
                    results.push(adjusted_interval);

                    if interval.metadata.query_id != current_target {
                        stack.push((interval.metadata.query_id, adjusted_start, adjusted_end));
                    }
                });
            }
        }

        results
    }
}

fn project_target_range_through_alignment(
    target_range: (i32, i32),
    mapping_record: (i32, i32, i32, i32),
    cigar_ops: &[CigarOp],
) -> (i32, i32) {
    let (target_start, target_end) = target_range;
    let (mapping_target_start, _, mapping_query_start, _) = mapping_record;

    let mut query_pos = mapping_query_start;
    let mut target_pos = mapping_target_start;

    let mut projected_start = -1;
    let mut projected_end = -1;

    for CigarOp { target_delta, query_delta } in cigar_ops {
        let target_delta_i32 = *target_delta as i32;
        let query_delta_i32 = *query_delta as i32;

        if target_pos >= target_start && projected_start == -1 {
            projected_start = query_pos + (target_start - target_pos).max(0);
        }

        if target_pos + target_delta_i32 > target_end {
            projected_end = query_pos + (target_end - target_pos).min(target_delta_i32);
            break;
        }

        target_pos += target_delta_i32;
        query_pos += query_delta_i32;
    }

    (projected_start, projected_end)
}

fn parse_cigar_to_delta(cigar: &str) -> Result<Vec<CigarOp>, ParseErr> {
    let mut ops = Vec::new();
    let mut num_buf = String::new();

    for c in cigar.chars() {
        if c.is_digit(10) {
            // Accumulate digit characters to build the length of the operation
            num_buf.push(c);
        } else {
            // When encountering an operation character, process the accumulated number and the operation
            match c {
                'M' | '=' | 'X' => { // 'M' for match/mismatch, '=' for match, 'X' for mismatch
                    let len = num_buf.parse::<usize>().map_err(|_| ParseErr::InvalidCigarFormat)?;
                    num_buf.clear(); // Reset the buffer for the next operation

                    // For 'M', '=', and 'X', both target and query positions advance
                    ops.extend(split_into_ops(len as i32, 1, 1));
                }
                'I' => { // Insertion to the query
                    let len = num_buf.parse::<usize>().map_err(|_| ParseErr::InvalidCigarFormat)?;
                    num_buf.clear();

                    // For an insertion, only the query position advances
                    ops.extend(split_into_ops(len as i32, 0, 1));
                }
                'D' => { // Deletion from the target
                    let len = num_buf.parse::<usize>().map_err(|_| ParseErr::InvalidCigarFormat)?;
                    num_buf.clear();

                    // For a deletion, only the target position advances
                    ops.extend(split_into_ops(len as i32, 1, 0));
                }
                _ => return Err(ParseErr::UnsupportedCigarOperation), // Unsupported operation
            }
        }
    }

    if !num_buf.is_empty() {
        // Handle any trailing numbers without an associated operation
        return Err(ParseErr::InvalidCigarFormat);
    }

    Ok(ops)
}

/// Splits large operations into multiple CigarOps if necessary, due to i32 limits.
fn split_into_ops(len: i32, target_delta: i32, query_delta: i32) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut remaining = len;

    while remaining > 0 {
        let delta = remaining.min(2147483647); // Max value for i32
        ops.push(CigarOp {
            target_delta: if target_delta > 0 { delta } else { 0 },
            query_delta: if query_delta > 0 { delta } else { 0 },
        });
        remaining -= delta;
    }

    ops
}
