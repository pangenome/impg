use std::collections::{HashMap, HashSet};
use coitrees::{COITree, Interval, IntervalTree};
use crate::paf::{PafRecord, ParseErr};
use crate::seqidx::SequenceIndex;

#[derive(Clone, Debug)]
#[derive(PartialEq, Eq, Hash)]
pub struct CigarOp {
    target_delta: i32,
    query_delta: i32,
}

#[derive(Clone, Debug)]
pub struct QueryMetadata {
    query_id: u32,
    cigar_ops: Vec<CigarOp>,
    target_start: i32,
    target_end: i32,
    query_start: i32,
    query_end: i32,
}

type QueryInterval = Interval<QueryMetadata>;
type TreeMap = HashMap<u32, COITree<QueryMetadata, u32>>;

pub struct Impg {
    trees: TreeMap,
}

impl Impg {
    /// Creates an Impg index from PAF records.
    ///
    /// # Examples
    ///
    /// ```
    /// let records = vec![paf::PafRecord { /* fields */ }];
    /// let seq_index = impg::SequenceIndex::new();
    /// // Add records to seq_index...
    /// let impg = Impg::from_paf_records(&records, &seq_index).unwrap();
    /// ```
    pub fn from_paf_records(records: &[PafRecord], seq_index: &SequenceIndex) -> Result<Self, ParseErr> {
        let mut intervals: HashMap<u32, Vec<Interval<QueryMetadata>>> = HashMap::new();

        for record in records {
            let cigar_ops = record.cigar.as_ref().map(|x| parse_cigar_to_delta(x)).transpose()?.unwrap_or_else(Vec::new);
            let query_id = seq_index.get_id(&record.query_name).expect("Query name not found in index");
            let target_id = seq_index.get_id(&record.target_name).expect("Target name not found in index");

            let query_metadata = QueryMetadata {
                query_id,
                cigar_ops,
                target_start: record.target_start as i32,
                target_end: record.target_end as i32,
                query_start: record.query_start as i32,
                query_end: record.query_end as i32,
            };

            intervals.entry(target_id).or_default().push(Interval {
                first: record.target_start as i32,
                last: record.target_end as i32,
                metadata: query_metadata,
            });
        }

        let trees: TreeMap = intervals.into_iter().map(|(target_id, interval_nodes)| {
            (target_id, COITree::new(interval_nodes.as_slice()))
        }).collect();

        Ok(Self { trees })
    }

    pub fn query(&self, target_id: u32, query_start: i32, query_end: i32) -> Vec<QueryInterval> {
        let mut results = Vec::new();
        if let Some(tree) = self.trees.get(&target_id) {
            tree.query(query_start, query_end, |interval| {
                let metadata = &interval.metadata;
                let (adjusted_start, adjusted_end) = project_target_range_through_alignment(
                    (interval.first, interval.last),
                    (metadata.target_start, metadata.target_end, metadata.query_start, metadata.query_end),
                    &metadata.cigar_ops
                );

                let adjusted_interval = QueryInterval {
                    first: adjusted_start,
                    last: adjusted_end,
                    metadata: metadata.clone(),
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
                    let metadata = &interval.metadata;
                    let (adjusted_start, adjusted_end) = project_target_range_through_alignment(
                        (interval.first, interval.last),
                        (metadata.target_start, metadata.target_end, metadata.query_start, metadata.query_end),
                        &metadata.cigar_ops
                    );

                    let adjusted_interval = QueryInterval {
                        first: adjusted_start,
                        last: adjusted_end,
                        metadata: metadata.clone(),
                    };
                    results.push(adjusted_interval);

                    if metadata.query_id != current_target {
                        stack.push((metadata.query_id, adjusted_end, adjusted_end));
                    }
                });
            }
        }

        results
    }
}

fn project_target_range_through_alignment(
    target_range: (i32, i32),
    record: (i32, i32, i32, i32),
    cigar_ops: &[CigarOp],
) -> (i32, i32) {
    let (target_start, target_end) = target_range;
    let (target_start, _, query_start, _) = record;

    let mut query_pos = query_start;
    let mut target_pos = target_start;

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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cigar_to_delta_basic() {
        let cigar = "10M5I5D";
        let expected_ops = vec![
            CigarOp { target_delta: 10, query_delta: 10 },
            CigarOp { target_delta: 0, query_delta: 5 },
            CigarOp { target_delta: 5, query_delta: 0 },
        ];
        let ops = parse_cigar_to_delta(cigar).unwrap();
        assert_eq!(ops, expected_ops);
    }

    #[test]
    fn test_parse_cigar_to_delta_invalid() {
        let cigar = "10M5X"; // 'X' is not a supported operation
        assert!(parse_cigar_to_delta(cigar).is_err());
    }

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
                cigar: Some("10M".to_string()),
            },
            // Add more test records as needed
        ];
        let records = parse_paf(reader).unwrap();
        assert_eq!(records, expected_records);
    }

    #[test]
    fn test_parse_paf_invalid() {
        let file = File::open("path/to/invalid/paf/file").unwrap();
        let reader = BufReader::new(file);
        let records = parse_paf(reader);
        assert!(records.is_err());
    }
}
