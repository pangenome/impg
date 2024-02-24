use std::collections::{HashMap, HashSet};
use coitrees::{COITree, Interval, IntervalTree};
use crate::paf::{PafRecord, ParseErr, Strand};
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
    strand: Strand,
}

type QueryInterval = Interval<QueryMetadata>;
type TreeMap = HashMap<u32, COITree<QueryMetadata, u32>>;

pub struct Impg {
    trees: TreeMap,
}

impl Impg {
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
                strand: record.strand,
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
                    (metadata.target_start, metadata.target_end, metadata.query_start, metadata.query_end, metadata.strand),
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
                        (metadata.target_start, metadata.target_end, metadata.query_start, metadata.query_end, metadata.strand),
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
    record: (i32, i32, i32, i32, Strand),
    cigar_ops: &[CigarOp],
) -> (i32, i32) {
    let (target_start, target_end, query_start, query_end, strand) = record;

    let mut query_pos = query_start;
    let mut target_pos = target_start;

    let mut projected_start = -1;
    let mut projected_end = -1;

    for CigarOp { target_delta, query_delta } in cigar_ops.iter() {
        if target_pos >= target_range.0 && projected_start == -1 {
            projected_start = query_pos;
        }

        if target_pos + target_delta > target_end {
            projected_end = query_pos + (target_end - target_pos).min(*query_delta);
            break;
        }

        target_pos += target_delta;
        query_pos += match strand {
            Strand::Forward => *query_delta,
            Strand::Reverse => -(*query_delta), // Reverse the query position increment for reverse strand
        };
    }

    if projected_end == -1 { // Ensure projected_end is set if the loop completes without setting it
        projected_end = query_pos;
    }

    assert!(projected_start != -1 && projected_end != -1, "Projection failed to calculate valid start and end positions");

    match strand {
        Strand::Forward => (projected_start, projected_end),
        Strand::Reverse => (query_start - projected_end, query_start - projected_start), // Adjust for reverse strand
    }
}

fn parse_cigar_to_delta(cigar: &str) -> Result<Vec<CigarOp>, ParseErr> {
    let mut ops = Vec::new();
    let mut num_buf = String::new();

    for c in cigar.chars() {
        if c.is_digit(10) {
            num_buf.push(c);
        } else {
            let len = num_buf.parse::<i32>().map_err(|_| ParseErr::InvalidCigarFormat)?;
            num_buf.clear(); // Reset the buffer for the next operation

            match c {
                'M' | '=' | 'X' => ops.push(CigarOp { target_delta: len, query_delta: len }),
                'I' => ops.push(CigarOp { target_delta: 0, query_delta: len }),
                'D' => ops.push(CigarOp { target_delta: len, query_delta: 0 }),
                _ => return Err(ParseErr::UnsupportedCigarOperation),
            }
        }
    }

    Ok(ops)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::BufReader;
    use crate::paf::parse_paf;

    #[test]
    fn test_project_target_range_through_alignment_forward() {
        let target_range = (100, 200);
        let record = (100, 200, 0, 100, Strand::Forward);
        let cigar_ops = vec![CigarOp { target_delta: 100, query_delta: 100 }];
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);

        assert_eq!(start, 0);
        assert_eq!(end, 100);
    }

    #[test]
    fn test_project_target_range_through_alignment_reverse() {
        let target_range = (100, 200);
        let record = (100, 200, 100, 200, Strand::Reverse);
        let cigar_ops = vec![CigarOp { target_delta: 100, query_delta: 100 }];
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);

        assert_eq!(start, 0);
        assert_eq!(end, 100);
    }

    #[test]
    fn test_project_target_range_through_alignment() {
        let cigar_ops = vec![
            CigarOp { target_delta: 10, query_delta: 10 },
            CigarOp { target_delta: 0, query_delta: 5 },
            CigarOp { target_delta: 5, query_delta: 0 },
        ];
        let (target_start, target_end, query_start, query_end, strand) = (0, 100, 10, 20, Strand::Forward);
        let expected_result = (10, 20);
        let result = project_target_range_through_alignment((0, 100), (target_start, target_end, query_start, query_end, strand), &cigar_ops);
        assert_eq!(result, expected_result);
    }

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
                strand: Strand::Forward,
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
