use coitrees::*;
use std::collections::HashMap;

// Assuming PafRecord and ParseErr are defined elsewhere in your project
use crate::paf::{PafRecord, ParseErr};

// Define QueryInterval to store the query name in its metadata field
type QueryInterval = Interval<String>;

// The TreeMap now stores COITrees with Interval metadata containing the QueryInterval
type TreeMap = HashMap<String, COITree<QueryInterval, u32>>;

pub struct Impg {
    trees: TreeMap,
}

impl Impg {
    pub fn from_paf_records(records: &[PafRecord]) -> Result<Self, ParseErr> {
        let mut intervals: HashMap<String, Vec<Interval<QueryInterval>>> = HashMap::new();

        // Collect intervals grouped by target name, each with a corresponding query range and name
        for record in records {
            let query_interval = Interval {
                first: record.query_start as i32,
                last: record.query_end as i32,
                metadata: record.query_name.clone(),  // Store the query name here
            };

            intervals.entry(record.target_name.clone())
                .or_default()
                .push(Interval {
                    first: record.target_start as i32,
                    last: record.target_end as i32,
                    metadata: query_interval,
                });
        }

        // Construct COITrees from the collected intervals
        let trees: TreeMap = intervals.into_iter()
            .map(|(seqname, interval_nodes)| {
                (seqname, COITree::new(&interval_nodes))
            })
            .collect();

        println!("built with {} trees", trees.len());

        Ok(Self { trees })
    }

    // Perform a direct query for overlaps on a specific target sequence
    // This returns the corresponding query intervals for the target range overlaps
    pub fn query(&self, target_name: &str, query_start: i32, query_end: i32) -> Vec<QueryInterval> {
        let mut results = Vec::new();
        println!("querying {}: {}-{}", target_name, query_start, query_end);
        if let Some(tree) = self.trees.get(target_name) {
            tree.query(query_start, query_end, |interval| {
                results.push(interval.metadata.clone());
            });
        }
        results
    }

    // Perform a transitive query for overlaps on a specific target sequence
    // This method finds all transitive overlaps for each target range, using the query name in the metadata
    pub fn query_transitive(&self, target_name: &str, query_start: i32, query_end: i32) -> Vec<QueryInterval> {
        let mut results = Vec::new();
        if let Some(tree) = self.trees.get(target_name) {
            let mut stack = vec![Interval { first: query_start, last: query_end, metadata: target_name.to_string() }];
            let mut visited = std::collections::HashSet::<(i32, i32, String)>::new();

            while let Some(interval) = stack.pop() {
                if let Some(query_tree) = self.trees.get(&interval.metadata) {
                    query_tree.query(interval.first, interval.last, |overlap| {
                        let other = overlap.metadata.clone();
                        let key = (other.first, other.last, other.metadata.clone());
                        if visited.insert(key) {
                            stack.push(overlap.metadata.clone());
                            results.push(overlap.metadata.clone());
                        }
                    });
                }
            }
        }
        results
    }
}

/// Parse a query string in the format "seqname:start-end" into a tuple
pub fn parse_query(query: &str) -> Result<(String, (i32, i32)), &'static str> {
    let parts: Vec<&str> = query.split(':').collect();
    if parts.len() != 2 {
        return Err("Invalid query format");
    }
    let seq_name = parts[0].to_string();
    let range_parts: Vec<&str> = parts[1].split('-').collect();
    if range_parts.len() != 2 {
        return Err("Invalid range format");
    }

    let start = range_parts[0].parse::<i32>().map_err(|_| "Invalid start range")?;
    let end = range_parts[1].parse::<i32>().map_err(|_| "Invalid end range")?;

    Ok((seq_name, (start, end)))
}
