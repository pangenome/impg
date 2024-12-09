use std::collections::{HashMap, VecDeque};
use coitrees::{BasicCOITree, Interval, IntervalTree};
use crate::impg::{AdjustedInterval, Impg};
pub struct Progressive<'a> {
    impg: &'a Impg,
    region_queue: VecDeque<Interval<u32>>,
    genome_masks: HashMap<u32, BasicCOITree<(), u32>>,
}

impl<'a> Progressive<'a> {
    pub fn new(impg: &'a Impg, initial_regions: Vec<Interval<u32>>) -> Self {
        let mut region_queue = VecDeque::new();
        region_queue.extend(initial_regions);
        
        Progressive {
            impg,
            region_queue,
            genome_masks: HashMap::new(),
        }
    }

    pub fn process(&mut self) -> Vec<Vec<AdjustedInterval>> {
        let mut all_results = Vec::new();
        
        while let Some(region) = self.region_queue.pop_front() {
            eprintln!("process: region {:?}", region);

            // Process the region and collect results
            let partition_results = self.process_region(&region);
            
            // Only add non-empty partitions
            if !partition_results.is_empty() {
                all_results.push(partition_results);
            } else {
                eprintln!("No results for region {:?}", region);
            }
        }

        all_results
    }

    fn process_region(&mut self, region: &Interval<u32>) -> Vec<AdjustedInterval> {
        // Query the impg index for transitive overlaps with masking filter
        let results = self.impg.query_transitive(
            region.metadata,
            region.first, 
            region.last,
            Some(&|interval: &Interval<u32>| !self.is_masked(interval))
        );

        // Update masks and enqueue new regions based on results
        for (query_interval, _, _) in results.iter() {
            eprintln!("\tprocess_region: query_interval {:?}", query_interval);
            // Update the genome mask for this interval
            self.update_genome_mask(query_interval);

            // Identify and enqueue potential new regions to explore
            self.enqueue_new_regions(query_interval, region);
        }

        // Return all discovered intervals for this region
        results
    }

    fn is_masked(&self, interval: &Interval<u32>) -> bool {
        self.genome_masks
            .get(&interval.metadata)
            .map_or(false, |tree| tree.query_count(interval.first, interval.last) > 0)
    }

    fn update_genome_mask(&mut self, interval: &Interval<u32>) {

    }

    fn identify_extended_regions(&self, interval: &Interval<u32>, original_region: &Interval<u32>) -> Vec<Interval<u32>> {
        let mut extended_regions = Vec::new();
                
        extended_regions
    }

    fn enqueue_new_regions(&mut self, interval: &Interval<u32>, original_region: &Interval<u32>) {

    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paf::PafRecord;

    #[test]
    fn test_progressive_basic() {
        // Create test data
        let records = vec![
            PafRecord {
                query_name: "seq1".to_string(),
                query_length: 100,
                query_start: 0,
                query_end: 50,
                target_name: "seq2".to_string(),
                target_length: 100,
                target_start: 0,
                target_end: 50,
                strand: crate::paf::Strand::Forward,
                cigar_offset: 0,
                cigar_bytes: 0,
            }
        ];

        let impg = Impg::from_paf_records(&records, "test.paf").unwrap();
        
        let initial_regions = vec![
            Interval {
                first: 0,
                last: 50,
                metadata: 0, // seq1's id
            }
        ];

        let mut progressive = Progressive::new(&impg, initial_regions);
        let results = progressive.process();

        assert!(!results.is_empty());
    }
}