use agc_rs::AGCFile;
use rustc_hash::FxHashMap;
use std::io::{self};
use std::sync::{Arc, Mutex};

// Wrapper to make AGC operations thread-safe
#[derive(Clone, Debug)]
struct ThreadSafeAgc {
    agc_files: Arc<Mutex<Vec<AGCFile>>>,
}

// Structure to manage AGC archives
#[derive(Debug)]
pub struct AgcIndex {
    agc_wrapper: ThreadSafeAgc,
    pub agc_paths: Vec<String>,
    sample_contig_to_agc: FxHashMap<String, usize>,
}

impl AgcIndex {
    fn new() -> Self {
        AgcIndex {
            agc_wrapper: ThreadSafeAgc {
                agc_files: Arc::new(Mutex::new(Vec::new())),
            },
            agc_paths: Vec::new(),
            sample_contig_to_agc: FxHashMap::default(),
        }
    }

    // Helper function to extract the short contig name (before first whitespace)
    fn extract_short_contig_name(full_name: &str) -> &str {
        // Split on any whitespace character and take the first part
        full_name
            .split(|c: char| c == ' ' || c == '\n' || c == '\r' || c == '\t')
            .next()
            .unwrap_or(full_name)
    }

    pub fn build_from_files(agc_files: &[String]) -> io::Result<Self> {
        let mut index = AgcIndex::new();

        for (agc_idx, agc_path) in agc_files.iter().enumerate() {
            index.agc_paths.push(agc_path.clone());

            let mut agc = AGCFile::new();
            if !agc.open(agc_path, true) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Failed to open AGC file: {}", agc_path),
                ));
            }

            // Get all samples in this AGC file
            let samples = agc.list_samples();

            for sample in samples {
                // Get all contigs for this sample
                let contigs = agc.list_contigs(&sample);

                for contig in contigs {
                    // Create a key that combines full contig name and sample name
                    let key = format!("{}@{}", contig, sample);
                    index.sample_contig_to_agc.insert(key.clone(), agc_idx);

                    // Also insert just the full contig name if it's unique
                    index.sample_contig_to_agc.entry(contig.clone()).or_insert(agc_idx);
                    
                    // Extract short contig name and create mappings
                    let short_contig = Self::extract_short_contig_name(&contig);
                    
                    // If short name differs from full name, also create mappings for short name
                    if short_contig != contig {
                        // Create key with short contig name and sample
                        let short_key = format!("{}@{}", short_contig, sample);
                        index.sample_contig_to_agc.entry(short_key).or_insert(agc_idx);
                        
                        // Also insert just the short contig name if it's unique
                        index.sample_contig_to_agc.entry(short_contig.to_string()).or_insert(agc_idx);
                    }
                }
            }

            index.agc_wrapper.agc_files.lock().unwrap().push(agc);
        }

        Ok(index)
    }

    fn parse_query(&self, seq_name: &str) -> (String, String, Option<usize>) {
        // Parse queries in the format:
        // - "contig@sample" -> (sample, contig, agc_idx)
        // - "contig" -> (sample, contig, agc_idx) if contig is unique

        if let Some((contig, sample)) = seq_name.split_once('@') {
            // Format: contig@sample
            let key = seq_name;
            let agc_idx = self.sample_contig_to_agc.get(key).copied();
            (sample.to_string(), contig.to_string(), agc_idx)
        } else {
            // Format: just contig name
            if let Some(&agc_idx) = self.sample_contig_to_agc.get(seq_name) {
                // Find which sample contains this contig
                let agc_files = self.agc_wrapper.agc_files.lock().unwrap();
                if let Some(agc) = agc_files.get(agc_idx) {
                    let samples = agc.list_samples();
                    for sample in samples {
                        let contigs = agc.list_contigs(&sample);
                        // Check both full names and short names
                        for contig in &contigs {
                            if contig == seq_name || Self::extract_short_contig_name(contig) == seq_name {
                                return (sample, contig.clone(), Some(agc_idx));
                            }
                        }
                    }
                }
            }
            ("".to_string(), seq_name.to_string(), None)
        }
    }

    pub fn fetch_sequence(&self, seq_name: &str, start: i32, end: i32) -> io::Result<Vec<u8>> {
        let (sample, contig, agc_idx) = self.parse_query(seq_name);

        let agc_idx = agc_idx.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{}' not found in any AGC file", seq_name),
            )
        })?;

        // AGC uses 0-based coordinates with inclusive end
        let mut agc_files = self.agc_wrapper.agc_files.lock().unwrap();
        let sequence = agc_files[agc_idx]
            .get_contig_sequence(&sample, &contig, start, end - 1)
            .map_err(|e| {
                io::Error::other(
                    format!(
                        "Failed to fetch sequence '{}@{}:{}:{}': {}",
                        contig, sample, start, end, e
                    ),
                )
            })?;

        Ok(sequence.into_bytes())
    }

    pub fn fetch_full_sequence(&self, seq_name: &str) -> io::Result<Vec<u8>> {
        let (sample, contig, agc_idx) = self.parse_query(seq_name);

        let agc_idx = agc_idx.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{}' not found in any AGC file", seq_name),
            )
        })?;

        let mut agc_files = self.agc_wrapper.agc_files.lock().unwrap();
        let sequence = agc_files[agc_idx]
            .get_full_contig(&sample, &contig)
            .map_err(|e| {
                io::Error::other(
                    format!(
                        "Failed to fetch full sequence '{}@{}': {}",
                        contig, sample, e
                    ),
                )
            })?;

        Ok(sequence.into_bytes())
    }

    /// Fetch multiple sequence subranges in a single batch operation.
    /// 
    /// This method groups requests by sequence/contig and fetches the min-max range once
    /// per group, then extracts individual subranges from the cached sequence.
    /// This reduces decompression overhead compared to individual fetches.
    /// 
    /// # Arguments
    /// * `requests` - A slice of (sequence_name, start, end) tuples
    /// 
    /// # Returns
    /// A vector of sequences in the same order as the input requests
    /// 
    /// # Performance Benefits
    /// - Reduces decompression from N operations to K operations (where K is unique sequences)
    /// - Fetches min-max range once per sequence/contig group
    /// - Extracts subranges from cached decompressed sequence
    pub fn fetch_sequences_batch(&self, requests: &[(String, i32, i32)]) -> io::Result<Vec<Vec<u8>>> {
        // Group requests by contig to minimize AGC calls
        let mut grouped_requests: FxHashMap<String, Vec<(usize, i32, i32)>> = 
            FxHashMap::default();
        
        // Parse and group requests by sample@contig
        for (idx, (seq_name, start, end)) in requests.iter().enumerate() {
            let (sample, contig, agc_idx) = self.parse_query(seq_name);
            if agc_idx.is_some() {
                let key = format!("{}@{}", contig, sample);
                grouped_requests.entry(key).or_default().push((idx, *start, *end));
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("Sequence '{}' not found", seq_name),
                ));
            }
        }
        
        // Pre-allocate result vector
        let mut results = vec![Vec::new(); requests.len()];
        
        // Lock AGC files once
        let mut agc_files = self.agc_wrapper.agc_files.lock().unwrap();
        
        // Fetch each group
        for (key, group) in grouped_requests {
            // Find the overall range for this contig
            let min_start = group.iter().map(|(_, s, _)| *s).min().unwrap();
            let max_end = group.iter().map(|(_, _, e)| *e).max().unwrap();
            
            // Parse key back to get sample and contig
            let (contig, sample) = key.split_once('@').unwrap();
            
            // Find which AGC file contains this contig
            let agc_idx = self.sample_contig_to_agc.get(&key).copied().unwrap();
            
            // Fetch the entire range once (AGC uses 0-based inclusive end)
            let full_sequence = agc_files[agc_idx]
                .get_contig_sequence(sample, contig, min_start, max_end - 1)
                .map_err(|e| {
                    io::Error::other(format!(
                        "Failed to fetch sequence '{}@{}:{}:{}': {}",
                        contig, sample, min_start, max_end, e
                    ))
                })?;
            
            let full_sequence_bytes = full_sequence.into_bytes();
            
            // Extract subranges for each request
            for (idx, start, end) in group {
                let offset_start = (start - min_start) as usize;
                let offset_end = (end - min_start) as usize;
                results[idx] = full_sequence_bytes[offset_start..offset_end].to_vec();
            }
        }
        
        Ok(results)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_query() {
        let index = AgcIndex::new();

        // Test contig@sample format
        let (sample, contig, _) = index.parse_query("chr1@sample1");
        assert_eq!(sample, "sample1");
        assert_eq!(contig, "chr1");

        // Test contig-only format
        let (sample, contig, _) = index.parse_query("chr1");
        assert_eq!(sample, "");
        assert_eq!(contig, "chr1");
    }

    #[test]
    fn test_fetch_sequences_batch_empty() {
        let index = AgcIndex::new();
        let requests = vec![];
        let result = index.fetch_sequences_batch(&requests);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 0);
    }

    #[test]
    fn test_batch_range_calculation() {
        // Test the min-max range calculation logic
        let requests = vec![
            ("chr1@sample1".to_string(), 100, 200),
            ("chr1@sample1".to_string(), 150, 300),
            ("chr1@sample1".to_string(), 50, 120),
        ];
        
        // Expected min_start = 50, max_end = 300
        // So we should fetch range 50:300, then extract:
        // - 100:200 -> relative 50:150
        // - 150:300 -> relative 100:250  
        // - 50:120 -> relative 0:70
        
        // This test validates the logic without requiring actual AGC files
        let min_start = requests.iter().map(|(_, start, _)| *start).min().unwrap();
        let max_end = requests.iter().map(|(_, _, end)| *end).max().unwrap();
        
        assert_eq!(min_start, 50);
        assert_eq!(max_end, 300);
        
        // Test relative coordinate calculation
        let (_, start1, end1) = &requests[0];
        let relative_start1 = (start1 - min_start) as usize;
        let relative_end1 = (end1 - min_start) as usize;
        assert_eq!(relative_start1, 50);
        assert_eq!(relative_end1, 150);
    }
}
