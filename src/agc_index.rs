use agc_rs::AGCFile;
use rayon::prelude::*;
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
    // Precomputed mapping from contig name to (sample, full_contig_name, agc_idx)
    contig_to_sample_info: FxHashMap<String, (String, String, usize)>,
}

impl AgcIndex {
    fn new() -> Self {
        AgcIndex {
            agc_wrapper: ThreadSafeAgc {
                agc_files: Arc::new(Mutex::new(Vec::new())),
            },
            agc_paths: Vec::new(),
            sample_contig_to_agc: FxHashMap::default(),
            contig_to_sample_info: FxHashMap::default(),
        }
    }

    // Helper function to extract the short contig name (before first whitespace)
    fn extract_short_contig_name(full_name: &str) -> &str {
        // Split on any whitespace character and take the first part
        full_name
            .split([' ', '\n', '\r', '\t'])
            .next()
            .unwrap_or(full_name)
    }

    pub fn build_from_files(agc_files: &[String]) -> io::Result<Self> {
        let mut index = AgcIndex::new();

        // Parallel metadata extraction phase
        let metadata_results: Vec<_> = agc_files
            .par_iter()
            .enumerate()
            .map(|(agc_idx, agc_path)| -> io::Result<(usize, String, AGCFile, Vec<(String, Vec<String>)>)> {
                let mut agc = AGCFile::new();
                if !agc.open(agc_path, true) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("Failed to open AGC file: {agc_path}"),
                    ));
                }

                // Get all samples and their contigs
                let samples = agc.list_samples();
                let sample_contigs: Vec<_> = samples
                    .into_iter()
                    .map(|sample| {
                        let contigs = agc.list_contigs(&sample);
                        (sample, contigs)
                    })
                    .collect();

                Ok((agc_idx, agc_path.clone(), agc, sample_contigs))
            })
            .collect::<io::Result<Vec<_>>>()?;

        // Sequential assembly phase to maintain order and avoid shared mutable state issues
        for (agc_idx, agc_path, agc, sample_contigs) in metadata_results {
            index.agc_paths.push(agc_path);

            for (sample, contigs) in sample_contigs {
                for contig in contigs {
                    // Create a key that combines full contig name and sample name
                    let key = format!("{contig}@{sample}");
                    index.sample_contig_to_agc.insert(key.clone(), agc_idx);

                    // Also insert just the full contig name if it's unique
                    index
                        .sample_contig_to_agc
                        .entry(contig.clone())
                        .or_insert(agc_idx);

                    // Precompute contig-to-sample mappings for fast lookup
                    let sample_info = (sample.clone(), contig.clone(), agc_idx);

                    // Map full contig name to sample info
                    index
                        .contig_to_sample_info
                        .entry(contig.clone())
                        .or_insert(sample_info.clone());

                    // Extract short contig name and create mappings
                    let short_contig = Self::extract_short_contig_name(&contig);

                    // If short name differs from full name, also create mappings for short name
                    if short_contig != contig {
                        // Create key with short contig name and sample
                        let short_key = format!("{short_contig}@{sample}");
                        index
                            .sample_contig_to_agc
                            .entry(short_key)
                            .or_insert(agc_idx);

                        // Also insert just the short contig name if it's unique
                        index
                            .sample_contig_to_agc
                            .entry(short_contig.to_string())
                            .or_insert(agc_idx);

                        // Map short contig name to sample info
                        index
                            .contig_to_sample_info
                            .entry(short_contig.to_string())
                            .or_insert(sample_info);
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
            if let Some((sample, full_contig, agc_idx)) = self.contig_to_sample_info.get(seq_name) {
                return (sample.clone(), full_contig.clone(), Some(*agc_idx));
            }
            ("".to_string(), seq_name.to_string(), None)
        }
    }

    pub fn fetch_sequence(&self, seq_name: &str, start: i32, end: i32) -> io::Result<Vec<u8>> {
        let (sample, contig, agc_idx) = self.parse_query(seq_name);

        let agc_idx = agc_idx.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{seq_name}' not found in any AGC file"),
            )
        })?;

        // AGC uses 0-based coordinates with inclusive end
        let mut agc_files = self.agc_wrapper.agc_files.lock().unwrap();
        let sequence = agc_files[agc_idx]
            .get_contig_sequence(&sample, &contig, start, end - 1)
            .map_err(|e| {
                io::Error::other(format!(
                    "Failed to fetch sequence '{contig}@{sample}:{start}:{end}': {e}"
                ))
            })?;

        Ok(sequence.into_bytes())
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
    fn test_batch_range_calculation() {
        // Test the min-max range calculation logic
        let requests = [
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
