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
}
