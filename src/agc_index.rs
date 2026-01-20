use ragc_core::{Decompressor, DecompressorConfig};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fmt;
use std::io::{self};
use std::sync::{Arc, Mutex};

// Structure to manage AGC archives using ragc-core
pub struct AgcIndex {
    decompressors: Arc<Mutex<Vec<Decompressor>>>,
    pub agc_paths: Vec<String>,
    // Interned strings - each unique string stored once
    interned_strings: FxHashMap<String, Arc<str>>,
    // Maps use Arc<str> for zero-cost sharing (String key for efficient lookup)
    sample_contig_to_agc: FxHashMap<String, usize>,
    // Precomputed mapping from contig name to (sample, full_contig_name, agc_idx)
    contig_to_sample_info: FxHashMap<String, (Arc<str>, Arc<str>, usize)>,
}

impl fmt::Debug for AgcIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AgcIndex")
            .field("agc_paths", &self.agc_paths)
            .field(
                "num_decompressors",
                &self.decompressors.lock().unwrap().len(),
            )
            .finish_non_exhaustive()
    }
}

impl AgcIndex {
    fn new() -> Self {
        AgcIndex {
            decompressors: Arc::new(Mutex::new(Vec::new())),
            agc_paths: Vec::new(),
            interned_strings: FxHashMap::default(),
            sample_contig_to_agc: FxHashMap::default(),
            contig_to_sample_info: FxHashMap::default(),
        }
    }

    /// Intern a string, returning an Arc<str> that can be shared
    fn intern(&mut self, s: &str) -> Arc<str> {
        if let Some(arc) = self.interned_strings.get(s) {
            Arc::clone(arc)
        } else {
            let arc: Arc<str> = s.into();
            self.interned_strings
                .insert(s.to_string(), Arc::clone(&arc));
            arc
        }
    }

    fn extract_short_contig_name(full_name: &str) -> &str {
        full_name.split_whitespace().next().unwrap_or(full_name)
    }

    pub fn build_from_files(agc_files: &[String]) -> io::Result<Self> {
        let mut index = AgcIndex::new();

        // Parallel metadata extraction phase
        let metadata_results: Vec<_> = agc_files
            .par_iter()
            .enumerate()
            .map(|(agc_idx, agc_path)| -> io::Result<(usize, String, Decompressor, Vec<(String, Vec<String>)>)> {
                let config = DecompressorConfig { verbosity: 0 };
                let mut decompressor = Decompressor::open(agc_path, config).map_err(|e| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("Failed to open AGC file: {agc_path}: {e}"),
                    )
                })?;

                // Get all samples and their contigs
                let samples = decompressor.list_samples();
                let sample_contigs: Vec<_> = samples
                    .into_iter()
                    .map(|sample| {
                        let contigs = decompressor.list_contigs(&sample).unwrap_or_default();
                        (sample, contigs)
                    })
                    .collect();

                Ok((agc_idx, agc_path.clone(), decompressor, sample_contigs))
            })
            .collect::<io::Result<Vec<_>>>()?;

        // Sequential assembly phase to maintain order and avoid shared mutable state issues
        for (agc_idx, agc_path, decompressor, sample_contigs) in metadata_results {
            index.agc_paths.push(agc_path);

            for (sample, contigs) in sample_contigs {
                let sample_arc = index.intern(&sample);

                for contig in contigs {
                    let contig_arc = index.intern(&contig);

                    // Key: contig@sample
                    let key = format!("{contig}@{sample}");
                    index.sample_contig_to_agc.insert(key, agc_idx);

                    // Key: contig alone (if unique)
                    index
                        .sample_contig_to_agc
                        .entry(contig.clone())
                        .or_insert(agc_idx);

                    // Contig -> (sample, contig, agc_idx) - values use Arc<str>
                    index
                        .contig_to_sample_info
                        .entry(contig.clone())
                        .or_insert((Arc::clone(&sample_arc), Arc::clone(&contig_arc), agc_idx));

                    // Handle short contig name if different
                    let short_contig = Self::extract_short_contig_name(&contig);
                    if short_contig != contig {
                        let short_key = format!("{short_contig}@{sample}");
                        index
                            .sample_contig_to_agc
                            .entry(short_key)
                            .or_insert(agc_idx);
                        index
                            .sample_contig_to_agc
                            .entry(short_contig.to_string())
                            .or_insert(agc_idx);
                        index
                            .contig_to_sample_info
                            .entry(short_contig.to_string())
                            .or_insert((Arc::clone(&sample_arc), Arc::clone(&contig_arc), agc_idx));
                    }
                }
            }

            index.decompressors.lock().unwrap().push(decompressor);
        }

        // Shrink to fit after building
        index.sample_contig_to_agc.shrink_to_fit();
        index.contig_to_sample_info.shrink_to_fit();

        // Clear the interner - we no longer need it after building
        index.interned_strings = FxHashMap::default();

        Ok(index)
    }

    fn parse_query(&self, seq_name: &str) -> (String, String, Option<usize>) {
        // Parse queries in the format:
        // - "contig@sample" -> (sample, contig, agc_idx)
        // - "contig" -> (sample, contig, agc_idx) if contig is unique

        if let Some((contig, sample)) = seq_name.split_once('@') {
            let agc_idx = self.sample_contig_to_agc.get(seq_name).copied();
            (sample.to_string(), contig.to_string(), agc_idx)
        } else if let Some((sample, full_contig, agc_idx)) =
            self.contig_to_sample_info.get(seq_name)
        {
            (sample.to_string(), full_contig.to_string(), Some(*agc_idx))
        } else {
            (String::new(), seq_name.to_string(), None)
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

        // Use efficient sub-sequence extraction (only decompresses needed segments)
        let mut decompressors = self.decompressors.lock().unwrap();
        let sequence = decompressors[agc_idx]
            .get_contig_range(&sample, &contig, start as usize, end as usize)
            .map_err(|e| {
                io::Error::other(format!(
                    "Failed to fetch sequence '{contig}@{sample}:{start}-{end}': {e}"
                ))
            })?;

        // Convert from numeric encoding (0-3) to ASCII (A,C,G,T)
        Ok(sequence
            .iter()
            .map(|&b| match b {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            })
            .collect())
    }

    pub fn get_sequence_length(&self, seq_name: &str) -> io::Result<usize> {
        let (sample, contig, agc_idx) = self.parse_query(seq_name);

        let agc_idx = agc_idx.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{seq_name}' not found in any AGC file"),
            )
        })?;

        let mut decompressors = self.decompressors.lock().unwrap();
        let length = decompressors[agc_idx]
            .get_contig_length(&sample, &contig)
            .map_err(|e| {
                io::Error::other(format!("Failed to get length for '{contig}@{sample}': {e}"))
            })?;

        Ok(length)
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
