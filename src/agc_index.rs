use ragc_core::{Decompressor, DecompressorConfig};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fmt;
use std::io::{self};
use std::sync::{Arc, Mutex};

// Structure to manage AGC archives using ragc-core
pub struct AgcIndex {
    /// Per-thread decompressor pools. Indexed by rayon::current_thread_index().
    /// Each slot is lazily initialized on first access by opening the AGC files.
    thread_decompressors: Vec<Mutex<Option<Vec<Decompressor>>>>,
    pub agc_paths: Vec<String>,
    // Maps use Arc<str> for zero-cost sharing (String key for efficient lookup)
    sample_contig_to_agc: FxHashMap<String, usize>,
    // Precomputed mapping from contig name to (sample, full_contig_name, agc_idx)
    contig_to_sample_info: FxHashMap<String, (Arc<str>, Arc<str>, usize)>,
}

impl fmt::Debug for AgcIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AgcIndex")
            .field("agc_paths", &self.agc_paths)
            .field("num_thread_slots", &self.thread_decompressors.len())
            .finish_non_exhaustive()
    }
}

impl AgcIndex {
    fn extract_short_contig_name(full_name: &str) -> &str {
        full_name.split_whitespace().next().unwrap_or(full_name)
    }

    pub fn build_from_files(agc_files: &[String]) -> io::Result<Self> {
        // Parallel metadata extraction phase — load only contig names, skip segment details
        let metadata_results: Vec<_> = agc_files
            .par_iter()
            .enumerate()
            .map(|(agc_idx, agc_path)| -> io::Result<(usize, Decompressor, Vec<(String, Vec<String>)>)> {
                let config = DecompressorConfig { verbosity: 0 };
                let mut decompressor = Decompressor::open(agc_path, config).map_err(|e| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("Failed to open AGC file: {agc_path}: {e}"),
                    )
                })?;

                // Get all samples and their contig names (skip loading segment details)
                let samples = decompressor.list_samples();
                let sample_contigs: Vec<_> = samples
                    .into_iter()
                    .map(|sample| {
                        let contigs = decompressor.list_contigs_names_only(&sample).unwrap_or_default();
                        (sample, contigs)
                    })
                    .collect();

                Ok((agc_idx, decompressor, sample_contigs))
            })
            .collect::<io::Result<Vec<_>>>()?;

        // Sequential assembly phase
        let mut interned_strings: FxHashMap<String, Arc<str>> = FxHashMap::default();
        let mut sample_contig_to_agc: FxHashMap<String, usize> = FxHashMap::default();
        let mut contig_to_sample_info: FxHashMap<String, (Arc<str>, Arc<str>, usize)> =
            FxHashMap::default();

        let mut intern = |s: &str| -> Arc<str> {
            if let Some(arc) = interned_strings.get(s) {
                Arc::clone(arc)
            } else {
                let arc: Arc<str> = s.into();
                interned_strings.insert(s.to_string(), Arc::clone(&arc));
                arc
            }
        };

        // Keep the metadata decompressors to reuse in thread slot 0
        let mut saved_decomps: Vec<Decompressor> = Vec::with_capacity(agc_files.len());
        for (agc_idx, decompressor, sample_contigs) in metadata_results {
            assert_eq!(agc_idx, saved_decomps.len());
            saved_decomps.push(decompressor);

            for (sample, contigs) in sample_contigs {
                let sample_arc = intern(&sample);

                for contig in contigs {
                    let contig_arc = intern(&contig);

                    // Key: contig@sample
                    let key = format!("{contig}@{sample}");
                    sample_contig_to_agc.insert(key, agc_idx);

                    // Key: contig alone (if unique)
                    sample_contig_to_agc
                        .entry(contig.clone())
                        .or_insert(agc_idx);

                    // Contig -> (sample, contig, agc_idx) - values use Arc<str>
                    contig_to_sample_info.entry(contig.clone()).or_insert((
                        Arc::clone(&sample_arc),
                        Arc::clone(&contig_arc),
                        agc_idx,
                    ));

                    // Handle short contig name if different
                    let short_contig = Self::extract_short_contig_name(&contig);
                    if short_contig != contig {
                        let short_key = format!("{short_contig}@{sample}");
                        sample_contig_to_agc.entry(short_key).or_insert(agc_idx);
                        sample_contig_to_agc
                            .entry(short_contig.to_string())
                            .or_insert(agc_idx);
                        contig_to_sample_info
                            .entry(short_contig.to_string())
                            .or_insert((Arc::clone(&sample_arc), Arc::clone(&contig_arc), agc_idx));
                    }
                }
            }
        }

        // Initialize per-thread decompressor slots (+1 for main-thread fallback)
        let num_slots = rayon::current_num_threads() + 1;
        let thread_decompressors: Vec<_> = (0..num_slots).map(|_| Mutex::new(None)).collect();

        // Reuse the metadata decompressors in thread slot 0 to avoid re-opening AGC files
        thread_decompressors[0]
            .lock()
            .unwrap()
            .replace(saved_decomps);

        // Shrink to fit after building
        sample_contig_to_agc.shrink_to_fit();
        contig_to_sample_info.shrink_to_fit();

        Ok(AgcIndex {
            thread_decompressors,
            agc_paths: agc_files.to_vec(),
            sample_contig_to_agc,
            contig_to_sample_info,
        })
    }

    fn parse_query(&self, seq_name: &str) -> (String, String, Option<usize>) {
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

        // Use per-thread decompressors to avoid mutex contention
        let thread_idx =
            rayon::current_thread_index().unwrap_or(self.thread_decompressors.len() - 1);
        let mut slot = self.thread_decompressors[thread_idx].lock().unwrap();
        let decomps = slot.get_or_insert_with(|| {
            self.agc_paths
                .iter()
                .map(|path| {
                    Decompressor::open(path, DecompressorConfig { verbosity: 0 })
                        .unwrap_or_else(|e| panic!("Failed to open AGC '{path}' for thread: {e}"))
                })
                .collect()
        });
        let sequence = decomps[agc_idx]
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

        let thread_idx =
            rayon::current_thread_index().unwrap_or(self.thread_decompressors.len() - 1);
        let mut slot = self.thread_decompressors[thread_idx].lock().unwrap();
        let decomps = slot.get_or_insert_with(|| {
            self.agc_paths
                .iter()
                .map(|path| {
                    Decompressor::open(path, DecompressorConfig { verbosity: 0 })
                        .unwrap_or_else(|e| panic!("Failed to open AGC '{path}' for thread: {e}"))
                })
                .collect()
        });
        let length = decomps[agc_idx]
            .get_contig_length(&sample, &contig)
            .map_err(|e| {
                io::Error::other(format!("Failed to get length for '{contig}@{sample}': {e}"))
            })?;

        Ok(length)
    }
}

#[cfg(test)]
mod tests {
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
