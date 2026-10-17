use ragc_core::{Decompressor, DecompressorConfig};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fmt;
use std::io::{self};
use std::sync::{Arc, Mutex};

// Structure to manage AGC archives using ragc-core
pub struct AgcIndex {
    /// Per-RAYON-thread decompressor pools. Indexed by
    /// rayon::current_thread_index(). Each slot is lazily initialized
    /// on first access by opening the AGC files.
    rayon_decompressors: Vec<Mutex<Option<Vec<Decompressor>>>>,
    /// Per-OS-thread decompressor pools for NON-rayon threads (the
    /// realign-score parallel seams run on scoped std::threads, not
    /// the rayon pool). Measured disease, the bind-seam profile's own
    /// pole: the committed code routed EVERY non-rayon thread to the
    /// ONE fallback slot behind a single Mutex, serializing all
    /// seam-thread AGC fetches (the parallel bind would have paid the
    /// serial decompression wall unchanged). Every OS thread now
    /// claims its OWN pool, lazily, keyed by thread id — decompressors
    /// are independent readers of immutable files, so which pool
    /// serves a fetch never changes the bytes returned.
    other_decompressors:
        Mutex<FxHashMap<std::thread::ThreadId, Arc<Mutex<Vec<Decompressor>>>>>,
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
            .field("num_rayon_slots", &self.rayon_decompressors.len())
            .finish_non_exhaustive()
    }
}

/// The calling thread's decompressor pool: a rayon thread's own
/// lazily-initialized slot, or a non-rayon (seam) thread's own
/// per-thread pool.
enum PoolRef<'a> {
    Rayon(&'a Mutex<Option<Vec<Decompressor>>>),
    Other(Arc<Mutex<Vec<Decompressor>>>),
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

        // Initialize per-rayon-thread decompressor slots (slot 0 reuses
        // the metadata decompressors to avoid re-opening AGC files)
        let num_slots = rayon::current_num_threads();
        let rayon_decompressors: Vec<_> = (0..num_slots).map(|_| Mutex::new(None)).collect();

        rayon_decompressors[0]
            .lock()
            .unwrap()
            .replace(saved_decomps);

        // Shrink to fit after building
        sample_contig_to_agc.shrink_to_fit();
        contig_to_sample_info.shrink_to_fit();

        Ok(AgcIndex {
            rayon_decompressors,
            other_decompressors: Mutex::new(FxHashMap::default()),
            agc_paths: agc_files.to_vec(),
            sample_contig_to_agc,
            contig_to_sample_info,
        })
    }

    /// The calling thread's decompressor pool: rayon threads use
    /// their own indexed slot (lazily initialized, exactly the
    /// committed code); every OTHER OS thread claims its OWN pool,
    /// lazily, keyed by thread id — never the ONE shared fallback
    /// slot the committed code put every non-rayon thread behind
    /// (the bind-seam profile's named serialization pole).
    fn open_decompressors(&self) -> Vec<Decompressor> {
        self.agc_paths
            .iter()
            .map(|path| {
                Decompressor::open(path, DecompressorConfig { verbosity: 0 })
                    .unwrap_or_else(|e| panic!("Failed to open AGC '{path}' for thread: {e}"))
            })
            .collect()
    }
    fn my_pool(&self) -> PoolRef<'_> {
        if let Some(index) = rayon::current_thread_index() {
            return PoolRef::Rayon(&self.rayon_decompressors[index]);
        }
        let thread_id = std::thread::current().id();
        let mut map = self.other_decompressors.lock().unwrap();
        PoolRef::Other(Arc::clone(
            map.entry(thread_id)
                .or_insert_with(|| Arc::new(Mutex::new(Vec::new()))),
        ))
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

        // Use per-thread decompressors to avoid mutex contention: rayon
        // threads take their own slot; every non-rayon OS thread claims
        // its own pool by thread id (the bind seam's scoped std::threads
        // would otherwise ALL queue on one shared fallback slot).
        let sequence = match self.my_pool() {
            PoolRef::Rayon(slot) => {
                let mut guard = slot.lock().unwrap();
                let decomps = guard.get_or_insert_with(|| self.open_decompressors());
                decomps[agc_idx].get_contig_range(
                    &sample,
                    &contig,
                    start as usize,
                    end as usize,
                )
            }
            PoolRef::Other(pool) => {
                let mut guard = pool.lock().unwrap();
                if guard.is_empty() {
                    *guard = self.open_decompressors();
                }
                guard[agc_idx].get_contig_range(
                    &sample,
                    &contig,
                    start as usize,
                    end as usize,
                )
            }
        }
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

        let length = match self.my_pool() {
            PoolRef::Rayon(slot) => {
                let mut guard = slot.lock().unwrap();
                let decomps = guard.get_or_insert_with(|| self.open_decompressors());
                decomps[agc_idx].get_contig_length(&sample, &contig)
            }
            PoolRef::Other(pool) => {
                let mut guard = pool.lock().unwrap();
                if guard.is_empty() {
                    *guard = self.open_decompressors();
                }
                guard[agc_idx].get_contig_length(&sample, &contig)
            }
        }
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
