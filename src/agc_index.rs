use ragc_core::{Decompressor, DecompressorConfig};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fmt;
use std::io::{self};
use std::sync::{Arc, Mutex};

// THE AGC'S OWN LAYOUT UNIT (the bind seek-floor anatomy's fix basis):
// the archive stores each contig as consecutive SEGMENTS (raw_length
// ~uniform); `get_contig_range` decompresses per overlapping segment
// (a pack read + a whole-pack zstd decompress + an LZ reference
// prepare + decode per segment, EVERY fetch — the per-fetch fixed
// cost the bind's ~6KB mean lane loads pay in full). The cumulative
// contribution table below reproduces `get_contig_range`'s own
// position mapping (contribution_i = raw_length_i for i = 0, else
// raw_length_i - kmer_length) from the archive's public descriptors,
// so a caller can fetch SEGMENT-ALIGNED crops — the same bytes for
// any requested sub-crop (a crop of fetch(a, b) IS fetch of the
// crop), with the per-segment fixed cost amortized over the whole
// segment.
pub struct AgcSegmentLayout {
    /// Cumulative segment start positions in contig coordinates
    /// (starts[0] = 0; starts[i + 1] = starts[i] + contribution_i).
    starts: Vec<u64>,
    /// Cumulative segment end positions (ends[i] = starts[i] +
    /// contribution_i; ends.last() = the contig's own length).
    ends: Vec<u64>,
    kmer_length: u32,
}

impl AgcSegmentLayout {
    /// The segment-aligned covering range for [start, end): the
    /// smallest union of whole segments covering the request. Every
    /// crop the caller serves out of a fetch of this range equals the
    /// direct fetch of that crop (the Sources fetch contract), so the
    /// alignment is a pure layout/prefetch choice — the archived
    /// segment is the decompression unit, and one aligned fetch
    /// decompresses each touched segment once instead of once per
    /// delta-fetch touching it.
    pub fn aligned_crop(&self, start: u64, end: u64) -> (u64, u64) {
        let contig_len = self.ends.last().copied().unwrap_or(0);
        if start >= end || self.starts.is_empty() || start >= contig_len {
            return (start, end);
        }
        // The segment containing `start` (starts is sorted; index of
        // the last start <= start).
        let i0 = self.starts.partition_point(|&s| s <= start) - 1;
        // The segment containing `end - 1`.
        let i1 = self
            .starts
            .partition_point(|&s| s < end)
            .saturating_sub(1)
            .min(self.ends.len() - 1);
        (self.starts[i0], self.ends[i1].max(end).min(*self.ends.last().unwrap()))
    }
}

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
    /// Per-contig cumulative segment boundary tables (the archive's own
    /// layout unit; built lazily from the public descriptors, cached).
    segment_layouts: Mutex<FxHashMap<String, Arc<AgcSegmentLayout>>>,
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
            segment_layouts: Mutex::new(FxHashMap::default()),
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

    /// The contig's cumulative segment layout (the archive's own
    /// decompression-unit boundaries, from the public descriptors —
    /// `get_contig_range`'s own position mapping, cached per contig).
    /// Returns None when the sequence is not in the archive.
    pub fn segment_layout(&self, seq_name: &str) -> Option<Arc<AgcSegmentLayout>> {
        if let Some(layout) = self.segment_layouts.lock().unwrap().get(seq_name) {
            return Some(Arc::clone(layout));
        }
        let (sample, contig, agc_idx) = self.parse_query(seq_name);
        let Some(agc_idx) = agc_idx else {
            return None;
        };
        let layout = match self.my_pool() {
            PoolRef::Rayon(slot) => {
                let mut guard = slot.lock().unwrap();
                let decomps = guard.get_or_insert_with(|| self.open_decompressors());
                Self::segment_layout_of(&mut decomps[agc_idx], &sample, &contig)
            }
            PoolRef::Other(pool) => {
                let mut guard = pool.lock().unwrap();
                if guard.is_empty() {
                    *guard = self.open_decompressors();
                }
                Self::segment_layout_of(&mut guard[agc_idx], &sample, &contig)
            }
        };
        let layout = Arc::new(layout?);
        self.segment_layouts
            .lock()
            .unwrap()
            .insert(seq_name.to_string(), Arc::clone(&layout));
        Some(layout)
    }
    fn segment_layout_of(
        decompressor: &mut Decompressor,
        sample: &str,
        contig: &str,
    ) -> Option<AgcSegmentLayout> {
        let segments = decompressor
            .get_contig_segments_desc(sample, contig)
            .ok()?;
        let kmer_length = decompressor.kmer_length;
        let mut starts = Vec::with_capacity(segments.len());
        let mut ends = Vec::with_capacity(segments.len());
        let mut position = 0u64;
        for (i, segment) in segments.iter().enumerate() {
            let raw_length = segment.raw_length as u64;
            let contribution = if i == 0 {
                raw_length
            } else {
                raw_length - kmer_length as u64
            };
            starts.push(position);
            position += contribution;
            ends.push(position);
        }
        Some(AgcSegmentLayout {
            starts,
            ends,
            kmer_length,
        })
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
    use super::*;
    use std::time::Instant;

    #[test]
    fn agc_segment_layout_aligned_crop_covers_the_request() {
        // The layout's own unit proof, on a synthetic cumulative table
        // (the same contribution formula get_contig_range uses): the
        // aligned crop must cover the request, stay inside the
        // contig, and widen only to whole segment boundaries.
        let layout = AgcSegmentLayout {
            // contributions 100, 90, 90 -> ends 100, 190, 280 (kmer 10)
            starts: vec![0, 100, 190],
            ends: vec![100, 190, 280],
            kmer_length: 10,
        };
        // Interior request inside one segment: aligns to that segment.
        assert_eq!(layout.aligned_crop(10, 20), (0, 100));
        // Request crossing a boundary: covers both whole segments.
        assert_eq!(layout.aligned_crop(95, 105), (0, 190));
        // Request spanning three segments.
        assert_eq!(layout.aligned_crop(0, 280), (0, 280));
        // start == a segment start: that segment's own start.
        assert_eq!(layout.aligned_crop(100, 150), (100, 190));
        // end == a segment end: that segment's own end.
        assert_eq!(layout.aligned_crop(150, 190), (100, 190));
        // Degenerate/empty ranges pass through unchanged (callers
        // never fetch them).
        assert_eq!(layout.aligned_crop(20, 20), (20, 20));
        assert_eq!(layout.aligned_crop(300, 400), (300, 400));
    }

    /// THE BIND SEEK-FLOOR ANATOMY RECEIPT (SS10.6.13 duty 2, the
    /// measured half — #[ignore], explicit-run at the campaign's AGC;
    /// skips with a note when the substrate is absent): the per-fetch
    /// cost of `fetch_sequence` at the bind's own request geometry.
    /// The bind's post-delta-extend lane loads average ~6KB (L1's
    /// receipt: 18,626-107,333 loads over ~111.4MB); the receipt's
    /// question: is the per-fetch cost ~constant (a fixed-cost floor —
    /// the per-segment pack read + whole-pack zstd decompress + LZ
    /// reference prepare + decode `get_contig_range` pays per touched
    /// segment, EVERY fetch) or ~proportional to bytes? The classes:
    /// (a) 900B random crops (the delta-extend's mean extend),
    /// (b) 6KB random crops (the lane-load mean), (c) sequential 64KB
    /// crops (the segment-aligned fix's own fetch class), (d) one
    /// 250KB crop (the territory-extent class). Deterministic LCG
    /// positions, single-threaded, warm page cache.
    #[test]
    #[ignore = "the campaign's AGC substrate (IMPG_BIND_ANATOMY_AGC, default the HPRC panel AGC) — explicit-run receipt"]
    fn bind_seek_floor_anatomy_receipt() {
        let agc_path = std::env::var("IMPG_BIND_ANATOMY_AGC").unwrap_or_else(|_| {
            "/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc".to_string()
        });
        let lane =
            std::env::var("IMPG_BIND_ANATOMY_LANE").unwrap_or_else(|_| "CHM13#0#chr1".into());
        if !std::path::Path::new(&agc_path).exists() {
            eprintln!("[bind-anatomy] AGC {agc_path} absent — skipping (explicit-run receipt)");
            return;
        }
        let index = AgcIndex::build_from_files(&[agc_path.clone()]).unwrap();
        let length = index.get_sequence_length(&lane).unwrap() as u64;
        let layout = index.segment_layout(&lane).expect("lane layout");
        eprintln!(
            "[bind-anatomy] lane {lane}: length {length}, {} segments, first raw {} (contribution {}), kmer {}",
            layout.starts.len(),
            layout.ends[0],
            layout.ends[0],
            layout.kmer_length,
        );
        assert!(layout.ends.last().copied().unwrap() <= length);
        // A deterministic LCG for reproducible positions.
        let mut lcg = 0x2545F4914F6CDD1Du64;
        let mut rand_below = |bound: u64| -> u64 {
            lcg = lcg.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (lcg >> 33) % bound.max(1)
        };
        let fetch = |start: u64, end: u64| index.fetch_sequence(&lane, start as i32, end as i32).unwrap();
        let verify = |seq: &[u8], start: u64, end: u64| {
            assert_eq!(seq.len(), (end - start) as usize, "fetch length at {start}");
            assert!(seq.iter().all(|&b| b"ACGTN".contains(&b)), "alphabet at {start}");
        };
        let receipt = |class: &str, fetches: usize, bytes: u64, nanos: u128| {
            eprintln!(
                "[bind-anatomy] {class}: {} fetches, {:.1} MB, {:.1} ms/fetch, {:.2} us/KB, {:.1} MB/s",
                fetches,
                bytes as f64 / (1024.0 * 1024.0),
                nanos as f64 / fetches as f64 / 1e6,
                nanos as f64 / (bytes / 1024).max(1) as f64 / 1e3,
                bytes as f64 / (nanos as f64 / 1e9) / (1024.0 * 1024.0),
            );
        };
        // (a) 900B random crops — the delta-extend mean.
        {
            let mut bytes = 0u64;
            let started = Instant::now();
            for _ in 0..200 {
                let start = rand_below(length - 10_000);
                let end = start + 900;
                let seq = fetch(start, end);
                verify(&seq, start, end);
                bytes += seq.len() as u64;
            }
            receipt("random 900B crops ", 200, bytes, started.elapsed().as_nanos());
        }
        // (b) 6KB random crops — the lane-load mean.
        {
            let mut bytes = 0u64;
            let started = Instant::now();
            for _ in 0..200 {
                let start = rand_below(length - 20_000);
                let end = start + 6_000;
                let seq = fetch(start, end);
                verify(&seq, start, end);
                bytes += seq.len() as u64;
            }
            receipt("random 6KB crops  ", 200, bytes, started.elapsed().as_nanos());
        }
        // (c) sequential 64KB crops from one random position (the
        // segment-aligned fetch class: adjacent whole segments).
        {
            let mut bytes = 0u64;
            let started = Instant::now();
            let base = rand_below(length - 5_000_000);
            for i in 0..64 {
                let start = base + i * 65_536;
                let end = start + 65_536;
                let seq = fetch(start, end);
                verify(&seq, start, end);
                bytes += seq.len() as u64;
            }
            receipt("sequential 64KB   ", 64, bytes, started.elapsed().as_nanos());
        }
        // (d) one 250KB crop — the territory-extent class.
        {
            let started = Instant::now();
            let start = rand_below(length.saturating_sub(300_000));
            let end = start + 250_000;
            let seq = fetch(start, end);
            verify(&seq, start, end);
            receipt("one 250KB crop    ", 1, seq.len() as u64, started.elapsed().as_nanos());
        }
        // (e) the same 900B crop repeated: the fixed-cost floor is
        // per-fetch (pack re-read + re-decompress + re-prepare), not
        // cached by the decompressor — the receipt the fix amortizes.
        {
            let started = Instant::now();
            let start = rand_below(length - 10_000);
            for _ in 0..50 {
                let seq = fetch(start, start + 900);
                verify(&seq, start, start + 900);
            }
            receipt("same 900B x50     ", 50, 50 * 900, started.elapsed().as_nanos());
        }
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
