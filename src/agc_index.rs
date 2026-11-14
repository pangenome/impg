use ragc_core::{Decompressor, DecompressorConfig};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fmt;
use std::io::{self};
use std::sync::atomic::{AtomicUsize, Ordering};
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

/// THE BOUNDED DECOMPRESSOR POOL'S DEFAULT: the K-sweep's operating
/// point. Each slot costs about an archive's decompression state, so
/// this number is a memory/time trade, measured - not a guess.
const DEFAULT_DECOMPRESSOR_SLOTS: usize = 8;

/// THE BOUNDED DECOMPRESSOR POOL'S SIZE: `IMPG_AGC_DECOMPRESSOR_SLOTS`
/// when set and parseable, else the default. The ambient thread count
/// is a CEILING, never a plan - see the field's own doc for why every
/// extra slot costs about an archive.
fn decompressor_slot_count() -> usize {
    let requested = std::env::var("IMPG_AGC_DECOMPRESSOR_SLOTS")
        .ok()
        .and_then(|value| value.trim().parse::<usize>().ok());
    bounded_slot_count(requested, rayon::current_num_threads())
}

/// The bound application, pure and testable: `None` (or a nonsense
/// zero) means the default; the ambient thread count is a CEILING on
/// the request, never a floor; the result is always >= 1.
fn bounded_slot_count(requested: Option<usize>, ceiling: usize) -> usize {
    requested
        .filter(|count| *count > 0)
        .unwrap_or(DEFAULT_DECOMPRESSOR_SLOTS)
        .min(ceiling.max(1))
        .max(1)
}

// Structure to manage AGC archives using ragc-core
pub struct AgcIndex {
    /// THE BOUNDED DECOMPRESSOR POOL (the memory-blowup fix): a fixed
    /// number of lazily-opened decompressor sets, shared round-robin by
    /// every caller - rayon threads and the realign-score seam's
    /// scoped std::threads alike.
    ///
    /// WHY BOUNDED: an AGC `Decompressor` holds the archive's own
    /// decompression state (the pack being expanded, the LZ reference
    /// state, the descriptor tables) in RAM, and that state measures
    /// about the archive's own size - 3.70 GB per slot against the
    /// 3.30 GB HPRC panel, the arithmetic of the routes derivation's
    /// measured 711 GB peak at 192 slots. The committed code opened ONE
    /// DECOMPRESSOR PER AMBIENT THREAD (`0..current_num_threads()`), so
    /// a 192-thread box paid 192 copies of the archive to avoid a mutex
    /// on an operation whose critical section IS the decompression.
    ///
    /// SHARING COSTS NO CORRECTNESS: decompressors are independent
    /// readers of immutable files, and whichever slot serves a fetch
    /// never changes the bytes returned (the previous per-thread design
    /// said exactly this in its own comment). Contention is bounded by
    /// the slot count, which is the tunable.
    decompressor_slots: Vec<Mutex<Option<Vec<Decompressor>>>>,
    /// Round-robin cursor over `decompressor_slots`.
    next_slot: AtomicUsize,
    /// How many slots have actually opened their decompressor sets
    /// (`<= decompressor_slots.len()`) - the memory bound's own
    /// regression receipt.
    slots_opened: AtomicUsize,
    pub agc_paths: Vec<String>,
    // Maps use Arc<str> for zero-cost sharing (String key for efficient lookup)
    sample_contig_to_agc: FxHashMap<String, usize>,
    // Precomputed mapping from contig name to (sample, full_contig_name, agc_idx)
    contig_to_sample_info: FxHashMap<String, (Arc<str>, Arc<str>, usize)>,
    /// Per-contig cumulative segment boundary tables (the archive's own
    /// layout unit; built lazily from the public descriptors, cached).
    segment_layouts: Mutex<FxHashMap<String, Arc<AgcSegmentLayout>>>,
    /// Every (sample, contig) in archive order.
    entries: Vec<(Arc<str>, Arc<str>)>,
}

impl fmt::Debug for AgcIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AgcIndex")
            .field("agc_paths", &self.agc_paths)
            .field("decompressor_slots", &self.decompressor_slots.len())
            .field("slots_opened", &self.slots_opened.load(Ordering::Relaxed))
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
        let mut entries: Vec<(Arc<str>, Arc<str>)> = Vec::new();
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
                    entries.push((Arc::clone(&sample_arc), Arc::clone(&contig_arc)));

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

        // THE BOUNDED SLOT COUNT (the memory-blowup fix's knob): the
        // ambient thread count is the CEILING, never the plan - each
        // slot costs about an archive's decompression state, so the
        // committed `0..current_num_threads()` opened 192 copies on a
        // 192-thread box (the measured 711 GB peak). Slot 0 reuses the
        // metadata decompressors to avoid re-opening the AGC files.
        let slot_count = decompressor_slot_count();
        let decompressor_slots: Vec<_> = (0..slot_count).map(|_| Mutex::new(None)).collect();

        decompressor_slots[0]
            .lock()
            .unwrap()
            .replace(saved_decomps);

        // Shrink to fit after building
        sample_contig_to_agc.shrink_to_fit();
        contig_to_sample_info.shrink_to_fit();

        Ok(AgcIndex {
            decompressor_slots,
            next_slot: AtomicUsize::new(0),
            slots_opened: AtomicUsize::new(1),
            agc_paths: agc_files.to_vec(),
            sample_contig_to_agc,
            contig_to_sample_info,
            segment_layouts: Mutex::new(FxHashMap::default()),
            entries,
        })
    }

    /// Open one decompressor set (one decompressor per AGC archive).
    /// Callers hold their slot's lock while using it, so a set is never
    /// aliased.
    fn open_decompressors(&self) -> Vec<Decompressor> {
        self.agc_paths
            .iter()
            .map(|path| {
                Decompressor::open(path, DecompressorConfig { verbosity: 0 })
                    .unwrap_or_else(|e| panic!("Failed to open AGC '{path}' for thread: {e}"))
            })
            .collect()
    }
    /// THE SLOT FOR THIS FETCH: round-robin over the bounded pool.
    /// Never per-thread state, so the process's decompressor footprint
    /// is bounded by the slot count no matter how many threads (rayon or
    /// scoped std::thread) ever touch the archive.
    fn decompressor_slot(&self) -> &Mutex<Option<Vec<Decompressor>>> {
        let slot = self.next_slot.fetch_add(1, Ordering::Relaxed) % self.decompressor_slots.len();
        &self.decompressor_slots[slot]
    }

    /// Borrow an initialised decompressor set from the pool for the
    /// duration of `f` (the slot's lock IS the borrow). A slot opens its
    /// decompressors on first use, under its own lock, so at most once.
    fn with_decompressors<T>(&self, f: impl FnOnce(&mut Vec<Decompressor>) -> T) -> T {
        let mut guard = self.decompressor_slot().lock().unwrap();
        if guard.is_none() {
            self.slots_opened.fetch_add(1, Ordering::Relaxed);
        }
        let decomps = guard.get_or_insert_with(|| self.open_decompressors());
        f(decomps)
    }

    /// The pool's size (the memory bound's own receipt field).
    pub fn decompressor_slot_count(&self) -> usize {
        self.decompressor_slots.len()
    }

    /// How many slots have actually opened their decompressor sets.
    pub fn decompressors_opened(&self) -> usize {
        self.slots_opened.load(Ordering::Relaxed)
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
        let layout = self.with_decompressors(|decomps| {
            Self::segment_layout_of(&mut decomps[agc_idx], &sample, &contig)
        });
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
        Some(AgcSegmentLayout { starts, ends })
    }

    /// One fetchable name per sequence, in archive order: the contig name
    /// when no other sample has it, else `contig@sample`.
    pub fn sequence_names(&self) -> Vec<String> {
        let mut contig_count: FxHashMap<&str, usize> = FxHashMap::default();
        for (_, contig) in &self.entries {
            *contig_count.entry(contig.as_ref()).or_default() += 1;
        }
        self.entries
            .iter()
            .map(|(sample, contig)| {
                if contig_count[contig.as_ref()] == 1 {
                    contig.to_string()
                } else {
                    format!("{contig}@{sample}")
                }
            })
            .collect()
    }

    pub fn fetch_sequence(&self, seq_name: &str, start: i32, end: i32) -> io::Result<Vec<u8>> {
        let requested = checked_range_len(seq_name, start, end)?;
        let (sample, contig, agc_idx) = self.parse_query(seq_name);

        let agc_idx = agc_idx.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{seq_name}' not found in any AGC file"),
            )
        })?;

        // THE BOUNDED POOL (see the field's doc): this fetch borrows one
        // of a fixed number of shared decompressor sets. No thread owns
        // a slot, so the archive's decompression state is held
        // `decompressor_slots` times at most - never once per thread.
        let sequence = self
            .with_decompressors(|decomps| {
                decomps[agc_idx].get_contig_range(&sample, &contig, start as usize, end as usize)
            })
            .map_err(|e| {
            io::Error::other(format!(
                "Failed to fetch sequence '{contig}@{sample}:{start}-{end}': {e}"
            ))
        })?;

        check_fetched_length(seq_name, requested, sequence.len(), start, end)?;

        Ok(numeric_to_ascii(&sequence))
    }

    pub fn get_sequence_length(&self, seq_name: &str) -> io::Result<usize> {
        let (sample, contig, agc_idx) = self.parse_query(seq_name);

        let agc_idx = agc_idx.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{seq_name}' not found in any AGC file"),
            )
        })?;

        let length = self
            .with_decompressors(|decomps| {
                decomps[agc_idx].get_contig_length(&sample, &contig)
            })
            .map_err(|e| {
            io::Error::other(format!("Failed to get length for '{contig}@{sample}': {e}"))
        })?;

        Ok(length)
    }
}

/// Validate an AGC range and return the number of bases it asks for.
fn checked_range_len(seq_name: &str, start: i32, end: i32) -> io::Result<usize> {
    crate::sequence_index::checked_range_len(seq_name, start, end, "the AGC archive")
}

/// AGC returns a short sequence when the requested range runs past the end of
/// the contig, and an archive appended with agc older than 3.2.2 returns wrong
/// ranges for the appended samples. Either way a short fetch shifts every
/// downstream coordinate, so it is an error, not a shorter answer.
fn check_fetched_length(
    seq_name: &str,
    requested: usize,
    fetched: usize,
    start: i32,
    end: i32,
) -> io::Result<()> {
    if fetched == requested {
        return Ok(());
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "AGC returned {fetched} bases for '{seq_name}:{start}-{end}', but {requested} were requested. \
             The range goes past the end of the contig, or the archive was appended with agc older than 3.2.2, \
             which returns wrong ranges for appended samples. Rebuild the archive with agc >= 3.2.2."
        ),
    ))
}

/// AGC numeric bases (0-3) to ASCII ACGT; any other code becomes N.
pub fn numeric_to_ascii(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fetched_range_must_have_the_requested_length() {
        // An exact fetch is fine.
        assert!(check_fetched_length("chr1@ref", 100, 100, 0, 100).is_ok());

        // A short fetch shifts every downstream coordinate, so it is an error
        // that names the sequence, the range and both lengths.
        let err = check_fetched_length("chr1@ref", 100, 40, 0, 100).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("chr1@ref:0-100"), "{msg}");
        assert!(msg.contains("40 bases"), "{msg}");
        assert!(msg.contains("100 were requested"), "{msg}");
        assert!(msg.contains("3.2.2"), "{msg}");

        // A long fetch is equally wrong.
        assert!(check_fetched_length("chr1@ref", 100, 101, 0, 100).is_err());
    }

    #[test]
    fn agc_range_length_rejects_reversed_and_negative_ranges() {
        assert_eq!(checked_range_len("chr1@ref", 10, 20).unwrap(), 10);
        assert_eq!(checked_range_len("chr1@ref", 7, 7).unwrap(), 0);

        for (start, end) in [(-1, 10), (20, 10), (0, -5)] {
            let err = checked_range_len("chr1@ref", start, end)
                .expect_err("invalid range must be rejected");
            assert!(
                err.to_string().contains("invalid range"),
                "unexpected error: {err}"
            );
        }
    }

    /// THE MEMORY BOUND IS THE POOL SIZE (the 711 GB regression's own
    /// guard): the ambient thread count may never inflate the pool, an
    /// explicit request may restore the old behaviour for a measured
    /// comparison, and a nonsense request degrades to the default.
    #[test]
    fn decompressor_pool_is_bounded_by_the_request_not_the_thread_count() {
        // the default at a 192-thread box: bounded, NOT 192
        assert_eq!(bounded_slot_count(None, 192), DEFAULT_DECOMPRESSOR_SLOTS);
        // an explicit request is honoured (the old behaviour, opt-in)
        assert_eq!(bounded_slot_count(Some(192), 192), 192);
        // ... but never above the ceiling
        assert_eq!(bounded_slot_count(Some(64), 16), 16);
        // zero/small ceiling/absent request all stay sane
        assert_eq!(bounded_slot_count(Some(0), 16), DEFAULT_DECOMPRESSOR_SLOTS);
        assert_eq!(bounded_slot_count(Some(4), 0), 1);
        assert_eq!(bounded_slot_count(None, 1), 1);
    }

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
}
