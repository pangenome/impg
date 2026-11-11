use rust_htslib::faidx;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::io::{self};
use std::path::Path;
use std::sync::Mutex;

/// Open a faidx reader. rust-htslib returns Ok with a NULL index when
/// htslib cannot load or build it (plain gzip, unreadable file), and the
/// first use then segfaults, so check the index first.
pub fn open_reader<P: AsRef<Path>>(path: P) -> io::Result<faidx::Reader> {
    let path = path.as_ref();
    let cpath =
        std::ffi::CString::new(path.as_os_str().as_encoded_bytes()).map_err(io::Error::other)?;
    let fai = unsafe { rust_htslib::htslib::fai_load(cpath.as_ptr()) };
    if fai.is_null() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "cannot load or build the FASTA index of '{}' (a .gz FASTA must be bgzip-compressed)",
                path.display()
            ),
        ));
    }
    unsafe { rust_htslib::htslib::fai_destroy(fai) };
    faidx::Reader::from_path(path).map_err(io::Error::other)
}

struct ReaderCache {
    entries: HashMap<usize, (faidx::Reader, u64)>,
    counter: u64,
    capacity: usize,
}

impl ReaderCache {
    fn new(capacity: usize) -> Self {
        Self {
            entries: HashMap::with_capacity(capacity),
            counter: 0,
            capacity,
        }
    }

    fn get_or_open(&mut self, fasta_idx: usize, fasta_path: &Path) -> &mut faidx::Reader {
        self.counter += 1;
        let counter = self.counter;

        if self.entries.contains_key(&fasta_idx) {
            let entry = self.entries.get_mut(&fasta_idx).unwrap();
            entry.1 = counter;
            return &mut entry.0;
        }

        // Evict least-recently-used reader if at capacity
        if self.entries.len() >= self.capacity {
            let lru_key = *self
                .entries
                .iter()
                .min_by_key(|(_, (_, ts))| *ts)
                .unwrap()
                .0;
            self.entries.remove(&lru_key);
        }

        let reader = open_reader(fasta_path)
            .unwrap_or_else(|e| panic!("Failed to open FASTA '{}': {}", fasta_path.display(), e));
        self.entries.insert(fasta_idx, (reader, counter));
        &mut self.entries.get_mut(&fasta_idx).unwrap().0
    }
}

fn get_soft_fd_limit() -> usize {
    std::fs::read_to_string("/proc/self/limits")
        .ok()
        .and_then(|contents| {
            contents
                .lines()
                .find(|l| l.starts_with("Max open files"))
                .and_then(|l| {
                    let token = l.split_whitespace().nth(3)?;
                    if token.eq_ignore_ascii_case("unlimited") {
                        Some(usize::MAX)
                    } else {
                        token.parse().ok()
                    }
                })
        })
        .unwrap_or(1024)
}

// Structure to manage multiple FASTA files
pub struct FastaIndex {
    pub fasta_paths: Vec<String>,
    pub path_key_to_fasta: FxHashMap<String, usize>,
    pub sequence_lengths: FxHashMap<String, usize>,
    thread_readers: Vec<Mutex<Option<ReaderCache>>>,
    readers_per_thread: usize,
}

impl std::fmt::Debug for FastaIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FastaIndex")
            .field("fasta_paths", &self.fasta_paths)
            .field("num_thread_slots", &self.thread_readers.len())
            .field("readers_per_thread", &self.readers_per_thread)
            .finish_non_exhaustive()
    }
}

impl FastaIndex {
    fn new() -> Self {
        FastaIndex {
            fasta_paths: Vec::new(),
            path_key_to_fasta: FxHashMap::default(),
            sequence_lengths: FxHashMap::default(),
            thread_readers: Vec::new(),
            readers_per_thread: 0,
        }
    }

    pub fn build_from_files(fasta_files: &[String]) -> io::Result<Self> {
        let mut index = FastaIndex::new();

        for (fasta_idx, fasta_path) in fasta_files.iter().enumerate() {
            index.fasta_paths.push(fasta_path.clone());

            // Read the .fai file to get sequence names
            let fai_path = format!("{fasta_path}.fai");

            // Try to open the .fai file, if it doesn't exist, try to create it
            let fai_content = match std::fs::read_to_string(&fai_path) {
                Ok(content) => content,
                Err(_) => {
                    // Try to create the index using rust-htslib
                    match open_reader(fasta_path) {
                        Ok(_) => {
                            // Index was created, now read it
                            std::fs::read_to_string(&fai_path)?
                        }
                        Err(e) => {
                            return Err(io::Error::other(format!(
                                "Failed to create FASTA index for '{fasta_path}': {e}"
                            )));
                        }
                    }
                }
            };

            // Parse the .fai file to get sequence names and lengths
            for line in fai_content.lines() {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() >= 2 {
                    let seq_name = fields[0];
                    if !seq_name.is_empty() {
                        if let Some(previous) = index
                            .path_key_to_fasta
                            .insert(seq_name.to_string(), fasta_idx)
                        {
                            if previous != fasta_idx {
                                // Later files silently won the lookup, so queries
                                // against this contig could read another sample's
                                // sequence.
                                log::warn!(
                                    "duplicate sequence name '{}' in '{}'; it shadows the record in '{}'",
                                    seq_name,
                                    fasta_path,
                                    index.fasta_paths[previous]
                                );
                            }
                        }

                        // Parse sequence length (second field in .fai file)
                        if let Ok(length) = fields[1].parse::<usize>() {
                            index.sequence_lengths.insert(seq_name.to_string(), length);
                        }
                    }
                }
            }
        }

        // Initialize per-thread reader cache slots
        let num_slots = rayon::current_num_threads() + 1;
        let fd_limit = get_soft_fd_limit();
        let reserved_fds = 64; // stdin/stdout/stderr, PAF, logs, etc.
        let available = fd_limit.saturating_sub(reserved_fds);
        index.readers_per_thread = (available / num_slots).max(1).min(fasta_files.len());
        index.thread_readers = (0..num_slots).map(|_| Mutex::new(None)).collect();

        Ok(index)
    }

    fn get_fasta_idx(&self, path_key: &str) -> Option<usize> {
        self.path_key_to_fasta.get(path_key).copied()
    }

    pub fn fetch_sequence(&self, seq_name: &str, start: i32, end: i32) -> io::Result<Vec<u8>> {
        let fasta_idx = self.get_fasta_idx(seq_name).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{seq_name}' not found in any FASTA file"),
            )
        })?;

        // Same rules as the AGC backend: a negative start or an inverted range
        // is a caller bug, an empty range is empty. `(end - 1) as usize` used to
        // underflow to `usize::MAX` when `end == 0`.
        if start < 0 || end < start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid range {start}-{end} for '{seq_name}' in the FASTA input"),
            ));
        }
        if end == start {
            return Ok(Vec::new());
        }

        let fasta_path = Path::new(&self.fasta_paths[fasta_idx]);

        // `thread_readers` is sized at build time from `rayon::current_num_threads()`.
        // Fetching from a later, larger pool gives a higher thread index, which used
        // to index out of bounds. Clamp into the last slot: sharing a slot only costs
        // some mutex contention, never correctness.
        let last_slot = self.thread_readers.len().saturating_sub(1);
        let thread_idx = rayon::current_thread_index()
            .map(|i| i.min(last_slot))
            .unwrap_or(last_slot);
        let mut slot = self.thread_readers[thread_idx].lock().unwrap();
        let cache = slot.get_or_insert_with(|| ReaderCache::new(self.readers_per_thread));
        let reader = cache.get_or_open(fasta_idx, fasta_path);

        // rust-htslib uses 0-based half-open coordinates internally
        // but fetch_seq expects 0-based inclusive end coordinate
        let mut seq_vec = reader
            .fetch_seq(seq_name, start as usize, (end - 1) as usize)
            .map_err(|e| io::Error::other(format!("Failed to fetch sequence for {seq_name}: {e}")))?
            .to_vec();
        seq_vec
            .iter_mut()
            .for_each(|byte| *byte = byte.to_ascii_uppercase());

        Ok(seq_vec)
    }

    pub fn get_sequence_length(&self, seq_name: &str) -> io::Result<usize> {
        self.sequence_lengths.get(seq_name).copied().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{seq_name}' not found"),
            )
        })
    }
}

#[cfg(test)]
mod open_reader_tests {
    use super::open_reader;
    use std::io::Write;

    #[test]
    fn plain_gzip_fasta_is_an_error_not_a_null_reader() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.fa.gz");
        let mut gz = flate2::write::GzEncoder::new(
            std::fs::File::create(&path).unwrap(),
            flate2::Compression::default(),
        );
        gz.write_all(b">s\nACGTACGT\n").unwrap();
        gz.finish().unwrap();
        let err = open_reader(&path).err().expect("plain gzip must fail");
        assert!(err.to_string().contains("bgzip"), "{err}");

        let plain = dir.path().join("x.fa");
        std::fs::write(&plain, b">s\nACGTACGT\n").unwrap();
        assert_eq!(open_reader(&plain).unwrap().n_seqs(), 1);
    }
}

#[cfg(test)]
mod fetch_tests {
    use super::FastaIndex;

    fn index() -> FastaIndex {
        let a = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/test_data/a.fa").to_string();
        FastaIndex::build_from_files(&[a]).unwrap()
    }

    #[test]
    fn ranges_follow_the_agc_rules() {
        let idx = index();
        assert_eq!(idx.fetch_sequence("chr1a", 0, 4).unwrap(), b"CTGA");
        assert!(idx.fetch_sequence("chr1a", 3, 3).unwrap().is_empty());
        // `(end - 1) as usize` used to wrap to usize::MAX here
        assert!(idx.fetch_sequence("chr1a", 0, 0).unwrap().is_empty());
        let err = idx.fetch_sequence("chr1a", 5, 2).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput, "{err}");
        assert!(idx.fetch_sequence("chr1a", -1, 2).is_err());
    }

    #[test]
    fn fetch_from_a_larger_pool_than_the_build_pool() {
        // Built in a 1-thread pool, fetched from 8 threads: thread index 7 used
        // to index past `thread_readers`
        let small = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let idx = small.install(index);
        let large = rayon::ThreadPoolBuilder::new().num_threads(8).build().unwrap();
        let seqs: Vec<Vec<u8>> = large.install(|| {
            use rayon::prelude::*;
            (0..64).into_par_iter().map(|_| idx.fetch_sequence("chr3a", 0, 9).unwrap()).collect()
        });
        assert!(seqs.iter().all(|s| s == b"AGTTTAGCT"));
    }
}
