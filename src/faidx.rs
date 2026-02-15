use rust_htslib::faidx;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::io::{self};
use std::path::Path;
use std::sync::Mutex;

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

        let reader = faidx::Reader::from_path(fasta_path).unwrap_or_else(|e| {
            panic!("Failed to open FASTA '{}': {}", fasta_path.display(), e)
        });
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
                    match faidx::Reader::from_path(fasta_path) {
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
                        index
                            .path_key_to_fasta
                            .insert(seq_name.to_string(), fasta_idx);

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

        let fasta_path = Path::new(&self.fasta_paths[fasta_idx]);

        let thread_idx = rayon::current_thread_index()
            .unwrap_or(self.thread_readers.len() - 1);
        let mut slot = self.thread_readers[thread_idx].lock().unwrap();
        let cache = slot.get_or_insert_with(|| ReaderCache::new(self.readers_per_thread));
        let reader = cache.get_or_open(fasta_idx, fasta_path);

        // rust-htslib uses 0-based half-open coordinates internally
        // but fetch_seq expects 0-based inclusive end coordinate
        let mut seq_vec = reader
            .fetch_seq(seq_name, start as usize, (end - 1) as usize)
            .map_err(|e| {
                io::Error::other(format!("Failed to fetch sequence for {seq_name}: {e}"))
            })?
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
