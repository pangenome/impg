use rust_htslib::faidx;
use rustc_hash::FxHashMap;
use std::io::{self};

// Structure to manage multiple FASTA files
#[derive(Debug)]
pub struct FastaIndex {
    pub fasta_paths: Vec<String>,
    pub path_key_to_fasta: FxHashMap<String, usize>,
}

impl FastaIndex {
    fn new() -> Self {
        FastaIndex {
            fasta_paths: Vec::new(),
            path_key_to_fasta: FxHashMap::default(),
        }
    }

    pub fn build_from_files(fasta_files: &[String]) -> io::Result<Self> {
        let mut index = FastaIndex::new();

        for (fasta_idx, fasta_path) in fasta_files.iter().enumerate() {
            index.fasta_paths.push(fasta_path.clone());

            // Read the .fai file to get sequence names
            let fai_path = format!("{}.fai", fasta_path);

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
                                "Failed to create FASTA index for '{}': {}",
                                fasta_path, e
                            )));
                        }
                    }
                }
            };

            // Parse the .fai file to get sequence names
            for line in fai_content.lines() {
                if let Some(seq_name) = line.split('\t').next() {
                    if !seq_name.is_empty() {
                        index
                            .path_key_to_fasta
                            .insert(seq_name.to_string(), fasta_idx);
                    }
                }
            }
        }

        Ok(index)
    }

    fn get_fasta_path(&self, path_key: &str) -> Option<&str> {
        self.path_key_to_fasta
            .get(path_key)
            .map(|&idx| self.fasta_paths[idx].as_str())
    }

    pub fn fetch_sequence(&self, seq_name: &str, start: i32, end: i32) -> io::Result<Vec<u8>> {
        let fasta_path = self.get_fasta_path(seq_name).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{}' not found in any FASTA file", seq_name),
            )
        })?;

        let reader = faidx::Reader::from_path(fasta_path).map_err(|e| {
            io::Error::other(format!("Failed to open FASTA file '{}': {}", fasta_path, e))
        })?;

        // rust-htslib uses 0-based half-open coordinates internally
        // but fetch_seq expects 0-based inclusive end coordinate
        let sequence = reader
            .fetch_seq(seq_name, start as usize, (end - 1) as usize)
            .map_err(|e| {
                io::Error::other(format!(
                    "Failed to fetch sequence '{}:{}:{}': {}",
                    seq_name, start, end, e
                ))
            })?;

        // Apply the fix for the rust_htslib memory leak bug
        // https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171
        let seq_vec = sequence.to_vec();
        // Free the memory allocated by htslib to prevent memory leak
        unsafe {
            libc::free(sequence.as_ptr() as *mut std::ffi::c_void);
        }

        Ok(seq_vec)
    }
}
