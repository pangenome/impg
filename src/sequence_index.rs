use std::io;
use std::path::Path;

#[cfg(feature = "agc")]
use crate::agc_index::AgcIndex;
use crate::faidx::FastaIndex;

// Trait for sequence fetching from different sources
pub trait SequenceIndex {
    fn fetch_sequence(&self, seq_name: &str, start: i32, end: i32) -> io::Result<Vec<u8>>;
}

// Enum to hold either FASTA or AGC index
#[derive(Debug)]
pub enum UnifiedSequenceIndex {
    Fasta(FastaIndex),
    #[cfg(feature = "agc")]
    Agc(AgcIndex),
}

impl UnifiedSequenceIndex {
    pub fn from_files(files: &[String]) -> io::Result<Self> {
        if files.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "No input files provided",
            ));
        }

        // Check if all files have the same extension
        let get_full_extension = |path: &str| -> String {
            let p = Path::new(path);
            let file_name = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
            
            // Handle compound extensions like .fa.gz, .fasta.gz, .fna.gz
            if file_name.ends_with(".fa.gz") {
                "fa.gz".to_string()
            } else if file_name.ends_with(".fasta.gz") {
                "fasta.gz".to_string()
            } else if file_name.ends_with(".fna.gz") {
                "fna.gz".to_string()
            } else {
                p.extension()
                    .and_then(|s| s.to_str())
                    .unwrap_or("")
                    .to_string()
            }
        };

        let first_ext = get_full_extension(&files[0]);

        let all_same_type = files.iter().all(|f| {
            get_full_extension(f) == first_ext
        });

        if !all_same_type {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Mixed file types not supported. All files must be either .fa/.fasta/.fna/.fa.gz/.fasta.gz/.fna.gz or .agc",
            ));
        }

        match first_ext.as_str() {
            "fa" | "fasta" | "fna" | "fa.gz" | "fasta.gz" | "fna.gz" => {
                let index = FastaIndex::build_from_files(files)?;
                Ok(UnifiedSequenceIndex::Fasta(index))
            }
            #[cfg(feature = "agc")]
            "agc" => {
                let index = AgcIndex::build_from_files(files)?;
                Ok(UnifiedSequenceIndex::Agc(index))
            }
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Unsupported file extension: {}", first_ext),
            )),
        }
    }
}

impl SequenceIndex for UnifiedSequenceIndex {
    fn fetch_sequence(&self, seq_name: &str, start: i32, end: i32) -> io::Result<Vec<u8>> {
        match self {
            UnifiedSequenceIndex::Fasta(index) => index.fetch_sequence(seq_name, start, end),
            #[cfg(feature = "agc")]
            UnifiedSequenceIndex::Agc(index) => index.fetch_sequence(seq_name, start, end),
        }
    }
}

// Implement for FastaIndex
impl SequenceIndex for FastaIndex {
    fn fetch_sequence(&self, seq_name: &str, start: i32, end: i32) -> io::Result<Vec<u8>> {
        FastaIndex::fetch_sequence(self, seq_name, start, end)
    }
}

// Implement for AgcIndex
#[cfg(feature = "agc")]
impl SequenceIndex for AgcIndex {
    fn fetch_sequence(&self, seq_name: &str, start: i32, end: i32) -> io::Result<Vec<u8>> {
        AgcIndex::fetch_sequence(self, seq_name, start, end)
    }
}