pub mod align;
pub mod graph;
pub mod lace;
pub mod partition;
pub mod refine;
pub mod similarity;

use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

use sweepga::aligner::Aligner;
use sweepga::fastga_integration::FastGAIntegration;
use sweepga::paf_filter::{FilterConfig, FilterMode, ScoringFunction};
use sweepga::wfmash_integration::WfmashIntegration;

/// Create an aligner backend based on the name ("wfmash" or "fastga").
pub fn create_aligner(
    aligner_name: &str,
    kmer_frequency: usize,
    num_threads: usize,
    min_aln_length: u64,
    map_pct_identity: Option<String>,
    temp_dir: Option<String>,
) -> io::Result<Box<dyn Aligner>> {
    create_aligner_adaptive(
        aligner_name,
        kmer_frequency,
        num_threads,
        min_aln_length,
        map_pct_identity,
        temp_dir,
        None,
        None,
        None,
        None,
        None,
    )
}

/// Create an aligner backend with adaptive wfmash parameters.
///
/// `segment_length`: wfmash segment length (-s). None = adaptive.
/// `avg_seq_len`: average input sequence length, used to adapt parameters.
/// `num_mappings`: wfmash -n flag. None = wfmash default (1).
pub fn create_aligner_adaptive(
    aligner_name: &str,
    kmer_frequency: usize,
    num_threads: usize,
    min_aln_length: u64,
    map_pct_identity: Option<String>,
    temp_dir: Option<String>,
    segment_length: Option<u64>,
    avg_seq_len: Option<u64>,
    sparsify: Option<f64>,
    num_mappings: Option<usize>,
    pairs_file: Option<PathBuf>,
) -> io::Result<Box<dyn Aligner>> {
    match aligner_name {
        "wfmash" => {
            let block_len = if min_aln_length > 0 {
                Some(min_aln_length)
            } else {
                None
            };
            let wfmash = WfmashIntegration::adaptive(
                num_threads,
                block_len,
                map_pct_identity,
                temp_dir,
                segment_length,
                avg_seq_len,
                sparsify,
                num_mappings,
                pairs_file,
            )
            .map_err(|e| io::Error::other(format!("Failed to create wfmash aligner: {e}")))?;
            Ok(Box::new(wfmash))
        }
        "fastga" => Ok(Box::new(FastGAIntegration::new(
            Some(kmer_frequency),
            num_threads,
            min_aln_length,
            temp_dir,
        ))),
        _ => Err(io::Error::other(format!(
            "Unknown aligner: {}. Valid options: wfmash, fastga",
            aligner_name
        ))),
    }
}

/// Round a value to a "nice" multiple based on its magnitude:
/// ≤500 → multiple of 50, ≤1000 → 100, ≤3000 → 200, >3000 → 500.
pub fn round_nice(v: u64) -> u64 {
    if v == 0 {
        return 0;
    }
    let step = if v <= 500 {
        50
    } else if v <= 1000 {
        100
    } else if v <= 3000 {
        200
    } else {
        500
    };
    ((v + step / 2) / step * step).max(step)
}

/// Parse filter mode string (e.g., "1:1", "many:many", "5:3") into FilterMode and limits.
pub fn parse_filter_mode(s: &str) -> (FilterMode, Option<usize>, Option<usize>) {
    let s_lower = s.to_lowercase();
    if s_lower == "many:many" || s_lower == "n:n" {
        return (FilterMode::ManyToMany, None, None);
    }

    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 2 {
        return (FilterMode::OneToOne, Some(1), Some(1));
    }

    let query_max = if parts[0] == "many" || parts[0] == "n" {
        None
    } else {
        parts[0].parse().ok()
    };

    let target_max = if parts[1] == "many" || parts[1] == "n" {
        None
    } else {
        parts[1].parse().ok()
    };

    match (query_max, target_max) {
        (Some(1), Some(1)) => (FilterMode::OneToOne, Some(1), Some(1)),
        (Some(1), _) => (FilterMode::OneToMany, Some(1), target_max),
        (_, Some(1)) => (FilterMode::OneToMany, query_max, Some(1)),
        (Some(q), Some(t)) => (FilterMode::ManyToMany, Some(q), Some(t)),
        (Some(q), None) => (FilterMode::ManyToMany, Some(q), None),
        (None, Some(t)) => (FilterMode::ManyToMany, None, Some(t)),
        (None, None) => (FilterMode::ManyToMany, None, None),
    }
}

/// Parameters for PAF alignment filtering, shared across all commands.
pub struct FilterParams {
    pub num_mappings: String,
    pub scaffold_jump: u64,
    pub scaffold_mass: u64,
    pub scaffold_filter: String,
    pub overlap: f64,
    pub min_identity: f64,
    pub scaffold_dist: u64,
    pub min_map_length: u64,
}

/// Build a `FilterConfig` with adaptive scaffold parameters.
///
/// Scaffold defaults (scaffold_mass=10kb, scaffold_jump=50kb) are tuned for
/// whole-genome alignments. For short sequences (e.g. 1kb excerpts from
/// `impg query -o fasta`), these thresholds would filter out every alignment.
/// This function clamps them based on `avg_seq_len`.
pub fn build_filter_config(params: &FilterParams, avg_seq_len: u64) -> FilterConfig {
    let (mapping_mode, mapping_per_query, mapping_per_target) =
        parse_filter_mode(&params.num_mappings);
    let (scaffold_mode, scaffold_per_query, scaffold_per_target) =
        parse_filter_mode(&params.scaffold_filter);

    let scaffold_mass = if avg_seq_len > 0 {
        round_nice(params.scaffold_mass.min(avg_seq_len * 3 / 5))
    } else {
        params.scaffold_mass
    };
    let scaffold_jump = if avg_seq_len > 0 {
        params.scaffold_jump.min(avg_seq_len * 10)
    } else {
        params.scaffold_jump
    };

    FilterConfig {
        chain_gap: 0,
        min_block_length: params.min_map_length,
        mapping_filter_mode: mapping_mode,
        mapping_max_per_query: mapping_per_query,
        mapping_max_per_target: mapping_per_target,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: scaffold_mode,
        scaffold_max_per_query: scaffold_per_query,
        scaffold_max_per_target: scaffold_per_target,
        overlap_threshold: params.overlap,
        sparsity: 1.0,
        no_merge: true,
        scaffold_gap: scaffold_jump,
        min_scaffold_length: scaffold_mass,
        scaffold_overlap_threshold: 0.5,
        scaffold_max_deviation: params.scaffold_dist,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: params.min_identity,
        min_scaffold_identity: params.min_identity,
    }
}

/// Resolve a list of files from either a direct `Vec` or a list file.
///
/// - If `files` is non-empty and `list` is `None`: validate all files exist and return them.
/// - If `files` is empty and `list` is `Some(path)`: read paths from the list file.
/// - If both are empty/None: return an error.
/// - If both are provided: return an error.
///
/// `label` is used in error messages (e.g. "FASTA", "GFA", "VCF").
pub fn resolve_file_list(
    files: Vec<String>,
    list: Option<String>,
    label: &str,
) -> io::Result<Vec<String>> {
    match (files.is_empty(), list) {
        (false, None) => {
            for file in &files {
                if !Path::new(file).exists() {
                    return Err(io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("{label} file '{file}' not found"),
                    ));
                }
            }
            Ok(files)
        }
        (true, Some(list_file)) => {
            let f = std::fs::File::open(&list_file)?;
            let reader = BufReader::new(f);
            let mut out = Vec::new();
            for line in reader.lines() {
                let line = line?;
                let trimmed = line.trim();
                if !trimmed.is_empty() && !trimmed.starts_with('#') {
                    if !Path::new(trimmed).exists() {
                        return Err(io::Error::new(
                            io::ErrorKind::NotFound,
                            format!("{label} file '{trimmed}' not found"),
                        ));
                    }
                    out.push(trimmed.to_string());
                }
            }
            if out.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("No valid {label} files found in list file: {list_file}"),
                ));
            }
            Ok(out)
        }
        (true, None) => {
            let label_lower = label.to_lowercase();
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Either --{label_lower}-files or --{label_lower}-list must be provided"),
            ))
        }
        (false, Some(_)) => {
            let label_lower = label.to_lowercase();
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Cannot specify both --{label_lower}-files and --{label_lower}-list"),
            ))
        }
    }
}

/// Count sequences and genomes across FASTA files (for k-mer frequency calculation).
///
/// Uses PanSN naming convention: for names like `SAMPLE#HAPLOTYPE#CONTIG`, genomes
/// are unique `SAMPLE#HAPLOTYPE` prefixes.
/// Returns `(num_sequences, num_genomes)`.
pub fn count_sequences_and_genomes(fasta_files: &[String]) -> io::Result<(usize, usize)> {
    use std::collections::HashSet;
    use std::fs::File;

    let mut seq_count = 0;
    let mut genome_prefixes: HashSet<String> = HashSet::new();

    for path in fasta_files {
        let file = File::open(path)?;
        let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
            io::Error::other(format!("Failed to open reader for '{}': {}", path, e))
        })?;
        let reader = BufReader::new(reader);

        for line in reader.lines() {
            let line: String = line?;
            if line.starts_with('>') {
                seq_count += 1;
                let name = line
                    .strip_prefix('>')
                    .unwrap_or("")
                    .split_whitespace()
                    .next()
                    .unwrap_or("");
                let parts: Vec<&str> = name.split('#').collect();
                let prefix = if parts.len() >= 2 {
                    format!("{}#{}", parts[0], parts[1])
                } else {
                    name.to_string()
                };
                genome_prefixes.insert(prefix);
            }
        }
    }

    let genome_count = genome_prefixes.len().max(1);
    Ok((seq_count, genome_count))
}
