pub mod align;
pub mod genome_infer;
pub mod genotype;
pub mod graph;
pub mod infer;
pub mod lace;
pub mod partition;
pub mod refine;
pub mod render;
pub mod similarity;
pub mod syng2gfa;

use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Floor for the haplotype-derived FastGA seed frequency, as the sweepga
/// help states. FastGA drops seeds with frequency >= -f, so -f equal to
/// the haplotype count drops every k-mer shared by all haplotypes.
pub const MIN_FASTGA_FREQUENCY: usize = 10;

/// wfmash -n for raw (unsparsified) all-vs-all: the larger side of
/// --num-mappings, or every sequence when unbounded or unfiltered.
pub fn wfmash_num_mappings(num_mappings: &str, no_filter: bool, sequence_count: usize) -> usize {
    let (_, query_limit, target_limit) = sweepga::library_api::parse_filter_mode(num_mappings);
    match (no_filter, query_limit, target_limit) {
        (false, Some(q), Some(t)) => q.max(t).max(1),
        _ => sequence_count.max(1),
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
pub fn count_sequences_and_genomes(fasta_files: &[String]) -> io::Result<(usize, usize)> {
    use std::fs::File;

    let paths: Vec<&Path> = fasta_files.iter().map(Path::new).collect();
    let genome_count = sweepga::pansn::count_haplotypes(&paths)
        .map_err(|e| io::Error::other(format!("PanSN haplotype count failed: {e}")))?;

    let mut seq_count = 0usize;
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
            }
        }
    }

    Ok((seq_count, genome_count))
}
