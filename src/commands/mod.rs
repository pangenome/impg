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

use std::fs::File;
use std::io::{self, BufRead, BufWriter, Read, Write};
use std::path::Path;

/// Create missing parent directories for `path` if non-empty.
pub fn ensure_parent_dir(path: &Path) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    Ok(())
}

/// Read all UTF-8 text from `input`, or stdin when `input == "-"`.
pub fn read_text_from_path_or_stdin(input: &str) -> io::Result<String> {
    let mut text = String::new();
    if input == "-" {
        io::stdin().read_to_string(&mut text)?;
    } else {
        File::open(input)?.read_to_string(&mut text)?;
    }
    Ok(text)
}

/// Write UTF-8 `text` to `output`, or stdout when `output == "-"`.
pub fn write_text_to_path_or_stdout(
    output: &str,
    text: &str,
    create_parent: bool,
) -> io::Result<()> {
    if output == "-" {
        let stdout = io::stdout();
        let mut out = BufWriter::with_capacity(1024 * 1024, stdout.lock());
        out.write_all(text.as_bytes())?;
        out.flush()?;
    } else {
        let output_path = Path::new(output);
        if create_parent {
            ensure_parent_dir(output_path)?;
        }
        let mut out = BufWriter::with_capacity(1024 * 1024, File::create(output_path)?);
        out.write_all(text.as_bytes())?;
        out.flush()?;
    }
    Ok(())
}

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

/// Paths in a list file: trimmed lines, skipping blank and `#` lines.
pub fn read_list_file(path: &str) -> io::Result<Vec<String>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| io::Error::new(e.kind(), format!("failed to read list file '{path}': {e}")))?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_string)
        .collect())
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
            let out = read_list_file(&list_file)?;
            if let Some(missing) = out.iter().find(|f| !Path::new(f).exists()) {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("{label} file '{missing}' not found"),
                ));
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
    let paths: Vec<&Path> = fasta_files.iter().map(Path::new).collect();
    let genome_count = sweepga::pansn::count_haplotypes(&paths)
        .map_err(|e| io::Error::other(format!("PanSN haplotype count failed: {e}")))?;

    let mut seq_count = 0usize;
    for path in fasta_files {
        for line in crate::fasta::open_reader(Path::new(path))?.lines() {
            let line: String = line?;
            if line.starts_with('>') {
                seq_count += 1;
            }
        }
    }

    Ok((seq_count, genome_count))
}

#[cfg(test)]
mod list_file_tests {
    #[test]
    fn read_list_file_skips_blank_and_comment_lines() {
        let dir = tempfile::tempdir().unwrap();
        let list = dir.path().join("files.txt");
        std::fs::write(&list, "  a.fa  \n\n# note\nb.fa\n   \n").unwrap();
        assert_eq!(
            super::read_list_file(list.to_str().unwrap()).unwrap(),
            ["a.fa", "b.fa"]
        );
        let err = super::read_list_file("/nonexistent/list.txt").unwrap_err();
        assert!(err.to_string().contains("/nonexistent/list.txt"), "{err}");
    }
}
