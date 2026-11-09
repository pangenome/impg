//! The TPA (tracepoint) path with a real TPA fixture.
//!
//! TPA is the format of the released HPRC r2 alignments, and it had no test at
//! all: no `.tpa` or `.1aln` file appeared in `tests/`. The fixture here is
//! generated from a PAF with a known answer, converted with the same `tpa` crate
//! impg reads with, so a projection that is wrong shows up as a wrong interval
//! rather than as a crash in a worker.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use tempfile::TempDir;

fn get_impg_binary() -> PathBuf {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return PathBuf::from(path);
    }
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    for candidate in [
        manifest_dir.join("target/release/impg"),
        manifest_dir.join("target/debug/impg"),
    ] {
        if candidate.exists() {
            return candidate;
        }
    }
    PathBuf::from("impg")
}

fn run(binary: &Path, dir: &Path, args: &[&str]) -> Output {
    Command::new(binary)
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to run {:?}: {e}", args))
}

/// Deterministic pseudo-random DNA, so the fixture is the same on every run.
fn dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            bases[((state >> 33) % 4) as usize]
        })
        .collect()
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => *other,
        })
        .collect()
}

/// Sequences that match the two PAF records exactly:
///
/// - `q` equals `t[0..1000]`, so record 1 is 1000 matches.
/// - `q2` is the reverse complement of `t[100..900]` with one base changed, so
///   record 2 is `400=1X399=` on the reverse strand.
fn write_sequences(path: &Path) {
    let target = dna(1200, 42);
    let query = target[..1000].to_vec();
    let mut query2 = revcomp(&target[100..900]);
    query2[400] = if query2[400] == b'A' { b'C' } else { b'A' };

    let mut file = File::create(path).unwrap();
    for (name, seq) in [("q", query), ("t", target), ("q2", query2)] {
        writeln!(file, ">{name}").unwrap();
        for chunk in seq.chunks(60) {
            writeln!(file, "{}", String::from_utf8_lossy(chunk)).unwrap();
        }
    }
}

/// The PAF records for those sequences.
fn write_paf(path: &Path) {
    let mut file = File::create(path).unwrap();
    writeln!(
        file,
        "q\t1000\t0\t1000\t+\tt\t1200\t0\t1000\t1000\t1000\t60\tcg:Z:1000="
    )
    .unwrap();
    writeln!(
        file,
        "q2\t800\t0\t800\t-\tt\t1200\t100\t900\t799\t800\t60\tcg:Z:400=1X399="
    )
    .unwrap();
}

/// Write the sequences and the PAF, and convert the PAF to TPA with the same
/// crate impg reads it with. Returns the TPA path.
fn write_fixtures(dir: &Path) -> PathBuf {
    write_sequences(&dir.join("seqs.fa"));
    let paf = dir.join("a.paf");
    write_paf(&paf);
    let tpa = dir.join("a.tpa");
    tpa::paf_to_tpa(
        paf.to_str().unwrap(),
        tpa.to_str().unwrap(),
        tpa::CompressionConfig::default(),
    )
    .expect("PAF to TPA conversion failed");
    tpa
}

fn data_rows(stdout: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(stdout)
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| l.to_string())
        .collect()
}

/// Index a file and return the index path argument.
fn index_file(binary: &Path, dir: &Path, alignment: &str, index: &str) {
    let output = run(binary, dir, &["index", "-a", alignment, "-i", index]);
    assert!(
        output.status.success(),
        "indexing {alignment} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Query with an explicit list of extra arguments.
fn query(binary: &Path, dir: &Path, extra: &[&str]) -> Output {
    let mut args = vec!["query"];
    args.extend_from_slice(extra);
    run(binary, dir, &args)
}

/// The TPA path and the PAF path must agree on the projected intervals.
#[test]
fn test_tpa_index_and_query_match_the_paf_path() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    let tpa_path = write_fixtures(work_dir);
    let tpa_arg = tpa_path.to_string_lossy().into_owned();

    index_file(&binary, work_dir, &tpa_arg, "tpa.impg");
    index_file(&binary, work_dir, "a.paf", "paf.impg");

    for target in ["t:0-1000", "t:100-900", "q2:0-800"] {
        let tpa_rows = data_rows(
            &query(
                &binary,
                work_dir,
                &[
                    "-a", &tpa_arg, "-i", "tpa.impg", "--sequence-files", "seqs.fa", "-r", target,
                    "-d", "0", "-o", "bedpe",
                ],
            )
            .stdout,
        );
        let paf_rows = data_rows(
            &query(
                &binary,
                work_dir,
                &[
                    "-a", "a.paf", "-i", "paf.impg", "--sequence-files", "seqs.fa", "-r", target,
                    "-d", "0", "-o", "bedpe",
                ],
            )
            .stdout,
        );

        assert!(!paf_rows.is_empty(), "{target}: the PAF path printed no row");
        assert_eq!(
            tpa_rows, paf_rows,
            "{target}: the TPA path and the PAF path must agree"
        );
    }

    // The same for paf output, which carries the projected CIGAR.
    let tpa_paf = data_rows(
        &query(
            &binary,
            work_dir,
            &[
                "-a",
                &tpa_arg,
                "-i",
                "tpa.impg",
                "--sequence-files",
                "seqs.fa",
                "-r",
                "q2:0-800",
                "-d",
                "0",
                "-o",
                "paf",
            ],
        )
        .stdout,
    );
    let paf_paf = data_rows(
        &query(
            &binary,
            work_dir,
            &[
                "-a",
                "a.paf",
                "-i",
                "paf.impg",
                "--sequence-files",
                "seqs.fa",
                "-r",
                "q2:0-800",
                "-d",
                "0",
                "-o",
                "paf",
            ],
        )
        .stdout,
    );
    assert_eq!(
        tpa_paf, paf_paf,
        "the projected paf row must not depend on the input format"
    );
}

/// Tracepoint input cannot be turned into CIGARs without the sequences, and
/// the run must say so instead of panicking.
#[test]
fn test_tpa_query_without_sequence_files_is_an_error() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    let tpa_path = write_fixtures(work_dir);
    let tpa_arg = tpa_path.to_string_lossy().into_owned();

    index_file(&binary, work_dir, &tpa_arg, "tpa.impg");

    for format in ["bed", "bedpe", "paf"] {
        let output = query(
            &binary,
            work_dir,
            &["-a", &tpa_arg, "-i", "tpa.impg", "-r", "t:0-1000", "-d", "0", "-o", format],
        );
        assert!(
            !output.status.success(),
            "-o {format} must require --sequence-files for tracepoint input"
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("Sequence files"),
            "-o {format}: the error must say that sequence files are required: {stderr}"
        );
        assert!(
            !stderr.contains("panicked"),
            "-o {format}: it must be an error, not a panic: {stderr}"
        );
    }
}

/// Approximate mode on TPA input: the release path, without sequence files.
#[test]
fn test_tpa_approximate_query_needs_no_sequence_files() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    let tpa_path = write_fixtures(work_dir);
    let tpa_arg = tpa_path.to_string_lossy().into_owned();

    index_file(&binary, work_dir, &tpa_arg, "tpa.impg");

    for format in ["bed", "bedpe"] {
        let approximate = query(
            &binary,
            work_dir,
            &[
                "-a",
                &tpa_arg,
                "-i",
                "tpa.impg",
                "-r",
                "t:0-1000",
                "-d",
                "0",
                "-o",
                format,
                "--approximate",
            ],
        );
        assert!(
            approximate.status.success(),
            "approximate -o {format} on TPA input failed: {}",
            String::from_utf8_lossy(&approximate.stderr)
        );
        assert!(
            !data_rows(&approximate.stdout).is_empty(),
            "approximate -o {format} printed nothing"
        );
    }

    // A range inside the alignment: the projected intervals must be identical
    // to the CIGAR path, which needs the sequences.
    let approximate_rows = data_rows(
        &query(
            &binary,
            work_dir,
            &[
                "-a", &tpa_arg, "-i", "tpa.impg", "-r", "t:100-900", "-d", "0", "-o", "bedpe",
                "--approximate",
            ],
        )
        .stdout,
    );
    let exact_rows = data_rows(
        &query(
            &binary,
            work_dir,
            &[
                "-a", &tpa_arg, "-i", "tpa.impg", "--sequence-files", "seqs.fa", "-r", "t:100-900",
                "-d", "0", "-o", "bedpe",
            ],
        )
        .stdout,
    );
    let intervals = |rows: &[String]| -> Vec<String> {
        rows.iter()
            .map(|row| {
                let f: Vec<&str> = row.split('\t').collect();
                format!("{}:{}-{}|{}:{}-{}", f[0], f[1], f[2], f[3], f[4], f[5])
            })
            .collect()
    };
    assert_eq!(
        intervals(&approximate_rows),
        intervals(&exact_rows),
        "inside the alignment the tracepoint projection must give the same intervals"
    );

    // The identity is not the same: for a Standard TPA record the tracepoint
    // path estimates it from the trace spacing, so a perfect alignment does not
    // read 1. It is a lower bound, never higher than the CIGAR path.
    let identity = |row: &str| -> f64 {
        let field = row
            .split('\t')
            .find(|f| f.starts_with("gi:f:"))
            .expect("bedpe row without gi:f");
        field.rsplit(':').next().unwrap().parse().unwrap()
    };
    for (approximate_row, exact_row) in approximate_rows.iter().zip(exact_rows.iter()) {
        let (a, e) = (identity(approximate_row), identity(exact_row));
        assert!(
            a <= e + 1e-9,
            "the tracepoint identity must not exceed the CIGAR identity: {a} > {e}"
        );
        assert!(
            a >= 0.9,
            "the tracepoint identity on an exact alignment is far too low: {a}"
        );
    }
}

