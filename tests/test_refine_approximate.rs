//! `impg refine --approximate` on tracepoint input.
//!
//! Approximate mode interpolates the query interval from the tracepoint
//! boundaries and never reads a CIGAR, so refine must accept it without
//! `--sequence-files`, with and without `--transitive`.
//!
//! The `.1aln` fixture costs one alignment run, so all tests in this file share
//! it through `FIXTURE_DIR`. Each test still gets its own directory for output.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::OnceLock;
use tempfile::TempDir;

/// Shared fixture directory, dropped when the test process exits.
static FIXTURE_DIR: OnceLock<Option<TempDir>> = OnceLock::new();

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

fn decompress_gz(gz_path: &Path, out_path: &Path) {
    let output = Command::new("zcat")
        .arg(gz_path)
        .output()
        .expect("zcat is required to unpack the test FASTA");
    assert!(output.status.success(), "zcat failed on {gz_path:?}");
    File::create(out_path)
        .unwrap()
        .write_all(&output.stdout)
        .unwrap();
}

/// Build the `.1aln` fixture once: yeast chrV, 7 strains.
///
/// `impg align --format 1aln` writes PAF for a single input FASTA, so the file is
/// converted with `PAFtoALN`, which ships next to the impg binary.
fn fixture() -> Option<&'static Path> {
    FIXTURE_DIR
        .get_or_init(|| {
            let dir = TempDir::new().ok()?;
            let binary = get_impg_binary();
            let work = dir.path();

            let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/test_data/yeast.chrV.fa.gz");
            decompress_gz(&fixture, &work.join("seqs.fa"));

            let output = run(
                &binary,
                work,
                &[
                    "align",
                    "--sequence-files",
                    "seqs.fa",
                    "-o",
                    "aln",
                    "--format",
                    "1aln",
                    "-t",
                    "2",
                ],
            );
            if !output.status.success() {
                eprintln!(
                    "skipping refine tests: impg align failed: {}",
                    String::from_utf8_lossy(&output.stderr)
                );
                return None;
            }

            if !work.join("aln/alignments.1aln").exists() {
                let paf_to_aln = binary.parent().unwrap().join("PAFtoALN");
                if !paf_to_aln.exists() {
                    eprintln!(
                        "skipping refine tests: {} is not next to the impg binary",
                        paf_to_aln.display()
                    );
                    return None;
                }
                let output = Command::new(&paf_to_aln)
                    .current_dir(work)
                    .args(["-T2", "aln/alignments.paf", "seqs.fa"])
                    .output()
                    .expect("failed to run PAFtoALN");
                if !output.status.success() {
                    eprintln!(
                        "skipping refine tests: PAFtoALN failed: {}",
                        String::from_utf8_lossy(&output.stderr)
                    );
                    return None;
                }
            }

            if work.join("aln/alignments.1aln").exists() {
                Some(dir)
            } else {
                eprintln!("skipping refine tests: no .1aln was produced");
                None
            }
        })
        .as_ref()
        .map(|dir| dir.path())
}

/// One refine output row: `#chrom start end name original.support new.support left.bp right.bp`.
#[derive(Debug, Clone, PartialEq)]
struct Row {
    chrom: String,
    start: i64,
    end: i64,
    name: String,
    new_support: i64,
}

fn parse_rows(stdout: &[u8]) -> Vec<Row> {
    String::from_utf8_lossy(stdout)
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert!(f.len() >= 8, "unexpected refine row: {l}");
            Row {
                chrom: f[0].to_string(),
                start: f[1].parse().unwrap(),
                end: f[2].parse().unwrap(),
                name: f[3].to_string(),
                new_support: f[5].parse().unwrap(),
            }
        })
        .collect()
}

/// Deterministic BED of windows spread over the reference, so the comparison
/// covers interior loci and the chromosome start.
fn write_target_bed(path: &Path, windows: usize) {
    let mut file = File::create(path).unwrap();
    let mut start: i64 = 0;
    for i in 0..windows {
        writeln!(file, "S288C#1#chrV\t{start}\t{}\tw{}", start + 2000, i + 1).unwrap();
        start += 18_000 + ((i as i64 * 997) % 5_000);
    }
}

/// Approximate mode interpolates the query interval from the tracepoint
/// boundaries, so it is close to the CIGAR path, not identical to it. This test
/// pins the agreement: how many rows match exactly, how far a coordinate may
/// move, and how far the support may move.
#[test]
fn test_refine_approximate_matches_the_cigar_path_on_many_windows() {
    let Some(fixture) = fixture() else { return };
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    let onealn = fixture.join("aln/alignments.1aln");
    let seqs = fixture.join("seqs.fa");

    let windows = 28;
    write_target_bed(&work_dir.join("targets.bed"), windows);

    let exact = run(
        &binary,
        work_dir,
        &[
            "refine",
            "-a",
            &onealn.to_string_lossy(),
            "--sequence-files",
            &seqs.to_string_lossy(),
            "-b",
            "targets.bed",
            "-d",
            "1000",
            "--span-bp",
            "100",
        ],
    );
    assert!(
        exact.status.success(),
        "exact refine failed: {}",
        String::from_utf8_lossy(&exact.stderr)
    );

    let approximate = run(
        &binary,
        work_dir,
        &[
            "refine",
            "-a",
            &onealn.to_string_lossy(),
            "-b",
            "targets.bed",
            "-d",
            "1000",
            "--span-bp",
            "100",
            "--approximate",
        ],
    );
    assert!(
        approximate.status.success(),
        "approximate refine failed without --sequence-files: {}",
        String::from_utf8_lossy(&approximate.stderr)
    );

    let exact_rows = parse_rows(&exact.stdout);
    let approximate_rows = parse_rows(&approximate.stdout);

    assert_eq!(exact_rows.len(), windows, "the exact run lost targets");
    assert_eq!(
        exact_rows.iter().map(|r| r.name.clone()).collect::<Vec<_>>(),
        approximate_rows
            .iter()
            .map(|r| r.name.clone())
            .collect::<Vec<_>>(),
        "the two modes must answer the same targets in the same order"
    );

    let extension_step = 1000;
    let mut identical = 0usize;
    let mut max_start_shift = 0i64;
    let mut max_support_delta = 0i64;

    for (e, a) in exact_rows.iter().zip(approximate_rows.iter()) {
        assert_eq!(e.chrom, a.chrom, "{}: chromosome changed", e.name);
        assert_eq!(
            e.end, a.end,
            "{}: approximate changed the right boundary by {} bp",
            e.name,
            a.end - e.end
        );

        let start_shift = (a.start - e.start).abs();
        let support_delta = (a.new_support - e.new_support).abs();
        max_start_shift = max_start_shift.max(start_shift);
        max_support_delta = max_support_delta.max(support_delta);

        assert!(
            start_shift <= extension_step,
            "{}: approximate moved the left boundary by {start_shift} bp, more than one extension step ({extension_step} bp)",
            e.name
        );
        assert!(
            start_shift % extension_step == 0,
            "{}: approximate moved the left boundary by {start_shift} bp, which is not a whole extension step",
            e.name
        );

        if *e == *a {
            identical += 1;
        }
    }

    let agreed = 100.0 * identical as f64 / windows as f64;
    eprintln!(
        "approximate vs CIGAR: {identical}/{windows} rows identical ({agreed:.1}%), \
         max left shift {max_start_shift} bp, max support delta {max_support_delta}"
    );

    assert!(
        agreed >= 85.0,
        "approximate mode agreed on only {identical}/{windows} rows ({agreed:.1}%)"
    );
    assert!(
        max_support_delta <= 2,
        "approximate mode changed the support by {max_support_delta}, more than this test allows"
    );
}

/// The point of approximate mode: no sequence files and no CIGAR read.
#[test]
fn test_refine_approximate_needs_no_sequence_files() {
    let Some(fixture) = fixture() else { return };
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    let onealn = fixture.join("aln/alignments.1aln");
    let seqs = fixture.join("seqs.fa");
    let target = "S288C#1#chrV:100000-101000";

    // The CIGAR path cannot run on tracepoint input without the sequences.
    let exact = run(
        &binary,
        work_dir,
        &[
            "refine",
            "-a",
            &onealn.to_string_lossy(),
            "-r",
            target,
            "-d",
            "1000",
            "--span-bp",
            "100",
        ],
    );
    assert!(
        !exact.status.success(),
        "the CIGAR path must still require --sequence-files for tracepoint input"
    );
    assert!(
        String::from_utf8_lossy(&exact.stderr).contains("Sequence files"),
        "the CIGAR path must say why it needs sequences"
    );

    let approximate = run(
        &binary,
        work_dir,
        &[
            "refine",
            "-a",
            &onealn.to_string_lossy(),
            "-r",
            target,
            "-d",
            "1000",
            "--span-bp",
            "100",
            "--approximate",
        ],
    );
    assert!(
        approximate.status.success(),
        "approximate mode must not need --sequence-files: {}",
        String::from_utf8_lossy(&approximate.stderr)
    );

    // Same answer as the CIGAR path on a clean locus.
    let exact = run(
        &binary,
        work_dir,
        &[
            "refine",
            "-a",
            &onealn.to_string_lossy(),
            "--sequence-files",
            &seqs.to_string_lossy(),
            "-r",
            target,
            "-d",
            "1000",
            "--span-bp",
            "100",
        ],
    );
    assert!(exact.status.success(), "exact refine failed");
    assert_eq!(
        parse_rows(&exact.stdout),
        parse_rows(&approximate.stdout),
        "the two modes must agree on a clean locus"
    );
}

#[test]
fn test_refine_approximate_with_transitive_still_works() {
    let Some(fixture) = fixture() else { return };
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    let onealn = fixture.join("aln/alignments.1aln");

    let output = run(
        &binary,
        work_dir,
        &[
            "refine",
            "-a",
            &onealn.to_string_lossy(),
            "-r",
            "S288C#1#chrV:100000-101000",
            "-d",
            "1000",
            "--span-bp",
            "100",
            "--approximate",
            "--transitive",
        ],
    );
    assert!(
        output.status.success(),
        "transitive approximate refine failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        !parse_rows(&output.stdout).is_empty(),
        "transitive approximate refine printed no data row"
    );
}

#[test]
fn test_refine_approximate_rejects_paf_input() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();

    let mut paf = File::create(work_dir.join("a.paf")).unwrap();
    writeln!(
        paf,
        "q\t1000\t0\t500\t+\tt\t1000\t0\t500\t500\t500\t60\tcg:Z:500="
    )
    .unwrap();

    let output = run(
        &binary,
        work_dir,
        &[
            "refine",
            "-a",
            "a.paf",
            "-r",
            "q:0-500",
            "-d",
            "1000",
            "--approximate",
        ],
    );
    assert!(!output.status.success(), "PAF input must be rejected");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("tracepoint alignment files"),
        "the error must name the accepted input: {stderr}"
    );
}
