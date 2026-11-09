//! `impg refine --approximate` on tracepoint input.
//!
//! Approximate mode interpolates the query interval from the tracepoint
//! boundaries and never reads a CIGAR, so refine must accept it without
//! `--sequence-files`, with and without `--transitive`.

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

fn sibling(binary: &Path, name: &str) -> PathBuf {
    binary.parent().unwrap().join(name)
}

fn run(binary: &Path, dir: &Path, args: &[&str]) -> Output {
    Command::new(binary)
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to run {} {:?}: {e}", binary.display(), args))
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

/// Build a `.1aln` file from the yeast chrV fixture.
///
/// `impg align --format 1aln` writes PAF for a single input FASTA, so the file
/// is converted with `PAFtoALN`, which ships next to the impg binary. Returns
/// None when the conversion tool is not present.
fn build_onealn(binary: &Path, dir: &Path) -> Option<PathBuf> {
    let fasta = dir.join("seqs.fa");
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/test_data/yeast.chrV.fa.gz");
    decompress_gz(&fixture, &fasta);

    let output = run(
        binary,
        dir,
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
    assert!(
        output.status.success(),
        "impg align failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // impg may write either the .1aln directly or a PAF to convert.
    let direct = dir.join("aln/alignments.1aln");
    if direct.exists() {
        return Some(direct);
    }

    let paf = dir.join("aln/alignments.paf");
    assert!(paf.exists(), "impg align wrote neither a .1aln nor a PAF");

    let paf_to_aln = sibling(binary, "PAFtoALN");
    if !paf_to_aln.exists() {
        eprintln!("skipping: {} is not next to the impg binary", paf_to_aln.display());
        return None;
    }
    let output = Command::new(&paf_to_aln)
        .current_dir(dir)
        .args(["-T2", "aln/alignments.paf", "seqs.fa"])
        .output()
        .expect("failed to run PAFtoALN");
    assert!(
        output.status.success(),
        "PAFtoALN failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let converted = dir.join("aln/alignments.1aln");
    assert!(
        converted.exists(),
        "PAFtoALN did not write aln/alignments.1aln"
    );
    Some(converted)
}

/// One refine output row: `#chrom start end name original.support new.support left.bp right.bp`.
#[derive(Debug, Clone, PartialEq)]
struct Row {
    chrom: String,
    start: i64,
    end: i64,
    name: String,
    new_support: i64,
    left_bp: i64,
    right_bp: i64,
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
                left_bp: f[6].parse().unwrap(),
                right_bp: f[7].parse().unwrap(),
            }
        })
        .collect()
}

/// Deterministic BED of windows spread over the reference, so the comparison
/// covers interior loci, the chromosome start and both strands of the
/// alignments that reach the reference.
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
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();

    let Some(onealn) = build_onealn(&binary, work_dir) else {
        return;
    };
    let onealn_arg = onealn
        .strip_prefix(work_dir)
        .unwrap_or(&onealn)
        .to_string_lossy()
        .into_owned();

    let windows = 28;
    write_target_bed(&work_dir.join("targets.bed"), windows);

    let exact = run(
        &binary,
        work_dir,
        &[
            "refine", "-a", &onealn_arg, "--sequence-files", "seqs.fa", "-b", "targets.bed", "-d",
            "1000", "--span-bp", "100",
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
            "refine", "-a", &onealn_arg, "-b", "targets.bed", "-d", "1000", "--span-bp", "100",
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

    // Every target is answered by both modes, in the same order.
    assert_eq!(
        exact_rows.len(),
        windows,
        "the exact run lost targets: {:?}",
        exact_rows.len()
    );
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

        // A coarser projection may pick a different flank step, never more.
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

    // The agreement this fixture shows today is 26/28 rows. Keep a margin so a
    // small change in the alignment does not fail the test, but catch a real
    // regression in the projection.
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
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();

    let Some(onealn) = build_onealn(&binary, work_dir) else {
        return;
    };
    let onealn_arg = onealn
        .strip_prefix(work_dir)
        .unwrap_or(&onealn)
        .to_string_lossy()
        .into_owned();

    // The CIGAR path cannot run on tracepoint input without the sequences.
    let exact = run(
        &binary,
        work_dir,
        &[
            "refine",
            "-a",
            &onealn_arg,
            "-r",
            "S288C#1#chrV:100000-101000",
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
            &onealn_arg,
            "-r",
            "S288C#1#chrV:100000-101000",
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
}

#[test]
fn test_refine_approximate_with_transitive_still_works() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();

    let Some(onealn) = build_onealn(&binary, work_dir) else {
        return;
    };
    let onealn_arg = onealn
        .strip_prefix(work_dir)
        .unwrap_or(&onealn)
        .to_string_lossy()
        .into_owned();

    let output = run(
        &binary,
        work_dir,
        &[
            "refine",
            "-a",
            &onealn_arg,
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
