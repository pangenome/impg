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

    // A range that starts before the alignment: the target interval must be
    // clipped to what the alignment covers, not copied from the request.
    let approximate_rows = data_rows(
        &query(
            &binary,
            work_dir,
            &[
                "-a", &tpa_arg, "-i", "tpa.impg", "-r", "t:0-1000", "-d", "0", "-o", "bedpe",
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
                "-a", &tpa_arg, "-i", "tpa.impg", "--sequence-files", "seqs.fa", "-r", "t:0-1000",
                "-d", "0", "-o", "bedpe",
            ],
        )
        .stdout,
    );
    assert_eq!(
        intervals(&approximate_rows),
        intervals(&exact_rows),
        "the target interval must be clipped to the alignment, so a range that \
         starts before it must not report homology there"
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

/// Per-file mode: two files whose targets share a local id and coordinates.
///
/// `t1` in `a.tpa` and `y` in `b.tpa` are both local id 1 over 0-1000, with
/// unrelated bases (the names are chosen so the ids match). The thread-local
/// target-sequence cache was keyed on `(id, start, end, strand)` only, so the
/// second query rebuilt `y`'s CIGAR from `t1`'s bases. Each query has its
/// 100 bp deletion in a different place, so the CIGAR in the PAF output shows
/// which bases were used.
#[test]
fn test_per_file_target_cache_does_not_cross_files() {
    let temp_dir = TempDir::new().unwrap();
    let dir = temp_dir.path();
    let binary = get_impg_binary();

    let t1 = dna(1000, 7);
    let q1: Vec<u8> = t1[..200].iter().chain(&t1[300..]).copied().collect();
    let t2 = dna(1000, 11);
    let q2: Vec<u8> = t2[..500].iter().chain(&t2[600..]).copied().collect();
    let mut fasta = File::create(dir.join("seqs.fa")).unwrap();
    for (name, seq) in [("q1", &q1), ("t1", &t1), ("x", &q2), ("y", &t2)] {
        writeln!(fasta, ">{name}\n{}", String::from_utf8_lossy(seq)).unwrap();
    }
    drop(fasta);
    for (name, line) in [
        ("a", "q1\t900\t0\t900\t+\tt1\t1000\t0\t1000\t900\t1000\t60\tcg:Z:200=100D700="),
        ("b", "x\t900\t0\t900\t+\ty\t1000\t0\t1000\t900\t1000\t60\tcg:Z:500=100D400="),
    ] {
        std::fs::write(dir.join(format!("{name}.paf")), format!("{line}\n")).unwrap();
        tpa::paf_to_tpa(
            dir.join(format!("{name}.paf")).to_str().unwrap(),
            dir.join(format!("{name}.tpa")).to_str().unwrap(),
            tpa::CompressionConfig::default(),
        )
        .expect("PAF to TPA conversion failed");
    }
    std::fs::write(dir.join("r.bed"), "t1\t0\t1000\ny\t0\t1000\n").unwrap();

    let out = query(
        &binary,
        dir,
        &[
            "-a", "a.tpa", "b.tpa", "--index-mode", "per-file", "--sequence-files", "seqs.fa",
            "-b", "r.bed", "-d", "0", "-o", "paf", "-t", "1",
        ],
    );
    assert!(out.status.success(), "query failed: {}", String::from_utf8_lossy(&out.stderr));
    let rows = data_rows(&out.stdout);
    let cigar = |q: &str| -> String {
        let row = rows.iter().find(|r| r.starts_with(&format!("{q}\t"))).unwrap_or_else(|| panic!("no {q} row: {rows:?}"));
        row.split('\t').find_map(|f| f.strip_prefix("cg:Z:")).unwrap_or_else(|| panic!("no cg:Z in {row}")).to_string()
    };
    assert_eq!(cigar("q1"), "200=100D700=", "{rows:?}");
    assert_eq!(cigar("x"), "500=100D400=", "y rebuilt from another file's bases: {rows:?}");
}

/// Bases a CIGAR consumes on each axis and its indel totals.
fn cigar_content(cg: &str) -> String {
    let (mut query, mut target, mut ins, mut del, mut len) = (0, 0, 0, 0, 0u64);
    for c in cg.chars() {
        match c.to_digit(10) {
            Some(d) => len = len * 10 + u64::from(d),
            None => {
                match c {
                    '=' | 'X' | 'M' => (query, target) = (query + len, target + len),
                    'I' => (query, ins) = (query + len, ins + len),
                    'D' => (target, del) = (target + len, del + len),
                    other => panic!("unexpected CIGAR op {other} in {cg}"),
                }
                len = 0;
            }
        }
    }
    format!("cg:query={query},target={target},I={ins},D={del}")
}

/// Every tracepoint type must project like the PAF it was converted from.
///
/// Before the fix `FastgaNoDiff` records panicked on `query_deltas[idx]`, and
/// `Variable` put the query axis into the target vector, so any record with
/// an indel projected wrong. `q3` has a 50 bp deletion and a 10 bp insertion,
/// so `Variable` carries both `(len, None)` and `(a, Some(b))` entries.
#[test]
fn test_every_tracepoint_type_matches_the_paf_path() {
    let temp_dir = TempDir::new().unwrap();
    let dir = temp_dir.path();
    let binary = get_impg_binary();

    let t = dna(1200, 42);
    let q3: Vec<u8> = t[..300]
        .iter()
        .chain(&t[350..700])
        .chain(b"ACGTTGCAAC".iter())
        .chain(&t[700..])
        .copied()
        .collect();
    let mut fasta = File::create(dir.join("seqs.fa")).unwrap();
    for (name, seq) in [("q3", &q3), ("t", &t)] {
        writeln!(fasta, ">{name}\n{}", String::from_utf8_lossy(seq)).unwrap();
    }
    drop(fasta);
    std::fs::write(
        dir.join("a.paf"),
        "q3\t1160\t0\t1160\t+\tt\t1200\t0\t1200\t1150\t1210\t60\tcg:Z:300=50D350=10I500=\n",
    )
    .unwrap();
    index_file(&binary, dir, "a.paf", "paf.impg");

    let rows = |aln: &str, idx: &str, region: &str, fmt: &str| -> Vec<String> {
        let out = query(
            &binary,
            dir,
            &["-a", aln, "-i", idx, "--sequence-files", "seqs.fa", "-r", region, "-d", "0", "-o", fmt],
        );
        assert!(out.status.success(), "{aln} {region}: {}", String::from_utf8_lossy(&out.stderr));
        // Drop gi:f/bi:f (TPA takes them from record metadata) and replace the
        // CIGAR by its content: with edit distance, FastGA reconstruction
        // always, an indel can sit anywhere in an equal-score run
        data_rows(&out.stdout)
            .iter()
            .map(|r| {
                r.split('\t')
                    .filter(|f| !f.starts_with("gi:f:") && !f.starts_with("bi:f:"))
                    .map(|f| match f.strip_prefix("cg:Z:") {
                        Some(cg) => cigar_content(cg),
                        None => f.to_string(),
                    })
                    .collect::<Vec<_>>()
                    .join("\t")
            })
            .collect()
    };

    use tpa::TracepointType as T;
    for tp_type in [T::Standard, T::Variable, T::Mixed, T::Fastga, T::FastgaNoDiff] {
        let name = format!("{tp_type:?}").to_lowercase();
        let tpa_file = format!("{name}.tpa");
        tpa::paf_to_tpa(
            dir.join("a.paf").to_str().unwrap(),
            dir.join(&tpa_file).to_str().unwrap(),
            tpa::CompressionConfig { tp_type, max_complexity: 100, ..tpa::CompressionConfig::default() },
        )
        .unwrap_or_else(|e| panic!("{name}: PAF to TPA failed: {e}"));
        index_file(&binary, dir, &tpa_file, &format!("{name}.impg"));
        // Region ends away from the indels, whose placement is not unique
        for region in ["t:0-1200", "t:100-250", "t:400-650", "t:800-1100", "q3:100-1100"] {
            for fmt in ["bedpe", "paf"] {
                let expected = rows("a.paf", "paf.impg", region, fmt);
                assert!(!expected.is_empty(), "{region}: the PAF path printed no row");
                assert_eq!(
                    rows(&tpa_file, &format!("{name}.impg"), region, fmt),
                    expected,
                    "{name} {region} -o {fmt}: TPA and PAF disagree"
                );
            }
        }
    }
}
