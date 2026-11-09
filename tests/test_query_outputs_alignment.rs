//! `impg query -o bedpe` and `-o paf` on the alignment backend.
//!
//! These two formats had no test on the alignment backend at all. They carry
//! the projected CIGAR, so they are what a script uses to build a small graph
//! or a BEDPE of projected intervals. The fixtures below are exact: every CIGAR
//! consumes exactly the lengths the PAF record declares.

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

/// One fixture with three exact records:
///
/// - `q2`/`t2`: 1000 bp of pure match, so every projection is a diagonal.
/// - `A`/`B`: `B` has 10 bases that `A` does not (`10D`), so the projected
///   query range is shorter than the target range and the CIGAR keeps the `D`.
/// - `top`/`bot`: `bot` is on the reverse strand.
fn write_fixture(dir: &Path) -> PathBuf {
    let paf = dir.join("a.paf");
    let mut file = File::create(&paf).unwrap();
    writeln!(
        file,
        "t2\t1000\t0\t1000\t+\tq2\t1000\t0\t1000\t1000\t1000\t60\tcg:Z:1000="
    )
    .unwrap();
    // 400 + 10 + 580 target bases, 400 + 580 query bases.
    writeln!(
        file,
        "A\t980\t0\t980\t+\tB\t990\t0\t990\t980\t990\t60\tcg:Z:400=10D580="
    )
    .unwrap();
    writeln!(
        file,
        "top\t1000\t2000\t2500\t-\tbot\t1000\t100\t600\t500\t500\t60\tcg:Z:500="
    )
    .unwrap();
    paf
}

fn write_index(binary: &Path, dir: &Path) {
    let output = run(binary, dir, &["index", "-a", "a.paf", "-i", "a.impg"]);
    assert!(
        output.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Run a query and return its data rows, each split into tab-separated fields.
fn query_rows(binary: &Path, dir: &Path, target: &str, format: &str) -> Vec<Vec<String>> {
    let output = run(
        binary,
        dir,
        &[
            "query", "-a", "a.paf", "-i", "a.impg", "-r", target, "-d", "0", "-o", format,
        ],
    );
    assert!(
        output.status.success(),
        "query {target} -o {format} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| l.split('\t').map(|f| f.to_string()).collect())
        .collect()
}

fn number(field: &str) -> f64 {
    field
        .split(':')
        .next_back()
        .and_then(|v| v.parse().ok())
        .unwrap_or_else(|| panic!("not a number: {field}"))
}

/// bedpe: query name, start, end, target name, start, end, name, score,
/// query strand, target strand, `gi:f`, `bi:f`.
#[test]
fn test_bedpe_on_the_alignment_backend() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    write_fixture(work_dir);
    write_index(&binary, work_dir);

    let rows = query_rows(&binary, work_dir, "t2:200-800", "bedpe");
    assert_eq!(rows.len(), 1, "expected one bedpe row, got {rows:?}");
    let row = &rows[0];
    assert!(
        row.len() >= 12,
        "bedpe row is shorter than 12 columns: {row:?}"
    );
    assert_eq!(row[0], "q2", "query name");
    assert_eq!((row[1].as_str(), row[2].as_str()), ("200", "800"), "query interval");
    assert_eq!(row[3], "t2", "target name");
    assert_eq!((row[4].as_str(), row[5].as_str()), ("200", "800"), "target interval");
    assert_eq!(row[6], "t2:200-800", "row name is the queried region");
    assert_eq!(row[7], "0", "score");
    assert_eq!((row[8].as_str(), row[9].as_str()), ("+", "+"), "strands");
    assert!((number(&row[10]) - 1.0).abs() < 1e-9, "gi:f is 1 on a match");
    assert!((number(&row[11]) - 1.0).abs() < 1e-9, "bi:f is 1 on a match");

    // A deletion in the target shortens the projected query range by 10 bp and
    // stays a `D` in the CIGAR of the projected row.
    let rows = query_rows(&binary, work_dir, "B:0-500", "bedpe");
    assert_eq!(rows.len(), 1, "expected one bedpe row for B:0-500, got {rows:?}");
    let row = &rows[0];
    assert_eq!(row[0], "A");
    assert_eq!(
        (row[1].as_str(), row[2].as_str()),
        ("0", "490"),
        "10 deleted target bases must not move the query end"
    );
    assert_eq!((row[4].as_str(), row[5].as_str()), ("0", "500"));

    // A reverse-strand alignment queried from the query side reports both
    // strands in the last two strand columns.
    let rows = query_rows(&binary, work_dir, "bot:100-600", "bedpe");
    assert_eq!(rows.len(), 1, "expected one bedpe row for bot:100-600, got {rows:?}");
    let row = &rows[0];
    assert_eq!(row[0], "top");
    assert_eq!((row[1].as_str(), row[2].as_str()), ("2000", "2500"));
    assert_eq!((row[8].as_str(), row[9].as_str()), ("-", "+"));
}

/// paf: the 12 PAF columns plus `gi:f`, `bi:f`, `cg:Z` (the projected CIGAR
/// slice) and `an:Z` (the queried region name).
#[test]
fn test_paf_on_the_alignment_backend() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    write_fixture(work_dir);
    write_index(&binary, work_dir);

    let rows = query_rows(&binary, work_dir, "t2:0-500", "paf");
    assert_eq!(rows.len(), 1, "expected one paf row, got {rows:?}");
    let row = &rows[0];
    assert!(row.len() >= 16, "paf row is shorter than 16 columns: {row:?}");

    assert_eq!(row[0], "q2", "query name");
    assert_eq!(row[1], "1000", "query length");
    assert_eq!((row[2].as_str(), row[3].as_str()), ("0", "500"), "query interval");
    assert_eq!(row[4], "+", "strand");
    assert_eq!(row[5], "t2", "target name");
    assert_eq!(row[6], "1000", "target length");
    assert_eq!((row[7].as_str(), row[8].as_str()), ("0", "500"), "target interval");
    assert_eq!(
        (row[9].as_str(), row[10].as_str()),
        ("500", "500"),
        "matching bases and alignment block length"
    );

    let tags: std::collections::HashMap<&str, &str> = row[12..]
        .iter()
        .map(|t| t.split_once(':').map(|(k, _)| k).map(|k| (k, t.as_str())).unwrap())
        .collect();
    assert_eq!(
        tags.get("cg").map(|v| v.rsplit(':').next().unwrap()),
        Some("500="),
        "the projected CIGAR"
    );
    assert_eq!(
        tags.get("an").map(|v| v.rsplit_once("an:Z:").unwrap().1),
        Some("t2:0-500"),
        "the region name"
    );
    assert!((number(tags["gi"]) - 1.0).abs() < 1e-9, "gi:f is 1");
    assert!((number(tags["bi"]) - 1.0).abs() < 1e-9, "bi:f is 1");

    // The CIGAR slice keeps the deletion and the query range shrinks by it.
    let rows = query_rows(&binary, work_dir, "B:0-500", "paf");
    assert_eq!(rows.len(), 1, "expected one paf row for B:0-500, got {rows:?}");
    let row = &rows[0];
    assert_eq!((row[2].as_str(), row[3].as_str()), ("0", "490"));
    assert_eq!((row[7].as_str(), row[8].as_str()), ("0", "500"));
    let cigar = row
        .iter()
        .find(|f| f.starts_with("cg:Z:"))
        .expect("paf row has no cg:Z tag");
    assert_eq!(
        cigar, "cg:Z:400=10D90=",
        "the projected CIGAR must keep the deletion"
    );
}

#[test]
fn test_query_outputs_reject_a_range_below_min_transitive_len() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    write_fixture(work_dir);
    write_index(&binary, work_dir);

    // 100 bp is below the default --min-transitive-len of 101. The run must
    // fail loudly instead of printing an empty result that looks like "no
    // homology here".
    for format in ["bedpe", "paf"] {
        let output = run(
            &binary,
            work_dir,
            &[
                "query", "-a", "a.paf", "-i", "a.impg", "-r", "t2:0-100", "-d", "0", "-o", format,
            ],
        );
        assert!(
            !output.status.success(),
            "a 100 bp query with -o {format} must fail, not print nothing"
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("below minimum"),
            "the error must explain the minimum length"
        );
    }
}

/// A reciprocal pair: the same alignment written in both directions, which the
/// released HPRC r2 PAF set contains 1,483 times.
///
/// The index is bidirectional already, so record `A -> B` and record `B -> A`
/// describe the same edge twice. Today both reach the output: bedpe and paf
/// print the row twice. Any script that counts rows, sums block lengths or
/// builds a BEDPE of projected intervals sees double.
///
/// The test is ignored because the fix needs a decision: drop the duplicate at
/// index time (costs a set of seen edges over 389 million records in the
/// release), or keep the index as it is and deduplicate at output time.
#[test]
#[ignore = "C8: reciprocal records produce duplicate bedpe and paf rows"]
fn test_reciprocal_records_do_not_duplicate_rows() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();

    let mut file = File::create(work_dir.join("r.paf")).unwrap();
    for (q, t) in [("A", "B"), ("B", "A")] {
        writeln!(
            file,
            "{q}\t1000\t0\t500\t+\t{t}\t1000\t0\t500\t500\t500\t60\tcg:Z:500="
        )
        .unwrap();
    }

    let output = run(binary.as_path(), work_dir, &["index", "-a", "r.paf", "-i", "r.impg"]);
    assert!(output.status.success(), "index failed");

    for format in ["bedpe", "paf"] {
        let output = run(
            binary.as_path(),
            work_dir,
            &[
                "query", "-a", "r.paf", "-i", "r.impg", "-r", "A:0-500", "-d", "0", "-o", format,
            ],
        );
        assert!(output.status.success(), "query -o {format} failed");
        let rows = String::from_utf8_lossy(&output.stdout)
            .lines()
            .filter(|l| !l.trim().is_empty())
            .count();
        assert_eq!(
            rows, 1,
            "-o {format} printed {rows} rows for one edge described twice"
        );
    }
}
