//! Records with the same endpoints must give the same output whatever the
//! order of the lines and of the files.
//!
//! Two records with the same endpoints can carry different CIGARs: an indel in
//! a homopolymer fits at several places with the same score. The index kept
//! one entry per edge and ignored the CIGAR, so the first record in input
//! order won and the projection moved with the line or file order. Now the
//! record with fewer edits wins, then the one whose first difference comes
//! first, read in the frame where the sequence with the smaller name is the
//! target, so both directions of a reciprocal pair use the same record.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

fn query(dir: &Path, args: &[&str]) -> String {
    let out = Command::new(PathBuf::from(env!("CARGO_BIN_EXE_impg")))
        .current_dir(dir)
        .args(["query", "-o", "paf", "--no-merge", "-t", "1", "--min-transitive-len", "1"])
        .args(args)
        .output()
        .unwrap();
    assert!(out.status.success(), "{args:?}: {}", String::from_utf8_lossy(&out.stderr));
    let mut rows: Vec<String> = String::from_utf8_lossy(&out.stdout).lines().map(str::to_string).collect();
    rows.sort();
    rows.join("\n")
}

/// Query `regions` for both line orders of `lines` in one file, and for both
/// file orders in single and per-file mode; every run must print the same rows.
fn same_in_every_order(lines: [&str; 2], regions: &[&str]) -> Vec<String> {
    let mut first: Option<Vec<String>> = None;
    for order in [[0, 1], [1, 0]] {
        let dir = TempDir::new().unwrap();
        fs::write(dir.path().join("both.paf"), format!("{}\n{}\n", lines[order[0]], lines[order[1]])).unwrap();
        fs::write(dir.path().join("x.paf"), format!("{}\n", lines[0])).unwrap();
        fs::write(dir.path().join("y.paf"), format!("{}\n", lines[1])).unwrap();
        let files = if order[0] == 0 { ["x.paf", "y.paf"] } else { ["y.paf", "x.paf"] };
        let mut outputs = Vec::new();
        for region in regions {
            let single = query(dir.path(), &["-a", "both.paf", "-r", region]);
            for mode in ["single", "per-file"] {
                let split = query(dir.path(), &["-a", files[0], files[1], "--index-mode", mode, "-r", region]);
                assert_eq!(split, single, "order {order:?}, {mode}, {region}");
            }
            outputs.push(single);
        }
        match &first {
            Some(want) => assert_eq!(&outputs, want, "line order changed the output"),
            None => first = Some(outputs),
        }
    }
    first.unwrap()
}

#[test]
fn reciprocal_pair_uses_the_record_with_fewer_edits_both_ways() {
    let ab = "A\t1000\t0\t500\t+\tB\t1000\t0\t500\t500\t500\t60\tcg:Z:500=";
    // B -> A with an indel pair at base 200: same endpoints, 2 edits
    let ba = "B\t1000\t0\t500\t+\tA\t1000\t0\t500\t499\t501\t60\tcg:Z:200=1I1D299=";
    let rows = same_in_every_order([ab, ba], &["B:0-500", "A:0-500"]);
    for (region, rows) in ["B:0-500", "A:0-500"].iter().zip(&rows) {
        assert!(rows.contains("cg:Z:500=") && rows.lines().count() == 1, "{region}: {rows}");
    }
}

#[test]
fn equal_edits_keep_the_earliest_difference_in_both_directions() {
    // A = ACGTTTTG, B = ACGTTTG: the extra T of A fits after ACG or after ACGTTT
    let ab = "A\t8\t0\t8\t+\tB\t7\t0\t7\t7\t8\t60\tcg:Z:3=1I4=";
    let ba = "B\t7\t0\t7\t+\tA\t8\t0\t8\t7\t8\t60\tcg:Z:6=1D1=";
    let rows = same_in_every_order([ab, ba], &["B:3-4", "A:4-5"]);
    // 3=1I4= puts B:3 on A:4, with the inserted A:3 at its left edge (6=1I1=
    // would give A:3-4, 1=), and A:4 back on B:3 (6=1D1= would give B:4-5)
    assert!(rows[0].starts_with("A\t8\t3\t5\t") && rows[0].contains("cg:Z:1I1=") && rows[0].lines().count() == 1, "B:3-4: {}", rows[0]);
    assert!(rows[1].starts_with("B\t7\t3\t4\t") && rows[1].contains("cg:Z:1=") && rows[1].lines().count() == 1, "A:4-5: {}", rows[1]);
}

#[test]
fn same_direction_records_keep_the_one_with_fewer_edits() {
    let one = "q\t1000\t0\t200\t+\tt\t1000\t0\t200\t200\t200\t60\tcg:Z:200=";
    let two = "q\t1000\t0\t200\t+\tt\t1000\t0\t200\t199\t201\t60\tcg:Z:100=1I1D99=";
    let rows = same_in_every_order([one, two], &["t:0-200"]);
    assert!(rows[0].contains("cg:Z:200=") && rows[0].lines().count() == 1, "{}", rows[0]);
}

#[test]
fn per_file_mode_ranks_entries_that_project_to_nothing_too() {
    // 4=1I1=1D4= wins (2 edits against 10) and t:5 is its deleted base: no
    // homolog. Per-file mode used to rank only the entries with a projection,
    // so it kept 10X and found one.
    let indel = "q\t100\t0\t10\t+\tt\t100\t0\t10\t8\t10\t60\tcg:Z:4=1I1=1D4=";
    let mismatch = "q\t100\t0\t10\t+\tt\t100\t0\t10\t0\t10\t60\tcg:Z:10X";
    let rows = same_in_every_order([indel, mismatch], &["t:5-6"]);
    assert!(!rows[0].lines().any(|l| l.starts_with("q\t")), "t:5-6 has no homolog: {}", rows[0]);
}
