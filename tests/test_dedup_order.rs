//! Records with the same endpoints must give the same output whatever the
//! order of the lines and of the files.
//!
//! The index kept one entry per edge (query, target, endpoints, strand) and
//! ignored the CIGAR, so the first record in sort order won: two records with
//! the same endpoints and different CIGARs projected differently by line
//! order, and in per-file mode the file order chose the CIGAR of a reciprocal
//! pair. Now the forward entry of a reciprocal pair wins, and two records of
//! the same direction both stay.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

fn query(dir: &Path, args: &[&str]) -> String {
    let out = Command::new(PathBuf::from(env!("CARGO_BIN_EXE_impg")))
        .current_dir(dir)
        .args(["query", "-o", "paf", "--no-merge", "-t", "1"])
        .args(args)
        .output()
        .unwrap();
    assert!(out.status.success(), "{args:?}: {}", String::from_utf8_lossy(&out.stderr));
    let mut rows: Vec<String> = String::from_utf8_lossy(&out.stdout).lines().map(str::to_string).collect();
    rows.sort();
    rows.join("\n")
}

const AB: &str = "A\t1000\t0\t500\t+\tB\t1000\t0\t500\t500\t500\t60\tcg:Z:500=";
// B -> A with an indel pair: the inverse of AB on the same endpoints, but not the same CIGAR
const BA: &str = "B\t1000\t0\t500\t+\tA\t1000\t0\t500\t499\t501\t60\tcg:Z:200=1I1D299=";

#[test]
fn reciprocal_pair_keeps_the_forward_record_in_any_order() {
    let mut seen = Vec::new();
    for (order, lines) in [("ab", [AB, BA]), ("ba", [BA, AB])] {
        let dir = TempDir::new().unwrap();
        fs::write(dir.path().join("both.paf"), format!("{}\n{}\n", lines[0], lines[1])).unwrap();
        fs::write(dir.path().join("ab.paf"), format!("{AB}\n")).unwrap();
        fs::write(dir.path().join("ba.paf"), format!("{BA}\n")).unwrap();
        let files: [&str; 2] = if order == "ab" { ["ab.paf", "ba.paf"] } else { ["ba.paf", "ab.paf"] };
        for region in ["B:0-500", "A:0-500"] {
            let single = query(dir.path(), &["-a", "both.paf", "-r", region]);
            for mode in ["single", "per-file"] {
                let split = query(dir.path(), &["-a", files[0], files[1], "--index-mode", mode, "-r", region]);
                assert_eq!(split, single, "{order} {mode} {region}");
            }
            // On B the record aligned to B is AB, on A it is BA
            let want = if region.starts_with('B') { "cg:Z:500=" } else { "cg:Z:200=1I1D299=" };
            assert!(single.contains(want) && single.lines().count() == 1, "{order} {region}: {single}");
            seen.push(single);
        }
    }
    assert_eq!(seen[0..2], seen[2..4], "line order changed the output");
}

#[test]
fn two_records_of_one_direction_both_stay() {
    let one = "q\t1000\t0\t200\t+\tt\t1000\t0\t200\t200\t200\t60\tcg:Z:200=";
    let two = "q\t1000\t0\t200\t+\tt\t1000\t0\t200\t199\t201\t60\tcg:Z:100=1I1D99=";
    let mut outputs = Vec::new();
    for lines in [[one, two], [two, one]] {
        let dir = TempDir::new().unwrap();
        fs::write(dir.path().join("x.paf"), format!("{}\n{}\n", lines[0], lines[1])).unwrap();
        let rows = query(dir.path(), &["-a", "x.paf", "-r", "t:0-200"]);
        assert_eq!(rows.lines().filter(|l| l.starts_with("q\t")).count(), 2, "{rows}");
        outputs.push(rows);
    }
    assert_eq!(outputs[0], outputs[1], "line order changed the output");
}
