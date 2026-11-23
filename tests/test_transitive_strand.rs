//! Strand of transitive query results.
//!
//! BED reports each homolog's strand relative to the queried region, composed
//! across alignment hops. Pairwise outputs (BEDPE, PAF) keep the strand of the
//! single alignment each row comes from.
//!
//! Fixture: B is the reverse complement of A and C is the reverse complement of
//! B, so C equals A. Alignments: B->A on '-' and C->B on '-'. Querying A:100-400
//! must give B:600-900 '-' and C:100-400 '+'.

use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

const LEN: usize = 1000;

fn paf_row(query: &str, target: &str, strand: char) -> String {
    format!("{query}\t{LEN}\t0\t{LEN}\t{strand}\t{target}\t{LEN}\t0\t{LEN}\t{LEN}\t{LEN}\t60\tcg:Z:{LEN}=\n")
}

/// Write one PAF per alignment so the same fixture runs in per-file mode.
fn write_chain(dir: &Path, hops: &[(&str, &str, char)]) -> Vec<String> {
    hops.iter()
        .enumerate()
        .map(|(i, (query, target, strand))| {
            let name = format!("hop{i}.paf");
            fs::write(dir.join(&name), paf_row(query, target, *strand)).unwrap();
            name
        })
        .collect()
}

fn query(dir: &Path, files: &[String], mode: &str, dfs: bool, format: &str) -> String {
    let mut args: Vec<String> = vec!["query".into(), "-a".into()];
    args.extend(files.iter().cloned());
    for arg in [
        "--index-mode",
        mode,
        "-r",
        "A:100-400",
        "-x",
        "-m",
        "2",
        "--no-merge",
        "-o",
        format,
    ] {
        args.push(arg.into());
    }
    if dfs {
        args.push("--transitive-dfs".into());
    }
    let out = Command::new(env!("CARGO_BIN_EXE_impg"))
        .current_dir(dir)
        .args(&args)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "impg {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

/// (name, start, end, strand) of each BED row
fn bed_rows(bed: &str) -> Vec<(String, u32, u32, char)> {
    bed.lines()
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 6, "BED row has {} columns: {line}", f.len());
            (
                f[0].to_string(),
                f[1].parse().unwrap(),
                f[2].parse().unwrap(),
                f[5].chars().next().unwrap(),
            )
        })
        .collect()
}

fn strands_of(rows: &[(String, u32, u32, char)], name: &str) -> Vec<char> {
    rows.iter().filter(|r| r.0 == name).map(|r| r.3).collect()
}

const MODES: [(&str, bool); 4] = [
    ("single", false),
    ("single", true),
    ("per-file", false),
    ("per-file", true),
];

#[test]
fn two_reverse_hops_compose_to_forward() {
    let tmp = TempDir::new().unwrap();
    let files = write_chain(tmp.path(), &[("B", "A", '-'), ("C", "B", '-')]);
    for (mode, dfs) in MODES {
        let rows = bed_rows(&query(tmp.path(), &files, mode, dfs, "bed"));
        let label = format!("--index-mode {mode}, dfs {dfs}: {rows:?}");
        assert!(rows.contains(&("B".into(), 600, 900, '-')), "{label}");
        assert!(rows.contains(&("C".into(), 100, 400, '+')), "{label}");
        assert_eq!(strands_of(&rows, "C"), vec!['+'], "{label}");
        let a = strands_of(&rows, "A");
        assert!(!a.is_empty() && a.iter().all(|&s| s == '+'), "{label}");
    }
}

#[test]
fn one_reverse_hop_in_a_chain_stays_reverse() {
    let tmp = TempDir::new().unwrap();
    // B is a forward copy of A, C the reverse complement of B
    let files = write_chain(tmp.path(), &[("B", "A", '+'), ("C", "B", '-')]);
    for (mode, dfs) in MODES {
        let rows = bed_rows(&query(tmp.path(), &files, mode, dfs, "bed"));
        let label = format!("--index-mode {mode}, dfs {dfs}: {rows:?}");
        assert_eq!(strands_of(&rows, "B"), vec!['+'], "{label}");
        assert_eq!(strands_of(&rows, "C"), vec!['-'], "{label}");
        assert!(rows.contains(&("C".into(), 600, 900, '-')), "{label}");
    }
}

#[test]
fn pairwise_outputs_keep_the_strand_of_each_hop() {
    let tmp = TempDir::new().unwrap();
    let files = write_chain(tmp.path(), &[("B", "A", '-'), ("C", "B", '-')]);
    for (mode, dfs) in MODES {
        let label = format!("--index-mode {mode}, dfs {dfs}");
        // BEDPE: the C row pairs C with B through the C->B alignment, which is '-'
        let bedpe = query(tmp.path(), &files, mode, dfs, "bedpe");
        let c_rows: Vec<&str> = bedpe.lines().filter(|l| l.starts_with("C\t")).collect();
        assert_eq!(c_rows.len(), 1, "{label}: {bedpe}");
        let f: Vec<&str> = c_rows[0].split('\t').collect();
        assert_eq!(
            (f[3], f[8], f[9]),
            ("B", "-", "+"),
            "{label}: {}",
            c_rows[0]
        );
        // PAF: the strand column of the C row must agree with its CIGAR hop
        let paf = query(tmp.path(), &files, mode, dfs, "paf");
        let c_rows: Vec<&str> = paf.lines().filter(|l| l.starts_with("C\t")).collect();
        assert_eq!(c_rows.len(), 1, "{label}: {paf}");
        let f: Vec<&str> = c_rows[0].split('\t').collect();
        assert_eq!((f[4], f[5]), ("-", "B"), "{label}: {}", c_rows[0]);
    }
}
