use std::io::{BufRead, BufReader, Write};
use std::process::Command;

/// Input: the first 30 kb of 3 yeast chrV haplotypes.
fn write_yeast_slices(path: &std::path::Path) -> usize {
    let file = std::fs::File::open("tests/test_data/yeast.chrV.fa.gz").expect("open yeast FASTA");
    let reader = BufReader::new(flate2::read::MultiGzDecoder::new(file));
    let mut records: Vec<(String, String)> = Vec::new();
    for line in reader.lines() {
        let line = line.expect("read yeast FASTA");
        if let Some(name) = line.strip_prefix('>') {
            records.push((name.split_whitespace().next().unwrap().to_string(), String::new()));
        } else if let Some((_, seq)) = records.last_mut() {
            seq.push_str(line.trim());
        }
    }
    let mut out = std::fs::File::create(path).expect("create FASTA");
    let mut total = 0;
    for (name, seq) in records.iter().take(3) {
        let slice = &seq[..30_000];
        writeln!(out, ">{name}\n{slice}").unwrap();
        total += slice.len();
    }
    total
}

/// Partitioned `impg graph` must keep every input base in the paths and
/// collapse homologous sequence: the old window loop duplicated intervals
/// (0 bp shared) and windowed only the first sequence (bases lost).
#[test]
fn graph_partitioned_pggb_keeps_paths_and_shares_sequence() {
    let dir = tempfile::tempdir().expect("create test temp dir");
    let fasta = dir.path().join("y.fa");
    let input_bp = write_yeast_slices(&fasta);
    let gfa_path = dir.path().join("out.gfa");
    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(["graph", "--sequence-files"])
        .arg(&fasta)
        .args(["--gfa-engine", "pggb:10000", "-d", "0", "-t", "2", "--temp-dir"])
        .arg(dir.path())
        .arg("-g")
        .arg(&gfa_path)
        .output()
        .expect("run impg graph");
    assert!(
        output.status.success(),
        "impg graph failed\nstderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let gfa = std::fs::read_to_string(&gfa_path).expect("read GFA");
    let mut seg_len = std::collections::HashMap::new();
    for line in gfa.lines().filter(|l| l.starts_with("S\t")) {
        let f: Vec<&str> = line.split('\t').collect();
        seg_len.insert(f[1].to_string(), f[2].len());
    }
    let graph_bp: usize = seg_len.values().sum();
    let path_bp: usize = gfa
        .lines()
        .filter(|l| l.starts_with("P\t"))
        .flat_map(|l| l.split('\t').nth(2).unwrap().split(',').collect::<Vec<_>>())
        .map(|step| seg_len[&step[..step.len() - 1]])
        .sum();
    assert_eq!(path_bp, input_bp, "paths must spell every input base");
    assert!(graph_bp < path_bp, "no shared sequence: graph {graph_bp} bp, paths {path_bp} bp");
}
