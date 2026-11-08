use std::io::{BufRead, BufReader, Write};
use std::process::Command;

/// 3 x 8 kb of yeast chrV (offset 100 kb).
fn write_slices(path: &std::path::Path) {
    let file = std::fs::File::open("tests/test_data/yeast.chrV.fa.gz").unwrap();
    let reader = BufReader::new(flate2::read::MultiGzDecoder::new(file));
    let mut records: Vec<(String, String)> = Vec::new();
    for line in reader.lines() {
        let line = line.unwrap();
        if let Some(name) = line.strip_prefix('>') {
            records.push((name.split_whitespace().next().unwrap().to_string(), String::new()));
        } else if let Some((_, seq)) = records.last_mut() {
            seq.push_str(line.trim());
        }
    }
    let mut out = std::fs::File::create(path).unwrap();
    for (name, seq) in records.iter().take(3) {
        writeln!(out, ">{name}\n{}", &seq[100_000..108_000]).unwrap();
    }
}

fn graph(dir: &std::path::Path, fasta: &std::path::Path, engine: &str, threads: &str, tag: &str) -> Vec<u8> {
    let gfa = dir.join(format!("{engine}_{tag}.gfa"));
    let out = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(["graph", "--sequence-files"])
        .arg(fasta)
        .args(["--gfa-engine", engine, "-t", threads, "--deterministic", "--temp-dir"])
        .arg(dir)
        .arg("-g")
        .arg(&gfa)
        .output()
        .expect("run impg graph");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    std::fs::read(&gfa).unwrap()
}

/// With --deterministic the GFA is the same bytes across runs and thread
/// counts. Without it, L-line order (pggb) and node order (poa) changed
/// from run to run, even with one thread.
#[test]
fn graph_output_is_reproducible_with_deterministic() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("y.fa");
    write_slices(&fasta);
    for engine in ["pggb", "poa"] {
        let first = graph(dir.path(), &fasta, engine, "4", "a");
        assert!(graph(dir.path(), &fasta, engine, "4", "b") == first, "{engine}: -t 4 runs differ");
        assert!(graph(dir.path(), &fasta, engine, "1", "c") == first, "{engine}: -t 1 differs from -t 4");
    }
}
