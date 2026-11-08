use std::process::{Command, Output};

/// Two near-identical 3 kb haplotypes, optionally soft-masked.
fn write_pair(dir: &std::path::Path, lower: bool, suffix: &str) -> std::path::PathBuf {
    let mut seq = Vec::new();
    let mut x: u64 = 11;
    for _ in 0..3000 {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        seq.push(b"ACGT"[(x >> 62) as usize]);
    }
    let mut other = seq.clone();
    other[1500] = if other[1500] == b'A' { b'C' } else { b'A' };
    if lower {
        seq.make_ascii_lowercase();
        other.make_ascii_lowercase();
    }
    let path = dir.join("pair.fa");
    let a = format!(">A#0#chr1{suffix}\n");
    let b = format!("\n>B#0#chr1{suffix}\n");
    std::fs::write(&path, [a.as_bytes(), &seq[..], b.as_bytes(), &other[..], b"\n"].concat())
        .unwrap();
    path
}

fn graph(dir: &std::path::Path, fasta: &std::path::Path, engine: &str) -> (Output, String) {
    let gfa = dir.join("out.gfa");
    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(["graph", "--sequence-files"])
        .arg(fasta)
        .args(["--gfa-engine", engine, "-t", "1", "-v", "1", "--temp-dir"])
        .arg(dir)
        .arg("-g")
        .arg(&gfa)
        .output()
        .expect("run impg graph");
    (output, std::fs::read_to_string(&gfa).unwrap_or_default())
}

/// `impg graph` without a window ran only the engine and dropped the
/// pipeline stages; `seqwish:smooth` now smooths.
#[test]
fn graph_without_window_runs_pipeline_stages() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = write_pair(dir.path(), false, "");
    let (out, gfa) = graph(dir.path(), &fasta, "seqwish:smooth");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "{stderr}");
    assert!(stderr.contains("smooth: "), "no smooth stage log:\n{stderr}");
    assert!(gfa.lines().any(|l| l.starts_with("P\t")), "no paths");
}

/// `cut-n=` clips query intervals; on whole FASTA input it is an error,
/// not a silent no-op.
#[test]
fn graph_without_window_rejects_terminal_n_clip() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = write_pair(dir.path(), false, "");
    let (out, _) = graph(dir.path(), &fasta, "cut-n=10:pggb");
    assert!(!out.status.success(), "cut-n= must be rejected without a window");
}

/// Soft-masked input kept `:start-end` names one base off (`:1-3001`).
#[test]
fn graph_pggb_keeps_names_of_soft_masked_input() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = write_pair(dir.path(), true, ":1-3000");
    let (out, gfa) = graph(dir.path(), &fasta, "pggb");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let mut names: Vec<&str> = gfa
        .lines()
        .filter_map(|l| l.strip_prefix("P\t"))
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    names.sort();
    assert_eq!(names, ["A#0#chr1:1-3000", "B#0#chr1:1-3000"]);
}

/// AGC input failed late with "stream did not contain valid UTF-8".
#[test]
fn graph_rejects_agc_input_up_front() {
    let dir = tempfile::tempdir().unwrap();
    let (out, _) = graph(dir.path(), std::path::Path::new("tests/test_data/test.agc"), "pggb");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success());
    assert!(stderr.contains("FASTA only"), "unclear error:\n{stderr}");
}
