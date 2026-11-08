use std::process::{Command, Output};

/// Two near-identical 3 kb haplotypes.
fn write_pair(dir: &std::path::Path) -> std::path::PathBuf {
    let mut seq = Vec::new();
    let mut x: u64 = 7;
    for _ in 0..3000 {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        seq.push(b"ACGT"[(x >> 62) as usize]);
    }
    let mut other = seq.clone();
    other[1500] = if other[1500] == b'A' { b'C' } else { b'A' };
    let path = dir.join("pair.fa");
    std::fs::write(&path, [b">A#0#chr1\n", &seq[..], b"\n>B#0#chr1\n", &other[..], b"\n"].concat())
        .unwrap();
    path
}

fn align(dir: &std::path::Path, extra: &[&str]) -> (Output, usize) {
    let out = dir.join(format!("out{}", extra.join("_").replace(['-', ' '], "")));
    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(["align", "--format", "paf", "--sequence-files"])
        .arg(write_pair(dir))
        .arg("-o")
        .arg(&out)
        .args(["-t", "1"])
        .args(extra)
        .output()
        .expect("run impg align");
    let lines = std::fs::read_to_string(out.join("alignments.paf"))
        .map(|s| s.lines().count())
        .unwrap_or(0);
    (output, lines)
}

/// `--min-aln-identity 90` means 0.9, as in impg graph; it was read as 90.0.
#[test]
fn align_min_aln_identity_percent_is_a_fraction() {
    let dir = tempfile::tempdir().unwrap();
    let (out, lines) = align(dir.path(), &["--min-aln-identity", "90"]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(lines > 0, "90 must keep the 99.97% identity alignments");
}

/// ANI presets need input alignments; they were silently taken as 0.
#[test]
fn align_min_aln_identity_ani_preset_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let (out, _) = align(dir.path(), &["--min-aln-identity", "ani50-2"]);
    assert!(!out.status.success(), "ani50-2 must be rejected by impg align");
}
