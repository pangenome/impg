use std::process::Command;

/// Two near-identical haplotypes: -f2 dropped every shared seed.
#[test]
fn align_two_haplotypes_finds_alignments() {
    let dir = tempfile::tempdir().expect("create test temp dir");
    let mut seq = Vec::new();
    let mut x: u64 = 7;
    for _ in 0..3000 {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        seq.push(b"ACGT"[(x >> 62) as usize]);
    }
    let mut mutated = seq.clone();
    mutated[1500] = if mutated[1500] == b'A' { b'C' } else { b'A' };
    let fasta = dir.path().join("pair.fa");
    std::fs::write(
        &fasta,
        [b">A#0#chr1\n", &seq[..], b"\n>B#0#chr1\n", &mutated[..], b"\n"].concat(),
    )
    .expect("write FASTA");

    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(["align", "--aligner", "fastga", "--format", "paf", "--sequence-files"])
        .arg(&fasta)
        .arg("-o")
        .arg(dir.path().join("out"))
        .args(["-t", "1"])
        .output()
        .expect("run impg align");
    assert!(
        output.status.success(),
        "impg align failed\nstderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let paf = std::fs::read_to_string(dir.path().join("out/alignments.paf")).expect("read PAF");
    assert!(paf.lines().count() > 0, "expected at least one alignment");
}
