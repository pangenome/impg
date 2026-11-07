use std::process::Command;

/// `impg align` with a bare relative FASTA name (no directory part).
#[test]
fn align_fastga_accepts_bare_relative_fasta_name() {
    let dir = tempfile::tempdir().expect("create test temp dir");
    let mut seq = String::new();
    let mut x: u64 = 7;
    for _ in 0..3000 {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        seq.push(b"ACGT"[(x >> 62) as usize] as char);
    }
    let mut mutated = seq.clone().into_bytes();
    mutated[1500] = if mutated[1500] == b'A' { b'C' } else { b'A' };
    let mutated = String::from_utf8(mutated).unwrap();
    std::fs::write(
        dir.path().join("pair.fa"),
        format!(">A#0#chr1\n{seq}\n>B#0#chr1\n{mutated}\n"),
    )
    .expect("write FASTA");

    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .current_dir(dir.path())
        .args(["align", "--aligner", "fastga", "--format", "paf"])
        .args(["--sequence-files", "pair.fa", "-o", "out", "-t", "1"])
        .output()
        .expect("run impg align");
    assert!(
        output.status.success(),
        "impg align failed with a bare relative FASTA name\nstderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        dir.path().join("out/alignments.paf").exists(),
        "impg align wrote no PAF"
    );
}
