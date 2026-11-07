use std::io::Write;
use std::process::Command;

/// `impg align` on a plain-gzip (not BGZF) FASTA used to segfault.
#[test]
fn align_plain_gzip_fasta_does_not_crash() {
    let dir = tempfile::tempdir().expect("create test temp dir");
    let path = dir.path().join("pair.fa.gz");
    let mut gz = flate2::write::GzEncoder::new(
        std::fs::File::create(&path).expect("create FASTA"),
        flate2::Compression::default(),
    );
    let seq: String = (0..2000).map(|i| b"ACGT"[(i * 7 + i / 3) % 4] as char).collect();
    write!(gz, ">A#0#chr1\n{seq}\n>B#0#chr1\n{seq}\n").unwrap();
    gz.finish().expect("finish gzip FASTA");

    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(["align", "--aligner", "fastga", "--format", "paf", "--sequence-files"])
        .arg(&path)
        .arg("-o")
        .arg(dir.path().join("out"))
        .args(["-t", "1"])
        .output()
        .expect("run impg align");
    assert!(
        output.status.code().is_some(),
        "impg align was killed by a signal\nstderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output.status.success(),
        "impg align failed\nstderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
