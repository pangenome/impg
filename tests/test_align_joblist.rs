use std::process::Command;

/// The joblist FastGA -f counts haplotypes, as the PAF path does; it
/// counted sequences (11 haplotypes x 2 contigs gave -f22).
#[test]
fn align_joblist_frequency_counts_haplotypes() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("h.fa");
    let mut x: u64 = 5;
    let mut text = String::new();
    for h in 0..11 {
        for c in 0..2 {
            text.push_str(&format!(">H{h}#1#ctg{c}\n"));
            for _ in 0..2000 {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                text.push(b"ACGT"[(x >> 62) as usize] as char);
            }
            text.push('\n');
        }
    }
    std::fs::write(&fasta, text).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(["align", "--format", "joblist", "--sequence-files"])
        .arg(&fasta)
        .arg("-o")
        .arg(dir.path().join("out"))
        .args(["-t", "1"])
        .output()
        .expect("run impg align");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let jobs = std::fs::read_to_string(dir.path().join("out/align_jobs.txt")).unwrap();
    assert!(jobs.contains(" -f11 "), "expected -f11 (11 haplotypes):\n{jobs}");
}
