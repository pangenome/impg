use std::process::Command;

fn query_gfa(dir: &std::path::Path, extra: &[&str]) -> (Option<i32>, String) {
    let seq: String = (0..500).map(|i| b"ACGT"[(i * 7 + i / 3) % 4] as char).collect();
    let fasta = dir.join("x.fa");
    std::fs::write(&fasta, format!(">a#1#c\n{seq}\n>b#1#c\n{seq}\n")).unwrap();
    let paf = dir.join("x.paf");
    std::fs::write(&paf, "b#1#c\t500\t0\t500\t+\ta#1#c\t500\t0\t500\t500\t500\t60\tcg:Z:500=\n").unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(["query", "-a"])
        .arg(&paf)
        .args(["-r", "a#1#c:0-500", "-d", "0", "-o", "gfa", "--sequence-files"])
        .arg(&fasta)
        .args(extra)
        .output()
        .expect("run impg query");
    (out.status.code(), String::from_utf8_lossy(&out.stderr).into_owned())
}

/// query -o gfa ran the aligner checks of impg graph only in graph;
/// now query, partition, graph and align share them.
#[test]
fn query_gfa_rejects_ignored_aligner_options() {
    let dir = tempfile::tempdir().unwrap();
    let (code, err) = query_gfa(dir.path(), &["--aligner", "wfmash", "--fastga-frequency", "5"]);
    assert_ne!(code, Some(0), "{err}");
    assert!(err.contains("--fastga-frequency is only supported with --aligner fastga"), "{err}");
    let (code, err) = query_gfa(dir.path(), &["--sparsify", "wfmash:auto"]);
    assert_ne!(code, Some(0), "{err}");
    assert!(err.contains("requires --aligner wfmash"), "{err}");
}
