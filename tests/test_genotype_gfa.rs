use std::path::{Path, PathBuf};
use std::process::Command;

fn impg_binary() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return Some(PathBuf::from(path));
    }

    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    for candidate in [
        manifest_dir.join("target/debug/impg"),
        manifest_dir.join("target/release/impg"),
    ] {
        if candidate.exists() {
            return Some(candidate);
        }
    }
    None
}

#[test]
fn genotype_cos_prebuilt_gfa_pack_calls_expected_heterozygote() {
    let Some(bin) = impg_binary() else {
        eprintln!("SKIP: impg binary not built");
        return;
    };
    let dir = tempfile::tempdir().unwrap();
    let graph = dir.path().join("locus.gfa");
    std::fs::write(
        &graph,
        "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
P\th1\t1+,2+,4+\t*
P\th2\t1+,3+,4+\t*
",
    )
    .unwrap();
    let pack = dir.path().join("sample.pack.tsv");
    std::fs::write(
        &pack,
        "\
#feature_space\tgfa-segment
#feature_id_mode\tsegment-name
#node_id\tcount
1\t2
2\t1
3\t1
4\t2
",
    )
    .unwrap();

    let output = Command::new(bin)
        .args([
            "genotype",
            "cos",
            "--graph",
            graph.to_str().unwrap(),
            "--pack",
            pack.to_str().unwrap(),
            "--graph-feature-id-mode",
            "segment-name",
            "--ploidy",
            "2",
            "--top-n",
            "1",
        ])
        .output()
        .expect("failed to run impg genotype cos --graph");
    assert!(
        output.status.success(),
        "impg genotype cos --graph failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("#impg genotype cos"), "{stdout}");
    assert!(stdout.contains("#feature_space\tgfa-segment"), "{stdout}");
    assert!(
        stdout.contains("#graph_feature_id_mode\tsegment-name"),
        "{stdout}"
    );
    assert!(
        stdout
            .lines()
            .any(|line| line.starts_with("1\tcos\t2\t") && line.contains("\th1,h2\t")),
        "{stdout}"
    );
}
