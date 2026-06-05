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

fn report_section_lines<'a>(report: &'a str, section: &str) -> Vec<&'a str> {
    let marker = format!("#section\t{section}");
    let mut in_section = false;
    let mut lines = Vec::new();
    for line in report.lines() {
        if line.starts_with("#section\t") {
            if in_section {
                break;
            }
            in_section = line == marker;
            continue;
        }
        if in_section && !line.trim().is_empty() {
            lines.push(line);
        }
    }
    lines
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

#[test]
fn genotype_cos_gfa_debug_report_exposes_lengths_repeats_and_scores() {
    let Some(bin) = impg_binary() else {
        eprintln!("SKIP: impg binary not built");
        return;
    };
    let dir = tempfile::tempdir().unwrap();
    let graph = dir.path().join("repeat_locus.gfa");
    std::fs::write(
        &graph,
        "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tCC
S\t3\tGGGG
S\t4\tT
P\tshort\t1+,2+,4+\t*
P\trepeat\t1+,2+,2+,4+\t*
P\tlong\t1+,3+,4+\t*
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
1\t1
2\t2
4\t1
",
    )
    .unwrap();
    let report = dir.path().join("genotype.report.tsv");

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
            "1",
            "--top-n",
            "3",
            "--debug-report",
            report.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg genotype cos --graph --debug-report");
    assert!(
        output.status.success(),
        "impg genotype cos --graph --debug-report failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout
            .lines()
            .any(|line| line.starts_with("1\tcos\t1\t") && line.contains("\trepeat\t")),
        "{stdout}"
    );

    let report = std::fs::read_to_string(&report).unwrap();
    assert!(report.contains("#impg genotype cos report"), "{report}");
    assert!(
        report.contains("contribution_model\traw"),
        "report should record the graph contribution model:\n{report}"
    );

    let universe_rows = report_section_lines(&report, "graph_feature_universe");
    assert!(
        universe_rows
            .iter()
            .any(|line| *line == "2\t2\t2\t2\t2.000000000"),
        "feature universe should expose segment name, length, raw count, and weight:\n{report}"
    );
    assert!(
        universe_rows
            .iter()
            .any(|line| *line == "3\t3\t4\t0\t0.000000000"),
        "feature universe should expose zero-evidence long allele node length:\n{report}"
    );

    let candidate_rows = report_section_lines(&report, "candidates");
    assert!(
        candidate_rows
            .iter()
            .any(|line| line.contains("\trepeat:0-6\t0\t6\t+\t4\t")
                && line.contains("\t3\t4\t3\t1\t1\t2\t")),
        "candidate summary should expose repeated-node counts for repeat path:\n{report}"
    );

    let feature_rows = report_section_lines(&report, "candidate_features");
    assert!(
        feature_rows
            .iter()
            .any(|line| *line == "0\trepeat\t2\t2\t2\t2.000000000\t2\t2.000000000\t4.000000000"),
        "candidate feature rows should expose repeated segment count and dot contribution:\n{report}"
    );

    let result_rows = report_section_lines(&report, "result_scores");
    assert!(
        result_rows.iter().any(|line| line.contains("\trepeat\t")),
        "result scores should include the final assignment:\n{report}"
    );
}
