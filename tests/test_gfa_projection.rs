use std::path::{Path, PathBuf};
use std::process::Command;

use impg::commands::genotype::{GraphContributionModel, GraphFeatureIdMode};
use impg::projection::converter::{
    project_gaf_to_gfa, GfaProjectionConfig, GfaProjectionOutputFormat,
};

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

fn tiny_graph() -> &'static str {
    "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
P\th1\t1+,2+,4+\t*
P\th2\t1+,3+,4+\t*
"
}

#[test]
fn gfa_gaf_projection_writes_expected_counts_metadata_and_repeated_visit_debug() {
    let dir = tempfile::tempdir().unwrap();
    let graph = dir.path().join("locus.gfa");
    let gaf = dir.path().join("reads.gaf");
    let projection = dir.path().join("sample.gfa.proj");
    std::fs::write(&graph, tiny_graph()).unwrap();
    std::fs::write(
        &gaf,
        "\
read_repeat\t3\t0\t3\t+\t>1>2>1\t3\t0\t3\t3\t3\t60
read_alt\t3\t0\t3\t+\t>1>3>4\t3\t0\t3\t3\t3\t60
",
    )
    .unwrap();

    let summary = project_gaf_to_gfa(&GfaProjectionConfig {
        gfa_path: &graph,
        gaf_path: &gaf,
        output_path: &projection,
        output_format: GfaProjectionOutputFormat::ProjectionBundle,
        feature_id_mode: GraphFeatureIdMode::SegmentName,
        contribution_model: GraphContributionModel::Raw,
        read_contributions_path: None,
    })
    .unwrap();

    assert_eq!(summary.feature_space, "gfa-segment");
    assert_eq!(summary.feature_id_mode, "segment-name");
    assert_eq!(summary.contribution_model, "raw");
    assert_eq!(summary.nonzero_features, 4);
    assert_eq!(summary.retained_records, 2);
    assert!(summary.graph_id.starts_with("gfa-fnv1a64:"));

    let pack = std::fs::read_to_string(projection.join("sample.pack.tsv")).unwrap();
    assert!(pack.contains("#feature_space\tgfa-segment"), "{pack}");
    assert!(
        pack.contains(&format!("#graph_id\t{}", summary.graph_id)),
        "{pack}"
    );
    assert!(pack.contains("#feature_id_mode\tsegment-name"), "{pack}");
    assert!(pack.contains("#graph_contribution_model\traw"), "{pack}");
    assert!(pack.contains("1\t3\n"), "{pack}");
    assert!(pack.contains("2\t1\n"), "{pack}");
    assert!(pack.contains("3\t1\n"), "{pack}");
    assert!(pack.contains("4\t1\n"), "{pack}");

    let debug = std::fs::read_to_string(projection.join("read-contributions.tsv")).unwrap();
    assert!(
        debug.contains(
            "read_repeat\t1\t3\t1\t+\t1\t2\t1\trepeated visit 2 to segment in read; counted again"
        ),
        "{debug}"
    );

    let manifest = std::fs::read_to_string(projection.join("manifest.json")).unwrap();
    assert!(
        manifest.contains("\"projection_method\": \"gaf-to-gfa\""),
        "{manifest}"
    );
    assert!(
        manifest.contains(&format!("\"graph_id\": \"{}\"", summary.graph_id)),
        "{manifest}"
    );
}

#[test]
fn gfa_gaf_projection_rejects_unknown_gaf_segments() {
    let dir = tempfile::tempdir().unwrap();
    let graph = dir.path().join("locus.gfa");
    let gaf = dir.path().join("reads.gaf");
    let pack = dir.path().join("sample.pack.tsv");
    std::fs::write(&graph, tiny_graph()).unwrap();
    std::fs::write(
        &gaf,
        "read_bad\t3\t0\t3\t+\t>1>missing>4\t3\t0\t3\t3\t3\t60\n",
    )
    .unwrap();

    let err = project_gaf_to_gfa(&GfaProjectionConfig {
        gfa_path: &graph,
        gaf_path: &gaf,
        output_path: &pack,
        output_format: GfaProjectionOutputFormat::PackTsv,
        feature_id_mode: GraphFeatureIdMode::SegmentName,
        contribution_model: GraphContributionModel::Raw,
        read_contributions_path: None,
    })
    .unwrap_err();

    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    assert!(
        err.to_string()
            .contains("GAF line 1 references unknown GFA segment 'missing'"),
        "{err}"
    );
}

#[test]
fn project_cli_bundle_can_feed_graph_genotype_without_pack_feature_space_override() {
    let Some(bin) = impg_binary() else {
        eprintln!("SKIP: impg binary not built");
        return;
    };
    let dir = tempfile::tempdir().unwrap();
    let graph = dir.path().join("locus.gfa");
    let gaf = dir.path().join("reads.gaf");
    let projection = dir.path().join("sample.gfa.proj");
    std::fs::write(&graph, tiny_graph()).unwrap();
    std::fs::write(
        &gaf,
        "\
read_h1\t3\t0\t3\t+\t>1>2>4\t3\t0\t3\t3\t3\t60
read_h2\t3\t0\t3\t+\t>1>3>4\t3\t0\t3\t3\t3\t60
",
    )
    .unwrap();

    let project_output = Command::new(&bin)
        .args([
            "project",
            "--gfa",
            graph.to_str().unwrap(),
            "--gaf",
            gaf.to_str().unwrap(),
            "-O",
            projection.to_str().unwrap(),
            "--graph-feature-id-mode",
            "segment-name",
        ])
        .output()
        .expect("failed to run impg project");
    assert!(
        project_output.status.success(),
        "impg project failed: {}",
        String::from_utf8_lossy(&project_output.stderr)
    );

    let genotype_output = Command::new(&bin)
        .args([
            "genotype",
            "cos",
            "--graph",
            graph.to_str().unwrap(),
            "--proj",
            projection.to_str().unwrap(),
            "--graph-feature-id-mode",
            "segment-name",
            "--ploidy",
            "2",
            "--top-n",
            "1",
        ])
        .output()
        .expect("failed to run impg genotype cos --graph --proj");
    assert!(
        genotype_output.status.success(),
        "impg genotype cos --graph --proj failed: {}",
        String::from_utf8_lossy(&genotype_output.stderr)
    );
    let stdout = String::from_utf8_lossy(&genotype_output.stdout);
    assert!(stdout.contains("#feature_space\tgfa-segment"), "{stdout}");
    assert!(
        stdout
            .lines()
            .any(|line| line.starts_with("1\tcos\t2\t") && line.contains("\th1,h2\t")),
        "{stdout}"
    );
}
