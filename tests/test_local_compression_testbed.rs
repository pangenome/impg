use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const REQUIRED_CLASSES: [&str; 13] = [
    "snp_bubble",
    "short_indel",
    "insertion_50_500bp",
    "alu_like_insertion",
    "adjacent_bubbles_compress_together",
    "bubble_split_by_fake_repeat_anchor",
    "repeated_motif_microtangle",
    "duplicated_flank_requires_path_context",
    "tandem_copy_number_loop_cyclic",
    "dispersed_repeat_glue_break_or_ignore",
    "inversion_like",
    "nested_bubbles_top_level_right",
    "nested_bubbles_top_level_wrong",
];

const REQUIRED_METHODS: [&str; 12] = [
    "local_syng_raw",
    "local_syng_crush_auto",
    "local_syng_crush_poa",
    "local_syng_crush_poasta",
    "local_syng_crush_sweepga",
    "top_flubble_nonoverlap_sweepga",
    "chunk_window_smooth_or_crush",
    "chunk_window_sweepga_seqwish",
    "whole_region_sweepga_seqwish",
    "pggb_control",
    "smoothxg_control",
    "pggb_plus_smoothxg_control",
];

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read_json(path: impl AsRef<Path>) -> Value {
    serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap()
}

fn parse_fasta(path: impl AsRef<Path>) -> BTreeMap<String, String> {
    let mut records = BTreeMap::new();
    let mut current_name: Option<String> = None;
    for raw in fs::read_to_string(path).unwrap().lines() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(name) = line.strip_prefix('>') {
            current_name = Some(name.to_string());
            records.insert(name.to_string(), String::new());
        } else {
            let name = current_name.as_ref().expect("sequence before FASTA header");
            records.get_mut(name).unwrap().push_str(line);
        }
    }
    records
}

fn parse_expected_paths(path: impl AsRef<Path>) -> BTreeMap<String, String> {
    let content = fs::read_to_string(path).unwrap();
    let mut lines = content.lines().map(str::to_string);
    assert_eq!(
        lines.next().as_deref(),
        Some("path_name\texpected_spelling")
    );
    lines
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let (name, spelling) = line.split_once('\t').unwrap();
            (name.to_string(), spelling.to_string())
        })
        .collect()
}

#[test]
fn local_compression_manifest_has_exact_fixture_paths_and_spellings() {
    let root = repo_root();
    let manifest_path = root.join("tests/test_data/local_compression/manifest.json");
    let manifest = read_json(&manifest_path);

    assert_eq!(
        manifest["schema"].as_str(),
        Some("local_compression_fixture_manifest_v1")
    );
    let required_classes: Vec<_> = manifest["required_fixture_classes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_str().unwrap())
        .collect();
    assert_eq!(required_classes, REQUIRED_CLASSES);

    let fixtures = manifest["fixtures"].as_array().unwrap();
    assert_eq!(fixtures.len(), REQUIRED_CLASSES.len());

    let mut seen_classes = BTreeSet::new();
    for entry in fixtures {
        let fixture_id = entry["fixture_id"].as_str().unwrap();
        let fixture_class = entry["fixture_class"].as_str().unwrap();
        seen_classes.insert(fixture_class.to_string());

        assert_eq!(
            entry["metadata_path"].as_str().unwrap(),
            format!("tests/test_data/local_compression/{fixture_id}/metadata.json")
        );
        assert_eq!(
            entry["input_fasta_path"].as_str().unwrap(),
            format!("tests/test_data/local_compression/{fixture_id}/input.fa")
        );
        assert_eq!(
            entry["expected_paths_path"].as_str().unwrap(),
            format!("tests/test_data/local_compression/{fixture_id}/expected_paths.tsv")
        );
        assert_eq!(
            entry["notes_path"].as_str().unwrap(),
            format!("tests/test_data/local_compression/{fixture_id}/notes.md")
        );

        let metadata_path = root.join(entry["metadata_path"].as_str().unwrap());
        let input_fasta_path = root.join(entry["input_fasta_path"].as_str().unwrap());
        let expected_paths_path = root.join(entry["expected_paths_path"].as_str().unwrap());
        let metadata = read_json(metadata_path);
        let fasta = parse_fasta(input_fasta_path);
        let expected_tsv = parse_expected_paths(expected_paths_path);

        assert_eq!(metadata["fixture_id"].as_str(), Some(fixture_id));
        assert_eq!(metadata["fixture_class"].as_str(), Some(fixture_class));
        assert_eq!(metadata["tier"], entry["tier"]);

        let expected_spellings: BTreeMap<String, String> = metadata["expected_path_spellings"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(name, spelling)| (name.clone(), spelling.as_str().unwrap().to_string()))
            .collect();
        assert_eq!(expected_tsv, expected_spellings);
        assert_eq!(fasta, expected_spellings);

        let input_names: BTreeSet<String> = metadata["input_sequence_names"]
            .as_array()
            .unwrap()
            .iter()
            .map(|name| name.as_str().unwrap().to_string())
            .collect();
        assert_eq!(input_names, expected_spellings.keys().cloned().collect());

        for record in metadata["input_sequences"].as_array().unwrap() {
            let name = record["name"].as_str().unwrap();
            assert_eq!(record["file"].as_str(), Some("input.fa"));
            assert_eq!(
                record["expected_spelling"].as_str(),
                expected_spellings.get(name).map(String::as_str)
            );
        }
    }

    assert_eq!(
        seen_classes,
        REQUIRED_CLASSES
            .iter()
            .map(|value| value.to_string())
            .collect()
    );
}

#[test]
fn local_compression_fast_runner_emits_complete_scoreboard_rows() {
    let root = repo_root();
    let out_dir = root.join("target/local-compression-testbed/cargo-test-snp");
    let _ = fs::remove_dir_all(&out_dir);

    let output = Command::new("python3")
        .current_dir(&root)
        .args([
            "scripts/local_compression_testbed.py",
            "run",
            "--profile",
            "fast",
            "--manifest",
            "tests/test_data/local_compression/manifest.json",
            "--fixtures",
            "snp_bubble_3path",
            "--out-dir",
            "target/local-compression-testbed/cargo-test-snp",
        ])
        .output()
        .expect("run local compression testbed runner");

    assert!(
        output.status.success(),
        "runner failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let scoreboard_json = read_json(out_dir.join("scoreboard.json"));
    let rows = scoreboard_json.as_array().unwrap();
    assert_eq!(rows.len(), REQUIRED_METHODS.len());

    let method_ids: BTreeSet<_> = rows
        .iter()
        .map(|row| row["method_id"].as_str().unwrap().to_string())
        .collect();
    assert_eq!(
        method_ids,
        REQUIRED_METHODS
            .iter()
            .map(|value| value.to_string())
            .collect()
    );

    let scoreboard_tsv = fs::read_to_string(out_dir.join("scoreboard.tsv")).unwrap();
    let header = scoreboard_tsv.lines().next().unwrap();
    for required_column in [
        "expected_topology_status",
        "graph_size_bytes",
        "path_replay_compression_ratio",
        "command_line",
        "command_log_path",
        "output_gfa_path",
        "render_status",
        "skipped_optional_tool_reason",
    ] {
        assert!(
            header.split('\t').any(|column| column == required_column),
            "missing TSV column {required_column}"
        );
    }

    for row in rows {
        assert_eq!(row["fixture_id"].as_str(), Some("snp_bubble_3path"));
        assert_eq!(row["profile"].as_str(), Some("fast"));
        assert_eq!(
            row["input_manifest_path"].as_str(),
            Some("tests/test_data/local_compression/manifest.json")
        );
        assert!(root
            .join(row["command_log_path"].as_str().unwrap())
            .exists());
        assert!(root.join(row["stdout_log_path"].as_str().unwrap()).exists());
        assert!(root.join(row["stderr_log_path"].as_str().unwrap()).exists());
        assert!(!row["command_line"].as_str().unwrap().is_empty());
        assert_eq!(row["render_status"].as_str(), Some("skipped"));
        assert!(!row["render_skip_reason"].as_str().unwrap().is_empty());

        let method = row["method_id"].as_str().unwrap();
        if method.ends_with("_control") {
            let command_line = row["command_line"].as_str().unwrap();
            assert!(
                command_line.contains("impg graph"),
                "control {method} must run through the internal impg graph CLI: {command_line}"
            );
            assert!(
                command_line.contains("--gfa-engine pggb"),
                "control {method} must select an internal impg GFA engine: {command_line}"
            );
            assert!(
                !command_line.contains(" smoothxg "),
                "control {method} must not invoke standalone smoothxg: {command_line}"
            );
            assert_ne!(
                row["command_status"].as_str(),
                Some("skipped"),
                "CI control {method} should execute through impg, not skip on external pggb/smoothxg discovery"
            );
            match row["command_status"].as_str().unwrap() {
                "pass" | "path_corrupt" => {
                    assert_eq!(row["tool_available"].as_str(), Some("true"));
                    assert!(root.join(row["output_gfa_path"].as_str().unwrap()).exists());
                    assert!(root
                        .join(row["normalized_gfa_path"].as_str().unwrap())
                        .exists());
                    assert!(root
                        .join(row["metrics_json_path"].as_str().unwrap())
                        .exists());
                    assert_eq!(
                        row["exact_path_preservation"].as_str(),
                        Some("pass"),
                        "control {method} should preserve exact input FASTA path names"
                    );
                    assert_eq!(row["hard_path_corruption"].as_bool(), Some(false));
                    assert!(row["expected_topology_assertion_id"].as_str().is_some());
                    assert!(row["expected_topology_message"].as_str().is_some());
                }
                "error" => {
                    assert_eq!(row["tool_available"].as_str(), Some("true"));
                    assert!(row["expected_topology_message"].as_str().is_some());
                    assert!(row["exit_code"].as_i64().is_some());
                }
                status => panic!("unexpected control command status {status}"),
            }
        } else {
            assert_eq!(row["command_status"].as_str(), Some("pass"));
            assert_eq!(row["exact_path_preservation"].as_str(), Some("pass"));
            assert_eq!(row["hard_path_corruption"].as_bool(), Some(false));
            assert!(root.join(row["output_gfa_path"].as_str().unwrap()).exists());
            assert!(root
                .join(row["normalized_gfa_path"].as_str().unwrap())
                .exists());
            assert!(root
                .join(row["metrics_json_path"].as_str().unwrap())
                .exists());
            assert!(row["expected_topology_assertion_id"].as_str().is_some());
            assert!(row["expected_topology_message"].as_str().is_some());
            assert!(row["graph_size_bytes"].as_u64().unwrap() > 0);
            assert!(row["path_count"].as_u64().unwrap() > 0);
            assert!(row["path_replay_compression_ratio"].as_f64().unwrap() >= 1.0);
        }
    }
}

#[test]
fn local_compression_fast_runner_supports_filtered_control_methods() {
    let root = repo_root();
    let out_dir = root.join("target/local-compression-testbed/cargo-test-filtered-control");
    let _ = fs::remove_dir_all(&out_dir);

    let output = Command::new("python3")
        .current_dir(&root)
        .args([
            "scripts/local_compression_testbed.py",
            "run",
            "--profile",
            "fast",
            "--manifest",
            "tests/test_data/local_compression/manifest.json",
            "--fixtures",
            "snp_bubble_3path",
            "--methods",
            "smoothxg_control",
            "--out-dir",
            "target/local-compression-testbed/cargo-test-filtered-control",
        ])
        .output()
        .expect("run filtered local compression testbed runner");

    assert!(
        output.status.success(),
        "filtered runner failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let rows = read_json(out_dir.join("scoreboard.json"));
    let rows = rows.as_array().unwrap();
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(row["method_id"].as_str(), Some("smoothxg_control"));
    assert_eq!(row["command_status"].as_str(), Some("pass"));
    assert_eq!(row["exact_path_preservation"].as_str(), Some("pass"));

    let report = fs::read_to_string(out_dir.join("report.md")).unwrap();
    assert!(report.contains("`smoothxg_control`"));
    assert!(!report.contains("`local_syng_raw` |"));
}

#[test]
fn local_compression_chunk_window_exposes_nested_parent_overmerge() {
    let root = repo_root();
    let out_dir = root.join("target/local-compression-testbed/cargo-test-nested-overmerge");
    let _ = fs::remove_dir_all(&out_dir);

    let output = Command::new("python3")
        .current_dir(&root)
        .args([
            "scripts/local_compression_testbed.py",
            "run",
            "--profile",
            "fast",
            "--manifest",
            "tests/test_data/local_compression/manifest.json",
            "--fixtures",
            "nested_top_level_wrong",
            "--methods",
            "top_flubble_nonoverlap_sweepga,chunk_window_smooth_or_crush",
            "--out-dir",
            "target/local-compression-testbed/cargo-test-nested-overmerge",
        ])
        .output()
        .expect("run nested overmerge local compression testbed case");

    assert!(
        output.status.success(),
        "nested overmerge runner failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let rows = read_json(out_dir.join("scoreboard.json"));
    let rows = rows.as_array().unwrap();
    assert_eq!(rows.len(), 2);

    let by_method: BTreeMap<_, _> = rows
        .iter()
        .map(|row| (row["method_id"].as_str().unwrap(), row))
        .collect();

    let chunk = by_method["chunk_window_smooth_or_crush"];
    assert_eq!(chunk["exact_path_preservation"].as_str(), Some("pass"));
    assert_eq!(chunk["hard_path_corruption"].as_bool(), Some(false));
    assert_eq!(chunk["expected_topology_status"].as_str(), Some("pass"));
    assert_eq!(chunk["candidate_count"].as_u64(), Some(2));
    assert_eq!(chunk["bubble_count"].as_u64(), Some(2));
    assert_eq!(chunk["flubble_count"].as_u64(), Some(2));

    let top_flubble = by_method["top_flubble_nonoverlap_sweepga"];
    assert_eq!(
        top_flubble["exact_path_preservation"].as_str(),
        Some("pass")
    );
    assert_eq!(
        top_flubble["expected_topology_status"].as_str(),
        Some("fail")
    );
    assert!(top_flubble["expected_topology_message"]
        .as_str()
        .unwrap()
        .contains("bubble_count: observed 1 below min 2"));
}

#[test]
fn local_compression_path_replay_compression_ratio() {
    let root = repo_root();
    let out_dir = root.join("target/local-compression-testbed/cargo-test-replay-ratio");
    let _ = fs::remove_dir_all(&out_dir);

    let output = Command::new("python3")
        .current_dir(&root)
        .args([
            "scripts/local_compression_testbed.py",
            "run",
            "--profile",
            "fast",
            "--manifest",
            "tests/test_data/local_compression/manifest.json",
            "--fixtures",
            "nested_top_level_wrong",
            "--methods",
            "local_syng_raw,chunk_window_smooth_or_crush",
            "--out-dir",
            "target/local-compression-testbed/cargo-test-replay-ratio",
        ])
        .output()
        .expect("run replay compression ratio local compression testbed case");

    assert!(
        output.status.success(),
        "replay ratio runner failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let rows = read_json(out_dir.join("scoreboard.json"));
    let rows = rows.as_array().unwrap();
    assert_eq!(rows.len(), 2);

    let by_method: BTreeMap<_, _> = rows
        .iter()
        .map(|row| (row["method_id"].as_str().unwrap(), row))
        .collect();

    let raw = by_method["local_syng_raw"];
    let chunk = by_method["chunk_window_smooth_or_crush"];
    assert_eq!(raw["exact_path_preservation"].as_str(), Some("pass"));
    assert_eq!(chunk["exact_path_preservation"].as_str(), Some("pass"));

    let raw_ratio = raw["path_replay_compression_ratio"].as_f64().unwrap();
    let chunk_ratio = chunk["path_replay_compression_ratio"].as_f64().unwrap();
    assert!(
        (raw_ratio - 1.0).abs() < 0.000001,
        "raw path-copy graph should have replay ratio 1.0, got {raw_ratio}"
    );
    assert!(
        chunk_ratio > raw_ratio,
        "chunk graph should show more replay compression than raw graph: raw={raw_ratio} chunk={chunk_ratio}"
    );
    assert_eq!(chunk["total_segment_bp"].as_u64(), Some(37));
    assert!(
        (chunk_ratio - (128.0 / 37.0)).abs() < 0.000001,
        "unexpected chunk replay ratio {chunk_ratio}"
    );

    let scoreboard_tsv = fs::read_to_string(out_dir.join("scoreboard.tsv")).unwrap();
    let header = scoreboard_tsv.lines().next().unwrap();
    assert!(header
        .split('\t')
        .any(|column| column == "path_replay_compression_ratio"));
}

#[test]
fn local_compression_chunk_window_sweepga_seqwish_nested_top_level_wrong() {
    let root = repo_root();
    let out_dir = root.join("target/local-compression-testbed/cargo-test-chunk-window-sweepga");
    let _ = fs::remove_dir_all(&out_dir);

    let output = Command::new("python3")
        .current_dir(&root)
        .args([
            "scripts/local_compression_testbed.py",
            "run",
            "--profile",
            "fast",
            "--manifest",
            "tests/test_data/local_compression/manifest.json",
            "--fixtures",
            "nested_top_level_wrong",
            "--methods",
            "chunk_window_sweepga_seqwish,pggb_control,smoothxg_control,pggb_plus_smoothxg_control",
            "--out-dir",
            "target/local-compression-testbed/cargo-test-chunk-window-sweepga",
        ])
        .output()
        .expect("run chunk-window SweepGA/seqwish local compression testbed case");

    assert!(
        output.status.success(),
        "chunk-window SweepGA/seqwish runner failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let rows = read_json(out_dir.join("scoreboard.json"));
    let rows = rows.as_array().unwrap();
    assert_eq!(rows.len(), 4);

    let by_method: BTreeMap<_, _> = rows
        .iter()
        .map(|row| (row["method_id"].as_str().unwrap(), row))
        .collect();

    let chunk = by_method["chunk_window_sweepga_seqwish"];
    assert_eq!(chunk["exact_path_preservation"].as_str(), Some("pass"));
    assert_eq!(chunk["hard_path_corruption"].as_bool(), Some(false));
    assert_eq!(chunk["expected_topology_status"].as_str(), Some("pass"));
    assert_eq!(chunk["candidate_count"].as_u64(), Some(2));
    assert_eq!(chunk["bubble_count"].as_u64(), Some(2));
    assert_eq!(chunk["flubble_count"].as_u64(), Some(2));
    assert!(chunk["path_replay_compression_ratio"].as_f64().unwrap() > 1.0);

    for control_id in [
        "pggb_control",
        "smoothxg_control",
        "pggb_plus_smoothxg_control",
    ] {
        let control = by_method[control_id];
        assert_eq!(
            control["exact_path_preservation"].as_str(),
            Some("pass"),
            "{control_id} should preserve exact paths"
        );
        assert_eq!(
            control["hard_path_corruption"].as_bool(),
            Some(false),
            "{control_id} should not corrupt paths"
        );
        assert_eq!(
            control["expected_topology_status"].as_str(),
            Some("fail"),
            "{control_id} remains a visible diagnostic topology-fail control row"
        );
        assert!(
            control["path_replay_compression_ratio"].as_f64().unwrap() >= 1.0,
            "{control_id} should report replay compression diagnostics"
        );
    }
}

#[test]
fn local_compression_chunk_window_sweepga_seqwish_is_resolver_distinct() {
    let root = repo_root();
    let out_dir = root.join("target/local-compression-testbed/cargo-test-chunk-window-distinct");
    let _ = fs::remove_dir_all(&out_dir);

    let output = Command::new("python3")
        .current_dir(&root)
        .args([
            "scripts/local_compression_testbed.py",
            "run",
            "--profile",
            "fast",
            "--manifest",
            "tests/test_data/local_compression/manifest.json",
            "--fixtures",
            "nested_top_level_wrong",
            "--methods",
            "chunk_window_smooth_or_crush,chunk_window_sweepga_seqwish",
            "--out-dir",
            "target/local-compression-testbed/cargo-test-chunk-window-distinct",
        ])
        .output()
        .expect("run resolver-distinct chunk-window local compression testbed case");

    assert!(
        output.status.success(),
        "resolver-distinct chunk-window runner failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let rows = read_json(out_dir.join("scoreboard.json"));
    let rows = rows.as_array().unwrap();
    assert_eq!(rows.len(), 2);

    let by_method: BTreeMap<_, _> = rows
        .iter()
        .map(|row| (row["method_id"].as_str().unwrap(), row))
        .collect();
    let smooth = by_method["chunk_window_smooth_or_crush"];
    let sweepga = by_method["chunk_window_sweepga_seqwish"];

    for row in [smooth, sweepga] {
        assert_eq!(row["exact_path_preservation"].as_str(), Some("pass"));
        assert_eq!(row["hard_path_corruption"].as_bool(), Some(false));
        assert_eq!(row["expected_topology_status"].as_str(), Some("pass"));
        assert_eq!(row["candidate_count"].as_u64(), Some(2));
        assert_eq!(row["bubble_count"].as_u64(), Some(2));
        assert_eq!(row["flubble_count"].as_u64(), Some(2));
    }

    let smooth_gfa = fs::read_to_string(root.join(smooth["normalized_gfa_path"].as_str().unwrap())).unwrap();
    let sweepga_gfa = fs::read_to_string(root.join(sweepga["normalized_gfa_path"].as_str().unwrap())).unwrap();
    assert_ne!(
        smooth_gfa, sweepga_gfa,
        "SweepGA/seqwish chunk row should exercise a resolver-distinct graph construction"
    );
    assert!(
        sweepga["total_segment_bp"].as_u64().unwrap() < smooth["total_segment_bp"].as_u64().unwrap(),
        "resolver-distinct row should reuse sequence across windows: smooth={} sweepga={}",
        smooth["total_segment_bp"],
        sweepga["total_segment_bp"]
    );
    assert!(
        sweepga["path_replay_compression_ratio"].as_f64().unwrap()
            > smooth["path_replay_compression_ratio"].as_f64().unwrap(),
        "resolver-distinct row should improve replay compression diagnostics: smooth={} sweepga={}",
        smooth["path_replay_compression_ratio"],
        sweepga["path_replay_compression_ratio"]
    );
}
