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

const REQUIRED_METHODS: [&str; 11] = [
    "local_syng_raw",
    "local_syng_crush_auto",
    "local_syng_crush_poa",
    "local_syng_crush_poasta",
    "local_syng_crush_sweepga",
    "top_flubble_nonoverlap_sweepga",
    "chunk_window_smooth_or_crush",
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
        assert!(root.join(row["command_log_path"].as_str().unwrap()).exists());
        assert!(root.join(row["stdout_log_path"].as_str().unwrap()).exists());
        assert!(root.join(row["stderr_log_path"].as_str().unwrap()).exists());
        assert_eq!(row["render_status"].as_str(), Some("skipped"));
        assert!(!row["render_skip_reason"].as_str().unwrap().is_empty());

        let method = row["method_id"].as_str().unwrap();
        if method.ends_with("_control") {
            assert_eq!(row["command_status"].as_str(), Some("skipped"));
            assert_eq!(
                row["skipped_optional_tool_reason"].as_str(),
                Some("profile_excludes_optional")
            );
            assert_eq!(row["exact_path_preservation"].as_str(), Some("not_run"));
        } else {
            assert_eq!(row["command_status"].as_str(), Some("pass"));
            assert_eq!(row["exact_path_preservation"].as_str(), Some("pass"));
            assert_eq!(row["hard_path_corruption"].as_bool(), Some(false));
            assert!(root.join(row["output_gfa_path"].as_str().unwrap()).exists());
            assert!(root.join(row["normalized_gfa_path"].as_str().unwrap()).exists());
            assert!(root.join(row["metrics_json_path"].as_str().unwrap()).exists());
            assert!(row["expected_topology_assertion_id"].as_str().is_some());
            assert!(row["expected_topology_message"].as_str().is_some());
            assert!(row["graph_size_bytes"].as_u64().unwrap() > 0);
            assert!(row["path_count"].as_u64().unwrap() > 0);
        }
    }
}
