//! Small native production CLI pipeline; no truth inputs or real panel loads.
use impg::genome_inference as genome;
use impg::syng::{SyncmerParams, SyngIndex};
use serde_json::{json, Value};
use std::{fs, path::Path, process::Command};
fn dna(n: usize, mut seed: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed & 3) as usize]
        })
        .collect()
}
fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}
fn run(args: &[&str], success: bool) {
    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(args)
        .output()
        .unwrap();
    assert_eq!(
        output.status.success(),
        success,
        "{args:?}\n{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}
fn read(path: &Path) -> Value {
    genome::read_json(path).unwrap()
}
fn rows(path: &Path) -> Vec<Value> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}
#[test]
fn partition_compiler_quantitative_rethread_reconstruct_and_fail_closed_cli() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let panel_path = root.join("panel.syng");
    let a = dna(2400, 17);
    let mut b = dna(2400, 29);
    b[1125..1275].copy_from_slice(&a[450..600]);
    let r = impg::graph::reverse_complement(&a[..1200]);
    let sequences = vec![
        ("A#0#chr1".into(), a.clone()),
        ("B#0#chr1".into(), b),
        ("R#0#reverse".into(), r),
    ];
    let mut panel = SyngIndex::build(SyncmerParams::default(), sequences.clone().into_iter());
    panel.save(p(&panel_path)).unwrap();
    let sources = root.join("sources.fa");
    let mut text = String::new();
    for (name, dna) in &sequences {
        text += &format!(">{name}\n{}\n", String::from_utf8_lossy(dna));
    }
    fs::write(&sources, &text).unwrap();
    let i = |name, start, end| json!({"path":name,"start":start,"end":end,"strand":null});
    let groups = root.join("groups.json");
    genome::write_json(&groups,&json!({"version":1,"groups":[
        {"id":"g0","scaffold":null,"occurrences":[i("A#0#chr1",0,1200),i("A#0#chr1",0,1200),i("B#0#chr1",0,1200),i("R#0#reverse",0,1200)]},
        {"id":"g1","scaffold":null,"occurrences":[i("A#0#chr1",1200,2400),i("B#0#chr1",1200,2400)]}]})).unwrap();
    let catalog = root.join("catalog");
    run(
        &[
            "genome-infer",
            "build-catalog",
            "--panel",
            p(&panel_path),
            "--groups",
            p(&groups),
            "--out-dir",
            p(&catalog),
        ],
        true,
    );
    let legacy = catalog.join("catalog.json");
    let compile = |out: &Path, scope: Option<&str>, length: &str, success| {
        let mut args = vec![
            "genome-infer",
            "build-partition-observations",
            "--panel",
            p(&panel_path),
            "--catalog",
            p(&legacy),
            "--sources",
            p(&sources),
            "--read-lengths",
            length,
            "--out-dir",
            p(out),
        ];
        if let Some(g) = scope {
            args.extend(["--feature-group", g]);
        }
        run(&args, success);
    };
    let observations = root.join("observations");
    compile(&observations, None, "150,220", true);
    let metadata = read(&observations.join("metadata.json"));
    assert_eq!(metadata["stats"]["sources_scanned"], 3);
    let stats = &metadata["stats"];
    let execution = &stats["profile_execution"];
    assert_eq!(stats["core_context_fetches"], stats["core_chunks"]);
    assert!(execution["reused_context_profiles"].as_u64().unwrap() > 200);
    assert!(
        execution["fallback_context_fetches"].as_u64().unwrap()
            < stats["core_context_fetches"].as_u64().unwrap()
    );
    assert!(execution["max_profile_context_bp"].as_u64().unwrap() <= 440);
    assert_eq!(
        execution["reused_context_profiles"].as_u64().unwrap()
            + execution["fallback_context_fetches"].as_u64().unwrap()
            + execution["no_start_profiles"].as_u64().unwrap(),
        stats["physical_shard_records_before_merge"]
            .as_u64()
            .unwrap()
    );
    assert_eq!(metadata["full_feature_scope"], true);
    eprintln!("compiler_stats {}", metadata["stats"]);
    let definitions = rows(&observations.join("definitions.jsonl"));
    let hits = rows(&observations.join("incidences.jsonl"));
    let old = genome::load_catalog(&legacy).unwrap();
    assert_eq!(definitions.len(), old.features.len());
    let mut recovered = 0;
    let mut physical = std::collections::BTreeSet::new();
    for h in &hits {
        let f = h["feature"].as_u64().unwrap() as usize;
        let source = h["source"].as_u64().unwrap() as usize;
        let start = h["start"].as_u64().unwrap();
        let end = h["end"].as_u64().unwrap();
        assert!(
            physical.insert((f, source, start, end)),
            "raw views or duplicate BED created a copy"
        );
        if !old.features[f]
            .locations
            .iter()
            .any(|l| l.source == source && l.start == start && l.end == end)
        {
            recovered += 1;
        }
        for q in h["contributions"].as_array().unwrap() {
            assert_eq!(q[0], q[1]);
        }
    }
    assert!(
        recovered > 0,
        "native reverse-view recovery must be exercised"
    );
    assert!(hits.iter().any(|h| h["context_nonlocal"] == true));
    assert!(
        hits.iter().any(|h| h["group"].is_null()),
        "new crossing physical support must be explicit"
    );
    let again = root.join("again");
    compile(&again, None, "220,150", true);
    assert_eq!(
        fs::read(observations.join("incidences.jsonl")).unwrap(),
        fs::read(again.join("incidences.jsonl")).unwrap()
    );
    compile(&observations, None, "150", false); // explicit resume unavailable
    compile(&root.join("unknown"), Some("missing"), "150", false);
    assert_eq!(
        read(&root.join("unknown/manifest.json"))["status"],
        "failed"
    );
    compile(&root.join("length-cap"), None, "1048577", false);
    let partial = root.join("partial");
    compile(&partial, Some("g0"), "150", true);
    let partial_meta = read(&partial.join("metadata.json"));
    assert_eq!(partial_meta["stats"]["sources_scanned"], 3);
    assert_eq!(partial_meta["stats"]["partitions_scanned"], 2);
    assert_eq!(partial_meta["full_feature_scope"], false);
    assert_eq!(
        rows(&partial.join("definitions.jsonl")).len(),
        old.features
            .iter()
            .filter(|f| f.owning_group == Some(0))
            .count()
    );
    let reads = root.join("reads.fa");
    let mut text = String::new();
    for start in 0..=a.len() - 150 {
        text += &format!(
            ">transient\n{}\n",
            String::from_utf8_lossy(&a[start..start + 150])
        );
    }
    fs::write(&reads, text).unwrap();
    let sample = root.join("sample");
    run(
        &[
            "genome-infer",
            "build-sample",
            "--panel",
            p(&panel_path),
            "--reads",
            p(&reads),
            "--out-dir",
            p(&sample),
        ],
        true,
    );
    let sample_path = sample.join("sample.membwt");
    let axis = root.join("axis.json");
    genome::write_json(&axis,&json!({"version":1,"coordinate_system":"fixture","intervals":[
        {"component":"chr1","start":0,"end":1200,"group":"g0","reference_occurrence":0,"reference_strand":"+"},
        {"component":"chr1","start":1200,"end":2400,"group":"g1","reference_occurrence":4,"reference_strand":"+"}]})).unwrap();
    let genotype = |directory: &Path, out: &Path, with_axis, success| {
        let mut args = vec![
            "genome-infer",
            "genotype-partitions",
            "--panel",
            p(&panel_path),
            "--observations",
            p(directory),
            "--sample",
            p(&sample_path),
            "--haploid-depth",
            "150",
            "--background",
            "0.1",
            "--max-mean-deviance",
            "10",
            "--allow-unvalidated-catalog",
            "--out-dir",
            p(out),
        ];
        if with_axis {
            args.extend(["--axis", p(&axis)]);
        }
        run(&args, success);
    };
    let calls = root.join("calls");
    genotype(&observations, &calls, true, true);
    let call_rows = read(&calls.join("calls.json"));
    assert_eq!(call_rows["model"], genome::observations::MODEL);
    assert_eq!(call_rows["catalog_payload_checksum"], Value::Null);
    assert!(
        rows(&calls.join("factor-ledger.jsonl"))
            .iter()
            .any(|f| f["original_owner"] == 0
                && f["owner"].is_null()
                && f["exclusion_reasons"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|r| r == "boundary-crossing")),
        "new global boundary support must invalidate the original local owner"
    );
    let g0 = &call_rows["calls"][0];
    let bundles = g0["bundles"].as_array().unwrap();
    let score = |name| {
        bundles.iter().find(|b| b["identity"] == name).unwrap()["score"]
            .as_f64()
            .unwrap()
    };
    assert!(
        (score("A#0") - score("R#0")).abs() < 1e-8,
        "equivalent observer predictions should tie"
    );
    assert_eq!(
        call_rows["observations"]["orientation_references"],
        json!([0, 4])
    );
    assert!(
        call_rows["observations"]["orientation_support"]
            .as_array()
            .unwrap()
            .len()
            <= old.occurrences.len()
    );
    let threads = read(&calls.join("threads.json"));
    assert_eq!(threads["resolved_intervals"], 2);
    assert!(
        threads["intervals"]
            .as_array()
            .unwrap()
            .iter()
            .flat_map(|r| r["states"].as_array().unwrap())
            .all(|s| s["identity"] != "B#0"),
        "DP resurrected unsupported source"
    );
    let rethread = root.join("rethread");
    run(
        &[
            "genome-infer",
            "thread-partitions",
            "--panel",
            p(&panel_path),
            "--observations",
            p(&observations),
            "--calls",
            p(&calls.join("calls.json")),
            "--axis",
            p(&axis),
            "--out-dir",
            p(&rethread),
        ],
        true,
    );
    let mut replayed = read(&rethread.join("threads.json"));
    // Match reconstruction's existing tolerance only for descriptive summed
    // segment scores after JSON float parsing; every state/selection stays exact.
    for (a, b) in threads["segments"]
        .as_array()
        .unwrap()
        .iter()
        .zip(replayed["segments"].as_array_mut().unwrap())
    {
        assert!((a["score"].as_f64().unwrap() - b["score"].as_f64().unwrap()).abs() < 1e-8);
        b["score"] = a["score"].clone();
    }
    assert_eq!(threads, replayed);
    let derived = read(&rethread.join("calls.json"));
    let mut original_scores = call_rows.clone();
    let mut derived_scores = derived.clone();
    original_scores
        .as_object_mut()
        .unwrap()
        .remove("observations");
    derived_scores
        .as_object_mut()
        .unwrap()
        .remove("observations");
    assert_eq!(
        original_scores, derived_scores,
        "rethread must not alter scores, statuses, rankings or rate summaries"
    );
    let original_bytes = fs::read(calls.join("calls.json")).unwrap();
    assert_eq!(
        derived["observations"]["derived_from_calls"]["fnv1a64"],
        format!("{:016x}", genome::checksum(&original_bytes))
    );
    assert_eq!(
        derived["observations"]["derived_from_calls"]["bytes"],
        original_bytes.len()
    );
    // Frozen-call headers and invariant provenance must be checked before
    // deriving any outputs; incompatible inputs must never be silently repaired.
    for (label, pointer, replacement) in [
        ("version", "/version", json!(999)),
        ("source-files", "/observations/source_files", json!([])),
        ("legacy-input", "/observations/legacy_input", json!({})),
        ("count-policy", "/count_policy", json!("unknown")),
        (
            "checksum-algorithm",
            "/artifact_checksum_algorithm",
            json!("unknown"),
        ),
        ("sample-checksum", "/sample_payload_checksum", json!("bad")),
        ("catalog-accepted", "/catalog_accepted", json!(true)),
        ("experimental", "/experimental", json!(false)),
        ("ploidy", "/ploidy", json!(2)),
        (
            "feature-groups",
            "/observations/feature_groups",
            json!(["foreign"]),
        ),
        (
            "compiler",
            "/observations/compiler_identity",
            json!("unknown"),
        ),
        (
            "full-scope",
            "/observations/full_feature_scope",
            json!(false),
        ),
        ("read-lengths", "/read_lengths", json!({"151": 1})),
        (
            "source-interval",
            "/source_occurrences/0/interval/start",
            json!(1),
        ),
        (
            "source-identity",
            "/source_occurrences/0/identity",
            json!("foreign#0"),
        ),
        ("group", "/calls/0/group", json!("foreign")),
        (
            "bundle-identity",
            "/calls/0/bundles/0/identity",
            json!("foreign#0"),
        ),
    ] {
        let mut incompatible = call_rows.clone();
        let field = incompatible.pointer_mut(pointer).unwrap();
        assert_ne!(*field, replacement);
        *field = replacement;
        let input = root.join(format!("incompatible-{label}.json"));
        let output = root.join(format!("incompatible-{label}"));
        genome::write_json(&input, &incompatible).unwrap();
        run(
            &[
                "genome-infer",
                "thread-partitions",
                "--panel",
                p(&panel_path),
                "--observations",
                p(&observations),
                "--calls",
                p(&input),
                "--axis",
                p(&axis),
                "--out-dir",
                p(&output),
            ],
            false,
        );
        assert!(!output.join("calls.json").exists());
        assert!(!output.join("threads.json").exists());
        assert_eq!(read(&input), incompatible);
    }
    let reconstruction = root.join("reconstruction");
    run(
        &[
            "genome-infer",
            "reconstruct",
            "--calls",
            p(&rethread.join("calls.json")),
            "--threads",
            p(&rethread.join("threads.json")),
            "--panel-names",
            p(&root.join("panel.syng.names")),
            "--sources",
            p(&sources),
            "--out-dir",
            p(&reconstruction),
        ],
        true,
    );
    let fasta = fs::read_to_string(reconstruction.join("reconstruction.fa")).unwrap();
    let spelled: String = fasta.lines().filter(|l| !l.starts_with('>')).collect();
    assert_eq!(spelled.as_bytes(), a);
    // Late axis selection starts from calls with NO orientation work. A new
    // reference is derived directly from reusable incidences, not transitively.
    let unthreaded = root.join("unthreaded");
    genotype(&observations, &unthreaded, false, true);
    let uncomputed = read(&unthreaded.join("calls.json"));
    assert_eq!(
        uncomputed["observations"]["orientation_references"],
        json!([])
    );
    assert_eq!(uncomputed["observations"]["orientation_support"], json!([]));
    let late_axis = root.join("late-axis.json");
    genome::write_json(&late_axis,&json!({"version":1,"coordinate_system":"late-fixture","intervals":[
        {"component":"first","start":0,"end":1200,"group":"g0","reference_occurrence":3,"reference_strand":"-"},
        {"component":"second","start":0,"end":1200,"group":"g1","reference_occurrence":4,"reference_strand":"+"}]})).unwrap();
    let late = root.join("late");
    run(
        &[
            "genome-infer",
            "thread-partitions",
            "--panel",
            p(&panel_path),
            "--observations",
            p(&observations),
            "--calls",
            p(&unthreaded.join("calls.json")),
            "--axis",
            p(&late_axis),
            "--out-dir",
            p(&late),
        ],
        true,
    );
    let late_calls = read(&late.join("calls.json"));
    assert_eq!(
        late_calls["observations"]["orientation_references"],
        json!([3, 4])
    );
    let mut before = uncomputed.clone();
    let mut after = late_calls.clone();
    before.as_object_mut().unwrap().remove("observations");
    after.as_object_mut().unwrap().remove("observations");
    assert_eq!(before, after);
    let identity = genome::PanelIdentity::read(p(&panel_path)).unwrap();
    let (registry, _) = genome::observations::load(&observations, &identity).unwrap();
    assert!(genome::observations::thread(
        &registry,
        &genome::read_json(&unthreaded.join("calls.json")).unwrap(),
        genome::read_json(&late_axis).unwrap(),
        10.0
    )
    .is_err());
    let late_reconstruction = root.join("late-reconstruction");
    run(
        &[
            "genome-infer",
            "reconstruct",
            "--calls",
            p(&late.join("calls.json")),
            "--threads",
            p(&late.join("threads.json")),
            "--panel-names",
            p(&root.join("panel.syng.names")),
            "--sources",
            p(&sources),
            "--out-dir",
            p(&late_reconstruction),
        ],
        true,
    );
    // The changed axis has separate components and an equivalent tied fragment.
    // Existing reconstruction emits resolved blocks before equivalent fragments;
    // these are independent records, not a promised concatenation order.
    let mut late_spelled: Vec<_> =
        genome::sequence_evaluation::read_fasta(&late_reconstruction.join("reconstruction.fa"))
            .unwrap()
            .into_values()
            .collect();
    let mut expected = vec![a[..1200].to_vec(), a[1200..].to_vec()];
    late_spelled.sort();
    expected.sort();
    assert_eq!(late_spelled, expected);
    genotype(&partial, &root.join("partial-calls"), false, true);
    genotype(&partial, &root.join("partial-axis"), true, false);
    assert!(!root.join("partial-axis/calls.json").exists());
    // Artifact-specific fingerprints and status, not the legacy header, gate use.
    use std::io::Write;
    fs::OpenOptions::new()
        .append(true)
        .open(observations.join("incidences.jsonl"))
        .unwrap()
        .write_all(b" \n")
        .unwrap();
    genotype(&observations, &root.join("tampered"), false, false);
    genotype(
        &root.join("unknown"),
        &root.join("failed-input"),
        false,
        false,
    );
    let mut wrong_panel = genome::PanelIdentity::read(p(&panel_path)).unwrap();
    wrong_panel.sidecars[0].2 = "0000000000000000".into();
    assert!(genome::observations::load(&again, &wrong_panel).is_err());
    let mut wrong_model = read(&again.join("manifest.json"));
    wrong_model["model"] = json!("foreign-model");
    genome::write_json(&again.join("manifest.json"), &wrong_model).unwrap();
    let identity = genome::PanelIdentity::read(p(&panel_path)).unwrap();
    assert!(genome::observations::load(&again, &identity).is_err());
    let wrong_length = root.join("wrong-length");
    compile(&wrong_length, None, "220", true);
    genotype(
        &wrong_length,
        &root.join("wrong-length-calls"),
        false,
        false,
    );
    eprintln!(
        "native recovered physical incidences={recovered}, compiled features={}, incidences={}",
        definitions.len(),
        hits.len()
    );
}
