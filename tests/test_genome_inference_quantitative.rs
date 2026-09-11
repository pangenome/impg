//! Actual reads -> reusable sample/catalog -> quantitative calls -> source threads.
//! Native fixture tests must run serialized.
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
fn path(p: &Path) -> &str {
    p.to_str().unwrap()
}
fn write(p: &Path, v: Value) {
    fs::write(p, serde_json::to_vec_pretty(&v).unwrap()).unwrap();
}
fn read(p: &Path) -> Value {
    genome::read_json(p).unwrap()
}
fn run(args: &[&str], success: bool) {
    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(args)
        .env("RAYON_NUM_THREADS", "4")
        .output()
        .unwrap();
    assert_eq!(
        output.status.success(),
        success,
        "{args:?}\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
fn interval(name: &str, start: u64, end: u64, strand: &str) -> Value {
    json!({"path":name,"start":start,"end":end,"strand":strand})
}
#[test]
fn quantitative_cli_bundles_threads_truth_independence_and_failures() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let panel_path = root.join("panel.syng");
    let mut sequences = vec![
        ("A#0#chr1".to_string(), dna(9000, 11)),
        ("B#0#foreign1".to_string(), dna(9000, 13)),
        ("A#0#chr2".to_string(), dna(9000, 17)),
        ("B#0#foreign2".to_string(), dna(9000, 19)),
        ("A#0#chr3".to_string(), dna(6000, 23)),
        ("B#0#reverse-contig".to_string(), dna(9000, 29)),
        ("A#0#copy1".to_string(), dna(3000, 31)),
        ("A#0#copy2".to_string(), dna(3000, 37)),
    ];
    let identical = sequences[2].1[..3000].to_vec();
    sequences[3].1[..3000].copy_from_slice(&identical);
    let mut panel = SyngIndex::build(SyncmerParams::default(), sequences.clone().into_iter());
    panel.save(path(&panel_path)).unwrap();
    let names: Vec<_> = sequences.iter().map(|(n, _)| n.as_str()).collect();
    let mut groups = Vec::new();
    for (g, start) in [("left-A", 0), ("middle-B", 3000), ("right-B", 6000)] {
        groups.push(json!({"id":g,"scaffold":null,"occurrences":[interval(names[0],start,start+3000,"+"), interval(names[1],start,start+3000,"+")]}));
    }
    groups.push(json!({"id":"tie","scaffold":null,"occurrences":[interval(names[2],0,3000,"+"),interval(names[3],0,3000,"+")]}));
    groups.push(json!({"id":"no-reads","scaffold":null,"occurrences":[interval(names[2],3000,6000,"+"),interval(names[3],3000,6000,"+")]}));
    groups.push(json!({"id":"repeat","scaffold":null,"occurrences":[interval(names[2],6000,7500,"+"),interval(names[2],7500,9000,"+")]}));
    groups.push(json!({"id":"reverse-left","scaffold":null,"occurrences":[interval(names[4],0,3000,"+"),interval(names[5],6000,9000,"-")]}));
    groups.push(json!({"id":"reverse-right","scaffold":null,"occurrences":[interval(names[4],3000,6000,"+"),interval(names[5],3000,6500,"-")]}));
    groups.push(json!({"id":"accessory-copies","scaffold":null,"occurrences":[interval(names[6],0,3000,"+"),interval(names[6],0,3000,"+"),interval(names[6],500,2500,"+"),interval(names[7],0,3000,"+")]}));
    let group_path = root.join("groups.json");
    write(&group_path, json!({"version":1,"groups":groups}));
    let reads = root.join("reads.fa");
    let mut text = String::new();
    for (source, begin, end) in [
        (0, 0, 3000),
        (1, 3000, 9000),
        (2, 0, 3000),
        (2, 6000, 9000),
        (5, 3000, 9000),
        (6, 0, 3000),
        (7, 0, 3000),
    ] {
        for start in (begin..=end - 500).step_by(100) {
            for _ in 0..2 {
                text += &format!(
                    ">transient\n{}\n",
                    String::from_utf8_lossy(&sequences[source].1[start..start + 500])
                );
            }
        }
    }
    fs::write(&reads, text).unwrap();
    let sample = root.join("sample");
    let catalog = root.join("catalog");
    run(
        &[
            "genome-infer",
            "build-sample",
            "--panel",
            path(&panel_path),
            "--reads",
            path(&reads),
            "--out-dir",
            path(&sample),
        ],
        true,
    );
    run(
        &[
            "genome-infer",
            "build-catalog",
            "--panel",
            path(&panel_path),
            "--groups",
            path(&group_path),
            "--out-dir",
            path(&catalog),
        ],
        true,
    );
    let cat = genome::load_catalog(&catalog.join("catalog.json")).unwrap();
    let reference = |group: &str, start: u64| {
        cat.occurrences
            .iter()
            .find(|o| {
                cat.groups[o.group].id == group
                    && o.interval.path.starts_with("A#0#")
                    && o.interval.start == start
            })
            .unwrap()
            .id
    };
    let axis_row = |component: &str, group: &str, start: u64, end: u64| json!({"component":component,"start":start,"end":end,"group":group,"reference_occurrence":reference(group,start),"reference_strand":"+"});
    let axis_value = json!({"version":1,"coordinate_system":"independent-fixture-reference-v1","intervals":[
        axis_row(names[0],"left-A",0,3000),axis_row(names[0],"middle-B",3000,6000),axis_row(names[0],"right-B",6000,9000),
        axis_row(names[2],"tie",0,3000),axis_row(names[2],"no-reads",3000,6000),axis_row(names[2],"repeat",6000,7500),axis_row(names[2],"repeat",7500,9000),
        axis_row(names[4],"reverse-left",0,3000),axis_row(names[4],"reverse-right",3000,6000)]});
    let axis = root.join("axis.json");
    write(&axis, axis_value.clone());
    let truth = root.join("truth.json");
    write(
        &truth,
        json!({"version":1,"coordinate_system":null,"groups":[
        {"group":"left-A","occurrences":[interval(names[0],0,3000,"+")]},
        {"group":"middle-B","occurrences":[interval(names[1],3000,6000,"+")]},
        {"group":"right-B","occurrences":[interval(names[1],6000,9000,"+")]},
        {"group":"accessory-copies","occurrences":[interval(names[6],0,3000,"+"),interval(names[7],0,3000,"+")]}]}),
    );
    let invoke = |out: &Path, truth_path: Option<&Path>, extra: &[&str], success| {
        let s = sample.join("sample.membwt");
        let c = catalog.join("catalog.json");
        let mut args = vec![
            "genome-infer",
            "genotype",
            "--panel",
            path(&panel_path),
            "--sample",
            path(&s),
            "--catalog",
            path(&c),
            "--haploid-depth",
            "10",
            "--axis",
            path(&axis),
            "--switch-penalty",
            "5",
            "--allow-unvalidated-catalog",
            "--out-dir",
            path(out),
        ];
        if let Some(t) = truth_path {
            args.extend(["--truth", path(t)]);
        }
        args.extend(extra);
        run(&args, success);
    };
    let out = root.join("quantitative");
    invoke(&out, Some(&truth), &[], true);
    let calls = read(&out.join("calls.json"));
    let threads = read(&out.join("threads.json"));
    assert_eq!(calls["model"], "haploid-bundle-composite-poisson-v1");
    assert_eq!(read(&out.join("manifest.json"))["model"], calls["model"]);
    let rows = calls["calls"].as_array().unwrap();
    let row = |g: &str| rows.iter().find(|r| r["group"] == g).unwrap();
    for (g, expected) in [
        ("left-A", "A#0"),
        ("middle-B", "B#0"),
        ("right-B", "B#0"),
        ("reverse-left", "B#0"),
        ("reverse-right", "B#0"),
    ] {
        assert_eq!(row(g)["status"], "informative", "{g}: {}", row(g));
        let b = row(g)["best_bundles"][0].as_u64().unwrap() as usize;
        assert_eq!(row(g)["bundles"][b]["identity"], expected);
    }
    assert_eq!(row("tie")["status"], "tied");
    assert_eq!(row("no-reads")["status"], "no-call-no-positive-features");
    assert_eq!(
        row("accessory-copies")["bundles"].as_array().unwrap().len(),
        1
    );
    assert_eq!(
        row("accessory-copies")["bundles"][0]["source_occurrences"]
            .as_array()
            .unwrap()
            .len(),
        4
    );
    assert_eq!(calls["feature_details_included"], false);
    assert!(calls["sample_payload_checksum"].is_string());
    assert!(calls["catalog_payload_checksum"].is_string());
    assert!(rows.iter().all(|r| r.get("factors").is_none()
        && r["bundles"]
            .as_array()
            .unwrap()
            .iter()
            .all(|b| b.get("physical_multiplicities").is_none())));
    let copy = &row("accessory-copies")["bundles"][0];
    assert_eq!(copy["physical_feature_locations"], copy["modeled_features"]);
    let debug = root.join("debug-details");
    invoke(&debug, None, &["--include-feature-details"], true);
    let mut detailed: genome::genotype::Genotypes =
        genome::read_json(&debug.join("calls.json")).unwrap();
    assert!(detailed.feature_details_included);
    assert!(detailed
        .calls
        .iter()
        .find(|c| c.group == "accessory-copies")
        .unwrap()
        .bundles[0]
        .physical_multiplicities
        .values()
        .all(|&m| m == 1));
    detailed.omit_feature_details();
    assert_eq!(serde_json::to_value(detailed).unwrap(), calls);
    assert_eq!(
        fs::read(debug.join("threads.json")).unwrap(),
        fs::read(out.join("threads.json")).unwrap()
    );
    let compact: genome::genotype::Genotypes = genome::read_json(&out.join("calls.json")).unwrap();
    assert!(
        !compact.feature_details_included && compact.calls.iter().all(|c| c.factors.is_empty())
    );
    assert_eq!(threads["repeated_axis_groups"], json!(["repeat"]));
    assert_eq!(threads["unscaffolded_groups"], json!(["accessory-copies"]));
    let tr = threads["intervals"].as_array().unwrap();
    assert_eq!(tr[3]["status"], "unresolved-tied-paths");
    assert_eq!(tr[4]["status"], "unresolved-no-call-no-positive-features");
    for r in &tr[5..7] {
        assert_eq!(r["status"], "unresolved-repeated-axis-group");
        assert!(r["states"].as_array().unwrap().is_empty());
        assert!(r["segment"].is_null());
    }
    assert_eq!(threads["segments"].as_array().unwrap().len(), 3);
    assert_eq!(threads["switches"].as_array().unwrap().len(), 1);
    assert_eq!(threads["switches"][0]["left_interval"], 0);
    assert_eq!(threads["switches"][0]["right_interval"], 1);
    assert_eq!(threads["switches"][0]["breakpoint_start"], 0);
    assert_eq!(threads["switches"][0]["breakpoint_end"], 6000);
    let blocks = threads["blocks"].as_array().unwrap();
    assert_eq!(blocks.len(), 3);
    assert_eq!(blocks[1]["interval_indices"], json!([1, 2]));
    assert_eq!(blocks[1]["source_gaps_bp"], json!([0]));
    assert_eq!(blocks[2]["strand"], "-");
    assert_eq!(blocks[2]["interval_indices"], json!([7, 8]));
    assert_eq!(blocks[2]["source_gaps_bp"], json!([-500]));
    assert_eq!(threads["resolved_intervals"], 5);
    assert_eq!(threads["axis_union_bp"], 24000);
    assert_eq!(threads["resolved_union_bp"], 15000);
    let eval = read(&out.join("evaluation.json"));
    assert_eq!(eval["threading"]["known_source_path_switches"], 1);
    assert_eq!(
        eval["threading"]["recovered_coarse_source_path_switches"],
        1
    );
    assert!(eval["bundle_compatibility"]
        .as_array()
        .unwrap()
        .iter()
        .all(|r| r["truth_set_contained_in_one_best_bundle"] == true));
    // Change truth, then remove truth: frozen quantitative calls and threads stay byte-identical.
    write(
        &truth,
        json!({"version":1,"coordinate_system":null,"groups":[{"group":"left-A","occurrences":[interval(names[1],0,3000,"+")]}]}),
    );
    for (name, t) in [("changed-truth", Some(truth.as_path())), ("no-truth", None)] {
        let next = root.join(name);
        invoke(&next, t, &[], true);
        for artifact in ["calls.json", "threads.json"] {
            assert_eq!(
                fs::read(out.join(artifact)).unwrap(),
                fs::read(next.join(artifact)).unwrap()
            );
        }
    }
    // Axis is not a candidate filter: omitting it leaves partition calls identical.
    let no_axis = root.join("no-axis");
    run(
        &[
            "genome-infer",
            "genotype",
            "--panel",
            path(&panel_path),
            "--sample",
            path(&sample.join("sample.membwt")),
            "--catalog",
            path(&catalog.join("catalog.json")),
            "--haploid-depth",
            "10",
            "--allow-unvalidated-catalog",
            "--out-dir",
            path(&no_axis),
        ],
        true,
    );
    assert_eq!(
        fs::read(out.join("calls.json")).unwrap(),
        fs::read(no_axis.join("calls.json")).unwrap()
    );
    assert!(!no_axis.join("threads.json").exists());
    // Error protocols after calls have been written must remove BOTH frozen result artifacts.
    let bad_truth = root.join("missing-truth");
    let failed = root.join("failed-truth");
    invoke(&failed, Some(&bad_truth), &[], false);
    assert_eq!(read(&failed.join("manifest.json"))["status"], "failed");
    assert!(!failed.join("calls.json").exists());
    assert!(!failed.join("threads.json").exists());
    for (name, extra) in [
        ("diploid", vec!["--ploidy", "2"]),
        ("bad-background", vec!["--background", "0"]),
        ("nan-background", vec!["--background", "NaN"]),
        ("bad-fit", vec!["--max-mean-deviance", "0"]),
    ] {
        let failed = root.join(name);
        invoke(&failed, None, &extra, false);
        assert_eq!(read(&failed.join("manifest.json"))["status"], "failed");
    }
    let poor = root.join("poor-fit");
    invoke(&poor, None, &["--max-mean-deviance", "0.0000001"], true);
    assert!(read(&poor.join("calls.json"))["calls"]
        .as_array()
        .unwrap()
        .iter()
        .any(|r| r["status"] == "poor-fit"));
    invoke(&out, None, &[], false); // Existing output directories are never reused.
    let mut malformed = axis_value.clone();
    malformed["intervals"][1]["component"] = names[2].into();
    write(&axis, malformed);
    let failed_axis = root.join("failed-axis");
    invoke(&failed_axis, None, &[], false);
    assert!(!failed_axis.join("calls.json").exists());
    // All emissions are used once: repeated-axis groups have no DP states; each other group occurs in one segment.
    let emitted: Vec<_> = tr
        .iter()
        .filter(|r| !r["segment"].is_null())
        .map(|r| r["axis"]["group"].as_str().unwrap())
        .collect();
    assert_eq!(
        emitted.len(),
        emitted
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
    );
    eprintln!("QUANTITATIVE FIXTURE {}", read(&out.join("manifest.json")));
}
