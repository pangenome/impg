//! Real executable bootstrap fixtures. Run native tests with --test-threads=1.
use impg::genome_inference::{self as genome, catalog, sample::SampleIndex};
use impg::syng::{SyncmerParams, SyngIndex};
use serde_json::{json, Value};
use std::fs;
use std::path::Path;
use std::process::Command;

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
fn command(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(args)
        .env("RAYON_NUM_THREADS", "4")
        .output()
        .unwrap()
}
fn ok(args: &[&str]) {
    let result = command(args);
    assert!(
        result.status.success(),
        "{:?}\n{}",
        args,
        String::from_utf8_lossy(&result.stderr)
    );
}
fn path(p: &Path) -> &str {
    p.to_str().unwrap()
}
fn json_file(p: &Path) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn write_json(p: &Path, value: Value) {
    fs::write(p, serde_json::to_vec_pretty(&value).unwrap()).unwrap();
}
fn interval(name: &str, start: u64, end: u64, strand: Value) -> Value {
    json!({"path":name,"start":start,"end":end,"strand":strand})
}
fn group(id: &str, component: &str, occurrences: Vec<Value>) -> Value {
    json!({"id":id,"scaffold":{"component":component,"start":0,"end":3000},"occurrences":occurrences})
}

#[test]
fn genome_infer_actual_cli_multichromosome_counts_catalog_calls_and_failures() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let panel_path = root.join("panel.syng");
    let specs = [
        ("A#0#chr1", 11),
        ("B#0#chr1", 13),
        ("A#0#chr2", 17),
        ("B#0#chr2", 17),
        ("A#0#chr3", 19),
        ("B#0#chr3", 23),
        ("A#0#copy1", 29),
        ("A#0#copy2", 31),
        ("C#0#shared1", 37),
        ("C#0#shared2", 37),
        ("D#0#outside", 37),
        ("E#0#split", 41),
        ("F#0#whole", 41),
    ];
    let mut sequences: Vec<_> = specs
        .iter()
        .map(|(name, seed)| (name.to_string(), dna(3000, *seed)))
        .collect();
    sequences.push(("R#0#tandem".into(), dna(1000, 43).repeat(3)));
    let mut panel = SyngIndex::build(SyncmerParams::default(), sequences.clone().into_iter());
    panel.save(path(&panel_path)).unwrap();
    let split_walk = panel.walk_path_range(11, 0, 3000).unwrap();
    let cut = split_walk[split_walk.len() / 2].1 + 1;
    let groups = root.join("groups.json");
    write_json(
        &groups,
        json!({"version":1,"groups":[
            group("unique","chr1",vec![interval(specs[0].0,0,3000,json!("+")),interval(specs[1].0,0,3000,json!("-"))]),
            group("tie","chr2",vec![interval(specs[2].0,0,3000,json!("+")),interval(specs[3].0,0,3000,json!("+"))]),
            group("no-reads","chr3",vec![interval(specs[4].0,0,3000,json!("+")),interval(specs[5].0,0,3000,json!("+"))]),
            group("two-copies","unplaced-explicit",vec![interval(specs[6].0,0,3000,json!("+")),interval(specs[7].0,0,3000,json!("+"))]),
            group("shared1","repeat1",vec![interval(specs[8].0,0,3000,Value::Null)]),
            group("shared2","repeat2",vec![interval(specs[9].0,0,3000,Value::Null)]),
            {"id":"split-left","scaffold":null,"occurrences":[interval(specs[11].0,0,cut,json!("-"))]},
            {"id":"split-right","scaffold":null,"occurrences":[interval(specs[11].0,cut,3000,json!("-"))]},
            {"id":"whole","scaffold":null,"occurrences":[interval(specs[12].0,0,3000,json!("+"))]},
            {"id":"tandem","scaffold":null,"occurrences":[interval("R#0#tandem",0,3000,json!("+"))]}
        ]}),
    );
    // Deterministic error-free reads from all target chromosomes with signal,
    // including BOTH nonidentical copies of one haploid assembly identity.
    let mut records = Vec::new();
    for id in [0, 2, 6, 7, 8, 11] {
        for start in (0..=2500).step_by(100) {
            let seq = sequences[id].1[start..start + 500].to_vec();
            records.push(seq.clone());
            records.push(impg::graph::reverse_complement(&seq));
        }
    }
    let reads = root.join("reads.fastq");
    let write_reads = |target: &Path, rs: &[Vec<u8>]| {
        let mut text = String::new();
        for (i, r) in rs.iter().enumerate() {
            text += &format!(
                "@transient-read-{i}\n{}\n+\n{}\n",
                String::from_utf8_lossy(r),
                "I".repeat(r.len())
            );
        }
        fs::write(target, text).unwrap();
    };
    write_reads(&reads, &records);
    let truth = root.join("truth.json");
    write_json(
        &truth,
        json!({"version":1,"coordinate_system":"fixture-haploid-scaffold","groups":[
            {"group":"unique","occurrences":[interval(specs[0].0,0,3000,json!("+"))]},
            {"group":"tie","occurrences":[interval(specs[2].0,0,3000,json!("+"))]},
            {"group":"no-reads","occurrences":[interval(specs[4].0,0,3000,json!("+"))]},
            {"group":"two-copies","occurrences":[interval(specs[6].0,0,3000,json!("+")),interval(specs[7].0,0,3000,json!("+"))]}
        ]}),
    );
    let run = root.join("run");
    ok(&[
        "genome-infer",
        "run",
        "--panel",
        path(&panel_path),
        "--reads",
        path(&reads),
        "--groups",
        path(&groups),
        "--truth",
        path(&truth),
        "--out-dir",
        path(&run),
        "--allow-unvalidated-catalog",
    ]);
    let calls = json_file(&run.join("calls.json"));
    let rows = calls["calls"].as_array().unwrap();
    let row = |id: &str| rows.iter().find(|v| v["group"] == id).unwrap();
    assert_eq!(row("unique")["status"], "unique-compatible-occurrence");
    assert_eq!(row("unique")["compatible_source_occurrences"], json!([0]));
    assert_eq!(row("tie")["status"], "ambiguous-compatible-occurrences");
    assert_eq!(row("tie")["compatible_source_occurrences"], json!([2, 3]));
    assert_eq!(row("no-reads")["status"], "no-call-no-positive-features");
    assert_eq!(
        row("two-copies")["status"],
        "conflict-inadequate-locus-or-evidence-model"
    );
    assert_eq!(row("shared1")["status"], "no-call-no-local-features");
    assert!(row("two-copies")["locus_and_copy_structure_unresolved"]
        .as_bool()
        .unwrap());
    assert_eq!(calls["model"], "presence-compatibility");
    let cat = genome::load_catalog(&run.join("catalog.json")).unwrap();
    assert_eq!(cat.occurrences.len(), 14);
    let tandem = cat.occurrences.last().unwrap();
    assert!(tandem
        .feature_multiplicities
        .values()
        .any(|&count| count >= 2));
    assert!(cat
        .features
        .iter()
        .any(|f| f.locations.iter().filter(|l| l.source == 13).count() >= 2));
    assert!(!cat.catalog_accepted);
    assert_eq!(cat.occurrences[1].interval.path, "B#0#chr1");
    assert_eq!(cat.occurrences[1].interval.strand.as_deref(), Some("-"));
    assert!(cat
        .features
        .iter()
        .any(|f| f.exclusion_reasons.contains(&"boundary-crossing".into())));
    assert!(cat
        .features
        .iter()
        .any(|f| f.exclusion_reasons.contains(&"outside-catalog".into())
            && f.exclusion_reasons
                .contains(&"shared-between-groups".into())));
    assert!(cat
        .features
        .iter()
        .filter(|f| f.exclusion_reasons.contains(&"boundary-crossing".into()))
        .all(|f| f.owning_group.is_none()));
    assert_eq!(
        calls["global_factors_used"].as_u64().unwrap() as usize,
        cat.features
            .iter()
            .filter(|f| f.owning_group.is_some())
            .count()
    );
    let reverse_link = cat
        .links
        .iter()
        .find(|l| cat.occurrences[l.from].source == 11)
        .unwrap();
    assert_eq!(cat.occurrences[reverse_link.from].interval.start, cut);
    assert_eq!(cat.occurrences[reverse_link.to].interval.start, 0);
    assert_eq!(reverse_link.source_gap_bp, 0);
    assert!(reverse_link.oriented_continuation);
    assert!(cat
        .links
        .iter()
        .all(|l| cat.occurrences[l.from].source == cat.occurrences[l.to].source));
    let identity = genome::PanelIdentity::read(path(&panel_path)).unwrap();
    let sample = SampleIndex::load(&run.join("sample.membwt"), &identity).unwrap();
    assert_eq!(sample.stats.reads, records.len() as u64);
    assert!(sample.stats.mem_records > 0);
    assert!(calls["global_feature_counts"]
        .as_array()
        .unwrap()
        .iter()
        .any(|c| c.as_u64().unwrap() > 1));
    let evaluation = json_file(&run.join("evaluation.json"));
    assert_eq!(evaluation["evaluated_scaffold_union_bp"], 12000);
    assert_eq!(evaluation["truth_compatible_scaffold_union_bp"], 6000);
    assert_eq!(evaluation["genotype_callable_bp"], Value::Null);
    assert_eq!(json_file(&run.join("manifest.json"))["status"], "succeeded");
    eprintln!(
        "FIXTURE MEASUREMENT {}",
        json_file(&run.join("manifest.json"))
    );

    // Exactly one new positive feature is enough to conflict. Magnitude is
    // intentionally ignored, not modeled as coverage or an error probability.
    let mut observed: std::collections::BTreeMap<Vec<u64>, u64> = cat
        .features
        .iter()
        .enumerate()
        .filter(|(i, _)| calls["global_feature_counts"][i].as_u64().unwrap() > 0)
        .map(|(_, f)| (f.tokens.clone(), 1))
        .collect();
    let erroneous = cat.occurrences[1]
        .feature_multiplicities
        .keys()
        .find(|f| {
            cat.features[**f].owning_group == Some(0)
                && !cat.occurrences[0].feature_multiplicities.contains_key(f)
        })
        .unwrap();
    let mut injected = SampleIndex::load(&run.join("sample.membwt"), &identity).unwrap();
    injected.counts = impg::sample_mem_bwt::WeightedBwt::build(&observed).unwrap();
    let magnitude_ignored = genome::calling::call(&cat, &injected, true, 1).unwrap();
    assert_eq!(
        magnitude_ignored
            .calls
            .iter()
            .find(|c| c.group == "unique")
            .unwrap()
            .status,
        "unique-compatible-occurrence"
    );
    observed.insert(cat.features[*erroneous].tokens.clone(), 1);
    injected.counts = impg::sample_mem_bwt::WeightedBwt::build(&observed).unwrap();
    assert_eq!(
        genome::calling::call(&cat, &injected, true, 1)
            .unwrap()
            .calls
            .iter()
            .find(|c| c.group == "unique")
            .unwrap()
            .status,
        "conflict-inadequate-locus-or-evidence-model"
    );

    // Reusable stages produce byte-identical counts independent of read order.
    records.reverse();
    let reordered = root.join("reordered.fastq");
    write_reads(&reordered, &records);
    let sample2 = root.join("sample2");
    ok(&[
        "genome-infer",
        "build-sample",
        "--panel",
        path(&panel_path),
        "--reads",
        path(&reordered),
        "--out-dir",
        path(&sample2),
    ]);
    assert_eq!(
        fs::read(run.join("sample.membwt")).unwrap(),
        fs::read(sample2.join("sample.membwt")).unwrap()
    );
    // Switching every input read orientation preserves the canonical collection.
    let single_orientation: Vec<_> = records.iter().step_by(2).cloned().collect();
    let single_reads = root.join("single.fastq");
    write_reads(&single_reads, &single_orientation);
    let single_out = root.join("single");
    ok(&[
        "genome-infer",
        "build-sample",
        "--panel",
        path(&panel_path),
        "--reads",
        path(&single_reads),
        "--out-dir",
        path(&single_out),
    ]);
    let flipped: Vec<_> = single_orientation
        .iter()
        .map(|r| impg::graph::reverse_complement(r))
        .collect();
    let flipped_reads = root.join("flipped.fastq");
    write_reads(&flipped_reads, &flipped);
    let flipped_out = root.join("flipped");
    ok(&[
        "genome-infer",
        "build-sample",
        "--panel",
        path(&panel_path),
        "--reads",
        path(&flipped_reads),
        "--out-dir",
        path(&flipped_out),
    ]);
    assert_eq!(
        fs::read(single_out.join("sample.membwt")).unwrap(),
        fs::read(flipped_out.join("sample.membwt")).unwrap()
    );
    let catalog2 = root.join("catalog2");
    ok(&[
        "genome-infer",
        "build-catalog",
        "--panel",
        path(&panel_path),
        "--groups",
        path(&groups),
        "--out-dir",
        path(&catalog2),
    ]);
    assert_eq!(
        fs::read(run.join("catalog.json")).unwrap(),
        fs::read(catalog2.join("catalog.json")).unwrap()
    );
    let calls2 = root.join("calls2");
    ok(&[
        "genome-infer",
        "call",
        "--panel",
        path(&panel_path),
        "--sample",
        path(&sample2.join("sample.membwt")),
        "--catalog",
        path(&catalog2.join("catalog.json")),
        "--out-dir",
        path(&calls2),
        "--allow-unvalidated-catalog",
    ]);
    assert_eq!(
        fs::read(run.join("calls.json")).unwrap(),
        fs::read(calls2.join("calls.json")).unwrap()
    );

    // A changed downstream truth cannot change a single byte of frozen calls.
    let mut altered = json_file(&truth);
    altered["groups"][0]["occurrences"][0] = interval(specs[1].0, 0, 3000, json!("-"));
    write_json(&truth, altered);
    let calls3 = root.join("calls3");
    ok(&[
        "genome-infer",
        "call",
        "--panel",
        path(&panel_path),
        "--sample",
        path(&sample2.join("sample.membwt")),
        "--catalog",
        path(&catalog2.join("catalog.json")),
        "--truth",
        path(&truth),
        "--out-dir",
        path(&calls3),
        "--allow-unvalidated-catalog",
    ]);
    assert_eq!(
        fs::read(run.join("calls.json")).unwrap(),
        fs::read(calls3.join("calls.json")).unwrap()
    );
    assert_eq!(
        json_file(&calls3.join("evaluation.json"))["groups"][0]
            ["truth_set_contained_in_compatible_occurrences"],
        false
    );

    for (label, extra) in [
        ("ack", vec![]),
        (
            "ploidy",
            vec!["--allow-unvalidated-catalog", "--ploidy", "2"],
        ),
    ] {
        let out = root.join(label);
        let mut args = vec![
            "genome-infer",
            "call",
            "--panel",
            path(&panel_path),
            "--sample",
            path(&sample2.join("sample.membwt")),
            "--catalog",
            path(&catalog2.join("catalog.json")),
            "--out-dir",
            path(&out),
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
        args.extend(extra.into_iter().map(String::from));
        let result = command(&args.iter().map(String::as_str).collect::<Vec<_>>());
        assert!(!result.status.success());
        assert_eq!(json_file(&out.join("manifest.json"))["status"], "failed");
        assert!(!out.join("calls.json").exists());
    }
    // Re-reserving an output never changes its accepted artifacts.
    assert!(!command(&[
        "genome-infer",
        "build-sample",
        "--panel",
        path(&panel_path),
        "--reads",
        path(&reads),
        "--out-dir",
        path(&run)
    ])
    .status
    .success());
    assert_eq!(json_file(&run.join("manifest.json"))["status"], "succeeded");
    let corrupt = root.join("corrupt.membwt");
    let mut bytes = fs::read(run.join("sample.membwt")).unwrap();
    bytes[30] ^= 1;
    fs::write(&corrupt, &bytes).unwrap();
    assert!(SampleIndex::load(&corrupt, &identity).is_err());
    fs::write(&corrupt, &bytes[..20]).unwrap();
    assert!(SampleIndex::load(&corrupt, &identity).is_err());
    let mut wrong = identity.clone();
    wrong.sidecars[1].2 = "different-dictionary".into();
    assert!(SampleIndex::load(&run.join("sample.membwt"), &wrong).is_err());
    let other_prefix = root.join("different-panel");
    let mut other_panel = SyngIndex::build(
        SyncmerParams::default(),
        vec![("Z#0#chr1".into(), dna(3000, 999))].into_iter(),
    );
    other_panel.save(path(&other_prefix)).unwrap();
    let other_identity = genome::PanelIdentity::read(path(&other_prefix)).unwrap();
    assert_ne!(identity, other_identity);
    assert!(SampleIndex::load(&run.join("sample.membwt"), &other_identity).is_err());
    let corrupt_cat = root.join("corrupt.json");
    let mut bad = json_file(&run.join("catalog.json"));
    bad["payload"]["occurrences"][0]["interval"]["end"] = json!(4000);
    write_json(&corrupt_cat, bad);
    assert!(genome::load_catalog(&corrupt_cat).is_err());
    // Truth is opened downstream: malformed truth leaves FAILED, no accepted calls/eval.
    let bad_truth = root.join("bad-truth.json");
    fs::write(&bad_truth, "{not json").unwrap();
    let failed = root.join("bad-truth-run");
    assert!(!command(&[
        "genome-infer",
        "call",
        "--panel",
        path(&panel_path),
        "--sample",
        path(&run.join("sample.membwt")),
        "--catalog",
        path(&run.join("catalog.json")),
        "--truth",
        path(&bad_truth),
        "--out-dir",
        path(&failed),
        "--allow-unvalidated-catalog"
    ])
    .status
    .success());
    assert_eq!(json_file(&failed.join("manifest.json"))["status"], "failed");
    assert!(!failed.join("calls.json").exists());
    assert!(!failed.join("evaluation.json").exists());

    // Valid individual scaffold coordinates can overflow their cross-component union.
    // Exercise evaluation failure AFTER calls are written and require their cleanup.
    let mut overflow_catalog = genome::load_catalog(&run.join("catalog.json")).unwrap();
    let mut overflow_truth = json_file(&truth);
    let evaluated_groups = overflow_truth["groups"].as_array_mut().unwrap();
    assert!(evaluated_groups.len() >= 2);
    evaluated_groups.truncate(2);
    for (i, t) in evaluated_groups.iter().enumerate() {
        let group = overflow_catalog
            .groups
            .iter_mut()
            .find(|g| Some(g.id.as_str()) == t["group"].as_str())
            .unwrap();
        group.scaffold = Some(catalog::Scaffold {
            component: format!("overflow-component-{i}"),
            start: 0,
            end: if i == 0 { u64::MAX } else { 1 },
        });
    }
    let overflow_catalog_path = root.join("overflow-catalog.json");
    genome::save_catalog(&overflow_catalog_path, &overflow_catalog).unwrap();
    let overflow_truth_path = root.join("overflow-truth.json");
    write_json(&overflow_truth_path, overflow_truth);
    let overflow_out = root.join("overflow-evaluation");
    let overflow_result = command(&[
        "genome-infer",
        "call",
        "--panel",
        path(&panel_path),
        "--sample",
        path(&run.join("sample.membwt")),
        "--catalog",
        path(&overflow_catalog_path),
        "--truth",
        path(&overflow_truth_path),
        "--out-dir",
        path(&overflow_out),
        "--allow-unvalidated-catalog",
    ]);
    assert!(!overflow_result.status.success());
    let overflow_manifest = json_file(&overflow_out.join("manifest.json"));
    assert_eq!(overflow_manifest["status"], "failed");
    assert!(overflow_manifest["error"]
        .as_str()
        .unwrap()
        .contains("scaffold union base-pair count overflow"));
    assert!(!overflow_out.join("calls.json").exists());
    assert!(!overflow_out.join("evaluation.json").exists());

    // Reserving an existing EMPTY output also fails without populating it.
    let empty_out = root.join("existing-empty-output");
    fs::create_dir(&empty_out).unwrap();
    assert!(!command(&[
        "genome-infer",
        "build-sample",
        "--panel",
        path(&panel_path),
        "--reads",
        path(&reads),
        "--out-dir",
        path(&empty_out),
    ])
    .status
    .success());
    assert_eq!(fs::read_dir(&empty_out).unwrap().count(), 0);

    // BED3 import preserves duplicate physical occurrences and leaves orientation unknown.
    let beds = root.join("beds");
    fs::create_dir(&beds).unwrap();
    fs::write(
        beds.join("one.bed"),
        "A#0#chr1\t0\t3000\nA#0#chr1\t0\t3000\n",
    )
    .unwrap();
    let imported = catalog::build(
        &panel,
        identity.clone(),
        catalog::import_beds(&beds).unwrap(),
    )
    .unwrap();
    assert_eq!(imported.occurrences.len(), 2);
    assert!(imported
        .occurrences
        .iter()
        .all(|o| o.interval.strand.is_none()));
    assert!(imported
        .features
        .iter()
        .all(|f| f.locations[0].containing_occurrences == [0, 1]));
    fs::write(beds.join("one.bed"), "A#0#chr1\t0\t3001\n").unwrap();
    assert!(catalog::build(
        &panel,
        identity.clone(),
        catalog::import_beds(&beds).unwrap()
    )
    .is_err());
    fs::write(beds.join("one.bed"), "unknown#0#name\t0\t3000\n").unwrap();
    assert!(catalog::build(
        &panel,
        identity.clone(),
        catalog::import_beds(&beds).unwrap()
    )
    .is_err());
    fs::write(beds.join("one.bed"), "A#0#chr1\t0\t1\n").unwrap();
    let unsupported = catalog::build(
        &panel,
        identity.clone(),
        catalog::import_beds(&beds).unwrap(),
    )
    .unwrap();
    assert_eq!(unsupported.occurrences.len(), 1);
    assert_eq!(unsupported.occurrences[0].fully_contained_anchors, 0);
    let unsupported_calls = genome::calling::call(&unsupported, &sample, true, 1).unwrap();
    assert_eq!(unsupported_calls.calls[0].unsupported_occurrences, 1);
    assert_eq!(
        unsupported_calls.calls[0].status,
        "no-call-no-local-features"
    );

    let linked: catalog::CatalogInput = serde_json::from_value(json!({"version":1,"groups":[{
        "id":"neighbors","scaffold":null,"occurrences":[
            interval(specs[0].0,0,1000,json!("+")), interval(specs[0].0,0,1000,json!("+")),
            interval(specs[0].0,900,1800,json!("+")), interval(specs[0].0,1900,2000,json!("+"))]
    }]}))
    .unwrap();
    let linked = catalog::build(&panel, identity.clone(), linked).unwrap();
    assert_eq!(
        linked
            .links
            .iter()
            .filter(|l| l.relation == "overlap" && l.to == 2)
            .count(),
        2
    );
    assert!(linked
        .links
        .iter()
        .any(|l| l.from == 2 && l.to == 3 && l.source_gap_bp == 100 && l.oriented_continuation));
    let mut incompatible = SampleIndex::load(&run.join("sample.membwt"), &identity).unwrap();
    incompatible.version += 1;
    let versioned = root.join("versioned.membwt");
    incompatible.save(&versioned).unwrap();
    assert!(SampleIndex::load(&versioned, &identity).is_err());
    let malformed_reads = root.join("malformed.fastq");
    fs::write(&malformed_reads, "@r\nACGT\n+\n!\n").unwrap();
    let failed_reads = root.join("failed-reads");
    assert!(!command(&[
        "genome-infer",
        "build-sample",
        "--panel",
        path(&panel_path),
        "--reads",
        path(&malformed_reads),
        "--out-dir",
        path(&failed_reads)
    ])
    .status
    .success());
    assert_eq!(
        json_file(&failed_reads.join("manifest.json"))["status"],
        "failed"
    );
    assert!(!failed_reads.join("sample.membwt").exists());
    #[cfg(unix)]
    {
        let symlink = root.join("symlink-output");
        std::os::unix::fs::symlink(&run, &symlink).unwrap();
        assert!(!command(&[
            "genome-infer",
            "build-sample",
            "--panel",
            path(&panel_path),
            "--reads",
            path(&reads),
            "--out-dir",
            path(&symlink)
        ])
        .status
        .success());
        assert_eq!(json_file(&run.join("manifest.json"))["status"], "succeeded");
    }
}
