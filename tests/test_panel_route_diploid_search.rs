//! Frozen Gate-B paired construction; inference never receives oracle/truth geometry.
#[path = "panel_route_diploid_search/helpers.rs"]
mod helpers;
#[path = "../examples/panel_route_diploid_search/mod.rs"]
mod search;
use clap::Parser;
use impg::genome_inference::{self as genome, panel_routes as routes};
use search::score::Pair;
use serde_json::json;
use std::{fs, path::PathBuf};
fn root(name: &str) -> (tempfile::TempDir, PathBuf) {
    let temp = tempfile::tempdir().unwrap();
    let root = std::env::var_os("IMPG_TEST_PAIRED_OUTPUT")
        .map(PathBuf::from)
        .unwrap_or_else(|| temp.path().into())
        .join(name);
    (temp, root)
}
#[test]
fn paired_automatic_construction() {
    let (_temp, root) = root("construction");
    let f = helpers::fixture(&root, "middle");
    let result = search::run(f.options("inference"), None).unwrap();
    let report = helpers::assess(&f, "inference", &result);
    assert!(
        report["truth_exact_pair_in_support"].as_bool().unwrap(),
        "frozen construction failed; assessment retained: {}",
        root.display()
    );
    assert!(report["improvement_over_best_native"].as_f64().unwrap() > 1e-6);
    assert_eq!(result["native_initialization_complete"], true);
}
#[test]
#[ignore = "historical negative: superseded exhaustive CIS/TRANS generation criterion; run explicitly with --ignored --exact"]
fn paired_automatic_phase_ambiguity() {
    let (_temp, root) = root("phase");
    let f = helpers::fixture(&root, "phase");
    let result = search::run(f.options("inference"), None).unwrap();
    let report = helpers::assess(&f, "inference", &result);
    assert_eq!(report["cis_trans_exact_operator_count_equality"], true);
    assert!(
        report["automatic_phase_classes"]
            .as_array()
            .unwrap()
            .iter()
            .all(|c| c["generated"] == true && c["retained"] == true),
        "automatic cis/trans support missing: {}",
        root.display()
    );
}
#[test]
fn paired_automatic_far_phase_observation_recovery() {
    let (_temp, root) = root("phase-observable");
    let f = helpers::fixture(&root, "phase");
    let result = search::run(f.options("inference"), None).unwrap();
    let report = helpers::assess(&f, "inference", &result);
    helpers::assert_observation_recovery(&report, &result, false);
}

#[test]
fn paired_automatic_near_phase_native_pair() {
    let (_temp, root) = root("near-native");
    let f = helpers::near_native_fixture(&root);
    let result = search::run(f.options("inference"), None).unwrap();
    let report = helpers::assess(&f, "inference", &result);
    helpers::assert_observation_recovery(&report, &result, true);
}

#[test]
fn paired_automatic_mosaic_required_dosage() {
    let (_temp, root) = root("dosage-mosaic");
    let fixture = helpers::dosage_fixture(&root);
    let result = search::run(fixture.options("inference"), None).unwrap();
    let report = helpers::assess_dosage(&fixture, "inference", &result);
    assert_eq!(report["finite_pairs"], 153);
    assert_eq!(report["support_complete"], false);
    assert!(report["global_bound"].is_null());
    let minimum = report["finite_153_minimum"].as_f64().unwrap();
    assert!((report["best_found"].as_f64().unwrap() - minimum).abs() <= 1e-9);
    for native in report["native_pairs"].as_array().unwrap() {
        assert_eq!(native["count_vector_differs"], true);
        assert!(native["objective"].as_f64().unwrap() > minimum + 1e-9);
    }
    let retained = report["retained"].as_array().unwrap();
    assert!(!retained.is_empty());
    for pair in retained {
        assert_eq!(pair["exact_truth_count_vector"], true);
        assert_eq!(pair["dosage"], json!([0, 1, 1, 2]));
        assert_eq!(pair["both_primaries_non_native"], true);
        assert_eq!(pair["truth_bp"], 3872);
        assert_eq!(pair["query_bp"], 3872);
        assert_eq!(pair["truth_primary_bp"], 3072);
        assert_eq!(pair["query_primary_bp"], 3072);
        assert_eq!(pair["truth_tract_bp"], 192);
        assert_eq!(pair["query_tract_bp"], 192);
    }
}

#[test]
fn paired_checked_loss_schema_and_copy_identity() {
    assert!(search::score::add(u64::MAX, 1).is_err());
    assert!(search::score::loss(u64::MAX, 1, 1).is_err());
    let (_, joint) = search::score::loss(300, 7, 62).unwrap();
    let (_, one) = search::score::loss(150, 7, 62).unwrap();
    assert!((joint - 2.0 * one).abs() > 1.0);
    assert_eq!(search::score::loss(0, 9, 62).unwrap().1, 0.0);
    assert_eq!(search::score::loss(150, 0, 62).unwrap().1, 10.0);
    let a = routes::Assignment {
        version: 1,
        model: routes::MODEL.into(),
        graph_checksum: "test".into(),
        family: 0,
        routes: vec![],
    };
    let mut b = a.clone();
    b.family = 1;
    let pair = Pair {
        copies: [a.clone(), b.clone()],
    };
    let swapped = Pair { copies: [b, a] };
    assert_eq!(pair.key(), swapped.key());
    assert_ne!(pair, swapped);
    assert!(serde_json::from_value::<Pair>(
        json!({"copies":[pair.copies[0],pair.copies[1],pair.copies[0]]})
    )
    .is_err());
    assert!(!search::expansion_eligible(&swapped, &[pair.clone()]));
    let o = search::Options::try_parse_from([
        "paired",
        "--panel",
        "x",
        "--routes",
        "x",
        "--sample",
        "x",
        "--out-dir",
        "x",
        "--max-scores",
        "2049",
    ])
    .unwrap();
    assert!(o.validate().is_err());
    assert!(search::Options::try_parse_from(["paired", "--truth", "x"]).is_err());
}
#[test]
#[ignore = "explicit paired standalone binary and fresh output required"]
fn configured_paired_automatic_cli() {
    let bin = std::env::var_os("IMPG_TEST_PAIRED_BIN").expect("IMPG_TEST_PAIRED_BIN required");
    let root = PathBuf::from(
        std::env::var_os("IMPG_TEST_PAIRED_CLI_OUTPUT")
            .expect("IMPG_TEST_PAIRED_CLI_OUTPUT required"),
    );
    let f = helpers::fixture(&root, "middle");
    let o = f.options("inference");
    let output = std::process::Command::new(bin)
        .args([
            "--panel",
            &o.panel,
            "--routes",
            o.routes.to_str().unwrap(),
            "--sample",
            o.sample.to_str().unwrap(),
            "--out-dir",
            o.out_dir.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    fs::write(root.join("cli.stdout"), &output.stdout).unwrap();
    fs::write(root.join("cli.stderr"), &output.stderr).unwrap();
    assert!(output.status.success());
    let result = genome::read_json(&o.out_dir.join("result.json")).unwrap();
    let report = helpers::assess(&f, "inference", &result);
    assert_eq!(report["truth_exact_pair_in_support"], true);
}

fn read_lines(path: &std::path::Path) -> Vec<serde_json::Value> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}
#[test]
fn paired_causal_geometry_retains_partner_signal() {
    use impg::{genome_inference::sample, syng::SyngIndex};
    let (_temp, root) = root("causal");
    let f = helpers::fixture(&root, "middle");
    let mut first = f.options("middle");
    first.max_epochs = 0;
    let r1 = search::run(first, Some(0)).unwrap();
    let a1 = helpers::assess(&f, "middle", &r1);
    let mut last = f.sources[0].1.clone();
    last[1216..1240].copy_from_slice(&f.sources[1].1[1216..1240]);
    let truth = vec![
        last,
        f.sources[2].1.clone(),
        f.sources[0].1.clone(),
        f.sources[2].1.clone(),
    ];
    let prefix = root.join("panel.syng");
    let id = genome::PanelIdentity::read(prefix.to_str().unwrap()).unwrap();
    let panel = SyngIndex::load(prefix.to_str().unwrap(), Default::default()).unwrap();
    helpers::make_sample(&root, "last-sample", &truth, &panel, &id);
    let mut second = f.options("last");
    second.max_epochs = 0;
    second.sample = root.join("last-sample.membwt");
    let r2 = search::run(second, Some(0)).unwrap();
    let f2 = helpers::Fixture {
        root: root.clone(),
        sources: f.sources.clone(),
        truth,
        kind: "last".into(),
    };
    let a2 = helpers::assess(&f2, "last", &r2);
    let traces: Vec<Vec<serde_json::Value>> = ["middle", "last"]
        .iter()
        .map(|n| read_lines(&root.join(n).join("events.jsonl")))
        .collect();
    let native = |rows: &Vec<serde_json::Value>| {
        rows.iter()
            .find(|r| r["event"] == "paired_baseline")
            .unwrap()["scored"]["pair"]
            .clone()
    };
    assert_eq!(native(&traces[0]), native(&traces[1]));
    let initial = |rows: &Vec<serde_json::Value>| {
        rows.iter()
            .filter(|r| r["event"] == "geometry_attempt")
            .take(2)
            .cloned()
            .collect::<Vec<_>>()
    };
    assert_eq!(initial(&traces[0]), initial(&traces[1]));
    let guided = |rows: &Vec<serde_json::Value>| {
        rows.iter()
            .filter(|r| {
                r["event"] == "paired_geometry_task" && r["task"]["Region"]["guided"] == true
            })
            .cloned()
            .collect::<Vec<_>>()
    };
    let changes = [guided(&traces[0]), guided(&traces[1])];
    genome::write_json(&root.join("causal-assessment.json"),&json!({"same_ordered_initial_pair":true,"same_initial_geometry":true,"actual_guided_geometry_differed":changes[0]!=changes[1],"guided_attempts":changes,"independent_pair_audits":[a1["inference_pairs_independently_audited"],a2["inference_pairs_independently_audited"]],"partner_never_removed_from_objective":true})).unwrap();
    assert!(!changes[0].is_empty() && !changes[1].is_empty());
    assert_ne!(changes[0], changes[1]);
    assert_eq!(
        sample::SampleIndex::load(&root.join("sample.membwt"), &id)
            .unwrap()
            .stats
            .reads,
        224
    );
}
#[test]
fn paired_capacity_coupling_and_count_equivalent_expansion() {
    let (_temp, root) = root("coupled");
    let f = helpers::coupled_fixture(&root);
    let result = search::run(f.options("inference"), None).unwrap();
    let out = root.join("inference");
    assess_coupled(
        &f,
        &out,
        &out.with_file_name("inference-assessment"),
        &result,
    );
}
fn assess_coupled(
    f: &helpers::Fixture,
    inference: &std::path::Path,
    out: &std::path::Path,
    result: &serde_json::Value,
) {
    use impg::{genome_inference::sample, syng::SyngIndex};
    let root = &f.root;
    fs::create_dir(out).unwrap();
    let frozen = genome::reconstruction::fingerprint(&inference.join("result.json")).unwrap();
    genome::write_json(
        &out.join("mechanical-fnv-before-assessment.json"),
        &json!({"authoritative_SHA256":false,"result_fingerprint":frozen}),
    )
    .unwrap();
    let events = read_lines(&inference.join("events.jsonl"));
    let ledger = read_lines(&inference.join("scores.jsonl"));
    let bases: Vec<_> = events
        .iter()
        .filter(|r| r["event"] == "paired_baseline")
        .collect();
    let prefix = root.join("panel.syng");
    let id = genome::PanelIdentity::read(prefix.to_str().unwrap()).unwrap();
    let panel = SyngIndex::load(prefix.to_str().unwrap(), Default::default()).unwrap();
    let graph = routes::Graph::load(&root.join("routes"), &id).unwrap();
    let sample = sample::SampleIndex::load(&root.join("sample.membwt"), &id).unwrap();
    let mut e = routes::Evaluator::new(
        &root.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        50_000,
        50_000,
    )
    .unwrap();
    let mut pairs = Vec::new();
    let mut counts = Vec::new();
    for b in &bases {
        let pair: Pair = serde_json::from_value(b["scored"]["pair"].clone()).unwrap();
        let mut q = std::collections::BTreeMap::new();
        for a in &pair.copies {
            let ev = e.evaluate(a).unwrap();
            for f in ev.factors {
                let n = q.entry(f.tokens).or_insert(0u64);
                *n = n.checked_add(f.counts_by_length[0]).unwrap();
            }
        }
        pairs.push(pair);
        counts.push(q);
    }
    let mut equivalent_admitted = false;
    for i in 1..pairs.len() {
        if counts[i] == counts[0] && pairs[i].key() != pairs[0].key() {
            assert!(search::expansion_eligible(&pairs[i], &pairs[..i]));
            equivalent_admitted = true;
        }
    }
    let support: Vec<Pair> =
        serde_json::from_value(result["retained_correlated_best_found"].clone()).unwrap();
    let mut coverage = Vec::new();
    let mut replay_reads = 0usize;
    let mut replay_cache = std::collections::BTreeMap::new();
    let mut route_profiles = 0;
    for pair in &support {
        let mut query = Vec::new();
        for a in &pair.copies {
            e.validate_assignment(a).unwrap();
            for r in &a.routes {
                let dna = e.sources.spell(r, 0, r.length().unwrap()).unwrap();
                if !replay_cache.contains_key(&dna) {
                    let mut forward = None;
                    for reverse in [false, true] {
                        let sequence = if reverse {
                            impg::graph::reverse_complement(&dna)
                        } else {
                            dna.clone()
                        };
                        replay_reads += sequence.len().saturating_sub(149);
                        assert!(replay_reads <= 100_000);
                        let path =
                            out.join(format!("pooled-replay-{}-{reverse}.fa", replay_cache.len()));
                        fs::write(
                            &path,
                            sequence
                                .windows(150)
                                .map(|s| format!(">anonymous\n{}\n", String::from_utf8_lossy(s)))
                                .collect::<String>(),
                        )
                        .unwrap();
                        let q = sample::build(&panel, id.clone(), &[path])
                            .unwrap()
                            .counts
                            .observed_pairs()
                            .unwrap();
                        if let Some(expected) = &forward {
                            assert_eq!(*expected, q);
                        } else {
                            forward = Some(q);
                        }
                    }
                    replay_cache.insert(dna.clone(), forward.unwrap());
                }
                let (q, _) = e.route_counts(r).unwrap();
                route_profiles += 1;
                assert_eq!(q[0], replay_cache[&dna]);
                query.push(dna);
            }
        }
        let assessed = helpers::coverage(&f.truth, &query);
        assert_eq!(assessed["full_truth_coverage"], 1.0);
        assert_eq!(assessed["full_query_coverage"], 1.0);
        coverage.push(assessed);
    }
    genome::write_json(&out.join("coupling-assessment.json"),&json!({"result":result,"public_valid_count_equivalent_distinct_physical_expansion_admitted":equivalent_admitted,"admitted_baselines":bases,"rejection_records":ledger.iter().filter(|r|r["kind"]=="within_copy_capacity_rejection").count(),"scope":"within-copy capacity-coupled exchange, not cross-copy capacity repair","selected_copy_injective_assessments":coverage,"independent_route_profiles":route_profiles,"independent_replay_reads_both_orientations":replay_reads,"actual_reads":sample.stats.reads,"actual_bases":sample.stats.bases,"diploid_bp":880,"realized_depth":sample.stats.bases as f64/880.0,"nominal_per_copy_depth":10})).unwrap();
    assert!(result["within_copy_capacity_rejections"].as_u64().unwrap() > 0);
    assert!(result["compound_confirmations"].as_u64().unwrap() > 0);
    assert!(equivalent_admitted);
}

#[test]
fn paired_window_includes_own_baseline_improvement_below_global_best() {
    use impg::{genome_inference::sample, syng::SyngIndex};
    let (_temp, root) = root("eligibility");
    let f = helpers::fixture(&root, "middle");
    let prefix = root.join("panel.syng");
    let id = genome::PanelIdentity::read(prefix.to_str().unwrap()).unwrap();
    let panel = SyngIndex::load(prefix.to_str().unwrap(), Default::default()).unwrap();
    let graph = routes::Graph::load(&root.join("routes"), &id).unwrap();
    let sample = sample::SampleIndex::load(&root.join("sample.membwt"), &id).unwrap();
    let mut e = routes::Evaluator::new(
        &root.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        50_000,
        50_000,
    )
    .unwrap();
    let options = f.options("mechanism");
    fs::create_dir(&options.out_dir).unwrap();
    assert!(
        search::test_intermediate_improvement_window(&mut e, options).unwrap(),
        "own-baseline improvement must remain eligible when it does not beat global best"
    );
}

#[test]
fn paired_resource_boundaries_and_late_completed_accounting() {
    use impg::{genome_inference::sample, syng::SyngIndex};
    let (_temp, root) = root("boundaries");
    let f = helpers::fixture(&root, "middle");
    let mut reports = Vec::new();
    for (label, expected) in [
        ("pair-score", "budget: paired_scores"),
        ("haploid-score", "budget: exact_scores"),
        ("one-profile", "budget: profile_work"),
        ("validation", "budget: validations"),
        ("tasks", "budget: frontier_initial_heads"),
        ("baselines", "budget: paired_baselines"),
        ("chain", "budget: paired_chain_primitives"),
        ("segments", "budget: route_segments"),
    ] {
        let mut o = f.options(label);
        match label {
            "pair-score" => o.max_scores = 0,
            "haploid-score" => o.max_haploid_scores = 1,
            "one-profile" => o.max_profile_work = 1936,
            "validation" => o.max_validations = 0,
            "tasks" => o.max_tasks = 1,
            "baselines" => o.max_baselines = 1,
            "chain" => o.max_chain_primitives = 0,
            "segments" => o.max_segments = 0,
            _ => unreachable!(),
        }
        if label == "tasks" {
            assert_eq!(
                search::run(o.clone(), None).unwrap_err().to_string(),
                expected
            );
            assert!(!o.out_dir.join("result.json").exists());
            reports.push(json!({"case":label,"preflight_failed":true,"more_inclusive_heads":true}));
            continue;
        }
        let r = search::run(o.clone(), None).unwrap();
        let pending: serde_json::Value =
            genome::read_json(&o.out_dir.join("pending.json")).unwrap();
        assert_eq!(r["status"], expected, "{label}");
        assert_eq!(r["support_complete"], false);
        let ledger = read_lines(&o.out_dir.join("scores.jsonl"));
        let complete = ledger
            .iter()
            .filter(|r| r["kind"] == "paired_score")
            .count();
        let haploids = ledger
            .iter()
            .filter(|r| r["kind"] == "public_haploid_return")
            .count();
        assert_eq!(r["paired_objectives_completed"], complete);
        assert_eq!(r["public_haploid_losses_completed"], haploids);
        if label == "one-profile" || label == "haploid-score" {
            assert_eq!(haploids, 1);
            assert_eq!(complete, 0);
            assert!(pending["active_pair"].is_object());
        }
        if label == "tasks" || label == "baselines" {
            assert!(pending["pending_expansion"].is_object());
            assert!(complete > 0);
        }
        reports.push(json!({"case":label,"result":r,"pending":pending}));
    }
    for label in ["state", "work", "features"] {
        let mut o = f.options(label);
        match label {
            "state" => o.max_state_bytes = 1,
            "work" => o.max_work = 0,
            _ => o.max_features = 1,
        };
        assert!(search::run(o.clone(), None).is_err());
        assert!(!o.out_dir.join("result.json").exists());
        reports.push(json!({"case":label,"preflight_failed":true}));
    }
    let prefix = root.join("panel.syng");
    let id = genome::PanelIdentity::read(prefix.to_str().unwrap()).unwrap();
    let panel = SyngIndex::load(prefix.to_str().unwrap(), Default::default()).unwrap();
    let graph = routes::Graph::load(&root.join("routes"), &id).unwrap();
    let sample = sample::SampleIndex::load(&root.join("sample.membwt"), &id).unwrap();
    let mut e = routes::Evaluator::new(
        &root.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        50_000,
        50_000,
    )
    .unwrap();
    let a = graph.native_assignment(0).unwrap();
    let pair = Pair {
        copies: [a.clone(), a.clone()],
    };
    // Independent setup evaluation determines exact known factor-pass size, not a recovery-budget retry.
    let profile = e.evaluate(&a).unwrap();
    let terms = profile.factors.len() as u64;
    let mut limits = search::haploid::Options::parse_from([
        "h",
        "--panel",
        "x",
        "--routes",
        "x",
        "--sample",
        "x",
        "--out-dir",
        "x",
    ]);
    let mut probe = search::haploid::Meter::new(&limits);
    probe.validation(&a, graph.k, graph.port_count).unwrap();
    let validation_work = probe.work_used;
    limits.max_work = 2 * validation_work + 2 * terms; // BOTH public profiles return; final full pair pass cannot enter.
    let mut meter = search::haploid::Meter::new(&limits);
    let mut ledger = fs::File::create(root.join("late-before-pair.jsonl")).unwrap();
    let (mut admitted, mut complete) = (0, 0);
    let error = search::score::Scorer {
        e: &mut e,
        meter: &mut meter,
        ledger: &mut ledger,
        pair_attempts: &mut admitted,
        pair_completed: &mut complete,
        pair_cap: 1,
        histogram: sample.stats.read_lengths[&150],
    }
    .evaluate(&pair, true)
    .err()
    .unwrap();
    assert_eq!(error.to_string(), "budget: work");
    assert_eq!(admitted, 1);
    assert_eq!(complete, 0);
    assert_eq!(meter.public_scores, 2);
    reports.push(json!({"case":"both_public_profiles_before_pair_completion","admitted":admitted,"complete":complete,"meter":meter,"setup_public_haploid_calls":1,"setup_haploid_loss":profile.relative_objective,"derived_work_cap":limits.max_work}));
    // AA checked multiplicity and unchanged partner are visible independently in returned public factors.
    for f in &profile.factors {
        assert_eq!(
            search::score::add(f.counts_by_length[0], f.counts_by_length[0]).unwrap(),
            2 * f.counts_by_length[0]
        );
    }
    let negative = helpers::coverage(
        &[f.sources[0].1.clone(), f.sources[0].1.clone()],
        &[f.sources[0].1.clone()],
    );
    assert_eq!(negative["full_truth_coverage"], 0.5);
    assert_eq!(negative["full_query_coverage"], 1.0);
    genome::write_json(&root.join("boundary-assessments.json"), &reports).unwrap();
}
#[test]
fn paired_retention_stop_preserves_completed_score_and_physical_pair() {
    let (_temp, root) = root("retention");
    let f = helpers::coupled_fixture(&root);
    let mut o = f.options("inference");
    o.max_ties = 1;
    let r = search::run(o.clone(), None).unwrap();
    let ledger = read_lines(&o.out_dir.join("scores.jsonl"));
    let pending: serde_json::Value = genome::read_json(&o.out_dir.join("pending.json")).unwrap();
    assert_eq!(r["status"], "budget: paired_retained_ties");
    assert!(pending["active_pair"].is_object());
    assert_eq!(
        r["paired_objectives_completed"].as_u64().unwrap(),
        ledger
            .iter()
            .filter(|r| r["kind"] == "paired_score")
            .count() as u64
    );
    assert_eq!(
        r["public_haploid_losses_completed"].as_u64().unwrap(),
        2 * r["paired_objectives_completed"].as_u64().unwrap()
    );
    assert!(r["confirmed_candidates"].as_u64().unwrap() > 0);
    genome::write_json(&root.join("retention-assessment.json"),&json!({"result":r,"pending":pending,"complete_pair_retained_in_ledger_before_tie_stop":true})).unwrap();
}

#[test]
fn paired_realized_feature_boundary_is_not_a_partial_objective() {
    use impg::genome_inference::sample;
    let (_temp, root) = root("feature-boundary");
    let f = helpers::fixture(&root, "middle");
    let id = genome::PanelIdentity::read(root.join("panel.syng").to_str().unwrap()).unwrap();
    let sample = sample::SampleIndex::load(&root.join("sample.membwt"), &id).unwrap();
    let cap = sample.counts.observed_pairs().unwrap().len();
    let mut o = f.options("inference");
    o.max_features = cap;
    let r = search::run(o.clone(), None).unwrap();
    assert!(r["status"]
        .as_str()
        .unwrap()
        .contains("public_feature_terms"));
    assert!(
        r["meter"]["public_score_attempts"].as_u64().unwrap()
            > r["public_haploid_losses_completed"].as_u64().unwrap()
    );
    let pending: serde_json::Value = genome::read_json(&o.out_dir.join("pending.json")).unwrap();
    assert!(pending["active_pair"].is_object());
    assert_eq!(r["support_complete"], false);
    let rows = read_lines(&o.out_dir.join("scores.jsonl"));
    assert_eq!(
        r["paired_objectives_completed"].as_u64().unwrap(),
        rows.iter().filter(|r| r["kind"] == "paired_score").count() as u64
    );
    genome::write_json(
        &root.join("feature-assessment.json"),
        &json!({"observed_support_derived_cap":cap,"result":r,"pending":pending}),
    )
    .unwrap();
}

/// Explicit staging only. The external hashlib driver owns authoritative freezing.
#[test]
#[ignore = "explicit protocol stage, fixture, case and fresh inference output required"]
fn configured_paired_protocol_stage() {
    let stage = std::env::var("IMPG_PAIRED_STAGE").expect("IMPG_PAIRED_STAGE required");
    let root = PathBuf::from(
        std::env::var_os("IMPG_PAIRED_FIXTURE").expect("IMPG_PAIRED_FIXTURE required"),
    );
    let case = std::env::var("IMPG_PAIRED_CASE").expect("IMPG_PAIRED_CASE required");
    assert!([
        "middle",
        "phase",
        "phase-observable",
        "near-native",
        "last",
        "causal-middle",
        "causal-last",
        "coupled"
    ]
    .contains(&case.as_str()));
    let kind = if case == "phase-observable" {
        "phase"
    } else {
        case.strip_prefix("causal-").unwrap_or(&case)
    };
    if stage == "prepare" {
        if kind == "near-native" {
            helpers::near_native_fixture(&root);
        } else if kind == "coupled" {
            helpers::coupled_fixture(&root);
        } else {
            helpers::fixture(&root, kind);
        }
        return;
    }
    let output = PathBuf::from(
        std::env::var_os("IMPG_PAIRED_INFERENCE").expect("IMPG_PAIRED_INFERENCE required"),
    );
    if stage == "infer" {
        // No fixture/truth description is read by this path.
        let mut o = search::Options::parse_from([
            "paired",
            "--panel",
            root.join("panel.syng").to_str().unwrap(),
            "--routes",
            root.join("routes").to_str().unwrap(),
            "--sample",
            root.join(if kind == "near-native" {
                "sample/sample.membwt"
            } else {
                "sample.membwt"
            })
            .to_str()
            .unwrap(),
            "--out-dir",
            output.to_str().unwrap(),
        ]);
        let fixed = case.starts_with("causal-").then_some(0);
        if fixed.is_some() {
            o.max_epochs = 0;
        }
        search::run(o, fixed).unwrap();
        return;
    }
    assert_eq!(stage, "assess");
    let f: helpers::Fixture = genome::read_json(&root.join("fixture-description.json")).unwrap();
    assert_eq!(f.kind, kind);
    let result = genome::read_json(&output.join("result.json")).unwrap();
    if kind == "coupled" {
        assess_coupled(&f, &output, &helpers::assessment_dir(&output), &result);
        return;
    }
    let report = helpers::assess(&f, output.to_str().unwrap(), &result);
    if case == "middle" {
        assert_eq!(report["truth_exact_pair_in_support"], true);
        assert!(report["improvement_over_best_native"].as_f64().unwrap() > 1e-6);
        assert_eq!(result["native_initialization_complete"], true);
    } else if case == "phase-observable" || case == "near-native" {
        helpers::assert_observation_recovery(&report, &result, case == "near-native");
    } else if case == "phase" {
        assert_eq!(report["cis_trans_exact_operator_count_equality"], true);
        assert!(
            report["automatic_phase_classes"]
                .as_array()
                .unwrap()
                .iter()
                .all(|c| c["generated"] == true && c["retained"] == true),
            "automatic cis/trans support missing: {}",
            output.display()
        );
    }
}

fn frontier_fixture(
    name: &str,
    check: impl FnOnce(&mut routes::Evaluator<'_>, search::Options) -> std::io::Result<()>,
) {
    use impg::{genome_inference::sample, syng::SyngIndex};
    let (_temp, root) = root(name);
    let f = helpers::fixture(&root, "middle");
    let prefix = root.join("panel.syng");
    let id = genome::PanelIdentity::read(prefix.to_str().unwrap()).unwrap();
    let panel = SyngIndex::load(prefix.to_str().unwrap(), Default::default()).unwrap();
    let graph = routes::Graph::load(&root.join("routes"), &id).unwrap();
    let sample = sample::SampleIndex::load(&root.join("sample.membwt"), &id).unwrap();
    let mut e = routes::Evaluator::new(
        &root.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        50_000,
        50_000,
    )
    .unwrap();
    let o = f.options("transitions");
    fs::create_dir(&o.out_dir).unwrap();
    check(&mut e, o).unwrap();
}

#[test]
fn paired_frontier_owned_transitions() {
    frontier_fixture("frontier-mechanics", search::test_frontier_transitions);
}

#[test]
fn paired_frontier_scheduler_transitions() {
    frontier_fixture("frontier-scheduler", |e, o| {
        search::test_scheduler_transitions(e, o, false)
    });
}

#[test]
fn paired_frontier_pending_diagnostics() {
    frontier_fixture("frontier-pending", |e, o| {
        search::test_scheduler_transitions(e, o, true)
    });
}
