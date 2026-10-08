//! Frozen synthetic B1 controls; truth is used only for read generation/assessment.
#[path = "../examples/panel_route_residual_search/mod.rs"]
mod search;
use clap::Parser;
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::{SyncmerParams, SyngIndex},
};
use serde_json::{json, Value};
use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};
fn dna(n: usize, mut s: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            b"ACGT"[(s & 3) as usize]
        })
        .collect()
}
fn write(p: &Path, v: &impl serde::Serialize) {
    genome::write_json(p, v).unwrap();
}
fn read(p: &Path) -> Value {
    genome::read_json(p).unwrap()
}
fn fasta(p: &Path, s: &[(String, Vec<u8>)]) {
    fs::write(
        p,
        s.iter()
            .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
}
fn command(args: &[&str]) {
    let o = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(args)
        .output()
        .unwrap();
    assert!(o.stdout.len() + o.stderr.len() < 1024 * 1024);
    assert!(
        o.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
}
fn p(p: &Path) -> &str {
    p.to_str().unwrap()
}
struct Fixture {
    root: PathBuf,
    sources: Vec<(String, Vec<u8>)>,
    truth: Vec<Vec<u8>>,
}
fn setup(root: &Path, coupled: bool) -> Fixture {
    setup_variant(root, coupled, false, false)
}
fn setup_variant(
    root: &Path,
    coupled: bool,
    reverse_second: bool,
    separate_families: bool,
) -> Fixture {
    fs::create_dir(root).unwrap();
    let (mut sources, truth) = if coupled {
        let a = dna(220, 113);
        let b = if reverse_second {
            impg::graph::reverse_complement(&a)
        } else {
            a.clone()
        };
        (
            vec![("A#0#one".into(), a.clone()), ("A#0#two".into(), b.clone())],
            vec![a, b],
        )
    } else {
        let a = dna(1536, 219);
        let mut b = a.clone();
        for start in [256, 704, 1216] {
            for x in &mut b[start..start + 24] {
                *x = if *x == b'A' { b'C' } else { b'A' };
            }
        }
        let mut truth = a.clone();
        truth[704..728].copy_from_slice(&b[704..728]);
        let passive = dna(400, 57);
        (
            vec![
                ("A#0#one".into(), a),
                ("B#0#one".into(), b),
                ("A#0#two".into(), passive.clone()),
                ("B#0#two".into(), dna(400, 58)),
            ],
            vec![truth, passive],
        )
    };
    if separate_families {
        sources[1].0 = "B#0#two".into();
        sources.push(("C#0#three".into(), sources[0].1.clone()));
    }
    let mut panel = SyngIndex::build(SyncmerParams::default(), sources.clone().into_iter());
    let prefix = root.join("panel.syng");
    panel.save(p(&prefix)).unwrap();
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    let fa = root.join("sources.fa");
    fasta(&fa, &sources);
    let groups = root.join("groups.json");
    write(
        &groups,
        &json!({"version":1,"groups":sources.iter().enumerate().map(|(i,(n,s))|json!({"id":format!("g{i}"),"scaffold":null,"occurrences":[{"path":n,"start":0,"end":s.len(),"strand":null}]})).collect::<Vec<_>>()}),
    );
    command(&[
        "genome-infer",
        "build-catalog",
        "--panel",
        p(&prefix),
        "--groups",
        p(&groups),
        "--out-dir",
        p(&root.join("catalog")),
    ]);
    command(&[
        "genome-infer",
        "build-panel-routes",
        "--panel",
        p(&prefix),
        "--catalog",
        p(&root.join("catalog/catalog.json")),
        "--sources",
        p(&fa),
        "--read-lengths",
        "150",
        "--core-bp",
        "37",
        "--out-dir",
        p(&root.join("routes")),
    ]);
    make_sample(root, "sample", &truth, &panel, &identity);
    write(
        &root.join("recipe.json"),
        &json!({"sources":sources.iter().map(|(n,s)|json!({"name":n,"bp":s.len()})).collect::<Vec<_>>(),"read_length":150,"step":15,"depth":10,"background":0.1,"truth_not_in_inference_arguments":true}),
    );
    Fixture {
        root: root.into(),
        sources,
        truth,
    }
}
fn make_sample(
    root: &Path,
    name: &str,
    molecules: &[Vec<u8>],
    panel: &SyngIndex,
    id: &genome::PanelIdentity,
) {
    let mut reads = Vec::new();
    for (m, s) in molecules.iter().enumerate() {
        let last = s.len() - 150;
        let mut starts = (0..=last).step_by(15).collect::<Vec<_>>();
        if starts.last() != Some(&last) {
            starts.push(last);
        }
        for i in starts {
            reads.push((format!("{m}-{i}"), s[i..i + 150].to_vec()));
        }
    }
    let fa = root.join(format!("{name}.fa"));
    fasta(&fa, &reads);
    sample::build(panel, id.clone(), &[fa])
        .unwrap()
        .save(&root.join(format!("{name}.membwt")))
        .unwrap();
}
fn options(f: &Fixture, name: &str) -> search::Options {
    search::Options::parse_from([
        "search",
        "--panel",
        p(&f.root.join("panel.syng")),
        "--routes",
        p(&f.root.join("routes")),
        "--sample",
        p(&f.root.join("sample.membwt")),
        "--out-dir",
        p(&f.root.join(name)),
    ])
}
fn lines(p: &Path) -> Vec<Value> {
    fs::read_to_string(p)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}
fn spell(a: &routes::Assignment, sources: &[(String, Vec<u8>)]) -> Vec<Vec<u8>> {
    a.routes
        .iter()
        .map(|r| {
            r.segments
                .iter()
                .flat_map(|p| {
                    let s = sources[p.source].1[p.start as usize..p.end as usize].to_vec();
                    if p.reverse {
                        impg::graph::reverse_complement(&s)
                    } else {
                        s
                    }
                })
                .collect()
        })
        .collect()
}
fn freeze(out: &Path) {
    let rows:Vec<_>=["scores.jsonl","events.jsonl","result.json","pending.json","accounting.jsonl","provenance.json","options.json","status.json"].iter().map(|n|json!({"file":n,"fingerprint":genome::reconstruction::fingerprint(&out.join(n)).unwrap()})).collect();
    write(&out.join("frozen-before-dna-assessment.json"), &rows);
}
fn assess(f: &Fixture, name: &str) {
    let out = f.root.join(name);
    freeze(&out);
    let result = read(&out.join("result.json"));
    let ledger = lines(&out.join("scores.jsonl"));
    let natives: Vec<_> = ledger.iter().filter(|r| r["kind"] == "native").collect();
    assert_eq!(natives.len(), 2);
    let best_native = natives
        .iter()
        .map(|n| n["public"]["relative_objective"].as_f64().unwrap())
        .fold(f64::INFINITY, f64::min);
    let best = result["best_public_objective"].as_f64().unwrap();
    let ties = result["retained_correlated_best_found"].as_array().unwrap();
    let dna_matches = ties
        .iter()
        .filter(|a| {
            spell(
                &serde_json::from_value::<routes::Assignment>((*a).clone()).unwrap(),
                &f.sources,
            ) == f.truth
        })
        .count();
    write(
        &out.join("assessment.json"),
        &json!({"best_native":best_native,"best_found":best,"best_found_improvement":best_native-best,"whole_dna_equal_best_found":dna_matches,"retained_best_assignments":ties.len(),"global_identifiability_certified":false}),
    );
    assert!(
        best < best_native - 1e-6,
        "new non-native score must beat all natives: {}",
        out.display()
    );
    assert!(
        !ties.is_empty() && dna_matches == ties.len(),
        "new complete best DNA differs; inspect {}",
        out.display()
    );
    assert!(
        ledger
            .iter()
            .filter(|r| r["kind"] == "candidate")
            .any(|r| r["guided"] == true),
        "guided neighborhoods must actually be evaluated"
    );
    let candidates: Vec<_> = ledger.iter().filter(|r| r["kind"] == "candidate").collect();
    assert!(
        candidates.iter().any(
            |r| r["delta"]["changes"]
                .as_array()
                .unwrap()
                .iter()
                .any(|c| c["observed"] == 0
                    && c["new_signal"].as_f64().unwrap() > 0.0
                    && c["old_counts"].as_array().unwrap().iter().all(|n| n == 0))
        ),
        "new zero-count realized features must be charged"
    );
    assert!(
        candidates.iter().any(|r| r["public"]["factors"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["positive_residual"] == true)),
        "unsupported positives remain explicit"
    );
    assert!(
        candidates
            .iter()
            .any(|r| r["public"]["assignment"]["routes"]
                .as_array()
                .unwrap()
                .iter()
                .any(
                    |r| r["segments"].as_array().unwrap().iter().any(|s| s["end"]
                        .as_u64()
                        .unwrap()
                        - s["start"].as_u64().unwrap()
                        < 150)
                )),
        "short/terminal pieces must be replayed"
    );
    let mut different_baselines_same_assignment = false;
    for (i, a) in candidates.iter().enumerate() {
        for b in &candidates[i + 1..] {
            if a["public"]["assignment"] == b["public"]["assignment"]
                && a["baseline"] != b["baseline"]
            {
                assert!(
                    (a["public"]["relative_objective"].as_f64().unwrap()
                        - b["public"]["relative_objective"].as_f64().unwrap())
                    .abs()
                        < 1e-9
                );
                if a["baseline_objective"] != b["baseline_objective"] {
                    assert_ne!(a["delta"]["exact_delta"], b["delta"]["exact_delta"]);
                    different_baselines_same_assignment = true;
                }
            }
        }
    }
    assert!(
        different_baselines_same_assignment,
        "changed baselines must recompute deltas despite profile reuse"
    );
    assert!(result["meter"]["route_cache_hits"].as_u64().unwrap() > 0);
}
#[test]
fn automatic_public_mosaic_and_fixed_baseline_feedback() {
    let t = tempfile::tempdir().unwrap();
    let retained = std::env::var_os("IMPG_TEST_RESIDUAL_OUTPUT").map(PathBuf::from);
    let root = retained.as_deref().unwrap_or(t.path()).join("mosaic");
    let f = setup(&root, false);
    search::run(options(&f, "automatic"), None).unwrap();
    assess(&f, "automatic");
    let id = genome::PanelIdentity::read(p(&f.root.join("panel.syng"))).unwrap();
    let panel = SyngIndex::load(p(&f.root.join("panel.syng")), Default::default()).unwrap();
    let mut last = f.sources[0].1.clone();
    last[1216..1240].copy_from_slice(&f.sources[1].1[1216..1240]);
    make_sample(
        &f.root,
        "last-sample",
        &[last, f.truth[1].clone()],
        &panel,
        &id,
    );
    let mut choices = Vec::new();
    for name in ["mechanism-middle", "mechanism-last"] {
        let mut o = options(&f, name);
        o.max_epochs = 0;
        if name.ends_with("last") {
            o.sample = f.root.join("last-sample.membwt");
        }
        search::run(o, Some(0)).unwrap();
        let events = lines(&f.root.join(name).join("events.jsonl"));
        let selected: Vec<_> = events
            .iter()
            .filter(|r| r["event"] == "refinement_choice" && r["batch"] == 8)
            .map(|r| r["children"].clone())
            .collect();
        assert!(!selected.is_empty());
        let attempted: Vec<_> = events
            .iter()
            .filter(|r| r["event"] == "geometry_attempt" && r["guided"] == true)
            .map(|r| r["region"].clone())
            .collect();
        assert!(!attempted.is_empty());
        choices.push((selected, attempted));
    }
    write(
        &f.root.join("mechanism-comparison.json"),
        &json!({"middle":choices[0],"last":choices[1],"baseline_family":0,"fixed_epoch":0}),
    );
    assert_ne!(
        choices[0].0, choices[1].0,
        "evidence must change selected refinement geometry"
    );
    assert_ne!(
        choices[0].1, choices[1].1,
        "evidence must change actual subsequent geometry attempts"
    );
}
#[test]
fn reciprocal_capacity_and_budget_stops() {
    let t = tempfile::tempdir().unwrap();
    let retained = std::env::var_os("IMPG_TEST_RESIDUAL_OUTPUT").map(PathBuf::from);
    let root = retained.as_deref().unwrap_or(t.path()).join("coupled");
    let f = setup(&root, true);
    let result = search::run(options(&f, "exchange"), None).unwrap();
    assert!(result["capacity_rejections"].as_u64().unwrap() > 0);
    assert!(result["compound_confirmations"].as_u64().unwrap() > 0);
    let rows = lines(&f.root.join("exchange/scores.jsonl"));
    let native = &rows.iter().find(|r| r["kind"] == "native").unwrap()["public"];
    let compound = rows
        .iter()
        .find(|r| r["kind"] == "candidate" && r["compound"] == true)
        .unwrap();
    assert_ne!(native["assignment"], compound["public"]["assignment"]);
    let count = |v: &Value| {
        v["factors"]
            .as_array()
            .unwrap()
            .iter()
            .map(|f| (f["tokens"].clone(), f["counts_by_length"].clone()))
            .collect::<Vec<_>>()
    };
    assert_eq!(count(native), count(&compound["public"]));
    for (name, change) in [
        ("scores", 0),
        ("validations", 1),
        ("profile", 2),
        ("tasks", 3),
        ("state", 4),
        ("work", 5),
        ("segments", 6),
        ("ties", 7),
        ("features", 8),
    ] {
        let mut o = options(&f, &format!("stop-{name}"));
        match change {
            0 => o.max_scores = 0,
            1 => o.max_validations = 0,
            2 => o.max_profile_work = 0,
            3 => o.max_tasks = 1,
            4 => o.max_state_bytes = 1,
            5 => o.max_work = 0,
            6 => o.max_segments = 0,
            7 => o.max_ties = 1,
            _ => o.max_features = 1,
        };
        let out = o.out_dir.clone();
        let r = search::run(o, None);
        let status = read(&out.join("status.json"));
        assert!(status["global_bound"].is_null());
        match r {
            Ok(r) => assert!(
                r["status"].as_str().unwrap().starts_with("budget:"),
                "{name}: {r}"
            ),
            Err(e) => assert!(e.to_string().starts_with("budget:"), "{name}: {e}"),
        }
    }
}
#[test]
fn reverse_crop_and_reciprocal_geometry() {
    let r = routes::Route {
        segments: vec![
            routes::Segment {
                source: 0,
                start: 10,
                end: 50,
                reverse: true,
            },
            routes::Segment {
                source: 1,
                start: 0,
                end: 20,
                reverse: false,
            },
        ],
    };
    assert_eq!(
        search::geometry::crop(&r, 5, 45),
        vec![
            routes::Segment {
                source: 0,
                start: 10,
                end: 45,
                reverse: true
            },
            routes::Segment {
                source: 1,
                start: 0,
                end: 5,
                reverse: false
            }
        ]
    );
    let region = search::geometry::Region {
        slot: 0,
        left: 20,
        right: 40,
        level: 2,
    };
    assert!(!region.refinements(60, 8).is_empty());
    assert!(region.refinements(60, 2).is_empty());
}
#[test]
fn public_reverse_shared_counts_and_unready_inventory_exploration() {
    let t = tempfile::tempdir().unwrap();
    let retained = std::env::var_os("IMPG_TEST_RESIDUAL_OUTPUT").map(PathBuf::from);
    let root = retained.as_deref().unwrap_or(t.path());
    let f = setup_variant(&root.join("reverse"), true, true, false);
    let mut o = options(&f, "exchange");
    o.max_scores = 64;
    let result = search::run(o, None).unwrap();
    assert!(
        result["unsupported_opposite_strand_reciprocals"]
            .as_u64()
            .unwrap()
            > 0
    );
    let identity = genome::PanelIdentity::read(p(&f.root.join("panel.syng"))).unwrap();
    let panel = SyngIndex::load(p(&f.root.join("panel.syng")), Default::default()).unwrap();
    let graph = routes::Graph::load(&f.root.join("routes"), &identity).unwrap();
    let sample = sample::SampleIndex::load(&f.root.join("sample.membwt"), &identity).unwrap();
    assert_eq!(graph.k % 2, 1);
    let base = graph.native_assignment(0).unwrap();
    let mut evaluator = routes::Evaluator::new(
        &f.root.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        50000,
        50000,
    )
    .unwrap();
    let mut meter = search::Meter::new(&options(&f, "replay-accounting"));
    fn patch(
        e: &mut routes::Evaluator<'_>,
        base: &routes::Assignment,
        left: u64,
        right: u64,
        m: &mut search::Meter,
    ) -> (u64, u64, routes::Segment, Value) {
        let (lo, p) = search::geometry::nearest(e, &base.routes[0], left, m)
            .unwrap()
            .unwrap();
        let (hi, q) = search::geometry::nearest(e, &base.routes[0], right, m)
            .unwrap()
            .unwrap();
        let l = search::geometry::bucket(e, &p.word, m).unwrap();
        let r = search::geometry::bucket(e, &q.word, m).unwrap();
        for i in l.0..l.1 {
            for j in r.0..r.1 {
                let a = search::geometry::read_port(&mut e.ports.global, e.graph.k, i, m).unwrap();
                let b = search::geometry::read_port(&mut e.ports.global, e.graph.k, j, m).unwrap();
                if a.source == 1 && b.source == 1 && a.reverse && b.reverse && a.cut > b.cut {
                    return (
                        lo,
                        hi,
                        routes::Segment {
                            source: 1,
                            start: b.cut,
                            end: a.cut,
                            reverse: true,
                        },
                        json!({"requested":[left,right],"snapped":[lo,hi],"entry":p,"exit":q,"donor_entry":a,"donor_exit":b,"hub_indices":[i,j]}),
                    );
                }
            }
        }
        panic!("RC-derived donor port pair missing")
    }
    fn old_swap(
        base: &routes::Assignment,
        lo: u64,
        hi: u64,
        donor: &routes::Segment,
    ) -> routes::Assignment {
        let mut a = base.clone();
        a.routes[0] =
            search::geometry::splice(&base.routes[0], lo, hi, vec![donor.clone()]).unwrap();
        let mut removed = search::geometry::crop(&base.routes[0], lo, hi);
        removed.reverse();
        for p in &mut removed {
            p.reverse = !p.reverse;
        }
        a.routes[1] =
            search::geometry::splice(&base.routes[1], donor.start, donor.end, removed).unwrap();
        a
    }
    // NEW causal replay: original root port/hub enumeration and queued successor
    // path. The historical failed object was not serialized and is not recovered.
    let length = base.routes[0].length().unwrap();
    let (lo, hi, donor, parent) = patch(&mut evaluator, &base, 0, length, &mut meter);
    let mut queue = std::collections::VecDeque::from([0, 1, 2]);
    let mut replay = Vec::new();
    let mut reciprocal_error = None;
    while let Some(task) = queue.pop_front() {
        if task == 1 {
            replay.push(json!({"task":task,"parent":"root-port-pair","outcome":"same-slot successor skipped"}));
            continue;
        }
        let mut a = base.clone();
        a.routes[0] =
            search::geometry::splice(&base.routes[0], lo, hi, vec![donor.clone()]).unwrap();
        if task == 2 {
            a = old_swap(&base, lo, hi, &donor);
        }
        let outcome = evaluator.validate_assignment(&a);
        let error = outcome.as_ref().err().map(|e| e.to_string());
        replay.push(json!({"task":task,"parent":"root-port-pair","other_slot":if task==2{Some(1)}else{None},"assignment":a,"error":error}));
        if task == 2 {
            reciprocal_error = error;
        }
    }
    write(
        &f.root.join("queued-opposite-strand-replay.json"),
        &json!({"provenance":"new causal replay, NOT recovered historical assignment","graph":graph.digest().unwrap(),"k":graph.k,"forward_cut_offset":graph.cut_offset,"reverse_cut_offset":graph.k-graph.cut_offset,"parent_task":parent,"queued_tasks":replay}),
    );
    assert_eq!(
        reciprocal_error.as_deref(),
        Some("undeclared outgoing switch")
    );
    assert!(matches!(
        search::geometry::reciprocal(&base, 0, lo, hi, &donor, 1).unwrap(),
        search::geometry::Reciprocal::UnsupportedOppositeStrand
    ));
    // Keep the coincidentally valid opposite-strand example: unsupported is NOT
    // a declaration of infeasibility under the unchanged public model.
    let (lo, hi, donor, parent) = patch(
        &mut evaluator,
        &base,
        length / 4,
        3 * length / 4,
        &mut meter,
    );
    let valid = old_swap(&base, lo, hi, &donor);
    evaluator.validate_assignment(&valid).unwrap();
    assert!(matches!(
        search::geometry::reciprocal(&base, 0, lo, hi, &donor, 1).unwrap(),
        search::geometry::Reciprocal::UnsupportedOppositeStrand
    ));
    write(
        &f.root.join("valid-opposite-strand-deferred.json"),
        &json!({"parent_task":parent,"assignment":valid,"public_feasible":true,"b1_constructor_support":false,"k":graph.k}),
    );
    // Additional positive case, NOT a replacement for the failed paired case:
    // RC-derived donor belongs to a separate inventory and is unoccupied.
    let f = setup_variant(&root.join("reverse-single"), true, true, true);
    let mut o = options(&f, "patches");
    o.max_scores = 64;
    o.max_epochs = 0;
    search::run(o, None).unwrap();
    let rows = lines(&f.root.join("patches/scores.jsonl"));
    let candidates: Vec<_> = rows.iter().filter(|r| r["kind"] == "candidate").collect();
    assert!(!candidates.is_empty());
    assert!(candidates
        .iter()
        .any(|r| r["public"]["assignment"]["routes"]
            .as_array()
            .unwrap()
            .iter()
            .any(|r| r["segments"]
                .as_array()
                .unwrap()
                .iter()
                .any(|s| s["reverse"] == true))));
    for r in candidates {
        assert!(r["delta"]["exact_delta"].as_f64().unwrap().abs() < 1e-9);
    }
    let f = setup_variant(&root.join("readiness"), true, false, true);
    let mut o = options(&f, "explore");
    o.max_scores = 16;
    o.max_epochs = 0;
    search::run(o, None).unwrap();
    let events = lines(&f.root.join("explore/events.jsonl"));
    let first = events
        .iter()
        .position(|r| r["event"] == "geometry_attempt")
        .unwrap();
    let initialized = events[..first]
        .iter()
        .filter(|r| r["event"] == "baseline")
        .count();
    assert!(
        initialized < 3 && initialized > 0,
        "ready exploration must precede all inventories being initialized"
    );
    let mut o = options(&f, "coarse-only");
    o.max_epochs = 0;
    o.max_level = 0;
    o.max_scores = 64;
    let result = search::run(o, None).unwrap();
    assert_eq!(result["refinement_tasks_generated"], 0);
    assert!(result["global_bound"].is_null());
    assert_eq!(result["correlated_support_complete"], false);
}
fn completed_score_records(out: &Path) -> usize {
    lines(&out.join("scores.jsonl"))
        .iter()
        .filter(|r| {
            matches!(
                r["kind"].as_str(),
                Some("native" | "public_confirmation" | "unconfirmed_exact_delta")
            )
        })
        .count()
}
#[test]
fn boundary_public_feature_caps_retain_native_and_candidate_work() {
    let t = tempfile::tempdir().unwrap();
    let retained = std::env::var_os("IMPG_TEST_RESIDUAL_BOUNDARY_OUTPUT").map(PathBuf::from);
    let root = retained.as_deref().unwrap_or(t.path()).join("public-caps");
    let f = setup(&root, false);
    for cap in [62, 84] {
        let mut o = options(&f, &format!("cap-{cap}"));
        o.max_features = cap;
        let out = o.out_dir.clone();
        let result = search::run(o, None);
        write(
            &out.join("boundary-observation.json"),
            &json!({"cap":cap,"error":result.as_ref().err().map(|e|e.to_string()),"has_result":out.join("result.json").exists(),"has_pending":out.join("pending.json").exists(),"completed_score_records":completed_score_records(&out)}),
        );
        let result =
            result.expect("public cap must finalize a bounded result, not generic failure");
        assert!(result["status"]
            .as_str()
            .unwrap()
            .starts_with("budget: public_feature_terms:"));
        let pending = read(&out.join("pending.json"));
        assert_eq!(
            result["meter"]["exact_scores"].as_u64().unwrap() as usize,
            completed_score_records(&out)
        );
        let meter = &result["meter"];
        assert_eq!(
            meter["score_attempts"].as_u64().unwrap(),
            meter["exact_scores"].as_u64().unwrap() + 1
        );
        assert_eq!(
            meter["public_score_attempts"].as_u64().unwrap(),
            meter["public_scores"].as_u64().unwrap() + 1
        );
        if cap == 62 {
            assert_eq!(meter["native_score_attempts"], 1);
            assert_eq!(meter["native_scores"], 0);
            assert_eq!(meter["completed_route_profiles_observed"], 0);
            assert_eq!(meter["public_evaluation_route_reservations"], 2);
            assert_eq!(result["native_initializations"], 0);
            assert_eq!(pending["active_native_family"], 0);
            assert_eq!(pending["pending_native_families"], json!([0, 1]));
        } else {
            assert_eq!(result["native_initializations"], 2);
            assert!(!pending["active_at_stop"].is_null());
            assert!(!pending["queued"].as_array().unwrap().is_empty());
            assert_eq!(pending["pending_native_families"], json!([]));
        }
        assert!(result["global_bound"].is_null());
    }
}
#[test]
fn boundary_only_specific_public_feature_errors_become_stops() {
    use std::io::{Error, ErrorKind};
    for message in [
        "native feature resource cap exhausted (no truncation)",
        "realized feature union exceeds resource cap; no partial objective",
    ] {
        let mapped = search::public_feature_limit(Error::new(ErrorKind::InvalidData, message));
        assert_eq!(
            mapped.to_string(),
            format!("budget: public_feature_terms: {message}")
        );
        let other_kind = search::public_feature_limit(Error::new(ErrorKind::Other, message));
        assert_eq!(other_kind.kind(), ErrorKind::Other);
        assert_eq!(other_kind.to_string(), message);
    }
    for message in [
        "undeclared outgoing switch",
        "switch DNA provenance mismatch",
        "nonfinite sparse full-change loss",
        "realized feature union exceeds resource cap; no partial objective: unrelated suffix",
    ] {
        let unchanged = search::public_feature_limit(Error::new(ErrorKind::InvalidData, message));
        assert_eq!(unchanged.kind(), ErrorKind::InvalidData);
        assert_eq!(unchanged.to_string(), message);
    }
}
#[test]
fn boundary_late_touched_work_never_counts_an_unfinished_delta() {
    let t = tempfile::tempdir().unwrap();
    let retained = std::env::var_os("IMPG_TEST_RESIDUAL_BOUNDARY_OUTPUT").map(PathBuf::from);
    let root = retained.as_deref().unwrap_or(t.path()).join("late-work");
    let f = setup(&root, false);
    // Locate the first candidate admission using only the two native scores.
    // No cap search or frozen recovery budget/recipe changes.
    let mut o = options(&f, "native-only");
    o.max_scores = 2;
    let stopped = search::run(o, None).unwrap();
    assert_eq!(stopped["status"], "budget: exact_scores");
    let pending = read(&f.root.join("native-only/pending.json"));
    let task = &pending["active_at_stop"]["Candidate"];
    assert_eq!(task["other"], 0);
    let base_id = task["base"].as_u64().unwrap() as usize;
    let base: routes::Assignment =
        serde_json::from_value(pending["immutable_baselines"][base_id]["assignment"].clone())
            .unwrap();
    let mut candidate = base.clone();
    let slot = task["region"]["slot"].as_u64().unwrap() as usize;
    candidate.routes[slot] = search::geometry::splice(
        &base.routes[slot],
        task["lo"].as_u64().unwrap(),
        task["hi"].as_u64().unwrap(),
        vec![serde_json::from_value(task["donor"].clone()).unwrap()],
    )
    .unwrap();
    let id = genome::PanelIdentity::read(p(&f.root.join("panel.syng"))).unwrap();
    let panel = SyngIndex::load(p(&f.root.join("panel.syng")), Default::default()).unwrap();
    let graph = routes::Graph::load(&f.root.join("routes"), &id).unwrap();
    let sample = sample::SampleIndex::load(&f.root.join("sample.membwt"), &id).unwrap();
    let mut evaluator = routes::Evaluator::new(
        &f.root.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        50000,
        50000,
    )
    .unwrap();
    evaluator.validate_assignment(&candidate).unwrap();
    let mut integration = 0u64;
    let mut touched = std::collections::BTreeSet::new();
    for (old, new) in base.routes.iter().zip(&candidate.routes) {
        if old != new {
            for route in [old, new] {
                let (q, _) = evaluator.route_counts(route).unwrap();
                for row in q {
                    integration += row.len() as u64;
                    touched.extend(row.keys().copied());
                }
            }
        }
    }
    assert!(touched.len() > 1);
    let mut validation = search::Meter::new(&options(&f, "validation-cost"));
    validation
        .validation(&candidate, graph.k, graph.port_count)
        .unwrap();
    let before = stopped["meter"]["work_used"].as_u64().unwrap();
    let cap = before + validation.work_used + integration + touched.len() as u64 - 1;
    write(
        &f.root.join("derived-late-cap.json"),
        &json!({"candidate_task":task,"work_before_validation":before,"validation":validation.work_used,"profile_integration_terms":integration,"touched_terms":touched.len(),"max_work":cap,"fails_at":"last touched feature in v1; whole touched pass admission in corrected code"}),
    );
    let mut o = options(&f, "late-boundary");
    o.max_work = cap;
    let result = search::run(o, None).unwrap();
    assert_eq!(result["status"], "budget: work");
    let pending = read(&f.root.join("late-boundary/pending.json"));
    assert_eq!(pending["active_at_stop"]["Candidate"], *task);
    let recorded = completed_score_records(&f.root.join("late-boundary"));
    write(
        &f.root.join("late-boundary/boundary-observation.json"),
        &json!({"meter":result["meter"],"completed_score_records":recorded,"active_retained":true}),
    );
    assert_eq!(
        result["meter"]["exact_scores"].as_u64().unwrap() as usize,
        recorded,
        "unfinished delta must not count as a completed score"
    );
    assert_eq!(result["meter"]["exact_delta_scores"], 0);
    assert_eq!(result["meter"]["score_attempts"], 2);
    assert_eq!(result["meter"]["explicit_profile_attempts"], 2);
    assert_eq!(result["meter"]["completed_route_profiles_observed"], 6);
    assert_eq!(
        result["meter"]["work_used"].as_u64().unwrap(),
        before + validation.work_used + integration
    );
}
#[test]
#[ignore = "requires IMPG_TEST_RESIDUAL_BIN executable and fresh IMPG_TEST_RESIDUAL_CLI_OUTPUT"]
fn configured_cli_automatic_mosaic() {
    let root = PathBuf::from(
        std::env::var_os("IMPG_TEST_RESIDUAL_CLI_OUTPUT")
            .expect("fresh configured output required"),
    );
    let f = setup(&root, false);
    let bin = std::env::var_os("IMPG_TEST_RESIDUAL_BIN").expect("standalone executable required");
    let status = Command::new(bin)
        .args([
            "--panel",
            p(&f.root.join("panel.syng")),
            "--routes",
            p(&f.root.join("routes")),
            "--sample",
            p(&f.root.join("sample.membwt")),
            "--out-dir",
            p(&f.root.join("cli")),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    assess(&f, "cli");
}
