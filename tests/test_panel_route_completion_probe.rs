//! Frozen meaningful non-native DNA gate. Truth never enters the proposal harness.
#[path = "../examples/panel_route_completion_probe/mod.rs"]
mod probe;
use clap::Parser;
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::{SyncmerParams, SyngIndex},
};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{BufRead, BufReader, Write},
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
fn p(p: &Path) -> &str {
    p.to_str().unwrap()
}
fn write(p: &Path, v: &Value) {
    fs::write(p, serde_json::to_vec_pretty(v).unwrap()).unwrap();
}
fn read(p: &Path) -> Value {
    genome::read_json(p).unwrap()
}
fn run(bin: &str, args: &[&str], log: &Path, ok: bool) {
    let mut f = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(log)
        .unwrap();
    writeln!(f, "{bin} {args:?}").unwrap();
    let o = Command::new(bin)
        .args(args)
        .stdout(f.try_clone().unwrap())
        .stderr(f)
        .status()
        .unwrap();
    assert!(
        fs::metadata(log).unwrap().len() < 64 * 1024 * 1024,
        "native log alert"
    );
    assert_eq!(o.success(), ok, "see {}", log.display());
}
fn fasta(path: &Path, seq: &[(String, Vec<u8>)]) {
    fs::write(
        path,
        seq.iter()
            .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
}
struct Fixture {
    root: PathBuf,
    panel: PathBuf,
    graph: PathBuf,
    sample: PathBuf,
    sources: Vec<(String, Vec<u8>)>,
}
fn fixture(root: &Path, sources: Vec<(String, Vec<u8>)>, truth: &[Vec<u8>]) -> Fixture {
    fs::create_dir(root).unwrap();
    let mut panel = SyngIndex::build(SyncmerParams::default(), sources.clone().into_iter());
    let prefix = root.join("panel.syng");
    panel.save(p(&prefix)).unwrap();
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    let fa = root.join("sources.fa");
    fasta(&fa, &sources);
    let groups = root.join("groups.json");
    write(
        &groups,
        &json!({"version":1,"groups":sources.iter().enumerate().map(|(i,(n,s))|json!({"id":format!("source-{i}"),"scaffold":null,"occurrences":[{"path":n,"start":0,"end":s.len(),"strand":null}]})).collect::<Vec<_>>()}),
    );
    let catalog = root.join("catalog");
    let bin = env!("CARGO_BIN_EXE_impg");
    run(
        bin,
        &[
            "genome-infer",
            "build-catalog",
            "--panel",
            p(&prefix),
            "--groups",
            p(&groups),
            "--out-dir",
            p(&catalog),
        ],
        &root.join("catalog.log"),
        true,
    );
    let graph = root.join("routes");
    run(
        bin,
        &[
            "genome-infer",
            "build-panel-routes",
            "--panel",
            p(&prefix),
            "--catalog",
            p(&catalog.join("catalog.json")),
            "--sources",
            p(&fa),
            "--read-lengths",
            "150",
            "--core-bp",
            "65536",
            "--out-dir",
            p(&graph),
        ],
        &root.join("graph.log"),
        true,
    );
    let mut reads = vec![];
    for (m, s) in truth.iter().enumerate() {
        let last = s.len() - 150;
        let mut starts = (0..=last).step_by(15).collect::<Vec<_>>();
        if starts.last() != Some(&last) {
            starts.push(last);
        }
        for start in starts {
            reads.push((format!("m{m}-s{start}"), s[start..start + 150].to_vec()));
        }
    }
    let reads_path = root.join("reads.fa");
    fasta(&reads_path, &reads);
    let sample_path = root.join("sample.membwt");
    sample::build(&panel, identity, &[reads_path])
        .unwrap()
        .save(&sample_path)
        .unwrap();
    Fixture {
        root: root.into(),
        panel: prefix,
        graph,
        sample: sample_path,
        sources,
    }
}
fn opts(f: &Fixture, snapshot: &Path, out: &Path) -> probe::Options {
    probe::Options::parse_from([
        "probe",
        "--panel",
        p(&f.panel),
        "--routes",
        p(&f.graph),
        "--sample",
        p(&f.sample),
        "--snapshot",
        p(snapshot),
        "--out-dir",
        p(out),
    ])
}
fn guided(f: &Fixture, bin: &str, out: &Path, extra: &[&str], ok: bool) {
    let mut args = vec![
        "genome-infer",
        "search-panel-routes-guided",
        "--panel",
        p(&f.panel),
        "--routes",
        p(&f.graph),
        "--sample",
        p(&f.sample),
        "--haploid-depth",
        "10",
        "--background",
        "0.1",
        "--max-work",
        "100000",
        "--max-evaluations",
        "10000",
        "--max-state-bytes",
        "268435456",
        "--max-optima",
        "20000",
        "--out-dir",
        p(out),
    ];
    args.extend(extra);
    run(bin, &args, &out.with_extension("log"), ok);
}
fn spell(a: &routes::Assignment, g: &routes::Graph, sources: &[(String, Vec<u8>)]) -> Vec<Vec<u8>> {
    a.routes
        .iter()
        .map(|r| {
            r.segments
                .iter()
                .flat_map(|s| {
                    let name = &g.lanes[s.source].name;
                    let seq = &sources.iter().find(|(n, _)| n == name).unwrap().1;
                    let v = seq[s.start as usize..s.end as usize].to_vec();
                    if s.reverse {
                        impg::graph::reverse_complement(&v)
                    } else {
                        v
                    }
                })
                .collect()
        })
        .collect()
}
fn whole_equal(a: &[Vec<u8>], b: &[Vec<u8>]) -> bool {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    a.sort();
    b.sort();
    a == b
}
/// Independent of route native profiles/replay: explicitly spell all molecules,
/// enumerate ALL L150 windows, and build a separate public MEM-BWT count index.
/// Only the fixed public read/count operator is shared, not evaluator counts/loss.
fn independent(
    f: &Fixture,
    g: &routes::Graph,
    panel: &SyngIndex,
    a: &routes::Assignment,
    e: &mut routes::Evaluator<'_>,
    label: &str,
) -> Value {
    let seq = spell(a, g, &f.sources);
    let windows = seq
        .iter()
        .enumerate()
        .flat_map(|(m, s)| {
            s.windows(150)
                .enumerate()
                .map(move |(i, w)| (format!("{m}-{i}"), w.to_vec()))
        })
        .collect::<Vec<_>>();
    let windows_path = f.root.join(format!("oracle-{label}.fa"));
    fasta(&windows_path, &windows);
    let id = genome::PanelIdentity::read(p(&f.panel)).unwrap();
    let q = sample::build(panel, id.clone(), &[windows_path])
        .unwrap()
        .counts
        .observed_pairs()
        .unwrap();
    let observed = sample::SampleIndex::load(&f.sample, &id)
        .unwrap()
        .counts
        .observed_pairs()
        .unwrap();
    let keys = q
        .keys()
        .chain(observed.keys())
        .copied()
        .collect::<BTreeSet<_>>();
    let mut loss = 0.0;
    let mut rows = vec![];
    let mut unsupported = 0;
    let mut realized_zero = 0;
    for key in keys {
        let count = q.get(&key).copied().unwrap_or(0);
        let observed = observed.get(&key).copied().unwrap_or(0);
        let signal = 10.0 * (count as f64 / 150.0);
        let term = signal - observed as f64 * (1.0 + signal / 0.1).ln();
        loss += term;
        unsupported += usize::from(count == 0 && observed > 0);
        realized_zero += usize::from(count > 0 && observed == 0);
        rows.push((key, count, observed));
    }
    let result = e.evaluate(a).unwrap();
    assert_eq!(result.factors.len(), rows.len());
    for (factor, (key, count, observed)) in result.factors.iter().zip(&rows) {
        assert_eq!(factor.tokens, *key);
        assert_eq!(factor.counts_by_length, vec![*count]);
        assert_eq!(factor.observed, *observed);
    }
    assert!(
        (loss - result.relative_objective).abs() <= 1e-8,
        "independent objective mismatch"
    );
    let public_routes = a
        .routes
        .iter()
        .map(|r| e.route_counts(r).unwrap().0[0].clone())
        .collect::<Vec<_>>();
    let mut aggregate = BTreeMap::new();
    for c in public_routes {
        for (k, n) in c {
            let v = aggregate.entry(k).or_insert(0u64);
            *v = v.checked_add(n).unwrap();
        }
    }
    assert_eq!(aggregate, q);
    json!({"method":"all spelled L150 windows -> independent public sample MEM-BWT; separate union/loss loop ln(1+s/beta)","independent_objective":loss,"evaluator_objective_bits":result.relative_objective.to_bits(),"per_length_exact_equal":true,"unsupported_positive_features":unsupported,"realized_zero_observation_features":realized_zero,"counts":rows,"admitted_starts":result.admitted_starts})
}
#[test]
#[ignore = "frozen meaningful mosaic gate; requires fresh IMPG_TEST_COMPLETION_OUTPUT and preserved IMPG_TEST_COMPLETION_V3"]
fn completion_probe_frozen_non_native_mosaic() {
    let root = PathBuf::from(
        std::env::var("IMPG_TEST_COMPLETION_OUTPUT").expect("fresh gate output required"),
    );
    let bin = std::env::var("IMPG_TEST_COMPLETION_V3").expect("frozen v3 binary required");
    assert!(Path::new(&bin).is_absolute() && Path::new(&bin).is_file());
    let a = dna(1200, 113);
    let mut b = a.clone();
    for start in [240, 600, 960] {
        for base in &mut b[start..start + 24] {
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
    }
    let mut mosaic = a.clone();
    mosaic[600..624].copy_from_slice(&b[600..624]);
    let second = dna(400, 57);
    let truth = vec![mosaic.clone(), second.clone()];
    let f = fixture(
        &root,
        vec![
            ("A#0#one".into(), a),
            ("B#0#one".into(), b),
            ("A#0#two".into(), second),
            ("B#0#two".into(), dna(400, 58)),
        ],
        &truth,
    );
    write(
        &root.join("truth-for-assessment-only.json"),
        &json!({"molecules":truth.iter().map(|s|String::from_utf8_lossy(s)).collect::<Vec<_>>(),"total_bp":1600}),
    );
    let snapshot = root.join("ordinary-snapshot");
    guided(&f, &bin, &snapshot, &[], true);
    let checkpoint_before = fs::read(snapshot.join("checkpoint.json")).unwrap();
    let ledger_before = fs::read(snapshot.join("evaluations.jsonl")).unwrap();
    let out = root.join("probe");
    probe::run(opts(&f, &snapshot, &out)).unwrap();
    // Durable output freeze precedes any spelling/selection of emitted assignments.
    let freeze = Command::new("sha256sum")
        .args([
            out.join("proposals.jsonl"),
            out.join("result.json"),
            out.join("provenance.json"),
        ])
        .output()
        .unwrap();
    assert!(freeze.status.success());
    fs::write(
        root.join("outputs-frozen-before-assessment.sha256"),
        freeze.stdout,
    )
    .unwrap();
    let id = genome::PanelIdentity::read(p(&f.panel)).unwrap();
    let g = routes::Graph::load(&f.graph, &id).unwrap();
    let panel = SyngIndex::load(p(&f.panel), SyncmerParams::default()).unwrap();
    let sample = sample::SampleIndex::load(&f.sample, &id).unwrap();
    let mut e = routes::Evaluator::new(&f.graph, &g, &panel, &sample, 10.0, 0.1, 50000000, 1000000)
        .unwrap();
    let cp: Value = serde_json::from_slice(&checkpoint_before).unwrap();
    let natives = cp["state"]["native"].as_array().unwrap();
    for n in natives {
        let a: routes::Assignment = serde_json::from_value(n["assignment"].clone()).unwrap();
        assert!(
            !whole_equal(&spell(&a, &g, &f.sources), &truth),
            "truth is native"
        );
    }
    let native_best = natives
        .iter()
        .map(|n| f64::from_bits(n["objective_bits"].as_u64().unwrap()))
        .min_by(f64::total_cmp)
        .unwrap();
    let mut previous_best = f64::INFINITY;
    let mut previous_truth = false;
    let mut previous_nonnative = 0;
    for line in ledger_before
        .split(|b| *b == b'\n')
        .filter(|l| !l.is_empty())
    {
        let row: Value = serde_json::from_slice(line).unwrap();
        if row["cross_source_route"] != true {
            continue;
        }
        let a: routes::Assignment =
            serde_json::from_value(row["scored"]["assignment"].clone()).unwrap();
        previous_nonnative += 1;
        previous_best = previous_best.min(row["scored"]["relative_objective"].as_f64().unwrap());
        previous_truth |= whole_equal(&spell(&a, &g, &f.sources), &truth);
    }
    let mut strict = 0;
    let mut exact = 0;
    let mut first = None;
    let mut minimum = f64::INFINITY;
    for line in BufReader::new(fs::File::open(out.join("proposals.jsonl")).unwrap()).lines() {
        let row: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let score = f64::from_bits(row["objective_bits"].as_u64().unwrap());
        minimum = minimum.min(score);
        if score < native_best {
            strict += 1;
            let a: routes::Assignment = serde_json::from_value(row["assignment"].clone()).unwrap();
            let spelled = spell(&a, &g, &f.sources);
            if whole_equal(&spelled, &truth) {
                exact += 1;
                if first.is_none() {
                    first = Some((a, row, spelled));
                }
            }
        }
    }
    write(
        &root.join("assessment-initial.json"),
        &json!({"native_best":native_best,"harness_best":minimum,"strict_better_visits":strict,"whole_truth_strict_visits":exact,"ordinary_prior_nonnative_visits":previous_nonnative,"ordinary_prior_best":if previous_best.is_finite(){Some(previous_best)}else{None},"ordinary_already_spelled_truth":previous_truth,"novel_improvement_over_ordinary_search_claimed":false}),
    );
    assert!(
        strict > 0,
        "no automatic proposal strictly beats EVERY native"
    );
    assert!(exact > 0, "no strict proposal has whole truth DNA");
    let (a, row, spelled) = first.unwrap();
    assert_eq!(e.validate_assignment(&a).unwrap(), a);
    assert_eq!(spelled.len(), 2);
    assert_eq!(spelled.iter().map(Vec::len).sum::<usize>(), 1600);
    write(
        &root.join("whole-dna-evidence.json"),
        &json!({"proposal":row,"spelled_molecules":spelled.iter().map(|s|String::from_utf8_lossy(s)).collect::<Vec<_>>(),"truth_molecules":truth.iter().map(|s|String::from_utf8_lossy(s)).collect::<Vec<_>>(),"exact_whole_multiset_equality":true,"complete_molecules":2,"truth_bp":1600,"query_bp":1600,"matched_truth_bp":1600,"matched_query_bp":1600,"truth_coverage":1.0,"query_coverage":1.0,"synthetic_test_memory_only":true}),
    );
    let mut oracle = vec![independent(&f, &g, &panel, &a, &mut e, "mosaic")];
    for (i, n) in natives.iter().enumerate() {
        let a: routes::Assignment = serde_json::from_value(n["assignment"].clone()).unwrap();
        oracle.push(independent(
            &f,
            &g,
            &panel,
            &a,
            &mut e,
            &format!("native-{i}"),
        ));
    }
    assert!(oracle
        .iter()
        .any(|r| r["unsupported_positive_features"].as_u64().unwrap() > 0));
    assert!(oracle
        .iter()
        .any(|r| r["realized_zero_observation_features"].as_u64().unwrap() > 0));
    write(&root.join("independent-oracle.json"), &json!(oracle));
    let resumed = root.join("unchanged-resume");
    guided(&f, &bin, &resumed, &["--resume-from", p(&snapshot)], true);
    assert_eq!(read(&resumed.join("checkpoint.json"))["state"], cp["state"]);
    assert_eq!(
        fs::read(resumed.join("evaluations.jsonl")).unwrap(),
        ledger_before
    );
    assert_eq!(
        fs::read(snapshot.join("checkpoint.json")).unwrap(),
        checkpoint_before
    );
    assert_eq!(
        fs::read(snapshot.join("evaluations.jsonl")).unwrap(),
        ledger_before
    );
    let mut zero = opts(&f, &snapshot, &root.join("storage-stop"));
    zero.max_state_bytes = 1;
    assert!(probe::run(zero).is_err());
    assert!(read(&root.join("storage-stop/failure.json"))["error"]
        .as_str()
        .unwrap()
        .contains("state-budget-exhausted"));
    integrity_stops(&f, &snapshot);
    write(
        &root.join("gate-passed.json"),
        &json!({"strict_native_improvement":true,"whole_dna":true,"independent_counts_and_objective":true,"unchanged_policy_resume":true,"old_input_unchanged":true}),
    );
}

#[test]
fn completion_probe_dead_tip_multidonor_deferred_and_distractor_service() {
    let temp = tempfile::tempdir().unwrap();
    let root = std::env::var_os("IMPG_TEST_COMPLETION_MECHANISM_OUTPUT")
        .map(PathBuf::from)
        .unwrap_or_else(|| temp.path().join("fixture"));
    let base = dna(400, 113);
    let mut private = base.clone();
    private[300..].copy_from_slice(&dna(100, 479));
    let f = fixture(
        &root,
        vec![
            ("A#0#one".into(), base.clone()),
            ("B#0#one".into(), base.clone()),
            ("C#0#one".into(), base.clone()),
            ("D#0#one".into(), private),
        ],
        &[base],
    );
    let identity = genome::PanelIdentity::read(p(&f.panel)).unwrap();
    let g = routes::Graph::load(&f.graph, &identity).unwrap();
    let panel = SyngIndex::load(p(&f.panel), SyncmerParams::default()).unwrap();
    let sample = sample::SampleIndex::load(&f.sample, &identity).unwrap();
    let mut e = routes::Evaluator::new(&f.graph, &g, &panel, &sample, 10.0, 0.1, 50000000, 1000000)
        .unwrap();
    probe::mechanism_checks(
        opts(&f, Path::new("not-a-checkpoint"), &root.join("mechanism")),
        &mut e,
    )
    .unwrap();
}
fn fnv(bytes: &[u8]) -> String {
    let hash = bytes.iter().fold(0xcbf29ce484222325u64, |h, b| {
        (h ^ *b as u64).wrapping_mul(0x100000001b3)
    });
    format!("fnv1a64:{}:{hash:016x}", bytes.len())
}
fn integrity_stops(f: &Fixture, snapshot: &Path) {
    let original_cp = fs::read(snapshot.join("checkpoint.json")).unwrap();
    let original_ledger = fs::read(snapshot.join("evaluations.jsonl")).unwrap();
    let mut results = vec![];
    for case in [
        "trailing",
        "truncated",
        "policy",
        "sample",
        "depth",
        "native-bits",
        "ledger-eof",
        "seal",
    ] {
        let dir = f.root.join(format!("bad-{case}"));
        fs::create_dir(&dir).unwrap();
        let parsed: Value = serde_json::from_slice(&original_cp).unwrap();
        let original_text = std::str::from_utf8(&original_cp).unwrap();
        // Preserve every unrelated byte, particularly numeric BTreeMap task order.
        // Serializing Value here would sort task IDs lexicographically and confound
        // every intended negative case with an earlier task-order rejection.
        let edit = match case {
            "policy" => Some((
                format!("\"version\": {}", parsed["bindings"]["policy"]["version"]),
                "\"version\": \"wrong\"".into(),
            )),
            "sample" => Some((
                format!("\"sample\": {}", parsed["bindings"]["sample"]),
                "\"sample\": \"wrong\"".into(),
            )),
            "depth" => Some((
                format!("\"depth_bits\": {}", parsed["bindings"]["depth_bits"]),
                format!("\"depth_bits\": {}", 11.0f64.to_bits()),
            )),
            "native-bits" => {
                let n = parsed["state"]["native"][0]["objective_bits"]
                    .as_u64()
                    .unwrap();
                Some((
                    format!("\"objective_bits\": {n}"),
                    format!("\"objective_bits\": {}", n ^ 1),
                ))
            }
            _ => None,
        };
        let mut cp = if let Some((from, to)) = edit {
            assert!(original_text.contains(&from));
            original_text.replacen(&from, &to, 1).into_bytes()
        } else {
            original_cp.clone()
        };
        if case == "trailing" {
            cp.extend(b" garbage");
        }
        if case == "truncated" {
            cp.pop();
        }
        let mut ledger = original_ledger.clone();
        if case == "ledger-eof" {
            ledger.extend(b"{\"visit\":");
        }
        fs::write(dir.join("checkpoint.json"), &cp).unwrap();
        fs::write(dir.join("evaluations.jsonl"), &ledger).unwrap();
        write(
            &dir.join("checkpoint-seal.json"),
            &json!({"version":1,"checkpoint":if case=="seal"{"wrong".into()}else{fnv(&cp)},"ledger":fnv(&ledger)}),
        );
        let mut o = opts(f, &dir, &f.root.join(format!("reject-{case}")));
        o.max_work = 0;
        let err = probe::run(o).expect_err(case).to_string();
        let expected = match case {
            "trailing" => "trailing JSON bytes",
            "truncated" => "premature JSON EOF",
            "policy" | "sample" | "depth" => "snapshot scientific/policy binding mismatch",
            "native-bits" => "native objective bits/ledger mismatch",
            "ledger-eof" => "missing JSON value",
            "seal" => "checkpoint/ledger seal mismatch",
            _ => unreachable!(),
        };
        results.push(json!({"case":case,"error":err,"expected":expected}));
        write(
            &f.root.join("integrity-and-stops-in-progress.json"),
            &json!(results),
        );
        assert_eq!(err, expected, "confounded integrity case {case}");
    }
    let mut zero = opts(f, snapshot, &f.root.join("work-zero"));
    zero.max_work = 0;
    probe::run(zero).unwrap();
    let base = read(&f.root.join("work-zero/result.json"))["logical_reserved_bytes"]
        .as_u64()
        .unwrap();
    for case in ["growth", "work", "evaluation", "distinct", "record"] {
        let out = f.root.join(format!("stop-{case}"));
        let mut o = opts(f, snapshot, &out);
        match case {
            "growth" => o.max_state_bytes = base + 1,
            "work" => o.max_work = 1,
            "evaluation" => o.max_evaluations = 0,
            "distinct" => o.max_distinct = 1,
            "record" => o.record_bytes = 64,
            _ => unreachable!(),
        }
        if case == "record" {
            let err = probe::run(o).unwrap_err().to_string();
            assert!(err.contains("input-record-byte-budget-exhausted"));
            results.push(json!({"case":case,"error":err}));
            continue;
        }
        probe::run(o).unwrap();
        let result = read(&out.join("result.json"));
        let expected = match case {
            "growth" => "state-budget-exhausted",
            "work" => "work-budget-exhausted",
            "evaluation" => "evaluation-budget-exhausted",
            "distinct" => "distinct-budget-exhausted",
            _ => unreachable!(),
        };
        assert_eq!(result["status"], expected);
        assert!(
            result["logical_reserved_bytes"].as_u64().unwrap()
                <= result["logical_limit_bytes"].as_u64().unwrap()
        );
        if case == "growth" {
            assert_eq!(result["counts"]["admitted"], 0);
            assert!(!result["pending_input"].is_null());
        }
        if case == "evaluation" {
            assert_eq!(result["counts"]["fresh_evaluations"], 0);
            assert!(!result["pending_chains"].as_array().unwrap().is_empty());
        }
        results.push(json!({"case":case,"status":result["status"],"pending_input":result["pending_input"],"pending_chains":result["pending_chains"]}));
    }
    assert_eq!(
        fs::read(snapshot.join("checkpoint.json")).unwrap(),
        original_cp
    );
    assert_eq!(
        fs::read(snapshot.join("evaluations.jsonl")).unwrap(),
        original_ledger
    );
    write(&f.root.join("integrity-and-stops.json"), &json!(results));
}

#[test]
fn completion_ready_family_causal_service_and_outcomes() {
    let temp = tempfile::tempdir().unwrap();
    let root = std::env::var_os("IMPG_TEST_READY_FAMILY_OUTPUT")
        .map(PathBuf::from)
        .unwrap_or_else(|| temp.path().join("fixture"));
    let base = dna(400, 113);
    let mut private = base.clone();
    private[300..].copy_from_slice(&dna(100, 479));
    let f = fixture(
        &root,
        vec![
            ("A#0#one".into(), base.clone()),
            ("B#0#one".into(), base.clone()),
            ("C#0#one".into(), base.clone()),
            ("D#0#one".into(), private),
        ],
        &[base],
    );
    let id = genome::PanelIdentity::read(p(&f.panel)).unwrap();
    let g = routes::Graph::load(&f.graph, &id).unwrap();
    let panel = SyngIndex::load(p(&f.panel), SyncmerParams::default()).unwrap();
    let sample = sample::SampleIndex::load(&f.sample, &id).unwrap();
    let mut e = routes::Evaluator::new(&f.graph, &g, &panel, &sample, 10.0, 0.1, 50000000, 1000000)
        .unwrap();
    probe::ready_family_checks(
        opts(
            &f,
            Path::new("mechanism-only-not-a-snapshot"),
            &root.join("service"),
        ),
        &mut e,
    )
    .unwrap();
}
