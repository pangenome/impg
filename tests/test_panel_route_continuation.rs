//! Continuation and source-order synthetic gates; frozen forward recipes remain unchanged.
#[path = "panel_route_source_pair_order/mod.rs"]
mod source_pair_order;
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
fn p(x: &Path) -> &str {
    x.to_str().unwrap()
}
fn read(x: &Path) -> Value {
    genome::read_json(x).unwrap()
}
fn write(x: &Path, v: &Value) {
    fs::write(x, serde_json::to_vec_pretty(v).unwrap()).unwrap();
}
fn run(binary: &str, args: &[&str], ok: bool) -> String {
    let o = Command::new(binary).args(args).output().unwrap();
    let text = format!(
        "{}{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    assert!(
        !text.contains("++match ") && !text.contains("++path "),
        "native debug dump"
    );
    assert!(text.len() <= 1048576, "diagnostic cap");
    assert_eq!(o.status.success(), ok, "{args:?}\n{text}");
    text
}
struct Fixture {
    root: PathBuf,
    panel: PathBuf,
    graph: PathBuf,
    samples: Vec<PathBuf>,
}
fn fixture(root: &Path, repeats: usize) -> Fixture {
    fixture_with_slots(root, repeats, 1)
}
fn fixture_with_slots(root: &Path, repeats: usize, slots: usize) -> Fixture {
    fixture_with_orientation(root, repeats, slots, false)
}
fn fixture_with_orientation(root: &Path, repeats: usize, slots: usize, reverse: bool) -> Fixture {
    fs::create_dir(root).expect("fresh synthetic fixture directory");
    let a = dna(400, 113);
    let mut b = a.clone();
    for i in [150, 260] {
        b[i] = if b[i] == b'A' { b'C' } else { b'A' };
    }
    let donor = if reverse {
        b.iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect::<Vec<_>>()
    } else {
        b.clone()
    };
    let mut seq = vec![
        ("A#0#native".to_string(), a.repeat(repeats)),
        ("B#0#native".to_string(), donor.repeat(repeats)),
    ];
    if slots == 2 {
        seq.extend([
            ("A#0#second".into(), dna(160, 57)),
            ("B#0#second".into(), dna(160, 58)),
        ]);
    }
    let mut panel = SyngIndex::build(SyncmerParams::default(), seq.clone().into_iter());
    let prefix = root.join("panel.syng");
    panel.save(p(&prefix)).unwrap();
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    let fasta = root.join("sources.fa");
    fs::write(
        &fasta,
        seq.iter()
            .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    let groups = root.join("groups.json");
    write(
        &groups,
        &json!({"version":1,"groups":seq.iter().enumerate().map(|(i,(n,s))|json!({"id":format!("source-{i}"),"scaffold":null,"occurrences":[{"path":n,"start":0,"end":s.len(),"strand":null}]})).collect::<Vec<_>>()}),
    );
    let catalog = root.join("catalog");
    let binary = env!("CARGO_BIN_EXE_impg");
    run(
        binary,
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
        true,
    );
    let graph = root.join("routes");
    run(
        binary,
        &[
            "genome-infer",
            "build-panel-routes",
            "--panel",
            p(&prefix),
            "--catalog",
            p(&catalog.join("catalog.json")),
            "--sources",
            p(&fasta),
            "--read-lengths",
            "150",
            "--core-bp",
            "65536",
            "--out-dir",
            p(&graph),
        ],
        true,
    );
    let mut samples = vec![];
    for (label, s) in [("a", a), ("b", b)] {
        let reads = root.join(format!("{label}.fa"));
        fs::write(
            &reads,
            s.windows(150)
                .map(|w| format!(">read\n{}\n", String::from_utf8_lossy(w)))
                .collect::<String>(),
        )
        .unwrap();
        let sample = sample::build(&panel, identity.clone(), &[reads]).unwrap();
        let path = root.join(format!("{label}.membwt"));
        sample.save(&path).unwrap();
        samples.push(path);
    }
    Fixture {
        root: root.into(),
        panel: prefix,
        graph,
        samples,
    }
}
fn guided(
    f: &Fixture,
    binary: &str,
    sample: usize,
    label: &str,
    extra: &[&str],
    ok: bool,
) -> PathBuf {
    let out = f.root.join(label);
    let mut args = vec![
        "genome-infer",
        "search-panel-routes-guided",
        "--panel",
        p(&f.panel),
        "--routes",
        p(&f.graph),
        "--sample",
        p(&f.samples[sample]),
        "--haploid-depth",
        "150",
        "--out-dir",
        p(&out),
    ];
    for (flag, value) in [
        ("--max-work", "20000"),
        ("--max-evaluations", "10000"),
        ("--max-state-bytes", "268435456"),
        ("--max-optima", "20000"),
    ] {
        if !extra.contains(&flag) {
            args.extend([flag, value]);
        }
    }
    args.extend_from_slice(extra);
    let text = run(binary, &args, ok);
    fs::write(f.root.join(format!("{label}.log")), text).unwrap();
    out
}
#[test]
#[ignore = "expensive frozen LONG-PRODUCER recipe; requires fresh IMPG_TEST_CONTINUATION_OUTPUT"]
fn guided_long_producer_frozen_gate() {
    let root = PathBuf::from(
        std::env::var("IMPG_TEST_CONTINUATION_OUTPUT")
            .expect("requested long gate needs IMPG_TEST_CONTINUATION_OUTPUT"),
    );
    assert!(root.is_absolute());
    let f = fixture(&root, 10000);
    let identity = genome::PanelIdentity::read(p(&f.panel)).unwrap();
    let g = routes::Graph::load(&f.graph, &identity).unwrap();
    let prerequisites = json!({"source_port_counts":g.lanes.iter().map(|l|l.port_count).collect::<Vec<_>>(),"virtual_ranges":g.lanes.iter().map(|l|l.port_count.next_power_of_two()).collect::<Vec<_>>(),"global_ports":g.port_count,"families":g.families.len(),"k":g.k});
    write(
        &root.join("prerequisites-before-search.json"),
        &prerequisites,
    );
    assert_eq!(g.families.len(), 2);
    assert!(g
        .lanes
        .iter()
        .all(|l| l.port_count.next_power_of_two() >= 200000));
    assert!(g.port_count > 4);
    let binary = env!("CARGO_BIN_EXE_impg");
    let mut orders = vec![];
    for sample in 0..2 {
        let mut baseline = None;
        for repeat in 0..2 {
            let out = guided(
                &f,
                binary,
                sample,
                &format!("sample-{sample}-{repeat}"),
                &[],
                true,
            );
            let r = read(&out.join("result.json"))["result"].clone();
            let ledger = fs::read(out.join("evaluations.jsonl")).unwrap();
            if let Some(ref expected) = baseline {
                assert_eq!(&ledger, expected);
            } else {
                baseline = Some(ledger.clone());
                orders.push(r["native_ranked_identities"].clone());
            }
            let mut first = std::collections::BTreeMap::new();
            for line in ledger.split(|&b| b == b'\n').filter(|s| !s.is_empty()) {
                let row: Value = serde_json::from_slice(line).unwrap();
                let a: routes::Assignment =
                    serde_json::from_value(row["scored"]["assignment"].clone()).unwrap();
                if row["cross_source_route"] == true {
                    assert!(a.routes.iter().any(|route| route
                        .segments
                        .iter()
                        .any(|s| s.source != route.segments[0].source && s.end > s.start)));
                    first
                        .entry(a.family)
                        .or_insert(row["work"].as_u64().unwrap());
                }
            }
            write(&out.join("first-mixed-work.json"), &json!(first));
            if repeat == 0 {
                for (&family, &work) in &first {
                    let stop = guided(
                        &f,
                        binary,
                        sample,
                        &format!("first-mixed-{sample}-{family}"),
                        &["--max-work", &work.to_string()],
                        true,
                    );
                    let metrics = read(&stop.join("result.json"))["result"].clone();
                    write(
                        &out.join(format!("first-mixed-family-{family}-metrics.json")),
                        &metrics,
                    );
                }
                let mut lo = 2;
                let mut hi = *first.values().min().unwrap();
                while lo < hi {
                    let mid = lo + (hi - lo) / 2;
                    let stop = guided(
                        &f,
                        binary,
                        sample,
                        &format!("donor-boundary-{sample}-{mid}"),
                        &["--max-work", &mid.to_string()],
                        true,
                    );
                    let r = read(&stop.join("result.json"));
                    if r["result"]["donor_transitions_examined"].as_u64().unwrap() > 0 {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }
                write(
                    &out.join("first-donor-examined-work.json"),
                    &json!({"work":lo,"meaning":"hub member examined; not itself a positive biological traversal"}),
                );
            }

            assert_eq!(first.len(), 2, "LONG-PRODUCER fixed gate failed: {r}");
            let cp = read(&out.join("checkpoint.json"));
            for family in 0..2 {
                assert!(
                    cp["state"]["tasks"]
                        .as_object()
                        .unwrap()
                        .values()
                        .any(|t| t["context"]["family"] == family
                            && t["op"].get("SourceScan").is_some()),
                    "source producer exhausted"
                );
            }
            eprintln!("LONG-PRODUCER sample={sample} repeat={repeat} first_mixed={first:?} mixed={} work={} peak_bytes={} peak_tasks={}",r["mixed_assignments_evaluated"],r["work"],r["peak_state_bytes"],r["peak_frontier_tasks"]);
        }
    }
    assert_eq!(orders[0][0], "A#0");
    assert_eq!(orders[1][0], "B#0");
}
fn hash(path: &Path) -> String {
    let b = fs::read(path).unwrap();
    format!("fnv1a64:{}:{:016x}", b.len(), genome::checksum(&b))
}
fn pin(root: &Path, outputs: &[PathBuf]) -> PathBuf {
    let path = root.join(format!("pins-{}.json", outputs.len()));
    write(
        &path,
        &json!({"version":1,"checkpoints":outputs.iter().map(|o|json!({"directory":fs::canonicalize(o).unwrap(),"checkpoint":hash(&o.join("checkpoint.json")),"ledger":hash(&o.join("evaluations.jsonl"))})).collect::<Vec<_>>()}),
    );
    path
}
#[test]
#[ignore = "legacy v1-to-v2 gate; requires IMPG_TEST_GUIDED_V1, IMPG_TEST_GUIDED_V2 and fresh IMPG_TEST_TRANSITION_OUTPUT"]
fn guided_exact_v1_conversion_and_resume() {
    let v2 = std::env::var("IMPG_TEST_GUIDED_V2")
        .expect("legacy conversion gate requires preserved v2 binary");
    let binary = v2.as_str();
    let old = std::env::var("IMPG_TEST_GUIDED_V1")
        .expect("requested transition gate needs IMPG_TEST_GUIDED_V1");
    let check = Command::new("sha256sum").arg(&old).output().unwrap();
    assert!(check.status.success());
    assert_eq!(
        String::from_utf8(check.stdout)
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap(),
        "6722ef903cb075f76a8b0f7d16fca9394952033d3aa80b1185195259f215fffb"
    );
    let root = PathBuf::from(
        std::env::var("IMPG_TEST_TRANSITION_OUTPUT")
            .expect("requested transition gate needs IMPG_TEST_TRANSITION_OUTPUT"),
    );
    assert!(root.is_absolute());
    let f = fixture(&root, 1);
    let mut ancestors: Vec<PathBuf> = vec![];
    for cap in [1, 2, 5, 32, 128, 512, 2048] {
        // Caps use a direct subprocess to avoid duplicate clap options.
        let out = f.root.join(format!("v1-{cap}"));
        let captext = cap.to_string();
        let mut args = vec![
            "genome-infer",
            "search-panel-routes-guided",
            "--panel",
            p(&f.panel),
            "--routes",
            p(&f.graph),
            "--sample",
            p(&f.samples[0]),
            "--haploid-depth",
            "150",
            "--max-work",
            &captext,
            "--max-evaluations",
            "10000",
            "--max-state-bytes",
            "268435456",
            "--max-optima",
            "20000",
            "--out-dir",
            p(&out),
        ];
        if let Some(prior) = ancestors.last() {
            args.extend(["--resume-from", p(prior), "--extend-budgets"]);
        }
        run(&old, &args, true);
        ancestors.push(out.clone());
        let declaration = pin(&root, &ancestors);
        let original = read(&out.join("checkpoint.json"));
        let oldledger = fs::read(out.join("evaluations.jsonl")).unwrap();
        guided(
            &f,
            binary,
            0,
            &format!("ordinary-reject-{cap}"),
            &["--resume-from", p(&out), "--extend-budgets"],
            false,
        );
        let converted = guided(
            &f,
            binary,
            0,
            &format!("converted-{cap}"),
            &[
                "--resume-from",
                p(&out),
                "--extend-budgets",
                "--convert-exact-v1-to-v2",
                "--v1-ancestry-manifest",
                p(&declaration),
            ],
            true,
        );
        let cp = read(&converted.join("checkpoint.json"));
        assert_eq!(cp["state"]["work"], original["state"]["work"]);
        assert_eq!(cp["state"]["evaluations"], original["state"]["evaluations"]);
        assert_eq!(
            fs::read(converted.join("evaluations.jsonl")).unwrap(),
            oldledger
        );
        for key in [
            "native",
            "fifo",
            "shallow",
            "next_id",
            "next_ready",
            "families",
            "incumbent",
            "support",
            "lost_support",
            "family_cursor",
            "quantum",
            "ranking",
        ] {
            assert_eq!(cp["state"][key], original["state"][key], "{key}");
        }
        for (id, t) in original["state"]["tasks"].as_object().unwrap() {
            let mut nt = cp["state"]["tasks"][id].clone();
            assert!(nt["prev"].is_null() && nt["next"].is_null());
            nt.as_object_mut().unwrap().remove("prev");
            nt.as_object_mut().unwrap().remove("next");
            assert_eq!(&nt, t);
        }
        let a = guided(
            &f,
            binary,
            0,
            &format!("resumed-a-{cap}"),
            &["--resume-from", p(&converted)],
            true,
        );
        let b = guided(
            &f,
            binary,
            0,
            &format!("resumed-b-{cap}"),
            &["--resume-from", p(&converted)],
            true,
        );
        assert_eq!(
            fs::read(a.join("evaluations.jsonl")).unwrap(),
            fs::read(b.join("evaluations.jsonl")).unwrap()
        );
        assert_eq!(
            read(&a.join("checkpoint.json"))["state"],
            read(&b.join("checkpoint.json"))["state"]
        );
        let path = out.join("checkpoint.json");
        let sealpath = out.join("checkpoint-seal.json");
        let bytes = fs::read(&path).unwrap();
        let sealbytes = fs::read(&sealpath).unwrap();
        let mut corrupted = original.clone();
        corrupted["state"]["donors"] = 999999.into();
        write(&path, &corrupted);
        let mut seal: Value = serde_json::from_slice(&sealbytes).unwrap();
        seal["checkpoint"] = hash(&path).into();
        write(&sealpath, &seal);
        guided(
            &f,
            binary,
            0,
            &format!("resealed-pin-reject-{cap}"),
            &[
                "--resume-from",
                p(&out),
                "--extend-budgets",
                "--convert-exact-v1-to-v2",
                "--v1-ancestry-manifest",
                p(&declaration),
            ],
            false,
        );
        fs::write(&path, bytes).unwrap();
        fs::write(&sealpath, sealbytes).unwrap();
    }
    let old_root = ancestors.last().unwrap();
    let declaration = root.join(format!("pins-{}.json", ancestors.len()));
    let path = old_root.join("checkpoint.json");
    let sealpath = old_root.join("checkpoint-seal.json");
    let original_bytes = fs::read(&path).unwrap();
    let original_seal = fs::read(&sealpath).unwrap();
    let original: Value = serde_json::from_slice(&original_bytes).unwrap();
    for key in [
        "policy",
        "machine",
        "adapter",
        "checkpoint",
        "cli",
        "version",
        "backend",
        "graph",
        "sample",
        "depth_bits",
        "background_bits",
        "tie_bits",
        "max_feature_terms",
        "cache_terms",
        "count_policy",
        "drop-task",
        "native-bits",
        "modes",
    ] {
        let mut bad = original.clone();
        match key {
            "policy" | "machine" | "adapter" | "checkpoint" | "cli" | "version" => {
                bad["bindings"]["policy"][key] = "altered".into()
            }
            "depth_bits" | "background_bits" | "tie_bits" | "max_feature_terms" | "cache_terms" => {
                bad["bindings"][key] = 1.into()
            }
            "drop-task" => {
                let id = bad["state"]["tasks"]
                    .as_object()
                    .unwrap()
                    .keys()
                    .next()
                    .unwrap()
                    .clone();
                bad["state"]["tasks"].as_object_mut().unwrap().remove(&id);
            }
            "native-bits" => bad["state"]["native"][0]["objective_bits"] = 0.into(),
            "modes" => bad["state"]["families"][0]["modes"][2] = 999999.into(),
            _ => bad["bindings"][key] = "altered".into(),
        }
        write(&path, &bad);
        let mut seal: Value = serde_json::from_slice(&original_seal).unwrap();
        seal["checkpoint"] = hash(&path).into();
        write(&sealpath, &seal);
        guided(
            &f,
            binary,
            0,
            &format!("corruption-{key}"),
            &[
                "--resume-from",
                p(old_root),
                "--extend-budgets",
                "--convert-exact-v1-to-v2",
                "--v1-ancestry-manifest",
                p(&declaration),
            ],
            false,
        );
        fs::write(&path, &original_bytes).unwrap();
        fs::write(&sealpath, &original_seal).unwrap();
    }
    guided(
        &f,
        binary,
        0,
        "missing-conversion-extension",
        &[
            "--resume-from",
            p(old_root),
            "--convert-exact-v1-to-v2",
            "--v1-ancestry-manifest",
            p(&declaration),
        ],
        false,
    );
    // Independent pins, interrupted writes and corrupted recognized boundaries fail closed.
    let missing = root.join("incomplete-pins.json");
    write(&missing, &json!({"version":1,"checkpoints":[]}));
    guided(
        &f,
        binary,
        0,
        "reject-incomplete-pins",
        &[
            "--resume-from",
            p(old_root),
            "--extend-budgets",
            "--convert-exact-v1-to-v2",
            "--v1-ancestry-manifest",
            p(&missing),
        ],
        false,
    );
    fs::write(&path, b"{\"version\":1,").unwrap();
    guided(
        &f,
        binary,
        0,
        "reject-interrupted-old",
        &[
            "--resume-from",
            p(old_root),
            "--extend-budgets",
            "--convert-exact-v1-to-v2",
            "--v1-ancestry-manifest",
            p(&declaration),
        ],
        false,
    );
    fs::write(&path, &original_bytes).unwrap();
    // A cumulative history must retain the exact earlier bytes, even when the
    // local file and its checksum agree. Check both conversion and later resume.
    for name in ["converted-2048", "resumed-a-2048"] {
        let directory = root.join(name);
        let ledger = directory.join("evaluations.jsonl");
        let sealpath = directory.join("checkpoint-seal.json");
        let original = fs::read_to_string(&ledger).unwrap();
        let original_seal = fs::read(&sealpath).unwrap();
        for change in ["objective", "flag", "spacing"] {
            let mut lines: Vec<String> = original.lines().map(str::to_owned).collect();
            if change == "spacing" {
                lines[2].insert(0, ' ');
            } else {
                let mut row: Value = serde_json::from_str(&lines[2]).unwrap();
                assert_eq!(row["kind"], "native-score-reuse");
                if change == "objective" {
                    row["scored"]["relative_objective"] = json!(0.0);
                } else {
                    row["cross_source_route"] = json!(true);
                }
                lines[2] = serde_json::to_string(&row).unwrap();
            }
            fs::write(&ledger, lines.join("\n") + "\n").unwrap();
            let mut seal: Value = serde_json::from_slice(&original_seal).unwrap();
            seal["ledger"] = hash(&ledger).into();
            write(&sealpath, &seal);
            guided(
                &f,
                binary,
                0,
                &format!("inconsistent-history-{name}-{change}"),
                &["--resume-from", p(&directory)],
                false,
            );
            fs::write(&ledger, &original).unwrap();
            fs::write(&sealpath, &original_seal).unwrap();
        }
    }
    let boundary = root.join("converted-2048");
    let boundarypath = boundary.join("checkpoint.json");
    let boundaryseal = boundary.join("checkpoint-seal.json");
    let bb = fs::read(&boundarypath).unwrap();
    let bs = fs::read(&boundaryseal).unwrap();
    for kind in ["receipt", "link"] {
        let mut cp: Value = serde_json::from_slice(&bb).unwrap();
        if kind == "receipt" {
            cp["transition"]["conversion_only"] = false.into();
        } else {
            cp["state"]["focus"][0] = 999999.into();
        }
        write(&boundarypath, &cp);
        let mut seal: Value = serde_json::from_slice(&bs).unwrap();
        seal["checkpoint"] = hash(&boundarypath).into();
        write(&boundaryseal, &seal);
        guided(
            &f,
            binary,
            0,
            &format!("reject-v2-{kind}"),
            &["--resume-from", p(&boundary)],
            false,
        );
        fs::write(&boundarypath, &bb).unwrap();
        fs::write(&boundaryseal, &bs).unwrap();
    }
    let resumed = root.join("resumed-a-2048");
    let resumed_again = guided(
        &f,
        binary,
        0,
        "v2-post-transition-second-epoch",
        &["--resume-from", p(&resumed)],
        true,
    );
    assert_eq!(
        read(&resumed.join("checkpoint.json"))["state"],
        read(&resumed_again.join("checkpoint.json"))["state"]
    );
    let near = guided(
        &f,
        &old,
        0,
        "v1-near-cap",
        &["--max-state-bytes", "32768"],
        true,
    );
    assert_eq!(
        read(&near.join("result.json"))["result"]["status"],
        "state-budget-exhausted"
    );
    let nearpins = root.join("near-cap-pins.json");
    write(
        &nearpins,
        &json!({"version":1,"checkpoints":[{"directory":fs::canonicalize(&near).unwrap(),"checkpoint":hash(&near.join("checkpoint.json")),"ledger":hash(&near.join("evaluations.jsonl"))}]}),
    );
    guided(
        &f,
        binary,
        0,
        "near-cap-no-growth",
        &[
            "--max-state-bytes",
            "32768",
            "--resume-from",
            p(&near),
            "--convert-exact-v1-to-v2",
            "--v1-ancestry-manifest",
            p(&nearpins),
        ],
        false,
    );
    let converted = guided(
        &f,
        binary,
        0,
        "near-cap-explicit-growth",
        &[
            "--resume-from",
            p(&near),
            "--extend-budgets",
            "--convert-exact-v1-to-v2",
            "--v1-ancestry-manifest",
            p(&nearpins),
        ],
        true,
    );
    assert_eq!(
        fs::read(near.join("evaluations.jsonl")).unwrap(),
        fs::read(converted.join("evaluations.jsonl")).unwrap()
    );
    let continued = guided(
        &f,
        binary,
        0,
        "near-cap-continued",
        &["--resume-from", p(&converted)],
        true,
    );
    assert!(
        read(&continued.join("result.json"))["result"]["work"]
            .as_u64()
            .unwrap()
            > read(&near.join("result.json"))["result"]["work"]
                .as_u64()
                .unwrap()
    );
}

#[test]
fn guided_v2_every_primitive_resume_state_and_trace() {
    let t = tempfile::tempdir().unwrap();
    let preserved = std::env::var_os("IMPG_TEST_V2_BOUNDARY_OUTPUT").map(PathBuf::from);
    let root = preserved.unwrap_or_else(|| t.path().join("boundaries"));
    let f = fixture_with_slots(&root, 1, 2);
    let binary = env!("CARGO_BIN_EXE_impg");
    let mut prior: Option<PathBuf> = None;
    let mut phases = std::collections::BTreeSet::new();
    let mut modes = std::collections::BTreeSet::new();
    for work in 1..=384 {
        let cap = work.to_string();
        let fresh = guided(
            &f,
            binary,
            0,
            &format!("fresh-{work}"),
            &["--max-work", &cap],
            true,
        );
        let cp = read(&fresh.join("checkpoint.json"));
        if let Some(ref previous) = prior {
            let resumed = guided(
                &f,
                binary,
                0,
                &format!("resume-{work}"),
                &[
                    "--max-work",
                    &cap,
                    "--resume-from",
                    p(previous),
                    "--extend-budgets",
                ],
                true,
            );
            assert_eq!(
                fs::read(fresh.join("evaluations.jsonl")).unwrap(),
                fs::read(resumed.join("evaluations.jsonl")).unwrap(),
                "ledger at work {work}"
            );
            assert_eq!(
                cp["state"],
                read(&resumed.join("checkpoint.json"))["state"],
                "state at work {work}"
            );
        }
        for task in cp["state"]["tasks"].as_object().unwrap().values() {
            let op = &task["op"];
            phases.insert(if op.is_string() {
                op.as_str().unwrap().to_string()
            } else {
                op.as_object().unwrap().keys().next().unwrap().clone()
            });
        }
        for family in cp["state"]["families"].as_array().unwrap() {
            modes.insert(family["work"].as_u64().unwrap() % 3);
        }
        prior = Some(fresh);
    }
    for phase in [
        "Start",
        "SourceBound",
        "SourceScan",
        "Check",
        "HubBound",
        "HubScan",
        "Child",
        "Close",
        "Probe",
        "ProbeCheck",
        "Evaluate",
    ] {
        assert!(
            phases.contains(phase),
            "missing primitive boundary {phase}: {phases:?}"
        );
    }
    assert_eq!(modes.len(), 3);
    write(
        &root.join("boundary-coverage.json"),
        &json!({"every_work_boundary":384,"modes":modes,"operations":phases,"comparison":"fresh versus resumed exact state and cumulative semantic ledger at every primitive"}),
    );
}
