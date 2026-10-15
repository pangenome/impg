//! New guided-policy CLI gates. Existing lexical fixtures/helpers remain untouched.
use super::*;
mod no_return;

pub(super) fn configured_748() -> String {
    let binary = std::env::var("IMPG_TEST_748_ROUTES")
        .expect("requested frozen-748 gate requires IMPG_TEST_748_ROUTES");
    assert!(Path::new(&binary).is_absolute() && Path::new(&binary).is_file());
    let output = Command::new("sha256sum")
        .arg(&binary)
        .output()
        .expect("explicit Linux gate needs sha256sum");
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout)
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap(),
        "ea22ec02df59715dd852d77056e0a757955d85db811bed24fa4b05eda55a0eb7"
    );
    binary
}
fn guided(prefix: &Path, graph: &Path, sample: &Path, out: &Path, extra: &[&str], ok: bool) {
    let mut args = vec![
        "genome-infer",
        "search-panel-routes-guided",
        "--panel",
        p(prefix),
        "--routes",
        p(graph),
        "--sample",
        p(sample),
        "--haploid-depth",
        "150",
        "--max-optima",
        "20000",
        "--out-dir",
        p(out),
    ];
    args.extend_from_slice(extra);
    run(env!("CARGO_BIN_EXE_impg"), &args, ok);
}
pub(super) fn resume_gates(root: &Path, prefix: &Path, graph: &Path, sample: &Path, whole: &Path) {
    let whole_ledger = fs::read(whole.join("evaluations.jsonl")).unwrap();
    // Base state fits, but the conservative first-operation reservation does not.
    // Even a zero-operation checkpoint must include retained state in its peak.
    let zero = root.join("zero-work-resource-stop");
    guided(
        prefix,
        graph,
        sample,
        &zero,
        &["--max-state-bytes", "10000"],
        true,
    );
    let zero_result = read(&zero.join("result.json"))["result"].clone();
    assert_eq!(zero_result["status"], "state-budget-exhausted");
    assert_eq!(zero_result["work"], 0);
    assert_eq!(zero_result["evaluations"], 0);
    assert!(
        zero_result["peak_state_bytes"].as_u64().unwrap()
            >= zero_result["retained_state_bytes"].as_u64().unwrap(),
        "zero-operation peak underreports retained occupancy: {zero_result}"
    );
    let zero_resumed = root.join("zero-work-resource-resumed");
    guided(
        prefix,
        graph,
        sample,
        &zero_resumed,
        &[
            "--max-state-bytes",
            "268435456",
            "--max-work",
            "5000000",
            "--max-evaluations",
            "100000",
            "--resume-from",
            p(&zero),
            "--extend-budgets",
        ],
        true,
    );
    assert_eq!(
        fs::read(zero_resumed.join("evaluations.jsonl")).unwrap(),
        whole_ledger
    );
    assert_eq!(
        read(&zero_resumed.join("checkpoint.json"))["state"],
        read(&whole.join("checkpoint.json"))["state"]
    );
    let mut parent: Option<std::path::PathBuf> = None;
    let mut phases = BTreeSet::new();
    for work in [1, 3, 4, 5, 6, 7, 8, 16, 32, 64, 128, 256, 512, 1024, 2048] {
        let out = root.join(format!("resume-{work}"));
        let work = work.to_string();
        let mut extra = vec!["--max-work", &work, "--max-evaluations", "100000"];
        if let Some(prior) = &parent {
            extra.extend(["--resume-from", p(prior), "--extend-budgets"]);
        }
        guided(prefix, graph, sample, &out, &extra, true);
        let ledger = fs::read(out.join("evaluations.jsonl")).unwrap();
        assert!(
            whole_ledger.starts_with(&ledger),
            "semantic prefix changed at {work}"
        );
        let cp = read(&out.join("checkpoint.json"));
        for task in cp["state"]["tasks"].as_object().unwrap().values() {
            let op = &task["op"];
            phases.insert(if op.is_string() {
                op.as_str().unwrap().to_string()
            } else {
                op.as_object().unwrap().keys().next().unwrap().clone()
            });
        }
        parent = Some(out);
    }
    for phase in [
        "SourceBound",
        "SourceScan",
        "HubBound",
        "HubScan",
        "Child",
        "Check",
        "Close",
        "Evaluate",
    ] {
        assert!(
            phases.contains(phase),
            "missing checkpoint boundary {phase}: {phases:?}"
        );
    }
    let parent = parent.unwrap();
    let final_out = root.join("resumed-complete");
    guided(
        prefix,
        graph,
        sample,
        &final_out,
        &[
            "--max-work",
            "5000000",
            "--max-evaluations",
            "100000",
            "--resume-from",
            p(&parent),
            "--extend-budgets",
        ],
        true,
    );
    assert_eq!(
        fs::read(final_out.join("evaluations.jsonl")).unwrap(),
        whole_ledger
    );
    assert_eq!(
        read(&final_out.join("checkpoint.json"))["state"],
        read(&whole.join("checkpoint.json"))["state"]
    );
    // Fresh output collision, missing explicit extension, and scientific-parameter changes fail.
    guided(prefix, graph, sample, &final_out, &[], false);
    guided(
        prefix,
        graph,
        sample,
        &root.join("no-extension"),
        &[
            "--max-work",
            "5000000",
            "--max-evaluations",
            "100000",
            "--resume-from",
            p(&parent),
        ],
        false,
    );
    let parameter_out = root.join("changed-parameter");
    let args = [
        "genome-infer",
        "search-panel-routes-guided",
        "--panel",
        p(prefix),
        "--routes",
        p(graph),
        "--sample",
        p(sample),
        "--haploid-depth",
        "151",
        "--out-dir",
        p(&parameter_out),
        "--resume-from",
        p(&parent),
        "--extend-budgets",
    ];
    run(env!("CARGO_BIN_EXE_impg"), &args, false);
    for name in ["checkpoint.json", "evaluations.jsonl"] {
        let path = parent.join(name);
        let original = fs::read(&path).unwrap();
        fs::write(&path, b"corrupt\n").unwrap();
        guided(
            prefix,
            graph,
            sample,
            &root.join(format!("corrupt-{name}")),
            &[
                "--max-work",
                "5000000",
                "--max-evaluations",
                "100000",
                "--resume-from",
                p(&parent),
                "--extend-budgets",
            ],
            false,
        );
        // An altered ancestor must also reject a sealed child.
        guided(
            prefix,
            graph,
            sample,
            &root.join(format!("ancestor-{name}")),
            &[
                "--max-work",
                "5000000",
                "--max-evaluations",
                "100000",
                "--resume-from",
                p(&final_out),
            ],
            false,
        );
        fs::write(&path, original).unwrap();
    }
    // Modify a binding while recomputing its corruption seal: binding checks are separate.
    for binding in [
        "policy",
        "backend",
        "graph",
        "sample",
        "depth_bits",
        "tie_bits",
        "max_feature_terms",
    ] {
        let path = parent.join("checkpoint.json");
        let original = fs::read(&path).unwrap();
        let sealpath = parent.join("checkpoint-seal.json");
        let seal = fs::read(&sealpath).unwrap();
        let mut cp: Value = serde_json::from_slice(&original).unwrap();
        if binding == "policy" {
            cp["bindings"][binding]["version"] = "changed".into();
        } else if binding.ends_with("bits") || binding == "max_feature_terms" {
            cp["bindings"][binding] = 1.into();
        } else {
            cp["bindings"][binding] = "changed".into();
        }
        write(&path, &cp);
        let changed = fs::read(&path).unwrap();
        let mut se: Value = serde_json::from_slice(&seal).unwrap();
        se["checkpoint"] = format!(
            "fnv1a64:{}:{:016x}",
            changed.len(),
            genome::checksum(&changed)
        )
        .into();
        write(&sealpath, &se);
        guided(
            prefix,
            graph,
            sample,
            &root.join(format!("binding-{binding}")),
            &[
                "--max-work",
                "5000000",
                "--max-evaluations",
                "100000",
                "--resume-from",
                p(&parent),
                "--extend-budgets",
            ],
            false,
        );
        fs::write(path, original).unwrap();
        fs::write(sealpath, seal).unwrap();
    }
    // Resource/evaluation stops retain initialization and all pending work.
    for (label, flag, cap) in [
        ("resource", "--max-state-bytes", "32768"),
        ("evaluation", "--max-evaluations", "1"),
    ] {
        let stopped = root.join(format!("stop-{label}"));
        guided(
            prefix,
            graph,
            sample,
            &stopped,
            &["--max-work", "5000000", flag, cap],
            true,
        );
        assert_eq!(
            read(&stopped.join("result.json"))["result"]["search_exhausted"],
            false
        );
        let resumed = root.join(format!("extend-{label}"));
        guided(
            prefix,
            graph,
            sample,
            &resumed,
            &[
                "--max-work",
                "5000000",
                "--max-evaluations",
                "100000",
                "--resume-from",
                p(&stopped),
                "--extend-budgets",
            ],
            true,
        );
        assert_eq!(
            fs::read(resumed.join("evaluations.jsonl")).unwrap(),
            whole_ledger
        );
        assert_eq!(
            read(&resumed.join("checkpoint.json"))["state"],
            read(&whole.join("checkpoint.json"))["state"]
        );
    }
    eprintln!("Resume semantic trace exact through phases {phases:?}; resource/evaluation/initialization and ancestry checks passed");
}
pub(super) fn paired_rank_and_bounded_gates(
    root: &Path,
    prefix: &Path,
    graph_dir: &Path,
    panel: &SyngIndex,
    graph: &routes::Graph,
    a: &[u8],
    b: &[u8],
) {
    // Predeclared engineering gate: 20,000 primitive work; 10,000 fresh calls;
    // 256 MiB logical retained state. Require real mixed completions in A AND B.
    let mut orders = Vec::new();
    for (label, sequence) in [("a", a), ("b", b)] {
        let reads = root.join(format!("rank-{label}.fa"));
        fs::write(
            &reads,
            sequence
                .windows(150)
                .map(|s| format!(">read\n{}\n", String::from_utf8_lossy(s)))
                .collect::<String>(),
        )
        .unwrap();
        let sample = sample::build(panel, graph.panel.clone(), &[reads]).unwrap();
        let sample_path = root.join(format!("rank-{label}.membwt"));
        sample.save(&sample_path).unwrap();
        let mut first = None;
        for repeat in 0..2 {
            let out = root.join(format!("rank-{label}-{repeat}"));
            guided(
                prefix,
                graph_dir,
                &sample_path,
                &out,
                &["--max-work", "20000"],
                true,
            );
            let r = read(&out.join("result.json"))["result"].clone();
            let order = r["native_ranked_identities"].clone();
            let ledger = fs::read(out.join("evaluations.jsonl")).unwrap();
            if repeat == 1 {
                assert_eq!(
                    ledger,
                    fs::read(root.join(format!("rank-{label}-0/evaluations.jsonl"))).unwrap()
                );
            }
            let first_serviced: Value = String::from_utf8(ledger)
                .unwrap()
                .lines()
                .map(|l| serde_json::from_str::<Value>(l).unwrap())
                .find(|row| row["kind"] != "native-initialization")
                .unwrap();
            assert_eq!(
                first_serviced["scored"]["assignment"]["family"],
                r["native_ranked_family_order"][0]
            );
            if let Some(previous) = &first {
                assert_eq!(previous, &order);
            } else {
                first = Some(order.clone());
            }
            assert!(
                r["family_statistics"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .filter(|f| f["mixed"].as_u64().unwrap() > 0)
                    .count()
                    >= 2,
                "no useful multi-family mixed exploration: {r}"
            );
            eprintln!("Bounded sample {label} repeat {repeat}: order={order}; mixed={}; switches={}; donors={}; peak_tasks={}; peak_bytes={}",r["mixed_assignments_evaluated"],r["switch_distribution"],r["mixed_donor_source_ids"],r["peak_frontier_tasks"],r["peak_state_bytes"]);
        }
        // A genuinely different sample payload (same graph/parameters) cannot resume.
        guided(
            prefix,
            graph_dir,
            &sample_path,
            &root.join(format!("reject-actual-sample-{label}")),
            &[
                "--max-work",
                "5000000",
                "--max-evaluations",
                "100000",
                "--resume-from",
                p(&root.join("searched")),
            ],
            false,
        );
        orders.push(first.unwrap());
    }
    assert_ne!(
        orders[0], orders[1],
        "samples must change native-ranked service order"
    );
    assert_eq!(orders[0][0], "A#0");
    assert_eq!(orders[1][0], "B#0");
    let words: BTreeSet<_> = ports(graph_dir, graph)
        .into_iter()
        .filter(|p| !p.reverse)
        .map(|p| p.word)
        .collect();
    no_return::no_simple_return_gate(root, words.into_iter().take(3).collect());
}
