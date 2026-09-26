use super::*;

fn conserved(old: &Value, new: &Value) {
    assert_eq!(old["bindings"]["backend"], new["bindings"]["backend"]);
    let mut expected = old["state"].clone();
    for task in expected["tasks"].as_object_mut().unwrap().values_mut() {
        if let Some(scan) = task["op"].get_mut("SourceScan") {
            let permutation = &mut scan["permutation"];
            permutation["old_prefix"] = permutation["ordinal"].clone();
            permutation["ordinal"] = 0.into();
        }
    }
    // These are the only accounting fields changed by additional cursor storage.
    expected["task_bytes"] = new["state"]["task_bytes"].clone();
    expected["peak_state_bytes"] = new["state"]["peak_state_bytes"].clone();
    assert_eq!(expected, new["state"]);
    assert_eq!(new["order_update"]["old_bindings"], old["bindings"]);
    assert_eq!(new["order_update"]["new_bindings"], new["bindings"]);
    assert_eq!(new["version"], 3);
}
fn update(f: &Fixture, parent: &Path, label: &str, cap: u64) -> PathBuf {
    let old = read(&parent.join("checkpoint.json"));
    let updated = guided(
        f,
        env!("CARGO_BIN_EXE_impg"),
        0,
        label,
        &[
            "--resume-from",
            p(parent),
            "--update-source-pair-order",
            "--max-work",
            &cap.to_string(),
        ],
        true,
    );
    conserved(&old, &read(&updated.join("checkpoint.json")));
    assert_eq!(
        fs::read(parent.join("evaluations.jsonl")).unwrap(),
        fs::read(updated.join("evaluations.jsonl")).unwrap()
    );
    assert_eq!(
        read(&updated.join("result.json"))["result"]["update_only"],
        true
    );
    updated
}
fn reverse_scan(cp: &Value) -> bool {
    cp["state"]["tasks"].as_object().unwrap().values().any(|t| {
        t["context"]["reverse"] == true
            && t["op"]["SourceScan"]["permutation"]["ordinal"]
                .as_u64()
                .is_some_and(|n| n > 0)
    })
}
fn reverse_mixed(ledger: &[u8]) -> Vec<Value> {
    ledger
        .split(|&b| b == b'\n')
        .filter(|s| !s.is_empty())
        .map(|s| serde_json::from_slice::<Value>(s).unwrap())
        .filter(|r| {
            r["cross_source_route"] == true
                && r["scored"]["assignment"]["routes"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|route| {
                        let segments = route["segments"].as_array().unwrap();
                        segments.iter().any(|s| {
                            s["reverse"] == true
                                && s["source"] != segments[0]["source"]
                                && s["end"].as_u64().unwrap() > s["start"].as_u64().unwrap()
                        })
                    })
        })
        .collect()
}

#[test]
#[ignore = "synthetic long reverse gate; needs fresh IMPG_TEST_PAIR_REVERSE_OUTPUT"]
fn source_pair_order_long_reverse_mixed_gate() {
    // Recipe and caps are fixed before graph construction/search outcomes. No topology edits.
    let root = PathBuf::from(
        std::env::var("IMPG_TEST_PAIR_REVERSE_OUTPUT").expect("fresh reverse output required"),
    );
    let f = fixture_with_orientation(&root, 1000, 1, true);
    let identity = genome::PanelIdentity::read(p(&f.panel)).unwrap();
    let g = routes::Graph::load(&f.graph, &identity).unwrap();
    write(
        &root.join("prerequisites-before-search.json"),
        &json!({
        "recipe":"seed113 400bp; B mutations150,260 then reverse complement; both repeated1000",
        "work_cap":200000,"evaluation_cap":10000,"state_cap":268435456,"optima_cap":20000,
        "port_counts":g.lanes.iter().map(|l| l.port_count).collect::<Vec<_>>() }),
    );
    assert_eq!(g.families.len(), 2);
    assert!(g.lanes.iter().all(|l| l.port_count >= 32768));
    for lane in &g.lanes {
        let raw = fs::read(f.graph.join(&lane.ports.path)).unwrap();
        let width = g.k as usize + 17;
        for pair in raw.chunks_exact(2 * width) {
            assert_eq!(
                &pair[g.k as usize..width - 1],
                &pair[width + g.k as usize..2 * width - 1]
            );
            assert_eq!(pair[width - 1], 0);
            assert_eq!(pair[2 * width - 1], 1);
        }
    }
    let binary = env!("CARGO_BIN_EXE_impg");
    let mut previous = None;
    for repeat in 0..2 {
        let out = guided(
            &f,
            binary,
            0,
            &format!("reverse-{repeat}"),
            &["--max-work", "200000"],
            true,
        );
        let ledger = fs::read(out.join("evaluations.jsonl")).unwrap();
        let mixed = reverse_mixed(&ledger);
        assert!(
            !mixed.is_empty(),
            "fixed reverse gate did not complete positive reverse donor traversal"
        );
        write(
            &out.join("reverse-mixed-evidence.json"),
            &json!({"count":mixed.len(),"first":mixed[0]}),
        );
        if let Some(prior) = previous {
            assert_eq!(ledger, prior);
        }
        previous = Some(ledger);
        assert!(
            reverse_scan(&read(&out.join("checkpoint.json"))),
            "long reverse work must remain live"
        );
    }
}

#[test]
#[ignore = "independent v2 update gate; needs IMPG_TEST_GUIDED_V2, IMPG_TEST_GUIDED_V1 and fresh IMPG_TEST_PAIR_UPDATE_OUTPUT"]
fn source_pair_order_v2_update_and_v1_ancestry() {
    let root = PathBuf::from(
        std::env::var("IMPG_TEST_PAIR_UPDATE_OUTPUT").expect("fresh update output required"),
    );
    let v2 = std::env::var("IMPG_TEST_GUIDED_V2").expect("preserved v2 binary required");
    let v1 = std::env::var("IMPG_TEST_GUIDED_V1").expect("preserved v1 binary required");
    let f = fixture_with_orientation(&root, 100, 1, true);
    let binary = env!("CARGO_BIN_EXE_impg");
    write(
        &root.join("update-predeclaration.json"),
        &json!({"recipe":"same seed113 reverse recipe; repeats100",
        "old_work_caps":[1,2,5,32,128,512,2048,8192],"resume_cap":20000,"evaluation_cap":10000,"state_cap":268435456,
        "requirements":["SourceBound","SourceScan","nonzero focused reverse cursor","v1-v2-v3 ancestry","update-only conservation","exact repeated/resumed continuation"]}),
    );
    let mut saw_reverse = false;
    let mut saw_focused_reverse = false;
    let mut ops = std::collections::BTreeSet::new();
    for cap in [1u64, 2, 5, 32, 128, 512, 2048, 8192] {
        let old = guided(
            &f,
            &v2,
            0,
            &format!("v2-{cap}"),
            &["--max-work", &cap.to_string()],
            true,
        );
        let cp = read(&old.join("checkpoint.json"));
        saw_reverse |= reverse_scan(&cp);
        for t in cp["state"]["tasks"].as_object().unwrap().values() {
            if let Some(op) = t["op"].as_object().and_then(|o| o.keys().next()) {
                ops.insert(op.clone());
            }
        }
        for id in cp["state"]["focus"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(Value::as_u64)
        {
            let t = &cp["state"]["tasks"][id.to_string()];
            saw_focused_reverse |= t["context"]["reverse"] == true
                && t["op"]["SourceScan"]["permutation"]["ordinal"]
                    .as_u64()
                    .is_some_and(|n| n > 0);
        }
        guided(
            &f,
            binary,
            0,
            &format!("reject-ordinary-v2-{cap}"),
            &["--resume-from", p(&old), "--max-work", &cap.to_string()],
            false,
        );
        let a = update(&f, &old, &format!("update-a-{cap}"), cap);
        let b = update(&f, &old, &format!("update-b-{cap}"), cap);
        assert_eq!(
            read(&a.join("checkpoint.json"))["state"],
            read(&b.join("checkpoint.json"))["state"]
        );
        guided(
            &f,
            binary,
            0,
            &format!("reject-repeat-update-{cap}"),
            &[
                "--resume-from",
                p(&a),
                "--update-source-pair-order",
                "--max-work",
                &cap.to_string(),
            ],
            false,
        );
        let direct = guided(
            &f,
            binary,
            0,
            &format!("direct-{cap}"),
            &["--resume-from", p(&a), "--extend-budgets"],
            true,
        );
        let halfway = guided(
            &f,
            binary,
            0,
            &format!("halfway-{cap}"),
            &[
                "--resume-from",
                p(&b),
                "--extend-budgets",
                "--max-work",
                "10000",
            ],
            true,
        );
        let resumed = guided(
            &f,
            binary,
            0,
            &format!("resumed-{cap}"),
            &["--resume-from", p(&halfway), "--extend-budgets"],
            true,
        );
        assert_eq!(
            read(&direct.join("checkpoint.json"))["state"],
            read(&resumed.join("checkpoint.json"))["state"]
        );
        assert_eq!(
            fs::read(direct.join("evaluations.jsonl")).unwrap(),
            fs::read(resumed.join("evaluations.jsonl")).unwrap()
        );
    }
    assert!(
        saw_reverse && saw_focused_reverse,
        "predeclared caps must include nonzero focused reverse SourceScan"
    );
    assert!(ops.contains("SourceBound") && ops.contains("SourceScan"));
    let first = guided(&f, &v1, 0, "v1-first", &["--max-work", "128"], true);
    let second = guided(
        &f,
        &v1,
        0,
        "v1-second",
        &[
            "--resume-from",
            p(&first),
            "--extend-budgets",
            "--max-work",
            "512",
        ],
        true,
    );
    let pins = pin(&root, &[first, second.clone()]);
    let converted = guided(
        &f,
        &v2,
        0,
        "v1-converted-v2",
        &[
            "--resume-from",
            p(&second),
            "--convert-exact-v1-to-v2",
            "--v1-ancestry-manifest",
            p(&pins),
            "--max-work",
            "512",
        ],
        true,
    );
    let ongoing = guided(
        &f,
        &v2,
        0,
        "v2-after-conversion",
        &[
            "--resume-from",
            p(&converted),
            "--extend-budgets",
            "--max-work",
            "8192",
        ],
        true,
    );
    let updated = update(&f, &ongoing, "v1-v2-v3-update", 8192);
    let a = guided(
        &f,
        binary,
        0,
        "ancestry-resume-a",
        &["--resume-from", p(&updated), "--extend-budgets"],
        true,
    );
    let b = guided(
        &f,
        binary,
        0,
        "ancestry-resume-b",
        &["--resume-from", p(&updated), "--extend-budgets"],
        true,
    );
    assert_eq!(
        read(&a.join("checkpoint.json"))["state"],
        read(&b.join("checkpoint.json"))["state"]
    );
    assert_eq!(
        fs::read(a.join("evaluations.jsonl")).unwrap(),
        fs::read(b.join("evaluations.jsonl")).unwrap()
    );
    write(
        &root.join("update-evidence.json"),
        &json!({"nonzero_reverse":saw_reverse,"nonzero_focused_reverse":saw_focused_reverse,"operations":ops,"v1_v2_v3_resume":true}),
    );
}
