//! Three successive donors sharing one distinct port per pair, no two-switch bridge.
use super::*;
pub(super) fn no_simple_return_gate(root: &Path, words: Vec<Vec<u8>>) {
    let root = root.join("no-simple-return");
    fs::create_dir(&root).unwrap();
    let mut sequences = Vec::new();
    for (i, (left, right)) in [(0, 2), (0, 1), (1, 2)].into_iter().enumerate() {
        let mut seq = dna(400, 1009 + i as u64 * 77);
        seq[90..90 + words[left].len()].copy_from_slice(&words[left]);
        seq[250..250 + words[right].len()].copy_from_slice(&words[right]);
        sequences.push((format!("{}#0#chain", ["A", "B", "C"][i]), seq));
    }
    let mut panel = SyngIndex::build(SyncmerParams::default(), sequences.clone().into_iter());
    let prefix = root.join("panel.syng");
    panel.save(p(&prefix)).unwrap();
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    let fasta = root.join("sources.fa");
    fs::write(
        &fasta,
        sequences
            .iter()
            .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    let groups = root.join("groups.json");
    write(
        &groups,
        &json!({"version":1,"groups":sequences.iter().enumerate().map(|(i,(n,s))|json!({"id":format!("g{i}"),"scaffold":null,"occurrences":[{"path":n,"start":0,"end":s.len(),"strand":null}]})).collect::<Vec<_>>()}),
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
        true,
    );
    let graph_dir = root.join("routes");
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
            p(&fasta),
            "--read-lengths",
            "150",
            "--out-dir",
            p(&graph_dir),
        ],
        true,
    );
    let graph = routes::Graph::load(&graph_dir, &identity).unwrap();
    let ports = ports(&graph_dir, &graph);
    let mut expected = BTreeSet::new();
    for family in 0..graph.families.len() {
        let target = graph.families[family].paths[0];
        let rows = independent_routes(&graph, &ports, target);
        if family == 0 {
            assert!(
                rows.iter().any(|r| r.segments.len() >= 4),
                "three-donor route missing"
            );
        }
        assert!(
            rows.iter().all(|r| r.segments.len() != 3),
            "fixture unexpectedly has a simple donor-return bridge"
        );
        for route in rows {
            expected.insert(routes::Assignment {
                version: 1,
                model: routes::MODEL.into(),
                graph_checksum: graph.digest().unwrap(),
                family,
                routes: vec![route],
            });
        }
    }
    let reads = root.join("reads.fa");
    fs::write(
        &reads,
        sequences
            .iter()
            .flat_map(|(_, s)| s.windows(150))
            .map(|s| format!(">read\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    let sample = sample::build(&panel, identity, &[reads]).unwrap();
    let sample_path = root.join("sample.membwt");
    sample.save(&sample_path).unwrap();
    let out = root.join("searched");
    guided(
        &prefix,
        &graph_dir,
        &sample_path,
        &out,
        &["--max-work", "100000"],
        true,
    );
    let result = read(&out.join("result.json"));
    assert_eq!(result["result"]["search_exhausted"], true);
    let actual: BTreeSet<routes::Assignment> = fs::read_to_string(out.join("evaluations.jsonl"))
        .unwrap()
        .lines()
        .map(|line| {
            let row: Value = serde_json::from_str(line).unwrap();
            serde_json::from_value(row["scored"]["assignment"].clone()).unwrap()
        })
        .collect();
    assert_eq!(actual, expected);
    assert!(
        result["result"]["maximum_switches_evaluated"]
            .as_u64()
            .unwrap()
            >= 3
    );
    eprintln!("No simple return bridge: {} exhaustive assignments, {} mixed complete assignments; three-donor construction independently matched",actual.len(),result["result"]["mixed_assignments_evaluated"]);
}
