//! Native public-sample oracle and real CLI contract. Run serialized.
use impg::{
    genome_inference::{self as genome, joint, sample},
    syng::{SyncmerParams, SyngIndex},
};
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, path::Path, process::Command};
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
fn p(p: &Path) -> &str {
    p.to_str().unwrap()
}
fn write(path: &Path, v: &Value) {
    fs::write(path, serde_json::to_vec_pretty(v).unwrap()).unwrap();
}
fn read(path: &Path) -> Value {
    genome::read_json(path).unwrap()
}
fn run(args: &[&str], ok: bool) {
    let o = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(args)
        .output()
        .unwrap();
    assert_eq!(
        o.status.success(),
        ok,
        "{args:?}\n{}",
        String::from_utf8_lossy(&o.stderr)
    );
}
fn piece(instance: &str, start: u64, end: u64) -> Value {
    json!({"instance":instance,"start":start,"end":end,"strand":"+"})
}
fn walk(id: &str, pieces: Vec<Value>) -> Value {
    let adj: Vec<_> = (0..pieces.len() - 1)
        .map(|n| json!({"from":n,"to":n+1,"kind":"abut"}))
        .collect();
    json!({"id":id,"topology":"linear","left_endpoint":"asserted-molecule-terminus","right_endpoint":"asserted-molecule-terminus","pieces":pieces,"adjacencies":adj})
}
fn oracle(
    panel: &SyngIndex,
    identity: &genome::PanelIdentity,
    path: &Path,
    sequence: &[u8],
) -> sample::SampleIndex {
    fs::write(
        path,
        sequence
            .windows(250)
            .map(|s| format!(">transient\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    sample::build(panel, identity.clone(), &[path.to_path_buf()]).unwrap()
}
#[test]
fn joint_cli_mixed_source_feature_closure_native_exhaustive_objective_and_legacy_rejection() {
    let temp = tempfile::tempdir().unwrap();
    let root = temp.path();
    let a = dna(900, 31);
    let b = dna(900, 39);
    let mut bridge = a[300..450].to_vec();
    bridge.extend_from_slice(&b[450..600]);
    let sequences = vec![
        ("A#0#chr".to_string(), a.clone()),
        ("B#0#chr".to_string(), b.clone()),
        ("C#0#bridge".to_string(), bridge),
    ];
    let prefix = root.join("panel.syng");
    let mut panel = SyngIndex::build(SyncmerParams::default(), sequences.clone().into_iter());
    panel.save(p(&prefix)).unwrap();
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    let source = root.join("source.fa");
    fs::write(
        &source,
        sequences
            .iter()
            .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    let mut mixed = a[..450].to_vec();
    mixed.extend_from_slice(&b[450..]);
    assert!(sequences.iter().all(|(_, s)| s != &mixed));
    let singleton = root.join("singletons.fa");
    let observed = oracle(&panel, &identity, &singleton, &mixed);
    let sample_path = root.join("sample.membwt");
    observed.save(&sample_path).unwrap();
    let native: Vec<_> = [a.clone(), b.clone(), mixed.clone()]
        .iter()
        .enumerate()
        .map(|(i, s)| oracle(&panel, &identity, &root.join(format!("oracle{i}.fa")), s))
        .collect();
    // Registry deliberately contains only original A definitions. Candidate junction
    // and B features must enter even when observed zero; registry IDs stay intact.
    let groups = root.join("groups.json");
    write(
        &groups,
        &json!({"version":1,"groups":[{"id":"a","scaffold":null,"occurrences":[{"path":"A#0#chr","start":0,"end":900,"strand":"+"}]}]}),
    );
    let catalog = root.join("catalog");
    run(
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
    let catalog_path = catalog.join("catalog.json");
    let registry = genome::load_catalog(&catalog_path).unwrap();
    let sources: Vec<_> = panel
        .name_map
        .path_to_name
        .iter()
        .enumerate()
        .map(|(id, name)| json!({"id":id,"name":name,"length":panel.name_map.path_to_length[id]}))
        .collect();
    let source_id = |name: &str| {
        panel
            .name_map
            .path_to_name
            .iter()
            .position(|n| n == name)
            .unwrap()
    };
    let instances:Vec<_>=(0..2).flat_map(|copy|["A","B"].into_iter().map(move |name|(copy,name))).map(|(copy,name)|json!({"id":format!("{name}{copy}"),"source":source_id(&format!("{name}#0#chr")),"start":0,"end":900})).collect();
    let layout = json!({"version":1,"model":joint::LAYOUT_MODEL,"panel":identity,"sources":sources,"instances":instances,
        "slots":[{"id":"slot0","alternatives":[walk("a",vec![piece("A0",0,900)]),walk("b",vec![piece("B0",0,900)]),
            walk("mosaic",vec![piece("A0",0,450),piece("B0",450,900)]),
            walk("same-physical-description-cut",vec![piece("A0",0,200),piece("A0",200,450),piece("B0",450,900)])]},
            {"id":"slot1","alternatives":[walk("a-copy",vec![piece("A1",0,900)]),walk("b-copy",vec![piece("B1",0,900)])]}]});
    let layout_path = root.join("layout.json");
    write(&layout_path, &layout);
    let compiled_dir = root.join("compiled");
    run(
        &[
            "genome-infer",
            "compile-joint-walks",
            "--panel",
            p(&prefix),
            "--layout",
            p(&layout_path),
            "--sources",
            p(&source),
            "--read-lengths",
            "250,1000",
            "--registry-catalog",
            p(&catalog_path),
            "--out-dir",
            p(&compiled_dir),
        ],
        true,
    );
    let profiles = compiled_dir.join("joint-profiles.json");
    let (compiled, digest) = joint::Compiled::load(&profiles, &identity).unwrap();
    assert_eq!(
        serde_json::to_value(&compiled.profiles[0][2].contributions).unwrap(),
        serde_json::to_value(&compiled.profiles[0][3].contributions).unwrap()
    );
    assert!(compiled
        .profiles
        .iter()
        .flatten()
        .all(|p| p.admitted_starts[1] == 0));
    let registry_tokens: std::collections::BTreeSet<_> =
        registry.features.iter().map(|f| f.tokens.clone()).collect();
    let a_counts = native[0].counts.observed_pairs().unwrap();
    let b_counts = native[1].counts.observed_pairs().unwrap();
    let m_counts = native[2].counts.observed_pairs().unwrap();
    assert!(
        m_counts
            .keys()
            .any(|t| !registry_tokens.contains(&t.to_vec())
                && !a_counts.contains_key(t)
                && !b_counts.contains_key(t)),
        "genuine junction creates a feature absent from either complete donor read profile"
    );
    assert!(
        a_counts
            .iter()
            .any(|(t, q)| *q > m_counts.get(t).copied().unwrap_or(0)),
        "junction destroys donor exposure"
    );
    for d in &compiled.definitions {
        for &id in &d.original_ids {
            assert_eq!(registry.features[id].tokens, d.tokens);
        }
    }
    for (ai, &oracle_id) in [0usize, 1, 2, 2].iter().enumerate() {
        let q = native[oracle_id].counts.observed_pairs().unwrap();
        let actual: BTreeMap<_, _> = compiled.profiles[0][ai]
            .contributions
            .iter()
            .map(|c| (c.tokens, c.totals[0][0]))
            .collect();
        assert_eq!(actual, q);
    }
    let solution = root.join("solution");
    run(
        &[
            "genome-infer",
            "solve-joint-walks",
            "--panel",
            p(&prefix),
            "--compiled",
            p(&profiles),
            "--sample",
            p(&sample_path),
            "--haploid-depth",
            "250",
            "--max-assignments",
            "8",
            "--out-dir",
            p(&solution),
        ],
        true,
    );
    let result = read(&solution.join("result.json"));
    let factors = read(&solution.join("factors.json"));
    let rows = factors.as_array().unwrap();
    assert!(rows.iter().any(|f| f["observed"] == 0));
    let mut brute = Vec::new();
    for (ai, &oi) in [0usize, 1, 2, 2].iter().enumerate() {
        for bi in 0..2 {
            let mut expected = 0.0;
            for f in rows {
                let t: Vec<u64> = serde_json::from_value(f["tokens"].clone()).unwrap();
                let c = observed.counts.count(&t).unwrap() as f64;
                let mu = 0.1
                    + native[oi].counts.count(&t).unwrap() as f64
                    + native[bi].counts.count(&t).unwrap() as f64;
                expected += mu - c * mu.ln();
            }
            brute.push((vec![ai, bi], expected));
            let assignment = root.join(format!("assignment-{ai}-{bi}.json"));
            write(
                &assignment,
                &json!({"version":1,"model":joint::MODEL,"compiled_checksum":digest,"choices":[ai,bi]}),
            );
            let out = root.join(format!("eval-{ai}-{bi}"));
            run(
                &[
                    "genome-infer",
                    "evaluate-joint-walks",
                    "--panel",
                    p(&prefix),
                    "--compiled",
                    p(&profiles),
                    "--sample",
                    p(&sample_path),
                    "--haploid-depth",
                    "250",
                    "--assignment",
                    p(&assignment),
                    "--out-dir",
                    p(&out),
                ],
                true,
            );
            assert!(
                (read(&out.join("result.json"))["result"]["objective"]
                    .as_f64()
                    .unwrap()
                    - expected)
                    .abs()
                    < 1e-7
            );
        }
    }
    let minimum = brute.iter().map(|(_, s)| *s).fold(f64::INFINITY, f64::min);
    assert!((result["result"]["incumbent"]["objective"].as_f64().unwrap() - minimum).abs() < 1e-7);
    let expected: Vec<_> = brute
        .iter()
        .filter(|(_, s)| s - minimum <= 1e-9)
        .map(|(a, _)| a.clone())
        .collect();
    let actual: Vec<Vec<usize>> = result["result"]["correlated_optima"]
        .as_array()
        .unwrap()
        .iter()
        .map(|a| serde_json::from_value(a["choices"].clone()).unwrap())
        .collect();
    assert_eq!(expected, actual);
    assert_eq!(result["result"]["global_objective_certified"], true);
    assert_eq!(result["result"]["sequence_emission_authorized"], false);
    let bounded = root.join("bounded");
    run(
        &[
            "genome-infer",
            "solve-joint-walks",
            "--panel",
            p(&prefix),
            "--compiled",
            p(&profiles),
            "--sample",
            p(&sample_path),
            "--haploid-depth",
            "250",
            "--max-assignments",
            "1",
            "--out-dir",
            p(&bounded),
        ],
        true,
    );
    assert_eq!(
        read(&bounded.join("result.json"))["result"]["status"],
        "budget-exhausted"
    );
    // New artifacts cannot flow through unary reconstruction/threading routes.
    let names = root.join("panel.syng.names");
    let result_path = solution.join("result.json");
    run(
        &[
            "genome-infer",
            "reconstruct",
            "--calls",
            p(&result_path),
            "--threads",
            p(&result_path),
            "--panel-names",
            p(&names),
            "--sources",
            p(&source),
            "--out-dir",
            p(&root.join("bad-reconstruction")),
        ],
        false,
    );
    let axis = root.join("axis.json");
    write(&axis, &json!({}));
    run(
        &[
            "genome-infer",
            "thread-partitions",
            "--panel",
            p(&prefix),
            "--observations",
            p(&compiled_dir),
            "--calls",
            p(&result_path),
            "--axis",
            p(&axis),
            "--out-dir",
            p(&root.join("bad-thread")),
        ],
        false,
    );
    let mut bad = layout.clone();
    bad["slots"][0]["alternatives"][0]["left_endpoint"] = json!("unknown");
    write(&layout_path, &bad);
    run(
        &[
            "genome-infer",
            "compile-joint-walks",
            "--panel",
            p(&prefix),
            "--layout",
            p(&layout_path),
            "--sources",
            p(&source),
            "--read-lengths",
            "250",
            "--out-dir",
            p(&root.join("bad-layout")),
        ],
        false,
    );
}
