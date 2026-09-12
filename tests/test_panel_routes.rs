//! Portable native/CLI exhaustive regression plus explicitly configured frozen-f7
//! and parent-Linux driver gates. The finite engine is independent of route replay.
mod panel_route_support;
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::{SyncmerParams, SyngIndex},
};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
    process::Command,
};
fn p(p: &Path) -> &str {
    p.to_str().unwrap()
}
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
fn write(p: &Path, v: &Value) {
    fs::write(p, serde_json::to_vec_pretty(v).unwrap()).unwrap();
}
fn read(p: &Path) -> Value {
    genome::read_json(p).unwrap()
}
fn run(binary: &str, args: &[&str], ok: bool) {
    let out = Command::new(binary).args(args).output().unwrap();
    assert_eq!(
        out.status.success(),
        ok,
        "{binary} {args:?}\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
}
#[derive(Clone)]
struct Port {
    word: Vec<u8>,
    source: usize,
    cut: u64,
    reverse: bool,
}
fn ports(root: &Path, g: &routes::Graph) -> Vec<Port> {
    let bytes = fs::read(root.join(&g.ports.path)).unwrap();
    bytes
        .chunks_exact(g.k as usize + 17)
        .map(|b| {
            let k = g.k as usize;
            let source = u64::from_le_bytes(b[k..k + 8].try_into().unwrap()) as usize;
            let anchor = u64::from_le_bytes(b[k + 8..k + 16].try_into().unwrap());
            let reverse = b[k + 16] == 1;
            Port {
                word: b[..k].to_vec(),
                source,
                cut: anchor + if reverse { g.k - g.k / 2 } else { g.k / 2 },
                reverse,
            }
        })
        .collect()
}
fn independent_routes(g: &routes::Graph, ports: &[Port], target: usize) -> Vec<routes::Route> {
    fn valid(segments: &[routes::Segment], p: &routes::Segment) -> bool {
        segments
            .iter()
            .all(|s| s.source != p.source || s.end <= p.start || p.end <= s.start)
    }
    fn visit(
        g: &routes::Graph,
        ports: &[Port],
        target: usize,
        source: usize,
        cut: u64,
        reverse: bool,
        segments: &mut Vec<routes::Segment>,
        out: &mut BTreeSet<routes::Route>,
        visits: &mut usize,
    ) {
        *visits += 1;
        assert!(*visits < 500000, "fixture is not a bounded oracle");
        if source == target && !reverse && cut < g.lanes[target].length {
            let tail = routes::Segment {
                source,
                start: cut,
                end: g.lanes[target].length,
                reverse: false,
            };
            if valid(segments, &tail) {
                let mut result = segments.clone();
                result.push(tail);
                out.insert(routes::Route { segments: result });
            }
        }
        for port in ports.iter().filter(|p| {
            p.source == source
                && p.reverse == reverse
                && if reverse { p.cut < cut } else { p.cut > cut }
        }) {
            let piece = routes::Segment {
                source,
                start: cut.min(port.cut),
                end: cut.max(port.cut),
                reverse,
            };
            if !valid(segments, &piece) {
                continue;
            }
            for donor in ports.iter().filter(|p| {
                p.word == port.word
                    && !(p.source == source && p.cut == port.cut && p.reverse == reverse)
            }) {
                segments.push(piece.clone());
                visit(
                    g,
                    ports,
                    target,
                    donor.source,
                    donor.cut,
                    donor.reverse,
                    segments,
                    out,
                    visits,
                );
                segments.pop();
            }
        }
    }
    let mut out = BTreeSet::new();
    visit(
        g,
        ports,
        target,
        target,
        0,
        false,
        &mut Vec::new(),
        &mut out,
        &mut 0,
    );
    out.into_iter().collect()
}
#[test]
fn automatic_cli_successive_donors_matches_portable_native_union() {
    route_fixture(env!("CARGO_BIN_EXE_impg"), false);
}

#[test]
#[ignore = "explicit Linux frozen-f7 gate; requires IMPG_TEST_F7_FINITE"]
fn frozen_panel_routes_f7_oracle() {
    route_fixture(&panel_route_support::configured_f7(), false);
}

#[test]
#[ignore = "explicit parent Linux driver gate; requires IMPG_TEST_PANEL_ROUTE_DRIVER=1 and IMPG_TEST_F7_FINITE"]
fn linux_parent_panel_route_driver() {
    assert_eq!(
        std::env::var("IMPG_TEST_PANEL_ROUTE_DRIVER").as_deref(),
        Ok("1"),
        "explicit driver gate requires IMPG_TEST_PANEL_ROUTE_DRIVER=1"
    );
    route_fixture(&panel_route_support::configured_f7(), true);
}

fn route_fixture(f7: &str, parent_driver: bool) {
    let t = tempfile::tempdir().unwrap();
    let preserved = if parent_driver {
        std::env::var_os("IMPG_TEST_PANEL_ROUTE_OUTPUT").map(std::path::PathBuf::from)
    } else {
        None
    };
    if let Some(path) = &preserved {
        assert!(
            path.is_absolute(),
            "configured fixture output must be absolute"
        );
        fs::create_dir(path).expect("configured fixture output must be fresh; retain failures");
        eprintln!("preserved synthetic driver fixture: {}", path.display());
    }
    let root = preserved.as_deref().unwrap_or_else(|| t.path());
    let binary = env!("CARGO_BIN_EXE_impg");
    let a = dna(400, 113);
    let mut b = a.clone();
    for pos in [150, 260] {
        b[pos] = if b[pos] == b'A' { b'C' } else { b'A' };
    }
    let c = dna(120, 59);
    let sequences = vec![
        ("A#0#native".to_string(), a.clone()),
        ("B#0#other-label".to_string(), b.clone()),
        ("C#0#short".to_string(), c.clone()),
    ];
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
        &json!({"version":1,"groups":sequences.iter().enumerate().map(|(i,(n,s))|json!({"id":format!("not-reference-chain-{i}"),"scaffold":null,"occurrences":[{"path":n,"start":0,"end":s.len(),"strand":null}]})).collect::<Vec<_>>()}),
    );
    let catalog = root.join("catalog");
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
    let graph_dir = root.join("routes");
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
            "150,500",
            "--core-bp",
            "37",
            "--out-dir",
            p(&graph_dir),
        ],
        true,
    );
    let graph = routes::Graph::load(&graph_dir, &identity).unwrap();
    assert_eq!(graph.lanes.len(), 3);
    assert_eq!(graph.families.len(), 3);
    let ports = ports(&graph_dir, &graph);
    let mut all = Vec::new();
    for family in 0..3 {
        let target = graph.families[family].paths[0];
        for route in independent_routes(&graph, &ports, target) {
            all.push((family, route));
        }
    }
    assert!(
        all.iter().any(|(_, r)| r.segments.len() >= 5),
        "fixture needs successive donor switches; {} routes, {} ports",
        all.len(),
        ports.len()
    );
    assert!(
        all.len() < 20000,
        "bounded oracle route count {}",
        all.len()
    );
    eprintln!(
        "Independent bounded graph: {} complete routes, {} oriented ports",
        all.len(),
        ports.len()
    );
    let mut mosaic = a.clone();
    mosaic[150] = b[150];
    assert!(sequences.iter().all(|(_, s)| s != &mosaic));
    let reads = root.join("reads.fa");
    let mut text = mosaic
        .windows(150)
        .map(|s| format!(">singleton\n{}\n", String::from_utf8_lossy(s)))
        .collect::<String>();
    let mut residual = c.clone();
    residual.extend(vec![b'N'; 30]);
    text += &format!(">sample-only\n{}\n", String::from_utf8_lossy(&residual));
    fs::write(&reads, text).unwrap();
    let sample = sample::build(&panel, identity.clone(), &[reads]).unwrap();
    let sample_path = root.join("sample.membwt");
    sample.save(&sample_path).unwrap();
    // F7 gets the COMPLETE union of all independently enumerated feasible routes,
    // not an incumbent-dependent feature set. One explicit slot per one-path family
    // alternative is sufficient for this bounded independent provider/search oracle.
    let layout = root.join("finite-layout.json");
    let alternatives:Vec<_>=all.iter().enumerate().map(|(i,(_,r))|json!({"id":format!("route-{i}"),"topology":"linear","left_endpoint":"asserted-molecule-terminus","right_endpoint":"asserted-molecule-terminus","pieces":r.segments.iter().map(|s|json!({"instance":format!("source-{}",s.source),"start":s.start,"end":s.end,"strand":if s.reverse{"-"}else{"+"}})).collect::<Vec<_>>(),"adjacencies":(0..r.segments.len()-1).map(|j|json!({"from":j,"to":j+1,"kind":"abut"})).collect::<Vec<_>>()})).collect();
    write(
        &layout,
        &json!({"version":1,"model":"explicit-physical-linear-layout-v1","panel":identity,"sources":graph.lanes.iter().map(|l|json!({"id":l.id,"name":l.name,"length":l.length})).collect::<Vec<_>>(),"instances":graph.lanes.iter().map(|l|json!({"id":format!("source-{}",l.id),"source":l.id,"start":0,"end":l.length})).collect::<Vec<_>>(),"slots":[{"id":"all-bounded-routes","alternatives":alternatives}]}),
    );
    let finite = root.join("finite");
    run(
        &f7,
        &[
            "genome-infer",
            "compile-joint-walks",
            "--panel",
            p(&prefix),
            "--layout",
            p(&layout),
            "--sources",
            p(&fasta),
            "--read-lengths",
            "150,500",
            "--registry-catalog",
            p(&catalog.join("catalog.json")),
            "--out-dir",
            p(&finite),
        ],
        true,
    );
    let compiled = read(&finite.join("joint-profiles.json"));
    let payload = &compiled["payload"];
    let mut union: BTreeSet<[u64; 3]> = sample
        .counts
        .observed_pairs()
        .unwrap()
        .keys()
        .copied()
        .collect();
    for d in payload["definitions"].as_array().unwrap() {
        union.insert(serde_json::from_value(d["tokens"].clone()).unwrap());
    }
    let beta = 0.1f64;
    let background: f64 = union
        .iter()
        .map(|t| beta - sample.counts.count(t).unwrap() as f64 * beta.ln())
        .sum();
    let mut evaluator = routes::Evaluator::new(
        &graph_dir, &graph, &panel, &sample, 150.0, beta, 1000000, 100000,
    )
    .unwrap();
    let mut expected = BTreeMap::new();
    let mut generated_zero = false;
    let mut positive_residual = false;
    for (i, (family, route)) in all.iter().enumerate() {
        let assignment = routes::Assignment {
            version: 1,
            model: routes::MODEL.into(),
            graph_checksum: graph.digest().unwrap(),
            family: *family,
            routes: vec![route.clone()],
        };
        let actual = evaluator.evaluate(&assignment).unwrap();
        let contributions = &payload["profiles"][0][i]["contributions"];
        let native: BTreeMap<[u64; 3], Vec<[u64; 2]>> = contributions
            .as_array()
            .unwrap()
            .iter()
            .map(|c| {
                (
                    serde_json::from_value(c["tokens"].clone()).unwrap(),
                    serde_json::from_value(c["totals"].clone()).unwrap(),
                )
            })
            .collect();
        let absolute: f64 = union
            .iter()
            .map(|t| {
                let q = native.get(t).map_or(0, |q| q[0][0]);
                let mu = beta + q as f64;
                mu - sample.counts.count(t).unwrap() as f64 * mu.ln()
            })
            .sum();
        assert!(
            (actual.relative_objective - (absolute - background)).abs() < 1e-7,
            "route {i}"
        );
        for f in &actual.factors {
            let q = native.get(&f.tokens).cloned().unwrap_or(vec![[0, 0]; 2]);
            assert_eq!(f.counts_by_length, vec![q[0][0], q[1][0]]);
            if f.observed == 0 && f.signal > 0.0 {
                generated_zero = true;
                assert_eq!(f.relative_loss, f.signal);
            }
            positive_residual |= f.positive_residual;
        }
        assert!(evaluator.lower_bound().unwrap() <= actual.relative_objective);
        expected.insert(
            serde_json::to_string(&actual.assignment).unwrap(),
            actual.relative_objective,
        );
    }
    assert!(generated_zero && positive_residual);
    let searched = root.join("searched");
    run(
        binary,
        &[
            "genome-infer",
            "search-panel-routes",
            "--panel",
            p(&prefix),
            "--routes",
            p(&graph_dir),
            "--sample",
            p(&sample_path),
            "--haploid-depth",
            "150",
            "--max-work",
            "5000000",
            "--max-evaluations",
            "100000",
            "--max-frontier",
            "1000",
            "--max-optima",
            "20000",
            "--out-dir",
            p(&searched),
        ],
        true,
    );
    let result = read(&searched.join("result.json"));
    let r = &result["result"];
    assert_eq!(r["native_initializations_complete"], true);
    assert_eq!(r["search_exhausted"], true, "{r}");
    assert_eq!(r["global_optimum_certified"], true);
    assert_eq!(r["correlated_optima_complete"], true);
    assert!(r["mixed_assignments_evaluated"].as_u64().unwrap() > 0);
    assert!(r["mixed_identity_assignments_evaluated"].as_u64().unwrap() > 0);
    assert!(
        r["non_native_assignments_evaluated"].as_u64().unwrap()
            >= r["mixed_assignments_evaluated"].as_u64().unwrap()
    );
    assert!(r["maximum_switches_evaluated"].as_u64().unwrap() >= 4);
    assert_eq!(r["sequence_emission_authorized"], false);
    let visited: BTreeSet<_> = fs::read_to_string(searched.join("evaluations.jsonl"))
        .unwrap()
        .lines()
        .map(|l| {
            let row: Value = serde_json::from_str(l).unwrap();
            let a: routes::Assignment =
                serde_json::from_value(row["scored"]["assignment"].clone()).unwrap();
            serde_json::to_string(&a).unwrap()
        })
        .collect();
    assert_eq!(visited, expected.keys().cloned().collect());
    let minimum = expected.values().copied().fold(f64::INFINITY, f64::min);
    assert!((r["incumbent"]["relative_objective"].as_f64().unwrap() - minimum).abs() < 1e-7);
    let ties: BTreeSet<_> = expected
        .iter()
        .filter(|(_, s)| **s - minimum <= 1e-9)
        .map(|(a, _)| a.clone())
        .collect();
    let actual: BTreeSet<_> = r["correlated_optima"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| {
            let a: routes::Assignment = serde_json::from_value(s["assignment"].clone()).unwrap();
            serde_json::to_string(&a).unwrap()
        })
        .collect();
    assert_eq!(actual, ties);
    // Frozen executable itself evaluates the full union and performs exact search.
    let finite_solved = root.join("finite-solved");
    run(
        &f7,
        &[
            "genome-infer",
            "solve-joint-walks",
            "--panel",
            p(&prefix),
            "--compiled",
            p(&finite.join("joint-profiles.json")),
            "--sample",
            p(&sample_path),
            "--haploid-depth",
            "150",
            "--max-assignments",
            "100000",
            "--max-optima",
            "20000",
            "--out-dir",
            p(&finite_solved),
        ],
        true,
    );
    let f7_result = read(&finite_solved.join("result.json"));
    assert!(
        (f7_result["result"]["incumbent"]["objective"]
            .as_f64()
            .unwrap()
            - background
            - minimum)
            .abs()
            < 1e-7
    );
    let evaluation = root.join("evaluated");
    run(
        binary,
        &[
            "genome-infer",
            "evaluate-panel-routes",
            "--panel",
            p(&prefix),
            "--routes",
            p(&graph_dir),
            "--sample",
            p(&sample_path),
            "--haploid-depth",
            "150",
            "--assignment",
            p(&searched.join("incumbent-assignment.json")),
            "--out-dir",
            p(&evaluation),
        ],
        true,
    );
    assert!(
        (read(&evaluation.join("result.json"))["result"]["relative_objective"]
            .as_f64()
            .unwrap()
            - minimum)
            .abs()
            < 1e-7
    );
    let capped = root.join("capped-optima");
    run(
        binary,
        &[
            "genome-infer",
            "search-panel-routes",
            "--panel",
            p(&prefix),
            "--routes",
            p(&graph_dir),
            "--sample",
            p(&sample_path),
            "--haploid-depth",
            "150",
            "--max-work",
            "5000000",
            "--max-evaluations",
            "100000",
            "--max-optima",
            "1",
            "--out-dir",
            p(&capped),
        ],
        true,
    );
    let cap_result = read(&capped.join("result.json"));
    assert_eq!(cap_result["result"]["global_optimum_certified"], true);
    assert_eq!(cap_result["result"]["correlated_optima_complete"], false);
    let bounded = root.join("bounded");
    run(
        binary,
        &[
            "genome-infer",
            "search-panel-routes",
            "--panel",
            p(&prefix),
            "--routes",
            p(&graph_dir),
            "--sample",
            p(&sample_path),
            "--haploid-depth",
            "150",
            "--max-work",
            "4",
            "--max-optima",
            "1",
            "--out-dir",
            p(&bounded),
        ],
        true,
    );
    assert_eq!(
        read(&bounded.join("result.json"))["result"]["global_optimum_certified"],
        false
    );
    run(
        binary,
        &[
            "genome-infer",
            "reconstruct",
            "--calls",
            p(&searched.join("result.json")),
            "--threads",
            p(&searched.join("result.json")),
            "--panel-names",
            p(&root.join("panel.syng.names")),
            "--sources",
            p(&fasta),
            "--out-dir",
            p(&root.join("bad-reconstruct")),
        ],
        false,
    );
    // Persisted native shards are checked, not trusted solely because graph.json
    // exists. Retry must still fail after an unsuccessful seal verification.
    let shard = graph_dir.join(&graph.lanes[0].profiles[0].totals.path);
    let before = fs::read(&shard).unwrap();
    fs::write(&shard, b"corrupt\n").unwrap();
    let mut bad =
        routes::Evaluator::new(&graph_dir, &graph, &panel, &sample, 150.0, beta, 1000000, 0)
            .unwrap();
    for _ in 0..2 {
        assert!(bad
            .evaluate(&graph.native_assignment(graph.lanes[0].family).unwrap())
            .is_err());
    }
    fs::write(shard, before).unwrap();
    if !parent_driver {
        return;
    }
    // Execute the parent driver only on this synthetic fixture. It independently
    // checks inventory, ownership, view/copy accounting and fresh output linkage.
    let inventory = root.join("inventory.json");
    write(
        &inventory,
        &json!({"path_count":graph.lanes.len(),"source_bp":graph.total_source_bp,"family_count":graph.families.len(),"zero_length_paths":[],"families":graph.families.iter().map(|f|json!({"identity":f.identity,"paths":f.paths.iter().map(|&s|json!({"source":s,"path":graph.lanes[s].name,"length":graph.lanes[s].length})).collect::<Vec<_>>()})).collect::<Vec<_>>()}),
    );
    let hash = |path: &Path| {
        let o = Command::new("sha256sum").arg(path).output().unwrap();
        assert!(o.status.success());
        String::from_utf8(o.stdout)
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .to_string()
    };
    let mut paths = vec![
        root.join("catalog/catalog.json"),
        inventory.clone(),
        Path::new(binary).to_path_buf(),
        Path::new(&f7).to_path_buf(),
        fasta.clone(),
        sample_path.clone(),
    ];
    for suffix in ["1gbwt", "1khash", "names", "meta", "spos", "pstep"] {
        paths.push(root.join(format!("panel.syng.{suffix}")));
    }
    paths.extend(
        graph
            .source_access
            .iter()
            .map(|b| Path::new(&b.path).to_path_buf()),
    );
    let frozen = root.join("frozen-inputs.json");
    write(
        &frozen,
        &serde_json::to_value(
            paths
                .iter()
                .map(|p| {
                    (
                        fs::canonicalize(p).unwrap().to_string_lossy().into_owned(),
                        hash(p),
                    )
                })
                .collect::<BTreeMap<_, _>>(),
        )
        .unwrap(),
    );
    let driver = root.join("parent-driver-synthetic");
    let binary_hash = hash(Path::new(binary));
    let catalog_path = root.join("catalog/catalog.json");
    let driver_args = [
        "scripts/panel-route-gate.py",
        "--execute",
        "--impg",
        binary,
        "--expected-binary-sha256",
        &binary_hash,
        "--f7-finite",
        &f7,
        "--panel",
        p(&prefix),
        "--catalog",
        p(&catalog_path),
        "--sources",
        p(&fasta),
        "--inventory",
        p(&inventory),
        "--frozen-inputs",
        p(&frozen),
        "--sample",
        "synthetic",
        p(&sample_path),
        "150",
        "--read-lengths",
        "150",
        "500",
        "--max-work",
        "5000000",
        "--max-evaluations",
        "100000",
        "--max-optima",
        "20000",
        "--out-dir",
        p(&driver),
    ];
    run("python3", &driver_args, true);
    let gate = read(&driver.join("generation-gate.json"));
    assert_eq!(gate["native_profiles_verified"], 6);
    assert_eq!(gate["native_starts_verified"], 502);
    assert_eq!(gate["native_totals_independently_verified"], true);
    assert!(gate["native_event_runs_verified"].as_u64().unwrap() > 0);
    eprintln!("independent generation gate: {gate}");

    let missing = root.join("frozen-inputs-without-fai.json");
    let mut declarations = read(&frozen);
    declarations.as_object_mut().unwrap().remove(
        &fs::canonicalize(&graph.source_access[0].path)
            .unwrap()
            .to_string_lossy()
            .into_owned(),
    );
    write(&missing, &declarations);
    let failed = root.join("parent-driver-missing-fai-declaration");
    let mut reject_args = driver_args;
    let i = reject_args
        .iter()
        .position(|s| *s == "--frozen-inputs")
        .unwrap();
    reject_args[i + 1] = p(&missing);
    let i = reject_args.iter().position(|s| *s == "--out-dir").unwrap();
    reject_args[i + 1] = p(&failed);
    let rejected = Command::new("python3").args(reject_args).output().unwrap();
    fs::write(
        root.join("missing-fai-declaration.stdout.log"),
        &rejected.stdout,
    )
    .unwrap();
    fs::write(
        root.join("missing-fai-declaration.stderr.log"),
        &rejected.stderr,
    )
    .unwrap();
    assert!(!rejected.status.success());
    assert!(
        String::from_utf8_lossy(&rejected.stderr).contains("missing frozen input SHA declaration")
    );
    assert_eq!(
        read(&failed.join("driver-status.json"))["stages"]
            .as_array()
            .unwrap()
            .len(),
        0
    );
    assert_eq!(
        read(&driver.join("driver-status.json"))["status"],
        "succeeded-operational-gates-only"
    );
    assert_eq!(
        read(&driver.join("generation-gate.json"))["native_paths"],
        3
    );
}
