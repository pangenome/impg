//! A bounded multi-molecule exchange against the independent finite native engine
//! (current build by default, explicitly configured frozen f7 in the ignored gate).
//! Atomic physical instances exactly match the finer canonical-span capacity.
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
fn write(p: &Path, v: Value) {
    fs::write(p, serde_json::to_vec_pretty(&v).unwrap()).unwrap()
}
fn read(p: &Path) -> Value {
    genome::read_json(p).unwrap()
}
fn run(bin: &str, args: &[&str]) {
    let o = Command::new(bin).args(args).output().unwrap();
    assert!(
        o.status.success(),
        "{args:?}\n{}",
        String::from_utf8_lossy(&o.stderr)
    );
}
#[test]
fn coupled_exchange_search_matches_portable_native_atomic_resource_oracle() {
    coupling_fixture(env!("CARGO_BIN_EXE_impg"));
}

#[test]
#[ignore = "explicit Linux frozen-f7 gate; requires IMPG_TEST_F7_FINITE"]
fn frozen_panel_route_coupled_f7_oracle() {
    coupling_fixture(&panel_route_support::configured_f7());
}

fn coupling_fixture(f7: &str) {
    let t = tempfile::tempdir().unwrap();
    let root = t.path();
    let bin = env!("CARGO_BIN_EXE_impg");
    let sequence = dna(220, 113);
    let sequences = vec![
        ("A#0#one".to_string(), sequence.clone()),
        ("A#0#two".to_string(), sequence.clone()),
    ];
    let mut panel = SyngIndex::build(SyncmerParams::default(), sequences.clone().into_iter());
    let prefix = root.join("p.syng");
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
        json!({"version":1,"groups":sequences.iter().enumerate().map(|(i,(name,_))|json!({"id":format!("g{i}"),"scaffold":null,"occurrences":[{"path":name,"start":0,"end":220,"strand":null}]})).collect::<Vec<_>>()}),
    );
    let catalog = root.join("catalog");
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
    );
    let out = root.join("routes");
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
            "220",
            "--core-bp",
            "31",
            "--out-dir",
            p(&out),
        ],
    );
    let g = routes::Graph::load(&out, &identity).unwrap();
    assert_eq!(g.families.len(), 1);
    // Random identical lanes have the same forward cuts and no reverse-equal
    // words. Enumerate all binary lane paths between these ordered cuts directly,
    // independently of the production DFS/hub iteration.
    let bytes = fs::read(out.join(&g.lanes[0].ports.path)).unwrap();
    let k = g.k as usize;
    let mut cuts = Vec::new();
    let mut forward = BTreeSet::new();
    let mut reverse = BTreeSet::new();
    for row in bytes.chunks_exact(k + 17) {
        let anchor = u64::from_le_bytes(row[k + 8..k + 16].try_into().unwrap());
        if row[k + 16] == 0 {
            cuts.push(anchor + g.cut_offset);
            forward.insert(row[..k].to_vec());
        } else {
            reverse.insert(row[..k].to_vec());
        }
    }
    assert!(forward.is_disjoint(&reverse));
    assert!(cuts.len() >= 2 && cuts.len() <= 10);
    let mut all = Vec::new();
    for target in 0..2 {
        let mut rows = Vec::new();
        for mask in 0..(1usize << cuts.len()) {
            let mut source = target;
            let mut start = 0;
            let mut pieces = Vec::new();
            for (i, &cut) in cuts.iter().enumerate() {
                let next = (mask >> i) & 1;
                if next != source {
                    pieces.push(routes::Segment {
                        source,
                        start,
                        end: cut,
                        reverse: false,
                    });
                    source = next;
                    start = cut;
                }
            }
            if source == target {
                pieces.push(routes::Segment {
                    source,
                    start,
                    end: 220,
                    reverse: false,
                });
                rows.push(routes::Route { segments: pieces });
            }
        }
        rows.sort();
        rows.dedup();
        all.push(rows);
    }
    let mut boundaries = vec![0];
    boundaries.extend(&cuts);
    boundaries.push(220);
    let atoms: Vec<_> = (0..2)
        .flat_map(|source| boundaries.windows(2).map(move |w| (source, w[0], w[1])))
        .collect();
    let instances: Vec<_> = atoms
        .iter()
        .map(|&(s, a, b)| json!({"id":format!("{s}:{a}:{b}"),"source":s,"start":a,"end":b}))
        .collect();
    let slots:Vec<_>=all.iter().enumerate().map(|(slot,rows)|json!({"id":format!("slot{slot}"),"alternatives":rows.iter().enumerate().map(|(i,r)|{let pieces:Vec<_>=r.segments.iter().flat_map(|p|atoms.iter().filter(move |&&(s,a,b)|s==p.source && a>=p.start && b<=p.end).map(|&(s,a,b)|json!({"instance":format!("{s}:{a}:{b}"),"start":a,"end":b,"strand":"+"}))).collect();json!({"id":format!("a{i}"),"topology":"linear","left_endpoint":"asserted-molecule-terminus","right_endpoint":"asserted-molecule-terminus","adjacencies":(0..pieces.len()-1).map(|j|json!({"from":j,"to":j+1,"kind":"abut"})).collect::<Vec<_>>(),"pieces":pieces})}).collect::<Vec<_>>()})).collect();
    let layout = root.join("layout.json");
    write(
        &layout,
        json!({"version":1,"model":"explicit-physical-linear-layout-v1","panel":identity,"sources":g.lanes.iter().map(|l|json!({"id":l.id,"name":l.name,"length":l.length})).collect::<Vec<_>>(),"instances":instances,"slots":slots}),
    );
    let sample = sample::build(&panel, identity.clone(), &[fasta.clone()]).unwrap();
    let samplepath = root.join("sample.membwt");
    sample.save(&samplepath).unwrap();
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
            "220",
            "--out-dir",
            p(&finite),
        ],
    );
    let fsol = root.join("fsol");
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
            p(&samplepath),
            "--haploid-depth",
            "220",
            "--max-assignments",
            "2000000",
            "--max-optima",
            "20000",
            "--out-dir",
            p(&fsol),
        ],
    );
    let searched = root.join("searched");
    run(
        bin,
        &[
            "genome-infer",
            "search-panel-routes",
            "--panel",
            p(&prefix),
            "--routes",
            p(&out),
            "--sample",
            p(&samplepath),
            "--haploid-depth",
            "220",
            "--max-work",
            "5000000",
            "--max-evaluations",
            "2000000",
            "--max-optima",
            "20000",
            "--out-dir",
            p(&searched),
        ],
    );
    let r = read(&searched.join("result.json"));
    let f = read(&fsol.join("result.json"));
    assert_eq!(r["result"]["global_optimum_certified"], true);
    assert_eq!(r["result"]["correlated_optima_complete"], true);
    assert_eq!(f["result"]["global_objective_certified"], true);
    assert_eq!(f["result"]["correlated_optima_complete"], true);
    let factors = read(&fsol.join("factors.json"));
    let beta = 0.1f64;
    let background: f64 = factors
        .as_array()
        .unwrap()
        .iter()
        .map(|f| beta - f["observed"].as_u64().unwrap() as f64 * beta.ln())
        .sum();
    assert!(
        (r["result"]["incumbent"]["relative_objective"]
            .as_f64()
            .unwrap()
            - (f["result"]["incumbent"]["objective"].as_f64().unwrap() - background))
            .abs()
            < 1e-8
    );
    let lookup: Vec<BTreeMap<routes::Route, usize>> = all
        .iter()
        .map(|rows| {
            rows.iter()
                .cloned()
                .enumerate()
                .map(|(i, r)| (r, i))
                .collect()
        })
        .collect();
    let actual: BTreeSet<Vec<usize>> = r["result"]["correlated_optima"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| {
            let a: routes::Assignment = serde_json::from_value(v["assignment"].clone()).unwrap();
            a.routes
                .iter()
                .enumerate()
                .map(|(s, r)| lookup[s][r])
                .collect()
        })
        .collect();
    let expected: BTreeSet<Vec<usize>> = f["result"]["correlated_optima"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| serde_json::from_value(v["choices"].clone()).unwrap())
        .collect();
    assert_eq!(actual, expected);
    assert!(actual.len() > 1);
    let marginals: Vec<BTreeSet<usize>> = (0..2)
        .map(|s| actual.iter().map(|a| a[s]).collect())
        .collect();
    assert!(
        marginals[0].len() * marginals[1].len() > actual.len(),
        "marginal Cartesian product must include physically conflicting copy exchanges"
    );
    assert!(r["result"]["mixed_assignments_evaluated"].as_u64().unwrap() > 0);
    assert_eq!(
        r["result"]["mixed_identity_assignments_evaluated"], 0,
        "exchanging source paths within A#0 is not a donor-identity change"
    );
    eprintln!("Coupled atomic-span oracle: {} x {} route alternatives; {} correlated feasible optimum assignments",all[0].len(),all[1].len(),actual.len());
}
