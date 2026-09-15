use super::search;
use clap::Parser;
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::{SyncmerParams, SyngIndex},
};
use search::score::Pair;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{BufRead, Read},
    path::{Path, PathBuf},
    process::Command,
};
pub fn dna(n: usize, mut s: u64) -> Vec<u8> {
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
fn write(p: &Path, v: &impl serde::Serialize) {
    genome::write_json(p, v).unwrap()
}
fn fasta(p: &Path, rows: &[(String, Vec<u8>)]) {
    fs::write(
        p,
        rows.iter()
            .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap()
}
fn cli(root: &Path, name: &str, args: &[&str]) {
    let o = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(args)
        .output()
        .unwrap();
    fs::write(root.join(format!("{name}.stdout")), &o.stdout).unwrap();
    fs::write(root.join(format!("{name}.stderr")), &o.stderr).unwrap();
    assert!(
        o.status.success(),
        "{name}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
}
#[derive(serde::Serialize, serde::Deserialize)]
pub struct Fixture {
    pub root: PathBuf,
    pub sources: Vec<(String, Vec<u8>)>,
    pub truth: Vec<Vec<u8>>,
    pub kind: String,
}
impl Fixture {
    pub fn options(&self, name: &str) -> search::Options {
        search::Options::parse_from([
            "paired",
            "--panel",
            p(&self.root.join("panel.syng")),
            "--routes",
            p(&self.root.join("routes")),
            "--sample",
            p(&self.root.join(if self.kind == "near-native" {
                "sample/sample.membwt"
            } else {
                "sample.membwt"
            })),
            "--out-dir",
            p(&self.root.join(name)),
        ])
    }
}
pub fn fixture(root: &Path, kind: &str) -> Fixture {
    fs::create_dir(root).unwrap();
    write(
        &root.join("recipe-before-operator.json"),
        &json!({"kind":kind,"development_B1_reuse":true,"seed":219,"main_bp":1536,"mutated_windows":[[256,280],[704,728],[1216,1240]],"passive_seeds":[57,58],"passive_bp":400,"read_length":150,"stride":15,"append_last":true,"alternating_global_rc":true,"nominal_per_copy_depth":10,"background_once":0.1}),
    );
    let a = dna(1536, 219);
    let mut b = a.clone();
    for start in [256, 704, 1216] {
        for x in &mut b[start..start + 24] {
            *x = if *x == b'A' { b'C' } else { b'A' };
        }
    }
    let passive = dna(400, 57);
    let sources = vec![
        ("A#0#one".into(), a.clone()),
        ("B#0#one".into(), b.clone()),
        ("A#0#two".into(), passive.clone()),
        ("B#0#two".into(), dna(400, 58)),
    ];
    let mutate = |start| {
        let mut s = a.clone();
        s[start..start + 24].copy_from_slice(&b[start..start + 24]);
        s
    };
    let truth = match kind {
        "middle" => vec![mutate(704), passive.clone(), a, passive],
        "last" => vec![mutate(1216), passive.clone(), a, passive],
        "phase" => vec![mutate(704), passive.clone(), mutate(1216), passive],
        _ => panic!("unknown frozen recipe"),
    };
    let mut panel = SyngIndex::build(SyncmerParams::default(), sources.clone().into_iter());
    let prefix = root.join("panel.syng");
    panel.save(p(&prefix)).unwrap();
    fasta(&root.join("sources.fa"), &sources);
    write(
        &root.join("groups.json"),
        &json!({"version":1,"groups":sources.iter().enumerate().map(|(i,(n,s))|json!({"id":format!("g{i}"),"scaffold":null,"occurrences":[{"path":n,"start":0,"end":s.len(),"strand":null}]})).collect::<Vec<_>>()}),
    );
    cli(
        root,
        "catalog",
        &[
            "genome-infer",
            "build-catalog",
            "--panel",
            p(&prefix),
            "--groups",
            p(&root.join("groups.json")),
            "--out-dir",
            p(&root.join("catalog")),
        ],
    );
    cli(
        root,
        "routes",
        &[
            "genome-infer",
            "build-panel-routes",
            "--panel",
            p(&prefix),
            "--catalog",
            p(&root.join("catalog/catalog.json")),
            "--sources",
            p(&root.join("sources.fa")),
            "--read-lengths",
            "150",
            "--core-bp",
            "37",
            "--out-dir",
            p(&root.join("routes")),
        ],
    );
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    make_sample(root, "sample", &truth, &panel, &identity);
    let fixture = Fixture {
        root: root.into(),
        sources,
        truth,
        kind: kind.into(),
    };
    write(&root.join("fixture-description.json"), &fixture);
    fixture
}
/// Exact accepted finite near recipe/order/simulation, without its supplied inference domain.
/// Keep synchronized with test_panel_route_diploid_diagnostic::fixture("near").
pub fn near_native_fixture(root: &Path) -> Fixture {
    let kind = "near";
    fs::create_dir(root).unwrap();
    // Freeze recipe artifact before any operator evaluation or sample scoring.
    genome::write_json(&root.join("recipe-before-operator.json"),&json!({"kind":kind,"seed":113,"length":600,"positions":match kind {"heterozygote"=>vec![300],"near"=>vec![220,260],"far"=>vec![180,420],_=>vec![]},"read_length":150,"start_stride":15,"alternating_rc":true,"nominal_per_copy_depth":10,"background":0.1})).unwrap();
    let base = dna(600, 113);
    let positions = match kind {
        "heterozygote" => vec![300],
        "near" => vec![220, 260],
        "far" => vec![180, 420],
        _ => vec![],
    };
    let mut sequences = Vec::new();
    for mask in 0..(1usize << positions.len()) {
        let mut s = base.clone();
        for (bit, &pos) in positions.iter().enumerate() {
            if mask & (1 << bit) != 0 {
                s[pos] = if s[pos] == b'A' { b'C' } else { b'A' };
            }
        }
        sequences.push((format!("H{mask}#0#chr"), s));
    }
    let mut panel = SyngIndex::build(SyncmerParams::default(), sequences.clone().into_iter());
    let prefix = root.join("panel.syng");
    panel.save(p(&prefix)).unwrap();
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
    genome::write_json(&groups,&json!({"version":1,"groups":sequences.iter().enumerate().map(|(i,(n,s))|json!({"id":format!("source-{i}"),"scaffold":null,"occurrences":[{"path":n,"start":0,"end":s.len(),"strand":null}]})).collect::<Vec<_>>()})).unwrap();
    let catalog = root.join("catalog");
    cli(
        root,
        "catalog",
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
    let graph_dir = root.join("routes");
    cli(
        root,
        "routes",
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
            "37",
            "--out-dir",
            p(&graph_dir),
        ],
    );
    // Sample generation alone uses truth. Inference gets only public artifacts.
    let truth_indices = [0, sequences.len() - 1];
    let reads = root.join("reads.fa");
    let mut text = String::new();
    let mut ordinal = 0;
    for i in truth_indices {
        for start in (0..=450).step_by(15) {
            let mut read = sequences[i].1[start..start + 150].to_vec();
            if ordinal % 2 == 1 {
                read = impg::graph::reverse_complement(&read);
            }
            text += &format!(">r{ordinal}\n{}\n", String::from_utf8_lossy(&read));
            ordinal += 1;
        }
    }
    fs::write(&reads, text).unwrap();
    let sample_dir = root.join("sample");
    cli(
        root,
        "sample",
        &[
            "genome-infer",
            "build-sample",
            "--panel",
            p(&prefix),
            "--reads",
            p(&reads),
            "--out-dir",
            p(&sample_dir),
        ],
    );
    let fixture = Fixture {
        root: root.into(),
        truth: truth_indices
            .iter()
            .map(|&i| sequences[i].1.clone())
            .collect(),
        sources: sequences,
        kind: "near-native".into(),
    };
    write(&root.join("fixture-description.json"), &fixture);
    fixture
}

pub fn dosage_fixture(root: &Path) -> Fixture {
    fs::create_dir(root).unwrap();
    let starts = [256usize, 576, 896, 1216];
    write(
        &root.join("recipe.json"),
        &json!({"seed":9917,"main_bp":1536,"passive_bp":400,"tracts":starts.map(|s|[s,s+24]),"truth_masks":[10,12],"read_length":150,"stride":15,"alternating_rc":true}),
    );
    let a = dna(1536, 9917);
    let mut b = a.clone();
    for start in starts {
        for x in &mut b[start..start + 24] {
            *x = if *x == b'A' { b'C' } else { b'A' };
        }
    }
    let passive = dna(400, 57);
    let sources = vec![
        ("A#0#main".into(), a.clone()),
        ("B#0#main".into(), b.clone()),
        ("A#0#passive".into(), passive.clone()),
        ("B#0#passive".into(), passive.clone()),
    ];
    let mosaic = |mask: usize| {
        let mut sequence = a.clone();
        for (bit, start) in starts.into_iter().enumerate() {
            if mask & (1 << bit) != 0 {
                sequence[start..start + 24].copy_from_slice(&b[start..start + 24]);
            }
        }
        sequence
    };
    let truth = vec![mosaic(10), passive.clone(), mosaic(12), passive];
    let mut panel = SyngIndex::build(SyncmerParams::default(), sources.clone().into_iter());
    let prefix = root.join("panel.syng");
    panel.save(p(&prefix)).unwrap();
    fasta(&root.join("sources.fa"), &sources);
    write(
        &root.join("groups.json"),
        &json!({"version":1,"groups":sources.iter().enumerate().map(|(i,(name,sequence))|json!({"id":format!("source-{i}"),"scaffold":null,"occurrences":[{"path":name,"start":0,"end":sequence.len(),"strand":null}]})).collect::<Vec<_>>()}),
    );
    cli(
        root,
        "catalog",
        &[
            "genome-infer",
            "build-catalog",
            "--panel",
            p(&prefix),
            "--groups",
            p(&root.join("groups.json")),
            "--out-dir",
            p(&root.join("catalog")),
        ],
    );
    cli(
        root,
        "routes",
        &[
            "genome-infer",
            "build-panel-routes",
            "--panel",
            p(&prefix),
            "--catalog",
            p(&root.join("catalog/catalog.json")),
            "--sources",
            p(&root.join("sources.fa")),
            "--read-lengths",
            "150",
            "--core-bp",
            "37",
            "--out-dir",
            p(&root.join("routes")),
        ],
    );
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    make_sample(root, "sample", &truth, &panel, &identity);
    let fixture = Fixture {
        root: root.into(),
        sources,
        truth,
        kind: "dosage-mosaic".into(),
    };
    write(&root.join("fixture-description.json"), &fixture);
    fixture
}

pub fn make_sample(
    root: &Path,
    name: &str,
    truth: &[Vec<u8>],
    panel: &SyngIndex,
    id: &genome::PanelIdentity,
) {
    let mut reads = Vec::new();
    for s in truth {
        let last = s.len() - 150;
        let mut starts = (0..=last).step_by(15).collect::<Vec<_>>();
        if starts.last() != Some(&last) {
            starts.push(last);
        }
        for start in starts {
            let mut read = s[start..start + 150].to_vec();
            if reads.len() % 2 == 1 {
                read = impg::graph::reverse_complement(&read);
            }
            reads.push((format!("anonymous{}", reads.len()), read));
        }
    }
    let fa = root.join(format!("{name}.fa"));
    fasta(&fa, &reads);
    sample::build(panel, id.clone(), &[fa])
        .unwrap()
        .save(&root.join(format!("{name}.membwt")))
        .unwrap();
}
fn canonical(s: &[u8]) -> Vec<u8> {
    s.to_vec().min(impg::graph::reverse_complement(s))
}
fn spell(pair: &Pair, sources: &[(String, Vec<u8>)]) -> Vec<Vec<u8>> {
    pair.copies
        .iter()
        .flat_map(|a| &a.routes)
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
pub fn coverage(truth: &[Vec<u8>], query: &[Vec<u8>]) -> Value {
    let mut unused: Vec<_> = query.iter().map(|q| Some(canonical(q))).collect();
    let mut matched = 0;
    let mut copies = 0;
    for t in truth {
        if let Some(i) = unused
            .iter()
            .position(|q| q.as_ref() == Some(&canonical(t)))
        {
            unused[i] = None;
            matched += t.len();
            copies += 1;
        }
    }
    let tb: usize = truth.iter().map(Vec::len).sum();
    let qb: usize = query.iter().map(Vec::len).sum();
    json!({"exact_matched_molecules":copies,"exact_matched_bp":matched,"truth_bp":tb,"query_bp":qb,"full_truth_coverage":matched as f64/tb as f64,"full_query_coverage":matched as f64/qb as f64,"metric":"injective exact-whole-DNA/RC, NOT alignment or QV"})
}
/// Each pass retains one raw buffer and one parsed row, never a ledger backlog.
fn ledger_rows(p: &Path) -> impl Iterator<Item = (Value, usize, usize)> {
    let mut reader = std::io::BufReader::new(fs::File::open(p).unwrap());
    let mut line = String::new();
    std::iter::from_fn(move || {
        line.clear();
        let bytes = reader.read_line(&mut line).unwrap();
        (bytes > 0).then(|| (serde_json::from_str(&line).unwrap(), line.capacity(), bytes))
    })
}
fn independent_loss(q: &BTreeMap<[u64; 3], u64>, observed: &BTreeMap<[u64; 3], u64>) -> f64 {
    let mut keys: BTreeSet<_> = q.keys().copied().collect();
    keys.extend(observed.keys().copied());
    assert!(keys.len() <= 50_000);
    let absolute: f64 = keys
        .iter()
        .map(|t| {
            let s = 10.0 * q.get(t).copied().unwrap_or(0) as f64 / 150.0;
            let c = observed.get(t).copied().unwrap_or(0) as f64;
            (0.1 + s) - c * (0.1 + s).ln() - (0.1 - c * 0.1f64.ln())
        })
        .sum();
    let relative: f64 = keys
        .iter()
        .map(|t| {
            let s = 10.0 * q.get(t).copied().unwrap_or(0) as f64 / 150.0;
            s - observed.get(t).copied().unwrap_or(0) as f64 * (s / 0.1).ln_1p()
        })
        .sum();
    assert!((relative - absolute).abs() <= 1e-8 * (1.0 + absolute.abs()));
    absolute
}
fn sum(a: &BTreeMap<[u64; 3], u64>, b: &BTreeMap<[u64; 3], u64>) -> BTreeMap<[u64; 3], u64> {
    let mut out = a.clone();
    for (t, n) in b {
        let x = out.entry(*t).or_default();
        *x = x.checked_add(*n).unwrap();
    }
    out
}
/// Oracle constructs independent witnessed subsets AFTER inference.
pub(crate) fn oracle_at(
    e: &mut routes::Evaluator<'_>,
    starts: &[usize],
) -> Vec<routes::Assignment> {
    let mut ports = Vec::new();
    for source in [0, 1] {
        let mut f = e.ports.source_file(e.graph, source).unwrap();
        for _ in 0..e.graph.lanes[source].port_count {
            let mut word = vec![0; e.graph.k as usize];
            f.read_exact(&mut word).unwrap();
            let (mut s, mut a, mut r) = ([0; 8], [0; 8], [0; 1]);
            f.read_exact(&mut s).unwrap();
            f.read_exact(&mut a).unwrap();
            f.read_exact(&mut r).unwrap();
            if r[0] == 0 {
                ports.push((
                    u64::from_le_bytes(s) as usize,
                    u64::from_le_bytes(a) + e.graph.cut_offset,
                    word,
                ));
            }
        }
    }
    let mut shared = Vec::new();
    for p in &ports {
        for q in &ports {
            if p.0 == 0 && q.0 == 1 && p.1 == q.1 && p.2 == q.2 {
                shared.push(p.1);
            }
        }
    }
    let flanks: Vec<_> = starts
        .iter()
        .map(|&lo| {
            let lo = lo as u64;
            (
                *shared.iter().filter(|&&p| p <= lo).max().unwrap(),
                *shared.iter().filter(|&&p| p >= lo + 24).min().unwrap(),
            )
        })
        .collect();
    let native = e.graph.native_assignment(0).unwrap();
    let slot = native
        .routes
        .iter()
        .position(|r| r.segments[0].source == 0)
        .unwrap();
    let mut domain = Vec::new();
    for mask in 0..(1usize << starts.len()) {
        let mut a = native.clone();
        let mut pieces = Vec::new();
        let mut start = 0;
        for (bit, &(lo, hi)) in flanks.iter().enumerate() {
            if mask & (1 << bit) != 0 {
                pieces.push(routes::Segment {
                    source: 0,
                    start,
                    end: lo,
                    reverse: false,
                });
                pieces.push(routes::Segment {
                    source: 1,
                    start: lo,
                    end: hi,
                    reverse: false,
                });
                start = hi;
            }
        }
        pieces.push(routes::Segment {
            source: 0,
            start,
            end: 1536,
            reverse: false,
        });
        a.routes[slot] = routes::Route { segments: pieces };
        domain.push(e.validate_assignment(&a).unwrap());
    }
    domain.push(e.graph.native_assignment(1).unwrap());
    domain
}
fn oracle(e: &mut routes::Evaluator<'_>) -> Vec<routes::Assignment> {
    oracle_at(e, &[256, 704, 1216])
}
pub fn assessment_dir(inference: &Path) -> PathBuf {
    std::env::var_os("IMPG_PAIRED_ASSESSMENT")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            inference.with_file_name(format!(
                "{}-assessment",
                inference.file_name().unwrap().to_str().unwrap()
            ))
        })
}
pub fn assess(f: &Fixture, name: &str, result: &Value) -> Value {
    let inference = f.root.join(name);
    let out = assessment_dir(&inference);
    fs::create_dir(&out).unwrap();
    let files = [
        "scores.jsonl",
        "events.jsonl",
        "pending.json",
        "result.json",
        "accounting.jsonl",
        "options.json",
        "provenance.json",
        "status.json",
    ];
    let frozen:Vec<_>=files.iter().map(|n|json!({"file":n,"fingerprint":genome::reconstruction::fingerprint(&inference.join(n)).unwrap()})).collect();
    write(
        &out.join("mechanical-fnv-before-assessment.json"),
        &json!({"authoritative_SHA256":false,"files":frozen}),
    );
    let ledger_path = inference.join("scores.jsonl");
    let prefix = f.root.join("panel.syng");
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    let panel = SyngIndex::load(p(&prefix), Default::default()).unwrap();
    let graph = routes::Graph::load(&f.root.join("routes"), &identity).unwrap();
    let options: Value = genome::read_json(&inference.join("options.json")).unwrap();
    let sample =
        sample::SampleIndex::load(Path::new(options["sample"].as_str().unwrap()), &identity)
            .unwrap();
    let observed = sample.counts.observed_pairs().unwrap();
    let mut e = routes::Evaluator::new(
        &f.root.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        50_000,
        50_000,
    )
    .unwrap();
    let mut profiles: BTreeMap<routes::Assignment, BTreeMap<[u64; 3], u64>> = BTreeMap::new();
    let near = f.kind == "near-native";
    let mut audit = AuditBudget {
        metadata_bytes: json_charge(result) + 1024 * 1024,
        ..Default::default()
    };
    let mut audited = 0usize;
    let mut public_losses = Vec::new();
    audit.ledger_passes += 1;
    for (row, capacity, bytes) in ledger_rows(&ledger_path) {
        audit.ledger_row(&row, capacity, bytes);
        if row["kind"] != "paired_score" {
            continue;
        }
        let pair: Pair = serde_json::from_value(row["scored"]["pair"].clone()).unwrap();
        for a in &pair.copies {
            if !profiles.contains_key(a) {
                audit.profile();
                let ev = e.evaluate(a).unwrap();
                public_losses.push(ev.relative_objective);
                let q = ev
                    .factors
                    .into_iter()
                    .filter_map(|f| {
                        (f.counts_by_length[0] > 0).then_some((f.tokens, f.counts_by_length[0]))
                    })
                    .collect();
                audit.retain(&q);
                profiles.insert(a.clone(), q);
            }
        }
        let q = sum(&profiles[&pair.copies[0]], &profiles[&pair.copies[1]]);
        audit.formulas(false);
        let objective = independent_loss(&q, &observed);
        let reported = row["scored"]["objective"].as_f64().unwrap();
        assert!((objective - reported).abs() < 1e-8 * (1.0 + objective.abs()));
        let mut publicq = BTreeMap::new();
        let mut public_keys = BTreeSet::new();
        for factor in row["factors"].as_array().unwrap() {
            let t: [u64; 3] = serde_json::from_value(factor["tokens"].clone()).unwrap();
            assert!(public_keys.insert(t), "duplicate factor");
            let n = factor["counts"].as_u64().unwrap();
            assert_eq!(n, q.get(&t).copied().unwrap_or(0));
            assert_eq!(
                factor["observed"].as_u64().unwrap(),
                observed.get(&t).copied().unwrap_or(0)
            );
            if n > 0 {
                publicq.insert(t, n);
            }
        }
        assert_eq!(q, publicq);
        let expected_keys: BTreeSet<_> = q.keys().chain(observed.keys()).copied().collect();
        assert_eq!(
            public_keys, expected_keys,
            "complete realized/observed universe"
        );
        audited += 1;
    }
    let domain = if near {
        assert_eq!(graph.families.len(), 4);
        (0..4)
            .map(|i| graph.native_assignment(i).unwrap())
            .collect()
    } else {
        oracle(&mut e)
    };
    let reference_pairs = if near {
        [[0, 3], [1, 2]]
    } else {
        [[0, 6], [2, 4]]
    };
    write(&out.join("assessment-oracle-domain.json"), &domain);
    let mut oracle_profiles = Vec::new();
    for a in &domain {
        audit.profile();
        let ev = e.evaluate(a).unwrap();
        public_losses.push(ev.relative_objective);
        let q: BTreeMap<_, _> = ev
            .factors
            .into_iter()
            .filter_map(|f| {
                (f.counts_by_length[0] > 0).then_some((f.tokens, f.counts_by_length[0]))
            })
            .collect();
        audit.retain(&q);
        oracle_profiles.push(q);
    }
    let mut oracle_rows = Vec::new();
    for i in 0..domain.len() {
        for j in i..domain.len() {
            let q = sum(&oracle_profiles[i], &oracle_profiles[j]);
            audit.formulas(true);
            let loss = independent_loss(&q, &observed);
            audit.retain(&q);
            oracle_rows
                .push(json!({"pair":[i,j],"objective":loss,"counts":q.iter().collect::<Vec<_>>()}));
        }
    }
    assert_eq!(oracle_rows.len(), if near { 10 } else { 45 });
    write(
        &out.join(if near {
            "assessment-oracle10.json"
        } else {
            "assessment-oracle45.json"
        }),
        &oracle_rows,
    );
    let [cis_indices, trans_indices] = reference_pairs;
    let cis = sum(
        &oracle_profiles[cis_indices[0]],
        &oracle_profiles[cis_indices[1]],
    );
    let trans = sum(
        &oracle_profiles[trans_indices[0]],
        &oracle_profiles[trans_indices[1]],
    );
    if !near {
        assert_eq!(cis, trans, "actual operator far-tract equivalence required");
    }
    let reference_loss = |indices: [usize; 2]| {
        oracle_rows
            .iter()
            .find(|r| r["pair"] == json!(indices))
            .unwrap()["objective"]
            .as_f64()
            .unwrap()
    };
    let phase_score_separation = reference_loss(trans_indices) - reference_loss(cis_indices);
    let pairs: Vec<Pair> =
        serde_json::from_value(result["retained_correlated_best_found"].clone()).unwrap();
    let mut assessments = Vec::new();
    let mut matched = false;
    let mut dna_classes = BTreeSet::new();
    // Public BWT replay cache: full spelled forward DNA, separately verifies RC; multiplicity never cached away.
    let mut replay_cache: BTreeMap<Vec<u8>, BTreeMap<[u64; 3], u64>> = BTreeMap::new();
    let mut replay_reads = 0usize;
    let mut validations = 0;
    for pair in pairs.iter().chain(
        [
            Pair {
                copies: [
                    domain[cis_indices[0]].clone(),
                    domain[cis_indices[1]].clone(),
                ],
            },
            Pair {
                copies: [
                    domain[trans_indices[0]].clone(),
                    domain[trans_indices[1]].clone(),
                ],
            },
        ]
        .iter(),
    ) {
        let mut pooled = BTreeMap::new();
        for a in &pair.copies {
            e.validate_assignment(a).unwrap();
            validations += 1;
            for r in &a.routes {
                let sequence = e.sources.spell(r, 0, r.length().unwrap()).unwrap();
                if !replay_cache.contains_key(&sequence) {
                    let mut forward = None;
                    for reverse in [false, true] {
                        let dna = if reverse {
                            impg::graph::reverse_complement(&sequence)
                        } else {
                            sequence.clone()
                        };
                        let n = dna.len().saturating_sub(149);
                        assert!(replay_reads + n <= 100_000, "independent replay cap");
                        replay_reads += n;
                        let fa = out.join(format!(
                            "independent-replay-{}-{reverse}.fa",
                            replay_cache.len()
                        ));
                        fasta(
                            &fa,
                            &dna.windows(150)
                                .enumerate()
                                .map(|(i, s)| (format!("anonymous{i}"), s.to_vec()))
                                .collect::<Vec<_>>(),
                        );
                        let replay = sample::build(&panel, identity.clone(), &[fa])
                            .unwrap()
                            .counts
                            .observed_pairs()
                            .unwrap();
                        if let Some(q) = &forward {
                            assert_eq!(*q, replay);
                        } else {
                            forward = Some(replay);
                        }
                    }
                    let counts = forward.unwrap();
                    audit.retain(&counts);
                    replay_cache.insert(sequence.clone(), counts);
                }
                let (q, _) = e.route_counts(r).unwrap();
                assert_eq!(q[0], replay_cache[&sequence]);
                pooled = sum(&pooled, &replay_cache[&sequence]);
            }
        }
        let profile_counts = if let Some(i) = reference_pairs.iter().position(|indices| {
            pair.copies == [domain[indices[0]].clone(), domain[indices[1]].clone()]
        }) {
            if i == 0 {
                cis.clone()
            } else {
                trans.clone()
            }
        } else {
            sum(&profiles[&pair.copies[0]], &profiles[&pair.copies[1]])
        };
        assert_eq!(
            pooled, profile_counts,
            "complete pooled replay including multiplicity"
        );
    }
    for pair in &pairs {
        let query = spell(pair, &f.sources);
        let cover = coverage(&f.truth, &query);
        matched |= cover["full_truth_coverage"] == 1.0 && cover["full_query_coverage"] == 1.0;
        let mut class: Vec<_> = query.iter().map(|s| canonical(s)).collect();
        class.sort();
        dna_classes.insert(class);
        let main_bp = if near { 600 } else { 1536 };
        let width = if near { 1 } else { 24 };
        let positions: &[usize] = if near { &[220, 260] } else { &[256, 704, 1216] };
        let donor = if near { 3 } else { 1 };
        let dosage: Vec<_> = positions.iter().map(|&start| json!({"tract":[start,start+width],"truth_modified_copies":f.truth.iter().filter(|s|s.len()==main_bp&&s[start..start+width]==f.sources[donor].1[start..start+width]).count(),"query_modified_copies":query.iter().filter(|s|s.len()==main_bp&&s[start..start+width]==f.sources[donor].1[start..start+width]).count()})).collect();
        let count_equal = sum(&profiles[&pair.copies[0]], &profiles[&pair.copies[1]]) == cis;
        let modified_truth: Vec<_> = f
            .truth
            .iter()
            .filter(|s| s.len() == main_bp && **s != f.sources[0].1)
            .cloned()
            .collect();
        let main_query: Vec<_> = query
            .iter()
            .filter(|s| s.len() == main_bp)
            .cloned()
            .collect();
        let modified = coverage(&modified_truth, &main_query);
        assessments.push(json!({"pair":pair,"coverage":cover,"local_tract_dosage":dosage,"modified_truth_vs_all_main_query_exact_coverage":modified,"modified_copy_full_bp_denominator":main_bp,"tract_bp_denominator_per_copy":width,"reference_integer_counts_equal":count_equal,"non_native_modified_truth":modified_truth.iter().all(|t|f.sources.iter().all(|(_,s)|canonical(t)!=canonical(s)))}));
    }
    let mut best_native = f64::INFINITY;
    audit.ledger_passes += 1;
    for (row, capacity, bytes) in ledger_rows(&ledger_path) {
        audit.ledger_row(&row, capacity, bytes);
        if row["kind"] == "paired_score" && row["native"] == true {
            best_native = best_native.min(row["scored"]["objective"].as_f64().unwrap());
        }
    }
    let best = result["best_paired_objective"].as_f64().unwrap();
    let class = |pair: &Pair| {
        let mut dna: Vec<_> = spell(pair, &f.sources)
            .iter()
            .map(|s| canonical(s))
            .collect();
        dna.sort();
        dna
    };
    let cis_class = class(&Pair {
        copies: [
            domain[cis_indices[0]].clone(),
            domain[cis_indices[1]].clone(),
        ],
    });
    let trans_class = class(&Pair {
        copies: [
            domain[trans_indices[0]].clone(),
            domain[trans_indices[1]].clone(),
        ],
    });
    let mut generated_phase = Vec::new();
    for (label, target) in [("cis", cis_class), ("trans", trans_class)] {
        let mut objectives = Vec::new();
        audit.ledger_passes += 1;
        for (row, capacity, bytes) in ledger_rows(&ledger_path) {
            audit.ledger_row(&row, capacity, bytes);
            if row["kind"] == "paired_score"
                && class(&serde_json::from_value(row["scored"]["pair"].clone()).unwrap()) == target
            {
                objectives.push(row["scored"]["objective"].as_f64().unwrap());
            }
        }
        generated_phase.push(json!({"class":label,"generated":!objectives.is_empty(),"generated_objectives":objectives,"retained":pairs.iter().any(|p|class(p)==target)}));
    }
    let report = json!({"automatic_phase_classes":generated_phase,"case":f.kind,"truth_exact_pair_in_support":matched,"best_native_pair":best_native,"best_found":best,"improvement_over_best_native":best_native-best,"selected_pair_assessments":assessments,"selected_DNA_classes":dna_classes.len(),"cis_trans_exact_operator_count_equality":cis == trans,"cis_trans_differing_count_features":cis.keys().chain(trans.keys()).copied().collect::<BTreeSet<_>>().iter().filter(|t|cis.get(*t)!=trans.get(*t)).count(),"wrong_phase_minus_correct_objective":phase_score_separation,"inference_pairs_independently_audited":audited,"independent_public_haploid_calls":public_losses.len(),"independent_haploid_losses":public_losses,"restricted_oracle_pairs":oracle_rows.len(),"audit_budget":audit,"oracle_minimum":oracle_rows.iter().map(|r|r["objective"].as_f64().unwrap()).fold(f64::INFINITY,f64::min),"replay_reads_total_both_orientations":replay_reads,"replay_DNA_cache_entries":replay_cache.len(),"independent_physical_validations":validations,"actual_reads":sample.stats.reads,"actual_bases":sample.stats.bases,"diploid_bp":f.truth.iter().map(Vec::len).sum::<usize>(),"realized_depth_per_diploid_base":sample.stats.bases as f64/f.truth.iter().map(Vec::len).sum::<usize>() as f64,"nominal_per_copy_depth":10,"oracle_not_inference_input":true,"full_route_support_complete":result["support_complete"],"phase_status":"uncertified full-route support; reference count distinction assessed separately","genotype_uniqueness_certified":false,"complete_observation_vectors_checked":true,"non_objective_pair_vector_comparisons":audited + pairs.len() + if near { 1 } else { 2 },"replay_pair_vector_comparisons":pairs.len() + 2});
    write(&out.join("assessment.json"), &report);
    for v in &frozen {
        assert_eq!(
            v["fingerprint"],
            genome::reconstruction::fingerprint(&inference.join(v["file"].as_str().unwrap()))
                .unwrap()
        );
    }
    report
}

pub fn assess_dosage(f: &Fixture, name: &str, result: &Value) -> Value {
    let inference = f.root.join(name);
    let out = inference.with_file_name(format!("{name}-assessment"));
    fs::create_dir(&out).unwrap();
    let prefix = f.root.join("panel.syng");
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    let panel = SyngIndex::load(p(&prefix), Default::default()).unwrap();
    let graph = routes::Graph::load(&f.root.join("routes"), &identity).unwrap();
    assert_eq!(graph.families.len(), 2);
    let sample = sample::SampleIndex::load(&f.root.join("sample.membwt"), &identity).unwrap();
    let observed = sample.counts.observed_pairs().unwrap();
    let mut e = routes::Evaluator::new(
        &f.root.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        50_000,
        50_000,
    )
    .unwrap();
    let starts = [256usize, 576, 896, 1216];
    let domain = oracle_at(&mut e, &starts);
    assert_eq!(domain.len(), 17);
    let profile = |e: &mut routes::Evaluator<'_>, a: &routes::Assignment| {
        e.evaluate(a)
            .unwrap()
            .factors
            .into_iter()
            .filter_map(|factor| {
                (factor.counts_by_length[0] > 0)
                    .then_some((factor.tokens, factor.counts_by_length[0]))
            })
            .collect::<BTreeMap<_, _>>()
    };
    let domain_profiles: Vec<_> = domain.iter().map(|a| profile(&mut e, a)).collect();
    let truth_counts = sum(&domain_profiles[10], &domain_profiles[12]);
    let mut finite = Vec::new();
    for i in 0..domain.len() {
        for j in i..domain.len() {
            let counts = sum(&domain_profiles[i], &domain_profiles[j]);
            finite.push(json!({"pair":[i,j],"objective":independent_loss(&counts,&observed),"truth_counts":counts==truth_counts}));
        }
    }
    assert_eq!(finite.len(), 153);
    let finite_minimum = finite
        .iter()
        .map(|row| row["objective"].as_f64().unwrap())
        .fold(f64::INFINITY, f64::min);
    let native_rows: Vec<_> = [[0usize, 0usize], [0, 16], [16, 16]]
        .into_iter()
        .map(|indices| {
            let row = finite
                .iter()
                .find(|row| row["pair"] == json!(indices))
                .unwrap();
            json!({"pair":indices,"objective":row["objective"],"count_vector_differs":row["truth_counts"]==false})
        })
        .collect();
    let retained: Vec<Pair> =
        serde_json::from_value(result["retained_correlated_best_found"].clone()).unwrap();
    assert!(!retained.is_empty(), "no retained paired hypothesis");
    let a = &f.sources[0].1;
    let b = &f.sources[1].1;
    let normalize_main = |sequence: &[u8]| {
        let reverse = impg::graph::reverse_complement(sequence);
        let mismatches = |candidate: &[u8]| candidate.iter().zip(a).filter(|(x, y)| x != y).count();
        if mismatches(sequence) <= mismatches(&reverse) {
            sequence.to_vec()
        } else {
            reverse
        }
    };
    let main_mask = |sequence: &[u8]| {
        let sequence = normalize_main(sequence);
        assert_eq!(sequence.len(), 1536);
        let mut expected = a.clone();
        let mut mask = 0usize;
        for (bit, start) in starts.into_iter().enumerate() {
            if sequence[start..start + 24] == b[start..start + 24] {
                expected[start..start + 24].copy_from_slice(&b[start..start + 24]);
                mask |= 1 << bit;
            } else {
                assert_eq!(sequence[start..start + 24], a[start..start + 24]);
            }
        }
        assert_eq!(
            sequence, expected,
            "query contains a non-panel tract allele"
        );
        mask
    };
    let mut selected = Vec::new();
    for pair in &retained {
        for copy in &pair.copies {
            e.validate_assignment(copy).unwrap();
        }
        let counts = sum(
            &profile(&mut e, &pair.copies[0]),
            &profile(&mut e, &pair.copies[1]),
        );
        let query = spell(pair, &f.sources);
        let mains: Vec<_> = query
            .iter()
            .filter(|sequence| sequence.len() == 1536)
            .collect();
        let passives: Vec<_> = query
            .iter()
            .filter(|sequence| sequence.len() == 400)
            .collect();
        assert_eq!((mains.len(), passives.len(), query.len()), (2, 2, 4));
        let masks = [main_mask(mains[0]), main_mask(mains[1])];
        let dosage: Vec<_> = (0..4)
            .map(|bit| {
                usize::from(masks[0] & (1 << bit) != 0) + usize::from(masks[1] & (1 << bit) != 0)
            })
            .collect();
        let both_mosaic = masks.iter().all(|&mask| mask != 0 && mask != 15);
        let mut selected_phase: Vec<_> = mains.iter().map(|sequence| canonical(sequence)).collect();
        selected_phase.sort();
        let mut truth_phase: Vec<_> = f
            .truth
            .iter()
            .filter(|sequence| sequence.len() == 1536)
            .map(|sequence| canonical(sequence))
            .collect();
        truth_phase.sort();
        selected.push(json!({
            "pair":pair,
            "objective":independent_loss(&counts,&observed),
            "exact_truth_count_vector":counts==truth_counts,
            "dosage":dosage,
            "both_primaries_non_native":both_mosaic,
            "exact_whole_phase":selected_phase==truth_phase,
            "coverage":coverage(&f.truth,&query),
            "truth_bp":f.truth.iter().map(Vec::len).sum::<usize>(),
            "query_bp":query.iter().map(Vec::len).sum::<usize>(),
            "truth_primary_bp":3072,
            "query_primary_bp":mains.iter().map(|sequence|sequence.len()).sum::<usize>(),
            "truth_tract_bp":192,
            "query_tract_bp":mains.len()*starts.len()*24,
        }));
    }
    let best = result["best_paired_objective"].as_f64().unwrap();
    let report = json!({
        "case":"dosage-mosaic",
        "retained":selected,
        "best_found":best,
        "finite_153_minimum":finite_minimum,
        "finite_pairs":finite.len(),
        "truth_dosage":[0,1,1,2],
        "native_pairs":native_rows,
        "phase_required":false,
        "phase_reason":"the two heterozygous tracts are separated beyond L150 linkage",
        "support_complete":result["support_complete"],
        "global_bound":result["global_bound"],
        "inference_received_truth_or_domain":false,
        "actual_reads":sample.stats.reads,
        "actual_bases":sample.stats.bases,
    });
    write(&out.join("assessment.json"), &report);
    write(&out.join("finite-153.json"), &finite);
    report
}

fn json_charge(value: &Value) -> usize {
    128 + match value {
        Value::String(s) => 2 * s.len(),
        Value::Array(a) => a.iter().map(json_charge).sum(),
        Value::Object(o) => o
            .iter()
            .map(|(k, v)| 128 + 2 * k.len() + json_charge(v))
            .sum(),
        _ => 0,
    }
}

/// Independent assessment units; finite reference formulas are separate from inference audit.
#[derive(Default, serde::Serialize)]
struct AuditBudget {
    metadata_bytes: usize,
    ledger_passes: usize,
    ledger_rows_read: usize,
    ledger_bytes_read: usize,
    peak_row_and_buffer_bytes: usize,
    public_haploid_calls: usize,
    inference_formula_evaluations: usize,
    reference_formula_evaluations: usize,
    retained_count_cells: usize,
    conservative_peak_bytes: usize,
}
impl AuditBudget {
    fn ledger_row(&mut self, row: &Value, capacity: usize, bytes: usize) {
        self.ledger_rows_read += 1;
        self.ledger_bytes_read += bytes;
        self.peak_row_and_buffer_bytes = self
            .peak_row_and_buffer_bytes
            .max(json_charge(row) + capacity + 8192);
        self.reserve(0);
    }
    fn reserve(&mut self, cells: usize) {
        let reserved = self.retained_count_cells + cells;
        // Full bounded temporary vectors plus conservative recursive JSON/metadata charge.
        let bytes =
            (reserved + 3 * 50_000) * 256 + self.metadata_bytes + self.peak_row_and_buffer_bytes;
        assert!(
            reserved <= 2_000_000 && bytes <= 128 * 1024 * 1024,
            "independent audit state cap"
        );
        self.conservative_peak_bytes = self.conservative_peak_bytes.max(bytes);
    }
    fn profile(&mut self) {
        assert!(self.public_haploid_calls < 4096);
        self.reserve(50_000);
        self.public_haploid_calls += 1;
    }
    fn retain(&mut self, q: &BTreeMap<[u64; 3], u64>) {
        self.reserve(q.len());
        self.retained_count_cells += q.len();
    }
    fn formulas(&mut self, reference: bool) {
        if reference {
            assert!(self.reference_formula_evaluations + 2 <= 90);
            self.reference_formula_evaluations += 2;
        } else {
            assert!(self.inference_formula_evaluations + 2 <= 8192);
            self.inference_formula_evaluations += 2;
        }
    }
}

pub fn assert_observation_recovery(report: &Value, result: &Value, near: bool) {
    assert_eq!(report["complete_observation_vectors_checked"], true);
    assert_eq!(result["support_complete"], false);
    assert!(result["global_bound"].is_null());
    assert_eq!(result["full_gate_b_complete"], false);
    assert_eq!(result["sequence_emission_authorized"], false);
    assert_eq!(report["genotype_uniqueness_certified"], false);
    assert!(
        (report["best_found"].as_f64().unwrap() - report["oracle_minimum"].as_f64().unwrap()).abs()
            <= 1e-9
    );
    let selected = report["selected_pair_assessments"].as_array().unwrap();
    assert!(!selected.is_empty());
    for pair in selected {
        assert_eq!(pair["pair"]["copies"].as_array().unwrap().len(), 2);
        assert_eq!(pair["reference_integer_counts_equal"], true);
        assert!(pair["local_tract_dosage"]
            .as_array()
            .unwrap()
            .iter()
            .all(|d| d["truth_modified_copies"] == d["query_modified_copies"]));
        if near {
            assert_eq!(pair["coverage"]["full_truth_coverage"], 1.0);
            assert_eq!(pair["coverage"]["full_query_coverage"], 1.0);
        }
    }
    assert_eq!(report["cis_trans_exact_operator_count_equality"], !near);
    if near {
        assert_eq!(report["truth_exact_pair_in_support"], true);
        assert!(
            report["wrong_phase_minus_correct_objective"]
                .as_f64()
                .unwrap()
                > 1e-9
        );
        assert!(
            report["cis_trans_differing_count_features"]
                .as_u64()
                .unwrap()
                > 0
        );
    }
}

pub fn coupled_fixture(root: &Path) -> Fixture {
    fs::create_dir(root).unwrap();
    write(
        &root.join("recipe-before-operator.json"),
        &json!({"seed":113,"bp":220,"native_molecules":2,"whole_copies":2,"same_family":true,"read_length":150,"stride":15,"append_last":true,"alternating_global_RC":true}),
    );
    let a = dna(220, 113);
    let sources = vec![("A#0#one".into(), a.clone()), ("A#0#two".into(), a.clone())];
    let truth = vec![a.clone(), a.clone(), a.clone(), a];
    let mut panel = SyngIndex::build(SyncmerParams::default(), sources.clone().into_iter());
    let prefix = root.join("panel.syng");
    panel.save(p(&prefix)).unwrap();
    fasta(&root.join("sources.fa"), &sources);
    write(
        &root.join("groups.json"),
        &json!({"version":1,"groups":sources.iter().enumerate().map(|(i,(n,s))|json!({"id":format!("g{i}"),"scaffold":null,"occurrences":[{"path":n,"start":0,"end":s.len(),"strand":null}]})).collect::<Vec<_>>()}),
    );
    cli(
        root,
        "catalog",
        &[
            "genome-infer",
            "build-catalog",
            "--panel",
            p(&prefix),
            "--groups",
            p(&root.join("groups.json")),
            "--out-dir",
            p(&root.join("catalog")),
        ],
    );
    cli(
        root,
        "routes",
        &[
            "genome-infer",
            "build-panel-routes",
            "--panel",
            p(&prefix),
            "--catalog",
            p(&root.join("catalog/catalog.json")),
            "--sources",
            p(&root.join("sources.fa")),
            "--read-lengths",
            "150",
            "--core-bp",
            "37",
            "--out-dir",
            p(&root.join("routes")),
        ],
    );
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    make_sample(root, "sample", &truth, &panel, &identity);
    let fixture = Fixture {
        root: root.into(),
        sources,
        truth,
        kind: "coupled".into(),
    };
    write(&root.join("fixture-description.json"), &fixture);
    fixture
}
