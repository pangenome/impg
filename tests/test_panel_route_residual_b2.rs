//! Frozen B2 controls and narrow service/accounting boundary regressions.
#[allow(dead_code)]
#[path = "../examples/panel_route_residual_search/mod.rs"]
mod search;
use clap::Parser;
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::{SyncmerParams, SyngIndex},
};
use serde_json::json;
use std::{
    collections::BTreeSet,
    fs,
    io::Read,
    path::{Path, PathBuf},
    process::Command,
};
const STRUCTURAL_CAP: usize = 100_000;
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
fn text(path: &Path) -> &str {
    path.to_str().unwrap()
}
fn write(path: &Path, value: &impl serde::Serialize) {
    genome::write_json(path, value).unwrap();
}
fn fasta(path: &Path, sequences: &[(String, Vec<u8>)]) {
    fs::write(
        path,
        sequences
            .iter()
            .map(|(name, sequence)| format!(">{name}\n{}\n", String::from_utf8_lossy(sequence)))
            .collect::<String>(),
    )
    .unwrap();
}
fn command(root: &Path, name: &str, args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(args)
        .output()
        .unwrap();
    fs::write(root.join(format!("{name}.stdout")), &output.stdout).unwrap();
    fs::write(root.join(format!("{name}.stderr")), &output.stderr).unwrap();
    assert!(output.status.success(), "{name}: {:?}", output.status);
}
fn graph(root: &Path, sequences: &[(String, Vec<u8>)]) -> (SyngIndex, routes::Graph) {
    fs::create_dir(root).unwrap();
    let prefix = root.join("panel.syng");
    let mut panel = SyngIndex::build(SyncmerParams::default(), sequences.to_vec().into_iter());
    panel.save(text(&prefix)).unwrap();
    fasta(&root.join("sources.fa"), sequences);
    write(
        &root.join("groups.json"),
        &json!({"version":1,"groups":sequences.iter().enumerate().map(|(i,(name,s))|json!({"id":format!("g{i}"),"scaffold":null,"occurrences":[{"path":name,"start":0,"end":s.len(),"strand":null}]})).collect::<Vec<_>>()}),
    );
    command(
        root,
        "build-catalog",
        &[
            "genome-infer",
            "build-catalog",
            "--panel",
            text(&prefix),
            "--groups",
            text(&root.join("groups.json")),
            "--out-dir",
            text(&root.join("catalog")),
        ],
    );
    command(
        root,
        "build-routes",
        &[
            "genome-infer",
            "build-panel-routes",
            "--panel",
            text(&prefix),
            "--catalog",
            text(&root.join("catalog/catalog.json")),
            "--sources",
            text(&root.join("sources.fa")),
            "--read-lengths",
            "150",
            "--core-bp",
            "37",
            "--out-dir",
            text(&root.join("routes")),
        ],
    );
    let identity = genome::PanelIdentity::read(text(&prefix)).unwrap();
    let graph = routes::Graph::load(&root.join("routes"), &identity).unwrap();
    (panel, graph)
}
#[derive(Clone, Debug, serde::Serialize)]
struct Port {
    source: usize,
    cut: u64,
    reverse: bool,
    word: Vec<u8>,
}
// Independent decoding, not the experimental lookup helper being tested.
fn ports(root: &Path, graph: &routes::Graph) -> Vec<Port> {
    let mut file = fs::File::open(root.join("routes").join(&graph.ports.path)).unwrap();
    (0..graph.port_count)
        .map(|_| {
            let mut word = vec![0; graph.k as usize];
            let (mut source, mut anchor, mut reverse) = ([0; 8], [0; 8], [0; 1]);
            file.read_exact(&mut word).unwrap();
            file.read_exact(&mut source).unwrap();
            file.read_exact(&mut anchor).unwrap();
            file.read_exact(&mut reverse).unwrap();
            assert!(reverse[0] <= 1);
            let reverse = reverse[0] == 1;
            Port {
                source: u64::from_le_bytes(source) as usize,
                cut: u64::from_le_bytes(anchor)
                    + if reverse {
                        graph.k - graph.cut_offset
                    } else {
                        graph.cut_offset
                    },
                reverse,
                word,
            }
        })
        .collect()
}
fn capacity(route: &routes::Route) -> bool {
    route.segments.iter().enumerate().all(|(i, a)| {
        route.segments[i + 1..]
            .iter()
            .all(|b| a.source != b.source || a.end <= b.start || b.end <= a.start)
    })
}
fn spell(route: &routes::Route, sequences: &[(String, Vec<u8>)]) -> Vec<u8> {
    route
        .segments
        .iter()
        .flat_map(|segment| {
            let s =
                sequences[segment.source].1[segment.start as usize..segment.end as usize].to_vec();
            if segment.reverse {
                impg::graph::reverse_complement(&s)
            } else {
                s
            }
        })
        .collect()
}
struct Oracle<'a> {
    graph: &'a routes::Graph,
    ports: &'a [Port],
    target: usize,
    states: usize,
    candidates: usize,
    routes: BTreeSet<routes::Route>,
}
impl Oracle<'_> {
    fn visit(&mut self, source: usize, entry: u64, reverse: bool, pieces: Vec<routes::Segment>) {
        self.states += 1;
        assert!(
            self.states < STRUCTURAL_CAP,
            "incomplete: structural state cap"
        );
        if source == self.target && !reverse && entry < self.graph.lanes[source].length {
            self.candidates += 1;
            assert!(
                self.candidates < STRUCTURAL_CAP,
                "incomplete: structural candidate cap"
            );
            let mut completed = pieces.clone();
            completed.push(routes::Segment {
                source,
                start: entry,
                end: self.graph.lanes[source].length,
                reverse,
            });
            self.routes.insert(
                routes::Route {
                    segments: completed,
                }
                .normalize(),
            );
        }
        if pieces.len() == 3 {
            return;
        }
        for exit in self.ports {
            if exit.source != source
                || exit.reverse != reverse
                || if reverse {
                    exit.cut >= entry
                } else {
                    exit.cut <= entry
                }
            {
                continue;
            }
            for next in self.ports {
                if next.word != exit.word
                    || (next.source == source && next.cut == exit.cut && next.reverse == reverse)
                {
                    continue;
                }
                let mut child = pieces.clone();
                child.push(routes::Segment {
                    source,
                    start: entry.min(exit.cut),
                    end: entry.max(exit.cut),
                    reverse,
                });
                self.visit(next.source, next.cut, next.reverse, child);
            }
        }
    }
}
#[test]
fn frozen_b2_no_return_structural_and_automatic() {
    no_return_gate(false);
}
#[test]
#[ignore = "requires explicit IMPG_TEST_RESIDUAL_B2_BIN and fresh IMPG_TEST_RESIDUAL_B2_CLI_OUTPUT"]
fn configured_b2_no_return_public_cli() {
    no_return_gate(true);
}
fn no_return_gate(cli: bool) {
    let temporary = tempfile::tempdir().unwrap();
    let retained = if cli {
        Some(PathBuf::from(
            std::env::var_os("IMPG_TEST_RESIDUAL_B2_CLI_OUTPUT")
                .expect("fresh explicit B2 CLI output required"),
        ))
    } else {
        std::env::var_os("IMPG_TEST_RESIDUAL_B2_OUTPUT").map(PathBuf::from)
    };
    let root = retained.as_deref().unwrap_or(temporary.path());
    if cli {
        fs::create_dir(root).unwrap();
    }
    let probe = root.join("word-probe");
    let (_, probe_graph) = graph(&probe, &[("Probe#0#words".into(), dna(400, 113))]);
    assert_eq!(probe_graph.k, 63);
    let words: Vec<_> = ports(&probe, &probe_graph)
        .into_iter()
        .filter(|p| !p.reverse)
        .map(|p| p.word)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .take(3)
        .collect();
    assert_eq!(words.len(), 3);
    let sequences: Vec<_> = [(0, 2), (0, 1), (1, 2)]
        .into_iter()
        .enumerate()
        .map(|(i, (left, right))| {
            let mut sequence = dna(400, 1009 + 77 * i as u64);
            sequence[90..153].copy_from_slice(&words[left]);
            sequence[250..313].copy_from_slice(&words[right]);
            (format!("{}#0#chain", ["A", "B", "C"][i]), sequence)
        })
        .collect();
    let chain = root.join("chain");
    let (panel, graph) = graph(&chain, &sequences);
    let all_ports = ports(&chain, &graph);
    write(&chain.join("public-ports.json"), &all_ports);
    let truth: Vec<_> = [
        &sequences[0].1[..121],
        &sequences[1].1[121..281],
        &sequences[2].1[121..281],
        &sequences[0].1[281..],
    ]
    .concat();
    let last = truth.len() - 150;
    let mut starts: Vec<_> = (0..=last).step_by(15).collect();
    if starts.last() != Some(&last) {
        starts.push(last);
    }
    fasta(
        &chain.join("reads.fa"),
        &starts
            .into_iter()
            .map(|start| (format!("read-{start}"), truth[start..start + 150].to_vec()))
            .collect::<Vec<_>>(),
    );
    let sample = sample::build(&panel, graph.panel.clone(), &[chain.join("reads.fa")]).unwrap();
    sample.save(&chain.join("sample.membwt")).unwrap();
    let mut evaluator = routes::Evaluator::new(
        &chain.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        50_000,
        50_000,
    )
    .unwrap();
    let joins: Vec<_> = [(0,121,1,121), (1,281,2,121), (2,281,0,281)].into_iter().map(|(a,x,b,y)| {
        let left = all_ports.iter().find(|p| p.source == a && p.cut == x && !p.reverse);
        let right = all_ports.iter().find(|p| p.source == b && p.cut == y && !p.reverse);
        json!({"left":left,"right":right,"supported":left.zip(right).is_some_and(|(l,r)| l.word == r.word)})
    }).collect();
    write(
        &chain.join("frozen-recipe.json"),
        &json!({"seeds":[1009,1086,1163],"probe_seed":113,"words":words,"anchors":[90,250],"truth_bp":truth.len(),"joins":joins,"optimization_run":false}),
    );
    assert!(
        joins.iter().all(|j| j["supported"] == true),
        "frozen joins absent: no recipe retry authorized"
    );
    let options = search::Options::parse_from([
        "search",
        "--panel",
        text(&chain.join("panel.syng")),
        "--routes",
        text(&chain.join("routes")),
        "--sample",
        text(&chain.join("sample.membwt")),
        "--out-dir",
        text(&chain.join("not-run")),
    ]);
    let base = graph.native_assignment(0).unwrap();
    let mut meter = search::Meter::new(&options);
    let mut reachable = Vec::new();
    for level in 0..=options.max_level {
        let n = 1u64 << level;
        for cell in 0..2 * n - 1 {
            let l = search::geometry::nearest(
                &mut evaluator,
                &base.routes[0],
                400 * cell / (2 * n),
                &mut meter,
            )
            .unwrap();
            let r = search::geometry::nearest(
                &mut evaluator,
                &base.routes[0],
                400 * (cell + 2) / (2 * n),
                &mut meter,
            )
            .unwrap();
            if l.as_ref().is_some_and(|(p, _)| *p == 121)
                && r.as_ref().is_some_and(|(p, _)| *p == 281)
            {
                reachable.push(json!({"level":level,"cell":cell,"left":l,"right":r}));
            }
        }
    }
    write(
        &chain.join("snapped-coverage.json"),
        &json!({"required":[121,281],"reachable":reachable,"meter":meter,"truth_coordinates_not_in_inference":true}),
    );
    // The original failing preflight is preserved outside the source tree.
    // Authorized follow-up: structural closure/oracle only, never optimization.
    let mut summaries = Vec::new();
    let mut all_rows = Vec::new();
    let mut target_count = 0;
    let mut one_donor_target_count = 0;
    for family in 0..graph.families.len() {
        let target = graph.families[family].paths[0];
        let mut oracle = Oracle {
            graph: &graph,
            ports: &all_ports,
            target,
            states: 0,
            candidates: 0,
            routes: BTreeSet::new(),
        };
        oracle.visit(target, 0, false, vec![]);
        let (mut feasible, mut rejected) = (0, 0);
        for route in &oracle.routes {
            let mut assignment = graph.native_assignment(family).unwrap();
            assignment.routes[0] = route.clone();
            let independent = capacity(route);
            let public = evaluator.validate_assignment(&assignment);
            assert_eq!(
                independent,
                public.is_ok(),
                "public geometry/capacity mismatch: {assignment:?} {public:?}"
            );
            if independent {
                feasible += 1;
            } else {
                rejected += 1;
            }
            let exact = independent && spell(route, &sequences) == truth;
            if exact {
                target_count += 1;
                if route.segments.len() <= 3 {
                    one_donor_target_count += 1;
                }
            }
            all_rows.push(json!({"assignment":assignment,"capacity":independent,"public_feasible":public.is_ok(),"error":public.err().map(|e|e.to_string()),"target_dna_equal":exact}));
        }
        summaries.push(json!({"family":family,"states":oracle.states,"complete_candidates":oracle.candidates,"normalized_routes":oracle.routes.len(),"feasible":feasible,"rejected":rejected,"exhausted_declared_domain":true}));
    }
    write(&chain.join("structural-oracle.json"), &all_rows);
    write(
        &chain.join("structural-summary.json"),
        &json!({"domain":"all <=3 noncontiguous switches, positive oriented public ports, native endpoints, all families; independent final canonical capacity","families":summaries,"target_count":target_count,"native_baseline_one_donor_target_count":one_donor_target_count,"optimization_run":false,"public_objectives_computed":0}),
    );
    assert!(target_count > 0);
    assert_eq!(one_donor_target_count, 0);
    refinement_closure(&chain, &mut evaluator, &base, &options, &all_rows);
    automatic_chain(&chain, &sequences, &truth, cli);
    // Already declared finite-profile audit, only AFTER automatic output freeze.
    assert!(all_rows.len() <= 1024);
    let profiles: Vec<_> = all_rows
        .iter()
        .filter(|r| r["public_feasible"] == true)
        .map(|r| {
            let a: routes::Assignment = serde_json::from_value(r["assignment"].clone()).unwrap();
            evaluator.evaluate(&a).unwrap()
        })
        .collect();
    write(
        &chain.join("post-inference-finite-profiles.json"),
        &profiles,
    );
}

fn refinement_closure(
    root: &Path,
    evaluator: &mut routes::Evaluator<'_>,
    base: &routes::Assignment,
    options: &search::Options,
    oracle: &[serde_json::Value],
) {
    use search::geometry::Region;
    let mut meter = search::Meter::new(options);
    // Fixed fields only: region, parent index, BTree key and queue index. This
    // provision is logical capacity accounting, not measured allocator/RSS use.
    meter.reserve(STRUCTURAL_CAP as u64 * 512).unwrap();
    let mut nodes: Vec<(Region, Option<usize>)> = Vec::with_capacity(STRUCTURAL_CAP);
    let mut keys = BTreeSet::new();
    for level in 0..=options.max_level {
        let n = 1u64 << level;
        for cell in 0..2 * n - 1 {
            let region = Region {
                slot: 0,
                left: 400 * cell / (2 * n),
                right: 400 * (cell + 2) / (2 * n),
                level,
            };
            if keys.insert((level, region.left, region.right)) {
                meter.work(1).unwrap();
                nodes.push((region, None));
            }
        }
    }
    let mut cursor = 0;
    let mut found = None;
    let mut stop = "exhausted_permitted_refinement_closure".to_string();
    while cursor < nodes.len() {
        let region = nodes[cursor].0.clone();
        let attempt = (|| -> std::io::Result<_> {
            meter.work(1)?;
            let left =
                search::geometry::nearest(evaluator, &base.routes[0], region.left, &mut meter)?;
            let right =
                search::geometry::nearest(evaluator, &base.routes[0], region.right, &mut meter)?;
            Ok((left, right))
        })();
        let (left, right) = match attempt {
            Ok(pair) => pair,
            Err(error) => {
                stop = error.to_string();
                break;
            }
        };
        if left.as_ref().is_some_and(|(p, _)| *p == 121)
            && right.as_ref().is_some_and(|(p, _)| *p == 281)
        {
            found = Some(cursor);
            stop = "conditional_witness_found_not_exhausted".into();
            break;
        }
        for child in region.refinements(400, options.max_level) {
            let key = (child.level, child.left, child.right);
            if !keys.contains(&key) {
                if nodes.len() == STRUCTURAL_CAP {
                    stop = "budget: structural_states".into();
                    break;
                }
                if let Err(error) = meter.work(1) {
                    stop = error.to_string();
                    break;
                }
                keys.insert(key);
                nodes.push((child, Some(cursor)));
            }
        }
        if stop.starts_with("budget:") {
            break;
        }
        cursor += 1;
    }
    let mut ancestry = Vec::new();
    let mut next = found;
    while let Some(index) = next {
        ancestry.push(index);
        next = nodes[index].1;
    }
    ancestry.reverse();
    let mut path = Vec::new();
    for index in ancestry {
        let (region, parent) = &nodes[index];
        let left =
            search::geometry::nearest(evaluator, &base.routes[0], region.left, &mut meter).unwrap();
        let right = search::geometry::nearest(evaluator, &base.routes[0], region.right, &mut meter)
            .unwrap();
        let candidates = oracle
            .iter()
            .filter(|row| {
                if row["public_feasible"] != true || row["assignment"]["family"] != base.family {
                    return false;
                }
                let a: routes::Assignment =
                    serde_json::from_value(row["assignment"].clone()).unwrap();
                let pieces = &a.routes[0].segments;
                pieces.len() >= 3
                    && left.as_ref().is_some_and(|(p, _)| pieces[0].end == *p)
                    && right
                        .as_ref()
                        .is_some_and(|(p, _)| pieces.last().unwrap().start == *p)
            })
            .count();
        path.push(json!({"index":index,"parent":parent,"region":region,"left":left,"right":right,"feasible_non_native_single_or_two_donor_candidates_at_pair":candidates}));
    }
    write(
        &root.join("refinement-closure.json"),
        &json!({
            "status":stop,"states_retained":nodes.len(),"states_completed":cursor,"active_state":nodes.get(cursor),
            "meter":meter,"path":path,"required_refinement_edges":path.len().saturating_sub(1),
            "conditional_reachability_only":true,"eight_confirmation_and_score_selection_prerequisites_tested":false,
            "optimization_run":false,"objectives_computed":0,
            "scope":"all coarse roots, BFS unique (level,left,right), permitted refinements only; stop on first witness or cap"
        }),
    );
}

fn sample_native(
    root: &Path,
    panel: &SyngIndex,
    graph: &routes::Graph,
    sequences: &[(String, Vec<u8>)],
) -> sample::SampleIndex {
    let mut reads = Vec::new();
    for (source, (_, dna)) in sequences.iter().enumerate() {
        let last = dna.len() - 150;
        let mut starts: Vec<_> = (0..=last).step_by(15).collect();
        if starts.last() != Some(&last) {
            starts.push(last);
        }
        for start in starts {
            reads.push((
                format!("{source}-{start}"),
                dna[start..start + 150].to_vec(),
            ));
        }
    }
    fasta(&root.join("reads.fa"), &reads);
    let sample = sample::build(panel, graph.panel.clone(), &[root.join("reads.fa")]).unwrap();
    sample.save(&root.join("sample.membwt")).unwrap();
    sample
}
fn test_options(root: &Path) -> search::Options {
    search::Options::parse_from([
        "search",
        "--panel",
        text(&root.join("panel.syng")),
        "--routes",
        text(&root.join("routes")),
        "--sample",
        text(&root.join("sample.membwt")),
        "--out-dir",
        text(&root.join("not-run")),
    ])
}
fn positioned_ports(route: &routes::Route, ports: &[Port]) -> Vec<(u64, Port)> {
    let mut rows = Vec::new();
    let mut offset = 0;
    let length = route.length().unwrap();
    for segment in &route.segments {
        for port in ports {
            if port.source != segment.source
                || port.reverse != segment.reverse
                || port.cut < segment.start
                || port.cut > segment.end
            {
                continue;
            }
            let position = offset
                + if segment.reverse {
                    segment.end - port.cut
                } else {
                    port.cut - segment.start
                };
            if position > 0 && position < length {
                rows.push((position, port.clone()));
            }
        }
        offset += segment.end - segment.start;
    }
    rows
}
fn exchange_domain(
    base: &routes::Assignment,
    ports: &[Port],
    opposite: bool,
) -> Vec<(usize, u64, u64, routes::Segment, usize)> {
    let mut rows = Vec::new();
    for (slot, route) in base.routes.iter().enumerate() {
        let positions = positioned_ports(route, ports);
        for (lo, left) in &positions {
            for (hi, right) in &positions {
                if lo >= hi {
                    continue;
                }
                for a in ports.iter().filter(|p| p.word == left.word) {
                    for b in ports.iter().filter(|p| p.word == right.word) {
                        if a.source != b.source
                            || a.reverse != b.reverse
                            || if a.reverse {
                                a.cut <= b.cut
                            } else {
                                a.cut >= b.cut
                            }
                        {
                            continue;
                        }
                        let donor = routes::Segment {
                            source: a.source,
                            start: a.cut.min(b.cut),
                            end: a.cut.max(b.cut),
                            reverse: a.reverse,
                        };
                        for (other, occupied) in base.routes.iter().enumerate() {
                            if other == slot {
                                continue;
                            }
                            if occupied.segments.iter().any(|p| {
                                p.source == donor.source
                                    && p.start <= donor.start
                                    && donor.end <= p.end
                                    && (p.reverse != donor.reverse) == opposite
                            }) {
                                assert!(rows.len() < STRUCTURAL_CAP, "incomplete exchange domain");
                                rows.push((slot, *lo, *hi, donor.clone(), other));
                            }
                        }
                    }
                }
            }
        }
    }
    rows
}
fn assignment_capacity(a: &routes::Assignment) -> bool {
    capacity(&routes::Route {
        segments: a.routes.iter().flat_map(|r| r.segments.clone()).collect(),
    })
}
fn oriented_oracle(
    root: &Path,
    e: &mut routes::Evaluator<'_>,
    base: &routes::Assignment,
    ports: &[Port],
    name: &str,
) -> serde_json::Value {
    let domain = exchange_domain(base, ports, true);
    let mut meter = search::Meter::new(&test_options(root));
    let (mut feasible, mut rejected, mut transported_feasible, mut transported_rejected) =
        (0, 0, 0, 0);
    let mut rows = Vec::new();
    for (slot, lo, hi, donor, other) in domain {
        let removed = search::geometry::crop(&base.routes[slot], lo, hi);
        let candidate = search::oriented::reciprocal_proposal(base, slot, lo, hi, &donor, other)
            .unwrap()
            .unwrap();
        let rejection = search::oriented::seam_rejection(e, &candidate, &mut meter).unwrap();
        let capacity = assignment_capacity(&candidate);
        let public = e.validate_assignment(&candidate);
        assert_eq!(
            rejection.is_none() && capacity,
            public.is_ok(),
            "{candidate:?} {rejection:?} {public:?}"
        );
        if public.is_ok() {
            feasible += 1;
            if removed.len() > 1 {
                transported_feasible += 1;
            }
        } else {
            rejected += 1;
            if removed.len() > 1 {
                transported_rejected += 1;
            }
        }
        rows.push(json!({"slot":slot,"lo":lo,"hi":hi,"donor":donor,"other":other,"removed":removed,"assignment":candidate,"seam_rejection":rejection,"independent_capacity":capacity,"public_feasible":public.is_ok(),"public_error":public.err().map(|e|e.to_string())}));
    }
    let summary = json!({"domain":"all matching actual recipient-port pairs and positive donor-port pairs wholly occupied in another single segment; opposite traversal only; every recipient slot","exhausted_declared_domain":true,"feasible":feasible,"rejected":rejected,"transported_interior_seam_feasible":transported_feasible,"transported_interior_seam_rejected":transported_rejected,"meter":meter,"profiles_computed":0,"objectives_computed":0});
    write(&root.join(format!("{name}-cases.json")), &rows);
    write(&root.join(format!("{name}-summary.json")), &summary);
    summary
}
#[test]
fn frozen_b2_oriented_exact_cut_structural_oracle() {
    let temporary = tempfile::tempdir().unwrap();
    let retained = std::env::var_os("IMPG_TEST_RESIDUAL_B2_OUTPUT").map(PathBuf::from);
    let root = retained.as_deref().unwrap_or(temporary.path());
    let a = dna(220, 113);
    for copies in [2, 3] {
        let root = root.join(format!("opposite-{copies}"));
        let mut sequences = vec![
            ("A#0#one".into(), a.clone()),
            ("A#0#two".into(), impg::graph::reverse_complement(&a)),
        ];
        if copies == 3 {
            sequences.push(("A#0#three".into(), a.clone()));
        }
        let (panel, graph) = graph(&root, &sequences);
        let ports = ports(&root, &graph);
        write(&root.join("public-ports.json"), &ports);
        let sample = sample_native(&root, &panel, &graph, &sequences);
        let mut e = routes::Evaluator::new(
            &root.join("routes"),
            &graph,
            &panel,
            &sample,
            10.0,
            0.1,
            50_000,
            50_000,
        )
        .unwrap();
        let native = graph.native_assignment(0).unwrap();
        assert_eq!(graph.k, 63);
        if copies == 2 {
            let summary = oriented_oracle(&root, &mut e, &native, &ports, "native");
            assert!(summary["feasible"].as_u64().unwrap() > 0);
            assert!(summary["rejected"].as_u64().unwrap() > 0);
            let mut options = test_options(&root);
            options.geometry_b2 = true;
            options.out_dir = root.join("automatic");
            search::run(options.clone(), None).unwrap();
            freeze_inference(&options.out_dir);
            let result: serde_json::Value =
                genome::read_json(&options.out_dir.join("result.json")).unwrap();
            let ledger = ledger(&options.out_dir);
            let native = &ledger.iter().find(|r| r["kind"] == "native").unwrap()["public"];
            let compounds: Vec<_> = ledger
                .iter()
                .filter(|r| r["kind"] == "candidate" && r["compound"] == true)
                .collect();
            assert!(!compounds.is_empty());
            let counts = |p: &serde_json::Value| {
                p["factors"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|f| (f["tokens"].clone(), f["counts_by_length"].clone()))
                    .collect::<Vec<_>>()
            };
            for row in &compounds {
                assert_eq!(counts(&row["public"]), counts(native));
                let a: routes::Assignment =
                    serde_json::from_value(row["public"]["assignment"].clone()).unwrap();
                assert!(assignment_capacity(&a));
                assert!(a.routes.iter().flat_map(|r| &r.segments).any(|p| p.reverse));
            }
            assert!(
                result["retained_correlated_best_found"]
                    .as_array()
                    .unwrap()
                    .len()
                    > 1
            );
            write(
                &options.out_dir.join("assessment.json"),
                &json!({"opposite_compounds":compounds.len(),"full_count_equivalence":true,"independent_capacity_passed":true,"identifiable_recovery_claim":false}),
            );
            let cases: Vec<serde_json::Value> =
                genome::read_json(&root.join("native-cases.json")).unwrap();
            let feasible: Vec<_> = cases
                .iter()
                .filter(|r| r["public_feasible"] == true)
                .collect();
            assert!(feasible.len() <= 1024);
            let profiles: Vec<_> = feasible
                .iter()
                .map(|r| {
                    let a: routes::Assignment =
                        serde_json::from_value(r["assignment"].clone()).unwrap();
                    e.evaluate(&a).unwrap()
                })
                .collect();
            write(&root.join("post-inference-finite-profiles.json"), &profiles);
        } else {
            let mut bases = BTreeSet::new();
            for (slot, lo, hi, donor, other) in exchange_domain(&native, &ports, false) {
                let proposal =
                    search::oriented::reciprocal_proposal(&native, slot, lo, hi, &donor, other)
                        .unwrap()
                        .unwrap();
                if proposal != native && e.validate_assignment(&proposal).is_ok() {
                    bases.insert(proposal);
                }
            }
            let base = bases.first().expect("no same-orientation seam baseline");
            write(&root.join("lexicographic-seam-baseline.json"), base);
            let summary = oriented_oracle(&root, &mut e, base, &ports, "interior");
            // Frozen first-baseline negative stays a negative; the separately
            // approved filtered extension below must supply its own witness.
            assert_eq!(summary["transported_interior_seam_feasible"], 0);
            interior_witness_extension(&root, &mut e, &native, &ports);
        }
    }
}

fn charged_exchange_domain(
    base: &routes::Assignment,
    ports: &[Port],
    opposite: bool,
    meter: &mut search::Meter,
    constructed: &mut usize,
) -> std::io::Result<Vec<(usize, u64, u64, routes::Segment, usize)>> {
    // Charge actual scalar membership/comparison loop visits, not backend CPU
    // instructions. Word comparisons and temporary copies get explicit k cost.
    let k = ports.first().map_or(0, |p| p.word.len()) as u64;
    let mut rows = Vec::new();
    for (slot, route) in base.routes.iter().enumerate() {
        meter.work(route.segments.len() as u64 * ports.len() as u64)?;
        let positions = positioned_ports(route, ports);
        meter.reserve(positions.capacity() as u64 * (128 + k) * 2)?;
        for (lo, left) in &positions {
            for (hi, right) in &positions {
                meter.work(1)?;
                if lo >= hi {
                    continue;
                }
                for a in ports {
                    meter.work(1 + k)?;
                    if a.word != left.word {
                        continue;
                    }
                    for b in ports {
                        meter.work(1 + k)?;
                        if b.word != right.word
                            || a.source != b.source
                            || a.reverse != b.reverse
                            || if a.reverse {
                                a.cut <= b.cut
                            } else {
                                a.cut >= b.cut
                            }
                        {
                            continue;
                        }
                        let donor = routes::Segment {
                            source: a.source,
                            start: a.cut.min(b.cut),
                            end: a.cut.max(b.cut),
                            reverse: a.reverse,
                        };
                        for (other, occupied) in base.routes.iter().enumerate() {
                            meter.work(1 + occupied.segments.len() as u64)?;
                            if other == slot {
                                continue;
                            }
                            if occupied.segments.iter().any(|p| {
                                p.source == donor.source
                                    && p.start <= donor.start
                                    && donor.end <= p.end
                                    && (p.reverse != donor.reverse) == opposite
                            }) {
                                if *constructed == STRUCTURAL_CAP {
                                    return Err(std::io::Error::other(
                                        "budget: total_constructed_candidates",
                                    ));
                                }
                                meter.reserve(256)?;
                                *constructed += 1;
                                rows.push((slot, *lo, *hi, donor.clone(), other));
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(rows)
}
fn interior_witness_extension(
    root: &Path,
    e: &mut routes::Evaluator<'_>,
    native: &routes::Assignment,
    ports: &[Port],
) {
    use std::io::Write;
    let mut meter = search::Meter::new(&test_options(root));
    let mut constructed = 0;
    let mut log = fs::File::create(root.join("extension-cases.jsonl")).unwrap();
    let (mut skipped, mut examined, mut witness) = (0, 0, None);
    let result = (|| -> std::io::Result<()> {
        let mut bases = BTreeSet::new();
        for (slot, lo, hi, donor, other) in
            charged_exchange_domain(native, ports, false, &mut meter, &mut constructed)?
        {
            meter.work(64)?;
            let a = search::oriented::reciprocal_proposal(native, slot, lo, hi, &donor, other)?
                .unwrap();
            meter.validation(&a, e.graph.k, e.graph.port_count)?;
            if a != *native && e.validate_assignment(&a).is_ok() {
                meter.reserve(4096)?;
                bases.insert(a);
            }
        }
        for (index, base) in bases.iter().enumerate() {
            let mut counterparts = true;
            for route in &base.routes {
                for pair in route.segments.windows(2) {
                    for (segment, cut) in [
                        (
                            &pair[0],
                            if pair[0].reverse {
                                pair[0].start
                            } else {
                                pair[0].end
                            },
                        ),
                        (
                            &pair[1],
                            if pair[1].reverse {
                                pair[1].end
                            } else {
                                pair[1].start
                            },
                        ),
                    ] {
                        meter.work(ports.len() as u64)?;
                        counterparts &= ports.iter().any(|p| {
                            p.source == segment.source
                                && p.cut == cut
                                && p.reverse != segment.reverse
                        });
                    }
                }
            }
            writeln!(
                log,
                "{}",
                json!({"kind":"baseline_filter","index":index,"assignment":base,"independent_opposite_counterparts_at_all_seams":counterparts})
            )?;
            if !counterparts {
                skipped += 1;
                continue;
            }
            examined += 1;
            for (slot, lo, hi, donor, other) in
                charged_exchange_domain(base, ports, true, &mut meter, &mut constructed)?
            {
                let removed = search::geometry::crop(&base.routes[slot], lo, hi);
                if removed.len() < 2 {
                    continue;
                }
                meter.work(128)?;
                let a = search::oriented::reciprocal_proposal(base, slot, lo, hi, &donor, other)?
                    .unwrap();
                let rejection = search::oriented::seam_rejection(e, &a, &mut meter)?;
                let capacity = assignment_capacity(&a);
                meter.validation(&a, e.graph.k, e.graph.port_count)?;
                let public = e.validate_assignment(&a);
                assert_eq!(rejection.is_none() && capacity, public.is_ok());
                assert_eq!(
                    a.routes.iter().map(|r| r.length().unwrap()).sum::<u64>(),
                    660
                );
                let row = json!({"kind":"opposite_candidate","baseline_index":index,"baseline":base,"slot":slot,"lo":lo,"hi":hi,"donor":donor,"other":other,"removed":removed,"assignment":a,"seam_rejection":rejection,"capacity":capacity,"public_feasible":public.is_ok(),"error":public.as_ref().err().map(|e|e.to_string())});
                writeln!(log, "{row}")?;
                if public.is_ok() {
                    witness = Some(row);
                    return Ok(());
                }
            }
        }
        Ok(())
    })();
    log.flush().unwrap();
    write(
        &root.join("extension-summary.json"),
        &json!({"status":result.as_ref().err().map(|e|e.to_string()).unwrap_or_else(||if witness.is_some(){"positive_witness_not_exhausted".into()}else{"filtered_domain_exhausted_without_witness".into()}),"constructed_total":constructed,"filtered_out_baselines":skipped,"examined_baselines":examined,"meter":meter,"witness":witness,"unfiltered_expanded_domain_exhausted":false,"profiles_computed":0,"objectives_computed":0}),
    );
    // This preserves an INCOMPLETE diagnostic stop, not the failed positive
    // capability gate. The original required-positive executable/log is retained.
    assert_eq!(result.unwrap_err().to_string(), "budget: work");
    assert!(witness.is_none());
}
fn ledger(root: &Path) -> Vec<serde_json::Value> {
    fs::read_to_string(root.join("scores.jsonl"))
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect()
}
fn freeze_inference(root: &Path) {
    let rows:Vec<_>=["scores.jsonl","events.jsonl","result.json","pending.json","accounting.jsonl","options.json","provenance.json","status.json"].iter().map(|name|json!({"file":name,"fingerprint":genome::reconstruction::fingerprint(&root.join(name)).unwrap()})).collect();
    write(&root.join("frozen-before-dna-assessment.json"), &rows);
}
fn automatic_chain(root: &Path, sequences: &[(String, Vec<u8>)], truth: &[u8], cli: bool) {
    let mut options = test_options(root);
    options.geometry_b2 = true;
    options.out_dir = root.join("automatic");
    if cli {
        let bin = std::env::var_os("IMPG_TEST_RESIDUAL_B2_BIN")
            .expect("explicit standalone B2 binary required");
        let output = Command::new(bin)
            .args([
                "--geometry-b2",
                "--panel",
                &options.panel,
                "--routes",
                text(&options.routes),
                "--sample",
                text(&options.sample),
                "--out-dir",
                text(&options.out_dir),
            ])
            .output()
            .unwrap();
        fs::write(root.join("automatic.stdout"), output.stdout).unwrap();
        fs::write(root.join("automatic.stderr"), output.stderr).unwrap();
        assert!(output.status.success());
    } else {
        search::run(options.clone(), None).unwrap();
    }
    freeze_inference(&options.out_dir);
    let result: serde_json::Value =
        genome::read_json(&options.out_dir.join("result.json")).unwrap();
    let ledger = ledger(&options.out_dir);
    let best_native = ledger
        .iter()
        .filter(|r| r["kind"] == "native")
        .map(|r| r["public"]["relative_objective"].as_f64().unwrap())
        .fold(f64::INFINITY, f64::min);
    let best = result["best_public_objective"].as_f64().unwrap();
    let ties = result["retained_correlated_best_found"].as_array().unwrap();
    let matches = ties
        .iter()
        .filter(|v| {
            let a: routes::Assignment = serde_json::from_value((*v).clone()).unwrap();
            assert!(assignment_capacity(&a));
            a.routes.len() == 1 && spell(&a.routes[0], sequences) == truth
        })
        .count();
    write(
        &options.out_dir.join("assessment.json"),
        &json!({"best_native":best_native,"best_found":best,"improvement":best_native-best,"retained_best":ties.len(),"exact_whole_dna_best":matches,"truth_oracle_objectives_computed":0,"full_b2_complete":false}),
    );
    assert!(best_native - best > 1.0);
    assert!(!ties.is_empty() && matches == ties.len());
    assert!(result["b2"]["chain_confirmations"].as_u64().unwrap() > 0);
    let recorded = ledger
        .iter()
        .filter(|r| {
            matches!(
                r["kind"].as_str(),
                Some("native" | "public_confirmation" | "unconfirmed_exact_delta")
            )
        })
        .count();
    assert_eq!(
        result["meter"]["exact_scores"].as_u64().unwrap() as usize,
        recorded
    );
    let pending: serde_json::Value =
        genome::read_json(&options.out_dir.join("pending.json")).unwrap();
    assert!(!pending["b2"]["queued_chains"]
        .as_array()
        .unwrap()
        .is_empty());
    assert_eq!(result["b2"]["full_b2_complete"], false);
}

fn service_fixture(root: &Path, families: usize) {
    let sequences: Vec<_> = (0..families)
        .map(|i| (format!("{}#0#one", (b'A' + i as u8) as char), dna(220, 113)))
        .collect();
    let (panel, graph) = graph(root, &sequences);
    sample_native(root, &panel, &graph, &sequences);
}
fn service_root(temporary: &Path, name: &str) -> PathBuf {
    std::env::var_os("IMPG_TEST_RESIDUAL_B2_SERVICE_OUTPUT")
        .map(PathBuf::from)
        .unwrap_or_else(|| temporary.into())
        .join(name)
}
fn json_lines(root: &Path, name: &str) -> Vec<serde_json::Value> {
    fs::read_to_string(root.join(name))
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect()
}
#[test]
fn b2_service_ordinary_geometry_before_native_completion() {
    let temporary = tempfile::tempdir().unwrap();
    let root = service_root(temporary.path(), "native-service");
    service_fixture(&root, 5);
    let mut options = test_options(&root);
    options.geometry_b2 = true;
    options.max_epochs = 0;
    options.max_scores = 3; // Three natives admitted, two inventories unfinished.
    options.out_dir = root.join("stopped");
    let result = search::run(options.clone(), None).unwrap();
    let pending: serde_json::Value =
        genome::read_json(&options.out_dir.join("pending.json")).unwrap();
    let events = json_lines(&options.out_dir, "events.jsonl");
    let first_geometry = events.iter().position(|r| r["event"] == "geometry_attempt");
    let natives_before_geometry = first_geometry.map(|index| {
        events[..index]
            .iter()
            .filter(|r| r["event"] == "baseline")
            .count()
    });
    write(
        &root.join("service-observation.json"),
        &json!({"result":result,"pending":pending,"first_geometry_event_index":first_geometry,"natives_before_first_geometry":natives_before_geometry,"derivation":"five inventories; max_scores=3; Coarse then Region need two ordinary-task turns, independent of chain turns"}),
    );
    assert_eq!(result["status"], "budget: exact_scores");
    assert_eq!(result["native_initializations"], 3);
    assert_eq!(result["native_initialization_complete"], false);
    assert_eq!(pending["active_native_family"], 3);
    assert_eq!(pending["pending_native_families"], json!([3, 4]));
    assert!(pending["active_at_stop"].is_null());
    assert!(!pending["queued"].as_array().unwrap().is_empty());
    assert!(pending["b2"]["active_chain"].is_null());
    assert_eq!(pending["b2"]["queued_chains"].as_array().unwrap().len(), 3);
    assert!(result["b2"]["chain_primitives_completed"].as_u64().unwrap() > 0);
    assert_eq!(result["meter"]["native_scores"], 3);
    assert_eq!(result["meter"]["exact_scores"], 3);
    assert!(
        natives_before_geometry.is_some_and(|n| n < 3),
        "ordinary geometry must receive service before native initialization finishes"
    );
}
#[test]
fn b2_service_post_confirmation_stop_keeps_chain_classification() {
    let temporary = tempfile::tempdir().unwrap();
    let root = service_root(temporary.path(), "chain-confirmation");
    service_fixture(&root, 2);
    let mut options = test_options(&root);
    options.geometry_b2 = true;
    options.max_epochs = 0;
    options.max_level = 0;
    // Two native ties + one ordinary root substitution per native. The next
    // distinct chain tie stops after public confirmation, not at score admission.
    options.max_ties = 4;
    options.out_dir = root.join("stopped");
    let result = search::run(options.clone(), None).unwrap();
    let pending: serde_json::Value =
        genome::read_json(&options.out_dir.join("pending.json")).unwrap();
    let events = json_lines(&options.out_dir, "events.jsonl");
    let ledger = ledger(&options.out_dir);
    let candidates: Vec<_> = ledger.iter().filter(|r| r["kind"] == "candidate").collect();
    let is_chain = |row: &&serde_json::Value| {
        events.iter().any(|event| {
            event["event"] == "chain_complete_geometry"
                && event["baseline"] == row["baseline"]
                && event["region"] == row["region"]
                && event["assignment"] == row["public"]["assignment"]
        })
    };
    let chain_records = candidates.iter().copied().filter(is_chain).count();
    let completed_records = ledger
        .iter()
        .filter(|r| {
            matches!(
                r["kind"].as_str(),
                Some("native" | "unconfirmed_exact_delta" | "public_confirmation")
            )
        })
        .count();
    write(
        &root.join("confirmation-observation.json"),
        &json!({"result":result,"pending":pending,"completed_chain_candidate_records":chain_records,"completed_candidate_records":candidates.len(),"completed_score_records":completed_records,"derivation":"max_level=0; two natives + two ordinary root ties fit max_ties=4; first new chain tie confirms then stops"}),
    );
    assert!(result["status"]
        .as_str()
        .unwrap()
        .starts_with("budget: retained_ties"));
    assert_eq!(result["native_initializations"], 2);
    assert_eq!(
        result["retained_correlated_best_found"]
            .as_array()
            .unwrap()
            .len(),
        4
    );
    assert_eq!(candidates.len() - chain_records, 2);
    assert_eq!(chain_records, 1);
    assert!(!pending["b2"]["active_chain"].is_null());
    assert!(!pending["b2"]["queued_chains"]
        .as_array()
        .unwrap()
        .is_empty());
    assert_eq!(
        result["confirmed_candidates"].as_u64().unwrap() as usize,
        candidates.len()
    );
    assert_eq!(
        result["meter"]["exact_scores"].as_u64().unwrap() as usize,
        completed_records
    );
    assert_eq!(result["meter"]["public_scores"], 2 + candidates.len());
    assert_eq!(
        result["meter"]["score_attempts"],
        result["meter"]["exact_scores"]
    );
    assert_eq!(
        result["b2"]["chain_confirmations"].as_u64().unwrap() as usize,
        chain_records,
        "a completed public/parity chain confirmation must survive downstream retention stop"
    );
}
