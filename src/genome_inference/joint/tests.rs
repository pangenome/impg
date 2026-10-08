use super::*;
use crate::{sample_mem_bwt::WeightedBwt, syng::SyncmerParams};
use serde_json::json;
use std::fs;
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
fn identity() -> PanelIdentity {
    PanelIdentity {
        checksum_algorithm: "fixture".into(),
        sidecars: vec![],
    }
}
fn tokens(node: i32) -> [u64; 3] {
    canonical(&encode_walk(&[(node, 0), (node + 1, 10)]).unwrap())
        .try_into()
        .unwrap()
}
fn walk(id: &str, instance: &str) -> layout::Walk {
    layout::Walk {
        id: id.into(),
        topology: "linear".into(),
        left_endpoint: "asserted-molecule-terminus".into(),
        right_endpoint: "asserted-molecule-terminus".into(),
        pieces: vec![layout::Piece {
            instance: instance.into(),
            start: 0,
            end: 100,
            strand: "+".into(),
        }],
        adjacencies: vec![],
    }
}
fn fixture() -> Compiled {
    let layout = layout::Layout {
        version: VERSION,
        model: LAYOUT_MODEL.into(),
        panel: identity(),
        sources: vec![layout::Source {
            id: 0,
            name: "A#0#chr".into(),
            length: 100,
        }],
        instances: ["copy0", "copy1"]
            .iter()
            .map(|id| layout::Instance {
                id: (*id).into(),
                source: 0,
                start: 0,
                end: 100,
            })
            .collect(),
        slots: (0..2)
            .map(|s| layout::Slot {
                id: format!("slot{s}"),
                alternatives: vec![
                    walk("present", &format!("copy{s}")),
                    walk("other", &format!("copy{s}")),
                ],
            })
            .collect(),
    };
    let p = Profile {
        contributions: vec![Contribution {
            tokens: tokens(1),
            totals: vec![[100, 100]],
        }],
        admitted_starts: vec![1],
        event_runs: 1,
        max_crop_bp: 100,
    };
    let z = Profile {
        contributions: vec![],
        ..p.clone()
    };
    Compiled {
        version: VERSION,
        model: MODEL.into(),
        compiler_identity: compiler_identity(),
        count_policy: COUNT_POLICY.into(),
        layout,
        read_lengths: vec![100],
        definitions: vec![
            Definition {
                tokens: tokens(1),
                original_ids: vec![7],
            },
            Definition {
                tokens: tokens(4),
                original_ids: vec![9],
            },
        ],
        profiles: vec![vec![p.clone(), z.clone()], vec![p, z]],
        source_files: vec![json!({"bytes":100,"fnv1a64":"0000000000000000"})],
        registry: None,
        context_complete: true,
        candidate_universe: "finite-explicit-layout-only".into(),
    }
}
fn sample_fixture() -> sample::SampleIndex {
    sample::SampleIndex {
        version: 1,
        panel: identity(),
        count_policy: COUNT_POLICY.into(),
        stats: sample::SampleStats {
            read_lengths: BTreeMap::from([(100, 1)]),
            ..Default::default()
        },
        counts: WeightedBwt::build(&BTreeMap::from([
            (tokens(1).to_vec(), 1),
            (tokens(8).to_vec(), 3),
        ]))
        .unwrap(),
    }
}
#[test]
fn coupled_objective_exhaustive_correlated_ties_zeros_residuals_and_budgets() {
    let c = fixture();
    let sample = sample_fixture();
    let p = Problem::new(&c, &sample, 1.0, 0.1).unwrap();
    assert_eq!(p.factors.len(), 3);
    assert_eq!(p.components, vec![vec![0, 1]]);
    assert!(p
        .factors
        .iter()
        .any(|f| f.observed == 0 && f.original_ids == vec![9]));
    assert!(p
        .factors
        .iter()
        .any(|f| f.observed == 3 && f.background_only && f.constant));
    let mut brute = Vec::new();
    for a in 0..2 {
        for b in 0..2 {
            let assignment = Assignment {
                version: VERSION,
                model: MODEL.into(),
                compiled_checksum: c.digest().unwrap(),
                choices: vec![a, b],
            };
            let e = evaluate(&p, &assignment).unwrap();
            let mu = 0.1 + f64::from(a == 0) + f64::from(b == 0);
            let objective = mu - mu.ln() + 0.1 + (0.1 - 3.0 * 0.1f64.ln());
            assert!((e.objective - objective).abs() < 1e-12);
            brute.push((vec![a, b], objective));
        }
    }
    let r = solve(&p, 4, 4, 1e-9).unwrap();
    let minimum = brute.iter().map(|(_, s)| *s).fold(f64::INFINITY, f64::min);
    assert!((r.incumbent.as_ref().unwrap().objective - minimum).abs() < 1e-12);
    assert_eq!(
        r.correlated_optima
            .iter()
            .map(|a| a.choices.clone())
            .collect::<Vec<_>>(),
        vec![vec![0, 1], vec![1, 0]]
    );
    assert!(r.global_objective_certified && r.correlated_optima_complete);
    assert!(!r.candidate_universe_complete && !r.sequence_emission_authorized);
    let r = solve(&p, 1, 4, 1e-9).unwrap();
    assert!(
        !r.search_exhausted
            && !r.global_objective_certified
            && !r.correlated_optima_complete
            && r.lower_bound.is_none()
    );
    let r = solve(&p, 4, 1, 1e-9).unwrap();
    assert!(r.global_objective_certified && !r.correlated_optima_complete);
    assert_eq!(r.correlated_optima.len(), 1);
    assert!(solve(&p, 0, 1, 0.0).is_err());
}
#[test]
fn physical_conflicts_connect_components_and_duplicate_descriptions_never_add_copies() {
    let mut c = fixture();
    c.profiles[1]
        .iter_mut()
        .for_each(|p| p.contributions.clear());
    for w in &mut c.layout.slots[1].alternatives {
        w.pieces[0].instance = "copy0".into();
    }
    let sample = sample_fixture();
    let p = Problem::new(&c, &sample, 1.0, 0.1).unwrap();
    assert_eq!(
        p.components,
        vec![vec![0, 1]],
        "physical conflicts connect even without shared factors"
    );
    let r = solve(&p, 4, 4, 0.0).unwrap();
    assert_eq!(r.status, "infeasible");
    assert_eq!(r.feasible_assignments, 0);
    let mut l = fixture().layout;
    let w = &mut l.slots[0].alternatives[0];
    w.pieces.push(w.pieces[0].clone());
    w.adjacencies.push(layout::Adjacency {
        from: 0,
        to: 1,
        kind: "abut".into(),
    });
    assert!(l.validate().is_err());
    // Moving a harmless cut reuses one resource but only disjoint coordinates.
    let w = &mut l.slots[0].alternatives[0];
    w.pieces[0].end = 30;
    w.pieces[1].start = 30;
    l.validate().unwrap();
    l.instances.push(l.instances[0].clone());
    assert!(l.validate().is_err());
}
#[test]
fn corrupt_layout_artifacts_and_assignment_provenance_fail_closed() {
    let c = fixture();
    let t = tempfile::tempdir().unwrap();
    let path = t.path().join("profiles.json");
    c.save(&path).unwrap();
    let (loaded, _) = Compiled::load(&path, &identity()).unwrap();
    assert_eq!(c.digest().unwrap(), loaded.digest().unwrap());
    for pointer in [
        "/payload/compiler_identity",
        "/payload/count_policy",
        "/payload/model",
        "/checksum",
    ] {
        let mut v: Value = read_json(&path).unwrap();
        *v.pointer_mut(pointer).unwrap() = json!("corrupt");
        let bad = t.path().join("bad.json");
        fs::write(&bad, serde_json::to_vec(&v).unwrap()).unwrap();
        assert!(Compiled::load(&bad, &identity()).is_err());
    }
    for kind in 0..6 {
        let mut c = fixture();
        match kind {
            0 => c.context_complete = false,
            1 => c.layout.slots[0].alternatives[0].topology = "circular".into(),
            2 => c.layout.slots[0].alternatives[0].left_endpoint = "unknown".into(),
            3 => c.profiles[0].pop().map(|_| ()).unwrap(),
            4 => c.profiles[0][0].admitted_starts[0] = 0,
            _ => c.profiles[0][0].contributions[0].totals[0][1] = 0,
        }
        assert!(c.validate().is_err());
    }
    let s = sample_fixture();
    let p = Problem::new(&c, &s, 1.0, 0.1).unwrap();
    let mut a = Assignment {
        version: 1,
        model: MODEL.into(),
        compiled_checksum: "wrong".into(),
        choices: vec![0, 1],
    };
    assert!(evaluate(&p, &a).is_err());
    a.compiled_checksum = c.digest().unwrap();
    a.choices.pop();
    assert!(evaluate(&p, &a).is_err());
    assert!(Problem::new(&c, &s, f64::INFINITY, 0.1).is_err());
    let mut s = sample_fixture();
    s.stats.read_lengths.insert(50, 1);
    assert!(Problem::new(&c, &s, 1.0, 0.1).is_err());
    let mut s = sample_fixture();
    let bytes = bincode::serde::encode_to_vec(&s.counts, bincode::config::standard()).unwrap();
    s.counts = bincode::serde::decode_from_slice(&bytes, bincode::config::standard())
        .unwrap()
        .0;
    assert!(
        Problem::new(&c, &s, 1.0, 0.1).is_err(),
        "incomplete support enumeration fails closed"
    );
}

/// Oracle recollects each literal singleton through sample::build's PUBLIC native
/// matcher (not raw-window restriction), builds the weighted BWT, then queries it.
fn oracle(
    panel: &SyngIndex,
    sequence: &[u8],
    length: usize,
    reverse: bool,
) -> BTreeMap<[u64; 3], u64> {
    if length > sequence.len() {
        return BTreeMap::new();
    }
    let t = tempfile::tempdir().unwrap();
    let path = t.path().join("singletons.fa");
    let mut text = String::new();
    for read in sequence.windows(length) {
        let r = if reverse {
            crate::graph::reverse_complement(read)
        } else {
            read.to_vec()
        };
        text += &format!(">transient\n{}\n", String::from_utf8_lossy(&r));
    }
    fs::write(&path, text).unwrap();
    let s = sample::build(panel, identity(), &[path]).unwrap();
    let support = s.counts.observed_pairs().unwrap();
    for (tokens, &count) in &support {
        assert_eq!(s.counts.count(tokens).unwrap(), count);
    }
    support
}
#[test]
fn native_fresh_walk_events_match_independent_singleton_bwt_across_multiple_joins_and_ends() {
    let sequence = dna(700, 113);
    let source_dna = vec![
        ("L#0#left".to_string(), sequence[..500].to_vec()),
        ("R#0#right".to_string(), sequence[200..].to_vec()),
    ];
    let panel = SyngIndex::build(SyncmerParams::default(), source_dna.clone().into_iter());
    let lengths = [90, 250, 700, 701];
    let full = replay::compile_walk(&panel, 700, &lengths, 37, |a, b| {
        Ok(sequence[a as usize..b as usize].to_vec())
    })
    .unwrap();
    let mut overlap = false;
    for (li, &l) in lengths.iter().enumerate() {
        let expected = oracle(&panel, &sequence, l as usize, false);
        assert_eq!(expected, oracle(&panel, &sequence, l as usize, true));
        let actual: BTreeMap<_, _> = full
            .contributions
            .iter()
            .filter(|c| c.totals[li][0] > 0)
            .map(|c| (c.tokens, c.totals[li][0]))
            .collect();
        assert_eq!(actual, expected, "L={l}");
        if l == 700 {
            overlap = actual.values().any(|&q| q >= 2);
        }
    }
    assert!(overlap, "native overlapping MEM count-two fixture");
    assert_eq!(full.admitted_starts, vec![611, 451, 1, 0]);
    assert!(full.max_crop_bp <= 37 + 700 - 1);
    // Genuine mixed sources, with three joins in a single read. Same original
    // source interval is never duplicated; ownership cut placement is irrelevant.
    let t = tempfile::tempdir().unwrap();
    let fasta = t.path().join("sources.fa");
    fs::write(
        &fasta,
        source_dna
            .iter()
            .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    let mut layout = layout::Layout {
        version: 1,
        model: LAYOUT_MODEL.into(),
        panel: identity(),
        sources: source_dna
            .iter()
            .enumerate()
            .map(|(id, (name, s))| layout::Source {
                id,
                name: name.clone(),
                length: s.len() as u64,
            })
            .collect(),
        instances: vec![
            layout::Instance {
                id: "l".into(),
                source: 0,
                start: 0,
                end: 500,
            },
            layout::Instance {
                id: "r".into(),
                source: 1,
                start: 0,
                end: 500,
            },
        ],
        slots: vec![layout::Slot {
            id: "molecule".into(),
            alternatives: vec![layout::Walk {
                id: "mosaic".into(),
                topology: "linear".into(),
                left_endpoint: "asserted-molecule-terminus".into(),
                right_endpoint: "asserted-molecule-terminus".into(),
                pieces: vec![
                    ("l", 0, 250),
                    ("r", 50, 80),
                    ("l", 280, 350),
                    ("r", 150, 500),
                ]
                .into_iter()
                .map(|(instance, start, end)| layout::Piece {
                    instance: instance.into(),
                    start,
                    end,
                    strand: "+".into(),
                })
                .collect(),
                adjacencies: (0..3)
                    .map(|i| layout::Adjacency {
                        from: i,
                        to: i + 1,
                        kind: "abut".into(),
                    })
                    .collect(),
            }],
        }],
    };
    layout.validate().unwrap();
    let sources =
        layout::Sources::open(&[fasta.to_str().unwrap().into()], &layout, &panel).unwrap();
    let w = &layout.slots[0].alternatives[0];
    assert_eq!(sources.spell(&layout, w, 0, 700).unwrap(), sequence);
    let mixed = replay::compile_walk(&panel, 700, &lengths, 29, |a, b| {
        sources.spell(&layout, w, a, b)
    })
    .unwrap();
    assert_eq!(
        serde_json::to_value(&full.contributions).unwrap(),
        serde_json::to_value(&mixed.contributions).unwrap()
    );
    // Reverse whole layout, rather than just reverse source strands in place.
    let w = &mut layout.slots[0].alternatives[0];
    w.pieces.reverse();
    for p in &mut w.pieces {
        p.strand = "-".into();
    }
    let w = &layout.slots[0].alternatives[0];
    let reversed = replay::compile_walk(&panel, 700, &lengths, 41, |a, b| {
        sources.spell(&layout, w, a, b)
    })
    .unwrap();
    assert_eq!(
        serde_json::to_value(&full.contributions).unwrap(),
        serde_json::to_value(&reversed.contributions).unwrap()
    );
    // Extension beyond source L's end changes old conditional exposure. A new
    // unrelated join destroys native source pairs; both are independently replayed.
    let old = replay::compile_walk(&panel, 500, &[250], 31, |a, b| {
        Ok(sequence[a as usize..b as usize].to_vec())
    })
    .unwrap();
    let extended: BTreeMap<_, _> = full
        .contributions
        .iter()
        .map(|c| (c.tokens, c.totals[1][0]))
        .collect();
    assert!(old
        .contributions
        .iter()
        .any(|c| extended[&c.tokens] > c.totals[0][0]));
    let mut changed = sequence[..350].to_vec();
    changed.extend(crate::graph::reverse_complement(&sequence[350..]));
    let broken = replay::compile_walk(&panel, 700, &[250], 23, |a, b| {
        Ok(changed[a as usize..b as usize].to_vec())
    })
    .unwrap();
    let actual: BTreeMap<_, _> = broken
        .contributions
        .iter()
        .map(|c| (c.tokens, c.totals[0][0]))
        .collect();
    assert_eq!(actual, oracle(&panel, &changed, 250, false));
    assert!(extended
        .iter()
        .any(|(t, q)| *q > actual.get(t).copied().unwrap_or(0)));
    // Native sample support beyond the sole admitted candidate is retained as
    // background-only residuals, not discarded for lack of a registry definition.
    let mut short_layout = layout.clone();
    let w = &mut short_layout.slots[0].alternatives[0];
    w.pieces = vec![layout::Piece {
        instance: "l".into(),
        start: 0,
        end: 500,
        strand: "+".into(),
    }];
    w.adjacencies.clear();
    let compiled = compile(
        &panel,
        identity(),
        short_layout,
        &[fasta.to_str().unwrap().into()],
        vec![250],
        None,
    )
    .unwrap();
    let reads = t.path().join("sample-native-residuals.fa");
    fs::write(
        &reads,
        sequence
            .windows(250)
            .map(|s| format!(">singleton\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    let sample = sample::build(&panel, identity(), &[reads]).unwrap();
    let problem = Problem::new(&compiled, &sample, 1.0, 0.1).unwrap();
    assert!(problem
        .factors
        .iter()
        .any(|f| f.observed > 0 && f.background_only && f.original_ids.is_empty()));
    let result = evaluate(
        &problem,
        &Assignment {
            version: VERSION,
            model: MODEL.into(),
            compiled_checksum: compiled.digest().unwrap(),
            choices: vec![0],
        },
    )
    .unwrap();
    for (f, mu) in problem.factors.iter().zip(result.rates) {
        if f.background_only {
            assert_eq!(mu, 0.1);
        }
        assert_eq!(sample.counts.count(&f.tokens).unwrap(), f.observed);
    }
}
