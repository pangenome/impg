use super::super::{catalog, sample};
use super::*;
use crate::syng::SyncmerParams;
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
        checksum_algorithm: "synthetic".into(),
        sidecars: vec![],
    }
}
fn setup(
    sequences: Vec<(String, Vec<u8>)>,
    lengths: Vec<u64>,
    core: u64,
) -> (tempfile::TempDir, SyngIndex, Graph, sample::SampleIndex) {
    let temp = tempfile::tempdir().unwrap();
    let panel = SyngIndex::build(SyncmerParams::default(), sequences.clone().into_iter());
    let mut text = String::new();
    for (n, s) in &sequences {
        text += &format!(">{n}\n{}\n", String::from_utf8_lossy(s));
    }
    let path = temp.path().join("sources.fa");
    fs::write(&path, text).unwrap();
    let groups = panel
        .name_map
        .path_to_name
        .iter()
        .enumerate()
        .flat_map(|(id, name)| {
            let n = panel.name_map.path_to_length[id];
            [
                catalog::GroupInput {
                    id: format!("different-label-{id}-left"),
                    scaffold: None,
                    occurrences: vec![
                        catalog::Interval {
                            path: name.clone(),
                            start: 0,
                            end: n / 2,
                            strand: None,
                        },
                        catalog::Interval {
                            path: name.clone(),
                            start: 0,
                            end: n / 2,
                            strand: None,
                        },
                    ],
                },
                catalog::GroupInput {
                    id: format!("different-label-{id}-right"),
                    scaffold: None,
                    occurrences: vec![catalog::Interval {
                        path: name.clone(),
                        start: n / 2,
                        end: n,
                        strand: Some("-".into()),
                    }],
                },
            ]
        })
        .collect();
    let catalog = catalog::build(
        &panel,
        identity(),
        catalog::CatalogInput { version: 1, groups },
    )
    .unwrap();
    let catpath = temp.path().join("catalog.json");
    super::super::save_catalog(&catpath, &catalog).unwrap();
    let out = temp.path().join("routes");
    super::super::with_output_model(&out, MODEL, || {
        build(
            &panel,
            identity(),
            &catpath,
            &[path.to_str().unwrap().into()],
            lengths.clone(),
            &out,
            core,
            100000,
        )
    })
    .unwrap();
    let graph = Graph::load(&out, &identity()).unwrap();
    let read = temp.path().join("read.fa");
    let n = lengths[0] as usize;
    let s = &sequences[0].1;
    fs::write(
        &read,
        format!(">r\n{}\n", String::from_utf8_lossy(&s[..n.min(s.len())])),
    )
    .unwrap();
    let sample = sample::build(&panel, identity(), &[read]).unwrap();
    (temp, panel, graph, sample)
}
#[test]
fn automatic_inventory_ownership_views_inversions_zero_anchors_and_storage_segmentation() {
    let a = dna(400, 77);
    let mut inversion = a.clone();
    inversion[100..300].copy_from_slice(&crate::graph::reverse_complement(&a[100..300]));
    let sequences = vec![
        ("A#0#chr".into(), a),
        ("A#0#zero".into(), vec![b'N'; 50]),
        ("B#0#different-label-inversion".into(), inversion),
    ];
    let (t, p, g, s) = setup(sequences, vec![90, 250, 500], 37);
    assert_eq!(g.lanes.len(), 3);
    assert_eq!(g.families.len(), 2);
    assert_eq!(g.total_source_bp, 850);
    assert!(g.lanes.iter().any(|l| l.port_count == 0 && l.length == 50));
    let mut e = Evaluator::new(
        &t.path().join("routes"),
        &g,
        &p,
        &s,
        1.0,
        0.1,
        100000,
        100000,
    )
    .unwrap();
    for f in 0..g.families.len() {
        let a = g.native_assignment(f).unwrap();
        e.validate_assignment(&a).unwrap();
        e.evaluate(&a).unwrap();
    }
    let out = t.path().join("recored");
    super::super::with_output_model(&out, MODEL, || {
        build(
            &p,
            identity(),
            &t.path().join("catalog.json"),
            &g.source_paths,
            g.read_lengths.clone(),
            &out,
            113,
            100000,
        )
    })
    .unwrap();
    let other = Graph::load(&out, &identity()).unwrap();
    assert_eq!(g.port_count, other.port_count);
    assert_eq!(g.ports.hash, other.ports.hash);
    assert_eq!(g.registry.hash, other.registry.hash);
    for (a, b) in g.lanes.iter().zip(&other.lanes) {
        assert_eq!(a.ports.hash, b.ports.hash);
        for (x, y) in a.profiles.iter().zip(&b.profiles) {
            assert_eq!(x.totals.hash, y.totals.hash);
        }
    }
    // Every physical anchor emits exactly two opposite observer ports, and the
    // reverse cut is p+k-floor(k/2), not p+floor(k/2) for odd k.
    let mut r =
        BufReader::new(File::open(t.path().join("routes").join(&g.lanes[0].ports.path)).unwrap());
    while let Some(f) = storage::Port::read(&mut r, g.k as usize).unwrap() {
        let rev = storage::Port::read(&mut r, g.k as usize).unwrap().unwrap();
        assert!(!f.reverse && rev.reverse);
        assert_eq!(f.anchor, rev.anchor);
        assert_eq!(f.word, crate::graph::reverse_complement(&rev.word));
        assert_eq!(rev.cut(g.k) - f.cut(g.k), g.k % 2);
    }
}
#[test]
fn native_index_plus_lazy_corrections_matches_public_singleton_bwt_all_orientations_and_ends() {
    let whole = dna(700, 113);
    let (t, p, g, s) = setup(
        vec![
            ("L#0#left".into(), whole[..500].to_vec()),
            ("R#0#right".into(), whole[200..].to_vec()),
        ],
        vec![90, 250, 700, 701],
        29,
    );
    let mut e = Evaluator::new(
        &t.path().join("routes"),
        &g,
        &p,
        &s,
        1.0,
        0.1,
        100000,
        100000,
    )
    .unwrap();
    let routes = vec![
        Route {
            segments: vec![
                Segment {
                    source: 0,
                    start: 0,
                    end: 250,
                    reverse: false,
                },
                Segment {
                    source: 1,
                    start: 50,
                    end: 80,
                    reverse: false,
                },
                Segment {
                    source: 0,
                    start: 280,
                    end: 350,
                    reverse: false,
                },
                Segment {
                    source: 1,
                    start: 150,
                    end: 500,
                    reverse: false,
                },
            ],
        },
        Route {
            segments: vec![
                Segment {
                    source: 0,
                    start: 100,
                    end: 300,
                    reverse: true,
                },
                Segment {
                    source: 1,
                    start: 20,
                    end: 500,
                    reverse: false,
                },
            ],
        },
        Route {
            segments: vec![Segment {
                source: 0,
                start: 100,
                end: 200,
                reverse: false,
            }],
        },
        Route {
            segments: vec![Segment {
                source: 0,
                start: 0,
                end: 500,
                reverse: true,
            }],
        },
    ];
    for (ri, route) in routes.iter().enumerate() {
        let sequence = e.sources.spell(route, 0, route.length().unwrap()).unwrap();
        let (q, corrected) = e.route_counts(route).unwrap();
        assert!(corrected > 0 || route.segments.len() == 1);
        for (li, &length) in g.read_lengths.iter().enumerate() {
            let mut text = String::new();
            if length <= sequence.len() as u64 {
                for read in sequence.windows(length as usize) {
                    text += &format!(">singleton\n{}\n", String::from_utf8_lossy(read));
                }
            }
            let expected = if text.is_empty() {
                BTreeMap::new()
            } else {
                let path = t.path().join(format!("oracle-{ri}-{li}.fa"));
                fs::write(&path, text).unwrap();
                let sample = sample::build(&p, identity(), &[path]).unwrap();
                let support = sample.counts.observed_pairs().unwrap();
                for (t, &n) in &support {
                    assert_eq!(sample.counts.count(t).unwrap(), n);
                }
                support
            };
            assert_eq!(q[li], expected, "route {ri}, L={length}");
            if ri == 0 && length == 700 {
                assert!(
                    q[li].values().any(|&q| q >= 2),
                    "native overlapping MEM count-two"
                );
            }
        }
    }
    let seg = Segment {
        source: 0,
        start: 10,
        end: 20,
        reverse: true,
    };
    assert_eq!(evaluate::native_domain(&seg, 0, 10, 5, 500), (6, 16));
    assert_eq!(evaluate::native_domain(&seg, 0, 10, 600, 500), (0, 0));
}
#[test]
fn resource_overlap_orientation_capacity_and_outward_bound() {
    let p = Segment {
        source: 0,
        start: 0,
        end: 10,
        reverse: false,
    };
    let mut q = p.clone();
    q.reverse = true;
    assert!(!evaluate::resources_feasible([p.clone(), q.clone()]));
    q.start = 10;
    q.end = 20;
    assert!(evaluate::resources_feasible([p, q]));
    let (t, p, g, s) = setup(
        vec![
            ("A#0#chr".into(), dna(180, 53)),
            ("B#0#chr".into(), dna(180, 59)),
        ],
        vec![90],
        37,
    );
    let mut e = Evaluator::new(
        &t.path().join("routes"),
        &g,
        &p,
        &s,
        10.0,
        0.1,
        100000,
        100000,
    )
    .unwrap();
    let lower = e.lower_bound().unwrap();
    for f in 0..2 {
        assert!(
            lower
                <= e.evaluate(&g.native_assignment(f).unwrap())
                    .unwrap()
                    .relative_objective
        );
    }
    let mut bad = g.native_assignment(0).unwrap();
    bad.graph_checksum = "corrupt".into();
    assert!(e.evaluate(&bad).is_err());
    let mut bad = g.native_assignment(0).unwrap();
    bad.routes[0].segments[0].end -= 1;
    assert!(e.evaluate(&bad).is_err());
    let out = t.path().join("search");
    fs::create_dir(&out).unwrap();
    let r = search(
        &mut e,
        SearchBudget {
            max_work: 1,
            max_evaluations: 10,
            max_frontier: 20,
            max_optima: 10,
            tie_epsilon: 1e-9,
        },
        &out,
    )
    .unwrap();
    assert!(
        !r.native_initializations_complete
            && !r.global_optimum_certified
            && !r.correlated_optima_complete
            && !r.sequence_emission_authorized
    );
}
#[test]
fn oriented_switches_disjoint_reuse_and_two_slot_counts_once_inside_log() {
    let sequence = dna(500, 113);
    let reverse = crate::graph::reverse_complement(&sequence);
    let (t, p, g, s) = setup(
        vec![
            ("A#0#one".into(), sequence.clone()),
            ("A#0#two".into(), sequence.clone()),
            ("B#0#reverse".into(), reverse),
        ],
        vec![90, 250],
        37,
    );
    let root = t.path().join("routes");
    let mut r = BufReader::new(File::open(root.join(&g.lanes[0].ports.path)).unwrap());
    let mut cuts = Vec::new();
    while let Some(port) = storage::Port::read(&mut r, g.k as usize).unwrap() {
        if !port.reverse {
            cuts.push(port.cut(g.k));
        }
    }
    assert!(cuts.len() >= 4);
    let mut e = Evaluator::new(&root, &g, &p, &s, 1.0, 0.1, 100000, 100000).unwrap();
    let family = g.lanes[0].family;
    let route = |source: usize, left: u64, right: u64| Route {
        segments: vec![
            Segment {
                source,
                start: 0,
                end: left,
                reverse: false,
            },
            Segment {
                source: 2,
                start: 500 - right,
                end: 500 - left,
                reverse: true,
            },
            Segment {
                source,
                start: right,
                end: 500,
                reverse: false,
            },
        ],
    };
    let mut assignment = g.native_assignment(family).unwrap();
    assignment.routes = vec![route(0, cuts[0], cuts[1]), route(1, cuts[2], cuts[3])];
    let result = e.evaluate(&assignment).unwrap();
    assert_eq!(result.assignment.routes.len(), 2);
    for (li, &length) in g.read_lengths.iter().enumerate() {
        let path = t.path().join(format!("two-copies-{li}.fa"));
        fs::write(
            &path,
            sequence
                .windows(length as usize)
                .map(|s| format!(">singleton\n{}\n", String::from_utf8_lossy(s)))
                .collect::<String>(),
        )
        .unwrap();
        let oracle = sample::build(&p, identity(), &[path]).unwrap();
        for f in &result.factors {
            assert_eq!(
                f.counts_by_length[li],
                2 * oracle.counts.count(&f.tokens).unwrap()
            );
            if li == 0 {
                let signal = f.counts_by_length[li] as f64 / 90.0;
                assert_eq!(signal, f.signal);
                assert!(
                    (f.relative_loss - (signal - f.observed as f64 * (signal / 0.1).ln_1p())).abs()
                        < 1e-12
                );
            }
        }
    }
    // Cache computations are reusable, but the fixed score and physical count
    // multiplicities cannot depend on which alternative was visited previously.
    e.evaluate(&g.native_assignment(g.lanes[2].family).unwrap())
        .unwrap();
    assert_eq!(
        e.evaluate(&assignment).unwrap().relative_objective,
        result.relative_objective
    );
    let mut wrong = assignment.clone();
    wrong.routes[0].segments[1].start += 1;
    assert!(
        e.evaluate(&wrong).is_err(),
        "odd-k reverse cut shift must be rejected"
    );
    assignment.routes[1] = route(1, cuts[0], cuts[1]);
    assert!(
        e.evaluate(&assignment).is_err(),
        "opposite traversal uses the same span capacity across genome"
    );
}
#[test]
fn repeat_hubs_keep_every_physical_port_and_budget_caps_do_not_change_rule() {
    let motif = dna(250, 113);
    let mut repeated = motif.clone();
    repeated.extend(&motif);
    let (t, p, g, s) = setup(
        vec![
            ("A#0#repeat".into(), repeated),
            ("B#0#single".into(), motif),
        ],
        vec![90],
        31,
    );
    let root = t.path().join("routes");
    let mut source = File::open(root.join(&g.lanes[0].ports.path)).unwrap();
    let port = storage::port_at(&mut source, g.k as usize, 0).unwrap();
    let mut all = File::open(root.join(&g.ports.path)).unwrap();
    let (lo, hi) = storage::bucket(&mut all, g.k as usize, g.port_count, &port.word).unwrap();
    let mut spans = BTreeSet::new();
    for i in lo..hi {
        let p = storage::port_at(&mut all, g.k as usize, i).unwrap();
        spans.insert((p.source, p.anchor, p.reverse));
    }
    assert!(spans.iter().filter(|(s, _, _)| *s == 0).count() >= 2);
    assert!(spans.iter().any(|(s, _, _)| *s == 1));
    assert_eq!(spans.len() as u64, hi - lo);
    let mut e = Evaluator::new(&root, &g, &p, &s, 1.0, 0.1, 100000, 0).unwrap();
    for (name, budget) in [
        (
            "frontier",
            SearchBudget {
                max_work: 100000,
                max_evaluations: 10000,
                max_frontier: 1,
                max_optima: 1,
                tie_epsilon: 1e-9,
            },
        ),
        (
            "evaluation",
            SearchBudget {
                max_work: 100000,
                max_evaluations: 1,
                max_frontier: 100,
                max_optima: 1,
                tie_epsilon: 1e-9,
            },
        ),
    ] {
        let out = t.path().join(name);
        fs::create_dir(&out).unwrap();
        let result = search(&mut e, budget, &out).unwrap();
        assert!(
            !result.search_exhausted
                && !result.global_optimum_certified
                && !result.correlated_optima_complete
        );
        assert!(result.generation_rule_complete);
        assert!(out.join("frontier.json").is_file());
        assert!(result.status.contains(name));
    }
    let small = t.path().join("small-cap");
    let result = super::super::with_output_model(&small, MODEL, || {
        build(
            &p,
            identity(),
            &t.path().join("catalog.json"),
            &g.source_paths,
            vec![90],
            &small,
            31,
            1,
        )
    });
    assert!(result.is_err());
    let manifest: serde_json::Value = read_json(&small.join("manifest.json")).unwrap();
    assert_eq!(manifest["status"], "failed");
    assert!(Graph::load(&small, &identity()).is_err());
}

#[test]
fn corrupt_graph_histogram_source_and_partial_support_fail_closed() {
    let (t, p, g, s) = setup(vec![("A#0#chr".into(), dna(300, 113))], vec![90], 37);
    for case in 0..7 {
        let mut bad = g.clone();
        match case {
            0 => bad.generation_complete = false,
            1 => bad.cut_offset += 1,
            2 => bad.port_count += 1,
            3 => bad.families[0].paths.clear(),
            4 => bad.lanes[0].profiles[0].length += 1,
            5 => bad.lanes[0].profiles[0].index.bytes += 1,
            _ => bad.compiler_identity = "unknown producer".into(),
        }
        assert!(bad.validate().is_err());
    }
    let root = t.path().join("routes");
    let bytes = bincode::serde::encode_to_vec(&s, bincode::config::standard()).unwrap();
    let (mut incomplete, _): (sample::SampleIndex, usize) =
        bincode::serde::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
    assert!(Evaluator::new(&root, &g, &p, &incomplete, 1.0, 0.1, 100000, 0).is_err());
    incomplete.counts.rebuild().unwrap();
    incomplete.stats.read_lengths.insert(91, 1);
    assert!(Evaluator::new(&root, &g, &p, &incomplete, 1.0, 0.1, 100000, 0).is_err());
    let source = Path::new(&g.source_paths[0]);
    let before = fs::read(source).unwrap();
    let mut changed = before.clone();
    let last = changed.len() - 2;
    changed[last] = if changed[last] == b'A' { b'C' } else { b'A' };
    fs::write(source, changed).unwrap();
    assert!(Evaluator::new(&root, &g, &p, &s, 1.0, 0.1, 100000, 0).is_err());
    fs::write(source, before).unwrap();
    let path = root.join("graph.json");
    let mut value: serde_json::Value = read_json(&path).unwrap();
    value["read_lengths"] = serde_json::json!([91]);
    fs::write(&path, serde_json::to_vec(&value).unwrap()).unwrap();
    assert!(Graph::load(&root, &identity()).is_err());
}

#[test]
fn review_fai_offset_swap_with_unchanged_fasta_must_fail_before_source_open() {
    let (t, p, g, s) = setup(
        vec![
            ("A#0#one".into(), dna(300, 113)),
            ("B#0#two".into(), dna(300, 59)),
        ],
        vec![90],
        37,
    );
    let source = Path::new(&g.source_paths[0]);
    let before = fs::read(source).unwrap();
    let fai = PathBuf::from(format!("{}.fai", source.display()));
    let original = fs::read_to_string(&fai).unwrap();
    let mut rows: Vec<Vec<String>> = original
        .lines()
        .map(|l| l.split('\t').map(String::from).collect())
        .collect();
    assert_eq!(rows[0][1], rows[1][1]);
    let offset = rows[0][2].clone();
    rows[0][2] = rows[1][2].clone();
    rows[1][2] = offset;
    fs::write(
        &fai,
        rows.iter().map(|r| r.join("\t") + "\n").collect::<String>(),
    )
    .unwrap();
    assert_eq!(fs::read(source).unwrap(), before);
    assert!(
        Evaluator::new(&t.path().join("routes"), &g, &p, &s, 1.0, 0.1, 100000, 0).is_err(),
        "changed FAI offsets were adopted for cached native profiles"
    );
    fs::write(&fai, &original).unwrap();
    fs::remove_file(&fai).unwrap();
    assert!(Evaluator::new(&t.path().join("routes"), &g, &p, &s, 1.0, 0.1, 100000, 0).is_err());
    assert!(
        !fai.exists(),
        "missing FAI was silently recreated during reuse"
    );
    fs::write(&fai, original).unwrap();
    Evaluator::new(&t.path().join("routes"), &g, &p, &s, 1.0, 0.1, 100000, 0).unwrap();
}

#[test]
fn review_single_source_repeat_is_not_a_mixed_source_assignment() {
    let motif = dna(250, 113);
    let mut repeated = motif.clone();
    repeated.extend(motif);
    let (t, p, g, s) = setup(vec![("A#0#repeat".into(), repeated)], vec![90], 31);
    let mut e = Evaluator::new(&t.path().join("routes"), &g, &p, &s, 1.0, 0.1, 100000, 0).unwrap();
    let mut cut = g.native_assignment(0).unwrap();
    let mut tail = cut.routes[0].segments[0].clone();
    tail.start = 250;
    cut.routes[0].segments[0].end = 250;
    cut.routes[0].segments.push(tail);
    assert_eq!(
        e.evaluate(&cut).unwrap().assignment.routes[0]
            .segments
            .len(),
        1,
        "harmless segmentation must normalize before switch metrics"
    );
    let out = t.path().join("search-repeat");
    fs::create_dir(&out).unwrap();
    let r = search(
        &mut e,
        SearchBudget {
            max_work: 100000,
            max_evaluations: 10000,
            max_frontier: 1000,
            max_optima: 1000,
            tie_epsilon: 1e-9,
        },
        &out,
    )
    .unwrap();
    assert!(
        r.maximum_switches_evaluated > 0,
        "fixture must explore a same-source repeat rearrangement"
    );
    assert_eq!(
        r.mixed_assignments_evaluated, 0,
        "same-source rearrangement passed the cross-source operational gate"
    );
    assert_eq!(r.mixed_identity_assignments_evaluated, 0);
    assert!(r.non_native_assignments_evaluated > 0);
    let native = fs::read_to_string(out.join("evaluations.jsonl"))
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_string();
    let row: serde_json::Value = serde_json::from_str(&native).unwrap();
    assert_eq!(row["switches"], 0);
    assert_eq!(row["cross_source_route"], false);
    assert_eq!(row["cross_identity_route"], false);
}

#[test]
fn review_bgzf_fai_and_gzi_are_bound_after_creation_and_not_recreated_on_reuse() {
    let (t, p, original, s) = setup(
        vec![
            ("A#0#one".into(), dna(300, 113)),
            ("B#0#two".into(), dna(300, 59)),
        ],
        vec![90],
        37,
    );
    let compressed = t.path().join("sources.fa.gz");
    {
        let mut writer = rust_htslib::bgzf::Writer::from_path(&compressed).unwrap();
        writer
            .write_all(&fs::read(&original.source_paths[0]).unwrap())
            .unwrap();
    }
    let out = t.path().join("bgzf-routes");
    super::super::with_output_model(&out, MODEL, || {
        build(
            &p,
            identity(),
            &t.path().join("catalog.json"),
            &[compressed.to_str().unwrap().into()],
            vec![90],
            &out,
            37,
            100000,
        )
    })
    .unwrap();
    let g = Graph::load(&out, &identity()).unwrap();
    assert_eq!(g.source_paths.len(), 1);
    assert_eq!(
        g.source_access
            .iter()
            .map(|b| b.path.clone())
            .collect::<Vec<_>>(),
        vec![
            format!("{}.fai", compressed.display()),
            format!("{}.gzi", compressed.display())
        ]
    );
    let data = fs::read(&compressed).unwrap();
    for binding in &g.source_access {
        let path = Path::new(&binding.path);
        let bytes = fs::read(path).unwrap();
        let mut corrupt = bytes.clone();
        corrupt[0] ^= 1;
        fs::write(path, corrupt).unwrap();
        assert!(Evaluator::new(&out, &g, &p, &s, 1.0, 0.1, 100000, 0).is_err());
        assert_eq!(fs::read(&compressed).unwrap(), data);
        fs::remove_file(path).unwrap();
        assert!(Evaluator::new(&out, &g, &p, &s, 1.0, 0.1, 100000, 0).is_err());
        assert!(
            !path.exists(),
            "missing access sidecar was recreated and adopted"
        );
        fs::write(path, bytes).unwrap();
    }
    Evaluator::new(&out, &g, &p, &s, 1.0, 0.1, 100000, 0).unwrap();
}
