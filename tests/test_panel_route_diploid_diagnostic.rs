//! Frozen finite fixed-copy recipes; assessment occurs after selection is written.
#[path = "../examples/panel_route_diploid_diagnostic/mod.rs"]
mod diagnostic;
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::{SyncmerParams, SyngIndex},
};
use serde_json::{json, Value};
use std::{fs, path::Path, process::Command};
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
fn cli(root: &Path, name: &str, args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(args)
        .output()
        .unwrap();
    fs::write(root.join(format!("{name}.stdout")), &output.stdout).unwrap();
    fs::write(root.join(format!("{name}.stderr")), &output.stderr).unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(!String::from_utf8_lossy(&output.stderr).contains("++match "));
}
fn p(p: &Path) -> &str {
    p.to_str().unwrap()
}
fn canonical(dna: &[u8]) -> Vec<u8> {
    let rc = impg::graph::reverse_complement(dna);
    if dna < rc.as_slice() {
        dna.to_vec()
    } else {
        rc
    }
}
/// Full molecule equality with an injective match; no alignment/partial QV claim.
fn assessment(truth: &[Vec<u8>], query: &[Vec<u8>]) -> Value {
    let mut unused: Vec<_> = query.iter().map(|s| Some(canonical(s))).collect();
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
    json!({"exact_matched_molecules":copies,"exact_matched_bases":matched,"truth_bases":tb,"query_bases":qb,"copy_weighted_exact_truth_coverage":matched as f64/tb as f64,"copy_weighted_exact_query_coverage":if qb>0 {matched as f64/qb as f64}else{0.0},"metric":"injective full-DNA/RC molecule multiset equality, no partial alignment claim"})
}
fn fixture(root: &Path, kind: &str) -> Value {
    fs::create_dir(root).unwrap();
    // Freeze recipe artifact before any operator evaluation or sample scoring.
    genome::write_json(&root.join("recipe.json"),&json!({"kind":kind,"seed":113,"length":600,"positions":match kind {"heterozygote"=>vec![300],"near"=>vec![220,260],"far"=>vec![180,420],_=>vec![]},"read_length":150,"start_stride":15,"alternating_rc":true,"nominal_per_copy_depth":10,"background":0.1})).unwrap();
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
    let graph = routes::Graph::load(&graph_dir, &identity).unwrap();
    assert_eq!(graph.lanes.len(), sequences.len());
    let assignments: Vec<_> = (0..graph.families.len())
        .map(|f| graph.native_assignment(f).unwrap())
        .collect();
    let domain = root.join("domain.json");
    genome::write_json(&domain, &assignments).unwrap();
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
    let sample_path = sample_dir.join("sample.membwt");
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
    let out = root.join("inference");
    let summary = diagnostic::run(&prefix, &graph_dir, &sample_path, &domain, &out).unwrap();
    // Nothing above consumes assessment/truth assignments. Freeze bytes now.
    let selection_bytes = fs::read(out.join("selection.json")).unwrap();
    let pairs_bytes = fs::read(out.join("pairs.json")).unwrap();
    let rows: Vec<diagnostic::PairRow> = serde_json::from_slice(&pairs_bytes).unwrap();
    let normalized: Vec<routes::Assignment> =
        genome::read_json(&out.join("normalized-domain.json")).unwrap();
    let sample = sample::SampleIndex::load(&sample_path, &identity).unwrap();
    let mut e = routes::Evaluator::new(
        &graph_dir,
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        diagnostic::MAX_FEATURES,
        diagnostic::MAX_FEATURES,
    )
    .unwrap();
    let molecules: Vec<Vec<Vec<u8>>> = normalized
        .iter()
        .map(|a| {
            a.routes
                .iter()
                .map(|r| e.sources.spell(r, 0, r.length().unwrap()).unwrap())
                .collect()
        })
        .collect();
    let truth: Vec<_> = truth_indices
        .iter()
        .map(|&i| sequences[i].1.clone())
        .collect();
    let pair_for = |source_pair: [usize; 2]| {
        let mut pair = source_pair.map(|source| {
            normalized
                .iter()
                .position(|a| a.routes[0].segments[0].source == source)
                .unwrap()
        });
        pair.sort();
        pair
    };
    let truth_pair = pair_for(truth_indices);
    let truth_row = rows.iter().find(|r| r.pair == truth_pair).unwrap();
    let independent_min = rows
        .iter()
        .map(|r| r.relative_objective)
        .fold(f64::INFINITY, f64::min);
    let support: Vec<_> = rows
        .iter()
        .filter(|r| (r.relative_objective - independent_min).abs() <= diagnostic::TIE)
        .map(|r| r.pair)
        .collect();
    assert_eq!(summary["correlated_optima"], json!(support));
    let features: Value = genome::read_json(&out.join("features.json")).unwrap();
    let tokens: Vec<[u64; 3]> = serde_json::from_value(features["tokens"].clone()).unwrap();
    let c: Vec<u64> = serde_json::from_value(features["observed_once"].clone()).unwrap();
    for (t, &n) in tokens.iter().zip(&c) {
        assert_eq!(sample.counts.count(t).unwrap(), n);
    }
    // Independent whole-pair public BWT check (not just sum of cached haploids).
    let replay = root.join("pooled-all-windows.fa");
    fs::write(
        &replay,
        truth
            .iter()
            .flat_map(|s| s.windows(150))
            .enumerate()
            .map(|(i, s)| format!(">r{i}\n{}\n", String::from_utf8_lossy(s)))
            .collect::<String>(),
    )
    .unwrap();
    let pooled = sample::build(&panel, identity.clone(), &[replay]).unwrap();
    for (t, &q) in tokens.iter().zip(&truth_row.counts_by_length_150) {
        assert_eq!(pooled.counts.count(t).unwrap(), q);
    }
    let mut phase = Value::Null;
    if positions.len() == 2 {
        let alternative = pair_for([1, 2]);
        let row = rows.iter().find(|r| r.pair == alternative).unwrap();
        let different = truth_row
            .counts_by_length_150
            .iter()
            .zip(&row.counts_by_length_150)
            .filter(|(a, b)| a != b)
            .count();
        if kind == "far" {
            assert_eq!(different, 0);
            assert_eq!(
                support.contains(&truth_pair),
                support.contains(&alternative)
            );
        }
        phase = json!({"coupling_pair":truth_pair,"repulsion_pair":alternative,"exact_distinguishing_features":different,"operator_phase_information":different>0,"coupling_optimal":support.contains(&truth_pair),"repulsion_optimal":support.contains(&alternative),"global_label_swap_is_symmetry":true,"distinct_pair_phase_is_not_label_symmetry":true});
    }
    let assessments:Vec<_>=support.iter().map(|pair|{let query:Vec<_>=pair.iter().flat_map(|&i|molecules[i].clone()).collect();json!({"pair":pair,"coverage":assessment(&truth,&query),"allele_dosage":positions.iter().map(|&pos|json!({"position":pos,"truth_mutated_copies":truth.iter().filter(|s|s[pos]!=base[pos]).count(),"query_mutated_copies":query.iter().filter(|s|s[pos]!=base[pos]).count()})).collect::<Vec<_>>()})}).collect();
    let missing = assessment(&[base.clone(), base.clone()], &[base.clone()]);
    assert_eq!(missing["copy_weighted_exact_truth_coverage"], 0.5);
    assert_eq!(missing["copy_weighted_exact_query_coverage"], 1.0);
    let swapped = diagnostic::FixedDiploid {
        copies: [
            normalized[truth_pair[1]].clone(),
            normalized[truth_pair[0]].clone(),
        ],
    }
    .validate(&mut e)
    .unwrap();
    let original = diagnostic::FixedDiploid {
        copies: [
            normalized[truth_pair[0]].clone(),
            normalized[truth_pair[1]].clone(),
        ],
    }
    .validate(&mut e)
    .unwrap();
    assert_eq!(swapped, original);
    e.depth = 20.0;
    assert!(diagnostic::analyze(&mut e, &sample, normalized.clone(), &out).is_err());
    e.depth = 10.0;
    e.background = 0.2;
    assert!(diagnostic::analyze(&mut e, &sample, normalized.clone(), &out).is_err());
    e.background = 0.1;
    assert!(serde_json::from_value::<diagnostic::FixedDiploid>(
        json!({"copies":[normalized[0],normalized[0],normalized[0]]})
    )
    .is_err());
    let report = json!({"case":kind,"truth_pair":truth_pair,"truth_in_finite_optima":support.contains(&truth_pair),"truth_objective":truth_row.relative_objective,"finite_minimum":independent_min,"phase":phase,"selected_copy_assessments":assessments,"missing_copy_negative":missing,"actual_reads":sample.stats.reads,"actual_bases":sample.stats.bases,"actual_bases_per_diploid_base":sample.stats.bases as f64/1200.0,"actual_bases_per_haploid_reference":sample.stats.bases as f64/600.0,"inference_used_truth":false,"selection_frozen_before_assessment":true,"representation_operator_correct":true});
    genome::write_json(&root.join("assessment.json"), &report).unwrap();
    assert_eq!(
        selection_bytes,
        fs::read(out.join("selection.json")).unwrap()
    );
    assert_eq!(pairs_bytes, fs::read(out.join("pairs.json")).unwrap());
    if kind == "aa" {
        assert_eq!(normalized.len(), 1);
        assert_eq!(summary["paired_hypotheses"], 1);
        assert!(summary["pair_route_cache_hits"].as_u64().unwrap() >= 2);
        assert_eq!(
            report["selected_copy_assessments"][0]["coverage"]["exact_matched_molecules"],
            2
        );
    }
    report
}
#[test]
fn frozen_fixed_diploid_public_cli_operator_and_copy_assessment() {
    let temp = tempfile::tempdir().unwrap();
    let retained = std::env::var_os("IMPG_TEST_DIPLOID_OUTPUT").map(std::path::PathBuf::from);
    if let Some(root) = &retained {
        fs::create_dir(root).unwrap();
    }
    let root = retained.as_deref().unwrap_or(temp.path());
    let reports: Vec<_> = ["aa", "heterozygote", "near", "far"]
        .iter()
        .map(|kind| fixture(&root.join(kind), kind))
        .collect();
    genome::write_json(&root.join("assessments.json"), &reports).unwrap();
}
#[test]
fn fixed_diploid_math_capacity_schema_and_resource_negatives() {
    assert!(diagnostic::add(u64::MAX, 1).is_err());
    assert!(diagnostic::joint_loss(&[u64::MAX], &[1], 62).is_err());
    assert!(diagnostic::joint_loss(&[1], &[1], u64::MAX).is_err());
    assert!(diagnostic::state_charge(diagnostic::MAX_CELLS + 1).is_err());
    assert!(diagnostic::state_charge(diagnostic::MAX_CELLS).is_err());
    assert!(diagnostic::state_charge(1000).unwrap() < diagnostic::MAX_BYTES);
    let joint = diagnostic::joint_loss(&[300, 150, 0], &[7, 0, 9], 62).unwrap();
    let independent = 20.0 - 7.0 * 201f64.ln() + 10.0;
    diagnostic::close(joint, independent).unwrap();
    let separate = 2.0 * diagnostic::joint_loss(&[150, 75, 0], &[7, 0, 9], 62).unwrap();
    assert!((joint - separate).abs() > 1.0); // Cannot add two haploid losses/C.
    assert_eq!(diagnostic::joint_loss(&[0], &[9], 62).unwrap(), 0.0);
    assert_eq!(diagnostic::joint_loss(&[150], &[0], 62).unwrap(), 10.0);
    let aa = assessment(&[b"ACG".to_vec(), b"ACG".to_vec()], &[b"CGT".to_vec()]);
    assert_eq!(aa["copy_weighted_exact_truth_coverage"], 0.5);
}

#[test]
#[ignore = "explicit unchanged finite513/coupled64x64 fixtures required"]
fn configured_fixed_copy_junction_terminal_reverse_and_cross_molecule_capacity() {
    use impg::genome_inference::joint;
    let output = std::path::PathBuf::from(
        std::env::var_os("IMPG_TEST_DIPLOID_CONTROL_OUTPUT").expect("control output required"),
    );
    fs::create_dir(&output).unwrap();
    let mut reports = Vec::new();
    for (name, var) in [
        ("finite", "IMPG_TEST_DIPLOID_FINITE"),
        ("coupled", "IMPG_TEST_DIPLOID_COUPLED"),
    ] {
        let root =
            std::path::PathBuf::from(std::env::var_os(var).expect("frozen fixture required"));
        let envelope: Value = genome::read_json(&root.join("finite/joint-profiles.json")).unwrap();
        let compiled: joint::Compiled =
            serde_json::from_value(envelope["payload"].clone()).unwrap();
        assert_eq!(envelope["checksum"], compiled.digest().unwrap());
        let sizes: Vec<_> = compiled
            .layout
            .slots
            .iter()
            .map(|s| s.alternatives.len())
            .collect();
        assert_eq!(
            sizes,
            if name == "finite" {
                vec![513]
            } else {
                vec![64, 64]
            }
        );
        let prefix = root.join(if name == "finite" {
            "panel.syng"
        } else {
            "p.syng"
        });
        let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
        let panel = SyngIndex::load(p(&prefix), Default::default()).unwrap();
        let graph = routes::Graph::load(&root.join("routes"), &identity).unwrap();
        let sample = sample::SampleIndex::load(&root.join("sample.membwt"), &identity).unwrap();
        let mut e = routes::Evaluator::new(
            &root.join("routes"),
            &graph,
            &panel,
            &sample,
            10.0,
            0.1,
            diagnostic::MAX_FEATURES,
            diagnostic::MAX_FEATURES,
        )
        .unwrap();
        let physical = |choices: &[usize]| {
            let routes: Vec<_> = choices
                .iter()
                .enumerate()
                .map(|(s, &a)| {
                    routes::Route {
                        segments: compiled.layout.slots[s].alternatives[a]
                            .pieces
                            .iter()
                            .map(|piece| routes::Segment {
                                source: compiled
                                    .layout
                                    .instances
                                    .iter()
                                    .find(|i| i.id == piece.instance)
                                    .unwrap()
                                    .source,
                                start: piece.start,
                                end: piece.end,
                                reverse: piece.strand == "-",
                            })
                            .collect(),
                    }
                    .normalize()
                })
                .collect();
            let family = graph
                .families
                .iter()
                .position(|f| {
                    f.paths.len() == routes.len()
                        && f.paths.iter().zip(&routes).all(|(&s, r)| {
                            r.segments[0].source == s && r.segments.last().unwrap().source == s
                        })
                })
                .unwrap();
            routes::Assignment {
                version: 1,
                model: routes::MODEL.into(),
                graph_checksum: graph.digest().unwrap(),
                family,
                routes,
            }
        };
        let mut feasible = Vec::new();
        let mut invalid_cross = None;
        for ordinal in 0..sizes.iter().product::<usize>() {
            let choices = if sizes.len() == 1 {
                vec![ordinal]
            } else {
                vec![ordinal / sizes[1], ordinal % sizes[1]]
            };
            let a = physical(&choices);
            match e.validate_assignment(&a) {
                Ok(a) => feasible.push(a),
                Err(error) => {
                    if invalid_cross.is_none() {
                        invalid_cross = Some((a, error.to_string()));
                    }
                }
            }
        }
        let a = if name == "finite" {
            feasible
                .iter()
                .find(|a| a.routes[0].segments.len() >= 5)
                .expect("existing multi-junction control")
                .clone()
        } else {
            feasible[0].clone()
        };
        let pair = diagnostic::FixedDiploid {
            copies: [a.clone(), a.clone()],
        }
        .validate(&mut e)
        .unwrap();
        let haploid = e.evaluate(&a).unwrap();
        let mut doubled = vec![diagnostic::Counts::new(); graph.read_lengths.len()];
        for copy in &pair.copies {
            for route in &copy.routes {
                let (q, _) = e.route_counts(route).unwrap();
                for (li, counts) in q.iter().enumerate() {
                    for (&t, &v) in counts {
                        let old = doubled[li].get(&t).copied().unwrap_or(0);
                        doubled[li].insert(t, diagnostic::add(old, v).unwrap());
                    }
                }
            }
        }
        let mut count_cells = 0;
        for (li, &length) in graph.read_lengths.iter().enumerate() {
            let reads = output.join(format!("{name}-L{length}-pooled.fa"));
            let mut text = String::new();
            for copy in &pair.copies {
                for route in &copy.routes {
                    let sequence = e.sources.spell(route, 0, route.length().unwrap()).unwrap();
                    for read in sequence.windows(length as usize) {
                        text += &format!(">anonymous\n{}\n", String::from_utf8_lossy(read));
                    }
                }
            }
            // L500 has no admitted starts here. The public FASTA reader rejects
            // empty files (FileTooShort); check the empty physical domain directly
            // rather than inventing a read or claiming an empty BWT was built.
            let empty = text.is_empty();
            fs::write(&reads, text).unwrap();
            if empty {
                assert!(pair
                    .copies
                    .iter()
                    .flat_map(|a| &a.routes)
                    .all(|r| r.length().unwrap() < length));
                assert!(doubled[li].is_empty());
            } else {
                let replay = sample::build(&panel, identity.clone(), &[reads]).unwrap();
                assert_eq!(doubled[li], replay.counts.observed_pairs().unwrap());
            }
            for factor in &haploid.factors {
                assert_eq!(
                    doubled[li].get(&factor.tokens).copied().unwrap_or(0),
                    2 * factor.counts_by_length[li]
                );
                count_cells += 1;
            }
        }
        let rejection = if let Some((bad, error)) = invalid_cross {
            assert!(error.contains("conflicting canonical source-span reuse"));
            assert!(diagnostic::FixedDiploid {
                copies: [a.clone(), bad]
            }
            .validate(&mut e)
            .is_err());
            Some(error)
        } else {
            None
        };
        if name == "coupled" {
            assert!(rejection.is_some());
        }
        let mut bad_route = a.routes[0].clone();
        bad_route.segments.push(bad_route.segments[0].clone());
        let route_error = e.route_counts(&bad_route).unwrap_err().to_string();
        assert!(route_error.contains("overlapping route profile resources"));
        // Harmless segmentation is one description, not an additional copy.
        let native = graph.native_assignment(a.family).unwrap();
        let mut split = native.clone();
        let whole = split.routes[0].segments[0].clone();
        let cut = (whole.start + whole.end) / 2;
        split.routes[0].segments = vec![
            routes::Segment {
                end: cut,
                ..whole.clone()
            },
            routes::Segment {
                start: cut,
                ..whole
            },
        ];
        assert_eq!(e.validate_assignment(&split).unwrap(), native);
        reports.push(json!({"fixture":name,"paired_assignment":pair,"haploid_admitted_starts":haploid.admitted_starts,"haploid_corrected_starts":haploid.corrected_starts,"all_length_double_count_cells_checked":count_cells,"exact_public_pooled_replay":true,"within_copy_cross_molecule_rejection":rejection,"within_route_repeat_rejection":route_error,"same_source_split_normalizes":true,"feasible_haploid_assignments":feasible.len()}));
    }
    genome::write_json(&output.join("controls.json"), &reports).unwrap();
}

#[test]
#[ignore = "explicit existing B2 opposite fixture, no new search"]
fn configured_fixed_copy_existing_reversed_interior_public_replay() {
    use std::io::BufRead;
    let root = std::path::PathBuf::from(
        std::env::var_os("IMPG_TEST_DIPLOID_OPPOSITE")
            .expect("existing B2 opposite fixture required"),
    );
    let out = std::path::PathBuf::from(
        std::env::var_os("IMPG_TEST_DIPLOID_OPPOSITE_OUTPUT").expect("fresh output required"),
    );
    fs::create_dir(&out).unwrap();
    let prefix = root.join("panel.syng");
    let identity = genome::PanelIdentity::read(p(&prefix)).unwrap();
    let panel = SyngIndex::load(p(&prefix), Default::default()).unwrap();
    let graph = routes::Graph::load(&root.join("routes"), &identity).unwrap();
    let sample = sample::SampleIndex::load(&root.join("sample.membwt"), &identity).unwrap();
    let mut e = routes::Evaluator::new(
        &root.join("routes"),
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        diagnostic::MAX_FEATURES,
        diagnostic::MAX_FEATURES,
    )
    .unwrap();
    let ledger = root.join("automatic/scores.jsonl");
    let mut selected = None;
    for (ordinal, line) in std::io::BufReader::new(fs::File::open(&ledger).unwrap())
        .lines()
        .enumerate()
    {
        let value: Value = serde_json::from_str(&line.unwrap()).unwrap();
        let Some(a) = value.get("assignment") else {
            continue;
        };
        let a: routes::Assignment = serde_json::from_value(a.clone()).unwrap();
        if !a.routes.iter().any(|r| {
            r.segments.len() >= 3
                && r.segments[1..r.segments.len() - 1]
                    .iter()
                    .any(|s| s.reverse && s.source != r.segments[0].source)
        }) {
            continue;
        }
        if let Ok(a) = e.validate_assignment(&a) {
            selected = Some((ordinal, a));
            break;
        }
    }
    let (ordinal, a) = selected.expect(
        "no existing publicly valid B2 reversed-interior assignment; do not invent geometry",
    );
    genome::write_json(&out.join("selected-before-paired-objective.json"),&json!({"ledger":ledger,"zero_based_ordinal":ordinal,"selection":"first qualifying public-valid reversed-interior/source-switch assignment, no score selection","assignment":a})).unwrap();
    let summary = diagnostic::analyze(&mut e, &sample, vec![a], &out).unwrap();
    assert_eq!(summary["paired_hypotheses"], 1);
    assert_eq!(summary["aa_exact_double_checks"], 1);
}
