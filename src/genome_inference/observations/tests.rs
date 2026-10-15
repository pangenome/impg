use super::super::sample::{collect_raw_views, collect_tagged_read, normalize_reverse};
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
#[test]
fn native_raw_crop_event_profiles_match_exhaustive_public_singletons() {
    let a = dna(470, 17);
    let b = dna(470, 23);
    let mut repeated = a[..180].to_vec();
    repeated.extend_from_slice(b"NNNNNNNNNNNNNNNNN");
    repeated.extend_from_slice(&a[..180]);
    let mut masked = a.clone();
    masked[181..188].fill(b'N');
    let panel = SyngIndex::build(
        SyncmerParams::default(),
        vec![
            ("A#0#a".into(), a.clone()),
            ("B#0#b".into(), b),
            ("R#0#r".into(), crate::graph::reverse_complement(&a)),
            ("A#0#repeat".into(), repeated.clone()),
        ]
        .into_iter(),
    );
    let temp = tempfile::tempdir().unwrap();
    let identity = PanelIdentity {
        checksum_algorithm: "fixture".into(),
        sidecars: vec![],
    };
    let k = panel.syncmer_length_bp() as u64;
    let mut checked = 0;
    let mut crops = 0;
    let mut singleton_checks = 0;
    let mut maximum_physical_multiplicity = 0;
    let mut multicopy_reads = 0;
    let mut exact_anchor_windows = 0;
    let mut anchors_at_read_end = 0;
    for source in [&a, &masked, &repeated] {
        let views = profile::raw_views(&panel, source).unwrap();
        let mut locations = BTreeMap::new();
        for walk in &views {
            for pair in walk.windows(2) {
                locations.insert(
                    (
                        canonical(&encode_walk(pair).unwrap()),
                        pair[0].1,
                        pair[1].1 + k,
                    ),
                    (),
                );
            }
        }
        // Native selection hashes K-1 bases, but exported anchors and the
        // production full-window guard require all K bases (63 here).
        for length in [k - 1, k, k + 1, 80, 150, 240, 350, 500] {
            let mut expected: BTreeMap<(Vec<u64>, u64, u64), u64> =
                locations.keys().map(|key| (key.clone(), 0)).collect();
            if length <= source.len() as u64 {
                for start in 0..=source.len() as u64 - length {
                    let read = &source[start as usize..(start + length) as usize];
                    let raw = profile::raw_views(&panel, read).unwrap();
                    for v in 0..2 {
                        assert_eq!(raw[v], profile::restrict(&views[v], start, length, k));
                        assert!(raw[v].iter().all(|&(_, p)| p + k <= length));
                        anchors_at_read_end +=
                            raw[v].iter().filter(|&&(_, p)| p + k == length).count();
                        if length < k {
                            assert!(raw[v].is_empty());
                        }
                        if length == k {
                            assert!(raw[v].len() <= 1 && raw[v].iter().all(|&(_, p)| p == 0));
                            if read.contains(&b'N') {
                                assert!(raw[v].is_empty());
                            }
                            exact_anchor_windows += usize::from(!raw[v].is_empty());
                        }
                        crops += 1;
                    }
                    let public = collect_tagged_read(&panel, read).unwrap();
                    let optimized = collect_raw_views(&panel, &raw[0], &raw[1], length).unwrap();
                    assert_eq!(public, optimized);
                    let reverse =
                        collect_tagged_read(&panel, &crate::graph::reverse_complement(read))
                            .unwrap();
                    let reverse: BTreeSet<_> = reverse
                        .iter()
                        .map(|w| normalize_reverse(w, length, k))
                        .collect();
                    assert_eq!(public.iter().cloned().collect::<BTreeSet<_>>(), reverse);
                    let mut coordinate_counts = BTreeMap::new();
                    let mut tagged_counts = BTreeMap::new();
                    for record in &public {
                        for pair in record.windows(2) {
                            let tokens = canonical(&encode_walk(pair).unwrap());
                            *expected
                                .get_mut(&(
                                    tokens.clone(),
                                    pair[0].1 + start,
                                    pair[1].1 + start + k,
                                ))
                                .expect("raw union missing native incidence") += 1;
                            *tagged_counts
                                .entry((tokens.clone(), pair[0].1, pair[1].1))
                                .or_insert(0u64) += 1;
                            *coordinate_counts.entry(tokens).or_insert(0u64) += 1;
                        }
                    }
                    maximum_physical_multiplicity = maximum_physical_multiplicity
                        .max(tagged_counts.values().copied().max().unwrap_or(0));
                    let mut copies = BTreeMap::new();
                    for (tokens, _, _) in tagged_counts.keys() {
                        *copies.entry(tokens).or_insert(0) += 1;
                    }
                    multicopy_reads += usize::from(copies.values().any(|&n| n > 1));
                    // Actual public sample/BWT oracle on a deterministic subset;
                    // all starts above use the unchanged production collector.
                    if start % 61 == 0 {
                        let path = temp.path().join("singleton.fa");
                        fs::write(
                            &path,
                            format!(">transient\n{}\n", String::from_utf8_lossy(read)),
                        )
                        .unwrap();
                        let sample = sample::build(&panel, identity.clone(), &[path]).unwrap();
                        for (tokens, _, _) in locations.keys() {
                            assert_eq!(
                                sample.counts.count(tokens).unwrap(),
                                coordinate_counts.get(tokens).copied().unwrap_or(0)
                            );
                        }
                        singleton_checks += 1;
                    }
                }
            }
            for ((tokens, start, end), expected) in expected {
                let actual = profile::total(
                    &panel,
                    &views,
                    0,
                    start,
                    end,
                    source.len() as u64,
                    length,
                    &tokens.try_into().unwrap(),
                )
                .unwrap()
                .0;
                assert_eq!(actual, expected, "span {start}-{end}, length {length}");
                checked += 1;
            }
        }
    }
    assert!(checked > 100 && crops > 1000 && singleton_checks > 20);
    assert!(multicopy_reads > 0);
    assert!(exact_anchor_windows > 0 && anchors_at_read_end > exact_anchor_windows);
    eprintln!("exact-K anchor windows={exact_anchor_windows}, anchors ending at read end={anchors_at_read_end}");
    eprintln!("max native q={maximum_physical_multiplicity}, multi-copy reads={multicopy_reads}");
    eprintln!("native profile totals={checked}, raw crop comparisons={crops}, public singleton indexes={singleton_checks}");
}

#[test]
fn fixed_ownership_normalizes_duplicate_descriptions_and_rejects_conflicts() {
    let identity = PanelIdentity {
        checksum_algorithm: "fixture".into(),
        sidecars: vec![],
    };
    let mut c = catalog::Catalog {
        version: 1,
        panel: identity,
        catalog_accepted: false,
        validation: String::new(),
        sources: vec![catalog::Source {
            id: 0,
            path: "A#0#a".into(),
            length: 100,
        }],
        groups: vec![],
        occurrences: vec![],
        features: vec![],
        links: vec![],
    };
    for (id, (start, end, group)) in [(0, 60, 0), (0, 60, 0), (20, 50, 0), (60, 100, 1)]
        .into_iter()
        .enumerate()
    {
        c.occurrences.push(catalog::Occurrence {
            id,
            group,
            source: 0,
            interval: catalog::Interval {
                path: "A#0#a".into(),
                start,
                end,
                strand: None,
            },
            fully_contained_anchors: 0,
            anchor_status: String::new(),
            feature_multiplicities: BTreeMap::new(),
        });
    }
    let cores = input::ownership(&c).unwrap();
    assert_eq!(cores[0].len(), 2);
    assert_eq!((cores[0][0].start, cores[0][0].end), (0, 60));
    c.occurrences[2].group = 1;
    assert!(input::ownership(&c).is_err());
    c.occurrences[2].group = 0;
    c.occurrences[3].interval.start = 61;
    assert!(input::ownership(&c).is_err());
    let mut maximum = u64::MAX;
    assert!(profile::add(&mut maximum, 1).is_err());
}

#[test]
fn overlapping_survivors_and_two_copy_attribution_are_integer_not_detection() {
    // The production collector's retained overlapping records, with the same
    // physical middle pair twice, and a different identical physical copy once.
    let records = [
        vec![(1, 0), (2, 5), (3, 10)],
        vec![(2, 5), (3, 10), (4, 15)],
        vec![(2, 40), (3, 45)],
    ];
    let tokens = canonical(&encode_walk(&[(2, 0), (3, 5)]).unwrap());
    let mut counts = BTreeMap::new();
    for record in records {
        for pair in record.windows(2) {
            if canonical(&encode_walk(pair).unwrap()) == tokens {
                *counts.entry((pair[0].1, pair[1].1 + 3)).or_insert(0) += 1;
            }
        }
    }
    assert_eq!(counts, BTreeMap::from([((5, 13), 2), ((40, 48), 1)]));
    assert_eq!(counts.values().sum::<u64>(), 3);
}

#[test]
fn native_overlapping_mem_profiles_preserve_double_counts() {
    let source = dna(700, 113);
    // No complete graph walk spans this synthetic observer context. Its two
    // surviving, overlapping MEMs share physical pairs but neither contains the
    // other. This is a native collector oracle, not an added source candidate.
    let panel = SyngIndex::build(
        SyncmerParams::default(),
        vec![
            ("L#0#left".into(), source[..500].to_vec()),
            ("R#0#right".into(), source[200..].to_vec()),
        ]
        .into_iter(),
    );
    let views = profile::raw_views(&panel, &source).unwrap();
    let records = collect_tagged_read(&panel, &source).unwrap();
    let k = panel.syncmer_length_bp() as u64;
    let mut expected = BTreeMap::new();
    for record in records {
        for pair in record.windows(2) {
            *expected
                .entry((
                    canonical(&encode_walk(pair).unwrap()),
                    pair[0].1,
                    pair[1].1 + k,
                ))
                .or_insert(0) += 1;
        }
    }
    assert!(
        expected.values().any(|&q| q >= 2),
        "fixture must exercise native overlapping records"
    );
    let temp = tempfile::tempdir().unwrap();
    let fasta = temp.path().join("overlap-singleton.fa");
    fs::write(
        &fasta,
        format!(">transient\n{}\n", String::from_utf8_lossy(&source)),
    )
    .unwrap();
    let sample = sample::build(
        &panel,
        PanelIdentity {
            checksum_algorithm: "fixture".into(),
            sidecars: vec![],
        },
        &[fasta],
    )
    .unwrap();
    let mut totals = BTreeMap::new();
    for ((tokens, _, _), &q) in &expected {
        *totals.entry(tokens.clone()).or_insert(0u64) += q;
    }
    for (tokens, q) in totals {
        assert_eq!(sample.counts.count(&tokens).unwrap(), q);
    }
    for ((tokens, start, end), q) in expected {
        assert_eq!(
            profile::total(
                &panel,
                &views,
                0,
                start,
                end,
                source.len() as u64,
                source.len() as u64,
                &tokens.try_into().unwrap()
            )
            .unwrap()
            .0,
            q
        );
    }
}

#[test]
fn native_partition_endpoints_recover_long_zero_exposure_across_empty_core() {
    let mut source = dna(160, 71);
    source.extend(vec![b'N'; 800]);
    source.extend(dna(160, 73));
    let sequences = vec![
        ("A#0#a".into(), source.clone()),
        ("B#0#b".into(), source.clone()),
    ];
    let panel = SyngIndex::build(SyncmerParams::default(), sequences.into_iter());
    let identity = PanelIdentity {
        checksum_algorithm: "fixture".into(),
        sidecars: vec![],
    };
    let interval = |path: &str, start, end| catalog::Interval {
        path: path.into(),
        start,
        end,
        strand: None,
    };
    let c = catalog::build(
        &panel,
        identity.clone(),
        catalog::CatalogInput {
            version: 1,
            groups: vec![
                catalog::GroupInput {
                    id: "whole".into(),
                    scaffold: None,
                    occurrences: vec![interval("A#0#a", 0, 1120), interval("B#0#b", 0, 400)],
                },
                catalog::GroupInput {
                    id: "empty-core".into(),
                    scaffold: None,
                    occurrences: vec![interval("B#0#b", 400, 800)],
                },
                catalog::GroupInput {
                    id: "end".into(),
                    scaffold: None,
                    occurrences: vec![interval("B#0#b", 800, 1120)],
                },
            ],
        },
    )
    .unwrap();
    let temp = tempfile::tempdir().unwrap();
    let legacy = temp.path().join("catalog.json");
    super::super::save_catalog(&legacy, &c).unwrap();
    let fasta = temp.path().join("source.fa");
    fs::write(
        &fasta,
        format!(
            ">A#0#a\n{}\n>B#0#b\n{}\n",
            String::from_utf8_lossy(&source),
            String::from_utf8_lossy(&source)
        ),
    )
    .unwrap();
    let out = temp.path().join("observations");
    super::super::with_output_model(&out, MODEL, || {
        build(
            &panel,
            identity.clone(),
            &legacy,
            &[fasta.to_str().unwrap().into()],
            vec![150],
            &[],
            &out,
        )
    })
    .unwrap();
    let mut reader = BufReader::new(File::open(out.join("incidences.jsonl")).unwrap());
    let mut long = Vec::new();
    while let Some(h) = next::<Incidence>(&mut reader).unwrap() {
        if h.end - h.start > 500 {
            long.push(h);
        }
    }
    assert!(long.iter().any(|h| h.source == 0 && h.group == Some(0)));
    assert!(long.iter().any(|h| h.source == 1 && h.group.is_none()));
    assert!(long.iter().all(|h| h.contributions == vec![[0, 0]]));
    let checkpoint: Value =
        super::super::read_json(&out.join("shards/partition-00000001.json")).unwrap();
    assert_eq!(checkpoint["completed"], true);
    assert_eq!(checkpoint["shards"], json!([]));
    let metadata: Metadata = super::super::read_json(&out.join("metadata.json")).unwrap();
    assert_eq!(metadata.stats["sources_scanned"], 2);
    assert_eq!(metadata.stats["partitions_scanned"], 3);
    assert_eq!(metadata.stats["core_chunks"], 4);
    for entry in fs::read_dir(out.join("shards")).unwrap() {
        let entry = entry.unwrap();
        if entry.path().extension().is_some_and(|ext| ext == "jsonl") {
            assert!(entry.metadata().unwrap().len() > 0);
        }
    }
}

#[test]
fn new_model_background_eligibility_and_reset_preserve_legacy_classification() {
    let identity = PanelIdentity {
        checksum_algorithm: "fixture".into(),
        sidecars: vec![],
    };
    let mut c = catalog::Catalog {
        version: 1,
        panel: identity.clone(),
        catalog_accepted: false,
        validation: String::new(),
        sources: vec![
            catalog::Source {
                id: 0,
                path: "A#0#a".into(),
                length: 20,
            },
            catalog::Source {
                id: 1,
                path: "B#0#b".into(),
                length: 20,
            },
        ],
        groups: vec![],
        occurrences: vec![],
        features: vec![],
        links: vec![],
    };
    for g in 0..2 {
        let mut ids = Vec::new();
        for source in 0..2 {
            let id = c.occurrences.len();
            ids.push(id);
            c.occurrences.push(catalog::Occurrence {
                id,
                group: g,
                source,
                interval: catalog::Interval {
                    path: c.sources[source].path.clone(),
                    start: g as u64 * 10,
                    end: g as u64 * 10 + 10,
                    strand: Some("+".into()),
                },
                fully_contained_anchors: 0,
                anchor_status: String::new(),
                feature_multiplicities: BTreeMap::new(),
            });
        }
        c.groups.push(catalog::Group {
            id: format!("g{g}"),
            scaffold: None,
            occurrences: ids,
        });
    }
    let sample = sample::SampleIndex {
        version: 1,
        panel: identity,
        count_policy: COUNT_POLICY.into(),
        stats: sample::SampleStats {
            read_lengths: BTreeMap::from([(10, 1)]),
            ..Default::default()
        },
        counts: crate::sample_mem_bwt::WeightedBwt::build(&BTreeMap::new()).unwrap(),
    };
    let parameters = genotype::Parameters {
        haploid_depth: 10.0,
        background: 0.1,
        max_mean_deviance: 10.0,
    };
    let mut calls = genotype::call(&c, &sample, parameters.clone(), true, 1).unwrap();
    calls.model = MODEL.into();
    calls.observations = Some(Provenance {
        metadata_checksum: "1234567890123456".into(),
        compiler_identity: compiler_identity(),
        full_feature_scope: true,
        feature_groups: vec![],
        legacy_input: json!({}),
        source_files: vec![],
        orientation_references: vec![0, 2],
        derived_from_calls: None,
        orientation_support: vec![],
    });
    for row in &mut calls.calls {
        row.retained_features = 1;
        row.positive_features = 1;
        row.observed_total = 1.0;
        row.background_only_score = 5.0;
        for b in &mut row.bundles {
            b.score = 5.0;
            b.positive_features = 1;
        }
    }
    // BG equality/preference means no eligible source despite group evidence.
    assert_eq!(
        classify(&calls.calls[0], &parameters).0,
        "ambiguous-background-equivalent"
    );
    calls.calls[0].bundles[0].score = 5.0 - genotype::TIE_EPSILON / 2.0;
    assert!(!eligible(&calls.calls[0], &calls.calls[0].bundles[0]));
    calls.calls[0].bundles[0].score = 6.0;
    calls.calls[0].bundles[1].score = 7.0;
    // A score beating BG still needs its OWN positive support.
    calls.calls[1].bundles[0].score = 0.0;
    calls.calls[1].bundles[0].positive_features = 0;
    calls.calls[1].bundles[1].score = 1.0;
    for row in &mut calls.calls {
        (row.status, row.best_bundles, row.score_gap) = classify(row, &parameters);
    }
    assert_eq!(calls.calls[1].best_bundles, vec![1]);
    assert_eq!(
        calls.calls[1].score_gap,
        Some(1.0),
        "gap must include the ineligible source"
    );
    let axis = || threading::Axis {
        version: 1,
        coordinate_system: "fixture".into(),
        intervals: (0..2)
            .map(|g| threading::AxisInterval {
                component: "chr".into(),
                start: g as u64 * 10,
                end: g as u64 * 10 + 10,
                group: format!("g{g}"),
                reference_occurrence: g * 2,
                reference_strand: "+".into(),
                orientations: BTreeMap::new(),
            })
            .collect(),
    };
    let threads = threading::thread(&c, &calls, axis(), 10000.0).unwrap();
    assert!(threads.intervals[0].states.is_empty());
    assert_eq!(threads.segments[0].first_interval, 1);
    assert_eq!(threads.intervals[1].states.len(), 1);
    assert_eq!(threads.intervals[1].states[0].identity, "B#0");
    calls.calls[1].bundles[0].positive_features = 1;
    calls.calls[1].bundles[0].score = 2.0;
    let row = &mut calls.calls[1];
    (row.status, row.best_bundles, row.score_gap) = classify(row, &parameters);
    assert_eq!(row.best_bundles, vec![1]);
    let alternatives = threading::thread(&c, &calls, axis(), 10000.0).unwrap();
    assert_eq!(
        alternatives.intervals[1].states.len(),
        2,
        "eligible non-minimum must survive"
    );
    let new_provenance = calls.observations.clone();
    calls.model = genotype::MODEL.into();
    calls.observations = None;
    calls.calls[0].status = "informative".into();
    calls.calls[1].status = "informative".into();
    let legacy = threading::thread(&c, &calls, axis(), 10000.0).unwrap();
    assert_eq!(legacy.intervals[0].states.len(), 2);
    assert_eq!(legacy.intervals[1].states.len(), 2);
    calls.model = MODEL.into();
    calls.observations = new_provenance;
    calls.calls[0].bundles[0].score = 5.0;
    calls.calls[0].bundles[1].score = 6.0;
    for row in &mut calls.calls {
        (row.status, row.best_bundles, row.score_gap) = classify(row, &parameters);
    }
    assert_eq!(calls.calls[0].status, "ambiguous-background-equivalent");
    assert!(calls.calls[0].best_bundles.is_empty());
    assert_eq!(calls.calls[0].score_gap, Some(1.0));
    let blocked = threading::thread(&c, &calls, axis(), 10000.0).unwrap();
    assert_eq!(
        blocked.intervals[0].status,
        "unresolved-ambiguous-background-equivalent"
    );
    assert!(
        blocked.intervals[0].states.is_empty()
            && blocked.intervals[0].optimal_states.is_empty()
            && blocked.intervals[0].segment.is_none()
    );
    assert_eq!(blocked.segments[0].first_interval, 1);
}

#[test]
fn occurrence_weighted_consumer_sums_copies_inside_one_rate_and_preserves_zero_factors() {
    let identity = PanelIdentity {
        checksum_algorithm: "fixture".into(),
        sidecars: vec![],
    };
    let sources = vec![
        catalog::Source {
            id: 0,
            path: "A#0#a".into(),
            length: 1000,
        },
        catalog::Source {
            id: 1,
            path: "B#0#b".into(),
            length: 1000,
        },
    ];
    let occurrences = [0, 1, 1]
        .into_iter()
        .enumerate()
        .map(|(id, source)| catalog::Occurrence {
            id,
            group: 0,
            source,
            interval: catalog::Interval {
                path: sources[source].path.clone(),
                start: 0,
                end: 1000,
                strand: None,
            },
            fully_contained_anchors: 0,
            anchor_status: String::new(),
            feature_multiplicities: BTreeMap::new(),
        })
        .collect();
    let c = catalog::Catalog {
        version: 1,
        panel: identity.clone(),
        catalog_accepted: false,
        validation: String::new(),
        sources,
        occurrences,
        groups: vec![catalog::Group {
            id: "g".into(),
            scaffold: None,
            occurrences: vec![0, 1, 2],
        }],
        features: vec![],
        links: vec![],
    };
    let tokens = canonical(&encode_walk(&[(1, 0), (2, 7)]).unwrap());
    let zero = canonical(&encode_walk(&[(3, 0), (4, 500)]).unwrap());
    let sample = sample::SampleIndex {
        version: 1,
        panel: identity,
        count_policy: COUNT_POLICY.into(),
        stats: sample::SampleStats {
            read_lengths: BTreeMap::from([(100, 2), (200, 1)]),
            ..Default::default()
        },
        counts: crate::sample_mem_bwt::WeightedBwt::build(&BTreeMap::from([(tokens.clone(), 9)]))
            .unwrap(),
    };
    let temp = tempfile::tempdir().unwrap();
    write_lines(
        &temp.path().join("definitions.jsonl"),
        [
            input::Definition {
                id: 0,
                tokens: tokens.try_into().unwrap(),
                original_owner: Some(0),
            },
            input::Definition {
                id: 1,
                tokens: zero.try_into().unwrap(),
                original_owner: Some(0),
            },
        ],
    )
    .unwrap();
    write_lines(
        &temp.path().join("incidences.jsonl"),
        [
            Incidence {
                feature: 0,
                source: 0,
                start: 100,
                end: 170,
                signed_mask: 1,
                group: Some(0),
                context_nonlocal: false,
                source_terminal_context: false,
                contributions: vec![[80, 80], [200, 200]],
            },
            Incidence {
                feature: 0,
                source: 1,
                start: 100,
                end: 170,
                signed_mask: 1,
                group: Some(0),
                context_nonlocal: false,
                source_terminal_context: false,
                contributions: vec![[30, 50], [100, 100]],
            },
            Incidence {
                feature: 0,
                source: 1,
                start: 300,
                end: 370,
                signed_mask: 1,
                group: Some(0),
                context_nonlocal: false,
                source_terminal_context: false,
                contributions: vec![[40, 40], [100, 100]],
            },
            Incidence {
                feature: 1,
                source: 0,
                start: 0,
                end: 563,
                signed_mask: 1,
                group: Some(0),
                context_nonlocal: false,
                source_terminal_context: false,
                contributions: vec![[0, 0], [0, 0]],
            },
        ],
    )
    .unwrap();
    let metadata = Metadata {
        version: 1,
        model: MODEL.into(),
        count_policy: COUNT_POLICY.into(),
        compiler_identity: compiler_identity(),
        capability: String::new(),
        catalog: c,
        legacy_input: json!({}),
        sources: vec![],
        syncmer_length_bp: 63,
        read_lengths: vec![100, 200],
        feature_groups: vec![],
        full_feature_scope: true,
        feature_count: 2,
        files: BTreeMap::new(),
        stats: json!({}),
    };
    let calls = call(
        &metadata,
        "1234567890123456".into(),
        temp.path(),
        &sample,
        "1234567890123456".into(),
        genotype::Parameters {
            haploid_depth: 10.0,
            background: 0.1,
            max_mean_deviance: 10.0,
        },
        temp.path(),
    )
    .unwrap();
    assert_eq!(calls.global_factors_used, 2);
    assert_eq!(calls.calls[0].status, "tied");
    for b in &calls.calls[0].bundles {
        assert!((b.expected_total - 9.2).abs() < 1e-10);
    }
    assert_eq!(calls.calls[0].bundles[0].physical_feature_locations, 2); // one positive + one long zero
    assert_eq!(calls.calls[0].bundles[1].physical_feature_locations, 2); // two real copies, not four BED descriptions
    assert_eq!(calls.calls[0].bundles[0].positive_features, 1);
    assert_eq!(calls.calls[0].bundles[1].positive_features, 1);
}

#[test]
fn deterministic_merge_unions_descriptions_without_adding_profiles_and_rejects_disagreement() {
    let temp = tempfile::tempdir().unwrap();
    let left = temp.path().join("left.jsonl");
    let right = temp.path().join("right.jsonl");
    let make = |feature, source, mask| Incidence {
        feature,
        source,
        start: 10,
        end: 80,
        signed_mask: mask,
        group: Some(0),
        context_nonlocal: false,
        source_terminal_context: false,
        contributions: vec![[7, 7]],
    };
    write_lines(&left, [make(0, 0, 1), make(1, 0, 1)]).unwrap();
    write_lines(&right, [make(0, 0, 2), make(0, 1, 1)]).unwrap();
    let a = temp.path().join("a.jsonl");
    let b = temp.path().join("b.jsonl");
    merge(&[left.clone(), right.clone()], &a).unwrap();
    merge(&[right.clone(), left.clone()], &b).unwrap();
    assert_eq!(fs::read(&a).unwrap(), fs::read(&b).unwrap());
    let mut reader = BufReader::new(File::open(&a).unwrap());
    let first = next::<Incidence>(&mut reader).unwrap().unwrap();
    assert_eq!(first.signed_mask, 3);
    assert_eq!(first.contributions, vec![[7, 7]]);
    assert_eq!(next::<Incidence>(&mut reader).unwrap().unwrap().source, 1);
    assert_eq!(next::<Incidence>(&mut reader).unwrap().unwrap().feature, 1);
    assert!(next::<Incidence>(&mut reader).unwrap().is_none());
    let mut bad = make(0, 0, 2);
    bad.contributions[0][0] = 8;
    write_lines(&right, [bad]).unwrap();
    assert!(merge(&[left, right], &temp.path().join("bad.jsonl")).is_err());
    assert!(next::<Incidence>(&mut std::io::Cursor::new(vec![b' '; 65537])).is_err());

    let empty = temp.path().join("zero-shards.jsonl");
    merge(&[], &empty).unwrap();
    assert_eq!(fs::metadata(&empty).unwrap().len(), 0);
    assert!(!empty.with_extension("incomplete").exists());
    let failed = temp.path().join("missing-input-output.jsonl");
    assert!(merge(&[temp.path().join("missing.jsonl")], &failed).is_err());
    assert!(!failed.exists());
}

#[test]
fn terminal_and_longer_than_source_contexts_are_nonlocal_but_profiles_stay_clipped() {
    let source = dna(200, 17);
    let panel = SyngIndex::build(
        SyncmerParams::default(),
        vec![("A#0#a".into(), source.clone())].into_iter(),
    );
    let identity = PanelIdentity {
        checksum_algorithm: "fixture".into(),
        sidecars: vec![],
    };
    let c = catalog::build(
        &panel,
        identity,
        catalog::CatalogInput {
            version: 1,
            groups: vec![catalog::GroupInput {
                id: "whole".into(),
                scaffold: None,
                occurrences: vec![catalog::Interval {
                    path: "A#0#a".into(),
                    start: 0,
                    end: 200,
                    strand: None,
                }],
            }],
        },
    )
    .unwrap();
    let temp = tempfile::tempdir().unwrap();
    let fasta = temp.path().join("source.fa");
    fs::write(
        &fasta,
        format!(">A#0#a\n{}\n", String::from_utf8_lossy(&source)),
    )
    .unwrap();
    let sources = input::Sources::open(&[fasta.to_str().unwrap().into()], &c).unwrap();
    let cores = input::ownership(&c).unwrap();
    let views = profile::raw_views(&panel, &source).unwrap();
    let k = panel.syncmer_length_bp() as u64;
    let pair = views[0].windows(2).find(|p| p[1].1 + k < 150).unwrap();
    let tokens: [u64; 3] = canonical(&encode_walk(pair).unwrap()).try_into().unwrap();
    let definition = input::Definition {
        id: 0,
        tokens,
        original_owner: Some(0),
    };
    let make = || Incidence {
        feature: 0,
        source: 0,
        start: pair[0].1,
        end: pair[1].1 + k,
        signed_mask: 1,
        group: None,
        context_nonlocal: false,
        source_terminal_context: false,
        contributions: vec![],
    };
    let mut terminal = make();
    profile::compile(
        &panel,
        &sources,
        &c,
        &cores,
        &[150],
        &definition,
        &mut terminal,
    )
    .unwrap();
    assert!(terminal.context_nonlocal && terminal.source_terminal_context);
    assert_eq!(terminal.group, Some(0));
    assert_eq!(
        terminal.contributions[0][0],
        profile::total(
            &panel,
            &views,
            0,
            terminal.start,
            terminal.end,
            200,
            150,
            &tokens
        )
        .unwrap()
        .0
    );
    let mut too_long = make();
    profile::compile(
        &panel,
        &sources,
        &c,
        &cores,
        &[250],
        &definition,
        &mut too_long,
    )
    .unwrap();
    assert!(too_long.context_nonlocal && too_long.source_terminal_context);
    assert_eq!(too_long.contributions, vec![[0, 0]]);
    let mut no_spanning_context = make();
    profile::compile(
        &panel,
        &sources,
        &c,
        &cores,
        &[30],
        &definition,
        &mut no_spanning_context,
    )
    .unwrap();
    assert!(!no_spanning_context.context_nonlocal && !no_spanning_context.source_terminal_context);
    assert_eq!(no_spanning_context.contributions, vec![[0, 0]]);
}

#[test]
fn shared_core_halos_match_bounded_native_fallback_without_per_hit_extraction() {
    let mut forward = dna(3400, 75);
    forward[650..750].fill(b'N');
    let repeated = forward[200..450].to_vec();
    forward[900..1150].copy_from_slice(&repeated);
    let reverse = crate::graph::reverse_complement(&forward);
    let sequences = vec![
        ("A#0#a".into(), forward.clone()),
        ("R#0#r".into(), reverse.clone()),
    ];
    let panel = SyngIndex::build(SyncmerParams::default(), sequences.into_iter());
    let identity = PanelIdentity {
        checksum_algorithm: "fixture".into(),
        sidecars: vec![],
    };
    let interval = |path: &str, start, end| catalog::Interval {
        path: path.into(),
        start,
        end,
        strand: None,
    };
    let catalog = catalog::build(
        &panel,
        identity,
        catalog::CatalogInput {
            version: 1,
            groups: vec![
                catalog::GroupInput {
                    id: "left".into(),
                    scaffold: None,
                    occurrences: vec![interval("A#0#a", 0, 1700), interval("R#0#r", 1700, 3400)],
                },
                catalog::GroupInput {
                    id: "right".into(),
                    scaffold: None,
                    occurrences: vec![interval("A#0#a", 1700, 3400), interval("R#0#r", 0, 1700)],
                },
            ],
        },
    )
    .unwrap();
    let definitions: Vec<_> = catalog
        .features
        .iter()
        .enumerate()
        .map(|(id, f)| input::Definition {
            id,
            tokens: f.tokens.clone().try_into().unwrap(),
            original_owner: f.owning_group,
        })
        .collect();
    let mut lookup: Vec<_> = (0..definitions.len()).collect();
    lookup.sort_by_key(|&f| definitions[f].tokens);
    let cores = input::ownership(&catalog).unwrap();
    let temp = tempfile::tempdir().unwrap();
    let fasta = temp.path().join("source.fa");
    fs::write(
        &fasta,
        format!(
            ">A#0#a\n{}\n>R#0#r\n{}\n",
            String::from_utf8_lossy(&forward),
            String::from_utf8_lossy(&reverse)
        ),
    )
    .unwrap();
    let sources = input::Sources::open(&[fasta.to_str().unwrap().into()], &catalog).unwrap();
    let lengths = [62, 63, 64, 150, 220];
    let halo = 220;
    let k = panel.syncmer_length_bp() as u64;
    let mut stats = profile::CompilationStats::default();
    let mut profiles = 0;
    let mut extractions = 0;
    let mut maximum_core_anchors = 0;
    let mut wrong_context_checked = false;
    for owner in cores.iter().flatten() {
        let mut start = owner.start;
        while start < owner.end {
            let end = (start + 400).min(owner.end);
            let crop_start = start.saturating_sub(halo);
            let crop_end = (end + halo).min(catalog.sources[owner.source].length);
            let dna = sources
                .fetch(&catalog, owner.source, crop_start, crop_end)
                .unwrap();
            let context = profile::SourceContext {
                source: owner.source,
                crop_start,
                crop_end,
                views: profile::raw_views(&panel, &dna).unwrap(),
            };
            extractions += 1;
            maximum_core_anchors =
                maximum_core_anchors.max(context.views[0].len() + context.views[1].len());
            let mut hits = BTreeMap::new();
            for view in &context.views {
                let owned: Vec<_> = view
                    .iter()
                    .map(|&(n, p)| (n, p + crop_start))
                    .filter(|&(_, p)| p >= start && p < end)
                    .collect();
                for pair in owned.windows(2) {
                    if let Some((f, h)) =
                        candidate(pair, owner.source, k, &definitions, &lookup).unwrap()
                    {
                        hits.insert(h.key(), (f, h));
                    }
                }
            }
            for (_, (f, mut shared)) in hits {
                let mut fallback = shared.clone();
                let expected = profile::compile(
                    &panel,
                    &sources,
                    &catalog,
                    &cores,
                    &lengths,
                    &definitions[f],
                    &mut fallback,
                )
                .unwrap();
                let actual = profile::compile_with_context(
                    &panel,
                    &sources,
                    &catalog,
                    &cores,
                    &lengths,
                    &definitions[f],
                    &mut shared,
                    Some(&context),
                    &mut stats,
                )
                .unwrap();
                assert_eq!(actual, expected);
                assert_eq!(
                    serde_json::to_value(&shared).unwrap(),
                    serde_json::to_value(&fallback).unwrap()
                );
                // total() itself must also restrict long input traces before its
                // event loop; supplying a full core gives the same exact totals.
                for (i, &length) in lengths.iter().enumerate() {
                    assert_eq!(
                        profile::total(
                            &panel,
                            &context.views,
                            crop_start,
                            shared.start,
                            shared.end,
                            catalog.sources[owner.source].length,
                            length,
                            &definitions[f].tokens
                        )
                        .unwrap()
                        .0,
                        shared.contributions[i][0]
                    );
                }
                if !wrong_context_checked && shared.end - shared.start < 150 {
                    let missing = profile::SourceContext {
                        source: shared.source,
                        crop_start: shared.start,
                        crop_end: shared.end,
                        views: [vec![], vec![]],
                    };
                    let mut hit = shared.clone();
                    hit.contributions.clear();
                    assert!(profile::compile_with_context(
                        &panel,
                        &sources,
                        &catalog,
                        &cores,
                        &lengths,
                        &definitions[f],
                        &mut hit,
                        Some(&missing),
                        &mut profile::CompilationStats::default()
                    )
                    .is_err());
                    wrong_context_checked = true;
                }
                profiles += 1;
            }
            start = end;
        }
    }
    assert!(profiles > 100 && stats.reused_context_profiles > extractions);
    assert_eq!(stats.fallback_context_fetches, 0);
    assert!(stats.max_profile_context_bp <= 2 * halo);
    assert!(stats.max_profile_context_anchors < maximum_core_anchors);
    assert!(wrong_context_checked);
    eprintln!("shared core oracle: {profiles} physical profiles, {extractions} core extractions, {stats:?}");
}

#[test]
fn positive_endpoint_fallback_matches_shared_core_profile_across_ownership_boundary() {
    let source = dna(1000, 79);
    let panel = SyngIndex::build(
        SyncmerParams::default(),
        vec![
            ("A#0#a".into(), source.clone()),
            ("B#0#b".into(), source.clone()),
        ]
        .into_iter(),
    );
    let k = panel.syncmer_length_bp() as u64;
    let walk = panel.walk_path_range(0, 0, 1000).unwrap();
    let pair = walk
        .windows(2)
        .find(|p| p[0].1 > 300 && p[1].1 < 700 && p[1].1 + k - p[0].1 < 150)
        .unwrap();
    let (start, end) = (pair[0].1, pair[1].1 + k);
    let border = (pair[0].1 + pair[1].1 + 1) / 2;
    let tokens = canonical(&encode_walk(pair).unwrap());
    let identity = PanelIdentity {
        checksum_algorithm: "fixture".into(),
        sidecars: vec![],
    };
    let interval = |path: &str, start, end| catalog::Interval {
        path: path.into(),
        start,
        end,
        strand: None,
    };
    let catalog = catalog::build(
        &panel,
        identity.clone(),
        catalog::CatalogInput {
            version: 1,
            groups: vec![
                catalog::GroupInput {
                    id: "whole".into(),
                    scaffold: None,
                    occurrences: vec![interval("A#0#a", 0, 1000)],
                },
                catalog::GroupInput {
                    id: "left".into(),
                    scaffold: None,
                    occurrences: vec![interval("B#0#b", 0, border)],
                },
                catalog::GroupInput {
                    id: "right".into(),
                    scaffold: None,
                    occurrences: vec![interval("B#0#b", border, 1000)],
                },
            ],
        },
    )
    .unwrap();
    let feature = catalog
        .features
        .iter()
        .position(|f| f.tokens == tokens)
        .unwrap();
    let temp = tempfile::tempdir().unwrap();
    let legacy = temp.path().join("catalog.json");
    super::super::save_catalog(&legacy, &catalog).unwrap();
    let fasta = temp.path().join("source.fa");
    fs::write(
        &fasta,
        format!(
            ">A#0#a\n{}\n>B#0#b\n{}\n",
            String::from_utf8_lossy(&source),
            String::from_utf8_lossy(&source)
        ),
    )
    .unwrap();
    let out = temp.path().join("observations");
    super::super::with_output_model(&out, MODEL, || {
        build(
            &panel,
            identity,
            &legacy,
            &[fasta.to_str().unwrap().into()],
            vec![150],
            &[],
            &out,
        )
    })
    .unwrap();
    let metadata: Metadata = super::super::read_json(&out.join("metadata.json")).unwrap();
    let stats = &metadata.stats["profile_execution"];
    assert!(stats["reused_context_profiles"].as_u64().unwrap() > 0);
    assert!(stats["fallback_context_fetches"].as_u64().unwrap() > 0);
    assert!(
        stats["fallback_context_bp_fetched"].as_u64().unwrap()
            <= 300 * stats["fallback_context_fetches"].as_u64().unwrap()
    );
    let mut reader = BufReader::new(File::open(out.join("incidences.jsonl")).unwrap());
    let mut matched = BTreeMap::new();
    while let Some(h) = next::<Incidence>(&mut reader).unwrap() {
        if h.feature == feature && h.start == start && h.end == end {
            matched.insert(h.source, h);
        }
    }
    assert_eq!(matched.len(), 2);
    assert!(matched[&0].contributions[0][0] > 0);
    assert_eq!(matched[&0].contributions, matched[&1].contributions);
    assert_eq!(matched[&0].group, Some(0));
    assert_eq!(matched[&1].group, None);
    assert!(!matched[&0].context_nonlocal && matched[&1].context_nonlocal);
    eprintln!("positive endpoint fallback: {}", metadata.stats);
}
