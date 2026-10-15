//! Axis-specific direct orientation evidence. The reusable incidence stream is
//! the registry; neither calling nor rethreading constructs an all-pairs matrix.
use super::*;

/// Max-end interval tree over original descriptions for one (group, source).
/// Candidate queries prune both starts beyond the span and subtrees ending before
/// it. Duplicate descriptions remain distinct IDs without scanning a whole path.
struct IntervalIndex {
    ids: Vec<usize>,
    size: usize,
    max_end: Vec<u64>,
}
impl IntervalIndex {
    fn new(mut ids: Vec<usize>, catalog: &catalog::Catalog) -> Self {
        ids.sort_unstable_by_key(|&id| (catalog.occurrences[id].interval.start, id));
        let size = ids.len().next_power_of_two();
        let mut max_end = vec![0; 2 * size];
        for (i, &id) in ids.iter().enumerate() {
            max_end[size + i] = catalog.occurrences[id].interval.end;
        }
        for i in (1..size).rev() {
            max_end[i] = max_end[2 * i].max(max_end[2 * i + 1]);
        }
        Self { ids, size, max_end }
    }
    fn containing(
        &self,
        catalog: &catalog::Catalog,
        start: u64,
        end: u64,
        mut visit: impl FnMut(usize),
    ) {
        let limit = self
            .ids
            .partition_point(|&id| catalog.occurrences[id].interval.start <= start);
        let mut stack = vec![(1, 0, self.size)];
        while let Some((node, lo, hi)) = stack.pop() {
            if lo >= limit || self.max_end[node] < end {
                continue;
            }
            if hi - lo == 1 {
                visit(self.ids[lo]);
            } else {
                let mid = lo + (hi - lo) / 2;
                stack.push((node * 2 + 1, mid, hi));
                stack.push((node * 2, lo, mid));
            }
        }
    }
}

pub fn axis_references(
    catalog: &catalog::Catalog,
    axis: &threading::Axis,
) -> io::Result<BTreeMap<usize, usize>> {
    let groups = threading::validate_axis(catalog, axis)?;
    let mut frequencies = BTreeMap::new();
    for &g in &groups {
        *frequencies.entry(g).or_insert(0) += 1;
    }
    Ok(axis
        .intervals
        .iter()
        .zip(groups)
        .filter(|(_, g)| frequencies[g] == 1)
        .map(|(a, g)| (g, a.reference_occurrence))
        .collect())
}

fn accumulate(
    signs: &BTreeMap<usize, BTreeMap<usize, u8>>,
    references: &BTreeMap<usize, usize>,
    support: &mut [u8],
) {
    for (&g, members) in signs {
        let Some(&x) = members.get(&references[&g]) else {
            continue;
        };
        for (&id, &y) in members {
            let same = x & y != 0;
            let reverse = (x & 1 != 0 && y & 2 != 0) || (x & 2 != 0 && y & 1 != 0);
            support[id] |= u8::from(same) | (u8::from(reverse) << 1);
        }
    }
}

/// At most one assessed reference per non-repeated axis group. Time per feature
/// is linear in that group's matching descriptions, not their Cartesian product;
/// persistent masks occupy one byte per original occurrence. Zero-support
/// references are explicitly retained as assessed, not silently uncomputed.
pub fn prepare_orientations(
    metadata: &Metadata,
    directory: &Path,
    calls: &mut genotype::Genotypes,
    axis: &threading::Axis,
) -> io::Result<()> {
    require(
        metadata.full_feature_scope && calls.model == MODEL,
        "partial/incompatible calls cannot derive axis orientations",
    )?;
    super::consumer::validate_frozen_calls(metadata, calls)?;
    let p = calls
        .observations
        .as_mut()
        .ok_or_else(|| invalid("missing observation provenance"))?;
    require(
        p.compiler_identity == metadata.compiler_identity && p.full_feature_scope,
        "incompatible orientation compiler/scope",
    )?;
    let catalog = &metadata.catalog;
    let references = axis_references(catalog, axis)?;
    let mut grouped: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for o in &catalog.occurrences {
        if references.contains_key(&o.group) {
            grouped.entry((o.group, o.source)).or_default().push(o.id);
        }
    }
    let index: BTreeMap<_, _> = grouped
        .into_iter()
        .map(|(key, ids)| (key, IntervalIndex::new(ids, catalog)))
        .collect();
    let mut support = vec![0; catalog.occurrences.len()];
    let mut signs: BTreeMap<usize, BTreeMap<usize, u8>> = BTreeMap::new();
    let mut feature = None;
    let mut previous = None;
    let mut reader = BufReader::new(File::open(directory.join("incidences.jsonl"))?);
    while let Some(h) = next::<Incidence>(&mut reader)? {
        require(
            previous.is_none_or(|key| key < h.key()),
            "unordered orientation incidence stream",
        )?;
        previous = Some(h.key());
        if feature != Some(h.feature) {
            accumulate(&signs, &references, &mut support);
            signs.clear();
            feature = Some(h.feature);
        }
        // Full containment in an original description is impossible for a
        // crossing location under the validated fixed ownership map.
        if h.signed_mask != 0 {
            if let Some(g) = h.group {
                if let Some(index) = index.get(&(g, h.source)) {
                    index.containing(catalog, h.start, h.end, |id| {
                        *signs.entry(g).or_default().entry(id).or_default() |= h.signed_mask;
                    });
                }
            }
        }
    }
    accumulate(&signs, &references, &mut support);
    p.orientation_references = references.values().copied().collect();
    p.orientation_references.sort_unstable();
    p.orientation_support = support
        .into_iter()
        .enumerate()
        .filter(|(_, mask)| *mask != 0)
        .map(|(id, mask)| (references[&catalog.occurrences[id].group], id, mask))
        .collect();
    p.orientation_support.sort_unstable();
    Ok(())
}

/// Locate exactly one top-level provenance value while validating JSON tokens.
/// IgnoredAny avoids allocating intervening fields and byte offsets let late
/// rethread preserve every original score/ranking/rate byte verbatim.
fn provenance_span(bytes: &[u8]) -> io::Result<(usize, usize)> {
    use serde::de::IgnoredAny;
    let skip = |mut i: usize| {
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        i
    };
    let mut i = skip(0);
    require(bytes.get(i) == Some(&b'{'), "calls must be a JSON object")?;
    i += 1;
    let mut found = None;
    loop {
        i = skip(i);
        if bytes.get(i) == Some(&b'}') {
            i = skip(i + 1);
            break;
        }
        let mut key = serde_json::Deserializer::from_slice(&bytes[i..]).into_iter::<String>();
        let name = key
            .next()
            .transpose()
            .map_err(io::Error::other)?
            .ok_or_else(|| invalid("missing calls field"))?;
        i = skip(i + key.byte_offset());
        require(bytes.get(i) == Some(&b':'), "missing calls field colon")?;
        i = skip(i + 1);
        let start = i;
        let mut value = serde_json::Deserializer::from_slice(&bytes[i..]).into_iter::<IgnoredAny>();
        value
            .next()
            .transpose()
            .map_err(io::Error::other)?
            .ok_or_else(|| invalid("missing calls value"))?;
        i += value.byte_offset();
        if name == "observations" {
            require(
                found.replace((start, i)).is_none(),
                "duplicate observation provenance",
            )?;
        }
        i = skip(i);
        match bytes.get(i) {
            Some(b',') => {
                i = skip(i + 1);
                require(bytes.get(i) != Some(&b'}'), "trailing calls comma")?;
            }
            Some(b'}') => {
                i = skip(i + 1);
                break;
            }
            _ => return Err(invalid("invalid calls object delimiter")),
        }
    }
    require(i == bytes.len(), "trailing calls content")?;
    found.ok_or_else(|| invalid("missing calls observation provenance"))
}

pub fn write_derived_calls(
    out: &Path,
    original_path: &Path,
    original: &[u8],
    calls: &mut genotype::Genotypes,
) -> io::Result<()> {
    let p = calls
        .observations
        .as_mut()
        .ok_or_else(|| invalid("missing derived call provenance"))?;
    p.derived_from_calls = Some(
        json!({"path":original_path,"bytes":original.len(),"fnv1a64":format!("{:016x}",super::super::checksum(original))}),
    );
    let (start, end) = provenance_span(original)?;
    let value = serde_json::to_vec_pretty(p).map_err(io::Error::other)?;
    let mut bytes = Vec::with_capacity(original.len() - (end - start) + value.len());
    bytes.extend_from_slice(&original[..start]);
    bytes.extend(value);
    bytes.extend_from_slice(&original[end..]);
    super::super::atomic_write(out, &bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn provenance_span_preserves_exact_numeric_lexemes_and_ignores_nested_fields() {
        let raw = br#" {"score":-15174.237639340325,"nested":{"observations":2},"observations":{"x":"}\\\""},"rate":1e-123} "#;
        let (a, b) = provenance_span(raw).unwrap();
        assert_eq!(&raw[a..b], br#"{"x":"}\\\""}"#);
        assert!(provenance_span(br#"{"observations":{},"observations":{}}"#).is_err());
        assert!(provenance_span(br#"{"observations":{},}"#).is_err());
    }
}

#[cfg(test)]
mod axis_tests {
    use super::super::consumer::Provenance;
    use super::*;
    fn fixture(members: usize) -> (Metadata, genotype::Genotypes) {
        let identity = PanelIdentity {
            checksum_algorithm: "fixture".into(),
            sidecars: vec![],
        };
        let sources: Vec<_> = (0..=members)
            .map(|id| catalog::Source {
                id,
                path: format!("S{id}#0#chr"),
                length: 100,
            })
            .collect();
        let occurrences = (0..=members)
            .map(|id| catalog::Occurrence {
                id,
                source: id,
                group: usize::from(id == members),
                interval: catalog::Interval {
                    path: sources[id].path.clone(),
                    start: 0,
                    end: 100,
                    strand: None,
                },
                fully_contained_anchors: 0,
                anchor_status: String::new(),
                feature_multiplicities: BTreeMap::new(),
            })
            .collect();
        let catalog = catalog::Catalog {
            version: 1,
            panel: identity.clone(),
            catalog_accepted: false,
            validation: String::new(),
            sources,
            occurrences,
            groups: vec![
                catalog::Group {
                    id: "many".into(),
                    scaffold: None,
                    occurrences: (0..members).collect(),
                },
                catalog::Group {
                    id: "zero".into(),
                    scaffold: None,
                    occurrences: vec![members],
                },
            ],
            features: vec![],
            links: vec![],
        };
        let sample = sample::SampleIndex {
            version: 1,
            panel: identity,
            count_policy: COUNT_POLICY.into(),
            stats: sample::SampleStats {
                read_lengths: BTreeMap::from([(100, 1)]),
                ..Default::default()
            },
            counts: crate::sample_mem_bwt::WeightedBwt::build(&BTreeMap::new()).unwrap(),
        };
        let mut calls = genotype::call(
            &catalog,
            &sample,
            genotype::Parameters {
                haploid_depth: 10.0,
                background: 0.1,
                max_mean_deviance: 10.0,
            },
            true,
            1,
        )
        .unwrap();
        calls.model = MODEL.into();
        calls.sample_payload_checksum = Some("1234567890123456".into());
        calls.artifact_checksum_algorithm = CHECKSUM_ALGORITHM.into();
        calls.observations = Some(Provenance {
            metadata_checksum: "1234567890123456".into(),
            compiler_identity: compiler_identity(),
            full_feature_scope: true,
            feature_groups: vec![],
            legacy_input: json!({}),
            source_files: vec![],
            orientation_references: vec![],
            orientation_support: vec![],
            derived_from_calls: None,
        });
        let metadata = Metadata {
            version: 1,
            model: MODEL.into(),
            count_policy: COUNT_POLICY.into(),
            compiler_identity: compiler_identity(),
            capability: String::new(),
            catalog,
            legacy_input: json!({}),
            sources: vec![],
            syncmer_length_bp: 63,
            read_lengths: vec![100],
            feature_groups: vec![],
            full_feature_scope: true,
            feature_count: 2,
            files: BTreeMap::new(),
            stats: json!({}),
        };
        (metadata, calls)
    }
    fn row(group: &str, reference: usize) -> threading::AxisInterval {
        threading::AxisInterval {
            component: group.into(),
            start: 0,
            end: 100,
            group: group.into(),
            reference_occurrence: reference,
            reference_strand: "+".into(),
            orientations: BTreeMap::new(),
        }
    }
    #[test]
    fn many_member_axis_evidence_is_linear_direct_and_records_assessed_zero() {
        let n = 2048;
        let (metadata, mut calls) = fixture(n);
        let temp = tempfile::tempdir().unwrap();
        let hit = |feature, source, mask| Incidence {
            feature,
            source,
            start: 10,
            end: 80,
            signed_mask: mask,
            group: Some(0),
            context_nonlocal: true,
            source_terminal_context: true,
            contributions: vec![[0, 0]],
        };
        let mut hits: Vec<_> = (0..n)
            .map(|id| {
                hit(
                    0,
                    id,
                    match id {
                        1 => 2,
                        2 => 3,
                        3 => 0,
                        _ => 1,
                    },
                )
            })
            .collect();
        // These members share another feature, but the requested reference does
        // not. It must not create transitive or majority orientation evidence.
        hits.extend([hit(1, 1, 1), hit(1, 2, 2)]);
        write_lines(&temp.path().join("incidences.jsonl"), hits).unwrap();
        let axis = threading::Axis {
            version: 1,
            coordinate_system: "fixture".into(),
            intervals: vec![row("many", 0), row("zero", n)],
        };
        assert!(
            threading::thread(
                &metadata.catalog,
                &calls,
                threading::Axis {
                    version: 1,
                    coordinate_system: "fixture".into(),
                    intervals: axis.intervals.clone()
                },
                0.0
            )
            .is_err(),
            "uncomputed provider must not mean zero support"
        );
        prepare_orientations(&metadata, temp.path(), &mut calls, &axis).unwrap();
        let p = calls.observations.as_ref().unwrap();
        assert_eq!(p.orientation_references, vec![0, n]);
        assert_eq!(p.orientation_support.len(), n - 1);
        assert!(p
            .orientation_support
            .iter()
            .all(|&(reference, _, _)| reference == 0));
        assert!(p.orientation_support.contains(&(0, 1, 2)));
        assert!(p.orientation_support.contains(&(0, 2, 3)));
        assert!(!p
            .orientation_support
            .iter()
            .any(|&(_, id, _)| id == 3 || id == n));
        threading::thread(&metadata.catalog, &calls, axis, 0.0).unwrap();
        let repeated = threading::Axis {
            version: 1,
            coordinate_system: "fixture".into(),
            intervals: vec![
                row("many", 0),
                threading::AxisInterval {
                    start: 100,
                    end: 200,
                    ..row("many", 0)
                },
                row("zero", n),
            ],
        };
        prepare_orientations(&metadata, temp.path(), &mut calls, &repeated).unwrap();
        assert_eq!(
            calls.observations.as_ref().unwrap().orientation_references,
            vec![n]
        );
        assert!(calls
            .observations
            .as_ref()
            .unwrap()
            .orientation_support
            .is_empty());
        let threads = threading::thread(&metadata.catalog, &calls, repeated, 0.0).unwrap();
        assert_eq!(
            threads.intervals[0].status,
            "unresolved-repeated-axis-group"
        );
        assert_eq!(
            threads.intervals[1].status,
            "unresolved-repeated-axis-group"
        );
    }
    #[test]
    fn interval_index_matches_containment_with_nested_duplicate_and_disjoint_descriptions() {
        let (mut metadata, _) = fixture(1000);
        for (id, o) in metadata.catalog.occurrences.iter_mut().enumerate() {
            o.interval.start = (id % 100) as u64;
            o.interval.end = o.interval.start + (id % 17 + 1) as u64;
        }
        let index = IntervalIndex::new((0..1000).collect(), &metadata.catalog);
        for (start, end) in [(0, 1), (50, 51), (45, 60), (99, 115), (110, 116)] {
            let mut found = BTreeSet::new();
            index.containing(&metadata.catalog, start, end, |id| {
                found.insert(id);
            });
            let expected = metadata.catalog.occurrences[..1000]
                .iter()
                .filter(|o| o.interval.start <= start && o.interval.end >= end)
                .map(|o| o.id)
                .collect();
            assert_eq!(found, expected);
        }
    }
    #[test]
    fn selected_partial_group_with_only_exclusions_is_assessed_not_outside_scope() {
        let (mut metadata, _) = fixture(2);
        metadata.full_feature_scope = false;
        metadata.feature_groups = vec!["many".into()];
        let temp = tempfile::tempdir().unwrap();
        let tokens = canonical(&encode_walk(&[(1, 0), (2, 7)]).unwrap());
        let other = canonical(&encode_walk(&[(3, 0), (4, 7)]).unwrap());
        write_lines(
            &temp.path().join("definitions.jsonl"),
            [
                input::Definition {
                    id: 0,
                    tokens: tokens.clone().try_into().unwrap(),
                    original_owner: Some(0),
                },
                input::Definition {
                    id: 1,
                    tokens: other.try_into().unwrap(),
                    original_owner: Some(0),
                },
            ],
        )
        .unwrap();
        write_lines(
            &temp.path().join("incidences.jsonl"),
            (0..2).map(|id| Incidence {
                feature: id,
                source: id,
                start: 10,
                end: 80,
                signed_mask: 1,
                group: Some(0),
                context_nonlocal: true,
                source_terminal_context: true,
                contributions: vec![[7, 7]],
            }),
        )
        .unwrap();
        let sample = sample::SampleIndex {
            version: 1,
            panel: metadata.catalog.panel.clone(),
            count_policy: COUNT_POLICY.into(),
            stats: sample::SampleStats {
                read_lengths: BTreeMap::from([(100, 1)]),
                ..Default::default()
            },
            counts: crate::sample_mem_bwt::WeightedBwt::build(&BTreeMap::from([(tokens, 9)]))
                .unwrap(),
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
        assert_eq!(calls.calls[0].retained_features, 0);
        assert_eq!(calls.calls[0].status, "no-call-no-local-features");
        assert_eq!(
            calls.calls[1].status,
            "no-call-outside-assessed-feature-scope"
        );
        assert_eq!(calls.excluded_features, 2);
        assert!(calls
            .observations
            .as_ref()
            .unwrap()
            .orientation_support
            .is_empty());
    }
    #[test]
    fn derived_calls_preserve_all_non_orientation_bytes_and_link_original_content() {
        let (_, mut calls) = fixture(2);
        let raw = String::from_utf8(serde_json::to_vec_pretty(&calls).unwrap())
            .unwrap()
            .replace("\"haploid_depth\": 10.0", "\"haploid_depth\": 1e1")
            .into_bytes();
        assert!(String::from_utf8_lossy(&raw).contains("1e1"));
        let (a, b) = provenance_span(&raw).unwrap();
        let temp = tempfile::tempdir().unwrap();
        let original = temp.path().join("original.json");
        let derived = temp.path().join("derived.json");
        fs::write(&original, &raw).unwrap();
        calls.observations.as_mut().unwrap().orientation_references = vec![0];
        calls.observations.as_mut().unwrap().orientation_support = vec![(0, 0, 1)];
        write_derived_calls(&derived, &original, &raw, &mut calls).unwrap();
        let output = fs::read(&derived).unwrap();
        let (c, d) = provenance_span(&output).unwrap();
        assert_eq!(&raw[..a], &output[..c]);
        assert_eq!(&raw[b..], &output[d..]);
        assert_eq!(fs::read(&original).unwrap(), raw);
        let link = calls
            .observations
            .as_ref()
            .unwrap()
            .derived_from_calls
            .as_ref()
            .unwrap();
        assert_eq!(link["fnv1a64"], fingerprint(&original).unwrap()["fnv1a64"]);
        assert_eq!(link["bytes"], raw.len());
    }
    #[test]
    fn reverse_index_epsilon_ties_reconstruct_but_duplicate_best_indices_fail() {
        let (metadata, mut calls) = fixture(2);
        calls.sample_payload_checksum = Some("1234567890123456".into());
        calls.artifact_checksum_algorithm = CHECKSUM_ALGORITHM.into();
        let provenance = calls.observations.as_mut().unwrap();
        provenance.orientation_references = vec![0, 2];
        provenance.orientation_support = vec![(0, 0, 1), (0, 1, 1), (2, 2, 1)];
        for row in &mut calls.calls {
            row.retained_features = 1;
            row.positive_features = 1;
            row.observed_total = 1.0;
            row.background_only_score = 5.0;
            for (index, bundle) in row.bundles.iter_mut().enumerate() {
                bundle.positive_features = 1;
                bundle.score = if index == 0 { 1.0 + 5e-9 } else { 1.0 };
            }
            (row.status, row.best_bundles, row.score_gap) = classify(row, &calls.parameters);
        }
        assert_eq!(calls.calls[0].status, "tied");
        assert_eq!(calls.calls[0].best_bundles, vec![1, 0]);
        let axis = threading::Axis {
            version: 1,
            coordinate_system: "fixture".into(),
            intervals: vec![row("many", 0), row("zero", 2)],
        };
        let threads = threading::thread(&metadata.catalog, &calls, axis, 10.0).unwrap();
        let lengths = metadata
            .catalog
            .sources
            .iter()
            .map(|s| (s.path.clone(), s.length))
            .collect();
        crate::genome_inference::reconstruction::validate(&calls, &threads, &lengths).unwrap();
        calls.calls[0].best_bundles = vec![1, 0, 0];
        assert!(
            crate::genome_inference::reconstruction::validate(&calls, &threads, &lengths).is_err()
        );
    }
    #[test]
    fn background_equivalent_status_blocks_reconstruction_replay() {
        let (metadata, mut calls) = fixture(2);
        calls.sample_payload_checksum = Some("1234567890123456".into());
        calls.artifact_checksum_algorithm = CHECKSUM_ALGORITHM.into();
        calls.observations.as_mut().unwrap().orientation_references = vec![0, 2];
        for row in &mut calls.calls {
            row.retained_features = 1;
            row.positive_features = 1;
            row.observed_total = 1.0;
            row.background_only_score = 5.0;
            for b in &mut row.bundles {
                b.score = 5.0;
            }
            (row.status, row.best_bundles, row.score_gap) = classify(row, &calls.parameters);
            assert_eq!(row.status, "ambiguous-background-equivalent");
            assert!(row.best_bundles.is_empty());
        }
        let axis = threading::Axis {
            version: 1,
            coordinate_system: "fixture".into(),
            intervals: vec![row("many", 0), row("zero", 2)],
        };
        let threads = threading::thread(&metadata.catalog, &calls, axis, 10.0).unwrap();
        assert!(threads
            .intervals
            .iter()
            .all(|r| r.states.is_empty() && r.optimal_states.is_empty() && r.segment.is_none()));
        let lengths = metadata
            .catalog
            .sources
            .iter()
            .map(|s| (s.path.clone(), s.length))
            .collect();
        crate::genome_inference::reconstruction::validate(&calls, &threads, &lengths).unwrap();
    }
}
