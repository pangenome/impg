use super::*;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Provenance {
    pub metadata_checksum: String,
    pub compiler_identity: String,
    pub full_feature_scope: bool,
    pub feature_groups: Vec<String>,
    pub legacy_input: Value,
    pub source_files: Vec<Value>,
    /// Only references requested by non-repeated axis groups. Empty until an
    /// axis-specific provider is derived from the checked incidence stream.
    pub orientation_references: Vec<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub derived_from_calls: Option<Value>,
    /// Direct shared-context OR masks, never an all-member-pairs matrix.
    pub orientation_support: Vec<(usize, usize, u8)>,
}

pub fn eligible(call: &genotype::GenotypeCall, bundle: &genotype::Bundle) -> bool {
    bundle.positive_features > 0
        && call.background_only_score - bundle.score > genotype::TIE_EPSILON
}
/// Shared by the new caller and reconstruction validation. BG tie/preference is
/// unresolved, never absence. Legacy classification is intentionally untouched.
pub fn classify(
    row: &genotype::GenotypeCall,
    parameters: &genotype::Parameters,
) -> (String, Vec<usize>, Option<f64>) {
    let mut sources: Vec<_> = (0..row.bundles.len()).collect();
    sources.sort_by(|&a, &b| row.bundles[a].score.total_cmp(&row.bundles[b].score));
    // Preserve the original source-versus-source diagnostic, including sources
    // that cannot enter the DP. Eligibility must not redefine this number.
    let gap = sources
        .get(1)
        .map(|&b| row.bundles[b].score - row.bundles[sources[0]].score);
    let minimum = sources
        .first()
        .map(|&b| row.bundles[b].score)
        .unwrap_or(f64::INFINITY);
    let order: Vec<_> = sources
        .into_iter()
        .filter(|&b| eligible(row, &row.bundles[b]))
        .collect();
    let status = if row.retained_features == 0 {
        "no-call-no-local-features"
    } else if row.observed_total == 0.0 {
        "no-call-no-positive-features"
    } else if minimum - row.background_only_score > genotype::TIE_EPSILON {
        "no-call-background-preferred"
    } else if (minimum - row.background_only_score).abs() <= genotype::TIE_EPSILON {
        "ambiguous-background-equivalent"
    } else if order.is_empty() {
        "no-call-background-preferred"
    } else if row.bundles[order[0]].mean_deviance > parameters.max_mean_deviance {
        "poor-fit"
    } else if order.get(1).is_some_and(|&b| {
        row.bundles[b].score - row.bundles[order[0]].score <= genotype::TIE_EPSILON
    }) {
        "tied"
    } else {
        "informative"
    };
    let best = if status.starts_with("no-call") || status == "ambiguous-background-equivalent" {
        Vec::new()
    } else {
        let score = row.bundles[order[0]].score;
        order
            .into_iter()
            .filter(|&b| row.bundles[b].score - score <= genotype::TIE_EPSILON)
            .collect()
    };
    (status.into(), best, gap)
}
fn deviance(x: f64, mu: f64) -> f64 {
    2.0 * (mu - x + if x == 0.0 { 0.0 } else { x * (x / mu).ln() })
}

pub fn call(
    metadata: &Metadata,
    metadata_checksum: String,
    directory: &Path,
    sample: &sample::SampleIndex,
    sample_checksum: String,
    parameters: genotype::Parameters,
    out: &Path,
) -> io::Result<genotype::Genotypes> {
    parameters.validate()?;
    require(
        metadata.version == FORMAT_VERSION
            && metadata.model == MODEL
            && metadata.count_policy == COUNT_POLICY
            && metadata.compiler_identity == compiler_identity(),
        "incompatible observation model/compiler",
    )?;
    require(
        metadata.catalog.panel == sample.panel && sample.count_policy == COUNT_POLICY,
        "sample/observations mismatch",
    )?;
    genotype::exposure(&sample.stats.read_lengths, 1)?;
    require(
        sample
            .stats
            .read_lengths
            .keys()
            .all(|&l| metadata.read_lengths.contains(&(l as u64))),
        "sample contains uncompiled read length",
    )?;
    let denominator: f64 = sample
        .stats
        .read_lengths
        .iter()
        .map(|(&l, &n)| l as f64 * n as f64)
        .sum();
    let mut calls = genotype::call(&metadata.catalog, sample, parameters.clone(), true, 1)?;
    calls.model = MODEL.into();
    calls.sample_payload_checksum = Some(sample_checksum);
    calls.catalog_payload_checksum = None;
    calls.artifact_checksum_algorithm = CHECKSUM_ALGORITHM.into();
    calls.feature_details_included = false;
    calls.assumptions = "Error-free uniform source-start occurrence-weighted composite Poisson, not calibrated posterior. Integer surviving-substring profiles averaged across input orientations and read-length histogram; physical rates summed inside one factor logarithm. Complete panel physical support, conservative boundary/shared/context-nonlocal exclusion; no novel-junction, noisy-read or nonlocal solver. Background tie/preference unresolved, not deletion. Partial feature scope is not genome-wide inference.".into();
    let catalog = &metadata.catalog;
    let identities = catalog
        .sources
        .iter()
        .map(|s| genotype::source_identity(&s.path))
        .collect::<io::Result<Vec<_>>>()?;
    let bundle_lookup: Vec<BTreeMap<_, _>> = calls
        .calls
        .iter()
        .map(|c| {
            c.bundles
                .iter()
                .enumerate()
                .map(|(b, x)| (x.identity.clone(), b))
                .collect()
        })
        .collect();
    let mut base_deviance = vec![0.0; calls.calls.len()];
    let mut definitions = BufReader::new(File::open(directory.join("definitions.jsonl"))?);
    let mut incidences = BufReader::new(File::open(directory.join("incidences.jsonl"))?);
    let mut head: Option<Incidence> = next(&mut incidences)?;
    let temporary = out.join("factor-ledger.incomplete");
    let ledger_path = out.join("factor-ledger.jsonl");
    let mut ledger = BufWriter::new(File::create(&temporary)?);
    let mut previous_feature = None;
    let mut previous_location = None;
    let mut features = 0;
    while let Some(definition) = next::<input::Definition>(&mut definitions)? {
        require(
            previous_feature.is_none_or(|old| old < definition.id),
            "unordered feature IDs",
        )?;
        previous_feature = Some(definition.id);
        features += 1;
        let count = sample.counts.count(&definition.tokens)?; // exactly once, including excluded/zero features
        let mut reasons = BTreeSet::new();
        let mut groups = BTreeSet::new();
        let mut rates: BTreeMap<(usize, usize), (f64, u64)> = BTreeMap::new();
        let mut locations = 0u64;
        let mut zero_profiles = 0u64;
        while head.as_ref().is_some_and(|h| h.feature == definition.id) {
            let h = head.take().unwrap();
            require(
                previous_location.is_none_or(|p| p < h.key()),
                "unordered/duplicate physical incidence",
            )?;
            previous_location = Some(h.key());
            require(
                h.source < catalog.sources.len()
                    && h.start < h.end
                    && h.end <= catalog.sources[h.source].length
                    && h.end - h.start == definition.tokens[1] / 2 + metadata.syncmer_length_bp
                    && h.contributions.len() == metadata.read_lengths.len()
                    && h.group.is_none_or(|g| g < catalog.groups.len())
                    && h.signed_mask <= 3,
                "malformed physical incidence",
            )?;
            profile::add(&mut locations, 1)?;
            zero_profiles += u64::from(h.contributions.iter().all(|q| *q == [0, 0]));
            if h.context_nonlocal {
                reasons.insert("context-nonlocal".to_string());
            }
            if h.source_terminal_context {
                reasons.insert("source-terminal-context".to_string());
            }
            if let Some(g) = h.group {
                groups.insert(g);
                let b = *bundle_lookup[g]
                    .get(&identities[h.source])
                    .ok_or_else(|| invalid("physical source missing original bundle"))?;
                let mut numerator = 0.0;
                for (&length, q) in metadata.read_lengths.iter().zip(&h.contributions) {
                    if let Some(&n) = sample.stats.read_lengths.get(&(length as usize)) {
                        numerator += n as f64 * (q[0] as f64 / 2.0 + q[1] as f64 / 2.0);
                    }
                }
                let rate = rates.entry((g, b)).or_default();
                rate.0 += numerator / denominator;
                profile::add(&mut rate.1, 1)?;
            } else {
                reasons.insert("boundary-crossing".to_string());
            }
            head = next(&mut incidences)?;
        }
        require(
            head.as_ref().is_none_or(|h| h.feature > definition.id),
            "incidence references absent feature ID",
        )?;
        if locations == 0 {
            reasons.insert("unobservable-no-physical-incidence".into());
        }
        if groups.len() > 1 {
            reasons.insert("shared-between-groups".into());
        }
        let owner = if reasons.is_empty() && groups.len() == 1 {
            groups.first().copied()
        } else {
            None
        };
        if let Some(g) = owner {
            let row = &mut calls.calls[g];
            let x = count as f64;
            let bg = parameters.background;
            row.retained_features += 1;
            row.positive_features += usize::from(count > 0);
            row.observed_total += x;
            row.background_only_score += bg - x * bg.ln();
            base_deviance[g] += deviance(x, bg);
            for ((_, b), (e, m)) in rates {
                let bundle = &mut row.bundles[b];
                let signal = parameters.haploid_depth * e;
                bundle.score += signal - x * (signal / bg).ln_1p();
                bundle.expected_total += signal;
                bundle.modeled_features += 1;
                bundle.positive_features += usize::from(count > 0 && e > 0.0);
                profile::add(&mut bundle.physical_feature_locations, m)?;
            }
            calls.global_factors_used += 1;
        } else {
            calls.excluded_features += 1;
            for reason in &reasons {
                *calls
                    .excluded_features_by_reason
                    .entry(reason.clone())
                    .or_default() += 1;
            }
        }
        json_line(
            &mut ledger,
            &json!({"feature":definition.id,"count":count,"original_owner":definition.original_owner,
            "owner":owner,"exclusion_reasons":reasons,"physical_locations":locations,"zero_profiles":zero_profiles}),
        )?;
    }
    require(
        head.is_none() && features == metadata.feature_count,
        "incomplete feature/incidence stream",
    )?;
    finish(&mut ledger, &temporary, &ledger_path)?;
    for (g, row) in calls.calls.iter_mut().enumerate() {
        for bundle in &mut row.bundles {
            bundle.mean_deviance = (base_deviance[g] + 2.0 * bundle.score).max(0.0)
                / row.retained_features.max(1) as f64;
            bundle.score += row.background_only_score;
            bundle.expected_total += parameters.background * row.retained_features as f64;
            require(
                [bundle.score, bundle.mean_deviance, bundle.expected_total]
                    .iter()
                    .all(|x| x.is_finite()),
                "nonfinite observation score",
            )?;
        }
        (row.status, row.best_bundles, row.score_gap) = classify(row, &parameters);
        if !metadata.full_feature_scope
            && row.retained_features == 0
            && !metadata.feature_groups.contains(&row.group)
        {
            row.status = "no-call-outside-assessed-feature-scope".into();
        }
    }
    calls.observations = Some(Provenance {
        metadata_checksum,
        compiler_identity: metadata.compiler_identity.clone(),
        full_feature_scope: metadata.full_feature_scope,
        feature_groups: metadata.feature_groups.clone(),
        legacy_input: metadata.legacy_input.clone(),
        source_files: metadata.sources.clone(),
        orientation_references: Vec::new(),
        derived_from_calls: None,
        orientation_support: Vec::new(),
    });
    Ok(calls)
}

/// Check immutable call headers against the already-verified registry before
/// replacing axis-specific provenance. The CLI separately checks the registry's
/// content checksum; no incompatible input is repaired by orientation derivation.
pub(super) fn validate_frozen_calls(
    metadata: &Metadata,
    calls: &genotype::Genotypes,
) -> io::Result<()> {
    let digest = |s: &str| s.len() == 16 && s.bytes().all(|b| b.is_ascii_hexdigit());
    require(
        metadata.version == FORMAT_VERSION
            && metadata.model == MODEL
            && metadata.count_policy == COUNT_POLICY
            && calls.version == FORMAT_VERSION
            && calls.model == MODEL
            && calls.panel == metadata.catalog.panel
            && calls.count_policy == COUNT_POLICY
            && calls.ploidy == 1
            && calls.experimental
            && !calls.catalog_accepted
            && calls.catalog_payload_checksum.is_none()
            && calls.artifact_checksum_algorithm == CHECKSUM_ALGORITHM
            && calls.sample_payload_checksum.as_deref().is_some_and(digest),
        "incompatible frozen partition call headers",
    )?;
    let provenance = calls
        .observations
        .as_ref()
        .ok_or_else(|| invalid("missing observation provenance"))?;
    require(
        digest(&provenance.metadata_checksum)
            && provenance.compiler_identity == metadata.compiler_identity
            && provenance.full_feature_scope == metadata.full_feature_scope
            && provenance.feature_groups == metadata.feature_groups
            && provenance.legacy_input == metadata.legacy_input
            && provenance.source_files == metadata.sources,
        "incompatible frozen observation provenance",
    )?;
    calls.parameters.validate()?;
    genotype::exposure(&calls.read_lengths, 1)?;
    let supported_lengths: BTreeSet<_> = metadata.read_lengths.iter().copied().collect();
    require(
        calls
            .read_lengths
            .keys()
            .all(|&length| supported_lengths.contains(&(length as u64)))
            && calls.calls.len() == metadata.catalog.groups.len()
            && calls
                .calls
                .iter()
                .zip(&metadata.catalog.groups)
                .all(|(row, group)| row.group == group.id)
            && calls.source_occurrences.len() == metadata.catalog.occurrences.len(),
        "incompatible frozen partition call scope",
    )?;
    let identities = metadata
        .catalog
        .sources
        .iter()
        .map(|s| genotype::source_identity(&s.path))
        .collect::<io::Result<Vec<_>>>()?;
    require(
        calls
            .source_occurrences
            .iter()
            .zip(&metadata.catalog.occurrences)
            .all(|(actual, expected)| {
                actual.id == expected.id
                    && metadata
                        .catalog
                        .groups
                        .get(expected.group)
                        .is_some_and(|g| actual.group == g.id)
                    && identities
                        .get(expected.source)
                        .is_some_and(|identity| actual.identity == *identity)
                    && actual.interval.path == expected.interval.path
                    && actual.interval.start == expected.interval.start
                    && actual.interval.end == expected.interval.end
                    && actual.interval.strand == expected.interval.strand
            }),
        "incompatible frozen source occurrence provenance",
    )
}

pub fn thread(
    metadata: &Metadata,
    calls: &genotype::Genotypes,
    axis: threading::Axis,
    penalty: f64,
) -> io::Result<threading::Threads> {
    require(
        metadata.full_feature_scope,
        "partial feature scope cannot be threaded/reconstructed",
    )?;
    require(
        metadata.model == MODEL && calls.model == MODEL,
        "incompatible partition threading model",
    )?;
    validate_frozen_calls(metadata, calls)?;
    for (g, row) in calls.calls.iter().enumerate() {
        for bundle in &row.bundles {
            require(
                bundle.score.is_finite()
                    && bundle.mean_deviance.is_finite()
                    && bundle.mean_deviance >= 0.0
                    && !bundle.source_occurrences.is_empty()
                    && bundle.source_occurrences.iter().all(|&id| {
                        metadata.catalog.occurrences.get(id).is_some_and(|o| {
                            o.group == g && calls.source_occurrences[id].identity == bundle.identity
                        })
                    }),
                "invalid partition bundle/state membership",
            )?;
        }
        let (status, best, _) = classify(row, &calls.parameters);
        require(
            row.status == status && row.best_bundles == best,
            "inconsistent frozen partition call classification",
        )?;
    }
    threading::thread(&metadata.catalog, calls, axis, penalty)
}
