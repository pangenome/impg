//! Restricted haploid composite Poisson model. Scores are not calibrated posteriors.
use super::{
    catalog::{Catalog, Interval},
    sample::SampleIndex,
    PanelIdentity, COUNT_POLICY, FORMAT_VERSION,
};
use crate::sample_mem_bwt::invalid;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::io;

pub const MODEL: &str = "haploid-bundle-composite-poisson-v1";
pub const TIE_EPSILON: f64 = 1e-8;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Parameters {
    pub haploid_depth: f64,
    pub background: f64,
    /// Descriptive working-model diagnostic cutoff, NOT a significance threshold.
    pub max_mean_deviance: f64,
}
impl Parameters {
    pub fn validate(&self) -> io::Result<()> {
        if [self.haploid_depth, self.background, self.max_mean_deviance]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "depth, background and max-mean-deviance must be finite and positive",
            ));
        }
        Ok(())
    }
}

/// PanSN sample/haplotype only; the entire remaining contig spelling is retained.
/// Non-PanSN names are rejected, never interpreted as chromosomes or haplotypes.
pub fn source_identity(path: &str) -> io::Result<String> {
    let parts: Vec<_> = path.splitn(3, '#').collect();
    if parts.len() != 3
        || parts
            .iter()
            .any(|s| s.trim().is_empty() || s.chars().any(char::is_control))
    {
        return Err(invalid(format!(
            "quantitative mode requires PanSN sample#haplotype#contig: {path}"
        )));
    }
    Ok(format!("{}#{}", parts[0], parts[1]))
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Factor {
    pub feature: usize,
    pub count: u64,
    pub span_bp: u64,
    pub exposure: f64,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct Bundle {
    pub identity: String,
    pub source_occurrences: Vec<usize>,
    /// Absent in compact outputs, NOT an assertion of zero multiplicity.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub physical_multiplicities: BTreeMap<usize, u64>,
    pub physical_feature_locations: u64,
    pub modeled_features: usize,
    pub positive_features: usize,
    pub score: f64,
    pub mean_deviance: f64,
    pub expected_total: f64,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct GenotypeCall {
    pub group: String,
    pub status: String,
    /// Absent in compact outputs; exact factors remain reproducible from artifacts.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub factors: Vec<Factor>,
    pub retained_features: usize,
    pub positive_features: usize,
    pub observed_total: f64,
    pub background_only_score: f64,
    pub bundles: Vec<Bundle>,
    pub best_bundles: Vec<usize>,
    pub score_gap: Option<f64>,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct SourceOccurrence {
    pub id: usize,
    pub group: String,
    pub identity: String,
    pub interval: Interval,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct Genotypes {
    pub version: u32,
    pub model: String,
    pub panel: PanelIdentity,
    pub count_policy: String,
    pub sample_payload_checksum: Option<String>,
    pub catalog_payload_checksum: Option<String>,
    pub artifact_checksum_algorithm: String,
    pub feature_details_included: bool,
    pub experimental: bool,
    pub catalog_accepted: bool,
    pub ploidy: usize,
    pub parameters: Parameters,
    pub assumptions: String,
    pub read_lengths: BTreeMap<usize, u64>,
    pub global_factors_used: usize,
    pub excluded_features: usize,
    pub excluded_features_by_reason: BTreeMap<String, usize>,
    pub source_occurrences: Vec<SourceOccurrence>,
    pub calls: Vec<GenotypeCall>,
}

impl Genotypes {
    /// Keep summaries/provenance while omitting reproducible derived feature maps.
    /// Deserialization defaults missing maps to empty; consumers must check the
    /// explicit feature_details_included flag before interpreting those maps.
    pub fn omit_feature_details(&mut self) {
        self.feature_details_included = false;
        for call in &mut self.calls {
            call.factors.clear();
            call.factors.shrink_to_fit();
            for bundle in &mut call.bundles {
                bundle.physical_multiplicities.clear();
            }
        }
    }
}

pub fn exposure(histogram: &BTreeMap<usize, u64>, span: u64) -> io::Result<f64> {
    let mut bases = 0.0;
    let mut opportunities = 0.0;
    for (&length, &n) in histogram {
        if length == 0 || n == 0 {
            return Err(invalid("invalid read-length histogram"));
        }
        bases += length as f64 * n as f64;
        opportunities += (length as u64).saturating_sub(span.saturating_sub(1)) as f64 * n as f64;
    }
    if bases == 0.0 || !bases.is_finite() || !opportunities.is_finite() || span == 0 {
        return Err(invalid("empty/nonfinite read-span exposure"));
    }
    Ok(opportunities / bases)
}
fn deviance(count: f64, lambda: f64) -> f64 {
    2.0 * (lambda - count
        + if count == 0.0 {
            0.0
        } else {
            count * (count / lambda).ln()
        })
}
fn correction(count: f64, signal: f64, background: f64) -> f64 {
    signal - count * (signal / background).ln_1p()
}
pub fn call(
    catalog: &Catalog,
    sample: &SampleIndex,
    parameters: Parameters,
    allow_unvalidated: bool,
    ploidy: usize,
) -> io::Result<Genotypes> {
    parameters.validate()?;
    if !allow_unvalidated {
        return Err(invalid("explicit --allow-unvalidated-catalog required"));
    }
    if ploidy != 1 {
        return Err(invalid("quantitative mode supports ploidy 1 only"));
    }
    catalog.validate()?;
    if catalog.panel != sample.panel || sample.count_policy != COUNT_POLICY {
        return Err(invalid("sample/catalog panel or count policy mismatch"));
    }
    exposure(&sample.stats.read_lengths, 1)?;
    let identities = catalog
        .sources
        .iter()
        .map(|s| source_identity(&s.path))
        .collect::<io::Result<Vec<_>>>()?;
    let mut calls = Vec::new();
    let mut occurrence_bundle = vec![0; catalog.occurrences.len()];
    for group in &catalog.groups {
        let mut members: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for &id in &group.occurrences {
            members
                .entry(identities[catalog.occurrences[id].source].clone())
                .or_default()
                .push(id);
        }
        let bundles = members
            .into_iter()
            .enumerate()
            .map(|(b, (identity, source_occurrences))| {
                for &id in &source_occurrences {
                    occurrence_bundle[id] = b;
                }
                Bundle {
                    identity,
                    source_occurrences,
                    physical_multiplicities: BTreeMap::new(),
                    physical_feature_locations: 0,
                    modeled_features: 0,
                    positive_features: 0,
                    score: 0.0,
                    mean_deviance: 0.0,
                    expected_total: 0.0,
                }
            })
            .collect();
        calls.push(GenotypeCall {
            group: group.id.clone(),
            status: String::new(),
            factors: Vec::new(),
            retained_features: 0,
            positive_features: 0,
            observed_total: 0.0,
            background_only_score: 0.0,
            bundles,
            best_bundles: Vec::new(),
            score_gap: None,
        });
    }
    let mut excluded = BTreeMap::new();
    let mut excluded_features = 0;
    let mut exposure_cache = BTreeMap::new();
    for (f, feature) in catalog.features.iter().enumerate() {
        let Some(g) = feature.owning_group else {
            excluded_features += 1;
            if feature.exclusion_reasons.is_empty() {
                *excluded.entry("no-global-owner".into()).or_default() += 1;
            }
            for reason in &feature.exclusion_reasons {
                *excluded.entry(reason.clone()).or_default() += 1;
            }
            continue;
        };
        if !feature.exclusion_reasons.is_empty() {
            return Err(invalid("owned feature has exclusion reasons"));
        }
        let span = feature
            .locations
            .first()
            .ok_or_else(|| invalid("owned feature without locations"))?
            .end
            - feature.locations[0].start;
        let e = match exposure_cache.entry(span) {
            std::collections::btree_map::Entry::Occupied(entry) => *entry.get(),
            std::collections::btree_map::Entry::Vacant(entry) => {
                *entry.insert(exposure(&sample.stats.read_lengths, span)?)
            }
        };
        let count = sample.counts.count(&feature.tokens)?;
        let row = &mut calls[g];
        row.factors.push(Factor {
            feature: f,
            count,
            span_bp: span,
            exposure: e,
        });
        row.retained_features += 1;
        row.positive_features += usize::from(count > 0);
        // One physical span counts once, regardless of duplicate/overlapping BED representations.
        let mut physical = BTreeSet::new();
        let mut multiplicities: BTreeMap<usize, u64> = BTreeMap::new();
        for location in &feature.locations {
            if location.end - location.start != span || location.containing_occurrences.is_empty() {
                return Err(invalid("inconsistent owned feature span/containment"));
            }
            let mut bundle = None;
            for &id in &location.containing_occurrences {
                let o = &catalog.occurrences[id];
                if o.group != g
                    || o.source != location.source
                    || location.start < o.interval.start
                    || location.end > o.interval.end
                {
                    return Err(invalid("invalid physical feature ownership"));
                }
                let b = occurrence_bundle[id];
                if bundle.is_some_and(|old| old != b) {
                    return Err(invalid("physical location has multiple source identities"));
                }
                bundle = Some(b);
            }
            if physical.insert((location.source, location.start, location.end)) {
                *multiplicities.entry(bundle.unwrap()).or_default() += 1;
            }
        }
        let x = count as f64;
        let bg = parameters.background;
        row.observed_total += x;
        row.background_only_score += bg - x * bg.ln();
        // Only sparse nonzero state corrections. Zero observations still incur signal cost.
        for (b, m) in multiplicities {
            let bundle = &mut row.bundles[b];
            bundle.physical_multiplicities.insert(f, m);
            bundle.physical_feature_locations = bundle
                .physical_feature_locations
                .checked_add(m)
                .ok_or_else(|| invalid("physical feature count overflow"))?;
            bundle.modeled_features += 1;
            bundle.positive_features += usize::from(count > 0);
            let signal = parameters.haploid_depth * e * m as f64;
            bundle.score += correction(x, signal, bg);
            bundle.expected_total += signal;
        }
    }
    for row in &mut calls {
        let base_deviance: f64 = row
            .factors
            .iter()
            .map(|f| deviance(f.count as f64, parameters.background))
            .sum();
        for bundle in &mut row.bundles {
            bundle.mean_deviance =
                (base_deviance + 2.0 * bundle.score).max(0.0) / row.factors.len().max(1) as f64;
            bundle.score += row.background_only_score;
            bundle.expected_total += parameters.background * row.factors.len() as f64;
            if [bundle.score, bundle.mean_deviance, bundle.expected_total]
                .iter()
                .any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "nonfinite count score; parameters exceed working numeric range",
                ));
            }
        }
        let mut order: Vec<_> = (0..row.bundles.len()).collect();
        order.sort_by(|&a, &b| row.bundles[a].score.total_cmp(&row.bundles[b].score));
        let best = row.bundles[order[0]].score;
        row.score_gap = order.get(1).map(|&b| row.bundles[b].score - best);
        row.status = if row.factors.is_empty() {
            "no-call-no-local-features"
        } else if row.observed_total == 0.0 {
            "no-call-no-positive-features"
        } else if row.bundles[order[0]].mean_deviance > parameters.max_mean_deviance {
            "poor-fit"
        } else if row.score_gap.is_some_and(|g| g <= TIE_EPSILON) {
            "tied"
        } else {
            "informative"
        }
        .into();
        if !row.status.starts_with("no-call") {
            row.best_bundles = order
                .into_iter()
                .filter(|&b| row.bundles[b].score - best <= TIE_EPSILON)
                .collect();
        }
    }
    Ok(Genotypes { version: FORMAT_VERSION, model: MODEL.into(), panel: catalog.panel.clone(),
        count_policy: sample.count_policy.clone(), sample_payload_checksum: None, catalog_payload_checksum: None,
        artifact_checksum_algorithm: "fnv1a64-payload-v1; sample=bincode-payload; catalog=compact-json-payload".into(),
        feature_details_included: true, experimental: true, catalog_accepted: false, ploidy,
        parameters, assumptions: "Working composite Poisson, not calibrated posterior: ideal read-span exposure; MEM overlap/selection, error background and feature dependence uncalibrated. Globally owned within-group factors only, each scored once; shared/nonlocal/boundary factors excluded, not solved. Source bundles do not certify homology.".into(),
        read_lengths: sample.stats.read_lengths.clone(), global_factors_used: calls.iter().map(|c| c.factors.len()).sum(),
        excluded_features, excluded_features_by_reason: excluded,
        source_occurrences: catalog.occurrences.iter().map(|o| SourceOccurrence { id: o.id,
            group: catalog.groups[o.group].id.clone(), identity: identities[o.source].clone(), interval: o.interval.clone() }).collect(), calls })
}

/// Truth is consumed strictly downstream of frozen count calls and threads.
pub fn evaluate(
    catalog: &Catalog,
    calls: &Genotypes,
    threads: Option<&super::threading::Threads>,
    truth: super::calling::Truth,
) -> io::Result<serde_json::Value> {
    // Reuse the established interval validation/compatibility evaluator, but preserve
    // its diagnostic labels rather than claiming certified biological accuracy.
    let diagnostic = super::calling::Calls {
        version: FORMAT_VERSION,
        model: MODEL.into(),
        experimental: true,
        catalog_accepted: false,
        interpretation: "best source-bundle compatibility".into(),
        global_feature_counts: Vec::new(),
        global_factors_used: calls.global_factors_used,
        excluded_features_by_reason: calls.excluded_features_by_reason.clone(),
        calls: calls
            .calls
            .iter()
            .enumerate()
            .map(|(g, c)| super::calling::DiagnosticCall {
                group: c.group.clone(),
                scaffold: catalog.groups[g].scaffold.clone(),
                status: c.status.clone(),
                compatible_source_occurrences: if c.status == "poor-fit" {
                    Vec::new()
                } else {
                    c.best_bundles
                        .iter()
                        .flat_map(|&b| c.bundles[b].source_occurrences.clone())
                        .collect()
                },
                positive_features: Vec::new(),
                retained_features: c.retained_features,
                unsupported_occurrences: 0,
                provisional: true,
                locus_and_copy_structure_unresolved: true,
            })
            .collect(),
    };
    let same = |a: &Interval, b: &Interval| {
        a.path == b.path && a.start == b.start && a.end == b.end && a.strand == b.strand
    };
    let bundle_rows: Vec<_> = truth.groups.iter().map(|t| {
        let call = calls.calls.iter().find(|c| c.group == t.group);
        let compatible: Vec<_> = call.map(|c| c.best_bundles.iter().copied().filter(|&b| {
            c.status != "poor-fit" && t.occurrences.iter().all(|i| c.bundles[b].source_occurrences.iter().any(|&id| same(i, &catalog.occurrences[id].interval)))
        }).collect()).unwrap_or_default();
        let truth_identities = t.occurrences.iter().map(|i| source_identity(&i.path)).collect::<io::Result<BTreeSet<_>>>();
        serde_json::json!({"group": t.group, "compatible_best_bundles": compatible,
            "truth_set_contained_in_one_best_bundle": !compatible.is_empty(),
            "truth_bundle_identities": truth_identities.ok(),
            "unique_bundle_truth_compatible": call.is_some_and(|c| c.best_bundles.len() == 1 && !compatible.is_empty())})
    }).collect();
    let thread_evaluation = threads.map(|threads| {
        let truth_by_group: BTreeMap<_, _> = truth.groups.iter().map(|t| (t.group.as_str(), t)).collect();
        let mut assessed = 0usize;
        let mut compatible = 0usize;
        let mut known_switches = 0usize;
        let mut recovered_switches = 0usize;
        let mut reported_assessable_switches = 0usize;
        for row in &threads.intervals {
            if let Some(t) = truth_by_group.get(row.axis.group.as_str()) {
                assessed += 1;
                if row.status == "resolved-experimental" {
                    let state = &row.states[row.optimal_states[0]];
                    if t.occurrences.iter().any(|i| same(i, &state.interval)) { compatible += 1; }
                }
            }
        }
        for (p, pair) in threads.intervals.windows(2).enumerate() {
            if pair[0].axis.component != pair[1].axis.component { continue; }
            let (Some(left), Some(right)) = (truth_by_group.get(pair[0].axis.group.as_str()), truth_by_group.get(pair[1].axis.group.as_str())) else { continue; };
            if left.occurrences.len() != 1 || right.occurrences.len() != 1 || threads.repeated_axis_groups.contains(&left.group) || threads.repeated_axis_groups.contains(&right.group) { continue; }
            let known = left.occurrences[0].path != right.occurrences[0].path;
            let reported = threads.switches.iter().any(|s| s.left_interval == p && s.right_interval == p+1 && s.reason == "source-path-change");
            known_switches += usize::from(known);
            reported_assessable_switches += usize::from(reported);
            recovered_switches += usize::from(known && reported);
        }
        serde_json::json!({"truth_assessed_intervals": assessed, "resolved_truth_source_compatible_intervals": compatible,
            "known_source_path_switches": known_switches, "recovered_coarse_source_path_switches": recovered_switches,
            "reported_assessable_source_path_switches": reported_assessable_switches,
            "scope": "Exact source-interval compatibility; switches evaluated only between adjacent same-component single-occurrence truth groups, excluding repeated-axis groups. No biological recombination or exact breakpoint claim."})
    });
    let mut result = super::calling::evaluate(catalog, &diagnostic, truth)?;
    result["evaluation"] =
        "experimental-best-source-bundle-compatibility-not-validated-locus-biology".into();
    result["bundle_compatibility"] = serde_json::json!(bundle_rows);
    result["threading"] = serde_json::json!(thread_evaluation);
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn read_span_exposure_and_sparse_poisson_match_dense_including_zeros() {
        let histogram = BTreeMap::from([(100, 2), (200, 1)]);
        assert!((exposure(&histogram, 150).unwrap() - 51.0 / 400.0).abs() < 1e-12);
        assert_eq!(exposure(&histogram, 201).unwrap(), 0.0);
        let bg: f64 = 0.1;
        for counts in [[0.0, 0.0, 0.0], [1.0, 10.0, 30.0]] {
            for signals in [[0.0, 8.0, 0.0], [16.0, 8.0, 24.0]] {
                let dense: f64 = counts
                    .iter()
                    .zip(signals)
                    .map(|(&x, s)| s + bg - x * (s + bg).ln())
                    .sum();
                let sparse: f64 = counts.iter().map(|&x| bg - x * bg.ln()).sum::<f64>()
                    + counts
                        .iter()
                        .zip(signals)
                        .map(|(&x, s)| correction(x, s, bg))
                        .sum::<f64>();
                assert!((dense - sparse).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn bundles_deduplicate_overlapping_bed_rows_but_keep_physical_copies_and_change_calls() {
        use super::super::{catalog::*, sample::SampleStats};
        use crate::sample_mem_bwt::{canonical, encode_walk, WeightedBwt};
        let identity = PanelIdentity {
            checksum_algorithm: "fixture".into(),
            sidecars: Vec::new(),
        };
        let sources: Vec<_> = ["A#0#chr1", "B#0#copy1", "B#0#copy2"]
            .iter()
            .enumerate()
            .map(|(id, name)| Source {
                id,
                path: name.to_string(),
                length: 100,
            })
            .collect();
        let occurrences: Vec<_> = [0, 0, 0, 1, 2]
            .iter()
            .enumerate()
            .map(|(id, &source)| Occurrence {
                id,
                group: 0,
                source,
                interval: Interval {
                    path: sources[source].path.clone(),
                    start: if id == 2 { 5 } else { 0 },
                    end: 50,
                    strand: None,
                },
                fully_contained_anchors: 2,
                anchor_status: "fixture".into(),
                feature_multiplicities: BTreeMap::from([(0, 1)]),
            })
            .collect();
        let tokens = canonical(&encode_walk(&[(1, 10), (2, 15)]).unwrap());
        let catalog = Catalog {
            version: FORMAT_VERSION,
            panel: identity.clone(),
            catalog_accepted: false,
            validation: "fixture".into(),
            sources,
            groups: vec![Group {
                id: "copy-number".into(),
                scaffold: None,
                occurrences: (0..5).collect(),
            }],
            occurrences,
            features: vec![Feature {
                tokens: tokens.clone(),
                owning_group: Some(0),
                exclusion_reasons: Vec::new(),
                locations: vec![
                    FeatureLocation {
                        source: 0,
                        start: 10,
                        end: 20,
                        signed_nodes: [1, 2],
                        containing_occurrences: vec![0, 1, 2],
                    },
                    FeatureLocation {
                        source: 0,
                        start: 10,
                        end: 20,
                        signed_nodes: [1, 2],
                        containing_occurrences: vec![0, 1, 2],
                    },
                    FeatureLocation {
                        source: 1,
                        start: 10,
                        end: 20,
                        signed_nodes: [1, 2],
                        containing_occurrences: vec![3],
                    },
                    FeatureLocation {
                        source: 2,
                        start: 10,
                        end: 20,
                        signed_nodes: [1, 2],
                        containing_occurrences: vec![4],
                    },
                ],
            }],
            links: Vec::new(),
        };
        for (count, winner) in [(9, 0), (18, 1)] {
            let sample = SampleIndex {
                version: FORMAT_VERSION,
                panel: identity.clone(),
                count_policy: COUNT_POLICY.into(),
                stats: SampleStats {
                    read_lengths: BTreeMap::from([(100, 1)]),
                    ..Default::default()
                },
                counts: WeightedBwt::build(&BTreeMap::from([(tokens.clone(), count)])).unwrap(),
            };
            let result = call(
                &catalog,
                &sample,
                Parameters {
                    haploid_depth: 10.0,
                    background: 0.1,
                    max_mean_deviance: 10.0,
                },
                true,
                1,
            )
            .unwrap();
            let row = &result.calls[0];
            assert_eq!(result.global_factors_used, 1);
            assert_eq!(row.bundles[0].source_occurrences, vec![0, 1, 2]);
            assert_eq!(row.bundles[1].source_occurrences, vec![3, 4]);
            assert_eq!(row.bundles[0].physical_multiplicities[&0], 1);
            assert_eq!(row.bundles[1].physical_multiplicities[&0], 2);
            assert_eq!(row.best_bundles, vec![winner]);
            assert_eq!(row.status, "informative");
        }
    }
    #[test]
    fn absolute_count_magnitude_and_copy_multiplicity_change_the_winner() {
        let score = |x, m| correction(x, 10.0 * m, 0.1);
        assert!(score(10.0, 1.0) < score(10.0, 2.0));
        assert!(score(20.0, 2.0) < score(20.0, 1.0));
        assert!(score(0.0, 1.0) < score(0.0, 2.0));
        assert_eq!(source_identity("A#0#copy#extra").unwrap(), "A#0");
        for s in ["chr1", "A#0", "A##chr1", "#0#chr1", "A#0#"] {
            assert!(source_identity(s).is_err());
        }
        for x in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(Parameters {
                haploid_depth: x,
                background: 0.1,
                max_mean_deviance: 10.0
            }
            .validate()
            .is_err());
        }
    }
}
