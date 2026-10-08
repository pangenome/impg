use super::{
    catalog::{Catalog, Interval, Scaffold},
    sample::SampleIndex,
    FORMAT_VERSION, MODEL,
};
use crate::sample_mem_bwt::invalid;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::io;

#[derive(Debug, Serialize, Deserialize)]
pub struct DiagnosticCall {
    pub group: String,
    pub scaffold: Option<Scaffold>,
    pub status: String,
    pub compatible_source_occurrences: Vec<usize>,
    pub positive_features: Vec<usize>,
    pub retained_features: usize,
    pub unsupported_occurrences: usize,
    pub provisional: bool,
    pub locus_and_copy_structure_unresolved: bool,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct Calls {
    pub version: u32,
    pub model: String,
    pub experimental: bool,
    pub catalog_accepted: bool,
    pub interpretation: String,
    pub global_feature_counts: Vec<u64>,
    pub global_factors_used: usize,
    pub excluded_features_by_reason: BTreeMap<String, usize>,
    pub calls: Vec<DiagnosticCall>,
}

pub fn call(
    catalog: &Catalog,
    sample: &SampleIndex,
    allow_unvalidated: bool,
    ploidy: usize,
) -> io::Result<Calls> {
    if ploidy != 1 {
        return Err(invalid(
            "only haploid bootstrap diagnostics are supported; no mixture/dosage model",
        ));
    }
    if !allow_unvalidated {
        return Err(invalid(
            "ownership groups are unvalidated: explicit --allow-unvalidated-catalog required",
        ));
    }
    catalog.validate()?;
    if catalog.panel != sample.panel {
        return Err(invalid("sample/catalog panel dictionary mismatch"));
    }
    let counts = catalog
        .features
        .iter()
        .map(|f| sample.counts.count(&f.tokens))
        .collect::<io::Result<Vec<_>>>()?;
    let mut by_group = vec![Vec::new(); catalog.groups.len()];
    let mut excluded = BTreeMap::new();
    for (id, f) in catalog.features.iter().enumerate() {
        if let Some(g) = f.owning_group {
            by_group[g].push(id);
        }
        for reason in &f.exclusion_reasons {
            *excluded.entry(reason.clone()).or_default() += 1;
        }
    }
    let mut calls = Vec::new();
    for (g, group) in catalog.groups.iter().enumerate() {
        let positive: Vec<_> = by_group[g]
            .iter()
            .copied()
            .filter(|&f| counts[f] > 0)
            .collect();
        let compatible: Vec<_> = if positive.is_empty() {
            Vec::new()
        } else {
            group
                .occurrences
                .iter()
                .copied()
                .filter(|&id| {
                    positive.iter().all(|f| {
                        catalog.occurrences[id]
                            .feature_multiplicities
                            .contains_key(f)
                    })
                })
                .collect()
        };
        let status = if by_group[g].is_empty() {
            "no-call-no-local-features"
        } else if positive.is_empty() {
            "no-call-no-positive-features"
        } else if compatible.is_empty() {
            "conflict-inadequate-locus-or-evidence-model"
        } else if compatible.len() == 1 {
            "unique-compatible-occurrence"
        } else {
            "ambiguous-compatible-occurrences"
        };
        calls.push(DiagnosticCall {
            group: group.id.clone(),
            scaffold: group.scaffold.clone(),
            status: status.into(),
            compatible_source_occurrences: compatible,
            positive_features: positive,
            retained_features: by_group[g].len(),
            unsupported_occurrences: group
                .occurrences
                .iter()
                .filter(|&&id| catalog.occurrences[id].fully_contained_anchors < 2)
                .count(),
            provisional: true,
            locus_and_copy_structure_unresolved: true,
        });
    }
    // Order only by explicitly supplied evaluation scaffold, never source-name
    // parsing or presumed homologous numeric coordinates. No transitions/stitching.
    calls.sort_by(|a, b| match (&a.scaffold, &b.scaffold) {
        (Some(x), Some(y)) => {
            (&x.component, x.start, x.end, &a.group).cmp(&(&y.component, y.start, y.end, &b.group))
        }
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => a.group.cmp(&b.group),
    });
    Ok(Calls { version: FORMAT_VERSION, model: MODEL.into(), experimental: true, catalog_accepted: false,
        interpretation: "source-occurrence compatibility only; no genotype, posterior, dosage, homology certification, phasing or assembly; positives ignore magnitude/exposure/absence/error".into(),
        global_feature_counts: counts, global_factors_used: by_group.iter().map(Vec::len).sum(), excluded_features_by_reason: excluded, calls })
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Truth {
    pub version: u32,
    /// Names the EXPLICIT group scaffold coordinate system; not the panel sum.
    pub coordinate_system: Option<String>,
    pub groups: Vec<TruthGroup>,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TruthGroup {
    pub group: String,
    pub occurrences: Vec<Interval>,
}
fn same_interval(a: &Interval, b: &Interval) -> bool {
    a.path == b.path && a.start == b.start && a.end == b.end && a.strand == b.strand
}
fn union_bp(scaffolds: Vec<&Scaffold>) -> io::Result<u64> {
    let mut components: BTreeMap<&str, Vec<(u64, u64)>> = BTreeMap::new();
    for s in scaffolds {
        components
            .entry(&s.component)
            .or_default()
            .push((s.start, s.end));
    }
    let mut total = 0u64;
    for ranges in components.values_mut() {
        ranges.sort();
        let mut end = 0;
        for &(s, e) in ranges.iter() {
            total = total
                .checked_add(e.saturating_sub(s.max(end)))
                .ok_or_else(|| invalid("scaffold union base-pair count overflow"))?;
            end = end.max(e);
        }
    }
    Ok(total)
}

/// Evaluation consumes frozen calls. Truth is not a parameter to discovery,
/// sample indexing, feature selection, compatibility or ordering.
pub fn evaluate(catalog: &Catalog, calls: &Calls, truth: Truth) -> io::Result<serde_json::Value> {
    if truth.version != FORMAT_VERSION || truth.groups.is_empty() {
        return Err(invalid("invalid/empty truth"));
    }
    let mut seen = BTreeSet::new();
    let mut rows = Vec::new();
    let mut statuses = BTreeMap::new();
    let mut evaluated = Vec::new();
    let mut supported = Vec::new();
    let mut bp_available = truth
        .coordinate_system
        .as_ref()
        .is_some_and(|s| !s.trim().is_empty());
    for t in &truth.groups {
        if !seen.insert(&t.group) || t.occurrences.is_empty() {
            return Err(invalid("duplicate/empty truth group"));
        }
        let group = catalog
            .groups
            .iter()
            .find(|g| g.id == t.group)
            .ok_or_else(|| invalid("unknown truth group"))?;
        let call = calls
            .calls
            .iter()
            .find(|c| c.group == t.group)
            .ok_or_else(|| invalid("missing frozen call"))?;
        for interval in &t.occurrences {
            let source = catalog
                .sources
                .iter()
                .find(|s| s.path == interval.path)
                .ok_or_else(|| invalid("unknown truth source"))?;
            if interval.start >= interval.end
                || interval.end > source.length
                || interval
                    .strand
                    .as_deref()
                    .is_some_and(|s| s != "+" && s != "-")
            {
                return Err(invalid("invalid truth interval"));
            }
        }
        let panel_contains = t.occurrences.iter().all(|i| {
            group
                .occurrences
                .iter()
                .any(|&id| same_interval(i, &catalog.occurrences[id].interval))
        });
        let truth_contained = t.occurrences.iter().all(|i| {
            call.compatible_source_occurrences
                .iter()
                .any(|&id| same_interval(i, &catalog.occurrences[id].interval))
        });
        let exact_unique = t.occurrences.len() == 1
            && call.compatible_source_occurrences.len() == 1
            && truth_contained;
        *statuses.entry(call.status.clone()).or_insert(0usize) += 1;
        if let Some(s) = &group.scaffold {
            evaluated.push(s);
            if truth_contained {
                supported.push(s);
            }
        } else {
            bp_available = false;
        }
        rows.push(serde_json::json!({"group": t.group, "status": call.status,
            "truth_set_contained_in_catalog": panel_contains, "truth_set_contained_in_compatible_occurrences": truth_contained,
            "exact_unique_source_occurrence_recovery": exact_unique}));
    }
    Ok(
        serde_json::json!({"version": FORMAT_VERSION, "evaluation": "provisional-source-occurrence-recovery-not-genotype-accuracy",
        "coordinate_system": truth.coordinate_system, "groups": rows, "status_counts": statuses,
        "evaluated_scaffold_union_bp": if bp_available { Some(union_bp(evaluated)?) } else { None },
        "truth_compatible_scaffold_union_bp": if bp_available { Some(union_bp(supported)?) } else { None },
        "genotype_callable_bp": serde_json::Value::Null,
        "bp_scope": "only explicit truth-evaluated scaffold intervals; not whole sample genome unless independently certified complete"}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scaffold_union_resets_chromosomes_and_deduplicates_overlaps() {
        let a = Scaffold {
            component: "chr1".into(),
            start: 0,
            end: 10,
        };
        let b = Scaffold {
            component: "chr1".into(),
            start: 5,
            end: 20,
        };
        let c = Scaffold {
            component: "chr2".into(),
            start: 0,
            end: 10,
        };
        assert_eq!(union_bp(vec![&a, &a, &b, &c]).unwrap(), 30);
    }
    #[test]
    fn scaffold_union_rejects_cross_component_overflow() {
        let a = Scaffold {
            component: "chr1".into(),
            start: 0,
            end: u64::MAX,
        };
        let b = Scaffold {
            component: "chr2".into(),
            start: 0,
            end: 1,
        };
        assert_eq!(union_bp(vec![]).unwrap(), 0);
        assert_eq!(union_bp(vec![&a, &a]).unwrap(), u64::MAX);
        let error = union_bp(vec![&a, &b]).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error
            .to_string()
            .contains("scaffold union base-pair count overflow"));
    }
    #[test]
    fn boundary_phase_counterexample_is_not_a_local_dosage_factor() {
        let count_ab = |chromosomes: [&str; 2]| {
            chromosomes
                .iter()
                .map(|s| s.as_bytes().windows(2).filter(|w| *w == b"AB").count())
                .sum::<usize>()
        };
        assert_eq!(count_ab(["AB", "ab"]), 1);
        assert_eq!(count_ab(["Ab", "aB"]), 0);
        let marginal = |chromosomes: [&str; 2]| {
            let mut l = chromosomes.map(|s| s.as_bytes()[0]);
            let mut r = chromosomes.map(|s| s.as_bytes()[1]);
            l.sort();
            r.sort();
            (l, r)
        };
        assert_eq!(marginal(["AB", "ab"]), marginal(["Ab", "aB"]));
        // This caller excludes crossing spans; it never asserts discrimination.
    }
}
