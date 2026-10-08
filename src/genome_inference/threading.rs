//! Chromosome-reset source-occurrence DP on an independently declared reference axis.
use super::{
    catalog::{Catalog, Interval},
    genotype::{Genotypes, TIE_EPSILON},
    FORMAT_VERSION,
};
use crate::sample_mem_bwt::invalid;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::io;

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Axis {
    pub version: u32,
    pub coordinate_system: String,
    pub intervals: Vec<AxisInterval>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AxisInterval {
    pub component: String,
    pub start: u64,
    pub end: u64,
    pub group: String,
    pub reference_occurrence: usize,
    /// Reference source traversal relative to increasing axis coordinates.
    pub reference_strand: String,
    /// Optional externally established orientations relative to this axis.
    #[serde(default)]
    pub orientations: BTreeMap<usize, String>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct State {
    pub bundle: usize,
    pub identity: String,
    /// Duplicate BED rows retained as provenance, not duplicate physical states.
    pub source_occurrences: Vec<usize>,
    pub interval: Interval,
    pub strand: Option<String>,
    pub orientation_evidence: String,
    pub emission: f64,
    pub mean_deviance: f64,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct ThreadInterval {
    pub axis: AxisInterval,
    pub status: String,
    pub states: Vec<State>,
    /// All states on an optimal path within the reported absolute score tolerance.
    pub optimal_states: Vec<usize>,
    pub segment: Option<usize>,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct Segment {
    pub first_interval: usize,
    pub last_interval: usize,
    pub score: f64,
    pub tied: bool,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct Block {
    pub component: String,
    pub source_path: String,
    pub strand: String,
    pub interval_indices: Vec<usize>,
    pub source_gaps_bp: Vec<i64>,
    pub reference_gaps_bp: Vec<i64>,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct Switch {
    pub component: String,
    pub left_interval: usize,
    pub right_interval: usize,
    pub from_path: String,
    pub to_path: String,
    pub reason: String,
    pub breakpoint_start: u64,
    pub breakpoint_end: u64,
    pub reference_gap_bp: i64,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct Threads {
    pub version: u32,
    pub model: String,
    pub panel: super::PanelIdentity,
    pub count_policy: String,
    pub sample_payload_checksum: Option<String>,
    pub catalog_payload_checksum: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub observation_metadata_checksum: Option<String>,
    pub experimental: bool,
    pub catalog_accepted: bool,
    pub coordinate_system: String,
    pub switch_penalty: f64,
    pub tie_epsilon: f64,
    pub interpretation: String,
    pub intervals: Vec<ThreadInterval>,
    pub segments: Vec<Segment>,
    pub blocks: Vec<Block>,
    pub switches: Vec<Switch>,
    pub unscaffolded_groups: Vec<String>,
    pub repeated_axis_groups: Vec<String>,
    pub resolved_intervals: usize,
    pub unresolved_intervals: usize,
    pub axis_union_bp: u64,
    pub resolved_union_bp: u64,
}
fn multiply_strands(a: &str, b: &str) -> String {
    if a == b { "+" } else { "-" }.into()
}
fn valid_strand(s: &str) -> bool {
    s == "+" || s == "-"
}
pub(super) fn validate_axis(catalog: &Catalog, axis: &Axis) -> io::Result<Vec<usize>> {
    if axis.version != FORMAT_VERSION
        || axis.coordinate_system.trim().is_empty()
        || axis.intervals.is_empty()
    {
        return Err(invalid("invalid/empty reference axis"));
    }
    let groups: BTreeMap<_, _> = catalog
        .groups
        .iter()
        .enumerate()
        .map(|(g, x)| (x.id.as_str(), g))
        .collect();
    let mut seen_components = BTreeSet::new();
    let mut previous: Option<&AxisInterval> = None;
    let mut result = Vec::new();
    for a in &axis.intervals {
        let &g = groups
            .get(a.group.as_str())
            .ok_or_else(|| invalid("unknown axis group"))?;
        let reference = catalog
            .occurrences
            .get(a.reference_occurrence)
            .ok_or_else(|| invalid("unknown reference occurrence"))?;
        if a.component.trim().is_empty()
            || a.component.chars().any(char::is_control)
            || a.start >= a.end
            || a.end > i64::MAX as u64
            || reference.group != g
            || !valid_strand(&a.reference_strand)
        {
            return Err(invalid("invalid reference-axis interval/strand/membership"));
        }
        for (&id, strand) in &a.orientations {
            if !valid_strand(strand)
                || !catalog.groups[g].occurrences.contains(&id)
                || (id == a.reference_occurrence && strand != &a.reference_strand)
            {
                return Err(invalid("invalid axis orientation override"));
            }
        }
        if let Some(p) = previous {
            if p.component == a.component {
                if a.start <= p.start
                    || a.end <= p.end
                    || catalog.occurrences[p.reference_occurrence].source != reference.source
                {
                    return Err(invalid("axis component must have increasing coordinates and one exact reference source path"));
                }
            } else if seen_components.contains(&a.component) {
                return Err(invalid(
                    "axis components must be contiguous, not interleaved",
                ));
            }
        }
        seen_components.insert(a.component.clone());
        previous = Some(a);
        result.push(g);
    }
    Ok(result)
}

/// Both signed two-node contexts must agree throughout all shared nonpalindromic
/// evidence. No majority voting, contig-name homology or coordinate guessing.
fn orientations(
    catalog: &Catalog,
    axis: &AxisInterval,
    group: usize,
    frozen: Option<&BTreeMap<(usize, usize), u8>>,
) -> BTreeMap<usize, (Option<String>, String)> {
    let reference = &catalog.occurrences[axis.reference_occurrence];
    let mut support: BTreeMap<usize, u8> = BTreeMap::new();
    for &f in reference.feature_multiplicities.keys() {
        let feature = &catalog.features[f];
        let reference_nodes: BTreeSet<_> = feature
            .locations
            .iter()
            .filter(|l| l.containing_occurrences.contains(&reference.id))
            .map(|l| l.signed_nodes)
            .collect();
        for location in &feature.locations {
            let mut mask = 0;
            for &nodes in &reference_nodes {
                let rc = [-nodes[1], -nodes[0]];
                if nodes == rc {
                    continue;
                } // Palindrome supplies no orientation.
                if nodes == location.signed_nodes {
                    mask |= 1;
                }
                if rc == location.signed_nodes {
                    mask |= 2;
                }
            }
            if mask != 0 {
                for &id in &location.containing_occurrences {
                    if catalog.occurrences[id].group == group {
                        *support.entry(id).or_default() |= mask;
                    }
                }
            }
        }
    }
    if let Some(frozen) = frozen {
        for (&(_, id), &mask) in frozen.range((reference.id, 0)..=(reference.id, usize::MAX)) {
            *support.entry(id).or_default() |= mask;
        }
    }
    catalog.groups[group]
        .occurrences
        .iter()
        .map(|&id| {
            let occurrence = &catalog.occurrences[id];
            let value = if id == reference.id
                || (occurrence.source == reference.source
                    && occurrence.interval.start == reference.interval.start
                    && occurrence.interval.end == reference.interval.end)
            {
                (
                    Some(axis.reference_strand.clone()),
                    "reference-identity".into(),
                )
            } else if let Some(s) = axis.orientations.get(&id) {
                (Some(s.clone()), "explicit-axis-input".into())
            } else if let (Some(s), Some(r)) =
                (&occurrence.interval.strand, &reference.interval.strand)
            {
                (
                    Some(multiply_strands(
                        &multiply_strands(s, r),
                        &axis.reference_strand,
                    )),
                    "explicit-catalog-relative-to-reference".into(),
                )
            } else {
                match support.get(&id).copied().unwrap_or(0) {
                    1 => (
                        Some(axis.reference_strand.clone()),
                        "consistent-signed-contexts".into(),
                    ),
                    2 => (
                        Some(multiply_strands("-", &axis.reference_strand)),
                        "consistent-signed-contexts".into(),
                    ),
                    3 => (None, "conflicting-signed-contexts".into()),
                    _ => (None, "unsupported-or-palindromic-contexts".into()),
                }
            };
            (id, value)
        })
        .collect()
}

/// Exact path identity and strict monotonicity of BOTH interval endpoints.
/// Overlaps/gaps are legal and emitted explicitly; no sequence is fabricated.
fn continuation(a: &State, b: &State) -> bool {
    a.interval.path == b.interval.path
        && a.strand.is_some()
        && a.strand == b.strand
        && match a.strand.as_deref() {
            Some("+") => b.interval.start > a.interval.start && b.interval.end > a.interval.end,
            Some("-") => b.interval.start < a.interval.start && b.interval.end < a.interval.end,
            _ => false,
        }
}
fn transition(a: &State, b: &State, penalty: f64) -> f64 {
    if continuation(a, b) {
        0.0
    } else {
        penalty
    }
}

/// Exact forward/backward min-sum DP. All candidates survive; no beam/pruning.
/// Row minima are removed for numerical stability, then restored in the score.
fn solve(rows: &[ThreadInterval], penalty: f64) -> io::Result<(f64, Vec<Vec<usize>>)> {
    let minima: Vec<f64> = rows
        .iter()
        .map(|r| {
            r.states
                .iter()
                .map(|s| s.emission)
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    let emissions: Vec<Vec<f64>> = rows
        .iter()
        .zip(&minima)
        .map(|(r, min)| r.states.iter().map(|s| s.emission - min).collect())
        .collect();
    let mut forward = vec![emissions[0].clone()];
    for i in 1..rows.len() {
        let scores = rows[i]
            .states
            .iter()
            .enumerate()
            .map(|(b, state)| {
                emissions[i][b]
                    + rows[i - 1]
                        .states
                        .iter()
                        .enumerate()
                        .map(|(a, prev)| forward[i - 1][a] + transition(prev, state, penalty))
                        .fold(f64::INFINITY, f64::min)
            })
            .collect();
        forward.push(scores);
    }
    let mut backward: Vec<Vec<f64>> = rows.iter().map(|r| vec![0.0; r.states.len()]).collect();
    for i in (0..rows.len() - 1).rev() {
        for (a, state) in rows[i].states.iter().enumerate() {
            backward[i][a] = rows[i + 1]
                .states
                .iter()
                .enumerate()
                .map(|(b, next)| {
                    transition(state, next, penalty) + emissions[i + 1][b] + backward[i + 1][b]
                })
                .fold(f64::INFINITY, f64::min);
        }
    }
    let best = forward
        .last()
        .unwrap()
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let score = best + minima.iter().sum::<f64>();
    if !score.is_finite() {
        return Err(invalid("nonfinite threading score"));
    }
    let optimal = rows
        .iter()
        .enumerate()
        .map(|(i, row)| {
            (0..row.states.len())
                .filter(|&s| (forward[i][s] + backward[i][s] - best).abs() <= TIE_EPSILON)
                .collect()
        })
        .collect();
    Ok((score, optimal))
}
fn union_bp<'a>(rows: impl Iterator<Item = &'a AxisInterval>) -> io::Result<u64> {
    let mut end_by_component = BTreeMap::new();
    let mut total = 0u64;
    for a in rows {
        let end = end_by_component.entry(&a.component).or_insert(0);
        total = total
            .checked_add(a.end.saturating_sub(a.start.max(*end)))
            .ok_or_else(|| invalid("axis coverage overflow"))?;
        *end = (*end).max(a.end);
    }
    Ok(total)
}

pub fn thread(
    catalog: &Catalog,
    calls: &Genotypes,
    axis: Axis,
    switch_penalty: f64,
) -> io::Result<Threads> {
    if !switch_penalty.is_finite() || switch_penalty < 0.0 {
        return Err(invalid("switch penalty must be finite and nonnegative"));
    }
    if calls.panel != catalog.panel
        || calls.calls.len() != catalog.groups.len()
        || calls.ploidy != 1
    {
        return Err(invalid("incompatible quantitative calls"));
    }
    let groups = validate_axis(catalog, &axis)?;
    let mut frequencies = BTreeMap::new();
    for &g in &groups {
        *frequencies.entry(g).or_insert(0usize) += 1;
    }
    let frozen = if calls.model == super::observations::MODEL {
        let p = calls
            .observations
            .as_ref()
            .ok_or_else(|| invalid("missing observation provenance"))?;
        if !p.full_feature_scope
            || !p.feature_groups.is_empty()
            || calls.catalog_payload_checksum.is_some()
            || !p
                .orientation_support
                .windows(2)
                .all(|w| (w[0].0, w[0].1) < (w[1].0, w[1].1))
        {
            return Err(invalid(
                "partial/incompatible observation calls cannot be threaded",
            ));
        }
        let expected: BTreeSet<_> = axis
            .intervals
            .iter()
            .zip(&groups)
            .filter(|(_, g)| frequencies[g] == 1)
            .map(|(a, _)| a.reference_occurrence)
            .collect();
        if p.orientation_references != expected.iter().copied().collect::<Vec<_>>() {
            return Err(invalid("axis reference orientation provider was not assessed; derive axis-specific calls first"));
        }
        let mut masks = BTreeMap::new();
        for &(a, b, mask) in &p.orientation_support {
            if a >= catalog.occurrences.len()
                || b >= catalog.occurrences.len()
                || !expected.contains(&a)
                || mask == 0
                || mask > 3
                || catalog.occurrences[a].group != catalog.occurrences[b].group
                || masks.insert((a, b), mask).is_some()
            {
                return Err(invalid("invalid frozen signed orientation evidence"));
            }
        }
        Some(masks)
    } else if calls.model == super::genotype::MODEL && calls.observations.is_none() {
        None
    } else {
        return Err(invalid("unknown quantitative threading model"));
    };
    let repeated_axis_groups = frequencies
        .iter()
        .filter(|(_, n)| **n > 1)
        .map(|(&g, _)| catalog.groups[g].id.clone())
        .collect();
    let unscaffolded_groups = catalog
        .groups
        .iter()
        .enumerate()
        .filter(|(g, _)| !frequencies.contains_key(g))
        .map(|(_, g)| g.id.clone())
        .collect();
    let mut intervals = Vec::new();
    for (a, g) in axis.intervals.into_iter().zip(groups) {
        let call = &calls.calls[g];
        if call.group != a.group {
            return Err(invalid("call/group ordering mismatch"));
        }
        let mut row = ThreadInterval {
            axis: a,
            status: String::new(),
            states: Vec::new(),
            optimal_states: Vec::new(),
            segment: None,
        };
        if frequencies[&g] > 1 {
            row.status = "unresolved-repeated-axis-group".into();
        } else if (calls.model == super::observations::MODEL
            && !matches!(call.status.as_str(), "informative" | "tied"))
            || (calls.model != super::observations::MODEL
                && (call.status.starts_with("no-call") || call.status == "poor-fit"))
        {
            row.status = format!("unresolved-{}", call.status);
        } else {
            let orientations = orientations(catalog, &row.axis, g, frozen.as_ref());
            let mut physical: BTreeMap<_, usize> = BTreeMap::new();
            for (b, bundle) in call.bundles.iter().enumerate() {
                if calls.model == super::observations::MODEL
                    && !super::observations::eligible(call, bundle)
                {
                    continue;
                }
                for &id in &bundle.source_occurrences {
                    let occurrence = &catalog.occurrences[id];
                    let (strand, evidence) = &orientations[&id];
                    let key = (
                        occurrence.source,
                        occurrence.interval.start,
                        occurrence.interval.end,
                        strand.clone(),
                    );
                    if let Some(&s) = physical.get(&key) {
                        row.states[s].source_occurrences.push(id);
                    } else {
                        physical.insert(key, row.states.len());
                        row.states.push(State {
                            bundle: b,
                            identity: bundle.identity.clone(),
                            source_occurrences: vec![id],
                            interval: occurrence.interval.clone(),
                            strand: strand.clone(),
                            orientation_evidence: evidence.clone(),
                            emission: bundle.score,
                            mean_deviance: bundle.mean_deviance,
                        });
                    }
                }
            }
        }
        intervals.push(row);
    }
    let mut segments = Vec::new();
    let mut start = 0;
    while start < intervals.len() {
        if intervals[start].states.is_empty() {
            start += 1;
            continue;
        }
        let mut end = start + 1;
        while end < intervals.len()
            && !intervals[end].states.is_empty()
            && intervals[end].axis.component == intervals[start].axis.component
        {
            end += 1;
        }
        let (score, optimal) = solve(&intervals[start..end], switch_penalty)?;
        let tied = optimal.iter().any(|v| v.len() > 1);
        for (row, states) in intervals[start..end].iter_mut().zip(optimal) {
            row.segment = Some(segments.len());
            row.status = if states.len() != 1 {
                "unresolved-tied-paths"
            } else if row.states[states[0]].strand.is_none() {
                "unresolved-orientation"
            } else if row.states[states[0]].mean_deviance > calls.parameters.max_mean_deviance {
                "unresolved-selected-poor-fit"
            } else {
                "resolved-experimental"
            }
            .into();
            row.optimal_states = states;
        }
        segments.push(Segment {
            first_interval: start,
            last_interval: end - 1,
            score,
            tied,
        });
        start = end;
    }
    let mut blocks: Vec<Block> = Vec::new();
    let mut switches = Vec::new();
    let mut previous: Option<usize> = None;
    for (i, row) in intervals.iter().enumerate() {
        if row.status != "resolved-experimental" {
            previous = None;
            continue;
        }
        let state = &row.states[row.optimal_states[0]];
        let mut extend = false;
        if let Some(p) = previous {
            let prev_row = &intervals[p];
            let prev = &prev_row.states[prev_row.optimal_states[0]];
            if prev_row.axis.component == row.axis.component && prev_row.segment == row.segment {
                let reference_gap = row.axis.start as i64 - prev_row.axis.end as i64;
                if continuation(prev, state) {
                    let block = blocks.last_mut().unwrap();
                    block.interval_indices.push(i);
                    block.reference_gaps_bp.push(reference_gap);
                    block
                        .source_gaps_bp
                        .push(if state.strand.as_deref() == Some("+") {
                            state.interval.start as i64 - prev.interval.end as i64
                        } else {
                            prev.interval.start as i64 - state.interval.end as i64
                        });
                    extend = true;
                } else {
                    switches.push(Switch {
                        component: row.axis.component.clone(),
                        left_interval: p,
                        right_interval: i,
                        from_path: prev.interval.path.clone(),
                        to_path: state.interval.path.clone(),
                        reason: if prev.interval.path != state.interval.path {
                            "source-path-change"
                        } else if prev.strand != state.strand {
                            "orientation-change-block-break-not-certified-inversion"
                        } else {
                            "nonmonotonic-source-block-break"
                        }
                        .into(),
                        breakpoint_start: prev_row.axis.start,
                        breakpoint_end: row.axis.end,
                        reference_gap_bp: reference_gap,
                    });
                }
            }
        }
        if !extend {
            blocks.push(Block {
                component: row.axis.component.clone(),
                source_path: state.interval.path.clone(),
                strand: state.strand.clone().unwrap(),
                interval_indices: vec![i],
                source_gaps_bp: Vec::new(),
                reference_gaps_bp: Vec::new(),
            });
        }
        previous = Some(i);
    }
    let resolved_intervals = intervals
        .iter()
        .filter(|r| r.status == "resolved-experimental")
        .count();
    Ok(Threads { version: FORMAT_VERSION, model: if calls.model == super::observations::MODEL { super::observations::THREAD_MODEL } else { "haploid-source-occurrence-min-sum-v1" }.into(),
        panel: calls.panel.clone(), count_policy: calls.count_policy.clone(), sample_payload_checksum: calls.sample_payload_checksum.clone(), catalog_payload_checksum: calls.catalog_payload_checksum.clone(), observation_metadata_checksum: calls.observations.as_ref().map(|p| p.metadata_checksum.clone()),
        experimental: true, catalog_accepted: false,
        coordinate_system: axis.coordinate_system, switch_penalty, tie_epsilon: TIE_EPSILON,
        interpretation: "Count-score emissions plus explicit continuation/switch cost, chromosome resets; no posterior or read-linkage phase. Optimal states are marginal alternatives, NOT freely combinable paths. Unknown orientations never establish continuity. Repeated-axis groups excluded entirely; accessory calls retained. Breakpoints span flanking partitions, not exact edges. Blocks preserve source/reference gaps and overlaps without sequence emission.".into(),
        axis_union_bp: union_bp(intervals.iter().map(|r| &r.axis))?, resolved_union_bp: union_bp(intervals.iter().filter(|r| r.status == "resolved-experimental").map(|r| &r.axis))?,
        unresolved_intervals: intervals.len() - resolved_intervals, resolved_intervals, intervals, segments, blocks, switches, unscaffolded_groups, repeated_axis_groups })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn state(path: &str, start: u64, end: u64, strand: Option<&str>, emission: f64) -> State {
        State {
            bundle: 0,
            identity: path.into(),
            source_occurrences: vec![0],
            interval: Interval {
                path: path.into(),
                start,
                end,
                strand: None,
            },
            strand: strand.map(str::to_string),
            orientation_evidence: "fixture".into(),
            emission,
            mean_deviance: 0.0,
        }
    }
    fn row(states: Vec<State>) -> ThreadInterval {
        ThreadInterval {
            axis: AxisInterval {
                component: "reference#0#chr1".into(),
                start: 0,
                end: 10,
                group: "g".into(),
                reference_occurrence: 0,
                reference_strand: "+".into(),
                orientations: BTreeMap::new(),
            },
            status: String::new(),
            states,
            optimal_states: Vec::new(),
            segment: None,
        }
    }
    #[test]
    fn signed_context_orientation_requires_consistency_not_majority() {
        use super::super::{catalog::*, PanelIdentity};
        use crate::sample_mem_bwt::{canonical, encode_walk};
        let sources: Vec<_> = (0..6)
            .map(|id| Source {
                id,
                path: format!("S{id}#0#full-contig"),
                length: 100,
            })
            .collect();
        let occurrences = (0..6)
            .map(|id| Occurrence {
                id,
                group: 0,
                source: id,
                interval: Interval {
                    path: sources[id].path.clone(),
                    start: 0,
                    end: 100,
                    strand: None,
                },
                fully_contained_anchors: 2,
                anchor_status: "fixture".into(),
                feature_multiplicities: if id == 0 {
                    BTreeMap::from([(0, 1), (1, 1)])
                } else {
                    BTreeMap::new()
                },
            })
            .collect();
        let location = |source, start, signed_nodes| FeatureLocation {
            source,
            start,
            end: start + 10,
            signed_nodes,
            containing_occurrences: vec![source],
        };
        let catalog = Catalog {
            version: FORMAT_VERSION,
            panel: PanelIdentity {
                checksum_algorithm: "fixture".into(),
                sidecars: Vec::new(),
            },
            catalog_accepted: false,
            validation: "fixture".into(),
            sources,
            occurrences,
            groups: vec![Group {
                id: "g".into(),
                scaffold: None,
                occurrences: (0..6).collect(),
            }],
            links: Vec::new(),
            features: vec![
                Feature {
                    tokens: canonical(&encode_walk(&[(1, 0), (2, 5)]).unwrap()),
                    owning_group: Some(0),
                    exclusion_reasons: Vec::new(),
                    locations: vec![
                        location(0, 10, [1, 2]),
                        location(1, 10, [1, 2]),
                        location(2, 10, [-2, -1]),
                        location(3, 10, [1, 2]),
                        location(3, 30, [1, 2]),
                        location(3, 50, [-2, -1]),
                    ],
                },
                Feature {
                    tokens: canonical(&encode_walk(&[(3, 0), (-3, 5)]).unwrap()),
                    owning_group: Some(0),
                    exclusion_reasons: Vec::new(),
                    locations: vec![location(0, 60, [3, -3]), location(4, 10, [3, -3])],
                },
            ],
        };
        let mut axis = row(Vec::new()).axis;
        let inferred = orientations(&catalog, &axis, 0, None);
        assert_eq!(inferred[&0].0.as_deref(), Some("+"));
        assert_eq!(inferred[&1].0.as_deref(), Some("+"));
        assert_eq!(inferred[&2].0.as_deref(), Some("-"));
        assert_eq!(inferred[&3].0, None); // two forward vs one reverse is NOT a vote.
        assert_eq!(inferred[&3].1, "conflicting-signed-contexts");
        assert_eq!(inferred[&4].0, None);
        assert_eq!(inferred[&5].0, None);
        axis.reference_strand = "-".into();
        assert_eq!(
            orientations(&catalog, &axis, 0, None)[&2].0.as_deref(),
            Some("+")
        );
        axis.orientations.insert(3, "-".into());
        assert_eq!(
            orientations(&catalog, &axis, 0, None)[&3].1,
            "explicit-axis-input"
        );
    }
    #[test]
    fn exact_dp_matches_exhaustive_paths_with_reverse_ties_gaps_and_path_changes() {
        for penalty in [0.0, 1.0, 5.0, 100.0] {
            let rows = vec![
                row(vec![
                    state("A#0#chr1", 10, 20, Some("+"), 0.0),
                    state("B#0#chrX", 90, 100, Some("-"), 1.0),
                ]),
                row(vec![
                    state("A#0#chr1", 18, 30, Some("+"), 2.0),
                    state("B#0#chrX", 70, 85, Some("-"), 0.0),
                    state("B#0#copy2", 50, 60, Some("-"), 0.0),
                ]),
                row(vec![
                    state("A#0#chr1", 31, 40, Some("+"), 0.0),
                    state("B#0#chrX", 50, 60, Some("-"), 1.0),
                    state("U#0#contig", 1, 10, None, 0.0),
                ]),
            ];
            let (score, optimal) = solve(&rows, penalty).unwrap();
            let mut exhaustive = Vec::new();
            for a in 0..rows[0].states.len() {
                for b in 0..rows[1].states.len() {
                    for c in 0..rows[2].states.len() {
                        let path = [a, b, c];
                        let s = (0..3)
                            .map(|i| rows[i].states[path[i]].emission)
                            .sum::<f64>()
                            + (1..3)
                                .map(|i| {
                                    transition(
                                        &rows[i - 1].states[path[i - 1]],
                                        &rows[i].states[path[i]],
                                        penalty,
                                    )
                                })
                                .sum::<f64>();
                        exhaustive.push((s, path));
                    }
                }
            }
            let best = exhaustive
                .iter()
                .map(|(s, _)| *s)
                .fold(f64::INFINITY, f64::min);
            assert!((score - best).abs() < 1e-10);
            for i in 0..3 {
                let expected: BTreeSet<_> = exhaustive
                    .iter()
                    .filter(|(s, _)| (s - best).abs() <= TIE_EPSILON)
                    .map(|(_, p)| p[i])
                    .collect();
                assert_eq!(
                    optimal[i].iter().copied().collect::<BTreeSet<_>>(),
                    expected
                );
            }
        }
        assert!(!continuation(
            &state("A", 10, 30, Some("+"), 0.0),
            &state("A", 11, 20, Some("+"), 0.0)
        ));
        assert!(!continuation(
            &state("A", 10, 20, None, 0.0),
            &state("A", 21, 30, None, 0.0)
        ));
        assert!(!continuation(
            &state("A", 10, 20, Some("+"), 0.0),
            &state("A", 1, 9, Some("+"), 0.0)
        ));
    }
}
