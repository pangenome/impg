//! Quarantined coverage-vector diagnostics (IMPG_COSINE_DIAG_OUTPUT and
//! IMPG_COSINE_EXHAUSTIVE_OUTPUT). Saved reads and product scores are not
//! modified; no cosine score feeds a sweep, DP, posterior, or product column.
use super::*;
use std::io::{self, BufWriter, Write};

/// Public panel material only: source BED/completion rows and adjacent
/// same-source stitches. Port-cut, attested-junction, and phasing-partial
/// constructs stay in their legacy lanes; they are not pure-route alleles.
pub(super) fn physical_material_row(row: &genome::SpanningTraversal) -> bool {
    if row.identity.starts_with("single:") {
        return row.segments.len() == 1;
    }
    row.identity.starts_with("stitch:")
        && row.segments.len() > 1
        && row.segments.windows(2).all(|pair| {
            let (left, right) = (&pair[0], &pair[1]);
            left.source == right.source
                && left.reverse == right.reverse
                && if left.reverse { left.start == right.end } else { left.end == right.start }
        })
}

fn merge_intervals(intervals: &mut Vec<(u64, u64)>) {
    intervals.sort_unstable();
    let mut end = 0;
    for index in 0..intervals.len() {
        let (lo, hi) = intervals[index];
        if end > 0 && intervals[end - 1].1 >= lo {
            intervals[end - 1].1 = intervals[end - 1].1.max(hi);
        } else {
            intervals[end] = (lo, hi);
            end += 1;
        }
    }
    intervals.truncate(end);
}

/// Material bins are disjoint path-coordinate intervals. Within one record,
/// overlapping MEM-subwalk placement spans cover each base ONCE, not once per
/// subwalk; the record then contributes its routing-established m/t share.
fn material_coverage_bins(
    alleles: &BTreeSet<usize>,
    candidates: &[genome::SpanningTraversal],
    path_of_source: &[usize],
    instances: &crate::InstanceStructure,
    locus: usize,
) -> (Vec<(usize, u64, u64, f64)>, usize) {
    let mut material: BTreeMap<usize, Vec<(u64, u64)>> = BTreeMap::new();
    for &allele in alleles {
        for segment in &candidates[allele].segments {
            if segment.start < segment.end {
                material.entry(path_of_source[segment.source])
                    .or_default().push((segment.start, segment.end));
            }
        }
    }
    for intervals in material.values_mut() { merge_intervals(intervals); }
    let mut records: Vec<u32> = instances.window_records[locus].values()
        .flat_map(|ids| ids.iter().copied()).collect();
    records.sort_unstable();
    records.dedup();
    let mut events: BTreeMap<usize, BTreeMap<u64, f64>> = BTreeMap::new();
    for &record in &records {
        let share = instances.record_shares[record as usize];
        let mut placed: BTreeMap<usize, Vec<(u64, u64)>> = BTreeMap::new();
        for spans in instances.record_spans[record as usize].values() {
            for &(path, lo, hi) in spans {
                if let Some(material_intervals) = material.get(&path) {
                    for &(start, end) in material_intervals {
                        if lo < end && start < hi {
                            placed.entry(path).or_default().push((lo.max(start), hi.min(end)));
                        }
                    }
                }
            }
        }
        for (path, mut intervals) in placed {
            merge_intervals(&mut intervals);
            let path_events = events.entry(path).or_default();
            for (lo, hi) in intervals {
                *path_events.entry(lo).or_default() += share;
                *path_events.entry(hi).or_default() -= share;
            }
        }
    }
    let mut bins = Vec::new();
    for (path, intervals) in material {
        for (start, end) in intervals {
            let mut at = start;
            let mut coverage = 0.0;
            if let Some(path_events) = events.get(&path) {
                for (&position, &delta) in path_events.range(start..=end) {
                    if position > at {
                        bins.push((path, at, position, coverage));
                    }
                    coverage += delta;
                    at = position;
                }
            }
            if at < end { bins.push((path, at, end, coverage)); }
        }
    }
    (bins, records.len())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn dump(
    path: &str,
    panel: &SyngIndex,
    sources: &routes::Sources,
    path_of_source: &[usize],
    ranges: &[Vec<genome::SpanningTraversal>],
    classing: &[LocusClassing],
    sweeps: &[LocusSweep],
    viable: &[Vec<bool>],
    truth_pieces: &[[Vec<(usize, u64, u64, bool, u32)>; 2]],
    selected_pairs: &[[usize; 2]],
    window_obs: &[HashMap<FeatureKey, f64>],
    instances: &crate::InstanceStructure,
    depth: f64,
    locus_offset: usize,
) -> io::Result<()> {
    let mut output = BufWriter::new(std::fs::File::create(path)?);
    for locus in 0..ranges.len() {
        let alleles = &ranges[locus];
        let truth = [0usize, 1].map(|copy| {
            let pieces = &truth_pieces[locus][copy];
            (!pieces.is_empty()).then(|| alleles.iter().position(|candidate| {
                candidate.segments.len() == pieces.len()
                    && candidate.segments.iter().zip(pieces).all(|(segment, piece)| {
                        (segment.source, segment.start, segment.end, segment.reverse)
                            == (piece.0, piece.1, piece.2, piece.3)
                    })
            })).flatten()
        });
        let (Some(first), Some(second)) = (truth[0], truth[1]) else {
            serde_json::to_writer(&mut output, &serde_json::json!({
                "locus": locus + locus_offset, "truth_pair_expressible": false
            }))?;
            writeln!(output)?;
            continue;
        };
        let mut representative = vec![None; classing[locus].profiles.len()];
        let mut representative_viable = vec![None; classing[locus].profiles.len()];
        for (allele, &class) in classing[locus].membership.iter().enumerate() {
            representative[class].get_or_insert(allele);
            if viable[locus][allele] {
                representative_viable[class].get_or_insert(allele);
            }
        }
        let mut top: Vec<([usize; 2], f64)> = Vec::new();
        let mut top_viable: Vec<([usize; 2], f64)> = Vec::new();
        for right in 0..representative.len() {
            for left in 0..=right {
                let loss = sweeps[locus].table[class_pair_index(left, right)];
                if !loss.is_finite() { continue; }
                let pair = [left, right];
                if top.len() < 8 || loss < top.last().unwrap().1 {
                    top.push((pair, loss));
                    top.sort_by(|a,b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
                    top.truncate(8);
                }
                if representative_viable[left].is_some() && representative_viable[right].is_some()
                    && (top_viable.len() < 8 || loss < top_viable.last().unwrap().1) {
                    top_viable.push((pair, loss));
                    top_viable.sort_by(|a,b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
                    top_viable.truncate(8);
                }
            }
        }
        let mut pairs = vec![("truth".to_owned(), [first, second], None),
                             ("selected".to_owned(), selected_pairs[locus], None)];
        for (rank, (classes, loss)) in top.iter().enumerate() {
            pairs.push((format!("old_objective_top_{}", rank + 1),
                        [representative[classes[0]].unwrap(), representative[classes[1]].unwrap()],
                        Some(*loss)));
        }
        for (rank, (classes, loss)) in top_viable.iter().enumerate() {
            pairs.push((format!("old_viable_top_{}", rank + 1),
                        [representative_viable[classes[0]].unwrap(), representative_viable[classes[1]].unwrap()],
                        Some(*loss)));
        }
        let used: BTreeSet<usize> = pairs.iter().flat_map(|(_, pair, _)| *pair).collect();
        let (bins, records) = material_coverage_bins(&used, alleles, path_of_source, instances, locus);
        let mut oracle_memo = HashMap::new();
        let mut expected = Vec::new();
        for &allele in &used {
            let profile = crate::oracle_allele_profile(panel, sources, &alleles[allele], &mut oracle_memo)?;
            expected.push(serde_json::json!({
                "allele": allele,
                "material_intervals": alleles[allele].segments.iter()
                    .map(|segment| (path_of_source[segment.source], segment.start, segment.end))
                    .collect::<Vec<_>>(),
                "feature_profile": profile.iter().map(|(key, &count)| (key, count))
                    .collect::<Vec<_>>(),
            }));
        }
        let mut observed: Vec<_> = window_obs[locus].iter().collect();
        observed.sort_by(|a,b| a.0.cmp(b.0));
        serde_json::to_writer(&mut output, &serde_json::json!({
            "locus": locus + locus_offset,
            "truth_pair_expressible": true,
            "truth_alleles_port_viable": viable[locus][first] && viable[locus][second],
            "truth_classes_port_viable": representative_viable[classing[locus].membership[first]].is_some()
                && representative_viable[classing[locus].membership[second]].is_some(),
            "depth_per_copy": depth,
            "record_count": records,
            "coverage_bins": bins,
            "observed_feature_mass": observed,
            "alleles": expected,
            "candidate_pairs": pairs.iter().map(|(label, pair, loss)| serde_json::json!({
                "label": label, "alleles": pair, "old_objective_loss": loss
            })).collect::<Vec<_>>(),
        }))?;
        writeln!(output)?;
    }
    output.flush()
}

#[cfg(test)]
mod tests {
    use super::{material_coverage_bins, material_overlap, merge_intervals, physical_material_row, CoverageSpace};
    use crate::{genome, FeatureKey, InstanceStructure, SourceRange};
    use std::collections::{BTreeSet, HashMap};

    #[test]
    fn cosine_probe_record_share_covers_union_once() {
        let mut record_spans = HashMap::new();
        record_spans.insert(vec![1u64] as FeatureKey, vec![(0usize, 0u64, 5u64)]);
        record_spans.insert(vec![2u64] as FeatureKey, vec![(0usize, 3u64, 10u64)]);
        let instances = InstanceStructure {
            record_spans: vec![record_spans],
            record_shares: vec![2.0],
            window_records: vec![HashMap::from([(vec![1], vec![0]), (vec![2], vec![0])])],
        };
        let ranges = [genome::SpanningTraversal {
            partition: 0,
            identity: "sample".into(),
            segments: vec![SourceRange {
                partition: 0, occurrence: 0, source: 0, start: 0, end: 10, reverse: false,
            }],
        }];
        let (bins, records) = material_coverage_bins(
            &BTreeSet::from([0]), &ranges, &[0], &instances, 0,
        );
        assert_eq!(records, 1);
        assert_eq!(bins, vec![(0, 0, 10, 2.0)]);
    }

    #[test]
    fn physical_material_row_rejects_constructed_port_and_partial_rows() {
        let base = SourceRange {
            partition: 0, occurrence: 0, source: 3, start: 0, end: 10, reverse: false,
        };
        let single = genome::SpanningTraversal::single(base.clone());
        assert!(physical_material_row(&single));
        let stitched = genome::SpanningTraversal {
            partition: 0, identity: "stitch:0:3:0-20:2".into(),
            segments: vec![base.clone(), SourceRange { start: 10, end: 20, ..base.clone() }],
        };
        assert!(physical_material_row(&stitched));
        assert!(!physical_material_row(&genome::SpanningTraversal {
            identity: "split:0:3".into(), ..stitched.clone()
        }));
        assert!(!physical_material_row(&genome::SpanningTraversal {
            identity: "cooc-partial:0:3".into(), ..single.clone()
        }));
        assert!(!physical_material_row(&genome::SpanningTraversal {
            segments: vec![base.clone(), SourceRange { source: 4, ..base }],
            ..stitched
        }));
    }

    #[test]
    fn cosine_factored_pairs_match_direct_base_vectors_and_collapse_duplicates() {
        let row = |id: &str, lo, hi| genome::SpanningTraversal {
            partition: 0, identity: id.into(),
            segments: vec![SourceRange { partition: 0, occurrence: 0, source: 0,
                start: lo, end: hi, reverse: false }],
        };
        let mut spans = HashMap::new();
        spans.insert(vec![1u64] as FeatureKey, vec![(0, 0, 5)]);
        let instances = InstanceStructure {
            record_spans: vec![spans], record_shares: vec![2.0],
            window_records: vec![HashMap::from([(vec![1], vec![0])])],
        };
        let rows = vec![row("a", 0, 5), row("b", 3, 8), row("a_duplicate", 0, 5)];
        let space = CoverageSpace::new(&rows, &[0], &instances, 0, 15.0);
        assert_eq!(space.rows.len(), 2);
        assert_eq!(material_overlap(&space.rows[0].material, &space.rows[1].material), 2.0);
        // Observation [2,2,2,2,2,0,0,0], pair expectation
        // [15,15,15,30,30,15,15,15] in source coordinates.
        let observed = [2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0];
        let expected = [15.0, 15.0, 15.0, 30.0, 30.0, 15.0, 15.0, 15.0];
        let dot = |a: &[f64; 8], b: &[f64; 8]| a.iter().zip(b).map(|(x,y)| x*y).sum::<f64>();
        let direct = dot(&observed, &expected)
            / (dot(&observed, &observed) * dot(&expected, &expected)).sqrt();
        assert!((space.cosine(0,1,15.0).unwrap() - direct).abs() < 1e-14);
    }

    #[test]
    fn cosine_probe_record_placement_intervals_do_not_double_count_subwalks() {
        let mut spans = vec![(11, 14), (10, 12), (20, 22), (14, 16)];
        merge_intervals(&mut spans);
        assert_eq!(spans, vec![(10, 16), (20, 22)]);
    }
}

// Exhaustive, assessment-side pilot. The physical coverage equivalence is
// EXACT interval equality, never the old MEM-profile/Poisson class equality.
type Material = Vec<(usize, u64, u64)>;

fn row_material(row: &genome::SpanningTraversal, path_of_source: &[usize]) -> Material {
    let mut intervals: BTreeMap<usize, Vec<(u64, u64)>> = BTreeMap::new();
    for segment in &row.segments {
        if segment.start < segment.end {
            intervals.entry(path_of_source[segment.source])
                .or_default().push((segment.start, segment.end));
        }
    }
    let mut result = Vec::new();
    for (path, mut ranges) in intervals {
        merge_intervals(&mut ranges);
        result.extend(ranges.into_iter().map(|(lo, hi)| (path, lo, hi)));
    }
    result
}

fn material_overlap(left: &Material, right: &Material) -> f64 {
    let (mut i, mut j, mut total) = (0, 0, 0u64);
    while i < left.len() && j < right.len() {
        let (path_a, lo_a, hi_a) = left[i];
        let (path_b, lo_b, hi_b) = right[j];
        if path_a == path_b { total += hi_a.min(hi_b).saturating_sub(lo_a.max(lo_b)); }
        if (path_a, hi_a) < (path_b, hi_b) { i += 1; } else { j += 1; }
    }
    total as f64
}

#[derive(Debug)]
struct CoverageRow {
    material: Material,
    members: Vec<usize>,
    observed_dot: f64,
    expected_norm: f64,
}

struct CoverageSpace {
    rows: Vec<CoverageRow>,
    observed_norm: f64,
}

impl CoverageSpace {
    fn new(
        candidates: &[genome::SpanningTraversal],
        path_of_source: &[usize],
        instances: &crate::InstanceStructure,
        locus: usize,
        depth: f64,
    ) -> Self {
        let mut groups: BTreeMap<Material, Vec<usize>> = BTreeMap::new();
        for (index, row) in candidates.iter().enumerate() {
            groups.entry(row_material(row, path_of_source)).or_default().push(index);
        }
        let all: BTreeSet<usize> = (0..candidates.len()).collect();
        let (bins, _) = material_coverage_bins(&all, candidates, path_of_source, instances, locus);
        let observed_norm = bins.iter()
            .map(|(_, lo, hi, observed)| (*hi - *lo) as f64 * observed * observed)
            .sum();
        let rows = groups.into_iter().map(|(material, members)| {
            let observed_dot = bins.iter().map(|&(path, lo, hi, observed)| {
                material.iter()
                    .filter(|&&(p, start, end)| p == path && start < hi && lo < end)
                    .map(|&(_, start, end)|
                        (hi.min(end) - lo.max(start)) as f64 * observed * depth)
                    .sum::<f64>()
            }).sum();
            let expected_norm = material.iter()
                .map(|&(_, lo, hi)| (hi - lo) as f64 * depth * depth)
                .sum();
            CoverageRow { material, members, observed_dot, expected_norm }
        }).collect();
        Self { rows, observed_norm }
    }

    fn cosine(&self, first: usize, second: usize, depth: f64) -> Option<f64> {
        let a = &self.rows[first];
        let b = &self.rows[second];
        let expected_norm = a.expected_norm + b.expected_norm
            + 2.0 * depth * depth * material_overlap(&a.material, &b.material);
        (self.observed_norm > 0.0 && expected_norm > 0.0)
            .then(|| (a.observed_dot + b.observed_dot)
                / (self.observed_norm * expected_norm).sqrt())
    }
}

fn truth_index(
    candidates: &[genome::SpanningTraversal],
    pieces: &[(usize, u64, u64, bool, u32)],
) -> Option<usize> {
    (!pieces.is_empty()).then(|| candidates.iter().position(|row| {
        row.segments.len() == pieces.len()
            && row.segments.iter().zip(pieces).all(|(segment, piece)|
                (segment.source, segment.start, segment.end, segment.reverse)
                    == (piece.0, piece.1, piece.2, piece.3))
    })).flatten()
}

/// Stream the entire diversity-bounded candidate-vector universe. A truth
/// pair is used only to ASSESS this separate receipt, not to construct rows,
/// observed mass, candidate scores, or production states.
pub(super) fn dump_exhaustive(
    path: &str,
    ranges: &[Vec<genome::SpanningTraversal>],
    path_of_source: &[usize],
    truth_pieces: &[[Vec<(usize, u64, u64, bool, u32)>; 2]],
    instances: &crate::InstanceStructure,
    depth: f64,
    locus_offset: usize,
) -> io::Result<()> {
    let mut report = BufWriter::new(std::fs::File::create(path)?);
    let mut competitors = BufWriter::new(std::fs::File::create(format!("{path}.competitors.jsonl"))?);
    for (locus, candidates) in ranges.iter().enumerate() {
        let truth = [0, 1].map(|copy| truth_index(candidates, &truth_pieces[locus][copy]));
        let space = CoverageSpace::new(candidates, path_of_source, instances, locus, depth);
        let mut id_to_row = vec![0usize; candidates.len()];
        for (index, row) in space.rows.iter().enumerate() {
            for &id in &row.members { id_to_row[id] = index; }
        }
        let pair_truth = match (truth[0], truth[1]) {
            (Some(a), Some(b)) => Some((id_to_row[a], id_to_row[b])),
            _ => None,
        };
        let truth_score = pair_truth.and_then(|(a,b)| space.cosine(a,b,depth));
        let mut better = 0u64;
        let mut ties = 0u64;
        let mut scored = 0u64;
        let mut best: Option<(f64, [usize;2])> = None;
        for second in 0..space.rows.len() {
            for first in 0..=second {
                let Some(score) = space.cosine(first, second, depth) else { continue; };
                scored += 1;
                if best.is_none_or(|(old, _)| score > old) { best = Some((score, [first, second])); }
                if let Some(truth_value) = truth_score {
                    if score > truth_value + 1e-12 {
                        better += 1;
                        let ids = [space.rows[first].members[0], space.rows[second].members[0]];
                        serde_json::to_writer(&mut competitors, &serde_json::json!({
                            "locus": locus + locus_offset,
                            "row_indices": ids,
                            "identities": [candidates[ids[0]].identity, candidates[ids[1]].identity],
                            "material": [&space.rows[first].material, &space.rows[second].material],
                            "cosine": score,
                            "truth_cosine": truth_value,
                        }))?;
                        writeln!(competitors)?;
                    } else if (score - truth_value).abs() <= 1e-12 { ties += 1; }
                }
            }
        }
        serde_json::to_writer(&mut report, &serde_json::json!({
            "locus": locus + locus_offset,
            "physical_rows": candidates.len(),
            "material_vectors": space.rows.len(),
            "eligible_pairs": scored,
            "truth_piece_presence": truth_pieces[locus].iter().map(|v| !v.is_empty()).collect::<Vec<_>>(),
            "truth_pair_expressible": pair_truth.is_some(),
            "truth_rows": truth,
            "truth_material_rows": pair_truth.map(|(a,b)| [a,b]),
            "truth_cosine": truth_score,
            "truth_rank": truth_score.map(|_| better + 1),
            "truth_tied_pairs": truth_score.map(|_| ties),
            "higher_cosine_competitors": better,
            "best_cosine": best.map(|(score, _)| score),
            "best_row_indices": best.map(|(_, indices)| indices
                .map(|index| space.rows[index].members[0])),
        }))?;
        writeln!(report)?;
    }
    report.flush()?;
    competitors.flush()
}

// ---------------------------------------------------------------------------
// Graph-space condensation (owner greenlight 2026-10-01). Both observed
// read evidence and per-candidate expected usage live on the panel graph's
// OWN shared coordinates: nodes are the syng syncmer segments (signed node
// ids, 1-based, sign = storage strand; identity frame-independent, shared
// across paths wherever sequence is conserved), and edges are the
// adjacencies between consecutive covered nodes along a placement path —
// so two conserved placements collapse onto the SAME node/edge features
// while diverged material keeps its own. Observed node mass: each
// contributing routed record's equal share m/t, once per graph node whose
// k-mer window is fully contained in the union of the record's routed
// placement spans (all features, both frames, merged per path); a record
// placing on several conserved paths votes a shared node ONCE. Observed
// edge mass: the adjacencies the records span — consecutive contained
// steps of one placement — same share, once per record per edge. A
// junction adjacency no single record spans would need the span index's
// chain-pair channel; the pure-route domain cannot contain one (admissible
// rows are `single:` or strictly adjacent same-source `stitch:`, so every
// spelled adjacency is same-path continuous), and the audited per-locus
// `multi_segment_rows`/`disjoint_material_rows` counts report that
// structure instead of assuming it. Expected usage: depth per covered
// copy with repeat-visit multiplicity; cosine factored per row exactly as
// the row-space pilot (`cosine_probe.rs` factorization), in TWO separate
// arms — nodes-only and nodes+edges — so the edge layer's contribution is
// a measured number. Assessment-side only; no thresholds, no constants.
// ---------------------------------------------------------------------------

/// A shared graph segment key: the syncmer node id, unsigned (the sign is
/// the storage orientation; the physical segment identity is the abs id).
type GraphNode = u32;
/// An unpacked shared adjacency: two consecutive covered node ids.
type GraphEdge = (u32, u32);
/// A shared adjacency key, packed for the O(rows^2) pair walk.
type PackedEdge = u64;

fn pack_edge(left: GraphNode, right: GraphNode) -> PackedEdge {
    (left as u64) << 32 | right as u64
}

/// One candidate row's graph usage: sorted distinct keys with repeat-visit
/// multiplicity (a row whose spelled material visits the same shared
/// segment/adjacency twice expects 2x depth there).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct GraphUsage {
    nodes: Vec<(GraphNode, u32)>,
    edges: Vec<(PackedEdge, u32)>,
}

/// The window-contained step subrange of `steps` (sorted by bp): steps
/// whose full k-mer window [bp, bp + k) lies inside [lo, hi). Containment,
/// not overlap, is the shared-segment coverage convention on both the
/// observed and expected sides.
fn contained_steps(steps: &[(u64, i32)], k: u64, lo: u64, hi: u64) -> std::ops::Range<usize> {
    let start = steps.partition_point(|&(bp, _)| bp < lo);
    let end = steps.partition_point(|&(bp, _)| bp.saturating_add(k) <= hi);
    start..end.max(start)
}

/// Sorted-key multiset intersection weight: Sum mult_a * mult_b over equal
/// keys (the <e_i, e_j> of the graph factorization, before depth scaling).
fn multiset_overlap(left: &[(u64, u32)], right: &[(u64, u32)]) -> f64 {
    let (mut i, mut j, mut total) = (0usize, 0usize, 0u64);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                total += left[i].1 as u64 * right[j].1 as u64;
                i += 1;
                j += 1;
            }
        }
    }
    total as f64
}

/// Per-locus graph row with the factored statistics. `node_dot`/`edge_dot`
/// are <x, e_i> over the observed mass; the norms are ||e_i||^2 per layer.
/// Node keys are widened to u64 (and edges packed) so the O(rows^2) pair
/// walk stays scalar.
struct GraphRow {
    nodes: Vec<(u64, u32)>,
    edges: Vec<(PackedEdge, u32)>,
    node_dot: f64,
    edge_dot: f64,
    node_norm: f64,
    edge_norm: f64,
    members: Vec<usize>,
}

struct GraphSpace {
    rows: Vec<GraphRow>,
    observed_node_norm: f64,
    observed_edge_norm: f64,
}

impl GraphSpace {
    /// One arm's pair score; `None` fails closed on a zero observed or
    /// expected norm (no silent award). Used for the assessment-side truth
    /// scores; the exhaustive pair walk inlines the same arithmetic with
    /// the overlaps shared between the two arms.
    fn cosine(&self, first: usize, second: usize, depth: f64, edges: bool) -> Option<f64> {
        let (a, b) = (&self.rows[first], &self.rows[second]);
        let node_overlap = multiset_overlap(&a.nodes, &b.nodes);
        let edge_overlap = if edges {
            multiset_overlap(&a.edges, &b.edges)
        } else {
            0.0
        };
        let (dot, expected_norm, observed_norm) = if edges {
            (
                a.node_dot + b.node_dot + a.edge_dot + b.edge_dot,
                a.node_norm + b.node_norm + a.edge_norm + b.edge_norm
                    + 2.0 * depth * depth * (node_overlap + edge_overlap),
                self.observed_node_norm + self.observed_edge_norm,
            )
        } else {
            (
                a.node_dot + b.node_dot,
                a.node_norm + b.node_norm + 2.0 * depth * depth * node_overlap,
                self.observed_node_norm,
            )
        };
        (observed_norm > 0.0 && expected_norm > 0.0)
            .then(|| dot / (observed_norm * expected_norm).sqrt())
    }
}

/// Multiset sum of two sorted usage vectors: a diplotype CLASS's total
/// usage (both copies merged; a double-copy class merges a row with
/// itself, doubling its multiplicities).
fn merged_multiset(left: &[(u64, u32)], right: &[(u64, u32)]) -> Vec<(u64, u32)> {
    let (mut i, mut j) = (0usize, 0usize);
    let mut merged: Vec<(u64, u32)> = Vec::with_capacity(left.len() + right.len());
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => {
                merged.push(left[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                merged.push(right[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                merged.push((left[i].0, left[i].1 + right[j].1));
                i += 1;
                j += 1;
            }
        }
    }
    merged.extend_from_slice(&left[i..]);
    merged.extend_from_slice(&right[j..]);
    merged
}

/// Symmetric multiset distance between two sorted usage vectors:
/// (number of keys whose multiplicities differ, total absolute
/// multiplicity difference). (0, 0) means bit-identical usage.
fn multiset_distance(left: &[(u64, u32)], right: &[(u64, u32)]) -> (u64, u64) {
    let (mut i, mut j, mut keys, mut mass) = (0usize, 0usize, 0u64, 0u64);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => {
                keys += 1;
                mass += left[i].1 as u64;
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                keys += 1;
                mass += right[j].1 as u64;
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                if left[i].1 != right[j].1 {
                    keys += 1;
                    mass += (left[i].1 as i64 - right[j].1 as i64).unsigned_abs();
                }
                i += 1;
                j += 1;
            }
        }
    }
    keys += (left.len() - i) as u64;
    mass += left[i..].iter().map(|&(_, mult)| mult as u64).sum::<u64>();
    keys += (right.len() - j) as u64;
    mass += right[j..].iter().map(|&(_, mult)| mult as u64).sum::<u64>();
    (keys, mass)
}

/// Deterministic FNV-1a hash over a row's exact sorted node+edge usage —
/// a receipt aid for comparing class signatures across entries. The
/// coalescing itself is exact structural equality, never this hash.
fn usage_hash(row: &GraphRow) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for &(key, mult) in row.nodes.iter().chain(row.edges.iter()) {
        for part in [key, mult as u64] {
            hash ^= part;
            hash = hash.wrapping_mul(0x100000001b3);
        }
    }
    hash
}

/// How many physical candidate-row pairs coalesce into ONE material class
/// (the owner's "effectively the same" made countable): a two-row class
/// merges members_a x members_b physical pairs; a double-copy class merges
/// a row's members with themselves, m*(m+1)/2 pairs, not m*m.
fn class_physical_pairs(same_row: bool, members_a: usize, members_b: usize) -> u64 {
    if same_row {
        members_a as u64 * (members_a as u64 + 1) / 2
    } else {
        members_a as u64 * members_b as u64
    }
}

/// Observed mass on the keys where two sorted usage multisets differ,
/// weighted by the multiplicity difference: the maximum |dot| contribution
/// the usage difference can carry. Exactly 0.0 means every differing
/// segment/adjacency holds NO observed read mass — no read breaks a tie
/// between the two usages.
fn differing_observed_mass(
    left: &[(u64, u32)],
    right: &[(u64, u32)],
    observed: &dyn Fn(u64) -> f64,
) -> f64 {
    let (mut i, mut j, mut total) = (0usize, 0usize, 0.0f64);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => {
                total += observed(left[i].0) * left[i].1 as f64;
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                total += observed(right[j].0) * right[j].1 as f64;
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                let delta = (left[i].1 as i64 - right[j].1 as i64).unsigned_abs();
                if delta > 0 {
                    total += observed(left[i].0) * delta as f64;
                }
                i += 1;
                j += 1;
            }
        }
    }
    while i < left.len() {
        total += observed(left[i].0) * left[i].1 as f64;
        i += 1;
    }
    while j < right.len() {
        total += observed(right[j].0) * right[j].1 as f64;
        j += 1;
    }
    total
}

/// The confidence in the emitted single-material draw under the forced
/// two-way normalization of the k bit-tied maximum hypotheses against the
/// best alternative a: p = s_win / (k*s_win + a). The delta and cluster
/// forms share this algebra — the delta form passes the bit-tied CLASS
/// count as k, the cluster form (the owner's ruling below) passes the
/// DIVERGENT-CLUSTER count. Within a bit-tied set the truth is uniform
/// (measured tie anatomy: the differing nodes/edges carry zero observed
/// read mass), so p = (k*s_win/(k*s_win + a)) * (1/k). None fails closed
/// when there is no similarity mass at all (s_win = 0 and a = 0) or no
/// called class.
fn qual_p(s_win: f64, called_classes: usize, alternative: Option<f64>) -> Option<f64> {
    if called_classes == 0 {
        return None;
    }
    let alternative = alternative.unwrap_or(0.0);
    let denominator = called_classes as f64 * s_win + alternative;
    (denominator > 0.0).then(|| s_win / denominator)
}

/// Phred-scaled QUAL from the delta-form confidence. None = unbounded
/// (p = 1: k = 1 with no alternative similarity — the called class is the
/// domain's only similarity); emitted as null, never clamped — no free
/// parameters, no cap.
fn qual_from_p(p: f64) -> Option<f64> {
    (p < 1.0).then(|| -10.0 * (1.0 - p).log10())
}

/// The median of a non-empty f64 slice (the average of the two middle
/// values for an even count). Inputs here are relative jumps in [0, 1]
/// by construction — finite, non-negative, NaN-free — so a total sort is
/// exact. The median is the derived "typical" jump the knee must dominate:
/// no threshold, no constant.
fn median_of(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).expect("finite relative jumps"));
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        values[middle]
    } else {
        (values[middle - 1] + values[middle]) / 2.0
    }
}

/// The KNEE of a sorted non-decreasing distance spectrum — the derived
/// cluster-cut rule, with NO tuning constants. A winner's near-identical
/// band is a dense run of small distances followed by the divergent tail,
/// so the band-to-tail transition is the MAXIMAL RELATIVE JUMP
/// r_i = (d_{i+1} - d_i)/d_{i+1} ∈ [0, 1] between consecutive distances
/// (0/0 := 0); ties resolve to the LARGEST index so the band absorbs the
/// full dense run, and the cut distance is d_{i*} — the last band point —
/// so the called cluster is every class within the cut. A spectrum whose
/// maximal relative jump does not exceed the MEDIAN relative jump is
/// scale-free (exact geometric growth: every jump equal) and HAS NO KNEE
/// — a measurement, not an error: the cut then stays at 0 (the
/// bit-identical innermost level) and the caller reports the nearest
/// strictly-positive distance as the fallback. The exactly-zero distances
/// (bit-identical material, or differing material carrying no observed
/// mass) are the innermost band and are always inside the cluster, so the
/// knee is derived over the STRICTLY-POSITIVE distances only — a jump out
/// of an exact zero is the maximal relative jump r = 1 by definition and
/// would otherwise pin the cut at the zero boundary and split the band.
/// With fewer than two positive distances there is no jump that can
/// dominate the median, so no knee is claimed.
struct KneeCut {
    has_knee: bool,
    cut: f64,
}

fn spectrum_knee(sorted: &[f64]) -> KneeCut {
    let positives: Vec<f64> = sorted.iter().copied().filter(|&d| d > 0.0).collect();
    if positives.len() < 2 {
        return KneeCut { has_knee: false, cut: 0.0 };
    }
    let jumps: Vec<f64> = positives
        .windows(2)
        .map(|pair| (pair[1] - pair[0]) / pair[1])
        .collect();
    let mut center = jumps.clone();
    let typical = median_of(&mut center);
    let (mut best, mut best_index) = (f64::NEG_INFINITY, 0usize);
    for (index, &jump) in jumps.iter().enumerate() {
        if jump >= best {
            best = jump;
            best_index = index;
        }
    }
    if best > typical {
        KneeCut { has_knee: true, cut: positives[best_index] }
    } else {
        KneeCut { has_knee: false, cut: 0.0 }
    }
}

/// The SIGNATURE-COSINE distance between two class usage signatures:
/// 1 - cosine of the RAW node+edge usage multisets (multiplicity vectors;
/// the per-copy depth is a common positive factor and cancels, so raw
/// multiplicities give the same cosine). The second view of material
/// distance — scale-free and evidence-independent — reported beside the
/// observed-mass distance, which is the one the cluster cut uses.
fn signature_cosine_distance(
    nodes_a: &[(u64, u32)],
    edges_a: &[(u64, u32)],
    nodes_b: &[(u64, u32)],
    edges_b: &[(u64, u32)],
) -> f64 {
    let norm = |nodes: &[(u64, u32)], edges: &[(u64, u32)]| {
        nodes.iter().map(|&(_, mult)| mult as f64 * mult as f64).sum::<f64>()
            + edges.iter().map(|&(_, mult)| mult as f64 * mult as f64).sum::<f64>()
    };
    let (norm_a, norm_b) = (norm(nodes_a, edges_a), norm(nodes_b, edges_b));
    if norm_a == 0.0 && norm_b == 0.0 {
        return 0.0; // two empty signatures are the same (empty) usage
    }
    if norm_a == 0.0 || norm_b == 0.0 {
        return 1.0; // an empty signature shares nothing with a nonempty one
    }
    let overlap = multiset_overlap(nodes_a, nodes_b) + multiset_overlap(edges_a, edges_b);
    1.0 - overlap / (norm_a * norm_b).sqrt()
}

/// The CLUSTER-FORM confidence in the emitted single-material draw (the
/// owner's ruling, 2026-10-01: the near-twin runner-ups are not real
/// alternatives — "there's 10 sequences that are almost identical and we
/// pick one of them; they shouldn't get penalized as if they're not" — so
/// the quality must measure the separation from the NEXT ACTUALLY
/// DIVERGENT cluster of sequences, not from the ninth twin; picking any
/// member of a near-identical band is a fine sequence call). The delta
/// form's formula is UNCHANGED but applied BETWEEN CLUSTERS:
/// p = s_win/(k*s_win + a), where s_win is the bit-identical maximum
/// similarity, k is the number of DIVERGENT clusters whose best candidates
/// bit-tie the maximum (the bit-tied called classes clustered by the same
/// near-identity relation — single-linkage components at the knee cut,
/// chains included — the winner's own cluster counts), and a is the best
/// similarity among classes OUTSIDE the near-identical called cluster(s) —
/// the best genuinely divergent rival. Members of the called cluster(s)
/// never supply a, never count in k, never lower QUAL; bit-identical
/// classes (distance 0) are innermost and always inside the cluster.
/// Derivation: the hypotheses are the k divergent tied-maximum clusters
/// (each carrying s_win) and the best divergent rival (carrying a); the
/// forced two-way normalization k*s_win vs a with uniform truth inside a
/// bit-tied cluster gives the confidence in the single emitted material
/// draw p = (k*s_win/(k*s_win + a))*(1/k) = s_win/(k*s_win + a) — the
/// delta form's algebra with clusters as the units.
struct ClusterQual {
    /// The knee distance of the winner's spectrum (None: no knee — a
    /// measurement; the cut then stays at the bit-identical level 0).
    knee: Option<f64>,
    /// The cluster cut actually applied (the knee distance, or 0 without a
    /// knee: the bit-identical innermost coalescing).
    cut: f64,
    shape: &'static str,
    /// Classes in the winner's near-identical cluster (the winner itself
    /// included).
    cluster_size: usize,
    /// The number of DIVERGENT clusters whose best candidates bit-tie the
    /// maximum.
    k: usize,
    /// Exclusion mask over `spectrum` (true = inside a called cluster:
    /// supplies no alternative, counts in no k).
    excluded: Vec<bool>,
    /// The best similarity among classes OUTSIDE the called cluster(s)
    /// (None iff every eligible class is inside; Some(0.0) when the best
    /// divergent rival carries no similarity).
    alternative: Option<f64>,
    p: Option<f64>,
}

fn cluster_form_qual(s_win: f64, spectrum: &[(f64, f64)], tied: &[(usize, &[f64])]) -> ClusterQual {
    let mut sorted: Vec<f64> = spectrum.iter().map(|&(distance, _)| distance).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
    let knee_cut = spectrum_knee(&sorted);
    let cut = knee_cut.cut;
    let mut excluded = vec![false; spectrum.len()];
    let mut band = 0usize;
    for (index, &(distance, _)) in spectrum.iter().enumerate() {
        if distance <= cut {
            excluded[index] = true;
            band += 1;
        }
    }
    for (_, band_distances) in tied {
        for (index, &distance) in band_distances.iter().enumerate() {
            if distance <= cut {
                excluded[index] = true;
            }
        }
    }
    // k: single-linkage components over the called classes (the winner is
    // node 0, tied class i is node i + 1). Two called classes are
    // near-identical iff their distance <= cut; chains count (a divergent
    // twin of the winner's twin is NOT a new cluster).
    let mut parent: Vec<usize> = (0..=tied.len()).collect();
    let mut root = |parent: &mut Vec<usize>, node: usize| {
        let mut current = node;
        while parent[current] != current {
            current = parent[current];
        }
        let mut walked = node;
        while parent[walked] != current {
            let next = parent[walked];
            parent[walked] = current;
            walked = next;
        }
        current
    };
    for (left, (index_left, _)) in tied.iter().enumerate() {
        // winner (0) vs tied class left + 1: the tied class's own
        // distance-to-winner entry in the spectrum.
        if spectrum[*index_left].0 <= cut {
            let (a, b) = (root(&mut parent, 0), root(&mut parent, left + 1));
            if a != b {
                parent[a] = b;
            }
        }
        for (right, (index_right, _)) in tied.iter().enumerate().skip(left + 1) {
            if tied[left].1[*index_right] <= cut {
                let (a, b) =
                    (root(&mut parent, left + 1), root(&mut parent, right + 1));
                if a != b {
                    parent[a] = b;
                }
            }
        }
    }
    let mut roots: Vec<usize> = (0..=tied.len()).map(|node| root(&mut parent, node)).collect();
    roots.sort_unstable();
    roots.dedup();
    let k = roots.len();
    let mut alternative: Option<f64> = None;
    for (index, &(_, score)) in spectrum.iter().enumerate() {
        if !excluded[index] {
            alternative = Some(alternative.map_or(score, |value| value.max(score)));
        }
    }
    let shape = if spectrum.is_empty() {
        "single_class"
    } else if knee_cut.has_knee {
        "knee"
    } else {
        "no_knee"
    };
    let p = qual_p(s_win, k, alternative);
    ClusterQual {
        knee: knee_cut.has_knee.then_some(cut),
        cut,
        shape,
        cluster_size: 1 + band,
        k,
        excluded,
        alternative,
        p,
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn dump_graph_exhaustive(
    path: &str,
    panel: &SyngIndex,
    ranges: &[Vec<genome::SpanningTraversal>],
    path_of_source: &[usize],
    truth_pieces: &[[Vec<(usize, u64, u64, bool, u32)>; 2]],
    instances: &crate::InstanceStructure,
    k: u64,
    depth: f64,
    locus_offset: usize,
) -> io::Result<()> {
    let mut report = BufWriter::new(std::fs::File::create(path)?);
    let mut competitors =
        BufWriter::new(std::fs::File::create(format!("{path}.competitors.jsonl"))?);
    // Classes within the stage-1 1e-12 window of the truth class (including
    // the truth class itself): the tie-evidence stream for the exact-tie
    // vs near-identical-signature question. Diagnostic only.
    let mut ties = BufWriter::new(std::fs::File::create(format!("{path}.ties.jsonl"))?);
    for (locus, candidates) in ranges.iter().enumerate() {
        // Per path: the candidate rows' merged material intervals and the
        // contributing records' merged placement spans (all features, both
        // frames, unioned per path — the record votes its share at a shared
        // segment once, however many subwalks or conserved paths cover it).
        struct PathWork {
            rows: Vec<(usize, u64, u64)>,
            records: Vec<(usize, u64, u64)>,
        }
        impl Default for PathWork {
            fn default() -> Self {
                PathWork { rows: Vec::new(), records: Vec::new() }
            }
        }
        let mut paths: std::collections::BTreeMap<usize, PathWork> =
            std::collections::BTreeMap::new();
        let mut row_materials: Vec<Material> = Vec::with_capacity(candidates.len());
        for (index, row) in candidates.iter().enumerate() {
            let material = row_material(row, path_of_source);
            for &(path, lo, hi) in &material {
                paths
                    .entry(path)
                    .or_default()
                    .rows
                    .push((index, lo, hi));
            }
            row_materials.push(material);
        }
        let mut records: Vec<u32> = instances.window_records[locus]
            .values()
            .flat_map(|ids| ids.iter().copied())
            .collect();
        records.sort_unstable();
        records.dedup();
        for &record in &records {
            let mut per_path: std::collections::BTreeMap<usize, Vec<(u64, u64)>> =
                std::collections::BTreeMap::new();
            for spans in instances.record_spans[record as usize].values() {
                for &(path, lo, hi) in spans {
                    per_path.entry(path).or_default().push((lo, hi));
                }
            }
            for (path, mut intervals) in per_path {
                merge_intervals(&mut intervals);
                for (lo, hi) in intervals {
                    paths
                        .entry(path)
                        .or_default()
                        .records
                        .push((record as usize, lo, hi));
                }
            }
        }
        // One bounding-range path walk per path; per-interval window
        // containment extracts nodes and the adjacencies between
        // consecutive contained steps.
        let mut row_nodes: Vec<std::collections::BTreeMap<GraphNode, u32>> =
            vec![std::collections::BTreeMap::new(); candidates.len()];
        let mut row_edges: Vec<std::collections::BTreeMap<PackedEdge, u32>> =
            vec![std::collections::BTreeMap::new(); candidates.len()];
        let mut record_nodes: HashMap<u32, std::collections::BTreeSet<GraphNode>> =
            HashMap::new();
        let mut record_edges: HashMap<u32, std::collections::BTreeSet<PackedEdge>> =
            HashMap::new();
        for (path, work) in &paths {
            let lo_min = work
                .rows
                .iter()
                .chain(work.records.iter())
                .map(|&(_, lo, _)| lo)
                .min()
                .unwrap_or(0);
            let hi_max = work
                .rows
                .iter()
                .chain(work.records.iter())
                .map(|&(_, _, hi)| hi)
                .max()
                .unwrap_or(0);
            if hi_max <= lo_min {
                continue;
            }
            let mut steps: Vec<(u64, i32)> = panel
                .walk_path_range(*path, lo_min, hi_max)?
                .into_iter()
                .map(|(node, bp)| (bp, node))
                .collect();
            steps.sort_unstable_by_key(|&(bp, _)| bp);
            let mut add_interval =
                |lo: u64, hi: u64, nodes: &mut dyn FnMut(GraphNode), edges: &mut dyn FnMut(GraphEdge)| {
                    let window = contained_steps(&steps, k, lo, hi);
                    let mut previous: Option<GraphNode> = None;
                    for &(_, node) in &steps[window] {
                        let key = node.unsigned_abs();
                        if let Some(left) = previous.take() {
                            edges((left, key));
                        }
                        previous = Some(key);
                        nodes(key);
                    }
                };
            for &(row_index, lo, hi) in &work.rows {
                let nodes = &mut row_nodes[row_index];
                let edges = &mut row_edges[row_index];
                add_interval(
                    lo,
                    hi,
                    &mut |key| *nodes.entry(key).or_default() += 1,
                    &mut |(left, right)| *edges.entry(pack_edge(left, right)).or_default() += 1,
                );
            }
            for &(record, lo, hi) in &work.records {
                let record = record as u32;
                let nodes = record_nodes.entry(record).or_default();
                let edges = record_edges.entry(record).or_default();
                add_interval(
                    lo,
                    hi,
                    &mut |key| {
                        nodes.insert(key);
                    },
                    &mut |(left, right)| {
                        edges.insert(pack_edge(left, right));
                    },
                );
            }
        }
        // The fixed observation universe: the nodes/edges covered by ANY
        // candidate row's spelled material (mass outside it is dropped,
        // exactly the row-space candidate-material convention).
        let mut universe_nodes: std::collections::BTreeSet<GraphNode> =
            std::collections::BTreeSet::new();
        let mut universe_edges: std::collections::BTreeSet<PackedEdge> =
            std::collections::BTreeSet::new();
        for nodes in &row_nodes {
            universe_nodes.extend(nodes.keys().copied());
        }
        for edges in &row_edges {
            universe_edges.extend(edges.keys().copied());
        }
        // Observed mass: per record, its routed share once per covered
        // universe segment/adjacency (the record-once dedup across paths
        // IS the near-twin collapse: one read, one shared segment, one
        // vote, whatever the number of conserved placements).
        let mut observed_nodes: HashMap<GraphNode, f64> = HashMap::new();
        let mut observed_edges: HashMap<PackedEdge, f64> = HashMap::new();
        let mut observed_node_mass = 0.0;
        let mut observed_edge_mass = 0.0;
        for &record in &records {
            let share = instances.record_shares[record as usize];
            if let Some(nodes) = record_nodes.get(&record) {
                for &node in nodes {
                    if universe_nodes.contains(&node) {
                        *observed_nodes.entry(node).or_default() += share;
                        observed_node_mass += share;
                    }
                }
            }
            if let Some(edges) = record_edges.get(&record) {
                for &edge in edges {
                    if universe_edges.contains(&edge) {
                        *observed_edges.entry(edge).or_default() += share;
                        observed_edge_mass += share;
                    }
                }
            }
        }
        let observed_node_norm: f64 =
            observed_nodes.values().map(|mass| mass * mass).sum();
        let observed_edge_norm: f64 =
            observed_edges.values().map(|mass| mass * mass).sum();
        // Exact-usage coalescing (a twin that spells the same segments and
        // adjacencies is ONE graph row; identities stay inspectable via
        // `members`).
        let mut groups: std::collections::BTreeMap<GraphUsage, Vec<usize>> =
            std::collections::BTreeMap::new();
        for (index, (nodes, edges)) in row_nodes.iter().zip(&row_edges).enumerate() {
            groups
                .entry(GraphUsage {
                    nodes: nodes.iter().map(|(&key, &mult)| (key, mult)).collect(),
                    edges: edges.iter().map(|(&key, &mult)| (key, mult)).collect(),
                })
                .or_default()
                .push(index);
        }
        let graph_rows = groups.len();
        let rows: Vec<GraphRow> = groups
            .into_iter()
            .map(|(usage, members)| {
                let node_dot = usage
                    .nodes
                    .iter()
                    .map(|&(key, mult)| {
                        observed_nodes.get(&key).copied().unwrap_or(0.0) * depth * mult as f64
                    })
                    .sum();
                let edge_dot = usage
                    .edges
                    .iter()
                    .map(|&(key, mult)| {
                        observed_edges.get(&key).copied().unwrap_or(0.0) * depth * mult as f64
                    })
                    .sum();
                let node_norm: f64 = depth
                    * depth
                    * usage
                        .nodes
                        .iter()
                        .map(|&(_, mult)| mult as f64 * mult as f64)
                        .sum::<f64>();
                let edge_norm: f64 = depth
                    * depth
                    * usage
                        .edges
                        .iter()
                        .map(|&(_, mult)| mult as f64 * mult as f64)
                        .sum::<f64>();
                GraphRow {
                    nodes: usage
                        .nodes
                        .into_iter()
                        .map(|(key, mult)| (key as u64, mult))
                        .collect(),
                    edges: usage.edges,
                    node_dot,
                    edge_dot,
                    node_norm,
                    edge_norm,
                    members,
                }
            })
            .collect();
        let space = GraphSpace {
            rows,
            observed_node_norm,
            observed_edge_norm,
        };
        let mut id_to_row = vec![0usize; candidates.len()];
        for (index, row) in space.rows.iter().enumerate() {
            for &id in &row.members {
                id_to_row[id] = index;
            }
        }
        let truth = [0, 1].map(|copy| truth_index(candidates, &truth_pieces[locus][copy]));
        let pair_truth = match (truth[0], truth[1]) {
            (Some(a), Some(b)) => Some((id_to_row[a], id_to_row[b])),
            _ => None,
        };
        let truth_nodes_score = pair_truth.and_then(|(a, b)| space.cosine(a, b, depth, false));
        let truth_combined_score = pair_truth.and_then(|(a, b)| space.cosine(a, b, depth, true));
        // Product QUAL state (material-class semantics, CLUSTER FORM). The
        // exhaustive walk below collects the classes at the bit-identical
        // maximum and stores every eligible class's similarity for the
        // post-walk cluster analysis (the material-distance spectrum from
        // the winner, the derived knee cut, the divergent-cluster count k
        // and the best genuinely divergent rival a). Exact-tie convention:
        // identical f64 score values — every score is finite, non-negative,
        // NaN-free, and -0.0 cannot arise (non-negative dots over positive
        // norms), so IEEE-754 equality IS bit identity here; NO epsilon
        // constant enters the product call.
        let truth_merged_nodes = pair_truth.map(|(a, b)| {
            merged_multiset(&space.rows[a].nodes, &space.rows[b].nodes)
        });
        let truth_merged_edges = pair_truth.map(|(a, b)| {
            merged_multiset(&space.rows[a].edges, &space.rows[b].edges)
        });
        let mut combined_total = 0.0f64;
        let mut called_score = f64::NEG_INFINITY;
        let mut called: Vec<[usize; 2]> = Vec::new();
        // Every eligible material class (its row pair and similarity), in
        // the deterministic enumeration order — the input of the winner's
        // distance spectrum. Diagnostic continuity: the total is still
        // emitted (`qual_similarity_total`) but no QUAL form uses it — it
        // was the falsified share form's denominator.
        let mut eligible_classes: Vec<(usize, usize, f64)> = Vec::new();
        let mut ties_streamed = 0u64;
        // Exhaustive assessment: every unordered graph-row pair, both arms.
        // The overlaps are computed once per pair and shared by both arms'
        // scores (the row-space factorization, one arithmetic pass).
        let mut counts = [(0u64, 0u64); 2];
        let mut eligible = [0u64; 2];
        let mut best: [Option<(f64, [usize; 2])>; 2] = [None, None];
        let mut competitor_entries = 0u64;
        for second in 0..graph_rows {
            for first in 0..=second {
                let (a, b) = (&space.rows[first], &space.rows[second]);
                let node_overlap = multiset_overlap(&a.nodes, &b.nodes);
                let edge_overlap = multiset_overlap(&a.edges, &b.edges);
                let nodes_score = (space.observed_node_norm > 0.0
                    && a.node_norm + b.node_norm + 2.0 * depth * depth * node_overlap > 0.0)
                    .then(|| {
                        (a.node_dot + b.node_dot)
                            / (space.observed_node_norm
                                * (a.node_norm + b.node_norm
                                    + 2.0 * depth * depth * node_overlap))
                                .sqrt()
                    });
                let expected_combined = a.node_norm
                    + b.node_norm
                    + a.edge_norm
                    + b.edge_norm
                    + 2.0 * depth * depth * (node_overlap + edge_overlap);
                let combined_score =
                    (space.observed_node_norm + space.observed_edge_norm > 0.0
                        && expected_combined > 0.0)
                        .then(|| {
                            (a.node_dot + b.node_dot + a.edge_dot + b.edge_dot)
                                / ((space.observed_node_norm + space.observed_edge_norm)
                                    * expected_combined)
                                    .sqrt()
                        });
                // Product QUAL accumulation: one value per distinct material
                // class; the called set is the classes at the bit-identical
                // maximum, and every eligible class's similarity is stored
                // for the post-walk cluster analysis.
                if let Some(score) = combined_score {
                    combined_total += score;
                    eligible_classes.push((first, second, score));
                    if score > called_score {
                        called_score = score;
                        called.clear();
                        called.push([first, second]);
                    } else if score == called_score {
                        called.push([first, second]);
                    }
                }
                // Tie-evidence stream: every class inside the stage-1 1e-12
                // window of the truth class (a diagnostic convention; the
                // product called set above is bit-exact only), carrying the
                // measured ULP delta and the node/edge usage distance versus
                // the truth class's total usage.
                if let (Some(truth_value), Some(score)) = (truth_combined_score, combined_score) {
                    if (score - truth_value).abs() <= 1e-12 {
                        let (ta, tb) = pair_truth.unwrap();
                        let merged_nodes = merged_multiset(&a.nodes, &b.nodes);
                        let merged_edges = merged_multiset(&a.edges, &b.edges);
                        let node_distance =
                            multiset_distance(&merged_nodes, truth_merged_nodes.as_ref().unwrap());
                        let edge_distance =
                            multiset_distance(&merged_edges, truth_merged_edges.as_ref().unwrap());
                        // The observed mass sitting on the differing segments:
                        // 0.0 means no read's evidence distinguishes the two
                        // classes (the tie is unbreakable by this sample).
                        let node_differing_mass = differing_observed_mass(
                            &merged_nodes,
                            truth_merged_nodes.as_ref().unwrap(),
                            &|key| {
                                observed_nodes.get(&(key as u32)).copied().unwrap_or(0.0)
                            },
                        );
                        let edge_differing_mass = differing_observed_mass(
                            &merged_edges,
                            truth_merged_edges.as_ref().unwrap(),
                            &|key| observed_edges.get(&key).copied().unwrap_or(0.0),
                        );
                        serde_json::to_writer(&mut ties, &serde_json::json!({
                            "locus": locus + locus_offset,
                            "row_indices": [first, second],
                            "identities": [candidates[a.members[0]].identity,
                                           candidates[b.members[0]].identity],
                            "nodes_cosine": nodes_score,
                            "combined_cosine": combined_score,
                            "truth_combined_cosine": truth_combined_score,
                            "combined_ulp_delta":
                                (score.to_bits() as i64 - truth_value.to_bits() as i64).abs(),
                            "combined_bit_exact": score.to_bits() == truth_value.to_bits(),
                            "shared_rows_with_truth": usize::from(first == ta)
                                + usize::from(first == tb)
                                + usize::from(second == ta)
                                + usize::from(second == tb),
                            "node_distance_vs_truth_class": [node_distance.0, node_distance.1],
                            "edge_distance_vs_truth_class": [edge_distance.0, edge_distance.1],
                            "node_differing_observed_mass_vs_truth": node_differing_mass,
                            "edge_differing_observed_mass_vs_truth": edge_differing_mass,
                            "row_member_counts": [a.members.len(), b.members.len()],
                            "is_truth_class": (first == ta && second == tb)
                                || (first == tb && second == ta),
                        }))?;
                        writeln!(ties)?;
                        ties_streamed += 1;
                    }
                }
                for (arm_index, score) in [nodes_score, combined_score].into_iter().enumerate() {
                    if let Some(score) = score {
                        eligible[arm_index] += 1;
                        if best[arm_index].is_none_or(|(old, _)| score > old) {
                            best[arm_index] = Some((score, [first, second]));
                        }
                    }
                }
                if let (Some(truth_value), Some(score)) = (truth_nodes_score, nodes_score) {
                    if score > truth_value + 1e-12 {
                        counts[0].0 += 1;
                    } else if (score - truth_value).abs() <= 1e-12 {
                        counts[0].1 += 1;
                    }
                }
                if let (Some(truth_value), Some(score)) = (truth_combined_score, combined_score) {
                    if score > truth_value + 1e-12 {
                        counts[1].0 += 1;
                    } else if (score - truth_value).abs() <= 1e-12 {
                        counts[1].1 += 1;
                    }
                }
                let better_nodes = truth_nodes_score
                    .zip(nodes_score)
                    .is_some_and(|(t, s)| s > t + 1e-12);
                let better_combined = truth_combined_score
                    .zip(combined_score)
                    .is_some_and(|(t, s)| s > t + 1e-12);
                if better_nodes || better_combined {
                    let ids = [space.rows[first].members[0], space.rows[second].members[0]];
                    let mut arms = Vec::new();
                    if better_nodes {
                        arms.push("nodes");
                    }
                    if better_combined {
                        arms.push("combined");
                    }
                    serde_json::to_writer(&mut competitors, &serde_json::json!({
                        "locus": locus + locus_offset,
                        "arms": arms,
                        "row_indices": ids,
                        "identities": [candidates[ids[0]].identity, candidates[ids[1]].identity],
                        "node_counts": [space.rows[first].nodes.len(), space.rows[second].nodes.len()],
                        "edge_counts": [space.rows[first].edges.len(), space.rows[second].edges.len()],
                        "nodes_cosine": nodes_score,
                        "combined_cosine": combined_score,
                        "truth_nodes_cosine": truth_nodes_score,
                        "truth_combined_cosine": truth_combined_score,
                    }))?;
                    writeln!(competitors)?;
                    competitor_entries += 1;
                }
            }
        }
        let mut identity_kinds: HashMap<&str, u64> = HashMap::new();
        for row in candidates.iter() {
            *identity_kinds
                .entry(row.identity.split(':').next().unwrap_or("?"))
                .or_default() += 1;
        }
        let multi_segment_rows = candidates
            .iter()
            .filter(|row| row.segments.len() > 1)
            .count();
        // Rows whose spelled material is NOT path-continuous: these would
        // be the only place a junction adjacency that no single record
        // spans could hide (measured audit of the chain-channel claim).
        let disjoint_material_rows = candidates
            .iter()
            .zip(&row_materials)
            .filter(|(row, material)| {
                row.segments.iter().map(|s| s.end.saturating_sub(s.start)).sum::<u64>()
                    > row_material_merged_bp(material)
            })
            .count();
        let ranks = [
            truth_nodes_score.map(|_| counts[0].0 + 1),
            truth_combined_score.map(|_| counts[1].0 + 1),
        ];
        let truth_in_called_set = pair_truth.map(|(ta, tb)| {
            called.iter().any(|&[first, second]| {
                (first == ta && second == tb) || (first == tb && second == ta)
            })
        });
        let first_called_nodes = called.first().map(|&[first, second]| {
            merged_multiset(&space.rows[first].nodes, &space.rows[second].nodes)
        });
        let first_called_edges = called.first().map(|&[first, second]| {
            merged_multiset(&space.rows[first].edges, &space.rows[second].edges)
        });
        // Product QUAL (material-class semantics, CLUSTER FORM — the
        // owner's ruling: the near-twin runner-ups are not real
        // alternatives; the quality measures the separation from the NEXT
        // ACTUALLY DIVERGENT cluster, not from the ninth twin). The called
        // set is unchanged (distinct classes at the bit-identical maximum).
        // The MATERIAL DISTANCE between two classes is the observed read
        // mass carried on the symmetric difference of their node+edge
        // usage signatures (the same routed-share units the cosines use),
        // with the signature-cosine distance as a second view. The cluster
        // cut is the derived knee of the winner's sorted distance
        // spectrum; k = the number of DIVERGENT clusters whose best
        // candidates bit-tie the maximum; a = the best similarity OUTSIDE
        // the near-identical called cluster(s); p = s_win/(k*s_win + a),
        // QUAL = -10*log10(1 - p). Members of the called cluster never
        // supply a, never count in k, never lower QUAL. Bit-identical
        // classes (distance 0) are innermost: always in the cluster.
        // Unbounded (p = 1) and massless loci emit null, never a clamp.
        // Member multiplicity enters NOWHERE.
        let winner = called.first().copied();
        // The winner's distance spectrum: (material distance, similarity,
        // signature-cosine distance, row pair) per eligible class other
        // than the winner, sorted by distance — the per-locus evidence of
        // where the dense near-identical band ends and where the cut
        // lands.
        let mut spectrum: Vec<(f64, f64, f64, usize, usize)> = Vec::new();
        if let (
            Some([winner_first, winner_second]),
            Some(winner_nodes),
            Some(winner_edges),
        ) = (
            winner,
            first_called_nodes.as_ref(),
            first_called_edges.as_ref(),
        ) {
            for &(first, second, score) in &eligible_classes {
                if first == winner_first && second == winner_second {
                    continue;
                }
                let nodes = merged_multiset(&space.rows[first].nodes, &space.rows[second].nodes);
                let edges = merged_multiset(&space.rows[first].edges, &space.rows[second].edges);
                let distance = differing_observed_mass(
                    &nodes,
                    winner_nodes,
                    &|key| observed_nodes.get(&(key as u32)).copied().unwrap_or(0.0),
                ) + differing_observed_mass(
                    &edges,
                    winner_edges,
                    &|key| observed_edges.get(&key).copied().unwrap_or(0.0),
                );
                let signature_distance =
                    signature_cosine_distance(&nodes, &edges, winner_nodes, winner_edges);
                spectrum.push((distance, score, signature_distance, first, second));
            }
        }
        spectrum.sort_by(|a, b| a.partial_cmp(b).expect("finite spectrum entries"));
        // Per other bit-tied called class: its spectrum index (for the k
        // components) and its distances to every spectrum member (for the
        // band exclusion).
        let mut tied_bands: Vec<(usize, Vec<f64>)> = Vec::new();
        for &[first, second] in called.iter().skip(1) {
            let index = spectrum
                .iter()
                .position(|&(_, _, _, a, b)| a == first && b == second)
                .expect("called class missing from the winner's spectrum");
            let nodes = merged_multiset(&space.rows[first].nodes, &space.rows[second].nodes);
            let edges = merged_multiset(&space.rows[first].edges, &space.rows[second].edges);
            let band: Vec<f64> = spectrum
                .iter()
                .map(|&(_, _, _, other_first, other_second)| {
                    let other_nodes = merged_multiset(
                        &space.rows[other_first].nodes,
                        &space.rows[other_second].nodes,
                    );
                    let other_edges = merged_multiset(
                        &space.rows[other_first].edges,
                        &space.rows[other_second].edges,
                    );
                    differing_observed_mass(
                        &nodes,
                        &other_nodes,
                        &|key| observed_nodes.get(&(key as u32)).copied().unwrap_or(0.0),
                    ) + differing_observed_mass(
                        &edges,
                        &other_edges,
                        &|key| observed_edges.get(&key).copied().unwrap_or(0.0),
                    )
                })
                .collect();
            tied_bands.push((index, band));
        }
        let tied_refs: Vec<(usize, &[f64])> = tied_bands
            .iter()
            .map(|(index, band)| (*index, band.as_slice()))
            .collect();
        let cluster = (!called.is_empty()).then(|| {
            let pairs: Vec<(f64, f64)> = spectrum
                .iter()
                .map(|&(distance, score, _, _, _)| (distance, score))
                .collect();
            cluster_form_qual(called_score, &pairs, &tied_refs)
        });
        let cluster_knee = cluster.as_ref().and_then(|state| state.knee);
        let cluster_shape = cluster.as_ref().map(|state| state.shape);
        let cluster_size = cluster.as_ref().map(|state| state.cluster_size);
        let cluster_k = cluster.as_ref().map(|state| state.k);
        let cluster_alternative = cluster.as_ref().and_then(|state| state.alternative);
        let cluster_excluded = cluster
            .as_ref()
            .map(|state| state.excluded.iter().filter(|excluded| **excluded).count());
        let confidence = cluster.as_ref().and_then(|state| state.p);
        let qual = confidence.and_then(qual_from_p);
        let qual_unbounded = confidence.is_some_and(|value| value >= 1.0);
        let qual_delta = cluster_alternative.map(|best_alternative| called_score - best_alternative);
        // The best genuinely divergent rival, NAMED: every class achieving a
        // outside the called clusters, with its distance and identities.
        let mut qual_alternative_classes: Vec<serde_json::Value> = Vec::new();
        if let (Some(state), Some(alternative)) = (&cluster, cluster_alternative) {
            for (index, ((distance, score, signature_distance, first, second), excluded)) in
                spectrum.iter().zip(&state.excluded).enumerate()
            {
                if !excluded && *score == alternative {
                    let (first, second) = (*first, *second);
                    qual_alternative_classes.push(serde_json::json!({
                        "spectrum_index": index,
                        "distance": distance,
                        "signature_distance": signature_distance,
                        "similarity": score,
                        "row_indices": [first, second],
                        "identities": [candidates[space.rows[first].members[0]].identity,
                                       candidates[space.rows[second].members[0]].identity],
                    }));
                }
            }
        }
        let spectrum_zero_distance = (!called.is_empty())
            .then(|| spectrum.iter().filter(|&&(distance, ..)| distance == 0.0).count());
        let nearest_rival_distance = spectrum
            .iter()
            .find(|&&(distance, ..)| distance > 0.0)
            .map(|&(distance, ..)| distance);
        // Second view (diagnostic only): the same derived knee on the
        // signature-cosine distance spectrum; the cluster cut uses the
        // observed-mass distance.
        let mut signature_sorted: Vec<f64> = spectrum
            .iter()
            .map(|&(_, _, signature, _, _)| signature)
            .collect();
        signature_sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite signature distances"));
        let signature_view = (!spectrum.is_empty()).then(|| spectrum_knee(&signature_sorted))
            .map(|knee| {
                let cluster_size = 1 + spectrum
                    .iter()
                    .filter(|&&(_, _, signature, _, _)| signature <= knee.cut)
                    .count();
                (knee.has_knee, knee.cut, cluster_size)
            });
        let qual_tied_class_bands: Vec<serde_json::Value> = called
            .iter()
            .skip(1)
            .zip(&tied_bands)
            .map(|(&[first, second], (index, band))| {
                serde_json::json!({
                    "row_indices": [first, second],
                    "spectrum_index": index,
                    "distance_to_winner": spectrum[*index].0,
                    "band_distances": band,
                })
            })
            .collect();
        let mut qual_called_physical_pairs = 0u64;
        let qual_called_classes: Vec<serde_json::Value> = called
            .iter()
            .enumerate()
            .map(|(position, &[first, second])| {
                let (row_a, row_b) = (&space.rows[first], &space.rows[second]);
                let physical =
                    class_physical_pairs(first == second, row_a.members.len(), row_b.members.len());
                qual_called_physical_pairs += physical;
                let node_distance = if position == 0 {
                    [0u64, 0u64]
                } else {
                    let merged = merged_multiset(&row_a.nodes, &row_b.nodes);
                    let distance = multiset_distance(&merged, first_called_nodes.as_ref().unwrap());
                    [distance.0, distance.1]
                };
                let edge_distance = if position == 0 {
                    [0u64, 0u64]
                } else {
                    let merged = merged_multiset(&row_a.edges, &row_b.edges);
                    let distance = multiset_distance(&merged, first_called_edges.as_ref().unwrap());
                    [distance.0, distance.1]
                };
                let node_differing_mass = if position == 0 {
                    0.0
                } else {
                    differing_observed_mass(
                        &merged_multiset(&row_a.nodes, &row_b.nodes),
                        first_called_nodes.as_ref().unwrap(),
                        &|key| observed_nodes.get(&(key as u32)).copied().unwrap_or(0.0),
                    )
                };
                let edge_differing_mass = if position == 0 {
                    0.0
                } else {
                    differing_observed_mass(
                        &merged_multiset(&row_a.edges, &row_b.edges),
                        first_called_edges.as_ref().unwrap(),
                        &|key| observed_edges.get(&key).copied().unwrap_or(0.0),
                    )
                };
                serde_json::json!({
                    "row_indices": [first, second],
                    "combined_cosine": (!called.is_empty()).then_some(called_score),
                    "row_member_counts": [row_a.members.len(), row_b.members.len()],
                    "physical_pair_members": physical,
                    "row_node_counts": [row_a.nodes.len(), row_b.nodes.len()],
                    "row_edge_counts": [row_a.edges.len(), row_b.edges.len()],
                    "row_usage_hashes": [usage_hash(row_a), usage_hash(row_b)],
                    "node_distance_to_first_called_class": node_distance,
                    "edge_distance_to_first_called_class": edge_distance,
                    "node_differing_observed_mass_to_first_called": node_differing_mass,
                    "edge_differing_observed_mass_to_first_called": edge_differing_mass,
                    "is_truth_class": pair_truth.is_some_and(|(ta, tb)| {
                        (first == ta && second == tb) || (first == tb && second == ta)
                    }),
                })
            })
            .collect();
        serde_json::to_writer(&mut report, &serde_json::json!({
            "locus": locus + locus_offset,
            "physical_rows": candidates.len(),
            "identity_kinds": identity_kinds.iter().map(|(kind, count)| (kind, count))
                .collect::<Vec<_>>(),
            "multi_segment_rows": multi_segment_rows,
            "disjoint_material_rows": disjoint_material_rows,
            "record_count": records.len(),
            "graph_rows": graph_rows,
            "node_universe_count": universe_nodes.len(),
            "edge_universe_count": universe_edges.len(),
            "observed_node_mass": observed_node_mass,
            "observed_edge_mass": observed_edge_mass,
            "observed_node_norm": observed_node_norm,
            "observed_edge_norm": observed_edge_norm,
            "eligible_pairs_nodes": eligible[0],
            "eligible_pairs_combined": eligible[1],
            "competitor_entries": competitor_entries,
            "truth_piece_presence": truth_pieces[locus].iter().map(|v| !v.is_empty())
                .collect::<Vec<_>>(),
            "truth_pair_expressible": pair_truth.is_some(),
            "truth_rows": truth,
            "truth_graph_rows": pair_truth.map(|(a, b)| [a, b]),
            "nodes_truth_cosine": truth_nodes_score,
            "nodes_truth_rank": ranks[0],
            "nodes_truth_tied_pairs": truth_nodes_score.map(|_| counts[0].1),
            "nodes_higher_competitors": truth_nodes_score.map(|_| counts[0].0),
            "nodes_best_cosine": best[0].map(|(score, _)| score),
            "nodes_best_row_indices": best[0].map(|(_, indices)| indices
                .map(|index| space.rows[index].members[0])),
            "combined_truth_cosine": truth_combined_score,
            "combined_truth_rank": ranks[1],
            "combined_truth_tied_pairs": truth_combined_score.map(|_| counts[1].1),
            "combined_higher_competitors": truth_combined_score.map(|_| counts[1].0),
            "combined_best_cosine": best[1].map(|(score, _)| score),
            "combined_best_row_indices": best[1].map(|(_, indices)| indices
                .map(|index| space.rows[index].members[0])),
            "rank_delta_combined_vs_nodes": match (ranks[0], ranks[1]) {
                (Some(a), Some(b)) => Some(b as i64 - a as i64),
                _ => None,
            },
            // Product QUAL block (nodes+edges arm; the stage-1 machinery),
            // CLUSTER FORM. The called set is the bit-identical maximum
            // (class signatures, not member route identities). The cluster
            // cut is the derived knee of the winner's material-distance
            // spectrum (`qual_knee_distance`, null = no knee — a
            // measurement); `qual_cluster_size` counts the classes within
            // the cut (the winner included); k = `qual_cluster_k` DIVERGENT
            // clusters whose best candidates bit-tie the maximum;
            // a = `qual_alternative_similarity` is the best similarity
            // OUTSIDE the near-identical called cluster(s) (null = none,
            // i.e. a = 0); p = s_win/(k*s_win + a),
            // QUAL = -10*log10(1 - p). Unbounded (p = 1) is null plus
            // `qual_unbounded` — never clamped. `qual_similarity_total`
            // remains a diagnostic emission only: the total was the
            // FALSIFIED share form's denominator and no QUAL form uses it.
            // The full sorted material-distance spectrum (with aligned
            // similarities and signature-cosine distances) is the
            // per-locus evidence of where the cut lands; the best divergent
            // rivals are named in `qual_alternative_classes`.
            "qual_similarity_total": (!called.is_empty()).then_some(combined_total),
            "qual_best_similarity": (!called.is_empty()).then_some(called_score),
            "qual_spectrum_shape": cluster_shape,
            "qual_knee_distance": cluster_knee,
            "qual_cluster_size": cluster_size,
            "qual_cluster_k": cluster_k,
            "qual_spectrum_classes": (!called.is_empty()).then_some(spectrum.len()),
            "qual_spectrum_zero_distance_classes": spectrum_zero_distance,
            "qual_nearest_rival_distance": nearest_rival_distance,
            "qual_excluded_class_count": cluster_excluded,
            "qual_distance_spectrum": spectrum.iter().map(|&(distance, ..)| distance)
                .collect::<Vec<_>>(),
            "qual_distance_spectrum_scores": spectrum.iter().map(|&(_, score, ..)| score)
                .collect::<Vec<_>>(),
            "qual_signature_distance_spectrum": spectrum
                .iter()
                .map(|&(_, _, signature, ..)| signature)
                .collect::<Vec<_>>(),
            "qual_signature_view_knee": signature_view.map(|(has_knee, _, _)| has_knee),
            "qual_signature_view_knee_distance": signature_view.map(|(_, cut, _)| cut),
            "qual_signature_view_cluster_size": signature_view.map(|(_, _, size)| size),
            "qual_tied_class_bands": qual_tied_class_bands,
            "qual_alternative_similarity": cluster_alternative,
            "qual_alternative_classes": qual_alternative_classes,
            "qual_delta_similarity": qual_delta,
            "qual_called_class_count": called.len(),
            "qual_called_classes": qual_called_classes,
            "qual_called_physical_pairs": qual_called_physical_pairs,
            "qual_truth_in_called_set": truth_in_called_set,
            "qual_p": confidence,
            "qual_unbounded": qual_unbounded,
            "qual": qual,
            "combined_ties_with_truth_streamed": ties_streamed,
        }))?;
        writeln!(report)?;
    }
    report.flush()?;
    competitors.flush()?;
    ties.flush()
}

/// Total bp of a row's merged per-path material (the path-continuity audit
/// above: a strictly adjacent same-source stitch merges to its full bp
/// count; a gapped or overlapping row does not).
fn row_material_merged_bp(material: &Material) -> u64 {
    material.iter().map(|&(_, lo, hi)| hi - lo).sum()
}

#[cfg(test)]
mod graph_tests {
    use super::{
        contained_steps, multiset_overlap, pack_edge, row_material_merged_bp, GraphRow,
        GraphSpace, GraphUsage,
    };

    fn dense_dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b).map(|(x, y)| x * y).sum()
    }

    /// The factored pair arithmetic must equal the direct dense-vector
    /// cosine over the combined node+edge space (and over the node space
    /// alone for the nodes-only arm) on a hand-checkable example with a
    /// repeat-visit multiplicity and a shared segment.
    #[test]
    fn graph_factored_pair_scores_match_direct_dense_vectors() {
        // Universe: nodes {1,2,3}, edge {(1,2)}; observed mass n1=2,
        // n2=3, e12=1. Row A: nodes {1:1,2:1}, edge {(1,2):1}; row B:
        // nodes {2:2,3:1} (node 2 visited twice), no edges.
        let observed_nodes = [(1u32, 2.0f64), (2, 3.0)];
        let observed_edges = [(pack_edge(1, 2), 1.0f64)];
        let depth = 15.0;
        let row = |nodes: Vec<(u32, u32)>, edges: Vec<((u32, u32), u32)>| GraphRow {
            node_dot: nodes.iter().map(|&(key, mult)| {
                observed_nodes.iter().find(|&&(n, _)| n == key)
                    .map_or(0.0, |&(_, m)| m) * depth * mult as f64
            }).sum(),
            edge_dot: edges.iter().map(|&((a, b), mult)| {
                observed_edges.iter().find(|&&(e, _)| e == pack_edge(a, b))
                    .map_or(0.0, |&(_, m)| m) * depth * mult as f64
            }).sum(),
            node_norm: depth * depth
                * nodes.iter().map(|&(_, mult)| mult as f64 * mult as f64).sum::<f64>(),
            edge_norm: depth * depth
                * edges.iter().map(|&(_, mult)| mult as f64 * mult as f64).sum::<f64>(),
            nodes: nodes.into_iter().map(|(key, mult)| (key as u64, mult)).collect(),
            edges: edges.into_iter().map(|((a, b), mult)| (pack_edge(a, b), mult)).collect(),
            members: vec![],
        };
        let space = GraphSpace {
            rows: vec![
                row(vec![(1, 1), (2, 1)], vec![((1, 2), 1)]),
                row(vec![(2, 2), (3, 1)], vec![]),
            ],
            observed_node_norm: observed_nodes.iter().map(|&(_, m)| m * m).sum(),
            observed_edge_norm: observed_edges.iter().map(|&(_, m)| m * m).sum(),
        };
        // Dense combined space (n1, n2, n3, e12).
        let observed = [2.0, 3.0, 0.0, 1.0];
        let expected_pair = [15.0, 45.0, 15.0, 15.0];
        let direct_combined = dense_dot(&observed, &expected_pair)
            / (dense_dot(&observed, &observed)
                * dense_dot(&expected_pair, &expected_pair)).sqrt();
        assert!((space.cosine(0, 1, depth, true).unwrap() - direct_combined).abs() < 1e-14);
        // Dense node space alone (n1, n2, n3).
        let observed_nodes_dense = [2.0, 3.0, 0.0];
        let expected_nodes_dense = [15.0, 45.0, 15.0];
        let direct_nodes = dense_dot(&observed_nodes_dense, &expected_nodes_dense)
            / (dense_dot(&observed_nodes_dense, &observed_nodes_dense)
                * dense_dot(&expected_nodes_dense, &expected_nodes_dense)).sqrt();
        assert!((space.cosine(0, 1, depth, false).unwrap() - direct_nodes).abs() < 1e-14);
        // The self-pair (diplotype of two identical copies) doubles the
        // expected mass: dense [30,90,30,30] vs the same observation.
        let expected_double = [30.0, 90.0, 30.0, 30.0];
        let direct_double = dense_dot(&observed, &expected_double)
            / (dense_dot(&observed, &observed)
                * dense_dot(&expected_double, &expected_double)).sqrt();
        assert!((space.cosine(0, 0, depth, true).unwrap() - direct_double).abs() < 1e-14);
    }

    #[test]
    fn contained_steps_requires_full_window_containment() {
        // k = 4; steps at bp 0, 4, 8, 12.
        let steps: Vec<(u64, i32)> = vec![(0, 1), (4, 2), (8, 3), (12, 4)];
        // [4, 12): the steps at 4 and 8 fit (8+4<=12), the step at 12 does
        // not (12+4>12), the step at 0 does not (0<4).
        let window = contained_steps(&steps, 4, 4, 12);
        assert_eq!(&steps[window], &[(4u64, 2i32), (8, 3)]);
        // [0, 16): everything.
        let window = contained_steps(&steps, 4, 0, 16);
        assert_eq!(&steps[window], &steps);
        // A window no step fits fails closed.
        let window = contained_steps(&steps, 4, 5, 6);
        assert!(window.is_empty());
    }

    #[test]
    fn multiset_overlap_multiplies_repeat_visits() {
        // Row A visits node 1 twice and node 2 once; row B visits node 1
        // once: overlap 2*1 on node 1, nothing elsewhere.
        let left = vec![(1u64, 2u32), (2, 1)];
        let right = vec![(1u64, 1u32)];
        assert_eq!(multiset_overlap(&left, &right), 2.0);
        assert_eq!(multiset_overlap(&right, &left), 2.0);
        assert_eq!(multiset_overlap(&left, &left), 5.0);
    }

    #[test]
    fn packed_edges_order_by_path_position() {
        // The traversal (7, 9) is a different adjacency than (9, 7).
        assert_ne!(pack_edge(7, 9), pack_edge(9, 7));
        assert_eq!(pack_edge(7, 9), (7u64 << 32) | 9);
    }

    #[test]
    fn graph_usage_rows_coalesce_only_on_exact_key_equality() {
        let twin = GraphUsage {
            nodes: vec![(1, 1), (2, 1)],
            edges: vec![(super::pack_edge(1, 2), 1)],
        };
        let rearranged = GraphUsage {
            nodes: vec![(1, 1), (2, 1)],
            edges: vec![(super::pack_edge(2, 1), 1)],
        };
        // Same segments, different adjacency: NOT coalesced (the edge layer
        // keeps rearranged spellings distinct).
        assert_ne!(twin, rearranged);
        // A repeat visit is a distinct usage vector too.
        let repeated = GraphUsage {
            nodes: vec![(1, 2), (2, 1)],
            edges: vec![(super::pack_edge(1, 2), 1)],
        };
        assert_ne!(twin, repeated);
        assert_eq!(twin, twin.clone());
    }

    /// A clean k-way material-class tie holding ALL the domain's similarity
    /// lands EXACTLY at the derived bound Q = -10*log10(1 - 1/k): with no
    /// alternative mass, p = s_win/(k*s_win) = 1/k bit-identically (the
    /// chosen s values are powers of two, so k*s_win is exact); an
    /// alternative carrying no similarity is the same value; any
    /// competitor mass strictly lowers p and QUAL below the bound.
    #[test]
    fn qual_all_mass_k_way_class_tie_lands_at_derived_bound() {
        for (k, s) in [(2usize, 0.5f64), (3, 0.25), (5, 0.03125), (9, 0.0078125)] {
            let p = super::qual_p(s, k, None).unwrap();
            assert_eq!(p, 1.0 / k as f64);
            let qual = super::qual_from_p(p).unwrap();
            assert_eq!(qual, -10.0 * (1.0 - 1.0 / k as f64).log10());
            // a = 0 as an explicit zero-valued alternative is identical to
            // no alternative at all (reduction iv's other half: classes
            // outside the called set enter only through their similarity).
            assert_eq!(super::qual_p(s, k, Some(0.0)).unwrap(), p);
            // Any competitor mass strictly LOWERS p, so measured Q sits
            // strictly BELOW the bound (the theorem direction).
            let diluted = super::qual_p(s, k, Some(s / 4.0)).unwrap();
            assert!(diluted < 1.0 / k as f64);
            assert!(super::qual_from_p(diluted).unwrap() < qual);
        }
    }

    /// QUAL is strictly monotone in the alternative a and in the called-set
    /// size k (p = s_win/(k*s_win + a) decreases in both); a near-twin
    /// alternative (a = s_win, k = 1) is the honest Q ~ 3.01 and a unique
    /// separated max (a << s_win) scores far above it — the reductions the
    /// falsified share form could not make.
    #[test]
    fn qual_monotonic_in_alternative_and_called_set_size() {
        let s = 0.4f64;
        // Monotone in a at k = 1, falling a raises QUAL.
        let near_tie = super::qual_from_p(super::qual_p(s, 1, Some(s)).unwrap()).unwrap();
        assert!((near_tie - 3.010299956639812).abs() < 1e-9); // -10*log10(1/2)
        let mut previous = near_tie;
        for a in [0.3, 0.2, 0.1, 0.01, 0.001, 0.00001] {
            let q = super::qual_from_p(super::qual_p(s, 1, Some(a)).unwrap()).unwrap();
            assert!(q > previous);
            previous = q;
        }
        // Reduction (i): a unique max with a clear gap (a = s/100) gives
        // p = 100/101, Q ~ 20 — high, versus the near-twin's ~3.01 (iii).
        let separated = super::qual_from_p(super::qual_p(s, 1, Some(s / 100.0)).unwrap()).unwrap();
        assert!((separated - (-10.0 * (1.0 - 100.0 / 101.0_f64).log10())).abs() < 1e-6);
        assert!(separated > 19.9);
        assert!(near_tie < 3.1);
        // Monotone in k at fixed a: a wider bit-tied called set lowers the
        // confidence in the emitted single-material draw (uniform truth
        // within the tied set).
        for k in 1..8usize {
            let wider = super::qual_p(s, k + 1, Some(0.05)).unwrap();
            assert!(wider < super::qual_p(s, k, Some(0.05)).unwrap());
        }
        // A k-way tie with any competitor mass lands strictly below the
        // all-mass bound -10*log10(1 - 1/k).
        let bound = -10.0 * (1.0 - 1.0 / 2.0_f64).log10();
        let diluted = super::qual_from_p(super::qual_p(s, 2, Some(s)).unwrap()).unwrap();
        assert!(diluted < bound);
        // Reduction (iv) is STRUCTURAL: member coalescing never lowers QUAL
        // because the formula's only inputs are (s_win, k, a) — k counts
        // DISTINCT classes and no member sum exists — so a called class
        // holding 12 route members and one holding 1 produce the same p.
        // Nothing to compute: the signature has no multiplicity parameter.
    }

    /// Unbounded (p = 1: k = 1 with no alternative similarity) and massless
    /// (s_win = 0, a = 0) loci fail closed to None — no clamped or invented
    /// value is ever emitted.
    #[test]
    fn qual_unbounded_and_empty_domain_fail_closed() {
        assert_eq!(super::qual_p(1.0, 1, None), Some(1.0));
        assert_eq!(super::qual_p(1.0, 1, Some(0.0)), Some(1.0));
        assert_eq!(super::qual_from_p(1.0), None);
        assert_eq!(super::qual_p(0.0, 3, None), None);
        assert_eq!(super::qual_p(0.0, 1, Some(0.0)), None);
        assert_eq!(super::qual_p(0.5, 0, Some(0.25)), None);
    }

    /// The class-multiplicity and signature-distance helpers: a double-copy
    /// class merges m members into m*(m+1)/2 physical pairs (never m*m),
    /// merged usage doubles self-multiplicity, and the symmetric distance is
    /// (0, 0) only on identical multisets.
    #[test]
    fn class_multiplicity_and_signature_distance_helpers() {
        assert_eq!(super::class_physical_pairs(true, 3, 3), 6);
        assert_eq!(super::class_physical_pairs(false, 3, 5), 15);
        let usage = vec![(1u64, 1u32), (2, 2)];
        assert_eq!(
            super::merged_multiset(&usage, &usage),
            vec![(1u64, 2u32), (2, 4)]
        );
        assert_eq!(
            super::merged_multiset(&usage, &[(1u64, 1u32), (3, 1)]),
            vec![(1u64, 2u32), (2, 2), (3, 1)]
        );
        assert_eq!(super::multiset_distance(&usage, &usage), (0, 0));
        assert_eq!(super::multiset_distance(&usage, &[(1u64, 1u32)]), (1, 2));
        assert_eq!(super::multiset_distance(&usage, &[(1u64, 2u32)]), (2, 3));
        // Observed mass on differing keys: zero when every differing key is
        // unobserved, nonzero otherwise, weighted by the multiplicity delta.
        let observed = |key: u64| if key == 1 { 2.5 } else { 0.0 };
        assert_eq!(
            // keys 2 and 3 differ (multiplicity 2 swapped); key 1 agrees —
            // every DIFFERING key is unobserved, so no read breaks the tie.
            super::differing_observed_mass(&usage, &[(1u64, 1u32), (3, 2)], &observed),
            0.0
        );
        assert_eq!(
            super::differing_observed_mass(&[(1u64, 2u32), (2, 1)], &[(1u64, 1u32), (2, 1)], &observed),
            2.5 // key 1 differs by one visit carrying 2.5 mass
        );
        assert_eq!(super::differing_observed_mass(&usage, &usage, &observed), 0.0);
        let row = |nodes: Vec<(u32, u32)>, edges: Vec<((u32, u32), u32)>| GraphRow {
            node_dot: 0.0,
            edge_dot: 0.0,
            node_norm: 0.0,
            edge_norm: 0.0,
            nodes: nodes.into_iter().map(|(key, mult)| (key as u64, mult)).collect(),
            edges: edges
                .into_iter()
                .map(|((a, b), mult)| (super::pack_edge(a, b), mult))
                .collect(),
            members: vec![],
        };
        let left = row(vec![(1, 1)], vec![((1, 2), 1)]);
        let twin = row(vec![(1, 1)], vec![((1, 2), 1)]);
        let different_edge = row(vec![(1, 1)], vec![((2, 1), 1)]);
        assert_eq!(super::usage_hash(&left), super::usage_hash(&twin));
        assert_ne!(super::usage_hash(&left), super::usage_hash(&different_edge));
    }

    /// The derived knee rule: the maximal relative jump between consecutive
    /// strictly-positive distances, with ties resolved to the largest
    /// index; a jump that does not exceed the median jump is scale-free
    /// and NO knee is claimed; the exactly-zero band never pins the cut.
    #[test]
    fn spectrum_knee_derived_rule() {
        // Dense near-identical band then a divergent tail: the cut is the
        // LAST band point before the maximal relative jump.
        let knee = super::spectrum_knee(&[0.0, 0.0, 0.001, 0.002, 5.0, 8.0, 12.0]);
        assert!(knee.has_knee);
        assert!((knee.cut - 0.002).abs() < 1e-12);
        // The zero band never pins the cut (bit-identical material is
        // innermost, always in the cluster): positives [0.5, 1, 500] give
        // jumps 0.5 and 0.998, so the knee is 1.0, not 0.5.
        let knee = super::spectrum_knee(&[0.0, 0.5, 1.0, 500.0]);
        assert!(knee.has_knee);
        assert!((knee.cut - 1.0).abs() < 1e-12);
        // Exact geometric growth is scale-free: every relative jump equals
        // the median, so no knee exists — a measurement, not an error.
        assert!(!super::spectrum_knee(&[1.0, 10.0, 100.0, 1000.0]).has_knee);
        // Fewer than two strictly-positive distances: no knee claimed.
        assert!(!super::spectrum_knee(&[0.0, 0.0]).has_knee);
        assert!(!super::spectrum_knee(&[0.0, 7.0]).has_knee);
        assert!(!super::spectrum_knee(&[]).has_knee);
        // Max-jump ties resolve to the LARGEST index so the band absorbs
        // the full dense run: [1, 2, 3, 30, 300] gives jumps
        // 0.5, 1/3, 0.9, 0.9 — the tie between 3→30 and 30→300 cuts at 30.
        let knee = super::spectrum_knee(&[1.0, 2.0, 3.0, 30.0, 300.0]);
        assert!(knee.has_knee);
        assert!((knee.cut - 30.0).abs() < 1e-12);
        // A flat all-zero spectrum (every class bit-identical): no knee,
        // cut 0 — the whole domain is one cluster by the distance rule.
        let knee = super::spectrum_knee(&[0.0, 0.0, 0.0]);
        assert!(!knee.has_knee);
        assert_eq!(knee.cut, 0.0);
    }

    /// The signature-cosine distance (the second view): 0 on identical
    /// signatures, 1 on disjoint ones, hand-computed in between; depth
    /// cancels because it is a common factor of both raw signatures.
    #[test]
    fn signature_cosine_distance_second_view() {
        let nodes = vec![(1u64, 1u32), (2, 1)];
        let edges = vec![(pack_edge(1, 2), 1u32)];
        assert_eq!(
            super::signature_cosine_distance(&nodes, &edges, &nodes, &edges),
            0.0
        );
        let other_nodes = vec![(3u64, 1u32)];
        let other_edges = vec![(pack_edge(3, 4), 1u32)];
        assert_eq!(
            super::signature_cosine_distance(&nodes, &edges, &other_nodes, &other_edges),
            1.0
        );
        // a = {1:1, 2:1} + edge(1,2) (||a||² = 3); b = {1:2, 2:1} + edge(1,2)
        // (||b||² = 6); overlap = 1*2 + 1*1 + 1*1 = 4.
        let doubled = vec![(1u64, 2u32), (2, 1)];
        let expected = 1.0 - 4.0 / (3.0f64 * 6.0f64).sqrt();
        assert!(
            (super::signature_cosine_distance(&nodes, &edges, &doubled, &edges) - expected).abs()
                < 1e-15
        );
        assert_eq!(super::signature_cosine_distance(&[], &[], &nodes, &edges), 1.0);
        assert_eq!(super::signature_cosine_distance(&[], &[], &[], &[]), 0.0);
    }

    /// The CLUSTER-FORM reductions, measured on synthetic domains (the
    /// owner's ruling made countable): (i) an all-near-identical domain
    /// with a far divergent tail rates HIGH; (ii) a bit-tie of divergent
    /// materials keeps the honest low bound; (iii) a divergent near-score
    /// rival supplies a and keeps Q low; (iv) bit-identical members are
    /// invisible. Plus the divergent-cluster k semantics (single-linkage
    /// chains) and monotonicity in a at cluster level.
    #[test]
    fn cluster_form_qual_owner_reductions() {
        let s = 0.5f64;
        // (i) The owner's ten-sequences case: nine near-identical twins
        // (small observed mass on differing material) form ONE cluster —
        // the knee absorbs them — and the far rival supplies a, so the
        // call rates HIGH.
        let mut spectrum = Vec::new();
        for _ in 0..9 {
            spectrum.push((0.001, 0.4999));
        }
        spectrum.push((50.0, 0.05));
        spectrum.push((80.0, 0.04));
        let cluster = super::cluster_form_qual(s, &spectrum, &[]);
        assert_eq!(cluster.cluster_size, 10); // the winner + nine twins
        assert_eq!(cluster.k, 1);
        assert_eq!(cluster.alternative, Some(0.05));
        assert!(cluster.p.unwrap() > 0.9); // s/(s + 0.05) = 0.909
        assert_eq!(cluster.shape, "knee");
        // (ii) A bit-tie of DIVERGENT materials: no knee on the two-point
        // spectrum, the tie partner sits outside the winner's cluster, so
        // k = 2 and p lands at the honest 1/2 bound (Q = 3.01).
        let spectrum = vec![(50.0, s)];
        let tied = vec![(0usize, vec![0.0f64])]; // the partner's band: itself
        let tied_refs: Vec<(usize, &[f64])> =
            tied.iter().map(|(i, b)| (*i, b.as_slice())).collect();
        let cluster = super::cluster_form_qual(s, &spectrum, &tied_refs);
        assert_eq!(cluster.cluster_size, 1);
        assert_eq!(cluster.k, 2);
        assert_eq!(cluster.alternative, None);
        assert_eq!(cluster.p, Some(0.5));
        // (iii) A divergent NEAR-SCORE rival: no knee, it stays outside the
        // cluster, supplies a ≈ s_win, and Q stays low.
        let cluster = super::cluster_form_qual(s, &[(50.0, 0.499)], &[]);
        assert_eq!(cluster.cluster_size, 1);
        assert_eq!(cluster.k, 1);
        assert_eq!(cluster.alternative, Some(0.499));
        assert!(cluster.p.unwrap() < 0.51);
        assert_eq!(cluster.shape, "no_knee");
        // (iv) Bit-identical members are invisible: distance 0 always
        // lands inside the cluster (no knee needed), never supplies a,
        // never counts in k.
        let cluster = super::cluster_form_qual(s, &[(0.0, 0.4), (0.0, 0.3), (5.0, 0.25)], &[]);
        assert_eq!(cluster.cluster_size, 3);
        assert_eq!(cluster.k, 1);
        assert_eq!(cluster.alternative, Some(0.25));
        // A bit-tied class INSIDE the winner's band does not count in k; a
        // chained twin (near the winner's twin, far from the winner) joins
        // the SAME cluster by single linkage.
        // spectrum: T1 at 0.5 (in band, knee cut 0.5), T2 at 100 (far),
        // rival at 150; d(T1, T2) = 0.4 chains T2 into the winner cluster.
        let spectrum = vec![(0.5, s), (100.0, s), (150.0, 0.1)];
        let tied = vec![
            (0usize, vec![0.0f64, 0.4, 30.0]), // T1's band
            (1usize, vec![0.4f64, 0.0, 30.0]), // T2's band
        ];
        let tied_refs: Vec<(usize, &[f64])> =
            tied.iter().map(|(i, b)| (*i, b.as_slice())).collect();
        let cluster = super::cluster_form_qual(s, &spectrum, &tied_refs);
        assert_eq!(cluster.k, 1); // winner—T1 in band; T2 chains through T1
        assert_eq!(cluster.cluster_size, 2); // T2 is outside the winner's band
        assert_eq!(cluster.alternative, Some(0.1)); // the rival supplies a
        // Monotonicity in a at cluster level: a closer divergent rival
        // (higher similarity) lowers p.
        let far = super::cluster_form_qual(s, &[(0.001, 0.4999), (50.0, 0.05), (80.0, 0.02)], &[]);
        let near = super::cluster_form_qual(s, &[(0.001, 0.4999), (50.0, 0.4), (80.0, 0.02)], &[]);
        assert!(far.p.unwrap() > near.p.unwrap());
        // A divergent bit-tie with a far tail: k = 2 clusters bit-tie the
        // max and the tail supplies a — p = s/(2s + a) stays low even
        // though the tail is far (the REAL ambiguity is the divergent
        // tie, exactly what must be kept).
        let spectrum = vec![(50.0, s), (200.0, 0.1)];
        let tied = vec![(0usize, vec![0.0f64, 200.0])];
        let tied_refs: Vec<(usize, &[f64])> =
            tied.iter().map(|(i, b)| (*i, b.as_slice())).collect();
        let cluster = super::cluster_form_qual(s, &spectrum, &tied_refs);
        assert_eq!(cluster.k, 2);
        assert_eq!(cluster.alternative, Some(0.1));
        assert!(cluster.p.unwrap() < 0.5);
    }
}
