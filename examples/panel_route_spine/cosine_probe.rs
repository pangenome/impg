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

/// The called material class's share of the domain's total similarity:
/// the class's OWN similarity value over the summed similarity of every
/// distinct material class, each summed once (deterministic pair-walk
/// order). Member multiplicity never enters — nine route identities over
/// one material class carry ONE class value, not nine. None fails closed
/// when there is no similarity mass at all.
fn qual_share(best: f64, total: f64) -> Option<f64> {
    (total > 0.0).then(|| best / total)
}

/// Phred-scaled QUAL from the class share. None = unbounded (share 1: the
/// called classes hold all the domain's similarity mass); emitted as
/// null, never clamped — no free parameters, no cap.
fn qual_from_share(share: f64) -> Option<f64> {
    (share < 1.0).then(|| -10.0 * (1.0 - share).log10())
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
        // Product QUAL state (material-class semantics). The exhaustive
        // walk below sums every distinct class's combined-arm similarity
        // exactly once and collects the classes at the bit-identical
        // maximum. Exact-tie convention: identical f64 score values — every
        // score is finite, non-negative, NaN-free, and -0.0 cannot arise
        // (non-negative dots over positive norms), so IEEE-754 equality IS
        // bit identity here; NO epsilon constant enters the product call.
        let truth_merged_nodes = pair_truth.map(|(a, b)| {
            merged_multiset(&space.rows[a].nodes, &space.rows[b].nodes)
        });
        let truth_merged_edges = pair_truth.map(|(a, b)| {
            merged_multiset(&space.rows[a].edges, &space.rows[b].edges)
        });
        let mut combined_total = 0.0f64;
        let mut called_score = f64::NEG_INFINITY;
        let mut called: Vec<[usize; 2]> = Vec::new();
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
                // class, summed in the deterministic enumeration order; the
                // called set is the classes at the bit-identical maximum.
                if let Some(score) = combined_score {
                    combined_total += score;
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
        // Product QUAL (material-class semantics): the called set is the
        // distinct material classes at the bit-identical maximum; the
        // share is the called class's own value over the total similarity
        // of every distinct class; QUAL = -10*log10(1 - share). Unbounded
        // (share 1) and empty-domain loci emit null, never a clamp.
        let share = (!called.is_empty())
            .then(|| qual_share(called_score, combined_total))
            .flatten();
        let qual = share.and_then(qual_from_share);
        let qual_unbounded = share.is_some_and(|value| value >= 1.0);
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
            // Product QUAL block (nodes+edges arm; the stage-1 machinery).
            // `qual_similarity_total` sums each distinct material class's
            // similarity exactly once; classes failing the zero-norm guard
            // carry no similarity value and contribute nothing. The called
            // set is the bit-identical maximum (class signatures, not
            // member route identities). Unbounded QUAL (share 1) is null
            // plus `qual_unbounded` — never clamped.
            "qual_similarity_total": (!called.is_empty()).then_some(combined_total),
            "qual_best_similarity": (!called.is_empty()).then_some(called_score),
            "qual_called_class_count": called.len(),
            "qual_called_classes": qual_called_classes,
            "qual_called_physical_pairs": qual_called_physical_pairs,
            "qual_truth_in_called_set": truth_in_called_set,
            "qual_share": share,
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
    /// lands EXACTLY at the derived bound Q = -10*log10(1 - 1/k), and extra
    /// zero-similarity classes do not dilute the share.
    #[test]
    fn qual_all_mass_k_way_class_tie_lands_at_derived_bound() {
        for (k, s) in [(2u64, 0.5f64), (3, 0.25), (5, 0.03125), (9, 0.0078125)] {
            let total = std::iter::repeat(s).take(k as usize).sum::<f64>();
            let share = super::qual_share(s, total).unwrap();
            assert_eq!(share, 1.0 / k as f64);
            let qual = super::qual_from_share(share).unwrap();
            assert_eq!(qual, -10.0 * (1.0 - 1.0 / k as f64).log10());
            // Zero-similarity classes contribute nothing to the total.
            let with_zeros = total + 0.0 + 0.0;
            assert_eq!(super::qual_share(s, with_zeros).unwrap(), share);
            // Any competitor mass strictly LOWERS the share, so measured Q
            // sits strictly BELOW the bound (the theorem direction).
            let diluted = super::qual_share(s, total + s / 4.0).unwrap();
            assert!(diluted < 1.0 / k as f64);
            assert!(super::qual_from_share(diluted).unwrap() < qual);
        }
    }

    /// QUAL rises monotonically with the share; a separated unique top
    /// scores far above a near-tie; the exact k=2 all-mass bound is the
    /// ceiling a 2-way tie can never exceed.
    #[test]
    fn qual_monotonic_in_share_and_separated_beats_near_tie() {
        let near_tie = super::qual_from_share(0.5).unwrap();
        let separated = super::qual_from_share(0.99).unwrap();
        assert!(near_tie < separated);
        assert!((near_tie - 3.010299956639812).abs() < 1e-9); // -10*log10(1/2)
        // share 0.99: exactly 1% of the domain mass sits outside the called
        // class, so Q is a hair above Q20 (1 - 0.99 is not exactly 0.01 in
        // binary; the expected value is derived, not hardcoded).
        assert!((separated - (-10.0 * (1.0 - 0.99f64).log10())).abs() < 1e-12);
        assert!(separated > 19.9);
        for (low, high) in [(0.1f64, 0.2), (0.2, 0.5), (0.5, 0.9), (0.9, 0.999)] {
            assert!(super::qual_from_share(low).unwrap() < super::qual_from_share(high).unwrap());
        }
        // A 2-way tie with any competitor mass lands strictly below the
        // all-mass bound -10*log10(1 - 1/2).
        let bound = -10.0 * (1.0 - 0.5f64).log10();
        let diluted = super::qual_from_share(super::qual_share(0.5, 1.0 + 0.5).unwrap()).unwrap();
        assert!(diluted < bound);
    }

    /// Unbounded (share 1) and no-mass (total 0) loci fail closed to None —
    /// no clamped or invented value is ever emitted.
    #[test]
    fn qual_unbounded_and_empty_domain_fail_closed() {
        assert_eq!(super::qual_share(1.0, 1.0), Some(1.0));
        assert_eq!(super::qual_from_share(1.0), None);
        assert_eq!(super::qual_share(0.0, 0.0), None);
        assert_eq!(super::qual_share(0.5, 0.0), None);
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
}
