//! Quarantined, read-only coverage-vector diagnostic (IMPG_COSINE_DIAG_OUTPUT).
//! No score here feeds a sweep, DP, posterior, or product column.
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
