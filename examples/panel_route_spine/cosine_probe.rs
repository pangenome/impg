//! Quarantined, read-only coverage-vector diagnostic (IMPG_COSINE_DIAG_OUTPUT).
//! No score here feeds a sweep, DP, posterior, or product column.
use super::*;
use std::io::{self, BufWriter, Write};

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
    use super::{material_coverage_bins, merge_intervals};
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
    fn cosine_probe_record_placement_intervals_do_not_double_count_subwalks() {
        let mut spans = vec![(11, 14), (10, 12), (20, 22), (14, 16)];
        merge_intervals(&mut spans);
        assert_eq!(spans, vec![(10, 16), (20, 22)]);
    }
}
