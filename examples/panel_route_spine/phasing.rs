//! CORRELATION PHASING (owner-directed exploration, 2026-09): combine the
//! Stage-1 local calls into chromosome fragments using PANEL CO-OCCURRENCE
//! correlations plus the routed boundary MEM evidence — with NO port-word
//! legality requirement on the phasing transitions.
//!
//! Local calling is solved (the spine's Stage 1 exhaustive sweep prefers the
//! truth's SK1 pieces at 9/9 tract loci under the same geometric-q tables);
//! the open problem is combining those local calls into coherent chromosome
//! fragments. The phasing information comes from the panel's other
//! haplotypes — which allele rows co-occur on the same panel chromosome
//! across loci — not from port-word row-chain continuity (which provably
//! breaks at indel-shifted syntenic grids: global SK1 row-chain
//! reachability 0/9).
//!
//! Model semantics (constants discipline: every quantity is data-derived
//! and reported; no tuning constants):
//!
//! - PANEL CO-OCCURRENCE: every panel chromosome is one source lane
//!   end-to-end (`path_of_source` maps each lane to exactly one panel
//!   path), so the co-occurrence statistic between a row pair at adjacent
//!   loci is the public same-source adjacency fact: exit and entry segments
//!   share source and orientation and are order-compatible in the source's
//!   own coordinates (forward: exit.end <= entry.start; reverse:
//!   exit.start >= entry.end). The gap (entry.start - exit.end, forward) is
//!   the indel shift between the two loci's syntenic row grids — reported
//!   per boundary as a distribution, never used as a threshold.
//! - TRANSITIONS: a CO-OCCURRING adjacency is charged the panel-attested
//!   REAL junction — the source's own contiguous material spanning the exit
//!   junction (gap-filled), through the same L149 event-seam machinery
//!   against the two adjacent partitions' routed shares (at gap 0 this is
//!   bit-identical to the port-legal seam). A NON-CO-OCCURRING adjacency
//!   is a NOVEL junction charged the plain juxtaposition seam. The
//!   correlation fact selects which boundary evidence applies; there is
//!   no additive switch penalty (that would be a Li-Stephens recombination
//!   constant, which the constants discipline forbids).
//! - STATES: the Stage-1 local pair choices within the REPORTED margin (the
//!   spine's local seam-swing margin, anchored at the locus's best
//!   port-viable pair; the all-native backbone pair is structurally
//!   retained), as ORDERED pairs — both homolog matchings carried as
//!   distinct states, so a whole donor stretch can sit on one homolog slot.
//! - EXACTNESS: plain min-plus over the ordered pair states with an
//!   admissible suffix-bound prune (per-locus observation floors plus
//!   per-boundary observation floors lower-bound every remaining local pair
//!   loss and transition charge; the native backbone chain under the
//!   phasing transitions is the incumbent). No width, no tie rotations.
//! - DOSAGE/CONFLICT LEGALITY is NOT enforced on phasing transitions; the
//!   selected chains' same-source overlap violations are measured and
//!   reported post hoc.
//! - OPTION A: port-word continuity failures are public structural facts;
//!   overlap-type failures among margin-retained member rows fire partial
//!   row generation (cut at the partner row's endpoint port, so the
//!   partial chains port-legally into it) plus the existing split
//!   machinery at the flanking loci over the involved rows. Gap-type
//!   failures need no partial under this model: the co-occurrence
//!   gap-filled transition carries them, and the domain's bridge rows
//!   remain the port-legal alternative (margin-dead under M1's additive
//!   per-partition charge — reported, not repaired here).
//! - HAPLOID ADMISSIBILITY (owner decision 2026-09-23, after the ploidy
//!   audit): the sample is HAPLOID (chrIII exists ONCE as the hybrid
//!   molecule at d10; truth-assignment carries exactly 17 routes, one per
//!   chromosome). The chain DP therefore runs TWO ploidy tracks and
//!   reports both: the diploid track (two real slots, the previous
//!   behavior) and the haploid track whose SECOND slot is legitimately
//!   EMPTY at every locus — a ploidy statement, not an allele. The empty
//!   slot predicts nothing, pays nothing, explains nothing; slot 1 keeps
//!   the full coverage-completeness obligation, so only scorable classes
//!   are haploid-admissible (the (c) hygiene applies to slot 1 unchanged)
//!   and the (empty, empty) pair is excluded because slot 1 can never be
//!   empty. Ploidy is CONSTANT along a chain: the tracks never mix (a
//!   slot-2 that is real anywhere is real everywhere) — a partial second
//!   copy would re-create the hallucination pressure the decision
//!   removes. The empty slot's local charge is the slot-1 allele's SINGLE
//!   merged-with-empty loss under the same merged-pair convention (counts
//!   add before exposure; the empty profile contributes nothing), and its
//!   boundary charge is zero — the honest haploid exposure at the sample's
//!   realized depth.

use super::{
    build_spine_boundary_draft, class_locus_alleles, class_pair_index, draft_feature_set,
    exhaustive_local_sweep, finalize_spine_boundary, scorable_classes, spine_viability,
    LocusClassing, LocusSweep, SpineBoundary, SpineBoundaryDraft,
};
use super::junction::JunctionSpanIndex;
use crate::*;

// ---------------------------------------------------------------------------
// Co-occurrence (public same-source adjacency).
// ---------------------------------------------------------------------------

/// The public co-occurrence relation between a left-locus row `left` and a
/// right-locus row `right`: `Some(gap)` when the rows share source and
/// orientation and are order-compatible in that source's own coordinate
/// frame (gap = the coordinate distance from the left row's physical exit
/// to the right row's physical entry, >= 0), `Some(negative)` when they
/// share source and orientation but OVERLAP (order-incompatible;
/// -gap = overlap bases), `None` when they cannot co-occur on any panel
/// chromosome (different source, or different orientation).
pub(in super) fn cooccurrence_gap(
    left: &genome::SpanningTraversal,
    right: &genome::SpanningTraversal,
) -> Option<i64> {
    let exit = left.segments.last().expect("nonempty traversal");
    let entry = right.segments.first().expect("nonempty traversal");
    if exit.source != entry.source || exit.reverse != entry.reverse {
        return None;
    }
    if exit.reverse {
        // Reverse traversal: the physical exit is at the exit segment's
        // start (coordinates decrease along the molecule); the next row's
        // physical entry is at its entry segment's end.
        Some(exit.start as i64 - entry.end as i64)
    } else {
        Some(entry.start as i64 - exit.end as i64)
    }
}

/// Per-boundary co-occurrence + port-legality diagnostic over the whole
/// public row grid (no profiles, no observations — pure structure).
#[allow(clippy::too_many_arguments)]
pub(in super) fn cooccurrence_diagnostic(
    boundary: usize,
    ranges_left: &[genome::SpanningTraversal],
    ranges_right: &[genome::SpanningTraversal],
    successors: &[Vec<(usize, f64)>],
    retained_left: &HashSet<usize>,
    retained_right: &HashSet<usize>,
    donor_source: Option<usize>,
) -> serde_json::Value {
    let legal = |left: usize, right: usize| {
        successors[left].iter().any(|&(next, _)| next == right)
    };
    let mut sources_both: HashSet<usize> = HashSet::new();
    let mut cooc_pairs = 0u64;
    let mut gap_zero = 0u64;
    let mut gaps: Vec<i64> = Vec::new();
    let mut overlap_pairs = 0u64;
    let mut overlaps: Vec<i64> = Vec::new();
    let mut cooc_legal = 0u64;
    let mut cooc_illegal = 0u64;
    let mut retained_member_failures: Vec<serde_json::Value> = Vec::new();
    let mut retained_overlap_failures: usize = 0;
    let mut retained_gap_failures: usize = 0;
    let mut donor_rows: Vec<serde_json::Value> = Vec::new();
    for (left, left_row) in ranges_left.iter().enumerate() {
        for (right, right_row) in ranges_right.iter().enumerate() {
            let Some(gap) = cooccurrence_gap(left_row, right_row) else {
                continue;
            };
            sources_both.insert(left_row.segments.last().expect("row").source);
            cooc_pairs += 1;
            if gap >= 0 {
                if gap == 0 {
                    gap_zero += 1;
                }
                gaps.push(gap);
            } else {
                overlap_pairs += 1;
                overlaps.push(-gap);
            }
            if legal(left, right) {
                cooc_legal += 1;
            } else {
                cooc_illegal += 1;
                if retained_left.contains(&left) && retained_right.contains(&right) {
                    if gap >= 0 {
                        retained_gap_failures += 1;
                    } else {
                        retained_overlap_failures += 1;
                        retained_member_failures.push(serde_json::json!({
                            "left": left,
                            "right": right,
                            "source": left_row.segments[0].source,
                            "overlap": -gap,
                        }));
                    }
                }
            }
            if donor_source.is_some_and(|source| {
                left_row.segments[0].source == source
                    || right_row.segments[0].source == source
            }) {
                donor_rows.push(serde_json::json!({
                    "left": left,
                    "left_segments": left_row.segments.iter().map(segment_json).collect::<Vec<_>>(),
                    "right": right,
                    "right_segments": right_row.segments.iter().map(segment_json).collect::<Vec<_>>(),
                    "gap": gap,
                    "port_legal": legal(left, right),
                }));
            }
        }
    }
    let summary = |mut values: Vec<i64>| -> serde_json::Value {
        values.sort_unstable();
        if values.is_empty() {
            return serde_json::Value::Null;
        }
        serde_json::json!({
            "count": values.len(),
            "min": values[0],
            "median": values[values.len() / 2],
            "max": values[values.len() - 1],
        })
    };
    serde_json::json!({
        "boundary": boundary,
        "sources_with_cooccurring_rows": sources_both.len(),
        "cooccurring_pairs": cooc_pairs,
        "gap_zero_pairs": gap_zero,
        "gap_distribution": summary(gaps),
        "overlap_pairs": overlap_pairs,
        "overlap_distribution": summary(overlaps),
        "cooccurring_port_legal": cooc_legal,
        "cooccurring_port_illegal": cooc_illegal,
        "retained_member_overlap_failures": retained_overlap_failures,
        "retained_member_gap_failures": retained_gap_failures,
        "retained_member_overlap_failure_pairs":
            retained_member_failures.iter().take(64).collect::<Vec<_>>(),
        "donor_source_rows": donor_rows,
    })
}

fn segment_json(segment: &SourceRange) -> serde_json::Value {
    serde_json::json!({
        "source": segment.source,
        "start": segment.start,
        "end": segment.end,
        "reverse": segment.reverse,
    })
}

// ---------------------------------------------------------------------------
// Option A: partial rows at structurally detected continuity failures.
// ---------------------------------------------------------------------------

/// The CANONICAL continuity failure per (source, boundary): among the
/// margin-retained single-segment forward row pairs of one source whose
/// grids OVERLAP across the boundary and have no port-legal link, the
/// minimal-overlap pair is the grid's canonical adjacency (the other
/// overlapping occurrence pairs are alternative alignments of the same
/// source's material — their partials would duplicate the canonical fix).
/// Overlap bases sort before row ids so the choice is total and
/// deterministic. Reported per boundary; never a threshold.
pub(in super) fn canonical_overlap_failures(
    boundary: usize,
    ranges_left: &[genome::SpanningTraversal],
    ranges_right: &[genome::SpanningTraversal],
    successors: &[Vec<(usize, f64)>],
    retained_left: &HashSet<usize>,
    retained_right: &HashSet<usize>,
) -> Vec<(usize, usize, usize, i64)> {
    let legal = |left: usize, right: usize| {
        successors[left].iter().any(|&(next, _)| next == right)
    };
    let mut canonical: BTreeMap<usize, (i64, usize, usize)> = BTreeMap::new();
    for (left, left_row) in ranges_left.iter().enumerate() {
        if !retained_left.contains(&left)
            || left_row.segments.len() != 1
            || left_row.segments[0].reverse
        {
            continue;
        }
        for (right, right_row) in ranges_right.iter().enumerate() {
            if !retained_right.contains(&right)
                || right_row.segments.len() != 1
                || right_row.segments[0].reverse
            {
                continue;
            }
            let Some(gap) = cooccurrence_gap(left_row, right_row) else {
                continue;
            };
            if gap >= 0 || legal(left, right) {
                continue;
            }
            let source = left_row.segments[0].source;
            let candidate = (gap, left, right);
            let replace = match canonical.get(&source) {
                Some(&current) => candidate < current,
                None => true,
            };
            if replace {
                canonical.insert(source, candidate);
            }
        }
    }
    canonical
        .into_iter()
        .map(|(source, (overlap, left, right))| (source, left, right, overlap))
        .collect()
}

/// Generate port-word-continuous partial rows for the CANONICAL overlap
/// failures (public structural facts): the exit row cut at the entry row's
/// start — or the entry row cut at the exit row's end — shares the
/// partner's endpoint port exactly, so the partial chains port-legally into
/// it. Cuts require a port at the cut (the row-endpoint ports); rows
/// without one are reported, not cut.
#[allow(clippy::too_many_arguments)]
pub(in super) fn overlap_failure_partials(
    boundary: usize,
    canonical: &[(usize, usize, usize, i64)],
    ranges_left: &[genome::SpanningTraversal],
    ranges_right: &[genome::SpanningTraversal],
    graph: &routes::Graph,
    ports: &mut routes::Ports,
) -> io::Result<(Vec<genome::SpanningTraversal>, Vec<genome::SpanningTraversal>, serde_json::Value)>
{
    let mut left_partials: BTreeMap<String, genome::SpanningTraversal> = BTreeMap::new();
    let mut right_partials: BTreeMap<String, genome::SpanningTraversal> = BTreeMap::new();
    let mut skipped_no_port = 0u64;
    let mut generated = 0u64;
    for &(_, left, right, overlap) in canonical {
        let left_range = &ranges_left[left].segments[0];
        let right_range = &ranges_right[right].segments[0];
        let source = left_range.source;
        let occurrence = boundary
            .checked_mul(10_000_000)
            .and_then(|value| value.checked_add(9_500_000))
            .and_then(|value| value.checked_add(source))
            .ok_or_else(|| invalid("co-occurrence partial occurrence overflow"))?;
        // Left partial: [left.start, right.start) at the left locus — its
        // exit port is exactly the right row's entry port.
        if left_range.start < right_range.start
            && ports.at_cut(graph, source, right_range.start, false)?.is_some()
        {
            let identity = format!(
                "cooc-partial:{}:{}:{}-{}",
                boundary, source, left_range.start, right_range.start
            );
            left_partials.insert(
                identity.clone(),
                genome::SpanningTraversal {
                    partition: left_range.partition,
                    identity,
                    segments: vec![SourceRange {
                        partition: left_range.partition,
                        occurrence,
                        source,
                        start: left_range.start,
                        end: right_range.start,
                        reverse: false,
                    }],
                },
            );
            generated += 1;
        } else {
            skipped_no_port += 1;
        }
        // Right partial: [left.end, right.end) at the right locus — its
        // entry port is exactly the left row's exit port.
        if left_range.end < right_range.end
            && ports.at_cut(graph, source, left_range.end, false)?.is_some()
        {
            let identity = format!(
                "cooc-partial-r:{}:{}:{}-{}",
                boundary, source, left_range.end, right_range.end
            );
            right_partials.insert(
                identity.clone(),
                genome::SpanningTraversal {
                    partition: right_range.partition,
                    identity,
                    segments: vec![SourceRange {
                        partition: right_range.partition,
                        occurrence: occurrence.wrapping_add(1),
                        source,
                        start: left_range.end,
                        end: right_range.end,
                        reverse: false,
                    }],
                },
            );
            generated += 1;
        } else {
            skipped_no_port += 1;
        }
    }
    let stats = serde_json::json!({
        "boundary": boundary,
        "canonical_failures": canonical.len(),
        "generated_partials": generated,
        "skipped_no_endpoint_port": skipped_no_port,
        "left_partial_count": left_partials.len(),
        "right_partial_count": right_partials.len(),
    });
    Ok((
        left_partials.into_values().collect(),
        right_partials.into_values().collect(),
        stats,
    ))
}

// ---------------------------------------------------------------------------
// Margin retention (the spine's reported local seam-swing margin, recomputed
// over the post-augmentation sweep and boundary tables).
// ---------------------------------------------------------------------------

/// The spine's LOCAL retention margin: a pair whose local loss exceeds the
/// locus's best port-viable pair by more than the total seam-evidence swing
/// at its two adjacent boundaries (over viable-viable links) cannot be
/// compensated by any boundary evidence. Data-derived and reported per locus;
/// the all-native backbone pair is structurally retained separately.
pub(in super) fn local_seam_swing_margins(
    sweeps: &[LocusSweep],
    successors: &[Vec<Vec<(usize, f64)>>],
    viable: &[Vec<bool>],
) -> Vec<f64> {
    let locus_count = sweeps.len();
    let extremes: Vec<(f64, f64)> = successors
        .iter()
        .enumerate()
        .map(|(boundary, lists)| {
            let mut min = f64::INFINITY;
            let mut max = f64::NEG_INFINITY;
            for (left, list) in lists.iter().enumerate() {
                if !viable[boundary][left] {
                    continue;
                }
                for &(right, score) in list {
                    if viable[boundary + 1][right] {
                        min = min.min(score);
                        max = max.max(score);
                    }
                }
            }
            (min, max)
        })
        .collect();
    (0..locus_count)
        .map(|locus| {
            let mut swing = 0.0f64;
            if locus > 0 {
                swing += (extremes[locus - 1].1 - extremes[locus - 1].0).max(0.0);
            }
            if locus + 1 < locus_count {
                swing += (extremes[locus].1 - extremes[locus].0).max(0.0);
            }
            sweeps[locus].min_viable_loss + swing
        })
        .collect()
}

/// Margin-retained class pairs at one locus (the Stage-1 reported margin;
/// the backbone class pair is structurally retained).
pub(in super) fn retained_class_pairs(
    sweep: &LocusSweep,
    classes: usize,
    backbone_class: usize,
    margin: f64,
) -> Vec<[usize; 2]> {
    let mut retained = Vec::new();
    for first in 0..classes {
        for second in first..classes {
            let loss = sweep.table[class_pair_index(first, second)];
            if loss <= margin || (first == backbone_class && second == backbone_class) {
                retained.push([first, second]);
            }
        }
    }
    retained
}

/// The retained member alleles at one locus (members of any margin-retained
/// class pair — the phasing-relevant allele universe).
pub(in super) fn retained_members(
    _locus: usize,
    locus_classes: &LocusClassing,
    sweep: &LocusSweep,
    membership: &[usize],
    backbone_allele: usize,
    margin: f64,
) -> HashSet<usize> {
    let classes = locus_classes.profiles.len();
    let backbone_class = membership[backbone_allele];
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); classes];
    for (allele, &class) in membership.iter().enumerate() {
        members[class].push(allele);
    }
    let mut retained = HashSet::new();
    for [first, second] in retained_class_pairs(sweep, classes, backbone_class, margin) {
        for class in [first, second] {
            for &allele in &members[class] {
                retained.insert(allele);
            }
        }
    }
    retained
}

/// Ordered pair states at one locus: every ordered physical allele pair of
/// the margin-retained class pairs (both homolog matchings carried), with
/// the pair's local class-table loss.
pub(in super) fn ordered_pair_states(
    _locus: usize,
    locus_classes: &LocusClassing,
    sweep: &LocusSweep,
    membership: &[usize],
    backbone_allele: usize,
    margin: f64,
) -> Vec<([usize; 2], f64)> {
    let classes = locus_classes.profiles.len();
    let backbone_class = membership[backbone_allele];
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); classes];
    for (allele, &class) in membership.iter().enumerate() {
        members[class].push(allele);
    }
    let mut states: Vec<([usize; 2], f64)> = Vec::new();
    for [first, second] in retained_class_pairs(sweep, classes, backbone_class, margin) {
        let loss = sweep.table[class_pair_index(first, second)];
        if first == second {
            let list = &members[first];
            for (offset, &left) in list.iter().enumerate() {
                for &right in &list[offset..] {
                    states.push(([left, right], loss));
                    if left != right {
                        states.push(([right, left], loss));
                    }
                }
            }
        } else {
            for &left in &members[first] {
                for &right in &members[second] {
                    states.push(([left, right], loss));
                    states.push(([right, left], loss));
                }
            }
        }
    }
    states
}

/// The EMPTY second-slot sentinel (haploid admissibility, owner decision
/// 2026-09-23): `route[locus][1] == EMPTY_SLOT2` means the chain's second
/// homolog slot is a legitimate ploidy statement — no allele, no profile,
/// no transitions, no charges. Never a `ranges` index.
pub(in super) const EMPTY_SLOT2: usize = usize::MAX;

/// The haploid track's per-window losses under the DIRECT record-once
/// site-map arithmetic (owner ruling, dp-evaluator-consistency, 2026-09-24):
/// ONE loss per retained ALLELE, evaluated by the SAME evaluator the
/// decompose implements — the oracle spelled allele profile (the
/// depth-calibrated spelled-occurrence counts, not the geometric sweep's
/// window-incidence multiplicity: the geometric signal scale pays
/// ~O(window) pure signal cost per feature regardless of the observed
/// share, which is what made short fragments win tract windows on signal
/// cost alone) charged by the decompose's exact helper
/// (`merged_single_loss_sample_with_omission`: the oracle profile against
/// the row's own owner-resolved record-once observed side — bit-identical
/// to `SiteObserved::site_map` at a single owner — plus the omission addend
/// over the window's in-window observed profile). The per-class restricted
/// interior-junction charge stays folded in (classes are charge-homogeneous,
/// so the class charge IS the allele's charge). Non-scorable alleles (the
/// empty-profile / sub-read-length hygiene (c)) stay INFINITY.
/// Computed BEFORE the margin retention so the retained candidate sets
/// reflect the new objective (a chain the retention excluded can never be
/// selected). Measured cost is reported in the run's stage walls.
/// One allele's exonation-channel inputs for the pairwise-coupling search
/// objective (supervisor ruling on genome/finalist-reranking): the allele's
/// oracle profile intersected with a window's observed-feature universe, as
/// SORTED INDICES into that window's observed feature list (the weighted
/// sums recover the per-window omission terms by index — one term table per
/// window, the same values the omission charge uses). `self` = this allele's
/// own window; `left`/`right` = the ADJACENT windows' (the pairwise
/// approximation of the cross-window exoneration channel).
pub(in super) struct HaploidAlleleExon {
    pub(in super) self_list: Vec<u32>,
    pub(in super) self_sum: f64,
    pub(in super) left: Option<(Vec<u32>, f64)>,
    pub(in super) right: Option<(Vec<u32>, f64)>,
    /// The INSTANCE-LEVEL covered sets (the corrected Ruling 2's coupling
    /// channel): per observed feature (the window's obs index), the covered
    /// instance ids (sorted, entry order of the window's instance list) and
    /// their mass. Empty maps = the feature-level fallback.
    pub(in super) self_instances: AlleleCoveredInstances,
    pub(in super) left_instances: Option<AlleleCoveredInstances>,
    pub(in super) right_instances: Option<AlleleCoveredInstances>,
}

/// Per (allele, window): the allele's instance-level covered sets, SORTED
/// by feature index (the pair arithmetic's merge join works on the sorted
/// vectors — the hash-loop form measured 1710s vs the legacy form's 59s).
/// `solo_credits[f]` = term(full_mass) - term(full_mass - masses[f]): the
/// exonation credit the allele alone gives a self row that covers nothing
/// of the feature (the pair-only overlap branch re-prices the remainder).
#[derive(Default)]
pub(in super) struct AlleleCoveredInstances {
    pub(in super) features: Vec<u32>,
    pub(in super) masses: Vec<f64>,
    pub(in super) solo_credits: Vec<f64>,
    pub(in super) ids: Vec<Vec<u32>>,
}

impl AlleleCoveredInstances {
    fn is_empty(&self) -> bool {
        self.features.is_empty()
    }
}

pub(in super) fn haploid_allele_losses(
    // The COMBINED per-locus allele list: the domain rows (indices
    // 0..domain_alleles, the classing/sweeps universe — charged exactly as
    // before, bit-identically) followed by the same-owner stitched spanning
    // candidates (indices domain_alleles.., charged by the owner-set
    // record-once form below).
    ranges: &[genome::SpanningTraversal],
    domain_alleles: usize,
    locus_classes: &LocusClassing,
    window_obs_map: &HashMap<FeatureKey, f64>,
    // The windows' observed-feature indices (this locus = window_obs_map's;
    // the neighbors' for the exonation channel).
    obs_index: &[crate::WindowObsIndex],
    // The WHOLE slice's observed maps (the instance-level covers' solo
    // credits price the TARGET window's own full masses).
    window_obs_all: &[HashMap<FeatureKey, f64>],
    locus: usize,
    routed_equal: &crate::RoutedObs,
    site: &crate::SiteObserved,
    component_locus_to_partition: &[u32],
    sample: &SampleSideBackgrounds,
    scorable: &[bool],
    model: &ScoreModel,
    panel: &SyngIndex,
    sources: &routes::Sources,
    // The instance-level exoneration's observed side (the whole slice;
    // this locus and its two neighbors are read) and the sources -> panel
    // path map (the spelled spans' coordinates).
    window_instances: Option<&crate::InstanceStructure>,
    path_of_source: &[usize],
) -> io::Result<(Vec<f64>, Vec<HaploidAlleleExon>)> {
    let alleles = ranges.len();
    ensure(
        locus_classes.membership.len() == domain_alleles,
        "haploid allele-loss membership cardinality mismatch",
    )?;
    ensure(
        scorable.len() == locus_classes.profiles.len(),
        "haploid scorable cardinality mismatch",
    )?;
    ensure(
        locus_classes.class_charges.len() == locus_classes.profiles.len(),
        "haploid class charge cardinality mismatch",
    )?;
    let mut oracle_memo: HashMap<String, Profile> = HashMap::new();
    // The instance-level exoneration's per-allele spelled spans (positions
    // only; memoized per traversal identity within this locus's pass — the
    // same memo granularity the oracle profiles keep).
    let mut allele_spans_memo: HashMap<String, HashMap<FeatureKey, Vec<crate::SpelledSpan>>> =
        HashMap::new();
    let mut piece_spans_memo: crate::SpelledMemo = HashMap::new();
    // The allele's instance-level covered sets at one window: per observed
    // feature, the covered instance ids (the window's instance entries'
    // order) and their mass. Coverage = a spelled span overlapping the
    // instance's placement span on the same panel path (the record-atomic
    // form: a record whose ANY in-window placement is covered is covered —
    // the established record-once share semantics the observed maps charge
    // with).
    let exon_instances = |spans: &HashMap<FeatureKey, Vec<crate::SpelledSpan>>,
                          window: usize|
     -> AlleleCoveredInstances {
        let obs = &obs_index[window];
        let mut scratch: Vec<(u32, Vec<u32>, f64)> = Vec::new();
        let (win_records, record_spans, record_shares) = match window_instances {
            Some(instances) => (
                instances.window_records.get(window),
                &instances.record_spans,
                &instances.record_shares,
            ),
            None => return AlleleCoveredInstances::default(),
        };
        for (feature, spans_f) in spans {
            let Some(&idx) = obs.index.get(feature) else {
                continue;
            };
            let Some(records) = win_records.and_then(|w| w.get(feature)) else {
                continue;
            };
            let mut ids: Vec<u32> = Vec::new();
            let mut mass = 0.0f64;
            for (id, record_index) in records.iter().enumerate() {
                let record_spans_f = record_spans
                    .get(*record_index as usize)
                    .and_then(|m| m.get(feature))
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]);
                if crate::record_covered(record_spans_f, spans_f) {
                    ids.push(id as u32);
                    mass += record_shares
                        .get(*record_index as usize)
                        .copied()
                        .unwrap_or(0.0);
                }
            }
            if !ids.is_empty() {
                scratch.push((idx, ids, mass));
            }
        }
        // Ascending feature-index order: the pair arithmetic's merge join
        // and the f64 sums' determinism both require it.
        scratch.sort_by_key(|(idx, _, _)| *idx);
        let features: Vec<u32> = scratch.iter().map(|(idx, _, _)| *idx).collect();
        let masses: Vec<f64> = scratch.iter().map(|(_, _, mass)| *mass).collect();
        let solo_credits: Vec<f64> = scratch
            .iter()
            .map(|(idx, _, mass)| {
                let full = window_obs_all[window][&obs.features[*idx as usize]];
                crate::omission_term_value(&obs.features[*idx as usize], full, sample)
                    - crate::omission_term_value(
                        &obs.features[*idx as usize],
                        (full - mass).max(0.0),
                        sample,
                    )
            })
            .collect();
        let ids: Vec<Vec<u32>> = scratch.into_iter().map(|(_, ids, _)| ids).collect();
        AlleleCoveredInstances {
            features,
            masses,
            solo_credits,
            ids,
        }
    };
    // The stitched alleles' multi-owner observed sides, cached per owner set
    // (the record-once site map over the union of the segments' owning
    // partitions — the established Fix-1 multi-owner form, the same
    // convention the decompose charges the truth's multi-owner rows with).
    let mut owner_set_memo: HashMap<Vec<u32>, HashMap<FeatureKey, f64>> = HashMap::new();
    let mut losses = vec![f64::INFINITY; alleles];
    let exon = |profile: &Profile, window: usize| -> (Vec<u32>, f64) {
        let obs = &obs_index[window];
        let mut list: Vec<u32> = Vec::new();
        for key in profile.keys() {
            if let Some(&i) = obs.index.get(key) {
                list.push(i);
            }
        }
        list.sort_unstable();
        // Ascending index order: the sum's f64 order is the list's order
        // (deterministic — hash iteration never enters a sum).
        let sum = list
            .iter()
            .map(|&i| obs.terms[i as usize])
            .sum::<f64>();
        (list, sum)
    };
    let mut exons: Vec<HaploidAlleleExon> = Vec::new();
    for allele in 0..alleles {
        let mut allele_profile: Option<Profile> = None;
        let loss = if allele < domain_alleles {
            let class = locus_classes.membership[allele];
            if !scorable[class] {
                exons.push(HaploidAlleleExon {
                    self_list: Vec::new(),
                    self_sum: 0.0,
                    left: None,
                    right: None,
                    self_instances: AlleleCoveredInstances::default(),
                    left_instances: None,
                    right_instances: None,
                });
                continue;
            }
            let owner = locus_classes.class_owners[class];
            let owner_obs =
                crate::owner_routed_obs(routed_equal, owner, component_locus_to_partition);
            let profile = crate::oracle_allele_profile(
                panel,
                sources,
                &ranges[allele],
                &mut oracle_memo,
            )?;
            let loss = crate::merged_single_loss_sample_with_omission(
                &profile,
                owner_obs,
                window_obs_map,
                sample,
                model,
            )? + locus_classes.class_charges[class];
            allele_profile = Some(profile);
            loss
        } else {
            // A stitched spanning candidate: the oracle profile of the whole
            // chain (segment profiles plus the interior co-occurring seams —
            // the established same-owner convention, identical to how the
            // truth's identical pieces are charged), charged against the
            // record-once map over the union of the segments' owning
            // partitions, plus the per-window omission addend over the
            // union's unpredicted remainder. No restricted interior charge:
            // every admitted chain's seams are gap-0 co-occurring (the
            // pooled form is merged into the profile). The (c) hygiene
            // applies unchanged: an empty profile or a sub-read-length chain
            // stays INFINITY.
            let traversal = &ranges[allele];
            let total_length: u64 = traversal
                .segments
                .iter()
                .map(|segment| segment.end.saturating_sub(segment.start))
                .sum();
            if total_length < crate::READ_LENGTH as u64 {
                exons.push(HaploidAlleleExon {
                    self_list: Vec::new(),
                    self_sum: 0.0,
                    left: None,
                    right: None,
                    self_instances: AlleleCoveredInstances::default(),
                    left_instances: None,
                    right_instances: None,
                });
                continue;
            }
            let profile =
                crate::oracle_allele_profile(panel, sources, traversal, &mut oracle_memo)?;
            if profile.is_empty() {
                exons.push(HaploidAlleleExon {
                    self_list: Vec::new(),
                    self_sum: 0.0,
                    left: None,
                    right: None,
                    self_instances: AlleleCoveredInstances::default(),
                    left_instances: None,
                    right_instances: None,
                });
                continue;
            }
            let mut owner_set: Vec<u32> = Vec::new();
            for segment in &traversal.segments {
                let owner = segment.partition as u32;
                if !owner_set.contains(&owner) {
                    owner_set.push(owner);
                }
            }
            let owner_obs: &HashMap<FeatureKey, f64> = if owner_set.len() == 1 {
                crate::owner_routed_obs(routed_equal, owner_set[0], component_locus_to_partition)
            } else {
                let key: Vec<u32> = owner_set
                    .iter()
                    .map(|&owner| {
                        crate::owner_universe_partition(owner, component_locus_to_partition)
                    })
                    .collect();
                owner_set_memo
                    .entry(key.clone())
                    .or_insert_with(|| site.site_map(key.iter().copied()))
            };
            let loss = crate::merged_single_loss_sample_with_omission(
                &profile,
                owner_obs,
                window_obs_map,
                sample,
                model,
            )?;
            allele_profile = Some(profile);
            loss
        };
        // The exonation channel's inputs (computed only for scorable
        // alleles; the loss is finite there by construction).
        let profile = allele_profile.expect("scorable allele carries its profile");
        let (self_list, self_sum) = exon(&profile, locus);
        // The instance-level covered sets (the corrected Ruling 2's
        // coupling channel): the allele's spelled anchors against each
        // window's observed instance entries.
        let allele_spans = crate::oracle_allele_spans(
            panel,
            sources,
            path_of_source,
            &ranges[allele],
            &mut allele_spans_memo,
            &mut piece_spans_memo,
        )?;
        let self_instances = exon_instances(&allele_spans, locus);
        let left_instances = (locus > 0).then(|| exon_instances(&allele_spans, locus - 1));
        let right_instances =
            (locus + 1 < obs_index.len()).then(|| exon_instances(&allele_spans, locus + 1));
        exons.push(HaploidAlleleExon {
            self_list,
            self_sum,
            left: (locus > 0).then(|| exon(&profile, locus - 1)),
            right: (locus + 1 < obs_index.len()).then(|| exon(&profile, locus + 1)),
            self_instances,
            left_instances,
            right_instances,
        });
        losses[allele] = loss;
    }
    Ok((losses, exons))
}

/// Haploid-track states at one locus: every margin-retained member allele
/// with an EMPTY second slot, charged by its OWN oracle record-once loss
/// (see `haploid_allele_losses`; the SAMPLE-SIDE cross-support backgrounds,
/// owner-approved family, genome-wide form 2026-09-24: beta_f = base +
/// Sum_r m_r*(1 - 1/t_r_genome)). Only scorable alleles are
/// haploid-admissible (the (c) hygiene applies to slot 1 unchanged: an
/// empty-profile or sub-read-length mini must not become a zero-cost
/// "haploid" winner — the empty slot is a ploidy statement, not an
/// allele).
pub(in super) fn ordered_haploid_states(
    // The per-allele oracle record-once losses (precomputed before the
    // retention; see `haploid_allele_losses`).
    haploid_loss: &[f64],
    retained: &HashSet<usize>,
) -> io::Result<Vec<([usize; 2], f64)>> {
    let mut states = Vec::new();
    for &allele in retained.iter() {
        let loss = haploid_loss[allele];
        if loss.is_finite() {
            states.push(([allele, EMPTY_SLOT2], loss));
        }
    }
    states.sort_by(|left, right| left.0.cmp(&right.0).then_with(|| {
        left.1.total_cmp(&right.1)
    }));
    Ok(states)
}

// ---------------------------------------------------------------------------
// Boundary transition costs (co-occurrence gap-filled vs novel junctions).
// ---------------------------------------------------------------------------

pub(in super) struct BoundaryTransitionCosts {
    /// Dense [retained_left x retained_right] single-homolog transition
    /// charge (u32::MAX dense index = not a retained member) — the DIPLOID
    /// track's matrix (panel-multiplicity backgrounds on co-occurring
    /// compositions; restricted charges are background-invariant).
    pub(in super) cost: Vec<f64>,
    /// The HAPLOID track's matrix over the same dense layout: co-occurring
    /// compositions charged with the SAMPLE-SIDE cross-support backgrounds
    /// (the haploid single-allele convention — Step B, owner-approved
    /// 2026-09-23); novel (restricted) entries are identical to `cost`.
    pub(in super) cost_haploid: Vec<f64>,
    pub(in super) left_index: Vec<u32>,
    pub(in super) right_index: Vec<u32>,
    pub(in super) right_count: u32,
    pub(in super) stats: serde_json::Value,
}

/// The transition charge between every retained-member allele pair at a
/// boundary: CO-OCCURRING adjacencies are charged the panel-attested real
/// (gap-filled) junction through the source's own material against the
/// adjacent partitions' pooled shares (a gap-0 co-occurring pair scores
/// bit-identically to the port-legal seam of the old chain model); all
/// other pairs are NOVEL junctions, paid by the reads that actually cross
/// them (the within-read adjacency restriction): the juxtaposition seam
/// profile charged against the crossing reads' equal-share routed mass at
/// the two adjacent partitions. Charges pool by seam composition (the seam
/// machinery's own junction granularity): a composition realized by any
/// co-occurring pair is POOLED (protecting the real junctions' bit-identity
/// — and the composition-granularity limit is reported as a measured
/// statistic); novel-only compositions are restricted.
/// Draft of one boundary's transition compositions: the retained-pair
/// composition enumeration, the composition seam profiles, and the
/// NOVEL-only restricted charges (the junction machinery — flat background,
/// unchanged). The POOLED composition losses are deferred to
/// `finalize_boundary_transition_costs`: their features must be in the
/// multiplicity-background scan before they are charged.
pub(in super) struct BoundaryTransitionDraft {
    pair_composition: Vec<u32>,
    pair_kind: Vec<u8>,
    composition_pooled: Vec<bool>,
    /// Per composition: the DISTINCT owning universe partitions of its
    /// realizing allele pairs (the owner-resolved observed sides).
    composition_owners: Vec<BTreeSet<u32>>,
    profiles: Vec<Profile>,
    novel_scores: Vec<(f64, u64)>,
    left_count: usize,
    right_count: usize,
    left_index: Vec<u32>,
    right_index: Vec<u32>,
    cooc_pairs: u64,
    cooc_gap_zero: u64,
    cooc_gap_positive: u64,
    novel_pairs: u64,
    overlap_novel: u64,
    empty_tail_pairs: u64,
    novel_riding_pooled: u64,
    novel_only: Vec<usize>,
}

impl BoundaryTransitionDraft {
    /// Every feature the boundary's compositions can charge (the pooled
    /// composition losses' scan-universe contribution).
    pub(in super) fn feature_set(&self) -> BTreeSet<FeatureKey> {
        let mut set = BTreeSet::new();
        for profile in &self.profiles {
            for key in profile.keys() {
                set.insert(key.clone());
            }
        }
        set
    }
}

#[allow(clippy::too_many_arguments)]
pub(in super) fn build_boundary_transition_draft(
    panel: &SyngIndex,
    sources: &routes::Sources,
    flank_memo: &FlankMemo,
    model: &ScoreModel,
    span: &JunctionSpanIndex,
    path_of_source: &[usize],
    // Per-allele OWNING universe partitions (the window-domain extension's
    // owner-resolved charging).
    left_owners: &[u32],
    right_owners: &[u32],
    component_locus_to_partition: &[u32],
    ranges_left: &[genome::SpanningTraversal],
    ranges_right: &[genome::SpanningTraversal],
    retained_left: &HashSet<usize>,
    retained_right: &HashSet<usize>,
) -> io::Result<BoundaryTransitionDraft> {
    let flank = READ_LENGTH - 1;
    let sorted_left: Vec<usize> = {
        let mut list: Vec<usize> = retained_left.iter().copied().collect();
        list.sort_unstable();
        list
    };
    let sorted_right: Vec<usize> = {
        let mut list: Vec<usize> = retained_right.iter().copied().collect();
        list.sort_unstable();
        list
    };
    let mut left_index = vec![u32::MAX; ranges_left.len()];
    let mut right_index = vec![u32::MAX; ranges_right.len()];
    for (dense, &allele) in sorted_left.iter().enumerate() {
        left_index[allele] = dense as u32;
    }
    for (dense, &allele) in sorted_right.iter().enumerate() {
        right_index[allele] = dense as u32;
    }

    // Distinct oriented flank byte strings, interned once (left = the
    // molecule-true exit tail; right = the novel entry head OR the
    // co-occurring source continuation across the junction).
    let mut left_flanks: Vec<Vec<u8>> = Vec::new();
    let mut right_flanks: Vec<Vec<u8>> = Vec::new();
    let mut left_flank_index: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut right_flank_index: HashMap<Vec<u8>, u32> = HashMap::new();
    let intern_left = |bytes: Vec<u8>,
                           left_flanks: &mut Vec<Vec<u8>>,
                           index: &mut HashMap<Vec<u8>, u32>|
     -> u32 {
        if let Some(&id) = index.get(&bytes) {
            return id;
        }
        let id = left_flanks.len() as u32;
        left_flanks.push(bytes.clone());
        index.insert(bytes, id);
        id
    };
    let intern_right = |bytes: Vec<u8>,
                            right_flanks: &mut Vec<Vec<u8>>,
                            index: &mut HashMap<Vec<u8>, u32>|
     -> u32 {
        if let Some(&id) = index.get(&bytes) {
            return id;
        }
        let id = right_flanks.len() as u32;
        right_flanks.push(bytes.clone());
        index.insert(bytes, id);
        id
    };
    let mut tail_ids: Vec<u32> = vec![u32::MAX; ranges_left.len()];
    let mut head_ids: Vec<u32> = vec![u32::MAX; ranges_right.len()];
    for &allele in &sorted_left {
        let (_, tail) = allele_endpoints(sources, flank_memo, &ranges_left[allele].segments, flank)?;
        tail_ids[allele] = intern_left(tail, &mut left_flanks, &mut left_flank_index);
    }
    for &allele in &sorted_right {
        let (head, _) = allele_endpoints(sources, flank_memo, &ranges_right[allele].segments, flank)?;
        head_ids[allele] = intern_right(head, &mut right_flanks, &mut right_flank_index);
    }

    // Compositions per retained pair: (left flank id, right flank id).
    let mut composition_owners: Vec<BTreeSet<u32>> = Vec::new();
    let mut composition_ids: HashMap<(u32, u32), u32> = HashMap::new();
    let mut compositions: Vec<(u32, u32)> = Vec::new();
    let mut pair_composition: Vec<u32> =
        vec![u32::MAX; sorted_left.len() * sorted_right.len()];
    // Dense-pair evidence kinds: 1 = co-occurring (gap-filled real junction),
    // 0 = novel junction (including same-source overlaps).
    let mut pair_kind: Vec<u8> = vec![0; sorted_left.len() * sorted_right.len()];
    // Per composition: pooled flag (some realizing pair is co-occurring)
    // and the realizing dense pairs (for the restricted charge's event
    // union over novel-only compositions).
    let mut composition_pooled: Vec<bool> = Vec::new();
    let mut composition_pairs: Vec<Vec<(u32, u32)>> = Vec::new();
    let mut continuation_ids: HashMap<(usize, u64, bool), u32> = HashMap::new();
    let mut cooc_pairs = 0u64;
    let mut cooc_gap_zero = 0u64;
    let mut cooc_gap_positive = 0u64;
    let mut novel_pairs = 0u64;
    let mut overlap_novel = 0u64;
    let mut empty_tail_pairs = 0u64;
    for (dense_left, &left) in sorted_left.iter().enumerate() {
        for (dense_right, &right) in sorted_right.iter().enumerate() {
            let gap = cooccurrence_gap(&ranges_left[left], &ranges_right[right]);
            let (left_id, right_id) = match gap {
                Some(gap) if gap >= 0 => {
                    cooc_pairs += 1;
                    if gap == 0 {
                        cooc_gap_zero += 1;
                    } else {
                        cooc_gap_positive += 1;
                    }
                    let exit = ranges_left[left]
                        .segments
                        .last()
                        .expect("nonempty traversal");
                    let (junction, forward) = if exit.reverse {
                        (exit.start, false)
                    } else {
                        (exit.end, true)
                    };
                    let lane_length = sources.lanes[exit.source].1;
                    let bytes = if forward {
                        if junction < lane_length {
                            sources.fetch(
                                exit.source,
                                junction,
                                (junction + flank as u64).min(lane_length),
                            )?
                        } else {
                            Vec::new()
                        }
                    } else {
                        let lo = junction.saturating_sub(flank as u64);
                        if lo < junction {
                            impg::graph::reverse_complement(&sources.fetch(
                                exit.source, lo, junction,
                            )?)
                        } else {
                            Vec::new()
                        }
                    };
                    let right_id = match continuation_ids.get(&(exit.source, junction, forward)) {
                        Some(&id) => id,
                        None => {
                            let id = intern_right(
                                bytes,
                                &mut right_flanks,
                                &mut right_flank_index,
                            );
                            continuation_ids.insert((exit.source, junction, forward), id);
                            id
                        }
                    };
                    (tail_ids[left], right_id)
                }
                _ => {
                    novel_pairs += 1;
                    if gap.is_some() {
                        overlap_novel += 1;
                    }
                    (tail_ids[left], head_ids[right])
                }
            };
            if left_flanks[left_id as usize].is_empty() {
                empty_tail_pairs += 1;
            }
            let composition = match composition_ids.get(&(left_id, right_id)) {
                Some(&id) => id,
                None => {
                    let id = compositions.len() as u32;
                    compositions.push((left_id, right_id));
                    composition_ids.insert((left_id, right_id), id);
                    composition_pooled.push(false);
                    composition_pairs.push(Vec::new());
                    composition_owners.push(BTreeSet::new());
                    id
                }
            };
            let cooccurring = gap.is_some_and(|value| value >= 0);
            if cooccurring {
                composition_pooled[composition as usize] = true;
            } else {
                composition_pairs[composition as usize].push((dense_left as u32, dense_right as u32));
            }
            composition_owners[composition as usize].insert(left_owners[left]);
            composition_owners[composition as usize].insert(right_owners[right]);
            pair_composition[dense_left * sorted_right.len() + dense_right] = composition;
            pair_kind[dense_left * sorted_right.len() + dense_right] = u8::from(cooccurring);
        }
    }
    let composition_count = compositions.len();
    let novel_only: Vec<usize> = (0..composition_count)
        .filter(|&id| !composition_pooled[id])
        .collect();
    let novel_riding_pooled: u64 = pair_kind
        .iter()
        .zip(pair_composition.iter())
        .filter(|(&kind, &composition)| {
            kind == 0 && composition_pooled[composition as usize]
        })
        .count() as u64;
    let composition_profiles: Vec<Profile> = (0..composition_count)
        .into_par_iter()
        .map(|id| {
            let (left_id, right_id) = compositions[id];
            let (profile, _) = genome::profile_event_seam(
                panel,
                &left_flanks[left_id as usize],
                &right_flanks[right_id as usize],
                READ_LENGTH,
                MAX_FEATURES,
            )?;
            Ok(profile)
        })
        .collect::<io::Result<_>>()?;
    let novel_scores: Vec<(f64, u64)> = (0..composition_count)
        .into_par_iter()
        .map(|id| {
            if composition_pooled[id] {
                return Ok((0.0, 0));
            }
            let profile = &composition_profiles[id];
            if profile.is_empty() {
                return Ok((0.0, 0));
            }
            let pairs: Vec<(SourceRange, SourceRange)> = composition_pairs[id]
                .iter()
                .map(|&(dense_left, dense_right)| {
                    (
                        ranges_left[sorted_left[dense_left as usize]]
                            .segments
                            .last()
                            .expect("nonempty traversal")
                            .clone(),
                        ranges_right[sorted_right[dense_right as usize]]
                            .segments
                            .first()
                            .expect("nonempty traversal")
                            .clone(),
                    )
                })
                .collect();
            let partitions = composition_owners[id]
                .iter()
                .map(|&owner| crate::owner_universe_partition(owner, component_locus_to_partition))
                .collect::<Vec<_>>();
            let outcome = span.restricted_charge(
                profile,
                &pairs,
                sources,
                path_of_source,
                &partitions,
                model,
            )?;
            Ok((outcome.charge, outcome.spanning_reads))
        })
        .collect::<io::Result<_>>()?;
    Ok(BoundaryTransitionDraft {
        pair_composition,
        pair_kind,
        composition_pooled,
        composition_owners,
        profiles: composition_profiles,
        novel_scores,
        left_count: sorted_left.len(),
        right_count: sorted_right.len(),
        left_index,
        right_index,
        cooc_pairs,
        cooc_gap_zero,
        cooc_gap_positive,
        novel_pairs,
        overlap_novel,
        empty_tail_pairs,
        novel_riding_pooled,
        novel_only,
    })
}

/// Transition-cost finalize: the POOLED composition losses under the
/// per-feature multiplicity backgrounds, the dense cost matrix, and the
/// boundary's transition statistics.
pub(in super) fn finalize_boundary_transition_costs(
    draft: BoundaryTransitionDraft,
    routed_equal: &crate::RoutedObs,
    component_locus_to_partition: &[u32],
    model: &ScoreModel,
    backgrounds: &FeatureBackgrounds,
    sample: &SampleSideBackgrounds,
) -> io::Result<BoundaryTransitionCosts> {
    let BoundaryTransitionDraft {
        pair_composition,
        pair_kind,
        composition_pooled,
        composition_owners,
        profiles,
        novel_scores,
        left_count,
        right_count,
        left_index,
        right_index,
        cooc_pairs,
        cooc_gap_zero,
        cooc_gap_positive,
        novel_pairs,
        overlap_novel,
        empty_tail_pairs,
        novel_riding_pooled,
        novel_only,
    } = draft;
    let composition_count = profiles.len();
    let mut composition_scores = novel_scores.clone();
    // The haploid track's composition scores: novel (restricted) entries
    // are background-invariant and identical; pooled (co-occurring)
    // compositions charge the sample-side cross-support backgrounds (the
    // haploid single-allele convention).
    let mut composition_scores_haploid = novel_scores;
    for id in 0..composition_count {
        if composition_pooled[id] {
            // Owner-resolved pooled charge: the composition's realizing
            // allele pairs' DISTINCT owning partitions' shares, summed per
            // feature (the boundary convention; the pre-extension model's
            // two locus partitions are the universal single case).
            let mut merged: HashMap<FeatureKey, f64> = HashMap::new();
            for owner in &composition_owners[id] {
                for (feature, share) in crate::owner_routed_obs(
                    routed_equal,
                    *owner,
                    component_locus_to_partition,
                ) {
                    *merged.entry(feature.clone()).or_default() += *share;
                }
            }
            composition_scores[id] = (
                profile_loss_boundary_multiplicity(
                    &profiles[id],
                    &merged,
                    &crate::EMPTY_ROUTED_OBS,
                    model,
                    backgrounds,
                )?,
                0,
            );
            composition_scores_haploid[id] = (
                profile_loss_boundary_sample(
                    &profiles[id],
                    &merged,
                    &crate::EMPTY_ROUTED_OBS,
                    sample,
                    model,
                )?,
                0,
            );
        }
    }
    let novel_spanning_compositions = novel_only
        .iter()
        .filter(|&&id| composition_scores[id].1 > 0)
        .count();
    let novel_spanning_reads: Vec<u64> = novel_only
        .iter()
        .map(|&id| composition_scores[id].1)
        .collect();
    let mut cost = vec![f64::INFINITY; left_count * right_count];
    let mut cost_haploid = vec![f64::INFINITY; left_count * right_count];
    for (index, &composition) in pair_composition.iter().enumerate() {
        cost[index] = composition_scores[composition as usize].0;
        cost_haploid[index] = composition_scores_haploid[composition as usize].0;
    }
    ensure(
        cost.iter().all(|value| value.is_finite())
            && cost_haploid.iter().all(|value| value.is_finite()),
        "phasing boundary cost matrix has unscored entries",
    )?;
    let distribution = |values: &mut Vec<f64>| -> serde_json::Value {
        values.sort_by(|left, right| left.total_cmp(right));
        if values.is_empty() {
            return serde_json::Value::Null;
        }
        serde_json::json!({
            "count": values.len(),
            "min": values[0],
            "median": values[values.len() / 2],
            "max": values[values.len() - 1],
            "negative_count": values.iter().filter(|&&value| value < 0.0).count(),
        })
    };
    let mut cooc_costs: Vec<f64> = Vec::new();
    let mut novel_costs: Vec<f64> = Vec::new();
    for (index, &kind) in pair_kind.iter().enumerate() {
        if kind == 1 {
            cooc_costs.push(cost[index]);
        } else {
            novel_costs.push(cost[index]);
        }
    }
    let stats = serde_json::json!({
        "retained_left": left_count,
        "retained_right": right_count,
        "retained_pairs": left_count * right_count,
        "cooccurring_pairs": cooc_pairs,
        "cooccurring_gap_zero": cooc_gap_zero,
        "cooccurring_gap_positive": cooc_gap_positive,
        "novel_pairs": novel_pairs,
        "novel_from_overlap": overlap_novel,
        "empty_tail_pairs": empty_tail_pairs,
        "distinct_compositions": composition_count,
        "pooled_compositions": composition_count - novel_only.len(),
        "novel_only_compositions": novel_only.len(),
        "novel_compositions_with_spanning_reads": novel_spanning_compositions,
        "novel_pairs_riding_pooled_compositions": novel_riding_pooled,
        "novel_spanning_reads_distribution": {
            "count": novel_spanning_reads.len(),
            "nonzero": novel_spanning_reads.iter().filter(|&&v| v > 0).count(),
            "max": novel_spanning_reads.iter().copied().max().unwrap_or(0),
        },
        "cooccurring_cost_distribution": distribution(&mut cooc_costs),
        "novel_cost_distribution": distribution(&mut novel_costs),
    });
    Ok(BoundaryTransitionCosts {
        cost,
        cost_haploid,
        left_index,
        right_index,
        right_count: right_count as u32,
        stats,
    })
}

// ---------------------------------------------------------------------------
// The exact correlation-phasing chain DP.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
pub(in super) struct PhState {
    score: f64,
    pair: [usize; 2],
    pred: u32,
    /// The state's OWN local loss (the table row's loss; the score is the
    /// best full-path value ending here). The finalist enumerator reconstructs
    /// chains and their totals from (score, loss, transitions).
    loss: f64,
}

/// One locus's ordered pair states, sorted by (dense first column, dense
/// second column) for the transition loop: within a first-column group the
/// second columns ascend, so the per-state cost-matrix loads stream instead
/// of scattering. Sorting changes only iteration order — the min-plus is
/// order-independent, so the DP result is unchanged.
pub(in super) struct LocusStateTable {
    /// Sorted states: (pair, first column, second column, local loss).
    pub(in super) rows: Vec<([usize; 2], u32, u32, f64)>,
    /// `group_start[column]..group_start[column + 1]` is the sorted row
    /// range whose first transition column is `column` (counting-sort
    /// boundaries over the dense right-member columns).
    pub(in super) group_start: Vec<usize>,
}

impl LocusStateTable {
    pub(in super) fn build(mut rows: Vec<([usize; 2], u32, u32, f64)>, columns: usize) -> Self {
        rows.sort_by(|left, right| {
            left.1
                .cmp(&right.1)
                .then_with(|| left.2.cmp(&right.2))
                .then_with(|| left.0.cmp(&right.0))
        });
        let mut group_start = vec![0usize; columns + 1];
        for &(_, first, _, _) in &rows {
            group_start[first as usize + 1] += 1;
        }
        for index in 1..group_start.len() {
            group_start[index] += group_start[index - 1];
        }
        Self { rows, group_start }
    }
}

pub(in super) struct PhasingDpOutcome {
    pub(in super) route: Vec<[usize; 2]>,
    pub(in super) best_score: f64,
    pub(in super) tie_count: usize,
    pub(in super) state_counts: Vec<usize>,
    pub(in super) candidate_counts: Vec<usize>,
    pub(in super) transitions: Vec<u64>,
    pub(in super) suffix_pruned: Vec<u64>,
    pub(in super) wall_seconds: f64,
    /// The full DP layers (score, pair, pred, own local loss per surviving
    /// state) — the k-best finalist enumeration's input.
    pub(in super) layers: Vec<Vec<PhState>>,
}

impl Clone for PhasingDpOutcome {
    fn clone(&self) -> Self {
        Self {
            route: self.route.clone(),
            best_score: self.best_score,
            tie_count: self.tie_count,
            state_counts: self.state_counts.clone(),
            candidate_counts: self.candidate_counts.clone(),
            transitions: self.transitions.clone(),
            suffix_pruned: self.suffix_pruned.clone(),
            wall_seconds: self.wall_seconds,
            layers: self.layers.clone(),
        }
    }
}

/// Plain min-plus over the ordered pair states. Exact: no width, no tie
/// rotations, no seed ladder; the only prune is the admissible suffix
/// bound (observation floors) against the native-backbone incumbent under
/// the phasing transition model — any pruned state provably lies in no
/// chain at or below the incumbent, and the optimum is at or below it.
/// `haploid` selects the ploidy track (owner decision 2026-09-23): every
/// state's second slot is the EMPTY sentinel, so the second boundary
/// charge is zero (a shared zero row replaces the second slot's cost
/// row) and the caller passes the HAPLOID suffix bound (the exact
/// sample-side per-locus and per-boundary minima — the panel pair floor
/// is not admissible for the sample-side charges). The track also selects
/// its own cost matrix: co-occurring transitions charge the sample-side
/// cross-support backgrounds under the haploid track, the panel
/// backgrounds under the diploid track (Step B, owner-approved
/// 2026-09-23).
pub(in super) fn run_phasing_chain_dp(
    locus_count: usize,
    per_locus_tables: &[LocusStateTable],
    boundary_costs: &[BoundaryTransitionCosts],
    suffix_from_locus: &[f64],
    incumbent: f64,
    haploid: bool,
    rss: &mut genome::PeriodicRssGuard,
) -> io::Result<PhasingDpOutcome> {
    let started = Instant::now();
    ensure(
        locus_count > 0 && per_locus_tables.len() == locus_count,
        "dp arity",
    )?;
    let mut layers: Vec<Vec<PhState>> = Vec::with_capacity(locus_count);
    let mut state_counts = Vec::with_capacity(locus_count);
    let mut candidate_counts = Vec::with_capacity(locus_count);
    let mut transitions = Vec::with_capacity(locus_count);
    let mut suffix_pruned = Vec::with_capacity(locus_count);

    // ---- initial layer.
    {
        let table = &per_locus_tables[0];
        let mut layer = Vec::new();
        for &(pair, _, _, loss) in &table.rows {
            if loss + suffix_from_locus[0] <= incumbent + BOUND_PRUNE_EPSILON {
                layer.push(PhState {
                    score: loss,
                    pair,
                    pred: u32::MAX,
                    loss,
                });
            }
        }
        ensure(!layer.is_empty(), "phasing DP initial layer is empty")?;
        state_counts.push(layer.len());
        candidate_counts.push(table.rows.len());
        transitions.push(table.rows.len() as u64);
        suffix_pruned.push((table.rows.len() - layer.len()) as u64);
        layers.push(layer);
    }
    for locus in 1..locus_count {
        rss.checkpoint("phasing_dp_hot_loop")?;
        let previous = layers.last().expect("layer");
        let table = &per_locus_tables[locus];
        let boundary = &boundary_costs[locus - 1];
        let right_count = boundary.right_count as usize;
        let mut best = vec![f64::INFINITY; table.rows.len()];
        let mut pred = vec![u32::MAX; table.rows.len()];
        let states = &table.rows;
        let group_start = &table.group_start;
        // The haploid track's empty second slot carries a zero boundary
        // charge (one shared zero row; the sentinel never indexes the
        // cost matrix).
        let zero_row = if haploid {
            Some(vec![0.0f64; right_count])
        } else {
            None
        };
        for (previous_index, state) in previous.iter().enumerate() {
            let base = state.score;
            let cost_matrix = if haploid {
                &boundary.cost_haploid
            } else {
                &boundary.cost
            };
            let row_first = &cost_matrix
                [boundary.left_index[state.pair[0]] as usize * right_count..][..right_count];
            let row_second: &[f64] = match (&zero_row, haploid) {
                (Some(zero), true) => zero,
                _ => {
                    &cost_matrix[boundary.left_index[state.pair[1]] as usize * right_count..]
                        [..right_count]
                }
            };
            // Per first-column group: hoist the constant prefix, then stream
            // the group's ascending second columns.
            for column in 0..right_count {
                let start = group_start[column];
                if start == group_start[column + 1] {
                    continue;
                }
                let prefix = base + row_first[column];
                let rows = &states[start..group_start[column + 1]];
                let best_slice = &mut best[start..group_start[column + 1]];
                let pred_slice = &mut pred[start..group_start[column + 1]];
                for (offset, &(_, _, second, loss)) in rows.iter().enumerate() {
                    let candidate = prefix + row_second[second as usize] + loss;
                    if candidate < best_slice[offset] {
                        best_slice[offset] = candidate;
                        pred_slice[offset] = previous_index as u32;
                    }
                }
            }
        }
        let mut layer: Vec<PhState> = Vec::new();
        let mut pruned = 0u64;
        for (next, &(pair, _, _, _)) in states.iter().enumerate() {
            if !best[next].is_finite() {
                pruned += 1;
                continue;
            }
            if best[next] + suffix_from_locus[locus] > incumbent + BOUND_PRUNE_EPSILON {
                pruned += 1;
                continue;
            }
            layer.push(PhState {
                score: best[next],
                pair,
                pred: pred[next],
                loss: states[next].3,
            });
        }
        ensure(!layer.is_empty(), "phasing DP layer {locus} is empty")?;
        eprintln!(
            "[phasing] dp layer {locus}: candidates {} states {} suffix_pruned {} ({:.1}s)",
            states.len(),
            layer.len(),
            pruned,
            started.elapsed().as_secs_f64()
        );
        state_counts.push(layer.len());
        candidate_counts.push(states.len());
        transitions.push(previous.len() as u64 * states.len() as u64);
        suffix_pruned.push(pruned);
        layers.push(layer);
    }

    // ---- terminal best, exact tie set, representative backtracking.
    let final_layer = layers.last().expect("final layer");
    let best = final_layer
        .iter()
        .map(|state| state.score)
        .fold(f64::INFINITY, f64::min);
    ensure(best.is_finite(), "phasing DP terminal layer has no finite state")?;
    let tied: Vec<usize> = final_layer
        .iter()
        .enumerate()
        .filter(|&(_, state)| state.score.to_bits() == best.to_bits())
        .map(|(index, _)| index)
        .collect();
    ensure(!tied.is_empty(), "phasing DP tie set is empty")?;
    let mut route = Vec::with_capacity(locus_count);
    let mut cursor = (tied[0], locus_count - 1);
    loop {
        let (index, locus) = cursor;
        let state = &layers[locus][index];
        route.push(state.pair);
        if locus == 0 {
            break;
        }
        cursor = (state.pred as usize, locus - 1);
    }
    route.reverse();
    ensure(route.len() == locus_count, "phasing DP backtrack arity")?;
    Ok(PhasingDpOutcome {
        route,
        best_score: best,
        tie_count: tied.len(),
        state_counts,
        candidate_counts,
        transitions,
        suffix_pruned,
        wall_seconds: started.elapsed().as_secs_f64(),
        layers,
    })
}

// ---------------------------------------------------------------------------
// The k-best finalist enumeration (supervisor ruling, genome/finalist-reranking).
// ---------------------------------------------------------------------------

/// A total-order wrapper for f64 heap keys (`total_cmp`; the DP totals are
/// finite by construction — the enumerator only materializes finite states).
#[derive(Clone, Copy, PartialEq, Debug)]
struct OrdF64(f64);
impl Eq for OrdF64 {}
impl PartialEq<f64> for OrdF64 {
    fn eq(&self, other: &f64) -> bool {
        self.0 == *other
    }
}
impl Ord for OrdF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}
impl PartialOrd for OrdF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// One node's lazy k-best stream: the k-th best full-path value ending at
/// the node, emitted in EXACTLY nondecreasing order. The candidate
/// multiset is { pathval(pred, j) + trans(pred -> node) + loss(node) } over
/// the node's finite-transition predecessors and their already-enumerated
/// path ranks; it is merged from two ordered sources — the predecessors
/// sorted by their OWN best (rank-0) value, and a heap of in-progress
/// (pred, rank >= 1) expansions whose values are materialized (recursively,
/// one locus up) before they enter the heap. Every emitted value records
/// its (pred, pred rank) backpointer, so a popped rank reconstructs its
/// whole chain.
#[derive(Default)]
struct NodeStream {
    /// (pred row at the previous locus, transition cost), sorted ascending
    /// by (pred's rank-0 value + transition) — the node's loss is constant
    /// and added at emission.
    preds: Vec<(u32, f64)>,
    pos: usize,
    /// In-progress expansions (value, pred row, pred rank >= 1), values
    /// precomputed at push time.
    heap: std::collections::BinaryHeap<
        std::cmp::Reverse<(OrdF64, u32, u32)>,
    >,
    /// Emitted (value, pred row, pred rank), in nondecreasing value order.
    out: Vec<(f64, u32, u32)>,
}

/// Lazy k-best enumeration over a haploid DP's layers. Chains are yielded
/// one at a time in EXACTLY nondecreasing total order (the same accumulation
/// order as the DP's forward pass, so a chain's enumerated value is
/// bit-identical to the DP's internal score of the same path). The
/// enumeration covers exactly the DP's surviving layers: states the DP's
/// admissible prune dropped lie on no chain at or below the prune threshold,
/// and the finalist stop threshold is at or below it (verified in-run).
pub(in super) struct FinalistEnumerator<'a> {
    layers: &'a [Vec<PhState>],
    boundary_costs: &'a [BoundaryTransitionCosts],
    /// streams[offset[locus] + row]: the node's lazy k-best stream.
    streams: Vec<NodeStream>,
    offsets: Vec<usize>,
    global: std::collections::BinaryHeap<
        std::cmp::Reverse<(OrdF64, u32, u32)>,
    >,
    exhausted: bool,
    /// Diagnostic: max |stream rank-0 value - the DP's own layer score|
    /// over materialized nodes (expected exactly 0.0 — same accumulation).
    pub(in super) max_seed_delta: f64,
}

impl<'a> FinalistEnumerator<'a> {
    pub(in super) fn new(
        layers: &'a [Vec<PhState>],
        boundary_costs: &'a [BoundaryTransitionCosts],
    ) -> Self {
        let mut offsets = Vec::with_capacity(layers.len() + 1);
        let mut total = 0usize;
        for layer in layers {
            offsets.push(total);
            total += layer.len();
        }
        offsets.push(total);
        Self {
            layers,
            boundary_costs,
            streams: (0..total).map(|_| NodeStream::default()).collect(),
            offsets,
            global: std::collections::BinaryHeap::new(),
            exhausted: false,
            max_seed_delta: 0.0f64,
        }
    }

    fn transition(&self, locus: usize, pred_row: usize, row: usize) -> f64 {
        let costs = &self.boundary_costs[locus - 1];
        let left = self.layers[locus - 1][pred_row].pair[0];
        let right = self.layers[locus][row].pair[0];
        costs.cost_haploid
            [costs.left_index[left] as usize * costs.right_count as usize
                + costs.right_index[right] as usize]
    }

    /// The node's predecessor list, built once on first touch: every
    /// previous-layer row with a finite transition, sorted ascending by
    /// (pred rank-0 value + transition).
    fn ensure_preds(&mut self, locus: usize, row: usize) {
        let id = self.offsets[locus] + row;
        if !self.streams[id].preds.is_empty() || locus == 0 {
            return;
        }
        let mut preds: Vec<(u32, f64)> = Vec::new();
        for pred_row in 0..self.layers[locus - 1].len() {
            let trans = self.transition(locus, pred_row, row);
            if trans.is_finite() {
                preds.push((pred_row as u32, trans));
            }
        }
        let previous_best = |p: u32| self.layers[locus - 1][p as usize].score;
        preds.sort_by(|left, right| {
            let left_key = previous_best(left.0) + left.1;
            let right_key = previous_best(right.0) + right.1;
            left_key.total_cmp(&right_key)
        });
        self.streams[id].preds = preds;
    }

    /// The node's path value for `pred` at `pred_rank` (the summand order
    /// matches the DP's forward pass: prefix + transition + own loss).
    fn candidate_value(&self, locus: usize, row: usize, pred_row: u32, trans: f64, pred_rank: usize) -> f64 {
        let prefix = if pred_rank == 0 {
            self.layers[locus - 1][pred_row as usize].score
        } else {
            self.streams[self.offsets[locus - 1] + pred_row as usize].out[pred_rank].0
        };
        prefix + trans + self.layers[locus][row].loss
    }

    /// Materialize the node's stream up to rank `k` (inclusive).
    fn advance(&mut self, locus: usize, row: usize, k: usize) {
        if locus == 0 {
            // The initial layer: the single rank-0 value IS the state's loss
            // (the DP's own seed). Only rank 0 exists.
            let id = self.offsets[0] + row;
            if self.streams[id].out.is_empty() {
                let state = &self.layers[0][row];
                self.streams[id].out.push((state.loss, u32::MAX, 0));
                self.streams[id].preds = Vec::new();
            }
            return;
        }
        self.ensure_preds(locus, row);
        let id = self.offsets[locus] + row;
        while self.streams[id].out.len() <= k {
            let head_value = {
                let stream = &self.streams[id];
                if stream.pos < stream.preds.len() {
                    let (pred, trans) = stream.preds[stream.pos];
                    Some(self.candidate_value(locus, row, pred, trans, 0))
                } else {
                    None
                }
            };
            let queued_value = self.streams[id]
                .heap
                .peek()
                .map(|entry| ((entry.0).0).0);
            let take_heap = match (head_value, queued_value) {
                (None, None) => break,
                (Some(_), None) => false,
                (None, Some(_)) => true,
                (Some(a), Some(b)) => b <= a,
            };
            if take_heap {
                let std::cmp::Reverse((value, pred, pred_rank)) =
                    self.streams[id].heap.pop().expect("peeked entry");
                self.streams[id].out.push((value.0, pred, pred_rank));
                self.queue_next(locus, row, pred, pred_rank + 1);
            } else {
                let (pred, trans) = {
                    let stream = &self.streams[id];
                    (stream.preds[stream.pos].0, stream.preds[stream.pos].1)
                };
                let value = self.layers[locus - 1][pred as usize].score + trans
                    + self.layers[locus][row].loss;
                let stream = &mut self.streams[id];
                stream.pos += 1;
                stream.out.push((value, pred, 0));
                self.queue_next(locus, row, pred, 1);
            }
        }
    }

    /// Queue the node's next expansion from `pred` at `pred_rank >= 1`:
    /// materialize the pred's rank first (one locus up — the recursion's
    /// only direction), then push the precomputed value. An EXHAUSTED pred
    /// stream (it cannot produce that rank) queues nothing.
    fn queue_next(&mut self, locus: usize, row: usize, pred: u32, pred_rank: u32) {
        self.advance(locus - 1, pred as usize, pred_rank as usize);
        let materialized = self.streams[self.offsets[locus - 1] + pred as usize]
            .out
            .len()
            > pred_rank as usize;
        if !materialized {
            return;
        }
        let value = self.candidate_value(
            locus,
            row,
            pred,
            self.transition(locus, pred as usize, row),
            pred_rank as usize,
        );
        self.streams[self.offsets[locus] + row]
            .heap
            .push(std::cmp::Reverse((OrdF64(value), pred, pred_rank)));
    }

    /// The next chain in exactly nondecreasing total order, with its total.
    pub(in super) fn next_chain(&mut self) -> Option<(f64, Vec<[usize; 2]>)> {
        if self.exhausted {
            return None;
        }
        let last = self.layers.len() - 1;
        if self.global.is_empty() {
            for (row, state) in self.layers[last].iter().enumerate() {
                self.global
                    .push(std::cmp::Reverse((OrdF64(state.score), row as u32, 0u32)));
            }
        }
        let std::cmp::Reverse((value, row, k)) = self.global.pop()?;
        // Materialize the chain's backpointers (the seed value is the DP's
        // own layer score; deeper ranks come from the node's stream).
        if k > 0 {
            self.advance(last, row as usize, k as usize);
        } else {
            self.advance(last, row as usize, 0);
        }
        let seed_delta = if k == 0 {
            (self.streams[self.offsets[last] + row as usize].out[0].0 - self.layers[last][row as usize].score)
                .abs()
        } else {
            0.0
        };
        if seed_delta > self.max_seed_delta {
            self.max_seed_delta = seed_delta;
        }
        // Reconstruct the chain.
        let mut route = Vec::with_capacity(self.layers.len());
        let mut cursor = (last, row as usize, k as usize);
        loop {
            let (locus, chain_row, chain_k) = cursor;
            let state = &self.layers[locus][chain_row];
            route.push(state.pair);
            if locus == 0 {
                break;
            }
            let (_, pred, pred_k) = self.streams[self.offsets[locus] + chain_row].out[chain_k];
            cursor = (locus - 1, pred as usize, pred_k as usize);
        }
        route.reverse();
        // Queue the final node's next rank.
        self.advance(last, row as usize, k as usize + 1);
        let next_k = k as usize + 1;
        if self.streams[self.offsets[last] + row as usize].out.len() > next_k {
            let value = self.streams[self.offsets[last] + row as usize].out[next_k].0;
            self.global
                .push(std::cmp::Reverse((OrdF64(value), row, next_k as u32)));
        }
        if self.global.is_empty() {
            self.exhausted = true;
        }
        Some((value.0, route))
    }
}

/// The weighted intersection of two sorted index lists over ONE window's
/// observed-feature space (the terms are that window's per-feature omission
/// values, identical for both lists).
fn coupling_overlap(list_a: &[u32], list_b: &[u32], terms: &[f64]) -> f64 {
    let mut total = 0.0f64;
    let (mut i, mut j) = (0usize, 0usize);
    while i < list_a.len() && j < list_b.len() {
        match list_a[i].cmp(&list_b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                total += terms[list_a[i] as usize];
                i += 1;
                j += 1;
            }
        }
    }
    total
}

/// The pairwise-coupling SEARCH objective's boundary matrices (supervisor
/// ruling on genome/finalist-reranking, option (b)): C_b(a, a') =
/// trans_b(a, a') - exon_pair(b; a, a'), where exon_pair prices the
/// ADJACENT-window exoneration channel with the SAME per-feature omission
/// terms the model of record charges: exon(b+1<-b) = the window-(b+1)
/// observed mass the window-b allele spells and the window-(b+1) allele does
/// not predict; exon(b<-b+1) = the mirror. C is an upper-bound form on the
/// model of record (the model's exoneration channel is the union over ALL
/// window pairs; the adjacent channel is its pairwise slice), so its
/// enumeration order is model-aware and its frontier supports the measured-
/// slack stop. ROLE SEPARATION (absolute, per the ruling): C exists ONLY for
/// search ordering — the selection and every reported score remain exact
/// model-of-record numbers from the re-ranker.
pub(in super) fn build_coupling_transition_costs(
    // The finalist states' allele ids, in the tables' dense order (the same
    // order the finalist transition matrices' index maps encode).
    sorted_left: &[usize],
    sorted_right: &[usize],
    // The tracks' finalist transition matrices (the trans component).
    costs: &BoundaryTransitionCosts,
    // The two loci's per-allele exonation-channel inputs (indexed by ALLELE).
    left_exons: &[HaploidAlleleExon],
    right_exons: &[HaploidAlleleExon],
    // The two windows' per-feature omission-term vectors.
    terms_left: &[f64],
    terms_right: &[f64],
    // The INSTANCE-LEVEL coupling (the corrected Ruling 2's pairwise
    // channel): when true, the exon terms price the UNCOVERED instance
    // masses (the neighbor's spelling credits only the instances it
    // actually covers at the matching positions); the per-window feature
    // keys, observed masses and per-feature instance shares below are the
    // term arithmetic's inputs. When false, the established feature-level
    // form (the exon lists' binary overlap).
    instance_coupling: bool,
    sample: Option<&SampleSideBackgrounds>,
    keys_left: &[FeatureKey],
    keys_right: &[FeatureKey],
    full_left: &[f64],
    full_right: &[f64],
    shares_left: &[Vec<f64>],
    shares_right: &[Vec<f64>],
) -> BoundaryTransitionCosts {
    let right_count = costs.right_count as usize;
    let pairs = sorted_left.len() * right_count;
    let mut cost = vec![0.0f64; pairs];
    // The instance-level marginal exonation credit: the neighbor covers
    // delta mass of the feature's instances that the self row does not
    // already cover; the credit = term(u_self) - term(u_self - delta) over
    // the UNCOVERED masses (the Poisson count's own argument — the same
    // nonlinearity the corrected chain-level addend prices with).
    // The instance-level marginal exonation credit (the sorted merge-join
    // form — the hash-loop form measured 1710s vs the legacy form's 59s):
    // the neighbor covers delta mass of the feature's instances that the
    // self row does not already cover; the credit =
    // term(u_self) - term(u_self - delta) over the UNCOVERED masses (the
    // Poisson count's own argument — the same nonlinearity the corrected
    // chain-level addend prices with). Features the neighbor covers
    // NOTHING of precompute to the allele's own solo credit (no pair work).
    let instance_exon = |covered: &AlleleCoveredInstances,
                         other: &AlleleCoveredInstances,
                         keys: &[FeatureKey],
                         full: &[f64],
                         shares: &[Vec<f64>]| {
        if covered.is_empty() {
            return 0.0f64;
        }
        let sample = sample.expect("the instance coupling implies the sample-side backgrounds");
        let mut exon = 0.0f64;
        let mut i = 0usize; // covered.features
        let mut j = 0usize; // other.features
        while i < covered.features.len() {
            let f = covered.features[i] as usize;
            if j >= other.features.len() || (other.features[j] as usize) > f {
                // The neighbor covers nothing of this feature: the solo
                // credit (precomputed per allele).
                exon += covered.solo_credits[i];
                i += 1;
            } else if (other.features[j] as usize) < f {
                j += 1;
            } else {
                // BOTH cover the feature: re-price the remainder.
                let ids_a = &covered.ids[i];
                let ids_ap = &other.ids[j];
                let u_ap = (full[f] - other.masses[j]).max(0.0);
                // mass(Sa \ Sa'): the merge join over the sorted id lists.
                let mut delta = 0.0f64;
                let mut x = 0usize;
                let mut y = 0usize;
                let shares_f = &shares[f];
                while x < ids_a.len() {
                    if y >= ids_ap.len() || ids_a[x] < ids_ap[y] {
                        delta += shares_f[ids_a[x] as usize];
                        x += 1;
                    } else if ids_a[x] > ids_ap[y] {
                        y += 1;
                    } else {
                        x += 1;
                        y += 1;
                    }
                }
                let u_joint = (u_ap - delta).max(0.0);
                exon += crate::omission_term_value(&keys[f], u_ap, sample)
                    - crate::omission_term_value(&keys[f], u_joint, sample);
                i += 1;
                j += 1;
            }
        }
        exon
    };
    for (dense_l, &a) in sorted_left.iter().enumerate() {
        let left_exon = &left_exons[a];
        for (dense_r, &b) in sorted_right.iter().enumerate() {
            let right_exon = &right_exons[b];
            let trans = costs.cost_haploid[dense_l * right_count + dense_r];
            // exon(b+1<-b): the window-(b+1) observed mass the LEFT allele
            // spells, minus what BOTH spell there.
            let exon_forward = if instance_coupling {
                match &left_exon.right_instances {
                    Some(covered) => instance_exon(
                        covered,
                        &right_exon.self_instances,
                        keys_right,
                        full_right,
                        shares_right,
                    ),
                    None => 0.0,
                }
            } else {
                match (&left_exon.right, &right_exon.self_list) {
                    (Some((list, sum)), self_list) => {
                        sum - coupling_overlap(list, self_list, terms_right)
                    }
                    _ => 0.0,
                }
            };
            // exon(b<-b+1): the mirror through the window-b terms.
            let exon_backward = if instance_coupling {
                match &right_exon.left_instances {
                    Some(covered) => instance_exon(
                        covered,
                        &left_exon.self_instances,
                        keys_left,
                        full_left,
                        shares_left,
                    ),
                    None => 0.0,
                }
            } else {
                match (&right_exon.left, &left_exon.self_list) {
                    (Some((list, sum)), self_list) => {
                        sum - coupling_overlap(list, self_list, terms_left)
                    }
                    _ => 0.0,
                }
            };
            let value = trans - exon_forward - exon_backward;
            cost[dense_l * right_count + dense_r] = value;
        }
    }
    BoundaryTransitionCosts {
        cost: cost.clone(),
        cost_haploid: cost,
        left_index: costs.left_index.clone(),
        right_index: costs.right_index.clone(),
        right_count: costs.right_count,
        stats: serde_json::Value::Null,
    }
}

// ---------------------------------------------------------------------------
// The correlation-phasing spine tail (Stages P1-P4 after the shared Stage 1).
// ---------------------------------------------------------------------------

/// Shared state handed over from `run_local_first_spine` after its Stage 1
/// (the exhaustive local sweep) completes.
pub(in super) struct PhasingShared<'a> {
    pub(in super) panel: &'a SyngIndex,
    pub(in super) sources: &'a routes::Sources,
    pub(in super) graph: &'a routes::Graph,
    pub(in super) ports: &'a mut routes::Ports,
    pub(in super) flank_memo: &'a FlankMemo,
    pub(in super) path_of_source: &'a [usize],
    pub(in super) axis_slice: &'a [genome::AxisInterval],
    pub(in super) territory: &'a [Vec<SourceRange>],
    pub(in super) partition_obs: &'a [HashMap<FeatureKey, f64>],
    /// The routed shares + the component locus -> universe partition map
    /// (the window-domain extension's owner-resolved charging).
    pub(in super) routed_equal: &'a crate::RoutedObs,
    pub(in super) component_locus_to_partition: &'a [u32],
    /// Per slice locus: the universe partition (the crossing reads' share
    /// target at every junction the phasing scores).
    pub(in super) locus_partitions: Vec<u32>,
    pub(in super) k: u64,
    pub(in super) locus_offset: usize,
    pub(in super) target_source: usize,
    pub(in super) model: &'a ScoreModel,
    pub(in super) sample_counts: &'a impg::sample_mem_bwt::WeightedBwt,
    pub(in super) reference_specs: &'a [(String, PathBuf, PathBuf)],
    pub(in super) truth_route_a: &'a routes::Route,
    pub(in super) truth_route_b: &'a routes::Route,
    pub(in super) rescore_directory: &'a std::path::Path,
    pub(in super) rss: &'a mut genome::PeriodicRssGuard,
    /// The within-read adjacency index (novel transitions are paid by
    /// crossing reads only).
    pub(in super) span: &'a JunctionSpanIndex,
    // Working state (owned by the phasing tail from here on).
    pub(in super) ranges: Vec<Vec<genome::SpanningTraversal>>,
    pub(in super) locus_classes: Vec<LocusClassing>,
    pub(in super) folded: Vec<LocusFolded>,
    pub(in super) loss_tables: Vec<Vec<Vec<f64>>>,
    pub(in super) boundaries: Vec<SpineBoundary>,
    pub(in super) successors_ref: Vec<Vec<Vec<(usize, f64)>>>,
    pub(in super) viable: Vec<Vec<bool>>,
    pub(in super) sweeps: Vec<LocusSweep>,
    pub(in super) margins: Vec<f64>,
    pub(in super) incumbent_spine: f64,
    pub(in super) backbone_chain: Vec<usize>,
    pub(in super) truth_pieces: Vec<[Vec<(usize, u64, u64, bool, u32)>; 2]>,
    pub(in super) truth_losses: Vec<f64>,
    pub(in super) donor_sources: BTreeSet<usize>,
    pub(in super) tract_loci: Vec<usize>,
    pub(in super) stage1_summary: serde_json::Value,
    /// The per-feature multiplicity backgrounds (owned by the spine; the
    /// phasing extends the scan for features only it charges).
    pub(in super) backgrounds: &'a mut FeatureBackgrounds,
    /// The sample-side cross-support backgrounds for the HAPLOID track
    /// (owner-approved family, genome-wide form 2026-09-24): the haploid
    /// single-allele table charges beta_f = base +
    /// Sum_r m_r*(1 - 1/t_r_genome) — the sample's own other-copy support
    /// over the genome-universe placement pass — instead of the panel
    /// multiplicity background (the panel's other contexts generate zero
    /// reads in a haploid sample; see genome/haploid-local-table.md and
    /// genome/genome-wide-cross-support.md for the measurements). Always
    /// present for phasing runs (the routing pass builds the genome-wide
    /// cross support whenever --spine-phasing is set).
    pub(in super) sample_backgrounds: Option<&'a SampleSideBackgrounds<'a>>,
    /// The record-once site-map builder (Fix 1's mixed-owner charging) and
    /// the windows' full observed profiles (Fix 2's omission universes).
    pub(in super) site: &'a crate::SiteObserved,
    pub(in super) window_obs: &'a [HashMap<FeatureKey, f64>],
    /// The windows' observed INSTANCE structure (the S2 record-level
    /// coverage's observed side). `None` keeps the feature-level
    /// chain-level form bit-identically.
    pub(in super) window_instances: Option<&'a crate::InstanceStructure>,
    /// The instance-probe diagnostic's feature list (the pooled
    /// decomposition's top truth-advantage features; empty = no probe
    /// report). Operational diagnostic input, not a model constant.
    pub(in super) instance_probe_features: Vec<FeatureKey>,
    /// Probe routes (repeatable --probe-route LABEL:PATH.json): previously
    /// measured chains mapped into the finalist tables — the re-rank's
    /// bit-identity anchors (the predecessor's selected chain and
    /// donor-probe chain) and the surrogate-rank measurement the ruling
    /// asks for. Their model totals come from the same evaluator as every
    /// finalist's.
    pub(in super) probe_specs: &'a [(String, std::path::PathBuf)],
}

/// The whole correlation-phasing flow: co-occurrence extraction, structural
/// detection + Option A partials, the phasing DP, the authoritative rescore,
/// and the tract verdict. Returns the complete `--spine` JSON value.
/// The haploid chain's CHAIN-LEVEL model-of-record total (Ruling 2's
/// evaluator: `score_reference_m1`'s haploid branch with the once-per-material
/// piece filter and the chain-level omission addend). This is THE finalist
/// re-ranking's exact per-chain evaluation AND the ladder's
/// `selected_chain_level_m1` — one call shape, so a finalist's re-rank total
/// is bit-identical to what the ladder reports when that finalist is
/// selected. `shared_piece_memo` shares the territory pieces' oracle profiles
/// across repeated evaluations (values are bit-identical either way).
#[allow(clippy::too_many_arguments)]
fn haploid_chain_level_m1(
    panel: &SyngIndex,
    sources: &routes::Sources,
    territory: &[Vec<SourceRange>],
    routed_equal: &crate::RoutedObs,
    site: &crate::SiteObserved,
    window_obs: &[HashMap<FeatureKey, f64>],
    // The windows' observed INSTANCE structure (the S2 record-level
    // coverage's observed side); `None` keeps the feature-level form
    // bit-identically.
    window_instances: Option<&crate::InstanceStructure>,
    component_locus_to_partition: &[u32],
    locus_count: usize,
    locus_offset: usize,
    model: &ScoreModel,
    span: &JunctionSpanIndex,
    path_of_source: &[usize],
    locus_partitions: &[u32],
    backgrounds: &mut FeatureBackgrounds,
    sample_backgrounds: &SampleSideBackgrounds,
    // The haploid chain: the slot-0 traversal at each locus.
    chain: &[&genome::SpanningTraversal],
    shared_piece_memo: Option<&mut HashMap<(usize, u64, u64), Profile>>,
    shared_spelled_memo: Option<&mut crate::SpelledMemo>,
    // OUT (the instance-probe diagnostic): the chain's spelled anchor
    // positions per feature.
    mut spelled_out: Option<&mut HashMap<FeatureKey, Vec<crate::SpelledSpan>>>,
) -> io::Result<f64> {
    let mut segments: Vec<routes::Segment> = Vec::new();
    for traversal in chain {
        for segment in &traversal.segments {
            if segment.start >= segment.end {
                continue;
            }
            segments.push(routes::Segment {
                source: segment.source,
                start: segment.start,
                end: segment.end,
                reverse: segment.reverse,
            });
        }
    }
    let route_a = routes::Route { segments };
    let route_b = routes::Route {
        segments: Vec::new(),
    };
    let (loss, _) = score_reference_m1(
        panel,
        sources,
        territory,
        [&route_a, &route_b],
        routed_equal,
        site,
        window_obs,
        window_instances,
        component_locus_to_partition,
        locus_count,
        locus_offset,
        model,
        Some((span, path_of_source, locus_partitions)),
        Some(&mut *backgrounds),
        Some(sample_backgrounds),
        true,
        shared_piece_memo,
        shared_spelled_memo,
        spelled_out,
    )?;
    Ok(loss)
}

pub(in super) fn run_correlation_phasing(
    shared: &mut PhasingShared,
) -> io::Result<serde_json::Value> {
    let started = Instant::now();
    let PhasingShared {
        panel,
        sources,
        graph,
        ports,
        flank_memo,
        path_of_source,
        axis_slice,
        territory,
        partition_obs,
        routed_equal,
        component_locus_to_partition,
        window_instances,
        locus_partitions,
        k,
        locus_offset,
        target_source,
        model,
        sample_counts,
        reference_specs,
        truth_route_a,
        truth_route_b,
        rescore_directory,
        rss,
        span,
        ranges,
        locus_classes,
        folded,
        loss_tables,
        boundaries,
        successors_ref,
        viable,
        sweeps,
        margins,
        incumbent_spine,
        backbone_chain,
        truth_pieces,
        truth_losses,
        donor_sources,
        tract_loci,
        stage1_summary,
        backgrounds,
        sample_backgrounds,
        site,
        window_obs,
        probe_specs,
        instance_probe_features,
    } = shared;
    // Option<&[...]> is Copy; the &mut destructure binds a reborrowable
    // reference, so copy it out once for the plain-Option call sites.
    let window_instances: Option<&crate::InstanceStructure> = *window_instances;
    let locus_count = ranges.len();
    ensure(locus_count > 1, "phasing needs at least two loci")?;
    let donor_source = donor_sources.iter().copied().next();
    // The haploid track's sample-side cross-support backgrounds (Step B):
    // phasing runs always build the pooled support, so this is present.
    let sample_backgrounds: &SampleSideBackgrounds = sample_backgrounds
        .as_ref()
        .ok_or_else(|| {
            invalid(
                "phasing requires the pooled sample support (the sample-side \
                 haploid background); the routing pass builds it under \
                 --spine-phasing"
            )
        })?;

    // ---------------------------------------------- P1: co-occurrence + retention.
    let cooc_started = Instant::now();
    let pre_retained: Vec<HashSet<usize>> = (0..locus_count)
        .map(|locus| {
            retained_members(
                locus,
                &locus_classes[locus],
                &sweeps[locus],
                &locus_classes[locus].membership,
                backbone_chain[locus],
                margins[locus],
            )
        })
        .collect();
    let cooccurrence_rows: Vec<serde_json::Value> = (0..locus_count - 1)
        .into_par_iter()
        .map(|boundary| {
            cooccurrence_diagnostic(
                boundary,
                &ranges[boundary],
                &ranges[boundary + 1],
                &successors_ref[boundary],
                &pre_retained[boundary],
                &pre_retained[boundary + 1],
                donor_source,
            )
        })
        .collect();
    let cooc_seconds = cooc_started.elapsed().as_secs_f64();
    eprintln!(
        "[phasing] co-occurrence diagnostic done over {} boundaries ({cooc_seconds:.2}s)",
        locus_count - 1
    );

    // ---------------------------------------- P2 (Option A): structural partials.
    // The firing rule is measured before it is trusted: the naive
    // all-pairs rule fired 4.8M partial rows on this slice (every foreign
    // source's misaligned grid overlaps somewhere among the margin-retained
    // chimera rows), so the rule is the CANONICAL per-source minimal-overlap
    // failure pair — bounded by the panel's source count per boundary, and
    // reported per boundary. The existing split machinery is NOT fired from
    // boundary failures: under the corrected continuity structure the
    // junction loci (12/26) have no port-word continuity failure to fire on
    // (the donor's grid butts exactly there), so a boundary-failure trigger
    // reaches only chimera territory — reported as a finding, with the
    // post-DP novel-junction detector spec'd as the natural public trigger
    // for mid-locus junction partials.
    let option_a_started = Instant::now();
    let mut partial_additions: Vec<Vec<genome::SpanningTraversal>> =
        vec![Vec::new(); locus_count];
    let mut option_a_stats: Vec<serde_json::Value> = Vec::new();
    for boundary in 0..locus_count - 1 {
        let canonical = canonical_overlap_failures(
            boundary,
            &ranges[boundary],
            &ranges[boundary + 1],
            &successors_ref[boundary],
            &pre_retained[boundary],
            &pre_retained[boundary + 1],
        );
        let (left_partials, right_partials, stats) = overlap_failure_partials(
            boundary,
            &canonical,
            &ranges[boundary],
            &ranges[boundary + 1],
            graph,
            ports,
        )?;
        let canonical_rows: Vec<serde_json::Value> = canonical
            .iter()
            .map(|&(source, left, right, overlap)| {
                serde_json::json!({
                    "source": source,
                    "left": left,
                    "right": right,
                    "overlap": overlap,
                })
            })
            .collect();
        option_a_stats.push(serde_json::json!({
            "boundary": boundary,
            "canonical_failures": canonical.len(),
            "left_partials": left_partials.len(),
            "right_partials": right_partials.len(),
            "canonical_failure_rows": canonical_rows,
            "detail": stats,
        }));
        for partial in left_partials {
            if !ranges[boundary]
                .iter()
                .any(|row| row.segments == partial.segments)
            {
                partial_additions[boundary].push(partial);
            }
        }
        for partial in right_partials {
            if !ranges[boundary + 1]
                .iter()
                .any(|row| row.segments == partial.segments)
            {
                partial_additions[boundary + 1].push(partial);
            }
        }
    }
    let option_a_seconds = option_a_started.elapsed().as_secs_f64();
    let augmented: Vec<usize> = (0..locus_count)
        .filter(|&locus| !partial_additions[locus].is_empty())
        .collect();
    let total_partials: usize = partial_additions.iter().map(|list| list.len()).sum();
    eprintln!(
        "[phasing] option A: {} partial rows at loci {:?} ({option_a_seconds:.2}s)",
        total_partials, augmented
    );

    // ---------------------------------------- P2b: rebuild at the augmented loci.
    let rebuild_started = Instant::now();
    let mut affected_boundaries: BTreeSet<usize> = BTreeSet::new();
    for &locus in &augmented {
        for partial in partial_additions[locus].iter() {
            if !ranges[locus]
                .iter()
                .any(|row| row.segments == partial.segments)
            {
                ranges[locus].push(partial.clone());
            }
        }
        if locus > 0 {
            affected_boundaries.insert(locus - 1);
        }
        if locus + 1 < locus_count {
            affected_boundaries.insert(locus);
        }
    }
    for &locus in &augmented {
        locus_classes[locus] = class_locus_alleles(
            panel,
            sources,
            flank_memo,
            path_of_source,
            &ranges[locus],
            locus,
            *k,
            span,
            locus_partitions[locus],
            routed_equal,
            component_locus_to_partition,
            model,
        )?;
    }
    // Boundary drafts at the affected boundaries.
    let mut rebuild_drafts: BTreeMap<usize, SpineBoundaryDraft> = BTreeMap::new();
    for &boundary in &affected_boundaries {
        let axis_window = &axis_slice[boundary..boundary + 2];
        let range_window = &ranges[boundary..boundary + 2];
        let links = genome::port_word_seams(axis_window, range_window, graph, ports)?;
        let links = links
            .into_iter()
            .next()
            .ok_or_else(|| invalid("phasing boundary link cardinality mismatch"))?;
        ensure(
            !links.is_empty(),
            "phasing boundary {boundary} lost all legal links after partials",
        )?;
        let left_owners: Vec<u32> =
            ranges[boundary].iter().map(crate::traversal_owner).collect();
        let right_owners: Vec<u32> =
            ranges[boundary + 1].iter().map(crate::traversal_owner).collect();
        rebuild_drafts.insert(
            boundary,
            build_spine_boundary_draft(
                panel,
                sources,
                flank_memo,
                model,
                span,
                path_of_source,
                &left_owners,
                &right_owners,
                routed_equal,
                component_locus_to_partition,
                &ranges[boundary],
                &ranges[boundary + 1],
                &links,
            )?,
        );
    }
    // Extend the multiplicity backgrounds for every feature the augmented
    // loci and rebuilt boundary compositions charge that the scan has not
    // measured yet.
    {
        let mut universe: BTreeSet<FeatureKey> = BTreeSet::new();
        for &locus in &augmented {
            for profile in &locus_classes[locus].profiles {
                for key in profile.keys() {
                    universe.insert(key.clone());
                }
            }
        }
        for draft in rebuild_drafts.values() {
            for key in draft_feature_set(draft) {
                universe.insert(key);
            }
        }
        let features: Vec<FeatureKey> = universe.into_iter().collect();
        backgrounds.scan_extend(panel, features, *k, model)?;
    }
    for &locus in &augmented {
        folded[locus] = LocusFolded::build_backgrounds(
            &locus_classes[locus].profiles,
            &locus_classes[locus].class_owners,
            routed_equal,
            component_locus_to_partition,
            backgrounds,
        )?;
        loss_tables[locus] = folded[locus]
            .loss_tables(model)
            .map(|(tables, _)| tables)?;
    }
    for (boundary, draft) in rebuild_drafts {
        let left_owners: Vec<u32> =
            ranges[boundary].iter().map(crate::traversal_owner).collect();
        let right_owners: Vec<u32> =
            ranges[boundary + 1].iter().map(crate::traversal_owner).collect();
        boundaries[boundary] = finalize_spine_boundary(
            draft,
            &ranges[boundary],
            &ranges[boundary + 1],
            &left_owners,
            &right_owners,
            routed_equal,
            component_locus_to_partition,
            model,
            backgrounds,
        )?;
        successors_ref[boundary] = boundaries[boundary].successors.clone();
    }
    *viable = spine_viability(ranges, successors_ref);
    for &locus in &augmented {
        sweeps[locus] = exhaustive_local_sweep(
            &folded[locus],
            &loss_tables[locus],
            &locus_classes[locus].membership,
            &viable[locus],
            Some(backbone_chain[locus]),
            model,
            &locus_classes[locus].class_charges,
            &scorable_classes(&locus_classes[locus], &ranges[locus]),
        )?;
    }
    *margins = local_seam_swing_margins(sweeps, successors_ref, viable);
    let rebuild_seconds = rebuild_started.elapsed().as_secs_f64();
    eprintln!(
        "[phasing] rebuild of {} loci done ({rebuild_seconds:.2}s)",
        augmented.len()
    );
    rss_probe(rss, "phasing_after_rebuild")?;

    // Same-owner stitched spanning candidates (supervisor ruling,
    // genome/stitching-omission-alignment 2026-09-24): chains of ADJACENT
    // same-source forward rows whose union contains the window's axis
    // interval. HAPLOID-ONLY augmentation: the stitched alleles join the
    // haploid track's candidate list, losses, retention, states, transition
    // matrices, route and rescores; the DIPLOID track's inputs (classing,
    // sweeps, folded tables, margins) stay bit-untouched by construction
    // (the frozen-track rule). The combined allele list = the domain rows
    // (indices 0..ranges.len(), the classing/sweeps universe) followed by
    // the stitched candidates; every haploid index below addresses it.
    let stitched_started = Instant::now();
    let stitched: Vec<Vec<genome::SpanningTraversal>> = (0..locus_count)
        .map(|locus| {
            genome::stitched_candidates(
                *locus_offset + locus,
                &ranges[locus],
                axis_slice[locus].start,
                axis_slice[locus].end,
            )
        })
        .collect();
    let stitched_count: usize = stitched.iter().map(Vec::len).sum();
    let stitched_seconds = stitched_started.elapsed().as_secs_f64();
    eprintln!(
        "[phasing] same-owner stitched spanning candidates: {} over {} loci ({stitched_seconds:.2}s)",
        stitched_count,
        stitched.len()
    );
    let combined: Vec<Vec<genome::SpanningTraversal>> = (0..locus_count)
        .map(|locus| {
            ranges[locus]
                .iter()
                .cloned()
                .chain(stitched[locus].iter().cloned())
                .collect()
        })
        .collect();

    // Shared reborrowings for the parallel stages (no further mutation of
    // the working state happens after the rebuild).
    let panel: &SyngIndex = panel;
    let sources: &routes::Sources = sources;
    let flank_memo: &FlankMemo = flank_memo;
    let model: &ScoreModel = model;
    let partition_obs: &[HashMap<FeatureKey, f64>] = partition_obs;
    let ranges: &[Vec<genome::SpanningTraversal>] = ranges;
    let locus_classes: &[LocusClassing] = locus_classes;
    let sweeps: &[LocusSweep] = sweeps;
    let margins: &[f64] = margins;
    let backbone_chain: &[usize] = backbone_chain;
    let donor_sources: &BTreeSet<usize> = donor_sources;
    let locus_offset: usize = *locus_offset;
    let target_source: usize = *target_source;
    let incumbent_spine: f64 = *incumbent_spine;

    // ---------------------------------------- P3: transition costs + floors.
    // The transition matrices are computed first; the boundary floors for
    // the DP's admissible suffix bound are then the EXACT minima of the
    // computed charge matrices (the achievable transition set), which is
    // admissible and strictly tighter than the old per-feature bound —
    // necessary here because a restricted (crossing-read-paid) charge's
    // per-feature minimum over unknown observed counts is not bounded
    // below, while the matrix minimum is exact.
    let costs_started = Instant::now();
    let membership_slices: Vec<Vec<usize>> = (0..locus_count)
        .map(|locus| locus_classes[locus].membership.clone())
        .collect();
    let retained_member_sets: Vec<HashSet<usize>> = (0..locus_count)
        .map(|locus| {
            retained_members(
                locus,
                &locus_classes[locus],
                &sweeps[locus],
                &membership_slices[locus],
                backbone_chain[locus],
                margins[locus],
            )
        })
        .collect();
    // The haploid track's omission-inclusive single-allele losses (Fix 2),
    // computed BEFORE the retention so the haploid candidate universe can
    // reflect the new objective: a donor-spelling class the DIPLOID pair
    // margins exclude is retained for the haploid track when its
    // omission-inclusive single loss is within the locus's own swing margin
    // of the best single loss (the single table's mirror of the pair
    // rule's threshold form — same derived margin, no new constants). The
    // DIPLOID track's retained sets stay exactly as they were
    // (bit-identical candidate universe); the boundary drafts enumerate the
    // UNION so the widened haploid states have transition costs, while the
    // diploid DP reads only its own states' matrix entries (values
    // unchanged).
    let haploid_oracle_started = Instant::now();
    // The windows' observed-feature indices (the exonation channel's term
    // tables; the same per-feature values the omission charge uses).
    let window_obs_index: Vec<crate::WindowObsIndex> = (0..locus_count)
        .map(|locus| crate::build_window_obs_index(&window_obs[locus], sample_backgrounds))
        .collect();
    // The instance-level coupling's per-window inputs (aligned with the obs
    // index): each feature's observed mass and its instance entries' shares
    // (entry order = the covered ids' order).
    let window_instance_tables: Vec<(Vec<f64>, Vec<Vec<f64>>)> = (0..locus_count)
        .map(|locus| {
            let idx = &window_obs_index[locus];
            let full: Vec<f64> = idx
                .features
                .iter()
                .map(|f| window_obs[locus][f])
                .collect();
            let shares: Vec<Vec<f64>> = idx
                .features
                .iter()
                .map(|f| {
                    window_instances
                        .map(|w| {
                            w.window_records
                                .get(locus)
                                .and_then(|m| m.get(f))
                                .map(|records| {
                                    records
                                        .iter()
                                        .map(|&r| {
                                            w.record_shares
                                                .get(r as usize)
                                                .copied()
                                                .unwrap_or(0.0)
                                        })
                                        .collect()
                                })
                                .unwrap_or_default()
                        })
                        .unwrap_or_default()
                })
                .collect();
            (full, shares)
        })
        .collect();
    let haploid_oracle: Vec<(Vec<f64>, Vec<HaploidAlleleExon>)> = (0..locus_count)
        .into_par_iter()
        .map(|locus| {
            let scorable = scorable_classes(&locus_classes[locus], &ranges[locus]);
            haploid_allele_losses(
                &combined[locus],
                ranges[locus].len(),
                &locus_classes[locus],
                &window_obs[locus],
                &window_obs_index,
                window_obs,
                locus,
                routed_equal,
                site,
                component_locus_to_partition,
                sample_backgrounds,
                &scorable,
                model,
                panel,
                sources,
                window_instances,
                path_of_source,
            )
        })
        .collect::<io::Result<_>>()?;
    let haploid_oracle_seconds = haploid_oracle_started.elapsed().as_secs_f64();
    eprintln!(
        "[phasing] haploid oracle allele losses: {:.2}s",
        haploid_oracle_seconds
    );
    let haploid_losses: Vec<Vec<f64>> = haploid_oracle
        .iter()
        .map(|(losses, _)| losses.clone())
        .collect();
    let haploid_exons: Vec<Vec<HaploidAlleleExon>> = haploid_oracle
        .into_iter()
        .map(|(_, exons)| exons)
        .collect();
    let retained_haploid_sets: Vec<HashSet<usize>> = (0..locus_count)
        .map(|locus| {
            let losses = &haploid_losses[locus];
            let best_single = losses
                .iter()
                .copied()
                .filter(|value| value.is_finite())
                .fold(f64::INFINITY, f64::min);
            // The single table's mirror of the pair rule's threshold form:
            // the pair rule retains pairs within the locus's SEAM-SWING
            // margin of its best viable pair (margin = min_viable_pair +
            // swing); the single rule retains singles within the SWING of
            // its own best single. Anchoring the singles' threshold at the
            // pair loss's absolute value (best_single + full margin)
            // mis-scales across the two objectives (measured on try 5: the
            // truth's own classes fell ~6k outside the threshold at the
            // tract loci while winning those loci by 12k-124k under the new
            // objective). swing = margin - min_viable_pair, both existing
            // derived quantities; no new constants.
            let swing = margins[locus] - sweeps[locus].min_viable_loss;
            let threshold = best_single + swing;
            let mut set = retained_member_sets[locus].clone();
            for (allele, &loss) in losses.iter().enumerate() {
                if loss.is_finite() && loss <= threshold {
                    set.insert(allele);
                }
            }
            set
        })
        .collect();
    let retained_union_sets: Vec<HashSet<usize>> = (0..locus_count)
        .map(|locus| {
            let mut set = retained_member_sets[locus].clone();
            for allele in &retained_haploid_sets[locus] {
                set.insert(*allele);
            }
            set
        })
        .collect();
    // ------------------------------------------------------------- finalists.
    // The k-best finalist machinery's universe (supervisor ruling,
    // genome/finalist-reranking): the DP keeps its per-window surrogate
    // objective AND its selection (the tracks above are untouched); the MODEL
    // OF RECORD re-ranks the enumerated finalists. The enumeration order is
    // the SURROGATE's own (the ruling's k-best DP paths); the model's
    // refund per enumerated finalist is measured (surrogate - model) and the
    // stop rule uses its running maximum: no un-enumerated chain (surrogate
    // at or above the frontier) can beat the re-ranked best unless its refund
    // exceeds every measured one — reported, not assumed. No new constants:
    // the finalist universe's widening uses the SAME established swing
    // statistic (margins - min_viable_loss) applied to the surrogate's own
    // loss table at the WIDENED threshold the misranking scale requires —
    // measured on the predecessor's run: the donor-best alleles sit
    // +29.5k..+66.0k above their loci's surrogate winners (the tract middle),
    // two orders above the surrogate swing (0.5k-12k), and the predecessor's
    // retention EXCLUDED the donor best alleles at loci 2,3,5,6,7,10,11,13,14.
    // The widening therefore admits, per locus, every scorable allele whose
    // surrogate loss is finite — the finalist universe is the full scorable
    // haploid candidate set (the per-locus surrogate narrowing is what the
    // re-ranking supersedes; the DP's own tracks keep their tables).
    let finalist_state_sets: Vec<HashSet<usize>> = (0..locus_count)
        .map(|locus| {
            let mut set = retained_haploid_sets[locus].clone();
            for (allele, &loss) in haploid_losses[locus].iter().enumerate() {
                if loss.is_finite() {
                    set.insert(allele);
                }
            }
            set
        })
        .collect();
    eprintln!(
        "[phasing] finalist states per locus: {:?}",
        finalist_state_sets
            .iter()
            .map(|set| set.len())
            .collect::<Vec<_>>()
    );
    let finalist_union_sets: Vec<HashSet<usize>> = finalist_state_sets.clone();
    let locus_owners_all: Vec<Vec<u32>> = (0..locus_count)
        .map(|locus| {
            combined[locus]
                .iter()
                .map(crate::traversal_owner)
                .collect::<Vec<_>>()
        })
        .collect();
    let transition_drafts: Vec<BoundaryTransitionDraft> = (0..locus_count - 1)
        .into_par_iter()
        .map(|boundary| {
            build_boundary_transition_draft(
                panel,
                sources,
                flank_memo,
                model,
                span,
                path_of_source,
                &locus_owners_all[boundary],
                &locus_owners_all[boundary + 1],
                component_locus_to_partition,
                &combined[boundary],
                &combined[boundary + 1],
                &retained_union_sets[boundary],
                &retained_union_sets[boundary + 1],
            )
        })
        .collect::<io::Result<_>>()?;
    // The FINALIST track's own boundary drafts over the finalist state sets
    // (the L DP's transitions and the re-ranked route's surrogate/boundary
    // lookups; the re-ranked chain's alleles are finalist states by
    // construction). Built separately from the tracks' drafts: a
    // composition's charge pools over its realizing pairs' owners and
    // junction events, so widening the pair universe can change a SHARED
    // composition's value — the tracks' matrices must stay bit-identical.
    let finalist_transition_started = Instant::now();
    let finalist_transition_drafts: Vec<BoundaryTransitionDraft> = (0..locus_count - 1)
        .into_par_iter()
        .map(|boundary| {
            build_boundary_transition_draft(
                panel,
                sources,
                flank_memo,
                model,
                span,
                path_of_source,
                &locus_owners_all[boundary],
                &locus_owners_all[boundary + 1],
                component_locus_to_partition,
                &combined[boundary],
                &combined[boundary + 1],
                &finalist_union_sets[boundary],
                &finalist_union_sets[boundary + 1],
            )
        })
        .collect::<io::Result<_>>()?;
    eprintln!(
        "[phasing] finalist boundary drafts: {:.2}s",
        finalist_transition_started.elapsed().as_secs_f64()
    );
    // Extend the multiplicity backgrounds for every pooled transition
    // composition feature the scan has not measured (the phasing boundary
    // compositions enumerate retained-pair juxtapositions beyond the spine
    // boundary links).
    {
        let mut universe: BTreeSet<FeatureKey> = BTreeSet::new();
        for draft in transition_drafts.iter().chain(finalist_transition_drafts.iter()) {
            for key in draft.feature_set() {
                universe.insert(key);
            }
        }
        let features: Vec<FeatureKey> = universe.into_iter().collect();
        backgrounds.scan_extend(panel, features, *k, model)?;
    }
    let backgrounds_shared: &FeatureBackgrounds = &**backgrounds;
    let transition_costs: Vec<BoundaryTransitionCosts> = transition_drafts
        .into_par_iter()
        .enumerate()
        .map(|(boundary, draft)| {
            finalize_boundary_transition_costs(
                draft,
                routed_equal,
                component_locus_to_partition,
                model,
                backgrounds_shared,
                sample_backgrounds,
            )
        })
        .collect::<io::Result<_>>()?;
    let finalist_transition_costs: Vec<BoundaryTransitionCosts> = finalist_transition_drafts
        .into_par_iter()
        .enumerate()
        .map(|(boundary, draft)| {
            finalize_boundary_transition_costs(
                draft,
                routed_equal,
                component_locus_to_partition,
                model,
                backgrounds_shared,
                sample_backgrounds,
            )
        })
        .collect::<io::Result<_>>()?;
    let boundary_floors: Vec<f64> = transition_costs
        .iter()
        .map(|costs| costs.cost.iter().copied().fold(f64::INFINITY, f64::min))
        .collect();
    let boundary_floors_haploid: Vec<f64> = transition_costs
        .iter()
        .map(|costs| {
            costs
                .cost_haploid
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    ensure(
        boundary_floors
            .iter()
            .chain(boundary_floors_haploid.iter())
            .all(|value| value.is_finite()),
        "phasing boundary floor over an unscored matrix",
    )?;
    let mut suffix_from_locus = vec![0.0f64; locus_count];
    {
        let mut local_suffix = 0.0f64;
        let mut boundary_suffix = 0.0f64;
        for locus in (0..locus_count).rev() {
            suffix_from_locus[locus] = local_suffix + 2.0 * boundary_suffix;
            local_suffix += sweeps[locus].floor;
            if locus > 0 {
                boundary_suffix += boundary_floors[locus - 1];
            }
        }
    }
    let costs_seconds = costs_started.elapsed().as_secs_f64();
    eprintln!(
        "[phasing] transition costs done: {} boundary matrices, {:.2}s",
        transition_costs.len(),
        costs_seconds
    );
    rss_probe(rss, "phasing_after_costs")?;

    // The native-backbone incumbent under the phasing transition model
    // (must equal the spine's incumbent up to f64 summation order: the
    // backbone's transitions are gap-0 co-occurrences).
    let mut incumbent_phasing = 0.0f64;
    for locus in 0..locus_count {
        incumbent_phasing += sweeps[locus]
            .native_pair_loss
            .ok_or_else(|| invalid("phasing backbone pair loss missing"))?;
    }
    for boundary in 0..locus_count - 1 {
        let costs = &transition_costs[boundary];
        let row = costs.left_index[backbone_chain[boundary]] as usize * costs.right_count as usize;
        let column = costs.right_index[backbone_chain[boundary + 1]] as usize;
        incumbent_phasing += 2.0 * costs.cost[row + column];
    }
    let incumbent_delta = incumbent_phasing - incumbent_spine;
    eprintln!(
        "[phasing] incumbent (phasing) {incumbent_phasing:.2} vs spine {incumbent_spine:.2} (delta {incumbent_delta:.4})"
    );

    // ---------------------------------------- P4: the phasing chain DP.
    // Both ploidy tracks (owner decision 2026-09-23): the diploid ordered
    // pair states, and the haploid states whose second slot is the EMPTY
    // sentinel over the same margin-retained member universe.
    let per_locus_builds: Vec<(LocusStateTable, LocusStateTable, Vec<f64>)> = (0..locus_count)
        .into_par_iter()
        .map(|locus| {
            let states = ordered_pair_states(
                locus,
                &locus_classes[locus],
                &sweeps[locus],
                &membership_slices[locus],
                backbone_chain[locus],
                margins[locus],
            );
            let haploid_states = ordered_haploid_states(
                &haploid_losses[locus],
                &retained_haploid_sets[locus],
            )?;
            let haploid_class_loss = haploid_losses[locus].clone();
            let diploid_table = if locus == 0 {
                LocusStateTable::build(
                    states
                        .into_iter()
                        .map(|(pair, loss)| (pair, 0, 0, loss))
                        .collect(),
                    1,
                )
            } else {
                let costs = &transition_costs[locus - 1];
                let right_index = &costs.right_index;
                LocusStateTable::build(
                    states
                        .into_iter()
                        .map(|(pair, loss)| {
                            let first = right_index[pair[0]];
                            let second = right_index[pair[1]];
                            (pair, first, second, loss)
                        })
                        .collect(),
                    costs.right_count as usize,
                )
            };
            // The haploid table's `second` column is unused by the haploid
            // DP (no second-slot boundary charge); it stays 0.
            let haploid_table = if locus == 0 {
                LocusStateTable::build(
                    haploid_states
                        .into_iter()
                        .map(|(pair, loss)| (pair, 0, 0, loss))
                        .collect(),
                    1,
                )
            } else {
                let costs = &transition_costs[locus - 1];
                let right_index = &costs.right_index;
                LocusStateTable::build(
                    haploid_states
                        .into_iter()
                        .map(|(pair, loss)| {
                            let first = right_index[pair[0]];
                            (pair, first, 0, loss)
                        })
                        .collect(),
                    costs.right_count as usize,
                )
            };
            Ok((diploid_table, haploid_table, haploid_class_loss))
        })
        .collect::<io::Result<_>>()?;
    let mut per_locus_tables = Vec::with_capacity(locus_count);
    let mut per_locus_haploid_tables = Vec::with_capacity(locus_count);
    let mut haploid_loss_tables = Vec::with_capacity(locus_count);
    for (diploid, haploid, losses) in per_locus_builds {
        per_locus_tables.push(diploid);
        per_locus_haploid_tables.push(haploid);
        haploid_loss_tables.push(losses);
    }
    // The haploid track's suffix bound: the EXACT per-locus minimum of the
    // per-allele oracle record-once losses plus the exact per-boundary minimum
    // of the haploid cost matrix (the panel pair floor is NOT admissible
    // for the sample-side charges — the sample-side credit can lie far
    // below the panel-beta box — so the haploid track carries its own
    // tight floors; the haploid DP pays one boundary charge per
    // boundary).
    let haploid_local_floors: Vec<f64> = haploid_loss_tables
        .iter()
        .map(|losses| {
            losses
                .iter()
                .copied()
                .filter(|value| value.is_finite())
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    ensure(
        haploid_local_floors.iter().all(|value| value.is_finite()),
        "phasing haploid local floor over unscored alleles",
    )?;
    let mut suffix_haploid_from_locus = vec![0.0f64; locus_count];
    {
        let mut local_suffix = 0.0f64;
        let mut boundary_suffix = 0.0f64;
        for locus in (0..locus_count).rev() {
            suffix_haploid_from_locus[locus] = local_suffix + boundary_suffix;
            local_suffix += haploid_local_floors[locus];
            if locus > 0 {
                boundary_suffix += boundary_floors_haploid[locus - 1];
            }
        }
    }
    // The native-backbone HAPLOID incumbent (slot 2 empty): the backbone
    // allele's single-class local losses plus its own boundary charges
    // once. Admissible-chain score — tightens the haploid track's prune
    // threshold alongside the diploid incumbent.
    let mut incumbent_haploid = 0.0f64;
    for locus in 0..locus_count {
        incumbent_haploid += haploid_loss_tables[locus][backbone_chain[locus]];
    }
    for boundary in 0..locus_count - 1 {
        let costs = &transition_costs[boundary];
        let row = costs.left_index[backbone_chain[boundary]] as usize * costs.right_count as usize;
        let column = costs.right_index[backbone_chain[boundary + 1]] as usize;
        incumbent_haploid += costs.cost_haploid[row + column];
    }
    ensure(
        incumbent_haploid.is_finite(),
        "haploid incumbent over an unscored backbone class",
    )?;
    // The omission charge (Fix 2) split the tracks' absolute scales: the
    // haploid class losses now carry the windows' omitted-observed costs,
    // so the DIPLOID incumbent (its model unchanged) is no longer a valid
    // prune threshold for the haploid track — pruning against it would
    // empty the haploid initial layer (measured: diploid 213,330.62 vs
    // haploid backbone 1,173,034.98 on the first omission run). The
    // haploid track prunes against its OWN model's incumbent (the backbone
    // chain's exact haploid score — admissible); the ploidy selection
    // below still compares the two tracks' optima.
    let haploid_prune_incumbent = incumbent_haploid;
    eprintln!(
        "[phasing] incumbents: diploid {incumbent_phasing:.2} haploid {incumbent_haploid:.2}"
    );
    let dp_diploid = run_phasing_chain_dp(
        locus_count,
        &per_locus_tables,
        &transition_costs,
        &suffix_from_locus,
        incumbent_phasing,
        false,
        rss,
    )?;
    eprintln!(
        "[phasing] diploid DP done: best {:.2}, states {:?}, {:.2}s",
        dp_diploid.best_score,
        dp_diploid.state_counts,
        dp_diploid.wall_seconds
    );
    rss_probe(rss, "phasing_after_dp_diploid")?;
    let dp_haploid = run_phasing_chain_dp(
        locus_count,
        &per_locus_haploid_tables,
        &transition_costs,
        &suffix_haploid_from_locus,
        haploid_prune_incumbent,
        true,
        rss,
    )?;
    eprintln!(
        "[phasing] haploid DP done: best {:.2}, states {:?}, {:.2}s",
        dp_haploid.best_score,
        dp_haploid.state_counts,
        dp_haploid.wall_seconds
    );
    rss_probe(rss, "phasing_after_dp_haploid")?;

    // -------------------------------- P4b: the k-best finalist enumeration +
    // the model-of-record re-rank (supervisor ruling, genome/finalist-reranking).
    // The DP keeps its per-window surrogate objective and its selection (the
    // guard below compares); the FINALIST set is enumerated from the L DP
    // (the surrogate's chain-decomposable part: the per-window loss minus the
    // allele's own in-window omission addend — the model's per-window refund
    // through the omission channel is bounded by that addend) and re-ranked
    // EXACTLY under the model of record (`haploid_chain_level_m1`, the
    // selected_chain_level_m1 evaluator). The re-ranked best becomes the
    // haploid track's representative.
    let finalists_started = Instant::now();
    let per_locus_finalist_tables: Vec<LocusStateTable> = (0..locus_count)
        .map(|locus| {
            let states = ordered_haploid_states(
                &haploid_losses[locus],
                &finalist_state_sets[locus],
            )?;
            if locus == 0 {
                Ok(LocusStateTable::build(
                    states
                        .into_iter()
                        .map(|(pair, loss)| (pair, 0, 0, loss))
                        .collect(),
                    1,
                ))
            } else {
                let costs = &finalist_transition_costs[locus - 1];
                let right_index = &costs.right_index;
                Ok(LocusStateTable::build(
                    states
                        .into_iter()
                        .map(|(pair, loss)| {
                            (pair, right_index[pair[0]], 0, loss)
                        })
                        .collect(),
                    costs.right_count as usize,
                ))
            }
        })
        .collect::<io::Result<_>>()?;
    // The pairwise-coupling search objective's boundary matrices (the
    // enumeration order; see build_coupling_transition_costs). The DP's local
    // losses stay the SURROGATE's (the tables above); the coupling lives in
    // the transitions.
    let coupling_started = Instant::now();
    let sorted_finalist_states: Vec<Vec<usize>> = per_locus_finalist_tables
        .iter()
        .map(|table| table.rows.iter().map(|row| row.0[0]).collect())
        .collect();
    // Operational probe-only mode (IMPG_FINALIST_SKIP_COUPLING=1): the
    // coupling matrices are SKIPPED (the enumerator has no order; the
    // probes' model totals and the instance-probe diagnostic do not need
    // them). Diagnostic guard, not a model path — a run that QUOTES
    // coupling numbers must not set it.
    let skip_coupling = std::env::var("IMPG_FINALIST_SKIP_COUPLING")
        .ok()
        .as_deref()
        == Some("1");
    let coupling_started = Instant::now();
    let coupling_transition_costs: Vec<BoundaryTransitionCosts> = if skip_coupling {
        eprintln!("[phasing] coupling transition matrices: SKIPPED (probe-only mode)");
        Vec::new()
    } else {
    (0..locus_count - 1)
        .into_par_iter()
        .map(|boundary| {
            let terms_left: Vec<f64> = (0..window_obs_index[boundary].features.len())
                .map(|i| window_obs_index[boundary].terms[i])
                .collect();
            let terms_right: Vec<f64> = (0..window_obs_index[boundary + 1].features.len())
                .map(|i| window_obs_index[boundary + 1].terms[i])
                .collect();
            build_coupling_transition_costs(
                &sorted_finalist_states[boundary],
                &sorted_finalist_states[boundary + 1],
                &finalist_transition_costs[boundary],
                &haploid_exons[boundary],
                &haploid_exons[boundary + 1],
                &terms_left,
                &terms_right,
                window_instances.is_some(),
                Some(sample_backgrounds),
                &window_obs_index[boundary].features,
                &window_obs_index[boundary + 1].features,
                &window_instance_tables[boundary].0,
                &window_instance_tables[boundary + 1].0,
                &window_instance_tables[boundary].1,
                &window_instance_tables[boundary + 1].1,
            )
        })
        .collect::<Vec<_>>()
    };
    eprintln!(
        "[phasing] coupling transition matrices: {:.2}s",
        coupling_started.elapsed().as_secs_f64()
    );
    let finalist_local_floors: Vec<f64> = per_locus_finalist_tables
        .iter()
        .map(|table| {
            table
                .rows
                .iter()
                .map(|&(_, _, _, loss)| loss)
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    ensure(
        finalist_local_floors.iter().all(|value| value.is_finite()),
        "finalist local floor over an empty finalist table",
    )?;
    let finalist_boundary_floors: Vec<f64> = if coupling_transition_costs.is_empty() {
        // The probe-only mode: the floors are unreachable (the enumeration
        // is skipped below); zeros keep the layout.
        eprintln!("[phasing] finalist boundary floors: SKIPPED (probe-only mode)");
        vec![0.0; locus_count.saturating_sub(1)]
    } else {
        coupling_transition_costs
            .iter()
            .map(|costs| {
                costs
                    .cost_haploid
                    .iter()
                    .copied()
                    .fold(f64::INFINITY, f64::min)
            })
            .collect()
    };
    ensure(
        finalist_boundary_floors.iter().all(|value| value.is_finite()),
        "finalist boundary floor over an unscored matrix",
    )?;
    let mut suffix_finalist_from_locus = vec![0.0f64; locus_count];
    {
        let mut local_suffix = 0.0f64;
        let mut boundary_suffix = 0.0f64;
        for locus in (0..locus_count).rev() {
            suffix_finalist_from_locus[locus] = local_suffix + boundary_suffix;
            local_suffix += finalist_local_floors[locus];
            if locus > 0 {
                boundary_suffix += finalist_boundary_floors[locus - 1];
            }
        }
    }
    // The FINALIST DP: the surrogate's local losses over the widened finalist
    // tables with the COUPLING transitions (the search objective C = the
    // surrogate minus the adjacent-window exoneration channel — an
    // upper-bound form on the model of record), NO incumbent prune (the
    // enumerator's own stop rule bounds the work). The k-best enumeration
    // walks its layers in C's own total order.
    let dp_finalist = if skip_coupling {
        // The probe-only mode: the finalist DP orders the enumeration the
        // probe-only run does not perform; the layers are unreachable.
        None
    } else {
        Some(run_phasing_chain_dp(
            locus_count,
            &per_locus_finalist_tables,
            &coupling_transition_costs,
            &suffix_finalist_from_locus,
            f64::INFINITY,
            true,
            rss,
        )?)
    };
    if let Some(outcome) = &dp_finalist {
        eprintln!(
            "[phasing] finalist DP done: best {:.2}, states {:?}, {:.2}s",
            outcome.best_score,
            outcome.state_counts,
            outcome.wall_seconds
        );
    }
    rss_probe(rss, "phasing_after_dp_finalist")?;

    // The probe routes (repeatable --probe-route LABEL:PATH.json): previously
    // measured chains mapped into the finalist tables — the re-rank's
    // bit-identity anchors (the predecessor's selected chain and
    // donor-probe chain) and the surrogate-rank measurement the ruling asks
    // for. Their model totals come from the SAME evaluator as every
    // finalist's.
    let parse_probe_route = |path: &std::path::Path| -> io::Result<Vec<usize>> {
        let value: serde_json::Value = read_json(path)?;
        let rows = value
            .as_array()
            .ok_or_else(|| invalid("probe route must be an array"))?;
        ensure(
            rows.len() == locus_count,
            "probe route cardinality mismatch",
        )?;
        let mut alleles = Vec::with_capacity(locus_count);
        for (locus, row) in rows.iter().enumerate() {
            let slots = row
                .as_array()
                .ok_or_else(|| invalid("probe route row must be an array"))?;
            ensure(
                slots.len() == 2 && slots[1].as_str() == Some("empty-slot2"),
                "probe route row must be a haploid slot-0/empty-slot2 pair",
            )?;
            let segments_value = slots[0]
                .get("segments")
                .and_then(serde_json::Value::as_array)
                .ok_or_else(|| invalid("probe route row lacks segments"))?;
            let mut probe_segments: Vec<(usize, u64, u64, bool)> = Vec::new();
            for segment in segments_value {
                probe_segments.push((
                    segment
                        .get("source")
                        .and_then(serde_json::Value::as_u64)
                        .ok_or_else(|| invalid("probe segment lacks source"))? as usize,
                    segment
                        .get("start")
                        .and_then(serde_json::Value::as_u64)
                        .ok_or_else(|| invalid("probe segment lacks start"))?,
                    segment
                        .get("end")
                        .and_then(serde_json::Value::as_u64)
                        .ok_or_else(|| invalid("probe segment lacks end"))?,
                    segment
                        .get("reverse")
                        .and_then(serde_json::Value::as_bool)
                        .ok_or_else(|| invalid("probe segment lacks reverse"))?,
                ));
            }
            let matched = combined[locus].iter().position(|candidate| {
                candidate.segments.len() == probe_segments.len()
                    && candidate
                        .segments
                        .iter()
                        .zip(probe_segments.iter())
                        .all(
                            |(segment, &(source, start, end, reverse))| {
                                segment.source == source
                                    && segment.start == start
                                    && segment.end == end
                                    && segment.reverse == reverse
                            },
                        )
            });
            alleles.push(matched.ok_or_else(|| {
                invalid(&format!(
                    "probe route row at locus {locus} matches no candidate allele"
                ))
            })?);
        }
        Ok(alleles)
    };
    let mut probe_states: Vec<(String, Option<Vec<usize>>)> = Vec::new();
    for (label, path) in probe_specs.iter() {
        let parsed = parse_probe_route(path).ok();
        if parsed.is_none() {
            eprintln!("[phasing] probe route {label}: MAPPED FAILED (skipped)");
        }
        probe_states.push((label.clone(), parsed));
    }
    // The chain's SEARCH-objective (coupling) total under the finalist
    // tables: the surrogate's local losses plus the COUPLING transitions
    // (path-order accumulation — the DP's own form).
    let chain_coupling = |alleles: &[usize]| -> f64 {
        if coupling_transition_costs.is_empty() {
            return 0.0; // the probe-only mode: no coupling matrices built
        }
        let mut total = 0.0f64;
        for (locus, &allele) in alleles.iter().enumerate() {
            total += haploid_losses[locus][allele];
            if locus > 0 {
                let costs = &coupling_transition_costs[locus - 1];
                total += costs.cost_haploid[costs.left_index[alleles[locus - 1]] as usize
                    * costs.right_count as usize
                    + costs.right_index[allele] as usize];
            }
        }
        total
    };
    // The chain's true SURROGATE total under the finalist tables (path-order
    // accumulation — the DP's own form).
    let chain_surrogate = |alleles: &[usize]| -> f64 {
        let mut surrogate = 0.0f64;
        for (locus, &allele) in alleles.iter().enumerate() {
            surrogate += haploid_losses[locus][allele];
            if locus > 0 {
                let costs = &finalist_transition_costs[locus - 1];
                surrogate += costs.cost_haploid[costs.left_index[alleles[locus - 1]] as usize
                    * costs.right_count as usize
                    + costs.right_index[allele] as usize];
            }
        }
        surrogate
    };
    let chain_in_finalist_union = |alleles: &[usize]| -> bool {
        alleles.iter().enumerate().all(|(locus, &allele)| {
            finalist_state_sets[locus].contains(&allele)
        })
    };
    // Evaluate one chain under the model of record (the shared piece memo
    // keeps the repeated evaluations cheap; values are memo-independent).
    let mut evaluate_chain = |alleles: &[usize],
                          memo: &mut HashMap<(usize, u64, u64), Profile>,
                          spelled_memo: &mut crate::SpelledMemo,
                          spelled_out: Option<
                              &mut HashMap<FeatureKey, Vec<crate::SpelledSpan>>,
                          >|
     -> io::Result<f64> {
        let chain: Vec<&genome::SpanningTraversal> = alleles
            .iter()
            .enumerate()
            .map(|(locus, &allele)| &combined[locus][allele])
            .collect();
        haploid_chain_level_m1(
            panel,
            sources,
            territory,
            routed_equal,
            site,
            window_obs,
            window_instances,
            component_locus_to_partition,
            locus_count,
            locus_offset,
            model,
            span,
            path_of_source,
            &locus_partitions,
            &mut **backgrounds,
            sample_backgrounds,
            &chain,
            Some(memo),
            Some(spelled_memo),
            spelled_out,
        )
    };
    let mut shared_piece_memo: HashMap<(usize, u64, u64), Profile> = HashMap::new();
    let mut shared_spelled_memo: crate::SpelledMemo = HashMap::new();
    // The probes first (their anchor values are reported verbatim).
    let mut probe_rows: Vec<serde_json::Value> = Vec::new();
    let mut probe_alleles: Vec<Option<Vec<usize>>> = Vec::new();
    for (label, parsed) in &probe_states {
        let row = match parsed {
            Some(alleles) if chain_in_finalist_union(alleles) => {
                let surrogate = chain_surrogate(alleles);
                let coupling = chain_coupling(alleles);
                let mut probe_spelled: HashMap<FeatureKey, Vec<crate::SpelledSpan>> =
                    HashMap::new();
                let spelled_out = (!instance_probe_features.is_empty())
                    .then_some(&mut probe_spelled as &mut _);
                let model = evaluate_chain(
                    alleles,
                    &mut shared_piece_memo,
                    &mut shared_spelled_memo,
                    spelled_out,
                )?;
                let instance_probe = (!instance_probe_features.is_empty()).then(|| {
                    crate::instance_probe_report(
                        window_obs,
                        window_instances.unwrap_or(&crate::InstanceStructure::default()),
                        &probe_spelled,
                        sample_backgrounds,
                        &instance_probe_features,
                    )
                });
                probe_alleles.push(Some(alleles.clone()));
                serde_json::json!({
                    "label": label,
                    "mapped": true,
                    "in_finalist_universe": true,
                    "surrogate_total": surrogate,
                    "coupling_search_total": coupling,
                    "coupling_minus_model": coupling - model,
                    "model_of_record_total": model,
                    "refund": surrogate - model,
                    "instance_probe": instance_probe,
                })
            }
            Some(alleles) => {
                let mut probe_spelled: HashMap<FeatureKey, Vec<crate::SpelledSpan>> =
                    HashMap::new();
                let spelled_out = (!instance_probe_features.is_empty())
                    .then_some(&mut probe_spelled as &mut _);
                let model = evaluate_chain(
                    alleles,
                    &mut shared_piece_memo,
                    &mut shared_spelled_memo,
                    spelled_out,
                )?;
                let instance_probe = (!instance_probe_features.is_empty()).then(|| {
                    crate::instance_probe_report(
                        window_obs,
                        window_instances.unwrap_or(&crate::InstanceStructure::default()),
                        &probe_spelled,
                        sample_backgrounds,
                        &instance_probe_features,
                    )
                });
                probe_alleles.push(Some(alleles.clone()));
                serde_json::json!({
                    "label": label,
                    "mapped": true,
                    "in_finalist_universe": false,
                    "surrogate_total": serde_json::Value::Null,
                    "model_of_record_total": model,
                    "refund": serde_json::Value::Null,
                    "instance_probe": instance_probe,
                })
            }
            None => {
                probe_alleles.push(None);
                serde_json::json!({ "label": label, "mapped": false })
            }
        };
        probe_rows.push(row);
    }
    // The enumeration + exact re-rank loop.
    let finalist_cap: usize = std::env::var("IMPG_FINALIST_CAP")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(2000);
    if skip_coupling {
        eprintln!("[phasing] finalist enumeration SKIPPED (probe-only mode)");
        // The probes' rows above carry the run's model-of-record evidence;
        // the enumeration/re-rank machinery needs the coupling matrices.
        let mut skipped_finalists = serde_json::json!({
            "skipped": "IMPG_FINALIST_SKIP_COUPLING=1 (probe-only mode)",
            "probe_rows": probe_rows,
            "reranked_best": serde_json::Value::Null,
        });
        skipped_finalists["coupling_skipped"] = serde_json::Value::Bool(true);
        return Ok(serde_json::json!({
            "model": "local-first-spine-v1-correlation-phasing",
            "finalists": skipped_finalists,
            "probe_only": true,
        }));
    }
    let mut enumerator = FinalistEnumerator::new(
            &dp_finalist
                .as_ref()
                .expect("the enumerator implies the finalist DP (probe-only mode skips it)")
                .layers,
            &coupling_transition_costs,
        );
    let mut w_best = f64::INFINITY;
    // The measured maximum refund (surrogate - model) over the enumerated
    // finalists: the stop rule's slack. An un-enumerated chain's surrogate is
    // at or above the frontier; it can beat W only with a refund larger than
    // every measured one (reported, not assumed).
    let mut max_refund = f64::NEG_INFINITY;
    // The containment target (design note (a)): the enumerated rank of the
    // best-known chain among the probes (the chain the k must CONTAIN).
    let mut containment_target_rank: Option<usize> = None;
    let mut best_probe_model = f64::INFINITY;
    let mut finalist_rows: Vec<(usize, f64, f64, f64, f64)> = Vec::new();
    let mut reranked_best: Option<(usize, f64, f64, Vec<[usize; 2]>)> = None;
    let mut surrogate_best_model = f64::INFINITY;
    let mut stopped_by = "exhausted";
    let mut rank = 0usize;
    let mut dry_run = std::env::var("IMPG_FINALIST_DRY_RUN").is_ok();
    loop {
        let Some((coupling_total, allele_route)) = enumerator.next_chain() else {
            stopped_by = "exhausted";
            break;
        };
        let allele_path: Vec<usize> = allele_route.iter().map(|pair| pair[0]).collect();
        for (index, probe) in probe_alleles.iter().enumerate() {
            if let Some(probe_allele_path) = probe {
                if *probe_allele_path == allele_path {
                    if let Some(row) = probe_rows.get_mut(index) {
                        row["enumerated_rank"] = serde_json::json!(rank);
                    }
                }
            }
        }
        if rank >= finalist_cap {
            stopped_by = "cap";
            break;
        }
        // The containment stop (design note (a)): the k-best prefix must
        // contain the donor-rich chain (the best-known probe chain) AND be
        // widened past it (its rank doubled) so a rival with a still-larger
        // refund has room to appear; every popped chain in between is
        // re-ranked exactly.
        if let Some(target) = containment_target_rank {
            if rank >= 2 * target {
                stopped_by = "containment_doubled";
                break;
            }
        }
        if dry_run {
            // Density measurement only: no model evaluation; W stays open so
            // the frontier runs to the cap (the run's log reports the
            // coupling frontier's drift and the probe ranks).
            finalist_rows.push((rank, coupling_total, f64::NAN, f64::NAN, f64::NAN));
            rank += 1;
            if rank % 200 == 0 {
                eprintln!(
                    "[phasing] finalists(dry): {} enumerated, coupling frontier {:.2} ({:.1}s)",
                    rank,
                    coupling_total,
                    finalists_started.elapsed().as_secs_f64()
                );
            }
            continue;
        }
        let surrogate_total = chain_surrogate(&allele_path);
        let model_total = evaluate_chain(&allele_path, &mut shared_piece_memo, &mut shared_spelled_memo, None)?;
        if rank == 0 {
            surrogate_best_model = model_total;
        }
        // Track the best-known probe chain's rank (the containment target).
        for (index, probe) in probe_alleles.iter().enumerate() {
            if let Some(probe_allele_path) = probe {
                if *probe_allele_path == allele_path {
                    if let Some(row) = probe_rows.get_mut(index) {
                        let model_value = row["model_of_record_total"].as_f64();
                        if model_value.is_some_and(|value| value < best_probe_model) {
                            best_probe_model = model_value.expect("checked");
                            containment_target_rank = Some(rank);
                        }
                    }
                }
            }
        }
        let refund = surrogate_total - model_total;
        if refund > max_refund {
            max_refund = refund;
        }
        // The search objective's error accounting (ruling requirement 5):
        // C - model = the model's LONG-RANGE exoneration channel the pairwise
        // approximation does not catch, plus the documented form residual.
        let coupling_minus_model = coupling_total - model_total;
        if model_total < w_best {
            w_best = model_total;
            reranked_best = Some((rank, surrogate_total, model_total, allele_route.clone()));
        }
        finalist_rows.push((rank, coupling_total, surrogate_total, model_total, coupling_minus_model));
        rank += 1;
        if rank % 50 == 0 {
            eprintln!(
                "[phasing] finalists: {} enumerated, best model {:.2}, coupling frontier {:.2}, max refund {:.2} ({:.1}s)",
                rank,
                w_best,
                coupling_total,
                max_refund,
                finalists_started.elapsed().as_secs_f64()
            );
        }
    }
    // The coupling frontier at the stop (the next un-enumerated chain's
    // search-objective total).
    let surrogate_frontier = enumerator.next_chain().map(|(value, _)| value);
    // Completeness check: could an allele OUTSIDE the finalist universe lie
    // on a chain that beats the re-ranked best? A chain through an excluded
    // allele needs a surrogate total at most W + max_refund (its model is its
    // surrogate minus a refund; refunds above the measured maximum are the
    // reported caveat). The chain's other loci contribute at least their
    // ALL-allele finite loss minima; the transitions contribute at least the
    // finalist matrices' minima (an excluded allele's own transitions were
    // never scored — documented edge; the transition scale is +/-0.5k, three
    // orders below the refund scale).
    let mut completeness_violations: Vec<serde_json::Value> = Vec::new();
    if w_best.is_finite() {
    {
        let local_minima: Vec<f64> = (0..locus_count)
            .map(|locus| {
                haploid_losses[locus]
                    .iter()
                    .copied()
                    .filter(|value| value.is_finite())
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        let boundary_min_total: f64 = finalist_boundary_floors.iter().copied().sum();
        for locus in 0..locus_count {
            let others: f64 = (0..locus_count)
                .filter(|&j| j != locus)
                .map(|j| local_minima[j])
                .sum::<f64>()
                + boundary_min_total;
            for (allele, &loss) in haploid_losses[locus].iter().enumerate() {
                if loss.is_finite() && !finalist_state_sets[locus].contains(&allele) {
                    if others + loss <= w_best + max_refund.max(0.0) {
                        completeness_violations.push(serde_json::json!({
                            "locus": locus,
                            "allele": allele,
                            "surrogate_loss": loss,
                        }));
                    }
                }
            }
        }
    }
    }
    let enumerator_seed_delta = enumerator.max_seed_delta;
    let finalists_seconds = finalists_started.elapsed().as_secs_f64();
    eprintln!(
        "[phasing] finalists done: {} enumerated ({}), best model {:.2}, {:.2}s",
        rank,
        stopped_by,
        w_best,
        finalists_seconds
    );
    let finalist_summary = {
        let count = finalist_rows.len();
        let shown: Vec<serde_json::Value> = finalist_rows
            .iter()
            .take(16)
            .chain(
                finalist_rows
                    .iter()
                    .rev()
                    .take(8)
                    .collect::<Vec<_>>()
                    .into_iter()
                    .rev(),
            )
            .map(|&(rank, coupling_total, surrogate_total, model_total, coupling_minus_model)| {
                serde_json::json!({
                    "rank_coupling": rank,
                    "coupling_search_total": coupling_total,
                    "surrogate_total": surrogate_total,
                    "model_of_record_total": model_total,
                    "coupling_minus_model": coupling_minus_model,
                })
            })
            .collect();
        let coupling_values: Vec<f64> = finalist_rows.iter().map(|row| row.1).collect();
        let surrogate_values: Vec<f64> = finalist_rows.iter().map(|row| row.2).collect();
        let model_values: Vec<f64> = finalist_rows.iter().map(|row| row.3).collect();
        let coupling_minus_model_values: Vec<f64> = finalist_rows.iter().map(|row| row.4).collect();
        // Ruling requirement 4: the TOP-SLICE verification — is the search
        // order enriched for good model totals?
        let top_slice: Vec<f64> = model_values.iter().copied().take(200).collect();
        let top_slice_summary = if top_slice.len() >= 200 {
            serde_json::json!({
                "count": 200,
                "model_min": top_slice.iter().copied().fold(f64::INFINITY, f64::min),
                "model_max": top_slice.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            })
        } else {
            serde_json::json!({ "count": top_slice.len() })
        };
        serde_json::json!({
            "role_separation": "the coupling objective is the SEARCH ordering ONLY (the k-best enumeration order); the selection and every reported score are exact model-of-record numbers from the re-ranker",
            "objective": "SEARCH: C = the surrogate's per-window losses minus the ADJACENT-window exoneration channel (per (boundary, allele-pair), priced with the model of record's own per-feature omission terms) - an upper-bound form on the model of record. MODEL: every enumerated finalist is re-ranked EXACTLY under the chain-level once-per-material model of record",
            "stop_rule": "enumerate the surrogate's k-best until the best-known probe chain's rank is contained AND doubled (design note (a)'s containment criterion, with a widening margin for rivals with still-larger refunds); every popped chain in between is re-ranked exactly; an un-enumerated chain can beat W only with a refund above every measured one (reported, not assumed)",
            "operational_cap": finalist_cap,
            "enumerated": count,
            "stopped_by": stopped_by,
            "containment_target_rank": containment_target_rank,
            "finalist_dp_best": dp_finalist.as_ref().expect("the finalists report implies the finalist DP").best_score,
            "enumerator_seed_delta_max": enumerator_seed_delta,
            "coupling_frontier_at_stop": surrogate_frontier,
            "w_best_model": if w_best.is_finite() { serde_json::json!(w_best) } else { serde_json::Value::Null },
            "max_refund_measured": if max_refund.is_finite() { serde_json::json!(max_refund) } else { serde_json::Value::Null },
            "gap_distribution": {
                "coupling_min": coupling_values.iter().copied().fold(f64::INFINITY, f64::min),
                "coupling_max": coupling_values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                "surrogate_min": surrogate_values.iter().copied().fold(f64::INFINITY, f64::min),
                "surrogate_max": surrogate_values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                "model_min": model_values.iter().copied().fold(f64::INFINITY, f64::min),
                "model_max": model_values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                "coupling_minus_model_min": coupling_minus_model_values.iter().copied().fold(f64::INFINITY, f64::min),
                "coupling_minus_model_max": coupling_minus_model_values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            },
            "top_slice_model_of_record": top_slice_summary,
            "reranked_best": reranked_best.as_ref().map(|&(rank, surrogate, model, _)| {
                serde_json::json!({
                    "rank_coupling": rank,
                    "surrogate_total": surrogate,
                    "model_of_record_total": model,
                })
            }).unwrap_or(serde_json::Value::Null),
            "reranked_best_is_search_best": reranked_best
                .as_ref()
                .map(|&(rank, _, _, _)| rank == 0)
                .unwrap_or(false),
            "completeness_check": {
                "convention": "an excluded allele could lie on a chain whose surrogate is at most \
                    W + max_refund when the all-allele loss minima plus the finalist boundary \
                    minima leave room; an excluded allele's own transitions were never scored \
                    (their matrix entries do not exist) — the check uses the finalist matrices' \
                    minima as their floor",
                "violations": completeness_violations,
                "ok": completeness_violations.is_empty(),
            },
            "per_finalist_shown": shown,
            "probe_routes": probe_rows,
            "wall_seconds": finalists_seconds,
        })
    };
    let reranked_route: Option<Vec<[usize; 2]>> = reranked_best
        .as_ref()
        .map(|&(_, _, _, ref route)| route.clone());
    let reranked_surrogate = reranked_best
        .as_ref()
        .map(|&(_, surrogate, _, _)| surrogate)
        .unwrap_or(dp_haploid.best_score);
    let reranked_model = reranked_best
        .as_ref()
        .map(|&(_, _, model, _)| model)
        .unwrap_or(f64::INFINITY);

    // Ploidy selection: the better track wins (identical tie handling to
    // the terminal-best rule; a bit-exact tie prefers the diploid track,
    // the previous behavior).
    // The ploidy comparison under the omission charge (owner ruling,
    // 2026-09-24): the haploid track's objective carries the omission costs
    // while the diploid track's internal DP stays on its frozen objective,
    // so the raw totals are cross-model. The comparison is made
    // model-consistent by adding the diploid best chain's omission addend
    // POST-HOC (exact for that chain: the per-locus in-window observed
    // features its two slots' profiles do not explain). The diploid DP's
    // own selection is untouched (documented approximation: which diploid
    // chain is best is decided under the frozen objective). Under the
    // finalist re-rank (design note (d)) the haploid side is the RE-RANKED
    // best chain's surrogate total — the same per-window objective family
    // the diploid side's post-hoc addend lands in; the model-of-record
    // ploidy numbers are reported alongside.
    let diploid_omission_addend: f64 = (0..locus_count)
        .map(|locus| -> io::Result<f64> {
            let classes = &locus_classes[locus];
            let membership = &classes.membership;
            let mut union = Profile::new();
            for slot in 0..2 {
                let allele = dp_diploid.route[locus][slot];
                let profile = &classes.profiles[membership[allele]];
                union = merge_profiles(&[&union, profile])?;
            }
            Ok(omission_charge_addend(&union, &window_obs[locus], sample_backgrounds))
        })
        .collect::<io::Result<Vec<f64>>>()?
        .into_iter()
        .sum();
    let diploid_effective = dp_diploid.best_score + diploid_omission_addend;
    let selected_ploidy =
        if (reranked_surrogate.total_cmp(&diploid_effective)).is_lt() {
            "haploid"
        } else {
            "diploid"
        };
    let dp = if selected_ploidy == "haploid" {
        dp_haploid.clone()
    } else {
        dp_diploid.clone()
    };
    eprintln!(
        "[phasing] DP done: best {:.2}, ploidy {selected_ploidy}, {:.2}s",
        dp.best_score,
        dp.wall_seconds
    );

    // ---------------------------------------- P5: authoritative rescore.
    let rescore_started = Instant::now();
    let route = match (&reranked_route, selected_ploidy) {
        (Some(route), "haploid") => route.clone(),
        _ => dp.route.clone(),
    };
    // The boundary cost set matching the selected route's allele universe:
    // the re-ranked haploid chain's alleles are finalist states (only the
    // finalist matrices define their transition entries); the diploid track's
    // and the DP-best chain's alleles are the old retained sets (the old
    // matrices). For allele pairs present in BOTH universes the finalist
    // matrix's entry can differ from the old one (a composition's realizing
    // pairs grew with the finalist universe) — the surrogate totals pair each
    // chain with its own universe's matrices, consistently.
    let selected_costs: &[BoundaryTransitionCosts] =
        if selected_ploidy == "haploid" && reranked_route.is_some() {
            &finalist_transition_costs
        } else {
            &transition_costs
        };
    // A selected allele as public JSON, or the empty-slot ploidy marker.
    let allele_or_empty_json = |locus: usize, slot: usize| -> serde_json::Value {
        if route[locus][slot] == EMPTY_SLOT2 {
            serde_json::json!("empty-slot2")
        } else {
            segment_list_json(&combined[locus][route[locus][slot]])
        }
    };
    // Junction census of the selected route (public structure of the
    // selected chains' transitions, with the crossing-read evidence of the
    // NOVEL junctions: spanning reads, the restricted charge, and the
    // pooled charge the OLD model would have given — the per-junction
    // before/after of the fix).
    let mut census_cooc_gap_zero = 0u64;
    let mut census_cooc_gap_positive = 0u64;
    let mut census_novel = 0u64;
    let mut census_novel_spanning = 0u64;
    let mut junction_rows: Vec<serde_json::Value> = Vec::new();
    for boundary in 0..locus_count - 1 {
        let mut kinds: [serde_json::Value; 2] = [serde_json::Value::Null, serde_json::Value::Null];
        for copy in 0..2 {
            // The haploid track's empty second slot: a ploidy statement —
            // no allele, no junction, no charge.
            if route[boundary][copy] == EMPTY_SLOT2 {
                kinds[copy] = serde_json::json!("empty-slot2");
                continue;
            }
            let left = &combined[boundary][route[boundary][copy]];
            let right = &combined[boundary + 1][route[boundary + 1][copy]];
            let kind = cooccurrence_gap(left, right);
            match kind {
                Some(gap) if gap >= 0 => {
                    if gap == 0 {
                        census_cooc_gap_zero += 1;
                        kinds[copy] = serde_json::json!("cooccurring-gap-zero");
                    } else {
                        census_cooc_gap_positive += 1;
                        kinds[copy] = serde_json::json!("cooccurring-gap");
                    }
                }
                _ => {
                    census_novel += 1;
                    // Crossing-read evidence for the selected novel junction.
                    let crossing = span.crossing_reads(
                        left.segments.last().expect("nonempty traversal"),
                        right.segments.first().expect("nonempty traversal"),
                        sources,
                        path_of_source,
                    )?;
                    if crossing.distinct_reads > 0 {
                        census_novel_spanning += 1;
                    }
                    let costs = &selected_costs[boundary];
                    let row =
                        costs.left_index[route[boundary][copy]] as usize * costs.right_count as usize;
                    let column = costs.right_index[route[boundary + 1][copy]] as usize;
                    let restricted = costs.cost[row + column];
                    let (left_tail, _) =
                        allele_endpoints(sources, flank_memo, &left.segments, READ_LENGTH - 1)?;
                    let (right_head, _) =
                        allele_endpoints(sources, flank_memo, &right.segments, READ_LENGTH - 1)?;
                    let (profile, _) = genome::profile_event_seam(
                        panel,
                        &left_tail,
                        &right_head,
                        READ_LENGTH,
                        MAX_FEATURES,
                    )?;
                    let pooled = profile_loss_boundary(
                        &profile,
                        &partition_obs[boundary],
                        &partition_obs[boundary + 1],
                        model,
                    )?;
                    kinds[copy] = serde_json::json!({
                        "kind": "novel",
                        "spanning_reads": crossing.distinct_reads,
                        "events": crossing.events,
                        "restricted_charge": restricted,
                        "pooled_charge_comparison": pooled,
                    });
                }
            }
        }
        junction_rows.push(serde_json::json!({
            "boundary": boundary,
            "slot0": [
                allele_or_empty_json(boundary, 0),
                allele_or_empty_json(boundary + 1, 0),
            ],
            "slot1": [
                allele_or_empty_json(boundary, 1),
                allele_or_empty_json(boundary + 1, 1),
            ],
            "slot0_kind": kinds[0],
            "slot1_kind": kinds[1],
        }));
    }
    // M1 oracle rescore: per-locus local pair terms over oracle (IR)
    // profiles plus the phasing boundary charges of the two slot chains.
    // Parametrized by the evaluated chain and its boundary cost set: the
    // DP-best chain's anchor rescore uses the tracks' own matrices; the
    // re-ranked selected chain's alleles may be finalist-universe additions,
    // covered only by the finalist matrices (the same per-window form, the
    // same accumulation order — the cost set changes only which entries are
    // defined).
    let mut surrogate_route_m1 = |route_rows: &[[usize; 2]],
                              costs_set: &[BoundaryTransitionCosts]|
     -> io::Result<(f64, f64)> {
    let mut oracle_memo: HashMap<String, Profile> = HashMap::new();
    let mut m1_local_total = 0.0f64;
    for locus in 0..locus_count {
        let first =
            oracle_allele_profile(panel, sources, &combined[locus][route_rows[locus][0]], &mut oracle_memo)?;
        // The haploid track's empty second slot contributes an empty
        // profile: the merged-pair loss degenerates to the single allele's
        // honest haploid charge.
        let second = if route_rows[locus][1] == EMPTY_SLOT2 {
            Profile::new()
        } else {
            oracle_allele_profile(panel, sources, &combined[locus][route_rows[locus][1]], &mut oracle_memo)?
        };
        // Extend the backgrounds for any oracle-profile feature beyond the
        // geometric candidate universe before charging (memoized re-hits).
        {
            let mut missing: Vec<FeatureKey> = Vec::new();
            for profile in [&first, &second] {
                for key in profile.keys() {
                    if !backgrounds.contains(key) {
                        missing.push(key.clone());
                    }
                }
            }
            if !missing.is_empty() {
                backgrounds.scan_extend(panel, missing, *k, model)?;
            }
        }
        // The haploid track's empty second slot is charged sample-side
        // (the same convention its DP optimized — the rescore and the DP
        // must be the same model); a diploid pair keeps the panel
        // convention.
        if route_rows[locus][1] == EMPTY_SLOT2 {
            // The omission charge (Fix 2): the rescore and the DP must be
            // the same model — the selected row keeps its current charge
            // against its owner-resolved observed side (the record-once map
            // over the allele's owner SET: the established Fix-1 form — a
            // single-owner allele's map is the precomputed partition map
            // bit-identically) and adds the per-window omission cost of the
            // window's observed DNA it does not spell (the DP's surrogate
            // layer keeps the per-window form — supervisor ruling (A)).
            let allele = &combined[locus][route_rows[locus][0]];
            let owner_set: BTreeSet<u32> = allele
                .segments
                .iter()
                .map(|segment| segment.partition as u32)
                .collect();
            let multi_owner_map = (owner_set.len() != 1).then(|| {
                site.site_map(
                    owner_set
                        .iter()
                        .map(|&owner| {
                            crate::owner_universe_partition(owner, component_locus_to_partition)
                        }),
                )
            });
            let owner_obs: &HashMap<FeatureKey, f64> = match &multi_owner_map {
                Some(map) => map,
                None => owner_routed_obs(
                    routed_equal,
                    *owner_set.iter().next().expect("nonempty owner set"),
                    component_locus_to_partition,
                ),
            };
            m1_local_total += crate::merged_single_loss_sample_with_omission(
                &first,
                owner_obs,
                &window_obs[locus],
                sample_backgrounds,
                model,
            )?;
        } else {
            let owner_first = traversal_owner(&combined[locus][route_rows[locus][0]]);
            let owner_second = traversal_owner(&combined[locus][route_rows[locus][1]]);
            if owner_first == owner_second {
                m1_local_total += merged_pair_loss_multiplicity(
                    &first,
                    &second,
                    owner_routed_obs(routed_equal, owner_first, component_locus_to_partition),
                    model,
                    backgrounds,
                )?;
            } else {
                // Fix 1 (record-once mixed-owner charging): the merged pair
                // charge's observed side attributes each record ONCE.
                let obs_site = site.site_map([
                    crate::owner_universe_partition(owner_first, component_locus_to_partition),
                    crate::owner_universe_partition(owner_second, component_locus_to_partition),
                ]);
                m1_local_total += merged_pair_loss_multiplicity(
                    &first,
                    &second,
                    &obs_site,
                    model,
                    backgrounds,
                )?;
                // Fix 2 (the omission charge, rescore scope): the diploid
                // pair's union profile pays the in-window observed features
                // neither slot explains.
                let pair_union = merge_profiles(&[&first, &second])?;
                m1_local_total += omission_charge_addend(
                    &pair_union,
                    &window_obs[locus],
                    sample_backgrounds,
                );
            }
        }
    }
    let haploid_selected = (0..locus_count).all(|locus| route_rows[locus][1] == EMPTY_SLOT2);
    let mut m1_oracle = m1_local_total;
    for boundary in 0..locus_count - 1 {
        let costs = &costs_set[boundary];
        let cost_matrix = if haploid_selected {
            &costs.cost_haploid
        } else {
            &costs.cost
        };
        for copy in 0..2 {
            // The empty second slot pays no boundary charge.
            if route_rows[boundary][copy] == EMPTY_SLOT2 {
                continue;
            }
            let row = costs.left_index[route_rows[boundary][copy]] as usize * costs.right_count as usize;
            let column = costs.right_index[route_rows[boundary + 1][copy]] as usize;
            m1_oracle += cost_matrix[row + column];
        }
    }
    Ok((m1_local_total, m1_oracle))
    };
    // The DP-best chain's anchor rescore (the tracks' own matrices): the
    // internal-vs-external self-check's external side — preserved bit-for-bit
    // from the predecessor's runs.
    let (_m1_local_total_dp_best, m1_oracle_dp_best) =
        surrogate_route_m1(&dp.route, &transition_costs)?;
    // The SELECTED chain's rescore (the ladder's selected_m1; the guard's
    // surrogate side) over the route's matching boundary cost set.
    let (m1_local_total, m1_oracle) = surrogate_route_m1(&route, selected_costs)?;
    // Pooled external rescore of the slot chains (the haploid track's
    // empty second slot spells nothing — external_rescore skips sequences
    // shorter than a read, so the pooled score is the honest single-
    // molecule one).
    let mut sequences = Vec::with_capacity(2);
    for copy in 0..2 {
        let mut sequence = Vec::new();
        for locus in 0..locus_count {
            if route[locus][copy] == EMPTY_SLOT2 {
                continue;
            }
            sequence.extend(allele_sequence(
                sources,
                &combined[locus][route[locus][copy]].segments,
            )?);
        }
        sequences.push(sequence);
    }
    let directory = rescore_directory.join("phasing-selected-pooled");
    let (pooled_loss, pooled_cost, pooled_runs) = genome::external_rescore(
        panel,
        [&sequences[0], &sequences[1]],
        sample_counts,
        model,
        &directory,
    )?;
    // References ladder (same-run scorer).
    let mut reference_scores = Vec::new();
    // The instance-probe diagnostic's spelled-map scratch (overwritten per
    // reference by the scorer's out param).
    let mut reference_spelled: HashMap<FeatureKey, Vec<crate::SpelledSpan>> = HashMap::new();
    for (index, (label, path_a, path_b)) in reference_specs.iter().enumerate() {
        let route_a: routes::Route = read_json(path_a)?;
        let route_b: routes::Route = read_json(path_b)?;
        let (loss, reference_census) = score_reference_m1(
            panel,
            sources,
            territory,
            [&route_a, &route_b],
            routed_equal,
            site,
            window_obs,
            window_instances,
            component_locus_to_partition,
            locus_count,
            locus_offset,
            model,
            Some((span, path_of_source, &locus_partitions)),
            Some(&mut **backgrounds),
            // A reference with an empty second route is a HAPLOID statement
            // and charges sample-side (the same convention as the haploid
            // track); diploid reference pairs keep the panel convention.
            Some(sample_backgrounds),
            // Group-shared material charged once (supervisor decision
            // 2026-09-23): the ladder's `m1_oracle_rescore` uses the
            // once-per-material convention; the legacy per-window number
            // is reported alongside during the transition.
            true,
            None,
            None,
            if instance_probe_features.is_empty() {
                None
            } else {
                Some(&mut reference_spelled)
            },
        )?;
        let instance_probe = (!instance_probe_features.is_empty()).then(|| {
            crate::instance_probe_report(
                window_obs,
                window_instances.unwrap_or(&crate::InstanceStructure::default()),
                &reference_spelled,
                sample_backgrounds,
                &instance_probe_features,
            )
        });
        let (loss_group_shared_legacy, _) = score_reference_m1(
            panel,
            sources,
            territory,
            [&route_a, &route_b],
            routed_equal,
            site,
            window_obs,
            // The legacy variant keeps the feature-level form (None).
            None,
            component_locus_to_partition,
            locus_count,
            locus_offset,
            model,
            Some((span, path_of_source, &locus_partitions)),
            Some(&mut **backgrounds),
            Some(sample_backgrounds),
            false,
            None,
            None,
            None,
        )?;
        let mut reference_sequences = Vec::with_capacity(2);
        for route_ref in [&route_a, &route_b] {
            let mut sequence = Vec::new();
            for segment in &route_ref.segments {
                if segment.start >= segment.end {
                    continue;
                }
                let mut part = sources.fetch(segment.source, segment.start, segment.end)?;
                if segment.reverse {
                    part = impg::graph::reverse_complement(&part);
                }
                sequence.extend(part);
            }
            reference_sequences.push(sequence);
        }
        let directory = rescore_directory.join(format!("phasing-reference-{index}-{label}"));
        let (pooled, pooled_reference_cost, pooled_reference_runs) = genome::external_rescore(
            panel,
            [&reference_sequences[0], &reference_sequences[1]],
            sample_counts,
            model,
            &directory,
        )?;
        reference_scores.push(serde_json::json!({
            "label": label,
            "m1_oracle_rescore": loss,
            "m1_oracle_rescore_group_shared_legacy": loss_group_shared_legacy,
            "m1_junction_census": reference_census,
            "instance_probe": instance_probe,
            "pooled_external_rescore": pooled,
            "pooled_mem_queries": pooled_reference_cost.mem_queries,
            "pooled_initial_runs": pooled_reference_runs,
        }));
        rss_probe(rss, &format!("phasing_reference_{label}"))?;
    }
    let rescore_seconds = rescore_started.elapsed().as_secs_f64();
    let reference_by_label = |label: &str| -> Option<f64> {
        reference_scores
            .iter()
            .find(|row| row["label"].as_str() == Some(label))
            .and_then(|row| row["m1_oracle_rescore"].as_f64())
    };
    let truth_m1 = reference_by_label("truth");
    let native2_m1 = reference_by_label("native2");
    // The HAPLOID truth (owner decision 2026-09-23): the single mosaic
    // route with an empty second slot — the generative fact on this
    // sample. Present when the run supplies a `haploid_truth` reference
    // pair (LABEL:truth-mosaic.json:empty-route.json).
    let haploid_truth_m1 = reference_by_label("haploid_truth");

    // The selected chain's CHAIN-LEVEL M1 (Ruling 2's model of record; the
    // M1 comparison's selected side): the selected route evaluated by the
    // SAME evaluator as the reference ladder — `score_reference_m1` with the
    // once-per-material piece filter and the chain-level omission addend.
    // The inline `m1_oracle` rescore above stays the DP's per-window
    // SURROGATE model (supervisor ruling (A): the DP keeps the per-window
    // objective with admissible bounds; its delta against this value is the
    // documented once-per-material refund, reported below). The finalist
    // re-rank evaluates every finalist with THIS SAME call shape
    // (`haploid_chain_level_m1`), so the selected chain's re-rank total is
    // bit-identical to this value by construction.
    let selected_chain: Vec<&genome::SpanningTraversal> = (0..locus_count)
        .map(|locus| &combined[locus][route[locus][0]])
        .collect();
    let selected_chain_level_m1 = haploid_chain_level_m1(
        panel,
        sources,
        territory,
        routed_equal,
        site,
        window_obs,
        window_instances,
        component_locus_to_partition,
        locus_count,
        locus_offset,
        model,
        span,
        path_of_source,
        &locus_partitions,
        &mut **backgrounds,
        sample_backgrounds,
        &selected_chain,
        None,
        Some(&mut shared_spelled_memo),
        None,
    )?;
    // The diploid best chain's chain-level M1 (supplementary, read-only): the
    // model-of-record ploidy comparison's diploid side — the diploid DP's own
    // route under the same scorer (the diploid branch; the track's selection
    // is untouched).
    let mut diploid_segments_a: Vec<routes::Segment> = Vec::new();
    let mut diploid_segments_b: Vec<routes::Segment> = Vec::new();
    for locus in 0..locus_count {
        for (slot, out) in [(0usize, &mut diploid_segments_a), (1usize, &mut diploid_segments_b)] {
            if dp_diploid.route[locus][slot] == EMPTY_SLOT2 {
                continue;
            }
            for segment in &combined[locus][dp_diploid.route[locus][slot]].segments {
                if segment.start >= segment.end {
                    continue;
                }
                out.push(routes::Segment {
                    source: segment.source,
                    start: segment.start,
                    end: segment.end,
                    reverse: segment.reverse,
                });
            }
        }
    }
    let diploid_best_route_a = routes::Route { segments: diploid_segments_a };
    let diploid_best_route_b = routes::Route { segments: diploid_segments_b };
    let (diploid_best_chain_level_m1, _diploid_best_chain_level_census) = score_reference_m1(
        panel,
        sources,
        territory,
        [&diploid_best_route_a, &diploid_best_route_b],
        routed_equal,
        site,
        window_obs,
        // The diploid branch keeps its per-window form (no instance
        // structure read).
        None,
        component_locus_to_partition,
        locus_count,
        locus_offset,
        model,
        Some((span, path_of_source, &locus_partitions)),
        Some(&mut **backgrounds),
        Some(sample_backgrounds),
        true,
        None,
        None,
        None,
    )?;

    // ---------------------------------------- P6: the tract verdict (assessment).
    let merge_intervals = |intervals: &[(usize, u64, u64)]| -> Vec<(usize, u64, u64)> {
        // Coverage is a SET union: sort by coordinate and merge OVERLAPPING
        // same-source intervals (the former adjacency-only merge double-
        // counted overlapping spellings — measured: the stitched chains'
        // donor coverage reported 126% when the strict union is 73.8%).
        let mut sorted: Vec<(usize, u64, u64)> = intervals.to_vec();
        sorted.sort_unstable();
        let mut merged: Vec<(usize, u64, u64)> = Vec::new();
        for (source, start, end) in sorted {
            match merged.last_mut() {
                Some(last) if last.0 == source && start <= last.2 => {
                    last.2 = last.2.max(end);
                }
                _ => merged.push((source, start, end)),
            }
        }
        merged
    };
    let truth_tract_intervals: Vec<(usize, u64, u64)> = merge_intervals(
        &truth_route_a
            .segments
            .iter()
            .chain(&truth_route_b.segments)
            .filter(|segment| donor_sources.contains(&segment.source))
            .map(|segment| (segment.source, segment.start, segment.end))
            .collect::<Vec<_>>(),
    );
    let truth_tract_bases: u64 = truth_tract_intervals
        .iter()
        .map(|&(_, start, end)| end - start)
        .sum();
    let mut donor_recovery = Vec::new();
    let mut donor_interval_rows = Vec::new();
    for copy in 0..2 {
        let raw: Vec<(usize, u64, u64)> = (0..locus_count)
            .flat_map(|locus| {
                // The empty second slot spells nothing.
                if route[locus][copy] == EMPTY_SLOT2 {
                    return Vec::new();
                }
                combined[locus][route[locus][copy]]
                    .segments
                    .iter()
                    .filter(|segment| donor_sources.contains(&segment.source))
                    .map(|segment| (segment.source, segment.start, segment.end))
                    .collect()
            })
            .collect();
        let intervals = merge_intervals(&raw);
        let total_bases: u64 = intervals.iter().map(|&(_, start, end)| end - start).sum();
        // STRICT set-union overlap (the predecessor's documented intent): the
        // intersection of [start,end) with [t_start,t_end) is
        // min(end,t_end) - max(start,t_start); the former
        // max(end,t_start) - min(start,t_end) form counted each interval's
        // FULL length whenever it merely touches the tract (the measured
        // 6,181-base over-run: [99069,102922) + [207743,210071)).
        let overlap: u64 = intervals
            .iter()
            .map(|&(_, start, end)| {
                truth_tract_intervals
                    .iter()
                    .map(|&(_, t_start, t_end)| {
                        end.min(t_end).saturating_sub(start.max(t_start))
                    })
                    .sum::<u64>()
            })
            .sum();
        donor_interval_rows.push(serde_json::json!({
            "copy": copy,
            "donor_intervals": intervals,
        }));
        donor_recovery.push(serde_json::json!({
            "copy": copy,
            "donor_intervals": intervals,
            "donor_bases": total_bases,
            "truth_tract_overlap_bases": overlap,
            "pct_of_truth_tract": if truth_tract_bases > 0 {
                100.0 * overlap as f64 / truth_tract_bases as f64
            } else {
                0.0
            },
        }));
    }
    // Per-locus structure comparison at the tract loci (up to swap).
    let sources_of_row = |row: &genome::SpanningTraversal| -> BTreeSet<usize> {
        row.segments.iter().map(|segment| segment.source).collect()
    };
    let donor_alleles_for = |locus: usize| -> Vec<usize> {
        combined[locus]
            .iter()
            .enumerate()
            .filter(|(_, row)| {
                row.segments
                    .iter()
                    .any(|segment| donor_sources.contains(&segment.source))
            })
            .map(|(allele, _)| allele)
            .collect()
    };
    let truth_sources_of = |locus: usize, copy: usize| -> BTreeSet<usize> {
        truth_pieces[locus][copy]
            .iter()
            .map(|&(source, _, _, _, _)| source)
            .collect()
    };
    let mut tract_rows = Vec::new();
    let mut structure_matches = 0usize;
    // Assessment side: the best (donor-row, native-row) pair at each tract
    // locus — the truth-structure pair at row granularity (the truth's own
    // mid-locus junction partials are not domain rows; this is the closest
    // expressible form), with its margin membership.
    let donor_native_pair = |locus: usize| -> Option<serde_json::Value> {
        let mut donor_alleles: Vec<usize> = Vec::new();
        let mut native_alleles: Vec<usize> = Vec::new();
        for (allele, row) in ranges[locus].iter().enumerate() {
            let sources: BTreeSet<usize> =
                row.segments.iter().map(|segment| segment.source).collect();
            if sources.iter().any(|source| donor_sources.contains(source)) {
                donor_alleles.push(allele);
            }
            if sources.len() == 1 && sources.contains(&target_source) {
                native_alleles.push(allele);
            }
        }
        let mut best: Option<(f64, usize, usize)> = None;
        for &donor in &donor_alleles {
            for &native in &native_alleles {
                let first = locus_classes[locus].membership[donor];
                let second = locus_classes[locus].membership[native];
                let loss = sweeps[locus].table[class_pair_index(first, second)];
                if best.is_none_or(|(value, _, _)| loss < value) {
                    best = Some((loss, donor, native));
                }
            }
        }
        best.map(|(loss, donor, native)| {
            serde_json::json!({
                "best_loss": loss,
                "within_margin": loss <= margins[locus],
                "donor_allele": donor,
                "native_allele": native,
            })
        })
    };
    let mut structure_matches_up_to_swap = 0usize;
    for &locus in tract_loci.iter() {
        let selected_first: BTreeSet<usize> = sources_of_row(&combined[locus][route[locus][0]]);
        let selected_second: BTreeSet<usize> = if route[locus][1] == EMPTY_SLOT2 {
            BTreeSet::new()
        } else {
            sources_of_row(&combined[locus][route[locus][1]])
        };
        let truth_first = truth_sources_of(locus, 0);
        let truth_second = truth_sources_of(locus, 1);
        let direct = selected_first == truth_first && selected_second == truth_second;
        let swapped = selected_first == truth_second && selected_second == truth_first;
        if direct || swapped {
            structure_matches_up_to_swap += 1;
        }
        if direct {
            structure_matches += 1;
        }
        // Ruling (d) verification (dp-evaluator-consistency): are the
        // decompose-level winners — the truth pieces and the donor pieces —
        // present, retained, and competitive under the per-allele oracle
        // record-once losses?
        let losses = &haploid_loss_tables[locus];
        let finite: Vec<(usize, f64)> = losses
            .iter()
            .copied()
            .enumerate()
            .filter(|&(_, value)| value.is_finite())
            .collect();
        let mut ranked = finite.clone();
        ranked.sort_by(|left, right| left.1.total_cmp(&right.1));
        let best_single = ranked.first().map(|&(_, value)| value);
        let swing = margins[locus] - sweeps[locus].min_viable_loss;
        let threshold = best_single.map(|best| best + swing);
        let truth_pieces_json: Vec<serde_json::Value> = truth_pieces[locus][0]
            .iter()
            .map(|&(source, start, end, reverse, _)| {
                // The piece's exact allele (same source interval and
                // orientation) and its oracle retention state.
                let exact = combined[locus].iter().position(|row| {
                    row.segments.len() == 1
                        && row.segments[0].source == source
                        && row.segments[0].start == start
                        && row.segments[0].end == end
                        && row.segments[0].reverse == reverse
                });
                match exact {
                    Some(allele) => serde_json::json!({
                        "piece": [source, start, end, reverse],
                        "allele": allele,
                        "oracle_loss": losses[allele],
                        "retained": retained_haploid_sets[locus].contains(&allele),
                        "rank": ranked.iter().position(|&(id, _)| id == allele),
                    }),
                    None => serde_json::json!({
                        "piece": [source, start, end, reverse],
                        "allele": null,
                    }),
                }
            })
            .collect();
        let donor_best = donor_alleles_for(locus)
            .into_iter()
            .filter_map(|allele| {
                let loss = losses[allele];
                if loss.is_finite() {
                    Some((loss, allele))
                } else {
                    None
                }
            })
            .min_by(|left, right| left.0.total_cmp(&right.0));
        let donor_best_json = donor_best.map(|(loss, allele)| {
            serde_json::json!({
                "allele": allele,
                "oracle_loss": loss,
                "retained": retained_haploid_sets[locus].contains(&allele),
                "rank": ranked.iter().position(|&(id, _)| id == allele),
                "segments": combined[locus][allele]
                    .segments
                    .iter()
                    .map(|segment| {
                        [segment.source, segment.start as usize, segment.end as usize]
                    })
                    .collect::<Vec<_>>(),
            })
        });
        let oracle_winners = serde_json::json!({
            "best_single_oracle": best_single,
            "swing": swing,
            "threshold": threshold,
            "scored_alleles": finite.len(),
            "truth_piece_alleles": truth_pieces_json,
            "best_donor_allele": donor_best_json,
        });
        tract_rows.push(serde_json::json!({
            "locus": locus,
            "full_locus": locus_offset + locus,
            "axis_interval": [axis_slice[locus].start, axis_slice[locus].end],
            "selected_slot0_sources": selected_first.iter().copied().collect::<Vec<_>>(),
            "selected_slot1_sources": selected_second.iter().copied().collect::<Vec<_>>(),
            "truth_copy0_sources": truth_first.iter().copied().collect::<Vec<_>>(),
            "truth_copy1_sources": truth_second.iter().copied().collect::<Vec<_>>(),
            "structure_matches_truth": direct,
            "structure_matches_truth_up_to_swap": direct || swapped,
            "truth_pair_loss": truth_losses[locus],
            "selected_pair_loss": if route[locus][1] == EMPTY_SLOT2 {
                serde_json::json!({
                    "ploidy": "haploid",
                    "haploid_loss": haploid_loss_tables[locus][route[locus][0]],
                })
            } else {
                serde_json::json!(sweeps[locus].table[class_pair_index(
                    locus_classes[locus].membership[route[locus][0]],
                    locus_classes[locus].membership[route[locus][1]],
                )])
            },
            "best_donor_native_pair": donor_native_pair(locus),
            "oracle_record_once_winners": oracle_winners,
        }));
    }
    // Post-hoc physical diagnostics: same-source overlap violations inside
    // each selected molecule (the phasing DP enforces no dosage legality).
    let molecule_overlap_violations = |copy: usize| -> serde_json::Value {
        let mut by_source: BTreeMap<usize, Vec<(u64, u64, bool)>> = BTreeMap::new();
        for locus in 0..locus_count {
            if route[locus][copy] == EMPTY_SLOT2 {
                continue;
            }
            for segment in &combined[locus][route[locus][copy]].segments {
                if segment.start < segment.end {
                    by_source
                        .entry(segment.source)
                        .or_default()
                        .push((segment.start, segment.end, segment.reverse));
                }
            }
        }
        let mut pairs = 0u64;
        let mut bases = 0u64;
        // Ruling-1 rerun diagnostic (genome/stitching-omission-alignment):
        // same-source overlapping spans split by CLASS. An identical-interval
        // same-orientation repeat is the harness's group-window structure —
        // the same piece legitimately assigned at several windows of a group
        // (the truth reference's own route has such repeats; the model of
        // record charges each piece once) — NOT a molecule-legality
        // violation. A violation is an ORIENTATION CONFLICT (the same DNA
        // used in both directions) or an overlap between DISTINCT segments
        // (a true self-overlap). The old double-pick (5553 fwd+rev) was the
        // conflict class.
        let mut conflict_pairs = 0u64;
        let mut conflict_bases = 0u64;
        let mut repeat_pairs = 0u64;
        let mut repeat_bases = 0u64;
        let mut rows = Vec::new();
        for (source, spans) in &mut by_source {
            spans.sort_unstable();
            for left in 0..spans.len() {
                for right in left + 1..spans.len() {
                    if spans[right].0 >= spans[left].1 {
                        break;
                    }
                    let overlap = spans[left].1.min(spans[right].1) - spans[right].0.max(spans[left].0);
                    if overlap > 0 {
                        pairs += 1;
                        bases += overlap;
                        let identical_repeat = spans[left].0 == spans[right].0
                            && spans[left].1 == spans[right].1
                            && spans[left].2 == spans[right].2;
                        if identical_repeat {
                            repeat_pairs += 1;
                            repeat_bases += overlap;
                        } else {
                            conflict_pairs += 1;
                            conflict_bases += overlap;
                        }
                        if rows.len() < 64 {
                            rows.push(serde_json::json!({
                                "source": source,
                                "left": [spans[left].0, spans[left].1],
                                "right": [spans[right].0, spans[right].1],
                                "overlap": overlap,
                            }));
                        }
                    }
                }
            }
        }
        serde_json::json!({
            "copy": copy,
            "overlapping_pairs": pairs,
            "overlap_bases": bases,
            "conflicting_pairs": conflict_pairs,
            "conflict_bases": conflict_bases,
            "identical_repeat_pairs": repeat_pairs,
            "identical_repeat_bases": repeat_bases,
            "examples": rows,
        })
    };
    let route_json: Vec<serde_json::Value> = (0..locus_count)
        .map(|locus| {
            serde_json::json!([
                if route[locus][0] == EMPTY_SLOT2 {
                    serde_json::json!("empty-slot2")
                } else {
                    traversal_json(&combined[locus][route[locus][0]])
                },
                if route[locus][1] == EMPTY_SLOT2 {
                    serde_json::json!("empty-slot2")
                } else {
                    traversal_json(&combined[locus][route[locus][1]])
                },
            ])
        })
        .collect();
    let self_delta = m1_oracle_dp_best - dp.best_score;
    let dp_beats_incumbent = (dp.best_score.total_cmp(&incumbent_phasing)).is_le();

    Ok(serde_json::json!({
        "model": "local-first-spine-v1-correlation-phasing",
        "stage1_sweep": stage1_summary,
        "phasing": {
            "co_occurrence": {
                "panel_structure": "each panel chromosome is a single source lane end-to-end; \
                    co-occurrence between rows at adjacent loci is the public \
                    same-source order-compatible adjacency (the panel provides no \
                    cross-source co-occurrence to correlate)",
                "boundaries": cooccurrence_rows,
                "wall_seconds": cooc_seconds,
            },
            "option_a_structural_partials": {
                "firing_rule": "canonical per-source minimal-overlap continuity failure \
                    among margin-retained single-segment forward rows; partials cut at the \
                    partner row's endpoint port (the existing split machinery is not fired \
                    from boundary failures — see the report)",
                "boundaries": option_a_stats,
                "augmented_loci": augmented,
                "total_partials": total_partials,
                "wall_seconds": option_a_seconds + rebuild_seconds,
            },
            "transition_costs": {
                "boundaries": transition_costs.iter().map(|costs| costs.stats.clone())
                    .collect::<Vec<_>>(),
                "boundary_observation_floors": boundary_floors,
                "wall_seconds": costs_seconds,
            },
            "chain_dp": {
                "retained_margins_local": margins,
                "suffix_from_locus": suffix_from_locus,
                "suffix_haploid_from_locus": suffix_haploid_from_locus,
                "selected_ploidy": selected_ploidy,
                "incumbent_native_backbone_phasing": incumbent_phasing,
                "incumbent_native_backbone_haploid": incumbent_haploid,
                "incumbent_native_backbone_spine": incumbent_spine,
                "incumbent_delta": incumbent_delta,
                "state_counts": dp.state_counts,
                "candidate_counts": dp.candidate_counts,
                "transitions": dp.transitions,
                "suffix_pruned": dp.suffix_pruned,
                "best_score_internal": dp.best_score,
                "diploid_best_score_internal": dp_diploid.best_score,
                "haploid_best_score_internal": dp_haploid.best_score,
                "finalist_l_dp_best_internal": dp_finalist.as_ref().expect("the finalists report implies the finalist DP").best_score,
                "terminal_tie_count": dp.tie_count,
                "dp_beats_incumbent": dp_beats_incumbent,
                "wall_seconds": dp.wall_seconds,
            },
            "finalists": finalist_summary,
            "rescore": {
                "selected_route": route_json,
                "selected_ploidy": selected_ploidy,
                "m1_oracle_rescore": m1_oracle,
                "m1_oracle_local_terms": m1_local_total,
                "selected_pooled_external_rescore": pooled_loss,
                "selected_pooled_mem_queries": pooled_cost.mem_queries,
                "selected_pooled_initial_runs": pooled_runs,
                "internal_vs_external_self_check": {
                    // The DP-BEST chain's anchor (the tracks' own matrices):
                    // preserved bit-for-bit from the predecessor's runs.
                    "internal_geometric": dp.best_score,
                    "external_m1_oracle": m1_oracle_dp_best,
                    "delta": self_delta,
                    "delta_percent": if dp.best_score != 0.0 {
                        100.0 * self_delta / dp.best_score
                    } else {
                        0.0
                    },
                },
                "selected_chain_self_check": {
                    // The SELECTED chain (the re-ranked finalist when the
                    // haploid track won): its surrogate total from the
                    // finalist tables vs its external per-window rescore.
                    "internal_geometric": if selected_ploidy == "haploid" && reranked_route.is_some() {
                        serde_json::json!(reranked_surrogate)
                    } else {
                        serde_json::json!(dp.best_score)
                    },
                    "external_m1_oracle": m1_oracle,
                },
                "references": reference_scores,
                "ladder": {
                    "truth_m1": truth_m1,
                    "native2_m1": native2_m1,
                    "haploid_truth_m1": haploid_truth_m1,
                    "selected_m1": m1_oracle,
                    "selected_chain_level_m1": selected_chain_level_m1,
                    "haploid_truth_minus_selected_chain_level": haploid_truth_m1
                        .map(|truth| truth - selected_chain_level_m1),
                    "truth_minus_selected": truth_m1.map(|truth| truth - m1_oracle),
                    "selected_minus_native2": native2_m1.map(|native| m1_oracle - native),
                    "haploid_truth_minus_selected": haploid_truth_m1
                        .map(|truth| truth - m1_oracle),
                    "truth_minus_native2": match (truth_m1, native2_m1) {
                        (Some(truth), Some(native)) => Some(truth - native),
                        _ => None,
                    },
                },
                "surrogate_vs_chain_level": {
                    // Supervisor guard (ruling (A)): the DP's per-window
                    // surrogate vs the chain-level once-per-material model of
                    // record. The refund is the documented once-per-material
                    // delta; the flag reports whether the surrogate's
                    // selected chain is ranked BELOW a chain-level rival
                    // (the haploid truth reference) — the documented
                    // surrogate limitation when true.
                    "selected_surrogate_m1": m1_oracle,
                    "selected_chain_level_m1": selected_chain_level_m1,
                    "once_per_material_refund": m1_oracle - selected_chain_level_m1,
                    "surrogate_selects_chain_level_suboptimal": haploid_truth_m1
                        .map(|truth| selected_chain_level_m1 > truth),
                    // The k-best finalist re-rank (supervisor ruling,
                    // genome/finalist-reranking): the surrogate-vs-model
                    // disagreement is now RESOLVED by construction — the
                    // selected chain IS the re-ranked best. The flag below
                    // reports whether the surrogate's own best chain lost the
                    // re-rank (the guard's original question, now measured
                    // over the enumerated finalist set).
                    "rerank_selected_model_of_record_total": if selected_ploidy == "haploid" && reranked_route.is_some() {
                        serde_json::json!(reranked_model)
                    } else {
                        serde_json::Value::Null
                    },
                    "rerank_selected_equals_surrogate_best": reranked_best
                        .as_ref()
                        .map(|&(rank, _, _, _)| rank == 0)
                        .unwrap_or(false),
                    "rerank_surrogate_best_model_of_record_total": if surrogate_best_model.is_finite() {
                        serde_json::json!(surrogate_best_model)
                    } else {
                        serde_json::Value::Null
                    },
                },
                "ploidy_comparison_model_of_record": {
                    // Supplementary (read-only): the re-ranked haploid best's
                    // chain-level total vs the diploid best chain's
                    // chain-level total (the same scorer; the diploid track's
                    // own selection untouched).
                    "haploid_reranked_best_model": if selected_ploidy == "haploid" && reranked_route.is_some() {
                        serde_json::json!(reranked_model)
                    } else {
                        serde_json::Value::Null
                    },
                    "diploid_best_chain_level_m1": diploid_best_chain_level_m1,
                },
                "junction_census": {
                    "cooccurring_gap_zero": census_cooc_gap_zero,
                    "cooccurring_gap_positive": census_cooc_gap_positive,
                    "novel": census_novel,
                    "novel_with_spanning_reads": census_novel_spanning,
                    "per_boundary": junction_rows,
                },
                "wall_seconds": rescore_seconds,
            },
            "tract_verdict": {
                "tract_loci": tract_loci,
                "tract_loci_count": tract_loci.len(),
                "per_locus": tract_rows,
                "structure_matches_truth_count": structure_matches,
                "structure_matches_truth_up_to_swap_count": structure_matches_up_to_swap,
                "donor_recovery": donor_recovery,
                "donor_intervals_per_copy": donor_interval_rows,
                "truth_tract_intervals": truth_tract_intervals,
                "truth_tract_bases": truth_tract_bases,
                "molecule_overlap_violations": [
                    molecule_overlap_violations(0),
                    molecule_overlap_violations(1),
                ],
            },
            "flags": {
                "exact_dp": "plain min-plus over margin-retained ordered pair states; the only \
                    prune is the admissible suffix bound (observation floors) against the \
                    native-backbone incumbent under the phasing transition model; TWO ploidy \
                    tracks (owner decision 2026-09-23, ploidy-verified sample): diploid (two \
                    real slots) and haploid (second slot legitimately empty at every locus — a \
                    ploidy statement, constant along the chain; the empty slot predicts \
                    nothing, pays nothing; only scorable classes are haploid-admissible); the \
                    better track is selected and both are reported",
                "no_port_word_legality": true,
                "dosage_legality": "not enforced on phasing transitions; the selected chains' \
                    same-source overlap violations are reported in tract_verdict",
                "transition_semantics": "co-occurring adjacencies are charged the panel-attested \
                    gap-filled real junction; all other pairs are charged the novel-junction \
                    juxtaposition seam; both use the same L149 event-seam machinery and the \
                    same adjacent-share charge, so gap-0 co-occurrences are bit-identical to \
                    the port-legal seams",
                "retained_states": "the Stage-1 reported local seam-swing margin over the \
                    best port-viable pair, plus structural retention of the backbone pair; \
                    non-augmented loci keep their pre-augmentation viability statistics \
                    (conservative, mirroring the spine Stage-3 precedent)",
                "co_occurrence_caveat": "for multi-segment rows whose exit segment is shorter \
                    than the L149 flank, the co-occurring left flank is the molecule-true \
                    oriented tail (which may cross the row's interior junction) while the \
                    right flank stays the source's own continuation",
                "no_tuning_constants": true,
                "target_source": target_source,
                "donor_sources": donor_sources.iter().copied().collect::<Vec<_>>(),
            },
            "wall_seconds": {
                "co_occurrence": cooc_seconds,
                "option_a_and_rebuild": option_a_seconds + rebuild_seconds,
                "transition_costs": costs_seconds,
                "chain_dp": dp.wall_seconds,
                "rescore": rescore_seconds,
                "total": started.elapsed().as_secs_f64(),
            },
        },
        "rss_peak_bytes": rss.peak_bytes(),
        "total_seconds": started.elapsed().as_secs_f64(),
    }))
}

fn segment_list_json(traversal: &genome::SpanningTraversal) -> serde_json::Value {
    serde_json::Value::Array(
        traversal.segments.iter().map(segment_json).collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The DP path's evaluator identity (owner ruling,
    /// dp-evaluator-consistency, 2026-09-24): the haploid track's states
    /// charge each retained ALLELE by its OWN per-allele record-once loss
    /// (`haploid_allele_losses`), never the class-level geometric loss of
    /// its membership class — the two indexings differ whenever alleles
    /// share a class under the oracle evaluation, and the class indexing
    /// is exactly the arithmetic that could not represent the
    /// decompose-level wins. Also locks the scorable hygiene: a retained
    /// allele whose loss is INFINITY (non-scorable) produces no state.
    #[test]
    fn haploid_states_charge_per_allele_record_once_losses() {
        let locus_classes = LocusClassing {
            profiles: vec![Profile::new(), Profile::new()],
            membership: vec![1, 0, 1],
            class_owners: vec![7, 9],
            class_charges: vec![0.0, 0.0],
            substage: [0.0; 4],
            split_seam_queries: 0,
            ir_extractions: 0,
        };
        // Per-allele oracle record-once losses: allele 0 (class 1) 10.0,
        // allele 1 (class 0) 5.0, allele 2 (class 1) INFINITY
        // (non-scorable). The CLASS-level indexing would charge allele 0
        // and allele 2 both with the class-1 entry and lose the
        // per-allele distinction.
        let losses = vec![10.0f64, 5.0, f64::INFINITY];
        let retained: HashSet<usize> = [0usize, 1, 2].into();
        let states = ordered_haploid_states(&losses, &retained).unwrap();
        // Allele 2's INFINITY loss is not a state (hygiene (c)); the
        // surviving states carry their OWN losses, sorted by allele id.
        assert_eq!(
            states,
            vec![([0usize, EMPTY_SLOT2], 10.0f64), ([1usize, EMPTY_SLOT2], 5.0f64)]
        );
        // A loss map where class-level indexing would invert the order:
        // allele 0 (class 1) cheaper than allele 1 (class 0) per allele,
        // while the class-1 geometric entry is the more expensive one.
        let losses = vec![3.0f64, 9.0];
        let retained: HashSet<usize> = [0usize, 1].into();
        let states = ordered_haploid_states(&losses, &retained).unwrap();
        assert_eq!(
            states,
            vec![([0usize, EMPTY_SLOT2], 3.0f64), ([1usize, EMPTY_SLOT2], 9.0f64)]
        );
    }

    /// The k-best finalist enumerator (supervisor ruling,
    /// genome/finalist-reranking): over a synthetic haploid DP it yields
    /// EVERY chain exactly once, in EXACTLY nondecreasing total order, with
    /// the DP's own accumulation order (the rank-0 value is bit-identical to
    /// the layer score; max_seed_delta stays 0.0), including bit-exact ties.
    #[test]
    fn finalist_enumerator_yields_all_chains_in_exact_total_order() {
        let empty_slot = [0usize, EMPTY_SLOT2];
        let losses: [Vec<f64>; 3] = [vec![10.0, 12.0, 30.0], vec![5.0, 5.0], vec![1.0, 1.0]];
        // Dense haploid transition matrices (3x2 then 2x2), all finite.
        let t0 = vec![1.5, 0.5, 2.0, 0.25, 0.0, 3.0];
        let t1 = vec![0.5, 2.5, 1.0, 0.0];
        let costs = vec![
            BoundaryTransitionCosts {
                cost: t0.clone(),
                cost_haploid: t0.clone(),
                left_index: vec![0, 1, 2],
                right_index: vec![0, 1],
                right_count: 2,
                stats: serde_json::Value::Null,
            },
            BoundaryTransitionCosts {
                cost: t1.clone(),
                cost_haploid: t1.clone(),
                left_index: vec![0, 1],
                right_index: vec![0, 1],
                right_count: 2,
                stats: serde_json::Value::Null,
            },
        ];
        // The DP's own layers (the same accumulation order the forward pass
        // uses: prefix + transition + own loss).
        let mut layers: Vec<Vec<PhState>> = Vec::new();
        layers.push(
            losses[0]
                .iter()
                .enumerate()
                .map(|(allele, &loss)| PhState {
                    score: loss,
                    pair: [allele, EMPTY_SLOT2],
                    pred: u32::MAX,
                    loss,
                })
                .collect(),
        );
        for locus in 1..3 {
            let previous = layers.last().expect("layer");
            let matrix = &costs[locus - 1].cost_haploid;
            let right_count = costs[locus - 1].right_count as usize;
            let mut layer = Vec::new();
            for (allele, &loss) in losses[locus].iter().enumerate() {
                let mut best = f64::INFINITY;
                for state in previous {
                    let candidate = state.score
                        + matrix[state.pair[0] * right_count + allele]
                        + loss;
                    if candidate < best {
                        best = candidate;
                    }
                }
                layer.push(PhState {
                    score: best,
                    pair: [allele, EMPTY_SLOT2],
                    pred: u32::MAX,
                    loss,
                });
            }
            layers.push(layer);
        }
        // Brute force: every chain, accumulated in path order.
        let mut brute: Vec<(f64, Vec<[usize; 2]>)> = Vec::new();
        for a0 in 0..3usize {
            for a1 in 0..2usize {
                for a2 in 0..2usize {
                    let mut total = losses[0][a0];
                    total += t0[a0 * 2 + a1];
                    total += losses[1][a1];
                    total += t1[a1 * 2 + a2];
                    total += losses[2][a2];
                    brute.push((total, vec![[a0, EMPTY_SLOT2], [a1, EMPTY_SLOT2], [a2, EMPTY_SLOT2]]));
                }
            }
        }
        brute.sort_by(|left, right| {
            left.0
                .total_cmp(&right.0)
                .then_with(|| left.1.cmp(&right.1))
        });
        let mut enumerator = FinalistEnumerator::new(&layers, &costs);
        let mut got: Vec<(f64, Vec<[usize; 2]>)> = Vec::new();
        while let Some((value, route)) = enumerator.next_chain() {
            got.push((value, route));
        }
        assert_eq!(got.len(), brute.len(), "every chain exactly once");
        for window in got.windows(2) {
            assert!(
                window[0].0.total_cmp(&window[1].0).is_le(),
                "exactly nondecreasing order"
            );
        }
        let mut got_sorted = got.clone();
        got_sorted.sort_by(|left, right| {
            left.0
                .total_cmp(&right.0)
                .then_with(|| left.1.cmp(&right.1))
        });
        assert_eq!(got_sorted, brute, "the same chain multiset and values");
        // The rank-0 value is the DP's own terminal best, bit-exact.
        let terminal_best = layers[2]
            .iter()
            .map(|state| state.score)
            .fold(f64::INFINITY, f64::min);
        assert_eq!(got[0].0.to_bits(), terminal_best.to_bits());
        assert_eq!(enumerator.max_seed_delta, 0.0f64);
        let _ = empty_slot;
    }

    /// The enumerator SKIPS chains through a forbidden (non-finite)
    /// transition, and its ordering survives the gap.
    #[test]
    fn finalist_enumerator_skips_forbidden_transitions() {
        let losses: [Vec<f64>; 2] = [vec![10.0, 12.0], vec![5.0, 7.0]];
        // 2x2 matrix with one forbidden entry (allele 0 -> allele 1).
        let t0 = vec![1.0, f64::INFINITY, 2.0, 0.5];
        let costs = vec![BoundaryTransitionCosts {
            cost: t0.clone(),
            cost_haploid: t0.clone(),
            left_index: vec![0, 1],
            right_index: vec![0, 1],
            right_count: 2,
            stats: serde_json::Value::Null,
        }];
        let layers = vec![
            losses[0]
                .iter()
                .enumerate()
                .map(|(allele, &loss)| PhState {
                    score: loss,
                    pair: [allele, EMPTY_SLOT2],
                    pred: u32::MAX,
                    loss,
                })
                .collect::<Vec<_>>(),
            losses[1]
                .iter()
                .enumerate()
                .map(|(allele, &loss)| PhState {
                    score: match allele {
                        0 => 10.0 + 1.0 + 5.0,
                        1 => 10.0 + 2.0 + 7.0,
                        _ => unreachable!(),
                    },
                    pair: [allele, EMPTY_SLOT2],
                    pred: 0,
                    loss,
                })
                .collect::<Vec<_>>(),
        ];
        let mut enumerator = FinalistEnumerator::new(&layers, &costs);
        let mut got: Vec<(f64, Vec<[usize; 2]>)> = Vec::new();
        while let Some((value, route)) = enumerator.next_chain() {
            got.push((value, route));
        }
        // 3 reachable chains (0->0, 1->0, 1->1); 0->1 is forbidden.
        assert_eq!(got.len(), 3, "got: {got:?}");
        for window in got.windows(2) {
            assert!(window[0].0.total_cmp(&window[1].0).is_le());
        }
        assert_eq!(got[0].0.to_bits(), (10.0f64 + 1.0 + 5.0).to_bits());
    }

    /// The pairwise-coupling search objective's arithmetic (supervisor
    /// ruling on genome/finalist-reranking, option (b)): C = trans minus the
    /// adjacent-window exoneration channel, priced with the model of
    /// record's own per-feature omission terms; the exoneration grows with
    /// theneighbor's spelling and shrinks with the allele's own prediction.
    #[test]
    fn coupling_transition_costs_subtract_adjacent_exoneration() {
        let terms_left = vec![10.0f64, 20.0, 40.0];
        let terms_right = vec![1.0f64, 2.0];
        let exon = |self_list: Vec<u32>, self_sum: f64, left: Option<(Vec<u32>, f64)>, right: Option<(Vec<u32>, f64)>| HaploidAlleleExon {
            self_list,
            self_sum,
            left,
            right,
            self_instances: AlleleCoveredInstances::default(),
            left_instances: None,
            right_instances: None,
        };
        // Locus b: allele 0 spells right-window features {0} (sum 1.0);
        // allele 1 spells {0,1} (sum 3.0).
        let left_exons = vec![
            exon(vec![], 0.0, None, Some((vec![0], 1.0))),
            exon(vec![], 0.0, None, Some((vec![0, 1], 3.0))),
        ];
        // Locus b+1: allele 0 predicts right-window NOTHING (self empty);
        // allele 1 predicts {0} (self sum 1.0).
        let right_exons = vec![
            exon(vec![], 0.0, Some((vec![1, 2], 60.0)), None),
            exon(vec![0], 1.0, Some((vec![2], 40.0)), None),
        ];
        let trans = vec![100.0f64, 100.0, 100.0, 100.0];
        let costs = BoundaryTransitionCosts {
            cost: trans.clone(),
            cost_haploid: trans.clone(),
            left_index: vec![0, 1],
            right_index: vec![0, 1],
            right_count: 2,
            stats: serde_json::Value::Null,
        };
        let coupling = build_coupling_transition_costs(
            &[0usize, 1],
            &[0usize, 1],
            &costs,
            &left_exons,
            &right_exons,
            &terms_left,
            &terms_right,
            // The legacy feature-level form (the test's own arithmetic).
            false,
            None,
            &[],
            &[],
            &[],
            &[],
            &[],
            &[],
        );
        // Pair (a=0, b'=0): exon(b+1<-b) = sum_right(0) - overlap({0},{}) =
        // 1.0; exon(b<-b+1) = sum_left(b'=0) - overlap({1,2},{}) = 60.0.
        assert_eq!(coupling.cost_haploid[0].to_bits(), (100.0f64 - 61.0).to_bits());
        // Pair (a=1, b'=0): exon(b+1<-b) = 3.0 - overlap({0,1},{}) = 3.0;
        // exon(b<-b+1) = 60.0 - overlap({1,2},{}) = 60.0.
        assert_eq!(coupling.cost_haploid[2].to_bits(), (100.0f64 - 63.0).to_bits());
        // Pair (a=0, b'=1): exon(b+1<-b) = 1.0 - overlap({0},{0}) = 0.0 (the
        // right allele predicts what the left spells - no exonation);
        // exon(b<-b+1) = sum_left(b'=1)=40.0 - overlap({2},{}) = 40.0.
        assert_eq!(coupling.cost_haploid[1].to_bits(), (100.0f64 - 40.0).to_bits());
        // Pair (a=1, b'=1): exon(b+1<-b) = 3.0 - overlap({0,1},{0}) = 2.0;
        // exon(b<-b+1) = 40.0 - overlap({2},{}) = 40.0.
        assert_eq!(coupling.cost_haploid[3].to_bits(), (100.0f64 - 42.0).to_bits());
        // The index maps are carried unchanged (the dense layout is the
        // finalist matrices').
        assert_eq!(coupling.left_index, costs.left_index);
        assert_eq!(coupling.right_count, costs.right_count);
    }

    /// The INSTANCE-LEVEL coupling (the corrected Ruling 2's pairwise
    /// channel): the neighbor's spelling credits only the instances it
    /// actually covers; the credit = term(u_self) - term(u_self - delta)
    /// over the UNCOVERED masses. The repeat-at-distinct-loci case: a
    /// spelling at one instance gives NO credit for the other instance's
    /// mass (the feature-level form credited the FULL term).
    #[test]
    fn instance_coupling_credits_only_covered_instances() {
        let model = ScoreModel {
            read_length: 150,
            histogram: 100,
            denominator: 15_000.0,
            depth: 10.0,
            background: 0.1,
        };
        let mut cross: HashMap<FeatureKey, f64> = HashMap::new();
        cross.insert(vec![7u64], 0.0);
        let sample = SampleSideBackgrounds::new(&cross, 0.1);
        // Window b+1 observes feature f (key index 0) with TWO instances at
        // DISTINCT positions: shares 3.0 and 1.0 (full mass 4.0).
        let keys_right = vec![vec![7u64]];
        let full_right = vec![4.0f64];
        let shares_right = vec![vec![3.0f64, 1.0]];
        let exon = |self_list: Vec<u32>, self_sum: f64, left: Option<(Vec<u32>, f64)>, right: Option<(Vec<u32>, f64)>| HaploidAlleleExon {
            self_list,
            self_sum,
            left,
            right,
            self_instances: AlleleCoveredInstances::default(),
            left_instances: None,
            right_instances: None,
        };
        // LEFT allele a spells NOTHING at window b (its right list empty);
        // its instance cover at window b+1 covers instance 0 only (share 3).
        let right_cover_a = AlleleCoveredInstances {
            features: vec![0u32],
            masses: vec![3.0f64],
            solo_credits: vec![
                crate::omission_term_value(&keys_right[0], 4.0, &sample)
                    - crate::omission_term_value(&keys_right[0], 1.0, &sample),
            ],
            ids: vec![vec![0u32]],
        };
        let mut left_exon_a = exon(vec![], 0.0, None, None);
        left_exon_a.right_instances = Some(right_cover_a);
        let left_exons = vec![left_exon_a];
        // RIGHT allele a' predicts nothing (self empty, no covers).
        let right_exons = vec![exon(vec![], 0.0, None, None)];
        let trans = vec![100.0f64];
        let costs = BoundaryTransitionCosts {
            cost: trans.clone(),
            cost_haploid: trans.clone(),
            left_index: vec![0],
            right_index: vec![0],
            right_count: 1,
            stats: serde_json::Value::Null,
        };
        let terms_right = vec![crate::omission_term_value(&keys_right[0], 4.0, &sample)];
        let coupling = build_coupling_transition_costs(
            &[0usize],
            &[0usize],
            &costs,
            &left_exons,
            &right_exons,
            &[],
            &terms_right,
            true,
            Some(&sample),
            &[],
            &keys_right,
            &[],
            &full_right,
            &[],
            &shares_right,
        );
        // exon(b+1<-b) = term(4.0) - term(1.0) (the spelled instance's
        // share 3.0 is covered; the OTHER instance's 1.0 stays uncovered) —
        // STRICTLY LESS than the feature-level form's term(4.0). The cost
        // = trans - exon (compared through the same subtraction).
        let expected = crate::omission_term_value(&keys_right[0], 4.0, &sample)
            - crate::omission_term_value(&keys_right[0], 1.0, &sample);
        assert_eq!(
            coupling.cost_haploid[0].to_bits(),
            (100.0f64 - expected).to_bits()
        );
        assert!(expected < terms_right[0]);
        // The mirror direction: a''s LEFT cover at window b; window b
        // observes nothing (empty keys) — no backward credit.
        // (full_left/shares_left empty; left_exon_a.self_instances empty.)
    }
}
