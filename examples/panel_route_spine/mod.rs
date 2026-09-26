//! Local-first spine (owner redesign 2026-09): exhaustive local sweep, exact
//! chain DP, targeted splits, authoritative rescore.
//!
//! The beam machinery and its accumulated constants (beam widths, tie
//! rotations, seed ladders, prefilter top-Ns, admission margins) are RETIRED
//! from this path. The quarantined `run_routed_dp` in the parent example is
//! untouched and remains the oracle comparison. The spine contains NO tuning
//! constants: every cut is either structural (exactness guards) or a
//! statistic derived from the run's own loss/seam/residual distributions and
//! reported with the result.
//!
//! Stages:
//! 1. EXHAUSTIVE LOCAL SWEEP — every locus, ALL allele pairs (unordered,
//!    with replacement) over the locus's panel (parent) allele set, scored
//!    under the routed M1 model (geometric q + routed observed support).
//!    Full ranking, margins, per-locus winner table, assessment-side truth
//!    comparison. No prefilter, no top-N, no beam.
//! 2. EXACT CHAIN DP — min-plus over locally-plausible pair states. Retention
//!    by an admissible, data-derived margin (a reported statistic); boundary
//!    evidence from the routed seam layer; both homolog matchings; dosage
//!    (same-source overlap) legality encoded in the state's windowed
//!    history. No width, no bound machinery beyond the exactness-preserving
//!    admissible-margin prune.
//! 3. TARGETED SPLITS — a residual detector (best local pair versus what
//!    the locus's observations can explain) selects loci; the existing split
//!    machinery is invoked AT THOSE LOCI ONLY. Selection is the natural gap
//!    of the residual distribution — a reported statistic, not a tuned cut.
//! 4. FINAL RESCORE — the selected pair through the authoritative oracle-q
//!    M1 rescore plus the unchanged pooled external rescore; internal
//!    versus external self-check.

use super::*;
pub(crate) mod junction;
mod phasing;

use junction::{JunctionSpanIndex, RestrictedCharge};
use phasing::{run_correlation_phasing, PhasingShared};

/// Triangular index into a per-locus exhaustive class-pair loss table
/// (`class_count * (class_count + 1) / 2` entries, `first <= second`).
pub(super) fn class_pair_index(first: usize, second: usize) -> usize {
    let (lo, hi) = (first.min(second), first.max(second));
    hi * (hi + 1) / 2 + lo
}

/// Exact-content interning of windowed conflict histories (each history
/// INCLUDES its layer's pair; the merge key `(pair, history_id)` is therefore
/// exact — no fingerprint collisions can merge distinct physical states).
fn intern_history(
    content: &[[usize; 2]],
    history_store: &mut Vec<Vec<[usize; 2]>>,
    history_map: &mut HashMap<Vec<[usize; 2]>, u32>,
) -> u32 {
    if let Some(&id) = history_map.get(content) {
        return id;
    }
    let id = history_store.len() as u32;
    history_store.push(content.to_vec());
    history_map.insert(content.to_vec(), id);
    id
}

// ---------------------------------------------------------------------------
// Per-locus classing (parent domain; reused after targeted splits).
// ---------------------------------------------------------------------------

pub(super) struct LocusClassing {
    pub(super) profiles: Vec<Profile>,
    pub(super) membership: Vec<usize>,
    /// Per-class OWNING universe partition (the window-domain extension's
    /// owner-resolved charging; the owning partition enters the class key,
    /// so classes stay owner-homogeneous; every class of the pre-extension
    /// model is owned by the locus's axis partition).
    pub(super) class_owners: Vec<u32>,
    /// Per-CLASS restricted interior-junction charge (the additive novel-
    /// junction term folded into the exhaustive class-pair table; zero for
    /// every class whose alleles carry no novel-only interior junction).
    /// Classes are charge-homogeneous by construction (the class key
    /// includes the charge), so the class-pair table with the charges
    /// folded in stays exact at allele granularity.
    pub(super) class_charges: Vec<f64>,
    /// (parents, sweeps, partials+seams, classing) sub-stage seconds.
    pub(super) substage: [f64; 4],
    pub(super) split_seam_queries: u64,
    pub(super) ir_extractions: u64,
}

/// Per-locus geometric classing of the spine's allele set. Faithful
/// extraction of the quarantined `run_routed_dp` classing block, reused by
/// the spine for the parent domain and again after targeted splits.
///
/// Orientation convention: a reverse single-segment allele shares its forward
/// twin's profile (exact for the self view — feature keys are RC-canonical
/// and window containment mirrors one-to-one under reverse complement, so a
/// sequence and its reverse complement produce identical profiles; the
/// inverted-repeat part is inherited from the forward twin's RC-view
/// extraction, matching the oracle's per-window both-orientation matching up
/// to the Gate-1 tie-window residual causes). Boundary seams and
/// reconstruction keep true orientation throughout.
pub(super) fn class_locus_alleles(
    panel: &SyngIndex,
    sources: &routes::Sources,
    flank_memo: &FlankMemo,
    path_of_source: &[usize],
    locus_ranges: &[genome::SpanningTraversal],
    locus: usize,
    k: u64,
    // The within-read adjacency index (the junction-evidence restriction:
    // novel interior junctions of 2-segment alleles are paid by crossing
    // reads only; co-occurring interior junctions keep the pooled seam
    // charge bit-identically).
    span: &JunctionSpanIndex,
    // The locus's universe partition (the crossing reads' share target for
    // the locus's own axis rows; extension-added rows charge their OWNING
    // partition — resolved per traversal below).
    partition: u32,
    // The routed shares (owner-resolved charging) and the component locus ->
    // universe partition map.
    routed_equal: &crate::RoutedObs,
    component_locus_to_partition: &[u32],
    model: &ScoreModel,
) -> io::Result<LocusClassing> {
    let parents_started = Instant::now();
    let mut parent_intervals: BTreeSet<(usize, u64, u64)> = BTreeSet::new();
    for traversal in locus_ranges {
        for segment in &traversal.segments {
            if segment.start < segment.end {
                parent_intervals.insert((segment.source, segment.start, segment.end));
            }
        }
    }
    let parent_list: Vec<(usize, u64, u64)> = parent_intervals.into_iter().collect();
    let parent_data: HashMap<(usize, u64, u64), (Vec<(i32, u64)>, Vec<Vec<(i32, u64)>>)> =
        parent_list
            .par_iter()
            .map(|&(source, start, end)| {
                let anchors =
                    contained_path_anchors(panel, path_of_source[source], start, end, k)?;
                let records = geometric_ir_records(panel, sources, source, start, end, k)?;
                Ok(((source, start, end), (anchors, records)))
            })
            .collect::<io::Result<_>>()?;
    let ir_extractions = parent_data.len() as u64;
    let parents_seconds = parents_started.elapsed().as_secs_f64();

    let sweeps_started = Instant::now();
    let parent_profiles: HashMap<(usize, u64, u64), Profile> = parent_data
        .par_iter()
        .map(|(key, (anchors, records))| {
            let profile = geometric_sweep_profile(
                anchors,
                records,
                key.2 - key.1,
                k,
                &format!("spine-locus-{locus}-interval-{}-{}-{}", key.0, key.1, key.2),
            )?;
            Ok((*key, profile))
        })
        .collect::<io::Result<_>>()?;
    let sweeps_seconds = sweeps_started.elapsed().as_secs_f64();

    let partials_started = Instant::now();
    let mut partial_keys: BTreeSet<(usize, u64, u64)> = BTreeSet::new();
    let mut junctions: BTreeSet<(Vec<u8>, Vec<u8>)> = BTreeSet::new();
    // Per seam composition: the realizing (exit, entry) segment pairs and
    // whether the composition is POOLED (some realizing pair is a
    // co-occurring real panel junction — those keep the pooled charge
    // bit-identically; the composition is the seam machinery's own
    // junction granularity, so the observed-side rule pools with it).
    let mut junction_pairs: BTreeMap<(Vec<u8>, Vec<u8>), (Vec<(SourceRange, SourceRange)>, bool)> =
        BTreeMap::new();
    for traversal in locus_ranges {
        if traversal.segments.len() == 2 {
            for segment in &traversal.segments {
                if segment.start < segment.end {
                    partial_keys.insert((segment.source, segment.start, segment.end));
                }
            }
            let junction = (
                segment_flank(
                    sources,
                    flank_memo,
                    &traversal.segments[0],
                    false,
                    READ_LENGTH - 1,
                )?,
                segment_flank(
                    sources,
                    flank_memo,
                    &traversal.segments[1],
                    true,
                    READ_LENGTH - 1,
                )?,
            );
            junctions.insert(junction.clone());
            let entry = junction_pairs.entry(junction).or_default();
            entry
                .0
                .push((traversal.segments[0].clone(), traversal.segments[1].clone()));
            let gap = junction::segment_pair_gap(&traversal.segments[0], &traversal.segments[1]);
            if gap.is_some_and(|value| value >= 0) {
                entry.1 = true;
            }
        }
    }
    let partials: HashMap<(usize, u64, u64), std::sync::Arc<Profile>> = partial_keys
        .par_iter()
        .map(|&key| {
            let profile = if let Some(full) = parent_profiles.get(&key) {
                full.clone()
            } else {
                partial_profile_from_parents(&parent_data, key, k)?
            };
            Ok((key, std::sync::Arc::new(profile)))
        })
        .collect::<io::Result<_>>()?;
    let seam_profiles_by_junction: HashMap<(Vec<u8>, Vec<u8>), std::sync::Arc<Profile>> =
        junctions
            .par_iter()
            .map(|junction| {
                let (seam, _) = genome::profile_event_seam(
                    panel,
                    &junction.0,
                    &junction.1,
                    READ_LENGTH,
                    MAX_FEATURES,
                )?;
                Ok((junction.clone(), std::sync::Arc::new(seam)))
            })
            .collect::<io::Result<_>>()?;
    let split_seam_queries = seam_profiles_by_junction.len() as u64;
    // Restricted charges for the NOVEL-ONLY interior junction compositions
    // (no realizing pair is a co-occurring real panel junction): the seam
    // profile is paid by the reads that actually cross the junction —
    // unioned over the composition's realizing segment pairs, against the
    // crossing reads' equal-share routed mass at this locus's partition.
    // Pooled compositions keep the merged seam profile (bit-identical to
    // the previous model for co-occurring interior junctions).
    let junction_outcomes: HashMap<(Vec<u8>, Vec<u8>), RestrictedCharge> = junction_pairs
        .par_iter()
        .filter(|(_, (pairs, pooled))| !*pooled && !pairs.is_empty())
        .map(|(junction, (pairs, _))| {
            let profile = seam_profiles_by_junction
                .get(junction)
                .ok_or_else(|| invalid("split junction seam missing"))?;
            // The junction's realizing pairs share one owning partition
            // (mixed-owner split admission is excluded); the crossing
            // reads' share target is that owner.
            let owner = pairs[0].0.partition as u32;
            let owner_partition = crate::owner_universe_partition(owner, component_locus_to_partition);
            let outcome = span.restricted_charge(
                profile,
                pairs,
                sources,
                path_of_source,
                &[owner_partition],
                model,
            )?;
            Ok((junction.clone(), outcome))
        })
        .collect::<io::Result<_>>()?;
    let novel_junction_count = junction_outcomes.len();
    let novel_junction_spanning = junction_outcomes
        .values()
        .filter(|outcome| outcome.spanning_reads > 0)
        .count() as u64;
    let partials_seconds = partials_started.elapsed().as_secs_f64();

    let classing_started = Instant::now();
    let mut class_profiles: Vec<Profile> = Vec::new();
    let mut class_charges: Vec<f64> = Vec::new();
    let mut class_owners: Vec<u32> = Vec::new();
    let mut class_buckets: HashMap<(u64, u64, u32), Vec<usize>> = HashMap::new();
    let mut membership = vec![0usize; locus_ranges.len()];
    for (chunk_base, chunk) in locus_ranges.chunks(2048).enumerate() {
        let built: Vec<(Profile, u64, f64, u32)> = chunk
            .par_iter()
            .map(|traversal| {
                // The candidate's OWNING universe partition: every segment
                // of a well-formed traversal shares one owner (single-
                // segment domain rows and the anchor-owned split candidates
                // alike).
                let owner = {
                    let first = traversal.segments[0].partition as u32;
                    if !traversal
                        .segments
                        .iter()
                        .all(|segment| segment.partition as u32 == first)
                    {
                        return Err(invalid(
                            "spine allele mixes candidate rows of distinct owning partitions",
                        ));
                    }
                    first
                };
                let mut charge = 0.0f64;
                let profile = if traversal.segments.len() == 1 {
                    let segment = &traversal.segments[0];
                    if segment.start == segment.end {
                        Profile::new()
                    } else {
                        parent_profiles
                            .get(&(segment.source, segment.start, segment.end))
                            .cloned()
                            .ok_or_else(|| invalid("spine allele interval missing parent data"))?
                    }
                } else {
                    ensure(
                        traversal.segments.len() == 2,
                        "spine allele has unsupported segment arity"
                    )?;
                    let left_segment = &traversal.segments[0];
                    let right_segment = &traversal.segments[1];
                    let mut parts: Vec<&Profile> = Vec::with_capacity(3);
                    if left_segment.start < left_segment.end {
                        let key = (left_segment.source, left_segment.start, left_segment.end);
                        parts.push(
                            partials
                                .get(&key)
                                .ok_or_else(|| invalid("split partial missing"))?,
                        );
                    }
                    if right_segment.start < right_segment.end {
                        let key = (right_segment.source, right_segment.start, right_segment.end);
                        parts.push(
                            partials
                                .get(&key)
                                .ok_or_else(|| invalid("split partial missing"))?,
                        );
                    }
                    let junction = (
                        segment_flank(
                            sources,
                            flank_memo,
                            left_segment,
                            false,
                            READ_LENGTH - 1,
                        )?,
                        segment_flank(
                            sources,
                            flank_memo,
                            right_segment,
                            true,
                            READ_LENGTH - 1,
                        )?,
                    );
                    match junction_outcomes.get(&junction) {
                        // Novel-only interior junction: the seam profile is
                        // NOT part of the allele's predicted profile; its
                        // restricted charge is the allele's additive term.
                        Some(outcome) => {
                            charge = outcome.charge;
                        }
                        // Pooled (panel-attested) interior junction: the
                        // seam profile stays merged, charged against the
                        // partition's pooled shares exactly as before.
                        None => {
                            parts.push(
                                seam_profiles_by_junction
                                    .get(&junction)
                                    .ok_or_else(|| invalid("split junction seam missing"))?,
                            );
                        }
                    }
                    merge_profiles(&parts)?
                };
                let signature = profile_signature_hash(&profile);
                Ok((profile, signature, charge, owner))
            })
            .collect::<io::Result<_>>()?;
        for (offset, (profile, signature, charge, owner)) in built.into_iter().enumerate() {
            // The class key includes the charge and the OWNING partition, so
            // classes stay charge- and owner-homogeneous and the class-pair
            // table with charges folded in remains exact at allele
            // granularity.
            let key = (signature, charge.to_bits(), owner);
            let class = if let Some(class) = class_buckets
                .get(&key)
                .and_then(|bucket| {
                    bucket
                        .iter()
                        .find(|&&class| {
                            class_profiles[class] == profile
                                && class_charges[class].to_bits() == charge.to_bits()
                                && class_owners[class] == owner
                        })
                        .copied()
                })
            {
                class
            } else {
                let class = class_profiles.len();
                class_profiles.push(profile);
                class_charges.push(charge);
                class_owners.push(owner);
                class_buckets.entry(key).or_default().push(class);
                class
            };
            membership[chunk_base * 2048 + offset] = class;
        }
    }
    let classing_seconds = classing_started.elapsed().as_secs_f64();
    if novel_junction_count > 0 {
        eprintln!(
            "[junction] locus {locus}: {} novel-only interior junction compositions, \
             {} with crossing reads",
            novel_junction_count, novel_junction_spanning
        );
    }
    Ok(LocusClassing {
        profiles: class_profiles,
        membership,
        class_owners,
        class_charges,
        substage: [
            parents_seconds,
            sweeps_seconds,
            partials_seconds,
            classing_seconds,
        ],
        split_seam_queries,
        ir_extractions,
    })
}

// ---------------------------------------------------------------------------
// Domain hygiene (owner decision (c)): per-class scorability.
// ---------------------------------------------------------------------------

/// A class is scorable iff its profile is nonempty AND every member allele
/// covers at least one read length of source (total segment length >= L):
/// empty-profile and sub-read-length mini-alleles are degenerate rows that
/// score exactly 0.0 under every model variant (no interior features, no
/// seams) — they stay in the domain as chaining pass-through structure but
/// are excluded from the scorable candidate set.
pub(super) fn scorable_classes(
    classing: &LocusClassing,
    ranges: &[genome::SpanningTraversal],
) -> Vec<bool> {
    let covered = |allele: usize| -> bool {
        ranges[allele]
            .segments
            .iter()
            .map(|segment| segment.end.saturating_sub(segment.start))
            .sum::<u64>()
            >= READ_LENGTH as u64
    };
    (0..classing.profiles.len())
        .map(|class| {
            !classing.profiles[class].is_empty()
                && (0..classing.membership.len())
                    .filter(|&allele| classing.membership[allele] == class)
                    .all(covered)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Stage 1: exhaustive local sweep.
// ---------------------------------------------------------------------------

pub(super) struct LocusSweep {
    /// Exhaustive class-pair loss table (triangular; every class pair of the
    /// locus, no prefilter).
    pub(super) table: Vec<f64>,
    pub(super) class_pairs: u64,
    pub(super) allele_pairs: u64,
    pub(super) viable_allele_pairs: u64,
    pub(super) best_loss: f64,
    /// Exact tie set at the best loss (class pairs; empty-set ties are
    /// singletons by construction).
    pub(super) best_class_pairs: Vec<[usize; 2]>,
    pub(super) second_loss: Option<f64>,
    /// Exact minimum pair loss over viable class pairs (admissible-bound
    /// input for the chain margin).
    pub(super) min_viable_loss: f64,
    pub(super) native_pair_loss: Option<f64>,
    pub(super) native_class_pair: Option<[usize; 2]>,
    /// Number of class pairs strictly better than the native pair.
    pub(super) native_rank: Option<u64>,
    /// Admissible per-locus floor: what the locus's own observations can
    /// explain at best (`locus_min_pair_lower_bound_folded`).
    pub(super) floor: f64,
    /// [min, q25, median, q75, max] of the class-pair loss distribution
    /// over the EXACTLY scored pairs (the admissibly-pruned beyond-margin
    /// mass is excluded; see `bound_pruned_pairs`).
    pub(super) quantiles: Vec<f64>,
    /// STEP-3 admissible-bound pruning: pairs whose derived lower bound
    /// provably exceeds the locus's retention margin (min viable pair loss
    /// + the adjacent seam-evidence swing — the DP's own retention rule),
    /// so their exact values are never computed; their table entries carry
    /// +INFINITY (the DP's non-retained path reads no value from them).
    pub(super) bound_pruned_pairs: u64,
    /// The retention margin the prune used (admissible, data-derived; equals
    /// the DP's `compute_local_margins` value bit-for-bit: the locus's
    /// min viable pair loss + the adjacent seam swings).
    pub(super) prune_margin: f64,
    /// Exact best over VIABLE class pairs (the chaining-relevant winner).
    pub(super) best_viable_loss: f64,
    pub(super) best_viable_class_pairs: Vec<[usize; 2]>,
}

/// The exhaustive local sweep at one locus: every class pair scored under the
/// routed M1 model, the full ranking's summary statistics, the exact best
/// tie set, and the native pair's position. Class-pair granularity is exact
/// over allele pairs (profile-identical alleles are one class — the
/// sufficiency principle); allele-pair counts are reported alongside.
/// Per-locus adjacent seam-evidence swings (both boundaries' [max-min]
/// over viable-linked class seams) — the STEP-3 sweep prune's admissible
/// margin input, identical to `compute_local_margins`' swing and the
/// phasing rebuild's `local_seam_swing_margins` swing.
pub(super) fn seam_swings_of(
    locus_count: usize,
    successors: &[Vec<Vec<(usize, f64)>>],
    viable: &[Vec<bool>],
) -> Vec<f64> {
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
            swing
        })
        .collect()
}

pub(super) fn exhaustive_local_sweep(
    folded: &LocusFolded,
    loss_tables: &[Vec<f64>],
    membership: &[usize],
    viable: &[bool],
    native_allele: Option<usize>,
    model: &ScoreModel,
    // STEP-3 admissible-bound input: the locus's adjacent seam-evidence
    // swing (both boundaries' [max-min] over viable-linked class seams —
    // the same quantity `compute_local_margins` uses). The sweep computes
    // the exact min viable pair loss first (bound-ordered branch-and-bound)
    // and then prunes every pair whose derived lower bound provably exceeds
    // min + swing (the DP's own retention rule) BEFORE scoring it.
    seam_swing: f64,
    // Per-class restricted interior-junction charges (charge-homogeneous
    // classes), folded into the class-pair table: the pair loss of two
    // alleles is the profile-pair loss plus each allele's additive novel-
    // junction charge, so the folded table stays exact at allele
    // granularity and every downstream statistic (best pairs, margins,
    // native rank, DP states) is charge-inclusive.
    class_charges: &[f64],
    // Domain hygiene (owner decision (c), 2026-09): per-class scorability.
    // A class is scorable iff its profile is nonempty AND every member
    // allele covers at least one read length of source (total segment
    // length >= L). Empty-profile and sub-read-length mini-alleles are
    // degenerate rows (they scored exactly 0 under the flat model too and
    // only ever lost because real pairs scored negative): they stay in the
    // domain as chaining pass-through structure but their pair entries are
    // +infinity, so they never win a local sweep, never enter a margin, and
    // never get selected.
    scorable: &[bool],
) -> io::Result<LocusSweep> {
    let classes = folded.prepared.len();
    ensure(classes > 0, "spine locus has no classes")?;
    ensure(scorable.len() == classes, "spine scorable cardinality mismatch")?;
    ensure(
        class_charges.len() == classes,
        "spine class charge cardinality mismatch"
    )?;
    ensure(
        seam_swing.is_finite() && seam_swing >= 0.0,
        "spine seam swing must be finite and nonnegative"
    )?;

    // ------------------------------------------------------------- per-class
    // admissible-bound summaries (STEP 3; derived entirely from the
    // per-feature loss arithmetic, no thresholds):
    //  - singleton loss S_c: the EXACT sum of the class's own per-feature
    //    terms (the terms the pair loss keeps for features the other class
    //    lacks);
    //  - union-minimum mass M_c and its positive part P_c: the per-feature
    //    global minima m_f over the locus's achievable count range and the
    //    maximal observed side (the loss is decreasing in observed support
    //    and concave in count, so the box minimum sits at an endpoint:
    //    m_f = min(0, c_f(2*max_f, 2*obs_f))), summed over the class's
    //    features (inclusion-exclusion over the pair's feature union);
    //  - interaction mass J_c: the per-feature UPPER bound on the shared-
    //    feature interaction I_f = c(qa,obsA) + c(qb,obsB) - c(qa+qb, .):
    //    every singleton term is maximized at the locus's MINIMAL observed
    //    side and every merged form is minimized at the MAXIMAL merged side
    //    (2*obs_f), so I_f is bounded by the achievable count-grid maximum
    //    I_f^up = max over qa,qb <= max_f of
    //    [c(qa,obs_min_f) + c(qb,obs_min_f) - c(qa+qb, 2*obs_f)].
    // Every quantity comes from `loss_fractional_entry`'s own form.
    let summaries_started = std::time::Instant::now();
    let per_count = model.histogram as f64 * model.depth / model.denominator;
    let mut m_by_fid = vec![0.0f64; folded.max_add.len()];
    let mut j_by_fid = vec![0.0f64; folded.max_add.len()];
    let mut obs_min_by_fid = vec![f64::INFINITY; folded.max_add.len()];
    for class in 0..classes {
        for &(_fid, _q, observed) in &folded.prepared[class] {
            if observed < obs_min_by_fid[_fid as usize] {
                obs_min_by_fid[_fid as usize] = observed;
            }
        }
    }
    for fid in 0..folded.max_add.len() as u32 {
        let max = folded.max_add[fid as usize];
        if max == 0 {
            continue;
        }
        let obs_max = folded.obs_by_fid[fid as usize];
        let obs_min = obs_min_by_fid[fid as usize];
        let entry = folded.entry(fid);
        let merged_floor = crate::loss_fractional_entry(
            model,
            2 * max,
            2.0 * obs_max,
            entry,
        )?;
        m_by_fid[fid as usize] = 0.0f64.min(merged_floor);
        // CLOSED-FORM interaction upper bound (derivation, straight from
        // c(q, obs) = q*s - w*obs*ln_1p(q*s/beta): floor the singleton
        // observed sides at obs_min (the loss is decreasing in observed
        // support), ceiling the merged side at 2*obs_max, and write
        // I(sa, sb) = w*[O*ln(1+u) - o*ln(1+u+v)] with u = (sa+sb)/beta,
        // v = sa*sb/beta^2 >= 0, O = 2*obs_max, o = obs_min. Along each
        // line sa+sb = t the maximum sits at v = 0 (one count zero), and
        // the resulting w*(O-o)*ln(1+u) is increasing in u, so the box
        // maximum is
        //   I_f^up = w*(2*obs_max - obs_min)*ln_1p(2*max_f*s_unit/beta).
        let beta = match entry {
            Some(entry) if entry.beta.is_finite() && entry.beta > 0.0 => entry.beta,
            _ => model.background,
        };
        let weight = match entry {
            Some(entry) => entry.weight,
            None => 1.0,
        };
        let interaction = weight
            * (2.0 * obs_max - obs_min)
                .max(0.0)
            * (2.0 * max as f64 * per_count / beta).ln_1p();
        j_by_fid[fid as usize] = interaction;
    }
    let grid_evaluations: u64 = folded
        .max_add
        .iter()
        .map(|&max| (max as u64 + 1) * (max as u64 + 1))
        .sum();
    let summaries_grid_seconds = summaries_started.elapsed().as_secs_f64();
    let summaries_started = std::time::Instant::now();
    let mut class_singleton = vec![0.0f64; classes];
    let mut class_min_mass = vec![0.0f64; classes];
    let mut class_min_pos = vec![0.0f64; classes];
    let mut class_interaction = vec![0.0f64; classes];
    for class in 0..classes {
        if !scorable[class] {
            continue;
        }
        let mut singleton = 0.0f64;
        let mut min_mass = 0.0f64;
        let mut min_pos = 0.0f64;
        let mut interaction = 0.0f64;
        for &(fid, q, observed) in &folded.prepared[class] {
            singleton += crate::folded_term(
                fid,
                q,
                observed,
                folded.entry(fid),
                loss_tables,
                model,
            )?;
            min_mass += m_by_fid[fid as usize];
            if m_by_fid[fid as usize] > 0.0 {
                min_pos += m_by_fid[fid as usize];
            }
            interaction += j_by_fid[fid as usize];
        }
        class_singleton[class] = singleton + class_charges[class];
        class_min_mass[class] = min_mass + class_charges[class];
        class_min_pos[class] = min_pos;
        class_interaction[class] = interaction;
    }

    let summaries_class_seconds = summaries_started.elapsed().as_secs_f64();
    let summaries_seconds = summaries_grid_seconds + summaries_class_seconds;
    eprintln!(
        "[sweep] summaries: grid evals {grid_evaluations} ({summaries_grid_seconds:.2}s), classes ({summaries_class_seconds:.2}s), total {summaries_seconds:.2}s"
    );
    let pass1_started = std::time::Instant::now();
    // Bound-ordered enumeration order: classes ascending by the singleton-
    // charged summary (the tighter anchor; best candidates first so the
    // branch-and-bound incumbent improves fast and rows terminate early).
    let mut order: Vec<usize> = (0..classes).filter(|&c| scorable[c]).collect();
    order.sort_by(|&a, &b| {
        class_singleton[a]
            .total_cmp(&class_singleton[b])
            .then(a.cmp(&b))
    });

    // The exact pair evaluation (bit-identical to the old every-pair sweep).
    let evaluate = |first: usize, second: usize| -> io::Result<f64> {
        let mut value = prepared_pair_loss_folded(
            &folded.prepared[first],
            &folded.prepared[second],
            folded.owners[first] == folded.owners[second],
            &folded.entries_by_fid,
            loss_tables,
            model,
        )?;
        value += class_charges[first] + class_charges[second];
        Ok(value)
    };
    // The two derived lower bounds (both admissible; take the max).
    let bounds = |first: usize, second: usize| -> f64 {
        let interaction = class_interaction[first].min(class_interaction[second]);
        let singleton_bound = class_singleton[first] + class_singleton[second] - interaction;
        let positive = class_min_pos[first].min(class_min_pos[second]);
        let union_bound = class_min_mass[first] + class_min_mass[second] - positive;
        singleton_bound.max(union_bound)
    };
    // Row-terminating (monotone) weakenings of the bounds for `first` fixed:
    // drop the `min` with the row partner's masses (<= the row's own).
    let row_bound_singleton = |first: usize, second: usize| -> f64 {
        class_singleton[first] + class_singleton[second] - class_interaction[first]
    };
    let row_bound_union = |first: usize, second: usize| -> f64 {
        class_min_mass[first] + class_min_mass[second] - class_min_pos[first]
    };

    let mut table = vec![f64::INFINITY; classes * (classes + 1) / 2];
    let mut scored = vec![false; classes * (classes + 1) / 2];
    for class in 0..classes {
        if !scorable[class] {
            continue;
        }
        // Unscorable partners keep +INFINITY (the old unscorable entry).
        for partner in (class + 1)..classes {
            if !scorable[partner] {
                scored[class_pair_index(class, partner)] = true;
            }
        }
    }

    let mut members: Vec<Vec<usize>> = vec![Vec::new(); classes];
    let mut viable_members: Vec<Vec<usize>> = vec![Vec::new(); classes];
    for (allele, &class) in membership.iter().enumerate() {
        members[class].push(allele);
        if viable[allele] {
            viable_members[class].push(allele);
        }
    }

    let order_seconds = pass1_started.elapsed().as_secs_f64();
    eprintln!("[sweep] order build: {order_seconds:.2}s");
    let pass1_started = std::time::Instant::now();
    // ------------------------------------------------- pass 1: the exact core
    // (best/second/min-viable/native + their exact tie sets), bound-ordered
    // branch-and-bound. The prune threshold is shared across the parallel
    // rows through one monotone atomic: it holds the running
    // min(second-best-so-far, min-viable-so-far), a skipped pair's bound
    // exceeding a read (which is >= the final value, the quantities only
    // decrease) proves the pair out of every consumer's exact set.
    #[derive(Clone)]
    struct Core {
        best_loss: f64,
        second_loss: Option<f64>,
        min_viable_loss: f64,
        best_viable_loss: f64,
        best_pairs: Vec<[usize; 2]>,
        viable_pairs: Vec<[usize; 2]>,
    }
    impl Core {
        fn new() -> Self {
            Self {
                best_loss: f64::INFINITY,
                second_loss: None,
                min_viable_loss: f64::INFINITY,
                best_viable_loss: f64::INFINITY,
                best_pairs: Vec::new(),
                viable_pairs: Vec::new(),
            }
        }
    }
    fn track(
        first: usize,
        second: usize,
        loss: f64,
        viable_pair: bool,
        core: &mut Core,
    ) {
        if viable_pair {
            if (core.min_viable_loss.total_cmp(&loss)).is_gt() {
                core.min_viable_loss = loss;
            }
            if (core.best_viable_loss.total_cmp(&loss)).is_gt() {
                core.best_viable_loss = loss;
                core.viable_pairs.clear();
                core.viable_pairs.push([first, second]);
            } else if loss.to_bits() == core.best_viable_loss.to_bits() {
                core.viable_pairs.push([first, second]);
            }
        }
        if (core.best_loss.total_cmp(&loss)).is_gt() {
            core.second_loss = Some(core.best_loss);
            core.best_loss = loss;
            core.best_pairs.clear();
            core.best_pairs.push([first, second]);
        } else if loss.to_bits() == core.best_loss.to_bits() {
            core.best_pairs.push([first, second]);
        } else if core
            .second_loss
            .is_none_or(|s| (s.total_cmp(&loss)).is_gt())
        {
            core.second_loss = Some(loss);
        }
    }
    let mut core = Core::new();
    let core_evaluations = AtomicU64::new(0);
    // Seed the incumbent with the native pair (the backbone's own class
    // pair — also the DP's structurally-retained pair).
    if let Some(native) = native_allele {
        let class = membership[native];
        let loss = evaluate(class, class)?;
        let index = class_pair_index(class, class);
        table[index] = loss;
        scored[index] = true;
        core_evaluations.fetch_add(1, Ordering::Relaxed);
        let viable_pair = !viable_members[class].is_empty();
        track(class, class, loss, viable_pair, &mut core);
    }
    // Two shared monotone thresholds: the running second-best (best/second
    // exactness, every pair) and the running min viable loss (margin
    // exactness, viable pairs). A pair is skip-SAFE iff its bound exceeds
    // every threshold its exactness depends on: non-viable pairs need
    // bound > second; viable pairs need bound > second AND bound >
    // min_viable (skipping must be safe for BOTH statistics). A racy read
    // only ever sees a staler (larger) value, costing work, never a wrong
    // skip.
    let shared_second = std::sync::atomic::AtomicU64::new(f64::INFINITY.to_bits());
    let shared_min_viable = std::sync::atomic::AtomicU64::new(f64::INFINITY.to_bits());
    let lower_shared = |cell: &std::sync::atomic::AtomicU64, value: f64| {
        let bits = value.to_bits();
        let mut current = cell.load(Ordering::Relaxed);
        while bits < current {
            match cell.compare_exchange_weak(
                current,
                bits,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(observed) => current = observed,
            }
        }
    };
    lower_shared(&shared_second, core.second_loss.unwrap_or(f64::INFINITY));
    lower_shared(&shared_min_viable, core.min_viable_loss);
    let results: Vec<io::Result<(Core, Vec<(usize, f64)>)>> = order
        .par_iter()
        .copied()
        .map(|first| {
            let mut local = Core::new();
            let mut computed: Vec<(usize, f64)> = Vec::new();
            for &second in order.iter() {
                if second < first {
                    continue;
                }
                let index = class_pair_index(first, second);
                if scored[index] {
                    continue;
                }
                let second_so_far = f64::from_bits(shared_second.load(Ordering::Relaxed));
                let min_viable_so_far =
                    f64::from_bits(shared_min_viable.load(Ordering::Relaxed));
                // Row termination: the monotone bounds grow with the
                // partner's anchors; past BOTH thresholds the whole
                // remainder of the row (mixed viability) is provably
                // excluded.
                if row_bound_singleton(first, second) > second_so_far
                    && row_bound_singleton(first, second) > min_viable_so_far
                    && row_bound_union(first, second) > second_so_far
                    && row_bound_union(first, second) > min_viable_so_far
                {
                    break;
                }
                let viable_pair = !viable_members[first].is_empty()
                    && !viable_members[second].is_empty();
                let skip_safe = if viable_pair {
                    bounds(first, second) > second_so_far
                        && bounds(first, second) > min_viable_so_far
                } else {
                    bounds(first, second) > second_so_far
                };
                if skip_safe {
                    continue;
                }
                let loss = evaluate(first, second)?;
                core_evaluations.fetch_add(1, Ordering::Relaxed);
                computed.push((index, loss));
                track(first, second, loss, viable_pair, &mut local);
                lower_shared(&shared_second, local.second_loss.unwrap_or(f64::INFINITY));
                lower_shared(&shared_min_viable, local.min_viable_loss);
            }
            Ok((local, computed))
        })
        .collect();
    for result in &results {
        let (_local, computed) = result.as_ref().expect("row sweep");
        for &(index, loss) in computed {
            table[index] = loss;
            scored[index] = true;
        }
    }
    // Merge the per-row cores into the sweep's exact statistics (the global
    // best is the min over row bests with the row tie sets; the global
    // second is the best of the merged runner-ups and non-winning row
    // bests; min-viable is the min over rows).
    let mut best_class_pairs: Vec<[usize; 2]> = Vec::new();
    let mut best_viable_class_pairs: Vec<[usize; 2]> = Vec::new();
    let mut second_candidates: Vec<f64> = Vec::new();
    for result in &results {
        let row = result.as_ref().map(|(local, _)| local).expect("row sweep");
        if (core.min_viable_loss.total_cmp(&row.min_viable_loss)).is_gt() {
            core.min_viable_loss = row.min_viable_loss;
        }
        if (core.best_viable_loss.total_cmp(&row.best_viable_loss)).is_gt() {
            core.best_viable_loss = row.best_viable_loss;
            best_viable_class_pairs.clear();
        }
        if row.best_viable_loss.to_bits() == core.best_viable_loss.to_bits() {
            best_viable_class_pairs.extend_from_slice(&row.viable_pairs);
        }
        if (core.best_loss.total_cmp(&row.best_loss)).is_gt() {
            core.best_loss = row.best_loss;
        }
        if let Some(row_second) = row.second_loss {
            second_candidates.push(row_second);
        }
        if row.best_loss.to_bits() != core.best_loss.to_bits() {
            second_candidates.push(row.best_loss);
        }
    }
    for result in &results {
        let row = result.as_ref().map(|(local, _)| local).expect("row sweep");
        if row.best_loss.to_bits() == core.best_loss.to_bits() {
            best_class_pairs.extend_from_slice(&row.best_pairs);
        }
    }
    if let Some(native) = native_allele {
        let class = membership[native];
        let loss = table[class_pair_index(class, class)];
        if loss.to_bits() == core.best_loss.to_bits()
            && !best_class_pairs.contains(&[class, class])
        {
            best_class_pairs.push([class, class]);
        }
        if !viable_members[class].is_empty()
            && loss.to_bits() == core.best_viable_loss.to_bits()
            && !best_viable_class_pairs.contains(&[class, class])
        {
            best_viable_class_pairs.push([class, class]);
        }
    }
    for value in second_candidates {
        if value.to_bits() == core.best_loss.to_bits() {
            continue;
        }
        if core
            .second_loss
            .is_none_or(|s| (s.total_cmp(&value)).is_gt())
        {
            core.second_loss = Some(value);
        }
    }
    let core_evaluations = core_evaluations.load(Ordering::Relaxed);
    // The old sweep's deterministic tie order: (first, second) id order.
    best_class_pairs.sort_unstable();
    best_class_pairs.dedup();
    best_viable_class_pairs.sort_unstable();
    best_viable_class_pairs.dedup();
    let best_loss = core.best_loss;
    let second_loss = core.second_loss;
    let min_viable_loss = core.min_viable_loss;
    let best_viable_loss = core.best_viable_loss;

    ensure(best_loss.is_finite(), "spine sweep found no scored pair")?;
    ensure(
        min_viable_loss.is_finite(),
        "spine locus has no viable class pair"
    )?;
    let pass1_seconds = pass1_started.elapsed().as_secs_f64();
    eprintln!("[sweep] pass 1 (exact core): {pass1_seconds:.2}s");
    let pass2_started = std::time::Instant::now();
    // The retention margin (identical to compute_local_margins' rule).
    let prune_margin = min_viable_loss + seam_swing;

    // --------------------------------------- pass 2: the margin fill (STEP 3).
    // Every pair whose derived lower bound provably exceeds the retention
    // margin is left at +INFINITY (the DP's own non-retained path reads no
    // value from it; the retention predicate, the stage-3 admission and
    // every margin consumer see exactly the retained/not-retained decision
    // the full table produced). Pairs inside the possible-retention range
    // get exact scores (parallel over rows; the threshold is fixed).
    let rows: Vec<io::Result<(Vec<(usize, f64)>, u64, u64)>> = order
        .par_iter()
        .copied()
        .map(|first| {
            let mut computed: Vec<(usize, f64)> = Vec::new();
            let mut evaluations = 0u64;
            let mut skips = 0u64;
            for &second in order.iter() {
                if second < first {
                    continue;
                }
                let index = class_pair_index(first, second);
                if scored[index] {
                    continue;
                }
                // Row termination against the fixed margin.
                if row_bound_singleton(first, second) > prune_margin + BOUND_PRUNE_EPSILON
                    && row_bound_union(first, second) > prune_margin + BOUND_PRUNE_EPSILON
                {
                    break;
                }
                if bounds(first, second) > prune_margin + BOUND_PRUNE_EPSILON {
                    skips += 1;
                    continue;
                }
                let loss = evaluate(first, second)?;
                evaluations += 1;
                computed.push((index, loss));
            }
            Ok((computed, evaluations, skips))
        })
        .collect();
    let mut margin_evaluations = 0u64;
    let mut bound_pruned_pairs = 0u64;
    for result in rows {
        let (computed, evaluations, skips) = result?;
        margin_evaluations += evaluations;
        bound_pruned_pairs += skips;
        for (index, loss) in computed {
            table[index] = loss;
            scored[index] = true;
        }
    }

    let pass2_seconds = pass2_started.elapsed().as_secs_f64();
    eprintln!("[sweep] pass 2 (margin fill): {pass2_seconds:.2}s");
    // ------------------------------------------------------- statistics.
    // The counting statistics are over the same pair population as before
    // (every scorable pair); the loss-distribution statistics are over the
    // EXACTLY scored pairs (the pruned mass is provably beyond every
    // consumer's decision range; reported separately).
    let mut class_pairs = 0u64;
    let mut allele_pairs = 0u64;
    let mut viable_allele_pairs = 0u64;
    let mut losses: Vec<f64> = Vec::new();
    for &first in &order {
        let m1 = members[first].len() as u64;
        if m1 == 0 {
            continue;
        }
        let v1 = viable_members[first].len() as u64;
        for &second in order.iter() {
            if second < first {
                continue;
            }
            let m2 = members[second].len() as u64;
            if m2 == 0 {
                continue;
            }
            class_pairs += 1;
            allele_pairs += if first == second {
                m1 * (m1 + 1) / 2
            } else {
                m1 * m2
            };
            let v2 = viable_members[second].len() as u64;
            if v1 > 0 && v2 > 0 {
                viable_allele_pairs += if first == second {
                    v1 * (v1 + 1) / 2
                } else {
                    v1 * v2
                };
            }
            let index = class_pair_index(first, second);
            if scored[index] {
                losses.push(table[index]);
            }
        }
    }
    let (native_pair_loss, native_class_pair, native_rank) = match native_allele {
        Some(native) => {
            let class = membership[native];
            let loss = table[class_pair_index(class, class)];
            let rank = losses
                .iter()
                .filter(|&&value| (value.total_cmp(&loss)).is_lt())
                .count() as u64;
            (Some(loss), Some([class, class]), Some(rank))
        }
        None => (None, None, None),
    };
    ensure(!losses.is_empty(), "spine sweep scored no pair")?;
    losses.sort_by(|a, b| a.total_cmp(b));
    let quantile = |fraction: f64| -> f64 {
        let index = ((losses.len() as f64 - 1.0) * fraction).round() as usize;
        losses[index.min(losses.len() - 1)]
    };
    let quantiles = vec![
        losses[0],
        quantile(0.25),
        quantile(0.5),
        quantile(0.75),
        losses[losses.len() - 1],
    ];
    eprintln!(
        "[sweep] bound prune: class pairs {class_pairs}, exact {} (core {core_evaluations} \
         + margin {margin_evaluations}), bound-pruned {bound_pruned_pairs}",
        losses.len()
    );
    Ok(LocusSweep {
        table,
        class_pairs,
        allele_pairs,
        viable_allele_pairs,
        best_loss,
        best_class_pairs,
        second_loss,
        min_viable_loss,
        native_pair_loss,
        native_class_pair,
        native_rank,
        floor: locus_min_pair_lower_bound_folded(folded, model)
            + 2.0 * class_charges.iter().copied().fold(0.0f64, f64::min).min(0.0),
        quantiles,
        best_viable_loss,
        best_viable_class_pairs,
        bound_pruned_pairs,
        prune_margin,
    })
}

// ---------------------------------------------------------------------------
// Native backbone chain (the incumbent U for the exact margin).
// ---------------------------------------------------------------------------

/// The same-source-contiguous native chain, scored through the class tables
/// and seam links. Its total is the incumbent every retained pair must be
/// able to beat — the admissible-margin input. Faithful extraction of the
/// quarantined `run_routed_dp` backbone block.
pub(super) fn native_backbone_chain(
    ranges: &[Vec<genome::SpanningTraversal>],
    successors: &[Vec<Vec<(usize, f64)>>],
    locus_classes: &[LocusClassing],
    routed_equal: &crate::RoutedObs,
    component_locus_to_partition: &[u32],
    target_source: usize,
    model: &ScoreModel,
    backgrounds: &FeatureBackgrounds,
) -> io::Result<Option<(Vec<usize>, Vec<f64>)>> {
    let locus_count = ranges.len();
    let first: Vec<Option<usize>> = (0..locus_count)
        .map(|locus| {
            ranges[locus].iter().position(|traversal| {
                traversal.segments.len() == 1
                    && traversal.segments[0].source == target_source
                    && !traversal.segments[0].reverse
                    && traversal.segments[0].start < traversal.segments[0].end
            })
        })
        .collect();
    let mut chain: Vec<usize> = Vec::new();
    let mut cumulative: Vec<f64> = Vec::new();
    let mut total = 0.0f64;
    for locus in 0..locus_count {
        let native = if locus == 0 {
            first[0]
        } else {
            chain.last().and_then(|&previous| {
                successors[locus - 1][previous]
                    .iter()
                    .filter(|&&(next, _)| first[locus] == Some(next))
                    .map(|&(next, score)| (next, score))
                    .next()
                    .map(|(next, _)| next)
            }).or_else(|| {
                chain.last().and_then(|&previous| {
                    successors[locus - 1][previous]
                        .iter()
                        .find(|&&(next, _)| {
                            ranges[locus][next].segments.len() == 1
                                && ranges[locus][next].segments[0].source == target_source
                                && !ranges[locus][next].segments[0].reverse
                                && ranges[locus][next].segments[0].start
                                    < ranges[locus][next].segments[0].end
                        })
                        .map(|&(next, _)| next)
                })
            })
        };
        let Some(native) = native else {
            return Ok(None);
        };
        let class = locus_classes[locus].membership[native];
        let profile = &locus_classes[locus].profiles[class];
        let owner = locus_classes[locus].class_owners[class];
        total += merged_pair_loss_multiplicity(
            profile,
            profile,
            crate::owner_routed_obs(routed_equal, owner, component_locus_to_partition),
            model,
            backgrounds,
        )?;
        if locus > 0 {
            let previous = *chain.last().expect("chain");
            let seam = successors[locus - 1][previous]
                .iter()
                .find(|&&(next, _)| next == native)
                .map(|&(_, score)| score)
                .unwrap_or(f64::NAN);
            total += 2.0 * seam;
        }
        chain.push(native);
        cumulative.push(total);
    }
    Ok(Some((chain, cumulative)))
}

// ---------------------------------------------------------------------------
// Assessment-side truth decomposition (geometric local model).
// ---------------------------------------------------------------------------

/// Per (locus, copy): the reference route's merged route-contiguous pieces
/// intersected with the locus's public BED territory rows (the
/// `score_reference_m1` decomposition).
pub(super) fn reference_local_piece_lists(
    territory: &[Vec<SourceRange>],
    route_pair: [&routes::Route; 2],
    locus_count: usize,
    locus_offset: usize,
) -> Vec<[Vec<(usize, u64, u64, bool, u32)>; 2]> {
    // Pieces carry the OWNING universe partition of the territory row they
    // intersect (the window-domain extension's owner-resolved charging;
    // every pre-extension territory row is owned by the locus's axis
    // partition, so the pieces are bit-identical there).
    let mut local_pieces: Vec<[Vec<(usize, u64, u64, bool, u32)>; 2]> =
        vec![[Vec::new(), Vec::new()]; locus_count];
    for (copy, route) in route_pair.iter().enumerate() {
        for segment in &route.segments {
            for locus in 0..locus_count {
                for interval in &territory[locus_offset + locus] {
                    if interval.source != segment.source {
                        continue;
                    }
                    let lo = segment.start.max(interval.start);
                    let hi = segment.end.min(interval.end);
                    if lo < hi {
                        local_pieces[locus][copy].push((
                            segment.source,
                            lo,
                            hi,
                            segment.reverse,
                            interval.partition as u32,
                        ));
                    }
                }
            }
        }
    }
    for locus in 0..locus_count {
        for copy in 0..2 {
            let raw_pieces = local_pieces[locus][copy].clone();
            let mut merged: Vec<(usize, u64, u64, bool, u32)> = Vec::new();
            for piece in raw_pieces {
                match merged.last_mut() {
                    Some(last)
                        if last.0 == piece.0
                            && last.3 == piece.3
                            && last.4 == piece.4
                            && last.2 >= piece.1
                            && last.2 <= piece.2 =>
                    {
                        last.2 = last.2.max(piece.2);
                    }
                    _ => merged.push(piece),
                }
            }
            local_pieces[locus][copy] = merged;
        }
    }
    local_pieces
}

/// Geometric profile of one arbitrary truth piece interval (standalone
/// extraction: contained path anchors plus the RC-view records, jointly
/// swept — the same machinery as a parent interval).
fn geometric_piece_profile(
    panel: &SyngIndex,
    sources: &routes::Sources,
    path_of_source: &[usize],
    source: usize,
    start: u64,
    end: u64,
    k: u64,
    memo: &mut HashMap<(usize, u64, u64), Profile>,
) -> io::Result<Profile> {
    if let Some(profile) = memo.get(&(source, start, end)) {
        return Ok(profile.clone());
    }
    let anchors = contained_path_anchors(panel, path_of_source[source], start, end, k)?;
    let records = geometric_ir_records(panel, sources, source, start, end, k)?;
    let profile = geometric_sweep_profile(
        &anchors,
        &records,
        end - start,
        k,
        &format!("spine-piece:{source}:{start}-{end}"),
    )?;
    memo.insert((source, start, end), profile.clone());
    Ok(profile)
}

/// Per-locus local M1 pair losses of a reference pair under the geometric
/// model (piece profiles plus the surviving L149 internal-seam machinery
/// between pieces within a locus), charged against the partition's routed
/// shares. Under the junction-evidence restriction: a piece-to-piece
/// interior junction that is a co-occurring real panel junction keeps the
/// merged (pooled) seam profile bit-identically; a NOVEL interior junction
/// (the reference's real recombination points) is paid by the reads that
/// actually cross it — its seam profile is charged restricted and added
/// outside the merged profile, and its crossing census is reported.
pub(super) fn spine_reference_local_losses(
    panel: &SyngIndex,
    sources: &routes::Sources,
    path_of_source: &[usize],
    piece_lists: &[[Vec<(usize, u64, u64, bool, u32)>; 2]],
    routed_equal: &crate::RoutedObs,
    component_locus_to_partition: &[u32],
    k: u64,
    model: &ScoreModel,
    span: &JunctionSpanIndex,
    mut backgrounds: Option<&mut FeatureBackgrounds>,
) -> io::Result<(Vec<f64>, Vec<serde_json::Value>)> {
    let locus_count = piece_lists.len();
    let mut piece_memo: HashMap<(usize, u64, u64), Profile> = HashMap::new();
    let mut losses = Vec::with_capacity(locus_count);
    let mut census: Vec<serde_json::Value> = Vec::new();
    // Features charged here that the background scan has not measured yet
    // (collected only under `Some(backgrounds)`); when any exist the scan is
    // extended and the scoring recurses once on the complete backgrounds.
    let mut missing: Vec<FeatureKey> = Vec::new();
    for locus in 0..locus_count {
        // Per copy: the piece stream decomposed into OWNER RUNS (maximal
        // same-owner stretches; interior co-occurring seams merge into the
        // run's profile exactly as the pre-extension copy profile did).
        let mut copy_runs: [Vec<(u32, Profile)>; 2] = [Vec::new(), Vec::new()];
        let mut copy_charges: [f64; 2] = [0.0; 2];
        for copy in 0..2 {
            let mut runs: Vec<(u32, Profile)> = Vec::new();
            let mut previous: Option<(Vec<u8>, (usize, u64, u64, bool, u32))> = None;
            for &(source, start, end, reverse, owner) in &piece_lists[locus][copy] {
                ensure(
                    !reverse,
                    "spine reference decomposition met a reverse piece (unsupported)"
                );
                let sequence = sources.fetch(source, start, end)?;
                let piece_profile = geometric_piece_profile(
                    panel,
                    sources,
                    path_of_source,
                    source,
                    start,
                    end,
                    k,
                    &mut piece_memo,
                )?;
                if let Some((prev, prev_piece)) = &previous {
                    let (seam, _) =
                        genome::profile_event_seam(panel, prev, &sequence, READ_LENGTH, MAX_FEATURES)?;
                    let piece = (source, start, end, reverse, owner);
                    let cooccurring = junction::segment_pair_gap(
                        &junction::junction_range(prev_piece.0, prev_piece.1, prev_piece.2, prev_piece.3, u32::MAX),
                        &junction::junction_range(source, start, end, reverse, u32::MAX),
                    )
                    .is_some_and(|gap| gap >= 0);
                    if prev_piece.4 == owner {
                        // Same-owner interior junction: exactly the
                        // pre-extension handling.
                        let run = runs.last_mut().expect("same-owner piece continues a run");
                        run.1 = merge_profiles(&[&run.1, &piece_profile])?;
                        if cooccurring {
                            run.1 = merge_profiles(&[&run.1, &seam])?;
                        } else {
                            let owner_partition = crate::owner_universe_partition(
                                owner,
                                component_locus_to_partition,
                            );
                            let outcome = span.restricted_charge(
                                &seam,
                                &[(
                                    junction::junction_range(prev_piece.0, prev_piece.1, prev_piece.2, prev_piece.3, owner_partition),
                                    junction::junction_range(source, start, end, reverse, owner_partition),
                                )],
                                sources,
                                path_of_source,
                                &[owner_partition],
                                model,
                            )?;
                            copy_charges[copy] += outcome.charge;
                            census.push(serde_json::json!({
                                "locus": locus,
                                "copy": copy,
                                "kind": "novel",
                                "left": [prev_piece.0, prev_piece.1, prev_piece.2, prev_piece.3],
                                "right": [source, start, end, reverse],
                                "spanning_reads": outcome.spanning_reads,
                                "events": outcome.events,
                                "span_features": outcome.span_feature_count,
                                "restricted_charge": outcome.charge,
                            }));
                        }
                    } else {
                        // Owner change inside the locus: the seam charges as
                        // a boundary across the two owners (pooled forms) or
                        // a restricted novel charge across both owners; the
                        // piece starts a new run.
                        if let Some(bg) = &backgrounds {
                            for key in seam.keys() {
                                if !bg.contains(key) {
                                    missing.push(key.clone());
                                }
                            }
                        }
                        let obs_prev = crate::owner_routed_obs(
                            routed_equal,
                            prev_piece.4,
                            component_locus_to_partition,
                        );
                        let obs_this = crate::owner_routed_obs(
                            routed_equal,
                            owner,
                            component_locus_to_partition,
                        );
                        if cooccurring {
                            let mut merged: HashMap<FeatureKey, f64> = HashMap::new();
                            for map in [obs_prev, obs_this] {
                                for (feature, share) in map {
                                    *merged.entry(feature.clone()).or_default() += *share;
                                }
                            }
                            let charge = match &backgrounds {
                                Some(bg) => {
                                    profile_loss_boundary_multiplicity_owned(
                                        &seam,
                                        &merged,
                                        &crate::EMPTY_ROUTED_OBS,
                                        false,
                                        model,
                                        bg,
                                    )?
                                }
                                None => profile_loss_boundary_owned(
                                    &seam,
                                    &merged,
                                    &crate::EMPTY_ROUTED_OBS,
                                    false,
                                    model,
                                )?,
                            };
                            copy_charges[copy] += charge;
                        } else {
                            let prev_partition = crate::owner_universe_partition(
                                prev_piece.4,
                                component_locus_to_partition,
                            );
                            let this_partition = crate::owner_universe_partition(
                                owner,
                                component_locus_to_partition,
                            );
                            let outcome = span.restricted_charge(
                                &seam,
                                &[(
                                    junction::junction_range(prev_piece.0, prev_piece.1, prev_piece.2, prev_piece.3, prev_partition),
                                    junction::junction_range(source, start, end, reverse, this_partition),
                                )],
                                sources,
                                path_of_source,
                                &[prev_partition, this_partition],
                                model,
                            )?;
                            copy_charges[copy] += outcome.charge;
                            census.push(serde_json::json!({
                                "locus": locus,
                                "copy": copy,
                                "kind": "novel",
                                "left": [prev_piece.0, prev_piece.1, prev_piece.2, prev_piece.3],
                                "right": [source, start, end, reverse],
                                "spanning_reads": outcome.spanning_reads,
                                "events": outcome.events,
                                "span_features": outcome.span_feature_count,
                                "restricted_charge": outcome.charge,
                                "owner_boundary": true,
                            }));
                        }
                        runs.push((owner, piece_profile));
                    }
                } else {
                    runs.push((owner, piece_profile));
                }
                previous = Some((sequence, (source, start, end, reverse, owner)));
            }
            copy_runs[copy] = runs;
        }
        // Per-locus charge: the copies' owner runs grouped by owning
        // partition; each owner's merged copy profiles charge that owner's
        // routed shares (single owner = the pre-extension merged charge,
        // bit-identically).
        let mut owners: BTreeSet<u32> = BTreeSet::new();
        for copy in 0..2 {
            for (owner, _) in &copy_runs[copy] {
                owners.insert(*owner);
            }
        }
        let mut loss = 0.0f64;
        for owner in owners {
            let obs = crate::owner_routed_obs(routed_equal, owner, component_locus_to_partition);
            let merged = |copy: usize| -> io::Result<Profile> {
                let mut profile = Profile::new();
                for (run_owner, run_profile) in &copy_runs[copy] {
                    if *run_owner == owner {
                        profile = merge_profiles(&[&profile, run_profile])?;
                    }
                }
                Ok(profile)
            };
            let first = merged(0)?;
            let second = merged(1)?;
            match &backgrounds {
                Some(bg) => {
                    for profile in [&first, &second] {
                        for key in profile.keys() {
                            if !bg.contains(key) {
                                missing.push(key.clone());
                            }
                        }
                    }
                    loss += merged_pair_loss_multiplicity(&first, &second, obs, model, bg)?;
                }
                None => loss += merged_pair_loss(&first, &second, obs, model)?,
            }
        }
        loss += copy_charges[0] + copy_charges[1];
        losses.push(loss);
    }
    if let Some(bg) = &mut backgrounds {
        if !missing.is_empty() {
            bg.scan_extend(panel, missing, k, model)?;
            return spine_reference_local_losses(
                panel,
                sources,
                path_of_source,
                piece_lists,
                routed_equal,
                component_locus_to_partition,
                k,
                model,
                span,
                Some(&mut **bg),
            );
        }
    }
    Ok((losses, census))
}

// ---------------------------------------------------------------------------
// Step-A measurement (owner-directed 2026-09-23): the per-feature
// decomposition of the HAPLOID SINGLE-ALLELE local table. NO MODEL CHANGE —
// this reads the same-run structures (routed shares, measured backgrounds,
// oracle/IR profiles — the rescore convention that produced the chain
// ladder) and decomposes, per locus and per feature, the local comparison
// between a previously selected chain's slot-0 rows and the truth
// reference's copy-0 (mosaic) pieces: which row predicts the feature, at
// what predicted signal s, with what observed routed share C, under what
// measured background beta_f, and the resulting per-feature loss
// contribution to each row.
// ---------------------------------------------------------------------------

/// FNV-1a digest of a feature key (a compact identity for the table; the
/// full key stays in the run's own structures).
fn decompose_feature_digest(feature: &FeatureKey) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for word in feature {
        for byte in word.to_le_bytes() {
            hash ^= byte as u64;
            hash = hash.wrapping_mul(0x100000001b3);
        }
    }
    hash
}

/// The predicted signal s_f of one feature count — the exact
/// `loss_fractional_entry` signal arithmetic (q * histogram * depth /
/// denominator).
fn decompose_signal(model: &ScoreModel, q: u64) -> io::Result<f64> {
    let product = q
        .checked_mul(model.histogram)
        .ok_or_else(|| invalid("partition exposure product overflow"))?;
    ensure(
        product <= 1 << 53,
        "partition exposure conversion precision limit",
    )?;
    Ok(product as f64 * model.depth / model.denominator)
}

/// Sequence of one raw segment (fetch + reverse-complement when reversed).
fn decompose_segment_sequence(
    sources: &routes::Sources,
    segment: &SourceRange,
) -> io::Result<Vec<u8>> {
    let sequence = sources.fetch(segment.source, segment.start, segment.end)?;
    if segment.reverse {
        Ok(impg::graph::reverse_complement(&sequence))
    } else {
        Ok(sequence)
    }
}

/// Oracle (IR) profile of one raw row — the rescore convention: a single
/// segment spells its own interior profile; two segments spell both
/// interiors plus the interior seam (exactly `oracle_allele_profile`'s
/// construction, including the seam as charged features). When
/// `contain_window` is set (DIAGNOSTIC-ONLY, env IMPG_DECOMPOSE_CONTAIN —
/// the owner decision package's row-extent probe), each segment is first
/// clipped to the axis window's length (a prefix clip in the row's own
/// coordinates): the measured counterfactual for a window-containment
/// domain contract. No charging path reads this.
fn decompose_row_profile(
    panel: &SyngIndex,
    sources: &routes::Sources,
    segments: &[SourceRange],
    contain_window: Option<u64>,
) -> io::Result<Profile> {
    let clip = |segment: &SourceRange| -> SourceRange {
        match contain_window {
            Some(window) => {
                let mut clipped = segment.clone();
                clipped.end = segment.end.min(segment.start.saturating_add(window));
                clipped
            }
            None => segment.clone(),
        }
    };
    match segments {
        [segment] => {
            let segment = clip(segment);
            Ok(oracle_segment_profile(
                panel,
                &decompose_segment_sequence(sources, &segment)?,
            )?)
        }
        [left, right] => {
            let left = clip(left);
            let right = clip(right);
            let left_sequence = decompose_segment_sequence(sources, &left)?;
            let right_sequence = decompose_segment_sequence(sources, &right)?;;
            let interior_left = oracle_segment_profile(panel, &left_sequence)?;
            let interior_right = oracle_segment_profile(panel, &right_sequence)?;
            let (seam, _) = genome::profile_event_seam(
                panel,
                &left_sequence,
                &right_sequence,
                READ_LENGTH,
                MAX_FEATURES,
            )?;
            merge_profiles(&[&interior_left, &interior_right, &seam])
        }
        // Same-owner stitched chains (supervisor ruling,
        // genome/stitching-omission-alignment): the established convention
        // at every interior seam — the segments' profiles plus each
        // consecutive pair's profile_event_seam (identical to the oracle
        // profile the rescore charges the selected allele with).
        _ => {
            let mut parts: Vec<Profile> = Vec::with_capacity(segments.len() + 1);
            let mut sequences: Vec<Vec<u8>> = Vec::with_capacity(segments.len());
            for segment in segments {
                let segment = clip(segment);
                let sequence = decompose_segment_sequence(sources, &segment)?;
                parts.push(oracle_segment_profile(panel, &sequence)?);
                sequences.push(sequence);
            }
            for pair in sequences.windows(2) {
                let (seam, _) = genome::profile_event_seam(
                    panel,
                    &pair[0],
                    &pair[1],
                    READ_LENGTH,
                    MAX_FEATURES,
                )?;
                parts.push(seam);
            }
            merge_profiles(&parts.iter().collect::<Vec<&Profile>>())
        }
    }
}

/// The Step-A deliverable: the measured per-feature decomposition of the
/// haploid single-allele local comparison at every locus of the slice.
/// Writes the full table (compact JSON) to `out_path` and returns a small
/// reconciliation summary for the run log. The truth side replicates
/// `score_reference_m1`'s per-locus copy-0 computation for the haploid
/// reference (merged pieces, oracle piece profiles, co-occurring seams
/// merged, novel interior junctions restricted) INCLUDING its copy-0
/// boundary seam charges, so the grand total reconciles against the
/// reference ladder's `m1_oracle_rescore` bit-identically; the selected side
/// replicates the rescore's per-locus single-row charge
/// (`merged_pair_loss_multiplicity` with the empty second slot), so its
/// total reconciles against `m1_oracle_local_terms`.
#[allow(clippy::too_many_arguments)]
pub(super) fn decompose_haploid_local_table(
    panel: &SyngIndex,
    sources: &routes::Sources,
    path_of_source: &[usize],
    axis_slice: &[genome::AxisInterval],
    // The slice's territory intervals (indexed by the slice-local locus): the
    // selected chain's per-locus pieces are the scorer's territory
    // intersection of its route segments (owners from the territory rows).
    territory_slice: &[Vec<SourceRange>],
    selected_path: &std::path::Path,
    truth_pieces: &[[Vec<(usize, u64, u64, bool, u32)>; 2]],
    partition_obs: &[HashMap<FeatureKey, f64>],
    partitions: &[u32],
    // The routed shares + the component locus -> universe partition map
    // (the window-domain extension's owner-resolved observed sides).
    routed_equal: &crate::RoutedObs,
    // The record-once site-map builder (Fix 1's mixed-owner charging) and
    // the windows' full observed profiles (Fix 2's omission universes,
    // indexed by the slice-local locus).
    site: &crate::SiteObserved,
    window_obs: &[HashMap<FeatureKey, f64>],
    component_locus_to_partition: &[u32],
    ranges: &[Vec<genome::SpanningTraversal>],
    locus_classes: &[LocusClassing],
    donor_sources: &BTreeSet<usize>,
    pooled_support: Option<&HashMap<FeatureKey, f64>>,
    // The haploid track's CHARGED genome-wide cross-support background per
    // feature (Sum_r m_r*(1 - 1/t_r_genome)); the sample-side columns of the
    // table charge THIS form, so the per-feature arithmetic is verifiable
    // against the run's own rescore. The legacy component-scoped
    // complementary-share columns are kept alongside for the old-vs-new
    // delta distribution.
    sample_cross_support: Option<&HashMap<FeatureKey, f64>>,
    locus_offset: usize,
    k: u64,
    model: &ScoreModel,
    span: &JunctionSpanIndex,
    backgrounds: &mut FeatureBackgrounds,
    out_path: &std::path::Path,
) -> io::Result<serde_json::Value> {
    let started = Instant::now();
    let locus_count = partition_obs.len();
    // DIAGNOSTIC-ONLY (owner decision package): the containment probe —
    // clip the selected chain's rows to each locus's axis window length.
    let contain_probe = std::env::var("IMPG_DECOMPOSE_CONTAIN").is_ok();
    ensure(
        truth_pieces.len() == locus_count,
        "decompose truth piece cardinality mismatch",
    )?;
    // Parse the selected chain's route JSON (the rescore `selected_route`
    // format: per locus [slot-0 {segments}, "empty-slot2"]). The haploid
    // single-allele table requires the empty second slot everywhere.
    let selected_json: serde_json::Value = read_json(selected_path)?;
    let selected_rows = selected_json
        .as_array()
        .ok_or_else(|| invalid("decompose selected route must be an array"))?;
    ensure(
        selected_rows.len() == locus_count,
        "decompose selected route cardinality mismatch",
    )?;
    let mut exotic_segments: Vec<Vec<SourceRange>> = Vec::with_capacity(locus_count);
    for row in selected_rows {
        let slots = row
            .as_array()
            .ok_or_else(|| invalid("decompose selected route row must be an array"))?;
        ensure(slots.len() == 2, "decompose selected route row arity")?;
        ensure(
            slots[1].as_str() == Some("empty-slot2"),
            "decompose selected route slot 1 must be empty-slot2 \
             (the haploid single-allele table)",
        )?;
        let segments_value = slots[0]
            .get("segments")
            .and_then(serde_json::Value::as_array)
            .ok_or_else(|| invalid("decompose selected row lacks segments"))?;
        let mut segments = Vec::with_capacity(segments_value.len());
        for segment in segments_value {
            segments.push(SourceRange {
                partition: 0,
                occurrence: 0,
                source: segment
                    .get("source")
                    .and_then(serde_json::Value::as_u64)
                    .ok_or_else(|| invalid("decompose segment lacks source"))?
                    as usize,
                start: segment
                    .get("start")
                    .and_then(serde_json::Value::as_u64)
                    .ok_or_else(|| invalid("decompose segment lacks start"))?,
                end: segment
                    .get("end")
                    .and_then(serde_json::Value::as_u64)
                    .ok_or_else(|| invalid("decompose segment lacks end"))?,
                reverse: segment
                    .get("reverse")
                    .and_then(serde_json::Value::as_bool)
                    .ok_or_else(|| invalid("decompose segment lacks reverse"))?,
            });
        }
        ensure(
            !segments.is_empty(),
            "decompose selected row has no segments"
        )?;
        exotic_segments.push(segments);
    }
    // Owner-resolved observed sides (the window-domain extension): the
    // truth row charges the union-sum over its pieces' owners. The SELECTED
    // side's observed side is the CHARGED pieces' owners (the scorer's
    // form) — built after pass 1, once the side's own shared-once filter
    // has run. Anchor-only loci reproduce the pre-extension single map
    // bit-for-bit.
    let mut obs_exotic: Vec<HashMap<FeatureKey, f64>> =
        vec![HashMap::new(); locus_count];
    let mut obs_truth: Vec<HashMap<FeatureKey, f64>> = Vec::with_capacity(locus_count);
    for locus in 0..locus_count {
        let mut truth_owners: BTreeSet<u32> = BTreeSet::new();
        for piece in &truth_pieces[locus][0] {
            truth_owners.insert(piece.4);
        }
        let truth_map: HashMap<FeatureKey, f64> = if truth_owners.is_empty() {
            partition_obs[locus].clone()
        } else {
            // Fix 1 (record-once mixed-owner charging): the truth row's
            // owner-resolved observed side attributes each record ONCE —
            // the union-sum over the owners' full maps counted a record
            // touching k of them k times (the 157k decompose-vs-scorer
            // residual's source).
            site.site_map(
                truth_owners
                    .iter()
                    .map(|&owner| crate::owner_universe_partition(owner, component_locus_to_partition)),
            )
        };
        obs_truth.push(truth_map);
    }
    // Pass 1 — build every charged structure and collect the features the
    // background scan has not measured yet (mirroring `score_reference_m1`'s
    // missing-collection, extended ONCE before any charging).
    let mut piece_memo: HashMap<(usize, u64, u64, bool), Profile> = HashMap::new();
    let mut exotic_profiles: Vec<Profile> = Vec::with_capacity(locus_count);
    let mut truth_profiles: Vec<Profile> = Vec::with_capacity(locus_count);
    let mut truth_interior_charge = vec![0.0f64; locus_count];
    let mut truth_interior_census: Vec<Vec<serde_json::Value>> =
        vec![Vec::new(); locus_count];
    // Per boundary (index = right locus): the co-occurring seam profile to
    // charge pooled, or the novel-junction census row already charged.
    let mut truth_boundary_pooled: Vec<Option<Profile>> = vec![None; locus_count];
    let mut truth_boundary_charge = vec![0.0f64; locus_count];
    // The sample-side pooled boundary charges (the reconciliation layer).
    let mut truth_boundary_charge_sample = vec![0.0f64; locus_count];
    // The restricted novel-junction charges are background-invariant, so the
    // same values enter both the panel-form and the sample-side grands
    // (accumulated across loci for the totals' reconciliation rows).
    let mut truth_interior_total = 0.0f64;
    let mut truth_boundary_restricted_total = 0.0f64;
    let mut truth_boundary_census: Vec<Vec<serde_json::Value>> = vec![Vec::new(); locus_count];
    // The selected (exotic) side's own charge layers (the scorer's exact
    // form; see the pass-1 construction): interior restricted charges
    // (background-invariant, pass 1), cross-locus boundary seams (pooled
    // seams collected for the pass-2 sample charge; novel restricted in
    // pass 1), and the per-locus census rows.
    let mut exotic_interior_charge = vec![0.0f64; locus_count];
    let mut exotic_interior_census: Vec<Vec<serde_json::Value>> = vec![Vec::new(); locus_count];
    let mut exotic_interior_total = 0.0f64;
    let mut exotic_boundary_pooled: Vec<Option<Profile>> = vec![None; locus_count];
    let mut exotic_boundary_charge = vec![0.0f64; locus_count];
    let mut exotic_boundary_restricted_total = 0.0f64;
    let mut exotic_boundary_census: Vec<Vec<serde_json::Value>> = vec![Vec::new(); locus_count];
    let mut exotic_boundary_left_piece: Vec<Option<(usize, u64, u64, bool, u32)>> =
        vec![None; locus_count];
    let mut exotic_boundary_right_piece: Vec<Option<(usize, u64, u64, bool, u32)>> =
        vec![None; locus_count];
    let mut exotic_boundary_charge_sample = vec![0.0f64; locus_count];
    let mut exotic_boundary_pooled_total = 0.0f64;
    let mut missing: Vec<FeatureKey> = Vec::new();
    let mut previous_tail: Vec<u8> = Vec::new();
    let mut previous_piece: Option<(usize, u64, u64, bool, u32)> = None;
    // The selected (exotic) side's own scorer-form state (Ruling 1 rerun):
    // the side's own shared-once filter bookkeeping, its cross-locus
    // boundary state, and the per-locus CHARGED pieces' owner sets (the
    // observed side the scorer charges the charged runs against — the
    // charged-owners form; the former raw-segments form over-credited the
    // repeated pieces' predictions at their later windows, measured
    // 52,816 on this run's chain).
    let mut exotic_charged_pieces: std::collections::HashSet<(usize, u64, u64, bool)> =
        std::collections::HashSet::new();
    let mut exotic_previous_tail: Vec<u8> = Vec::new();
    let mut exotic_previous_piece: Option<(usize, u64, u64, bool, u32)> = None;
    let mut exotic_charged_owners: Vec<BTreeSet<u32>> = Vec::with_capacity(locus_count);
    // GROUP-SHARED MATERIAL CHARGED ONCE: the decompose must replicate
    // `score_reference_m1`'s once-per-material convention exactly (the
    // reference ladder runs with `charge_group_shared_once = true`), or the
    // truth side charges pieces the scorer charged at their first window —
    // a structural disagreement at the shared-piece loci (measured: a
    // ~414k decompose-vs-ladder residual on try 5). Same bookkeeping as
    // the scorer's per-copy `charged_pieces` set, copy 0.
    let mut charged_pieces: std::collections::HashSet<(usize, u64, u64, bool)> =
        std::collections::HashSet::new();
    // The truth boundary's two pieces per boundary (the owner-resolved
    // pooled boundary charge's observed sides): the previous locus's LAST
    // piece and this locus's FIRST piece (copy 0).
    let mut boundary_left_piece: Vec<Option<(usize, u64, u64, bool, u32)>> =
        vec![None; locus_count];
    let mut boundary_right_piece: Vec<Option<(usize, u64, u64, bool, u32)>> =
        vec![None; locus_count];
    for locus in 0..locus_count {
        // Selected (exotic) side: the SCORER's exact copy-0 construction
        // (the model of record the M1 comparison charges) — the route's
        // segments intersected with the locus's territory intervals (owners
        // from the territory rows), merged into route-contiguous pieces,
        // the side's OWN shared-once filter, oracle piece profiles,
        // co-occurring seams merged, novel interior junctions
        // restricted-charged (the cross-locus boundary seams below). The
        // former form here — `decompose_row_profile`'s segments+L149-seams
        // row profile against the raw segments' owners, unfiltered —
        // decomposed the DP's SURROGATE per-row charge and diverges from
        // the scorer wherever a piece repeats across windows (the filter
        // empties the later charged profile; the charged-owners observed
        // side narrows): measured 52,816 on this run's chain. The
        // containment probe keeps the old row profile (its own diagnostic
        // form, documented).
        let exotic = if contain_probe {
            let window_len = axis_slice[locus].end - axis_slice[locus].start;
            // The probe's own diagnostic form; the charged owners are the
            // raw segments' domain-row resolutions (the pre-stitching
            // fallback), so the downstream columns stay well-defined.
            let mut owners: BTreeSet<u32> = BTreeSet::new();
            for segment in &exotic_segments[locus] {
                if let Some(owner) = ranges[locus]
                    .iter()
                    .find(|row| {
                        row.segments.len() == 1
                            && row.segments[0].source == segment.source
                            && row.segments[0].start == segment.start
                            && row.segments[0].end == segment.end
                            && row.segments[0].reverse == segment.reverse
                    })
                    .map(crate::traversal_owner)
                {
                    owners.insert(crate::owner_universe_partition(
                        owner,
                        component_locus_to_partition,
                    ));
                }
            }
            exotic_charged_owners.push(owners);
            exotic_previous_tail = Vec::new();
            exotic_previous_piece = None;
            decompose_row_profile(
                panel,
                sources,
                &exotic_segments[locus],
                Some(window_len),
            )?
        } else {
            let mut raw_pieces: Vec<(usize, u64, u64, bool, u32)> = Vec::new();
            for segment in &exotic_segments[locus] {
                for interval in &territory_slice[locus] {
                    if interval.source != segment.source {
                        continue;
                    }
                    let lo = segment.start.max(interval.start);
                    let hi = segment.end.min(interval.end);
                    if lo < hi {
                        raw_pieces.push((
                            segment.source,
                            lo,
                            hi,
                            segment.reverse,
                            interval.partition as u32,
                        ));
                    }
                }
            }
            // The scorer's merge rule (same source and orientation with no
            // gap and the same owning partition).
            let mut merged_pieces: Vec<(usize, u64, u64, bool, u32)> = Vec::new();
            for piece in raw_pieces {
                match merged_pieces.last_mut() {
                    Some(last)
                        if last.0 == piece.0
                            && last.3 == piece.3
                            && last.4 == piece.4
                            && last.2 >= piece.1
                            && last.2 <= piece.2 =>
                    {
                        last.2 = last.2.max(piece.2);
                    }
                    _ => merged_pieces.push(piece),
                }
            }
            // The side's own shared-once filter in route order.
            let pieces: Vec<(usize, u64, u64, bool, u32)> = merged_pieces
                .iter()
                .filter(|piece| {
                    exotic_charged_pieces.insert((piece.0, piece.1, piece.2, piece.3))
                })
                .copied()
                .collect();
            let mut owners: BTreeSet<u32> = BTreeSet::new();
            let mut profile = Profile::new();
            let mut previous: Option<(Vec<u8>, (usize, u64, u64, bool, u32))> = None;
            let mut head: Option<Vec<u8>> = None;
            let mut head_piece: Option<(usize, u64, u64, bool, u32)> = None;
            let mut tail: Vec<u8> = Vec::new();
            for &(source, start, end, reverse, piece_owner) in &pieces {
                owners.insert(piece_owner);
                let sequence = sources.fetch(source, start, end)?;
                let sequence = if reverse {
                    impg::graph::reverse_complement(&sequence)
                } else {
                    sequence
                };
                let piece_profile = match piece_memo.get(&(source, start, end, reverse)) {
                    Some(profile) => profile.clone(),
                    None => {
                        let profile = oracle_segment_profile(panel, &sequence)?;
                        piece_memo.insert((source, start, end, reverse), profile.clone());
                        profile
                    }
                };
                profile = merge_profiles(&[&profile, &piece_profile])?;
                if let Some((previous_sequence, previous_seg)) = &previous {
                    let (seam, _) = genome::profile_event_seam(
                        panel,
                        previous_sequence,
                        &sequence,
                        READ_LENGTH,
                        MAX_FEATURES,
                    )?;
                    let cooccurring = junction::segment_pair_gap(
                        &junction::junction_range(
                            previous_seg.0,
                            previous_seg.1,
                            previous_seg.2,
                            previous_seg.3,
                            u32::MAX,
                        ),
                        &junction::junction_range(source, start, end, reverse, u32::MAX),
                    )
                    .is_some_and(|gap| gap >= 0);
                    if cooccurring {
                        profile = merge_profiles(&[&profile, &seam])?;
                    } else {
                        let piece_partition = crate::owner_universe_partition(
                            piece_owner,
                            component_locus_to_partition,
                        );
                        let prev_partition = crate::owner_universe_partition(
                            previous_seg.4,
                            component_locus_to_partition,
                        );
                        let outcome = span.restricted_charge(
                            &seam,
                            &[(
                                junction::junction_range(
                                    previous_seg.0,
                                    previous_seg.1,
                                    previous_seg.2,
                                    previous_seg.3,
                                    prev_partition,
                                ),
                                junction::junction_range(
                                    source, start, end, reverse, piece_partition,
                                ),
                            )],
                            sources,
                            path_of_source,
                            &[prev_partition, piece_partition],
                            model,
                        )?;
                        exotic_interior_charge[locus] += outcome.charge;
                        exotic_interior_total += outcome.charge;
                        exotic_interior_census[locus].push(serde_json::json!({
                            "kind": "interior",
                            "left": [previous_seg.0, previous_seg.1, previous_seg.2, previous_seg.3],
                            "right": [source, start, end, reverse],
                            "spanning_reads": outcome.spanning_reads,
                            "events": outcome.events,
                            "span_features": outcome.span_feature_count,
                            "restricted_charge": outcome.charge,
                        }));
                    }
                }
                if head.is_none() {
                    head = Some(sequence[..sequence.len().min(READ_LENGTH - 1)].to_vec());
                    head_piece = Some((source, start, end, reverse, piece_owner));
                }
                tail = sequence[sequence.len().saturating_sub(READ_LENGTH - 1)..].to_vec();
                previous = Some((sequence, (source, start, end, reverse, piece_owner)));
            }
            exotic_charged_owners.push(owners);
            // Cross-locus boundary seam (the scorer's boundary machinery):
            // pooled seams collected for the pass-2 sample charge; novel
            // junctions restricted-charged here (background-invariant).
            if locus > 0 && !exotic_previous_tail.is_empty() {
                exotic_boundary_left_piece[locus] = exotic_previous_piece;
                exotic_boundary_right_piece[locus] = head_piece;
                if let (Some(head), Some(head_seg), Some(prev_seg)) =
                    (&head, head_piece, exotic_previous_piece)
                {
                    let (seam, _) = genome::profile_event_seam(
                        panel,
                        &exotic_previous_tail,
                        head,
                        READ_LENGTH,
                        MAX_FEATURES,
                    )?;
                    let cooccurring = junction::segment_pair_gap(
                        &junction::junction_range(
                            prev_seg.0,
                            prev_seg.1,
                            prev_seg.2,
                            prev_seg.3,
                            u32::MAX,
                        ),
                        &junction::junction_range(
                            head_seg.0,
                            head_seg.1,
                            head_seg.2,
                            head_seg.3,
                            u32::MAX,
                        ),
                    )
                    .is_some_and(|gap| gap >= 0);
                    if cooccurring {
                        for key in seam.keys() {
                            if !backgrounds.contains(key) {
                                missing.push(key.clone());
                            }
                        }
                        exotic_boundary_pooled[locus] = Some(seam);
                    } else {
                        let prev_partition = crate::owner_universe_partition(
                            prev_seg.4,
                            component_locus_to_partition,
                        );
                        let head_partition = crate::owner_universe_partition(
                            head_seg.4,
                            component_locus_to_partition,
                        );
                        let outcome = span.restricted_charge(
                            &seam,
                            &[(
                                junction::junction_range(
                                    prev_seg.0,
                                    prev_seg.1,
                                    prev_seg.2,
                                    prev_seg.3,
                                    prev_partition,
                                ),
                                junction::junction_range(
                                    head_seg.0,
                                    head_seg.1,
                                    head_seg.2,
                                    head_seg.3,
                                    head_partition,
                                ),
                            )],
                            sources,
                            path_of_source,
                            &[prev_partition, head_partition],
                            model,
                        )?;
                        exotic_boundary_charge[locus] += outcome.charge;
                        exotic_boundary_restricted_total += outcome.charge;
                        exotic_boundary_census[locus].push(serde_json::json!({
                            "kind": "boundary",
                            "left": [prev_seg.0, prev_seg.1, prev_seg.2, prev_seg.3],
                            "right": [head_seg.0, head_seg.1, head_seg.2, head_seg.3],
                            "spanning_reads": outcome.spanning_reads,
                            "events": outcome.events,
                            "span_features": outcome.span_feature_count,
                            "restricted_charge": outcome.charge,
                        }));
                    }
                }
            }
            exotic_previous_tail = tail;
            exotic_previous_piece = pieces.last().copied();
            profile
        };
        for key in exotic.keys() {
            if !backgrounds.contains(key) {
                missing.push(key.clone());
            }
        }
        exotic_profiles.push(exotic);
        // Truth mosaic row: `score_reference_m1`'s copy-0 computation
        // (with its once-per-material piece filter — see the charged_pieces
        // preamble).
        let pieces: Vec<(usize, u64, u64, bool, u32)> = truth_pieces[locus][0]
            .iter()
            .copied()
            .filter(|piece| {
                charged_pieces.insert((piece.0, piece.1, piece.2, piece.3))
            })
            .collect();
        let pieces = &pieces;
        let mut profile = Profile::new();
        let mut previous: Option<(Vec<u8>, (usize, u64, u64, bool, u32))> = None;
        let mut head: Option<Vec<u8>> = None;
        let mut head_piece: Option<(usize, u64, u64, bool, u32)> = None;
        let mut tail: Vec<u8> = Vec::new();
        for &(source, start, end, reverse, piece_owner) in pieces {
            let sequence = sources.fetch(source, start, end)?;
            let sequence = if reverse {
                impg::graph::reverse_complement(&sequence)
            } else {
                sequence
            };
            let piece_profile = match piece_memo.get(&(source, start, end, reverse)) {
                Some(profile) => profile.clone(),
                None => {
                    let profile = oracle_segment_profile(panel, &sequence)?;
                    piece_memo.insert((source, start, end, reverse), profile.clone());
                    profile
                }
            };
            profile = merge_profiles(&[&profile, &piece_profile])?;
            if let Some((previous_sequence, previous_seg)) = &previous {
                let (seam, _) = genome::profile_event_seam(
                    panel,
                    previous_sequence,
                    &sequence,
                    READ_LENGTH,
                    MAX_FEATURES,
                )?;
                let piece = (source, start, end, reverse, piece_owner);
                let cooccurring = junction::segment_pair_gap(
                    &junction::junction_range(
                        previous_seg.0,
                        previous_seg.1,
                        previous_seg.2,
                        previous_seg.3,
                        u32::MAX,
                    ),
                    &junction::junction_range(source, start, end, reverse, u32::MAX),
                )
                .is_some_and(|gap| gap >= 0);
                if cooccurring {
                    profile = merge_profiles(&[&profile, &seam])?;
                } else {
                    let piece_partition = crate::owner_universe_partition(
                        piece_owner,
                        component_locus_to_partition,
                    );
                    let prev_partition = crate::owner_universe_partition(
                        previous_seg.4,
                        component_locus_to_partition,
                    );
                    let outcome = span.restricted_charge(
                        &seam,
                        &[(
                            junction::junction_range(
                                previous_seg.0,
                                previous_seg.1,
                                previous_seg.2,
                                previous_seg.3,
                                prev_partition,
                            ),
                            junction::junction_range(
                                source, start, end, reverse, piece_partition,
                            ),
                        )],
                        sources,
                        path_of_source,
                        &[prev_partition, piece_partition],
                        model,
                    )?;
                    truth_interior_charge[locus] += outcome.charge;
                    truth_interior_total += outcome.charge;
                    truth_interior_census[locus].push(serde_json::json!({
                        "kind": "interior",
                        "left": [previous_seg.0, previous_seg.1, previous_seg.2, previous_seg.3],
                        "right": [source, start, end, reverse],
                        "spanning_reads": outcome.spanning_reads,
                        "events": outcome.events,
                        "span_features": outcome.span_feature_count,
                        "restricted_charge": outcome.charge,
                    }));
                }
            }
            if head.is_none() {
                head = Some(sequence[..sequence.len().min(READ_LENGTH - 1)].to_vec());
                head_piece = Some((source, start, end, reverse, piece_owner));
            }
            tail = sequence[sequence.len().saturating_sub(READ_LENGTH - 1)..].to_vec();
            previous = Some((sequence, (source, start, end, reverse, piece_owner)));
        }
        for key in profile.keys() {
            if !backgrounds.contains(key) {
                missing.push(key.clone());
            }
        }
        // Boundary seam against the previous locus (copy 0 only — the
        // haploid reference's second slot is empty).
        if locus > 0 && !previous_tail.is_empty() {
            boundary_left_piece[locus] = previous_piece;
            boundary_right_piece[locus] = head_piece;
            if let (Some(head), Some(head_seg), Some(prev_seg)) =
                (&head, head_piece, previous_piece)
            {
                let (seam, _) = genome::profile_event_seam(
                    panel,
                    &previous_tail,
                    head,
                    READ_LENGTH,
                    MAX_FEATURES,
                )?;
                let cooccurring = junction::segment_pair_gap(
                    &junction::junction_range(
                        prev_seg.0, prev_seg.1, prev_seg.2, prev_seg.3, u32::MAX,
                    ),
                    &junction::junction_range(
                        head_seg.0, head_seg.1, head_seg.2, head_seg.3, u32::MAX,
                    ),
                )
                .is_some_and(|gap| gap >= 0);
                if cooccurring {
                    for key in seam.keys() {
                        if !backgrounds.contains(key) {
                            missing.push(key.clone());
                        }
                    }
                    truth_boundary_pooled[locus] = Some(seam);
                } else {
                    let prev_partition = crate::owner_universe_partition(
                        prev_seg.4,
                        component_locus_to_partition,
                    );
                    let head_partition = crate::owner_universe_partition(
                        head_seg.4,
                        component_locus_to_partition,
                    );
                    let outcome = span.restricted_charge(
                        &seam,
                        &[(
                            junction::junction_range(
                                prev_seg.0,
                                prev_seg.1,
                                prev_seg.2,
                                prev_seg.3,
                                prev_partition,
                            ),
                            junction::junction_range(
                                head_seg.0,
                                head_seg.1,
                                head_seg.2,
                                head_seg.3,
                                head_partition,
                            ),
                        )],
                        sources,
                        path_of_source,
                        &[prev_partition, head_partition],
                        model,
                    )?;
                    truth_boundary_charge[locus] += outcome.charge;
                    truth_boundary_restricted_total += outcome.charge;
                    truth_boundary_census[locus].push(serde_json::json!({
                        "kind": "boundary",
                        "left": [prev_seg.0, prev_seg.1, prev_seg.2, prev_seg.3],
                        "right": [head_seg.0, head_seg.1, head_seg.2, head_seg.3],
                        "spanning_reads": outcome.spanning_reads,
                        "events": outcome.events,
                        "span_features": outcome.span_feature_count,
                        "restricted_charge": outcome.charge,
                    }));
                }
            }
        }
        previous_tail = tail;
        previous_piece = pieces.last().copied();
        truth_profiles.push(profile);
    }
    if !missing.is_empty() {
        backgrounds.scan_extend(panel, missing, k, model)?;
    }
    // Pass 2 — charge everything against the completed backgrounds and
    // emit the per-feature table.
    // Ruling 2 (the chain-level once-per-material omission): the chains'
    // explained-feature sets — the union of the charged per-window profiles'
    // features (exactly the definition the scorer's chain-level omission
    // uses; the decompose's profiles ARE the scorer's charged profiles by
    // construction, so the exoneration agrees bit-for-bit). A feature the
    // chain spells at ANY window is exonerated from every window's omission
    // addend; the per-feature rows below carry the exonerations (their
    // sample_* columns read 0 where exonerated).
    let truth_explained: std::collections::HashSet<FeatureKey> = truth_profiles
        .iter()
        .flat_map(|profile| profile.keys().cloned())
        .collect();
    let exotic_explained: std::collections::HashSet<FeatureKey> = exotic_profiles
        .iter()
        .flat_map(|profile| profile.keys().cloned())
        .collect();
    // The selected side's observed side: the CHARGED pieces' owners (the
    // scorer's charged-runs form; the charged-owners map narrows exactly
    // where the shared-once filter removed a repeated piece).
    for locus in 0..locus_count {
        obs_exotic[locus] = if exotic_charged_owners[locus].is_empty() {
            partition_obs[locus].clone()
        } else {
            site.site_map(exotic_charged_owners[locus].iter().map(|&owner| {
                crate::owner_universe_partition(owner, component_locus_to_partition)
            }))
        };
    }
    let variant = match multiplicity_variant() {
        MultiplicityVariant::GentleBeta => "gentle-beta",
        MultiplicityVariant::Attribution => "attribution",
    };
    let mut loci_json: Vec<serde_json::Value> = Vec::with_capacity(locus_count);
    let mut exotic_local_total = 0.0f64;
    let mut truth_local_total = 0.0f64;
    let mut truth_grand_total = 0.0f64;
    // The sample-side background — the CHARGED genome-wide cross-support
    // form (owner decision 2026-09-24): beta_f = base +
    // Sum_r m_r*(1 - 1/t_r_genome). The LEGACY component-scoped
    // complementary-share form is accumulated alongside for the old-vs-new
    // delta (measurement only).
    let mut sample_exotic_total = 0.0f64;
    let mut sample_truth_total = 0.0f64;
    let mut sample_acceptance_violations = 0u64;
    let mut sample_acceptance_violation_sum = 0.0f64;
    let mut legacy_exotic_total = 0.0f64;
    let mut legacy_truth_total = 0.0f64;
    let mut legacy_acceptance_violations = 0u64;
    let mut legacy_acceptance_violation_sum = 0.0f64;
    let mut unpredicted_count = 0u64;
    let mut unpredicted_share = 0.0f64;
    // The sample-side backgrounds handle for the omission charge's helper
    // totals (Fix 2) and the sample-side pooled boundary charges.
    let sample_bg = sample_cross_support
        .map(|cross| crate::SampleSideBackgrounds::new(cross, model.background));
    let mut sample_boundary_pooled_total = 0.0f64;
    let mut sample_side_exonerated_exotic = 0.0f64;
    let mut sample_side_exonerated_truth = 0.0f64;
    for locus in 0..locus_count {
        let obs = &partition_obs[locus];
        // The owner-resolved observed sides (see the Pass-1 preamble): the
        // exotic row charges its domain allele's owner, the truth row its
        // pieces' owners (Fix 1: the record-once aggregation). The CHARGED
        // sample-side columns use the WINDOW'S full observed profile
        // (obs_window, Fix 2's omission universe) for BOTH rows.
        let obs_exotic = &obs_exotic[locus];
        let obs_truth = &obs_truth[locus];
        let obs_window: &HashMap<FeatureKey, f64> = &window_obs[locus];
        let exotic = &exotic_profiles[locus];
        let truth = &truth_profiles[locus];
        // Charge the boundary co-occurring seam pooled (after the scan
        // extension), exactly `score_reference_m1`'s boundary branch —
        // owner-resolved: the boundary pieces' owning partitions' shares.
        if let Some(seam) = &truth_boundary_pooled[locus] {
            let mut boundary_map: HashMap<FeatureKey, f64> = HashMap::new();
            if let Some(prev_seg) = &boundary_left_piece[locus] {
                for (feature, share) in crate::owner_routed_obs(
                    routed_equal,
                    prev_seg.4,
                    component_locus_to_partition,
                ) {
                    *boundary_map.entry(feature.clone()).or_default() += *share;
                }
            }
            if let Some(head_seg) = &boundary_right_piece[locus] {
                for (feature, share) in crate::owner_routed_obs(
                    routed_equal,
                    head_seg.4,
                    component_locus_to_partition,
                ) {
                    *boundary_map.entry(feature.clone()).or_default() += *share;
                }
            }
            if boundary_map.is_empty() {
                boundary_map = partition_obs[locus].clone();
            }
            truth_boundary_charge[locus] += profile_loss_boundary_multiplicity(
                seam,
                &boundary_map,
                &crate::EMPTY_ROUTED_OBS,
                model,
                backgrounds,
            )?;
            // The sample-side pooled boundary charge (the boundary layer
            // keeps the baselined sample form — the omission charge applies
            // at the per-window candidate charges only): needed so the
            // sample-side truth grand reconciles with the reference
            // ladder's haploid truth rescore.
            if let Some(sample) = &sample_bg {
                let charge = profile_loss_boundary_sample(
                    seam,
                    &boundary_map,
                    &crate::EMPTY_ROUTED_OBS,
                    sample,
                    model,
                )?;
                truth_boundary_charge_sample[locus] += charge;
                sample_boundary_pooled_total += charge;
            }
        }
        // The SELECTED side's pooled cross-locus boundary seam, charged in
        // the scorer's exact two-map sample form (the model of record the
        // selected chain's chain-level M1 charges).
        if let (Some(seam), Some(sample)) = (
            &exotic_boundary_pooled[locus],
            sample_bg.as_ref(),
        ) {
            if let (Some(prev_seg), Some(head_seg)) = (
                &exotic_boundary_left_piece[locus],
                &exotic_boundary_right_piece[locus],
            ) {
                let obs_left = crate::owner_routed_obs(
                    routed_equal,
                    prev_seg.4,
                    component_locus_to_partition,
                );
                let obs_right = crate::owner_routed_obs(
                    routed_equal,
                    head_seg.4,
                    component_locus_to_partition,
                );
                let owners_equal = prev_seg.4 == head_seg.4;
                let charge = if owners_equal {
                    crate::profile_loss_boundary_sample(
                        seam,
                        obs_left,
                        &crate::EMPTY_ROUTED_OBS,
                        sample,
                        model,
                    )?
                } else {
                    crate::profile_loss_boundary_sample(
                        seam,
                        obs_left,
                        obs_right,
                        sample,
                        model,
                    )?
                };
                exotic_boundary_charge_sample[locus] += charge;
                exotic_boundary_pooled_total += charge;
            }
        }
        // Per-feature union walk over the UNION of the two rows' profiles
        // AND the window's full observed profile (Fix 2's omission universe:
        // a feature observed in the window but predicted by neither row is
        // now a charged row — each row's omission term). Key order is the
        // merged BTreeMap order (deterministic; the profiles and the union
        // are BTreeMaps).
        let mut union: BTreeMap<&FeatureKey, (u64, u64)> = BTreeMap::new();
        for (key, &q) in exotic.iter() {
            union.entry(key).or_insert((0, 0)).0 += q;
        }
        for (key, &q) in truth.iter() {
            union.entry(key).or_insert((0, 0)).1 += q;
        }
        for key in obs_window.keys() {
            union.entry(key).or_insert((0, 0));
        }
        let mut records: Vec<serde_json::Value> = Vec::new();
        let mut sum_exotic = 0.0f64;
        let mut sum_truth = 0.0f64;
        let mut sum_legacy_exotic = 0.0f64;
        let mut sum_legacy_truth = 0.0f64;
        let mut sample_row_exotic = 0.0f64;
        let mut sample_row_truth = 0.0f64;
        // Ruling 2 diagnostic: the per-window omission mass the chain-level
        // once-per-material form exonerates (the chain spells the feature
        // somewhere else along the route).
        let mut exonerated_exotic = 0.0f64;
        let mut exonerated_truth = 0.0f64;
        let mut stats = DecomposeCategoryStats::default();
        for (key_ref, &(q_exotic, q_truth)) in union.iter() {
            let key: &FeatureKey = key_ref;
            let observed = obs.get(key).copied().unwrap_or(0.0);
            // Owner-resolved per-side observed support: the exotic row's
            // term charges its owner's share, the truth row's term its
            // pieces' owners' union-sum (identical to `observed` at every
            // anchor-only locus).
            let observed_exotic = obs_exotic.get(key).copied().unwrap_or(0.0);
            let observed_truth = obs_truth.get(key).copied().unwrap_or(0.0);
            let entry = backgrounds.entry(key);
            let beta = entry
                .map(|entry| {
                    if entry.beta.is_finite() && entry.beta > 0.0 {
                        entry.beta
                    } else {
                        model.background
                    }
                })
                .unwrap_or(model.background);
            let paths = entry.map(|entry| entry.paths).unwrap_or(0);
            let s_exotic = decompose_signal(model, q_exotic)?;
            let s_truth = decompose_signal(model, q_truth)?;
            let term_exotic = if q_exotic > 0 {
                loss_fractional_entry(model, q_exotic, observed_exotic, entry)?
            } else {
                0.0
            };
            let term_truth = if q_truth > 0 {
                loss_fractional_entry(model, q_truth, observed_truth, entry)?
            } else {
                0.0
            };
            sum_exotic += term_exotic;
            sum_truth += term_truth;
            stats.record(q_exotic > 0, q_truth > 0, observed_truth, term_exotic, term_truth);
            // The sample-side background columns: the LEGACY component-
            // scoped complementary share (kept for the old-vs-new delta
            // distribution) and the CHARGED genome-wide cross-support form
            // (the haploid track's model of record — reconciles with the
            // run's own sample-side charges).
            let (
                pooled,
                complement_legacy,
                beta_legacy,
                legacy_exotic,
                legacy_truth,
                genome_cross,
                beta_sample,
                sample_exotic,
                sample_truth,
            ) = match (pooled_support, sample_cross_support) {
                (Some(pooled_map), Some(cross_map)) => {
                    let pooled = pooled_map.get(key).copied().unwrap_or(0.0);
                    let complement_legacy = (pooled - observed).max(0.0);
                    let beta_legacy = model.background + complement_legacy;
                    let legacy_term = |signal: f64, observed: f64| {
                        if signal > 0.0 {
                            signal - observed * (signal / beta_legacy).ln_1p()
                        } else {
                            0.0
                        }
                    };
                    let genome_cross = cross_map.get(key).copied().unwrap_or(0.0);
                    let beta_sample = model.background + genome_cross;
                    // THE OMISSION CHARGE (Fix 2): the row's per-feature
                    // charge is its CURRENT relative term when it predicts
                    // the feature (its own owner-resolved observed side —
                    // bit-identical to the established form) and the
                    // pure-background omission cost
                    // beta_f - C_f*ln(beta_f) — with C_f the WINDOW'S
                    // record-once routed share — when it does not but the
                    // window observes it. A row omitting an observed
                    // feature therefore pays for explaining its observed
                    // reads as pure background: ignoring observed DNA is no
                    // longer free.
                    let observed_window = obs_window.get(key).copied().unwrap_or(0.0);
                    // The same Poisson-complete omission term the charging
                    // helper uses (including the count's log-factorial).
                    let omission_term = beta_sample
                        - observed_window * beta_sample.ln()
                        + crate::ln_gamma_observation(observed_window);
                    let sample_exotic = if q_exotic > 0 {
                        s_exotic - observed_exotic * (s_exotic / beta_sample).ln_1p()
                    } else if observed_window > 0.0 && !exotic_explained.contains(key) {
                        omission_term
                    } else {
                        0.0
                    };
                    let sample_truth = if q_truth > 0 {
                        s_truth - observed_truth * (s_truth / beta_sample).ln_1p()
                    } else if observed_window > 0.0 && !truth_explained.contains(key) {
                        omission_term
                    } else {
                        0.0
                    };
                    if q_exotic == 0 && observed_window > 0.0 && exotic_explained.contains(key) {
                        exonerated_exotic += omission_term;
                    }
                    if q_truth == 0 && observed_window > 0.0 && truth_explained.contains(key) {
                        exonerated_truth += omission_term;
                    }
                    (
                        pooled,
                        complement_legacy,
                        beta_legacy,
                        legacy_term(s_exotic, observed_exotic),
                        legacy_term(s_truth, observed_truth),
                        genome_cross,
                        beta_sample,
                        sample_exotic,
                        sample_truth,
                    )
                }
                _ => (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            };
            sample_row_exotic += sample_exotic;
            sample_row_truth += sample_truth;
            legacy_exotic_total += legacy_exotic;
            legacy_truth_total += legacy_truth;
            sum_legacy_exotic += legacy_exotic;
            sum_legacy_truth += legacy_truth;
            // Step-B acceptance arithmetic (the regression guard, now under
            // the charged genome-wide form): a correctly-predicted observed
            // feature (the row predicts it and the sample realizes it) must
            // not net positive against the sample-side background. The
            // legacy counters track the same test under the old form.
            if q_truth > 0 && observed_truth > 0.0 && sample_truth > 0.0 {
                sample_acceptance_violations += 1;
                sample_acceptance_violation_sum += sample_truth;
            }
            if q_truth > 0 && observed_truth > 0.0 && legacy_truth > 0.0 {
                legacy_acceptance_violations += 1;
                legacy_acceptance_violation_sum += legacy_truth;
            }
            records.push(serde_json::json!([
                decompose_feature_digest(key),
                key.len(),
                q_exotic,
                q_truth,
                observed,
                observed_exotic,
                observed_truth,
                paths,
                beta,
                s_exotic,
                s_truth,
                term_exotic,
                term_truth,
                pooled,
                complement_legacy,
                beta_legacy,
                legacy_exotic,
                legacy_truth,
                genome_cross,
                beta_sample,
                sample_exotic,
                sample_truth,
            ]));
        }
        // Observed features predicted by neither row (they cost both rows
        // exactly 0 — the deletion-bias control column).
        let mut neither_count = 0u64;
        let mut neither_share = 0.0f64;
        for (key, &share) in obs.iter() {
            if !exotic.contains_key(key) && !truth.contains_key(key) {
                neither_count += 1;
                neither_share += share;
            }
        }
        unpredicted_count += neither_count;
        unpredicted_share += neither_share;
        // Domain-allele context of the exotic row (validation: the chain's
        // states are domain alleles).
        let domain_allele = ranges[locus].iter().position(|traversal| {
            traversal.segments.len() == exotic_segments[locus].len()
                && traversal.segments.iter().zip(exotic_segments[locus].iter()).all(
                    |(segment, selected)| {
                        segment.source == selected.source
                            && segment.start == selected.start
                            && segment.end == selected.end
                            && segment.reverse == selected.reverse
                    },
                )
        });
        let segments_json: Vec<serde_json::Value> = exotic_segments[locus]
            .iter()
            .map(|segment| {
                serde_json::json!({
                    "source": segment.source,
                    "start": segment.start,
                    "end": segment.end,
                    "reverse": segment.reverse,
                })
            })
            .collect();
        let pieces_json: Vec<serde_json::Value> = truth_pieces[locus][0]
            .iter()
            .map(|&(source, start, end, reverse, owner)| {
                serde_json::json!([source, start, end, reverse, owner])
            })
            .collect();
        // The charged sample-side TOTALS come from the same helper the
        // rescore uses (bit-reconciliation by construction); the per-feature
        // rows above are the exact per-feature terms (their walk-order sums
        // are reported as the row sums).
        let mut sample_locus_exotic = 0.0f64;
        let mut sample_locus_truth = 0.0f64;
        if let Some(sample) = &sample_bg {
            sample_locus_exotic = crate::merged_single_loss_sample_with_omission_except(
                exotic,
                obs_exotic,
                &exotic_explained,
                obs_window,
                sample,
                model,
            )?;
            sample_locus_truth = crate::merged_single_loss_sample_with_omission_except(
                truth,
                obs_truth,
                &truth_explained,
                obs_window,
                sample,
                model,
            )?;
            sample_exotic_total += sample_locus_exotic;
            sample_truth_total += sample_locus_truth;
            sample_side_exonerated_exotic += exonerated_exotic;
            sample_side_exonerated_truth += exonerated_truth;
        }
        let truth_local = sum_truth + truth_interior_charge[locus];
        exotic_local_total += sum_exotic;
        truth_local_total += truth_local;
        truth_grand_total += truth_local + truth_boundary_charge[locus];
        let is_tract = truth_pieces[locus]
            .iter()
            .flat_map(|pieces| pieces.iter().map(|&(source, _, _, _, _)| source))
            .any(|source| donor_sources.contains(&source));
        loci_json.push(serde_json::json!({
            "locus": locus,
            "full_locus": locus_offset + locus,
            "axis_interval": [axis_slice[locus].start, axis_slice[locus].end],
            "is_tract": is_tract,
            "exotic": {
                "segments": segments_json,
                "domain_allele": domain_allele,
                "class": domain_allele
                    .map(|allele| locus_classes[locus].membership[allele]),
                "features": exotic.len(),
                "sum_terms": sum_exotic,
                "interior_restricted_charge": exotic_interior_charge[locus],
                "interior_census": exotic_interior_census[locus],
                "boundary_charge": exotic_boundary_charge[locus],
                "boundary_census": exotic_boundary_census[locus],
                "boundary_pooled": exotic_boundary_charge_sample[locus],
                "grand_including_boundary": sample_locus_exotic
                    + exotic_interior_charge[locus]
                    + exotic_boundary_charge[locus]
                    + exotic_boundary_charge_sample[locus],
            },
            "truth": {
                "pieces": pieces_json,
                "features": truth.len(),
                "sum_terms": sum_truth,
                "interior_restricted_charge": truth_interior_charge[locus],
                "interior_census": truth_interior_census[locus],
                "boundary_charge": truth_boundary_charge[locus],
                "boundary_census": truth_boundary_census[locus],
                "local": truth_local,
            },
            "delta_truth_minus_exotic_local": truth_local - sum_exotic,
            "sample_side": {
                // The CHARGED genome-wide cross-support locals under the
                // omission charge (Fix 2: unbaselined over the window's full
                // observed profile — helper totals, bit-reconciling with the
                // rescore), the per-feature row-walk sums, and the LEGACY
                // component-scoped complementary-share locals, per locus.
                "exotic_local": sample_locus_exotic,
                "truth_local": sample_locus_truth,
                "delta_truth_minus_exotic": sample_locus_truth - sample_locus_exotic,
                "row_walk_sums": {
                    "exotic": sample_row_exotic,
                    "truth": sample_row_truth,
                },
                "boundary_pooled": truth_boundary_charge_sample[locus],
                "omission_exonerated": {
                    "exotic": exonerated_exotic,
                    "truth": exonerated_truth,
                },
                "legacy": {
                    "exotic_local": sum_legacy_exotic,
                    "truth_local": sum_legacy_truth,
                    "delta_truth_minus_exotic": sum_legacy_truth - sum_legacy_exotic,
                },
            },
            "categories": stats.json(),
            "unpredicted_observed": {
                "count": neither_count,
                "sum_share": neither_share,
            },
            // [digest, feature_token_len (odd; anchors = (len+1)/2),
            //  q_exotic, q_truth, C (axis view), observed_exotic (owner),
            //  observed_truth (owner, record-once), paths, beta,
            //  s_exotic, s_truth, term_exotic, term_truth, pooled,
            //  complement_legacy, beta_legacy, legacy_term_exotic,
            //  legacy_term_truth, genome_cross_support, beta_sample,
            //  sample_term_exotic, sample_term_truth (both charged over the
            //  WINDOW'S record-once observed share, unbaselined — Fix 2)]
            "features": records,
        }));
    }
    let out = serde_json::json!({
        "what": "per-feature decomposition of the haploid single-allele local \
                 table: selected chain slot-0 rows vs the truth mosaic copy-0 \
                 pieces, per locus (rescore convention: routed shares, measured \
                 backgrounds, oracle profiles; the sample-side columns charge \
                 the omission charge — Fix 2 — over the window's full \
                 record-once observed profile, and the owner-resolved truth \
                 column aggregates records once — Fix 1)",
        "variant": variant,
        "containment_probe": contain_probe,
        "model": {
            "base": model.background,
            "per_count": model.histogram as f64 * model.depth / model.denominator,
            "depth": model.depth,
            "histogram": model.histogram,
            "denominator": model.denominator,
        },
        "loci": loci_json,
        "totals": {
            "exotic_local": exotic_local_total,
            "truth_local": truth_local_total,
            "truth_grand_including_boundary": truth_grand_total,
            "truth_boundary_total": truth_grand_total - truth_local_total,
            "sample_side_exotic_local": sample_exotic_total,
            "sample_side_truth_local": sample_truth_total,
            "sample_side_margin_truth_minus_exotic": sample_truth_total - sample_exotic_total,
            "sample_side_boundary_pooled_total": sample_boundary_pooled_total,
            "sample_side_exotic_grand_including_boundary": sample_exotic_total
                + exotic_interior_total
                + exotic_boundary_restricted_total
                + exotic_boundary_pooled_total,
            "sample_side_omission_exonerated_exotic": sample_side_exonerated_exotic,
            "sample_side_omission_exonerated_truth": sample_side_exonerated_truth,
            "sample_side_truth_grand_including_boundary": sample_truth_total
                + sample_boundary_pooled_total
                + truth_interior_total
                + truth_boundary_restricted_total,
            "sample_side_acceptance_violations": sample_acceptance_violations,
            "sample_side_acceptance_violation_sum": sample_acceptance_violation_sum,
            "sample_side_legacy_exotic_local": legacy_exotic_total,
            "sample_side_legacy_truth_local": legacy_truth_total,
            "sample_side_legacy_margin_truth_minus_exotic": legacy_truth_total
                - legacy_exotic_total,
            "sample_side_legacy_acceptance_violations": legacy_acceptance_violations,
            "sample_side_legacy_acceptance_violation_sum": legacy_acceptance_violation_sum,
            "unpredicted_observed_count": unpredicted_count,
            "unpredicted_observed_share": unpredicted_share,
        },
        "sample_side_background": {
            "definition": "beta_f = base + Sum_r m_r*(1 - 1/t_r_genome) over \
                            every multiset record r realizing f; m_r = the \
                            record's genome-wide multiplicity, t_r_genome = its \
                            touched count across ALL components' partitions \
                            (the genome-universe placement pass). The observed \
                            shares stay component-scoped. Owner-decided \
                            candidate family 1 (the genome-wide touched set), \
                            2026-09-24.",
            "charged_here": "the haploid track's sample-side charges under \
                             the omission charge (Fix 2, owner-approved \
                             2026-09-24): every haploid-track candidate is \
                             charged over the WINDOW'S full observed profile \
                             (the record-once routed shares over the window's \
                             owning partitions), unbaselined: term = \
                             (s_f + beta_f) - C_f*ln(s_f + beta_f); an omitted \
                             observed feature contributes beta_f - \
                             C_f*ln(beta_f). Boundary seams keep the \
                             established baselined boundary convention. The \
                             legacy component-scoped complementary-share \
                             columns are kept for the old-vs-new delta \
                             distribution",
        },
    });
    let bytes = serde_json::to_vec(&out).map_err(io::Error::other)?;
    std::fs::write(out_path, &bytes)?;
    Ok(serde_json::json!({
        "out_path": out_path.display().to_string(),
        "variant": variant,
        "sample_side_truth_grand_including_boundary": sample_truth_total
            + sample_boundary_pooled_total
            + truth_interior_total
            + truth_boundary_restricted_total,
        "exotic_local_total": exotic_local_total,
        "truth_local_total": truth_local_total,
        "truth_grand_total": truth_grand_total,
        "per_locus": (0..locus_count)
            .map(|locus| {
                let exotic = out["loci"][locus]["exotic"]["sum_terms"].as_f64().unwrap_or(0.0);
                let truth = out["loci"][locus]["truth"]["local"].as_f64().unwrap_or(0.0);
                serde_json::json!({
                    "locus": locus,
                    "full_locus": locus_offset + locus,
                    "is_tract": out["loci"][locus]["is_tract"].as_bool().unwrap_or(false),
                    "exotic_local": exotic,
                    "truth_local": truth,
                    "delta": truth - exotic,
                })
            })
            .collect::<Vec<_>>(),
        "wall_seconds": started.elapsed().as_secs_f64(),
    }))
}

/// Per-locus category aggregates of the decomposition table (which row
/// predicts the feature, summed observed share and loss contributions).
#[derive(Default)]
struct DecomposeCategoryStats {
    both: DecomposeCategoryAggregate,
    exotic_only: DecomposeCategoryAggregate,
    truth_only: DecomposeCategoryAggregate,
}

#[derive(Default)]
struct DecomposeCategoryAggregate {
    count: u64,
    sum_share: f64,
    sum_term_exotic: f64,
    sum_term_truth: f64,
    negative_truth: u64,
    positive_truth: u64,
}

impl DecomposeCategoryStats {
    fn record(
        &mut self,
        predicted_exotic: bool,
        predicted_truth: bool,
        share: f64,
        term_exotic: f64,
        term_truth: f64,
    ) {
        let aggregate = match (predicted_exotic, predicted_truth) {
            (true, true) => &mut self.both,
            (true, false) => &mut self.exotic_only,
            (false, true) => &mut self.truth_only,
            (false, false) => return,
        };
        aggregate.count += 1;
        aggregate.sum_share += share;
        aggregate.sum_term_exotic += term_exotic;
        aggregate.sum_term_truth += term_truth;
        if term_truth < 0.0 {
            aggregate.negative_truth += 1;
        } else if term_truth > 0.0 {
            aggregate.positive_truth += 1;
        }
    }

    fn json(&self) -> serde_json::Value {
        let encode = |aggregate: &DecomposeCategoryAggregate| {
            serde_json::json!({
                "count": aggregate.count,
                "sum_share": aggregate.sum_share,
                "sum_term_exotic": aggregate.sum_term_exotic,
                "sum_term_truth": aggregate.sum_term_truth,
                "truth_terms_negative": aggregate.negative_truth,
                "truth_terms_positive": aggregate.positive_truth,
            })
        };
        serde_json::json!({
            "both": encode(&self.both),
            "exotic_only": encode(&self.exotic_only),
            "truth_only": encode(&self.truth_only),
        })
    }
}



// ---------------------------------------------------------------------------
// Stage 3 detector: the natural gap of a residual distribution.
// ---------------------------------------------------------------------------

/// Parameter-free demand rule: among the loci's positive statistic values,
/// the largest ratio gap between consecutive sorted values is the
/// distribution's natural break; loci above it are selected. A largest gap
/// of exactly 1 (no gap) selects nothing. The full sorted distribution and
/// every gap are reported with the selection — the rule is a reported
/// statistic, not a tuned constant.
pub(super) fn natural_gap_selection(values: &[(usize, f64)]) -> (Vec<usize>, serde_json::Value) {
    let mut sorted: Vec<(usize, f64)> =
        values.iter().copied().filter(|&(_, v)| v > 0.0).collect();
    sorted.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let mut best: Option<(usize, f64)> = None;
    for i in 0..sorted.len().saturating_sub(1) {
        let ratio = sorted[i + 1].1 / sorted[i].1;
        if best.is_none_or(|(_, best_ratio)| ratio > best_ratio) {
            best = Some((i, ratio));
        }
    }
    let selected: Vec<usize> = match best {
        Some((i, ratio)) if ratio > 1.0 => {
            sorted[i + 1..].iter().map(|&(locus, _)| locus).collect()
        }
        _ => Vec::new(),
    };
    let rows: Vec<serde_json::Value> = sorted
        .windows(2)
        .enumerate()
        .map(|(i, window)| {
            serde_json::json!({
                "locus": window[0].0,
                "value": window[0].1,
                "next_locus": window[1].0,
                "next_value": window[1].1,
                "ratio_gap": window[1].1 / window[0].1,
                "_index": i,
            })
        })
        .collect();
    let largest_gap = best.map(|(i, ratio)| serde_json::json!({
        "below_locus": sorted[i].0,
        "below_value": sorted[i].1,
        "above_locus": sorted.get(i + 1).map(|&(l, _)| l),
        "above_value": sorted.get(i + 1).map(|&(_, v)| v),
        "ratio": ratio,
    }));
    (
        selected,
        serde_json::json!({
            "sorted": sorted
                .iter()
                .map(|&(locus, value)| serde_json::json!({"locus": locus, "value": value}))
                .collect::<Vec<_>>(),
            "gaps": rows,
            "largest_ratio_gap": largest_gap,
        }),
    )
}

// ---------------------------------------------------------------------------
// Stage 2: exact chain DP.
// ---------------------------------------------------------------------------

pub(super) struct SpineState {
    pub(super) score: f64,
    pub(super) pair: [usize; 2],
    pub(super) pred: Option<u32>,
    /// Windowed conflict history INCLUDING this pair (interned exactly).
    pub(super) history_id: u32,
}

pub(super) struct SpineDpOutcome {
    pub(super) route: Vec<[usize; 2]>,
    pub(super) best_score: f64,
    /// Terminal states whose score is bit-identical to the best (the exact
    /// tie set, retained as per-locus pair sets below).
    pub(super) tie_count: usize,
    pub(super) tie_pair_sets: Vec<Vec<[usize; 2]>>,
    pub(super) state_counts: Vec<usize>,
    pub(super) retained_pairs: Vec<u64>,
    pub(super) transitions: Vec<u64>,
    pub(super) merged_drops: Vec<u64>,
    pub(super) bound_pruned: Vec<u64>,
    pub(super) margin_pruned: Vec<u64>,
}

/// The EXACT min-plus chain DP. No width, no tie rotations, no seed ladder:
/// every pair state within the admissible margin is kept; per-state merging
/// is exact on (pair, windowed history); the only prune is the
/// exactness-preserving admissible check against the incumbent backbone
/// total (any pruned edge provably sits in no chain at or below the
/// incumbent, and the optimum is at or below it).
#[allow(clippy::too_many_arguments)]
pub(super) fn run_exact_chain_dp(
    ranges: &[Vec<genome::SpanningTraversal>],
    membership: &[Vec<usize>],
    sweeps: &[LocusSweep],
    successors: &[Vec<Vec<(usize, f64)>>],
    viable: &[Vec<bool>],
    margins: &[f64],
    suffix_from_next: &[f64],
    incumbent: f64,
    conflict_window: usize,
    backbone: &[usize],
    rss: &mut genome::PeriodicRssGuard,
) -> io::Result<SpineDpOutcome> {
    let locus_count = ranges.len();
    ensure(locus_count > 0, "empty spine DP domain")?;
    ensure(backbone.len() == locus_count, "backbone cardinality mismatch")?;
    // The all-native backbone pair is structurally retained at every locus
    // (the incumbent chain that anchors admissibility and guarantees a
    // completable chain exists); every other pair is retained iff its loss
    // is within the local seam-swing margin of the locus's viable best.
    let backbone_class = |locus: usize| membership[locus][backbone[locus]];
    let class_pair_retained = |locus: usize, c1: usize, c2: usize, loss: f64| -> bool {
        let cb = backbone_class(locus);
        loss <= margins[locus] || (c1 == cb && c2 == cb)
    };
    let mut history_store: Vec<Vec<[usize; 2]>> = vec![Vec::new()];
    let mut history_map: HashMap<Vec<[usize; 2]>, u32> = HashMap::new();
    history_map.insert(Vec::new(), 0);

    // Viable members per class (per locus).
    let viable_members: Vec<Vec<Vec<usize>>> = (0..locus_count)
        .map(|locus| {
            let classes = membership[locus]
                .iter()
                .copied()
                .max()
                .map(|class| class + 1)
                .unwrap_or(0);
            let mut members = vec![Vec::new(); classes];
            for (allele, &class) in membership[locus].iter().enumerate() {
                if viable[locus][allele] {
                    members[class].push(allele);
                }
            }
            members
        })
        .collect();

    // Retained physical pair counts per locus (margin statistics; the
    // structurally-retained backbone pair is counted once when it sits
    // outside the margin).
    let retained_pairs: Vec<u64> = (0..locus_count)
        .map(|locus| {
            let classes = viable_members[locus].len();
            let mut count = 0u64;
            for c1 in 0..classes {
                if viable_members[locus][c1].is_empty() {
                    continue;
                }
                for c2 in c1..classes {
                    if viable_members[locus][c2].is_empty() {
                        continue;
                    }
                    if !class_pair_retained(
                        locus,
                        c1,
                        c2,
                        sweeps[locus].table[class_pair_index(c1, c2)],
                    ) {
                        continue;
                    }
                    let (m1, m2) = (
                        viable_members[locus][c1].len() as u64,
                        viable_members[locus][c2].len() as u64,
                    );
                    count += if c1 == c2 { m1 * (m1 + 1) / 2 } else { m1 * m2 };
                }
            }
            count
        })
        .collect();

    let mut layers: Vec<Vec<SpineState>> = Vec::with_capacity(locus_count);
    let mut state_counts: Vec<usize> = Vec::with_capacity(locus_count);
    let mut transitions: Vec<u64> = Vec::with_capacity(locus_count);
    let mut merged_drops: Vec<u64> = Vec::with_capacity(locus_count);
    let mut bound_pruned: Vec<u64> = Vec::with_capacity(locus_count);
    let mut margin_pruned: Vec<u64> = Vec::with_capacity(locus_count);

    // ---- initial layer (locus 0): every retained viable pair.
    {
        let classes = viable_members[0].len();
        let mut layer: Vec<SpineState> = Vec::new();
        let mut index: HashMap<([usize; 2], u32), usize> = HashMap::new();
        let mut pruned_bound = 0u64;
        let mut pruned_margin = 0u64;
        for c1 in 0..classes {
            if viable_members[0][c1].is_empty() {
                continue;
            }
            for c2 in c1..classes {
                if viable_members[0][c2].is_empty() {
                    continue;
                }
                let loss = sweeps[0].table[class_pair_index(c1, c2)];
                let (m1, m2) = (
                    viable_members[0][c1].len() as u64,
                    viable_members[0][c2].len() as u64,
                );
                let feasible = if c1 == c2 { m1 * (m1 + 1) / 2 } else { m1 * m2 };
                if !class_pair_retained(0, c1, c2, loss) {
                    pruned_margin += feasible;
                    continue;
                }
                if loss + suffix_from_next[0] > incumbent + BOUND_PRUNE_EPSILON {
                    pruned_bound += feasible;
                    continue;
                }
                for &a in &viable_members[0][c1] {
                    for &b in &viable_members[0][c2] {
                        let pair = [a.min(b), a.max(b)];
                        let content = vec![pair];
                        let history_id =
                            intern_history(&content, &mut history_store, &mut history_map);
                        match index.get(&(pair, history_id)) {
                            Some(&slot) => {
                                if layer[slot].score > loss {
                                    layer[slot].score = loss;
                                }
                            }
                            None => {
                                index.insert((pair, history_id), layer.len());
                                layer.push(SpineState {
                                    score: loss,
                                    pair,
                                    pred: None,
                                    history_id,
                                });
                            }
                        }
                    }
                }
            }
        }
        ensure(!layer.is_empty(), "spine DP initial layer is empty")?;
        state_counts.push(layer.len());
        transitions.push(layer.len() as u64);
        merged_drops.push(0);
        bound_pruned.push(pruned_bound);
        margin_pruned.push(pruned_margin);
        layers.push(layer);
    }

    // ---- transitions.
    for locus in 1..locus_count {
        let previous = layers.last().expect("layer");
        let next_count = ranges[locus].len();
        let mut spans_by_source: BTreeMap<usize, Vec<(u64, u64, usize)>> = BTreeMap::new();
        for (member, traversal) in ranges[locus].iter().enumerate() {
            for segment in &traversal.segments {
                if segment.start < segment.end {
                    spans_by_source
                        .entry(segment.source)
                        .or_default()
                        .push((segment.start, segment.end, member));
                }
            }
        }
        for spans in spans_by_source.values_mut() {
            spans.sort();
        }
        // Folded successor groups per PREVIOUS-locus left allele: members
        // grouped by (next class, seam score) — every member of a group pair
        // contributes the same local constant and seam terms, so the exact
        // DP computes the transition once per group pair and enumerates the
        // physical member pairs (both homolog matchings arise by iterating
        // all group-pair combinations, exactly as in the quarantined beam).
        let grouped: Vec<Vec<(usize, f64, Vec<usize>)>> = successors[locus - 1]
            .iter()
            .map(|list| {
                let mut groups: Vec<(usize, f64, Vec<usize>)> = Vec::new();
                for &(next, seam) in list {
                    if !viable[locus][next] {
                        continue;
                    }
                    let class = membership[locus][next];
                    match groups.iter_mut().find(|group| {
                        group.0 == class && group.1.to_bits() == seam.to_bits()
                    }) {
                        Some(group) => group.2.push(next),
                        None => groups.push((class, seam, vec![next])),
                    }
                }
                groups.sort_by(|a, b| {
                    a.0.cmp(&b.0).then_with(|| a.1.to_bits().cmp(&b.1.to_bits()))
                });
                for group in &mut groups {
                    group.2.sort_unstable();
                }
                groups
            })
            .collect();
        let mut layer: Vec<SpineState> = Vec::new();
        let mut index: HashMap<([usize; 2], u32), usize> = HashMap::new();
        let mut transitions_locus = 0u64;
        let mut merged_locus = 0u64;
        let mut pruned_bound = 0u64;
        let mut pruned_margin = 0u64;
        for (state_index, state) in previous.iter().enumerate() {
            rss.checkpoint("spine_dp_hot_loop")?;
            let history: Vec<[usize; 2]> = history_store[state.history_id as usize].clone();
            let first_locus = locus - history.len();
            let mut eligible = [vec![true; next_count], vec![true; next_count]];
            for (offset, previous_pair) in history.iter().enumerate() {                let previous_locus = first_locus + offset;
                for copy in 0..2 {
                    for segment in &ranges[previous_locus][previous_pair[copy]].segments {
                        let Some(spans) = spans_by_source.get(&segment.source) else {
                            continue;
                        };
                        let prefix =
                            spans.partition_point(|&(start, _, _)| start < segment.end);
                        for &(_start, end, member) in &spans[..prefix] {
                            if !eligible[copy][member] || end <= segment.start {
                                continue;
                            }
                            if ranges[locus][member].segments.iter().any(|right| {
                                segment.source == right.source
                                    && segment.start < right.end
                                    && right.start < segment.end
                            }) {
                                eligible[copy][member] = false;
                            }
                        }
                    }
                }
            }
            let member_stats = |copy: usize, members: &[usize]| -> (u64, usize) {
                let mut count = 0u64;
                let mut first = usize::MAX;
                for &member in members {
                    if eligible[copy][member] {
                        count += 1;
                        if first == usize::MAX {
                            first = member;
                        }
                    }
                }
                (count, first)
            };
            let first_groups = &grouped[state.pair[0]];
            let second_groups = &grouped[state.pair[1]];
            for (class_first, first_seam, first_members) in first_groups {
                let (first_count, _) = member_stats(0, first_members);
                if first_count == 0 {
                    continue;
                }
                for (class_second, second_seam, second_members) in second_groups {
                    let (second_count, _) = member_stats(1, second_members);
                    let feasible = first_count.saturating_mul(second_count);
                    if feasible == 0 {
                        continue;
                    }
                    let pair_loss =
                        sweeps[locus].table[class_pair_index(*class_first, *class_second)];
                    if !class_pair_retained(locus, *class_first, *class_second, pair_loss) {
                        pruned_margin += feasible;
                        continue;
                    }
                    let score = state.score + first_seam + second_seam + pair_loss;
                    if score + suffix_from_next[locus] > incumbent + BOUND_PRUNE_EPSILON {
                        pruned_bound += feasible;
                        continue;
                    }
                    for &first in first_members {
                        if !eligible[0][first] {
                            continue;
                        }
                        for &second in second_members {
                            if !eligible[1][second] {
                                continue;
                            }
                            transitions_locus += 1;
                            let pair = [first.min(second), first.max(second)];
                            let keep =
                                history.len().saturating_sub(conflict_window.saturating_sub(1));
                            let mut content: Vec<[usize; 2]> = Vec::with_capacity(conflict_window);
                            content.extend_from_slice(&history[keep..]);
                            content.push(pair);
                            let history_id = intern_history(
                                &content,
                                &mut history_store,
                                &mut history_map,
                            );
                            match index.get(&(pair, history_id)) {
                                Some(&slot) => {
                                    if layer[slot].score > score {
                                        layer[slot].score = score;
                                        layer[slot].pred = Some(state_index as u32);
                                    } else {
                                        merged_locus += 1;
                                    }
                                }
                                None => {
                                    index.insert((pair, history_id), layer.len());
                                    layer.push(SpineState {
                                        score,
                                        pair,
                                        pred: Some(state_index as u32),
                                        history_id,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
        ensure(!layer.is_empty(), "spine DP layer {locus} is empty")?;
        eprintln!(
            "[spine] dp layer {locus}: states {} transitions {} merged {} bound_pruned {} margin_pruned {}",
            layer.len(),
            transitions_locus,
            merged_locus,
            pruned_bound,
            pruned_margin
        );
        state_counts.push(layer.len());
        transitions.push(transitions_locus);
        merged_drops.push(merged_locus);
        bound_pruned.push(pruned_bound);
        margin_pruned.push(pruned_margin);
        layers.push(layer);
    }

    // ---- terminal best, exact tie set, representative backtracking.
    let final_layer = layers.last().expect("final layer");
    let best = final_layer
        .iter()
        .map(|state| state.score)
        .fold(f64::INFINITY, f64::min);
    ensure(best.is_finite(), "spine DP terminal layer has no finite state")?;
    let tied: Vec<usize> = final_layer
        .iter()
        .enumerate()
        .filter(|&(_, state)| state.score.to_bits() == best.to_bits())
        .map(|(index, _)| index)
        .collect();
    ensure(!tied.is_empty(), "spine DP tie set is empty")?;
    let backtrack = |start: usize| -> Vec<[usize; 2]> {
        let mut route = Vec::with_capacity(locus_count);
        let mut cursor = (start, locus_count - 1);
        loop {
            let (state_index, locus) = cursor;
            let state = &layers[locus][state_index];
            route.push(state.pair);
            match state.pred {
                Some(previous) if locus > 0 => cursor = (previous as usize, locus - 1),
                _ => break,
            }
        }
        route.reverse();
        route
    };
    let route = backtrack(tied[0]);
    let mut pair_sets: Vec<HashSet<[usize; 2]>> = vec![HashSet::new(); locus_count];
    for &terminal in &tied {
        for (locus, pair) in backtrack(terminal).iter().enumerate() {
            pair_sets[locus].insert(*pair);
        }
    }
    let tie_pair_sets: Vec<Vec<[usize; 2]>> = pair_sets
        .into_iter()
        .map(|mut set| {
            let mut pairs: Vec<[usize; 2]> = set.drain().collect();
            pairs.sort_unstable();
            pairs
        })
        .collect();
    Ok(SpineDpOutcome {
        route,
        best_score: best,
        tie_count: tied.len(),
        tie_pair_sets,
        state_counts,
        retained_pairs,
        transitions,
        merged_drops,
        bound_pruned,
        margin_pruned,
    })
}

// ---------------------------------------------------------------------------
// Boundary construction (shared by the initial pass and Stage-3 rebuilds).
// ---------------------------------------------------------------------------

pub(super) struct SpineBoundary {
    /// Per left allele: (right allele, seam loss) legal links.
    pub(super) successors: Vec<Vec<(usize, f64)>>,
    pub(super) seam_by_pair: HashMap<(usize, usize), std::sync::Arc<Profile>>,
    pub(super) compositions: usize,
}

/// Draft of one boundary's seam compositions: profiles, realizing pairs,
/// pooled/novel kinds, and the NOVEL-only restricted charges (the junction
/// machinery — flat background, unchanged). The POOLED composition losses
/// are deferred to `finalize_spine_boundary`: their features must be in
/// the multiplicity-background scan before they are charged, so the draft
/// is built first, its features join the scan universe, and the finalize
/// step charges them under the per-feature backgrounds.
pub(super) struct SpineBoundaryDraft {
    pub(super) profiles: Vec<std::sync::Arc<Profile>>,
    pub(super) composition_pooled: Vec<bool>,
    /// Restricted charge + spanning-read count per composition (valid where
    /// !composition_pooled).
    pub(super) novel_scores: Vec<(f64, u64)>,
    pub(super) link_composition: Vec<u32>,
    pub(super) links: Vec<(usize, usize)>,
    pub(super) pooled_count: usize,
    pub(super) novel_spanning: usize,
}

/// Endpoint-class boundary construction, phase 1: L149 flanks per allele
/// (orientation-true), deduplicated seam compositions over the legal links,
/// their profiles, and the novel-only restricted charges.
#[allow(clippy::too_many_arguments)]
pub(super) fn build_spine_boundary_draft(
    panel: &SyngIndex,
    sources: &routes::Sources,
    flank_memo: &FlankMemo,
    model: &ScoreModel,
    span: &JunctionSpanIndex,
    path_of_source: &[usize],
    // Per-allele OWNING universe partitions (the window-domain extension's
    // owner-resolved charging; the pre-extension model's every allele is
    // owned by its locus's axis partition).
    left_owners: &[u32],
    right_owners: &[u32],
    routed_equal: &crate::RoutedObs,
    component_locus_to_partition: &[u32],
    left_ranges: &[genome::SpanningTraversal],
    right_ranges: &[genome::SpanningTraversal],
    links: &[(usize, usize)],
) -> io::Result<SpineBoundaryDraft> {
    let left_ends: Vec<Vec<u8>> = left_ranges
        .par_iter()
        .map(|traversal| {
            allele_endpoints(sources, flank_memo, &traversal.segments, 149).map(|(_, tail)| tail)
        })
        .collect::<io::Result<_>>()?;
    let right_starts: Vec<Vec<u8>> = right_ranges
        .par_iter()
        .map(|traversal| {
            allele_endpoints(sources, flank_memo, &traversal.segments, 149).map(|(head, _)| head)
        })
        .collect::<io::Result<_>>()?;
    // Seam compositions over the legal links (the seam machinery's own
    // junction granularity), each with its realizing (exit, entry) segment
    // pairs and its kind: a composition is POOLED when any realizing link
    // is a co-occurring real panel junction (those keep the pooled charge
    // bit-identically); novel-only compositions are paid by the reads that
    // actually cross them.
    let mut composition_ids: HashMap<(Vec<u8>, Vec<u8>), u32> = HashMap::new();
    let mut composition_keys: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut composition_pairs: Vec<Vec<(SourceRange, SourceRange)>> = Vec::new();
    let mut composition_pooled: Vec<bool> = Vec::new();
    // Per composition: the DISTINCT owning universe partitions of its
    // realizing links (the crossing reads' share targets; the pre-extension
    // model's compositions realize only the two locus partitions).
    let mut composition_owners: Vec<BTreeSet<u32>> = Vec::new();
    let mut link_composition: Vec<u32> = Vec::with_capacity(links.len());
    for &(left, right) in links {
        let key = (left_ends[left].clone(), right_starts[right].clone());
        let id = match composition_ids.get(&key) {
            Some(&id) => id,
            None => {
                let id = composition_keys.len() as u32;
                composition_ids.insert(key.clone(), id);
                composition_keys.push(key);
                composition_pairs.push(Vec::new());
                composition_pooled.push(false);
                composition_owners.push(BTreeSet::new());
                id
            }
        };
        composition_pairs[id as usize].push((
            left_ranges[left]
                .segments
                .last()
                .expect("nonempty traversal")
                .clone(),
            right_ranges[right]
                .segments
                .first()
                .expect("nonempty traversal")
                .clone(),
        ));
        composition_owners[id as usize].insert(left_owners[left]);
        composition_owners[id as usize].insert(right_owners[right]);
        if phasing::cooccurrence_gap(&left_ranges[left], &right_ranges[right])
            .is_some_and(|gap| gap >= 0)
        {
            composition_pooled[id as usize] = true;
        }
        link_composition.push(id);
    }
    let compositions = composition_keys.len();
    let profiles: Vec<std::sync::Arc<Profile>> = composition_keys
        .par_iter()
        .map(|(left_end, right_start)| {
            let (profile, _) =
                genome::profile_event_seam(panel, left_end, right_start, READ_LENGTH, MAX_FEATURES)?;
            Ok(std::sync::Arc::new(profile))
        })
        .collect::<io::Result<_>>()?;
    // Novel-only restricted charges (the junction machinery, unchanged);
    // pooled compositions carry a placeholder until the finalize step.
    let novel_scores: Vec<(f64, u64)> = (0..compositions)
        .into_par_iter()
        .map(|id| {
            if composition_pooled[id] {
                Ok((0.0, 0))
            } else {
                let partitions = composition_owners[id]
                    .iter()
                    .map(|&owner| {
                        crate::owner_universe_partition(owner, component_locus_to_partition)
                    })
                    .collect::<Vec<_>>();
                let outcome = span.restricted_charge(
                    &profiles[id],
                    &composition_pairs[id],
                    sources,
                    path_of_source,
                    &partitions,
                    model,
                )?;
                Ok((outcome.charge, outcome.spanning_reads))
            }
        })
        .collect::<io::Result<_>>()?;
    let pooled_count = composition_pooled.iter().filter(|&&pooled| pooled).count();
    let novel_spanning = novel_scores
        .iter()
        .zip(composition_pooled.iter())
        .filter(|(&(_, spanning), &pooled)| !pooled && spanning > 0)
        .count();
    Ok(SpineBoundaryDraft {
        profiles,
        composition_pooled,
        novel_scores,
        link_composition,
        links: links.to_vec(),
        pooled_count,
        novel_spanning,
    })
}

/// Endpoint-class boundary construction, phase 2: the POOLED composition
/// losses under the per-feature multiplicity backgrounds, and the successor
/// assembly (a composition realized by any co-occurring pair keeps the
/// pooled charge; novel-only compositions keep their restricted charge).
pub(super) fn finalize_spine_boundary(
    draft: SpineBoundaryDraft,
    left_ranges: &[genome::SpanningTraversal],
    right_ranges: &[genome::SpanningTraversal],
    left_owners: &[u32],
    right_owners: &[u32],
    routed_equal: &crate::RoutedObs,
    component_locus_to_partition: &[u32],
    model: &ScoreModel,
    backgrounds: &FeatureBackgrounds,
) -> io::Result<SpineBoundary> {
    let SpineBoundaryDraft {
        profiles,
        composition_pooled,
        novel_scores,
        link_composition,
        links,
        pooled_count,
        novel_spanning,
    } = draft;
    let compositions = profiles.len();
    let mut composition_scores = novel_scores;
    for id in 0..compositions {
        if composition_pooled[id] {
            // The pooled seam charge's observed sides: the SUM over the
            // DISTINCT owning partitions of the composition's realizing
            // links (the pre-extension model's two locus partitions are the
            // universal single case, reproduced bit-for-bit).
            let mut owners: BTreeSet<u32> = BTreeSet::new();
            for (&(left, right), &composition) in
                links.iter().zip(link_composition.iter())
            {
                if composition as usize == id {
                    owners.insert(left_owners[left]);
                    owners.insert(right_owners[right]);
                }
            }
            let mut merged: HashMap<FeatureKey, f64> = HashMap::new();
            for owner in owners {
                let map = crate::owner_routed_obs(
                    routed_equal,
                    owner,
                    component_locus_to_partition,
                );
                for (feature, share) in map {
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
        }
    }
    if compositions > 0 {
        eprintln!(
            "[junction] boundary: {} compositions ({} pooled, {} novel-only, \
             {} novel with crossing reads)",
            compositions,
            pooled_count,
            compositions - pooled_count,
            novel_spanning
        );
    }
    let mut successors = vec![Vec::new(); left_ranges.len()];
    let mut by_pair = HashMap::with_capacity(links.len());
    for (index, &(left, right)) in links.iter().enumerate() {
        let composition = link_composition[index] as usize;
        let profile = std::sync::Arc::clone(&profiles[composition]);
        successors[left].push((right, composition_scores[composition].0));
        by_pair.insert((left, right), profile);
    }
    Ok(SpineBoundary {
        successors,
        seam_by_pair: by_pair,
        compositions,
    })
}

/// The set of features any seam composition at a boundary can predict — the
/// boundary-explainable universe of the orphan-share diagnostic. Draft
/// variant: the same set before the boundary is finalized.
pub(super) fn draft_feature_set(draft: &SpineBoundaryDraft) -> HashSet<FeatureKey> {
    let mut set = HashSet::new();
    for profile in &draft.profiles {
        for key in profile.keys() {
            set.insert(key.clone());
        }
    }
    set
}

/// The set of features any seam composition at a boundary can predict — the
/// boundary-explainable universe of the orphan-share diagnostic.
pub(super) fn boundary_feature_set(
    seam_by_pair: &HashMap<(usize, usize), std::sync::Arc<Profile>>,
) -> HashSet<FeatureKey> {
    let mut set = HashSet::new();
    for profile in seam_by_pair.values() {
        for key in profile.keys() {
            set.insert(key.clone());
        }
    }
    set
}

// ---------------------------------------------------------------------------
// Viability (topological legality: forward/backward reachability through
// the legal links, plus per-allele local span feasibility).
// ---------------------------------------------------------------------------

pub(super) fn spine_viability(
    ranges: &[Vec<genome::SpanningTraversal>],
    successors: &[Vec<Vec<(usize, f64)>>],
) -> Vec<Vec<bool>> {
    let locus_count = ranges.len();
    let mut forward: Vec<Vec<bool>> = ranges.iter().map(|l| vec![true; l.len()]).collect();
    for locus in (0..locus_count - 1).rev() {
        for (allele, list) in successors[locus].iter().enumerate() {
            forward[locus][allele] = list.iter().any(|&(next, _)| forward[locus + 1][next]);
        }
    }
    let mut backward: Vec<Vec<bool>> = ranges.iter().map(|l| vec![false; l.len()]).collect();
    if locus_count > 0 {
        backward[0].iter_mut().for_each(|value| *value = true);
    }
    for locus in 1..locus_count {
        for (allele, list) in successors[locus - 1].iter().enumerate() {
            if backward[locus - 1][allele] {
                for &(next, _) in list {
                    backward[locus][next] = true;
                }
            }
        }
    }
    (0..locus_count)
        .map(|locus| {
            ranges[locus]
                .iter()
                .enumerate()
                .map(|(allele, traversal)| {
                    forward[locus][allele]
                        && backward[locus][allele]
                        && spans_feasible_local(&traversal.segments)
                })
                .collect()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Multiplicity-background census (the panplexity companion diagnostic):
// reported statistics only — the per-feature multiplicity distribution,
// the derived beta_f distribution, and a cheap Shannon-entropy linguistic-
// complexity proxy per locus window (the panel-context diversity of the
// features each territory window charges). NOTHING is filtered: the
// Poisson discount is the whole model change.
// ---------------------------------------------------------------------------

fn multiplicity_histogram(paths: u32) -> u32 {
    match paths {
        0 => 0,
        1 => 1,
        2 => 2,
        3..=4 => 3,
        5..=8 => 4,
        9..=16 => 5,
        17..=32 => 6,
        33..=64 => 7,
        65..=128 => 8,
        129..=256 => 9,
        _ => 10,
    }
}

pub(super) fn multiplicity_census(
    locus_classes: &[LocusClassing],
    backgrounds: &FeatureBackgrounds,
    model: &ScoreModel,
    universe_size: usize,
    background_seconds: f64,
) -> serde_json::Value {
    let per_count = model.histogram as f64 * model.depth / model.denominator;
    // Global multiplicity histogram + attribution-weight and background
    // distributions over every scanned feature.
    let mut histogram = [0u64; 11];
    let mut betas: Vec<f64> = Vec::with_capacity(backgrounds.len());
    let mut weights: Vec<f64> = Vec::with_capacity(backgrounds.len());
    let mut multiplicity: Vec<u32> = Vec::with_capacity(backgrounds.len());
    for (_, entry) in backgrounds.iter_entries() {
        histogram[multiplicity_histogram(entry.paths) as usize] += 1;
        betas.push(entry.beta);
        weights.push(entry.weight);
        multiplicity.push(entry.paths);
    }
    betas.sort_by(|a, b| a.total_cmp(b));
    weights.sort_by(|a, b| a.total_cmp(b));
    let quantile = |values: &[f64], fraction: f64| -> f64 {
        if values.is_empty() {
            return 0.0;
        }
        let index = ((values.len() as f64 - 1.0) * fraction).round() as usize;
        values[index.min(values.len() - 1)]
    };
    multiplicity.sort_unstable();
    // Per-locus window diversity: the Shannon entropy (nats) of the locus's
    // feature multiplicity distribution, the effective panel-context count
    // exp(H), and the locus's hardest-discounted features.
    let per_locus: Vec<serde_json::Value> = (0..locus_classes.len())
        .map(|locus| {
            let mut features: BTreeSet<&FeatureKey> = BTreeSet::new();
            for profile in &locus_classes[locus].profiles {
                for key in profile.keys() {
                    features.insert(key);
                }
            }
            let mut by_multiplicity: BTreeMap<u32, u64> = BTreeMap::new();
            let mut rows: Vec<serde_json::Value> = Vec::new();
            let mut top: Vec<(&FeatureKey, FeatureBackgroundEntry)> = features
                .iter()
                .filter_map(|&key| backgrounds.entry(key).map(|entry| (key, entry)))
                .collect();
            for (_, entry) in &top {
                *by_multiplicity.entry(entry.paths).or_default() += 1;
            }
            top.sort_by(|a, b| b.1.paths.cmp(&a.1.paths));
            for (feature, entry) in top.iter().take(8) {
                rows.push(serde_json::json!({
                    "feature": feature,
                    "paths": entry.paths,
                    "incidence": entry.incidence,
                    "attribution_weight": entry.weight,
                    "beta": entry.beta,
                    "per_copy_expected_s": per_count * entry.incidence as f64
                        / entry.paths.max(1) as f64,
                }));
            }
            let total = features.len() as f64;
            let entropy: f64 = by_multiplicity
                .values()
                .map(|&count| {
                    let p = count as f64 / total;
                    -p * p.ln()
                })
                .sum();
            let unique = by_multiplicity.get(&1).copied().unwrap_or(0);
            serde_json::json!({
                "locus": locus,
                "features": features.len(),
                "multiplicity_entropy_nats": entropy,
                "effective_panel_contexts": entropy.exp(),
                "unique_fraction": if total > 0.0 { unique as f64 / total } else { 0.0 },
                "top_attributed_features": rows,
            })
        })
        .collect();
    serde_json::json!({
        "model": "attribution: loss_f = s - (C_f/m_f)*ln(1+s/base); (b) diagnostic: beta_f = base + (m_f-1)*per_count",
        "variant": format!("{:?}", multiplicity_variant()),
        "base_background": model.background,
        "per_count": per_count,
        "universe_features": universe_size,
        "scanned_features": backgrounds.len(),
        "scan_wall_seconds": background_seconds,
        "scans": backgrounds.scans,
        "attribution_weight_quantiles": {
            "min": quantile(&weights, 0.0),
            "q25": quantile(&weights, 0.25),
            "median": quantile(&weights, 0.5),
            "q75": quantile(&weights, 0.75),
            "max": quantile(&weights, 1.0),
        },
        "beta_quantiles": {
            "min": quantile(&betas, 0.0),
            "q25": quantile(&betas, 0.25),
            "median": quantile(&betas, 0.5),
            "q75": quantile(&betas, 0.75),
            "max": quantile(&betas, 1.0),
        },
        "multiplicity_histogram": {
            "absent": histogram[0],
            "1": histogram[1],
            "2": histogram[2],
            "3-4": histogram[3],
            "5-8": histogram[4],
            "9-16": histogram[5],
            "17-32": histogram[6],
            "33-64": histogram[7],
            "65-128": histogram[8],
            "129-256": histogram[9],
            ">256": histogram[10],
        },
        "per_locus_window_diversity": per_locus,
    })
}

// ---------------------------------------------------------------------------
// The spine (stages 1-4).
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
pub(super) fn run_local_first_spine(
    panel: &SyngIndex,
    sources: &routes::Sources,
    graph: &routes::Graph,
    routes_dir: &PathBuf,
    path_of_source: &[usize],
    axis_slice: &[genome::AxisInterval],
    mut ranges: Vec<Vec<genome::SpanningTraversal>>,
    territory: &[Vec<SourceRange>],
    k: u64,
    routed_equal: &RoutedObs,
    // The record-once site-map builder (Fix 1's mixed-owner charging) and
    // the windows' full observed profiles (Fix 2's omission universes,
    // indexed by the slice-local locus).
    site: &crate::SiteObserved,
    window_obs: &[HashMap<FeatureKey, f64>],
    // The windows' observed INSTANCE structure (the instance-level
    // exoneration's observed side; `None` keeps the established
    // feature-level chain-level form bit-identically).
    window_instances: Option<&crate::InstanceStructure>,
    // The instance-probe diagnostic's feature list (empty = no probe
    // report; an operational diagnostic input, not a model constant).
    instance_probe_features: &[FeatureKey],
    component_locus_to_partition: &[u32],
    locus_offset: usize,
    target_source: usize,
    _target_length: u64,
    model: &ScoreModel,
    sample_counts: &impg::sample_mem_bwt::WeightedBwt,
    reference_specs: &[(String, PathBuf, PathBuf)],
    sweep_only: bool,
    // Correlation phasing (owner-directed exploration): replace the
    // port-word-viability chain DP with the co-occurrence + boundary-evidence
    // phasing DP (no port-word legality on transitions).
    phasing: bool,
    // Step-A measurement (owner-directed 2026-09-23): when present, the
    // per-feature haploid single-allele local table between the previously
    // selected chain's slot-0 rows (first path) and the truth reference's
    // copy-0 pieces is decomposed and written (second path). No model
    // change — a measurement over the run's own structures.
    decompose: Option<&(PathBuf, PathBuf)>,
    // The LEGACY sample-side pooled support per feature (component-scoped,
    // each slice-touching record's full multiplicity per feature) —
    // DECOMPOSE-TABLE DIAGNOSTIC COLUMNS ONLY (the old-vs-new background
    // delta distribution); no charging path reads it.
    pooled_support: Option<&HashMap<FeatureKey, f64>>,
    // The haploid track's CHARGED background (owner decision 2026-09-24,
    // the genome-wide cross-support form): per feature,
    // Sum_r m_r*(1 - 1/t_r_genome) over every multiset record realizing it,
    // t_r_genome = the record's touched count across ALL components'
    // partitions (the genome-universe placement pass). The haploid
    // single-allele charges use beta_f = base + this value; the observed
    // shares stay component-scoped.
    sample_cross_support: Option<&HashMap<FeatureKey, f64>>,
    // The within-read adjacency index (the junction-evidence restriction:
    // novel junctions are paid by crossing reads only, at every seam the
    // spine scores — Stage-1 interior junctions, Stage-2 boundary links,
    // the reference ladder, and the phasing transitions).
    span: &JunctionSpanIndex,
    rescore_directory: &std::path::Path,
    rss: &mut genome::PeriodicRssGuard,
    // Probe routes (the finalist re-rank's bit-identity anchors and the
    // surrogate-rank measurement; see PhasingShared::probe_specs).
    probe_specs: &[(String, PathBuf)],
) -> io::Result<serde_json::Value> {
    let started = Instant::now();
    let locus_count = ranges.len();
    ensure(locus_count > 0, "empty spine domain")?;
    ensure(
        axis_slice.len() == locus_count,
        "axis/spine domain cardinality mismatch",
    )?;
    let flank_memo: FlankMemo = std::sync::Mutex::new(HashMap::new());
    let mut ports = routes::Ports::open_without_global_verification(routes_dir, graph)?;
    let truth_spec = reference_specs
        .iter()
        .find(|(label, _, _)| label == "truth")
        .cloned()
        .ok_or_else(|| invalid("--spine needs a --reference-pair labeled truth"))?;

    let partition_obs: Vec<HashMap<FeatureKey, f64>> = (0..locus_count)
        .map(|locus| {
            let partition = component_locus_to_partition[locus_offset + locus];
            ensure(partition != u32::MAX, "spine locus outside the routing universe")?;
            Ok(routed_equal.get(&partition).cloned().unwrap_or_default())
        })
        .collect::<io::Result<_>>()?;

    // ------------------------------------------------ Stage 1a: classing.
    let classing_started = Instant::now();
    let mut locus_classes: Vec<LocusClassing> = Vec::with_capacity(locus_count);
    let mut classing_substage = [0.0f64; 4];
    let mut split_seam_queries = 0u64;
    let mut ir_extractions = 0u64;
    for locus in 0..locus_count {
        let classing = class_locus_alleles(
            panel,
            sources,
            &flank_memo,
            path_of_source,
            &ranges[locus],
            locus,
            k,
            span,
            component_locus_to_partition[locus_offset + locus],
            routed_equal,
            component_locus_to_partition,
            model,
        )?;
        for (slot, value) in classing.substage.iter().enumerate() {
            classing_substage[slot] += value;
        }
        split_seam_queries += classing.split_seam_queries;
        ir_extractions += classing.ir_extractions;
        locus_classes.push(classing);
        rss_probe(rss, &format!("spine_classing_locus_{locus}"))?;
    }
    let mut class_counts: Vec<usize> = locus_classes.iter().map(|l| l.profiles.len()).collect();
    let classing_seconds = classing_started.elapsed().as_secs_f64();
    eprintln!(
        "[spine] classing done: {:?} classes (subtotals parents {:.2}s sweeps {:.2}s \
         partials+seams {:.2}s classing {:.2}s)",
        class_counts,
        classing_substage[0],
        classing_substage[1],
        classing_substage[2],
        classing_substage[3]
    );

    // ------------------------------------------------ boundaries (drafts).
    let boundary_started = Instant::now();
    let mut boundary_drafts: Vec<SpineBoundaryDraft> = Vec::with_capacity(locus_count - 1);
    for boundary in 0..locus_count - 1 {
        let axis_window = &axis_slice[boundary..boundary + 2];
        let range_window = &ranges[boundary..boundary + 2];
        let links = genome::port_word_seams(axis_window, range_window, graph, &mut ports)?;
        let links = links
            .into_iter()
            .next()
            .ok_or_else(|| invalid("spine boundary link cardinality mismatch"))?;
        ensure(!links.is_empty(), "spine boundary {boundary} has no legal links")?;
        let left_owners: Vec<u32> = ranges[boundary]
            .iter()
            .map(crate::traversal_owner)
            .collect();
        let right_owners: Vec<u32> = ranges[boundary + 1]
            .iter()
            .map(crate::traversal_owner)
            .collect();
        let draft = build_spine_boundary_draft(
            panel,
            sources,
            &flank_memo,
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
        )?;
        boundary_drafts.push(draft);
        rss_probe(rss, &format!("spine_boundary_{boundary}"))?;
    }
    let draft_seconds = boundary_started.elapsed().as_secs_f64();

    // ------------------------------------------------ multiplicity backgrounds.
    // The interior charging universe: every classed locus profile feature
    // plus every boundary seam composition feature, scanned once against
    // the panel's paths (distinct paths containing the feature and its
    // per-path L150-window incidence) to derive the per-feature Poisson
    // backgrounds beta_f = base + (m_f - 1) * per_copy_expected_s_f.
    let background_started = Instant::now();
    let mut backgrounds = FeatureBackgrounds::new(model.background);
    let background_universe_size: usize;
    {
        let mut universe: BTreeSet<FeatureKey> = BTreeSet::new();
        for locus in 0..locus_count {
            for profile in &locus_classes[locus].profiles {
                for key in profile.keys() {
                    universe.insert(key.clone());
                }
            }
        }
        for draft in &boundary_drafts {
            for key in draft_feature_set(draft) {
                universe.insert(key);
            }
        }
        background_universe_size = universe.len();
        let features: Vec<FeatureKey> = universe.into_iter().collect();
        backgrounds.scan_extend(panel, features, k, model)?;
        eprintln!(
            "[backgrounds] multiplicity scan: {} universe features, {:.2}s",
            background_universe_size,
            background_started.elapsed().as_secs_f64()
        );
    }
    let background_seconds = background_started.elapsed().as_secs_f64();
    rss_probe(rss, "spine_background_scan")?;

    // ------------------------------------------------ boundaries (finalize).
    let finalize_started = Instant::now();
    let mut boundaries: Vec<SpineBoundary> = Vec::with_capacity(locus_count - 1);
    for (boundary, draft) in boundary_drafts.into_iter().enumerate() {
        let left_owners: Vec<u32> =
            ranges[boundary].iter().map(crate::traversal_owner).collect();
        let right_owners: Vec<u32> =
            ranges[boundary + 1].iter().map(crate::traversal_owner).collect();
        let boundary_data = finalize_spine_boundary(
            draft,
            &ranges[boundary],
            &ranges[boundary + 1],
            &left_owners,
            &right_owners,
            routed_equal,
            component_locus_to_partition,
            model,
            &backgrounds,
        )?;
        boundaries.push(boundary_data);
    }
    let boundary_seconds = draft_seconds + finalize_started.elapsed().as_secs_f64();
    eprintln!(
        "[spine] boundaries done: compositions {:?}",
        boundaries.iter().map(|b| b.compositions).collect::<Vec<_>>()
    );

    // ------------------------------------------------ viability.
    let mut successors_ref: Vec<Vec<Vec<(usize, f64)>>> = boundaries
        .iter()
        .map(|b| b.successors.clone())
        .collect();
    let mut viable = spine_viability(&ranges, &successors_ref);

    // ------------------------------------------------ folded tables + Stage 1b sweep.
    let sweep_started = Instant::now();
    let mut folded: Vec<LocusFolded> = (0..locus_count)
        .into_par_iter()
        .map(|locus| {
            LocusFolded::build_backgrounds(
                &locus_classes[locus].profiles,
                &locus_classes[locus].class_owners,
                routed_equal,
                component_locus_to_partition,
                &backgrounds,
            )
        })
        .collect::<io::Result<_>>()?;
    let mut loss_tables: Vec<Vec<Vec<f64>>> = (0..locus_count)
        .into_par_iter()
        .map(|locus| folded[locus].loss_tables(model).map(|(tables, _)| tables))
        .collect::<io::Result<_>>()?;

    // ------------------------------------------------ native backbone + incumbent.
    let backbone = native_backbone_chain(
        &ranges,
        &successors_ref,
        &locus_classes,
        routed_equal,
        component_locus_to_partition,
        target_source,
        model,
        &backgrounds,
    )?;
    let (backbone_chain, backbone_cumulative) = backbone
        .ok_or_else(|| invalid("spine native backbone chain is broken"))?;
    let incumbent = *backbone_cumulative.last().expect("backbone total");
    ensure(incumbent.is_finite(), "nonfinite spine incumbent");
    eprintln!("[spine] native backbone incumbent {incumbent:.2}");

    // Exhaustive per-locus sweeps (sequential over loci; STEP-3: the
    // admissible-bound branch-and-bound with margin-sentinel pruning inside).
    let seam_swings = seam_swings_of(locus_count, &successors_ref, &viable);
    let mut sweeps: Vec<LocusSweep> = Vec::with_capacity(locus_count);
    for locus in 0..locus_count {
        sweeps.push(exhaustive_local_sweep(
            &folded[locus],
            &loss_tables[locus],
            &locus_classes[locus].membership,
            &viable[locus],
            Some(backbone_chain[locus]),
            model,
            seam_swings[locus],
            &locus_classes[locus].class_charges,
            &scorable_classes(&locus_classes[locus], &ranges[locus]),
        )?);
        rss_probe(rss, &format!("spine_sweep_locus_{locus}"))?;
    }
    let sweep_seconds = sweep_started.elapsed().as_secs_f64();
    let total_class_pairs: u64 = sweeps.iter().map(|s| s.class_pairs).sum();
    let total_allele_pairs: u64 = sweeps.iter().map(|s| s.allele_pairs).sum();
    eprintln!(
        "[spine] exhaustive sweep done: {total_class_pairs} class pairs, \
         {total_allele_pairs} allele pairs, {:.2}s",
        sweep_seconds
    );

    // ------------------------------------------------ truth decomposition (assessment side).
    let truth_started = Instant::now();
    let truth_route_a: routes::Route = read_json(&truth_spec.1)?;
    let truth_route_b: routes::Route = read_json(&truth_spec.2)?;
    let truth_pieces = reference_local_piece_lists(
        territory,
        [&truth_route_a, &truth_route_b],
        locus_count,
        locus_offset,
    );
    let reference_partitions: Vec<u32> = (0..locus_count)
        .map(|locus| component_locus_to_partition[locus_offset + locus])
        .collect();
    let (truth_losses, truth_junction_census) = spine_reference_local_losses(
        panel,
        sources,
        path_of_source,
        &truth_pieces,
        routed_equal,
        component_locus_to_partition,
        k,
        model,
        span,
        Some(&mut backgrounds),
    )?;
    let truth_seconds = truth_started.elapsed().as_secs_f64();
    // The truth pair's donor sources (non-target sources — SK1 here).
    let donor_sources: BTreeSet<usize> = truth_route_a
        .segments
        .iter()
        .chain(&truth_route_b.segments)
        .map(|segment| segment.source)
        .filter(|&source| source != target_source)
        .collect();

    // ------------------------------------------------ Step-A measurement (owner-directed
    // 2026-09-23): the per-feature haploid single-allele local table,
    // a previously selected chain's slot-0 rows vs the truth mosaic's
    // copy-0 pieces (no model change; full table to the requested file,
    // reconciliation summary into the run log).
    let decompose_summary = match decompose {
        Some((selected_path, out_path)) => Some(decompose_haploid_local_table(
            panel,
            sources,
            path_of_source,
            axis_slice,
            &territory[locus_offset..locus_offset + locus_count],
            selected_path,
            &truth_pieces,
            &partition_obs,
            &reference_partitions,
            routed_equal,
            site,
            window_obs,
            component_locus_to_partition,
            &ranges,
            &locus_classes,
            &donor_sources,
            pooled_support,
            sample_cross_support,
            locus_offset,
            k,
            model,
            span,
            &mut backgrounds,
            out_path,
        )?),
        None => None,
    };
    let decompose_seconds = match &decompose_summary {
        Some(summary) => summary["wall_seconds"].as_f64().unwrap_or(0.0),
        None => 0.0,
    };
    // The sample-side cross-support backgrounds for the haploid track
    // (the genome-wide form): built from the per-feature genome-wide
    // cross-support map whenever the routing pass produced it (phasing runs
    // and decomposition runs).
    let sample_backgrounds = sample_cross_support
        .map(|cross| SampleSideBackgrounds::new(cross, model.background));

    // ------------------------------------------------ margins (admissible, data-derived).
    // Seam-score extremes per boundary over links whose BOTH sides are
    // viable: the boundary-evidence variation a chain can actually trade
    // local loss against (non-viable alleles never enter a chain, so their
    // link scores are not compensatory evidence).
    let seam_extremes = |successors: &[Vec<Vec<(usize, f64)>>],
                         viable: &[Vec<bool>]|
     -> Vec<(f64, f64)> {
        successors
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
            .collect()
    };
    // LOCAL margin (the retention rule the spine is specified with): a pair
    // whose local loss exceeds the locus's viable best by more than the
    // total seam-evidence swing at its two adjacent boundaries cannot be
    // compensated by any boundary evidence — retained pairs are those
    // within that swing of the local best. The margin is computed from the
    // run's own seam-score distribution and reported per locus; the
    // all-native backbone pair is structurally retained (the incumbent
    // chain that anchors admissibility and feasibility).
    let compute_local_margins = |min_pair: &[f64],
                                 successors: &[Vec<Vec<(usize, f64)>>],
                                 viable: &[Vec<bool>]|
     -> Vec<f64> {
        let extremes = seam_extremes(successors, viable);
        (0..locus_count)
            .map(|locus| {
                let mut swing = 0.0f64;
                if locus > 0 {
                    let (min, max) = extremes[locus - 1];
                    swing += (max - min).max(0.0);
                }
                if locus + 1 < locus_count {
                    let (min, max) = extremes[locus];
                    swing += (max - min).max(0.0);
                }
                min_pair[locus] + swing
            })
            .collect()
    };
    // Admissible GLOBAL margin (reported diagnostic + exactness reference):
    // any pair outside it provably lies in no chain at or below the native
    // backbone incumbent.
    let compute_margins = |min_pair: &[f64], min_seam: &[f64]| -> io::Result<(Vec<f64>, Vec<f64>)> {
        ensure(
            min_seam.iter().all(|&value| value.is_finite()),
            "a spine boundary has no legal link scores",
        )?;
        let mut suffix_from_next = vec![0.0f64; locus_count];
        for locus in (0..locus_count.saturating_sub(1)).rev() {
            suffix_from_next[locus] =
                min_pair[locus + 1] + 2.0 * min_seam[locus] + suffix_from_next[locus + 1];
        }
        let min_pair_total: f64 = min_pair.iter().sum();
        let min_seam_total: f64 = min_seam.iter().sum();
        let margins: Vec<f64> = (0..locus_count)
            .map(|locus| {
                incumbent - (min_pair_total - min_pair[locus]) - 2.0 * min_seam_total
            })
            .collect();
        Ok((margins, suffix_from_next))
    };
    let min_pair: Vec<f64> = sweeps.iter().map(|s| s.min_viable_loss).collect();
    let min_seam: Vec<f64> = successors_ref
        .iter()
        .map(|boundary| {
            boundary
                .iter()
                .flat_map(|list| list.iter().map(|&(_, score)| score))
                .reduce(f64::min)
                .unwrap_or(f64::INFINITY)
        })
        .collect();
    let (mut margins_admissible, mut suffix_from_next) = compute_margins(&min_pair, &min_seam)?;
    for locus in 0..locus_count {
        ensure(
            margins_admissible[locus] >= sweeps[locus].native_pair_loss.unwrap_or(f64::INFINITY) - 1e-9,
            "admissible margin dropped the native backbone pair (internal error)"
        )?;
    }
    // The retention margin is the LOCAL seam-swing margin; the backbone
    // pair is structurally retained per-pair inside the DP (not by raising
    // the margin floor — a floor at the native pair's loss would retain
    // every pair better than native, which is nearly everything at loci
    // where native fits poorly).
    let mut margins = compute_local_margins(&min_pair, &successors_ref, &viable);

    // ------------------------------------------------ Stage 1 winner table.
    let member_of_class = |locus: usize, class: usize| -> Option<usize> {
        locus_classes[locus]
            .membership
            .iter()
            .position(|&member| member == class)
    };
    let sources_of_pieces = |pieces: &[(usize, u64, u64, bool, u32)]| -> BTreeSet<usize> {
        pieces.iter().map(|&(source, _, _, _, _)| source).collect()
    };
    let allele_segments_json = |traversal: &genome::SpanningTraversal| -> Vec<serde_json::Value> {
        traversal
            .segments
            .iter()
            .map(|segment| {
                serde_json::json!({
                    "source": segment.source,
                    "start": segment.start,
                    "end": segment.end,
                    "reverse": segment.reverse,
                })
            })
            .collect()
    };
    // Per copy: the domain allele index whose segments exactly equal the
    // truth's merged local pieces (None when the truth piece is not a
    // parent-domain allele — junction partials, merged multi-row spans).
    let truth_pair_alleles = |locus: usize| -> [Option<usize>; 2] {
        let mut found: [Option<usize>; 2] = [None, None];
        for copy in 0..2 {
            let pieces: Vec<(usize, u64, u64, bool)> = truth_pieces[locus][copy]
                .iter()
                .map(|&(source, start, end, reverse, _owner)| (source, start, end, reverse))
                .collect();
            if pieces.is_empty() {
                continue;
            }
            found[copy] = ranges[locus].iter().position(|traversal| {
                traversal.segments.len() == pieces.len()
                    && traversal
                        .segments
                        .iter()
                        .zip(pieces.iter())
                        .all(|(segment, piece)| {
                            segment.source == piece.0
                                && segment.start == piece.1
                                && segment.end == piece.2
                                && segment.reverse == piece.3
                        })
            });
        }
        found
    };
    let mut stage1_loci: Vec<serde_json::Value> = Vec::with_capacity(locus_count);
    let mut tract_loci: Vec<usize> = Vec::new();
    for locus in 0..locus_count {
        let best_pair = sweeps[locus].best_class_pairs[0];
        let winner_alleles: Vec<Option<usize>> = best_pair
            .iter()
            .map(|&class| member_of_class(locus, class))
            .collect();
        let winner_sources: BTreeSet<usize> = winner_alleles
            .iter()
            .flatten()
            .flat_map(|&allele| {
                ranges[locus][allele]
                    .segments
                    .iter()
                    .map(|segment| segment.source)
            })
            .collect();
        let truth_sources: BTreeSet<usize> = truth_pieces[locus]
            .iter()
            .flat_map(|pieces| sources_of_pieces(pieces))
            .collect();
        let allele_matches_pieces =
            |allele: usize, pieces: &[(usize, u64, u64, bool)]| -> bool {
                ranges[locus][allele].segments.len() == pieces.len()
                    && ranges[locus][allele]
                        .segments
                        .iter()
                        .zip(pieces.iter())
                        .all(|(segment, piece)| {
                            segment.source == piece.0
                                && segment.start == piece.1
                                && segment.end == piece.2
                                && segment.reverse == piece.3
                        })
            };
        let truth_pieces_a: Vec<(usize, u64, u64, bool)> = truth_pieces[locus][0]
            .iter()
            .map(|&(source, start, end, reverse, _owner)| (source, start, end, reverse))
            .collect();
        let truth_pieces_b: Vec<(usize, u64, u64, bool)> = truth_pieces[locus][1]
            .iter()
            .map(|&(source, start, end, reverse, _owner)| (source, start, end, reverse))
            .collect();
        // The winner pair (unordered) is the truth pair iff it matches the
        // truth's two piece lists in some assignment.
        let winner_is_truth = match (winner_alleles[0], winner_alleles[1]) {
            (Some(a), Some(b)) => {
                (allele_matches_pieces(a, &truth_pieces_a)
                    && allele_matches_pieces(b, &truth_pieces_b))
                    || (allele_matches_pieces(a, &truth_pieces_b)
                        && allele_matches_pieces(b, &truth_pieces_a))
            }
            _ => false,
        };
        let truth_alleles = truth_pair_alleles(locus);
        let truth_pair_viable = truth_alleles
            .iter()
            .all(|allele| allele.is_some_and(|allele| viable[locus][allele]));
        let truth_link_status: Vec<serde_json::Value> = (0..2)
            .map(|copy| match truth_alleles[copy] {
                None => serde_json::Value::Null,
                Some(allele) => {
                    let has_in = locus == 0
                        || successors_ref[locus - 1]
                            .iter()
                            .any(|list| list.iter().any(|&(next, _)| next == allele));
                    let has_out = locus + 1 == locus_count
                        || !successors_ref[locus][allele].is_empty();
                    serde_json::json!({
                        "allele": allele,
                        "viable": viable[locus][allele],
                        "has_in_links": has_in,
                        "has_out_links": has_out,
                    })
                }
            })
            .collect();
        let native_pair_loss = sweeps[locus].native_pair_loss;
        let truth_loss = truth_losses[locus];
        let truth_rank = sweeps[locus]
            .table
            .iter()
            .filter(|&&value| (value.total_cmp(&truth_loss)).is_lt())
            .count() as u64;
        let is_tract_locus = truth_sources
            .iter()
            .any(|source| donor_sources.contains(source));
        if is_tract_locus {
            tract_loci.push(locus);
        }
        let donor_in_winner = winner_sources
            .iter()
            .any(|&source| source != target_source);
        let truth_donor_in_winner = winner_sources
            .iter()
            .any(|source| donor_sources.contains(source));
        let best_viable_pair = sweeps[locus].best_viable_class_pairs[0];
        let viable_winner_alleles: Vec<Option<usize>> = best_viable_pair
            .iter()
            .map(|&class| member_of_class(locus, class))
            .collect();
        let viable_winner_sources: BTreeSet<usize> = viable_winner_alleles
            .iter()
            .flatten()
            .flat_map(|&allele| {
                ranges[locus][allele]
                    .segments
                    .iter()
                    .map(|segment| segment.source)
            })
            .collect();
        let truth_donor_in_viable_winner = viable_winner_sources
            .iter()
            .any(|source| donor_sources.contains(source));
        let native_in_viable_winner = viable_winner_sources
            .iter()
            .any(|&source| source == target_source);
        stage1_loci.push(serde_json::json!({
            "locus": locus,
            "full_locus": locus_offset + locus,
            "axis_interval": [axis_slice[locus].start, axis_slice[locus].end],
            "alleles": ranges[locus].len(),
            "classes": class_counts[locus],
            "viable_alleles": viable[locus].iter().filter(|&&v| v).count(),
            "class_pairs": sweeps[locus].class_pairs,
            "allele_pairs": sweeps[locus].allele_pairs,
            "viable_allele_pairs": sweeps[locus].viable_allele_pairs,
            "best_loss": sweeps[locus].best_loss,
            "best_class_pair_ties": sweeps[locus].best_class_pairs,
            "best_tie_count": sweeps[locus].best_class_pairs.len(),
            "winner_pair_alleles": winner_alleles
                .iter()
                .map(|allele| allele.map(|a| allele_segments_json(&ranges[locus][a])))
                .collect::<Vec<_>>(),
            "second_best_loss": sweeps[locus].second_loss,
            "margin_to_second": sweeps[locus].second_loss.map(|s| s - sweeps[locus].best_loss),
            "native_pair_loss": native_pair_loss,
            "native_pair_rank": sweeps[locus].native_rank,
            "truth_pair_loss": truth_loss,
            "truth_pair_rank": truth_rank,
            "truth_pair_in_domain":
                [truth_alleles[0].is_some(), truth_alleles[1].is_some()],
            "truth_pair_viable": truth_pair_viable,
            "truth_allele_link_status": truth_link_status,
            "winner_sources": winner_sources.into_iter().collect::<Vec<_>>(),
            "truth_sources": truth_sources.into_iter().collect::<Vec<_>>(),
            "winner_is_truth_pair": winner_is_truth,
            "winner_contains_donor": donor_in_winner,
            "winner_contains_truth_donor": truth_donor_in_winner,
            "best_viable_loss": sweeps[locus].best_viable_loss,
            "best_viable_class_pair_ties": sweeps[locus].best_viable_class_pairs,
            "viable_winner_sources": viable_winner_sources.into_iter().collect::<Vec<_>>(),
            "viable_winner_contains_truth_donor": truth_donor_in_viable_winner,
            "viable_winner_contains_native": native_in_viable_winner,
            "truth_vs_native_margin": native_pair_loss.map(|native| native - truth_loss),
            "admissible_floor": sweeps[locus].floor,
            "loss_residual": sweeps[locus].best_loss - sweeps[locus].floor,
            "loss_quantiles_class_pairs": sweeps[locus].quantiles,
            "retained_margin_local": margins[locus],
            "retained_margin_admissible": margins_admissible[locus],
            "min_viable_pair_loss": sweeps[locus].min_viable_loss,
            "is_truth_tract_locus": is_tract_locus,
        }));
    }

    // Orphan-share diagnostic (observed support no parent interior and no
    // adjacent boundary seam can predict) — per locus.
    let orphan_started = Instant::now();
    let compute_orphan_shares = || -> io::Result<Vec<(usize, f64)>> {
        let mut shares = Vec::with_capacity(locus_count);
        for locus in 0..locus_count {
            let mut class_set: HashSet<FeatureKey> = HashSet::new();
            for profile in &locus_classes[locus].profiles {
                for key in profile.keys() {
                    class_set.insert(key.clone());
                }
            }
            let mut explainable = class_set;
            if locus > 0 {
                for key in boundary_feature_set(&boundaries[locus - 1].seam_by_pair) {
                    explainable.insert(key);
                }
            }
            if locus + 1 < locus_count {
                for key in boundary_feature_set(&boundaries[locus].seam_by_pair) {
                    explainable.insert(key);
                }
            }
            let mut total = 0.0f64;
            let mut orphan = 0.0f64;
            for (key, &value) in &partition_obs[locus] {
                total += value;
                if !explainable.contains(key) {
                    orphan += value;
                }
            }
            let share = if total > 0.0 { orphan / total } else { 0.0 };
            shares.push((locus, share));
        }
        Ok(shares)
    };
    let orphan_shares = compute_orphan_shares()?;
    let orphan_seconds = orphan_started.elapsed().as_secs_f64();
    for (index, row) in stage1_loci.iter_mut().enumerate() {
        if let Some(object) = row.as_object_mut() {
            object.insert(
                "orphan_observed_share".into(),
                serde_json::json!(orphan_shares[index].1),
            );
        }
    }

    // ---- the tract verdict (model-vs-search question).
    let tract_rows: Vec<serde_json::Value> = tract_loci
        .iter()
        .map(|&locus| {
            let row = &stage1_loci[locus];
            serde_json::json!({
                "locus": locus,
                "axis_interval": row["axis_interval"],
                "winner_is_truth": row["winner_is_truth_pair"],
                "winner_contains_donor": row["winner_contains_donor"],
                "truth_pair_loss": row["truth_pair_loss"],
                "native_pair_loss": row["native_pair_loss"],
                "best_loss": row["best_loss"],
                "truth_vs_native_margin": row["truth_vs_native_margin"],
                "truth_pair_viable": row["truth_pair_viable"],
                "best_viable_loss": row["best_viable_loss"],
                "viable_winner_contains_truth_donor":
                    row["viable_winner_contains_truth_donor"],
                "viable_winner_contains_native": row["viable_winner_contains_native"],
                "in_c19_beam_native_region":
                    axis_slice[locus].start < 194421 && axis_slice[locus].end > 113881,
            })
        })
        .collect();
    let tract_winner_truth = tract_loci
        .iter()
        .filter(|&&locus| stage1_loci[locus]["winner_is_truth_pair"].as_bool().unwrap_or(false))
        .count();
    let tract_winner_donor = tract_loci
        .iter()
        .filter(|&&locus| stage1_loci[locus]["winner_contains_donor"].as_bool().unwrap_or(false))
        .count();
    let tract_winner_truth_donor = tract_loci
        .iter()
        .filter(|&&locus| {
            stage1_loci[locus]["winner_contains_truth_donor"]
                .as_bool()
                .unwrap_or(false)
        })
        .count();
    let tract_viable_winner_truth_donor = tract_loci
        .iter()
        .filter(|&&locus| {
            stage1_loci[locus]["viable_winner_contains_truth_donor"]
                .as_bool()
                .unwrap_or(false)
        })
        .count();
    let tract_truth_beats_native = tract_loci
        .iter()
        .filter(|&&locus| {
            stage1_loci[locus]["truth_vs_native_margin"]
                .as_f64()
                .is_some_and(|margin| margin > 0.0)
        })
        .count();
    let tract_truth_pair_in_domain = tract_loci
        .iter()
        .filter(|&&locus| {
            stage1_loci[locus]["truth_pair_in_domain"].as_array().is_some_and(|flags| {
                flags.len() == 2
                    && flags[0].as_bool().unwrap_or(false)
                    && flags[1].as_bool().unwrap_or(false)
            })
        })
        .count();
    let tract_truth_pair_viable = tract_loci
        .iter()
        .filter(|&&locus| {
            stage1_loci[locus]["truth_pair_viable"]
                .as_bool()
                .unwrap_or(false)
        })
        .count();
    let left_region: Vec<usize> = tract_loci
        .iter()
        .copied()
        .filter(|&locus| axis_slice[locus].start < 194421 && axis_slice[locus].end > 113881)
        .collect();
    let left_region_truth = left_region
        .iter()
        .filter(|&&locus| stage1_loci[locus]["winner_is_truth_pair"].as_bool().unwrap_or(false))
        .count();
    let left_region_donor = left_region
        .iter()
        .filter(|&&locus| stage1_loci[locus]["winner_contains_truth_donor"].as_bool().unwrap_or(false))
        .count();
    let left_region_viable_donor = left_region
        .iter()
        .filter(|&&locus| {
            stage1_loci[locus]["viable_winner_contains_truth_donor"]
                .as_bool()
                .unwrap_or(false)
        })
        .count();
    let left_region_truth_beats_native = left_region
        .iter()
        .filter(|&&locus| {
            stage1_loci[locus]["truth_vs_native_margin"]
                .as_f64()
                .is_some_and(|margin| margin > 0.0)
        })
        .count();
    let tract_verdict = serde_json::json!({
        "tract_loci": tract_loci,
        "tract_loci_count": tract_loci.len(),
        "per_locus": tract_rows,
        "winner_is_truth_count": tract_winner_truth,
        "winner_contains_donor_count": tract_winner_donor,
        "winner_contains_truth_donor_count": tract_winner_truth_donor,
        "viable_winner_contains_truth_donor_count": tract_viable_winner_truth_donor,
        "truth_beats_native_loci_count": tract_truth_beats_native,
        "truth_pair_in_domain_count": tract_truth_pair_in_domain,
        "truth_pair_viable_count": tract_truth_pair_viable,
        "c19_beam_native_region_loci": left_region,
        "c19_region_winner_is_truth_count": left_region_truth,
        "c19_region_winner_contains_truth_donor_count": left_region_donor,
        "c19_region_viable_winner_contains_truth_donor_count": left_region_viable_donor,
        "c19_region_truth_beats_native_count": left_region_truth_beats_native,
    });

    let mut stage1_summary = serde_json::json!({
        "loci": stage1_loci,
        "totals": {
            "class_pairs": total_class_pairs,
            "allele_pairs": total_allele_pairs,
            "tract_loci": tract_loci.len(),
            "tract_winner_is_truth": tract_winner_truth,
            "tract_winner_contains_donor": tract_winner_donor,
            "tract_winner_contains_truth_donor": tract_winner_truth_donor,
            "tract_viable_winner_contains_truth_donor": tract_viable_winner_truth_donor,
            "tract_truth_beats_native": tract_truth_beats_native,
            "tract_truth_pair_in_domain": tract_truth_pair_in_domain,
            "tract_truth_pair_viable": tract_truth_pair_viable,
        },
        "tract_verdict": tract_verdict,
        "multiplicity_background": multiplicity_census(
            &locus_classes,
            &backgrounds,
            model,
            background_universe_size,
            background_seconds,
        ),
        "wall_seconds": {
            "classing": classing_seconds,
            "boundaries": boundary_seconds,
            "background_scan": background_seconds,
            "sweep": sweep_seconds,
            "truth_decomposition": truth_seconds,
            "orphan_diagnostic": orphan_seconds,
            "decompose_local_table": decompose_seconds,
            "total": classing_seconds + boundary_seconds + sweep_seconds + truth_seconds,
        },
        "classing_substage_seconds": {
            "parents": classing_substage[0],
            "sweeps": classing_substage[1],
            "partials_seams": classing_substage[2],
            "classing": classing_substage[3],
        },
        "ir_extractions": ir_extractions,
        "split_seam_queries": split_seam_queries,
        "incumbent_native_backbone": incumbent,
        "truth_interior_junction_census": truth_junction_census,
    });
    if let (Some(summary), Some(object)) = (
        &decompose_summary,
        stage1_summary.as_object_mut(),
    ) {
        object.insert("haploid_local_decomposition".into(), summary.clone());
    }

    if sweep_only {
        return Ok(serde_json::json!({
            "stage1_sweep": stage1_summary,
            "rss_peak_bytes": rss.peak_bytes(),
            "total_seconds": started.elapsed().as_secs_f64(),
        }));
    }

    // ---------------------------------------- correlation phasing (the new tail).
    if phasing {
        let summary = run_correlation_phasing(&mut PhasingShared {
            panel,
            sources,
            graph,
            ports: &mut ports,
            flank_memo: &flank_memo,
            path_of_source,
            axis_slice,
            territory,
            partition_obs: &partition_obs,
            routed_equal,
            component_locus_to_partition,
            locus_partitions: reference_partitions.clone(),
            k,
            locus_offset,
            target_source,
            model,
            sample_counts,
            reference_specs: &reference_specs,
            truth_route_a: &truth_route_a,
            truth_route_b: &truth_route_b,
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
            incumbent_spine: incumbent,
            backbone_chain,
            truth_pieces,
            truth_losses,
            donor_sources,
            tract_loci,
            stage1_summary,
            backgrounds: &mut backgrounds,
            sample_backgrounds: sample_backgrounds.as_ref(),
            site,
            window_obs,
            window_instances,
            instance_probe_features: instance_probe_features.to_vec(),
            probe_specs,
        })?;
        return Ok(summary);
    }

    // ------------------------------------------------ Stage 3: targeted splits.
    // The residual detector (best local pair versus the admissible
    // observations-can-explain floor) selects loci through the natural gap
    // of its own distribution; the split machinery runs AT THOSE LOCI ONLY,
    // with full cut enumeration there (no rank checkpoints, no top-N). A
    // locus is split at most once; rounds repeat while the detector demands
    // fresh loci.
    let splits_started = Instant::now();
    let conflict_window = conflict_window_of(&ranges).max(1);
    let mut split_done = vec![false; locus_count];
    let mut split_candidate_counts = vec![0usize; locus_count];
    let mut split_admitted_counts = vec![0usize; locus_count];
    let mut split_rounds: Vec<serde_json::Value> = Vec::new();
    let mut min_pair = sweeps.iter().map(|s| s.min_viable_loss).collect::<Vec<f64>>();
    let mut min_seam = successors_ref
        .iter()
        .map(|boundary| {
            boundary
                .iter()
                .flat_map(|list| list.iter().map(|&(_, score)| score))
                .reduce(f64::min)
                .unwrap_or(f64::INFINITY)
        })
        .collect::<Vec<f64>>();
    let mut margins = margins;
    let mut suffix_from_next = suffix_from_next;
    let mut orphan_shares = orphan_shares;
    let mut seam_feature_sets: Vec<HashSet<FeatureKey>> = boundaries
        .iter()
        .map(|boundary| boundary_feature_set(&boundary.seam_by_pair))
        .collect();
    loop {
        let residuals: Vec<(usize, f64)> = (0..locus_count)
            .map(|locus| {
                (
                    locus,
                    (sweeps[locus].best_loss - sweeps[locus].floor).max(0.0),
                )
            })
            .collect();
        let (selected, residual_diag) = natural_gap_selection(&residuals);
        let (orphan_selected, orphan_diag) = natural_gap_selection(&orphan_shares);
        let fresh: Vec<usize> = selected
            .iter()
            .copied()
            .filter(|&locus| !split_done[locus])
            .collect();
        split_rounds.push(serde_json::json!({
            "residual_rule": {"selected": selected, "diagnostics": residual_diag},
            "orphan_share_rule": {"selected": orphan_selected, "diagnostics": orphan_diag},
            "fresh_split_loci": fresh,
        }));
        eprintln!(
            "[spine] stage3 detector: residual selected {:?}, orphan selected {:?}, fresh {:?}",
            selected, orphan_selected, fresh
        );
        if fresh.is_empty() {
            break;
        }
        for &locus in &fresh {
            split_done[locus] = true;
            let full_locus = locus_offset + locus;
            // Admission: the forward single-segment alleles of the
            // margin-retained class pairs at this locus (the locally
            // plausible set — evidence-derived, no top-N).
            let classes = locus_classes[locus].profiles.len();
            let mut members_by_class: Vec<Vec<usize>> = vec![Vec::new(); classes];
            for (allele, &class) in locus_classes[locus].membership.iter().enumerate() {
                if viable[locus][allele] {
                    members_by_class[class].push(allele);
                }
            }
            let cb = locus_classes[locus].membership[backbone_chain[locus]];
            let mut admitted: Vec<usize> = Vec::new();
            for c1 in 0..classes {
                if members_by_class[c1].is_empty() {
                    continue;
                }
                for c2 in c1..classes {
                    if members_by_class[c2].is_empty() {
                        continue;
                    }
                    let retained = sweeps[locus].table[class_pair_index(c1, c2)]
                        <= margins[locus]
                        || (c1 == cb && c2 == cb);
                    if !retained {
                        continue;
                    }
                    for &class in &[c1, c2] {
                        for &allele in &members_by_class[class] {
                            let candidate = &ranges[locus][allele];
                            if candidate.segments.len() == 1
                                && !candidate.segments[0].reverse
                                && candidate.segments[0].start < candidate.segments[0].end
                                && !admitted.contains(&allele)
                            {
                                admitted.push(allele);
                            }
                        }
                    }
                }
            }
            admitted.sort_unstable();
            admitted.dedup();
            split_admitted_counts[locus] = admitted.len();
            let round_started = Instant::now();
            let split = genome::split_candidates(
                full_locus,
                &ranges[locus],
                &admitted,
                graph,
                &mut ports,
            )?;
            split_candidate_counts[locus] = split.len();
            eprintln!(
                "[spine] stage3 locus {locus}: admitted {} alleles -> {} split candidates ({:.2}s)",
                admitted.len(),
                split.len(),
                round_started.elapsed().as_secs_f64()
            );
            ranges[locus].extend(split);
        }
        // Re-derive at the fresh loci: classing, boundary drafts, background
        // extension for the new features, folded tables, boundary finalize.
        for &locus in &fresh {
            locus_classes[locus] = class_locus_alleles(
                panel,
                sources,
                &flank_memo,
                path_of_source,
                &ranges[locus],
                locus,
                k,
                span,
                component_locus_to_partition[locus_offset + locus],
                routed_equal,
                component_locus_to_partition,
                model,
            )?;
            class_counts[locus] = locus_classes[locus].profiles.len();
        }
        let mut affected: BTreeSet<usize> = BTreeSet::new();
        for &locus in &fresh {
            if locus > 0 {
                affected.insert(locus - 1);
            }
            if locus + 1 < locus_count {
                affected.insert(locus);
            }
        }
        let mut fresh_drafts: BTreeMap<usize, SpineBoundaryDraft> = BTreeMap::new();
        for &boundary in &affected {
            let axis_window = &axis_slice[boundary..boundary + 2];
            let range_window = &ranges[boundary..boundary + 2];
            let links = genome::port_word_seams(axis_window, range_window, graph, &mut ports)?;
            let links = links
                .into_iter()
                .next()
                .ok_or_else(|| invalid("spine boundary link cardinality mismatch"))?;
            ensure(
                !links.is_empty(),
                "spine boundary {boundary} lost all legal links after splits",
            )?;
            let left_owners: Vec<u32> = ranges[boundary]
                .iter()
                .map(crate::traversal_owner)
                .collect();
            let right_owners: Vec<u32> = ranges[boundary + 1]
                .iter()
                .map(crate::traversal_owner)
                .collect();
            fresh_drafts.insert(
                boundary,
                build_spine_boundary_draft(
                    panel,
                    sources,
                    &flank_memo,
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
        // Extend the multiplicity backgrounds for every feature the fresh
        // loci and rebuilt boundary compositions charge that the initial
        // scan has not measured (split partials can produce new junction
        // seam subwalks; parents' features are already known).
        {
            let mut universe: BTreeSet<FeatureKey> = BTreeSet::new();
            for &locus in &fresh {
                for profile in &locus_classes[locus].profiles {
                    for key in profile.keys() {
                        universe.insert(key.clone());
                    }
                }
            }
            for draft in fresh_drafts.values() {
                for key in draft_feature_set(draft) {
                    universe.insert(key);
                }
            }
            let features: Vec<FeatureKey> = universe.into_iter().collect();
            backgrounds.scan_extend(panel, features, k, model)?;
        }
        for &locus in &fresh {
            folded[locus] = LocusFolded::build_backgrounds(
                &locus_classes[locus].profiles,
                &locus_classes[locus].class_owners,
                routed_equal,
                component_locus_to_partition,
                &backgrounds,
            )?;
            loss_tables[locus] = folded[locus]
                .loss_tables(model)
                .map(|(tables, _)| tables)?;
        }
        for (boundary, draft) in fresh_drafts {
            let left_owners: Vec<u32> =
                ranges[boundary].iter().map(crate::traversal_owner).collect();
            let right_owners: Vec<u32> =
                ranges[boundary + 1].iter().map(crate::traversal_owner).collect();
            let boundary_data = finalize_spine_boundary(
                draft,
                &ranges[boundary],
                &ranges[boundary + 1],
                &left_owners,
                &right_owners,
                routed_equal,
                component_locus_to_partition,
                model,
                &backgrounds,
            )?;
            boundaries[boundary] = boundary_data;
            successors_ref[boundary] = boundaries[boundary].successors.clone();
            seam_feature_sets[boundary] = boundary_feature_set(&boundaries[boundary].seam_by_pair);
        }
        viable = spine_viability(&ranges, &successors_ref);
        let seam_swings = seam_swings_of(locus_count, &successors_ref, &viable);
        for &locus in &fresh {
            sweeps[locus] = exhaustive_local_sweep(
                &folded[locus],
                &loss_tables[locus],
                &locus_classes[locus].membership,
                &viable[locus],
                Some(backbone_chain[locus]),
                model,
                seam_swings[locus],
                &locus_classes[locus].class_charges,
                &scorable_classes(&locus_classes[locus], &ranges[locus]),
            )?;
        }
        // Orphan shares at the fresh loci only (seam sets changed).
        for &locus in &fresh {
            let mut explainable: HashSet<FeatureKey> = HashSet::new();
            for profile in &locus_classes[locus].profiles {
                for key in profile.keys() {
                    explainable.insert(key.clone());
                }
            }
            if locus > 0 {
                for key in &seam_feature_sets[locus - 1] {
                    explainable.insert(key.clone());
                }
            }
            if locus + 1 < locus_count {
                for key in &seam_feature_sets[locus] {
                    explainable.insert(key.clone());
                }
            }
            let mut total = 0.0f64;
            let mut orphan = 0.0f64;
            for (key, &value) in &partition_obs[locus] {
                total += value;
                if !explainable.contains(key) {
                    orphan += value;
                }
            }
            orphan_shares[locus].1 = if total > 0.0 { orphan / total } else { 0.0 };
        }
        // Admissible-margin inputs recompute (min pair losses may have
        // dropped at the split loci; min seam scores at affected boundaries).
        min_pair = sweeps.iter().map(|s| s.min_viable_loss).collect();
        min_seam = successors_ref
            .iter()
            .map(|boundary| {
                boundary
                    .iter()
                    .flat_map(|list| list.iter().map(|&(_, score)| score))
                    .reduce(f64::min)
                    .unwrap_or(f64::INFINITY)
            })
            .collect();
        let recomputed = compute_margins(&min_pair, &min_seam)?;
        margins_admissible = recomputed.0;
        suffix_from_next = recomputed.1;
        margins = compute_local_margins(&min_pair, &successors_ref, &viable);
        rss_probe(rss, "spine_stage3_round")?;
    }
    let splits_seconds = splits_started.elapsed().as_secs_f64();

    // ------------------------------------------------ Stage 2: the exact chain DP.
    let dp_started = Instant::now();
    let membership_slices: Vec<Vec<usize>> = locus_classes
        .iter()
        .map(|classing| classing.membership.clone())
        .collect();
    let dp = run_exact_chain_dp(
        &ranges,
        &membership_slices,
        &sweeps,
        &successors_ref,
        &viable,
        &margins,
        &suffix_from_next,
        incumbent,
        conflict_window,
        &backbone_chain,
        rss,
    )?;
    let dp_seconds = dp_started.elapsed().as_secs_f64();
    eprintln!(
        "[spine] exact DP done: best {:.2}, states {:?}, {:.2}s",
        dp.best_score,
        dp.state_counts,
        dp_seconds
    );

    // ------------------------------------------------ Stage 4: final rescore.
    let rescore_started = Instant::now();
    let route = dp.route.clone();
    let seam_score = |boundary: usize, left: usize, right: usize| -> Option<f64> {
        successors_ref[boundary][left]
            .iter()
            .find(|&&(next, _)| next == right)
            .map(|&(_, score)| score)
    };
    let mut matchings: Vec<bool> = Vec::with_capacity(locus_count - 1);
    let mut identity_missing = 0usize;
    for boundary in 0..locus_count - 1 {
        let (l0, l1) = (route[boundary][0], route[boundary][1]);
        let (r0, r1) = (route[boundary + 1][0], route[boundary + 1][1]);
        let identity = seam_score(boundary, l0, r0)
            .zip(seam_score(boundary, l1, r1))
            .map(|(a, b)| a + b);
        let swap = seam_score(boundary, l0, r1)
            .zip(seam_score(boundary, l1, r0))
            .map(|(a, b)| a + b);
        let (use_swap, _total) = match (identity, swap) {
            (Some(i), Some(s)) => {
                if (s.total_cmp(&i)).is_lt() {
                    (true, s)
                } else {
                    (false, i)
                }
            }
            (Some(i), None) => (false, i),
            (None, Some(s)) => (true, s),
            (None, None) => return Err(invalid("selected route has an illegal boundary")),
        };
        if identity.is_none() {
            identity_missing += 1;
        }
        matchings.push(use_swap);
    }
    // Molecules: per-copy chains through the chosen matchings.
    let mut molecules: Vec<Vec<usize>> = vec![vec![route[0][0]], vec![route[0][1]]];
    for locus in 1..locus_count {
        for copy in 0..2 {
            let allele = if matchings[locus - 1] ^ (copy == 1) {
                route[locus][1]
            } else {
                route[locus][0]
            };
            molecules[copy].push(allele);
        }
    }
    // M1 oracle rescore: per-locus local pair terms plus per-boundary seam
    // terms, under the run's identity-matching convention (continuity with
    // the quarantined path's ladder) and under the best matching. The seam
    // terms are the boundary link scores themselves (the fixed-rule
    // charges: co-occurring junctions pooled bit-identically, novel
    // junctions paid by crossing reads).
    let mut oracle_memo: HashMap<String, Profile> = HashMap::new();
    let mut m1_local_total = 0.0f64;
    for locus in 0..locus_count {
        let first =
            oracle_allele_profile(panel, sources, &ranges[locus][route[locus][0]], &mut oracle_memo)?;
        let second =
            oracle_allele_profile(panel, sources, &ranges[locus][route[locus][1]], &mut oracle_memo)?;
        // Extend the backgrounds for any oracle-profile feature beyond the
        // geometric candidate universe before charging (memoized: the
        // scoring loop below re-hits the same profiles).
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
                let k = panel.syncmer_length_bp() as u64;
                backgrounds.scan_extend(panel, missing, k, model)?;
            }
        }
        let owner_first = crate::traversal_owner(&ranges[locus][route[locus][0]]);
        let owner_second = crate::traversal_owner(&ranges[locus][route[locus][1]]);
        if owner_first == owner_second {
            m1_local_total += merged_pair_loss_multiplicity(
                &first,
                &second,
                crate::owner_routed_obs(routed_equal, owner_first, component_locus_to_partition),
                model,
                &backgrounds,
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
                &backgrounds,
            )?;
        }
    }
    let mut m1_identity = m1_local_total;
    let mut m1_best = m1_local_total;
    for boundary in 0..locus_count - 1 {
        for copy in 0..2 {
            m1_identity += seam_score(boundary, route[boundary][copy], route[boundary + 1][copy])
                .ok_or_else(|| invalid("selected route boundary lacks a seam score"))?;
            m1_best += seam_score(
                boundary,
                molecules[copy][boundary],
                molecules[copy][boundary + 1],
            )
            .ok_or_else(|| invalid("selected molecule boundary lacks a seam score"))?;
        }
    }
    // Pooled external rescore of the two molecule sequences.
    let mut sequences = Vec::with_capacity(2);
    for copy in 0..2 {
        let mut sequence = Vec::new();
        for locus in 0..locus_count {
            sequence.extend(allele_sequence(
                sources,
                &ranges[locus][molecules[copy][locus]].segments,
            )?);
        }
        sequences.push(sequence);
    }
    let directory = rescore_directory.join("selected-pooled");
    let (pooled_loss, pooled_cost, pooled_runs) = genome::external_rescore(
        panel,
        [&sequences[0], &sequences[1]],
        sample_counts,
        model,
        &directory,
    )?;
    // References ladder (same-run scorer).
    let mut reference_scores = Vec::new();
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
            Some((span, path_of_source, &reference_partitions)),
            Some(&mut backgrounds),
            // The exact-chain-DP path is the DIPLOID pair ladder: panel
            // convention everywhere.
            None,
            // Group-shared material charged once (supervisor decision
            // 2026-09-23).
            true,
            None,
            None,
            None,
        )?;
        let mut reference_sequences = Vec::with_capacity(2);
        for route in [&route_a, &route_b] {
            let mut sequence = Vec::new();
            for segment in &route.segments {
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
        let directory = rescore_directory.join(format!("reference-{index}-{label}"));
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
            "m1_junction_census": reference_census,
            "pooled_external_rescore": pooled,
            "pooled_mem_queries": pooled_reference_cost.mem_queries,
            "pooled_initial_runs": pooled_reference_runs,
        }));
        rss_probe(rss, &format!("spine_reference_{label}"))?;
    }
    // Donor recovery per molecule (assessment side: the truth's non-target
    // sources define the tract; recovery is the selected molecules' donor
    // material, in donor coordinates).
    let merge_intervals = |intervals: &[(usize, u64, u64)]| -> Vec<(usize, u64, u64)> {
        let mut merged: Vec<(usize, u64, u64)> = Vec::new();
        for &(source, start, end) in intervals {
            match merged.last_mut() {
                Some(last) if last.0 == source && last.2 == start => last.2 = last.2.max(end),
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
    for copy in 0..2 {
        let raw: Vec<(usize, u64, u64)> = (0..locus_count)
            .flat_map(|locus| {
                ranges[locus][molecules[copy][locus]]
                    .segments
                    .iter()
                    .filter(|segment| donor_sources.contains(&segment.source))
                    .map(|segment| (segment.source, segment.start, segment.end))
            })
            .collect();
        let intervals = merge_intervals(&raw);
        let total_bases: u64 = intervals.iter().map(|&(_, start, end)| end - start).sum();
        let overlap: u64 = intervals
            .iter()
            .map(|&(_, start, end)| {
                truth_tract_intervals
                    .iter()
                    .map(|&(_, t_start, t_end)| end.max(t_start).saturating_sub(start.min(t_end)))
                    .sum::<u64>()
            })
            .sum();
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
    let rescore_seconds = rescore_started.elapsed().as_secs_f64();
    let reference_by_label = |label: &str| -> Option<f64> {
        reference_scores
            .iter()
            .find(|row| row["label"].as_str() == Some(label))
            .and_then(|row| row["m1_oracle_rescore"].as_f64())
    };
    let truth_m1 = reference_by_label("truth");
    let native2_m1 = reference_by_label("native2");
    let self_delta = m1_best - dp.best_score;

    let route_json: Vec<serde_json::Value> = (0..locus_count)
        .map(|locus| {
            serde_json::json!([
                traversal_json(&ranges[locus][route[locus][0]]),
                traversal_json(&ranges[locus][route[locus][1]]),
            ])
        })
        .collect();
    let molecules_json: Vec<Vec<serde_json::Value>> = (0..2)
        .map(|copy| {
            (0..locus_count)
                .map(|locus| traversal_json(&ranges[locus][molecules[copy][locus]]))
                .collect::<Vec<_>>()
        })
        .collect();

    Ok(serde_json::json!({
        "stage1_sweep": stage1_summary,
        "stage3_targeted_splits": {
            "rounds": split_rounds,
            "split_loci": (0..locus_count).filter(|&locus| split_done[locus]).collect::<Vec<_>>(),
            "split_admitted_counts": split_admitted_counts,
            "split_candidate_counts": split_candidate_counts,
            "total_split_candidates": split_candidate_counts.iter().sum::<usize>(),
            "wall_seconds": splits_seconds,
        },
        "stage2_exact_chain": {
            "conflict_window": conflict_window,
            "retention_margins_local": margins,
            "admissible_margins_global": margins_admissible,
            "viable_seam_extremes": seam_extremes(&successors_ref, &viable),
            "min_pair": min_pair,
            "min_seam": min_seam,
            "suffix_from_next": suffix_from_next,
            "retained_pairs": dp.retained_pairs,
            "state_counts": dp.state_counts,
            "transitions": dp.transitions,
            "merged_drops": dp.merged_drops,
            "bound_pruned": dp.bound_pruned,
            "margin_pruned": dp.margin_pruned,
            "best_score_internal": dp.best_score,
            "terminal_tie_count": dp.tie_count,
            "tie_pair_sets": dp.tie_pair_sets,
            "wall_seconds": dp_seconds,
        },
        "stage4_rescore": {
            "selected_route": route_json,
            "molecules": molecules_json,
            "boundary_matchings_swap": matchings,
            "identity_matching_missing_boundaries": identity_missing,
            "m1_oracle_rescore_identity": m1_identity,
            "m1_oracle_rescore_best_matching": m1_best,
            "selected_pooled_external_rescore": pooled_loss,
            "selected_pooled_mem_queries": pooled_cost.mem_queries,
            "selected_pooled_initial_runs": pooled_runs,
            "internal_vs_external_self_check": {
                "internal_geometric": dp.best_score,
                "external_m1_oracle": m1_best,
                "delta": self_delta,
                "delta_percent": if dp.best_score != 0.0 {
                    100.0 * self_delta / dp.best_score
                } else {
                    0.0
                },
            },
            "references": reference_scores,
            "ladder": {
                "truth_m1": truth_m1,
                "native2_m1": native2_m1,
                "selected_m1_identity": m1_identity,
                "selected_m1_best_matching": m1_best,
                "truth_minus_selected": truth_m1.map(|truth| truth - m1_identity),
                "selected_minus_native2": native2_m1.map(|native| m1_identity - native),
                "truth_minus_native2": match (truth_m1, native2_m1) {
                    (Some(truth), Some(native)) => Some(truth - native),
                    _ => None,
                },
            },
            "donor_recovery": donor_recovery,
            "truth_tract_intervals": truth_tract_intervals,
            "truth_tract_bases": truth_tract_bases,
            "wall_seconds": rescore_seconds,
        },
        "flags": {
            "exact_dp": true,
            "exactness_caveats": [
                "the DP is an exact min-plus over the margin-retained state space; \
                 the retention margin is the LOCAL seam-swing bound (pairs whose local \
                 loss exceeds the viable best by more than the total seam-evidence \
                 swing at the two adjacent boundaries are excluded by design, and the \
                 exclusion counts are reported); the all-native backbone pair is \
                 structurally retained at every locus", 
                "edge pruning is admissible against the native-backbone incumbent: any \
                 pruned edge provably lies in no chain at or below the incumbent, and \
                 the optimum is at or below it",
                "same-state equal-score merges keep one predecessor representative: the \
                 optimal VALUE and the terminal tie sets are exact; the pre-terminal route \
                 multiplicity of tied histories is folded, and reported via merged_drops",
                "reverse single-segment alleles share their forward twin's local profile \
                 (RC-canonical keys make the self view exact; the inverted-repeat part is \
                 inherited from the forward twin's extraction)",
            ],
            "no_tuning_constants": true,
            "guard_events": "none (the RSS guard fails closed; reaching this report means \
                             no guard event fired)",
        },
        "stage_wall_seconds": {
            "classing": classing_seconds,
            "boundaries": boundary_seconds,
            "sweep": sweep_seconds,
            "truth_decomposition": truth_seconds,
            "orphan_diagnostic": orphan_seconds,
            "targeted_splits": splits_seconds,
            "exact_chain_dp": dp_seconds,
            "rescore": rescore_seconds,
            "total": started.elapsed().as_secs_f64(),
        },
        "rss_peak_bytes": rss.peak_bytes(),
        "total_seconds": started.elapsed().as_secs_f64(),
    }))
}
