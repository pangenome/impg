//! Syng-seeded homology query with ends-only projection refinement.
//!
//! Pipeline per hop:
//!   1. `query_region_with_anchors` — syng returns padded syncmer-
//!      resolution intervals per (target, strand), with anchor positions.
//!      No K-cap on anchor collection; chain-level positional cap
//!      bounds memory downstream.
//!   2. `chain_anchors` — flatten per-visit anchors into PAF-equivalent
//!      exact-match records and delegate chain/scaffold filtering to
//!      SweepGA's plane-sweep scaffold filter. This keeps the off-diagonal
//!      short-hit handling in one battle-tested implementation.
//!   3. `dedupe_strand_overlaps` — per path, overlapping +/- chains
//!      resolved by majority anchor count.
//!   4. `refine_ends_only` — per chain, two ends-free BiWFA calls:
//!      backward from the first anchor to the query region start, and
//!      forward from the last anchor to the query region end. Each is
//!      bounded by `extend_budget` on the target axis and by the
//!      neighbor chain's target edge. The interior is trusted (every
//!      interior anchor is a syncmer match); we do not align through
//!      it. Real-CIGAR output would require per-anchor-pair interior
//!      gap-BiWFA and is left for when PAF-style output is needed.
//!   5. Multihop outer loop.

use std::cell::RefCell;
use std::io;

use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus, Distance, MemoryMode};
use log::warn;
use rustc_hash::FxHashSet;

use rayon::prelude::*;

use crate::graph::reverse_complement;
use crate::sequence_index::{SequenceIndex as _, UnifiedSequenceIndex};
use crate::syng::{
    Anchor, HomologousInterval, HomologousIntervalWithAnchors, SyngIndex, SyngSeedFilter,
};

/// Minimum fraction of query bases that must align at 'M'/'=' ops for
/// the refined homolog to be accepted. Below this, the "homolog" is
/// treated as syncmer noise and the refinement falls back to syng's
/// padded bounds. Permissive by default (0.3) — only filters obvious
/// non-homology.
const MIN_ALIGNMENT_IDENTITY: f64 = 0.3;

thread_local! {
    /// Edit-distance aligner for per-homolog refinement. `MemoryMode::Low`
    /// (piggyback BiWFA) benchmarked faster than `High` on our short
    /// per-projection windows despite the extra per-step allocations.
    static EDGE_ALIGNER: RefCell<Option<AffineWavefronts>> =
        const { RefCell::new(None) };
}

fn with_edge_aligner<F, R>(f: F) -> R
where
    F: FnOnce(&mut AffineWavefronts) -> R,
{
    EDGE_ALIGNER.with(|cell| {
        let mut opt = cell.borrow_mut();
        if opt.is_none() {
            *opt = Some(Distance::Edit.create_aligner(None, Some(&MemoryMode::Low)));
        }
        f(opt.as_mut().unwrap())
    })
}

/// Chain exact syncmer anchors with SweepGA's scaffold-aware plane sweep.
///
/// Syncmer anchors are represented as 1:1 exact PAF-style records and
/// passed to `sweepga::paf_filter::PafFilter`. The syng query code must
/// not maintain a parallel plane-sweep/scaffold implementation: SweepGA's
/// scaffold filter is where off-diagonal short-hit suppression belongs.
pub fn chain_anchors(
    hits: Vec<HomologousIntervalWithAnchors>,
    extend_budget: u64,
    syncmer_len: u64,
) -> Vec<HomologousIntervalWithAnchors> {
    let max_query_hop = extend_budget.saturating_mul(3).max(extend_budget);
    let drift_per_hop = extend_budget;
    let max_chain_span = extend_budget.saturating_mul(10).max(extend_budget);
    chain_anchors_with_limits(
        hits,
        syncmer_len,
        max_query_hop,
        drift_per_hop,
        max_chain_span,
    )
}

/// Collinear-chain shared syncmer anchors per `(target path, strand)` by
/// delegating to SweepGA's scaffold filter.
///
/// `max_query_hop` is passed as SweepGA's scaffold gap. The drift/span
/// parameters are retained in the signature for older callers but are no
/// longer interpreted locally; avoiding an impg-side clone of SweepGA's
/// scaffold algorithm is the point of this path.
pub fn chain_anchors_with_limits(
    hits: Vec<HomologousIntervalWithAnchors>,
    syncmer_len: u64,
    max_query_hop: u64,
    _drift_per_hop: u64,
    _max_chain_span: u64,
) -> Vec<HomologousIntervalWithAnchors> {
    chain_anchors_with_sweepga_scaffold_mass(hits, syncmer_len, max_query_hop, syncmer_len)
}

/// SweepGA-backed anchor chaining with an explicit scaffold mass.
///
/// `min_scaffold_length` is in bp of exact syncmer seed support, matching
/// SweepGA's scaffold-mass semantics. Callers that already know the desired
/// anchor-count threshold should pass `threshold * syncmer_len` here so
/// short off-diagonal seed runs are removed before interval refinement.
pub fn chain_anchors_with_sweepga_scaffold_mass(
    hits: Vec<HomologousIntervalWithAnchors>,
    syncmer_len: u64,
    scaffold_gap: u64,
    min_scaffold_length: u64,
) -> Vec<HomologousIntervalWithAnchors> {
    if hits.is_empty() {
        return hits;
    }

    use std::collections::HashMap;
    use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, RecordMeta, ScoringFunction};

    let syncmer_len = syncmer_len.max(1);
    let scaffold_gap = scaffold_gap.max(syncmer_len);
    let min_scaffold_length = min_scaffold_length.max(syncmer_len);
    let per_target: HashMap<(String, char), Vec<Anchor>> = hits
        .into_par_iter()
        .fold(HashMap::new, |mut buckets, hit| {
            buckets
                .entry((hit.genome, hit.strand))
                .or_insert_with(Vec::new)
                .extend(hit.anchors);
            buckets
        })
        .reduce(HashMap::new, |mut left, right| {
            for (key, anchors) in right {
                left.entry(key).or_insert_with(Vec::new).extend(anchors);
            }
            left
        });

    if per_target.is_empty() {
        return Vec::new();
    }

    let filter_config = FilterConfig {
        chain_gap: 0,
        min_block_length: syncmer_len,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::ManyToMany,
        scaffold_max_per_query: None,
        scaffold_max_per_target: None,
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: true,
        scaffold_gap,
        min_scaffold_length,
        scaffold_overlap_threshold: 0.5,
        scaffold_max_deviation: 0,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let mut out: Vec<HomologousIntervalWithAnchors> = per_target
        .into_par_iter()
        .flat_map(|((genome, strand), mut anchors)| {
        anchors.sort_by(|a, b| {
            a.query_pos
                .cmp(&b.query_pos)
                .then(a.target_pos.cmp(&b.target_pos))
        });
        anchors.dedup_by(|a, b| {
            a.query_pos == b.query_pos && a.target_pos == b.target_pos && a.node_id == b.node_id
        });
        if anchors.is_empty() {
            return Vec::new();
        }
        let records: Vec<RecordMeta> = anchors
            .iter()
            .enumerate()
            .map(|(rank, anchor)| RecordMeta {
                rank,
                query_name: "__syng_query#0#query".to_string(),
                target_name: genome.clone(),
                query_start: anchor.query_pos,
                query_end: anchor.query_pos.saturating_add(syncmer_len),
                target_start: anchor.target_pos,
                target_end: anchor.target_pos.saturating_add(syncmer_len),
                block_length: syncmer_len,
                identity: 1.0,
                matches: syncmer_len,
                alignment_length: syncmer_len,
                strand,
                chain_id: None,
                chain_status: sweepga::mapping::ChainStatus::Unassigned,
                discard: false,
                overlapped: false,
            })
            .collect();
        let passing = match PafFilter::new(filter_config.clone())
            .with_keep_self(true)
            .with_scaffolds_only(true)
            .apply_filters(records)
        {
            Ok(passing) => passing,
            Err(err) => {
                warn!("SweepGA scaffold filtering failed for syng anchors on {genome}/{strand}: {err}");
                return Vec::new();
            }
        };
        let mut grouped: HashMap<String, Vec<Anchor>> = HashMap::new();
        for (rank, meta) in passing {
            let Some(anchor) = anchors.get(rank).copied() else {
                continue;
            };
            let chain_id = meta.chain_id.unwrap_or_else(|| format!("rank_{rank}"));
            grouped.entry(chain_id).or_default().push(anchor);
        }
        let mut bucket_out = Vec::with_capacity(grouped.len());
        for (_chain_id, mut chain_anchors) in grouped {
            chain_anchors.sort_by(|a, b| {
                a.query_pos
                    .cmp(&b.query_pos)
                    .then(a.target_pos.cmp(&b.target_pos))
            });
            let start = chain_anchors.iter().map(|a| a.target_pos).min().unwrap_or(0);
            let end = chain_anchors
                .iter()
                .map(|a| a.target_pos.saturating_add(syncmer_len))
                .max()
                .unwrap_or(start);
            bucket_out.push(HomologousIntervalWithAnchors {
                genome: genome.clone(),
                start,
                end,
                strand,
                anchors: chain_anchors,
            });
        }
        bucket_out
    })
    .collect();
    out.sort_by(|a, b| {
        a.genome
            .cmp(&b.genome)
            .then(a.strand.cmp(&b.strand))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.anchors.len().cmp(&b.anchors.len()))
    });
    out
}

fn query_scaled_chain_limits(
    query_range_len: u64,
    extend_budget: u64,
    syncmer_len: u64,
    merge_distance: i32,
) -> (u64, u64, u64) {
    let merge_gap = merge_distance.max(0) as u64;
    let chain_slop = extend_budget.max(merge_gap);
    let chain_gap = query_range_len
        .saturating_add(chain_slop)
        .max(extend_budget.saturating_mul(3))
        .max(syncmer_len);
    let drift_budget = query_range_len
        .saturating_add(chain_slop)
        .max(chain_slop)
        .max(syncmer_len);
    // Boundary realignment budget is not an SV-size budget. Keep the target
    // span cap query-scaled so large insertions, deletions, and local copy
    // changes can survive chaining. The user merge distance is also part of
    // the same semantic gap model, so it participates here before refinement
    // and before the next transitive hop.
    let sv_span_slop = query_range_len
        .saturating_mul(3)
        .max(50_000)
        .max(extend_budget.saturating_mul(10))
        .max(merge_gap);
    let target_span_cap = query_range_len
        .saturating_add(sv_span_slop)
        .max(syncmer_len);
    (chain_gap, drift_budget, target_span_cap)
}

/// Linear projection from an anchor to a query position — coordinate
/// arithmetic only. Used as a fallback when pairwise alignment is
/// unavailable.
///
/// * `'+'` strand: `t_anchor + (q_pos - q_anchor)`
/// * `'-'` strand: `t_anchor + syncmer_len - (q_pos - q_anchor)`
///
/// Clamped to `[0, target_len]`.
fn project_query_to_target(
    q_pos: u64,
    anchor: Anchor,
    strand: char,
    syncmer_len: u64,
    target_len: u64,
) -> u64 {
    let q_anchor = anchor.query_pos;
    let t_anchor = anchor.target_pos;
    let projected: i128 = match strand {
        '+' => t_anchor as i128 + (q_pos as i128 - q_anchor as i128),
        '-' => t_anchor as i128 + syncmer_len as i128 - (q_pos as i128 - q_anchor as i128),
        _ => t_anchor as i128,
    };
    if projected < 0 {
        0
    } else if projected as u128 > target_len as u128 {
        target_len
    } else {
        projected as u64
    }
}

/// Linear-projection fallback for one merged homolog. Used when the
/// pairwise alignment path is unavailable (no sequence index, fetch
/// failure, alignment failure).
fn refine_by_linear_projection(
    homolog: &HomologousIntervalWithAnchors,
    query_start: u64,
    query_end: u64,
    syncmer_len: u64,
    target_len: u64,
) -> (u64, u64) {
    if homolog.anchors.is_empty() {
        return (homolog.start, homolog.end);
    }
    let mut anchors = homolog.anchors.clone();
    anchors.sort_by_key(|a| a.query_pos);
    let leftmost = anchors[0];
    let rightmost = *anchors.last().unwrap();
    let t_for_qs = project_query_to_target(
        query_start,
        leftmost,
        homolog.strand,
        syncmer_len,
        target_len,
    );
    let t_for_qe = project_query_to_target(
        query_end,
        rightmost,
        homolog.strand,
        syncmer_len,
        target_len,
    );
    let (start, end) = if t_for_qs <= t_for_qe {
        (t_for_qs, t_for_qe)
    } else {
        (t_for_qe, t_for_qs)
    };
    (start.max(homolog.start), end.min(homolog.end))
}

/// Align one end of a chain with a pinned sub-query and one-sided
/// ends-free target. Returns `(skip_outer_bp, cigar)` on success.
///
/// * `outer_end` = `Outer::Left` → text prefix is free (the target
///   window's outer edge is its start).
/// * `outer_end` = `Outer::Right` → text suffix is free.
///
/// The sub-query is pinned End2End; the aligner can only extend or
/// truncate the target on `outer_end`. On completion, `skip_outer_bp`
/// is the leading (Left) or trailing (Right) `I` count, i.e. how
/// many target bp at the outer edge were not covered by the
/// alignment. Caller translates that back into a forward-strand
/// refined boundary.
#[derive(Copy, Clone)]
enum Outer {
    Left,
    Right,
}

#[allow(clippy::too_many_arguments)]
use std::sync::atomic::{AtomicU64, Ordering};
static PROF_FETCH_NS: AtomicU64 = AtomicU64::new(0);
static PROF_ALIGN_NS: AtomicU64 = AtomicU64::new(0);
static PROF_FETCH_COUNT: AtomicU64 = AtomicU64::new(0);
static PROF_ALIGN_COUNT: AtomicU64 = AtomicU64::new(0);

fn project_end_align(
    sub_query: &[u8],
    target_fwd_start: u64,
    target_fwd_end: u64,
    strand: char,
    genome: &str,
    outer: Outer,
    sequence_index: &UnifiedSequenceIndex,
    want_cigar: bool,
) -> Option<(u64, Option<Vec<u8>>)> {
    if sub_query.is_empty() || target_fwd_end <= target_fwd_start {
        return None;
    }
    let t_fetch = std::time::Instant::now();
    let t_bytes_fwd = sequence_index
        .fetch_sequence(genome, target_fwd_start as i32, target_fwd_end as i32)
        .ok()?;
    PROF_FETCH_NS.fetch_add(t_fetch.elapsed().as_nanos() as u64, Ordering::Relaxed);
    PROF_FETCH_COUNT.fetch_add(1, Ordering::Relaxed);
    if t_bytes_fwd.is_empty() {
        return None;
    }
    let t_bytes_oriented = if strand == '-' {
        reverse_complement(&t_bytes_fwd)
    } else {
        t_bytes_fwd
    };
    let t_len = t_bytes_oriented.len() as u64;
    // Outer end in oriented view = outer end in forward view for '+',
    // flipped for '-'. The BiWFA call is always on the oriented
    // sequence, so we map the desired forward outer to the oriented
    // outer here.
    let oriented_outer = match (strand, outer) {
        ('+', Outer::Left) | ('-', Outer::Right) => Outer::Left,
        ('+', Outer::Right) | ('-', Outer::Left) => Outer::Right,
        _ => return None,
    };
    let (text_begin_free, text_end_free) = match oriented_outer {
        Outer::Left => (t_len as i32, 0i32),
        Outer::Right => (0i32, t_len as i32),
    };

    let t_align = std::time::Instant::now();
    let res = with_edge_aligner(|aligner| {
        unsafe {
            lib_wfa2::bindings::wfa::wavefront_aligner_set_alignment_free_ends(
                aligner.aligner_mut(),
                0,
                0,
                text_begin_free,
                text_end_free,
            );
        }
        let status = aligner.align(sub_query, &t_bytes_oriented);
        if !matches!(status, AlignmentStatus::Completed) {
            return None;
        }
        let cigar = aligner.cigar();
        let leading_i = cigar.iter().take_while(|&&op| op == b'I').count();
        let trailing_i = cigar.iter().rev().take_while(|&&op| op == b'I').count();
        let core_end = cigar.len().saturating_sub(trailing_i);
        let matches = if core_end > leading_i {
            cigar[leading_i..core_end]
                .iter()
                .filter(|&&op| op == b'M' || op == b'=')
                .count()
        } else {
            0
        };
        let denom = sub_query.len();
        if denom == 0 || (matches as f64) / (denom as f64) < MIN_ALIGNMENT_IDENTITY {
            return None;
        }
        let skip_oriented_outer_bp = match oriented_outer {
            Outer::Left => leading_i as u64,
            Outer::Right => trailing_i as u64,
        };
        if skip_oriented_outer_bp >= t_len {
            return None;
        }
        Some((
            skip_oriented_outer_bp,
            if want_cigar {
                Some(cigar.to_vec())
            } else {
                None
            },
        ))
    });
    PROF_ALIGN_NS.fetch_add(t_align.elapsed().as_nanos() as u64, Ordering::Relaxed);
    PROF_ALIGN_COUNT.fetch_add(1, Ordering::Relaxed);
    res
}

/// Ends-only projection of one chain: two small ends-free BiWFAs,
/// one from the chain's first anchor backward to the query region
/// start, one from the last anchor forward to the query region end.
/// Each is bounded on the target axis by `extend_budget` and by the
/// neighbor chain's target edge on the same `(genome, strand)`.
/// Returns refined forward-strand `(ts, te)`. Interior anchors are
/// trusted — no CIGAR is emitted. Real CIGAR output would require
/// per-anchor-pair interior BiWFA and is left as a future extension.
#[allow(clippy::too_many_arguments)]
fn refine_ends_only(
    query_bytes_full: &[u8],
    query_region_start: u64,
    query_region_end: u64,
    extend_budget: u64,
    hit: &HomologousIntervalWithAnchors,
    prev_target_end: u64,
    next_target_start: u64,
    sequence_index: &UnifiedSequenceIndex,
    target_len: u64,
    syncmer_len: u64,
) -> Option<(u64, u64)> {
    if hit.anchors.is_empty() {
        return None;
    }
    let mut a_first = hit.anchors[0];
    let mut a_last = hit.anchors[0];
    for &a in &hit.anchors {
        if a.query_pos < a_first.query_pos {
            a_first = a;
        }
        if a.query_pos > a_last.query_pos {
            a_last = a;
        }
    }

    // Sub-query slice helper (clip to [query_region_start, query_region_end]
    // then index into query_bytes_full).
    let slice_sub_q = |lo: u64, hi: u64| -> Option<&[u8]> {
        let lo = lo.max(query_region_start);
        let hi = hi.min(query_region_end);
        if hi <= lo {
            return None;
        }
        let a = (lo - query_region_start) as usize;
        let b = (hi - query_region_start) as usize;
        if b > query_bytes_full.len() {
            return None;
        }
        Some(&query_bytes_full[a..b])
    };

    // Target-span computation is a function of how far the query
    // position sits from the anchor (`query_gap`). Target extends
    // just that distance plus a slop budget. The first attempt uses
    // a small "homologous short-circuit" slop (~5% of gap); if that
    // alignment fails we retry with the full `extend_budget`. Both
    // caps are per-projection and independent of tile size.
    // Target window = min(query_gap + slop, query_gap + extend_budget,
    // 2 × extend_budget). Caps target-side BiWFA work per call.
    let target_span_for_slop = |query_gap: u64, slop: u64| -> u64 {
        query_gap
            .saturating_add(slop)
            .min(extend_budget.saturating_mul(2))
            .max(syncmer_len)
    };
    let small_slop = |query_gap: u64| -> u64 { (query_gap / 20).max(32).min(extend_budget) };
    // Parallel-walk projection: for each edge, do BiWFA on at most
    // `extend_budget` bp of query closest to the anchor (where drift is
    // informative and alignment is cheap), then linearly extrapolate
    // the remainder of the gap along the anchor's diagonal. BiWFA
    // cost per call is bounded; long edges get the BiWFA refinement at
    // the near edge plus colinear extension for the far reach.
    //
    // Tiny gaps (< SKIP_BIWFA_MIN_GAP) skip BiWFA entirely — linear
    // drift over <64 bp is negligible.
    const SKIP_BIWFA_MIN_GAP: u64 = 64;
    let back_query_gap = a_first.query_pos.saturating_sub(query_region_start);
    let back_outer = match hit.strand {
        '+' => Outer::Left,
        '-' => Outer::Right,
        _ => return None,
    };
    // BiWFA covers the last `biwfa_back_gap` bp before A_first; the
    // remaining `linear_back_ext` bp are linear-extrapolated to qs.
    let biwfa_back_gap = back_query_gap.min(extend_budget);
    let linear_back_ext = back_query_gap - biwfa_back_gap;
    let back_biwfa_q_lo = a_first
        .query_pos
        .saturating_sub(biwfa_back_gap)
        .max(query_region_start);
    let back_biwfa_q_hi = (a_first.query_pos + syncmer_len).min(query_region_end);
    // Shift backward_bounds to be centered on biwfa_back_gap rather
    // than the full back_query_gap.
    let backward_bounds_bwg = |span: u64| -> (u64, u64) {
        match hit.strand {
            '+' => (
                a_first.target_pos.saturating_sub(span).max(prev_target_end),
                (a_first.target_pos + syncmer_len).min(target_len),
            ),
            '-' => (
                a_first.target_pos.min(target_len),
                (a_first.target_pos + syncmer_len + span)
                    .min(next_target_start)
                    .min(target_len),
            ),
            _ => (0, 0),
        }
    };
    let linear_back_from = |anchor_t: u64| -> u64 {
        // Linearly extend anchor_t by linear_back_ext along the
        // a_first diagonal, then clip to neighbor/target bounds.
        match hit.strand {
            '+' => anchor_t
                .saturating_sub(linear_back_ext)
                .max(prev_target_end),
            '-' => (anchor_t + linear_back_ext)
                .min(next_target_start)
                .min(target_len),
            _ => anchor_t,
        }
    };
    let pure_linear_back = || -> u64 {
        let t = project_query_to_target(
            query_region_start,
            a_first,
            hit.strand,
            syncmer_len,
            target_len,
        );
        match hit.strand {
            '+' => t.max(prev_target_end),
            '-' => t.min(next_target_start),
            _ => t,
        }
    };
    let try_back = |slop: u64| -> Option<u64> {
        if back_query_gap < SKIP_BIWFA_MIN_GAP {
            return Some(pure_linear_back());
        }
        let span = target_span_for_slop(biwfa_back_gap, slop);
        let (ts, te) = backward_bounds_bwg(span);
        if te <= ts {
            return None;
        }
        let sq = slice_sub_q(back_biwfa_q_lo, back_biwfa_q_hi)?;
        let (skip, _cigar) = project_end_align(
            sq,
            ts,
            te,
            hit.strand,
            &hit.genome,
            back_outer,
            sequence_index,
            false,
        )?;
        let t_at_biwfa_lo = match hit.strand {
            '+' => ts + skip,
            '-' => te.saturating_sub(skip),
            _ => return None,
        };
        Some(linear_back_from(t_at_biwfa_lo))
    };
    let back_result = try_back(small_slop(biwfa_back_gap)).or_else(|| try_back(extend_budget));

    // ── Forward projection ──
    let fwd_query_gap =
        query_region_end.saturating_sub((a_last.query_pos + syncmer_len).min(query_region_end));
    let fwd_outer = match hit.strand {
        '+' => Outer::Right,
        '-' => Outer::Left,
        _ => return None,
    };
    let biwfa_fwd_gap = fwd_query_gap.min(extend_budget);
    let linear_fwd_ext = fwd_query_gap - biwfa_fwd_gap;
    let fwd_biwfa_q_lo = a_last.query_pos.max(query_region_start);
    let fwd_biwfa_q_hi = (a_last.query_pos + syncmer_len + biwfa_fwd_gap).min(query_region_end);
    let forward_bounds_bwg = |span: u64| -> (u64, u64) {
        match hit.strand {
            '+' => (
                a_last.target_pos.min(target_len),
                (a_last.target_pos + syncmer_len + span)
                    .min(next_target_start)
                    .min(target_len),
            ),
            '-' => (
                a_last.target_pos.saturating_sub(span).max(prev_target_end),
                (a_last.target_pos + syncmer_len).min(target_len),
            ),
            _ => (0, 0),
        }
    };
    let linear_fwd_from = |anchor_t: u64| -> u64 {
        match hit.strand {
            '+' => (anchor_t + linear_fwd_ext)
                .min(next_target_start)
                .min(target_len),
            '-' => anchor_t.saturating_sub(linear_fwd_ext).max(prev_target_end),
            _ => anchor_t,
        }
    };
    let pure_linear_fwd = || -> u64 {
        let t = project_query_to_target(
            query_region_end,
            a_last,
            hit.strand,
            syncmer_len,
            target_len,
        );
        match hit.strand {
            '+' => t.min(next_target_start),
            '-' => t.max(prev_target_end),
            _ => t,
        }
    };
    let try_fwd = |slop: u64| -> Option<u64> {
        if fwd_query_gap < SKIP_BIWFA_MIN_GAP {
            return Some(pure_linear_fwd());
        }
        let span = target_span_for_slop(biwfa_fwd_gap, slop);
        let (ts, te) = forward_bounds_bwg(span);
        if te <= ts {
            return None;
        }
        let sq = slice_sub_q(fwd_biwfa_q_lo, fwd_biwfa_q_hi)?;
        let (skip, _cigar) = project_end_align(
            sq,
            ts,
            te,
            hit.strand,
            &hit.genome,
            fwd_outer,
            sequence_index,
            false,
        )?;
        let t_at_biwfa_hi = match hit.strand {
            '+' => te.saturating_sub(skip),
            '-' => ts + skip,
            _ => return None,
        };
        Some(linear_fwd_from(t_at_biwfa_hi))
    };
    let fwd_result = try_fwd(small_slop(biwfa_fwd_gap)).or_else(|| try_fwd(extend_budget));

    // `back_result` / `fwd_result` already hold the refined forward-
    // strand target position for whichever tier succeeded. Compose
    // into (ts, te) per strand.
    let (refined_ts, refined_te) = match hit.strand {
        '+' => (
            back_result.unwrap_or(a_first.target_pos),
            fwd_result.unwrap_or((a_last.target_pos + syncmer_len).min(target_len)),
        ),
        '-' => (
            fwd_result.unwrap_or(a_last.target_pos),
            back_result.unwrap_or((a_first.target_pos + syncmer_len).min(target_len)),
        ),
        _ => return None,
    };

    if refined_te <= refined_ts {
        return None;
    }

    Some((refined_ts, refined_te))
}

/// Per-path strand dedupe: if a `+` and a `-` interval on the same
/// target path overlap on forward-strand coordinates, keep the majority-
/// anchor strand; drop the minority as random-syncmer-coincidence noise.
/// Non-overlapping `+`/`-` on the same path stay separate (real
/// inversion on a distinct region of the contig).
///
/// Bucket by genome first so the O(K²) cross-strand check runs per
/// target only, not across all genomes. Without this, rDNA queries
/// that emit millions of chains would block for hours on the global
/// O(N²) compare.
fn dedupe_strand_overlaps(
    hits: Vec<HomologousIntervalWithAnchors>,
) -> Vec<HomologousIntervalWithAnchors> {
    if hits.len() <= 1 {
        return hits;
    }
    use std::collections::HashMap;
    let mut by_genome: HashMap<String, Vec<HomologousIntervalWithAnchors>> = HashMap::new();
    for h in hits {
        by_genome.entry(h.genome.clone()).or_default().push(h);
    }
    let mut out: Vec<HomologousIntervalWithAnchors> = Vec::new();
    for (_, group) in by_genome {
        if group.len() <= 1 {
            out.extend(group);
            continue;
        }
        let n = group.len();
        let mut keep = vec![true; n];
        for i in 0..n {
            if !keep[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !keep[j] {
                    continue;
                }
                let a = &group[i];
                let b = &group[j];
                // Same-genome invariant by construction; just check
                // cross-strand overlap.
                if a.strand == b.strand || a.start >= b.end || b.start >= a.end {
                    continue;
                }
                if a.anchors.len() >= b.anchors.len() {
                    keep[j] = false;
                } else {
                    keep[i] = false;
                    break;
                }
            }
        }
        out.extend(
            group
                .into_iter()
                .zip(keep)
                .filter_map(|(iv, k)| if k { Some(iv) } else { None }),
        );
    }
    out
}

/// Default cap on the adaptive minimum anchor count per chain for emission.
///
/// The actual threshold is derived from query length, syncmer density, and an
/// assumed identity in [`effective_min_chain_anchors_for_syncmer`], then capped
/// here. A modest cap avoids dropping legitimate CNV/SV loci whose SweepGA
/// scaffold chains are intentionally local rather than one locus-spanning
/// alignment.
pub const DEFAULT_MIN_CHAIN_ANCHORS: usize = 20;

/// Default identity used to estimate how many exact syncmers should survive in
/// a true homologous chain.
pub const DEFAULT_MIN_CHAIN_IDENTITY: f64 = 0.95;

/// Fraction of expected exact shared syncmers required for emission.
///
/// We do not require the full expectation because SNP clustering, local
/// divergence, structural variation, and syncmer-selection edge effects all
/// lower observed anchor counts. The value is high enough to reject weak
/// paralog/repeat chains in AMY-scale loci while preserving short-query
/// sensitivity.
pub const DEFAULT_MIN_CHAIN_EXPECTED_FRACTION: f64 = 0.10;

/// Effective syng chain anchor threshold for a specific query span.
///
/// `requested_min` is treated as a user/CLI cap. The returned value scales with
/// the query span and the expected exact-syncmer survival at 95% identity. `0`
/// disables the anchor-count filter.
pub fn effective_min_chain_anchors_for_syncmer(
    query_range_len: u64,
    syncmer_len: u64,
    smer_len: u64,
    requested_min: usize,
) -> usize {
    if requested_min == 0 {
        return 0;
    }

    let density = closed_syncmer_density(syncmer_len, smer_len);
    let exact_survival = DEFAULT_MIN_CHAIN_IDENTITY.powf(syncmer_len as f64);
    let expected_shared =
        query_range_len as f64 * density * exact_survival * DEFAULT_MIN_CHAIN_EXPECTED_FRACTION;
    let required = expected_shared.ceil() as usize;
    requested_min.min(required.max(1))
}

pub fn effective_min_chain_anchors(
    query_range_len: u64,
    syncmer_len: u64,
    requested_min: usize,
) -> usize {
    effective_min_chain_anchors_for_syncmer(query_range_len, syncmer_len, 8, requested_min)
}

fn closed_syncmer_density(syncmer_len: u64, smer_len: u64) -> f64 {
    if syncmer_len == 0 || smer_len == 0 || smer_len > syncmer_len {
        return 1.0;
    }
    let window_count = syncmer_len - smer_len + 1;
    2.0 / (window_count as f64 + 1.0)
}

/// Run one hop of syng-seeded homology query.
pub fn one_hop(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    one_hop_ext(
        syng_index,
        query_name,
        query_start,
        query_end,
        padding,
        0,
        DEFAULT_EXTEND_BUDGET_BP,
        DEFAULT_MIN_CHAIN_ANCHORS,
        0.0,
        sequence_index,
    )
}

/// [`one_hop`] with tunable source-side widening and per-chain
/// extension budget.
#[allow(clippy::too_many_arguments)]
pub fn one_hop_ext(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    query_extension: u64,
    extend_budget: u64,
    min_chain_anchors: usize,
    min_chain_fraction: f64,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    one_hop_ext_visited(
        syng_index,
        query_name,
        query_start,
        query_end,
        padding,
        query_extension,
        extend_budget,
        min_chain_anchors,
        min_chain_fraction,
        None,
        sequence_index,
    )
}

/// Core single-hop entry: chain anchors via mutual-best-buddy, then
/// ends-only projection per chain.
///
/// Chains are dropped at emission when either:
///   - anchor count is below the query-scaled effective `min_chain_anchors`
///   - query-extent / query_range_len < `min_chain_fraction`
///     (default 0.0 in the library helper: no length filter)
#[allow(clippy::too_many_arguments)]
pub fn one_hop_ext_visited(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    query_extension: u64,
    extend_budget: u64,
    min_chain_anchors: usize,
    min_chain_fraction: f64,
    visited_nodes: Option<&mut FxHashSet<u32>>,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    one_hop_ext_visited_with_seed_filter(
        syng_index,
        query_name,
        query_start,
        query_end,
        padding,
        query_extension,
        extend_budget,
        min_chain_anchors,
        min_chain_fraction,
        visited_nodes,
        SyngSeedFilter::default(),
        -1,
        sequence_index,
    )
}

/// [`one_hop_ext_visited`] with a high-frequency syncmer seed filter.
#[allow(clippy::too_many_arguments)]
pub fn one_hop_ext_visited_with_seed_filter(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    query_extension: u64,
    extend_budget: u64,
    min_chain_anchors: usize,
    min_chain_fraction: f64,
    visited_nodes: Option<&mut FxHashSet<u32>>,
    seed_filter: SyngSeedFilter,
    merge_distance: i32,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    use std::time::Instant;
    let prof = std::env::var("SYNG_PROFILE").is_ok();
    let t0 = Instant::now();
    let hits = syng_index.query_region_with_anchors_ext_visited_seed_filtered(
        query_name,
        query_start,
        query_end,
        padding,
        query_extension,
        visited_nodes,
        seed_filter,
    )?;
    if prof {
        eprintln!("PROF emit={:?} hits={}", t0.elapsed(), hits.len());
    }
    let syncmer_len = (syng_index.params.w + syng_index.params.k) as u64;
    let query_range_len = query_end.saturating_sub(query_start);
    let t1 = Instant::now();
    let total_anchors: usize = hits.iter().map(|h| h.anchors.len()).sum();
    let (chain_gap, drift_budget, target_span_cap) =
        query_scaled_chain_limits(query_range_len, extend_budget, syncmer_len, merge_distance);
    // Filter chains by a length-scaled anchor-count model. The default
    // estimates how many exact syncmers should survive at high identity and
    // requires a fraction of that, so short queries stay sensitive while
    // long repeat-dense loci need multiple bounded GBWT-walk seed islands.
    let effective_min = effective_min_chain_anchors_for_syncmer(
        query_range_len,
        syncmer_len,
        syng_index.params.k as u64,
        min_chain_anchors,
    );
    let hits = chain_anchors_with_sweepga_scaffold_mass(
        hits,
        syncmer_len,
        chain_gap,
        (effective_min as u64).saturating_mul(syncmer_len),
    );
    let pre_filter = hits.len();
    // Chain query-extent = a_last.query_pos + syncmer_len - a_first.query_pos.
    // Threshold in bp; 0 means "no extent filter".
    let min_chain_extent_bp = (query_range_len as f64 * min_chain_fraction.max(0.0)) as u64;
    let hits: Vec<HomologousIntervalWithAnchors> = hits
        .into_iter()
        .filter(|h| {
            if h.anchors.len() < effective_min {
                return false;
            }
            if min_chain_extent_bp == 0 {
                return true;
            }
            let mut qmin = u64::MAX;
            let mut qmax = 0u64;
            for a in &h.anchors {
                qmin = qmin.min(a.query_pos);
                qmax = qmax.max(a.query_pos);
            }
            let extent = qmax.saturating_sub(qmin).saturating_add(syncmer_len);
            extent >= min_chain_extent_bp
        })
        .collect();
    if prof {
        eprintln!(
            "PROF chain_anchors={:?} total_anchors={} chains={} after_filter={} (min_anchors={} eff={} min_frac={} min_bp={} chain_gap={} drift={} span_cap={})",
            t1.elapsed(),
            total_anchors,
            pre_filter,
            hits.len(),
            min_chain_anchors,
            effective_min,
            min_chain_fraction,
            min_chain_extent_bp,
            chain_gap,
            drift_budget,
            target_span_cap,
        );
    }
    let t2 = Instant::now();
    let hits = dedupe_strand_overlaps(hits);
    if prof {
        eprintln!("PROF dedupe={:?} after={}", t2.elapsed(), hits.len());
    }

    // Fetch the full query bytes ONCE per hop; projection extracts
    // the sub-ranges its anchors + extension budget cover.
    let t3 = Instant::now();
    let query_bytes_cache: Option<Vec<u8>> = sequence_index
        .fetch_sequence(query_name, query_start as i32, query_end as i32)
        .ok();
    if prof {
        eprintln!("PROF fetch_query={:?}", t3.elapsed());
    }

    // Attach per-chain target-axis neighbors (prev_end, next_start)
    // within each (genome, strand) bucket. These bound extension.
    let t4 = Instant::now();
    let hits_with_neighbors = attach_target_neighbors(hits);
    if prof {
        eprintln!("PROF neighbors={:?}", t4.elapsed());
    }

    // Projection is independent per chain; drive it on rayon's global
    // work-stealing pool. Singleton chains (one anchor) skip BiWFA —
    // there's nothing to refine between two unanchored ends, so
    // linear projection gives the same answer for free. Multi-anchor
    // chains get two ends-free BiWFA calls, one per extension side.
    //
    // Thread-local `EDGE_ALIGNER` means each rayon worker reuses its
    // own aligner instance — no cross-thread sharing. The global
    // rayon pool is sized to `num_cpus::get()` by default (overridable
    // via `RAYON_NUM_THREADS`), so CPU use stays bounded regardless of
    // how many chains this hop emits.
    let qb_slice = query_bytes_cache.as_deref();
    let t5 = Instant::now();
    let singleton_count = hits_with_neighbors
        .iter()
        .filter(|(h, _, _)| h.anchors.len() < 2)
        .count();
    let out: Vec<HomologousInterval> = hits_with_neighbors
        .par_iter()
        .with_min_len(128)
        .filter_map(|(hit, prev_end, next_start)| {
            let target_len = sequence_index
                .get_sequence_length(&hit.genome)
                .map(|l| l as u64)
                .unwrap_or(u64::MAX);

            let (refined_start, refined_end) = if hit.anchors.len() < 2 {
                refine_by_linear_projection(hit, query_start, query_end, syncmer_len, target_len)
            } else {
                let projected = qb_slice.and_then(|qb| {
                    refine_ends_only(
                        qb,
                        query_start,
                        query_end,
                        extend_budget,
                        hit,
                        *prev_end,
                        *next_start,
                        sequence_index,
                        target_len,
                        syncmer_len,
                    )
                });
                projected.unwrap_or_else(|| {
                    refine_by_linear_projection(
                        hit,
                        query_start,
                        query_end,
                        syncmer_len,
                        target_len,
                    )
                })
            };
            if refined_end <= refined_start {
                return None;
            }
            Some(HomologousInterval {
                genome: hit.genome.clone(),
                start: refined_start,
                end: refined_end,
                strand: hit.strand,
                cigar: None,
            })
        })
        .collect();
    if prof {
        let fetch_ns = PROF_FETCH_NS.swap(0, Ordering::Relaxed);
        let align_ns = PROF_ALIGN_NS.swap(0, Ordering::Relaxed);
        let fetch_n = PROF_FETCH_COUNT.swap(0, Ordering::Relaxed);
        let align_n = PROF_ALIGN_COUNT.swap(0, Ordering::Relaxed);
        eprintln!(
            "PROF refine={:?} chains={} singletons={} out={} fetch_sum={:.2}s (n={}) align_sum={:.2}s (n={})",
            t5.elapsed(),
            hits_with_neighbors.len(),
            singleton_count,
            out.len(),
            fetch_ns as f64 / 1e9,
            fetch_n,
            align_ns as f64 / 1e9,
            align_n,
        );
    }
    Ok(out)
}

/// For each cluster, attach the `(prev_target_end, next_target_start)`
/// of its signature-sorted neighbors within the same
/// `(genome, strand)`. A cluster with no predecessor gets `prev_end = 0`;
/// with no successor, `next_start = u64::MAX`.
fn attach_target_neighbors(
    hits: Vec<HomologousIntervalWithAnchors>,
) -> Vec<(HomologousIntervalWithAnchors, u64, u64)> {
    use std::collections::HashMap;
    let mut groups: HashMap<(String, char), Vec<HomologousIntervalWithAnchors>> = HashMap::new();
    for h in hits {
        groups
            .entry((h.genome.clone(), h.strand))
            .or_default()
            .push(h);
    }
    let mut out: Vec<(HomologousIntervalWithAnchors, u64, u64)> = Vec::new();
    for (_, mut group) in groups {
        group.sort_by_key(|h| h.start);
        let bounds: Vec<(u64, u64)> = group.iter().map(|h| (h.start, h.end)).collect();
        let n = group.len();
        for (i, h) in group.into_iter().enumerate() {
            let prev_end = if i > 0 { bounds[i - 1].1 } else { 0 };
            let next_start = if i + 1 < n { bounds[i + 1].0 } else { u64::MAX };
            out.push((h, prev_end, next_start));
        }
    }
    out
}

/// Default sub-query/target extension budget used by the thin
/// [`query_transitive`] / [`query_transitive_ext`] wrappers when no
/// explicit value is passed. Matches the default on the CLI's
/// `--syng-extend-budget`.
pub const DEFAULT_EXTEND_BUDGET_BP: u64 = 1_000;

/// Run a multihop syng-seeded homology query.
pub fn query_transitive(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    max_depth: u16,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    query_transitive_ext(
        syng_index,
        query_name,
        query_start,
        query_end,
        padding,
        max_depth,
        0,
        DEFAULT_EXTEND_BUDGET_BP,
        DEFAULT_MIN_CHAIN_ANCHORS,
        0.0,
        sequence_index,
    )
}

/// [`query_transitive`] with extension tunables plumbed through every hop.
///
/// * `query_extension`: widens syncmer discovery on the source side at
///   each frontier step, so BFS can follow conserved blocks whose
///   endpoints fall just outside each frontier region.
/// * `extend_budget`: per-chain bp budget for ends-free target-side
///   extension in projection; see [`one_hop_ext_visited`].
/// * `min_chain_fraction`: drop chains whose query-extent is less than
///   `fraction × query_range_len`. 0.0 = no filter.
#[allow(clippy::too_many_arguments)]
/// Multi-hop BFS that accumulates anchor-bearing chains across all
/// discovered regions. Unlike [`query_transitive_ext`] this never
/// refines via BiWFA — it just collects raw per-hop anchor chains,
/// which is what the syng-native graph engine needs to derive pair
/// anchors across transitively-connected members.
///
/// Each returned chain carries its source-side anchor positions on the
/// *specific frontier region* it was discovered from (not the original
/// seed). Callers that need to compose chains from different hops must
/// track the frontier source themselves; for rDNA-like clusters where
/// single-hop misses most paralogs, just exposing the fuller set of
/// chains already recovers most pair-anchor opportunities.
pub fn query_transitive_with_anchors(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    max_depth: u16,
) -> io::Result<Vec<HomologousIntervalWithAnchors>> {
    let mut visited: FxHashSet<(String, u64, u64, char)> = FxHashSet::default();
    visited.insert((query_name.to_string(), query_start, query_end, '+'));
    let mut visited_nodes: FxHashSet<u32> = FxHashSet::default();
    let mut frontier: Vec<(String, u64, u64)> =
        vec![(query_name.to_string(), query_start, query_end)];
    let mut all_chains: Vec<HomologousIntervalWithAnchors> = Vec::new();
    for _ in 0..max_depth {
        if frontier.is_empty() {
            break;
        }
        let mut next_frontier: Vec<(String, u64, u64)> = Vec::new();
        for (q_name, q_start, q_end) in &frontier {
            let hits = match syng_index.query_region_with_anchors_ext_visited(
                q_name,
                *q_start,
                *q_end,
                padding,
                0,
                Some(&mut visited_nodes),
            ) {
                Ok(h) => h,
                Err(e) => {
                    warn!(
                        "syng transitive-with-anchors hop failed for {}:{}-{}: {}",
                        q_name, q_start, q_end, e
                    );
                    continue;
                }
            };
            for hit in hits {
                let key = (hit.genome.clone(), hit.start, hit.end, hit.strand);
                if visited.insert(key) {
                    next_frontier.push((hit.genome.clone(), hit.start, hit.end));
                    all_chains.push(hit);
                }
            }
        }
        frontier = next_frontier;
    }
    Ok(all_chains)
}

pub fn query_transitive_ext(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    max_depth: u16,
    query_extension: u64,
    extend_budget: u64,
    min_chain_anchors: usize,
    min_chain_fraction: f64,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    let merge_distance = extend_budget.min(i32::MAX as u64) as i32;
    query_transitive_ext_with_seed_filter(
        syng_index,
        query_name,
        query_start,
        query_end,
        padding,
        max_depth,
        query_extension,
        extend_budget,
        min_chain_anchors,
        min_chain_fraction,
        SyngSeedFilter::default(),
        merge_distance,
        sequence_index,
    )
}

/// [`query_transitive_ext`] with a high-frequency syncmer seed filter
/// applied independently at each BFS hop.
#[allow(clippy::too_many_arguments)]
pub fn query_transitive_ext_with_seed_filter(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    max_depth: u16,
    query_extension: u64,
    extend_budget: u64,
    min_chain_anchors: usize,
    min_chain_fraction: f64,
    seed_filter: SyngSeedFilter,
    merge_distance: i32,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    let mut visited: FxHashSet<(String, u64, u64, char)> = FxHashSet::default();
    visited.insert((query_name.to_string(), query_start, query_end, '+'));
    // Shared across all hops: any syncmer node consumed at an earlier hop
    // need not be re-entered. This collapses quadratic BFS cost on repeat-
    // dense regions (same node hit from many overlapping frontier ranges).
    let mut visited_nodes: FxHashSet<u32> = FxHashSet::default();

    let mut frontier: Vec<(String, u64, u64)> =
        vec![(query_name.to_string(), query_start, query_end)];
    let mut all_hits: Vec<HomologousInterval> = Vec::new();

    for depth in 0..max_depth {
        if frontier.is_empty() {
            break;
        }
        // Extension only widens the initial query. Subsequent hops already
        // land inside syncmer-connected neighborhoods, so extending them
        // compounds BFS fan-out without adding new reach.
        let hop_extension = if depth == 0 { query_extension } else { 0 };
        let mut next_frontier: Vec<(String, u64, u64)> = Vec::new();
        for (q_name, q_start, q_end) in &frontier {
            let hits = match one_hop_ext_visited_with_seed_filter(
                syng_index,
                q_name,
                *q_start,
                *q_end,
                padding,
                hop_extension,
                extend_budget,
                min_chain_anchors,
                min_chain_fraction,
                Some(&mut visited_nodes),
                seed_filter,
                merge_distance,
                sequence_index,
            ) {
                Ok(h) => h,
                Err(e) => {
                    warn!(
                        "syng-transitive hop {} failed for {}:{}-{}: {}",
                        depth, q_name, q_start, q_end, e
                    );
                    continue;
                }
            };
            let hits = merge_nearby_on_same_path(hits, merge_distance);
            for hit in hits {
                let key = (hit.genome.clone(), hit.start, hit.end, hit.strand);
                if visited.insert(key) {
                    next_frontier.push((hit.genome.clone(), hit.start, hit.end));
                    all_hits.push(hit);
                }
            }
        }
        frontier = next_frontier;
    }

    Ok(merge_nearby_on_same_path(all_hits, merge_distance))
}

/// Collapse nearby intervals that share the same
/// `(genome, strand)`. Different transitive paths can reach the same
/// homolog and produce intervals that differ by a few bp at the boundaries
/// (BiWFA refinement is not idempotent across hop depths), leaving tens of
/// near-duplicate rows per haplotype for a single query.
///
/// `merge_distance` is the user's `-d` query gap model. With `-d >= 0`, merge
/// same-path intervals separated by at most that many bp before the next
/// transitive hop. With `--no-merge` (`merge_distance < 0`), preserve the
/// historical strict-overlap collapse used only to remove duplicate hits.
fn merge_nearby_on_same_path(
    hits: Vec<HomologousInterval>,
    merge_distance: i32,
) -> Vec<HomologousInterval> {
    use std::collections::HashMap;
    let prof = std::env::var("SYNG_PROFILE").is_ok();
    let input_count = hits.len();
    let gap = merge_distance.max(0) as u64;
    // Merging loses CIGAR alignment context since merged coordinates
    // don't map back to either source CIGAR. Group by (genome, strand)
    // with (start, end, cigar); if exactly one entry ends up unmerged
    // in its group slot, keep its cigar; otherwise set None.
    let mut groups: HashMap<(String, char), Vec<(u64, u64, Option<Vec<u8>>)>> = HashMap::new();
    for h in hits {
        groups
            .entry((h.genome, h.strand))
            .or_default()
            .push((h.start, h.end, h.cigar));
    }
    let mut out: Vec<HomologousInterval> = Vec::new();
    let mut merged_groups = 0usize;
    let mut merge_examples: Vec<String> = Vec::new();
    for ((genome, strand), mut ivs) in groups {
        let group_input_count = ivs.len();
        ivs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut iter = ivs.into_iter();
        let Some(mut cur) = iter.next() else {
            continue;
        };
        let mut group_output_count = 0usize;
        let mut cur_merged = false;
        let mut group_example_parts: Vec<String> = if prof && group_input_count > 1 {
            vec![format!("{}({}):", genome, strand)]
        } else {
            Vec::new()
        };
        for next in iter {
            let should_merge = if merge_distance >= 0 {
                next.0 <= cur.1.saturating_add(gap)
            } else {
                next.0 < cur.1
            };
            if should_merge {
                if prof && merge_examples.len() < 20 {
                    group_example_parts.push(format!("{}-{}~{}-{}", cur.0, cur.1, next.0, next.1));
                }
                if next.1 > cur.1 {
                    cur.1 = next.1;
                }
                cur_merged = true;
            } else {
                out.push(HomologousInterval {
                    genome: genome.clone(),
                    start: cur.0,
                    end: cur.1,
                    strand,
                    cigar: if cur_merged { None } else { cur.2.take() },
                });
                group_output_count += 1;
                cur = next;
                cur_merged = false;
            }
        }
        out.push(HomologousInterval {
            genome,
            start: cur.0,
            end: cur.1,
            strand,
            cigar: if cur_merged { None } else { cur.2 },
        });
        group_output_count += 1;
        if group_output_count < group_input_count {
            merged_groups += 1;
            if prof && merge_examples.len() < 20 {
                merge_examples.push(group_example_parts.join(" "));
            }
        }
    }
    out.sort_by(|a, b| {
        a.genome
            .cmp(&b.genome)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.strand.cmp(&b.strand))
    });
    if prof {
        eprintln!(
            "PROF merge_same_path input={} output={} merged_groups={} gap={} examples={}",
            input_count,
            out.len(),
            merged_groups,
            merge_distance,
            merge_examples.join(" | ")
        );
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_anchor(q: u64, t: u64) -> Anchor {
        Anchor {
            query_pos: q,
            target_pos: t,
            node_id: 0,
        }
    }

    #[test]
    fn project_forward_strand_linear() {
        let a = mk_anchor(100, 200);
        assert_eq!(project_query_to_target(50, a, '+', 63, 10_000), 150);
        assert_eq!(project_query_to_target(400, a, '+', 63, 10_000), 500);
    }

    #[test]
    fn project_reverse_strand_linear() {
        let a = mk_anchor(100, 200);
        assert_eq!(project_query_to_target(50, a, '-', 63, 10_000), 313);
        assert_eq!(project_query_to_target(150, a, '-', 63, 10_000), 213);
    }

    #[test]
    fn project_clamps_to_target_bounds() {
        let a = mk_anchor(100, 200);
        assert_eq!(project_query_to_target(0, a, '+', 63, 10_000), 100);
        assert_eq!(project_query_to_target(100_000, a, '+', 63, 300), 300);
    }

    #[test]
    fn linear_projection_uses_leftmost_and_rightmost_anchors() {
        let hit = HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0,
            end: 10_000,
            strand: '+',
            anchors: vec![
                mk_anchor(550, 550),
                mk_anchor(1000, 1000),
                mk_anchor(1500, 1500),
                mk_anchor(2450, 2440),
            ],
        };
        let (s, e) = refine_by_linear_projection(&hit, 500, 2500, 63, 10_000);
        assert_eq!(s, 500);
        assert_eq!(e, 2490);
    }

    #[test]
    fn linear_projection_reverse_strand_reports_forward_bounds() {
        let hit = HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 1900,
            end: 3100,
            strand: '-',
            anchors: vec![mk_anchor(500, 3000), mk_anchor(1500, 2000)],
        };
        let (s, e) = refine_by_linear_projection(&hit, 400, 1600, 63, 10_000);
        assert_eq!(s, 1963);
        assert_eq!(e, 3100);
    }

    // syncmer_len matches default (w + k = 10 + 20 = 30 in yeast tests;
    // 63 in the pipeline). Tests use 63 for consistency with the
    // production default.
    const TEST_SYNCMER_LEN: u64 = 63;

    #[test]
    fn chain_anchors_empty_and_single() {
        assert!(chain_anchors(Vec::new(), 5_000, TEST_SYNCMER_LEN).is_empty());
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0,
            end: 100,
            strand: '+',
            anchors: vec![mk_anchor(10, 10)],
        }];
        let out = chain_anchors(hits, 5_000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 1);
    }

    #[test]
    fn chain_anchors_colinear_forms_one_chain() {
        // Four anchors on the primary diagonal within the queried range:
        // mutual-best links pair (i, i+1) for every i.
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0,
            end: 500,
            strand: '+',
            anchors: vec![
                mk_anchor(0, 0),
                mk_anchor(50, 50),
                mk_anchor(100, 100),
                mk_anchor(150, 150),
            ],
        }];
        let out = chain_anchors(hits, 500, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 4);
    }

    #[test]
    fn chain_anchors_range_constraint_rejects_out_of_range_pair() {
        // Two anchors on the same diagonal but 10 kb apart on query
        // axis. query_range_len = 1000 < 10_000 → no edge can form;
        // each anchor is its own singleton chain.
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0,
            end: 20_000,
            strand: '+',
            anchors: vec![mk_anchor(0, 0), mk_anchor(10_000, 10_000)],
        }];
        let out = chain_anchors(hits, 1_000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn chain_anchors_with_limits_can_span_large_bubble_flanks() {
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0,
            end: 90_000,
            strand: '+',
            anchors: vec![mk_anchor(0, 0), mk_anchor(80_000, 80_000)],
        }];
        let out = chain_anchors_with_limits(hits, TEST_SYNCMER_LEN, 90_000, 90_000, 90_000);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 2);
        assert_eq!(out[0].start, 0);
        assert_eq!(out[0].end, 80_000 + TEST_SYNCMER_LEN);
    }

    #[test]
    fn sweepga_scaffold_mass_filters_short_off_diagonal_hits() {
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g#0#chr1".into(),
            start: 0,
            end: 10_000,
            strand: '+',
            anchors: vec![
                // Long primary scaffold: five exact syncmer hits.
                mk_anchor(0, 1000),
                mk_anchor(100, 1100),
                mk_anchor(200, 1200),
                mk_anchor(300, 1300),
                mk_anchor(400, 1400),
                // Short off-diagonal coincidence: close enough to merge
                // locally, but below SweepGA scaffold mass.
                mk_anchor(0, 5000),
                mk_anchor(100, 5100),
            ],
        }];
        let out = chain_anchors_with_sweepga_scaffold_mass(
            hits,
            TEST_SYNCMER_LEN,
            500,
            3 * TEST_SYNCMER_LEN,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 5);
        assert_eq!(out[0].start, 1000);
        assert_eq!(out[0].end, 1400 + TEST_SYNCMER_LEN);
    }

    #[test]
    fn query_scaled_chain_limits_do_not_treat_alignment_budget_as_sv_budget() {
        let query_len = 83_361;
        let extend_budget = 1_000;
        let (_gap, _drift, span_cap) =
            query_scaled_chain_limits(query_len, extend_budget, TEST_SYNCMER_LEN, 10_000);
        assert!(
            span_cap > query_len + 100_000,
            "C4-scale queries need enough target span to retain local CNV/SV alleles"
        );
        assert!(
            span_cap > query_len + extend_budget * 2,
            "the boundary realignment budget must not cap target SV span"
        );
    }

    #[test]
    fn query_scaled_chain_limits_use_merge_distance_before_refinement() {
        let query_len = 1_000;
        let extend_budget = 100;
        let merge_distance = 10_000;
        let (gap, drift, span_cap) =
            query_scaled_chain_limits(query_len, extend_budget, TEST_SYNCMER_LEN, merge_distance);
        assert!(
            gap >= query_len + merge_distance as u64,
            "-d should widen the pre-refinement chain gap"
        );
        assert!(
            drift >= query_len + merge_distance as u64,
            "-d should widen the pre-refinement drift budget"
        );
        assert!(
            span_cap >= query_len + merge_distance as u64,
            "-d should contribute to the pre-refinement target span cap"
        );
    }

    #[test]
    fn effective_min_chain_anchors_scales_with_query_span() {
        assert_eq!(effective_min_chain_anchors(0, TEST_SYNCMER_LEN, 1_000), 1);
        assert_eq!(effective_min_chain_anchors(63, TEST_SYNCMER_LEN, 1_000), 1);
        assert_eq!(
            effective_min_chain_anchors(1_000, TEST_SYNCMER_LEN, 1_000),
            1
        );
        assert_eq!(
            effective_min_chain_anchors(10_000, TEST_SYNCMER_LEN, 1_000),
            2
        );
        assert_eq!(
            effective_min_chain_anchors(50_000, TEST_SYNCMER_LEN, 1_000),
            7
        );
        assert_eq!(
            effective_min_chain_anchors(100_000, TEST_SYNCMER_LEN, 1_000),
            14
        );
        assert_eq!(
            effective_min_chain_anchors(400_000, TEST_SYNCMER_LEN, 1_000),
            56
        );
        assert_eq!(effective_min_chain_anchors(50_000, TEST_SYNCMER_LEN, 3), 3);
        assert_eq!(effective_min_chain_anchors(50_000, TEST_SYNCMER_LEN, 0), 0);
    }

    #[test]
    fn merge_nearby_on_same_path_uses_user_gap_distance() {
        let hits = vec![
            HomologousInterval {
                genome: "g".into(),
                start: 0,
                end: 100,
                strand: '+',
                cigar: None,
            },
            HomologousInterval {
                genome: "g".into(),
                start: 150,
                end: 220,
                strand: '+',
                cigar: None,
            },
            HomologousInterval {
                genome: "g".into(),
                start: 400,
                end: 500,
                strand: '+',
                cigar: None,
            },
        ];
        let out = merge_nearby_on_same_path(hits, 50);
        assert_eq!(out.len(), 2);
        assert_eq!((out[0].start, out[0].end), (0, 220));
        assert_eq!((out[1].start, out[1].end), (400, 500));
    }

    #[test]
    fn merge_nearby_on_same_path_no_merge_keeps_gapped_intervals() {
        let hits = vec![
            HomologousInterval {
                genome: "g".into(),
                start: 0,
                end: 100,
                strand: '+',
                cigar: None,
            },
            HomologousInterval {
                genome: "g".into(),
                start: 100,
                end: 200,
                strand: '+',
                cigar: None,
            },
            HomologousInterval {
                genome: "g".into(),
                start: 190,
                end: 250,
                strand: '+',
                cigar: None,
            },
        ];
        let out = merge_nearby_on_same_path(hits, -1);
        assert_eq!(out.len(), 2);
        assert_eq!((out[0].start, out[0].end), (0, 100));
        assert_eq!((out[1].start, out[1].end), (100, 250));
    }

    #[test]
    fn chain_anchors_absorbs_within_drift_tolerance() {
        // Three anchors: (0, 0) [sig 0], (100, 100) [sig 0], (100, 500)
        // [sig 400]. With extend_budget=1000 the per-hop drift tolerance
        // equals 1000, so the 400bp diagonal jump from (100,100) to
        // (100,500) is admitted. Plane-sweep assembles one chain covering
        // all three anchors — intentionally more permissive than the old
        // mutual-best pipeline, which would split them.
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0,
            end: 1000,
            strand: '+',
            anchors: vec![mk_anchor(0, 0), mk_anchor(100, 100), mk_anchor(100, 500)],
        }];
        let out = chain_anchors(hits, 1000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 3);
    }

    #[test]
    fn chain_anchors_tandem_emits_parallel_chains() {
        // 3 query copies × 3 target copies at period 100. Mutual-best
        // pairs same-copy anchors on each of 5 diagonals (shift
        // -2, -1, 0, +1, +2). Primary diagonal (shift 0): three
        // anchors forming a chain of length 3. Off-diagonals have
        // fewer anchors and form shorter chains or singletons.
        let mut anchors = Vec::new();
        for qi in 0..3u64 {
            for tj in 0..3u64 {
                anchors.push(Anchor {
                    query_pos: qi * 100,
                    target_pos: 1000 + tj * 100,
                    node_id: 1,
                });
            }
        }
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 1000,
            end: 1300,
            strand: '+',
            anchors,
        }];
        let out = chain_anchors(hits, 500, TEST_SYNCMER_LEN);
        // Plane-sweep with extend_budget=500 admits within-drift moves;
        // all 9 anchors of the 3×3 grid share query_pos columns at
        // 0/100/200 with sig range within the drift window, so they
        // assemble into fewer, longer chains than mutual-best would.
        // Invariant: at least one chain exists and the longest chain
        // covers every query-column (length ≥ 3).
        assert!(!out.is_empty());
        let max_len = out.iter().map(|iv| iv.anchors.len()).max().unwrap();
        assert!(max_len >= 3);
    }

    #[test]
    fn chain_anchors_respects_strand_boundaries() {
        // +/- never co-chain.
        let hits = vec![
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 0,
                end: 100,
                strand: '+',
                anchors: vec![mk_anchor(0, 0)],
            },
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 50,
                end: 150,
                strand: '-',
                anchors: vec![mk_anchor(50, 150)],
            },
        ];
        let out = chain_anchors(hits, 10_000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 2);
    }

    fn mk_hit(
        genome: &str,
        start: u64,
        end: u64,
        strand: char,
        n_anchors: usize,
    ) -> HomologousIntervalWithAnchors {
        HomologousIntervalWithAnchors {
            genome: genome.into(),
            start,
            end,
            strand,
            anchors: (0..n_anchors)
                .map(|i| mk_anchor(start + i as u64, start + i as u64))
                .collect(),
        }
    }

    #[test]
    fn dedupe_keeps_majority_strand_on_overlap() {
        let hits = vec![
            mk_hit("X", 100, 500, '+', 50),
            mk_hit("X", 200, 400, '-', 3),
        ];
        let out = dedupe_strand_overlaps(hits);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, '+');
    }

    #[test]
    fn dedupe_preserves_distinct_real_inversion() {
        let hits = vec![
            mk_hit("X", 100, 500, '+', 50),
            mk_hit("X", 1000, 1200, '-', 40),
        ];
        let out = dedupe_strand_overlaps(hits);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn dedupe_is_per_path_not_cross_path() {
        let hits = vec![
            mk_hit("Y", 100, 500, '+', 50),
            mk_hit("Z", 200, 400, '-', 50),
        ];
        let out = dedupe_strand_overlaps(hits);
        assert_eq!(out.len(), 2);
    }
}
