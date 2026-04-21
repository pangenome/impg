//! Syng-seeded homology query with ends-only projection refinement.
//!
//! Pipeline per hop:
//!   1. `query_region_with_anchors` — syng returns padded syncmer-
//!      resolution intervals per (target, strand), with anchor positions.
//!      No K-cap on anchor collection; chain-level positional cap
//!      bounds memory downstream.
//!   2. `chain_anchors` — flatten per-visit anchors to per-(target,
//!      strand) anchor lists; greedy most-similar-diagonal assignment
//!      with both a signature tolerance (`merge_distance`) and a
//!      query-axis positional cap (`positional_cap_multiplier ×
//!      extend_budget`). Defeats the "same-diagonal-far-apart"
//!      pathology of pure diagonal clustering.
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
use crate::syng::{Anchor, HomologousInterval, HomologousIntervalWithAnchors, SyngIndex};

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

/// Mutual-best-buddy anchor chaining with range-bounded
/// diagonal-preferring scoring.
///
/// Caps are derived from `extend_budget` rather than the queried range
/// so chain shape (and therefore projection endpoints) are independent
/// of tile size. Liftover of any query position is a function of local
/// alignment evidence, not how the user chopped the query.
///
/// For a pair `(A, B)` with `B` downstream of `A` on query AND target
/// axes:
///
/// * **Hard constraints** (scale with `extend_budget`):
///   - `0 < dq ≤ 3 × extend_budget`       (max hop on query axis)
///   - `0 < dt ≤ 3 × extend_budget`       (max hop on target axis)
///   - `|dq − dt| ≤ 3 × extend_budget`    (max diagonal drift per hop)
/// * **Score** (lower is better):
///   `W_sig × |dq − dt| + max(dq, dt)` with `W_sig ≫ 1` — diagonal
///   match dominates; 2D proximity is the tiebreaker.
///
/// Chains form by mutual-best-buddy edges; union-find assembles
/// connected components. Chain-total span is capped at
/// `10 × extend_budget` during union so stacked small drifts don't
/// produce unbounded target spans.
pub(crate) fn chain_anchors(
    hits: Vec<HomologousIntervalWithAnchors>,
    extend_budget: u64,
    syncmer_len: u64,
) -> Vec<HomologousIntervalWithAnchors> {
    if hits.is_empty() {
        return hits;
    }

    use std::collections::HashMap;
    let mut per_path: HashMap<(String, char), Vec<Anchor>> = HashMap::new();
    for h in hits {
        per_path
            .entry((h.genome.clone(), h.strand))
            .or_default()
            .extend(h.anchors.into_iter());
    }

    // Iterative union-find with path compression.
    fn uf_find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    fn uf_union(parent: &mut [usize], rank: &mut [u32], x: usize, y: usize) {
        let rx = uf_find(parent, x);
        let ry = uf_find(parent, y);
        if rx == ry {
            return;
        }
        if rank[rx] < rank[ry] {
            parent[rx] = ry;
        } else if rank[rx] > rank[ry] {
            parent[ry] = rx;
        } else {
            parent[ry] = rx;
            rank[rx] += 1;
        }
    }

    const W_SIG: i64 = 1_000_000;
    // All hop caps scale with `extend_budget`, not the tile size.
    // `3 × budget` is a generous per-hop bound that admits normal
    // indel drift and small TE insertions while rejecting
    // wildly-discordant cross-paralog jumps.
    let hop_cap = extend_budget.saturating_mul(3).max(extend_budget);
    let max_dq = hop_cap;
    let max_dt = hop_cap;
    let sig_diff_abs_ceiling = hop_cap;

    let mut out: Vec<HomologousIntervalWithAnchors> = Vec::new();
    for ((genome, strand), mut anchors) in per_path {
        if anchors.is_empty() {
            continue;
        }
        anchors.sort_by(|a, b| {
            a.query_pos
                .cmp(&b.query_pos)
                .then(a.target_pos.cmp(&b.target_pos))
        });
        anchors.dedup_by(|a, b| a.query_pos == b.query_pos && a.target_pos == b.target_pos);
        let n = anchors.len();

        let mut best_fwd: Vec<Option<(usize, i64)>> = vec![None; n];
        let mut best_bwd: Vec<Option<(usize, i64)>> = vec![None; n];

        for i in 0..n {
            let a = anchors[i];
            for j in (i + 1)..n {
                let b = anchors[j];
                let dq = b.query_pos.saturating_sub(a.query_pos);
                if dq == 0 {
                    continue;
                }
                if dq > max_dq {
                    // Anchors are sorted by query_pos; further j can only
                    // have larger dq → abort inner loop.
                    break;
                }
                let dt_signed: i64 = match strand {
                    '+' => b.target_pos as i64 - a.target_pos as i64,
                    '-' => a.target_pos as i64 - b.target_pos as i64,
                    _ => 0,
                };
                if dt_signed <= 0 {
                    continue;
                }
                let dt = dt_signed as u64;
                if dt > max_dt {
                    continue;
                }
                let sig_diff_abs = (dq as i64 - dt as i64).unsigned_abs();
                if sig_diff_abs > sig_diff_abs_ceiling {
                    continue;
                }
                let prox = dq.max(dt) as i64;
                let s = (sig_diff_abs as i64)
                    .saturating_mul(W_SIG)
                    .saturating_add(prox);
                if best_fwd[i].map_or(true, |(_, bs)| s < bs) {
                    best_fwd[i] = Some((j, s));
                }
                if best_bwd[j].map_or(true, |(_, bs)| s < bs) {
                    best_bwd[j] = Some((i, s));
                }
            }
        }

        // Per-root target-axis span cap. Scaled by `extend_budget`, not
        // the queried range, so chain shape and projection endpoints
        // are independent of tile size. `10 × budget` admits TE/MEI
        // insertions plus normal homology extent; bigger chains would
        // indicate stacked discordant hops that shouldn't co-locate.
        let max_chain_span = extend_budget.saturating_mul(10).max(extend_budget);
        let mut parent: Vec<usize> = (0..n).collect();
        let mut rank: Vec<u32> = vec![0; n];
        let mut span_lo: Vec<u64> = anchors.iter().map(|a| a.target_pos).collect();
        let mut span_hi: Vec<u64> = anchors
            .iter()
            .map(|a| a.target_pos + syncmer_len)
            .collect();
        for i in 0..n {
            if let Some((j, _)) = best_fwd[i] {
                if let Some((k, _)) = best_bwd[j] {
                    if k == i {
                        let ri = uf_find(&mut parent, i);
                        let rj = uf_find(&mut parent, j);
                        if ri != rj {
                            let lo = span_lo[ri].min(span_lo[rj]);
                            let hi = span_hi[ri].max(span_hi[rj]);
                            if hi.saturating_sub(lo) <= max_chain_span {
                                uf_union(&mut parent, &mut rank, i, j);
                                let r = uf_find(&mut parent, i);
                                span_lo[r] = lo;
                                span_hi[r] = hi;
                            }
                        }
                    }
                }
            }
        }

        let mut by_root: HashMap<usize, Vec<Anchor>> = HashMap::new();
        for i in 0..n {
            let r = uf_find(&mut parent, i);
            by_root.entry(r).or_default().push(anchors[i]);
        }

        let mut raw_chains: Vec<(u64, u64, Vec<Anchor>)> = Vec::new();
        for (_, mut chain) in by_root {
            chain.sort_by(|a, b| a.query_pos.cmp(&b.query_pos));
            let mut tmin = u64::MAX;
            let mut tmax = 0u64;
            for a in &chain {
                tmin = tmin.min(a.target_pos);
                tmax = tmax.max(a.target_pos + syncmer_len);
            }
            raw_chains.push((tmin, tmax, chain));
        }

        // Post-merge pass: two chains on the same (genome, strand) that
        // are (a) adjacent on both axes with small bilateral gaps and
        // (b) colinear (gaps on query and target match within
        // tolerance) belong together. Best-buddy can fragment an
        // otherwise-clean chain when local anchor density creates an
        // edge where A's best forward ≠ B's best backward; this pass
        // heals those cuts.
        let merge_gap_cap: u64 = (extend_budget / 4).max(128);
        // Sort by target_start for + strand, target_end for - strand —
        // within a colinear chain, consecutive anchors have increasing
        // target_pos on +, decreasing on -. Either sort yields
        // physically adjacent chains for bilateral-gap checks.
        raw_chains.sort_by_key(|c| c.0);
        let mut merged_chains: Vec<(u64, u64, Vec<Anchor>)> = Vec::new();
        for (cs, ce, canchors) in raw_chains {
            if let Some(last) = merged_chains.last_mut() {
                // Query-axis bounds per chain
                let last_q_hi = last.2.last().map(|a| a.query_pos).unwrap_or(0);
                let cur_q_lo = canchors.first().map(|a| a.query_pos).unwrap_or(0);
                let last_q_lo = last.2.first().map(|a| a.query_pos).unwrap_or(0);
                let cur_q_hi = canchors.last().map(|a| a.query_pos).unwrap_or(0);
                let (q_gap, t_gap): (i64, i64) = match strand {
                    '+' => (
                        cur_q_lo as i64 - last_q_hi as i64,
                        cs as i64 - last.1 as i64,
                    ),
                    '-' => {
                        // Sort by target_start asc puts "later-on-query"
                        // chains FIRST on '-' (target decreasing with
                        // query). Use that: cur is upstream of last on
                        // query, so q_gap measured from cur_q_hi to
                        // last_q_lo.
                        (
                            last_q_lo as i64 - cur_q_hi as i64,
                            cs as i64 - last.1 as i64,
                        )
                    }
                    _ => (i64::MAX, i64::MAX),
                };
                let colinear =
                    (q_gap - t_gap).unsigned_abs() <= merge_gap_cap
                        && q_gap.unsigned_abs() <= merge_gap_cap
                        && t_gap.unsigned_abs() <= merge_gap_cap
                        && q_gap >= -(merge_gap_cap as i64)
                        && t_gap >= -(merge_gap_cap as i64);
                // Also enforce the chain-total span cap on the merge.
                let merged_span = ce.max(last.1).saturating_sub(cs.min(last.0));
                if colinear && merged_span <= max_chain_span {
                    last.0 = last.0.min(cs);
                    last.1 = last.1.max(ce);
                    last.2.extend(canchors);
                    last.2.sort_by(|a, b| a.query_pos.cmp(&b.query_pos));
                    continue;
                }
            }
            merged_chains.push((cs, ce, canchors));
        }

        for (tmin, tmax, chain) in merged_chains {
            out.push(HomologousIntervalWithAnchors {
                genome: genome.clone(),
                start: tmin,
                end: tmax,
                strand,
                anchors: chain,
            });
        }
    }
    out
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
        '-' => {
            t_anchor as i128 + syncmer_len as i128 - (q_pos as i128 - q_anchor as i128)
        }
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
        query_start, leftmost, homolog.strand, syncmer_len, target_len,
    );
    let t_for_qe = project_query_to_target(
        query_end, rightmost, homolog.strand, syncmer_len, target_len,
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
            if want_cigar { Some(cigar.to_vec()) } else { None },
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
    let small_slop = |query_gap: u64| -> u64 {
        (query_gap / 20).max(32).min(extend_budget)
    };
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
    let back_biwfa_q_lo = a_first.query_pos.saturating_sub(biwfa_back_gap).max(query_region_start);
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
            '+' => anchor_t.saturating_sub(linear_back_ext).max(prev_target_end),
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
            sq, ts, te, hit.strand, &hit.genome, back_outer, sequence_index, false,
        )?;
        let t_at_biwfa_lo = match hit.strand {
            '+' => ts + skip,
            '-' => te.saturating_sub(skip),
            _ => return None,
        };
        Some(linear_back_from(t_at_biwfa_lo))
    };
    let back_result = try_back(small_slop(biwfa_back_gap))
        .or_else(|| try_back(extend_budget));

    // ── Forward projection ──
    let fwd_query_gap = query_region_end
        .saturating_sub((a_last.query_pos + syncmer_len).min(query_region_end));
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
            '+' => (anchor_t + linear_fwd_ext).min(next_target_start).min(target_len),
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
            sq, ts, te, hit.strand, &hit.genome, fwd_outer, sequence_index, false,
        )?;
        let t_at_biwfa_hi = match hit.strand {
            '+' => te.saturating_sub(skip),
            '-' => ts + skip,
            _ => return None,
        };
        Some(linear_fwd_from(t_at_biwfa_hi))
    };
    let fwd_result = try_fwd(small_slop(biwfa_fwd_gap))
        .or_else(|| try_fwd(extend_budget));

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

/// Default minimum anchor count per chain for emission. Chains with
/// fewer anchors (singletons = 1) are filtered out — they have no
/// mutual-best colinear partner and are the weakest possible evidence.
pub const DEFAULT_MIN_CHAIN_ANCHORS: usize = 2;

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
        None,
        sequence_index,
    )
}

/// Core single-hop entry: chain anchors via mutual-best-buddy, then
/// ends-only projection per chain.
///
/// Chains with fewer than `min_chain_anchors` anchors are dropped at
/// emission (default 2: singletons filtered as weak evidence).
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
    visited_nodes: Option<&mut FxHashSet<u32>>,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    use std::time::Instant;
    let prof = std::env::var("SYNG_PROFILE").is_ok();
    let t0 = Instant::now();
    let hits = syng_index.query_region_with_anchors_ext_visited(
        query_name, query_start, query_end, padding, query_extension, visited_nodes,
    )?;
    if prof {
        eprintln!(
            "PROF emit={:?} hits={}",
            t0.elapsed(),
            hits.len()
        );
    }
    let syncmer_len = (syng_index.params.w + syng_index.params.k) as u64;
    let query_range_len = query_end.saturating_sub(query_start);
    let t1 = Instant::now();
    let total_anchors: usize = hits.iter().map(|h| h.anchors.len()).sum();
    let hits = chain_anchors(hits, extend_budget, syncmer_len);
    let pre_filter = hits.len();
    // Filter chains by anchor count. Default min=2 drops singletons
    // (one-anchor chains with no mutual-best colinear partner — mostly
    // noise in repeat-dense regions). Users who want to keep
    // singletons set --syng-min-chain-anchors 1.
    //
    // Auto-relax for short queries: a query range shorter than
    // 2 × syncmer_len can't physically host a chain of 2 distinct
    // syncmer anchors, so filtering would drop everything. Fall back
    // to min=1 when the range is too small for the requested minimum.
    let effective_min = if query_range_len < syncmer_len.saturating_mul(min_chain_anchors as u64) {
        1
    } else {
        min_chain_anchors
    };
    let hits: Vec<HomologousIntervalWithAnchors> = hits
        .into_iter()
        .filter(|h| h.anchors.len() >= effective_min)
        .collect();
    if prof {
        eprintln!(
            "PROF chain_anchors={:?} total_anchors={} chains={} after_filter={} (min={} eff={})",
            t1.elapsed(),
            total_anchors,
            pre_filter,
            hits.len(),
            min_chain_anchors,
            effective_min,
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
    let singleton_count = hits_with_neighbors.iter().filter(|(h, _, _)| h.anchors.len() < 2).count();
    let out: Vec<HomologousInterval> = hits_with_neighbors
        .par_iter()
        .with_min_len(128)
        .filter_map(|(hit, prev_end, next_start)| {
            let target_len = sequence_index
                .get_sequence_length(&hit.genome)
                .map(|l| l as u64)
                .unwrap_or(u64::MAX);

            let (refined_start, refined_end) = if hit.anchors.len() < 2 {
                refine_by_linear_projection(
                    hit, query_start, query_end, syncmer_len, target_len,
                )
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
                        hit, query_start, query_end, syncmer_len, target_len,
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
        groups.entry((h.genome.clone(), h.strand)).or_default().push(h);
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
#[allow(clippy::too_many_arguments)]
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
            let hits = match one_hop_ext_visited(
                syng_index,
                q_name,
                *q_start,
                *q_end,
                padding,
                hop_extension,
                extend_budget,
                min_chain_anchors,
                Some(&mut visited_nodes),
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

    Ok(merge_overlapping_on_same_path(all_hits))
}

/// Collapse strictly overlapping intervals that share the same
/// `(genome, strand)`. Different transitive paths can reach the same
/// homolog and produce intervals that differ by a few bp at the boundaries
/// (BiWFA refinement is not idempotent across hop depths), leaving tens of
/// near-duplicate rows per haplotype for a single query. Strictly touching
/// or gapped intervals are *not* merged — those represent genuinely adjacent
/// homologs (e.g. tandem arrays) that should stay separate.
fn merge_overlapping_on_same_path(hits: Vec<HomologousInterval>) -> Vec<HomologousInterval> {
    use std::collections::HashMap;
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
    for ((genome, strand), mut ivs) in groups {
        ivs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut iter = ivs.into_iter();
        let Some(mut cur) = iter.next() else {
            continue;
        };
        let mut cur_merged = false;
        for next in iter {
            if next.0 < cur.1 {
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
    }
    out.sort_by(|a, b| {
        a.genome
            .cmp(&b.genome)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.strand.cmp(&b.strand))
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_anchor(q: u64, t: u64) -> Anchor {
        Anchor { query_pos: q, target_pos: t, node_id: 0 }
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
    fn chain_anchors_mutual_best_rejects_non_colinear_triangle() {
        // Three anchors: (0, 0) [sig 0], (100, 100) [sig 0], (100, 500)
        // [sig 400]. (0,0)'s best forward buddy is (100, 100) — same
        // diagonal. (100, 100) and (100, 500) have dq=0 → invalid.
        // (0, 0)'s forward options: (100, 100) on sig 0, (100, 500) on
        // sig 400. Mutual-best picks (100, 100). The off-diagonal
        // anchor (100, 500) is a singleton.
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0,
            end: 1000,
            strand: '+',
            anchors: vec![
                mk_anchor(0, 0),
                mk_anchor(100, 100),
                mk_anchor(100, 500),
            ],
        }];
        let out = chain_anchors(hits, 1000, TEST_SYNCMER_LEN);
        // One chain of 2 (the colinear pair) + one singleton (the
        // off-diagonal anchor).
        assert_eq!(out.len(), 2);
        let with_two = out.iter().find(|iv| iv.anchors.len() == 2).unwrap();
        assert_eq!(with_two.anchors[0].query_pos, 0);
        assert_eq!(with_two.anchors[1].query_pos, 100);
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
        // Three diagonals (shift -1, 0, +1) have ≥2 anchors and pair up.
        // Others are singletons. Exact count depends on tiebreak but
        // there should be at least one chain with 3 anchors (primary).
        let max_len = out.iter().map(|iv| iv.anchors.len()).max().unwrap();
        assert_eq!(max_len, 3);
    }

    #[test]
    fn chain_anchors_respects_strand_boundaries() {
        // +/- never co-chain.
        let hits = vec![
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 0, end: 100, strand: '+',
                anchors: vec![mk_anchor(0, 0)],
            },
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 50, end: 150, strand: '-',
                anchors: vec![mk_anchor(50, 150)],
            },
        ];
        let out = chain_anchors(hits, 10_000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 2);
    }

    fn mk_hit(genome: &str, start: u64, end: u64, strand: char, n_anchors: usize) -> HomologousIntervalWithAnchors {
        HomologousIntervalWithAnchors {
            genome: genome.into(),
            start, end, strand,
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
