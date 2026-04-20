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
//!      it unless `--syng-emit-cigar` is set (approximate M-fill
//!      across the interior; interior gap-BiWFA is a follow-up).
//!   5. Multihop outer loop.

use std::cell::RefCell;
use std::io;

use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus, Distance, MemoryMode};
use log::warn;
use rustc_hash::FxHashSet;

use crate::graph::reverse_complement;
use crate::sequence_index::{SequenceIndex as _, UnifiedSequenceIndex};
use crate::syng::{Anchor, HomologousInterval, HomologousIntervalWithAnchors, SyngIndex};

/// Default multiplier applied to `extend_budget` to obtain the query-axis
/// positional cap used in `chain_anchors`. Two anchors cluster into the
/// same chain only if they are within `merge_distance` in signature AND
/// within `positional_cap = multiplier × extend_budget` bp on the query
/// axis.
pub const DEFAULT_POSITIONAL_CAP_MULTIPLIER: u64 = 4;

/// Minimum fraction of query bases that must align at 'M'/'=' ops for
/// the refined homolog to be accepted. Below this, the "homolog" is
/// treated as syncmer noise and the refinement falls back to syng's
/// padded bounds. Permissive by default (0.3) — only filters obvious
/// non-homology.
const MIN_ALIGNMENT_IDENTITY: f64 = 0.3;

thread_local! {
    /// Edit-distance aligner for per-homolog refinement.
    ///
    /// `MemoryMode::Low` = piggyback bidirectional WFA: O(s · log s)
    /// working memory vs `High`'s full O(s² + n·m) DP matrix. Traceback
    /// (which we need for leading / trailing-I counts in
    /// `refine_homolog_by_alignment`) is preserved in all modes — only the
    /// intermediate storage shrinks. On our hot path (thousands of per-hit
    /// refinements, small windows) `High` was the primary reason a single
    /// query could OOM at 77 GB RSS when a merged repeat cluster pushed
    /// any one target window into the Mb range.
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

/// Co-linearity signature for an anchor on a given strand.
///
/// Two anchors on the same biological homology (modulo indels) have the
/// same signature modulo the local indel burden; anchors on different
/// homologies (paralog, inversion elsewhere) have very different
/// signatures.
///
/// * `'+'` strand:  `signature = target_pos - query_pos`  (constant under
///   collinear homology — increases equally on both axes).
/// * `'-'` strand:  `signature = query_pos + target_pos`  (constant under
///   RC homology — target decreases as query increases, so the sum is
///   stable).
fn anchor_signature(a: Anchor, strand: char) -> i64 {
    match strand {
        '+' => a.target_pos as i64 - a.query_pos as i64,
        '-' => a.query_pos as i64 + a.target_pos as i64,
        _ => 0,
    }
}

/// Greedy per-anchor chaining with most-similar-diagonal assignment
/// and a query-axis positional cap.
///
/// Flattens all per-visit hits' anchors into per-(target, strand)
/// anchor lists, sorts each list by `query_pos`, and walks. For each
/// anchor `A`, candidate chains are those where:
///
/// * `|anchor_signature(A) − running_mean_sig(chain)| ≤ merge_distance`
///   (diagonal tolerance — absorbs small-indel drift)
/// * `A.query_pos − chain.last_anchor.query_pos ≤ positional_cap`
///   (locality — prevents two anchors on the "same virtual diagonal"
///   but far apart in genomic coordinates from collapsing into one
///   chain; that's the paralog-across-chromosome pathology of pure
///   signature clustering.)
///
/// If no candidate qualifies, a new chain starts with `A`. If multiple
/// qualify, the chain with the **most similar** running mean signature
/// wins — this keeps the correct repeat copy when several copies
/// overlap on the query axis.
///
/// Each emitted `HomologousIntervalWithAnchors` represents one chain
/// with `(start, end)` set to the anchor-supported target span
/// `[min(target_pos), max(target_pos) + syncmer_len]`.
pub(crate) fn chain_anchors(
    hits: Vec<HomologousIntervalWithAnchors>,
    merge_distance: u64,
    positional_cap: u64,
    syncmer_len: u64,
) -> Vec<HomologousIntervalWithAnchors> {
    if hits.is_empty() {
        return hits;
    }
    // `merge_distance = 0` is the CLI default and carries the classical
    // bedtools "overlap-only" semantics: don't *gate* chaining on
    // signature, just rely on the positional cap. Under that setting
    // every anchor on the same (target, strand) within `positional_cap`
    // joins one chain, independent of diagonal drift — which is what
    // users get when they don't pass `-d` explicitly. Positive values
    // apply signature tolerance in the natural way.
    let sig_tol: i64 = if merge_distance == 0 {
        i64::MAX
    } else {
        merge_distance as i64
    };

    use std::collections::HashMap;
    let mut per_path: HashMap<(String, char), Vec<Anchor>> = HashMap::new();
    for h in hits {
        per_path
            .entry((h.genome.clone(), h.strand))
            .or_default()
            .extend(h.anchors.into_iter());
    }

    struct Chain {
        anchors: Vec<Anchor>,
        sig_sum: i128,
        sig_count: i64,
        last_query_pos: u64,
    }

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

        let mut chains: Vec<Chain> = Vec::new();
        for a in anchors {
            let a_sig = anchor_signature(a, strand);
            let mut best: Option<(usize, i64)> = None;
            for (i, c) in chains.iter().enumerate() {
                // Chain is always non-empty by construction.
                let chain_mean: i64 = (c.sig_sum / c.sig_count as i128) as i64;
                let sig_diff = (a_sig - chain_mean).abs();
                if sig_diff > sig_tol {
                    continue;
                }
                let pos_diff = a.query_pos.saturating_sub(c.last_query_pos);
                if pos_diff > positional_cap {
                    continue;
                }
                best = match best {
                    Some((_, bd)) if bd <= sig_diff => best,
                    _ => Some((i, sig_diff)),
                };
            }
            match best {
                Some((i, _)) => {
                    let c = &mut chains[i];
                    c.anchors.push(a);
                    c.sig_sum += a_sig as i128;
                    c.sig_count += 1;
                    c.last_query_pos = a.query_pos;
                }
                None => chains.push(Chain {
                    anchors: vec![a],
                    sig_sum: a_sig as i128,
                    sig_count: 1,
                    last_query_pos: a.query_pos,
                }),
            }
        }

        for mut c in chains {
            c.anchors.sort_by(|a, b| {
                a.query_pos
                    .cmp(&b.query_pos)
                    .then(a.target_pos.cmp(&b.target_pos))
            });
            let mut tmin = u64::MAX;
            let mut tmax = 0u64;
            for a in &c.anchors {
                tmin = tmin.min(a.target_pos);
                tmax = tmax.max(a.target_pos + syncmer_len);
            }
            out.push(HomologousIntervalWithAnchors {
                genome: genome.clone(),
                start: tmin,
                end: tmax,
                strand,
                anchors: c.anchors,
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
    let t_bytes_fwd = sequence_index
        .fetch_sequence(genome, target_fwd_start as i32, target_fwd_end as i32)
        .ok()?;
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

    with_edge_aligner(|aligner| {
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
    })
}

/// Ends-only projection of one chain: two small ends-free BiWFAs,
/// one from the chain's first anchor backward to the query region
/// start, one from the last anchor forward to the query region end.
/// Each is bounded on the target axis by `extend_budget` and by the
/// neighbor chain's target edge on the same `(genome, strand)`.
/// Returns refined forward-strand `(ts, te)` and, when `emit_cigar`
/// is set, an approximate per-segment CIGAR formed by concatenating
/// the backward end-CIGAR + `M`×(interior query bp) + forward
/// end-CIGAR. The interior is trusted (every interior anchor is a
/// syncmer match); interior gap-BiWFA is a follow-up enhancement.
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
    emit_cigar: bool,
) -> Option<(u64, u64, Option<Vec<u8>>)> {
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

    // ── Backward projection: sub-query covers [query_region_start,
    //    a_first.query_pos + syncmer_len); its OUTER edge (toward
    //    query_region_start) is what we're trying to place on target.
    //    Forward-strand target span to fetch:
    //      '+': [a_first.target_pos − extend_budget  ..  a_first.target_pos + syncmer_len]
    //      '-': [a_first.target_pos  ..  a_first.target_pos + syncmer_len + extend_budget]
    //    (on '-', going backward on query is going forward on forward target.)
    let back_q_lo = query_region_start;
    let back_q_hi = (a_first.query_pos + syncmer_len).min(query_region_end);
    let (back_ts, back_te) = match hit.strand {
        '+' => (
            a_first
                .target_pos
                .saturating_sub(extend_budget)
                .max(prev_target_end),
            (a_first.target_pos + syncmer_len).min(target_len),
        ),
        '-' => (
            a_first.target_pos.min(target_len),
            (a_first.target_pos + syncmer_len + extend_budget)
                .min(next_target_start)
                .min(target_len),
        ),
        _ => return None,
    };
    let back_result = slice_sub_q(back_q_lo, back_q_hi).and_then(|sq| {
        // Backward projection: outer forward edge = Left for '+', Right for '-'.
        let outer = match hit.strand {
            '+' => Outer::Left,
            '-' => Outer::Right,
            _ => return None,
        };
        project_end_align(
            sq,
            back_ts,
            back_te,
            hit.strand,
            &hit.genome,
            outer,
            sequence_index,
            emit_cigar,
        )
    });

    // ── Forward projection: sub-query covers [a_last.query_pos,
    //    query_region_end); outer edge is toward query_region_end.
    let fwd_q_lo = a_last.query_pos.max(query_region_start);
    let fwd_q_hi = query_region_end;
    let (fwd_ts, fwd_te) = match hit.strand {
        '+' => (
            a_last.target_pos.min(target_len),
            (a_last.target_pos + syncmer_len + extend_budget)
                .min(next_target_start)
                .min(target_len),
        ),
        '-' => (
            a_last
                .target_pos
                .saturating_sub(extend_budget)
                .max(prev_target_end),
            (a_last.target_pos + syncmer_len).min(target_len),
        ),
        _ => return None,
    };
    let fwd_result = slice_sub_q(fwd_q_lo, fwd_q_hi).and_then(|sq| {
        let outer = match hit.strand {
            '+' => Outer::Right,
            '-' => Outer::Left,
            _ => return None,
        };
        project_end_align(
            sq,
            fwd_ts,
            fwd_te,
            hit.strand,
            &hit.genome,
            outer,
            sequence_index,
            emit_cigar,
        )
    });

    // Compose refined forward-strand bounds.
    //
    // Forward-strand outer edge mapping:
    //   '+' backward → leading bp on forward target from back_ts inward
    //       ⇒ refined start = back_ts + skip
    //   '+' forward → trailing bp on forward target from fwd_te inward
    //       ⇒ refined end = fwd_te - skip
    //   '-' backward (oriented outer = Right on forward view): trailing bp
    //       from back_te inward ⇒ refined end = back_te - skip
    //   '-' forward (oriented outer = Left on forward view): leading bp
    //       from fwd_ts inward ⇒ refined start = fwd_ts + skip
    let (refined_ts, refined_te, back_cigar, fwd_cigar) = match hit.strand {
        '+' => {
            let ts = match back_result.as_ref() {
                Some((skip, _)) => back_ts + *skip,
                None => a_first.target_pos,
            };
            let te = match fwd_result.as_ref() {
                Some((skip, _)) => fwd_te.saturating_sub(*skip),
                None => (a_last.target_pos + syncmer_len).min(target_len),
            };
            (
                ts,
                te,
                back_result.and_then(|(_, c)| c),
                fwd_result.and_then(|(_, c)| c),
            )
        }
        '-' => {
            let te = match back_result.as_ref() {
                Some((skip, _)) => back_te.saturating_sub(*skip),
                None => (a_first.target_pos + syncmer_len).min(target_len),
            };
            let ts = match fwd_result.as_ref() {
                Some((skip, _)) => fwd_ts + *skip,
                None => a_last.target_pos,
            };
            (
                ts,
                te,
                back_result.and_then(|(_, c)| c),
                fwd_result.and_then(|(_, c)| c),
            )
        }
        _ => return None,
    };

    if refined_te <= refined_ts {
        return None;
    }

    let combined_cigar = if emit_cigar {
        let mut cig: Vec<u8> = Vec::new();
        if let Some(c) = back_cigar {
            cig.extend_from_slice(&c);
        }
        // Interior fill: trust every interior anchor is a syncmer
        // match; emit 'M' across the anchor-supported query span.
        // This is an approximation — exact interior ops require
        // per-gap BiWFA (tracked as future work).
        let interior = a_last
            .query_pos
            .saturating_sub(a_first.query_pos + syncmer_len) as usize;
        cig.extend(std::iter::repeat(b'M').take(interior));
        if let Some(c) = fwd_cigar {
            cig.extend_from_slice(&c);
        }
        Some(cig)
    } else {
        None
    };

    Some((refined_ts, refined_te, combined_cigar))
}

/// Per-path strand dedupe: if a `+` and a `-` interval on the same
/// target path overlap on forward-strand coordinates, keep the majority-
/// anchor strand; drop the minority as random-syncmer-coincidence noise.
/// Non-overlapping `+`/`-` on the same path stay separate (real
/// inversion on a distinct region of the contig).
fn dedupe_strand_overlaps(
    hits: Vec<HomologousIntervalWithAnchors>,
) -> Vec<HomologousIntervalWithAnchors> {
    if hits.len() <= 1 {
        return hits;
    }
    let n = hits.len();
    let mut keep = vec![true; n];
    for i in 0..n {
        if !keep[i] {
            continue;
        }
        for j in (i + 1)..n {
            if !keep[j] {
                continue;
            }
            let a = &hits[i];
            let b = &hits[j];
            if a.genome != b.genome
                || a.strand == b.strand
                || a.start >= b.end
                || b.start >= a.end
            {
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
    hits.into_iter()
        .zip(keep)
        .filter_map(|(iv, k)| if k { Some(iv) } else { None })
        .collect()
}

/// Run one hop of syng-seeded homology query.
///
/// `merge_distance` = bedtools `-d`. 0 = overlap-only merge.
pub fn one_hop(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    merge_distance: u64,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    one_hop_ext(
        syng_index,
        query_name,
        query_start,
        query_end,
        padding,
        merge_distance,
        0,
        DEFAULT_EXTEND_BUDGET_BP,
        DEFAULT_POSITIONAL_CAP_MULTIPLIER,
        false,
        sequence_index,
    )
}

/// [`one_hop`] with tunable source-side widening and per-cluster
/// extension budget.
#[allow(clippy::too_many_arguments)]
pub fn one_hop_ext(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    merge_distance: u64,
    query_extension: u64,
    extend_budget: u64,
    positional_cap_multiplier: u64,
    emit_cigar: bool,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    one_hop_ext_visited(
        syng_index,
        query_name,
        query_start,
        query_end,
        padding,
        merge_distance,
        query_extension,
        extend_budget,
        positional_cap_multiplier,
        emit_cigar,
        None,
        sequence_index,
    )
}

/// Core single-hop entry: chain anchors, then ends-only projection
/// per chain.
///
/// Anchors are chained by `chain_anchors` under two caps:
/// `merge_distance` (diagonal tolerance) and `positional_cap =
/// positional_cap_multiplier × extend_budget` (query-axis locality).
/// Each chain produces a single `HomologousInterval` whose forward-
/// strand `(start, end)` is refined by two ends-free BiWFA calls —
/// one backward from the first anchor, one forward from the last —
/// each bounded on the target axis by `extend_budget` and by the
/// neighbor chain's target edge.
#[allow(clippy::too_many_arguments)]
pub fn one_hop_ext_visited(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    merge_distance: u64,
    query_extension: u64,
    extend_budget: u64,
    positional_cap_multiplier: u64,
    emit_cigar: bool,
    visited_nodes: Option<&mut FxHashSet<u32>>,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    let hits = syng_index.query_region_with_anchors_ext_visited(
        query_name, query_start, query_end, padding, query_extension, visited_nodes,
    )?;
    let syncmer_len = (syng_index.params.w + syng_index.params.k) as u64;
    let positional_cap = positional_cap_multiplier.saturating_mul(extend_budget);
    let hits = chain_anchors(hits, merge_distance, positional_cap, syncmer_len);
    let hits = dedupe_strand_overlaps(hits);

    // Fetch the full query bytes ONCE per hop; projection extracts
    // the sub-ranges its anchors + extension budget cover.
    let query_bytes_cache: Option<Vec<u8>> = sequence_index
        .fetch_sequence(query_name, query_start as i32, query_end as i32)
        .ok();

    // Attach per-chain target-axis neighbors (prev_end, next_start)
    // within each (genome, strand) bucket. These bound extension.
    let hits_with_neighbors = attach_target_neighbors(hits);

    let mut out: Vec<HomologousInterval> = Vec::with_capacity(hits_with_neighbors.len());
    for (hit, prev_end, next_start) in &hits_with_neighbors {
        let target_len = sequence_index
            .get_sequence_length(&hit.genome)
            .map(|l| l as u64)
            .unwrap_or(u64::MAX);

        let projected = query_bytes_cache.as_deref().and_then(|qb| {
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
                emit_cigar,
            )
        });
        let (refined_start, refined_end, cigar) = match projected {
            Some(t) => t,
            None => {
                let (s, e) = refine_by_linear_projection(
                    hit, query_start, query_end, syncmer_len, target_len,
                );
                (s, e, None)
            }
        };
        if refined_end <= refined_start {
            continue;
        }
        out.push(HomologousInterval {
            genome: hit.genome.clone(),
            start: refined_start,
            end: refined_end,
            strand: hit.strand,
            cigar,
        });
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
    merge_distance: u64,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    query_transitive_ext(
        syng_index,
        query_name,
        query_start,
        query_end,
        padding,
        max_depth,
        merge_distance,
        0,
        DEFAULT_EXTEND_BUDGET_BP,
        DEFAULT_POSITIONAL_CAP_MULTIPLIER,
        false,
        sequence_index,
    )
}

/// [`query_transitive`] with all tunables plumbed through every hop.
///
/// * `query_extension`: widens syncmer discovery on the source side at
///   each frontier step, so BFS can follow conserved blocks whose
///   endpoints fall just outside each frontier region.
/// * `extend_budget`: per-chain bp budget for ends-free target-side
///   extension in projection; see [`one_hop_ext_visited`].
/// * `positional_cap_multiplier`: scales `extend_budget` into the
///   query-axis locality cap for [`chain_anchors`].
/// * `emit_cigar`: attach a per-segment CIGAR to each emitted
///   `HomologousInterval` (approximate in the interior — exact for
///   the two ends-free projection runs).
#[allow(clippy::too_many_arguments)]
pub fn query_transitive_ext(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    max_depth: u16,
    merge_distance: u64,
    query_extension: u64,
    extend_budget: u64,
    positional_cap_multiplier: u64,
    emit_cigar: bool,
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
                merge_distance,
                hop_extension,
                extend_budget,
                positional_cap_multiplier,
                emit_cigar,
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
        assert!(chain_anchors(Vec::new(), 1000, 4000, TEST_SYNCMER_LEN).is_empty());
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0,
            end: 100,
            strand: '+',
            anchors: vec![mk_anchor(10, 10)],
        }];
        let out = chain_anchors(hits, 1000, 4000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 1);
    }

    #[test]
    fn chain_anchors_joins_matching_diagonal_within_positional_cap() {
        // Colinear anchors (sig=0) spaced 50 bp apart → one chain.
        let hits = vec![
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 0, end: 100, strand: '+',
                anchors: vec![mk_anchor(0, 0), mk_anchor(50, 50)],
            },
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 200, end: 300, strand: '+',
                anchors: vec![mk_anchor(200, 200), mk_anchor(250, 250)],
            },
        ];
        let out = chain_anchors(hits, 150, 4000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 4);
    }

    #[test]
    fn chain_anchors_splits_when_positional_cap_exceeded() {
        // Two anchor clusters on the SAME diagonal (sig=0) but 10 kb
        // apart on the query axis. positional_cap=4000 → they must
        // split into two chains even though the signatures match
        // perfectly. This is the test for the "same-diagonal, far-
        // apart" pathology.
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0, end: 20_000, strand: '+',
            anchors: vec![
                mk_anchor(0, 0),
                mk_anchor(100, 100),
                mk_anchor(10_000, 10_000),
                mk_anchor(10_100, 10_100),
            ],
        }];
        let out = chain_anchors(hits, 1000, 4000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn chain_anchors_splits_when_signature_gap_exceeds_merge_distance() {
        // Two anchor clusters in the same query locality but on
        // different diagonals (sig=1000 vs sig=7275).
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0, end: 20_000, strand: '+',
            anchors: vec![
                Anchor { query_pos: 100, target_pos: 1100, node_id: 1 },
                Anchor { query_pos: 200, target_pos: 1200, node_id: 2 },
                Anchor { query_pos: 300, target_pos: 7575, node_id: 3 },
                Anchor { query_pos: 400, target_pos: 7675, node_id: 4 },
            ],
        }];
        let out = chain_anchors(hits, 3000, 4000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 2);
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
        let out = chain_anchors(hits, 10_000, 10_000, TEST_SYNCMER_LEN);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn chain_anchors_most_similar_diagonal_tiebreak() {
        // Three anchors: two on sig=1000, one on sig=1100. All within
        // positional cap. With merge_distance=200, the third could
        // join either chain — but "most similar" picks the sig=1000
        // chain if it's closer. Here we set up so the tiebreak is
        // clear.
        let hits = vec![HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 0, end: 2000, strand: '+',
            anchors: vec![
                Anchor { query_pos: 100, target_pos: 1100, node_id: 1 },
                Anchor { query_pos: 200, target_pos: 1200, node_id: 2 },
                Anchor { query_pos: 300, target_pos: 1400, node_id: 3 },
            ],
        }];
        let out = chain_anchors(hits, 200, 4000, TEST_SYNCMER_LEN);
        // All three match diagonal within ±200 of sig=1000; positional
        // cap permits. Result: one chain.
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].anchors.len(), 3);
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
