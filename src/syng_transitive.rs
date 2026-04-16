//! Syng-seeded transitive query via boundary realignment.
//!
//! See `notes/SYNG_TRANSITIVE_DESIGN.md` for the algorithm description and
//! rationale. In short: syng returns homologous intervals at
//! syncmer-resolution plus a set of shared-syncmer anchors; we use BiWFA
//! (lib_wfa2 `MemoryMode::Ultralow`) on small windows between flanking
//! anchors to snap each fuzzy edge to base-pair precision, then iterate.

use std::cell::RefCell;
use std::io;

use lib_wfa2::affine_wavefront::{AffineWavefronts, Distance, MemoryMode};
use log::{debug, warn};
use rustc_hash::FxHashSet;

use crate::graph::reverse_complement;
use crate::sequence_index::{SequenceIndex as _, UnifiedSequenceIndex};
use crate::syng::{Anchor, HomologousInterval, HomologousIntervalWithAnchors, SyngIndex};

/// Maximum size (bp) of an anchor gap we'll BiWFA-realign across. If the
/// chosen flanking anchors are farther apart than this on the query side,
/// fall back to the syncmer-resolution coordinate for that edge — the
/// cost/benefit of local realignment inverts for very large windows.
const MAX_REALIGN_WINDOW_BP: u64 = 2048;

/// Extra bp collected outside the query region when enumerating shared
/// syncmers, so that each boundary has an anchor on either side. Without
/// this, the left edge (qStart) has no anchor with `query_pos <= qStart`
/// because all in-range anchors sit AT or ABOVE qStart. 512bp picks up a
/// comfortable number of flanking syncmers at default (k=8, w=55) syncmer
/// density — roughly 8 syncmers per side — without inflating the syng
/// hit graph meaningfully.
const ANCHOR_FLANK_BP: u64 = 512;

thread_local! {
    /// Dedicated BiWFA (Ultralow memory) edit-distance aligner for boundary
    /// realignment. Kept separate from the `with_aligner` pool in
    /// `src/impg.rs` — that pool uses `MemoryMode::High`, which would be
    /// wasteful (and unnecessary) for short boundary windows.
    static BIWFA_EDIT_ALIGNER: RefCell<Option<AffineWavefronts>> =
        const { RefCell::new(None) };
}

fn with_biwfa_edit_aligner<F, R>(f: F) -> R
where
    F: FnOnce(&mut AffineWavefronts) -> R,
{
    BIWFA_EDIT_ALIGNER.with(|cell| {
        let mut opt = cell.borrow_mut();
        if opt.is_none() {
            *opt = Some(Distance::Edit.create_aligner(None, Some(&MemoryMode::Ultralow)));
        }
        f(opt.as_mut().unwrap())
    })
}

/// Project a query-coordinate offset through a CIGAR onto the target.
///
/// Walks the CIGAR op-by-op, consuming `query_offset` bases from the query
/// axis. Returns the corresponding number of bases consumed on the target
/// axis. Returns `None` if the CIGAR does not cover `query_offset` bases.
///
/// CIGAR is the raw `&[u8]` returned by `AffineWavefronts::cigar()`:
/// per-base ops `M`, `X`, `=`, `I` (insertion in query relative to
/// target — consumes query only), `D` (deletion — consumes target only).
fn project_query_offset_via_cigar(cigar: &[u8], query_offset: u64) -> Option<u64> {
    let mut q: u64 = 0;
    let mut t: u64 = 0;
    for &op in cigar {
        if q >= query_offset {
            break;
        }
        match op {
            b'M' | b'X' | b'=' => {
                q += 1;
                t += 1;
            }
            b'I' => {
                q += 1;
            }
            b'D' => {
                t += 1;
            }
            _ => {
                // Unknown op — bail out conservatively.
                return None;
            }
        }
    }
    if q >= query_offset { Some(t) } else { None }
}

/// Pick flanking anchors for a query edge.
///
/// Returns `(left, right)` where `left.query_pos <= edge` is the rightmost
/// such anchor, and `right.query_pos > edge` is the leftmost such anchor.
/// Either side may be `None` if the edge lies outside the anchor set.
fn flanking_anchors(anchors: &[Anchor], edge: u64) -> (Option<Anchor>, Option<Anchor>) {
    let left = anchors
        .iter()
        .filter(|a| a.query_pos <= edge)
        .max_by_key(|a| a.query_pos)
        .copied();
    let right = anchors
        .iter()
        .filter(|a| a.query_pos > edge)
        .min_by_key(|a| a.query_pos)
        .copied();
    (left, right)
}

/// Resolve one query edge to a precise target forward-strand position using BiWFA.
///
/// `edge` is the query-side coordinate we want to project. `(left, right)`
/// are flanking anchors on the query axis: `left.query_pos <= edge <
/// right.query_pos`. Anchors store each syncmer's START position on both
/// query and target forward strands. Under forward homology the query
/// syncmer at q_pos aligns base-for-base to the target syncmer at t_pos;
/// under reverse (RC) homology, a[q_pos..q_pos+k] = RC(b[t_pos..t_pos+k]),
/// meaning a[q_pos] (the FIRST base of the query syncmer) corresponds to
/// b[t_pos + k - 1] (the LAST base of the target syncmer).
///
/// `syncmer_len` is `k` — used to account for the syncmer width under RC.
///
/// Returns the target forward-strand projection of `edge`, or `None` if no
/// precise resolution is possible.
fn resolve_edge_via_biwfa(
    query_name: &str,
    target_name: &str,
    edge: u64,
    left: Option<Anchor>,
    right: Option<Anchor>,
    sequence_index: &UnifiedSequenceIndex,
    target_strand: char,
    syncmer_len: u64,
) -> Option<u64> {
    let (l, r) = match (left, right) {
        (Some(l), Some(r)) => (l, r),
        _ => return None,
    };

    if r.query_pos <= l.query_pos {
        return None;
    }
    let gap_bp = r.query_pos - l.query_pos;
    if gap_bp > MAX_REALIGN_WINDOW_BP {
        debug!(
            "syng-transitive: anchor gap {}bp > {}bp cap on {}:{}, falling back to syncmer-res edge",
            gap_bp, MAX_REALIGN_WINDOW_BP, query_name, edge
        );
        return None;
    }

    // Compute the target-forward window and the anchor-point that query_pos
    // maps to on the target forward strand, per strand.
    //
    // For '+': the query-axis point q_pos corresponds base-for-base to
    //   target-axis point t_pos (syncmer starts aligned). The between-anchors
    //   window on target is [t_l, t_r) of length (t_r - t_l). The target
    //   anchor-point for the LEFT anchor is l.target_pos.
    //
    // For '-': under RC, a[q_pos..q_pos+k] = RC(b[t_pos..t_pos+k]). The
    //   first base of query syncmer a[q_pos] corresponds to b[t_pos + k - 1]
    //   (i.e. the LAST base of the target syncmer). In exclusive half-open
    //   convention, the target-forward position "just past" a[q_pos] is
    //   t_pos + k. The between-anchors target window (inclusive of the
    //   right-anchor's matching region) is [t_r + k, t_l + k), length
    //   (t_l - t_r). Under RC, target forward coord decreases as query
    //   forward coord increases, so the window's forward bases read in RC
    //   align to query forward bases. The anchor-point for l under RC is
    //   l.target_pos + k.
    let (t_fwd_start, t_fwd_end, t_anchor_point) = match target_strand {
        '+' => {
            if r.target_pos <= l.target_pos {
                return None;
            }
            (l.target_pos, r.target_pos, l.target_pos)
        }
        '-' => {
            if l.target_pos <= r.target_pos {
                return None;
            }
            (r.target_pos + syncmer_len, l.target_pos + syncmer_len, l.target_pos + syncmer_len)
        }
        _ => return None,
    };

    let q_slice = sequence_index
        .fetch_sequence(query_name, l.query_pos as i32, r.query_pos as i32)
        .ok()?;
    let t_slice_fwd = sequence_index
        .fetch_sequence(target_name, t_fwd_start as i32, t_fwd_end as i32)
        .ok()?;
    if q_slice.is_empty() || t_slice_fwd.is_empty() {
        return None;
    }

    let t_slice_oriented = if target_strand == '-' {
        reverse_complement(&t_slice_fwd)
    } else {
        t_slice_fwd
    };

    let q_offset_in_window = edge.saturating_sub(l.query_pos);
    if q_offset_in_window == 0 {
        return Some(t_anchor_point);
    }
    let q_slice_len = q_slice.len() as u64;
    if q_offset_in_window >= q_slice_len {
        // Edge past the right anchor in query — return the right anchor's
        // target-forward anchor point.
        let r_anchor_point = match target_strand {
            '+' => r.target_pos,
            '-' => r.target_pos + syncmer_len,
            _ => return None,
        };
        return Some(r_anchor_point);
    }

    let projected_in_oriented = with_biwfa_edit_aligner(|aligner| {
        aligner.align(&q_slice, &t_slice_oriented);
        let cigar = aligner.cigar();
        project_query_offset_via_cigar(cigar, q_offset_in_window)
    })?;

    // Convert the oriented-slice offset back to target forward coord.
    let refined = match target_strand {
        '+' => t_anchor_point.saturating_add(projected_in_oriented),
        '-' => t_anchor_point.saturating_sub(projected_in_oriented),
        _ => return None,
    };
    Some(refined)
}

/// Snap the two fuzzy edges of one homolog to base-pair precision on the
/// target forward strand.
///
/// `query_start` / `query_end` are the original query coordinates on
/// `query_name`. The homolog's `anchors` carry the shared-syncmer positions
/// that tie the query edges to the target; `homolog.strand` tells us whether
/// the homology is forward (`+`) or reverse-complement (`-`).
///
/// Returns `(start, end)` on the target forward strand — always start <= end
/// regardless of strand. For RC homology the projected target positions come
/// out reversed (qStart → larger target coord, qEnd → smaller) so the final
/// ordering is swapped here.
///
/// Edges with insufficient anchor support fall back to the syncmer-resolution
/// padded values from the homolog.
fn refine_boundaries(
    query_name: &str,
    query_start: u64,
    query_end: u64,
    homolog: &HomologousIntervalWithAnchors,
    sequence_index: &UnifiedSequenceIndex,
    syncmer_len: u64,
) -> (u64, u64) {
    if homolog.anchors.is_empty() {
        return (homolog.start, homolog.end);
    }

    // Each query edge projects to a position on the target forward strand.
    let (l_left, l_right) = flanking_anchors(&homolog.anchors, query_start);
    let t_for_qs = resolve_edge_via_biwfa(
        query_name,
        &homolog.genome,
        query_start,
        l_left,
        l_right,
        sequence_index,
        homolog.strand,
        syncmer_len,
    );

    let (r_left, r_right) = flanking_anchors(&homolog.anchors, query_end);
    let t_for_qe = resolve_edge_via_biwfa(
        query_name,
        &homolog.genome,
        query_end,
        r_left,
        r_right,
        sequence_index,
        homolog.strand,
        syncmer_len,
    );

    // For forward homology: qs projects to target_start, qe to target_end.
    // For RC homology: qs projects to target_END (larger coord), qe to
    // target_START (smaller coord). Use min/max to get the correct forward
    // strand bounds regardless.
    let (mut refined_start, mut refined_end) = match (t_for_qs, t_for_qe) {
        (Some(a), Some(b)) => (a.min(b), a.max(b)),
        // If one edge resolves and the other doesn't, use the resolved edge
        // on the correct side and fall back to the padded bound on the
        // other. For '+' strand: qs→start, qe→end. For '-' strand: qs→end,
        // qe→start.
        (Some(a), None) => match homolog.strand {
            '+' => (a, homolog.end),
            '-' => (homolog.start, a),
            _ => (homolog.start, homolog.end),
        },
        (None, Some(b)) => match homolog.strand {
            '+' => (homolog.start, b),
            '-' => (b, homolog.end),
            _ => (homolog.start, homolog.end),
        },
        (None, None) => (homolog.start, homolog.end),
    };

    // Clamp against the sequence length and keep start <= end.
    let target_len = sequence_index
        .get_sequence_length(&homolog.genome)
        .unwrap_or(usize::MAX) as u64;
    refined_start = refined_start.min(target_len);
    refined_end = refined_end.min(target_len);
    if refined_end < refined_start {
        std::mem::swap(&mut refined_start, &mut refined_end);
    }
    (refined_start, refined_end)
}

/// Distance-merge anchored intervals in place: any two intervals on the same
/// (genome, strand) whose target-axis gap is `<= merge_distance` are merged
/// and their anchor sets are unioned. Mirrors bedtools `merge -d`. A distance
/// of 0 is a no-op (raw output from `query_region_with_anchors` is already
/// overlap-merged).
fn distance_merge_anchored(
    mut hits: Vec<HomologousIntervalWithAnchors>,
    merge_distance: u64,
) -> Vec<HomologousIntervalWithAnchors> {
    if hits.len() <= 1 || merge_distance == 0 {
        return hits;
    }
    hits.sort_by(|a, b| {
        a.genome
            .cmp(&b.genome)
            .then(a.strand.cmp(&b.strand))
            .then(a.start.cmp(&b.start))
    });
    let mut merged: Vec<HomologousIntervalWithAnchors> = Vec::with_capacity(hits.len());
    for mut iv in hits {
        if let Some(last) = merged.last_mut() {
            if last.genome == iv.genome
                && last.strand == iv.strand
                && iv.start <= last.end.saturating_add(merge_distance)
            {
                last.end = last.end.max(iv.end);
                last.anchors.append(&mut iv.anchors);
                continue;
            }
        }
        merged.push(iv);
    }
    for iv in &mut merged {
        iv.anchors
            .sort_by(|a, b| a.query_pos.cmp(&b.query_pos).then(a.target_pos.cmp(&b.target_pos)));
        iv.anchors
            .dedup_by(|a, b| a.query_pos == b.query_pos && a.target_pos == b.target_pos);
    }
    merged
}

/// Run one hop of boundary-realignment transitive query.
///
/// Returns refined intervals (base-pair-precise endpoints) for all homologs
/// of the input query region. The syng query is internally widened by
/// [`ANCHOR_FLANK_BP`] on each side so that each edge of `[query_start,
/// query_end)` has shared-syncmer anchors on both sides for BiWFA
/// projection. Refinement uses the original (un-widened) edges.
///
/// `merge_distance` controls bedtools-style distance merging of padded
/// syncmer hits before realignment (honours `-d` from the CLI). 0 = only
/// merge overlapping intervals.
pub fn one_hop(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    merge_distance: u64,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    let widened_start = query_start.saturating_sub(ANCHOR_FLANK_BP);
    let query_total_len = sequence_index
        .get_sequence_length(query_name)
        .map(|l| l as u64)
        .unwrap_or(u64::MAX);
    let widened_end = query_end
        .saturating_add(ANCHOR_FLANK_BP)
        .min(query_total_len);
    let hits = syng_index.query_region_with_anchors(
        query_name,
        widened_start,
        widened_end,
        padding,
    )?;
    let hits = distance_merge_anchored(hits, merge_distance);
    let syncmer_len = (syng_index.params.w + syng_index.params.k) as u64;
    let mut out: Vec<HomologousInterval> = Vec::with_capacity(hits.len());
    for hit in &hits {
        let (refined_start, refined_end) = refine_boundaries(
            query_name,
            query_start,
            query_end,
            hit,
            sequence_index,
            syncmer_len,
        );
        // Skip hits whose refined interval collapses (e.g. the anchor set
        // only straddled one of the two edges and the fallback on the other
        // put it on the wrong side).
        if refined_end <= refined_start {
            continue;
        }
        out.push(HomologousInterval {
            genome: hit.genome.clone(),
            start: refined_start,
            end: refined_end,
            strand: hit.strand,
        });
    }
    Ok(out)
}

/// Run a multihop boundary-realignment transitive query.
///
/// At each hop, refined intervals from the previous hop's output are used
/// as new syng seeds. Visited `(genome, start, end)` tuples are deduped to
/// prevent cycles. Terminates after `max_depth` hops or when no new
/// intervals are discovered.
///
/// `merge_distance` is the bedtools-style `-d` merge distance applied to
/// padded syncmer hits before realignment, at every hop.
///
/// Because each hop snaps boundaries to base-pair precision via BiWFA,
/// slop does NOT compound across hops.
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
    let mut visited: FxHashSet<(String, u64, u64, char)> = FxHashSet::default();
    visited.insert((query_name.to_string(), query_start, query_end, '+'));

    let mut frontier: Vec<(String, u64, u64)> =
        vec![(query_name.to_string(), query_start, query_end)];
    let mut all_hits: Vec<HomologousInterval> = Vec::new();

    for depth in 0..max_depth {
        if frontier.is_empty() {
            break;
        }
        let mut next_frontier: Vec<(String, u64, u64)> = Vec::new();
        for (q_name, q_start, q_end) in &frontier {
            let hits = match one_hop(
                syng_index,
                q_name,
                *q_start,
                *q_end,
                padding,
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

    Ok(all_hits)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn project_cigar_simple_match() {
        // All matches — target offset == query offset.
        let cigar = b"MMMMMMMMMM";
        assert_eq!(project_query_offset_via_cigar(cigar, 0), Some(0));
        assert_eq!(project_query_offset_via_cigar(cigar, 5), Some(5));
        assert_eq!(project_query_offset_via_cigar(cigar, 10), Some(10));
    }

    #[test]
    fn project_cigar_with_insertion() {
        // I = query-only. Target offset lags query offset by the number of I's.
        // CIGAR: MM I MM  → query: 5, target: 4
        let cigar = b"MMIMM";
        assert_eq!(project_query_offset_via_cigar(cigar, 2), Some(2));
        assert_eq!(project_query_offset_via_cigar(cigar, 3), Some(2)); // consumed I
        assert_eq!(project_query_offset_via_cigar(cigar, 5), Some(4));
    }

    #[test]
    fn project_cigar_with_deletion() {
        // D = target-only. Target offset exceeds query offset by the number of D's.
        // CIGAR: MM D MM  → query: 4, target: 5
        // We project the MINIMAL target offset at a given query offset — so at
        // query=2 the answer is 2 (just before the D; the D is target-only and
        // is not yet "consumed" by a query base). At query=3 we must consume
        // the D before the next M, giving target=4. At query=4, target=5.
        let cigar = b"MMDMM";
        assert_eq!(project_query_offset_via_cigar(cigar, 2), Some(2));
        assert_eq!(project_query_offset_via_cigar(cigar, 3), Some(4));
        assert_eq!(project_query_offset_via_cigar(cigar, 4), Some(5));
    }

    #[test]
    fn project_cigar_short_circuits_past_end() {
        let cigar = b"MMM";
        assert_eq!(project_query_offset_via_cigar(cigar, 5), None);
    }

    #[test]
    fn flanking_anchors_normal() {
        let anchors = vec![
            Anchor { query_pos: 100, target_pos: 200, node_id: 1 },
            Anchor { query_pos: 300, target_pos: 450, node_id: 2 },
            Anchor { query_pos: 500, target_pos: 700, node_id: 3 },
        ];
        let (l, r) = flanking_anchors(&anchors, 400);
        assert_eq!(l.map(|a| a.query_pos), Some(300));
        assert_eq!(r.map(|a| a.query_pos), Some(500));
    }

    #[test]
    fn flanking_anchors_before_first() {
        let anchors = vec![
            Anchor { query_pos: 100, target_pos: 200, node_id: 1 },
        ];
        let (l, r) = flanking_anchors(&anchors, 50);
        assert!(l.is_none());
        assert_eq!(r.map(|a| a.query_pos), Some(100));
    }

    #[test]
    fn flanking_anchors_after_last() {
        let anchors = vec![
            Anchor { query_pos: 100, target_pos: 200, node_id: 1 },
        ];
        let (l, r) = flanking_anchors(&anchors, 500);
        assert_eq!(l.map(|a| a.query_pos), Some(100));
        assert!(r.is_none());
    }
}
