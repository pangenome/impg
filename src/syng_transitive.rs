//! Syng-seeded transitive query via anchor-based boundary refinement.
//!
//! For closely-related genomes (the intended use case) the boundaries of
//! homologous regions can be projected from the shared-syncmer anchors
//! alone — linear extrapolation from the innermost shared anchor to the
//! query edge is accurate to within the local indel burden, which for
//! such genomes is a handful of bp. No per-edge realignment or sequence
//! fetching is required.
//!
//! Pipeline per hop:
//!   1. `query_region_with_anchors` → raw padded syncmer intervals +
//!      per-hit shared-syncmer anchor positions.
//!   2. bedtools-style distance merge (honouring `-d`) → one merged
//!      interval per `(target, strand)` homology.
//!   3. For each merged interval, project the user's query edges onto
//!      the target forward strand by linear extrapolation from the
//!      innermost (leftmost-q and rightmost-q) anchors.
//!   4. Multihop: refined intervals feed the next hop's syng seed.
//!
//! See `notes/SYNG_TRANSITIVE_DESIGN.md`.

use std::cell::RefCell;
use std::io;

use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus, Distance, MemoryMode};
use log::warn;
use rustc_hash::FxHashSet;

use crate::graph::reverse_complement;
use crate::sequence_index::{SequenceIndex as _, UnifiedSequenceIndex};
use crate::syng::{Anchor, HomologousInterval, HomologousIntervalWithAnchors, SyngIndex};

/// Upper bound on the query-side outer span `q_anchor - qs` (or `qe -
/// q_anchor`) we'll realign. In practice the outer span is typically <
/// one syncmer gap (~30bp at default params); the cap is there to avoid
/// pathological cases (sparse anchors in divergent regions).
const EDGE_ALIGN_CAP_BP: u64 = 4096;

/// Extra bp fetched on the target side beyond the query-side window, to
/// accommodate net insertions / deletions in the inter-anchor region. Sized
/// to comfortably exceed typical indel burden over ≤4kb of query-outer span.
const EDGE_ALIGN_TARGET_BUFFER_BP: u64 = 256;

thread_local! {
    /// Edit-distance aligner for edge refinement. High-memory mode (small
    /// windows — 100–4kb — so memory isn't a concern, but EndsFree needs
    /// the full traceback that Ultralow's BiWFA doesn't always provide).
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
            *opt = Some(Distance::Edit.create_aligner(None, Some(&MemoryMode::High)));
        }
        f(opt.as_mut().unwrap())
    })
}

/// Align a small outer window (at a merged-homolog edge) against the
/// corresponding target window with the ANCHOR end fixed and the OUTER end
/// free on the target side (to absorb inter-anchor indels). Returns the
/// number of leading target bases the alignment skipped — equivalently,
/// the position in the oriented target window where query position 0
/// starts aligning.
///
/// `q_slice` and `t_slice_oriented` are both oriented so the anchor end is
/// at their right (last) end. `target_begin_buffer` is the number of
/// extra target bases fetched at the outer end, past where the query
/// window would land under zero-indel linear projection.
fn align_edge_to_anchor(
    q_slice: &[u8],
    t_slice_oriented: &[u8],
    target_begin_buffer: u64,
) -> Option<u64> {
    if q_slice.is_empty() || t_slice_oriented.is_empty() {
        return None;
    }
    // WFA2 requires text_begin_free <= |text|. Cap to the slice length;
    // if the target is shorter than the nominal buffer we still allow it
    // to skip up to its full length.
    let text_begin_free = target_begin_buffer.min(t_slice_oriented.len() as u64) as i32;
    with_edge_aligner(|aligner| {
        unsafe {
            lib_wfa2::bindings::wfa::wavefront_aligner_set_alignment_free_ends(
                aligner.aligner_mut(),
                0,                // pattern_begin_free: qs is fixed
                0,                // pattern_end_free: anchor end fixed
                text_begin_free,  // text_begin_free: outer end on target can skip
                0,                // text_end_free: anchor end fixed
            );
        }
        let status = aligner.align(q_slice, t_slice_oriented);
        if !matches!(status, AlignmentStatus::Completed) {
            return None;
        }
        let cigar = aligner.cigar();
        // WFA2 CIGAR convention (opposite of PAF/SAM):
        //   M / = / X  — both pattern and text advance
        //   I          — TEXT advances only  (skipped target base)
        //   D          — PATTERN advances only (skipped query base)
        // Leading 'I' ops = target bases skipped before the alignment
        // actually begins matching query.
        let skipped = cigar.iter().take_while(|&&op| op == b'I').count() as u64;
        Some(skipped)
    })
}

/// Co-linearity signature for an anchor on a given strand.
///
/// Two anchors on the same biological homology (modulo indels) will have
/// the same signature modulo the local indel burden; anchors on different
/// homologies (paralog, inversion elsewhere) will have very different
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

/// Median anchor signature for an interval (representative co-linearity
/// "level" of the interval's homology).
fn interval_signature(iv: &HomologousIntervalWithAnchors) -> Option<i64> {
    if iv.anchors.is_empty() {
        return None;
    }
    let mut sigs: Vec<i64> = iv.anchors.iter().map(|&a| anchor_signature(a, iv.strand)).collect();
    sigs.sort();
    Some(sigs[sigs.len() / 2])
}

/// Distance-merge anchored intervals in place: two intervals on the same
/// `(genome, strand)` merge when BOTH:
///   1. Their target-axis gap is `<= merge_distance` (bedtools `-d` style), AND
///   2. Their median anchor co-linearity signatures are within
///      `merge_distance / 10` of each other.
///
/// Criterion #2 prevents merging across paralogs / distinct homologies that
/// happen to sit within `-d` of each other but have very different
/// signatures (e.g. one at delta=1047 and one at delta=7275 on the same
/// contig — clearly two separate biological regions, not fragments of one).
///
/// `merge_distance = 0` is a no-op.
fn distance_merge_anchored(
    mut hits: Vec<HomologousIntervalWithAnchors>,
    merge_distance: u64,
) -> Vec<HomologousIntervalWithAnchors> {
    if hits.len() <= 1 || merge_distance == 0 {
        return hits;
    }
    let signature_tolerance = (merge_distance / 10).max(1) as i64;
    hits.sort_by(|a, b| {
        a.genome
            .cmp(&b.genome)
            .then(a.strand.cmp(&b.strand))
            .then(a.start.cmp(&b.start))
    });
    let mut merged: Vec<HomologousIntervalWithAnchors> = Vec::with_capacity(hits.len());
    let mut merged_sigs: Vec<Option<i64>> = Vec::with_capacity(hits.len());
    for mut iv in hits {
        let iv_sig = interval_signature(&iv);
        let can_merge = {
            if let (Some(last), Some(last_sig)) =
                (merged.last(), merged_sigs.last().and_then(|s| *s))
            {
                last.genome == iv.genome
                    && last.strand == iv.strand
                    && iv.start <= last.end.saturating_add(merge_distance)
                    && iv_sig
                        .map(|s| (s - last_sig).abs() <= signature_tolerance)
                        .unwrap_or(false)
            } else {
                false
            }
        };
        if can_merge {
            let last = merged.last_mut().unwrap();
            last.end = last.end.max(iv.end);
            last.anchors.append(&mut iv.anchors);
            // Update the cached signature for the (now larger) last interval.
            *merged_sigs.last_mut().unwrap() = interval_signature(last);
        } else {
            merged.push(iv);
            merged_sigs.push(iv_sig);
        }
    }
    for iv in &mut merged {
        iv.anchors
            .sort_by(|a, b| a.query_pos.cmp(&b.query_pos).then(a.target_pos.cmp(&b.target_pos)));
        iv.anchors
            .dedup_by(|a, b| a.query_pos == b.query_pos && a.target_pos == b.target_pos);
    }
    merged
}

/// Project a single query position onto the target forward strand via
/// linear extrapolation from a reference anchor.
///
/// For anchor `(q_anchor, t_anchor)` on strand `s`, the target-forward
/// coordinate corresponding to query base `q_pos`:
///
/// * `s = '+'`:  `t_anchor + (q_pos - q_anchor)` — increases with query.
/// * `s = '-'`:  `t_anchor + syncmer_len - (q_pos - q_anchor)` — decreases
///   with query. The `+ syncmer_len` accounts for the anchor's target
///   position pointing at the *start* of the syncmer, while under RC the
///   query base `q_anchor` corresponds to the syncmer's *last* target
///   base (exclusive: `t_anchor + syncmer_len`).
///
/// Underflow/overflow is clamped to `[0, target_len]`.
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

/// Refine one edge of a merged homolog by local realignment between the
/// innermost shared-syncmer anchor and the user's query boundary.
///
/// `is_left_edge = true` for qs projection (innermost = leftmost-q anchor);
/// `false` for qe (innermost = rightmost-q anchor).
///
/// Returns the target forward-strand coordinate where the query edge lands.
/// Falls back to linear projection when the outer span is zero, too large
/// for realignment, or a sequence fetch fails.
fn refine_edge_via_realignment(
    homolog: &HomologousIntervalWithAnchors,
    query_name: &str,
    query_edge: u64,
    inner_anchor: Anchor,
    is_left_edge: bool,
    syncmer_len: u64,
    target_len: u64,
    sequence_index: &UnifiedSequenceIndex,
) -> u64 {
    let linear = project_query_to_target(
        query_edge,
        inner_anchor,
        homolog.strand,
        syncmer_len,
        target_len,
    );
    let q_anchor = inner_anchor.query_pos;
    let t_anchor = inner_anchor.target_pos;

    // Outer span on query (from edge to the innermost anchor).
    let outer_span = if is_left_edge {
        q_anchor.saturating_sub(query_edge)
    } else {
        query_edge.saturating_sub(q_anchor + syncmer_len)
    };
    if outer_span == 0 || outer_span > EDGE_ALIGN_CAP_BP {
        return linear;
    }

    // Build query and target windows oriented so the anchor sits at the
    // RIGHT end of both (so the outer end is at position 0 and the
    // alignment can use the anchor as a fixed right-edge tether).
    //
    // For the LEFT edge: q_window = query[qs .. q_anchor + k]
    //   '+' strand: t_window = target[t_anchor + k - (q_window_len) - buffer .. t_anchor + k]
    //   '-' strand: t_window_fwd = target[t_anchor .. t_anchor + q_window_len + buffer],
    //               then RC it — the anchor ends up at the right, as intended.
    //
    // For the RIGHT edge the orientation is flipped: we want the anchor at
    // the LEFT end of both windows (so outer is on the right). The BiWFA call
    // with `text_end_free` would handle that symmetrically — but we can also
    // simply reverse both sequences and reuse the same left-edge routine.
    let buffer = EDGE_ALIGN_TARGET_BUFFER_BP;
    let q_window_len = outer_span + syncmer_len;

    let (q_slice, t_slice_fwd_range, rc_target) = if is_left_edge {
        let q_start = query_edge;
        let q_end = q_anchor + syncmer_len;
        let (t_start, t_end, rc) = match homolog.strand {
            '+' => {
                let t_end = t_anchor + syncmer_len;
                let t_start = t_end.saturating_sub(q_window_len + buffer);
                (t_start, t_end, false)
            }
            '-' => {
                let t_start = t_anchor;
                let t_end = t_anchor + syncmer_len + outer_span + buffer;
                (t_start, t_end.min(target_len), true)
            }
            _ => return linear,
        };
        ((q_start, q_end), (t_start, t_end), rc)
    } else {
        // RIGHT edge: flip both windows so the anchor ends up on the right.
        // Query window: query[q_anchor .. qe], reversed for alignment.
        let q_start = q_anchor;
        let q_end = query_edge;
        let (t_start, t_end, rc) = match homolog.strand {
            '+' => {
                // Anchor starts at t_anchor. Outward = rightward on target.
                let t_start = t_anchor;
                let t_end = t_anchor + q_window_len + buffer;
                (t_start, t_end.min(target_len), false)
            }
            '-' => {
                // Under RC: rightward on query = leftward on target forward.
                let t_end = t_anchor + syncmer_len;
                let t_start = t_end.saturating_sub(q_window_len + buffer);
                (t_start, t_end, true)
            }
            _ => return linear,
        };
        ((q_start, q_end), (t_start, t_end), rc)
    };

    let q_bytes = match sequence_index.fetch_sequence(
        query_name,
        q_slice.0 as i32,
        q_slice.1 as i32,
    ) {
        Ok(b) => b,
        Err(_) => return linear,
    };
    let t_bytes_fwd = match sequence_index.fetch_sequence(
        &homolog.genome,
        t_slice_fwd_range.0 as i32,
        t_slice_fwd_range.1 as i32,
    ) {
        Ok(b) => b,
        Err(_) => return linear,
    };

    // Orient both slices so the anchor is on the RIGHT end of both.
    let (q_oriented, t_oriented) = if is_left_edge {
        let t_o = if rc_target { reverse_complement(&t_bytes_fwd) } else { t_bytes_fwd.clone() };
        (q_bytes.clone(), t_o)
    } else {
        // Right edge: reverse both so the anchor lands at the right.
        let mut q_rev = q_bytes.clone();
        q_rev.reverse();
        let t_o = if rc_target {
            // Target was fetched with anchor at right in forward coords; RC for '-' strand
            // alignment, then reverse again so anchor lands right in oriented frame? Let's
            // think: for '-' strand RIGHT edge, anchor end on target-forward is at t_end
            // = t_anchor + syncmer_len. We fetched [t_end - (q_window_len + buffer), t_end).
            // RC puts the anchor at the LEFT of oriented. Reversing afterwards (making it
            // the RC's reverse = original forward) puts anchor back at right... confusing.
            // Simpler: for the right edge, build the slice already-anchor-right:
            //  '+': target[t_anchor .. t_anchor + q_window_len + buffer).
            //       reverse this so anchor (at left in forward) ends up on right.
            //  '-': target[t_end - (q_window_len + buffer) .. t_end).
            //       RC + reverse = complement only (identity on indexing).
            let mut t_rev = t_bytes_fwd.clone();
            t_rev.reverse();
            // Now the LEFT-most byte is the anchor's last base on target forward, i.e.
            // the base corresponding to q_anchor under RC. Complement each byte to get
            // the RC'd oriented slice (anchor on right now).
            for b in &mut t_rev {
                *b = match *b {
                    b'A' | b'a' => b'T',
                    b'T' | b't' => b'A',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    _ => *b,
                };
            }
            t_rev
        } else {
            // '+' strand right-edge: target is forward, anchor at LEFT of fetched slice.
            // Reverse the bytes to put the anchor at right (so we can reuse left-edge-style
            // EndsFree with free text_begin).
            let mut t_rev = t_bytes_fwd.clone();
            t_rev.reverse();
            t_rev
        };
        (q_rev, t_o)
    };

    // Align (anchor-end-fixed = right end; outer-end-free = left end on target).
    let skipped = match align_edge_to_anchor(&q_oriented, &t_oriented, buffer) {
        Some(s) => s,
        None => return linear,
    };

    // Translate `skipped` back to a target forward-strand position.
    // In oriented space, q_oriented[0] aligns to t_oriented[skipped].
    // t_oriented = positions in fetched target slice, possibly RC'd and/or reversed.
    // Mapping back:
    //   LEFT edge:
    //     '+': oriented = t_bytes_fwd. position `skipped` in oriented = position `skipped`
    //          in fwd slice. Target forward coord = t_slice_start + skipped.
    //     '-': oriented = RC(t_bytes_fwd). position `skipped` = RC of fwd[len-1-skipped].
    //          Target forward coord (exclusive) = t_slice_end - skipped.
    //   RIGHT edge (both sequences were reversed to put anchor at right):
    //     '+': oriented = reverse(t_bytes_fwd). position `skipped` = fwd[len-1-skipped].
    //          Target forward coord (exclusive) = t_slice_end - skipped.
    //     '-': oriented = reverse(RC(t_bytes_fwd)) = complement(t_bytes_fwd).
    //          Position `skipped` in oriented = complement(fwd[skipped]), which
    //          sits at target forward coord t_slice_start + skipped.
    let (t_start, t_end) = t_slice_fwd_range;
    let refined = if is_left_edge {
        match homolog.strand {
            '+' => t_start + skipped,
            '-' => t_end.saturating_sub(skipped),
            _ => linear,
        }
    } else {
        match homolog.strand {
            '+' => t_end.saturating_sub(skipped),
            '-' => t_start + skipped,
            _ => linear,
        }
    };
    refined.min(target_len)
}

/// Project the user's query edges onto the target forward strand for one
/// merged homolog. When `sequence_index` is provided, each edge is refined
/// by a local BiWFA realignment between the innermost anchor and the query
/// boundary — this captures indels that sit in the inter-anchor gap (which
/// linear projection from syncmer anchors cannot see). Otherwise (no
/// sequence access), falls back to pure linear extrapolation.
///
/// "Innermost" = leftmost-by-query_pos and rightmost-by-query_pos anchors
/// in the merged homolog.
fn refine_boundaries(
    homolog: &HomologousIntervalWithAnchors,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    syncmer_len: u64,
    target_len: u64,
    sequence_index: Option<&UnifiedSequenceIndex>,
) -> (u64, u64) {
    if homolog.anchors.is_empty() {
        return (homolog.start, homolog.end);
    }

    let mut anchors = homolog.anchors.clone();
    anchors.sort_by_key(|a| a.query_pos);
    let leftmost = anchors[0];
    let rightmost = *anchors.last().unwrap();

    let (t_for_qs, t_for_qe) = match sequence_index {
        Some(seq_idx) => (
            refine_edge_via_realignment(
                homolog, query_name, query_start, leftmost, true,
                syncmer_len, target_len, seq_idx,
            ),
            refine_edge_via_realignment(
                homolog, query_name, query_end, rightmost, false,
                syncmer_len, target_len, seq_idx,
            ),
        ),
        None => (
            project_query_to_target(query_start, leftmost, homolog.strand, syncmer_len, target_len),
            project_query_to_target(query_end, rightmost, homolog.strand, syncmer_len, target_len),
        ),
    };

    let (start, end) = if t_for_qs <= t_for_qe {
        (t_for_qs, t_for_qe)
    } else {
        (t_for_qe, t_for_qs)
    };
    (start.max(homolog.start), end.min(homolog.end))
}

/// Per-path strand dedupe: for each target path (= one haplotype's one
/// contig), if a `+` interval and a `-` interval overlap on the target
/// forward strand, keep only the one with more anchor support. The minority
/// strand is treated as noise from random-coincidence syncmer matches (short
/// k-mers whose canonical hash also hits at unrelated positions in the
/// opposite orientation).
///
/// Non-overlapping `+` / `-` intervals on the same path represent genuinely
/// distinct biological homologies (e.g. a collinear region plus an inversion
/// elsewhere on the same contig) and both are preserved.
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
            // Overlap on same path with opposite strands — majority wins.
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
/// `merge_distance` controls bedtools-style `-d` distance merging of padded
/// syncmer hits before edge projection. 0 = only merge overlapping intervals.
pub fn one_hop(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    merge_distance: u64,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    use crate::sequence_index::SequenceIndex as _;

    let hits = syng_index.query_region_with_anchors(
        query_name, query_start, query_end, padding,
    )?;
    let hits = distance_merge_anchored(hits, merge_distance);
    // Collapse strand-duplication noise: if a single path has a `+` and a
    // `-` interval overlapping on the target forward strand, the minority-
    // anchor strand is random syncmer coincidence — drop it.
    let hits = dedupe_strand_overlaps(hits);
    let syncmer_len = (syng_index.params.w + syng_index.params.k) as u64;

    let mut out: Vec<HomologousInterval> = Vec::with_capacity(hits.len());
    for hit in &hits {
        let target_len = sequence_index
            .get_sequence_length(&hit.genome)
            .map(|l| l as u64)
            .unwrap_or(u64::MAX);
        let (refined_start, refined_end) = refine_boundaries(
            hit,
            query_name,
            query_start,
            query_end,
            syncmer_len,
            target_len,
            Some(sequence_index),
        );
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

/// Run a multihop syng-seeded homology query.
///
/// At each hop, refined intervals from the previous hop's output are used as
/// new syng seeds. Visited `(genome, start, end, strand)` tuples are deduped
/// to prevent cycles. Terminates after `max_depth` hops or when no new
/// intervals are discovered.
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

    fn mk_anchor(q: u64, t: u64) -> Anchor {
        Anchor { query_pos: q, target_pos: t, node_id: 0 }
    }

    #[test]
    fn project_forward_strand_linear() {
        // Anchor at (100, 200) on '+': qs=50 → t=150, qe=400 → t=500.
        let a = mk_anchor(100, 200);
        assert_eq!(project_query_to_target(50, a, '+', 63, 10_000), 150);
        assert_eq!(project_query_to_target(400, a, '+', 63, 10_000), 500);
    }

    #[test]
    fn project_reverse_strand_linear() {
        // Anchor at (100, 200) on '-' with syncmer_len=63:
        //   qs=50  → t = 200 + 63 - (50 - 100)  = 263 + 50 = 313
        //   qe=150 → t = 200 + 63 - (150 - 100) = 263 - 50 = 213
        // I.e. increasing query → decreasing target (RC).
        let a = mk_anchor(100, 200);
        assert_eq!(project_query_to_target(50, a, '-', 63, 10_000), 313);
        assert_eq!(project_query_to_target(150, a, '-', 63, 10_000), 213);
    }

    #[test]
    fn project_clamps_to_target_bounds() {
        let a = mk_anchor(100, 200);
        // Extrapolation underflow → 0
        assert_eq!(project_query_to_target(0, a, '+', 63, 10_000), 100);
        // Extrapolation past target_len → clamp
        assert_eq!(project_query_to_target(100_000, a, '+', 63, 300), 300);
    }

    #[test]
    fn refine_uses_leftmost_and_rightmost_anchors() {
        // Homolog with several anchors; only leftmost & rightmost should drive the
        // refined edges.
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
        // No sequence_index → pure linear projection path.
        let (s, e) = refine_boundaries(&hit, "q", 500, 2500, 63, 10_000, None);
        // Left: leftmost anchor (550, 550), qs=500 → 550 + (500-550) = 500.
        // Right: rightmost anchor (2450, 2440), qe=2500 → 2440 + (2500-2450) = 2490.
        assert_eq!(s, 500);
        assert_eq!(e, 2490);
    }

    #[test]
    fn refine_reverse_strand_reports_forward_bounds() {
        // Under RC: as query increases, target decreases.
        // Anchors: (500, 3000), (1500, 2000) with syncmer_len=63.
        let hit = HomologousIntervalWithAnchors {
            genome: "g".into(),
            start: 1900,
            end: 3100,
            strand: '-',
            anchors: vec![mk_anchor(500, 3000), mk_anchor(1500, 2000)],
        };
        // qs=400: anchor (500, 3000) → t = 3000 + 63 - (400-500) = 3163.
        // qe=1600: anchor (1500, 2000) → t = 2000 + 63 - (1600-1500) = 1963.
        let (s, e) = refine_boundaries(&hit, "q", 400, 1600, 63, 10_000, None);
        // Forward-strand interval: min..max, clamped to [1900, 3100].
        assert_eq!(s, 1963);
        assert_eq!(e, 3100); // 3163 clamped to homolog.end=3100
    }

    #[test]
    fn distance_merge_zero_is_noop() {
        let hits = vec![
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 0, end: 100, strand: '+',
                anchors: vec![mk_anchor(0, 0)],
            },
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 200, end: 300, strand: '+',
                anchors: vec![mk_anchor(200, 200)],
            },
        ];
        let merged = distance_merge_anchored(hits, 0);
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn distance_merge_joins_within_d() {
        let hits = vec![
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 0, end: 100, strand: '+',
                anchors: vec![mk_anchor(0, 0)],
            },
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 200, end: 300, strand: '+',
                anchors: vec![mk_anchor(250, 250)],
            },
        ];
        let merged = distance_merge_anchored(hits, 150);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].end, 300);
        assert_eq!(merged[0].anchors.len(), 2);
    }

    #[test]
    fn distance_merge_refuses_when_signatures_diverge() {
        // Two intervals on the same (genome, strand), within -d on target
        // axis, but with very different co-linearity signatures — these are
        // paralogs, not fragments of one homology. Must NOT merge.
        // Anchors for first: delta=1047 (target - query). Second: delta=7275.
        let hits = vec![
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 409000, end: 410200, strand: '+',
                anchors: vec![
                    Anchor { query_pos: 408100, target_pos: 409147, node_id: 1 },
                    Anchor { query_pos: 408500, target_pos: 409547, node_id: 2 },
                    Anchor { query_pos: 409000, target_pos: 410047, node_id: 3 },
                ],
            },
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 416300, end: 417400, strand: '+',
                anchors: vec![
                    Anchor { query_pos: 409100, target_pos: 416375, node_id: 4 },
                    Anchor { query_pos: 409500, target_pos: 416775, node_id: 5 },
                    Anchor { query_pos: 409900, target_pos: 417175, node_id: 6 },
                ],
            },
        ];
        let merged = distance_merge_anchored(hits, 10_000);
        assert_eq!(
            merged.len(),
            2,
            "Signatures 1047 vs 7275 differ far beyond tolerance — must stay separate"
        );
    }

    #[test]
    fn distance_merge_joins_when_signatures_close() {
        // Same strand, within -d, and signatures within tolerance (within
        // a few bp of each other — real fragmentation of the same homology).
        let hits = vec![
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 100, end: 300, strand: '+',
                anchors: vec![
                    Anchor { query_pos: 100, target_pos: 200, node_id: 1 },
                    Anchor { query_pos: 200, target_pos: 300, node_id: 2 },
                ],
            },
            HomologousIntervalWithAnchors {
                genome: "g".into(),
                start: 500, end: 700, strand: '+',
                anchors: vec![
                    Anchor { query_pos: 400, target_pos: 500, node_id: 3 },
                    Anchor { query_pos: 500, target_pos: 600, node_id: 4 },
                ],
            },
        ];
        // Both intervals have signature delta=100.
        let merged = distance_merge_anchored(hits, 10_000);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].anchors.len(), 4);
    }

    #[test]
    fn distance_merge_respects_strand_boundaries() {
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
        let merged = distance_merge_anchored(hits, 10_000);
        // Different strands must not merge at this stage (cross-strand
        // dedupe happens separately).
        assert_eq!(merged.len(), 2);
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
        // On path X, a dense `+` interval at [100, 500) (50 anchors) overlaps
        // a sparse `-` interval at [200, 400) (3 anchors — random coincidence).
        // Dedupe should keep only the `+`.
        let hits = vec![
            mk_hit("X", 100, 500, '+', 50),
            mk_hit("X", 200, 400, '-', 3),
        ];
        let out = dedupe_strand_overlaps(hits);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, '+');
        assert_eq!((out[0].start, out[0].end), (100, 500));
    }

    #[test]
    fn dedupe_preserves_distinct_real_inversion() {
        // On path X: collinear `+` at [100, 500) AND real inversion `-` at
        // [1000, 1200) (non-overlapping target coords). Both kept.
        let hits = vec![
            mk_hit("X", 100, 500, '+', 50),
            mk_hit("X", 1000, 1200, '-', 40),
        ];
        let out = dedupe_strand_overlaps(hits);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn dedupe_is_per_path_not_cross_path() {
        // Two different paths: Y '+', Z '-'. They don't overlap logically
        // because they're different sequences. Both kept.
        let hits = vec![
            mk_hit("Y", 100, 500, '+', 50),
            mk_hit("Z", 200, 400, '-', 50),
        ];
        let out = dedupe_strand_overlaps(hits);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn dedupe_tie_prefers_first_wins() {
        // Equal anchor counts on overlap → keep the one we encountered first
        // (stable choice). Not a critical property but the algorithm must be
        // deterministic.
        let hits = vec![
            mk_hit("X", 100, 500, '+', 10),
            mk_hit("X", 200, 400, '-', 10),
        ];
        let out = dedupe_strand_overlaps(hits);
        assert_eq!(out.len(), 1);
        // `>=` in the comparator means the `+` (first) wins ties.
        assert_eq!(out[0].strand, '+');
    }
}
