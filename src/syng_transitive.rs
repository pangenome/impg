//! Syng-seeded homology query with per-homolog pairwise refinement.
//!
//! Pipeline per hop:
//!   1. `query_region_with_anchors` — syng returns padded syncmer-
//!      resolution intervals per (target, strand), with anchor positions.
//!   2. `distance_merge_anchored` — bedtools-style `-d` merge, with a
//!      co-linearity-signature guard that prevents merging paralogs.
//!   3. `dedupe_strand_overlaps` — per path, overlapping +/- intervals
//!      resolved by majority anchor count.
//!   4. `refine_homolog_by_alignment` — per merged homolog, BiWFA
//!      EndsFree between query bytes and a generously-padded target
//!      window; the CIGAR projects qs/qe onto target forward strand
//!      with bp precision. Linear-projection fallback when sequence
//!      fetch or alignment fails / scores too low.
//!   5. Multihop outer loop.
//!
//! See `notes/SYNG_OPTION3_PAIRWISE_REFINE.md` for the design rationale:
//! why anchor-only refinement was insufficient on RC homology, and how
//! per-homolog alignment restores PAF-quality boundary precision.

use std::cell::RefCell;
use std::io;

use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus, Distance, MemoryMode};
use log::warn;
use rustc_hash::FxHashSet;

use crate::graph::reverse_complement;
use crate::sequence_index::{SequenceIndex as _, UnifiedSequenceIndex};
use crate::syng::{Anchor, HomologousInterval, HomologousIntervalWithAnchors, SyngIndex};

/// Padding (bp) added on each side of syng's padded homolog bounds when
/// fetching target bytes for pairwise refinement. Must comfortably exceed
/// the worst anchor-to-boundary gap observed — ~1500bp on RC homologs in
/// yeast235 validation — so the real alignment endpoint lies within the
/// fetched window even when syng's bounds are loose.
const HOMOLOG_FETCH_PAD_BP: u64 = 2048;

/// Minimum fraction of query bases that must align at 'M'/'=' ops for
/// the refined homolog to be accepted. Below this, the "homolog" is
/// treated as syncmer noise and the refinement falls back to syng's
/// padded bounds. Permissive by default (0.3) — only filters obvious
/// non-homology.
const MIN_ALIGNMENT_IDENTITY: f64 = 0.3;

thread_local! {
    /// Edit-distance aligner for per-homolog refinement. High memory mode
    /// — windows are small (a few kb per side) and EndsFree semantics
    /// need the full traceback.
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

fn interval_signature(iv: &HomologousIntervalWithAnchors) -> Option<i64> {
    if iv.anchors.is_empty() {
        return None;
    }
    let mut sigs: Vec<i64> = iv.anchors.iter().map(|&a| anchor_signature(a, iv.strand)).collect();
    sigs.sort();
    Some(sigs[sigs.len() / 2])
}

/// Distance-merge anchored intervals: two intervals on the same
/// `(genome, strand)` merge when BOTH target-axis gap ≤ `merge_distance`
/// AND median anchor co-linearity signatures within `merge_distance / 10`.
///
/// The signature guard prevents merging paralogs that happen to sit
/// within `-d` of each other but belong to different biological regions.
pub(crate) fn distance_merge_anchored(
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

/// BiWFA-align query bytes against a padded target window for one merged
/// homolog. Returns refined forward-strand bounds, or None on
/// fetch/alignment failure or below-threshold identity.
///
/// The alignment is configured EndsFree on the target (both ends) and
/// End2End on the query — so the full query must align somewhere inside
/// the padded target window. Leading and trailing `I` ops (WFA2 CIGAR
/// convention: `I` = text-only advance) measure how much of the target
/// flanks fall outside the query's extent, giving us the refined
/// boundaries directly.
fn refine_homolog_by_alignment(
    query_bytes: &[u8],
    hit: &HomologousIntervalWithAnchors,
    sequence_index: &UnifiedSequenceIndex,
    target_len: u64,
) -> Option<(u64, u64)> {
    let fetch_start = hit.start.saturating_sub(HOMOLOG_FETCH_PAD_BP);
    let fetch_end = hit.end.saturating_add(HOMOLOG_FETCH_PAD_BP).min(target_len);
    if fetch_end <= fetch_start {
        return None;
    }
    let t_bytes_fwd = sequence_index
        .fetch_sequence(&hit.genome, fetch_start as i32, fetch_end as i32)
        .ok()?;
    if t_bytes_fwd.is_empty() || query_bytes.is_empty() {
        return None;
    }

    let t_bytes_oriented = if hit.strand == '-' {
        reverse_complement(&t_bytes_fwd)
    } else {
        t_bytes_fwd
    };
    let t_len = t_bytes_oriented.len() as u64;

    let result = with_edge_aligner(|aligner| {
        unsafe {
            lib_wfa2::bindings::wfa::wavefront_aligner_set_alignment_free_ends(
                aligner.aligner_mut(),
                0,                // pattern (query) must begin at 0
                0,                // pattern (query) must end at qe
                t_len as i32,     // text (target) can skip any prefix
                t_len as i32,     // text (target) can skip any suffix
            );
        }
        let status = aligner.align(query_bytes, &t_bytes_oriented);
        if !matches!(status, AlignmentStatus::Completed) {
            return None;
        }
        let cigar = aligner.cigar();
        let leading_i = cigar.iter().take_while(|&&op| op == b'I').count();
        let trailing_i = cigar.iter().rev().take_while(|&&op| op == b'I').count();
        // Identity gate: count M / = ops in the aligned portion only.
        let core_slice_end = cigar.len().saturating_sub(trailing_i);
        let matches = if core_slice_end > leading_i {
            cigar[leading_i..core_slice_end]
                .iter()
                .filter(|&&op| op == b'M' || op == b'=')
                .count()
        } else {
            0
        };
        let identity = matches as f64 / query_bytes.len() as f64;
        if identity < MIN_ALIGNMENT_IDENTITY {
            return None;
        }
        Some((leading_i as u64, trailing_i as u64))
    });

    let (leading_i, trailing_i) = result?;
    if leading_i + trailing_i >= t_len {
        // Degenerate: CIGAR was entirely flanking skips — no actual
        // alignment content.
        return None;
    }

    // Translate oriented-CIGAR offsets back to target forward-strand
    // coordinates.
    //
    // Oriented target window layout:
    //   '+' : identical to fetched forward slice;
    //         oriented[i] ↔ fwd[fetch_start + i].
    //   '-' : reverse-complement of fetched forward slice;
    //         oriented[i] (exclusive) ↔ fwd[fetch_end - i] (exclusive).
    //
    // Query aligns starting at oriented[leading_i] and ending at
    // oriented[t_len - trailing_i]. So:
    //   '+' : fwd_start = fetch_start + leading_i
    //         fwd_end   = fetch_end - trailing_i
    //   '-' : oriented[leading_i]        → fwd position (fetch_end - leading_i)  = forward END of homology
    //         oriented[t_len - trailing_i] → fwd position (fetch_start + trailing_i) = forward START
    let (refined_start, refined_end) = match hit.strand {
        '+' => (
            fetch_start + leading_i,
            fetch_end.saturating_sub(trailing_i),
        ),
        '-' => (
            fetch_start + trailing_i,
            fetch_end.saturating_sub(leading_i),
        ),
        _ => return None,
    };
    if refined_end <= refined_start {
        return None;
    }
    Some((refined_start, refined_end))
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
        sequence_index,
    )
}

/// [`one_hop`] with a `query_extension` bp source-side widening for
/// syncmer discovery. Lets syng reach homologs whose shared-syncmer
/// support ends just outside the user's declared query range; BiWFA
/// refinement still pairs against the original query bytes.
pub fn one_hop_ext(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    merge_distance: u64,
    query_extension: u64,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    let hits = syng_index.query_region_with_anchors_ext(
        query_name, query_start, query_end, padding, query_extension,
    )?;
    let hits = distance_merge_anchored(hits, merge_distance);
    let hits = dedupe_strand_overlaps(hits);
    let syncmer_len = (syng_index.params.w + syng_index.params.k) as u64;

    // Fetch query bytes ONCE per hop — the same slice is compared against
    // every homolog's target window. Avoids 270× redundant AGC reads.
    let query_bytes_cache: Option<Vec<u8>> = sequence_index
        .fetch_sequence(query_name, query_start as i32, query_end as i32)
        .ok();

    let mut out: Vec<HomologousInterval> = Vec::with_capacity(hits.len());
    for hit in &hits {
        let target_len = sequence_index
            .get_sequence_length(&hit.genome)
            .map(|l| l as u64)
            .unwrap_or(u64::MAX);

        // Preferred path: per-homolog pairwise alignment.
        let aligned = query_bytes_cache.as_deref().and_then(|qb| {
            refine_homolog_by_alignment(qb, hit, sequence_index, target_len)
        });
        let (refined_start, refined_end) = aligned.unwrap_or_else(|| {
            refine_by_linear_projection(hit, query_start, query_end, syncmer_len, target_len)
        });
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
        sequence_index,
    )
}

/// [`query_transitive`] with an optional `query_extension` plumbed through
/// every hop. The extension widens syncmer discovery on the source side at
/// each frontier step, so the BFS can follow conserved blocks whose
/// endpoints fall just outside each frontier region.
pub fn query_transitive_ext(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    max_depth: u16,
    merge_distance: u64,
    query_extension: u64,
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
        // Extension only widens the initial query. Subsequent hops already
        // land inside syncmer-connected neighborhoods, so extending them
        // compounds BFS fan-out without adding new reach.
        let hop_extension = if depth == 0 { query_extension } else { 0 };
        let mut next_frontier: Vec<(String, u64, u64)> = Vec::new();
        for (q_name, q_start, q_end) in &frontier {
            let hits = match one_hop_ext(
                syng_index,
                q_name,
                *q_start,
                *q_end,
                padding,
                merge_distance,
                hop_extension,
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
    let mut groups: HashMap<(String, char), Vec<(u64, u64)>> = HashMap::new();
    for h in hits {
        groups.entry((h.genome, h.strand)).or_default().push((h.start, h.end));
    }
    let mut out: Vec<HomologousInterval> = Vec::new();
    for ((genome, strand), mut ivs) in groups {
        ivs.sort_unstable();
        let mut iter = ivs.into_iter();
        let Some(mut cur) = iter.next() else {
            continue;
        };
        for next in iter {
            if next.0 < cur.1 {
                if next.1 > cur.1 {
                    cur.1 = next.1;
                }
            } else {
                out.push(HomologousInterval {
                    genome: genome.clone(),
                    start: cur.0,
                    end: cur.1,
                    strand,
                });
                cur = next;
            }
        }
        out.push(HomologousInterval {
            genome,
            start: cur.0,
            end: cur.1,
            strand,
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
    fn distance_merge_joins_within_d_with_matching_signatures() {
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
        let merged = distance_merge_anchored(hits, 150);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].end, 300);
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
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn distance_merge_refuses_when_signatures_diverge() {
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
