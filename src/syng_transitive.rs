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

use std::io;

use log::warn;
use rustc_hash::FxHashSet;

use crate::sequence_index::UnifiedSequenceIndex;
use crate::syng::{Anchor, HomologousInterval, HomologousIntervalWithAnchors, SyngIndex};

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

/// Project the user's query edges onto the target forward strand for one
/// merged homolog, using linear extrapolation from the innermost anchors.
///
/// "Innermost" = leftmost-by-query_pos and rightmost-by-query_pos anchors
/// in the merged homolog. These anchors sit somewhere inside the
/// homology region; the user's query edges `qs` / `qe` may sit anywhere
/// relative to them (inside the homology, at the boundary, or outside
/// in padding). Linear extrapolation gives a coordinate to within the
/// local indel burden — which for closely-related genomes is a handful
/// of bp, matching the precision of PAF-based queries for the same case.
fn refine_boundaries(
    homolog: &HomologousIntervalWithAnchors,
    query_start: u64,
    query_end: u64,
    syncmer_len: u64,
    target_len: u64,
) -> (u64, u64) {
    if homolog.anchors.is_empty() {
        // No anchor data → fall back to padded syncmer-resolution bounds.
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

    // For '+': t_for_qs <= t_for_qe. For '-': reversed. Sort for a
    // forward-strand half-open interval output.
    let (start, end) = if t_for_qs <= t_for_qe {
        (t_for_qs, t_for_qe)
    } else {
        (t_for_qe, t_for_qs)
    };

    // Clamp to the padded syncmer-resolution bounds. The padded bounds are
    // always a valid superset of the homology; projection outside them is
    // an extrapolation too far (e.g. qs far outside the homology region).
    (start.max(homolog.start), end.min(homolog.end))
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
    let syncmer_len = (syng_index.params.w + syng_index.params.k) as u64;

    let mut out: Vec<HomologousInterval> = Vec::with_capacity(hits.len());
    for hit in &hits {
        let target_len = sequence_index
            .get_sequence_length(&hit.genome)
            .map(|l| l as u64)
            .unwrap_or(u64::MAX);
        let (refined_start, refined_end) =
            refine_boundaries(hit, query_start, query_end, syncmer_len, target_len);
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
        let (s, e) = refine_boundaries(&hit, 500, 2500, 63, 10_000);
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
        let (s, e) = refine_boundaries(&hit, 400, 1600, 63, 10_000);
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
        // Different strands must not merge.
        assert_eq!(merged.len(), 2);
    }
}
