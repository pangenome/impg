//! Syng-seeded homology query with per-homolog pairwise refinement.
//!
//! Pipeline per hop:
//!   1. `query_region_with_anchors` — syng returns padded syncmer-
//!      resolution intervals per (target, strand), with anchor positions.
//!   2. `cluster_by_signature` — group hits into homology blocks by
//!      co-linearity signature (`target_pos − query_pos` on `+`, the
//!      sum on `-`). Each block is one colinear chain; paralog copies
//!      and insertions that shift signatures by > `merge_distance`
//!      fall into separate clusters.
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

fn interval_signature(iv: &HomologousIntervalWithAnchors) -> Option<i64> {
    if iv.anchors.is_empty() {
        return None;
    }
    let mut sigs: Vec<i64> = iv.anchors.iter().map(|&a| anchor_signature(a, iv.strand)).collect();
    sigs.sort();
    Some(sigs[sigs.len() / 2])
}

/// Signature-space clustering of anchored syncmer hits into homology
/// blocks.
///
/// Each hit carries a co-linearity signature (see [`anchor_signature`]):
///
/// * `'+'` strand: `signature = target_pos − query_pos`
/// * `'-'` strand: `signature = target_pos + query_pos`
///
/// Along a single colinear homology the signature is invariant modulo
/// the local indel burden. Paralog copies of the same sequence appear
/// at structurally different signatures — a tandem array at period P
/// has copies whose signatures differ by P; two paralogs elsewhere on
/// the same chromosome have signatures that differ by the paralog
/// offset. That's the right physical basis for block identity.
///
/// Algorithm: group by `(genome, strand)`, sort by signature, walk the
/// sorted list and start a new cluster whenever the signature gap
/// between consecutive hits exceeds `merge_distance`. Each cluster
/// becomes one `HomologousIntervalWithAnchors` spanning the union of
/// its members' target bounds and carrying the union of their anchors.
///
/// `merge_distance` is the user's `-d` / `--merge-distance` flowing
/// through unchanged. It's a single knob with a clean physical meaning:
/// the maximum signature excursion tolerated within one block. Typical
/// within-block jitter from small indels is tens of bp; paralog
/// separations are kb-scale on yeast; so `-d` up to ~1 kb reliably
/// merges within-homology chains without collapsing paralog copies.
pub(crate) fn cluster_by_signature(
    hits: Vec<HomologousIntervalWithAnchors>,
    merge_distance: u64,
) -> Vec<HomologousIntervalWithAnchors> {
    if hits.len() <= 1 || merge_distance == 0 {
        return hits;
    }
    let tol = merge_distance as i64;

    // Partition into (genome, strand) buckets; within each, sort by
    // signature and walk — clustering in one pass per path/strand.
    use std::collections::HashMap;
    let mut buckets: HashMap<(String, char), Vec<HomologousIntervalWithAnchors>> = HashMap::new();
    for h in hits {
        buckets
            .entry((h.genome.clone(), h.strand))
            .or_default()
            .push(h);
    }

    let mut out: Vec<HomologousIntervalWithAnchors> = Vec::new();
    for ((genome, strand), mut group) in buckets {
        // Hits with no anchors have no signature and can't be clustered;
        // pass them through unchanged.
        group.sort_by_key(|h| interval_signature(h).unwrap_or(i64::MIN));

        let mut cur: Option<HomologousIntervalWithAnchors> = None;
        let mut cur_sig: Option<i64> = None;
        for mut iv in group {
            let iv_sig = interval_signature(&iv);
            let start_new = match (cur_sig, iv_sig) {
                (Some(cs), Some(is)) => (is - cs).abs() > tol,
                // Never collapse across the "no-signature" boundary.
                _ => cur.is_some(),
            };
            if start_new {
                if let Some(mut done) = cur.take() {
                    finalize_cluster_anchors(&mut done);
                    out.push(done);
                }
            }
            match cur.as_mut() {
                Some(accum) => {
                    accum.end = accum.end.max(iv.end);
                    accum.start = accum.start.min(iv.start);
                    accum.anchors.append(&mut iv.anchors);
                    cur_sig = interval_signature(accum);
                }
                None => {
                    cur = Some(HomologousIntervalWithAnchors {
                        genome: genome.clone(),
                        start: iv.start,
                        end: iv.end,
                        strand,
                        anchors: std::mem::take(&mut iv.anchors),
                    });
                    cur_sig = iv_sig;
                }
            }
        }
        if let Some(mut done) = cur.take() {
            finalize_cluster_anchors(&mut done);
            out.push(done);
        }
    }
    out
}

/// Sort anchors within a finished cluster by `(query_pos, target_pos)`
/// and drop exact duplicates. Keeps projection fast and leaves the
/// signature-stable invariant intact.
fn finalize_cluster_anchors(iv: &mut HomologousIntervalWithAnchors) {
    iv.anchors
        .sort_by(|a, b| a.query_pos.cmp(&b.query_pos).then(a.target_pos.cmp(&b.target_pos)));
    iv.anchors
        .dedup_by(|a, b| a.query_pos == b.query_pos && a.target_pos == b.target_pos);
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

/// Padding (bp) added on each side of a cluster's anchor-covered query
/// span when extracting the sub-query bytes to feed BiWFA. Gives the
/// aligner a small margin of query sequence beyond the anchors to
/// refine the exact boundary.
const SUBQUERY_FLANK_PAD_BP: u64 = 200;

/// BiWFA-align one cluster's anchor-covered query span against a padded
/// target window. Returns refined forward-strand target bounds, or None
/// on fetch / alignment failure or below-threshold identity.
///
/// Key: we align only the *sub-query* covered by the cluster's anchors
/// (leftmost-anchor to rightmost-anchor + small flank padding), NOT the
/// user's full query region. A cluster of syncmer hits on one homology
/// block covers a specific contiguous stretch of query; forcing the
/// 21 kb full query to squeeze into a 4 kb target window was the reason
/// a single refinement took ~250 ms.
///
/// Alignment is configured EndsFree on the target (both ends) and
/// End2End on the sub-query — the sub-query must align, but the target
/// window can skip arbitrary prefix / suffix. Leading and trailing `I`
/// ops in the CIGAR (WFA2 convention: `I` = text-only advance) measure
/// how much of the target flanks fall outside the sub-query's extent,
/// giving us the refined boundaries directly.
fn refine_homolog_by_alignment(
    query_bytes_full: &[u8],
    query_region_start: u64,
    hit: &HomologousIntervalWithAnchors,
    sequence_index: &UnifiedSequenceIndex,
    target_len: u64,
) -> Option<(u64, u64)> {
    if hit.anchors.is_empty() {
        return None;
    }
    let syncmer_end_offset: u64 = 1; // anchors point at syncmer starts; include one base beyond to bound range

    // Sub-query span from cluster anchors (anchor query_pos values are
    // absolute query coordinates; subtract query_region_start to index
    // into query_bytes_full).
    let mut q_lo: u64 = u64::MAX;
    let mut q_hi: u64 = 0;
    for a in &hit.anchors {
        q_lo = q_lo.min(a.query_pos);
        q_hi = q_hi.max(a.query_pos.saturating_add(syncmer_end_offset));
    }
    let sub_lo = q_lo
        .saturating_sub(SUBQUERY_FLANK_PAD_BP)
        .max(query_region_start);
    let sub_hi = q_hi
        .saturating_add(SUBQUERY_FLANK_PAD_BP)
        .min(query_region_start + query_bytes_full.len() as u64);
    if sub_hi <= sub_lo {
        return None;
    }
    let sub_q_start_idx = (sub_lo - query_region_start) as usize;
    let sub_q_end_idx = (sub_hi - query_region_start) as usize;
    if sub_q_end_idx <= sub_q_start_idx || sub_q_end_idx > query_bytes_full.len() {
        return None;
    }
    let sub_query_bytes = &query_bytes_full[sub_q_start_idx..sub_q_end_idx];

    let fetch_start = hit.start.saturating_sub(HOMOLOG_FETCH_PAD_BP);
    let fetch_end = hit.end.saturating_add(HOMOLOG_FETCH_PAD_BP).min(target_len);
    if fetch_end <= fetch_start {
        return None;
    }
    let t_bytes_fwd = sequence_index
        .fetch_sequence(&hit.genome, fetch_start as i32, fetch_end as i32)
        .ok()?;
    if t_bytes_fwd.is_empty() || sub_query_bytes.is_empty() {
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
                0,                // pattern (sub-query) must begin at 0
                0,                // pattern (sub-query) must end at qe
                t_len as i32,     // text (target) can skip any prefix
                t_len as i32,     // text (target) can skip any suffix
            );
        }
        let status = aligner.align(sub_query_bytes, &t_bytes_oriented);
        if !matches!(status, AlignmentStatus::Completed) {
            return None;
        }
        let cigar = aligner.cigar();
        let leading_i = cigar.iter().take_while(|&&op| op == b'I').count();
        let trailing_i = cigar.iter().rev().take_while(|&&op| op == b'I').count();
        let core_slice_end = cigar.len().saturating_sub(trailing_i);
        let matches = if core_slice_end > leading_i {
            cigar[leading_i..core_slice_end]
                .iter()
                .filter(|&&op| op == b'M' || op == b'=')
                .count()
        } else {
            0
        };
        let identity = matches as f64 / sub_query_bytes.len() as f64;
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
    one_hop_ext_visited(
        syng_index,
        query_name,
        query_start,
        query_end,
        padding,
        merge_distance,
        query_extension,
        None,
        sequence_index,
    )
}

/// As [`one_hop_ext`], but threading a shared `visited_nodes` set through
/// syng so that syncmer nodes already consumed by an earlier BFS hop are
/// skipped at lookup time — avoiding re-running BiWFA refinement for
/// homologs that transitive-level dedup would drop after the fact.
pub fn one_hop_ext_visited(
    syng_index: &SyngIndex,
    query_name: &str,
    query_start: u64,
    query_end: u64,
    padding: u64,
    merge_distance: u64,
    query_extension: u64,
    visited_nodes: Option<&mut FxHashSet<u32>>,
    sequence_index: &UnifiedSequenceIndex,
) -> io::Result<Vec<HomologousInterval>> {
    let hits = syng_index.query_region_with_anchors_ext_visited(
        query_name, query_start, query_end, padding, query_extension, visited_nodes,
    )?;
    let hits = cluster_by_signature(hits, merge_distance);
    let hits = dedupe_strand_overlaps(hits);
    let syncmer_len = (syng_index.params.w + syng_index.params.k) as u64;

    // Fetch the full query bytes ONCE per hop. Per-cluster BiWFA
    // extracts just the sub-range its anchors cover, so this single
    // allocation is reused across every refinement without per-hit
    // AGC reads.
    let query_bytes_cache: Option<Vec<u8>> = sequence_index
        .fetch_sequence(query_name, query_start as i32, query_end as i32)
        .ok();

    let mut out: Vec<HomologousInterval> = Vec::with_capacity(hits.len());
    for hit in &hits {
        let target_len = sequence_index
            .get_sequence_length(&hit.genome)
            .map(|l| l as u64)
            .unwrap_or(u64::MAX);

        // Every cluster gets its own BiWFA pass — the sub-query bytes
        // covered by its anchors against its own small target window,
        // a fast alignment regardless of how wide the user's query
        // region is.
        let aligned = query_bytes_cache.as_deref().and_then(|qb| {
            refine_homolog_by_alignment(qb, query_start, hit, sequence_index, target_len)
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
    fn cluster_by_signature_zero_is_noop() {
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
        let merged = cluster_by_signature(hits, 0);
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn cluster_by_signature_joins_matching_signatures() {
        // Two hits, signatures both 0 (target == query). Any positive
        // merge_distance should cluster them together.
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
        let merged = cluster_by_signature(hits, 150);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 0);
        assert_eq!(merged[0].end, 300);
        assert_eq!(merged[0].anchors.len(), 4);
    }

    #[test]
    fn cluster_by_signature_respects_strand_boundaries() {
        // +/- on the same path never co-cluster — they come from
        // different strand semantics (distinct signatures).
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
        let merged = cluster_by_signature(hits, 10_000);
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn cluster_by_signature_splits_when_sig_gap_exceeds_merge_distance() {
        // Hit A: target-query == 1047. Hit B: target-query == 7275.
        // Sig gap = 6228. merge_distance = 3000 < 6228 → separate clusters.
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
        let merged = cluster_by_signature(hits, 3_000);
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn cluster_by_signature_merges_when_sig_gap_within_merge_distance() {
        // Same two hits as above but merge_distance = 10_000 > 6228 sig
        // gap → user said "merge within 10 kb", so they cluster.
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
        let merged = cluster_by_signature(hits, 10_000);
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].anchors.len(), 6);
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
