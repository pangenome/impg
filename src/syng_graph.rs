//! Syng-native graph engine primitives.
//!
//! (end-of-file) also contains `build_gfa_syng_native_from_sequences`,
//! which stitches the pairwise primitives into a full GFA via the existing
//! seqwish induction pipeline.
//!
//! This module provides the pieces needed to turn a set of homologous
//! sequences (e.g., members of a partition) into a variation graph without
//! running an external pairwise aligner (wfmash, fastGA). The key operation
//! is BiWFA between pairs, wrapped in a PAF emitter so the result can flow
//! through the existing seqwish induction pipeline.
//!
//! Staging plan (see `GfaEngine::SyngNative` in `lib.rs`):
//!   - v0 (this file): pairwise BiWFA → PAF line (tested in isolation).
//!   - v1: all-pairs driver + temp FASTA + seqwish integration.
//!   - v2: sweepga-knn-graph sparsification + anchor-seeded gap-only BiWFA
//!         + strand groom + indel left-align.
//!
//! BiWFA returns a per-base CIGAR byte string (e.g., `b"MMMIMMM"`); we
//! compact it to the `NumOp` form that PAF's `cg:Z:` tag expects.

use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus, Distance, MemoryMode};
use std::cell::RefCell;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

static DEBUG_GRAPH_BUILD_ID: AtomicUsize = AtomicUsize::new(0);

#[derive(Clone, Copy, Debug, Default)]
struct PafExactMatchRunStats {
    records: usize,
    parsed_records: usize,
    records_with_cigar: usize,
    name_mismatches: usize,
    length_mismatches: usize,
    max_exact_match_run: u64,
}

thread_local! {
    /// Gap-affine aligner reused across pair calls on the same thread.
    /// `MemoryMode::Low` (piggyback BiWFA) is faster on our size range
    /// (5 kb × low-divergence) than High mode per syng_transitive's
    /// benchmarking.
    static PAIR_ALIGNER: RefCell<Option<AffineWavefronts>> =
        const { RefCell::new(None) };
}

/// Gap-affine scoring that matches minimap2's defaults (asm5-ish):
/// match=0, mismatch=4, gap_open=6, gap_extend=2. Good starting point for
/// low-divergence within-species alignment. We align on global (end2end)
/// boundaries — the caller is responsible for picking sequences that
/// should align end-to-end (e.g., partition members of similar length).
const GAP_AFFINE_MISMATCH: i32 = 4;
const GAP_AFFINE_GAP_OPEN: i32 = 6;
const GAP_AFFINE_GAP_EXTEND: i32 = 2;

fn with_pair_aligner<F, R>(f: F) -> R
where
    F: FnOnce(&mut AffineWavefronts) -> R,
{
    PAIR_ALIGNER.with(|cell| {
        let mut opt = cell.borrow_mut();
        if opt.is_none() {
            *opt = Some(
                Distance::GapAffine {
                    mismatch: GAP_AFFINE_MISMATCH,
                    gap_opening: GAP_AFFINE_GAP_OPEN,
                    gap_extension: GAP_AFFINE_GAP_EXTEND,
                }
                .create_aligner(None, Some(&MemoryMode::Low)),
            );
        }
        f(opt.as_mut().unwrap())
    })
}

fn complement_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'U' => b'A',
        other => other,
    }
}

fn bases_equal_for_seqwish(query: u8, target: u8) -> bool {
    let query = query.to_ascii_uppercase();
    let target = target.to_ascii_uppercase();
    query == target && query != b'N'
}

fn paf_cigar<'a>(line_fields: &'a [&'a str]) -> Option<&'a str> {
    line_fields
        .iter()
        .skip(12)
        .find_map(|field| field.strip_prefix("cg:Z:"))
}

fn paf_line_max_exact_match_run(
    fields: &[&str],
    seq_by_name: &std::collections::HashMap<&str, &[u8]>,
) -> Option<(u64, bool)> {
    if fields.len() < 12 {
        return None;
    }
    let qseq = *seq_by_name.get(fields[0])?;
    let tseq = *seq_by_name.get(fields[5])?;
    let qlen = fields[1].parse::<usize>().ok()?;
    let qstart = fields[2].parse::<usize>().ok()?;
    let qend = fields[3].parse::<usize>().ok()?;
    let tlen = fields[6].parse::<usize>().ok()?;
    let tstart = fields[7].parse::<usize>().ok()?;
    let tend = fields[8].parse::<usize>().ok()?;
    if qlen != qseq.len()
        || tlen != tseq.len()
        || qstart > qend
        || tstart > tend
        || qend > qseq.len()
        || tend > tseq.len()
    {
        return None;
    }

    let Some(cigar) = paf_cigar(fields) else {
        // No base-level CIGAR: use the PAF block length as a conservative
        // upper bound. seqwish itself requires a CIGAR to split exact runs, but
        // this preserves old behaviour for PAF producers that emit single-run
        // alignments without cg:Z.
        return fields
            .get(10)
            .and_then(|value| value.parse::<u64>().ok())
            .map(|len| (len, false));
    };

    let query_is_rev = fields[4] == "-";
    let mut q_pos = if query_is_rev {
        qend.checked_sub(1)? as isize
    } else {
        qstart as isize
    };
    let q_step: isize = if query_is_rev { -1 } else { 1 };
    let mut t_pos = tstart;
    let mut digits = 0u64;
    let mut saw_op = false;
    let mut run = 0u64;
    let mut max_run = 0u64;

    for op in cigar.bytes() {
        if op.is_ascii_digit() {
            digits = digits.saturating_mul(10).saturating_add((op - b'0') as u64);
            continue;
        }
        if digits == 0 {
            return None;
        }
        saw_op = true;
        let len = digits as usize;
        digits = 0;
        match op {
            b'M' | b'=' | b'X' => {
                for _ in 0..len {
                    if q_pos < 0 || q_pos as usize >= qseq.len() || t_pos >= tseq.len() {
                        return None;
                    }
                    let mut query_base = qseq[q_pos as usize];
                    if query_is_rev {
                        query_base = complement_base(query_base);
                    }
                    let target_base = tseq[t_pos];
                    if bases_equal_for_seqwish(query_base, target_base) {
                        run += 1;
                        max_run = max_run.max(run);
                    } else {
                        run = 0;
                    }
                    q_pos += q_step;
                    t_pos += 1;
                }
            }
            b'I' => {
                run = 0;
                q_pos += q_step.saturating_mul(len as isize);
            }
            b'D' | b'N' => {
                run = 0;
                t_pos = t_pos.saturating_add(len);
            }
            b'S' => {
                run = 0;
                q_pos += q_step.saturating_mul(len as isize);
            }
            b'H' | b'P' => {
                run = 0;
            }
            _ => {
                run = 0;
            }
        }
    }
    if digits != 0 || !saw_op {
        return None;
    }
    Some((max_run, true))
}

fn filtered_paf_exact_match_run_stats(
    seqs: &[(String, Vec<u8>)],
    paf_path: &std::path::Path,
) -> std::io::Result<PafExactMatchRunStats> {
    use std::io::BufRead;

    let seq_by_name: std::collections::HashMap<&str, &[u8]> = seqs
        .iter()
        .map(|(name, seq)| (name.as_str(), seq.as_slice()))
        .collect();
    let file = std::fs::File::open(paf_path)?;
    let reader = std::io::BufReader::new(file);
    let mut stats = PafExactMatchRunStats::default();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        stats.records += 1;
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 12 {
            continue;
        }
        let names_match =
            seq_by_name.contains_key(fields[0]) && seq_by_name.contains_key(fields[5]);
        if !names_match {
            stats.name_mismatches += 1;
            continue;
        }
        let lengths_match = fields[1]
            .parse::<usize>()
            .ok()
            .zip(seq_by_name.get(fields[0]).map(|seq| seq.len()))
            .is_some_and(|(paf_len, seq_len)| paf_len == seq_len)
            && fields[6]
                .parse::<usize>()
                .ok()
                .zip(seq_by_name.get(fields[5]).map(|seq| seq.len()))
                .is_some_and(|(paf_len, seq_len)| paf_len == seq_len);
        if !lengths_match {
            stats.length_mismatches += 1;
            continue;
        }
        if let Some((max_run, has_cigar)) = paf_line_max_exact_match_run(&fields, &seq_by_name) {
            stats.parsed_records += 1;
            if has_cigar {
                stats.records_with_cigar += 1;
            }
            stats.max_exact_match_run = stats.max_exact_match_run.max(max_run);
        }
    }

    Ok(stats)
}

fn gfa_line_stats(gfa: &str) -> (usize, usize, usize, usize) {
    let mut segments = 0usize;
    let mut segment_bp = 0usize;
    let mut links = 0usize;
    let mut paths = 0usize;
    for line in gfa.lines() {
        if line.starts_with("S\t") {
            segments += 1;
            if let Some(seq) = line.split('\t').nth(2) {
                segment_bp += seq.len();
            }
        } else if line.starts_with("L\t") {
            links += 1;
        } else if line.starts_with("P\t") || line.starts_with("W\t") {
            paths += 1;
        }
    }
    (segments, segment_bp, links, paths)
}

#[derive(Default)]
pub(crate) struct PafStageMetrics {
    records: usize,
    directed_pairs: usize,
    undirected_pairs: usize,
    aligned_block_bp: u64,
    covered_sequences: usize,
    zero_coverage_sequences: usize,
    min_coverage_fraction: f64,
    mean_coverage_fraction: f64,
}

pub(crate) fn write_path_spelling_debug(
    dir: &std::path::Path,
    seqs: &[(String, Vec<u8>)],
) -> std::io::Result<()> {
    let mut out = String::from("sequence\tlength\n");
    for (name, seq) in seqs {
        out.push_str(name);
        out.push('\t');
        out.push_str(&seq.len().to_string());
        out.push('\n');
    }
    std::fs::write(dir.join("path_spellings.tsv"), out)
}

pub(crate) fn write_paf_stage_debug(
    dir: &std::path::Path,
    stage: &str,
    seqs: &[(String, Vec<u8>)],
    paf_path: &std::path::Path,
) -> std::io::Result<PafStageMetrics> {
    use std::io::BufRead;

    let seq_index = seqs
        .iter()
        .enumerate()
        .map(|(idx, (name, seq))| (name.as_str(), (idx, seq.len())))
        .collect::<std::collections::HashMap<_, _>>();
    let mut coverage_diff = seqs
        .iter()
        .map(|(_, seq)| vec![0i32; seq.len().saturating_add(1)])
        .collect::<Vec<_>>();
    let mut query_records = vec![0usize; seqs.len()];
    let mut target_records = vec![0usize; seqs.len()];
    let mut neighbors = vec![std::collections::HashSet::<usize>::new(); seqs.len()];
    let mut directed_pairs = std::collections::HashSet::<(usize, usize)>::new();
    let mut undirected_pairs = std::collections::HashSet::<(usize, usize)>::new();
    let mut metrics = PafStageMetrics::default();

    let file = std::fs::File::open(paf_path)?;
    let reader = std::io::BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 11 {
            continue;
        }
        metrics.records += 1;
        metrics.aligned_block_bp = metrics
            .aligned_block_bp
            .saturating_add(fields[10].parse::<u64>().unwrap_or(0));
        let Some(&(q_idx, q_len)) = seq_index.get(fields[0]) else {
            continue;
        };
        let Some(&(t_idx, t_len)) = seq_index.get(fields[5]) else {
            continue;
        };
        query_records[q_idx] += 1;
        target_records[t_idx] += 1;
        directed_pairs.insert((q_idx, t_idx));
        let pair = if q_idx <= t_idx {
            (q_idx, t_idx)
        } else {
            (t_idx, q_idx)
        };
        undirected_pairs.insert(pair);
        if q_idx != t_idx {
            neighbors[q_idx].insert(t_idx);
            neighbors[t_idx].insert(q_idx);
        }

        let mark = |diff: &mut [i32], len: usize, start: &str, end: &str| {
            let start = start.parse::<usize>().unwrap_or(0).min(len);
            let end = end.parse::<usize>().unwrap_or(0).min(len);
            if start < end {
                diff[start] += 1;
                diff[end] -= 1;
            }
        };
        mark(&mut coverage_diff[q_idx], q_len, fields[2], fields[3]);
        mark(&mut coverage_diff[t_idx], t_len, fields[7], fields[8]);
    }

    metrics.directed_pairs = directed_pairs.len();
    metrics.undirected_pairs = undirected_pairs.len();
    let mut coverage_rows =
        String::from("stage\tsequence\tlength\tcovered_bp\tcoverage_fraction\tquery_records\ttarget_records\tpair_degree\n");
    let mut total_fraction = 0.0;
    metrics.min_coverage_fraction = if seqs.is_empty() { 0.0 } else { 1.0 };
    for (idx, (name, seq)) in seqs.iter().enumerate() {
        let mut active = 0i32;
        let mut covered = 0usize;
        for delta in coverage_diff[idx].iter().take(seq.len()) {
            active += *delta;
            if active > 0 {
                covered += 1;
            }
        }
        let fraction = if seq.is_empty() {
            0.0
        } else {
            covered as f64 / seq.len() as f64
        };
        if covered > 0 {
            metrics.covered_sequences += 1;
        } else {
            metrics.zero_coverage_sequences += 1;
        }
        metrics.min_coverage_fraction = metrics.min_coverage_fraction.min(fraction);
        total_fraction += fraction;
        coverage_rows.push_str(&format!(
            "{stage}\t{name}\t{}\t{covered}\t{fraction:.6}\t{}\t{}\t{}\n",
            seq.len(),
            query_records[idx],
            target_records[idx],
            neighbors[idx].len()
        ));
    }
    if !seqs.is_empty() {
        metrics.mean_coverage_fraction = total_fraction / seqs.len() as f64;
    }

    std::fs::write(dir.join(format!("{stage}.coverage.tsv")), coverage_rows)?;
    let summary = format!(
        concat!(
            "stage\trecords\tdirected_pairs\tundirected_pairs\taligned_block_bp\t",
            "sequences\tcovered_sequences\tzero_coverage_sequences\t",
            "min_coverage_fraction\tmean_coverage_fraction\n",
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\n"
        ),
        stage,
        metrics.records,
        metrics.directed_pairs,
        metrics.undirected_pairs,
        metrics.aligned_block_bp,
        seqs.len(),
        metrics.covered_sequences,
        metrics.zero_coverage_sequences,
        metrics.min_coverage_fraction,
        metrics.mean_coverage_fraction
    );
    std::fs::write(dir.join(format!("{stage}.metrics.tsv")), summary)?;
    Ok(metrics)
}

/// Compact a byte CIGAR (one char per aligned position) into PAF-style
/// run-length form: `b"MMMIMMM"` → `"3M1I3M"`.
pub fn compact_cigar(cigar: &[u8]) -> String {
    if cigar.is_empty() {
        return String::new();
    }
    let mut out = String::new();
    let mut prev = cigar[0];
    let mut count = 1u32;
    for &op in &cigar[1..] {
        if op == prev {
            count += 1;
        } else {
            out.push_str(&count.to_string());
            out.push(prev as char);
            prev = op;
            count = 1;
        }
    }
    out.push_str(&count.to_string());
    out.push(prev as char);
    out
}

/// Walk a byte CIGAR once and return `(matches, mismatches, insertions, deletions)`.
/// Counts per-base, not per-op (so `"MMI"` → `matches=2, ins=1`).
pub fn cigar_stats(cigar: &[u8]) -> (u32, u32, u32, u32) {
    let mut m = 0u32;
    let mut x = 0u32;
    let mut i = 0u32;
    let mut d = 0u32;
    for &op in cigar {
        match op {
            b'M' | b'=' => m += 1,
            b'X' => x += 1,
            b'I' => i += 1,
            b'D' => d += 1,
            _ => {}
        }
    }
    (m, x, i, d)
}

/// Align two sequences end-to-end with BiWFA and return a single PAF
/// line with `cg:Z:` tag. Returns `None` if either sequence is empty or
/// the alignment does not complete.
///
/// PAF convention: sequence `a` is the query, `b` is the target. Both
/// are aligned forward-strand — callers must reverse-complement one
/// side before calling if the pair is `'-'` strand (e.g., per syng
/// chain strand in the partition-reference frame).
pub fn pairwise_biwfa_paf(
    a_name: &str,
    a_seq: &[u8],
    b_name: &str,
    b_seq: &[u8],
) -> Option<String> {
    if a_seq.is_empty() || b_seq.is_empty() {
        return None;
    }
    let cigar_bytes = with_pair_aligner(|aligner| {
        let status = aligner.align(a_seq, b_seq);
        if !matches!(status, AlignmentStatus::Completed) {
            return None;
        }
        Some(aligner.cigar().to_vec())
    })?;
    // WFA2 CIGAR convention (per WFA2-lib/alignment/cigar.c:387-392):
    //   I (insertion) consumes TEXT  — the second argument `b`
    //   D (deletion)  consumes PATTERN — the first argument `a`
    // PAF/SAM convention:
    //   I (insertion) consumes QUERY  — i.e. `a` when we call align(a=query, b=target)
    //   D (deletion)  consumes TARGET — i.e. `b`
    // So we must swap I↔D when converting WFA2 → PAF.
    let wfa_to_paf: Vec<u8> = cigar_bytes
        .iter()
        .map(|&op| match op {
            b'I' => b'D',
            b'D' => b'I',
            other => other,
        })
        .collect();

    // Left-align indels in the query frame. Canonicalizes microsat and
    // homopolymer indel placements → stops seqwish from seeing the same
    // event as two different variants. `a` is the query by our
    // convention, so the CIGAR is already in the query frame.
    let paf_cigar_bytes = crate::syng_graph_norm::left_align_indels(&wfa_to_paf, a_seq, b_seq);

    let (matches, mismatches, ins, del) = cigar_stats(&paf_cigar_bytes);
    let block_len = matches + mismatches + ins + del;
    if block_len == 0 {
        return None;
    }
    let a_len = a_seq.len();
    let b_len = b_seq.len();

    // Validate that the CIGAR consumes exactly a_len from query and
    // b_len from target. If it doesn't, the PAF line would break seqwish's
    // position math — drop it.
    let q_consumed = (matches + mismatches + ins) as usize;
    let t_consumed = (matches + mismatches + del) as usize;
    if q_consumed != a_len || t_consumed != b_len {
        log::debug!(
            "syng_graph: dropping pair {}/{} — cigar consumes q={} t={} but lens are q={} t={}",
            a_name,
            b_name,
            q_consumed,
            t_consumed,
            a_len,
            b_len
        );
        return None;
    }
    let cigar_str = compact_cigar(&paf_cigar_bytes);
    Some(format!(
        "{}\t{}\t{}\t{}\t+\t{}\t{}\t{}\t{}\t{}\t{}\t255\tcg:Z:{}\n",
        a_name, a_len, 0, a_len, b_name, b_len, 0, b_len, matches, block_len, cigar_str,
    ))
}

/// Build a PAF string from all pairs of sequences (i < j). Rayon-parallel
/// over pairs. Intended for small partitions or as a correctness baseline;
/// larger partitions should use `sparse_pairs_paf` which selects a
/// tree-kNN + stranger-joining + random-sample subset of pairs via
/// `sweepga::knn_graph::extract_tree_pairs`.
pub fn all_pairs_paf(seqs: &[(String, Vec<u8>)]) -> String {
    use rayon::prelude::*;
    let pairs: Vec<(usize, usize)> = (0..seqs.len())
        .flat_map(|i| ((i + 1)..seqs.len()).map(move |j| (i, j)))
        .collect();
    pairs
        .par_iter()
        .filter_map(|&(i, j)| pairwise_biwfa_paf(&seqs[i].0, &seqs[i].1, &seqs[j].0, &seqs[j].1))
        .collect::<Vec<_>>()
        .concat()
}

/// Build a PAF string from a sparsified pair set chosen via mash-distance
/// kNN + k-farthest ("stranger-joining") + deterministic-hashed random
/// subset. Delegates to `sweepga::knn_graph::extract_tree_pairs` for
/// selection and runs BiWFA in parallel on the selected pairs.
///
/// Sensible defaults for a partitioning use case: `k_nearest=3` (MST-ish
/// connectivity backbone), `k_farthest=1` (captures novel edges across
/// the cluster), `random_fraction=0.01` (completeness/robustness).
pub fn sparse_pairs_paf(
    seqs: &[(String, Vec<u8>)],
    k_nearest: usize,
    k_farthest: usize,
    random_fraction: f64,
) -> String {
    use rayon::prelude::*;
    if seqs.len() < 2 {
        return String::new();
    }
    let raw: Vec<Vec<u8>> = seqs.iter().map(|(_, s)| s.clone()).collect();
    let pairs = sweepga::knn_graph::extract_tree_pairs(
        &raw,
        k_nearest,
        k_farthest,
        random_fraction,
        &sweepga::knn_graph::MashParams::default(),
    );
    pairs
        .par_iter()
        .filter_map(|&(i, j)| pairwise_biwfa_paf(&seqs[i].0, &seqs[i].1, &seqs[j].0, &seqs[j].1))
        .collect::<Vec<_>>()
        .concat()
}

/// Anchor-seeded BiWFA: given two sequences and their shared syncmer
/// anchors (pairs of `(a_pos, b_pos)` on each), build a full-pair PAF
/// CIGAR by using each anchor as a match block and running BiWFA only
/// on the inter-anchor gaps.
///
/// `anchors` must be sorted by `a_pos` ascending, colinear on `b_pos`
/// (same strand, no crossings), and have `syncmer_len > 0`.
///
/// Returns None if:
///   - sequences are empty
///   - anchors are not colinear on `b_pos`
///   - any gap BiWFA fails to consume exactly the expected query/target bp
///
/// This is the v2.3 speedup over full-pair BiWFA: each BiWFA call
/// operates on ~30-300 bp gaps instead of the full 5 kb, cutting per-pair
/// cost by 10-50× on low-divergence input.
pub fn anchor_seeded_biwfa_paf(
    a_name: &str,
    a_seq: &[u8],
    b_name: &str,
    b_seq: &[u8],
    anchors: &[(usize, usize)],
    syncmer_len: usize,
) -> Option<String> {
    anchor_seeded_biwfa_paf_strand(a_name, a_seq, b_name, b_seq, anchors, syncmer_len, '+')
}

/// Strand-aware variant. Caller passes `paf_strand` ('+' or '-') which
/// becomes column 5 of the emitted PAF line. When `paf_strand == '-'`
/// the caller has already reverse-complemented `b_seq` and translated
/// anchor target positions to the rev-comp frame; the CIGAR is computed
/// on (a_seq, rev_comp_b) but the PAF reports forward-frame target
/// coords with strand '-' so seqwish reads it correctly.
pub fn anchor_seeded_biwfa_paf_strand(
    a_name: &str,
    a_seq: &[u8],
    b_name: &str,
    b_seq: &[u8],
    anchors: &[(usize, usize)],
    syncmer_len: usize,
    paf_strand: char,
) -> Option<String> {
    if a_seq.is_empty() || b_seq.is_empty() || syncmer_len == 0 {
        return None;
    }
    let a_len = a_seq.len();
    let b_len = b_seq.len();

    // Verify colinearity and bounds.
    for w in anchors.windows(2) {
        if w[1].0 < w[0].0 || w[1].1 < w[0].1 {
            return None;
        }
    }
    for &(ap, bp) in anchors {
        if ap + syncmer_len > a_len || bp + syncmer_len > b_len {
            return None;
        }
    }

    // Walk: (0..first_anchor.a_pos, 0..first_anchor.b_pos) gap → BiWFA,
    // then anchor match block, then gap between anchors, ... then final
    // gap (last_anchor..a_len, last_anchor..b_len) → BiWFA.
    let mut paf_cigar: Vec<u8> = Vec::new();

    let push_gap = |paf_cigar: &mut Vec<u8>,
                    a_start: usize,
                    a_end: usize,
                    b_start: usize,
                    b_end: usize|
     -> bool {
        let a_gap = &a_seq[a_start..a_end];
        let b_gap = &b_seq[b_start..b_end];
        if a_gap.is_empty() && b_gap.is_empty() {
            return true;
        }
        if a_gap.is_empty() {
            for _ in 0..b_gap.len() {
                paf_cigar.push(b'D');
            }
            return true;
        }
        if b_gap.is_empty() {
            for _ in 0..a_gap.len() {
                paf_cigar.push(b'I');
            }
            return true;
        }
        // Align the gap pair with BiWFA; translate WFA→PAF CIGAR on the
        // fly (WFA's I consumes text=b=target, WFA's D consumes pattern
        // =a=query; PAF swaps those).
        let wfa_cig = with_pair_aligner(|aligner| {
            let status = aligner.align(a_gap, b_gap);
            if !matches!(status, AlignmentStatus::Completed) {
                return None;
            }
            Some(aligner.cigar().to_vec())
        });
        let wfa_cig = match wfa_cig {
            Some(c) => c,
            None => return false,
        };
        // Verify the gap CIGAR consumes what we expect after I/D swap.
        let mut q_consumed = 0usize;
        let mut t_consumed = 0usize;
        for &op in &wfa_cig {
            match op {
                b'M' | b'=' | b'X' => {
                    q_consumed += 1;
                    t_consumed += 1;
                    paf_cigar.push(op);
                }
                b'D' => {
                    // WFA D consumes pattern (a=query) → PAF I
                    q_consumed += 1;
                    paf_cigar.push(b'I');
                }
                b'I' => {
                    // WFA I consumes text (b=target) → PAF D
                    t_consumed += 1;
                    paf_cigar.push(b'D');
                }
                _ => {}
            }
        }
        if q_consumed != a_gap.len() || t_consumed != b_gap.len() {
            return false;
        }
        true
    };

    // Gap before first anchor (or whole span if no anchors).
    let first = anchors.first().copied().unwrap_or((a_len, b_len));
    if !push_gap(&mut paf_cigar, 0, first.0, 0, first.1) {
        return None;
    }

    // For each anchor: emit syncmer_len M ops, then the gap to the
    // next anchor (if any).
    for i in 0..anchors.len() {
        let (ap, bp) = anchors[i];
        for _ in 0..syncmer_len {
            paf_cigar.push(b'M');
        }
        let a_after = ap + syncmer_len;
        let b_after = bp + syncmer_len;
        let (next_a, next_b) = if i + 1 < anchors.len() {
            let (nap, nbp) = anchors[i + 1];
            // Next anchor must start at or after the current one ended.
            if nap < a_after || nbp < b_after {
                // Overlapping / out-of-order anchor — skip to end.
                return None;
            }
            (nap, nbp)
        } else {
            (a_len, b_len)
        };
        if !push_gap(&mut paf_cigar, a_after, next_a, b_after, next_b) {
            return None;
        }
    }

    // Left-align the final CIGAR.
    let paf_cigar = crate::syng_graph_norm::left_align_indels(&paf_cigar, a_seq, b_seq);
    let (matches, mismatches, ins, del) = cigar_stats(&paf_cigar);
    let block_len = matches + mismatches + ins + del;
    if block_len == 0 {
        return None;
    }
    // Final sanity check.
    let q_consumed = (matches + mismatches + ins) as usize;
    let t_consumed = (matches + mismatches + del) as usize;
    if q_consumed != a_len || t_consumed != b_len {
        log::debug!(
            "anchor_seeded: dropping {}/{}: q={} t={} vs a_len={} b_len={}",
            a_name,
            b_name,
            q_consumed,
            t_consumed,
            a_len,
            b_len
        );
        return None;
    }
    let cigar_str = compact_cigar(&paf_cigar);
    Some(format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tcg:Z:{}\n",
        a_name, a_len, 0, a_len, paf_strand, b_name, b_len, 0, b_len, matches, block_len, cigar_str,
    ))
}

/// Use SweepGA's scaffold-aware plane sweep to turn raw shared syncmer
/// positions for one sequence pair into supported colinear anchor chains.
///
/// The anchors are exact PAF-style records in the oriented pair coordinate
/// system. SweepGA builds scaffold chains from those records, filters the
/// scaffold chains with the scaffold plane sweep, then returns the member
/// records of the scaffold chains that survived. We keep all returned scaffold
/// members, not just the largest inner chain.
fn plane_sweep_pair_syncmer_anchors(
    a_name: &str,
    b_name: &str,
    anchors: &[(usize, usize)],
    syncmer_len: usize,
) -> Vec<(usize, usize)> {
    if anchors.is_empty() || syncmer_len == 0 {
        return Vec::new();
    }

    use std::collections::BTreeMap;
    use sweepga::mapping::ChainStatus;
    use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, RecordMeta, ScoringFunction};

    let syncmer_len_u64 = syncmer_len as u64;
    let mut sorted = anchors.to_vec();
    sorted.sort_unstable();
    sorted.dedup();

    let records: Vec<RecordMeta> = sorted
        .iter()
        .enumerate()
        .map(|(rank, &(ap, bp))| RecordMeta {
            rank,
            query_name: a_name.to_string(),
            target_name: b_name.to_string(),
            query_start: ap as u64,
            query_end: ap as u64 + syncmer_len_u64,
            target_start: bp as u64,
            target_end: bp as u64 + syncmer_len_u64,
            block_length: syncmer_len_u64,
            identity: 1.0,
            matches: syncmer_len_u64,
            alignment_length: syncmer_len_u64,
            strand: '+',
            chain_id: None,
            chain_status: ChainStatus::Unassigned,
            discard: false,
            overlapped: false,
        })
        .collect();

    let filter_config = FilterConfig {
        chain_gap: 0,
        min_block_length: syncmer_len_u64,
        mapping_filter_mode: FilterMode::ManyToMany,
        mapping_max_per_query: None,
        mapping_max_per_target: None,
        plane_sweep_secondaries: 0,
        scaffold_filter_mode: FilterMode::OneToOne,
        scaffold_max_per_query: Some(1),
        scaffold_max_per_target: Some(1),
        overlap_threshold: 0.95,
        sparsity: 1.0,
        no_merge: true,
        scaffold_gap: 50_000,
        min_scaffold_length: syncmer_len_u64.saturating_mul(5),
        scaffold_overlap_threshold: 0.5,
        scaffold_max_deviation: syncmer_len_u64.saturating_mul(10),
        prefix_delimiter: '#',
        skip_prefix: true,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: 0.0,
        min_scaffold_identity: 0.0,
    };

    let passing = match PafFilter::new(filter_config)
        .with_keep_self(true)
        .with_scaffolds_only(true)
        .apply_filters(records)
    {
        Ok(passing) => passing,
        Err(err) => {
            log::debug!("syng_graph: SweepGA anchor plane sweep failed: {err}");
            return Vec::new();
        }
    };

    let mut by_chain: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
    for (rank, meta) in passing {
        let Some(&anchor) = sorted.get(rank) else {
            continue;
        };
        if meta.chain_status == ChainStatus::Scaffold {
            let chain_id = meta.chain_id.unwrap_or_else(|| format!("rank_{rank}"));
            by_chain.entry(chain_id).or_default().push(anchor);
        }
    }

    let mut anchors = Vec::new();
    for (_chain_id, mut chain) in by_chain {
        chain.sort_unstable();
        chain.dedup();
        anchors.extend(chain);
    }
    anchors.sort_unstable();
    anchors.dedup();

    let mut out = Vec::with_capacity(anchors.len());
    let mut last_a_end = 0usize;
    let mut last_b_end = 0usize;
    for (ap, bp) in anchors {
        if out.is_empty() || (ap >= last_a_end && bp >= last_b_end) {
            last_a_end = ap.saturating_add(syncmer_len);
            last_b_end = bp.saturating_add(syncmer_len);
            out.push((ap, bp));
        }
    }
    out
}

/// Build a GFA from a set of named sequences using pairwise BiWFA for
/// alignments, feeding the resulting PAF through the existing seqwish
/// induction pipeline in `crate::commands::graph::induce_graph_from_alignment`.
///
/// v0 uses naive all-pairs. Sparsification via
/// `sweepga::knn_graph::extract_tree_pairs_from_matrix` lands in v1.
/// Strand grooming and indel left-align land in v2. Returned GFA is
/// pre-normalization — callers apply `normalize_and_sort` if desired.
pub fn build_gfa_syng_native_from_sequences(
    seqs: &[(String, Vec<u8>)],
    config: &crate::commands::graph::GraphBuildConfig,
) -> std::io::Result<String> {
    if seqs.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }
    // Pair selection. Below a small cutoff, all-pairs is cheap and
    // avoids the mash-sketch cost. Above, delegate to sweepga's
    // tree-kNN sparsifier (k_near=3 backbone, k_far=1 novelty,
    // random=0.01 completeness).
    const ALL_PAIRS_CUTOFF: usize = 16;
    let paf_content = if seqs.len() <= ALL_PAIRS_CUTOFF {
        all_pairs_paf(seqs)
    } else {
        sparse_pairs_paf(seqs, 3, 1, 0.01)
    };
    build_gfa_from_paf_and_sequences(seqs, &paf_content, config)
}

/// Feed an arbitrary pre-generated PAF through the seqwish induction
/// pipeline. Exposed so callers with their own alignment source
/// (e.g. `build_paf_anchor_seeded`) can reuse the tail of the pipeline.
pub fn build_gfa_from_paf_and_sequences(
    seqs: &[(String, Vec<u8>)],
    paf_content: &str,
    config: &crate::commands::graph::GraphBuildConfig,
) -> std::io::Result<String> {
    use std::io::{BufWriter, Write};
    if seqs.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }
    let raw_paf_records = paf_content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .count();
    let debug_dir = std::env::var_os("IMPG_CRUSH_DEBUG_DIR").map(|root| {
        let id = DEBUG_GRAPH_BUILD_ID.fetch_add(1, Ordering::Relaxed);
        let dir = PathBuf::from(root).join(format!("graph_build_{id:04}"));
        let _ = std::fs::create_dir_all(&dir);
        dir
    });

    // Write combined FASTA to a temp file. One entry per input sequence;
    // seqwish indexes by the header name, so names must match the PAF's
    // query/target columns exactly.
    let combined_fasta = tempfile::Builder::new()
        .suffix(".fa")
        .tempfile()
        .map_err(|e| std::io::Error::other(format!("Failed to create temp FASTA: {}", e)))?;
    {
        let mut w = BufWriter::new(&combined_fasta);
        for (name, seq) in seqs {
            writeln!(w, ">{}", name)?;
            w.write_all(seq)?;
            writeln!(w)?;
        }
        w.flush()?;
    }
    if let Some(dir) = &debug_dir {
        let _ = std::fs::copy(combined_fasta.path(), dir.join("combined.fa"));
        let _ = write_path_spelling_debug(dir, seqs);
    }

    let mut paf_file = tempfile::Builder::new()
        .suffix(".paf")
        .tempfile()
        .map_err(|e| std::io::Error::other(format!("Failed to create temp PAF: {}", e)))?;
    paf_file.write_all(paf_content.as_bytes())?;
    paf_file.flush()?;
    if let Some(dir) = &debug_dir {
        let _ = std::fs::write(dir.join("raw.paf"), paf_content);
        let _ = write_paf_stage_debug(dir, "raw", seqs, paf_file.path());
    }

    let total_bases = seqs.iter().map(|(_, seq)| seq.len() as u64).sum::<u64>();
    let avg_seq_len = if seqs.is_empty() {
        0
    } else {
        total_bases / seqs.len() as u64
    };
    let filtered_paf = crate::commands::graph::filter_generated_paf(paf_file, avg_seq_len, config)?;
    if let Some(dir) = &debug_dir {
        let _ = std::fs::copy(filtered_paf.path(), dir.join("filtered.paf"));
        let _ = std::fs::copy(filtered_paf.path(), dir.join("final.paf"));
        let _ = write_paf_stage_debug(dir, "filtered", seqs, filtered_paf.path());
        let _ = write_paf_stage_debug(dir, "final", seqs, filtered_paf.path());
    }
    let paf_match_stats = filtered_paf_exact_match_run_stats(seqs, filtered_paf.path())?;

    // Optional adaptive seqwish floor. Local replacement induction defaults to
    // `min_match_len == 1`; callers that explicitly request adaptive behaviour
    // start from their configured floor and allow this wrapper to lower it to
    // observed local sequence/CIGAR evidence before invoking seqwish.
    let min_seq_len = seqs.iter().map(|(_, s)| s.len() as u64).min().unwrap_or(0);
    let mut effective_config = config.clone();
    if effective_config.adaptive_min_match_len {
        if min_seq_len > 0 && min_seq_len < effective_config.min_match_len {
            log::info!(
                "crush adaptive min-match: clamping seqwish min_match_len {} -> {} (shortest input traversal is {} bp)",
                effective_config.min_match_len,
                min_seq_len,
                min_seq_len,
            );
            effective_config.min_match_len = min_seq_len;
        }
        if paf_match_stats.max_exact_match_run > 0
            && paf_match_stats.max_exact_match_run < effective_config.min_match_len
        {
            log::info!(
                "crush adaptive min-match: clamping seqwish min_match_len {} -> {} (filtered PAF has {} record(s), {} with cg:Z, longest exact run seqwish can index is {} bp)",
                effective_config.min_match_len,
                paf_match_stats.max_exact_match_run,
                paf_match_stats.records,
                paf_match_stats.records_with_cigar,
                paf_match_stats.max_exact_match_run,
            );
            effective_config.min_match_len = paf_match_stats.max_exact_match_run;
        }
    }
    if paf_match_stats.records > 0 && paf_match_stats.max_exact_match_run == 0 {
        log::warn!(
            "crush seqwish min-match diagnostic: filtered PAF has {} record(s) but no indexable exact match run; names_missing={}, length_mismatches={}",
            paf_match_stats.records,
            paf_match_stats.name_mismatches,
            paf_match_stats.length_mismatches,
        );
    }
    log::debug!(
        "crush seqwish input summary: sequences={}, raw_paf_records={}, filtered_paf_records={}, parsed_filtered_records={}, max_exact_match_run={}, effective_min_match_len={}",
        seqs.len(),
        raw_paf_records,
        paf_match_stats.records,
        paf_match_stats.parsed_records,
        paf_match_stats.max_exact_match_run,
        effective_config.min_match_len,
    );

    // Hand off to the shared seqwish induction pipeline.
    let num_seqs = seqs.len();
    let num_genomes = seqs
        .iter()
        .filter_map(|(n, _)| n.split('#').next())
        .collect::<std::collections::HashSet<_>>()
        .len();
    let aln_result = crate::commands::graph::AlignmentResult {
        combined_fasta,
        filtered_paf,
        num_sequences: num_seqs,
        num_genomes,
    };
    let mut gfa_buf: Vec<u8> = Vec::new();
    crate::commands::graph::induce_graph_from_alignment(
        aln_result,
        &mut gfa_buf,
        &effective_config,
    )?;
    let gfa_text = String::from_utf8(gfa_buf).map_err(|e| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Invalid UTF-8 in GFA: {}", e),
        )
    })?;
    if let Some(dir) = &debug_dir {
        let (segments, segment_bp, links, paths) = gfa_line_stats(&gfa_text);
        let summary = format!(
            concat!(
                "sequences\t{}\n",
                "total_sequence_bp\t{}\n",
                "raw_paf_records\t{}\n",
                "filtered_paf_records\t{}\n",
                "parsed_filtered_paf_records\t{}\n",
                "filtered_paf_records_with_cigar\t{}\n",
                "paf_name_mismatches\t{}\n",
                "paf_length_mismatches\t{}\n",
                "max_exact_match_run\t{}\n",
                "adaptive_min_match_len\t{}\n",
                "effective_min_match_len\t{}\n",
                "seqwish_segments\t{}\n",
                "seqwish_segment_bp\t{}\n",
                "seqwish_links\t{}\n",
                "seqwish_paths\t{}\n"
            ),
            seqs.len(),
            total_bases,
            raw_paf_records,
            paf_match_stats.records,
            paf_match_stats.parsed_records,
            paf_match_stats.records_with_cigar,
            paf_match_stats.name_mismatches,
            paf_match_stats.length_mismatches,
            paf_match_stats.max_exact_match_run,
            effective_config.adaptive_min_match_len,
            effective_config.min_match_len,
            segments,
            segment_bp,
            links,
            paths,
        );
        let _ = std::fs::write(dir.join("seqwish.gfa"), gfa_text.as_bytes());
        let _ = std::fs::write(dir.join("summary.tsv"), summary);
    }
    Ok(gfa_text)
}

/// Anchor-seeded driver: re-queries syng from a seed interval to recover
/// anchor chains across the partition's members, then builds PAF via
/// `anchor_seeded_biwfa_paf` on pairs where we have usable anchors,
/// falling back to `pairwise_biwfa_paf` on pairs where anchors aren't
/// derivable.
///
/// `members` are pre-fetched, strand-groomed sequences keyed by their
/// `(chrom, fwd_start, fwd_end, strand)` interval. Anchor positions
/// from syng are in full-chromosome coordinates on the target side,
/// so we translate to fetched-sequence-local coords by offsetting
/// by `fwd_start` (for '+' strand members; '-' strand members are
/// skipped for anchor-seeded and fall back to full-pair).
///
/// This is conservative but correct:
/// - '+' strand members with fresh chains overlapping our interval → use anchors
/// - '-' strand members or members with no usable chain → full-pair BiWFA
pub fn build_paf_anchor_seeded(
    members: &[(String, Vec<u8>, Member)],
    seed_chrom: &str,
    seed_start: u64,
    seed_end: u64,
    syng_index: &crate::syng::SyngIndex,
    syng_padding: u64,
    max_depth: u16,
    k_near: usize,
    k_far: usize,
    random_fraction: f64,
) -> String {
    use rayon::prelude::*;
    use std::collections::HashMap;

    if members.len() < 2 {
        return String::new();
    }

    let t_total = std::time::Instant::now();
    let t_query = std::time::Instant::now();
    let chains = match crate::syng_transitive::query_transitive_with_anchors(
        syng_index,
        seed_chrom,
        seed_start,
        seed_end,
        syng_padding,
        max_depth.max(1),
    ) {
        Ok(cs) => cs,
        Err(e) => {
            log::debug!(
                "syng re-query failed: {}; falling back to full-pair BiWFA",
                e
            );
            return sparse_pairs_paf(
                &members
                    .iter()
                    .map(|(n, s, _)| (n.clone(), s.clone()))
                    .collect::<Vec<_>>(),
                k_near,
                k_far,
                random_fraction,
            );
        }
    };
    let dt_query = t_query.elapsed();

    let syncmer_len = (syng_index.params.w + syng_index.params.k) as usize;

    // Index chains by their (genome, start, end) — we match a partition
    // member to its chain by overlap: the member's interval should
    // contain (or be contained by) the chain's interval. For simplicity
    // we require exact match on (genome, start, end); members with
    // no exact-match chain fall back to full-pair.
    // Index chains by genome; multiple paralog chains may land on the
    // same genome, so we keep a list and pick the best-overlapping one
    // per member at lookup time.
    let mut chains_by_genome: HashMap<String, Vec<&crate::syng::HomologousIntervalWithAnchors>> =
        HashMap::new();
    for c in &chains {
        chains_by_genome
            .entry(c.genome.clone())
            .or_default()
            .push(c);
    }

    // Diagnostic counters — emitted at end.
    let counter_anchor = std::sync::atomic::AtomicUsize::new(0);
    let counter_fallback = std::sync::atomic::AtomicUsize::new(0);
    let counter_skipped = std::sync::atomic::AtomicUsize::new(0);

    let t_pairs = std::time::Instant::now();
    // Pair selection uses mash distance from raw sequences. We tested
    // replacing this with an anchor-count Jaccard derived from the
    // fresh chains (free, ~0.01s instead of mash's ~1.5s) but that
    // selected more pairs lacking shared anchors, pushing them to
    // full-pair BiWFA fallback and growing align cost by ~3s — net
    // negative. The mash sketch cost is buying us better pair quality
    // (pairs whose sequences are similar tend to share anchors).
    let raw: Vec<Vec<u8>> = members.iter().map(|(_, s, _)| s.clone()).collect();
    let pairs = sweepga::knn_graph::extract_tree_pairs(
        &raw,
        k_near,
        k_far,
        random_fraction,
        &sweepga::knn_graph::MashParams::default(),
    );
    let dt_pairs = t_pairs.elapsed();

    let t_align = std::time::Instant::now();
    let lines: Vec<String> = pairs
        .par_iter()
        .filter_map(|&(i, j)| {
            let (a_name, a_seq, a_meta) = &members[i];
            let (b_name, b_seq, b_meta) = &members[j];

            // For syng-partitioned input, member.strand is always '+' since
            // partition emits forward intervals; this is just a guard.
            if a_meta.strand == '+' && b_meta.strand == '+' {
                let ac = best_overlapping_chain(
                    chains_by_genome.get(&a_meta.chrom),
                    a_meta.fwd_start,
                    a_meta.fwd_end,
                );
                let bc = best_overlapping_chain(
                    chains_by_genome.get(&b_meta.chrom),
                    b_meta.fwd_start,
                    b_meta.fwd_end,
                );
                if let (Some(ac), Some(bc)) = (ac, bc) {
                    // Pair-relative strand: '+' if both chains agree, '-' if
                    // they disagree. Same-strand pairs (both '+' or both '-')
                    // align directly because both members are in the same
                    // relative orientation w.r.t. each other (even if both
                    // are flipped w.r.t. the seed). Mixed-strand pairs need
                    // one side reverse-complemented.
                    let pair_strand_relative = if ac.strand == bc.strand { '+' } else { '-' };

                    let a_start_i = a_meta.fwd_start as i64;
                    let b_start_i = b_meta.fwd_start as i64;
                    let a_len = a_seq.len() as i64;
                    let b_len = b_seq.len() as i64;
                    let scl = syncmer_len as i64;

                    // A's anchors → A-local fetched coords.
                    let mut a_by_node: HashMap<u32, i64> = HashMap::new();
                    for a in &ac.anchors {
                        a_by_node.insert(a.node_id, a.target_pos as i64 - a_start_i);
                    }
                    // B's anchors → B-local coords; if pair is '-', use
                    // rev-comp position.
                    let mut pair_anchors: Vec<(usize, usize)> = Vec::new();
                    for b in &bc.anchors {
                        if let Some(&a_local) = a_by_node.get(&b.node_id) {
                            let b_fwd_local = b.target_pos as i64 - b_start_i;
                            let b_local = if pair_strand_relative == '+' {
                                b_fwd_local
                            } else {
                                // After rev-comp(B): the syncmer occupying
                                // [b_fwd_local, b_fwd_local + scl) on forward
                                // ends up at [b_len - b_fwd_local - scl,
                                // b_len - b_fwd_local) on the rev-comp.
                                b_len - b_fwd_local - scl
                            };
                            if a_local < 0
                                || b_local < 0
                                || a_local + scl > a_len
                                || b_local + scl > b_len
                            {
                                continue;
                            }
                            pair_anchors.push((a_local as usize, b_local as usize));
                        }
                    }
                    let filtered_anchors = plane_sweep_pair_syncmer_anchors(
                        a_name,
                        b_name,
                        &pair_anchors,
                        syncmer_len,
                    );
                    if !filtered_anchors.is_empty() {
                        // For mixed-strand pairs, rev-comp B before BiWFA
                        // and emit PAF strand '-'. Same-strand pairs use
                        // forward B and strand '+'.
                        let b_oriented_owned: Vec<u8>;
                        let b_oriented: &[u8] = if pair_strand_relative == '-' {
                            b_oriented_owned = crate::graph::reverse_complement(b_seq);
                            &b_oriented_owned
                        } else {
                            b_seq
                        };
                        if let Some(line) = anchor_seeded_biwfa_paf_strand(
                            a_name,
                            a_seq,
                            b_name,
                            b_oriented,
                            &filtered_anchors,
                            syncmer_len,
                            pair_strand_relative,
                        ) {
                            counter_anchor.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            return Some(line);
                        }
                    } else if !pair_anchors.is_empty() {
                        counter_skipped.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        return None;
                    }
                }
            }
            // Fallback: full-pair BiWFA on forward sequences (PAF '+').
            counter_fallback.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            pairwise_biwfa_paf(a_name, a_seq, b_name, b_seq)
        })
        .collect();

    let dt_align = t_align.elapsed();
    let dt_total = t_total.elapsed();
    let anchor_used = counter_anchor.load(std::sync::atomic::Ordering::Relaxed);
    let fallback_used = counter_fallback.load(std::sync::atomic::Ordering::Relaxed);
    let skipped = counter_skipped.load(std::sync::atomic::Ordering::Relaxed);
    log::info!(
        "syng_graph: PAF total {:.2}s = query {:.2}s + pair-select {:.2}s + align {:.2}s | anchor-seeded {} / fallback {} / skipped {} pairs ({} members, {} chains)",
        dt_total.as_secs_f64(),
        dt_query.as_secs_f64(),
        dt_pairs.as_secs_f64(),
        dt_align.as_secs_f64(),
        anchor_used, fallback_used, skipped, members.len(), chains.len(),
    );

    lines.concat()
}

/// Pick the chain whose [start, end] has the largest overlap with [member_start, member_end].
/// Returns None if no chain provided or no overlap.
fn best_overlapping_chain<'a>(
    chains: Option<&'a Vec<&'a crate::syng::HomologousIntervalWithAnchors>>,
    member_start: u64,
    member_end: u64,
) -> Option<&'a crate::syng::HomologousIntervalWithAnchors> {
    let chains = chains?;
    let mut best: Option<(&crate::syng::HomologousIntervalWithAnchors, u64)> = None;
    for &c in chains {
        let ov_start = c.start.max(member_start);
        let ov_end = c.end.min(member_end);
        if ov_end > ov_start {
            let overlap = ov_end - ov_start;
            if best.map_or(true, |(_, b)| overlap > b) {
                best = Some((c, overlap));
            }
        }
    }
    best.map(|(c, _)| c)
}

/// Compact per-member metadata used by `build_paf_anchor_seeded`.
#[derive(Clone, Debug)]
pub struct Member {
    pub chrom: String,
    pub fwd_start: u64,
    pub fwd_end: u64,
    pub strand: char,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compact_cigar_basic() {
        assert_eq!(compact_cigar(b""), "");
        assert_eq!(compact_cigar(b"MMM"), "3M");
        assert_eq!(compact_cigar(b"MMMIMMM"), "3M1I3M");
        assert_eq!(compact_cigar(b"MIDXM"), "1M1I1D1X1M");
        assert_eq!(compact_cigar(b"MMMMMMMMMM"), "10M");
    }

    #[test]
    fn cigar_stats_basic() {
        assert_eq!(cigar_stats(b"MMMIMMM"), (6, 0, 1, 0));
        assert_eq!(cigar_stats(b"MMXIMD"), (3, 1, 1, 1));
        assert_eq!(cigar_stats(b""), (0, 0, 0, 0));
    }

    #[test]
    fn biwfa_identical_sequences() {
        let a = b"ACGTACGTACGT".to_vec();
        let line = pairwise_biwfa_paf("a", &a, "b", &a).unwrap();
        // All matches — 12M in CIGAR, matches == block_len == 12.
        assert!(line.contains("\t12\t12\t"));
        assert!(line.contains("cg:Z:12M"));
    }

    #[test]
    fn biwfa_with_substitution() {
        let a = b"ACGTACGTACGT".to_vec();
        let b = b"ACGTATGTACGT".to_vec(); // one mismatch at position 5
        let line = pairwise_biwfa_paf("a", &a, "b", &b).unwrap();
        // 11 matches, 1 mismatch, 12 total block.
        assert!(line.contains("\t11\t12\t"), "unexpected PAF: {line}");
    }

    #[test]
    fn biwfa_with_insertion() {
        // b has an extra base inserted
        let a = b"ACGTACGT".to_vec();
        let b = b"ACGTAACGT".to_vec();
        let line = pairwise_biwfa_paf("a", &a, "b", &b).unwrap();
        // PAF reports query length 8, target length 9.
        assert!(
            line.starts_with("a\t8\t0\t8\t+\tb\t9\t0\t9\t"),
            "unexpected PAF: {line}"
        );
        // CIGAR should contain a 1D (deletion on query = extra base in target).
        assert!(
            line.contains("1D") || line.contains("1I"),
            "no indel op: {line}"
        );
    }

    #[test]
    fn biwfa_empty_returns_none() {
        assert!(pairwise_biwfa_paf("a", b"", "b", b"ACGT").is_none());
        assert!(pairwise_biwfa_paf("a", b"ACGT", "b", b"").is_none());
    }

    #[test]
    fn anchor_seeded_identical_with_one_anchor() {
        // Two identical 12 bp sequences, one anchor of length 4 in the middle.
        // CIGAR should be 12M and PAF should report a=b=12, matches=12.
        let a = b"ACGTACGTACGT";
        let b = b"ACGTACGTACGT";
        let anchors = vec![(4, 4)]; // 4-bp anchor at positions 4..8
        let line = anchor_seeded_biwfa_paf("q", a, "t", b, &anchors, 4).unwrap();
        assert!(line.contains("\t12\t12\t"), "unexpected: {line}");
        assert!(line.contains("cg:Z:12M"), "unexpected cigar: {line}");
    }

    #[test]
    fn anchor_seeded_with_insertion_in_gap() {
        // Same query/target but target has one extra A in the first gap.
        let a: Vec<u8> = b"ACGTACGTACGT".to_vec();
        let b: Vec<u8> = b"ACGTAACGTACGT".to_vec(); // one extra A at position 4
                                                    // Anchor at (4,5): aligns a[4..8]=ACGT to b[5..9]=ACGT
        let anchors = vec![(4, 5)];
        let line = anchor_seeded_biwfa_paf("q", &a, "t", &b, &anchors, 4).unwrap();
        // query 12, target 13, one deletion from query side.
        assert!(
            line.starts_with("q\t12\t0\t12\t+\tt\t13\t0\t13"),
            "unexpected: {line}"
        );
    }

    #[test]
    fn anchor_seeded_zero_anchors_falls_through_to_full_biwfa() {
        // With no anchors, this is just a full-pair BiWFA.
        let a = b"ACGTACGTACGT";
        let b = b"ACGTACGTACGT";
        let line = anchor_seeded_biwfa_paf("q", a, "t", b, &[], 4).unwrap();
        assert!(line.contains("cg:Z:12M"));
    }

    #[test]
    fn anchor_seeded_rejects_non_colinear() {
        let a = b"ACGTACGTACGT";
        let b = b"ACGTACGTACGT";
        let anchors = vec![(4, 8), (8, 4)]; // b_pos goes backward
        assert!(anchor_seeded_biwfa_paf("q", a, "t", b, &anchors, 2).is_none());
    }

    #[test]
    fn pair_anchor_plane_sweep_keeps_supported_diagonal() {
        let mut anchors: Vec<(usize, usize)> = (0..6).map(|i| (i * 200, i * 200)).collect();
        anchors.push((1200, 100_000));

        let filtered = plane_sweep_pair_syncmer_anchors("q", "t", &anchors, 100);

        let expected: Vec<(usize, usize)> = (0..6).map(|i| (i * 200, i * 200)).collect();
        assert_eq!(filtered, expected);
    }

    #[test]
    fn pair_anchor_plane_sweep_keeps_indel_shifted_chain() {
        let anchors = vec![
            (0, 0),
            (200, 200),
            (400, 430),
            (600, 630),
            (800, 830),
            (1_000, 1_030),
        ];

        let filtered = plane_sweep_pair_syncmer_anchors("q", "t", &anchors, 100);

        assert_eq!(filtered, anchors);
    }

    #[test]
    fn pair_anchor_plane_sweep_keeps_multiple_scaffold_chains() {
        let mut anchors: Vec<(usize, usize)> = (0..5).map(|i| (i * 200, i * 200)).collect();
        anchors.extend((0..5).map(|i| (100_000 + i * 200, 100_000 + i * 200)));

        let filtered = plane_sweep_pair_syncmer_anchors("q", "t", &anchors, 100);

        assert_eq!(filtered, anchors);
    }

    #[test]
    fn build_gfa_three_identical_sequences() {
        // End-to-end: three identical short sequences should produce a
        // GFA with exactly one segment walked by three paths.
        let seqs = vec![
            ("sample1#0#chr1:0-12".to_string(), b"ACGTACGTACGT".to_vec()),
            ("sample2#0#chr1:0-12".to_string(), b"ACGTACGTACGT".to_vec()),
            ("sample3#0#chr1:0-12".to_string(), b"ACGTACGTACGT".to_vec()),
        ];
        let config = crate::commands::graph::GraphBuildConfig {
            num_threads: 1,
            ..Default::default()
        };
        let gfa = build_gfa_syng_native_from_sequences(&seqs, &config).unwrap();
        // Must have H line, at least one S line, and 3 P lines.
        assert!(gfa.contains("H\tVN:Z"), "missing header: {gfa}");
        let p_count = gfa.lines().filter(|l| l.starts_with("P\t")).count();
        assert_eq!(p_count, 3, "expected 3 paths, got {p_count}:\n{gfa}");
    }

    #[test]
    fn all_pairs_three_sequences() {
        let seqs = vec![
            ("a".to_string(), b"ACGTACGT".to_vec()),
            ("b".to_string(), b"ACGTACGT".to_vec()),
            ("c".to_string(), b"ACGTACGT".to_vec()),
        ];
        let paf = all_pairs_paf(&seqs);
        // 3 pairs: (a,b), (a,c), (b,c).
        assert_eq!(paf.lines().count(), 3);
    }

    #[test]
    fn seqwish_tail_does_not_rescue_records_dropped_by_explicit_filter() {
        let seq = b"ACGT".repeat(25);
        let seqs = vec![("q".to_string(), seq.clone()), ("t".to_string(), seq)];
        let paf = "q\t100\t0\t100\t+\tt\t100\t0\t100\t100\t100\t255\tcg:Z:100M\n";
        let config = crate::commands::graph::GraphBuildConfig {
            num_threads: 1,
            show_progress: false,
            no_filter: false,
            num_mappings: "many:many".to_string(),
            scaffold_filter: "many:many".to_string(),
            scaffold_mass: 0,
            min_map_length: 101,
            min_match_len: 1,
            adaptive_min_match_len: false,
            ..Default::default()
        };

        let gfa = build_gfa_from_paf_and_sequences(&seqs, paf, &config).unwrap();
        let mut support = std::collections::HashMap::<String, usize>::new();
        for line in gfa.lines().filter(|line| line.starts_with("P\t")) {
            if let Some(walk) = line.split('\t').nth(2) {
                for step in walk.split(',').filter(|step| !step.is_empty()) {
                    let node = step.trim_end_matches(|c| c == '+' || c == '-');
                    *support.entry(node.to_string()).or_default() += 1;
                }
            }
        }
        let shared_nodes = support.values().filter(|&&count| count > 1).count();
        assert_eq!(
            shared_nodes, 0,
            "the wrapper must not add raw PAF records back after min_map_length filtering:\n{gfa}"
        );
    }
}
