//! Syng-native graph engine primitives.
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
    let (matches, mismatches, ins, del) = cigar_stats(&cigar_bytes);
    let block_len = matches + mismatches + ins + del;
    if block_len == 0 {
        return None;
    }
    let cigar_str = compact_cigar(&cigar_bytes);
    let a_len = a_seq.len();
    let b_len = b_seq.len();
    Some(format!(
        "{}\t{}\t{}\t{}\t+\t{}\t{}\t{}\t{}\t{}\t{}\t255\tcg:Z:{}\n",
        a_name, a_len, 0, a_len,
        b_name, b_len, 0, b_len,
        matches, block_len,
        cigar_str,
    ))
}

/// Build a PAF string from all pairs of sequences (i < j). Simplest
/// possible sparsification: none. Intended for small partitions or as
/// a correctness baseline; real use goes through
/// `sweepga::knn_graph::extract_tree_pairs_from_matrix` with a distance
/// matrix derived from syng anchor counts.
pub fn all_pairs_paf(seqs: &[(String, Vec<u8>)]) -> String {
    let mut out = String::new();
    for i in 0..seqs.len() {
        for j in (i + 1)..seqs.len() {
            if let Some(line) = pairwise_biwfa_paf(
                &seqs[i].0, &seqs[i].1, &seqs[j].0, &seqs[j].1,
            ) {
                out.push_str(&line);
            }
        }
    }
    out
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
        assert!(line.starts_with("a\t8\t0\t8\t+\tb\t9\t0\t9\t"), "unexpected PAF: {line}");
        // CIGAR should contain a 1D (deletion on query = extra base in target).
        assert!(line.contains("1D") || line.contains("1I"), "no indel op: {line}");
    }

    #[test]
    fn biwfa_empty_returns_none() {
        assert!(pairwise_biwfa_paf("a", b"", "b", b"ACGT").is_none());
        assert!(pairwise_biwfa_paf("a", b"ACGT", "b", b"").is_none());
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
}
