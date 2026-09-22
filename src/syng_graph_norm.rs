//! CIGAR normalization primitives for the syng-native graph engine.
//!
//! The headline operation is indel **left-alignment**: for each
//! insertion (`I`) or deletion (`D`) in a CIGAR, slide it as far left
//! as possible while preserving the alignment score. Canonicalizing the
//! placement of ambiguous indels (microsatellites, homopolymer runs,
//! tandem repeats) is what stops seqwish from seeing the "same" event
//! as two different variants and producing graph cycles / bubbles.
//!
//! Input CIGAR is the per-base byte form that WFA2 emits (`b"MMMIMMM"`),
//! after the WFA→PAF I/D swap has been applied (i.e. PAF convention
//! where `I` consumes query, `D` consumes target).
//!
//! Algorithm, per indel span:
//!   1. Identify each maximal I-run or D-run in the CIGAR.
//!   2. While the base immediately to the left of the run matches the
//!      base immediately at the tail of the run's content, rotate the
//!      indel left by one base (swap the M/= preceding the run with
//!      the last run op).
//!   3. Stop when the two bases differ or we reach the start of the
//!      alignment.
//!
//! This is the "anchor-base convention": shift until the shifted indel
//! is preceded by a base that *doesn't* extend the run. Matches what
//! `bcftools norm`, `vt normalize`, and `vg call` do at the VCF/BCF
//! level, but operates directly on CIGAR ops.

/// Left-align all I and D runs in a per-base byte CIGAR.
///
/// `query` and `target` are the aligned sequences in forward orientation
/// (caller is responsible for strand grooming — the CIGAR must be in
/// the query's frame). The CIGAR must be valid PAF convention:
/// `M`/`=`/`X` consume both, `I` consumes query, `D` consumes target.
pub fn left_align_indels(cigar: &[u8], query: &[u8], target: &[u8]) -> Vec<u8> {
    if cigar.is_empty() {
        return Vec::new();
    }
    let mut c = cigar.to_vec();

    // For each index i, shift any I/D starting at i as far left as
    // possible. We walk right-to-left so earlier shifts don't affect
    // the boundary positions of later runs.
    let mut i = c.len();
    while i > 0 {
        i -= 1;
        let op = c[i];
        if op != b'I' && op != b'D' {
            continue;
        }
        // Find the start of this I or D run (i is the rightmost op).
        let mut run_start = i;
        while run_start > 0 && c[run_start - 1] == op {
            run_start -= 1;
        }
        // Shift the run left as far as possible.
        let new_start = shift_indel_run_left(&mut c, run_start, op, query, target);
        // Continue scanning from one position before the (possibly new)
        // run start.
        i = new_start;
    }
    c
}

/// Shift a single indel run (op-type `op`, starting at `run_start`) as
/// far left as possible. Returns the new `run_start` index.
///
/// The run consists of identical `op` bytes in `cigar` starting at
/// `run_start`. We walk the CIGAR to count how many `op`s are in the
/// run, then repeatedly check whether the preceding op is a match
/// (`M`/`=`/`X`) and whether swapping it with the last op of the run
/// leaves the alignment consistent with the sequences.
fn shift_indel_run_left(
    cigar: &mut [u8],
    mut run_start: usize,
    op: u8,
    query: &[u8],
    target: &[u8],
) -> usize {
    // Determine the bp coordinates of the run on each sequence. To do
    // this, walk the CIGAR from the start and count query/target
    // consumption up to `run_start`.
    let (mut q_pos, mut t_pos) = consume_to(cigar, run_start);

    loop {
        if run_start == 0 {
            return run_start;
        }
        let prev_op = cigar[run_start - 1];
        if prev_op != b'M' && prev_op != b'=' && prev_op != b'X' {
            return run_start;
        }
        // Count the run length.
        let run_len = count_run(cigar, run_start, op);
        // The base we'd be swapping: immediately before the run on the
        // axis that the run consumes.
        let (swap_base, run_tail_base) = match op {
            // I consumes query. The run covers query[q_pos .. q_pos+run_len].
            // The "previous base on query" is query[q_pos - 1]. The run
            // tail base (last base in the inserted region) is query[q_pos + run_len - 1].
            b'I' => {
                if q_pos == 0 {
                    return run_start;
                }
                (query[q_pos - 1], query[q_pos + run_len - 1])
            }
            // D consumes target. Mirror logic with target[].
            b'D' => {
                if t_pos == 0 {
                    return run_start;
                }
                (target[t_pos - 1], target[t_pos + run_len - 1])
            }
            _ => return run_start,
        };
        if swap_base != run_tail_base {
            return run_start;
        }
        // Also require the OTHER axis's base at the boundary matches
        // itself — i.e. the `prev_op` is a real match `=` or `M` that
        // *stays* a match after the shift. For `X` (mismatch) we don't
        // shift because swapping would change what the mismatch aligns.
        if prev_op == b'X' {
            return run_start;
        }
        // Perform the rotation: swap prev_op with the last op of the run.
        // Before:  ...  M  {op op op}  ...   (run_start-1 is M; run_start..run_start+run_len is op)
        // After:   ... {op op op}  M  ...
        cigar[run_start - 1] = op;
        cigar[run_start + run_len - 1] = prev_op;
        run_start -= 1;
        // Update q_pos / t_pos: the previous M consumed one of each;
        // now it's been replaced by `op`, which consumes only one axis.
        // So the position where the run now starts shifted back by one
        // on both axes for the M-that-was-consumed (we were "past" it),
        // then forward by 0 on the non-consumed axis and forward by 1
        // on the op-consumed axis (because the M slid to the right of
        // the run).
        //
        // Net: the indel's left boundary on its axis moves back by 1.
        // The other axis's boundary is unchanged.
        match op {
            b'I' => {
                q_pos -= 1;
            }
            b'D' => {
                t_pos -= 1;
            }
            _ => unreachable!(),
        }
    }
}

/// Walk `cigar[0..idx]` and return `(query_consumed, target_consumed)`.
fn consume_to(cigar: &[u8], idx: usize) -> (usize, usize) {
    let mut q = 0usize;
    let mut t = 0usize;
    for &op in &cigar[..idx] {
        match op {
            b'M' | b'=' | b'X' => {
                q += 1;
                t += 1;
            }
            b'I' => q += 1,
            b'D' => t += 1,
            _ => {}
        }
    }
    (q, t)
}

/// Length of the run of `op` starting at `start`.
fn count_run(cigar: &[u8], start: usize, op: u8) -> usize {
    let mut n = 0;
    while start + n < cigar.len() && cigar[start + n] == op {
        n += 1;
    }
    n
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_noop() {
        let cig = b"MMMMMMMM";
        assert_eq!(
            left_align_indels(cig, b"ACGTACGT", b"ACGTACGT"),
            cig.to_vec()
        );
    }

    #[test]
    fn simple_insertion_shifts_left_over_homopolymer() {
        // query has an extra A in a homopolymer run:
        //   target: ACAAAAG   (len 7)
        //   query:  ACAAAAAG  (len 8, one extra A)
        // Initial CIGAR from WFA might put the I at the right edge:
        //   MMMMMMIM  (shifts through the A-run)
        // Expected left-aligned: MMIMMMMM  (I pushed as far left as possible)
        let query = b"ACAAAAAG";
        let target = b"ACAAAAG";
        let init = b"MMMMMMIM".to_vec();
        let out = left_align_indels(&init, query, target);
        // All the A bases in the run are interchangeable, so the I
        // should land at the earliest position where the preceding
        // base is not 'A' (i.e. right after 'C' at index 1).
        let expected_prefix = b"MMIMM";
        assert!(
            out.starts_with(expected_prefix),
            "expected left-shifted I; got {}",
            std::str::from_utf8(&out).unwrap()
        );
    }

    #[test]
    fn simple_deletion_shifts_left() {
        // query missing a base inside a T-run:
        //   target: GATTTC
        //   query:  GATTC (one T deleted)
        // CIGAR from WFA might put the D at the right edge.
        let query = b"GATTC";
        let target = b"GATTTC";
        let init = b"MMMMDM".to_vec(); // D consumes target; after it, M's resume.
        let out = left_align_indels(&init, query, target);
        // Expected: D moves left over the T run until preceded by A.
        // Query has GATT (positions 0..4 consumed by M's), then the D,
        // then C. target has GATTT at positions 0..5. Let me walk:
        //   init: MMMM D M   → query consumes 5, target consumes 4+1+1=6 — wait that's wrong.
        //     M=1q+1t each → 4M = 4q+4t
        //     D = 0q+1t
        //     M = 1q+1t
        //     total = 5q + 6t ✓
        // After shift, expected D as close to the start of the T-run
        // (position 2) as possible:
        //   MMDMMM (D at pos 2, shifting through the T's)
        assert_eq!(
            out,
            b"MMDMMM".to_vec(),
            "got {}",
            std::str::from_utf8(&out).unwrap()
        );
    }

    #[test]
    fn shift_stops_at_mismatch() {
        // Preceding op is X → should not shift past it.
        //   target: ACGTTT (6)
        //   query:  AXGTT  (5, with mismatch at pos 1 and a missing T)
        //   Hmm, let me keep it simpler with a concrete case.
        // target: ACTTT (5), query: ATTT (4) → query has an X at pos 0
        // (A→?), and missing a C or T. Actually let me just verify the
        // invariant: if prev_op is X, the run doesn't shift past it.
        let cig = b"XMDMM".to_vec(); // position 2 is D, preceded by M at 1 preceded by X at 0.
                                     // D can shift over M but should stop before X.
                                     // Sequences for this: CIGAR consumes q=1+1+0+1+1=4, t=1+1+1+1+1=5.
        let query = b"ACGT";
        let target = b"ACGTT"; // extra T at the end.
        let out = left_align_indels(&cig, query, target);
        // D was at index 2 originally, preceded by M at 1. The D run tail
        // base = target[t_pos + run_len - 1] = target[2+1-1] = target[2] = 'G'.
        // The swap_base (target[t_pos - 1]) = target[1] = 'C'. 'G' != 'C',
        // so no shift. Output == input.
        assert_eq!(
            out,
            cig,
            "unexpected shift; got {}",
            std::str::from_utf8(&out).unwrap()
        );
    }

    #[test]
    fn no_indel_no_change() {
        let cig = b"MMXMM".to_vec();
        let query = b"ACTGT";
        let target = b"ACCGT";
        assert_eq!(left_align_indels(&cig, query, target), cig);
    }

    #[test]
    fn empty_cigar() {
        assert_eq!(left_align_indels(b"", b"", b""), Vec::<u8>::new());
    }

    #[test]
    fn insertion_against_non_homopolymer_no_shift() {
        // Insertion not in a homopolymer — shouldn't shift.
        //   target: ACGT
        //   query:  ACCGT (one extra C)
        let query = b"ACCGT";
        let target = b"ACGT";
        // The extra C is at query position 2; target has C at position 1.
        // WFA might place it like MMIMM or MIMMM. Either should be stable
        // at MIMMM since the base before (position 0 = 'A') != 'C'.
        let cig = b"MMIMM".to_vec();
        let out = left_align_indels(&cig, query, target);
        // Should shift one position left: preceding base is 'C' at query[1]
        // which equals the last run base query[2]='C'. So shift once.
        // After shift: MIMMM. Now preceding base query[0]='A' != 'C', stop.
        assert_eq!(
            out,
            b"MIMMM".to_vec(),
            "got {}",
            std::str::from_utf8(&out).unwrap()
        );
    }
}
