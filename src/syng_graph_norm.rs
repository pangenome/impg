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

    // Maintain exact (q_pos, t_pos) consumed by c[0..i] in O(L) total instead
    // of re-scanning from 0 on every indel run.
    let (mut q_pos, mut t_pos) = consume_to(&c, c.len());

    // For each index i, shift any I/D starting at i as far left as
    // possible. We walk right-to-left so earlier shifts don't affect
    // the boundary positions of later runs.
    let mut i = c.len();
    while i > 0 {
        i -= 1;
        let op = c[i];
        match op {
            b'I' => q_pos = q_pos.saturating_sub(1),
            b'D' => t_pos = t_pos.saturating_sub(1),
            b'M' | b'=' | b'X' => {
                q_pos = q_pos.saturating_sub(1);
                t_pos = t_pos.saturating_sub(1);
                continue;
            }
            _ => continue,
        }
        // Find the start of this I or D run (i is the rightmost op).
        let mut run_start = i;
        while run_start > 0 && c[run_start - 1] == op {
            run_start -= 1;
            match op {
                b'I' => q_pos = q_pos.saturating_sub(1),
                b'D' => t_pos = t_pos.saturating_sub(1),
                _ => {}
            }
        }
        let mut run_len = i - run_start + 1;
        // Include any contiguous `op` bytes immediately to the right of `i`
        // that coalesced here when a later run shifted left.
        while run_start + run_len < c.len() && c[run_start + run_len] == op {
            run_len += 1;
        }
        // Shift the run left as far as possible; updates (q_pos, t_pos)
        // in place so they still reflect consumption of c[0..new_start].
        let new_start =
            shift_indel_run_left(&mut c, run_start, run_len, op, query, target, &mut q_pos, &mut t_pos);
        // Continue scanning from one position before the (possibly new)
        // run start.
        i = new_start;
    }
    c
}

/// Shift a single indel run (op-type `op`, length `run_len`, starting at
/// `run_start`) as far left as possible. Returns the new `run_start` index
/// and updates `q_pos`/`t_pos` to match `consume_to(cigar, new_start)`.
fn shift_indel_run_left(
    cigar: &mut [u8],
    mut run_start: usize,
    run_len: usize,
    op: u8,
    query: &[u8],
    target: &[u8],
    q_pos: &mut usize,
    t_pos: &mut usize,
) -> usize {
    if run_len == 0 {
        return run_start;
    }

    loop {
        if run_start == 0 {
            return run_start;
        }
        let prev_op = cigar[run_start - 1];
        // Only shift over true match ops (not X mismatches or another indel).
        if prev_op != b'M' && prev_op != b'=' {
            return run_start;
        }
        if *q_pos == 0 || *t_pos == 0 {
            return run_start;
        }
        // The base we'd be swapping: immediately before the run on the
        // axis that the run consumes, against the last base of the run.
        let (swap_base, run_tail_base) = match op {
            b'I' => {
                assert!(
                    *q_pos + run_len <= query.len(),
                    "left_align_indels: CIGAR exceeds query length {} (insertion ends at {})",
                    query.len(),
                    *q_pos + run_len
                );
                (query[*q_pos - 1], query[*q_pos + run_len - 1])
            }
            b'D' => {
                assert!(
                    *t_pos + run_len <= target.len(),
                    "left_align_indels: CIGAR exceeds target length {} (deletion ends at {})",
                    target.len(),
                    *t_pos + run_len
                );
                (target[*t_pos - 1], target[*t_pos + run_len - 1])
            }
            _ => return run_start,
        };
        if swap_base != run_tail_base {
            return run_start;
        }
        // Perform the rotation: swap prev_op with the last op of the run.
        // Before:  ...  M  {op op op}  ...   (run_start-1 is M; run_start..run_start+run_len is op)
        // After:   ... {op op op}  M  ...
        cigar[run_start - 1] = op;
        cigar[run_start + run_len - 1] = prev_op;
        run_start -= 1;
        // prev_op (M/=) consumed 1 on both axes and has moved to the right
        // of the indel run, so c[0..run_start] now consumes 1 less on both.
        *q_pos -= 1;
        *t_pos -= 1;
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

    #[test]
    fn multiple_indels_track_coordinates() {
        //   target: ACAAAATGGG  (len 10)
        //   query:  ACAAAAATGGGG (len 12): extra A + extra G
        //   raw:    MMMMMMIMMMMI (len 12: q=6+1+4+1=12, t=6+4=10)
        let target = b"ACAAAATGGG";
        let query = b"ACAAAAATGGGG";
        let cig = b"MMMMMMIMMMMI";
        let out = left_align_indels(cig, query, target);
        assert_eq!(out, b"MMIMMMMMIMMM");
    }

    #[test]
    #[should_panic(expected = "CIGAR exceeds query length")]
    fn overlong_insertion_cigar_panics() {
        left_align_indels(b"MMMI", b"AA", b"AA");
    }

    #[test]
    #[should_panic(expected = "CIGAR exceeds target length")]
    fn overlong_deletion_cigar_panics() {
        left_align_indels(b"MMMD", b"AA", b"AA");
    }

    #[test]
    fn coalescing_indel_runs_shift_left_together() {
        // Two separated insertions that coalesce into II= and then shift across
        // the leading match together to form II==.
        assert_eq!(left_align_indels(b"=I=I", b"ATAA", b"AA"), b"II==");
        assert_eq!(left_align_indels(b"=I=I=I", b"AAAAAA", b"AAA"), b"III===");
        assert_eq!(left_align_indels(b"=D=D", b"AA", b"ATAA"), b"DD==");
    }
}
