//! THE biWFA ALIGNMENT HELPER (receipt-side assessment machinery for
//! the called-vs-truth sequence-QV stage).
//!
//! The likegt sequence-QV pattern (likegt/src/sequence_qv.rs, the
//! owner's named reference for this stage): a gap-affine WFA2
//! End2End alignment with match 0 / mismatch 4 / gap-open 6 /
//! gap-extend 2, Medium memory, no heuristic, AlignmentScope::
//! Alignment so the CIGAR is available; the per-base error is
//! (mismatches + insertions + deletions) / (total CIGAR columns).
//!
//! Protocol: one request per stdin line, "seqA\tseqB" (the sequences
//! are uppercase panel windows; they carry no tabs); one response per
//! request on stdout:
//!   "matches\tmismatches\tins\tdels"        on success
//!   "FAIL:<reason>"                        on any failure (never
//!        silently approximated: a failed alignment is reported, the
//!        caller names the locus and excludes it from the aggregate
//!        honestly)
//!
//! The CIGAR is walked against BOTH sequences byte-by-byte: match /
//! mismatch columns are decided by comparing the actual bases (exact
//! regardless of whether the library emits 'M', '=' or 'X'), and the
//! walk must consume both sequences exactly (End2End) — any
//! inconsistency is a FAIL, never a guess.

use std::io::{self, BufRead, Write};

use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, AlignmentStatus, MemoryMode,
};

fn main() {
    // The likegt penalties exactly (BIWFA_INTEGRATION.md /
    // sequence_qv.rs): mismatch 4, gap-opening 6, gap-extension 2,
    // match score 0, Medium memory, End2End span, Alignment scope
    // (CIGAR), no heuristic (with_penalties_and_memory_mode sets
    // wf_heuristic_strategy_wf_heuristic_none at this rev).
    let mut aligner = AffineWavefronts::with_penalties_and_memory_mode(
        0, 4, 6, 2, MemoryMode::Medium,
    );
    aligner.set_alignment_scope(AlignmentScope::Alignment);
    aligner.set_alignment_span(AlignmentSpan::End2End);

    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    for line in stdin.lock().lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                writeln!(out, "FAIL:stdin-read:{e}").unwrap();
                continue;
            }
        };
        let mut fields = line.splitn(2, '\t');
        let a = fields.next().unwrap_or("");
        let b = match fields.next() {
            Some(b) => b,
            None => {
                writeln!(out, "FAIL:malformed-request").unwrap();
                continue;
            }
        };
        if a.is_empty() || b.is_empty() {
            writeln!(out, "FAIL:empty-sequence").unwrap();
            continue;
        }
        let status = aligner.align(a.as_bytes(), b.as_bytes());
        let cigar = aligner.cigar();
        let mut matches = 0usize;
        let mut mismatches = 0usize;
        let mut ins = 0usize;
        let mut dels = 0usize;
        let mut i = 0usize; // consumed from a
        let mut j = 0usize; // consumed from b
        let mut ok = matches!(status, AlignmentStatus::Completed);
        for &op in cigar {
            match op {
                b'M' | b'=' | b'X' => {
                    if i >= a.len() || j >= b.len() {
                        ok = false;
                        break;
                    }
                    // Decide the column by the actual bases: exact
                    // for any op encoding the library chooses.
                    if a.as_bytes()[i] == b.as_bytes()[j] {
                        matches += 1;
                    } else {
                        mismatches += 1;
                    }
                    i += 1;
                    j += 1;
                }
                b'I' | b'D' => {
                    // A gap column: one indel column regardless of
                    // which side carries the gap (the QV convention
                    // counts gap columns in both the numerator and
                    // the denominator, side-agnostic). WFA2's
                    // convention (measured, and the same reading as
                    // likegt's comment): 'I' consumes the TEXT
                    // (b), 'D' consumes the PATTERN (a).
                    if op == b'I' {
                        ins += 1;
                        j += 1;
                    } else {
                        dels += 1;
                        i += 1;
                    }
                    if i > a.len() || j > b.len() {
                        ok = false;
                        break;
                    }
                }
                _ => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok || i != a.len() || j != b.len() {
            writeln!(
                out,
                "FAIL:cigar-walk:consumed_a={i}/{},consumed_b={j}/{}",
                a.len(),
                b.len()
            )
            .unwrap();
            out.flush().unwrap();
            continue;
        }
        writeln!(out, "{matches}\t{mismatches}\t{ins}\t{dels}").unwrap();
        out.flush().unwrap();
    }
}
