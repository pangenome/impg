//! Regression tests for crush failures documented in docs/crush-audit.md.
//!
//! Each test corresponds to a numbered failure in that document and asserts
//! the graph-level invariant the failure violates. Failures 1-4 are marked
//! `#[ignore]` to keep CI green; they are intentionally RED on HEAD (90ba74f).
//!
//! ## To run and verify they are red on HEAD:
//!
//! ```
//! source /home/erikg/impg/env.sh
//! cargo test --test test_crush_integration -- --ignored --nocapture 2>&1
//! ```
//!
//! ## Real data used:
//!
//! - `tests/test_data/crush/c4_slice_1500_3000.gfa` — 2942 segments / 64k bp
//!   extracted from the C4A blunt pangenome GFA (chr6:31891045-32123783, 465
//!   haplotypes, seqwish-k=311). Real chromosomal data, NOT synthetic.
//! - `/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/` — full C4A blunt
//!   GFA and known-good round-1 output (pre-0af1a4c baseline). Tests that
//!   require the full GFA skip automatically if the file is not present.

use impg::resolution::{path_sequences, resolve_gfa_bubbles, ResolutionConfig};
use std::collections::HashSet;
use std::path::Path;
use std::process::Command;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Locate the compiled `impg` binary.
fn impg_binary() -> Option<std::path::PathBuf> {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return Some(std::path::PathBuf::from(path));
    }
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    for candidate in [
        manifest_dir.join("target/release/impg"),
        manifest_dir.join("target/debug/impg"),
        // Fall back to the main repo binary on this machine
        std::path::PathBuf::from("/home/erikg/impg/target/release/impg"),
    ] {
        if candidate.exists() {
            return Some(candidate);
        }
    }
    None
}

/// Count S-lines and sum their sequence lengths in a GFA text.
fn gfa_segment_stats(gfa: &str) -> (usize, usize) {
    let mut count = 0usize;
    let mut bp = 0usize;
    for line in gfa.lines() {
        if line.starts_with("S\t") {
            let mut parts = line.splitn(4, '\t');
            parts.next(); // "S"
            parts.next(); // id
            if let Some(seq) = parts.next() {
                bp += seq.len();
                count += 1;
            }
        }
    }
    (count, bp)
}

/// Count S-lines and sum bp in a GFA file without reading it all into memory.
fn gfa_file_segment_stats(path: &str) -> (usize, usize) {
    use std::io::BufRead;
    let file = std::fs::File::open(path).expect("failed to open GFA file");
    let reader = std::io::BufReader::new(file);
    let mut count = 0usize;
    let mut bp = 0usize;
    for line in reader.lines() {
        let line = line.expect("failed to read line");
        if line.starts_with("S\t") {
            let mut parts = line.splitn(4, '\t');
            parts.next(); // "S"
            parts.next(); // id
            if let Some(seq) = parts.next() {
                bp += seq.len();
                count += 1;
            }
        }
    }
    (count, bp)
}

/// Return the set of sequences present in the S-lines.
fn gfa_segment_seqs(gfa: &str) -> HashSet<String> {
    gfa.lines()
        .filter(|l| l.starts_with("S\t"))
        .filter_map(|l| {
            let mut p = l.splitn(4, '\t');
            p.next(); // "S"
            p.next(); // id
            p.next().map(|s| s.to_owned()) // sequence
        })
        .collect()
}

/// Count duplicate sequences in a GFA file.
fn gfa_file_duplicate_seqs(path: &str) -> usize {
    use std::io::BufRead;
    let file = std::fs::File::open(path).expect("failed to open GFA file");
    let reader = std::io::BufReader::new(file);
    let mut seqs: HashSet<String> = HashSet::new();
    let mut total = 0usize;
    for line in reader.lines() {
        let line = line.expect("failed to read line");
        if line.starts_with("S\t") {
            let mut parts = line.splitn(4, '\t');
            parts.next(); // "S"
            parts.next(); // id
            if let Some(seq) = parts.next() {
                if !seqs.insert(seq.to_owned()) {
                    total += 1;
                }
            }
        }
    }
    total
}

// ---------------------------------------------------------------------------
// Failure 2 — No round-level quality gate: rounds that grow score are accepted
//
// docs/crush-audit.md §"Failure 2"
//
// Removed in commit 0af1a4c ("Accept validated crush replacements"):
//   - `pub const DEFAULT_MAX_ROUND_SCORE_GROWTH: f64 = 0.02`
//   - `fn round_quality_decision` and `enum RoundQualityDecision`
//   - The rollback that emitted "crush round N: rejecting M replacement(s): ..."
//
// Invariant violated: after accepting a round, segment-bp must not grow by
// more than 10% (the deleted threshold was 2%, we allow 10% here as a generous
// bound — the actual HEAD growth is 18.7% with one-round sweepga or 78.8% with
// the full canonical command including polish).
//
// Expected failure on HEAD (90ba74f): round 1 of the C4 command grows
// segment-bp by 18.7% (389 316 → 461 998) which exceeds the 10% bound.
// ---------------------------------------------------------------------------

#[test]
#[ignore = "RED on HEAD 90ba74f: round 1 grows segment-bp 18.7% (> 10% threshold); uses full C4 GFA (2 min)"]
fn c4_round1_sweep_quality_gate_rejects_score_growth() {
    // docs/crush-audit.md §Failure2 — maps to 0af1a4c removing round_quality_decision
    let Some(bin) = impg_binary() else {
        eprintln!("SKIP: impg binary not found");
        return;
    };
    let gfa_path =
        "/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa";
    if !Path::new(gfa_path).exists() {
        eprintln!("SKIP: {} not found", gfa_path);
        return;
    }
    let out_path = std::env::temp_dir().join("c4_round1_sweep_quality_gate.gfa");

    let (before_segs, before_bp) = gfa_file_segment_stats(gfa_path);
    eprintln!("Input: {} segments / {} bp", before_segs, before_bp);

    // Parameters match the canonical C4 command from docs/c4-crush-handoff.md.
    // max-iterations=1 → exactly round 1 (the round that fails the quality gate).
    let run = Command::new(&bin)
        .args([
            "crush",
            "--gfa", gfa_path,
            "--output", out_path.to_str().unwrap(),
            "--method", "sweepga",
            "--sweepga-aligner", "fastga",
            "--seqwish-k", "311",
            "--min-traversal-len", "5000",
            "--max-iterations", "1",
            "-v", "2",
        ])
        .output()
        .expect("failed to run impg crush");

    let stderr = String::from_utf8_lossy(&run.stderr);
    eprintln!("stderr (tail):\n{}", stderr.lines().rev().take(5).collect::<Vec<_>>().into_iter().rev().collect::<Vec<_>>().join("\n"));
    assert!(
        run.status.success(),
        "impg crush failed unexpectedly: {}",
        stderr
    );

    let (after_segs, after_bp) = gfa_file_segment_stats(out_path.to_str().unwrap());
    let growth = (after_bp as f64 - before_bp as f64) / before_bp as f64;
    eprintln!(
        "After round 1: {} segments / {} bp (Δsegs={:+}, Δbp={:+}, bp_growth={:.1}%)",
        after_segs, after_bp,
        after_segs as i64 - before_segs as i64,
        after_bp as i64 - before_bp as i64,
        growth * 100.0,
    );

    std::fs::remove_file(&out_path).ok();

    // The deleted quality gate (DEFAULT_MAX_ROUND_SCORE_GROWTH=0.02) would have
    // rejected this round. We assert a generous 10% bound — still clearly violated
    // by the 18.7% observed growth on HEAD.
    assert!(
        growth <= 0.10,
        "Failure 2: segment-bp grew {:.1}% in round 1 (before={}, after={}) — \
         the quality gate (DEFAULT_MAX_ROUND_SCORE_GROWTH=0.02) was removed in 0af1a4c; \
         this round should have been rejected",
        growth * 100.0,
        before_bp,
        after_bp
    );
}

// ---------------------------------------------------------------------------
// Failure 3 — Segment deduplication missing across replacements
//
// docs/crush-audit.md §"Failure 3"
//
// In `render_rewritten_graph` (src/resolution.rs:2810-2900), every
// `Replacement(plan_idx, node_idx)` node is assigned a fresh segment id via
// `next_unused_segment_id` with no sequence-deduplication. Two replacements
// that produce byte-identical segments emit two separate `S` lines.
//
// Invariant violated: no two S-lines may have the same sequence after crush.
//
// Expected failure on HEAD (90ba74f): the C4 round-1 output has tandem-
// paralog bubbles producing sequence-identical segments across parallel
// replacements.
// ---------------------------------------------------------------------------

#[test]
#[ignore = "RED on HEAD 90ba74f: multi-replacement round emits duplicate S-line sequences; uses full C4 GFA"]
fn c4_round1_render_emits_no_duplicate_segment_sequences() {
    // docs/crush-audit.md §Failure3 — render_rewritten_graph :2851 lacks seq->id dedup
    let Some(bin) = impg_binary() else {
        eprintln!("SKIP: impg binary not found");
        return;
    };
    let gfa_path =
        "/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa";
    if !Path::new(gfa_path).exists() {
        eprintln!("SKIP: {} not found", gfa_path);
        return;
    }
    let out_path = std::env::temp_dir().join("c4_round1_no_dup_seqs.gfa");

    let (before_segs, _) = gfa_file_segment_stats(gfa_path);

    let run = Command::new(&bin)
        .args([
            "crush",
            "--gfa", gfa_path,
            "--output", out_path.to_str().unwrap(),
            "--method", "sweepga",
            "--sweepga-aligner", "fastga",
            "--seqwish-k", "311",
            "--min-traversal-len", "5000",
            "--max-iterations", "1",
            "-v", "1",
        ])
        .output()
        .expect("failed to run impg crush");

    assert!(
        run.status.success(),
        "impg crush failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );

    let out_str = out_path.to_str().unwrap();
    let (after_segs, _) = gfa_file_segment_stats(out_str);
    let duplicates = gfa_file_duplicate_seqs(out_str);

    eprintln!(
        "Before: {} segs | After: {} segs | Duplicate seqs: {}",
        before_segs, after_segs, duplicates
    );

    std::fs::remove_file(&out_path).ok();

    // Invariant: no two S-lines should spell the same sequence after crush.
    // The failure is in render_rewritten_graph (src/resolution.rs:2851) which
    // assigns fresh IDs to replacement segments without deduplicating sequences.
    assert_eq!(
        duplicates, 0,
        "Failure 3: {} segment(s) share a sequence with another segment after round 1; \
         render_rewritten_graph does not deduplicate across parallel replacements \
         (see src/resolution.rs:2851)",
        duplicates
    );
}

// ---------------------------------------------------------------------------
// Failure 4 — Round-over-round candidate inflation
//
// docs/crush-audit.md §"Failure 4"
//
// Invariant violated: after 2 rounds, total segment-bp must not grow by more
// than 20% from the original input. The baseline grows only +8.4% in round 1
// then STOPS. On HEAD, after 2 rounds segment-bp reaches ~740k (+90%).
//
// Expected failure on HEAD (90ba74f): 2 rounds grow segment-bp >20% because
// round 2 starts from a bloated round-1 graph (+18.7%) and further inflates it.
// ---------------------------------------------------------------------------

#[test]
#[ignore = "RED on HEAD 90ba74f: segment-bp grows >20% across 2 rounds (candidate inflation); 15+ min"]
fn c4_round2_segment_bp_does_not_exceed_round1() {
    // docs/crush-audit.md §Failure4 — consequence of missing quality gate at :611-616
    let Some(bin) = impg_binary() else {
        eprintln!("SKIP: impg binary not found");
        return;
    };
    let gfa_path =
        "/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa";
    if !Path::new(gfa_path).exists() {
        eprintln!("SKIP: {} not found", gfa_path);
        return;
    }
    let out_path = std::env::temp_dir().join("c4_round2_inflation.gfa");

    let (_, before_bp) = gfa_file_segment_stats(gfa_path);

    // Run TWO rounds to observe the inflation cascade
    let run = Command::new(&bin)
        .args([
            "crush",
            "--gfa", gfa_path,
            "--output", out_path.to_str().unwrap(),
            "--method", "sweepga",
            "--sweepga-aligner", "fastga",
            "--seqwish-k", "311",
            "--min-traversal-len", "5000",
            "--max-iterations", "2",
            "-v", "1",
        ])
        .output()
        .expect("failed to run impg crush (2 rounds)");

    let stderr = String::from_utf8_lossy(&run.stderr);
    assert!(
        run.status.success(),
        "impg crush failed: {}",
        stderr
    );

    let (after_segs, after_bp) = gfa_file_segment_stats(out_path.to_str().unwrap());
    let growth = (after_bp as f64 - before_bp as f64) / before_bp as f64;
    eprintln!(
        "After 2 rounds: {} segs / {} bp (growth={:.1}%)",
        after_segs, after_bp,
        growth * 100.0,
    );

    std::fs::remove_file(&out_path).ok();

    // The baseline (pre-0af1a4c) kept segment-bp nearly flat: 389k → 422k (+8.4%)
    // in round 1, then REJECTED round 2. On HEAD, segment-bp reaches ~460k after
    // round 1 and continues to grow in round 2. We allow a generous 20% bound.
    assert!(
        growth <= 0.20,
        "Failure 4: segment-bp grew {:.1}% across 2 rounds (before={}, after={}) — \
         candidate inflation caused by accepting quality-degrading rounds; \
         the quality gate removed in 0af1a4c would have stopped this",
        growth * 100.0,
        before_bp,
        after_bp
    );
}

// ---------------------------------------------------------------------------
// Failure 1 — Canonical command does not finish within 30 minutes
//
// docs/crush-audit.md §"Failure 1"
//
// Invariant violated: the crush pipeline should complete within 360 seconds
// on the C4 canonical command. The baseline completed in 4:39.
//
// Expected failure on HEAD (90ba74f): the CLI hangs at round 3 candidate 8/8;
// this test will hit the 360s timeout and fail with a non-zero exit status.
// ---------------------------------------------------------------------------

#[test]
#[ignore = "RED on HEAD 90ba74f: canonical C4 command hangs >30 min; test uses 360s budget"]
fn c4_canonical_command_completes_within_budget() {
    // docs/crush-audit.md §Failure1 — hang at round 3 7/8 → 8/8 never finishes
    let Some(bin) = impg_binary() else {
        eprintln!("SKIP: impg binary not found");
        return;
    };
    let gfa_path =
        "/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa";
    if !Path::new(gfa_path).exists() {
        eprintln!("SKIP: {} not found", gfa_path);
        return;
    }

    let out_path = std::env::temp_dir().join("c4_canonical_test_out.gfa");

    // Known-good baseline completed in 4:39 (279s). We allow 360s (6 min).
    // On HEAD the command never finishes within 30 minutes.
    let run = Command::new("timeout")
        .args([
            "--signal=SIGTERM",
            "--kill-after=30s",
            "360s",
            bin.to_str().unwrap(),
            "crush",
            "--gfa", gfa_path,
            "--output", out_path.to_str().unwrap(),
            "--method", "sweepga",
            "--sweepga-aligner", "fastga",
            "--seqwish-k", "311",
            "--min-traversal-len", "5000",
            "--max-iterations", "until-done",
            "--polish-method", "smooth",
            "--polish-rounds", "until-done",
            "--polish-max-traversal-len", "10000",
            "--polish-max-median-traversal-len", "1000",
            "-v", "1",
        ])
        .output()
        .expect("failed to run impg crush under timeout");

    let stderr = String::from_utf8_lossy(&run.stderr);
    eprintln!("exit status: {}", run.status);
    eprintln!("stderr (last 5 lines):\n{}",
        stderr.lines().rev().take(5).collect::<Vec<_>>().into_iter().rev().collect::<Vec<_>>().join("\n"));

    std::fs::remove_file(&out_path).ok();

    assert!(
        run.status.success(),
        "Failure 1: impg crush did not complete within 360s on the C4 canonical command \
         (exit={:?}). On HEAD (90ba74f) round 3 hangs because rounds 1 and 2 were accepted \
         despite segment-bp growth (Failure 2), inflating subsequent candidates (Failure 4). \
         The known-good baseline (pre-0af1a4c) completed in 4:39.",
        run.status.code()
    );
}

// ---------------------------------------------------------------------------
// Companion test: small C4 slice with auto method should not degrade quality
//
// This test uses the committed test data slice and PASSES on HEAD because
// auto-SPOA on small bubbles (<500 bp) actually improves quality. It confirms
// the crush API is callable and path sequences are preserved.
// ---------------------------------------------------------------------------

#[test]
fn c4_slice_auto_crush_preserves_path_sequences() {
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let gfa_path = manifest_dir.join("tests/test_data/crush/c4_slice_1500_3000.gfa");
    if !gfa_path.exists() {
        eprintln!("SKIP: test data slice not found at {:?}", gfa_path);
        return;
    }

    let gfa = std::fs::read_to_string(&gfa_path).expect("failed to read C4 slice GFA");
    let before_seqs: std::collections::HashMap<_, _> =
        path_sequences(&gfa).unwrap().into_iter().collect();

    let resolved = resolve_gfa_bubbles(&gfa, &ResolutionConfig::default())
        .expect("resolve_gfa_bubbles failed on C4 slice");

    let after_seqs: std::collections::HashMap<_, _> =
        path_sequences(&resolved.gfa).unwrap().into_iter().collect();

    for (name, before_seq) in &before_seqs {
        let after_seq = after_seqs.get(name).expect("path disappeared after crush");
        assert_eq!(
            before_seq, after_seq,
            "path {} changed sequence after crush",
            name
        );
    }

    let (before_segs, before_bp) = gfa_segment_stats(&gfa);
    let (after_segs, after_bp) = gfa_segment_stats(&resolved.gfa);
    eprintln!(
        "slice crush: {} paths preserved, {} resolved, {} bailed",
        before_seqs.len(),
        resolved.stats.resolved,
        resolved.stats.bailed
    );
    eprintln!(
        "slice: segs {} → {}, bp {} → {}",
        before_segs, after_segs, before_bp, after_bp
    );
}
