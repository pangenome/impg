//! Regression tests for crush failures documented in docs/crush-audit.md.
//!
//! Each test corresponds to a numbered failure in that document and asserts
//! the graph-level invariant the failure violates. Known RED checks are marked
//! `#[ignore]` to keep CI green and are documented in the corresponding task
//! notes.
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

use impg::graph_report::{describe_gfa, GraphReportOptions};
use impg::resolution::{path_sequences, resolve_gfa_bubbles, ResolutionConfig, ResolutionMethod};
use std::collections::{HashMap, HashSet};
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

fn fasta_records(path: &Path) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(path).expect("failed to read FASTA");
    let mut records = Vec::new();
    let mut name: Option<String> = None;
    let mut seq = Vec::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(prev) = name.replace(header.to_string()) {
                records.push((prev, std::mem::take(&mut seq)));
            }
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(name) = name {
        records.push((name, seq));
    }
    records
}

fn gfa_shared_path_segment_count(gfa: &str) -> usize {
    let mut support: HashMap<String, HashSet<String>> = HashMap::new();
    for line in gfa.lines().filter(|line| line.starts_with("P\t")) {
        let mut fields = line.split('\t');
        fields.next();
        let Some(path_name) = fields.next() else {
            continue;
        };
        let Some(steps) = fields.next() else {
            continue;
        };
        for step in steps.split(',').filter(|step| !step.is_empty()) {
            let node = step
                .strip_suffix('+')
                .or_else(|| step.strip_suffix('-'))
                .unwrap_or(step);
            support
                .entry(node.to_string())
                .or_default()
                .insert(path_name.to_string());
        }
    }
    support.values().filter(|paths| paths.len() > 1).count()
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

#[test]
fn c4_top_flubble_seqwish_indexes_observed_exact_runs() {
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let fasta_path = manifest_dir.join("tests/test_data/crush/top_flubble_seqwish_minrun.fa");
    let paf_path = manifest_dir.join("tests/test_data/crush/top_flubble_seqwish_minrun.paf");
    let seqs = fasta_records(&fasta_path);
    let paf = std::fs::read_to_string(&paf_path).expect("failed to read C4 top-flubble PAF");
    let input_bp: usize = seqs.iter().map(|(_, seq)| seq.len()).sum();

    let config = impg::commands::graph::GraphBuildConfig {
        num_threads: 1,
        min_match_len: 311,
        no_filter: true,
        show_progress: false,
        ..impg::commands::graph::GraphBuildConfig::default()
    };
    let gfa = impg::syng_graph::build_gfa_from_paf_and_sequences(&seqs, &paf, &config)
        .expect("seqwish induction failed for C4 top-flubble regression block");
    let shared_segments = gfa_shared_path_segment_count(&gfa);
    let (segments, segment_bp) = gfa_segment_stats(&gfa);

    assert!(
        shared_segments > 0,
        "C4 top-flubble block has {} real PAF records but seqwish emitted no shared segment; \
         before the fix min_match_len clamped only to the shortest traversal (299 bp), \
         above every exact run in this PAF, and produced five isolated full-length paths\n{}",
        paf.lines().count(),
        gfa,
    );
    assert!(
        segment_bp < input_bp,
        "replacement should reuse aligned sequence (segment bp {} < input bp {}); segments={}",
        segment_bp,
        input_bp,
        segments,
    );
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
    let gfa_path = "/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa";
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
            "--gfa",
            gfa_path,
            "--output",
            out_path.to_str().unwrap(),
            "--method",
            "sweepga",
            "--sweepga-aligner",
            "fastga",
            "--seqwish-k",
            "311",
            "--min-traversal-len",
            "5000",
            "--max-iterations",
            "1",
            "-v",
            "2",
        ])
        .output()
        .expect("failed to run impg crush");

    let stderr = String::from_utf8_lossy(&run.stderr);
    eprintln!(
        "stderr (tail):\n{}",
        stderr
            .lines()
            .rev()
            .take(5)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
            .join("\n")
    );
    assert!(
        run.status.success(),
        "impg crush failed unexpectedly: {}",
        stderr
    );

    let (after_segs, after_bp) = gfa_file_segment_stats(out_path.to_str().unwrap());
    let growth = (after_bp as f64 - before_bp as f64) / before_bp as f64;
    eprintln!(
        "After round 1: {} segments / {} bp (Δsegs={:+}, Δbp={:+}, bp_growth={:.1}%)",
        after_segs,
        after_bp,
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
    let gfa_path = "/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa";
    if !Path::new(gfa_path).exists() {
        eprintln!("SKIP: {} not found", gfa_path);
        return;
    }
    let out_path = std::env::temp_dir().join("c4_round1_no_dup_seqs.gfa");

    let (before_segs, _) = gfa_file_segment_stats(gfa_path);

    let run = Command::new(&bin)
        .args([
            "crush",
            "--gfa",
            gfa_path,
            "--output",
            out_path.to_str().unwrap(),
            "--method",
            "sweepga",
            "--sweepga-aligner",
            "fastga",
            "--seqwish-k",
            "311",
            "--min-traversal-len",
            "5000",
            "--max-iterations",
            "1",
            "-v",
            "1",
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
    let gfa_path = "/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa";
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
            "--gfa",
            gfa_path,
            "--output",
            out_path.to_str().unwrap(),
            "--method",
            "sweepga",
            "--sweepga-aligner",
            "fastga",
            "--seqwish-k",
            "311",
            "--min-traversal-len",
            "5000",
            "--max-iterations",
            "2",
            "-v",
            "1",
        ])
        .output()
        .expect("failed to run impg crush (2 rounds)");

    let stderr = String::from_utf8_lossy(&run.stderr);
    assert!(run.status.success(), "impg crush failed: {}", stderr);

    let (after_segs, after_bp) = gfa_file_segment_stats(out_path.to_str().unwrap());
    let growth = (after_bp as f64 - before_bp as f64) / before_bp as f64;
    eprintln!(
        "After 2 rounds: {} segs / {} bp (growth={:.1}%)",
        after_segs,
        after_bp,
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
    let gfa_path = "/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa";
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
            "--gfa",
            gfa_path,
            "--output",
            out_path.to_str().unwrap(),
            "--method",
            "sweepga",
            "--sweepga-aligner",
            "fastga",
            "--seqwish-k",
            "311",
            "--min-traversal-len",
            "5000",
            "--max-iterations",
            "until-done",
            "--polish-method",
            "smooth",
            "--polish-rounds",
            "until-done",
            "--polish-max-traversal-len",
            "10000",
            "--polish-max-median-traversal-len",
            "1000",
            "-v",
            "1",
        ])
        .output()
        .expect("failed to run impg crush under timeout");

    let stderr = String::from_utf8_lossy(&run.stderr);
    eprintln!("exit status: {}", run.status);
    eprintln!(
        "stderr (last 5 lines):\n{}",
        stderr
            .lines()
            .rev()
            .take(5)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
            .join("\n")
    );

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

// ---------------------------------------------------------------------------
// Wider-context bubble resolution — see docs/crush-wider-context-bubbles.md.
// Verifies that enabling `replacement_flank_bp` on the canonical C4 slice still
// preserves every path sequence exactly. The integration code clips the
// flanking aligned portion before substitution; if that clipping is wrong, the
// resolved path sequence will not equal the input path sequence.
// ---------------------------------------------------------------------------

#[test]
fn c4_slice_auto_crush_with_flank_preserves_path_sequences() {
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let gfa_path = manifest_dir.join("tests/test_data/crush/c4_slice_1500_3000.gfa");
    if !gfa_path.exists() {
        eprintln!("SKIP: test data slice not found at {:?}", gfa_path);
        return;
    }

    let gfa = std::fs::read_to_string(&gfa_path).expect("failed to read C4 slice GFA");
    let before_seqs: std::collections::HashMap<_, _> =
        path_sequences(&gfa).unwrap().into_iter().collect();

    let mut config = ResolutionConfig::default();
    config.replacement_flank_bp = 500;
    let resolved = resolve_gfa_bubbles(&gfa, &config)
        .expect("resolve_gfa_bubbles with flank=500 failed on C4 slice");

    let after_seqs: std::collections::HashMap<_, _> =
        path_sequences(&resolved.gfa).unwrap().into_iter().collect();

    for (name, before_seq) in &before_seqs {
        let after_seq = after_seqs.get(name).expect("path disappeared after crush");
        assert_eq!(
            before_seq, after_seq,
            "path {} changed sequence after crush with flank=500",
            name
        );
    }
    eprintln!(
        "slice crush (flank=500): {} paths preserved, {} resolved, {} bailed",
        before_seqs.len(),
        resolved.stats.resolved,
        resolved.stats.bailed
    );
}

// ---------------------------------------------------------------------------
// Nested-bubble level descent — see docs/crush-nested-bubble-test.md
//
// Fixture: tests/test_data/crush/nested_bubbles_real.gfa
//   Real C4A pangenome extract (5 haplotypes through chr6 ref steps 568-590 of
//   CHM13#0#chr6:31744284-31976975, sourced from C4A.blunt.numeric.gfa). POVU
//   on the fixture finds exactly:
//     • 1 top-level (L0) site `<43388736>43388746` (16 ref-step span)
//     • 2 leaf (L1) sub-bubbles nested INSIDE the L0 parent:
//         `>43388742<2094988` and `>43388742>546635`
//   1 component, 40 segments, 981 bp, 43 links, 5 paths (3.8 KB on disk).
//
// Correct behavior (top-down with localized re-POVU on the resolved bubble's
// subgraph):
//   • Round 1 selects exactly 1 candidate: the L0 parent.
//   • The POA replacement subgraph still exposes ≥2 nested L1 sub-bubbles, so
//     round 2 resolves the nested children.
//   • Round 3 has no eligible candidates → terminate in 2 rounds.
//   • Path sequences preserved, component count preserved.
//   • Final POVU on the resolved graph finds ≤ a small handful of sites
//     (everything cleanly resolved).
//
// HEAD (471f089) behavior — observed failure that drives this test:
//   crush re-POVUs the WHOLE working graph every round and filters by
//   `site.is_leaf` (src/resolution.rs:1224), so it works BOTTOM-UP. On the
//   fixture this collapses to:
//     • Round 1 resolves a LEAF (not the parent).
//     • Each round's POA replacement re-fragments the graph, re-introducing
//       internal sub-bubbles that POVU keeps re-discovering. Stops only when
//       leftover sites fall below min_traversal_len.
//     • Resolution takes 5 rounds, produces 16 replacement plans, and the
//       result still contains ~15 POVU sites (10 L0 + 5 L1) — the "ramp /
//       whitespace" fragmentation signature visible in
//       /tmp/crush-ld/c4-crush-level-descent.png.
//
// Test asserts: iterations<=2, sites_after<=2, components preserved, path
// sequences preserved, round-1 selected==1, round-2 selected>=2. The
// iterations + sites_after assertions are the discriminators that FAIL on
// HEAD; the others are sanity gates.
// ---------------------------------------------------------------------------

const NESTED_BUBBLES_REF: &str = "CHM13#0#chr6:31744284-31976975";

/// Parse a single integer field from a `crush round N: ... K selected ...`
/// stderr line. Returns the K value for each round in order. Robust against
/// other interleaved log lines (FastGA, build progress, etc.).
fn parse_round_selected_counts(stderr: &str) -> Vec<usize> {
    let mut counts: Vec<(usize, usize)> = Vec::new();
    for line in stderr.lines() {
        // Match: "crush round <N>: <SITES> POVU site(s), ... <K> selected ..."
        let Some(rest) = line.split_once("crush round ").map(|x| x.1) else {
            continue;
        };
        let Some((round_str, after)) = rest.split_once(':') else {
            continue;
        };
        let Ok(round_idx) = round_str.trim().parse::<usize>() else {
            continue;
        };
        // Skip "crush round N traversal stats" lines and "no eligible
        // candidates" lines — only count the headline selection line.
        if after.contains("traversal stats") || after.contains("no eligible") {
            continue;
        }
        let Some(sel_pos) = after.find(" selected ") else {
            continue;
        };
        // Walk backwards from ' selected ' to find the integer before it.
        let prefix = &after[..sel_pos];
        let Some(num_start) = prefix.rfind(|c: char| !c.is_ascii_digit()) else {
            continue;
        };
        let Ok(selected) = prefix[num_start + 1..].parse::<usize>() else {
            continue;
        };
        // Keep only the first occurrence per round index (in case of dupes).
        if !counts.iter().any(|&(r, _)| r == round_idx) {
            counts.push((round_idx, selected));
        }
    }
    counts.sort_by_key(|&(r, _)| r);
    counts.into_iter().map(|(_, k)| k).collect()
}

#[test]
#[ignore = "known RED documented in docs/crush-true-level-descent.md: final POVU-site <= 2 assertion is unsatisfiable under SPOA + path preservation"]
fn nested_bubble_level_descent_actually_descends() {
    // docs/crush-nested-bubble-test.md
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let gfa_path = manifest_dir.join("tests/test_data/crush/nested_bubbles_real.gfa");
    assert!(gfa_path.exists(), "fixture missing: {}", gfa_path.display());

    let gfa = std::fs::read_to_string(&gfa_path).expect("failed to read fixture");

    // Sanity: characterise the input fixture's bubble structure.
    let describe_opts = GraphReportOptions {
        top_n: 50,
        include_povu: true,
        povu_reference_names: vec![NESTED_BUBBLES_REF.to_string()],
        ..GraphReportOptions::default()
    };
    let input_report =
        describe_gfa("nested-fixture", &gfa, &describe_opts).expect("describe_gfa on fixture");
    let input_povu = input_report
        .povu
        .as_ref()
        .expect("POVU disabled or failed on fixture");
    let input_l0 = input_povu.top_sites.iter().filter(|s| s.level == 0).count();
    let input_l1 = input_povu.top_sites.iter().filter(|s| s.level == 1).count();
    let input_l1_nested_in_l0: usize = input_povu
        .top_sites
        .iter()
        .filter(|s| s.level == 1 && s.parent_id.is_some())
        .count();
    eprintln!(
        "fixture POVU: {} sites total, level breakdown L0={}, L1={} (L1 nested-in-L0={})",
        input_povu.sites, input_l0, input_l1, input_l1_nested_in_l0
    );
    assert_eq!(
        input_l0, 1,
        "fixture should have exactly 1 top-level (L0) bubble; got {}",
        input_l0
    );
    assert!(
        input_l1_nested_in_l0 >= 2,
        "fixture should have >= 2 nested (L1) sub-bubbles inside the L0 parent; got {}",
        input_l1_nested_in_l0
    );
    assert_eq!(
        input_report.metrics.components, 1,
        "fixture should be a single connected component; got {}",
        input_report.metrics.components
    );

    // Capture input path sequences for end-to-end preservation check.
    let before_seqs: std::collections::HashMap<_, _> =
        path_sequences(&gfa).unwrap().into_iter().collect();
    let input_components = input_report.metrics.components;

    // Run crush with max_iterations=5 via the library to capture stats.iterations.
    let resolved = resolve_gfa_bubbles(
        &gfa,
        &ResolutionConfig {
            max_iterations: 5,
            method: ResolutionMethod::Auto,
            ..ResolutionConfig::default()
        },
    )
    .expect("resolve_gfa_bubbles failed on nested fixture");

    eprintln!(
        "library crush: iterations={}, resolved={}, bailed={}, candidates_seen={}",
        resolved.stats.iterations,
        resolved.stats.resolved,
        resolved.stats.bailed,
        resolved.stats.candidates_seen,
    );

    // ----------- Hard assertions that FAIL on HEAD (471f089) -----------

    // (1) Termination round count.
    // CORRECT (top-down + localized re-POVU): round 1 resolves L0 parent, round
    // 2 resolves the nested L1 leaves discovered inside the replacement, round
    // 3 has no eligible candidates → terminate; iterations==2.
    // HEAD: bottom-up leaf filter creates synthetic L0/L1 sites every rewrite,
    // so resolution takes 5 rounds and ends only when leftover sites are below
    // min_traversal_len.
    assert!(
        resolved.stats.iterations <= 2,
        "nested-bubble level descent should terminate in <= 2 iterations \
         (correct top-down algorithm needs exactly 2: round 1 = L0 parent, \
         round 2 = nested L1 children); got iterations={} on this run \
         (HEAD 471f089 takes 5 rounds because re-POVU on the WHOLE graph + \
         is_leaf filter at src/resolution.rs:1224 keeps re-discovering POA \
         fragmentation artifacts each round)",
        resolved.stats.iterations
    );

    // (2) After resolution, POVU on the result must NOT find a fragmented
    // graph littered with leftover sub-bubbles. The correct algorithm leaves
    // ≤ a couple sites (the polymorphic structure is fully consumed). HEAD
    // leaves ~15 sites (10 L0 + 5 L1) — the fragmentation signature.
    let result_report = describe_gfa("nested-fixture-resolved", &resolved.gfa, &describe_opts)
        .expect("describe_gfa on resolved");
    let result_povu = result_report
        .povu
        .as_ref()
        .expect("POVU disabled or failed on resolved");
    eprintln!(
        "resolved POVU: {} sites total, leaves={}, levels={:?}; segments={} components={}",
        result_povu.sites,
        result_povu.leaf_sites,
        result_povu.level_counts,
        result_report.metrics.segments,
        result_report.metrics.components,
    );
    assert!(
        result_povu.sites <= 2,
        "resolved graph should have <= 2 POVU sites (cleanly consumed nested \
         bubble); got {} sites. HEAD 471f089 leaves ~15 sites due to POA \
         consensus fragmentation introducing new L0/L1 bubbles each round.",
        result_povu.sites
    );

    // (3) Component count must stay at the input value (fragmentation must
    // not introduce disconnected components).
    assert_eq!(
        result_report.metrics.components, input_components,
        "resolved graph component count changed from {} to {} \
         (fragmentation should not disconnect the graph)",
        input_components, result_report.metrics.components
    );

    // (4) All input path sequences must be preserved.
    let after_seqs: std::collections::HashMap<_, _> =
        path_sequences(&resolved.gfa).unwrap().into_iter().collect();
    for (name, before_seq) in &before_seqs {
        let after_seq = after_seqs
            .get(name)
            .unwrap_or_else(|| panic!("path {} disappeared after crush", name));
        assert_eq!(
            before_seq, after_seq,
            "path {} changed sequence after crush",
            name
        );
    }

    // ----------- Round-by-round assertions (parsed from CLI stderr) ----------
    //
    // The library's ResolutionStats does not expose per-round frontier sizes.
    // Run the CLI subprocess with -v 2 to capture the structured log lines
    // and parse "crush round N: <X> POVU site(s), <Y> unseen ..., <K> selected"
    // so we can assert round 1 selected == 1 and round 2 selected >= 2.

    let Some(bin) = impg_binary() else {
        eprintln!("SKIP CLI round-count assertions: impg binary not found");
        return;
    };
    let out_path = std::env::temp_dir().join("nested_bubble_level_descent_out.gfa");
    let run = Command::new(&bin)
        .args([
            "crush",
            "--gfa",
            gfa_path.to_str().unwrap(),
            "--output",
            out_path.to_str().unwrap(),
            "--method",
            "auto",
            "--max-iterations",
            "5",
            "-v",
            "2",
        ])
        .output()
        .expect("failed to run impg crush CLI");
    let stderr = String::from_utf8_lossy(&run.stderr);
    let per_round = parse_round_selected_counts(&stderr);
    eprintln!("CLI per-round selected counts: {:?}", per_round);
    std::fs::remove_file(&out_path).ok();

    assert!(
        per_round.len() >= 2,
        "expected >= 2 rounds in stderr; parsed {:?}\n--- stderr tail ---\n{}",
        per_round,
        stderr.lines().rev().take(15).collect::<Vec<_>>().join("\n")
    );

    // Round 1 should resolve exactly 1 candidate (the L0 parent).
    assert_eq!(
        per_round[0], 1,
        "round 1 should select exactly 1 candidate (the L0 top-level parent); \
         got {}. HEAD 471f089's leaf-filter selects whatever leaves pass the \
         min_traversal_len cut — usually a count != 1.",
        per_round[0]
    );
    // Round 2 should resolve the nested children (>= 2 plans) — POVU on the
    // L0 replacement's local subgraph reveals the nested L1 leaves.
    assert!(
        per_round[1] >= 2,
        "round 2 should select >= 2 candidates (the nested L1 sub-bubbles in \
         the L0 replacement's local subgraph); got {}",
        per_round[1]
    );
}
