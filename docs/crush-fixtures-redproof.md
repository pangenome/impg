# crush-fixtures red-proof

Captured test output demonstrating each `#[ignore]` test FAILS on HEAD `90ba74f`.

Run environment: `/home/erikg/impg/.wg-worktrees/agent-52`, 2026-05-24.

---

## Test 1 — Failure 2: No round-level quality gate

**Test:** `c4_round1_sweep_quality_gate_rejects_score_growth`

**Assertion:** segment-bp growth after round 1 ≤ 10%

**Observed on HEAD:** 18.7% growth (before=389,316 bp → after=461,998 bp)

**Root cause in code:** commit `0af1a4c` removed `DEFAULT_MAX_ROUND_SCORE_GROWTH=0.02`,
`fn round_quality_decision`, and the round-level rollback from `src/resolution.rs`.
The round loop at line 616 now unconditionally accepts `graph = next_graph`.

**Run output (key lines):**
```
Input: 18048 segments / 389316 bp
[INFO  impg::resolution] crush round 1: resolved 2/8 replacement(s) in 35.45s;
  rewrite+validate 14.42s; total 75.44s;
  quality score=217455065, segments=18048, segment-bp=389316, ...
  -> score=151542250, segments=18206, segment-bp=461998, ...
After round 1: 18206 segments / 461998 bp (Δsegs=+158, Δbp=+72682, bp_growth=18.7%)

thread 'c4_round1_sweep_quality_gate_rejects_score_growth' panicked at tests/test_crush_integration.rs:205:5:
Failure 2: segment-bp grew 18.7% in round 1 (before=389316, after=461998) —
the quality gate (DEFAULT_MAX_ROUND_SCORE_GROWTH=0.02) was removed in 0af1a4c;
this round should have been rejected

test result: FAILED. 0 passed; 1 failed; 0 ignored; finished in 86.85s
```

---

## Test 2 — Failure 3: Duplicate segment sequences in render output

**Test:** `c4_round1_render_emits_no_duplicate_segment_sequences`

**Assertion:** zero S-lines sharing a sequence after round 1

**Observed on HEAD:** 5,474 duplicate sequences (before=18,048 segs → after=18,206 segs)

**Root cause in code:** `render_rewritten_graph` at `src/resolution.rs:2851` calls
`next_unused_segment_id()` to assign fresh IDs to every replacement segment without
deduplicating sequences across parallel replacements. Tandem-paralog bubbles in C4
produce byte-identical sequences from independent replacements.

**Run output (key lines):**
```
Before: 18048 segs | After: 18206 segs | Duplicate seqs: 5474

thread 'c4_round1_render_emits_no_duplicate_segment_sequences' panicked at tests/test_crush_integration.rs:286:5:
assertion `left == right` failed: Failure 3: 5474 segment(s) share a sequence with
another segment after round 1; render_rewritten_graph does not deduplicate across
parallel replacements (see src/resolution.rs:2851)
  left: 5474
 right: 0

test result: FAILED. 0 passed; 1 failed; 0 ignored; finished in 80.93s
```

---

## Test 3 — Failure 4: Round-over-round candidate inflation

**Test:** `c4_round2_segment_bp_does_not_exceed_round1`

**Assertion:** segment-bp growth after 2 rounds ≤ 20%

**Observed on HEAD:** 33.6% growth (before=389,316 bp → after=520,224 bp across 2 rounds)

**Root cause in code:** Round 1 is accepted despite 18.7% growth (Failure 2). Round 2
then sees a larger, more-complex graph and introduces additional inflation. The cascade
continues indefinitely with no quality gate to halt it.

**Run output (key lines):**
```
After 2 rounds: 18360 segs / 520224 bp (growth=33.6%)

thread 'c4_round2_segment_bp_does_not_exceed_round1' panicked at tests/test_crush_integration.rs:362:5:
Failure 4: segment-bp grew 33.6% across 2 rounds (before=389316, after=520224) —
candidate inflation caused by accepting quality-degrading rounds; the quality gate
removed in 0af1a4c would have stopped this

test result: FAILED. 0 passed; 1 failed; 0 ignored; finished in 163.80s
```

---

## Test 4 — Failure 1: Canonical command does not finish within budget

**Test:** `c4_canonical_command_completes_within_budget`

**Assertion:** `impg crush ... --max-iterations until-done` completes within 360 seconds

**Observed on HEAD:** exit code 124 (killed by `timeout 360s`); process was still in
round 4 when killed (progress log showed `6/6 accepted` for round 4 just before kill)

**Root cause:** Failures 2 and 4 cascade — each accepted round inflates the graph,
causing subsequent rounds to see larger bubbles with more traversals. The known-good
baseline (pre-`0af1a4c`) completed the same command in 4:39 (279s); on HEAD the
command never finishes within 30 minutes.

**Run output (key lines):**
```
[INFO  impg::resolution] crush round 4: replacement build progress 5/6 (accepted)
[INFO  impg::smooth] [smooth] Smoothed 20 blocks (5 passthrough)
[INFO  impg::smooth] [smooth] Lacing complete
[INFO  impg::resolution] crush round 4: replacement build progress 6/6 (accepted)
exit status: exit status: 124

thread 'c4_canonical_command_completes_within_budget' panicked at tests/test_crush_integration.rs:434:5:
Failure 1: impg crush did not complete within 360s on the C4 canonical command
(exit=Some(124)). On HEAD (90ba74f) round 3 hangs because rounds 1 and 2 were accepted
despite segment-bp growth (Failure 2), inflating subsequent candidates (Failure 4).
The known-good baseline (pre-0af1a4c) completed in 4:39.

test result: FAILED. 0 passed; 1 failed; 0 ignored; finished in 360.24s
```

---

## Companion test — passes on HEAD

**Test:** `c4_slice_auto_crush_preserves_path_sequences` (no `#[ignore]`)

**Input:** `tests/test_data/crush/c4_slice_1500_3000.gfa` (2942 segs / 64,348 bp, real C4A slice)

**Result:** PASS — 465 paths preserved, 147 bubbles resolved, 0 bailed, bp reduced
64,348 → 58,003 (–9.9%)

**Output:**
```
slice crush: 465 paths preserved, 147 resolved, 0 bailed
slice: segs 2942 → 2984, bp 64348 → 58003
test c4_slice_auto_crush_preserves_path_sequences ... ok
test result: ok. 1 passed; 0 failed; finished in 14.95s
```
