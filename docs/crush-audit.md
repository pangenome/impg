# Crush behavioral audit on real C4 data

**Branch / HEAD:** `90ba74f` (worktree `wg/agent-49/crush-audit`).
**Reproducer command:** the `current_auto` block in `docs/c4-crush-handoff.md`
(GRCh38#0#chr6:31891045-32123783, `-d 50k`, `method=sweepga,aligner=fastga`,
`seqwish-k=311`, `polish-rounds=until-done`, polish defaults to `smooth`).
**Reproducer artifacts:** `data/c4_handoff_repro_20260524T142743Z/current_auto.nosort.{stderr,stdout}`
(this worktree). Inputs are the real syng+agc at
`/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.{syng,agc}`, 32 threads.
**Run wall budget:** the canonical command was capped at 1800 s (30 min) with
`timeout --signal=SIGTERM --kill-after=60s 1800s`. The process exited with code 124
at `14:58:01Z` (the kernel SIGTERM from `timeout`); no output GFA was emitted.
**Known-good baseline for comparison:**
`/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.parent5k.sweepga_allvsall_fastga.k311.done.nosort.{gfa,stderr}`,
generated 2026-05-23T21:44Z from the same syng query with the same canonical command,
predating commit `0af1a4c` ("Accept validated crush replacements"). Wall 4:39, max RSS 55 GB.

This audit is observation only. No code change is committed.

## Round-by-round trajectory observed on current HEAD

Each row is taken verbatim from the `crush round N: resolved ...` summary line in
`current_auto.nosort.stderr` (numbers in the source log are unspaced; spaces here are for
readability only).

| round | wall  | accepted | score (before → after)        | segments (Δ) | seg-bp (Δ)        | links (Δ) | path-steps (Δ)    | ws-p99 (Δ)          | ws-max (Δ)          | ws-long≥10k (Δ) |
|-------|-------|----------|--------------------------------|--------------|-------------------|-----------|--------------------|----------------------|----------------------|------------------|
| 0 (start) | -    | -        | 217 461 559                    | 18 048       | 389 354           | 20 933    | 4 591 855          | 177 799              | 388 845              | 208 494          |
| 1     | 79.05s| 8/8      | 217 461 559 → **357 006 340**  | 18 048 → 20 919 (+2 871) | 389 354 → **696 283** (+78.8%) | 20 933 → 25 323 (+4 390) | 4 591 855 → 4 333 110 (-5.6%) | 177 799 → **316 426** (+77.9%) | 388 845 → **691 267** (+77.8%) | 208 494 → 145 844 (-30.1%) |
| 2     | 793.14s| 8/8     | 357 006 340 → 208 208 844      | 20 919 → 21 701 (+782)   | 696 283 → 740 381 (+6.3%)      | 25 323 → 26 353 (+1 030) | 4 333 110 → 4 182 971 (-3.5%) | 316 426 → 172 155 (-45.6%)     | 691 267 → 734 444 (+6.2%)      | 145 844 → 139 029 (-4.7%) |
| 3     | n/a (timeout) | 7/8 then SIGTERM | (round 3 never completed) | | | | | | | |

For comparison the known-good baseline reaches its final state in *one accepted
replacement*, with `score` improving by -30.8%, `seg-bp` growing by only +8.4% (+32 653
bp), and the next round (round 2) rejected outright on a +6.18% score growth.

Net effect after 25:50 of accepted work on HEAD: `score` 217 461 559 → 208 208 844
(-4.3%, much worse than baseline's 1-round -30.8%), `seg-bp` 389 354 → 740 381
(+90.2%), `segments` 18 048 → 21 701 (+20.2%), `ws-max` 388 845 → 734 444 (+88.9%);
then the third round hung.

## Pipeline stages

For pinning failures the stages of the crush pipeline are:

| stage | code |
|-------|------|
| bubble identification (flubble detect, POVU)        | `find_candidate_frontier`, src/resolution.rs:458 |
| aligner selection by bubble size                    | `candidate_selection_method`, src/resolution.rs:838 |
| aligner invocation                                  | `build_replacement_with_method`, src/resolution.rs:1544 (sweepga at :1891, allwave at :1797) |
| replacement integration into working graph          | round loop, src/resolution.rs:454-632 (no quality gate); `apply_replacement_frontier`, :1452 |
| path rewriting through replaced bubble              | `apply_replacement_frontier` step rewrite, src/resolution.rs:1480-1528 |
| graph emission/serialization                        | `render_rewritten_graph`, src/resolution.rs:2810 |

## Failures

### Failure 1 — the canonical command does not finish on real C4 within 30 minutes; final GFA is never written

(a) **Input slice.** The exact canonical command from `docs/c4-crush-handoff.md` §Baseline
(`current_auto` variant):

```bash
impg query -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O .../current_auto.nosort -v 1
```

(b) **Expected behavior.** The known-good baseline (predates commit `0af1a4c`) ran the
identical command in 4:39 wall, 55 GB RSS, and emitted a 47 MB sorted GFA that the
report (`...Ygs.report.md`) characterises as `18 089 segments / 21 018 links / 465 paths
/ 4 590 745 path steps / 422 007 segment bp`.

(c) **Actual behavior.** The same command on HEAD does not finish within a 30 min wall
budget. After 30 min the `timeout` wrapper returns exit code 124 and no output GFA is
written (`current_auto.nosort.stdout` is 0 bytes; the in-progress `nosort.gfa` file is
not created because the rounds never converge to a final write). Last log line in
`current_auto.nosort.stderr`:

```
[2026-05-24T14:52:30Z INFO  impg::resolution] crush round 3: replacement build progress 7/8 (accepted)
```

Round-by-round wall costs grew sharply:
- round 1: 79.05s (matches the baseline scale)
- round 2: 793.14s (10x worse) — the per-replacement build on the largest round-2
  candidate (`traversals=433, max-len=42 544, median-len=42 430, total-len=18 370 524`)
  dominated this round.
- round 3: 7/8 done in 519s, then ≥330s on the 8th before SIGTERM. The biggest
  round-3 candidate is `traversals=440, max-len=42 628, median-len=42 513,
  total-len=18 702 925`.

The proximate cause of the hang is that each accepted round grows the working graph in
the bubble (see Failures 2-4), so the next POVU pass selects bigger candidates, and
sweepga's `build_sweepga_seqwish_replacement` followed by the smooth polish on a
~18 MB candidate dominates wall time.

(d) **Pipeline stage where divergence first appears.** Replacement integration into
working graph — the round loop in `src/resolution.rs:454-632`. The loop has no
acceptance gate other than per-path sequence equality (see Failure 2) and no
runtime/work budget per round.

(e) **Test that should have caught it and why it didn't.**

- `src/resolution.rs:3483` `first_crush_round_accepts_path_step_growth_without_visual_tail_regression`
  is the only test that exercises the multi-round loop, but on a 100-bp synthetic
  two-path GFA with `StarBiwfa`. It cannot reach the sweepga + smooth path that
  dominates wall time on real data, and it has no wall-clock or growth assertion.
- `tests/test_syng_integration.rs:207` `test_crush_cli_resolves_blunt_gfa` exercises
  `impg crush` end-to-end but on a 7-bp insertion bubble. It also has no wall-clock
  or growth assertion.
- `tests/test_pipeline_integration.rs:59` `test_full_pipeline` is `#[ignore]`-gated,
  runs `index → partition → graph → lace` on 7-strain yeast chrV, and never goes
  through the `crush` stage at all.

No test in the repo runs the canonical C4 command (or anything close to its scale)
under a wall-clock budget, so the regression past the 4:39 baseline is not caught
by CI.

### Failure 2 — there is no round-level quality gate; rounds that explode the graph are accepted unconditionally

(a) **Input slice.** Same canonical command. The relevant per-round summary is round 1
on HEAD vs the same round 1 on the baseline.

(b) **Expected behavior.** When a round's collective replacement would grow the
visual-tail score substantially, reject the round and keep the pre-round graph. The
known-good baseline shows this exact behavior at `2026-05-23T21:44:20Z`:

```
crush round 2: rejecting 3 replacement(s): visual-tail score grew by 0.0618, above allowed 0.0200;
  before score=150573617, segments=18089, segment-bp=422007, ...
  after  score=159879124, segments=18538, segment-bp=708222, ...
```

The reference is the implementation history: commit `259689b` ("Guard crush quality
before accepting replacements", 2026-05-23) introduced a per-replacement quality guard
with test `first_crush_round_rejects_quality_regression`. Commit `b163321` ("Score crush
rounds by visual tail quality", 2026-05-23) added `round_quality_decision()` with
`DEFAULT_MAX_ROUND_SCORE_GROWTH = 0.02`. Both were present in the baseline binary.

(c) **Actual behavior.** Round 1 on HEAD (`current_auto.nosort.stderr:14:33:56Z`)
accepts all 8/8 replacements unconditionally. Quantitative comparison:

|                        | baseline rd1 (accepted)       | baseline rd2 (rejected)     | HEAD rd1 (accepted)         |
|------------------------|-------------------------------:|----------------------------:|------------------------------:|
| replacements accepted  | 1/8                            | 0/3                         | **8/8**                       |
| score growth           | -30.8% (improve)               | +6.18% (rejected)           | **+64.2% (accepted)**         |
| segments added         | +41                            | +449                        | **+2 871**                    |
| segment bp added       | +32 653                        | +286 215                    | **+306 929**                  |
| links added            | +85                            | +901                        | **+4 390**                    |
| ws-total change        | -40%                           | +9.3%                       | **+15.2%**                    |
| ws-p99 change          | -31%                           | +3.5%                       | **+77.9%**                    |
| ws-max change          | +7%                            | +65.9%                      | **+77.8%**                    |

Every regression metric the baseline guard would have triggered on is exceeded by the
HEAD round-1 result, but HEAD has no guard.

(d) **Pipeline stage where divergence first appears.** Replacement integration —
src/resolution.rs:454-632. Specifically src/resolution.rs:611-616:

```rust
let resolved_count = plans.len();
let rewrite_start = Instant::now();
let next_graph = apply_replacement_frontier(&graph, &plans)?;
let rewrite_elapsed = rewrite_start.elapsed();
let after_quality = graph_quality(&next_graph);
graph = next_graph;                       // unconditional accept
```

`before_quality` and `after_quality` are still computed and logged (the `summary()` at
line 977), but the code that compared them and rolled the graph back was removed in
commit `0af1a4c` ("Accept validated crush replacements", 2026-05-23T23:47:54Z, +343
-376 lines in src/resolution.rs). That commit removed:

- `pub max_round_score_growth: f64` (ResolutionConfig field) and
  `pub const DEFAULT_MAX_ROUND_SCORE_GROWTH: f64 = 0.02;`,
- `fn round_quality_decision` and `enum RoundQualityDecision { Reject { reason }, ... }`,
- `fn replacement_method_needs_quality_guard` and the per-replacement quality guard for
  pairwise (sweepga/allwave) methods,
- the round-level rollback that produced
  `crush round N: rejecting M replacement(s): visual-tail score grew by ...`.

The orphan comment that survives at src/resolution.rs:1030-1034 now reads:

```rust
// This is a visual-tail score, not a raw complexity count. Crushing can
// legitimately increase path steps while improving the graph the user sees:
// fewer long empty jumps, tighter high-depth node reuse, and shorter
// segment sequence. This is logged for diagnostics only; acceptance is
// based on exact replacement path validation.
```

The only surviving acceptance criterion is `path_sequences_equal`
(src/resolution.rs:1536, 2796), i.e. "every input path still spells the same
sequence". A replacement plan can quintuple segment count, double segment bp, and
double `ws-p99` and the round still passes.

(e) **Test that should have caught it and why it didn't.**

- `src/resolution.rs:3482` `first_crush_round_accepts_path_step_growth_without_visual_tail_regression`
  is the renamed-and-inverted descendant of the original
  `first_crush_round_rejects_quality_regression`. History (via `git log -G`):
  added by `259689b`; renamed and its assertion flipped by `b163321`; survived
  `0af1a4c`. It now builds a 100-bp two-path synthetic GFA, runs `StarBiwfa`, and
  asserts only `resolved.stats.resolved == 1`. It never runs `sweepga`, never tests
  the multi-round growth path, and never asserts an upper bound on any quality
  metric. Removing the guard does not change the result of this test, which is
  why CI stayed green through commit `0af1a4c`.
- The tests that *would* have caught the gate removal — `round_quality_rejects_large_complexity_growth`,
  `round_quality_accepts_visual_tail_improvement_despite_score_growth`,
  `direct_poa_replacements_bypass_local_quality_gate` — were deleted by the same
  commit (`0af1a4c`).
- `tests/test_syng_integration.rs:207` `test_crush_cli_resolves_blunt_gfa` is the
  only end-to-end CLI test for `impg crush`; it asserts only that the binary exits 0
  and stderr contains `crush: 1 resolved`, on a 7-bp insertion bubble.

No surviving test exercises a multi-replacement round on a real graph and asserts an
upper bound on segment count, segment bp, links, or ws-* metrics.

### Failure 3 — replacement integration multiplies segment bp because cross-replacement segments are never deduplicated

(a) **Input slice.** Same canonical command. Round 1 accepts 8 sweepga/seqwish
replacements in parallel covering 8 different POVU sites. The C4 locus has tandem-
paralog structure (the report on the baseline graph at
`C4A.parent5k.sweepga_allvsall_fastga.k311.done.Ygs.report.md` flags
`duplicate_sequence_frac` and `local_repeat_contexts` as warnings), so independent
sites have many short locally-identical sequences.

(b) **Expected behavior.** When two replacements emit segments with the same DNA
sequence the combined working graph should share them or at minimum the total segment
bp should not grow super-linearly in the number of accepted replacements. The
known-good baseline shows the desired scale: +41 segments and +32 653 segment bp for
the 1 accepted replacement.

(c) **Actual behavior.** Round 1 of HEAD: +2 871 segments and +306 929 segment bp for
8 accepted replacements. That is ~70x more segments and ~9x more added segment bp
*per accepted replacement* than the baseline saw. The input to round 1 had 389 354 bp
of unique segment sequence across 18 048 segments; the output has 696 283 bp across
20 919 segments, i.e. 78.8% more segment bp.

(d) **Pipeline stage where divergence first appears.** Graph emission —
`render_rewritten_graph`, src/resolution.rs:2810-2900. The function iterates over
`out_paths`, collects unique `OutNode` values (`OutNode::Original(idx)` for
original-graph indices and `OutNode::Replacement(plan_idx, node_idx)` for
per-replacement nodes), and for every `Replacement(plan_idx, node_idx)` node it
calls `next_unused_segment_id(&mut used_ids, &mut next_id)` (line 2851) and emits a
fresh `S` line with the replacement's sequence. Two replacements that produced
byte-identical segment sequences end up as two separate `S` lines with two different
numeric IDs; there is no `seq -> id` deduplication in this function. The only sequence
collapsing that exists is *within* one replacement (seqwish does that as part of
building the replacement graph) and across replacements via `gfaffix` only if the
caller runs gfaffix afterward — which the `nosort` output stage explicitly skips.

(e) **Test that should have caught it and why it didn't.**

- `src/resolution.rs:3482` `first_crush_round_accepts_path_step_growth_without_visual_tail_regression`
  uses one replacement of length ≤100 bp on a 2-path GFA; even if every emitted
  segment were duplicated, bp growth would be a handful of bp. No segment-count
  invariant is asserted.
- `src/resolution.rs:3252` `auto_routes_small_bubbles_to_spoa_and_larger_bubbles_to_allwave`
  checks the routing decision but not the post-rewrite graph shape.
- `tests/test_syng_integration.rs:207` `test_crush_cli_resolves_blunt_gfa` runs
  end-to-end but with one bubble; the segment count assertion only requires `\nS\t`
  to appear in the output.

No test exercises a many-replacement round on a real graph and asserts that segment
bp grows by O(replacement work) rather than O(replacement work × replacement count),
and no test asserts the absence of duplicate-sequence segments emitted by
`render_rewritten_graph`.

### Failure 4 — round-over-round candidate inflation: accepting round 1 makes round 2 see the same site as a much larger bubble; sweepga + smooth then dominates wall time

(a) **Input slice.** Same canonical command. The bubble in question is the
highest-traversal candidate in the first POVU pass. Round-1 logs:

```
crush round 1: 2441 POVU site(s), 72 unseen polymorphic candidate(s), 8 selected, 0 selection-guarded
  candidate 1: traversals=312  max=31 478  median=25 155  total=5 132 443   root-span=110
  candidate 4: traversals=437  max=42 362  median=157     total=2 893 782   root-span=157
```

After round 1 accepts all 8 replacements, POVU re-runs on the rewritten working graph:

```
crush round 2: 2340 POVU site(s), 61 unseen polymorphic candidate(s), 8 selected, 0 selection-guarded
  candidate 1: traversals=433  max=42 544  median=42 430  total=18 370 524  root-span=42 394
  candidate 8: traversals=312  max=31 536  median=25 281  total=7 966 585   root-span=31 376
```

And after round 2:

```
crush round 3: 2648 POVU site(s), 49 unseen polymorphic candidate(s), 8 selected, 0 selection-guarded
  candidate 1: traversals=440  max=42 628  median=42 513  total=18 702 925  root-span=42 478
  candidate 2: traversals=384  max=31 691  median=25 439  total=9 116 260   root-span=327
```

The relevant change is the jump in `median-len` (round-1 candidates with `median-len=157`
because most traversals were short come back in round 2 with `median-len=42 430`) and
`total-len` (round-1 sum across all 8 selected candidates is 14 585 526 bp of traversal
sequence; round-2 selected candidates already sum to ~35 MB of traversal sequence; the
round-3 selected candidates sum to ~36 MB).

(b) **Expected behavior.** A useful round should monotonically reduce the work the
next round sees. The baseline confirms this: baseline round 1 accepted 1/8, the
round-1 working graph had `seg-bp=422 007 / ws-total=14.84 Gbp`, and POVU on that
graph found 3 candidates in round 2 which the guard then rejected. The
next-round candidate budget did not blow up. In particular: in the baseline, after
round 1 the total selection-eligible traversal sequence is bounded by `seg-bp`
(422 kb), so even the worst-case selected candidate is bounded. On HEAD that
invariant is broken because `seg-bp` itself grew (Failure 3) and the bubble
boundaries widened.

(c) **Actual behavior.** HEAD round 1 selects 8 candidates with combined traversal
sequence ~15 MB; after acceptance, round 2 selects 8 candidates totaling ~35 MB with
the largest single candidate at 18.4 MB. Smooth+SPOA on 18 MB of traversal sequence
over 433 paths takes 12+ minutes wall on this hardware, and the same shape candidate
comes back in round 3 at 18.7 MB. The round-3 8th replacement build had not finished
after 5+ minutes when the 30-minute SIGTERM hit (Failure 1).

(d) **Pipeline stage where divergence first appears.** Replacement integration
causes bubble boundaries to widen: when an accepted replacement re-aligns short
traversals against the longest traversal (sweepga all-vs-all, `seqwish-k=311`), the
resulting compacted replacement collapses sequence-identical anchors that previously
fenced the bubble. The next POVU pass on the rewritten graph sees a wider flubble at
the same path region. Concretely the proximate divergence is the lack of a
round-acceptance gate (Failure 2) — without that gate the loop has no
monotone-progress condition. The structural cause is at the boundary between
`apply_replacement_frontier` (src/resolution.rs:1452) and the next
`find_candidate_frontier` (src/resolution.rs:458) iteration.

(e) **Test that should have caught it and why it didn't.** No test in the repo runs
>1 round on a real or even synthetic non-trivial graph and asserts that round N+1's
selected budget is smaller than round N's. The candidates one might expect to
include such an assertion (`auto_resolves_indel_ladder_without_changing_path_sequences`
at src/resolution.rs:3578, `resolves_maximal_eligible_nested_site_first` at :3673,
`resolves_independent_sites_in_one_round` at :3722) all assert single-round
resolution and never re-run POVU after acceptance. The closest historical guard —
the deleted `round_quality_rejects_large_complexity_growth` test — only checked
that growth above a threshold caused rejection, not that subsequent rounds saw
non-inflated candidates.

## Notes on what was *not* a failure

- Aligner correctness: the aligner-level pieces (sweepga / FastGA / SPOA / POASTA)
  produce the per-replacement output expected by the rest of the pipeline; their
  per-replacement results pass the `validate_replacement_paths` byte-equality check
  at src/resolution.rs:2436.
- Flubble detection: POVU site counts are consistent run-over-run (2441 → 2340 →
  2648) and the selected candidates' `root-span` and `traversals` are internally
  consistent with the working graph at that point in time. POVU is correctly
  reporting the bubble shape it sees; the bubble shape just keeps changing because
  the working graph keeps changing (Failure 4).
- Path-sequence preservation: every accepted round preserves every input path
  sequence exactly (the `path_sequences_equal` check at src/resolution.rs:1536
  enforces this on every round). The graph grows and white-space widens, but no
  path is corrupted.
