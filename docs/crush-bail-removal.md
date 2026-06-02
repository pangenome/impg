# crush — bail/rejection mechanism removal

**Task:** `crush-remove-bails`
**Date:** 2026-05-25
**Branch:** `wg/agent-83/crush-remove-bails`
**Author:** agent-83

## What was asked

Remove ALL bail / rejection mechanisms in `impg crush`. Per user direction:
the algorithm should run POA-style local alignment on every bubble at every
level until no bubbles remain. No quality-based bailing. `path_sequences_equal`
remains the only correctness invariant (paths must spell their input sequences
exactly).

## Investigation

### Stale claim in the task description

The task notes that the prior `crush-add-gfaffix-2` run reported:

```
crush round 1: rejecting 33/33 replacement(s); score grew 15.3% (threshold 10.0%)
```

I could not reproduce this message at the task's base commit (`b7b7e92`).
A `git log -S "score grew"` shows the round-quality reject logic was actually
deleted in commit `0af1a4c` ("Accept validated crush replacements"), which
predates `b7b7e92`. The "score grew" message in `docs/crush-gfaffix-run.md` is
a stale paste from older logs — the live binary at `b7b7e92` already prints
`0 bailed` on the slice.

Baseline confirmation (before any of this task's edits, using the existing
`/home/erikg/impg/target/release/impg` binary at commit `b7b7e92`):

```
impg crush -g tests/test_data/crush/c4_slice_1500_3000.gfa \
  --method auto --max-iterations 10
…
crush: 2027 resolved, 0 bailed, 2027 candidates seen across 10 rounds
```

So the *round-quality* gate from `0af1a4c` is genuinely gone. The
mechanisms removed by THIS task are the remaining dormant guards plus a
small number of cost caps in the pairwise-induction backends.

### Bail mechanisms found in current source (pre-removal)

| # | `src/resolution.rs` line(s) | Mechanism | Status |
|---|------|------|------|
| 1 | 1271–1275 | Selection-guard: `if !candidate.within_budget { frontier.bailed.push(candidate); continue; }` | Dormant under default config (`max_bubble_span=0`) but live code path |
| 2 | 491–494 | Round-loop counter that increments `stats.bailed` per selection-guarded candidate | Dormant counter |
| 3 | 496–504 | Early `break` in the round loop with log `"stopping because remaining candidate(s) are outside selection guards"` | Dormant break |
| 4 | 879–900 | `candidate_within_selection_budget(…)` helper | Sole call site is line 1224 |
| 5 | 341 | `within_budget: bool` field on `BubbleCandidate` | Used only by the guards above |
| 6 | 366 | `bailed: Vec<BubbleCandidate>` field on `CandidateFrontier` | Used only by the guards above |
| 7 | 1874–1880 | AllWave `max_pair_alignments` cap → returns `Err`, replacement dropped | Live cost cap (default 10k) |
| 8 | 1907–1913 | AllWave `max_replacement_paf_bytes` cap → returns `Err` | Live cost cap (default 64 MiB) |
| 9 | 1975–1981 | SweepGA `max_pair_alignments` cap → returns `Err` | Live cost cap (default 10k) |
| 10 | 1982–1988 | SweepGA `max_replacement_paf_bytes` cap → returns `Err` | Live cost cap (default 64 MiB) |

Items #1–#6 are dormant under default config but they are the dead-code
infrastructure for quality-bails of the family the task asked us to delete,
so they were removed alongside the live ones.

`path_sequences_equal` (`src/resolution.rs:2818`) is **kept** — that is the
correctness invariant, not a bail.

Build-time error propagation (e.g. BiWFA CIGAR over-consumption checks,
unsupported scoring params for POASTA, validation that the replacement
graph contains the expected per-traversal paths) is **kept** — those
indicate genuine algorithm bugs, not quality decisions.

## What was removed

In `src/resolution.rs`:

1. **`within_budget` field** removed from `BubbleCandidate` (was line 341).
2. **`bailed` vec** removed from `CandidateFrontier` (was line 366).
3. **Selection-guard loop branch** removed from `find_candidate_frontier`
   (was lines 1271–1275).
4. **Per-candidate `within_budget` initialization** removed at the build
   site of `BubbleCandidate` (was lines 1223–1224, 1240).
5. **Round-loop bailed counter and early-break** removed
   (was lines 491–494 and 496–504). The break condition on line 460 (no
   selected candidates remaining at all) is kept — this is the natural
   termination when `find_candidate_frontier` returns nothing.
6. **`candidate_within_selection_budget` helper** removed (was lines 879–900).
7. **AllWave pair cap** removed (was lines 1874–1880).
8. **AllWave PAF byte cap** removed (was lines 1907–1913).
9. **SweepGA pair cap** removed (was lines 1975–1981).
10. **SweepGA PAF byte cap** removed (was lines 1982–1988).

The `max_pair_alignments` and `max_replacement_paf_bytes` config fields and
their CLI flags are **kept** for back-compat (CLI smoke tests assert on the
parsed values), but the algorithm no longer reads them.

`stats.bailed` is still emitted in the run summary; it now only counts
build-time failures (empty replacement or backend error) — not quality
decisions, since none remain.

### Round-loop log line cleanup

The discovery summary log no longer mentions `selection-guarded`:

```text
// before
"crush round {}: {} POVU site(s), {} unseen polymorphic candidate(s), {} selected, {} selection-guarded in {:.2?}"
// after
"crush round {}: {} POVU site(s), {} unseen polymorphic candidate(s), {} selected in {:.2?}"
```

The traversal-stats log dropped its second half (`format_candidate_length_summary("selection-guarded", &frontier.bailed)`) for the same reason.

## Tests changed (and why)

| Test | Before | After |
|---|---|---|
| `direct_budgets_do_not_filter_pairwise_induction_candidates` | Asserted `candidate_within_selection_budget(...)` returned expected booleans for a mix of method/span/budget cases. | **Deleted.** The helper no longer exists. |
| `allwave_pair_alignment_cap_bails_before_seqwish` → renamed to `allwave_processes_candidate_regardless_of_pair_alignment_cap` | Asserted `resolved == 0, bailed == 1` when `max_pair_alignments=1`. | Now asserts `resolved == 1, bailed == 0` — the cap is no longer enforced, so the candidate is processed. |
| `direct_poa_bails_out_when_candidate_exceeds_median_traversal_budget` → renamed to `direct_poa_processes_candidate_regardless_of_median_traversal_budget` | Asserted `resolved == 0, bailed >= 1` when `max_median_traversal_len=4` cut off the candidate. | Now asserts `resolved == 1, bailed == 0`. POA runs on every candidate now. |
| `direct_poa_bails_out_when_candidate_exceeds_budget` → renamed to `direct_poa_processes_candidate_regardless_of_max_traversal_len` | Asserted `resolved == 0, bailed >= 1, gfa == input` when `max_traversal_len=4`. | Now asserts `resolved == 1, bailed == 0`. |
| `over_budget_candidates_are_not_marked_seen` → renamed to `all_candidates_processed_regardless_of_budget` | Asserted that over-budget candidates remain unseen across rounds (`resolved == 1, bailed >= 2`). | Now asserts `resolved >= 2, bailed == 0` — every candidate is processed in round 1 and signed off. |

The `BubbleCandidate { within_budget: true, … }` initialization was also
removed from 6 unit-test fixtures (a pure mechanical edit; the field no
longer exists on the struct).

All renamed tests had their assertion semantics inverted to match the new
"always process" behavior — none were silently weakened.

## Validation

### `cargo build --release`

```
Finished `release` profile [optimized] target(s) in 27.53s
```

(5 warnings, all pre-existing in vendored gfaffix — none from this change.)

### `cargo test --release --all`

```
test result: ok. 261 passed; 0 failed; 0 ignored
test result: ok. 56 passed; 0 failed; 0 ignored
test result: ok. 7 passed; 0 failed; 0 ignored
test result: ok. 25 passed; 0 failed; 0 ignored
test result: ok. 8 passed; 0 failed; 0 ignored
test result: ok. 10 passed; 0 failed; 0 ignored
(+ several smaller suites all green)
```

### Slice run — `tests/test_data/crush/c4_slice_1500_3000.gfa`

```
impg crush --gfa tests/test_data/crush/c4_slice_1500_3000.gfa \
  --method auto --max-iterations 10 -o /tmp/new_slice.gfa
```

Result: exit 0, 10 rounds completed. Run summary:

```
crush: 2027 resolved, 0 bailed, 2027 candidates seen across 10 rounds
```

`grep -E "bailed|reject|grew|stopping" /tmp/new_slice.stderr`:

```
[INFO  impg] crush: 2027 resolved, 0 bailed, 2027 candidates seen across 10 rounds
```

The only hit is the `0 bailed` in the summary counter. **No round-level
"bailed", "rejecting", "score grew", "stopping because …" lines.**

Per-round log no longer mentions "selection-guarded" (confirmed across
all 10 rounds — the new format is `… N selected in T`, with no second
clause):

```
crush round 1: 313 POVU site(s), 313 unseen polymorphic candidate(s), 147 selected in 527.93ms
crush round 2: 506 POVU site(s), 506 unseen polymorphic candidate(s), 163 selected in 569.82ms
...
crush round 10: 597 POVU site(s), 462 unseen polymorphic candidate(s), 202 selected in 505.22ms
```

`path_sequences_equal` still gates the rewrite (`src/resolution.rs:1557`):
the validate call is exercised on every round (`rewrite+validate Tms` in
each round's summary) and the run exits 0, confirming the invariant held.

Slice numbers match the pre-change baseline exactly (2027 / 0 / 2027) —
no candidate ever hit any of the dormant guards under default config,
which is why the slice never bailed in the first place. The point of
this task is to delete the *machinery*, not to change the slice's
numerical outcome.

### Full C4 GRCh38 — 5-min wall budget

```
timeout 300 impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O /tmp/crush-remove-bails/C4A.crush_full.nosort \
  -v 1
```

Result: **exit 124** (timeout, as expected per task spec). Progress at
cutoff:

```
crush round 1: 2441 POVU site(s), 72 unseen polymorphic candidate(s), 8 selected in 3.48s
crush round 1: building replacement 1/8 with Sweepga; traversals=312, max-len=31478, …
…
crush round 1: replacement build progress 1/8 (accepted)
crush round 1: replacement build progress 2/8 (accepted)
crush round 1: replacement build progress 3/8 (accepted)
crush round 1: replacement build progress 4/8 (accepted)
crush round 1: replacement build progress 5/8 (accepted)
crush round 1: replacement build progress 6/8 (accepted)
```

6 of 8 replacements built in 5 minutes; replacements 7 and 8 (which had
the largest `total-len` values, including a 5,132,443-bp Sweepga
all-vs-all between 312 long traversals) were still building when
`timeout` fired.

`grep -E "bailed|reject|grew|stopping|emitted unchanged" /tmp/crush-remove-bails/c4full.stderr`:

```
(no output — zero matches across all 446,368 stderr lines)
```

This is the F3-F4 cascade — long-running Sweepga inductions on
megabase-scale traversal sets — exactly the behavior the task description
predicted ("May hang due to the F3-F4 cascade (Phase 6 will fix this in
the next task). If it hangs, that's the cascade; document the last-emitted
line"). It is **not** a bail mechanism; nothing is being rejected. The
algorithm is simply slow on these inputs and would, given enough wall
time, accept replacements 7 and 8 just like 1–6.

Last-emitted log line in the C4 run at the 5-minute mark:

```
[2026-05-25T00:06:14Z INFO  impg::resolution] crush round 1: replacement build progress 6/8 (accepted)
```

(The trailing burst of `[sweepga::plane_sweep_scaffold]` chromosome-pair
lines and `[transclosure]` notices in `c4full.stderr` come from the
SweepGA backend's own progress reporting for the still-building 7th
replacement — not from crush itself.)

## Files modified

- `src/resolution.rs` — bail removal + test updates.

No CLI surface changed; the deprecated `--max-pair-alignments` and
`--max-replacement-paf-bytes` flags are still parsed (so existing
invocations keep working) but the algorithm no longer reads their values.

## Follow-ups (out of scope for this task)

- The unused config fields `max_pair_alignments` and
  `max_replacement_paf_bytes`, plus their CLI flags, could be deleted in
  a future cleanup. They are kept here so CLI smoke tests
  (`gfa:syng:crush,…,max-pair-alignments=20k,max-paf-bytes=128m,…`) keep
  parsing.
- The F3-F4 cascade is **`crush-level-descent`'s** scope (the downstream
  consumer task listed in this task's spec). The cascade is exposed but
  not addressed here.
