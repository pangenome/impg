# crush-layer-a-validate

Validation run for the Layer A fixes described in `docs/crush-spec-audit.md`.
The two fixes under test are:

1. **Fresh-ID-range fix (invariant 4b)**: `next_id` should start at
   `original.segments.len()` instead of `1`.
2. **Three-tier routing fix (Phase 2)**: auto mode should route to sPOA / POASTA /
   sweepga by median traversal length, not just sPOA / Allwave by max length.

**Verdict: both fixes are absent. Both hard gates FAIL.**

---

## Validation environment

- Branch: `wg/agent-75/crush-layer-a-2`
- Binary: `/home/erikg/impg/target/release/impg` (built 2026-05-24T18:25 from
  commit `176136f`; two docs-only commits above `4545cf8` — source code unchanged)
- Input: `tests/test_data/crush/c4_slice_1500_3000.gfa`
  (2 942 segments, 465 paths, 64 348 bp; segment IDs in `[50064, 272213932]`)

---

## Run 1: slice, method=auto, max-iterations=1

```bash
time impg crush -g tests/test_data/crush/c4_slice_1500_3000.gfa \
  --method auto --max-iterations 1
```

**Wall time:** 21.66 s (45.12 s user / 1.11 s sys, 213% CPU)

### Stderr summary

```
crush round 1: 381 POVU site(s), 292 unseen polymorphic candidate(s),
  33 selected, 0 selection-guarded in 604.52ms
crush round 1 traversal stats: selected n=33,
  max-len median/max=116/20988, median-len median/max=116/20840,
  p90-len median/max=116/20877, traversals max=465,
  total max=9416057, step-savings max=469971, root-span max=166
crush round 1: rejecting 33/33 replacement(s); score grew 15.3%
  (threshold 10.0%); ws-p99 33508->31140, ws-max 64200->75498; total 21.56s
crush: 0 resolved, 33 bailed, 33 candidates seen across 1 rounds
```

All 33 selected candidates were **rejected** because the visual-tail score grew
15.3% (threshold 10.0%). The output GFA is bit-for-bit identical to the input.

### Graph metrics (output = input, all replacements rejected)

| metric | value |
|---|---|
| total segments | 2 942 |
| total segment bp | 64 348 |
| min segment ID | 50064 |
| max segment ID | 272 213 932 |
| segments with ID ≥ n\_orig (2942) | 2 942 (all original IDs are ≥ 2942) |
| new replacement IDs | **0** (no replacements accepted) |

Because all replacements were rejected, the ID-allocation code was never
invoked. Invariant 4b cannot be confirmed or denied from this run alone.

### Routing counts (auto mode, round 1, 33 selected)

Using the routing function at `src/resolution.rs:838-855`:
- Condition for POA: `max_len ≤ auto_spoa_max_traversal_len` (2000 bp)
- Otherwise: Allwave

| tier | spec threshold | code behaviour | count (estimated) |
|---|---|---|---|
| sPOA | median ≤ 1 kb | max\_len ≤ 2000 bp | ~32 (max-len 116 bp) |
| POASTA | 1–10 kb median | **not reachable** from auto | 0 |
| sweepga | > 10 kb median | **not reachable** from auto | 0 |
| Allwave | — (spec has no allwave tier) | max\_len > 2000 bp | ~1 (max-len 20988 bp) |

All 33 were ultimately rejected, so accepted-replacement tier counts are all 0.

---

## Run 2: prior audit data (slice-auto-1round, agent-62)

Agent-62 (`docs/crush-spec-audit.md`) ran the same command at `4545cf8` and
produced 147/147 accepted replacements. The output was retained at
`/tmp/crush-audit-62/slice_auto_1round.gfa` and is used here to test the
ID-allocation invariants (since the current run produced 0 new IDs).

### Graph metrics (agent-62 run, 147 replacements accepted)

| metric | input | output | delta |
|---|---|---|---|
| total segments | 2 942 | 2 984 | +42 |
| total segment bp | 64 348 | 60 987 | −3 361 |
| new IDs | — | 1 205 | — |
| dropped IDs | — | 1 163 | — |
| survivor IDs | — | 1 779 | — |

New ID range: `[1, 1205]`

### Invariant 4b check (agent-62 output)

| invariant | result |
|---|---|
| 4b: all new IDs ≥ n\_orig (2942) | **FAIL** — 0/1205 new IDs ≥ 2942 |
| new IDs below n\_orig | **1205/1205** (range [1, 1205]) |
| source: `render_rewritten_graph:2860` | `let mut next_id = 1usize;` — unfixed |

All 1205 replacement IDs fall in `[1, 1205]`, all below `n_orig = 2942`. The fix
(`let mut next_id = original.segments.len() + 1;` or `max_original_id + 1`) has
**not been applied**.

### Invariant 4a check (agent-62 output, 0 collisions)

No collision between new IDs `[1, 1205]` and survivor IDs (which are all ≥ 50064).
The `used_ids` set guard at `src/resolution.rs:2856-2859` prevents operational
collisions. **PASS** (operational, not structural).

### Invariant 2 check

Both the current run and the agent-62 run exited with code 0. The per-path
byte-equality check at `src/resolution.rs:1572` (`path_sequences_equal`) runs
after every replacement batch; accepted batches passed it. **PASS**.

---

## Full C4 GRCh38 run (5-minute budget)

The canonical command from `docs/c4-crush-handoff.md` requires:
- `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng` (does not exist in this environment)
- `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc` (does not exist in this environment)

The full GRCh38 run cannot be executed; the required data files are absent. The
handoff document (`docs/c4-crush-handoff.md`) records the baseline metrics from
the known-good run at `data/c4_crush_eval_20260523T140141Z/`:
18 089 segments, 422 007 bp, 465 paths.

---

## Hard validation gates

| gate | result | evidence |
|---|---|---|
| Invariant 4b holds on slice (new IDs ≥ n\_orig) | **FAIL** | 1205/1205 new IDs < 2942; code still has `next_id = 1usize` at `resolution.rs:2860` |
| Invariant 4a holds (no ID collision) | PASS | 0 collisions in agent-62 run |
| Invariant 2 holds (path-equality validation passes) | PASS | exit 0 on both runs; `path_sequences_equal` accepted all batches |
| Routing counts reported with numbers per tier | PASS | sPOA ~32, Allwave ~1, POASTA 0, sweepga 0 |
| Routing produces NON-ZERO counts for POASTA/sweepga | **FAIL** | neither tier is reachable from `method=auto`; routing at `resolution.rs:843-851` returns only `Poa` or `Allwave` |
| `docs/crush-layer-a-validate.md` committed | — (this file) |

---

## Root-cause summary

Layer A fixes have NOT been applied to the codebase.

### Fix 1: fresh-ID-range (1 LOC)

- **File:** `src/resolution.rs:2860`
- **Current:** `let mut next_id = 1usize;`
- **Required:** `let mut next_id = original_max_int_id + 1;` (or `original.segments.len() + 1`)
- **Status:** UNFIXED

### Fix 2: three-tier routing (≈30 LOC)

- **File:** `src/resolution.rs:838-855` (`candidate_selection_method`)
- **Current:** two tiers — POA if `max_len ≤ 2000`, else Allwave
- **Required:** three tiers — sPOA if `median_len ≤ 1000`, POASTA if `median_len ≤ 10000`, sweepga above; both POASTA and sweepga must be reachable from `method=auto`
- **Status:** UNFIXED — POASTA and sweepga are unreachable from auto mode

---

## Invariants 1 and 3 (Layer B / Layer C scope)

Reported per task spec; not expected to pass.

**Invariant 1** (no two segments share a sequence): the agent-62 run produced 948
duplicate-sequence segments out of 2984 (31.8%). Unfixed; requires Phase 3/4
cross-plan dedup from `docs/crush-architecture-spec.md`.

**Invariant 3** (each round shrinks the graph): current run had 0 net change
(all replacements rejected). Agent-62 run: segments 2942 → 2984 (+42, +1.4%
increase); bp 64348 → 60987 (−5.2% decrease). Segment count growth is the
cascade failure described in `docs/crush-audit.md`. Layer C scope.
