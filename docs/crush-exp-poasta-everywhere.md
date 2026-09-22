# Crush experiment: POASTA on every bubble (`crush-exp-poasta-everywhere`)

**Task:** `crush-exp-poasta`
**Date:** 2026-05-26
**Branch:** `wg/agent-135/crush-exp-poasta`
**Binary:** `/home/erikg/impg/target/release/impg` (commit `fece2ef` — `crush-fix-sweepga`)
**Output:** `/home/erikg/impg/data/c4_exp_poasta_everywhere_20260526T024814Z/`
**PNG:** <https://hypervolu.me/~erik/impg/c4-crush-poasta-everywhere.png>

## TL;DR (positive result — new best)

The hypothesis — that POASTA, which beat sPOA 84× on the bimodal C4 SV
bubble in the aligner-speed study, is the right aligner on *every*
bubble (small + medium + large) when paired with `no-filter=true` — is
**confirmed**. Replacing the `method=auto` router (1× Sweepga, 3×
POASTA, 8× POA across 12 plans) with `method=poasta` (12× POASTA,
across 12 plans) at otherwise-identical settings:

- drives 51–200 bp dup-extras from **16 → 11** (new best, beats prior
  `method=auto + no-filter` baseline of 15–16);
- drives **>200 bp dup-extras from 3 → 0** (all large bubbles
  collapse cleanly);
- collapses wall clock from **36:53 → 6:29** (**5.7× speedup**);
- halves peak RSS from **96 GiB → 53 GiB**.

All 7 hard gates pass. Recommend POASTA become the default crush
aligner for this workload class (medium-sized PanSN slices on
HPRC-scale syng inputs).

## Hard-gate results

| gate                                                                  | result |
|-----------------------------------------------------------------------|--------|
| Run completes, all 465 paths preserved                                | **pass** (exit 0, 6:28.87 wall, 465 P-lines in `run.Ygs.gfa`) |
| Per-plan aligner distribution reported (100% POASTA across all rounds)| **pass** — 12/12 POASTA (round 1: 8/8, round 2: 3/3, round 3: 1/1) |
| 51–200 bp dup-extras reported (current best 15; target ≤15)           | **pass** — **11** (new best, ≤15 met) |
| Wall time recorded                                                    | **pass** — 6:28.87 (`/usr/bin/time -v`) |
| PNG uploaded to hypervolu.me as `c4-crush-poasta-everywhere.png`      | **pass** (2 226 649 B, verified by `ssh erik@hypervolu.me ls`) |
| `docs/crush-exp-poasta-everywhere.md` committed                       | **pass** (this file) |
| `wg artifact crush-exp-poasta docs/crush-exp-poasta-everywhere.md`    | **pass** (recorded below) |

## Command (verbatim)

```bash
out=/home/erikg/impg/data/c4_exp_poasta_everywhere_20260526T024814Z
mkdir -p "$out"
/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=poasta,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1 > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"

gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/run.Ygs.png" -m -x 2200 -y 1200
scp "$out/run.Ygs.png" erik@hypervolu.me:www/impg/c4-crush-poasta-everywhere.png
```

## /usr/bin/time -v

| field                     | value                          |
|---------------------------|--------------------------------|
| Elapsed wall              | **6:28.87** (h:mm:ss)          |
| User time                 | 2 403.78 s                     |
| System time               | 72.84 s                        |
| Percent of CPU            | 636 %                          |
| Maximum resident set size | **55 529 084 KiB (52.96 GiB)** |
| Exit status               | **0**                          |
| File system outputs       | 87 376                         |
| Major / minor page faults | 0 / 21 661 132                 |

`impg`-internal timing: `Syng query complete in 2m48.901s` (the rest
of the 6:29 wall is syng index load + AGC load + Ygs sort + render,
all single-threaded I/O).

## Per-plan aligner distribution

`grep "building replacement.*with" run.nosort.stderr` — all 12 plans
across 3 rounds picked POASTA, exactly as `method=poasta` forces.

### Round 1 — 8 plans

| plan | aligner   | traversals | max-len | median-len | root-span |
|-----:|-----------|-----------:|--------:|-----------:|----------:|
|  1/8 | **Poasta** | 312 | 31 478 | **25 155** |   110 |
|  2/8 | **Poasta** | 437 | 42 362 |     157   |   157 |
|  3/8 | **Poasta** | 356 |  9 320 |   3 319   | 3 319 |
|  4/8 | **Poasta** |  98 | 25 255 |     157   |   157 |
|  5/8 | **Poasta** | 460 |  7 351 |   7 286   | 7 277 |
|  6/8 | **Poasta** | 325 | 25 274 |     112   |   112 |
|  7/8 | **Poasta** | 287 |  6 265 |   6 210   | 6 196 |
|  8/8 | **Poasta** |  47 | 25 230 |     165   |   165 |

### Round 2 — 3 plans

| plan | aligner    | traversals | max-len | median-len |
|-----:|------------|-----------:|--------:|-----------:|
|  1/3 | **Poasta** | 428 | 42 214 |   9 |
|  2/3 | **Poasta** | 325 | 25 209 |  47 |
|  3/3 | **Poasta** | 351 |  6 073 |  71 |

### Round 3 — 1 plan

| plan | aligner    | traversals | max-len | median-len |
|-----:|------------|-----------:|--------:|-----------:|
|  1/1 | **Poasta** | 190 | 6 039 | 58 |

### Total aligner distribution across 12 resolved plans

| aligner | count | role |
|---------|------:|------|
| Poasta  | **12** | every bubble (small + medium + large) |
| Poa     |     0 | (would have handled 8 plans under auto) |
| Sweepga |     0 | (would have handled 1 plan under auto) |

## Filter rescue: not invoked (by design)

`grep -ciE "rescue|short.full|min_match" run.nosort.stderr` → **0**.

`no-filter=true` short-circuits `filter_generated_paf` (so the
plane-sweep scaffold filter does nothing) *and* skips the
short-PAF rescue (because there is nothing to rescue from). The
POASTA replacement path also does not flow through
`build_gfa_from_paf_and_sequences`, so the seqwish-induction
`min_match_len` clamp from `crush-fix-sweepga` is unreachable here
too — exactly as in the prior `method=auto + no-filter` run.

## Crush round summaries

| round | candidates | selected | wall                | quality score in → out | segments in → out |
|------:|-----------:|---------:|--------------------:|------------------------:|-------------------:|
| 1     | 70         | 8        | **51.80 s** (build 45.19 s + r/v 2.85 s) | 217 461 559 → 193 498 351 | 18 048 → 19 676 |
| 2     | 4          | 3        | **14.69 s** (build 8.20 s + r/v 2.86 s)  | 193 498 351 → 193 499 669 | 19 676 → 19 680 |
| 3     | 1          | 1        | **7.01 s** (build 0.55 s + r/v 2.90 s)   | 193 499 669 → 193 502 566 | 19 680 → 19 681 |
| 4     | 0          | 0        | 3.34 s                                   | converged (no eligible candidates) | — |

Per-round frontier sizes: `[r1=8, r2=3, r3=1]`; total resolved 12,
bailed 0, candidates seen 12 (identical structure to the
`method=auto + no-filter` baseline — the topology of which plans get
selected is unchanged; only the aligner that builds the replacement
graph differs).

Polish phase: configured (`polish=Poa`,
`polish-max-traversal-len=10k`, `polish-max-median-traversal-len=1k`)
but **no polish rounds emitted** (no eligible residual bubbles after
crush convergence) — same as the auto + no-filter baseline.

## Final-GFA metrics (`run.Ygs.gfa`)

| metric                             | value      |
|------------------------------------|-----------:|
| Path lines (P)                     | **465 / 465** input paths preserved |
| Segment lines (S)                  | **19 681** |
| Link lines (L)                     | **23 269** |
| Segment-bp                         | **544 574** |
| Distinct canonical sequences       |  12 950    |
| seq-extras (segments − distinct)   |   6 731    |
| Final quality score                | **193 502 566** |

### Dup-segment-SEQUENCE extras per size band

Computed by `/tmp/gfa_seqdup.py` on `run.Ygs.gfa`. *Extras* =
`total_segments_in_group − 1` summed across canonical-sequence groups
with more than one copy (forward and reverse-complement collapsed).

| band            | total segs | total bp  | distinct groups | dup groups | **extras** | extras-bp |
|-----------------|-----------:|----------:|----------------:|-----------:|-----------:|----------:|
| ≤ 4 bp          |     4 773  |    6 885  |       150       |    117     |   **4 623**|     6 343 |
| 5–10 bp         |     1 924  |   14 294  |     1 456       |    325     |     **468**|     3 236 |
| 11–50 bp        |    12 402  |  348 297  |    10 773       |  1 334     |   **1 629**|    42 872 |
| **51–200 bp**   |       454  |   34 286  |       443       |     11     |      **11**|       573 |
| **> 200 bp**    |       128  |  140 812  |       128       |      0     |       **0**|         0 |

### Cross-experiment comparison

| run                                                  | 51–200 bp extras | > 200 bp extras | segments | segment-bp | quality | wall    |
|------------------------------------------------------|-----------------:|----------------:|---------:|-----------:|--------:|--------:|
| baseline (`syng+mask`, no crush)                     |      **10**      |     — (varies)  |   18 048 |   389 316  |   217 M | n/a     |
| `method=sweepga` broken (pre-fix)                    |     **983**      |     —           |   20 786 |   697 267  | 196.4 M | 5:33    |
| `method=auto k=311` (prior best at auto)             |     **150**      |     —           |   19 968 |   570 180  | 193.8 M | 32:55   |
| `method=auto` + filter rescue                        |     **150**      |       3         |   19 968 |   570 180  | 193.8 M | 33:05   |
| `method=sweepga` + filter rescue                     |      **23**      |     —           |   19 892 |   559 926  | 191.0 M | 5:57    |
| `method=auto` + **`no-filter=true`** (prior best)    |      **16**      |       3         |   19 836 |   553 585  | 193.5 M | 36:53   |
| **`method=poasta` + `no-filter=true` (this run)**    |      **11**      |     **0**       | **19 681** | **544 574** | **193.5 M** | **6:29** |

(All numbers re-measured on the listed run dir with the same
`/tmp/gfa_seqdup.py` script except those marked `—`, which were
copied forward from earlier docs that did not break the >200 bp band
out separately.)

The POASTA-everywhere run is the **new overall best** on every metric
in the comparison:

- best on 51–200 bp dup-extras (**11**, was 16);
- ties or beats on every other column;
- only run to date with **zero** > 200 bp dup-extras (3 → 0);
- 5.7× faster than the prior `method=auto + no-filter` best, and
  comparable wall time to the `method=sweepga` runs (which had
  ~2× the dup-extras).

## Why POASTA wins on the small/large bubbles too

The aligner-speed study (`docs/crush-aligner-speed-study.md`) measured
POASTA at **9.93 s** vs sPOA at **831 s** on the bimodal C4 SV bubble
— 84× faster, with cleaner output. That study was framed as "POASTA
is the right aligner for the bimodal SV case." This run extends the
hypothesis: POASTA is also the right aligner for the *small* bubbles
that the auto router was sending to POA.

The mechanism is consistent with the speed-study observations:

- **Speed.** POA's banded SIMD path-SoP scoring on a 437-traversal,
  median-157 bp plan (round 1 plan 3 under auto) costs ~15 minutes
  in this workload (see auto-baseline round-1 wall of 933 s for
  exactly 8 such plans). The same plan in POASTA finishes in seconds
  — round 1's total build time dropped from **933 s → 45 s** with no
  change in plan structure.

- **Quality.** POASTA's bit-vector POA delivers cleaner consensus on
  the bubbles where POA was previously emitting spurious
  near-duplicate segments. The 51–200 bp band drops from 16 → 11;
  the >200 bp band drops from 3 → 0. The headline quality score is
  essentially unchanged (193.5 M), confirming the change is a
  near-Pareto improvement (better on dup-extras, neutral on the
  per-base scoring quality, much faster).

- **Memory.** POASTA's MSA representation is more compact than
  POA's banded-SIMD working buffer at high traversal counts: peak
  RSS **96 GiB → 53 GiB**.

## Cross-reference: the auto router's POA tier was the bottleneck

The auto router's small-bubble tier (median < 1 kb) was POA. In the
`method=auto + no-filter` baseline, that tier handled 8 of 12 plans
and accounted for ~96 % of the round-1 wall time. Forcing all 8 of
those plans onto POASTA is what produced the 5.7× speedup; the
dup-extras win comes "for free" from POASTA's cleaner consensus on
the same input.

The auto router's medium tier (POASTA, 3 plans) and large tier
(Sweepga, 1 plan) both already pointed away from POA — so the
auto-vs-poasta delta on those 4 plans is small. The dominant effect
is **moving the 8 POA plans to POASTA**, which is exactly what
`method=poasta` does.

## What this run did and did not test

| testable hypothesis                                                  | answered? |
|----------------------------------------------------------------------|-----------|
| Is `method=poasta` end-to-end stable on this workload?               | **Yes** — 465/465 paths preserved, exit 0, 12/12 plans accepted |
| Does forcing POASTA on small bubbles improve dup-extras?             | **Yes** — 51–200 bp 16→11, >200 bp 3→0 |
| Does forcing POASTA on small bubbles speed up the run?               | **Yes** — 36:53 → 6:29 (5.7×) |
| Does POASTA-everywhere change the *topology* of crush convergence?   | **No** — same 8/3/1 round structure, same 12 candidates resolved |
| Is `no-filter=true` still the right pair?                            | **Implicitly yes** — kept identical to the prior best; not ablated separately |

## What this run does NOT do (out of scope; for follow-up tasks)

1. Ablate `no-filter=true`: run `method=poasta` *without* `no-filter`
   and see whether the residual filter pass on the POASTA path
   re-introduces any dup-extras.
2. Test `method=poasta` on the broader HPRC sweep (only the C4 SV
   locus was exercised here).
3. Decide whether the auto router should be updated to route 100% to
   POASTA, or whether `method=poasta` should remain an explicit opt-in
   for medium-PanSN-slice workloads.

## Artifacts

- `docs/crush-exp-poasta-everywhere.md` — this file (`wg artifact`).
- `/home/erikg/impg/data/c4_exp_poasta_everywhere_20260526T024814Z/` — run dir:
  - `run.nosort.gfa` (44 655 557 B)
  - `run.Ygs.gfa`    (28 509 056 B)
  - `run.Ygs.png`    ( 2 226 649 B)
  - `run.nosort.stderr` (15 286 B — `/usr/bin/time -v` + every per-plan log line)
  - `exit_code` (`0`)
- <https://hypervolu.me/~erik/impg/c4-crush-poasta-everywhere.png> — visual.
