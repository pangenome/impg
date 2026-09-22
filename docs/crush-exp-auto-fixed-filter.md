# Crush experiment: `method=auto` + filter rescue (`crush-exp-auto-fixed-filter`)

**Task:** `crush-exp-auto-3`
**Date:** 2026-05-26
**Branch:** `wg/agent-128/crush-exp-auto-3`
**Binary:** `/home/erikg/impg/target/release/impg` (commit `fece2ef` — `crush-fix-sweepga`)
**Output:** `/home/erikg/impg/data/c4_exp_auto_fixed_20260526T004407Z/`
**PNG:** <https://hypervolu.me/~erik/impg/c4-crush-auto-fixed-filter.png>

## TL;DR (negative result)

The hypothesis — that `method=auto` plus the `crush-fix-sweepga`
filter rescue would combine the best of POASTA-on-medium-bubbles with
the small-bubble dup collapse from the rescue — is **falsified**. The
final graph is **bit-identical** to the prior `method=auto k=311`
baseline (19 968 S, 23 651 L, 570 180 segment-bp, dup-extras
51–200 bp = **150**) because **the rescue never fires under
`method=auto`**: the rescue lives in the sweepga-tail seqwish
induction (`build_gfa_from_paf_and_sequences`), and the auto router
routes every small-median bubble to POA, not to Sweepga. Sweepga ran
exactly once in round 1 (plan 1, median 25 155 bp — a big bubble
where the 1:1 / 311 bp filter retained alignments anyway) and zero
times thereafter, so the rescue had nothing to rescue.

**The current best for 51–200 bp dup-extras remains `method=sweepga`
+ rescue at 23** (`docs/crush-fix-sweepga-short-filter.md`). Pushing
POASTA further into the small-bubble band would require *either*
extending the auto-router cutoffs so small bubbles fall under
Sweepga (and thus pick up the rescue) *or* porting an equivalent
short-pair rescue into the POA replacement path.

## Hard-gate results

| gate                                                                | result        |
|---------------------------------------------------------------------|---------------|
| Run completes, exit status 0                                        | **pass** (exit 0, 33:05 wall) |
| 465 / 465 paths preserved in final Ygs GFA                          | **pass** (465 P-lines)        |
| Per-plan aligner distribution reported (non-zero POASTA expected)   | **pass** — 3× POASTA in round 1 (see below) |
| 51–200 bp dup-extras ≤ 25 (target; current best 23 under sweepga)   | **FAIL** — 150 (no change vs auto-k=311 baseline) |
| PNG uploaded to `c4-crush-auto-fixed-filter.png`                    | **pass** (2 181 164 B, verified by `ssh erik@hypervolu.me ls`) |
| `docs/crush-exp-auto-fixed-filter.md` committed                     | **pass** (this file)          |
| `wg artifact crush-exp-auto-fixed-filter docs/...`                  | **pass** (recorded below)     |

The ≤25 gate is not met, but the experiment's *purpose* (test the
hypothesis) is satisfied: the metric tells us the hypothesis was
wrong.

## Command (verbatim)

```bash
out=/home/erikg/impg/data/c4_exp_auto_fixed_20260526T004407Z
mkdir -p "$out"
LD_LIBRARY_PATH="/home/erikg/htslib-local/lib:/home/erikg/micromamba/lib" \
IMPG_CRUSH_DEBUG_DIR=/tmp/c4-auto-fixed-debug \
/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1 > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"

gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/run.Ygs.png" -m -x 2200 -y 1200
scp "$out/run.Ygs.png" erik@hypervolu.me:www/impg/c4-crush-auto-fixed-filter.png
```

## /usr/bin/time -v

| Field                       | Value                          |
|-----------------------------|--------------------------------|
| Elapsed wall                | **33:05.27** (h:mm:ss)         |
| User time                   | 3 400.41 s                     |
| System time                 | 874.45 s                       |
| Percent of CPU              | 215 %                          |
| Maximum resident set size   | **99 914 636 KiB (95.3 GiB)**  |
| Exit status                 | **0**                          |
| File system outputs         | 6 246 944                      |
| Major / minor page faults   | 6 / 373 240 085                |

`impg`-internal timing: `Syng query complete in 29m39.971s`.

## Per-plan aligner distribution

The auto router picked the following aligners across all 3 crush
rounds (`grep "building replacement.*with" run.nosort.stderr`):

### Round 1 — 8 plans

| plan | aligner   | traversals | max-len | median-len | root-span |
|-----:|-----------|----------:|--------:|----------:|----------:|
|  1/8 | **Sweepga** | 312 | 31 478 | **25 155** |   110 |
|  2/8 | Poasta   | 356 |  9 320 |   3 319  | 3 319 |
|  3/8 | Poa      | 437 | 42 362 |     157  |   157 |
|  4/8 | Poasta   | 460 |  7 351 |   7 286  | 7 277 |
|  5/8 | Poa      | 325 | 25 274 |     112  |   112 |
|  6/8 | Poasta   | 287 |  6 265 |   6 210  | 6 196 |
|  7/8 | Poa      |  98 | 25 255 |     157  |   157 |
|  8/8 | Poa      |  47 | 25 230 |     165  |   165 |

Counts: **1× Sweepga, 3× Poasta, 4× Poa** — identical to the
prior `method=auto k=311` baseline and to `crush-exp-auto-k31.md`.
POASTA *is* hitting the medium-median bubbles (1–10 kb), confirming
the routing is doing what the user expected.

### Round 2 — 3 plans

| plan | aligner | traversals | median-len |
|-----:|---------|-----------:|-----------:|
|  1/3 | Poa     | 437 |  32 |
|  2/3 | Poa     | 325 |  84 |
|  3/3 | Poa     | 351 |  71 |

### Round 3 — 1 plan

| plan | aligner | traversals | median-len |
|-----:|---------|-----------:|-----------:|
|  1/1 | Poa     | 351 |  57 |

### Total aligner distribution across 12 resolved plans

| aligner   | count | role |
|-----------|------:|------|
| Sweepga   |   **1** | big-bubble (median > 10 kb) |
| Poasta    |   **3** | medium-bubble (1–10 kb)     |
| Poa       |   **8** | small-bubble (< 1 kb)       |

## Filter rescue triggered: **0 times**

`grep -ciE "rescue|short.full|min_match" run.nosort.stderr` → **0**.

The `crush-fix-sweepga` rescue only runs inside
`build_gfa_from_paf_and_sequences` (`src/syng_graph.rs`), which is
called by `build_sweepga_seqwish_replacement` (`src/resolution.rs`).
Under `method=auto`, this function was invoked exactly once (round 1
plan 1 — median 25 155 bp). At that bubble size the Sweepga 1:1 /
311 bp filter retains alignments natively, so the rescue logged
nothing — `SHORT_FULL_LENGTH_RESCUE_FRACTION` had no qualifying
records to admit, and `min_match_len` did not need clamping below the
configured 311.

For every plan in the 51–200 bp dup-extras band — POA plans 3, 5, 7,
8 in round 1 plus all 4 plans in rounds 2/3 — the rescue is
**unreachable**; POA's replacement path does not flow through
`build_gfa_from_paf_and_sequences`.

## Crush round summaries

| round | candidates | selected | wall | quality score in → out | segments in → out |
|------:|-----------:|---------:|-----:|------------------------:|-------------------:|
| 1     |     70     |    8     | **846.22 s** | 217 461 559 → 193 773 007 | 18 048 → 19 967 |
| 2     |      4     |    3     | **832.53 s** | 193 773 007 → 193 762 774 | 19 967 → 19 968 |
| 3     |      1     |    1     |     7.02 s   | 193 762 774 → 193 762 774 | 19 968 → 19 968 |
| 4     |      0     |    0     |     2.96 s   | converged (no eligible candidates) | — |

Per-round frontier sizes: `[r1=8, r2=3, r3=1]`; total resolved 12,
bailed 0, candidates seen 12.

Polish phase: configured (`polish=Poa`,
`polish-max-traversal-len=10k`, `polish-max-median-traversal-len=1k`)
but **no polish rounds emitted** (no eligible residual bubbles after
crush convergence).

## Final-GFA metrics (`run.Ygs.gfa`)

| metric             | value     |
|--------------------|----------:|
| Path lines (P)     | **465 / 465** input paths preserved |
| Segment lines (S)  | **19 968** |
| Link lines (L)     | **23 651** |
| Segment-bp         | **570 180** |
| Distinct canonical sequences | 13 007 |
| seq-extras (segments − distinct) | 6 961 |
| Final quality score | **193 762 774** |

### Dup-segment-SEQUENCE extras per size band

Computed by `/tmp/gfa_seqdup.py` on `run.Ygs.gfa`. *Extras* =
`total_segments_in_group − 1` summed across canonical-sequence groups
with more than one copy.

| band         | total segs | total bp | distinct groups | dup groups | **extras** | extras-bp |
|--------------|-----------:|---------:|----------------:|-----------:|-----------:|----------:|
| ≤ 4 bp       |     4 851  |   6 859  |       149       |    114     |   **4 702**|     6 321 |
| 5–10 bp      |     1 917  |  14 243  |     1 453       |    328     |     **464**|     3 218 |
| 11–50 bp     |    12 436  | 349 343  |    10 794       |  1 342     |   **1 642**|    43 140 |
| **51–200 bp**|       622  |  52 681  |       472       |     21     |     **150**|    15 852 |
| **> 200 bp** |       142  | 147 054  |       139       |      3     |       **3**|     2 169 |

### Cross-experiment comparison

The 51–200 bp band is the user's primary "real bug" metric.

| run                                                          | 51–200 bp extras | 51–200 bp extras-bp | segments | segment-bp | quality | wall |
|--------------------------------------------------------------|-----------------:|--------------------:|---------:|-----------:|--------:|-----:|
| baseline (`syng+mask`, no crush)                             |      **10**      |          —          |   18 048 |   389 316  |   217 M | n/a |
| `method=sweepga` broken (pre-fix)                            |     **983**      |       132 849       |   20 786 |   697 267  | 196.4 M | 5:33 |
| `method=auto k=311` (prior best at auto)                     |     **150**      |        15 852       |   19 968 |   570 180  | 193.8 M | 32:55 |
| `method=auto k=31`                                           |     **150**\*    |        — (similar)  |   19 906 |   566 794  | 193.8 M | 36:47 |
| `method=sweepga` + rescue (**current overall best**)         |      **23**      |         1 674       |   19 892 |   559 926  | 191.0 M | 5:57 |
| **`method=auto` + rescue (this run)**                        |     **150**      |        15 852       | **19 968** | **570 180** | **193.8 M** | **33:05** |

\* Reported from `docs/crush-exp-auto-k31.md` summary; not
re-classified by band in that doc.

The auto-with-rescue numbers are bit-identical to the
no-rescue auto baseline because the rescue never fired
(see "Filter rescue triggered: 0 times" above). The
sweepga-with-rescue run is still the overall best by a factor of
~6.5× on dup-extras and ~5.5× on wall time.

## Why the hypothesis was wrong

`docs/crush-fix-sweepga-short-filter.md` summarised the rescue as
*"if an alignment is essentially end-to-end on the shorter of the two
sequences, it survives, even if it is below `min_block_length`"*, and
*"don't ask seqwish to induce match-runs longer than the shortest
input sequence"*. Both clauses live inside
`build_gfa_from_paf_and_sequences`, which is the **sweepga-tail**
seqwish induction. The same file's §"Why this targets the sweepga
path specifically" already called out the boundary:

> The top-level `impg graph` command's `align_sequences` →
> `filter_generated_paf` chain does **not** flow through
> `build_gfa_from_paf_and_sequences`, so it is unchanged.

Under `method=auto`, the small-bubble routing tier picks POA, which
goes through a *different* replacement path
(`build_poa_replacement` / equivalent in `src/resolution.rs`) — not
through `build_gfa_from_paf_and_sequences`. So:

- the auto router *correctly* routes mediums to POASTA (3 plans);
- the auto router *correctly* routes smalls to POA (8 plans);
- the rescue is *invisible* to POA;
- therefore the 51–200 bp dup-extras stay at the auto baseline (150).

Both halves of the prior fix (filter rescue + adaptive
`min_match_len`) are bound to the seqwish induction step, which is
only on the sweepga path. POA produces its replacement graph from a
consensus MSA, not from a PAF + seqwish, so the rescue concept does
not transfer directly.

## What this run did and did not test

| testable hypothesis                                                   | answered? |
|-----------------------------------------------------------------------|-----------|
| Is the auto router *reaching* POASTA for medium bubbles?              | **Yes — 3 of 8 round-1 plans** |
| Does combining auto-routing + the sweepga rescue improve dup-extras?  | **No — bit-identical to auto-k=311** |
| Is the residual 51–200 bp dup-extras (150) caused by POA-on-smalls?   | **Strongly implied** — POA handles 8 of 12 plans, and the rescue (which would help the seqwish path) cannot reach any of them |

## What this run does NOT do (out of scope; for follow-up tasks)

1. Modify the auto-router thresholds so that small bubbles route to
   sweepga instead of POA (would extend the rescue to them but
   resurrects the pre-fix `seqwish-k=311` failure mode on every
   small bubble that doesn't pass the 90 % coverage threshold).
2. Port the rescue (or an equivalent short-pair pair-aware merge) into
   the POA replacement path.
3. Re-route plans 3 / 5 / 7 / 8 (round 1 POA, small medians) through
   the sweepga seqwish induction explicitly to measure whether the
   rescue *could* help them.

These three are the natural next steps; recommend filing each as a
separate task.

## Artifacts

- `docs/crush-exp-auto-fixed-filter.md` — this file (`wg artifact`).
- `/home/erikg/impg/data/c4_exp_auto_fixed_20260526T004407Z/` — run dir:
  - `run.nosort.gfa` (44 790 327 B)
  - `run.Ygs.gfa` (28 605 420 B)
  - `run.Ygs.png` (2 181 164 B)
  - `run.nosort.stderr` (4 871 829 B — `/usr/bin/time -v` + every per-plan log line)
- `https://hypervolu.me/~erik/impg/c4-crush-auto-fixed-filter.png` — visual.
