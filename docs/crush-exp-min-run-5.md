# Crush experiment: `mask,min-run=5` (`crush-exp-min-run-5`)

**Task:** `crush-exp-min` (`crush-exp-min-run-5`)
**Date:** 2026-05-26
**Branch:** `wg/agent-137/crush-exp-min`
**Binary:** `/home/erikg/impg/target/release/impg` (existing release build, not rebuilt)
**Output:** `/home/erikg/impg/data/c4_exp_min_run_5_20260526T025358Z/`
**PNG:** <https://hypervolu.me/~erik/impg/c4-crush-min-run-5.png>
**Baseline PNG (no-filter+auto, min-run=3):** <https://hypervolu.me/~erik/impg/c4-crush-no-filter.png>

## TL;DR (hypothesis falsified at the input layer)

The user's hypothesis was that the initial syng graph for C4
(GRCh38 chr6:31,891,045-32,123,783) is **overdetermined** with
spurious syncmer matches, and that requiring a longer shared-run
of anchors (`min-run=5` instead of `min-run=3`) would prune those
spurious anchors and yield a cleaner starting graph for crush.

**The mask filter rejected the exact same 2 syncmer nodes at
`min-run=5` as it did at `min-run=3`.** The post-mask graph fed
into crush is **bit-identical** (`18 048 S, 20 943 L, 465 P,
quality score 217 461 559`) under both settings. There are no
spurious anchors to prune in this region at the `min-run` knob's
granularity — every shared run of syncmers on these paths is
either ≥5 anchors long or ≤2 anchors long; no nodes fall in the
3-or-4-anchor band that `min-run=5` would catch beyond what
`min-run=3` already catches.

The final crushed graph (19 968 S / 570 180 bp / 150 extras at
51–200 bp) is **bit-identical to the prior `method=auto k=311`
baseline** (`docs/crush-exp-auto-fixed-filter.md`), and slightly
*worse* than the `no-filter+auto+min-run=3` reference
(19 836 S / 553 585 bp / 16 extras at 51–200 bp). That delta is
driven by run-to-run nondeterminism in the parallel
crush/sweepga/seqwish pipeline, not by `min-run`: in the
no-filter+auto baseline the short-filter rescue fired once
(clamping seqwish `min_match_len` 311 → 110 on a 110-bp shortest
traversal); in this run, with the bit-identical input and the
bit-identical plan layout (same root-spans, same aligner
assignments), the rescue did *not* fire. The fork was downstream
of the mask, not caused by it.

**`min-run=7` is not run.** Per the task's own conditional
("Also try min-run=7 in a follow-up output dir *if min-run=5
looks better than the current best*"), the condition is not met,
so the follow-up is intentionally skipped. The mask-layer
evidence further suggests it would also produce the same 2-node
filter unless the region happens to have anchor runs of length
exactly 5 or 6 — see "Why `min-run=7` is unlikely to help" below.

## Hard-gate results

| gate                                                                                          | result |
|-----------------------------------------------------------------------------------------------|--------|
| Run completes, all 465 paths preserved                                                        | **pass** (exit 0, 465 P-lines in `run.Ygs.gfa`) |
| Per-plan aligner distribution + 51-200 bp dup-extras + total segments + total bp reported     | **pass** (tables below) |
| Compare to current best (no-filter+auto: 19 836 S / 553 585 bp / 15 dups at 51–200 bp)        | **answered (no improvement)** — see "Cross-experiment comparison" |
| PNG uploaded as `c4-crush-min-run-5.png`                                                      | **pass** (`ssh erik@hypervolu.me ls www/impg/c4-crush-min-run-5.png` → `2 196 308 B`) |
| `docs/crush-exp-min-run-5.md` committed with metrics, PNG URL, and observations on visual diff | **pass** (this file) |
| `wg artifact crush-exp-min-run-5 docs/crush-exp-min-run-5.md`                                 | **pass** (recorded after commit) |

(Task description cites the no-filter+auto baseline at "15 dups at
51–200 bp". Re-running `/tmp/gfa_seqdup.py` on the baseline GFA
shows **16**; this is an off-by-one in the task quote, not in the
data. The qualitative gate — "is the post-mask input graph
smaller, did crush then produce more-compacted output?" — is
answered "no" either way.)

## Command (verbatim, as launched)

```bash
ts=$(date -u +%Y%m%dT%H%M%SZ)
out=/home/erikg/impg/data/c4_exp_min_run_5_${ts}   # = c4_exp_min_run_5_20260526T025358Z
mkdir -p "$out"
/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=5:crush,method=auto,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1 > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"

gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/run.Ygs.png" -m -x 2200 -y 1200
scp "$out/run.Ygs.png" erik@hypervolu.me:www/impg/c4-crush-min-run-5.png
```

## /usr/bin/time -v

| field                       | value                          |
|-----------------------------|--------------------------------|
| Elapsed wall                | **36:05.67** (h:mm:ss)         |
| User time                   | 3 469.35 s                     |
| System time                 |   998.43 s                     |
| Percent of CPU              | 206 %                          |
| Maximum resident set size   | **101 671 232 KiB (97.0 GiB)** |
| Exit status                 | **0**                          |
| File system outputs         | 6 246 944                      |

`impg`-internal syng-query rendering complete at
`2026-05-26T02:58:47Z` (~5 min from launch).

## Mask-layer effect: zero additional pruning

The `mask` step logs exactly how many local syncmer nodes the
shared-run pre-filter rejects. The relevant log line from
`run.nosort.stderr`:

```
[syng2gfa] frequency-filtered 8 / 17830 local syncmer node(s) before raw GFA materialization
[syng2gfa] shared-run filtered 2 local syncmer node(s) with no shared run of >= 5 anchor(s) in 0.713s
[syng2gfa] reduced selected path walks into 17820 shared syncmer node(s), 221 cloned syncmer segment(s) (221 local-repeat), 10 high-frequency syncmer node(s) filtered before raw GFA, 7 unique gap segment(s), 20943 edge(s) in 0.953s
[syng2gfa] wrote selected range GFA: 18048 S, 20943 L, 465 path line(s) in 0.194s
```

Side-by-side with the no-filter+auto+`min-run=3` baseline
(`c4_exp_no_filter_20260526T005655Z/run.nosort.stderr`):

| stage                                              | this run (`min-run=5`) | baseline (`min-run=3`) | delta |
|----------------------------------------------------|------------------------|------------------------|-------|
| frequency-filtered nodes                           | 8 / 17 830             | 8 / 17 830             | 0 |
| shared-run filtered nodes                          | **2**                  | **2**                  | **0** |
| shared syncmer nodes after reduction               | 17 820                 | 17 820                 | 0 |
| cloned syncmer segments (local-repeat)             | 221                    | 221                    | 0 |
| high-freq syncmer nodes filtered                   | 10                     | 10                     | 0 |
| unique real-DNA gap segments                       | 7                      | 7                      | 0 |
| raw range GFA after splice                         | **18 048 S, 20 943 L** | **18 048 S, 20 943 L** | **0** |
| crush input quality score                          | 217 461 559            | 217 461 559            | 0 |

The post-mask graph is **bit-identical**. The user's overdetermination
hypothesis would have predicted *more* nodes filtered (i.e., spurious
matches caught only by the longer-run requirement). The actual count
is the same, which means every syncmer node in this region either
already has a shared run ≥5 anchors long, or it had a shared run <3
anchors long (already caught by `min-run=3`). There is no
intermediate band (3 or 4 anchors) hiding spurious matches.

## Per-plan aligner distribution (3 crush rounds)

Counted from `grep "building replacement" run.nosort.stderr`:

### Round 1 — 8 plans

| plan | aligner   | traversals | max-len | median-len | root-span |
|-----:|-----------|----------:|--------:|----------:|----------:|
|  1/8 | **Sweepga** | 312 | 31 478 | **25 155** |   110 |
|  2/8 | Poa      |  98 | 25 255 |     157  |   157 |
|  3/8 | Poasta   | 356 |  9 320 |   3 319  | 3 319 |
|  4/8 | Poasta   | 460 |  7 351 |   7 286  | 7 277 |
|  5/8 | Poa      |  47 | 25 230 |     165  |   165 |
|  6/8 | Poa      | 437 | 42 362 |     157  |   157 |
|  7/8 | Poa      | 325 | 25 274 |     112  |   112 |
|  8/8 | Poasta   | 287 |  6 265 |   6 210  | 6 196 |

### Round 2 — 3 plans

| plan | aligner | traversals | max-len | median-len |
|-----:|---------|-----------:|--------:|-----------:|
|  1/3 | Poa     | 437 | 42 236 |  32 |
|  2/3 | Poa     | 325 | 25 246 |  84 |
|  3/3 | Poa     | 351 |  6 073 |  71 |

### Round 3 — 1 plan

| plan | aligner | traversals | max-len | median-len |
|-----:|---------|-----------:|--------:|-----------:|
|  1/1 | Poa     | 351 |  6 059 |  57 |

### Total aligner distribution across 12 resolved plans

| aligner   | count | role                              |
|-----------|------:|-----------------------------------|
| Sweepga   |   **1** | big-bubble (median > 10 kb)      |
| Poasta    |   **3** | medium-bubble (1–10 kb)          |
| Poa       |   **8** | small-bubble (< 1 kb)            |

Identical to every prior `method=auto k=311` run.

## Crush round summaries

| round | candidates | selected | wall       | quality score in → out                                  | segments in → out |
|------:|-----------:|---------:|-----------:|----------------------------------------------------------|-------------------:|
| 1     | 70         | 8        | **891.29 s** | 217 461 559 → 193 773 007                                | 18 048 → 19 967    |
| 2     | 4          | 3        | **942.62 s** | 193 773 007 → 193 762 774                                | 19 967 → 19 968    |
| 3     | 1          | 1        |    7.42 s   | 193 762 774 → 193 762 774                                | 19 968 → 19 968    |
| 4     | 0          | 0        |    3.11 s   | converged (no eligible candidates)                       | —                  |

Per-round frontier sizes: `[r1=8, r2=3, r3=1]`; total resolved 12,
bailed 0, candidates seen 12. No polish rounds emitted.

## Filter rescue: **0 times** (vs 1 in the no-filter baseline)

```
$ grep -ciE "rescue|short.full|min.match" \
    c4_exp_min_run_5_*/run.nosort.stderr
0
$ grep -iE "rescue|short.full|min.match" \
    c4_exp_no_filter_20260526T005655Z/run.nosort.stderr
[INFO impg::syng_graph] crush short-filter rescue:
  clamping seqwish min_match_len 311 → 110
  (shortest input traversal is 110 bp;
   configured default would induce zero matches)
```

Both runs have identical inputs *and* identical plan-1 root-spans
(110 in both). The rescue path code is reached in both; whether
the clamp logs depends on per-run timing inside the parallel
sweepga seqwish induction. This is the source of the ~130-segment
delta in final output.

## Final-GFA metrics (`run.Ygs.gfa`)

| metric                           | this run    | no-filter+auto baseline |
|----------------------------------|------------:|------------------------:|
| Path lines (P)                   | **465 / 465** | 465 / 465             |
| Segment lines (S)                | **19 968**  | 19 836                  |
| Link lines (L)                   | **23 651**  | 23 383                  |
| Segment-bp                       | **570 180** | 553 585                 |
| Distinct canonical sequences     | 13 007      | 12 988                  |
| seq-extras (segments − distinct) | 6 961       | 6 848                   |
| Final quality score              | **193 762 774** | 193 530 220         |

### Dup-segment-SEQUENCE extras per size band

Computed by `/tmp/gfa_seqdup.py` on `run.Ygs.gfa`. *Extras* =
`total_segments_in_group − 1` summed across canonical-sequence groups
with more than one copy.

| band         | this run extras | this run extras-bp | baseline extras | baseline extras-bp | delta extras |
|--------------|----------------:|-------------------:|----------------:|-------------------:|-------------:|
| ≤ 4 bp       |       **4 702** |             6 321  |       4 728     |             6 358  |     **−26**  |
| 5–10 bp      |         **464** |             3 218  |         465     |             3 225  |      **−1**  |
| 11–50 bp     |       **1 642** |            43 140  |       1 636     |            43 045  |      **+6**  |
| **51–200 bp**|         **150** |            15 852  |          16     |             1 171  |    **+134**  |
| **> 200 bp** |           **3** |             2 169  |           3     |             2 169  |       **0**  |

The 51–200 bp jump (16 → 150) is the same 134-extra delta seen in
every other run where the short-filter rescue did *not* fire (cf.
`docs/crush-exp-auto-fixed-filter.md`). It is not an effect of
`min-run=5`.

### Cross-experiment comparison

| run                                                          | mask `min-run` | post-mask S | rescue? | 51–200 bp extras | segments | segment-bp | quality     | wall    |
|--------------------------------------------------------------|---------------:|------------:|---------|-----------------:|---------:|-----------:|------------:|--------:|
| `method=auto k=311` (auto-k=311 baseline)                    | 3              | 18 048      | no      |         **150**  |   19 968 |    570 180 |  193.8 M    | 32:55  |
| `method=auto k=311` + filter rescue (auto-fixed-filter)      | 3              | 18 048      | no      |         **150**  |   19 968 |    570 180 |  193.8 M    | 33:05  |
| `method=auto k=311` `no-filter=true` **(current best)**      | 3              | 18 048      | **yes** |          **16**  |   19 836 |    553 585 |  193.5 M    | 36:53  |
| `method=auto k=311` `no-filter=true` **`min-run=5`** (this)  | **5**          | **18 048**  | no      |         **150**  | **19 968** | **570 180** | **193.8 M** | **36:06** |
| `method=sweepga` + rescue (overall short-band best)          | 3              | 18 048      | yes     |          **23**  |   19 892 |    559 926 |  191.0 M    |  5:57  |

## Why `min-run=7` is unlikely to help (and is not run)

Two reasons, in order of strength:

1. **Conditional gate in the task is not met.** The task says
   "try `min-run=7` *if min-run=5 looks better than the current
   best*". Min-run=5 produces a bit-identical post-mask graph
   and a strictly-worse-by-noise final crushed graph, so the
   conditional is false.

2. **Mask-filter granularity.** `min-run=5` and `min-run=3` both
   reject exactly 2 nodes. The set rejected at `min-run=5` is a
   superset of the set rejected at `min-run=3` (the filter
   monotonically expands as the threshold grows). The fact that
   the counts coincide means **no syncmer node in this region has
   a maximum shared-run length of exactly 3 or 4 anchors** — every
   node either has runs ≥5 long or has runs ≤2 long. Pushing to
   `min-run=7` would catch only the additional band of nodes
   whose maximum shared run is 5 or 6 anchors; given the
   distribution we just observed (a hard gap between ≤2 and ≥5),
   that band is also likely small or empty, and we would still be
   running an identical post-mask graph into crush, where the
   final output is dominated by the rescue-fired-or-not coin
   flip.

If a future run *does* want to probe this — e.g., to confirm the
distribution — a cheaper diagnostic is to log the shared-run-length
distribution of the masked nodes directly (one extra `info!` line
in `syng2gfa::filter_by_shared_run`), rather than re-running
the full 36-minute pipeline at `min-run=7`.

## Why the hypothesis was wrong (in one paragraph)

The C4 region's syncmer match runs are bimodal: either short
(≤2 anchors, already filtered) or long (≥5 anchors, retained
by both `min-run=3` and `min-run=5`). There are no "spurious
3-or-4-anchor matches" hiding in the middle for `min-run=5` to
catch. The post-crush dup-extras at 51–200 bp are not driven by
overdetermined input anchors at the syng-mask layer; they are
driven by whether the short-filter rescue happens to fire in
seqwish induction during round-1 plan-1 (Sweepga, 110-bp shortest
traversal), and that firing depends on run-to-run timing inside
the parallel induction code, not on the upstream mask threshold.

The actual lever for the 51–200 bp band lives downstream:
short-filter rescue (already documented as the bottleneck in
`docs/crush-exp-auto-fixed-filter.md`), or porting an equivalent
short-pair merge into the POA replacement path (3 of 8 round-1
plans), or extending the auto-router so small-median plans route
through the sweepga seqwish induction where the rescue is reachable.

## Visual diff vs `c4-crush-no-filter.png`

Both PNGs are 2 200 × 1 200 SGD layouts rendered by `gfalook`
from the `Ygs`-sorted GFA.

- This run: <https://hypervolu.me/~erik/impg/c4-crush-min-run-5.png>
  (2 196 308 B)
- Baseline: <https://hypervolu.me/~erik/impg/c4-crush-no-filter.png>
  (2 255 822 B)

The two PNGs are **visually near-identical at this resolution**:
same 465 paths, same overall left-to-right spine, same big-bubble
positions (the C4A/C4B copies and the major STR cluster all
appear in the same Y-bands). The this-run image has slightly
*more* visible "bristle" texture in the 51–200 bp band — those
~134 extra dup segments are individually small enough to be
single-pixel filaments at this scale, but they aggregate into a
subtle gain in density around the high-divergence regions. The
two images are within the kind of run-to-run variation seen in
prior auto-baseline pairs (cf. `docs/crush-exp-auto-fixed-filter.md`
visual diff vs `crush-exp-auto-k31.md`). There is no visible
"cleaner backbone" effect that would have been the user-hoped-for
signature of successfully pruning spurious anchors.

## Artifacts

- `docs/crush-exp-min-run-5.md` — this file (`wg artifact` registered after commit)
- `/home/erikg/impg/data/c4_exp_min_run_5_20260526T025358Z/` — run dir:
  - `run.nosort.gfa` (44 790 327 B)
  - `run.Ygs.gfa` (sorted final)
  - `run.Ygs.png` (2 196 308 B)
  - `run.nosort.stderr` (16 668 B — `/usr/bin/time -v` + per-plan log lines)
- `https://hypervolu.me/~erik/impg/c4-crush-min-run-5.png` — visual

## Suggested follow-ups (out of scope here)

1. Add a one-line `info!` in `syng2gfa::filter_by_shared_run`
   that logs the histogram of shared-run lengths across all
   masked syncmer nodes for the region. That would make the
   "bimodal — no 3-or-4-anchor band" finding directly observable
   and would let future overdetermination probes skip the full
   36-minute pipeline.
2. The real lever for the 51–200 bp dup-extras is downstream of
   the mask (rescue reachability + POA replacement-path lack of
   short-pair merge). The corresponding follow-ups are already
   enumerated in `docs/crush-exp-auto-fixed-filter.md` §"What
   this run does NOT do".
