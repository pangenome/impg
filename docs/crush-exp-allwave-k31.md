# crush experiment — method=allwave, seqwish-k=31

Run of `impg query -o gfa:…crush,method=allwave,…seqwish-k=31…` on the
canonical C4/GRCh38 query. Recorded for the `crush-exp-allwave-k31` task,
sibling of `crush-exp-auto-k{31,51}` and `crush-exp-sweepga-k31`.

## Inputs

- Region: `GRCh38#0#chr6:31891045-32123783` (C4 SV locus on chr6)
- Index: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng`
- Sequences: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc`
- Distance threshold: `-d 50k`
- Threads: 32
- impg: `target/release/impg` v0.4.1, built on this branch (no rebuild for
  this task)
- Output directory: `/home/erikg/impg/data/c4_exp_allwave_k31_20260525T202440Z/`

## Recipe (verbatim, as run)

```bash
out=/home/erikg/impg/data/c4_exp_allwave_k31_20260525T202440Z
/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=allwave,min-traversal-len=5k,\
max-rounds=until-done,seqwish-k=31,max-pair-alignments=0,max-paf-bytes=0,\
polish-rounds=until-done,polish-max-traversal-len=10k,\
polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1 > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"
gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/run.Ygs.png" -m -x 2200 -y 1200
scp "$out/run.Ygs.png" erik@hypervolu.me:www/impg/c4-crush-allwave-k31.png
```

## Result — top line

| | |
|---|---|
| **Exit status** | `0` |
| **Wall (impg)** | **40 min 11 s** (`Elapsed (wall clock)` = `40:11.26`) |
| **Syng query wall** | 36 min 37.900 s (impg's own end-to-end timer) |
| **Crush phase wall** | 35 min 4 s (sum of round totals: 2023.18 s + 81.73 s) |
| **Max RSS** | **106.7 GB** (`Maximum resident set size` = 106 730 548 kB) |
| **CPU utilisation** | 538 % (≈ 5.4 cores effective on a 32-thread budget) |
| **User CPU** | 11 945.85 s |
| **System CPU** | 1 042.79 s |
| **Output GFA (`run.nosort.gfa`)** | 43 886 837 B, **465 paths preserved** |
| **Sorted GFA (`run.Ygs.gfa`, `gfasort -p Ygs`)** | 28 123 673 B |
| **Rendered PNG (`run.Ygs.png`)** | 1 913 940 B, 2200×1200, uploaded as `c4-crush-allwave-k31.png` |
| **Crush resolved** | 9 sites (round 1: 8/8; round 2: 1/1) |
| **Crush bailed** | 0 |

All 465 paths are preserved in `run.nosort.gfa` (`grep -c '^P\t'` = `465`),
satisfying the "all 465 paths preserved" gate.

## Per-round timing — and the bottleneck plan

Crush dispatches the 8 round-1 plans in parallel under rayon. The 32 threads
are shared by 8 plans (each AllWave job rayon-parallelises its internal pair
schedule). The 8 plans accept-times below are wall-clock from a single
`building 1..8/8` dispatch line at `20:55:20Z`:

| accept # | wall (relative to 20:55:20Z) |
|---------:|----:|
| 1 | 25 s |
| 2 | 3 min 15 s |
| 3 | 3 min 20 s |
| 4 | 3 min 23 s |
| 5 | 3 min 23 s |
| 6 | 3 min 38 s |
| 7 | 6 min 4 s |
| **8** | **33 min 37 s** |

i.e. **seven of the eight plans complete inside 7 minutes; one plan
single-handedly drags round 1 to 33.6 minutes**. Without per-plan ID logging
we cannot pin which input was the laggard, but by elimination it is one of
the bimodal plans (plans 2 / 6 / 7 / 8, all `median ∈ [112, 165]`, `max ∈
[25 230, 42 362]`). Plan 1 (long-uniform: `median=25 155, max=31 478, count=312`)
is more likely the round-1 winner that returned at 25 s, since seqwish-k=31
on long, uniformly-divergent traversals is favourable for AllWave.

Round-2 plan 1/1 (`traversals=236, max-len=42237, median-len=33, total-len=2 832 930`)
is the same SV-locus shape (the residual short-allele/long-allele bubble)
and resolved in 75 s. Round 3 found no eligible candidates → terminated.

## Quality-score progression

The crush quality score (working-set integrated divergence; lower = better)
moved as:

```
round 1:  217 461 559  →  204 893 738   (Δ = −12 567 821, −5.78 %)
round 2:  204 893 738  →  205 009 314   (Δ = +115 576,   +0.06 %)
round 3:  no eligible candidates
```

The +0.06 % bump in round 2 is within rounding (the rewrite added a few
edges and 100 segments). Net improvement: **−5.7 % vs the input working set**.

Working-set bookkeeping (`ws-…`, larger = worse — total span × per-edge cost):

| metric | input | post-r1 | post-r2 |
|---:|---:|---:|---:|
| segments | 18 048 | 20 541 | 20 641 |
| segment-bp | 389 354 | 673 543 | 679 663 |
| links | 20 933 | 24 807 | 25 040 |
| path-steps | 4 591 855 | 4 196 122 | 4 192 232 |
| ws-total | 24 769 855 962 | 19 361 250 168 | 19 509 188 868 |
| ws-p99 | 177 799 | 170 062 | 170 062 |
| ws-max | 388 845 | 664 177 | 670 088 |
| ws-long (≥10 kb) | 208 494 | 140 713 | 140 685 |

The working set grew (segments, segment-bp, links — the graph got bigger)
while the path-traversal-weighted cost (`ws-total`, `ws-long≥10kb`) fell.
That is exactly the expected crush trade-off: spend extra segments to
factor out the long path-step penalty.

## Polish phase

No polish work was emitted: with `polish-max-traversal-len=10k` and
`polish-max-median-traversal-len=1k`, every remaining bubble after round 2
was either above the per-traversal cap or the median cap, so polish was a
no-op. Round 3 (`no eligible candidates`) confirms the structure is no
longer cuttable under the configured thresholds.

## Comparison with the prior C4 method-auto profile

Reference: `docs/crush-aligner-speed-study.md` characterises the same
canonical query under `method=auto, seqwish-k=311`. Salient deltas:

| | auto, k=311 (prior study) | allwave, k=31 (this run) |
|---|---|---|
| Round-1 wall | 838.34 s (~14 min) | 2 023.18 s (~33.6 min) |
| Round-1 outlier plan wall | 831.92 s (sPOA on plan 2) | ≈ 32 min (one of the bimodal plans) |
| Plans built in round 1 | 8 (sPOA / POASTA / SweepGA mix routed by median) | 8 (AllWave for every plan, by `method=allwave`) |
| Aligner | per-plan `auto_method_by_median(...)` | uniform AllWave |
| Seqwish k | 311 | 31 |
| Max RSS | (~99 GB for sPOA process-wide, study-time) | 106.7 GB |
| Paths preserved | 465 | 465 |

Two things matter here:

1. **Forcing AllWave on every plan does not avoid the bimodal SV
   bottleneck**, it merely moves it from sPOA to AllWave. AllWave on the
   bimodal C4 plan is faster than sPOA's pathological 832 s
   (`crush-aligner-speed-study.md` § "Slow plan — input characterization"),
   but it is still the rate-limiter — one plan takes >30 min while every
   other plan finishes in <7 min.

2. **`seqwish-k=31` is more sensitive than `k=311`** — it accepts more,
   shorter co-linear matches when inducing the GFA from the PAF. That makes
   AllWave's induced graph larger (more segments) but more
   topologically reduced (`ws-long≥10kb` dropped from 208k to 141k in r1).
   On this fixture the **`-5.7 %` total quality improvement is comparable
   to what the prior `method=auto` runs achieved**, but at ~2.4× the wall
   time and roughly the same RSS ceiling.

In short: **`method=allwave, seqwish-k=31` finishes, preserves all 465
paths, and improves quality, but does not unblock the C4 SV plan.** The
real win (POASTA's 9.93 s on the same fixture, from the aligner-speed-study)
required routing — not a uniform switch.

## Artifacts

- `data/c4_exp_allwave_k31_20260525T202440Z/run.nosort.gfa` — raw crush GFA, 465 paths
- `data/c4_exp_allwave_k31_20260525T202440Z/run.nosort.stderr` — full impg stderr (1.68 MB)
- `data/c4_exp_allwave_k31_20260525T202440Z/run.Ygs.gfa` — `gfasort -p Ygs` (Y, groom, topo-sort)
- `data/c4_exp_allwave_k31_20260525T202440Z/run.Ygs.png` — `gfalook` render
  (2200×1200, 20 641 segments, 465 paths, 25 040 edges; mirrored to
  `https://hypervolu.me/~erik/impg/c4-crush-allwave-k31.png`)

## Hard-gate checklist

- [x] Real run completed, exit code recorded (`Exit status: 0`), all 465 paths preserved (`grep -c '^P\t' run.nosort.gfa` = 465)
- [x] PNG uploaded as `c4-crush-allwave-k31.png` (1 913 940 B, scp to `erik@hypervolu.me:www/impg/c4-crush-allwave-k31.png`)
- [x] `docs/crush-exp-allwave-k31.md` committed (this file)
- [x] `wg artifact crush-exp-allwave-k31 docs/crush-exp-allwave-k31.md` recorded
