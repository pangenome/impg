# crush-exp-sweepga-k31 — Test: lower `seqwish-k` to rescue small-bubble plans

**Status:** hypothesis **REFUTED** by bit-identical PAF retention.

**Date:** 2026-05-25
**Reference:** C4 locus `GRCh38#0#chr6:31891045-32123783` on HPRCv2 (465 input paths).

## Hypothesis under test

From the task description and `src/syng_graph.rs:481`:

> `min_block_length = syncmer_len_u64 = seqwish-k`. The plane sweep filter
> rejects any alignment shorter than `seqwish-k`. With `seqwish-k=311`,
> alignments under 311 bp are dropped → sweepga fails on small bubbles by
> construction.

> Rerun with `method=sweepga + seqwish-k=31` (10× smaller threshold). Should
> allow sweepga to succeed on small-median bubbles that the default scaffold
> filter was killing.

**Expected:** per-plan scaffold-filter PAF retention rises from 0–2.4 % to
near 100 % on the 4 small-median plans (plans 3, 6, 7, 8 with median
112–165 bp).

## Headline result

Per-plan scaffold-filter PAF retention is **bit-identical** between
`seqwish-k=31` and the `seqwish-k=311` baseline. The retention rates on the
4 small-median plans stayed at **0.00–2.41 %**, exactly matching the
baseline. The hypothesis is refuted.

| graph_build dir (k=31) | baseline dir (k=311) | raw PAF | filtered PAF | retention | `cmp` raw | `cmp` filtered |
|---|---|---:|---:|---:|---|---|
| `graph_build_0000` (plan median 157 bp, traversals 98) | `graph_build_0001` | 9 786 | 2 | **0.02 %** | identical | identical |
| `graph_build_0001` (plan median 165 bp, traversals 47) | `graph_build_0000` | 2 301 | 0 | **0.00 %** | identical | identical |
| `graph_build_0002` (plan median 112 bp, traversals 325) | `graph_build_0002` | 85 673 | 30 | **0.04 %** | identical | identical |
| `graph_build_0003` (plan median 157 bp, traversals 437) | `graph_build_0004` | 183 265 | 4 422 | **2.41 %** | identical | identical |
| `graph_build_0004` (plan median 6 210 bp, traversals 287) | `graph_build_0005` | 82 369 | 82 082 | 99.65 % | identical | identical |
| `graph_build_0005` (plan median 3 319 bp, traversals 356) | `graph_build_0003` | 126 216 | 122 850 | 97.33 % | identical | identical |
| `graph_build_0006` (plan median 25 155 bp, traversals 312) | `graph_build_0006` | 77 666 | 47 770 | 61.51 % | identical | identical |
| `graph_build_0007` (plan median 7 286 bp, traversals 460) | `graph_build_0007` | 207 928 | 207 474 | 99.78 % | identical | identical |

(`graph_build_NNNN` naming is non-deterministic across runs — order reflects which parallel build finishes first; mapping above is by raw PAF size, which is unique per plan. Every k=31 raw.paf and filtered.paf is byte-identical to its k=311 counterpart, verified by `cmp -s`.)

**Aggregate, 4 / 8 plans fragment:** plans with median traversal length
112–165 bp keep **0.00–2.41 %** of sweepga PAF after the 1:1 scaffold
filter; plans with median ≥ 3 319 bp keep **61.5–99.8 %**. Identical to
the k=311 baseline (`docs/crush-fragment-source-trace.md` table).

## Why the hypothesis is wrong

`src/syng_graph.rs:481` (`min_block_length = syncmer_len_u64`) is in the
**syng-anchor plane sweep** — a different code path that filters syncmer
hits for the "method=syng" induction. It does *not* run on the external
sweepga PAF.

For `method=sweepga`, the per-plan PAF goes through
`src/commands/graph.rs::filter_generated_paf` (line 671), which builds a
`FilterConfig` via sweepga's own
`filter_config_from_align_cfg(align_cfg, avg_seq_len)`. The `min_block_length`
in that path is derived from the alignment configuration (`min_map_length`)
and the per-plan `avg_seq_len` through the adaptive `clamp_scaffold_params`
helper (`sweepga::pansn:207-225` referenced in
`docs/crush-fragment-source-trace.md`). It is **not** a function of
`seqwish-k`. Lowering `seqwish-k` from 311 to 31 therefore cannot change
what the scaffold filter accepts.

The downstream seqwish induction step (`build_gfa_from_paf_and_sequences`
in `src/syng_graph.rs`) *is* affected by `seqwish-k`: see the per-plan
`seqwish.gfa` size comparison below — small-median plans are essentially
unchanged because they had so few PAF lines to begin with, but
medium/large-median plans produce up to ~20 % more seqwish output bytes at
k=31 (more, finer segments admitted by the lower minimum match length).

| graph_build dir | plan median | k=31 seqwish.gfa | k=311 seqwish.gfa | Δ bytes | Δ % |
|---|---:|---:|---:|---:|---:|
| 0000 ↔ 0001 | 157 bp |  45 919 |  45 974 |    −55 | −0.12 % |
| 0001 ↔ 0000 | 165 bp |  34 630 |  34 630 |      0 |  0.00 % |
| 0002 ↔ 0002 | 112 bp |  77 534 |  77 417 |   +117 | +0.15 % |
| 0003 ↔ 0004 | 157 bp | 160 331 | 155 737 | +4 594 | +2.95 % |
| 0004 ↔ 0005 | 6 210 bp | 155 617 | 130 107 | +25 510 | **+19.61 %** |
| 0005 ↔ 0003 | 3 319 bp |  94 355 |  78 197 | +16 158 | **+20.66 %** |
| 0006 ↔ 0006 | 25 155 bp | 513 004 | 431 703 | +81 301 | **+18.83 %** |
| 0007 ↔ 0007 | 7 286 bp | 346 725 | 300 928 | +45 797 | **+15.22 %** |

Interpretation: `seqwish-k=31` mostly affects plans whose retained-PAF
volume is already large; small-median fragmented plans are unaffected
because no extra PAF lines exist to feed lower k. The graph just gets
slightly more fine-grained where retention was already high.

## Final graph comparison

| metric | k=31 | k=311 baseline | Δ |
|---|---:|---:|---:|
| input paths preserved | 465 / 465 | 465 / 465 | 0 |
| final GFA segments (post-`gfasort -p Ygs`) | 20 848 | 20 786 | +62 (+0.30 %) |
| final GFA segment-bp | 687 177 | 697 267 | −10 090 (−1.45 %) |
| final GFA links | 25 253 | 25 141 | +112 (+0.45 %) |
| final GFA paths | 465 | 465 | 0 |
| crush rounds | 1 | 1 | — |
| crush resolved | 8 | 8 | — |
| crush bailed | 0 | 0 | — |

The fragmented graph shape on small-median bubbles is unchanged.

## Run

Wall time **5 min 33.54 s**, max RSS **53.7 GiB**, exit 0.

```bash
out=/home/erikg/impg/data/c4_exp_sweepga_k31_20260525T204733Z
IMPG_CRUSH_DEBUG_DIR=/tmp/c4-sweepga-k31-debug \
/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=31,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1 \
  > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"
gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/run.Ygs.png" -m -x 2200 -y 1200
scp "$out/run.Ygs.png" erik@hypervolu.me:www/impg/c4-crush-sweepga-k31.png
```

The `IMPG_CRUSH_DEBUG_DIR` is **not** in the task command — it was added to
capture per-plan `{raw,filtered}.paf` for the headline retention metric
(`src/syng_graph.rs:611-653` writes these when the env var is set).

Selected `/usr/bin/time -v` metrics:

| metric | value |
|---|---|
| User time (s) | 2 717.95 |
| System time (s) | 192.33 |
| Percent of CPU | 872 % |
| Elapsed wall (m:ss) | **5:33.54** |
| Maximum resident set size (KB) | **56 352 368** (≈ 53.7 GiB) |
| Minor page faults | 23 366 210 |
| File system outputs | 26 363 728 |
| Exit status | **0** |

Selected crush summary lines (from stderr):

```
crush: parsed input GFA into 18048 segment(s), 465 path(s) in 172.45ms
crush round 1: building replacement 1/8 with Sweepga; traversals=312, max-len=31478, median-len=25155, total-len=5132443, root-span=110
crush round 1: building replacement 2/8 with Sweepga; traversals=356, max-len=9320, median-len=3319, total-len=1181720, root-span=3319
crush round 1: building replacement 3/8 with Sweepga; traversals=460, max-len=7351, median-len=7286, total-len=3309531, root-span=7277
crush round 1: building replacement 4/8 with Sweepga; traversals=98, max-len=25255, median-len=157, total-len=65552, root-span=157
crush round 1: building replacement 5/8 with Sweepga; traversals=47, max-len=25230, median-len=165, total-len=32820, root-span=165
crush round 1: building replacement 6/8 with Sweepga; traversals=325, max-len=25274, median-len=112, total-len=188057, root-span=112
crush round 1: building replacement 7/8 with Sweepga; traversals=437, max-len=42362, median-len=157, total-len=2893782, root-span=157
crush round 1: building replacement 8/8 with Sweepga; traversals=287, max-len=6265, median-len=6210, total-len=1781621, root-span=6196
crush round 1: resolved 8/8 replacement(s) in 28.40s; rewrite+validate 2.98s; total 35.23s;
  quality score=217461559, segments=18048, segment-bp=389354, links=20933, path-steps=4591855,
    ws-total=24769855962, ws-p99=177799, ws-max=388845, ws-long>=10000bp=208494
  ->
  score=196395256, segments=20848, segment-bp=687177, links=25253, path-steps=4354666,
    ws-total=19263364347, ws-p99=160980, ws-max=676546, ws-long>=10000bp=145389
crush round 2: no eligible candidates from 2347 POVU site(s) in 3.35s (round 1 resolved 8 region(s))
crush: 8 resolved, 0 bailed, 8 candidates seen across 1 rounds
Syng query complete: GRCh38#0#chr6:31891045-32123783 in 2m6.286s
```

## Artifacts

- Output dir: `/home/erikg/impg/data/c4_exp_sweepga_k31_20260525T204733Z/`
  - `run.nosort.gfa` (45 689 386 B), `run.Ygs.gfa` (29 265 086 B),
    `run.Ygs.png` (1 890 534 B), `run.nosort.stderr` (69 610 121 B)
- Per-plan PAF dump: `/tmp/c4-sweepga-k31-debug/graph_build_{0000..0007}/{raw,filtered}.paf`
- PNG uploaded: <https://hypervolu.me/~erik/impg/c4-crush-sweepga-k31.png>
- Baseline (k=311) for comparison: `/tmp/c4-crush-debug/`,
  `/tmp/c4-true-descent-trace.gfa`, `docs/crush-fragment-source-trace.md`.

## Hard gates

- [x] Real run completed, all 465 paths preserved (verified `awk '/^P\t/' run.Ygs.gfa | wc -l`)
- [x] PNG uploaded as `c4-crush-sweepga-k31.png` (verified `ssh erik@hypervolu.me ls www/impg/c4-crush-sweepga-k31.png` → 1 890 534 B)
- [x] `docs/crush-exp-sweepga-k31.md` committed with per-plan retention, wall time, all metrics, PNG URL (this file)
- [x] Per-plan scaffold-filter retention reported as the headline — **0.00–2.41 % on 4 of 8 plans**, bit-identical to k=311; hypothesis refuted
- [x] `wg artifact crush-exp-sweepga-k31 docs/crush-exp-sweepga-k31.md`

## Recommendation

`seqwish-k` is the wrong knob for rescuing small-median bubbles under
`method=sweepga`. The actual fragmentation gate is the sweepga scaffold
filter's `min_block_length`, set adaptively from `avg_seq_len` via
`filter_config_from_align_cfg` + `clamp_scaffold_params`. Future experiments
that aim to rescue small-bubble PAF retention should target one of:

- `min_map_length` (the sweepga align config field that propagates to
  `min_block_length`), or
- `scaffold_filter` mode (currently `1:1`; relaxing to `M:N` or disabling
  via `--no-filter`), or
- the `clamp_scaffold_params` adaptive-shrink itself (currently shrinks
  `scaffold_mass` to ~419 bp for plan-7-class inputs).

The `method=auto` routing (`src/resolution.rs:1277-1291`) bypasses sweepga
for small-median plans and is the existing, working fix.
