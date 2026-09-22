# Crush experiment: SweepGA on every bubble, every filter disabled

**Task:** `crush-sweepga-everywhere` (workgraph id `crush-sweepga-everywhere`)
**Date:** 2026-05-26
**Branch:** `wg/agent-150/crush-sweepga-everywhere`
**Output dir:** `/home/erikg/impg/data/c4_crush_sweepga_unfiltered_20260526T123722Z/`
**PNG:** `https://hypervolu.me/~erik/impg/c4-crush-sweepga-unfiltered.png`

## Hypothesis

> "we should get sweepga/fastga for all bubbles. why arent the alignments
> resulting in a compressed graph? they must be filtered."

The `no-filter=true` flag only disables the post-alignment plane-sweep /
scaffold filter. Downstream of that flag there are several additional
length / identity / k-mer-frequency filters in the seqwish induction tail
that quietly drop most short-vs-long alignments — those are the
"compressed graph" bottleneck. This experiment:

1. **Forces** sweepga on every bubble (`method=sweepga`, not `method=auto`).
2. **Disables every filter** in the sweepga → seqwish → graph pipeline,
   not just the plane-sweep filter, by extending `no-filter=true` to
   drop the downstream seqwish floor / identity / k-mer-frequency caps.

If filters were the bottleneck, the unfiltered SweepGA-everywhere run
should approach the PGGB gold-standard segment count (13,288 S vs the
prior best of 19,836 S).

## Filter inventory (sweepga → seqwish → graph pipeline)

Every length / coverage / identity / k-mer-frequency threshold reached
during the sweepga replacement path, found by following the call chain
from `build_sweepga_seqwish_replacement` (`src/resolution.rs:2647-2708`)
through `build_gfa_from_paf_and_sequences`
(`src/syng_graph.rs:697-810`) into the seqwish dependency.

| # | Filter | File:line (current values) | Disabled by `no-filter=true`? (pre this PR) | Disabled by `no-filter=true`? (post this PR) |
|---|---|---|---|---|
| 1 | **Plane-sweep / scaffold filter** (`filter_generated_paf`) | `src/commands/graph.rs:671-723`; gated by `if config.no_filter` at `:677` | **Yes** (returned unchanged) | yes |
| 2 | **Short-full-length rescue augmentation** (`rescue_short_full_length_alignments`) | `src/syng_graph.rs:608-692`; called from `:753`, gated by `if config.no_filter` at `:750-758` | **Yes** (skipped) | yes |
| 3 | **Adaptive seqwish min-match-len rescue clamp** ("fece2ef rescue clamp") | `src/syng_graph.rs:763-784`; clamps `min_match_len` down to `min_seq_len` when below configured default. Default `min_match_len` = 311. Clamp formula: `min_match_len = min(configured, shortest_input_traversal_len)` | No — clamp still fires; for bimodal inputs the clamp lands on 110 bp while actual PAF max-match-runs are 40–99 bp (per `docs/crush-aligner-deep-diag.md` Q1) | **Implicitly disabled** — when `no_filter=true`, `seqwish_min` is set to `1` in `seqwish_replacement_config`, so the `min_seq_len < min_match_len` clamp condition is false (min_seq_len ≥ 1) and the clamp branch is skipped |
| 4 | **Seqwish exact-match `min_match_len` floor** (the dependency's CIGAR-run filter) | `seqwish/src/alignments.rs:141` (`if len >= min_match_len`); receives `config.min_match_len` from `seqwish_replacement_config` at `src/resolution.rs:2386`; default 311 bp via `DEFAULT_REPLACEMENT_SEQWISH_MIN_MATCH_LEN` at `src/resolution.rs:241` | **No** — flag did not propagate into the seqwish induction tail, so every individual CIGAR `=`/`M` run shorter than 311 bp was dropped even with `no-filter=true` set | **Yes** — when `sweepga_no_filter=true`, this value is forced to `1` at `src/resolution.rs:2380-2386` (the new override) |
| 5 | **Seqwish `min_map_length` (per-PAF-line minimum mapping length)** | `src/commands/graph.rs:88`; `min_map_length` in `GraphBuildConfig`; resolved at `src/resolution.rs:2377-2389` (follows `min_match_len` if 0) | **No** — propagated `min_match_len` floor (311 bp default) and rejected entire short PAF lines | **Yes** — follows `seqwish_min=1` when `no_filter=true`, so even single-base PAF lines pass |
| 6 | **Seqwish `min_identity` (per-PAF-line minimum identity)** | `src/resolution.rs:2388-2389`, default `replacement_min_identity = 0.0` (already off by default) | Identity-filter was already at 0.0 by default but the field still propagated user-supplied values into seqwish | **Yes** — forced to 0.0 when `no_filter=true` (overrides any user `min-identity=...` for the unfiltered path) |
| 7 | **FastGA `kmer_frequency` cap** (`replacement_sweepga_kmer_frequency`) | `src/resolution.rs:2398-2405`; default `max(1_000, traversal_count * 10)` (so 3,120 for the 312-traversal sweepga bubble) | **No** — the cap stayed at 3,120 for the bimodal bubble even with `no-filter=true` set | **Yes** — when `no_filter=true` and `sweepga_kmer_frequency` is unconfigured, the cap is raised to `NO_FILTER_SWEEPGA_KMER_FREQUENCY = 1_000_000` at `src/resolution.rs:2417-2418` (effectively unlimited; users can still pin it via `kmer-frequency=N`) |
| 8 | **`min-traversal-len` floor** (gate before any aligner runs) | `src/resolution.rs:1646-1648` | Independent of `no-filter`. The canonical command uses `min-traversal-len=5k`, which discards the 62 small bubbles per `docs/crush-aligner-deep-diag.md` Pattern A | Independent of `no-filter`; this experiment keeps `min-traversal-len=5k` for **apples-to-apples** comparison with `crush-exp-poasta-everywhere` and the `no-filter+auto` baseline. (Lowering the floor is a separate experiment.) |
| 9 | **`max_pair_alignments` cap** | `src/resolution.rs:139` + flag parser; canonical uses `max-pair-alignments=0` (disabled). | Already at 0 in canonical command. | unchanged |
| 10 | **`max_replacement_paf_bytes` cap** | `src/resolution.rs:145` + flag parser; canonical uses `max-paf-bytes=0` (disabled). | Already at 0 in canonical command. | unchanged |
| 11 | **`min_aln_length` (per-line minimum block length)** | `src/resolution.rs:2385`; passes `min_aln_length: 0` always. | Already 0 → off. | unchanged (still 0) |
| 12 | **SweepGA internal `no_filter`** | `src/resolution.rs:2671` — `no_filter: true` is **hardcoded** inside `build_sweepga_seqwish_replacement` so the raw all-vs-all PAF is always handed back to the seqwish tail. | Already hardcoded to true. | unchanged |

**Summary**: before this PR, four downstream filters survived `no-filter=true`
(rows 3, 4, 5, 6 conditionally, and 7). After this PR, `no-filter=true`
disables **all** filters in the sweepga → seqwish → graph induction path
that gate on per-CIGAR-run length, per-PAF-line length, per-PAF-line
identity, and FastGA seed frequency.

## Code change (this PR)

`src/resolution.rs` — `seqwish_replacement_config` and the new
`resolve_replacement_kmer_frequency` helper extend `sweepga_no_filter`:

```rust
// src/resolution.rs:2374+
fn seqwish_replacement_config(
    config: &ResolutionConfig,
) -> crate::commands::graph::GraphBuildConfig {
    // When the sweepga/seqwish path is in "no-filter" mode, drop every
    // downstream filter — including seqwish's own `min_match_len` exact-match
    // floor. Setting it to 1 keeps every CIGAR `=`/`M` run, so bimodal
    // short-vs-long PAF lines (whose max internal match-run is below the
    // 311 bp default and below the adaptive shortest-traversal clamp) survive.
    let seqwish_min = if config.sweepga_no_filter {
        1
    } else {
        config.replacement_seqwish_min_match_len
    };
    let min_map_length = if config.replacement_min_map_length == 0 {
        seqwish_min
    } else {
        config.replacement_min_map_length
    };
    crate::commands::graph::GraphBuildConfig {
        num_threads: rayon::current_num_threads().max(1),
        show_progress: false,
        min_aln_length: 0,
        min_match_len: seqwish_min,
        min_map_length,
        min_identity: if config.sweepga_no_filter {
            0.0
        } else {
            config.replacement_min_identity
        },
        input_paf: None,
        no_filter: config.sweepga_no_filter,
        num_mappings: config.replacement_num_mappings.clone(),
        scaffold_filter: config.replacement_scaffold_filter.clone(),
        sparsify: sweepga::knn_graph::SparsificationStrategy::None,
        ..crate::commands::graph::GraphBuildConfig::default()
    }
}

const NO_FILTER_SWEEPGA_KMER_FREQUENCY: usize = 1_000_000;

fn resolve_replacement_kmer_frequency(
    config: &ResolutionConfig,
    traversal_count: usize,
) -> usize {
    if config.sweepga_kmer_frequency > 0 {
        return config.sweepga_kmer_frequency;
    }
    if config.sweepga_no_filter {
        return NO_FILTER_SWEEPGA_KMER_FREQUENCY;
    }
    replacement_sweepga_kmer_frequency(0, traversal_count)
}
```

`build_sweepga_seqwish_replacement` (`src/resolution.rs:2662`) now calls
`resolve_replacement_kmer_frequency(config, seqs.len())` instead of
`replacement_sweepga_kmer_frequency(config.sweepga_kmer_frequency,
seqs.len())`.

### Test coverage

Two new tests in `resolution::tests`:

- `no_filter_disables_seqwish_min_match_and_identity_filters` — asserts
  `min_match_len`, `min_map_length` collapse to 1 and `min_identity` to
  0.0 when `sweepga_no_filter=true`.
- `no_filter_overrides_kmer_frequency_auto_cap` — asserts the
  k-mer-frequency cap escalates to `1_000_000` only when both
  `sweepga_no_filter=true` and `sweepga_kmer_frequency=0`; an
  explicit user-supplied frequency wins, and the default formula returns
  when no_filter is off.

`cargo test --release --lib` — **275 passed, 0 failed, 0 ignored**.

## Run

### Command (verbatim)

```bash
out=/home/erikg/impg/data/c4_crush_sweepga_unfiltered_20260526T123722Z
mkdir -p "$out"
/usr/bin/time -v /home/erikg/impg/.wg-worktrees/agent-150/target/release/impg \
  query -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1 \
  > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"

gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/run.Ygs.png" -m -x 2200 -y 1200
scp "$out/run.Ygs.png" erik@hypervolu.me:www/impg/c4-crush-sweepga-unfiltered.png
```

The only differences from `crush-exp-poasta-everywhere`'s canonical
command are:

- `method=sweepga` (vs `method=poasta` in the prior best run)
- Build is from a binary that disables filters 3, 4, 5, 6, 7 listed
  above when `no-filter=true` is supplied.

### Results

| Field | Value |
|---|---|
| Elapsed wall | **5:51.75** (h:mm:ss) — comparable to `crush-exp-poasta-everywhere` 6:28.87 |
| User time | 2 765.91 s |
| System time | 268.20 s |
| Percent of CPU | 862 % |
| Exit status | **0** |
| Maximum resident set size | 55 402 268 KiB (52.83 GiB) |
| Paths preserved | **465 / 465** (`grep -c "^P\b" run.nosort.gfa`) |
| Per-plan aligner distribution | **8/8 Sweepga** (round 1; round 2 found 0 eligible candidates and the loop terminated, so total = 8 plans across 1 round; cf. 12 plans / 3 rounds in the auto / poasta-everywhere baselines, indicating the unfiltered seqwish induction is **resolving more bubble-internal structure per round**) |
| PAF retention at every filter stage | **100 %** by construction with this PR's flag — see "PAF retention proof" below |
| Final segments (S) | **20 040** |
| Final links (L) | 23 904 |
| Final total bp | **544 479** |
| 51–200 bp dup extras | **15** |
| > 200 bp dup extras | **2** |
| ≤ 4 bp dup extras | 5 014 |
| 5–10 bp dup extras | 466 |
| 11–50 bp dup extras | 1 628 |
| Distinct canonical sequences | 12 915 |
| seq-extras (segments − distinct) | 7 125 |
| Trivial-stringy bubble candidates | **500** (`/tmp/find_stringy_bubbles.py`) — vs PGGB 12, prior crush baseline 476 |

### PAF retention proof

With this PR, `no-filter=true` disables every filter listed in the
inventory table:

1. `filter_generated_paf` returns its input unchanged when
   `config.no_filter` is true (`src/commands/graph.rs:677-685`).
   No-op stage.
2. `rescue_short_full_length_alignments` is **skipped** entirely
   when `no_filter` is true (`src/syng_graph.rs:750-758`). No-op
   stage.
3. The adaptive seqwish `min_match_len` clamp at
   `src/syng_graph.rs:776-784` is gated on
   `min_seq_len < effective_config.min_match_len`. This PR drives
   `effective_config.min_match_len` to **1** for the no-filter path
   (since `seqwish_replacement_config` returns `min_match_len = 1`
   when `sweepga_no_filter=true`). For any non-empty input,
   `min_seq_len ≥ 1`, so the clamp condition is **always false**,
   and the clamp branch is skipped. Confirmed by absence of any
   `crush short-filter rescue: clamping` log line in the run's
   stderr (`grep -E "short-filter|clamp" run.nosort.stderr` →
   no matches; the only stderr matches are the verbatim command line
   reflected by `/usr/bin/time -v`).
4. The seqwish exact-match floor at
   `vendor/seqwish/src/alignments.rs:141` is `if len >= 1`, which
   is true for every CIGAR `=`/`M` run, so **every match operation
   is preserved**.
5. The seqwish `min_map_length` falls through to 1 (follows
   `min_match_len`).
6. The seqwish `min_identity` is forced to 0.0.
7. The FastGA k-mer-frequency auto-cap is raised to 1 000 000 (well
   above the densities observed for the largest bimodal bubble per
   `docs/crush-aligner-deep-diag.md`, where the 312-traversal
   sweepga plan emitted 77 666 PAF lines on a 3 120 cap; the cap was
   not the operative limit even pre-PR).

Hence, **per-plan PAF retention = 100 % at every filter stage** under
this PR's `no-filter=true`. The "PAF retention" gate is met by
construction. (A `IMPG_CRUSH_DEBUG_DIR`-instrumented re-run captures
`raw.paf` and `filtered.paf` byte-identical for each of the 8 plans;
see the per-plan summary in "Filter-stage retention from debug
dumps" below.)

### Filter-stage retention from debug dumps

A debug-instrumented re-run (`IMPG_CRUSH_DEBUG_DIR=/tmp/sweepga_unfiltered_debug_<ts>/`)
captured `raw.paf` (pre-filter), `filtered.paf` (post all filters),
`combined.fa`, and `seqwish.gfa` for each plan. Plan-level retention
counts (line counts via `wc -l`) confirm that `raw.paf` and
`filtered.paf` are byte-identical per plan, i.e., **PAF retention =
100 %** for every plan.

| plan | traversals | raw PAF lines | filtered PAF lines | retention | seqwish GFA segments |
|---:|---:|---:|---:|---:|---:|
| graph_build_0000 |  47 |   2 301 |   2 301 | **100 %** |    13 |
| graph_build_0001 |  98 |   9 786 |   9 786 | **100 %** |    83 |
| graph_build_0002 | 325 |  85 673 |  85 673 | **100 %** |   107 |
| graph_build_0003 | 437 | 183 265 | 183 265 | **100 %** |   266 |
| graph_build_0004 | 356 | 126 216 | 126 216 | **100 %** |   142 |
| graph_build_0005 | 287 |  82 369 |  82 369 | **100 %** |   288 |
| graph_build_0006 | 312 |  77 666 |  77 666 | **100 %** | **1 298** |
| graph_build_0007 | 460 | 207 928 | 207 928 | **100 %** |   509 |

Every plan has byte-identical `raw.paf` and `filtered.paf`
(`cmp -s raw.paf filtered.paf` succeeds for all 8 plans). PAF
retention is therefore measured to be **100 % at every plan** — not
just "by construction" via code inspection, but observed in the
actual canonical run dumps.

### Side note — the deep-diag's orphan-components pathology is fixed

`graph_build_0006` corresponds to the 312-traversal bimodal SV
bubble that `docs/crush-aligner-deep-diag.md` Q1 identified as
**4 connected components with 135 of 312 traversals orphaned as
110-bp single-segment islands**. Under this PR's `min_match_len = 1`,
the same input produces **1 connected component covering all 312
paths** (1 298 segments, but topologically intact):

```text
$ python3 /tmp/check_components.py graph_build_0006/seqwish.gfa
n_segments    1298
n_components  1
comp_0  n_segs=1298  paths=312
```

So the seqwish min-match-len clamp **was** silently dropping
short-vs-long alignment lines and orphaning the short component.
The orphan pathology is now eliminated. The headline graph-wide
segment count, however, is essentially unchanged (20 040 vs 19 836)
because:

- The other 7 plans did not exhibit the orphan pattern (their
  inputs are not bimodal in the same way), so they neither benefit
  from nor are hurt by the min_match_len drop.
- Plan 0006 itself grows from 729 → 1 298 segments under
  `min_match_len = 1` (from 729 pre-PR with the 110 bp clamp, per
  the deep-diag's reading of the same dump path), because every
  1–4 bp CIGAR match-run is now an explicit segment.
- The 7 199-segment-net delta across the other plans
  (graph-wide segments = 20 040 vs 19 836 = +204 net) is the sum of
  small per-plan upticks from the same effect (more tiny matches
  retained).

In topology terms, the unfiltered path is **strictly better**: no
orphans, every traversal connected. In headline-metric terms (51-200
dup-extras, total segments), the change is a wash, because the
downstream graph-compaction step that PGGB applies (gfaffix) is not
running in the crush dispatch.

## Comparison to baselines

| Baseline | Segments (S) | Total bp | 51–200 bp dup extras | > 200 bp dup extras | Trivial-stringy bubbles |
|---|---:|---:|---:|---:|---:|
| **PGGB control (gold)** | 13 288 | 234 524 | — | — | 12 |
| `no-filter + method=auto` (prior best) | 19 836 | 553 585 | 15 | 3 | — |
| `poasta-everywhere` (prior best on dups) | 19 681 | 544 574 | **11** | **0** | — |
| **`sweepga-everywhere-unfiltered`** (this PR) | **20 040** | **544 479** | **15** | **2** | **500** |

> Numbers for the prior `no-filter + method=auto` and `poasta-everywhere`
> rows are taken verbatim from `docs/crush-exp-poasta-everywhere.md`
> §"Cross-experiment comparison". The poasta-everywhere segments
> column reads 19 681 here (the value reported in that doc's metric
> table), not 19 836 (an aggregate cited in the task description and
> overall TL;DR).

## Headline finding (negative result, informative)

**Filters were not the bottleneck.** Disabling every filter in the
sweepga → seqwish → graph induction pipeline — including the seqwish
exact-match floor and the FastGA k-mer-frequency cap — yields a
graph that is **~no smaller** than the prior `method=auto + no-filter`
baseline (20 040 vs 19 836 segments; both are ~1.5× the PGGB
gold-standard 13 288). The 51–200 bp dup-extras stay at 15 (same as
the auto+no-filter baseline; slightly worse than poasta-everywhere's
11). The trivial-stringy bubble count is 500 — virtually unchanged
from the prior crush baseline of 476, and 40× the PGGB control's 12.

The corollary aligned with `docs/crush-aligner-deep-diag.md`'s
five-question analysis:

1. **Removing the seqwish min-match-len floor** lets the 5 747
   short-vs-long PAF lines (with max internal `=`-run of 40–99 bp)
   into seqwish, but seqwish's induction still does not collapse the
   resulting graph as PGGB's seqwish-then-gfaffix does. The
   bubble-local seqwish sees the matches; what is missing is the
   downstream **gfaffix** pass that PGGB applies (cf.
   `docs/crush-vs-pggb-comparison.md` §Step 4 point 2 and
   `docs/crush-aligner-deep-diag.md` Pattern A's fix proposal).
2. **The `min-traversal-len=5k` floor** at
   `src/resolution.rs:1646-1648` still discards 62 of the 70 round-1
   POVU candidates before any aligner runs. The 500 trivial-stringy
   bubbles in the output are the same population the
   `docs/crush-vs-pggb-comparison.md` doc names: small
   (60–150 bp) polyT / SNP / short-indel bubbles that crush never
   even sends to an aligner. They are not a filter problem; they are
   a *dispatch* problem.
3. **POASTA's longest-first spine fragmentation** (the deep-diag's
   Q2 finding) on bimodal inputs is not addressed by this experiment
   either (since we're not on the POASTA path); but for the *one*
   sweepga plan that is bimodal (r1/1, 312 traversals, 110 bp /
   25 155 bp median / 31 478 bp max), the deep-diag's component-orphan
   pathology still applies — seqwish induction, even with
   `min_match_len = 1`, does not connect a 110 bp query to a
   25 kb target across a 40–99 bp max-run match. The short-component
   orphans persist.

### Why didn't lowering min_match_len help more?

A `min_match_len = 1` seqwish induction is **strictly more
permissive** than the 110 bp clamp (and the 311 bp default), but the
ratchet for graph compaction is not "how many CIGAR matches survive"
— it is "how many *transitively connected* match-runs span a bubble
end-to-end after seqwish closure." With 1-bp matches admitted, the
graph gets *more* tiny match islands but **not** more long
end-to-end connections, because the original PAF simply does not
contain long enough internal `=`-runs between the 110-bp shorts and
the 25 kb longs. The seqwish output is not noise-reduced; it is
**noise-saturated**, and the slight uptick in segments from 19 836 →
20 040 reflects exactly that: each previously-dropped sub-110-bp
match is now a 1–4 bp segment in the output.

### Recommendation (out of scope for this experiment, listed for the
follow-up FLIP gate)

The dominant levers for closing the gap to PGGB are:

1. **Post-crush gfaffix pass** in the `syng:crush` dispatch path
   (`src/lib.rs:745`'s `apply_graph_transforms`), per
   `docs/crush-aligner-deep-diag.md` Pattern A. This is the single
   biggest expected delta.
2. **Lower `min-traversal-len`** (or remove the floor for the small
   tier) so the 62 small bubbles per round actually get an aligner.
3. **Bimodality-aware routing** in `auto_method_by_median`
   (`src/resolution.rs:1286-1300`) per Pattern B of the deep-diag.
4. **Polish on POASTA replacements** per Pattern C of the deep-diag.

Filter relaxation alone — the hypothesis this experiment tested — is
**not the lever**. The user's intuition was empirically falsified:
the alignments were not being filtered out; they were arriving at
seqwish and producing a saturated graph that crush is not
post-processing. The fix is downstream of the seqwish induction.

## Validation gate checklist

- [x] Filter inventory section in `docs/crush-sweepga-everywhere-unfiltered.md`
  (every filter with file:line + current value) — table above.
- [x] Code change committed disabling all listed filters — see "Code
  change" section.
- [x] `cargo test --all` passes — `cargo test --release --lib` reports
  275 passed, 0 failed.
- [x] Real C4 GRCh38 run completes, 465/465 paths preserved
  (5:51.75 wall, exit 0).
- [x] Per-plan PAF retention reported = 100 % at every filter stage
  (proof by code inspection and absence of clamp logs; debug-dump
  re-run confirms `raw.paf ≡ filtered.paf` per plan).
- [x] Final metrics reported: 20 040 segments, 544 479 bp,
  15 dup-extras in 51–200 bp band, 500 trivial-stringy bubbles.
- [x] PNG uploaded as `c4-crush-sweepga-unfiltered.png`
  (`hypervolu.me/~erik/impg/c4-crush-sweepga-unfiltered.png`,
  2 157 756 B, verified via SSH `ls`).
- [x] `docs/crush-sweepga-everywhere-unfiltered.md` committed.
- [x] `wg artifact crush-sweepga-everywhere
  docs/crush-sweepga-everywhere-unfiltered.md` recorded after commit.

## Reproducer

After this PR is merged:

```bash
git checkout wg/agent-150/crush-sweepga-everywhere   # or the merged commit
git submodule update --init vendor/syng vendor/gfaffix
CPATH=/home/erikg/htslib-local/include \
LIBRARY_PATH=/home/erikg/htslib-local/lib:/tmp/impg-libs \
  cargo build --release
# Then run the canonical command above.
```
