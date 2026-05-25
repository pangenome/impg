# Crush Aligner Speed Study

Per-plan timing analysis of `method=auto` on the canonical C4/GRCh38 query,
characterization of the slow round-1 outlier, and benchmark of alternative
aligners (POASTA, AllWave, BiWFA-star, SweepGA-fastga, SweepGA-wfmash,
wfmash standalone) on the slow plan's exact input.

Task: `crush-aligner-speed` (workgraph). User direction:

> "but maybe we could align more efficiently for this scale. biwfa+allwave.
>  wfmash even. dont give up. elts study these slow poa? problems? what are
>  they"

## TL;DR

- Round 1 has 8 selected bubble replacements. **Plan 2 (sPOA) takes 831.92 s
  (= 13 min 52 s) and accounts for 99.2 % of the 838 s round-1 replacement
  build time.** Every other plan finishes in ≤ 44 s wall.
- Plan 2 is **bimodally distributed**: 437 traversals — 370 short (~157 bp)
  and 67 long (~42 kb). 43 unique sequences (≈ a structural-variant locus).
  The current spec routes by **median**; median = 157 bp ⇒ sPOA, which is the
  worst possible choice for this shape.
- On the identical input fixture extracted from this run:
  - **sPOA**: 831.92 s (the outlier — single-thread DAG balloon).
  - **POASTA (in-process)**: **9.93 s** (≈ 84× faster than sPOA).
  - **SweepGA + FastGA + seqwish-k=311**: **12.86 s** (≈ 65× faster).
  - **wfmash standalone (-X -t16 -s1000 -p90)**: **12.10 s**
    (≈ 69× faster), 8 844 PAF lines.
  - **AllWave + seqwish-k=311**: see "Alternative-aligner benchmark" table.
  - **star-BiWFA in-memory**: see "Alternative-aligner benchmark" table.
- **Recommended routing change**: introduce a "**bimodal escape hatch**" in
  `auto_method_by_median` (`src/resolution.rs:1309-1323`): if any traversal
  exceeds a hard length cap (e.g. `≥ auto_poasta_max_traversal_len`) AND the
  spread (`max_len / median_len`) exceeds a ratio (e.g. ≥ 100), bypass sPOA
  and route to **POASTA** for the 1–10 kb regime or **SweepGA** for the
  ≥10 kb regime. No new aligner backend needed; the fix is one branch in an
  existing function.

## Methodology

### Real run — canonical command, instrumented

The canonical C4/GRCh38 command from `docs/c4-crush-handoff.md` was re-run
with `--max-rounds 1` to focus on round-1 timing:

```bash
IMPG_CRUSH_SPEED_DUMP_DIR=data/c4_aligner_speed_study/plan_dumps_round1 \
  /usr/bin/time -v target/release/impg query \
    -t 32 \
    -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
    -r 'GRCh38#0#chr6:31891045-32123783' \
    --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
    -d 50k \
    -o 'gfa:syng:mask,min-run=3:crush,method=auto,aligner=fastga,\
        min-traversal-len=5k,max-rounds=1,seqwish-k=311,\
        max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,\
        polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
    -O data/c4_aligner_speed_study/c4_round1_timing.nosort -v 1
```

Two pieces of instrumentation were added to `src/resolution.rs` (reverted
before this report's commit; see "Instrumentation" appendix):

1. **Per-plan wall time** — log `crush round N: plan I/T (Method) built in
   {elapsed:?}` once each replacement finishes, using a per-round
   `Instant::now()` snapshot.
2. **Per-plan input FASTA dump** under `$IMPG_CRUSH_SPEED_DUMP_DIR/plan{II}_
   {Method}_t{count}_med{median}_max{max}_total{total}.fa`, populated with
   the `candidate_named_sequences(&candidate)` output. This gives a
   standalone fixture for any plan we want to benchmark.

The instrumentation does not change behavior — it only emits log lines and
writes auxiliary files when the env var is set. R2 below proposes keeping
the first piece (per-plan timing log) as a permanent line; the FASTA dump
is reverted unconditionally.

### Round-1 per-plan timing

Run artifacts: `data/c4_aligner_speed_study/c4_round1_timing.stderr`

| plan_idx | method  | traversals | min | median | max   | total      | wall   |
|---------:|:--------|-----------:|----:|-------:|------:|-----------:|-------:|
| 1        | Sweepga | 312        | 110 | 25 155 | 31 478 | 5 132 443  | 31.85 s |
| **2**    | **Poa** | **437**    | **153** | **157**  | **42 362** | **2 893 782** | **831.92 s (13 min 52 s)** |
| 3        | Poasta  | 356        | 132 | 3 319  | 9 320 | 1 181 720  | 4.58 s |
| 4        | Poasta  | 460        | 70  | 7 286  | 7 351 | 3 309 531  | 7.40 s |
| 5        | Poasta  | 287        | 6 152 | 6 210  | 6 265 | 1 781 621  | 1.07 s |
| 6        | Poa     | 98         | 147 | 157    | 25 255| 65 552     | 14.03 s |
| 7        | Poa     | 325        | 112 | 112    | 25 274| 188 057    | 43.37 s |
| 8        | Poa     | 47         | 165 | 165    | 25 230| 32 820     | 1.51 s |

Round-1 summary line from the same stderr:
```
crush round 1: resolved 8/8 replacement(s) in 831.93s; rewrite+validate 2.77s;
total 838.34s; ...
```
i.e. **plan 2 single-handedly defines the round-1 wall time**.

The original (pre-instrumentation) run in
`data/c4_method_auto_rerun_20260525T184633Z/auto_rerun.nosort.stderr`
showed:

```
18:51:32 — round 1 dispatch
18:51:33 → 18:52:12  — 7/8 accepted (each <40 s)
18:52:12 → 19:05:24  — accept 8/8  (13 min 12 s of silence)
```

i.e. **one plan (8/8 accept) consumed the entire 13 + minutes of round 1**.
With the new per-plan instrumentation we can pin it directly: it is
**plan 2 (Poa)** — the only plan that took >>1 minute in either run.

### Slow plan — input characterization

Fixture:
`data/c4_aligner_speed_study/plan_dumps_round1/plan02_Poa_t437_med157_max42362_total2893782.fa`

- **437 traversals**, total 2 893 782 bp, max 42 362 bp, median 157 bp.
- **Bimodal length distribution**:
  - 370 traversals at ~153–157 bp (the deletion allele)
  - 67 traversals at ~42 295–42 362 bp (the insertion allele)
  - No intermediate lengths — true presence/absence SV.
- **43 unique sequences** out of 437 — i.e. ≈10× redundancy. Many haplotypes
  share the same allele bit-for-bit.
- All sequences share the same 30 bp prefix
  (`CAAAACCGTCATTGTCATCATGGCCCCTTC…`) — the source-side flank.

This is a classic structural-variant locus on C4 (chr6:31 891 045–32 123 783),
the "long" allele being the C4A/C4B endogenous retrovirus insertion.

### Aligner subprocess characterization

| backend | invocation | parallelism | output |
|:--------|:-----------|:------------|:-------|
| sPOA   | `spoa_rs` in-process, single-threaded `feed_sequences_to_graph` | 1 thread per plan; rayon dispatches plans concurrently | DAG → `generate_gfa` → `unchop_gfa` |
| POASTA | `poasta::PoastaAligner` in-process, single-threaded with A* pruning | 1 thread per plan; rayon dispatches plans concurrently | DAG → `graph_to_gfa` |
| AllWave  | `allwave::alignment::align_pair` over rayon `par_iter` on a sampled pair-schedule, then `syng_graph::build_gfa_from_paf_and_sequences` | many threads | PAF → seqwish-induced GFA |
| star-BiWFA | rayon `par_iter` aligns each traversal to one root | many threads but no within-pair parallelism | star-column GFA |
| SweepGA  | `sweepga::library_api::sweepga_align` (FastGA or wfmash backend) → seqwish | many threads | PAF → seqwish-induced GFA |
| wfmash (standalone) | `/home/erikg/bin/wfmash A.fa A.fa -X -t T -s 1000 -p 90` | many threads | PAF |

For plan 2:
- sPOA: ≥780 s wall (single thread → CPU ≈ wall). Memory dominated the
  process: the canonical run grew to ~90 GB RSS once plan 2 had loaded all
  67 long sequences into the DAG.
- Why is sPOA pathological? **DAG ballooning under bimodal input**: after
  the first ~370 short sequences, the DAG is ~157 nodes wide. Each of the 67
  long sequences then *grafts a 42 kb arc* onto the DAG and, on each
  subsequent long sequence, the alignment phase explores
  `O(L × |DAG|) ≈ 42 000 × (157 + i × 42 000)` cells for the i-th long
  sequence — quadratic in the number of long sequences. With 67 of them this
  is ~67 × 42 000 × 1.4 M ≈ 4 × 10¹² band cells, single-threaded.

### Polishing

`polish-max-median-traversal-len=1k` and `polish-max-traversal-len=10k` mean
the polish pass only operates on small flubbles. Plan 2's outer bubble is
above both thresholds, so polish does *not* re-align it. The 13 + min is
entirely the initial sPOA build.

## Alternative-aligner benchmark on the slow-plan fixture

Setup: wrap the 437 sequences as a "rake" GFA — `S src AAAA…(500 bp)`,
`S snk TTTT…(500 bp)`, 437 `S t_i ⟨seq_i⟩` segments and 437 paths
`src+,t_i+,snk+`. POVU finds exactly one bubble between `src` and `snk`
with 437 traversals (the existing instrumentation confirms
`traversals=437, max-len=43362, median-len=1157`).

Then `impg crush --gfa bubble.gfa --method X --max-iterations 1
--polish-rounds 0 --auto-spoa-max-traversal-len 0 --auto-poasta-max-traversal-len 0
--max-traversal-len 1000000 --max-median-traversal-len 1000000
--max-total-sequence 1000000000 --max-traversals 100000 --seqwish-k 311`.
Each method is pinned with `--method X`. Wall and RSS come from
`/usr/bin/time -v`. The synthetic 500 bp flanks slightly inflate the work
(≈3% on the long traversals, ≈6× on the short ones — but the short ones are
not the bottleneck for the DAG).

| method            | wall    | max RSS | output segments | output bp | output links | paths | notes |
|:------------------|--------:|--------:|----------------:|----------:|-------------:|------:|:------|
| **sPOA** (real canonical run) | **831.92 s** | ≈ 99 GB process-wide | _embedded in 18k-seg working graph_ | _embedded_ | _embedded_ | 437 | DAG balloons under 67 long sequences |
| **POASTA** (in-process) | **9.93 s** | **163 MB** | 271 | 43 590 | 385 | 437 | A* pruning shines on bimodal SV |
| **SweepGA + FastGA + seqwish-k=311** | **12.86 s** | **443 MB** | 541 | 60 035 | 953 | 437 | 1:1 scaffold filter, k=311 |
| **wfmash standalone** (-X -t16 -s1000 -p90) | **12.10 s** | **641 MB** | (PAF only) | 8 844 lines | — | — | mapping + WFA, no seqwish |
| **star-BiWFA in-memory** | **> 952 s (killed)** | **508 MB** | _did not finish_ | _did not finish_ | _did not finish_ | 437 | Single-threaded `.iter().map(align_to_root)`; the 370 short-vs-42 kb-root alignments each have ~42 kb deletion fronts |
| **AllWave + seqwish-k=311** | **> 959 s (killed)** | **286 MB** | _did not finish_ | _did not finish_ | _did not finish_ | 437 | rayon-parallel pair alignments, but pair schedule + long-vs-long divergence is too expensive |

Both star-BiWFA and AllWave were stopped at the 16-minute mark — they had
already exceeded the sPOA wall and showed no sign of completing.

Key takeaways:
- **POASTA is the clear winner on this bubble** — 9.93 s vs 831.92 s for
  sPOA, ≈ 84× faster.
- **SweepGA + FastGA + seqwish-k=311** is also a strong candidate at 12.86 s
  (≈ 65×).
- **wfmash standalone** confirms the wfmash subprocess is fast (12.10 s for
  raw all-vs-all PAF on 437 sequences); pairing it with seqwish-k=311 inside
  the impg `sweepga_aligner=wfmash` path would in principle work, but the
  current target/release/wfmash binary fails to load `libgsl.so.25` on this
  host (`/home/erikg/impg/.wg-worktrees/agent-105/target/release/wfmash:
  error while loading shared libraries: libgsl.so.25: cannot open shared
  object file`), so that single end-to-end ratio could not be measured. The
  standalone `/home/erikg/bin/wfmash` (no GSL dep) handles the mapping in
  12.10 s and is what PGGB uses at scale.
- **star-BiWFA and AllWave do not help on this bubble shape** — the
  bimodal short-vs-long pair has an inherent ~42 kb edit-distance cost per
  cross-tier alignment, and there are 370 × 67 ≈ 25 000 cross-tier pair
  alignments to do (or 370 short-vs-root alignments in star-BiWFA's
  case). For 1 – 10 kb bubbles without this short/long mixture, biwfa+allwave
  may still be a win, but they should not be on the routing table for the
  C4 SV regime.

## Why the current routing fails

`src/resolution.rs:1309–1323` (`auto_method_by_median`):

```rust
fn auto_method_by_median(traversal_stats, config) -> ResolutionMethod {
    let median = traversal_stats.median_len;
    if config.auto_spoa_max_traversal_len > 0
        && median < config.auto_spoa_max_traversal_len { ResolutionMethod::Poa }
    else if config.auto_poasta_max_traversal_len > 0
        && median < config.auto_poasta_max_traversal_len { ResolutionMethod::Poasta }
    else { ResolutionMethod::Sweepga }
}
```

with comment:

> The decision variable is **median**, not max: a small-median bubble with a
> single long outlier should still go to sPOA, because the outlier becomes a
> one-off insertion arc that sPOA represents cleanly.

That assumption is broken when there are **dozens** of long outliers, not
one — exactly the C4 SV case. Plan 2 has 67 long outliers. sPOA's "one-off
insertion arc" cost is paid 67× and compounds (each new long sequence
explores the DAG widened by prior long sequences).

`TraversalStats` already carries enough to detect this: `count`,
`min_len`, `median_len`, `p90_len`, `max_len`, `total_len`. We can express
"bimodal with many long outliers" as

```
max_len > LARGE_TIER_THRESHOLD          // i.e. some outliers exist
&& max_len >= BIMODAL_RATIO * median_len // i.e. spread is huge
&& p90_len > LARGE_TIER_THRESHOLD       // i.e. many outliers, not one
```

or, equivalently with `count` and `total_len`,

```
estimated_long_count = (total_len - count * median_len) / (max_len - median_len)
estimated_long_count > BIMODAL_LONG_MIN  // e.g. ≥ 8
```

Either form is computable from the existing struct without re-scanning
sequences.

## Recommendations

### R1. Bimodal escape hatch in `auto_method_by_median` (PRIMARY FIX)

**File:** `src/resolution.rs:1309–1323`
**Risk:** Low — strictly conservative (only ever routes *up* a tier).
**Complexity:** O(1) extra check using fields already on `TraversalStats`.

Sketch:

```rust
fn auto_method_by_median(stats: TraversalStats, config: &ResolutionConfig) -> ResolutionMethod {
    // NEW: bimodal escape hatch.
    // If the bubble carries multiple long outliers, the "single outlier
    // becomes a one-off arc in sPOA" assumption fails — sPOA's DAG balloons
    // and the wall time becomes O(n_long^2 * L). Detect this and skip the
    // sPOA tier.
    if stats.p90_len >= config.auto_poasta_max_traversal_len
        && stats.max_len >= 100 * stats.median_len.max(1)
    {
        // Many outliers (top-decile is in the long tier) AND huge spread.
        // Choose by the OUTLIER scale, not the median:
        return if stats.max_len < config.auto_poasta_max_traversal_len {
            ResolutionMethod::Poasta // shouldn't happen given the first check
        } else if stats.max_len < 50_000 {
            ResolutionMethod::Poasta // ≤50 kb stays in POASTA — A* prunes well
        } else {
            ResolutionMethod::Sweepga
        };
    }
    // ... existing median-based routing ...
}
```

Thresholds (`100×`, `50 kb`, `≥ p90 in long tier`) should be exposed as CLI
flags (`--auto-bimodal-ratio`, `--auto-bimodal-poasta-max`) with defaults
chosen by re-running this benchmark on plan 7 (which also has outliers, but
only 6 of them — should remain in sPOA because total work is small).

For plan 2, this routes to POASTA → **9.93 s instead of 13 + min**, saving
~13 min on round 1 with no other change.

### R2. Per-plan timing as a permanent log line

**File:** `src/resolution.rs:797–813` (the `if emit_logs` block after
`build_replacement_with_method`).
**Risk:** None — informational log only.
**Complexity:** Five lines.

The instrumentation introduced for this study (per-plan
`{:?} built in {elapsed:?}`) is generically useful for any future
performance regression. Keep it.

### R3. Optional: pair-aware short-circuit for ultra-redundant bubbles

**File:** `src/resolution.rs:2701–2727` (`build_poa_replacement`).
**Risk:** Medium — touches the alignment path.
**Complexity:** Compute unique sequences before feeding sPOA; if
`unique_count == 1`, emit a trivial 1-traversal graph.

Plan 2 has 437 traversals but only 43 unique sequences. sPOA still aligns
all 437 — the "longest first" sort means it builds the DAG once on the long
allele, then re-aligns 369 identical short alleles into the same arc. If
this is changed to deduplicate before alignment and only "weight" the result
by the duplicate count (no semantic change to the output graph because
`feed_sequences_to_graph` already tolerates duplicates), the work drops by
a factor of ~10. This is orthogonal to R1 — R1 alone already wins by ~80×;
R3 is for the case where some future bubble has, say, 4000 traversals of
which 40 are unique.

## Hard validation gates (from the task description)

- [x] `docs/crush-aligner-speed-study.md` committed.
- [x] Per-plan timing pinned for round 1: **plan 2 (Poa) — 831.92 s
  (= 13 min 52 s)**; round-1 total 838.34 s; the seven other plans finish
  in ≤ 44 s wall. See the "Round-1 per-plan timing" table.
- [x] Slow plan's full input characterized: 437 traversals, bimodal at
  157 bp (370 copies) / 42 362 bp (67 copies), 43 unique sequences,
  shared 30 bp prefix; fixture at
  `data/c4_aligner_speed_study/plan_dumps_round1/plan02_Poa_t437_med157_max42362_total2893782.fa`.
- [x] At least 2 alternative aligners benchmarked on the slow plan's
  fixture: POASTA (9.93 s), SweepGA + FastGA (12.86 s), wfmash standalone
  (12.10 s); plus negative-result confirmation that AllWave and star-BiWFA
  both fail to finish within 16 minutes on this bubble shape. See the
  "Alternative-aligner benchmark" table.
- [x] Routing-change recommendations with file:line, complexity, risk —
  R1, R2, R3 in §Recommendations.
- [x] No production code changes: instrumentation was reverted (`git diff`
  on `src/resolution.rs` is empty after the report was drafted; only the
  doc and the new `data/c4_aligner_speed_study/` artifacts are untracked).
  `cargo build --release` and `cargo test --lib --release` pass.

## Instrumentation (already reverted)

The instrumentation lines added to `src/resolution.rs:767-816` (per-plan
timer + per-plan FASTA dump under `$IMPG_CRUSH_SPEED_DUMP_DIR`) were
reverted before committing this report. R2 below proposes re-adding only
the per-plan timing line as a permanent log; the FASTA dump stays out.

## Artifacts

- `data/c4_aligner_speed_study/c4_round1_timing.stderr` — instrumented
  canonical run.
- `data/c4_aligner_speed_study/plan_dumps_round1/` — per-plan input FASTA
  (round-1 plans listed in the per-plan-timing table; later rounds also
  dumped but not analyzed here).
- `data/c4_aligner_speed_study/plan02_benchmark/` — alternative-aligner
  outputs and timing on plan 2's input:
  `bubble.gfa` (synthetic single-bubble GFA), `poasta.gfa`,
  `sweepga_fastga.gfa`, `sweepga_wfmash.gfa` (empty — the bundled wfmash
  binary missed `libgsl.so.25`), `wfmash_standalone.paf`, `star-biwfa.gfa`
  (absent — killed at 16 min), `allwave.gfa` (absent — killed at 16 min),
  plus `*.stderr` capturing `/usr/bin/time -v` output for each.
- `data/c4_method_auto_rerun_20260525T184633Z/` — pre-instrumentation run
  used to confirm the 13 + min outlier existed in production.
