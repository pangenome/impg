# Crush experiment: SweepGA many-to-many + systematic POASTA polish

**Task:** `crush-sweepga-many-to-many-poasta-polish` (workgraph id `crush-sweepga-many`)
**Date:** 2026-05-26
**Branch:** `wg/agent-163/crush-sweepga-many`
**Output dir:** `/home/erikg/impg/data/c4_crush_many_to_many_poasta_polish_20260526T185552Z/`
**PNG:** `https://hypervolu.me/~erik/impg/c4-crush-many-to-many-poasta-polish.png`

## User direction

> we should use many-to-many filtering. No scaffold filtering, no plane
> sweep filtering, just like go straight for it. Take whatever FastGA
> gives us. Then that gets crushed into a graph. Plus a systematic
> application of POASTA afterward.

## Hypothesis

The prior `crush-sweepga-everywhere-unfiltered` (`crush-sweepga-everywhere`,
agent-150) PR already disabled every length / identity / k-mer-frequency
filter in the SweepGA → seqwish → graph induction path when
`no-filter=true` is set. After that PR, the per-plan PAF retention is
100 % at every filter stage. The only remaining candidate "bottleneck"
inside SweepGA that the user direction implies could still be active is
the **mapping / scaffold filter mode** — when `no_filter=false` the
default modes are `"1:1"` (one-to-one); the user wants `many:many`.

The other half of the direction is to **polish with POASTA**, not sPOA.
The polish loop already runs at every iteration of the level-descent —
it just dispatches by `polish_method`, which defaults to `Poa` (sPOA),
not `Poasta`. Setting `polish-method=poasta` routes every sub-bubble
inside what SweepGA produced through POASTA.

## Filter inventory (every knob, file:line, current value)

This is the complete set of length / identity / mode / mass knobs that
gate alignment retention in the SweepGA replacement path. Most are
already disabled in the **prior** `crush-sweepga-everywhere-unfiltered`
PR when `no-filter=true` is set; the two columns below show whether
each knob is disabled in the canonical `no-filter=true` invocation.

| # | Knob | Where it lives (file:line) | Current default | Effect when `no-filter=true` (this PR) |
|---|---|---|---|---|
| 1 | **SweepGA library `no_filter`** | `src/resolution.rs:2796` (set to `true` hard-coded) | `true` (impg always passes `true`) | Skips SweepGA's internal post-alignment plane-sweep / scaffold filter. Unchanged by this PR. |
| 2 | **SweepGA library `num_mappings` (filter mode)** | `src/resolution.rs:2797` (this PR pins `"many:many"`) | `"many:many"` (sweepga lib default) | Pinned to `"many:many"` so the documented intent survives if the filter is ever re-enabled. **New in this PR.** |
| 3 | **SweepGA library `scaffold_filter` (mode)** | `src/resolution.rs:2798` | `"many:many"` (this PR) | Pinned to `"many:many"`. **New in this PR.** |
| 4 | **SweepGA library `scaffold_mass`** | `src/resolution.rs:2799` | `0` (this PR) | Pinned to `0`. **New in this PR.** |
| 5 | **GraphBuildConfig `no_filter`** | `src/resolution.rs:2495` (`config.sweepga_no_filter`) | `false` | When `true`, `filter_generated_paf` (`src/commands/graph.rs:677`) early-returns its input unchanged. |
| 6 | **GraphBuildConfig `num_mappings` (filter mode)** | `src/resolution.rs:2496`; sourced from `replacement_num_mappings` default `"1:1"` (`src/resolution.rs:300`) | `"1:1"` | This PR forces this to `"many:many"` in `seqwish_replacement_config` when `sweepga_no_filter=true` (the `filter_generated_paf` step is also skipped, so this is defensive coupling). |
| 7 | **GraphBuildConfig `scaffold_filter` (mode)** | `src/resolution.rs:2497`; sourced from `replacement_scaffold_filter` default `"1:1"` (`src/resolution.rs:301`) | `"1:1"` | This PR forces `"many:many"` when `sweepga_no_filter=true`. |
| 8 | **GraphBuildConfig `scaffold_mass` (min scaffold chain length)** | `src/resolution.rs:2498`; sourced from `replacement_scaffold_mass` default `0` (`src/resolution.rs:264`) | `0` | Already 0 (set by `crush-scaffold-mass` / agent-160). This PR forces it to 0 when `sweepga_no_filter=true` for defensive coupling. |
| 9 | **Short-full-length rescue augmentation** (`rescue_short_full_length_alignments`) | `src/syng_graph.rs:608-692`; called at `:753`, gated by `if !config.no_filter` at `:750-758` | enabled when `no_filter=false` | Skipped when `no_filter=true`. |
| 10 | **Adaptive seqwish `min_match_len` clamp** | `src/syng_graph.rs:763-784` | clamps `min_match_len` down to `min_seq_len` | Implicitly skipped when `no_filter=true` because `seqwish_min` is forced to `1` (`src/resolution.rs:2455-2459`), so `min_seq_len < 1` is false for any non-empty input. |
| 11 | **Seqwish `min_match_len` floor** | `vendor/seqwish/src/alignments.rs:141`; sourced from `replacement_seqwish_min_match_len` default `311` (`src/resolution.rs:258`) | `311` | Forced to `1` when `sweepga_no_filter=true`. |
| 12 | **Seqwish `min_map_length` (per-PAF-line min)** | `src/commands/graph.rs:88`; sourced from `replacement_min_map_length` default `0` (`src/resolution.rs:259`) | `0` → follows `min_match_len` | Falls through to `seqwish_min=1` when `no_filter=true`. |
| 13 | **Seqwish `min_identity`** | `src/resolution.rs:2473`; sourced from `replacement_min_identity` default `0.0` (`src/resolution.rs:260`) | `0.0` | Forced to `0.0` when `no_filter=true`. |
| 14 | **FastGA `kmer_frequency` cap (`replacement_sweepga_kmer_frequency`)** | `src/resolution.rs:2486-2511`; default `max(1000, traversal_count*10)` for the auto path | `max(1000, n*10)` (≈3120 for the 312-traversal bubble) | Raised to `NO_FILTER_SWEEPGA_KMER_FREQUENCY = 1_000_000` when `sweepga_no_filter=true` and `sweepga_kmer_frequency=0` (the unconfigured default). |
| 15 | **`min-traversal-len` floor** (gate before any aligner runs) | `src/resolution.rs:1646-1648` | Independent of `no-filter`. Canonical command uses `min-traversal-len=5k`, which discards 62 small bubbles per `docs/crush-aligner-deep-diag.md` Pattern A | Independent of `no-filter`; kept at `5k` for apples-to-apples comparison. |
| 16 | **`max_pair_alignments` cap** | `src/resolution.rs:139` + flag parser; canonical uses `max-pair-alignments=0` | `10_000` (CLI default) | Canonical command sets `max-pair-alignments=0` (disabled). |
| 17 | **`max_replacement_paf_bytes` cap** | `src/resolution.rs:151` + flag parser; canonical uses `max-paf-bytes=0` | `64 MiB` (CLI default) | Canonical command sets `max-paf-bytes=0` (disabled). |
| 18 | **`min_aln_length` (per-line min block length)** | `src/resolution.rs:2772`; passes `min_aln_length: 0` always | `0` | Always 0 — off. |
| 19 | **Polish-loop tier limits** | `src/resolution.rs:248-253`: `polish_max_traversal_len=10k`, `polish_max_median_traversal_len=1k` | per defaults | Sub-bubbles within these limits are polished by POASTA when `polish-method=poasta`. |

**Summary.** Before this PR, `no-filter=true` already disabled #1, #5,
#9, #10, #11, #12, #13, #14. The mapping / scaffold filter modes (#2,
#3, #6, #7) and scaffold mass (#4, #8) were *unused-but-not-pinned*
(`filter_generated_paf` early-returned without consulting them). This
PR pins them all to `many:many` / `0` so the documented "fully
unfiltered + many-to-many" semantics survive future refactors.

## Code change

`src/resolution.rs` — two changes:

1. **`build_sweepga_seqwish_replacement`** (`src/resolution.rs:2796-2802`):
   the `SweepgaAlignConfig` now passes `num_mappings: "many:many"`,
   `scaffold_filter: "many:many"`, and `scaffold_mass: 0` alongside the
   already-present `no_filter: true`. SweepGA's library default is
   already `many:many`, so this is a no-op at the current sweepga rev
   (`ddd31d3`) — it is defensive coupling for any future change that
   re-enables sweepga's library filter.

2. **`seqwish_replacement_config`** (`src/resolution.rs:2465-2483`):
   when `sweepga_no_filter=true`, the `num_mappings`,
   `scaffold_filter`, and `scaffold_mass` fields in the
   `GraphBuildConfig` are now forced to `"many:many"`, `"many:many"`,
   and `0` respectively. `filter_generated_paf` already early-returns
   on `no_filter=true`, so this is also defensive coupling.

The log line in `build_sweepga_seqwish_replacement`
(`src/resolution.rs:2826-2842`) is updated to print the actual filter
state instead of the stale "scaffold filter=1:1" string.

## Validation

- `cargo build --release` succeeds.
- `cargo test --lib` 281 passed, 0 failed.
- `cargo test --bin impg` 56 passed, 0 failed.
- `cargo test --test test_crush_integration` has one pre-existing
  failure (`nested_bubble_level_descent_actually_descends`) that
  reproduces identically on the pre-change tree
  (verified via `git stash`); it is not introduced by this PR. The
  test failure is annotated in its `panic` message as a known issue at
  `HEAD 471f089`.

## C4 GRCh38 run

### Command

```bash
out=/home/erikg/impg/data/c4_crush_many_to_many_poasta_polish_20260526T185552Z
/usr/bin/time -v target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-method=poasta,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1
```

Only differences from the prior `crush-sweepga-everywhere-unfiltered`
canonical command:

- `aligner=fastga` (was implicit; sweepga default is `wfmash`, but
  `crush-sweepga-everywhere`'s docs used `fastga`)
- **`polish-method=poasta`** (was implicit `poa`)

### Results

| Field | Value |
|---|---|
| Elapsed wall | **5:24.84** |
| User time | 2 736.69 s |
| System time | 390.48 s |
| Percent of CPU | 962 % |
| Exit status | **0** |
| Maximum resident set size | 55 381 148 KiB (52.81 GiB) |
| Paths preserved | **465 / 465** (`grep -c "^P\b" run.nosort.gfa`) |
| Per-plan aligner distribution | **8/8 Sweepga** (round 1; round 2 found 0 eligible candidates and the loop terminated) — same as `crush-sweepga-everywhere-unfiltered` |
| Per-plan PAF retention | **100 %** at every filter stage (by construction — `no-filter=true` skips every #5–#14 in the inventory table, and #1–#4 are also pinned `many:many`/`0` by this PR) |
| Final segments (S) | **19 958** |
| Final links (L) | 23 822 |
| Final total bp | **544 641** |
| Total path steps | 4 522 594 |
| Trivial-stringy bubble candidates | **499** (`/tmp/find_stringy_bubbles.py` heuristic) — vs PGGB 12, prior `crush-sweepga-everywhere-unfiltered` 500, prior `crush-scaffold-mass-zero` 500 |

### Comparison vs prior baselines

| Config | S | L | bp | Polish | Filter |
|---|---|---|---|---|---|
| **This PR** (sweepga, no-filter=true, **polish=poasta**) | 19 958 | 23 822 | 544 641 | POASTA | none |
| `crush-sweepga-everywhere-unfiltered` (agent-150) | 20 040 | 23 904 | 544 479 | sPOA | none |
| `crush-scaffold-mass-zero` (agent-160) | 20 040 | 23 904 | 544 479 | sPOA | scaffold-mass=0 |
| `crush-hierarchical-sweepga` (agent-152) | 20 326 | 24 028 | — | sPOA | inherited |
| `crush-exp-hybrid` (auto + sweepga top-tier) | 19 843 | 23 391 | — | sPOA | none |
| **PGGB gold standard** | 13 288 | — | ~234 000 | (PGGB pipeline) | — |

Switching the polish from sPOA → POASTA reduced segments by **82** (0.4 %).
Pinning the many-to-many semantics inside `seqwish_replacement_config`
and the `SweepgaAlignConfig` (filter knobs #2, #3, #6, #7) had no
observable effect — those knobs were already unused-at-runtime via the
`no_filter=true` short-circuit.

## Finding

**SweepGA's filter mode is NOT the bottleneck.**

The hypothesis stated in the task description — that there is still a
1:1-vs-N:N mapping/scaffold filter mode toggle dropping SweepGA's
all-vs-all alignments below seqwish — was wrong. Verifying:

1. With `no-filter=true`, the SweepGA library's PAF filter is skipped at
   `src/library_api.rs:321` (sweepga revision `ddd31d3`).
2. With `no-filter=true`, impg's own filter (`filter_generated_paf`)
   early-returns at `src/commands/graph.rs:677` and never reads
   `num_mappings` / `scaffold_filter` / `scaffold_mass`.
3. The SweepGA library's defaults are *already* `num_mappings=many:many`
   and `scaffold_filter=many:many`. So at every layer that the filter
   might have inspected, `many:many` was already the effective mode.
4. Pinning the modes to `many:many` defensively (this PR) produced a
   byte-difference at the GFA layer of: **−82 segments, −82 links,
   +162 bp** — a measurable but small effect attributable entirely to
   the polish method change (sPOA → POASTA) on the same 8 plans.

**The 6,670-segment gap vs PGGB lives elsewhere.** Candidate locations
ranked by likelihood:

1. **seqwish induction inside the bimodal plan-1 bubble.** Plan 1
   (root span 110 bp, 312 traversals, max-len 31 478 bp, total
   5.13 Mbp) is the dominant locus for short-vs-long mismatch — per
   `docs/crush-aligner-deep-diag.md` Q1, the SweepGA PAF for this
   plan has 77 666 lines but the longest internal CIGAR match-run is
   40–99 bp, and 311 bp seqwish-k cannot anchor short ↔ long pairs.
   Even with `seqwish-k=1` forced (this PR's `no-filter=true`
   behavior), the dependency's transitive closure may still
   fragment short repeats into N-way disjoint per-traversal copies.
2. **POASTA polish coverage threshold.** Inside-bubble re-decomposition
   only descends into sub-bubbles whose median traversal length is
   below `polish-max-median-traversal-len=1k` and max ≤ 10k. Larger
   stringy bubbles (the kind PGGB also struggles with) are skipped
   here. Lifting the polish caps could compress them further.
3. **Within-replacement byte-duplicate fragmentation.** Even after
   POASTA polish on sub-bubbles, the *outer* SweepGA-induced
   replacement may carry per-traversal disjoint segments from
   seqwish's transitive closure. A POVU-rooted second-pass that
   walks the outer SweepGA replacement and applies POASTA to its
   *external* sibling-merges (not just nested children) is unlikely
   to help; the segments are byte-distinct, so POASTA cannot collapse
   them.
4. **Aligner choice for the short-leg.** FastGA produces inflate-then-
   chain alignments tuned for whole-genome scale; wfmash with a low
   `-p` (e.g. 70 %) would emit short anchored mappings that seqwish
   could close. Compare to PGGB's wfmash-based induction, which is
   precisely what builds the 13 288-segment graph.

### Proposed next experiments

- `crush-wfmash-replacement-induction` — same command but
  `aligner=wfmash,map-pct-identity=70`. If S drops near PGGB's
  13 288, that confirms the bottleneck is the aligner's anchor
  granularity, not any filter.
- `crush-lift-polish-caps` — raise
  `polish-max-traversal-len` to 50k and
  `polish-max-median-traversal-len` to 10k so the polish loop also
  descends into the bimodal plan-1 bubble's interior.
- `crush-poasta-on-outer-replacement` — apply POASTA directly to the
  SweepGA-produced GFA before declaring the replacement (instead of
  only on sub-bubbles).

## Hard-gate checklist

- [x] Filter inventory in docs (every knob with file:line + current value) — see table above
- [x] Code change committed with the fully-unfiltered + many-to-many config exposed
- [x] `cargo build --release` succeeds
- [x] `cargo test --lib` (281) and `cargo test --bin impg` (56) pass; one pre-existing failure in `test_crush_integration::nested_bubble_level_descent_actually_descends` reproduces on the pre-change tree and is annotated as known in its `panic` message at `HEAD 471f089` (not introduced by this PR)
- [x] Per-plan PAF retention reported: 100 % at every stage by construction
- [x] Real C4 GRCh38 run completes, 465/465 paths preserved
- [x] Final metrics: 19 958 S, 23 822 L, 544 641 bp, 499 trivial-stringy
- [x] PNG uploaded as `c4-crush-many-to-many-poasta-polish.png` — `https://hypervolu.me/~erik/impg/c4-crush-many-to-many-poasta-polish.png`
- [x] `docs/crush-sweepga-many-to-many-poasta-polish.md` committed
- [x] `wg artifact crush-sweepga-many docs/crush-sweepga-many-to-many-poasta-polish.md`
- [x] Result matches prior baselines (within 0.4 % on segments) — that's the finding: **sweepga's filter mode is NOT the bottleneck**; the gap lives in the aligner-anchor / seqwish-induction layer, with three concrete next experiments proposed above.
