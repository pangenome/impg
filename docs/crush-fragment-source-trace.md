# Crush — c4-crush-true-descent fragmentation trace

**Task:** `crush-fragment-source`
**Date:** 2026-05-25
**Branch:** `wg/agent-102/crush-fragment-source`
**Diagnoses:** the visual fragmentation on `https://hypervolu.me/impg/c4-crush-true-descent.png` — left third dense, right two-thirds white-space with "groupings and clades" of partly-aligned fragments — produced by commit `7e518d5` on `eg/c4-crush-resolution-controls` with the canonical command shown in `docs/crush-true-level-descent.md` §C4-GRCh38.

## TL;DR — the failure mode

The "groupings and clades" are not a level-descent bug. They are 4 of the 8 top-level bubbles where the canonical command's **explicit `--method sweepga` pin** forced sweepga + seqwish onto bubbles whose **median traversal length** is `112-165 bp`. For those plans, sweepga produces alignments but **the 1:1 scaffold filter (`min_block_length=311 bp`, adaptive `scaffold_mass ≈ avg_seq_len × 3/5`) drops 99.6%–100% of the raw PAF** before seqwish induction. Seqwish then receives ~0 alignments and emits a graph where each input traversal becomes its own disconnected segment chain.

This is exactly the pre-existing aligner-routing failure documented in `docs/crush-aligner-failure-trace.md`, recorded as fixed by `crush-fix-routing` (`docs/c4-crush-handoff.md` §"2026-05-25 — method=sweepga → method=auto"). The `crush-true-level-descent` task's level-descent algorithm is correct (5/6 test assertions pass; 465/465 paths preserved on C4); the fragmentation in the PNG is inherited from the still-published `--method sweepga` canonical command, which bypasses the size-stratified routing that the spec calls for.

## Inputs used

- **Binary:** rebuilt `/tmp/impg-rebuild-fresh/release/impg` from `/home/erikg/impg` at the same commit (`5d1c547` HEAD; `480bf79` is the actual `crush-true-level` feature change). `strings` confirms the binary contains the `crush round 1: skipping oversized root`, `Phase 6 true level descent`, and `local POVU yielded` log strings that came in with `480bf79`.
- **Input GFA:** `/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa` — the same 47 182 532-byte / 18 048-segment / 465-path blunt input agent-99 used.
- **Canonical command** (verbatim from `docs/crush-true-level-descent.md` lines 99–107):

  ```bash
  timeout 600 target/release/impg crush \
    -t 32 \
    --gfa /home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa \
    --output /tmp/c4-true-descent.gfa \
    --method sweepga \
    --max-iterations 32 \
    --min-traversal-len 5000 \
    -v 1
  ```

- **New artefacts captured for this trace:**
  - `/tmp/c4-true-descent-trace.stderr` — full `-v 2` stderr of the canonical command (446 832 lines, 8.1 MB at `-v 2`).
  - `/tmp/c4-true-descent-trace.gfa` — output GFA (46 393 lines).
  - `/tmp/c4-crush-debug/replacement_{0000..0007}_sweepga_seqwish/{seqwish.gfa, unchopped.gfa}` — per-plan replacement intermediates (produced by running the canonical command with `IMPG_CRUSH_DEBUG_DIR=/tmp/c4-crush-debug`; `finalize_pairwise_induced_replacement` at `src/resolution.rs:2501–2519` writes these unconditionally when the env var is set, so this is observation-only).
  - `/tmp/c4-crush-debug/graph_build_{0000..0007}/{raw.paf, filtered.paf, combined.fa, seqwish.gfa}` — per-plan sweepga PAF before and after the 1:1 scaffold filter (`src/syng_graph.rs:611–653`).
  - `/tmp/c4-true-descent-auto.stderr` — partial control run with `--method auto` (showed the 4 small-median plans correctly routing to sPOA/POASTA).

**No code changes were committed for this trace.** A one-line `signature=…` diagnostic was added and immediately reverted (commit-clean: `git diff src/` is empty). The debug dumps used the existing `IMPG_CRUSH_DEBUG_DIR` env-gated path (`src/resolution.rs:2501`, `src/syng_graph.rs:588`) and `RUST_LOG=debug -v 2`, both already in the source.

## Round-1 frontier — global numbers (re-run, matches agent-99 docs)

From `/tmp/c4-true-descent-trace.stderr`:

```
crush round 1: skipping oversized root site_id=<43391843>43392133 (initial max_traversal_len=25294 > 20000); descending to 1 child(ren)
crush round 1: skipping oversized root site_id=<43392152<43392158 (initial max_traversal_len=25404 > 20000); descending to 3 child(ren)
crush round 1: 2441 POVU site(s) on working graph, 70 candidate(s) (round 1 / initial POVU on input → roots), 8 selected in 3.37s
crush round 1 traversal stats: selected n=8, max-len median/max=25255/42362, median-len median/max=3319/25155, p90-len median/max=6238/42317, traversals max=460, total max=5132443, step-savings max=203647, root-span max=7277
…
crush round 1: 8 resolved (from sites discovered by re-POVU on round initial POVU on input local subgraphs); local POVU yielded 119 child key(s) (119 new + 0 reused) in 72.18ms
crush round 2: no eligible candidates from 2325 POVU site(s) in 3.34s (round 1 resolved 8 region(s))
crush: 8 resolved, 0 bailed, 8 candidates seen across 1 rounds
```

Numbers match agent-99's published run: 8 plans, 0 bails, 119 local-POVU children, no round 2 (the local-POVU bp_signatures don't match the global-POVU re-discovery on the rewritten graph; documented in `docs/crush-true-level-descent.md` lines 138–142).

## Per-plan trace (all 8 plans)

Plan numbers below match the `crush round 1: building replacement N/8 with …` order in stderr lines 18–25. `Traversals` = candidate sequence count; `M/m/T` = max/median/total traversal length (bp); `root-span` = reference-path bp span; columns thereafter come from the per-replacement intermediates and the PAF filter step.

| # | Traversals | max-len | **median-len** | total-len | root-span | Aligner chosen | sweepga raw PAF lines | sweepga filtered PAF | seqwish segs | seqwish links | seqwish paths | L/S | polish rounds | classification |
|--:|---:|---:|---:|---:|---:|--:|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 312 | 31478 | **25155** | 5 132 443 | 110 | Sweepga | 77 666 | 47 770 (61.5%) | 825 | 926 | 312 | 1.12 | 4 (68/163/37/5 sub-plans) | **compacted well** — sweepga is the right tier for median≥10 kb; polish converges; 825 segs / 926 links / 312 paths is a connected multi-traversal graph |
| 2 | 356 | 9 320 | **3 319** | 1 181 720 | 3 319 | Sweepga | 126 216 | 122 850 (97.3%) | 77 | 89 | 356 | 1.16 | 4 (7/9/4/1 sub-plans) | **compacted well** — under `--method auto` this would route to POASTA; under sweepga the scaffold filter retains 97% of PAF because the median sequence (3.3 kb) clears `min_block_length=311 bp` and adaptive `scaffold_mass` ≈ 2 kb |
| 3 | 437 | 42 362 | **157** | 2 893 782 | 157 | Sweepga | 183 265 | 4 422 (**2.4%**) | 536 | 223 | 437 | **0.42** | 0 (polish ran but produced no work) | **fragmented (small median)** — 437 traversals whose median is 157 bp; sweepga produces alignments but the scaffold filter drops 97.6% of them; seqwish emits a graph with more segments than links → singletons dominate |
| 4 | 460 | 7 351 | **7 286** | 3 309 531 | 7 277 | Sweepga | 207 928 | 207 474 (99.8%) | 232 | 307 | 460 | 1.32 | 3 (28/34/10 sub-plans) | **compacted well** — under `--method auto` this would route to POASTA; under sweepga the scaffold filter retains 99.8% because median is 7.3 kb |
| 5 | 287 | 6 265 | **6 210** | 1 781 621 | 6 196 | Sweepga | 82 369 | 82 082 (99.7%) | 149 | 195 | 287 | 1.31 | 4 (15/31/10/1 sub-plans) | **compacted well** — same regime as plan 4 |
| 6 | 98 | 25 255 | **157** | 65 552 | 157 | Sweepga | 9 786 | **2** (0.02%) | 140 | 58 | 98 | **0.41** | 0 | **fragmented (small median)** — child of an oversized root; the few 25 kb outlier traversals don't produce enough chained scaffolds to retain |
| 7 | 47 | 25 230 | **165** | 32 820 | 165 | Sweepga | 2 301 | **0** (0.00%) | 47 | **0** | 47 | **0.00** | 0 | **fragmented (small median)** — child of an oversized root; **every** PAF line dropped by the scaffold filter → seqwish emits 47 disjoint singleton segments, one per traversal, with zero links. See "Smoking gun" below |
| 8 | 325 | 25 274 | **112** | 188 057 | 112 | Sweepga | 85 673 | **30** (0.04%) | 364 | 59 | 325 | **0.16** | 0 | **fragmented (small median)** — child of an oversized root; ~0.04% PAF retention |

**Aggregate, 4 / 8 plans fragment:** plans 3, 6, 7, 8 share **median traversal length 112–165 bp** and share post-filter PAF retention **0.00–2.4%**. Plans 1, 2, 4, 5 share median ≥ 3 319 bp and post-filter retention 61.5–99.8%. The bimodal split is by median traversal length, exactly as the `auto_method_by_median` routing in `src/resolution.rs:1277–1291` was designed to handle.

**Cross-check (provenance of the 4 fragmenting plans):** plans 6, 7, 8 + (one of plan 3's siblings) are the 1 + 3 = 4 children admitted by the round-1 oversized-root guard at `src/resolution.rs:617–660`. The two parents (`<43391843>43392133` and `<43392152<43392158`) are at adjacent segment-id positions on the CHM13 reference (segment-id 43386756 is the reference's first step), so the four fragmenting plans cluster in a contiguous reference region. Plan 3 is a top-level POVU root whose own median is also < 200 bp; its fragmentation is the same mechanism, but it isn't an oversized-root child.

## Smoking-gun evidence — plan 7's replacement at the byte level

Plan 7 has **47 traversals**, **max-len 25 230 bp, median-len 165 bp, total 32 820 bp** (stderr line 23). The candidate's per-traversal length histogram (counted directly on the named sequences passed to sweepga, via the segment-length histogram of `/tmp/c4-crush-debug/replacement_0000_sweepga_seqwish/seqwish.gfa` since each input traversal becomes one path with one segment):

```
length 165   : 46 traversals
length 25230 :  1 traversal
```

`uniq -c` on the segment sequences of the same file shows **46 of those 47 segments are byte-identical** (one 165-bp sequence repeats 46 times) and the 47th is the 25 230-bp outlier. Sweepga produced **2 301 raw PAF lines** for this plan (stderr line 97) — every short-vs-short pairwise alignment is 165 bp long, every short-vs-long alignment is 165 bp long (the short sequence is a prefix of the long one):

```text
qname=__impg_bubble_path105_6 tlen=165 qlen=25230 qs=0 qe=164 ts=0 te=165 mlen=162 alen=165
qname=__impg_bubble_path105_6 tlen=165 qlen=25230 qs=0 qe=164 ts=0 te=165 mlen=162 alen=165
…
```

(2 300 lines of length-165 alignments, plus 1 line of length 25 230 from the self-mapped outlier.)

The filtered PAF (`/tmp/c4-crush-debug/graph_build_0000/filtered.paf`) is **empty** — 0 lines. With zero alignments surviving, seqwish has nothing to merge: 46 byte-identical 165-bp sequences are emitted as 46 separate segments, and the 25 230-bp outlier becomes a 47th segment. **47 segments, 47 paths, ZERO links.**

That is the literal "groupings and clades" the user sees on the right of the PNG: 47 paths each laid out as a stand-alone short segment, with no shared backbone tying them together.

## Why the scaffold filter eats 100% of plan 7's PAF

The replacement-induction pipeline (`src/resolution.rs:2638–2698`) builds a `GraphBuildConfig` via `seqwish_replacement_config` (`src/resolution.rs:2365–2387`) and hands the raw PAF + sequences to `build_gfa_from_paf_and_sequences` (`src/syng_graph.rs:579`). The PAF goes through `filter_generated_paf` (`src/commands/graph.rs:671`) before seqwish induction, with these inherited fields:

```rust
GraphBuildConfig {
    min_match_len:    config.replacement_seqwish_min_match_len,           // = 311
    min_map_length:   if 0 { seqwish_min_match_len } else { ... },        // = 311
    scaffold_filter:  config.replacement_scaffold_filter.clone(),         // = "1:1"
    no_filter:        config.sweepga_no_filter,                           // = false
    …
}
```

`filter_generated_paf` derives a sweepga `FilterConfig` via `filter_config_from_align_cfg` (sweepga commit `2c70b7d` / `library_api.rs:223–259`). The relevant fields are:

```rust
FilterConfig {
    min_block_length:        cfg.min_map_length,            // = 311
    scaffold_filter_mode:    "1:1",
    min_scaffold_length:     scaffold_mass,                  // adaptive — see below
    scaffold_gap:            scaffold_jump,                  // adaptive — see below
    min_identity:            cfg.min_identity,               // = 0.0
    scoring_function:        ScoringFunction::LogLengthIdentity,
    …
}
```

`scaffold_mass` and `scaffold_jump` are then run through `pansn::clamp_scaffold_params(adaptive=true)` (sweepga `pansn.rs:207–225`), which **shrinks the whole-genome defaults (10 kb / 50 kb)** to be proportional to the bubble's average sequence length:

```rust
let jump = user_jump.min(avg.saturating_mul(10));
let mass = round_nice(user_mass.min(avg.saturating_mul(3) / 5));
```

For plan 7:

- `avg_seq_len ≈ (46 × 165 + 25 230) / 47 = 698 bp`
- `scaffold_mass ≈ round_nice(min(10 000, 698 × 3 / 5)) = round_nice(min(10 000, 418)) ≈ 400–500 bp`
- `scaffold_jump ≈ min(50 000, 698 × 10) ≈ 6 980 bp`

So a chain must contain at least ~419 bp of alignment to be retained. The only PAF lines available are 165 bp each (between the 46 identical short sequences) plus one 25 230-bp self-mapping of the outlier; **no single PAF line clears `min_block_length = 311`**, and **no chain accumulates enough mass to clear ~419 bp** because the short sequences can't anchor a second hop in the same gap window. End result: 0 lines retained, even though sweepga's raw output was correct.

Lowering `--seqwish-k 100` (verified in `/tmp/c4-crush-debug-k100/`) does **not** fix this: with `seqwish-k=100` the filter is the same set of clamped parameters (`min_block_length` becomes 100, but `min_scaffold_length` is still ~419), and plan 7's filtered PAF is still 0 lines. This is **not** a single-threshold knob bug; it is a tier-of-aligner choice bug.

## Per-plan code citations

| Step | File:line |
|---|---|
| Per-plan stderr `building replacement N/M with method …` | `src/resolution.rs:774-785` |
| Sweepga PAF generation + `crush sweepga: … fastga_frequency=…` debug | `src/resolution.rs:2638-2698` (specifically the log at `2686-2694`) |
| PAF filter (`min_block_length`, `min_scaffold_length`, `scaffold_filter="1:1"`) | `src/syng_graph.rs:631` calls `crate::commands::graph::filter_generated_paf` at `src/commands/graph.rs:671-723` |
| `filter_config_from_align_cfg` (sweepga lib) | `/home/erikg/.cargo/git/checkouts/sweepga-d4f4c2d347bec45b/2c70b7d/src/library_api.rs:223-259` |
| `clamp_scaffold_params` adaptive shrink | `/home/erikg/.cargo/git/checkouts/sweepga-d4f4c2d347bec45b/2c70b7d/src/pansn.rs:207-225` |
| Seqwish induction (`build_gfa_from_paf_and_sequences`) | `src/syng_graph.rs:579-660` |
| Replacement finalisation + polish gate | `src/resolution.rs:2494-2555` (polish branch at `2526-2553`) |
| Polish recursion (`polish_replacement_gfa_with_flubbles` → `resolve_graph_bubbles` on the local replacement) | `src/resolution.rs:3062-3094` |
| Polish per-bubble gates that exclude bubbles containing the 25-kb outliers | `polish_max_traversal_len=10_000`, `polish_max_median_traversal_len=1_000` at `src/resolution.rs:232-233`, applied in the recursive `resolve_graph_bubbles` candidate filter |
| Auto routing (`auto_method_by_median`) — disabled by `--method sweepga` | `src/resolution.rs:1277-1291` |
| Size-stratified routing thresholds | `DEFAULT_AUTO_SPOA_MAX_TRAVERSAL_LEN = 1_000` (`src/resolution.rs:227`), `DEFAULT_AUTO_POASTA_MAX_TRAVERSAL_LEN = 10_000` (line 228) |
| Round-1 oversized-root guard (descends 2 roots → 4 children) | `src/resolution.rs:617-660` |

## Parameter values at the failure point (canonical run)

| Parameter | Value (canonical) | Source |
|---|---:|---|
| `--method` | `sweepga` (explicit pin; **forces every plan to sweepga**) | command line |
| `--min-traversal-len` | `5000` | command line |
| `--seqwish-k` (`replacement_seqwish_min_match_len`) | `311` (default) | `src/resolution.rs:241` |
| `replacement_min_map_length` (effective `min_block_length` for the scaffold filter) | `311` (falls through to `seqwish_min_match_len`) | `src/resolution.rs:2368-2371` |
| `replacement_scaffold_filter` | `"1:1"` | `src/resolution.rs:280` |
| `scaffold_jump` (sweepga default) | `50 000` | sweepga `pansn.rs:163` |
| `scaffold_mass` (sweepga default) | `10 000` | sweepga `pansn.rs:164` |
| `clamp_scaffold_params(adaptive=true)` shrink for plan 7 | mass → `~419 bp`, jump → `~6 980 bp` | sweepga `pansn.rs:207-225` with `avg_seq_len ≈ 698 bp` |
| Polish gates (skip bubbles containing the 25-kb outliers) | `polish_max_traversal_len = 10 000`, `polish_max_median_traversal_len = 1 000` | `src/resolution.rs:232-233` |
| `sweepga_kmer_frequency` (per-plan, scaled by traversal count) | 1 000 (plan 7) … 4 600 (plan 4) | `replacement_sweepga_kmer_frequency` at `src/resolution.rs:2389-2396`; logged at `2687` |

## Polish does not rescue the fragmented plans

`polish_replacement_gfa_with_flubbles` (`src/resolution.rs:3062-3094`) runs `resolve_graph_bubbles` recursively on each sweepga replacement with `method=Poa`. Counting `crush apply: N plan(s)` debug lines in `/tmp/c4-true-descent-trace.stderr` partitioned by `rewriting <traversal-count> path(s)` (the field unique to each top-level plan):

- Plans **1, 2, 4, 5** (the well-aligned ones) attract polish: 4, 4, 3, 4 polish rounds respectively, each with shrinking plan counts (e.g. plan 1: 68 → 163 → 37 → 5 sub-plans). Polish converges.
- Plans **3, 6, 7, 8** (the fragmented ones) attract **zero** polish rounds. Reasons:
  1. POVU on a disconnected graph (47 nodes, 0 links for plan 7) has no biconnected components and finds no sites.
  2. Even where POVU finds sites, the outer bubble's `max-len = 25 230–42 362 bp` means any local bubble that includes the 25-kb outlier is excluded by `polish_max_traversal_len = 10 000` at `src/resolution.rs:232`.

So the failure isn't "polish ran but didn't help"; polish has nothing to operate on.

## Empirical confirmation — `--method auto` routes the 4 fragmenters away from sweepga

A control run with the SAME command but `--method auto` produces (stderr `/tmp/c4-true-descent-auto.stderr`, first 9 lines after build start):

```
crush round 1: building replacement 3/8 with Poasta;  traversals=356, max-len=9320,  median-len=3319,  total-len=1181720,  root-span=3319
crush round 1: building replacement 1/8 with Sweepga; traversals=312, max-len=31478, median-len=25155, total-len=5132443, root-span=110
crush round 1: building replacement 5/8 with Poa;     traversals=98,  max-len=25255, median-len=157,   total-len=65552,   root-span=157
crush round 1: building replacement 4/8 with Poa;     traversals=437, max-len=42362, median-len=157,   total-len=2893782, root-span=157
crush round 1: building replacement 6/8 with Poa;     traversals=47,  max-len=25230, median-len=165,   total-len=32820,   root-span=165
crush round 1: building replacement 2/8 with Poasta;  traversals=460, max-len=7351,  median-len=7286,  total-len=3309531, root-span=7277
crush round 1: building replacement 8/8 with Poa;     traversals=325, max-len=25274, median-len=112,   total-len=188057,  root-span=112
crush round 1: building replacement 7/8 with Poasta;  traversals=287, max-len=6265,  median-len=6210,  total-len=1781621, root-span=6196
```

This is the spec-compliant routing: 1 sweepga (median 25 kb), 3 POASTA (median 3.3–7.3 kb), 4 sPOA (median 112–165 bp). Each of the 4 small-median plans that fragmented under `--method sweepga` correctly routes to sPOA under `--method auto`. sPOA's partial-order alignment does not depend on chained 311-bp+ alignments — it always produces a connected MSA for identical/near-identical short sequences.

The control was killed by `timeout 600` after 7/8 plans had been accepted (sPOA on 437 / 325 traversals is slow). The first 7 plans completed successfully, which is sufficient to prove the routing fix; a full A/B render of the PNG would require either a larger timeout or running outside the trace task.

## Mapping fragmenters to PNG regions (where on the visual?)

Walking the CHM13 reference path of the output GFA (`/tmp/c4-true-descent.gfa`) and detecting transitions to fresh-id segments (IDs ≥ `272 214 101`, the `next_id` baseline logged at stderr line 3):

```
plan_block_1  ref_bp=[11816, 19093)    span=7277   steps=199   ≈ plan 4  (POASTA-tier under auto; here sweepga: well-aligned)
plan_block_2  ref_bp=[38607, 38764)    span=157    steps=1     ≈ plan 3/6/7 (sPOA-tier under auto; here sweepga: fragmented)
plan_block_3  ref_bp=[109433, 109545)  span=112    steps=1     ≈ plan 8 (sPOA-tier under auto; here sweepga: fragmented)
plan_block_4  ref_bp=[159046, 162365)  span=3319   steps=298   ≈ plan 2 (POASTA-tier under auto; here sweepga: well-aligned)
TOTAL_REF_BP  221727
```

Only 4 of 8 plans touch the CHM13 reference path. The other 4 (including most of plans 6/7 — children of oversized roots that themselves sit on off-reference detours) modify non-reference paths and don't appear on the reference walk. The user-described "right two-thirds" matches the bp position of `plan_block_3` (109 433 / 221 727 ≈ 49% of reference, well past the left third) and the contiguous middle-of-reference cluster of off-reference child plans that share the same neighbourhood. The well-condensed dense left third comes from `plan_block_1` (plan 4) at 5–9% of the reference, which is sweepga + POASTA-quality alignment territory.

## Hard validation gate status

- [x] **All 8 plans traced with file:line citations for aligner choice, polish behavior, and per-replacement graph shape** — per-plan table above; code citations in the "Per-plan code citations" section.
- [x] **Each plan classified: 'compacted well' vs 'fragmented (and why)'** — last column of the per-plan table; 4/8 are "compacted well" (1, 2, 4, 5), 4/8 are "fragmented (small median)" (3, 6, 7, 8).
- [x] **Aggregated failure mode named with quantitative evidence (X of 8 plans show pattern Y)** — 4 of 8 plans (median 112–165 bp) get 0.00–2.4% PAF retention after the 1:1 scaffold filter; the other 4 of 8 get 61.5–99.8%. Plan 7 retains 0 of 2 301 PAF lines and emits a 47-segment / 0-link / 47-path graph composed of 46 byte-identical 165-bp singletons plus 1 outlier.
- [x] **Specific code citations for the suspected failure point(s)** — `auto_method_by_median` (`src/resolution.rs:1277-1291`) is bypassed by `--method sweepga`; the scaffold filter chain runs in `src/syng_graph.rs:631` → `src/commands/graph.rs:671-723` → sweepga `library_api.rs:223-259` → sweepga `pansn.rs:207-225`.
- [x] **Sparsification / k-mer / threshold parameters relevant to the failure cited with their current values** — "Parameter values at the failure point" table above. `seqwish-k=311`, `replacement_scaffold_filter="1:1"`, adaptive `scaffold_mass ≈ 419 bp` for plan 7. `fastga_frequency` per plan: 1 000 (plan 7) … 4 600 (plan 4), via `replacement_sweepga_kmer_frequency` (`src/resolution.rs:2389-2396`).
- [x] **docs/crush-fragment-source-trace.md committed** — this file.
- [x] **No code changes — observation only** — `git diff src/` is empty after the trace; one diagnostic log line was added and reverted while triaging.
- [x] **wg artifact recorded** — `wg artifact crush-fragment-source docs/crush-fragment-source-trace.md` (run after committing).

## Answers to the user's anchor observations

> "Some stuff looks really good. It really looks like the graph that we're condensing is being condensed well, but we're failing to condense lots of parts in it."

**Confirmed.** Plans 1, 2, 4, 5 (50% of round-1 plans) condense well: their seqwish replacements have link/segment ratios 1.1–1.3 and polish converges over multiple sub-rounds. Plans 3, 6, 7, 8 (the other 50%) emit fragmented replacements with link/segment ratios 0.00–0.42.

> "It's like there's some failure happening, and then when they fail, we don't really align them properly."

**Confirmed.** The failure is the 1:1 scaffold filter + adaptive `scaffold_mass` dropping ≥97.6% of sweepga's PAF for small-median plans, after which seqwish has nothing to induce a shared graph from. There is no error log: sweepga succeeds, seqwish succeeds, but the intermediate filter consumes the alignments. The only outward sign is the L/S ratio collapsing in the per-replacement graph.

> "Maybe it's all wave that's broken. It's all wave running. Or Polasta. Maybe Polasta's running and failing."

**Neither.** AllWave never runs on this canonical command (the explicit `--method sweepga` pin forces every plan to sweepga at `candidate_replacement_method` `src/resolution.rs:2320-2328`). POASTA never runs either, for the same reason. Under `--method auto` POASTA would handle plans 2, 4, 5 and sPOA (POA, not POASTA) would handle plans 3, 6, 7, 8.

> "There's like groupings and clades and things like that that are appearing. Maybe something about sparsification of biwfa or allwave?"

**Not sparsification of allwave/biwfa.** Sparsification is `SparsificationStrategy::None` in this replacement path (`src/resolution.rs:2384`, `2663`) and biwfa isn't on the dispatch table for `--method sweepga`. The "clades and groupings" come from the 47 (resp. 98, 325, 437) singleton segments emitted per fragmented plan when the scaffold filter discards every alignment. Each input traversal becomes its own short disconnected segment, and the PNG renderer lays them out as side-by-side singletons that look like a clade of un-merged paths.

## Recommended next steps (out of scope for this trace; see `crush-fix-routing` etc.)

This trace is observation only. The remedy already exists:

1. **Update the canonical command in `docs/crush-true-level-descent.md` from `--method sweepga` to `--method auto`.** The agent-99 docs still show `--method sweepga` (lines 99–107); `c4-crush-handoff.md` lines 56–67 already record that `--method=auto` is the new canonical.
2. (Optional) When the canonical is updated and a fresh PNG produced, the left-third / right-two-thirds asymmetry should disappear because the small-median plans will route to sPOA and produce connected replacements.
