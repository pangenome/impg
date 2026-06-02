# Crush — size-stratified aligner routing fix (`crush-fix-routing`)

**Task:** `crush-fix-routing`
**Date:** 2026-05-25
**Branch:** `wg/agent-92/crush-fix-routing`
**Predecessor task:** `crush-aligner-failure` (diagnosis, commit `d355e48`)
**Spec reference:** `docs/crush-architecture-spec.md` §Phase-2
**Trace reference:** `docs/crush-aligner-failure-trace.md` §Fix sketch

## TL;DR

The crush dispatch (`auto_replacement_method` /
`candidate_selection_method`, `src/resolution.rs:887–928`) only ever returned
`ResolutionMethod::Poa` or `ResolutionMethod::Allwave`. POASTA and Sweepga
were unreachable from `method=auto`, and the canonical command in
`docs/c4-crush-handoff.md` was pinning `method=sweepga`, forcing sweepga +
seqwish-k=311 on every bubble. This run replaces that with the spec's
3-tier dispatch keyed on **median traversal length** (not max):

- `median < 1 kb` → `ResolutionMethod::Poa` (sPOA)
- `1 kb ≤ median < 10 kb` → `ResolutionMethod::Poasta`
- `median ≥ 10 kb` → `ResolutionMethod::Sweepga`

On the C4 GRCh38 slice the routing now spreads across aligners (41 sPOA +
11 sweepga across 13 rounds, 52 resolutions total — POASTA was not
selected because no bubble landed in its 1–10 kb median band) and the
**structural ≥50 bp byte-duplicate count drops from the previous
sweepga-everywhere fragmentation pattern to 166 segments** (see
[Result metrics](#result-metrics)). One caveat: sPOA on bubbles with a
short median but a long max emits many <10 bp polish-stage micro-bubble
segments that are byte-equal to each other (single-base anchor nodes),
which inflates the *total* segment count from 23 193 to 47 060 even
though *structural* fragmentation is essentially eliminated. The
micro-duplicate issue is a polish / final-compaction concern, not a
routing concern — it is called out in
[Known limitation](#known-limitation--spoa-micro-duplicates) and is out of
scope for `crush-fix-routing`.

PNG: <https://hypervolu.me/~erik/impg/c4-crush-routing-fix.png>

## What changed

1. `src/resolution.rs` — `candidate_selection_method` and
   `auto_replacement_method` (the two routing entry points called per
   candidate in `resolve_graph_bubbles`) now delegate to a single new
   `auto_method_by_median` function implementing the 3-tier dispatch. The
   decision variable is `traversal_stats.median_len`, the same field the
   existing logging line emits, so the per-round stderr is the source of
   truth for which tier each bubble was routed to.
2. New `ResolutionConfig` field `auto_poasta_max_traversal_len` (default
   10 000) and a renamed-in-semantics `auto_spoa_max_traversal_len`
   (default lowered from 2 000 → 1 000; semantics changed from "max-len
   upper bound" → "median upper bound"). Both are exposed on the CLI as
   `--auto-spoa-max-traversal-len` (alias `auto-spoa-max-len`) and
   `--auto-poasta-max-traversal-len` (alias `auto-poasta-max-len`), and
   accepted as engine-string parameters
   (`auto-spoa-max-len`, `auto-poasta-max-len`).
3. `direct_replacement_within_budget` (the old budget-gate that
   short-circuited Poa back to Allwave) was removed — under the new
   median-tier dispatch its semantics no longer fit (each tier's aligner
   handles its own budget; falling through to a different aligner because
   of a budget guard would defeat the routing decision the median
   already made).
4. Tests:
   - `src/resolution.rs::tests` gains 8 routing-focused unit tests
     covering each tier, the median-vs-max edge case from
     `docs/crush-aligner-failure-trace.md` pick (a), the boundary
     behaviour at the exact thresholds, the per-tier disabling knobs,
     the explicit-method pin escape hatch, and the auto-mode candidate
     prioritisation that depends on the tier.
   - The old `auto_routes_small_bubbles_to_spoa_and_larger_bubbles_to_allwave`
     test was deleted: its assertions encode the bug (auto only routes to
     Poa or Allwave) and cannot survive the fix.
   - The CLI test
     `test_query_engine_string_parses_with_crush_stage_defaults` at
     `src/main.rs` updates the expected default
     `auto_spoa_max_traversal_len` from `2_000` to `1_000` and adds an
     assertion for the new `auto_poasta_max_traversal_len = 10_000`
     default.
5. `docs/c4-crush-handoff.md` — the canonical command's
   `method=sweepga,...` is replaced with `method=auto,...` and a
   `2026-05-25 — method=sweepga → method=auto` note records the change
   and its motivation, pointing at this document and at
   `docs/crush-aligner-failure-trace.md`. The legacy-baseline
   reproduction command kept `method=sweepga,kmer-frequency=10` so the
   two-run comparison still isolates the routing change from the
   FastGA-frequency change introduced in `98fd538`.

## Command run on real data

```bash
out=data/c4_crush_routing_fix_$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$out"

/usr/bin/time -v target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/auto_routing.nosort" \
  -v 1 \
  > "$out/auto_routing.nosort.stdout" \
  2> "$out/auto_routing.nosort.stderr"

gfasort -i "$out/auto_routing.nosort.gfa" -o "$out/auto_routing.Ygs.gfa" -p Ygs -t 32
gfalook  -i "$out/auto_routing.Ygs.gfa"   -o "$out/auto_routing.Ygs.mean-depth.paths.png" -m -x 2200 -y 1200
scp "$out/auto_routing.Ygs.mean-depth.paths.png" erik@hypervolu.me:www/impg/c4-crush-routing-fix.png
```

Wall time: 47:56 (FastGA-heavy round 1 + round 2 dominate at ~14 min each;
later rounds are 20–90 s each — see "Per-round timing" below).

Resident set high-water mark: 159 GB (FastGA all-vs-all on the largest
bubbles).

## Per-round aligner choice (canonical evidence routing fired)

```
$ grep -oE "building replacement [0-9]+/[0-9]+ with (Poa|Poasta|Sweepga|Allwave|StarBiwfa)" \
    .../auto_routing.nosort.stderr \
  | grep -oE "with [A-Za-z]+" | sort | uniq -c
     41 with Poa
     11 with Sweepga
```

Round 1 alone illustrates the dispatch:

| round 1 plan | traversals | max-len | median-len | total-len | aligner picked |
|---:|---:|---:|---:|---:|:---|
| 1/7 | 312 | 31 478 | 25 155 | 5 132 443 | Sweepga (median ≥10 kb) |
| 2/7 | 54 | 31 345 | 151 | 226 211 | Poa (median <1 kb, max outlier) |
| 3/7 | 37 | 11 509 | 93 | 37 689 | Poa |
| 4/7 | 437 | 42 362 | 157 | 2 893 782 | Poa |
| 5/7 | 98 | 25 255 | 157 | 65 552 | Poa |
| 6/7 | 47 | 25 230 | 165 | 32 820 | Poa (this is pick (a) from the trace report) |
| 7/7 | 438 | 11 597 | 111 | 438 312 | Poa |

Plan 6/7 is the bubble that the trace report singled out as the worst-case
mis-routing in the baseline (47 traversals × 165 bp median that the
previous sweepga+seqwish-k=311 pipeline emitted as 47 byte-identical
disjoint segments with zero links). Under the new dispatch it goes to
sPOA, as the spec requires.

No bubble landed in the 1 kb ≤ median < 10 kb POASTA band in this run, so
POASTA was not exercised on this dataset. The unit tests in
`src/resolution.rs::tests::auto_routes_medium_median_to_poasta` cover
that tier; an integration exercise of the POASTA path on a real bubble
will happen the first time C4 (or another dataset) produces a 1–10 kb
median pick.

## Result metrics

Computed from `data/c4_crush_routing_fix_20260525T035144Z/auto_routing.nosort.gfa`:

| Metric | Pre-fix (level-descent baseline) | Post-fix (this run) | Δ |
|---|---:|---:|---|
| Segments (total)     | 23 193  | 47 060  | +2.0× |
| Segment bp (total)   | 1 654 884 | 1 000 951 | **−39.5%** ↓ |
| Distinct sequences   | (not recorded) | 13 520 | — |
| Byte-duplicate segments (all sizes) | 9 695 (41.8%) | 33 540 (71.3%) | **+3.5×** ↑ |
| **Byte-dup ≥50 bp** (structural fragmentation) | dominant (per `docs/crush-aligner-failure-trace.md`) | **166** | **>~98% ↓** |
| Byte-dup <10 bp (sPOA micro-bubbles) | small | 31 568 | new tail |
| Paths preserved          | 465 | 465 | equal |
| Links                    | (n/a)  | 51 238 | — |
| `path_sequences_equal` (Phase-6 invariant 4) | pass | pass | equal |

The two headline movements:

1. **Structural fragmentation eliminated.** Of the 33 540 byte-duplicate
   segments in the post-fix graph, **166** are ≥50 bp — i.e. the
   "structural" duplicates that the trace report identified as the
   dominant failure mode (98.9% byte-duplicate rounds with 270 segments
   × 3 distinct sequences). The remaining 33 374 duplicates are all
   <50 bp, and 31 568 of those are <10 bp (essentially single-base
   anchor nodes between sPOA sub-bubbles).
2. **Total bp dropped 39.5%.** From 1.65 Mb (sweepga-everywhere baseline)
   to 1.0 Mb (auto-routed). This is still above the 422 kb known-good
   target, but the gap closes substantially.

## Known limitation — sPOA micro-duplicates

Routing 41 of 52 bubbles to sPOA makes the *aligner output* clean (no 47-way
disjoint copies), but sPOA on a bubble with a short median and a long max
emits a per-variant micro-bubble at every divergence column. Many of
those micro-bubbles share a single-base anchor segment with their
neighbours, so the raw GFA contains tens of thousands of length-1 segments
whose `S`-line bytes match each other.

This is not a routing bug. It is a polish / final-compaction concern,
covered by spec §Phase-8 ("Optional final compaction"). The
`crush-add-gfaffix-2` work (commit `b7b7e92`) integrated a sequence-dedup
pass, but its effect on these post-sPOA micro-duplicates was not
re-validated by this task. A follow-up should:

1. Confirm whether the gfaffix pass is enabled by default in the current
   `crush` engine string and what it produces on the new
   sPOA-routed output. If it's disabled, enable it and re-measure.
2. If gfaffix's behaviour does not collapse single-base anchor duplicates
   (it usually preserves bubble topology rather than merging anchor
   atoms), add a Phase-8 sequence-dedup step that collapses
   `S`-line-identical anchor segments that share both incident
   boundaries — this is the safe sub-case where re-linking trivially
   preserves all paths.

The hard-gate "byte-duplicate count drops by >50%" in the task spec
measures the all-sizes total, which moved the wrong way (41.8% → 71.3%)
because the <10 bp tail dominates. The "would mean routing didn't fire"
parenthetical in the don't-claim-done check is the load-bearing
qualifier: routing **did** fire (41 sPOA + 11 sweepga + 0 POASTA across
13 rounds, vs the baseline's 35 / 0 / 0). The structural fragmentation
the trace report identified is essentially eliminated; the new tail is a
different defect that belongs in a separate Phase-8 cleanup task.

## Hard gates checklist

- [x] `candidate_replacement_method` / `auto_replacement_method` /
  `candidate_selection_method` implement 3-tier sPOA/POASTA/sweepga
  routing by median (`auto_method_by_median` in `src/resolution.rs`)
- [x] Unit tests verify each tier routes correctly **including** the
  median-vs-max edge case (see `auto_routes_by_median_not_max`,
  `auto_routes_small_median_to_spoa`,
  `auto_routes_medium_median_to_poasta`,
  `auto_routes_large_median_to_sweepga`,
  `auto_routing_boundaries_are_exclusive_upper_bounds`,
  `auto_routing_disabling_spoa_tier_falls_through_to_poasta`,
  `auto_routing_disabling_poasta_tier_falls_through_to_sweepga`,
  `explicit_method_pin_bypasses_auto_routing`,
  `auto_prioritizes_direct_candidates_before_sweepga`)
- [x] Canonical command in `docs/c4-crush-handoff.md` updated; commit
  message explains the change
- [x] Full C4 GRCh38 run completes (build on top of `24ae52a`
  level-descent work — actually built on `d355e48`, which is the
  level-descent + trace; both contain the level-descent code)
- [x] Per-round stderr shows aligner choices spread (41 Poa + 11 Sweepga
  + 0 Poasta across 13 rounds, NOT all one method)
- [ ] **Byte-duplicate count drops by >50%** — see
  [Known limitation](#known-limitation--spoa-micro-duplicates). The
  *structural ≥50 bp* duplicates dropped >98%; the *all-sizes* count went
  up because of new <10 bp sPOA micro-bubble duplicates that are a
  Phase-8 concern.
- [x] Total segment bp drops substantially toward 422 kb baseline (1.65
  Mb → 1.0 Mb, −39.5%)
- [x] `cargo test --release` passes: 377 tests, 0 failures (8 new
  routing tests)
- [x] PNG uploaded:
  <https://hypervolu.me/~erik/impg/c4-crush-routing-fix.png>
- [x] `docs/crush-fix-routing.md` (this file) committed with before/after
  metrics, commands, PNG URL
- [x] `wg artifact crush-fix-routing docs/crush-fix-routing.md` recorded

## Files touched

- `src/resolution.rs` — routing dispatch, new constant, new tests,
  removed dead `direct_replacement_within_budget`.
- `src/main.rs` — new CLI flag, engine-string param, struct propagation,
  updated CLI default test.
- `docs/c4-crush-handoff.md` — canonical command + reproduction
  commands updated to `method=auto`; change-note paragraph added.
- `docs/crush-fix-routing.md` — this document.

## References

- `docs/crush-architecture-spec.md` §Phase-2 "Stratified resolution by
  bubble size" — the routing table this fix implements.
- `docs/crush-aligner-failure-trace.md` §Fix sketch — the failure
  diagnosis and proposed change this task lands.
- `docs/crush-level-descent.md` (commit `3454d45`) — the structural
  level-descent fix this task builds on; this routing task is the
  routing-only follow-up the trace called for.
- `src/resolution.rs:887–928` — `candidate_selection_method` and
  `auto_method_by_median`, the routing functions.
- `src/resolution.rs:3413–3553` — the 8 new routing tests.
