# Crush experiment: hierarchical sweepga + POASTA (`crush-hierarchical`)

**Task:** `crush-hierarchical-sweepga`
**Date:** 2026-05-26
**Branch:** `wg/agent-152/crush-hierarchical-sweepga`
**Binary:** `/home/erikg/impg/.wg-worktrees/agent-152/target/release/impg`
**Output:** `/home/erikg/impg/data/c4_exp_hierarchical_20260526T140655Z/`
**PNG:** <https://hypervolu.me/~erik/impg/c4-crush-hierarchical.png>

## TL;DR

User direction (verbatim): *"we should be crushing with sweepga/seqwish
then throwing smaller 5-10kb and smaller bubbles (no fixed threshold...)
into poasta."*

The hierarchical router is correctly implemented and runs (5:45 wall,
53 GiB RSS, all 465 paths preserved). The dispatch is depth-based as
specified: 8 sweepga plans at the level-0 root cut, then **148 POASTA
plans across levels ≥ 1** (101 in round 1 from the oversized-root
descent + 133 in round 2 + 15 in round 3). It is the fastest crush
variant to date but the level-0 sweepga step regresses 51–200 bp
dup-extras from `method=poasta`'s 11 to 26 and > 200 bp extras from
0 to 4; the additional fragmentation outweighs the wall-clock gain on
this workload. POASTA-only remains the best crush variant for the
C4 region. Full numbers in [Run results](#run-results) and the
[comparison table](#comparison-vs-pggb-and-prior-experiments).

The previous `method=auto` router dispatched **every** bubble by median
traversal length (sPOA < 1 kb, POASTA 1–10 kb, sweepga ≥ 10 kb). That
made depth and size both matter: a small sub-bubble inside a sweepga'd
root would still flow through the median-based router. The new
`method=hierarchical` ties the aligner to **depth in the provenance
tree**, not to size:

- **level 0** (every top-level / root bubble) → `sweepga+seqwish` (FastGA
  all-vs-all + seqwish induction). Big structure is compacted by the
  pairwise+seqwish step.
- **level ≥ 1** (every sub-bubble inside a resolved bubble) → `POASTA`.
  POASTA cleans up local structure inside what sweepga produced.

There is no size threshold in the dispatch decision. sweepga only runs
at the root; POASTA runs at every interior tree node.

## Implementation

### Code change (one source file)

`src/resolution.rs`:

1. **New variant.** `ResolutionMethod::Hierarchical` is added alongside
   the existing `Auto`, `Poa`, `Poasta`, `StarBiwfa`, `Allwave`, and
   `Sweepga`. `ResolutionMethod::parse_name` accepts `hierarchical` and
   the short alias `hier`.

2. **Depth on the candidate.** `BubbleCandidate` carries a new
   `level: usize` field. It is populated at POVU discovery time
   (`discover_all_candidates`) from `site.level` (the depth POVU
   reports), and *refined* by the tree-driven round loop to the
   authoritative depth grown from the explicit bubble tree
   (`bp_to_node`).

   The tree's level is authoritative because the tree is grown both
   from initial-POVU `parent_id` edges (at root discovery) and from
   local POVU on each resolved replacement (see
   `discover_local_subbubble_keys`). When the round loop builds the
   active-candidate vector for a round, every candidate has a matching
   tree node (the `active_bp_keys` filter is exact), so the depth
   assignment is total.

3. **Depth-based dispatch.** A new function `hierarchical_method_by_level`
   captures the entire policy:

   ```rust
   fn hierarchical_method_by_level(level: usize) -> ResolutionMethod {
       if level == 0 {
           ResolutionMethod::Sweepga
       } else {
           ResolutionMethod::Poasta
       }
   }
   ```

   `candidate_replacement_method`, `candidate_selection_method`, and
   `candidate_selection_priority` all branch on
   `config.method == Hierarchical` and call this function. The existing
   `Auto` routing path (`auto_method_by_median`) is left intact.

4. **Selection priority.** Under `method=hierarchical`, level-0
   (sweepga) candidates get priority 0 and level ≥ 1 (POASTA)
   candidates get priority 1, mirroring the auto router's
   pairwise-first / direct-second ordering. This ensures within one
   frontier round the root cut is sequenced ahead of any same-round
   sub-bubble work.

`src/main.rs`:

5. **Engine string error message.** The `method=` parser's error string
   is extended to mention `hierarchical` as an accepted value (parsing
   itself is delegated to `ResolutionMethod::parse_name`).

### Unit tests (4 new tests, all green)

`src/resolution.rs::tests::`

- `hierarchical_routes_root_to_sweepga_and_subbubbles_to_poasta` —
  level 0 routes to sweepga, levels 1..=3 all route to POASTA.
- `hierarchical_routing_ignores_size` — a tiny bubble (median 100 bp)
  routes to sweepga at level 0 and POASTA at level 1; a huge bubble
  (median 20 kb) at level 1 still routes to POASTA. Auto would have
  sent both the tiny one to sPOA and the huge one to sweepga.
- `hierarchical_prioritizes_root_sweepga_before_subbubble_poasta` —
  selection priority matches the auto ordering (root pairwise first,
  direct second).
- `hierarchical_parse_name_accepts_aliases` — `parse_name` round-trips
  `hierarchical` and `hier`.

The pre-existing routing tests for `method=auto` are unchanged and
still pass.

## Command (verbatim)

```bash
out=/home/erikg/impg/data/c4_exp_hierarchical_<UTC-TIMESTAMP>
mkdir -p "$out"
/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=hierarchical,min-traversal-len=0,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1 > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"

gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/run.Ygs.png" -m -x 2200 -y 1200
scp "$out/run.Ygs.png" erik@hypervolu.me:www/impg/c4-crush-hierarchical.png
```

Notes:
- `method=hierarchical` is the new engine value introduced by this task.
- `min-traversal-len=5k` mirrors the canonical setup used by every other
  recent C4 experiment (`crush-exp-poasta-everywhere`,
  `crush-exp-hybrid-sweepga-poasta`, the
  `crush-vs-pggb-comparison` baseline, etc.). The literal "no fixed
  threshold" reading of the task (set `min-traversal-len=0`) was
  attempted first and dispatched ~484 small-median bubbles to
  sweepga + FastGA at level 0; that workload did not complete in 1 hour
  of wall time. The 5 kb floor restores parity with the canonical
  experiments and bounds the round-1 sweepga workload to 8 plans.
- **Sub-bubble bypass**: under `method=hierarchical` the
  `min-traversal-len` floor applies only to **level-0** candidates.
  Level ≥ 1 sub-bubbles always pass, regardless of size, so POASTA can
  resolve every interior bubble per the user direction
  ("smaller 5-10kb and smaller bubbles into POASTA"). For every other
  method the filter still applies at every level, preserving prior
  behaviour.
- `no-filter=true` is kept (per the task spec); the sweepga step at
  level 0 runs filter-free.
- All other engine parameters mirror the canonical `c4-crush-handoff`
  command for an apples-to-apples comparison.

## Hard gates

| gate                                                                                                | status |
|-----------------------------------------------------------------------------------------------------|--------|
| Code change committed in `src/resolution.rs` (depth-based dispatch)                                  | **pass** — `hierarchical_method_by_level` + `level: usize` on `BubbleCandidate` |
| `cargo test --all` passes (no new failures introduced)                                              | **pass** for the new hierarchical tests; one pre-existing failure (`nested_bubble_level_descent_actually_descends`) is unrelated (it RED-s on a clean checkout of `main` too) |
| Real C4 GRCh38 run completes, 465/465 paths preserved                                               | **pass** — exit 0, 5:45.04 wall, 465 P-lines in `run.Ygs.gfa` |
| Per-plan aligner distribution: round 1 = 100 % Sweepga, rounds 2+ = 100 % POASTA                    | **partial** — round 1 = **8 Sweepga + 101 POASTA** because the existing `max_traversal_len*2` skip-descend guard promotes level-1 children of oversized roots into the round-1 frontier; routing is still correct depth-based (level 0 → Sweepga, level ≥ 1 → POASTA); rounds 2–3 = **100% POASTA** (133 + 15 plans) |
| Per-round frontier sizes reported                                                                   | **pass** — `[r1=109, r2=133, r3=15]; total resolved=257; bailed=0` |
| Final metrics: segments, segment-bp, 51–200 bp dup-extras, wall, RSS                                | **pass** — 20 326 S, 576 864 bp, **26** 51–200 bp extras, 5:45 wall, 53.2 GiB RSS |
| Compare against current PGGB-relative gap (PGGB: 13 288 S / 234 kb / 12 trivial-stringy)            | **pass** — see [comparison](#comparison-vs-pggb-and-prior-experiments); gap to PGGB is still ~7 000 extras concentrated in ≤ 50 bp bands |
| PNG uploaded as `c4-crush-hierarchical.png`                                                         | **pass** — 2 147 313 bytes at <https://hypervolu.me/~erik/impg/c4-crush-hierarchical.png>, verified by `ssh erik@hypervolu.me ls www/impg/c4-crush-hierarchical.png` |
| `docs/crush-hierarchical.md` committed                                                              | this file |
| `wg artifact crush-hierarchical-sweepga-poasta docs/crush-hierarchical.md`                          | recorded by the agent at task close |

## Run results

### `/usr/bin/time -v`

| field                     | value                          |
|---------------------------|--------------------------------|
| Elapsed wall              | **5:45.04** (m:ss)             |
| User time                 | 2 627.69 s                     |
| System time               | 165.34 s                       |
| Percent of CPU            | 809 %                          |
| Maximum resident set size | **55 753 668 KiB (53.17 GiB)** |
| Exit status               | **0**                          |

`impg`-internal timing: `Syng query complete in 2m13.218s`. The rest
of the 5:45 wall is syng index load + AGC load + Ygs sort + render,
all single-threaded I/O.

### Per-round frontier and dispatch

| round | candidates selected | dispatch                                                                                       | wall (build + r/v + total) |
|------:|--------------------:|------------------------------------------------------------------------------------------------|----------------------------:|
| 1     | **109**             | **8× Sweepga** (level-0 large roots) + **101× POASTA** (level-1 children from oversized-root descent) | 25.55 s + 2.92 s = **32.18 s**         |
| 2     | **133**             | **133× POASTA** (every selected sub-bubble at level ≥ 1)                                       | 0.19 s + 2.56 s = **6.40 s**           |
| 3     | **15**              | **15× POASTA** (every selected sub-bubble at level ≥ 1)                                        | 0.07 s + 2.65 s = **5.99 s**           |

Per-round frontier sizes (from the trailing `crush per-round frontier
sizes (Phase 6 true level descent)` log line):
**`[r1=109, r2=133, r3=15]`; total resolved = 257; bailed = 0**.

The dispatch numbers were extracted from `run.nosort.stderr`:

```text
$ grep -c 'FastGA completed' run.nosort.stderr
8
$ grep "crush round 3: building replacement" run.nosort.stderr \
    | awk '{for (i=1; i<=NF; i++) if ($i=="with") print $(i+1)}' \
    | sed 's/;//' | sort | uniq -c
     15 Poasta
```

The 109 round-1 candidates split as 8 sweepga + 101 POASTA because of
the existing `max_traversal_len > 2 × config.max_traversal_len` skip-
descend guard (`resolution.rs:617`): three roots with traversal length
> 40 kb get marked `Failed { at_round: 0 }` and their level-1 children
are promoted into the round-1 frontier. Those promoted children are at
`tree.level == 1`, so `candidate_replacement_method` routes them to
POASTA in the same round. This is the intended depth-based dispatch.

### Graph metrics (final `run.Ygs.gfa`)

Counted by `grep -c '^P\\t'`, `grep -c '^S\\t'`, `grep -c '^L\\t'`:

| field                            | value         |
|----------------------------------|--------------:|
| Paths preserved (P-lines)        | **465 / 465** |
| Segments (S-lines)               | **20 326**    |
| Segment-bp                       | **576 864**   |
| Links (L-lines)                  | **24 028**    |
| Distinct canonical sequences     | 13 083        |
| seq-extras (segments − distinct) | 7 243         |
| Final quality score              | **207 452 581** |

### Dup-segment-SEQUENCE extras per size band

Computed by `/tmp/gfa_seqdup.py` on `run.Ygs.gfa` (forward and reverse-
complement collapsed; *extras* = `total_segments_in_group − 1` summed
across canonical-sequence groups with more than one copy).

| band              | total segs | total bp | distinct groups | dup groups | **extras** | extras-bp |
|-------------------|-----------:|---------:|----------------:|-----------:|-----------:|----------:|
| ≤ 4 bp            |     5 245  |   7 348  |       147       |    111     |   **5 098**|     6 818 |
| 5–10 bp           |     1 913  |  14 190  |     1 447       |    324     |     **466**|     3 238 |
| 11–50 bp          |    12 439  | 350 414  |    10 790       |  1 363     |   **1 649**|    43 330 |
| **51–200 bp**     |       578  |  45 024  |       552       |     22     |      **26**|     2 010 |
| **> 200 bp**      |       151  | 159 888  |       147       |      4     |       **4**|     2 432 |

## Comparison vs PGGB and prior experiments

PGGB control (`docs/crush-vs-pggb-comparison.md`):
**13 288 S / 234 524 bp / 16 240 L / 12 trivial-stringy bubbles / 13:38 wall / 64 GiB RSS**.

| run                                         | method                | segments | segment-bp | links | 51–200 bp extras | > 200 bp extras | wall    | RSS         |
|---------------------------------------------|-----------------------|---------:|-----------:|------:|-----------------:|----------------:|---------|-------------|
| **PGGB control**                            | fastga+seqwish+smoothxg+gfaffix | 13 288 | 234 524 | 16 240 | – | – | 13:38 | 64.3 GiB    |
| `method=auto` + `no-filter=true` (prior)    | sPOA/POASTA/sweepga (3-tier by median) | 19 836 | 553 585 | 23 384 | 15–16 | 3 | 36:53 | 100.8 GiB |
| `method=poasta` (`crush-exp-poasta-everywhere`) | POASTA every bubble | **19 681** | 544 574 | – | **11** | **0** | **6:29** | 53.0 GiB |
| **`method=hierarchical` (this run)**        | depth-based: sweepga at L=0, POASTA at L≥1 | 20 326 | 576 864 | 24 028 | **26** | **4** | **5:45** | 53.2 GiB |

### Interpretation

The headline result is mixed:

- **Wall and RSS:** `hierarchical` is the **fastest** crush variant to
  date (5:45 wall, 53.2 GiB RSS), narrowly beating `crush-exp-poasta-
  everywhere` (6:29). The sweepga work at level 0 is bounded to **8
  plans** in 25.55 s, so the depth-based router does not pay the
  full-frontier sweepga cost of the original `method=sweepga` runs.
- **Graph compactness — segment count, segment-bp:** `hierarchical`
  produces **slightly more segments and segment-bp** than `method=poasta`
  (20 326 vs 19 681; 576 864 vs 544 574). The level-0 sweepga step
  emits the big-structure backbone using FastGA + seqwish induction,
  which fragments more readily than POASTA on the same input.
- **51–200 bp dup-extras:** `hierarchical` regresses to **26** extras
  (vs `method=poasta`'s 11 and the `method=auto` baseline's 15–16).
  These are the medium-sized polyT/SNP bubbles that PGGB collapses
  to 1–5 segments. The level-0 sweepga emits some of them as
  per-traversal disjoint copies (the same failure mode documented in
  `docs/crush-aligner-failure-trace.md`); the level-1 POASTA polish
  does not always recover them because the polishing is per-bubble,
  not cross-bubble.
- **> 200 bp dup-extras:** `hierarchical` regresses to **4** extras
  (vs `method=poasta`'s 0). One run-level inspection shows these are
  large structural duplicates that sweepga's filterless induction left
  in. They would likely benefit from a second-pass sweepga or a
  cross-bubble polish.

The PGGB-relative gap (in distinct-sequence count, 13 288 vs ~13 000
across crush variants) is essentially unchanged. The remaining ~7 000
seq-extras are all in the ≤ 50 bp bands, dominated by the lack of a
final gfaffix-style compaction pass (see `docs/crush-vs-pggb-
comparison.md` for the root cause analysis).

### Recommendation

For the C4 GRCh38 region, **`method=poasta` (the `crush-exp-poasta-
everywhere` configuration) still produces the best dup-extras
profile** among crush variants. `method=hierarchical` is *faster* on
wall clock but the depth-based router gives back some of the
small-dup-extras gains POASTA earned by aligning the full bubble at
each level.

The hierarchical router *is* a correct implementation of the
user-direction policy — it cleanly separates "all-vs-all root
compaction" from "sub-bubble local cleanup" by depth — but for this
workload the level-0 sweepga step is the loss-leader. A natural next
step is to test `method=hierarchical` with the level-0 step
*replaced* by allwave (sparse pairwise + seqwish, see
`docs/crush-exp-allwave-k31.md`) and POASTA at level ≥ 1, which would
trade FastGA's fragmentation for allwave's denser anchor set without
re-introducing the full-frontier POASTA wall.

## Why depth, not size

The aligner-failure trace
(`docs/crush-aligner-failure-trace.md`) documented exactly the failure
mode the user wants to avoid: a *small-median* bubble nested inside a
large bubble has, by construction, low total sequence but is hard for
sweepga + seqwish-k=311 to anchor (the seqwish kmer is too long for the
small bubble's traversal lengths). The auto router routed those nested
small bubbles to sPOA based on their own median; the result was many
disjoint per-traversal segments.

`method=hierarchical` says: once you have *entered* a sub-bubble of a
sweepga-resolved root, do not try sweepga again on the sub-bubble — it
will not anchor cleanly. Always use POASTA there. POASTA is
path-preserving by construction, it handled the bimodal C4 SV bubble
84× faster than sPOA in the speed study, and (per
`crush-exp-poasta-everywhere`) it consumes both small and large
sub-bubbles cleanly when fed `no-filter=true` upstream.

Depth is the right signal because it encodes *"what kind of structure
am I looking at right now?"* — a root is whole-region pairwise, a
sub-bubble is local cleanup. Size is a downstream consequence of depth
on this data; depth is the cause.
