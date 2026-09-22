# Nested-bubble level descent regression test

This document captures the design and the HEAD-failure evidence for the
regression test `nested_bubble_level_descent_actually_descends` in
`tests/test_crush_integration.rs`.

The test pins down the user-visible behavior shown in
`/tmp/crush-ld/c4-crush-level-descent.png` and described by the user:

> "There's this ramp in the underline. It's incredibly clear. ... Are we
> getting a new bubble tree and putting it inside the bubble we've crashed?
> And resolving that too?"

The white-space "ramp" is the visible signature of crush re-POVU'ing the
WHOLE working graph each round and selecting LEAVES (via the `site.is_leaf`
filter at `src/resolution.rs:1224`) instead of resolving each nested-bubble
parent first and re-POVU'ing on its LOCAL replacement subgraph.

## Fixture: tests/test_data/crush/nested_bubbles_real.gfa

Real C4A pangenome extract — NOT synthetic.

- **Provenance.** Sliced from `tests/test_data/crush/c4_slice_1500_3000.gfa`
  (itself a real C4A blunt-pangenome extract:
  `CHM13#0#chr6:31744284-31976975`, 465 haplotypes, seqwish-k=311, derived
  from `/home/erikg/impg/data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.gfa`).
- **Selection range.** Reference path step indices 568-590 of the CHM13
  reference walk (a ≈22-step window around the POVU site
  `<43388736>43388746`).
- **Paths kept (5).** CHM13#0#chr6 plus four haplotypes selected to cover
  five DISTINCT path signatures through the parent bubble, so that POVU
  exposes the nested structure on the slim subgraph:
  - `CHM13#0#chr6:31744284-31976975` (reference allele)
  - `HG03927#1#CM086636.1:31887369-32120054`
  - `HG02622#1#JAHAOO020000001.1:32050104-32302570`
  - `NA18945#1#CM101613.1:31897691-32124049`
  - `HG02647#1#CM086559.1:31866435-32092803`
- **Subgraph extraction.** All segments visited by the 5 chosen path slices
  were kept; links between any two retained segments were kept; sequences
  copied verbatim from the source GFA (real biological bp, no edits).
- **Size on disk.** 3 857 bytes (3.8 KB), 40 segments, 981 bp of segment
  sequence, 43 links, 5 paths. One connected component.
- **POVU on the fixture (reference `CHM13#0#chr6:31744284-31976975`).** 3
  sites, 2 leaves, levels `[L0: 1, L1: 2]`:
  - `<43388736>43388746` — L0, leaf=false, span 16 ref-steps (the
    TOP-LEVEL bubble)
  - `>43388742<2094988` — L1, leaf=true, span 4 ref-steps, parent
    `<43388736>43388746` (nested sub-bubble 1)
  - `>43388742>546635` — L1, leaf=true, span 3 ref-steps, parent
    `<43388736>43388746` (nested sub-bubble 2)

This matches the task contract: exactly one top-level bubble containing
≥2 nested sub-bubbles, real-data biological sequence, small enough to
`cat` or `ls`.

## What the test asserts

`nested_bubble_level_descent_actually_descends` runs `resolve_gfa_bubbles`
on the fixture with `max_iterations=5, method=Auto, …default`, then also
spawns the `impg crush` CLI with `-v 2` to capture per-round counts. It
checks:

1. **Fixture sanity.** The input has exactly 1 L0 site and ≥2 L1 sites
   nested inside that L0, and 1 connected component.
2. **Termination round count.** `stats.iterations <= 2`. The correct
   top-down algorithm needs exactly two rounds: round 1 resolves the L0
   parent, round 2 resolves the L1 children discovered by re-POVU on the
   replacement's LOCAL subgraph; round 3 has no eligible candidates and
   terminates. ← **Primary discriminator.**
3. **No fragmentation.** A fresh POVU pass on the resolved graph finds
   ≤ 2 sites total. HEAD's repeated whole-graph POA-fragmentation
   re-discovery leaves ~15 sites (10 L0 + 5 L1). ← **Quantitative ramp
   check.**
4. **No new connected components.** `metrics.components` is unchanged
   from input.
5. **Path sequence preservation.** Every input path's sequence appears
   verbatim in the resolved graph (sanity gate that should pass on both
   HEAD and any fix).
6. **Per-round selected counts (parsed from CLI stderr `-v 2`).** Round 1
   should select exactly 1 candidate (the L0 parent); round 2 should
   select ≥ 2 candidates (the nested L1 leaves).

## Behavior on HEAD (471f089) — captured failure output

`HEAD = 471f089a1f4eef355eaa84e1459365a2b7471c3b` ("feat: crush-fix-routing (agent-92)")

Running the test reproduces this failure:

```
$ source ./env.sh && cargo test --release --test test_crush_integration \
    nested_bubble_level_descent_actually_descends -- --nocapture

running 1 test
fixture POVU: 3 sites total, level breakdown L0=1, L1=2 (L1 nested-in-L0=2)
library crush: iterations=4, resolved=6, bailed=0, candidates_seen=6

thread 'nested_bubble_level_descent_actually_descends' panicked at
  tests/test_crush_integration.rs:675:
  nested-bubble level descent should terminate in <= 2 iterations
  (correct top-down algorithm needs exactly 2: round 1 = L0 parent,
   round 2 = nested L1 children); got iterations=4 on this run
  (HEAD 471f089 takes 5 rounds because re-POVU on the WHOLE graph +
   is_leaf filter at src/resolution.rs:1224 keeps re-discovering
   POA fragmentation artifacts each round)
test nested_bubble_level_descent_actually_descends ... FAILED
```

The library run reports `iterations=4` (with the default `max_iterations=5`
cap applied; the natural round count is 5 — see the CLI section below).
Either way `iterations > 2`, so the discriminator assertion fires.

### Observed CLI trace (HEAD)

`impg crush --gfa nested_bubbles_real.gfa --method auto --max-iterations 5 -v 2`:

| round | POVU sites | selected (frontier) | notes                                                       |
|------:|-----------:|--------------------:|-------------------------------------------------------------|
|     1 |          3 |                   1 | leaf with median-len 406, NOT the L0 parent                  |
|     2 |         16 |                   5 | re-POVU on whole rewritten graph: 5 new replacement targets |
|     3 |         16 |                   6 | further whole-graph fragmentation re-detected               |
|     4 |         16 |                   4 |                                                             |
|     5 |         16 |                   1 |                                                             |
|     6 |         15 |                   0 | no eligible candidates → terminate                          |

Final result: 16 replacements made across 5 rounds; POVU on the resolved
graph still finds 15 sites (10 L0 + 5 L1) — the fragmentation signature
that the test pins down via the `result_povu.sites <= 2` assertion.

### What the corrected behavior should look like

| round | POVU sites | selected (frontier) | notes                                                                   |
|------:|-----------:|--------------------:|-------------------------------------------------------------------------|
|     1 |          3 |                   1 | the L0 parent `<43388736>43388746` is resolved                          |
|     2 |          ≥2 (LOCAL re-POVU)               |                   ≥2 | re-POVU on the L0 replacement subgraph exposes ≥2 nested L1 leaves         |
|     3 |          0 |                   0 | no eligible candidates → terminate, iterations==2                       |

POVU on the resolved graph should find ≤ 2 sites (cleanly consumed).

## How to reproduce

```bash
source /home/erikg/impg/env.sh   # Guix GCC + glibc 2.35
git checkout 471f089             # current HEAD
cargo build --release --test test_crush_integration
cargo test --release --test test_crush_integration \
    nested_bubble_level_descent_actually_descends -- --nocapture
```

Expected: the test PANICS with the `iterations <= 2` message shown above.

## Hand-off

This test is the regression gate for `crush-true-level-descent` (the
downstream task that fixes the algorithm). Any candidate fix must keep
the assertions in this test green AND keep the existing
`phase6_level_descent_terminates_below_max_iterations` and
`c4_slice_auto_crush_preserves_path_sequences` tests green.
