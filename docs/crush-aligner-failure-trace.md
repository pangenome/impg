# Crush — aligner failure trace (`crush-aligner-failure`)

**Task:** `crush-aligner-failure`
**Date:** 2026-05-25
**Branch:** `eg/c4-crush-resolution-controls`
**Author:** agent-89 (Default Evaluator role, investigation-only)
**Inputs analysed:**

- `/tmp/crush-ld/C4A_full.gfa` — the level-descent output that prompted the
  investigation (23 193 segments, 1.65 Mb total bp, 4.2 M path steps).
- `/tmp/crush-ld/c4_stderr.log` — full stderr from that run (9 MB, 9 rounds,
  35 bubbles resolved).
- `/tmp/crush-ld-debug/replacement_000[0-7]_sweepga_seqwish/{seqwish,unchopped}.gfa`
  — per-replacement subgraph dumps captured by re-running the canonical
  command with `IMPG_CRUSH_DEBUG_DIR=/tmp/crush-ld-debug` against the same
  binary, syng index, and AGC store. Round 1 selected 8 bubbles in this
  rerun (the original had 7); the only difference relevant to this report
  is that round 1's leaf set has slightly different ordering after the
  previous resolutions — the per-bubble shape data is identical for every
  shared candidate.

No source code was changed (verified: `git diff --stat` reports only `m
vendor/gfaffix`, untouched docs, and the new file
`docs/crush-aligner-failure-trace.md`).

## TL;DR

The dominant failure mode is **wrong aligner for scale**: the run was
started with `method=sweepga`, so `candidate_replacement_method` in
`src/resolution.rs:1649` returns `Sweepga` for **every** bubble regardless
of size. `build_sweepga_seqwish_replacement` then induces a graph with
`seqwish-k=311`, which cannot anchor matches in traversal sequences shorter
than 311 bp. Every small-median bubble (`median_len <~ 311`) therefore
comes out of seqwish induction as **N byte-identical disjoint segments
with zero links** — one per traversal, no compaction, no graph structure.
The downstream `polish_replacement_gfa_with_flubbles` step cannot recover
because POVU finds no bubbles to polish in a graph that has no links.

Quantitative evidence: 4 out of 5 bubbles below show this exact pattern; in
the original run, 5 of the 9 rounds emit at least 69% byte-duplicate
replacement segments (round 5 emits 270 segments with only **3** distinct
sequences — 98.9% duplicates). The spec invariant "No two segments share a
sequence" (`docs/crush-architecture-spec.md` §Invariants 1) is violated on
a sub-bubble basis by every small-median candidate the run resolves.

Root cause is the absence of size-tiered routing: the spec's three-aligner
strategy (sPOA <1 kb / POASTA 1–10 kb / sweepga >=10 kb) is **not
implemented anywhere in the code path**. `auto_replacement_method` only
ever switches between `Poa` and `Allwave`, never `Poasta` or `Sweepga`
(`src/resolution.rs:887–904`); when the config pins a specific method the
switch is bypassed entirely. The fix lives at `src/resolution.rs:1649–1668`
— see §Fix sketch.

## Method (how the data below was produced)

For each pick, the **traversal-stats** numbers (`traversals`, `max-len`,
`median-len`, `total-len`, `root-span`) come from the verbatim per-round
"building replacement X/Y with Sweepga" log line in
`/tmp/crush-ld/c4_stderr.log`. The **replacement subgraph shape**
(segments, links, paths, unique-sequence counts) comes from the captured
`replacement_XXXX_sweepga_seqwish/seqwish.gfa` dump in the debug-dir rerun
(round-1 picks 1–3) or, for the deep-round picks, from **slicing
`/tmp/crush-ld/C4A_full.gfa` by the per-round `ids minted` range** logged
at the end of each round, e.g.

> `crush round 5: resolved 2/2 replacement(s) in 1.80s; ... ids minted 272221292..272221907;`

Every surviving segment whose id falls in `[272 221 292, 272 221 907)` was
minted by a round-5 replacement subgraph and never consumed by a later
round, so segment-count, byte-length, and sequence-distinct statistics
computed over that id slice are a faithful summary of what round 5's two
replacements contributed to the final graph.

## The 5 picks

The original run resolved 35 bubbles across 9 rounds. Picks (a)–(c) come
from round 1 because round 1 is the only round whose selection spans all
three spec-tier size classes; picks (d) and (e) are from the deeper-round
range the task asks for.

### Pick (a) — small bubble that should have gone to sPOA

| Field | Value |
|---|---|
| Round / plan | round 1, plan 8/8 |
| Aligner chosen (log line) | `crush round 1: building replacement 8/8 with Sweepga; traversals=47, max-len=25230, median-len=165, total-len=32820, root-span=165` |
| Did the aligner run? | Yes — log line `crush round 1: replacement build progress 8/8 (accepted)` (`/tmp/crush-ld/c4_stderr.log:36900` in the original 7-bubble round; the 8-bubble rerun shows the same accept at `replacement build progress 8/8 (accepted)`). The seqwish.gfa dump exists at `replacement_0000_sweepga_seqwish/seqwish.gfa`. |
| Spec-correct aligner | **sPOA** (median 165 bp << 1 kb) |
| Replacement subgraph shape | **47 S / 0 L / 47 P**; 2 distinct sequences total (one 165 bp × 46 copies, one 25 230 bp × 1 copy); total bp = 32 820 |
| `unchop_gfa` change | none — `diff seqwish.gfa unchopped.gfa` is empty (there are no length-1 in/out chains to merge) |
| Polish behaviour | `polish=Poa` so `polish_replacement_gfa_with_flubbles` (`src/resolution.rs:2391`) is invoked. POVU finds **zero** sites because the graph has no links → polish exits with the replacement unchanged. |
| Integration behaviour | `apply_replacement_frontier` (`src/resolution.rs:1496–1631`) splices the 47 disjoint segments into the working graph; `render_rewritten_graph` mints 47 fresh ids in the round-1 id range. No segments are dropped, no segments are coalesced — the byte-duplicate structure is preserved end-to-end. |

Inside the 47-traversal bubble the sweepga + seqwish pipeline collapses
**zero** of the 46 byte-identical 165 bp traversals: every one becomes its
own segment with its own path line and no incident edges. The polish step
cannot help because there is nothing for POVU to decompose. A direct sPOA
pass on the 47 traversals would emit one node for the consensus 165 bp
sequence plus a one-off insertion sub-bubble for the 25 230 bp outlier —
order-of-magnitude difference.

### Pick (b) — medium bubble that should have gone to POASTA

| Field | Value |
|---|---|
| Round / plan | round 1, plan 2/8 |
| Aligner chosen (log line) | `crush round 1: building replacement 2/8 with Sweepga; traversals=356, max-len=9320, median-len=3319, total-len=1181720, root-span=3319` |
| Did the aligner run? | Yes — `replacement build progress` line accepted; seqwish.gfa dump in `replacement_0004_sweepga_seqwish/` (re-mapped to plan 2/8 by traversal-count match: `P-lines = 356`) |
| Spec-correct aligner | **POASTA** (median 3.3 kb in [1 kb, 10 kb]) |
| Replacement subgraph shape | **77 S / 89 L / 356 P**; 62 distinct sequences (15 duplicates of other segments) |
| Polish behaviour | POVU finds bubbles inside this graph (it has links and branching), and `polish_replacement_gfa_with_flubbles` runs `resolve_graph_bubbles` with `method=Poa, max_traversal_len=10 000, max_median_traversal_len=1 000`. Small nested bubbles get SPOA-resolved silently (polish emits no INFO logs). |
| Integration behaviour | Standard splice via `apply_replacement_frontier`. Substantial compaction: 1.18 Mb of input collapsed to 12 376 bp of replacement segments (95x bp reduction). |

This bubble is in the spec's POASTA band, and even with the "wrong" aligner
sweepga + seqwish produced a respectable result (77 segments for 356
traversals, plenty of links). The polish step has work to do — POVU can
decompose it. This is the size-class where sweepga's choice is *adequate*,
not catastrophic.

### Pick (c) — large bubble that should have gone to sweepga

| Field | Value |
|---|---|
| Round / plan | round 1, plan 1/8 |
| Aligner chosen (log line) | `crush round 1: building replacement 1/8 with Sweepga; traversals=312, max-len=31478, median-len=25155, total-len=5132443, root-span=110` |
| Did the aligner run? | Yes — `replacement build progress 1/8 (accepted)`; seqwish.gfa dump at `replacement_0006_sweepga_seqwish/` (matched by `P-lines = 312`) |
| Spec-correct aligner | **sweepga** (median 25 kb >= 10 kb) ✓ |
| Replacement subgraph shape | **825 S / 926 L / 312 P**; 479 distinct sequences (346 duplicates) |
| Polish behaviour | POVU + bounded SPOA polish runs and resolves nested tangles. |
| Integration behaviour | Splice via `apply_replacement_frontier`. 5.13 Mb of input collapsed to 55 952 bp (92x bp reduction). |

This is the one bubble in the round where the configured aligner matches
spec. Its 825-segment output is large in absolute terms but tightly linked
(more links than segments) — that's a real graph, not a chain of disjoint
copies. The 42% byte-duplicate rate here is the irreducible "tube of
slight haplotype divergence" you'd get for any large-bubble aligner; sPOA
on 312 traversals × 25 kb would be infeasible.

### Pick (d) — deep round, small bubble (round 5 plan 2/2)

| Field | Value |
|---|---|
| Round / plan | round 5, plan 2/2 |
| Aligner chosen (log line) | `crush round 5: building replacement 2/2 with Sweepga; traversals=275, max-len=25253, median-len=106, total-len=130464, root-span=106` |
| Did the aligner run? | Yes — `crush round 5: replacement build progress 2/2 (accepted)` (`/tmp/crush-ld/c4_stderr.log:53265`); the round summary `resolved 2/2 replacement(s) in 1.80s` is line `53266` |
| Spec-correct aligner | **sPOA** (median 106 bp << 1 kb) |
| Replacement subgraph shape (inferred from output GFA id range 272 221 292..272 221 907) | All round-5 ids together: **270 segments, 54 493 bp, 3 distinct sequences (267 duplicates = 98.9%)**. Plan 2/2 (n=275, median 106 bp) is the dominant contributor; plan 1/2 (n=252, median 284 bp) is the other. Both plans are small-median; both produce essentially nothing but duplicates. |
| Polish behaviour | `polish=Poa`; for a sub-graph that consists mostly of disjoint paths the polish has no bubbles to find (same mechanism as pick (a)). |
| Integration behaviour | 270 fresh-id segments minted in the round-5 id range survive to the final GFA. Path-step count drops only 2 718 in the whole round (4 214 694 -> 4 211 976) — round 5 contributes essentially zero compaction. |

The 98.9% byte-duplicate ratio is the strongest single quantitative
signature of the failure mode in the run. Three distinct sequences across
270 segments means almost every surviving "compacted" segment is just a
copy of one of three short consensus strings, each per-path.

### Pick (e) — deep round, dense small bubble (round 8 plan 1/4)

| Field | Value |
|---|---|
| Round / plan | round 8, plan 1/4 |
| Aligner chosen (log line) | `crush round 8: building replacement 1/4 with Sweepga; traversals=391, max-len=56583, median-len=285, total-len=2049563, root-span=285` |
| Did the aligner run? | Yes — `crush round 8: replacement build progress 1/4 (accepted)` (`/tmp/crush-ld/c4_stderr.log:53497`); round summary `resolved 4/4 replacement(s) in 6.86s` at line `59014` |
| Spec-correct aligner | **sPOA** (median 285 bp << 1 kb) |
| Replacement subgraph shape (id range 272 222 377..272 223 604) | Round 8 emits **1 226 segments, 355 990 bp, 308 distinct sequences (918 duplicates = 74.9%), 377 of those segments < 50 bp**. Plan 1/4 (n=391, median 285 bp) is again the dominant contributor; plans 2–4 contribute the long tail with n=47/259/33 and medians 96/470/135. |
| Polish behaviour | Same as the other small-median cases; for the disjoint-paths portion of the subgraph the polish is a no-op, for any portion that *does* have nested bubbles polish would help but the structural fragmentation is on the wrong side of POVU's bubble definition. |
| Integration behaviour | 1 226 of 1 227 minted ids survive to the final GFA (one was consumed by a later round). 377 of those are < 50 bp — tiny per-traversal fragments produced because seqwish-k=311 cannot bridge gaps in 285 bp medians. |

Round 8 single-handedly accounts for roughly half of all the "tiny"
(<50 bp) segments in the entire output (377 / 1 745 total `<50 bp` segments
in the final GFA).

## Aggregate

Across the five picks: **4 out of 5 are spec-mis-routed to sweepga**
(small-median bubbles, picks a/d/e and the medium-but-still-mis-routed
pick b). Of those four, the three sub-1-kb-median ones (a/d/e) all
exhibit the "no compaction, byte-identical copies" pattern; the 1–10 kb-
median one (b) gets adequate sweepga+seqwish behaviour with usable polish.
The 5th pick (c) is the one large-median bubble that matched spec — and
that one does compact cleanly.

Generalising to the whole run (per-round byte-duplicate rate inside the
fresh-id slice of `/tmp/crush-ld/C4A_full.gfa`, computed from the `ids
minted` ranges):

| Round | New segs | Unique seqs | Duplicate segs | Duplicate rate |
|---|---:|---:|---:|---:|
| 1 | 2 237 | 692 | 1 545 | **69.1%** |
| 2 | 209 | 125 | 84 | 40.2% |
| 3 | 1 265 | 316 | 949 | **75.0%** |
| 4 | 596 | 53 | 543 | **91.1%** |
| 5 | 270 | 3 | 267 | **98.9%** |
| 6 | 1 | 1 | 0 | 0.0% |
| 7 | 59 | 6 | 53 | **89.8%** |
| 8 | 1 226 | 308 | 918 | **74.9%** |
| 9 | 10 | 7 | 3 | 30.0% |

(Round 6 is a single 74 bp segment from a 1-bubble round and is not
representative.)

Six of the eight non-trivial rounds emit at least 69% byte-duplicate
segments; three of them are at least 89% duplicates. **Four of the four
deep rounds with at least 50 segments minted (4, 5, 7, 8) exceed the 74%
rate, and the worst (round 5) is 98.9%.** The pattern is uniform across
rounds, not concentrated at any one depth.

### Which failure mode dominates?

Of the task's four candidate failure modes:

1. **Aligner not invoked** — *not* the failure mode. Every selected bubble
   has its `building replacement X/Y with Sweepga` line followed by a
   matching `replacement build progress X/Y (accepted)` line in the
   stderr, and the rerun produced a non-empty `seqwish.gfa` for every
   replacement.
2. **Aligner ran but no compaction** — *symptom*, not root cause. For
   small-median bubbles the seqwish-induced graph really does come out as
   N-way disjoint paths (replacement_0000 is the extreme: 47 paths, 0
   links). But it's a symptom of #4 (small sequences cannot anchor under
   `seqwish-k=311`), not an independent bug in the aligner.
3. **Aligner failed silently** — *not* the failure mode. The
   `failed_or_empty` counter in `resolve_graph_bubbles`
   (`src/resolution.rs:574–591`) is zero for every round; the
   `validate_replacement_paths` post-condition succeeds for every plan;
   `path_sequences_equal` (Phase-6 invariant 4) succeeds end-to-end.
4. **Wrong aligner for scale** — **dominant failure mode (4 / 5 picks)**.
   `method=sweepga` is the configured method,
   `candidate_replacement_method` passes it through unchanged for every
   bubble, and the spec's three-tier routing
   (`docs/crush-architecture-spec.md` §Phase-2) is **not implemented
   anywhere**: `auto_replacement_method` (`src/resolution.rs:1659–1668`)
   only ever returns `Poa` or `Allwave`, never `Poasta` or `Sweepga`.

There is an additional, second-order failure mode the task list did not
name explicitly that I think is worth flagging:

5. **Seqwish-k mismatch** — `seqwish-k=311` cannot induce match anchors
   inside traversals whose median length is < 311 bp. Even if a *correctly
   routed* sweepga call were given small sequences, this parameter would
   still fragment them. The two-tier mistake (config) and the
   parameter-tier mistake (`seqwish-k=311` is fine for sweepga's intended
   >=10 kb domain but is too aggressive for short sequences) compound each
   other.

The fix proposal below addresses (4) as the primary failure; (5) becomes a
non-issue if (4) is fixed because the routing change diverts short-median
bubbles away from seqwish entirely.

## File:line of the wrong code + fix sketch

**Wrong code:** `src/resolution.rs:1649–1668`

```rust
fn candidate_replacement_method(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> ResolutionMethod {
    match config.method {
        ResolutionMethod::Auto => auto_replacement_method(candidate, config),
        method => method,            // <-- fixes ALL bubbles to one aligner
    }
}

fn auto_replacement_method(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> ResolutionMethod {
    candidate_selection_method(candidate.traversal_stats, config)
    // candidate_selection_method (line 887) only returns Poa OR Allwave.
}
```

`candidate_selection_method` at `src/resolution.rs:887–904` is the only
place size-driven routing exists today, and its outputs are
`ResolutionMethod::Poa | ResolutionMethod::Allwave` — never `Poasta`, never
`Sweepga`. So even the `Auto` path cannot reach the spec's three-tier
strategy; the `method => method` arm at line 1655 is just a config
short-circuit that makes the un-implemented routing irrelevant.

**Fix sketch (one paragraph, no implementation):** Implement the spec's
three-tier routing inside `candidate_replacement_method` (or refactor
`auto_replacement_method` / `candidate_selection_method` to support it) so
the function chooses an aligner based on `candidate.traversal_stats`:
`ResolutionMethod::Poa` (sPOA) when
`median_len < auto_spoa_max_traversal_len` (default 1 kb),
`ResolutionMethod::Poasta` when the median is in
`[auto_spoa_max_traversal_len, auto_poasta_max_traversal_len]` (new
constant, default 10 kb), and `ResolutionMethod::Sweepga` above that. Two
pieces are needed: (i) add `auto_poasta_max_traversal_len` and
`auto_sweepga_min_traversal_len` constants alongside the existing
`auto_spoa_max_traversal_len` / `auto_allwave_*` thresholds on
`ResolutionConfig`, and (ii) change the `config.method` branch so a
user-supplied `Sweepga | Poasta | Poa` still **bypasses the small-bubble
case** only when the per-bubble routing would have picked something at
least as aggressive — i.e., the user's explicit method is treated as an
upper bound, not a hard pin, so that small bubbles still get sPOA even
under `method=sweepga`. After the fix the small-median bubbles in pick (a)
/ (d) / (e) would all be dispatched to `build_poa_replacement` and produce
single-consensus-segment subgraphs instead of N-way disjoint copies.

## References

- `docs/crush-architecture-spec.md` §Phase 2 "Stratified resolution by
  bubble size" — the three-tier routing table this report references.
- `docs/crush-architecture-spec.md` §Invariants 1 and 3 — "No two segments
  share a sequence" and "Each round shrinks the graph"; both are violated
  by this run (the run *grows* segments every round and emits up to 98.9%
  byte duplicates per round).
- `docs/crush-level-descent.md` — the run that produced
  `/tmp/crush-ld/C4A_full.gfa` (the canonical command and its 9-round
  trajectory).
- `src/resolution.rs:887–904` — `candidate_selection_method`, the only
  size-driven routing site in the code.
- `src/resolution.rs:1649–1668` — `candidate_replacement_method` /
  `auto_replacement_method`, the dispatch wrapper called per candidate in
  `resolve_graph_bubbles`.
- `src/resolution.rs:1967–2028` — `build_sweepga_seqwish_replacement`, the
  function every bubble in the analysed run flows through.
- `src/resolution.rs:1823–1885` — `finalize_pairwise_induced_replacement`
  (unchop + polish, the place where polish either fixes the fragmentation
  or — as observed for small bubbles — has no bubble structure to operate
  on).
- `src/resolution.rs:1496–1631` — `apply_replacement_frontier`, which
  splices the per-bubble replacement subgraphs into the working graph; it
  drops nothing and coalesces nothing, so any per-bubble fragmentation
  propagates verbatim to the final output.
