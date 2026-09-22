# Crush — Architectural Specification

This is the canonical specification for what the `crush` operation in impg should do. It exists because three days of agents working on crush have been flailing without a clear spec, each guessing at the algorithm and producing inconsistent or buggy results. **No code change to crush should be proposed without referencing this document.**

## Goal

Take a syng graph (bubbly topology, ~minutes-to-produce) and emit a **reasonably compacted variation graph** suitable for downstream analysis and genotyping.

**Fast** — minutes end-to-end, achieved by **avoiding all-to-all alignment** in favor of local per-bubble condensation. This is the whole point of crush vs the PGGB pipeline.

PGGB is the quality oracle for comparison, but byte-for-byte PGGB equivalence is **not** the goal.

## Input

A syng graph: bubbly topology where most edges and segments are correct, with **possible spurious runs of matching synchromers**, sometimes 5–10 in a row. These are mostly fine but we may want to prune more aggressively than the current 3-at-a-time policy.

## Algorithm

### Phase 1 — Bubble detection
POVU / flubble detection (Lean-proven correct). For each iteration: identify **top-level non-overlapping flubbles**.

### Phase 2 — Stratified resolution by bubble size

Route each non-overlapping bubble by **median traversal length**:

| Median length | Aligner | Why |
|---|---|---|
| < ~1 kb | **sPOA** | Linearizes consensus; short tandem repeats and small recurrent motifs become clean linear representations |
| ~1 kb – ~10 kb | **POASTA** (possibly via allwave first) | Fast POA for medium scale; still produces near-linear output |
| ≥ ~10 kb | **sweepga** | Above 10kb we expect non-linear patterns (gene conversion, non-allelic homologous recombination). Can't linearize — sweepga handles non-linear topology |

The aligner is chosen by **what the output is allowed to look like** at each scale: small bubbles must linearize; large bubbles preserve non-linear topology.

### Phase 3 — Local bubble graph construction

For each bubble, independently:

1. Extract the paths **inside** the bubble (between boundary nodes)
2. Pass those sequences to the appropriate aligner
3. Aligner returns a local compacted graph
4. **Allocate fresh node IDs starting at `n + current_graph_size`** so the new local graph's IDs are outside the existing range and cannot collide with surviving segments

**Boundary handling (open question to nail in implementation):** what counts as "inside" the bubble. The boundary nodes themselves are shared with the rest of the graph — they MUST NOT be re-emitted with fresh IDs. We may need to realign content right up to (but not including) the boundary segments. To be confirmed against POVU's site definitions.

### Phase 4 — Strike + link

For each resolved bubble:
1. **Strike (remove)** the old inner segments from the parent graph
2. **Insert** the new local-aligner segments
3. **Add edges** connecting boundary nodes to the first/last segments of each path through the new local graph

Old inner nodes are GONE. New inner nodes have IDs in the n+size range. Boundaries are preserved with original IDs.

### Phase 5 — Batched lacing

Process the **whole non-overlapping batch** in one round:
1. Compute all resolutions in parallel (rayon over the batch)
2. Apply strike+link as a single batch update on the parent graph
3. Sort the resulting graph

Non-overlap is what makes batching safe. Non-overlap is **guaranteed by the level hierarchy** if we do top-level first.

### Phase 6 — Level descent (iterate)

After resolving level-1 bubbles:
- Inside each resolved bubble there may be smaller bubbles → resolve at level-2
- **Do not re-resolve the same level-1 bubble** — that's wasted work
- Continue descending until no resolvable bubbles remain

### Phase 7 — POA normalization of small recurrent motifs

After sweepga (or allwave) on a bubble, the local graph may contain tight high-copy nodes from short tandem repeats. These should be **re-aligned locally with sPOA / POASTA** to linearize. This is the "linearize small recurrent motifs" step and is what produces clean compact output for downstream genotyping.

### Phase 8 — (Optional) final compaction

Run a gfaffix-style sequence-dedup pass on the final emit. If Phases 1–7 are correct, this should be a no-op — duplicate-sequence segments shouldn't exist by construction.

## Parameters

- **seqwish-k**: `311` has been working well. Smaller values are possible but require more POA/POASTA normalization downstream. Trade-off between aligner sensitivity and downstream cleanup.
- **Synchromer pruning at input**: current `min-run=3`; consider going further (5–10) for spurious-match suppression.

## Performance goal

End-to-end **minutes** on a syng graph for a single locus like C4 GRCh38. Achievable because:
- No all-to-all alignment (that was PGGB's bottleneck)
- Per-bubble work is local and rayon-parallelizable
- Level hierarchy bounds iteration depth

## Invariants that fall out of the algorithm

A correct implementation produces output that satisfies, **by construction**:

1. **No two segments share a sequence.** If two segments have identical bytes, that's evidence of a bug — bubbles overlap, boundary sequences were absorbed into a replacement, or the local aligner emitted redundant nodes.
2. **Every input path is preserved exactly** (per-path byte equality).
3. **Each round shrinks the graph** (segment count and segment bp both ≤ pre-round). If a round grows the graph, the algorithm is fighting itself — either bubble selection is wrong or strike-and-replace is leaving orphans.
4. **No segment ID is re-used** between original survivors and new replacement nodes (new IDs in `[n+size, ...)`).

## Non-goals

- PGGB byte-for-byte equivalence
- All-to-all sweepga (defeats the purpose)
- Matching the historical 4:39 baseline exactly (that baseline is itself underaligned per `docs/c4-crush-handoff.md`)
