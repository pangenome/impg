# Design: Syng-seeded transitive queries via implicit local graphs

**PR:** pangenome/impg#162 (feature/syng-integration)
**Date:** 2026-04-16
**Status:** Design — not yet implemented

## Context

Syng's `query_region` returns homologous intervals at syncmer resolution
(~63bp granularity with default k=8, w=55), padded by a user-specified
amount. These intervals are approximate — syncmer boundaries don't align
with biological breakpoints, and the slop compounds at every transitive
hop. Without cigar-based coordinate projection, multi-hop queries
explode to chromosome-scale intervals.

The existing impg codebase already solves precise coordinate projection
through alignment graphs via `ImpgIndex::query_transitive_{dfs,bfs}`,
which walks PAF cigars to project intervals through an implicit graph.
The missing piece is using syng as the *seed phase* that replaces
precomputed all-vs-all alignments with fast approximate one-hop lookups,
followed by targeted local realignment to produce the cigars needed for
precise projection.

## Architecture

### Per-hop pipeline

Each transitive hop follows this sequence:

```
syng one-hop (fast, approximate)
    → fetch sequences from AGC
        → pairwise align (FastGA/sweepga) → PAFs with cigars
            → load PAFs into in-memory ImpgIndex (implicit graph)
                → re-query through implicit graph (precise boundaries)
                    → output OR feed into next hop
```

**There is no explicit graph.** No seqwish, no GFA, no graph induction.
The implicit graph IS the set of PAFs — impg's core data structure.

### Step-by-step

#### 1. Syng one-hop expansion

Call `SyngIndex::query_region(genome, start, end, padding)` via the
FastLocate fast path (O(k log r) per syncmer node). Returns a set of
`HomologousInterval { genome, start, end, strand }` — rough coordinates
at syncmer resolution with padding slop.

These intervals are the *seed* — they tell us WHERE to look, not the
precise boundaries.

#### 2. Sequence fetch

For each HomologousInterval, fetch the actual sequence bytes from AGC
(or FASTA via the UnifiedSequenceIndex). This produces a set of
`(name, subsequence)` pairs ready for alignment.

Include the query region itself in the fetch set — it participates in
the all-vs-all alignment.

#### 3. Local pairwise alignment

Run FastGA or sweepga on the fetched sequences to produce PAF records
with full CIGAR strings. This is targeted alignment — only the regions
syng identified as potentially homologous, not a whole-genome all-vs-all.

The alignment scope is small (one query region + its one-hop homologs),
so this should be fast even on 2-vCPU CI runners.

**Key:** the aligner must produce PAFs with cigars (not just approximate
mappings). impg's coordinate projection requires base-level alignment
information.

#### 4. Build in-memory ImpgIndex

Parse the PAF records into the existing `Impg` / `ImpgIndex` structures.
This is the same code path used when loading PAFs from disk, except the
PAFs come from the local alignment in step 3 instead of a file.

The result is an implicit local graph — a set of pairwise alignments
that impg can query through using cigar-based interval projection.

#### 5. Re-query through the implicit graph

Call `impg_index.query()` or `query_transitive_dfs()` on the in-memory
index with the ORIGINAL query coordinates. The cigar-based projection
gives precise start/end on each target genome — base-pair resolution,
not syncmer resolution.

The output of this step is the final result for a single-hop query, OR
the input seed for the next transitive hop.

#### 6. Transitive iteration (optional)

For multi-hop queries (depth > 1), take the precise intervals from
step 5 and feed them back into step 1 as new query regions. Each hop:

```
precise intervals from previous hop
    → syng one-hop on each interval
        → fetch + align + re-query
            → new precise intervals
```

Because step 5 tightens the boundaries via cigars, the intervals fed
into the next hop are PRECISE, not sloppy. This prevents the
compounding-slop explosion that would occur if we chained syng queries
directly without realignment.

**Termination:** same as existing `query_transitive` — stop at
`max_depth` or when no new intervals are discovered.

## Key design decisions

### Why not build an explicit graph?

Running seqwish + smoothxg + gfaffix per hop is:
- Slow (seqwish transclosure is O(n²) on dense input)
- Unnecessary — impg's core query engine works on PAFs, not GFA
- Lossy — graph induction discards alignment detail that cigars preserve

The implicit graph (PAFs with cigars) is the natural input to impg's
existing coordinate projection machinery.

### Why realign at each hop instead of just expanding syng intervals?

Without cigars, each hop adds ~syncmer_len slop to interval boundaries.
After k hops: `slop ≈ k × (syncmer_len + 2 × padding)`. For k=3 with
63bp syncmers and 120bp padding: slop ≈ 900bp per boundary per
direction. That's manageable for small k but grows linearly.

More critically: without alignment, there's no way to know WHERE within
a syncmer interval the actual homology boundary falls. Two adjacent
syncmers that share a node might have a breakpoint between them that
only alignment can resolve.

Realigning at each hop costs O(local_region_size²) for the aligner but
gives base-pair precision, keeping subsequent hops tight.

### Padding strategy

The syng `padding` parameter controls how much slop to add around
syncmer hits in step 1. Too little → misses homologs whose breakpoints
fall between syncmers. Too much → fetches excess sequence, slowing
alignment.

Recommended default: `padding = syncmer_len` (63bp). This ensures that
any homology extending one syncmer beyond the query boundaries is
captured. The realignment in step 3 will tighten the result regardless.

### Integration with existing code

The existing `SyngImpgWrapper` in `src/lib.rs` already implements the
`ImpgIndex` trait, bridging syng queries into impg's query interface.
The transitive pipeline can be implemented as a new method on
`SyngImpgWrapper` (or a new struct) that:

1. Calls `self.syng_index.query_region()` for the one-hop seed
2. Fetches sequences via the existing `UnifiedSequenceIndex`
3. Calls the aligner (FastGA/sweepga) — may need a new internal API
   for in-process alignment rather than subprocess invocation
4. Builds a temporary `Impg` from the resulting PAFs
5. Calls `impg.query()` on the temporary index for precise projection

Steps 3-5 are already available as library functions in impg (used by
the partition and GFA pipelines). The new part is wiring them into a
loop seeded by syng.

## Files likely modified

- `src/lib.rs` — new `SyngTransitiveQuery` struct or method on
  `SyngImpgWrapper` implementing the hop pipeline
- `src/main.rs` — wire the transitive pipeline into the `query --syng`
  CLI path, replacing the current direct `query_region` call
- `src/syng.rs` — possibly expose additional query primitives
  (e.g., raw node-level hits before interval merging)
- Tests — parity test comparing syng-seeded transitive query against
  alignment-based transitive query on shared test data

## Open questions

1. **In-process vs subprocess alignment.** Currently FastGA/sweepga are
   invoked as subprocesses. For per-hop local alignment on small regions,
   subprocess overhead dominates. An in-process alignment API (calling
   the aligner as a library function) would be much faster. Does
   sweepga/FastGA expose a Rust library API, or do we need to add one?

2. **Alignment filtering.** The current pipeline applies overlap/identity
   filters to PAFs before loading into impg. Should the per-hop local
   alignment use the same filters, or relax them (since the regions are
   already pre-selected by syng as likely homologous)?

3. **Memory management.** Each hop creates a temporary in-memory
   ImpgIndex. For deep transitive queries on large regions, the
   accumulated PAFs could grow. Should we discard intermediate indexes
   after extracting the precise intervals, or keep them for debugging?

4. **Parallelism.** Per-hop alignment is embarrassingly parallel across
   independent query intervals. The existing `par_iter` in partition
   could be reused. But the transitive iteration itself is sequential
   (hop k+1 depends on hop k's output). Is there a useful fan-out
   within a single hop?

5. **Strand handling.** Syng currently only reports `strand: '+'`.
   Reverse-complement homology is captured via syng's internal RC path
   insertion, but the strand information is lost when collapsing via
   `unsigned_abs()`. Should the transitive pipeline recover strand from
   the alignment cigars instead?
