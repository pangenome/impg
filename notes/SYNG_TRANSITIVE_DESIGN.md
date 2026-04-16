# Design: Syng-seeded transitive queries via boundary realignment

**PR:** pangenome/impg#162 (feature/syng-integration)
**Date:** 2026-04-16
**Status:** Design — not yet implemented

## Context

Syng's `query_region` returns homologous intervals at syncmer resolution
(~63bp granularity with default k=8, w=55), padded by a user-specified
amount. These intervals are approximate — syncmer boundaries don't align
with biological breakpoints, and the slop compounds at every transitive
hop. Without base-pair-precise coordinate resolution, multi-hop queries
explode to chromosome-scale intervals.

The naive solution is to do full-region realignment at each hop
(FastGA/WFMASH over all pairs of fetched homolog sequences, build an
in-memory Impg, project through the implicit graph). That is expensive:
`O(region²)` per pair, dominated by subprocess and index setup overhead
when regions are small. It also computes much more than we need —
interior CIGARs are not used for transitive projection, only endpoints.

This design takes the leaner path: only resolve the two fuzzy boundaries
of each homolog, using small local BiWFA alignments anchored to shared
syncmer positions. Directly analogous to how impg's TPA/tracepoint
approximate mode does targeted realignment around specific query
coordinates rather than reconstructing full CIGARs.

## Architecture

### Per-hop pipeline

Each transitive hop follows this sequence:

```
syng one-hop query (returns homologs + shared-syncmer anchors)
    → for each homolog:
        → find nearest shared anchor to each query edge
            → fetch small window on query & target around each edge
                → BiWFA align window
                    → trace CIGAR to project edge → precise coord
                        → output (genome, start_refined, end_refined)
```

**There is no in-memory Impg.** No PAF, no graph induction. Syng is the
homology oracle; BiWFA is the coordinate snapper. Nothing else.

### Step-by-step

#### 1. Syng one-hop expansion with anchors

Extend `SyngIndex::query_region` (or add a sibling primitive) to return
both the homolog interval and the set of shared-syncmer anchors that
produced it:

```rust
pub struct HomologousIntervalWithAnchors {
    pub genome: String,
    pub start: u64,           // syncmer-resolution + padding
    pub end: u64,             // syncmer-resolution + padding
    pub strand: char,
    /// (query_pos, target_pos, node_id) for each shared syncmer,
    /// sorted by query_pos. At least one anchor per homolog.
    pub anchors: Vec<(u64, u64, i64)>,
}
```

Anchors are already computed during the path walk inside `query_region`
— they are discarded when `HomologousInterval` is built. We just need to
preserve them.

#### 2. Pick anchors flanking each query edge

For the left edge (`qStart`):
- Find the rightmost anchor with `query_pos <= qStart` — call it
  `L = (q_L, t_L)`.
- Find the leftmost anchor with `query_pos > qStart` — call it
  `R = (q_R, t_R)`.

Fallback cases:
- **No anchor with `query_pos <= qStart`**: the query edge sits before
  the first shared syncmer. Use `L = (0, 0)` (the path origin acts as
  an implicit anchor) or fall back to the syncmer-resolution `start`
  unchanged.
- **No anchor with `query_pos > qStart`**: rare; query edge is beyond
  the last anchor for this homolog. Same fallback.

Symmetric for the right edge.

#### 3. Fetch small windows

For the left edge, fetch:
- `query_slice = query_seq[q_L .. q_R]` — typically one inter-syncmer
  gap, so ≤ ~2 × syncmer_len bp (~120bp at defaults)
- `target_slice = target_seq[t_L .. t_R]`

Both slices come from the `UnifiedSequenceIndex`. Fetches are tiny
(a few hundred bp), so AGC or FASTA random-access cost is negligible.

For reverse-strand homologs: `target_slice` is reverse-complemented
before alignment, and anchor positions on the target side are already
expressed in forward-strand coordinates by syng (confirm during
implementation — see open questions).

#### 4. BiWFA align

Call into `lib_wfa2` using the existing `with_aligner` thread-local
pattern at `src/impg.rs:53`, but with `MemoryMode::Ultralow` (BiWFA)
instead of the current default `High`. BiWFA gives O(s) memory and O(ns)
time — perfect for windows in the 100-500bp range.

```rust
with_biwfa_aligner(&Distance::Edit, |aligner| {
    aligner.align(&query_slice, &target_slice);
    parse_cigar_to_delta(aligner.cigar())
})
```

Output: CIGAR ops covering the `[q_L, q_R)` × `[t_L, t_R)` window.

#### 5. Trace CIGAR to project edge

Walk the CIGAR consuming `qStart - q_L` query bases. The corresponding
number of target bases gives the target position: `t_refined = t_L +
consumed_target_bases`. This is the precise base-pair projection of
`qStart` onto the target.

Symmetric for `qEnd` on the right edge.

Output: `(genome, t_refined_start, t_refined_end, strand)`.

#### 6. Multihop iteration

For `max_depth > 1`: collect refined intervals from hop k, feed each
back into step 1 as a new query at hop k+1. Deduplicate against a
visited set keyed on `(genome, start, end)`. Terminate when
`depth == max_depth` or no new intervals are discovered.

Because each hop snaps boundaries precisely, slop does **not** compound.
Hop-k intervals are just as tight as hop-0 intervals.

## Key design decisions

### Why boundary-only instead of full-region realignment?

- **Cost.** Full-region realignment is O(region²) per pair. Boundary
  realignment is O(slop²) per edge per pair — typically 100× cheaper on
  kb-scale regions, more on longer ones.
- **Right-sized.** We only need endpoint coordinates to project
  intervals. Interior CIGAR is not consumed anywhere in the transitive
  pipeline.
- **No subprocess, no temp FASTAs.** BiWFA is a direct library call. On
  tiny regions the FastGA/WFMASH subprocess setup would dominate
  wall-clock.
- **Matches existing impg machinery.** The TPA approximate-mode path
  already does exactly this — syncmers are the tracepoints, BiWFA is
  the re-aligner. We reuse `with_aligner` and the Ultralow memory mode.

### What we give up

The ability to project a sub-range of a current hop's interval through
the current hop's "edge" (the alignment we just computed). We never
store interior CIGARs. For transitive iteration this is fine — each hop
re-queries syng, which has its own full anchor graph. The "only need
endpoints" property holds at every hop.

If a caller later needs interior CIGARs (e.g. for GFA induction from
syng-seeded regions), that is a separate pipeline that *can* opt into
full-region realignment downstream. It's orthogonal to this design.

### Why BiWFA specifically

`lib_wfa2` exposes `MemoryMode::Ultralow`, which selects wfa2-lib's
BiWFA path — O(s) memory, O(ns) time. For 100-500bp windows with ≤10%
divergence, s is ≤50 and n is ≤500, so alignment finishes in
microseconds with a few KB of memory. No risk of the `O(s²)` memory
blowups that full WFA can hit on divergent pairs.

`lib_wfa2` is already a direct dep of impg, already linked in, and the
thread-local aligner pool pattern at `src/impg.rs:43-69` is already
used in production for TPA tracepoint reconstruction. We add one more
thread-local for `Ultralow` and one helper to mirror `with_aligner`.

### Anchor-selection edge cases

1. **No anchor before `qStart`.** Use path origin `(0, 0)` — syng
   guarantees both paths share their origin as an implicit anchor
   (to be confirmed — see open q). If that's not safe, fall back to
   returning the syncmer-resolution `start` unchanged for that edge.

2. **No anchor after `qEnd`.** Symmetric — use path end (or fall back).

3. **`qStart == anchor_pos` exactly.** `q_L == qStart`, so
   `consumed = 0` and `t_refined == t_L`. Correct; no alignment needed.

4. **Gap between `q_L` and `q_R` is too large.** E.g., query falls in a
   long syncmer-desert. If `q_R - q_L` exceeds a threshold (say 2kb),
   decline to realign and return the syncmer-resolution edge — the
   cost/benefit inverts for large windows. Configurable cap.

5. **Anchor on wrong strand.** Syng currently always reports
   `strand: '+'`; RC homologs are handled internally via the node sign.
   For boundary realignment we need the target-side anchor positions to
   be expressed in a consistent frame with the fetched (possibly
   RC'd) target slice. See open q 3.

### Integration with existing code

- `SyngImpgWrapper::query_transitive_seeded(target_id, range_start,
  range_end, max_depth, padding, sequence_index) -> Vec<HomologousInterval>`
  implements the full per-hop pipeline + multihop loop.
- `SyngImpgWrapper::query` (the trait method on `ImpgIndex`) routes to
  `query_transitive_seeded` with `max_depth = 1` by default.
- `--syng-raw` CLI flag preserves the current debug behavior (raw
  syncmer-resolution intervals, no realignment).

## Files modified

- `src/syng.rs` — expose anchors from `query_region`, either by
  augmenting `HomologousInterval` or by adding a new
  `query_region_with_anchors` primitive.
- `src/lib.rs` — `SyngImpgWrapper::query_transitive_seeded`, update
  trait method implementations to route through the new pipeline by
  default.
- `src/impg.rs` — add `Ultralow`-memory BiWFA aligner to the
  thread-local pool + `with_biwfa_aligner` helper (or extend
  `with_aligner` to take a memory mode).
- `src/main.rs` — add `--syng-raw` flag to the query subcommand; route
  default `--syng` to the new pipeline.
- `tests/test_syng_integration.rs` — parity test comparing
  syng-realign transitive output against PAF-based
  `query_transitive_dfs` on a shared fixture.

## Status of open questions

### Resolved

1. **Anchor surface from syng.** Added `query_region_with_anchors`
   returning `Vec<HomologousIntervalWithAnchors>`; `query_region` is
   now a thin wrapper that drops anchors. See `src/syng.rs`.

2. **Strand / RC homology.** Implemented: `query_region_with_anchors`
   queries both GBWT node orientations (`2*N` and `2*N+1`), tags
   anchors per (query_orient XOR target_orient), and emits forward
   and RC hits as separate intervals with `strand='+'` / `strand='-'`.
   Also fixed `build_fast_locate` to preserve orientation in the GBWT
   encoding (previously collapsed via `unsigned_abs()`, which silently
   dropped every RC visit). Boundary realignment handles `strand='-'`
   by RC'ing the target slice before BiWFA and projecting offsets in
   the correct direction.

### Known limitations (carried forward)

3. **Sparse RC-shared syncmers.** Syng's closed-syncmer selection is
   not RC-symmetric: a window on forward strand picks a minimizer
   k-mer that is generally NOT the minimizer of the same window on
   the reverse strand. So when genome A and genome B share a region
   under RC homology, only a small fraction of syncmers (empirically
   ~3–5 per kb in random sequence, vs ~35 per kb on matched forward
   strand) are shared between them. Consequences:
   - Padded-syncmer-resolution RC homologs often fragment into
     multiple intervals with inter-anchor gaps that exceed the
     padding.
   - Boundary realignment has fewer anchor options per edge; some
     edges fall back to the syncmer-resolution bound.
   - `test_syng_rc_homolog_end_to_end` validates via interval
     overlap (≥200bp with the expected window) and base-content
     alignment (≥30bp exact-match run when RC'd), rather than exact
     coordinates — the coordinate precision achievable depends on
     the specific sequence.

   Not a bug in impg's pipeline — inherent to the syng scheme. A
   future improvement could use RC-symmetric syncmer selection (e.g.,
   open syncmers with a rotational-symmetric hash). Out of scope.

4. **Syncmer-desert cap.** What threshold for `q_R - q_L` before we
   stop realigning? Default proposal: `MAX_REALIGN_WINDOW_BP = 2048`.
   Could become configurable if needed.

5. **Multihop frontier deduplication.** Keying on
   `(genome, start, end)` is correct for exact-interval dedup, but
   near-duplicates (off by 1bp after realignment noise) could cause
   redundant work. Consider interval-tree dedup with a small tolerance
   instead of a hashmap.

6. **Path-origin as implicit anchor.** Not relied on — the actual
   fallback is "return the syncmer-resolution padded edge" when a
   flank is missing, which is conservative and safe regardless of
   index-build assumptions.
