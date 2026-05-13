# Design: Syng-seeded transitive queries via anchor-based projection

**PR:** pangenome/impg#162 (feature/syng-integration)
**Date:** 2026-04-16
**Status:** Implemented

## Context

Syng's `query_region` returns homologous intervals at syncmer resolution
(~63bp granularity with default k=8, w=55), padded by a user-specified
amount. These intervals are approximate — syncmer boundaries don't align
with biological breakpoints, and on sparse-anchor paths (divergent
homology, RC matches) a single biological region can fragment into
multiple padded intervals.

For the target use case — closely-related genomes (pangenomes of the
same species) — the homology structure is simple and mostly collinear.
The shared-syncmer anchors inside each homology already encode the
alignment structure at syncmer resolution; base-pair-precise boundaries
can be obtained by linear extrapolation from the innermost anchors to
the user's query edges. This matches what PAF-based queries give you
for the same case: precise coordinates without heavyweight realignment.

## Architecture

### Per-hop pipeline

```
syng query_region_with_anchors (returns per-homolog shared-syncmer anchors)
  → bedtools-style distance merge within each (target, strand) group
    → per merged homolog, linear-project query edges to target forward
      strand using innermost anchors
        → output (genome, start_refined, end_refined, strand)
          → (optionally) feed into next hop's syng seed
```

**No alignment. No sequence fetches for refinement.** For closely-related
genomes the anchor grid is dense enough that linear extrapolation gives
coordinate precision within the local indel burden (typically single-digit
bp). Divergent or RC cases with sparse anchors are less precise but still
consistent.

### Step-by-step

#### 1. Anchored syng query

`SyngIndex::query_region_with_anchors` returns `Vec<HomologousIntervalWithAnchors>`
— one entry per (target path, strand) that shares at least one syncmer
with the padded query region, with `anchors: Vec<Anchor>` giving the
per-hit `(query_pos, target_pos, node_id)` shared-syncmer positions.

Both forward- and reverse-orientation GBWT nodes are queried (see
`build_fast_locate` for orientation-preserving GBWT encoding). RC
homologs are emitted with `strand='-'`; forward with `strand='+'`.
Intervals on the same target path but different strands are kept as
separate homologs — they describe distinct biology.

#### 2. Distance merge

In `syng_transitive::distance_merge_anchored`, each `(target, strand)`
group is merged using a bedtools `-d` semantic: two padded intervals
whose target-axis gap is `≤ merge_distance` are coalesced and their
anchor sets unioned. Honours the CLI's `-d` / `--min-distance-between-ranges`
(default 0 = overlap-only, matching the existing PAF-based query
convention).

This collapses the syncmer-sparsity fragmentation: where padded hits
don't quite overlap but sit within `-d` bp of each other on the target,
they merge into a single biological homolog.

#### 3. Linear projection per merged homolog

For each merged homolog, the user's query edges `qs`, `qe` are projected
onto the target forward strand using the merged anchor list:

- **Left edge** (`qs`): project from `leftmost = min(anchors by query_pos)`.
- **Right edge** (`qe`): project from `rightmost = max(anchors by query_pos)`.

Projection formula:
- `'+'` strand: `t(q_pos) = t_anchor + (q_pos - q_anchor)`
- `'-'` strand: `t(q_pos) = t_anchor + syncmer_len - (q_pos - q_anchor)`

The `+ syncmer_len` on the RC case accounts for the anchor's target
position pointing at the *start* of its syncmer on the forward strand,
while under RC the query base `q_anchor` corresponds to the *last* base
of the target syncmer. Results are clamped to `[0, target_len]` and to
the padded syncmer-resolution bounds (refinement can only tighten, not
extend, the homology window).

Accuracy is bounded by the indel burden in `[qs, q_anchor_leftmost]` and
`[q_anchor_rightmost, qe]`: for closely-related genomes these are typically
small gaps containing at most a few indels per kb.

#### 4. Multihop iteration

For `max_depth > 1`: refined intervals from hop k seed hop k+1. Visited
`(genome, start, end, strand)` tuples are deduped to prevent cycles.
Because boundaries are projected linearly (not via compounding padding),
slop does NOT grow across hops.

## Key design decisions

### Why linear projection instead of realignment?

Earlier iterations used BiWFA on anchor-flanked windows per edge per
homolog. On real pangenome data (e.g. yeast235 with ~313 haplotypes
per region) this produced:

- Thousands of AGC sequence fetches per query — the dominant cost.
- Subtle coordinate errors from window selection edge cases.
- ~10-minute wall time for a 2kb query over the yeast pangenome.

Linear projection from innermost anchors:

- Zero sequence fetches during refinement.
- Correct to within the local indel burden for closely-related genomes.
- <1s wall time on the same query.

If future use cases need tighter boundaries across divergent or RC-heavy
homology, a validation-alignment step can be added downstream without
changing this core path.

### Why not validate the projection with a quick alignment?

For closely-related genomes the indel burden is typically <5% per kb, so
linear projection error is bounded by ~50bp per kb. That's comparable to
what PAF-based queries produce from their stored CIGARs (the CIGARs are
from alignments that themselves have some error). Adding a validation
alignment would double the wall time for marginal gain.

### Why keep RC as a separate homolog, not try to merge strands?

A target genome can have both forward and RC-matching regions relative
to the query — these are different biology (a direct duplicate vs. an
inversion). Reporting them as separate intervals with explicit `strand`
gives downstream tools (GFA construction, untangle) the information to
handle them correctly.

## Integration with existing code

- `src/syng.rs`: `query_region_with_anchors` exposes the anchor graph
  from syng; `Anchor` and `HomologousIntervalWithAnchors` types.
- `src/syng_transitive.rs`: `distance_merge_anchored`, `project_query_to_target`,
  `refine_boundaries`, `one_hop`, `query_transitive`.
- `src/main.rs`: syng-index query routes through `query_transitive` for all
  output formats (bed, fasta, gbwt, gfa) when `--syng-raw` isn't set.
  Passes `query.transitive_opts.min_distance_between_ranges` as the
  merge distance.

## Known limitations

1. **Indel burden in outermost gaps.** Linear projection accuracy is
   bounded by the indel burden in the unaligned region between the
   innermost anchor and the query edge. For closely-related genomes
   (<5% divergence) this is single-digit bp; for more divergent cases
   precision degrades. A future validation alignment pass could tighten
   this but is not necessary for the stated use case.

2. **RC syncmer sparsity.** Closed syncmers are not RC-symmetric — a
   window on forward strand picks a minimizer k-mer that is generally
   NOT the minimizer of the same window on reverse strand. So RC
   homologs have ~10× fewer shared syncmers per kb than forward
   homologs. Combined with `-d`, users querying RC-heavy regions may
   need a larger `-d` value to collapse the fragmentation. Out of scope
   to fix in the syncmer scheme.

3. **Syng walk_path first-syncmer offset.** `GbwtPathStart` now carries
   `first_syncmer_pos` so `walk_path` returns absolute sequence
   coordinates (previously returned relative-to-first-syncmer, silently
   shifting every query). Indexes built before this fix load
   `first_syncmer_pos = 0` and retain the old behaviour — rebuild the
   index for correct anchor positions.
