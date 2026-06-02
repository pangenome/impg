# Syng GFA scaffold-chain filtering

This note documents the conversion-time filtering used by syng-to-GFA range
materialization. It is intentionally a conversion-layer rule: the `.syng` index
and the syng query collector still expose raw syncmer hits, while the local GFA
writer decides which selected-range occurrences may become shared graph
topology.

## Goal

The raw syng walk can contain exact syncmer hits that are true colinear anchors
in one selected path pair and unrelated off-chain coincidences elsewhere. A
node-level filter cannot distinguish those cases: once one occurrence of a
syncmer node is supported, every other occurrence of that node is also left
shared. The GFA writer now filters at occurrence level.

An occurrence is the path-local index of a retained syncmer after node-level
frequency masking has removed high-copy seeds. The identity is:

```text
(selected_path_index, post-frequency-mask_syncmer_occurrence_index)
```

This occurrence identity is local to one conversion call. It is not stored in
the syng index and is not exposed as a persistent graph identifier.

## Algorithm

For `gfa:syng:mask` and equivalent local syng GFA paths, conversion proceeds in
this order:

1. Collect the selected ranges as oriented local syncmer walks. Forward ranges
   use range-relative forward coordinates. Reverse ranges use reverse-oriented
   range coordinates, so colinear anchors are still increasing in the emitted
   path coordinate system.
2. Select explicit high-frequency syncmer nodes using `top=` and/or
   `max-occ=`. In the default occurrence-level policy, selected nodes are not
   removed globally. Their occurrences are private-split only when they lack
   support from `freq-run` consecutive syncmers or `freq-span` bp of repeated
   exact local sequence. `legacy-freq-mask=true` restores the older node-level
   removal behavior for comparisons.
3. Build occurrence candidates from the post-frequency-selection walks. Only
   syncmer nodes with more than one remaining occurrence are candidates for
   scaffold filtering; singletons can be emitted as ordinary private topology
   because no other path occurrence shares their node. The converter also
   records a compact local node spectrum over the remaining shared syncmers:
   total local occurrences, carrying path/haplotype count,
   occurrence-per-carrying-path ratio, maximum copies on any one path, and
   maximum path-local positional span for repeated copies. Nodes in the
   dispersed high-copy tail are treated as a second high-frequency source when
   they satisfy all calibrated defaults: at least `64` local occurrences,
   occurrence-per-carrying-path ratio at least `2.0`, at least two copies on one
   carrying path, and maximum path-local span at least the scaffold gap budget
   (`1 kb` by default). These nodes are not blindly private-split; their
   occurrences use the same `freq-run` / `freq-span` rescue as explicit
   high-frequency selections.
4. Build bounded scaffold candidates from repeated windows of `min-run`
   consecutive shared syncmer occurrences. The old exhaustive all-pairs
   expansion was:

   ```text
   all occurrences of node X on path A x all occurrences of node X on path B
   ```

   That is too large for C4-scale high-copy syncmers. The current conversion
   first indexes path-local shared-node windows whose adjacent anchors are
   within the normal scaffold gap. A pairwise anchor is only materialized when
   the full `min-run` node window is observed on at least two selected paths.
   Each repeated window directly supports its member occurrences. This
   preserves the configured consecutive-run semantics even when syncmer anchors
   overlap in base coordinates. Pathologically dense repeated windows, or
   otherwise moderate windows that would exceed the global candidate-anchor
   budget in aggregate, stop there rather than allocating every path-pair
   combination.
5. For repeated windows that remain within the candidate budget, also convert
   matching occurrences into exact PAF-style syng anchors:

   ```text
   query_pos  = occurrence coordinate on the lower-index path
   target_pos = occurrence coordinate on the higher-index path
   node_id    = absolute syncmer node id
   length     = syncmer length
   identity   = 1.0
   ```

6. Pass those bounded anchors to
   `syng_transitive::chain_anchors_with_sweepga_scaffold_mass`. This is the same
   SweepGA-backed scaffold-chain path used by syng query collection. The GFA
   conversion does not reimplement SweepGA's diagonal-offset, drift, or
   plane-sweep logic; it only prunes anchors that cannot satisfy the configured
   `min-run` scaffold context.
7. Use `mask,min-run=` as the scaffold mass threshold: the minimum scaffold
   length is `min_run * syncmer_len`. The scaffold gap follows the query
   collection default `DEFAULT_EXTEND_BUDGET_BP`. The default therefore remains
   conservative, and existing `:mask` / `:filter` syntax continues to tune the
   support requirement.
8. Mark every occurrence whose pairwise anchor is returned as a member of a
   retained SweepGA scaffold chain. These occurrences remain shared syncmer
   topology.
9. Optionally apply the existing `sequence-k=` exact-span rule as an additional
   occurrence-level support source. This keeps previous exact sequence-span
   behavior but no longer blesses all occurrences of a node.
10. Any remaining shared-node occurrence that is not supported by its applicable
   policy is split into a private occurrence during GFA materialization.

## Node vs occurrence semantics

The conversion has three distinct decisions:

| Decision | Granularity | Effect |
| --- | --- | --- |
| Frequency selection (`top=`, `max-occ=`) | Node for selection, occurrence for output | Selects explicit high-frequency nodes; unsupported occurrences are split/private, while occurrences rescued by `freq-run` or `freq-span` remain shared. |
| Spectrum scaffold-glue source | Node for selection, occurrence for output | Selects dispersed high-copy local-spectrum tail nodes and applies the same `freq-run` / `freq-span` rescue before private splitting. |
| Scaffold-chain support (`min-run=`) | Occurrence | Keeps non-spectrum shared syncmer occurrences that are members of retained SweepGA scaffold chains. |
| Local repeat context rescue/split | Occurrence | Clones rare repeated local contexts when a near-single-copy node appears in a minor context. |

The important change is the middle row. A syncmer node can now have both shared
and private occurrences in the same output GFA:

```text
path A: ... x y Z ... Z ...
path B: ... x y Z ...
path C: ........ Z ...
```

If `x,y,Z` form a retained SweepGA scaffold chain between `A` and `B`, those two
`Z` occurrences remain shared. The extra `A` occurrence and the singleton `C`
occurrence are not chain members, so they are emitted as private syncmer
segments instead of reusing node `Z`.

## Relation to syng query collection

Syng query collection already turns raw syncmer hits into exact PAF-style
records and delegates chain/scaffold filtering to SweepGA through
`chain_anchors_with_sweepga_scaffold_mass`. Conversion now reuses that helper
for bounded local selected-path candidates. This keeps off-diagonal and
indel-shifted colinear handling in one implementation once an occurrence pair
has enough local context to be a plausible scaffold member:

- indel-shifted but colinear runs can remain one scaffold chain;
- short off-chain singleton hits fail the scaffold mass threshold;
- no custom diagonal test in the GFA writer decides chain membership;
- high-copy singleton coincidences are pruned before SweepGA instead of being
  expanded into an all-pairs anchor set.

No SweepGA filtering semantics were changed for this integration.

## Path reconstruction

Filtering happens before shared topology is materialized, but path spelling is
preserved:

- in the default occurrence policy, explicit high-frequency and
  spectrum-selected scaffold-glue syncmers are not deleted, only cloned/private
  when unsupported, so they keep the same syncmer sequence on the path;
- in legacy node-level frequency masking, removed anchors are bridged with the
  selected input sequence when sequence files are available;
- exact blunt output still trims syncmer slices by the same overlap/clip logic
  used before this change;
- optional N cutting is the only intended path-spelling change, and only for
  gap DNA according to the configured `cut-n` policy.
