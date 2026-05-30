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
2. Apply the existing high-frequency syncmer mask. This remains node-level:
   a node rejected by `top=` or `max-occ=` is removed from all local walks, and
   adjacent retained anchors are bridged with source sequence or `N` fill as
   before.
3. Build occurrence candidates from the post-frequency-mask walks. Only syncmer
   nodes with more than one remaining occurrence are candidates for scaffold
   filtering; singletons can be emitted as ordinary private topology because no
   other path occurrence shares their node.
4. For every pair of selected paths that share candidate nodes, convert matching
   occurrences into exact PAF-style syng anchors:

   ```text
   query_pos  = occurrence coordinate on the lower-index path
   target_pos = occurrence coordinate on the higher-index path
   node_id    = absolute syncmer node id
   length     = syncmer length
   identity   = 1.0
   ```

5. Pass those anchors to
   `syng_transitive::chain_anchors_with_sweepga_scaffold_mass`. This is the same
   SweepGA-backed scaffold-chain path used by syng query collection. The GFA
   conversion does not implement its own diagonal-offset, drift, or plane-sweep
   logic.
6. Use `mask,min-run=` as the scaffold mass threshold: the minimum scaffold
   length is `min_run * syncmer_len`. The scaffold gap follows the query
   collection default `DEFAULT_EXTEND_BUDGET_BP`. The default therefore remains
   conservative, and existing `:mask` / `:filter` syntax continues to tune the
   support requirement.
7. Mark every occurrence whose pairwise anchor is returned as a member of a
   retained SweepGA scaffold chain. These occurrences remain shared syncmer
   topology.
8. Optionally apply the existing `sequence-k=` exact-span rule as an additional
   occurrence-level support source. This keeps previous exact sequence-span
   behavior but no longer blesses all occurrences of a node.
9. Any remaining shared-node occurrence that is not SweepGA-scaffold-supported
   or `sequence-k`-supported is split into a private occurrence during GFA
   materialization.

## Node vs occurrence semantics

The conversion has three distinct decisions:

| Decision | Granularity | Effect |
| --- | --- | --- |
| Frequency mask (`top=`, `max-occ=`) | Node | Removes all local occurrences of rejected high-copy syncmer nodes and bridges them by sequence. |
| Scaffold-chain support (`min-run=`) | Occurrence | Keeps only the exact syncmer occurrences that are members of retained SweepGA scaffold chains. |
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
for local selected paths. This keeps off-diagonal and indel-shifted colinear
handling in one implementation:

- indel-shifted but colinear runs can remain one scaffold chain;
- short off-chain singleton hits fail the scaffold mass threshold;
- no custom diagonal test in the GFA writer decides chain membership.

No SweepGA filtering semantics were changed for this integration.

## Path reconstruction

Filtering happens before shared topology is materialized, but path spelling is
preserved:

- frequency-masked nodes are removed from the anchor walk and the gap between
  neighboring retained anchors is filled from the selected input sequence when
  sequence files are available;
- scaffold-filtered off-chain syncmers are not deleted, only cloned/private, so
  they keep the same syncmer sequence on the path;
- exact blunt output still trims syncmer slices by the same overlap/clip logic
  used before this change;
- optional N cutting is the only intended path-spelling change, and only for
  gap DNA according to the configured `cut-n` policy.
