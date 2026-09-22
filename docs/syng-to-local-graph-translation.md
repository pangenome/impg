# Syng Evidence To Local Graph Translation

IMPG has two useful feature spaces that must interoperate:

```text
syng-syncmer-node       global implicit syncmer graph used by impg map/pack
gfa-segment            local explicit variation graph from syng/poa/seqwish/pggb/crush
```

Short-read evidence is cheap to collect in syng space. Local genotyping and
graph smoothing often need explicit local graph nodes. The render bundle must
therefore include a translation that projects syng-node evidence onto the local
graph without remapping reads.

## Required Invariant

The root coordinate system is still source sequence coordinates. Translation
must be path-step aware:

```text
syng node occurrence
  -> source sequence + bp position + orientation
  -> rendered input path + path step offset
  -> local graph node + orientation
```

A plain `syng_node_id -> gfa_node_id` map is not sufficient. Syncmers can occur
in repeats, local graph nodes can be shared by many paths, and graph smoothing
can replace a path interval with new nodes. The projection key must include the
rendered path context and source-coordinate overlap.

## Bundle Schema Extension

Keep the existing `RenderedPathRecord` and `StepTranslationRecord`. Add a second
table for cross-feature projection:

```rust
pub struct SyngToGraphProjectionRecord {
    pub rendered_path_id: u32,
    pub rendered_step: u32,
    pub local_feature_id: u32,
    pub local_orientation: char,
    pub source_start: u64,
    pub source_end: u64,
    pub syng_node_id: u32,
    pub syng_orientation: char,
    pub syng_source_start: u64,
    pub syng_source_end: u64,
}
```

Interpretation:

- `local_feature_id` is the local GFA segment/node ID on the rendered path.
- `source_start..source_end` is the source interval covered by that local graph
  step on this rendered path.
- `syng_node_id` is the global syng syncmer node ID.
- `syng_source_start..syng_source_end` is the source interval covered by that
  syncmer occurrence. It is normally `k=63` bp, adjusted for source strand.

The same local graph node may have many projection records, one per rendered
path context. That is intentional.

Manifest additions:

```json
{
  "feature_space": "gfa-segment",
  "source_feature_space": "syng-syncmer-node",
  "syng_to_graph": "syng-to-graph.bin",
  "syng_to_graph_tsv": "syng-to-graph.tsv"
}
```

## Build Algorithm

For every rendered path in a local graph bundle:

1. Get the rendered source interval from `RenderedPathRecord`.
2. Walk the source interval in the original syng index:

   ```rust
   full_syng.walk_path_range(source_path_id, source_start, source_end)
   ```

3. Parse the local GFA path walk and its per-step source offsets. This is
   already done by `collect_gfa_step_samples`.
4. Sweep both sorted interval streams:

   ```text
   syng occurrence intervals on source coordinates
   local graph step intervals on source coordinates
   ```

5. Emit a projection record for each overlap.

For exact syng-native bundles this table is identity-like: a syncmer node maps
to the matching syng GFA node on the same rendered path. For POA/seqwish/pggb
and future `:crush`, the table maps each syncmer occurrence onto the local graph
step or steps whose path-coordinate span overlaps that syncmer.

## Projection Of Pack Evidence

Given a syng pack:

```text
syng_node_id -> read_count
```

and `syng_to_graph` records for a local bundle, produce a local graph pack:

```text
local_feature_id -> projected_count
```

Initial scoring should be simple and deterministic:

```text
for each projection record:
  contribution = syng_count[syng_node_id] * overlap_bp / syncmer_len
  local_count[local_feature_id] += contribution
```

For genotyping, path-context records can also build haplotype feature vectors:

```text
rendered_path_id -> vector(local_feature_id weights)
```

This lets `genotype cos` compare read evidence to local graph haplotypes even
when the reads were only mapped once to the whole syng index.

## Ambiguity Rules

Projection is allowed to be many-to-many, but it must be explicit:

- Repeated syncmer occurrences produce multiple records with different
  `rendered_path_id` / source coordinates.
- A syncmer spanning a graph-node boundary contributes proportionally to both
  local nodes.
- Reverse-strand rendered intervals store source coordinates in forward source
  coordinates and set orientation fields explicitly.
- No deduplication across paths is done in the translation table. Downstream
  pack projection may choose path-aware or node-aggregate semantics.

## Implementation Plan

1. Extend `RenderTranslationTables` with `syng_to_graph: Vec<SyngToGraphProjectionRecord>`.
2. Add binary and TSV read/write support with backward-compatible versioning or
   bump the render bundle version.
3. In `render_local_graph`, build this table from the original full syng index,
   not the rendered FASTA alone.
4. Add `impg pack-project --render-bundle <bundle> --pack <syng.pack>` or wire
   automatic projection into `genotype cos --render-bundle` when the bundle
   feature space is `gfa-segment`.
5. Test on synthetic bubbles:
   - SNP bubble;
   - insertion/deletion bubble;
   - repeated syncmer copied twice;
   - reverse-strand interval;
   - a syncmer crossing a local graph node boundary.

The implementation must preserve Pan-SN names through `RenderedPathRecord`.
That is what lets the projected evidence support diploid/path-aware genotyping
instead of only anonymous node coverage.
