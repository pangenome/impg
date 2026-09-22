# Syng Position Query Index

## Problem

The current repaired HPRCv2 syng index has correct regular positional sidecars:

```text
HPRC_r2_assemblies_0.6.1.syng.spos
HPRC_r2_assemblies_0.6.1.syng.pstep
```

Those sidecars fix query correctness for unsampled syncmer positions by
walking the GBWT to the nearest sampled checkpoint. They do not yet make
the first query as fast as it should be on a human-scale graph.

There were two remaining latency costs:

1. `SyngIndex::load` read the multi-GB `.pstep` payload
   into memory.
2. The first arbitrary target-position locate lazily builds an in-memory
   checkpoint hash from every sampled `.pstep` checkpoint:

```text
(signed_syncmer_node, absolute_gbwt_occurrence_rank) -> (path_id, bp_pos)
```

On HPRCv2 this means hundreds of millions of checkpoints. Building that
hash is the "materialize all positions" behavior that can dominate a
small query.

## Current Format

`.syng.pstep` is path-major:

```text
path_id -> sampled checkpoints ordered by path bp position
```

This is good for jumping into one source path range, because the query
path can binary-search or scan its own checkpoints and then walk forward
through the GBWT.

It is not good for target occurrence location. For a shared syncmer node,
query needs to enumerate GBWT occurrences and repeatedly ask:

```text
does this current (signed_node, occurrence_rank) equal a sampled checkpoint?
```

Today that question is answered by building the whole in-memory hash.

## Implemented Fix

`impg syng-repair` and normal `impg syng` construction now write
`.syng.spos` as a persistent checkpoint lookup sidecar:

```text
{prefix}.syng.spos
```

The sidecar is occurrence-major:

```text
signed_syncmer_node -> sorted occurrence_rank checkpoints
                         -> path_id, bp_pos
```

At query time, lookup becomes:

1. Memory-map `.syng.spos`.
2. Find the signed-node group by binary search or a small resident
   sampled dictionary.
3. Binary-search the group's occurrence ranks.
4. If present, return `(path_id, bp_pos)`.
5. If absent, advance the GBWT by one step and try again, up to the
   position sample rate.

This preserves the current correctness model but removes the global
first-query hash build.

## File Layout

The v1 `.syng.spos` layout is:

```text
header:
  magic, version, sample_rate, checkpoint_count

node table:
  signed_node ids with at least one checkpoint
  byte offsets or record offsets for each signed_node group

payload:
  per signed_node:
    fixed-width occurrence_rank, path_id, bp_pos records
```

The fixed-width payload is intentionally simple: it can be memory-mapped
and binary-searched directly by `(signed_node, occurrence_rank)` without
decoding unrelated nodes or building a global hash. The payload can later
be block-compressed or Elias-Fano encoded if disk footprint becomes the
dominant concern.

## Load Behavior

`SyngIndex::load` no longer eagerly reads all positional payloads:

- load `.1gbwt`, `.1khash`, `.names`, and `.meta` as today;
- memory-map `.pstep` and `.spos`;
- do not build the checkpoint hash at all when `.spos` is present.

For exact syng query, `.pstep` and `.spos` are the critical pair:
`.pstep` enters a path by coordinate, and `.spos` resolves sampled syncmer
occurrences to path coordinates.

## Repair And Distribution

The existing HPRCv2 repair does not need to be rerun for correctness, but
it needs one more repair/build step to rewrite `.syng.spos` into the new
occurrence-major format:

```bash
impg syng-repair -a HPRC_r2_assemblies_0.6.1.syng \
  --position-sample-rate 256 \
  --force \
  -t 32
```

That repair writes `.pstep` and `.spos` together. The positional sidecars
should then be uploaded together to S3 so other users do not pay the
repair cost:

```text
HPRC_r2_assemblies_0.6.1.syng.spos
HPRC_r2_assemblies_0.6.1.syng.pstep
```

Older `impg` binaries will not understand the new `.spos` sidecar. The
failure mode should be explicit: ask the user to update `impg` or rebuild
sidecars with the current binary.

## Testing

Correctness tests:

- exact locate of sampled and unsampled syncmer occurrences;
- query-region anchors match the old in-memory checkpoint hash;
- reverse-complement occurrences resolve to the same coordinates and
  correct strand semantics;
- short paths still sample first and terminal syncmers;
- normal query/partition/map loads require the new `.spos`; repair loads can
  tolerate a missing or stale `.spos` and rebuild it from `.pstep`.

Performance tests:

- HPRCv2 tiny query should not build a hundreds-of-millions-entry hash;
- first-query wall time should be dominated by GBWT/khash load, not
  positional checkpoint materialization;
- repeated query batches should not grow memory after the first query;
- report load time, mmap setup time, and checkpoint lookup count under
  `SYNG_EMIT_PROFILE=1`.
