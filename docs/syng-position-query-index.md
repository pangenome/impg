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

There are two remaining latency costs:

1. `SyngIndex::load` reads the multi-GB `.spos` and `.pstep` byte payloads
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

## Required Fix

Add a persistent checkpoint lookup sidecar, built by `impg syng-repair`
and by normal `impg syng` construction:

```text
{prefix}.syng.ckpt   # tentative name
```

The sidecar should be occurrence-major:

```text
signed_syncmer_node -> sorted occurrence_rank checkpoints
                         -> path_id, bp_pos
```

At query time, lookup becomes:

1. Memory-map `.syng.ckpt`.
2. Find the signed-node group by binary search or a small resident
   sampled dictionary.
3. Binary-search the group's occurrence ranks.
4. If present, return `(path_id, bp_pos)`.
5. If absent, advance the GBWT by one step and try again, up to the
   position sample rate.

This preserves the current correctness model but removes the global
first-query hash build.

## File Layout Direction

A simple v1 `.syng.ckpt` can be:

```text
header:
  magic, version, sample_rate, checkpoint_count

node table:
  signed_node ids with at least one checkpoint
  byte offsets or record offsets for each signed_node group

payload:
  per signed_node:
    delta-coded occurrence_rank
    delta-coded path_id or raw path_id
    delta-coded bp_pos within path
```

This is already enough to be queryable without decoding unrelated nodes.
The payload can later become block-compressed, Elias-Fano encoded, or
partitioned by signed node. The important first step is the access shape:
node-local lookup without whole-index materialization.

## Load Behavior

`SyngIndex::load` should stop eagerly reading all positional payloads:

- load `.1gbwt`, `.1khash`, `.names`, and `.meta` as today;
- load only small sidecar metadata immediately;
- memory-map `.pstep`, `.spos`, and `.ckpt` payloads;
- do not load `.spos` unless a command actually needs sampled
  node-to-position records;
- do not build the checkpoint hash at all when `.ckpt` is present.

For exact syng query, `.pstep` and `.ckpt` are the critical pair.
`.spos` is useful as a sampled node-to-position sidecar, but it is not the
final answer for exact query positions.

## Repair And Distribution

The existing HPRCv2 repair does not need to be rerun for correctness, but
it will need one more repair/build step after `.syng.ckpt` exists:

```bash
impg syng-repair -a HPRC_r2_assemblies_0.6.1.syng \
  --position-sample-rate 256 \
  --force \
  -t 32
```

That repair should write `.spos`, `.pstep`, and `.ckpt` together. The
three positional sidecars should then be uploaded together to S3 so other
users do not pay the repair cost:

```text
HPRC_r2_assemblies_0.6.1.syng.spos
HPRC_r2_assemblies_0.6.1.syng.pstep
HPRC_r2_assemblies_0.6.1.syng.ckpt
```

Older `impg` binaries will not understand the new sidecar. The failure
mode should be explicit: either ignore `.ckpt` and fall back to the old
hash build, or error with a message asking the user to update `impg`.

## Testing

Correctness tests:

- exact locate of sampled and unsampled syncmer occurrences;
- query-region anchors match the old in-memory checkpoint hash;
- reverse-complement occurrences resolve to the same coordinates and
  correct strand semantics;
- short paths still sample first and terminal syncmers;
- missing `.ckpt` falls back or errors according to the chosen policy.

Performance tests:

- HPRCv2 tiny query should not build a hundreds-of-millions-entry hash;
- first-query wall time should be dominated by GBWT/khash load, not
  positional checkpoint materialization;
- repeated query batches should not grow memory after the first query;
- report load time, mmap setup time, and checkpoint lookup count under
  `SYNG_EMIT_PROFILE=1`.
