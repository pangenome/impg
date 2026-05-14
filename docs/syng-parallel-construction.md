# Syng Parallel Construction Design

## Goal

Build a syng index from many genomes without being limited by one serial
`kmerHashAdd` / `syngBWTpathAdd` stream for the whole input. The construction
must produce a normal syng index:

- `{prefix}.1khash`
- `{prefix}.1gbwt`
- `{prefix}.syng.names`
- `{prefix}.syng.spos`
- `{prefix}.syng.meta`

The output does not need to assign the same numeric syncmer node ids as the
legacy one-pass build. It must assign ids consistently across all files and
produce equivalent query/map behavior.

## Mathematical Model

For each genome path `p`, syncmer extraction produces an ordered stream:

```text
S_p = [(x_1, strand_1, pos_1), ..., (x_n, strand_n, pos_n)]
```

where `x_i` is the canonical full syncmer sequence, `strand_i` records whether
the occurrence matched the canonical orientation or reverse complement, and
`pos_i` is the bp coordinate on `p`.

The global syncmer dictionary is:

```text
K = union_p { x_i : (x_i, strand_i, pos_i) in S_p }
```

Construction chooses one deterministic bijection:

```text
phi: K -> {1, 2, ..., |K|}
```

Every path is then encoded as signed global node ids:

```text
P_p = [(strand_i * phi(x_i), pos_i)]
```

The GBWT and sampled-position sidecar are functions of the ordered collection
of `P_p` streams. A parallel build is correct if all shards use the same `phi`
and the final path order is the chosen global path order.

## Production Algorithm

The clean end-state is a two-pass parallel dictionary plus per-node reduce:

1. Enumerate all input paths and assign stable global path ids.
2. In parallel, extract canonical syncmers and build shard-local unique sets.
3. Sort/deduplicate the union and assign global ids.
4. In parallel, stream inputs again and map syncmers through the global
   dictionary.
5. For each shard, emit per-node occurrence records and sampled-position
   records using global node ids and global path coordinates.
6. Reduce per-node occurrence records into the final GBWT node structures.
7. Sort/group sampled-position records and write `.syng.spos`.

Step 6 is the only part not exposed by the current syng C API. The existing
API has one mutable `syngBWTpathAdd` stream, but no `syngBWTmerge` or builder
from per-node occurrence lists. That API should be added only after the
dictionary split is tested.

## Implemented First Step

The first implementation stage is intentionally conservative:

1. Build the global syncmer dictionary in parallel.
2. Construct a syng `KmerHash` from that frozen dictionary.
3. Replay paths through the existing serial GBWT writer, using dictionary
   lookup instead of adding new syncmers.

This validates the global-id model and allows tests to compare the new two-pass
path against the legacy one-pass path. It also removes hash-table growth and
syncmer insertion from the serial GBWT phase. The final GBWT insertion remains
serial until a true per-node reduce API is added.

This path is exposed as `impg syng --parallel-dictionary` for both FASTA and
AGC inputs.

## Correctness Requirements

- Syncmer parameters are read from `.syng.meta` and written with every build.
- Packed canonical syncmers use the same canonical orientation as syng's
  `kmerHashAdd`.
- The frozen dictionary rejects duplicate packed syncmers.
- Replay fails if any extracted syncmer is absent from the frozen dictionary.
- `.syng.spos` sampling uses global path ids and positions, so it remains
  merge-compatible.

## Next Step: True GBWT Reduce

The next implementation layer uses a new internal representation:

```text
node_id -> ordered incoming/outgoing occurrence lists
```

For a fixed global path order, every directed edge occurrence receives a
deterministic rank within that edge. Each node can then be reduced
independently:

- directory counts come from incoming edge ranks on that side;
- BWT run lists are ordered by the opposite side's incoming edge ranks;
- the run symbols are edge ids on the side being built.

This is exposed as `impg syng --parallel-reduce`. Rust extracts ranked path
occurrences and reduces each node independently; a C helper in `syngbwt3.c`
installs the resulting edge directories and run lists directly into syng's
`NodeSide` / `Rskip` representation. The implementation is currently in-memory:
the next scaling step is to spill shard occurrence records and externally sort
by `(node_id, side)` and `(directed_edge, occurrence_order)` for HPRC-scale
construction without retaining all path occurrences at once.
