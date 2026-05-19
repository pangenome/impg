# Graph Pipeline DSL

IMPG graph output is moving toward a typed graph-construction pipeline instead
of one-off CLI flags for every engine and transformation.

The intended grammar is:

```text
gfa:<stage>[,<key=value>...][:<stage>[,<key=value>...]...]
vcf:<stage>[,<key=value>...][:<stage>[,<key=value>...]...]
```

Stages are separated by `:`. Parameters belong to the stage immediately before
them and are separated by commas:

```text
syng,k=63,s=8,seed=7:blunt
seqwish,min-match-len=70
pggb,window=20k
syng,k=63,s=8:crush,max-span=10k,max-traversals=10k
```

This is a graph-pipeline DSL, not a general shell pipeline. Each stage is a
typed graph producer or graph transform with declared inputs, outputs, and
parameters.

## Current Runtime Support

The parser accepts the staged grammar and runtime execution supports a small
set of typed producers plus the exact path-preserving `crush` transform.

Supported executable producers today:

```text
pggb
seqwish
poa
syng
```

Supported syng graph modes:

```text
syng
syng:blunt
syng:raw
syng,mode=blunt
syng,mode=raw
syng,k=63,s=8,seed=7:blunt
syng:blunt,k=63,s=8,seed=7
```

Supported legacy alignment prefixes remain:

```text
gfa:wfmash:seqwish
gfa:fastga:pggb
gfa:sweepga:seqwish
```

Supported partition/window forms remain:

```text
gfa:pggb:10k
gfa:pggb,window=10k
gfa:seqwish:20k
gfa:sweepga:fastga:pggb,window=20k
```

`crush` is a generic blunt-graph transform:

```text
input blunt graph
  -> POVU native flubble decomposition
  -> bottom-up leaf-site ordering
  -> exact path-preserving local replacement
  -> recompute site discovery after each replacement
  -> stop when no bounded sites remain
```

For syng, `crush` implies blunt input:

```text
syng:crush == syng:blunt -> crush
```

For seqwish and pggb, `crush` operates on the already-blunt local graph:

```text
seqwish:crush
pggb:crush
```

The important design point is that `crush` is not syng-specific. Syng, seqwish,
pggb, high-k seqwish, and future GBZ/handlegraph-backed renderers should all
feed a common path-embedded blunt graph into the same transform.

`syng:raw:crush` is rejected because the resolver requires blunt `0M` links.

## Stage Parameters

Each stage owns its own parameter namespace. The top-level CLI should parse the
pipeline into stages and then let each stage validate its own parameters.

Examples:

```text
syng,k=63,s=8,seed=7
seqwish,min-match-len=70,sparse-factor=0.001
pggb,window=20k
crush,max-span=10k,max-traversal-len=10k,max-traversals=10k,method=poa
```

Unknown parameters should be errors, not warnings. Silent ignoring would make
graph parameterization hard to reproduce.

## POVU Decomposition Boundary

`src/resolution.rs` does not implement its own flubble finder. It passes the
current blunt GFA to `povu-rs` and calls `NativeGfa::decompose_flubbles`, then
uses POVU's site IDs, parent/child relationships, entry/exit steps, and
bottom-up leaf ordering to schedule exact path-preserving replacement.

The crush loop is iterative because each exact replacement changes the local
graph topology and therefore can create, remove, or expose nested POVU sites.
One iteration means one candidate attempt, not one SPOA alignment step:

1. Decompose the current graph into flubbles/bubbles/tangles.
2. Build a parent/child hierarchy.
3. Select deepest non-overlapping sites.
4. Extract observed path traversals through each site.
5. Replace bounded sites exactly.
6. Recompute decomposition on the changed graph.

For debugging or batch processing of an already rendered blunt GFA, the same
transform is exposed directly:

```bash
impg crush -g local.blunt.gfa -o local.crushed.gfa
```

The defaults are sized for human panels (`512` candidate attempts and `10k`
path traversals per candidate). Traversal count is deliberately a high safety
rail. The main alignment budgets are `max-span`, `max-traversal-len`, and
`max-total-seq`, because a common allele represented by many haplotypes should
not fail only because many paths traverse it.
