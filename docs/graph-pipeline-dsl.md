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
syng,k=63,s=8:crush,max-span=10k,max-traversals=128
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
  -> POVU-style collinear path-site discovery
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
crush,max-span=10k,max-traversal-len=10k,max-traversals=128,method=poa
```

Unknown parameters should be errors, not warnings. Silent ignoring would make
graph parameterization hard to reproduce.

## POVU Requirement For Fuller Crush

The current `src/resolution.rs` implementation uses the same reference-path
collinear-match site model as POVU's native Rust VCF extractor, then processes
leaf sites first and validates exact path-sequence preservation after every
replacement. This is enough to make `:crush` executable and testable.

The fuller scheduler should be POVU/flubble-driven when `povu-rs` exposes the
needed hierarchy API:

1. Decompose the current graph into flubbles/bubbles/tangles.
2. Build a parent/child hierarchy.
3. Select deepest non-overlapping sites.
4. Extract observed path traversals through each site.
5. Replace bounded sites exactly.
6. Recompute decomposition on the changed graph.

At the current pinned `povu-rs` revision, IMPG has native GFA-to-VCF support,
but not a public flubble hierarchy API. The current `:crush` implementation
therefore does not consume a public POVU hierarchy yet. We still need a
POVU-side API that exposes site IDs, parent/child relationships, entry/exit
handles, path-supported traversal boundaries, and enough node/edge identity to
schedule non-overlapping bottom-up batches.
