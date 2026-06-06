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
syng,k=63,s=8:crush,max-traversal-len=1k,max-traversals=10k
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
syng:mask,top=0.001,max-occ=500
syng:mask,top=0.001,max-occ=500,freq-run=10,freq-span=1k
syng:mask,legacy-freq-mask=true
syng:filter,top=0.001,max-occ=500
syng:cut-ns
syng:filter,cut-ns=true,cut-n-min-run=10
syng:nomask
syng:nofilter
```

For syng GFA extraction, `syng` defaults to a small raw-layer high-frequency
syncmer filter: the top 0.05% most frequent local syncmer nodes are selected
before bluntification, but the default policy is occurrence-level. Unsupported
high-frequency occurrences are private-split, while occurrences in supported
high-frequency runs or exact spans stay shared. `freq-run` / `high-freq-run`
controls the minimum run length in syncmers and defaults to `10`; `freq-span` /
`high-freq-span` controls exact-sequence span rescue in bp and defaults to
`1000`. The same `freq-run` / `freq-span` rescue applies to the
spectrum-selected dispersed scaffold-glue source. Use
`freq-run-aware=false` or `legacy-freq-mask=true` for the historical node-level
removal behavior. Rare repeated-copy local syncmer contexts are still split
rather than emitted as one globally shared graph node. This is a
graph-materialization policy, separate from the query seed filter.

`syng:cut-ns` is an optional materialization policy for assembly gaps. It drops
N-runs from fetched inter-syncmer DNA and splits emitted GFA paths at those
breaks, preventing large path-private `N` segments from dominating graph
layouts. `cut-n-min-run` controls the minimum ambiguous run length and defaults
to `1`.

Supported legacy alignment prefixes remain:

```text
gfa:wfmash:seqwish
gfa:fastga:pggb
gfa:sweepga:seqwish
```

Supported query-sequence preprocessing stages:

```text
gfa:cut-n=100:pggb
vcf:cut-n=100:seqwish
```

`cut-n=<bp>` is explicit and has no bare default. It clips only terminal
`N`/`n` runs from each query-extracted sequence when the run length is at
least the requested threshold before graph construction. Internal N-runs and
terminal runs shorter than the threshold are left unchanged. If both ends clip
away the full extracted sequence, that sequence is omitted from the local graph.
This is separate from syng's `syng:cut-ns` materialization policy, which cuts
N-runs from fetched inter-syncmer gap DNA.

Supported partition/window forms remain:

```text
gfa:pggb:10k
gfa:pggb,window=10k
gfa:seqwish:20k
gfa:sweepga:fastga:pggb,window=20k
```

The same staged spelling can ask an alignment-derived graph engine to run the
path-preserving crush machinery after graph induction:

```text
gfa:pggb:crush
gfa:seqwish:crush
```

For these alignment engines, `crush` runs on the induced blunt GFA and the final
graph is emitted with the default self-loop run normalization plus `Ygs`
sorting/finalization.

`crush` is a generic blunt-graph transform:

```text
input blunt graph
  -> POVU native flubble decomposition
  -> maximal eligible non-overlapping site frontier
  -> parallel exact path-preserving local replacement
  -> single graph rewrite for the frontier
  -> recompute site discovery after each round
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
crush,method=auto,max-median-traversal-len=1k,max-traversal-len=10k
crush,method=biwfa,max-median-traversal-len=10k,polish-max-median-traversal-len=256
crush,method=allwave,k-nearest=5,k-farthest=2,random-fraction=0.05,mash-k=17
crush,method=sweepga,aligner=fastga,kmer-frequency=42,no-filter=true
```

Unknown parameters should be errors, not warnings. Silent ignoring would make
graph parameterization hard to reproduce.

## POVU Decomposition Boundary

`src/resolution.rs` does not implement its own flubble finder. It passes the
current blunt GFA to `povu-rs` and calls `NativeGfa::decompose_flubbles`, then
uses POVU's site IDs, parent/child relationships, entry/exit steps, and
root-path spans to schedule exact path-preserving replacement.

The crush loop is iterative because each exact replacement changes the local
graph topology and therefore can create, remove, or expose nested POVU sites.
One round resolves a frontier of independent sites in parallel; it is not one
SPOA alignment step:

1. Decompose the current graph into flubbles/bubbles/tangles.
2. Build a parent/child hierarchy.
3. Build candidates for all observed path-supported sites.
4. Select a maximal non-overlapping frontier: choose the largest eligible site
   under the budgets, and descend only when a containing site is too large.
5. Resolve the selected sites in parallel.
6. Rewrite the graph once and validate exact path-sequence preservation.
7. Recompute decomposition on the changed graph.

For debugging or batch processing of an already rendered blunt GFA, the same
transform is exposed directly:

```bash
impg crush -g local.blunt.gfa -o local.crushed.gfa
```

The defaults are sized for human panels (`1` frontier round and `10k` path
traversals per candidate). Traversal count is deliberately a high safety rail.
The main
alignment budgets are `max-median-traversal-len` (default `1k`),
`max-traversal-len` (default `10k`), and `max-total-seq`, because a common
allele represented by many haplotypes should not fail only because many paths
traverse it. `max-span` is optional and disabled by default; when set, it caps
the span on the POVU root path, currently the first GFA path, so it is a rooted
coordinate guard rather than the main runtime budget.

`method=auto` uses direct global/end-to-end SPOA only for bubbles whose longest
traversal is at most `auto-spoa-max-len` (default `2k`) and whose traversal
count / total sequence size are still inside the direct replacement budgets.
Short but very high-copy bubbles are routed to the scalable pairwise induction
path instead of being selection-guarded. Remaining bounded bubbles use pairwise
graph induction, with seqwish using a high exact-match length by default
(`seqwish-k=311`) so human repeats are not glued through short off-diagonal
matches before the small SPOA polish pass. This is currently a human-repeat
default chosen just above Alu length; the long-term default should be derived
from the expected identity / repeat model for the local sequence set. Direct
SPOA/POASTA replacements are exact path-sequence validated and are not rejected
by the local graph-layout quality heuristic; pairwise-induced replacements still
must pass local quality guards. All methods pass the round-level visual-tail
guard: the score is dominated by long path white-space bridges (`ws-p99` and
`ws-max`) and only lightly penalizes total white space, path steps, and link
count, so useful local crushing is not rolled back merely because it splits
paths more finely. `method=poasta` constructs POASTA in
global/end-to-end mode with `GapAffine2Piece` costs from `poa-scoring`, and impg
normalizes POASTA's clipped `W` walk GFA export before exact sequence
validation. It remains an explicit experimental method because POASTA 0.1.0's
public API does not give impg a stable exact numeric score contract for
independently asserting long-gap two-piece affine costs; current tests validate
structural wiring plus exact path preservation.

Pairwise-induced replacements are filtered before seqwish induction. `seqwish-k`
sets the minimum exact-match scale used by seqwish transitive closure, and
`min-map-length` / `replacement-min-map-length` sets the minimum PAF mapping
length kept by the SweepGA filter. A value of `min-map-length=0` follows
`seqwish-k`, so the default graph-induction mapping scale is also `311`.
`min-identity` / `replacement-min-identity` can require a minimum mapping
identity before a pairwise alignment is allowed to seed local graph induction.
For `method=sweepga`, `kmer-frequency=0` is the default and resolves per
candidate to at least `1000` and otherwise `10 * traversal_count`. This differs
from whole-genome FastGA defaults deliberately: local bubble crushing needs
repeated seeds to survive, while the downstream plane-sweep/scaffold filter and
`seqwish-k` control graph glue.

Syng-native graph extraction has an earlier raw-layer mask before bluntification
or crush. `gfa:syng:mask,min-run=3,freq-run=10:crush` private-splits locally
shared syncmer occurrences that never appear in a supported shared run, and
private-splits high-frequency occurrences from explicit frequency selection or
the spectrum-glue source unless they are rescued by the high-frequency
`freq-run` or `freq-span` policy. Private splits preserve the selected path
spelling while preventing isolated shared syncmers from becoming graph glue.
`legacy-freq-mask=true` restores the older behavior that removed
selected high-frequency syncmer nodes and bridged them in cis using source
sequence.

`method=biwfa` is the coarse condenser path. It aligns every selected POVU
bubble traversal end-to-end against the longest traversal with BiWFA, induces a
path-preserving in-memory column graph from those root alignments, then runs one
small-scale SPOA polish pass over the induced replacement graph. This avoids a
seqwish transitive-closure step inside bounded bubble replacement. The two
scales are intentionally separate: `max-median-traversal-len` controls which
biological bubble traversals may be induced by BiWFA, while
`polish-max-median-traversal-len` controls the tiny STR/indel tangles we are
willing to clean inside that induced graph.

`method=allwave` is the many-sequence BiWFA path. It gives every selected
POVU bubble traversal to AllWave at once, uses AllWave's tree/kNN pair
selection and strand-specific Mash orientation detection, emits oriented PAF,
and induces the replacement graph through seqwish. This is useful when a bubble
has many homologous traversals and root-star alignment is too biased. The
replacement is still exact-path validated after graph induction.

`method=sweepga` uses SweepGA as the replaceable alignment provider before the
same seqwish induction and validation step. It is intentionally separate from
`method=allwave`: SweepGA has different filtering/scaling semantics, and its
knobs can be used to leave harder bubbles in the graph for a later small SPOA
crush pass. The shared pair-sampling knobs are `k-nearest`, `k-farthest`,
`random-fraction`, and `mash-k`; SweepGA-specific knobs include `aligner`,
`kmer-frequency`, `min-aln-length`, `map-pct-identity`, and `no-filter`.
