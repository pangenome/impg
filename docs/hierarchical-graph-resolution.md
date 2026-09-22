# Hierarchical Graph Resolution

IMPG can already build local sequence graphs from pangenome intervals, including
graphs induced from sparse de Bruijn or syncmer structures such as syng. Once
those graphs are blunt variation graphs, the next missing control surface is
not primarily VCF formatting. It is graph resolution: finding the regions of a
graph that should be locally realigned, smoothed, compressed, simplified, or
only summarized at a coarser topology-preserving scale.

This note proposes a hierarchical bubble decomposition layer after local graph
construction. The layer should identify clean bubbles and more general flubbles,
organize them into larger sites, and choose an operation for each region:
lossless pass-through, local partial order realignment, bidirectional-WFA pair
compression, representative traversal collapse, symbolic coarse representation,
VCF emission, or genotype/deconvolution against selected haplotypes.

The important shift is that VCF is only one view over the resolved graph. The
same hierarchy should also drive bubble-guided smoothing for pggb/seqwish/syng
graphs, compression of poor local alignments, graph cleanup, and local
haplotype inference.

## Problem

Current graph construction can leave us with the wrong local resolution. De
Bruijn-style graph induction, syncmer graphs, poor pairwise alignments, or
under-smoothed seqwish output may all produce tangles that are topologically
correct at a large scale but too noisy or too fragmented inside local bubbles.
Conversely, trying to fully align very large, high-divergence bubbles can be
expensive, lossy, or misleading.

The core issues are:

- Local pangenome graphs are often nested: small variants sit inside insertion
  alleles, repeat copies, or alternative structural haplotypes.
- Existing smoothxg-style smoothing chooses blocks from a 1D graph sort and
  heuristic block decomposition. That works, but it is not biologically or
  topologically targeted: the blocks are not necessarily bubbles.
- Bubble decomposition gives exact local work units. A bubble has common
  anchors at both ends, so local alignment has a clear reason to align these
  sequences together.
- Small bubbles can be locally realigned into a clean partial order graph.
  Larger bubbles may need compression or representative traversal selection
  rather than full base-level smoothing.
- VCF emission, simplified GFA, and genotyping need the same notion of
  resolution, or the outputs will disagree.
- Genotyping needs representative haplotypes that can explain read coverage and
  linkage, not just syntactically valid VCF alleles.

The desired behavior is a tunable path:

```text
local GFA
  -> normalize/clean graph
  -> detect bubbles/flubbles and cyclic tangles
  -> cluster into hierarchical flubbles/sites
  -> choose per-site operation by size/divergence/topology
  -> smooth, compress, simplify, or preserve each site
  -> emit resolved GFA plus optional VCF/sidecars
  -> genotype/deconvolve sample evidence against representatives
```

## Goals

- Support multiple resolution levels from atomic bubbles to large SV/CNV sites.
- Use bubble/flubble decomposition as the primary unit of graph smoothing and
  compression.
- Use POVU-style bubble, flubble, and cycle structure to find tight tangles
  that should be resolved or summarized as one local unit.
- Produce locally consistent partial order graphs for bubbles below a practical
  sequence/complexity threshold.
- Preserve coarse sparse-de-Bruijn/syng topology above that threshold instead
  of forcing a fragile full alignment.
- Keep source path identity and path intervals through every transformation.
- Select representative traversals that cover common haplotype space while
  bounding allele count.
- Emit VCF records at the requested resolution as a downstream view, with clear
  provenance back to graph sites and lower-level bubbles.
- Simplify or clean GFA graphs with the same site hierarchy used for local
  smoothing and compression.
- Expose representative haplotypes as feature vectors for genotyping,
  deconvolution, and local haplotype inference.

Non-goals for the first implementation:

- Perfect decomposition of every cyclic assembly graph.
- A new whole-genome variant caller independent of IMPG's local graph and pack
  abstractions.
- Replacing POVU immediately. The first target can wrap or preprocess GFA for
  POVU and then progressively move more logic into IMPG.
- Reimplementing smoothxg wholesale. We should learn from its local POA
  smoothing, but the work units should be graph-topological bubbles, not
  heuristic chunks from a 1D sorted graph.

## Hierarchical Model

Use a hierarchy of typed graph regions:

```text
level 0: atomic graph features
  segment, oriented segment, edge, path step, syncmer node

level 1: atomic bubbles
  single-entry/single-exit alternative traversals, usually small variants

level 2: compound flubbles / variant sites
  adjacent, nested, or weakly interacting bubbles that should often be operated
  on together, such as an insertion with internal SNPs or a deletion plus
  nearby breakpoint polymorphism

level 3: complex loci / blocks
  large SV, CNV, duplicated, or repetitive regions where representative
  traversal selection, compression, or topology preservation is required
```

"Bubble" should be reserved for clean source/sink subgraphs when possible.
"Flubble" is useful for an almost-bubble: a local component with recognizable
entry/exit anchors and path traversals, but with internal nesting, weak cycles,
copy-number variation, or multiple plausible boundaries. Resolution operations
may use either, but the model should record boundary confidence and
decomposition type.

The hierarchy is usually a tree along a chosen reference path, but the data
model should allow a DAG. Repetitive loci can cause one lower-level bubble to
participate in multiple larger candidate sites, and a strict tree would force
early arbitrary choices.

## Proposed Data Model

The decomposition layer should produce a portable sidecar object next to the
GFA/VCF, even if the first implementation keeps it internal.

```text
GraphResolutionIndex
  graph_id
  graph_build_engine
  reference_paths
  feature_space
  normalization_steps
  sites[]
  traversals[]
  representatives[]
  provenance_maps[]
```

Each site describes a graph region at one level:

```text
Site
  id
  level
  kind = atomic_bubble | flubble | cnv_block | repeat_block | component
  parent_ids
  child_ids
  entry_handles
  exit_handles
  node_set
  edge_set
  reference_path
  reference_interval
  boundary_confidence
  acyclic
  copy_variable
  path_count
  traversal_count
  notes
```

Each traversal is an observed or inferred path through a site:

```text
Traversal
  id
  site_id
  oriented_handles
  sequence_hash
  sequence_length
  source_path_intervals
  child_site_alleles
  copy_number_vector
  frequency_or_path_support
  feature_vector_id
```

Each representative is a selected allele or polytype at a chosen resolution:

```text
Representative
  id
  site_id
  level
  method = exact | medoid | k_medoids | centroid | k_centroids | polytype | reference
  member_traversal_ids
  representative_traversal_id
  assigned_weight
  max_member_distance
  feature_vector_id
  emitted_allele_id
  deconvolution_children
```

Important invariants:

- A representative must map back to one or more graph traversals.
- Smoothed subgraphs, compressed paths, VCF alleles, simplified GFA paths, and
  genotype candidates must all reference the same site and representative IDs.
- The reference allele is explicit, even when it is not the most common
  traversal.
- Source path identity is never discarded. It may be summarized, but the full
  path membership must remain recoverable for audit and training data.
- Exact graph rendering does not need a full old-handle-to-new-handle
  translation table as long as all output paths are embedded correctly in the
  new graph. Path coordinates can be recovered from those emitted paths.
- Approximate, representative, or otherwise lossy graph rewrites must produce
  enough coordinate/provenance metadata to explain how old path intervals map
  to new graph paths, representative paths, or symbolic blocks.

## Decomposition Algorithm

The first implementation should be deterministic and conservative.

1. Parse and normalize the GFA.

   Build a handle graph view with segment lengths, edges, path walks, and path
   step positions. Optionally compact linear chains, remove exact duplicate
   segments, and annotate tips or low-support graph pieces. If the rewritten
   graph still embeds exact source paths, the path embeddings are the primary
   coordinate system. A separate old-to-new handle provenance map is only
   required once the operation changes path sequence, collapses alternatives,
   or emits synthetic representatives.

2. Choose anchor paths and coordinate frames.

   Prefer the requested reference path/range when `query` has one. Also keep
   all source paths from syng/Pan-SN metadata. Site coordinates should be
   reported in reference coordinates when available and graph-local coordinates
   otherwise.

3. Detect atomic bubbles and cyclic tangles.

   Start with superbubble-like detection on the directed bidirected graph after
   orienting around the reference path. Boundaries should be graph handles with
   multiple path-consistent traversals between them. Reject or downgrade sites
   with ambiguous entry/exit, too many internal cycles, or insufficient path
   support. Also record dense local cycles and knots as tangle candidates even
   when they are not clean bubbles; these are candidates for compression or
   coarse preservation.

4. Enumerate path-supported traversals.

   Use existing graph paths as the default traversal source. Do not enumerate
   all possible walks in a repetitive graph. For small acyclic bubbles, optional
   recombination across child bubbles can be allowed, but it should be gated by
   `--max-traversals` and clearly labeled as inferred rather than observed.

5. Build larger flubbles/sites.

   Cluster atomic bubbles when any of the following holds:

   - nested containment in the graph
   - adjacent reference spans within `--site-merge-distance`
   - high path co-occurrence or linkage disequilibrium across source paths
   - shared internal nodes or edges
   - breakpoints of one SV overlap the reference span of another bubble
   - copy-number traversals cannot be represented by independent child bubbles

   Build candidate parents bottom-up, then score whether a parent should exist.
   Useful parent scores include reference span, graph node count, path count,
   entropy of child allele combinations, traversal count, and copy-number
   dispersion.

6. Assign resolution levels and operations.

   Each site can be operated on at multiple levels:

   ```text
   atomic       resolve child bubbles independently
   site         operate on a compound flubble as one local graph block
   representative collapse to selected traversal representatives
   block        preserve larger sparse topology with symbolic/sidecar summary
   ```

   The default should be conservative: atomic for simple sites, site-level for
   obvious nested SV/indel structures, local smoothing for tractable bubbles,
   and representative/block mode only when full local smoothing would be too
   expensive or would collapse important large-scale topology.

7. Select representatives.

   Collapse duplicate traversals first, then select exact or approximate
   representatives according to the requested resolution and bounds.

8. Materialize the resolved view.

   Replace each operated-on site with either a smoothed subgraph, a compressed
   representative subgraph, or the original topology plus annotations. Then
   lace the transformed sites back into the surrounding graph using the site
   entry/exit handles.

## Bubble-Guided Smoothing And Compression

This hierarchy should become the block selector for smoothing. The current
smoothxg-style path in IMPG inherits the pggb idea: sort the graph into one
dimension, cut blocks with heuristics, run POA per block, and lace the result
back. That is useful, but it is not the cleanest model once we have a bubble
decomposition. A clean bubble already says: these traversals share a start and
end and should be compared as alternatives in this local context.

The proposed smoothing engine should be bubble-guided:

```text
input graph
  -> POVU / internal bubble+flubble/cycle decomposition
  -> classify each site by span, traversal count, divergence, cycles, alignability
  -> choose operation
  -> run local POA / BiWFA compression / representative collapse / passthrough
  -> lace the replacement graph back between the same boundary handles
  -> repeat on the changed graph until no eligible bubbles remain or budget is hit
  -> write resolved graph plus site metadata and optional coordinate provenance
```

For syng, this is a direct extension of the existing model. The syncmer graph
already gives a sparse topology over path-supported homology. Resolving a small
bubble to a tighter local graph and lacing it back compresses that sparse graph
without needing to materialize or globally realign everything. After one bubble
is replaced, the surrounding graph changes; that replacement may become a
fragment inside a larger bubble, which can then be considered by the next
iteration. The algorithm is therefore naturally bottom-up and hierarchical.

Suggested operation classes:

```text
passthrough
  keep the original subgraph; annotate why it was not touched

local-poa
  extract all path-supported traversal sequences through a small bubble and
  build a local partial order graph

pair-compress
  align traversal pairs with bidirectional WFA inside the bubble anchors and
  induce a compact local graph, seqwish-style but using the known bubble
  context rather than global mapping

representative
  select medoids/k-medoids/polytypes and collapse low-value redundant
  traversals onto them, preserving provenance

symbolic-block
  keep large-scale topology and emit a symbolic summary/sidecar for downstream
  VCF or genotyping
```

Initial thresholds can be simple and explicit:

```text
local-poa if:
  reference_span <= 1000 bp
  max_traversal_length <= 3000 bp
  traversal_count <= 128
  graph is acyclic or weakly cyclic after path support filtering

pair-compress if:
  max_traversal_length <= 10000 bp
  traversal_count is too high for all-vs-all POA but pair selection is bounded
  entry/exit anchors are clean

representative/block if:
  traversal_count, copy-number entropy, cycle count, or span exceeds budget
```

The exact numbers should be configurable and tuned empirically. The guiding
principle is that we want base-level local consistency up to a practical kernel
size, perhaps 1-3 kb by default, and then we deliberately retain coarser
topology above that scale rather than pretending a forced alignment is more
truthful.

### Bottom-Up Bubble Compression

The core loop can be simple:

```text
while true:
  find current bubbles/flubbles/tangles
  rank smallest and cleanest candidates first
  for each candidate within budget:
    extract observed path traversals
    build a replacement graph with the selected operation
    verify exact path embedding for exact modes
    lace replacement between candidate entry/exit handles
  stop when no candidate changes the graph or global budget is reached
```

This makes the hierarchy emerge from the graph instead of being fixed at the
start. A local replacement can simplify a nested structure and expose the next
larger bubble. The next pass can then operate on that larger unit, using the
already-compressed internal graph as one of its traversed fragments.

Budgeting is part of correctness, not just performance. The engine should have
hard stop conditions:

- maximum bubble reference span, for example 10 kb by default
- maximum traversal length, total extracted sequence, and traversal count
- maximum cycle/tangle score
- maximum pairwise-alignment count
- maximum wall time or memory per site
- maximum graph growth ratio after replacement

When a site exceeds budget, the correct first approximation is to leave its
current topology in place, optionally after compressing eligible internal
components. A 10 kb bubble with no shared syncmers in the middle is likely a
real structural difference or a region where base-level homology is not worth
forcing. Preserving the sparse syng/de-Bruijn topology is often more honest
than spending unbounded time trying to align it.

The compression schedule should avoid all-vs-all alignment by default. Useful
bounded schedules include:

- guide-tree or spanning-tree alignments over traversal sketches
- reference-plus-nearest-neighbor alignments
- random chords between traversal clusters
- extra alignments for high-support or long-distance outliers
- iterative refinement only while replacement quality improves

This is enough to get most of the benefit without quadratic blowups. Exact
small sites can use denser alignment, while large or repetitive sites should
fall back to representative or passthrough modes.

### Resolution Policy

The policy should be explicit about what scale is being represented:

```text
below the local alignment kernel
  make a clean local partial order graph when the traversal set is bounded and
  alignable

near the kernel boundary
  compress with bounded pairwise alignments or choose representatives, keeping
  exact path membership

above the kernel boundary
  preserve sparse de Bruijn / syng / large-SV topology and expose it as a
  block-level object
```

This avoids the main failure mode of global smoothing: forcing base-level
homology into regions where the input graph is really expressing larger-scale
topology. A syng or sparse de Bruijn graph can be exactly the right object at
large scale while still being too coarse or fragmented inside a small bubble.
Hierarchical resolution should improve the small local pieces without
destroying the larger topology that made the sparse graph useful.

The alignability model can start as a deterministic classifier over:

- maximum and median traversal length
- traversal count and duplicate sequence count
- length dispersion across traversals
- cycle/tangle count inside the site
- anchor uniqueness and boundary confidence
- repeat/self-similarity sketch score
- pairwise sketch distance between traversals
- path support and Pan-SN source diversity

The output is an operation choice plus a reason string. That reason should be
emitted in `resolution.json`, because users will need to know whether a site was
smoothed, pair-compressed, represented by medoids, or deliberately left as
coarse topology.

### BiWFA Inside Bubbles

Bidirectional WFA is attractive here because the bubble boundaries are known.
We do not need to discover whether two sequences belong together; the graph
topology has already asserted that they are alternative traversals between the
same anchors. That lets us run pairwise compression in a strongly constrained
local coordinate frame:

- extract traversal sequences from entry to exit
- select pairs by path support, distance, or k-nearest traversal sketches
- run BiWFA for each selected pair
- left-normalize / compress indel-equivalent alignments locally
- feed the alignments to seqwish/cseqwish-like induction or a smaller internal
  alignment-to-graph builder
- stitch the induced subgraph back between entry and exit

The main risk is local over-compression of indel/repeat structure. We should
therefore preserve path-step provenance and keep the operation bounded to
regions where the alignment is credible. For STRs, satellites, or high-copy
CNVs, representative or symbolic block modes may be more honest than forcing
base-level compression.

The important distinction from whole-graph alignment is that the bubble gives
the coordinate frame. We are not asking whether two sequences should align
somewhere in the graph; we are asking how to express alternatives between the
same entry and exit handles. That makes pair compression a local graph
resolution operation rather than a discovery operation.

## Representative Traversals, Centroids, And Polytypes

Representative selection should optimize coverage of common haplotype space
under allele-count and divergence constraints.

A practical objective:

```text
maximize covered_path_weight
penalize representative_count
penalize member_distance_to_representative
penalize loss of copy-number-distinguishing features
force include reference traversal
force include long or high-impact outliers above distance threshold
```

Distance can start simple:

- edit distance between traversal sequences for small sites
- Jaccard or weighted cosine distance over graph feature vectors for large sites
- child-allele Hamming distance for nested bubble vectors
- copy-number vector distance for duplicated regions
- breakpoint distance for SV alleles

Representative modes:

- `exact`: one allele per unique observed traversal sequence or child-allele
  vector.
- `medoid`: choose the observed traversal with minimum weighted distance to
  all assigned traversals.
- `k-medoids`: choose up to `k` observed traversals. This is safer than a true
  centroid because every emitted allele exists in the graph.
- `centroid` / `k-centroids`: compute one or more consensus traversals or
  sequences. This is useful for graph cleaning but should be marked as
  synthetic if emitted.
- `polytype`: represent a repeated or CNV locus by a compact type vector, such
  as copy count plus ordered or unordered child allele composition.

Polytypes are important for CNV and repeats. A two-copy allele and a three-copy
allele may share most nodes but imply different dosage and should not collapse
unless the requested resolution explicitly allows it. A polytype should be able
to say:

```text
site=C4_like_block
copy_number=3
copy_units=[A, B, B]
representative_path=sample42#1#chr6:...
assigned_paths=...
```

Exact and medoid representatives can become paths in a simplified graph. They
can also become VCF sequence alleles or symbolic alleles when the selected site
is being exported as variants. Centroids and polytypes may be useful for graph
cleaning or compression, but should be marked as synthetic if emitted as
sequence paths because they may not correspond to a real source haplotype.

## Resolved Graph Output

The primary output of hierarchical resolution should be a graph, not a VCF.
VCF, genotype candidates, and traversal FASTA are views over the same resolved
site model.

Possible graph outputs:

```text
resolved.gfa
  graph after bubble-guided smoothing/compression

resolution.json
  site hierarchy, operations, thresholds, and provenance maps

coordinate_map.tsv / coordinate_map.bin
  optional old path interval -> new path interval, representative assignment,
  or symbolic block assignment for approximate/lossy operations

site_traversals.fa
  traversal sequences used for local POA/BiWFA/representative selection
```

For exact path-preserving rendering, `resolved.gfa` and `resolution.json` are
sufficient. The emitted paths are the coordinate system. A coordinate map is
only needed when the output graph stops being an exact embedding of the source
paths, for example after representative collapse, centroid consensus,
symbolic-block emission, or any approximation that changes sequence or path
membership.

Resolution levels:

```text
fine
  aggressively smooth small bubbles into local POA graphs

balanced
  smooth tractable bubbles, compress moderate sites, preserve large topology

coarse
  collapse redundant traversals into representatives/polytypes

topology
  preserve the input sparse graph topology with only annotations and cleanup
```

This is the graph-compression analogue of variant calling. The caller chooses
where the graph should represent base-level homology and where it should
represent only coarser equivalence classes of haplotype paths.

## VCF Emission At Chosen Resolution

The VCF writer should consume the resolved hierarchy instead of rediscovering
variants from raw graph topology. It is a downstream exporter from graph
resolution, not the reason the hierarchy exists.

Proposed modes:

```text
--vcf-resolution atomic
  one record per atomic bubble when it can be linearized

--vcf-resolution site
  one record per selected flubble/site; alleles are observed traversals when
  bounded and symbolic alleles otherwise

--vcf-resolution representative
  one record per selected site with representative traversal alleles only

--vcf-resolution block
  one record per complex locus, primarily symbolic, with sidecar traversal
  metadata for downstream genotyping/deconvolution
```

For linearizable sites, emit ordinary REF/ALT sequences. For larger or
ambiguous sites, use symbolic alleles and INFO fields:

```text
IMPG_SITE=site_id
IMPG_LEVEL=2
IMPG_KIND=flubble
IMPG_CHILDREN=child_site_ids
IMPG_REP=representative_ids
IMPG_METHOD=k_medoids
IMPG_MAX_DIST=...
IMPG_PATH_SUPPORT=...
IMPG_CN=...
IMPG_GRAPH=graph_id
```

When a site cannot be represented honestly as one VCF record, the writer should
split it at the highest safe level and record why:

```text
IMPG_SPLIT_REASON=cycle|too_many_alleles|ambiguous_boundary|missing_ref_path
```

VCF emission should be deterministic:

- sort sites by reference path, start, end, level, and site ID
- keep allele order stable: reference, exact high-support alleles, selected
  representatives, then symbolic residual alleles
- record all thresholds and versioned graph IDs in the header

## Graph Simplification And Cleaning

The hierarchy is also a graph simplification plan.

At low resolution, simplification can:

- compact linear chains that are not site boundaries
- remove tips not supported by any source path or sample evidence
- merge exact duplicate traversal sequences inside a site
- annotate low-confidence segments without deleting them

At representative resolution, simplification can:

- replace each selected site with representative traversal paths
- collapse member traversals onto their representative when allowed
- preserve child site annotations as nested metadata
- write a coordinate map from original path intervals to simplified path
  intervals, representative paths, or symbolic site records when exact path
  embeddings are no longer present

This is closely related to assembly graph simplification, but the objective is
different. Assembly cleaning often tries to recover one assembly graph or remove
errors. IMPG should preserve population alleles and source path identity while
making the graph easier to emit, view, and genotype at a requested resolution.

The simplifier should therefore have two classes of operations:

```text
lossless-ish normalization
  compaction, duplicate collapse, path-preserving cleanup

lossy resolution reduction
  representative-only site replacement, rare traversal parking,
  symbolic residual alleles
```

Lossy operations must require explicit resolution settings and must produce a
provenance sidecar. Exact smoothing can avoid a separate coordinate map if it
emits correct full paths through the resolved graph.

## Genotyping And Deconvolution

Representatives should become genotype candidates in the same feature-pack
model used elsewhere in IMPG.

For each representative:

```text
representative traversal -> graph feature vector
sample reads/alignments   -> evidence feature vector
genotype combination      -> sum of representative vectors
score                     -> cosine, count likelihood, alignment likelihood,
                             or learned local scorer
```

This fits the existing architecture:

- syng packs can score representatives by syncmer-node support when the site is
  derived from syng paths
- local GFA/GBZ alignment can score representatives by graph-node or path-step
  support
- bubble-level feature packs can score by child allele support
- read-walk evidence can disambiguate representatives with similar aggregate
  coverage but different ordering or phase

The bottom-up compression hierarchy also gives a natural feature hierarchy for
genotyping. Fine bubbles can contribute exact node or sequence support. Larger
compressed bubbles can contribute support for their replacement paths. Large
passthrough blocks can contribute coarse topology, copy-number, and read-walk
features without pretending that they have a trustworthy base-level alignment.

Deconvolution is the inverse view: after choosing a representative genotype,
estimate which child bubbles or member traversals explain the evidence.

Useful outputs:

```text
top representative genotype
assigned child allele probabilities
copy-number estimate per subsite
residual unexplained evidence
member traversal posterior weights
read-link support between child choices
```

This matters when a VCF is emitted coarsely but downstream analysis still wants
fine structure. A sample may genotype as representative `R3` for a large CNV
site, while deconvolution reports that the evidence prefers child alleles
`b1=A,b2=DEL,b3=copy2:B` inside that representative cluster.

## CLI And Config Knobs

The first user-facing knobs should be few and stable:

```text
impg query -o gfa:syng \
  --graph-resolution balanced \
  --bubble-smoothing local-poa \
  --max-poa-span 1000 \
  --max-poa-traversal-len 3000 \
  --bubble-compression biwfa \
  --emit-resolution-index out.resolution.json
```

Potential options:

```text
--graph-resolution fine|balanced|coarse|topology
--site-operation auto|passthrough|local-poa|pair-compress|representative|symbolic
--bubble-smoothing none|local-poa|biwfa|auto
--bubble-min-anchor-bp N
--bubble-max-span N
--max-poa-span N
--max-poa-traversal-len N
--max-poa-traversals N
--max-biwfa-span N
--max-biwfa-pairs N
--site-merge-distance N
--site-max-span N
--max-site-traversals N
--max-site-alleles N
--representative-method exact|medoid|k-medoids|centroid|k-centroids|polytype
--representative-k N
--representative-max-distance D
--include-rare-long-alleles N
--min-path-support N
--preserve-reference yes|no
--emit-simplified-gfa PATH
--emit-resolution-index PATH
--emit-traversal-fasta PATH
--vcf-resolution atomic|site|representative|block
```

Defaults should favor auditability:

- include the reference traversal
- use observed traversals rather than synthetic centroids for graph paths or
  VCF alleles
- preserve large sparse graph topology when local alignment budgets are
  exceeded
- avoid lossy representative collapse unless allele/traversal limits are hit or
  the user requests representative/block resolution
- write clear warnings for ambiguous sites

## Validation And Adversarial Tests

Start with synthetic graph fixtures where the expected hierarchy is known.

Tests should cover:

- simple SNP/indel bubbles emitted identically at atomic and site resolution
- adjacent independent SNPs that remain split unless merge distance or linkage
  requires a parent site
- small bubbles that are smoothed into a locally consistent POA graph
- small bubbles with poor original alignment that improve after local POA
- moderate bubbles compressed by BiWFA pair induction without losing path
  membership
- large bubbles that deliberately pass through as coarse topology rather than
  forced POA
- insertion with internal SNPs, emitted as child bubbles at atomic resolution
  and one insertion-site record at site resolution
- nested deletion and inversion-like path where boundaries are ambiguous
- tandem duplication with one-copy, two-copy, and three-copy traversals
- paralogous swapped-copy haplotypes that have similar node coverage but
  different ordered read-walk support
- repeat tangle with too many possible walks, where only observed path
  traversals are used
- rare long insertion that must not be collapsed into a common short centroid
- graph with unsupported tips, where cleanup removes or annotates tips only
  under explicit settings
- missing or multiple reference traversals, requiring symbolic output or split
  records
- deterministic output under path order permutations

Validation metrics:

```text
site boundary precision/recall against hand-labeled fixtures
resolved graph node/edge count versus input
sequence/path preservation through emitted graph paths
coordinate-map correctness when lossy rewrites are requested
allele reconstruction identity for emitted sequence alleles
path coverage by representatives
maximum and mean traversal-to-representative distance
VCF record count versus resolution
genotype accuracy under simulated reads
copy-number accuracy under coverage gradients
residual unexplained evidence after deconvolution
runtime and memory on large local graphs
```

Real-data validation should use HPRC-like local regions and known difficult
sites, but the adversarial synthetic set is more important for preventing
regressions in decomposition logic.

## Relation To Other Work Areas

PGGB / smoothxg:

- smoothxg provides the precedent for local POA smoothing and graph lacing
- the proposed difference is block choice: use bubble/flubble topology instead
  of sorted-graph chunk heuristics
- the target is a locally consistent partial order graph up to a configurable
  kernel size, while preserving larger sparse topology above that size
- the useful lesson is not the exact chunking heuristic, but the smoothing
  machinery: extract local sequences, align them into a partial order graph,
  and lace the result back with path provenance intact

Assembly graph simplification:

- shares graph cleanup operations such as tip trimming, chain compaction, and
  duplicate path collapse
- differs because IMPG must retain population alleles, source path membership,
  and selectable resolution rather than producing one cleaned assembly

Variant calling:

- VCF emission becomes a view over graph hierarchy and representative traversal
  choice
- the same graph can emit fine small-variant records or coarse SV/CNV records
  depending on analysis needs

Local haplotype inference:

- representative selection defines candidate local haplotypes
- genotyping and deconvolution score those candidates against sample evidence
- read linkage and panel path continuity can later stitch local calls into
  phased mosaics

Learned or simulation-trained scoring:

- the hierarchy gives natural training labels: site boundaries, child allele
  vectors, representative assignments, and genotype states
- simulated reads over known traversals can train scoring functions to choose
  the correct representative genotype and deconvolve child structure
- learned scorers should remain optional and auditable, with deterministic
  heuristic scoring as the baseline

## Open Risks

- Graph decomposition may be unstable in cyclic or highly repetitive regions.
  The model needs explicit boundary confidence and safe fallback to block-level
  symbolic output.
- Representative collapse can hide rare clinically or functionally important
  alleles. Long, high-impact, or high-distance outliers should force their own
  representative by default.
- VCF cannot faithfully encode all graph traversal structure. Sidecar metadata
  is not optional for block/polytype mode.
- Synthetic centroids may not correspond to real haplotypes. Prefer medoids for
  emitted VCF alleles unless the user explicitly requests centroid consensus.
- Coverage-only scoring can confuse paralogous copies with similar node
  content. Read-walk order, edge support, and local realignment should be used
  where available.
- Source path and coordinate provenance can become large. The sidecar format
  needs compression and stable IDs.
- Resolution choices may affect downstream association or annotation. VCF
  headers and sidecars must record all thresholds and assignment methods.

## Suggested Implementation Phases

1. Add an internal resolution index for local GFAs with bubble/flubble
   detection, cycle/tangle marking, path-supported traversal extraction, and
   JSON debug output. Start with POVU output where possible and keep the data
   model independent enough to replace pieces later.
2. Add deterministic site clustering over adjacent/nested bubbles, tight
   tangles, and cyclic local components with stable site IDs.
3. Add the deterministic alignability classifier and operation assignment:
   passthrough, local POA, BiWFA pair-compression, representative, or symbolic
   block.
4. Add graph output passthrough and exact path-preserving rendering, proving
   that the hierarchy can round-trip a graph without changing path sequences or
   path coordinates.
5. Add bubble-guided local POA smoothing for small clean bubbles and lace the
   smoothed subgraph back through entry/exit handles.
6. Add the bottom-up iterative compression loop, rerunning decomposition after
   replacements and stopping when no eligible bounded sites remain.
7. Add BiWFA pair-compression for moderate bubbles using bounded alignment
   schedules such as guide trees, spanning trees, nearest neighbors, and random
   chords instead of default all-vs-all alignment.
8. Add exact and medoid representative selection with feature vectors for
   larger or repetitive sites.
9. Add coordinate/provenance sidecars for representative, symbolic, centroid,
   or other lossy modes.
10. Route `query -o vcf` through the hierarchy for `atomic` and `site` modes,
   keeping POVU as the sequence-allele backend where possible.
11. Add representative/block VCF modes with symbolic alleles and sidecar
   traversal metadata.
12. Connect representatives to `genotype`/`infer` as candidate haplotypes in a
   declared feature space.
13. Add simulation-trained or learned scoring only after deterministic fixtures
    and real-region validation are stable.

The key design constraint is that resolution is not a VCF formatting option.
It is a shared graph model used by bubble-guided smoothing, graph compression,
graph simplification, VCF emission, and sample evidence interpretation.
