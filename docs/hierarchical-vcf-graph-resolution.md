# Hierarchical VCF And Graph Resolution

IMPG can already build local sequence graphs from pangenome intervals, including
syng-derived graphs, and route `query -o vcf` through POVU-style GFA-to-VCF
conversion. The missing control surface is resolution: the same local graph may
contain single-base substitutions, nested indels, multi-kilobase SVs, CNV
structure, and repetitive tangles. A flat GFA-to-VCF pass tends to either emit
too many small records or lose the structure needed to type larger sites.

This note proposes a hierarchical bubble decomposition layer between local GFA
construction and VCF/genotyping. The layer should identify fine bubbles,
organize them into larger flubbles/sites, select representative traversals at
configurable resolution, and preserve enough provenance to simplify graphs,
emit VCF, and match sample evidence back to representative haplotypes.

## Problem

Current GFA-to-VCF conversion needs controllable clustering of variant graph
structure, especially in large SV, CNV, and repetitive regions.

The core issues are:

- Local pangenome graphs are often nested: small variants sit inside insertion
  alleles, repeat copies, or alternative structural haplotypes.
- A "one record per smallest bubble" VCF is too fragmented for SV/CNV sites and
  loses haplotype context across nearby dependent choices.
- A "one record per graph component" VCF can become untypeable if it enumerates
  every observed path in a repetitive or copy-variable locus.
- Graph simplification and VCF emission need the same notion of resolution, or
  the simplified graph and the VCF records will disagree.
- Genotyping needs representative haplotypes that can explain read coverage and
  linkage, not just syntactically valid VCF alleles.

The desired behavior is a tunable path:

```text
local GFA
  -> normalize/clean graph
  -> detect atomic bubbles
  -> cluster into hierarchical flubbles/sites
  -> select representative traversals at requested resolution
  -> emit VCF and optional simplified GFA
  -> genotype/deconvolve sample evidence against representatives
```

## Goals

- Support multiple resolution levels from atomic bubbles to large SV/CNV sites.
- Keep source path identity and path intervals through every transformation.
- Select representative traversals that cover common haplotype space while
  bounding allele count.
- Emit VCF records at the requested resolution with clear provenance back to
  graph sites and lower-level bubbles.
- Simplify or clean GFA graphs with the same site hierarchy used for VCF.
- Expose representative haplotypes as feature vectors for genotyping,
  deconvolution, and local haplotype inference.

Non-goals for the first implementation:

- Perfect decomposition of every cyclic assembly graph.
- A new whole-genome variant caller independent of IMPG's local graph and pack
  abstractions.
- Replacing POVU immediately. The first target can wrap or preprocess GFA for
  POVU and then progressively move more logic into IMPG.

## Hierarchical Model

Use a hierarchy of typed graph regions:

```text
level 0: atomic graph features
  segment, oriented segment, edge, path step, syncmer node

level 1: atomic bubbles
  single-entry/single-exit alternative traversals, usually small variants

level 2: compound flubbles / variant sites
  adjacent, nested, or weakly interacting bubbles that should often be emitted
  together, such as an insertion with internal SNPs or a deletion plus nearby
  breakpoint polymorphism

level 3: complex loci / blocks
  large SV, CNV, duplicated, or repetitive regions where representative
  traversal selection is required
```

"Bubble" should be reserved for clean source/sink subgraphs when possible.
"Flubble" is useful for an almost-bubble: a local component with recognizable
entry/exit anchors and path traversals, but with internal nesting, weak cycles,
copy-number variation, or multiple plausible boundaries. VCF emission may use
either, but the model should record boundary confidence and decomposition type.

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
- VCF alleles, simplified GFA paths, and genotype candidates must all reference
  the same site and representative IDs.
- The reference allele is explicit, even when it is not the most common
  traversal.
- Source path identity is never discarded. It may be summarized, but the full
  path membership must remain recoverable for audit and training data.

## Decomposition Algorithm

The first implementation should be deterministic and conservative.

1. Parse and normalize the GFA.

   Build a handle graph view with segment lengths, edges, path walks, and path
   step positions. Optionally compact linear chains, remove exact duplicate
   segments, and annotate tips or low-support graph pieces, but keep a
   provenance map from old handles to new handles.

2. Choose anchor paths and coordinate frames.

   Prefer the requested reference path/range when `query` has one. Also keep
   all source paths from syng/Pan-SN metadata. Site coordinates should be
   reported in reference coordinates when available and graph-local coordinates
   otherwise.

3. Detect atomic bubbles.

   Start with superbubble-like detection on the directed bidirected graph after
   orienting around the reference path. Boundaries should be graph handles with
   multiple path-consistent traversals between them. Reject or downgrade sites
   with ambiguous entry/exit, too many internal cycles, or insufficient path
   support.

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

6. Assign resolution levels.

   Each site can be emitted at multiple levels:

   ```text
   atomic       emit child bubbles independently
   site         emit a compound flubble as one VCF record
   representative emit only selected traversal representatives
   block        emit a larger complex locus with symbolic alleles or sidecar
   ```

   The default should be conservative: atomic for simple sites, site-level for
   obvious nested SV/indel structures, and representative mode only when allele
   explosion would otherwise occur.

7. Select representatives.

   Collapse duplicate traversals first, then select exact or approximate
   representatives according to the requested resolution and bounds.

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

For VCF, exact and medoid representatives can usually become sequence alleles
or symbolic alleles. Centroids and polytypes may require symbolic alleles plus
sidecar annotations because standard VCF cannot fully encode arbitrary graph
traversal structure.

## VCF Emission At Chosen Resolution

The VCF writer should consume the hierarchy instead of rediscovering variants
from raw graph topology.

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
- write a translation table from original handles/path steps to simplified
  handles/path steps

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
provenance sidecar.

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
impg query -o vcf:syng \
  --vcf-resolution site \
  --site-merge-distance 1000 \
  --max-site-alleles 16 \
  --representative-method k-medoids \
  --representative-k 8 \
  --emit-resolution-index out.resolution.json
```

Potential options:

```text
--vcf-resolution atomic|site|representative|block
--bubble-min-anchor-bp N
--bubble-max-span N
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
```

Defaults should favor auditability:

- include the reference traversal
- use observed traversals rather than synthetic centroids for VCF alleles
- avoid lossy representative collapse unless allele/traversal limits are hit or
  the user requests representative/block resolution
- write clear warnings for ambiguous sites

## Validation And Adversarial Tests

Start with synthetic graph fixtures where the expected hierarchy is known.

Tests should cover:

- simple SNP/indel bubbles emitted identically at atomic and site resolution
- adjacent independent SNPs that remain split unless merge distance or linkage
  requires a parent site
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

1. Add an internal resolution index for local GFAs with atomic bubble detection,
   path-supported traversal extraction, and JSON debug output.
2. Add site clustering over adjacent/nested bubbles and deterministic site IDs.
3. Add exact and medoid representative selection with feature vectors.
4. Route `query -o vcf` through the hierarchy for `atomic` and `site` modes,
   keeping POVU as the sequence-allele backend where possible.
5. Add representative/block VCF modes with symbolic alleles and sidecar
   traversal metadata.
6. Add simplified GFA emission using the same hierarchy and provenance maps.
7. Connect representatives to `genotype`/`infer` as candidate haplotypes in a
   declared feature space.
8. Add simulation-trained or learned scoring only after deterministic fixtures
   and real-region validation are stable.

The key design constraint is that resolution is not just a VCF formatting
option. It is a shared graph model used by VCF emission, graph simplification,
and sample evidence interpretation.
