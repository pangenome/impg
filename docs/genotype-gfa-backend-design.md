# Backend-Neutral GFA Genotype Evidence Design

This document designs the path from the current `syng-syncmer-node` genotype
evidence path to backend-neutral GFA and variation-graph genotyping. It is
scoped as a design artifact: no current syng semantics should change until the
GFA backend is added behind typed feature-space checks.

The current audit in
[`genotype-evidence-audit.md`](genotype-evidence-audit.md) is the baseline:
`impg genotype cos` consumes a syng index or syng-native render bundle plus a
syng pack/projection, scores unsigned syncmer-node count vectors, and rejects
local render bundles whose feature space is `gfa-segment`. The new path must
extend the feature-vector model without making syng packs look like GFA packs.

## Decision Summary

- The default GFA feature namespace should be **render-bundle translation
  feature IDs** with `feature_space = gfa-segment`. Each internal dense
  `feature_id` maps to one original GFA segment name, sequence length, and
  optional source-coordinate/path-step metadata. Raw GFA segment IDs are
  user-facing aliases, not the machine identity.
- Path-step-aware translation is required for evidence projection and source
  compatibility, but path steps are not the first scoring namespace. If segment
  aggregation is later insufficient for repeat-heavy loci, add an explicit
  `gfa-path-step` feature space rather than overloading `gfa-segment`.
- Bubble alleles are a later derived feature space. They should be built from
  the same graph/translation tables after node-level scoring is inspectable.
- Read evidence should enter GFA genotyping through a typed projection object
  that contains a pack-like vector, read-contribution debug tables, graph
  identity, feature-space identity, and projection parameters.
- Variable-length graph nodes must be length-normalized for scoring. For GFA
  features, cosine should use depth-like weights:
  `sample_weight_f = observed_aligned_bp_f / scoring_length_f` and
  `candidate_weight_f = candidate_traversed_bp_f / scoring_length_f`. The
  current syng score remains unchanged.
- The minimal CLI should start with a projection boundary:
  `impg project --gfa/--render-bundle ... --gaf/--syng-pack/--bam ... -O
  sample.proj`, then `impg genotype cos --render-bundle ... --proj sample.proj`.
  `impg map --gfa` can become a convenience wrapper once graph-read mapping is
  implemented. Direct `genotype cos --gfa graph.gfa --pack sample.pack` should
  wait until packs are typed enough to reject syng/GFA mismatches reliably.

## Implementation Status

`impg genotype cos` now has an initial backend-neutral graph path for
`gfa-segment` / `variation-graph-node` evidence while preserving the existing
syng syncmer-node path.

Implemented graph source modes:

- `--graph <local.gfa>` loads a prebuilt local GFA directly.
- `--render-bundle <bundle>` loads render bundles whose manifest declares
  `feature_space = gfa-segment` or `variation-graph-node` and contains
  `graph_gfa`.
- `-a <panel.syng> -r <REF:start-end> --sequence-files ... --gfa-engine ...`
  builds a dynamic GFA through the same `dispatch_gfa_engine` path used by
  `impg query -o gfa`, then genotypes that emitted graph.
- A GFA produced by `impg query -o gfa...` can be handed back through
  `--graph` with identical downstream semantics to an independently built GFA.

Implemented graph feature behavior:

- GFA `S` records define the node universe. Length comes from sequence length
  or `LN:i:` when the sequence is `*`.
- GFA `P` paths and `W` walks become candidate haplotypes. W-line names are
  imported as Pan-SN-like `sample#hap#seqid` names.
- `--target-path PATH[:start-end]` selects a path interval. Without an
  interval, each graph path is scored as a full-path candidate. For interval
  targets, the same path-coordinate interval is clipped across candidate paths.
- `--graph-feature-id-mode auto|dense|segment-name` controls integer feature
  IDs. `auto` uses numeric segment names when all segment names are unique
  positive `u32`s, otherwise dense import order.
- `--graph-contribution-model raw` scores raw node traversal counts.
  `--graph-contribution-model length-normalized` divides sample counts by node
  length and scores partial candidate traversals as covered-bp divided by node
  length.

Pack compatibility is explicit. A graph pack must declare `feature_space` in
one of these ways:

```text
#feature_space    gfa-segment
#feature_id_mode  segment-name
```

or a sidecar:

```text
sample.pack.tsv.meta.tsv
sample.pack.tsv.metadata.tsv
```

with `key<TAB>value` rows, or the CLI override:

```bash
--pack-feature-space gfa-segment
```

If a pack declares `graph_id`, it must match the loaded graph ID. If it
declares `feature_id_mode`, it must match the effective GFA import mode. Pack
feature IDs not present in the loaded graph are rejected.

The debug report from `--emit-report` now has graph sections for source,
dynamic build command, graph ID, feature universe, node lengths, sample raw
counts and weights, candidate counts, contribution model, and node-level
dot/norm decomposition.

Remaining limitation: this implementation consumes a pack-like vector over
graph feature IDs. It does not map BAM, CRAM, FASTQ, or GAF reads into graph
packs. That projection layer should be implemented as a separate typed
converter so read-walk evidence, graph identity, feature namespace, and
normalization policy are preserved rather than faked inside genotyping.

## Current Constraints To Preserve

The syng backend has deliberately simple semantics:

- Feature IDs are unsigned syng syncmer nodes.
- `impg map -o pack` deduplicates node IDs per read before incrementing counts.
- Candidate vectors count syncmer-node occurrences in candidate path intervals.
- Cosine uses raw counts with no node-length, read-length, mapping-quality, or
  orientation component.
- Projection bundles may carry GAF read walks, but `genotype cos` consumes only
  the pack path.

Those semantics are documented behavior. The GFA backend should use shared
scoring machinery where possible, but it must be a separate feature-space
adapter with explicit metadata and error messages. A current syng `.pack` file
must not be accepted as `gfa-segment` evidence just because both vectors use
integer-looking feature IDs.

## Feature Namespace

### Recommended Default: Rendered Graph Feature IDs

Use the render bundle or graph import step to assign a dense stable feature ID:

```text
feature_id: u32
feature_space: gfa-segment
graph_id: digest of normalized graph topology, segment sequences, paths, and
          translation-table identity

feature table row:
  feature_id
  original_segment_name
  segment_sequence_length
  scoring_length
  graph_component / local graph id
  optional source intervals observed on rendered paths
```

This matches the existing render-bundle direction in
[`render-gbz-translation-design.md`](render-gbz-translation-design.md), where a
rendered graph is a view over the source sequence namespace and the feature
space is declared in the manifest. The current local graph bundle already uses
`feature_space = gfa-segment`; the missing piece is making those segment
features consumable by genotype.

Reasons not to use raw segment IDs as the machine namespace:

- GFA segment names are arbitrary strings and may be renumbered by sorting,
  chopping, unchopping, crush, or external graph tools.
- Dense integer IDs are needed for compact vector storage and random access.
- The same external segment name could occur in unrelated imported graphs.
- A typed `graph_id + feature_space_id + feature_id` key makes pack/graph
  compatibility checkable.

Reasons not to use path-step translated syng nodes as the GFA namespace:

- They would keep the scorer in syng space and fail to genotype graphs whose
  useful evidence comes from local graph alignment.
- A syng node occurrence must be path-step aware to avoid repeat mistakes, but
  that is projection metadata, not the local graph feature identity.
- For local POA/seqwish/pggb/crush graphs, several syng syncmers may overlap one
  graph segment, and one syncmer may span multiple graph segments.

Reasons not to start with bubble alleles:

- Bubble decomposition is graph-transform dependent.
- Nested and overlapping bubbles need a hierarchy before they are safe as a
  scoring namespace.
- Node-level vectors are easier to debug against GAF paths, rendered paths, and
  the existing validation suite.

### Optional Future Namespaces

The first implementation should make feature-space selection explicit so later
namespaces are additive:

```text
gfa-segment       default node/segment aggregate support
gfa-path-step     occurrence-aware support keyed by rendered_path_id + step
bubble-allele     derived allele support over a validated bubble hierarchy
haplotype-segment  larger copied blocks or representatives
```

The scorer should require candidate vectors and sample vectors to have exactly
the same `graph_id` and `feature_space_id`.

## Projection Object And Pack-Like Vector

The existing `sample.proj/` idea should generalize from syng to typed evidence:

```text
sample.gfa.proj/
  manifest.json
  features.tsv.zst              optional debug feature table
  sample.vector.tsv.zst         normalized vector debug export
  sample.pack                   compact vector, once a typed binary exists
  reads.gaf.zst                 optional original or normalized read walks
  read-contrib.tsv.zst          optional per-read feature contributions
  projection-summary.tsv        optional filters and aggregate counters
```

The manifest should declare at least:

```json
{
  "format": "impg-projection",
  "version": 2,
  "graph_id": "...",
  "feature_space": "gfa-segment",
  "feature_space_id": "...",
  "projection_method": "gaf-graph-walk",
  "sample": "...",
  "read_source": "...",
  "weighting": "aligned-bp-per-scoring-length",
  "graph": "graph.gfa",
  "render_bundle": "locus.impg-gbz",
  "pack": "sample.pack",
  "read_contrib": "read-contrib.tsv.zst"
}
```

The in-memory type should be graph-backend neutral:

```text
TypedCoverage
  graph_id
  feature_space_id
  feature_space_kind
  counts/weights: feature_id -> FeatureObservation

FeatureObservation
  reads: u64
  aligned_bp: u64
  mapq_weighted_bp: f64
  weight: f64
```

The current syng `pack::Coverage` can be adapted into this type with
`feature_space_kind = syng-syncmer-node` and integer weights, but the on-disk
syng pack format should not be changed in the GFA slice.

## Evidence Input Paths

All evidence inputs should produce the same two layers:

```text
read/alignment records
  -> per-read feature contribution rows
  -> aggregate TypedCoverage vector
```

The per-read layer is essential because the audit found that current syng
debugging requires manual correlation between GAF paths, pack counts, and
genotype output. GFA should not repeat that gap.

### GAF To GFA Features

The first implemented GFA projection supports GAF records whose path field is a
graph walk over segment names:

```text
read1  qlen  qstart  qend  +  >s1>s2<s3  tlen  tstart  tend  ...
```

Implemented projection steps:

1. Resolve each segment name in the GAF path to `feature_id`.
2. Use segment lengths and GAF target start/end to compute per-step graph-bp
   overlap. The current converter uses graph-span overlap and raw traversal
   count deltas; `cg` / `cs` base-level refinement is still future work.
3. Emit one contribution row per retained read-feature overlap:

```text
read_name  read_ordinal  step_index  segment_name  orientation
feature_id  segment_visit_in_read  count_delta  explanation
```

4. Aggregate by feature into a typed pack TSV:

```text
#feature_space              gfa-segment
#graph_id                   gfa-fnv1a64:...
#feature_id_mode            segment-name|dense
#graph_contribution_model   raw|length-normalized
feature_id                  count
```

Unlike syng pack construction, GFA projection should not deduplicate repeated
visits to the same feature within a read by default. If a graph alignment walks
the same segment twice, both traversals carry copy-number evidence. The
read-contribution table makes repeated visits inspectable with
`segment_visit_in_read` and an explanation field, for example `repeated visit 2
to segment in read; counted again`.

Future projection refinements can add primary/supplementary filtering,
minimum mapQ, multimapping weight splitting, and base-level aligned-bp weights
without changing the typed pack compatibility checks.

### Syng Pack To Local GFA Features

For render bundles built from syng candidates, syng evidence can be projected
without remapping reads using the path-step-aware records designed in
[`syng-to-local-graph-translation.md`](syng-to-local-graph-translation.md):

```text
syng_node_id + rendered_path_id + source interval
  -> local gfa feature_id + orientation + source overlap
```

Projection should use overlap-weighted contributions:

```text
for each syng_to_graph projection record:
  syng_count = syng_pack[syng_node_id]
  contribution_bp = syng_count * overlap_bp
  gfa_observed_bp[local_feature_id] += contribution_bp

sample_weight_f = gfa_observed_bp_f / scoring_length_f
```

This is a bridge for existing syng packs, not a change to syng pack meaning.
The projection manifest must say `projection_method = syng-pack-to-gfa` and
must record the source syng graph identity and target GFA graph identity.

### BAM/CRAM To GFA Features

BAM/CRAM support should be translation-backed rather than raw coordinate
guessing:

- If the BAM reference names match source sequence names in a render bundle,
  project aligned reference intervals through the bundle's source-to-rendered
  path-step translation.
- If the BAM is aligned to rendered haplotype FASTA paths, resolve those path
  names through the rendered path table.
- If the BAM is aligned only to one linear reference path, it can support
  features on that path but cannot directly support alternate-only segments
  without local realignment.

The projection output is still the same per-read contribution table and
aggregate vector. BAM projection should record CIGAR-derived aligned bp,
mapping quality, and optional base-quality/damage summaries for later scorers.

### FASTA/FASTQ To GFA Features

FASTA/FASTQ require a graph-read mapper. Do not make `genotype cos` map reads
internally. The intended flow is:

```text
reads.fa/fq
  -> graph mapper producing GAF/read walks
  -> impg project --gfa/--render-bundle --gaf reads.gaf
  -> sample.gfa.proj
  -> impg genotype cos --proj sample.gfa.proj
```

`impg map --gfa graph.gfa -q reads.fq -o proj -O sample.gfa.proj` can be added
later as a convenience wrapper around mapper selection plus the same projector.
Until mapper semantics are chosen, accepting externally produced GAF is the
smallest inspectable slice.

## Node Lengths And Scoring

The syng scorer can ignore node length because syncmer features are fixed-size
anchors and the implemented semantics are raw counts. GFA nodes can range from
one base to many kilobases, so raw node counts would make a 1 bp SNP segment and
a 10 kb invariant segment comparable in a way that is biologically wrong.

GFA cosine should use length-normalized entries:

```text
scoring_length_f = segment sequence length, or a validated mappability-adjusted
                   length in a later implementation

sample_weight_f = mapq_weighted_aligned_bp_f / scoring_length_f

candidate_weight_f = sum(candidate path-step bp overlap with feature f)
                     / scoring_length_f

genotype_weight_f = candidate_weight_f(h1) + candidate_weight_f(h2) + ...

score = cosine(sample_weight, genotype_weight)
```

Interpretation:

- A haplotype traversing a full segment once contributes approximately `1.0`.
- A diploid genotype traversing it twice contributes approximately `2.0`.
- A long segment needs proportionally more aligned bases to reach the same
  depth-like weight.
- Partial boundary nodes contribute by the number of candidate-overlapped bases.
- Repeated traversal of the same segment by one candidate contributes more than
  `1.0`, preserving copy-number evidence.

Implementation requirements:

- Reject or quarantine GFA segments with unknown length (`S` sequence `*`
  without an `LN` tag or equivalent bundle length table).
- Store both raw `aligned_bp` and normalized `weight` so debug output can
  explain length effects.
- Keep `scoring_length` in the feature table. Later mappability or mask-aware
  correction can change `scoring_length` without changing feature identity.
- Keep the syng backend on its existing raw count vectors. Any length-normalized
  syng experiment should be a separate feature space or scoring mode.

## Candidate Vectors On GFA

The backend-neutral genotype layer needs the same operations described in
[`pangenome-genotyping-roadmap.md`](pangenome-genotyping-roadmap.md):

```text
GraphBackend
  graph_id()
  feature_space()
  candidates(target range or graph subregion)
  candidate_features(candidate) -> feature vector
  candidate_walk(candidate) -> ordered oriented feature walk
  fetch_sequence(candidate) -> DNA sequence, when available
```

For a render bundle with local GFA paths:

1. Translate the requested source `-r seq:start-end` through the namespace and
   rendered path tables into graph path-step intervals.
2. Collect candidate haplotype subwalks from rendered paths that overlap or
   span the target interval, preserving Pan-SN identity.
3. Build candidate feature vectors by sweeping the candidate path steps and
   adding `step_overlap_bp / scoring_length_f` to each feature.
4. Keep an oriented candidate walk for read-link/infer extensions, even if the
   first GFA cosine score is unoriented.

For an imported external GFA without render translation:

- Candidates can be graph paths or named walks in graph coordinates.
- `-r` must refer to a graph path coordinate unless a separate source
  translation table is supplied.
- Output can report graph path names, but source Pan-SN fields are available
  only when path names parse as Pan-SN or the import provided metadata.

## Render Bundle And Pan-SN Compatibility

GFA genotyping must preserve the same source namespace invariant as render and
infer:

```text
source sequence interval
  <-> rendered path step
  <-> GFA feature_id / segment name
  <-> evidence vector feature
  <-> genotype output haplotype and region
```

Compatibility rules:

- Use render-bundle `SequenceNamespace` and `RenderedPathRecord` as the
  authority for Pan-SN sample, haplotype, contig, full path name, and source
  interval.
- Do not infer biological sample/haplotype fields from arbitrary non-Pan-SN
  names. Store them as unparsed source names.
- `genotype` output for GFA should report both source intervals and graph
  feature-space metadata:

```text
#feature_space  gfa-segment
#graph_id       ...
#projection     sample.gfa.proj
#weighting      aligned-bp-per-scoring-length
```

- Render-bundle target projection should remain the preferred way to genotype
  syng-derived local graphs. It gives `genotype` enough information to translate
  source ranges, preserve Pan-SN naming, and validate evidence compatibility.
- Syng-to-local projection records must be path-step aware. A direct
  `syng_node_id -> gfa_segment` map is not valid because repeats and graph
  smoothing make the same node/segment occur in multiple source contexts.

## CLI Shape

### First Slice

The first slice now includes a projection command for external GFA/GAF evidence:

```bash
# External graph aligner produced GAF over graph.gfa segment names.
impg project \
  --gfa graph.gfa \
  --gaf reads.gaf \
  -O sample.gfa.proj

# Syng evidence projected through a local render bundle.
impg project \
  --render-bundle locus.impg-gbz \
  --syng-pack sample.syng.pack \
  -O sample.gfa.proj

# Genotype a typed local graph bundle.
impg genotype cos \
  --render-bundle locus.impg-gbz \
  --proj sample.gfa.proj \
  -r HG00000#0#chr6:31972057-32055418
```

For a pure external GFA with paths but no source translation:

```bash
impg genotype cos \
  --gfa graph.gfa \
  --proj sample.gfa.proj \
  -r graph_path:10000-20000
```

The `--proj` path resolves the bundle's typed pack and relies on pack metadata
for feature-space, graph ID, feature ID mode, and contribution-model checks.
Direct `--gfa --pack` also works for typed pack TSVs, but `--proj` preserves the
read-contribution table and manifest alongside the aggregate counts.

### Later Convenience

After graph-read mapping is implemented:

```bash
impg map \
  --gfa graph.gfa \
  -q reads.fq.gz \
  -o proj \
  -O sample.gfa.proj
```

This command should produce the same projection format as `impg project`. It is
a convenience entry point, not a separate evidence model.

### What To Avoid Initially

Avoid this as the first interface:

```bash
impg genotype cos --gfa graph.gfa --pack sample.pack
```

The existing `--pack` flag is strongly associated with syng packs, and current
binary packs have no graph identity or feature-space metadata. A direct
`--gfa --pack` mode is acceptable only after typed packs/projections can reject:

- syng pack against GFA graph;
- GFA pack against a different graph digest;
- same graph with a different feature namespace;
- same graph after segment renumbering without a matching translation table.

## Debuggability Requirements

The audit lists several missing human-debug artifacts for the current syng path.
The GFA backend should add them as part of the first implementation, then later
syng can reuse the same debug surface.

Required debug exports:

```text
features.tsv
  feature_id  segment_name  length  scoring_length  source_context_count

read-contrib.tsv
  read_name  read_ordinal  step_index  segment_name  orientation
  feature_id  segment_visit_in_read  count_delta  explanation

candidate-features.tsv
  candidate_id  path_name  source_interval  feature_id  segment_name
  traversed_bp  normalized_weight

score-residuals.tsv
  rank  genotype  feature_id  sample_weight  genotype_weight  residual

candidate-filters.tsv
  candidate_id  path_name  reason  anchors/features/span/top_k_status

projection-summary.tsv
  reads_seen  alignments_seen  alignments_retained  reads_retained
  aligned_bp  filtered_by_reason...
```

These tables should be optional (`--debug-evidence-dir` or manifest debug
exports) but covered by tests. They directly answer the audit gaps: per-read
pack contribution, candidate vectors beside sample counts, node-length effects,
feature residuals, and candidate filter reasons.

## Staged Implementation Plan

### Stage 1: Typed Feature Tables And Coverage Types

- Add a graph-backend-neutral `FeatureSpace`/`TypedCoverage` model in memory.
- Add a GFA feature table builder that assigns dense `feature_id`s, records
  original segment names, validates lengths, and computes a `graph_id`.
- Add compatibility checks: `graph_id`, `feature_space`, and
  `feature_space_id` must match before scoring.
- Keep current syng commands and on-disk syng pack format unchanged.

Validation:

- Unit tests for stable GFA feature ID assignment and graph digest changes when
  segment sequence/path topology changes.
- Unit tests rejecting unknown segment lengths.
- Unit tests rejecting feature-space/graph mismatches.

### Stage 2: Projection Builders

- Implement GAF-to-GFA projection into `TypedCoverage`.
- Implement syng-pack-to-local-GFA projection using path-step-aware
  `syng_to_graph` records from render bundles.
- Emit projection manifests plus `read-contrib.tsv` and `features.tsv` debug
  exports.
- Keep BAM/CRAM as design-supported but implement after GAF/syng projection
  unless a downstream task needs it first.

Validation:

- Synthetic SNP bubble: reads over each allele produce allele-specific segment
  weights.
- Insertion/deletion bubble: length-normalized weights call the correct allele
  instead of favoring the longer node by raw count.
- Repeated segment/path case: a read traversing the same segment twice records
  two contributions and aggregates expected copy support.
- Reverse-orientation GAF path: feature IDs aggregate correctly and orientation
  is preserved in the read-contribution table.
- Syng projection boundary case: a syncmer spanning two graph nodes splits by
  overlap.

### Stage 3: GFA Candidate Backend For `genotype cos`

- Factor the current syng-specific compute path so the scorer can consume
  backend-provided candidates and typed vectors.
- Implement candidate extraction from render-bundle local GFA paths first.
- Build candidate vectors with `traversed_bp / scoring_length`.
- Add GFA output headers and optional `candidate-features.tsv` and
  `score-residuals.tsv`.

Validation:

- Existing syng genotype and infer tests continue to pass unchanged.
- `genotype cos --render-bundle local-gfa --proj sample.gfa.proj` calls a known
  synthetic genotype.
- Mismatched projection/graph gives an explicit error before scoring.
- `candidate-top-k`, ploidy, top-n, and combination budget behavior match the
  existing scorer on a small deterministic backend-neutral fixture.

### Stage 4: CLI Integration

- Add `impg project --gfa graph.gfa --gaf reads.gaf -O sample.gfa.proj`.
- Add `impg project --render-bundle bundle --syng-pack sample.pack -O
  sample.gfa.proj`.
- Add `genotype cos --render-bundle <gfa-segment bundle> --proj <typed proj>`.
- Add `genotype cos --gfa <graph.gfa> --proj <typed proj>` only when graph path
  coordinate candidate extraction is ready.
- Keep direct `--gfa --pack` out of the first release unless the pack is typed.

Validation:

- CLI tests for help text and expected errors.
- CLI tests proving a syng pack cannot be used against `--gfa`.
- CLI tests proving a GFA projection cannot be used against a syng index.
- CLI tests for compressed debug/projection outputs where supported.

### Stage 5: Mapping And BAM Extensions

- Add BAM/CRAM projection through render-bundle source/path translation.
- Add `impg map --gfa` as a wrapper around a chosen graph alignment backend plus
  the same projection builder.
- Add read-link support for GFA candidate walks in `infer`.

Validation:

- BAM path-coordinate fixture with CIGAR clipping, insertion, deletion, and
  reverse-strand alignments.
- FASTQ-to-GFA smoke test once mapper selection is deterministic enough for CI.
- GFA read-link/infer test where read order distinguishes two equal cosine
  candidates.

## Test Strategy

The test suite should extend, not replace, the current syng suite documented in
the audit.

### Must Not Regress

- All current `tests/test_syng_integration.rs` genotype, map, projection, and
  infer tests should continue to pass.
- Current pack TSV/binary output remains byte/field compatible.
- Current syng render-bundle genotype/infer support remains limited to
  `syng-syncmer-node` until the GFA branch is explicitly implemented.

### New Unit Tests

- Feature table import from GFA `S`, `P`, and `W` records.
- Pan-SN path parsing through imported graph paths and render bundles.
- Length-normalized sample/candidate vector calculations.
- Graph/feature-space compatibility errors.
- Per-read contribution aggregation, including repeated feature visits.

### New Integration Tests

- Tiny two-haplotype SNP graph with GAF evidence for homozygous and
  heterozygous calls.
- Insertion/deletion graph where raw node counts would be misleading but
  length-normalized weights recover the expected call.
- Repeated-copy graph matching the audit's repeated-node concern, with
  `read-contrib.tsv` proving how repeated visits are counted.
- Render-bundle local graph test preserving Pan-SN source names in genotype
  output.
- Syng-pack projection to local GFA with a known synthetic bundle and
  overlap-split syncmer.
- Negative CLI tests for pack/projection feature-space mismatch.

### Validation Matrix

Once the first GFA backend works, add a small known-truth simulation matrix:

```text
coverage:        low, medium, high
read length:     short, long
error rate:      exact, noisy
ploidy:          1, 2, 3
copy number:     single, duplicated, triplicated
allele type:     SNP, insertion, deletion, inversion-like orientation case
repeat content:  none, shared flank, repeated internal segment
decoys:          absent, close unsampled decoy
```

The matrix should emit both calls and debug artifacts. The assertion should not
only be "top genotype is correct"; it should also check that:

- sample and candidate feature spaces match;
- nonzero sample features are explainable from `read-contrib.tsv`;
- candidate vector entries are explainable from `candidate-features.tsv`;
- length-normalized weights differ from raw aligned bp where node lengths vary;
- filter decisions explain missing candidates.

### Real-Data Smoke Path

After synthetic coverage is stable, add a C4/HPRCv2-style local render smoke
test when the dataset path is stable. This should use the same output contract:
projection manifest, typed feature pack, candidate features, score residuals,
and Pan-SN-preserving genotype rows.

## Open Questions

- Should the compact typed pack store `f64` weights directly, or store raw
  integer bp/count fields plus a manifest-defined normalization transform? The
  safer first implementation is raw fields plus a debug TSV; compact storage
  can follow once the normalization policy settles.
- Should `gfa-path-step` be implemented before BAM projection? Segment-level
  support is simpler, but BAM/source projection may expose repeats where
  path-step features are more faithful.
- How much GAF CIGAR parsing is required for the first slice? Coarse path-span
  overlap is inspectable and useful, but CIGAR-aware aligned bp will matter for
  noisy and indel-heavy reads.
- Should invariant flank nodes be downweighted or masked for cosine? The first
  backend should match current "all locus features" behavior, then add an
  explicit informative-feature mask only after debug residuals show the need.

## Bottom Line

Treat GFA genotyping as a typed graph-feature projection path:

```text
GFA/render bundle
  -> stable feature table and graph identity
reads/alignments/syng pack
  -> per-read or per-source projection records
  -> length-normalized typed coverage vector
candidate graph paths
  -> length-normalized typed candidate vectors
genotype cos
  -> same combination scorer, guarded by feature-space compatibility
```

This gets IMPG to backend-neutral COSIGT-style graph genotyping while keeping
current syng packs and syncmer-node scoring intact.
