# Pangenome Genotyping Roadmap

IMPG should become a flexible genotyping and genome-inference system on top of
pangenomes. The reason to build it here is that the hard problem is not a single
local scorer. The hard problem is making the whole workflow implicit, scalable,
auditable, and interchangeable across graph representations.

The long-term target is:

```text
panel sequences / pangenome graph
  -> implicit graph backend
  -> sample evidence projection
  -> local candidate subwalks
  -> local genotype scoring
  -> recombination/copying inference
  -> inferred phased haplotype mosaics
```

Syng is the first scalable backend. It is not the only backend.

## Sequence Namespace

The coordinate authority is a single namespace of named source sequences:

```text
source_sequence_id : [0, source_length)
```

Those source sequences may be chromosome-length assemblies, contigs, scaffolds,
fragments, reads, partial assemblies, pooled consensus sequences, or other
collection members. Graphs, GBWTs, syng indexes, alignments, packs, local
renders, and inferred mosaics are all views over this sequence coordinate
space.

This is the core synthesis: IMPG is a population-scale genotyping and
translation engine over source sequence coordinates. Homology graphs describe
relationships among those sequences; they do not replace the source namespace.

## Path Identity And Pan-SN

All backend and render work must preserve path identity. When source sequences
use Pan-SN names, the structured path identity is:

```text
sample#haplotype#contig
```

Every candidate interval, rendered graph path, GBWT path, evidence projection,
genotype result, and inferred mosaic segment must be able to recover the source
sequence id, sample, haplotype, contig, full path name, and source interval when
those fields exist. Non-Pan-SN names should be supported as explicit unparsed
source sequence names, but they must not be silently interpreted as biological
sample/haplotype metadata.

This is a scaling requirement, not just output formatting. The copying model
needs stable haplotype identity, cohort outputs need sample grouping, and local
graph renders need a compact translation back to source panel paths.

## Core Abstractions

The central invariant is that every backend must expose the same small set of
operations:

```text
GraphBackend
  graph_id()
  feature_space()
  candidates(target range or subgraph) -> candidate subwalks/haplotypes
  candidate_features(candidate) -> feature vector
  candidate_walk(candidate) -> ordered oriented feature walk
  fetch_sequence(candidate) -> DNA sequence, when available
  project_evidence(reads/alignments) -> pack + optional read-walk/link evidence
```

The current syng implementation already has most of these pieces, but they are
still named and wired as syng-specific internals:

```text
feature_space         syng-syncmer-node
candidate discovery   syng query over path ranges
candidate features    syncmer-node traversal counts
candidate walk        signed syncmer nodes with bp positions
sample evidence       pack from impg map
read linkage          GAF syncmer walks + syng GBWT MEMs
```

Variation graphs should fit the same abstraction:

```text
feature_space         variation-graph-node, path-step, bubble-allele, or segment
candidate discovery   paths through a local graph/subgraph
candidate features    node/path-step/bubble traversal counts
candidate walk        ordered oriented graph walk
sample evidence       graph alignment, local realignment, MEM projection, or pack
read linkage          read GAF/GAM/BAM projection through ordered graph walks
```

## Graph Construction Backends

We should not block on one graph builder. IMPG can consume or invoke multiple
construction paths.

### Syng Backend

Use the existing syng index for whole-panel implicit graph operations:

```text
FASTA/AGC -> impg syng -> .1gbwt/.1khash/.spos/.pstep/.names/.meta
FASTQ     -> impg map -o proj -> sample.pack + reads.gaf.zst
```

This path is best for large cohorts and whole-genome scale because the graph is
implicit and compressed.

### Local Variation Graph Backend

Build a local graph from candidate sequences extracted around a target region:

```text
target range
  -> candidate haplotype sequences from syng query / panel paths
  -> pggb, seqwish/smoothxg, or an internal sequence-graph pipeline
  -> local GFA/GBZ
  -> project reads or packs to local graph features
```

This is the right place to support LikeGT/COSIGT-style variation-graph
genotyping, local realignment, ancient DNA likelihoods, and bubble-level models.

### Direct Sequence Backend

For some loci, a graph may be unnecessary:

```text
candidate haplotype FASTA
  -> align reads/MEMs to candidate haplotypes
  -> candidate sequence support features
  -> genotype / copy-number / mosaic inference
```

This is useful for quick experiments, ancient DNA, and difficult loci where
local haplotype realignment is more sensitive than graph-node projection.

## Evidence Projections

Evidence projection is separate from scoring. Every projection should produce a
typed feature vector and, when possible, read-level linkage.

Supported and planned projection types:

```text
syncmer-pack
  FASTQ -> syng syncmer matches -> node-count pack

syncmer-walk
  FASTQ -> signed syncmer GAF with query positions -> GBWT MEM read evidence

variation-graph-align
  BAM/CRAM/GAF/GAM -> graph-node support and read walks

haplotype-align
  FASTQ/BAM/CRAM -> local alignment to candidate haplotype sequences

mem-project
  reads -> MEMs on graph/path backend -> node/subwalk support

hybrid-local
  global syng map gathers candidates/read subsets, then local graph alignment
```

The projection object should evolve beyond the current `sample.proj/` directory:

```text
sample.proj/
  manifest.json
  sample.pack
  reads.gaf.zst
  optional local-read index / read locator
  optional alignment summaries
  optional scorer-ready feature tensors
```

The manifest must identify:

```text
graph backend and graph identity
feature space
projection method
sample/read source
parameters
software versions
```

## Scoring Methods

Scoring should be pluggable and small-locus friendly.

### Cosine (`cos`)

Current production method:

```text
sample vector x
candidate genotype vector g = h1 + h2 + ...
score = cosine(x, g)
```

Strengths:

- simple and auditable
- works with pack counts
- catches copy-number dosage when features are informative
- graph-backend neutral if features match

Limitations:

- ignores base qualities and damage models
- mostly unoriented unless read-walk evidence is layered in
- does not model sampling variance explicitly

### Count Likelihood

A direct probabilistic count model:

```text
x_f ~ Poisson(depth * g_f + error_f)
or
x_f ~ NegativeBinomial(depth * g_f, dispersion)
```

This should improve dosage/CNV calls and support calibrated likelihoods.

### Alignment Likelihood

LikeGT/COSIGT-style graph or haplotype alignment scoring:

```text
read alignments -> per-read likelihood over candidate haplotypes/subwalks
genotype likelihood = product/sum over read likelihoods
```

This is essential for ancient DNA, short damaged reads, and regions where
syncmer matching loses sensitivity.

### GBWT / Panel Imputation

Use panel path structure as a prior over copied haplotype mosaics:

```text
emission = local evidence score
transition = probability/cost of copying same or nearby panel haplotype
state = one or more panel paths/subwalks at a genomic block
```

This is the principled version of the current beam stitcher.

### Learned Scoring

Learned local emission scoring is feasible because the inputs are small per
locus. The first learned model should score local candidate diplotypes, not own
the whole recombination problem.

Recommended first model:

```text
input:
  observed feature vector
  candidate haplotype feature vectors
  optional read-link/read-walk summaries
  optional depth/error/context features

output:
  local score for candidate genotype / candidate subwalk pair
```

Training data comes from simulation:

```text
choose panel haplotypes or recombinant mosaics
simulate reads with coverage/error/damage parameters
project reads through selected backend
train scorer to recover generating genotype/subwalks
```

Then genome inference remains an HMM/beam/Viterbi problem:

```text
learned scorer = local emission
copying model  = transition
decoder        = best phased mosaic
```

This keeps the model interpretable and prevents a small neural network from
having to learn the entire recombination state space.

## Recombination And Copying Model

The biological object we want is a pair of phased haplotypes represented as
copied subwalks from the panel graph, with possible recombinations and CNVs.

A practical model:

```text
block j state:
  phase 0 copies candidate subwalk a_j
  phase 1 copies candidate subwalk b_j

emission:
  score(evidence_j | a_j, b_j)

transition:
  same source haplotype continuation is cheap
  switching source haplotype costs recombination/copy switch
  read-walk links reward phase-consistent adjacent states
  GBWT/path priors reward population-supported transitions
```

Initial decoders:

```text
beam search over top local genotype states
Viterbi over a bounded candidate lattice
```

Later decoders:

```text
multi-resolution blocks
adaptive block boundaries from graph ambiguity/read support
CNV-aware states with variable copy count
population prior from GBWT/haplotype panel
```

The current `infer --stitch beam` is an early version of this. It already has
local emissions, switch penalties, and read-walk link rewards. The next step is
to make the state space and transition model backend-neutral and explicit.

## Local Graph Workflow For Hard Loci

For C4/C4A/C4B-like loci and other structurally complex regions:

```text
1. use syng/HPRCv2 to gather homologous candidate haplotypes
2. extract candidate sequences from AGC/FASTA
3. build a local variation graph with pggb or the internal sequence pipeline
4. project sample evidence to the local graph
5. run multiple scorers:
   - syng pack/cos
   - local graph node/cos
   - alignment likelihood
   - learned local emission
6. decode a phased copied mosaic
7. emit:
   - local genotype table
   - mosaic TSV
   - diplotype FASTA
   - diplotype GFA
   - debug graph/evidence bundle
```

This lets us keep the whole-genome index implicit while using explicit local
graphs only where they add value.

The render and translation layer is described in
[`render-gbz-translation-design.md`](render-gbz-translation-design.md). That
document is the implementation anchor for Pan-SN-aware GBZ-style bundles,
succinct translation tables, syng-native renders, local pggb/seqwish renders,
and tiled whole-genome rendering.

## Implementation Phases

### Phase 1: Make Current Syng Model Explicit

- Rename internal concepts around generic evidence/candidate/scorer terms.
- Keep current CLI behavior stable.
- Expose a backend-neutral local call structure.
- Keep syng candidate discovery as the first backend.
- Preserve `cos` and read-walk/GBWT-MEM evidence behavior.

### Phase 2: Backend Interface

- Add a `GraphBackend`-like trait or module boundary.
- Move syng-specific candidate discovery behind `SyngBackend`.
- Define generic candidate subwalk, feature vector, and read-link evidence
  structs.
- Keep existing `impg genotype cos` and `impg infer` as wrappers over the
  generic pipeline.

### Phase 3: Local Variation Graph Backend

- Add local graph construction from candidate sequences.
- Support pggb/seqwish external invocation first, then internal sequence graph
  construction where practical.
- Project reads/alignments to local graph nodes.
- Produce graph-node packs and read walks.
- Run existing `cos` scorer on variation graph nodes.

### Phase 4: Alignment Evidence

- Add BAM/CRAM/GAF input ingestion.
- Collect reads by target partition/reference interval.
- Realign locally to candidate haplotypes or local variation graph.
- Produce pack-like support and read-level likelihood summaries.

### Phase 5: Probabilistic And Learned Scorers

- Add count-likelihood scorer.
- Add alignment-likelihood scorer.
- Add simulator and scorer training/evaluation harness.
- Add first learned emission scorer over local feature tensors.

### Phase 6: Copying/Imputation Decoder

- Replace ad hoc stitch scoring with a formal state/transition model.
- Add GBWT/panel priors where available.
- Support adaptive block boundaries and multi-resolution inference.
- Emit phased copied subwalks and confidence annotations.

### Phase 7: Validation

- Synthetic adversarial tests remain mandatory.
- Add real C4/C4A/C4B debug dataset if available or construct it from HPRCv2.
- Add HG002 truth-backed local validations.
- Add HPRCv2 whole-genome smoke tests for pack/proj/infer.
- Track accuracy, runtime, memory, and output stability.

## Testing Matrix

Synthetic tests should cover every backend/scorer combination we claim:

```text
backends:
  syng
  local variation graph
  direct haplotype sequence

projection:
  syncmer pack
  syncmer walk / GBWT MEM
  graph alignment
  haplotype alignment

events:
  SNP/short indel
  insertion
  deletion
  inversion
  nested SV
  duplicated copy
  triplicated copy
  swapped paralogs
  low coverage
  noisy reads
  ancient-DNA-like damage
  recombination inside a target
  recombination between targets

outputs:
  local genotype table
  mosaic TSV
  diplotype FASTA
  diplotype GFA
```

The current syng tests already cover many of these for syncmer pack/projection
evidence. The missing coverage is variation-graph projection, alignment-based
evidence, learned scoring, and truth-backed real loci.

## Guiding Principle

Do not make a special-purpose C4 genotyper, a special-purpose syng genotyper,
and a special-purpose variation-graph genotyper. Build one inference system
with interchangeable backends and scorers.

The scalable core is implicit:

```text
large panel graph stays compressed
sample evidence is compact and typed
local explicit graphs are materialized only when needed
```

That is the reason IMPG is the right place for this work.
