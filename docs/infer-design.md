# IMPG Genome Inference Design

`impg infer` is the genome-scale layer above `partition`, `map`, `pack`, and
`genotype`. Its job is to turn sample evidence into allele calls over genomic
ranges, and eventually into a coherent inferred diploid or polyploid genome.

The key design rule is that every evidence backend produces a pack-like vector
in a declared feature space. Genotyping and inference consume those vectors;
they should not care whether the evidence came from native syncmer mapping,
BWA-aln reads, local realignment, GAF, or MEM hits.

## Core Pipeline

```text
reads / alignments
  -> evidence projection
  -> typed feature pack
  -> partition/range selection
  -> local candidate haplotypes
  -> local genotype scoring
  -> optional genome-wide imputation / stitching
  -> allele calls over ranges
```

For modern WGS, the first-class evidence path is native:

```text
FASTQ + impg map -a panel.syng -o packbin -> sample.packbin
```

For ancient DNA or very short degraded reads, the intended path is:

```text
BWA-aln/MEM BAM/CRAM
  -> collect reads by partition/reference interval
  -> local realign against candidate haplotypes
  -> local pack
```

The second path should still produce the same abstraction: a feature support
vector with enough metadata to say how it was made.

## Inputs Drive Operation

`infer` should choose its work pattern from the inputs:

1. `--target-range`: type exactly one coordinate range.
2. `--target-bed`: type each coordinate range from a BED-like file.
3. `--partitions`: type each ready partition from an `impg partition` BED.
4. No explicit ranges/partitions: run partition-style discovery first.

Discovery mode must require `-d/--merge-distance` (or an explicit no-merge
mode if we add one later). This is the same biological knob as in
`query`/`partition`: it defines the maximum internal gap/SV that one query hop
can absorb into one typed range.

The initial implementation supports pack-based evidence. BWA/local-realignment
is a backend, not a separate genotyping model.

## Partitions

Partitioning defines where to type. It is not a syncmer-ID windowing scheme.

The first implementation reuses `impg partition` semantics:

```bash
impg partition -a panel.syng -w 1000000 -d 10000 -o bed --output-folder parts
impg infer -a panel.syng -p sample.packbin --partitions parts/partitions.bed
```

`infer` may also run the partition step internally:

```bash
impg infer -a panel.syng -p sample.packbin --window-size 1000000 -d 10000
```

Internal partition discovery writes temporary partitions, parses them back, then
types each partition. This keeps the first implementation aligned with existing
partition behavior and avoids inventing a parallel windowing model.

## Feature Packs

Current pack feature space:

```text
syng-syncmer-node: node_id -> count
```

Future pack feature spaces should include:

```text
local-haplotype-segment
variation-graph-node
reference-interval
mem-hit
```

Pack metadata should eventually record:

```text
feature_space
sample
graph/index identity
evidence_backend = impg-map | bwa-aln | local-align | gaf | mem
weighting / normalization
```

The important invariant is compatibility: candidate allele vectors and sample
evidence vectors must live in the same feature space.

## Scoring

The first scoring method is `cos`, the same cosine scorer used by
`impg genotype cos`.

For a local typing range:

```text
candidate haplotype h -> vector v_h
genotype G=(h1,h2,...) -> vector g = v_h1 + v_h2 + ...
sample evidence -> vector x
score = cosine(x, g)
```

Later scoring methods can include count likelihoods, damage-aware ancient-DNA
likelihoods, learned scoring functions, or GBWT/panel-based imputers.

## Genome-Wide Inference

The first implementation emits independent local calls per range/partition.

The next layer should stitch those local calls into a haplotype mosaic:

```text
state_j = haplotype tuple in window j
emission_j(state_j) = local evidence score
transition(state_j-1, state_j) = panel/GBWT continuity + switch penalty
best genome = argmax sum_j emission_j + transition_j
```

This can start as beam search over the top local genotype states, then move to
more principled HMM/Viterbi inference.

## Output

The first output is a TSV of allele calls over ranges:

```text
#impg infer
#evidence_backend pack
#score cos
#feature_space syng-syncmer-node
#rank  partition  chrom  start  end  method  ploidy  similarity  qv  haplotypes  regions  candidate_anchors  candidate_span_fractions
```

Each row is a local genotype state. Adjacent rows with compatible inferred
haplotypes can later be merged into longer genome segments.

## Testing Strategy

Tests should cover:

- explicit single range
- BED/range list input
- ready partition BED input
- internal partition discovery requiring `-d`
- pack, compressed pack, and packbin evidence
- short-read simulation against a panel with a decoy haplotype
- invalid CLI combinations
- output compression and deterministic row shape

The synthetic graphs are not substitutes for real population graph validation.
They should progressively add SNPs, indels, nested SVs, repeats/paralogy,
coverage variation, sequencing error, and recombination/mosaic samples.
