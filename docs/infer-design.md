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
FASTQ + impg map -a panel.syng -o pack -> sample.pack
FASTQ + impg map -a panel.syng -o proj -> sample.proj/
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

The initial implementation supports pack-based evidence and projection bundles.
BWA/local-realignment is a backend, not a separate genotyping model.

## Partitions

Partitioning defines where to type. It is not a syncmer-ID windowing scheme.

The first implementation reuses `impg partition` semantics:

```bash
impg partition -a panel.syng -w 1000000 -d 10000 -o bed --output-folder parts
impg infer -a panel.syng --proj sample.proj --partitions parts/partitions.bed
```

`infer` may also run the partition step internally:

```bash
impg infer -a panel.syng -p sample.pack --window-size 1000000 -d 10000
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

## Projections

A projection is the native sample evidence object produced by `impg map -o
proj`. The first directory-backed format contains:

```text
sample.proj/
  manifest.json
  sample.pack
  reads.gaf.zst
```

`sample.pack` is the aggregate support vector. `reads.gaf.zst` is the per-read
syncmer walk layer used by stitched inference to supply linkage/phase evidence.
Commands should prefer `--proj` when both layers are available, while still
accepting `--pack` for large cohort-scale aggregate workflows.

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

The first implementation emits independent local calls per range/partition and
can optionally stitch retained local states into a phased mosaic.

The next layer should stitch those local calls into a haplotype mosaic:

```text
state_j = haplotype tuple in window j
emission_j(state_j) = local evidence score
transition(state_j-1, state_j) = panel/GBWT continuity + switch penalty
best genome = argmax sum_j emission_j + transition_j
```

This can start as beam search over the top local genotype states, then move to
more principled HMM/Viterbi inference.

Current stitching is a beam-search haplotype-copying model over retained local
calls. Each local genotype is expanded into ordered ploidy states. Same-path
continuation is cheap, haplotype switches pay `--switch-penalty`, and read walks
from a projection bundle add transition reward when the same read supports the
specific previous-candidate -> current-candidate phase link.

The read-link reward is intentionally simple and auditable:

```text
candidate feature = syng syncmer node
read hits candidate = shared read-walk nodes >= --min-read-link-anchors
transition support = adjacent candidate hits from the same read walk
transition reward = --read-link-weight * 10 * log10(1 + normalized_anchor_support)
```

Ambiguous reads divide their support across all candidate links they match, so
shared sequence does not create artificial certainty. The mosaic TSV reports
`transition_cost`, `read_link_reads`, `read_link_anchors`, and
`read_link_reward` per phase.

This is now a real read-link imputer across the current partition lattice. The
same state/transition machinery can also be run within a large target by adding
`--phase-block-size N`, which splits each target/partition into fixed-width
internal phase blocks before local genotyping and stitching. This lets `infer`
discover copied IBD segments inside a broad partition and then carry phase
between adjacent partitions with the same model. The fixed grid is an initial
lattice; future work can choose block boundaries from read density, graph
breakpoints, or local ambiguity instead of a constant width.

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

With `--stitch beam`, side outputs can be requested:

```bash
impg infer -a panel.syng --proj sample.proj --partitions parts/partitions.bed \
  --stitch beam --emit-mosaic sample.mosaic.tsv \
  --sequence-files panel.agc \
  --emit-fasta sample.haps.fa --emit-gfa sample.diplotype.gfa
```

The FASTA writer concatenates phase tracks and inserts `N`s at uncertain joins
unless `--strict-stitch` is set. The GFA writer emits one segment per selected
haplotype interval and path lines for the inferred phases.

## Testing Strategy

Tests should cover:

- explicit single range
- BED/range list input
- ready partition BED input
- internal partition discovery requiring `-d`
- pack, pack-tsv, compressed pack-tsv, compatibility packbin, and projection evidence
- short-read simulation against a panel with a decoy haplotype
- invalid CLI combinations
- output compression and deterministic row shape
- stitched recombinant phasing where read links disambiguate copied haplotype
  segments across both explicit partitions and fixed internal phase blocks
- read-link control settings, including disabled read-link rewards
- CNV/repeated-syncmer paths where a duplicated-copy haplotype creates repeated
  nodes in read GAF walks and must not collapse into a single-copy or unrelated
  allele
- triplicated copy-number paths where three-copy evidence beats one-copy,
  two-copy, and unrelated triplicated decoys
- nested insertion+deletion haplotypes under sparse noisy reads
- swapped paralogous-copy haplotypes where read links must avoid collapsing to
  homo-copy decoys

The synthetic graphs are not substitutes for real population graph validation.
They should progressively add real HPRCv2 truth-backed validation, richer
coverage gradients, sequencing error models, and orientation-confounded
inversion/paralogy decoys. The latter is a known model pressure point because
the current pack score is mostly unoriented syncmer-node support; fully
resolving it likely requires using oriented read-walk evidence in local emission
as well as transition scoring.

The broader roadmap for turning this into a backend-neutral pangenome genotyping
and phased genome inference system is in
[`pangenome-genotyping-roadmap.md`](pangenome-genotyping-roadmap.md).
The GBZ-style render bundle, translation layer, and Pan-SN path identity
requirements are described in
[`render-gbz-translation-design.md`](render-gbz-translation-design.md).
