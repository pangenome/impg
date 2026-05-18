# Render, GBZ, And Translation Design

IMPG is a translation system over implicit pangenomes. The scalable object is
not a single materialized graph. The scalable object is the ability to move
between source sequence coordinates, graph features, evidence projections, and
inferred haplotypes without losing identity.

The root namespace is named biological sequences plus coordinates:

```text
source_sequence_id : [0, source_length)
```

A source sequence can be a chromosome-length assembly path, contig, scaffold,
fragment, read, partial assembly, pooled consensus, or any other sequence
collection member. Everything else in IMPG is a view over this coordinate
substrate:

```text
syng index             view of syncmer homology over source sequences
explicit pggb/seqwish  view of local sequence homology over source intervals
GBWT                   compressed walks through one graph view
pack/projection        sample evidence projected onto one feature view
inferred haplotype     copied/recombined sequence intervals from the namespace
```

This is what lets the same system operate on comparative genomes, huge panels
of very homologous haplotypes, fragmentary assemblies, short reads, ancient DNA,
and mixed collections. The graph is a representation of relationships among
source sequence intervals, not the only coordinate authority.

The target render product is a GBZ-style bundle:

```text
render bundle
  graph topology + node sequences
  GBWT haplotype path index
  source-sequence/coordinate/feature translation tables
  manifest with graph identity, parameters, and software versions
  optional evidence projections
```

The near-term format can be a directory or simple container. The design target
should remain a coherent GBZ-style object, not unrelated side files.

## Sequence Namespace And Pan-SN

Every render bundle must define one sequence namespace. Pan-SN is the preferred
structured naming convention inside that namespace when the source data supports
sample/haplotype/contig identity.

Panel path names must preserve Pan-SN identity wherever the source data follows
that convention:

```text
sample#haplotype#contig
```

Every render, local graph, GBWT path, translation table, genotype result,
mosaic output, FASTA output, and GFA/GBZ path must be able to recover:

```text
source sequence id
sample id
haplotype id / phase id
contig or chromosome id
source sequence/path name
source interval
```

This matters because sample/haplotype grouping is the biological identity of
the panel. A path interval copied from `HG00000#1#chr6` is not just a string
label; it is one phased haplotype of one sample. Recombination/copying models,
diploid output, validation, and cohort summaries all depend on keeping this
identity exact.

Non-Pan-SN names must still be first-class source sequences, but they should be
represented as explicit unknown or single-field path metadata rather than
silently parsed into incorrect sample/haplotype fields. Fragmentary sequences
and reads may have no sample/haplotype/contig decomposition at all; they still
live in the same sequence namespace and can still be projected onto graph
features.

## Translation

Translation is the bookkeeping that connects the feature spaces:

```text
source sequence interval
  <-> syng syncmer walk
  <-> extracted DNA interval
  <-> rendered graph path steps
  <-> rendered GBWT path id
  <-> evidence feature vector
  <-> inferred haplotype/mosaic segment
```

The translation layer must be succinct. A human-readable TSV export is useful
for debugging, but the primary representation should be compact tables with
integer IDs, sampled offsets, and block-local arrays.

Core tables:

```text
String table
  sample names
  haplotype names
  contig names
  full source sequence/path names
  collection names

Source sequence table
  source_sequence_id
  full_name_id
  sample_id or null
  haplotype_id or null
  contig_id or null
  collection_id or null
  source_length

Rendered path table
  rendered_path_id
  source_sequence_id
  source_start
  source_end
  render_block_id
  gbwt_path_id

View table
  view_id
  kind: syng | local-variation-graph | read-projection | inferred-mosaic
  feature_space_id
  graph_id or null
  parameter_set_id

Feature-space table
  feature_space_id
  kind: syng-syncmer-node | variation-graph-node | haplotype-segment | ...
  view_id

Step translation table
  rendered_path_id
  sampled rendered step offset
  source bp position
  optional syng node offset
```

The step translation table should be sampled, not fully expanded by default.
Dense exports can be produced for debugging or small local graphs.

## Why GBZ-Style

GFA is useful but incomplete. It gives explicit topology and node sequences,
but it is not the compact haplotype path index.

GBWT is useful but incomplete. It gives compressed path walks, but it needs the
graph topology and metadata to make those walks interpretable.

A GBZ-style bundle is the natural render target because it packages the graph,
the path collection, and the identity/translation layer together.

Near-term bundle layout:

```text
render.impg-gbz/
  manifest.json
  graph.gfa.zst or graph.bin
  paths.gbwt
  translation.bin
  translation.tsv.zst          optional debug export
  path-names.tsv.zst           optional debug export
  evidence/                    optional projection files
```

Longer-term output may be a true GBZ-compatible file if we want direct
interoperability with GBZ tooling. The IMPG-specific translation tables remain
needed even when the graph and GBWT are GBZ-compatible.

## Render Modes

### Syng-Native Render

Render the syng graph itself:

```text
syng .1gbwt + .1khash + .names + .spos/.pstep
  -> explicit syncmer-node graph
  -> GBWT paths over syncmer nodes
  -> translation to source path intervals and syncmer positions
```

This is the most direct way to inspect or export the implicit syng graph. The
feature space remains `syng-syncmer-node`.

### Local Sequence Graph Render

Use syng as the global implicit oracle, then materialize an explicit local
sequence graph only for the target region:

```text
syng query target
  -> candidate source sequence intervals with optional Pan-SN metadata
  -> fetch interval DNA from AGC/FASTA
  -> pggb/seqwish/internal graph build
  -> thread input haplotypes as graph paths
  -> write graph-step translation tables
```

This is the first target for C4/C4A/C4B-like loci.

Implemented local render engines:

```bash
impg render -a panel.syng -r sample#0#chr6:31972057-32055418 \
  --sequence-files panel.fa -O c4.poa.impg-gbz --engine poa
impg render -a panel.syng -r sample#0#chr6:31972057-32055418 \
  --sequence-files panel.fa -O c4.seqwish.impg-gbz --engine seqwish
impg render -a panel.syng -r sample#0#chr6:31972057-32055418 \
  --sequence-files panel.fa -O c4.pggb.impg-gbz --engine pggb
```

The current local graph bundle writes `rendered.fa`, `graph.gfa`,
`namespace.json`, `translation.bin`, and `translation.tsv`. Its feature space is
`gfa-segment`. A true local GBWT over the explicit graph paths is still a
separate pending step; the rendered FASTA path names and translation tables are
already organized so that GBWT path IDs can be attached without changing the
source namespace model.

Syng-native bundles can be used as local genotyping views by mapping reads to
the bundle's `paths` syng index and passing the bundle to `genotype cos`:

```bash
impg map -a c4.syng-native.impg-gbz/paths -q reads.fq.gz -o pack \
  -O c4.local.pack
impg genotype cos --render-bundle c4.syng-native.impg-gbz --pack c4.local.pack
impg infer --render-bundle c4.syng-native.impg-gbz --pack c4.local.pack
```

When `-r/--target-range` is omitted, `genotype cos` and `infer` use the source
target from the render manifest and project it onto the containing rendered path.
A different source-coordinate subrange can be supplied with `-r`; it is
projected through the same translation table before scoring. `infer
--render-bundle` currently supports this single-target mode; BED and partition
inputs should be projected through the bundle as a deliberate batch operation in
a later slice.

### Whole-Genome Render

Whole-genome rendering should be tiled:

```text
partition genome
  -> render each tile
  -> remap node IDs into disjoint ranges
  -> merge graph chunks
  -> merge GBWT path collections
  -> merge translation tables
```

Node ID overlap is mechanical. The required invariant is that merged paths and
translation records still recover the original source sequence identity,
coordinates, and Pan-SN metadata when present.

## Proposed CLI

Start regional, then scale to partitions and whole genome:

```bash
impg render \
  -a panel.syng \
  --range sample#0#chr6:31900000-32100000 \
  --sequence-files HPRC.agc \
  --engine syng-native|seqwish|pggb \
  --emit-gbz c4.render.impg-gbz
```

Partitioned render:

```bash
impg render \
  -a panel.syng \
  --partitions parts.bed \
  --sequence-files HPRC.agc \
  --engine seqwish \
  --emit-gbz hprc.render.impg-gbz
```

Debug exports should be optional:

```text
--emit-gfa
--emit-gbwt
--emit-translation-tsv
--keep-intermediates
```

## Relationship To Genotyping And Inference

`impg genotype` and `impg infer` should consume typed render/projection
objects instead of reimplementing path lookup for each backend.

The intended flow:

```text
render bundle
  -> backend-neutral candidate subwalks
  -> backend-neutral feature vectors
  -> pack/projection evidence in the same feature space
  -> local scorer
  -> copying/imputation decoder
  -> Pan-SN-aware mosaic output
```

Syng packs and syng read walks remain first-class evidence. Local graph
alignment, haplotype alignment, and MEM projection become additional evidence
projection methods into declared feature spaces.

Population-scale genotyping is then a translation and inference problem:

```text
sample evidence
  -> project into one or more feature views
  -> translate features to candidate source sequence intervals
  -> score local copied haplotype states
  -> decode a phased mosaic over the source sequence namespace
```

This is why the design has to stay implicit and succinct. At cohort scale, the
system cannot depend on fully materializing every graph view, every path step,
or every evidence relationship. It should materialize local explicit graphs
only when a locus needs that representation.

## Implementation Steps

1. Define `SequenceNamespace`, `SourceSequenceId`, `SourceInterval`, and
   `PathIdentity` in one shared module.
2. Implement Pan-SN parsing/roundtrip behavior as optional metadata on source
   sequences, not as a requirement for membership in the namespace.
3. Define render bundle manifest and succinct translation table schemas.
4. Implement regional syng-native render with translation tables.
5. Implement regional local sequence graph render using POA, seqwish, and pggb,
   with threaded input haplotype paths and graph-step translation tables.
6. Teach `genotype`/`infer` to consume render bundles as a backend.
7. Add tiled render and merge.
8. Add true GBZ-compatible packaging if needed for interoperability.
9. Add local GBWT construction for explicit graph paths.

## Tests

Required tests:

- Pan-SN parsing and non-Pan-SN fallback.
- Source namespace roundtrip for chromosome assemblies, contigs, fragments,
  reads, and pooled/collection sequences.
- Source path identity survives syng build/load/render.
- Rendered path identity survives local graph construction.
- Translation maps source interval to rendered path steps and back.
- GBWT path IDs recover the correct source Pan-SN path metadata.
- Diploid sample grouping remains correct when two haplotypes share a contig.
- Duplicate contig names across samples do not collide.
- Fragmentary and partial sequences can be projected without Pan-SN metadata.
- Tile merge preserves source path identity and coordinate ranges.
- Genotype/infer outputs retain Pan-SN names after consuming a render bundle.
- C4/HPRCv2 local render smoke test once the real dataset path is stable.

The tests should include small synthetic Pan-SN FASTA inputs and at least one
real HPRCv2/C4 smoke path when available.
