# impg: implicit pangenome graph

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/impg/README.html)

## Why impg?

Studying variation at specific loci across populations or species usually
means either building an expensive whole-genome graph or falling back to
a single reference. `impg` takes a third path: it treats all-vs-all
pairwise alignments as an *implicit pangenome graph* and projects
target ranges through the alignment network to extract only the
homologous sequences you need. Query regions across hundreds of
genomes in seconds, walk transitive alignments, partition a cohort into
comparable loci, refine regions to maximize sample support — all
without ever materializing a graph structure.

## What does it do?

At its core, `impg` lifts ranges from a target sequence (the reference
in a given alignment) into the queries aligned onto it. It outputs
BED / BEDPE / PAF — ready to feed FASTA extraction, multiple sequence
alignment, or a graph builder like `pggb` or `minigraph-cactus` — and
can also emit GFA directly by chaining sweepga + seqwish + smoothxg-style
smoothing.

`impg` also embeds [syng](https://github.com/richarddurbin/syng) as a second backend: build a syncmer GBWT index once from a FASTA or AGC, then `query` / `partition` against it without ever running an aligner. Shared syncmers act as anchors; transitive chains define homology. Same output formats as the alignment path. `impg map` projects short reads onto a syng index in GAF or PAF form.

## How does it work?

`impg` uses [coitrees](https://github.com/dcjones/coitrees) (cache-oblivious
interval trees) for fast range lookup, and stores CIGAR strings as
compact deltas. The result is fast, memory-efficient projection of
sequence ranges through alignment networks.

## Install

```bash
# Bioconda
conda install -c bioconda impg

# Source
git clone --recursive https://github.com/pangenome/impg.git
cd impg
cargo install --force --path .
```

The source install places `impg`, `gfaffix`, and the companion aligner
binaries (`wfmash`, `FastGA`) into `~/.cargo/bin/`.

### Docker

```bash
docker pull pangenome/impg
docker run pangenome/impg
```

### Troubleshooting

On older glibc systems (e.g. Debian Buster) a plain `cargo build` can
fail because wfmash needs modern CMake / GCC / glibc. Use Guix's
toolchain:

```bash
source ./env.sh
cargo build --release
```

See `.guix/` for the Guix build recipe (`guix build -L .guix/modules
--file=guix.scm`). For `libclang` link errors, set `LIBCLANG_PATH` to
your LLVM install (see `env -i … LIBCLANG_PATH=…`).

## Quick start

```bash
impg query -a cerevisiae.pan.paf.gz -r S288C#1#chrI:50000-100000 -d 1000 -x
```

- `-a` — alignment file (PAF / 1ALN / TPA). PAF must use `=`/`X` CIGAR ops
  (from `wfmash` or `minimap2 --eqx`).
- `-r` — target range, `seq:start-end`.
- `-d` — merge query-gathered ranges separated by up to this many bp.
  This is also the largest internal gap/SV one query hop can absorb into
  one reported interval.
- `-x` — walk the transitive closure: find everything aligned to the
  initial result, recursively.

Example output (BED):

```
S288C#1#chrI         50000   100000
DBVPG6044#1#chrI     35335    85288
Y12#1#chrI           36263    86288
DBVPG6765#1#chrI     36166    86150
YPS128#1#chrI        47080    97062
UWOPS034614#1#chrI   36826    86817
SK1#1#chrI           52740   102721
```

## Commands

All commands accept `-a` (alignment files, mixed PAF/1ALN/TPA) or
`--alignment-list` (text file, one per line), `-t` / `--threads`, and
`-v 0|1|2` for verbosity. Every command has a `--help` with the
exhaustive flag list — this section covers the flags you'll actually turn.

### `query` — project a range through alignments

```bash
# A single range
impg query -a aln.paf -r chr1:1000-2000 -d 100

# Transitive closure (depth 2 by default)
impg query -a aln.paf -r chr1:1000-2000 -d 100 -x -m 3

# Many regions from a BED, mixed PAF + 1ALN
impg query -a f1.paf f2.1aln -b regions.bed -d 100

# One local graph per BED4 row, named from column 4
impg query -a pan.syng -b regions.bed -d 100k -o gfa:syng:crush \
           --sequence-files genomes.agc -O graphs/

# Output formats: auto | bed | bedpe | paf | gfa | vcf | maf | fasta | fasta+paf | fasta-aln
impg query -a aln.paf -r chr1:1000-2000 -d 100 -o bed
impg query -a aln.paf -r chr1:1000-2000 -d 100 -o gfa --sequence-files genomes.fa
impg query -a aln.paf -r chr1:1000-2000 -d 100 -o vcf --sequence-files genomes.fa
impg query -a aln.1aln -r chr1:1000-2000 -d 100 -o fasta --sequence-files *.fa \
           --reverse-complement

# Filter / shape the result
impg query -a aln.paf -r chr1:1000-2000 --min-identity 0.9 -l 5000 -d 1000

# Restrict to a sequence whitelist (also filters transitive intermediates)
impg query -a aln.paf -r chr1:1000-2000 -d 100 -x --subset-sequence-list seqs.txt

# Fast approximate mode (.1aln only; bed/bedpe output)
impg query -a aln.1aln -r chr1:1000-2000 -d 100 --approximate

# Alignment-free path: -a accepts a syng index prefix (see `syng` below)
impg query -a pan.syng -r chr1:1000-2000 -d 100 --sequence-files genomes.fa
```

GFA / MAF / FASTA outputs need `--sequence-files` (FASTA or AGC
archive) or `--sequence-list`. See [GFA engines](#gfa-engines) for
engine selection and partitioned builds. With `-b` and graph-like outputs
(`gfa`, `vcf`, or `gbwt`), `-O` is treated as an output directory and
each BED row is written separately using the BED column 4 name as the
file stem.

### `graph` — build a pangenome graph from FASTA

Runs alignment + seqwish + (optional) smoothing, no pre-computed
alignment needed.

```bash
# Default pipeline: pggb (align → seqwish → smooth → gfaffix)
impg graph --sequence-files genomes.fa -g output.gfa -t 16

# Partitioned mode for large inputs (aligns once, then builds per-window)
impg graph --sequence-files genomes.fa -g output.gfa --gfa-engine pggb:10000

# Reuse an existing PAF instead of aligning
impg graph --sequence-files genomes.fa -g output.gfa --paf-file aln.paf

# Batch alignment to cap per-batch RAM (wfmash) or disk (FastGA)
impg graph --sequence-files genomes.fa -g output.gfa --batch-bytes 2G
```

`query -o gfa` and `graph` share the same engine code and flags — the
only difference is where the sequences come from (IMPG index +
sequence files for `query`; FASTAs directly for `graph`).

### `partition` — split the cohort into windowed loci

```bash
# 1Mb windows, single BED output
impg partition -a aln.paf -w 1000000 -d 100000

# One FASTA per partition (for downstream pipelines)
impg partition -a aln.1aln -w 1000000 -d 100000 -o fasta --sequence-files *.fa \
               --separate-files --output-folder partitions/

# Selection strategies pick the next starting sequence
impg partition -a aln.paf -w 1000000 -d 100000 --selection-mode longest     # default
impg partition -a aln.paf -w 1000000 -d 100000 --selection-mode sample      # PanSN sample
impg partition -a aln.paf -w 1000000 -d 100000 --selection-mode haplotype   # PanSN haplotype

# Start from a fixed list of sequences
impg partition -a aln.paf -w 1000000 -d 100000 --starting-sequences-file seqs.txt

# GFA output per partition; engines: pggb | seqwish | poa
impg partition -a aln.paf -w 1000000 -d 100000 -o gfa --gfa-engine pggb \
               --sequence-files *.fa --separate-files --output-folder gfas/

# Fully partitioned pipeline: build → lace → one gfaffix pass
impg partition -a aln.paf -w 100000 -d 100000 -o gfa --gfa-engine pggb:10000 \
               --sequence-files *.fa --output-folder results/
```

### `refine` — tighten a locus to maximize sample support

Explores asymmetric left/right expansions around each range, picking
the smallest window that keeps the most sequences, samples, or
haplotypes fully spanning it. Useful for anchoring loci outside
structural variants.

```bash
impg refine -a aln.paf -r chr1:1000-2000
impg refine -a aln.paf -b loci.bed --span-bp 2000 -d 200000

# Maximize PanSN samples / haplotypes instead of raw sequence count
impg refine -a aln.paf -r chr1:1000-2000 --pansn-mode sample
impg refine -a aln.paf -r chr1:1000-2000 --pansn-mode haplotype

# Cap expansion distance
impg refine -a aln.paf -r chr1:1000-2000 --max-extension 0.90    # 90% of locus
impg refine -a aln.paf -r chr1:1000-2000 --max-extension 50000   # 50kb absolute

# Emit the spanning-entity list alongside the refined BED
impg refine -a aln.paf -r chr1:1000-2000 --support-output support.bed
```

### `similarity` — pairwise similarity / distance within a region

```bash
impg similarity -a aln.paf -r chr1:1000-2000 --sequence-files *.fa
impg similarity -a aln.1aln -b regions.bed --sequence-files *.fa --distances

# Group by PanSN prefix
impg similarity -a aln.paf -r chr1:1000-2000 --sequence-files *.fa \
                --delim '#' --delim-pos 2     # sample#haplotype

# PCA / MDS on the distance matrix
impg similarity -a aln.paf -r chr1:1000-2000 --sequence-files *.fa \
                --pca --pca-components 3 --pca-measure cosine
```

### `lace` — combine many per-window graphs into one

```bash
# GFAs (auto-detected)
impg lace -f gfa1.gfa gfa2.gfa gfa3.gfa -o combined.gfa

# From a list file, fill inter-window gaps with sequence
impg lace -l gfa_list.txt -o combined.gfa --fill-gaps 1 --sequence-files ref.fa

# VCFs
impg lace -f *.vcf -o combined.vcf --reference ref.fa
```

Path names must follow `NAME:START-END` (e.g.
`HG002#1#chr20:1000-2000`); the coordinates drive reassembly. `NAME`
may contain `:` — the last `:` is the separator.

Recommended post-processing:

```bash
gfaffix combined.gfa -o combined.fix.gfa &> /dev/null
odgi unchop -i combined.fix.gfa -o - -t 16 | \
  odgi sort  -i -              -o - -p gYs -t 16 | \
  odgi view  -i -              -g > combined.final.gfa
```

### `index` — build a reusable IMPG index

```bash
# Single combined index
impg index -a aln.paf -i aln.impg

# Mixed PAF + 1ALN + TPA
impg index -a f1.paf f2.1aln f3.tpa -i all.impg

# Per-file index (faster incremental rebuilds for large cohorts)
impg index --alignment-list files.txt --index-mode per-file
```

`--index-mode auto` (default) picks per-file when ≥ 100 files are
listed, single otherwise. impg warns when the index is older than its
input alignments; `-f/--force-reindex` rebuilds.

Bgzipped PAFs work natively (`impg` reads `.paf.gz`); optionally
`bgzip -r alignments.paf.gz` creates a `.gzi` sidecar to speed up the
first read.

### `stats` — summarize alignments

```bash
impg stats -a aln.paf
impg stats -a f1.paf f2.1aln
```

### `syng` and `map` — alignment-free syncmer backend

`impg syng` builds a [syng](https://github.com/richarddurbin/syng) index from FASTA or AGC. Six sidecars are written under one prefix (`<prefix>.1khash` dictionary, `.1gbwt` GBWT, `.syng.names`, `.syng.pstep` sampled path-step checkpoints, `.syng.spos` sampled syncmer occurrence positions, `.syng.meta` parameters — auto-loaded on read). Any later `impg query` / `partition` / `map` can then point `-a` at the prefix (or any sidecar) and skip pairwise alignment. Exact target coordinates are located by walking GBWT occurrences forward to the next sampled checkpoint and resolving it through `.syng.spos`.

Parameters follow the syng paper: `--smer-length` (`s`, default 8) and `--syncmer-length` (`k`, must be odd, default 63). Position sidecars use a regular per-path syncmer-step grid plus the terminal syncmer: `--position-sample-rate 256` samples steps `0, 256, 512, ...` and the final step on each path. `--parallel-dictionary` adds a deterministic prepass for large inputs.

`impg map` projects FASTA/FASTQ queries onto a syng index via shared syncmers. The default output is GAF (per-read syncmer-node walks); pass `-o paf` for projected genome coordinates, `-o pack` for a compact binary node support vector, `-o pack-tsv` for a human-readable TSV support vector, or `-o proj` for a sample projection bundle containing both `sample.pack` and `reads.gaf.zst`. Text map output written with `-O` is compressed automatically when the filename ends in `.zst` or `.zstd`; `pack` uses internal block zstd compression for random access by node ID. `packbin` remains accepted as a compatibility alias for compact pack output.

Use `impg syng-repair -a <prefix> --position-sample-rate <N> --force` to rebuild or resample `.syng.pstep` and `.syng.spos` from an existing `.1gbwt` / `.1khash` syng index without re-reading the original sequences.

End-to-end walkthrough using ODGI's C4 test GFA (90 HPRC haplotypes, ~6.9 Mb total):

```bash
# 1. Get the C4 GFA from odgi/test/, dump paths to FASTA
curl -sL -O https://raw.githubusercontent.com/pangenome/odgi/master/test/chr6.C4.gfa
odgi paths -i chr6.C4.gfa -f > chr6.C4.fa
samtools faidx chr6.C4.fa

# 2. Build the syng index (~80 ms on this dataset, 4 threads)
impg syng -f chr6.C4.fa -o c4.syng \
          --syncmer-length 63 --smer-length 8 --syncmer-seed 7 \
          --position-sample-rate 256 -t 4
# writes c4.syng.1khash / .1gbwt / .syng.names / .syng.spos / .syng.pstep / .syng.meta

# 3. Query a 10 kb grch38 sub-window onto all 90 haplotypes
impg query -a c4.syng -r 'grch38#chr6:31972046-32055647:0-10000' \
           -d 150 --sequence-files chr6.C4.fa | head -3
# HG00673#1#JAHBBZ010000030.1:31835924-31919525  0  10000  grch38...:0-10000  .  +
# HG01123#2#JAGYYY010000050.1:31954985-32038586  0  10001  grch38...:0-10000  .  +
# HG01243#1#JAHEOY010000117.1:3252171-3329407    0   9999  grch38...:0-10000  .  +

# 4. Whole grch38 C4 path (length from .fai) → 36 hits, FASTA out
LEN=$(awk -F'\t' '$1=="grch38#chr6:31972046-32055647"{print $2}' chr6.C4.fa.fai)
impg query -a c4.syng -r "grch38#chr6:31972046-32055647:0-${LEN}" \
           -d 150 --sequence-files chr6.C4.fa -o fasta > c4.homologs.fa

# 5. Map a 2 kb probe back onto the index
samtools faidx chr6.C4.fa 'grch38#chr6:31972046-32055647:5000-7000' > probe.fa
impg map -a c4.syng -q probe.fa | cut -f1-6 | head -1
# grch38...:5000-7000  2001  16  1979  +  <264>265>266<267<268>...
impg map -a c4.syng -q probe.fa -o pack-tsv | head
# #node_id  count
# 264      1
# 265      1
impg map -a c4.syng -q probe.fa -o pack-tsv -O probe.pack.tsv.zst
impg map -a c4.syng -q probe.fa -o pack -O probe.pack \
         --pack-compression-level 12
impg map -a c4.syng -q probe.fa -o proj -O probe.proj
impg map -a c4.syng -q probe.fa -o paf
# grch38...:5000-7000  2001  1302  1833  +  HG02109#1#...  77232  6300  6831  126 531 0 an:i:2 sk:i:63
```

`pack` stores a dense u8 count vector over syng node IDs, zstd-compressed
in independently addressable blocks, plus a small overflow table for node
counts above 255. `--pack-compression-level` defaults to 12; level 19 is a
compact cohort/archive setting that is slower to write but just as fast to
read in practice.

A `proj` directory is the richer sample projection object:

```text
sample.proj/
  manifest.json
  sample.pack
  reads.gaf.zst
```

Use `pack` when aggregate graph support is enough or when storing many
samples. Use `proj` when inference should also have read-walk linkage evidence
available.

`impg genotype` (`impg gt`) is the namespace for graph-based genotyping
methods. The first listed method is `cos`: cosine genotyping over graph-feature
coverage, following the COSIGT/LikeGT scoring model:

```bash
impg genotype cos -a c4.syng -p sample.pack \
  -r 'grch38#chr6:31972046-32055647:0-10000' \
  --ploidy 2 --top-n 20
```

`cosigt` is accepted as an alias for `cos` to make the paper/tool lineage
clear. `cos` extracts haplotype candidates from the requested reference path
range, builds traversal-count vectors over the candidate syncmer nodes, and ranks
ploidy-sized haplotype combinations against the sample `pack`
coverage vector by cosine similarity. The default candidate mode is
`spanning`: candidates must have shared anchors spanning the requested
reference interval. Use `--candidate-mode overlapping` to score each gathered
candidate interval independently. See `docs/genotype-architecture.md` for the
graph-feature evidence model this is built around.

`impg infer` lifts the same pack/cos scoring path over ranges or partitions and
emits allele calls over genomic intervals:

```bash
impg infer -a hprc.syng --proj sample.proj \
  --partitions partitions.bed --top-n 1 -O sample.infer.tsv.zst
```

Inputs choose the operation mode: `-r/--target-range` types one range,
`--target-bed` types a BED-like range list, `--partitions` consumes ready
`impg partition` BED output, and omitting all three runs internal syng partition
discovery, which requires `-d/--merge-distance`. `infer` can consume `--pack`
for aggregate graph support or `--proj` for a bundle that also carries read
walks. Add `--stitch beam --emit-mosaic out.tsv --emit-fasta haps.fa
--emit-gfa diplotype.gfa --sequence-files panel.fa-or.agc` to stitch local
calls into phased mosaic paths and materialize sequence. When `--proj` or
`--gaf` read walks are available, stitched inference rewards phase transitions
supported by the same reads (`--read-link-weight`,
`--min-read-link-anchors`) while `--switch-penalty` keeps unnecessary panel
crossovers expensive. Add `--phase-block-size N` to split broad targets into
internal phase blocks before stitching, which lets the same model infer copied
segments inside a partition. See `docs/infer-design.md`.

`-r` splits on the **last** `:`. Path names from `odgi paths -f`
already contain coordinates (`grch38#chr6:31972046-32055647`), so a
whole-sequence query still needs an explicit trailing sub-range
(`:0-${LEN}`) — otherwise impg parses `grch38#chr6` as the sequence
and the bare `31972046-32055647` as the range, and the genome lookup
fails.

Default query mode runs BiWFA boundary refinement and so requires
`--sequence-files`. Pass `--syng-raw` for the raw syncmer-resolution
pass-through. Tune via `--syng-padding`, `--syng-min-chain-anchors`
(lower -> more paralog hits; default is an adaptive cap derived from
query length, syncmer density, and expected exact-syncmer survival at
95% identity), `--syng-min-chain-fraction` (default `0.5`; set `0` for
exploratory local-chain discovery), and the seed-frequency filters
`--syng-seed-max-occurrences` / `--syng-seed-drop-top-fraction`.
By default syng query drops the top 0.05% most frequent query-local seed
syncmers, seeds ranges from bounded exact GBWT walk seeds of five
syncmers (`--syng-seed-walk-anchors`; set `3` for more sensitive seeds),
and does not use a fixed absolute occurrence cap; set
`--syng-seed-max-occurrences` only when a hard panel-specific ceiling is
desired. High-frequency syncmers are ignored only while seeding candidate
ranges; once a range survives, downstream scoring can still recover all
syncmers by walking that path range. `query -o gbwt` emits a
region-specific sub-GBWT. The syng prefix passed to `-a` is resolved
relative to cwd, so either `cd` to the index directory or pass an absolute
path. For a focused syng-backed local GFA recipe, see
[`docs/syng-gfa-query.md`](docs/syng-gfa-query.md).

For graph outputs, `query -o gfa|vcf --render-graph` also renders the final
1D graph with `gfalook`. The default image is PNG beside the `-O` output prefix
(`<prefix>.png`); `--render-graph-output` overrides the image path, and
`--render-graph-format svg` or a `.svg` suffix emits SVG. Add
`--render-graph-depth` to pass `-m` to `gfalook` for mean-depth coloring. With
`-b regions.bed` graph output, `-O` is a directory for the graph files and the
render output is written per BED row using the sanitized BED column 4 name.

If you already have a local GFA, call variants directly with POVU via
`impg gfa2vcf -g local.gfa -o local.vcf -r ref_path`. This is the same
GFA-to-VCF conversion used internally by `query -o vcf` and `partition -o vcf`.

## GFA engines

The `graph`, `query -o gfa`, and `partition -o gfa` commands share
one set of engine implementations, selected via `--gfa-engine`:

| Engine | Pipeline | Use for |
|---|---|---|
| `pggb` (default) | sweepga + seqwish + smoothxg-style smoothing + gfaffix | smoothed variation graphs |
| `seqwish` | sweepga + seqwish + gfaffix | raw (unsmoothed) graphs |
| `poa` | single-pass SPOA | small regions, quick MSA-based output |
| `syng` / `syng:blunt` | regional syng syncmer graph with exact zero-overlap path materialization | source-spelling syng graph output from syng indexes |
| `syng:raw` | regional syng syncmer overlap graph | explicit native overlap graph for debugging or overlap-aware consumers |
| `syng-local` / `syng-local:blunt` | extract query-selected sequences, build a fresh local syng graph, then exact zero-overlap materialization | experimental regional syncmer parameter sweeps |

For syng-index queries, `--gfa-engine syng` defaults to `syng:blunt`.
The compact form `-o gfa:syng:blunt,k=63,s=8,seed=7` is accepted as
shorthand for `-o gfa --gfa-engine syng:blunt,k=63,s=8,seed=7`; the
`k/s/seed` tail is checked against the loaded syng index.
Blunt syng output is sequence-preserving when the graph can fetch source DNA
for non-syncmer spans via `--sequence-files` or the local temporary FASTA used
by `syng-local`. Without source sequence files, missing gap DNA is emitted as
`N`, and those paths are explicitly not source-preserving. Use `syng:raw` when
you want the native syng overlap graph; raw path spelling requires an
overlap-aware parser and is not the same as concatenating S-line DNA.
Use `syng-local` when the local graph should be rebuilt from the extracted
regional sequences with its own syncmer scheme, for example
`-o gfa:syng-local:blunt,k=127,s=16,seed=7:crush`; in this mode `k/s/seed`
select the local rebuild parameters rather than asserting the global index.
Unlike `syng`, `syng-local` does not apply the syng frequency mask unless
`:mask` is requested explicitly.
Syng GFA extraction selects the top 0.05% most frequent local syncmer nodes and
private-splits unsupported high-frequency occurrences before raw or blunt graph
materialization. Occurrences in supported high-frequency runs or exact spans
stay shared (`freq-run=10` and `freq-span=1k` by default). The same
run/span rescue is applied to spectrum-selected dispersed scaffold-glue nodes,
so there is no second private-split path that bypasses the high-frequency
support policy. It also splits rare repeated-copy local syncmer contexts by
default, so high-copy repeats and single-syncmer repeat loops do not become
global graph glue. Use
`-o gfa:syng:nomask` to disable this, or
`-o gfa:syng:mask,top=0.001,max-occ=500,freq-run=10,freq-span=1k:crush` to
tune it. Use `freq-run-aware=false` or `legacy-freq-mask=true` to reproduce the
older node-level frequency removal for debugging.
Add `:cut-ns`, for example `-o gfa:syng:cut-ns:crush`, to drop assembly N-runs
from fetched gap DNA and split graph paths at those breaks.

For graph builds from query-extracted sequences, add an explicit terminal
N-run clipping stage before the engine: `-o gfa:cut-n=100:pggb`. This clips
only leading and trailing `N`/`n` runs whose length is at least the requested
threshold before graph construction; internal N-runs and shorter terminal
N-runs are preserved. There is no default clipping threshold, and omitting
`cut-n=<bp>` leaves extracted sequences unchanged.

The same shorthand works for the other engines: `-o gfa:pggb`,
`-o gfa:seqwish`, `-o gfa:poa`, and `-o gfa:syng`. Alignment-backed graph
builds may also include the aligner prefix, for example
`-o gfa:wfmash:seqwish`, `-o gfa:fastga:pggb`, or
`-o gfa:sweepga:seqwish`; this is equivalent to setting `--aligner` and
`--gfa-engine` separately. Add a generic `:crush` stage to run exact
path-preserving flubble/motif crush resolution over the induced blunt GFA,
then emit the finalized graph with default self-loop normalization and `Ygs`
sorting, for example `-o gfa:pggb:crush` or `-o gfa:seqwish:crush`.

For a C4/HPRCv2-style local render where terminal assembly N blocks would
otherwise become noisy graph tips:

```bash
impg query -a ~/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:<C4-range> \
  --sequence-files ~/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 100000 \
  -o gfa:cut-n=100:pggb \
  -O c4.cutn100.pggb.gfa
```

For a C4/HPRCv2-style local render where terminal assembly N blocks would
otherwise become noisy graph tips:

```bash
impg query -a ~/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:<C4-range> \
  --sequence-files ~/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 100000 \
  -o gfa:cut-n=100:pggb \
  -O c4.cutn100.pggb.gfa
```

VCF output uses the same graph engines and then converts the resulting local
GFA through POVU. Use `-o vcf --gfa-engine <engine>` or the shorthand
`-o vcf:<engine>`, for example `-o vcf:syng`.

### Partitioned mode

Append `:WINDOW` to any engine to build per-window and lace:

```bash
impg query    -a aln.paf -r chr1:0-500000 -o gfa \
              -d 1000 --gfa-engine pggb:10000 --sequence-files *.fa -O out
impg graph    --sequence-files *.fa -g out.gfa --gfa-engine seqwish:10000
impg partition -a aln.paf -w 100000 -d 100000 -o gfa --gfa-engine pggb:10000 \
               --sequence-files *.fa --output-folder results/
```

Window size is in bp (≥ 1000). Partitioned mode is the recommended
approach for large regions — it caps peak memory and runs one final
gfaffix pass over the laced graph.

### Tuning

The flags below are available on all three GFA-producing commands.
Defaults match pggb's conventions; only tune if the default graph
doesn't meet your need.

```bash
# Seqwish induction
--min-match-len 23            # minimum transitive-match length
--transclose-batch 10000000   # batch size (reduce for lower memory)
--sparse-factor 0.0           # drop this fraction of input matches
--disk-backed                 # use disk-backed interval trees
--repeat-max / --min-repeat-dist

# Smoothxg-style smoothing (pggb only)
--target-poa-length 700,1100  # one pass per value
--max-node-length 100
--poa-padding-fraction 0.001

# Alignment filtering (sweepga, seqwish + pggb only)
--no-filter                   # skip post-alignment filtering
--num-mappings many:many      # plane-sweep cardinality
--scaffold-jump 50000         # scaffold chaining gap (0 = off)
--scaffold-mass 10000         # min scaffold chain length
--overlap 0.95
--min-aln-identity 0.9

# Aligner backend
--aligner wfmash              # default; alt: fastga
--sparsify auto               # wfmash-only; pair-selection heuristic
--map-pct-identity 90         # wfmash -p value
--fastga-frequency / --fastga-frequency-multiplier   # fastga-only

# Temp files (can be large)
--temp-dir /scratch/tmp       # explicit path
--temp-dir ramdisk            # → /dev/shm on Linux
```

Combining `--aligner fastga` with `--sparsify` or `--aligner wfmash`
with `--fastga-frequency` is rejected at parse time.

## Common options

- `-a / --alignment-files` — one or more PAF/1ALN/TPA files (can be `.gz`).
- `--alignment-list` — text file, one alignment path per line.
- `-i / --index` — existing IMPG index.
- `-f / --force-reindex` — rebuild even if the index is up-to-date.
- `-t / --threads` — default `4`.
- `-d / --merge-distance` — required for query/partition/refine/similarity:
  merge query-gathered ranges within this gap (bp). This is the largest
  internal gap/SV one query hop can absorb into one reported interval.
- `--no-merge` — disable merging.
- `--consider-strandness` — keep strands separate during merge.
- `--subset-sequence-list` — restrict results to listed sequences.
- `--unidirectional` — disable bidirectional alignment interpretation.

Sequence-requiring outputs (GFA/MAF/FASTA, `similarity`, `lace
--fill-gaps`) take `--sequence-files` (FASTA or AGC) or
`--sequence-list`.

## Tutorial: yeast pangenome graph

```bash
FASTA="cerevisiae.fa.gz"
PAF="cerevisiae.paf"
THREADS=16

# 1. Index
impg index -a "$PAF" -i yeast.impg -t "$THREADS"

# 2. Partition into 100kb windows, one FASTA per window
mkdir -p partitions gfas
impg partition -i yeast.impg -w 100000 \
    -d 100000 --sequence-files "$FASTA" -o fasta \
    --separate-files --output-folder partitions -t "$THREADS"

# 3. Build per-partition GFAs in parallel
ls partitions/*.fasta | xargs -P 4 -I {} bash -c '
    f="{}"; base=$(basename "$f" .fasta)
    impg graph --sequence-files "$f" -g "gfas/${base}.gfa" -t 4
'

# 4. Lace, filling inter-window gaps with reference sequence
find gfas -name "*.gfa" -size +0 | sort -V > gfa_list.txt
impg lace --file-list gfa_list.txt --sequence-files "$FASTA" \
    -o yeast.gfa --fill-gaps 2 -t "$THREADS"

# 5. Post-process with odgi
odgi build  -g yeast.gfa          -o yeast.og       -t "$THREADS"
odgi sort   -i yeast.og           -o yeast.sort.og  -O -p Ygs -t "$THREADS"
odgi layout -i yeast.sort.og      -o yeast.lay      -t "$THREADS"
odgi viz    -i yeast.sort.og      -o yeast.viz.png  -x 4000 -y 1000 -s '#'
odgi draw   -i yeast.sort.og      -c yeast.lay      -p yeast.draw.png
```

For modern inputs, you can replace steps 1–4 with a single
`impg graph --sequence-files "$FASTA" -g yeast.gfa --gfa-engine pggb:100000`.

## Visualizing FASTA alignments

`scripts/faln2html.py` renders the `fasta-aln` output into an
interactive HTML MSA using [react-msa](https://github.com/GMOD/JBrowseMSA)
or [ProSeqViewer](https://github.com/BioComputingUP/ProSeqViewer).

```bash
impg query -a aln.paf -r chr1:1000-2000 -d 100 -o fasta-aln --sequence-files *.fa \
  | python scripts/faln2html.py -i - -o alignment.html [--tool proseqviewer]
```

## Authors

Andrea Guarracino <aguarra1@uthsc.edu> · Bryce Kille
<brycekille@gmail.com> · Erik Garrison <erik.garrison@gmail.com>

## License

MIT.
