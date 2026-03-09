# impg: implicit pangenome graph

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/impg/README.html)

## Why impg?

Studying genomic variation at specific loci—disease genes, regulatory elements, structural variants—across populations or species traditionally requires either building expensive whole-genome graphs or using reference-based methods that miss variation. `impg` solves this by treating all-vs-all pairwise alignments as an *implicit pangenome graph*, rapidly projecting target ranges through the alignment network to extract only the homologous sequences you need. Query regions across many genomes in seconds. Perform transitive searches to discover connected sequences. Partition genomes into comparable loci. Refine regions to maximize sample coverage—all without constructing explicit graph structure. This makes pangenome-scale comparative genomics fast and practical.

## Usage

Here's a basic example:

```bash
impg query -a cerevisiae.pan.paf.gz -r S288C#1#chrI:50000-100000 -x
```

- `-a` specifies the path to the alignment file in PAF, 1ALN, or TPA format. PAF files must use CIGAR strings with `=` for matches and `X` for mismatches (e.g., from `wfmash` or `minimap2 --eqx`).
- `-r` defines the target range in the format of `seq_name:start-end`
- `-x` requests a *transitive closure* of the matches. That is, for each collected range, we then find what sequence ranges are aligned onto it. This is done progressively until we've closed the set of alignments connected to the initial target range.

Depending on your alignments, this might result in the following BED file:

```txt
S288C#1#chrI        50000  100000
DBVPG6044#1#chrI    35335  85288
Y12#1#chrI          36263  86288
DBVPG6765#1#chrI    36166  86150
YPS128#1#chrI       47080  97062
UWOPS034614#1#chrI  36826  86817
SK1#1#chrI          52740  102721
```

## Installation

You need Rust (`cargo`) installed. Then:

```bash
git clone https://github.com/pangenome/impg.git
cd impg
cargo install --force --path .
```

### Troubleshooting

If you encounter issues related to `libclang` during the build process, you may need to set specific environment variables to point to your LLVM installation. 

```shell
env -i HOME="$HOME" PATH="/usr/local/bin:/usr/bin:/bin:$HOME/.cargo/bin" LIBCLANG_PATH="/usr/lib/llvm-7/lib" BINDGEN_EXTRA_CLANG_ARGS="-I/usr/lib/llvm-7/lib/clang/7.0.1/include" bash -c 'cargo build --release'
```

Alternatively, install from Bioconda:

```bash
conda install -c bioconda impg
```

## Building with GNU Guix

**NOTE**: All paths are relative to the repository root. If you are not working
in the repository root, please update the paths to fit your work scenario.

To build `impg` with guix without making it available, run the following
command:

```sh
guix build -L .guix/modules --file=guix.scm
```

The `-L` option adds the `.guix/modules` directory to the front of the guile
load path. The `--file` option points to the `guix.scm` at the root of the
repository.

To build and "install" `impg` with guix, run:

```sh
guix install -L .guix/modules --file=guix.scm
```

## Commands

### Query

Query overlaps in the alignment:

```bash
# Query a single region
impg query -a alignments.paf -r chr1:1000-2000

# Query a whole sequence by name
impg query -a alignments.paf -r chr1

# Query multiple regions from a BED file (mix PAF and .1aln)
impg query -a file1.paf file2.1aln -b regions.bed

# Enable transitive overlap search
impg query -a alignments.paf -r chr1:1000-2000 -x

# Set maximum transitive depth (default: 2)
impg query -a alignments.1aln -r chr1:1000-2000 -x -m 3

# Filter by minimum gap-compressed identity
impg query -a alignments.paf alignments.1aln -r chr1:1000-2000 --min-identity 0.9

# Output formats (auto/bed/bedpe/paf/gfa/maf/fasta/fasta+paf/fasta-aln)
impg query -a alignments.paf -r chr1:1000-2000 -o bed
impg query -a alignments.1aln -r chr1:1000-2000 -o bedpe
impg query -a file1.paf file2.1aln -b chr1:1000-2000 -o paf

# Write output to file instead of stdout (using -O / --output-prefix)
impg query -a alignments.paf -r chr1:1000-2000 -o bed -O results       # creates results.bed

# gfa/maf/fasta output requires sequence files (--sequence-files or --sequence-list)
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --sequence-files ref.fa genomes.fa
impg query -a alignments.1aln -r chr1:1000-2000 -o maf --sequence-list fastas.txt
impg query -a file1.paf file2.1aln -r chr1:1000-2000 -o fasta --sequence-files *.fa

# fasta+paf combines FASTA and PAF output
impg query -a alignments.paf -r chr1:1000-2000 -o fasta+paf --sequence-files *.fa

# fasta-aln outputs POA-based FASTA alignment
impg query -a alignments.1aln -r chr1:1000-2000 -o fasta-aln --sequence-files *.fa

# Works with AGC archives too
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --sequence-files genomes.agc

# fasta output with reverse complement for reverse strand sequences
impg query -a alignments.1aln -r chr1:1000-2000 -o fasta --sequence-files *.fa --reverse-complement

# Merge nearby regions (default: 0)
impg query -a file1.paf file2.1aln -r chr1:1000-2000 -d 1000

# Filter results by minimum length
impg query -a alignments.paf -r chr1:1000-2000 -l 5000

# Use DFS instead of BFS for transitive search (slower but fewer overlapping results)
impg query -a alignments.paf -r chr1:1000-2000 -x --transitive-dfs

# Fast approximate mode for .1aln files (bed/bedpe only)
impg query -a alignments.1aln -r chr1:1000-2000 --approximate

# Fast approximate mode with transitive queries (requires --min-transitive-len > trace_spacing)
impg query -a alignments.1aln -r chr1:1000-2000 --approximate -x --min-transitive-len 101

# Restrict results to sequences listed in a file (one name per line)
# Also filters intermediate steps during transitive queries
impg query -a alignments.paf -r chr1:1000-2000 -x --subset-sequence-list sequences.txt

# Transform coordinates back to original sequences when using subsequence inputs (seq_name:start-end)
impg query -a alignments.paf -r chr1:1000-2000 -o paf --original-sequence-coordinates

# Disable merging entirely
impg query -a alignments.paf -r chr1:1000-2000 --no-merge

# Force processing large regions (>10kbp) with gfa/maf output
impg query -a alignments.paf -r chr1:1000-50000 -o gfa --sequence-files *.fa --force-large-region

# Sparsify wfmash mappings for GFA output (auto or explicit fraction)
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --sequence-files *.fa --sparsify auto
```

#### GFA engine selection

When outputting GFA (`-o gfa`), use `--engine` to choose the graph construction algorithm (default: `pggb`):

| Engine | Algorithm | Best for |
|--------|-----------|----------|
| `pggb` | seqwish + smoothxg-style smoothing + gfaffix normalization | Default. Smoothed, normalized variation graphs. |
| `recursive` | Recursive sweepga-guided chunking + POA per chunk + lacing | Fast, compact graphs for mid-size regions. |
| `seqwish` | All-vs-all sweepga + transitive closure graph induction | Raw (unsmoothed) variation graphs. |
| `poa` | Single-pass partial order alignment (SPOA) | Small regions, quick MSA-based output. |

```bash
# Choose a GFA engine (default: pggb)
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --engine pggb --sequence-files *.fa
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --engine recursive --sequence-files *.fa
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --engine seqwish --sequence-files *.fa
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --engine poa --sequence-files *.fa

# Disable alignment filtering for seqwish (faster but may produce less clean graphs)
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --engine seqwish -N --sequence-files *.fa

# Tune recursive engine parameters
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --engine recursive \
    --recursive-poa-threshold 1000 --recursive-chunk-size 5000 --recursive-padding 100 \
    --sequence-files *.fa

# Custom POA scoring (match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2)
impg query -a alignments.paf -r chr1:1000-2000 -o gfa --poa-scoring 5,4,6,2,24,1 --sequence-files *.fa

# Tune seqwish graph induction parameters (seqwish and pggb engines)
impg query -a alignments.paf -r chr1:1000-50000 -o gfa --engine seqwish --sequence-files *.fa \
    --min-match-len 50 --transclose-batch 5000000 --sparse-factor 0.5

# Tune smoothxg-style smoothing (pggb engine, default: two passes at 700,1100 bp)
impg query -a alignments.paf -r chr1:1000-50000 -o gfa --engine pggb --sequence-files *.fa \
    --target-poa-length 700,1100 --max-node-length 200 --poa-padding-fraction 0.001
```

#### Alignment visualizations

The `scripts/faln2html.py` tool converts FASTA alignments into interactive HTML visualizations that can be viewed in any web browser. It supports [react-msa](https://github.com/GMOD/JBrowseMSA) and [ProSeqViewer](https://github.com/BioComputingUP/ProSeqViewer) as MSA viewers.

```bash
# Visualize FASTA alignments in the browser (pipe directly to visualization script)
impg query -a alignments.paf -r chr1:1000-2000 -o fasta-aln --sequence-files *.fa | \
  python scripts/faln2html.py -i - -o alignment.html

# Choose visualization tool (reactmsa or proseqviewer)
impg query -a alignments.paf -r chr1:1000-2000 -o fasta-aln --sequence-files *.fa | \
  python scripts/faln2html.py -i - -o alignment.html --tool proseqviewer
```

### Partition

Partition the alignment into smaller pieces:

```bash
# Basic partitioning with 1Mb windows (outputs single partitions.bed file with partition number in 4th column)
impg partition -a alignments.paf -w 1000000

# Output separate files for each partition
impg partition -a alignments.1aln -w 1000000 --separate-files

# Specify output folder for partition files (directory will be created if it doesn't exist)
impg partition -a file1.paf file2.1aln -w 1000000 --output-folder results

# Start from specific sequences (one per line)
impg partition -a alignments.paf -w 1000000 --starting-sequences-file seqs.txt

# Merge nearby intervals within partitions
impg partition -a file1.paf file2.1aln -w 1000000 -d 10000

# Selection strategies for next sequence
impg partition -a alignments.paf -w 1000000 --selection-mode longest        # longest missing region
impg partition -a alignments.1aln -w 1000000 --selection-mode total          # most total missing
impg partition -a alignments.paf -w 1000000 --selection-mode sample         # by sample (PanSN)
impg partition -a alignments.1aln -w 1000000 --selection-mode haplotype      # by haplotype (PanSN)

# Control transitive search depth and minimum region size
impg partition -a file1.paf file2.1aln -w 1000000 -m 2 --min-transitive-len 10000

# Approximate mode with partition (requires --min-transitive-len > trace_spacing)
impg partition -a alignments.1aln -w 1000000 --approximate --min-transitive-len 101 -o bed
# Output as GFA, MAF or FASTA requires sequence files and --separate-files flag
impg partition -a alignments.paf -w 1000000 -o gfa --sequence-files *.fa --separate-files --output-folder gfa_partitions
impg partition -a alignments.1aln -w 1000000 -o maf --sequence-list fastas.txt --separate-files --output-folder maf_partitions
impg partition -a file1.paf file2.1aln -w 1000000 -o fasta --sequence-files *.fa --separate-files --output-folder fasta_partitions
# Works with AGC archives too
impg partition -a alignments.paf -w 1000000 -o gfa --sequence-files genomes.agc --separate-files --output-folder gfa_partitions

# GFA engine selection (same engines as query: pggb, recursive, seqwish, poa; default: pggb)
impg partition -a alignments.paf -w 1000000 -o gfa --engine pggb --sequence-files *.fa --separate-files
impg partition -a alignments.paf -w 1000000 -o gfa --engine seqwish --sequence-files *.fa --separate-files
impg partition -a alignments.paf -w 1000000 -o gfa --engine seqwish -N --sequence-files *.fa --separate-files  # no filtering

# Tune seqwish graph induction parameters (seqwish and pggb engines)
impg partition -a alignments.paf -w 1000000 -o gfa --engine seqwish --sequence-files *.fa --separate-files \
    --min-match-len 50 --transclose-batch 5000000

# Tune smoothxg-style smoothing (pggb engine, default: two passes at 700,1100 bp)
impg partition -a alignments.paf -w 1000000 -o gfa --engine pggb --sequence-files *.fa --separate-files \
    --target-poa-length 700,1100 --max-node-length 200
```

### Similarity

Compute pairwise similarity between sequences in a region:

```bash
# Basic similarity computation
impg similarity -a alignments.paf -r chr1:1000-2000 --sequence-files ref.fa genomes.fa

# Query multiple regions from a BED file (it produces multiple similarity matrices)
impg similarity -a alignments.1aln -b regions.bed --sequence-files *.fa

# Output distances instead of similarities
impg similarity -a file1.paf file2.1aln -r chr1:1000-2000 --sequence-files *.fa --distances

# Include all pairs (even those with zero similarity)
impg similarity -a alignments.paf -r chr1:1000-2000 --sequence-files *.fa --all

# Restrict analysis to sequences listed in a file (one name per line)
# Entries may be full contig names or sample identifiers (e.g., HG00097 or HG00097_hap1)
impg similarity -a alignments.1aln -r chr1:1000-2000 --sequence-files *.fa --subset-sequence-list sequences.txt

# Group sequences by delimiter (e.g., for PanSN naming, "sample#haplotype#chr" -> "sample")
impg similarity -a alignments.paf -r chr1:1000-2000 --sequence-files *.fa --delim '#'

# Use 2nd occurrence of delimiter for grouping (e.g., for PanSN naming, "sample#haplotype#chr" -> "sample#haplotype")
impg similarity -a alignments.1aln -r chr1:1000-2000 --sequence-files *.fa --delim '#' --delim-pos 2

# Perform PCA/MDS dimensionality reduction
impg similarity -a file1.paf file2.1aln -r chr1:1000-2000 --sequence-files *.fa --pca

# Specify number of PCA components (default: 2)
impg similarity -a alignments.paf -r chr1:1000-2000 --sequence-files *.fa --pca --pca-components 3

# Choose similarity measure for PCA distance matrix (jaccard/cosine/dice, default: jaccard)
impg similarity -a alignments.1aln -r chr1:1000-2000 --sequence-files *.fa --pca --pca-measure cosine

# PCA with adaptive polarization using previous regions
impg similarity -a file1.paf file2.1aln -b regions.bed --sequence-files *.fa --pca --polarize-n-prev 3

# PCA with sample-guided polarization
impg similarity -a alignments.paf -b regions.bed --sequence-files *.fa --pca --polarize-guide-samples sample1,sample2

# Engine selection (poa or recursive; default: poa)
impg similarity -a alignments.paf -r chr1:1000-2000 --sequence-files *.fa --engine recursive
```

### Refine

Refine loci to maximize sample support:

```bash
# Refine a single region to maximize the number of sequences spanning both ends
impg refine -a alignments.paf -r chr1:1000-2000

# Refine many regions from a BED file
impg refine -a alignments.paf -b loci.bed

# Allow merging within 200 kb and require at least 2 kb coverage near each end
impg refine -a alignments.paf -r chr1:1000-2000 -d 200000 --span-bp 2000

# Expand up to 90% of the locus length on each side (default: 0.5)
impg refine -a alignments.paf -r chr1:1000-2000 --max-extension 0.90

# Or cap the search to an absolute flank size
impg refine -a alignments.paf -r chr1:1000-2000 --max-extension 50000

# Maximize PanSN sample or haplotype counts instead of sequence counts
impg refine -a alignments.paf -r chr1:1000-2000 --pansn-mode sample
impg refine -a alignments.paf -r chr1:1000-2000 --pansn-mode haplotype

# Capture the supporting entities in a separate BED file
impg refine -a alignments.paf -r chr1:1000-2000 --support-output refine_support.bed

# Control extension step size (default: 1000 bp)
impg refine -a alignments.paf -r chr1:1000-2000 --extension-step 500

# Exclude regions from entity counting via blacklist
impg refine -a alignments.paf -r chr1:1000-2000 --blacklist-bed excluded.bed

# Restrict to specific sequences (filters transitive steps too)
impg refine -a alignments.paf -r chr1:1000-2000 --subset-sequence-list samples.txt

# Works with .1aln files too (requires --sequence-files)
impg refine -a alignments.1aln --sequence-files sequences.fa -r chr1:1000-2000

# Fast approximate mode for .1aln files (requires --min-transitive-len > trace_spacing)
impg refine -a alignments.1aln -r chr1:1000-2000 --approximate --min-transitive-len 101
```

When `--support-output` is provided, the tool emits a BED file listing every sequence/sample/haplotype that spans the refined region: `sequence	start	end	region-name`.

`impg refine` explores asymmetric left/right expansions around each target region to find the smallest window that maximizes the number of sequences, samples, or haplotypes. Keeping start/end alignment anchors outside structural variants helps avoid selecting loci that terminate inside large insertions or deletions.

### Stats

Print alignment statistics:

```bash
# Get statistics for PAF file
impg stats -a alignments.paf

# Get statistics for .1aln file
impg stats -a alignments.1aln

# Get statistics for combined PAF and .1aln files
impg stats -a file1.paf file2.1aln
```

### Lace

Combine multiple GFA or VCF files:

```bash
# Combine multiple GFA files (auto-detects format)
impg lace -f file1.gfa file2.gfa file3.gfa -o combined.gfa

# Combine multiple VCF files
impg lace -f file1.vcf file2.vcf file3.vcf -o combined.vcf

# Use a list file containing file paths
impg lace -l files.txt -o combined.gfa

# Explicitly specify input format (gfa, vcf, auto)
impg lace -f *.gfa -o combined.gfa --format gfa

# Fill gaps between contiguous path segments (GFA only)
impg lace -f *.gfa -o combined.gfa --fill-gaps 1 # Fill with N's
impg lace -f *.gfa -o combined.gfa --fill-gaps 1 --sequence-files sequence.fa # Fill with sequences

# Fill all gaps, including start and end gaps (GFA only, requires sequence files)
impg lace -f *.gfa -o combined.gfa --fill-gaps 2 --sequence-files sequence.fa

# Control output compression
impg lace -f *.gfa -o combined.gfa.gz --compress gzip
impg lace -f *.gfa -o combined.gfa.bgz --compress bgzip
impg lace -f *.gfa -o combined.gfa.zst --compress zstd

# Use reference for VCF contig validation
impg lace -f *.vcf -o combined.vcf --reference reference.fa

# Use custom temporary directory
impg lace -f *.gfa -o combined.gfa --temp-dir /tmp/lace_work
```

#### Path Name Format

The command expects path names in the format:

```
NAME:START-END
```

Example: `HG002#1#chr20:1000-2000`

The command uses these coordinates to:
1. Identify which sequences belong together
2. Order the sequences correctly
3. Detect and handle overlaps or gaps

Note: `NAME` can contain ':' characters. When parsing coordinates, the command uses the last occurrence of ':' to separate the name from the coordinate range.

#### Post-processing recommendations

After combining the GFA files, the resulting graph will already have compacted node IDs ranging from `1` to the total number of nodes. However, it is strongly recommended to perform post-processing steps using **[ODGI](https://github.com/pangenome/odgi)** to unchop and sort the graph.

```bash
odgi unchop -i combined.gfa -o - -t 16 | \
    odgi sort -i - -o - -p gYs -t 16 | \
    odgi view -i - -g > combined.final.gfa
```

If overlaps were present, and then trimmed during the merging process, it's advisable to run **[GFAffix](https://github.com/marschall-lab/GFAffix)** before the ODGI pipeline to remove redundant nodes introduced by the overlap trimming.

```bash
gfaffix combined.gfa -o combined.fix.gfa &> /dev/null

odgi unchop -i combined.fix.gfa -o - -t 16 | \
    odgi sort -i - -o - -p gYs -t 16 | \
    odgi view -i - -g > combined.final.gfa
```

### Graph

Build a pangenome graph directly from FASTA sequences (no pre-computed alignments needed):

```bash
# Build graph from FASTA files (default engine: seqwish)
impg graph --fasta-files sequences.fa -g output.gfa -t 16

# Build from multiple FASTA files
impg graph --fasta-files file1.fa file2.fa file3.fa -g output.gfa

# Use a list file containing FASTA paths
impg graph --fasta-list fasta_files.txt -g output.gfa

# Write to stdout
impg graph --fasta-files sequences.fa -g - | odgi build -g - -o output.og
```

#### Engine selection

Use `--engine` to choose the graph construction algorithm (default: `pggb`):

```bash
# Pggb: seqwish + smoothxg-style smoothing + gfaffix normalization (default)
impg graph --fasta-files sequences.fa -g output.gfa --engine pggb

# Seqwish: sweepga alignment + transitive closure graph induction (raw, unsmoothed)
impg graph --fasta-files sequences.fa -g output.gfa --engine seqwish

# Recursive: sweepga-guided chunking + POA per chunk + lacing
impg graph --fasta-files sequences.fa -g output.gfa --engine recursive

# POA: single-pass partial order alignment (fastest for small inputs)
impg graph --fasta-files sequences.fa -g output.gfa --engine poa
```

All engines produce sorted, unchopped GFA with consistent path names.

#### Alignment backend options

The seqwish, pggb, and recursive engines run all-vs-all alignment internally. The aligner can be selected with `--aligner` (default: `wfmash`; alternative: `fastga`). PanSN-formatted sequence names are handled automatically, and the k-mer frequency is scaled by the number of genomes.

> **Aligner compatibility:**
> - `--sparsify` is **wfmash-only** — combining `--aligner fastga` with `--sparsify` is an error.
> - `--frequency` / `-f` (k-mer frequency) is **fastga-only** — combining `--aligner wfmash` with `--frequency` is an error.

```bash
# Choose aligner backend (default: wfmash)
impg graph --fasta-files sequences.fa -g output.gfa --aligner fastga

# Adjust k-mer frequency multiplier for fastga (default: 10x number of genomes)
impg graph --fasta-files sequences.fa -g output.gfa --aligner fastga -f 5

# Set explicit k-mer frequency for fastga
impg graph --fasta-files sequences.fa -g output.gfa --aligner fastga --frequency 100

# Filter alignments by minimum length
impg graph --fasta-files sequences.fa -g output.gfa --min-alignment-length 500

# Sparsify wfmash mappings (auto or explicit fraction; wfmash only)
impg graph --fasta-files sequences.fa -g output.gfa --sparsify auto

# Skip pre-computed alignment step with a PAF file
impg graph --fasta-files sequences.fa -g output.gfa --paf-file alignments.paf
```

#### Seqwish graph induction options

These options control the transitive closure step used by the `seqwish` and `pggb` engines. They are available in `graph`, `query -o gfa`, and `partition -o gfa`.

```bash
# Minimum match length for alignments (default: 23)
impg graph --fasta-files sequences.fa -g output.gfa --min-match-len 50

# Batch size for transitive closure (default: 10000000; reduce for lower memory)
impg graph --fasta-files sequences.fa -g output.gfa --transclose-batch 5000000

# Use sparse factor to reduce alignment density (0.0 = keep all; default: 0.0)
impg graph --fasta-files sequences.fa -g output.gfa --sparse-factor 0.5

# Use disk-backed interval trees for very large datasets (slower but lower memory)
impg graph --fasta-files *.fa -g output.gfa --disk-backed

# Repeat filtering (suppress paths through repeats that appear too many times)
impg graph --fasta-files sequences.fa -g output.gfa --repeat-max 1000 --min-repeat-dist 500
```

#### Smoothxg-style smoothing options (pggb engine)

These options control the per-block POA smoothing step in the `pggb` engine. They are available in `graph`, `query -o gfa`, and `partition -o gfa`.

The `--target-poa-length` parameter accepts a comma-separated list of values, one per smoothing pass (matching pggb's `-G` flag). The default `700,1100` runs two passes: first resolving variation up to ~700 bp, then a second pass resolving up to ~1100 bp. Each pass feeds its output into the next, progressively smoothing longer variation.

```bash
# Two-pass smoothing with custom lengths (default: "700,1100")
impg graph --fasta-files sequences.fa -g output.gfa --engine pggb --target-poa-length 700,1100

# Single-pass smoothing
impg graph --fasta-files sequences.fa -g output.gfa --engine pggb --target-poa-length 700

# Three-pass smoothing for very diverse regions
impg graph --fasta-files sequences.fa -g output.gfa --engine pggb --target-poa-length 700,1100,1500

# Maximum node length before chopping (default: 100)
impg graph --fasta-files sequences.fa -g output.gfa --engine pggb --max-node-length 200

# POA padding fraction of average block length (default: 0.001)
impg graph --fasta-files sequences.fa -g output.gfa --engine pggb --poa-padding-fraction 0.01
```

#### Alignment filtering options (seqwish engine)

By default, `impg graph` applies sweepga's plane-sweep filtering to produce clean alignments with scaffold-based chaining. This removes spurious cross-chromosome alignments (e.g., from repetitive elements like TEs, telomeres, or rDNA).

```bash
# Disable filtering entirely (much faster, but may produce less clean graphs)
impg graph --fasta-files sequences.fa -g output.gfa -N

# Control mapping cardinality (default: many:many)
impg graph --fasta-files sequences.fa -g output.gfa --num-mappings 1:1    # strictest
impg graph --fasta-files sequences.fa -g output.gfa --num-mappings n:n    # most permissive

# Scaffold-based filtering (chains alignments along diagonal)
impg graph --fasta-files sequences.fa -g output.gfa --scaffold-jump 50000   # Max gap between chained alignments
impg graph --fasta-files sequences.fa -g output.gfa --scaffold-mass 10000   # Min total aligned bases in scaffold
impg graph --fasta-files sequences.fa -g output.gfa --scaffold-filter 1:1   # Scaffold cardinality

# Overlap and identity thresholds
impg graph --fasta-files sequences.fa -g output.gfa --overlap 0.95      # Max overlap between alignments (0.0-1.0)
impg graph --fasta-files sequences.fa -g output.gfa --min-identity 0.9  # Min alignment identity (0.0-1.0)
```

#### Recursive and POA engine options

```bash
# Tune recursive engine parameters
impg graph --fasta-files sequences.fa -g output.gfa --engine recursive \
    --recursive-poa-threshold 1000 --recursive-chunk-size 5000

# Custom POA scoring (match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2)
impg graph --fasta-files sequences.fa -g output.gfa --engine poa --poa-scoring 5,4,6,2,24,1
```

#### Temporary files

The aligner writes large intermediate files during graph construction. By default these go to `$TMPDIR` or the current working directory.

```bash
# Use a specific directory for temp files
impg graph --fasta-files sequences.fa -g output.gfa --temp-dir /scratch/tmp

# Use RAM-backed storage for faster I/O (requires enough free RAM)
impg graph --fasta-files sequences.fa -g output.gfa --temp-dir ramdisk
```

The `ramdisk` shortcut maps to `/dev/shm` on Linux. Use it only when you have sufficient free memory, as temp files can be several GB for large inputs.

#### query vs graph

Both `query -o gfa` and `graph` share the same engine implementations and parameter sets. The difference is only in how sequences are obtained:

- **`query -o gfa`**: extracts subsequences from an IMPG index + sequence files for a target region, then builds a graph.
- **`graph`**: reads FASTA files directly and builds a graph (no alignments or index needed).

All seqwish graph induction options (`--min-match-len`, `--transclose-batch`, `--sparse-factor`, `--disk-backed`, `--repeat-max`, `--min-repeat-dist`) and smoothxg smoothing options (`--target-poa-length`, `--max-node-length`, `--poa-padding-fraction`) are available in both `graph` and `query -o gfa` (and `partition -o gfa`).

### Align

Generate pairwise alignment jobs with sparsification strategies for large cohorts:

```bash
# Default: giant-component sparsification (99% edge probability)
impg align --fasta-files sequences.fa -o alignments -t 16

# No sparsification (all-vs-all)
impg align --fasta-files sequences.fa -o alignments -p none

# Random sparsification (keep 50% of pairs)
impg align --fasta-files sequences.fa -o alignments -p random:0.5

# Tree-based sparsification
impg align --fasta-files sequences.fa -o alignments -p tree:3:1:0.1

# Output as PAF or 1ALN instead of joblists
impg align --fasta-files sequences.fa -o alignments --format paf
impg align --fasta-files sequences.fa -o alignments --format 1aln

# Choose aligner backend (default: wfmash)
impg align --fasta-files sequences.fa -o alignments --aligner fastga

# fastga-specific: adjust k-mer frequency multiplier (ignored by wfmash)
impg align --fasta-files sequences.fa -o alignments --aligner fastga -f 5
impg align --fasta-files sequences.fa -o alignments --aligner fastga --frequency 100

# With alignment filtering
impg align --fasta-files sequences.fa -o alignments --num-mappings 1:1 --scaffold-filter 1:1
```

### Index

Create an IMPG index from alignment files:

```bash
# Index a single PAF file
impg index -a alignments.paf

# Index a single .1aln file
impg index -a alignments.1aln

# Index multiple alignment files (PAF and .1aln mixed)
impg index -a file1.paf file2.1aln file3.paf

# Create index with custom name
impg index -a alignments.paf -i custom.impg

# Index from a list of alignment files (can mix formats)
impg index --alignment-list alignment_files.txt
```

#### Indexing Modes

Use `--index-mode` to control how indices are built (`auto` by default):

| Mode | Description |
|------|-------------|
| `auto` | Single index when < 100 files, per-file when >= 100 |
| `single` | Always create a single combined `.impg` index |
| `per-file` | One `.impg` per alignment file (e.g., `data.paf.impg`) |

```bash
# Single combined index (explicit or via -i)
impg index -a file1.paf file2.1aln -i combined.impg
impg query -i combined.impg -r chr1:0-1000

# Per-file index
impg index --alignment-list files.txt --index-mode per-file -t 32
impg query --alignment-list files.txt --index-mode per-file -r chr1:0-1000
```

Both modes work with PAF, 1ALN, and TPA files (can be mixed in `--alignment-list`).

**When to use per-file indexing:**
- Incremental updates (only rebuild changed alignment files)
- Many alignment files (auto mode switches at 100 files)

**Stale index detection:** impg warns if alignment files are modified after index creation. Use `-f/--force-reindex` to rebuild.

**Note on compressed files**: `impg` works directly with bgzip-compressed PAF files (`.paf.gz`, `.paf.bgz`). For large files, creating a GZI index can speed up initial index creation:

```bash
bgzip -r alignments.paf.gz  # Creates alignments.paf.gz.gzi (optional)
```

### Common options

All commands support these options:
- `-a, --alignment-files`: One or more paths to alignment files in PAF, 1ALN, or TPA format (can be mixed). Files can be gzipped or uncompressed.
- `--alignment-list`: Path to a plain-text file listing one alignment path per line (PAF, 1ALN, or TPA files can be mixed).
- `-i, --index`: Path to an existing IMPG index file.
- `-f, --force-reindex`: Always regenerate the IMPG index even if it already exists.
- `-t, --threads`: Number of threads (default: 4)
- `--unidirectional`: Disable bidirectional alignment interpretation.
- `-v, --verbose`: Verbosity level (0=error/silent, 1=info with progress bar, 2=debug)

### Sequence file options

For GFA/MAF/FASTA output and similarity computation:

- `--sequence-files`: List of sequence files (FASTA or AGC)
- `--sequence-list`: Text file listing sequence files (FASTA or AGC) (one per line)
- `--poa-scoring`: POA scoring parameters as `match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2` (default: `5,4,6,2,24,1`)
- `--reverse-complement`: Reverse complement sequences on the reverse strand (for FASTA output)

### Merging behaviour

- `-d, --merge-distance <INT>`: Merge nearby hits within this distance (bp).
- `--no-merge`: Disable merging entirely for all output formats.
- `--consider-strandness`: Keep forward and reverse strands separate when merging. By default, strands are merged for BED/GFA/MAF outputs and kept separate for FASTA/FASTA-ALN.

## What does `impg` do?

At its core, `impg` lifts over ranges from a target sequence (used as reference) into the queries (the other sequences aligned to the sequence used as reference) described in alignments.
In effect, it lets us pick up homologous loci from all genomes mapped onto our specific target region.
This is particularly useful when you're interested in comparing a specific genomic region across different individuals, strains, or species in a pangenomic or comparative genomic setting.
The output is provided in BED, BEDPE and PAF formats, making it straightforward to use to extract FASTA sequences for downstream use in multiple sequence alignment (like `mafft`) or pangenome graph building (e.g., `pggb` or `minigraph-cactus`).

## How does it work?

`impg` uses [`coitrees`](https://github.com/dcjones/coitrees) (Cache Oblivious Interval Trees) to provide efficient range lookup over the input alignments.
CIGAR strings are converted to a compact delta encoding.
This approach allows for fast and memory-efficient projection of sequence ranges through alignments.

## Tutorial: Building a Yeast Pangenome Graph

This tutorial walks through building a complete pangenome graph from 7 *S. cerevisiae* strains using `impg`'s partition-based workflow. The approach partitions the pangenome into manageable pieces, builds graphs for each partition using sweepga + seqwish, and laces them together into a final graph.

### Prerequisites

- `impg` compiled with all dependencies
- `odgi` for visualization
- Yeast pangenome FASTA (PanSN naming: `SAMPLE#HAPLOTYPE#CONTIG`)
- Pre-computed alignments in PAF format (e.g., from `wfmash`)

### Step 1: Create the IMPG Index

Build an index from your alignment file:

```bash
# Create working directory
mkdir -p yeast_pangenome && cd yeast_pangenome

# Build the index from alignments
impg index -a cerevisiae.paf -i yeast.impg -t 16
```

### Step 2: Partition the Pangenome

Divide the pangenome into ~100kb regions:

```bash
# Partition into FASTA files for graph construction
impg partition -i yeast.impg \
    -w 100000 \
    --sequence-files cerevisiae.fa.gz \
    -o fasta \
    --separate-files \
    --output-folder partitions \
    -t 16
```

This creates one FASTA file per partition (e.g., `partitions/partition0.fasta`, `partitions/partition1.fasta`, ...).

### Step 3: Build Graphs for Each Partition

Use `impg graph` to build a GFA for each partition. By default, `impg graph` uses the `pggb` engine (seqwish + smoothing + gfaffix normalization) with 1:1 alignment filtering and scaffold-based chaining to remove spurious cross-chromosome alignments from repetitive elements:

```bash
# Create output directory
mkdir -p gfas

# Build graphs in parallel (default filtering: 1:1 + scaffolding)
ls partitions/*.fasta | xargs -P 4 -I {} bash -c '
    f="{}"; base=$(basename "$f" .fasta)
    impg graph --fasta-files "$f" -g "gfas/${base}.gfa" -t 4
    echo "Done: $base"
'
```

The k-mer frequency is automatically calculated based on the number of genomes (unique SAMPLE#HAPLOTYPE prefixes), not the total number of sequences. For 7 yeast strains, this means `-f 10` results in a frequency of 70 (10 × 7 genomes).

### Step 4: Lace Partition Graphs Together

Combine all partition GFAs into a single pangenome graph:

```bash
# Create list of valid (non-empty) GFA files
find gfas -name "*.gfa" -size +0 | sort -V > gfa_list.txt

# Lace together with gap filling
impg lace \
    --file-list gfa_list.txt \
    --sequence-files cerevisiae.fa.gz \
    -o yeast_pangenome.gfa \
    --fill-gaps 2 \
    -t 16
```

The `--fill-gaps 2` option fills gaps between partitions with the original sequence.

### Step 5: Post-process and Visualize with ODGI

Convert to ODGI format and sort for optimal visualization:

```bash
# Build ODGI graph
odgi build -g yeast_pangenome.gfa -o yeast_pangenome.og -t 16

# Sort the graph
odgi sort -i yeast_pangenome.og -o yeast_pangenome_sorted.og -O -p Ygs -t 16

# Generate 2D layout
odgi layout -i yeast_pangenome_sorted.og -o yeast_pangenome.lay -t 16 -P

# Create visualizations
# Linear view (paths as rows, colored by sample)
odgi viz -i yeast_pangenome_sorted.og -o yeast_pangenome_viz.png -x 4000 -y 1000 -s '#'

# 2D graph drawing
odgi draw -i yeast_pangenome_sorted.og -c yeast_pangenome.lay -p yeast_pangenome_draw.png -w 4000 -H 2000
```

### Step 6: Verify the Graph

Check graph statistics:

```bash
odgi stats -i yeast_pangenome_sorted.og -S
```

Expected output for 7 yeast strains with 16 chromosomes each:
- ~35 Mbp total length
- ~440,000 nodes
- ~616,000 edges
- 112 paths (7 strains × 16 chromosomes)

### Complete Pipeline Script

```bash
#!/bin/bash
set -euo pipefail

FASTA="cerevisiae.fa.gz"
PAF="cerevisiae.paf"
THREADS=16
WINDOW=100000

# Step 1: Index
impg index -a "$PAF" -i yeast.impg -t "$THREADS"

# Step 2: Partition
mkdir -p partitions gfas
impg partition -i yeast.impg -w "$WINDOW" \
    --sequence-files "$FASTA" -o fasta \
    --separate-files --output-folder partitions -t "$THREADS"

# Step 3: Build partition graphs
ls partitions/*.fasta | xargs -P 4 -I {} bash -c '
    f="{}"; base=$(basename "$f" .fasta)
    impg graph --fasta-files "$f" -g "gfas/${base}.gfa" -t 4
'

# Step 4: Lace
find gfas -name "*.gfa" -size +0 | sort -V > gfa_list.txt
impg lace --file-list gfa_list.txt --sequence-files "$FASTA" \
    -o yeast_pangenome.gfa --fill-gaps 2 -t "$THREADS"

# Step 5: ODGI post-processing
odgi build -g yeast_pangenome.gfa -o yeast_pangenome.og -t "$THREADS"
odgi sort -i yeast_pangenome.og -o yeast_pangenome_sorted.og -O -p Ygs -t "$THREADS"
odgi layout -i yeast_pangenome_sorted.og -o yeast_pangenome.lay -t "$THREADS"
odgi viz -i yeast_pangenome_sorted.og -o yeast_pangenome_viz.png -x 4000 -y 1000 -s '#'
odgi draw -i yeast_pangenome_sorted.og -c yeast_pangenome.lay -p yeast_pangenome_draw.png

echo "Done! Check yeast_pangenome_viz.png and yeast_pangenome_draw.png"
```

### Exploration: Effect of Partition Size

The partition window size (`-w`) affects the resulting graph structure. Smaller partitions create more alignment subproblems but may fragment complex structural variants, while larger partitions allow better representation of larger variations but increase memory usage per partition.

Here's a comparison using the 7-strain yeast pangenome with different partition sizes:

| Partition Size | Partitions | Non-empty GFAs | Graph Length | Nodes | Edges | Steps |
|---------------|------------|----------------|--------------|-------|-------|-------|
| 10kb | 1,537 | 1,533 | 63.0 Mb | 164,582 | 221,966 | 243,526 |
| 50kb | 368 | 367 | 59.5 Mb | 206,982 | 280,259 | 298,275 |
| 100kb | 227 | 227 | 59.1 Mb | 221,360 | 299,952 | 325,202 |

**Observations:**

- **Smaller partitions (10kb)** produce graphs with fewer nodes and edges but more total sequence length due to gap filling between many partition boundaries
- **Larger partitions (50-100kb)** create more complex graphs with better variant representation but require more memory per partition
- All produce 112 paths (7 strains x 16 chromosomes)

**Recommendations:**

- Use **10kb partitions** for small genomes or when memory is constrained
- Use **50-100kb partitions** for better structural variant representation
- For human-scale genomes, consider **100kb-500kb partitions**

To experiment with partition sizes:

```bash
# Try different window sizes
for WINDOW in 10000 50000 100000; do
    mkdir -p partitions_${WINDOW} gfas_${WINDOW}

    impg partition -i index.impg -w $WINDOW \
        --sequence-files sequences.fa -o fasta \
        --separate-files --output-folder partitions_${WINDOW}

    ls partitions_${WINDOW}/*.fasta | xargs -P 4 -I {} bash -c '
        f="{}"; base=$(basename "$f" .fasta)
        impg graph --fasta-files "$f" -g "gfas_${WINDOW}/${base}.gfa" -t 4
    '

    find gfas_${WINDOW} -name "*.gfa" -size +0 | sort -V > gfa_list_${WINDOW}.txt
    impg lace --file-list gfa_list_${WINDOW}.txt --sequence-files sequences.fa \
        -o pangenome_${WINDOW}.gfa --fill-gaps 2
done
```

## Authors

Andrea Guarracino <aguarra1@uthsc.edu> \
Bryce Kille <brycekille@gmail.com> \
Erik Garrison <erik.garrison@gmail.com>

## License

MIT
