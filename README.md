# impg: implicit pangenome graph

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/impg/README.html)

## Why impg?

Studying genomic variation at specific loci—disease genes, regulatory elements, structural variants—across populations or species traditionally requires either building expensive whole-genome graphs or using reference-based methods that miss variation. `impg` solves this by treating all-vs-all pairwise alignments as an *implicit pangenome graph*, rapidly projecting target ranges through the alignment network to extract only the homologous sequences you need. Query regions across many genomes in seconds. Perform transitive searches to discover connected sequences. Partition genomes into comparable loci. Refine regions to maximize sample coverage—all without constructing explicit graph structure. This makes pangenome-scale comparative genomics fast and practical.

## Usage

Here's a basic example:

```bash
impg query -a cerevisiae.pan.paf.gz -r S288C#1#chrI:50000-100000 -x
```

- `-a` specifies the path to the alignment file in PAF or .1aln format. PAF files must use CIGAR strings with `=` for matches and `X` for mismatches (e.g., from `wfmash` or `minimap2 --eqx`).
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

To install without AGC support (using only FASTA files):

```bash
git clone https://github.com/pangenome/impg.git
cd impg
cargo install --force --path . --no-default-features
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

# Set maximum transitive depth (default: 0 = unlimited)
impg query -a alignments.1aln -r chr1:1000-2000 -x -m 3

# Filter by minimum gap-compressed identity
impg query -a alignments.paf alignments.1aln -r chr1:1000-2000 --min-identity 0.9

# Output formats (auto/bed/bedpe/paf/gfa/maf/fasta/fasta+paf/fasta-aln)
impg query -a alignments.paf -r chr1:1000-2000 -o bed
impg query -a alignments.1aln -r chr1:1000-2000 -o bedpe
impg query -a file1.paf file2.1aln -b chr1:1000-2000 -o paf

# Write output to file instead of stdout (using -O / --output-prefix)
impg query -p alignments.paf -r chr1:1000-2000 -o bed -O results       # creates results.bed

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

# Use DFS for transitive search (slower but fewer overlapping results)
impg query -a alignments.paf -r chr1:1000-2000 --transitive-dfs

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
impg similarity -a alignments.paf -r chr1:1000-2000 --sequence-files *.fa -a

# Restrict analysis to sequences listed in a file (one name per line)
# Entries may be full contig names or sample identifiers (e.g., HG00097 or HG00097_hap1)
impg similarity -a alignments.1aln -r chr1:1000-2000 --sequence-files *.fa --subset-sequence-list sequences.txt
# Show the progress bar
impg similarity -a file1.paf file2.1aln -r chr1:1000-2000 --sequence-files *.fa --progress-bar

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

**Combined index** (default): Creates a single `.impg` file for all alignments.
```bash
impg index -a file1.paf file2.1aln -i combined.impg
impg query -i combined.impg -r chr1:0-1000
```

**Per-file index** (`--per-file-index`): Creates one `.impg` per alignment file (e.g., `data.paf.impg`).
```bash
impg index --alignment-list files.txt --per-file-index -t 32
impg query --alignment-list files.txt --per-file-index -r chr1:0-1000
```

Both modes work with PAF and .1aln files (can be mixed in `--alignment-list`).

**When to use per-file indexing:**
- Incremental updates (only rebuild changed alignment files)
- Many alignment files

**Stale index detection:** impg warns if alignment files are modified after index creation. Use `-f/--force-reindex` to rebuild.

**Note on compressed files**: `impg` works directly with bgzip-compressed PAF files (`.paf.gz`, `.paf.bgz`). For large files, creating a GZI index can speed up initial index creation:

```bash
bgzip -r alignments.paf.gz  # Creates alignments.paf.gz.gzi (optional)
```

### Common options

All commands support these options:
- `-a, --alignment-files`: One or more paths to alignment files in PAF or .1aln format (can be mixed). Files can be gzipped or uncompressed.
- `--alignment-list`: Path to a plain-text file listing one alignment path per line (PAF or .1aln files can be mixed).
- `-i, --index`: Path to an existing IMPG index file.
- `-f, --force-reindex`: Always regenerate the IMPG index even if it already exists.
- `-t, --threads`: Number of threads (default: 4)
- `-v, --verbose`: Verbosity level (0=error, 1=info, 2=debug)

### Sequence file options

For GFA/MAF/FASTA output and similarity computation:

- `--sequence-files`: List of sequence files (FASTA or AGC*)
- `--sequence-list`: Text file listing sequence files (FASTA or AGC*) (one per line)
- `--poa-scoring`: POA scoring parameters as `match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2` (default: `1,4,6,2,26,1`)
- `--reverse-complement`: Reverse complement sequences on the reverse strand (for FASTA output)

*AGC files are only supported in the full installation (default features). For FASTA-only support, install with `--no-default-features`.

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

## Authors

Andrea Guarracino <aguarra1@uthsc.edu> \
Bryce Kille <brycekille@gmail.com> \
Erik Garrison <erik.garrison@gmail.com>

## License

MIT
