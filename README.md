# impg: implicit pangenome graph

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/impg/README.html)

Pangenome graphs and whole genome multiple alignments are powerful tools, but they are expensive to build and manipulate.
Often, we would like to be able to break a small piece out of a pangenome without constructing the whole thing.
`impg` lets us do this by projecting sequence ranges through many-way (e.g. all-vs-all) pairwise alignments built by tools like `wfmash` and `minimap2`.

## Usage

Here's a basic example:

```bash
impg query -p cerevisiae.pan.paf.gz -r S288C#1#chrI:50000-100000 -x
```

- `-p` specifies the path to the PAF file. Your alignments must use `wfmash` default or `minimap2 --eqx` type CIGAR strings which have `=` for matches and `X` for mismatches.
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

Alternatively, install from Bioconda:

```bash
conda install -c bioconda impg
```

## Commands

### Query

Query overlaps in the alignment:

```bash
# Query a single region
impg query -p alignments.paf -r chr1:1000-2000 

# Query multiple regions from a BED file
impg query -p alignments.paf -b regions.bed

# Enable transitive overlap search
impg query -p alignments.paf -r chr1:1000-2000 -x

# Set maximum transitive depth (default: 0 = unlimited)
impg query -p alignments.paf -r chr1:1000-2000 -x -m 3

# Filter by minimum gap-compressed identity
impg query -p alignments.paf -r chr1:1000-2000 --min-identity 0.9

# Output formats (auto/bed/bedpe/paf/gfa/maf)
impg query -p alignments.paf -r chr1:1000-2000 -o bedpe
impg query -p alignments.paf -b regions.bed -o paf

# GFA/MAF/FASTA output requires FASTA files (--fasta-files or --fasta-list)
impg query -p alignments.paf -r chr1:1000-2000 -o gfa --fasta-files ref.fa genomes.fa
impg query -p alignments.paf -r chr1:1000-2000 -o maf --fasta-list fastas.txt
impg query -p alignments.paf -r chr1:1000-2000 -o fasta --fasta-files *.fa

# FASTA output with reverse complement for reverse strand sequences
impg query -p alignments.paf -r chr1:1000-2000 -o fasta --fasta-files *.fa --reverse-complement

# Merge nearby regions (default: 0)
impg query -p alignments.paf -r chr1:1000-2000 -d 1000

# Use DFS for transitive search (slower but fewer overlapping results)
impg query -p alignments.paf -r chr1:1000-2000 --transitive-dfs
```

### Partition

Partition the alignment into smaller pieces:

```bash
# Basic partitioning with 1Mb windows
impg partition -p alignments.paf -w 1000000

# Specify output folder for partition files (directory will be created if it doesn't exist)
impg partition -p alignments.paf -w 1000000 --output-folder results

# Start from specific sequences (one per line)
impg partition -p alignments.paf -w 1000000 --starting-sequences-file seqs.txt

# Merge nearby intervals within partitions
impg partition -p alignments.paf -w 1000000 -d 10000

# Selection strategies for next sequence
impg partition -p alignments.paf -w 1000000 --selection-mode longest        # longest missing region
impg partition -p alignments.paf -w 1000000 --selection-mode total          # most total missing
impg partition -p alignments.paf -w 1000000 --selection-mode sample         # by sample (PanSN)
impg partition -p alignments.paf -w 1000000 --selection-mode haplotype      # by haplotype (PanSN)

# Control transitive search depth and minimum sizes
impg partition -p alignments.paf -w 1000000 -m 2 -l 10000

# Output as GFA, MAF or FASTA requires FASTA files (--fasta-files or --fasta-list)
impg partition -p alignments.paf -w 1000000 -o gfa --fasta-files *.fa --output-folder gfa_partitions
impg partition -p alignments.paf -w 1000000 -o maf --fasta-list fastas.txt --output-folder maf_partitions
impg partition -p alignments.paf -w 1000000 -o fasta --fasta-files *.fa --output-folder fasta_partitions
```

### Similarity

Compute pairwise similarity between sequences in a region:

```bash
# Basic similarity computation
impg similarity -p alignments.paf -r chr1:1000-2000 --fasta-files ref.fa genomes.fa

# Query multiple regions from a BED file (it produces multiple similarity matrices)
impg similarity -p alignments.paf -b regions.bed --fasta-files *.fa

# Output distances instead of similarities
impg similarity -p alignments.paf -r chr1:1000-2000 --fasta-files *.fa --distances

# Include all pairs (even those with zero similarity)
impg similarity -p alignments.paf -r chr1:1000-2000 --fasta-files *.fa -a

# Group sequences by delimiter (e.g., for PanSN naming, "sample#haplotype#chr" -> "sample")
impg similarity -p alignments.paf -r chr1:1000-2000 --fasta-files *.fa --delim '#'

# Use 2nd occurrence of delimiter for grouping (e.g., for PanSN naming, "sample#haplotype#chr" -> "sample#haplotype")
impg similarity -p alignments.paf -r chr1:1000-2000 --fasta-files *.fa --delim '#' --delim-pos 2

# Perform PCA/MDS dimensionality reduction
impg similarity -p alignments.paf -r chr1:1000-2000 --fasta-files *.fa --pca

# Specify number of PCA components (default: 2)
impg similarity -p alignments.paf -r chr1:1000-2000 --fasta-files *.fa --pca --pca-components 3

# Choose similarity measure for PCA distance matrix (jaccard/cosine/dice, default: jaccard)
impg similarity -p alignments.paf -r chr1:1000-2000 --fasta-files *.fa --pca --pca-measure cosine

# PCA with adaptive polarization using previous regions
impg similarity -p alignments.paf -b regions.bed --fasta-files *.fa --pca --polarize-n-prev 3

# PCA with sample-guided polarization
impg similarity -p alignments.paf -b regions.bed --fasta-files *.fa --pca --polarize-guide-samples sample1,sample2
```

### Stats

Print alignment statistics:

```bash
impg stats -p alignments.paf
```

### Index

Create an IMPG index from PAF files:

```bash
# Index a single PAF file
impg index -p alignments.paf

# Index multiple PAF files
impg index -p file1.paf file2.paf file3.paf

# Create index with custom name
impg index -p alignments.paf -i custom.impg

# Index from a list of PAF files
impg index --paf-list paf_files.txt
```


### Common options

All commands support these options:
- `-p, --paf-files`: One or more paths to PAF files (gzipped or uncompressed).
- `--paf-list`: Path to a plain-text file listing one PAF path per line.
- `-i, --index`: Path to an existing IMPG index file.
- `-f, --force-reindex`: Always regenerate the IMPG index even if it already exists.
- `-t, --threads`: Number of threads (default: 4)
- `-v, --verbose`: Verbosity level (0=error, 1=info, 2=debug)

### FASTA options

For GFA/MAF/FASTA output and similarity computation:

- `--fasta-files`: List of FASTA files
- `--fasta-list`: Text file listing FASTA files (one per line)
- `--poa-scoring`: POA scoring parameters as `match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2` (default: `1,4,6,2,26,1`)
- `--reverse-complement`: Reverse complement sequences on the reverse strand (for FASTA output)

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
