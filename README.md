# impg: implicit pangenome graph

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/impg/README.html)

Pangenome graphs and whole genome multiple alignments are powerful tools, but they are expensive to build and manipulate.
Often, we would like to be able to break a small piece out of a pangenome without constructing the whole thing.
`impg` lets us do this by projecting sequence ranges through many-way (e.g. all-vs-all) pairwise alignments built by tools like `wfmash` and `minimap2`.

## Using `impg`

Here's a basic example:

```bash
impg query -p cerevisiae.pan.paf.gz -r S288C#1#chrI:50000-100000 -x
```

- `-p` specifies the path to the PAF file. Your alignments must use `wfmash` default or `minimap2 --eqx` type CIGAR strings which have `=` for matches and `X` for mismatches. The `M` positional match character is not allowed.
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

You need Rust (cargo) installed. Then:

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

# Output formats (auto/bed/bedpe/paf)
impg query -p alignments.paf -r chr1:1000-2000 -o bedpe
impg query -p alignments.paf -b regions.bed -o paf

# Merge nearby regions (default: 0)
impg query -p alignments.paf -r chr1:1000-2000 -d 1000
```

### Partition

Partition the alignment into smaller pieces:

```bash
# Basic partitioning with 1Mb windows
impg partition -p alignments.paf -w 1000000

# Start from specific sequences (one per line)
impg partition -p alignments.paf -w 1000000 --starting-sequences-file seqs.txt

# Merge nearby intervals within partitions
impg partition -p alignments.paf -w 1000000 -d 10000

# Selection strategies for next sequence
impg partition -p alignments.paf -w 1000000 --selection-mode longest        # longest missing region
impg partition -p alignments.paf -w 1000000 --selection-mode total          # most total missing
impg partition -p alignments.paf -w 1000000 --selection-mode sample         # by sample (PanSN)
impg partition -p alignments.paf -w 1000000 --selection-mode haplotype,#    # by haplotype

# Control transitive search depth and minimum sizes
impg partition -p alignments.paf -w 1000000 -m 2 -l 10000
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


### Common Options

All commands support these options:
- `-p, --paf-files`: One or more paths to PAF files (gzipped or uncompressed).
- `--paf-list`: Path to a plain-text file listing one PAF path per line.
- `-i, --index`: Path to an existing IMPG index file.
- `-f, --force-reindex`: Always regenerate the IMPG index even if it already exists.
- `-t, --num-threads`: Number of threads (default: 4)
- `-v, --verbose`: Verbosity level (0=error, 1=info, 2=debug)

## What does `impg` do?

At its core, `impg` lifts over ranges from a target sequence (used as reference) into the queries (the other sequences aligned to the sequence used as reference) described in alignments.
In effect, it lets us pick up homologous loci from all genomes mapped onto our specific target region.
This is particularly useful when you're interested in comparing a specific genomic region across different individuals, strains, or species in a pangenomic or comparative genomic setting.
The output is provided in BED, BEDPE and PAF formats, making it straightforward to use to extract FASTA sequences for downstream use in multiple sequence alignment (like `mafft`) or pangenome graph building (e.g., `pggb` or `minigraph-cactus`).

## How does it work?

`impg` uses `coitrees` (implicit interval trees) to provide efficient range lookup over the input alignments.
CIGAR strings are converted to a compact delta encoding.
This approach allows for fast and memory-efficient projection of sequence ranges through alignments.

## Authors

Erik Garrison <erik.garrison@gmail.com>
Andrea Guarracino <aguarra1@uthsc.edu>
Bryce Kille <brycekille@gmail.com>

## License

MIT
