# impg: implicit pangenome graph

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/impg/README.html)

Pangenome graphs and whole genome multiple alignments are powerful tools, but they are expensive to build and manipulate.
Often, we would like to be able to break a small piece out of a pangenome without constructing the whole thing.
`impg` lets us do this by projecting sequence ranges through many-way (e.g. all-vs-all) pairwise alignments built by tools like `wfmash` and `minimap2`.

## Using `impg`

Getting started with `impg` is straightforward. Here's a basic example of how to use the command-line utility:

```bash
impg query -p cerevisiae.pan.paf.gz -r S288C#1#chrI:50000-100000 -x
```

Your alignments must use `wfmash` default or `minimap2 --eqx` type CIGAR strings which have `=` for matches and `X` for mismatches. The `M` positional match character is not allowed.

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

In this example, `-p` specifies the path to the PAF file, `-r` defines the target range in the format of `seq_name:start-end`, and `-x` requests a *transitive closure* of the matches.
That is, for each collected range, we then find what sequence ranges are aligned onto it.
This is done progressively until we've closed the set of alignments connected to the initial target range.

## Installation

To compile and install `impg` from source, you'll need a recent rust build toolchain and cargo.

1. Clone the repository:
   ```bash
   git clone https://github.com/pangenome/impg.git
   ```
2. Navigate to the `impg` directory:
   ```bash
   cd impg
   ```
3. Compile the tool:
   ```bash
   cargo install --force --path .
   ```
## Commands

`impg` provides three main commands:

### Query
Query overlaps in the alignment:
```bash
# Query a single region
impg query -p alignments.paf -r chr1:1000-2000 

# Query multiple regions from a BED file
impg query -p alignments.paf -b regions.bed

# Enable transitive overlap search
impg query -p alignments.paf -r chr1:1000-2000 -x

# Output in PAF format
impg query -p alignments.paf -r chr1:1000-2000 -P
```

### Partition
Partition the alignment into smaller pieces:
```bash
impg partition -p alignments.paf -w 1000000 -s chr1 -d 10000 -l 5000
```
- `-w`: Window size for partitioning
- `-s`: Prefix of sequence names to start partitioning from
- `-d`: Maximum distance to merge intervals in each partition
- `-l`: Minimum length for intervals in each partition (this can lead to overlapping partitions)

### Stats
Print alignment statistics:
```bash
impg stats -p alignments.paf
```

### Common Options

All commands support these options:
- `-p, --paf-file`: Path to PAF file (gzipped or uncompressed)
- `-t, --num-threads`: Number of threads (default: 1)
- `-I, --force-reindex`: Force regeneration of index
- `-v, --verbose`: Verbosity level (0=error, 1=info, 2=debug)

## What does `impg` do?

At its core, `impg` lifts over ranges from a target sequence (used as reference) into the queries (the other sequences aligned to the sequence used as reference) described in alignments.
In effect, it lets us pick up homologous loci from all genomes mapped onto our specific target region.
This is particularly useful when you're interested in comparing a specific genomic region across different individuals, strains, or species in a pangenomic or comparative genomic setting.
The output is provided in BED, BEDPE and PAF formats, making it straightforward to use to extract FASTA sequences for downstream use in multiple sequence alignment (like `mafft`) or pangenome graph building (e.g., `pggb` or `minigraph-cactus`).

## How does it work?

`impg` uses coitrees (implicit interval trees) to provide efficient range lookup over the input alignments.
CIGAR strings are converted to a compact delta encoding.
This approach allows for fast and memory-efficient projection of sequence ranges through alignments.

## Authors

Erik Garrison <erik.garrison@gmail.com>
Andrea Guarracino <aguarra1@uthsc.edu>
Bryce Kille <brycekille@gmail.com>

## License

MIT
