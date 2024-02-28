# impg: implicit pangenome graph

Pangenome graphs and whole genome multiple alignments are powerful tools, but they are expensive to build and manipulate.
Often, we would like to be able to break a small piece out of a pangenome without constructing the whole thing.
`impg` lets us do this by projecting sequence ranges through many-way (e.g. all-vs-all) pairwise alignments built by tools like `wfmash` and `minimap2`.

## What does `impg` do?

At its core, `impg` lifts over ranges from a target sequence into the other genomes described in alignments.
In effect, it lets us pick up homologous loci from all genomes mapped onto our specific target region.
This is particularly useful when you're interested in comparing a specific genomic region across different individuals, strains, or species in a pangenomic or comparative genomic setting.
The output is provided in BED format, making it straightforward to use to extract FASTA sequences for downstream use in multiple sequence alignment (like `mafft`) or pangenome graph building (e.g., `pggb` or `minigraph-cactus`).

## How does it work?

`impg` uses coitrees (implicit interval trees) to provide efficient range lookup over the input alignments.
CIGAR strings are converted to a compact delta encoding.
This approach allows for fast and memory-efficient projection of sequence ranges through alignments.

## Using `impg`

Getting started with `impg` is straightforward. Here's a basic example of how to use the command-line utility:

```bash
impg -p cerevisiae.pan.paf.gz -q S288C#1#chrI:50000-100000 -x
```

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

In this example, `-p` specifies the path to the PAF file (compressed with zstd in this case), `-q` defines the query in the format of `target:start-end`, and `-x` requests a *transitive closure* of the matches.
That is, for each collected range, we then find what sequence ranges are aligned onto it.
This is done progressively until we've closed the set of alignments connected to the initial query range.

### Installation

To compile and install `impg` from source, you'll need a recent rust build toolchain and cargo.

1. Clone the repository:
   ```bash
   git clone https://github.com/ekg/impg.git
   ```
2. Navigate to the `impg` directory:
   ```bash
   cd impg
   ```
3. Compile the tool (requires rust build tools):
   ```bash
   cargo install --force --path .
   ```

## Author

Erik Garrison <erik.garrison@gmail.com>

## License

MIT
