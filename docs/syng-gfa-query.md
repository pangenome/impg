# Querying a Syng Index to GFA

This is the short path for making a local variation graph from an existing
syng index.

Use `impg query -a <syng-prefix> -o gfa`. The `-a` argument is the syng index
prefix, not a GFA file, and `-o gfa` selects GFA output. There is no
`impg query -g gfa`; `--gfa-engine` only chooses the GFA construction engine.
`impg` uses the syng GBWT/syncmer index to find homologous intervals, fetches
the matching sequences, and then builds a local GFA with the selected GFA
engine.

## Inputs

You need:

```text
panel.syng.1gbwt
panel.syng.1khash
panel.syng.names
panel.syng.spos
panel.syng.pstep
panel.syng.meta
panel.fa or panel.agc
```

The sequence file is required for `-o gfa`, because the query finds intervals
from syng but the GFA engine needs the actual DNA strings. Sequence names in
the FASTA/AGC must match the path names stored in the syng index.

## Single Region

```bash
impg query \
  -a panel.syng \
  -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k \
  -x \
  -o gfa \
  --sequence-files panel.agc \
  --force-large-region \
  -O c4.query
```

This writes `c4.query.gfa`.

Important flags:

- `-d 50k` is required. It merges query-gathered ranges separated by at most
  this distance and sets the largest one-hop gap/SV absorbed into one interval.
- `-x` enables transitive query. For duplicated loci, this is usually what you
  want.
- `--sequence-files` can be FASTA or AGC.
- `--force-large-region` is needed for large GFA/MAF requests.
- `-O` is an output prefix; the `.gfa` suffix is added by `impg`.

## BED Batch

Batching is preferred when querying many regions because the syng index is
loaded once.

```bash
impg query \
  -a panel.syng \
  -b regions.bed \
  -d 50k \
  -x \
  -o gfa \
  --sequence-files panel.agc \
  --force-large-region \
  -O regions
```

BED rows use the syng path name in column 1:

```text
GRCh38#0#chr6	31982056	32035418	C4
```

If a path name itself contains `:`, `-r` still works because `impg` splits
target ranges on the last colon:

```bash
impg query \
  -a c4.syng \
  -r 'grch38#chr6:31972046-32055647:0-10000' \
  -d 150 \
  -o gfa \
  --sequence-files chr6.C4.fa \
  -O c4.10kb
```

## GFA Engine

The default engine is `pggb`. You can choose another engine with
`--gfa-engine`:

```bash
# Raw seqwish graph
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o gfa --gfa-engine seqwish \
  --sequence-files panel.agc --force-large-region -O c4.seqwish

# Partitioned graph build, useful for broad regions
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o gfa --gfa-engine pggb:10000 \
  --sequence-files panel.agc --force-large-region -O c4.partitioned

# Syng-native blunt graph (default when selecting syng)
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o gfa --gfa-engine syng \
  --sequence-files panel.agc --force-large-region -O c4.syng-blunt

# Same thing with explicit mode and syncmer-parameter assertion
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o gfa:syng:blunt,k=63,s=8,seed=7 \
  --sequence-files panel.agc --force-large-region -O c4.syng-blunt

# Native syng overlap graph, for debugging
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o gfa --gfa-engine syng:raw \
  --sequence-files panel.agc --force-large-region -O c4.syng-raw
```

Available engines:

```text
pggb          default; sweepga + seqwish + smoothing + gfaffix
seqwish       unsmoothed graph
poa           small-region MSA graph
pggb:10000    partitioned mode with 10 kb windows
syng          syng syncmer graph, defaults to syng:blunt
syng:blunt    syng graph processed through pangenome/bluntg; links are 0M
syng:raw      native syng overlap graph
```

## Raw Syng Graph Dump

`impg syng2gfa` is different from `impg query -o gfa`.

```bash
impg syng2gfa -a panel.syng --sequence-files panel.agc -o panel.syncmers.gfa
```

This dumps the whole syng syncmer graph: one GFA segment per syncmer plus gap
segments between syncmers. By default this writes blunt graph output via
pangenome/bluntg. Use `--gfa-mode raw` when you want to inspect native syng
overlaps. Use `impg query -o gfa` when you want a local graph for a genomic
region.
