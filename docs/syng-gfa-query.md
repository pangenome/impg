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

Add `--render-graph` to also write a 1D `gfalook` rendering of the final graph:

```bash
impg query \
  -a panel.syng \
  -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k \
  -x \
  -o gfa:syng:crush \
  --sequence-files panel.agc \
  --force-large-region \
  --render-graph \
  --render-graph-depth \
  -O c4.syng-crush
```

This writes `c4.syng-crush.gfa` and `c4.syng-crush.png`. PNG is the default;
use `--render-graph-format svg` or `--render-graph-output c4.svg` for SVG.
`--render-graph-depth` passes `-m` to `gfalook`, giving the mean-depth 1D
rendering that makes path-depth architecture easier to see.
The renderer consumes the final graph after the configured graph transforms,
so syng output is rendered after the default `Ygs` sort unless `:nosort` is
explicitly requested.

Important flags:

- `-d 50k` is required. It merges query-gathered ranges separated by at most
  this distance and sets the largest one-hop gap/SV absorbed into one interval.
- `-x` enables transitive query. For duplicated loci, this is usually what you
  want.
- `--sequence-files` can be FASTA or AGC.
- `--force-large-region` is needed for large GFA/MAF requests.
- `-O` is an output prefix; the `.gfa` suffix is added by `impg`.
- Add `:cut-ns` to the syng engine, for example `-o gfa:syng:cut-ns:crush`,
  when fetched gap DNA contains assembly N-runs that should break paths rather
  than become path-private N segments. `cut-n-min-run=10` raises the minimum
  N-run length.

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

With `--render-graph`, images are written per BED row as well. By default they
go beside the graph outputs under `regions`; `--render-graph-output renders`
uses a separate render directory. Filenames are sanitized from BED column 4.

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

# Syng-native crush graph, sorted with the default gfasort Ygs pipeline
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o gfa:syng:crush \
  --sequence-files panel.agc --force-large-region -O c4.syng-crush

# VCF calls from the same local syng graph, using POVU
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o vcf:syng \
  --sequence-files panel.agc --force-large-region -O c4.syng-blunt

# Same thing with explicit mode and syncmer-parameter assertion
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o gfa:syng:blunt,k=63,s=8,seed=7 \
  --sequence-files panel.agc --force-large-region -O c4.syng-blunt

# Native syng overlap graph, for debugging
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o gfa --gfa-engine syng:raw \
  --sequence-files panel.agc --force-large-region -O c4.syng-raw

# Faster visualization-oriented order; use :nosort to preserve construction order
impg query -a panel.syng -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 50k -x -o gfa:syng:crush:sort,pipeline=Yg \
  --sequence-files panel.agc --force-large-region -O c4.syng-crush.yg

# Syng-native blunt graph for every discovered partition
impg partition -a panel.syng -w 1m -d 50k -o gfa --gfa-engine syng \
  --sequence-files panel.agc --separate-files --output-folder parts.syng-gfa
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
syng:crush    syng:blunt plus exact path-preserving bubble resolution
syng:...:sort,pipeline=Ygs
              final gfasort ordering; Ygs is the default for syng GFA output
syng:...:nosort
              disable final gfasort ordering
syng:...:nomask
              disable raw-layer high-frequency syncmer filtering
syng:...:mask,top=0.001,max-occ=500
              filter high-frequency syncmers using a top-fraction and/or
              whole-index occurrence threshold; rare repeated-copy local
              contexts are split as per-occurrence clones
```

`query` and `partition` differ only in how they gather local intervals. Once
the intervals are selected, both commands call the same GFA engine dispatcher.
For `--gfa-engine syng`, the dispatcher extracts the selected source-forward
strings, builds a regional syng index with the loaded syncmer parameters, then
emits either `syng:blunt` or `syng:raw` GFA.

The `-o gfa:<spec>` shorthand is generic. `-o gfa:pggb`, `-o gfa:seqwish`,
`-o gfa:poa`, and `-o gfa:syng` are equivalent to `-o gfa --gfa-engine
<engine>`. For alignment-backed engines, the spec can also name the alignment
backend: `-o gfa:wfmash:seqwish`, `-o gfa:fastga:pggb`, or
`-o gfa:sweepga:seqwish`. Add `:crush` after the graph engine to induce the
GFA, run the shared exact path-preserving crush pass, and then emit the default
self-loop-normalized plus `Ygs`-sorted final graph; supported spellings include
`-o gfa:pggb:crush` and `-o gfa:seqwish:crush`.

The shorthand is parsed as a staged graph pipeline:

```text
gfa:<stage>[,<key=value>...][:<stage>[,<key=value>...]...]
```

For example, `gfa:syng,k=63,s=8,seed=7:blunt` and the legacy
`gfa:syng:blunt,k=63,s=8,seed=7` are both accepted. Graph transforms use the
same staged syntax.
`gfa:syng:crush,method=auto,max-median-traversal-len=1k` emits blunt syng GFA
and then runs exact path-preserving bubble resolution; see
`docs/graph-pipeline-dsl.md`.
Crush methods include `poa`, `biwfa`, `allwave`, and `sweepga`; for example,
`gfa:syng:crush,method=allwave,k-nearest=5,k-farthest=2` uses AllWave's
many-sequence pair selection and Mash orientation detection inside each
selected POVU bubble.
For syng GFA output, `impg` then runs an internal gfasort pass by default using
pipeline `Ygs` (path-guided ordering, grooming, then topological sort). The pipeline
can be changed with a sort stage, for example
`gfa:syng:crush:sort,pipeline=Yg`, or disabled with `gfa:syng:crush:nosort`.

Syng GFA extraction frequency-filters high-copy syncmer nodes by default before
raw syng GFA is handed to bluntg. Filtered syncmers are removed from the raw
topology and bridged by sequence, so Alu-like or other repeated syncmers do not
become shared graph glue. It also splits rare repeated-copy local syncmer
contexts, which catches single-syncmer loops where an otherwise single-copy
syncmer appears a second time in one path. Use `gfa:syng:nomask` or
`gfa:syng:nofilter` to disable this, or tune it with a stage such as
`gfa:syng:filter,top=0.001,max-occ=500:crush`.

The transform can also be run on an existing blunt GFA:

```bash
impg crush -g local.blunt.gfa -o local.crushed.gfa
```

VCF output uses the same engine dispatch and passes the local GFA to POVU:
`-o vcf:syng`, `-o vcf:pggb`, `-o vcf:seqwish`, and `-o vcf:poa` are
equivalent to `-o vcf --gfa-engine <engine>`.

If the graph already exists, call POVU directly through impg:

```bash
impg gfa2vcf -g c4.syng-crush.gfa -o c4.syng-crush.vcf -r ref_path_name
```

`-r/--reference-name` is a path/name hint for POVU and can be repeated. This
uses the same Rust POVU conversion as `query -o vcf`.

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
