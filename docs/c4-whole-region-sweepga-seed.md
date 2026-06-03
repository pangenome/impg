# C4 whole-region SweepGA seed graph evaluation

Date: 2026-06-03

Output directory:

`/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z`

## Question

This run tests whether syng should be used only to collect C4 locus sequences,
with the initial local graph built from the whole collected region using
SweepGA/FastGA plus seqwish. No syng/flubble topology was used as the first
local graph.

## Inputs

- Syng prefix: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng`
- AGC: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc`
- Region: `GRCh38#0#chr6:31891045-32123783`
- Merge distance: `50k`
- Threads: `32`
- CLI used for the experiment: installed `impg 0.4.1`

FASTA collection produced `c4_whole_region.fa` with 465 records and
106,710,959 bytes. FASTA headers retained the PanSN/range/strand names from
query output, for example `CHM13#0#chr6:31744284-31976975(+)`.

## Commands

The one syng-derived step was sequence collection:

```bash
impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r 'GRCh38#0#chr6:31891045-32123783' \
  -d 50k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o fasta \
  --reverse-complement \
  -O /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region \
  -t 32 -v 1
```

The three initial graph modes were built from the collected FASTA with:

```bash
impg graph \
  --sequence-files c4_whole_region.fa \
  --gfa-engine seqwish \
  --fastga \
  --num-mappings <MODE> \
  --scaffold-filter <MODE> \
  --scaffold-jump 0 \
  -g graphs/<label>.initial.gfa \
  -t 32 -v 1
```

Modes run:

- `many_many_scaffold0`: `<MODE>=many:many`
- `one_many_scaffold0`: `<MODE>=1:many`
- `ninety_ninety_scaffold0`: `<MODE>=90:90`

The `90:90` syntax was accepted by the parser. On this input it kept the same
600,404 mappings as `many:many`.

Each initial graph was then polished with the small-bubble-only POASTA pass:

```bash
impg crush \
  -g graphs/<label>.initial.gfa \
  -o graphs/<label>.poasta2.gfa \
  --method poasta \
  --max-iterations 2 \
  --max-traversal-len 10k \
  --max-median-traversal-len 1k \
  --max-total-sequence 1m \
  --threads 32 -v 1
```

POASTA succeeded for all three modes; no POA fallback was needed. For every
initial and polished graph I also ran:

```bash
impg graph-report -g <graph.gfa> -o reports/<label>.graph-report.tsv --format tsv -t 32 -v 1
gfasort -i <graph.gfa> -o sorted/<label>.Ygs.gfa -p Ygs -t 32 -v 1
gfalook -i sorted/<label>.Ygs.gfa -o renders/<label>.Ygs.mean-depth.png -m -x 3200 -y 1800 -a 3 -t 32 -v 1
```

Complete command and timing logs are in `commands.md`, `runtime_summary.tsv`,
and `logs/` under the output directory.

## Runtime

Major commands, from `/usr/bin/time -v`:

| command | wall | max RSS GiB | status |
|---|---:|---:|---:|
| FASTA query | 5:07.80 | 51.57 | 0 |
| graph `many_many_scaffold0` | 5:21.33 | 9.49 | 0 |
| graph `one_many_scaffold0` | 5:16.64 | 9.11 | 0 |
| graph `ninety_ninety_scaffold0` | 5:20.51 | 9.50 | 0 |
| POASTA `many_many_scaffold0` | 0:15.59 | 1.32 | 0 |
| POASTA `one_many_scaffold0` | 0:15.24 | 1.31 | 0 |
| POASTA `ninety_ninety_scaffold0` | 0:15.70 | 1.32 | 0 |

All graph-report, sort, and render commands exited 0. Their timings are in
`runtime_summary.tsv`; each graph-report took about 7 seconds, each Ygs sort
about 10-11 seconds, and each gfalook render about 2 seconds.

## Metrics

Key `impg graph-report --format tsv` metrics:

| graph | segments | total segment bp | bp-weighted coverage | singleton bp | path jump p99 | white-space p99 bp | white-space bridges >=1 kb |
|---|---:|---:|---:|---:|---:|---:|---:|
| `many_many_scaffold0.initial` | 7,431 | 234,157 | 455.63 | 2,921 | 28 | 168 | 28,025 |
| `many_many_scaffold0.poasta2` | 7,886 | 248,237 | 429.79 | 4,038 | 3,427 | 110,838 | 223,826 |
| `one_many_scaffold0.initial` | 7,411 | 234,828 | 454.33 | 2,922 | 5 | 66 | 8,790 |
| `one_many_scaffold0.poasta2` | 7,860 | 248,941 | 428.57 | 4,152 | 3,564 | 113,389 | 214,252 |
| `ninety_ninety_scaffold0.initial` | 7,431 | 234,157 | 455.63 | 2,921 | 23 | 144 | 26,157 |
| `ninety_ninety_scaffold0.poasta2` | 7,886 | 248,237 | 429.79 | 4,038 | 3,427 | 110,838 | 223,826 |

The best run by the graph-report stress metrics is
`one_many_scaffold0.initial`: it has the lowest path jump p99, lowest
white-space p99, and lowest count of >=1 kb white-space bridges. The two
unbounded modes are nearly identical; `90:90` was effectively unbounded here.

The POASTA polish preserved path spelling but made this graph-report objective
worse: it increased segments, singleton bp, path jump p99, white-space p99, and
long white-space bridge counts in all three modes. The configured small-bubble
pass should therefore not be accepted as an improvement for this run.

## Validation

Direct source validation was attempted with:

```bash
/home/erikg/impg/target/release/examples/validate_gfa_path_sources \
  <graph.gfa> \
  c4_whole_region.fa
```

This validator could not parse the exact preserved FASTA/GFA path names because
the coordinate suffix includes strand, e.g. `CHM13#0#chr6:31744284-31976975(+)`.
For every completed graph it reported:

- `paths=465`
- `unparsable_paths=465`
- `spelling_mismatches=0`
- exit status `1`

Because the direct source validator could not perform source fetch validation
on these exact names, I used the documented fallback:

```bash
/home/erikg/impg/target/release/examples/compare_gfa_paths \
  graphs/<label>.initial.gfa \
  graphs/<label>.poasta2.gfa
```

All three initial-vs-polished comparisons passed:

| mode | expected paths | observed paths | missing | extra | spelling mismatches |
|---|---:|---:|---:|---:|---:|
| `many_many_scaffold0` | 465 | 465 | 0 | 0 | 0 |
| `one_many_scaffold0` | 465 | 465 | 0 | 0 | 0 |
| `ninety_ninety_scaffold0` | 465 | 465 | 0 | 0 | 0 |

I also compared exact GFA P-line path-name sets against the collected FASTA
headers. Every initial and polished graph had `missing=0` and `extra=0`.

Validation limitation: exact path names are validated against the collected
FASTA headers, and exact path spellings are validated from each initial GFA to
its polished GFA. The current source validator did not validate graph spellings
directly against the collected FASTA because it does not parse the preserved
`(+)/(-)` suffix form.

## Renders

All six completed graphs were sorted with `gfasort -p Ygs -t 32` and rendered
with `gfalook -m`.

Uploaded stable HyperVolume render URLs:

- https://hypervolu.me/~erik/impg/c4-whole-region-sweepga-seed-one-many-initial.png
- https://hypervolu.me/~erik/impg/c4-whole-region-sweepga-seed-one-many-poasta2.png

Visual judgment from the mean-depth renders and the report metrics is that the
unpolished `one_many_scaffold0.initial` graph is the cleanest product from this
run. The polished one-many render is path-preserving, but the long-jump and
white-space metrics show that the POASTA pass stretched/disconnected the
layout relative to the initial seed graph.

## Baseline comparison

Relevant C4 baselines from `docs/local-syng-parameter-sweep.md` and the latest
WG FastGA 8-round report:

| graph | segments | total segment bp | bp-weighted coverage | singleton bp | white-space bridges >=1 kb |
|---|---:|---:|---:|---:|---:|
| syng-local `k=63,s=8` C4 | 9,037 | 261,748 | 193.11 | 23,085 | 15,684 |
| global syng mask+crush C4 | 9,062 | 263,218 | 192.48 | 26,333 | 13,496 |
| PGGB control C4 | 7,170 | 89,342 | 595.08 | 523 | 9,607 |
| WG FastGA 8-round `/home/erikg/impg/data/c4_blocker_05b_fastga_wg_20260603T084706Z/graph-report.tsv` | 15,521 | 371,375 | 287.28 | 29,704 | 423,508 |
| this run: `one_many_scaffold0.initial` | 7,411 | 234,828 | 454.33 | 2,922 | 8,790 |
| this run: `one_many_scaffold0.poasta2` | 7,860 | 248,941 | 428.57 | 4,152 | 214,252 |

`one_many_scaffold0.initial` beats the current syng-derived C4 path on the
main stress indicators: fewer segments, much lower singleton bp, higher
bp-weighted coverage, and fewer long white-space bridges than both the local
syng `k=63,s=8` and global syng mask+crush baselines. It also beats the latest
WG FastGA 8-round result by a wide margin on segment count, segment bp,
singleton bp, and long white-space bridges.

It does not fully beat the PGGB control: PGGB still has much lower stored
segment bp and singleton bp. However, this whole-region FastGA/seqwish seed has
slightly fewer >=1 kb white-space bridges than the PGGB C4 control in the
local-syng sweep table.

## Conclusion

The experiment supports using syng only for C4 sequence collection and building
the first local graph from the whole collected region with FastGA/SweepGA plus
seqwish. The best graph here is `one_many_scaffold0.initial`, not any
syng-topology seed and not the small-bubble-polished derivative.

The configured two-iteration POASTA polish is path-preserving but not
quality-improving for this graph. It should be treated as a failed polish
choice for this specific C4 whole-region seed rather than as evidence against
the whole-region FastGA/seqwish seed itself.

No SweepGA upstream filtering behavior was changed for this task, and no code
changes were required.
