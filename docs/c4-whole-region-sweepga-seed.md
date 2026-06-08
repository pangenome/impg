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

After `one_many_scaffold0.initial` was selected as the best initial graph, I
ran the missing SPOA-only small-bubble pass. This deliberately used
`--method poa`, not POASTA:

```bash
impg crush \
  --gfa /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_scaffold0.initial.gfa \
  --output /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_scaffold0.spoa5.gfa \
  --method poa \
  --max-iterations 5 \
  --max-traversal-len 10k \
  --max-median-traversal-len 1k \
  --max-total-sequence 1m \
  --max-traversals 10k \
  --threads 32 \
  -v 1
```

The SPOA pass reached the fifth configured round and then stopped changing:
rounds 1-4 resolved 759, 618, 157, and 1 sites, and round 5 found no eligible
candidates. The full log is
`logs/one_many_scaffold0.spoa5.crush.log`.

The follow-up low-min-match seed graphs used the same collected FASTA and the
same FastGA/seqwish path, with only `--min-match-len 1` added:

```bash
impg graph \
  --sequence-files c4_whole_region.fa \
  --gfa-engine seqwish \
  --fastga \
  --num-mappings <MODE> \
  --scaffold-filter <MODE> \
  --scaffold-jump 0 \
  --min-match-len 1 \
  -g graphs/<label>.initial.gfa \
  -t 32 -v 1
```

Follow-up modes run:

- `one_many_minmatch1_scaffold0`: `<MODE>=1:many`
- `many_many_minmatch1_scaffold0`: `<MODE>=many:many`

The strict 1:1 low-min-match follow-up was run later in a separate timestamped
artifact directory:

`/home/erikg/impg/data/c4_whole_region_sweepga_seed_one_one_minmatch1_20260603T155138Z`

It used the same collected FASTA and the same FastGA/seqwish path, with strict
filters on both mapping axes:

```bash
/home/erikg/.cargo/bin/impg graph \
  --sequence-files /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa \
  --gfa-engine seqwish \
  --fastga \
  --num-mappings 1:1 \
  --scaffold-filter 1:1 \
  --scaffold-jump 0 \
  --min-match-len 1 \
  --temp-dir /tmp/c4_451_20260603T155138Z \
  --threads 32 \
  -v 1 \
  --output /home/erikg/impg/data/c4_whole_region_sweepga_seed_one_one_minmatch1_20260603T155138Z/graphs/one_one_minmatch1_scaffold0.initial.gfa
```

The short `--temp-dir` avoided the long temporary FastGA path form that had
triggered `FAtoGDB`'s buffer overflow in the first low-min-match attempt.

The first low-min-match attempt was launched from the output directory and
failed before alignment: `FAtoGDB` aborted with `*** buffer overflow detected
***` while preparing the temporary GDB from
`/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/.tmp6YkA6c.fa`.
The captured stderr is
`logs/one_many_minmatch1_scaffold0.graph.attempt1.stderr.log`; the timed exit
status was 1. This was not a parser rejection of `--min-match-len 1`. Rerunning
from the WG worktree, so the temporary FastGA paths lived under the shorter
`/home/erikg/impg/.wg-worktrees/agent-447`, succeeded for both modes.

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
| SPOA `one_many_scaffold0.spoa5` | 0:36.45 | 1.54 | 0 |
| graph `one_many_minmatch1_scaffold0` attempt 1 | 0:01.40 | 0.01 | 1 |
| graph `one_many_minmatch1_scaffold0` rerun | 5:35.78 | 10.58 | 0 |
| graph `many_many_minmatch1_scaffold0` | 6:17.99 | 11.25 | 0 |
| graph `one_one_minmatch1_scaffold0` | 5:16.12 | 10.46 | 0 |

All graph-report, sort, and render commands exited 0. Their timings are in
`runtime_summary.tsv`; the SPOA graph-report took 7.51 seconds, the Ygs sort
took 11.77 seconds, and the gfalook render took 1.84 seconds. The two
min-match-1 graph-report commands took 9.99 and 10.95 seconds, their Ygs sorts
took 15.67 and 20.86 seconds, and their gfalook renders took 1.97 and 2.01
seconds.

For `one_one_minmatch1_scaffold0.initial`, graph-report took 9.94 seconds,
`gfasort -p Ygs` took 16.14 seconds, `gfalook -m` took 2.00 seconds, and the
PNG upload took 1.29 seconds. All exited 0.

## Metrics

Key `impg graph-report --format tsv` metrics:

| graph | segments | total segment bp | bp-weighted coverage | singleton bp | path jump p99 | white-space p99 bp | white-space bridges >=1 kb |
|---|---:|---:|---:|---:|---:|---:|---:|
| `many_many_scaffold0.initial` | 7,431 | 234,157 | 455.63 | 2,921 | 28 | 168 | 28,025 |
| `many_many_scaffold0.poasta2` | 7,886 | 248,237 | 429.79 | 4,038 | 3,427 | 110,838 | 223,826 |
| `one_many_scaffold0.initial` | 7,411 | 234,828 | 454.33 | 2,922 | 5 | 66 | 8,790 |
| `one_many_scaffold0.poasta2` | 7,860 | 248,941 | 428.57 | 4,152 | 3,564 | 113,389 | 214,252 |
| `one_many_scaffold0.spoa5` | 8,387 | 250,010 | 426.74 | 4,609 | 3,445 | 112,207 | 214,719 |
| `ninety_ninety_scaffold0.initial` | 7,431 | 234,157 | 455.63 | 2,921 | 23 | 144 | 26,157 |
| `ninety_ninety_scaffold0.poasta2` | 7,886 | 248,237 | 429.79 | 4,038 | 3,427 | 110,838 | 223,826 |
| `one_many_minmatch1_scaffold0.initial` | 8,186 | 232,093 | 459.68 | 2,808 | 6,687 | 189,754 | 256,597 |
| `many_many_minmatch1_scaffold0.initial` | 8,154 | 230,664 | 462.53 | 2,795 | 19 | 24 | 34,018 |
| `one_one_minmatch1_scaffold0.initial` | 8,186 | 232,100 | 459.67 | 2,810 | 16 | 111 | 32,989 |

The best run by the graph-report stress metrics is
`one_many_scaffold0.initial`: it has the lowest path jump p99, lowest
white-space p99, and lowest count of >=1 kb white-space bridges. The two
unbounded modes are nearly identical; `90:90` was effectively unbounded here.

The POASTA polish preserved path spelling but made this graph-report objective
worse: it increased segments, singleton bp, path jump p99, white-space p99, and
long white-space bridge counts in all three modes. The configured small-bubble
pass should therefore not be accepted as an improvement for this run.

The SPOA-only pass also preserved path spelling, and it did exhaust the small
local candidate queue by round 5, but it did not fix the local-loop concern
without reproducing the long-jump regression. Relative to
`one_many_scaffold0.initial`, SPOA increased path jump p99 from 5 to 3,445,
white-space p99 from 66 bp to 112,207 bp, and >=1 kb white-space bridges from
8,790 to 214,719. Relative to `one_many_scaffold0.poasta2`, SPOA is only
slightly lower on path jump p99 and white-space p99, while singleton bp and
long white-space bridge count are higher. The result is therefore
path-preserving but not quality-improving.

The low-min-match follow-up did not improve the seed choice. Lowering
`--min-match-len` to 1 changed graph topology and increased segment count in
both modes. `one_many_minmatch1_scaffold0.initial` has slightly lower stored
segment bp than `one_many_scaffold0.initial`, but it regresses the topology
metrics badly: path jump p99 rises from 5 to 6,687, white-space p99 rises from
66 bp to 189,754 bp, and >=1 kb white-space bridges rise from 8,790 to 256,597.
`many_many_minmatch1_scaffold0.initial` has the best white-space p99 in the
table at 24 bp and lower singleton bp than the current best, but it still has
more segments, higher path jump p99, and 34,018 long white-space bridges. Its
render also shows more large bottom jump structure than
`one_many_scaffold0.initial`. Since neither min-match-1 graph is visibly and
topologically better than the current best seed, no SPOA polish was run for
either follow-up graph.

The strict 1:1 min-match-1 follow-up is a large improvement over
`one_many_minmatch1_scaffold0.initial`, but not a decisive improvement over
`many_many_minmatch1_scaffold0.initial`. Relative to one-many min-match-1,
strict 1:1 keeps the same segment count while reducing path jump p99 from
6,687 to 16, white-space p99 from 189,754 bp to 111 bp, and >=1 kb
white-space bridges from 256,597 to 32,989. Relative to many-many min-match-1,
strict 1:1 improves path jump p99 from 19 to 16 and long bridge count from
34,018 to 32,989, but it worsens white-space p99 from 24 bp to 111 bp and has
more segments, total segment bp, and singleton bp. The uploaded mean-depth
render should therefore be read as a clear visual/metric improvement over
one-many min-match-1 and a mixed, not overall-better, result versus many-many
min-match-1.

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
  graphs/<label>.polished.gfa
```

All three initial-vs-POASTA comparisons passed, and the required
initial-vs-SPOA comparison also passed:

| mode | expected paths | observed paths | missing | extra | spelling mismatches |
|---|---:|---:|---:|---:|---:|
| `many_many_scaffold0` | 465 | 465 | 0 | 0 | 0 |
| `one_many_scaffold0` | 465 | 465 | 0 | 0 | 0 |
| `ninety_ninety_scaffold0` | 465 | 465 | 0 | 0 | 0 |
| `one_many_scaffold0.spoa5` | 465 | 465 | 0 | 0 | 0 |

I also compared exact GFA P-line path-name sets against the collected FASTA
headers. Every initial, POASTA-polished, SPOA-polished, and min-match-1 graph
had `missing=0` and `extra=0`; the two min-match-1 graphs each had exactly 465
P-line path names.

| graph | P-line path names | missing FASTA names | extra GFA names |
|---|---:|---:|---:|
| `one_many_minmatch1_scaffold0.initial` | 465 | 0 | 0 |
| `many_many_minmatch1_scaffold0.initial` | 465 | 0 | 0 |
| `one_one_minmatch1_scaffold0.initial` | 465 | 0 | 0 |

Validation limitation: exact path names are validated against the collected
FASTA headers, and exact path spellings are validated from each initial GFA to
its polished GFA. The current source validator did not validate graph spellings
directly against the collected FASTA because it does not parse the preserved
`(+)/(-)` suffix form.

The strict 1:1 min-match-1 run used the same direct FASTA-header-to-GFA-path
set comparison. It observed 465 FASTA headers, 465 GFA paths, zero missing
names, zero extra names, and status `PASS`; see
`/home/erikg/impg/data/c4_whole_region_sweepga_seed_one_one_minmatch1_20260603T155138Z/validation/path_name_set_summary.tsv`.

## Renders

All original seven graphs and the three min-match-1 follow-up graphs were sorted
with `gfasort -p Ygs -t 32` and rendered with `gfalook -m`.

Uploaded stable HyperVolume render URLs:

- https://hypervolu.me/~erik/impg/c4-whole-region-sweepga-seed-one-many-initial.png
- https://hypervolu.me/~erik/impg/c4-whole-region-sweepga-seed-one-many-poasta2.png
- https://hypervolu.me/~erik/impg/c4-whole-region-sweepga-seed-one-many-spoa5.png
- https://hypervolu.me/~erik/impg/c4-whole-region-sweepga-seed-many-many-initial.png
- https://hypervolu.me/~erik/impg/c4-whole-region-sweepga-seed-one-many-minmatch1-scaffold0-initial.png
- https://hypervolu.me/~erik/impg/c4-whole-region-sweepga-seed-many-many-minmatch1-scaffold0-initial.png
- https://hypervolu.me/impg/c4-whole-region-sweepga-seed-one-one-minmatch1-scaffold0-initial.png

Visual judgment from the mean-depth renders and the report metrics is that the
unpolished `one_many_scaffold0.initial` graph is the cleanest product from this
run. The existing `many_many_scaffold0.initial` render is close, but its path
jump p99 and long white-space bridge count are worse than one-many. The
min-match-1 renders do not overturn that result: one-many min-match-1 has a
large long-jump regression, and many-many min-match-1 adds enough bottom jump
structure and long bridge count that it is not a safer seed for small-loop
polish. Strict 1:1 min-match-1 removes the severe one-many min-match-1
long-jump structure, but its mixed metrics versus many-many min-match-1 do not
make it a better seed than the original one-many scaffold0 graph.

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
| this run: `one_many_scaffold0.spoa5` | 8,387 | 250,010 | 426.74 | 4,609 | 214,719 |
| this run: `one_many_minmatch1_scaffold0.initial` | 8,186 | 232,093 | 459.68 | 2,808 | 256,597 |
| this run: `many_many_minmatch1_scaffold0.initial` | 8,154 | 230,664 | 462.53 | 2,795 | 34,018 |
| this run: `one_one_minmatch1_scaffold0.initial` | 8,186 | 232,100 | 459.67 | 2,810 | 32,989 |

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
syng-topology seed, not a many-many seed, not a lower-min-match seed including
the strict 1:1 min-match-1 follow-up, and not a small-bubble-polished
derivative.

The configured two-iteration POASTA polish and the five-iteration SPOA polish
are path-preserving but not quality-improving for this graph. SPOA exhausted
eligible small-bubble candidates by round 5, but it did not repair local loops
without introducing the POASTA-style long-jump and white-space regression. Both
polish choices should be treated as failed polish choices for this specific C4
whole-region seed rather than as evidence against the whole-region
FastGA/seqwish seed itself.

No SweepGA upstream filtering behavior was changed for this task, and no code
changes were required.
