# C4 Small-Bubble Threshold Sweep

Task: `c4-threshold-sweep`

Run root:

```text
/home/erikg/impg/.wg-worktrees/agent-450/data/c4_threshold_sweep_20260603T154441Z
```

Input seed:

```text
/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa
```

Binary:

```text
/home/erikg/.cargo/bin/impg
```

Local release rebuild was attempted after initializing submodules, but the task
message at 2026-06-03 15:44 UTC explicitly directed this run to bypass local
build repair and use `/home/erikg/.cargo/bin/impg`. The failed rebuild path was
blocked by `wfmash-rs` CMake selecting `/usr/bin/cc` under the Guix glibc
loader, producing a glibc mismatch. No more build time was spent after that
instruction.

## Validation

All seven requested matrix entries completed with exit status 0 for `crush`,
`graph-report`, `gfasort -p Ygs`, `gfalook -m`, and `scp`.

Path validation passed for every output:

- 465 paths present.
- No missing or extra path names relative to the seed.
- 0 path spelling mismatches relative to the seed GFA path sequences.
- 0 spelling mismatches relative to `c4_whole_region.fa`.

The installed `gfasort` exposes the requested pipeline flag as `-p/--pipeline`,
so the runs used `-p Ygs`.

## Outputs

Machine-readable summaries:

```text
metrics.tsv
summary.tsv
path-validation.tsv
uploads.tsv
manifest.json
seed.graph-report.tsv
```

Contact sheet:

```text
c4-smallbubble-contact-sheet.png
```

Uploaded render URLs:

```text
https://hypervolu.me/~erik/impg/c4-smallbubble-poa-median100.png
https://hypervolu.me/~erik/impg/c4-smallbubble-poa-median500.png
https://hypervolu.me/~erik/impg/c4-smallbubble-poa-median1000.png
https://hypervolu.me/~erik/impg/c4-smallbubble-poa-median2000.png
https://hypervolu.me/~erik/impg/c4-smallbubble-poasta-median1000.png
https://hypervolu.me/~erik/impg/c4-smallbubble-poasta-median2500.png
https://hypervolu.me/~erik/impg/c4-smallbubble-poasta-median5000.png
```

Remote upload confirmation:

```text
-rw-r--r-- 1 erik erik 1.1M Jun  3 15:48 www/impg/c4-smallbubble-poa-median100.png
-rw-r--r-- 1 erik erik 1.1M Jun  3 15:49 www/impg/c4-smallbubble-poa-median500.png
-rw-r--r-- 1 erik erik 1.1M Jun  3 15:50 www/impg/c4-smallbubble-poa-median1000.png
-rw-r--r-- 1 erik erik 1.1M Jun  3 15:51 www/impg/c4-smallbubble-poa-median2000.png
-rw-r--r-- 1 erik erik 1.2M Jun  3 15:52 www/impg/c4-smallbubble-poasta-median1000.png
-rw-r--r-- 1 erik erik 1.2M Jun  3 15:53 www/impg/c4-smallbubble-poasta-median2500.png
-rw-r--r-- 1 erik erik 1.1M Jun  3 15:54 www/impg/c4-smallbubble-poasta-median5000.png
```

## Baseline

Seed graph-report baseline:

| Graph | S | L | bp | bp-wtd depth | singleton bp | seg white-space bp | path ws p99/max | path ws bridges |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| seed | 8,186 | 11,303 | 232,093 | 459.684 | 2,808 | 16,967,514 | 189,754 / 217,619 | 256,597 |

## Summary Table

`summary.tsv` contains the full narrow table. The compact view below keeps the
same fields requested by the task: runtime, max RSS, resolved sites per round,
graph size, node/path-depth metrics, path whitespace and jump metrics, and a
visual/metric recommendation.

| Run | Wall | RSS GiB | Resolved per round | Graph S/L/bp | Node depth mean/bp-wtd/p50 | Path depth p50/p95/max | Whitespace and jumps | Recommendation |
| --- | ---: | ---: | --- | --- | --- | --- | --- | --- |
| poa median 100 | 0:33.76 | 1.81 | r1=748/748, r2=597/597, r3=151/151 | 9,079 / 12,339 / 240,589 | 533.45 / 443.45 / 464 | 464 / 931 / 161,984 | seg-ws=20,865,810; path-ws p99/max/bridges=100,506/240,054/431,597; link jump p99/max=8,147/9,016; path jump p99/max=3,682/9,016 | Best representative. All POA thresholds were identical. Safer than POASTA by bp-depth and whitespace, but still worse than seed on size/depth/singletons. |
| poa median 500 | 0:34.90 | 1.82 | r1=748/748, r2=597/597, r3=151/151 | 9,079 / 12,339 / 240,589 | 533.45 / 443.45 / 464 | 464 / 931 / 161,984 | seg-ws=20,865,810; path-ws p99/max/bridges=100,506/240,054/431,597; link jump p99/max=8,147/9,016; path jump p99/max=3,682/9,016 | Same output as other POA thresholds. |
| poa median 1000 | 0:35.41 | 1.83 | r1=748/748, r2=597/597, r3=151/151 | 9,079 / 12,339 / 240,589 | 533.45 / 443.45 / 464 | 464 / 931 / 161,984 | seg-ws=20,865,810; path-ws p99/max/bridges=100,506/240,054/431,597; link jump p99/max=8,147/9,016; path jump p99/max=3,682/9,016 | Same output as other POA thresholds. |
| poa median 2000 | 0:33.63 | 1.81 | r1=748/748, r2=597/597, r3=151/151 | 9,079 / 12,339 / 240,589 | 533.45 / 443.45 / 464 | 464 / 931 / 161,984 | seg-ws=20,865,810; path-ws p99/max/bridges=100,506/240,054/431,597; link jump p99/max=8,147/9,016; path jump p99/max=3,682/9,016 | Same output as other POA thresholds. |
| poasta median 1000 | 0:35.67 | 1.77 | r1=748/748, r2=629/629, r3=161/161, r4=3/3 | 8,573 / 11,833 / 241,841 | 539.46 / 441.16 / 464 | 464 / 931 / 161,984 | seg-ws=21,447,990; path-ws p99/max/bridges=103,177/241,388/432,042; link jump p99/max=7,695/8,510; path jump p99/max=3,810/8,510 | Not preferred. Fewer S/L than POA, but worse total bp, bp-weighted depth, singleton bp, and whitespace. |
| poasta median 2500 | 0:35.18 | 1.76 | r1=748/748, r2=629/629, r3=161/161, r4=3/3 | 8,573 / 11,833 / 241,841 | 539.46 / 441.16 / 464 | 464 / 931 / 161,984 | seg-ws=21,447,990; path-ws p99/max/bridges=103,177/241,388/432,042; link jump p99/max=7,695/8,510; path jump p99/max=3,810/8,510 | Same output as other POASTA thresholds. |
| poasta median 5000 | 0:34.58 | 1.78 | r1=748/748, r2=629/629, r3=161/161, r4=3/3 | 8,573 / 11,833 / 241,841 | 539.46 / 441.16 / 464 | 464 / 931 / 161,984 | seg-ws=21,447,990; path-ws p99/max/bridges=103,177/241,388/432,042; link jump p99/max=7,695/8,510; path jump p99/max=3,810/8,510 | Same output as other POASTA thresholds. |

## Recommendation

Do not promote these small-bubble crush outputs as a clear improvement over the
visually strong seed. The no-admission-gate sweep does preserve all paths and
reduces path white-space p99, but it increases graph size and singleton bp,
lowers bp-weighted node depth, and increases total segment white-space and long
white-space bridge counts.

If a representative output is needed, use `poa_median100` or any POA row, since
all POA thresholds are byte/metric equivalent here and POA has better
bp-weighted depth and whitespace than POASTA. POASTA visually looks slightly
more compact in the central rendered block and has fewer segments/links than
POA, but the metric regression is larger, so I do not recommend POASTA from this
sweep.

The threshold values did not change outputs within either method in this matrix.
