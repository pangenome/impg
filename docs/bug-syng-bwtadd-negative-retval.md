# Bug: `impg syng` dies with `syngBWTadd -ve output retVal -3` on two HPRC r2 haplotypes

Found 2026-11-23 on gemini while building a syng index of 40 HPRC release 2 haplotypes for read partitioning.

## Summary

- `impg syng` stops with `FATAL ERROR: syngBWTadd -ve output retVal -3 callCount 433514526 k -32263223` and exit status 255.
- The crash is deterministic: every crashing input fails at the same node, `k = -32263223`, and inputs with the same sequences before the crash contig also print the same `callCount`.
- The smallest input that crashes is NA20346 hap1 (100 sequences, 2,949,496,588 bp) followed by one contig of hap2, `NA20346#2#JBIRDS010000003.1` (52,272,629 bp, 0 N): 101 sequences, 3,001,769,217 bp. The crash comes while indexing that contig.
- The crash depends on context:
  - `NA20346#2#JBIRDS010000003.1` alone passes;
  - the 17 hap2 sequences up to and including it, without hap1, pass;
  - NA20346 hap2 alone (108 sequences, the crash contig included) passes;
  - a 10-haplotype set of 30.6 Gb passes.
- It is not the total size. The crash comes after 5.1 Gb and 215 M syncmers, while 30.6 Gb passed with the same binary.
- `k` is negative, so the failing call is in the negative-orientation branch of `syngBWTadd`. `docs/bug-syng-segfault.md` also points at that branch (`k = -324`), in a segfault that also depended on context.

## Environment

- impg main a34a167 (2026-11-22, "Merge pull request #328 from pangenome/fix/drop-build-env"), with `vendor/syng` at d7e4e3c. Built on gemini with GCC 14 (`/tgen_labs/guarracino/condatools/_buildenv_gcc14`), `cargo build --release`. Binary at `/tgen_labs/guarracino/condatools/impg/a34a167/bin/impg`, run with `LD_LIBRARY_PATH=/tgen_labs/guarracino/condatools/_buildenv_gcc14/lib`.
- Syncmer parameters are the defaults, from the log: `paper k=63, paper s=8, internal k=8, w=55, seed=7`.
- Input sequences come from `agc getset` (agc 3.2.1) on `/tgen_labs/guarracino/alejandro/project-2-its/stage-3-hprc_scan/HPRC_r2_assemblies_0.6.1.agc`. HPRC release 2 has 465 haplotypes plus GRCh38, and no HG002.

## Inputs and results

All files are in `/scratch/aguarracino/tmp/its-reads/syng_bug/` on gemini, except the 40-haplotype input. Each run is `impg syng -f X.fa -o X -t 1 -v 1`.

| input | sequences | bp | result | wall | max RSS |
|---|---:|---:|---|---:|---:|
| `pan40.fa` (40 haplotypes, `-t 32`) | 3,401 | 121,367,231,852 bytes of FASTA | crash after 115 sequences, 5,102,438,437 bp, 214,995,528 syncmers | 3:06 | 11.1 GB |
| `both.fa` (NA20346 hap1 + hap2) | 208 | 5,998,862,381 | crash, same message | 214.6 s | 11.1 GB |
| `upto.fa` (hap1 + hap2 up to `JBIRDS010000003.1`) | 117 | 5,154,711,066 | crash, same message | 232.8 s | 11.1 GB |
| `h1c.fa` (hap1 + `JBIRDS010000003.1` only) | 101 | 3,001,769,217 | crash, `callCount 251800108 k -32263223` | 207.8 s | 10.2 GB |
| `h2c.fa` (the 17 hap2 sequences of `upto.fa`, no hap1) | 17 | 2,205,214,478 | passes, exit 0 | 152.5 s | 9.4 GB |
| `hap2.fa` (NA20346 hap2 alone) | 108 | | passes, exit 0 | 205.8 s | 10.6 GB |
| `one.fa` (`NA20346#2#JBIRDS010000003.1` alone) | 1 | 52,272,629 (0 N) | passes, exit 0 | 3.3 s | 0.35 GB |
| `hprc10.fa` (HG00097, HG00099, HG00126, HG00128, HG00133, both haplotypes; `-t 32`) | | 30.6 Gb | passes | 15:52 | 13.7 GB |

- `pan40.fa` is at `/scratch/aguarracino/tmp/its-reads/hprc10x/pan40/pan40.fa`, with its haplotype list in `sets.txt`. NA20346 hap1 and hap2 are its first two sets. Its log is `/scratch/aguarracino/tmp/its-reads/hprc10x/hx_pan40.41310649.log`.
- `hprc10.fa` is at `/scratch/aguarracino/tmp/its-reads/syng_pilot/hprc10.fa`.
- The last log lines before the crash in `pan40` are:

```
Indexed NA20346#2#JBIRDS010000002.1: 140364940 bp, 5791077 syncmers, indexed=true, elapsed 3.989s, total 116 sequences / 5102438437 bp / 214995528 syncmers at +3m3.
Processing NA20346#2#JBIRDS010000003.1 (52272629 bp)
FATAL ERROR: syngBWTadd -ve output retVal -3 callCount 433514526 k -32263223
```

## Reproduce (on gemini)

```shell
T=/tgen_labs/guarracino/condatools; export LD_LIBRARY_PATH=$T/_buildenv_gcc14/lib
cd /scratch/aguarracino/tmp/its-reads/syng_bug
$T/impg/a34a167/bin/impg syng -f h1c.fa -o h1c -t 1 -v 1       # 3.5 min, 10 GB; exit 255, callCount 251800108 k -32263223
# to rebuild upto.fa from the archive:
AGC=/tgen_labs/guarracino/alejandro/project-2-its/stage-3-hprc_scan/HPRC_r2_assemblies_0.6.1.agc
$T/agc/3.2.1/bin/agc getset $AGC NA20346_hap1_hprc_r2_v1.0.1 NA20346_hap2_hprc_r2_v1.0.1 > both.fa
$T/samtools/1.24/bin/samtools faidx both.fa
grep '#1#' both.fa.fai | cut -f1 > h1c.list; echo 'NA20346#2#JBIRDS010000003.1' >> h1c.list
$T/samtools/1.24/bin/samtools faidx both.fa -r h1c.list > h1c.fa
```

The scripts are `syng_bug.sh` and `syng_bug2.sh` in the same folder.

## Where it dies

`vendor/syng/syngbwt3.c`, `syngBWTadd`, the check at the end:

```c
if (retVal < 0) die ("syngBWTadd -ve output retVal %d callCount %llu k %d", retVal, callCount, k) ;
```

- `k < 0`, so the call takes the `else` branch ("have to come from .out == -in and leave on .in == -out").
- In that branch `retVal` is either `j` or the value of `rsAddSyng (&n->in.rs, j, -out, outOff)`.
- −3 is an error code, not a rank: `RS_ERROR_INDEX_OVERFLOW` in `vendor/syng/rskip.h`. `addDirect` in `rskip.c` returns it when it is asked to insert at a position past the end of the list (`index %u > sum %d`). So `j` was larger than the number of paths that leave node 32263223 on its `in` side.

## Cause

- A simple node side stores its offset in a `U16`. When the edge that creates a node is longer than 65,535 bp (a stretch with no syncmer), `nodeCreateOnePath` (added by syng `dc18190`) builds that side as an Rskip instead.
- It gave that Rskip one run and a directory count of 1, and set the count of the other side to 1. A run counts a path that leaves through the side, a directory count a path that enters through it. So the node records a second path, in the opposite direction, that never existed.
- The neighbour nodes do not know that path. Later paths that use the same edges get a rank one too high. Simple nodes pass a rank on unchanged, so the error travels along the path to the next node where paths branch. There the insert position is past the end and the build dies. If the shifted rank is still in range, the build finishes and the GBWT walks the wrong nodes, with no error.
- In NA20346 hap1 alone (old binary, `fix/h1_old.1gbwt`) there are 3 edges longer than 65,535 bp: 82,217, 597,338 and 207,138 bp. At both ends of an edge the counts must be equal. On the 597,338 bp edge they are 2 at one end and 1 at the other (node 93558372 `e 93558373 597338 2`, node 93558373 `E 93558372 597338 1`), and the same on the 207,138 bp edge.
- Crash node 32263223 has 2 edges in and 2 edges out. One branch is an unbranched chain of 85 nodes, 93558288 to 93558372, that ends at the 597,338 bp edge (`fix/h1_old.nbhd.txt`). Hap1 alone builds and walks correctly. The hap2 contig `JBIRDS010000003.1` uses that chain, and its rank overflows at 32263223.
- The crash contig alone has no edge longer than 65,535 bp (`one.1gbwt`), which is why it passes on its own.

## Fix

`vendor/syng/syngbwt3.c`, `nodeCreateOnePath` and `rskipBuildSingleDynamicSyng`. The new node gets only the path that created it.

| node side | path enters by this side | path leaves by this side |
|---|---|---|
| simple | count 1 | count 0 |
| Rskip | directory count 1, no run | run of 1, directory count 0 |

Test `test_long_edge_paths_walk_back_exactly` in `tests/test_syng_startcount.rs`:
- It adds sequences forward and reverse, as `SyngIndex` does, and walks every path back, both in memory and after a write and read of the `.1gbwt`.
- Its two cases are the smallest ones found by a random search. The first is one sequence with syncmers 2, 1, 2, 1 at positions 0, 20, 70,020 and 70,040. The second is two sequences: (1,0) (1,20) (−2,70,020), and (1,0) (−2,70,000) (1,70,020).
- On gemini (job 41325710), with the fix all 9 tests of the file pass. With the fix reverted, the test dies with `FATAL ERROR: syngBWTadd -ve output retVal -3 callCount 6 k -2` (exit 255).

Random test of the C code alone (2,000 graphs, 30 sequences of 2 to 9 syncmers, each added forward and reverse, 5 nodes):

| edge offsets | old: walks right | old: crash `-3` | old: walks wrong | fixed: walks right |
|---|---:|---:|---:|---:|
| 10, 20 bp | 2,000 | 0 | 0 | 2,000 |
| 10, 20 or 65,535 bp | 2,000 | 0 | 0 | 2,000 |
| 10, 20 or 65,536 bp | 997 | 59 | 944 | 2,000 |
| 10, 20 or 70,000 bp | 997 | 59 | 944 | 2,000 |

## Results with the fix (gemini)

- Fixed binary: `/scratch/aguarracino/tmp/its-reads/syng_bug/fix/impg.fix` (a34a167 plus the fix, md5 `264d50959892`). Outputs and logs: `/scratch/aguarracino/tmp/its-reads/syng_bug/fix/`.
- Check of an index: copy its 6 files and run `impg syng-repair -a COPY --force`, which walks every GBWT path again. If its `.syng.pstep`, `.syng.spos` and `.syng.names` are byte-identical to those of the build, every path walks back as it was added.

| input | old binary | fixed binary | re-walk of the index |
|---|---|---|---|
| `one.fa` | exit 0 | exit 0, all 6 files byte-identical to old | identical (old and fixed) |
| `h2c.fa` | exit 0 | exit 0, all 6 files byte-identical to old | identical (old and fixed) |
| `h1.fa` (hap1 alone) | exit 0, 266.8 s, 10.2 GB | not run | identical (old) |
| `h1c.fa` | crash after 207.8 s | exit 0, 182.0 s, 10.2 GB | identical (fixed) |
| `both.fa` | crash after 214.6 s | exit 0, 329.0 s, 12.0 GB | identical (fixed) |
| `hprc10.fa` (`-t 32`) | exit 0, 15:52 | exit 0, 1,349.9 s, 13.7 GB | old: `FATAL ERROR: shouldn't get here in rsFind` after 538 of 765 paths. Fixed: identical, 110.1 s |

- `hprc10`: the old and the fixed index differ only in `.1gbwt` (first difference at byte 2,534,605,571). `.1khash` and `.syng.meta` are identical. The fixed `.1gbwt` has 86 edge records with an offset above 65,535 bp.
- The old pilot index `/scratch/aguarracino/tmp/its-reads/syng_pilot/hprc10.*` built without error but cannot be walked. Anything made from it needs a rebuild.
