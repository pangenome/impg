# Bug: `impg syng` segfaults (signal 11) on some inputs

Found 2026-09-22/23 on gemini while building a chr20 syncmer index for `impg infer`.

## Summary

- `impg syng` dies with signal 11, with no Rust panic or message, on a whole GRCh38 chr20, on a chr20 panel, and on single pieces of some haplotypes.
- Changing the syncmer settings makes some inputs pass and others fail, so the crash depends on which syncmers get picked.
- Two patterns stand out:
  - the one input tested that holds long N runs (whole GRCh38 chr20) crashes with every setting tried, and the same sequence cut at its N runs passes with two of them;
  - one centromeric satellite piece crashes, or runs for over 10 min, even with no N.
- The crash also depends on context. The same 648,981 bp piece crashes when indexed alone with `--syncmer-seed 13`, but passes inside a set of 85 pieces with the same setting.
- I have no backtrace: gemini has no gdb, and the release binary has no debug symbols. A build of HEAD 50a14f3 on gemini failed; `docs/bug-agc-appended-sample-range.md` has the other issue from the same session.

## Environment

- impg: `/tgen_labs/guarracino/condatools/impg/latest/bin/impg`, `impg 0.5.0`. A 5052c41 build (`/tgen_labs/guarracino/condatools/impg/5052c41/bin/impg`) also segfaults on the 35 Mb piece below.
- All test inputs are in `/tgen_labs/guarracino/genome-infer/syngtest/` on gemini:
  - `one.fa`: GRCh38 chr20, named `GRCh38#0#chr20:0-64444167`. 64,444,167 bp, of which 499,910 are N.
  - `nsplit.fa`: `one.fa` cut at every N run, keeping the 85 pieces of at least 1 kb. 63,944,064 bp, 0 N. Pieces are named `GRCh38#0#chr20:start-end`.
  - `last.fa`: `GRCh38#0#chr20:26646554-27295535`, a centromeric satellite piece. 648,981 bp, 0 N, identical to the same-named record in `nsplit.fa` (md5 checked).
  - `min.fa`: bases 100001-500000 of the panel piece `HG03579#2#CM089297.1:32240917-67462154`. 400,000 bp, 0 N.
- The full inputs, both from impg query on the fixed `human579+GRCh38.agc`:
  - `/scratch/aguarracino/genome-infer/chr20/panel.fa`: 24,079 sequences, 39,805,776,171 bp, cut at N runs;
  - the two arm panels in `/scratch/aguarracino/genome-infer/chr20_0_25500000/` and `chr20_30500000_64444167/`.

## Reproduce (on gemini)

```shell
I=/tgen_labs/guarracino/condatools/impg/latest/bin/impg
cd /tgen_labs/guarracino/genome-infer/syngtest
for f in one.fa nsplit.fa last.fa min.fa; do
  for a in "" "--syncmer-seed 13" "--smer-length 12" "--syncmer-length 127"; do
    timeout 600 $I syng $a -f $f -o x -t 1 > x.log 2>&1; echo "$f [$a] rc=$?"
  done
done
# rc 139 = segfault, rc 124 = killed by the 600 s timeout
```

## Observed (one thread, 600 s limit)

| input | default | `--syncmer-seed 13` | `--smer-length 12` | `--syncmer-length 127` |
| --- | --- | --- | --- | --- |
| `one.fa`, GRCh38 chr20 with 499,910 N | segfault | segfault | segfault | segfault |
| `nsplit.fa`, the same cut at N runs (85 pieces) | segfault | ok (5.3 s, 417,040 KB) | ok | segfault |
| `last.fa`, 649 kb satellite piece | segfault | segfault | timeout | timeout |
| `min.fa`, 400 kb of HG03579 | segfault | ok | ok | ok |

More observations:
- Default settings, one thread, on sub-ranges (1-based) of the same HG03579 piece:

| sub-range | result |
| --- | --- |
| 1-100,000 | ok |
| 1-200,000 | ok |
| 1-400,000 | ok |
| 1-500,000 | segfault |
| 1-600,000 | segfault |
| 100,001-500,000 (`min.fa`) | segfault |
| 100,001-600,000 | segfault |
| 200,001-600,000 | ok |
| 300,001-600,000 | ok |
| 400,001-600,000 | ok |
| 500,001-1,500,000 | ok |

- On `min.fa`, `--syncmer-length 31` and `--position-sample-rate 16` also segfault.
- On bases 1-1,000,000, `--parallel-dictionary` and `--position-sample-rate 64` both segfault.
- On a 35,221,237 bp N-free piece of the same contig, the segfault stays with `ulimit -s unlimited`, and also with the 5052c41 build.
- The whole chr20 panel, cut at N runs (24,079 sequences, `-t 16`), crashes with all four settings within 16 s.
- The two arm panels, which stay off the centromere (chr20:0-25,500,000 and chr20:30,500,000-64,444,167), both index fine with `--syncmer-seed 13`: 6 min 13 s and 10 min 27 s, 180 GB allowed.
- With `-t 16`, the last `Processing` log line before the crash is not a reliable pointer to the sequence that caused it.

## Possibly related

- Before the AGC fix, a panel whose GRCh38 record was a wrong 19,245,237 bp sequence (see `docs/bug-agc-appended-sample-range.md`) made syng run for over 1 hour on that record.
- Not yet checked whether `impg infer` works on the arm indexes. That runs next.

## Root cause and fix (found 2026-09-23)

Reproduced on a machine with gdb + valgrind using `last.fa` and `min.fa`. Backtrace of the
release-equivalent debug build:

```
addDirect           vendor/syng/rskip.c:1071   <- node[rs.dynamic->start].before, start = 0x5555
rsAddSyng           vendor/syng/rskip.c:1676
syngBWTadd          vendor/syng/syngbwt3.c:389   (negative-orientation branch, k = -324)
syngBWTpathAdd      vendor/syng/syngbwt3.c:540
add_sequence_internal  src/syng.rs:2646
```

Valgrind on the same input pins the origin:

```
Conditional jump or move depends on uninitialised value(s)
   at addDirect (rskip.c:1050)
Uninitialised value was created by a stack allocation at buildDynamic (rskip.c:608)
```

Root cause, in `buildDynamic` (`vendor/syng/rskip.c`):

- `symMax` starts at `-1` and is only set inside the run loop, so with `nRun == 0` it stays `-1`;
- `rs.dynamic->start = sTop[symMax]` then reads `sTop[-1]`, i.e. uninitialised stack memory
  (the crash showed `start = 21845 = 0x5555`, `max = 129`, `free = 128` - exactly what
  `buildDynamic` produces for `nRun = 0`);
- `addDirect` only treats `start == 0` as the empty-structure case, so the garbage `start`
  sends it to `node[21845]`, outside the array, and it segfaults.

The `nRun == 0` path is reachable: a `LINEAR_SYNG` rskip can hold many directory entries but no
runs (`convertToRskip` builds one with `nRun = 0`; directory entries accumulate via
`rsDirAddSyng` while the out side stays SIMPLE and short-circuits `rsAddSyng`). Once the
directory no longer fits in `MAX_LINEAR`, `rsDirAddSyng` converts to DYNAMIC with `nRun = 0`,
and the next `rsAddSyng` walks the garbage `start`. This is why the crash depends on syncmer
settings and on which (sync, offset) pairs a sequence produces.

Fix (merged in `pangenome/syng` main via PR #4, https://github.com/pangenome/syng/pull/4):

```c
rs.dynamic->start = nRun ? sTop[symMax] : 0 ; // nRun == 0 leaves no columns - start 0 marks empty
```

With the fix, `min.fa` and `last.fa` pass under all four settings from the matrix above
(default, `--syncmer-seed 13`, `--smer-length 12`, `--syncmer-length 127`), and the syng test
suite (`test_syng_integration`, `test_syng_startcount`) passes. `one.fa`/`nsplit.fa` were not
re-run locally (only `last.fa`/`min.fa` were copied off gemini); the N-run cases failed with the
same signature, so they are expected to be fixed too, but confirm on gemini.

Note: the earlier 600 s timeouts on `last.fa` with default settings disappeared as well - the
corrupted `start` was also sending `addDirect` into pathological walks before crashing.

## Original investigation notes

- A debug build with a backtrace on `last.fa` (smallest crash that fails with seed 13) and on `min.fa` (smallest default-setting crash).
- A guard that reports the sequence and position being indexed, instead of a bare segfault.
