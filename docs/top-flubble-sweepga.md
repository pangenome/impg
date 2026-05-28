# Top-level flubble SweepGA crush

## Summary

`method=top-flubble-sweepga` (alias `chain-povu-sweepga`) resolves only POVU
level-0 flubbles as replacement boundaries. Descendant and neighboring bubble
content inside a top-level region is included in the same local SweepGA/seqwish
alignment problem, and no compression-quality guard is allowed to keep the
original region. The only acceptance checks are validity checks: replacement GFA
parses, path names are preserved, and path spellings are preserved.

The real C4 GRCh38 acceptance run completed successfully and preserved all
paths, but it does not fix the under-alignment/explosion problem. It proves
SweepGA does produce alignments for large top-level regions, but most small
regions still receive no backend alignments under the wfmash-safe floor, and the
final graph is much larger than all guarded baselines.

## Implementation

- Added `ResolutionMethod::TopFlubbleSweepga`.
- Accepted method names include `top-flubble-sweepga`,
  `top-level-flubble-sweepga`, and `chain-povu-sweepga`.
- The mode runs one POVU discovery pass, selects only `povu_level == 0`
  regions, materializes those non-overlapping regions, and applies all accepted
  replacements in a single frontier rewrite.
- Descendant POVU sites are counted for evidence but not selected as separate
  replacement jobs.
- Replacement building uses SweepGA/seqwish in local contained-region mode.
- For top-level mode, the code trusts SweepGA's own replacement-tier PAF filter
  and disables the seqwish-tail filter to avoid applying a second wrapper-side
  filtering decision after SweepGA has already made its documented decision.
- wfmash replacement jobs no longer withhold short non-empty records in impg;
  wfmash/MashMap own any minimum-length semantics and should fail explicitly if
  an input record is unsupported.
- Top-level replacement jobs are built sequentially from smaller to larger
  regions. This avoids concurrently running several large wfmash/seqwish jobs,
  which caused a 30 minute CPU-bound plateau with no completed regions.

## Final C4 command

Output directory:

```bash
/home/erikg/impg/data/c4_top_flubble_sweepga_wfmash_seq_20260527T172621Z
```

Command:

```bash
TMPDIR=/home/erikg/impg/data/c4_top_flubble_sweepga_wfmash_seq_20260527T172621Z/tmp \
RAYON_NUM_THREADS=8 \
LD_LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/home/erikg/micromamba/lib:${LD_LIBRARY_PATH:-} \
RUST_LOG=impg=info \
/usr/bin/time -v -o /home/erikg/impg/data/c4_top_flubble_sweepga_wfmash_seq_20260527T172621Z/time.txt \
  impg query \
    -t 32 \
    --temp-dir /home/erikg/impg/data/c4_top_flubble_sweepga_wfmash_seq_20260527T172621Z/tmp \
    -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
    -r 'GRCh38#0#chr6:31891045-32123783' \
    --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
    -d 50k \
    -o 'gfa:syng:mask,min-run=3:crush,method=top-flubble-sweepga,aligner=wfmash,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=0:nosort' \
    -O /home/erikg/impg/data/c4_top_flubble_sweepga_wfmash_seq_20260527T172621Z/run.nosort \
    -v 1
```

The command exited with status 0.

## C4 validation

Top-level decomposition and containment:

```text
POVU sites: 2441
Top-level flubbles: 1559
Descendant/non-root sites: 878
Sequence-variable top-level regions selected: 484
No cross-top-level merging: yes
Top-level size distribution: n=1559, min=30, p50=124, p90=171, max=42420, total=1213152
Selected traversal stats: max-len median/max=132/42362, median-len median/max=129/25155,
  p90-len median/max=129/42317, traversals max=465, total max=5132443
```

Replacement evidence:

```text
Applied replacements: 484/484
Round count: 1
Block count: 484
Regions with emitted alignments: 25/484
Total emitted alignment records: 2899229
Total aligned bp: 4202270125
Total PAF bytes after SweepGA filtering/dedup: 496727803
Build time: 295.21s
Rewrite plus validation time: 2.81s
No wfmash minimum-segment errors in final run
```

Path preservation:

```text
Reference no-crush C4 paths: 465
Output paths: 465
Missing path names: 0
Extra path names: 0
Path spelling mismatches versus no-crush reference: 0
Duplicate output path names: 0
```

Final graph metrics:

```text
Segments: 221342
Links: 397240
Paths: 465
Segment bp: 26037441
Trivial-stringy candidates: 398
Wall time: 9:56.16
Max RSS: 55212712 KB
```

The trivial-stringy count was computed with `/tmp/find_stringy_bubbles.py`.

PNG:

```text
Local: /home/erikg/impg/data/c4_top_flubble_sweepga_wfmash_seq_20260527T172621Z/c4-top-flubble-sweepga.png
Remote: www/impg/c4-top-flubble-sweepga.png
ssh ls: -rw-r--r-- 1 erik erik 865K May 27 17:59 www/impg/c4-top-flubble-sweepga.png
```

The graph was rendered directly from the final unsorted GFA. `gfasort -p Ygs`
was stopped after more than 15 minutes on the 221342 segment graph while still
in topo-sort; rendering the actual final GFA directly produced the required PNG.

## Comparison

| Run | Segments | Segment bp | Trivial-stringy |
| --- | ---: | ---: | ---: |
| top-level SweepGA, this run | 221342 | 26037441 | 398 |
| no-guard POASTA explosion | 193318 | 70400000 | 99 |
| guarded chain-POVU | 19433 | 357353 | 265 |
| neighbor-merge | 18761 | 461241 | 72 |
| iterative-multi baseline | 16547 | 396438 | 103 |

## Conclusion

Top-level SweepGA crush is valid and path-preserving on real C4, and the logs
show real SweepGA output for large top-level regions rather than a pure
zero-alignment concatenation path. However, this experiment does not fix the C4
under-alignment/explosion problem. It reduces segment bp compared with the
no-guard POASTA explosion, but it increases segment count and produces far more
trivial-stringy candidates than the useful baselines. The output is also far
larger than the pre-crush no-sort graph (`18048 S / 389354 bp`).

The practical failure mode is mixed: large top-level regions can align and
collapse aggressively, but many shorter top-level regions are below the safe
wfmash aligner-input floor and therefore pass through as validated,
path-preserving stringy replacements. Because quality metrics are diagnostic
only for this task, those replacements are applied instead of falling back.
