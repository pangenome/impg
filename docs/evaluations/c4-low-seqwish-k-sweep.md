# Full C4 low seqwish-k sweep plus POA polish

Task: `c4-low-seqwish-k-sweep`

Date: 2026-06-05 UTC

Output root:

```text
/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z
```

Driver:

```bash
scripts/c4-low-seqwish-k-sweep.py --out-dir /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z
```

Input FASTA:

```text
/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa
```

## Provenance

All required full-C4 seed floors were built:

- `--min-match-len 79`
- `--min-match-len 127`
- `--min-match-len 191`
- `--min-match-len 255`
- `--min-match-len 311`

The first graph, K=79, ran the full requested FastGA/SweepGA plus seqwish recipe and used a short temp path, `/tmp/c4lk79`, to avoid the `FAtoGDB` buffer overflow seen in an earlier K=311 attempt under `c4_hard_seqwish_k_20260605T135038Z`.

K=79 saved:

```text
/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/debug/k79/raw.paf
```

The remaining seed floors reused that raw PAF with `impg graph --paf-file`, while preserving the same FASTA, `--num-mappings 1:many`, `--scaffold-filter 1:many`, `--scaffold-jump 0`, and seqwish engine settings. This keeps the alignment evidence fixed and changes only the seqwish min exact-run floor. The commands are recorded in:

```text
/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/commands.md
```

## Validation

All required validation criteria were met:

- Full-C4 seed graphs were produced for K=79, 127, 191, 255, and 311.
- Full-C4 seed PNGs and POA-polished PNGs were rendered with `gfasort -p Ygs` and `gfalook -m`.
- All PNGs were uploaded to `https://hypervolu.me/~erik/impg/`.
- For every 1kb POA-polished graph, `compare_gfa_paths seed polished` reported 465 expected paths, 465 observed paths, 0 missing paths, 0 extra paths, and 0 spelling mismatches.
- Optional 2kb POA polish was also run for K=311 and K=255, the two best seed floors by the driver ranking over direct self-loops, adjacent same-node repeats, path whitespace, and K.

Path validation summary:

| graph | expected | observed | missing | extra | spelling mismatches |
| --- | ---: | ---: | ---: | ---: | ---: |
| `c4.k79.poa1kb` | 465 | 465 | 0 | 0 | 0 |
| `c4.k127.poa1kb` | 465 | 465 | 0 | 0 | 0 |
| `c4.k191.poa1kb` | 465 | 465 | 0 | 0 | 0 |
| `c4.k255.poa1kb` | 465 | 465 | 0 | 0 | 0 |
| `c4.k311.poa1kb` | 465 | 465 | 0 | 0 | 0 |
| `c4.k255.poa2kb` | 465 | 465 | 0 | 0 | 0 |
| `c4.k311.poa2kb` | 465 | 465 | 0 | 0 | 0 |

## Scoreboard

Full machine-readable results are in:

```text
/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/scoreboard.tsv
```

Key metrics:

| K | X | segs | links | path steps | seg bp | self-loop edges | adjacent same-node repeats | singleton bp | bp-weighted coverage | ws p99 | ws max | duplicate frac | local repeat ctx | resolved | bailed | build/crush wall | PNG |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 79 | seed | 7275 | 9884 | 3280404 | 237559 | 72 | 314738 | 3359 | 449.107270 | 65 | 34091 | 0.665842 | 1007 |  |  | 5:07.25 | https://hypervolu.me/~erik/impg/c4.k79.seed.png |
| 127 | seed | 7245 | 9847 | 3249356 | 239236 | 72 | 311665 | 3943 | 445.959112 | 37 | 34124 | 0.663216 | 1012 |  |  | 3:41.24 | https://hypervolu.me/~erik/impg/c4.k127.seed.png |
| 191 | seed | 7222 | 9814 | 3201532 | 241441 | 66 | 286005 | 4249 | 441.886316 | 855 | 224363 | 0.660620 | 1018 |  |  | 3:40.80 | https://hypervolu.me/~erik/impg/c4.k191.seed.png |
| 255 | seed | 7180 | 9755 | 3157212 | 246092 | 61 | 271129 | 5539 | 433.534914 | 666 | 228564 | 0.656685 | 1026 |  |  | 3:47.46 | https://hypervolu.me/~erik/impg/c4.k255.seed.png |
| 311 | seed | 7151 | 9716 | 3129613 | 250580 | 60 | 263737 | 5531 | 425.770109 | 767 | 232746 | 0.653195 | 1043 |  |  | 3:09.44 | https://hypervolu.me/~erik/impg/c4.k311.seed.png |
| 79 | 1kb | 8250 | 11025 | 3520043 | 245792 | 67 | 306081 | 4998 | 434.064062 | 114255 | 245235 | 0.727030 | 1122 | 1576 | 0 | 0:31.46 | https://hypervolu.me/~erik/impg/c4.k79.poa1kb.png |
| 127 | 1kb | 8268 | 11045 | 3507551 | 246050 | 67 | 303008 | 5021 | 433.608917 | 114602 | 245456 | 0.727504 | 1116 | 1568 | 0 | 0:31.76 | https://hypervolu.me/~erik/impg/c4.k127.poa1kb.png |
| 191 | 1kb | 8310 | 11095 | 3464424 | 247909 | 62 | 277811 | 5413 | 430.357405 | 115733 | 247351 | 0.727437 | 1115 | 1563 | 0 | 0:31.52 | https://hypervolu.me/~erik/impg/c4.k191.poa1kb.png |
| 255 | 1kb | 8408 | 11167 | 3505268 | 250445 | 59 | 264705 | 6397 | 425.999617 | 115923 | 249854 | 0.729186 | 1111 | 1533 | 0 | 0:31.59 | https://hypervolu.me/~erik/impg/c4.k255.poa1kb.png |
| 311 | 1kb | 8001 | 10749 | 3318893 | 253182 | 57 | 252715 | 6371 | 421.394388 | 121543 | 252604 | 0.714036 | 1113 | 1512 | 0 | 0:29.79 | https://hypervolu.me/~erik/impg/c4.k311.poa1kb.png |
| 311 | 2kb | 8160 | 10928 | 3374388 | 251692 | 57 | 252715 | 6350 | 423.889015 | 120462 | 251134 | 0.718750 | 1112 | 1517 | 0 | 0:33.26 | https://hypervolu.me/~erik/impg/c4.k311.poa2kb.png |
| 255 | 2kb | 8419 | 11195 | 3488995 | 250395 | 59 | 264705 | 6363 | 426.084682 | 115974 | 249837 | 0.728828 | 1114 | 1525 | 0 | 0:32.15 | https://hypervolu.me/~erik/impg/c4.k255.poa2kb.png |

All self-loop length bucket diagnostics are preserved in `scoreboard.tsv`. Every row still has 1bp direct self-loop bucket entries, but the count falls as K increases: seed K=79/127 have 72 direct self-loop edges, seed K=255 has 61, and seed K=311 has 60.

## Interpretation

K=311 is not too large for this full-C4 recipe by the graph-report stress metrics. The lower floors K=79 and K=127 improve seed path white-space p99 and max substantially, but they carry more graph-artifact evidence: more direct self-loop edges, more adjacent same-node path repeats, more path steps, more segments, more duplicate sequence fraction, and higher bp-weighted node coverage. In other words, lower K gives a visually/layout-smoother seed according to white-space p99, but it does not preserve architecture better once self-loop/repeat diagnostics are included.

K=79 and K=127 do avoid the severe K=1 artifact in the earlier baseline in the sense that direct self-loop edges are lower than the K=1 run, but they do not remove the 1bp self-loop class. They still produce 68 one-bp self-loop nodes/edges in the seed graph, compared with 57 at K=311.

The 1kb POA polish preserves spelling for all K values, resolves roughly 1512 to 1576 local candidates with zero bails, and slightly reduces the self-loop/repeat counts. It also makes the path white-space p99 much worse for every K, jumping from tens/hundreds of bp in the seed graphs to roughly 114 to 122 kb after polish. Therefore the 1kb POA result is path-preserving and locally simplifying, but it is not a layout/white-space improvement.

The best pre-polish tradeoff is K=311 if artifact/repeat compactness is weighted ahead of seed white-space p99. It has the fewest segments, links, path steps, direct self-loop edges, adjacent same-node repeats, and duplicate sequence fraction among the required seed floors. If seed white-space p99 alone is the objective, K=127 wins, but that is not the best overall architecture tradeoff.

The best post-1kb POA tradeoff is also K=311. It has the fewest direct self-loop edges, adjacent same-node repeats, segments, path steps, and duplicate sequence fraction among the 1kb polished outputs, with exact path spelling preservation. Its white-space p99 is the worst of the 1kb set by about 5 to 7 kb, so this conclusion depends on prioritizing architecture/repeat diagnostics over layout p99.

The optional 2kb polish did not change the ranking. K=311 2kb slightly improves K=311 1kb white-space p99 and max, and reduces segment bp, but increases segments and path steps. K=255 2kb is essentially similar to K=255 1kb. Neither 2kb output is a decisive improvement over K=311 1kb.

## Answer

K=311 is not too large for C4 in this explicit full-C4 seqwish-k sweep. Smaller floors such as K=79 or K=127 can make the seed render/layout look less jumpy by white-space p99, but they do not preserve the C4 architecture better once direct self-loops, adjacent same-node repeats, duplicate sequence fraction, and compactness are included. The best overall tradeoff before POA is K=311. After 1kb POA, K=311 remains the best artifact/repeat tradeoff and preserves spelling exactly, although POA substantially worsens path white-space for every K.
