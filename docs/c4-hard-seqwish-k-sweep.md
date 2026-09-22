# Full C4 hard seqwish min-match sweep before 1 kb POA polish

Task: `c4-hard-seqwish-k-sweep`

Date: 2026-06-05

## Summary

Ran the requested full-C4 SweepGA/FastGA + seqwish seed sweep at
`--min-match-len` 311, 511, and 1001 with the same broad filtering as the
liked graph:

```bash
--gfa-engine seqwish --fastga --num-mappings 1:many --scaffold-filter 1:many --scaffold-jump 0 -t 32 -v 1
```

Then ran the requested POA-only 1 kb polish:

```bash
impg crush -g SEED -o POLISHED --method poa --max-iterations 5 \
  --max-traversal-len 1k --max-median-traversal-len 1k \
  --max-total-sequence 1m --max-traversals 10k -t 32 -v 1
```

The raw FastGA PAF was generated once during the k=311 build and saved at:

`/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/debug/one_many_minmatch311_scaffold0/raw.paf`

The k=511 and k=1001 builds reused that raw PAF through `--paf-file`; each
re-applied the same 1:many filtering and kept the same 477,003 of 600,404
mappings. This keeps alignment evidence constant and isolates the seqwish
min-match floor.

Output directory:

`/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z`

Full generated scoreboard:

`/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/metrics.scoreboard.tsv`

## Verdict

Harder seqwish floors reduce the 1 bp repeat-unit self-loop artifact at the
source, but they do not remove it. The baseline min-match-1 seed has 152 direct
self-loop edges, mostly 1 bp nodes; k=311 reduces that to 60, k=511 to 55, and
k=1001 to 50. Adjacent same-node repeat steps likewise drop from 617,670 to
263,737 / 243,911 / 196,423.

k=1001 underaligns too much for this C4 graph. It increases seed singleton bp
to 40,603, drops bp-weighted node coverage to 283.39, and inflates total
segment bp to 376,474. After POA, it remains path-preserving but has 41,008
singleton bp and path white-space p99 of 140,823 bp.

The best floor among this sweep is k=311. k=511 buys only five fewer seed
self-loop edges than k=311, but doubles singleton bp and did not exhaust the
1 kb POA frontier after three complete 5-round passes. k=1001 is a clear
underalignment point.

Best current recipe: k=311 hard seqwish floor plus the 1 kb POA polish. Use
explicit `normalize-self-loops` only as an optional final reporting/consumer
variant when loop-free repeat runs are required. Normalization removes all
direct self-loop and adjacent same-node repeat steps, but it worsens
order/layout-sensitive white-space metrics and should not be conflated with
seed quality.

## Scoreboard excerpt

The full scoreboard TSV includes segments, links, paths, path steps, segment
bp, direct self-loop edges, adjacent same-node/same-step repeats, singleton bp,
bp-weighted coverage, path white-space p99/max, duplicate sequence fraction,
local repeat contexts, POVU sites, path validation, and wall/RSS timing for
graph build, POA, sort, and render.

| Graph | Segments | Links | Paths | Path steps | Segment bp | Direct self-loops | Adjacent repeat steps | Singleton bp | bp-weighted coverage | White-space p99 / max | Duplicate frac | Local repeat contexts | Build / POA wall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| min-match-1 seed baseline | 8,186 | 11,303 | 465 | 4,632,725 | 232,093 | 152 | 617,670 | 2,808 | 459.68 | 189,754 / 217,619 | 0.7046 | 572 / 997 | existing |
| min-match-1 POA 1 kb baseline | 9,079 | 12,339 | 465 | 4,843,152 | 240,589 | 144 | 604,675 | 4,055 | 443.45 | 100,506 / 240,054 | 0.7507 | 657 / 1,106 | existing |
| min-match-1 POA normalized | 9,416 | 12,994 | 465 | 4,238,477 | 243,782 | 0 | 0 | 4,490 | 437.64 | 200,983 / 240,054 | 0.7589 | 644 / 1,088 | existing |
| k=311 seed | 7,151 | 9,716 | 465 | 3,129,613 | 250,580 | 60 | 263,737 | 5,531 | 425.77 | 858 / 232,746 | 0.6532 | 621 / 1,043 | 281.36 s |
| k=311 POA 1 kb | 8,001 | 10,749 | 465 | 3,318,893 | 253,182 | 57 | 252,715 | 6,371 | 421.39 | 121,543 / 252,604 | 0.7140 | 667 / 1,113 | 281.36 / 30.49 s |
| k=311 POA normalized | 8,191 | 11,110 | 465 | 3,066,178 | 255,657 | 0 | 0 | 6,769 | 417.31 | 187,958 / 252,604 | 0.7194 | 651 / 1,097 | +2.88 s norm |
| k=511 seed | 7,067 | 9,612 | 465 | 3,011,518 | 275,754 | 55 | 243,911 | 10,956 | 386.90 | 402 / 44,466 | 0.6411 | 668 / 1,092 | 194.06 s |
| k=511 POA 1 kb | 9,195 | 12,135 | 465 | 3,371,167 | 296,607 | 49 | 216,066 | 14,275 | 359.70 | 133,723 / 294,781 | 0.7302 | 770 / 1,243 | 194.06 / 99.19 s |
| k=511 POA normalized | 9,351 | 12,429 | 465 | 3,155,101 | 298,719 | 0 | 0 | 14,642 | 357.16 | 202,604 / 294,781 | 0.7337 | 758 / 1,234 | +2.93 s norm |
| k=1001 seed | 6,919 | 9,421 | 465 | 2,749,461 | 376,474 | 50 | 196,423 | 40,603 | 283.39 | 921 / 73,158 | 0.6195 | 736 / 1,161 | 162.14 s |
| k=1001 POA 1 kb | 9,799 | 12,540 | 465 | 3,214,924 | 386,628 | 44 | 169,188 | 41,008 | 275.95 | 140,823 / 383,602 | 0.7462 | 758 / 1,200 | 162.14 / 58.46 s |
| k=1001 POA normalized | 9,932 | 12,793 | 465 | 3,045,736 | 388,500 | 0 | 0 | 41,380 | 274.62 | 248,734 / 383,602 | 0.7483 | 753 / 1,197 | +2.94 s norm |

## POA convergence

| Variant | POA passes | Converged? | Note |
| --- | ---: | --- | --- |
| k=311 | 1 | yes | Round 5 reported no eligible candidates after 1,512 resolved candidates. |
| k=511 | 3 | no | Continued for three complete 5-round passes. The third pass still selected 2 candidates in round 5 after resolving 1,515 more candidates in that pass. I stopped there because the frontier was not simply draining; each continuation re-opened hundreds of candidates, and the graph metrics were already worse than k=311. |
| k=1001 | 2 | yes | First pass hit the 5-round limit with 1 selected candidate; second pass reached no eligible candidates in round 5 after 1,089 more resolved candidates. |

## Path validation

Every POA-polished output and every optional self-loop-normalized output
preserved all 465 seed paths exactly:

```text
expected_paths  465
observed_paths  465
missing_paths   0
extra_paths     0
spelling_mismatches  0
```

Validation logs are in:

`/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/validation/`

## PNG URLs

All URLs below returned HTTP 200 on 2026-06-05.

| Variant | PNG URL |
| --- | --- |
| k=311 seed | http://hypervolu.me/~erik/impg/c4-hard-seqwish-one_many_minmatch311_scaffold0.initial.Ygs.mean-depth.png |
| k=311 POA 1 kb | http://hypervolu.me/~erik/impg/c4-hard-seqwish-one_many_minmatch311_scaffold0.poa1kb.Ygs.mean-depth.png |
| k=311 POA normalized | http://hypervolu.me/~erik/impg/c4-hard-seqwish-one_many_minmatch311_scaffold0.poa1kb.selfloop-normalized.Ygs.mean-depth.png |
| k=511 seed | http://hypervolu.me/~erik/impg/c4-hard-seqwish-one_many_minmatch511_scaffold0.initial.Ygs.mean-depth.png |
| k=511 POA 1 kb | http://hypervolu.me/~erik/impg/c4-hard-seqwish-one_many_minmatch511_scaffold0.poa1kb.Ygs.mean-depth.png |
| k=511 POA normalized | http://hypervolu.me/~erik/impg/c4-hard-seqwish-one_many_minmatch511_scaffold0.poa1kb.selfloop-normalized.Ygs.mean-depth.png |
| k=1001 seed | http://hypervolu.me/~erik/impg/c4-hard-seqwish-one_many_minmatch1001_scaffold0.initial.Ygs.mean-depth.png |
| k=1001 POA 1 kb | http://hypervolu.me/~erik/impg/c4-hard-seqwish-one_many_minmatch1001_scaffold0.poa1kb.Ygs.mean-depth.png |
| k=1001 POA normalized | http://hypervolu.me/~erik/impg/c4-hard-seqwish-one_many_minmatch1001_scaffold0.poa1kb.selfloop-normalized.Ygs.mean-depth.png |

## Artifacts

- Commands: `/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/commands.md`
- Runtime summary: `/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/runtime_summary.tsv`
- Scoreboard: `/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/metrics.scoreboard.tsv`
- Continuation notes: `/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/continuation_notes.tsv`
- Reports: `/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/reports/`
- Rendered PNGs: `/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/renders/`
- Path validation: `/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/validation/`

## Interpretation

The high floor does address the suspected seqwish source mechanism, but only
partially. The residual hard-floor loops remain mostly 1 bp repeat-unit loops:

- k=311 seed: 57 one-base self-loop nodes, plus one 2 bp, one 3 bp, and one
  111 bp loop node.
- k=511 seed: 52 one-base self-loop nodes, plus one 2 bp, one 3 bp, and one
  111 bp loop node.
- k=1001 seed: 47 one-base self-loop nodes, plus two 3 bp and one 111 bp loop
  node.

The source-level loop reduction from k=311 to k=511 and k=1001 is modest
relative to the cost. k=311 is the useful inflection point: it cuts seed
self-loop edges by 60.5% versus min-match 1 while keeping singleton bp in the
same order of magnitude as the previous seed family. k=511 and k=1001 continue
lowering loop counts, but the singleton/coverage cost accelerates and POA
frontier behavior becomes less stable.

So the C4 answer is not "make seqwish k about 1 kb." A floor around 311 is the
best tested value here; if the residual short-unit loops are unacceptable, use
the explicit path-preserving self-loop normalizer as a reporting/final-output
option rather than raising the seqwish floor to 1 kb.
