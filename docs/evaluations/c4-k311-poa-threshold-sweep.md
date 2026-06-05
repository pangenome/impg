# Full C4 k311 Seed Plus POA Threshold Sweep

Task: `c4-k311-poa-threshold-sweep`
Date: 2026-06-05

## Summary

Ran the requested full-C4 experiment with seqwish induction fixed at
`--min-match-len 311`, then swept direct POA crush traversal ceilings:
`500`, `1000`, `2000`, `5000`, and `10000` bp.

Output directory:

`/home/erikg/impg/data/c4_k311_poa_threshold_20260605T140018Z`

The k311 seed was reused from the concurrent `c4-hard-seqwish-k-sweep`
run after verifying command provenance and path count. Its graph-build log
shows the requested recipe:

```bash
/home/erikg/.cargo/bin/impg graph --sequence-files /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa --gfa-engine seqwish --fastga --num-mappings 1:many --scaffold-filter 1:many --scaffold-jump 0 --min-match-len 311 --temp-dir /tmp/c4k_135300Z/k311 --debug-dir /home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/debug/one_many_minmatch311_scaffold0 -g /home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/graphs/one_many_minmatch311_scaffold0.initial.gfa -t 32 -v 1
```

Seed build timing from that provenance log:

- Wall time: `4:41.36`
- Max RSS: `8,697,416 KB`
- Seed path count: `465`

## Artifacts

- Runner script: `scripts/run-c4-k311-poa-threshold-sweep.sh`
- Reused source seed:
  `/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/graphs/one_many_minmatch311_scaffold0.initial.gfa`
- Copied task seed:
  `/home/erikg/impg/data/c4_k311_poa_threshold_20260605T140018Z/graphs/c4.k311.seed.gfa`
- Raw scoreboard:
  `/home/erikg/impg/data/c4_k311_poa_threshold_20260605T140018Z/scoreboard.tsv`
- Sorted-order scoreboard:
  `/home/erikg/impg/data/c4_k311_poa_threshold_20260605T140018Z/sorted_scoreboard.tsv`
- Commands:
  `/home/erikg/impg/data/c4_k311_poa_threshold_20260605T140018Z/commands.md`
- Runtime summary:
  `/home/erikg/impg/data/c4_k311_poa_threshold_20260605T140018Z/runtime_summary.tsv`
- Crush convergence status:
  `/home/erikg/impg/data/c4_k311_poa_threshold_20260605T140018Z/poa_sweep_status.tsv`
- Upload ledger:
  `/home/erikg/impg/data/c4_k311_poa_threshold_20260605T140018Z/uploaded_urls.tsv`
- Path validation logs:
  `/home/erikg/impg/data/c4_k311_poa_threshold_20260605T140018Z/validation/`

Uploaded PNG URLs:

| Output | PNG |
| --- | --- |
| seed | `http://hypervolu.me/~erik/impg/c4-k311-poa-threshold-c4.k311.seed.Ygs.mean-depth.png` |
| X=500 | `http://hypervolu.me/~erik/impg/c4-k311-poa-threshold-c4.k311.poa500.Ygs.mean-depth.png` |
| X=1000 | `http://hypervolu.me/~erik/impg/c4-k311-poa-threshold-c4.k311.poa1000.Ygs.mean-depth.png` |
| X=2000 | `http://hypervolu.me/~erik/impg/c4-k311-poa-threshold-c4.k311.poa2000.Ygs.mean-depth.png` |
| X=5000 | `http://hypervolu.me/~erik/impg/c4-k311-poa-threshold-c4.k311.poa5000.Ygs.mean-depth.png` |
| X=10000 | `http://hypervolu.me/~erik/impg/c4-k311-poa-threshold-c4.k311.poa10000.Ygs.mean-depth.png` |

## Sweep Result

Every threshold converged in one `--max-iterations 5` crush call. Each run
resolved candidates for four rounds, then the fifth discovery reported no
eligible candidates. There were no bails.

| X | Resolved | Bailed | Rounds | Crush wall | Crush RSS |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 500 | 1,461 | 0 | 4 | 29.55 s | 1,463,628 KB |
| 1000 | 1,512 | 0 | 4 | 34.09 s | 1,599,296 KB |
| 2000 | 1,517 | 0 | 4 | 34.07 s | 1,598,156 KB |
| 5000 | 1,517 | 0 | 4 | 32.47 s | 1,611,552 KB |
| 10000 | 1,518 | 0 | 4 | 35.75 s | 1,625,512 KB |

## Path Validation

All five polished outputs preserve the seed paths exactly:

| X | expected paths | observed paths | missing | extra | spelling mismatches |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 500 | 465 | 465 | 0 | 0 | 0 |
| 1000 | 465 | 465 | 0 | 0 | 0 |
| 2000 | 465 | 465 | 0 | 0 | 0 |
| 5000 | 465 | 465 | 0 | 0 | 0 |
| 10000 | 465 | 465 | 0 | 0 | 0 |

## Sorted-Order Metrics

`path_white_space_*` depends on segment order. The raw `scoreboard.tsv`
contains graph-report metrics on the unsorted POA GFAs, which is useful for
provenance but overstates visual whitespace because replacement nodes are
appended by crush. The table below uses graph reports on the `gfasort -p Ygs`
GFAs, matching the uploaded PNG order.

| X | Segments | Links | Path steps | Segment bp | Self-loop L edges | Adjacent same-node path steps | Singleton bp | BP-weighted node coverage | White p99/max | Duplicate seq frac | Local repeat contexts |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| seed | 7,151 | 9,716 | 3,129,613 | 250,580 | 60 | 263,737 | 5,531 | 425.770109 | 234/36,303 | 0.653195 | 621/1,043 |
| 500 | 8,149 | 10,895 | 3,380,078 | 254,017 | 57 | 252,715 | 6,634 | 420.009188 | 139/25,170 | 0.719597 | 671/1,114 |
| 1000 | 8,001 | 10,749 | 3,318,893 | 253,182 | 57 | 252,715 | 6,371 | 421.394388 | 141/25,170 | 0.714036 | 667/1,113 |
| 2000 | 8,160 | 10,928 | 3,374,388 | 251,692 | 57 | 252,715 | 6,350 | 423.889015 | 139/25,170 | 0.718750 | 667/1,112 |
| 5000 | 8,160 | 10,928 | 3,374,388 | 251,692 | 57 | 252,715 | 6,350 | 423.889015 | 139/25,170 | 0.718750 | 667/1,112 |
| 10000 | 8,230 | 11,023 | 3,369,192 | 258,174 | 57 | 251,738 | 6,355 | 413.246392 | 139/31,652 | 0.718955 | 675/1,125 |

## Interpretation

The hard k311 seed already avoids the worst one-base induction artifact, but
it still has 60 direct self-loop L edges. POA crush only reduces that to 57
at every ceiling tested; it does not eliminate self-loops.

Raising X from 500 bp to 1 kb is useful. It resolves 51 additional candidates,
keeps path spelling exact, reduces singleton bp relative to 500, and produces
fewer segments, links, path steps, and duplicate-sequence nodes than the 500 bp
output.

Raising X from 1 kb to 2 kb has only a small gain: five additional candidates,
slightly better singleton bp and BP-weighted coverage, and the same sorted
white-space max. The trade-off is more segments, links, path steps, and a
higher duplicate sequence fraction than 1 kb.

X=5 kb is metric-identical to X=2 kb in this run, so there is no benefit to
using 5 kb for this seed under the direct POA budget.

X=10 kb resolves one additional candidate beyond 2/5 kb, but it worsens the
larger structural metrics: more segments, links, total segment bp, local repeat
contexts, worse BP-weighted coverage, and worse sorted white-space max. It is
not worth using as the default POA-only ceiling.

## Recommendation

Use `k311 + POA X=1kb` as the conservative default for this hard seed. It gets
nearly all of the direct-bubble cleanup available in this sweep, stays fast,
preserves all 465 paths, and avoids the extra fragmentation seen at larger
ceilings.

Use `X=2kb` only if the priority is the small additional singleton/coverage
improvement over 1 kb and the extra segments/path steps are acceptable. There
is no evidence from this run that `X=5kb` or `X=10kb` is needed.

No source code was modified for this experiment.
