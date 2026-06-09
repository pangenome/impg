# Fuller C4 localized polish after exact-run fix

Date: 2026-06-09

Task: `run-fuller-c4`

Generated artifacts are outside git under:

```text
/home/erikg/impg/data/run_fuller_c4_20260609T160000Z
```

The run used branch `wg/agent-688/run-fuller-c4` at `08f12cf`, which includes
the `diagnose-exact-c4` exact-run scan fix. I rebuilt
`target/release/impg` from the worktree after initializing `vendor/syng` and
`vendor/gfaffix`.

## Verdict

Expanded localized-polish budgets do not resolve the remaining C4 open-bubble
signature. They do apply more exact-path-preserving replacements and reduce
some repeat/self-loop diagnostics, but the final candidate-chunk count stalls at
about 21.9k and graph-shape metrics move in mixed or negative directions:
segments, segment bp, singleton bp, segment white-space, and replay/compression
worsen as the chunk budget increases.

No exact path corruption was observed. Exact path preservation remained the only
hard rejection gate in the localized-polish logs:

```text
hard_gate=exact_path_corruption, metrics_gate=false
```

Both expanded runs ended by chunk budget, not convergence.

## Render uploads

Both URLs were verified with `curl -I` as `HTTP/1.1 200 OK` and
`Content-Type: image/png`.

| run | PNG URL | HTTP content length |
|---|---|---:|
| moderate, 2 iterations, 32 chunks/iteration, 64 total chunks | http://hypervolu.me/~erik/impg/c4.fuller-localized.moderate-2x32-64.Ygs.gfalook-m.png | 554,210 |
| strong, 2 iterations, 128 chunks/iteration, 256 total chunks | http://hypervolu.me/~erik/impg/c4.fuller-localized.strong-2x128-256.Ygs.gfalook-m.png | 714,916 |

Only these PNGs were uploaded to `erik@hypervolu.me:www/impg/`.

## Configurations

Common C4 query settings matched the diagnostic tuned run:

```text
-a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng
-b /home/erikg/impg/data/run_fuller_c4_20260609T160000Z/c4.bed
-d 100k
--sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc
--num-mappings 1:many
--scaffold-filter 1:many
--scaffold-jump 0
--min-match-len 311
-o gfa:syng-local:localized,...:nosort
-t 16 -v 1
```

Localized-polish configs:

| label | localized config |
|---|---|
| diagnostic one-chunk baseline | `iterations=1,max-chunks-per-iteration=1,max-total-chunks=1,max-total-bp=200k,max-runtime-secs=300,flank=1k,method=poa,polish-rounds=0` |
| moderate | `iterations=2,max-chunks-per-iteration=32,max-total-chunks=64,max-total-bp=2m,max-runtime-secs=1800,flank=1k,method=poa,polish-rounds=0` |
| strong | `iterations=2,max-chunks-per-iteration=128,max-total-chunks=256,max-total-bp=8m,max-runtime-secs=3600,flank=1k,method=poa,polish-rounds=0` |

One caveat: the moderate run collected 53,133,182 bp of spelled paths, while
the diagnostic baseline and strong run collected 53,165,920 bp. Exact validation
was still against each run's own collected local sequences.

## Runtime

`total wall` and `max RSS` are from `/usr/bin/time -v` around the full query.
`post-index query wall` is the `Syng query complete ... in` log value, which
excludes syng index loading. `localized elapsed` is the localized-polish loop
time after local seed induction.

| run | exit | total wall | post-index query wall | local seed elapsed | localized elapsed | max RSS KB |
|---|---:|---:|---:|---:|---:|---:|
| diagnostic one-chunk baseline | 0 | 11:00.11 | 7:29.598 | 300.730s | 67.907s | 49,943,524 |
| moderate 2x32/64 | 0 | 13:52.12 | 10:06.947 | 214.749s | 315.521s | 49,864,908 |
| strong 2x128/256 | 0 | 23:24.41 | 19:59.315 | 217.634s | 905.025s | 49,960,968 |

The exact-run scan fix held in both expanded runs: after PAF filtering, logs
handed the filtered PAF directly to seqwish without the prior expensive exact
CIGAR/base replay. Seqwish induction remained about 152-155s in the expanded
runs.

## Dirty-region movement

The table uses localized-polish diagnostic counters from stderr and
`localized_polish_summary.tsv`. `after candidates` is the post-iteration
candidate-chunk count from the final after-scan.

| run | initial dirty sites | initial candidates | selected chunks | selected bp | applied chunks | skipped/non-applied chunks | after candidates | final status | exact paths |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| diagnostic one-chunk baseline | 425,704 | 23,906 | 1 | 2,024 | 1 | 0 | 23,551 | `chunk-budget-exhausted` | pass |
| moderate 2x32/64 | 425,418 | 23,890 | 64 | 154,803 | 20 | 44 | 21,951 | `chunk-budget-exhausted` | pass |
| strong 2x128/256 | 425,704 | 23,906 | 256 | 621,928 | 75 | 181 | 21,951 | `chunk-budget-exhausted` | pass |

Observations:

- The moderate pass reduced candidates from 23,890 to 21,951.
- The stronger pass reached the same final candidate count, 21,951, despite
  selecting 4x as many chunks and applying 75 replacements.
- Many non-applied chunks were skipped because shared flank anchors were
  ambiguous or missing. No selected chunk reported path-invalid.
- Dirty sites fell after the first expanded iteration, but final post-run dirty
  site count is not logged by the current after-summary. Candidate chunks are
  the available final dirty-region count.

## Ygs graph metrics

I ran `gfasort -p Ygs` for both expanded outputs, reused the diagnostic
baseline's saved Ygs GFA, and then ran `target/release/impg graph-report
--format tsv --povu` on all three Ygs GFAs.

| metric | diagnostic one-chunk | moderate 2x32/64 | strong 2x128/256 |
|---|---:|---:|---:|
| segments | 3,851 | 4,000 | 4,158 |
| links | 5,220 | 5,391 | 5,560 |
| paths | 466 | 466 | 466 |
| path steps | 1,784,781 | 1,750,053 | 1,753,671 |
| total segment bp | 97,960 | 103,754 | 110,122 |
| singleton bp | 2,531 | 3,055 | 4,863 |
| bp-weighted path-depth / replay compression | 542.730911 | 512.107312 | 482.791086 |
| path-depth median | 457 | 455 | 455 |
| path-depth p95 | 933 | 933 | 933 |
| path-depth max | 16,075 | 16,075 | 16,075 |
| path white-space p99 bp | 210 | 262 | 262 |
| path white-space max bp | 26,933 | 29,096 | 23,107 |
| white-space bridges >= 1 kb | 5,987 | 6,908 | 5,528 |
| segment white-space total bp | 5,656,162 | 8,333,591 | 11,269,943 |
| direct self-loop edges | 14 | 11 | 10 |
| adjacent same-node path steps | 44,368 | 24,642 | 24,606 |
| self-loop repeat runs | 9,125 | 5,899 | 5,881 |
| POVU sites | 1,510 | 1,501 | 1,505 |
| POVU leaf sites | 1,452 | 1,426 | 1,414 |

Quality movement:

- The expanded budgets improve the small-loop/repeat counters substantially:
  direct self-loop edges go 14 -> 11 -> 10, and adjacent same-node path steps
  fall from 44,368 to about 24.6k.
- The path-step count improves relative to the one-chunk baseline, but most of
  that improvement arrives in the moderate run; the stronger run does not
  materially improve it further.
- Graph condensation gets worse by the replay/compression proxy: 542.73 ->
  512.11 -> 482.79. Segment count, segment bp, singleton bp, and segment
  white-space all increase with budget.
- The stronger run reduces long white-space bridge count and max path
  white-space relative to moderate, but the p99 remains worse than the
  one-chunk baseline and segment white-space is highest in the strong run.

## Conclusion

Another pass/chunk budget does not currently resolve the remaining visually
trivial C4 bubbles. The larger budget applies more exact-path-preserving local
SPOA replacements, but it mostly trades some self-loop/repeat cleanup for higher
segment count, lower replay compression, more singleton sequence, and more
segment white-space. The persistent candidate count and frequent ambiguous-flank
skips point to selection/resolver boundary semantics, not simply insufficient
chunk budget.

## Artifact policy

Generated GFA, sorted GFA, PNG, report TSV, debug, and log artifacts remain
outside the repository under:

```text
/home/erikg/impg/data/run_fuller_c4_20260609T160000Z
```

This commit includes only this report.
