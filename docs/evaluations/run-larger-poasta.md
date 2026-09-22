# Run larger-POASTA C4 internal crush experiment

Date: 2026-06-01
Task: `run-larger-poasta`

## Inputs

Starting GFA:

`/home/erikg/impg/data/c4_multibubble_impl3_autoroute_20260601T1412Z/run.nosort.gfa`

Successful output directory:

`/home/erikg/impg/data/run_larger_poasta_20260601T1548Z/`

Successful PNG:

https://hypervolu.me/~erik/impg/run-larger-poasta-c4-20260601T1548Z.png

## Prior run confirmation

The existing run log confirms the requested skipped-candidate counters and run
settings:

| field | value |
| --- | ---: |
| `total-bp-cap-rejected` | 4,524 |
| `max-len-cap-rejected` | 2,254 |
| `median-len-cap-rejected` | 10,023 |
| `overlap_deferred` | 170 |
| `polish-rounds` | 0 |
| `max-rounds` | 1 |

The prior command used the current multi-bubble engine string shape:

`gfa:syng:mask,top=0.001,freq-run=10,freq-span=1000:crush,method=multi-bubble,window-mode=combined,max-window-sites=8,window-target-bp=30k,candidate-limit=192,max-rounds=1,min-match-length=off,max-pair-alignments=0,max-replacement-paf-bytes=0,sweepga-no-filter=true,polish-rounds=0:nosort`

## Successful run

Binary:

`/home/erikg/.cargo/bin/impg` (`impg 0.4.1`, installed at
2026-06-01 14:25 UTC by the C4 multi-bubble implementation work).

Engine:

`gfa:syng:mask,top=0.001,freq-run=10,freq-span=1000:crush,method=multi-bubble,window-mode=combined,max-window-sites=8,window-target-bp=30k,candidate-limit=192,max-rounds=2,min-match-length=off,max-traversal-len=50k,max-median-traversal-len=10k,max-total-sequence=10m,auto-2tier=true,auto-poasta-max-len=50k,max-pair-alignments=0,max-replacement-paf-bytes=0,sweepga-no-filter=true,polish-method=poasta,polish-rounds=until-done,polish-max-traversal-len=50k,polish-max-median-traversal-len=10k,polish-max-total-sequence=10m:nosort`

The successful log shows the alias was accepted and normalized internally as
`method=iterative-multi-level`, with
`replacement-routing=multi-site-auto-poasta-or-sweepga/single-site-auto`. Both
rounds logged `replacement build concurrency=4 (global threads=32,
candidates=192)`. The successful log contains no `FastGA`, `GIXmake`,
`FAtoGDB`, or `SweepGA` lines. Because admitted candidates had max median
lengths below 10 kb and `auto-poasta-max-len=50k`, the built candidates routed
through POASTA rather than SweepGA.

## Wrong-binary attempts

Two wrong-binary attempts were kept under the data directory for provenance:

| directory | binary | result |
| --- | --- | --- |
| `run_larger_poasta_20260601T1535Z` | `/home/erikg/impg/target/release/impg` | Failed after syng index load: `method=multi-bubble` alias unsupported. |
| `run_larger_poasta_20260601T1541Z` | `/home/erikg/impg/target/release/impg` | Relaunched with `method=iterative-multi-level`, then stopped because the log showed `replacement-routing=multi-site-sweepga/single-site-auto` and nested `FastGA`/`GIXmake -T32` work. |

This confirmed that `/home/erikg/impg/target/release/impg` predates the
multi-bubble alias and the bounded auto-POASTA replacement routing needed for
this experiment.

## Wall time and memory

| stage | wall | max RSS |
| --- | ---: | ---: |
| Full query, including syng index load | 13:40.32 | 59,124,608 KB |
| Syng query body reported by impg | 10:12.316 | included above |
| Initial GBWT load | about 3:10 | included above |
| AGC interval collection before blunt render | about 1:20 | included above |
| Exact blunt syng GFA render | 16.092 s | included above |
| Crush round 1 total | 264.45 s | included above |
| Crush round 1 discovery | 23.29 s | included above |
| Crush round 1 replacement build | 238.55 s | included above |
| Crush round 1 rewrite and validate | 2.04 s | included above |
| Crush round 2 total | 251.14 s | included above |
| Crush round 2 discovery | 23.12 s | included above |
| Crush round 2 replacement build | 225.55 s | included above |
| Crush round 2 rewrite and validate | 1.87 s | included above |
| `impg graph-report` | 0:08.68 | 1,041,840 KB |
| `gfasort -p Ygs` | 0:46.56 | 25,740,232 KB |
| `gfalook -m` | 0:02.04 | 184,320 KB |

## Crush results

| round | generated | considered | applied | overlap deferred | failed/empty | total cap rejects | max-len rejects | median-len rejects |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 17,464 | 192 | 15 | 177 | 0 | 67 | 193 | 1,682 |
| 2 | 16,124 | 192 | 11 | 181 | 0 | 21 | 185 | 2,170 |

Per-round accepted candidates: `[r1=15, r2=11]`; total resolved: 26; total
candidates seen: 384.

The raised POASTA budget reduced length-cap skipping substantially relative to
the starting run, but it produced worse final graph shape. Round 1 especially
expanded the graph before round 2 partially recovered whitespace.

| metric | starting run | larger-POASTA final |
| --- | ---: | ---: |
| segments | 17,771 | 18,954 |
| links | 20,677 | 22,265 |
| paths | 465 | 465 |
| path steps | 4,108,043 | 3,690,330 |
| total segment bp | 421,076 | 566,439 |
| mean node coverage | 231.165551 | 194.699272 |
| bp-weighted node coverage | 253.373439 | 188.351215 |
| singleton bp | 77,622 | 69,405 |
| high-coverage bp | 203,812 | 196,404 |
| path whitespace p99 | 134,504 | 180,373 |
| path whitespace max | 416,904 | 562,589 |
| whitespace bridges >= threshold | 163,077 | 153,037 |
| duplicate sequence fraction | 0.449271 | 0.501108 |

The internal crush quality log shows the same pattern:

| point | score | segments | segment bp | ws p99 | ws max |
| --- | ---: | ---: | ---: | ---: | ---: |
| before round 1 | 160,294,797 | 18,970 | 425,742 | 130,878 | 422,205 |
| after round 1 | 229,721,656 | 19,649 | 544,514 | 198,321 | 540,634 |
| after round 2 | 204,829,447 | 18,954 | 566,439 | 174,454 | 562,594 |

## Path validation

Exact path-spelling validation against the starting GFA passed:

| field | value |
| --- | ---: |
| expected paths | 465 |
| observed paths | 465 |
| missing paths | 0 |
| extra paths | 0 |
| spelling mismatches | 0 |

Validation files:

- `/home/erikg/impg/data/run_larger_poasta_20260601T1548Z/path-validation-vs-start.tsv`
- `/home/erikg/impg/data/run_larger_poasta_20260601T1548Z/path-validation-vs-start.json`

## Conclusion

The larger internal POASTA budget did what it was asked to do mechanically:
more candidates passed length caps, two rounds completed with bounded
replacement build concurrency, no SweepGA activity was observed in the
successful run, and exact path spellings were preserved.

It is not an improvement for this C4 topology. Compared with the previous
single-round run, the larger-POASTA final graph has more segments, more links,
more stored segment bp, worse duplicate-sequence fraction, and much worse p99
and max path whitespace. It does reduce path steps, singleton bp, and the count
of long whitespace bridges, but the total shape regression is large enough that
these parameters should not be promoted as the next default.

Recommended next experiment: keep the installed/current binary and the bounded
replacement pool, but either restore the tighter `max-median-traversal-len=1k`
for multi-site sibling windows or add a round-level shape acceptance gate before
allowing raised-POASTA candidates to commit.
