# compare-c4-cap64

Date: 2026-06-04

Task: compare the existing C4 cap64 `iterative-multi-level` run against a variant
that disables the POASTA middle tier by routing median traversal length `<15k`
through global POA/SPOA and `>=15k` through SweepGA.

## Summary

The completed baseline POASTA-middle-tier run finished in `4:10.66` wall time
with `4,712,032 KB` max RSS.

The first run using the suggested equal thresholds
`--auto-spoa-max-traversal-len 15k --auto-poasta-max-traversal-len 15k`
finished in `4:08.07` wall time with `4,756,200 KB` max RSS, but it did not
actually disable the POASTA middle tier for multi-site windows in the current
routing code. The log still contains `20` `with Poasta` candidate builds and
`23` `method=Poasta` records. This run is retained as a diagnostic control and
has GFA/sorted/report/render artifacts.

After applying the routing guard that treats equal/inverted SPOA and POASTA
ceilings as an intentional middle-tier collapse, the corrected global-POA run
confirmed the requested routing but timed out after `30:01.66` with status
`124` and `17,545,328 KB` max RSS. It reached only build progress `44/64` and
did not write the final GFA, so no sorted GFA, graph-report TSV, or corrected
render exists for that timed-out run.

## Runs

### Baseline: SPOA 2k, POASTA 15k

Directory:
`data/c4_cap64_diag_20260604T154104Z`

Command shape:

```bash
target/release/impg crush \
  --gfa data/c4_grouped_multi_bubble_seed_top1_20260604T121002Z/sorted/combined_compact_12000bp_4sites_top1_r1.Ygs.gfa \
  --output data/c4_cap64_diag_20260604T154104Z/graphs/multilevel_cap64_spoa2k_poasta15k_r1.gfa \
  --method iterative-multi-level \
  --window-mode combined \
  --window-target-bp 12000 \
  --max-window-sites 4 \
  --candidate-limit 64 \
  --max-iterations 1 \
  --max-span 0 \
  --max-traversal-len 100k \
  --max-median-traversal-len 50k \
  --max-total-sequence 500m \
  --max-traversals 100k \
  --auto-spoa-max-traversal-len 2k \
  --auto-poasta-max-traversal-len 15k \
  --poa-scoring 1,4,6,2,26,1 \
  --min-match-length off \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --sweepga-no-filter true \
  --polish-rounds 0 \
  --threads 32 \
  --verbose 1
```

Result:

| metric | value |
| --- | ---: |
| exit status | 0 |
| wall time | 4:10.66 |
| max RSS | 4,712,032 KB |
| built candidates | 63/64 |
| applied candidates | 23 |
| applied source mix | complete-homologous-window=5, top-level=18 |

Graph-report highlights from
`data/c4_cap64_diag_20260604T154104Z/reports/multilevel_cap64_spoa2k_poasta15k_r1.graph-report.tsv`:

| metric | value |
| --- | ---: |
| status | REVIEW |
| segments | 7405 |
| links | 10122 |
| paths | 465 |
| path_steps | 3619762 |
| total_segment_bp | 242307 |
| node_coverage_mean | 488.826739 |
| bp-weighted coverage | 440.307024 |
| singleton_bp | 2928 |
| duplicate_sequence_frac | 0.675219 |
| local_repeat_context_nodes | 568 |
| path_white_space_bp_p99 | 375 |
| path_white_space_bp_max | 40106 |
| path_white_space_bridges_ge_threshold | 5890 |

### Diagnostic: Suggested Equal Thresholds Before Routing Fix

Directory:
`data/c4_cap64_spoa15k_no_poasta_20260604T180927Z`

This used the requested threshold flags:

```bash
--auto-spoa-max-traversal-len 15k --auto-poasta-max-traversal-len 15k
```

but it was run before the source routing guard was restored. In that code path,
multi-site windows whose median length nominally selected `Poa` were still
promoted to the scalable POASTA tier because `auto_poasta_max_traversal_len > 0`.
The run is therefore not a valid global-POA comparison, but it is useful as a
control showing why the code change was needed.

Result:

| metric | value |
| --- | ---: |
| exit status | 0 |
| wall time | 4:08.07 |
| max RSS | 4,756,200 KB |
| built candidates | 63/64 |
| applied candidates | 23 |
| build-method log counts | `41 with Poa`, `20 with Poasta`, `2 with Sweepga` |
| method-record log counts | `59 method=Poa`, `23 method=Poasta`, `4 method=Sweepga` |

Artifacts produced:

| artifact | path |
| --- | --- |
| GFA | `data/c4_cap64_spoa15k_no_poasta_20260604T180927Z/graphs/multilevel_cap64_spoa15k_no_poasta_r1.gfa` |
| sorted GFA | `data/c4_cap64_spoa15k_no_poasta_20260604T180927Z/sorted/multilevel_cap64_spoa15k_no_poasta_r1.Ygs.gfa` |
| graph-report TSV | `data/c4_cap64_spoa15k_no_poasta_20260604T180927Z/reports/multilevel_cap64_spoa15k_no_poasta_r1.graph-report.tsv` |
| render | `data/c4_cap64_spoa15k_no_poasta_20260604T180927Z/renderings/c4-cap64-spoa15k-suggested-r1.Ygs.mean-depth.png` |
| run logs | `data/c4_cap64_spoa15k_no_poasta_20260604T180927Z/logs/` |

Compared with the baseline graph report, core counts match exactly for
segments, links, paths, path steps, total segment bp, coverage, singleton bp,
duplicate sequence fraction, and local repeat contexts. The graph-report
differences are in path jump / whitespace fields:

| metric | baseline | diagnostic |
| --- | ---: | ---: |
| path_jump_p99 | 32 | 23 |
| path_white_space_bp_p99 | 375 | 117 |
| path_white_space_bridges_ge_threshold | 5890 | 6006 |

### Corrected Global POA/SPOA Middle-Tier Collapse

Directory:
`data/c4_cap64_global_spoa15k_20260604T181607Z`

The corrected run used the same seed graph and the requested equal thresholds:

```bash
timeout 1800s target/release/impg crush \
  --gfa data/c4_grouped_multi_bubble_seed_top1_20260604T121002Z/sorted/combined_compact_12000bp_4sites_top1_r1.Ygs.gfa \
  --output data/c4_cap64_global_spoa15k_20260604T181607Z/graphs/multilevel_cap64_global_spoa15k_r1.gfa \
  --method iterative-multi-level \
  --window-mode combined \
  --window-target-bp 12000 \
  --max-window-sites 4 \
  --candidate-limit 64 \
  --max-iterations 1 \
  --max-span 0 \
  --max-traversal-len 100k \
  --max-median-traversal-len 50k \
  --max-total-sequence 500m \
  --max-traversals 100k \
  --auto-spoa-max-traversal-len 15k \
  --auto-poasta-max-traversal-len 15k \
  --poa-scoring 1,4,6,2,26,1 \
  --min-match-length off \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --sweepga-no-filter true \
  --polish-rounds 0 \
  --threads 32 \
  --verbose 1
```

Result:

| metric | value |
| --- | ---: |
| exit status | 124 |
| timeout | 1800s |
| wall time | 30:01.66 |
| max RSS | 17,545,328 KB |
| CPU | 360% |
| final build progress | 44/64 |
| final GFA written | no |

Routing was confirmed:

| evidence | log line |
| --- | --- |
| `<15k` complete-homologous windows use POA | `building candidate 2/64 ... median-len=14839 ... with Poa` |
| `>=15k` windows use SweepGA | `building candidate 4/64 ... median-len=15005 ... with Sweepga` |
| no POASTA records in corrected run | method-record counts: `44 method=Poa`; build-method counts: `46 with Poa`, `2 with Sweepga` |

Bottleneck:

The bottleneck is the global POA build for high-fanout, near-15kb
complete-homologous windows. Examples:

| candidate | method | median length | timing evidence |
| --- | --- | ---: | --- |
| 2/64 | Poa | 14839 | started 18:16:33; built at 18:27:14 |
| 43/64 | Poa | 14988 | started 18:16:59; built at 18:27:46 |
| 46/64 | Poa | 14957 | started 18:27:46; built at 18:38:34 |
| 47/64 | Poa | 14939 | started 18:38:34; still in progress when the timeout killed the run |
| 48/64 | Poa | 14949 | started 18:38:41; still in progress when the timeout killed the run |

Because the corrected run timed out before graph rewrite, no corrected GFA,
sorted GFA, graph-report TSV, full-overlay render, or pathname render could be
produced.

## Code Note

The code change in `src/resolution.rs` makes equal or inverted POA/POASTA
ceilings explicit: when `auto_spoa_max_traversal_len >=
auto_poasta_max_traversal_len`, the multi-site window override no longer
promotes POA-selected windows into POASTA. This lets the command flags express
the requested middle-tier collapse.

The same source edit also adds the missing `ResolutionMethod::Abpoa` arm to the
debug replacement-build label helper so release builds remain exhaustive.

Validation:

```bash
cargo test --lib
cargo build --release --bin impg
cargo install --path .
```

All passed after adding an assertion for the equal-threshold case.

## Conclusion

The POASTA middle tier is necessary for this C4 cap64 workload at a 15kb
boundary. With global POA handling `<15k` complete-homologous windows, the run
exceeded the baseline wall time by more than 7x, used about 3.7x the baseline
max RSS, and still reached only `44/64` built candidates before timeout.
