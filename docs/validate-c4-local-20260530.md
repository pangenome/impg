# Validate C4 local syng filtering and crush impact

Date: 2026-05-30

Task: `validate-c4-local`

Branch commit tested: `4e067f6` (`wg/agent-328/validate-c4-local`)

Binary: `target/release/impg`, built in this worktree at
2026-05-30 17:33:25 UTC after initializing the `vendor/syng` and
`vendor/gfaffix` submodules. Version output: `impg 0.4.1`.

Output directory:

`/home/erikg/impg/data/validate_c4_local_20260530T173344Z`

## Result

The current default local syng GFA filter did not complete on the C4 locus in a
30 minute validation window. The run reached the new occurrence-level
scaffold-context evaluation after the usual local rebuild and frequency mask,
then remained CPU-bound until interrupted. No GFA, graph report, or PNG was
produced, so there were no new PNGs to upload and the crush run was not started:
the crush command would have to pass through the same default-filter conversion
gate before any crush stage could run.

This is a validation failure for the landed default behavior, but it is also a
clear operational result: on this C4 panel, the new occurrence/scaffold filter is
currently a conversion-time bottleneck before topology metrics can be measured.

## Command

```bash
/usr/bin/time -v target/release/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31982056-32035418 \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o 'gfa:syng-local:blunt,k=63,s=8,seed=7:nosort' \
  -O /home/erikg/impg/data/validate_c4_local_20260530T173344Z/c4-local-syng-k63-default-filter-blunt-only.gfa \
  --render-graph \
  --render-graph-output /home/erikg/impg/data/validate_c4_local_20260530T173344Z/c4-local-syng-k63-default-filter-blunt-only.png \
  --render-graph-depth \
  --render-graph-width 2200 \
  --render-graph-height 1200 \
  --describe-graph \
  --graph-report-output /home/erikg/impg/data/validate_c4_local_20260530T173344Z/c4-local-syng-k63-default-filter-blunt-only.report.tsv \
  --graph-report-format tsv \
  -t 32 -v 1
```

Captured artifacts:

- `/home/erikg/impg/data/validate_c4_local_20260530T173344Z/c4-local-syng-k63-default-filter-blunt-only.stderr`
- `/home/erikg/impg/data/validate_c4_local_20260530T173344Z/c4-local-syng-k63-default-filter-blunt-only.stdout`

Expected but not produced:

- `c4-local-syng-k63-default-filter-blunt-only.gfa`
- `c4-local-syng-k63-default-filter-blunt-only.report.tsv`
- `c4-local-syng-k63-default-filter-blunt-only.png`

## Timings

From the interrupted stderr:

| Stage | Time / status |
| --- | ---: |
| Cold-start elapsed before interrupt | 30:07.82 |
| Syng index load | ~3:06, from 17:34:15 to 17:37:21 UTC |
| Query interval collection | ~1:18, from 17:37:21 to 17:38:39 UTC |
| Local regional syng GBWT rebuild | 0.960 s |
| Temporary FASTA write | 0.049 s |
| Temporary sequence index build | 0.159 s |
| Regional syng index reload | 0.009 s |
| Local syncmer walk collection for masking | 0.074 s |
| Frequency node mask | 3 / 7015 local syncmer nodes in 0.001 s |
| Scaffold-context occurrence evaluation | Did not complete; ran for >25:40 before interrupt |
| Syng2GFA conversion | Not completed |
| Render time | Not reached |

`/usr/bin/time -v` reported 4,078.10 user seconds, 287.02 system seconds,
241% CPU, and 80,416,540 kB maximum RSS. The process was interrupted with
SIGINT after it had emitted:

```text
[syng2gfa] frequency-filtered 3 / 7015 local syncmer node(s) before raw GFA materialization
```

No subsequent `scaffold-context evaluated ...` line was emitted.

## Baseline metrics

Existing baseline reports are under:

`/home/erikg/impg/data/local_syng_guardless_explore_20260530T154204Z`

| Run | segments | total_segment_bp | bp-weighted node coverage | singleton_bp | segment_white_space_bp_fraction | path_white_space_bp_p99 | path_white_space_bp_max | duplicate_sequence_frac | local_repeat_context_nodes | local_repeat_context_occurrences |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| old unmasked local k63 blunt-only | 7,486 | 186,828 | 284.396247 | 37,095 | 0.523969 | 63,632 | 182,857 | 0.405691 | 617 | 763 |
| default node-level mask blunt-only (`cb13376`) | 7,700 | 193,715 | 274.285326 | 43,543 | 0.540225 | 64,672 | 189,067 | 0.435584 | 558 | 686 |
| previous guardless crush control | 8,758 | 213,534 | 248.827737 | 29,174 | 0.548046 | 292 | 71,483 | 0.471683 | 598 | 768 |

The old node-level default mask filtered only 3 / 7015 local syncmer nodes and
did not move C4 toward PGGB-like condensation. The current occurrence-level
default still starts with the same 3 / 7015 node frequency mask, but the new
scaffold-context occurrence pass did not complete, so no current topology
metrics are available to compare against the baseline table.

## Interpretation

This run does not support a claim that the current default filter makes the C4
local syng graph more PGGB-like. It also cannot show whether the final topology
is still underaligned, because final topology was not emitted. The actionable
finding is narrower and stronger: with 466 selected C4 intervals and 2,105,221
local syncmer steps, the occurrence-level scaffold filter is too expensive in
its current default form to reach the reporting/rendering stage.

The failure occurs before crush. Therefore the requested "current best crush"
validation is blocked by the same syng-to-GFA conversion stage, not by crush
itself.

Recommended next step: fix or bound the occurrence-level scaffold-context pass
for high-copy local syng walks, then rerun this exact validation command and the
matching crush command.
