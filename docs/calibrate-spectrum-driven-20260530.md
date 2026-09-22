# Calibrate spectrum-driven syng scaffold-glue mask

Date: 2026-05-30

Task: `calibrate-spectrum-driven`

## Code rule

The previous C4 guard split any shared local syncmer whose local occurrence
count exceeded `ceil(selected_paths * 1.25)`. That bounded the scaffold-context
stage but treated selected path count as a proxy for local copy structure.

The replacement rule records the post-frequency-mask local shared-node spectrum
and selects only dispersed high-copy tail nodes for private splitting. A node is
removed from scaffold support when all of these are true:

- total local occurrences >= `64`;
- occurrence-per-carrying-path ratio >= `2.0`;
- maximum copies on any one carrying path >= `2`;
- maximum path-local span between copies on one path >= `1,000 bp`.

The rule is still provisional, but the default is now tied to local spectra
rather than a fixed selected-path copy factor. Selected nodes are unfolded into
private occurrences; no path sequence is deleted.

## Commands

C4:

```bash
/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31982056-32035418 \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o 'gfa:syng-local:blunt,k=63,s=8,seed=7:nosort' \
  -O /home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/c4_default_spectrum/c4-local-syng-k63-spectrum-default.gfa \
  --render-graph \
  --render-graph-output /home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/c4_default_spectrum/c4-local-syng-k63-spectrum-default.png \
  --render-graph-depth --render-graph-width 2200 --render-graph-height 1200 \
  --describe-graph \
  --graph-report-output /home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/c4_default_spectrum/c4-local-syng-k63-spectrum-default.report.tsv \
  --graph-report-format tsv \
  -t 32 -v 1
```

Non-C4 control:

```bash
/usr/bin/time -v /home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/impg-agent-332-spectrum query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r CHM13#0#chr6:50000000-50100000 \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o 'gfa:syng-local:blunt,k=63,s=8,seed=7:nosort' \
  -O /home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/boring_chr6_spectrum_rerun/boring-chr6-local-syng-k63-spectrum-default.gfa \
  --render-graph \
  --render-graph-output /home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/boring_chr6_spectrum_rerun/boring-chr6-local-syng-k63-spectrum-default.png \
  --render-graph-depth --render-graph-width 2200 --render-graph-height 1200 \
  --describe-graph \
  --graph-report-output /home/erikg/impg/data/calibrate_spectrum_driven_20260530T203220Z/boring_chr6_spectrum_rerun/boring-chr6-local-syng-k63-spectrum-default.report.tsv \
  --graph-report-format tsv \
  -t 32 -v 1
```

The control was rerun from a copied branch binary because the shared
`/home/erikg/impg/target/release/impg` binary was overwritten by another
worktree between validation runs.

## C4 spectrum

Overall local shared-node spectrum:

| Metric | Value |
| --- | ---: |
| Selected paths | 466 |
| Shared nodes | 5,688 |
| Shared-node occurrences | 2,100,031 |
| Max occurrences on one node | 1,214 |
| Max carrying paths on one node | 466 |
| Max copies per path | 5 |
| Max path-local span | 98,328 bp |

Occurrence-per-carrying-path ratio bins:

| Ratio bin | Nodes | Occurrences |
| --- | ---: | ---: |
| <=1.0 | 3,844 | 1,017,602 |
| <=1.25 | 530 | 122,592 |
| <=2.0 | 1,034 | 701,669 |
| <=4.0 | 280 | 258,168 |
| <=8.0 | 0 | 0 |
| >8.0 | 0 | 0 |

Max copies per path bins:

| Max copies/path | Nodes | Occurrences |
| --- | ---: | ---: |
| 1 | 3,844 | 1,017,602 |
| 2 | 456 | 67,265 |
| 3-4 | 1,386 | 1,014,694 |
| 5-8 | 2 | 470 |
| >8 | 0 | 0 |

Path-local dispersion bins:

| Max path span | Nodes | Occurrences |
| --- | ---: | ---: |
| 0 bp | 3,844 | 1,017,602 |
| <1 kb | 2 | 470 |
| 1-10 kb | 29 | 19,840 |
| 10-100 kb | 1,813 | 1,062,119 |
| >100 kb | 0 | 0 |

Threshold candidates:

| Candidate | Nodes split | Occurrences split |
| --- | ---: | ---: |
| Old selected-path factor >1.25 | 1,070 | 897,742 |
| Dispersed ratio >=1.25 | 1,223 | 959,033 |
| Dispersed ratio >=2.0 | 379 | 352,120 |
| Dispersed ratio >=3.0 | 0 | 0 |
| Dispersed ratio >=4.0 | 0 | 0 |

The selected default is `dispersed ratio >=2.0`. It keeps the broad C4
moderate-copy spectrum shared for downstream graph resolution while still
splitting a clearly separated high-copy tail.

## Non-C4 control spectrum

Control locus: `CHM13#0#chr6:50000000-50100000`.

Overall local shared-node spectrum:

| Metric | Value |
| --- | ---: |
| Selected paths | 465 |
| Shared nodes | 6,951 |
| Shared-node occurrences | 1,936,844 |
| Max occurrences on one node | 466 |
| Max carrying paths on one node | 465 |
| Max copies per path | 7 |
| Max path-local span | 115,588 bp |

Occurrence-per-carrying-path ratio bins:

| Ratio bin | Nodes | Occurrences |
| --- | ---: | ---: |
| <=1.0 | 3,071 | 225,468 |
| <=1.25 | 3,873 | 1,711,095 |
| <=2.0 | 7 | 281 |
| <=4.0 | 0 | 0 |
| <=8.0 | 0 | 0 |
| >8.0 | 0 | 0 |

Max copies per path bins:

| Max copies/path | Nodes | Occurrences |
| --- | ---: | ---: |
| 1 | 3,071 | 225,468 |
| 2 | 3,879 | 1,711,111 |
| 3-4 | 0 | 0 |
| 5-8 | 1 | 265 |
| >8 | 0 | 0 |

Path-local dispersion bins:

| Max path span | Nodes | Occurrences |
| --- | ---: | ---: |
| 0 bp | 3,071 | 225,468 |
| <1 kb | 3 | 273 |
| 1-10 kb | 0 | 0 |
| 10-100 kb | 0 | 0 |
| >100 kb | 3,877 | 1,711,103 |

Threshold candidates:

| Candidate | Nodes split | Occurrences split |
| --- | ---: | ---: |
| Old selected-path factor >1.25 | 0 | 0 |
| Dispersed ratio >=1.25 | 0 | 0 |
| Dispersed ratio >=2.0 | 0 | 0 |
| Dispersed ratio >=3.0 | 0 | 0 |
| Dispersed ratio >=4.0 | 0 | 0 |

This control has many nodes with long path-local span, but their
occurrence-per-carrying-path ratios stay near one. The spectrum rule therefore
does not split any additional scaffold glue in this non-C4 local graph.

## Graph metrics

| Run | Segments | Segment bp | Singleton bp | Path whitespace p99 | Path whitespace max | Long links >=1 kb | Render |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| C4 old selected-path factor 1.25 | 910,115 | 24,423,312 | 24,305,964 | 19,345,196 | 24,422,683 | 148,161 | produced |
| C4 spectrum dispersed ratio >=2.0 | 359,381 | 9,688,629 | 9,550,083 | 8,583,211 | 9,647,318 | 272,987 | produced |
| Non-C4 spectrum dispersed ratio >=2.0 | 9,591 | 239,486 | 72,272 | 90,927 | 234,413 | 80,909 | produced |

The old C4 row is the validated artifact from
`docs/fix-c4-syng-20260530.md`:

```text
/home/erikg/impg/data/fix_c4_syng_20260530T1938Z/fixed3/c4-local-syng-k63-default-filter-blunt-only.report.tsv
```

The spectrum rule substantially reduces C4 private clone inflation relative to
the old fixed factor: segments drop by 550,734, stored segment bp drops by
14.7 Mb, and singleton bp drops by 14.8 Mb. The long-link bridge count is worse
on the raw syng-local blunt graph, but the p99 and maximum whitespace lengths
are smaller. Because sequence is preserved as private/unfolded occurrences,
SweepGA/POA/crush remain responsible for resolving repeat interiors rather
than the conversion mask deleting them.

## Recommendation

Use the `>=2.0` dispersed occurrence-per-carrying-path rule as the provisional
default. It is validated on C4 as a large reduction in over-splitting compared
with the old fixed selected-path factor, and on the non-C4 chr6 control as a
zero-additional-split rule despite many long-span but near-single-copy shared
nodes.

This is not a final biological copy-number classifier. It is a conversion-time
scaffold-glue guard calibrated to remove dispersed high-copy local syng tail
nodes while preserving compact or near-single-copy shared structure for
downstream graph resolution.
