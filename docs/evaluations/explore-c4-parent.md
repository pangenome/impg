# Explore C4 Parent-First Large Crush

Task: `explore-c4-parent`

Date: 2026-06-01

## Scope

This experiment reran HPRCv2 C4 on
`GRCh38#0#chr6:31891045-32123783` with `window-mode=outward` so the
multi-bubble resolver would generate parent-sized/outward residual windows
rather than only the 8-10 kb sibling-run windows used by the previous larger
POASTA run.

Input data:

- syng prefix: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng`
- AGC: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc`
- query distance: `-d 50k`
- threads: `-t 32`

Output root:

`/home/erikg/impg/data/explore_c4_parent_20260601T212246Z`

Binary provenance:

- used `/home/erikg/.cargo/bin/impg`
- `impg 0.4.1`
- installed timestamp: `2026-06-01 14:25:06 +0000`
- worktree HEAD and source checkout HEAD:
  `53af68430ef3b5261122e69eddd06245505d3986`
- `cargo build --release` in the isolated worktree was attempted first but
  failed in the vendored `wfmash-rs` build because `htslib/faidx.h` is missing.
  I did not fall back to an unverified stale target binary; the installed
  binary was checked for the multi-bubble/window/replacement flags, and the
  run logs below show `mode=Outward` plus replacement-routing telemetry.

Lineage:

- baseline/current best:
  `/home/erikg/impg/data/c4_multibubble_impl3_autoroute_20260601T1412Z/run.nosort.gfa`
- prior larger POASTA comparator:
  `/home/erikg/impg/data/run_larger_poasta_20260601T1548Z/run.nosort.gfa`
- this run mutates the larger POASTA configuration by switching from
  `window-mode=combined` to `window-mode=outward` and by forcing the SweepGA
  variant with `auto-poasta-max-len=0` plus `replacement-*=1-to-1`.

## Run Outcomes

| run | engine | status | wall time | max RSS |
| --- | --- | --- | ---: | ---: |
| `parent_poasta_outward` | full parent/outward POASTA/auto route | terminated before GFA | 17:56.19 | 74,337,904 KB |
| `parent_sweepga_outward` | forced SweepGA, initial `1:1` spelling | CLI error before crush | 3:19.14 | 35,604,620 KB |
| `parent_sweepga_outward_1to1` | forced SweepGA, `1-to-1`, `auto-poasta-max-len=0` | terminated before GFA | 21:42.93 | 62,911,736 KB |
| `parent_sweepga_outward_single_site` | forced SweepGA, `max-window-sites=1` diagnostic | terminated before GFA | 8:28.84 | 57,889,524 KB |
| `parent_outward_admission_control` | attempted no-build control | terminated before GFA | about 8:33 | about 62,031,496 KB observed |
| `masked_syng_control` | no-crush masked syng control | success | 4:58.83 query | 57,790,516 KB |

Structured run table:

`/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/run-outcomes.tsv`

The initial forced SweepGA spelling used
`replacement-num-mappings=1:1,replacement-scaffold-filter=1:1`. The engine DSL
split on `:`, so it parsed as `replacement-num-mappings=1` and failed. The
successful spelling for reaching replacement routing is `1-to-1`.

`candidate-limit=0` is not a dry-run/no-build setting in this binary. It means
unlimited, so the admission-control attempt entered SweepGA/transclosure on the
first large outward windows and was stopped.

## Admission Evidence

The required parent/outward configurations did reach outward candidate
generation. The full POASTA and forced SweepGA `1-to-1` runs both logged:

- `mode=Outward`
- `replacement-routing=multi-site-auto-poasta-or-sweepga/single-site-auto`
- `2479 POVU site(s), 2475 polymorphic site(s)`
- `72295 generated candidate(s), 32 considered after cap`
- `sources outward-residual-window=72295`
- `total-bp-cap-rejected=0`
- `max-len-cap-rejected=0`
- `median-len-cap-rejected=0`

The considered outward windows were parent-sized:

- largest considered root span: `60682 bp`
- largest considered traversal max length: `126158 bp`
- largest considered median traversal length: `54316 bp`
- largest considered total sequence: `25519322 bp`

The largest discovered residual sites included multiple parent residual chunks
above 30 kb:

- `root_span=32834bp`
- `root_span=32809bp`
- `root_span=32797bp`
- `root_span=32793bp`
- `root_span=32787bp`

The single-site forced-SweepGA diagnostic also admitted a direct >30 kb
candidate:

- candidate `#6`
- `source=outward-residual-window`
- `sites=1`
- `root_span=32787bp`
- `max=32787`
- `median=26418`
- `total=2363564`

This satisfies the admission question: the cap was no longer blocking
parent-sized chunks. The blocker moved downstream into replacement building.

## Replacement Routing

POASTA parent/outward:

- first four builds were routed `with Poasta`
- candidate totals were about 23.8-25.5 Mb each
- no completion was logged after the first four large POASTA replacements
- the run was stopped after 17:56.19 with 74.3 GB max RSS

Forced SweepGA parent/outward:

- first four builds were routed `with Sweepga`
- FastGA/FAtoGDB/GIXmake/ALNtoPAF logs were present
- PAF outputs for the first batch were about 195-269 MB
- transclosure logged `453 edges from 102831 total pairs (227x reduction)`
- the run was stopped after 21:42.93 with 62.9 GB max RSS

Forced SweepGA single-site diagnostic:

- first four single-site outward candidates were routed `with Sweepga`
- transclosure logged `166 edges from 13861 total pairs (83x reduction)`
- the admitted >30 kb residual candidate was present in the considered set,
  but build ordering started with higher step-savings candidates and no GFA was
  produced before the run was stopped.

## Successful Control Output

No parent/outward replacement run produced a GFA. To still complete the
path-validation, graph-report, sort, render, and upload steps on a successful C4
artifact from the same source/query/mask, I produced a no-crush masked syng
control:

`/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa`

Control post-processing:

- graph report:
  `/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/graph-report.tsv`
- sorted GFA:
  `/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.Ygs.gfa`
- PNG:
  `/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/masked-syng-control.Ygs.gfalook-m.png`
- uploaded PNG:
  `erik@hypervolu.me:www/impg/explore-c4-parent-masked-syng-control-20260601T212246Z.png`

Remote upload was verified with `ssh erik@hypervolu.me ls -lh`, showing a 2.8
MB PNG at that path.

## Path Validation

Path validation used the repo's built release example:

`/home/erikg/impg/target/release/examples/compare_gfa_paths`

Comparison:

- expected: current best/start GFA
- observed: this run's masked syng control GFA

Result:

| field | value |
| --- | ---: |
| expected paths | 465 |
| observed paths | 465 |
| missing paths | 0 |
| extra paths | 0 |
| spelling mismatches | 0 |

Validation files:

- `/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/path-validation-vs-start.tsv`
- `/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/path-validation-vs-start.json`

## Metrics

Full metrics comparison:

`/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/metrics-comparison.tsv`

Selected graph-report metrics:

| run | segments | links | path steps | segment bp | ws p99 | ws max | duplicate sequence fraction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| masked syng control | 18,970 | 21,863 | 4,593,198 | 425,742 | 130,878 | 422,205 | 0.468213 |
| current best impl3 autoroute | 17,771 | 20,677 | 4,108,043 | 421,076 | 134,504 | 416,904 | 0.449271 |
| larger POASTA | 18,954 | 22,265 | 3,690,330 | 566,439 | 180,373 | 562,589 | 0.501108 |

The masked syng control is not a candidate improvement. It is included to
record a successful output from the same source/query/mask and to anchor the
path/render validation. The current best remains better than the no-crush
control on segments, links, path steps, and duplicate sequence fraction.

## Recommendation

`window-mode=outward` does what the task needed at admission time. It generates
and considers outward residual windows, including parent-sized residual chunks
above 30 kb. The length and total-sequence caps used here are not the remaining
blocker.

Neither parent/outward POASTA nor forced SweepGA improves C4 graph shape in the
current implementation because neither produced a replacement GFA at this
scale:

- POASTA stalls on several 24-25 Mb multi-parent windows.
- SweepGA reaches alignment but stalls in transclosure with very large PAFs.
- The single-site SweepGA diagnostic proves >30 kb chunks are admitted but also
  shows the current build/order path can spend time on smaller high-step-saving
  windows before the direct parent residual.

Next code change:

Add an explicit admission-only/dry-run mode and a replacement-build budget
boundary for outward residual windows. The boundary should log and skip or split
candidates before expensive builder work when estimated PAF/transclosure cost is
too high. A useful follow-up mode would prioritize isolated parent residual
sites above 30 kb separately from 50-60 kb aggregate outward windows, so the
resolver can test true parent chunks without first attempting four 24-25 Mb
aggregate replacement builds.
