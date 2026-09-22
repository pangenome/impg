# C4 Blocker 03: Full Rerun Scoreboard

Task: `c4-blocker-03-full-rerun-scoreboard`

Branch/commit tested: `wg/agent-422/c4-blocker-03-full-rerun-scoreboard` at
`47fef15eb9c3a47b9b9c19c348a955cfefbfe692`.

Result: **regressed, not fixed**. The full C4 crush completed and preserved all
465 paths, but the final graph moved farther away from PGGB-like quality than
both the integrated C4 baseline and the blocker-02 representative slice.

## Run

Input GFA:

`/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa`

Output directory:

`/home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z`

Command:

```bash
/usr/bin/time -v /home/erikg/.cargo/bin/impg crush \
  --gfa /home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa \
  --output /home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z/full_c4_after_residual_routing.gfa \
  --method iterative-multi-level \
  --window-mode largest \
  --window-target-bp 30k \
  --max-window-sites 8 \
  --candidate-limit 1 \
  --max-iterations 8 \
  --polish-method poasta \
  --polish-rounds 0 \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --threads 32 \
  -v 2
```

Completion:

- Exit status: `0`
- Wall time: `3:05.08`
- Maximum resident set size: `2446448` KB
- Accepted candidates: `8/8`
- Per-round accepted candidates: `[r1=1, r2=1, r3=1, r4=1, r5=1, r6=1, r7=1, r8=1]`

No graph-quality metric was used as an application gate. Metrics below are
diagnostic only; exact path preservation remains the hard acceptance invariant.

## Artifacts

- GFA:
  `/home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z/full_c4_after_residual_routing.gfa`
- Crush stderr log:
  `/home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z/crush.stderr.log`
- Crush stdout log:
  `/home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z/crush.stdout.log`
- Graph-report TSV:
  `/home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z/full_c4_after_residual_routing.graph-report.tsv`
- Path-validation TSV:
  `/home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z/full_c4_after_residual_routing.path-validation.tsv`
- gfalook `-m` PNG:
  `/home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z/full_c4_after_residual_routing.gfalook-m.png`
- Scoreboard TSV:
  `/home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z/c4_blocker_03_scoreboard.tsv`
- Metric deltas TSV:
  `/home/erikg/impg/data/c4_blocker_03_full_rerun_20260603T041449Z/c4_blocker_03_metric_deltas.tsv`

The PNG copy succeeded:

`erik@hypervolu.me:www/impg/c4-blocker-03-full-rerun-scoreboard.gfalook-m.png`

## Path Validation

External validation against the pre-crush C4 input:

```text
expected_paths	465
observed_paths	465
missing_paths	0
extra_paths	0
spelling_mismatches	0
```

Visible path names remain semantic sample/haplotype-compatible names; the run did
not use sidecars to recover names.

## Scoreboard

The PGGB control row uses the existing 465-path control report at
`/home/erikg/impg/data/c4_blocker_01_poasta_scale_20260603T0300Z/before/pggb_control.graph-report.tsv`,
whose source is
`/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.nosort.gfa`.
The blocker-02 row is the requested representative slice, not a full
multi-round C4 graph.

| label | segments | total segment bp | singleton bp | bp-weighted coverage | path whitespace p99 | path whitespace max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PGGB control, 465 paths | 13288 | 234524 | 2890 | 454.919215 | 14 | 219917 |
| Integrated C4 baseline | 19026 | 462940 | 78435 | 230.460695 | 160890 | 459430 |
| Blocker-02 largest1 slice | 19162 | 488855 | 78142 | 218.243598 | 172012 | 485345 |
| Blocker-03 full rerun | 20674 | 532484 | 101312 | 200.361840 | 196939 | 529093 |

Current full rerun deltas:

| metric | vs integrated baseline | vs blocker-02 slice | current / PGGB |
| --- | ---: | ---: | ---: |
| segments | +8.66% | +7.89% | 1.556x |
| total segment bp | +15.02% | +8.92% | 2.270x |
| singleton bp | +29.17% | +29.65% | 35.060x |
| bp-weighted coverage | -13.06% | -8.19% | 0.440x |
| path whitespace p99 | +22.41% | +14.49% | 14067x |
| path whitespace max | +15.16% | +9.01% | 2.406x |

Decision: the actual full C4 problem is **not fixed**. It is improved only in
the operational sense that the full run now completes quickly enough with broad
residual SweepGA routing and exact path preservation. Quality is **regressed**
relative to the integrated C4 baseline.

## Diagnosis

The main failure mode is repeated broad residual churn around the same root
coordinate interval after successful replacements. The run no longer stalls on
the blocker-02 small-POA case, but it keeps rebuilding near-identical broad
top-level sites:

| round | method | traversals | median len | total len | root span | interval |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| 1 | SweepGA | 117 | 26418 | 2363564 | 32787 | `>272218963..>272219340` |
| 2 | Poasta | 143 | 38 | 348244 | 32776 | `>272218968..>272219341` |
| 3 | SweepGA | 256 | 26523 | 7245209 | 32891 | `>272218974..>272219346` |
| 4 | SweepGA | 269 | 26631 | 7243526 | 32999 | `>272218973..>272218980` |
| 5 | SweepGA | 388 | 26752 | 7301649 | 33120 | `>272218972..>272218984` |
| 6 | SweepGA | 388 | 26800 | 7320273 | 33168 | `>272218971..>272218985` |
| 7 | SweepGA | 388 | 26884 | 7352865 | 33252 | `>272218970..>272218986` |
| 8 | SweepGA | 247 | 26937 | 7319463 | 33304 | `>272218968..>272218987` |

Rounds 5-7 are especially diagnostic: the same 388-traversal broad residual is
rebuilt with nearly identical root spans and replacement size (`897` segments,
about `33136` bp). The exact-path gate is working, but the frontier keeps
choosing structurally overlapping residuals instead of converging toward a
PGGB-like graph.

## Next Blocker

Created exactly one next WG blocker task:

`c4-blocker-04-stop` - **C4 blocker 04: stop repeated broad residual churn**

The blocker asks for structural residual-novelty diagnostics and a frontier
selection change so a just-resolved broad residual cannot be immediately
reselected as a near-identical top-level target unless it is structurally novel
or explicitly expanded/settled as one residual. The task explicitly forbids
graph-quality metric application gates, SweepGA filtering duplication, sidecar
name recovery, and synthetic visible IDs.
