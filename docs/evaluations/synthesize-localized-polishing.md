# Synthesis: Localized Polishing Results And Decision

Date: 2026-06-08
Task: `synthesize-localized-polishing`

## Decision

**Iterate. Do not land this as solved, and do not stop/rethink the whole
architecture yet.**

The workgroup implemented the intended production shape for localized polishing
over SYNG-collected local sequence sets, but validation did not show
PGGB/SmoothXG-like quality on C4 or on any non-C4 SVR locus. The integrated C4
production runs timed out before emitting a seed GFA, so the new localized
output could not be exact-path checked and could not be compared against the
PGGB/SmoothXG control as an output graph.

The next decision is therefore not "tune quality metrics harder." The hard
blocker is making `gfa:syng-local:localized` produce an exact-path-preserving
C4-scale seed/output under budget. Graph-quality metrics should remain
diagnostics after that invariant is observable.

## Final Task Graph Order

The final implementation and validation order came from
`quality-pass-for` and was executed as:

1. `quality-pass-for` - commit `98fef15`
2. `audit-smoothxg-pggb` - commit `a457ee5`
3. `design-localized-smoothxg` - commit `dda9d3f`
4. `implement-explicit-local` - commit `8096460`
5. `implement-graph-badness` - commit `f011d21`
6. `implement-localized-resolver` - commit `2882dce`
7. `integrate-iterative-localized` - commit `b0912ce`
8. `validate-localized-polishing` - commit `ee84b87`
9. `synthesize-localized-polishing` - this report

The audit/design pair ran in parallel. The implementation tasks were serialized
because they touch shared graph output, replacement, resolver, metric, and
path-preservation surfaces.

## What Was Implemented

| Task | Commit | Implemented result |
| --- | --- | --- |
| `quality-pass-for` | `98fef15` | Tightened the workgroup around localized SYNG sequence-set smoothing, exact path preservation as the only hard gate, diagnostic metrics, explicit validation, and no-big-artifact policy. |
| `audit-smoothxg-pggb` | `a457ee5` | Audited existing IMPG PGGB/SmoothXG mechanics and mapped them to `SYNG sequence collection -> region seed induction -> localized dirty-region smoothing -> path-preserving lacing`. |
| `design-localized-smoothxg` | `dda9d3f` | Wrote the implementable design for `syng-local:localized`, including seed options, dirty detectors, chunk selection, resolver tiers, iteration rules, and validation plan. |
| `implement-explicit-local` | `8096460` | Added the explicit local seed induction driver in `src/local_seed.rs`, wired local sequence collection and whole-region SweepGA/FastGA plus seqwish-style seed induction through existing graph output conventions, and added path-name/path-sequence tests. |
| `implement-graph-badness` | `f011d21` | Added `src/graph_badness.rs` with diagnostic graph-quality metrics and dirty-region/chunk detection for singleton/low-depth sequence, white-space/underalignment, path jumps, and loop/self-loop signals. |
| `implement-localized-resolver` | `2882dce` | Added localized dirty-chunk resolver primitives in `src/resolution.rs`, including flanked replacement/trim-back behavior and exact path reconstruction checks. |
| `integrate-iterative-localized` | `b0912ce` | Added `src/localized_polish.rs` and CLI/parser wiring for `gfa:syng-local:localized,...`, with iteration budgets, chunk budgets, resolver settings, diagnostic debug reports, convergence handling, and path-validation gates. |
| `validate-localized-polishing` | `ee84b87` | Ran real C4 production attempts and non-C4 SVR diagnostics, then committed `docs/evaluations/validate-localized-polishing.md` and `docs/evaluations/validate-localized-polishing-metrics.tsv`. |

Design/audit references:

- `docs/designs/localized-smoothxg-polishing.md`
- `docs/evaluations/audit-smoothxg-pggb.md`
- `docs/evaluations/quality-pass-for-localized-smoothing.md`

Validation artifacts:

- Metrics table: `docs/evaluations/validate-localized-polishing-metrics.tsv`
- Validation report: `docs/evaluations/validate-localized-polishing.md`
- C4/local generated artifact root:
  `/home/erikg/impg/data/validate_localized_polishing_20260608T225223Z`
- Current C4 control artifact root:
  `/home/erikg/impg/data/validate_current_c4_20260608T125053Z`
- PGGB/SmoothXG C4 control GFA:
  `/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.Ygs.gfa`
- C4 uploaded visualizations:
  `http://hypervolu.me/~erik/impg/c4.validate-current-c4.k311.seed.repro.Ygs.png`,
  `http://hypervolu.me/~erik/impg/c4.validate-current-c4.k311.poa1kb.repro.Ygs.png`,
  `http://hypervolu.me/~erik/impg/c4.validate-current-c4.pggb-control.Ygs.png`
- Non-C4 SVR local visualization paths:
  `/home/erikg/impg/data/cosigt_svr_hprcv2_crush_20260521T114027Z/renders_biwfa_1d/AMY1A.1d.png`,
  `/home/erikg/impg/data/cosigt_svr_hprcv2_crush_20260521T114027Z/renders_biwfa_1d/CYP2D7.1d.png`,
  `/home/erikg/impg/data/cosigt_svr_hprcv2_crush_20260521T114027Z/renders_biwfa_1d/HLA-A.1d.png`,
  `/home/erikg/impg/data/cosigt_svr_hprcv2_crush_20260521T114027Z/renders_biwfa_1d/KIR2DL1.1d.png`

## Performance And Quality Change

The implemented path is runnable on synthetic integration coverage and is wired
into the production graph-output syntax, but the real C4 production path did
not finish seed induction:

| Run | Status | Runtime | Max RSS | Last observed stage |
| --- | --- | ---: | ---: | --- |
| C4 `syng-local:localized` default | timeout, no output | `20:04.53` | `47,119,692 kB` | FastGA emitted `541,784` unique raw PAF records, `541,318` inter-sequence records |
| C4 `syng-local:localized` tuned k311/1:many | timeout, no output | `15:04.60` | `46,565,888 kB` | SweepGA/FastGA filtering reduced `541,388` inter-sequence records to `431,600` mappings |

Because no integrated localized GFA was emitted, exact path preservation for the
new output is `not_testable_no_output`. That is the decisive validation result.

The current C4 control that does emit output remains exact-path preserving but
does not move toward the PGGB/SmoothXG quality target:

| C4 graph | Status | Segments | Segment bp | Replay/compression | Singleton bp | Low-depth bp | Path steps | White-space bp total | Exact path status |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| k311 seed | pass | `7,346` | `253,140` | `421.464` | `5,925` | `52,671` | `2,865,437` | `24,352,099` | reference |
| current `poa1kb` | pass | `8,142` | `258,628` | `412.521` | `7,492` | `56,791` | `3,065,010` | `26,794,743` | `465` paths, `0` missing/extra/spelling mismatches |
| PGGB/SmoothXG control | pass | `13,288` | `234,524` | `454.919` | `2,890` | `37,477` | `5,538,879` | `17,165,717` | not name-comparable to SYNG seed |

Metric interpretation:

- Exact path corruption remains the hard invariant. The current `poa1kb`
  control passes it, but the integrated localized path has no output to check.
- The current `poa1kb` control worsens relative to its k311 seed on segments,
  segment bp, replay/compression ratio, singleton bp, low-depth bp, path steps,
  and total white-space proxy.
- The PGGB/SmoothXG control is not perfect: it has more segments, more path
  steps, more white-space bridges by count, and self-loop signals. It is still
  the better C4 quality target overall because it has lower segment bp, much
  lower singleton bp, lower low-depth bp, higher replay/compression, and lower
  total white-space proxy.

This is why the conclusion is not based on a single metric. The current local
control preserves paths but is diagnostically below the PGGB/SmoothXG target,
and the integrated localized route has not reached a comparable output.

## C4 And Non-C4 Goal Check

**C4:** the localized SYNG sequence-set smoothing goal was **not met**. The
new `syng-local:localized` route timed out before a seed GFA, so it did not
produce a C4 graph that could be exact-path checked or compared to the
PGGB/SmoothXG control. The older exact-preserving `poa1kb` local-compression
control is still not PGGB/SmoothXG-like quality.

**Non-C4 SVR loci:** the goal was also **not met** on non-C4 SVR loci. The
validation task ran current-schema diagnostics on existing COSIGT/HGSVC/SVR
GFAs for `AMY1A`, `CYP2D7`, `HLA-A`, and `KIR2DL1`, but these were not
integrated `syng-local:localized` production runs and had no source comparator
for exact-path validation.

| Locus | Segments | Segment bp | Replay/compression | Singleton bp | Low-depth bp | White-space bp total | Self-loop repeat runs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `AMY1A` | `61,605` | `1,154,327` | `236.897` | `242,351` | `580,843` | `294,595,633` | `38` |
| `CYP2D7` | `24,289` | `391,057` | `242.604` | `60,306` | `183,874` | `89,705,883` | `0` |
| `HLA-A` | `74,012` | `1,227,678` | `104.214` | `130,519` | `907,951` | `445,626,460` | `4` |
| `KIR2DL1` | `52,928` | `822,137` | `165.220` | `130,694` | `463,241` | `210,805,297` | `1` |

These non-C4 rows are useful diagnostics: the low-depth and white-space signals
are not C4-only. They are not evidence that localized SYNG smoothing reached
PGGB/SmoothXG-like quality on those loci.

## Remaining Blocker

The remaining blocker is **C4-scale local seed induction inside
`gfa:syng-local:localized`**. The production path reaches the expensive
SweepGA/FastGA PAF stage, consumes roughly `46-47 GB` RSS, and times out before
emitting a seed GFA. Until that step produces a seed/output graph, the workgroup
cannot evaluate the actual hard invariant for integrated localized smoothing:
same path names, no missing/extra paths, and no spelling mismatches.

This is a runtime/scale blocker first, not a proof that the localized polishing
architecture corrupts paths or cannot improve quality.

## Next Single Highest-Value Action

Make the `syng-local:localized` seed induction path budget-aware enough to emit
an exact-path-checkable C4-scale seed/output under the validation budget, then
rerun `validate-localized-polishing`.

Concretely, the next task should target the PAF/seed-induction bottleneck:
stream or chunk the SweepGA/FastGA mapping records, cap/filter pair work before
large in-memory accumulation, and preserve enough provenance to compare every
final path spelling against the SYNG-collected local sequence set. Only after a
C4 integrated output exists should the team spend more effort interpreting
singleton, low-depth, replay/compression, or white-space movement.

## Artifact Policy Check

This synthesis commits only this markdown report. It does not add generated
GFA, compressed GFA, PNG, PDF, bulky logs, per-fixture command artifacts, or any
large binary artifact.
