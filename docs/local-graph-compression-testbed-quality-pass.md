# Local Graph Compression Testbed Quality Pass

Date: 2026-06-07
Task: `quality-pass-local-graph-compression-testbed`

## Scope

This quality pass tightened the downstream WG task batch for the local graph
compression testbed before implementation starts. The target batch is:

1. `design-local-graph-compression-testbed`
2. `implement-local-compression-fixtures`
3. `implement-local-compression-testbed-runner`
4. `evaluate-local-compression-testbed`

The dependency order is now explicit as:

```text
design -> fixtures -> runner -> scoreboard/report
```

## Downstream Contract

The edited downstream descriptions require workers to produce a tiny,
adversarial, reproducible fixture suite with exact expected path spellings and
fixture-specific topology/compression expectations. Fixture metadata must cover
the fixture id, class, paths, expected path spellings, expected topology or
metric ranges, known failure mode, CI/local tier, and optional render hints.

The required fixture categories are:

- SNP bubble
- Short indel
- 50-500 bp insertion
- Alu-like insertion
- Adjacent bubbles that should compress together
- Bubble split by a fake repeat anchor
- Repeated motif or microtangle
- Duplicated flank requiring path context
- Tandem copy-number loop that should remain cyclic
- Dispersed repeat glue that should be broken or ignored
- Inversion-like case
- Nested bubbles where top-level compression is right
- Nested bubbles where top-level compression is wrong

The runner task now requires a method matrix covering raw/local SYNG, SYNG plus
flubble/local crush variants, highest-level non-overlapping flubble windows,
SmoothXG-style chunk windows, whole-region SweepGA/seqwish, and optional
PGGB/SmoothXG controls when installed and bounded.

The scoreboard contract requires exact path preservation, expected-topology
pass/fail, graph size, node/path-depth distributions, white-space proxy,
self-loop and repeat-loop counts, flubble/bubble counts, long links, runtime,
command status, skip reason, output paths, and command log paths.

## Guardrail

The downstream specs now state that no hidden filtering, candidate suppression,
quality guard, or silent acceptance gate is allowed. The only permitted hard
rejection is exact path spelling corruption. All other quality problems must
remain visible in TSV/JSON scoreboard rows and in the final report.

## Quality Assessment

Underspecification flag: false after this pass.

Dimension scores:

- Dependency clarity: 1.00
- Fixture specificity: 0.95
- Method matrix coverage: 0.95
- Metric and artifact coverage: 1.00
- Correctness acceptance criteria: 0.95
- Usability acceptance criteria: 0.95

Overall readiness grade: 0.97

Confidence: 0.86. The remaining uncertainty is implementation-path specific:
the design worker may choose exact paths or helper APIs that differ from the
examples named in the descriptions. The task contracts now require those
choices to be written into the design and consumed by later tasks instead of
remaining chat-only context.

## Verification Notes

WG descriptions were updated in place for the design, fixture, runner, and
scoreboard/report tasks. The downstream validation sections now force both
correctness and usability by requiring exact path checks, topology assertions,
complete TSV/JSON scoreboards, per-fixture notes, optional failure renders, and
clear reproduction details.
