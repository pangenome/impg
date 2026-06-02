# Crush POA/POASTA pass-through audit

Task: `crush-poasta-pass-through-audit`

## Finding

POA/POASTA replacement acceptance is path-validation gated. For the standard
crush path, each selected candidate is built with
`build_replacement_with_method` and then kept unless the builder errors or
returns an empty graph (`src/resolution.rs:1158-1167`). The frontier rewrite is
then applied by `apply_replacement_frontier`, which rejects the whole rewrite
only if path spellings differ from the input graph
(`src/resolution.rs:6146-6149`). Per-replacement validation similarly compares
each replacement path sequence against the candidate traversal sequence
(`src/resolution.rs:7792-7817`).

No graph-quality, compression-ratio, objective, or size-improvement check gates
an already path-valid POA/POASTA replacement. Quality summaries are logged
before/after a rewrite, compression ratio is diagnostic-only, and multi-level
objective values rank/report candidates rather than veto built replacements.
The code comments for `multi_level_objective`, `IterativeMultiLevel`, and
`CoverageMultiBubble` were updated to match this behavior.

## Candidate Eligibility Budgets

The remaining non-acceptance controls are candidate eligibility or construction
budgets:

- CLI parsing accepts `method=poa`, `method=poasta`, `method=auto`, and
  `polish-method=poasta` through the crush parser (`src/main.rs:2647-2760`,
  `src/main.rs:2957-2969`, `src/main.rs:8674-8681`).
- Discovery requires a POVU-supported bubble with at least two uniquely
  anchored traversals (`src/resolution.rs:4911-4941`).
- `min_traversal_len` can filter legacy one-bubble candidates before
  materialization, but tree-driven and multi-bubble modes explicitly keep the
  full POVU topology for later path-validation checks
  (`src/resolution.rs:4963-4978`).
- Candidate materialization skips only traversal-identical sites because there
  is no replacement to perform (`src/resolution.rs:5884-5907`).
- Non-overlap scheduling chooses which path-valid candidates can be applied in
  the same round; deferred overlaps are scheduling, not quality rejection.
- Pairwise graph-induction methods still have construction budgets such as
  selected pair limits, PAF byte limits, PAF filters, and backend failures. When
  those produce an error or empty graph, no replacement exists to accept.
- Multi-level window caps (`candidate-limit`, max window sites, target bp) bound
  which candidate windows are built. Objective floors are diagnostic only.

## Regression

`direct_poa_and_poasta_accept_topology_changing_replacements` builds a tiny
shared-prefix bubble:

- input segments: `A`, `CTT`, `CGG`, `T`
- accepted POA/POASTA replacement segmentation includes shared `AC`
- path spellings remain `ref=ACTTT` and `alt=ACGGT`

The test sets an intentionally unreachable compression threshold and asserts
that both direct SPOA (`method=poa`) and POASTA still accept the replacement,
preserve path spellings, report one resolved replacement, and perform no retry
or metric-based swap.

## Logging

Replacement build progress now reports `accepted`, `empty`, `failed`, or
`path-invalid`. Round summaries separately count `failed`, `empty`, and
`path-invalid` outcomes (`src/resolution.rs:1171-1178`,
`src/resolution.rs:1186-1235`, `src/resolution.rs:7899-7918`).
