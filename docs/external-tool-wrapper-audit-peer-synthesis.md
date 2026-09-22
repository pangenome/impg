# External Tool Wrapper Audit Peer Synthesis

Task: `synthesize-external-wrapper`

Date: 2026-05-28

## Inputs Reviewed

- Primary audit artifact: `docs/external-tool-wrapper-audit.md`.
- Supporting C4 follow-up artifact: `docs/top-flubble-sweepga.md`.
- Code and regression evidence cited by the primary audit:
  `src/commands/align.rs`, `src/graph.rs`, `src/resolution.rs`,
  `src/smooth.rs`, `src/syng_graph.rs`, and
  `tests/test_crush_integration.rs`.

At synthesis time, no completed peer review artifact was available. The graph
had started `peer-audit-aligner` and `peer-audit-smoothing`, but their expected
outputs, `docs/external-tool-wrapper-audit-aligner-review.md` and
`docs/external-tool-wrapper-audit-render-review.md`, did not yet exist in their
worktrees and the tasks were not done. This synthesis therefore cites no peer
review findings as completed evidence.

## Synthesis

The primary audit establishes the policy that wrappers may marshal inputs and
outputs, invoke tools with documented parameters, parse outputs, validate
correctness invariants, and log diagnostics, but must not silently add filters,
fallbacks, substitutions, name mutations, or minimum-length/identity thresholds
outside explicit tool parameters (`docs/external-tool-wrapper-audit.md:7`,
`docs/external-tool-wrapper-audit.md:19`).

The audit's decisions are coherent with that invariant. It reports removal of
hidden behavior in the alignment, local graph induction, replacement, polishing,
smoothing, POVU reference-hint, and per-pair alignment wrappers
(`docs/external-tool-wrapper-audit.md:27`). The wrapper table then classifies the
remaining integration points as retained, fixed, explicit-parameter behavior,
marshal/invoke/parse behavior, or correctness validation
(`docs/external-tool-wrapper-audit.md:40`). The C4 fragment regressions and full
validation list cover the key failure modes that motivated the audit, including
PAF rescue after explicit filtering, POVU reference-hint fallback, hidden wfmash
input floors, and C4 fragment seqwish induction under local graph defaults
(`docs/external-tool-wrapper-audit.md:76`).

The top-level flubble artifact is consistent with the audit result. It states
that `top-flubble-sweepga` accepts only validity checks, not quality guards:
replacement GFA must parse and preserve path names and spellings
(`docs/top-flubble-sweepga.md:5`). It also records that top-level mode trusts
SweepGA's own replacement-tier PAF filter and disables the second seqwish-tail
filter (`docs/top-flubble-sweepga.md:29`), and that the C4 run preserved all
path names and spellings (`docs/top-flubble-sweepga.md:99`). The remaining C4
quality problem, many short top-level regions passing through as valid but
stringy replacements, is not itself a hidden wrapper semantic because the task's
policy intentionally made quality diagnostic-only and removed fallback guards
(`docs/top-flubble-sweepga.md:146`).

## Follow-Up Decision

One unresolved semantic risk remains and was converted into a concrete
implementation task: `fix-syng-query`.

The primary audit explicitly leaves the shared `ImpgIndex` trait as a remaining
guardrail because it returns vectors rather than `Result`, preventing every
syng query-layer failure from propagating (`docs/external-tool-wrapper-audit.md:107`).
The code still has two places where backend/filter failures are represented as
empty semantic output:

- `SyngImpgWrapper::query_via_syng_raw` logs `query_region` failures and returns
  an empty vector (`src/lib.rs:304`, `src/lib.rs:310`, `src/lib.rs:317`).
- `SyngImpgWrapper::query_via_syng_filtered` logs
  `query_region_with_anchors_ext` failures and returns an empty vector
  (`src/lib.rs:338`, `src/lib.rs:346`, `src/lib.rs:354`).
- The `ImpgIndex` implementation then delegates `query`, `query_with_cache`,
  `query_transitive_dfs`, and `query_transitive_bfs` to those vector-returning
  helpers (`src/lib.rs:417`, `src/lib.rs:434`, `src/lib.rs:458`,
  `src/lib.rs:477`).
- Syng anchor chaining also logs SweepGA scaffold-filter failures and returns an
  empty vector (`src/syng_transitive.rs:209`, `src/syng_transitive.rs:215`).

Those cases make a wrapper/backend failure indistinguishable from a successful
query with no hits. That is the same class of hidden semantic behavior the audit
rules disallow: the caller receives a valid-looking empty result rather than an
error that explains the query could not be completed. The follow-up task
`fix-syng-query` was created with file scope for `src/impg_index.rs`,
`src/lib.rs`, `src/syng_transitive.rs`, affected `ImpgIndex` implementers and
call sites, plus focused tests. Its validation requires failing tests first,
error propagation for syng backend and scaffold-filter failures, preservation of
successful no-hit empty results, real CLI-path coverage where user-visible
behavior is affected, and full build/test/install validation.

No other follow-up implementation task is needed from the available evidence.
The audit's fixed wrapper behaviors have direct regression coverage
(`docs/external-tool-wrapper-audit.md:78`), the C4 fragment validation preserved
fixture paths while retaining shared graph induction evidence
(`docs/external-tool-wrapper-audit.md:97`), and no completed peer artifact added
another missed behavior.

## Validation Notes

- Synthesis cites the primary `audit-external-tool` artifact and supporting C4
  artifact.
- No completed peer review artifact was available to cite.
- The unresolved syng query error-as-empty behavior was converted into
  `fix-syng-query` with concrete file scope and validation.
- This task modified only `docs/external-tool-wrapper-audit-peer-synthesis.md`
  and created the workgraph follow-up task; no production code was modified.
