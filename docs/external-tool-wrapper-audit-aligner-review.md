# External Tool Wrapper Audit: Aligner Review

Task: `peer-audit-aligner`

Date: 2026-05-28

Reviewed artifact: `docs/external-tool-wrapper-audit.md`

## Summary

The committed audit is directionally strong: it identifies the prior hidden PAF
rescue, removed the impg-side wfmash length floor, tightened smoothing fallback
behavior, and added C4 fragment coverage for the seqwish tail. The central
invariant in `docs/external-tool-wrapper-audit.md:7-25` is the right standard for
these wrappers.

I would score the audit at **0.78 / 1.00** with **0.72 confidence**. It covers the
direct wrapper call sites well, but it under-classifies several behaviors where
the wrapper result is still changed or suppressed by the surrounding graph
induction path.

Dimension scores:

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Source coverage | 0.84 | SweepGA/FastGA, wfmash, seqwish, SPOA/POASTA, AllWave, smooth, and syng-native WFA paths are present in the audit table. |
| Concrete evidence | 0.86 | The audit cites files and line numbers, and the cited fixes are mostly present in current source. |
| Missed-semantics classification | 0.62 | Several residual behaviors are not classified as explicit parameters, correctness validation, or hidden semantic changes. |
| Regression coverage | 0.75 | Added tests cover the PAF filter rescue and C4 seqwish-tail regressions, but not wrapper-error propagation, zero-PAF AllWave, ignored sparse-pairs, or syng anchor-filter dropout. |
| Production-code scope | 1.00 | This peer review modifies documentation only. |

Rubric underspecification flag: **not blocking**. The task gives concrete
acceptance criteria but no numeric grading rubric; the score above is calibrated
against the stated wrapper invariant and validation checklist.

## Confirmed Good Coverage

- `src/commands/graph.rs:150-157` delegates graph construction through
  `align_sequences` and `induce_graph_from_alignment`, and
  `src/commands/graph.rs:731-970` combines FASTA, invokes/copies PAF, applies the
  shared filter, and returns the temp files. This is correctly classified as
  marshal/invoke/parse plus explicit parameters.
- `src/commands/graph.rs:675-727` applies the shared SweepGA PAF filter and
  returns early on `no_filter`, matching the audit's explicit-parameter
  classification.
- `src/commands/align.rs:1028-1095` now errors out of generic per-pair
  alignment when one pair invocation fails, rather than returning a partial PAF.
  This matches the audit decision at `docs/external-tool-wrapper-audit.md:38`.
- `src/syng_graph.rs:796-923` writes exact FASTA names, runs only
  `filter_generated_paf`, optionally lowers `min_match_len` only under
  `adaptive_min_match_len`, and then calls the shared seqwish tail. The regression
  at `src/syng_graph.rs:1409-1440` directly protects against re-adding raw PAF
  records after explicit filter removal.
- `src/resolution.rs:4218-4276` makes replacement seqwish filtering explicit:
  `sweepga_no_filter` sets `min_match_len=1`, disables identity filtering, and
  forces many-to-many/scaffold-mass-zero semantics. The tests at
  `src/resolution.rs:6988-7076` cover those default and no-filter cases.
- `src/resolution.rs:4402-4450`, `src/resolution.rs:4453-4511`, and
  `src/smooth.rs:424-487` now return errors rather than substituting unpolished,
  direct-POASTA, or passthrough outputs from failed smoothable blocks.
- `src/smooth.rs:1051-1070` and `src/smooth.rs:2606-2624` now fail an explicit
  POVU reference-hint miss rather than falling back to the first graph path.
- `tests/test_crush_integration.rs:427-462` and
  `tests/test_crush_integration.rs:465-580` are useful C4 fragment guards for
  path spelling, shared segment evidence, and unfolded/duplicated replacement
  signatures.

## Findings

### 1. Replacement-wrapper errors still become original-graph retention in the outer resolver

Classification: **hidden semantic change** and **correctness validation gap**.

The audit states that polish/smooth/AllWave/SweepGA wrapper failures now error
instead of substituting fallback outputs. That is true inside the immediate
wrapper functions, but the normal resolver catches those `Err` values and keeps
going:

- `src/resolution.rs:1002-1011` builds a replacement with `?`, but
  `src/resolution.rs:1034-1053` records `Err` as `stats.bailed` and does not
  propagate it.
- If every selected replacement fails, `src/resolution.rs:1056-1068` continues
  to the next round with no applied plan, leaving the current graph unchanged.
- The chain-greedy path has the same shape at `src/resolution.rs:1385-1395` and
  `src/resolution.rs:1418-1435`.
- Top-flubble SweepGA has an explicit zero-alignment identity fallback:
  `src/resolution.rs:1923-1934` skips a zero-PAF region and leaves the original
  interval intact, while `src/resolution.rs:1942-1952` returns the current graph
  if every selected top-level region emitted zero PAF records.

This means a wrapper failure can still silently become "no replacement applied"
at the workflow level. That is not the same fallback shape as the removed
unpolished/direct-SPOA substitutions, but it is still a fallback output under the
audit invariant's wording in `docs/external-tool-wrapper-audit.md:19-25`.

Recommended follow-up: split "candidate not eligible" from "wrapper invocation
or correctness validation failed." For explicit `method=sweepga`, `method=wfmash`,
`method=allwave`, `method=poasta`, or `method=chain-povu`, wrapper failure should
either propagate or require an explicit user parameter such as
`allow-failed-replacement-passthrough`. Add tests that force a wrapper failure and
assert the command fails rather than reporting only `stats.bailed`.

### 2. `sweepga_sparse_pairs` is documented but ignored in replacement induction

Classification: **explicit parameter violation** with a **hidden semantic
change**.

`ResolutionConfig` documents `sweepga_sparse_pairs` as a control that can force
sparse pair dispatch for SweepGA/FastGA graph induction:

- `src/resolution.rs:167-175` describes `sweepga_sparse_pairs` as user-facing
  replacement-induction behavior.
- The CLI parser accepts `sparse-pairs=true`; see the parse expectation at
  `src/main.rs:13238-13262`.

The replacement builder then ignores it:

- `src/resolution.rs:4709-4712` logs that sparse-pairs is ignored.
- `src/resolution.rs:4727-4748` constructs the SweepGA align config with
  `SparsificationStrategy::None`.

The audit's SweepGA/seqwish row at `docs/external-tool-wrapper-audit.md:56`
classifies the path as explicit user parameters, but this particular parameter
does not affect invocation semantics. That is a missed behavior because it can
change pair coverage, runtime, and resulting PAF density while making the user's
explicit request a no-op.

Recommended follow-up: either remove/deprecate the parameter from replacement
induction with a clear parse-time error, or wire it through to the same
pair-selection semantics used by `src/commands/align.rs:801-865` and
`src/commands/align.rs:937-1026`.

### 3. AllWave guards zero selected pairs, but not zero emitted PAF

Classification: **hidden semantic change** and **correctness validation gap**.

The audit correctly notes that AllWave now rejects fewer than two non-empty
traversals and zero scheduled pairs:

- `src/resolution.rs:4557-4565` rejects candidates with fewer than two non-empty
  traversals.
- `src/resolution.rs:4585-4590` rejects an empty pair schedule.

However, after AllWave runs, empty PAF output is still accepted:

- `src/resolution.rs:4593-4608` runs each selected pair and drops empty PAF
  strings.
- `src/resolution.rs:4611-4617` converts zero unique PAF lines into an empty PAF
  string.
- `src/resolution.rs:4626-4628` passes that empty PAF into the seqwish tail and
  final replacement validation.

Path preservation can still pass for an unfolded graph with no shared induced
structure, so this is not fully covered by the current correctness validation.
The SweepGA path has a zero-PAF guard at `src/resolution.rs:4769-4783`; the
AllWave path should have an equivalent guard or should explicitly document that
empty-PAF AllWave replacement is allowed.

Recommended follow-up: fail when `pairs` is non-empty but `lines` is empty, and
add a regression that asserts AllWave does not convert zero alignment evidence
into a path-preserving but unshared replacement graph.

### 4. Syng-native WFA silently drops selected pair evidence

Classification: **correctness validation that currently acts as hidden
filtering**.

The audit classifies syng-native BiWFA output as correctness validation. The
source does validate CIGAR consumption, but the validation result is represented
as `None`, and callers concatenate only successful PAF lines:

- `src/syng_graph.rs:335-397` returns `None` on empty input, incomplete WFA
  status, zero block length, or CIGAR consumption mismatch.
- `src/syng_graph.rs:410-419` uses `filter_map` for all-pairs PAF and drops any
  pair returning `None`.
- `src/syng_graph.rs:430-452` does the same for sparse-pair PAF.
- `src/syng_graph.rs:785-790` then feeds the remaining PAF, possibly empty or
  partial, into seqwish.

Rejecting invalid CIGARs is valid correctness validation, but silently dropping a
selected pair is also a graph-induction semantic change. The audit should have
called out whether WFA failure/invalid CIGAR is expected to fail the wrapper or
is an explicit "best-effort PAF" policy.

Recommended follow-up: return a `Result` with pair-level error context from the
WFA PAF builders, or at least count selected pairs versus emitted PAF lines and
error when the gap exceeds an explicit configured tolerance.

### 5. Anchor-seeded syng WFA uses hard-coded SweepGA scaffold filtering and can skip rather than fallback

Classification: **hidden semantic change** unless promoted to an explicit
algorithm parameter.

The audit says anchor-seeded syng BiWFA uses a SweepGA anchor scaffold filter and
that fallback to full-pair BiWFA affects performance/coverage. The implementation
has two additional semantics that should be recorded:

- `src/syng_graph.rs:696-718` hard-codes the anchor filter to one-to-one
  scaffold mode, `min_scaffold_length = syncmer_len * 5`,
  `scaffold_gap = 50_000`, `overlap_threshold = 0.95`, and
  `scaffold_max_deviation = syncmer_len * 10`.
- `src/syng_graph.rs:720-729` returns an empty anchor set if the SweepGA anchor
  filter itself errors.
- In the anchor driver, `src/syng_graph.rs:1144-1175` returns `None` when raw
  pair anchors exist but every anchor is removed by the plane-sweep filter. That
  path does not fall back to full-pair BiWFA; fallback happens only later at
  `src/syng_graph.rs:1179-1181` when there were no usable anchors for the pair.

That behavior can reduce PAF coverage based on fixed local thresholds rather
than only marshalling/querying/validating WFA output. If these thresholds are
part of the syng-native engine definition, the audit should classify them as
explicit documented defaults. As written, they look like hidden wrapper-side
filtering.

Recommended follow-up: expose these thresholds in the syng-native graph config,
or make the filter failure/all-filtered case fall back to full-pair BiWFA unless
the anchor set is provably invalid for correctness reasons.

### 6. `impg align` wfmash percent identity is hard-coded at invocation time

Classification: **explicit parameter gap** and possible **hidden semantic
change**.

The graph/query paths carry `map_pct_identity` through `GraphBuildConfig`:
`src/commands/graph.rs:100-101`, `src/commands/graph.rs:136`, and
`src/main.rs:9394-9399`. The direct `impg align` path does not have the same
field in `AlignConfig`:

- `src/commands/align.rs:23-49` defines `AlignConfig` without
  `map_pct_identity`.
- `src/main.rs:8436-8460` builds that config from CLI arguments without carrying
  `aln.sw.map_pct_identity`.
- `src/commands/align.rs:533-546` and `src/commands/align.rs:566-582` invoke
  wfmash/SweepGA with `Some("90".to_string())`.

The audit's aligner-factory row at `docs/external-tool-wrapper-audit.md:50`
states that documented backend options pass through to SweepGA integrations.
That is true for graph/query construction but not for this direct align command
path.

Recommended follow-up: add `map_pct_identity` to `AlignConfig`, populate it from
`aln.sw.map_pct_identity`, and use it in the direct and batched align invocations.

## Coverage Gaps To Add

- A wrapper-failure integration test that proves explicit replacement methods
  fail instead of merely incrementing `stats.bailed` and preserving the original
  graph.
- An AllWave test with selected pairs but no emitted PAF records, expecting an
  error rather than an unfolded replacement.
- A SweepGA replacement test proving `sweepga_sparse_pairs=true` changes
  invocation semantics or is rejected as unsupported.
- A syng-native WFA test that selected pair count and emitted PAF count cannot
  silently diverge without an explicit best-effort setting.
- A direct `impg align` config test proving `map-pct-identity` reaches the
  wfmash integration.

## Bottom Line

The original audit removed several concrete hidden transformations and improved
the C4 seqwish-tail regression suite. The remaining risk is one layer higher:
some wrapper failures, empty alignment evidence, and hard-coded pair/anchor
filters still become path-preserving but semantically different graph outputs.
Those should be classified explicitly before treating the wrapper invariant as
fully enforced.
