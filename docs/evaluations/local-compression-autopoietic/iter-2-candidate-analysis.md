# Local Compression Candidate Analysis After Iteration 1

Date: 2026-06-07

Task: `analyze-next-local`

## Scope

This report compares the remaining local-compression algorithm candidates after
iteration 1 of `autopoietic-local-compression-loop`. It is a decision report
only: no runner defaults, fixture metadata, or testbed code were changed.

Source anchors:

- Candidate method definitions and fixed parameters are in
  `scripts/local_compression_testbed.py:122`, `:129`, `:136`, `:143`, and
  internal controls at `:150`, `:157`, `:164`.
- Current chunk-window behavior discovers same-length variant runs, merges runs
  only when the invariant gap is below `CHUNK_WINDOW_MIN_INVARIANT_GAP_BP`, and
  falls back to `compact_bubble` when only one window remains
  (`scripts/local_compression_testbed.py:927`, `:950`, `:966`, `:972`).
- Iteration 1 established that `chunk_window_smooth_or_crush` is the only method
  passing `nested_top_level_wrong`, while internal controls preserve exact paths
  but fail current topology assertions
  (`docs/evaluations/local-compression-autopoietic/iter-1.md:68` and `:99`).
- The fast-profile control mapping and coverage are recorded in
  `docs/evaluations/local-compression-testbed-fast/report.md:35` and `:45`.

## Baseline Signal

The post-iteration-1 fast profile has 13 fixtures, 11 methods, and 143 rows.
There are 11 graph-producing CI rows per method in the comparison below; the two
local-tier fixtures remain skipped by fast-profile policy.

| Method or control | Candidate hypothesis | Exact paths on graph rows | Topology pass/fail on graph rows | Mean graph bytes | Mean segments | Mean path steps | Mean white-space proxy | Mean runtime seconds |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `chunk_window_smooth_or_crush` | Sorted same-length chunk windows can split independent local events. | 11/11 | 11/0 | 670.1 | 5.364 | 12.364 | 2.545 | 0.000364 |
| `top_flubble_nonoverlap_sweepga` | Highest-level non-overlapping flubble windows avoid nested overmerge. | 11/11 | 10/1 | 588.5 | 4.909 | 10.364 | 0.545 | 0.000347 |
| `local_syng_crush_sweepga` | SweepGA/seqwish local resolver at current `k=19` improves local topology. | 11/11 | 10/1 | 588.5 | 4.909 | 10.364 | 0.545 | 0.000391 |
| `whole_region_sweepga_seqwish` | Whole-region seqwish induction can resolve the local window without special boundaries. | 11/11 | 10/1 | 588.5 | 4.909 | 10.364 | 0.545 | 0.000351 |
| `pggb_control` | Internal PGGB-style compactness/path-spelling control. | 11/11 | 0/11 | 308.5 | 3.727 | 3.727 | 0.000 | 1.433 |
| `smoothxg_control` | Internal SmoothXG-style compactness/path-spelling control. | 11/11 | 0/11 | 308.5 | 3.727 | 3.727 | 0.000 | 0.03354 |
| `pggb_plus_smoothxg_control` | Internal PGGB plus SmoothXG-style compactness/path-spelling control. | 11/11 | 0/11 | 308.5 | 3.727 | 3.727 | 0.000 | 1.423 |

Interpretation: exact path spelling no longer separates these candidates. The
controls are the compactness reference, not the current topology winner: each
control is much smaller than the local candidates but fails all 11 synthetic
topology assertions. The local candidates pass most topology checks, but only
sorted chunk windows pass the promoted nested-parent adversarial case.

## Candidate Hypotheses

### 1. Sorted Chunk Windows

Hypothesis: preserve exact paths while splitting separable same-length local
variant runs into independent bounded windows.

Evidence for:

- It is the only candidate with 11/11 topology passes on graph-producing CI
  rows.
- It uniquely passes `nested_top_level_wrong` with `candidate_count=2`,
  `bubble_count=2`, and `flubble_count=2`; the scoreboard row is line 140 in
  `docs/evaluations/local-compression-testbed-fast/scoreboard.tsv`.
- The promoted fixture is exactly the failure mode targeted by iteration 1:
  two independent events inside an overbroad parent. The variant runs are split
  by a 5 bp invariant gap, and the current strict `< 5` merge rule keeps them
  separate.

Evidence against:

- It is less compact than every internal control: 670.1 mean bytes and 12.364
  mean path steps versus 308.5 mean bytes and 3.727 mean path steps for each
  control.
- It is also larger than the one-parent local placeholder family. On
  `nested_top_level_wrong`, the chunk row emits 931 bytes and 20 path steps
  versus 668 bytes and 12 path steps for the flubble/SweepGA placeholder rows.
- Its success is currently a topology-shape success, not evidence of
  PGGB/SmoothXG-like compression quality.

Assessment: strongest next boundary policy, but not yet a compression-quality
winner. It should be used to fix candidate-window selection before resolver
parameter sweeps.

### 2. Highest-Level Non-Overlapping Flubble Windows

Hypothesis: selecting only highest-level non-overlapping flubbles should avoid
overlapping-window artifacts while retaining compact local topology.

Evidence for:

- It preserves exact paths on all 11 graph-producing rows.
- It ties the compact placeholder family on mean size, path steps, white-space
  proxy, and runtime, and it passes 10/11 topology assertions.
- Its mean graph is closer than chunk windows to the internal controls, though
  it is still larger than them.

Evidence against:

- It fails the only fixture added specifically to expose top-level overmerge:
  `nested_top_level_wrong`.
- On that fixture, it still reports one candidate, one bubble, and one flubble
  at line 139 of `scoreboard.tsv`; this is the same one-parent boundary failure
  as `local_syng_crush_sweepga` and `whole_region_sweepga_seqwish`.
- Because it fails before the resolver can matter, it is not a good standalone
  next implementation target.

Assessment: keep as a comparison arm, but do not prioritize until its boundary
selection is split by a stricter local-window detector.

### 3. SweepGA/Seqwish `k` and Min-Match Sweeps

Hypothesis: changing seqwish/SweepGA sensitivity, such as lowering `k` or adding
min-match variants, could recover topology without a new boundary selector.

Evidence for:

- Current SweepGA/seqwish rows preserve exact paths on graph-producing fast
  rows.
- Both `local_syng_crush_sweepga` and `whole_region_sweepga_seqwish` are
  SYNG-fast in the synthetic profile, so a future sweep would be cheap to run
  once the tested parameter space is meaningful.

Evidence against:

- The current fixed `seqwish_k=19` rows fail exactly where the boundary selector
  emits one overbroad candidate. On `nested_top_level_wrong`,
  `local_syng_crush_sweepga`, `top_flubble_nonoverlap_sweepga`, and
  `whole_region_sweepga_seqwish` all have `candidate_count=1`,
  `bubble_count=1`, and `flubble_count=1` at lines 138, 139, and 141 of
  `scoreboard.tsv`.
- The whole-region row also fails the promoted nested fixture, so simply giving
  seqwish more context did not solve this class of topology assertion.
- A `k` or min-match sweep can tune alignment sensitivity, but it is unlikely to
  produce the missing two-window decomposition when the candidate generator
  presents only one window.

Assessment: postpone broad `k` and min-match sweeps. They should happen after
the next candidate generator emits the right window count on the adversarial
fixture.

## Control Comparison

The internal PGGB/SmoothXG controls are valuable because they define the
compactness pressure missing from the local topology pass/fail assertions:

- All three controls preserve exact path names and spellings on every
  graph-producing CI row.
- All three controls fail current synthetic topology assertions on every
  graph-producing CI row.
- All three controls are much smaller than local candidates: 308.5 mean graph
  bytes and 3.727 mean path steps, with zero white-space proxy, compared with
  588.5 to 670.1 mean graph bytes and 10.364 to 12.364 mean path steps for the
  local compact/window candidates.

Therefore the next local candidate should not be chosen solely by topology pass
count. It must preserve the boundary win from sorted chunk windows while moving
toward the controls on graph size, path replay depth, and white-space.

## Decision

Recommended next change: implement one chunked SweepGA/seqwish candidate that
uses the existing sorted chunk-window boundary detector with
`min_invariant_gap_bp=5`, keeps the resolver at the current `seqwish_k=19`, and
does not introduce a min-match sweep yet.

Rationale:

- This is one implementation change: replace the one-parent candidate window
  fed to the SweepGA/seqwish local candidate with the already validated sorted
  chunk-window selector.
- It preserves the only demonstrated fix for `nested_top_level_wrong`: two
  windows across the 5 bp invariant gap.
- It avoids spending the next iteration on a broad resolver parameter sweep
  before the window selector is capable of emitting the expected candidate
  count.
- It creates a direct test of whether the SweepGA/seqwish resolver can improve
  compression-like metrics once boundary selection is no longer the dominant
  failure.

Deferred options for the next iteration:

- Do not tune `min_invariant_gap_bp` upward in the next pass; the adversarial
  fixture's separating invariant gap is exactly 5 bp, and the current strict
  merge rule is what keeps the two windows separated.
- Do not prioritize highest-level flubble windows as a standalone algorithm
  candidate until they inherit or match the sorted chunk-window split.
- Do not run a broad `k` or min-match sweep until the candidate generator passes
  the nested-parent window-count check.

Expected validation for the follow-up implementation:

- `nested_top_level_wrong` remains exact-path `pass` and topology `pass` with
  `candidate_count=2`, `bubble_count=2`, and `flubble_count=2`.
- Mean graph/path metrics are compared against the three internal controls, not
  only against the local placeholder family.
- Runner defaults and existing control rows remain visible; no hidden quality
  gate is added.
