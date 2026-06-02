# Peer Audit: Smoothing and Renderer Wrapper Review

Task: `peer-audit-smoothing`

Date: 2026-05-28

Reviewed artifact: `docs/external-tool-wrapper-audit.md`.

## Evaluation Summary

Overall grade for the committed `audit-external-tool` work: **0.84 / 1.00**.

Confidence: **0.74**. The task has concrete validation bullets, but no numeric
rubric, so the score is calibrated against the stated acceptance criteria and
software-audit norms rather than an explicit weighting scheme.

Dimension scores:

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Wrapper scope coverage | 0.84 | The audit covers POASTA, SPOA, smoothing, gfaffix, gfasort, gfalook, and render bundles, but under-describes several renderer and post-smooth normalization behaviors. |
| Concrete source evidence | 0.78 | Most references are concrete, but at least one cited line is stale or wrong: `docs/external-tool-wrapper-audit.md:61` points SPOA graph generation partly at `src/resolution.rs:5956`, which is currently POASTA error wrapping, not the SPOA replacement wrapper. |
| Behavior classification | 0.82 | The audit uses the requested categories well for the major crush/smoothing fixes, but misses some invoke, parse, explicit-parameter, and correctness-validation details in renderer helpers. |
| Hidden semantic change detection | 0.86 | I did not find an unreported production hidden-semantic fallback in the reviewed POASTA/SPOA/smoothing paths. The main remaining risk is under-documented post-smooth gfaffix coupling and missing render-output validation. |
| Validation discipline | 0.90 | Regression coverage cited in `docs/external-tool-wrapper-audit.md:76-95` targets the high-risk wrapper fixes; renderer-specific validation is thinner. |

Rubric underspecification flag: **partial**. The task validation defines required
document properties, but not how to weight coverage gaps, citation defects, or
low-risk renderer omissions in a numerical score.

## Scope Reviewed

I checked the committed audit against these concrete code areas:

- POASTA and SPOA replacement wrappers: `src/resolution.rs:4871`,
  `src/resolution.rs:4907`, `src/resolution.rs:4965`,
  `src/resolution.rs:5925`.
- Crush pairwise induction, smoothing/condenser, and fallback removal:
  `src/resolution.rs:4218`, `src/resolution.rs:4402`,
  `src/resolution.rs:4453`, `src/resolution.rs:4553`,
  `src/resolution.rs:4696`.
- smoothxg-style block smoothing and POVU-guided block placement:
  `src/smooth.rs:206`, `src/smooth.rs:298`, `src/smooth.rs:424`,
  `src/smooth.rs:1035`, `src/smooth.rs:1191`,
  `src/smooth.rs:2016`, `src/smooth.rs:2043`.
- gfaffix, gfasort, and graph transform dispatch:
  `src/graph.rs:732`, `src/graph.rs:767`, `src/graph.rs:844`,
  `src/lib.rs:772`.
- Render bundle and gfalook helpers: `src/commands/render.rs:45`,
  `src/commands/render.rs:133`, `src/commands/render.rs:321`,
  `src/commands/render.rs:453`, `src/main.rs:9888`,
  `src/main.rs:10067`.

## Findings

### 1. Wrong SPOA line reference in the audit table

Classification: **citation defect, not a missed wrapper behavior**.

`docs/external-tool-wrapper-audit.md:61` lists SPOA graph generation as
`src/graph.rs:310`, `src/resolution.rs:5956`, and `src/smooth.rs:2043`.
`src/graph.rs:310-320` does define the SPOA engine builders, and
`src/smooth.rs:2043-2259` is the per-block SPOA smoother. However,
`src/resolution.rs:5956` is currently:

- `src/resolution.rs:5955-5956`: `poasta_to_io_error`, returning
  `"POASTA replacement failed: {err}"`.

The crush SPOA replacement wrapper is actually:

- `src/resolution.rs:4871-4897`: `build_poa_replacement`, which builds a
  global SPOA graph, unchops it, orders replacement paths, and validates
  replacement path preservation.

Impact: low-to-moderate for the audit artifact. It does not indicate a hidden
semantic change in production code, but it weakens the "concrete files and line
numbers" evidence for the SPOA row.

### 2. Post-smooth gfaffix is coupled to `:smooth`, not independently selected

Classification: **explicit parameter / coupled stage**. If left undocumented,
it can look like a **hidden semantic change** because gfaffix modifies graph
topology after smoothing.

The audit states in `docs/external-tool-wrapper-audit.md:66` that graph engine
dispatch transforms are applied only when selected by `EngineOpts`. The broad
statement is mostly right, but it hides one important coupling:

- `src/lib.rs:785-815` runs `smooth::smooth_gfa` when
  `engine_opts.smooth_after_crush` is set.
- `src/lib.rs:822-824` then runs `graph::run_gfaffix` unconditionally inside
  that smooth branch.
- `src/lib.rs:827-837` applies the optional gfasort pipeline after that, or
  skips only sorting when `graph_sort_pipeline` is `None`.

So `gfa:syng:crush:smooth:nosort` still performs gfaffix normalization after
smooth; `:nosort` disables the gfasort stage, not gfaffix. The wrapper
implementation is consistent with the comment at `src/lib.rs:822-823` that this
matches pggb's post-smooth step, and `src/graph.rs:844-898` errors if gfaffix is
missing or exits non-zero. The audit should make this coupling explicit.

Impact: moderate audit-coverage gap. I do not see a production fallback bug
here, but a reader could incorrectly infer that gfaffix is independently
selected or disabled by `:nosort`.

### 3. gfalook wrapper checks process status but not output artifact integrity

Classification: **correctness validation gap**.

The audit row at `docs/external-tool-wrapper-audit.md:68` correctly says the
gfalook wrapper skips only graphs with no segment records and treats missing
binary/non-zero status as errors. The code path is:

- `src/main.rs:10074-10080`: skip render if no `S` records exist.
- `src/main.rs:10082-10095`: marshal GFA through either the already-written
  output path or a temporary `.gfa`.
- `src/main.rs:10098-10114`: invoke `gfalook` from `PATH` with `-i`, `-o`,
  width, height, path-height, thread count, and optional `-m`.
- `src/main.rs:10115-10130`: map missing binary/start errors and non-zero exit
  status to `io::Error`.

What is not validated: after a zero exit status, the wrapper does not check that
`render_path` exists, is non-empty, or matches the requested image type. That is
not a marshal/invoke/parse issue; it is a missing correctness validation on a
user-requested side-effect.

Impact: low-to-moderate. A faulty or incompatible `gfalook` could return success
without producing a usable image, and `impg` would report success. This does not
change graph semantics, but it is a renderer wrapper validation gap.

### 4. gfalook binary resolution is an invoke behavior not called out

Classification: **invoke**.

The gfaffix wrapper intentionally resolves only a sibling binary:

- `src/graph.rs:844-860`: `run_gfaffix` finds `gfaffix` next to the current
  executable and does not search `PATH`.

The gfalook wrapper does the opposite:

- `src/main.rs:10098`: `Command::new("gfalook")` uses normal executable lookup.

This may be a reasonable renderer policy, but it is a real invoke behavior and
differs from gfaffix's version-drift guard. The audit mentions gfaffix's sibling
lookup in `docs/external-tool-wrapper-audit.md:64`, but does not explicitly
contrast gfalook's `PATH` lookup in `docs/external-tool-wrapper-audit.md:68`.

Impact: low. It is not a hidden semantic transform, but a reproducibility risk
for renderer output.

### 5. Local render bundles have additional explicit defaults and parse validations

Classification: **explicit parameter**, **parse**, and **correctness validation**.

The render bundle row at `docs/external-tool-wrapper-audit.md:67` says local
graph rendering builds a FASTA, invokes the selected local graph engine, and
records translation metadata. That is accurate but incomplete.

Additional explicit parameters:

- `src/commands/render.rs:327-331`: local render engine `poa` uses fixed POA
  scoring `(5, 4, 6, 2, 24, 1)`.
- `src/commands/render.rs:336-345`: local render engine `pggb` uses fixed
  smoothing parameters `[700, 1100]`, max node length `100`, and padding
  fraction `0.001`.

Additional parse/correctness validation:

- `src/commands/render.rs:453-485`: translation sampling parses `S`, `P`, and
  `W` lines from the emitted local graph.
- `src/commands/render.rs:487-507`: every rendered interval must have a path
  and every path step must reference a known segment.
- `src/commands/render.rs:532-550`: path-name matching allows an exact match or
  one unique metadata-suffixed match; ambiguity returns no match and later
  errors.
- `src/commands/render.rs:600-608`: render translation currently requires GFA
  segment IDs to parse as `u32`.

These are reasonable render-bundle semantics, but they should be classified in
the audit because they are stricter than generic "marshal/invoke/parse".

Impact: low-to-moderate. The numeric segment-ID requirement can reject otherwise
valid GFA from a future engine, but it fails explicitly and therefore is a
correctness validation, not a hidden semantic change.

### 6. gfasort "empty pipeline" behavior needs a CLI/API qualifier

Classification: **explicit parameter** plus **correctness validation**.

`docs/external-tool-wrapper-audit.md:65` says empty pipeline and trivial graph
passthrough are identity cases. That is true for the lower-level helper:

- `src/graph.rs:772-775`: `sort_gfa_pipeline` returns the input unchanged if
  `pipeline.trim()` is empty.
- `src/graph.rs:781-785`: trivial graphs with zero or one nodes are returned
  unchanged.

However, the user-facing graph-engine parser rejects an empty sort pipeline:

- `src/main.rs:3159-3164`: `validate_gfasort_pipeline_param` errors when the
  parsed `pipeline` value is empty.

Impact: low. This is not a production bug, but the audit should distinguish the
library helper's identity case from the CLI's explicit-parameter validation.

## Confirmed Major Decisions

I did not find an unreported fallback substitution in the core smoothing and
replacement wrappers reviewed here.

- Pairwise polish no longer substitutes unpolished replacements:
  `src/resolution.rs:4434-4447`.
- Chain-POVU smooth-to-POASTA rejects empty traversal and smooth failures rather
  than falling back to direct POASTA: `src/resolution.rs:4457-4462` and
  `src/resolution.rs:4496-4510`.
- AllWave rejects too few non-empty traversals and zero selected pairs rather
  than falling back to SPOA: `src/resolution.rs:4557-4565` and
  `src/resolution.rs:4585-4591`.
- SweepGA replacement induction errors on zero PAF records, trusts SweepGA's
  filter only in the top-flubble tail mode, and treats no-shared replacement
  segments as a warning: `src/resolution.rs:4769-4783`,
  `src/resolution.rs:4808-4817`, and `src/resolution.rs:4821-4828`.
- wfmash and FastGA replacement inputs no longer receive an impg-side non-empty
  length floor above one base: `src/resolution.rs:4841-4851`.
- POASTA replacement and block wrappers validate exact paths:
  `src/resolution.rs:4907-4955` and `src/resolution.rs:4965-5038`.
- SPOA replacement validates paths after graph generation and unchop:
  `src/resolution.rs:4871-4897`.
- Smoothable blocks now error instead of falling back to passthrough after
  aligner failure or empty smoothing output: `src/smooth.rs:424-479`.
- Explicit POVU reference hints fail instead of falling back to graph path 0:
  `src/smooth.rs:1051-1069` and `src/smooth.rs:1213-1232`.

## Recommendation

No production code should be changed under `peer-audit-smoothing`; the assigned
file scope is documentation only. The highest-value follow-up would be a small
renderer task that adds post-success output validation to `render_graph_image`
and updates the audit table to call out `PATH` lookup, post-smooth gfaffix, and
render-bundle translation constraints.
