# Audit and Validate: C4 SPOA Replacement Lacing Dedup

Task: `audit-and-validate`
Date: 2026-06-04

## Verdict

The audited failure mode is confirmed: the replacement graphs were built and
carried into `apply_replacement_frontier`, but the old render/lace-in path
treated each `OutNode::Replacement(plan_idx, node_idx)` as identity-distinct and
therefore minted duplicate segment IDs for byte-identical replacement segment
sequences.

The current chat-agent patch addresses that mechanism. `render_rewritten_graph`
now projects replacement segment sequence bytes through a sequence-to-ID map,
reuses an already emitted ID when the sequence is already present, and only mints
a fresh ID for the first occurrence of a sequence. Targeted tests and a saved C4
slice rerun pass the hard path-sequence preservation gate.

Calibrated patch-readiness grade: `0.86` with `0.82` confidence. The core
renderer bug is fixed and validated on the right path, but I am not assigning a
higher score because the full C4 true-SPOA visual rerun did not complete during
this audit and the output can still contain duplicate segment sequences inherited
from surviving original nodes. That residual duplication is not a failure of the
replacement seq-to-ID projection and was not treated as a hard gate.

Rubric underspecification flag: `false` for the engineering decision, because
the task gives one hard gate (path-sequence preservation) and one specific
mechanism to audit (seq-to-ID replacement projection). There is no detailed
numeric rubric, so the grade above is calibrated against the stated acceptance
criteria rather than against a formal point schedule.

## Code Evidence

Replacement graphs are built and laced in as replacement nodes:

- `src/resolution.rs:7509-7522` maps each accepted replacement path step to
  `OutNode::Replacement(plan_idx, step.node)`.
- `src/resolution.rs:7548-7581` rewrites each affected path by splicing those
  replacement steps into the path sequence.
- `src/resolution.rs:7584-7609` passes the replacement graphs and rewritten
  paths into `render_rewritten_graph`.
- `src/resolution.rs:7630-7634` keeps exact path-sequence validation as the hard
  post-render gate.

The seq-to-ID projection is present in the renderer:

- `src/resolution.rs:9968-9969` creates `id_by_node` and the
  `id_by_replacement_seq` sequence projection map.
- `src/resolution.rs:10069-10079` emits surviving original nodes and seeds the
  sequence projection map with their sequences.
- `src/resolution.rs:10081-10093` handles replacement nodes: if the replacement
  segment sequence already has an ID, the renderer reuses it and skips emitting a
  duplicate `S` line; otherwise it mints the next fresh ID.
- `src/resolution.rs:10098-10124` reconstructs links from the projected IDs, so
  paths point at the canonical emitted segment IDs.

The patch also covers the insertion-order corner case needed for path-start
replacement lacing:

- `src/resolution.rs:10008-10064` tracks pending replacements and inserts them
  before the first surviving original node when a replacement starts a path.
- `src/resolution.rs:10294-10318` tests that path-start replacement ordering
  preserves path spelling.
- `src/resolution.rs:10320-10355` tests that independent replacements with the
  same sequence share one emitted segment while preserving both path spellings.

## Validation

Initial validation found one compile blocker in the uncommitted patch:
`cargo check` failed with `E0425` because the renderer referenced
`original_ordered` before declaration. The current working tree fixes that by
using `original_node_count` before `ordered_nodes` is moved.

Passed validation commands:

```bash
cargo check
cargo test --lib replacement_at_path_start_is_ordered_before_first_surviving_original
cargo test --lib independent_replacements_with_same_sequence_share_one_emitted_segment
cargo test --test test_crush_integration c4_slice_auto_crush_preserves_path_sequences -- --nocapture
cargo run --bin impg -- crush --gfa tests/test_data/crush/c4_slice_1500_3000.gfa --method auto --max-iterations 1 --output /tmp/audit-and-validate/c4_slice_auto1_seqdedup.gfa
cargo run --example compare_gfa_paths -- tests/test_data/crush/c4_slice_1500_3000.gfa /tmp/audit-and-validate/c4_slice_auto1_seqdedup.gfa
cargo build
cargo test --lib
```

Results:

- `cargo check`: passed, with pre-existing warnings.
- Targeted path-start ordering test: passed.
- Targeted duplicate replacement projection test: passed.
- C4 slice integration test: passed; `465` paths preserved, `147` resolved,
  `0` bailed, segments `2942 -> 2366`, segment bp `64348 -> 57281`.
- Saved one-round C4 slice CLI rerun: completed; `147` resolved, `0` bailed,
  IDs minted `272213933..272214520`.
- Path comparison for the saved rerun: `465` expected paths, `465` observed
  paths, `0` missing, `0` extra, `0` spelling mismatches.
- Duplicate-sequence count, audit-only: input slice had `2942` segments,
  `2241` unique sequences, `701` duplicate extras across `242` duplicate
  groups; the saved one-round output had `2366` segments, `2036` unique
  sequences, `330` duplicate extras across `76` duplicate groups.
- `cargo build`: passed, with pre-existing warnings.
- `cargo test --lib`: passed, `377` passed, `0` failed.

## C4 Interpretation

The small C4 rerun supports the seq-to-ID hypothesis. Earlier audit material
reported the pre-projection slice-auto one-round output as `2984` segments and
`58003` segment bp with substantial duplicate replacement emission. The current
patch emits `2366` segments and `57281` segment bp for the same slice while
preserving all path spellings. That is not a quality-gate claim; it is evidence
that the old visual/no-op symptom was at least partly caused by duplicate
replacement segment lacing, and that projecting `sequence -> emitted segment ID`
removes that specific no-op mechanism.

A full true-SPOA C4 visual rerun was not used as final evidence in this audit.
The directory `data/c4_spoa_dedup_probe_20260604T000000Z` existed during the
audit but contained only logs and no graph/report/render artifacts at inspection
time. The small committed C4 slice and targeted renderer tests were therefore
the practical validation surface.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Root-cause audit | 0.95 | The code path from replacement graph construction to render/lace-in was traced with line references, and the old failure mode matches prior C4 audit evidence. |
| Patch mechanism | 0.90 | The seq-to-ID projection is implemented in the renderer and reconnects paths/links through canonical emitted IDs. |
| Targeted tests | 0.90 | Focused tests cover both duplicate replacement sequences and path-start insertion ordering. |
| C4 validation | 0.78 | A real C4 slice rerun validates path preservation and shows reduced duplicate extras; full true-SPOA visual validation remains pending. |
| Constraint fidelity | 0.92 | No quality/metric gate was added; path spelling remained the only hard C4 acceptance gate. |
| Reporting transparency | 0.90 | Commands, counts, residual risks, and the incomplete full-rerun evidence are recorded explicitly. |

Overall: `0.86`.
