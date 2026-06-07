# Local Compression Iteration 2 Synthesis

Date: 2026-06-07

Task: `synthesize-next-local`

## Source Reports

This synthesis integrates exactly two predecessor analyses:

- Algorithm candidate analysis:
  `docs/evaluations/local-compression-autopoietic/iter-2-candidate-analysis.md`.
  Its decision is to implement one chunked SweepGA/seqwish candidate by feeding
  the existing sorted chunk-window selector, with `min_invariant_gap_bp=5`, into
  the current `seqwish_k=19` resolver (`iter-2-candidate-analysis.md:160-191`).
- Metric blind-spot analysis:
  `docs/evaluations/local-compression-autopoietic/next-metric-blind-spots.md`.
  Its highest-value metric recommendation is
  `path_replay_compression_ratio = total_spelled_path_bp / total_segment_bp`,
  introduced first as a diagnostic rather than a hard gate
  (`next-metric-blind-spots.md:89-118`).

## Integrated Finding

Iteration 1 proved a boundary-selection point, not a compression-quality point.
The sorted chunk-window branch is the only current fast-profile method that
passes `nested_top_level_wrong` with two candidates, two bubbles, and two
flubbles, but it is larger and deeper than both the one-parent placeholder
family and the internal PGGB/SmoothXG controls. The candidate analysis therefore
argues that the next resolver experiment must inherit the sorted chunk-window
boundary before spending work on broad `k` or min-match sweeps.

The metric analysis adds the missing evaluation pressure: exact path spelling
and binary topology pass/fail no longer rank the rows well enough. Current
controls preserve paths and are compact, but fail synthetic topology assertions;
current local methods satisfy more topology checks, but can be over-expanded.
`path_replay_compression_ratio` is the smallest useful diagnostic because it is
computed from GFA path replay and segment bases rather than trusting
`# testbed_metrics` comments.

## Next Hypothesis

If the next loop implements a metric-instrumented chunked SweepGA/seqwish
candidate, then the fast profile should keep the iteration-1 boundary win on
`nested_top_level_wrong` while producing enough replay-compression evidence to
tell whether the resolver improves sharing or merely preserves a topology-shaped
expanded graph.

Expected support for the hypothesis:

- The new candidate preserves exact paths on graph-producing fast rows.
- On `nested_top_level_wrong`, it reports `candidate_count=2`,
  `bubble_count=2`, and `flubble_count=2`.
- The fast scoreboard includes `path_replay_compression_ratio` for graph rows.
- The replay-ratio table is compared against the one-parent rows and the three
  internal controls, with the ratio treated as diagnostic in this iteration.

## Recommended Next Loop Action

Implement one metric-instrumented chunked SweepGA/seqwish candidate:
`chunk_window_sweepga_seqwish`.

The action is one implementation task with one purpose: evaluate the best known
boundary selector under the current SweepGA/seqwish resolver while adding the
diagnostic metric needed to judge the result. The implementation should reuse
the sorted chunk-window selector from `chunk_window_smooth_or_crush`, keep
`min_invariant_gap_bp=5`, keep `seqwish_k=19`, and add
`path_replay_compression_ratio` to the scoreboard as a reported diagnostic. It
should not change runner defaults or hide existing local/control rows.

Rationale:

- The candidate report identifies the sorted chunk-window boundary as the only
  current fix for the nested-parent overmerge fixture, and specifically
  recommends using that boundary before resolver parameter sweeps.
- The metric report shows that topology pass/fail alone can over-reward
  synthetic structure, so the next algorithm row needs replay-compression
  evidence in the same run.
- Combining the candidate and metric in one implementation gives the next
  synthesis a direct answer: whether chunked SweepGA/seqwish preserves the
  two-window topology win and improves sequence sharing relative to the
  one-parent rows and path-copy-like controls.

## Validation Commands For Next Implementation

The next implementation task should add focused tests first, then run these
commands after implementation:

```bash
python3 -m py_compile scripts/local_compression_testbed.py
```

```bash
cargo test --test test_local_compression_testbed local_compression_path_replay_compression_ratio -- --nocapture
```

```bash
cargo test --test test_local_compression_testbed local_compression_chunk_window_sweepga_seqwish_nested_top_level_wrong -- --nocapture
```

```bash
python3 scripts/local_compression_testbed.py run --profile fast --manifest tests/test_data/local_compression/manifest.json --out-dir docs/evaluations/local-compression-testbed-fast
```

```bash
python3 - <<'PY'
import csv
from pathlib import Path

scoreboard = Path("docs/evaluations/local-compression-testbed-fast/scoreboard.tsv")
rows = list(csv.DictReader(scoreboard.open(), delimiter="\t"))
assert rows, "scoreboard has no rows"
assert "path_replay_compression_ratio" in rows[0], "missing replay-compression metric"

target = [
    row
    for row in rows
    if row["fixture_id"] == "nested_top_level_wrong"
    and row["method_id"] == "chunk_window_sweepga_seqwish"
]
assert len(target) == 1, f"expected one target row, found {len(target)}"
row = target[0]
assert row["exact_path_preservation"] == "pass", row
assert row["expected_topology_status"] == "pass", row
assert row["candidate_count"] == "2", row
assert row["bubble_count"] == "2", row
assert row["flubble_count"] == "2", row
ratio = float(row["path_replay_compression_ratio"])
assert ratio > 1.0, row

comparison_methods = {
    "local_syng_crush_sweepga",
    "top_flubble_nonoverlap_sweepga",
    "chunk_window_smooth_or_crush",
    "chunk_window_sweepga_seqwish",
    "pggb_control",
    "smoothxg_control",
    "pggb_plus_smoothxg_control",
}
for comparison in rows:
    if (
        comparison["fixture_id"] == "nested_top_level_wrong"
        and comparison["method_id"] in comparison_methods
    ):
        print(
            comparison["method_id"],
            comparison["exact_path_preservation"],
            comparison["expected_topology_status"],
            comparison["candidate_count"],
            comparison["bubble_count"],
            comparison["flubble_count"],
            comparison["path_replay_compression_ratio"],
        )
PY
```

```bash
cargo build
```

```bash
cargo test
```

```bash
git diff --check
```

## Current Task Validation

- This synthesis cites both predecessor reports by path and decision lines.
- It recommends exactly one next loop action:
  `chunk_window_sweepga_seqwish`.
- It includes explicit validation commands for the next implementation task.
- This task only writes a synthesis report under
  `docs/evaluations/local-compression-autopoietic/`; it does not modify runner
  defaults, fixtures, generated scoreboards, or testbed code.
