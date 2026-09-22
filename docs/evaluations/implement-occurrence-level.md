# Evaluation: implement-occurrence-level

Task: `implement-occurrence-level`
Evaluator: `agent-322`
Date: 2026-05-30

## Grade

Overall score: **0.20 / 1.00**
Confidence: **0.78**
Rubric underspecified: **false**

The task did not have an explicit weighted rubric, but the implementation notes and validation checklist are specific enough to grade against. The low score is driven by failure to implement the core requested behavior: occurrence-level shared-run filtering. The current code still filters at syncmer-node level and then clones every occurrence of a weak node.

Status caveat: at evaluation time, `wg show implement-occurrence-level` reported the task as `in-progress`, and the branch head before this evaluation artifact matched `origin/eg/c4-crush-resolution-controls`. I did not find a task-specific implementation commit for `implement-occurrence-level`; this grade is therefore based on the available branch state and workgraph record.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Core occurrence-level filtering | 0.00 | No stable occurrence identity is computed for support. `SharedContextFilterResult` still stores `weak_nodes: FxHashSet<u32>`, and support is represented as node sets. |
| Use of existing local walks | 0.45 | Existing `LocalSyncmerWalkStep` walks are collected for masking, but only `node` and `pos` are used to produce node-level support, not occurrence-level decisions. |
| Raw/blunt private materialization integration | 0.35 | Raw and exact-blunt writers can materialize private/cloned per-occurrence segments, but they are driven by a node-level `split_syncmers` set, so supported occurrences of a weakly isolated node would also be cloned. |
| Frequency masking and local-repeat preservation | 0.60 | Existing frequency masking and local-repeat clone logic appears preserved, but this is legacy behavior rather than a task-specific occurrence-level implementation. |
| `:nomask` behavior | 0.35 | The existing disabled mask path likely prevents the shared-context filter from running, but no new or adjusted task-specific test proves this. |
| Required tests | 0.05 | Existing tests assert node-level behavior, including "clones all weak node occurrences"; the required test where only an isolated/off-diagonal occurrence is cloned was not added. |
| Validation evidence | 0.00 | Required validation commands could not be completed in this worktree because vendored source files are missing. |

## Evidence

The current implementation remains node-level:

- `src/commands/syng2gfa.rs:1265` defines `weak_shared_run_syncmers(walks: &[Vec<u32>], ...) -> FxHashSet<u32>`, so the support result is a set of node ids, not occurrences.
- `src/commands/syng2gfa.rs:1330` defines `SharedContextFilterResult { weak_nodes: FxHashSet<u32>, ... }`.
- `src/commands/syng2gfa.rs:1370` and `src/commands/syng2gfa.rs:1386` return node sets for run-supported and sequence-supported syncmers.
- `src/commands/syng2gfa.rs:1452` computes `weak_nodes` as `shared_nodes.difference(&supported)`, which means one supported occurrence of a node globally blesses that node.
- `src/commands/syng2gfa.rs:1579` defines `split_syncmer_clone_occurrences(work, split_syncmers: &FxHashSet<u32>)`; `src/commands/syng2gfa.rs:1586` then clones every raw occurrence whose absolute node id is in `split_syncmers`.
- `src/commands/syng2gfa.rs:2276` and `src/commands/syng2gfa.rs:2444` assign `split_syncmers` from `context_filter.weak_nodes`, preserving the node-level path into both exact-blunt and raw range emission.
- `src/commands/syng2gfa.rs:3416` has an existing test named `test_sequence_context_split_clones_all_weak_node_occurrences`, which explicitly asserts the behavior this task was meant to replace.

This conflicts directly with the acceptance requirement:

> Compute weak occurrences as `(node, path/range-local position or another stable occurrence identity)`, not just weak nodes.

It also conflicts with the required test scenario:

> node X appears inside a valid run in one pair and as an isolated/off-diagonal occurrence elsewhere; only the isolated occurrence is cloned/private.

The present code would either keep all occurrences shared if the node is globally supported by a repeated run, or clone all occurrences if the node is globally weak. It cannot express "clone only this occurrence of node X."

## Validation Attempted

- `cargo fmt --check` failed before formatting checks because `vendor/gfaffix/src/main.rs` is missing.
- `source ./env.sh && cargo test test_sequence_context_split_clones_all_weak_node_occurrences --lib -- --nocapture` failed during the `impg` build script because multiple `vendor/syng/*.c` sources are missing, including `vendor/syng/syngbwt3.c`.

I did not run the broader requested test matrix after these workspace-level failures, because the same missing vendored files block compilation.

## Calibration

A score of 0.20 gives credit for pre-existing infrastructure that would be useful for an occurrence-level implementation: local walks are collected, raw path work tracks per-syncmer occurrences, and raw/blunt emitters can already materialize private clone segments. However, the actual acceptance criterion is the local run/sequence-k occurrence filter, and that is absent. The implementation would still allow the original reported failure mode: a node supported in one good repeated run can bless unrelated off-diagonal occurrences elsewhere.
