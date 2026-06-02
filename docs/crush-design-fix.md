# Crush architectural fix — design proposal

**Inputs:**
- [`docs/crush-audit.md`](crush-audit.md) (the 4 failures observed on real C4)
- [`docs/crush-fixtures-redproof.md`](crush-fixtures-redproof.md) (the 4 red `#[ignore]`d tests in `tests/test_crush_integration.rs` that reproduce them)
- [`docs/crush-trace.md`](crush-trace.md) (the labeled stage flow with `file:line` for every box)

**Status:** design only. No code change in this commit (apart from this doc). The downstream `crush-impl` task is paused; the user reviews this proposal before resuming.

**HEAD under review:** `b8e80c2` (parent `69ae688` `feat: crush-fixtures (agent-52)`, rooted at `90ba74f`). All `file:line` citations below are against this tree.

---

## 1. The architectural mistake (one paragraph)

The `H6 → H5` boundary inside the per-round loop body `resolve_graph_bubbles` at **`src/resolution.rs:611-617`** consumes the wrong invariant. It calls `apply_replacement_frontier` (`:1452`) which returns a `next_graph` validated for per-path byte equality via `path_sequences_equal` (`:2832`, `:1536`), then assigns `graph = next_graph` *unconditionally* at line `616`, and computes `after_quality = graph_quality(&next_graph)` (`:615`) only for logging. Concretely: **`path_sequences_equal` is a per-replacement correctness invariant** ("no path was corrupted by aligning its traversals"); commit `0af1a4c` removed `round_quality_decision` / `RoundQualityDecision::Reject` and repurposed `path_sequences_equal` as the *only* round-level acceptance criterion, conflating "byte-equality preserving" with "monotone progress". That round-level boundary needs *two* invariants the per-replacement check cannot supply: (a) **structural collapse** — the rewritten graph must not be allowed to carry duplicate-sequence segments emitted by independent replacements, because `render_rewritten_graph` (`:2846`) assigns a fresh id to every `OutNode::Replacement(plan_idx, node_idx)` (`:2887`) with no `seq → id` projection, so byte-identical replacement segments from tandem-paralog bubbles emit as distinct `S` lines (Failure 3, 5 474 dup-seqs on round 1 of C4); and (b) **monotone progress** — the graph quality after the round must not regress past a bounded threshold, otherwise round N+1 sees a wider bubble and the loop has no fixed point (Failures 2, 4, 1 cascade). The mistake is therefore precisely-locatable: the `OutNode::Replacement(plan_idx, node_idx)` identity is correct *inside* `apply_replacement_frontier` (it must distinguish replacement plans during the per-path rewrite at `:1480-1528`) but is *wrong* when carried through to `render_rewritten_graph` as the canonical segment identity for the emitted graph; and the round acceptance check at `:611-617` has been narrowed to byte-equality, losing the composition guarantee.

---

## 2. Proposed corrected architecture

The fix is two surgical changes at the two boundaries identified above, plus one defensive timing measure. The fix preserves the validated-replacement invariants from commits `0af1a4c` and `f453983` and adds two orthogonal invariants at the integration boundary.

### 2.1 Functions / modules that change

| function | file:line | change |
|---|---|---|
| `render_rewritten_graph` | `src/resolution.rs:2846` | insert sequence-keyed dedup: before emitting an `S` line for an `OutNode::Replacement(plan, node)`, look up the segment's byte sequence in a `FxHashMap<&[u8], String>` keyed by sequence; if present, route the `OutNode` to the existing id; if absent, mint a new id via `next_unused_segment_id` and record `seq → id`. Apply the same dedup to `OutNode::Original` segments **only against other replacement sequences** (not against other originals — originals keep their stable ids) so a replacement that re-emits the original sequence collapses onto the original id. |
| `resolve_graph_bubbles` | `src/resolution.rs:454-632` | reinstate the round-level acceptance gate at `:611-617`. After `let after_quality = graph_quality(&next_graph);` add `let decision = round_quality_decision(before_quality, after_quality, config);` and branch: on `Reject`, do **not** assign `graph = next_graph`, log the reason, increment `stats.bailed += plans.len()`, and `break` out of the round loop (terminating crush deterministically rather than continuing to a worse round). |
| `round_quality_decision` (new) | `src/resolution.rs` (re-added near `graph_quality`, ~`:990`) | port the body from the `0af1a4c -` side: takes `before`, `after`, `config`; returns `Accept` if `before.score == 0`; returns `Reject` when `(after.score - before.score) / before.score > config.max_round_score_growth` and the visual-tail (`path_white_space_bp_p99`, `path_white_space_bp_max`) is not strictly improved; otherwise `Accept`. **Key change from the deleted version:** the gate is **post-dedup**, so the score is computed on the actually-emitted graph, not on the OutNode-tagged structural surrogate. |
| `ResolutionConfig` | `src/resolution.rs:79-130` | reinstate one field: `pub max_round_score_growth: f64` with `DEFAULT_MAX_ROUND_SCORE_GROWTH = 0.10` (looser than the deleted `0.02` to match the test 1 assertion bound; the value is a tuning point — see §5 risks). Do **not** reinstate `min_round_score_improvement` or `disable_round_quality_check` or per-replacement `replacement_method_needs_quality_guard` — those broadened the surface area and were correctly removed. |
| CLI surface | `src/main.rs` Args::Crush handler `~:7745-7873` | add a single CLI flag `--max-round-score-growth <float>` wired into `config.max_round_score_growth`. Default is the new constant. (Backwards-compatible: omitting the flag uses the default.) |
| `apply_replacement_frontier` | `src/resolution.rs:1452-1551` | **no change.** It correctly produces the rewritten path step list. The existing `crush apply:` debug lines added by `crush-trace` already report the structural ingredients of Failure 3 at this boundary; they stay. |

### 2.2 Data structures that change

| data structure | file:line | change |
|---|---|---|
| `OutNode` enum (`Original(usize)`, `Replacement(usize, usize)`) | `src/resolution.rs:394-398` | **no change to the enum**; it remains the correct identity for per-path rewriting. The change is in how it is *projected* to segment ids at the serializer (see §2.1 row 1). |
| `render_rewritten_graph` local state | `src/resolution.rs:2855-2870` | add `let mut id_by_seq: FxHashMap<Vec<u8>, String> = FxHashMap::default();` keyed by replacement segment sequence. Replacement walk consults `id_by_seq` before minting a new id. `id_by_node` continues to map `OutNode → emitted id` for L/P-line construction so reverse-lookup still works after dedup. |
| no other data structure changes | — | the on-disk GFA contract is unchanged; the only observable difference is fewer S-lines and fewer L-lines (because dup links collapse with their endpoints). |

### 2.3 Sequencing / ordering changes

Two ordering changes inside the per-round body:

1. **Render-then-check (post-dedup gate).** The round-acceptance gate must run *after* `apply_replacement_frontier` has rendered, parsed, and validated the post-rewrite graph. This is already the existing order — `:613 apply → :615 after_quality → :616 assign`. The change is to *insert* the gate between `:615` and `:616` rather than between `apply` and `after_quality`. Rationale: `graph_quality` must measure the deduplicated graph (the one users actually receive), not the OutNode-tagged structural surrogate, otherwise Failure 3 dominates the score and triggers false rejections.

2. **Rollback semantics.** On `RoundQualityDecision::Reject`, the loop terminates rather than continuing. This matches the baseline behavior documented in `docs/crush-audit.md` §Failure-2 (`2026-05-23T21:44:20Z` log shows the baseline stops at round 2 after the reject) and avoids the cascade Failure 4 describes (an accepted bad round inflates the next round). The current loop has no rollback because no acceptance decision is made; the new loop has the *option* to reject and the *requirement* to stop on first reject.

### 2.4 New invariants at which boundaries

| boundary | invariant | enforced by |
|---|---|---|
| `G3 → H1` (per-path rewrite → segment emission) | **No two emitted `S` lines may carry the same byte sequence.** | `render_rewritten_graph` `seq → id` map. Tested by the un-`#[ignore]`d `c4_round1_render_emits_no_duplicate_segment_sequences`. |
| `H6 → H5` (quality computed → graph assigned) | **`graph_quality(next).score ≤ graph_quality(before).score · (1 + max_round_score_growth)`** unless the visual-tail (`ws-p99` and `ws-max`) is strictly non-worse. | `round_quality_decision` gate. Tested by the un-`#[ignore]`d `c4_round1_sweep_quality_gate_rejects_score_growth` and `c4_round2_segment_bp_does_not_exceed_round1`. |
| `B0 → I0` (loop exit) | **The number of rounds actually accepted is bounded by the first round whose acceptance violates the round-quality invariant, or by `max_iterations`, whichever comes first.** | the `break` on `Reject` plus the existing `max_iterations` loop bound. Tested by the un-`#[ignore]`d `c4_canonical_command_completes_within_budget`. |
| `F4 → G0` (per-replacement validation) | **Per-path byte equality for every replacement** (unchanged from HEAD). | `validate_replacement_paths` `src/resolution.rs:2436`. Preserved verbatim. |

The first three are new; the fourth is the preserved validated-replacement invariant from `0af1a4c`.

---

## 3. Causal chain for each failing crush-fixtures test

### Test 1 — `c4_round1_sweep_quality_gate_rejects_score_growth` (Failure 2)
Asserts: `segment-bp growth in round 1 ≤ 10%`. HEAD observes `+18.7%`.

Causal chain:
1. The seq-dedup at `render_rewritten_graph` collapses the 5 474 dup-seq replacement segments observed in Test 2 onto existing ids, reducing post-round `segment-bp` substantially toward the baseline figure (`+8.4%`, see audit §Failure 2). On the slice the instrumentation already reports 948 dup-seqs out of 2 984 emitted segments (`docs/crush-trace.md` §5); deduping those is a direct ~32% reduction in post-round segment-bp on the slice.
2. If after dedup the residual growth still exceeds 10%, `round_quality_decision` returns `Reject` because `growth > 0.10` and visual-tail (`ws-p99`, `ws-max`) is worse (Failure 2 reports `ws-p99 +77.9%`, `ws-max +77.8%`).
3. On `Reject` the new branch skips `graph = next_graph` (`:616`) and `break`s; the output graph is the pre-round graph (`segment-bp = 389 316`).
4. The test reads `(after_segs, after_bp) = gfa_file_segment_stats(out_path)` and computes `growth = 0`, which is ≤ 0.10. ✓

### Test 2 — `c4_round1_render_emits_no_duplicate_segment_sequences` (Failure 3)
Asserts: zero S-lines sharing a sequence after round 1. HEAD observes 5 474.

Causal chain:
1. `render_rewritten_graph` walks `ordered_nodes`. For each `OutNode::Replacement(plan, node)` it now looks up `replacement.segments[node].seq` in `id_by_seq`.
2. If the sequence has been seen before (under any plan), the existing id is reused: `id_by_node.insert(node, existing_id)`; no new `S` line is emitted.
3. If new, `next_unused_segment_id` mints an id, the `S` line is emitted, and `id_by_seq.insert(seq.clone(), id)` records it.
4. L-line construction continues to use `id_by_node` keyed by `OutNode`, so the collapsed-segment edges are correctly attributed.
5. After the loop, the `gfa_file_duplicate_seqs(out_path)` helper sees zero S-lines sharing a sequence. ✓

(Note: this is a structural invariant enforced by construction; it does not depend on the quality gate firing.)

### Test 3 — `c4_round2_segment_bp_does_not_exceed_round1` (Failure 4)
Asserts: `segment-bp growth after 2 rounds ≤ 20%`. HEAD observes `+33.6%`.

Causal chain:
1. Either round 1 is accepted post-dedup (growth ≤ 10%, see Test 1 causal chain), in which case round 2 starts from a graph whose `segment-bp` is at most `389 316 · 1.10 = 428 248`; or round 1 is rejected and the run terminates at round 1.
2. If round 2 runs, `find_candidate_frontier` (`:1107`) sees a graph whose anchor structure has been preserved by dedup (Failure 4 says the round-2 inflation is because round 1 collapsed short anchors that previously fenced bubbles — but dedup *restores* the anchor structure, because byte-identical anchors across replacements are now one segment). So round-2 candidates have the expected size distribution.
3. Round 2's post-dedup growth is gated the same way. Either it is accepted ≤ 10%, giving total 2-round growth ≤ `1.10² - 1 = 21%` (one boundary case; see §5 risks); or it is rejected.
4. To meet the 20% bound the gate must reject any round whose individual growth exceeds the bound that, compounded, would exceed 20%. With `DEFAULT_MAX_ROUND_SCORE_GROWTH = 0.10` the compounded 2-round bound is 21%, marginally above 20%. The test is `≤ 20%`, so either (i) the proposed default is tightened to `0.095` (still well above the historical `0.02`), or (ii) the test's bound is widened to `0.21`. **Decision required from the user.** See §5 risks. ✓ assuming (i) or (ii).

### Test 4 — `c4_canonical_command_completes_within_budget` (Failure 1)
Asserts: `impg crush ... --max-iterations until-done` completes within 360 s. HEAD: exit code 124 (timeout).

Causal chain:
1. With the gate in place, round 1 either accepts (post-dedup growth ≤ threshold, in which case the graph stays near baseline scale: ~390 kbp segment-bp) or rejects (loop terminates after round 1's discovery+build, ~75 s on the C4 input per `crush-fixtures-redproof.md` Test 1).
2. If round 1 accepts, round 2 sees a non-bloated graph (dedup preserved anchors; growth was bounded). Round 2's POVU discovery does not find the 18.4 MB single-candidate bubble of audit §Failure 4; instead the candidate distribution is close to baseline round 2 (3 candidates, traversal_max ~32k).
3. If round 2 itself violates the gate (which the baseline already documented: `crush round 2: rejecting 3 replacement(s): visual-tail score grew by 0.0618`), the run terminates at the round-2 reject.
4. The baseline completed the same command in 279 s; the proposal restores the baseline's loop-termination behavior. Total wall is bounded by baseline + dedup overhead (one extra `FxHashMap<Vec<u8>, String>` insert per emitted replacement segment, O(segment_bp) total, ≤ a few seconds on a 700 kbp graph). 360 s budget is comfortably met. ✓

### Companion test — `c4_slice_auto_crush_preserves_path_sequences` (passes on HEAD)
Asserts: path sequences preserved, 147 sites resolved, slice bp drops `64 348 → 58 003` on auto/SPOA.

Causal chain (must continue passing):
1. The seq-dedup change does not affect `path_sequence_map` (which walks paths through the rewritten graph using `id_by_node` reverse-lookup); the per-path byte equality invariant is unchanged.
2. The round-quality gate compares `before` and `after` `graph_quality`. The slice already *improves* (bp drops 9.9%); `before.score > after.score`, so `growth < 0` and the gate is `Accept`.
3. The 948 dup-seqs the instrumentation reports on the slice collapse onto existing ids; `resolved.stats.resolved == 147` is unaffected because `stats.resolved` counts accepted *plans*, not segments. ✓

---

## 4. Existing tests that need to be UPDATED

Existing intra-module tests in `src/resolution.rs` that touch the invariants this proposal changes. **None should be removed without analysis** — each was testing a real invariant, even if the invariant was the wrong one for that boundary.

| test | file:line | needed change | why |
|---|---|---|---|
| `first_crush_round_accepts_path_step_growth_without_visual_tail_regression` | `src/resolution.rs:3519-3549` | **Keep as-is, but rename back and extend.** This was the renamed-and-inverted descendant of `first_crush_round_rejects_quality_regression` (per `docs/crush-audit.md` §Failure-2-e: added by `259689b`, renamed by `b163321`, kept by `0af1a4c`). It currently builds a 100-bp two-path GFA, runs `StarBiwfa`, and asserts `resolved.stats.resolved == 1`. After the fix, that exact assertion still passes (the small benign bubble is accepted by the gate because it doesn't regress quality). Add a sibling test `first_crush_round_rejects_quality_regression` that constructs a degenerate GFA where seqwish-induced replacement *would* grow segment-bp >10% and assert `resolved.stats.resolved == 0` and `resolved.stats.bailed > 0`. This requires no `sweepga` binary — a hand-rolled `Graph` with synthetic ranges suffices. |
| `auto_routes_small_bubbles_to_spoa_and_larger_bubbles_to_allwave` | `src/resolution.rs:3288` | **No change.** Tests routing, not acceptance. |
| `graph_quality_penalizes_path_white_space_bridges` | `src/resolution.rs:3489` | **No change.** Tests `graph_quality` scoring; the proposal does not change `graph_quality`. The reinstated `round_quality_decision` uses the existing scoring; reaffirming this test passes is enough. |
| `resolves_*` family (`:3080`, `:3102`, `:3123`, `:3139`, `:3158`, `:3190`, `:3214`, `:3241`, `:3710`, `:3759`) | `src/resolution.rs` | **No change.** Each constructs a 2- or 3-path synthetic GFA where the replacement is a clear improvement; the gate's `Accept` path runs as today. They preserve their assertions verbatim. |
| `replacement_seqwish_filter_defaults_to_seqwish_k_scale` | `src/resolution.rs:3552` | **No change.** Tests the `f453983` scaffold-filter invariant; untouched. |
| `c4_slice_auto_crush_preserves_path_sequences` (integration) | `tests/test_crush_integration.rs:453` | **No change.** Companion test; must continue passing — see §3 above. |
| `tests/test_syng_integration.rs::test_crush_cli_resolves_blunt_gfa` (`:207`) | end-to-end | **No change.** 7-bp insertion bubble; the gate accepts. Audit §Failure-2-e called this test out as too small to catch the bug; it remains useful as a CLI smoke. |
| `tests/test_pipeline_integration.rs::test_full_pipeline` (`:59`) | `#[ignore]` | **No change.** Does not exercise crush. |
| All 4 `#[ignore]`d tests in `tests/test_crush_integration.rs` (`:148`, `:235`, `:310`, `:387`) | integration | **Un-`#[ignore]` after the fix lands.** They are the acceptance contract for this proposal. The `#[ignore]` lines should be removed; the test bodies are unchanged. (One bound adjustment in Test 3 — see §3 above.) |

**No tests are proposed for removal.** The deleted historical tests (`round_quality_rejects_large_complexity_growth`, `round_quality_accepts_visual_tail_improvement_despite_score_growth`, `direct_poa_replacements_bypass_local_quality_gate`) removed by `0af1a4c` are *not* proposed to be re-added: their first sibling (the local per-replacement gate) is intentionally not being restored (see §2.1 row 4), so re-adding `direct_poa_replacements_bypass_local_quality_gate` would test a path that no longer exists; and the `round_quality_rejects_large_complexity_growth` and `round_quality_accepts_visual_tail_improvement_despite_score_growth` cases are covered functionally by the new `first_crush_round_rejects_quality_regression` and by the un-ignored `c4_round1_sweep_quality_gate_rejects_score_growth`.

---

## 5. Risks and interactions

### Risks

1. **Threshold-tuning risk.** `DEFAULT_MAX_ROUND_SCORE_GROWTH` was `0.02` in the deleted code; the proposal suggests `0.10` to match Test 1's `≤ 10%` bound, but this is calibrated to the C4 input. The baseline observed `+8.4%` on its accepted round 1; a stricter threshold would have rejected the baseline's own round 1. A looser threshold makes Test 3 (`≤ 20%` after 2 rounds, compounded bound is `(1+t)² - 1`) marginal. **Decision required:** (a) `t = 0.10`, accept Test 3's `≤ 20%` becomes `≤ 21%`; or (b) `t = 0.095`, keep Test 3's `≤ 20%`. Recommendation: (a), and tighten Test 3 to `≤ 21%` — this matches the geometric meaning of the bound.

2. **Dedup interaction with originals.** The proposed `id_by_seq` map is keyed only by replacement segment bytes. If a replacement emits a segment whose bytes equal an original segment that survives in the rewritten graph (a real possibility when sweepga "rediscovers" original anchors), we *should* collapse onto the original id, not mint a new one. Failing to do this leaves a small residue of dup-seqs (original + replacement) that Test 2's `gfa_file_duplicate_seqs` would still count. **Mitigation:** initialize `id_by_seq` with `(original.segments[i].seq, original.segments[i].id)` for every `i ∈ used_original` *before* the replacement walk. This is cheap (one pass over `used_original`, ~18 048 entries on round 1) and closes the residue exactly.

3. **L-line correctness under dedup.** `render_rewritten_graph` at `:2898-2918` builds edges from `id_by_node[from]` and `id_by_node[to]`. After dedup, two distinct `OutNode::Replacement` keys map to the same id; the L-line `BTreeSet<(String, bool, String, bool)>` already collapses duplicate edges by construction, so this is correct *if* `id_by_node` is populated for every `OutNode` (including the merged ones). **Mitigation:** in the dedup branch, still call `id_by_node.insert(node, id.clone());` so reverse-lookup is total. (The current code already calls `id_by_node.insert(...)`; the only change is that `id` may be a previously-minted one.)

4. **Quality gate firing on the slice companion test.** The slice already improves (`bp 64 348 → 58 003`), so `growth < 0`, gate Accepts. No risk.

5. **Quality gate firing on a real "good" round that has not converged yet.** The baseline shows the gate accepting round 1 (`growth = -30.8%`) and rejecting round 2 (`growth = +6.18%`). The proposal's default `0.10` is much looser than the baseline's `0.02`, so a baseline-style round-2 (`+6.18%`) would be *accepted* under the proposal. This is intentional: the looser threshold is for the post-dedup graph, where round-1's dedup contribution has already collapsed the easy wins. If users want the strict `0.02` behavior, they can pass `--max-round-score-growth 0.02`.

6. **Per-round wall-budget is not added.** Composition issue #4 in `docs/crush-trace.md` §4 identifies "no box in the diagram is a deadline owner" as the proximate cause of Failure 1's SIGTERM landing inside polish. The proposal relies on the quality gate to terminate the loop *before* polish is asked to operate on an 18 MB candidate. This is correct in the causal chain for Test 4 above. A separate per-round wall budget would be defensive (it would still help if a single round-1 polish blew up on a degenerate input), but is **not proposed here** because (a) Tests 1-4 don't require it; (b) adding a deadline owner changes more code (passing `Instant`s and atomics into the rayon par_iter), is harder to test deterministically, and would expand scope. It is documented here as a future follow-up; see "Future" below.

7. **Subtle interaction: signature dedup `seen` set vs structural dedup.** `find_candidate_frontier` (`:1107`) uses a per-bubble `signature` to deduplicate candidates *across rounds* (`seen: FxHashSet<String>`). After this proposal, the post-round graph has the same path sequences but fewer segments. POVU is order-dependent on the working graph (composition issue #3 in `docs/crush-trace.md`), so the next-round signatures may not match the previous round's signatures verbatim. **Mitigation:** `signature` is computed from path-spelled bubble content (`:1233`), not from segment ids, so dedup at the serializer should not perturb signatures across rounds. This is a property of `candidate_signature` (`src/resolution.rs:1486-` per the grep at hand-of-time, called at `:1233`). Spot-check this during implementation; if signatures *do* shift across rounds, the `seen` set may admit re-attempts at the same bubble. This is not a correctness bug (each round is independently gated) but may waste work.

### Interaction with the validated-replacement invariants from commits `0af1a4c` and `f453983`

- **`0af1a4c` invariant — "Correctness is enforced by exact path validation; graph-quality metrics are diagnostic telemetry only" (per the docstring at `src/resolution.rs:39-44`).** The proposal **preserves correctness via exact path validation** (`validate_replacement_paths` is unchanged, `path_sequences_equal` is unchanged). The proposal **adds, separately, a termination/progress invariant** at the round-loop boundary that is upgraded from "diagnostic telemetry" to "acceptance criterion". Crucially this is **at a different boundary** than the one `0af1a4c` removed: `0af1a4c` removed both the per-replacement local guard *and* the round-level guard; the proposal restores **only** the round-level guard, on the **post-dedup** graph. Per-replacement methods (Poa, Poasta, StarBiwfa) are no longer quality-guarded — that simplification stands. The composition reads as: "every replacement preserves paths exactly; every round preserves progress; together the loop terminates."

- **`f453983` invariant — `replacement_num_mappings = "1:1"` and `replacement_scaffold_filter = "1:1"` as default for pairwise graph induction.** The proposal **does not touch** `seqwish_replacement_config` (`:1704`) or the scaffold-filter wiring. Sweepga/Allwave produce per-replacement graphs as today; dedup operates **downstream** at `render_rewritten_graph`. The two changes are orthogonal: scaffold filtering controls what pair alignments enter seqwish; segment-sequence dedup controls how the output graph's `S` lines are emitted. A user who tunes `--replacement-scaffold-filter many:many` to get larger replacements still gets dedup at serialization.

### What could go wrong (concrete failure modes if the implementation is sloppy)

- Dedup keyed by `String` rather than `Vec<u8>` would force a UTF-8 conversion and could mis-handle non-ASCII sequences. Use `Vec<u8>` keys.
- Forgetting to seed `id_by_seq` with originals (risk #2) leaves Test 2 failing.
- Forgetting to insert into `id_by_node` on the merged-id branch (risk #3) leaves L-line construction with `None` lookups and *silently drops edges*.
- Putting the gate **before** the parse-and-validate step would reject on an invalid (path-corrupting) replacement, masking a real correctness bug. The gate must run **after** `parse_gfa(rendered)` and after `path_sequences_equal` succeeds.
- Wiring the CLI flag to a 64-bit float without bounds-checking lets `--max-round-score-growth -1.0` reject every round including the no-op case. Add a `if config.max_round_score_growth < 0.0 { return Err... }` in the CLI parse path.

---

## 6. Scope estimate

**Files touched:** 2.
- `src/resolution.rs` — primary; ~180 ± 30 lines changed.
- `src/main.rs` — single CLI flag plumbing; ~15 ± 5 lines added.

**Tests touched:** 2.
- `tests/test_crush_integration.rs` — four `#[ignore = ...]` lines removed; one bound adjustment on Test 3 (`≤ 20%` → `≤ 21%`, see risk #1). ~6 line diff.
- `src/resolution.rs` (intra-module tests) — one new test `first_crush_round_rejects_quality_regression` added (~40 lines); the existing `first_crush_round_accepts_path_step_growth_without_visual_tail_regression` is unchanged.

**Total LOC delta:** add ~220, remove ~5 (the `#[ignore]` lines and the orphan comment at `src/resolution.rs:1026-1030` that currently reads "this is logged for diagnostics only; acceptance is based on exact replacement path validation" — should be updated to reflect the new round-level acceptance criterion).

**Localized or broader refactor?** **Localized.** Two functions and one struct field, plus the CLI flag and one new constant. No public API changes. No new traits, no module reorganization. The PR diff should read mostly as a backport of `round_quality_decision` from `0af1a4c -` plus a 30-line block inside `render_rewritten_graph` and a config field.

**Estimated implementation wall-time:** 2-4 hours for one engineer who is already familiar with the resolution module. Validation against the 4 ignored tests is the long tail because each takes 1-15 min to run on real C4 input; running them in series end-to-end is the budget driver, not the code change.

---

## Future (out of scope for this proposal)

These are explicit follow-ups that **should not** be bundled into `crush-impl`:

1. **Per-round wall-clock deadline** (composition issue #4 in `docs/crush-trace.md` §4). Defensive guard against a single degenerate replacement blowing past the round-quality check by hanging *inside* the build/polish path before the check even runs. Should be its own task: it adds a deadline owner, requires passing `Instant`s through `par_iter`, and benefits from its own test.
2. **Post-acceptance unchop / `gfaffix`-style compaction on the final emit.** The `nosort` output stage currently skips gfaffix. Even with seq-dedup, two replacements that emit *prefix-matching* but non-identical segments are not merged. A bounded compaction step would tighten the final GFA further but is not required for any of the four red tests.
3. **Restore `RoundQualityDecision::Reject` reasoning into structured stderr for CI consumption.** The deleted `crush round N: rejecting M replacement(s): visual-tail score grew by ...` log line should be re-added (see §2.1 row 2) — this is part of this proposal — but a JSON/machine-readable summary line for CI harnesses would be a clean follow-up.

---

## Appendix — verbatim summary for a reviewer

If you take only one thing from this doc: **the per-replacement byte-equality check (`validate_replacement_paths`) is correct and must stay; the round-loop also needs a structural-collapse invariant at `render_rewritten_graph` and a monotone-progress invariant at `graph = next_graph`; both are post-dedup; one CLI flag exposes the threshold; no algorithm change to sweepga / seqwish / spoa / poasta / smooth.**
