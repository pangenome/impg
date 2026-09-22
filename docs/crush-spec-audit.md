# Crush implementation vs. canonical specification — audit

**Spec:** [`docs/crush-architecture-spec.md`](crush-architecture-spec.md) (HEAD = `4545cf8`).
**Audit HEAD:** `4545cf8` (worktree `wg/agent-62/crush-spec-audit`). Source tree differs
from `main` only in `docs/`; the live `src/resolution.rs` audited here is identical
to the one traced in [`docs/crush-trace.md`](crush-trace.md) and behaviourally
described in [`docs/crush-audit.md`](crush-audit.md). No code change was made by
this audit; the pre-existing debug log lines inside `apply_replacement_frontier`
(committed by `wg/agent-55/crush-trace`) were used as-is to surface invariant 1
on the canonical slice — no new instrumentation was added.

**Real C4 fixtures used for invariant checks (slice scale, per task instructions —
no full GRCh38 perf run):**

| name | command | output |
|---|---|---|
| slice-auto-1round   | `impg crush -g tests/test_data/crush/c4_slice_1500_3000.gfa --method auto    --max-iterations 1` | `/tmp/crush-audit-62/slice_auto_1round.gfa` (5 388 459 bytes, 2 984 segments, 58 003 bp) |
| slice-sweepga-1round| `impg crush -g tests/test_data/crush/c4_slice_1500_3000.gfa --method sweepga --max-iterations 1 --seqwish-k 311 --min-traversal-len 100` | `/tmp/crush-audit-62/slice_sweepga_1round.gfa` (4 289 segments, 178 933 bp) |
| slice-auto-3round   | `impg crush -g tests/test_data/crush/c4_slice_1500_3000.gfa --method auto    --max-iterations 3` | `/tmp/crush-audit-62/slice_auto_3round.gfa` (2 773 segments, 47 509 bp) |

All three were produced by the binary built from this worktree
(`target/release/impg` at HEAD `4545cf8`, built 2026-05-24T19:57Z; `git submodule
update --init`, `CPATH`/`LIBRARY_PATH` pointed at the gnu_old htslib and jemalloc
stores to satisfy the wfmash transitive build).

## Phase-by-phase audit

For each phase the audit answers, in order:

1. Present in the code?
2. If yes, where (file:line entry point + key functions)?
3. What does the code actually do (plain language)?
4. Where does it diverge from the spec?
5. Is the divergence a localized bug, a partial implementation, or "phase entirely missing"?

### Phase 1 — Bubble detection (POVU)

1. **Present:** yes.
2. **Where:** `find_candidate_frontier` `src/resolution.rs:1107`, called once per round at `src/resolution.rs:458`. The POVU side is in the `povu-rs` dep: `NativeGfa::parse` is invoked at `src/resolution.rs:1142` and `decompose_flubbles` at `src/resolution.rs:1152` (impl: `povu-rs/src/native_gfa.rs:446`).
3. **What it does:** renders the current working graph to a GFA string (`render_graph` `src/resolution.rs:1133`), parses it as a POVU `NativeGfa`, then asks POVU to decompose flubbles using path 0 as the root. POVU returns `FlubbleDecomposition.sites` which is a `Vec<FlubbleSite>` carrying every site at every level (`{id, parent_id, level, reference_start_step, reference_end_step, start, end, is_leaf}`) sorted by reference start/end (`povu-rs/src/native_gfa.rs:472`). Crush then materialises one `BubbleCandidate` per site that has ≥2 distinct path anchors (`unique_anchor_range` `:1395`), computes traversal statistics, and runs greedy non-overlapping selection against an `occupied_by_path` interval set (`:1270-1281`).
4. **Divergence from spec:** the spec says "identify **top-level non-overlapping flubbles**". The code does not consume POVU's `level` or `is_leaf` field — it consumes **all** sites and produces non-overlap by greedy selection on a per-path step interval set (`candidate_conflicts_with_occupied` `:1325`). This happens to produce a valid non-overlapping frontier but the priority order is `Allwave-before-Poa` then `estimated_step_savings` then `unique_steps` then `root_step_span` (`:1254-1262`), not "top-level first". A site nested inside another can be selected if its parent has lower priority. The path-step `occupied_by_path` set is also indexed by step rather than by reference coordinate, so two distinct sites that share a path-step run can each be rejected as overlapping with each other even if neither is parent of the other.
5. **Classification:** **partial implementation.** The output is a non-overlapping frontier, which the spec requires; but the frontier is not "top-level" in the level-hierarchy sense and the greedy selection's interaction with priority means the selected frontier is not deterministic in the same sense the spec implies. Localized fix.

### Phase 2 — Stratified resolution by bubble size

1. **Present:** yes, but the routing function does not implement the spec.
2. **Where:** routing is decided per-candidate by `candidate_selection_method` `src/resolution.rs:838-855`, dispatched through `candidate_replacement_method` `:1595` and the per-task call site `:519`. Default thresholds: `DEFAULT_AUTO_SPOA_MAX_TRAVERSAL_LEN = 2000` `:223`, `DEFAULT_MAX_TRAVERSAL_LEN = 10_000` `:218`, `DEFAULT_MAX_MEDIAN_TRAVERSAL_LEN = 1_000` `:220`, `DEFAULT_MAX_TOTAL_SEQUENCE = 1_000_000` `:221`, `DEFAULT_MAX_TRAVERSALS = 10_000` `:222`.
3. **What it does:** for `method=auto`,
   ```
   if  auto_spoa_max_traversal_len > 0
   AND traversal_stats.max_len ≤ auto_spoa_max_traversal_len     // default 2 kb
   AND direct_replacement_within_budget(...)                     // count, max-len, median-len, total-len budgets
   then ResolutionMethod::Poa
   else ResolutionMethod::Allwave
   ```
   For any other `method=...` (poa, poasta, star-biwfa, allwave, sweepga) the user's choice is honoured for every bubble regardless of size — there is no auto-routing.
4. **Divergence from spec:** four distinct departures from the spec table.
   - **Decision variable.** Spec says route by **median** traversal length. Code routes by **max** traversal length (`traversal_stats.max_len`, `:845`). For a bubble whose paths run 100 / 100 / 100 / 10 000 bp the spec routes to sPOA (median 100 bp); the code routes to Allwave (max 10 000 bp). The audit verified this on the slice: round 1 selected n=147 candidates with `max-len median/max = 129/498` and `median-len median/max = 127/406`, all under the 2 kb cutoff, so all 147 went to SPOA — but on a bubble with one long outlier traversal, max and median diverge and the code would mis-route.
   - **Thresholds.** Spec calls for two cutoffs: ~1 kb (sPOA→POASTA) and ~10 kb (POASTA→sweepga). Code uses one cutoff (`auto_spoa_max_traversal_len`, default 2 kb) producing two tiers.
   - **Aligner choice above the cutoff.** Spec calls for **POASTA** at 1–10 kb and **sweepga** above 10 kb. Code defaults to **Allwave/seqwish** above 2 kb (`ResolutionMethod::Allwave` at `:850`), never to POASTA in auto mode and never to sweepga in auto mode. Sweepga and POASTA are reachable only by explicit `--method sweepga` / `--method poasta`, in which case **every** bubble goes through the chosen aligner regardless of size.
   - **Auto path is binary.** Spec is three-tier (sPOA | POASTA | sweepga); code is two-tier (SPOA | Allwave).
5. **Classification:** **partial implementation.** The skeleton (per-candidate dispatch, per-method backend, parallel fan-out) is in place. The actual routing function is wrong both in the decision variable (max vs median) and in the codomain (SPOA / Allwave vs SPOA / POASTA / sweepga). Localized fix at `candidate_selection_method` + thresholds (≤30 LOC), plus a small amount of wiring to make POASTA and sweepga reachable from auto.

### Phase 3 — Local bubble graph construction (and fresh ID allocation)

1. **Present:** yes for the per-bubble extract + aligner step; **no for fresh-ID-at-`n+size`**.
2. **Where:**
   - extract paths inside the bubble: `materialize_candidate_sequences` `src/resolution.rs:1417`, fed by `candidate_named_sequences` `:1616` and `range_sequence` `:807`.
   - pass to aligner: `build_replacement_with_method` `:1580`, which dispatches to `build_poa_replacement` `:2004`, `build_poasta_replacement` `:2033`, `build_biwfa_inmemory_replacement` `:2137`, `build_allwave_seqwish_replacement` `:1833`, or `build_sweepga_seqwish_replacement` `:1927`.
   - aligner returns a local graph: each backend produces a per-replacement `Graph` (`finalize_pairwise_induced_replacement` `:1769` for the pairwise-induced ones).
   - **fresh ID allocation:** the actually-used function is `next_unused_segment_id` `src/resolution.rs:2938`, called from `render_rewritten_graph` `:2887`. There is **no** `next_unused_segment_id`-style allocator at `n + current_graph_size`.
3. **What it does:** for each bubble:
   - Builds `(headers, seqs)` where `headers[i] = "__impg_bubble_path{path_idx}_{i}"` and `seqs[i]` is the sequence of the path's traversal `[begin_step, end_step)` including both boundary anchor segments (`unique_anchor_range` returns `(i, j+1)` where `i, j` are the anchor step indices).
   - Feeds those sequences to the aligner; aligner returns a local GFA which is re-parsed into a `Graph`, then `unchop_gfa` (gfaffix-style) compaction is applied (`:1795`), then the bounded polish runs (`:1803`).
   - At emit time (`render_rewritten_graph` `:2846`), the renderer iterates `out_paths`, collects unique `OutNode`s, and for each `OutNode::Replacement(plan_idx, node_idx)` calls `next_unused_segment_id(&mut used_ids, &mut next_id)` `:2887`. `next_id` starts at `1usize` `:2860` and `used_ids` is seeded only with the **survivor** original IDs `:2856-2859`. The loop in `next_unused_segment_id` `:2938-2946` increments `next_id` until it finds one not in `used_ids`.
4. **Divergence from spec:**
   - **Fresh ID range.** Spec: "Allocate fresh node IDs starting at `n + current_graph_size` so the new local graph's IDs are outside the existing range and cannot collide with surviving segments." Code: starts at `1` and skips survivor IDs by string set lookup. Numerically verified on slice-auto-1round: input GFA has 2 942 segments with integer IDs in `[50064, 272213932]`; output has 1 205 new-replacement IDs all in `[1, 1205]`. All 1 205 are below the original-max ID (`272213932`); the spec's structural invariant is fully broken even though the operational invariant (no collision with survivor IDs) is preserved by the `used_ids` set membership test.
   - **Boundary handling.** Spec: "The boundary nodes themselves are shared with the rest of the graph — they MUST NOT be re-emitted with fresh IDs." Code: `unique_anchor_range` `:1412` returns `(i, j+1)` where `i`, `j` are the indices of the entry and exit anchor steps **inclusive**; `apply_replacement_frontier` `:1503-1508` replaces `path.steps[begin_step..end_step]` (i.e. the anchors too) with replacement steps. The aligner is given sequences that include the anchor sequences, and the aligner's output graph contains its own version of those anchor sequences as replacement nodes with fresh IDs. So boundary anchors **are** re-emitted with fresh IDs, in violation of the spec invariant. This is half of the proximate cause of invariant 1 violations (the same anchor sequence appears once per surrounding plan that touches it).
5. **Classification:** the per-bubble extract + aligner + per-bubble-graph step is **partial implementation**. The fresh-ID-allocation strategy is a **localized bug** (1-line fix to seed `next_id` at `original.segments.len()` or `max(original_int_ids)+1`). The boundary-shared-with-parent rule is a **partial implementation / phase missing** — the code has no concept of "boundary node preserved with original ID"; the path rewrite at `:1503` and the candidate-range definition at `:1412` would both need to change so that the replacement only replaces `[begin_step+1, end_step-1]` (or the equivalent open-interval semantics in the aligner-input view).

### Phase 4 — Strike + link

1. **Present:** yes, implicitly via path rewriting and the re-render round trip.
2. **Where:** `apply_replacement_frontier` `src/resolution.rs:1452-1577`. Path step rewrite at `:1494-1528`; render at `:1551` (`render_rewritten_graph` `:2846`); re-parse at `:1552`; per-path byte-equality validation at `:1572` (`path_sequences_equal` `:2832`).
3. **What it does:**
   - For each plan, builds `replacements_by_path: FxHashMap<usize, Vec<PathReplacement>>` mapping each path index to the per-path replacement intervals `[begin_step, end_step)` and the replacement step list `OutNode::Replacement(plan_idx, step_idx)` (`:1455-1478`).
   - Sorts per-path replacement intervals and rejects any overlapping pair within one path (`:1480-1488`).
   - Per-path linear scan: walk the path's `steps`; when `i == replacement.begin_step`, emit the replacement's steps verbatim and jump `i` to `replacement.end_step`; otherwise emit the original step as `OutNode::Original(node)` and tag the surviving node id in `used_original`.
   - `render_rewritten_graph` walks the rewritten paths and emits **only** the nodes (`OutNode`) seen on at least one path. Original nodes that no longer appear on any path are dropped from the output GFA simply by not being emitted; their old IDs end up in neither `used_ids` nor any S-line.
   - L-lines are reconstructed from consecutive path-step pairs (`:2898-2917`) so cross-bubble edges to/from the surviving boundary path-step are re-issued.
4. **Divergence from spec:**
   - **Strike works**, but indirectly: old inner segments are dropped because the renderer's output graph is path-defined, not because there is an explicit removal step. (You see this in slice-auto-1round: input 2 942 segments, output 2 984 segments; 2 942 − 1 779 = 1 163 original segments were dropped, replaced by 1 205 fresh-ID nodes.) The strike is **not** a localized bug — it is operationally correct.
   - **Boundary preservation is missing** (see Phase 3): "old inner nodes are GONE" per the spec, but the code's "inner" range also includes the boundary anchors, so boundary nodes are stricken alongside actual inner nodes, then re-emitted as fresh-ID replacements.
   - **No seq → id deduplication across plans.** Two replacements that emit byte-identical sequences (typical when a bubble's flanking anchor sequence appears in another sibling bubble's replacement) each get a fresh `next_id`, producing two distinct S-lines. The per-replacement parse/render at `:1552-1571` carries the dup-seq pathology through the parse/validate round trip (the duplicate-segment debug log at `:1563-1570` fires before validation, then validation passes because path sequences are still preserved byte-for-byte).
   - The spec's phrasing implies an explicit "remove inner nodes from parent graph" step; the code's "drop by not being emitted" approach is operationally equivalent for the strike but loses any chance to assert that the strike actually happened (no `assert!(original_id not in next.segments)` survives the re-render).
5. **Classification:** **partial implementation.** The strike-and-link semantics are present in spirit; what is broken is (a) boundary inclusion in the strike range (a Phase 3 / Phase 4 boundary issue) and (b) the absence of seq→id deduplication at emit time. Both are localized fixes inside `apply_replacement_frontier` + `render_rewritten_graph` (~30-50 LOC total). The `crush-trace.md` § 4 item 2 ("`OutNode::Replacement(plan_idx, node_idx)` carried all the way to the serializer without being projected through `seq → id`") names exactly this gap.

### Phase 5 — Batched lacing

1. **Present:** yes.
2. **Where:** the per-bubble build is rayon-parallel at `src/resolution.rs:515-563` (`par_iter().map(|candidate| build_replacement_with_method(...))`). The strike+link applies all accepted plans as a single batch at `src/resolution.rs:611-617` (`apply_replacement_frontier(&graph, &plans)`).
3. **What it does:** each round (a) computes the candidate frontier (POVU + non-overlap selection), (b) builds all replacements in parallel, (c) collects accepted `ReplacementPlan`s, (d) calls `apply_replacement_frontier` once with the whole `Vec<ReplacementPlan>`. `apply_replacement_frontier` per-path linearly rewrites each path through all of its replacements at once and then re-renders the whole graph once. There is no per-replacement sequential re-render.
4. **Divergence from spec:** the spec also says "Sort the resulting graph". The code does **not** sort the rewritten graph after a round; the final `render_graph` at `src/resolution.rs:636` emits paths and segments in path-discovery order (via the `ordered_nodes` walk at `render_rewritten_graph:2862-2870`, which preserves first-path-step-of-occurrence ordering). The downstream `:nosort` pipeline modifier (visible in `docs/c4-crush-handoff.md`) intentionally skips the optimization sort. Whether this divergence matters depends on downstream consumers; for invariants 1–4 it has no effect.
5. **Classification:** **present, mostly correct.** The non-overlap precondition the spec relies on is enforced both at selection time (`occupied_by_path` in `find_candidate_frontier`) and at apply time (the per-path overlap rejection at `:1482-1487`). Sort-after-batch is not present but is not a correctness issue.

### Phase 6 — Level descent (iterate)

1. **Present:** **no — what is present is "re-run POVU on the rewritten graph", not "descend into the resolved bubble".**
2. **Where:** the round loop `src/resolution.rs:454-632`. Each iteration calls `find_candidate_frontier(&graph, config, &seen)` on the **current working graph** (which is the post-rewrite graph from the previous round). The `seen: FxHashSet<String>` (`:436`) accumulates `candidate_signature` strings (`:1430`) so the same signature is not re-selected.
3. **What it does:** every round = a fresh POVU pass over the entire working graph, followed by candidate selection, build, and apply. The `seen` set's signature is `"{root_path_name}:{root_pos_begin}-{root_pos_end}:path_idx:begin-end|..."` (`:1430-1450`) — coordinates in the **previous-iteration** working graph. After the rewrite, those step coordinates shift, and a "site at the same biological position" reappears with a different signature, so it can be re-selected.
4. **Divergence from spec:** the spec says "Inside each resolved bubble there may be smaller bubbles → resolve at level-2. **Do not re-resolve the same level-1 bubble — that's wasted work**." The code does not descend into a specific resolved bubble; it runs the full POVU decomposition on the whole graph again, relying on `seen` for dedup. Two consequences:
   - When the previous round's resolution **widened** a bubble's path-coordinate boundaries (Failure 3 / Failure 4 in [`docs/crush-audit.md`](crush-audit.md), which is exactly this), the same biological site re-enters the frontier as a much bigger candidate (round-2 candidate 1 on real C4: traversals=433, median-len=42 430, root-span=42 394 — vs round-1 median-len 157).
   - There is no notion of "inside the resolved bubble"; descent and re-resolution are conflated. The level information POVU computes (`FlubbleSite.level`, `FlubbleSite.parent_id` in `povu-rs/src/native_gfa.rs:459-468`) is fetched but never used by `find_candidate_frontier`.
5. **Classification:** **partial implementation / wrong shape.** Iteration exists; descent semantics do not. The `seen` set protects against trivial re-resolution but not against the spec's failure mode (re-resolving a bubble whose boundary widened). Localized-but-non-trivial fix: filter `decomposition.sites` to `level == 0` (or use the parent_id to walk a tree) and only descend after a successful resolution; or alternatively keep the current "re-decompose" model but stop the loop when the per-round monotone-progress invariant fails (which is `crush-audit.md` Failure 2 territory). 50-150 LOC depending on which model is chosen.

### Phase 7 — POA normalization of small recurrent motifs

1. **Present:** yes, but only **inside the per-replacement pipeline** — not as a separate pass over the output graph.
2. **Where:** `polish_replacement_gfa` `src/resolution.rs:2349`, invoked at `:1803` from `finalize_pairwise_induced_replacement` `:1769`. Two backends:
   - `polish_replacement_gfa_with_flubbles` `:2365` — re-decomposes the per-replacement GFA with POVU, runs SPOA or POASTA on each inner flubble recursively (it builds a fresh `ResolutionConfig` and calls `resolve_graph_bubbles` on the per-replacement graph at `:2396`).
   - `polish_replacement_gfa_with_smooth` `:2399` — runs smoothxg-style block smoothing.
3. **What it does:** after AllWave/seqwish or sweepga/seqwish builds a per-replacement graph and `unchop_gfa` compacts it, the polish pass either (a) recursively crushes the small inner flubbles with SPOA (`polish_method=poa` — the CLI default at `src/main.rs:4754`), or (b) runs smoothxg on the per-replacement graph (`polish_method=smooth` — what the canonical C4 command in `docs/c4-crush-handoff.md` uses). The polished replacement then goes through `validate_replacement_paths` before being accepted.
4. **Divergence from spec:** the spec phrasing is "After sweepga (or allwave) on a bubble, the local graph may contain tight high-copy nodes from short tandem repeats. These should be **re-aligned locally with sPOA / POASTA** to linearize." The code's `polish_method=poa` matches the spec exactly. `polish_method=smooth` (the canonical-command choice) **does not** match the spec — smooth is a different normalization algorithm — but it is still a "local re-alignment after pairwise induction" step at the right place in the pipeline.
5. **Classification:** **present, mostly correct.** Under default settings (`polish_method=poa`) the implementation matches the spec. Under the canonical command's settings (`polish_method=smooth`) the implementation uses smoothxg in place of SPOA, which is a behavioural rather than architectural divergence. The trace's Failure 1 wall-clock blow-up lives in `polish_replacement_gfa_with_smooth` `:2399` when smoothxg is fed an 18 MB / 440-path bubble — that is a Phase 2 routing problem (the bubble should not have been routed to sweepga + smooth at all), not a Phase 7 problem.

### Phase 8 — (Optional) final compaction

1. **Present:** yes per-replacement (`unchop_gfa` at `:1795`); **no** as a separate final-emit pass on the output of `resolve_graph_bubbles`.
2. **Where:** `unchop_gfa` is called inside `finalize_pairwise_induced_replacement` `:1795`. The final emit at `src/resolution.rs:636` is `render_graph(&graph)` with no further sequence-dedup.
3. **What it does:** every per-replacement build runs gfaffix-style compaction (`unchop_gfa`) on its own replacement GFA before the polish pass and before the final replacement integration. No equivalent sweep runs on the post-`apply_replacement_frontier` working graph or on the final emit.
4. **Divergence from spec:** the spec says "Run a gfaffix-style sequence-dedup pass on the final emit. If Phases 1–7 are correct, this should be a no-op — duplicate-sequence segments shouldn't exist by construction." The code never runs the dedup pass at the final-emit boundary; instead the `:nosort` output stage explicitly **skips** any downstream gfaffix. Because Phases 3 and 4 leak duplicate-sequence segments across plans (see Invariant 1 below), the missing final pass means duplicate-sequence segments make it into the on-disk output.
5. **Classification:** **phase entirely missing at the final-emit boundary**. The per-replacement `unchop_gfa` does not substitute for it because compaction is local to one replacement, not cross-replacement. A localized add: call `unchop_gfa` (or seq-keyed dedup) once at `:636` before `render_graph`. ≤30 LOC, but only useful as a band-aid — the spec is explicit that "if Phases 1–7 are correct, this should be a no-op." A working final-emit gfaffix pass would mask Failure 3 from `docs/crush-audit.md`, not fix it.

## Invariants — graph-level check on real C4 slice

Numbers below come from the three crush runs listed at the top, computed by a
small Python pass over the emitted GFAs that counts unique segment sequences,
duplicate-sequence segments, and the integer-id range of original survivors vs
new replacements. The full source of the check is the inline script that produced
the table; results were copied verbatim.

| invariant | slice-auto-1round | slice-sweepga-1round | slice-auto-3round | holds? |
|---|---:|---:|---:|---|
| 1. no two segments share a sequence (dup-extras = 0) | **948 / 2 984** | **2 056 / 4 289** | **1 160 / 2 773** | **FAIL** (32% / 48% / 42% of output segments are duplicate-sequence) |
| 2. every input path is preserved exactly | enforced by `path_sequences_equal` `:1572` — runs accepted | runs accepted | runs accepted | **PASS** (the resolved-gfa step would return `io::Error` otherwise; all three runs exited 0) |
| 3. each round shrinks the graph (segments ≤ pre-round AND segment-bp ≤ pre-round) | seg 2 942 → 2 984 (**+1.4%**), bp 64 348 → 58 003 (−9.9%) — segments grew; bp shrank | seg 2 942 → 4 289 (**+45.8%**), bp 64 348 → 178 933 (**+178%**) — both grew | round-by-round: r1 +1.4%/−9.9%; r2 −2.4%/−10.4%; r3 −4.8%/−8.6% — r1 segments grew, r2-r3 OK | **FAIL** on round 1 segment count for all three runs; **FAIL** completely under sweepga; matches the `docs/crush-audit.md` real-C4 result (+79% bp on round 1 of full C4) |
| 4a. no segment ID re-used between original survivors and new replacements | 0 collisions (1 205 new ids vs 1 163 dropped originals; `used_ids` set check at `:2856-2859` enforces this) | 0 collisions | 0 collisions | **PASS** |
| 4b. new IDs in `[n + current_graph_size, ...)` (n = original segment count = 2 942; original ID max = 272 213 932) | all 1 205 new ids in `[1, 1205]` — all 1 205 are below original-max (272 213 932) | all 1 385 new ids in `[1, 1385]` | all 1 734 new ids in `[1, 1734]` | **FAIL** — the structural property the spec names is fully broken on every run; the spec's operational consequence (no collision) is preserved by the `used_ids` membership test |

Invariant 3 is the headline failure — and it's the same one the existing
[`docs/crush-audit.md`](crush-audit.md) reports on full C4 with sweepga (+78.8%
bp on round 1, +88.9% on ws-max). On the slice with `method=auto` (so every
bubble below 2 kb max-len is routed to SPOA), invariant 3 fails on segment
count but passes on bp; this confirms that the per-round growth is **not**
purely a sweepga-routing artefact — the per-plan emit pathway (Phase 3/4
duplicate-sequence emission) inflates segment count even in pure-SPOA mode.

Invariant 1 dovetails with invariant 3: 948 / 2 984 = 31.8% duplicate
sequences on the slice even under method=auto means that ~1/3 of the round-1
emit is sequence-redundant — every duplicate segment is, by construction, a
candidate cross-plan dedup target. The 178% bp growth on slice-sweepga-1round
is dominated by the same effect amplified by sweepga's coarser pairwise output
(2 056 duplicate-extras on 4 289 segments = 48% of the output is sequence
duplicate).

Invariant 4a is the only invariant that fully holds on every run, and it
holds **only because** the `used_ids` membership test at `next_unused_segment_id:2942`
explicitly prevents collisions. The structural property the spec actually
names (4b: new ids start at `n + current_graph_size`) is fully broken on
every run — a 1-line fix at `render_rewritten_graph:2860` (seed
`next_id = original.segments.len()` or `next_id = max_original_int_id + 1`)
would make this hold without changing any other behaviour.

## Scope estimate

Per-phase classification with a rough LOC estimate. "LOC" is changed lines in
`src/resolution.rs` unless noted, ignoring tests.

| phase | divergence | class | rough LOC |
|---|---|---|---:|
| Phase 1 | greedy non-overlap selection is not top-level; level/parent fields unused | (a) localized fix | ~20 LOC in `find_candidate_frontier` to filter by `site.level == 0` (or walk `parent_id`) before the non-overlap selection |
| Phase 2 | routing uses `max_len` not `median_len`; only 2 tiers, not 3; never selects POASTA or sweepga in `auto` | (a) localized fix | ~30 LOC in `candidate_selection_method` + new thresholds + tests; CLI flag wiring |
| Phase 3 — extract+aligner | path extraction includes boundary anchors | (b) partial-impl-needs-completion | ~30-60 LOC across `unique_anchor_range`, `materialize_candidate_sequences`, `apply_replacement_frontier`. Requires defining boundary semantics with POVU site-definition (spec open question) |
| Phase 3 — fresh ID range | `next_id = 1` instead of `n + size`; structural invariant 4b broken | (a) localized fix | 1 LOC at `render_rewritten_graph:2860` (`let mut next_id = original.segments.len() + 1;`) |
| Phase 4 — strike | works indirectly via path-rewrite + re-render; no explicit assertion | (a) localized fix (only adds an assertion or explicit removal pass) | ~10 LOC if we want explicit, 0 LOC if we keep implicit |
| Phase 4 — cross-plan seq dedup | duplicate-sequence segments emitted (Failure 3 in [`docs/crush-audit.md`](crush-audit.md)) | (b) partial-impl-needs-completion | ~40-80 LOC in `render_rewritten_graph` to add a `BTreeMap<Vec<u8>, String>` from sequence to canonical id and route every `OutNode` through it. Path L-line reconstruction needs adjustment for orientation |
| Phase 5 — batched lacing | sort-after-batch missing | (a) localized fix | ~5-15 LOC if we want to add sort; otherwise (0 LOC) leave as-is |
| Phase 6 — level descent | re-runs POVU on the rewritten graph instead of descending into resolved bubbles | (c) phase-missing-needs-implementation | ~80-150 LOC: track resolved-bubble path-ranges, restrict next-round POVU to those ranges (or refactor the loop to iterate over `FlubbleSite.parent_id`/`level`) |
| Phase 7 — POA post-normalization | present for `polish_method=poa`; canonical command uses `polish_method=smooth` instead | (b) partial-impl-needs-completion (decision-level) | ~5-10 LOC to change the canonical command default; alignment with spec is mostly a policy/CLI choice |
| Phase 8 — final compaction | no final-emit dedup pass; depended on Phases 1-7 being correct, which they aren't | (c) phase-missing-needs-implementation **only as a band-aid** | ~20-30 LOC to add `unchop_gfa` (or seq-keyed dedup) at `resolve_graph_bubbles` `:636`; spec says this should be a no-op once upstream is fixed, so the band-aid is acceptable but should not substitute for the Phase 3/4 fix |

**Rough total**, treating Phase 8 as a band-aid that should be deferred until
Phases 3/4 are correct: **~200-400 LOC** across `src/resolution.rs` plus a
handful of new tests, concentrated in `find_candidate_frontier`,
`candidate_selection_method`, `apply_replacement_frontier`, `render_rewritten_graph`,
and the round loop in `resolve_graph_bubbles`. The fresh-ID-range and
routing-by-median fixes are 1-30 LOC each and can land independently of the
larger Phase 4 / Phase 6 refactors.

## Notes on what was explicitly *not* changed by this audit

- **No new instrumentation.** The two `log::debug!` lines inside
  `apply_replacement_frontier` at `src/resolution.rs:1543-1549` and
  `:1563-1570` were committed by `wg/agent-55/crush-trace` and remain
  unchanged. The invariant-1 dup-segment count surfaced by the second of
  those lines (`crush apply: ... U unique segment sequence(s), D duplicate-sequence segment(s)`)
  was verified to be the same as the independent Python-pass count over the
  emitted GFA — 948 in both cases — so the instrumentation is correct and
  load-bearing for this audit.
- **No code change.** This audit is observation-only per the task description.
  The fixes named above are not proposed here — the downstream
  [`flip-crush-spec-audit`](#) task is expected to take this divergence list as
  input for a redesign.
- **Slice-scale only.** Per the task brief, invariant checks were run on the
  committed C4 slice (`tests/test_data/crush/c4_slice_1500_3000.gfa`), not on
  the full GRCh38 C4 input. The slice reproduces the same architectural
  pathologies (invariant 1 fails at 32-48% dup rate; invariant 3 fails on
  segment count; invariant 4b fails fully) without needing a 30-minute wall
  budget. The full-C4 numbers in [`docs/crush-audit.md`](crush-audit.md) are
  consistent with the slice numbers reported here.
