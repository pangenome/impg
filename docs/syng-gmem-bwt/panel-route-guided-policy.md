# Guided panel-route policy v3 (experimental, no sequence emission)

`genome-infer search-panel-routes-guided` is a **separate command-owned policy**.
The existing `search-panel-routes` lexical DFS, graph format, compiler identity,
loader seals, native counting, complete-assignment evaluator, normalization and
lower bound are unchanged. Existing graphs are accepted only through normal
`Graph::load` and `Evaluator::new` verification, never migrated or re-identified.

## Parent release validation

V3 passed **682 ordinary Rust tests**, **nine additional configured scientific,
update and driver gates**, a repeated every-primitive resume gate, and **six
Python tests**. Parent checks verified the source package, patch apply/reversal,
file modes and all 96 unchanged protected backend files/gitlinks. Independent
source review found no issues. This establishes tested enumeration and resume
behavior, not real-panel genome recovery, QV or coverage.

## Source pair order and explicit v2 update

V3 changes **only source-port enumeration**: for a virtual `2^b` domain,
`index = 2 * reverse_(b-1)(ordinal / 2) + ordinal % 2` (0/1 handled directly).
Adjacent source positions are interleaved while anchor-pair groups remain spread.
Plain bit reversal on paired forward/reverse source records visits only even
positions for half the domain; the paired order offers both orientations in its
first two positions. This is ordering, not candidate/evidence filtering or an
assumption that adjacent records are homologous copies. Arbitrary bases and odd
lengths retain full coverage. Hub-member order remains plain bit reversal.

Ordinary v3 resume requires an exact v3 policy identity. To preserve an existing
compatible `fair-sample-ranked-live-continuation-v2` checkpoint, explicitly run:

```
impg genome-infer search-panel-routes-guided \
  --panel PANEL --routes ROUTES --sample SAMPLE \
  --haploid-depth DEPTH --background 0.1 \
  --resume-from V2_OUTPUT --update-source-pair-order \
  --max-work OLD_WORK --max-evaluations OLD_EVALUATIONS \
  --max-state-bytes AUTHORIZED_BYTES --max-optima OLD_OPTIMA \
  --out-dir NEW_UPDATE_OUTPUT
```

This is **update-only**, never search or native rescoring. Preserve all prior paths.
Use `--extend-budgets` for any cap increase, including space for new cursor fields.
Then resume `NEW_UPDATE_OUTPUT` normally with explicitly increased work if desired.
A repeated update of v3 is rejected; updating the same immutable v2 parent twice
produces equivalent states. Actual old/new source identities are recorded honestly
in `order_update`; all non-policy scientific/backend bindings remain exact.

Schema3 source scans retain `{len, ordinal, old_prefix}`. The update preserves each
legacy scan's visited ordinal boundary in `old_prefix`, sets its paired cursor to
zero, and changes no other task payload, live link, FIFO/shallow index, focus,
clock, native objective bit, counter or ledger byte. Task storage is recomputed;
the old occupancy is recorded and the peak never decreases. Each paired cursor
position is charged normally, including virtual positions and indices whose
inverse old ordinal is below `old_prefix`. Thus the old visited set and newly
yielded set are disjoint and together cover the original domain. Pending bounds
create fresh paired scans (`old_prefix=0`). No native work is restarted.

V2 is accepted only for explicit updates and historical ancestry. Existing
v1-to-v2 receipts remain readable using their recorded historical source identity
and ordinary conservation checks, not comparison to today's transform source.
Direct v1 conversion is rejected in v3: use the preserved v2 executable for the
historical conversion below, then update its v2 output. No new source archive,
signature, or authentication machinery is introduced.

Synthetic opt-in gates in `tests/panel_route_source_pair_order/mod.rs` require
fresh `IMPG_TEST_PAIR_REVERSE_OUTPUT` / `IMPG_TEST_PAIR_UPDATE_OUTPUT`; the latter
also requires preserved `IMPG_TEST_GUIDED_V1` and `IMPG_TEST_GUIDED_V2` binaries.
The reverse recipe is seed113, 400bp, mutations150/260 in B followed by reverse
complement, repeated1000, with fixed 200,000 work / 10,000 evaluations / 256MiB
logical-state caps. It requires actual complete mixed assignments with positive
reverse donor spans, not merely a reverse task. The update recipe repeats100 and
checks caps1/2/5/32/128/512/2048/8192, nonzero focused reverse work, exact state and
ledger conservation, repeated updates from equivalent parents, resumed continuation,
and genuine v1→v2→v3 ancestry. Existing forward/coupled/finite gates remain in place.
Synthetic success does **not** establish real best-family mixed recovery, whole-genome
reconstruction, QV, coverage, or sequence-emission authorization.

## Historical v2 parent validation

The corrected v2 passed **678 ordinary Rust tests** (16 explicitly ignored),
**eight configured oracle/driver/long-producer/transition tests**, **six Python
tests**, and the near-cap backlog continuation gate. Parent source checks confirm
all 96 protected backend files/modes and gitlinks are unchanged. Independent
review's remaining history-consistency finding was reproduced, corrected and
re-reviewed clean: conversion retains its prior ledger, and each ordinary resume
preserves the exact parent prefix while allowing new records to be appended.
These are normal scientific-data consistency checks, not cryptographic assurances.

The prior real-panel v1 search failed to construct mixed assignments; successful
synthetic v2 tests do not establish real-panel recovery or authorize sequence
emission. The next milestone is useful continued search on that preserved problem.

## CLI and explicit resume

```
impg genome-infer search-panel-routes-guided \
  --panel PANEL --routes ROUTES --sample SAMPLE \
  --haploid-depth DEPTH --background 0.1 \
  --max-work 100000 --max-evaluations 10000 \
  --max-state-bytes 268435456 --max-optima 1000 --out-dir NEW_OUTPUT

impg genome-infer search-panel-routes-guided \
  --panel PANEL --routes ROUTES --sample SAMPLE \
  --haploid-depth DEPTH --background 0.1 \
  --resume-from PREVIOUS_OUTPUT --extend-budgets \
  --max-work 200000 --max-evaluations 20000 \
  --max-state-bytes 268435456 --max-optima 1000 --out-dir ANOTHER_NEW_OUTPUT
```

All four resource caps are **cumulative**, cannot decrease, and increases require
`--extend-budgets`. Specify the old values explicitly when they differ from CLI
defaults. `--tie-epsilon` (default `1e-9`), depth, background, feature cap and cache
cap are immutable bindings. Fresh output is atomically reserved; even an existing
empty directory fails. A budget stop is a successful, explicitly incomplete run,
not exhaustion. Reaching the evaluation cap stops at the next *fresh* evaluation,
not at recorded native-score reuse or a port operation. Initialization may stop
partway; no native assignment is rescored on resume. Loader/sample verification
and evaluator cache reconstruction still run normally.

## Scheduler and construction

* Every family's complete native assignment is scored once. Initialization cursor,
  bit-exact scores and root tasks persist. Turns begin only after initialization;
  ascending native objective determines family order, with identity ties.
* All live families rotate in **32 primitive-work** quanta. Within each family the
  service cycle is FIFO, shallowest committed switch depth (creation-ID ties),
  live focused continuation. A task has one owner, two exact ordered indices, and
  optional reciprocal previous/next links in its family's single live return chain.
  An empty focus seeds from FIFO on the third share. Unfinished
  tasks move to the FIFO tail; modes, IDs, requeue clocks and remaining quantum
  persist. FIFO has a fixed one-third service share. No stale index entries exist.
* Each source-bound probe, source ordinal, hub-bound probe, hub-member ordinal,
  pairwise feasibility comparison and traversal transition is one resumable
  primitive. Scanner and hub tasks are independent. Source-paired and hub bit-reversed
  ordinal orders spread visits over indexed intervals; virtual out-of-range positions and
  failed proposals also consume work. A `u128` cursor safely represents a virtual
  `2^64` domain. Record reads always seek to a checked absolute offset, including
  when interleaved with evaluator operations on shared handles.
* Only completed routes and committed positive prefixes reserve canonical spans.
  Both orientations share the same source-forward capacity. Legal forward closure
  at the paired native endpoint creates ordinary next-slot work **and** a lazy
  native-tail completion probe. A probe conflict destroys only the probe. No
  hypothetical native tails constrain ordinary coupled exchanges.
* Children are immediately runnable, without switch-depth limits, shallow-stream
  exhaustion, or a successful single-donor return bridge. Complete candidates use
  the unchanged whole-genome evaluator. No partial likelihood, prior, donor score,
  frequency mask, annealing or altered physical-copy/count model is introduced.
* Exact backend-normalized native completions reuse their recorded complete score.
  Visits and reuse are separate from fresh evaluations. Non-native memoization is
  deliberately absent; non-native duplicates may consume evaluations, but neither
  physical rates nor correlated optima are duplicated.

Work is a reproducible **search operation count**, not CPU time. Native evaluation,
complete validation/counting/scoring, graph/sample/source verification, cache work,
serialization, ordered-index maintenance and I/O are not bounded constant-time
operations. Fresh evaluator calls are atomic and separately counted. There is no
claim that 32 work units take equal wall time across families.

## Durable artifacts and provenance

`result.json` retains the existing backend/graph/sample/depth/background envelope;
its result contains the new policy identity, ranking, native scores, cumulative
statistics, per-family work/modes/mixed completions, switch distribution, visited
source IDs and actual mixed-route donor source IDs (excluding each paired endpoint
source), completion-probe outcomes, occupancy, incumbent and correlated support.
`incumbent-assignment.json` contains coordinates only, not spelled sequence. Guided
search does not make a hidden extra evaluator call to write incumbent diagnostics.

`evaluations.jsonl` is the cumulative semantic ledger: visit, fresh-call count,
work, stable task ID, initialization/fresh/reuse kind, normalized assignment,
assignment checksum, objective and switch/mixing flags. It excludes timings and
cold-cache diagnostics. Objective state is serialized as IEEE-754 **integer bits**,
so JSON floating-point parsing cannot perturb ranking or resumed tie decisions.

`checkpoint.json` schema version 3 contains (v2 ancestry is retained):

* `bindings`: actual separately named source checksums for policy, state machine,
  adapter, checkpoint and CLI, policy version, normally computed backend compiler
  identity, graph and sample checksums, count policy, scientific parameter bits,
  feature and cache caps;
* `parent`: canonical immutable parent output path and expected checkpoint/ledger
  identities (or null); `budgets`: cumulative, explicitly extended cap history;
* `state`: initialization/ranking, tasks and all operation phases (lookup bounds,
  permutations, committed segments, pending children, closures, probes and exact
  evaluations), both scheduler indices, family heads, live links and counters, native-score table,
  incumbent/support, sticky support loss and logical occupancy accounting.

`checkpoint-seal.json` binds the exact checkpoint bytes and cumulative ledger by
length and FNV-1a64, matching the repository's corruption-checksum convention.
These are **not cryptographic authenticity signatures**. Resume verifies both
seals, bindings and every immutable ancestor's expected checkpoint/ledger identity;
it copies the ledger prefix to the new directory and verifies the copied bytes
before proceeding. Ledger and checkpoint files are synced; checkpoint installation
streams through a command-owned buffered atomic JSON writer (byte-identical pretty
JSON, no trailing newline), and the directory is synced on Unix. Serialization
failure leaves an `.incomplete` file, never an installed successful checkpoint. Preserve
all ancestors and their original paths. Missing seals/ancestors and changed bindings
fail; resuming a partially written checkpoint is not supported. File permissions
are not an immutable-storage enforcement mechanism: callers must keep parents
unchanged. Failed/running output manifests never mean successful completion.

## Logical storage contract (not RSS)

`--max-state-bytes` caps a deterministic retained-state occupancy measure:

* compact serialized bytes of each owned task (including context, pending operation,
  segment bindings and permutations), plus **256 units per task** for map ownership
  two ordered indices and live-link allowances;
* compact serialized bytes of every native score, incumbent and retained correlated
  score, plus **64 units per score** (each actual stored copy is charged);
* **4096 fixed units**, **512 per family** for scheduler/statistic storage, and
  **64 per occupied switch-histogram entry, visited-source ID or mixed-donor ID**;
* compact serialized bytes of bindings, cumulative budget history and parent link;
  conversion-boundary receipt metadata and historical accounting are also charged.

One unit is one logical byte/allowance unit, **not a measured resident byte**.
This does not cap allocator overhead, verified graph/sample/evaluator memory or
caches, transient operation/checkpoint/report buffers, or operating-system caches.
Those unchanged evaluator limits remain separate. Peak task count and core-state
high-water mark persist; reported total occupancy also charges current checkpoint
metadata (which naturally differs between resumed and uninterrupted invocations).

Before consuming work, a conservative growth reservation is checked:
`8 * (selected_task_weight + native_assignment_JSON_bytes + 16*k) + 8192`.
Initialization has zero selected-task weight. This bounds retained task/score
copies, a native tail, new indices and counters for one operation; it intentionally
may stop well below the cap. Insufficient reservation leaves the operation and
fairness counters untouched. An insufficient cap for already restored/base state
plus checkpoint metadata fails. No unexplored task is evicted, truncated or spilled.
Support alone may be capped by `--max-optima`, with **sticky** incompleteness even
after a later budget extension or better incumbent.

## Guarantees and gates

Generation-rule completeness, initialization completion, exact evaluation,
construction exhaustion, global optimum certification, correlated support and
biological completeness are separate. Certification requires all initialization
and every scanner, lookup, child, closure, probe and evaluation to drain. Exhaustion
certifies the finite retained domain, not biological topology. The backend lower
bound is unchanged. Biological completeness and sequence-emission authorization
remain **false**, even for exhausted fixtures.

Historical v1 parent release acceptance passed **664 ordinary Rust tests** across 25 suites
(12 explicitly ignored), followed by **six configured oracle/driver tests** and
**six Python tests**. All 96 protected backend source hashes/modes and gitlinks
were unchanged. Independent review found one zero-operation peak-occupancy
reporting defect; it was reproduced, fixed and re-reviewed with no remaining
findings. A resource stop before the first operation now includes its retained
base state in the peak, and extending that checkpoint preserves uninterrupted
ledger/core-state parity. These are source/regression gates, not real-panel
mosaic-recovery acceptance.

Portable gates (no private binaries, external hash tools or real inputs):

```
cargo test --offline --locked --release --lib panel_route_search_policy -- --test-threads=1
cargo test --offline --locked --release \
  --test test_panel_route_policy --test test_panel_route_policy_coupling -- --test-threads=1
```

New tests retain the existing tests unchanged and independently reuse their small
fixture recipe. They check verified adapter parity, invalid/truncated records,
shared-seek interleaving, permutation boundaries, native reuse, all 513 finite
assignment/count/objective/support results against the independent finite engine,
64-by-64 simultaneous copy exchanges with 64 correlated optima, a three-donor chain
with no simple donor-return bridge, and multi-checkpoint semantic ledger/state
parity at initialization, source lookup, hub lookup/scan, child, closure, evaluation,
probe and resource boundaries. Corruption, binding changes, ancestor changes,
output collisions, implicit extensions and sticky support loss are exercised.
Fresh CLI processes reject native `++match`/`++path` dumps and enforce a 1 MiB small
fixture diagnostic cap, without truncating any biological output.

Predeclared useful-work gate: **20,000 work, 10,000 fresh evaluations, 256 MiB logical
state** on paired synthetic samples. It must produce actual mixed complete
assignments in at least two endpoint families and deterministically reverse native
service ranking, not merely change the incumbent. Observed: A/B ranking reverses,
263 fresh mixed completions, 2/4-switch assignments, and fewer than 500 live tasks
at peak (under 0.4 MiB measured logical state). Defaults remain engineering caps,
not validated production capacity guarantees. No default was tuned against real
truth. Large-repeat frontier growth, duplicated non-native completions, atomic
verification/scoring cost and slow endpoint closure remain limitations.

Explicit opt-in Linux gates (requested missing configuration fails):

```
IMPG_TEST_748_ROUTES=/absolute/frozen-748-binary \
IMPG_TEST_POLICY_OUTPUT=/absolute/NEW-fixture-directory \
cargo test --offline --locked --release --test test_panel_route_policy \
  guided_frozen_748_backend_reuse -- --ignored --exact --test-threads=1

IMPG_TEST_F7_FINITE=/absolute/frozen-f7-binary \
cargo test --offline --locked --release --test test_panel_route_policy \
  guided_frozen_f7_finite_oracle -- --ignored --exact --test-threads=1

IMPG_TEST_F7_FINITE=/absolute/frozen-f7-binary \
cargo test --offline --locked --release --test test_panel_route_policy_coupling \
  guided_frozen_coupled_f7_oracle -- --ignored --exact --test-threads=1
```

The frozen-748 producer is SHA-checked before it builds a new synthetic graph;
the candidate CLI loads it normally. Frozen f7 is independently SHA-checked by the
existing configured helper. Optional absolute `IMPG_TEST_POLICY_OUTPUT` and
`IMPG_TEST_POLICY_COUPLING_OUTPUT` directories preserve portable synthetic fixtures
and checkpoints; they must be fresh. Parent-only backend SHA/mode/gitlink audit and
real-graph useful-mixing/resource gates remain required before production acceptance.
No real-data or biological recovery claim follows from these tests.


## V2 continuation and exact-v1 transition

Source scanners and hub scanners hand an active focus to their emitted child,
retaining themselves as live returns. A waiting producer's emission is ordinary
work and does not preempt focus. Waiting nonproducer transformations replace
only their current chain position. Every share performs the same bounded local
link repair; consumed middle/tail tasks are unlinked immediately, not tombstoned.
There is no fourth index, depth cap, stale-entry scan, or discarded construction.

Unchanged native `Start` follows source-bound construction rather than closure.
Changed geometry (including same-source repeats, orientation changes, and earlier
completed mixed routes) follows legal closure feasibility first, retaining its
source-bound alternative behind it. Changed closure follows its completion probe
with ordinary next-slot construction immediately behind. Both probe success and
conflict return to that ordinary work; coupled exchanges are never pruned by the
hypothetical native tail. Native-only completion receives one service when seeded
but is not pinned through the completion chain. All checks, lookups and evaluator
calls retain their original one-primitive charging and native-score reuse rules.

The following is historical v2 behavior: invoke the **preserved v2 binary**, not
v3. Ordinary v2 resume rejects v1. Conversion is a separate **conversion-only** opt-in:

```
impg genome-infer search-panel-routes-guided \
  --panel PANEL --routes ROUTES --sample SAMPLE \
  --haploid-depth DEPTH --background 0.1 \
  --resume-from EXACT_V1_OUTPUT --convert-exact-v1-to-v2 \
  --v1-ancestry-manifest OWNER_FROZEN_PINS.json \
  --max-work OLD_WORK --max-evaluations OLD_EVALUATIONS \
  --max-state-bytes AUTHORIZED_BYTES --max-optima OLD_OPTIMA \
  --out-dir NEW_CONVERSION_OUTPUT
```

Any cap increase additionally requires `--extend-budgets`. Conversion performs
normal graph/sample/source/evaluator construction but **no search or rescoring**.
It verifies exact supported a520 v1 source identity from the five byte-identical
archived inputs under `supported_v1/` (including the entire CLI), not a trusted
hash string supplied by the checkpoint. Those archived inputs, the transition,
schema and validation modules all enter actual v2 policy identity independently.
The backend identity continues to be computed and checked normally.

The separately frozen owner declaration is mandatory and must precede any
suspected corruption. Its schema is:

```json
{"version":1,"checkpoints":[{"directory":"/canonical/old/output",
 "checkpoint":"fnv1a64:BYTE_LENGTH:16_HEX_DIGITS",
 "ledger":"fnv1a64:BYTE_LENGTH:16_HEX_DIGITS"}]}
```

Include the selected root and **every** old ancestor, with no duplicates or
extras. These pins must be independently verified/frozen by the owner, not
created from suspect inputs during conversion. Preserve this declaration at its
original path. The recognized receipt binds its path and exact content identity.
Existing FNV/length identities are noncryptographic corruption checks, **not
signatures or proof of authenticity**. Self-consistent resealed unlogged work
cannot be authenticated by structural checks alone; the independent pins close
that corruption gap only while the owner keeps the declaration trustworthy.

The sealed checkpoint embeds a versioned conversion receipt containing old/new
bindings/schema, actual transform source, old root/ledger identity, byte-identical
ledger-prefix identity, budget history/authorized caps, old mode/counter baseline,
old occupancy evidence and streaming converted-state fingerprint. Stronger
non-scoring validation rejects malformed ownership, duplicate indices, counters,
native ledger discrepancies, geometry/cursors and links. These additional checks
were not part of the historical v1 loader.

Conversion preserves tasks/IDs/operation payloads, FIFO/shallow order and clocks,
native objective bits, incumbent/support, sticky loss, cumulative counters and
family rank/cursor/quantum. Only obsolete newest-ready indexing is removed and
focus heads/links are initialized empty. `historical_accounting` preserves old
third-share modes and v1 occupancy; `focused_services_by_family` reports only
post-conversion focused service. New occupancy/peak use v2 accounting rather than
relabeling old high-water values. Converted occupancy **plus the required next
operation reservation** must fit the explicitly authorized cap. Failure never
silently grows caps or evicts work.

Subsequent v2 ancestry crosses policies only through the recognized receipt,
re-verifying its exact old pinned ancestry and conversion conservation. Every
other epoch binding must match exactly. Keep ancestors and failed/interrupted
outputs; incomplete writes are not successful transitions. No v1/v2 trace equality
is claimed. No real checkpoint conversion or practical real-panel acceptance has
been authorized or demonstrated by these synthetic tests.
