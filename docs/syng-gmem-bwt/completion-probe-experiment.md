# Standalone completion-probe experiment

This **one-shot, budget-bounded proposal subset** reads an immutable compatible
v3 checkpoint. It does not advance, replace, prune, resume, or import scores into
production search. In particular a failed hypothetical native tail has no effect
on ordinary coupled next-molecule work. No sequence output, biological support,
optimum, genome recovery, QV, or real-data coverage is authorized.

## Usage

Build with existing dependencies and example auto-discovery:

```sh
cargo build --release --offline --locked --example panel_route_completion_probe
# Use a NEW output directory. These are synthetic/example input paths.
target/release/examples/panel_route_completion_probe \
  --panel fixture/panel.syng --routes fixture/routes \
  --sample fixture/sample.membwt --snapshot fixture/ordinary-snapshot \
  --out-dir fixture/probe
```

Controls are fixed at depth 10, background 0.1, and exactly read length 150.
The snapshot must have complete native initialization, current v3 policy identity,
and matching graph, sample, count policy, and controls. Evaluator feature/cache
limits are taken from those bindings, not silently changed. Normal public graph,
panel, sample, source/FAI/GZI, port, and native profile verification remains active.
There is no expected sequence, desired coordinate, family rank, or donor list input.

Defaults (explicit command-line options may lower/choose budgets for a separate
experiment; the acceptance recipe never tunes them after outcomes):

| Resource | Default |
|---|---:|
| charged work | 5,000,000 |
| fresh evaluator calls | 100,000 |
| logical reservation limit | 134,217,728 bytes |
| active task producers | 16 |
| active closure chains | 16 |
| bounded compact input record | 1,048,576 bytes |
| distinct normalized physical assignments | 200,000 |

Each producer retains at most one reconstructed entry and one target-interval
cursor. Full queues apply backpressure; there is no eviction or reservoir.
Discovery and active service alternate when both can run; active service
alternates producer work and round-robin closure chains. Producer work alternates
source/bound RR with ready-member grants when both can run. Total producer
occupancy includes both source and member-ready queues; the width is unchanged.
Neither discovery nor an entire source/member range must finish before closure
chains receive service.

Revision `experimental-snapshot-completion-ready-family-v2` groups producers with
`Members { next < end }` into owned per-family FIFO queues. At a free chain slot,
a single ordered index chooses the **least member-attempt-granted currently-ready
family**, breaking ties only by saved complete-native objective order and stable
family ID. No native rescore, partial score, source priority, or family exclusion
is used. Each grant inspects one member and increments that family's counter,
including rejected members and chains that later fail capacity or reuse a native
or cached score. A resource-failed attempt also consumes its grant without advancing
the retained member. There is at most **one in-flight chain per family**; ownership
ends only on actual chain termination, not on a resource stop. Unready families
are absent from the eligible index and cannot impose an all-family barrier.

Ready ownership transfers preserve cursors and producer payloads. A charged
exhausted-member step returns to ordinary source service. Only the affected family
key is updated at enqueue/grant/retirement; there are no full-pool scans or
family-wide rekey loops in service. This bounds repeated service to one ready
family but does not promise every family a completion within arbitrary budgets.
Every rejected record, eligibility span/route boundary, prefix position, binary-bound comparison,
source ordinal (including virtual/rejected positions), member, blocked producer
turn, native-tail append, and whole-assignment capacity comparison costs a turn.
An evaluator call is atomic and separately capped; work is not CPU instruction
count or a wall-time guarantee.

## Input and proposals

`stream.rs` first makes a bounded syntax/EOF/seal preflight pass. It retains no
candidate/task index. Its bytes/time are reported separately, not called discovery
work. Native baseline assignments and score bits are checked against the graph
and streamed initialization ledger rows without rescoring. Non-native ledger,
incumbent, support, and saved Evaluate/Probe payloads are **never** proposal seeds.
The reader understands the current task schema; discarded scheduling indices are
syntax-checked, not replayed or claimed to be a replacement v3 resume validator.
Snapshot parent/budget/transition provenance is reported, not migrated. The
original executable remains the authority for ordinary ancestry/resume validation.

A second reader retains a task cursor at the task-map offset. Each task inspection
is charged. Native-only contexts consume discovery/eligibility service but not
producer slots. Eligibility checks one span or route boundary per turn, including
same-source gaps/jumps, completed native extent, and live-cut discontinuity. Each
molecule uses its own native source/extent and resets at zero; contiguous same-source
storage splits remain native-only. These checks do not replace final validation. Positive
HubBound/HubScan contexts receive direct service. Earlier entry positions are
recovered one at a time from committed segments, keeping exactly the prefix before
that entry and earlier complete molecules, never its later detour.

The current production read-only fixed-width adapter is reused through a Rust
module path; no production implementation is copied/modified or called as a
scheduler. Positive alternative exits use its paired source ordering. Target hub
bounds compare `(full_DNA_word, decoded_numeric_source_id)` in the existing global
sorted records, not little-endian byte strings. No homology/repeat placement,
first-exit/widest-only restriction, frequency filtering, or truth-derived ordering
is assumed. A positive target suffix and incremental remaining native routes form
a whole assignment; every pair of canonical physical spans is then checked across
all molecules/orientations. Only complete feasible assignments reach the unchanged
public evaluator and its real full-word/sequence checks.

## Accounting and output

Logical payload reservations are compact wire length times 16 plus 4096 bytes of
node/cursor slack. Whole-chain reservations include the entire possible native
tail before construction. Producer reservations include both retained port words.
Queues, native metadata, cache keys, score-tie indices, and fixed metadata are
reserved before retention. Ready-service metadata still reserves 256 bytes
per native family and `(8 * size_of(Producer) + 256)` bytes per producer slot,
before allocating family rows, FIFO storage and eligible keys. Each owned ready
FIFO is now a standard `LinkedList<Producer>`: one owned node per element, with
removed nodes deallocated by `pop_front`. Empty or partly drained families retain
no historical backing capacity. The total live producer width therefore bounds
aggregate owned list nodes across **all** families, independently of churn.

A conservative node allowance is `size_of(Producer) + 2 * pointer_size +
2 * max_alignment` (two links plus alignment padding). The unchanged per-slot
reservation covers this node storage, a locally moved producer, and even an extra
transient node. Producer payload buffers remain separately reserved and are moved,
not cloned, during FIFO transfers. Global source/chain deques retain their original
fixed-width reservations. Freed allocator caches/size-class overhead are not owned
list storage and are not an RSS claim.

The previous per-family `VecDeque` representation was incorrect for this bound:
`pop_front` retained each family's high-water capacity, which could accumulate
across many families despite bounded live length. No cap or scheduling change is
used to correct this; only the owned ready-FIFO representation changes.
An additional `16 * record_bytes + 135168` reservation
bounds parsing/current-input transient payloads and reader buffers before input is
opened. Overlarge records explicitly fail rather than allocating an unbounded
whole JSON document. A resource stop preserves/report pending cursors in the
one-shot result; it is **not** a resumable sidecar.

These are conservative logical reservations, **not RSS**. Public graph/panel/sample
loading, native source handles, evaluator cache and atomic evaluator scratch are
outside this accounting; measure RSS/time externally. Streamed coordinate results
consume disk, not a retained result list. Huge snapshots require a full sequential
validation read before charged discovery, with no finite-budget discovery promise.

Outputs are exclusively coordinate/score/accounting artifacts:

- `options.json`: invocation and caps;
- `provenance.json`: actual snapshot seals, scientific/policy bindings, saved native
  score bits, ancestry/budgets, preflight accounting, explicit method/service revision;
- `proposals.jsonl`: every evaluated/reused complete assignment, objective bits,
  bit-preserved native deltas, task/entry provenance, input position/unfinished flag,
  fresh/cache/native visit classification, distinct-physical and score-tie flags;
- `result.json`: explicit stop, inspected/skipped/admitted and stage counters,
  pending source/member-ready producers/chains/input, per-family grant/in-flight/ready
  counters, logical occupancy and output checksum;
- `failure.json`: input/initial resource failure (no partial scientific result).

Caching only reuses computation for exactly identical normalized physical
assignments. Equal DNA or equal score never collapses different physical resources.
Saved native baselines and native reuse visits are reported separately. Outputs
are flushed before post-run snapshot/ledger mutation checks and final result.

## Tests and frozen scientific gate

`tests/test_panel_route_completion_probe.rs` contains the explicitly configured,
default-ignored `completion_probe_frozen_non_native_mosaic` gate. It requires
`IMPG_TEST_COMPLETION_OUTPUT` (fresh) and `IMPG_TEST_COMPLETION_V3` (preserved normal
v3 executable). It freezes coordinate outputs before test-only whole-DNA spelling.
The recipe combines only B's middle 24-base variant block with A's flanking alleles
and A's independent second molecule. Actual L150 MEM-BWT reads are made every 15
bases (including the last possible start); every native genome differs from truth.
The ordinary snapshot uses 100,000 work / 10,000 evaluations / 256 MiB /
20,000 optima. No expected route/context is injected into that end-to-end gate.

Independent assessment enumerates **all** L150 windows of complete spelled
molecules into a separate public sample MEM-BWT, checks exact integer counts
against public route counts/factors, and computes the union objective in a separate
loop. This is independent of native profile/route replay and evaluator loss, not
an independent implementation of the shared public read/count operator. Unsupported
positive features and realized zero-observation penalties remain present. Existing
configured finite f7 gates separately retain their independent implementation oracle.

Portable mechanism tests use **explicitly labelled injected contexts**, not the
end-to-end recipe: a current tip with no target return, earlier entries retaining
two distinct donor interiors, deferred positive exits, and 10,000 native distractors
still undiscovered when closures finish. They also check whole-genome capacity
before scoring, numeric IDs across the 255/256 boundary, streaming EOF/record bounds,
and reservation failure without growth. The configured mosaic gate adds corrupt
copy binding/seal/EOF/native-bit checks, storage/work/evaluation/distinct stops,
original-byte preservation, and unchanged-policy resume. All existing production
Rust/Python tests and configured finite/coupled/no-return/reverse/ancestry gates
remain unchanged.

The first frozen scientific run recovered exact whole synthetic truth and strictly
beat all native baselines. Its ordinary input snapshot had **already** recovered
that same truth/score. Thus this establishes standalone synthetic correctness, not
novel improvement over ordinary search. The injected mechanism cases establish
service of earlier/deferred prefixes separately. No real-data result is implied.

The portable `completion_ready_family_causal_service_and_outcomes` test freezes
three already-admitted, target-member-ready producers (A, A, D), a live B source
producer and 10,000 native distractors under two fresh evaluations. A has a strictly
better complete native score than D. Isolated original producer-RR selection uses
both fresh calls on A; ready-family service evaluates A and D under the identical
caps and unchanged member/tail/capacity/evaluator primitives. The test checks local
queue/cursor conservation, stopped-chain ownership, and equal grant accounting
for fresh/cache/native/capacity-rejected/member-rejected outcomes. Set
`IMPG_TEST_READY_FAMILY_OUTPUT` to a fresh path to retain this mechanism evidence.
It injects readiness for a scheduling regression only; the original normally
generated 1600bp DNA gate and its thresholds/caps remain unchanged.

The frozen `ready_family_storage_churn_releases_nodes_not_historical_capacity`
regression fills/drains32 families at width16 for four passes, then performs128
bounded-width rotations with partially drained and concurrently nonempty families.
It records **actual old VecDeque capacities**, reproduces their reservation
violation, and checks new aggregate node/transfer bounds using LinkedList's node
ownership/deallocation guarantees (not a length-as-capacity assumption). FIFO and
local key conservation are checked throughout. Set `IMPG_TEST_READY_STORAGE_OUTPUT`
to a fresh JSON filename to preserve this evidence. Existing failed-reservation,
ready-cursor, chain-ownership, causal and frozen scientific gates are unchanged.
Old/new ready-service executable proposal ledgers must also be byte-identical on
the same frozen synthetic snapshot: the storage correction is not a policy change.
