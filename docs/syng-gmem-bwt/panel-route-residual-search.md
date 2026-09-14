# Gate B1: automatic two-ended proposals with exact global feedback

This standalone experiment implements the **approved B1 subset**, not all of
[Gate B](structured-poisson-optimization.md). Existing tracked code, production
CLI, operator, dependencies, Gate A and completion prototypes are unchanged.
There is no installation, warm start, truth/assignment input, sequence emission,
fractional genome, global pricing bound or complete-support claim.

## Inputs and candidate construction

`examples/panel_route_residual_search.rs` accepts fresh public panel/routes/sample
artifacts and explicit resource limits. It initializes complete native assignments
for every supplied endpoint inventory, interleaving initialization with already
ready work. Native scores and their full route-profile cost are charged. Budgets
may stop initialization; `native_initialization_complete` distinguishes this.
There is no all-family exhaustion barrier or native-score donor filter.

For a complete immutable baseline, each molecule has a lazy coordinate schedule:
level l has n=2^l, and cells `[m*j/(2*n), m*(j+2)/(2*n)]` for j=0..2*n-2.
These are dyadic cells plus half-cell shifts. Ends snap to actual oriented ports
using binary search in the public verified source index, inspecting nearby
records of the required orientation. Entry and exit full-DNA words define hub
ranges in the verified global index. Their product is traversed **incrementally**,
not materialized as a patch/edge matrix. Two members on one donor lane with equal
orientation and positive traversal define a donor segment. A complete molecule
is constructed from the original prefix, donor and original suffix. Endpoints,
join DNA and global canonical capacity must pass public complete validation.
All donor interiors/strands remain eligible; a poor whole-native donor can supply
a beneficial tract. This is a sampled, capped neighborhood, not a full oracle.

Candidate construction and scoring are separate. A rejected single patch still
has its reciprocal successor queued. When the donor span is occupied wholly
inside one segment of another molecule, B1 can replace that span by the removed
subwalk, simultaneously patching both molecules. No hypothetical single genome
is treated as the feasibility test for the final reciprocal assignment. The
**final compound assignment** must validate before profiling or scoring.

### Important oriented-coupling limitation

The reciprocal constructor supports **same traversal orientation only**.
Independent reverse-donor single patches remain eligible. Opposite-strand
reciprocal attempts generate an explicit `unsupported_opposite_strand_reciprocal`
event/counter, not biological infeasibility. At odd k, forward cut=floor(k/2) and
reverse cut=k-floor(k/2); reversing a cut subwalk need not produce declared ports
at those exact canonical coordinates. No base is shifted to force validity.

The derived seed113 reverse fixture has k=63 from its actual graph. A NEW causal
replay of the original root enumeration and queued single/reciprocal successors
records graph-derived ports/cuts, parent/task IDs, candidate assignment and the
exact public error `undeclared outgoing switch`. It is **not a recovered historical
failed object**: the original failing implementation did not serialize that object.
The original executable/logs and mistaken reconstructions are preserved with this
limitation. A middle-region opposite swap is separately shown **publicly feasible**
yet unsupported by B1. Thus the restriction must not be generalized to all
opposite-strand swaps or to the overall candidate model.

Oriented two-ended coupled construction and two-donor/no-simple-return interior
chains remain unfinished Gate-B work. Passing the unchanged old no-simple-return
regression protects old functionality; it does not show this generator discovers
those solutions. A focused geometry-enabling experiment is needed before claiming
full Gate-B acceptance.

## What feedback changes

Every eight publicly confirmed candidates on one immutable baseline, the two
lowest exact whole-genome scores select endpoint-refinement neighborhoods.
Selection includes non-improving candidates: coarse exploration continues
independently, and worsening immediate steps do not erase future work. Refinement
moves either endpoint by half the parent interval width, then re-snaps to actual
ports. The next generated/attempted coordinates, not merely final rankings, are
recorded. This is direct evidence-responsive candidate construction.

A strictly improved complete family incumbent can create a new baseline and
coarse stream, up to the declared expansion cap. All old tasks retain their
original complete-baseline IDs/counts; no cursor or sparse delta is relabeled.
There are no imports of saved scored assignments, old ledgers or snapshots.
A profile may be reused by normalized complete route; a score/delta is never
cached across changed baselines. All public scores remain in the streamed ledger.
Distinct tied physical assignments are retained, never interpreted as independently
combinable marginals. If in-memory tie storage fills, the experiment stops; the
newly scored assignment remains in the ledger, and support is explicitly incomplete.

## Public scoring seam and arithmetic

The unchanged public APIs are `Evaluator::validate_assignment`, `::route_counts`
and `::evaluate`. Complete route profiles include unioned junction/terminal
corrections and both native orientations. `route_counts` alone is not complete
genome feasibility; it is called for a proposal only after final validation.

Baseline integer counts come from complete public factors. For changed molecules,
checked integer old profiles are subtracted and checked new profiles added before
histogram exposure conversion. Sparse touched-feature union includes newly realized
zero-count features; every observed sample-positive is retained. Shared molecules
contribute inside **one** logarithm. Fixed depth10/background0.1 and unchanged MEM
occurrence semantics apply. Checked integer histogram products precede f64 exposure;
conversion is rejected beyond the declared exact integer range.

Once exact changed-molecule profiles have been computed, the sparse exact loss
change is nearly the same traversal as a residual dot product:

`sum[(t-s) - C*ln1p((t-s)/(0.1+s))]`.

Therefore B1 uses the **exact global delta** for guidance, rather than calling a
redundant gradient pass cheap pricing. Diagnostic `g dot delta` is accumulated in
the same sparse pass, but is not a candidate-acceptance rule or pricing certificate.
**Every exact sparse delta counts as a complete score evaluation.** Each proposal
also receives an independently counted public complete confirmation. Public integer
counts and objective must agree; comparison tolerance is1e-7 scaled, while discrete
score ties use absolute1e-9. Ordinary f64 arithmetic is not a rigorous certificate.

The ledger separately records native scores, every computed unconfirmed exact
delta, every returned public confirmation, and parity-checked candidates. If a
budget stops between delta and confirmation, the unconfirmed score is retained
but not accepted as a result. Unexpected validation errors record only the current
task/candidate/error and fail; they are not swallowed as capacity rejection.
The complete known touched-feature pass is admitted atomically before charging
its score attempt; completed delta counters are committed only after construction.
A work stop cannot count an unfinished delta as a completed objective.

## Frozen controls and limits

The recipe/caps were approved/frozen before outcomes. Automatic control:

- xorshift seed219 DNA1536 for A; B copies A with24-bp mutations beginning at
  256,704,1216 (`A→C`, other bases→A);
- passive400-bp paths from seeds57/58 for A/B respectively;
- sample truth is A with **only** the middle704..728 B tract, plus A's passive path;
- public L150 reads at starts step15, including the last valid start;
- CLI sees only panel/routes/sample/budgets. Inference files are fingerprinted
  before any DNA comparison. This is a synthetic development control, not blind
  genome validation or proof of global biological identifiability.

The fixed-native mechanism holds graph, legal geometry and complete baseline A
constant while counts favor the middle versus last tract; both selected child
coordinates and actual subsequent geometry attempts must differ. This injected
helper is separate from the automatic CLI and has no corresponding CLI family or
assignment input. The220-bp identical two-copy test is a physical/capacity mechanism,
not identifiable sequence recovery. Derived reverse/readiness fixtures test
operator/orientation and interleaving laws; they do not replace the frozen recovery
control. All prior DNA1600 and Gate-A recipes remain unchanged regressions.

Defaults and hard maxima:2m work units;2m profiled-base-length units;2048 exact
score admissions **including native/delta/public confirmation attempts**;1024
validation admissions;128 pending tasks;64 retained ties;64 segments per route;50k feature terms;
128MiB logical accounted state;4 incumbent expansions per family;dyadic/refinement
level8. Lower caps are configurable and explicitly tested. No outcome-based cap
or recipe extensions occurred.

Profile work is precisely `sum complete route lengths * compiled read-length count`
for every admitted explicit `route_counts` call and every route reserved for an
admitted `evaluate`, including cache hits, native initialization and confirmations.
An aborted public evaluation keeps the complete prepaid profile charge; it is
not refunded or reported as a completed objective. This is conservative admission
accounting, not a count of actual internal base traversals on failure.
It is not merely corrected-junction length. `accounting.jsonl` preserves every
profile charge and state reservation. The initial two-native control costs3872
units. A conservative first-eight native-baseline candidate cycle plus those
natives costs at most68,512 units, below2m, checked before implementation.

Port decoding, task operations, sparse feature loops and validation/copy upper
charges consume work units. Validation charges segment-pair and indexed-port/DNA
lookup upper counts; it does not pretend to expose exact backend instruction
counts. State reservations prepay cleanup work (`ceil(bytes/64)`) and account for
retained vector capacities, queue/support slots, baseline maps, public profile/cache
workspace and temporary transfers. Cleanup credits dominate some small runs and
are disclosed, not hidden as useful search effort. State/profile caps are logical
accounting, **not an RSS/instruction/time guarantee**. Public provenance/native
loading and sealed-file I/O remain real costs; input byte counts, call/cache-hit
counts, native correction starts and wall/RSS are separately retained. No new
whole-domain RAM port index or eager patch-pair matrix is constructed.

`score_attempts` enforces the unchanged score-admission budget; `exact_scores`,
`exact_delta_scores`, `public_scores` and `native_scores` count only complete
returned objectives. `public_score_attempts`/`native_score_attempts` count admitted
public calls. `explicit_profile_attempts` counts entered `route_counts` calls;
`public_evaluation_route_reservations` counts routes prepaid for `evaluate`.
`completed_route_profiles_observed` counts successful explicit returns and routes
in successful evaluation returns only: the public API does **not** expose internal
profile completion in a failed evaluation. It is not a total completed-profile
claim when evaluations abort. Cache hits include exposed hits on failed calls;
corrected starts include successful returns only. `validations` is a budget
admission count, not a successful-validation count (a later profile-admission
stop can occur before the reserved validation is entered).

Every cap denotes a resource limit, never biological infeasibility. Once the engine
is initialized, local budget errors and exactly the public InvalidData messages
`native feature resource cap exhausted (no truncation)` and
`realized feature union exceeds resource cap; no partial objective` finalize a
bounded stop. Unrelated public errors remain failures; input/preflight failures
before engine creation remain explicit failures, not invented pending searches.
Pending candidate work/baselines, the active native family and **all** uncompleted
native family IDs are retained in `pending.json`, including a family whose public
call failed. Native cursors advance only after initialization finishes. There is
**no resume contract** and no partial objective/profile is fabricated.
The wider finer/oriented/chain domain remains unassessed even if the capped declared
schedule ends. `global_bound`/`global_gap` stay null, support completeness and
sequence-emission permission stay false. Stagnation is not a genome optimum.

## Execution and evidence

Use the existing offline native environment, taskset CPUs252–255/nice10, at most
four jobs/threads and serialized tests. No install or staging:

```sh
cargo test --offline --locked --release -j4 --test test_panel_route_residual_search -- --test-threads=1
cargo build --offline --locked --release -j4 --example panel_route_residual_search
cargo run --offline --locked --release -j4 --example panel_route_residual_search -- \
  --panel /fresh/panel.syng --routes /fresh/routes --sample /fresh/sample.membwt \
  --out-dir /new/output
IMPG_TEST_RESIDUAL_BIN=/absolute/built/panel_route_residual_search \
IMPG_TEST_RESIDUAL_CLI_OUTPUT=/new/frozen-cli-fixture \
cargo test --offline --locked --release -j4 --test test_panel_route_residual_search \
  configured_cli_automatic_mosaic -- --ignored --exact --test-threads=1
```

Optional `IMPG_TEST_RESIDUAL_OUTPUT` preserves portable fixture directories beneath
an existing fresh root; `IMPG_TEST_RESIDUAL_BOUNDARY_OUTPUT` preserves the focused
feature-cap and analytically derived late-work boundary fixtures. Missing requested private configuration fails, never falls
back to a path. Event-driven log-size alerts preserve failures and all buffered
diagnostics; no file-size limit is applied to scientific artifacts.

Evidence packages: shared target `experiments/genome-mem-bwt-pipeline/candidate-copy-repair/structured-poisson-gate-b-worker-final-v1`
remains immutable; sibling `structured-poisson-gate-b-worker-final-v2` contains the
review corrections, red/green boundary cases, narrow v1 diff and fresh regressions.
Each contains source/modes/hashes, executable, frozen contract, approved geometry
limitation, causal before/after mechanism traces, all score/cost ledgers, preserved
failures, full portable/configured/Gate-A/DNA1600/Python regression commands/results,
and independent verification. The v2 audit independently checks **score-record
arithmetic and emitted profile/state reservation sums**, not all costs. Primitive
work, validation admissions, actual call/cache-hit counts, opaque failed-call
profile progress and allocation-bound properties are **not independently
reconstructed**. Checking that reported meters fit configured caps is not proof of
those charges or an allocation bound. The old v1 all-costs flag was too broad.
Old ordinary/prototype baselines are not imported
or assessed on the new1536-bp control; old DNA1600 comparisons remain separate.

**Exit:** B1 can be reviewed for automatic tiny non-native improvement/recovery,
evidence-responsive refinement, exact complete score parity and supported compound
capacity behavior. Full Gate B and genome-scale use are not accepted by these
results. Next: reviewed enabling geometry experiments for oriented coupling and
two-donor/no-simple-return construction, then separately approved development-data
experiments; Gate C still requires new frozen blind challenges.
