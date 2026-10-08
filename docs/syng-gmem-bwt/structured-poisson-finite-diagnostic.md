# Gate A: standalone finite Poisson diagnostic

This addition implements **only Gate A** of [the committed plan](structured-poisson-optimization.md). It does not change production, the observation operator, dependencies, or the completion prototype. It emits no sequence or genotype. Gate B remains **NO-GO pending independent review and parent authorization**; these results do not demonstrate genome-scale recovery or tractable genome-wide pricing.

## Public seam and finite domain

`examples/panel_route_convex_diagnostic.rs` accepts a retained **synthetic** fixture from the unchanged `test_panel_route_policy` finite513 or `test_panel_route_policy_coupling` 64x64 generator. Its CLI rejects other layout/source-length/read-length shapes before native recompilation. This is a small diagnostic, not an input route for biological data.

The historical compiled envelope checksum is verified, without spoofing its compiler identity. The current public `joint::compile` recompiles the complete layout using the source-bound public panel graph and full registry. Every alternative's integer orientation totals and admitted starts must exactly equal the historical public finite oracle. Public graph loading/evaluation verifies current compiler, sample, panel, source/index and native-profile bindings.

Every Cartesian choice vector is enumerated, including conflicts. An independent all-pairs source-forward overlap check must agree with `Layout::feasible`, **and** the public panel evaluator's complete assignment validation. A mismatch fails: distinct explicit instance IDs can assert copies that the canonical source-span route model does not permit. No capacity model is weakened. The coupled fixture uses atomic spans, making these two declared models agree. The finite513 generator exhausts its small source-port traversal domain; this diagnostic's own bound scope is precisely the supplied, exhausted explicit layout, not arbitrary unprovided routes.

For every feasible choice, orientation totals must agree. Selected physical counts are checked-added in `u64` **before** exposure, including checked histogram multiplication/summation. Conversion requires integers <=2^53. Exposure is computed independently as `10 * sum(H*q) / sum(L*H)`. Each complete assignment is evaluated by both public `joint::evaluate` and `panel_routes::Evaluator`; exact per-length counts, observed counts, signals, and objectives must agree. The feature universe is the union of compiled definitions (including all imported registry definitions and all generated profiles) and complete sample-positive support; every public realized/sample-positive factor is checked against that union. Registry-only zero/zero terms contribute zero. Realized zeros pay signal cost, and unsupported positives remain explicit.

Legacy joint/f7 absolute NLL is compared after subtracting the **same complete union's** background constant `sum(beta - C*ln(beta))`. It is not compared directly to the background-relative objective. Shared physical contributions enter one logarithm. The portable L150 public replay test also compares expected per-start integer counts to counts from actual singleton sample reads, tests forward/RC and duplicate-description invariance, checks distinct-source two-copy multiplicity, and rejects a deliberately incompatible capacity assertion.

## Frozen recipes and resource policy

Recipes/caps were written to the evidence package's `frozen-contract.json` before implementation/outcomes. Existing DNA1600bp completion and all prior scientific recipes remain unchanged.

| Case | Declared domain | Actual diagnostic histogram | Fixed depth/background |
|---|---|---|---|
| finite513 | all 513 supplied complete routes, three sources 400/400/120 bp | L=[150,500], H=[252,0] | 10 / 0.1 |
| coupled64x64 | all 4096 products, 64 canonically feasible assignments, two identical 220-bp sources | L=[220], H=[2] | 10 / 0.1 |
| portable public | seed113 220-bp source duplicated physically; forward/reverse/duplicate-description/distinct-source alternatives, plus derived two-copy checks | L=[150], H=[71] | 10 / 0.1 |
| analytic | explicitly declared scalar/vector columns, not biological candidates | none | beta0.1 |

L220 is the **unchanged coupled regression**, not a replacement for the L150 primary control. The old configured regressions keep their original scoring settings; this diagnostic and additional frozen-f7 solves use depth10.

Limits: 4096 Cartesian products; 1024 feasible columns/public complete evaluations per evaluator; 4096 features; 4,194,304 feature-column cells; 128 hull iterations; 64 deterministic derivative-bisection line-search steps per iteration; 600,000,000 charged optimization scalar operations; 128 MiB logical numeric state; 64 MiB input JSON files. A cap fails closed, retains a failed/incomplete status, and withholds full-domain bounds/support; it never silently prices a prefix as exhaustive. Complete feasible enumeration precedes optimization.

Pricing costs `columns * features` scalar products **per iterate**, not one cheap genome oracle call. The work counter conservatively charges that product plus `(2*64+12)*features` each iterate. Public replay, input decoding, exhaustive feasibility checks, profile construction and identity checks are separate work, not hidden inside that counter. The numeric-state accounting counts dense numeric cells and retained trace vectors, not Rust allocator/JSON representation overhead or native-reader memory; it is **not an RSS cap**. Actual `/usr/bin/time -v` RSS/wall time are retained separately. The command wrapper enforces bounded logs; no runtime/memory scalability claim follows from these tiny cases.

## Mathematics and outputs

The loss, residual and full-change formulas are independent code in `math.rs`:

- `F(s) = sum(s - C*ln1p(s/beta))`;
- `g = 1 - C/(beta+s)`;
- exact full delta `sum((t-s) - C*ln1p((t-s)/(beta+s)))`.

At each iteration, every complete column is priced by `g dot column`, deterministic index order breaks pricing ties, and a bounded derivative-bisection line search moves along the hull edge. Current signal, residual, all column prices, raw gap, numerical tangent lower bound, convex mixture weights, objective and next step are retained. The best integral genome and all correlated public-score ties are selected separately by exhausted complete evaluation, never by interpreting mixture weights as copies, ploidy, support or genomes.

Each reported tangent bound is checked against the exhausted public discrete minimum. Known analytic tests also check it against the known continuous minimum. The public fixtures happen to reach integral hull optima; **the fractional integrality-gap demonstration is analytic**, C=1 with columns0/2 and optimum signal0.9. Another analytic example has a negative derivative from0 toward100 but a worsening finite step. Restricted-column optimization is explicitly not a global lower bound. Incomplete pricing returns `null` gap/bound and no complete correlated support, even if its apparent candidate gap is zero. The coupled public fixture demonstrates infeasible single moves and distinct count-equivalent complete exchanges; the improving two-move vector example is labeled analytic. Unchanged no-simple-return and DNA completion gates are regressions, not new optimizer discovery claims.

Three deliberately separate numerical tolerances:

1. Absolute **1e-9** public complete-score tie epsilon, used only for correlated discrete support.
2. Formula/bound comparison: **1e-7 * (1 + max(abs(a), abs(b)))**. This cannot merge discrete optima or establish epsilon-optimality.
3. Numerical hull stopping: raw gap <= **1e-9 * (1 + abs(F))**. Derivative checks use central h=1e-6 at the midpoint of each tested full change.

These are **ordinary f64 diagnostics**, not interval/logarithm-directed rounding or rigorous real-arithmetic certificates. Gap values are not clamped. A tiny negative raw gap or bound excess within comparison tolerance remains visible. Exhaustive minimum checks do not convert f64 arithmetic into rigorous certification. No result is a genome-wide bound, biological identifiability statement, posterior or calibrated confidence.

Per fixture output:

- `compiled-public.json`, `features.json`: exact public profiles and complete union/observations/histogram/background shift;
- `feasibility.json`: every product, both independent feasibility decisions, public rejection reasons;
- `hypotheses.json`: every feasible choice, normalized physical routes, full integer profiles, public/independent scores and actual discrepancies;
- `identities.json`: every column's exact change from column0, public score delta, midpoint derivative/finite difference, raw discrepancies;
- `relaxation.json`: all pricing/residual/bound/mixture traces and discrete support;
- `summary.json`, `status.json`, `frozen-profile-comparison.json`: scope, counts, timing and provenance.

Integer count-equivalence groups retain all column indices; inspect physical routes to distinguish descriptions from distinct assignments. **Equal scores alone do not establish count or DNA equality.** Complete correlated optima are never replaced by Cartesian products of marginals.

## Commands

Use the existing offline native environment, CPUs252–255/nice10, four jobs/threads, serialized tests. No install or staging.

```sh
cargo test --offline --locked --release -j4 --test test_panel_route_convex_diagnostic -- --test-threads=1
cargo build --offline --locked --release -j4 --example panel_route_convex_diagnostic
# Both arguments must name small synthetic fixtures / fresh outputs.
cargo run --offline --locked --release -j4 --example panel_route_convex_diagnostic -- \
  --fixture /absolute/fresh/finite-final --out /absolute/new/finite-diagnostic
# Requested configuration is mandatory, never a fallback to historical/private paths.
IMPG_TEST_CONVEX_FINITE=/absolute/fresh/finite-final \
IMPG_TEST_CONVEX_COUPLED=/absolute/fresh/coupled-final \
IMPG_TEST_CONVEX_OUTPUT=/absolute/new/diagnostics \
cargo test --offline --locked --release -j4 --test test_panel_route_convex_diagnostic \
  configured_finite_and_coupled_public_pricing -- --ignored --exact --test-threads=1
```

Optional `IMPG_TEST_CONVEX_ANALYTIC_OUTPUT` and `IMPG_TEST_CONVEX_PORTABLE_OUTPUT` preserve fresh directories from portable tests; default tests require no private paths. The worker evidence package contains the fully configured executable wrappers, source snapshots/diffs/modes/hashes, protected-file/installed-binary verification, independent trace verification, measured logs and failures:

`target/experiments/genome-mem-bwt-pipeline/candidate-copy-repair/structured-poisson-finite-worker-final-v1` in the shared target tree.
