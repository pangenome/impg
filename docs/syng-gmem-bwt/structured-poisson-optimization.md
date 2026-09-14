# Whole-genome mosaic inference: structured Poisson optimization and completion plan

Status: **Gate A and bounded B1 accepted; full Gate B incomplete, 2026-09-14**. The observation operator and complete-genome scorer exist; reliable genome-scale mosaic recovery does not yet exist. [Gate A's executable finite diagnostic](structured-poisson-finite-diagnostic.md) was independently reviewed and parent-reproduced (commit `9d8cc24`): public counts, exhaustive feasible optima/support and numerical bounds agreed; 699 portable Rust tests and the configured gates passed. This is finite numerical evidence, not a genome-scale certificate. [B1's score-guided generator](panel-route-residual-search.md), commit `594a6f3`, now constructs new complete mosaics and changes subsequent refinement geometry in response to sample counts. Independent review and parent reruns passed, including 706 portable tests and the original resource-boundary red/green cases. It recovered the exact frozen 1,936bp two-molecule synthetic truth, improving the objective by 67.5373 over every native. This is not whole-yeast recovery. Next: oriented coupled and two-donor/no-simple-return construction, then a separately declared practical genome-scale resource/development experiment.

## Goal and scope

Infer complete, physically feasible haplotype-mosaic genomes from sample MEM-BWT occurrence counts against the syng panel. Use all panel donors and the supported topology inventories; do not reduce inference to choosing a reference. Preserve correlated alternatives, source coordinates, orientation, continuity and canonical-span capacity. Assess reconstructed hypotheses with **both sequence accuracy and complete truth/query coverage**.

The current panel contains 235 sample–haplotype assembly identities. These provide candidate native path inventories/endpoint pairings, not certified chromosome biology or ploidy. Novel topology/karyotype, noisy reads, stock/ploidy validation and statistical calibration require separate gates. Completion of the initial error-free, supported-topology workflow must not be advertised as completion of those broader claims.

Related implementation contracts: [finite joint walks](joint-walks.md), [panel routes](panel-routes.md), [guided policy](panel-route-guided-policy.md), and [completion prototype](completion-probe-experiment.md).

## Where we are: biological result, not service counters

The standalone completion service now generates and globally scores complete donor substitutions, including a 451-kb replacement in a native SK1 control. Native control scores remained best among tested proposals. This established useful candidate construction, not mosaic recovery.

The first blind whole-yeast mosaic challenge then failed recovery:

- 12,233,923 bp in 17 candidate molecules; one internal 104,821-bp SK1 tract replaces 113,840 bp of S288C chrIII.
- Recipe, sampling and budgets were frozen first. Inference received the unchanged panel and sample counts, not truth breakpoints/assignment. The hybrid whole molecule is absent from native panel molecules in either orientation.
- Fresh ordinary v3 search: 6m work, 1,000 evaluations, 1 GiB logical state. Standalone probe: 5m work, 256 fresh evaluations, 128 MiB logical state, 256 producers/16 chains. No extensions or outcome-based recipe tuning.
- Ordinary search stayed at the best native. The probe found a slightly better-scoring but incorrect mosaic.

| Frozen hypothesis | Relative objective, lower better | Alignment QV | Full truth coverage | Full query coverage |
|---|---:|---:|---:|---:|
| Native / ordinary incumbent | -10,473,857.958090 | 36.3255 | 99.7428% | 99.6721% |
| Probe representative | -10,474,553.566036 | 36.5136 | 99.5438% | 99.6771% |
| Truth oracle, evaluated only afterward | -10,506,105.570280 | — | — | — |

On the modified molecule, the probe had QV20.3610 and only88.4361% selected truth coverage, versus native QV20.5092 and95.7575%. The whole-genome QV obscures the missed tract because most of the genome is unchanged. Missing alignment is explicit, not automatically classified as a base error.

The truth scores31,552.004244 better than our best found genome. This establishes a better feasible solution that the search missed; it does **not** prove truth is the global objective optimum or uniquely identifiable.

Evidence: `/home/erikg/yeast/genome-blind-mosaic-20260914T121834Z/`, authoritative assessment `assessment-completion-v2/report.json`. The first truth-oracle request lacked version/model/graph-checksum serialization fields; its failure is preserved in `assessment-v1`. A separate envelope added only those fields, preserving frozen geometry, truth, inference, selection and alignments. It did not alter the scientific experiment.

This challenge is now a **development case**. A later run on it is not a new blind validation.

## 1. Mathematical model: convex loss over a discrete feasible set

For feature f, let C_f be its observed MEM **occurrence** count, e_f(G) the exact physical exposure of complete genome G, and s_f(G)=d e_f(G). With fixed d=10 and beta=0.1, retain the existing background-relative objective:

\[
R(G)=F(s(G))=\sum_f\left[s_f(G)-C_f\log(1+s_f(G)/\beta)\right].
\]

Up to constants on the fixed feature universe, this is Poisson negative log-likelihood / generalized KL divergence between C and beta+s. Background-relative evaluation leaves uninstantiated zero/zero contributions at symbolic zero. All realized zero-count features pay their signal cost; unsupported positives remain explicit. Shared physical contributions add **inside one logarithm**.

The derivatives are

\[
g_f(s)=1-\frac{C_f}{\beta+s_f},\qquad
H_{ff}(s)=\frac{C_f}{(\beta+s_f)^2}\ge0,
\]

with zero off-diagonal Hessian entries in signal coordinates. Thus F is convex in predicted signal. Independently unconstrained feature minima occur at s_f=max(C_f-beta,0), but such independently fitted rates generally do not describe any genome. An independent-feature bound can consequently be extremely weak.

MEMs overlap and their counts are correlated. Convexity is an algebraic fact about this objective; it does not imply an exact independent-Poisson generative likelihood, calibrated confidence, or applicability of statistical guarantees from independent observations.

### Discrete genome constraints

Let S={s(G): G is a legal complete genome under the declared candidate model}. The real optimization is min over S, not over arbitrary nonnegative vectors.

One conceptual formulation is s=A x, with integer variables selecting exact complete-molecule alternatives, endpoint-inventory choices and capacity constraints. Columns must include their full junction and terminal contexts. Assigning fixed profiles to arbitrary reference intervals and ignoring new junctions would be incorrect. Alternatively, complete-genome profiles themselves are columns; this gives a clean conceptual convex hull even when a compact molecule formulation is unavailable.

A is **implicit**: no eager genome matrix, read-by-haplotype matrix, artificial read placements or all-pairs expansion. Expected profiles continue to come from the exact public operator, with checked integer per-length counts before exposure conversion.

## 2. Residual weights guide proposals, not biological validity

A proposed complete change has exact signal difference Delta s. Its exact loss difference is

\[
F(s+\Delta s)-F(s)
=\sum_f\left[\Delta s_f-C_f\log\left(1+\frac{\Delta s_f}{\beta+s_f}\right)\right].
\]

The first-order quantity g(s)^T Delta s supplies a proposal objective:

- negative weights indicate underpredicted observations;
- positive weights indicate excess predicted signal;
- C_f=0 gives weight1, preserving the cost of unsupported signal;
- all weights refer to the current **complete genome**, including repeats/shared observations.

A negative first-order prediction does not guarantee improvement for a finite discrete move. Always compute the unchanged exact complete-genome objective before reporting an improvement. Convexity gives a tangent lower bound for an exact full signal change; it does not justify pruning partially specified routes using a donor-local likelihood or an approximate context profile.

This adds no switch penalty, family prior, annealing, fitted background/depth, feature masking or count-policy change.

## 3. Conditional gradient / column generation and valid bounds

Consider the relaxation min F(s) over conv(S). Fractional combinations are optimization devices, **not inferred ploidy, averaged genomes or permission to emit a consensus**.

For a point s in that convex hull, solve the linear pricing problem

\[
\widehat G\in\arg\min_{G\text{ legal}} g(s)^T s(G).
\]

An exact full-domain oracle gives gap

\[
\delta=g(s)^T(s-s(\widehat G)),\qquad
F(s)-\delta\le\min_{z\in\operatorname{conv}(S)}F(z)\le\min_{G}R(G).
\]

A certified lower bound ell on the pricing minimum also yields the valid tangent bound F(s)-g(s)^T s+ell. A merely found pricing candidate gives an **upper**, not lower, bound on that pricing minimum and cannot certify this gap.

Critical distinctions:

- A restricted-master minimum over only generated columns is **not** a lower bound on the full integer-genome problem in general. Missing columns can improve it.
- Exhaustive pricing certifies only the explicitly exhausted finite universe. No finite-fixture certificate becomes a genome-wide certificate.
- Even an exact convex-relaxation solution may have an integrality gap. The best legal genome remains a separate incumbent.
- Floating-point claims must state their arithmetic/tolerance scope; no rigorous real-number certificate is claimed from ordinary f64 calculations alone.

Linear pricing may become a weighted path/flow problem if sufficient read-context state is represented. Canonical-span capacity, multiple molecules and topology selection can still make it difficult. We will not assume an ordinary shortest-path or min-cost-flow oracle solves the actual resource-constrained problem.

## 4. What this explains about optimization behavior

**Local traps:** convex observation loss does not make single-edit genomic search globally convergent. Two compatible changes can help together while either alone hurts or is infeasible. Coupled/no-simple-return fixtures remain mandatory; a greedy-only hill climber is insufficient.

**Non-identifiability:** s(G1)=s(G2) makes the genomes indistinguishable to this objective, even if their sequences or donor origins differ. Near-equal rate vectors can also provide weak discrimination. A Hessian/conditioning diagnostic is not a calibrated uncertainty estimate, particularly with correlated MEM counts. Equal scores alone do not prove equal rates or DNA.

**HMM limitation:** standard copying-HMM/Viterbi models assume a local emission structure we do not have. Unlocalized MEM counts and repeat/copy sharing couple positions across the genome. Any dynamic program must represent those dependencies rather than invent locus assignments or independent local evidence.

## 5. Delivery gates, in order

### Gate A — executable finite-oracle mathematics (ACCEPTED, finite numerical scope)

Add an isolated experimental diagnostic; keep production and the accepted completion prototype unchanged.

1. Use explicit synthetic cases plus existing frozen finite/coupled fixtures with complete feasible enumeration. Preserve all correlated alternatives and actual resource constraints.
2. Obtain exact integer counts/exposures from the unchanged public machinery; independently reconstruct the fixed objective/gradient and exact score differences. Include zero counts, unsupported positives, junctions, reverse views and duplicate descriptions versus physical copies.
3. Exhaustively price all feasible columns in the small domain. Implement a bounded convex-hull conditional-gradient diagnostic, with deterministic line search, explicit tolerances and separate discrete incumbent.
4. Check tangent/gap bounds against the exhaustive discrete minimum and known analytic cases, at every tested iterate. Include a restricted-column counterexample showing why its optimum is not automatically a global lower bound, and a heuristic-pricing example that withholds certification.
5. Show a compound-move case, exact count-equivalent alternatives, and a finite-step overshoot where negative derivative does not mean improvement.
6. Freeze fixture recipes/caps before outcomes. Report wall time, oracle effort and failures as well as objective traces. No claim that the continuous iterate is a legal genome.

Exit: reviewed executable evidence establishes the formulas, the scope of every bound and the failure cases. If the mathematics/representation does not match the actual operator, correct that mismatch before genome-scale implementation.

### Gate B — structured complete-genome proposal/pricing experiment (IN PROGRESS; B1 accepted)

**B1 exit:** automatic single-donor replacements, score-dependent refinement and supported same-orientation exchanges passed the frozen tiny recovery/causal/physical-resource tests. Exact sparse complete-objective deltas drive feedback; the diagnostic linear residual price is not advertised as a cheaper oracle. Seven retained best physical assignments spell the exact test genome; support remains incomplete. Opposite-strand reciprocal construction and two-donor/no-simple-return generation remain unfinished. The frozen B1 resource caps cannot initialize a whole yeast genome; no genome-scale or full Gate-B acceptance follows from this result.

Use Gate A to choose the smallest scalable candidate generator: coarse two-ended replacements, boundary refinement and compound changes, with sufficient junction context and globally valid capacity. Residual weights prioritize **complete feasible changes**; exact full-genome evaluation accepts their reported scores.

Keep multiple alternatives and exploratory/coupled work. Do not use single-change failures to prune ordinary exchange paths. All235 endpoint inventories and all donor interiors remain eligible. No truth coordinates, correct-donor hints or saved scored assignments silently reinterpreted as a different proposal policy.

Explicitly specify initialization, neighborhood, exploration, retained state, work charging, cache reuse, stopping, and any warm-start contract before implementation. Preserve current searches/prototypes as immutable baselines; no silent production/checkpoint update. Avoid dense matrices and hidden all-pool scans. Exact pricing may be infeasible at scale: report heuristic pricing as heuristic and withhold unsupported bounds.

Exit: finite regressions still pass; the generator finds materially better complete hypotheses on development data under predeclared budgets. If coarse moves do not provide useful guidance or pricing dominates runtime, investigate that specific failure rather than adding budget or claiming expected convergence.

### Gate C — fresh blind whole-genome recovery

Freeze NEW sample-independent recipes, budgets, candidate-generation code and assessment rules. Include native controls and increasing mosaic difficulty: single tract, multiple tracts, more than one molecule, orientation and resource-coupled alternatives. Keep known non-identifiability separate from search failure.

Only panel/count inputs enter inference. Freeze hypotheses and correlated alternatives before truth assessment. Score the truth oracle only afterward. Report native versus inferred loss, exact whole-molecule equality where applicable, full truth/query coverage, alignment errors/QV, modified-molecule/tract performance, unsupported signal and resource limits. No whole-genome average may conceal failure on the altered part.

Exact finite identifiable fixtures must recover their optimum/support. Genome-scale acceptance thresholds and ambiguity handling must be declared per challenge before outcomes; current data must not be relabeled blind after tuning. Coverage loss is not repaired merely because aligned QV increases. Do not collapse score-tied physical origins into a certified winner.

### Gate D — supported inference and reconstruction product

After recovery gates: explicitly reviewed production integration and output/emission contract. Separate best-found, objective-bound scope, feasibility, support completeness, sequence ambiguity and calibrated confidence. Unknown/background preference is not deletion; no invented sequence, fixed-N joins, forced donor, majority orientation or tied-copy winner emission.

Provide reproducible counts-to-inference-to-assessment commands with practical time/memory limits and reports. Preserve unchanged-policy resumes and make any policy update explicit; never reinterpret an old cursor ordinal. Install production via cargo install only after acceptance. Regular trusted-local checks suffice; no migration/security framework expansion.

Exit: a documented, usable workflow for the declared supported scope, backed by blind recovery and full-coverage evidence—not a collection of successful internal service tests.

### Gate E — broader biological validation

Held-out assemblies/aliases, real noisy reads, stock identity/ploidy, dosage, calibration and novel topology require explicit independent experiments and models where needed. Track each as passed, unsupported or blocked. Do not declare general whole-genome inference complete while these remain implicitly assumed.

## Execution and decision discipline

- Sole source/build writer -> fresh read-only review -> parent acceptance/reruns/publication. Parent alone runs real data. Existing protected backend/dependency files and production policy remain unchanged in the next experimental slice.
- Use existing offline dependencies, at most4 threads/jobs, CPUs252–255/nice10, serialized tests, managed processes and bounded logs. Stream large inputs.
- Cache computation, not copy multiplicity; retain distinct normalized physical assignments and correlated alternatives.
- Report actual loss, valid bound scope, sequence results, elapsed time/RSS and logical retained state separately. Visits, fresh evaluations, physical assignments and DNA outcomes are different quantities.
- Preserve failed attempts. Correct administrative/schema mistakes transparently without changing frozen biological recipes or inference selections.
- Each gate ends with an explicit go/no-go and the smallest next experiment. A negative result is useful; indefinite scheduler polishing is not the project goal.

## References

- Boyd and Vandenberghe, *Convex Optimization*, Poisson estimation example: <https://web.cvxr.com/cvx/examples/cvxbook/Ch07_statistical_estim/html/counting_problem_poisson.html>.
- Jaggi (2013), *Revisiting Frank-Wolfe: Projection-Free Sparse Convex Optimization*: <https://proceedings.mlr.press/v28/jaggi13.html>.
- These references support the optimization framework. They do not establish our biological identifiability, independent-observation assumptions, or a tractable full-genome pricing oracle.
