# Partition-local diplotype inference and physical chaining

Status: **finite experimental implementation; not a production CLI**, 2026-09-15.

This document is the working plan for replacing genome-wide mosaic frontier search
with bounded local genotyping followed by physical haplotype chaining. It is meant
to keep the next implementation focused on the genomics algorithm rather than on
making the existing 1.5 kb prototype more elaborate.

The immediate target is an experimental, deterministic end-to-end recovery test.
It is not yet a production CLI design, a confidence model, a ploidy caller, or a
whole-yeast success claim.

## Decision

Do **not** continue trying to discover a 12 Mb mosaic by applying increasingly
clever complete-genome edits to a global incumbent.

Instead:

1. treat each existing syng catalog group (a homologous subpangenome) as a
   genotyping locus;
2. preserve every complete oriented source interval occurrence in that group as
   a physical allele, including coordinate and copy multiplicity;
3. genotype unordered pairs of those intervals from their summed variable-length
   maximal-MEM-subwalk count profiles;
4. retain all relevant physical alternatives, including count-equivalent ones;
5. chain local diplotypes through legal panel continuations, considering both
   homolog matchings at each boundary;
6. reconstruct two complete source-coordinate routes and sequences; and
7. rescore both complete sequences with the same variable-length MEM-subwalk
   profile and joint Poisson loss used locally.

The local model proposes genomes. It does not replace complete physical validation
or the global score.

## Relationship to LikeGT

The reference implementation is the canonical checkout at `/home/erikg/likegt`
(`https://github.com/pangenome/likegt`). Its core genotype kernel is
`src/commands/geno.rs::evaluate_genotypes`:

- represent each candidate haplotype by an evidence vector;
- enumerate combinations with replacement for the requested ploidy;
- add the selected haplotype vectors elementwise;
- score the combined vector against the sample; and
- rank the combinations.

We should reuse that simple exhaustive **local genotype construction**, not its
particular evidence model.

Our feature space is stronger than LikeGT's node-coverage vector. For every
retained maximal exact syng-anchor record, it contains every canonical odd-length
node-to-node contiguous subwalk: single-node, adjacent-node, intermediate, and
full-record patterns. Every occurrence position contributes, multiplied by the
maximal-record multiplicity. Direct candidate accumulation must equal
`WeightedBwt::count(key)`, including RC orbits, palindromes, and repeated patterns.
The features retain ordered graph context, exact inter-anchor spacing, orientation
orbit, observed zeros, unsupported positives, and junction changes. They must
never be collapsed to node depth or adjacent `[node,gap,node]` factors.

The differences are therefore substantive:

| LikeGT regional kernel | This design |
|---|---|
| graph path | complete oriented catalog source-interval occurrence |
| node coverage vector | sparse MEM-BWT occurrence-count profile |
| cosine similarity | existing background-relative Poisson loss |
| independent regional result | phase-aware state in a physical whole-genome chain |
| one path identity per vector | score-profile reuse may cover many distinct physical continuations |

The analogy ends at local combination enumeration. Physical chaining, boundary
MEM evidence, source-span capacity, exact route reconstruction, and global joint
rescoring are additional requirements here.

## Non-negotiable inference invariants

1. A diploid hypothesis contains exactly two independently public-valid complete
   `panel_routes::Assignment`s.
2. Counts from the two copies are added before exposure conversion and before the
   single logarithm:

   \[
   q_f=q_f(G_1)+q_f(G_2),\qquad
   R_f=s_f-C_f\log(1+s_f/\beta).
   \]

   The observation and background occur once per feature. Independently scored
   haploid losses are never added.
3. Each local allele is one complete oriented source interval from the catalog
   group, not a graph route or arbitrary fragment that merely touches a locus.
4. Source coordinates, orientation, entry and exit geometry, legal seams, and
   copy multiplicity remain explicit.
5. Canonical source-span conflicts remain illegal within one copy. Reusing the
   same source span across the two copies is legal dosage.
6. Equal count vectors do not make physical states interchangeable. States with
   different continuations, source-span ledgers, or reconstructed sequence remain
   distinct.
7. L150 linkage is used where the MEM-BWT retained it. Unlinked distant phase is
   preserved as ambiguity rather than invented.
8. No truth tract, donor hint, closest-reference choice, sample family, read
   placement, or assessment-only finite domain enters inference.
9. Complete paired rescoring uses the same fixed variable-length feature universe
   and arithmetic as local/transition scoring and remains authoritative.
10. Exhaustion, support completeness, and optimality are claimed only for a domain
    actually enumerated. A cap stop is not such a claim.

## Initial evidence universe

The implementation uses canonical variable-length subwalk keys, sample counts
queried from `WeightedBwt`, the L150 histogram, depth, background, and existing
Poisson arithmetic. For one read length, `denominator=L*H`, exposure is `q*H /
denominator`, and signal is therefore `depth*q/L`; the histogram must not remain in
the final scale. The stride-1 finite control uses depth 150, while the production-like
stride-15 route control uses depth 10. It does not create a surrogate based on base
mismatches, node presence, donor labels, or read placement.

The universe is frozen from all public catalog-allele interior profiles and every
publicly legal adjacent-seam profile before genotyping. Observations are point
queries against the supplied aggregate `WeightedBwt`; raw reads and truth-derived
profiles never enter inference or feature-cap accounting. Unsupported candidate
positives remain in the score. Longer sample-only keys cannot be enumerated from
the BWT and are reported as an audit limitation; numerically their candidate-zero
background-relative contribution is zero. `max_features` fails closed with no
truncation, and truth never selects feature length or content.

Singleton-node maximal MEM/subwalks are valid one-token features and must not be
noise-filtered. A globally unique singleton can provide strong local
presence/dosage evidence, but it provides no junction or phase linkage and one
observation is not treated as definitive. The experimental audit reports exact
bins by node count/token length (singleton, pair, and three-or-more-node), summed
observed and predicted multiplicity, and empirical local-versus-shared incidence.

### Production projection versus retained information

This distinction is critical. L150 and the weighted MEM-BWT are not the limitation
exposed by the central phase counterexample: `WeightedBwt::count` can count any
valid odd-length node-to-node pattern. The current production route registry and
scorer deliberately materialize only `[node,gap,node]` triples. That downstream
projection can make two hypotheses equal after the BWT has retained a distinguishing
three-node/five-token (or longer) pattern.

The isolated experimental module proves the minimal extension seam without
rewriting production scoring:

1. use `Vec<u64>` as the registered feature key;
2. register at least all public candidate-derived five-token patterns, and in the
   experiment every node-to-node subwalk through the full maximal record;
3. replay complete ranged starts and complete reconstructed sequences through
   `mem_records::canonical_mem_records`;
4. obtain each observation with `WeightedBwt::count(key)`;
5. dynamically retain candidate-positive/sample-zero terms rather than dropping
   unsupported positives; and
6. add both copy counts before the existing exposure conversion and single
   background-relative Poisson logarithm.

Production `panel_routes` remains triple-projected in this slice. Moving this
`Vec<u64>` registry/replay seam into the authoritative whole-route scorer requires
an explicit reviewed follow-up; a local surrogate must not claim phase that the
production global score cannot yet distinguish.

Every registered factor is classified exactly once for proposal scoring:

### 1. Partition-owned emission factor

A factor may be a local emission only when the existing physical ownership and
full read-context checks establish that its prediction is determined entirely by
one partition allele. This includes zero-count contexts. For a local diplotype
`{a,b}`, compute its integer count as

\[
q_f(a,b)=q_f(a)+q_f(b)
\]

and apply one existing loss term using the sample's single `C_f`.

### 2. Adjacent-boundary transition factor

A factor may be transition evidence only when its complete support envelope is
determined by the joined choices on exactly two adjacent partitions. Its predicted
count is computed for the actual oriented seam and for each homolog matching.
These factors distinguish, for example, `A-A + B-B` from `A-B + B-A` when an L150
MEM context spans the boundary.

A feature that spans more than one boundary, has source-terminal uncertainty, or
cannot be attributed to one particular adjacency is not forced into this class.

### 3. Deferred global factor

Shared, repeated, nonlocal, multi-boundary, or otherwise context-dependent factors
are omitted from local/transition scores. They are evaluated only by the exact
whole-genome scorer after reconstruction.

Deferred does not mean absent or zero. It means the local proposal score is
partial. The final global score can reorder locally preferred chains.

### Accounting audit

Maintain a feature-ID classification ledger. The owned and boundary sets must be
disjoint, and their union must be disjoint from the deferred set. A factor may
contribute at most one background/observation term to the local decomposition.
The ledger must report why each deferred factor was not localized.

## Candidate and state model

The initial Rust types may evolve, but the implementation must preserve these
logical fields.

### Partition allele

```text
PartitionAllele
    partition/group ID
    catalog occurrence ID
    source ID and half-open source coordinates
    explicit orientation (both retained when the catalog strand is unknown)
    exact oriented interval sequence
    sparse variable-length MEM-subwalk profile for wholly interior L150 starts
```

Two alleles with identical sequence or counts may share cached profile storage.
They remain separate physical alleles if their occurrence, source coordinates,
orientation, multiplicity, or legal continuations differ. Routes are reconstructed
from chosen source intervals after chaining; they do not define the partition or
its alleles.

### Local diplotype

```text
LocalDiplotype
    unordered allele IDs [min(a,b), max(a,b)]
    summed sparse integer owned-factor profile
    local emission loss
    two physical entry keys
    two physical exit keys
```

Homozygous combinations are included. Enumeration is the diploid
`combinations_with_replacement(2)` pattern used by LikeGT.

The diplotype is unordered for local genotyping. It becomes temporarily ordered
when a transition assigns its two alleles to the two growing homologs.

### Chaining state

```text
ChainState
    ordered route prefix for copy 0
    ordered route prefix for copy 1
    current exit key for each copy
    within-copy canonical source-span ledger for each copy
    accumulated owned-emission and boundary-transition score
    exact local count/signature information needed for ties and audit
    backpointer or bounded reconstruction record
```

At each boundary evaluate both matchings:

```text
previous copy 0 -> next allele 0, previous copy 1 -> next allele 1
previous copy 0 -> next allele 1, previous copy 1 -> next allele 0
```

Reject a matching if either join is physically illegal or if extending either
copy violates within-copy source-span capacity. Copy labels are bookkeeping; final
pairs are normalized only after reconstruction.

A DP key must include the sufficient physical continuation and resource state.
It is forbidden to merge states solely because they have the same local count
vector or score. Profile equivalence may cache arithmetic, never discard distinct
continuations.

## Allele enumeration

The first implementation reuses the existing public catalog rather than treating a
graph route as a locus.

1. Read each catalog group and every homologous source interval occurrence in it.
2. Preserve the exact half-open coordinates, orientation, occurrence identity,
   and multiplicity. If strand is unspecified, retain both orientations.
3. Fetch and orient that exact sequence subrange; reject invalid or incomplete
   source bindings.
4. Profile every complete L150 start wholly within the subrange directly through
   `mem_records::canonical_mem_records` and variable-length subwalk accumulation.
5. Enumerate public physically compatible adjacent interval pairs and profile only
   complete L150 windows crossing their concatenated seam.
6. Enforce this slice's decomposition precondition: every typed interval has
   length at least L150, so one complete start crosses at most one boundary.
7. Cache profile arithmetic without deduplicating physical occurrences.

The initial deterministic fixture may use a declared linear partition axis. A
whole-genome extension must obtain adjacency from public axis/route/port structure,
not from a chosen closest sample. Rearrangements or alternate topology that cannot
be represented by the initial linear chain remain explicitly unsupported rather
than silently coerced.

## Local diplotype scoring

For each partition with `n` physical alleles, enumerate `n(n+1)/2` unordered pairs.
Do not call the complete whole-genome scorer for every local pair.

1. Obtain each allele's sparse integer profile once.
2. Add the two profiles with checked integer arithmetic.
3. Evaluate all partition-owned factors, including observed positive factors with
   zero prediction and predicted positive factors with zero observation.
4. Apply the same histogram exposure conversion and `beta` as the exact scorer.
5. Store the local loss, exact count signature, and physical allele pair.
6. Retain all exact count-equivalent physical diplotypes needed by distinct
   continuation/resource states.

Many of the 235 panel identities may have the same local allele profile. It is
valid and important to evaluate an identical profile-pair loss once and reuse that
number. Every physical allele-pair member must still be available to chaining.
This is score caching, not closest-reference selection and not physical-state
pruning.

Floating loss equality is not used as a substitute for integer count equality.
Use exact sparse count signatures for count equivalence and the existing numerical
tolerance only when comparing computed objectives.

## Boundary scoring and phase

A transition combines physical validity with the boundary factor subset.

For each predecessor state, next local diplotype, and the two homolog matchings:

1. join the actual oriented source intervals for each copy;
2. compute or retrieve the maximal-MEM-subwalk profile of complete L150 windows
   crossing that concatenated seam;
3. add counts across both copies;
4. score every factor assigned to that boundary once;
5. update each copy's physical source-span ledger; and
6. retain all best/tied states with distinct sufficient continuation signatures.

Boundary evidence is the mechanism for local phase. If the two matchings produce
different count vectors, the better supported matching must be recoverable. If
they produce the same vector and remain physically valid, both are legitimate
alternatives.

Do not add a switch penalty, donor persistence prior, family prior, or arbitrary
phase preference. Physical compatibility and observed MEM counts are the only
transition criteria in this slice.

## Reconstruction and authoritative scoring

A completed chain is only a proposal until it passes all of the following:

1. concatenate the selected oriented interval sequences for each copy;
2. retain their exact source-coordinate route descriptions and validate every
   physical join and within-copy source-span constraint;
3. replay all complete L150 starts on each reconstructed sequence through the same
   maximal-MEM/subwalk profiler;
4. add integer counts across copies before one Poisson loss over the fixed feature
   universe; and
5. retain global-score ties and exact count-equivalent alternatives within the
   declared state limits.

Report both the decomposed proposal score and exact global score. They are expected
to differ whenever deferred factors exist. Selection of final hypotheses is by the
exact global score, not by the local DP score.

## First deterministic red/green experiment

Build the smallest fixture that exercises local dosage, boundary phase, and full
route reconstruction together.

### Required shape

- Four ordered partitions.
- At least two complete spanning alleles, A and B, in every partition.
- Native whole routes A-A-A-A and B-B-B-B.
- Diploid truth has local B dosage `[0,1,1,2]`.
- The two truth copies are both non-native mosaics.
- The central two heterozygous partitions admit both local phase matchings.
- An L150 boundary context makes one central matching distinguishable in the
  linked control.

A suitable abstract truth is:

```text
partition       1  2  3  4
copy 0          A  A  B  B
copy 1          A  B  A  B
B dosage        0  1  1  2
```

The alternative central phase is:

```text
copy 0          A  A  A  B
copy 1          A  B  B  B
```

The inference input contains only the panel, routes, sample MEM-BWT, partition
structure, and ordinary options. Truth paths and tract coordinates exist only in
assessment code.

### Controls

1. **Linked positive:** boundary-spanning MEM counts differ between the two phase
   matchings. Every retained optimum must recover the supported matching, exact
   dosage, and two complete non-native mosaics.
2. **Unlinked ambiguity:** remove the informative boundary linkage without
   changing the local diplotypes. Both count-equivalent phase classes must remain
   represented if both are constructed within an exhaustive finite run.
3. **Joint-count negative:** independently summing haploid losses must disagree
   with the accepted joint formula on at least one heterozygous state.
4. **Physical negative:** a locally attractive allele sequence with an illegal
   continuation or within-copy reused source span must be rejected.
5. **Deferred-factor negative:** a nonlocal/shared factor must not be charged in
   two partition emissions or silently assigned zero.

### Finite oracle

Assessment enumerates every legal whole haplotype in this small panel and every
unordered diploid pair. It independently computes complete occurrence-count
vectors and the joint objective. This finite oracle is assessment-only.

The green test requires:

- exact agreement with the finite minimum within the existing tolerance;
- exact truth count-vector recovery in the linked control;
- dosage `[0,1,1,2]`;
- two public-valid non-native complete assignments;
- correct linked phase;
- retained count-equivalent alternatives in the unlinked exhaustive control;
- no factor-ID overlap among local, boundary, and deferred ledgers; and
- exact whole-genome rescoring agreement with an independent formula.

## Implementation order

### Milestone 0 — inventory and classification audit

Before implementing DP, write a small read-only probe over the fixture that emits:

- partitions and their physical source intervals;
- candidate spanning traversals and entry/exit keys;
- factor classification as local, one-boundary, or deferred; and
- allele/profile equivalence-class sizes.

Stop and correct the representation if any expected local or linked-boundary
factor cannot be classified without truth knowledge.

### Milestone 1 — red fixture and exhaustive oracle

Add the four-partition fixture, independent finite oracle, and failing test before
adding the partition search. The expected failure is absence of the new search
result, not an artificially weakened biological assertion.

### Milestone 2 — allele profiles and local exhaustive genotyping

Implement spanning-allele extraction, sparse profile caching, combinations with
replacement, joint local loss, exact count signatures, and physical alternative
retention. Verify every partition's expected genotype against its exhaustive local
oracle before adding phase chaining.

### Milestone 3 — two-matching physical DP

Implement legal joins, per-copy resource ledgers, both phase matchings, boundary
factor scoring, and backtracking. Start exhaustive on the tiny fixture. Do not add
a beam until the exact state behavior is understood.

### Milestone 4 — complete reconstruction and global rescore

Produce two complete assignments, validate them publicly, rescore their complete
sequences with the experimental variable-length universe and joint formula, and
compare against the separately named finite whole-pair assessment oracle. The
production triple scorer remains unchanged and is insufficient for the long-linkage
witness. The old global frontier
may be run only as a baseline comparison; it is not extended to make the test pass.

### Milestone 5 — scaling ladder

After the finite test passes and receives review, test approximately 10 kb, 100 kb,
and 1 Mb deterministic mosaics before retrying whole yeast. At every size report:

- partitions, physical alleles, and unique allele profiles per partition;
- physical and unique-profile diplotype counts;
- transition candidates and legal matchings;
- retained physical states and equivalence classes;
- owned, boundary, and deferred feature counts;
- profile work, local score work, exact global scores, validations, wall time,
  peak RSS, and logical state;
- exact sequence/dosage/phase recovery and retained ambiguity; and
- whether search/support was exhaustive or capped.

Only optimize the measured dominant cost. Likely safe optimizations include sparse
profile addition, exact profile-pair score caching, shared boundary-profile caching,
and backpointer storage. They must not erase physical continuations.

## Resource discipline

The reviewed finite prototype keeps its original ceilings unchanged: 4M work,
8M profile work, 2,048 complete paired scores, 4,096 public haploid scores,
128 MiB logical state, 50,000 features per materialized structure, and 64
segments. The separate genome experiment has a documented 600,000,000
MEM-query-run ceiling. This corrects the earlier prototype-scale mistake of
applying the finite 8M meter globally: one exhaustive public pass is estimated at
about 500M event runs, while comparable existing machinery processed 2.17B events
in about 70 minutes. Genome accounting reports each locus, each boundary, and the
global total; exceeding 600M stops incomplete. It does not alter production or
the finite module.

Use at most four threads, CPUs 252–255 at nice 10, Cargo
`--offline --locked --release -j4`, and serialized tests. A bounded stop reports
incomplete search. It never becomes evidence of infeasibility or exhaustive
support.

## Genome-scale experimental path

`examples/panel_route_diploid_search/genome.rs` loads the public reference axis
and BED3 catalog groups, preserves occurrence, source coordinates, multiplicity,
and orientation, and permits contiguous same-source traversals plus cross-source
seams whose verified production exit/entry port words match exactly. A bounded,
truth-free refinement can also replace one forward interval with a two-segment
spanning traversal: an occurrence-local prefix joined to another occurrence's
suffix at equal interior production port words. Blocks are expanded from reported
owned-feature margins, with selected classes retained in the top-K scope. Cuts are
ranked by the complete composed traversal profile (prefix interior + suffix
interior + internal seam) jointly with the locus's public native allele, retaining
top-N per ordered occurrence pair plus every exact best-loss tie. The physical
beam is configurable (`--beam-width`, default 1024) and remains subject to the
128 MiB logical-state cap. Each boundary is filled globally into one bounded
physical-ledger layer; there are no per-state successor quotas or local-loss early
stops once cumulative scores are in use. Physical boundary links sharing exact
oriented L149 endpoint sequences reference
one composed seam profile, so expansion never profiles physical member pairs.
Either truncation makes the run explicitly `bounded_refined`, never exact. A
configurable RSS guard defaults to 16 GiB and emits stage/allocation diagnostics
before aborting. Interior and seam profiles use event-compressed
complete L150 replay. One full canonical
maximal-MEM query is made per invariant run and every odd-length node-to-node
subwalk contribution is multiplied by the run length. Sealed append-only caches
are keyed by physical occurrence and orientation. A parity regression compares
event compression against exhaustive per-window replay.

Feature incidence is classified in deterministic shards, so no materialized
incidence or scoring structure exceeds 50,000 keys. The bounded merger carries a
persistent cumulative joint-count profile as an `Arc` parent plus an interned
sparse locus/seam delta. Every transition applies the exact Poisson difference
`loss(q_old + delta) - loss(q_old)`, including shared and Deferred features.
Identical cumulative profiles with identical conflict-window suffixes merge only
for search, after exact profile equality verification; physical histories remain
distinct tie backpointers. Logical-state charging counts each reachable shared
node and delta once rather than charging duplicated complete maps.

The merger is exact min-plus DP over ordered physical exit pairs, with arithmetic
cached by exact interior-profile and seam-signature classes. Boundary arithmetic
adds both copy counts before one Poisson loss,
and both homolog matchings arise from the two ordered continuations. Backpointers
retain equal-score predecessor sets; a 2,048-alternative overflow is explicit
`complete=false`, never silent pruning. Before chaining, the audit reports each
exit-signature `E_i`, `E_i^2 E_{i+1}^2`, and the maximum locus distance `W` of
candidate same-source overlap. `W=0` makes the current physical exit pair a
sufficient resource state. For nonzero `W`, the experimental path retains exact
in-chain overlap eligibility in a conflict-window history and applies the configured
beam only when required; it always reports bounded/incomplete and counts both
dropped cutoff ties and dropped non-ties. Internally overlapping two-segment
traversals are rejected before admission to the initial or successor layer.

The incidence pre-pass stores only two hashes and compact location state per key.
A detected primary-hash collision is wholly deferred. Profiles are loaded one
locus at a time from an offset-indexed sealed JSONL cache, then discarded. Sample
counts are queried lazily and retained in 256 bounded hash shards. Complete-copy
rescoring emits at most 50,000-key sorted runs, merges them with fan-in 16, sums
the two copies per key, and performs one final Poisson loss.

The chrI exhaustive profile preflight completed all 8,944 oriented alleles in
11,172,346 MEM-query runs (80,158,038 integrated windows), only 1.40% above the
previous 8M early-stop average per allele and therefore below the owner’s 2x
review threshold. The occurrence-keyed chrI audit measured 11,225,933 executed
runs including 2,871 seams. Its exit-class counts were
`[2,2,3,22,30,26,49,30,24,23,25,52,31,30,17,3,10,15,34,2,51]`, the largest
`E_i^2 E_{i+1}^2` was 2,598,544, and `W=0`. The eight-partition 80 kb control now
runs through the exact DP and completes without beam pruning. The streaming chrI
gate initially completed in 894.38 seconds with 1,794,012 KiB peak RSS. Preparing
owned profile classes once and scoring their Cartesian products on three Rayon
workers reduced the same exact run to 59.60 seconds with 1,904,004 KiB peak RSS:
input 5.95 s, cache index 3.41 s, streamed interior incidence 12.70 s, seams
0.06 s, exact DP 26.74 s, and external final rescore 10.65 s. The optimized run
retained one exact physical pair, used
721,174 feature hashes with 20,234 deferred and zero detected collisions, and
rescored 860 bounded runs / 61,892 MEM-query runs to
`-1286625.6161305453` for chrI. Both reconstructed routes are the complete public
S288C chrI source interval and pass coordinate-seam/span validation. A
current-source route-artifact rebuild completed. Current-source evaluation of the
unchanged S288C family-216 assignment is `-10525297.667494176`; the historical
artifact objective is retired because that artifact records an unreviewed WIP
compiler identity. Lanes 0--3 and S288C lanes 9564--9567 have byte-identical
native totals across the artifacts, and registry/ownership seals are also
identical, localizing the change to current evaluator semantics rather than
native-profile compilation.

The resumed chromosome-reset whole-axis run completed all 17 components in
4:00:06 with 3,508,952 KiB peak RSS and 538,106,401 new profile-query runs, below
the 600M ceiling. Ten `W=0` components used exact min-plus DP. Seven components
(chrIII, chrIV, chrX, chrXII, chrXIII, chrXV, chrXVI) used the approved bounded
conflict-window policy (`W=2..82`) and therefore remain explicitly incomplete,
even though completion-viability filtering made every dropped-tie and dropped
non-tie counter zero in this run. One external sorted-run merge over all 17
chromosomes and both copies produced full-subwalk loss `-34233513.67185511`.
Matched native and assessment-only truth duplicate-pair rescoring on the rebuilt
artifact produced `-34260135.38842692` and `-34133457.58802463`, respectively.
The chain did not recover the chrIII SK1 tract: both selected chrIII copies are
native S288C source 9564 over `[0,341580)`. Moreover, only chain copy 1 passes the
authoritative production Evaluator; copy 0 fails native endpoint pairing because
several full public-axis donor routes are not legal production route assignments.
The external score is thus an experimental physical-copy score, not a validated
production pair objective.

## The sufficiency principle (logged 2026-09, from the owner's formulation)

The controlling design rule for all search and scoring layers, distilled from the
campaign's measured failures and fixes:

**We never need to score all possible pairs — only a sufficient set. Two
candidates are interchangeable until some evidence can distinguish them, and the
MEM structure itself defines where evidence exists at all.**

The evidence is symmetric MEM-finding: haplotype-haplotype MEMs are what the
pangenome graph is *built from* (syng), and read-haplotype MEMs are what the
sample index stores (`sample.membwt`). They are the same object class. Their
intersection with a decision defines the sufficient scoring set: a feature can
discriminate a pair only if (a) the candidates' predicted counts differ on it
AND (b) it has observations. Features with no observations, or identical
predictions across candidates, cancel in every comparison — scoring them is
pure waste. Candidates with identical profiles are one class — enumerating
their members separately is pure waste. Pairs provably unable to beat the
incumbent are decided — expanding them is pure waste.

The campaign's three measured collapses were this one principle applied at
three layers, each with its measured win:

1. **Feature sufficiency** (the telescoped extend): features with q_old = 0
   contribute per-feature *constants* — per-state arithmetic is needed only on
   the overlap set (features the state has touched). Measured: 240 µs/extend →
   40–60 ns per delta feature.
2. **Candidate/class sufficiency** (byte-identical class dedup; top-N +
   exact-best-ties admission; tie-sets retained as sets, never physically
   enumerated): measured 38M physical seams → 1,204 endpoint classes; 4.88B
   transition instances → 0.67M folded edges.
3. **Decision sufficiency** (seed incumbent + admissible bounds): check pairs
   only until optimality is *proven* — the search's job is proof, not
   discovery (the slice's greedy incumbent was already optimal to 5e-9).

The corresponding **failure mode**, also measured three times: any new path
that re-instantiates search in unfolded form (per-state × per-feature ×
per-physical-pair arithmetic) reintroduces the explosion — the composition
explosion, the 5.3B-transition DP, and the routed-DP c14 run (492M dropped
non-ties, 42 min) were all the same disease arriving through fresh code.

**Structural obligation:** the folded search must be a shared engine that
consumes (evidence, candidates) as inputs. No integration path may instantiate
its own DP over physical pairs. Any new evidence layer (including the MEM-routed
collection) must feed the shared folded engine, never re-implement scoring.

## The block model (logged 2026-09, from the owner's formulation; companion to the sufficiency principle)

The evidence architecture, stated in one object hierarchy:

1. **Blocks** are the maximal shared segments between haplotypes — which is
   what the pangenome graph already encodes. The panel's compression IS the
   block set; a "deeply covering set of MEMs or blocks" over the haplotype set
   is the graph's own segment/occurrence structure. No new structure to
   invent: the territory index and geometric-q derivation already consume
   panel segments this way, and candidates are paths over blocks.
2. **Read MEMs against blocks** are what `sample.membwt` already stores
   (read-vs-panel maximal MEM records, canonical subwalks). The missing base
   data type is a first-class **read-to-block index**: every read mapped to
   the blocks it matches, with offsets and orientations, built ONCE per
   sample. It replaces the per-run sample-record re-derivation (the measured
   5.8–13 s API gap), compresses evidence finding to block-local lookups,
   and — because a read spanning two blocks carries their adjacency with
   offsets — it also stores the within-read MEM adjacency that the
   evidence-topology scout independently proved is required for junction
   detection. One base type fixes three known problems: the API gap, evidence
   compression, and junction observability.
3. **Candidates are paths over blocks**; predicted counts derive from geometry
   (Gate 1 validated: geometric q vs the profile_read oracle, every residual
   diff classified, decision-level effect <= 2.1% with direction unchanged).
4. **The search runs on block-level sufficiency** through the shared folded
   engine (the structural obligation above): per-block observed support,
   per-class-pair constants, O(overlap) extends, proof-by-incumbent. No path
   may re-instantiate unfolded search.

**The one open design choice**: block granularity — deep coverage (few large
blocks) maximizes compression but loses variant resolution at boundaries;
shallow (many small blocks) the reverse. It is a measurable tradeoff (score
resolution vs block count), not an open research problem; the k63 anchor
length sets the floor and the current BED/segment structure is one point in
the space.

**Status (2026-09)**: pieces 1, 3, 4 exist and are validated (routed slice
25.2 s zero-cache; routing 0.4 s; geometric profiles 4.5 s; exit-cut decision
preserved in all model variants). Piece 2 (the read-to-block index as a base
type) is queued. The routed-DP port of the folded engine is in flight; the
42-minute c14 run was the unfolded-search failure mode, not a property of the
block model.

## Explicit non-goals for this slice

- No production diploid CLI.
- No raw-read placement or reconstruction of discarded read identities.
- No closest-reference or sample-family selection.
- No CNV/ploidy inference.
- No switch penalty or copying prior.
- No calibrated posterior/confidence.
- No forced resolution of L150-unlinked phase.
- No conversion to PGGB or graph smoothing prerequisite.
- No claim that local proposal scores fully decompose the global repeat/shared
  likelihood.
- No more blind-challenge or sealing infrastructure.

## Stop conditions

Pause implementation and diagnose rather than adding complexity if:

1. a local allele cannot be defined as a complete physical spanning traversal;
2. factor ownership would count one observation/background more than once;
3. a boundary score requires hidden read placement or truth coordinates;
4. state deduplication would merge different physical continuations;
5. reconstruction cannot yield two independently valid complete assignments;
6. the linked finite oracle prefers truth but the DP does not construct it; or
7. the exact global scorer reverses the local result because an allegedly local
   factor was misclassified.

The intended next result is narrow and concrete: a reviewed deterministic test in
which LikeGT-style exhaustive local diplotype genotyping over the stronger MEM-BWT
feature space, followed by phase-aware physical chaining, recovers a multi-partition
mosaic that complete-genome frontier search did not efficiently discover.
