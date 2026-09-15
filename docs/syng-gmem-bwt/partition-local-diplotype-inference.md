# Partition-local diplotype inference and physical chaining

Status: **implementation guide; not yet implemented or validated**, 2026-09-15.

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

1. treat each bounded syng partition as a genotyping locus;
2. enumerate complete physical traversals spanning that partition;
3. genotype unordered pairs of those traversals from their summed MEM-BWT count
   profiles;
4. retain all relevant physical alternatives, including count-equivalent ones;
5. chain local diplotypes through legal panel continuations, considering both
   homolog matchings at each boundary;
6. reconstruct two complete `panel_routes::Assignment`s; and
7. use the unchanged exact whole-genome paired loss as the authoritative final
   score.

The local model proposes genomes. It does not replace complete-route validation or
the global score.

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

Our feature space is stronger than LikeGT's node-coverage vector. The current
route scorer queries canonical node-gap-node substrings from a weighted MEM-BWT
whose records originate as maximal exact syng-anchor subwalks. The scored features
retain ordered graph context, exact inter-anchor spacing, orientation orbit,
occurrence multiplicity, observed zeros, unsupported positives, and junction or
terminal changes induced by a physical route. They must never be collapsed to
plain node depth.

The differences are therefore substantive:

| LikeGT regional kernel | This design |
|---|---|
| graph path | complete physical partition-spanning traversal |
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
3. Each local allele is a complete spanning traversal, not an arbitrary fragment
   that merely touches a partition.
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
9. The existing complete paired scorer remains unchanged and authoritative.
10. Exhaustion, support completeness, and optimality are claimed only for a domain
    actually enumerated. A cap stop is not such a claim.

## Initial evidence universe

The first implementation must use the same registered feature tokens, sample
counts, read-length histogram, depth, background, and arithmetic as the public
route evaluator. It must not create a separate surrogate based on base mismatches,
node presence, donor labels, or read placement.

The presently registered route factors are canonical `[node, gap, node]` token
triples. They are exact ordered-and-spaced substring queries against the weighted
MEM-BWT, not node-depth observations. Supporting longer registered subwalk factors
later is compatible with this architecture, but is not required for the first
slice.

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
    stable physical allele ID
    complete oriented Route through the partition
    entry continuation/port key
    exit continuation/port key
    canonical source spans consumed by this allele
    sparse integer count profile for owned factors
```

Two alleles with identical sequence or counts may share cached profile storage.
They remain separate physical alleles if their ports, source coordinates,
orientation, or legal continuations differ.

For the first linear fixture, an allele may be one source segment. The abstraction
must nevertheless be a `Route`, so a later internally branching partition can use
a multi-segment spanning traversal without changing the diplotype interface.

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

The first implementation should reuse existing public panel data rather than build
a second graph model.

1. Read fixed partition ownership/group intervals and public route/port data.
2. For every physical source occurrence capable of entering and exiting the
   partition, construct the oriented spanning `Route` and its boundary keys.
3. Normalize harmless adjacent storage pieces without changing physical copy
   count.
4. Reject nonspanning fragments, ambiguous/incomplete boundary geometry, invalid
   orientation seams, and paths that fail the existing public validation rules.
5. Deduplicate only exact duplicate physical descriptions.
6. Cache allele profiles by exact route and, separately, cache score arithmetic by
   exact count vector. Do not deduplicate the physical alleles themselves by
   profile.

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

1. join the actual oriented route pieces for each copy;
2. compute or retrieve the exact boundary-seam count profile;
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

1. coalesce exactly contiguous same-source/same-orientation segments;
2. build exactly two complete `panel_routes::Assignment`s;
3. run `Evaluator::validate_assignment` independently on both copies;
4. evaluate the pair through the unchanged joint scorer, adding integer route
   counts across copies before one loss; and
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

Produce two complete assignments, validate them publicly, run the unchanged paired
scorer, and compare against the finite whole-pair oracle. The old global frontier
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

The existing paired-search ceilings remain unchanged: 4M work, 8M profile work,
2,048 complete paired scores, 4,096 public haploid scores, 128 MiB logical state,
50,000 features, and 64 segments, together with the other existing host limits.
Local partial-score operations require explicit accounting; they are not hidden
complete scores. Exact global rescoring consumes the existing paired/public-score
budgets normally.

Use at most four threads, CPUs 252–255 at nice 10, Cargo
`--offline --locked --release -j4`, and serialized tests. A bounded stop reports
incomplete search. It never becomes evidence of infeasibility or exhaustive
support.

If whole-yeast local combination counts exceed the work ceiling, first exploit
exact profile equivalence for arithmetic reuse while retaining all physical
members. Do not respond by choosing a nearest donor, dropping count-equivalent
continuations, raising the cap, or weakening the feature space.

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
