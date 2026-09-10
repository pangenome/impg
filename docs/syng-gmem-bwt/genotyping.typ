#set document(
  title: "Local genotyping with a syng-gMEM-BWT",
  description: "Proposed sample evidence index and local haplotype-mixture model for impg.",
)
#set page(paper: "a4", margin: (x: 22mm, y: 20mm), numbering: "1")
#set text(size: 10.5pt)
#set par(justify: true)
#set heading(numbering: "1.")
#show heading.where(level: 1): set block(above: 1.3em, below: 0.6em)
#show link: set text(fill: rgb("245886"))
#show raw.where(block: true): set block(breakable: false)

#let occ = math.op("occ")
#let select = math.op("select")
#let rc = math.op("rc")
#let argmax = math.op("arg max", limits: true)

#align(center)[
  #text(size: 21pt, weight: "bold")[Local genotyping with a syng-gMEM-BWT]
  #v(0.5em)
  #text(size: 12pt)[A compact sample index for contextual haplotype coverage]
  #v(0.5em)
  Proposed design for impg
]

#block(fill: rgb("f0f4f7"), inset: 10pt, width: 100%)[
  *Status: specification, not an implemented genotyper.*
  impg already provides a syng panel index, anchor mapping, GBWT MEM queries,
  node-coverage genotyping, and experimental mosaic inference. The sample
  gMEM-BWT and contextual scorer below are proposed additions. The design stays
  entirely in syng; explicit variation graphs and bluntification are optional.
]

= Objective

Map a sample once to the syng pangenome, collect its graph maximal exact matches
(gMEMs), and index their strings in a collection BWT. Query this sample index
with candidate haplotypes at any locus to obtain contextual matching counts.
Infer the local genotype as a mixture of haplotype copy-number vectors.

```text
reads -> syng gMEM collection -> sample gMEM-BWT
                                      ^
                                      | matching-count queries
panel syng GBWT -> locus haplotypes ---+
                                      |
                         local haplotype-mixture scoring
                                      |
                         optional chromosome mosaics
```

The durable sample object is a compact collection of observed matched strings,
not a per-read alignment database. Read names, original alignments, individual
panel placements, and a read-by-haplotype matrix are not required. Store the
panel fingerprint, syncmer parameters, counting convention, and aggregate
read-length/error metadata needed to interpret counts. Retaining raw reads
elsewhere permits remapping, but is not a query dependency.

= Panel and sample representations

== Syng haplotype walks

Represent a haplotype as

$ H = (v_1, delta_1, v_2, dots, delta_(n-1), v_n), $

where $v_i$ is an oriented syncmer node and $delta_i$ is the distance between
consecutive syncmer starts. A contiguous subwalk includes both node identities
and these distances. With syncmer length $k$, its sequence span is

$ s_f = k + sum_(i=1)^(m-1) delta_i, $

for a feature containing $m$ nodes. These are matches in the *anchor-and-spacing
representation*. Matching spacing across an uncovered gap does not prove that
all intervening nucleotides match.

== Sample gMEM collection

Let the mapped sample yield a multiset

$ cal(M) = {(M_t, w_t)}_(t=1)^T, $

where $M_t$ is a gMEM string and $w_t$ is its multiplicity. Identical strings may
be collapsed. Maximality is defined relative to the mapped query trace and the
panel GBWT, under the mapper's node/spacing rules. It need not remain maximal
relative to a particular candidate haplotype. Indexing all substrings of each
stored gMEM preserves access to its shorter compatible contexts.

Build a collection BWT/FM-index over a lossless tokenization of these strings.
Separators prevent matches across gMEM boundaries. Search patterns begin and end
on node tokens; spacing tokens belong to their corresponding internal
transitions. No nucleotide graph needs to be materialized.

A full collection BWT counts occurrences by suffix-interval width. If duplicate
strings are collapsed, the index must instead support a *weighted interval sum*,
with each suffix carrying its string's multiplicity. Weighting only complete
MEM records, without a corresponding substring-count mechanism, is insufficient.
Run-length compression and weighted document storage are implementation choices;
compression ratios and query costs must be measured.

Define the sample count primitive as

$ C(f) = sum_(t=1)^T w_t occ(f, M_t). $

Here $occ$ counts contiguous occurrences, identifying a feature with its reverse
complement. Query both orientations, counting a reverse-complement-palindromic
feature only once. Repeated occurrences at distinct positions remain distinct.
Do not concatenate adjacent or overlapping gMEMs into an unobserved string.

The unit is a *gMEM-substring occurrence*, not a distinct-read count. Overlapping
MEMs from one query may contribute multiple observations. Freeze the extraction
and overlap convention, and reproduce it when estimating exposure and validating
counts. This design deliberately does not preserve molecule identity or linkage
between separate gMEMs.

= Locus candidates and contextual features

== Partitions are not genotype or recombination boundaries

Computational partitions organize source-coordinate occurrences, not disjoint
sets of graph nodes. Refine them into homologous inference chunks with validated
spanning candidates, orientation, boundaries and source-continuation links.
Seed retention establishes source accounting, not homology or callable evidence.
Recombination can occur inside these units. Keep core ownership separate from
overlapping matching context, and count each predicted genomic occurrence once.

Many panel homologs are alternative candidates, not sample copy number. Tandem
repeat multiplicity belongs inside spanning candidate walks; dispersed copies
retain their separate occurrence identities. Candidate equivalence classes must
preserve source members and boundary compatibility, not just local feature
vectors. The companion `partition-inference-contract.md` records the current
yeast evidence, bounded scheduling repair and implementation gates.

== Candidate genotypes

For a locus $W$, retrieve homologous panel subwalks

$ cal(H)_W = {H_1, dots, H_K}. $

Use source-coordinate homology from impg, not identical numerical intervals on
each haplotype. Candidates may differ in length, orientation, and repeat count.
A diploid genotype is an unordered pair with repetition:

$ G = (H_a, H_b), quad 1 <= a <= b <= K. $

There are $K(K+1)/2$ such pairs, including homozygotes. For ploidy $p$, use
multisets of $p$ haplotypes. A haplotype's internal repeat multiplicity is
separate from the number of chromosome copies in the genotype.

== A common feature universe

A feature is a contiguous oriented syng subwalk. Select a shared feature set
from the candidate walks:

$ cal(F)_W = select(union_(h=1)^K op("subwalks")(H_h)). $

A practical first policy samples a bounded set of context lengths and start
positions. Adaptive extension can distinguish candidates or separate a locus
from off-locus repeats. Freeze selection rules across the compared genotypes;
avoid searching all possible substrings, and retain zero-count sample features.
Record selected features and weights so scores can be reproduced.

For each $f in cal(F)_W$, query the sample count and candidate multiplicity:

$ x_f = C(f), quad a_(h f) = occ(f, H_h). $

Under a diploid genotype, expected feature copy number is

$ g_(G f) = a_(a f) + a_(b f). $

One-node features measure node support; two-node features add orientation and
spacing; longer features measure haplotype context. A haplotype matching profile
can be written as

$ Q_h(i, ell) = C(H_h[i:i+ell]), $

where the slice denotes $ell$ consecutive nodes and their internal spacings.
This profile is a query interface, not a license to independently score each
haplotype in a different feature space. Compare all genotype mixtures against
the same $x$ and $cal(F)_W$. Merely adding individual haplotype scores does not
jointly explain shared observations.

= Observation opportunity and genomic ambiguity

== Exposure

Long contexts have fewer opportunities to fit in reads and survive MEM
extraction. Let $e_f$ be the expected indexed occurrence count per genomic copy
of $f$, per unit haploid sequencing depth. For ideal error-free reads of fixed
length $R$, with complete retention of each feature occurrence,

$ e_f approx max(R - s_f + 1, 0) / R. $

This is an illustrative sampling limit, not the exposure of the current mapper.
Actual exposure depends on the read-length distribution, error process, syncmer
selection, gMEM fragmentation, and overlap convention. Estimate it using the
same mapping pipeline on controlled simulations. Maximality and flanking
context may require a candidate-specific exposure $e_(h f)$; use the simplified
shared $e_f$ only when that approximation is supported. Features with zero or
unreliably estimated exposure cannot support exposure-normalized scoring.

== Whole-pangenome counts

$C(f)$ counts the entire sample, not just reads originating in $W$. Off-locus
copies must not be silently interpreted as local dosage. Either extend to
locus-specific contexts, estimate the off-locus contribution, or use a shared
count factor in joint inference. In a conservative first implementation,
restrict local scoring to features whose panel occurrences are contained in
the locus, while retaining repeats *within* the locus as informative copy number.
Panel specificity is not a guarantee against unrepresented off-locus copies.

Let $beta_f$ represent expected off-locus/error background. Estimate or constrain
it independently of each candidate score; allowing arbitrary per-feature
background to explain every mismatch would destroy identifiability.

== Shared repeat factors across chunks

For chunk $j$, let $n_(j h)$ be the selected chromosomal multiplicity of candidate
$h$, with $sum_h n_(j h) = 2$ in a diploid. Let $a_(j h f)$ count feature
occurrences assigned to that candidate's core under an explicit ownership rule;
overlapping retrieval context must not duplicate a genomic occurrence. A joint
extension is

$ mu_f = d sum_(j=1)^J sum_(h in cal(H)_j)
  e_(j h f) n_(j h) a_(j h f) + beta_f. $

This expression applies only to occurrences whose existence and modeled
exposure are determined by the owning chunk state. Unique ownership alone does
not make boundary-spanning context local. A feature crossing a donor-switch
boundary needs occurrence/exposure terms conditioned on the joint phased states
on both sides (or all states spanned), or an explicit exclusion/approximation.

For example, diploid mosaics `AB | ab` and `Ab | aB` have identical per-chunk
multiplicities for left candidates `A/a` and right candidates `B/b`, but feature
`AB` occurs once in the first and zero times in the second. Fixed candidate-local
terms cannot distinguish them. Include this two-chunk counterexample in
mathematical validation of boundary-feature scoring.

Here exposure may depend on the candidate and must be calibrated. One global
$C(f)$ supplies one shared evidence factor, rather than an independent emission
for every chunk containing $f$. This does not make overlapping features
statistically independent or identify copy placement when contexts cannot
distinguish it. The expression is a proposed model, not an implemented or
calibrated caller.

= Local genotype scoring

== Contextual COSIGT baseline

For features with valid shared exposure, let $y_f = x_f / e_f$ and choose
nonnegative fixed feature weights $omega_f$. A weighted cosine score is

$ S_("cos")(G) =
  frac(sum_(f in cal(F)_W) omega_f y_f g_(G f),
    sqrt(sum_(f in cal(F)_W) omega_f y_f^2)
    sqrt(sum_(f in cal(F)_W) omega_f g_(G f)^2)). $

This baseline assumes negligible or separately handled background. It replaces
node-only features with contextual subwalk features, retaining COSIGT's
haplotype-vector summation. It recovers a node-feature formulation when only
one-node features are used and evidence units agree; it is not automatically
identical to COSIGT's aligned-base coverage pipeline.

Return no call if the weighted sample norm vanishes. Exclude each individual
genotype with zero weighted norm as unscorable before ranking, and return no
call if none remain. Retain tied or indistinguishable genotype classes rather
than assigning certainty through arbitrary ordering. Cosine is scale-invariant: proportional candidate vectors
cannot be distinguished, and absolute CN cannot be recovered from direction
alone.

== Count-based extension

With depth $d$ per haploid copy, define

$ mu_f(G,d) = d e_f g_(G f) + beta_f. $

With candidate-specific exposure, replace $e_f g_(G f)$ by
$e_(a f) a_(a f) + e_(b f) a_(b f)$. A weighted Poisson composite score is

$ S_W(G,d) = sum_(f in cal(F)_W) omega_f
  [x_f log mu_f(G,d) - mu_f(G,d) - log(x_f !)] + log P(G). $

Omit features with $omega_f = 0$ before evaluating terms or checking rates;
masked features contribute neither penalties nor zero-rate rejection. On the
active feature set, use $0 log 0 = 0$ and reject a zero predicted rate with a
positive count, or use a justified positive error floor. $P(G)$ is an explicit genotype prior; a
uniform prior gives an evidence-only comparison. The call is

$ hat(G)_W = argmax_G max_(d >= 0) S_W(G,d). $

Estimate depth globally or constrain it locally for copy-number inference.
Freely profiling depth cannot distinguish proportional genotype expectations.
A negative-binomial model can accommodate extra count dispersion.

Overlapping subwalks, nested contexts, and multiple MEMs from one read create
correlated counts. Consequently this is initially a *composite score*, not a
calibrated joint likelihood or genotype posterior. Select/thin features and
weight redundant context scales; validate uncertainty empirically. Report
rankings, score margins, feature support, and no-call conditions. Do not convert
a cosine residual or an uncalibrated composite score into claimed error odds.

= Windowed recombination inference

Infer locally rather than choosing two panel donors for an entire chromosome.
For windows $W_1, dots, W_J$, let

$ Z_j = (h_j^(1), h_j^(2)) $

be an ordered local haplotype state. A copying-model extension selects

$ hat(Z)_(1:J) = argmax_(Z_(1:J))
  [sum_(j=1)^J S_j(Z_j)
   + sum_(j=2)^J log P(Z_j | Z_(j-1))]. $

The transition model favors compatible source continuation and permits donor
switches, preferably using genomic distance or a recombination map. If local
scores are composite scores, this remains a penalized mosaic objective; their
scale relative to transition penalties must be calibrated before a probabilistic
interpretation is warranted. The additive local-score form also requires
consistent evidence allocation or local factors. Shared dispersed-repeat
counts can couple distant chunks; ordinary local-state dynamic programming
then requires joint treatment or an explicit approximation, not repeated use
of the same count as independent evidence.

Recombination can occur inside a chosen window. Permit subdivision and retain
multiple local states rather than freezing a single early call. Reset inference
at chromosomes/components and preserve unsupported boundaries and phase blocks.
Phase labels are not parental-origin labels.

The same sample BWT serves every window. Assign evidence factors consistently:
reusing the same global feature count in neighboring windows does not create
independent observations. Long gMEMs can supply cross-boundary context, but
linkage between separate gMEMs or mates is unavailable in this representation.
Sequence output additionally requires strand-aware monotonicity, overlap
resolution, and explicit gap handling; the current impg mosaic prototype does
not yet establish these whole-genome guarantees.

= Algorithm and implementation sequence

```text
BUILD SAMPLE INDEX
  map reads against the fixed syng panel
  emit gMEM strings under a declared counting convention
  aggregate identical strings with multiplicities
  construct a boundary-aware, weighted collection BWT

GENOTYPE LOCUS
  retrieve homologous candidate haplotype walks
  select a shared, bounded contextual feature set
  query observed counts from the sample BWT
  count candidate occurrences and establish exposure/background
  enumerate local genotype mixtures and rank their scores
  emit supported calls, ties, uncertainty diagnostics, or no call

OPTIONAL MOSAIC
  combine retained local states chromosome by chromosome
  account for shared evidence and unsupported boundaries
```

First implement the sample count primitive and matching-coverage queries, then
the contextual cosine baseline. Validate the count model before adding global
inference. No per-locus remapping, graph-to-graph evidence translation, or
durable per-read placement table is needed for this sequence of work.

= Validation and existing evidence

Compare BWT queries with exhaustive substring counting on a small explicit MEM
multiset, including duplicates, reverse complements, palindromes, repeated
positions, spacing differences, separators, and overlapping-MEM conventions.
Verify weighted and uncollapsed indexes agree. Test candidate features with no
observations, no-call behavior, proportional-vector ties, depth changes, and
shared off-locus repeats.

For genotyping, freeze HPRCv2 locus candidates and simulated reads. Compare
node-only and contextual models on ordinary loci, C4/RCCX, and amylase/CNV
examples. Separate hold-0, candidate-only hold-out, and strict panel/index
hold-out; strict hold-out requires rebuilding small training-only indexes, not
merely removing truth from the final candidate list. Evaluate actual genotype
or sequence agreement, CN error, phase/switch error, no-call rate, index size,
construction time, and locus-query throughput. Compare representations and
scorers separately.

A completed impg 0.5.0 dictionary-mapping smoke test on the local HPRCv2 index
retained 1,867 of 2,000 HG002 reads, with a median of three anchors per retained
read. Wall time was 24.27 seconds, mostly loading, and peak RSS was 12.33 GiB.
This checks the existing mapper/dataset connection only: it does not validate a
gMEM-BWT, genotype accuracy, full coordinate queries, or genome-scale throughput.

Relevant existing source entry points are `SyngIndex::gbwt_mems_for_walk` in
`src/syng.rs`, local cosine scoring in `src/commands/genotype.rs`, and experimental
MEM-assisted mosaic scoring in `src/commands/infer.rs`. They are foundations,
not implementations of this sample-index design.

= References and build

- #link("https://doi.org/10.1186/s13059-026-04242-4")[Bolognini et al., COSIGT, Genome Biology (2026).]
- #link("https://davidebolo1993.github.io/cosigtdoc/workflow/workflow.html")[COSIGT workflow and coverage semantics.]
- #link("https://github.com/pangenome/likegt")[LikeGT: Rust coverage-based genotyping toolkit.]
- Repository notes: `docs/syng-gmem-bwt/benchmark.md` and
  `docs/syng-gmem-bwt/evidence-notes.md`. Earlier read-provenance proposals in those
  notes are superseded by this compact sample-index design.

This document has no external Typst package dependencies. From the repository
root, using Typst 0.14.2:

```sh
mkdir -p target/docs
typst compile docs/syng-gmem-bwt/genotyping.typ \
  target/docs/syng-gmem-bwt-genotyping.pdf
```
