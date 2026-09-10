# Partition discovery, high-copy states, and genome inference

Status: execution-plan clarification. The new sample gMEM-BWT scorer and reliable
whole-genome mosaic inference are not implemented. Partitioning is the immediate
repair task, not evidence that genotyping has already succeeded.

## Three different boundaries

1. **Computational partitions** organize a cover of physical source-path
   intervals. They are not disjoint sets of graph nodes.
2. **Inference chunks** contain validated homologous, spanning candidate walks,
   with explicit orientation, coordinates, repeated occurrences and links.
3. **Recombination boundaries** are inferred transitions. They can occur inside
   either of the above and require refinement rather than forced edge switches.

A seed such as `S288C#0#chrI:100000-110000` discovers intervals at different
coordinates and lengths on other paths. Independent numerical slicing of every
assembly is not homology discovery. A seed-only partition can account for source
bases without supplying any additional homolog or callable genotype evidence.

Assemblies supply the possible paths through repeats. We do not reconstruct
arbitrary paths from disconnected node coverage. Candidate validation must
separate complete spanning walks, partial paths, unknown sequence, and missing
candidates. BED source coverage alone establishes none of those properties.

## Latest partition evidence

The prerequisite fixes are separate draft PRs:

- [#241](https://github.com/pangenome/impg/pull/241): exact collinear seed-block
  normalization and legitimate seed retention, fixing fragmentation and a
  no-progress loop without relaxing scaffold thresholds.
- [#242](https://github.com/pangenome/impg/pull/242), implementation `c707e2e`:
  full-original-window ambiguity exclusion, explicit metadata-v2 empty paths,
  and safe rejection of unloadable seed-free native graphs. Independent review,
  413 library tests, 39 integration tests and all three Rust CI platforms passed.

The original index contains ambiguity-derived fake exact seeds and is retained
only for diagnosis. Fresh construction is mandatory; legacy load/resave or
position-sidecar repair does not clean it. The replacement is at
`~/yeast/syng-k63-s8-seed7-acgt-only/`, including its matching frozen binary,
source patch, input/binary/index hashes, and build manifest. It preserves all
9,901 source IDs, names and lengths (3,336,976,856 bp). Older binaries reject v2.

The hard queries `S288C#0#chrXII:450000-460000` and
`S288C#0#chrXII:460000-470000` now take 0.123s and 0.015s for raw lookup plus
chaining, excluding loading. Fully contained stored windows overlapping
ambiguity fall to zero. Neighboring N-free repeat controls keep their stored
step counts, unique-node counts and occurrence totals. This fixes fabricated
seeds, not the biological difficulty of real repeats.

The subsequent unchanged-parameter full-panel run **timed out at 1,200s**:

| Measurement | Result |
|---|---:|
| Written partitions | 3,876 |
| Source intervals | 328,536 |
| Represented source paths | 9,339 |
| Emitted source bases | 2,959,027,068 (88.6739%) |
| Uncovered source bases | 377,949,788 |

All emitted intervals are in bounds; summed interval length equals coordinate
union length. The wrapper's exit0 means it successfully recorded timeout/QC,
not that partitioning completed. No inference catalog is accepted.

Artifacts: `~/yeast/partition-acgt-only-w10k-d1k-full-1200s/` contains the exact
command, run log, report and partial BEDs. Parameters remain window10kb,
merge1kb, anchors5, fraction0.5, haplotype selection, 16 threads, no singleton
rehoming. Starting S288C nuclear sequences seed discovery; they do not scope it.

## Next repair: discovery scheduling

### Confirmed source behavior, not yet a complete performance diagnosis

At `c707e2e`, `select_and_window_sequences` in `src/commands/partition.rs` selects
sample/haplotype groups using missing bases but generates windows across their
whole paths. The drain loop calls the backend for each queued window without a
fresh missing-core check. The syng wrapper's BFS/DFS methods in `src/lib.rs`
ignore the supplied masks and dispatch full queries; normal partition masking
happens after querying/chaining. Other backends must be checked separately.

One late 10kb query at `ANM#4#block75_contig3:0-10000` contributed only 211 new
source bases. This is evidence of low incremental coverage, **not proof that
its flanking context or all its discovery work was unnecessary**. We have not
measured what fraction of total runtime is redundant querying versus legitimate
high-copy matching, chaining, masking, or writing.

### Bounded repair sequence

1. Add focused per-query counters/timings: queued and dispatched queries,
   current missing core bp, query/context span, raw hits/anchors, passing chains,
   query/chaining/masking time, newly assigned source bp, and repeat visits.
   Avoid per-base/per-read retained state or huge debug dumps.
2. Reproduce the scheduling issue on small multi-path fixtures: a stale queued
   window covered by an earlier query; a mostly covered selected group; a tiny
   residual needing surrounding context; and repeated copies at distinct source
   coordinates. Also test directional/asymmetric discovery.
3. Separate an **unassigned core** from its **query context**. Generate coverage
   work from remaining source intervals, preserving useful flanks and the
   declared threshold denominator. Do not silently shrink a 10kb query to a
   tiny core and thereby change eligibility.
4. Revalidate work when dequeued. Suppress demonstrably redundant work, but
   **do not assume that a covered source window cannot discover new homologs**.
   A covered B can still reveal C that the earlier A query missed. Preserve an
   explicit discovery obligation or demonstrate an equivalent later query.
   Exact-repeat suppression likewise needs a proof appropriate to backend and
   mask behavior, not a generic visited-node rule.
5. Keep source ownership occurrence-aware and monotone. Distinct repeat copies,
   orientations and coordinates survive. Successful seed processing consumes
   its legitimate remaining source interval; backend errors remain errors.
6. Measure the bounded real-data replay, then rerun the unchanged full panel.
   If genuine high-copy expansion dominates, isolate it separately before
   changing matching internals. Do not hide it with frequency masks, weaker
   thresholds, dropped sequences, or a larger timeout as the first response.

A greedy scheduler change may alter partition grouping. Validate source union,
known homolog/copy recall and context-sensitive fixtures, not only speed or
byte-identical BED output. Scope syng-specific assumptions to syng; do not
silently change PAF/BFS/DFS semantics. Check minimum-size omissions explicitly.

**Acceptance:** fewer demonstrably redundant backend calls; no lost source
occurrences or required homology in adversarial fixtures; errors and seed
progress preserved; bounded yeast completion and resource measurements. Source
coverage, validated homology and callable coverage remain separate gates.

## High-copy regions are normal states

Three cases need different handling:

- Many panel assemblies share a locus: these are alternative candidate
  haplotypes, not hundreds of copies in the sample.
- A spanning haplotype contains a tandem array: its walk carries the internal
  repeat multiplicity. A diploid call still selects two chromosomal haplotypes,
  not an independently chosen state for every repeat copy.
- Similar copies occur at dispersed loci: retain their location/occurrence
  identities. Shared features couple inference across those loci; do not assign
  the same observed count independently to every partition.

Candidates indistinguishable under available features can be represented by
classes, but source members and continuation compatibility must survive. Equal
local feature vectors alone do not justify discarding different boundary links.
If the data cannot resolve copy placement or arrangement, report that ambiguity.
Assembly context gives possible configurations, not unlimited read linkage.

## Count and mosaic contracts

Two indexes, different roles:

- Panel syng GBWT: discover candidate regions and supply assembly walks.
- Sample gMEM-BWT: weighted contextual substring counts across observed MEMs.

For feature f and weighted MEM collection {(M_t, w_t)},

    C(f) = sum_t w_t * occ(f, M_t)

counts MEM-substring occurrences, not distinct reads. Full signed-node order,
spacing, collection boundaries, overlap convention and RC/palindrome rules must
be fixed. No durable read-ID or read-by-haplotype matrix is required.

Let n[j,h] be the selected chromosomal multiplicity of candidate h at chunk j.
For a diploid, sum_h n[j,h] = 2. Let a[j,h,f] count genomic occurrences owned by
that candidate's core, with each occurrence assigned once even when retrieved
through overlapping context. A proposed shared-feature expectation is

    mu[f] = d * sum_j sum_h e[j,h,f] * n[j,h] * a[j,h,f] + beta[f]

This expression is restricted to occurrences whose existence and modeled
exposure are determined by the owning chunk state. Ownership prevents duplicate
counting; it does not make boundary-spanning features depend on only one chunk.
Features spanning donor switches require occurrence/exposure terms conditioned
on the joint phased states they span, or an explicit exclusion/approximation.

The mathematical test plan must include this counterexample: left candidates
A/a and right candidates B/b give the same per-chunk multiplicities in diploid
mosaics `AB | ab` and `Ab | aB`, but feature AB occurs once in the first and zero
times in the second. A candidate-local count formula must not claim to distinguish
them; a joint boundary-feature model must reproduce the exhaustive counts.

Exposure e is calibrated through the actual extraction pipeline and beta is
constrained background. Candidate-dependent exposure is allowed; this is a
model extension, not a calibrated implementation. Overlapping features also
remain statistically dependent. Pure cosine cannot identify absolute dosage
when candidate vectors differ only by scale.

One global C(f) gets one shared evidence factor, not a fresh independent count
for each chunk. Ordinary additive local-score/Viterbi inference is appropriate
only when evidence has been allocated consistently or the relevant factors are
local. Dispersed repeats can create nonlocal coupling requiring joint treatment
or an explicit approximation. A catalog order alone does not remove this issue.

Mosaic transitions use actual source-coordinate continuation, orientation and
chromosome identity, while allowing donor switches. Retain alternatives and
refine chunks for internal switches. Reset at chromosomes, trim overlaps,
handle gaps explicitly, and do not infer sequence absent from the panel. Phase
beyond a single stored MEM is panel-model evidence, not preserved mate linkage.

## Execution order

1. Repair and independently review partition scheduling; complete bounded yeast
   discovery without hiding remaining source sequence.
2. Validate chunk homology, spanning candidates and source-coordinate links;
   reconcile this contract with the formal specification as tests settle it.
3. Implement and exhaustively test the sample count primitive and contextual
   scorer, then chromosome-safe mosaics and the frozen truth regimes.

No new inference machinery or mandatory graph projection is a prerequisite for
step1. The unresolved stored-path reverse-complement sketch asymmetry remains a
separate candidate-recall blocker, not something a scheduler repair fixes.
