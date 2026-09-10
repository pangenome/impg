# Yeast whole-pangenome inference plan

## Objective and scope

Infer assembly-supported genome mosaics from short reads using one compact
sample syng-gMEM-BWT. The local inference units are **homologous pangenome
chunks recovered through impg**, not independently sliced coordinates on every
genome and not disjoint sets of graph nodes.

Start with `/home/erikg/yeast/yeast235.agc`. Yeast is the end-to-end development
panel; HPRCv2 remains a later scale test. Work stays in syng. Spanning assembly
paths provide the haplotypes through repeats, including repeated occurrences
and copy number; bluntification, graph projection and de novo repeat resolution
are not prerequisites.

The current mapper, panel GBWT and local cosine scorer are available. The sample
gMEM-BWT and contextual genotype scorer are **not implemented yet**. Each stage
below has an explicit gate; dataset preparation is not genotype validation.

## Execution state

**Current priority:** repair remaining partition scheduling, then validate the
chunk catalog before implementing sample inference. The
[partition/inference contract and repair plan](partition-inference-contract.md)
records the current source audit, high-copy semantics, shared-count math and
acceptance gates. It supersedes the initial pilot status below:

- Partition/scaffold repair #241 and ambiguity repair #242 are reviewed draft
  PRs, unmerged; Rust CI passed. The ambiguity-safe rebuilt index is
  `~/yeast/syng-k63-s8-seed7-acgt-only/` (metadata v2); retain the old index only
  for diagnosis. Fresh rebuilding is required, not legacy sidecar repair.
- With both fixes and unchanged thresholds, the 1,200s full-panel run timed out
  after 3,876 partitions and 2,959,027,068 emitted bp (88.6739%). Bounds and
  non-overlap QC passed; 377,949,788 source bp remain uncovered. Output:
  `~/yeast/partition-acgt-only-w10k-d1k-full-1200s/`. This is not an accepted
  inference catalog. Diagnose redundant scheduling versus genuine repeat cost
  before increasing budgets or changing policy.

### Earlier baseline and data provenance

- Work branch: `work/yeast-gmem-bwt-inference`, initially based on the reviewed
  specification branch at `b2692a1`. Documentation PR #239 remains separate.
- Panel build completed with the installed impg 0.5.0 binary (upstream Rust
  source `7bd50a8`), k=63, s=8, seed=7, position sampling 256, 16 threads and a
  parallel dictionary prepass: 4m50.51s, 9,967,384 KiB peak RSS, approximately
  448 MiB of index sidecars. The archive contains 9,901 paths / 3,336,976,856 bp,
  with 167 sample-name tokens and 235 sample/haplotype identities, not 235
  independently verified strains. S288C, SK1 and Y12 are present by name.
- Data/build outputs: `/home/erikg/yeast/syng-k63-s8-seed7/`. The input archive is
  not modified. The build manifest records its SHA-256 and command parameters.
- The [baseline harness](yeast-baseline.md) passed 29 fixture tests and a real
  four-window query on `S288C#0#chrI:100000-140000`: 960 candidate intervals in
  8.44s query time, including reverse-strand and multiple same-path occurrences.
  Outputs: `~/yeast/baseline-inventory/` and
  `~/yeast/baseline-pilot-chrI-4windows/`, with a successful hardened rerun in
  `~/yeast/baseline-pilot-chrI-4windows-hardened/`. These are bounded candidate
  diagnostics, not validated homology, spanning haplotypes or genotype results.
- The full-panel partition discovery pilot **hit its 600-second timeout** at
  `~/yeast/partition-pilot-w10k-d1k/` (exit 124). It started from S288C nuclear
  chromosomes. Its last progress record reported partition202 and 16,190,787
  total partitioned bp (0.4852%); no final result files were emitted. This run
  is incomplete and not accepted as a catalog. Diagnose partition discovery
  and source accounting against explicit queries before increasing its budget.
- Public SK1/Y12 reads are downloaded under `~/yeast/reads/`: all four paired
  FASTQ files passed their ENA expected-size and MD5 checks. The frozen file
  list is `ena-manifest.tsv`; verification is recorded in
  `download-report.json`. Assembly-truth correspondence remains unverified.
  No new inference accuracy has been measured.

## 1. Freeze and inventory the panel

**Deliverables:** reproducible build command, input checksum, syng sidecars,
path/sample inventory and basic build resource measurements.

1. Build the full panel once at k=63/s=8/seed=7.
2. Inspect actual sequence names and lengths. Verify archive sample count rather
   than assuming the filename describes complete, haploid or diploid assemblies.
3. Locate S288C, SK1 and Y12 if present. Match source assemblies to public read
   accessions using provenance, not strain-name similarity alone.
4. Classify nuclear, mitochondrial, unplaced and ambiguous source sequences;
   retain all in the inventory. Begin chromosome mosaics on nuclear sequences,
   with other components reported separately rather than silently discarded.
5. Smoke-test coordinate queries and source-sequence retrieval, including a
   reverse-strand homolog where available.

**Gate:** index loads; names/lengths are coherent; query intervals are in bounds;
source spelling can be recovered. Build logs and input identity are retained.

## 2. Establish the pangenome chunk catalog

Three objects must remain distinct:

- **Computational partitions** cover source-path intervals and organize work.
- **Inference chunks** contain homologous, preferably spanning candidate paths.
- **Recombination boundaries** are inferred biological transitions, not assumed
  to coincide with either initial partition or chunk boundaries.

Use `impg partition` with the syng backend as a discovery scaffold. Validate it
before accepting the result as an inference catalog. The current command's
`--starting-sequences-file` controls the initial seeds, not the final scope:
it continues over remaining source intervals. A small chromosome pilot should
use explicit `query` seeds, not pretend that this option limits partitioning.

### Pilot and full catalog

Before expanding the catalog, execute the bounded scheduling repair in the
[partition/inference contract](partition-inference-contract.md#next-repair-discovery-scheduling).
Keep missing core intervals separate from query context; a covered seed can
still reveal a previously undiscovered homolog. A skip optimization must
preserve that discovery obligation, not just pass a coverage-union test.

1. Pilot explicit 10 kb seed intervals across one verified nuclear chromosome,
   with a small selection of ordinary and repeated regions. Query the full panel
   for homologous source intervals; do not cut every candidate at the same
   numeric coordinates.
2. Audit partition discovery separately, initially with bounded time/resources.
   Compare candidate groups and source interval coverage with the query pilot.
   Vary chain-support/extent filters if raw union-find closure merges unrelated
   regions; report effects on candidate recall. These are homology-discovery
   controls, not blanket repeat masks.
3. For each selected seed/partition, query/refine spanning haplotype candidates
   and retain boundary anchors, source path, start/end, strand and eligibility
   status. A repeated node may appear in several chunks or several occurrences
   of a candidate. Preserve those occurrences.
4. Keep core intervals for scoring and flanking context for matching/boundary
   assessment. Do not count overlapping evidence independently in every core.
5. Build source-coordinate predecessor/successor relations between chunk
   occurrences. Strand, chromosome identity, overlaps and unrepresented gaps
   remain explicit; adjacency is not inferred merely from matching node IDs.
6. Scale the validated policy across the full panel, including sequence absent
   from a chosen reference. Record intervals excluded by minimum-size or
   candidate-support rules and uncovered source bases.

**Catalog artifacts (proposed schemas):**

- `sequences.tsv`: source IDs/names, lengths, parsed sample/haplotype identity.
- `chunks.tsv`: chunk IDs, seed/core/flank intervals and discovery provenance.
- `candidates.tsv`: chunk ID, source interval, orientation, boundary support,
  spanning/partial status, and occurrence identity.
- `links.tsv`: compatible source-coordinate continuation between occurrences.
- `catalog-qc.json`: coverage, gaps, overlaps, candidate counts, partial paths,
  anomalously large groups and failures. A list of BED intervals alone is not
  proof of validated homology, source spelling or spanning support.

**Gate:** ordinary chunks have consistent homology/spanning candidates; source
interval coverage is accounted for; no false cross-chromosome joins; repeated
copies stay attached to their spanning haplotype context. Diagnose truncation
or absent assemblies as missing candidates, not a need to resolve the graph.

## 3. Prepare three truth regimes

### In-graph

Simulate deterministic reads from known panel paths, then use real reads from
verified matching assemblies. Run haploid strains individually and combine two
read sets at controlled per-haplotype depth for a pseudo-diploid. Truth remains
in the panel. Measure local path/equivalence-class recovery and artificial
switches across chromosomes.

### Recombinant

Create explicit haploid and diploid mosaics of panel source intervals using
validated homologous joins, including switches inside initial chunks. Simulate
reads from the resulting sequences, not by mixing reads independently per
window (which would omit junction-spanning reads). Retain true breakpoints and
source provenance. Do not add whole recombinant chromosomes to the training
panel; test subwindow inference rather than whole-donor recognition.

### Independent assembly / strict hold-out

Use reads with an independent truth assembly absent from the panel. For
controlled tests, remove the target strain's complete assembly/haplotype set
and rebuild a training-only index and catalog. Truth must not seed candidate
selection or feature/index construction. Keep candidate-only exclusion as a
separate, clearly labelled easier test.

**Verified public read candidates** (ENA project PRJNA340312):

| Strain | Run | Layout |
|---|---|---|
| SK1 | SRR4074258 | paired Illumina WGS |
| Y12 | SRR4074358 | paired Illumina WGS |
| S288C | SRR4074255 | paired Illumina WGS |

SK1 and Y12 were used as haploid inputs to an established pseudo-diploid yeast
benchmark. Verify correspondence with our archive before treating them as exact
in-graph truth. SK1 plus Y12 FASTQs total approximately 4.7 GB compressed.
Use ENA file checksums and record assembly versions, stock differences and
sample ploidy. Download only the selected runs. Store data under `~/yeast/`,
not in Git.

Start at 150/151 bp with fixed-seed subsamples at 5x, 15x and 30x total nuclear
coverage, where total sequenced nuclear bases divided by haploid genome size
is the stated depth (an equal diploid mixture gives half that per haplotype).
Separate technical read effects from biological copy-number/aneuploidy tests.

**Gate:** exact source/assembly truth manifest; fixed reads and sampling seeds;
no hidden candidate/index leakage. Public metadata/assembly claims are checked,
not inferred from a filename.

## 4. Implement the sample gMEM-BWT count primitive

1. Reuse `SyngIndex::gbwt_mems_for_walk` to collect signed-node/spacing gMEM
   strings under a declared extraction/overlap convention.
2. Aggregate identical strings with multiplicities. Build a collection BWT that
   respects MEM boundaries and supports multiplicity-weighted substring counts.
   No durable read IDs, placement tables or read-by-haplotype matrix are needed.
3. Define reverse-complement counting and start/end token alignment. Test that
   queries cannot cross MEM boundaries or lose spacing information.
4. Expose observed subwalk counts and candidate matching-coverage profiles.
   Start with a correct small implementation; choose run-length/other compressed
   representations using measured yeast index size and query throughput.
5. Store panel fingerprint, syncmer parameters, extraction/counting policy and
   aggregate sampling metadata with the sample index.

**Gate:** counts exactly match exhaustive counting on toy MEM collections,
including duplicates, palindromes, repeated occurrences, changed distances,
zero matches and overlapping MEM conventions. Weighted and uncollapsed
collections agree. Measure real yeast build size/time/RSS before optimizing.

## 5. Local contextual genotyping over chunks

1. Derive one shared, bounded contextual feature universe from each chunk's
   eligible candidate haplotypes, preserving zero-count observations. Distinguish
   alternative panel homologs from within-haplotype copy number and dispersed
   copies. Retain occurrence identities and boundary compatibility when grouping
   locally equivalent candidates.
2. Query sample MEM-BWT counts and candidate occurrence multiplicities.
3. Implement contextual COSIGT as the initial scorer. Compare node-only and
   contextual features on the same chunk candidates and evidence conventions;
   compare conventional COSIGT/LikeGT separately when graph/evidence units differ.
4. Estimate exposure versus context span through the actual simulation/mapping
   pipeline. Retain copy-number signal within repeated haplotype paths. Extend
   context or model background when the same pattern has off-locus copies.
5. Retain ranked/equivalent candidate combinations, no-call conditions and
   diagnostics; do not treat geometric QV as genotype posterior confidence.
6. Test smaller/subdivided chunks when one combination fails to explain the
   whole interval. Record candidate recall before attributing errors to scoring.

**Gate:** in-graph recovery and contextual-order tests pass; actual genotype
accuracy is measured, not a high-cosine proxy. Hold-out error is decomposed into
panel representation limits, candidate discovery loss and scorer error.

## 6. Whole-genome mosaic inference and output

Use the chunk occurrence links as a chromosome-specific inference scaffold.
Local states are ordered haplotype combinations; transitions favor source
continuation and allow donor switches. Recombination may occur inside initial
chunks, so permit refinement rather than forcing every switch to their edges.

Before emitting genome sequences, fix the existing stitcher's chromosome reset,
strand-aware continuation, overlap trimming and gap handling. Do not bridge
failed chunks silently or concatenate chromosomes. The sample MEM-BWT supports
context within each stored MEM, not mate/read linkage between different MEMs;
phase beyond that comes from the panel model and must retain uncertainty.
Shared global feature counts must not become independent observations merely
because two chunks query them. The shared-feature expectation in the
[contract](partition-inference-contract.md#count-and-mosaic-contracts) counts each
predicted genomic occurrence once despite overlapping query context. Nonlocal
repeat factors require joint treatment or an explicit approximation; ordinary
additive local-score dynamic programming is not automatically sufficient.

**Outputs:** local calls, per-chromosome haplotype mosaics, support/phase blocks,
sequence FASTA where joins are justified, and explicit unknown/gap intervals.
No claim of de novo recovery of sequence absent from the panel.

**Gate:** whole-genome in-graph and recombinant truth recovered to measured
accuracy; chromosome identities preserved; no fabricated sequence continuity.

## Evaluation and stopping rules

Report by whole genome and by chunk class:

- Callable fraction and uncovered truth/source sequence.
- Unordered genotype/equivalence-class accuracy for in-graph truth.
- Sequence edit/error rate and structural/CN agreement against independent truth.
- Switch errors and breakpoint uncertainty, preserving chromosome-wide phase
  consistency rather than rematching phases independently in every window.
- Candidate recall and missing/partial spanning assembly rates.
- Index build time/size/RSS, mapping time and all-chunk query/scoring throughput.

For hold-out, compute an evaluation-only best attainable panel mosaic under the
same chunk and transition constraints. Separate its representation error from
additional inference error. Do not train on this oracle or present best local
matches with unlimited switches as the same constrained oracle.

Stop at the smallest failing layer: source retrieval, chunk homology, count
semantics, candidate recall, local score, then mosaic assembly. Keep a failing
toy/small-chromosome case before attempting whole-panel parameter sweeps.
Do not make broad default changes merely to improve one yeast locus.

## References

- [Formal algorithm](genotyping.typ).
- [COSIGT/LikeGT and HPRCv2 audit](benchmark.md).
- [ENA PRJNA340312](https://www.ebi.ac.uk/ena/browser/view/PRJNA340312).
- [SK1/Y12 pseudo-diploid benchmark](https://pmc.ncbi.nlm.nih.gov/articles/PMC6022571/).
