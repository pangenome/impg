# Bounded syng partition discovery repair

## Two distinct failures

Parent-run k63/s8/seed7 yeast probes found that irregular dense exact self
anchors fragmented under the alignment-oriented SweepGA best-buddy heuristic.
For `S288C#0#chrI:100000-110000` and `S288C#0#chrI:110000-120000`,
452/404 exact self anchors became 69/59 chains, with largest spans
1096/1238bp. A .5 per-chain query-span fraction
therefore discarded even self. Increasing scaffold gap did not fix this.
Separately, a deterministic 10kb singleton partition repeated the same missing
interval indefinitely when all query hits were filtered out. The earlier
600-second full-panel timeout alone does not establish that same loop.

## Repair and invariants

- Partition discovery explicitly retains the actual source seed interval before
  normal merge/boundary/mask accounting. Empty successful discovery can emit a
  legitimate seed singleton, not invented homolog support. Already claimed
  bases are still subtracted. Backend errors still propagate.
- Before the existing SweepGA helper, normalize overlapping or touching exact
  collinear anchors into gap-free exact blocks. Bucket by path and strand,
  then signed diagonal (`t-q` forward; `t+q` reverse, for fixed k). Separate
  diagonals/copies remain separate blocks. No unobserved gap is merged here.
- Each PAF-like record uses its interval-union length; it does not count
  overlapping anchors times k as sequence length. Every retained block expands
  to its original anchors, including original node IDs and multiplicity.
- SweepGA still owns chaining, scaffold filtering, and chain IDs. Its pinned
  implementation computes scaffold `total_length = q_max - q_min`, which can
  include gaps. This **scaffold span is not exact covered-base mass**. That
  upstream behavior and all parameter defaults remain unchanged.

Focused tests cover irregular dense forward/reverse blocks, memberships,
multiple paths/strands/copies, touching seeds, one-base gaps, union lengths,
empty successful discovery, masks, and duplicate self accounting. Integration
checks cover filtered self/forward discovery, forced-empty filtered partition
output (exactly the source singleton), and normal multi-path partition output.

## Parent validation and bounded panel replay

All 409 library tests and 38 integration tests across `test_syng_integration`,
`test_genotype_gfa`, `test_gfa_projection`, and
`test_genotype_validation_suite` passed on the repaired native release build.
Independent source review accepted the bounded repair with the limitations
below; that is not an inference-catalog acceptance gate.

The original deterministic `toy#0#chr1` 10kb case now completes in 0.033s with
anchors=5/fraction=.5 and emits exactly the source interval. The original
binary timed out after 3s; debug logging showed 99 selections of the identical
missing interval in 0.2s. Raw and fraction-zero controls also still complete.

With the same 1kb scaffold gap and .5 fraction, before partition masking:

| Query | Passing paths before | Passing paths after |
|---|---:|---:|
| `S288C#0#chrI:100000-110000` | 28 | 256 |
| `S288C#0#chrI:110000-120000` | 0 | 235 |

The exact self anchors now form one chain spanning 10018bp and 10046bp,
respectively, without thinning away anchors or lowering thresholds. Path counts
include source occurrences and are not numbers of validated spanning haplotypes.
The four-window refined CLI pilot also completes: 961 preliminary intervals,
5.48s query time in this run, versus 960 intervals before. These timings are
single-run observations, not a general throughput benchmark.

A fresh 60s full-panel discovery diagnostic used the earlier 10kb windows,
merge distance 1kb, anchors 5/fraction .5, and S288C nuclear starting sequences.
It emitted 186 partial partitions containing 44,138 in-bounds source intervals
on 3,679 paths: 397,232,361 source bp (11.90% of the names-table panel total).
Their coordinate-union size equals their summed lengths, so no source bases
are double-assigned in these partial outputs. The run hit its explicit budget
and is **incomplete, not an accepted catalog**. Unlike the earlier 600s pilot,
it used `--separate-files` to preserve inspectable intermediate BEDs. This is
source-coordinate accounting, not a homology or callable-coverage measurement.

Artifacts remain outside Git:

- `target/experiments/syng-partition-diagnosis/after-*.tsv` and
  `termination-toy/after-results.json`: isolated before/after controls.
- `parent-full-library.log` and `parent-integration.log` in that directory:
  test output.
- `~/yeast/baseline-pilot-chrI-4windows-seed-block-repair/`: full pilot argv,
  executable fingerprint, raw BED and preliminary diagnostics.
- `~/yeast/partition-pilot-w10k-d1k-seed-block-repair-60s/`: bounded discovery
  argv, log, partial BEDs, executable SHA256 and coordinate-QC report.

### Duplicate-anchor API caveat

The old helper deduplicated identical `(query_pos, target_pos, node_id)` triples;
normalization retains supplied memberships, including their multiplicity.
Current discovery callers already deduplicate these occurrences, and review
found no concrete production regression. Nevertheless, deliberately duplicated
triples supplied to this public helper are repeated provenance, **not independent
biological support**; downstream anchor-count filters remain count-sensitive.
Different coordinates, paths and orientations must not be conflated with such
duplicate discovery records.

## Unresolved backend limitations (not fixed or validated here)

### Reverse-complement stored-path discovery

The deterministic fixture retained in
`test_syng_filtered_wrapper_retains_dense_self_and_forward_homologs` contains
10kb self, identical forward copy, and exact reverse-complement copy. With
`SyncmerParams { k: 8, w: 55, seed: 7 }`, querying self `[0,10000)` produced
369 anchors each for self/forward but only 15 reverse anchors across 13 isolated
hits, **before normalization**. All reverse anchors had `q+t=9937`; there was
insufficient connected support for the .5 fraction filter. The wrapper test
asserts only self/forward retention; synthetic dense reverse normalization
passes, but reverse end-to-end discovery is **not validated**.

`src/syng.rs::matched_syncmers_best_query_orientation_impl` explicitly notes
that closed-syncmer extraction is not symmetric under reverse complementation.
The mapper tries both read orientations; stored-path queries use stored
sketches. The observed scarcity is consistent with this pre-existing issue.
No C/backend workaround, threshold reduction, or reverse inference acceptance
is implied by this repair.

Reproduction: the test generator starts a u32 state at 7, repeatedly computes
`state = state * 1103515245 + 12345` with wrapping arithmetic, then maps
`(state >> 16) % 4` to A/C/G/T for 10000 bases. Build three paths named
`self#0#chr1`, `copy#0#chr1`, and `reverse#0#chr1` with the sequences above;
call `query_region_with_anchors_ext("self#0#chr1", 0, 10000, 0, 0)` and inspect
per-hit anchors. Parent scratch evidence is in
`target/experiments/syng-partition-diagnosis/worker-wrapper-raw.log`.

### Entirely seedless short paths

Adding `>short#0#chr1` with `ACGTACGTACGTACGTACGT` to a FASTA containing the
10kb generated path, building with `impg syng -f fixture.fa -o idx
--syncmer-length 63 --smer-length 8 --syncmer-seed 7
--position-sample-rate 1 -t 1`, and loading/querying the short path `[0,20)`
returns `InvalidData: No GBWT path start info for path 1 — index may need
rebuilding`. Its name exists, but no GBWT start exists. This is an error, not
an empty successful hit list; it remains an error before seed retention.
Whole-panel completion including these paths is not established.

A random10kb + N10kb control was unsuitable as a seed-free fixture: the C
indexing path produced repeated seeds across the N run. Querying `[10100,20000)`
returned 19775 anchors in that attempt. This observation is not a validation of
ambiguous-base handling. Scratch evidence: `worker-anchorless-control.log` in
the same diagnosis directory.

## Parameter and validation boundaries

The requested minimum anchor count remains an adaptive cap (5 becomes 2 at
10kb/k63/s8). CLI wrapper filtering is still enabled only when the requested
anchor count is positive; a positive fraction alone does not enable it.
Those control semantics need separate review, not silent reinterpretation.
No completed whole-panel run, inference catalog acceptance, or global
sensitivity/specificity claim is part of this bounded repair. Reverse sketch
orientation, entirely seedless paths, ambiguous-base behavior and resource
adaptation still need follow-up before whole-pangenome inference is accepted.
