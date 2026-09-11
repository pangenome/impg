# Source-supported sequence blocks and independent alignment accounting

These are new **experimental executable stages**, downstream of frozen quantitative
calls and threads. They do not implement substring cover, rebuild the catalog,
change genotype/DP parameters, infer recombinant splice positions, or certify a
complete phased chromosome assembly. Catalog acceptance remains false.

## Commands

After the existing `genome-infer genotype --axis ...` command:

```sh
IMPG=/home/erikg/impg/target/release/impg
# Every output directory must be NEW; its parent must exist.
"$IMPG" genome-infer reconstruct \
  --calls "$QUANTITATIVE/calls.json" --threads "$QUANTITATIVE/threads.json" \
  --panel-names /home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng.names \
  --sources /home/erikg/yeast/yeast235.agc \
  --gap-policy split --out-dir "$RECONSTRUCTION"

# Truth is an input ONLY to the independent evaluation stage.
"$IMPG" genome-infer align-sequence \
  --query "$RECONSTRUCTION/reconstruction.fa" --truth "$WHOLE_TRUTH_FASTA" \
  --wfmash /home/erikg/.cargo/bin/wfmash --threads 4 --out-dir "$ALIGNMENT"

# Alternatively, replay an independently supplied base-level PAF, without an aligner:
"$IMPG" genome-infer evaluate-sequence \
  --query "$RECONSTRUCTION/reconstruction.fa" --truth "$WHOLE_TRUTH_FASTA" \
  --paf "$BASE_LEVEL_PAF" --out-dir "$EVALUATION"
# Optional: --aligner-metadata metadata.json (recorded, not treated as authenticated).
```

Use the exact `.names` sidecar path present on disk (a prefix ending in `.syng`
uses `PREFIX.names`, not `PREFIX.syng.names`). No graph or feature catalog is
loaded. FASTA source files can be supplied instead of AGC; each source path must
be available in exactly one declared file and have the panel's exact length.
FASTA inputs may be compressed and must have unique, nonempty sequence names and
nonempty IUPAC DNA records. **Concatenated XZ streams are not supported** by the
current decoder: they fail with `corrupt xz stream`, rather than silently counting
only the first stream. Decompress such files completely to plain FASTA first.
Single-stream XZ is supported; a regression checks each stream separately and the
combined file's rejection. Full PanSN contig names are retained, not suffix-matched.
Reconstruction does not create FASTA indexes beside input files.

Truth headers such as `S288C#0#chrI:0-230218` can be consumed **literally**. No
coordinate-suffix normalization is necessary for this evaluator. Names in the PAF
must exactly match the first whitespace-delimited FASTA header tokens.

## Reconstruction validation and trust boundary

`src/genome_inference/reconstruction.rs` validates:

- Format/model, haploid scope, count policy, checksum algorithm, experimental
  flags, matching panel identity and nonempty sample/catalog payload digests.
- The supplied names file's byte length/content checksum against the calls' panel
  identity; every source occurrence ID, group, bundle identity, membership,
  coordinate, strand, name and length. Bundle membership must be complete.
- Local status/best-bundle consistency with the frozen scores and fit cutoff.
- The **entire** candidate state set, physical occurrence membership and emissions,
  axis validity, all optimal marginals, selected statuses, DP segments, blocks,
  gaps, switches, coverage and genotype-only/repeated-axis group lists. It creates
  only small source/group metadata and reuses the existing threading function to
  replay the original DP, not a second model or catalog rebuild.
- Reference, explicit-axis and explicit-catalog orientation consistency. Signed
  context evidence is consumed as frozen threading evidence; revalidating the
  underlying anchor contexts would require the feature catalog. Segment score
  sums alone allow `1e-8 + abs(score)*1e-12` round-trip tolerance; selected states,
  block structure and other fields must agree exactly.

All input files and emitted sequences receive FNV-1a64 content fingerprints,
including the whole declared AGC/FASTA files (streamed hashing), calls, threads and
names. Source-file bindings record the exact source path, full length and file
index. These are corruption/provenance fingerprints, **not authentication**. Names
and lengths cannot prove that an arbitrary same-length FASTA is the original
panel sequence; declared source files remain a trust input. The original
sample/catalog digest equality links frozen artifacts but is not a recomputation
of their counts. No truth is accepted or accessed by reconstruction.

## Spelling policy and version-1 outputs

- `reconstruction.fa`: deterministic `block000001` identifiers; uppercase oriented
  full donor spans, reverse-complemented with IUPAC complements when needed.
- `provenance.json`: model/scope, input fingerprints and source bindings,
  `blocks`, `group_ledger`, `gap_ledger`, cross-block reuse diagnostics and totals.
  Each block has exact source extraction coordinates/orientation, original source
  states (bundle, identity, occurrence IDs, scores/evidence), axis memberships and
  span, sequence hash and emitted/scored/bridged/suppressed bp counts.
- `unresolved.json`: unresolved axis rows and reasons, unscored axis gaps,
  unscaffolded groups and repeated-axis groups. These categories do not disappear
  because an independent local genotype exists.
- `unresolved.bed`: zero-based half-open **axis** coordinates with group/reason;
  it is not a claim of precisely mapped missing donor or truth intervals.
- `manifest.json`: authoritative running/succeeded/failed stage status. Success
  means the computation completed, **not** that reconstruction coverage passed.

The default `split` policy breaks at either a positive donor gap or a positive
reference-axis gap. `--gap-policy copy-source` can copy intervening sequence only
inside one existing same-path, known-orientation DP block. Those bp are explicitly
`source_bridged_imputed_bp`; continuity is not genotype evidence for the gap.
Unscored axis gaps remain labelled even when source-bridged. Source and axis spans
are separate quantities; they need not have equal lengths.

Within one DP block, same-source coordinate overlap is spelled **once**, including
across reference-gap splits. Such fragments preserve their break, original source
span and axis memberships, but trim the already emitted physical source overlap
and report `suppressed_within_block_overlap_bp`. This trim does **not** invent an
exact source-to-axis nucleotide correspondence. For a reverse traversal the high
end, not the low end, is trimmed. Duplicate physical BED rows keep all occurrence
IDs without emitting duplicate sequence.

Different source paths, chromosomes, unresolved rows, unknown orientations and
nonmonotonic blocks are never joined with Ns or invented splices. Across
**disconnected** blocks, the same donor template may be reused. Such overlaps are
listed in `cross_block_reused_source_overlaps` as copy/placement ambiguity, **not**
certified independent physical copies. Both literal selected block hypotheses
remain in the whole-output evaluation denominator. Donor reuse alone cannot
establish sample copy number; globally suppressing it would hide excess sequence.
Distinct physical copies on different source coordinates/paths are never removed
merely because their strings are identical.

For an optimal tied row, a separate `sequence-equivalent-unphased-fragment` is
emitted only if **every** optimal state is oriented, passes the fit cutoff, and
spells exactly the same complete oriented nucleotide string. All alternatives
remain provenance. Such fragments are never concatenated into an arbitrary path
of independently selected optimal marginals. Nonidentical/unknown/poor-fit ties
stay unresolved. Ties involving a simultaneous multicopy bundle stay explicitly
`multicopy-placement-tie-not-single-allele`: a general multicopy allele/placement
model is deferred rather than collapsing equal copies or emitting alternatives
as simultaneously present. Even for uniquely placed multicopy bundles, the
ledger says which source placement was spelled and that remaining bundle copies
are genotype-only, not reconstructed.

Empty resolved FASTA is a valid no-call. Unscaffolded, repeated-axis and unresolved
sequence is never counted as emitted; the **whole truth**, not called intervals,
remains the downstream denominator.

## Native alignment settings and raw evidence

`align-sequence` requires the inspected **wfmash b55cf75** option contract. It
records `--version`, the executable's byte length/content fingerprint, arguments,
input/working-copy fingerprints and working directory in `aligner.json`.
Equivalent direct invocation on its reserved working copies is:

```sh
cd "$ALIGNMENT"
/home/erikg/.cargo/bin/wfmash alignment-truth.fa alignment-query.fa \
  -t 4 -H 0 -f -n 10 -S 10 -s 1000 -l 1000 -p 90 -k 19 \
  > alignments.paf 2> aligner.stderr.log
```

There is **no `-m` approximate mode**, no `-4` internal one-to-one filter, no `-X`
self-name skipping and no silently inherited frequency mask. `-H 0` disables
most-frequent-kmer exclusion; `-f` disables mapping filtering. `-n 10` and `-S 10`
explicitly request mappings per query/reference pair for ordinary/short queries.
The inspected help says `-S` applies to queries shorter than segment length.
`k19`, identity threshold90%, segment1000 and minimum homology block1000 remain
real ascertainment limits. Very short, repetitive, ambiguous or divergent blocks
may be missed. These are evaluation settings, not changes to syng sketching,
partitioning or inference. Other defaults are pinned by the tool fingerprint;
other aligners/versions can supply external PAFs to `evaluate-sequence`.

Raw PAF and stderr are retained. Native temporary/index files are isolated in the
new output directory. Working FASTAs preserve every base/name/length, converting
only to uppercase uncompressed IUPAC records; no interval trimming or truth-based
parameter adjustment occurs. Empty query/truth skips alignment explicitly and is
still evaluated. External PAF evaluation retains a copy as `alignments.paf`.

## Independent CIGAR accounting

`src/genome_inference/sequence_evaluation.rs` loads **entire** FASTAs and replays
every PAF alignment against their actual bases, including reverse query traversal.
It checks names, reported lengths, bounds, strands, MAPQ range, one `cg:Z:` tag,
positive operation lengths and exact query/target/column consumption. `M` is
resolved by comparing bases; `=`/`X` must agree with assessed ACGT bases. `I` and
`D` count bases, not events. Optional terminal `S`/`H` must agree with the full
oriented query's clipped coordinates; internal clipping and unsupported operations
fail. Without explicit clipping, flanks remain visible through PAF coordinates
and whole-query unaligned bp. PAF's reported matches are retained, never trusted
instead of replay. Unknown/IUPAC columns are unassessed, including gaps containing
unknown bases, not credited as exact ACGT matches.

The selected one-to-one accounting is deliberately conservative and deterministic:

1. Sort by decreasing `min(query span, truth span)`, decreasing MAPQ, then PAF line.
   **Do not sort by local identity.**
2. Reject an entire alignment if **any** consumed query or truth coordinate
   overlaps a previously selected alignment. Do not truncate/reinvent a CIGAR.
3. Credit each selected query/truth base at most once. Retain rejected alignments,
   overlap reasons, raw span sums, raw independent coordinate unions and selected
   independent coordinate unions. This is not a globally optimal assignment and
   can under-credit partially overlapping or ambiguously mapped sequence.

Coverage is selected alignment-span union (including insertion/deletion spans),
not exact-match coverage. Both complete FASTA lengths, all unaligned bp, total and
selected unknown bp, per-sequence accounting and alternative counts are reported.
Split alignments retain query order and diagnose collinear splits, noncollinear
order, orientation and target-path changes. High local identity never hides these
structural diagnostics. Multiple reconstructed copies cannot both receive credit
for a single truth copy. A greedy decision can miss a recoverable alternative;
this limitation is explicit rather than called globally optimal copy recovery.

Let `M`, `X`, `I`, `D` be selected assessed ACGT matches, substitutions, inserted bp
and deleted bp. Ambiguous columns are excluded from these counts, but not from
whole-genome length/coverage denominators:

```
assessed_columns = M + X + I + D
error_columns = X + I + D
identity = M / assessed_columns
alignment_qv = -10 * log10(error_columns / assessed_columns)
```

No assessed columns gives unavailable identity/QV. Zero observed errors gives
`alignment_qv:null`, `alignment_qv_status:"zero-observed-errors"`, with explicit
assessed bp and coverage, **not infinite/perfect-genome QV**. Empty query has
undefined query coverage and zero coverage of nonempty truth. This QV is
alignment-derived, **not calibrated Merqury or assembly-wide QV**. Missing
alignment coverage is visible, not automatically labelled inference error.

## Reproducible small executable fixture and checks

```sh
export CFLAGS=-I/home/erikg/.cache/impg/native-build-deps/root/usr/include
export LIBRARY_PATH=/home/erikg/.cache/impg/native-build-deps/root/usr/lib/x86_64-linux-gnu
export LD_LIBRARY_PATH="$LIBRARY_PATH"
export CARGO_TARGET_DIR=/home/erikg/impg/target
export CARGO_BUILD_JOBS=4 RAYON_NUM_THREADS=4
export IMPG_TEST_WFMASH=/home/erikg/.cargo/bin/wfmash
# Optional: existing directory in which to retain the tiny fixture, including raw PAF.
export IMPG_TEST_SEQUENCE_ARTIFACT_DIR=/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/sequence-evaluation
cargo test --offline --release -j4 --lib genome_inference -- --test-threads=1
cargo test --offline --release -j4 --test test_genome_inference_quantitative -- --test-threads=1 --nocapture
cargo test --offline --release -j4 --test test_genome_inference -- --test-threads=1
```

The quantitative integration test runs real sample/catalog/genotype executables,
then reconstruction, independent exact-PAF evaluation and (when the explicit
`IMPG_TEST_WFMASH` environment variable is set) the real pinned native aligner.
It asserts exact four-block strings, reverse overlap deduplication, a
sequence-equivalent donor tie, unscaffolded copies and chromosome/unresolved
breaks. It emits 18,000 bp and assesses 18,000 exact matches against **30,000 whole
truth bp**: query coverage100%, truth coverage60%, missing truth12,000 bp. That is
an intentional missing-coverage result, not a perfect genome. Changing truth
leaves reconstruction FASTA and provenance byte-identical. Native smoke is not
silently claimed when the variable is absent.

Unit fixtures additionally give exact mismatch/indel counts (`M4,X1,I1,D1`, seven
assessed columns), reverse clipping and IUPAC accounting, positive donor gap
policies, forward/reverse reference-gap overlap trimming, identical versus
different-sequence ties, multicopy ambiguity, duplicate physical rows versus
distinct physical copies, disconnected template reuse with reduced one-to-one
query coverage, split/noncollinear/chromosome mappings, empty output, malformed
CIGAR/artifacts, new-directory reservation and injected partial-output cleanup.
Final focused validation passed 444 library tests, 85 CLI unit tests, the existing
bootstrap integration test and the extended quantitative/native integration test,
with four build jobs and serialized test execution. The earlier empty-FASTA
niffler failure was fixed and rerun; empty and one-base plain FASTAs now bypass
five-byte compression detection. Full-panel yeast execution, installation,
biological acceptance and held-out rebuilds remain parent-owned follow-up gates.
