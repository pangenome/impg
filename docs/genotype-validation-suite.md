# Syng Genotype Evidence Validation Suite

This report describes the focused validation suite added for syng
pack/genotype/infer evidence allocation. The tests are intentionally small,
truth-known simulations inspired by LikeGT/COSIGT-style count-vector checks:
the sample evidence vector is controlled, the candidate haplotypes are known,
and simple cases must not silently produce strange assignments.

## Test Matrix

| Test | Fixture | Evidence path | What is proven |
| --- | --- | --- | --- |
| `test_pack_tsv_matches_independent_syncmer_vector_for_fasta_and_fastq` | Two-haplotype panel with allele-specific reads plus one read shorter than the syng syncmer length | `impg map -o pack-tsv` against FASTQ and FASTA query files | The CLI pack TSV exactly equals an independently accumulated vector from `SyngMatcher::matched_syncmers_in_sequence`; FASTA and FASTQ inputs produce identical vectors; short reads below syncmer length do not contribute phantom evidence. |
| `test_genotype_truth_known_homozygote_heterozygote_and_decoy` | A/B haplotypes with a third unsampled decoy haplotype sharing the same flanks | `impg genotype cos` over independently written pack TSVs | An easy A/A sample ranks A/A first; an A/B sample ranks A/B first; the decoy is available as a candidate but not selected in the top call. |
| `test_repeated_node_pack_dedup_and_cnv_counting_counterfactual` | Single-copy and duplicated-copy CNV haplotypes, with one full duplicated read | `impg map -o pack-tsv` plus direct cosine scoring over a minimal CNV feature model | Current pack semantics are distinct syng nodes per retained read. Repeated syncmer-node occurrences are present in the read, but the emitted pack vector deduplicates them per read. A counterfactual per-occurrence sample vector flips a one-copy/two-copy ranking, making the semantic choice explicit. |
| `test_projection_bundle_gaf_read_walk_evidence_rewards_recombinant_stitching` | Four haplotypes A/B plus recombinant C/D across two blocks | `impg map -o proj`, compressed GAF read walks, `impg infer --stitch beam` | Projection bundles contain oriented syng GAF read walks; read-link rewards are nonzero across the second partition; beam stitching selects the sampled recombinant haplotypes C and D. |
| `star_segments_use_ln_tags_for_candidate_intervals_and_weights` | GFA with `S` records using `*` sequence plus `LN:i` length tags | GFA graph candidate extraction with `length-normalized` contribution model | GFA path coordinates and candidate weights use segment lengths from `LN:i` when sequence is absent. |
| `repeated_gfa_path_steps_are_counted_in_candidate_vectors` | GFA path that traverses one segment twice | GFA graph candidate extraction with raw contribution model | Repeated GFA path steps increase candidate raw counts and scoring weights. |
| `genotype_cos_gfa_debug_report_exposes_lengths_repeats_and_scores` | Three-path GFA with unequal segment lengths and one repeated segment path | `impg genotype cos --graph --pack --debug-report` | The graph debug report exposes segment lengths, raw sample counts, repeated candidate segment counts, dot contributions, and the final assignment. |

The syng tests live in `tests/test_genotype_validation_suite.rs`. They build
small syng indexes in-process and then exercise the real CLI for the user-facing
map, genotype, and infer behavior. The graph tests live in
`src/commands/genotype.rs` unit tests and `tests/test_genotype_gfa.rs`.

## Current Counting Semantics

`impg map -o pack` and `impg map -o pack-tsv` currently count each distinct
syng syncmer node at most once per retained read. Repeated occurrences of the
same node inside one read are not added multiple times to the sample pack.

This is deliberate current behavior, not an accidental side effect:

- `build_syng_map_pack_chunk` collects matched node IDs for one read, sorts
  them, deduplicates them, applies `--min-anchors` to the distinct-node count,
  and then increments each retained node by one.
- Candidate haplotype vectors still count path occurrences. A candidate that
  traverses a syncmer node twice has candidate count `2` for that node.
- The new CNV counterfactual demonstrates that per-read deduplication versus
  per-occurrence sample counting can change genotype order on a small
  one-copy/two-copy feature model. Under current semantics the deduplicated
  sample vector ranks the single-copy model first; under per-occurrence sample
  counting the duplicated model ranks first.

The report emitted by `impg genotype cos --emit-report` already states this as
`sample_pack_counting_semantics = distinct_nodes_per_read`.

## Node Length Assessment

Node length is not an extra weighting term for the current
`syng-syncmer-node` feature space. A syng node is a fixed-length syncmer under
the index metadata, with length `syncmer_k + syncmer_w` in the current internal
parameter names. Because all syng syncmer-node features have the same sequence
length within one index, weighting the cosine vectors by node length would be a
constant factor and would not change rankings.

Node lengths will matter for future non-syng backends:

- `variation-graph-node` features can have different segment lengths, so raw
  node counts and base-pair coverage are not equivalent.
- Local GFA backends may include gap segments, cloned repeat segments, or
  bluntified fragments whose lengths encode important dosage or support.
- A future backend should record whether its sample vector is in node-count,
  occurrence-count, base-pair, or alignment-likelihood units. Mixing those
  units without explicit metadata would make COSIGT-like scores difficult to
  interpret.

For the current syng path, the relevant ambiguity is counting unit
(`distinct_nodes_per_read` versus per-occurrence), not node length.

## Commands and Observed Outputs

The focused suite was run with the same local build environment needed by the
existing project tests:

```bash
source ./env.sh
export C_INCLUDE_PATH=/home/erikg/htslib-local/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/home/erikg/htslib-local/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=/home/erikg/htslib-local/lib:$LIBRARY_PATH
export PKG_CONFIG_PATH=/home/erikg/htslib-local/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/home/erikg/htslib-local/lib:/home/erikg/.cargo/lib:${LD_LIBRARY_PATH:-}
export LDFLAGS="-L/home/erikg/htslib-local/lib -Wl,-rpath,/home/erikg/htslib-local/lib -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu -lcurl -L/home/erikg/.cargo/lib -Wl,-rpath,/home/erikg/.cargo/lib"
export CMAKE_PREFIX_PATH=/home/erikg/htslib-local
export CMAKE_C_COMPILER="$CC"
export CMAKE_CXX_COMPILER="$CXX"
unset CFLAGS
cargo test --test test_genotype_validation_suite -- --test-threads=1
```

Observed output:

```text
running 4 tests
test test_genotype_truth_known_homozygote_heterozygote_and_decoy ... ok
test test_pack_tsv_matches_independent_syncmer_vector_for_fasta_and_fastq ... ok
test test_projection_bundle_gaf_read_walk_evidence_rewards_recombinant_stitching ... ok
test test_repeated_node_pack_dedup_and_cnv_counting_counterfactual ... ok

test result: ok. 4 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```

Broader validation commands for this branch:

```bash
cargo build
cargo test -p impg --lib commands::genotype::tests::star_segments_use_ln_tags_for_candidate_intervals_and_weights
cargo test -p impg --lib commands::genotype::tests::repeated_gfa_path_steps_are_counted_in_candidate_vectors
cargo test --test test_genotype_gfa -- genotype_cos_gfa_debug_report_exposes_lengths_repeats_and_scores --nocapture
cargo test --test test_gfa_projection -- --nocapture
cargo test --test test_genotype_validation_suite -- --test-threads=1
cargo test --test test_syng_integration -- --test-threads=1
cargo test --all -- --test-threads=1
cargo install --path .
```

Observed branch validation:

- `cargo build`: passed.
- `cargo test --test test_genotype_validation_suite -- --test-threads=1`: passed,
  4 tests.
- Existing genotype/infer tests:
  `cargo test --test test_syng_integration -- --test-threads=1` passed, 27
  tests.
- `cargo test --all -- --test-threads=1`: passed. The run included 320 lib
  unit tests, 76 main CLI unit tests, the new 4-test focused suite, 27 syng
  integration tests, 8 syng start-count tests, 10 transitive integrity tests,
  and the existing integration targets. Repository-marked long/red C4 tests
  remained ignored.
- `cargo install --path .`: passed and replaced `/home/erikg/.cargo/bin/impg`
  from this worktree.

## Before/After Behavior

No scoring redesign was made. The new suite did not expose a production code
bug requiring a fix. The only implementation change is additional validation
coverage and documentation.

Before this suite, a regression in syng pack-vector allocation, FASTA/FASTQ
query parsing, decoy handling, short-read filtering, repeated-node CNV
semantics, or projection GAF read-link use could pass without a direct
truth-known assertion. After this suite, those cases have focused deterministic
coverage.

## Remaining Gaps

- The suite uses tiny synthetic panels. It is meant to catch simple impossible
  assignments, not to validate population-scale calibration.
- Pack semantics are documented as `distinct_nodes_per_read`; choosing
  per-occurrence or base-pair-weighted packs for CNV-sensitive genotyping would
  be a behavioral change and should be introduced as an explicit new mode.
- Current syng pack counts ignore base qualities and mapping ambiguity beyond
  exact syncmer membership. Future evidence backends should make weighting and
  uncertainty metadata explicit.
- Non-syng GFA/variation-graph feature vectors still need their own length and
  unit validation once those backends become first-class genotype inputs.
