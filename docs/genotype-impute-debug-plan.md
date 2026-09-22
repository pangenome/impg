# Genotype/Impute Debug Plan

This note records the current state of `impg genotype`, `impg map`,
`impg project`, and `impg infer` for locus-specific genotyping and syng-backed
imputation-style stitching. It is intentionally concrete: what works now, what
does not, what evidence units are being scored, and which debug artifacts a
human can inspect.

## Current Input Modes

`impg genotype cos` has two implemented scoring backends.

| User input | Directly accepted by `genotype cos` today? | Current path |
| --- | --- | --- |
| Syng index prefix / sidecar | Yes | `-a/--index panel.syng -r PATH:start-end` builds syng syncmer-node candidate vectors and scores them against `--pack` or syng `--proj`. |
| Syng-native render bundle | Yes | `--render-bundle` is routed to the syng backend when the manifest `feature_space` is `syng-syncmer-node`. |
| Prebuilt GFA | Yes | `--graph local.gfa --pack sample.graph.pack.tsv` uses the GFA segment backend. `P` paths and `W` walks become candidate haplotypes. |
| GFA-segment render bundle | Yes | `--render-bundle` is routed to the graph backend when the manifest `feature_space` is `gfa-segment` or `variation-graph-node`. |
| Syng-derived local GFA built on demand | Yes, indirectly | Supplying `-a panel.syng -r PATH:start-end --sequence-files ... --gfa-engine ... --pack graph.pack.tsv` builds a temporary query GFA and then uses the graph backend. |
| Binary/text pack | Yes | `--pack` is the aggregate coverage/support vector. TSV may be plain or compressed through niffler; binary uses the `IMPGPKB1` format. |
| Projection bundle | Yes | `--proj` resolves to its `sample.pack` path. For GFA projection bundles the graph metadata is also checked. |
| GAF | Not as a direct genotype evidence argument | For GFA graphs, first run `impg project --gfa local.gfa --gaf reads.gaf -O sample.gfa.proj` or `-o pack-tsv`; then genotype the resulting pack/projection. For syng, `impg map -o proj` includes `reads.gaf.zst`, but local `genotype cos` only uses the pack component. |
| FASTA / FASTQ reads | Not directly | Use `impg map -a panel.syng -q reads.fa/fq -o pack|pack-tsv|proj` first. |
| BAM / CRAM | No | There is no BAM/CRAM reader in the genotype evidence path. Convert or align elsewhere, then project into a supported pack/projection format. |
| Raw FASTA panel sequences | Not as sample evidence | FASTA/AGC panel sequence files can build syng indexes, query local graphs, and materialize graph output. They are not read-evidence inputs to `genotype cos`. |

The practical boundary is that `genotype cos` scores a typed feature vector. It
does not align reads itself. Read/query sequences become syng feature evidence
through `impg map`; GAF graph walks become GFA-segment evidence through
`impg project`.

## GFA And Syng-Derived Graph Boundaries

The syng backend operates in `syng-syncmer-node` feature space. Candidate
haplotypes are intervals on indexed syng paths found by querying a target
range. The candidate vector is an unoriented count vector over unsigned syng
syncmer node IDs.

The GFA backend operates in `gfa-segment` or `variation-graph-node` feature
space. It accepts:

- `--graph <local.gfa>`;
- graph render bundles with graph-node feature space;
- dynamic local GFAs constructed from a syng query when `--sequence-files` and
  a graph engine are supplied.

GFA `S` lines define features. Segment IDs are selected with
`--graph-feature-id-mode`:

- `auto`: numeric positive `S` names use `segment-name`, otherwise dense import
  order;
- `segment-name`: require every `S` name to be a unique positive `u32`;
- `dense`: feature IDs are `1..N` in `S`-line import order.

GFA `P` paths and `W` walks define candidate haplotypes. `S` lines with real
sequence use sequence length. `S` lines with `*` must carry `LN:i:<length>`.
Target intervals are path-coordinate intervals over those segment lengths.

Dynamic syng-derived local GFA genotyping is therefore graph genotyping, not
syng syncmer genotyping. The graph pack must be in the graph's feature space.
A syng pack cannot be silently reused as a GFA-segment pack because graph
genotype checks feature-space metadata and, when present, graph ID,
feature-ID mode, and contribution model.

## Evidence Allocation Semantics

### Syng Pack Evidence

`impg map -o pack` and `impg map -o pack-tsv` currently count distinct syng
nodes per retained read:

```text
for each read/query record:
  collect matched syng syncmer node IDs
  sort and deduplicate node IDs within that read
  require at least --min-anchors distinct IDs
  add +1 to each retained distinct node
```

Consequences:

- Repeated occurrences of the same syng node inside one read do not increase
  that node's pack count.
- The pack count is not weighted by read length, base quality, mapping quality,
  node length, or orientation.
- `--min-anchors` for pack output means distinct node count after per-read
  deduplication.
- Syng GAF output is occurrence-preserving and oriented; a syng projection
  bundle can therefore show repeated GAF steps that collapsed in the pack.

Syng genotype candidate vectors are different from the sample pack vector:
candidate intervals count every path occurrence of a syng node. This mismatch
is now documented and covered by
`test_repeated_node_pack_dedup_and_cnv_counting_counterfactual`.

### GFA Projection Evidence

`impg project --gfa local.gfa --gaf reads.gaf` projects oriented GAF graph
walks onto GFA segments. It counts every contributing segment visit in the
aligned GAF path interval. Repeated visits to the same segment in one read are
counted again, and the projection bundle writes `read-contributions.tsv` so a
human can audit each read step.

For graph scoring, `--graph-contribution-model raw` uses raw visit counts.
`--graph-contribution-model length-normalized` divides sample counts by segment
length and gives each candidate a partial-coverage weight of
covered-bp / segment-length for target intervals. This makes unequal GFA node
lengths explicit instead of hidden in a raw node-count vector.

## Debug Artifacts Available Today

The following intermediate products can be emitted now:

- `impg map -o pack-tsv`: sorted text support vector over syng node IDs.
- `impg map -o pack`: compact binary support vector with retained-record and
  raw-anchor metadata.
- `impg map -o gaf`: syng read syncmer walks, including orientation and query
  syncmer positions.
- `impg map -o proj`: syng projection bundle containing `manifest.json`,
  `sample.pack`, and `reads.gaf.zst`.
- `impg project -o pack-tsv`: typed GFA graph pack TSV with graph metadata.
- `impg project -o proj`: GFA projection bundle containing `manifest.json`,
  `sample.pack.tsv`, a GAF copy, and `read-contributions.tsv`.
- `impg genotype cos --emit-report <path>` or `--debug-report <path>`:
  sectioned TSV report for one genotype call.
- `impg infer --emit-mosaic`, `--emit-fasta`, `--emit-gfa`: stitched local
  call products when beam stitching is enabled or requested.

The genotype debug report is the main human-readable locus artifact. For syng
it records input metadata, pack counting semantics, selected locus sample
counts, candidate node counts, overlap summaries, top genotype scores, and
per-node dot contributions. For GFA it additionally records graph source,
graph ID, feature-ID mode, contribution model, segment lengths, sample raw
counts, sample weights, candidate segment counts/weights, and result feature
decomposition.

Remaining debug gaps:

- Syng pack output still has no per-read contribution table analogous to GFA
  `read-contributions.tsv`.
- Candidate filter rejects from `--min-anchors`, `--min-span-fraction`, and
  `--candidate-top-k` are not emitted as a separate rejected-candidates table.
- Infer read-walk evidence does not yet dump GBWT MEM intervals and filtered
  read/candidate hits as a compact debug table.
- BAM/CRAM evidence projection is absent.

## GBWT Imputation Status

There are three distinct concepts that should not be conflated:

1. `genotype cos`: local genotype scoring over coverage vectors. This is a
   COSIGT/LikeGT-style cosine scorer. It is not a recombinant GBWT imputer.
2. syng GBWT lookup: syng indexes are GBWT-backed and are used to find
   homologous path intervals, walk path ranges, and compute MEMs over signed
   syncmer walks.
3. `infer --stitch beam`: a current imputation-like phase/mosaic layer over
   local genotype calls. It can use syng GAF read walks and syng GBWT MEMs to
   reward read-supported transitions, then emit mosaic TSV/FASTA/GFA products.

What is implemented today is local coverage-vector genotyping plus a beam
stitcher over retained local states. The stitcher is a useful early
haplotype-copying/imputation surface, but it is not yet a formal GBWT panel
prior or full recombinant path decoder. The intended direction is to preserve
the current pack/cos evidence model while adding explicit panel transition
priors, richer read likelihoods, and debuggable MEM/read-link tables.

## Focused Tests Added Or Extended

The current focused coverage includes:

- `test_pack_tsv_matches_independent_syncmer_vector_for_fasta_and_fastq`
- `test_genotype_truth_known_homozygote_heterozygote_and_decoy`
- `test_repeated_node_pack_dedup_and_cnv_counting_counterfactual`
- `test_projection_bundle_gaf_read_walk_evidence_rewards_recombinant_stitching`
- `gfa_gaf_projection_writes_expected_counts_metadata_and_repeated_visit_debug`
- `project_cli_bundle_can_feed_graph_genotype_without_pack_feature_space_override`
- `star_segments_use_ln_tags_for_candidate_intervals_and_weights`
- `repeated_gfa_path_steps_are_counted_in_candidate_vectors`
- `genotype_cos_gfa_debug_report_exposes_lengths_repeats_and_scores`

These tests cover FASTA/FASTQ syng mapping, known homozygote/heterozygote
calls, decoys, repeated syng nodes, recombinant read-link stitching, GFA GAF
projection, repeated GFA segment visits, explicit GFA segment lengths, and the
human-readable graph genotype report.

## Validation Commands

Run the focused workstream checks with:

```bash
cargo test -p impg --lib commands::genotype::tests::star_segments_use_ln_tags_for_candidate_intervals_and_weights
cargo test -p impg --lib commands::genotype::tests::repeated_gfa_path_steps_are_counted_in_candidate_vectors
cargo test --test test_genotype_gfa -- genotype_cos_gfa_debug_report_exposes_lengths_repeats_and_scores --nocapture
cargo test --test test_gfa_projection -- --nocapture
cargo test --test test_genotype_validation_suite -- --test-threads=1
```

Broader validation for this branch should include:

```bash
cargo build
cargo test
cargo install --path .
```

If full `cargo test` is too slow in a local iteration, the focused commands
above cover the genotype/impute evidence surfaces touched by this workstream.
