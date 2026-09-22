# IMPG Genotype Architecture

`impg genotype` is organized around graph-feature evidence rather than a
particular graph representation.

The shared model is:

1. Choose a locus.
2. Extract candidate haplotypes or paths that cover the locus.
3. Represent each candidate as a vector over graph features.
4. Represent the sample as a coverage/support vector over the same feature IDs.
5. Score ploidy-sized candidate combinations.

For `impg genotype cos`, the scoring step is COSIGT/LikeGT-style cosine
similarity over graph-feature coverage. `cosigt` is accepted as an alias for
`cos` to make that lineage explicit:

```text
sample vector S
candidate haplotype vectors H_i
genotype vector G = H_i + H_j + ...
score = cosine(S, G)
```

The command now has two concrete `cos` backends:

1. **Syng syncmer-node backend**:
   - graph features are syng syncmer node IDs;
   - sample evidence comes from `impg map -o pack` or a projection bundle from
     `impg map -o proj`;
   - candidates are path intervals gathered by querying a reference path range;
   - candidate vectors are syncmer-node traversal counts through those
     intervals.

2. **GFA segment backend**:
   - graph features are imported GFA `S` segments;
   - graph sources are `--graph <local.gfa>`, local graph render bundles with
     `feature_space = gfa-segment`, or a dynamic query graph built from
     `-a/--index`, `-r/--target-range`, `--sequence-files`, and
     `--gfa-engine`;
   - GFA `P` paths and `W` walks become genotype haplotype candidates;
   - candidate vectors count segment traversals, with
     `--graph-contribution-model raw` or `length-normalized`;
  - sample evidence comes from a typed graph pack or a projection bundle
    produced by `impg project --gfa graph.gfa --gaf reads.gaf`;
  - graph evidence is guarded by explicit feature-space, graph ID, feature ID
    mode, and contribution-model metadata.

These are intentionally only the first backends, not the whole model. Other
projection methods should be able to feed the same typed coverage abstraction:

- MEM/syncmer matching to graph features for high-quality WGS
- read alignment to local variation graph nodes
- read alignment to extracted haplotypes for short or damaged ancient DNA
- condensed graph or bubble-level features

`pack` is the compact binary coverage-vector format used by the syng path.
`pack-tsv` is the human-readable export of the same abstraction, and `packbin`
remains a compatibility alias for compact pack output. For GFA graph
genotyping, a pack must declare `feature_space = gfa-segment` or
`variation-graph-node` in TSV header metadata, a `.meta.tsv` / `.metadata.tsv`
sidecar, or the CLI override `--pack-feature-space`. Typed GFA packs produced
by `impg project` also declare `graph_id`, `feature_id_mode`, and
`graph_contribution_model`; `genotype cos --graph ... --proj ...` resolves the
bundle pack and checks the same metadata. If a graph ID is declared, it must
match the loaded graph. This prevents a syng-node pack from being silently
treated as graph-node evidence.

GFA feature IDs are selected by `--graph-feature-id-mode`:

- `auto`: use numeric segment names when every `S` name is a unique positive
  `u32`; otherwise use dense import order;
- `segment-name`: require numeric `S` names and use those as pack feature IDs;
- `dense`: assign feature IDs `1..N` in `S`-line order.

The implemented GFA projection slice accepts GAF records whose path field is an
oriented graph walk over GFA segment names, such as `>segA<segB>segC`. It
counts every contributing segment visit in the aligned graph interval; unlike
the syng pack builder, repeated visits to the same segment within one read are
not deduplicated. Projection bundles include `sample.pack.tsv`,
`manifest.json`, a copy of the source GAF, and `read-contributions.tsv` so the
aggregate pack can be audited back to individual read-walk steps.

## Evidence Debug Reports

`impg genotype cos` accepts `--emit-report <path>` (alias:
`--debug-report <path>`) to write an optional sectioned TSV text report without
changing the normal genotype TSV on stdout or `-O`. The report is intended for
human audits and focused regression tests of one locus.

For syng, the report records the resolved syng prefix, pack path, target range,
candidate mode, ploidy, `top_n`, `candidate_top_k`, syng padding/extension,
candidate anchor/span filters, and pack metadata. It also states the current
pack counting semantics explicitly: packs produced by `impg map -o pack` count
each distinct syng node at most once per retained read, so repeated node
occurrences within a single read do not increase that node's sample count.

Report sections expose the selected locus sample vector, candidate vectors,
candidate overlap summaries, top genotype score rows, and feature-level score
decomposition. The `result_scores` section uses the same similarity, QV, dot,
sample norm, genotype norm, haplotypes, regions, anchor counts, and span
fraction formatting as the primary genotype TSV. The `result_features` section
shows `node_id`, `sample_count`, per-haplotype candidate counts,
`genotype_count`, and each node's dot contribution, making repeated-node and
ploidy dosage effects visible.

For GFA graphs, the same report flag records graph source, graph source path,
dynamic graph build command when present, graph ID, feature space,
feature-ID mode, contribution model, node lengths, sample raw counts and
weights, candidate counts, and feature-level dot/norm contributions by segment.
The normal TSV result shape is preserved: rank, method, ploidy, similarity, QV,
dot, sample norm, genotype norm, haplotypes, regions, candidate support counts,
and span fractions.

Common GFA examples:

```bash
impg project \
  --gfa local.gfa \
  --gaf reads.gaf \
  -O sample.gfa.proj

impg genotype cos \
  --graph local.gfa \
  --proj sample.gfa.proj \
  --target-path REF_PATH:0-10000

impg genotype cos \
  --graph local.gfa \
  --pack sample.graph.pack.tsv \
  --pack-feature-space gfa-segment \
  --target-path REF_PATH:0-10000

impg genotype cos \
  --render-bundle locus.impg-gbz \
  --pack sample.graph.pack.tsv \
  --pack-feature-space gfa-segment

impg genotype cos \
  -a panel.syng \
  -r REF:1000-2000 \
  --sequence-files panel.agc \
  --gfa-engine syng:crush \
  --pack sample.graph.pack.tsv \
  --pack-feature-space gfa-segment
```

See [`pangenome-genotyping-roadmap.md`](pangenome-genotyping-roadmap.md) for
the longer-term plan: backend-neutral pangenome genotyping, local variation
graph backends, alignment evidence, learned local scorers, and phased copying
inference.

For the current genotype/impute workstream status, input-mode matrix, evidence
counting semantics, debug artifact inventory, and GBWT-imputation boundary, see
[`genotype-impute-debug-plan.md`](genotype-impute-debug-plan.md).

For rendered graph products, the genotype layer must preserve source path
identity through Pan-SN-aware translation tables. See
[`render-gbz-translation-design.md`](render-gbz-translation-design.md).
