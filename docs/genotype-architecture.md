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

The current backend is syng:

- graph features are syng syncmer node IDs
- sample evidence comes from `impg map -o pack` or a projection bundle from
  `impg map -o proj`
- candidates are path intervals gathered by querying a reference path range
- candidate vectors are syncmer-node traversal counts through those intervals

This is intentionally the first backend, not the whole model. Other projection
methods should be able to feed the same `pack` abstraction:

- MEM/syncmer matching to graph features for high-quality WGS
- read alignment to local variation graph nodes
- read alignment to extracted haplotypes for short or damaged ancient DNA
- condensed graph or bubble-level features

`pack` is the compact binary coverage-vector format used by the current syng
path. `pack-tsv` is the human-readable export of the same abstraction, and
`packbin` remains a compatibility alias for compact pack output. Projection
bundles add metadata and read-walk evidence around the pack so packs cannot be
accidentally combined with the wrong graph or feature namespace.

See [`pangenome-genotyping-roadmap.md`](pangenome-genotyping-roadmap.md) for
the longer-term plan: backend-neutral pangenome genotyping, local variation
graph backends, alignment evidence, learned local scorers, and phased copying
inference.

For rendered graph products, the genotype layer must preserve source path
identity through Pan-SN-aware translation tables. See
[`render-gbz-translation-design.md`](render-gbz-translation-design.md).
