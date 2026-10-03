# syng-gMEM-BWT genotyping

**Proposed algorithm, not an implemented impg command.**

Build a compact sample BWT over syng graph MEM strings and their multiplicities.
Query it with candidate haplotype walks to obtain contextual matching coverage,
then score haplotype mixtures locally. Reuse the same sample index across loci;
no per-locus remapping, per-read alignment database, or bluntified variation
graph is required.

## Documents

- **[Formal specification](genotyping.typ)** — the authoritative proposal:
  weighted substring counts, shared contextual features, exposure and off-locus
  background, contextual COSIGT/count scoring, and windowed mosaic inference.
- [COSIGT/LikeGT audit and HPRCv2 pilot](benchmark.md) — source findings,
  completed dictionary-mapping smoke test, and proposed truth-based comparisons.
- [Existing impg capabilities and earlier alternatives](evidence-notes.md) —
  evidence semantics, graph-translation limitations, and whole-genome caveats.

The supporting notes preserve earlier read-by-candidate/provenance proposals
for context. Those storage proposals are superseded by the compact MEM-BWT
specification; they are not requirements for the algorithm.

## Build the specification

From the repository root, with **Typst 0.14.2**:

```bash
mkdir -p target/docs
typst compile docs/syng-gmem-bwt/genotyping.typ \
  target/docs/syng-gmem-bwt-genotyping.pdf
```

No external Typst packages are needed. Generated files stay under the ignored
`target/` directory. The `Typst specification` GitHub Actions workflow compiles
the document and uploads the PDF as the `syng-gmem-bwt-genotyping` artifact.

## First implementation milestone

Implement the sample count primitive and matching-coverage query, validate them
against exhaustive counting on small MEM multisets, then compare contextual
cosine with node-only coverage on fixed-truth HPRCv2 loci. The count-likelihood
and chromosome-mosaic extensions require additional calibration and validation.
