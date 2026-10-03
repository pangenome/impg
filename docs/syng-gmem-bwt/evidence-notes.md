# Syng mapping: from anchor support to diplotype inference

Status: proposed next work, based on code at `6dd2e73` (upstream main
`7bd50a8` plus the local Guix environment fix). No new mapper or inference
behavior is implemented by this note.

**Design update:** the authoritative proposal is now
[`genotyping.typ`](genotyping.typ): a compact
sample gMEM-BWT queried by local haplotype walks. The read-provenance and
per-read-table proposals below are retained as earlier alternatives, not
implementation requirements. Variation-graph projection is optional, not on
the genotyping critical path. See [README.md](README.md) in this directory
for the Typst build command.

## What exists

- `impg map -o pack`: unsigned syncmer-node support, counted once per distinct
  node per retained read. This is read support, not aligned-base depth.
- `impg map -o gaf`: oriented matched syncmer occurrences in query order,
  including repeats and `qp:B:I` query positions. The lightweight mapper only
  loads the dictionary, not GBWT topology. These are anchor traces, not
  necessarily contiguous graph walks or base-level alignments. MAPQ is zero;
  the reported path length is anchor count times syncmer length.
- `impg map -o proj`: both products. Currently scans reads twice; pack and GAF
  use distinct-node versus raw-occurrence `--min-anchors` thresholds.
- `impg genotype cos`: local syng node-multiplicity scoring, or explicit GFA
  segment/path scoring with evidence already in that graph's namespace.
- `impg project --gfa ... --gaf ...`: destination-GFA GAF to segment evidence.
  It does not translate global syng anchor IDs into local GFA segment IDs.
- `impg infer`: local calls over ranges/partitions, optionally followed by a
  beam over ordered ploidy states with read-walk emissions and linkage rewards.
  It can emit mosaic TSV, FASTA and GFA. Local retained states originate from
  cosine scoring; later walk support cannot rescue already-pruned states.
- Syng raw/blunt and other local GFA builders exist. Render bundles preserve
  source interval/path-step provenance, but do not yet implement the full
  global-syncmer-occurrence to post-bluntification segment-span relation.

Primary code: `src/main.rs` (`write_syng_map_gaf`,
`build_syng_map_gaf_chunk`, `build_syng_map_pack_chunk`), `src/syng.rs`
(`SyngMatcher`), `src/commands/{genotype,infer,render}.rs`,
`src/projection/converter.rs`, and `src/render_bundle.rs`.

## Evidence model: preserve reads, derive summaries

Use the richer projection object as the evidence source, rather than making
coverage the only durable product. Keep three distinct layers:

1. **Node support:** unique-read support, occurrence counts, and aligned bases
   are different statistics and must be named separately.
2. **Oriented edge support:** evidence for traversing actual graph adjacencies.
   A bidirected edge `(u,v)` is equivalent to `(-v,-u)`. Parallel transitions
   with distinct sequence/distance context may need occurrence-level identity.
3. **Read subwalks:** query offsets, oriented occurrences, alternative placements,
   gaps, molecule identity and confidence. These retain phase and repeat-copy
   context that neither node nor edge marginals can recover.

Do not count consecutive matched anchors as an edge unless adjacency is
validated. Missing anchors may represent errors, unindexed sequence, a skipped
subwalk, or a novel junction. Label anchor links separately from verified edges;
retain unresolved alternatives rather than silently choosing a panel path.

First implementation slice: expose read-by-candidate haplotype matching
statistics in homologous windows, using the existing GBWT MEM and candidate
subwalk machinery. Keep node-coverage cosine as the baseline; edge summaries
can follow rather than gate the read-level scorer. Also converge pack and GAF
onto a shared retained-read evidence stream. When adding topology-aware output,
split/extend traces into validated subwalks and derive edges only from those;
retain dictionary-only mode as explicitly sparse anchor evidence. Declare
graph/index identity, feature space and contribution semantics in versioned
metadata. See [benchmark.md](benchmark.md) for the COSIGT/LikeGT audit,
completed HPRCv2 mapping pilot and proposed fixed-truth comparison.

## Translation to blunt/local variation graphs

The required relation is:

```
global syng occurrence (path, step, source interval, orientation)
  -> source sequence interval
  -> final rendered path occurrence
  -> oriented local segment span(s)
```

Build translation against the final graph after splitting, bluntification,
normalization and sorting, not by assuming node IDs survive these operations.
Source spelling must be verified. Overlapping raw GFA paths and unknown segment
sequence need explicit handling or rejection, not concatenated S-line lengths.
Non-syncmer gap sequence has no direct syncmer observation; do not manufacture
support for it by distributing counts over all nearby sequence.

`docs/syng-to-local-graph-translation.md` sketches the missing overlap table.
Its simple additive projection is a prototype scoring rule, not a general
coverage-preserving mapping: one global count copied onto every panel path can
multiply evidence with panel size. Keep alternative occurrence assignments
explicit and normalize across mutually exclusive placement hypotheses. A read
may legitimately support several consecutive features; conservation applies to
alternative placement weight, not the sum of every feature count.

Aggregate node packs cannot recover occurrence assignment, edge order, or phase
once discarded. Use read subwalks for ambiguous projection, with local alignment
as a fallback when base-resolution GAF/CIGAR or novel alleles are required.

## Local and global inference

For local inference, retain cosine as a transparent baseline, then evaluate
read/fragment likelihoods over candidate haplotype pairs. Nodes, edges and
subwalks derived from the same reads are correlated, not independent likelihood
terms. Include depth, error and ambiguous placements; expose uncertainty and
no-call outcomes. Cosine's QV transform is not a calibrated genotype posterior,
and scale invariance prevents absolute copy-number inference from cosine alone.

For chromosome inference, the existing beam is a useful starting point for a
panel-copying model: local haplotype-pair states, evidence emissions, continuity
and distance-aware recombination transitions, plus molecule-spanning phase
constraints. Infer mosaics of panel haplotypes, not one panel donor pair for the
entire genome. Novel sequence/copy-number states need explicit additional
modeling; syng-to-source mappings alone are not sample genotype evidence.

Before whole-genome sequence claims, fix current stitching boundaries:

- Reset inference and output tracks at chromosomes/components and explicit
  unsupported boundaries; do not silently bridge failed calls.
- Enforce strand-aware source monotonicity. Current `phase_transition` treats
  `curr.start <= prev.end` as continuation, even for backward jumps.
- Resolve source overlaps and actual gaps. Current FASTA output concatenates
  full intervals and inserts a fixed ten Ns for uncertain/recombine joins.
- Preserve phase blocks and ambiguity. Short reads do not establish chromosome-
  wide phase without linking evidence or statistical priors; phase labels are
  interchangeable across disconnected blocks/chromosomes and are not parental
  origin labels without additional evidence.

## Validation gates

- Same nodes/different edges; same node and edge marginals/different long walks.
- Reverse complements, missing anchors, nonexistent adjacency, repeated nodes,
  paralogous placements, and reads rejected by the shared retention rule.
- SNP/indel bubbles, syncmers spanning multiple final segments, repeated source
  occurrences, graph split/merge/renumber invariance, and duplicate panel paths
  not inflating evidence.
- Two chromosomes, reverse-strand tracks, overlapping windows, missing calls,
  uncertain joins and exact source-spelling round trips.
- Held-out haplotypes and recombinant truth; vary depth, errors, read length,
  copy number and panel composition. Compare node-only, node+edge and read-walk
  models using genotype error, switch error, CN error, no-call rate and runtime.

Prioritize evidence semantics and read-by-haplotype window scoring; treat
chromosome-safe stitching as a prerequisite to whole-genome evaluation. Existing
synthetic infer tests establish useful behavior but not chromosome-scale
accuracy, calibrated uncertainty, or genome-scale resource use.

Baseline validation at this revision: all 36 tests passed across
`test_syng_integration`, `test_genotype_gfa`, `test_gfa_projection` and
`test_genotype_validation_suite`, using `cargo test --release --locked` with
one test thread and `CARGO_BIN_EXE_impg` pointing to the installed binary.
