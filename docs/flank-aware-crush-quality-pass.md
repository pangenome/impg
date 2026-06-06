# Flank-Aware Crush Quality Pass

Date: 2026-06-06
Task: `quality-pass-flank-aware-crush`

This quality pass tightened the downstream WG tasks for flank-aware crush before
implementation starts. The task graph records were updated for:

- `design-flank-aware-crush`
- `implement-flank-aware-crush`
- `validate-flank-aware-crush-c4`

## Canonical Model

The downstream tasks now specify one shared model:

1. Identify the original bubble, flubble, nested-bubble, or motif-local target
   interval before adding context.
2. Collect common left and right flanks from the same path/graph occurrence as
   each target traversal. Requested flank length and actual flank length are
   recorded separately because path ends, contig ends, Ns, fragmented
   haplotypes, and missing context can shorten the usable flank.
3. Normalize traversal orientation before resolution. Reverse-orientation
   traversals must carry enough metadata to trim and lace the result back with
   the original path orientation and names intact.
4. Run the selected resolver in global/end-to-end mode on the flanked
   sequences. SPOA, abPOA, POASTA, SweepGA, and future generic graph-output
   paths are required to share the same extraction/trimming contract.
5. Trim the resolved replacement exactly at the recorded left/right flank
   boundaries before lacing it back into the graph. Added flank context is not
   part of the replacement interval.
6. Verify exact full-graph path-spelling preservation and path-name stability
   after lacing.

## Required Edge Cases

The implementation and validation tasks now require coverage for:

- repeated contexts where the same flank sequence appears in multiple
  positions;
- short, zero-length, and one-sided flanks at path or contig ends;
- nested bubbles/flubbles and motif-local regions;
- Ns, fragmented haplotypes, missing traversals, and inconsistent path entry or
  exit;
- reverse-orientation traversals and reverse-complement trim coordinates;
- stable sample/haplotype/contig path names and replacement/source names.

## Guard Policy

No downstream task authorizes hidden filtering, compression gates, graph-quality
guards, or score-based replacement filtering. A candidate may be skipped or
failed only for hard path-corruption risks, such as ambiguous or unmappable trim
boundaries, inconsistent orientation metadata, missing traversal interval data,
or observed path-spelling corruption after a trial replacement.

## Validation Expectations

The implementation task now requires synthetic correctness tests for flank-only
alignment, repeated motifs, short/path-end flanks, orientation, trim exclusion,
path-name stability, exact path-spelling preservation, and hard
path-corruption rejection.

The C4 validation task now requires durable measurement artifacts: a report
under `docs/evaluations/`, command lines, input/output graph paths, metrics
tables, exact path-preservation and candidate-accounting results, rendered PNG
URLs, and GFA/zst paths when produced. It also requires a synthetic or
small-locus sanity run in addition to C4.
