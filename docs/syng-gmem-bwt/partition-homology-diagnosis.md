# Yeast partition homology and graph-view diagnosis

This is diagnostic evidence, not a policy change or a validated inference catalog.
The complete source partitioning run remains uninterrupted, without a wall-clock
cutoff. Local PGGB graph construction is an optional diagnostic view, not a new
requirement for the sample MEM-BWT method.

## Inputs and scope

Three saved groups from the static sample64 1200-second snapshot:

| Partition | Owned `S288C#0#chrI` interval | Occurrences | Panel identities |
|---|---:|---:|---:|
| 1 | 10058–20036 | 6 | 6 |
| 10 | 100054–110021 | 256 | 206 |
| 16 | 160031–170031 | 574 | 110 |

All source-forward BED intervals were extracted with their complete source names
and coordinates. The internal PGGB engine used default FastGA alignment plus
SweepGA filtering, seqwish, smoothing and gfaffix. All 836 graph paths reproduce
the extracted source sequences exactly, and all path adjacencies exist. Those
checks are necessary but, as the controls below demonstrate, insufficient for
preserving useful homology representation.

Source commit: `295bca975f9c584cb59919451d11b8a18671fe2f`; frozen executable SHA256:
`6ea56ee21eed2f627efad079a2ac52e38fa284c688e0a96cc31acffd7013c726`.
Jobs were serial, restricted to four CPU cores and four threads.

- [Original GFA Look renderings and graph downloads](http://hypervolu.me/~erik/yeast/partition-graphs-20260910T234952Z.html)
- [Diagnosis, matrices, fixtures, stage measurements and checks](http://hypervolu.me/~erik/yeast/partition-graphs-20260910T234952Z-diagnosis.html)
- Local evidence: `~/yeast/partition-graphs-20260910T234952Z/sharing-diagnosis/`.
- Run scripts: `target/experiments/yeast-partition-renderings/` in the source checkout.

## 1. Reverse-oriented homology is lost during smoothing

The partition10 weighted node-sharing matrix separates 28 occurrences from the
main group despite one weakly connected graph component. These occurrences have
90.6–99.4% source-base coverage by exact ACGT 63mers found in the reference in
either orientation, but only about 2% of their bases use reference graph nodes
in the published PGGB graph.

Stage isolation:

1. Raw FastGA PAF contains 12,736 cross-group alignments, all reverse-strand.
2. SweepGA retains 12,488 of these alignments.
3. The unsmoothed seqwish graph gives the 28 paths 98.4–99.4% shared-reference-node
   coverage and 175,656 reverse-oriented graph bp in total.
4. A direct `smooth_gfa` replay reduces shared-reference-node coverage to
   1.94–2.37%, before gfaffix is invoked. The published full PGGB output similarly
   has 2.01–2.20% sharing, with no reverse graph steps in these 28 paths.
5. Every source sequence remains exact across the direct smoothing comparison.

A minimal orientation control removes discovery and alignment entirely. Start
with one deterministic, N-free 2000bp sequence node and two paths through it:

| Two paths | Shared-node coverage before smoothing | After default smoothing |
|---|---:|---:|
| forward / forward | 100% | 100% |
| forward / reverse | 100% | 1.5% |

Both path spellings remain exact. This is a homology-representation defect, not
sequence loss. `src/smooth.rs::smooth_block` extracts each range in its own
path-traversal direction and feeds those strings into POA without normalizing
opposing block orientations and retaining that orientation through reconstruction.
Consequently true reverse homology becomes largely separate forward-oriented
sequence, potentially joined by incidental short matches.

Do not interpret the published smoothed graph's connectivity or `-z` coloring as
validated evidence of homology or absence of inversions. This smoothing problem
is separate from the already-known stored-path reverse-complement sketch recall
limitation. Sorting/layout can vary between replays; the claims above concern
spelled sequences, alignment direction and node-sharing coverage, not stable IDs
or byte-identical layouts.

## 2. Automatic FastGA frequency selection causes a separate dropout

In partition1, `CBM#1#chrI_1:14075-21029` has no nodes in common with the reference
and occupies its own graph component. However, 91.3% of its source bases are
covered by exact shared 63mers. Its original passing discovery chain contains
221 anchors on 219 distinct syncmer nodes, with 6013bp query seed-union coverage.

The unchanged default replay has only a self-alignment for this sequence in the
raw PAF: the lost connections do not originate in SweepGA filtering or smoothing.

A single controlled counterfactual changes only FastGA's cutoff on the same
six-sequence input:

| Cutoff | CBM raw cross-path alignments | Retained by SweepGA | Seqwish components |
|---|---:|---:|---:|
| automatic 6 | 0 | 0 | 2 |
| diagnostic 7 | 10 | 10 | 1 |

The recovered alignments span the entire 6954bp CBM interval at about 99%
identity. Both graphs still reproduce all six input sequences exactly.

Cutoff7 is not an adopted default, and this experiment does not justify a
blanket cutoff change. The native FastGA self-adaptamer code has strict
frequency comparisons; a panel-identity count is also not generally the same
as a bound on occurrence multiplicity in repetitive inputs.

Partition16 contains a main component of 566 paths and six other components
containing eight paths. Those eight have 85.6–93.6% exact-63mer source coverage
against the reference, but no raw alignments connect their groups to the main
component. Their dropout is localized before filtering; the precise frequency
mechanism has not been separately proven for every one of them.

## 3. Anchor extent and matched coverage are different

The partition wrapper's 50% criterion uses the envelope from the first to last
query anchor, not the union of bases actually covered by exact seeds. The
requested minimum5 is adaptively2 for these 10kb queries, although the observed
retained chains have many more distinct nodes:

| Partition | Matching passing chains | Minimum distinct nodes | Minimum query seed-union bp |
|---|---:|---:|---:|
| 1 | 6 | 219 | 6013 |
| 10 | 256 | 35 | 1205 |
| 16 | 591 | 30 | 1320 |

Multiple passing chains can overlap one retained/merged BED occurrence. All 836
BED occurrences overlap a reproduced passing chain. This is not a full replay
of historical greedy ownership masking or proof of full-span support for each
owned interval.

Example: `ASN#0#chrX:363741-369117` has 30 distinct seed nodes and a5377bp query
extent, satisfying the5000bp extent threshold, but only1320bp of exact query
seed-union support. That does not establish that the interior is unrelated; it
establishes that passing the threshold does not prove50%-covered homology.

The trace directly compared156,645 anchor occurrences with source nucleotides,
with zero mismatches;722 boundary-clipped anchors were not checked. The
single-syncmer-bridge hypothesis is not supported for these inspected groups.

Do not blindly require50% seed-union coverage instead: the reverse-oriented
homologs in partition10 also have sparse stored-sketch support. Validation must
use independent sequence alignment and supported flanks, rather than turn an
observed diagnostic difference into an untested rejection rule.

## Diagnostic metrics and future regression gates

The retained matrix uses node-length-weighted Jaccard on distinct node IDs.
Source-base sharing separately counts all source traversal bp using nodes
present in the reference. Independent exact-63mer coverage ignores N-containing
windows and checks both orientations. It is not percent identity, a collinear
alignment, dosage, or evidence of unique placement; repeats can match multiple
positions. Full matrix row identities and per-occurrence metrics are retained.

Before accepting graph-derived homology evidence:

- Add an orientation regression requiring the exact forward/reverse control to
  retain shared homology, not only its spelled sequences.
- Retain the six-path CBM fixture as an automatic-frequency/alignment-completeness
  regression. The diagnostic scripts currently reproduce the problem; they are
  not assertions that the implementation has been repaired.
- Compare raw versus filtered alignment neighbors and unsmoothed versus smoothed
  sharing, including weakly connected subgroups rather than components alone.
- Track distinct seeds, union coverage, anchor envelope, collinear nucleotide
  coverage and flank support independently.
- Validate shared cores and spanning candidates before treating a computational
  partition as one inference chunk. Preserve repeat occurrences, coordinates,
  continuations and uncertainty rather than dropping outliers to make graphs
  appear cleaner.

No smoothing fix, frequency-policy change, stronger partition threshold or
catalog acceptance is claimed by this diagnosis.
