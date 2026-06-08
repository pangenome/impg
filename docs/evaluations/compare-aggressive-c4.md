# Aggressive C4 motif-window POA vs abPOA

Task: `compare-aggressive-c4`
Date: 2026-06-05

## Summary

Global SPOA/POA is the better resolver for this C4 input. The 10 kbp and 25 kbp motif-window POA runs produced identical path-preserving graphs, reducing direct self-loop edges from 57 to 24 and adjacent same-node path steps from 252,715 to 40,952 without large graph expansion.

abPOA was verified in global mode, but it is not better here. The 10 kbp abPOA run preserved paths and resolved more candidates, but it expanded the graph from 8,160 to 177,970 segments and from 3.37M to 79.72M path steps while leaving 27 direct self-loop edges and 41,956 adjacent same-node path steps. The 25 kbp abPOA variant was stopped after the 10 kbp result showed this expansion pattern; the abPOA render attempt was also stopped after `gfasort -p Ygs` ran for 15:20.98 without producing a sorted GFA.

The run directory is `/home/erikg/impg/data/c4_aggressive_motif_matrix_20260605T210000Z`.

## Setup

Input GFA:

`/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.poa2kb.gfa`

Runner:

`scripts/run-c4-aggressive-motif-matrix.py`

Common settings:

- `--window-mode motif`
- `--candidate-limit 512`
- `--max-iterations 8`
- `--threads 32`
- `--poa-scoring 1,4,4,2,24,1`
- POA windows: 10 kbp and 25 kbp
- abPOA windows: 10 kbp completed; 25 kbp stopped after 89.94s

abPOA mode was made explicit via the wrapper log at:

`/home/erikg/impg/data/c4_aggressive_motif_matrix_20260605T210000Z/reports/abpoa.invocations.tsv`

The recorded abPOA argv starts with:

`-m 0 -M 1 -X 4 -O 4,24 -E 2,1 -r 3 ...`

That confirms global alignment mode (`-m 0`) and the expected two-piece/convex gap penalties.

## Metrics

Metrics are from `graph-report` and `compare_gfa_paths`. Runtime and RSS are from `/usr/bin/time -v`.

| run | path preserving | resolved | bailed | candidates seen | runtime s | max RSS KB | segments | links | path steps | segment bp | singleton bp | direct self-loops | adjacent same-node steps | self-loop repeat runs | bp-weighted cov | path whitespace p99/max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| input | yes |  |  |  |  |  | 8,160 | 10,928 | 3,374,388 | 251,692 | 6,350 | 57 | 252,715 | 28,031 | 423.889015 | 120,462 / 251,134 |
| POA motif 10k | yes | 1,587 | 510 | 2,566 | 178.45 | 1,511,644 | 10,201 | 13,492 | 3,460,617 | 298,971 | 5,347 | 24 | 40,952 | 4,245 | 356.855595 | 153,011 / 298,597 |
| POA motif 25k | yes | 1,587 | 510 | 2,566 | 169.72 | 1,573,816 | 10,201 | 13,492 | 3,460,617 | 298,971 | 5,347 | 24 | 40,952 | 4,245 | 356.855595 | 153,011 / 298,597 |
| abPOA motif 10k | yes | 2,538 | 396 | 3,314 | 1,513.49 | 24,171,868 | 177,970 | 181,256 | 79,717,063 | 283,273 | 4,927 | 27 | 41,956 | 4,669 | 376.631285 | 0 / 282,899 |
| abPOA motif 25k | no output | partial |  |  | 89.94 | 12,822,992 |  |  |  |  |  |  |  |  |  |  |

All completed outputs passed `compare_gfa_paths` against the input:

| output | expected paths | observed paths | missing | extra | spelling mismatches |
|---|---:|---:|---:|---:|---:|
| POA 10k | 465 | 465 | 0 | 0 | 0 |
| POA 25k | 465 | 465 | 0 | 0 | 0 |
| abPOA 10k | 465 | 465 | 0 | 0 | 0 |

## Remaining Loops

We are not crushing everything. POA removes most of the original repeat signal, but 24 direct self-loop edges and 40,952 adjacent same-node path steps remain. abPOA does not solve this residual either; it leaves 27 direct self-loop edges and 41,956 adjacent same-node path steps while creating a much larger graph.

Top remaining POA motifs:

| node | seq | direct loops | adjacent steps | repeat runs | max run | paths | evidence |
|---|---|---:|---:|---:|---:|---:|---|
| 6326 | A | 1 | 8,627 | 465 | 20 | 465 | emitted in crush log |
| 2108 | A | 1 | 8,267 | 465 | 20 | 465 | emitted in crush log |
| 1599 | T | 1 | 7,422 | 501 | 19 | 465 | emitted in crush log |
| 961 | A | 1 | 6,334 | 465 | 17 | 465 | emitted in crush log |
| 672 | A | 1 | 2,568 | 428 | 7 | 428 | emitted in crush log |
| 6871 | T | 1 | 1,959 | 516 | 10 | 460 | not emitted in top-candidate log |

Top remaining abPOA motifs are nearly the same: 6326/A, 2108/A, 1599/T, 961/A, 672/A, 6871/T, 6717/T, 2567/C, 1174/T, and 6674/T dominate the remaining list.

The main reason these persist is not lack of path preservation. The resolver hit `max-iterations` with residual candidates still appearing. In the POA 25k final round, 105 candidates were generated, 27 built, 22 applied, 78 failed or were empty, 21 were below objective in diagnostics, and 5 were overlap-deferred. Many surviving one-base loop nodes were emitted in candidate logs, but the final rounds mostly selected sparse offshoot windows rather than fully eliminating the high-frequency 1 bp repeat nodes.

## Window Size

Larger POA windows did not help. The 10 kbp and 25 kbp POA runs had identical generated/applied counts and identical graph metrics. The 25 kbp run was slightly faster in this execution, but there is no quality gain.

Larger abPOA windows were not worth continuing. The 10 kbp abPOA run already took 25:13.49, used 24,171,868 KB RSS, wrote a 575 MB GFA, and expanded path steps to 79.7M. The 25 kbp abPOA run was stopped after 89.94s once it had entered the same expansion path.

## Renders

POA render completed and was uploaded:

- Local PNG: `/home/erikg/impg/data/c4_aggressive_motif_matrix_20260605T210000Z/renders/c4.k311.poa.motif25k.cand512.iter8.Ygs.gfalook-m.png`
- URL: `http://hypervolu.me/~erik/impg/c4.k311.poa.motif25k.cand512.iter8.Ygs.gfalook-m.png`
- URL check: HTTP 200

abPOA render did not complete:

- Input GFA: `/home/erikg/impg/data/c4_aggressive_motif_matrix_20260605T210000Z/graphs/c4.k311.abpoa.motif10k.cand512.iter8.gfa`
- Requested sorted GFA: `/home/erikg/impg/data/c4_aggressive_motif_matrix_20260605T210000Z/graphs/c4.k311.abpoa.motif10k.cand512.iter8.Ygs.gfa`
- Requested PNG: `/home/erikg/impg/data/c4_aggressive_motif_matrix_20260605T210000Z/renders/c4.k311.abpoa.motif10k.cand512.iter8.Ygs.gfalook-m.png`
- Status: `gfasort -p Ygs` stopped after 15:20.98 with exit 143 and no sorted output.

## Decision

Use global SPOA/POA, not abPOA, for this C4 motif-local condensation path. abPOA global mode preserves paths, but on this workload it behaves like an over-fragmenting local/base-level expander rather than a useful condenser. Larger motif windows do not improve POA quality in this matrix, and larger abPOA windows are not justified by the 10 kbp behavior.

The practical next lever is not "make the window larger." It is candidate/objective control for the persistent 1 bp repeat motifs: either add a targeted high-frequency loop resolver or change final-round candidate prioritization so repeated one-base self-loop nodes are not left behind after sparse offshoot cleanup.
