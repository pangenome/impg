# C4 Graph-Quality Optimization Loop

Run root: `/home/erikg/impg/data/c4_graph_quality_20260611T1323Z`

Corrected target/control used: `/home/erikg/impg/data/c4_control_provenance_20260610T154843Z/graphs/current_53kb_pggb/C4_GRCh38_53kb.gfa`.
The historical `PGGB_control` was not used as a gold target.

Current best PNG:
`http://hypervolu.me/~erik/impg/c4-graph-quality-q_seed_sensitive_pggb-20260611T1323Z.Ygs.gfalook-m.png`

Committed machine-readable artifacts:

- `docs/evaluations/c4-graph-quality-scoreboard.tsv`
- `docs/evaluations/c4-graph-quality-left-edge-recheck.tsv`
- `scripts/c4-graph-quality-scoreboard.py`

## Verdict

Best candidate: `q_seed_sensitive_pggb`.

This is still PGGB-style construction, but with query extraction made more sensitive at the left edge by disabling high-copy seed dropping and requiring three seed-walk anchors:

```bash
impg query ... --syng-seed-drop-top-fraction 0 --syng-seed-walk-anchors 3 -o gfa:pggb
```

It preserves exact path spelling (`466/466`, `0` missing, `0` extra, `0` spelling mismatches), fixes all ten diagnosed left-edge paths, and improves most current-control target metrics:

| trial | exact | segments | total_segment_bp | singleton_bp | sparse_bp | ws_frac | ws_total | p99 | max | bridges | povu_sites | left fixed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| control_current_pggb | pass | 7252 | 89353 | 525 | 4691 | 0.063563 | 2646685 | 3 | 33593 | 3102 | 1618 | 0/10 |
| q_seed_sensitive_pggb | pass | 7214 | 89215 | 513 | 4601 | 0.055115 | 2291342 | 3 | 33621 | 4043 | 1591 | 10/10 |
| q_boundary_wide_pggb | pass | 4780 | 57962 | 403 | 5034 | 0.089077 | 2405994 | 12 | 33607 | 6279 | 1248 | 4/10 |
| crush_coverage_repeat_aware | pass | 7707 | 90711 | 584 | 7735 | 0.098803 | 4176552 | 65 | 36492 | 12319 | 1665 | 0/10 |
| crush_motif_local_right | pass | 7252 | 89353 | 525 | 4691 | 0.063563 | 2646685 | 5 | 33593 | 3102 | 1618 | 0/10 |
| crush_outward_guarded | pass | 6054 | 92194 | 2223 | 6285 | 0.078853 | 3387716 | 20 | 33593 | 3571 | 1615 | 0/10 |

`q_seed_sensitive_pggb` is therefore better than the current control on segment count, total segment bp, singleton bp, sparse-coverage bp, segment whitespace fraction/total, and POVU site counts. It is worse on `path_white_space_bridges_ge_threshold` and trivially worse on max path whitespace (`33621` vs `33593`), so the right-side repeated/shared artifact is not fully solved.

## Left-Edge Recheck

The left-edge issue was a source extraction artifact, not a rendering artifact. The seed-sensitive query trial changed the extracted source spellings for all ten suspect paths and moved their first graph node to the beginning of the layout (`54-60 bp` instead of about `33 kb`). The missing left side is now present in the extracted source spellings and in the graph layout.

| suspect path | q_seed_sensitive extracted path | start delta | first graph bp | status |
|---|---|---:|---:|---|
| HG01952#2#CM087999.1:31995790-32075638 | HG01952#2#CM087999.1:31963104-32075581 | -32686 | 54 | left side present |
| HG01960#2#JBHIHN010000011.1:32078636-32158486 | HG01960#2#JBHIHN010000011.1:32045953-32158429 | -32683 | 58 | left side present |
| HG02027#2#JBHDTY010000008.1:5188422-5268271 | HG02027#2#JBHDTY010000008.1:5155736-5268214 | -32686 | 60 | left side present |
| HG02132#2#CM086908.1:31997546-32077395 | HG02132#2#CM086908.1:31964860-32077338 | -32686 | 60 | left side present |
| HG02178#2#CM089938.1:31943310-32055920 | HG02178#2#CM089938.1:31910627-32055863 | -32683 | 54 | left side present |
| HG03710#2#CM086882.1:31974825-32061063 | HG03710#2#CM086882.1:31942142-32061006 | -32683 | 57 | left side present |
| HG03874#1#CM089432.1:31915349-32027959 | HG03874#1#CM089432.1:31882666-32027902 | -32683 | 57 | left side present |
| NA18565#2#CM094342.1:31929266-32009115 | NA18565#2#CM094342.1:31896580-32009058 | -32686 | 60 | left side present |
| NA18945#2#CM101630.1:32018016-32097865 | NA18945#2#CM101630.1:31985330-32097808 | -32686 | 55 | left side present |
| NA18960#2#CM101559.1:31867170-31947019 | NA18960#2#CM101559.1:31834484-31946962 | -32686 | 60 | left side present |

The direct wide-boundary trial only fixed four of ten and shortened the overall graph substantially, so it is not accepted.

## Right-Side Repeated/Shared Region

Three occurrence-aware post-control smoothing/construction variants were tested against the diagnosed right-side repeated/shared region while preserving exact path spelling:

- `crush_coverage_repeat_aware`: 207 replacements across 2 rounds, exact spelling preserved, but severe metric regression (`ws_frac 0.098803`, `p99 65`, `bridges 12319`).
- `crush_motif_local_right`: no candidates generated, exact spelling preserved, effectively a no-op.
- `crush_outward_guarded`: 6 replacements from 32 candidates, exact spelling preserved, reduced segment count (`6054`) but regressed singleton/sparse/whitespace metrics (`singleton_bp 2223`, `sparse_bp 6285`, `ws_frac 0.078853`).

Conclusion: post-control occurrence-aware smoothing did not beat the seed-sensitive PGGB-style query candidate. The best graph fixes source extraction and remains path-correct, but the right-side repeated/shared condensation still needs more selective occurrence-aware construction or replacement gating before it should be accepted as solved.

## Notes

- All scored trials passed exact path validation with zero spelling mismatches.
- Heavy generated GFA, PAF, PNG, and log artifacts remain under the run root and are not committed.
- A broad outward stress run was stopped after it held a large FastGA replacement alignment without producing a complete scored row; the committed scoreboard uses the bounded `crush_outward_guarded` variant instead.
