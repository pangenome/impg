# C4 Repeat-Collapse Filter Sensitivity Sweep

Run root: `/home/erikg/impg/data/c4_repeat_filter_sweep_20260612T1508Z`

Driver: `scripts/c4-graph-quality-scoreboard.py --preset repeat-filter-sweep`

Committed machine-readable outputs:

- `docs/evaluations/c4-graph-quality-repeat-filter-sweep.tsv`
- `docs/evaluations/c4-graph-quality-repeat-filter-sweep-left-edge-recheck.tsv`

The sweep reran `q_seed_sensitive_pggb` as the same-run baseline and compared
only variants derived from that setting. The old historical `PGGB_control` was
not used as an optimization target. The same-run baseline is the comparison
anchor below; the earlier June 11 row remains provenance for the original
candidate but is not mixed into this matrix.

Command used:

```bash
python scripts/c4-graph-quality-scoreboard.py \
  --preset repeat-filter-sweep \
  --out-dir /home/erikg/impg/data/c4_repeat_filter_sweep_20260612T1508Z \
  --threads 16 \
  --timeout 5400 \
  --render-labels q_seed_sensitive_pggb,q_seed_nm11_sf11_pggb,q_seed_nm1many_sf11_pggb,q_seed_nm11_sf11_nojump_pggb,q_seed_minid98_pggb,q_seed_minid99_pggb
```

Identity filtering was run through the first-class `impg query
--min-aln-identity` option. No hidden wrapper PAF filtering was used.

## Verdict

Keep `q_seed_sensitive_pggb` as the current default/reference. The stricter
global SweepGA/scaffold filters did not convincingly reduce the questionable
right-side repeat-collapse shape, and the small reductions they did produce
came with worse graph-quality metrics. The 99% identity filter produced the
lowest bridge count in the matrix, but it was an obvious regression: it created
much more singleton/sparse sequence and visibly broader vertical gaps in the
right-side graph.

Recommended next experiment: do not make global `1:1` scaffold filtering or a
global 99% identity threshold the default. Instead, keep the seed-sensitive
source extraction fix and investigate targeted, occurrence-aware right-side
repeat controls: for example, copy-aware/local repeat gating or right-region
diagnostics that distinguish homologous shared blocks from spurious
many-to-many glue before graph induction.

## Summary Table

| trial | exact | filtered PAF | segments | total bp | bp weighted cov | singleton bp | sparse bp | seg ws frac | seg ws total | p99 | max | bridges | povu sites | left fixed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| q_seed_sensitive_pggb | pass | 538859 | 6557 | 89237 | 598.805350 | 511 | 4615 | 0.055249 | 2297508 | 3 | 33630 | 3124 | 1592 | 10/10 |
| q_seed_nm11_sf11_pggb | pass | 359661 | 6738 | 89315 | 598.282405 | 506 | 4696 | 0.056021 | 2331652 | 2 | 33637 | 3097 | 1598 | 10/10 |
| q_seed_nm1many_sf11_pggb | pass | 360900 | 6620 | 89334 | 598.155159 | 507 | 4692 | 0.055953 | 2329297 | 2 | 33637 | 3558 | 1593 | 10/10 |
| q_seed_nm11_sf11_nojump_pggb | pass | 390880 | 6694 | 89332 | 598.168551 | 514 | 4684 | 0.055903 | 2327183 | 2 | 33619 | 3102 | 1590 | 10/10 |
| q_seed_minid98_pggb | pass | 538728 | 6644 | 89226 | 598.879172 | 507 | 4611 | 0.055177 | 2294233 | 2 | 33619 | 3590 | 1594 | 10/10 |
| q_seed_minid99_pggb | pass | 537914 | 5994 | 93063 | 574.187303 | 2875 | 6980 | 0.078426 | 3401143 | 3 | 33621 | 2948 | 1644 | 10/10 |

All rows had exact path validation `466/466` with zero missing paths, zero
extra paths, and zero spelling mismatches. The detailed ten-path left-edge
recheck is committed in
`docs/evaluations/c4-graph-quality-repeat-filter-sweep-left-edge-recheck.tsv`.
Every variant placed all ten suspect left-edge paths at the beginning of the
graph layout, with first-node offsets from 54 to 60 bp, so none of the stricter
filters reintroduced the prior left-edge truncation.

## Right-Side Repeat-Collapse Inspection

The visual and metric evidence agree that global stricter filtering is not a
clean fix.

`q_seed_nm11_sf11_pggb` applies `--num-mappings 1:1 --scaffold-filter 1:1`.
It removes about one third of the final PAF records relative to baseline
(`538859` to `359661`) and slightly lowers bridge count (`3124` to `3097`),
but the rendered right side still has the same large shared/collapsed block
shape. It also raises segments (`6557` to `6738`), sparse coverage bp
(`4615` to `4696`), segment whitespace fraction (`0.055249` to `0.056021`),
and total segment whitespace (`2297508` to `2331652`). This is the best
strict scaffold-filter diagnostic by bridge count, but not a recommended
default.

`q_seed_nm1many_sf11_pggb` applies `--num-mappings 1:many
--scaffold-filter 1:1`. It also keeps the left-edge fix, but bridge count
worsens to `3558` and whitespace metrics are worse than baseline. This argues
against treating strict scaffold filtering alone as the remedy.

`q_seed_nm11_sf11_nojump_pggb` applies `--num-mappings 1:1
--scaffold-filter 1:1 --scaffold-jump 0`. It improves p99 whitespace from `3`
to `2`, lowers max path whitespace from `33630` to `33619`, and lowers bridge
count only slightly (`3124` to `3102`). The cost is more segments, more total
bp, more singleton bp, more sparse bp, and more total segment whitespace. The
PNG looks only modestly reorganized on the right side, not convincingly
uncollapsed.

`q_seed_minid98_pggb` uses `--min-aln-identity 0.98`. It barely changes PAF
yield (`538859` to `538728`) and slightly improves singleton/sparse/segment
whitespace metrics, but bridge count worsens to `3590`. A 98% global identity
filter is therefore not addressing the right-side collapse signal.

`q_seed_minid99_pggb` uses `--min-aln-identity 0.99`. It has the lowest bridge
count (`2948`), but this is not a good graph: total segment bp jumps to `93063`,
bp-weighted coverage drops to `574.187303`, singleton bp jumps to `2875`,
sparse bp jumps to `6980`, segment whitespace fraction jumps to `0.078426`,
and total segment whitespace jumps to `3401143`. The PNG visibly shows broad
extra vertical gaps/fragmentation on the right side. Treat this as an obvious
identity-threshold regression, not a solution.

## Uploaded PNGs

All rendered PNG uploads exited with status 0.

| trial | role | PNG |
|---|---|---|
| q_seed_sensitive_pggb | same-run baseline | http://hypervolu.me/~erik/impg/c4-graph-quality-q_seed_sensitive_pggb-c4_repeat_filter_sweep_20260612T1508Z.Ygs.gfalook-m.png |
| q_seed_nm11_sf11_pggb | best strict scaffold-filter diagnostic by bridge count | http://hypervolu.me/~erik/impg/c4-graph-quality-q_seed_nm11_sf11_pggb-c4_repeat_filter_sweep_20260612T1508Z.Ygs.gfalook-m.png |
| q_seed_nm1many_sf11_pggb | scaffold-only regression | http://hypervolu.me/~erik/impg/c4-graph-quality-q_seed_nm1many_sf11_pggb-c4_repeat_filter_sweep_20260612T1508Z.Ygs.gfalook-m.png |
| q_seed_nm11_sf11_nojump_pggb | strictest scaffold/no-jump diagnostic | http://hypervolu.me/~erik/impg/c4-graph-quality-q_seed_nm11_sf11_nojump_pggb-c4_repeat_filter_sweep_20260612T1508Z.Ygs.gfalook-m.png |
| q_seed_minid98_pggb | moderate identity filter | http://hypervolu.me/~erik/impg/c4-graph-quality-q_seed_minid98_pggb-c4_repeat_filter_sweep_20260612T1508Z.Ygs.gfalook-m.png |
| q_seed_minid99_pggb | obvious identity-threshold regression | http://hypervolu.me/~erik/impg/c4-graph-quality-q_seed_minid99_pggb-c4_repeat_filter_sweep_20260612T1508Z.Ygs.gfalook-m.png |

## Validation Notes

- Baseline plus five stricter/identity variants were evaluated.
- Every variant passed exact path validation and the ten-path left-edge
  recheck.
- PNGs were uploaded for all variants.
- Heavy generated GFA, PAF, PNG, and log artifacts remain under the run root
  and are not committed.
