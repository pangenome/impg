# Validate localized polishing on C4 and SVR loci

Date: 2026-06-08

Task: `validate-localized-polishing`

Run directory for generated artifacts:

`/home/erikg/impg/data/validate_localized_polishing_20260608T225223Z`

Committed metric table:

`docs/evaluations/validate-localized-polishing-metrics.tsv`

## Verdict

The integrated `syng-local:localized` production path was exercised on the
real C4 locus, but both full-C4 attempts timed out before emitting a seed GFA.
This is a runtime/scale blocker in local seed induction, not exact path
corruption:

- default local seed settings timed out after `20:04.53`, max RSS
  `47,119,692 kB`, after FastGA emitted `541,784` raw PAF records;
- tuned C4/k311-style seed settings timed out after `15:04.60`, max RSS
  `46,565,888 kB`, after plane-sweep filtering reduced `541,388` inter-sequence
  PAF records to `431,600` mappings.

No localized-polish output GFA was produced in either run, so exact path
preservation for the new integrated output was not testable. The prior current
C4 local-compression control remains exact-path preserving, but it is still not
PGGB/SmoothXG-quality by diagnostic metrics. The current C4 `poa1kb` output is
best classified as **still not solved**: it preserves paths, but it worsens
segment count, segment bp, singleton bp, low-depth bp, path steps, and total
white-space proxy relative to its seed while staying below the PGGB control on
singleton/low-depth and replay/compression metrics.

## C4 Commands

Default integrated localized attempt:

```bash
source ./env.sh >/dev/null
/usr/bin/time -v -o /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/logs/c4.localized.time.txt \
  /usr/bin/timeout 1200s target/debug/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -b /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/c4.bed \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o 'gfa:syng-local:localized,iterations=1,max-chunks-per-iteration=1,max-total-chunks=1,max-total-bp=200k,max-runtime-secs=600,flank=1k,method=poa,polish-rounds=0,debug-dir=/home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/debug:nosort' \
  -O /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/graphs \
  --render-graph \
  --render-graph-output /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/renders \
  --render-graph-mean-depth \
  --describe-graph \
  --graph-report-output /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/reports \
  --graph-report-format tsv \
  -t 16 -v 1
```

Result: exit `124` from `/usr/bin/timeout`. Last log stage:

```text
[local seed] SweepGA/FastGA alignment emitted 541784 unique raw PAF record(s), 541318 inter-sequence record(s) in 65.244s
```

Tuned integrated localized attempt, using the known C4 k311-style seed options:

```bash
source ./env.sh >/dev/null
/usr/bin/time -v -o /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/logs/c4.localized.tuned.time.txt \
  /usr/bin/timeout 900s target/debug/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -b /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/c4.bed \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  --num-mappings 1:many \
  --scaffold-filter 1:many \
  --scaffold-jump 0 \
  --min-match-len 311 \
  --debug-dir /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/debug/c4_tuned_pipeline \
  -o 'gfa:syng-local:localized,iterations=1,max-chunks-per-iteration=1,max-total-chunks=1,max-total-bp=200k,max-runtime-secs=300,flank=1k,method=poa,polish-rounds=0,debug-dir=/home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/debug/c4_tuned_localized:nosort' \
  -O /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/graphs_tuned \
  --render-graph \
  --render-graph-output /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/renders_tuned \
  --render-graph-mean-depth \
  --describe-graph \
  --graph-report-output /home/erikg/impg/data/validate_localized_polishing_20260608T225223Z/reports_tuned \
  --graph-report-format tsv \
  -t 16 -v 1
```

Result: exit `124` from `/usr/bin/timeout`. Last log stage:

```text
[sweepga] Summary: 541388 -> 431600 mappings (79.7% kept)
```

## C4 Comparison

The comparison rows use current local files already available under
`/home/erikg/impg/data/validate_current_c4_20260608T125053Z` and
`/home/erikg/impg/data/c4_pggb_control_20260526T025439Z`.

| run | status | segments | segment bp | replay/compression | singleton bp | low-depth bp | white-space bridges | white-space bp total | self-loop runs | runtime |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| C4 k311 seed | pass | 7,346 | 253,140 | 421.464 | 5,925 | 52,671 | 5,958 | 24,352,099 | 0 | 3:18.19 |
| C4 poa1kb current | pass | 8,142 | 258,628 | 412.521 | 7,492 | 56,791 | 4,900 | 26,794,743 | 0 | 0:33.07 |
| C4 PGGB/SmoothXG control | pass | 13,288 | 234,524 | 454.919 | 2,890 | 37,477 | 22,851 | 17,165,717 | 15 | 13:37.90 |

Interpretation:

- Current `poa1kb` is exact-path preserving but does not improve over its k311
  seed by the diagnostic metrics above. It increases segments, segment bp,
  singleton bp, low-depth bp, path steps, and total white-space proxy.
- PGGB/SmoothXG has more segments and more long white-space bridges by count,
  but stores fewer segment bp, has much lower singleton and low-depth bp, and
  has a higher replay/compression ratio. It remains the better C4-quality
  target overall.
- The integrated localized path has not yet reached an output graph on full C4
  within the tested budgets, so it is not approaching PGGB/SmoothXG on real C4
  yet.

Exact path preservation:

- Current C4 seed -> `poa1kb`: `expected_paths=465`, `observed_paths=465`,
  `missing_paths=0`, `extra_paths=0`, `spelling_mismatches=0`.
- PGGB control is not exact-name comparable to the seed in this validation:
  `missing_paths=465`, `extra_paths=465`, `spelling_mismatches=0`.
- No exact path corruption was observed in an integrated localized output,
  because no integrated localized output was emitted.

## Non-C4 SVR Panel

Local non-C4 SVR material was available from the existing COSIGT/HGSVC/SVR run:

`/home/erikg/impg/data/cosigt_svr_hprcv2_crush_20260521T114027Z`

Checked inputs:

- BED panel: `cosigt_svr_hprcv2_bounded5.bed`
- GFAs: `graphs_biwfa/AMY1A.gfa`, `CYP2D7.gfa`, `HLA-A.gfa`, `KIR2DL1.gfa`
- Existing renders: `renders_biwfa_1d/*.png` and `renders_biwfa_2d/*.png`

I ran current-schema `impg graph-report --format tsv` diagnostics on all four
non-C4 SVR GFAs. These are diagnostics over existing local SVR graph artifacts,
not integrated `syng-local:localized` production runs, because the full C4
production path had already demonstrated seed-induction timeout at the same
panel scale.

| locus | segments | segment bp | replay/compression | singleton bp | low-depth bp | white-space bridges | white-space bp total | self-loop runs | runtime |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| AMY1A | 61,605 | 1,154,327 | 236.897 | 242,351 | 580,843 | 313,437 | 294,595,633 | 38 | 3:21.13 |
| CYP2D7 | 24,289 | 391,057 | 242.604 | 60,306 | 183,874 | 25,496 | 89,705,883 | 0 | 1:14.03 |
| HLA-A | 74,012 | 1,227,678 | 104.214 | 130,519 | 907,951 | 66,416 | 445,626,460 | 4 | 2:25.82 |
| KIR2DL1 | 52,928 | 822,137 | 165.220 | 130,694 | 463,241 | 243,167 | 210,805,297 | 1 | 2:13.65 |

This panel confirms the diagnostic signals are not C4-only: all four non-C4
SVR graphs have substantial low-depth bp and white-space proxy, and AMY1A/KIR2DL1
also expose self-loop/repeat-run signals.

## Controls And Visualizations

PGGB/SmoothXG C4 control:

`/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.Ygs.gfa`

Available visualization URLs from the current C4 validation upload:

- `http://hypervolu.me/~erik/impg/c4.validate-current-c4.k311.seed.repro.Ygs.png`
- `http://hypervolu.me/~erik/impg/c4.validate-current-c4.k311.poa1kb.repro.Ygs.png`
- `http://hypervolu.me/~erik/impg/c4.validate-current-c4.pggb-control.Ygs.png`

Local non-C4 visualization paths:

- `/home/erikg/impg/data/cosigt_svr_hprcv2_crush_20260521T114027Z/renders_biwfa_1d/AMY1A.1d.png`
- `/home/erikg/impg/data/cosigt_svr_hprcv2_crush_20260521T114027Z/renders_biwfa_1d/CYP2D7.1d.png`
- `/home/erikg/impg/data/cosigt_svr_hprcv2_crush_20260521T114027Z/renders_biwfa_1d/HLA-A.1d.png`
- `/home/erikg/impg/data/cosigt_svr_hprcv2_crush_20260521T114027Z/renders_biwfa_1d/KIR2DL1.1d.png`

## Artifact Policy

Generated GFAs, PNGs, and logs remain outside the repository under
`/home/erikg/impg/data/validate_localized_polishing_20260608T225223Z` and the
pre-existing data directories listed above. This commit includes only this
concise report and the small TSV metrics table.
