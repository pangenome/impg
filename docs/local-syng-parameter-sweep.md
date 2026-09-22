# Local syng parameter sweep

Date: 2026-05-29

Branch: `wg/agent-287/evaluate-local-syng`

Output directory:

`/home/erikg/impg/data/local_syng_parameter_sweep_20260529T193452Z`

## Question

This run tested whether C4 residual fragmentation is caused by the global
HPRCv2 syng index using the paper/default `k=63,s=8` syncmer graph. The
alternate path implemented here is explicit: `gfa:syng-local` extracts the
query-selected sequences, builds a fresh regional syng index with requested
local syncmer parameters, renders it, then applies the normal graph transforms
such as `:crush` and sort. `gfa:syng` semantics are unchanged.

`syng-local` now uses the same default conversion mask as `gfa:syng`: top
0.05% local syncmer frequency filtering, rare-repeat context splitting, and a
short consecutive shared-syncmer run requirement. Add `:nomask` only when the
raw local syng topology is the experiment. `k/s/seed` remain local rebuild
parameters for `syng-local` and assertions against the loaded global index for
`gfa:syng`.

## Test Matrix

Input index and sequences:

- Syng index: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng`
- AGC sequences: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc`
- BED:
  - `GRCh38#0#chr6:31982056-32035418` as `C4_GRCh38_53kb`
  - `CHM13#0#chr6:50000000-50100000` as `boring_chr6_100kb`

Main local sweep:

```bash
impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -b /home/erikg/impg/data/local_syng_parameter_sweep_20260529T193452Z/loci.bed \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o 'gfa:syng-local:blunt,k=<K>,s=<S>,seed=7:crush' \
  -O <run>/graphs \
  --render-graph --render-graph-output <run>/renders \
  --describe-graph --graph-report-output <run>/reports \
  --graph-report-format tsv \
  -t 16 -v 1
```

Local parameter sets:

- `k=63,s=8,seed=7`
- `k=127,s=16,seed=7`
- `k=191,s=24,seed=7`
- `k=311,s=31,seed=7`

Controls:

- `gfa:syng:mask:crush`
- `gfa:pggb`

The PGGB control completed C4, then failed on `boring_chr6_100kb` with:

```text
smooth block 16 produced no GFA; refusing passthrough fallback
```

The C4 PGGB output is included below. The failed lower-locus PGGB run is kept
in `pggb_control/stderr.log` and `pggb_control/time.txt`.

## Implementation

Code changes:

- Added `GfaEngine::SyngLocal`.
- Added `write_syng_region_gfa_from_sequences_with_params`, which builds a
  fresh regional `SyngIndex::new(params).build_region_gbwt(...)` from the
  extracted sequence records instead of filtering the global syng GBWT.
- Added `gfa:syng-local` / `gfa:local-syng` parser support.
- Made `k/s/seed` select the local rebuild parameters for `syng-local`, while
  preserving assertion semantics for `syng`.
- `syng-local` originally defaulted unmasked; current behavior is masked by
  default to match `gfa:syng`, with `:nomask` as the explicit raw-topology mode.

Parser/regression coverage added:

- `test_gfa_engine_syng_local_mode_and_rebuild_params`
- `test_gfa_engine_syng_local_rejects_even_total_length`
- `syng_local_rebuild_uses_requested_syncmer_length`

## Metrics

The long-link count is `path_white_space_bridges_ge_threshold` from
`impg graph-report`, i.e. path white-space bridges at the report threshold
(default 1 kb). Runtime and memory are per two-locus command, except PGGB,
whose command exited after C4 and the lower-locus failure.

| run | locus | segments | segment bp | bp-weighted node path coverage | singleton bp | long links >=1 kb | wall | max RSS GiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| syng-local k=63,s=8 | C4 | 9,037 | 261,748 | 193.11 | 23,085 | 15,684 | 6:05.81 | 72.07 |
| syng-local k=127,s=16 | C4 | 8,350 | 299,351 | 167.13 | 47,694 | 16,415 | 6:14.09 | 54.59 |
| syng-local k=191,s=24 | C4 | 8,017 | 364,876 | 136.16 | 93,969 | 22,092 | 6:09.43 | 44.47 |
| syng-local k=311,s=31 | C4 | 7,322 | 487,442 | 101.51 | 119,908 | 69,328 | 5:09.34 | 43.78 |
| global syng mask+crush | C4 | 9,062 | 263,218 | 192.48 | 26,333 | 13,496 | 6:03.19 | 72.30 |
| PGGB control | C4 | 7,170 | 89,342 | 595.08 | 523 | 9,607 | 12:35.56 | 47.61 |
| syng-local k=63,s=8 | boring chr6 | 11,081 | 203,434 | 216.75 | 22,885 | 568 | 6:05.81 | 72.07 |
| syng-local k=127,s=16 | boring chr6 | 8,816 | 267,993 | 163.99 | 50,056 | 11,681 | 6:14.09 | 54.59 |
| syng-local k=191,s=24 | boring chr6 | 8,534 | 327,982 | 130.58 | 74,336 | 21,865 | 6:09.43 | 44.47 |
| syng-local k=311,s=31 | boring chr6 | 7,420 | 468,720 | 92.07 | 135,254 | 39,630 | 5:09.34 | 43.78 |
| global syng mask+crush | boring chr6 | 11,097 | 225,594 | 195.62 | 45,113 | 1,946 | 6:03.19 | 72.30 |
| PGGB control | boring chr6 | n/a | n/a | n/a | n/a | n/a | 12:35.56 | 47.61 |

Source tables:

- `metrics.tsv`
- `path_validation.tsv`

## Path Validation

Path-name validation:

- C4: `syng-local k=127/191/311` and PGGB had identical path-name sets.
- C4: `syng-local k=63` and `global syng mask+crush` had identical path-name
  sets, but one `HG01960#2` coordinate suffix differed from the PGGB/source
  set (`32078636-32158486` vs `32045898-32158486`).
- boring chr6: all completed syng-local and global syng outputs had identical
  path-name sets.

Strict path-spelled DNA validation:

- PGGB C4 spelled 53,165,920 path bp, matching the extracted sequence panel for
  that run.
- The syng GFA outputs spelled shorter paths in S-line DNA. For example,
  unmasked local C4 path bp ranged from 49,478,677 to 50,547,131, and global
  syng C4 spelled 50,662,998 bp.
- Therefore strict source-sequence byte validation is not satisfied by the
  current syng GFA materialization path. `crush` preserves the path strings it
  is given internally, but the syng/blunt materialization is still not a
  source-sequence-preserving GFA renderer.

Follow-up on 2026-05-30: `syng:blunt` and `syng-local:blunt` now bypass the
generic bluntg trimming pass and materialize exact zero-overlap paths directly
from the syng walk plus source gap DNA. `syng:raw` remains the explicit native
overlap graph mode; concatenating raw S-line DNA without interpreting overlaps
is not a sequence-preserving validation mode. Blunt syng output is exact only
when source sequence files are available for non-syncmer spans; otherwise the
documented `N` gap fill remains an explicit non-source-preserving fallback.

## Render URLs

- https://hypervolu.me/~erik/impg/local-syng-sweep-local_nomask_k63_s8-C4_GRCh38_53kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-local_nomask_k63_s8-boring_chr6_100kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-local_nomask_k127_s16-C4_GRCh38_53kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-local_nomask_k127_s16-boring_chr6_100kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-local_nomask_k191_s24-C4_GRCh38_53kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-local_nomask_k191_s24-boring_chr6_100kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-local_nomask_k311_s31-C4_GRCh38_53kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-local_nomask_k311_s31-boring_chr6_100kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-global_mask_crush-C4_GRCh38_53kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-global_mask_crush-boring_chr6_100kb.png
- https://hypervolu.me/~erik/impg/local-syng-sweep-pggb_control-C4_GRCh38_53kb.png

Upload was verified with:

```bash
ssh erik@hypervolu.me 'ls -lh www/impg/local-syng-sweep-*.png'
```

## Interpretation

Larger local syncmers reduce segment count on C4, but the reduction is not a
quality win. As `k` grows from 63 to 311, C4 segment count drops from 9,037 to
7,322, but stored segment bp increases from 261,748 to 487,442, singleton bp
increases from 23,085 to 119,908, bp-weighted path coverage drops from 193 to
102, and long-link count rises from 15,684 to 69,328. The visual result is a
sparser graph with longer uncovered path jumps, not a PGGB-like condensation.

The lower-complexity chr6 locus shows the over-fragmentation more clearly:
`k=311,s=31` lowers the segment count to 7,420, but segment bp grows to 468,720,
singleton bp grows to 135,254, coverage drops to 92, and long-link count grows
to 39,630. The larger syncmer definitions are too sparse for this regional
rendering path.

The best local syng result is effectively the default `k=63,s=8`; it is close
to the direct global syng control on C4 and slightly better than global on the
boring locus for singleton bp and long links. It still does not close the gap
to PGGB on C4: PGGB has fewer segments, far less stored sequence bp, much lower
singleton bp, and much higher bp-weighted coverage.

## Recommendation

Keep `gfa:syng-local` as an explicit experimental engine. Do not replace
`gfa:syng` semantics, and do not promote larger local syncmers as defaults.

Recommended permanent syntax:

- Keep `gfa:syng-local:blunt,k=<K>,s=<S>,seed=<SEED>:crush` for experiments.
- Use `:nomask` explicitly on `syng-local` only for raw-topology experiments.
- Keep `gfa:syng:mask:crush` as the global-index comparison/control syntax.

Recommended follow-up before promotion:

- Fix strict source-sequence path spelling for syng GFA materialization.
- Investigate the PGGB lower-locus failure in `smooth_gfa`: `smooth block 16
  produced no GFA; refusing passthrough fallback`.
