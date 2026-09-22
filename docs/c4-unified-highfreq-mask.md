# C4 unified high-frequency run/span mask

Date: 2026-05-31

Task: `unify-high-frequency`

Output directory:

```text
/home/erikg/impg/data/c4_unified_highfreq_mask
```

## Policy

The high-frequency mask now uses one occurrence-level rule for both explicit
frequency selections and spectrum-selected dispersed scaffold-glue nodes:

```text
private split = selected high-frequency source AND not rescued by freq-run/freq-span
```

Selection sources:

- explicit frequency mask: `top=` and/or `max-occ=`;
- spectrum glue: dispersed high-copy local syncmers selected from the shared
  node spectrum.

Rescue knobs:

- `freq-run` / `high-freq-run`: minimum supported run length in syncmers;
- `freq-span` / `high-freq-span`: minimum repeated exact sequence span in bp.

The defaults are `freq-run=10` and `freq-span=1000`. The older
`min-run` / `sequence-k` knobs still control the generic non-spectrum shared
context filter. `legacy-freq-mask=true` remains available for old node-level
frequency removal. Crush replacement acceptance remains path-preservation based;
compression and white-space metrics are logged as diagnostics only.

## Command

The best row used the exact engine string written to
`unified_top001_run10_span1k/engine.txt`:

```text
gfa:syng:mask,top=0.001,freq-run=10,freq-span=1000,freq-run-aware=true:crush,method=auto,auto-2tier=true,auto-poasta-max-len=10k,max-traversal-len=10k,max-median-traversal-len=1k,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-method=poasta,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort
```

Full query command:

```bash
/home/erikg/impg/target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,top=0.001,freq-run=10,freq-span=1000,freq-run-aware=true:crush,method=auto,auto-2tier=true,auto-poasta-max-len=10k,max-traversal-len=10k,max-median-traversal-len=1k,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-method=poasta,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O /home/erikg/impg/data/c4_unified_highfreq_mask/unified_top001_run10_span1k/run.nosort \
  -v 1
```

The sweep helper invocation was:

```bash
python3 scripts/c4-highfreq-mask-crush-sweep.py \
  --out-dir /home/erikg/impg/data/c4_unified_highfreq_mask \
  --impg /home/erikg/impg/target/release/impg \
  --threads 32 \
  --configs unified_top001_run10_span1k,unified_maxocc1000_run10_span1k \
  --render-top 1
```

## Results

The unified top-fraction row is the best new row. It improves graph
condensation over the previous best `maxocc1000_run10_span1k`, but still does
not reach the PGGB control graph.

| graph | segments | links | total bp | bp-weighted depth | singleton bp | TV | EMD | bp ratio vs PGGB | score |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PGGB target | 7,170 | 8,664 | 89,342 | 595.083 | 523 | - | - | 1.00x | - |
| previous best `maxocc1000_run10_span1k` | 24,168 | 29,224 | 691,137 | 154.368 | 205,569 | 0.677260 | 0.298795 | 7.74x | 1,580.890 |
| new `unified_maxocc1000_run10_span1k` | 24,168 | 29,224 | 691,137 | 154.368 | 205,569 | 0.677260 | 0.298795 | 7.74x | 1,580.890 |
| new best `unified_top001_run10_span1k` | 19,686 | 22,855 | 585,021 | 182.369 | 97,547 | 0.628427 | 0.280768 | 6.55x | 1,481.414 |

Mask counters for the new best row:

| source | selected nodes | occurrences | rescued | run-supported | span-supported | private split |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| explicit HF (`top=0.001`) | 17 | 5,191 | 5,186 | 20 | 5,182 | 5 |
| spectrum glue | 401 | 372,905 | 372,901 | 365,676 | 372,727 | 4 |
| total scaffold context | - | 4,585,120 shared | - | 4,212,124 generic scaffold-supported | 0 generic sequence-supported | 95 |

The important change is that the spectrum-glue path no longer blindly
private-splits the 372,905 selected occurrences. Under `freq-run=10` and
`freq-span=1000`, 372,901 of them are rescued and remain shared.

The absolute-cap row is unchanged relative to the previous best because the
explicit `max-occ=1000` source had already removed most dispersed glue before
spectrum selection; only 3 spectrum-glue nodes remained and none had local
run/span support.

## Validation

Code validation:

- `rustfmt --check src/commands/syng2gfa.rs src/main.rs` passed.
- `cargo fmt --check` still reports pre-existing formatting diffs inside the
  initialized `vendor/gfaffix/src/main.rs` submodule. That vendor file was not
  modified for this task.
- `cargo test --lib syng2gfa::tests -- --nocapture` passed: 24 syng2gfa tests,
  including explicit HF split/rescue and spectrum-glue split/rescue cases.
- `cargo test --bin impg test_gfa_output_format_accepts_syng -- --nocapture`
  passed: 8 parser tests.
- `cargo build --release --bin impg` passed.
- `cargo test` passed: 346 library tests, 78 main tests, integration tests, and
  doctests; the repository's annotated ignored C4 tests remained ignored.

C4 validation:

- Both rows exited with status 0 and graph-report succeeded.
- Both rows have 465 paths and `path_count_ok=true`.
- Best-row sorted-vs-nosort path spelling validation:

  ```text
  expected_paths	465
  observed_paths	465
  missing_paths	0
  extra_paths	0
  spelling_mismatches	0
  ```

- Best-row AGC source validation:

  ```text
  paths	465
  forward_matches	465
  reverse_complement_matches	0
  unparsable_paths	0
  spelling_mismatches	0
  ```

- POA/POASTA/aggressive crush pass-through behavior is intact. The C4 logs show
  replacement build progress entries ending in `(accepted)`, and compression
  diagnostics report `threshold=disabled (diagnostic only); no alternate
  aligner attempted and no replacement was rejected by ratio`.

## PNG

Best-row sorted Ygs `gfalook -m` render:

- local:
  `/home/erikg/impg/data/c4_unified_highfreq_mask/unified_top001_run10_span1k/unified_top001_run10_span1k.Ygs.gfalook-m.png`
- uploaded:
  https://hypervolu.me/~erik/impg/c4-highfreq-mask-crush-sweep-unified_top001_run10_span1k.png

The remote upload was verified with:

```text
-rw-r--r-- 1 erik erik 2.2M May 31 21:39 www/impg/c4-highfreq-mask-crush-sweep-unified_top001_run10_span1k.png
```
