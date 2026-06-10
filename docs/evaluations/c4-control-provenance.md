# C4 Control Provenance

Task: `rebuild-and-verify`
Date: 2026-06-10 UTC

## Decision

The old `PGGB_control` is historical only. It is not an external PGGB gold
graph, and it is not a valid C4 graph-quality optimization target.

The corrected target for the next C4 graph-quality optimization loop is:

```text
/home/erikg/impg/data/c4_control_provenance_20260610T154843Z/graphs/current_53kb_pggb/C4_GRCh38_53kb.gfa
```

This is the current same-range `impg query -o gfa:pggb` output, validated and
rendered in this task. It should be called the current impg PGGB-style control
candidate, not external PGGB. The current `gfa:syng-local` output is retained
as a current SYNG-local reference/control, but not as the optimization target.

The target metric row is committed in
`docs/evaluations/c4-control-provenance-metrics.tsv`.

## Current Run Scope

Run root:

```text
/home/erikg/impg/data/c4_control_provenance_20260610T154843Z
```

Current code and binary provenance:

| Item | Value |
|---|---|
| source checkout | `/home/erikg/impg` |
| git commit | `f9da204099b79261630f11b3c7be0429acc729a1` |
| impg binary | `/home/erikg/impg/target/release/impg` |
| impg version | `impg 0.4.1` |
| HPRCv2 syng | `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng` |
| HPRCv2 AGC | `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc` |
| BED | `/home/erikg/impg/data/c4_control_provenance_20260610T154843Z/c4_53kb.bed` |
| BED row | `GRCh38#0#chr6 31982056 32035418 C4_GRCh38_53kb` |

Generated GFA, PNG, PAF, FASTA, validation GFA, and log artifacts remain under
the run root and are not committed.

## Current Commands

### Current impg PGGB-style control candidate

Query command, from `logs/current_53kb_pggb.query.time.txt`:

```bash
/usr/bin/timeout 2400s target/release/impg query -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng -b /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/c4_53kb.bed -d 100k --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc --debug-dir /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/debug/current_53kb_pggb_query -o gfa:pggb -O /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/graphs/current_53kb_pggb -t 32 -v 1
```

The command ran from `/home/erikg/impg`, so `target/release/impg` resolves to
`/home/erikg/impg/target/release/impg`. Exit status was 0, wall time was
`9:06.62`, and max RSS was `59,699,612 KB`.

Post-query artifact commands:

```bash
gfasort -i /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/graphs/current_53kb_pggb/C4_GRCh38_53kb.gfa -o /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/sorted/current_53kb_pggb.C4_GRCh38_53kb.Ygs.gfa -p Ygs -t 32 -v 1
/home/erikg/impg/target/release/impg graph-report -g /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/sorted/current_53kb_pggb.C4_GRCh38_53kb.Ygs.gfa -o /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/reports/current_53kb_pggb.C4_GRCh38_53kb.Ygs.graph-report.tsv --format tsv --povu --top 20 -t 32 -v 1
/home/erikg/impg/target/release/examples/compare_gfa_paths /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/validation/current_53kb_pggb.source.gfa /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/sorted/current_53kb_pggb.C4_GRCh38_53kb.Ygs.gfa
gfalook -i /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/sorted/current_53kb_pggb.C4_GRCh38_53kb.Ygs.gfa -o /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/renders/current_53kb_pggb.C4_GRCh38_53kb.Ygs.gfalook-m.png -m -x 3200 -y 1800 -a 3 -t 32 -v 1
scp -o BatchMode=yes -o ConnectTimeout=20 /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/renders/current_53kb_pggb.C4_GRCh38_53kb.Ygs.gfalook-m.png erik@hypervolu.me:www/impg/c4-control-provenance-current-53kb-pggb-20260610T154843Z.Ygs.gfalook-m.png
```

The extracted source spelling table has 467 lines including the header, so the
source graph has 466 expected paths.

### Current SYNG-local reference

Query command, from `logs/current_53kb_syng_local.query.time.txt`:

```bash
/usr/bin/timeout 2400s /home/erikg/impg/target/release/impg query -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng -b /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/c4_53kb.bed -d 100k --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc --debug-dir /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/debug/current_53kb_syng_local_query -o gfa:syng-local -O /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/graphs/current_53kb_syng_local -t 32 -v 1
```

Exit status was 0, wall time was `6:40.50`, and max RSS was `59,008,976 KB`.
The run collected 6,617 local sequences, 26,195,409 bp, 2,128,665 raw PAF
records, and 1,289,044 filtered PAF records. The debug `path_spellings.tsv`
has 6,618 lines including the header.

Post-query artifact commands:

```bash
gfasort -i /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/graphs/current_53kb_syng_local/C4_GRCh38_53kb.gfa -o /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/sorted/current_53kb_syng_local.C4_GRCh38_53kb.Ygs.gfa -p Ygs -t 32 -v 1
/home/erikg/impg/target/release/impg graph-report -g /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/sorted/current_53kb_syng_local.C4_GRCh38_53kb.Ygs.gfa -o /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/reports/current_53kb_syng_local.C4_GRCh38_53kb.Ygs.graph-report.tsv --format tsv --povu --top 20 -t 32 -v 1
/home/erikg/impg/target/release/examples/compare_gfa_paths /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/validation/current_53kb_syng_local.source.gfa /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/sorted/current_53kb_syng_local.C4_GRCh38_53kb.Ygs.gfa
gfalook -i /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/sorted/current_53kb_syng_local.C4_GRCh38_53kb.Ygs.gfa -o /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/renders/current_53kb_syng_local.C4_GRCh38_53kb.Ygs.gfalook-m.png -m -x 3200 -y 1800 -a 3 -t 32 -v 1
scp -o BatchMode=yes -o ConnectTimeout=20 /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/renders/current_53kb_syng_local.C4_GRCh38_53kb.Ygs.gfalook-m.png erik@hypervolu.me:www/impg/c4-control-provenance-current-53kb-syng-local-20260610T154843Z.Ygs.gfalook-m.png
```

## Current Metrics

These rows come from the current graph-report TSV files under the run root.

| Metric | current_53kb_pggb | current_53kb_syng_local |
|---|---:|---:|
| status | REVIEW | REVIEW |
| paths | 466 | 6,617 |
| segments | 7,252 | 5,594 |
| links | 8,751 | 2,397 |
| path_steps | 3,575,249 | 858,019 |
| total_segment_bp | 89,353 | 6,522,327 |
| node_coverage_bp_weighted_mean | 594.643515 | 4.016267 |
| singleton_bp | 525 | 6,483,157 |
| components | 1 | 3,803 |
| largest_component_frac | 1.000000 | 0.318734 |
| path_white_space_bp_p99 | 3 | 18 |
| path_white_space_bp_max | 33,593 | 33,926 |
| path_white_space_bridges_ge_threshold | 3,102 | 703 |
| segment_occupancy_bp_fraction | 0.936437 | 0.000607 |
| segment_white_space_bp_fraction | 0.063563 | 0.999393 |
| segment_white_space_bp_total | 2,646,685 | 43,132,054,664 |
| sparse_coverage_segment_bp | 4,691 | 6,522,327 |
| path_depth_median | 456 | 1 |
| path_depth_p95 | 933 | 929 |
| path_depth_max | 1,866 | 2,700 |
| duplicate_sequence_frac | 0.820877 | 0.864498 |
| local_repeat_context_nodes | 355 | 94 |
| local_repeat_context_occurrences | 559 | 154 |
| direct_self_loop_edges | 0 | 0 |
| adjacent_same_node_path_steps | 0 | 0 |
| adjacent_same_step_path_steps | 0 | 0 |
| self_loop_repeat_runs | 0 | 0 |
| povu_sites | 1,618 | 0 |
| povu_leaf_sites | 1,603 | 0 |

Interpretation:

- `current_53kb_pggb` is the usable target vector for the next optimization
  loop because it uses the current code, same 53 kb BED, HPRCv2 syng/AGC input,
  466 source paths, and exact path validation.
- `current_53kb_syng_local` is useful as a current reference for the SYNG-local
  query path, but it is not a target: it has 6,617 fragment/source paths,
  3,803 components, 6,483,157 singleton bp, and 0.999393 segment whitespace
  fraction.
- Low `path_white_space_bridges_ge_threshold` is not sufficient by itself.
  Future scoring should use the metric vector above plus the visual renders.

## Path Validation

| Control | Expected paths | Observed paths | Missing | Extra | Spelling mismatches | Log |
|---|---:|---:|---:|---:|---:|---|
| current_53kb_pggb | 466 | 466 | 0 | 0 | 0 | `/home/erikg/impg/data/c4_control_provenance_20260610T154843Z/validation/current_53kb_pggb.compare_gfa_paths.stdout.log` |
| current_53kb_syng_local | 6,617 | 6,617 | 0 | 0 | 0 | `/home/erikg/impg/data/c4_control_provenance_20260610T154843Z/validation/current_53kb_syng_local.compare_gfa_paths.stdout.log` |

## Uploaded Renders

Both URLs returned HTTP 200 with `Content-Type: image/png` on 2026-06-10.

| Control | PNG URL | Content-Length |
|---|---|---:|
| current_53kb_pggb | http://hypervolu.me/~erik/impg/c4-control-provenance-current-53kb-pggb-20260610T154843Z.Ygs.gfalook-m.png | 188,033 |
| current_53kb_syng_local | http://hypervolu.me/~erik/impg/c4-control-provenance-current-53kb-syng-local-20260610T154843Z.Ygs.gfalook-m.png | 1,066,730 |

## Historical Controls

These files remain useful context, but they must be labeled as historical and
must not be used as the active optimization target.

| Historical label | Artifact / report | Provenance | Status | Reason |
|---|---|---|---|---|
| `PGGB_control` | `/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.nosort.gfa`, sorted as `pggb.Ygs.gfa`; documented in `docs/crush-vs-pggb-comparison.md` | `impg query -o 'gfa:pggb'` on 2026-05-26, from `/home/erikg/impg/target/release/impg`, syng query `GRCh38#0#chr6:31891045-32123783`, `-d 50k`, HPRCv2 syng/AGC, 465 path ranges | invalid as target; historical only | It was an internal impg query path, not external PGGB. It used a different larger region than the current 53 kb BED and has 465 paths rather than 466. Older reports also disagree on some whitespace metrics, which is another reason not to promote it as a gold target. |
| `onechunk_baseline` | `/home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z/graphs_tuned/C4_GRCh38_53kb.gfa`; documented in `docs/evaluations/diagnose-exact-c4.md` and `docs/evaluations/diagnose-c4-query-2.md` | Current-at-the-time localized `impg query -o "gfa:syng-local:localized,...,max-chunks-per-iteration=1,...,method=poa,polish-rounds=0"` against `c4.bed`, `-d 100k`, HPRCv2 syng/AGC; exit 0, 466 paths | valid historical diagnostic; not target | It reproduced a one-chunk localized path and exposed query/end-effect issues. It was not a clean target control and its raw/sorted reports have metric drift around whitespace fields. |
| `syng_k311_seed` | Original sweep `/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.seed.gfa`; reproduced row `syng_k311_seed_repro` in `docs/evaluations/validate-current-c4-metrics.tsv` | `impg graph --gfa-engine seqwish --fastga --num-mappings 1:many --scaffold-filter 1:many --scaffold-jump 0 --min-match-len 311` using the C4 whole-region FASTA and historical k79 raw PAF | valid historical seed/reproduction; not target | Useful for k-sweep and seed comparisons, but it is not the same extraction/build path as the current `impg query -o gfa:pggb` control and its reproduced graph has 465 paths. |
| `fix_c4_near` | `/home/erikg/impg/data/fix_c4_near_20260610T105108Z_crush_manymany_large_final/C4_GRCh38_53kb.crush-sweepga-manymany.large.gfa`; report `docs/evaluations/fix-c4-alignment-yield.md` | `impg crush` with `--method sweepga`, one iteration, many:many replacement defaults, input `/home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z/graphs_tuned/C4_GRCh38_53kb.gfa`; exact path validation passed for 466 paths | valid historical alignment-yield proof; not target | It proved replacement alignment and lacing can apply path-preservingly. The report explicitly treats graph-quality movement as diagnostic only, not as a target metric. |

Historical commit/version evidence:

| Historical label | Recoverable commit evidence | Binary/version evidence |
|---|---|---|
| `PGGB_control` | Report commit `ddee2ad9447551979c0ca8babb3877202f2e968c`; task commit `68021d08351dd864d1db2390636a002d9583dcc5`; artifact log timestamps begin 2026-05-26T02:54:52Z | Artifact records `/home/erikg/impg/target/release/impg`; `impg --version` was not recorded in the artifact, so the exact binary version is unknown. |
| `onechunk_baseline` | Code/report lineage commits `f7e41df20414e2ab4b91ab7223c8b3b53504bca1` and `08f12cf7612d01c9b6be2f3f6f7f468ae9cc0082` | Artifact records `target/release/impg`; `impg --version` was not recorded in the artifact. |
| `syng_k311_seed` | Original sweep commits `8bcb58e52897e685cc3b7de7e2964fd910e2781c` and `8861af82cbaa5351fdfb8ca5ad6352c4fd2b7a82`; reproduced validation commits `b1410e522603641e11df779bc592702064578305` and `8b5081df883f5516a9bca7642d2f1ab005764831`; reproduced report branch head `e3c2494` | Original commands record `impg` from PATH; reproduced commands record `target/release/impg`; exact `impg --version` was not recorded in those artifacts. |
| `fix_c4_near` | Code/report lineage commits `b0759d864f52f2127c3dfd7ba1cc2d10fbe0931f` and `f9da204099b79261630f11b3c7be0429acc729a1` | Artifact records `target/release/impg`; current checkout at `f9da204099b79261630f11b3c7be0429acc729a1` reports `impg 0.4.1`. |

Historical `PGGB_control` command, exactly as documented:

```bash
out=/home/erikg/impg/data/c4_pggb_control_20260526T025439Z
/usr/bin/time -v /home/erikg/impg/target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:pggb' \
  -O "$out/pggb.nosort" -v 1 > "$out/pggb.nosort.stdout" 2> "$out/pggb.nosort.stderr"

gfasort -i "$out/pggb.nosort.gfa" -o "$out/pggb.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/pggb.Ygs.gfa" -o "$out/pggb.Ygs.png" -m -x 2200 -y 1200
scp "$out/pggb.Ygs.png" erik@hypervolu.me:www/impg/c4-pggb-control.png
```

Historical `onechunk_baseline` command, exactly as documented:

```bash
source ./env.sh >/dev/null
out=/home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z
IMPG_CRUSH_DEBUG_DIR="$out/debug/internal_tail" \
/usr/bin/time -v -o "$out/logs/internal.release.tuned.time.txt" \
  /usr/bin/timeout 900s target/release/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -b "$out/c4.bed" \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  --num-mappings 1:many \
  --scaffold-filter 1:many \
  --scaffold-jump 0 \
  --min-match-len 311 \
  --debug-dir "$out/debug/c4_tuned_pipeline" \
  -o "gfa:syng-local:localized,iterations=1,max-chunks-per-iteration=1,max-total-chunks=1,max-total-bp=200k,max-runtime-secs=300,flank=1k,method=poa,polish-rounds=0,debug-dir=$out/debug/c4_tuned_localized:nosort" \
  -O "$out/graphs_tuned" \
  --render-graph \
  --render-graph-output "$out/renders_tuned" \
  --render-graph-mean-depth \
  --describe-graph \
  --graph-report-output "$out/reports_tuned" \
  --graph-report-format tsv \
  -t 16 -v 1 \
  > "$out/logs/internal.release.tuned.stdout.log" \
  2> "$out/logs/internal.release.tuned.stderr.log"
printf '%s\n' "$?" > "$out/logs/internal.release.tuned.exitcode"
```

Historical original `syng_k311_seed` command, exactly as recorded in
`/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/commands.md`:

```bash
/usr/bin/time -v -o /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/logs/c4.k311.seed.graph.time.txt impg graph --sequence-files /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa --paf-file /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/debug/k79/raw.paf --gfa-engine seqwish --fastga --num-mappings 1:many --scaffold-filter 1:many --scaffold-jump 0 --min-match-len 311 --temp-dir /tmp/c4lk311 --debug-dir /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/debug/k311 -g /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.seed.gfa -t 32 -v 1
```

Historical reproduced `syng_k311_seed` command, exactly as documented in
`docs/evaluations/validate-current-c4.md`:

```bash
target/release/impg graph \
  --sequence-files /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa \
  --paf-file /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/debug/k79/raw.paf \
  --gfa-engine seqwish \
  --fastga \
  --num-mappings 1:many \
  --scaffold-filter 1:many \
  --scaffold-jump 0 \
  --min-match-len 311 \
  --temp-dir /home/erikg/impg/data/validate_current_c4_latest/tmp/c4.k311.seed \
  --debug-dir /home/erikg/impg/data/validate_current_c4_latest/debug/c4.k311.seed \
  -g /home/erikg/impg/data/validate_current_c4_latest/graphs/c4.k311.seed.repro.gfa \
  -t 32 -v 1
```

Historical `fix_c4_near` command, exactly as documented:

```bash
IMPG_CRUSH_DEBUG_DIR="$out/debug/crush" \
/usr/bin/time -v -o "$out/logs/crush.time.txt" \
  /usr/bin/timeout 600s target/release/impg crush \
    -g /home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z/graphs_tuned/C4_GRCh38_53kb.gfa \
    -o "$out/C4_GRCh38_53kb.crush-sweepga-manymany.large.gfa" \
    --method sweepga \
    --max-iterations 1 \
    --max-span 7k \
    --min-traversal-len 5k \
    --max-traversal-len 7k \
    --max-median-traversal-len 7k \
    --max-total-sequence 2m \
    --max-traversals 10k \
    --polish-rounds 0 \
    -t 16 -v 1
```

## External PGGB Feasibility

Standalone `pggb` and `smoothxg` were not installed in this environment:

```text
command -v pggb     -> not found
command -v smoothxg -> not found
```

Therefore a true external/current PGGB + SmoothXG build was not feasible in
this pass. The extracted 466-sequence FASTA for the current same-range PGGB-style
candidate exists at:

```text
/home/erikg/impg/data/c4_control_provenance_20260610T154843Z/debug/current_53kb_pggb_query/combined.fa
```

If external PGGB/SmoothXG is installed later, build from that FASTA or a
freshly re-extracted same-commit equivalent, then require graph-report,
`compare_gfa_paths`, `gfalook -m`, upload, and this same provenance labeling
before promoting it.

## Superseding c4-graph-quality

The paused `c4-graph-quality` work should be superseded with this corrected
target row:

| Field | Corrected target value |
|---|---:|
| target graph | `/home/erikg/impg/data/c4_control_provenance_20260610T154843Z/graphs/current_53kb_pggb/C4_GRCh38_53kb.gfa` |
| target class | current impg `gfa:pggb` PGGB-style control candidate |
| external PGGB gold? | no |
| exact paths | pass, 466 expected / 466 observed / 0 mismatches |
| segments | 7,252 |
| links | 8,751 |
| paths | 466 |
| path_steps | 3,575,249 |
| total_segment_bp | 89,353 |
| node_coverage_bp_weighted_mean | 594.643515 |
| singleton_bp | 525 |
| sparse_coverage_segment_bp | 4,691 |
| segment_white_space_bp_fraction | 0.063563 |
| segment_white_space_bp_total | 2,646,685 |
| path_white_space_bp_p99 | 3 |
| path_white_space_bp_max | 33,593 |
| path_white_space_bridges_ge_threshold | 3,102 |
| direct_self_loop_edges | 0 |
| adjacent_same_node_path_steps | 0 |
| adjacent_same_step_path_steps | 0 |
| povu_sites | 1,618 |
| povu_leaf_sites | 1,603 |
| visual URL | http://hypervolu.me/~erik/impg/c4-control-provenance-current-53kb-pggb-20260610T154843Z.Ygs.gfalook-m.png |

Future optimization should compare candidates against the full vector above and
the render, not against one old whitespace counter. The old `PGGB_control`,
`onechunk_baseline`, `syng_k311_seed`, and `fix_c4_near` rows can be cited as
historical references only.

## Validation Notes

- Current `gfa:pggb` query exit code: 0.
- Current `gfa:syng-local` query exit code: 0.
- `compare_gfa_paths` passed for both current controls.
- `gfalook -m` PNGs were uploaded and verified by HTTP HEAD.
- Only this Markdown report and the small TSV summary are intended for commit.
- Generated GFA, PNG, FASTA, PAF, validation GFA, and log files remain outside
  git under `/home/erikg/impg/data/c4_control_provenance_20260610T154843Z`.
