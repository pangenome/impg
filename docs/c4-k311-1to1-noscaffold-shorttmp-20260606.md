# C4 k311 1:1 no-scaffold short-temp tracking

Task: `track-active-c4`

Run directory:
`/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z`

## Verdict

The short-temp run completed successfully and produced seed, polished, sorted,
report, render, upload, and path-preservation artifacts.

`--num-mappings 1:1 --scaffold-jump 0` partially improves the seed relative to
the existing k311 `1:many` scaffold0 seed: it eliminates direct self loops and
adjacent same-step path repeats, and it sharply lowers worst jump/whitespace
metrics. It does not clearly improve overall graph quality after SPOA2kb
default-selfloop polishing. Compared with the existing verified
SPOA2kb/default-selfloop output, the 1:1 polished graph has nearly identical
jump metrics and self-loop status, but worse whitespace fraction, whitespace
bridge count, singleton bp, sparse coverage bp, and local repeat contexts.

Recommendation: keep this run as a diagnostic artifact, but do not promote the
1:1 no-scaffold configuration as an improvement over the existing
SPOA2kb/default-selfloop polished output.

## Exact Command

The already-running shell command used `RUN_ID=20260606T070148Z`, wrote to
`OUT=/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z`,
and used `/tmp/c4_1to1_20260606T070148Z` as the short temp root.

```bash
set -euo pipefail
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
OUT="/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_${RUN_ID}"
TMP_ROOT="/tmp/c4_1to1_${RUN_ID}"
mkdir -p "$OUT"/graphs "$OUT"/reports "$OUT"/sorted "$OUT"/renders "$OUT"/logs "$OUT"/validation "$OUT"/debug "$TMP_ROOT"
FASTA="/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa"
IMPG="/home/erikg/impg/target/release/impg"
GFASORT="/home/erikg/.cargo/bin/gfasort"
GFALOOK="/home/erikg/.cargo/bin/gfalook"
COMPARE="/home/erikg/impg/target/release/examples/compare_gfa_paths"
SEED="$OUT/graphs/c4.k311.1to1.noscaffold.seed.gfa"
POLISHED="$OUT/graphs/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.gfa"
SEED_SORTED="$OUT/sorted/c4.k311.1to1.noscaffold.seed.Ygs.gfa"
POLISHED_SORTED="$OUT/sorted/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.Ygs.gfa"
SEED_REPORT="$OUT/reports/c4.k311.1to1.noscaffold.seed.Ygs.graph-report.tsv"
POLISHED_REPORT="$OUT/reports/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.Ygs.graph-report.tsv"
SEED_PNG="$OUT/renders/c4.k311.1to1.noscaffold.seed.Ygs.mean-depth.png"
POLISHED_PNG="$OUT/renders/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.Ygs.mean-depth.png"

/usr/bin/time -v "$IMPG" graph --sequence-files "$FASTA" --gfa-engine seqwish --fastga --num-mappings 1:1 --scaffold-jump 0 --min-match-len 311 --temp-dir "$TMP_ROOT" --debug-dir "$OUT/debug/seed" -g "$SEED" -t 32 -v 1 >"$OUT/logs/seed.graph.stdout.log" 2>"$OUT/logs/seed.graph.stderr.log"
/usr/bin/time -v "$IMPG" crush -g "$SEED" -o "$POLISHED" --method poa --max-iterations 5 --max-traversal-len 2000 --max-median-traversal-len 2000 --max-total-sequence 1m --max-traversals 10k -t 32 -v 1 >"$OUT/logs/crush.stdout.log" 2>"$OUT/logs/crush.stderr.log"
/usr/bin/time -v "$COMPARE" "$SEED" "$POLISHED" >"$OUT/validation/compare_gfa_paths.stdout.log" 2>"$OUT/validation/compare_gfa_paths.stderr.log"
/usr/bin/time -v "$GFASORT" -i "$SEED" -o "$SEED_SORTED" -p Ygs -t 32 -v 1 >"$OUT/logs/seed.gfasort.stdout.log" 2>"$OUT/logs/seed.gfasort.stderr.log"
/usr/bin/time -v "$GFASORT" -i "$POLISHED" -o "$POLISHED_SORTED" -p Ygs -t 32 -v 1 >"$OUT/logs/polished.gfasort.stdout.log" 2>"$OUT/logs/polished.gfasort.stderr.log"
/usr/bin/time -v "$IMPG" graph-report -g "$SEED_SORTED" -o "$SEED_REPORT" --format tsv --povu --top 20 -t 32 -v 1 >"$OUT/logs/seed.graph-report.stdout.log" 2>"$OUT/logs/seed.graph-report.stderr.log"
/usr/bin/time -v "$IMPG" graph-report -g "$POLISHED_SORTED" -o "$POLISHED_REPORT" --format tsv --povu --top 20 -t 32 -v 1 >"$OUT/logs/polished.graph-report.stdout.log" 2>"$OUT/logs/polished.graph-report.stderr.log"
/usr/bin/time -v "$GFALOOK" -i "$SEED_SORTED" -o "$SEED_PNG" -m -x 3200 -y 1800 -a 3 -t 32 -v 1 >"$OUT/logs/seed.gfalook.stdout.log" 2>"$OUT/logs/seed.gfalook.stderr.log"
/usr/bin/time -v "$GFALOOK" -i "$POLISHED_SORTED" -o "$POLISHED_PNG" -m -x 3200 -y 1800 -a 3 -t 32 -v 1 >"$OUT/logs/polished.gfalook.stdout.log" 2>"$OUT/logs/polished.gfalook.stderr.log"
/usr/bin/time -v scp -o BatchMode=yes -o ConnectTimeout=20 "$SEED_PNG" erik@hypervolu.me:www/impg/c4-k311-1to1-noscaffold-seed-${RUN_ID}.Ygs.mean-depth.png >"$OUT/logs/seed.upload.stdout.log" 2>"$OUT/logs/seed.upload.stderr.log"
/usr/bin/time -v scp -o BatchMode=yes -o ConnectTimeout=20 "$POLISHED_PNG" erik@hypervolu.me:www/impg/c4-k311-1to1-noscaffold-spoa2k-default-selfloop-${RUN_ID}.Ygs.mean-depth.png >"$OUT/logs/polished.upload.stdout.log" 2>"$OUT/logs/polished.upload.stderr.log"
printf '%s\n' "$OUT" > /home/erikg/impg/data/latest_c4_k311_1to1_noscaffold_shorttmp_dir.txt
{
  printf 'seed\thttp://hypervolu.me/~erik/impg/c4-k311-1to1-noscaffold-seed-%s.Ygs.mean-depth.png\n' "$RUN_ID"
  printf 'spoa2k_default_selfloop\thttp://hypervolu.me/~erik/impg/c4-k311-1to1-noscaffold-spoa2k-default-selfloop-%s.Ygs.mean-depth.png\n' "$RUN_ID"
} > "$OUT/uploaded_urls.tsv"
rm -rf "$TMP_ROOT"
```

## Produced Artifacts

Seed GFA:
`/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/graphs/c4.k311.1to1.noscaffold.seed.gfa`

Polished GFA:
`/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/graphs/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.gfa`

Ygs-sorted seed GFA:
`/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/sorted/c4.k311.1to1.noscaffold.seed.Ygs.gfa`

Ygs-sorted polished GFA:
`/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/sorted/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.Ygs.gfa`

Seed graph-report TSV:
`/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/reports/c4.k311.1to1.noscaffold.seed.Ygs.graph-report.tsv`

Polished graph-report TSV:
`/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/reports/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.Ygs.graph-report.tsv`

Seed PNG:
`/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/renders/c4.k311.1to1.noscaffold.seed.Ygs.mean-depth.png`

Polished PNG:
`/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/renders/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.Ygs.mean-depth.png`

Uploaded PNG URLs:

| artifact | URL |
| --- | --- |
| seed | http://hypervolu.me/~erik/impg/c4-k311-1to1-noscaffold-seed-20260606T070148Z.Ygs.mean-depth.png |
| spoa2k default-selfloop | http://hypervolu.me/~erik/impg/c4-k311-1to1-noscaffold-spoa2k-default-selfloop-20260606T070148Z.Ygs.mean-depth.png |

## Path Preservation

Log paths:

| log | path |
| --- | --- |
| stdout | `/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/validation/compare_gfa_paths.stdout.log` |
| stderr/time | `/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/validation/compare_gfa_paths.stderr.log` |

Summary:

| metric | value |
| --- | ---: |
| expected_paths | 465 |
| observed_paths | 465 |
| missing_paths | 0 |
| extra_paths | 0 |
| spelling_mismatches | 0 |

## Runtime and Exit Status

All short-temp steps exited with status 0. Total summed wall-clock runtime was
361.91 seconds, approximately 6:01.91.

| step | elapsed | max RSS KB | exit |
| --- | ---: | ---: | ---: |
| seed graph | 4:45.65 | 8689464 | 0 |
| crush polish | 0:33.92 | 1705596 | 0 |
| path compare | 0:00.45 | 307200 | 0 |
| seed sort | 0:08.86 | 156868 | 0 |
| polished sort | 0:10.25 | 172176 | 0 |
| seed graph-report | 0:07.77 | 980828 | 0 |
| polished graph-report | 0:08.65 | 1058588 | 0 |
| seed render | 0:01.79 | 112640 | 0 |
| polished render | 0:01.81 | 116736 | 0 |
| seed upload | 0:01.38 | 8192 | 0 |
| polished upload | 0:01.38 | 8192 | 0 |

## Graph-Report Metrics

Comparison baselines:

| label | graph-report TSV |
| --- | --- |
| 1:many scaffold0 seed | `/home/erikg/impg/data/c4_hard_seqwish_k_20260605T135300Z/reports/one_many_minmatch311_scaffold0.initial.graph-report.tsv` |
| SPOA2kb/default-selfloop verified | `/home/erikg/impg/data/c4_spoa2k_default_selfloop_verify_20260606T065905Z/reports/c4.k311.seed.spoa2k.default-selfloop.Ygs.graph-report.tsv` |

Selected Ygs graph-report metrics:

| artifact | status | failures | warnings | segments | links | path_steps | self_loop_edges | adjacent_same_step_path_steps | self_loop_runs | link_jump_max | path_jump_max | whitespace_p99_bp | whitespace_max_bp | whitespace_bridges | occupancy_frac | whitespace_frac | duplicate_frac | local_repeat_contexts | singleton_nodes | singleton_bp | povu_sites | povu_leaf_sites |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 1:many scaffold0 seed | REVIEW | internal_tips>0,common_end_frac,link_jump_max_frac | duplicate_sequence_frac,local_repeat_contexts,direct_self_loop_edges,adjacent_same_step_path_repeats | 7151 | 9716 | 3129613 | 60 | 263737 | 29415 | 6669 | 6669 | 858 | 232746 | 30956 | 0.798984 | 0.201016 | 0.653195 | 621 nodes / 1043 occurrences | 940 | 5531 | 2745 | 2651 |
| 1:1 no-scaffold seed | REVIEW | internal_tips>0,common_end_frac,link_jump_max_frac | duplicate_sequence_frac,local_repeat_contexts | 7346 | 10094 | 2857231 | 0 | 0 | 0 | 2012 | 2012 | 230 | 37388 | 5793 | 0.791122 | 0.208878 | 0.660632 | 622 nodes / 1049 occurrences | 959 | 6082 | 2671 | 2575 |
| verified SPOA2kb/default-selfloop | REVIEW | internal_tips>0,common_end_frac | duplicate_sequence_frac,local_repeat_contexts | 8350 | 11289 | 3121673 | 0 | 0 | 0 | 1134 | 1134 | 141 | 25254 | 6740 | 0.791178 | 0.208822 | 0.724072 | 651 nodes / 1096 occurrences | 1031 | 6748 | 2829 | 2646 |
| 1:1 no-scaffold SPOA2kb/default-selfloop | REVIEW | internal_tips>0,common_end_frac | duplicate_sequence_frac,local_repeat_contexts | 8301 | 11147 | 3112253 | 0 | 0 | 0 | 1135 | 1135 | 159 | 25250 | 7198 | 0.781303 | 0.218697 | 0.724130 | 666 nodes / 1115 occurrences | 1040 | 7643 | 2865 | 2617 |

Interpretation:

- Against the 1:many scaffold0 seed, the 1:1 seed is better on self-loop and
  repeat-run metrics: direct self-loop edges drop from 60 to 0,
  adjacent same-step path steps drop from 263737 to 0, and self-loop repeat
  runs drop from 29415 to 0.
- Against the 1:many scaffold0 seed, the 1:1 seed also reduces the worst
  path/link jump from 6669 to 2012 and reduces max path whitespace from
  232746 bp to 37388 bp.
- The 1:1 seed still has `REVIEW` status and the same failure classes for
  internal tips, common end fraction, and link jump max fraction. It also has
  more segments, links, singleton nodes, duplicate sequence fraction, and
  segment whitespace fraction than the 1:many seed.
- Against the verified SPOA2kb/default-selfloop output, the 1:1 polished graph
  is not better overall. It is slightly smaller in segments, links, and path
  steps, but has worse occupancy/whitespace fractions, higher whitespace bridge
  count, higher sparse/singleton burden, and more local repeat contexts.

## Failure Logs

The tracked short-temp run had no failed steps; every `/usr/bin/time -v` log
reported exit status 0.

The earlier long-temp no-scaffold attempts failed before alignment, during
FastGA GDB preparation:

| attempt | failing log | observed failure |
| --- | --- | --- |
| `/home/erikg/impg/data/c4_k311_1to1_noscaffold_20260606T065737Z` | `logs/seed.graph.stderr.log` | `fastga alignment failed: Failed to prepare GDB: FAtoGDB failed with code None`; stderr included `*** buffer overflow detected ***: terminated`; exit status 1 |
| `/home/erikg/impg/data/c4_k311_1to1_noscaffold_retry_20260606T065905Z` | `logs/seed.graph.stderr.log` | `fastga alignment failed: Failed to prepare GDB: FAtoGDB failed with code None`; exit status 1 |

Both failures occurred immediately after the logs reported
`[graph::align] ... Running fastga alignment` and before a usable alignment or
seed GFA was produced.
