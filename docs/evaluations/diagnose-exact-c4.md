# Diagnose exact C4 seed-induction stall

Date: 2026-06-09

Task: `diagnose-exact-c4`

All generated FASTA/PAF/GFA/render/log artifacts were kept outside git under:

```text
/home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z
```

## Summary

The C4 seed path is not bottlenecked by SweepGA/FastGA or seqwish. The release
reproduction shows:

- SweepGA/FastGA invocation: 60.803s.
- SweepGA filter: 3.423s.
- seqwish induction after handoff: 153.753s.
- Hidden IMPG wrapper scan between PAF filtering and seqwish handoff: 78.282s.

The concrete stall boundary is after PAF filtering and before seqwish handoff in
`src/syng_graph.rs`, where IMPG unconditionally replayed every filtered PAF
CIGAR against the FASTA to compute `max_exact_match_run`. For this C4 input that
scan covered 431,600 filtered PAF records and 25,953.2 Mb of aligned block bp.
In release it took 78.282s. In the prior timeout, the same work was running in
`target/debug/impg` with no graph-stage progress logs, so it appeared as a seed
induction stall after the filter.

I implemented the small fix: the full exact-run scan now runs only when
`adaptive_min_match_len` is enabled. The normal C4 path only counts filtered PAF
records before handing the PAF to seqwish.

## Commands

### Internal C4 reproduction

This reproduces the `validate-localized-polishing` C4 tuned seed command with
the same syng index, AGC input, BED interval, syng-local localized output, C4
limits, and tuning parameters. The run used `target/release/impg` so the
boundary timings are comparable to the known-fast standalone release command.

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

Result:

- Exit code: 0.
- Wall time: 11:00.11.
- Max RSS: 49,943,524 KB.
- Query log complete: `2026-06-09T09:24:06Z`.
- Graph: `$out/graphs_tuned/C4_GRCh38_53kb.gfa`.

### Standalone comparison using identical internal artifacts

This compares the internal seed tail against standalone `impg graph` using the
exact FASTA and raw PAF emitted by the internal run:

- `$out/debug/internal_tail/graph_build_0000/combined.fa`
- `$out/debug/internal_tail/graph_build_0000/raw.paf`

```bash
source ./env.sh >/dev/null
out=/home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z
/usr/bin/time -v -o "$out/logs/standalone.release.same_artifacts.time.txt" \
  /usr/bin/timeout 900s target/release/impg graph \
  --sequence-files "$out/debug/internal_tail/graph_build_0000/combined.fa" \
  --paf-file "$out/debug/internal_tail/graph_build_0000/raw.paf" \
  --gfa-engine seqwish \
  --fastga \
  --num-mappings 1:many \
  --scaffold-filter 1:many \
  --scaffold-jump 0 \
  --min-match-len 311 \
  --temp-dir "$out/tmp/standalone_same" \
  --debug-dir "$out/debug/standalone_same" \
  -g "$out/graphs_standalone/standalone.same-artifacts.gfa" \
  -t 16 -v 1 \
  > "$out/logs/standalone.release.same_artifacts.stdout.log" \
  2> "$out/logs/standalone.release.same_artifacts.stderr.log"
printf '%s\n' "$?" > "$out/logs/standalone.release.same_artifacts.exitcode"
```

Result:

- Exit code: 0.
- Wall time: 2:47.77.
- Max RSS: 4,526,708 KB.
- Graph: `$out/graphs_standalone/standalone.same-artifacts.gfa`.

This confirms that the same PAF filters and induces with seqwish normally when
run through the standalone graph path.

## Internal timing table

Timestamps are UTC from
`$out/logs/internal.release.tuned.stderr.log`; wall and RSS are from
`$out/logs/internal.release.tuned.time.txt`.

| Boundary | Timestamp(s) | Wall time | Evidence |
|---|---:|---:|---|
| Syng index load | 09:13:20 to 09:16:36 | ~196s | Read 272,213,874 vertices and 683,539,246 edges |
| AGC index setup | 09:16:36 to 09:16:37 | ~1s | `Building AGC index for 1 file` |
| Syng query interval collection | 09:16:37 to 09:17:22 | ~45s | `Collected 466 merged intervals` |
| Local sequence extraction | 09:17:22 | 0.042s | 466 sequences, 53,165,920 bp |
| SweepGA/FastGA invocation | 09:17:22 to 09:18:23 | 60.803s | PAF path returned by FastGA wrapper |
| Raw PAF load/sort/dedup | 09:18:23 to 09:18:24 | 0.950s | 348,857,330 bytes, 541,854 unique records |
| Seed FASTA write | 09:18:24 | 0.044s | 466 sequences, 53,165,920 bp |
| Raw PAF temp write | 09:18:24 | 0.276s | 541,854 records, 348,857,330 bytes |
| SweepGA PAF filter | 09:18:25 to 09:18:29 | 3.423s | 541,388 to 431,600 mappings, 79.7% kept |
| IMPG exact-run PAF scan | 09:18:29 to 09:19:48 | 78.282s | 431,600 parsed records, `max_exact_match_run=145272` |
| Seqwish handoff | 09:19:48 | immediate after scan | 431,600 records, `min_match_len=311`, 16 threads |
| Seqwish sequence index | 09:19:48 to 09:19:49 | 0.095s | Indexed 466 sequences |
| Seqwish alignment indexing | 09:19:49 to 09:20:05 | 16.721s cumulative | Alignments indexed |
| Seqwish transclosure | 09:20:05 to 09:22:20 | 151.782s cumulative, ~135.061s stage | Transitive closure complete |
| Seqwish compaction | 09:22:20 | 0.036s stage | 3,850 nodes |
| Seqwish link derivation | 09:22:20 | 0.251s stage | 5,219 links |
| GFA generation/sort | 09:22:20 to 09:22:22 | 1.494s stage | Graph building complete at 153.563s |
| Seqwish induction return | 09:22:22 | 153.753s | 10,303,414 GFA bytes |
| Local seed completion | 09:22:23 | 300.730s from local seed start | 466 paths, path validation pass |
| Localized handoff/dirty-site scan | 09:22:23 to 09:23:22 | ~59s | First localized polish iteration log |
| Localized resolver chunk | 09:23:22 to 09:24:01 | 67.907s total polish | 1 chunk attempted, 1 applied, path validation pass |
| Rendering/report | 09:24:01 to 09:24:06 | ~5s | PNG and TSV written outside git |
| Whole process | time(1) | 11:00.11 | Max RSS 49,943,524 KB |

## Standalone timing table

Standalone used the internal run's exact `combined.fa` and `raw.paf`.

| Boundary | Timestamp(s) | Wall time | Evidence |
|---|---:|---:|---|
| FASTA count | 09:26:24 | 0.060s | 466 sequences in 466 genomes |
| FASTA combine | 09:26:24 | 0.126s | Average sequence length 114,089 bp |
| Provided PAF handoff | 09:26:25 | 0.242s | Used internal raw PAF |
| SweepGA PAF filter | 09:26:26 to 09:26:30 | 3.568s | 541,388 to 431,600 mappings |
| Seqwish sequence index | 09:26:31 | 0.102s | Indexed 466 sequences |
| Seqwish alignment indexing | 09:26:32 to 09:26:49 | 17.185s cumulative | Alignments indexed |
| Seqwish transclosure | 09:26:49 to 09:29:03 | 151.864s cumulative, ~134.679s stage | Transitive closure complete |
| Seqwish compaction | 09:29:03 | 0.030s stage | 3,850 nodes |
| Seqwish link derivation | 09:29:03 | 0.241s stage | 5,219 links |
| GFA generation/sort | 09:29:03 to 09:29:05 | 1.500s stage | Graph building complete at 153.635s |
| gfaffix/load/normalization | 09:29:05 to 09:29:07+ | included in full wall | 95 shared prefixes collapsed; self-loop normalization ran |
| Whole process | time(1) | 2:47.77 | Max RSS 4,526,708 KB |

The standalone and internal seed tail agree on the important boundaries:
SweepGA filtering is ~3.5s and seqwish graph induction is ~153s. The internal
path's extra pre-handoff exact-run scan is not present in the standalone graph
path.

## Root cause

The precise cause is IMPG wrapper logic, not SweepGA/FastGA and not seqwish.

Before handing the filtered PAF to seqwish, `build_gfa_from_paf_and_sequences`
called `filtered_paf_exact_match_run_stats` unconditionally. That function
parsed every filtered PAF row, looked up both path sequences, replayed `cg:Z`
CIGAR operations base by base, compared query and target bases, and computed
the longest exact run. This is only needed to adaptively lower
`min_match_len`, but the C4 syng-local command was not using adaptive min-match.

The prior timeout was made misleading by two wrapper/config choices:

1. Validation used `target/debug/impg`, while the known-fast standalone command
   used `target/release/impg`. The full CIGAR/base replay is especially costly
   in debug.
2. Query-side engine setup forced `pipeline.show_progress=false`, so graph
   stage logs were suppressed unless called through another path. The tail also
   only honored `IMPG_CRUSH_DEBUG_DIR`, not the configured `--debug-dir`.

The internal query also carries the full syng index in memory while standalone
`impg graph --paf-file` does not. That explains the RSS gap
(49,943,524 KB internal vs 4,526,708 KB standalone), but it is not the seed
stall boundary.

## Fix implemented

Code changes:

- `src/syng_graph.rs`
  - Added timing/progress logs around FASTA writing, raw PAF writing, PAF
    filtering, filtered PAF stats, seqwish handoff, and seqwish return.
  - Made graph-build debug artifacts honor `--debug-dir` as well as
    `IMPG_CRUSH_DEBUG_DIR`.
  - Changed filtered PAF exact-run scanning to run only when
    `adaptive_min_match_len=true`; otherwise it only performs a cheap nonempty
    record count before seqwish.
- `src/local_seed.rs`
  - Added separate timing logs for SweepGA/FastGA invocation and raw PAF
    load/sort/dedup.
- `src/main.rs`
  - Plumbed `-v`/verbose into query-side graph engine `show_progress`.
  - Added coverage in the existing syng-local CLI parser test.

Targeted tests added/updated:

- `syng_graph::tests::exact_run_scan_is_adaptive_only`
- `syng_graph::tests::non_adaptive_filtered_stats_counts_records_without_replaying_cigars`
- `tests::test_gfa_engine_accepts_syng_local_localized_polish_stage`

## Validation results

Completed during this task:

```text
cargo fmt
cargo test exact_run_scan_is_adaptive_only -- --nocapture
cargo test non_adaptive_filtered_stats_counts_records_without_replaying_cigars -- --nocapture
cargo test test_gfa_engine_accepts_syng_local_localized_polish_stage -- --nocapture
cargo build --release --bin impg --bin gfaffix
cargo build
cargo test
cargo install --path .
git diff --cached --check
staged artifact/blob scan
```

The first attempted focused test command used multiple Cargo test filters and
failed with Cargo's usage error; the three tests above were rerun individually
and passed.

`cargo build`, `cargo test`, and `cargo build --release --bin impg --bin
gfaffix` passed with pre-existing warnings only. The staged diff check passed,
and the staged artifact/blob scan found no generated GFA/GFA.zst/PNG/PDF/log
artifacts and no staged blob over 1 MiB.

## Artifact policy

Generated graph/render/log artifacts remain outside the repository under
`/home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z`.

Largest diagnostic files outside git:

| Artifact | Size |
|---|---:|
| `debug/internal_tail/graph_build_0000/raw.paf` | 348,857,330 bytes |
| `debug/internal_tail/graph_build_0000/filtered.paf` | 284,155,587 bytes |
| `debug/internal_tail/graph_build_0000/final.paf` | 284,155,587 bytes |
| `debug/internal_tail/graph_build_0000/combined.fa` | 53,186,044 bytes |
| `debug/internal_tail/graph_build_0000/seqwish.gfa` | 10,303,414 bytes |
| `graphs_tuned/C4_GRCh38_53kb.gfa` | 10,304,399 bytes |
| `graphs_standalone/standalone.same-artifacts.gfa` | 10,186,451 bytes |

None of these are staged or intended for commit.
