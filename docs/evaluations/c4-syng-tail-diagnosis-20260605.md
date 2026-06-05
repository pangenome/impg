# C4 syng 3-prime tail diagnosis

Task: `diagnose-and-eliminate`

Date: 2026-06-05

Artifact directory: `data/c4_syng_tail_diagnosis_20260605T145815Z`

## Reproducer

Current-code baseline reproduced the reported tail with:

```bash
target/debug/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  -d 50k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o fasta \
  --reverse-complement \
  -O data/c4_syng_tail_diagnosis_20260605T145815Z/d50k \
  -t 32 \
  -v 1
```

The baseline FASTA matched the prior observation: 465 records, median 226,371 bp,
max 298,173 bp, 15 records >265 kb, and 4 records >285 kb.

All matrix commands are recorded in:

- `data/c4_syng_tail_diagnosis_20260605T145815Z/commands.sh`
- `data/c4_syng_tail_diagnosis_20260605T145815Z/fraction_sweep/commands.sh`

The driver used for reproducible commands and metrics is
`scripts/c4_syng_tail_diagnosis.py`.

## Length Metrics

| setting | records | median | p95 | p98 | max | >265k | >285k | note |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `-d 50k` | 465 | 226,371 | 259,184 | 265,429 | 298,173 | 15 | 4 | reproduced tail |
| `-d 10k` | 465 | 226,371 | 259,184 | 265,429 | 298,173 | 15 | 4 | unchanged |
| `-d 1k` | 465 | 226,371 | 259,184 | 265,429 | 298,173 | 15 | 4 | unchanged |
| `-d 0` | 474 | 226,370 | 259,184 | 265,429 | 298,173 | 15 | 4 | tail remains |
| old `--no-merge` | 474 | 226,370 | 259,184 | 265,429 | 298,173 | 15 | 4 | old code still merged overlaps |
| `-d 50k --syng-extend-budget 0` | 465 | 226,373 | 259,184 | 265,429 | 298,173 | 15 | 4 | BiWFA extension budget is not the lever |
| `-d 50k --syng-min-chain-fraction 0.02` | 6 | 232,736 | 232,736 | 232,736 | 232,736 | 0 | 0 | drops expected haplotypes |
| `-d 50k --syng-min-chain-fraction 0.05` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | drops all |
| `-d 50k --syng-min-chain-fraction 0.10` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | drops all |
| `-d 50k --syng-min-chain-fraction 0.20` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | drops all |
| `-d 50k --syng-min-chain-fraction 0.50` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | drops all |
| `-d 50k --syng-min-chain-fraction 0.80` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | drops all |
| `-d 50k --syng-raw` | 53,530 | 312 | 13,358 | 23,698 | 232,940 | 0 | 0 | raw syncmer hits, not haplotype FASTA |
| fixed `--no-merge` | 4,586 | 44,738 | 85,610 | 102,560 | 165,158 | 0 | 0 | sidecar isolated as fragments |

Per-record length tables are in `*.lengths.tsv`; run stderr profiles are in
`*.stderr.log`.

## Stage Diagnosis

The extra 3-prime tail is introduced by default syng boundary-realignment
fragment merging, specifically `merge_nearby_on_same_path` in
`src/syng_transitive.rs`.

Evidence:

- Raw syncmer FASTA with `--syng-raw` has no high tail: max 232,940 bp, zero
  records >265 kb. This bypasses anchor chaining, BiWFA end projection, and
  same-path syng-transitive merging for FASTA output.
- Baseline boundary realignment emits 53,530 raw hits, SweepGA chains them to
  10,293 scaffold chains, filters to 4,586 chains, refines 4,586 intervals, then
  same-path merging collapses them to 466 intervals in the first hop. The final
  merge sees 465 intervals and does not merge further.
- `-d 10k`, `-d 1k`, `-d 0`, and old `--no-merge` all keep the same high tail.
  Therefore the problem is not the positive 50 kb merge distance.
- `--syng-extend-budget 0` keeps the same high tail, and the profile shows the
  BiWFA align time is small relative to interval fetch/refinement. Therefore
  BiWFA end extension is not the root lever.
- This repro does not use user transitive search (`-x` is absent). The syng path
  sets `syng_max_depth = 1`, so transitive hop collection is not the cause.
- `syng_intervals_to_merged_query_intervals` is not on the syng FASTA path. The
  FASTA branch writes the intervals returned by
  `query_transitive_ext_with_seed_filter` directly.

The exact collapse profile for the baseline is:

```text
PROF chain_anchors ... chains=10293 after_filter=4586 ... chain_gap=282738 ...
PROF refine ... chains=4586 ... out=4586 ...
PROF merge_same_path input=4586 output=466 merged_groups=466 gap=50000 ...
PROF merge_same_path input=465 output=465 merged_groups=0 gap=50000 ...
```

With `-d 0` and old `--no-merge`, overlapping projected fragments still collapse
into the same high-tail records. That is why changing merge distance alone did
not help.

## Fix

I fixed the concrete CLI semantics bug in `merge_nearby_on_same_path`: a
negative merge distance now returns sorted intervals unchanged. Before this
change, `--no-merge` still merged overlapping syng fragments.

The fixed diagnostic command is:

```bash
target/debug/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --no-merge \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -o fasta \
  --reverse-complement \
  -O data/c4_syng_tail_diagnosis_20260605T145815Z/no_merge_fixed \
  -t 32 \
  -v 1
```

Result: 4,586 records, median 44,738 bp, max 165,158 bp, zero records >265 kb
and zero >285 kb. This does not produce one merged FASTA per haplotype, but it
does preserve the underlying C4 fragments and clearly isolates the sidecar
instead of gluing it into the 15 high-tail haplotypes.

No tested existing setting produced a 465-record merged haplotype FASTA while
also removing the tail. `--syng-min-chain-fraction` is not viable for this
region: even 0.02 keeps only 6 records because valid C4 reconstruction depends
on many local chains whose query span is below 2% of the full query.

## Validation

Commands run:

```bash
bash -lc 'source ./env.sh >/tmp/impg-env.log && \
  CMAKE_C_COMPILER="$CC" CMAKE_CXX_COMPILER="$CXX" \
  CFLAGS="-I/gnu/store/2r5ryq2ibvy44jkz9diar0fvf5cm7q4c-profile/include" \
  LDFLAGS="-L/gnu/store/2r5ryq2ibvy44jkz9diar0fvf5cm7q4c-profile/lib" \
  cargo test syng_transitive::tests::merge_nearby_on_same_path_no_merge_keeps_touching_and_overlapping_intervals'

bash -lc 'source ./env.sh >/tmp/impg-env.log && \
  CMAKE_C_COMPILER="$CC" CMAKE_CXX_COMPILER="$CXX" \
  CFLAGS="-I/gnu/store/2r5ryq2ibvy44jkz9diar0fvf5cm7q4c-profile/include" \
  LDFLAGS="-L/gnu/store/2r5ryq2ibvy44jkz9diar0fvf5cm7q4c-profile/lib" \
  cargo build'

bash -lc 'source ./env.sh >/tmp/impg-env.log && \
  CMAKE_C_COMPILER="$CC" CMAKE_CXX_COMPILER="$CXX" \
  CFLAGS="-I/gnu/store/2r5ryq2ibvy44jkz9diar0fvf5cm7q4c-profile/include" \
  LDFLAGS="-L/gnu/store/2r5ryq2ibvy44jkz9diar0fvf5cm7q4c-profile/lib" \
  cargo test'

bash -lc 'mkdir -p target/release/build/wfmash-rs-61a2de665a016fe9/out && \
  cp /home/erikg/impg/target/release/build/wfmash-rs-61a2de665a016fe9/out/wfmash \
     target/release/build/wfmash-rs-61a2de665a016fe9/out/wfmash && \
  chmod +x target/release/build/wfmash-rs-61a2de665a016fe9/out/wfmash && \
  source ./env.sh >/tmp/impg-env.log && \
  CMAKE_C_COMPILER="$CC" CMAKE_CXX_COMPILER="$CXX" \
  CFLAGS="-I/gnu/store/2r5ryq2ibvy44jkz9diar0fvf5cm7q4c-profile/include" \
  LDFLAGS="-L/gnu/store/2r5ryq2ibvy44jkz9diar0fvf5cm7q4c-profile/lib" \
  cargo install --path .'
```

Results:

- Targeted regression test passed.
- `cargo build` passed.
- Full `cargo test` passed: 390 library tests, 83 main tests, syng/query
  integration tests, and the rest of the non-ignored integration tests.
- `cargo install --path .` passed and replaced the installed `impg` binary from
  this worktree.
- Build/test commands passed after initializing `vendor/syng` and `vendor/gfaffix`
  submodules and seeding the generated `wfmash-rs` cache from the existing local
  target. The default build environment initially failed because htslib headers
  and submodule files were not visible until `env.sh` and submodule
  initialization were used.
