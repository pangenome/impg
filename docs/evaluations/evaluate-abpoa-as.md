# Evaluate abPOA as C4 Crush POA Resolver Tier

Task: `evaluate-abpoa-as`

Date: 2026-06-04

## Recommendation

Add abPOA as an explicit experimental crush resolver method, but do not replace the
current SPOA/POASTA automatic local tier yet.

The hard gate passed: explicit abPOA preserved every C4 slice path sequence exactly
and rejected a synthetic corruption test. However, on the small C4 slice it expanded
topology substantially compared with both the input slice and the referenced cap64
diagnostic run. That is acceptable for a manually selected alternative, but not a
good default tier replacement without more work on graph shaping/polishing.

## Implementation

- Added `ResolutionMethod::Abpoa` plus CLI/pipeline parsing for `method=abpoa` and
  `--abpoa-bin` / `abpoa-bin=` in `src/resolution.rs:93`,
  `src/resolution.rs:7813`, `src/main.rs:2756`, and `src/main.rs:5501`.
- abPOA external IO uses the original semantic source path token as the FASTA
  header and expected GFA path name. It does not call the local synthetic
  uniquifier used by other POA paths; duplicate semantic names fail before tool
  invocation instead of being rewritten. See `src/resolution.rs:8094` and
  `src/resolution.rs:8134`.
- abPOA invocation is global alignment and emits GFA with `-r 3`; scoring is wired
  from the existing crush POA scoring tuple. See `src/resolution.rs:8994`.
- Path sequence preservation remains the only acceptance gate in the resolver:
  abPOA output is parsed, reordered by exact source path names, and validated
  against expected traversal sequences before integration. Quality/compression
  metrics are only logged diagnostics.

## Tests

- Unit path-preservation test with a fake abPOA binary verifies semantic path
  names are passed through exactly and no synthetic local IDs appear in external
  FASTA/GFA IO: `src/resolution.rs:10743`.
- Unit corruption test verifies changed abPOA output is hard-gated and leaves the
  original graph unchanged: `src/resolution.rs:10783`.
- CLI parser coverage verifies pipeline syntax accepts `method=abpoa` and
  `abpoa-bin=`: `src/main.rs:14601`.
- C4 slice smoke test uses a real abPOA binary when `IMPG_ABPOA_BIN` is set and
  asserts exact preservation for all output paths: `tests/test_crush_integration.rs:1081`.

## Evaluation Setup

- abPOA source: `https://github.com/yangao07/abPOA.git`
- abPOA commit: `b29b9f4e8641650597d6bf214313aef5910a4303`
- abPOA version: `1.5.6`
- Binary used locally: `/tmp/abPOA-evaluate-abpoa-as/bin/abpoa`
- Slice input: `tests/test_data/crush/c4_slice_1500_3000.gfa`
- Referenced baseline diagnostics:
  `/home/erikg/impg/data/c4_cap64_diag_20260604T154104Z`

The explicit abPOA artifact run used:

```bash
target/debug/impg crush \
  -g tests/test_data/crush/c4_slice_1500_3000.gfa \
  -o data/evaluate_abpoa_as_20260604T181500Z/c4_slice_abpoa.gfa \
  --method abpoa \
  --abpoa-bin /tmp/abPOA-evaluate-abpoa-as/bin/abpoa \
  --max-iterations 1 \
  --max-traversal-len 10k \
  --max-median-traversal-len 10k \
  --max-total-sequence 1m \
  --max-traversals 10k \
  -t 4 -v 1
```

## Results

Primary artifacts:

- `data/evaluate_abpoa_as_20260604T181500Z/c4_slice_abpoa.gfa`
- `data/evaluate_abpoa_as_20260604T181500Z/c4_slice_abpoa.graph-report.tsv`
- `data/evaluate_abpoa_as_20260604T181500Z/c4_slice_abpoa.gfalook-compressed.png`
- `data/evaluate_abpoa_as_20260604T181500Z/c4_slice_abpoa.time.txt`
- `data/evaluate_abpoa_as_20260604T181500Z/c4_slice_abpoa.path-preservation.tsv`
- `data/evaluate_abpoa_as_20260604T181500Z/abpoa_vs_c4_cap64_diag.summary.tsv`

Path preservation:

- `PASS`: 465 paths before, 465 paths after, 0 missing, 0 extra, 0 changed.
- Resolver log: 147 resolved, 0 bailed, 147 candidates seen across 1 round.
- Timed run: wall `1:07.77`, max RSS `1,360,772 KB`.

Selected graph-report comparison:

| Scenario | Segments | Links | Path Steps | Segment bp | POVU Sites | POVU Leaves | Duplicate Seq Frac | Runtime | Max RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| C4 slice input | 2,942 | 3,319 | 697,500 | 64,348 | 313 | 292 | 0.320530 | n/a | n/a |
| C4 slice abPOA | 21,529 | 22,150 | 8,764,488 | 58,031 | 546 | 453 | 0.933485 | 1:07.77 | 1,360,772 KB |
| cap64 diag spoa2k/poasta15k | 7,405 | 10,122 | 3,619,762 | 242,307 | 0 | 0 | 0.675219 | 4:10.66 | 4,712,032 KB |

The cap64 diagnostic is not the same input as the small slice run; it is included
because the task requested comparison against
`data/c4_cap64_diag_20260604T154104Z`. The useful signal is directional:
abPOA is exact and fast enough on the small slice, but it produces many more
segments/path steps and leaves more residual POVU sites than the cap64 diagnostic
row reports for the larger scheduled run.

## Decision

Keep the existing SPOA/POASTA automatic routing. Add abPOA as an explicit
resolver for manual experiments and follow-up work, not as the default C4 local
tier. The exact-path hard gate is satisfied, but the graph-shape result is not
yet strong enough to justify replacing POASTA in auto mode.
