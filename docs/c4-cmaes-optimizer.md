# C4 CMA-ES Optimizer

`scripts/c4-crush-cmaes.py` is an external CMA-ES study runner for the C4
syng plus crush graph-building problem. It does not add any quality gate inside
`impg`: every candidate is run as an ordinary `impg crush` or `impg query`
command, `impg graph-report` is run afterwards, and this wrapper alone marks
trials invalid or assigns objective values.

The primary optimizer is the installed Python `cma` package. If importing
`cma` fails, the script falls back to a deterministic seeded random sampler so a
study can still produce auditable trial artifacts, but the expected production
path is ask/tell CMA-ES.

## Modes

`--mode crush-only` optimizes `impg crush` on an existing blunt GFA:

```bash
scripts/c4-crush-cmaes.py \
  --mode crush-only \
  --input-gfa input.gfa \
  --method-family sweepga \
  --study-dir data/c4_cmaes/sweepga \
  --max-trials 80 --population-size 8 --jobs 2 --threads 16
```

For every crush-only trial the wrapper validates that the output GFA has the
same path names and exact spelled path sequences as the input GFA. Missing
paths, extra paths, or sequence spelling differences make the trial invalid.

`--mode full-query` optimizes full `impg query` GFA generation over syng engine
strings:

```bash
scripts/c4-crush-cmaes.py \
  --mode full-query \
  --index /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  --target-range 'GRCh38#0#chr6:31982056-32035418' \
  --expected-path-count 465 \
  --syng-engine syng-local \
  --method-family coverage-multi-bubble \
  --study-dir data/c4_cmaes/full_coverage \
  --max-trials 80 --population-size 8 --jobs 2 --threads 16
```

Full-query trials build an output string of the form:

```text
gfa:syng-local:blunt,k=<K>,s=<S>,seed=<SEED>:mask,top=<F>,max-occ=<N>,min-run=<R>,sequence-k=<K2>:crush,...:sort,pipeline=Ygs
```

or:

```text
gfa:syng:blunt:mask,top=<F>,max-occ=<N>,min-run=<R>,sequence-k=<K2>:crush,...:sort,pipeline=Ygs
```

The hard constraints for full-query are:

- the primary `impg query` command exits successfully;
- exactly one output GFA is found for scoring;
- `impg graph-report --format tsv` exits successfully;
- the output path count and report path count both match
  `--expected-path-count`.

Quality metrics never gate a candidate. They only contribute to the external
objective after hard constraints pass.

## Top-Level Method Family

The method family is a study-level argument, not a hidden continuous
categorical dimension. Each study uses exactly one `--method-family`:

- `auto`
- `sweepga`
- `chain-povu`
- `iterative-multi-level`
- `coverage-multi-bubble`

Run separate studies to compare method families. This keeps CMA-ES operating on
continuous/numeric neighborhoods rather than trying to learn discontinuous
method switches.

## CMA Vector Decode

CMA-ES operates over normalized continuous vectors bounded to `[0, 1]`. Each
dimension has an explicit decode function:

- `linear integer [a, b]`: `round(a + x * (b - a))`
- `linear odd integer [a, b]`: linear integer decode, adjusted to the nearest
  odd value inside the range
- `log integer [a, b]`: `round(exp(log(a) + x * (log(b) - log(a))))`
- `linear float [a, b]`: `a + x * (b - a)`
- `log float [a, b]`: `exp(log(a) + x * (log(b) - log(a)))`

The decoded common crush parameters are normalized before command generation so
`max_median_traversal_len <= max_traversal_len` and
`polish_max_median_traversal_len <= polish_max_traversal_len`.

## Parameter Spaces

All method families share these crush dimensions:

| Dimension | Decode | impg parameter |
| --- | --- | --- |
| `max_rounds` | linear integer `[1, 4]` | `max-rounds` |
| `min_traversal_len` | log integer `[1, 20000]` | `min-traversal-len` |
| `max_traversal_len` | log integer `[1000, 120000]` | `max-traversal-len` |
| `max_median_traversal_len` | log integer `[200, 30000]` | `max-median-traversal-len` |
| `max_traversals` | log integer `[64, 25000]` | `max-traversals` |
| `polish_rounds` | linear integer `[1, 4]` | `polish-rounds` |
| `polish_max_traversal_len` | log integer `[500, 50000]` | `polish-max-traversal-len` |
| `polish_max_median_traversal_len` | log integer `[100, 10000]` | `polish-max-median-traversal-len` |
| `replacement_flank_bp` | log integer `[1, 1000]` | `replacement-flank-bp` |

Full-query studies additionally optimize syng mask dimensions:

| Dimension | Decode | Engine parameter |
| --- | --- | --- |
| `mask_top_fraction` | log float `[0.00005, 0.003]` | `mask,top` |
| `mask_max_occ` | log integer `[100, 2500]` | `mask,max-occ` |
| `mask_min_run` | linear integer `[2, 8]` | `mask,min-run` |
| `mask_sequence_k` | linear odd integer `[63, 311]` | `mask,sequence-k` |

Full-query `--syng-engine syng-local` also optimizes local rebuild dimensions:

| Dimension | Decode | Engine parameter |
| --- | --- | --- |
| `syng_k` | linear odd integer `[63, 311]` | `syng-local,k` |
| `syng_s` | linear integer `[8, 31]`, clamped below `syng_k / 2` | `syng-local,s` |

### `auto`

| Dimension | Decode | impg parameter |
| --- | --- | --- |
| `auto_poasta_max_traversal_len` | log integer `[1000, 60000]` | `auto-poasta-max-len` |

The auto family intentionally optimizes the stable POASTA routing cutoff while
leaving the smaller internal auto routing defaults to `impg`.

### `sweepga`

| Dimension | Decode | impg parameter |
| --- | --- | --- |
| `k_nearest` | linear integer `[1, 8]` | `k-nearest` |
| `k_farthest` | linear integer `[0, 4]` | `k-farthest` |
| `pair_trees` | linear integer `[1, 4]` | `pair-trees` |
| `random_fraction` | log float `[0.0001, 0.05]` | `random-fraction` |
| `mash_k` | linear integer `[11, 21]` | `mash-k` |
| `seqwish_k` | linear odd integer `[51, 501]` | `seqwish-k` |
| `min_match_length` | log integer `[31, 501]` | `min-match-length` |
| `replacement_min_map_length` | log integer `[50, 2000]` | `replacement-min-map-length` |
| `replacement_min_identity` | linear float `[0.0, 0.99]` | `replacement-min-identity` |
| `max_pair_alignments` | log integer `[500, 100000]` | `max-pair-alignments` |
| `max_replacement_paf_bytes` | log integer `[1000000, 512000000]` | `max-replacement-paf-bytes` |
| `sweepga_kmer_frequency` | linear integer `[0, 100]` | `kmer-frequency` in engine strings, `sweepga-kmer-frequency` for `impg crush` |

### `chain-povu`

| Dimension | Decode | impg parameter |
| --- | --- | --- |
| `chain_target_bp` | log integer `[2000, 100000]` | `chain-target-bp` |

### `iterative-multi-level`

| Dimension | Decode | impg parameter |
| --- | --- | --- |
| `window_target_bp` | log integer `[2000, 120000]` | `window-target-bp` |
| `max_window_sites` | linear integer `[2, 20]` | `max-window-sites` |
| `candidate_limit` | log integer `[16, 768]` | `candidate-limit` |
| `min_objective_delta` | linear integer `[0, 10]` | `min-objective-delta` |

The study-level `--multi-level-window-mode` selects the fixed window source set
(`sibling`, `sliding`, `outward`, or `combined`). The generated crush stage uses
`objective=size` and does not enable repeat-aware boundaries.

### `coverage-multi-bubble`

`coverage-multi-bubble` uses the same four dimensions as
`iterative-multi-level`, but the generated crush stage uses
`objective=coverage` and `repeat-aware-boundaries=true`.

## Objective

The objective is minimized. Invalid trials receive:

```text
1e12 + trial_id
```

Valid trials are scored from raw `impg graph-report --format tsv` fields:

```text
score =
  total_segment_bp / bp_scale
+ 2 * singleton_bp / total_segment_bp
+ 10 * path_white_space_bp_p99 / total_segment_bp
+ 5 * path_white_space_bp_max / total_segment_bp
+ path_white_space_bridges_ge_threshold / paths
+ duplicate_sequence_frac
+ local_repeat_context_nodes / segments
- node_coverage_bp_weighted_mean / paths
```

`bp_scale` defaults to the input GFA segment bp in crush-only mode. In
full-query mode it defaults to `1,000,000` unless overridden by
`--objective-bp-scale`.

The coverage term is a reward: higher bp-weighted node coverage reduces the
score after normalization by path count. All raw graph-report fields are still
recorded in each trial log; the objective terms are an additional derived view.

## Audit and Resume Artifacts

Each trial gets a dedicated directory:

```text
trial-000123/
  command.sh
  decoded.json
  stdout.log
  stderr.log
  time.txt
  output.gfa
  graph-report.tsv
  graph-report.stdout.log
  graph-report.stderr.log
  graph-report.time.txt
  path-validation.json
  trial.json
```

The study root also contains:

```text
study.json
trials.jsonl
```

`trials.jsonl` is the resumable audit log. On resume, the script reconstructs
the CMA-ES state by replaying every complete generation in trial-log order. A
crash in the middle of a generation leaves those partial trial directories and
JSONL rows auditable, but only complete generations are replayed into CMA-ES.
The next run continues with fresh candidates and monotonically increasing trial
IDs.

The `command.sh` file is the exact impg command sequence for the trial. Timing
and RSS are captured with `/usr/bin/time -v` when available; stdout/stderr are
kept separate for the primary impg command and graph-report command.

## Validation Commands

Implementation smoke checks are intentionally tiny:

```bash
python -m py_compile scripts/c4-crush-cmaes.py
scripts/c4-crush-cmaes.py --help
scripts/c4-crush-cmaes.py \
  --mode crush-only \
  --input-gfa tests/test_data/crush/small_insertion.gfa \
  --method-family sweepga \
  --max-trials 1 \
  --population-size 2 \
  --study-dir /tmp/c4-cmaes-dry-run \
  --dry-run
```

Do not use implementation validation to launch long C4 sweeps. Real studies
should be scheduled explicitly with bounded `--max-trials`, `--jobs`,
`--threads`, and `--timeout`.
