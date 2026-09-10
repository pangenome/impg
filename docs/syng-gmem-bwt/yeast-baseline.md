# Yeast names inventory and explicit-query baseline

`scripts/yeast_inference_baseline.py` is a Python 3.9+ **stdlib-only** stage 1/2
experiment harness. It inventories syng source names and runs a bounded batch
of coordinate queries against an existing full-panel syng index plus AGC.
It does not build an index, parse/open AGC sequence contents, run partition,
retrieve source spelling, implement a sample gMEM-BWT, or genotype anything.

All rows are **preliminary candidate intervals**, not proven spanning
haplotypes or validated homology. JSON metadata and candidate TSV eligibility
fields explicitly leave homology/spanning checks pending. Source spelling,
chromosome identity, ploidy, and assembly provenance also remain pending.
Nothing here certifies a nuclear chromosome, an archive sample count, panel
completeness, or passage of the full stage 1/2 gates.

## Inventory first

Use the actual build prefix and its names sidecar, not a strain-name guess:

```sh
python3 scripts/yeast_inference_baseline.py inventory \
  --names /path/to/panel.syng.names \
  --out-dir /path/to/new-inventory
```

Outputs are `sequences.tsv` and `inventory.json`: numeric source IDs, exact
complete names, positive lengths, conservative parsed PanSN identities,
all source records/base totals, parsed sample and sample/haplotype token
counts, and unparsed-name/legacy-format ambiguity diagnostics. No path is
classified as nuclear, mitochondrial, unplaced, haploid, or diploid by its name.
Source coverage here means the names-table inventory, not independent AGC
verification or complete biological assembly coverage.

Supported names tables have 3, 6, or 7 tab-separated columns, without a header.
Only ID/name/length are interpreted; optional native anchor columns are not
validated. IDs must fit unsigned 32-bit integers; lengths must be positive
unsigned 64-bit integers. Duplicate numeric IDs, duplicate exact names, empty
names/control characters, malformed rows, and invalid lengths are errors.
Blank lines are ignored. Names are never split on spaces or normalized.
PanSN parsing requires the **whole name** to have exactly three nonempty
`sample#haplotype#contig` components and no whitespace. Description-bearing
names remain exact, unparsed source names; even parsed tokens are not evidence
of chromosome identity or biological ploidy. Six-column tables are flagged for
missing first-syncmer offsets. Arbitrary unique IDs are inventoried, but the
pilot requires IDs `0..N-1` in row order, matching impg's loader behavior.

## Bounded explicit-query pilot

After independently identifying the intended source path, copy its complete
name from the inventory and quote it. The parent-run inventory confirmed the
exact name `S288C#0#chrI`; this example uses that spelling without claiming
independent chromosome or biological validation by this harness:

```sh
python3 scripts/yeast_inference_baseline.py pilot \
  --index-prefix /path/to/panel.syng \
  --agc /path/to/yeast235.agc \
  --impg /path/to/impg \
  --reference 'S288C#0#chrI' \
  --start 0 --window-size 10000 --max-windows 3 \
  --threads 4 --merge-distance 0 --timeout 300 \
  --out-dir /path/to/new-pilot
```

The reference must match exactly; there is no fuzzy strain/chromosome lookup.
All coordinates are zero-based, half-open. Consecutive seeds start at `start`
and advance by `window-size`, stopping at `max-windows` or source end. The last
window is clipped to the source end, never dropped or extended out of bounds.
`start` must lie within the source. Seeds are uniquely labelled BED4 rows
(`seed_000000`, etc.). Only this source is seeded; the query searches the
supplied full index, without coordinate slicing other source paths.

The harness requires `.1gbwt`, `.1khash`, and the names/spos/pstep/meta sidecars
of the current complete build. For prefix `panel.syng`, sidecars are
`panel.syng.names`, etc.; for prefix `panel`, they are `panel.syng.names`, etc.
It does not inspect native sidecar contents or prove they belong to the AGC.
Use the parent's frozen build manifest/checksums to establish that linkage.

After one `impg --version`, the harness invokes **one** batched command:

```text
impg query -a PREFIX -b seeds.bed --sequence-files AGC -o bed -t THREADS -d DISTANCE --consider-strandness
```

`--consider-strandness` prevents the default cross-strand merge policy from
obscuring orientations. Other query parameters use the recorded executable's
defaults. impg may still merge nearby/overlapping intervals according to its
query policy and merge distance. The harness preserves every **emitted** row,
including repeated identical intervals and occurrences on the same source
under different labels. It cannot reconstruct occurrences already merged by
impg, or infer alignment, repeat context, boundary anchors, or spanning status
from BED intervals.

**Seed absence among refined hits is expected, not a failure.** In
`src/syng_transitive.rs`, `query_transitive_ext_with_seed_filter` marks the exact
seed tuple visited before traversal and emits only newly visited hits. The
harness therefore preserves `seeds.bed` separately and never requires the seed
path among hits or fabricates an extra emitted row. A future inference catalog
must explicitly retain its seed candidate rather than accidentally excluding it
when importing refined query output. Summary metadata records this policy.

A parent-run smoke of `S288C#0#chrI:100000-110000` returned 240 intervals on 240
paths but only 210 sample/haplotype naming identities, with lengths 476–10024 bp
and median 9999 bp. These parent observations are not harness measurements:
path count is not haplotype count, and none of the 240 intervals is thereby
proven spanning. The observed BED shape is the six fields
`path`, `start`, `end`, `seed label`, `.`, `+` (with tabs between fields).

### Artifacts and failure behavior

- `seeds.bed`, `inventory.json`, `sequences.tsv`: exact inputs and inventory.
- `argv.json`, `command.sh`, `manifest.json`: resolved executable/input paths,
  argument array, shell-quoted command, working directory, version, parameters,
  per-command duration/return status and input identity metadata.
- `version.txt`, `version.stderr.log`, `query.stderr.log`: subprocess outputs.
- `candidates.bed`: raw query stdout, including partial output if a command fails.
- `candidates.tsv`: one occurrence ID per validated BED6 row, preserving its
  label, source ID/name, start/end, strand, length, score and raw line number.
- `summary.json`: each seed's row count, exact unique paths/count, strand counts,
  min/max/sum/mean candidate length, exact-seed-interval count and diagnostics.
  Every source has seeded/candidate coordinate-union bases and uncovered bases;
  unions ignore orientation and do not double-count repeated rows. Candidate
  counts/length sums **do** retain multiplicity. These are coordinate coverage
  statistics, not callable or homologous coverage.
- `failure.json`: when a subprocess or output validation fails, every seed has
  a failure diagnostic and unknown (`null`), not falsely zero, candidate count.
  Preflight errors are reported to stderr before creating a result directory.

Failed/timed-out commands, unknown names/labels, malformed BED6, bad strand or
score, empty/reversed/out-of-bounds intervals, and excess candidate rows return
exit code 2. No summary is accepted from a failed or partial query. A successful
empty BED is valid and yields `no_candidates` for every affected seed, with
null length extrema/mean. Seeds with only their own exact intervals get
`only_exact_seed_intervals`. Neither diagnostic proves biological absence;
the presence or absence of seed intervals is not independent homology evidence.

All subcommands refuse nonempty output directories (and symlink directories).
An exclusively created `.baseline-reserved` marker prevents concurrent harness
runs from sharing an initially empty directory; it remains after success or
failure. JSON is atomically replaced, and candidate TSV output is completed
before the summary is published. Handled publication failures discard accepted
result artifacts. Require `manifest.json` status `succeeded` as well as the
summary for pilot acceptance; abrupt termination can leave staging files or a
running manifest, which are not accepted results.
Use a fresh directory for reruns; never treat preserved failed raw BED as an
accepted result. `command.sh` is an audit/replay command, not the guarded pilot:
executing it directly redirects over the raw files and does not apply the
harness timeout. Prefer rerunning the Python CLI into a fresh directory.

### Resource and reproducibility limits

`max-windows`, threads, and subprocess timeout are explicit. Timeout applies
separately to version and query, not to all Python preprocessing/postprocessing.
`--max-candidate-rows` defaults to 1,000,000 and bounds accepted rows held in
memory; lower it for small pilots. Raw stdout/stderr are streamed to disk, not
kept in subprocess memory, but their disk size is **not capped**. The timeout
terminates the direct subprocess, not an arbitrary wrapper's process tree.
Names/seed tables and accepted candidate records are held in memory; there is
no OS-level memory limit. This is a bounded experiment, not a hostile-input
sandbox or a whole-panel scaling benchmark.

Names, seeds, harness and executable get SHA-256 identifiers. Large index/AGC
files get absolute paths, byte sizes and mtimes only, avoiding a second full
archive scan. The caller environment is inherited, not hermetically captured.
The frozen build/input checksum manifest must accompany published runs;
matching paths alone cannot establish byte-identical inputs. No real yeast
query or biological validation is claimed by the harness tests.

## Summarize an existing labelled BED

```sh
python3 scripts/yeast_inference_baseline.py summarize \
  --names /path/to/panel.syng.names \
  --seeds /path/to/pilot/seeds.bed \
  --bed /path/to/pilot/candidates.bed \
  --max-candidate-rows 100000 \
  --out-dir /path/to/new-summary
```

The seed file must be nonempty BED4 with unique labels and in-bounds intervals.
Candidate input is BED6 with exactly six tab-separated columns, score `.` or
an integer 0–1000, and strand `+`/`-`. Blank lines and `#` comments are ignored
(exact known source names beginning with `#` are still treated as records).
All candidate labels must occur in the seed file. External BED has
`command_status: not_verified_external_bed`: syntax validation cannot establish
that the producing command succeeded. Keep and inspect its original manifest.

## Real-data smoke and focused tests

The parent also ran the four-seed pilot on
`S288C#0#chrI:100000-140000` against the full yeast index. It returned 960
preliminary intervals: 240, 235, 245 and 240 rows per seed. The third seed
included eight reverse-strand rows; the fourth retained multiple intervals
on some paths (240 rows on 236 paths). The hardened output-publication version
was rerun successfully into
`~/yeast/baseline-pilot-chrI-4windows-hardened/`. These checks establish CLI
integration and interval parsing, not candidate completeness or spanning
support. No real dataset is required by the fixture CI workflow.


```sh
python3 -m unittest discover -s tests -p 'test_yeast_inference_baseline.py' -v
```

The 29 tests use tiny text fixtures and a fake executable only. They include
interleaved output reservations, injected TSV/JSON/manifest write failures,
and preparation-stage diagnostics, in addition to names formats and
ambiguity, duplicates/length validation, clipped tail windows, exact occurrence
and label retention, strand/length/coverage statistics, malformed BED and bounds,
no-hits/self-only diagnostics, valid refined output without the seed,
command/version/timeout failures, replay artifacts,
row limits, external summaries, and overwrite refusal. No Rust/core code, AGC
parser, datasets, installation, network access, or real impg commands are needed.
