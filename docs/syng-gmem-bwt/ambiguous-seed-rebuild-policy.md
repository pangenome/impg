# Ambiguous windows require a full syng rebuild

Construction previously substituted unknown bases with numeric A for C safety,
then inserted those windows as exact seeds. In the diagnosed
`S288C#0#chrXII:450000-460000` window, 8,389 source Ns accompanied 8,514
stored steps; 8,382 contained seed
windows overlapped ambiguity. This created a high-copy fake exact seed, not
biological repeat evidence. The neighboring N-free repeat windows are not masked.

Fresh construction now accepts a selected full `k`-base window only when **every
original byte** is ASCII A/C/G/T (either case) or numeric 0–3. N/n, other IUPAC
symbols, numeric 4+, and all other bytes invalidate the window. Temporary numeric
zero substitution remains solely to protect C's four-element lookup tables.
A constant-space streaming check uses the same alphabet as query matching.
Nothing is deleted or concatenated: first offsets, absolute positions, distances
across gaps, names, and lengths remain source-based. Sequential, dictionary
prepass/replay, FASTA/AGC streaming, and regional builds share this rule. The C
sampler, reverse-complement selection policy, defaults, and real repeats remain
unchanged.

## Yeast rebuild and hard-window replay

A fresh k63/s8/seed7 parallel-dictionary build completed in 5m08.68s, using
9,826,368 KiB peak RSS. The index sidecars total 444.08 MiB. All 9,901 source
IDs, full names and lengths match the old inventory (3,336,976,856 bp); no
source sequence was deleted. Stored steps decreased from 131,823,304 to
129,915,680 after excluding invalid windows.

| Source query | Stored steps: old → new | Ambiguity-overlapping steps: old → new | New raw lookup + chaining |
|---|---:|---:|---:|
| `S288C#0#chrXII:450000-460000` | 8514 → 70 | 8382 → 0 | 0.123 s |
| `S288C#0#chrXII:460000-470000` | 9058 → 35 | 8961 → 0 | 0.015 s |

The invalid-window check here covers fully contained stored seed windows.
Index loading is excluded from these single-query timings. In the old full-panel
run, processing the first interval took roughly 684s between partition progress
records; that is a different timing boundary, not a calibrated speedup ratio.
The query-count × panel-occurrence product fell from roughly 31/33 billion to
18,768/4,330, respectively. Panel occurrence counts include both orientations.

N-free controls at `S288C#0#chrXII:430000-440000`,
`S288C#0#chrXII:440000-450000`, and `S288C#0#chrXII:470000-480000` retained
identical stored-step counts, unique-node counts and occurrence totals. Node IDs
can change after a fresh dictionary build; old sample/index products are not
interchangeable with this index.

The hard windows still contain over 80% unknown source bases. With the unchanged
50% span filter they have no passing spanning chains; partition seed retention
can account for their source intervals as singletons. This does not invent
homology or recover missing assembly sequence. Full-panel completion and
biological chunk validation are separate gates.

Artifacts remain outside Git:

- `~/yeast/syng-k63-s8-seed7-acgt-only/`: complete index, metadata v2,
  original-input checksum, frozen build binary and executable checksum,
  build-source patch, command, resource log and sidecar checksums.
- `target/experiments/syng-ambiguity-fix/rebuilt-composition.tsv` and
  `rebuilt-hard-query.tsv`: the isolated replay and profiling logs alongside.
- `parent-native-validation.log` in the same directory: the parent rerun of
  all 413 library and 39 integration tests, all passing. Independent source
  review accepted the bounded fix with the compatibility limits below.

The original index at `~/yeast/syng-k63-s8-seed7/` remains unchanged for diagnosis.
The rebuilt directory's `impg-build-snapshot` is the matching tested binary;
an older installed impg will reject the new metadata version.

## Explicit empty-path compatibility

The original `.meta` v1 and all-zero `.names` rows cannot distinguish known empty
paths from unavailable legacy path starts. Therefore this fix uses **metadata
v2** for fresh builds (a reviewed compatibility expansion). The existing 7-column
names shape is unchanged; in v2, a well-formed zero-count row records a known
empty path. Unknown starts are saved as the already-supported 3-column form.
Supplied numeric fields are parsed strictly. v1 long zero-count rows remain
unknown and error on range queries/sidecar rebuilding; legacy paths shorter than
`k` are safely known empty from length. Positive-count path and checkpoint errors
are not converted to empty results. Mixed valid/empty indexes preserve empty
names/lengths through saves, loads, queries, and sampled-sidecar rebuilding,
allowing singleton source accounting by partition seed retention.

This binary reads v1 and v2; old binaries reject v2 as an unsupported version.
Loading, appending to, repairing sidecars for, or re-saving v1 preserves its v1
identity: **none of these operations cleans old fake seeds**. Build from original
input into a **new prefix** to adopt the new extraction policy. Do not relabel a
legacy metadata file as v2. Freeze input hashes, binary hash, parameters, and this
policy in the build manifest; metadata alone is not proof of provenance.

## Seed-free native graphs

All-N/too-short input builds and queries as empty in memory, without fake nodes.
Native `syngBWTwrite` emits no vertices for an empty graph, but `syngBWTread`
terminates with `no Vertex objects in .1gbwt file or can't locate to them`.
Consequently save (including regional save) returns `InvalidInput` for a genuinely
seed-free graph **before creating or truncating any files**. A preloaded unused
dictionary does not count as a graph. Mixed valid/empty indexes are supported.
Previously written empty native files still must not be loaded: this bounded fix
does not change the C reader or add general native-corruption recovery.

Deterministic offline tests cover the full byte alphabet, unfiltered ACGT sampler
controls, long N gaps, repeats, dictionary replay, regional and sidecar round
trips, v1/v2 distinctions, malformed positive starts, seed-free rejection, and
mixed FASTA/AGC CLI builds. No full yeast rebuild is part of these tests.
