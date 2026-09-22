# Genotype Evidence Pipeline Audit

This document audits the current `impg genotype cos` evidence path and the
related `impg infer` path as implemented in this tree. It is written for
debugging simple cases where an apparently obvious assignment scores oddly.

Code references use current source line numbers.

## Short Answer

`impg genotype cos` is a cosine scorer over typed graph-feature vectors. Today
it has two concrete backends:

- A syng syncmer-node backend fed by a syng index prefix or syng-native render
  bundle plus `--pack` or syng `--proj`.
- A graph-node backend fed by `--graph <local.gfa>`, a graph render bundle, or
  a syng-derived dynamic query GFA plus a typed graph pack or GFA projection
  bundle.

Both backends accept pack files in plain TSV, zstd-compressed TSV, or current
binary pack format, because `pack::read` dispatches on binary magic and
otherwise opens a niffler reader for TSV.

`impg genotype cos` does not directly consume FASTA, FASTQ, BAM, or CRAM as
sample evidence. FASTA/FASTQ reads must first go through `impg map` to produce
a syng pack/projection. Raw GAF read walks are not a direct genotype evidence
argument either; for GFA graphs they must first go through
`impg project --gfa local.gfa --gaf reads.gaf`, and for syng projections the
GAF component is used by `infer` read-link stitching rather than local
`genotype cos` scoring.

`impg map` currently builds this evidence from a syng index plus a query FASTA
or FASTQ file. It can emit `gaf`, `paf`, binary `pack`, text `pack-tsv`, or
`proj` (`src/main.rs:5577`, `src/main.rs:5583`, `src/main.rs:5587`). There is
no BAM, CRAM, or graph-GFA ingestion backend in `impg map`; GAF is an output
read-walk format here, not an input format for building pack coverage.

## CLI Surface Today

### Genotype

Representative commands that work today:

```bash
impg genotype cos -a panel.syng -p sample.pack -r sampleA#0#chr1:100000-102000
impg genotype cos -a panel.syng -p sample.pack.tsv.zst -r sampleA#0#chr1:100000-102000
impg genotype cos -a panel.syng --proj sample.proj -r sampleA#0#chr1:100000-102000
impg genotype cos --render-bundle locus.impg-gbz --pack local.sample.pack --ploidy 2
impg genotype cos --graph local.gfa --pack sample.graph.pack.tsv --pack-feature-space gfa-segment
impg genotype cos --graph local.gfa --proj sample.gfa.proj --target-path REF:0-10000
impg gt cosigt -a panel.syng -p sample.pack -r sampleA#0#chr1:100000-102000
```

The `genotype` command has the alias `gt`, and `cos` has the hidden alias
`cosigt` (`src/main.rs:4275`, `src/main.rs:4277`, `src/main.rs:4929`). The
implementation still emits `#alias\tcosigt` in output metadata
(`src/commands/genotype.rs:438`, `src/commands/genotype.rs:442`).

The graph source is `--index/-a`, `--graph`, `--render-bundle`, or a dynamic
query graph built from `--index`, `--target-range`, and `--sequence-files`.
Render bundles whose manifest feature space is `syng-syncmer-node` route to
the syng backend; bundles whose feature space is `gfa-segment` or
`variation-graph-node` route to the graph backend.

The required sample evidence is `--pack/-p` or `--proj`. In local genotype
scoring, `--proj` resolves to a pack path. Syng projection GAF paths are not
used by `genotype cos`; GFA projection bundles are used as typed graph evidence
because their manifest records the graph pack, graph ID, feature-ID mode,
contribution model, source GAF copy, and read-contribution table.

### Map

Representative commands that build evidence:

```bash
impg map -a panel.syng -q reads.fq -o pack -O sample.pack
impg map -a panel.syng -q reads.fq -o pack-tsv -O sample.pack.tsv.zst
impg map -a panel.syng -q reads.fq -o gaf -O reads.gaf.zst
impg map -a panel.syng -q reads.fq -o proj -O sample.proj
impg map -a panel.syng -q reads.fq -o paf -O reads.paf
```

`map` accepts one query file, documented as FASTA or FASTQ (`src/main.rs:5583`).
The parsers detect `>` or `@` after opening the file with niffler, so common
compression handled by niffler is accepted (`src/main.rs:416`,
`src/main.rs:433`, `src/main.rs:435`, `src/main.rs:962`, `src/main.rs:982`,
`src/main.rs:990`). Anything else errors as "not FASTA or FASTQ"
(`src/main.rs:438`, `src/main.rs:993`).

The output format gate is explicit: `gaf`, `paf`, `proj`, binary pack aliases,
or text pack aliases (`src/main.rs:8724`, `src/main.rs:8732`;
`src/pack.rs:25`, `src/pack.rs:29`). `proj` requires `-O` because it writes a
directory (`src/main.rs:8737`). Binary pack and projection validate zstd
compression level and block size (`src/main.rs:8743`).

There are two syng mapping modes:

- `gaf`, `pack-tsv`, binary `pack`, and `proj` load the lightweight
  `SyngMatcher`, which only requires `.syng.meta` and `.1khash` and is intended
  for read-to-node membership queries (`src/syng.rs:2114`, `src/syng.rs:2144`;
  `src/main.rs:8774`, `src/main.rs:8785`, `src/main.rs:8796`,
  `src/main.rs:8814`).
- `paf` loads the full `SyngIndex` and maps query syncmer anchors to indexed
  paths (`src/main.rs:8833`, `src/main.rs:8835`, `src/syng.rs:4392`).

Missing backends are as important as existing ones: no `impg map` branch accepts
BAM, CRAM, GAM, existing GAF, or arbitrary graph GFA as sample input. Current
pack evidence is built only by scanning read/query sequences for syng syncmer
nodes.

For graph GAF evidence, use `impg project` rather than `impg map`:

```bash
impg project --gfa local.gfa --gaf reads.gaf -O sample.gfa.proj
impg genotype cos --graph local.gfa --proj sample.gfa.proj
```

## Intermediate Products

Existing products:

- `pack-tsv`: human-readable `#node_id\tcount` table, sorted by node id
  (`src/pack.rs:63`). It can be stdout or `.zst/.zstd` if written through the
  map output writer (`src/main.rs:688`, `src/main.rs:691`).
- Binary pack: magic `IMPGPKB1`, dense byte blocks, overflow records for counts
  above 255, metadata for universe nodes, retained reads, and raw syncmer
  anchors (`src/pack.rs:5`, `src/pack.rs:73`, `src/pack.rs:151`,
  `src/pack.rs:154`, `src/pack.rs:156`).
- Projection bundle: directory with `manifest.json`, `sample.pack`, and
  `reads.gaf.zst` (`src/main.rs:1428`, `src/main.rs:1439`,
  `src/main.rs:1440`, `src/main.rs:1476`; `src/projection.rs:37`).
- GAF read walks: `impg map -o gaf` or projection GAF. Paths are signed syng
  syncmer node walks and include `qp:B:I` query syncmer positions
  (`src/main.rs:516`, `src/main.rs:533`, `src/main.rs:539`).
- GFA projection bundle: directory with `manifest.json`, typed
  `sample.pack.tsv`, a source GAF copy, and `read-contributions.tsv`, produced
  by `impg project --gfa ... --gaf ...`.
- Genotype result TSV: top cosine genotype combinations with similarity, QV,
  dot, sample norm, genotype norm, haplotypes, regions, candidate anchors, and
  candidate span fractions (`src/commands/genotype.rs:462`,
  `src/commands/genotype.rs:487`).
- Genotype debug report: `impg genotype cos --emit-report <path>` or
  `--debug-report <path>` writes a sectioned TSV with input metadata, coverage
  summaries, selected sample features, candidate features, top result scores,
  and per-feature score decomposition. The graph backend also reports segment
  lengths, raw counts, normalized weights, graph ID, and contribution model.
- Infer local call TSV: one or more target ranges with local genotype-style
  calls and `PASS`/`NO_CALL` status (`src/commands/infer.rs:448`,
  `src/commands/infer.rs:475`).
- Infer mosaic TSV, and optional FASTA/GFA, when stitching or emit options are
  requested (`src/commands/infer.rs:1267`, `src/commands/infer.rs:1336`,
  `src/commands/infer.rs:1380`).

Missing or weak for human debugging:

- No per-read syng pack contribution table. `reads.gaf.zst` shows every
  syncmer walk, but the syng pack path deduplicates nodes per read before
  counting, so the GAF is not a direct explanation of pack counts.
- No report of candidates removed by `--min-anchors`, `--min-span-fraction`, or
  `--candidate-top-k`.
- No read-length or effective-depth diagnostics. The GFA report exposes segment
  lengths; syng feature lengths are uniform within an index and are not a
  weighting term.
- Infer mosaic exposes read-link/read-emission weights, but not the underlying
  GAF records, GBWT MEM intervals, matched candidate occurrences, or filtered
  candidates that produced those weights.

## Read To Node Coverage

The pack builder is the core "read to coverage" path. For each query record,
`build_syng_map_pack_chunk` does this (`src/main.rs:1196`):

```text
syncmers = matched syng syncmers in the read
node_ids = abs node ids from syncmers
sort node_ids
deduplicate node_ids within this read
if distinct node_ids < --min-anchors: skip the read
for each distinct node id: count[node_id] += 1
```

The exact code clears `node_ids`, extends it from `syncmer.node_id`, sorts and
deduplicates, then applies the `min_anchors` filter to the deduplicated length
(`src/main.rs:1204`, `src/main.rs:1211`, `src/main.rs:1214`,
`src/main.rs:1215`). Each retained distinct node gets `+1`
(`src/main.rs:1219`). The `anchors` statistic counts raw syncmers before
deduplication (`src/main.rs:1222`).

Consequences:

- Node IDs are deduplicated per read before pack counting.
- Repeated occurrences of the same node inside one read do not increase that
  node's pack count.
- Counts are not weighted by node length, syncmer span, read length, base
  coverage, or mapping quality.
- `--min-anchors` for pack output means minimum distinct syng node IDs in the
  read. A read with many repeated occurrences of one node can still fail the
  filter.
- Orientation is ignored in pack counts because the stored field is `node_id`,
  the unsigned absolute syncmer node (`src/syng.rs:1933`,
  `src/syng.rs:2030`).

The GAF path is different. GAF emission retains every matched syncmer occurrence
in query order, including repeated node IDs, and writes orientation as `>` or
`<` from the signed node (`src/main.rs:533`, `src/main.rs:535`). Its
`--min-anchors` check uses raw `syncmers.len()` in the GAF builder
(`src/main.rs:777`, `src/main.rs:778`), not distinct node count. That means a
projection bundle can contain a GAF read walk with repeats that the accompanying
pack has collapsed for coverage.

Read orientation selection happens before both GAF and pack construction.
Syng scans the input read and its reverse complement, keeps the orientation with
more syncmer matches, and negates/repositions reverse-complement syncmers when
that orientation wins (`src/syng.rs:2070`, `src/syng.rs:2082`,
`src/syng.rs:2108`).

## Candidate Vectors And Cosine Scoring

Candidate collection starts from the user target range. `compute_syng_cosigt`
parses `seq:start-end`, validates coordinates, then calls
`collect_syng_candidates` (`src/commands/genotype.rs:379`,
`src/commands/genotype.rs:393`). Candidate discovery is based on syng
homologous intervals from `query_region_with_anchors_ext`
(`src/commands/genotype.rs:165`, `src/syng.rs:4508`).

In syng, that query path is walked over the requested source interval, converted
to a map from absolute node id to query positions and orientation, and then
located across indexed paths (`src/syng.rs:4604`, `src/syng.rs:4608`,
`src/syng.rs:4633`, `src/syng.rs:4643`). Located hits are grouped by target
path and strand (`src/syng.rs:5074`, `src/syng.rs:5078`,
`src/syng.rs:5160`), then overlapping intervals on the same path/strand are
merged and anchors are sorted/deduplicated (`src/syng.rs:5221`,
`src/syng.rs:5297`, `src/syng.rs:5313`).

Candidate modes:

- `spanning` groups all hits by `(path_name, strand)`, expands to min start/max
  end, sorts/deduplicates anchors, and requires both `--min-anchors` and
  `--min-span-fraction` (`src/commands/genotype.rs:201`,
  `src/commands/genotype.rs:220`, `src/commands/genotype.rs:226`,
  `src/commands/genotype.rs:233`).
- `overlapping` scores each query-gathered hit independently and only applies
  `--min-anchors` plus non-empty interval checks; it computes span fraction but
  does not filter on `--min-span-fraction` in that branch
  (`src/commands/genotype.rs:175`, `src/commands/genotype.rs:178`,
  `src/commands/genotype.rs:181`).

For each candidate interval, `syng_candidate_features` walks the indexed path
range and counts every syncmer step by unsigned node id
(`src/commands/genotype.rs:116`, `src/commands/genotype.rs:138`,
`src/commands/genotype.rs:139`). Repeated node occurrences in the candidate
range therefore increase the candidate feature count. The candidate also keeps
an oriented signed walk for infer read-walk scoring (`src/commands/genotype.rs:28`,
`src/commands/genotype.rs:140`), but cosine scoring uses the unoriented absolute
feature vector (`src/commands/genotype.rs:33`, `src/commands/genotype.rs:145`).

The candidate's `strand` is retained for reporting and stitching
(`src/commands/genotype.rs:25`, `src/commands/infer.rs:1108`). It does not make
cosine features signed. A candidate on `-` still contributes counts by absolute
node id, while FASTA emission later reverse-complements sequence if needed
(`src/commands/infer.rs:1321`).

Ranking and combination scoring:

- Build a locus feature universe from all candidates (`src/commands/genotype.rs:274`;
  `src/genotyping.rs:96`).
- Compute a one-haplotype cosine similarity for each candidate using sample
  counts restricted to that universe (`src/commands/genotype.rs:279`,
  `src/commands/genotype.rs:281`; `src/genotyping.rs:119`).
- Sort candidate preselection by single similarity, then anchor count, path
  name, and start (`src/commands/genotype.rs:288`). `--candidate-top-k` then
  truncates candidates when greater than zero (`src/commands/genotype.rs:295`).
- Recompute selected features from the retained candidates and reject zero
  sample norm over those features (`src/commands/genotype.rs:299`,
  `src/commands/genotype.rs:304`, `src/commands/genotype.rs:305`).
- Enumerate nondecreasing combinations with replacement for the requested
  ploidy (`src/genotyping.rs:210`, `src/genotyping.rs:231`,
  `src/genotyping.rs:233`). This permits homozygous combinations like
  candidate 0 + candidate 0.
- Sum candidate feature counts into a genotype vector
  (`src/genotyping.rs:147`, `src/genotyping.rs:149`).
- Score cosine as `dot(sample, genotype) / (||sample|| * ||genotype||)`
  (`src/genotyping.rs:154`, `src/genotyping.rs:162`,
  `src/genotyping.rs:167`).
- Compute QV as `999` for similarity at least 1, `0` for similarity at most 0,
  otherwise `-10 * log10(1 - similarity)` (`src/genotyping.rs:169`).
- Sort results by similarity descending, dot descending, and combination index
  order, then truncate to `--top-n` (`src/genotyping.rs:201`,
  `src/commands/genotype.rs:323`).

There is no node length weighting, read length normalization, per-read
likelihood, mapping quality, or orientation-aware cosine component in the
current genotype score.

## Infer Versus Genotype

`infer` reuses genotype's local cosine scorer. `compute_local_call_sets` creates
a `SyngCosigtQuery` for each target and calls `genotype::compute_syng_cosigt`
(`src/commands/infer.rs:411`, `src/commands/infer.rs:420`,
`src/commands/infer.rs:432`). Therefore local infer calls have the same pack
coverage, candidate collection, candidate top-k, combination search, cosine, and
QV behavior as `genotype cos`.

What infer adds:

- Multiple targets from one range, BED, partition BED, or internal partition
  discovery (`src/main.rs:4963`, `src/main.rs:4967`, `src/main.rs:4971`,
  `src/main.rs:4975`; `src/commands/infer.rs:356`).
- Optional phase-block splitting before local calls (`src/commands/infer.rs:378`,
  `src/commands/infer.rs:1440`).
- Beam stitching over ordered local genotype states (`src/commands/infer.rs:591`,
  `src/commands/infer.rs:1181`).
- Optional GAF read-walk evidence for stitched emissions/transitions
  (`src/main.rs:4959`, `src/main.rs:7982`, `src/commands/infer.rs:1447`,
  `src/commands/infer.rs:1454`).
- Optional mosaic, FASTA, and GFA outputs (`src/main.rs:5092`,
  `src/main.rs:5096`, `src/main.rs:5100`; `src/commands/infer.rs:1487`,
  `src/commands/infer.rs:1490`, `src/commands/infer.rs:1505`).

If `--proj` is used, infer automatically resolves both `sample.pack` and
`reads.gaf.zst` when the manifest contains GAF (`src/main.rs:7972`,
`src/main.rs:7982`; `src/projection.rs:102`). If `--pack` is used, `--gaf` can
be supplied separately. Local call scoring still uses pack coverage; GAF read
walks are only built into evidence when stitching or an emit option requests a
mosaic/sequence path (`src/commands/infer.rs:1454`,
`src/commands/infer.rs:1459`). The local infer header advertises read walks
when a GAF path is present (`src/commands/infer.rs:455`), but the GAF does not
change the local cosine rows unless the stitched mosaic path is active.

Read-walk evidence has two modes:

- With `qp:B:I` query positions, infer builds signed read walk steps, computes
  syng GBWT MEMs on the forward read walk and on a reversed walk, and scores
  candidate overlaps from those MEMs (`src/commands/infer.rs:682`,
  `src/commands/infer.rs:713`, `src/commands/infer.rs:731`,
  `src/commands/infer.rs:1071`, `src/commands/infer.rs:1072`,
  `src/commands/infer.rs:1075`; `src/syng.rs:4016`).
- Without query positions, infer falls back to whole-walk orientation hits using
  a longest-increasing-subsequence length over candidate positions, checking
  both read orientations (`src/commands/infer.rs:831`,
  `src/commands/infer.rs:844`, `src/commands/infer.rs:871`,
  `src/commands/infer.rs:1083`).

The candidate walk index is signed and occurrence-aware: every candidate
`oriented_walk` step is inserted by signed node, with its position in the
candidate walk (`src/commands/infer.rs:1002`, `src/commands/infer.rs:1012`).
For MEM mode, each MEM accumulates best overlap per candidate and adds anchor
counts to `walk_counts` (`src/commands/infer.rs:958`,
`src/commands/infer.rs:965`, `src/commands/infer.rs:988`). Hits are filtered by
`--min-read-link-anchors`, grouped by call, and used for both candidate support
and adjacent-call links (`src/commands/infer.rs:749`,
`src/commands/infer.rs:756`, `src/commands/infer.rs:801`,
`src/commands/infer.rs:770`). Multiple hits in one call split read/support
weight across candidates (`src/commands/infer.rs:818`,
`src/commands/infer.rs:825`), and link weights are split across adjacent hit
pairs (`src/commands/infer.rs:781`, `src/commands/infer.rs:793`).

Rewards are simple log-scaled anchor weights:
`read_link_weight * 10 * log10(1 + anchor_weight)` (`src/commands/infer.rs:1094`).
Stitching adds local QV plus unique candidate read-emission rewards
(`src/commands/infer.rs:601`, `src/commands/infer.rs:615`,
`src/commands/infer.rs:624`), subtracts same-path/recombine transition costs,
and adds read-link rewards across adjacent local calls
(`src/commands/infer.rs:1102`, `src/commands/infer.rs:1146`,
`src/commands/infer.rs:1167`, `src/commands/infer.rs:1218`).

## Likely Failure Modes In Weird Simple Assignments

1. Pack dedup versus candidate repeated counts.

   Sample pack coverage deduplicates a node within each read
   (`src/main.rs:1211`, `src/main.rs:1214`), but candidate vectors count every
   occurrence of a repeated node in the candidate range
   (`src/commands/genotype.rs:138`, `src/commands/genotype.rs:139`). A
   candidate with repeated copies can therefore have a genotype vector that
   expects multiplicity the sample vector may not record when reads span
   repeated occurrences of the same syncmer node. Projection GAF will still show
   those repeats because GAF emission writes every syncmer occurrence
   (`src/main.rs:533`).

2. No node-length or read-length normalization.

   A retained read contributes `+1` to each distinct node it contains,
   independent of read length and independent of how much sequence each node
   represents (`src/main.rs:1219`). Candidate and genotype norms are pure count
   vector norms (`src/genotyping.rs:154`, `src/genotyping.rs:158`). Long regions,
   short regions, high-density syncmer regions, and repeated-node regions are
   all compared in raw syncmer-node count space.

3. `--min-anchors` can mean different things in pack and GAF.

   Pack `--min-anchors` uses distinct node count after per-read dedup
   (`src/main.rs:1211`, `src/main.rs:1215`). GAF `--min-anchors` uses raw
   syncmer occurrence count (`src/main.rs:777`, `src/main.rs:778`). A read can
   appear in GAF but not add pack coverage if its repeated path has too few
   distinct nodes.

4. Orientation is mostly absent from genotype.

   Candidate discovery records strand and infer read-walk evidence uses signed
   walks, but genotype cosine features are unsigned node IDs
   (`src/commands/genotype.rs:139`, `src/commands/genotype.rs:145`). A simple
   inversion-like or orientation-ambiguous case can therefore look identical to
   cosine unless candidate range collection or span filters separate it.

5. Candidate preselection can hide the right ploidy combination.

   `--candidate-top-k` filters by single-haplotype cosine before combination
   search (`src/commands/genotype.rs:281`, `src/commands/genotype.rs:295`).
   The correct combination can be impossible to recover if one member is weak as
   a singleton but necessary in combination.

6. Spanning versus overlapping changes filters materially.

   Spanning mode applies `--min-span-fraction`; overlapping mode does not
   (`src/commands/genotype.rs:178`, `src/commands/genotype.rs:233`). A debug run
   should compare both modes before concluding the evidence itself is wrong.

7. Syng projection can make GAF look relevant to local genotype when it is not.

   `sample.proj` carries both `sample.pack` and `reads.gaf.zst`, but genotype
   only consumes the pack path for syng local cosine scoring. GAF read-walk
   order affects infer stitching, not local `genotype cos` rows.

8. Graph and syng projection bundles have different semantics.

   GFA projection bundles are typed graph evidence and include
   `read-contributions.tsv`; syng projection bundles are syng pack evidence
   plus optional GAF read walks for infer. They share the `--proj` flag but not
   the same feature space or read-link behavior.

## Current Tests

Unit tests in `src/genotyping.rs` cover stable user-visible names, combination
search with replacement genotypes, and max-combination budget errors
(`src/genotyping.rs:263`, `src/genotyping.rs:271`, `src/genotyping.rs:290`).

Integration tests in `tests/test_syng_integration.rs` cover most of the CLI
surface:

- Render bundle genotype/infer on syng-native bundles
  (`tests/test_syng_integration.rs:416`, `tests/test_syng_integration.rs:451`).
- `impg map` GAF default, pack-tsv, compressed pack-tsv, binary pack,
  projection bundle, reverse-complement GAF, PAF, khash-only GAF/pack, and
  deterministic parallel GAF output (`tests/test_syng_integration.rs:960`,
  `tests/test_syng_integration.rs:1035`, `tests/test_syng_integration.rs:1231`,
  `tests/test_syng_integration.rs:1277`, `tests/test_syng_integration.rs:1310`,
  `tests/test_syng_integration.rs:1350`, `tests/test_syng_integration.rs:1405`).
- The map pack test specifically asserts one output row per distinct GAF node
  and count `1` for a single read (`tests/test_syng_integration.rs:1062`,
  `tests/test_syng_integration.rs:1083`).
- Genotype CLI permutations across `genotype`/`gt`, `cos`/`cosigt`,
  pack-tsv, compressed pack-tsv, binary pack, spanning/overlapping modes,
  `--proj`, ploidy/top-n options, output compression, and expected errors
  (`tests/test_syng_integration.rs:1463`, `tests/test_syng_integration.rs:1647`,
  `tests/test_syng_integration.rs:1769`, `tests/test_syng_integration.rs:1819`,
  `tests/test_syng_integration.rs:1928`).
- Short-read simulated heterozygote with an unsampled decoy
  (`tests/test_syng_integration.rs:1978`, `tests/test_syng_integration.rs:2044`,
  `tests/test_syng_integration.rs:2084`).
- Infer from projection over explicit range, target BED, partition BED,
  discovery mode, and mosaic/FASTA/GFA emission
  (`tests/test_syng_integration.rs:2147`, `tests/test_syng_integration.rs:2245`,
  `tests/test_syng_integration.rs:2297`, `tests/test_syng_integration.rs:2352`,
  `tests/test_syng_integration.rs:2449`).
- Read-walk stitching rewards, disabled read-link weight, and phase-block
  stitching (`tests/test_syng_integration.rs:2496`,
  `tests/test_syng_integration.rs:2604`, `tests/test_syng_integration.rs:2662`,
  `tests/test_syng_integration.rs:2726`).
- Repeated-copy and triplicated-copy infer cases
  (`tests/test_syng_integration.rs:2809`, `tests/test_syng_integration.rs:2917`,
  `tests/test_syng_integration.rs:3047`).
- Noisy low-coverage nested SV, read-walk order decoy, and paralogous swapped
  copies (`tests/test_syng_integration.rs:3207`,
  `tests/test_syng_integration.rs:3395`, `tests/test_syng_integration.rs:3591`).

Syng unit tests cover the GBWT MEM helper used by infer read-walk scoring,
including offset-sensitive full-length MEMs (`src/syng.rs:5694`,
`src/syng.rs:5717`, `src/syng.rs:5734`).

## Validation Gaps

Compared with a LikeGT/COSIGT-style validation suite, the current tests are
useful but not systematic. Main gaps:

- No matrix of known-truth simulations varying coverage, read length, error
  rate, ploidy, copy number, allele divergence, repeat content, and decoy
  density.
- No explicit genotype test for the repeated-node mismatch where a long read's
  repeated GAF nodes collapse in pack coverage but candidate vectors retain copy
  multiplicity.
- No direct tests for rejecting or ignoring FASTA/FASTQ/BAM/CRAM as genotype
  evidence inputs beyond the current CLI shape.
- No benchmark against external LikeGT/COSIGT expected calls or score
  distributions.
- No tests for read-length normalization because the implementation has no such
  normalization today.
- Candidate filter-rejection debug tables are still missing, even though the
  current debug report exposes selected sample vectors, candidate vectors, and
  score decomposition.

## Practical Debug Recipe

For a small weird assignment today:

1. Build both pack and GAF from the same reads.

   ```bash
   impg map -a panel.syng -q reads.fq -o pack-tsv -O sample.pack.tsv
   impg map -a panel.syng -q reads.fq -o gaf -O reads.gaf
   ```

2. Check whether nodes repeated in `reads.gaf` appear only once per read in the
   pack table. That is expected.

3. Run genotype with no candidate prefilter and relaxed span filter.

   ```bash
   impg genotype cos -a panel.syng -p sample.pack.tsv \
     -r sampleA#0#chr1:100000-102000 \
     --candidate-top-k 0 --top-n 20 --min-span-fraction 0 \
     --emit-report genotype.report.tsv
   ```

4. Compare `spanning` and `overlapping`.

   ```bash
   impg genotype cos -a panel.syng -p sample.pack.tsv \
     -r sampleA#0#chr1:100000-102000 --candidate-mode overlapping \
     --candidate-top-k 0 --top-n 20 --min-anchors 1 \
     --emit-report genotype.overlap.report.tsv
   ```

5. If read order or phase should matter, use infer with a projection or explicit
   GAF and inspect the mosaic read-link/read-emission columns.

   ```bash
   impg map -a panel.syng -q reads.fq -o proj -O sample.proj
   impg infer -a panel.syng --proj sample.proj \
     --partitions partitions.bed --stitch beam --emit-mosaic mosaic.tsv \
     --candidate-top-k 0 --top-n 20 --read-link-weight 5
   ```

The key remaining limitation is per-read syng pack provenance: `--emit-report`
now exposes selected sample vectors, candidate vectors, and feature-level score
decomposition, but syng pack counts still need manual correlation back to
`reads.gaf.zst` when repeated read-walk nodes collapsed during pack creation.
