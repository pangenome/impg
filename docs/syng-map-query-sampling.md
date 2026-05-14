# Syng Mapping, Query, and Position Sampling

This note documents the current syng-backed `impg query` / `impg map`
contract and the positional sampling tradeoff we need to resolve.

## Current Map Contract

`impg map` is currently a syncmer-anchor mapper, not a base-level aligner.

For a query sequence, it extracts canonical syng syncmers, looks up matching
syncmer nodes in the syng index, projects sampled hits onto indexed paths via
`.syng.spos`, groups anchors by target path and strand, and chains colinear
anchors. The output is therefore anchor support over a projected interval.

Consequences:

- A read with an internal deletion can still produce one PAF/GAF record if
  anchors on both sides remain colinear. That is not automatically wrong.
- A single output record does not imply a gap-free base-level alignment.
- The gap/deletion is only visible indirectly today, through query span,
  target span, anchor count, and chain structure. We should avoid emitting a
  CIGAR or equivalent unless a real refinement/alignment step has run.
- GAF output should be interpreted as a syncmer-node support path. PAF output
  should be interpreted as a projected target interval from syncmer anchors.

The next mapper tests should pin these cases explicitly: exact reads,
internal deletions, large deletions that should split, reverse-complement
reads, repeats, and reads whose source/target hits fall between sampled
positions.

## Current Query Contract

Indexed-path queries use `.syng.pstep` to find a sampled checkpoint on the
source path, resume GBWT traversal from that checkpoint, collect source
syncmer nodes overlapping the requested interval, and project those nodes to
target path coordinates through `.syng.spos`.

This is the intended single code path for syng-backed query, partition, and
future mapping improvements:

1. Source path/range to source syncmer nodes via `.syng.pstep` plus GBWT walk.
2. Source syncmer nodes to target occurrences via `.syng.spos`.
3. Anchor grouping/chaining/filtering once, shared by callers.

Correctness requires that a query landing between checkpoints still walks from
the nearest usable checkpoint to the exact requested interval and includes
boundary-overlapping syncmers.

## Current Sampling Rule

The current sidecar sampling is deterministic pseudorandom sampling, not
fixed-period sampling.

For each path syncmer occurrence, the builder hashes:

- syncmer node id
- path id
- step index on that path
- base-pair position on that path
- seed

It keeps the occurrence when the low `sample_shift` bits of the mixed hash are
zero. With the default `sample_shift = 8`, the expected sampling rate is
approximately one checkpoint per 256 syncmer occurrences.

This was a reasonable locate-style design for projected MEM/syncmer hits: the
samples are deterministic, reproducible, and spread pseudorandomly across
nodes, paths, and positions instead of being synchronized to the same step
phase on every path.

However, it is not the same as fixed-width path sampling.

## Tradeoff

Pseudorandom sampling is good when the lookup is "given a syncmer/MEM hit,
find representative path positions" because it avoids phase bias and gives an
unbiased subset of occurrences.

It is awkward for "given path coordinate X, jump close to X" because the gap
between checkpoints is geometric rather than bounded. With `sample_shift = 8`,
the expected distance is about 256 syncmer steps, but a particular interval can
land after a much longer unsampled run. Queries remain correct if we always
walk from the checkpoint to the requested range, but latency has a tail.

Fixed-period path-step sampling would give a hard bound: sample every
`2^sample_shift` syncmer steps per path, so a path-coordinate query needs to
scan at most that many syncmer steps from the preceding checkpoint. That is the
better shape for `.syng.pstep` as a path-position checkpoint index.

## Likely Direction

We should separate the concepts if both are needed:

- `.syng.pstep`: fixed-period per-path checkpoints for bounded
  path-coordinate query lookup.
- `.syng.spos`: either keep deterministic pseudorandom occurrence sampling for
  locate-style projection, or derive/store whatever target occurrence coverage
  the mapper/query contract actually requires.

If query and map need complete target occurrence projection, `.syng.spos`
cannot be sparse without accepting missed anchors. In that case it should be
complete or backed by another exact lookup path. Sparse sampling is appropriate
only when the caller explicitly wants sampled representative hits.

Before changing the file format, add tests that distinguish:

- fixed-period checkpoint distance guarantees on `.syng.pstep`
- pseudorandom sampling distribution, if retained for `.syng.spos`
- query correctness for intervals between checkpoints
- map behavior for deletions and reverse-complement reads
