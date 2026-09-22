# Succinct Read-Syncmer Index Design

This note describes the production direction for `impg read-index`: a
succinct inverted index from syng syncmer nodes to read ordinals. The
read-walk GBWT prototype is useful as an experiment, but it is not the
right storage model for short-read support queries because it pays GBWT
path overhead for hundreds of millions of tiny read walks.

## Required Operation

The main query is:

```text
syng subgraph / set of syncmer nodes -> reads supporting those nodes
```

The storage model should therefore be node-major:

```text
syncmer node -> sorted read ordinals
```

This is the standard inverted-index problem: a dictionary of terms
syncmer nodes, each with a monotone posting list of document/read IDs.

## Implementation Status

As of May 2026, `impg read-index` has a builder for the first prototype:

```text
impg read-index -a panel.syng -q reads.fq.gz -o sample
```

This writes:

```text
sample.r2s.meta
sample.r2s.sample
sample.r2s.post
```

That prototype is a build artifact, not yet the production query format.
It proves that node-major read postings can be built, but it has not
fixed the fast-load / immediate-query problem yet:

- there is no public `.r2s` query command yet;
- there is no memory-mapped or zero-copy reader yet;
- the builder still collects all `(syncmer_node, read_ordinal)` pairs and
  sorts them before writing;
- read IDs are still FASTQ-order ordinals, so posting-list deltas are
  much less compressible than they should be;
- the sparse `.r2s.sample` file can seek near a node block, but it is
  not the final observed-node dictionary / offset table we want.

So the correction we still need is `.r2s2`: an indexable, node-major
posting-list format with a real loader/query path. The goal is that a
subgraph query loads only metadata and small dictionaries up front, then
decodes only the posting lists for requested syncmer nodes.

This is separate from the syng positional sidecars (`.syng.spos` and
`.syng.pstep`). Those sidecars already use the regular sampled
path-step scheme and are used by syng query/partition/map. They solve
path-coordinate lookup inside the panel index; they do not solve
read-to-syncmer support lookup.

## Current Prototype

The current `.r2s` prototype writes each node block as:

```text
node_delta, posting_count, read_delta_1, read_delta_2, ...
```

with byte varints and a sparse node-offset sample. This is correct, but
not yet succinct.

Measured on one corrected HG002 FASTQ shard against HPRCv2:

```text
reads              1,595,085
retained reads     1,553,547
postings           5,258,318
observed nodes     4,865,284
.r2s.post          25.7 MB
.r2s.sample        0.3 MB
bits/posting       39.2
```

For this shard, a simple `zstd -3` stream compression of `.r2s.post`
reduces 25.7 MB to 9.7 MB. That proves the current byte stream still
contains substantial structure, but a generic compressor is not an
indexable data structure.

## Key Problem

The current read ordinals follow FASTQ order. For a given syncmer node,
the reads supporting that node are therefore almost random over the read
ID universe. Delta coding random read IDs is expensive.

Succinct list encodings help, but the larger win is read-ID reassignment:
assign nearby IDs to reads that share syncmer nodes. This is a known
inverted-index compression technique, usually called document identifier
reassignment.

For this domain the natural reassignment is genomic/syng locality:

1. Map each read to its matched syng syncmer nodes.
2. Choose a stable ordering key, such as the canonical positional key of
   the read's first or best anchor using the syng positional index.
3. Physically store the read collection, or a read locator permutation,
   in this order.
4. Build node postings over the reassigned read IDs.

After reassignment, postings for adjacent syncmer nodes should become
clustered and often run-like, because reads covering the same locus get
nearby IDs.

## Recommended Format

The production `.r2s` format should be an adaptive, partitioned inverted
index:

```text
.r2s.meta      JSON metadata and codec version
.r2s.nodes     observed syng node IDs, Elias-Fano encoded
.r2s.offsets   monotone block offsets, Elias-Fano encoded
.r2s.post      compressed posting-list payloads
.r2s.reads     optional read-order/locator metadata
```

Posting-list payloads should be chosen per list or per partition:

```text
singleton       inline one read ID
tiny list       packed delta varints or Simple8b/StreamVByte
clustered list  partitioned Elias-Fano
dense/repeat    Roaring-like containers or masked repetitive-node channel
```

The first implementation can use fixed-size partitions, then move to
space-optimized variable partitions once correctness and query behavior
are stable.

## Why Partitioned Elias-Fano

Elias-Fano represents a monotone sequence of `m` integers from universe
`u` in roughly:

```text
m * ceil(log2(u / m)) + 2m bits
```

and supports efficient access/search. Partitioned Elias-Fano improves on
plain Elias-Fano by splitting lists into chunks so local clustering is
encoded with a smaller local universe. This matches our expected data
after read-ID reassignment: local clusters of reads supporting nearby
syncmer nodes.

Relevant literature:

- Vigna, "Quasi-Succinct Indices": Elias-Fano based inverted indexes with
  efficient search operations.
- Ottaviano and Venturini, "Partitioned Elias-Fano Indexes": two-level
  partitioned lists with strong compression/time tradeoffs.
- Ding, Attenberg, and Suel, "Scalable techniques for document identifier
  assignment in inverted indexes": read/document ID reassignment improves
  index compression.
- Chambi et al., "Better bitmap performance with Roaring bitmaps":
  Roaring is a strong representation for set operations on dense or
  containerized integer sets.
- Lemire, Boytsov, and Kurz, "SIMD Compression and the Intersection of
  Sorted Integers": SIMD codecs are excellent when decode speed is more
  important than minimum size.

## Expected Size Direction

The current corrected prototype scales to roughly 10 GB per 30x genome.
That is a baseline, not the target.

Expected improvement layers:

```text
block/indexable compression over current stream       ~2x
read-ID reassignment by syng position                 additional large win
partitioned Elias-Fano over clustered postings        additional win
adaptive singleton/tiny/dense encodings               cleanup win
repeat-node masking or separate repetitive channel    workload-dependent win
```

The important point is that "succinct postings over FASTQ-order read IDs"
is not enough. A competitive format must combine a succinct monotone-list
codec with read-ID reassignment.

## Implementation Plan

1. Add a benchmark command that builds several candidate encodings on the
   same read-hit stream and reports bytes/posting, build time, and lookup
   time.
2. Add read-ID reassignment:
   - derive a canonical syng-position key per retained read;
   - sort retained reads by that key;
   - produce a read-order locator so output can recover original file
     ordinal or readpack offset.
3. Implement `.r2s2` with observed-node dictionary plus offset table.
4. Implement adaptive posting codecs:
   - singleton inline;
   - tiny packed list;
   - fixed-partition Elias-Fano;
   - optional Roaring-style container for dense lists.
5. Implement query over `.r2s2`:
   - translate query syncmer nodes to observed-node entries;
   - decode/iterate postings;
   - union or count read ordinals;
   - return read ordinals or extract reads through the locator.
6. Validate on the HG002 shard, then one full 30x sample.

## Non-Goals

The production read-support index should not store read bases and should
not store short reads as GBWT paths by default. A read-walk GBWT can stay
as a separate experimental branch for exploring path-level read queries,
but it is not the compact support index.

The `.r2s` / future `.s2r` read-syncmer postings are working indexes for
read retrieval and subgraph-to-read support queries. They are not the
primary long-term cohort archive for per-sample coverage. For storage at
tens of thousands of samples, the compact archive target is the
node-coverage vector (`impg map -o pack`): exact per-node counts,
internally block-compressed with zstd, random-accessible by syncmer node
ID, and cheap to regenerate from reads when a read-level support index is
needed.

## References

- Sebastiano Vigna. "Quasi-Succinct Indices."
  https://arxiv.org/abs/1206.4300
- Giuseppe Ottaviano and Rossano Venturini. "Partitioned Elias-Fano
  Indexes." https://pages.di.unipi.it/rossano/assets/pdf/papers/SIGIR14.pdf
- Shuai Ding, Josh Attenberg, and Torsten Suel. "Scalable techniques for
  document identifier assignment in inverted indexes."
  https://doi.org/10.1145/1772690.1772723
- Samy Chambi, Daniel Lemire, Owen Kaser, and Robert Godin. "Better
  bitmap performance with Roaring bitmaps." https://arxiv.org/abs/1402.6407
- Daniel Lemire, Leonid Boytsov, and Nathan Kurz. "SIMD Compression and
  the Intersection of Sorted Integers." https://arxiv.org/abs/1401.6399
- Giulio Ermanno Pibiri and Rossano Venturini. "Techniques for Inverted
  Index Compression." https://arxiv.org/abs/1908.10598
