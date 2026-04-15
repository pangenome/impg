//! FastLocate: r-index locate structure for a classical GBWT.
//!
//! Port of jltsiren/gbwt's C++ FastLocate (`fast_locate.{h,cpp}`) to Rust on
//! top of the `gbz` crate (which, despite the crate name on crates.io, *is*
//! `jltsiren/gbwt-rs`). Builds a succinct locate structure whose space is
//! O(r) where r is the total number of BWT runs across all records.
//!
//! Algorithm: Gagie, Navarro, Prezza (JACM 2020), as adapted to multi-string
//! BWTs in Sirén, Garrison, Novak, Paten, Durbin (Bioinformatics 2020).
//!
//! Text positions are packed as `seq_id * max_length + seq_offset`. GBWT
//! indexes the *reverse* paths, so `seq_offset` is the distance to the
//! sequence's endmarker (i.e. distance to the end of the path when walked
//! forward from the GBWT's first position).

use gbz::GBWT;
use gbz::bwt::Pos;
use gbz::support::RLEIter;
use gbz::ENDMARKER;
use simple_sds::ops::PredSucc;
use simple_sds::serialize::Serialize as SdsSerialize;
use simple_sds::sparse_vector::{SparseBuilder, SparseVector};
use std::io::{self, Read, Write};

/// Packed text position = `seq_id * max_length + seq_offset`.
pub type PackedPos = u64;

/// r-index locate structure built over a classical GBWT.
///
/// The structure stores, per BWT run:
///   - one *head* sample: the packed text-position of the first visit in
///     the run (indexed into `samples` by global run id)
///   - one *tail* sample: the packed text-position of the last visit in
///     the run; the set of all tail samples is stored as a sparse bitvector
///     (`last`) with a parallel `last_to_run[]` giving, for a rank-of-1,
///     the global run id whose tail landed at that position.
///
/// `locate_next(prev)` then implements the φ⁻¹ jump from one visit of a
/// node to the next (in text-position order) in `O(log r)`.
pub struct FastLocate {
    /// `samples[global_run_id]` = packed text-position at the START of that
    /// BWT run (head sample), sorted by `global_run_id`.
    samples: Vec<PackedPos>,
    /// Sparse bitvector over the text-position universe (`sequences * max_length`),
    /// with a set bit at every run's *tail* packed position.
    last: SparseVector,
    /// `last_to_run[k]` = global run id of the `k`-th set bit of `last`
    /// (where `k` is the 0-based rank in text-position order).
    last_to_run: Vec<u32>,
    /// `comp_to_run[record_id]` = first global run id belonging to that
    /// record. Lets us compute `globalRunId(node, local_run) = comp_to_run[record_id(node)] + local_run`.
    comp_to_run: Vec<u32>,
    /// Total BWT runs across all records (== `samples.len()`).
    total_runs: usize,
    /// Packing universe: every path has strictly fewer than `max_length` nodes.
    max_length: u64,
    /// Number of sequences in the source GBWT. Kept for future use when
    /// wiring a bp-offset sidecar keyed by `(seq_id, forward_node_idx)`.
    #[allow(dead_code)]
    sequences: usize,
}

impl FastLocate {
    /// Packs a `(seq_id, seq_offset)` pair into a single text position.
    #[inline]
    fn pack(&self, seq_id: usize, seq_offset: usize) -> PackedPos {
        seq_id as u64 * self.max_length + seq_offset as u64
    }

    /// Sequence id component of a packed position.
    #[inline]
    pub fn seq_id(&self, p: PackedPos) -> usize {
        (p / self.max_length) as usize
    }

    /// Sequence offset component of a packed position.
    ///
    /// Note: this is the *reverse* offset — distance to the end of the
    /// sequence — because GBWT indexes reverse paths.
    #[inline]
    pub fn seq_offset(&self, p: PackedPos) -> usize {
        (p % self.max_length) as usize
    }

    /// Number of BWT runs stored.
    pub fn total_runs(&self) -> usize {
        self.total_runs
    }

    /// Value of `max_length` used for packing (caller must preserve it to
    /// interpret packed positions).
    pub fn max_length(&self) -> u64 {
        self.max_length
    }

    /// Builds the FastLocate over `gbwt` by walking every sequence once and
    /// collecting head/tail samples at BWT run boundaries.
    ///
    /// Mirrors `FastLocate::FastLocate(const GBWT& source)` in the C++ source.
    pub fn build(gbwt: &GBWT) -> Self {
        let effective = gbwt.effective_size();
        if effective == 0 || gbwt.sequences() == 0 {
            // Empty GBWT — return a zero-state FastLocate that satisfies the
            // invariants used by `decompress_da` (which checks has_node first).
            let empty_last = SparseVector::try_from_iter(std::iter::empty::<usize>())
                .expect("SparseVector::try_from_iter(empty)");
            return FastLocate {
                samples: Vec::new(),
                last: empty_last,
                last_to_run: Vec::new(),
                comp_to_run: Vec::new(),
                total_runs: 0,
                max_length: 1,
                sequences: 0,
            };
        }

        // Phase 1: count runs per record and build comp_to_run as a prefix sum.
        //
        // In gbz 0.6.0 the BWT always has `effective_size` records numbered
        // `0..effective_size`. Record 0 is the endmarker record regardless of
        // `alphabet_offset`; we derive its run count from grouping the
        // sequences by their first non-endmarker node (not from iterating
        // bwt.record(0)), because record 0's internal layout is implementation-
        // specific. Records `1..effective_size` correspond to real nodes,
        // one per node id in `[first_node, alphabet_size)`.
        let mut comp_to_run: Vec<u32> = Vec::with_capacity(effective);
        let mut total_runs: usize = 0;
        for comp in 0..effective {
            comp_to_run.push(total_runs as u32);
            if comp == 0 {
                // Endmarker slot — count via the start() grouping.
                total_runs += count_endmarker_runs(gbwt);
                continue;
            }
            let Some((edges_bytes, bwt_bytes)) = gbwt.as_ref().compressed_record(comp) else {
                continue;
            };
            let Some((edges, _)) = decompress_edges(edges_bytes) else {
                continue;
            };
            let sigma = edges.len();
            let mut runs: usize = 0;
            // Important: runs whose successor is ENDMARKER contribute one
            // logical run PER POSITION, not one per physical RLE run. This
            // mirrors the C++ `LFLoop` convention: each endmarker entry is
            // effectively a distinct run because each ends a different
            // sequence.
            for run in RLEIter::with_sigma(bwt_bytes, sigma) {
                if edges[run.value].node == ENDMARKER {
                    runs += run.len;
                } else {
                    runs += 1;
                }
            }
            total_runs += runs;
        }
        // Endmarker runs always live at `comp_to_run[0]`.
        let endmarker_base: u32 = comp_to_run[0];

        // Phase 2: endmarker runs — one per group of consecutive sequences
        // sharing the same first node.
        let sequences = gbwt.sequences();
        let mut endmarker_runs: Vec<u32> = vec![0; sequences];
        if sequences > 0 {
            let mut run_id: u32 = 0;
            let mut prev_first = first_node_of_seq(gbwt, 0);
            endmarker_runs[0] = run_id;
            for (i, slot) in endmarker_runs.iter_mut().enumerate().skip(1) {
                let curr_first = first_node_of_seq(gbwt, i);
                if curr_first == ENDMARKER || curr_first != prev_first {
                    run_id += 1;
                    prev_first = curr_first;
                }
                *slot = run_id;
            }
        }

        // Allocate sample buffers.
        let mut head_samples: Vec<SampleRecord> = Vec::with_capacity(total_runs);
        let mut tail_samples: Vec<SampleRecord> = Vec::with_capacity(total_runs);

        // Phase 3: per-sequence walk.
        let mut max_length: u64 = 1;
        for i in 0..sequences {
            let mut seq_offset: usize = 0;
            let em_run = endmarker_runs[i];
            let em_global = endmarker_base + em_run;

            // Head sample for the endmarker run if this is the first sequence
            // of its group.
            if i == 0 || em_run != endmarker_runs[i - 1] {
                head_samples.push(SampleRecord { seq_id: i as u32, seq_offset: seq_offset as u64, run_id: em_global });
            }
            // Tail sample for the endmarker run if this is the last sequence
            // of its group.
            if i + 1 >= sequences || em_run != endmarker_runs[i + 1] {
                tail_samples.push(SampleRecord { seq_id: i as u32, seq_offset: seq_offset as u64, run_id: em_global });
            }

            // Walk the sequence, tracking offset and capturing head/tail at
            // each run boundary within the records visited.
            let Some(mut curr) = gbwt.start(i) else {
                // Empty sequence — only the endmarker sample applies.
                seq_offset = 1;
                let total = seq_offset as u64;
                // Flip seq_offsets for samples of this sequence.
                flip_seq_offsets(&mut head_samples, i as u32, total);
                flip_seq_offsets(&mut tail_samples, i as u32, total);
                if total > max_length { max_length = total; }
                continue;
            };
            seq_offset += 1;

            while curr.node != ENDMARKER {
                let record_id = gbwt.node_to_record(curr.node);
                let Some((edges_bytes, bwt_bytes)) = gbwt.as_ref().compressed_record(record_id) else {
                    break;
                };
                let Some((edges, _)) = decompress_edges(edges_bytes) else {
                    break;
                };
                let sigma = edges.len();
                let Some(lf_info) = lf_with_run(&edges, bwt_bytes, sigma, curr.offset) else {
                    break;
                };

                let global_rid = comp_to_run[record_id] + lf_info.run_idx as u32;
                let is_head = curr.offset == lf_info.run_start;
                let is_tail = curr.offset == lf_info.run_end;
                if is_head {
                    head_samples.push(SampleRecord {
                        seq_id: i as u32,
                        seq_offset: seq_offset as u64,
                        run_id: global_rid,
                    });
                }
                if is_tail {
                    tail_samples.push(SampleRecord {
                        seq_id: i as u32,
                        seq_offset: seq_offset as u64,
                        run_id: global_rid,
                    });
                }

                curr = lf_info.next_pos;
                seq_offset += 1;
            }

            // Flip seq_offsets so they mean "distance to end of sequence."
            let total = seq_offset as u64;
            flip_seq_offsets(&mut head_samples, i as u32, total);
            flip_seq_offsets(&mut tail_samples, i as u32, total);
            if total > max_length {
                max_length = total;
            }
        }

        // Phase 4: store samples[] sorted by run_id, and last / last_to_run
        // sorted by text position.
        head_samples.sort_by_key(|s| s.run_id);
        assert_eq!(
            head_samples.len(),
            total_runs,
            "head_samples count mismatch: expected {} runs, got {}",
            total_runs,
            head_samples.len()
        );

        // Build samples[] in run-id order.
        let mut samples: Vec<PackedPos> = Vec::with_capacity(total_runs);
        {
            let tmp = FastLocate {
                samples: Vec::new(),
                last: SparseVector::try_from_iter(std::iter::empty::<usize>()).unwrap(),
                last_to_run: Vec::new(),
                comp_to_run: comp_to_run.clone(),
                total_runs,
                max_length,
                sequences,
            };
            for s in &head_samples {
                samples.push(tmp.pack(s.seq_id as usize, s.seq_offset as usize));
            }
        }

        // Sort tail_samples by (seq_id, seq_offset).
        tail_samples.sort_by(|a, b| (a.seq_id, a.seq_offset).cmp(&(b.seq_id, b.seq_offset)));
        assert_eq!(
            tail_samples.len(),
            total_runs,
            "tail_samples count mismatch: expected {} runs, got {}",
            total_runs,
            tail_samples.len()
        );

        // Build the sparse bitvector over the text-position universe.
        let universe = sequences as u64 * max_length;
        let mut builder = SparseBuilder::new(universe as usize, total_runs)
            .expect("SparseBuilder::new failed");
        let mut last_to_run: Vec<u32> = Vec::with_capacity(total_runs);
        for s in &tail_samples {
            let packed = s.seq_id as u64 * max_length + s.seq_offset;
            builder.set(packed as usize);
            last_to_run.push(s.run_id);
        }
        let last = SparseVector::try_from(builder).expect("SparseVector::try_from failed");

        FastLocate {
            samples,
            last,
            last_to_run,
            comp_to_run,
            total_runs,
            max_length,
            sequences,
        }
    }

    /// Enumerates every visit of `node` as a `(seq_id, seq_offset_from_end)`
    /// pair. This mirrors C++ `decompressDA`.
    ///
    /// The returned `seq_offset` is the distance to the end of the path (not
    /// the forward node index). Callers that need forward bp positions must
    /// look them up through a per-sequence sidecar keyed by the *forward*
    /// node index = `path_node_count - 1 - seq_offset_from_end`.
    pub fn decompress_da(&self, gbwt: &GBWT, node: usize) -> Vec<(usize, usize)> {
        if !gbwt.has_node(node) {
            return Vec::new();
        }
        let record_id = gbwt.node_to_record(node);
        let Some(rec) = gbwt.as_ref().record(record_id) else {
            return Vec::new();
        };
        let len = rec.len();
        if len == 0 {
            return Vec::new();
        }

        let first_run = self.comp_to_run[record_id] as usize;
        let mut result: Vec<(usize, usize)> = Vec::with_capacity(len);
        let mut curr = self.samples[first_run];
        result.push((self.seq_id(curr), self.seq_offset(curr)));
        for _ in 1..len {
            curr = self.locate_next(curr);
            result.push((self.seq_id(curr), self.seq_offset(curr)));
        }
        result
    }

    /// Returns the number of visits to `node` (== `record.len()`).
    pub fn node_count(&self, gbwt: &GBWT, node: usize) -> usize {
        if !gbwt.has_node(node) {
            return 0;
        }
        let record_id = gbwt.node_to_record(node);
        gbwt.as_ref()
            .record(record_id)
            .map(|r| r.len())
            .unwrap_or(0)
    }

    /// The r-index φ⁻¹ step: given the packed text position of one visit,
    /// return the packed text position of the "next" visit (the one
    /// immediately after it in the BWT permutation cycle for that node).
    fn locate_next(&self, prev: PackedPos) -> PackedPos {
        // Find the predecessor 1-bit in `last` at or before `prev`.
        let mut iter = self.last.predecessor(prev as usize);
        let (rank, pred_pos) = iter.next().expect("locate_next: no predecessor in `last`");
        let run = self.last_to_run[rank] as usize;
        // samples[run + 1] is the head of the NEXT run in global run-id order.
        // Distance from the predecessor tail to `prev` is preserved across φ⁻¹.
        self.samples[run + 1] + (prev - pred_pos as u64)
    }

    // ─── serialization ──────────────────────────────────────────────
    //
    // Custom binary format (not compatible with the C++ `.ri` SDSL layout).
    // We mix simple-sds `Serialize` for the `SparseVector` with our own
    // little-endian framing for the `Vec<...>` fields and scalar header.
    //
    // Layout:
    //   [u64 MAGIC]
    //   [u64 version]
    //   [u64 total_runs]
    //   [u64 max_length]
    //   [u64 sequences]
    //   [u64 samples.len()]    [u64 * samples.len()]
    //   [u64 last_to_run.len()] [u32 * last_to_run.len()]
    //   [u64 comp_to_run.len()] [u32 * comp_to_run.len()]
    //   [SparseVector last via simple-sds Serialize]

    const MAGIC: u64 = 0x494D50_46415354; // "IMPFAST"
    const VERSION: u64 = 1;

    pub fn save<W: Write>(&self, w: &mut W) -> io::Result<()> {
        write_u64(w, Self::MAGIC)?;
        write_u64(w, Self::VERSION)?;
        write_u64(w, self.total_runs as u64)?;
        write_u64(w, self.max_length)?;
        write_u64(w, self.sequences as u64)?;
        write_u64_slice(w, &self.samples)?;
        write_u32_slice(w, &self.last_to_run)?;
        write_u32_slice(w, &self.comp_to_run)?;
        self.last
            .serialize(w)
            .map_err(|e| io::Error::other(format!("SparseVector::serialize: {}", e)))?;
        Ok(())
    }

    pub fn load<R: Read>(r: &mut R) -> io::Result<Self> {
        let magic = read_u64(r)?;
        if magic != Self::MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("FastLocate: bad magic 0x{:x}", magic),
            ));
        }
        let version = read_u64(r)?;
        if version != Self::VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("FastLocate: unsupported version {}", version),
            ));
        }
        let total_runs = read_u64(r)? as usize;
        let max_length = read_u64(r)?;
        let sequences = read_u64(r)? as usize;
        let samples = read_u64_vec(r)?;
        let last_to_run = read_u32_vec(r)?;
        let comp_to_run = read_u32_vec(r)?;
        let last = SparseVector::load(r)
            .map_err(|e| io::Error::other(format!("SparseVector::load: {}", e)))?;
        if samples.len() != total_runs {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FastLocate: samples.len() != total_runs",
            ));
        }
        Ok(Self {
            samples,
            last,
            last_to_run,
            comp_to_run,
            total_runs,
            max_length,
            sequences,
        })
    }
}

fn write_u64<W: Write>(w: &mut W, v: u64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_u64<R: Read>(r: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn write_u64_slice<W: Write>(w: &mut W, v: &[u64]) -> io::Result<()> {
    write_u64(w, v.len() as u64)?;
    let mut buf = vec![0u8; v.len() * 8];
    for (i, &x) in v.iter().enumerate() {
        buf[i * 8..(i + 1) * 8].copy_from_slice(&x.to_le_bytes());
    }
    w.write_all(&buf)
}

fn read_u64_vec<R: Read>(r: &mut R) -> io::Result<Vec<u64>> {
    let n = read_u64(r)? as usize;
    let mut buf = vec![0u8; n * 8];
    r.read_exact(&mut buf)?;
    let mut v = Vec::with_capacity(n);
    for chunk in buf.chunks_exact(8) {
        v.push(u64::from_le_bytes(chunk.try_into().unwrap()));
    }
    Ok(v)
}

fn write_u32_slice<W: Write>(w: &mut W, v: &[u32]) -> io::Result<()> {
    write_u64(w, v.len() as u64)?;
    let mut buf = vec![0u8; v.len() * 4];
    for (i, &x) in v.iter().enumerate() {
        buf[i * 4..(i + 1) * 4].copy_from_slice(&x.to_le_bytes());
    }
    w.write_all(&buf)
}

fn read_u32_vec<R: Read>(r: &mut R) -> io::Result<Vec<u32>> {
    let n = read_u64(r)? as usize;
    let mut buf = vec![0u8; n * 4];
    r.read_exact(&mut buf)?;
    let mut v = Vec::with_capacity(n);
    for chunk in buf.chunks_exact(4) {
        v.push(u32::from_le_bytes(chunk.try_into().unwrap()));
    }
    Ok(v)
}

//----------------------------- helpers ---------------------------------

#[derive(Clone, Copy, Debug)]
struct SampleRecord {
    seq_id: u32,
    seq_offset: u64,
    run_id: u32,
}

/// Flip `seq_offset` of all samples matching `seq_id` so they encode the
/// distance to the end of the sequence instead of the distance from the
/// start. Called once per sequence after the walk finishes.
fn flip_seq_offsets(samples: &mut [SampleRecord], seq_id: u32, total: u64) {
    // Iterate in reverse because only the tail of the Vec contains samples
    // belonging to this sequence (walks append sequentially).
    for s in samples.iter_mut().rev() {
        if s.seq_id != seq_id {
            break;
        }
        s.seq_offset = total - 1 - s.seq_offset;
    }
}

/// First node of sequence `seq_id` (== `ENDMARKER` if empty).
fn first_node_of_seq(gbwt: &GBWT, seq_id: usize) -> usize {
    gbwt.start(seq_id).map(|p| p.node).unwrap_or(ENDMARKER)
}

/// Count runs in the endmarker record by grouping consecutive sequences by
/// their first non-endmarker node. Matches the C++ endmarker-run numbering.
fn count_endmarker_runs(gbwt: &GBWT) -> usize {
    let n = gbwt.sequences();
    if n == 0 {
        return 0;
    }
    let mut runs: usize = 1;
    let mut prev = first_node_of_seq(gbwt, 0);
    for i in 1..n {
        let curr = first_node_of_seq(gbwt, i);
        if curr == ENDMARKER || curr != prev {
            runs += 1;
            prev = curr;
        }
    }
    runs
}

/// Output of `lf_with_run`: one LF step through a record plus run info.
#[derive(Debug)]
struct LfInfo {
    /// Next position in the sequence (the result of `record.lf(offset)`).
    next_pos: Pos,
    /// First offset belonging to the run that contained `offset`.
    run_start: usize,
    /// Last offset belonging to that run (inclusive).
    run_end: usize,
    /// Index of that run within the record (0-based).
    run_idx: usize,
}

/// Walks a record's RLE-encoded BWT one step, mirroring `Record::lf` but
/// also returning the containing-run bounds and the **logical** run index.
///
/// When the run's successor is ENDMARKER, each position is treated as its
/// own logical run (matching the C++ `LFLoop` convention), so `run_idx` is
/// the per-position run id and `run_start == run_end == i`.
fn lf_with_run(edges: &[Pos], bwt_bytes: &[u8], sigma: usize, i: usize) -> Option<LfInfo> {
    let mut edges_local = edges.to_vec();
    let mut offset: usize = 0;
    let mut runs_seen: usize = 0;
    for run in RLEIter::with_sigma(bwt_bytes, sigma) {
        // Update local edge offset (used to compute the next position).
        // Advance runs_seen by the logical run count for this physical run.
        let is_endmarker = edges_local[run.value].node == ENDMARKER;
        runs_seen += if is_endmarker { run.len } else { 1 };
        if offset + run.len > i {
            let (run_start, run_end, run_idx) = if is_endmarker {
                // Each endmarker position is its own logical run.
                let local_idx_in_run = i - offset;
                (i, i, runs_seen - (run.len - local_idx_in_run))
            } else {
                (offset, offset + run.len - 1, runs_seen - 1)
            };
            edges_local[run.value].offset += i - offset;
            return Some(LfInfo {
                next_pos: edges_local[run.value],
                run_start,
                run_end,
                run_idx,
            });
        }
        edges_local[run.value].offset += run.len;
        offset += run.len;
    }
    None
}

/// Decode the edge-list prefix of a compressed-record byte slice. This
/// duplicates `gbz::bwt::Record::decompress_edges` since that function is
/// only reachable via owning a `Record` — here we want to decode directly
/// from the byte slice returned by `BWT::compressed_record`.
fn decompress_edges(bytes: &[u8]) -> Option<(Vec<Pos>, usize)> {
    // gbz re-exports Record::decompress_edges as an associated fn.
    gbz::bwt::Record::decompress_edges(bytes)
}

//------------------------------ tests ----------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use gbz::GBWTBuilder;

    /// Build a tiny BIDIRECTIONAL GBWT from a few hard-coded paths and return
    /// `(gbwt, forward_paths)`.
    ///
    /// The C++ FastLocate tests use bidirectional GBWTs and so do we — the
    /// algorithm's predecessor-query invariants depend on the extra reverse
    /// paths guaranteeing that every head sample has a preceding tail sample
    /// in packed text position order.
    fn tiny_gbwt() -> (GBWT, Vec<Vec<usize>>) {
        // Paths of encoded node ids. The builder treats them as arbitrary
        // `usize` sequences; bidirectional=true inserts each path and its
        // reverse. We pick even-numbered encoded ids to mimic the
        // "Node::encode(i, false) = 2*i" convention in the jltsiren tests.
        let paths: Vec<Vec<usize>> = vec![
            vec![2, 4, 6, 8, 10],
            vec![2, 4, 6, 12, 14],
            vec![16, 4, 6, 8, 10],
            vec![2, 4, 18, 20],
        ];
        let mut builder = GBWTBuilder::new(true, false, 1024);
        for p in &paths {
            builder.insert(p, None).expect("insert failed");
        }
        let gbwt = builder.build().expect("build failed");
        (gbwt, paths)
    }

    /// Ground truth via `gbwt.sequence(i)`: for every sequence in the GBWT,
    /// walk it and record (seq_id, forward_node_idx, forward_seq_len) for
    /// each visit. Returns a map `node → Vec<(seq_id, seq_offset_from_end)>`
    /// which is what `decompress_da` should match after sorting.
    fn naive_visits_from_gbwt(gbwt: &GBWT) -> std::collections::HashMap<usize, Vec<(usize, usize)>> {
        let mut out: std::collections::HashMap<usize, Vec<(usize, usize)>> =
            Default::default();
        for seq_id in 0..gbwt.sequences() {
            let Some(iter) = gbwt.sequence(seq_id) else {
                continue;
            };
            let nodes: Vec<usize> = iter.collect();
            let total = nodes.len() + 1; // +1 for endmarker, matches walk seq_offset
            for (fwd_idx, &node) in nodes.iter().enumerate() {
                // forward offset in walk = fwd_idx + 1
                // flipped = total - 1 - (fwd_idx + 1) = total - 2 - fwd_idx
                let seq_off_from_end = total - 2 - fwd_idx;
                out.entry(node).or_default().push((seq_id, seq_off_from_end));
            }
        }
        for v in out.values_mut() {
            v.sort();
        }
        out
    }

    fn assert_fast_locate_matches(gbwt: &GBWT, fl: &FastLocate) {
        let truth = naive_visits_from_gbwt(gbwt);
        for (&node, expected_visits) in &truth {
            if !gbwt.has_node(node) {
                continue;
            }
            let got_raw = fl.decompress_da(gbwt, node);
            let mut got = got_raw.clone();
            got.sort();
            let mut expected_sorted = expected_visits.clone();
            expected_sorted.sort();
            assert_eq!(
                got, expected_sorted,
                "decompress_da mismatch for node {}: got {:?}, expected {:?}",
                node, got, expected_sorted
            );
        }
    }

    #[test]
    fn fast_locate_matches_naive_walk_tiny() {
        let (gbwt, _paths) = tiny_gbwt();
        let fl = FastLocate::build(&gbwt);
        assert!(fl.total_runs() > 0, "should have at least one run");
        assert_fast_locate_matches(&gbwt, &fl);
    }

    #[test]
    fn fast_locate_matches_naive_walk_larger() {
        // Longer, more varied paths with heavy sharing.
        let paths: Vec<Vec<usize>> = vec![
            vec![2, 4, 6, 8, 10, 12, 14, 16],
            vec![2, 4, 6, 8, 10, 18, 20, 22],
            vec![2, 4, 6, 24, 26, 12, 14, 16],
            vec![2, 28, 30, 8, 10, 12, 14, 16],
            vec![32, 4, 6, 8, 10, 12, 14, 34],
            vec![2, 4, 36, 8, 10, 12, 38, 16],
            vec![2, 4, 6, 8, 40],
            vec![42, 44, 6, 8, 10, 12, 46],
        ];
        let mut builder = GBWTBuilder::new(true, false, 4096);
        for p in &paths {
            builder.insert(p, None).expect("insert failed");
        }
        let gbwt = builder.build().expect("build failed");
        let fl = FastLocate::build(&gbwt);
        assert!(fl.total_runs() > 0);
        assert_fast_locate_matches(&gbwt, &fl);
    }

    #[test]
    fn fast_locate_single_path() {
        // Edge case: one path, no sharing, every run is length 1.
        let paths: Vec<Vec<usize>> = vec![vec![2, 4, 6, 8, 10]];
        let mut builder = GBWTBuilder::new(true, false, 1024);
        for p in &paths {
            builder.insert(p, None).expect("insert failed");
        }
        let gbwt = builder.build().expect("build failed");
        let fl = FastLocate::build(&gbwt);
        assert_fast_locate_matches(&gbwt, &fl);
    }

    #[test]
    fn fast_locate_save_load_roundtrip() {
        let (gbwt, _) = tiny_gbwt();
        let fl = FastLocate::build(&gbwt);
        let mut buf: Vec<u8> = Vec::new();
        fl.save(&mut buf).expect("save failed");
        let mut cursor = std::io::Cursor::new(buf);
        let loaded = FastLocate::load(&mut cursor).expect("load failed");
        assert_eq!(loaded.samples, fl.samples);
        assert_eq!(loaded.last_to_run, fl.last_to_run);
        assert_eq!(loaded.comp_to_run, fl.comp_to_run);
        assert_eq!(loaded.total_runs, fl.total_runs);
        assert_eq!(loaded.max_length, fl.max_length);
        // Functional check: decompress_da results match.
        assert_fast_locate_matches(&gbwt, &loaded);
    }

    #[test]
    fn fast_locate_identical_paths() {
        // Edge case: two identical paths. Every run has length 2, and all
        // visits share predecessors — stresses multi-visit run logic.
        let paths: Vec<Vec<usize>> = vec![vec![2, 4, 6, 8], vec![2, 4, 6, 8]];
        let mut builder = GBWTBuilder::new(true, false, 1024);
        for p in &paths {
            builder.insert(p, None).expect("insert failed");
        }
        let gbwt = builder.build().expect("build failed");
        let fl = FastLocate::build(&gbwt);
        assert_fast_locate_matches(&gbwt, &fl);
    }
}
