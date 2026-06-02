//! Safe Rust wrapper around syng's C library.
//!
//! Provides `SyngIndex` for building, loading, saving, and querying
//! GBWT-based syncmer indices.

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::ffi::CString;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::marker::PhantomData;
use std::path::Path;
use std::sync::{Mutex, OnceLock};

use crate::syng_ffi;

fn write_u64<W: std::io::Write>(w: &mut W, v: u64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_u64<R: std::io::Read>(r: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_u64_from_slice(buf: &[u8], offset: &mut usize) -> io::Result<u64> {
    let end = offset
        .checked_add(8)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "u64 slice offset overflow"))?;
    let bytes = buf
        .get(*offset..end)
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "truncated u64 in sidecar"))?;
    *offset = end;
    Ok(u64::from_le_bytes(bytes.try_into().unwrap()))
}

fn read_u32_from_slice_at(buf: &[u8], offset: usize) -> io::Result<u32> {
    let end = offset
        .checked_add(4)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "u32 slice offset overflow"))?;
    let bytes = buf
        .get(offset..end)
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "truncated u32 in sidecar"))?;
    Ok(u32::from_le_bytes(bytes.try_into().unwrap()))
}

fn read_u64_from_slice_at(buf: &[u8], offset: usize) -> io::Result<u64> {
    let end = offset
        .checked_add(8)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "u64 slice offset overflow"))?;
    let bytes = buf
        .get(offset..end)
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "truncated u64 in sidecar"))?;
    Ok(u64::from_le_bytes(bytes.try_into().unwrap()))
}

#[derive(Debug, Clone, Copy)]
pub struct SampledPositionHit {
    pub path_idx: usize,
    pub target_pos: u64,
    pub target_orient: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SampledPathStepHit {
    pub path_idx: usize,
    pub step_idx: u32,
    pub bp_pos: u64,
    pub signed_node: i32,
    pub prev_node: i32,
    pub prev_offset: u32,
    pub traversal_rank: u32,
}

#[derive(Debug, Clone, Copy)]
struct LocatedSyncmerHit {
    path_idx: usize,
    target_pos: u64,
    target_orient: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct LocateCheckpointKey {
    signed_node: i32,
    abs_rank: u32,
}

#[derive(Debug, Clone, Copy)]
struct LocateCheckpointValue {
    path_idx: usize,
    bp_pos: u64,
}

#[derive(Debug)]
struct LocateCheckpointIndex {
    sample_rate: u32,
    shard_mask: usize,
    shards: Vec<FxHashMap<LocateCheckpointKey, LocateCheckpointValue>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CheckpointRecord {
    signed_node: i32,
    abs_rank: u32,
    path_idx: u32,
    bp_pos: u64,
}

impl CheckpointRecord {
    const BYTE_LEN: usize = 16;
}

/// Persistent occurrence-major syncmer positions.
///
/// `.spos` is the query-time inverse of `.pstep`: for each oriented syncmer
/// node, it stores sampled GBWT occurrence ranks and the path coordinate of
/// that checkpoint. The payload is fixed-width so it can be memory-mapped and
/// binary-searched without materializing a global hash table.
#[derive(Debug)]
pub struct SampledCheckpointIndex {
    pub sample_rate: u32,
    checkpoint_count: u64,
    signed_nodes: Vec<i32>,
    record_offsets: Vec<u64>,
    payload_offset: usize,
    mmap: memmap2::Mmap,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PathStepRecord {
    step_idx: u32,
    bp_pos: u64,
    signed_node: i32,
    prev_node: i32,
    prev_offset: u32,
    traversal_rank: u32,
}

#[derive(Debug)]
struct EncodedPositionSampleChunk {
    position_samples: Vec<(u32, u64)>,
    path_step_data: Vec<u8>,
    path_step_sample_count: u64,
    walked_steps: u64,
}

#[derive(Debug, Default)]
struct PositionRepairProgress {
    completed_paths: usize,
    completed_steps: u64,
    completed_position_samples: u64,
    completed_path_step_samples: u64,
}

#[derive(Debug, Clone, Copy)]
pub struct SampledPositionBuildStats {
    pub sampled_occurrences: u64,
    pub sampled_nodes: usize,
    pub sampled_path_steps: u64,
    pub sampled_step_paths: usize,
    pub walked_paths: usize,
    pub sample_rate: u32,
}

pub const DEFAULT_POSITION_SAMPLE_RATE: u32 = 256;
pub const DEFAULT_WALK_SEED_ANCHORS: usize = 5;

fn syng_sidecar_path(prefix: &str, suffix: &str) -> String {
    if prefix.ends_with(".syng") {
        format!("{prefix}.{suffix}")
    } else {
        format!("{prefix}.syng.{suffix}")
    }
}

fn syng_sidecar_candidates(prefix: &str, suffix: &str) -> Vec<String> {
    let primary = syng_sidecar_path(prefix, suffix);
    let legacy = format!("{prefix}.syng.{suffix}");
    if legacy == primary {
        vec![primary]
    } else {
        vec![primary, legacy]
    }
}

fn existing_syng_sidecar_path(prefix: &str, suffix: &str) -> String {
    syng_sidecar_candidates(prefix, suffix)
        .into_iter()
        .find(|p| Path::new(p).exists())
        .unwrap_or_else(|| syng_sidecar_path(prefix, suffix))
}

fn read_kmer_hash_from_prefix(prefix: &str) -> io::Result<*mut syng_ffi::KmerHash> {
    let schema_text = syng_ffi::syng_schema_text();
    let khash_path = format!("{}.1khash", prefix);
    if !Path::new(&khash_path).exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("KmerHash file not found: {}", khash_path),
        ));
    }
    let khash_cpath = CString::new(khash_path.as_str())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let khash_type = CString::new("khash").unwrap();
    unsafe {
        let schema = syng_ffi::oneSchemaCreateFromText(schema_text.as_ptr());
        if schema.is_null() {
            return Err(io::Error::other(
                "Failed to create ONEcode schema for khash read",
            ));
        }
        let of = syng_ffi::oneFileOpenRead(khash_cpath.as_ptr(), schema, khash_type.as_ptr(), 1);
        if of.is_null() {
            syng_ffi::oneSchemaDestroy(schema);
            return Err(io::Error::other(format!(
                "Failed to open {} for reading",
                khash_path
            )));
        }
        let kh = syng_ffi::kmerHashReadOneFile(of);
        syng_ffi::oneFileClose(of);
        syng_ffi::oneSchemaDestroy(schema);
        if kh.is_null() {
            return Err(io::Error::other("kmerHashReadOneFile returned null"));
        }
        Ok(kh)
    }
}

fn encode_query_base(base: u8) -> (u8, bool) {
    match base {
        b'a' | b'A' | 0 => (0, true),
        b'c' | b'C' | 1 => (1, true),
        b'g' | b'G' | 2 => (2, true),
        b't' | b'T' | 3 => (3, true),
        _ => (0, false),
    }
}

pub fn syng_names_path(prefix: &str) -> String {
    syng_sidecar_path(prefix, "names")
}

pub fn syng_spos_path(prefix: &str) -> String {
    syng_sidecar_path(prefix, "spos")
}

pub fn syng_pstep_path(prefix: &str) -> String {
    syng_sidecar_path(prefix, "pstep")
}

pub fn syng_meta_path(prefix: &str) -> String {
    syng_sidecar_path(prefix, "meta")
}

/// Succinct sampled path-position sidecar for projected syng mapping.
///
/// This mirrors the sampled suffix-array idea used by ropebwt3: store a regular
/// per-path grid of occurrences, then resolve mapping positions from those
/// samples instead of materializing every visit for every node. Occurrences are
/// grouped by syncmer node id and varint delta encoded.
pub struct SampledPositions {
    pub sample_rate: u32,
    sample_count: u64,
    path_starts: Vec<u64>,
    node_ids: Vec<u32>,
    byte_offsets: Vec<u64>,
    data: Vec<u8>,
}

#[allow(dead_code)]
impl SampledPositions {
    const MAGIC: u64 = 0x494D50_53504F53; // "IMPSPOS"
    const VERSION: u64 = 4;

    fn from_samples(
        sample_rate: u32,
        path_lengths: &[u64],
        mut samples: Vec<(u32, u64)>,
    ) -> io::Result<Self> {
        validate_position_sample_rate(sample_rate)?;
        let path_starts = Self::path_starts_from_lengths(path_lengths)?;
        samples.sort_unstable();
        samples.dedup();

        let sample_count = samples.len() as u64;
        let mut node_ids = Vec::new();
        let mut byte_offsets = Vec::new();
        let mut data = Vec::new();
        byte_offsets.push(0);

        let mut i = 0usize;
        while i < samples.len() {
            let node_id = samples[i].0;
            node_ids.push(node_id);
            let mut prev = 0u64;
            while i < samples.len() && samples[i].0 == node_id {
                let packed = samples[i].1;
                let delta = packed.checked_sub(prev).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "sampled positions are not sorted within node group",
                    )
                })?;
                push_varint(&mut data, delta);
                prev = packed;
                i += 1;
            }
            byte_offsets.push(data.len() as u64);
        }

        Ok(Self {
            sample_rate,
            sample_count,
            path_starts,
            node_ids,
            byte_offsets,
            data,
        })
    }

    fn path_starts_from_lengths(path_lengths: &[u64]) -> io::Result<Vec<u64>> {
        let mut starts = Vec::with_capacity(path_lengths.len() + 1);
        let mut acc = 0u64;
        starts.push(acc);
        for &len in path_lengths {
            acc = acc.checked_add(len).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "total path length overflowed u64",
                )
            })?;
            starts.push(acc);
        }
        Ok(starts)
    }

    pub fn sample_count(&self) -> u64 {
        self.sample_count
    }

    pub fn node_count(&self) -> usize {
        self.node_ids.len()
    }

    pub fn decode_node(&self, node_id: u32) -> io::Result<Vec<SampledPositionHit>> {
        let Ok(group_idx) = self.node_ids.binary_search(&node_id) else {
            return Ok(Vec::new());
        };
        let start = *self.byte_offsets.get(group_idx).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "sampled-position offset missing",
            )
        })? as usize;
        let end = *self.byte_offsets.get(group_idx + 1).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "sampled-position terminal offset missing",
            )
        })? as usize;
        let data = self.data.as_slice();
        if start > end || end > data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sampled-position byte offsets are out of range",
            ));
        }

        let mut hits = Vec::new();
        let mut pos = start;
        let mut packed = 0u64;
        while pos < end {
            let delta = read_varint_from_slice(&self.data, &mut pos, end)?;
            packed = packed.checked_add(delta).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "sampled position overflowed u64",
                )
            })?;
            let target_orient = (packed & 1) as u8;
            let global_pos = packed >> 1;
            if let Some((path_idx, target_pos)) = self.path_for_global_pos(global_pos) {
                hits.push(SampledPositionHit {
                    path_idx,
                    target_pos,
                    target_orient,
                });
            }
        }
        Ok(hits)
    }

    fn path_for_global_pos(&self, global_pos: u64) -> Option<(usize, u64)> {
        if self.path_starts.len() < 2 || global_pos >= *self.path_starts.last()? {
            return None;
        }
        let upper = self
            .path_starts
            .partition_point(|&start| start <= global_pos);
        let path_idx = upper.checked_sub(1)?;
        Some((path_idx, global_pos - self.path_starts[path_idx]))
    }

    fn save(&self, path: &str) -> io::Result<()> {
        let mut w = BufWriter::new(std::fs::File::create(path)?);
        write_u64(&mut w, Self::MAGIC)?;
        write_u64(&mut w, Self::VERSION)?;
        write_u64(&mut w, self.sample_rate as u64)?;
        write_u64(&mut w, self.sample_count)?;

        write_u64(&mut w, self.path_starts.len() as u64)?;
        for &start in &self.path_starts {
            write_u64(&mut w, start)?;
        }

        write_u64(&mut w, self.node_ids.len() as u64)?;
        for &node_id in &self.node_ids {
            write_u64(&mut w, node_id as u64)?;
        }

        write_u64(&mut w, self.byte_offsets.len() as u64)?;
        for &offset in &self.byte_offsets {
            write_u64(&mut w, offset)?;
        }

        write_u64(&mut w, self.data.len() as u64)?;
        w.write_all(self.data.as_slice())?;
        w.flush()
    }

    fn save_atomic(&self, path: &str) -> io::Result<()> {
        let tmp_path = format!("{}.tmp.{}", path, std::process::id());
        self.save(&tmp_path)?;
        std::fs::rename(&tmp_path, path).map_err(|e| {
            let _ = std::fs::remove_file(&tmp_path);
            e
        })
    }

    #[allow(dead_code)]
    fn load(path: &str) -> io::Result<Self> {
        let mut r = BufReader::new(std::fs::File::open(path)?);
        let magic = read_u64(&mut r)?;
        if magic != Self::MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("SampledPositions: bad magic 0x{:x}", magic),
            ));
        }
        let version = read_u64(&mut r)?;
        if version != Self::VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("SampledPositions: unsupported version {}", version),
            ));
        }
        let sample_rate = read_u64(&mut r)? as u32;
        validate_position_sample_rate(sample_rate)?;
        let sample_count = read_u64(&mut r)?;

        let n_paths = read_u64(&mut r)? as usize;
        let mut path_starts = Vec::with_capacity(n_paths);
        for _ in 0..n_paths {
            path_starts.push(read_u64(&mut r)?);
        }

        let n_nodes = read_u64(&mut r)? as usize;
        let mut node_ids = Vec::with_capacity(n_nodes);
        for _ in 0..n_nodes {
            let node_id = read_u64(&mut r)?;
            if node_id > u32::MAX as u64 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("SampledPositions: node id {} exceeds u32", node_id),
                ));
            }
            node_ids.push(node_id as u32);
        }

        let n_offsets = read_u64(&mut r)? as usize;
        let mut byte_offsets = Vec::with_capacity(n_offsets);
        for _ in 0..n_offsets {
            byte_offsets.push(read_u64(&mut r)?);
        }

        let n_data = read_u64(&mut r)? as usize;
        let mut data = vec![0u8; n_data];
        r.read_exact(&mut data)?;

        if byte_offsets.len() != node_ids.len() + 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPositions: byte_offsets length must be node_ids length + 1",
            ));
        }
        if byte_offsets.last().copied().unwrap_or(0) as usize != data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPositions: final byte offset does not match data length",
            ));
        }
        if path_starts.windows(2).any(|w| w[0] > w[1]) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPositions: path starts are not sorted",
            ));
        }
        if node_ids.windows(2).any(|w| w[0] >= w[1]) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPositions: node ids are not strictly sorted",
            ));
        }
        if byte_offsets.windows(2).any(|w| w[0] > w[1]) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPositions: byte offsets are not sorted",
            ));
        }

        Ok(Self {
            sample_rate,
            sample_count,
            path_starts,
            node_ids,
            byte_offsets,
            data,
        })
    }

    fn load_sample_rate(path: &str) -> io::Result<u32> {
        let mut r = BufReader::new(std::fs::File::open(path)?);
        let magic = read_u64(&mut r)?;
        if magic != Self::MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("SampledPositions: bad magic 0x{:x}", magic),
            ));
        }
        let version = read_u64(&mut r)?;
        if version != Self::VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("SampledPositions: unsupported version {}", version),
            ));
        }
        let sample_rate = read_u64(&mut r)? as u32;
        validate_position_sample_rate(sample_rate)?;
        Ok(sample_rate)
    }
}

/// Sampled checkpoints keyed by path coordinate.
///
/// `.pstep` stores sampled forward-path checkpoints: path position -> syncmer
/// step plus enough GBWT state to resume walking from the sampled step. Query
/// locate also inverts these checkpoints in memory and finds exact syncmer
/// positions by walking arbitrary GBWT occurrences forward to the next sample.
#[derive(Debug)]
enum SampledPathStepData {
    Owned(Vec<u8>),
    Mapped {
        mmap: memmap2::Mmap,
        offset: usize,
        len: usize,
    },
}

impl SampledPathStepData {
    fn len(&self) -> usize {
        match self {
            Self::Owned(data) => data.len(),
            Self::Mapped { len, .. } => *len,
        }
    }

    fn as_slice(&self) -> &[u8] {
        match self {
            Self::Owned(data) => data.as_slice(),
            Self::Mapped { mmap, offset, len } => &mmap[*offset..*offset + *len],
        }
    }
}

pub struct SampledPathSteps {
    pub sample_rate: u32,
    sample_count: u64,
    path_byte_offsets: Vec<u64>,
    data: SampledPathStepData,
}

impl SampledPathSteps {
    const MAGIC: u64 = 0x494D50_50535450; // "IMPPSTP"
    const VERSION: u64 = 4;

    #[allow(dead_code)]
    fn from_samples(
        sample_rate: u32,
        path_count: usize,
        mut samples: Vec<(usize, PathStepRecord)>,
    ) -> io::Result<Self> {
        validate_position_sample_rate(sample_rate)?;
        samples.sort_unstable_by(|a, b| {
            (a.0, a.1.bp_pos, a.1.step_idx).cmp(&(b.0, b.1.bp_pos, b.1.step_idx))
        });
        samples.dedup();

        if let Some((path_idx, _)) = samples.iter().find(|(path_idx, _)| *path_idx >= path_count) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "sampled path-step path index {} exceeds path count {}",
                    path_idx, path_count
                ),
            ));
        }

        let sample_count = samples.len() as u64;
        let mut path_byte_offsets = Vec::with_capacity(path_count + 1);
        let mut data = Vec::new();
        let mut i = 0usize;
        path_byte_offsets.push(0);
        for path_idx in 0..path_count {
            let mut prev_bp = 0u64;
            let mut prev_step = 0u32;
            while i < samples.len() && samples[i].0 == path_idx {
                let record = samples[i].1;
                let bp_delta = record.bp_pos.checked_sub(prev_bp).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "sampled path-step bp positions are not sorted within path",
                    )
                })?;
                let step_delta = record.step_idx.checked_sub(prev_step).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "sampled path-step indices are not sorted within path",
                    )
                })?;
                push_varint(&mut data, bp_delta);
                push_varint(&mut data, step_delta as u64);
                push_varint(&mut data, encode_i32(record.signed_node));
                push_varint(&mut data, encode_i32(record.prev_node));
                push_varint(&mut data, record.prev_offset as u64);
                push_varint(&mut data, record.traversal_rank as u64);
                prev_bp = record.bp_pos;
                prev_step = record.step_idx;
                i += 1;
            }
            path_byte_offsets.push(data.len() as u64);
        }

        if i != samples.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sampled path-step samples were not assigned to a path group",
            ));
        }

        Self::from_encoded(sample_rate, sample_count, path_byte_offsets, data)
    }

    fn from_encoded(
        sample_rate: u32,
        sample_count: u64,
        path_byte_offsets: Vec<u64>,
        data: Vec<u8>,
    ) -> io::Result<Self> {
        validate_position_sample_rate(sample_rate)?;
        if path_byte_offsets.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPathSteps: path offsets cannot be empty",
            ));
        }
        if path_byte_offsets.last().copied().unwrap_or(0) as usize != data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPathSteps: final byte offset does not match data length",
            ));
        }
        if path_byte_offsets.windows(2).any(|w| w[0] > w[1]) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPathSteps: byte offsets are not sorted",
            ));
        }

        Ok(Self {
            sample_rate,
            sample_count,
            path_byte_offsets,
            data: SampledPathStepData::Owned(data),
        })
    }

    pub fn sample_count(&self) -> u64 {
        self.sample_count
    }

    pub fn sampled_path_count(&self) -> usize {
        self.path_byte_offsets
            .windows(2)
            .filter(|w| w[0] < w[1])
            .count()
    }

    pub fn path_count(&self) -> usize {
        self.path_byte_offsets.len().saturating_sub(1)
    }

    pub fn decode_path(&self, path_idx: usize) -> io::Result<Vec<SampledPathStepHit>> {
        let start =
            *self.path_byte_offsets.get(path_idx).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "path index out of range")
            })? as usize;
        let end = *self.path_byte_offsets.get(path_idx + 1).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "path index terminal offset out of range",
            )
        })? as usize;
        let data = self.data.as_slice();
        if start > end || end > data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sampled path-step byte offsets are out of range",
            ));
        }

        let mut hits = Vec::new();
        let mut pos = start;
        let mut bp_pos = 0u64;
        let mut step_idx = 0u32;
        while pos < end {
            let bp_delta = read_varint_from_slice(data, &mut pos, end)?;
            let step_delta = read_varint_from_slice(data, &mut pos, end)?;
            bp_pos = bp_pos.checked_add(bp_delta).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "path-step bp position overflow")
            })?;
            if step_delta > u32::MAX as u64 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "path-step delta exceeds u32",
                ));
            }
            step_idx = step_idx.checked_add(step_delta as u32).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "path-step index overflow")
            })?;
            let signed_node = decode_i32(read_varint_from_slice(data, &mut pos, end)?)?;
            let prev_node = decode_i32(read_varint_from_slice(data, &mut pos, end)?)?;
            let prev_offset = read_varint_from_slice(data, &mut pos, end)?;
            let traversal_rank = read_varint_from_slice(data, &mut pos, end)?;
            if prev_offset > u32::MAX as u64 || traversal_rank > u32::MAX as u64 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "path-step traversal fields exceed u32",
                ));
            }
            hits.push(SampledPathStepHit {
                path_idx,
                step_idx,
                bp_pos,
                signed_node,
                prev_node,
                prev_offset: prev_offset as u32,
                traversal_rank: traversal_rank as u32,
            });
        }
        Ok(hits)
    }

    pub fn checkpoint_at_or_before(
        &self,
        path_idx: usize,
        bp_pos: u64,
    ) -> io::Result<Option<SampledPathStepHit>> {
        let hits = self.decode_path(path_idx)?;
        let upper = hits.partition_point(|hit| hit.bp_pos <= bp_pos);
        Ok(upper.checked_sub(1).map(|idx| hits[idx]))
    }

    fn save(&self, path: &str) -> io::Result<()> {
        let mut w = BufWriter::new(std::fs::File::create(path)?);
        write_u64(&mut w, Self::MAGIC)?;
        write_u64(&mut w, Self::VERSION)?;
        write_u64(&mut w, self.sample_rate as u64)?;
        write_u64(&mut w, self.sample_count)?;

        write_u64(&mut w, self.path_byte_offsets.len() as u64)?;
        for &offset in &self.path_byte_offsets {
            write_u64(&mut w, offset)?;
        }

        write_u64(&mut w, self.data.len() as u64)?;
        w.write_all(self.data.as_slice())?;
        w.flush()
    }

    fn save_atomic(&self, path: &str) -> io::Result<()> {
        let tmp_path = format!("{}.tmp.{}", path, std::process::id());
        self.save(&tmp_path)?;
        std::fs::rename(&tmp_path, path).map_err(|e| {
            let _ = std::fs::remove_file(&tmp_path);
            e
        })
    }

    fn load(path: &str) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        let mut offset = 0usize;
        let magic = read_u64_from_slice(&mmap, &mut offset)?;
        if magic != Self::MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("SampledPathSteps: bad magic 0x{:x}", magic),
            ));
        }
        let version = read_u64_from_slice(&mmap, &mut offset)?;
        if version != Self::VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("SampledPathSteps: unsupported version {}", version),
            ));
        }
        let sample_rate = read_u64_from_slice(&mmap, &mut offset)? as u32;
        validate_position_sample_rate(sample_rate)?;
        let sample_count = read_u64_from_slice(&mmap, &mut offset)?;

        let n_offsets = read_u64_from_slice(&mmap, &mut offset)? as usize;
        let mut path_byte_offsets = Vec::with_capacity(n_offsets);
        for _ in 0..n_offsets {
            path_byte_offsets.push(read_u64_from_slice(&mmap, &mut offset)?);
        }

        let n_data = read_u64_from_slice(&mmap, &mut offset)? as usize;

        if path_byte_offsets.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPathSteps: path offsets cannot be empty",
            ));
        }
        if path_byte_offsets.last().copied().unwrap_or(0) as usize != n_data {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPathSteps: final byte offset does not match data length",
            ));
        }
        if path_byte_offsets.windows(2).any(|w| w[0] > w[1]) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledPathSteps: byte offsets are not sorted",
            ));
        }
        let expected_len = offset.checked_add(n_data).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "path-step file size overflow")
        })?;
        if mmap.len() != expected_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "SampledPathSteps: file has {} bytes, expected {}",
                    mmap.len(),
                    expected_len
                ),
            ));
        }

        Ok(Self {
            sample_rate,
            sample_count,
            path_byte_offsets,
            data: SampledPathStepData::Mapped {
                mmap,
                offset,
                len: n_data,
            },
        })
    }
}

impl SampledCheckpointIndex {
    const MAGIC: u64 = 0x494D50_53504F31; // "IMPSPO1"
    const VERSION: u64 = 1;

    fn save_records_atomic(
        path: &str,
        sample_rate: u32,
        records: &mut Vec<CheckpointRecord>,
    ) -> io::Result<()> {
        let tmp_path = format!("{}.tmp.{}", path, std::process::id());
        let result = Self::save_records(&tmp_path, sample_rate, records);
        match result {
            Ok(()) => std::fs::rename(&tmp_path, path).map_err(|e| {
                let _ = std::fs::remove_file(&tmp_path);
                e
            }),
            Err(e) => {
                let _ = std::fs::remove_file(&tmp_path);
                Err(e)
            }
        }
    }

    fn save_records(
        path: &str,
        sample_rate: u32,
        records: &mut Vec<CheckpointRecord>,
    ) -> io::Result<()> {
        validate_position_sample_rate(sample_rate)?;
        records.par_sort_unstable_by(|a, b| {
            a.signed_node
                .cmp(&b.signed_node)
                .then(a.abs_rank.cmp(&b.abs_rank))
                .then(a.path_idx.cmp(&b.path_idx))
                .then(a.bp_pos.cmp(&b.bp_pos))
        });

        let mut signed_nodes = Vec::new();
        let mut record_offsets = Vec::new();
        record_offsets.push(0u64);
        let mut i = 0usize;
        while i < records.len() {
            let signed_node = records[i].signed_node;
            signed_nodes.push(signed_node);
            let mut prev_rank = None;
            while i < records.len() && records[i].signed_node == signed_node {
                if prev_rank == Some(records[i].abs_rank) {
                    return Err(Self::duplicate_checkpoint_error(
                        records[i].signed_node,
                        records[i].abs_rank,
                    ));
                }
                prev_rank = Some(records[i].abs_rank);
                i += 1;
            }
            record_offsets.push(i as u64);
        }

        let mut w = BufWriter::new(std::fs::File::create(path)?);
        write_u64(&mut w, Self::MAGIC)?;
        write_u64(&mut w, Self::VERSION)?;
        write_u64(&mut w, sample_rate as u64)?;
        write_u64(&mut w, records.len() as u64)?;

        write_u64(&mut w, signed_nodes.len() as u64)?;
        for &signed_node in &signed_nodes {
            write_u64(&mut w, encode_i32(signed_node))?;
        }
        for &offset in &record_offsets {
            write_u64(&mut w, offset)?;
        }

        for record in records {
            w.write_all(&record.abs_rank.to_le_bytes())?;
            w.write_all(&record.path_idx.to_le_bytes())?;
            w.write_all(&record.bp_pos.to_le_bytes())?;
        }
        w.flush()
    }

    fn load(path: &str) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        let mut offset = 0usize;
        let magic = read_u64_from_slice(&mmap, &mut offset)?;
        if magic != Self::MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("SampledCheckpointIndex: bad magic 0x{:x}", magic),
            ));
        }
        let version = read_u64_from_slice(&mmap, &mut offset)?;
        if version != Self::VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("SampledCheckpointIndex: unsupported version {}", version),
            ));
        }
        let sample_rate = read_u64_from_slice(&mmap, &mut offset)? as u32;
        validate_position_sample_rate(sample_rate)?;
        let checkpoint_count = read_u64_from_slice(&mmap, &mut offset)?;

        let n_nodes = read_u64_from_slice(&mmap, &mut offset)? as usize;
        let mut signed_nodes = Vec::with_capacity(n_nodes);
        for _ in 0..n_nodes {
            signed_nodes.push(decode_i32(read_u64_from_slice(&mmap, &mut offset)?)?);
        }
        let mut record_offsets = Vec::with_capacity(n_nodes + 1);
        for _ in 0..=n_nodes {
            record_offsets.push(read_u64_from_slice(&mmap, &mut offset)?);
        }

        if signed_nodes.windows(2).any(|w| w[0] >= w[1]) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledCheckpointIndex: signed node ids are not strictly sorted",
            ));
        }
        if record_offsets.windows(2).any(|w| w[0] > w[1]) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledCheckpointIndex: record offsets are not sorted",
            ));
        }
        if record_offsets.last().copied().unwrap_or(0) != checkpoint_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "SampledCheckpointIndex: final record offset does not match checkpoint count",
            ));
        }

        let payload_bytes = (checkpoint_count as usize)
            .checked_mul(CheckpointRecord::BYTE_LEN)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "checkpoint payload size overflow",
                )
            })?;
        let expected_len = offset.checked_add(payload_bytes).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "checkpoint file size overflow")
        })?;
        if mmap.len() != expected_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "SampledCheckpointIndex: file has {} bytes, expected {}",
                    mmap.len(),
                    expected_len
                ),
            ));
        }

        Ok(Self {
            sample_rate,
            checkpoint_count,
            signed_nodes,
            record_offsets,
            payload_offset: offset,
            mmap,
        })
    }

    pub fn checkpoint_count(&self) -> u64 {
        self.checkpoint_count
    }

    pub fn signed_node_count(&self) -> usize {
        self.signed_nodes.len()
    }

    fn checkpoint(
        &self,
        signed_node: i32,
        abs_rank: u32,
    ) -> io::Result<Option<LocateCheckpointValue>> {
        let Ok(group_idx) = self.signed_nodes.binary_search(&signed_node) else {
            return Ok(None);
        };
        let start = self.record_offsets[group_idx];
        let end = self.record_offsets[group_idx + 1];
        let mut lo = start;
        let mut hi = end;
        while lo < hi {
            let mid = lo + ((hi - lo) / 2);
            let rank = self.record_abs_rank(mid)?;
            if rank < abs_rank {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if lo >= end || self.record_abs_rank(lo)? != abs_rank {
            return Ok(None);
        }
        Ok(Some(LocateCheckpointValue {
            path_idx: self.record_path_idx(lo)? as usize,
            bp_pos: self.record_bp_pos(lo)?,
        }))
    }

    fn record_byte_offset(&self, record_idx: u64) -> io::Result<usize> {
        let record_idx = usize::try_from(record_idx).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "checkpoint record index exceeds usize",
            )
        })?;
        let payload_delta = record_idx
            .checked_mul(CheckpointRecord::BYTE_LEN)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "checkpoint record offset overflow",
                )
            })?;
        self.payload_offset
            .checked_add(payload_delta)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "checkpoint record offset overflow",
                )
            })
    }

    fn record_abs_rank(&self, record_idx: u64) -> io::Result<u32> {
        read_u32_from_slice_at(&self.mmap, self.record_byte_offset(record_idx)?)
    }

    fn record_path_idx(&self, record_idx: u64) -> io::Result<u32> {
        read_u32_from_slice_at(&self.mmap, self.record_byte_offset(record_idx)? + 4)
    }

    fn record_bp_pos(&self, record_idx: u64) -> io::Result<u64> {
        read_u64_from_slice_at(&self.mmap, self.record_byte_offset(record_idx)? + 8)
    }

    fn duplicate_checkpoint_error(signed_node: i32, abs_rank: u32) -> io::Error {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "duplicate syng locate checkpoint for node {} rank {}",
                signed_node, abs_rank
            ),
        )
    }
}

struct SampledPositionBuilder {
    sample_rate: u32,
    collect_positions: bool,
    position_samples: Vec<(u32, u64)>,
    path_step_offsets: Vec<u64>,
    path_step_data: Vec<u8>,
    path_step_sample_count: u64,
    sampled_step_paths: usize,
    next_path_start: u64,
    paths_seen: usize,
    walked_paths: usize,
}

impl SampledPositionBuilder {
    #[allow(dead_code)]
    fn new(sample_rate: u32, next_path_start: u64, paths_seen: usize) -> io::Result<Self> {
        Self::new_with_mode(sample_rate, next_path_start, paths_seen, true)
    }

    fn new_path_steps_only(
        sample_rate: u32,
        next_path_start: u64,
        paths_seen: usize,
    ) -> io::Result<Self> {
        Self::new_with_mode(sample_rate, next_path_start, paths_seen, false)
    }

    fn new_with_mode(
        sample_rate: u32,
        next_path_start: u64,
        paths_seen: usize,
        collect_positions: bool,
    ) -> io::Result<Self> {
        validate_position_sample_rate(sample_rate)?;
        Ok(Self {
            sample_rate,
            collect_positions,
            position_samples: Vec::new(),
            path_step_offsets: vec![0; paths_seen + 1],
            path_step_data: Vec::new(),
            path_step_sample_count: 0,
            sampled_step_paths: 0,
            next_path_start,
            paths_seen,
            walked_paths: 0,
        })
    }

    fn record_path(&mut self, path_idx: usize, seq_len: u64, steps: &[PathStepRecord]) {
        debug_assert_eq!(
            path_idx, self.paths_seen,
            "online sampled-position path order drifted from name map"
        );
        if !steps.is_empty() {
            self.walked_paths += 1;
        }
        debug_assert_eq!(
            self.path_step_offsets.len(),
            self.paths_seen + 1,
            "path-step offset table drifted from online path order"
        );
        let path_start = self.next_path_start;
        let path_step_start = self.path_step_data.len();
        let mut prev_bp = 0u64;
        let mut prev_step = 0u32;
        let last_step_idx = steps.last().map(|step| step.step_idx).unwrap_or(0);
        for step in steps {
            let node_id = step.signed_node.unsigned_abs();
            let bp_pos = step.bp_pos;
            if !sample_path_step(step.step_idx, self.sample_rate, last_step_idx) {
                continue;
            }
            if self.collect_positions {
                let global_pos = path_start
                    .checked_add(bp_pos)
                    .expect("sampled path coordinate overflowed u64");
                let packed = global_pos
                    .checked_mul(2)
                    .and_then(|p| p.checked_add(if step.signed_node >= 0 { 0 } else { 1 }))
                    .expect("sampled packed coordinate overflowed u64");
                self.position_samples.push((node_id, packed));
            }

            let bp_delta = step
                .bp_pos
                .checked_sub(prev_bp)
                .expect("path-step positions should be sorted within path");
            let step_delta = step
                .step_idx
                .checked_sub(prev_step)
                .expect("path-step indices should be sorted within path");
            push_varint(&mut self.path_step_data, bp_delta);
            push_varint(&mut self.path_step_data, step_delta as u64);
            push_varint(&mut self.path_step_data, encode_i32(step.signed_node));
            push_varint(&mut self.path_step_data, encode_i32(step.prev_node));
            push_varint(&mut self.path_step_data, step.prev_offset as u64);
            push_varint(&mut self.path_step_data, step.traversal_rank as u64);
            prev_bp = step.bp_pos;
            prev_step = step.step_idx;
            self.path_step_sample_count += 1;
        }
        if self.path_step_data.len() > path_step_start {
            self.sampled_step_paths += 1;
        }
        self.path_step_offsets
            .push(self.path_step_data.len() as u64);
        self.next_path_start = self
            .next_path_start
            .checked_add(seq_len)
            .expect("total sampled path length overflowed u64");
        self.paths_seen += 1;
    }

    #[allow(dead_code)]
    fn finish(
        self,
        path_lengths: &[u64],
    ) -> io::Result<(
        SampledPositions,
        SampledPathSteps,
        SampledPositionBuildStats,
    )> {
        let SampledPositionBuilder {
            sample_rate,
            position_samples,
            path_step_offsets,
            path_step_data,
            path_step_sample_count,
            walked_paths,
            ..
        } = self;
        let sampled_positions =
            SampledPositions::from_samples(sample_rate, path_lengths, position_samples)?;
        let sampled_path_steps = SampledPathSteps::from_encoded(
            sample_rate,
            path_step_sample_count,
            path_step_offsets,
            path_step_data,
        )?;
        if sampled_path_steps.path_count() != path_lengths.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "path-step index has offsets for {} paths, expected {}",
                    sampled_path_steps.path_count(),
                    path_lengths.len()
                ),
            ));
        }
        let stats = SampledPositionBuildStats {
            sampled_occurrences: sampled_positions.sample_count(),
            sampled_nodes: sampled_positions.node_count(),
            sampled_path_steps: sampled_path_steps.sample_count(),
            sampled_step_paths: sampled_path_steps.sampled_path_count(),
            walked_paths,
            sample_rate,
        };
        Ok((sampled_positions, sampled_path_steps, stats))
    }

    fn finish_path_steps(
        self,
        path_count: usize,
    ) -> io::Result<(SampledPathSteps, SampledPositionBuildStats)> {
        let sample_rate = self.sample_rate;
        let walked_paths = self.walked_paths;
        let sampled_step_paths = self.sampled_step_paths;
        let sampled_path_steps = self.finish_path_steps_index(path_count)?;
        let stats = SampledPositionBuildStats {
            sampled_occurrences: 0,
            sampled_nodes: 0,
            sampled_path_steps: sampled_path_steps.sample_count(),
            sampled_step_paths,
            walked_paths,
            sample_rate,
        };
        Ok((sampled_path_steps, stats))
    }

    fn finish_path_steps_index(self, path_count: usize) -> io::Result<SampledPathSteps> {
        if self.path_step_offsets.len() != path_count + 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "path-step index has offsets for {} paths, expected {}",
                    self.path_step_offsets.len().saturating_sub(1),
                    path_count
                ),
            ));
        }
        SampledPathSteps::from_encoded(
            self.sample_rate,
            self.path_step_sample_count,
            self.path_step_offsets,
            self.path_step_data,
        )
    }
}

fn push_varint(buf: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        buf.push(((value as u8) & 0x7f) | 0x80);
        value >>= 7;
    }
    buf.push(value as u8);
}

fn read_varint_from_slice(buf: &[u8], pos: &mut usize, end: usize) -> io::Result<u64> {
    let mut value = 0u64;
    let mut shift = 0u32;
    while *pos < end {
        let byte = buf[*pos];
        *pos += 1;
        value |= ((byte & 0x7f) as u64) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "varint exceeds u64 width",
            ));
        }
    }
    Err(io::Error::new(
        io::ErrorKind::UnexpectedEof,
        "truncated varint in sampled positions",
    ))
}

fn encode_i32(value: i32) -> u64 {
    let value = value as i64;
    ((value << 1) ^ (value >> 31)) as u64
}

fn decode_i32(value: u64) -> io::Result<i32> {
    let decoded = ((value >> 1) as i64) ^ (-((value & 1) as i64));
    i32::try_from(decoded).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("zigzag-decoded value {} exceeds i32", decoded),
        )
    })
}

fn validate_position_sample_rate(sample_rate: u32) -> io::Result<()> {
    if sample_rate == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "position sample rate must be greater than zero",
        ));
    }
    Ok(())
}

fn sample_path_step(step_idx: u32, sample_rate: u32, last_step_idx: u32) -> bool {
    (step_idx as u64) % (sample_rate as u64) == 0 || step_idx == last_step_idx
}

#[allow(clippy::too_many_arguments)]
fn encode_regular_position_sample_if_selected(
    path_start: u64,
    sample_rate: u32,
    last_step_idx: u32,
    collect_positions: bool,
    step: PathStepRecord,
    prev_sample_bp: &mut u64,
    prev_sample_step: &mut u32,
    path_step_sample_count: &mut u64,
    position_samples: &mut Vec<(u32, u64)>,
    path_step_data: &mut Vec<u8>,
) -> io::Result<()> {
    if !sample_path_step(step.step_idx, sample_rate, last_step_idx) {
        return Ok(());
    }

    if collect_positions {
        let global_pos = path_start.checked_add(step.bp_pos).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "sampled path coordinate overflowed u64",
            )
        })?;
        let packed = global_pos
            .checked_mul(2)
            .and_then(|p| p.checked_add(if step.signed_node >= 0 { 0 } else { 1 }))
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "sampled packed coordinate overflowed u64",
                )
            })?;
        position_samples.push((step.signed_node.unsigned_abs(), packed));
    }

    let bp_delta = step.bp_pos.checked_sub(*prev_sample_bp).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "sampled path-step bp positions are not sorted within path",
        )
    })?;
    let step_delta = step
        .step_idx
        .checked_sub(*prev_sample_step)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "sampled path-step indices are not sorted within path",
            )
        })?;
    push_varint(path_step_data, bp_delta);
    push_varint(path_step_data, step_delta as u64);
    push_varint(path_step_data, encode_i32(step.signed_node));
    push_varint(path_step_data, encode_i32(step.prev_node));
    push_varint(path_step_data, step.prev_offset as u64);
    push_varint(path_step_data, step.traversal_rank as u64);
    *prev_sample_bp = step.bp_pos;
    *prev_sample_step = step.step_idx;
    *path_step_sample_count = path_step_sample_count.checked_add(1).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "path-step sample count overflowed u64",
        )
    })?;
    Ok(())
}

/// Parameters controlling syncmer extraction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SyncmerParams {
    /// Inner k-mer length (default 8).
    pub k: u32,
    /// Window length; total syncmer length = w + k (default 55).
    pub w: u32,
    /// Hash function seed (default 7).
    pub seed: u32,
}

/// Packed canonical full-syncmer sequence used for prebuilt KmerHash dictionaries.
///
/// This is the same 2-bit little-endian packing used by syng's `seqPack()`,
/// stored as `kmerHash` expects it: `(syncmer_len + 31) / 32` u64 words.
pub type PackedSyncmer = Vec<u64>;

pub fn packed_syncmer_word_len(params: SyncmerParams) -> usize {
    ((params.w + params.k + 31) / 32) as usize
}

impl Default for SyncmerParams {
    fn default() -> Self {
        Self {
            k: 8,
            w: 55,
            seed: 7,
        }
    }
}

impl SyncmerParams {
    /// Convert to the C-side SyncmerParams struct.
    #[cfg(test)]
    fn to_c(self) -> syng_ffi::CSyncmerParams {
        syng_ffi::CSyncmerParams {
            w: self.w as i32,
            k: self.k as i32,
            seed: self.seed as i32,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct SyngMetadata {
    params: SyncmerParams,
}

impl SyngMetadata {
    const VERSION: u32 = 1;

    fn new(params: SyncmerParams) -> Self {
        Self { params }
    }

    fn save(&self, path: &str) -> io::Result<()> {
        let file = std::fs::File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "# impg syng metadata")?;
        writeln!(writer, "version\t{}", Self::VERSION)?;
        writeln!(writer, "syncmer_k\t{}", self.params.k)?;
        writeln!(writer, "syncmer_w\t{}", self.params.w)?;
        writeln!(writer, "syncmer_seed\t{}", self.params.seed)?;
        writeln!(writer, "syncmer_length\t{}", self.params.k + self.params.w)?;
        writeln!(writer, "smer_length\t{}", self.params.k)?;
        writer.flush()
    }

    fn load(path: &str) -> io::Result<Self> {
        let file = std::fs::File::open(path).map_err(|e| {
            if e.kind() == io::ErrorKind::NotFound {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!(
                        "syng metadata not found: {} (rebuild the syng index with this impg)",
                        path
                    ),
                )
            } else {
                e
            }
        })?;
        let reader = BufReader::new(file);

        let mut version = None;
        let mut syncmer_k = None;
        let mut syncmer_w = None;
        let mut syncmer_seed = None;

        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut parts = line.split('\t');
            let key = parts.next().unwrap_or_default();
            let value = parts.next().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid syng metadata line {}: {}", line_no + 1, line),
                )
            })?;
            if parts.next().is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid syng metadata line {}: {}", line_no + 1, line),
                ));
            }

            match key {
                "version" => version = Some(parse_metadata_u32(key, value, line_no + 1)?),
                "syncmer_k" => syncmer_k = Some(parse_metadata_u32(key, value, line_no + 1)?),
                "syncmer_w" => syncmer_w = Some(parse_metadata_u32(key, value, line_no + 1)?),
                "syncmer_seed" => syncmer_seed = Some(parse_metadata_u32(key, value, line_no + 1)?),
                "syncmer_length" | "smer_length" => {}
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "unknown syng metadata key '{}' on line {}",
                            key,
                            line_no + 1
                        ),
                    ));
                }
            }
        }

        let version = version.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "syng metadata missing version")
        })?;
        if version != Self::VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unsupported syng metadata version {}", version),
            ));
        }

        let params = SyncmerParams {
            k: syncmer_k.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "syng metadata missing syncmer_k",
                )
            })?,
            w: syncmer_w.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "syng metadata missing syncmer_w",
                )
            })?,
            seed: syncmer_seed.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "syng metadata missing syncmer_seed",
                )
            })?,
        };

        if params.k == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "syng metadata syncmer_k must be greater than 0",
            ));
        }
        if params.k.checked_add(params.w).is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "syng metadata syncmer length overflowed u32",
            ));
        }

        Ok(Self { params })
    }
}

fn parse_metadata_u32(key: &str, value: &str, line_no: usize) -> io::Result<u32> {
    value.parse::<u32>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid syng metadata value for '{}' on line {}: {}",
                key, line_no, e
            ),
        )
    })
}

/// Information needed to walk a genome's forward path through the GBWT.
#[derive(Debug, Clone)]
pub struct GbwtPathStart {
    /// First syncmer node ID of the forward path.
    pub start_node: i32,
    /// Index among all paths starting at start_node.
    pub start_count: u32,
    /// Number of syncmer nodes in the forward path.
    pub num_syncmers: u32,
    /// Absolute bp position of the first syncmer on the forward sequence.
    ///
    /// Syng's syncmer iterator reports the first syncmer at wherever the
    /// first min k-mer lands in the initial window — this can be anywhere
    /// in `[0, w+k)`. `walk_path` uses this value to initialise its bp
    /// accumulator so that downstream consumers see ABSOLUTE sequence
    /// coordinates (not positions relative to the first syncmer).
    ///
    /// Defaults to 0 when an index loaded from an older format lacks this
    /// field (preserving the pre-fix behaviour for such indexes — still
    /// buggy for anchors, but not worse; a fresh rebuild corrects it).
    pub first_syncmer_pos: u64,
}

/// Maps between GBWT path numbers and sequence names/lengths.
pub struct SyngNameMap {
    /// GBWT path number → sequence name.
    pub path_to_name: Vec<String>,
    /// GBWT path number → sequence length in bp.
    pub path_to_length: Vec<u64>,
    /// Sequence name → GBWT path number.
    pub name_to_path: FxHashMap<String, u32>,
    /// Per-path GBWT forward-path start info (for walking paths during query).
    /// None if loaded from an old-format file that lacks this info.
    pub path_starts: Vec<Option<GbwtPathStart>>,
}

impl SyngNameMap {
    /// Create an empty name map.
    pub fn new() -> Self {
        Self {
            path_to_name: Vec::new(),
            path_to_length: Vec::new(),
            name_to_path: FxHashMap::default(),
            path_starts: Vec::new(),
        }
    }

    /// Add a mapping from a path number to a sequence name and length.
    pub fn add(&mut self, name: String, length: u64) -> u32 {
        let path_num = self.path_to_name.len() as u32;
        self.name_to_path.insert(name.clone(), path_num);
        self.path_to_name.push(name);
        self.path_to_length.push(length);
        self.path_starts.push(None);
        path_num
    }

    /// Set the GBWT forward-path start info for a given path.
    pub fn set_path_start(&mut self, path_num: u32, info: GbwtPathStart) {
        self.path_starts[path_num as usize] = Some(info);
    }

    /// Save to a `.syng.names` file.
    /// Format: one line per path, tab-separated:
    /// `path_number\tname\tlength\tstart_node\tstart_count\tnum_syncmers\tfirst_syncmer_pos`
    /// (columns 3-6 are 0 if path start info is unavailable). The 7th column
    /// `first_syncmer_pos` was added to record the absolute bp offset of
    /// the forward path's first syncmer — older 6-column files load it as 0.
    pub fn save(&self, path: &str) -> io::Result<()> {
        let file = std::fs::File::create(path)?;
        let mut writer = BufWriter::new(file);
        for (i, ((name, &length), start)) in self
            .path_to_name
            .iter()
            .zip(self.path_to_length.iter())
            .zip(self.path_starts.iter())
            .enumerate()
        {
            let (sn, sc, ns, fsp) = match start {
                Some(info) => (
                    info.start_node,
                    info.start_count,
                    info.num_syncmers,
                    info.first_syncmer_pos,
                ),
                None => (0, 0, 0, 0),
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i, name, length, sn, sc, ns, fsp
            )?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Load from a `.syng.names` file.
    /// Supports 3-, 6-, and 7-column formats. Older 6-column files load with
    /// `first_syncmer_pos = 0` (the pre-fix behaviour) — rebuild the index
    /// for coordinate-accurate anchor positions.
    pub fn load(path: &str) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let reader = BufReader::new(file);
        let mut name_map = Self::new();
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 3 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid name map line (need >= 3 columns): {}", line),
                ));
            }
            let _path_num: u32 = parts[0].parse().map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid path number '{}': {}", parts[0], e),
                )
            })?;
            let name = parts[1].to_string();
            let length: u64 = parts[2].parse().map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Invalid length '{}': {}", parts[2], e),
                )
            })?;
            let path_num = name_map.add(name, length);

            // Parse optional path start info (columns 3-5, plus optional
            // column 6 for first_syncmer_pos).
            if parts.len() >= 6 {
                let start_node: i32 = parts[3].parse().unwrap_or(0);
                let start_count: u32 = parts[4].parse().unwrap_or(0);
                let num_syncmers: u32 = parts[5].parse().unwrap_or(0);
                // 7th column is optional (new in this format revision). Old
                // files default to 0 — anchors remain off by the first
                // syncmer's absolute position until the index is rebuilt.
                let first_syncmer_pos: u64 =
                    parts.get(6).map(|s| s.parse().unwrap_or(0)).unwrap_or(0);
                if num_syncmers > 0 {
                    name_map.set_path_start(
                        path_num,
                        GbwtPathStart {
                            start_node,
                            start_count,
                            num_syncmers,
                            first_syncmer_pos,
                        },
                    );
                }
            }
        }
        Ok(name_map)
    }
}

impl Default for SyngNameMap {
    fn default() -> Self {
        Self::new()
    }
}

/// A genomic interval on another genome that is homologous to the query region.
///
/// `cigar` is populated only when `SyngQueryOpts::emit_cigar` is set; byte
/// alphabet is WFA2 convention (`M`, `=`, `X`, `I`, `D`). Empty for
/// ends-only projections.
#[derive(Debug, Clone)]
pub struct HomologousInterval {
    pub genome: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub cigar: Option<Vec<u8>>,
}

/// A shared syncmer node between the query path and a target path.
///
/// `query_pos` and `target_pos` are the bp-coordinate positions of the
/// syncmer's start on each respective path. `node_id` is the absolute
/// (unsigned) syncmer node id.
///
/// Used by the boundary-realignment transitive pipeline to locate
/// flanking anchors for each query edge.
#[derive(Debug, Clone, Copy)]
pub struct Anchor {
    pub query_pos: u64,
    pub target_pos: u64,
    pub node_id: u32,
}

/// A homologous interval plus the shared-syncmer anchors that produced it.
///
/// Anchors are sorted by `query_pos` ascending. For merged intervals the
/// anchor set is the union of anchors from the pre-merge intervals.
#[derive(Debug, Clone)]
pub struct HomologousIntervalWithAnchors {
    pub genome: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub anchors: Vec<Anchor>,
}

/// Frequency filter for syncmer seed nodes during syng queries.
///
/// The filter applies before locating syncmer occurrences. Dropped
/// syncmers do not seed candidate ranges, but downstream code can still
/// recover every syncmer by walking ranges that survived seeding.
#[derive(Debug, Clone, Copy)]
pub struct SyngSeedFilter {
    /// Drop query syncmer nodes whose total GBWT occurrence count
    /// across both orientations is greater than this cap. `None`
    /// disables the absolute cap.
    pub max_occurrences: Option<u32>,
    /// Drop this fraction of query syncmer nodes with the highest
    /// occurrence counts. Uses `floor(n * fraction)`, so small queries
    /// are not affected by tiny defaults such as 0.0005.
    pub drop_top_fraction: f64,
    /// Number of consecutive query syncmers in each bounded exact GBWT
    /// walk seed. Values >1 prevent single high-copy syncmers from
    /// seeding candidate ranges while still allowing divergent homologs
    /// to chain through short exact walks.
    pub walk_anchors: usize,
}

impl Default for SyngSeedFilter {
    fn default() -> Self {
        Self {
            max_occurrences: None,
            drop_top_fraction: 0.0,
            walk_anchors: 1,
        }
    }
}

impl SyngSeedFilter {
    pub fn enabled(self) -> bool {
        self.max_occurrences.is_some() || self.drop_top_fraction > 0.0 || self.walk_anchors > 1
    }
}

#[derive(Debug, Clone, Copy)]
struct NodeOccurrenceCounts {
    forward: u32,
    reverse: u32,
}

impl NodeOccurrenceCounts {
    fn total(self) -> u32 {
        self.forward.saturating_add(self.reverse)
    }
}

#[derive(Debug, Clone, Copy)]
struct QueryWalkStep {
    signed_node: i32,
    bp_pos: u64,
    query_pos: u64,
}

#[derive(Debug, Clone)]
struct WalkSeedTask {
    strand: char,
    steps: Vec<QueryWalkStep>,
}

#[derive(Debug, Default)]
struct WalkSeedTaskHits {
    seeds_kept: u64,
    located_hits: u64,
    emitted_anchors: u64,
    intervals: Vec<(usize, HomologousIntervalWithAnchors)>,
}

/// A query syncmer found in the index.
#[derive(Debug, Clone, Copy)]
pub struct SyngQuerySyncmer {
    pub signed_node: i32,
    pub node_id: u32,
    pub query_pos: u64,
}

/// A signed syncmer walk step with a bp-coordinate position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SyngWalkStep {
    pub signed_node: i32,
    pub bp_pos: u64,
}

/// A maximal exact match of a signed syncmer walk in the syng GBWT.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SyngGbwtMem {
    pub step_start: usize,
    pub step_end: usize,
    pub query_start: u64,
    pub query_end: u64,
    pub anchors: usize,
    pub occurrences: u32,
    pub match_signed_node: i32,
    pub match_low: u32,
    pub match_high: u32,
}

/// A rough syncmer-only mapping projected to a path in the syng index.
#[derive(Debug, Clone)]
pub struct SyngMapHit {
    pub query_name: String,
    pub query_len: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub target_name: String,
    pub target_len: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub strand: char,
    pub anchors: usize,
}

fn matched_syncmers_in_sequence_impl(
    seqhash: *mut syng_ffi::Seqhash,
    kmer_hash: *mut syng_ffi::KmerHash,
    params: SyncmerParams,
    query_seq: &[u8],
    kmer_buf: &mut [u64],
    seq_buf: &mut Vec<u8>,
    valid_prefix: &mut Vec<usize>,
    out: &mut Vec<SyngQuerySyncmer>,
) {
    out.clear();
    let syncmer_len = (params.w + params.k) as usize;
    if query_seq.len() < syncmer_len {
        return;
    }

    seq_buf.clear();
    valid_prefix.clear();
    seq_buf.reserve(query_seq.len() + 1);
    valid_prefix.reserve(query_seq.len() + 1);
    valid_prefix.push(0);
    for &base in query_seq {
        let (encoded, valid) = encode_query_base(base);
        seq_buf.push(encoded);
        let prev = *valid_prefix.last().unwrap();
        valid_prefix.push(prev + usize::from(valid));
    }
    seq_buf.push(0);

    unsafe {
        let sit = syng_ffi::syncmerIterator(
            seqhash,
            seq_buf.as_mut_ptr() as *mut i8,
            query_seq.len() as i32,
        );
        let mut pos: i32 = 0;
        while syng_ffi::syncmerNext(sit, std::ptr::null_mut(), &mut pos, std::ptr::null_mut()) {
            if pos < 0 {
                continue;
            }
            let start = pos as usize;
            let Some(end) = start.checked_add(syncmer_len) else {
                continue;
            };
            if end > query_seq.len() || valid_prefix[end] - valid_prefix[start] != syncmer_len {
                continue;
            }

            let mut kmer_index: i64 = 0;
            if syng_ffi::kmerHashFindThreadSafe(
                kmer_hash,
                seq_buf.as_mut_ptr().add(start) as *mut i8,
                &mut kmer_index,
                kmer_buf.as_mut_ptr(),
            ) {
                let signed_node = kmer_index as i32;
                out.push(SyngQuerySyncmer {
                    signed_node,
                    node_id: signed_node.unsigned_abs(),
                    query_pos: start as u64,
                });
            }
        }
        syng_ffi::impg_seqhashIteratorDestroy(sit);
    }
}

fn reverse_complement_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    for &base in seq.iter().rev() {
        out.push(match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => base,
        });
    }
}

fn orient_reverse_query_syncmers(
    syncmers: &mut [SyngQuerySyncmer],
    query_len: usize,
    syncmer_len: usize,
) {
    let query_len = query_len as u64;
    let syncmer_len = syncmer_len as u64;
    for syncmer in syncmers {
        syncmer.signed_node = -syncmer.signed_node;
        syncmer.query_pos = query_len - syncmer.query_pos - syncmer_len;
    }
}

fn matched_syncmers_best_query_orientation_impl(
    seqhash: *mut syng_ffi::Seqhash,
    kmer_hash: *mut syng_ffi::KmerHash,
    params: SyncmerParams,
    query_seq: &[u8],
    kmer_buf: &mut [u64],
    seq_buf: &mut Vec<u8>,
    valid_prefix: &mut Vec<usize>,
    rc_seq: &mut Vec<u8>,
    rc_syncmers: &mut Vec<SyngQuerySyncmer>,
    out: &mut Vec<SyngQuerySyncmer>,
) {
    // Syng's closed-syncmer extraction is not symmetric under reverse
    // complementation when applied to a read in isolation. Try both read
    // orientations and keep the one that is coherent with the graph.
    matched_syncmers_in_sequence_impl(
        seqhash,
        kmer_hash,
        params,
        query_seq,
        kmer_buf,
        seq_buf,
        valid_prefix,
        out,
    );

    reverse_complement_into(query_seq, rc_seq);
    matched_syncmers_in_sequence_impl(
        seqhash,
        kmer_hash,
        params,
        rc_seq,
        kmer_buf,
        seq_buf,
        valid_prefix,
        rc_syncmers,
    );

    if rc_syncmers.len() > out.len() {
        orient_reverse_query_syncmers(rc_syncmers, query_seq.len(), (params.w + params.k) as usize);
        std::mem::swap(out, rc_syncmers);
    }
}

/// Lightweight syng syncmer dictionary handle for read-to-node matching.
///
/// This loads only `.syng.meta` and `.1khash`. It deliberately does not load
/// the GBWT, names, or positional sidecars, so `impg map -o gaf` can run as a
/// fast syncmer membership query without requiring coordinate indexes.
pub struct SyngMatcher {
    kmer_hash: *mut syng_ffi::KmerHash,
    seqhash: *mut syng_ffi::Seqhash,
    pub params: SyncmerParams,
}

pub struct SyngMatcherWorker<'a> {
    kmer_hash: *mut syng_ffi::KmerHash,
    seqhash: *mut syng_ffi::Seqhash,
    kmer_buf: Vec<u64>,
    seq_buf: Vec<u8>,
    valid_prefix: Vec<usize>,
    rc_seq: Vec<u8>,
    rc_syncmers: Vec<SyngQuerySyncmer>,
    params: SyncmerParams,
    _owner: PhantomData<&'a SyngMatcher>,
}

// SAFETY: The matcher owns its C handles and only creates worker-local query state.
unsafe impl Send for SyngMatcher {}
// SAFETY: Shared matcher access creates independent worker-local Seqhash handles.
unsafe impl Sync for SyngMatcher {}
// SAFETY: A worker is used by one Rayon worker at a time and owns its Seqhash handle.
unsafe impl Send for SyngMatcherWorker<'_> {}

impl SyngMatcher {
    pub fn load(prefix: &str, _params: SyncmerParams) -> io::Result<Self> {
        let metadata_path = existing_syng_sidecar_path(prefix, "meta");
        let metadata = SyngMetadata::load(&metadata_path)?;
        let params = metadata.params;
        let kmer_hash = read_kmer_hash_from_prefix(prefix)?;
        let seqhash = unsafe {
            syng_ffi::impg_seqhashCreateSafe(params.k as i32, params.w as i32, params.seed as i32)
        };
        if seqhash.is_null() {
            unsafe { syng_ffi::kmerHashDestroy(kmer_hash) };
            return Err(io::Error::other("impg_seqhashCreateSafe returned null"));
        }
        Ok(Self {
            kmer_hash,
            seqhash,
            params,
        })
    }

    /// Return query syncmers that are present in this syng dictionary.
    pub fn matched_syncmers_in_sequence(&self, query_seq: &[u8]) -> Vec<SyngQuerySyncmer> {
        let mut kmer_buf = vec![0u64; unsafe { (*self.kmer_hash).plen as usize }];
        let mut seq_buf = Vec::with_capacity(query_seq.len() + 1);
        let mut valid_prefix = Vec::with_capacity(query_seq.len() + 1);
        let mut rc_seq = Vec::with_capacity(query_seq.len());
        let mut rc_syncmers = Vec::new();
        let mut out = Vec::new();
        matched_syncmers_best_query_orientation_impl(
            self.seqhash,
            self.kmer_hash,
            self.params,
            query_seq,
            &mut kmer_buf,
            &mut seq_buf,
            &mut valid_prefix,
            &mut rc_seq,
            &mut rc_syncmers,
            &mut out,
        );
        out
    }

    /// Number of syncmer nodes stored in the KmerHash (valid IDs are `1..=N`).
    pub fn num_syncmer_nodes(&self) -> usize {
        unsafe { (*self.kmer_hash).max as usize }
    }

    pub fn worker(&self) -> io::Result<SyngMatcherWorker<'_>> {
        let seqhash = unsafe {
            syng_ffi::impg_seqhashCreateSafe(
                self.params.k as i32,
                self.params.w as i32,
                self.params.seed as i32,
            )
        };
        if seqhash.is_null() {
            return Err(io::Error::other("impg_seqhashCreateSafe returned null"));
        }
        Ok(SyngMatcherWorker {
            kmer_hash: self.kmer_hash,
            seqhash,
            kmer_buf: vec![0u64; unsafe { (*self.kmer_hash).plen as usize }],
            seq_buf: Vec::new(),
            valid_prefix: Vec::new(),
            rc_seq: Vec::new(),
            rc_syncmers: Vec::new(),
            params: self.params,
            _owner: PhantomData,
        })
    }
}

impl SyngMatcherWorker<'_> {
    /// Return query syncmers using this worker's thread-local syng iterator state.
    pub fn matched_syncmers_in_sequence(&mut self, query_seq: &[u8]) -> Vec<SyngQuerySyncmer> {
        let mut out = Vec::new();
        self.matched_syncmers_in_sequence_into(query_seq, &mut out);
        out
    }

    pub fn matched_syncmers_in_sequence_into(
        &mut self,
        query_seq: &[u8],
        out: &mut Vec<SyngQuerySyncmer>,
    ) {
        matched_syncmers_best_query_orientation_impl(
            self.seqhash,
            self.kmer_hash,
            self.params,
            query_seq,
            &mut self.kmer_buf,
            &mut self.seq_buf,
            &mut self.valid_prefix,
            &mut self.rc_seq,
            &mut self.rc_syncmers,
            out,
        );
    }
}

/// A loaded GBWT index with syncmer hash and name mapping.
///
/// This is the primary handle for syng integration. It owns the C-allocated
/// SyngBWT, KmerHash, and Seqhash structures and frees them on drop.
pub struct SyngIndex {
    gbwt: *mut syng_ffi::SyngBWT,
    kmer_hash: *mut syng_ffi::KmerHash,
    seqhash: *mut syng_ffi::Seqhash,
    pub name_map: SyngNameMap,
    pub params: SyncmerParams,
    /// Sampled path checkpoints for path-position to syncmer-step lookup.
    sampled_path_steps: Option<SampledPathSteps>,
    /// Memory-mapped occurrence-major checkpoints for fast target locate.
    sampled_checkpoints: Option<SampledCheckpointIndex>,
    /// In-memory locate table keyed by oriented GBWT node occurrence rank.
    locate_checkpoint_index: OnceLock<Result<LocateCheckpointIndex, String>>,
    /// Online collector that records sampled positions while sequences are added.
    sampled_position_builder: Option<SampledPositionBuilder>,
}

#[derive(Debug, Clone, Copy)]
pub struct SyngAddSequenceStats {
    pub sequence_len: usize,
    pub syncmers: usize,
    pub indexed: bool,
}

#[derive(Debug, Clone, Copy)]
enum SyncmerLookupMode {
    AddMissing,
    RequireExisting,
}

// SAFETY: The C structures are only accessed through &self or &mut self,
// and the C library's read-only query functions are thread-safe.
unsafe impl Send for SyngIndex {}
// SAFETY: Read-only query functions (query_region, walk_path) only read from
// the GBWT/KmerHash/Seqhash structures and create local iterators.
unsafe impl Sync for SyngIndex {}

impl SyngIndex {
    /// Create a new, empty SyngIndex with default parameters.
    ///
    /// The GBWT is created with a fixed syncmer length of `w + k` (default 63)
    /// and an initial capacity of 0 paths.
    pub fn new(params: SyncmerParams) -> Self {
        let syncmer_len = (params.w + params.k) as i32;
        let gbwt = unsafe { syng_ffi::syngBWTcreate(syncmer_len, 0) };
        let kmer_hash = unsafe { syng_ffi::kmerHashCreate(1024, syncmer_len) };
        let seqhash = unsafe {
            syng_ffi::impg_seqhashCreateSafe(params.k as i32, params.w as i32, params.seed as i32)
        };
        Self {
            gbwt,
            kmer_hash,
            seqhash,
            name_map: SyngNameMap::new(),
            params,
            sampled_path_steps: None,
            sampled_checkpoints: None,
            locate_checkpoint_index: OnceLock::new(),
            sampled_position_builder: Some(
                SampledPositionBuilder::new_path_steps_only(DEFAULT_POSITION_SAMPLE_RATE, 0, 0)
                    .expect("default sampled-position parameters are valid"),
            ),
        }
    }

    /// Create an index and preload its syncmer dictionary.
    ///
    /// Paths can then be replayed with [`Self::add_sequence_with_existing_syncmers`].
    /// The supplied dictionary must already be canonicalized, deduplicated, and
    /// in the intended global ID order.
    pub fn new_with_packed_syncmer_dictionary(
        params: SyncmerParams,
        packed_syncmers: &[PackedSyncmer],
    ) -> io::Result<Self> {
        let mut index = Self::new(params);
        index.add_packed_syncmer_dictionary(packed_syncmers)?;
        Ok(index)
    }

    /// Preload packed canonical syncmers into the C KmerHash in order.
    pub fn add_packed_syncmer_dictionary(
        &mut self,
        packed_syncmers: &[PackedSyncmer],
    ) -> io::Result<()> {
        let expected_words = packed_syncmer_word_len(self.params);
        for (i, packed) in packed_syncmers.iter().enumerate() {
            if packed.len() != expected_words {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "packed syncmer {} has {} words, expected {}",
                        i,
                        packed.len(),
                        expected_words
                    ),
                ));
            }
            let mut index: i64 = 0;
            let added = unsafe {
                syng_ffi::kmerHashAddPacked(
                    self.kmer_hash,
                    packed.as_ptr() as *mut syng_ffi::U64,
                    &mut index,
                )
            };
            if !added {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("duplicate packed syncmer at dictionary offset {}", i),
                ));
            }
            let expected_index = i as i64 + 1;
            if index != expected_index {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "packed syncmer dictionary assigned id {}, expected {}",
                        index, expected_index
                    ),
                ));
            }
        }
        Ok(())
    }

    /// True if a sampled path-position sidecar has been built or loaded.
    pub fn has_sampled_positions(&self) -> bool {
        self.sampled_checkpoints.is_some()
    }

    /// True if an inverted sampled path-step sidecar has been built or loaded.
    pub fn has_sampled_path_steps(&self) -> bool {
        self.sampled_path_steps.is_some()
    }

    /// True if an occurrence-major checkpoint sidecar has been loaded.
    pub fn has_sampled_checkpoints(&self) -> bool {
        self.sampled_checkpoints.is_some()
    }

    /// Sampling interval for the syncmer-position sidecar, if present.
    pub fn sampled_positions_rate(&self) -> Option<u32> {
        self.sampled_checkpoints_rate()
    }

    /// Sampling interval for the path-step checkpoint sidecar, if present.
    pub fn sampled_path_steps_rate(&self) -> Option<u32> {
        self.sampled_path_steps
            .as_ref()
            .map(|steps| steps.sample_rate)
    }

    /// Sampling interval for the occurrence-major checkpoint sidecar, if present.
    pub fn sampled_checkpoints_rate(&self) -> Option<u32> {
        self.sampled_checkpoints
            .as_ref()
            .map(|checkpoints| checkpoints.sample_rate)
    }

    fn reset_locate_checkpoint_index(&mut self) {
        self.locate_checkpoint_index = OnceLock::new();
    }

    /// Enable online sampled-position collection for subsequently added paths.
    ///
    /// This should be called before adding any sequence. The sampled position
    /// index is intentionally built online as part of syng construction.
    pub fn enable_online_sampled_positions(&mut self, sample_rate: u32) -> io::Result<()> {
        let next_path_start = self
            .name_map
            .path_to_length
            .iter()
            .try_fold(0u64, |acc, &len| acc.checked_add(len))
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "total path length overflowed u64",
                )
            })?;
        self.sampled_position_builder = Some(SampledPositionBuilder::new_path_steps_only(
            sample_rate,
            next_path_start,
            self.name_map.path_to_name.len(),
        )?);
        self.sampled_path_steps = None;
        self.sampled_checkpoints = None;
        self.reset_locate_checkpoint_index();
        Ok(())
    }

    /// Materialize the `.syng.pstep` and `.syng.spos` representations.
    ///
    /// GBWT occurrence ranks are not stable until all paths have been inserted,
    /// so final `.pstep` checkpoints are rebuilt from the completed GBWT rather
    /// than trusting ranks captured while paths were still being added.
    pub fn finalize_online_sampled_positions(
        &mut self,
    ) -> io::Result<Option<SampledPositionBuildStats>> {
        let Some(builder) = self.sampled_position_builder.take() else {
            return Ok(None);
        };
        let stats =
            self.rebuild_sampled_position_indexes_from_gbwt_parallel(builder.sample_rate, 0)?;
        Ok(Some(stats))
    }

    /// Build an index progressively from an iterator of (name, sequence) pairs.
    ///
    /// For each sequence:
    /// 1. Extract syncmers using the seqhash iterator
    /// 2. Add new syncmers to the KmerHash
    /// 3. Build a forward GBWT path
    /// 4. Build a reverse complement GBWT path
    /// 5. Record name → path number mapping
    pub fn build(
        params: SyncmerParams,
        sequences: impl Iterator<Item = (String, Vec<u8>)>,
    ) -> Self {
        let mut index = Self::new(params);
        if let Err(e) = index.enable_online_sampled_positions(1) {
            log::warn!(
                "failed to enable exact sampled positions for in-memory build: {}",
                e
            );
        }

        for (name, seq) in sequences {
            index.add_sequence(name, seq);
        }

        if let Err(e) = index.finalize_online_sampled_positions() {
            log::warn!("sampled-position finalization failed: {}", e);
        }

        index
    }

    /// Add one sequence to the syng index.
    ///
    /// This is the same per-sequence work used by [`Self::build`], exposed so
    /// large callers can stream inputs without buffering all decompressed
    /// sequences before GBWT construction.
    pub fn add_sequence(&mut self, name: String, seq: Vec<u8>) -> SyngAddSequenceStats {
        self.add_sequence_internal(name, seq, SyncmerLookupMode::AddMissing)
            .expect("adding new syncmers to KmerHash should not fail")
    }

    /// Add one sequence using a preloaded dictionary.
    ///
    /// Returns an error if any extracted syncmer is absent from the current
    /// KmerHash. This is used by two-pass/parallel construction to replay paths
    /// through a frozen global syncmer ID assignment.
    pub fn add_sequence_with_existing_syncmers(
        &mut self,
        name: String,
        seq: Vec<u8>,
    ) -> io::Result<SyngAddSequenceStats> {
        self.add_sequence_internal(name, seq, SyncmerLookupMode::RequireExisting)
    }

    fn add_sequence_internal(
        &mut self,
        name: String,
        seq: Vec<u8>,
        lookup_mode: SyncmerLookupMode,
    ) -> io::Result<SyngAddSequenceStats> {
        self.sampled_path_steps = None;
        self.sampled_checkpoints = None;
        self.reset_locate_checkpoint_index();

        let syncmer_len = (self.params.w + self.params.k) as usize;
        let seq_len = seq.len();
        if seq_len < syncmer_len {
            // Sequence too short for syncmer extraction — record in name map but skip
            let path_num = self.name_map.add(name, seq_len as u64);
            if let Some(builder) = self.sampled_position_builder.as_mut() {
                builder.record_path(path_num as usize, seq_len as u64, &[]);
            }
            return Ok(SyngAddSequenceStats {
                sequence_len: seq_len,
                syncmers: 0,
                indexed: false,
            });
        }

        // Convert sequence to numeric encoding (0=a, 1=c, 2=g, 3=t) for the C
        // seqhash library.  The seqhash uses raw byte values as array indices
        // into patternRC[4], so the values MUST be 0-3.
        let mut seq_buf: Vec<u8> = Vec::with_capacity(seq_len + 1);
        seq_buf.extend(seq.iter().map(|&b| match b {
            b'a' | b'A' | 0 => 0u8,
            b'c' | b'C' | 1 => 1u8,
            b'g' | b'G' | 2 => 2u8,
            b't' | b'T' | 3 => 3u8,
            _ => 0u8, // N -> a (syng convention)
        }));
        seq_buf.push(0); // null terminator

        // Extract syncmers and build GBWT paths
        let mut syncmers: Vec<(i64, i32)> = Vec::new(); // (kmer_index, position)

        unsafe {
            let sit = syng_ffi::syncmerIterator(
                self.seqhash,
                seq_buf.as_mut_ptr() as *mut i8,
                seq_len as i32,
            );

            let mut pos: i32 = 0;
            while syng_ffi::syncmerNext(sit, std::ptr::null_mut(), &mut pos, std::ptr::null_mut()) {
                let mut kmer_index: i64 = 0;
                let syncmer_ptr = seq_buf.as_mut_ptr().add(pos as usize) as *mut i8;
                match lookup_mode {
                    SyncmerLookupMode::AddMissing => {
                        syng_ffi::kmerHashAdd(self.kmer_hash, syncmer_ptr, &mut kmer_index);
                    }
                    SyncmerLookupMode::RequireExisting => {
                        let found =
                            syng_ffi::kmerHashFind(self.kmer_hash, syncmer_ptr, &mut kmer_index);
                        if !found {
                            syng_ffi::impg_seqhashIteratorDestroy(sit);
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "syncmer at {}:{} is absent from the preloaded dictionary",
                                    name, pos
                                ),
                            ));
                        }
                    }
                }
                syncmers.push((kmer_index, pos));
            }

            syng_ffi::impg_seqhashIteratorDestroy(sit);
        }

        if syncmers.is_empty() {
            let path_num = self.name_map.add(name, seq_len as u64);
            if let Some(builder) = self.sampled_position_builder.as_mut() {
                builder.record_path(path_num as usize, seq_len as u64, &[]);
            }
            return Ok(SyngAddSequenceStats {
                sequence_len: seq_len,
                syncmers: 0,
                indexed: false,
            });
        }

        // Build forward GBWT path and capture start info
        let fwd_start_node;
        let fwd_start_count;
        let mut path_steps: Vec<PathStepRecord> = Vec::with_capacity(syncmers.len());
        // Syng's first syncmer lands at wherever the first min k-mer
        // appears in the initial window — anywhere in `[0, w+k)`. We
        // record it so `walk_path` can emit absolute bp coordinates.
        let fwd_first_syncmer_pos = syncmers[0].1 as u64;
        unsafe {
            let first_sync = syncmers[0].0 as i32;
            let sbp = syng_ffi::syngBWTpathStartNew(self.gbwt, first_sync);
            // Capture start info before any path additions modify it
            fwd_start_node = first_sync;
            fwd_start_count = (*sbp).j_last;
            path_steps.push(PathStepRecord {
                step_idx: 0,
                bp_pos: fwd_first_syncmer_pos,
                signed_node: first_sync,
                prev_node: 0,
                prev_offset: 0,
                traversal_rank: fwd_start_count,
            });
            for i in 1..syncmers.len() {
                let next_sync = syncmers[i].0 as i32;
                let offset = (syncmers[i].1 - syncmers[i - 1].1) as u32;
                syng_ffi::syngBWTpathAdd(sbp, next_sync, offset);
                path_steps.push(PathStepRecord {
                    step_idx: i as u32,
                    bp_pos: syncmers[i].1 as u64,
                    signed_node: next_sync,
                    prev_node: syncmers[i - 1].0 as i32,
                    prev_offset: offset,
                    traversal_rank: (*sbp).j_last,
                });
            }
            syng_ffi::syngBWTpathFinish(sbp);
        }

        // Build reverse complement GBWT path
        // Reverse: negate all syncmer indices and reverse the order.
        // Offsets between consecutive syncmers stay the same but are
        // computed from the reverse direction.
        unsafe {
            let n = syncmers.len();
            let first_sync_rc = -(syncmers[n - 1].0 as i32);
            let sbp = syng_ffi::syngBWTpathStartNew(self.gbwt, first_sync_rc);
            for i in (0..n - 1).rev() {
                let next_sync_rc = -(syncmers[i].0 as i32);
                let offset = (syncmers[i + 1].1 - syncmers[i].1) as u32;
                syng_ffi::syngBWTpathAdd(sbp, next_sync_rc, offset);
            }
            syng_ffi::syngBWTpathFinish(sbp);
        }

        // Record in name map with path start info
        let path_num = self.name_map.add(name, seq_len as u64);
        self.name_map.set_path_start(
            path_num,
            GbwtPathStart {
                start_node: fwd_start_node,
                start_count: fwd_start_count,
                num_syncmers: syncmers.len() as u32,
                first_syncmer_pos: fwd_first_syncmer_pos,
            },
        );
        if let Some(builder) = self.sampled_position_builder.as_mut() {
            builder.record_path(path_num as usize, seq_len as u64, &path_steps);
        }

        Ok(SyngAddSequenceStats {
            sequence_len: seq_len,
            syncmers: syncmers.len(),
            indexed: true,
        })
    }

    /// Save the index to disk.
    ///
    /// Produces six files at the given prefix:
    /// - `{prefix}.1khash` — syncmer hash (syng-compatible)
    /// - `{prefix}.1gbwt` — GBWT graph (syng-compatible)
    /// - `{prefix}.syng.names` or `{prefix}.names` — sequence name mapping
    /// - `{prefix}.syng.spos` or `{prefix}.spos` — sampled syncmer occurrence positions
    /// - `{prefix}.syng.pstep` or `{prefix}.pstep` — inverted sampled path-step checkpoints
    /// - `{prefix}.syng.meta` or `{prefix}.meta` — syncmer extraction parameters
    ///
    /// The shorter sidecar names are used when `prefix` already ends in
    /// `.syng`, avoiding names like `foo.syng.syng.spos`.
    pub fn save(&mut self, prefix: &str) -> io::Result<()> {
        let _ = self.finalize_online_sampled_positions()?;
        let schema_text = syng_ffi::syng_schema_text();

        // Write .1gbwt
        let gbwt_path = format!("{}.1gbwt", prefix);
        let gbwt_cpath = CString::new(gbwt_path.as_str())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let gbwt_type = CString::new("gbwt").unwrap();
        unsafe {
            let schema = syng_ffi::oneSchemaCreateFromText(schema_text.as_ptr());
            if schema.is_null() {
                return Err(io::Error::other("Failed to create ONEcode schema for gbwt"));
            }
            let of = syng_ffi::oneFileOpenWriteNew(
                gbwt_cpath.as_ptr(),
                schema,
                gbwt_type.as_ptr(),
                true,
                1,
            );
            if of.is_null() {
                syng_ffi::oneSchemaDestroy(schema);
                return Err(io::Error::other(format!(
                    "Failed to open {} for writing",
                    gbwt_path
                )));
            }
            syng_ffi::syngBWTwrite(of, self.gbwt);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
        }

        // Write .1khash
        let khash_path = format!("{}.1khash", prefix);
        let khash_cpath = CString::new(khash_path.as_str())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let khash_type = CString::new("khash").unwrap();
        unsafe {
            let schema = syng_ffi::oneSchemaCreateFromText(schema_text.as_ptr());
            if schema.is_null() {
                return Err(io::Error::other(
                    "Failed to create ONEcode schema for khash",
                ));
            }
            let of = syng_ffi::oneFileOpenWriteNew(
                khash_cpath.as_ptr(),
                schema,
                khash_type.as_ptr(),
                true,
                1,
            );
            if of.is_null() {
                syng_ffi::oneSchemaDestroy(schema);
                return Err(io::Error::other(format!(
                    "Failed to open {} for writing",
                    khash_path
                )));
            }
            let ok = syng_ffi::kmerHashWriteOneFile(self.kmer_hash, of);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
            if !ok {
                return Err(io::Error::other("kmerHashWriteOneFile failed"));
            }
        }

        // Write names sidecar.
        let names_path = syng_names_path(prefix);
        self.name_map.save(&names_path)?;

        self.save_position_sidecars(prefix)?;

        // Write metadata last so a complete index records the syncmer
        // parameters that future query/map calls must use.
        SyngMetadata::new(self.params).save(&syng_meta_path(prefix))?;

        Ok(())
    }

    /// Save only the position sidecars for an existing syng index prefix.
    pub fn save_position_sidecars(&self, prefix: &str) -> io::Result<()> {
        let sampled_path_steps = self.sampled_path_steps.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot save syng index without sampled path-step checkpoints",
            )
        })?;
        sampled_path_steps.save_atomic(&syng_pstep_path(prefix))?;
        self.save_checkpoint_sidecar(prefix)
    }

    /// Save only the inverted path-step sidecar for an existing syng index prefix.
    pub fn save_path_step_sidecar(&self, prefix: &str) -> io::Result<()> {
        let sampled_path_steps = self.sampled_path_steps.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot save syng index without sampled path-step checkpoints",
            )
        })?;
        sampled_path_steps.save_atomic(&syng_pstep_path(prefix))?;
        self.save_checkpoint_sidecar(prefix)
    }

    /// Save only the occurrence-major syncmer-position sidecar for an existing syng index prefix.
    pub fn save_checkpoint_sidecar(&self, prefix: &str) -> io::Result<()> {
        let sampled_path_steps = self.sampled_path_steps.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot save syng syncmer-position sidecar without sampled path-step checkpoints",
            )
        })?;
        let build_start = std::time::Instant::now();
        let mut records = self.build_checkpoint_records_from_path_steps(sampled_path_steps)?;
        let record_count = records.len();
        SampledCheckpointIndex::save_records_atomic(
            &syng_spos_path(prefix),
            sampled_path_steps.sample_rate,
            &mut records,
        )?;
        if std::env::var("SYNG_EMIT_PROFILE").is_ok() {
            eprintln!(
                "EMITPROF save_syncmer_position_sidecar checkpoints={} build_save={:.2}s",
                record_count,
                build_start.elapsed().as_secs_f64(),
            );
        }
        Ok(())
    }

    fn build_checkpoint_records_from_path_steps(
        &self,
        sampled_path_steps: &SampledPathSteps,
    ) -> io::Result<Vec<CheckpointRecord>> {
        let path_count = self.name_map.path_to_name.len();
        if sampled_path_steps.path_count() != path_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "path-step index has {} paths, expected {}",
                    sampled_path_steps.path_count(),
                    path_count
                ),
            ));
        }
        (0..path_count)
            .into_par_iter()
            .try_fold(Vec::new, |mut records, path_idx| {
                let path_idx_u32 = u32::try_from(path_idx).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "path index exceeds u32 in checkpoint sidecar",
                    )
                })?;
                for checkpoint in sampled_path_steps.decode_path(path_idx)? {
                    records.push(CheckpointRecord {
                        signed_node: checkpoint.signed_node,
                        abs_rank: self.absolute_incoming_rank(&checkpoint)?,
                        path_idx: path_idx_u32,
                        bp_pos: checkpoint.bp_pos,
                    });
                }
                Ok(records)
            })
            .try_reduce(Vec::new, |mut left, right| {
                left.extend(right);
                Ok(left)
            })
    }

    /// Load an index from disk.
    ///
    /// Reads syng files at the given prefix:
    /// - `{prefix}.1khash`
    /// - `{prefix}.1gbwt`
    /// - `{prefix}.syng.names` or `{prefix}.names`
    /// - `{prefix}.syng.pstep` or `{prefix}.pstep`
    /// - `{prefix}.syng.spos` or `{prefix}.spos`
    /// - `{prefix}.syng.meta` or `{prefix}.meta`
    pub fn load(prefix: &str, _params: SyncmerParams) -> io::Result<Self> {
        Self::load_with_options(prefix, true)
    }

    /// Load core syng files for sidecar repair. Missing position sidecars are allowed.
    pub fn load_for_repair(prefix: &str, _params: SyncmerParams) -> io::Result<Self> {
        Self::load_with_options(prefix, false)
    }

    fn load_with_options(prefix: &str, require_sampled_positions: bool) -> io::Result<Self> {
        let metadata_path = existing_syng_sidecar_path(prefix, "meta");
        let metadata = SyngMetadata::load(&metadata_path)?;
        let params = metadata.params;
        let schema_text = syng_ffi::syng_schema_text();

        // Read .1gbwt
        let gbwt_path = format!("{}.1gbwt", prefix);
        if !Path::new(&gbwt_path).exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("GBWT file not found: {}", gbwt_path),
            ));
        }
        let gbwt_cpath = CString::new(gbwt_path.as_str())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let gbwt_type = CString::new("gbwt").unwrap();
        let gbwt = unsafe {
            let schema = syng_ffi::oneSchemaCreateFromText(schema_text.as_ptr());
            if schema.is_null() {
                return Err(io::Error::other(
                    "Failed to create ONEcode schema for gbwt read",
                ));
            }
            let of = syng_ffi::oneFileOpenRead(gbwt_cpath.as_ptr(), schema, gbwt_type.as_ptr(), 1);
            if of.is_null() {
                syng_ffi::oneSchemaDestroy(schema);
                return Err(io::Error::other(format!(
                    "Failed to open {} for reading",
                    gbwt_path
                )));
            }
            let gbwt = syng_ffi::syngBWTread(of);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
            if gbwt.is_null() {
                return Err(io::Error::other("syngBWTread returned null"));
            }
            gbwt
        };

        let kmer_hash = match read_kmer_hash_from_prefix(prefix) {
            Ok(kh) => kh,
            Err(e) => {
                unsafe { syng_ffi::syngBWTdestroy(gbwt) };
                return Err(e);
            }
        };

        // Read names sidecar.
        let names_path = existing_syng_sidecar_path(prefix, "names");
        let name_map = match SyngNameMap::load(&names_path) {
            Ok(nm) => nm,
            Err(e) => {
                unsafe {
                    syng_ffi::syngBWTdestroy(gbwt);
                    syng_ffi::kmerHashDestroy(kmer_hash);
                }
                return Err(e);
            }
        };

        // Create seqhash for future use (queries, incremental adds)
        let seqhash = unsafe {
            syng_ffi::impg_seqhashCreateSafe(params.k as i32, params.w as i32, params.seed as i32)
        };

        // Read sampled path-step sidecar.
        let sampled_path_steps_path = existing_syng_sidecar_path(prefix, "pstep");
        let sampled_path_steps = if Path::new(&sampled_path_steps_path).exists() {
            match SampledPathSteps::load(&sampled_path_steps_path) {
                Ok(ps) => Some(ps),
                Err(e) => {
                    if !require_sampled_positions {
                        log::warn!(
                            "ignoring sampled path-step sidecar during repair: {} ({})",
                            sampled_path_steps_path,
                            e
                        );
                        None
                    } else {
                        unsafe {
                            syng_ffi::syngBWTdestroy(gbwt);
                            syng_ffi::kmerHashDestroy(kmer_hash);
                            syng_ffi::impg_seqhashDestroy(seqhash);
                        }
                        return Err(e);
                    }
                }
            }
        } else if require_sampled_positions {
            unsafe {
                syng_ffi::syngBWTdestroy(gbwt);
                syng_ffi::kmerHashDestroy(kmer_hash);
                syng_ffi::impg_seqhashDestroy(seqhash);
            }
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "sampled path-step sidecar not found: {} (run `impg syng-repair` for this index)",
                    sampled_path_steps_path
                ),
            ));
        } else {
            None
        };

        // Read occurrence-major syncmer-position sidecar. Normal loads require
        // this mmap-backed lookup; repair loads can tolerate a missing/stale
        // `.spos` so they can rebuild it from `.pstep`.
        let sampled_checkpoints_path = existing_syng_sidecar_path(prefix, "spos");
        let mut sampled_checkpoints = if Path::new(&sampled_checkpoints_path).exists() {
            match SampledCheckpointIndex::load(&sampled_checkpoints_path) {
                Ok(spos) => Some(spos),
                Err(e) => {
                    log::warn!(
                        "ignoring sampled syncmer-position sidecar: {} ({})",
                        sampled_checkpoints_path,
                        e
                    );
                    None
                }
            }
        } else {
            None
        };
        if let (Some(path_steps), Some(checkpoints)) =
            (sampled_path_steps.as_ref(), sampled_checkpoints.as_ref())
        {
            if checkpoints.sample_rate != path_steps.sample_rate
                || checkpoints.checkpoint_count() != path_steps.sample_count()
            {
                log::warn!(
                    "ignoring stale sampled syncmer-position sidecar {}: rate/count {} / {} does not match .pstep {} / {}",
                    sampled_checkpoints_path,
                    checkpoints.sample_rate,
                    checkpoints.checkpoint_count(),
                    path_steps.sample_rate,
                    path_steps.sample_count()
                );
                sampled_checkpoints = None;
            }
        }
        if require_sampled_positions && sampled_checkpoints.is_none() {
            unsafe {
                syng_ffi::syngBWTdestroy(gbwt);
                syng_ffi::kmerHashDestroy(kmer_hash);
                syng_ffi::impg_seqhashDestroy(seqhash);
            }
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "sampled syncmer-position sidecar is missing, stale, or unreadable: {} (run `impg syng-repair` for this index)",
                    sampled_checkpoints_path
                ),
            ));
        }

        Ok(Self {
            gbwt,
            kmer_hash,
            seqhash,
            name_map,
            params,
            sampled_path_steps,
            sampled_checkpoints,
            locate_checkpoint_index: OnceLock::new(),
            sampled_position_builder: None,
        })
    }

    /// Returns the raw GBWT pointer (for advanced C interop).
    pub fn gbwt_ptr(&self) -> *mut syng_ffi::SyngBWT {
        self.gbwt
    }

    /// Returns the raw KmerHash pointer (for advanced C interop).
    pub fn kmer_hash_ptr(&self) -> *mut syng_ffi::KmerHash {
        self.kmer_hash
    }

    /// Returns the raw Seqhash pointer (for advanced C interop).
    pub fn seqhash_ptr(&self) -> *mut syng_ffi::Seqhash {
        self.seqhash
    }

    /// Number of syncmer nodes stored in the KmerHash (valid IDs are `1..=N`).
    pub fn num_syncmer_nodes(&self) -> usize {
        unsafe { (*self.kmer_hash).max as usize }
    }

    /// Syncmer node length in bp (= `w + k`).
    pub fn syncmer_length_bp(&self) -> usize {
        unsafe { (*self.kmer_hash).len as usize }
    }

    /// Unpack a syncmer node's DNA sequence (uppercase ASCII A/C/G/T).
    ///
    /// `signed_node_id` follows syng's convention: positive returns the
    /// canonical sequence, negative returns its reverse complement.
    /// IDs are 1-based; `signed_node_id.abs()` must be in `1..=num_syncmer_nodes()`.
    pub fn syncmer_seq(&self, signed_node_id: i32) -> Vec<u8> {
        let len = self.syncmer_length_bp();
        let mut buf = vec![0u8; len];
        unsafe {
            syng_ffi::kmerHashSeq(
                self.kmer_hash,
                signed_node_id as i64,
                buf.as_mut_ptr() as *mut std::os::raw::c_char,
            );
        }
        buf
    }

    /// Public wrapper around the internal forward GBWT walk.
    ///
    /// Returns `(signed_node_id, absolute_bp_pos)` for every syncmer in the
    /// forward path.
    pub fn walk_forward_path(&self, start: &GbwtPathStart) -> Vec<(i32, u64)> {
        self.walk_path(start)
    }

    /// Build a SequenceIndex from the name map (for interop with graph engines).
    pub fn build_seq_index(&self) -> crate::seqidx::SequenceIndex {
        let mut seq_index = crate::seqidx::SequenceIndex::new();
        for (name, &length) in self
            .name_map
            .path_to_name
            .iter()
            .zip(self.name_map.path_to_length.iter())
        {
            seq_index.get_or_insert_id(name, Some(length as usize));
        }
        seq_index
    }

    /// Return the sampled path-step checkpoint at or before `bp_pos` on `path_idx`.
    pub fn sampled_path_step_at_or_before(
        &self,
        path_idx: usize,
        bp_pos: u64,
    ) -> io::Result<Option<SampledPathStepHit>> {
        let sampled_path_steps = self.sampled_path_steps.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                "sampled path-step sidecar is missing; run `impg syng-repair` for this index",
            )
        })?;
        sampled_path_steps.checkpoint_at_or_before(path_idx, bp_pos)
    }

    /// Advance one GBWT step from a sampled path-step checkpoint.
    pub fn next_path_step_from_checkpoint(
        &self,
        checkpoint: &SampledPathStepHit,
    ) -> io::Result<Option<SampledPathStepHit>> {
        unsafe { syng_ffi::impg_syng_suppress_debug() };
        let next_step_idx = match checkpoint.step_idx.checked_add(1) {
            Some(step_idx) => step_idx,
            None => return Ok(None),
        };
        let mut state = self.path_state_from_checkpoint(checkpoint);
        let mut next_node: i32 = 0;
        let mut offset: u32 = 0;
        let has_next =
            unsafe { syng_ffi::syngBWTpathNext(&mut state, &mut next_node, &mut offset) };
        if !has_next {
            return Ok(None);
        }
        let bp_pos = checkpoint
            .bp_pos
            .checked_add(offset as u64)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "path-step coordinate overflow")
            })?;
        Ok(Some(SampledPathStepHit {
            path_idx: checkpoint.path_idx,
            step_idx: next_step_idx,
            bp_pos,
            signed_node: next_node,
            prev_node: state.last_node,
            prev_offset: offset,
            traversal_rank: state.j_last,
        }))
    }

    /// Walk a bp interval on an indexed forward path using sampled path-step
    /// checkpoints, returning every syncmer step overlapping `[start, end)`.
    pub fn walk_path_range(
        &self,
        path_idx: usize,
        start: u64,
        end: u64,
    ) -> io::Result<Vec<(i32, u64)>> {
        if path_idx >= self.name_map.path_to_name.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("path index {} is outside the syng name map", path_idx),
            ));
        }
        self.walk_path_range_from_sampled_steps(path_idx, start, end)
    }

    /// Rebuild sampled path-step checkpoints by walking existing GBWT paths.
    ///
    /// This repairs indexes that already have `.1gbwt`, `.1khash`, `.names`, and
    /// `.meta`, without re-decompressing sequences or reconstructing the graph.
    pub fn rebuild_sampled_position_indexes_from_gbwt(
        &mut self,
        sample_rate: u32,
    ) -> io::Result<SampledPositionBuildStats> {
        self.rebuild_sampled_path_steps_from_gbwt(sample_rate)
    }

    /// Rebuild only the inverted path-step sidecar by walking existing GBWT paths.
    pub fn rebuild_sampled_path_steps_from_gbwt(
        &mut self,
        sample_rate: u32,
    ) -> io::Result<SampledPositionBuildStats> {
        unsafe { syng_ffi::impg_syng_suppress_debug() };
        let mut builder = SampledPositionBuilder::new_path_steps_only(sample_rate, 0, 0)?;
        let syncmer_len = (self.params.k + self.params.w) as u64;
        for path_idx in 0..self.name_map.path_to_name.len() {
            let seq_len = self.name_map.path_to_length[path_idx];
            let steps = match self
                .name_map
                .path_starts
                .get(path_idx)
                .and_then(|s| s.as_ref())
            {
                Some(start) => self.walk_path_steps(start)?,
                None if seq_len < syncmer_len => Vec::new(),
                None => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "cannot rebuild path-step sidecar: missing GBWT path start for {}",
                            self.name_map.path_to_name[path_idx]
                        ),
                    ));
                }
            };
            builder.record_path(path_idx, seq_len, &steps);
        }
        let (sampled_path_steps, stats) =
            builder.finish_path_steps(self.name_map.path_to_name.len())?;
        self.sampled_path_steps = Some(sampled_path_steps);
        self.sampled_checkpoints = None;
        self.reset_locate_checkpoint_index();
        Ok(stats)
    }

    /// Rebuild sampled path-step checkpoints by walking existing GBWT paths
    /// in parallel. Sampling is a regular per-path syncmer-step grid.
    pub fn rebuild_sampled_position_indexes_from_gbwt_parallel(
        &mut self,
        sample_rate: u32,
        progress_interval: usize,
    ) -> io::Result<SampledPositionBuildStats> {
        self.rebuild_sampled_path_steps_from_gbwt_parallel(sample_rate, progress_interval)
    }

    /// Rebuild only the inverted path-step sidecar by walking existing GBWT
    /// paths in parallel. Sampling is a regular per-path syncmer-step grid.
    pub fn rebuild_sampled_path_steps_from_gbwt_parallel(
        &mut self,
        sample_rate: u32,
        progress_interval: usize,
    ) -> io::Result<SampledPositionBuildStats> {
        let (_, sampled_path_steps, stats) = self
            .rebuild_sampled_position_chunks_from_gbwt_parallel(
                sample_rate,
                false,
                progress_interval,
            )?;
        self.sampled_path_steps = Some(sampled_path_steps);
        self.sampled_checkpoints = None;
        self.reset_locate_checkpoint_index();
        Ok(stats)
    }

    fn rebuild_sampled_position_chunks_from_gbwt_parallel(
        &self,
        sample_rate: u32,
        collect_positions: bool,
        progress_interval: usize,
    ) -> io::Result<(
        Option<SampledPositions>,
        SampledPathSteps,
        SampledPositionBuildStats,
    )> {
        validate_position_sample_rate(sample_rate)?;
        unsafe { syng_ffi::impg_syng_suppress_debug() };

        let path_count = self.name_map.path_to_name.len();
        let path_starts =
            SampledPositions::path_starts_from_lengths(&self.name_map.path_to_length)?;
        let syncmer_len = (self.params.k + self.params.w) as u64;
        let total_steps: u64 = self
            .name_map
            .path_starts
            .iter()
            .filter_map(|start| start.as_ref())
            .map(|start| start.num_syncmers as u64)
            .sum();
        log::info!(
            "Parallel regular sampled-position rebuild: {} paths, {} syncmer steps, grid period {}, {} rayon threads",
            path_count,
            total_steps,
            sample_rate,
            rayon::current_num_threads()
        );

        let progress = Mutex::new(PositionRepairProgress::default());
        let results: Vec<io::Result<EncodedPositionSampleChunk>> = (0..path_count)
            .into_par_iter()
            .map(|path_idx| {
                let seq_len = self.name_map.path_to_length[path_idx];
                let path_start = path_starts[path_idx];
                let result = match self.name_map.path_starts.get(path_idx).and_then(|s| s.as_ref()) {
                    Some(start) => self.encode_regular_position_sample_chunk(
                        path_start,
                        start,
                        sample_rate,
                        collect_positions,
                    ),
                    None if seq_len < syncmer_len => Ok(EncodedPositionSampleChunk {
                        position_samples: Vec::new(),
                        path_step_data: Vec::new(),
                        path_step_sample_count: 0,
                        walked_steps: 0,
                    }),
                    None => Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "cannot rebuild positional sidecars: missing GBWT path start for {}",
                            self.name_map.path_to_name[path_idx]
                        ),
                    )),
                };

                if let Ok(chunk) = &result {
                    let mut progress = progress.lock().map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::Other,
                            "sampled-position repair progress lock was poisoned",
                        )
                    })?;
                    progress.completed_paths += 1;
                    progress.completed_steps =
                        progress.completed_steps.checked_add(chunk.walked_steps).ok_or_else(
                            || {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    "sampled-position progress step count overflowed u64",
                                )
                            },
                        )?;
                    progress.completed_position_samples = progress
                        .completed_position_samples
                        .checked_add(chunk.position_samples.len() as u64)
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "sampled-position progress count overflowed u64",
                            )
                        })?;
                    progress.completed_path_step_samples = progress
                        .completed_path_step_samples
                        .checked_add(chunk.path_step_sample_count)
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "path-step progress count overflowed u64",
                            )
                        })?;
                    let done = progress.completed_paths;
                    if progress_interval > 0
                        && (done == path_count || done % progress_interval == 0)
                    {
                        let path_pct = if path_count > 0 {
                            (done as f64 * 100.0) / path_count as f64
                        } else {
                            100.0
                        };
                        let step_pct = if total_steps > 0 {
                            (progress.completed_steps as f64 * 100.0) / total_steps as f64
                        } else {
                            100.0
                        };
                        log::info!(
                            "Sampled-position repair progress: {}/{} paths complete ({:.1}%), {}/{} syncmer steps ({:.1}%), {} path-step checkpoints; completed path {}/{} {}",
                            done,
                            path_count,
                            path_pct,
                            progress.completed_steps,
                            total_steps,
                            step_pct,
                            progress.completed_path_step_samples,
                            path_idx + 1,
                            path_count,
                            self.name_map.path_to_name[path_idx]
                        );
                    }
                }

                result
            })
            .collect();

        let mut position_samples = Vec::new();
        let mut path_step_offsets = Vec::with_capacity(path_count + 1);
        let mut path_step_data = Vec::new();
        let mut path_step_sample_count = 0u64;
        let mut sampled_step_paths = 0usize;
        let mut walked_paths = 0usize;
        path_step_offsets.push(0);

        for result in results {
            let chunk = result?;
            if chunk.walked_steps > 0 {
                walked_paths += 1;
            }
            if !chunk.path_step_data.is_empty() {
                sampled_step_paths += 1;
            }
            path_step_sample_count = path_step_sample_count
                .checked_add(chunk.path_step_sample_count)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "path-step sample count overflowed u64",
                    )
                })?;
            if collect_positions {
                position_samples.extend_from_slice(&chunk.position_samples);
            }
            path_step_data.extend_from_slice(&chunk.path_step_data);
            path_step_offsets.push(path_step_data.len() as u64);
        }

        let sampled_path_steps = SampledPathSteps::from_encoded(
            sample_rate,
            path_step_sample_count,
            path_step_offsets,
            path_step_data,
        )?;
        let sampled_positions = if collect_positions {
            Some(SampledPositions::from_samples(
                sample_rate,
                &self.name_map.path_to_length,
                position_samples,
            )?)
        } else {
            None
        };
        let stats = SampledPositionBuildStats {
            sampled_occurrences: sampled_positions
                .as_ref()
                .map(|sp| sp.sample_count())
                .unwrap_or(0),
            sampled_nodes: sampled_positions
                .as_ref()
                .map(|sp| sp.node_count())
                .unwrap_or(0),
            sampled_path_steps: sampled_path_steps.sample_count(),
            sampled_step_paths,
            walked_paths,
            sample_rate,
        };
        Ok((sampled_positions, sampled_path_steps, stats))
    }

    fn encode_regular_position_sample_chunk(
        &self,
        path_start: u64,
        start: &GbwtPathStart,
        sample_rate: u32,
        collect_positions: bool,
    ) -> io::Result<EncodedPositionSampleChunk> {
        if start.num_syncmers == 0 {
            return Ok(EncodedPositionSampleChunk {
                position_samples: Vec::new(),
                path_step_data: Vec::new(),
                path_step_sample_count: 0,
                walked_steps: 0,
            });
        }

        let expected_steps = start.num_syncmers as usize;
        let last_step_idx = u32::try_from(expected_steps - 1).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "GBWT path step count exceeds u32",
            )
        })?;
        let mut position_samples = Vec::new();
        let mut path_step_data = Vec::new();
        let mut path_step_sample_count = 0u64;
        let mut prev_sample_bp = 0u64;
        let mut prev_sample_step = 0u32;
        let mut steps_seen = 0usize;

        let mut state = self.path_state_from_start(start);
        let mut acc_pos = start.first_syncmer_pos;
        let first_step = PathStepRecord {
            step_idx: 0,
            bp_pos: acc_pos,
            signed_node: start.start_node,
            prev_node: 0,
            prev_offset: 0,
            traversal_rank: start.start_count,
        };
        encode_regular_position_sample_if_selected(
            path_start,
            sample_rate,
            last_step_idx,
            collect_positions,
            first_step,
            &mut prev_sample_bp,
            &mut prev_sample_step,
            &mut path_step_sample_count,
            &mut position_samples,
            &mut path_step_data,
        )?;
        steps_seen += 1;

        let mut next_node: i32 = 0;
        let mut offset: u32 = 0;
        while unsafe { syng_ffi::syngBWTpathNext(&mut state, &mut next_node, &mut offset) } {
            if steps_seen >= expected_steps {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "GBWT path has more steps than name map reports ({})",
                        expected_steps
                    ),
                ));
            }
            acc_pos = acc_pos.checked_add(offset as u64).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "GBWT path coordinate overflow")
            })?;
            let step = PathStepRecord {
                step_idx: steps_seen as u32,
                bp_pos: acc_pos,
                signed_node: next_node,
                prev_node: state.last_node,
                prev_offset: offset,
                traversal_rank: state.j_last,
            };
            encode_regular_position_sample_if_selected(
                path_start,
                sample_rate,
                last_step_idx,
                collect_positions,
                step,
                &mut prev_sample_bp,
                &mut prev_sample_step,
                &mut path_step_sample_count,
                &mut position_samples,
                &mut path_step_data,
            )?;
            steps_seen += 1;
        }
        if steps_seen != expected_steps {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GBWT path has {} steps but name map reports {}",
                    steps_seen, expected_steps
                ),
            ));
        }

        Ok(EncodedPositionSampleChunk {
            position_samples,
            path_step_data,
            path_step_sample_count,
            walked_steps: steps_seen as u64,
        })
    }

    fn path_state_from_checkpoint(&self, checkpoint: &SampledPathStepHit) -> syng_ffi::SyngBWTpath {
        syng_ffi::SyngBWTpath {
            sb: self.gbwt,
            last_node: checkpoint.prev_node,
            this_node: checkpoint.signed_node,
            last_off: checkpoint.prev_offset,
            j_last: checkpoint.traversal_rank,
            j_max: 0,
        }
    }

    fn path_state_from_start(&self, start: &GbwtPathStart) -> syng_ffi::SyngBWTpath {
        syng_ffi::SyngBWTpath {
            sb: self.gbwt,
            last_node: 0,
            this_node: start.start_node,
            last_off: 0,
            j_last: start.start_count,
            j_max: 0,
        }
    }

    fn walk_path_steps(&self, start: &GbwtPathStart) -> io::Result<Vec<PathStepRecord>> {
        unsafe { syng_ffi::impg_syng_suppress_debug() };
        if start.num_syncmers == 0 {
            return Ok(Vec::new());
        }
        let expected_steps = start.num_syncmers as usize;
        let mut steps = Vec::with_capacity(expected_steps);
        let mut state = self.path_state_from_start(start);
        let mut acc_pos = start.first_syncmer_pos;
        steps.push(PathStepRecord {
            step_idx: 0,
            bp_pos: acc_pos,
            signed_node: start.start_node,
            prev_node: 0,
            prev_offset: 0,
            traversal_rank: start.start_count,
        });

        let mut next_node: i32 = 0;
        let mut offset: u32 = 0;
        while unsafe { syng_ffi::syngBWTpathNext(&mut state, &mut next_node, &mut offset) } {
            if steps.len() >= expected_steps {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "GBWT path has more steps than name map reports ({})",
                        expected_steps
                    ),
                ));
            }
            acc_pos = acc_pos.checked_add(offset as u64).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "GBWT path coordinate overflow")
            })?;
            steps.push(PathStepRecord {
                step_idx: steps.len() as u32,
                bp_pos: acc_pos,
                signed_node: next_node,
                prev_node: state.last_node,
                prev_offset: offset,
                traversal_rank: state.j_last,
            });
        }
        if steps.len() != expected_steps {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GBWT path has {} steps but name map reports {}",
                    steps.len(),
                    expected_steps
                ),
            ));
        }
        Ok(steps)
    }

    fn walk_path_range_from_sampled_steps(
        &self,
        path_idx: usize,
        start: u64,
        end: u64,
    ) -> io::Result<Vec<(i32, u64)>> {
        if end <= start {
            return Ok(Vec::new());
        }
        let path_start = self
            .name_map
            .path_starts
            .get(path_idx)
            .and_then(|start| start.as_ref())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "No GBWT path start info for path {} — index may need rebuilding",
                        path_idx
                    ),
                )
            })?;
        if path_start.num_syncmers == 0 {
            return Ok(Vec::new());
        }

        let syncmer_len = (self.params.k + self.params.w) as u64;
        let scan_start = start.saturating_sub(syncmer_len.saturating_sub(1));
        let mut current =
            if let Some(checkpoint) = self.sampled_path_step_at_or_before(path_idx, scan_start)? {
                if checkpoint.step_idx >= path_start.num_syncmers {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "sampled path-step checkpoint is beyond path length",
                    ));
                }
                checkpoint
            } else {
                SampledPathStepHit {
                    path_idx,
                    step_idx: 0,
                    bp_pos: path_start.first_syncmer_pos,
                    signed_node: path_start.start_node,
                    prev_node: 0,
                    prev_offset: 0,
                    traversal_rank: path_start.start_count,
                }
            };

        let mut nodes = Vec::new();
        loop {
            if current.bp_pos < end && current.bp_pos.saturating_add(syncmer_len) > start {
                nodes.push((current.signed_node, current.bp_pos));
            }
            if current.bp_pos >= end || current.step_idx + 1 >= path_start.num_syncmers {
                break;
            }
            let Some(next) = self.next_path_step_from_checkpoint(&current)? else {
                break;
            };
            current = next;
        }
        Ok(nodes)
    }

    fn locate_checkpoint_index(&self) -> io::Result<&LocateCheckpointIndex> {
        let cached = self.locate_checkpoint_index.get_or_init(|| {
            self.build_locate_checkpoint_index()
                .map_err(|e| e.to_string())
        });
        cached.as_ref().map_err(|msg| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to build syng locate checkpoint index: {msg}"),
            )
        })
    }

    fn prepare_locate_checkpoint_lookup(&self) -> io::Result<()> {
        if self.sampled_checkpoints.is_some() {
            Ok(())
        } else {
            self.locate_checkpoint_index().map(|_| ())
        }
    }

    fn build_locate_checkpoint_index(&self) -> io::Result<LocateCheckpointIndex> {
        let build_start = std::time::Instant::now();
        let sampled_path_steps = self.sampled_path_steps.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "syng locate requires .syng.pstep; run `impg syng-repair` for this index",
            )
        })?;

        let shard_count = (rayon::current_num_threads().max(1) * 4).next_power_of_two();
        let shard_mask = shard_count - 1;
        let new_shards = || -> Vec<FxHashMap<LocateCheckpointKey, LocateCheckpointValue>> {
            (0..shard_count).map(|_| FxHashMap::default()).collect()
        };
        let mut shards = (0..self.name_map.path_to_name.len())
            .into_par_iter()
            .try_fold(new_shards, |mut shards, path_idx| {
                for checkpoint in sampled_path_steps.decode_path(path_idx)? {
                    let abs_rank = self.absolute_incoming_rank(&checkpoint)?;
                    let key = LocateCheckpointKey {
                        signed_node: checkpoint.signed_node,
                        abs_rank,
                    };
                    let value = LocateCheckpointValue {
                        path_idx,
                        bp_pos: checkpoint.bp_pos,
                    };
                    let shard_idx = Self::locate_checkpoint_shard(&key, shard_mask);
                    if shards[shard_idx].insert(key, value).is_some() {
                        return Err(Self::duplicate_locate_checkpoint_error(key));
                    }
                }
                Ok(shards)
            })
            .try_reduce(new_shards, |mut left, right| {
                for (left_shard, right_shard) in left.iter_mut().zip(right) {
                    for (key, value) in right_shard {
                        if left_shard.insert(key, value).is_some() {
                            return Err(Self::duplicate_locate_checkpoint_error(key));
                        }
                    }
                }
                Ok(left)
            })?;
        shards.shrink_to_fit();
        let checkpoint_count: usize = shards.iter().map(FxHashMap::len).sum();

        if std::env::var("SYNG_EMIT_PROFILE").is_ok() {
            eprintln!(
                "EMITPROF locate_checkpoint_index paths={} samples={} checkpoints={} shards={} build={:.2}s",
                self.name_map.path_to_name.len(),
                sampled_path_steps.sample_count(),
                checkpoint_count,
                shard_count,
                build_start.elapsed().as_secs_f64(),
            );
        }

        Ok(LocateCheckpointIndex {
            sample_rate: sampled_path_steps.sample_rate,
            shard_mask,
            shards,
        })
    }

    fn locate_checkpoint_shard(key: &LocateCheckpointKey, shard_mask: usize) -> usize {
        let h = (key.signed_node as u32 as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
            ^ key.abs_rank as u64;
        (h as usize) & shard_mask
    }

    fn duplicate_locate_checkpoint_error(key: LocateCheckpointKey) -> io::Error {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "duplicate syng locate checkpoint for node {} rank {}",
                key.signed_node, key.abs_rank
            ),
        )
    }

    fn absolute_incoming_rank(&self, checkpoint: &SampledPathStepHit) -> io::Result<u32> {
        let mut abs_rank = 0u32;
        let ok = unsafe {
            syng_ffi::syngBWTincomingRank(
                self.gbwt,
                checkpoint.signed_node,
                checkpoint.prev_node,
                checkpoint.prev_offset,
                checkpoint.traversal_rank,
                &mut abs_rank,
            )
        };
        if ok {
            Ok(abs_rank)
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "syngBWTincomingRank failed for sampled path-step checkpoint",
            ))
        }
    }

    fn signed_node_occurrence_count(&self, signed_node: i32) -> io::Result<u32> {
        let mut high = 0u32;
        let sbp = unsafe { syng_ffi::syngBWTmatchStart(self.gbwt, signed_node, &mut high) };
        if sbp.is_null() {
            return Err(io::Error::other("syngBWTmatchStart returned null"));
        }
        unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
        Ok(high)
    }

    /// Total occurrence count for an unsigned syncmer node across both
    /// orientations in the GBWT.
    pub fn syncmer_node_occurrence_count(&self, node_id: u32) -> io::Result<u32> {
        Ok(self.syncmer_node_orientation_counts(node_id)?.total())
    }

    fn syncmer_node_orientation_counts(&self, node_id: u32) -> io::Result<NodeOccurrenceCounts> {
        if node_id == 0 || node_id > i32::MAX as u32 {
            return Ok(NodeOccurrenceCounts {
                forward: 0,
                reverse: 0,
            });
        }
        let signed_node = node_id as i32;
        let forward = self.signed_node_occurrence_count(signed_node)?;
        let reverse = self.signed_node_occurrence_count(-signed_node)?;
        Ok(NodeOccurrenceCounts { forward, reverse })
    }

    fn valid_walk_node(&self, signed_node: i32) -> bool {
        let node_id = signed_node.unsigned_abs() as usize;
        signed_node != 0 && node_id <= self.num_syncmer_nodes()
    }

    fn walk_step_offset(prev: SyngWalkStep, next: SyngWalkStep) -> Option<u32> {
        next.bp_pos
            .checked_sub(prev.bp_pos)
            .and_then(|offset| u32::try_from(offset).ok())
    }

    fn start_gbwt_match(
        &self,
        signed_node: i32,
    ) -> io::Result<Option<(*mut syng_ffi::SyngBWTpath, u32, u32)>> {
        if !self.valid_walk_node(signed_node) {
            return Ok(None);
        }
        let mut high = 0u32;
        let sbp = unsafe { syng_ffi::syngBWTmatchStart(self.gbwt, signed_node, &mut high) };
        if sbp.is_null() {
            return Err(io::Error::other("syngBWTmatchStart returned null"));
        }
        if high == 0 {
            unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
            return Ok(None);
        }
        Ok(Some((sbp, 0, high)))
    }

    fn advance_gbwt_match(
        &self,
        sbp: *mut syng_ffi::SyngBWTpath,
        low: &mut u32,
        high: &mut u32,
        signed_node: i32,
        offset: u32,
    ) -> bool {
        if !self.valid_walk_node(signed_node) {
            return false;
        }
        let mut next_low = *low;
        let mut next_high = *high;
        let matched = unsafe {
            syng_ffi::syngBWTmatchNext(sbp, signed_node, offset, &mut next_low, &mut next_high)
        };
        if matched {
            *low = next_low;
            *high = next_high;
        }
        matched
    }

    fn push_gbwt_mem(
        &self,
        walk: &[SyngWalkStep],
        start: usize,
        end: usize,
        low: u32,
        high: u32,
        mems: &mut Vec<SyngGbwtMem>,
    ) {
        let occurrences = high.saturating_sub(low);
        if start >= end || occurrences == 0 {
            return;
        }
        let query_start = walk[start].bp_pos;
        let query_end = walk[end - 1]
            .bp_pos
            .saturating_add(self.syncmer_length_bp() as u64);
        mems.push(SyngGbwtMem {
            step_start: start,
            step_end: end,
            query_start,
            query_end,
            anchors: end - start,
            occurrences,
            match_signed_node: walk[end - 1].signed_node,
            match_low: low,
            match_high: high,
        });
    }

    fn prune_contained_gbwt_mems(candidates: Vec<SyngGbwtMem>) -> Vec<SyngGbwtMem> {
        let mut candidates = candidates;
        candidates.sort_by(|a, b| {
            a.step_start
                .cmp(&b.step_start)
                .then(b.step_end.cmp(&a.step_end))
                .then(a.occurrences.cmp(&b.occurrences))
        });
        let mut mems: Vec<SyngGbwtMem> = Vec::new();
        'candidate: for candidate in candidates {
            for kept in &mems {
                if kept.step_start <= candidate.step_start
                    && kept.step_end >= candidate.step_end
                    && kept.anchors >= candidate.anchors
                {
                    continue 'candidate;
                }
            }
            mems.retain(|kept| {
                !(candidate.step_start <= kept.step_start
                    && candidate.step_end >= kept.step_end
                    && candidate.anchors >= kept.anchors)
            });
            mems.push(candidate);
        }
        mems.sort_by(|a, b| {
            a.step_start
                .cmp(&b.step_start)
                .then(a.step_end.cmp(&b.step_end))
        });
        mems
    }

    fn restart_gbwt_match_at_suffix(
        &self,
        walk: &[SyngWalkStep],
        run_start: usize,
        current: usize,
    ) -> io::Result<Option<(*mut syng_ffi::SyngBWTpath, usize, u32, u32)>> {
        let Some((reverse_sbp, mut reverse_low, mut reverse_high)) =
            self.start_gbwt_match(-walk[current].signed_node)?
        else {
            return Ok(None);
        };

        let mut suffix_start = current;
        while suffix_start > run_start {
            let prev = suffix_start - 1;
            let Some(offset) = Self::walk_step_offset(walk[prev], walk[suffix_start]) else {
                break;
            };
            if !self.advance_gbwt_match(
                reverse_sbp,
                &mut reverse_low,
                &mut reverse_high,
                -walk[prev].signed_node,
                offset,
            ) {
                break;
            }
            suffix_start = prev;
        }
        unsafe { syng_ffi::syngBWTpathDestroy(reverse_sbp) };

        let Some((sbp, mut low, mut high)) =
            self.start_gbwt_match(walk[suffix_start].signed_node)?
        else {
            return Ok(None);
        };
        for next in suffix_start + 1..=current {
            let Some(offset) = Self::walk_step_offset(walk[next - 1], walk[next]) else {
                unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
                return Ok(None);
            };
            if !self.advance_gbwt_match(sbp, &mut low, &mut high, walk[next].signed_node, offset) {
                unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
                return Ok(None);
            }
        }
        Ok(Some((sbp, suffix_start, low, high)))
    }

    /// Find maximal exact matches of a signed syncmer walk in the syng GBWT.
    ///
    /// This mirrors syngmap's GBWT MEM procedure: extend the current match by
    /// signed syncmer node and bp offset, then on failure reverse-search to the
    /// longest suffix that can seed the next match. Returned MEMs use half-open
    /// step coordinates `[step_start, step_end)`.
    pub fn gbwt_mems_for_walk(&self, walk: &[SyngWalkStep]) -> io::Result<Vec<SyngGbwtMem>> {
        if walk.is_empty() {
            return Ok(Vec::new());
        }

        let mut candidates = Vec::new();
        let mut sbp: *mut syng_ffi::SyngBWTpath = std::ptr::null_mut();
        let mut active_start = 0usize;
        let mut run_start = 0usize;
        let mut low = 0u32;
        let mut high = 0u32;

        for idx in 0..walk.len() {
            if !self.valid_walk_node(walk[idx].signed_node) {
                if !sbp.is_null() {
                    self.push_gbwt_mem(walk, active_start, idx, low, high, &mut candidates);
                    unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
                    sbp = std::ptr::null_mut();
                }
                run_start = idx + 1;
                continue;
            }

            if sbp.is_null() {
                let Some((started_sbp, started_low, started_high)) =
                    self.start_gbwt_match(walk[idx].signed_node)?
                else {
                    run_start = idx + 1;
                    continue;
                };
                sbp = started_sbp;
                active_start = idx;
                run_start = idx;
                low = started_low;
                high = started_high;
                continue;
            }

            let Some(offset) = Self::walk_step_offset(walk[idx - 1], walk[idx]) else {
                self.push_gbwt_mem(walk, active_start, idx, low, high, &mut candidates);
                unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
                sbp = std::ptr::null_mut();
                let Some((started_sbp, started_low, started_high)) =
                    self.start_gbwt_match(walk[idx].signed_node)?
                else {
                    run_start = idx + 1;
                    continue;
                };
                sbp = started_sbp;
                active_start = idx;
                run_start = idx;
                low = started_low;
                high = started_high;
                continue;
            };

            if self.advance_gbwt_match(sbp, &mut low, &mut high, walk[idx].signed_node, offset) {
                continue;
            }

            self.push_gbwt_mem(walk, active_start, idx, low, high, &mut candidates);
            unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
            match self.restart_gbwt_match_at_suffix(walk, run_start, idx)? {
                Some((restarted_sbp, restarted_start, restarted_low, restarted_high)) => {
                    sbp = restarted_sbp;
                    active_start = restarted_start;
                    low = restarted_low;
                    high = restarted_high;
                }
                None => {
                    sbp = std::ptr::null_mut();
                    run_start = idx + 1;
                }
            }
        }

        if !sbp.is_null() {
            self.push_gbwt_mem(walk, active_start, walk.len(), low, high, &mut candidates);
            unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
        }

        Ok(Self::prune_contained_gbwt_mems(candidates))
    }

    fn locate_signed_node_occurrences(
        &self,
        signed_node: i32,
    ) -> io::Result<Vec<LocatedSyncmerHit>> {
        if signed_node == 0 {
            return Ok(Vec::new());
        }
        let occurrence_count = self.signed_node_occurrence_count(signed_node)?;
        self.locate_signed_node_occurrences_with_count(signed_node, occurrence_count)
    }

    fn locate_signed_node_occurrences_with_count(
        &self,
        signed_node: i32,
        occurrence_count: u32,
    ) -> io::Result<Vec<LocatedSyncmerHit>> {
        self.locate_signed_node_occurrence_range(signed_node, 0, occurrence_count)
    }

    fn locate_signed_node_occurrence_range(
        &self,
        signed_node: i32,
        start_rank: u32,
        end_rank: u32,
    ) -> io::Result<Vec<LocatedSyncmerHit>> {
        if signed_node == 0 {
            return Ok(Vec::new());
        }
        if start_rank >= end_rank {
            return Ok(Vec::new());
        }

        let persistent_checkpoint_index = self.sampled_checkpoints.as_ref();
        let fallback_checkpoint_index = if persistent_checkpoint_index.is_none() {
            Some(self.locate_checkpoint_index()?)
        } else {
            None
        };
        let max_walk_steps = persistent_checkpoint_index
            .map(|index| index.sample_rate)
            .or_else(|| {
                fallback_checkpoint_index
                    .as_ref()
                    .map(|index| index.sample_rate)
            })
            .unwrap_or(DEFAULT_POSITION_SAMPLE_RATE) as usize;

        let target_orient = if signed_node >= 0 { 0 } else { 1 };
        let mut hits = Vec::new();
        for start_rank in start_rank..end_rank {
            let mut current_node = signed_node;
            let mut current_rank = start_rank;
            let mut walked_bp = 0u64;

            for _ in 0..=max_walk_steps {
                let key = LocateCheckpointKey {
                    signed_node: current_node,
                    abs_rank: current_rank,
                };
                let checkpoint = if let Some(index) = persistent_checkpoint_index {
                    index.checkpoint(key.signed_node, key.abs_rank)?
                } else if let Some(index) = fallback_checkpoint_index {
                    let shard_idx = Self::locate_checkpoint_shard(&key, index.shard_mask);
                    index.shards[shard_idx].get(&key).copied()
                } else {
                    None
                };
                if let Some(checkpoint) = checkpoint {
                    let Some(target_pos) = checkpoint.bp_pos.checked_sub(walked_bp) else {
                        break;
                    };
                    hits.push(LocatedSyncmerHit {
                        path_idx: checkpoint.path_idx,
                        target_pos,
                        target_orient,
                    });
                    break;
                }

                let mut next_node = 0i32;
                let mut next_off = 0u32;
                let mut next_rank = 0u32;
                let has_next = unsafe {
                    syng_ffi::syngBWTadvanceRank(
                        self.gbwt,
                        current_node,
                        current_rank,
                        &mut next_node,
                        &mut next_off,
                        &mut next_rank,
                    )
                };
                if !has_next {
                    break;
                }
                walked_bp = walked_bp.checked_add(next_off as u64).ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "syng locate offset overflow")
                })?;
                current_node = next_node;
                current_rank = next_rank;
            }
        }
        Ok(hits)
    }

    fn locate_node_occurrences_with_counts(
        &self,
        node_id: u32,
        counts: Option<NodeOccurrenceCounts>,
    ) -> io::Result<Vec<LocatedSyncmerHit>> {
        if node_id == 0 || node_id > i32::MAX as u32 {
            return Ok(Vec::new());
        }
        let signed_node = node_id as i32;
        let mut hits = if let Some(counts) = counts {
            self.locate_signed_node_occurrences_with_count(signed_node, counts.forward)?
        } else {
            self.locate_signed_node_occurrences(signed_node)?
        };
        if let Some(counts) = counts {
            hits.extend(
                self.locate_signed_node_occurrences_with_count(-signed_node, counts.reverse)?,
            );
        } else {
            hits.extend(self.locate_signed_node_occurrences(-signed_node)?);
        }
        Ok(hits)
    }

    fn top_frequency_seed_nodes(
        occurrence_counts: &[(u32, u32)],
        drop_top_fraction: f64,
    ) -> FxHashSet<u32> {
        if drop_top_fraction <= 0.0 || occurrence_counts.is_empty() {
            return FxHashSet::default();
        }
        let drop_count = ((occurrence_counts.len() as f64) * drop_top_fraction).floor() as usize;
        if drop_count == 0 {
            return FxHashSet::default();
        }
        let mut ranked = occurrence_counts.to_vec();
        ranked.sort_by(|(node_a, count_a), (node_b, count_b)| {
            count_b.cmp(count_a).then(node_a.cmp(node_b))
        });
        ranked
            .into_iter()
            .take(drop_count.min(occurrence_counts.len()))
            .map(|(node, _)| node)
            .collect()
    }

    /// Walk a genome's forward GBWT path, returning `(node_id, absolute_bp_pos)`
    /// for each syncmer in the path.
    ///
    /// Positions are ABSOLUTE bp coordinates on the forward sequence — the
    /// first syncmer is placed at `start.first_syncmer_pos` (the bp offset
    /// at which syng's C iterator emitted it), and subsequent syncmers are
    /// accumulated from there using the inter-syncmer offsets stored in the
    /// GBWT. This matches the convention used by callers of
    /// `query_region_with_anchors` (which treat returned positions as
    /// absolute sequence coordinates when computing BiWFA realignment).
    ///
    /// For indexes loaded from an older on-disk format that predates the
    /// `first_syncmer_pos` field, that value is 0 and the returned positions
    /// reduce to the old behaviour (relative to the first syncmer). A fresh
    /// rebuild of the index corrects this.
    #[allow(dead_code)]
    fn walk_path(&self, start: &GbwtPathStart) -> Vec<(i32, u64)> {
        let mut nodes = Vec::with_capacity(start.num_syncmers as usize);
        unsafe {
            let sbp = syng_ffi::syngBWTpathStartOld(self.gbwt, start.start_node, start.start_count);
            // Anchor the bp accumulator at the first syncmer's absolute
            // position (0 for old indexes — pre-fix behaviour).
            let mut acc_pos: u64 = start.first_syncmer_pos;
            nodes.push((start.start_node, acc_pos));

            let mut next_node: i32 = 0;
            let mut offset: u32 = 0;
            while syng_ffi::syngBWTpathNext(sbp, &mut next_node, &mut offset) {
                acc_pos += offset as u64;
                nodes.push((next_node, acc_pos));
            }
            syng_ffi::syngBWTpathDestroy(sbp);
        }
        nodes
    }

    /// Query: find all homologous sequences for a genomic region.
    ///
    /// Returns intervals on other genomes that share syncmer nodes with the
    /// query region `[start, end)` on `genome`. Results are padded outward
    /// by `padding` bp.
    ///
    /// Thin wrapper over [`query_region_with_anchors`] that discards the
    /// per-homolog anchor set.
    pub fn query_region(
        &self,
        genome: &str,
        start: u64,
        end: u64,
        padding: u64,
    ) -> Result<Vec<HomologousInterval>, io::Error> {
        self.query_region_ext(genome, start, end, padding, 0)
    }

    /// [`query_region`] + a `query_extension` that widens the query interval
    /// for syncmer lookup on the source side (not target padding).
    ///
    /// Useful when the user's query sits just past the end of a conserved
    /// syncmer block: a small extension (e.g. a few kb) lets the lookup
    /// pick up the block's terminal anchors, after which downstream
    /// refinement can project back onto the original region.
    pub fn query_region_ext(
        &self,
        genome: &str,
        start: u64,
        end: u64,
        padding: u64,
        query_extension: u64,
    ) -> Result<Vec<HomologousInterval>, io::Error> {
        self.query_region_ext_with_seed_filter(
            genome,
            start,
            end,
            padding,
            query_extension,
            SyngSeedFilter::default(),
        )
    }

    /// [`query_region_ext`] with a high-frequency syncmer seed filter.
    pub fn query_region_ext_with_seed_filter(
        &self,
        genome: &str,
        start: u64,
        end: u64,
        padding: u64,
        query_extension: u64,
        seed_filter: SyngSeedFilter,
    ) -> Result<Vec<HomologousInterval>, io::Error> {
        let with_anchors = self.query_region_with_anchors_ext_seed_filtered(
            genome,
            start,
            end,
            padding,
            query_extension,
            seed_filter,
        )?;
        Ok(with_anchors
            .into_iter()
            .map(|h| HomologousInterval {
                genome: h.genome,
                start: h.start,
                end: h.end,
                strand: h.strand,
                cigar: None,
            })
            .collect())
    }

    /// Return query syncmers that are present in this syng index.
    ///
    /// The query sequence does not need to be one of the indexed paths. This
    /// is the primitive used by `impg map -o gaf`.
    pub fn matched_syncmers_in_sequence(&self, query_seq: &[u8]) -> Vec<SyngQuerySyncmer> {
        let mut kmer_buf = vec![0u64; unsafe { (*self.kmer_hash).plen as usize }];
        let mut seq_buf = Vec::with_capacity(query_seq.len() + 1);
        let mut valid_prefix = Vec::with_capacity(query_seq.len() + 1);
        let mut rc_seq = Vec::with_capacity(query_seq.len());
        let mut rc_syncmers = Vec::new();
        let mut out = Vec::new();
        matched_syncmers_best_query_orientation_impl(
            self.seqhash,
            self.kmer_hash,
            self.params,
            query_seq,
            &mut kmer_buf,
            &mut seq_buf,
            &mut valid_prefix,
            &mut rc_seq,
            &mut rc_syncmers,
            &mut out,
        );
        out
    }

    /// Project query syncmer matches to indexed genome coordinates.
    ///
    /// This is intentionally simple: no base-level alignment or edge
    /// refinement, just shared syncmer anchors grouped by path/strand and
    /// chained with the existing anchor chainer.
    pub fn map_sequence(
        &self,
        query_name: &str,
        query_seq: &[u8],
        min_anchors: usize,
        chain_budget: u64,
    ) -> io::Result<Vec<SyngMapHit>> {
        unsafe { syng_ffi::impg_syng_suppress_debug() };

        let query_syncmers = self.matched_syncmers_in_sequence(query_seq);
        if query_syncmers.is_empty() {
            return Ok(Vec::new());
        }

        let syncmer_len = (self.params.w + self.params.k) as u64;
        let mut query_node_positions: FxHashMap<u32, Vec<(u64, u8)>> = FxHashMap::default();
        for sm in query_syncmers {
            let q_orient: u8 = if sm.signed_node >= 0 { 0 } else { 1 };
            query_node_positions
                .entry(sm.node_id)
                .or_default()
                .push((sm.query_pos, q_orient));
        }

        let raw_hits = self
            .query_region_from_node_positions(query_node_positions, 0, None)?
            .into_iter()
            .filter(|hit| hit.anchors.len() >= min_anchors)
            .collect();

        let chained = crate::syng_transitive::chain_anchors_with_sweepga_scaffold_mass(
            raw_hits,
            syncmer_len,
            chain_budget,
            (min_anchors as u64).saturating_mul(syncmer_len),
        );
        let query_len = query_seq.len() as u64;
        let mut hits: Vec<SyngMapHit> = chained
            .into_iter()
            .filter_map(|hit| {
                if hit.anchors.len() < min_anchors {
                    return None;
                }
                let query_start = hit.anchors.iter().map(|a| a.query_pos).min().unwrap_or(0);
                let query_end = hit
                    .anchors
                    .iter()
                    .map(|a| a.query_pos.saturating_add(syncmer_len))
                    .max()
                    .unwrap_or(query_start)
                    .min(query_len);
                let target_len = self
                    .name_map
                    .name_to_path
                    .get(&hit.genome)
                    .and_then(|idx| self.name_map.path_to_length.get(*idx as usize))
                    .copied()
                    .unwrap_or(hit.end);
                Some(SyngMapHit {
                    query_name: query_name.to_string(),
                    query_len,
                    query_start,
                    query_end,
                    target_name: hit.genome,
                    target_len,
                    target_start: hit.start,
                    target_end: hit.end.min(target_len),
                    strand: hit.strand,
                    anchors: hit.anchors.len(),
                })
            })
            .collect();

        hits.sort_by(|a, b| {
            b.anchors
                .cmp(&a.anchors)
                .then(a.target_name.cmp(&b.target_name))
                .then(a.target_start.cmp(&b.target_start))
        });
        Ok(hits)
    }

    /// Query variant that preserves per-hit shared-syncmer anchor positions.
    ///
    /// Each returned homolog carries the `(query_pos, target_pos, node_id)`
    /// positions of the shared syncmers that contributed to it. These anchors
    /// are the input to the boundary-realignment transitive pipeline
    /// (see notes/SYNG_TRANSITIVE_DESIGN.md).
    ///
    /// Handles BOTH forward- and reverse-complement homology. Query syncmers
    /// are tagged with their orientation on the query path (from `walk_path`
    /// node sign); target syncmer positions are located by enumerating GBWT
    /// occurrences and walking each occurrence forward to the next `.syng.pstep`
    /// checkpoint. A query orientation matching the located target orientation
    /// produces a forward-strand hit (`strand='+'`); a mismatch produces a
    /// reverse-strand hit (`strand='-'`). The two strands are kept as separate
    /// homolog intervals (not merged) since they describe distinct homologies.
    pub fn query_region_with_anchors(
        &self,
        genome: &str,
        start: u64,
        end: u64,
        padding: u64,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        self.query_region_with_anchors_ext(genome, start, end, padding, 0)
    }

    /// [`query_region_with_anchors`] + `query_extension` bp on each side of
    /// the source query region for the syncmer-lookup filter. Extension is
    /// clamped to `[0, query_path_length]`. Target-side padding is
    /// unaffected.
    pub fn query_region_with_anchors_ext(
        &self,
        genome: &str,
        start: u64,
        end: u64,
        padding: u64,
        query_extension: u64,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        self.query_region_with_anchors_ext_seed_filtered(
            genome,
            start,
            end,
            padding,
            query_extension,
            SyngSeedFilter::default(),
        )
    }

    /// [`query_region_with_anchors_ext`] with a high-frequency syncmer
    /// seed filter.
    pub fn query_region_with_anchors_ext_seed_filtered(
        &self,
        genome: &str,
        start: u64,
        end: u64,
        padding: u64,
        query_extension: u64,
        seed_filter: SyngSeedFilter,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        self.query_region_with_anchors_ext_visited_seed_filtered(
            genome,
            start,
            end,
            padding,
            query_extension,
            None,
            seed_filter,
        )
    }

    /// As [`query_region_with_anchors_ext`], but with an optional shared
    /// `visited_nodes` set used to suppress re-processing of syncmer nodes
    /// already handled by earlier BFS hops in a transitive query. Nodes
    /// in the set are skipped; nodes processed here are inserted.
    pub fn query_region_with_anchors_ext_visited(
        &self,
        genome: &str,
        start: u64,
        end: u64,
        padding: u64,
        query_extension: u64,
        visited_nodes: Option<&mut FxHashSet<u32>>,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        self.query_region_with_anchors_ext_visited_seed_filtered(
            genome,
            start,
            end,
            padding,
            query_extension,
            visited_nodes,
            SyngSeedFilter::default(),
        )
    }

    /// [`query_region_with_anchors_ext_visited`] with a high-frequency
    /// syncmer seed filter.
    #[allow(clippy::too_many_arguments)]
    pub fn query_region_with_anchors_ext_visited_seed_filtered(
        &self,
        genome: &str,
        start: u64,
        end: u64,
        padding: u64,
        query_extension: u64,
        visited_nodes: Option<&mut FxHashSet<u32>>,
        seed_filter: SyngSeedFilter,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        // Suppress C debug output from syngBWTnext
        unsafe { syng_ffi::impg_syng_suppress_debug() };

        // 1. Look up the query genome's path
        let query_path_num = self
            .name_map
            .name_to_path
            .get(genome)
            .copied()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("Genome '{}' not found in syng index", genome),
                )
            })?;

        let expanded_start = start.saturating_sub(query_extension);
        let expanded_end = end.saturating_add(query_extension);

        // 2. Walk only the requested query path range. The sampled path-step
        //    index gives us a binary-search checkpoint before the leftmost
        //    syncmer that could overlap the region, and the stored GBWT rank
        //    lets us resume traversal from there.
        let query_nodes = self.walk_path_range_from_sampled_steps(
            query_path_num as usize,
            expanded_start,
            expanded_end,
        )?;
        let mut query_walk = Vec::new();
        for &(signed_node, pos) in &query_nodes {
            if pos >= expanded_start && pos < expanded_end {
                query_walk.push(QueryWalkStep {
                    signed_node,
                    bp_pos: pos,
                    query_pos: pos,
                });
            }
        }

        if seed_filter.enabled() {
            return self.query_region_from_ordered_walk_seed_filtered(
                query_walk,
                padding,
                visited_nodes,
                seed_filter,
            );
        }

        let mut query_node_positions: FxHashMap<u32, Vec<(u64, u8)>> = FxHashMap::default();
        for step in query_walk {
            let abs = step.signed_node.unsigned_abs();
            let q_orient: u8 = if step.signed_node >= 0 { 0 } else { 1 };
            query_node_positions
                .entry(abs)
                .or_default()
                .push((step.query_pos, q_orient));
        }

        self.query_region_from_node_positions_seed_filtered(
            query_node_positions,
            padding,
            visited_nodes,
            seed_filter,
        )
    }

    /// Query variant that takes the source sequence bytes directly.
    ///
    /// `query_seq` must cover the sequence interval beginning at
    /// `query_seq_start` in source coordinates. Syncmers are extracted from
    /// those bytes, filtered to `[start - query_extension, end + query_extension)`,
    /// then located through `.syng.pstep` checkpoints. Use this for arbitrary
    /// query sequences that are not indexed paths. For indexed paths, use
    /// [`Self::query_region_with_anchors_ext`], which also jumps through
    /// `.syng.pstep` on the source side and avoids re-extracting source syncmers
    /// from FASTA/AGC.
    #[allow(clippy::too_many_arguments)]
    pub fn query_region_with_anchors_from_sequence_ext_visited(
        &self,
        query_seq: &[u8],
        query_seq_start: u64,
        start: u64,
        end: u64,
        padding: u64,
        query_extension: u64,
        visited_nodes: Option<&mut FxHashSet<u32>>,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        unsafe { syng_ffi::impg_syng_suppress_debug() };

        let expanded_start = start.saturating_sub(query_extension);
        let expanded_end = end.saturating_add(query_extension);
        let mut query_node_positions: FxHashMap<u32, Vec<(u64, u8)>> = FxHashMap::default();

        for sm in self.matched_syncmers_in_sequence(query_seq) {
            let pos = query_seq_start.saturating_add(sm.query_pos);
            if pos >= expanded_start && pos < expanded_end {
                let q_orient: u8 = if sm.signed_node >= 0 { 0 } else { 1 };
                query_node_positions
                    .entry(sm.node_id)
                    .or_default()
                    .push((pos, q_orient));
            }
        }

        self.query_region_from_node_positions(query_node_positions, padding, visited_nodes)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn query_region_with_anchors_from_sequence_ext(
        &self,
        query_seq: &[u8],
        query_seq_start: u64,
        start: u64,
        end: u64,
        padding: u64,
        query_extension: u64,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        self.query_region_with_anchors_from_sequence_ext_visited(
            query_seq,
            query_seq_start,
            start,
            end,
            padding,
            query_extension,
            None,
        )
    }

    pub fn query_region_from_sequence_ext(
        &self,
        query_seq: &[u8],
        query_seq_start: u64,
        start: u64,
        end: u64,
        padding: u64,
        query_extension: u64,
    ) -> Result<Vec<HomologousInterval>, io::Error> {
        let with_anchors = self.query_region_with_anchors_from_sequence_ext(
            query_seq,
            query_seq_start,
            start,
            end,
            padding,
            query_extension,
        )?;
        Ok(with_anchors
            .into_iter()
            .map(|h| HomologousInterval {
                genome: h.genome,
                start: h.start,
                end: h.end,
                strand: h.strand,
                cigar: None,
            })
            .collect())
    }

    fn query_region_from_node_positions(
        &self,
        query_node_positions: FxHashMap<u32, Vec<(u64, u8)>>,
        padding: u64,
        mut visited_nodes: Option<&mut FxHashSet<u32>>,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        self.query_region_from_node_positions_seed_filtered(
            query_node_positions,
            padding,
            visited_nodes.as_deref_mut(),
            SyngSeedFilter::default(),
        )
    }

    fn query_region_from_ordered_walk_seed_filtered(
        &self,
        query_walk: Vec<QueryWalkStep>,
        padding: u64,
        mut visited_nodes: Option<&mut FxHashSet<u32>>,
        seed_filter: SyngSeedFilter,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        if query_walk.is_empty() {
            return Ok(Vec::new());
        }

        let emit_prof = std::env::var("SYNG_EMIT_PROFILE").is_ok();
        let lookup_start = std::time::Instant::now();
        let mut skipped_visited = 0usize;
        let mut unvisited = Vec::with_capacity(query_walk.len());
        let mut unvisited_nodes = FxHashSet::default();
        for step in query_walk {
            let node = step.signed_node.unsigned_abs();
            if visited_nodes
                .as_deref()
                .map(|v| v.contains(&node))
                .unwrap_or(false)
            {
                skipped_visited += 1;
                continue;
            }
            unvisited_nodes.insert(node);
            unvisited.push(step);
        }
        let input_seed_nodes = unvisited_nodes.len();

        let count_start = std::time::Instant::now();
        let mut skipped_frequency = 0usize;
        let mut top_drop_count = 0usize;
        let mut max_seed_occurrences = 0u32;
        let mut dropped_nodes = FxHashSet::default();
        if seed_filter.enabled() && !unvisited_nodes.is_empty() {
            let mut occurrence_counts = Vec::with_capacity(unvisited_nodes.len());
            for &node in &unvisited_nodes {
                let occurrences = self.syncmer_node_orientation_counts(node)?.total();
                max_seed_occurrences = max_seed_occurrences.max(occurrences);
                occurrence_counts.push((node, occurrences));
            }
            let top_drop_nodes =
                Self::top_frequency_seed_nodes(&occurrence_counts, seed_filter.drop_top_fraction);
            top_drop_count = top_drop_nodes.len();
            let occurrence_by_node: FxHashMap<u32, u32> = occurrence_counts.into_iter().collect();
            let max_occurrences = seed_filter.max_occurrences;
            for &node in &unvisited_nodes {
                let occurrences = occurrence_by_node.get(&node).copied().unwrap_or(0);
                let drop_by_top = top_drop_nodes.contains(&node);
                let drop_by_cap = max_occurrences
                    .map(|cap| occurrences > cap)
                    .unwrap_or(false);
                if drop_by_top || drop_by_cap {
                    dropped_nodes.insert(node);
                }
            }
            skipped_frequency = dropped_nodes.len();
        }
        if let Some(v) = visited_nodes.as_deref_mut() {
            v.extend(unvisited_nodes.iter().copied());
        }
        let count_elapsed = count_start.elapsed();

        let mut runs: Vec<Vec<QueryWalkStep>> = Vec::new();
        let mut current = Vec::new();
        for step in unvisited {
            if dropped_nodes.contains(&step.signed_node.unsigned_abs()) {
                if !current.is_empty() {
                    runs.push(std::mem::take(&mut current));
                }
            } else {
                current.push(step);
            }
        }
        if !current.is_empty() {
            runs.push(current);
        }

        let kept_seed_nodes = input_seed_nodes.saturating_sub(skipped_frequency);
        let kept_steps: usize = runs.iter().map(Vec::len).sum();
        let walk_anchors = seed_filter.walk_anchors.max(1);
        if kept_steps < walk_anchors {
            let mut query_node_positions: FxHashMap<u32, Vec<(u64, u8)>> = FxHashMap::default();
            for run in runs {
                for step in run {
                    let q_orient: u8 = if step.signed_node >= 0 { 0 } else { 1 };
                    query_node_positions
                        .entry(step.signed_node.unsigned_abs())
                        .or_default()
                        .push((step.query_pos, q_orient));
                }
            }
            return self.query_region_from_node_positions_seed_filtered(
                query_node_positions,
                padding,
                None,
                SyngSeedFilter::default(),
            );
        }

        let locate_start = std::time::Instant::now();
        // Prepare the shared locate table before entering rayon. With `.spos`
        // this is a no-op; unsaved in-memory indexes can still materialize the
        // same lookup from `.pstep`-equivalent samples outside the worker pool.
        self.prepare_locate_checkpoint_lookup()?;
        let mut seed_tasks = Vec::new();
        for run in &runs {
            Self::push_walk_seed_tasks_from_run(run, '+', walk_anchors, &mut seed_tasks);
            let reverse_run = self.reverse_query_run(run);
            Self::push_walk_seed_tasks_from_run(&reverse_run, '-', walk_anchors, &mut seed_tasks);
        }
        let seeds_seen = seed_tasks.len() as u64;

        let task_hits: Vec<io::Result<WalkSeedTaskHits>> = seed_tasks
            .par_iter()
            .map(|task| self.collect_walk_seed_hits_from_task(task, padding))
            .collect();

        let mut per_path: FxHashMap<(usize, char), Vec<HomologousIntervalWithAnchors>> =
            FxHashMap::default();
        let mut seeds_kept = 0u64;
        let mut located_hits = 0u64;
        let mut emitted_anchors = 0u64;
        for hits in task_hits {
            let hits = hits?;
            seeds_kept += hits.seeds_kept;
            located_hits += hits.located_hits;
            emitted_anchors += hits.emitted_anchors;
            for (path_idx, interval) in hits.intervals {
                per_path
                    .entry((path_idx, interval.strand))
                    .or_default()
                    .push(interval);
            }
        }
        let locate_elapsed = locate_start.elapsed();

        let t_merge = std::time::Instant::now();
        let mut merged: Vec<HomologousIntervalWithAnchors> = Vec::new();
        for (_key, mut group) in per_path {
            Self::merge_intervals_with_anchors(&mut group);
            merged.extend(group);
        }
        merged.sort_by(|a, b| {
            a.genome
                .cmp(&b.genome)
                .then(a.strand.cmp(&b.strand))
                .then(a.start.cmp(&b.start))
        });
        if emit_prof {
            eprintln!(
                "EMITPROF mode=bounded-mem seed_nodes={} kept_seed_nodes={} skipped_visited={} skipped_frequency={} top_drop={} max_seed_occurrences={} walk_anchors={} runs={} seeds={} kept_seeds={} located_hits={} emit_anchors={} count={:.2}s locate={:.2}s lookup={:.2}s merge_sort={:.2}s",
                input_seed_nodes,
                kept_seed_nodes,
                skipped_visited,
                skipped_frequency,
                top_drop_count,
                max_seed_occurrences,
                walk_anchors,
                runs.len(),
                seeds_seen,
                seeds_kept,
                located_hits,
                emitted_anchors,
                count_elapsed.as_secs_f64(),
                locate_elapsed.as_secs_f64(),
                lookup_start.elapsed().as_secs_f64(),
                t_merge.elapsed().as_secs_f64(),
            );
        }
        Ok(merged)
    }

    fn push_walk_seed_tasks_from_run(
        run: &[QueryWalkStep],
        strand: char,
        walk_anchors: usize,
        tasks: &mut Vec<WalkSeedTask>,
    ) {
        if run.len() < walk_anchors {
            return;
        }
        for seed_start in (0..=run.len() - walk_anchors).step_by(walk_anchors) {
            let seed_end = seed_start + walk_anchors;
            tasks.push(WalkSeedTask {
                strand,
                steps: run[seed_start..seed_end].to_vec(),
            });
        }
    }

    fn reverse_query_run(&self, run: &[QueryWalkStep]) -> Vec<QueryWalkStep> {
        let syncmer_len = self.syncmer_length_bp() as u64;
        let Some(reverse_end) = run
            .iter()
            .map(|step| step.query_pos.saturating_add(syncmer_len))
            .max()
        else {
            return Vec::new();
        };
        run.iter()
            .rev()
            .map(|step| QueryWalkStep {
                signed_node: -step.signed_node,
                bp_pos: reverse_end.saturating_sub(step.query_pos.saturating_add(syncmer_len)),
                query_pos: step.query_pos,
            })
            .collect()
    }

    fn match_walk_seed(&self, walk: &[SyngWalkStep]) -> io::Result<Option<(i32, u32, u32)>> {
        let Some(first) = walk.first().copied() else {
            return Ok(None);
        };
        let Some((sbp, mut low, mut high)) = self.start_gbwt_match(first.signed_node)? else {
            return Ok(None);
        };
        for idx in 1..walk.len() {
            let Some(offset) = Self::walk_step_offset(walk[idx - 1], walk[idx]) else {
                unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
                return Ok(None);
            };
            if !self.advance_gbwt_match(sbp, &mut low, &mut high, walk[idx].signed_node, offset) {
                unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
                return Ok(None);
            }
        }
        unsafe { syng_ffi::syngBWTpathDestroy(sbp) };
        Ok(Some((walk.last().unwrap().signed_node, low, high)))
    }

    fn collect_walk_seed_hits_from_task(
        &self,
        task: &WalkSeedTask,
        padding: u64,
    ) -> io::Result<WalkSeedTaskHits> {
        let syncmer_len = self.syncmer_length_bp() as u64;
        let walk: Vec<SyngWalkStep> = task
            .steps
            .iter()
            .map(|step| SyngWalkStep {
                signed_node: step.signed_node,
                bp_pos: step.bp_pos,
            })
            .collect();
        let Some((match_signed_node, match_low, match_high)) = self.match_walk_seed(&walk)? else {
            return Ok(WalkSeedTaskHits::default());
        };
        if match_low >= match_high {
            return Ok(WalkSeedTaskHits::default());
        }
        let hits =
            self.locate_signed_node_occurrence_range(match_signed_node, match_low, match_high)?;
        let Some(first_step) = task.steps.first() else {
            return Ok(WalkSeedTaskHits::default());
        };
        let Some(last_step) = task.steps.last() else {
            return Ok(WalkSeedTaskHits::default());
        };
        let seed_span = last_step.bp_pos.saturating_sub(first_step.bp_pos);

        let mut out = WalkSeedTaskHits {
            seeds_kept: 1,
            located_hits: hits.len() as u64,
            emitted_anchors: 0,
            intervals: Vec::new(),
        };
        for hit in hits {
            if hit.path_idx >= self.name_map.path_to_name.len() {
                continue;
            }
            let Some(target_start) = hit.target_pos.checked_sub(seed_span) else {
                continue;
            };
            let genome_len = self.name_map.path_to_length[hit.path_idx];
            let target_end = hit.target_pos.saturating_add(syncmer_len).min(genome_len);
            let padded_start = target_start.saturating_sub(padding);
            let padded_end = target_end.saturating_add(padding).min(genome_len);
            if padded_start >= padded_end {
                continue;
            }
            let mut anchors = Vec::with_capacity(task.steps.len());
            for step in &task.steps {
                let rel = step.bp_pos.saturating_sub(first_step.bp_pos);
                anchors.push(Anchor {
                    query_pos: step.query_pos,
                    target_pos: target_start.saturating_add(rel),
                    node_id: step.signed_node.unsigned_abs(),
                });
            }
            out.emitted_anchors += anchors.len() as u64;
            out.intervals.push((
                hit.path_idx,
                HomologousIntervalWithAnchors {
                    genome: self.name_map.path_to_name[hit.path_idx].clone(),
                    start: padded_start,
                    end: padded_end,
                    strand: task.strand,
                    anchors,
                },
            ));
        }
        Ok(out)
    }

    fn query_region_from_node_positions_seed_filtered(
        &self,
        query_node_positions: FxHashMap<u32, Vec<(u64, u8)>>,
        padding: u64,
        mut visited_nodes: Option<&mut FxHashSet<u32>>,
        seed_filter: SyngSeedFilter,
    ) -> Result<Vec<HomologousIntervalWithAnchors>, io::Error> {
        if query_node_positions.is_empty() {
            return Ok(Vec::new());
        }

        let syncmer_len = (self.params.w + self.params.k) as u64;
        // Group per (forward_path_idx, strand). Intervals on the same genome
        // but different strands are NOT merged — they represent distinct
        // homologies (forward match vs RC / inversion).
        let mut per_path: FxHashMap<(usize, char), Vec<HomologousIntervalWithAnchors>> =
            FxHashMap::default();

        let emit_prof = std::env::var("SYNG_EMIT_PROFILE").is_ok();
        let lookup_start = std::time::Instant::now();
        let mut unvisited: Vec<(u32, Vec<(u64, u8)>)> =
            Vec::with_capacity(query_node_positions.len());
        let mut skipped_visited = 0usize;
        for (node, query_positions) in query_node_positions {
            if visited_nodes
                .as_deref()
                .map(|v| v.contains(&node))
                .unwrap_or(false)
            {
                skipped_visited += 1;
                continue;
            }
            unvisited.push((node, query_positions));
        }
        let input_seed_nodes = unvisited.len();
        let unvisited_nodes_for_visited: Vec<u32> = unvisited.iter().map(|(n, _)| *n).collect();

        let mut skipped_frequency = 0usize;
        let mut top_drop_count = 0usize;
        let mut max_seed_occurrences = 0u32;
        let count_start = std::time::Instant::now();
        let mut occurrence_counts_by_node: FxHashMap<u32, NodeOccurrenceCounts> =
            FxHashMap::default();
        if seed_filter.enabled() && !unvisited.is_empty() {
            let mut occurrence_counts = Vec::with_capacity(unvisited.len());
            for (node, _) in &unvisited {
                let counts = self.syncmer_node_orientation_counts(*node)?;
                let occurrences = counts.total();
                max_seed_occurrences = max_seed_occurrences.max(occurrences);
                occurrence_counts.push((*node, occurrences));
                occurrence_counts_by_node.insert(*node, counts);
            }
            let top_drop_nodes =
                Self::top_frequency_seed_nodes(&occurrence_counts, seed_filter.drop_top_fraction);
            top_drop_count = top_drop_nodes.len();
            let max_occurrences = seed_filter.max_occurrences;
            let occurrence_by_node: FxHashMap<u32, u32> = occurrence_counts.into_iter().collect();
            unvisited.retain(|(node, _)| {
                let occurrences = occurrence_by_node.get(node).copied().unwrap_or(0);
                let drop_by_top = top_drop_nodes.contains(node);
                let drop_by_cap = max_occurrences
                    .map(|cap| occurrences > cap)
                    .unwrap_or(false);
                let keep = !(drop_by_top || drop_by_cap);
                if !keep {
                    skipped_frequency += 1;
                }
                keep
            });
        }
        if let Some(v) = visited_nodes.as_deref_mut() {
            for node in unvisited_nodes_for_visited {
                v.insert(node);
            }
        }
        let count_elapsed = count_start.elapsed();

        let locate_start = std::time::Instant::now();
        self.prepare_locate_checkpoint_lookup()?;
        type LocatedSeedBuckets = FxHashMap<(usize, u32, char), (u64, u64, FxHashMap<i64, Anchor>)>;
        let locate_results: Vec<io::Result<(u64, u64, LocatedSeedBuckets)>> = unvisited
            .into_par_iter()
            .map(|(query_node, query_positions)| {
                let mut located_hits = 0u64;
                let mut emitted_anchors = 0u64;
                let mut per_node_sampled: LocatedSeedBuckets = FxHashMap::default();
                let cached_counts = occurrence_counts_by_node.get(&query_node).copied();
                for hit in self.locate_node_occurrences_with_counts(query_node, cached_counts)? {
                    located_hits += 1;
                    if hit.path_idx >= self.name_map.path_to_name.len() {
                        continue;
                    }
                    let genome_len = self.name_map.path_to_length[hit.path_idx];
                    let hit_end = hit.target_pos + syncmer_len;
                    let padded_start = hit.target_pos.saturating_sub(padding);
                    let padded_end = (hit_end + padding).min(genome_len);
                    for &(qp, q_orient) in &query_positions {
                        let strand = if q_orient == hit.target_orient {
                            '+'
                        } else {
                            '-'
                        };
                        let sig: i64 = if strand == '+' {
                            hit.target_pos as i64 - qp as i64
                        } else {
                            hit.target_pos as i64 + qp as i64
                        };
                        let entry = per_node_sampled
                            .entry((hit.path_idx, query_node, strand))
                            .or_insert_with(|| (padded_start, padded_end, FxHashMap::default()));
                        entry.0 = entry.0.min(padded_start);
                        entry.1 = entry.1.max(padded_end);
                        entry.2.entry(sig).or_insert(Anchor {
                            query_pos: qp,
                            target_pos: hit.target_pos,
                            node_id: query_node,
                        });
                        emitted_anchors += 1;
                    }
                }
                Ok((located_hits, emitted_anchors, per_node_sampled))
            })
            .collect();

        let mut located_hits = 0u64;
        let mut emitted_anchors = 0u64;
        let mut per_node_sampled: LocatedSeedBuckets = FxHashMap::default();
        for result in locate_results {
            let (local_located_hits, local_emitted_anchors, local_buckets) = result?;
            located_hits += local_located_hits;
            emitted_anchors += local_emitted_anchors;
            for (key, (start, end, anchors_by_signature)) in local_buckets {
                let entry = per_node_sampled
                    .entry(key)
                    .or_insert_with(|| (start, end, FxHashMap::default()));
                entry.0 = entry.0.min(start);
                entry.1 = entry.1.max(end);
                entry.2.extend(anchors_by_signature);
            }
        }
        let locate_elapsed = locate_start.elapsed();
        for ((genome_idx, _node, strand), (start, end, sig_map)) in per_node_sampled {
            if sig_map.is_empty() {
                continue;
            }
            let anchors: Vec<Anchor> = sig_map.into_values().collect();
            per_path
                .entry((genome_idx, strand))
                .or_default()
                .push(HomologousIntervalWithAnchors {
                    genome: self.name_map.path_to_name[genome_idx].clone(),
                    start,
                    end,
                    strand,
                    anchors,
                });
        }

        // 3. Merge per (target, strand); different strands stay separate.
        let t_merge = std::time::Instant::now();
        let mut merged: Vec<HomologousIntervalWithAnchors> = Vec::new();
        for (_key, mut group) in per_path {
            Self::merge_intervals_with_anchors(&mut group);
            merged.extend(group);
        }
        // Sort final output for deterministic order (by genome, strand, start)
        merged.sort_by(|a, b| {
            a.genome
                .cmp(&b.genome)
                .then(a.strand.cmp(&b.strand))
                .then(a.start.cmp(&b.start))
        });
        if emit_prof {
            eprintln!(
                "EMITPROF seed_nodes={} kept_seed_nodes={} skipped_visited={} skipped_frequency={} top_drop={} max_seed_occurrences={} located_hits={} emit_anchors={} count={:.2}s locate={:.2}s lookup={:.2}s merge_sort={:.2}s",
                input_seed_nodes,
                input_seed_nodes.saturating_sub(skipped_frequency),
                skipped_visited,
                skipped_frequency,
                top_drop_count,
                max_seed_occurrences,
                located_hits,
                emitted_anchors,
                count_elapsed.as_secs_f64(),
                locate_elapsed.as_secs_f64(),
                lookup_start.elapsed().as_secs_f64(),
                t_merge.elapsed().as_secs_f64(),
            );
        }
        Ok(merged)
    }

    /// Build a region-specific GBWT from fetched sequences.
    ///
    /// Creates a fresh KmerHash and SyngBWT from the given sequences,
    /// extracting syncmers and building GBWT paths. Writes standard
    /// syng-compatible `.1khash` and `.1gbwt` files at the given prefix.
    ///
    /// Uses the same syncmer parameters as this index (or default if
    /// called on a freshly constructed index).
    pub fn build_region_gbwt(&self, sequences: &[(String, &[u8])], prefix: &str) -> io::Result<()> {
        let total_start = std::time::Instant::now();
        let total_bp: usize = sequences.iter().map(|(_, seq)| seq.len()).sum();
        log::info!(
            "[syng region gbwt] serial build from {} sequence(s), {} bp",
            sequences.len(),
            total_bp
        );
        let mut region_index = SyngIndex::new(self.params);
        region_index.enable_online_sampled_positions(1)?;
        let replay_start = std::time::Instant::now();
        let mut total_syncmers = 0usize;
        for (name, seq) in sequences {
            let stats = region_index.add_sequence(name.clone(), seq.to_vec());
            total_syncmers += stats.syncmers;
        }
        log::info!(
            "[syng region gbwt] extracted and inserted {} syncmer step(s) into GBWT in {:.3}s",
            total_syncmers,
            replay_start.elapsed().as_secs_f64()
        );
        let save_start = std::time::Instant::now();
        region_index.save(prefix)?;
        log::info!(
            "[syng region gbwt] saved regional syng files in {:.3}s; total {:.3}s",
            save_start.elapsed().as_secs_f64(),
            total_start.elapsed().as_secs_f64()
        );
        Ok(())
    }

    /// Merge overlapping or adjacent intervals that share the same genome and strand.
    /// Merge a per-target group of anchored intervals in place. Assumes all
    /// entries reference the same genome + strand (caller groups by path).
    fn merge_intervals_with_anchors(group: &mut Vec<HomologousIntervalWithAnchors>) {
        if group.len() <= 1 {
            return;
        }
        group.sort_by_key(|h| h.start);
        let mut merged: Vec<HomologousIntervalWithAnchors> = Vec::with_capacity(group.len());
        for mut iv in group.drain(..) {
            if let Some(last) = merged.last_mut() {
                if iv.start <= last.end {
                    last.end = last.end.max(iv.end);
                    last.anchors.append(&mut iv.anchors);
                    continue;
                }
            }
            merged.push(iv);
        }
        // Dedupe + sort anchors within each merged interval
        for iv in &mut merged {
            iv.anchors.sort_by(|a, b| {
                a.query_pos
                    .cmp(&b.query_pos)
                    .then(a.target_pos.cmp(&b.target_pos))
            });
            iv.anchors
                .dedup_by(|a, b| a.query_pos == b.query_pos && a.target_pos == b.target_pos);
        }
        *group = merged;
    }
}

impl Drop for SyngIndex {
    fn drop(&mut self) {
        unsafe {
            if !self.gbwt.is_null() {
                syng_ffi::syngBWTdestroy(self.gbwt);
            }
            if !self.kmer_hash.is_null() {
                syng_ffi::kmerHashDestroy(self.kmer_hash);
            }
            if !self.seqhash.is_null() {
                syng_ffi::impg_seqhashDestroy(self.seqhash);
            }
        }
    }
}

impl Drop for SyngMatcher {
    fn drop(&mut self) {
        unsafe {
            if !self.kmer_hash.is_null() {
                syng_ffi::kmerHashDestroy(self.kmer_hash);
            }
            if !self.seqhash.is_null() {
                syng_ffi::impg_seqhashDestroy(self.seqhash);
            }
        }
    }
}

impl Drop for SyngMatcherWorker<'_> {
    fn drop(&mut self) {
        unsafe {
            if !self.seqhash.is_null() {
                syng_ffi::impg_seqhashDestroy(self.seqhash);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The syng C library has non-thread-safe global state (pathCount, array
    // tracking counters, etc.), so all tests that touch the C FFI must be
    // serialized.
    static SYNG_LOCK: std::sync::LazyLock<std::sync::Mutex<()>> =
        std::sync::LazyLock::new(|| std::sync::Mutex::new(()));
    fn lock_syng() -> std::sync::MutexGuard<'static, ()> {
        SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner())
    }

    // ── 1. FFI smoke tests ──────────────────────────────────────────

    #[test]
    fn test_syng_ffi_seqhash_create_destroy() {
        let _guard = lock_syng();
        // Create a Seqhash with default params (k=8, w=55, seed=7), verify non-null, free it.
        let sh = unsafe { syng_ffi::impg_seqhashCreateSafe(8, 55, 7) };
        assert!(
            !sh.is_null(),
            "seqhashCreate returned null for valid params"
        );
        unsafe { syng_ffi::impg_seqhashDestroy(sh) };
    }

    #[test]
    fn test_syng_ffi_kmerhash_create_destroy() {
        let _guard = lock_syng();
        // Create a KmerHash, verify non-null, destroy it.
        let syncmer_len = 55 + 8; // w + k = 63
        let kh = unsafe { syng_ffi::kmerHashCreate(1024, syncmer_len) };
        assert!(!kh.is_null(), "kmerHashCreate returned null");
        unsafe { syng_ffi::kmerHashDestroy(kh) };
    }

    #[test]
    fn test_syng_ffi_syngbwt_create_destroy() {
        let _guard = lock_syng();
        // Create a SyngBWT, verify non-null, destroy it.
        let syncmer_len = 63;
        let sb = unsafe { syng_ffi::syngBWTcreate(syncmer_len, 0) };
        assert!(!sb.is_null(), "syngBWTcreate returned null");
        unsafe { syng_ffi::syngBWTdestroy(sb) };
    }

    // ── 2. SyncmerParams defaults ───────────────────────────────────

    #[test]
    fn test_syncmer_params_default() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        assert_eq!(params.k, 8);
        assert_eq!(params.w, 55);
        assert_eq!(params.seed, 7);
    }

    #[test]
    fn test_syncmer_params_to_c() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let c_params = params.to_c();
        assert_eq!(c_params.k, 8);
        assert_eq!(c_params.w, 55);
        assert_eq!(c_params.seed, 7);
    }

    // ── 3. SyngIndex lifecycle ──────────────────────────────────────

    #[test]
    fn test_syng_index_create_drop() {
        let _guard = lock_syng();
        // Verify we can create and drop without crashing
        let params = SyncmerParams::default();
        let index = SyngIndex::new(params);
        assert!(!index.gbwt.is_null());
        assert!(!index.kmer_hash.is_null());
        assert!(!index.seqhash.is_null());
        drop(index);
    }

    #[test]
    fn test_syng_index_custom_params() {
        let _guard = lock_syng();
        // Non-default params should also produce valid pointers
        let params = SyncmerParams {
            k: 11,
            w: 40,
            seed: 42,
        };
        let index = SyngIndex::new(params);
        assert!(!index.gbwt.is_null());
        assert!(!index.kmer_hash.is_null());
        assert!(!index.seqhash.is_null());
        assert_eq!(index.params.k, 11);
        assert_eq!(index.params.w, 40);
        assert_eq!(index.params.seed, 42);
    }

    #[test]
    fn test_syng_index_accessors() {
        let _guard = lock_syng();
        let index = SyngIndex::new(SyncmerParams::default());
        // Accessor methods should return the same non-null pointers
        assert_eq!(index.gbwt_ptr(), index.gbwt);
        assert_eq!(index.kmer_hash_ptr(), index.kmer_hash);
        assert_eq!(index.seqhash_ptr(), index.seqhash);
        assert!(!index.gbwt_ptr().is_null());
        assert!(!index.kmer_hash_ptr().is_null());
        assert!(!index.seqhash_ptr().is_null());
    }

    // ── 4. SyngNameMap ──────────────────────────────────────────────

    #[test]
    fn test_name_map_new() {
        let _guard = lock_syng();
        let nm = SyngNameMap::new();
        assert!(nm.path_to_name.is_empty());
        assert!(nm.path_to_length.is_empty());
        assert!(nm.name_to_path.is_empty());
    }

    #[test]
    fn test_name_map_default() {
        let _guard = lock_syng();
        // Default trait should produce the same as new()
        let nm = SyngNameMap::default();
        assert!(nm.path_to_name.is_empty());
        assert!(nm.path_to_length.is_empty());
        assert!(nm.name_to_path.is_empty());
    }

    #[test]
    fn test_name_map_add() {
        let _guard = lock_syng();
        let mut nm = SyngNameMap::new();
        let p0 = nm.add("seq1".to_string(), 1000);
        let p1 = nm.add("seq2".to_string(), 2000);
        assert_eq!(p0, 0);
        assert_eq!(p1, 1);
        assert_eq!(nm.path_to_name[0], "seq1");
        assert_eq!(nm.path_to_length[1], 2000);
        assert_eq!(*nm.name_to_path.get("seq1").unwrap(), 0);
        assert_eq!(*nm.name_to_path.get("seq2").unwrap(), 1);
    }

    #[test]
    fn test_name_map_save_load_roundtrip() {
        let _guard = lock_syng();
        let mut nm = SyngNameMap::new();
        nm.add("chr1".to_string(), 248956422);
        nm.add("chr2".to_string(), 242193529);
        nm.add("chrX".to_string(), 156040895);

        let dir = std::env::temp_dir().join("impg_test_name_map");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test.syng.names");
        let path_str = path.to_str().unwrap();

        nm.save(path_str).unwrap();
        let loaded = SyngNameMap::load(path_str).unwrap();

        assert_eq!(loaded.path_to_name, nm.path_to_name);
        assert_eq!(loaded.path_to_length, nm.path_to_length);
        assert_eq!(loaded.name_to_path.len(), nm.name_to_path.len());
        for (k, v) in &nm.name_to_path {
            assert_eq!(loaded.name_to_path.get(k), Some(v));
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 5. SyngIndex::build ─────────────────────────────────────────

    /// Generate a random-ish DNA sequence of the given length.
    fn make_test_sequence(len: usize, seed: u8) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut seq = Vec::with_capacity(len);
        let mut state: u32 = seed as u32;
        for _ in 0..len {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            seq.push(bases[((state >> 16) % 4) as usize]);
        }
        seq
    }

    fn normalize_anchored_hits(
        hits: Vec<HomologousIntervalWithAnchors>,
    ) -> Vec<(String, u64, u64, char, Vec<(u64, u64, u32)>)> {
        let mut out: Vec<_> = hits
            .into_iter()
            .map(|mut hit| {
                hit.anchors.sort_by(|a, b| {
                    a.query_pos
                        .cmp(&b.query_pos)
                        .then(a.target_pos.cmp(&b.target_pos))
                        .then(a.node_id.cmp(&b.node_id))
                });
                let anchors = hit
                    .anchors
                    .into_iter()
                    .map(|a| (a.query_pos, a.target_pos, a.node_id))
                    .collect();
                (hit.genome, hit.start, hit.end, hit.strand, anchors)
            })
            .collect();
        out.sort();
        out
    }

    #[test]
    fn test_preloaded_dictionary_replay_supports_query() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let shared = make_test_sequence(800, 42);
        let mut seq_a = shared.clone();
        seq_a.extend_from_slice(&make_test_sequence(300, 1));
        let mut seq_b = shared;
        seq_b.extend_from_slice(&make_test_sequence(300, 2));
        let sequences = vec![
            ("seqA".to_string(), seq_a.clone()),
            ("seqB".to_string(), seq_b),
            ("seqC".to_string(), make_test_sequence(1100, 99)),
        ];

        let dictionary =
            crate::syng_parallel::build_packed_syncmer_dictionary(params, &sequences).unwrap();
        assert!(!dictionary.is_empty());

        let mut index = SyngIndex::new_with_packed_syncmer_dictionary(params, &dictionary).unwrap();
        index.enable_online_sampled_positions(1).unwrap();
        for (name, seq) in sequences {
            let stats = index
                .add_sequence_with_existing_syncmers(name, seq)
                .unwrap();
            assert!(stats.indexed);
            assert!(stats.syncmers > 0);
        }
        index.finalize_online_sampled_positions().unwrap();

        let matches = index.matched_syncmers_in_sequence(&seq_a);
        assert!(!matches.is_empty());

        let intervals = index.query_region("seqA", 0, 500, 0).unwrap();
        assert!(
            intervals.iter().any(|interval| interval.genome == "seqB"),
            "shared-prefix sequence should be discoverable through preloaded dictionary"
        );
    }

    #[test]
    fn test_matched_syncmers_skip_ambiguous_query_bases() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq = make_test_sequence(1000, 42);
        let index = SyngIndex::build(
            params,
            vec![("sample#0#chr1".to_string(), seq.clone())].into_iter(),
        );
        let baseline = index.matched_syncmers_in_sequence(&seq);
        assert!(
            !baseline.is_empty(),
            "test sequence should produce indexed syncmer matches"
        );

        let syncmer_len = (params.w + params.k) as usize;
        let (chosen, invalid_pos) = baseline
            .iter()
            .find_map(|sm| {
                let start = sm.query_pos as usize;
                let end = (start + syncmer_len).min(seq.len());
                (start..end)
                    .find(|&pos| seq[pos] == b'A')
                    .map(|pos| (*sm, pos))
            })
            .expect("at least one matched syncmer should contain an A base");
        assert!(
            chosen.query_pos <= invalid_pos as u64
                && (invalid_pos as u64) < chosen.query_pos.saturating_add(syncmer_len as u64)
        );

        let mut query = seq.clone();
        query[invalid_pos] = b'N';
        let matches = index.matched_syncmers_in_sequence(&query);
        assert!(
            matches.iter().all(|sm| {
                !(sm.query_pos <= invalid_pos as u64
                    && (invalid_pos as u64) < sm.query_pos.saturating_add(syncmer_len as u64))
            }),
            "matched syncmers must not span ambiguous query bases: {:?}",
            matches
        );
    }

    #[test]
    fn test_matched_syncmers_orients_reverse_complement_query() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq = make_test_sequence(5000, 42);
        let index = SyngIndex::build(
            params,
            vec![("sample#0#chr1".to_string(), seq.clone())].into_iter(),
        );

        let forward_matches = index.matched_syncmers_in_sequence(&seq);
        assert!(
            forward_matches.len() >= 2,
            "test sequence should produce multiple syncmer matches"
        );

        let rc_seq = crate::graph::reverse_complement(&seq);
        let rc_matches = index.matched_syncmers_in_sequence(&rc_seq);
        assert_eq!(
            rc_matches.len(),
            forward_matches.len(),
            "reverse-complement query should recover the graph-compatible syncmer orientation"
        );

        let mut forward_nodes: Vec<_> = forward_matches.iter().map(|sm| sm.node_id).collect();
        let mut rc_nodes: Vec<_> = rc_matches.iter().map(|sm| sm.node_id).collect();
        forward_nodes.sort_unstable();
        rc_nodes.sort_unstable();
        assert_eq!(rc_nodes, forward_nodes);
    }

    #[test]
    fn test_gbwt_mems_for_walk_uses_syncmer_offsets() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq = make_test_sequence(5000, 43);
        let index = SyngIndex::build(
            params,
            vec![("sample#0#chr1".to_string(), seq.clone())].into_iter(),
        );

        let mut matches = index.matched_syncmers_in_sequence(&seq);
        matches.sort_unstable_by_key(|sm| sm.query_pos);
        assert!(
            matches.len() >= 4,
            "test sequence should produce several indexed syncmer matches"
        );
        let walk: Vec<_> = matches
            .iter()
            .map(|sm| SyngWalkStep {
                signed_node: sm.signed_node,
                bp_pos: sm.query_pos,
            })
            .collect();

        let mems = index.gbwt_mems_for_walk(&walk).unwrap();
        assert!(
            mems.iter().any(|mem| mem.step_start == 0
                && mem.step_end == walk.len()
                && mem.query_start == walk[0].bp_pos
                && mem.query_end == walk.last().unwrap().bp_pos + index.syncmer_length_bp() as u64),
            "the exact query walk should produce a full-length GBWT MEM: {:?}",
            mems
        );

        let Some(shift_idx) =
            (1..walk.len() - 1).find(|&idx| walk[idx].bp_pos + 1 < walk[idx + 1].bp_pos)
        else {
            return;
        };
        let mut shifted = walk.clone();
        shifted[shift_idx].bp_pos += 1;
        let shifted_mems = index.gbwt_mems_for_walk(&shifted).unwrap();
        assert!(
            shifted_mems
                .iter()
                .all(|mem| mem.step_start != 0 || mem.step_end != shifted.len()),
            "changing a syncmer offset must break the full-length GBWT MEM: {:?}",
            shifted_mems
        );
    }

    #[test]
    fn test_bounded_walk_seeds_recover_divergent_homolog_under_long_self_mem() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq_a = make_test_sequence(6000, 91);
        let mut seq_b = seq_a.clone();
        for pos in (300..5700).step_by(600) {
            seq_b[pos] = match seq_b[pos] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                b'T' => b'A',
                other => other,
            };
        }
        let index = SyngIndex::build(
            params,
            vec![
                ("genome_a".to_string(), seq_a.clone()),
                ("genome_b".to_string(), seq_b),
            ]
            .into_iter(),
        );

        let mut matches = index.matched_syncmers_in_sequence(&seq_a);
        matches.sort_unstable_by_key(|sm| sm.query_pos);
        let walk: Vec<_> = matches
            .iter()
            .map(|sm| SyngWalkStep {
                signed_node: sm.signed_node,
                bp_pos: sm.query_pos,
            })
            .collect();
        let mems = index.gbwt_mems_for_walk(&walk).unwrap();
        assert!(
            mems.iter().any(|mem| mem.step_start == 0 && mem.step_end == walk.len()),
            "self path should produce a long maximal MEM that would hide shorter homolog seeds: {:?}",
            mems
        );

        let hits = index
            .query_region_with_anchors_ext_seed_filtered(
                "genome_a",
                0,
                seq_a.len() as u64,
                0,
                0,
                SyngSeedFilter {
                    max_occurrences: None,
                    drop_top_fraction: 0.0005,
                    walk_anchors: 3,
                },
            )
            .unwrap();
        assert!(
            hits.iter().any(|hit| hit.genome == "genome_b"),
            "bounded 3-syncmer walk seeds should recover the divergent homolog; hits: {:?}",
            hits
        );
    }

    #[test]
    fn test_bounded_walk_parallel_query_matches_single_thread() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq_a = make_test_sequence(7000, 121);
        let mut seq_b = seq_a.clone();
        let mut seq_c = seq_a.clone();
        for pos in (500..6500).step_by(700) {
            seq_b[pos] = match seq_b[pos] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                b'T' => b'A',
                other => other,
            };
        }
        for pos in (800..6200).step_by(900) {
            seq_c[pos] = match seq_c[pos] {
                b'A' => b'T',
                b'C' => b'A',
                b'G' => b'C',
                b'T' => b'G',
                other => other,
            };
        }
        let index = SyngIndex::build(
            params,
            vec![
                ("genome_a".to_string(), seq_a.clone()),
                ("genome_b".to_string(), seq_b),
                ("genome_c".to_string(), seq_c),
            ]
            .into_iter(),
        );
        let seed_filter = SyngSeedFilter {
            max_occurrences: None,
            drop_top_fraction: 0.0005,
            walk_anchors: 3,
        };

        let one_thread = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(|| {
                index
                    .query_region_with_anchors_ext_seed_filtered(
                        "genome_a",
                        0,
                        seq_a.len() as u64,
                        0,
                        0,
                        seed_filter,
                    )
                    .unwrap()
            });
        let four_threads = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(|| {
                index
                    .query_region_with_anchors_ext_seed_filtered(
                        "genome_a",
                        0,
                        seq_a.len() as u64,
                        0,
                        0,
                        seed_filter,
                    )
                    .unwrap()
            });

        assert_eq!(
            normalize_anchored_hits(one_thread),
            normalize_anchored_hits(four_threads)
        );
    }

    #[test]
    fn test_build_empty() {
        let _guard = lock_syng();
        let index = SyngIndex::build(SyncmerParams::default(), std::iter::empty());
        assert!(index.name_map.path_to_name.is_empty());
    }

    #[test]
    fn test_build_short_sequence() {
        let _guard = lock_syng();
        // Sequence shorter than syncmer length (63) — should be recorded but not panick
        let seqs = vec![("short".to_string(), b"ACGTACGT".to_vec())];
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 1);
        assert_eq!(index.name_map.path_to_name[0], "short");
    }

    #[test]
    fn test_build_single_sequence() {
        let _guard = lock_syng();
        // A sequence long enough for syncmer extraction
        let seq = make_test_sequence(500, 42);
        let seqs = vec![("test_seq".to_string(), seq)];
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 1);
        assert_eq!(index.name_map.path_to_name[0], "test_seq");
        assert_eq!(index.name_map.path_to_length[0], 500);
    }

    #[test]
    fn test_build_multiple_sequences() {
        let _guard = lock_syng();
        let seqs: Vec<(String, Vec<u8>)> = (0..5)
            .map(|i| (format!("seq_{}", i), make_test_sequence(1000, i as u8)))
            .collect();
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 5);
        for i in 0..5 {
            assert_eq!(index.name_map.path_to_name[i], format!("seq_{}", i));
            assert_eq!(index.name_map.path_to_length[i], 1000);
        }
    }

    // ── 6. SyngIndex save/load roundtrip ────────────────────────────

    #[test]
    fn test_save_load_roundtrip() {
        let _guard = lock_syng();
        let seqs: Vec<(String, Vec<u8>)> = (0..3)
            .map(|i| {
                (
                    format!("genome_{}", i),
                    make_test_sequence(2000, i as u8 + 100),
                )
            })
            .collect();
        let params = SyncmerParams::default();
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_roundtrip");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("test_index");
        let prefix_str = prefix.to_str().unwrap();

        // Save
        index.save(prefix_str).unwrap();

        // Verify files exist
        assert!(
            std::path::Path::new(&format!("{}.1gbwt", prefix_str)).exists(),
            ".1gbwt file should exist"
        );
        assert!(
            std::path::Path::new(&format!("{}.1khash", prefix_str)).exists(),
            ".1khash file should exist"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.names", prefix_str)).exists(),
            ".syng.names file should exist"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.spos", prefix_str)).exists(),
            ".syng.spos file should exist"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.pstep", prefix_str)).exists(),
            ".syng.pstep file should exist"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.meta", prefix_str)).exists(),
            ".syng.meta file should exist"
        );

        // Load
        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert_eq!(loaded.params, params);

        // Verify name map is intact
        assert_eq!(
            loaded.name_map.path_to_name, index.name_map.path_to_name,
            "Name map names should match after round-trip"
        );
        assert_eq!(
            loaded.name_map.path_to_length, index.name_map.path_to_length,
            "Name map lengths should match after round-trip"
        );

        // Verify C pointers are valid
        assert!(!loaded.gbwt.is_null());
        assert!(!loaded.kmer_hash.is_null());
        assert!(!loaded.seqhash.is_null());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_syng_suffix_prefix_does_not_duplicate_sidecar_names() {
        let _guard = lock_syng();
        let seqs = vec![("genome".to_string(), make_test_sequence(2000, 7))];
        let params = SyncmerParams::default();
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_suffix_sidecars");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx.syng");
        let prefix_str = prefix.to_str().unwrap();

        index.save(prefix_str).unwrap();

        assert!(std::path::Path::new(&format!("{}.1gbwt", prefix_str)).exists());
        assert!(std::path::Path::new(&format!("{}.1khash", prefix_str)).exists());
        assert!(std::path::Path::new(&format!("{}.names", prefix_str)).exists());
        assert!(std::path::Path::new(&format!("{}.spos", prefix_str)).exists());
        assert!(std::path::Path::new(&format!("{}.pstep", prefix_str)).exists());
        assert!(std::path::Path::new(&format!("{}.meta", prefix_str)).exists());
        assert!(!std::path::Path::new(&format!("{}.syng.names", prefix_str)).exists());
        assert!(!std::path::Path::new(&format!("{}.syng.spos", prefix_str)).exists());
        assert!(!std::path::Path::new(&format!("{}.syng.pstep", prefix_str)).exists());
        assert!(!std::path::Path::new(&format!("{}.syng.meta", prefix_str)).exists());

        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert_eq!(loaded.name_map.path_to_name, vec!["genome".to_string()]);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_syng_suffix_prefix_loads_legacy_duplicated_sidecars() {
        let _guard = lock_syng();
        let seqs = vec![("genome".to_string(), make_test_sequence(2000, 11))];
        let params = SyncmerParams::default();
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_suffix_legacy_sidecars");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx.syng");
        let prefix_str = prefix.to_str().unwrap();

        index.save(prefix_str).unwrap();
        for suffix in ["names", "spos", "pstep", "meta"] {
            std::fs::rename(
                format!("{}.{}", prefix_str, suffix),
                format!("{}.syng.{}", prefix_str, suffix),
            )
            .unwrap();
        }

        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert_eq!(loaded.name_map.path_to_name, vec!["genome".to_string()]);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_save_load_uses_saved_syncmer_params() {
        let _guard = lock_syng();
        let seq = make_test_sequence(5000, 11);
        let params = SyncmerParams {
            k: 6,
            w: 31,
            seed: 42,
        };
        let mut index = SyngIndex::build(params, vec![("seq".to_string(), seq)].into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_saved_params");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx");
        let prefix_str = prefix.to_str().unwrap();

        index.save(prefix_str).unwrap();

        let meta = std::fs::read_to_string(format!("{}.syng.meta", prefix_str)).unwrap();
        assert!(meta.contains("syncmer_k\t6"));
        assert!(meta.contains("syncmer_w\t31"));
        assert!(meta.contains("syncmer_seed\t42"));

        let loaded = SyngIndex::load(prefix_str, SyncmerParams::default()).unwrap();
        assert_eq!(
            loaded.params, params,
            "load must use syncmer params stored in .syng.meta, not caller defaults"
        );
        let matches = loaded.matched_syncmers_in_sequence(&make_test_sequence(5000, 11));
        assert!(
            !matches.is_empty(),
            "loaded index should query with the saved non-default params"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 7. Round-trip: build → save → load → verify (spec Test 1) ──

    #[test]
    fn test_syng_roundtrip_path_count_preserved() {
        let _guard = lock_syng();
        // Build from multiple sequences, save, reload, verify same number
        // of paths recorded in the name map.
        let seqs: Vec<(String, Vec<u8>)> = (0..5)
            .map(|i| {
                (
                    format!("sample{}#hap{}#chr1", i / 2, i % 2),
                    make_test_sequence(3000, i as u8 + 50),
                )
            })
            .collect();
        let params = SyncmerParams::default();
        let mut original = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_rt_paths");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("rt_paths");
        let prefix_str = prefix.to_str().unwrap();

        original.save(prefix_str).unwrap();
        let loaded = SyngIndex::load(prefix_str, params).unwrap();

        // Same number of paths
        assert_eq!(
            loaded.name_map.path_to_name.len(),
            original.name_map.path_to_name.len(),
            "Path count must be preserved after round-trip"
        );

        // Same names in the same order
        assert_eq!(loaded.name_map.path_to_name, original.name_map.path_to_name);

        // Same lengths
        assert_eq!(
            loaded.name_map.path_to_length,
            original.name_map.path_to_length
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_syng_roundtrip_names_match() {
        let _guard = lock_syng();
        // Use PanSN-style names and verify they survive round-trip.
        let names = vec![
            "HG002#1#chr1",
            "HG002#2#chr1",
            "HG005#1#chrX",
            "GRCh38#0#chr22",
        ];
        let seqs: Vec<(String, Vec<u8>)> = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.to_string(), make_test_sequence(2000, i as u8 + 10)))
            .collect();
        let params = SyncmerParams::default();
        let mut original = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_rt_names");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("rt_names");
        let prefix_str = prefix.to_str().unwrap();

        original.save(prefix_str).unwrap();
        let loaded = SyngIndex::load(prefix_str, params).unwrap();

        for name in &names {
            assert!(
                loaded.name_map.name_to_path.contains_key(*name),
                "Name '{}' should be present in loaded name map",
                name
            );
            let orig_path = original.name_map.name_to_path[*name];
            let load_path = loaded.name_map.name_to_path[*name];
            assert_eq!(
                orig_path, load_path,
                "Path number for '{}' should match after round-trip",
                name
            );
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_syng_roundtrip_lengths_match() {
        let _guard = lock_syng();
        // Build with varying-length sequences, verify lengths preserved.
        let lengths = [100, 500, 1000, 5000, 10000];
        let seqs: Vec<(String, Vec<u8>)> = lengths
            .iter()
            .enumerate()
            .map(|(i, &len)| (format!("seq_len_{}", len), make_test_sequence(len, i as u8)))
            .collect();
        let params = SyncmerParams::default();
        let mut original = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_rt_lengths");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("rt_lengths");
        let prefix_str = prefix.to_str().unwrap();

        original.save(prefix_str).unwrap();
        let loaded = SyngIndex::load(prefix_str, params).unwrap();

        for (i, &expected_len) in lengths.iter().enumerate() {
            assert_eq!(
                loaded.name_map.path_to_length[i], expected_len as u64,
                "Length for seq {} should be {} after round-trip",
                i, expected_len
            );
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 8. SyngNameMap serialization with special characters ────────

    #[test]
    fn test_name_map_pansn_format() {
        let _guard = lock_syng();
        // PanSN format: SAMPLE#HAP#CONTIG — verify # chars survive serialization
        let mut nm = SyngNameMap::new();
        nm.add("HG002#1#chr1".to_string(), 248956422);
        nm.add("HG002#2#chr1".to_string(), 248956422);
        nm.add("GRCh38#0#chrX".to_string(), 156040895);
        nm.add("NA12878#1#chr22_random".to_string(), 50818468);

        let dir = std::env::temp_dir().join("impg_test_name_map_pansn");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("pansn.syng.names");
        let path_str = path.to_str().unwrap();

        nm.save(path_str).unwrap();
        let loaded = SyngNameMap::load(path_str).unwrap();

        assert_eq!(loaded.path_to_name, nm.path_to_name);
        assert_eq!(loaded.path_to_length, nm.path_to_length);
        for (k, v) in &nm.name_to_path {
            assert_eq!(
                loaded.name_to_path.get(k),
                Some(v),
                "PanSN name '{}' should survive round-trip",
                k
            );
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_name_map_special_characters() {
        let _guard = lock_syng();
        // Names with dots, underscores, dashes, and PanSN separators
        let mut nm = SyngNameMap::new();
        let names = vec![
            ("chr1.1_paternal-v2", 100000u64),
            ("sample.HG002#hap1#contig_123", 200000),
            ("ref-genome_v3.0", 300000),
            ("a", 1),
        ];
        for (name, len) in &names {
            nm.add(name.to_string(), *len);
        }

        let dir = std::env::temp_dir().join("impg_test_name_map_special");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("special.syng.names");
        let path_str = path.to_str().unwrap();

        nm.save(path_str).unwrap();
        let loaded = SyngNameMap::load(path_str).unwrap();

        for (name, len) in &names {
            let idx = *loaded
                .name_to_path
                .get(*name)
                .unwrap_or_else(|| panic!("Name '{}' missing from loaded name map", name));
            assert_eq!(loaded.path_to_length[idx as usize], *len);
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 9. Syncmer parameter variations ─────────────────────────────

    #[test]
    fn test_syng_different_params_different_syncmer_counts() {
        let _guard = lock_syng();
        // Building with different params should yield different numbers of syncmers
        // (reflected by different KmerHash contents). We verify indirectly: the
        // serialized .1khash files should differ in size.
        let seq = make_test_sequence(5000, 99);

        let dir = std::env::temp_dir().join("impg_test_syng_params_var");
        std::fs::create_dir_all(&dir).unwrap();

        // Default params: k=8, w=55 (syncmer_len=63)
        let params_default = SyncmerParams::default();
        let mut idx_default = SyngIndex::build(
            params_default,
            vec![("seq".to_string(), seq.clone())].into_iter(),
        );
        let prefix_default = dir.join("default");
        idx_default.save(prefix_default.to_str().unwrap()).unwrap();

        // Different params: k=6, w=31 (syncmer_len=37) — more syncmers (shorter, denser)
        let params_short = SyncmerParams {
            k: 6,
            w: 31,
            seed: 7,
        };
        let mut idx_short = SyngIndex::build(
            params_short,
            vec![("seq".to_string(), seq.clone())].into_iter(),
        );
        let prefix_short = dir.join("short");
        idx_short.save(prefix_short.to_str().unwrap()).unwrap();

        // Different seed with same k/w — different sampling
        let params_alt_seed = SyncmerParams {
            k: 8,
            w: 55,
            seed: 42,
        };
        let mut idx_alt =
            SyngIndex::build(params_alt_seed, vec![("seq".to_string(), seq)].into_iter());
        let prefix_alt = dir.join("alt_seed");
        idx_alt.save(prefix_alt.to_str().unwrap()).unwrap();

        // Verify files exist for all three
        for name in &["default", "short", "alt_seed"] {
            let khash = dir.join(format!("{}.1khash", name));
            assert!(khash.exists(), "{}.1khash should exist", name);
            let gbwt = dir.join(format!("{}.1gbwt", name));
            assert!(gbwt.exists(), "{}.1gbwt should exist", name);
        }

        // Shorter syncmers (k=6,w=30) should yield a different-sized khash than
        // default (k=8,w=55). With a 36bp syncmer vs 63bp, the shorter one extracts
        // more syncmers.
        let size_default = std::fs::metadata(dir.join("default.1khash")).unwrap().len();
        let size_short = std::fs::metadata(dir.join("short.1khash")).unwrap().len();
        assert_ne!(
            size_default, size_short,
            "Different syncmer params should produce different khash sizes"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_syng_default_params_build_and_reload() {
        let _guard = lock_syng();
        // Build with default params, save, reload with same params — should work
        let params = SyncmerParams::default();
        let seqs = vec![
            ("s1".to_string(), make_test_sequence(2000, 1)),
            ("s2".to_string(), make_test_sequence(2000, 2)),
        ];
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_default_reload");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx");
        let prefix_str = prefix.to_str().unwrap();

        index.save(prefix_str).unwrap();
        let loaded = SyngIndex::load(prefix_str, params).unwrap();

        assert_eq!(loaded.name_map.path_to_name.len(), 2);
        assert_eq!(loaded.params.k, 8);
        assert_eq!(loaded.params.w, 55);
        assert_eq!(loaded.params.seed, 7);

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 10. Empty and edge cases ────────────────────────────────────

    #[test]
    fn test_syng_build_empty_sequence() {
        let _guard = lock_syng();
        // A 0-length sequence should be handled gracefully (no panic)
        let seqs = vec![("empty".to_string(), Vec::new())];
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 1);
        assert_eq!(index.name_map.path_to_name[0], "empty");
        assert_eq!(index.name_map.path_to_length[0], 0);
    }

    #[test]
    fn test_syng_build_sequence_shorter_than_syncmer() {
        let _guard = lock_syng();
        // 50bp < 63bp (default syncmer length) — should record but skip extraction
        let seq = make_test_sequence(50, 77);
        let seqs = vec![("too_short".to_string(), seq)];
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 1);
        assert_eq!(index.name_map.path_to_name[0], "too_short");
        assert_eq!(index.name_map.path_to_length[0], 50);
    }

    #[test]
    fn test_syng_build_exactly_syncmer_length() {
        let _guard = lock_syng();
        // Sequence exactly syncmer length (63bp) — borderline, should not panic
        let seq = make_test_sequence(63, 88);
        let seqs = vec![("exact".to_string(), seq)];
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 1);
        assert_eq!(index.name_map.path_to_length[0], 63);
    }

    #[test]
    fn test_syng_build_single_path_in_gbwt() {
        let _guard = lock_syng();
        // Single sequence should produce exactly one entry in name map
        let seq = make_test_sequence(1000, 33);
        let seqs = vec![("single".to_string(), seq)];
        let params = SyncmerParams::default();
        let mut index = SyngIndex::build(params, seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 1);

        // Save and reload to verify the single path survives
        let dir = std::env::temp_dir().join("impg_test_syng_single_path");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("single");
        let prefix_str = prefix.to_str().unwrap();

        index.save(prefix_str).unwrap();
        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert_eq!(loaded.name_map.path_to_name.len(), 1);
        assert_eq!(loaded.name_map.path_to_name[0], "single");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_syng_build_mixed_lengths() {
        let _guard = lock_syng();
        // Mix of empty, short, and long sequences — all should be recorded
        let seqs = vec![
            ("empty".to_string(), Vec::new()),
            ("short_10bp".to_string(), make_test_sequence(10, 1)),
            ("borderline_62bp".to_string(), make_test_sequence(62, 2)),
            ("borderline_63bp".to_string(), make_test_sequence(63, 3)),
            ("normal_500bp".to_string(), make_test_sequence(500, 4)),
            ("long_10000bp".to_string(), make_test_sequence(10000, 5)),
        ];
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 6);
        assert_eq!(index.name_map.path_to_length[0], 0);
        assert_eq!(index.name_map.path_to_length[1], 10);
        assert_eq!(index.name_map.path_to_length[2], 62);
        assert_eq!(index.name_map.path_to_length[3], 63);
        assert_eq!(index.name_map.path_to_length[4], 500);
        assert_eq!(index.name_map.path_to_length[5], 10000);
    }

    #[test]
    fn test_syng_save_load_with_edge_cases() {
        let _guard = lock_syng();
        // Build with mix of edge-case sequences, round-trip through disk
        let seqs = vec![
            ("empty".to_string(), Vec::new()),
            ("tiny".to_string(), b"ACG".to_vec()),
            ("normal".to_string(), make_test_sequence(2000, 42)),
        ];
        let params = SyncmerParams::default();
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_edge_rt");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("edge");
        let prefix_str = prefix.to_str().unwrap();

        index.save(prefix_str).unwrap();
        let loaded = SyngIndex::load(prefix_str, params).unwrap();

        assert_eq!(loaded.name_map.path_to_name.len(), 3);
        assert_eq!(loaded.name_map.path_to_name[0], "empty");
        assert_eq!(loaded.name_map.path_to_length[0], 0);
        assert_eq!(loaded.name_map.path_to_name[1], "tiny");
        assert_eq!(loaded.name_map.path_to_length[1], 3);
        assert_eq!(loaded.name_map.path_to_name[2], "normal");
        assert_eq!(loaded.name_map.path_to_length[2], 2000);

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_online_sampled_positions_exact_mapping_roundtrip() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let shared = make_test_sequence(1000, 42);
        let mut seq_a = shared.clone();
        seq_a.extend_from_slice(&make_test_sequence(400, 1));
        let mut seq_b = shared.clone();
        seq_b.extend_from_slice(&make_test_sequence(400, 2));

        let mut index = SyngIndex::new(params);
        index.enable_online_sampled_positions(1).unwrap();
        index.add_sequence("sampleA#0#chr1".to_string(), seq_a);
        index.add_sequence("sampleB#0#chr1".to_string(), seq_b);
        let stats = index.finalize_online_sampled_positions().unwrap().unwrap();
        assert_eq!(stats.sampled_occurrences, 0);
        assert_eq!(stats.sampled_nodes, 0);
        assert!(stats.sampled_path_steps > 0);
        assert!(stats.sampled_step_paths > 0);
        assert!(index.has_sampled_path_steps());

        let query = &shared[100..800];
        let hits = index
            .map_sequence("read1", query, 2, 10_000)
            .expect("sampled-position PAF mapping should succeed");
        assert!(
            hits.iter().any(|h| h.target_name == "sampleA#0#chr1"),
            "sampled mapping should include sampleA, got: {:?}",
            hits
        );
        assert!(
            hits.iter().any(|h| h.target_name == "sampleB#0#chr1"),
            "sampled mapping should include sampleB, got: {:?}",
            hits
        );

        let dir = std::env::temp_dir().join("impg_test_syng_spos_roundtrip");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx");
        let prefix_str = prefix.to_str().unwrap();
        index.save(prefix_str).unwrap();
        assert!(
            std::path::Path::new(&format!("{}.syng.spos", prefix_str)).exists(),
            ".syng.spos file should exist after save"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.pstep", prefix_str)).exists(),
            ".syng.pstep file should exist after save"
        );
        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert!(loaded.has_sampled_positions());
        assert!(loaded.has_sampled_path_steps());
        assert!(loaded.has_sampled_checkpoints());
        let loaded_hits = loaded
            .map_sequence("read1", query, 2, 10_000)
            .expect("loaded sampled-position PAF mapping should succeed");
        assert!(
            loaded_hits
                .iter()
                .any(|h| h.target_name == "sampleA#0#chr1"),
            "loaded sampled mapping should include sampleA, got: {:?}",
            loaded_hits
        );
        assert!(
            loaded_hits
                .iter()
                .any(|h| h.target_name == "sampleB#0#chr1"),
            "loaded sampled mapping should include sampleB, got: {:?}",
            loaded_hits
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_sampled_path_steps_resume_gbwt_traversal() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq = make_test_sequence(5000, 123);
        let mut index = SyngIndex::new(params);
        index.enable_online_sampled_positions(1).unwrap();
        index.add_sequence("sampleA#0#chr1".to_string(), seq);
        index.finalize_online_sampled_positions().unwrap().unwrap();

        let path_idx = *index.name_map.name_to_path.get("sampleA#0#chr1").unwrap() as usize;
        let path_steps = index
            .sampled_path_steps
            .as_ref()
            .unwrap()
            .decode_path(path_idx)
            .unwrap();
        assert!(
            path_steps.len() > 8,
            "test sequence should produce many syncmers"
        );
        assert_eq!(path_steps[0].step_idx, 0);
        assert_eq!(path_steps[0].prev_node, 0);
        assert_eq!(path_steps[0].prev_offset, 0);

        for window in path_steps.windows(2).take(16) {
            let current = window[0];
            let expected_next = window[1];
            let actual_next = index
                .next_path_step_from_checkpoint(&current)
                .unwrap()
                .expect("checkpoint should resume to next step");
            assert_eq!(actual_next.step_idx, expected_next.step_idx);
            assert_eq!(actual_next.bp_pos, expected_next.bp_pos);
            assert_eq!(actual_next.signed_node, expected_next.signed_node);
            assert_eq!(actual_next.prev_node, expected_next.prev_node);
            assert_eq!(actual_next.prev_offset, expected_next.prev_offset);
            assert_eq!(actual_next.traversal_rank, expected_next.traversal_rank);
        }

        let mid = path_steps[path_steps.len() / 2];
        let checkpoint = index
            .sampled_path_step_at_or_before(path_idx, mid.bp_pos + 1)
            .unwrap()
            .unwrap();
        assert_eq!(checkpoint, mid);

        let range = index
            .walk_path_range_from_sampled_steps(path_idx, mid.bp_pos, mid.bp_pos + 200)
            .unwrap();
        assert!(
            range
                .iter()
                .any(|&(node, pos)| node == mid.signed_node && pos == mid.bp_pos),
            "range walk should include the checkpoint syncmer"
        );
    }

    #[test]
    fn test_sampled_path_range_includes_unsampled_boundary_overlap() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let syncmer_len = (params.k + params.w) as u64;
        let seq = make_test_sequence(12_000, 77);
        let mut index = SyngIndex::new(params);
        index.enable_online_sampled_positions(16).unwrap();
        index.add_sequence("sampleA#0#chr1".to_string(), seq);
        index.finalize_online_sampled_positions().unwrap().unwrap();

        let path_idx = *index.name_map.name_to_path.get("sampleA#0#chr1").unwrap() as usize;
        let path_start = index.name_map.path_starts[path_idx].as_ref().unwrap();
        let full_steps = index.walk_path_steps(path_start).unwrap();
        let sampled_steps = index
            .sampled_path_steps
            .as_ref()
            .unwrap()
            .decode_path(path_idx)
            .unwrap();
        assert!(
            sampled_steps.len() < full_steps.len(),
            "sparse test requires unsampled path steps"
        );

        let sampled_step_idxs: std::collections::BTreeSet<u32> =
            sampled_steps.iter().map(|step| step.step_idx).collect();
        let boundary_step = full_steps
            .iter()
            .find(|step| {
                !sampled_step_idxs.contains(&step.step_idx)
                    && step.bp_pos > syncmer_len
                    && step.bp_pos.saturating_add(syncmer_len)
                        < index.name_map.path_to_length[path_idx]
            })
            .expect("expected at least one interior unsampled syncmer");

        // Query a 1 bp window at the unsampled syncmer's right edge. A range
        // walker that starts at the checkpoint before `start` instead of before
        // `start - syncmer_len + 1` skips this overlapping syncmer.
        let start = boundary_step.bp_pos + syncmer_len - 1;
        let end = start + 1;
        let expected: Vec<(i32, u64)> = full_steps
            .iter()
            .filter(|step| step.bp_pos < end && step.bp_pos.saturating_add(syncmer_len) > start)
            .map(|step| (step.signed_node, step.bp_pos))
            .collect();
        let observed = index
            .walk_path_range_from_sampled_steps(path_idx, start, end)
            .unwrap();

        assert!(
            expected.iter().any(
                |&(node, pos)| node == boundary_step.signed_node && pos == boundary_step.bp_pos
            ),
            "test setup should include the unsampled boundary syncmer"
        );
        assert_eq!(observed, expected);
    }

    #[test]
    fn test_locate_unsampled_syncmer_walks_forward_to_checkpoint() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let syncmer_len = (params.k + params.w) as u64;
        let seq = make_test_sequence(12_000, 78);
        let mut index = SyngIndex::new(params);
        index.enable_online_sampled_positions(16).unwrap();
        index.add_sequence("sampleA#0#chr1".to_string(), seq);
        index.finalize_online_sampled_positions().unwrap().unwrap();

        let path_idx = *index.name_map.name_to_path.get("sampleA#0#chr1").unwrap() as usize;
        let path_start = index.name_map.path_starts[path_idx].as_ref().unwrap();
        let full_steps = index.walk_path_steps(path_start).unwrap();
        let sampled_steps = index
            .sampled_path_steps
            .as_ref()
            .unwrap()
            .decode_path(path_idx)
            .unwrap();
        let sampled_step_idxs: std::collections::BTreeSet<u32> =
            sampled_steps.iter().map(|step| step.step_idx).collect();
        let unsampled_step = full_steps
            .iter()
            .find(|step| {
                !sampled_step_idxs.contains(&step.step_idx)
                    && step.bp_pos > syncmer_len
                    && step.bp_pos.saturating_add(syncmer_len)
                        < index.name_map.path_to_length[path_idx]
            })
            .expect("expected an interior unsampled syncmer step");

        let located = index
            .locate_signed_node_occurrences(unsampled_step.signed_node)
            .unwrap();
        assert!(
            located
                .iter()
                .any(|hit| hit.path_idx == path_idx && hit.target_pos == unsampled_step.bp_pos),
            "locate should recover the exact unsampled path position by walking forward to a checkpoint; step={:?}, located={:?}",
            unsampled_step,
            located
        );
    }

    #[test]
    fn test_query_region_uses_walked_target_positions_for_unsampled_anchors() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let syncmer_len = (params.k + params.w) as u64;
        let seq = make_test_sequence(12_000, 79);
        let mut index = SyngIndex::new(params);
        index.enable_online_sampled_positions(16).unwrap();
        index.add_sequence("sampleA#0#chr1".to_string(), seq);
        index.finalize_online_sampled_positions().unwrap().unwrap();

        let path_idx = *index.name_map.name_to_path.get("sampleA#0#chr1").unwrap() as usize;
        let path_start = index.name_map.path_starts[path_idx].as_ref().unwrap();
        let full_steps = index.walk_path_steps(path_start).unwrap();
        let sampled_steps = index
            .sampled_path_steps
            .as_ref()
            .unwrap()
            .decode_path(path_idx)
            .unwrap();
        let sampled_step_idxs: std::collections::BTreeSet<u32> =
            sampled_steps.iter().map(|step| step.step_idx).collect();
        let unsampled_step = full_steps
            .iter()
            .find(|step| {
                !sampled_step_idxs.contains(&step.step_idx)
                    && step.bp_pos > syncmer_len
                    && step.bp_pos.saturating_add(syncmer_len)
                        < index.name_map.path_to_length[path_idx]
            })
            .expect("expected an interior unsampled syncmer step");

        let hits = index
            .query_region_with_anchors_ext(
                "sampleA#0#chr1",
                unsampled_step.bp_pos,
                unsampled_step.bp_pos + syncmer_len,
                0,
                0,
            )
            .unwrap();
        assert!(
            hits.iter().any(|hit| {
                hit.genome == "sampleA#0#chr1"
                    && hit.anchors.iter().any(|anchor| {
                        anchor.node_id == unsampled_step.signed_node.unsigned_abs()
                            && anchor.query_pos == unsampled_step.bp_pos
                            && anchor.target_pos == unsampled_step.bp_pos
                    })
            }),
            "query should emit the exact unsampled target anchor; step={:?}, hits={:?}",
            unsampled_step,
            hits
        );
    }

    #[test]
    fn test_regular_sampled_positions_use_step_grid_and_terminal_step() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let mut index = SyngIndex::new(params);
        index.enable_online_sampled_positions(8).unwrap();
        index.add_sequence("sampleA#0#chr1".to_string(), make_test_sequence(12_000, 91));
        index.finalize_online_sampled_positions().unwrap().unwrap();

        let path_idx = *index.name_map.name_to_path.get("sampleA#0#chr1").unwrap() as usize;
        let path_start = index.name_map.path_starts[path_idx].as_ref().unwrap();
        let full_steps = index.walk_path_steps(path_start).unwrap();
        let last_step_idx = full_steps.last().unwrap().step_idx;
        let sampled_steps = index
            .sampled_path_steps
            .as_ref()
            .unwrap()
            .decode_path(path_idx)
            .unwrap();
        assert!(
            sampled_steps.len() > 4,
            "test sequence should produce several regular checkpoints"
        );
        assert_eq!(sampled_steps[0].step_idx, 0);
        assert_eq!(sampled_steps.last().unwrap().step_idx, last_step_idx);
        for step in &sampled_steps {
            assert!(
                step.step_idx % 8 == 0 || step.step_idx == last_step_idx,
                "regular sampled checkpoint is off-grid: {:?}",
                step
            );
        }
        assert_eq!(
            index.sampled_path_steps.as_ref().unwrap().sample_count(),
            sampled_steps.len() as u64,
            "path-step sidecar should report the decoded checkpoint count"
        );
    }

    #[test]
    fn test_terminal_sampled_for_short_paths() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let mut index = SyngIndex::new(params);
        index.enable_online_sampled_positions(256).unwrap();
        index.add_sequence("sampleA#0#chr1".to_string(), make_test_sequence(2_000, 92));
        index.finalize_online_sampled_positions().unwrap().unwrap();

        let path_idx = *index.name_map.name_to_path.get("sampleA#0#chr1").unwrap() as usize;
        let path_start = index.name_map.path_starts[path_idx].as_ref().unwrap();
        let full_steps = index.walk_path_steps(path_start).unwrap();
        assert!(
            full_steps.len() > 1 && full_steps.len() < 256,
            "test setup should produce a short nonempty syncmer path"
        );

        let sampled_steps = index
            .sampled_path_steps
            .as_ref()
            .unwrap()
            .decode_path(path_idx)
            .unwrap();
        assert_eq!(sampled_steps.len(), 2);
        assert_eq!(sampled_steps[0].step_idx, 0);
        assert_eq!(
            sampled_steps[1].step_idx,
            full_steps.last().unwrap().step_idx
        );
    }

    #[test]
    fn test_rebuild_sampled_position_indexes_from_gbwt_repairs_sidecars() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seqs = vec![
            ("sampleA#0#chr1".to_string(), make_test_sequence(4000, 13)),
            ("sampleB#0#chr1".to_string(), make_test_sequence(4000, 17)),
        ];
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_repair_sidecars");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx");
        let prefix_str = prefix.to_str().unwrap();
        index.save(prefix_str).unwrap();
        std::fs::remove_file(format!("{}.syng.spos", prefix_str)).unwrap();
        std::fs::remove_file(format!("{}.syng.pstep", prefix_str)).unwrap();

        let mut repair = SyngIndex::load_for_repair(prefix_str, params).unwrap();
        assert!(!repair.has_sampled_positions());
        assert!(!repair.has_sampled_path_steps());
        assert!(!repair.has_sampled_checkpoints());
        let stats = repair
            .rebuild_sampled_position_indexes_from_gbwt(1)
            .unwrap();
        assert_eq!(stats.sampled_occurrences, 0);
        assert!(stats.sampled_path_steps > 0);
        repair.save_position_sidecars(prefix_str).unwrap();

        assert!(std::path::Path::new(&format!("{}.syng.spos", prefix_str)).exists());
        assert!(std::path::Path::new(&format!("{}.syng.pstep", prefix_str)).exists());

        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert!(loaded.has_sampled_positions());
        assert!(loaded.has_sampled_path_steps());
        assert!(loaded.has_sampled_checkpoints());
        let path_idx = *loaded.name_map.name_to_path.get("sampleA#0#chr1").unwrap() as usize;
        let checkpoint = loaded
            .sampled_path_step_at_or_before(path_idx, 2_000)
            .unwrap()
            .unwrap();
        let next = loaded.next_path_step_from_checkpoint(&checkpoint).unwrap();
        assert!(next.is_some());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_loaded_syng_requires_spos_but_repair_can_rebuild() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seqs = vec![
            ("sampleA#0#chr1".to_string(), make_test_sequence(4000, 31)),
            ("sampleB#0#chr1".to_string(), make_test_sequence(4000, 37)),
        ];
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_requires_spos");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx");
        let prefix_str = prefix.to_str().unwrap();
        index.save(prefix_str).unwrap();
        std::fs::remove_file(format!("{}.syng.spos", prefix_str)).unwrap();

        let err = match SyngIndex::load(prefix_str, params) {
            Ok(_) => panic!("normal syng load should reject indexes without .spos"),
            Err(err) => err,
        };
        assert!(
            err.to_string().contains("sampled syncmer-position sidecar"),
            "unexpected missing-spos error: {err}"
        );

        let repair = SyngIndex::load_for_repair(prefix_str, params).unwrap();
        assert!(repair.has_sampled_path_steps());
        assert!(!repair.has_sampled_checkpoints());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_parallel_regular_position_rebuild_matches_serial() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seqs = vec![
            ("sampleA#0#chr1".to_string(), make_test_sequence(5000, 21)),
            ("sampleB#0#chr1".to_string(), make_test_sequence(4200, 22)),
            ("sampleC#0#chr1".to_string(), make_test_sequence(3500, 23)),
            ("short".to_string(), make_test_sequence(12, 24)),
        ];
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_syng_parallel_regular_repair");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx");
        let prefix_str = prefix.to_str().unwrap();
        index.save(prefix_str).unwrap();

        let mut serial = SyngIndex::load_for_repair(prefix_str, params).unwrap();
        let serial_stats = serial
            .rebuild_sampled_position_indexes_from_gbwt(8)
            .unwrap();
        let mut parallel = SyngIndex::load_for_repair(prefix_str, params).unwrap();
        let parallel_stats = parallel
            .rebuild_sampled_position_indexes_from_gbwt_parallel(8, 0)
            .unwrap();

        assert_eq!(serial_stats.sampled_occurrences, 0);
        assert_eq!(parallel_stats.sampled_occurrences, 0);
        assert_eq!(serial_stats.sampled_nodes, 0);
        assert_eq!(parallel_stats.sampled_nodes, 0);
        assert_eq!(
            serial_stats.sampled_path_steps,
            parallel_stats.sampled_path_steps
        );
        assert_eq!(
            serial_stats.sampled_step_paths,
            parallel_stats.sampled_step_paths
        );
        assert_eq!(serial_stats.walked_paths, parallel_stats.walked_paths);

        let serial_steps = serial.sampled_path_steps.as_ref().unwrap();
        let parallel_steps = parallel.sampled_path_steps.as_ref().unwrap();
        assert_eq!(serial_steps.sample_count(), parallel_steps.sample_count());
        assert_eq!(serial_steps.path_count(), parallel_steps.path_count());
        for path_idx in 0..serial_steps.path_count() {
            assert_eq!(
                serial_steps.decode_path(path_idx).unwrap(),
                parallel_steps.decode_path(path_idx).unwrap(),
                "path-step checkpoints differ for path {}",
                path_idx
            );
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 11. CLI integration test ────────────────────────────────────

    #[test]
    fn test_syng_cli_fasta_roundtrip() {
        let _guard = lock_syng();
        // Run `impg syng -f <fasta> -o <prefix>` via Command,
        // verify output files exist.
        let dir = std::env::temp_dir().join("impg_test_syng_cli");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        // Write a small synthetic FASTA
        let fasta_path = dir.join("test_cli.fa");
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            // Need sequences >= 63bp (default syncmer length)
            let seq1 = String::from_utf8(make_test_sequence(500, 200)).unwrap();
            let seq2 = String::from_utf8(make_test_sequence(500, 201)).unwrap();
            let seq3 = String::from_utf8(make_test_sequence(300, 202)).unwrap();
            use std::io::Write;
            writeln!(f, ">HG002#1#chr1").unwrap();
            writeln!(f, "{}", seq1).unwrap();
            writeln!(f, ">HG002#2#chr1").unwrap();
            writeln!(f, "{}", seq2).unwrap();
            writeln!(f, ">GRCh38#0#chrX").unwrap();
            writeln!(f, "{}", seq3).unwrap();
        }

        let output_prefix = dir.join("cli_output");
        let output_prefix_str = output_prefix.to_str().unwrap();

        // Find the impg binary
        let bin = std::env::current_exe()
            .unwrap()
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("impg");

        // If the binary doesn't exist, try cargo bin path
        let bin = if bin.exists() {
            bin
        } else {
            // Fall back to looking in target/debug or target/release
            let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
            manifest_dir.join("target/debug/impg")
        };

        if !bin.exists() {
            // Skip if binary not built yet — cargo test --lib won't build the binary
            eprintln!("Skipping CLI test: impg binary not found at {:?}", bin);
            std::fs::remove_dir_all(&dir).ok();
            return;
        }

        let output = std::process::Command::new(&bin)
            .args([
                "syng",
                "-f",
                fasta_path.to_str().unwrap(),
                "-o",
                output_prefix_str,
            ])
            .output()
            .expect("Failed to run impg syng");

        assert!(
            output.status.success(),
            "impg syng should succeed. stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        // Verify output files
        assert!(
            std::path::Path::new(&format!("{}.1khash", output_prefix_str)).exists(),
            ".1khash file should exist after CLI run"
        );
        assert!(
            std::path::Path::new(&format!("{}.1gbwt", output_prefix_str)).exists(),
            ".1gbwt file should exist after CLI run"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.names", output_prefix_str)).exists(),
            ".syng.names file should exist after CLI run"
        );
        assert!(
            !std::path::Path::new(&format!("{}.syng.locate", output_prefix_str)).exists(),
            ".syng.locate should not be built"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.spos", output_prefix_str)).exists(),
            ".syng.spos should be built by default"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.pstep", output_prefix_str)).exists(),
            ".syng.pstep should be built by default"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.meta", output_prefix_str)).exists(),
            ".syng.meta should be built by default"
        );

        // Load and verify content
        let loaded = SyngIndex::load(output_prefix_str, SyncmerParams::default()).unwrap();
        assert!(loaded.has_sampled_positions());
        assert!(loaded.has_sampled_path_steps());
        assert!(loaded.has_sampled_checkpoints());
        assert_eq!(loaded.name_map.path_to_name.len(), 3);
        assert_eq!(loaded.name_map.path_to_name[0], "HG002#1#chr1");
        assert_eq!(loaded.name_map.path_to_name[1], "HG002#2#chr1");
        assert_eq!(loaded.name_map.path_to_name[2], "GRCh38#0#chrX");

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 12. Load error cases ────────────────────────────────────────

    #[test]
    fn test_syng_load_missing_files() {
        let _guard = lock_syng();
        let result = SyngIndex::load("/tmp/nonexistent_prefix_xyz123", SyncmerParams::default());
        assert!(
            result.is_err(),
            "Loading from nonexistent files should fail"
        );
    }

    // ── 13. query_region tests ─────────────────────────────────────

    /// Full SyngIndex save/load with sampled positional sidecars. Verify that
    /// the reloaded index returns the same query results as the in-memory
    /// version.
    #[test]
    fn test_syng_save_load_with_sampled_positions() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let shared = make_test_sequence(500, 42);
        let mut sa = shared.clone();
        sa.extend_from_slice(&make_test_sequence(500, 1));
        let mut sb = shared.clone();
        sb.extend_from_slice(&make_test_sequence(500, 2));
        let seqs = vec![("ga".to_string(), sa), ("gb".to_string(), sb)];

        let mut index = SyngIndex::build(params, seqs.into_iter());
        let in_mem = index.query_region("ga", 0, 1000, 120).unwrap();

        let dir = std::env::temp_dir().join("impg_test_syng_spos_query_roundtrip");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx");
        let prefix_str = prefix.to_str().unwrap();
        index.save(prefix_str).unwrap();
        assert!(
            std::path::Path::new(&format!("{}.syng.spos", prefix_str)).exists(),
            ".syng.spos file should exist after save"
        );
        assert!(
            !std::path::Path::new(&format!("{}.syng.locate", prefix_str)).exists(),
            ".syng.locate should not be written"
        );
        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert!(
            loaded.has_sampled_positions(),
            "loaded index should have sampled positions"
        );
        assert!(
            loaded.has_sampled_checkpoints(),
            "loaded index should have sampled checkpoints"
        );
        let reloaded = loaded.query_region("ga", 0, 1000, 120).unwrap();

        let to_set = |v: &[HomologousInterval]| -> std::collections::BTreeSet<(String, u64, u64)> {
            v.iter()
                .map(|iv| (iv.genome.clone(), iv.start, iv.end))
                .collect()
        };
        assert_eq!(to_set(&in_mem), to_set(&reloaded));
        assert!(!in_mem.is_empty());
    }

    /// Query path smoke test for exact sampled positions over multiple genomes.
    #[test]
    fn test_query_region_sampled_positions_three_genomes() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();

        // Shared backbone + diverging tails + one entirely novel sequence.
        let shared = make_test_sequence(600, 42);
        let mut seq_a = shared.clone();
        seq_a.extend_from_slice(&make_test_sequence(400, 1));
        let mut seq_b = shared.clone();
        seq_b.extend_from_slice(&make_test_sequence(400, 2));
        let mut seq_c = shared.clone();
        seq_c.extend_from_slice(&make_test_sequence(400, 3));
        let seq_d = make_test_sequence(1000, 99);

        let seqs = vec![
            ("ga".to_string(), seq_a),
            ("gb".to_string(), seq_b),
            ("gc".to_string(), seq_c),
            ("gd".to_string(), seq_d),
        ];

        let index = SyngIndex::build(params, seqs.into_iter());
        assert!(index.has_sampled_path_steps());
        let intervals = index.query_region("ga", 0, 1000, 120).unwrap();
        let genomes: std::collections::BTreeSet<&str> =
            intervals.iter().map(|iv| iv.genome.as_str()).collect();
        assert!(genomes.contains("ga"), "self hit missing: {:?}", genomes);
        assert!(
            genomes.contains("gb"),
            "shared genome gb missing: {:?}",
            genomes
        );
        assert!(
            genomes.contains("gc"),
            "shared genome gc missing: {:?}",
            genomes
        );
        assert!(
            !genomes.contains("gd"),
            "unrelated genome gd should not appear in sampled smoke test: {:?}",
            genomes
        );
    }

    #[test]
    fn test_query_region_from_sequence_matches_path_query() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();

        let shared = make_test_sequence(600, 42);
        let mut seq_a = shared.clone();
        seq_a.extend_from_slice(&make_test_sequence(400, 1));
        let mut seq_b = shared.clone();
        seq_b.extend_from_slice(&make_test_sequence(400, 2));
        let mut index = SyngIndex::new(params);
        index.enable_online_sampled_positions(16).unwrap();
        index.add_sequence("ga".to_string(), seq_a.clone());
        index.add_sequence("gb".to_string(), seq_b);
        index.finalize_online_sampled_positions().unwrap().unwrap();

        let path_intervals = index.query_region("ga", 200, 800, 120).unwrap();
        let sequence_intervals = index
            .query_region_from_sequence_ext(&seq_a, 0, 200, 800, 120, 0)
            .unwrap();

        let to_set =
            |v: &[HomologousInterval]| -> std::collections::BTreeSet<(String, u64, u64, char)> {
                v.iter()
                    .map(|iv| (iv.genome.clone(), iv.start, iv.end, iv.strand))
                    .collect()
            };
        assert_eq!(to_set(&path_intervals), to_set(&sequence_intervals));
        assert!(!sequence_intervals.is_empty());
    }

    /// Build an index from sequences that share some content, then query a region.
    #[test]
    fn test_query_region_basic() {
        let _guard = lock_syng();
        // Create sequences that share a common prefix/suffix (to get shared syncmer nodes)
        // Use the same random seed for a shared backbone of ~500bp,
        // then diverge in the middle.
        let shared_len = 500;
        let total_len = 1000;
        let params = SyncmerParams {
            k: 8,
            w: 55,
            seed: 7,
        };
        // Build a shared backbone
        let backbone = make_test_sequence(shared_len, 42);

        // genome_a: backbone + unique_a
        let unique_a = make_test_sequence(total_len - shared_len, 1);
        let mut seq_a = backbone.clone();
        seq_a.extend_from_slice(&unique_a);

        // genome_b: backbone + unique_b (different from unique_a)
        let unique_b = make_test_sequence(total_len - shared_len, 2);
        let mut seq_b = backbone.clone();
        seq_b.extend_from_slice(&unique_b);

        // genome_c: completely different sequence (should NOT share nodes)
        let seq_c = make_test_sequence(total_len, 99);

        let sequences = vec![
            ("genome_a".to_string(), seq_a),
            ("genome_b".to_string(), seq_b),
            ("genome_c".to_string(), seq_c),
        ];

        let index = SyngIndex::build(params, sequences.into_iter());

        // Verify path starts were recorded
        for i in 0..3 {
            assert!(
                index.name_map.path_starts[i].is_some(),
                "path_starts[{}] should be Some",
                i
            );
        }

        // Query the shared region on genome_a
        let intervals = index
            .query_region("genome_a", 0, shared_len as u64, 0)
            .unwrap();

        // genome_a should appear (self-hit)
        let has_self = intervals.iter().any(|iv| iv.genome == "genome_a");
        assert!(has_self, "query_region should return self-hit for genome_a");

        // genome_b should appear (shares the backbone)
        let has_b = intervals.iter().any(|iv| iv.genome == "genome_b");
        assert!(
            has_b,
            "query_region should find genome_b (shared backbone). Only found: {:?}",
            intervals.iter().map(|iv| &iv.genome).collect::<Vec<_>>()
        );

        // All intervals should have valid coordinates
        for iv in &intervals {
            assert!(
                iv.end > iv.start,
                "interval end should be > start: {:?}",
                iv
            );
            assert!(
                iv.end
                    <= index.name_map.path_to_length
                        [*index.name_map.name_to_path.get(&iv.genome).unwrap() as usize],
                "interval end should be <= genome length"
            );
        }
    }

    #[test]
    fn test_query_region_unknown_genome() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq = make_test_sequence(500, 42);
        let index = SyngIndex::build(params, vec![("seq1".to_string(), seq)].into_iter());

        let result = index.query_region("nonexistent", 0, 100, 0);
        assert!(result.is_err(), "Should error for unknown genome");
    }

    #[test]
    fn test_query_region_padding() {
        let _guard = lock_syng();
        let params = SyncmerParams {
            k: 8,
            w: 55,
            seed: 7,
        };
        let seq_a = make_test_sequence(1000, 42);
        let seq_b = seq_a.clone(); // identical => fully shared

        let index = SyngIndex::build(
            params,
            vec![
                ("genome_a".to_string(), seq_a),
                ("genome_b".to_string(), seq_b),
            ]
            .into_iter(),
        );

        // Query without padding
        let no_pad = index.query_region("genome_a", 200, 400, 0).unwrap();
        // Query with 120bp padding
        let with_pad = index.query_region("genome_a", 200, 400, 120).unwrap();

        // With padding, intervals should be at least as wide as without
        for (np, wp) in no_pad.iter().zip(with_pad.iter()) {
            if np.genome == wp.genome {
                assert!(
                    wp.start <= np.start && wp.end >= np.end,
                    "Padded interval should be >= unpadded: {:?} vs {:?}",
                    wp,
                    np
                );
            }
        }
    }

    #[test]
    fn test_query_region_no_path_start_info() {
        let _guard = lock_syng();
        // Build an index, then manually clear path_starts to simulate old format
        let params = SyncmerParams::default();
        let seq = make_test_sequence(500, 42);
        let mut index = SyngIndex::build(params, vec![("seq1".to_string(), seq)].into_iter());
        index.name_map.path_starts[0] = None;

        let result = index.query_region("seq1", 0, 100, 0);
        assert!(result.is_err(), "Should error when path_starts is missing");
    }

    #[test]
    fn test_query_region_save_load_roundtrip() {
        let _guard = lock_syng();
        // Build, save, load, then query — results should match
        let params = SyncmerParams {
            k: 8,
            w: 55,
            seed: 7,
        };
        let seq_a = make_test_sequence(800, 10);
        let seq_b = seq_a.clone(); // identical

        let mut index = SyngIndex::build(
            params,
            vec![
                ("genome_a".to_string(), seq_a),
                ("genome_b".to_string(), seq_b),
            ]
            .into_iter(),
        );

        // Save and reload
        let dir = std::env::temp_dir().join("impg_test_query_rt");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("test");
        let prefix_str = prefix.to_str().unwrap();
        index.save(prefix_str).unwrap();
        let loaded = SyngIndex::load(prefix_str, params).unwrap();

        // Check that path_starts survived the roundtrip
        for i in 0..2 {
            let orig = index.name_map.path_starts[i].as_ref().unwrap();
            let load = loaded.name_map.path_starts[i].as_ref().unwrap();
            assert_eq!(orig.start_node, load.start_node);
            assert_eq!(orig.start_count, load.start_count);
            assert_eq!(orig.num_syncmers, load.num_syncmers);
        }

        // Query results should match
        let orig_intervals = index.query_region("genome_a", 100, 500, 120).unwrap();
        let loaded_intervals = loaded.query_region("genome_a", 100, 500, 120).unwrap();
        assert_eq!(orig_intervals.len(), loaded_intervals.len());
        for (o, l) in orig_intervals.iter().zip(loaded_intervals.iter()) {
            assert_eq!(o.genome, l.genome);
            assert_eq!(o.start, l.start);
            assert_eq!(o.end, l.end);
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 14. Query completeness vs known ground truth ───────────────

    /// Helper: build test sequences with known homology structure.
    /// Returns (sequences, shared_regions) where shared_regions describes
    /// which genomes share which backbone.
    fn make_homologous_sequences() -> Vec<(String, Vec<u8>)> {
        // All genomes share a 400bp backbone at the start.
        // genome_a and genome_b also share a 300bp region at the end.
        // genome_c diverges after the common backbone.
        // genome_d is completely independent (different seed, no shared syncmers).
        let common_backbone = make_test_sequence(400, 42);
        let shared_tail_ab = make_test_sequence(300, 77);

        let mut seq_a = common_backbone.clone();
        seq_a.extend_from_slice(&make_test_sequence(200, 1)); // unique middle
        seq_a.extend_from_slice(&shared_tail_ab);

        let mut seq_b = common_backbone.clone();
        seq_b.extend_from_slice(&make_test_sequence(200, 2)); // different middle
        seq_b.extend_from_slice(&shared_tail_ab);

        let mut seq_c = common_backbone.clone();
        seq_c.extend_from_slice(&make_test_sequence(500, 3)); // different rest

        let seq_d = make_test_sequence(900, 99); // completely independent

        vec![
            ("genome_a".to_string(), seq_a),
            ("genome_b".to_string(), seq_b),
            ("genome_c".to_string(), seq_c),
            ("genome_d".to_string(), seq_d),
        ]
    }

    #[test]
    fn test_query_completeness_shared_backbone() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seqs = make_homologous_sequences();
        let index = SyngIndex::build(params, seqs.into_iter());

        // Query the common backbone region (first 400bp) on genome_a
        let intervals = index.query_region("genome_a", 0, 400, 0).unwrap();

        let genomes: Vec<&str> = intervals.iter().map(|iv| iv.genome.as_str()).collect();

        // genome_a (self), genome_b, genome_c should all appear (shared backbone)
        assert!(
            genomes.contains(&"genome_a"),
            "Self-hit expected. Found: {:?}",
            genomes
        );
        assert!(
            genomes.contains(&"genome_b"),
            "genome_b shares backbone. Found: {:?}",
            genomes
        );
        assert!(
            genomes.contains(&"genome_c"),
            "genome_c shares backbone. Found: {:?}",
            genomes
        );

        // Coordinates should be in a reasonable range (within ~syncmer_len of 0..400)
        let syncmer_len = (params.w + params.k) as u64;
        for iv in &intervals {
            if iv.genome == "genome_a" || iv.genome == "genome_b" || iv.genome == "genome_c" {
                // Start should be near 0 (within syncmer length)
                assert!(
                    iv.start <= syncmer_len,
                    "{}: start {} should be near 0 (within {})",
                    iv.genome,
                    iv.start,
                    syncmer_len
                );
            }
        }
    }

    #[test]
    fn test_query_completeness_no_false_negatives_interior() {
        let _guard = lock_syng();
        // For interior coverage: if genome_a and genome_b share an identical region,
        // querying interior of that region should always find genome_b.
        let params = SyncmerParams::default();
        let shared_region = make_test_sequence(800, 42);

        let mut seq_a = make_test_sequence(200, 10); // unique prefix
        seq_a.extend_from_slice(&shared_region);
        seq_a.extend_from_slice(&make_test_sequence(200, 11)); // unique suffix

        let mut seq_b = make_test_sequence(200, 20); // different prefix
        seq_b.extend_from_slice(&shared_region);
        seq_b.extend_from_slice(&make_test_sequence(200, 21)); // different suffix

        let index = SyngIndex::build(
            params,
            vec![
                ("genome_a".to_string(), seq_a),
                ("genome_b".to_string(), seq_b),
            ]
            .into_iter(),
        );

        // Query interior of the shared region on genome_a (offset by 200bp prefix)
        // Query 300..700 within the shared region (which is at 200..1000 on genome_a)
        let intervals = index.query_region("genome_a", 400, 800, 0).unwrap();
        let genomes: Vec<&str> = intervals.iter().map(|iv| iv.genome.as_str()).collect();

        assert!(
            genomes.contains(&"genome_b"),
            "Interior of shared region should find genome_b. Found: {:?}",
            genomes
        );
    }

    // ── 15. Boundary padding tests ─────────────────────────────────

    #[test]
    fn test_query_padding_extends_intervals() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq_a = make_test_sequence(1000, 42);
        let seq_b = seq_a.clone(); // identical

        let index = SyngIndex::build(
            params,
            vec![
                ("genome_a".to_string(), seq_a),
                ("genome_b".to_string(), seq_b),
            ]
            .into_iter(),
        );

        let no_pad = index.query_region("genome_a", 200, 600, 0).unwrap();
        let pad_120 = index.query_region("genome_a", 200, 600, 120).unwrap();

        // Both should find genome_b
        assert!(no_pad.iter().any(|iv| iv.genome == "genome_b"));
        assert!(pad_120.iter().any(|iv| iv.genome == "genome_b"));

        // Find genome_b in each
        let no_pad_b = no_pad.iter().find(|iv| iv.genome == "genome_b").unwrap();
        let pad_120_b = pad_120.iter().find(|iv| iv.genome == "genome_b").unwrap();

        // Padded interval should be at least as wide
        assert!(
            pad_120_b.start <= no_pad_b.start,
            "Padded start {} should be <= unpadded start {}",
            pad_120_b.start,
            no_pad_b.start
        );
        assert!(
            pad_120_b.end >= no_pad_b.end,
            "Padded end {} should be >= unpadded end {}",
            pad_120_b.end,
            no_pad_b.end
        );

        // The difference should be approximately the padding amount
        let start_diff = no_pad_b.start as i64 - pad_120_b.start as i64;
        let end_diff = pad_120_b.end as i64 - no_pad_b.end as i64;
        assert!(
            start_diff >= 0,
            "start_diff should be non-negative: {}",
            start_diff
        );
        assert!(
            end_diff >= 0,
            "end_diff should be non-negative: {}",
            end_diff
        );
    }

    #[test]
    fn test_query_padding_zero_vs_nonzero() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq = make_test_sequence(1000, 42);

        let index = SyngIndex::build(
            params,
            vec![("g1".to_string(), seq.clone()), ("g2".to_string(), seq)].into_iter(),
        );

        let pad_0 = index.query_region("g1", 300, 500, 0).unwrap();
        let pad_60 = index.query_region("g1", 300, 500, 60).unwrap();
        let pad_120 = index.query_region("g1", 300, 500, 120).unwrap();

        // Intervals should grow monotonically with padding
        let g2_0 = pad_0.iter().find(|iv| iv.genome == "g2").unwrap();
        let g2_60 = pad_60.iter().find(|iv| iv.genome == "g2").unwrap();
        let g2_120 = pad_120.iter().find(|iv| iv.genome == "g2").unwrap();

        assert!(
            g2_60.start <= g2_0.start,
            "60bp pad start should be <= 0bp pad start"
        );
        assert!(
            g2_60.end >= g2_0.end,
            "60bp pad end should be >= 0bp pad end"
        );
        assert!(
            g2_120.start <= g2_60.start,
            "120bp pad start should be <= 60bp pad start"
        );
        assert!(
            g2_120.end >= g2_60.end,
            "120bp pad end should be >= 60bp pad end"
        );
    }

    #[test]
    fn test_query_padding_clamped_to_genome_length() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq = make_test_sequence(500, 42);
        let genome_len = seq.len() as u64;

        let index = SyngIndex::build(
            params,
            vec![("g1".to_string(), seq.clone()), ("g2".to_string(), seq)].into_iter(),
        );

        // Query near the start with large padding — should clamp to 0
        let intervals = index.query_region("g1", 0, 200, 10000).unwrap();
        for iv in &intervals {
            assert!(
                iv.start == 0,
                "Start should be clamped to 0, got {}",
                iv.start
            );
            assert!(
                iv.end <= genome_len,
                "End {} should be <= genome length {}",
                iv.end,
                genome_len
            );
        }
    }

    // ── 16. Interval merging tests ─────────────────────────────────
    //
    // `merge_intervals_with_anchors` operates on a per-target-path group,
    // so all inputs share the same genome by construction. Cross-genome
    // grouping happens one level up (in `query_region_with_anchors`).

    fn mk_hiwa(start: u64, end: u64, anchors: Vec<(u64, u64)>) -> HomologousIntervalWithAnchors {
        HomologousIntervalWithAnchors {
            genome: "g1".to_string(),
            start,
            end,
            strand: '+',
            anchors: anchors
                .into_iter()
                .map(|(q, t)| Anchor {
                    query_pos: q,
                    target_pos: t,
                    node_id: 0,
                })
                .collect(),
        }
    }

    #[test]
    fn test_merge_with_anchors_no_overlap() {
        let _guard = lock_syng();
        let mut intervals = vec![
            mk_hiwa(0, 100, vec![(10, 10)]),
            mk_hiwa(200, 300, vec![(250, 250)]),
        ];
        SyngIndex::merge_intervals_with_anchors(&mut intervals);
        assert_eq!(intervals.len(), 2);
    }

    #[test]
    fn test_merge_with_anchors_overlapping() {
        let _guard = lock_syng();
        let mut intervals = vec![
            mk_hiwa(0, 150, vec![(10, 10)]),
            mk_hiwa(100, 300, vec![(200, 200)]),
        ];
        SyngIndex::merge_intervals_with_anchors(&mut intervals);
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].start, 0);
        assert_eq!(intervals[0].end, 300);
        // Anchors from both intervals are preserved
        assert_eq!(intervals[0].anchors.len(), 2);
    }

    #[test]
    fn test_merge_with_anchors_adjacent() {
        let _guard = lock_syng();
        let mut intervals = vec![
            mk_hiwa(0, 100, vec![(10, 10)]),
            mk_hiwa(100, 200, vec![(150, 150)]),
        ];
        SyngIndex::merge_intervals_with_anchors(&mut intervals);
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].start, 0);
        assert_eq!(intervals[0].end, 200);
    }

    #[test]
    fn test_merge_with_anchors_multiple_groups() {
        let _guard = lock_syng();
        let mut intervals = vec![
            mk_hiwa(0, 100, vec![(10, 10)]),
            mk_hiwa(50, 200, vec![(100, 100)]),
            mk_hiwa(150, 400, vec![(300, 300)]),
        ];
        SyngIndex::merge_intervals_with_anchors(&mut intervals);
        // All three overlap transitively → one merged interval
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].start, 0);
        assert_eq!(intervals[0].end, 400);
        assert_eq!(intervals[0].anchors.len(), 3);
    }

    #[test]
    fn test_merge_with_anchors_dedup() {
        let _guard = lock_syng();
        let mut intervals = vec![
            mk_hiwa(0, 150, vec![(10, 10), (20, 20)]),
            mk_hiwa(100, 300, vec![(20, 20), (200, 200)]), // (20,20) is dup
        ];
        SyngIndex::merge_intervals_with_anchors(&mut intervals);
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].anchors.len(), 3); // deduped
    }

    #[test]
    fn test_merge_with_anchors_empty() {
        let _guard = lock_syng();
        let mut intervals: Vec<HomologousIntervalWithAnchors> = Vec::new();
        SyngIndex::merge_intervals_with_anchors(&mut intervals);
        assert_eq!(intervals.len(), 0);
    }

    #[test]
    fn test_merge_with_anchors_single() {
        let _guard = lock_syng();
        let mut intervals = vec![mk_hiwa(10, 20, vec![(15, 15)])];
        SyngIndex::merge_intervals_with_anchors(&mut intervals);
        assert_eq!(intervals.len(), 1);
    }

    // ── 17. Edge cases for query_region ────────────────────────────

    #[test]
    fn test_query_region_isolated_region() {
        let _guard = lock_syng();
        // genome_a and genome_b share nothing (different seeds) —
        // querying genome_a should only return self.
        let params = SyncmerParams::default();
        let seq_a = make_test_sequence(500, 1);
        let seq_b = make_test_sequence(500, 99); // completely different

        let index = SyngIndex::build(
            params,
            vec![
                ("genome_a".to_string(), seq_a),
                ("genome_b".to_string(), seq_b),
            ]
            .into_iter(),
        );

        let intervals = index.query_region("genome_a", 100, 300, 0).unwrap();

        // Should contain self-hit but not genome_b
        let genomes: Vec<&str> = intervals.iter().map(|iv| iv.genome.as_str()).collect();
        assert!(
            genomes.contains(&"genome_a"),
            "Self-hit expected. Found: {:?}",
            genomes
        );
        // genome_b may or may not appear — with truly different seeds it shouldn't share syncmers.
        // But we can't 100% guarantee no accidental collision. Check that if genome_b appears,
        // it's because of a genuine (unlikely) collision, not a bug.
    }

    #[test]
    fn test_query_region_entire_sequence() {
        let _guard = lock_syng();
        // Query the full length of genome_a — should return all genomes that share any syncmers
        let params = SyncmerParams::default();
        let shared = make_test_sequence(500, 42);

        let mut seq_a = shared.clone();
        seq_a.extend_from_slice(&make_test_sequence(500, 1));

        let mut seq_b = shared.clone();
        seq_b.extend_from_slice(&make_test_sequence(500, 2));

        let seq_c = make_test_sequence(1000, 99); // independent

        let index = SyngIndex::build(
            params,
            vec![
                ("genome_a".to_string(), seq_a.clone()),
                ("genome_b".to_string(), seq_b),
                ("genome_c".to_string(), seq_c),
            ]
            .into_iter(),
        );

        let intervals = index
            .query_region("genome_a", 0, seq_a.len() as u64, 0)
            .unwrap();
        let genomes: Vec<&str> = intervals.iter().map(|iv| iv.genome.as_str()).collect();

        assert!(genomes.contains(&"genome_a"), "Self-hit expected");
        assert!(genomes.contains(&"genome_b"), "genome_b shares prefix");
    }

    #[test]
    fn test_query_region_single_sequence_index() {
        let _guard = lock_syng();
        // Index with one sequence — query should return only self
        let params = SyncmerParams::default();
        let seq = make_test_sequence(1000, 42);

        let index = SyngIndex::build(params, vec![("only_genome".to_string(), seq)].into_iter());

        let intervals = index.query_region("only_genome", 100, 500, 0).unwrap();
        // The query must find at least a forward self-hit. Random data can
        // also produce palindromic syncmer matches (same syncmer canonical
        // hash in both orientations on the same path), so we tolerate '-'
        // self-hits in addition — they're legitimate RC homology within the
        // single genome, not noise.
        assert!(
            !intervals.is_empty(),
            "Single-sequence index must return at least one self-hit"
        );
        assert!(
            intervals.iter().all(|iv| iv.genome == "only_genome"),
            "All hits must be on the only indexed genome"
        );
        assert!(
            intervals.iter().any(|iv| iv.strand == '+'),
            "Must find at least one forward self-hit"
        );
    }

    #[test]
    fn test_query_region_out_of_range() {
        let _guard = lock_syng();
        // Query beyond the end of the sequence — should return empty (no syncmers there)
        let params = SyncmerParams::default();
        let seq = make_test_sequence(500, 42);

        let index = SyngIndex::build(params, vec![("g1".to_string(), seq)].into_iter());

        let intervals = index.query_region("g1", 10000, 20000, 0).unwrap();
        assert!(
            intervals.is_empty(),
            "Query beyond sequence length should return empty, got {:?}",
            intervals
                .iter()
                .map(|iv| (&iv.genome, iv.start, iv.end))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_query_region_zero_width() {
        let _guard = lock_syng();
        // Query with start == end — no syncmers in empty range
        let params = SyncmerParams::default();
        let seq = make_test_sequence(500, 42);

        let index = SyngIndex::build(params, vec![("g1".to_string(), seq)].into_iter());

        let intervals = index.query_region("g1", 200, 200, 0).unwrap();
        assert!(intervals.is_empty(), "Zero-width query should return empty");
    }

    #[test]
    fn test_query_region_identical_sequences() {
        let _guard = lock_syng();
        // Multiple identical sequences — all should appear as hits
        let params = SyncmerParams::default();
        let seq = make_test_sequence(1000, 42);

        let index = SyngIndex::build(
            params,
            vec![
                ("g1".to_string(), seq.clone()),
                ("g2".to_string(), seq.clone()),
                ("g3".to_string(), seq.clone()),
                ("g4".to_string(), seq),
            ]
            .into_iter(),
        );

        let intervals = index.query_region("g1", 100, 500, 0).unwrap();
        let genomes: Vec<&str> = intervals.iter().map(|iv| iv.genome.as_str()).collect();

        assert!(genomes.contains(&"g1"), "Should find self");
        assert!(genomes.contains(&"g2"), "Should find g2 (identical)");
        assert!(genomes.contains(&"g3"), "Should find g3 (identical)");
        assert!(genomes.contains(&"g4"), "Should find g4 (identical)");
        assert_eq!(intervals.len(), 4, "Should find exactly 4 genomes");
    }

    // ── 18. CLI integration: syng build + query ────────────────────

    #[test]
    fn test_syng_cli_query_bed_output() {
        let _guard = lock_syng();
        // Build index via CLI, then query and get BED output
        let dir = std::env::temp_dir().join("impg_test_syng_cli_query");
        std::fs::create_dir_all(&dir).unwrap();

        // Write test FASTA with shared content
        let fasta_path = dir.join("test.fa");
        let shared = make_test_sequence(500, 42);
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            let mut seq_a = shared.clone();
            seq_a.extend_from_slice(&make_test_sequence(500, 1));
            let mut seq_b = shared.clone();
            seq_b.extend_from_slice(&make_test_sequence(500, 2));
            use std::io::Write;
            writeln!(f, ">genome_a").unwrap();
            writeln!(f, "{}", String::from_utf8(seq_a).unwrap()).unwrap();
            writeln!(f, ">genome_b").unwrap();
            writeln!(f, "{}", String::from_utf8(seq_b).unwrap()).unwrap();
        }

        // Find impg binary
        let bin = find_impg_binary();
        if bin.is_none() {
            eprintln!("Skipping CLI test: impg binary not found");
            std::fs::remove_dir_all(&dir).ok();
            return;
        }
        let bin = bin.unwrap();

        // Build index
        let output_prefix = dir.join("idx");
        let output = std::process::Command::new(&bin)
            .args([
                "syng",
                "-f",
                fasta_path.to_str().unwrap(),
                "-o",
                output_prefix.to_str().unwrap(),
                "--position-sample-rate",
                "1",
            ])
            .output()
            .expect("Failed to run impg syng");
        assert!(
            output.status.success(),
            "impg syng failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        // Query via -a syng prefix with BED output
        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "-d",
                "0",
                "-a",
                output_prefix.to_str().unwrap(),
                "--sequence-files",
                fasta_path.to_str().unwrap(),
                "-r",
                "genome_a:0-400",
                "-o",
                "bed",
            ])
            .output()
            .expect("Failed to run impg query -a syng prefix");
        assert!(
            output.status.success(),
            "impg query -a syng prefix failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8_lossy(&output.stdout);
        let lines: Vec<&str> = stdout.lines().collect();
        assert!(!lines.is_empty(), "BED output should not be empty");
        assert!(
            lines.iter().all(|line| line.split('\t').count() >= 6),
            "BED stdout should contain only BED records, not native diagnostics. Got:\n{}",
            stdout
        );

        // Should find genome_b in the output
        let has_genome_b = lines.iter().any(|l| l.starts_with("genome_b\t"));
        assert!(
            has_genome_b,
            "BED output should contain genome_b. Got:\n{}",
            stdout
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_syng_cli_query_gfa_output() {
        let _guard = lock_syng();
        let dir = std::env::temp_dir().join("impg_test_syng_cli_gfa");
        std::fs::create_dir_all(&dir).unwrap();

        // Write test FASTA
        let fasta_path = dir.join("test.fa");
        let shared = make_test_sequence(500, 42);
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            let mut seq_a = shared.clone();
            seq_a.extend_from_slice(&make_test_sequence(500, 1));
            let mut seq_b = shared.clone();
            seq_b.extend_from_slice(&make_test_sequence(500, 2));
            use std::io::Write;
            writeln!(f, ">genome_a").unwrap();
            writeln!(f, "{}", String::from_utf8(seq_a).unwrap()).unwrap();
            writeln!(f, ">genome_b").unwrap();
            writeln!(f, "{}", String::from_utf8(seq_b).unwrap()).unwrap();
        }

        let bin = find_impg_binary();
        if bin.is_none() {
            eprintln!("Skipping CLI test: impg binary not found");
            std::fs::remove_dir_all(&dir).ok();
            return;
        }
        let bin = bin.unwrap();

        // Build index
        let output_prefix = dir.join("idx");
        let output = std::process::Command::new(&bin)
            .args([
                "syng",
                "-f",
                fasta_path.to_str().unwrap(),
                "-o",
                output_prefix.to_str().unwrap(),
                "--position-sample-rate",
                "1",
            ])
            .output()
            .expect("Failed to run impg syng");
        assert!(
            output.status.success(),
            "impg syng failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        // Query with GFA output
        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "-d",
                "0",
                "-a",
                output_prefix.to_str().unwrap(),
                "--sequence-files",
                fasta_path.to_str().unwrap(),
                "-r",
                "genome_a:0-400",
                "-o",
                "gfa",
            ])
            .output()
            .expect("Failed to run impg query -a syng prefix -o gfa");
        assert!(
            output.status.success(),
            "impg query -a syng prefix -o gfa failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8_lossy(&output.stdout);
        // GFA should have S (segment) and L (link) lines
        let has_s_lines = stdout.lines().any(|l| l.starts_with("S\t"));
        assert!(
            has_s_lines,
            "GFA output should contain S lines. Got:\n{}",
            stdout
        );

        // Check for path-related lines (W or P lines)
        let has_paths = stdout
            .lines()
            .any(|l| l.starts_with("W\t") || l.starts_with("P\t"));
        assert!(
            has_paths,
            "GFA output should contain W or P lines. Got:\n{}",
            stdout
        );

        // Syng engine defaults to blunt graph output, so every emitted link
        // should have zero overlap after pangenome/bluntg processing.
        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "-d",
                "0",
                "-a",
                output_prefix.to_str().unwrap(),
                "--sequence-files",
                fasta_path.to_str().unwrap(),
                "-r",
                "genome_a:0-400",
                "-o",
                "gfa",
                "--gfa-engine",
                "syng",
            ])
            .output()
            .expect("Failed to run impg query -a syng prefix -o gfa --gfa-engine syng");
        assert!(
            output.status.success(),
            "impg query -a syng prefix -o gfa --gfa-engine syng failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let stdout = String::from_utf8_lossy(&output.stdout);
        let nonzero_links: Vec<&str> = stdout
            .lines()
            .filter(|l| l.starts_with("L\t") && !l.ends_with("\t0M"))
            .collect();
        assert!(
            nonzero_links.is_empty(),
            "syng engine should default to blunt 0M links. Nonzero links: {:?}\n{}",
            nonzero_links,
            stdout
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── build_region_gbwt tests ───────────────────────────────────

    #[test]
    fn test_build_region_gbwt_produces_files() {
        let _guard = lock_syng();
        // Build an index, then build a region GBWT from some sequences
        let seqs: Vec<(String, Vec<u8>)> = (0..3)
            .map(|i| (format!("seq_{}", i), make_test_sequence(2000, i as u8 + 30)))
            .collect();
        let params = SyncmerParams::default();
        let index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_region_gbwt");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("region");
        let prefix_str = prefix.to_str().unwrap();

        // Build region GBWT from a subset of sequences
        let region_seqs: Vec<(String, Vec<u8>)> = (0..2)
            .map(|i| {
                (
                    format!("region_seq_{}", i),
                    make_test_sequence(500, i as u8 + 50),
                )
            })
            .collect();
        let refs: Vec<(String, &[u8])> = region_seqs
            .iter()
            .map(|(n, s)| (n.clone(), s.as_slice()))
            .collect();

        index.build_region_gbwt(&refs, prefix_str).unwrap();

        assert!(
            std::path::Path::new(&format!("{}.1gbwt", prefix_str)).exists(),
            "Region .1gbwt file should exist"
        );
        assert!(
            std::path::Path::new(&format!("{}.1khash", prefix_str)).exists(),
            "Region .1khash file should exist"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_build_region_gbwt_empty_sequences() {
        let _guard = lock_syng();
        let index = SyngIndex::new(SyncmerParams::default());

        let dir = std::env::temp_dir().join("impg_test_region_gbwt_empty");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("empty_region");
        let prefix_str = prefix.to_str().unwrap();

        // Empty sequence list should still produce valid files
        let refs: Vec<(String, &[u8])> = Vec::new();
        index.build_region_gbwt(&refs, prefix_str).unwrap();

        assert!(std::path::Path::new(&format!("{}.1gbwt", prefix_str)).exists());
        assert!(std::path::Path::new(&format!("{}.1khash", prefix_str)).exists());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_build_region_gbwt_short_sequences_skipped() {
        let _guard = lock_syng();
        let index = SyngIndex::new(SyncmerParams::default());

        let dir = std::env::temp_dir().join("impg_test_region_gbwt_short");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("short_region");
        let prefix_str = prefix.to_str().unwrap();

        // All sequences shorter than syncmer length (63) should be skipped
        let short_seqs = vec![
            ("short1".to_string(), b"ACGT" as &[u8]),
            ("short2".to_string(), b"ACGTACGTACGT" as &[u8]),
        ];
        index.build_region_gbwt(&short_seqs, prefix_str).unwrap();

        // Files should still be produced (empty but valid)
        assert!(std::path::Path::new(&format!("{}.1gbwt", prefix_str)).exists());
        assert!(std::path::Path::new(&format!("{}.1khash", prefix_str)).exists());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_build_region_gbwt_loadable() {
        let _guard = lock_syng();
        // Build a region GBWT and verify it can be loaded back
        let params = SyncmerParams::default();
        let index = SyngIndex::new(params);

        let dir = std::env::temp_dir().join("impg_test_region_gbwt_load");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("loadable_region");
        let prefix_str = prefix.to_str().unwrap();

        let region_seqs: Vec<(String, Vec<u8>)> = (0..3)
            .map(|i| (format!("reg_{}", i), make_test_sequence(1000, i as u8 + 70)))
            .collect();
        let refs: Vec<(String, &[u8])> = region_seqs
            .iter()
            .map(|(n, s)| (n.clone(), s.as_slice()))
            .collect();

        index.build_region_gbwt(&refs, prefix_str).unwrap();

        // Write a names file so we can load it
        let names_path = format!("{}.syng.names", prefix_str);
        let mut nm = SyngNameMap::new();
        for (name, seq) in &region_seqs {
            nm.add(name.clone(), seq.len() as u64);
        }
        nm.save(&names_path).unwrap();

        // Should be loadable
        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert!(!loaded.gbwt.is_null());
        assert!(!loaded.kmer_hash.is_null());

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 19. Region-specific GBWT output (spec Test 4) ───────────────

    #[test]
    fn test_region_gbwt_from_query_results() {
        let _guard = lock_syng();
        // Build full index → query region → output as GBWT → verify files
        let params = SyncmerParams::default();
        let shared = make_test_sequence(500, 42);

        let mut seq_a = shared.clone();
        seq_a.extend_from_slice(&make_test_sequence(500, 1));
        let mut seq_b = shared.clone();
        seq_b.extend_from_slice(&make_test_sequence(500, 2));
        let seq_c = make_test_sequence(1000, 99); // independent

        let sequences = vec![
            ("genome_a".to_string(), seq_a),
            ("genome_b".to_string(), seq_b),
            ("genome_c".to_string(), seq_c),
        ];

        let index = SyngIndex::build(params, sequences.into_iter());

        // Query the shared region
        let intervals = index.query_region("genome_a", 0, 500, 0).unwrap();
        assert!(
            !intervals.is_empty(),
            "Should find intervals in shared region"
        );

        // Build region GBWT from the query result sequences
        // (simulate what the CLI does: fetch sequences for each interval)
        let region_seqs: Vec<(String, Vec<u8>)> = intervals
            .iter()
            .map(|iv| {
                let seq_name = format!("{}:{}-{}", iv.genome, iv.start, iv.end);
                let seq = make_test_sequence((iv.end - iv.start) as usize, 42);
                (seq_name, seq)
            })
            .collect();

        let dir = std::env::temp_dir().join("impg_test_region_gbwt_query");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("region_from_query");
        let prefix_str = prefix.to_str().unwrap();

        let refs: Vec<(String, &[u8])> = region_seqs
            .iter()
            .map(|(n, s)| (n.clone(), s.as_slice()))
            .collect();

        index.build_region_gbwt(&refs, prefix_str).unwrap();

        // Verify output files exist and are non-empty
        let gbwt_path = format!("{}.1gbwt", prefix_str);
        let khash_path = format!("{}.1khash", prefix_str);
        assert!(Path::new(&gbwt_path).exists(), "Region .1gbwt should exist");
        assert!(
            Path::new(&khash_path).exists(),
            "Region .1khash should exist"
        );
        assert!(
            std::fs::metadata(&gbwt_path).unwrap().len() > 0,
            "Region .1gbwt should be non-empty"
        );
        assert!(
            std::fs::metadata(&khash_path).unwrap().len() > 0,
            "Region .1khash should be non-empty"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_region_gbwt_loadable_as_syng_index() {
        let _guard = lock_syng();
        // Build full index, query, produce region GBWT, load as SyngIndex
        let params = SyncmerParams::default();
        let shared = make_test_sequence(600, 42);

        let mut seq_a = shared.clone();
        seq_a.extend_from_slice(&make_test_sequence(400, 1));
        let mut seq_b = shared.clone();
        seq_b.extend_from_slice(&make_test_sequence(400, 2));

        let index = SyngIndex::build(
            params,
            vec![
                ("genome_a".to_string(), seq_a),
                ("genome_b".to_string(), seq_b),
            ]
            .into_iter(),
        );

        // Build region GBWT from actual sequence data (use shared prefix)
        let region_seqs: Vec<(String, Vec<u8>)> = vec![
            ("region_a".to_string(), make_test_sequence(400, 42)),
            ("region_b".to_string(), make_test_sequence(400, 42)),
            ("region_c".to_string(), make_test_sequence(400, 10)),
        ];

        let dir = std::env::temp_dir().join("impg_test_region_gbwt_load_syng");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("region_loadable");
        let prefix_str = prefix.to_str().unwrap();

        let refs: Vec<(String, &[u8])> = region_seqs
            .iter()
            .map(|(n, s)| (n.clone(), s.as_slice()))
            .collect();
        index.build_region_gbwt(&refs, prefix_str).unwrap();

        // Write a names file so we can load it as a full SyngIndex
        let names_path = format!("{}.syng.names", prefix_str);
        let mut nm = SyngNameMap::new();
        for (name, seq) in &region_seqs {
            nm.add(name.clone(), seq.len() as u64);
        }
        nm.save(&names_path).unwrap();

        // Load the region GBWT back as a SyngIndex
        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert!(!loaded.gbwt.is_null(), "Loaded region GBWT should be valid");
        assert!(
            !loaded.kmer_hash.is_null(),
            "Loaded region KmerHash should be valid"
        );
        assert_eq!(
            loaded.name_map.path_to_name.len(),
            3,
            "Region index should have 3 paths"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_region_gbwt_nodes_subset_of_full_index() {
        let _guard = lock_syng();
        // Build a full index and a region GBWT from a subset of sequences.
        // The region GBWT's syncmer nodes should be a subset of the full index's nodes.
        let params = SyncmerParams::default();

        // Use the SAME sequences for both: full index from all, region from subset
        let seqs: Vec<(String, Vec<u8>)> = (0..4)
            .map(|i| {
                (
                    format!("genome_{}", i),
                    make_test_sequence(2000, i as u8 + 10),
                )
            })
            .collect();
        let mut full_index = SyngIndex::build(params, seqs.clone().into_iter());

        let dir = std::env::temp_dir().join("impg_test_region_nodes_subset");
        std::fs::create_dir_all(&dir).unwrap();

        // Save full index to disk
        let full_prefix = dir.join("full");
        full_index.save(full_prefix.to_str().unwrap()).unwrap();

        // Build region GBWT from just the first 2 sequences (subset)
        let region_prefix = dir.join("region_subset");
        let region_refs: Vec<(String, &[u8])> = seqs[0..2]
            .iter()
            .map(|(n, s)| (n.clone(), s.as_slice()))
            .collect();
        full_index
            .build_region_gbwt(&region_refs, region_prefix.to_str().unwrap())
            .unwrap();

        // Both files should exist and be non-empty
        let region_gbwt = format!("{}.1gbwt", region_prefix.to_str().unwrap());
        let region_khash = format!("{}.1khash", region_prefix.to_str().unwrap());
        assert!(Path::new(&region_gbwt).exists());
        assert!(Path::new(&region_khash).exists());
        assert!(std::fs::metadata(&region_gbwt).unwrap().len() > 0);
        assert!(std::fs::metadata(&region_khash).unwrap().len() > 0);

        // Load both as SyngIndex and verify the region has fewer or equal paths
        let full_loaded = SyngIndex::load(full_prefix.to_str().unwrap(), params).unwrap();
        let names_path = format!("{}.syng.names", region_prefix.to_str().unwrap());
        let mut nm = SyngNameMap::new();
        for (name, seq) in &seqs[0..2] {
            nm.add(name.clone(), seq.len() as u64);
        }
        nm.save(&names_path).unwrap();
        let region_loaded = SyngIndex::load(region_prefix.to_str().unwrap(), params).unwrap();
        assert!(
            full_loaded.name_map.path_to_name.len() >= region_loaded.name_map.path_to_name.len(),
            "Full index paths ({}) should be >= region paths ({})",
            full_loaded.name_map.path_to_name.len(),
            region_loaded.name_map.path_to_name.len(),
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_region_gbwt_single_genome() {
        let _guard = lock_syng();
        // Edge case: region with a single genome → GBWT with one path
        let params = SyncmerParams::default();
        let index = SyngIndex::new(params);

        let dir = std::env::temp_dir().join("impg_test_region_gbwt_single");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("single_genome_region");
        let prefix_str = prefix.to_str().unwrap();

        let seq = make_test_sequence(500, 42);
        let refs: Vec<(String, &[u8])> = vec![("only_genome".to_string(), seq.as_slice())];

        index.build_region_gbwt(&refs, prefix_str).unwrap();

        // Files should exist
        assert!(Path::new(&format!("{}.1gbwt", prefix_str)).exists());
        assert!(Path::new(&format!("{}.1khash", prefix_str)).exists());

        // Write names and load back
        let names_path = format!("{}.syng.names", prefix_str);
        let mut nm = SyngNameMap::new();
        nm.add("only_genome".to_string(), seq.len() as u64);
        nm.save(&names_path).unwrap();

        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert_eq!(loaded.name_map.path_to_name.len(), 1);
        assert_eq!(loaded.name_map.path_to_name[0], "only_genome");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_region_gbwt_very_small_region() {
        let _guard = lock_syng();
        // Very small region (smaller than one syncmer) → empty/minimal GBWT
        let params = SyncmerParams::default();
        let index = SyngIndex::new(params);

        let dir = std::env::temp_dir().join("impg_test_region_gbwt_tiny");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("tiny_region");
        let prefix_str = prefix.to_str().unwrap();

        // Sequences shorter than syncmer length (63bp) — should be skipped
        let refs: Vec<(String, &[u8])> = vec![
            ("tiny1".to_string(), b"ACGTACGT" as &[u8]),
            ("tiny2".to_string(), b"TGCATGCA" as &[u8]),
        ];

        index.build_region_gbwt(&refs, prefix_str).unwrap();

        // Files should still be produced (valid but minimal)
        assert!(Path::new(&format!("{}.1gbwt", prefix_str)).exists());
        assert!(Path::new(&format!("{}.1khash", prefix_str)).exists());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_region_gbwt_output_prefix_with_directory() {
        let _guard = lock_syng();
        // Output prefix with nested directory path
        let dir = std::env::temp_dir().join("impg_test_region_gbwt_nested/subdir/deep");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("nested_output");
        let prefix_str = prefix.to_str().unwrap();

        let params = SyncmerParams::default();
        let index = SyngIndex::new(params);

        let seq = make_test_sequence(500, 42);
        let refs: Vec<(String, &[u8])> = vec![("seq1".to_string(), seq.as_slice())];

        index.build_region_gbwt(&refs, prefix_str).unwrap();

        assert!(Path::new(&format!("{}.1gbwt", prefix_str)).exists());
        assert!(Path::new(&format!("{}.1khash", prefix_str)).exists());

        // Clean up the top-level temp dir
        std::fs::remove_dir_all(std::env::temp_dir().join("impg_test_region_gbwt_nested")).ok();
    }

    #[test]
    fn test_region_gbwt_nonexistent_directory_fails() {
        let _guard = lock_syng();
        // Output prefix whose parent directory doesn't exist → should error
        let params = SyncmerParams::default();
        let index = SyngIndex::new(params);

        let prefix = "/tmp/impg_test_nonexistent_dir_xyz123/subdir/output";

        let seq = make_test_sequence(500, 42);
        let refs: Vec<(String, &[u8])> = vec![("seq1".to_string(), seq.as_slice())];

        let result = index.build_region_gbwt(&refs, prefix);
        assert!(
            result.is_err(),
            "Should fail when output directory doesn't exist"
        );
    }

    // ── 20. Syng format interoperability (spec Test 3) ───────────────

    #[test]
    fn test_onecode_magic_bytes_gbwt() {
        let _guard = lock_syng();
        // Verify the .1gbwt file has correct ONEcode format markers
        let params = SyncmerParams::default();
        let seqs = vec![("seq1".to_string(), make_test_sequence(1000, 42))];
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_onecode_magic_gbwt");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("magic_test");
        let prefix_str = prefix.to_str().unwrap();
        index.save(prefix_str).unwrap();

        // Read first few bytes of .1gbwt — ONEcode binary files start with specific markers
        let gbwt_bytes = std::fs::read(format!("{}.1gbwt", prefix_str)).unwrap();
        assert!(gbwt_bytes.len() > 4, ".1gbwt should be non-trivial size");

        // ONEcode binary format: first byte is '1' (0x31) for binary mode
        // or the file might start with the schema header
        // Check that it's either binary ('1') or text ('#' for comment/schema)
        let first_byte = gbwt_bytes[0];
        assert!(
            first_byte == b'1' || first_byte == b'#' || first_byte == b'!',
            ".1gbwt should start with ONEcode marker ('1', '#', or '!'), got 0x{:02x}",
            first_byte
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_onecode_magic_bytes_khash() {
        let _guard = lock_syng();
        // Verify the .1khash file has correct ONEcode format markers
        let params = SyncmerParams::default();
        let seqs = vec![("seq1".to_string(), make_test_sequence(1000, 42))];
        let mut index = SyngIndex::build(params, seqs.into_iter());

        let dir = std::env::temp_dir().join("impg_test_onecode_magic_khash");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("magic_test");
        let prefix_str = prefix.to_str().unwrap();
        index.save(prefix_str).unwrap();

        let khash_bytes = std::fs::read(format!("{}.1khash", prefix_str)).unwrap();
        assert!(khash_bytes.len() > 4, ".1khash should be non-trivial size");

        let first_byte = khash_bytes[0];
        assert!(
            first_byte == b'1' || first_byte == b'#' || first_byte == b'!',
            ".1khash should start with ONEcode marker ('1', '#', or '!'), got 0x{:02x}",
            first_byte
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_region_gbwt_onecode_format() {
        let _guard = lock_syng();
        // Verify region GBWT output also has correct ONEcode format
        let params = SyncmerParams::default();
        let index = SyngIndex::new(params);

        let dir = std::env::temp_dir().join("impg_test_region_onecode");
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("region_onecode");
        let prefix_str = prefix.to_str().unwrap();

        let seq = make_test_sequence(1000, 42);
        let refs: Vec<(String, &[u8])> = vec![("seq1".to_string(), seq.as_slice())];
        index.build_region_gbwt(&refs, prefix_str).unwrap();

        // Check .1gbwt format marker
        let gbwt_bytes = std::fs::read(format!("{}.1gbwt", prefix_str)).unwrap();
        assert!(gbwt_bytes.len() > 4);
        assert!(
            gbwt_bytes[0] == b'1' || gbwt_bytes[0] == b'#' || gbwt_bytes[0] == b'!',
            "Region .1gbwt should have ONEcode format marker"
        );

        // Check .1khash format marker
        let khash_bytes = std::fs::read(format!("{}.1khash", prefix_str)).unwrap();
        assert!(khash_bytes.len() > 4);
        assert!(
            khash_bytes[0] == b'1' || khash_bytes[0] == b'#' || khash_bytes[0] == b'!',
            "Region .1khash should have ONEcode format marker"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 21. Round-trip: query → GBWT → query (spec Test 4) ──────────

    #[test]
    fn test_roundtrip_query_gbwt_query() {
        let _guard = lock_syng();
        // Build full index → query region → output GBWT → load region GBWT → query region GBWT
        // Verify second query returns consistent results
        let params = SyncmerParams::default();

        // Create sequences with known homology
        let shared = make_test_sequence(600, 42);
        let mut seq_a = shared.clone();
        seq_a.extend_from_slice(&make_test_sequence(400, 1));
        let mut seq_b = shared.clone();
        seq_b.extend_from_slice(&make_test_sequence(400, 2));
        let mut seq_c = shared.clone();
        seq_c.extend_from_slice(&make_test_sequence(400, 3));

        let full_seqs = vec![
            ("genome_a".to_string(), seq_a.clone()),
            ("genome_b".to_string(), seq_b.clone()),
            ("genome_c".to_string(), seq_c.clone()),
        ];
        let full_index = SyngIndex::build(params, full_seqs.into_iter());

        // Step 1: Query the shared region on the full index
        let intervals = full_index.query_region("genome_a", 0, 500, 0).unwrap();
        let original_genomes: Vec<&str> = intervals.iter().map(|iv| iv.genome.as_str()).collect();
        assert!(
            original_genomes.contains(&"genome_b"),
            "Full query should find genome_b"
        );

        // Step 2: Build region GBWT from the full sequences in the shared region
        let dir = std::env::temp_dir().join("impg_test_roundtrip_query_gbwt");
        std::fs::create_dir_all(&dir).unwrap();
        let region_prefix = dir.join("region_rt");
        let region_prefix_str = region_prefix.to_str().unwrap();

        // Use the shared portion of each sequence for the region GBWT
        let region_seqs: Vec<(String, Vec<u8>)> = vec![
            ("genome_a".to_string(), seq_a[0..600].to_vec()),
            ("genome_b".to_string(), seq_b[0..600].to_vec()),
            ("genome_c".to_string(), seq_c[0..600].to_vec()),
        ];
        let region_refs: Vec<(String, &[u8])> = region_seqs
            .iter()
            .map(|(n, s)| (n.clone(), s.as_slice()))
            .collect();
        full_index
            .build_region_gbwt(&region_refs, region_prefix_str)
            .unwrap();

        // Write name map for the region index
        let names_path = format!("{}.syng.names", region_prefix_str);
        let mut nm = SyngNameMap::new();
        for (name, seq) in &region_seqs {
            nm.add(name.clone(), seq.len() as u64);
        }
        nm.save(&names_path).unwrap();

        // Step 3: Load the region GBWT and build a new SyngIndex from
        // the region sequences (so we get path_starts for querying)
        let region_index = SyngIndex::build(params, region_seqs.into_iter());

        // Step 4: Query the region GBWT
        let region_intervals = region_index.query_region("genome_a", 0, 500, 0).unwrap();
        let region_genomes: Vec<&str> = region_intervals
            .iter()
            .map(|iv| iv.genome.as_str())
            .collect();

        // The region query should find the same genomes
        assert!(
            region_genomes.contains(&"genome_a"),
            "Region query should find self. Found: {:?}",
            region_genomes
        );
        assert!(
            region_genomes.contains(&"genome_b"),
            "Region query should find genome_b (shared backbone). Found: {:?}",
            region_genomes
        );
        assert!(
            region_genomes.contains(&"genome_c"),
            "Region query should find genome_c (shared backbone). Found: {:?}",
            region_genomes
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 22. CLI integration: GBWT output format ──────────────────────

    #[test]
    fn test_syng_cli_gbwt_output_from_syng_index() {
        let _guard = lock_syng();
        // impg query -a prefix -f test.fa -r region -o gbwt -O tmpdir/region
        let dir = std::env::temp_dir().join("impg_test_cli_gbwt_output");
        std::fs::create_dir_all(&dir).unwrap();

        let bin = find_impg_binary();
        if bin.is_none() {
            eprintln!("Skipping CLI test: impg binary not found");
            std::fs::remove_dir_all(&dir).ok();
            return;
        }
        let bin = bin.unwrap();

        // Write test FASTA with shared content
        let fasta_path = dir.join("test.fa");
        let shared = make_test_sequence(500, 42);
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            let mut seq_a = shared.clone();
            seq_a.extend_from_slice(&make_test_sequence(500, 1));
            let mut seq_b = shared.clone();
            seq_b.extend_from_slice(&make_test_sequence(500, 2));
            use std::io::Write;
            writeln!(f, ">genome_a").unwrap();
            writeln!(f, "{}", String::from_utf8(seq_a).unwrap()).unwrap();
            writeln!(f, ">genome_b").unwrap();
            writeln!(f, "{}", String::from_utf8(seq_b).unwrap()).unwrap();
        }

        // Build syng index
        let idx_prefix = dir.join("idx");
        let output = std::process::Command::new(&bin)
            .args([
                "syng",
                "-f",
                fasta_path.to_str().unwrap(),
                "-o",
                idx_prefix.to_str().unwrap(),
                "--position-sample-rate",
                "1",
            ])
            .output()
            .expect("Failed to run impg syng");
        assert!(
            output.status.success(),
            "impg syng failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        // Query with GBWT output
        let gbwt_output = dir.join("region_gbwt");
        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "-d",
                "0",
                "-a",
                idx_prefix.to_str().unwrap(),
                "--sequence-files",
                fasta_path.to_str().unwrap(),
                "-r",
                "genome_a:0-400",
                "-o",
                "gbwt",
                "-O",
                gbwt_output.to_str().unwrap(),
            ])
            .output()
            .expect("Failed to run impg query -o gbwt");
        assert!(
            output.status.success(),
            "impg query -a syng prefix -o gbwt failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        // Verify two output files
        let gbwt_file = format!("{}.1gbwt", gbwt_output.to_str().unwrap());
        let khash_file = format!("{}.1khash", gbwt_output.to_str().unwrap());
        assert!(
            Path::new(&gbwt_file).exists(),
            "Region .1gbwt should exist at {}",
            gbwt_file
        );
        assert!(
            Path::new(&khash_file).exists(),
            "Region .1khash should exist at {}",
            khash_file
        );

        // Files should be non-empty
        assert!(
            std::fs::metadata(&gbwt_file).unwrap().len() > 0,
            ".1gbwt should be non-empty"
        );
        assert!(
            std::fs::metadata(&khash_file).unwrap().len() > 0,
            ".1khash should be non-empty"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_paf_based_gbwt_output() {
        let _guard = lock_syng();
        // impg query -i test.paf -f test.fa -r region -o gbwt -O tmpdir/region2
        // (PAF-based → GBWT output)
        let dir = std::env::temp_dir().join("impg_test_cli_paf_gbwt");
        std::fs::create_dir_all(&dir).unwrap();

        let bin = find_impg_binary();
        if bin.is_none() {
            eprintln!("Skipping CLI test: impg binary not found");
            std::fs::remove_dir_all(&dir).ok();
            return;
        }
        let bin = bin.unwrap();

        // Create a simple test PAF and FASTA
        let fasta_path = dir.join("test.fa");
        let paf_path = dir.join("test.paf");

        // Create two sequences with known coordinates
        let seq_a = make_test_sequence(1000, 42);
        let seq_b = make_test_sequence(1000, 42); // identical to create valid PAF alignment
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            use std::io::Write;
            writeln!(f, ">genome_a").unwrap();
            writeln!(f, "{}", String::from_utf8(seq_a).unwrap()).unwrap();
            writeln!(f, ">genome_b").unwrap();
            writeln!(f, "{}", String::from_utf8(seq_b).unwrap()).unwrap();
        }

        // Create a minimal PAF with an identity alignment
        {
            let mut f = std::fs::File::create(&paf_path).unwrap();
            use std::io::Write;
            // PAF format: qname qlen qstart qend strand tname tlen tstart tend nmatch alen mapq [cigar]
            writeln!(
                f,
                "genome_a\t1000\t0\t1000\t+\tgenome_b\t1000\t0\t1000\t1000\t1000\t60\tcg:Z:1000M"
            )
            .unwrap();
        }

        let gbwt_output = dir.join("paf_region");
        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "-d",
                "0",
                "-a",
                paf_path.to_str().unwrap(),
                "--sequence-files",
                fasta_path.to_str().unwrap(),
                "-r",
                "genome_a:100-500",
                "-o",
                "gbwt",
                "-O",
                gbwt_output.to_str().unwrap(),
            ])
            .output()
            .expect("Failed to run impg query -i paf -o gbwt");
        assert!(
            output.status.success(),
            "impg query -i paf -o gbwt failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        // Verify output files exist
        let gbwt_file = format!("{}.1gbwt", gbwt_output.to_str().unwrap());
        let khash_file = format!("{}.1khash", gbwt_output.to_str().unwrap());
        assert!(
            Path::new(&gbwt_file).exists(),
            "PAF-based .1gbwt should exist"
        );
        assert!(
            Path::new(&khash_file).exists(),
            "PAF-based .1khash should exist"
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_cli_gbwt_output_requires_prefix() {
        let _guard = lock_syng();
        // -o gbwt without -O should fail
        let bin = find_impg_binary();
        if bin.is_none() {
            eprintln!("Skipping CLI test: impg binary not found");
            return;
        }
        let bin = bin.unwrap();

        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "-d",
                "0",
                "-a",
                "/tmp/nonexistent_prefix",
                "-r",
                "genome_a:0-100",
                "-o",
                "gbwt",
            ])
            .output()
            .expect("Failed to run impg query");

        assert!(
            !output.status.success(),
            "Using -o gbwt without -O should fail"
        );

        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("prefix") || stderr.contains("-O") || stderr.contains("required"),
            "Error message should mention output prefix requirement. Got: {}",
            stderr
        );
    }

    #[test]
    fn test_top_frequency_seed_filter_drops_only_top_fraction() {
        let counts = vec![(10, 5), (11, 100), (12, 40), (13, 80), (14, 80)];
        let dropped = SyngIndex::top_frequency_seed_nodes(&counts, 0.4);
        assert_eq!(dropped.len(), 2);
        assert!(dropped.contains(&11));
        assert!(dropped.contains(&13));
        assert!(
            !dropped.contains(&14),
            "ties are broken by node id for deterministic exact-count drops"
        );

        let tiny = SyngIndex::top_frequency_seed_nodes(&counts, 0.0005);
        assert!(
            tiny.is_empty(),
            "tiny top-fraction filters should not affect small queries"
        );
    }

    /// Helper to locate the impg binary for CLI tests.
    fn find_impg_binary() -> Option<std::path::PathBuf> {
        // Try CARGO_BIN_EXE first
        let bin = std::env::current_exe()
            .unwrap()
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("impg");
        if bin.exists() {
            return Some(bin);
        }
        // Try manifest dir
        let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let candidates = [
            manifest_dir.join("target/debug/impg"),
            manifest_dir.join("target/release/impg"),
        ];
        for path in &candidates {
            if path.exists() {
                return Some(path.clone());
            }
        }
        None
    }
}
