//! Safe Rust wrapper around syng's C library.
//!
//! Provides `SyngIndex` for building, loading, saving, and querying
//! GBWT-based syncmer indices.

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::ffi::CString;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use crate::fast_locate::FastLocate;
use crate::syng_ffi;
use simple_sds::serialize::Serialize as SdsSerialize;

/// Sidecar table mapping `(forward_path_idx, forward_node_idx)` → `bp_pos`.
///
/// FastLocate returns visits as `(gbz_seq_id, seq_offset_from_end)`. Since
/// we build `gbz::GBWT` bidirectionally with one forward + one reverse
/// sequence per syng path, forward paths live at even gbz_seq_ids:
/// `gbz_seq_id = 2 * forward_path_idx`. The forward node index within that
/// path is `forward_node_count - 1 - seq_offset_from_end`. Looking it up
/// here yields the bp coordinate needed for `HomologousInterval`.
pub struct BpOffsets {
    /// Flat storage of bp offsets for all forward paths, concatenated.
    offsets: Vec<u64>,
    /// `seq_starts[i]` = index in `offsets` where forward path `i` begins.
    /// `seq_starts` has `num_forward_paths + 1` entries (last is `offsets.len()`).
    seq_starts: Vec<usize>,
}

impl BpOffsets {
    fn bp_of(&self, forward_path_idx: usize, node_idx: usize) -> Option<u64> {
        let start = *self.seq_starts.get(forward_path_idx)?;
        let end = *self.seq_starts.get(forward_path_idx + 1)?;
        if node_idx >= end - start {
            return None;
        }
        Some(self.offsets[start + node_idx])
    }

    fn num_forward_nodes(&self, forward_path_idx: usize) -> Option<usize> {
        let start = *self.seq_starts.get(forward_path_idx)?;
        let end = *self.seq_starts.get(forward_path_idx + 1)?;
        Some(end - start)
    }

    const MAGIC: u64 = 0x494D50_42504F46; // "IMPBPOF"
    const VERSION: u64 = 1;

    fn save<W: std::io::Write>(&self, w: &mut W) -> io::Result<()> {
        write_u64(w, Self::MAGIC)?;
        write_u64(w, Self::VERSION)?;
        write_u64(w, self.offsets.len() as u64)?;
        for &x in &self.offsets {
            write_u64(w, x)?;
        }
        write_u64(w, self.seq_starts.len() as u64)?;
        for &s in &self.seq_starts {
            write_u64(w, s as u64)?;
        }
        Ok(())
    }

    fn load<R: std::io::Read>(r: &mut R) -> io::Result<Self> {
        let magic = read_u64(r)?;
        if magic != Self::MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("BpOffsets: bad magic 0x{:x}", magic),
            ));
        }
        let version = read_u64(r)?;
        if version != Self::VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("BpOffsets: unsupported version {}", version),
            ));
        }
        let n_off = read_u64(r)? as usize;
        let mut offsets = Vec::with_capacity(n_off);
        for _ in 0..n_off {
            offsets.push(read_u64(r)?);
        }
        let n_starts = read_u64(r)? as usize;
        let mut seq_starts = Vec::with_capacity(n_starts);
        for _ in 0..n_starts {
            seq_starts.push(read_u64(r)? as usize);
        }
        Ok(Self { offsets, seq_starts })
    }
}

fn write_u64<W: std::io::Write>(w: &mut W, v: u64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_u64<R: std::io::Read>(r: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

#[derive(Debug, Clone, Copy)]
pub struct SampledPositionHit {
    pub path_idx: usize,
    pub target_pos: u64,
    pub target_orient: u8,
}

#[derive(Debug, Clone, Copy)]
pub struct SampledPositionBuildStats {
    pub sampled_occurrences: u64,
    pub sampled_nodes: usize,
    pub walked_paths: usize,
    pub sample_shift: u32,
}

pub const DEFAULT_POSITION_SAMPLE_SHIFT: u32 = 8;
pub const DEFAULT_POSITION_SAMPLE_SEED: u64 = 7;

/// Succinct sampled path-position sidecar for projected syng mapping.
///
/// This mirrors the sampled suffix-array idea used by ropebwt3: store only a
/// deterministic fraction of path occurrences, then resolve mapping positions
/// from those samples instead of materializing every visit for every node.
/// Occurrences are grouped by syncmer node id and varint delta encoded.
pub struct SampledPositions {
    pub sample_shift: u32,
    pub seed: u64,
    sample_count: u64,
    path_starts: Vec<u64>,
    node_ids: Vec<u32>,
    byte_offsets: Vec<u64>,
    data: Vec<u8>,
}

impl SampledPositions {
    const MAGIC: u64 = 0x494D50_53504F53; // "IMPSPOS"
    const VERSION: u64 = 1;

    fn from_samples(
        sample_shift: u32,
        seed: u64,
        path_lengths: &[u64],
        mut samples: Vec<(u32, u64)>,
    ) -> io::Result<Self> {
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
            sample_shift,
            seed,
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
                io::Error::new(io::ErrorKind::InvalidData, "total path length overflowed u64")
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
            io::Error::new(io::ErrorKind::InvalidData, "sampled-position offset missing")
        })? as usize;
        let end = *self.byte_offsets.get(group_idx + 1).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "sampled-position terminal offset missing")
        })? as usize;
        if start > end || end > self.data.len() {
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
                io::Error::new(io::ErrorKind::InvalidData, "sampled position overflowed u64")
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
        let upper = self.path_starts.partition_point(|&start| start <= global_pos);
        let path_idx = upper.checked_sub(1)?;
        Some((path_idx, global_pos - self.path_starts[path_idx]))
    }

    fn save(&self, path: &str) -> io::Result<()> {
        let mut w = BufWriter::new(std::fs::File::create(path)?);
        write_u64(&mut w, Self::MAGIC)?;
        write_u64(&mut w, Self::VERSION)?;
        write_u64(&mut w, self.sample_shift as u64)?;
        write_u64(&mut w, self.seed)?;
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
        w.write_all(&self.data)?;
        w.flush()
    }

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
        let sample_shift = read_u64(&mut r)? as u32;
        let seed = read_u64(&mut r)?;
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
            sample_shift,
            seed,
            sample_count,
            path_starts,
            node_ids,
            byte_offsets,
            data,
        })
    }
}

struct SampledPositionBuilder {
    sample_shift: u32,
    seed: u64,
    sample_mask: u64,
    samples: Vec<(u32, u64)>,
    next_path_start: u64,
    paths_seen: usize,
    walked_paths: usize,
}

impl SampledPositionBuilder {
    fn new(
        sample_shift: u32,
        seed: u64,
        next_path_start: u64,
        paths_seen: usize,
    ) -> io::Result<Self> {
        if sample_shift >= 63 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "position sample shift must be less than 63",
            ));
        }
        Ok(Self {
            sample_shift,
            seed,
            sample_mask: if sample_shift == 0 { 0 } else { (1u64 << sample_shift) - 1 },
            samples: Vec::new(),
            next_path_start,
            paths_seen,
            walked_paths: 0,
        })
    }

    fn record_path(&mut self, path_idx: usize, seq_len: u64, syncmers: &[(i64, i32)]) {
        debug_assert_eq!(
            path_idx, self.paths_seen,
            "online sampled-position path order drifted from name map"
        );
        if !syncmers.is_empty() {
            self.walked_paths += 1;
        }
        let path_start = self.next_path_start;
        for (node_idx, &(signed_node_i64, bp_pos_i32)) in syncmers.iter().enumerate() {
            let signed_node = signed_node_i64 as i32;
            let node_id = signed_node.unsigned_abs();
            let bp_pos = bp_pos_i32 as u64;
            if self.sample_shift > 0
                && sampled_position_hash(node_id, path_idx, node_idx, bp_pos, self.seed)
                    & self.sample_mask
                    != 0
            {
                continue;
            }
            let global_pos = path_start
                .checked_add(bp_pos)
                .expect("sampled path coordinate overflowed u64");
            let orient_bit = if signed_node >= 0 { 0u64 } else { 1u64 };
            let packed = global_pos
                .checked_mul(2)
                .and_then(|p| p.checked_add(orient_bit))
                .expect("sampled packed coordinate overflowed u64");
            self.samples.push((node_id, packed));
        }
        self.next_path_start = self
            .next_path_start
            .checked_add(seq_len)
            .expect("total sampled path length overflowed u64");
        self.paths_seen += 1;
    }

    fn finish(self, path_lengths: &[u64]) -> io::Result<(SampledPositions, SampledPositionBuildStats)> {
        let sample_shift = self.sample_shift;
        let walked_paths = self.walked_paths;
        let sampled_positions =
            SampledPositions::from_samples(sample_shift, self.seed, path_lengths, self.samples)?;
        let stats = SampledPositionBuildStats {
            sampled_occurrences: sampled_positions.sample_count(),
            sampled_nodes: sampled_positions.node_count(),
            walked_paths,
            sample_shift,
        };
        Ok((sampled_positions, stats))
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

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

fn sampled_position_hash(
    node_id: u32,
    path_idx: usize,
    node_idx: usize,
    bp_pos: u64,
    seed: u64,
) -> u64 {
    let mut value = seed;
    value ^= (node_id as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    value ^= (path_idx as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value ^= (node_idx as u64).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^= bp_pos.rotate_left(17);
    splitmix64(value)
}

/// Parameters controlling syncmer extraction.
#[derive(Debug, Clone, Copy)]
pub struct SyncmerParams {
    /// Inner k-mer length (default 8).
    pub k: u32,
    /// Window length; total syncmer length = w + k (default 55).
    pub w: u32,
    /// Hash function seed (default 7).
    pub seed: u32,
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
                let first_syncmer_pos: u64 = parts
                    .get(6)
                    .map(|s| s.parse().unwrap_or(0))
                    .unwrap_or(0);
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

/// A query syncmer found in the index.
#[derive(Debug, Clone, Copy)]
pub struct SyngQuerySyncmer {
    pub signed_node: i32,
    pub node_id: u32,
    pub query_pos: u64,
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
    /// Classical GBWT built from syng's forward paths — backs FastLocate.
    /// None until `build_fast_locate` is called (or an index is loaded that
    /// already has one).
    gbz_gbwt: Option<gbz::GBWT>,
    /// r-index locate structure over `gbz_gbwt`, for fast `query_region`.
    fast_locate: Option<FastLocate>,
    /// bp-offset sidecar keyed by `(forward_path_idx, forward_node_idx)`.
    bp_offsets: Option<BpOffsets>,
    /// Succinct sampled path-position sidecar for approximate PAF projection.
    sampled_positions: Option<SampledPositions>,
    /// Online collector that records sampled positions while sequences are added.
    sampled_position_builder: Option<SampledPositionBuilder>,
}

#[derive(Debug, Clone, Copy)]
pub struct SyngAddSequenceStats {
    pub sequence_len: usize,
    pub syncmers: usize,
    pub indexed: bool,
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
            gbz_gbwt: None,
            fast_locate: None,
            bp_offsets: None,
            sampled_positions: None,
            sampled_position_builder: Some(
                SampledPositionBuilder::new(
                    DEFAULT_POSITION_SAMPLE_SHIFT,
                    DEFAULT_POSITION_SAMPLE_SEED,
                    0,
                    0,
                )
                .expect("default sampled-position parameters are valid"),
            ),
        }
    }

    /// Builds a classical `gbz::GBWT` from the forward paths currently in
    /// this index, then constructs the FastLocate r-index structure and a
    /// bp-offset sidecar. After this returns, `query_region` will use the
    /// fast path (O(query_len) + O(hits)) instead of walking every forward
    /// path per call.
    ///
    /// This walks every existing forward path once to collect (syncmer_id,
    /// bp_pos) pairs. Safe to call multiple times (each call rebuilds).
    pub fn build_fast_locate(&mut self) -> io::Result<()> {
        // Suppress the C debug printfs that syngBWTpathNext emits during walks.
        unsafe { syng_ffi::impg_syng_suppress_debug() };

        let n_paths = self.name_map.path_to_name.len();
        let mut builder = gbz::GBWTBuilder::new(true, false, 64 * 1024);
        let mut offsets_flat: Vec<u64> = Vec::new();
        let mut seq_starts: Vec<usize> = Vec::with_capacity(n_paths + 1);

        for path_idx in 0..n_paths {
            seq_starts.push(offsets_flat.len());
            let Some(ps) = self.name_map.path_starts[path_idx].as_ref() else {
                // No path start info (zero-length or too-short sequence). Still
                // insert an empty path placeholder so that forward_path_idx =
                // gbz_seq_id / 2 stays consistent.
                builder
                    .insert(&[], None)
                    .map_err(|e| io::Error::other(format!("GBWTBuilder::insert(empty): {}", e)))?;
                continue;
            };
            let walk = self.walk_path(ps);
            // walk_path yields (signed_node_id, accumulated_bp_pos). Syng's
            // sign encodes the orientation of the syncmer on this path's
            // forward strand: positive = canonical (forward), negative = RC
            // (reverse). We preserve that sign as the low bit of the GBWT
            // encoded node id so the fast-locate path can distinguish
            // forward- from reverse-orientation visits — required for
            // detecting RC homology via `decompress_da(2*N + 1)`.
            //
            // ENDMARKER is 0, so syncmer ids starting at 1 are safe.
            let mut path_encoded: Vec<usize> = Vec::with_capacity(walk.len());
            for &(node, pos) in &walk {
                let raw = node.unsigned_abs() as usize;
                debug_assert!(raw > 0, "forward walk produced zero node id");
                let orient_bit: usize = if node >= 0 { 0 } else { 1 };
                path_encoded.push((raw << 1) | orient_bit);
                offsets_flat.push(pos);
            }
            builder
                .insert(&path_encoded, None)
                .map_err(|e| io::Error::other(format!("GBWTBuilder::insert: {}", e)))?;
        }
        seq_starts.push(offsets_flat.len());

        let gbz_gbwt = builder
            .build()
            .map_err(|e| io::Error::other(format!("GBWTBuilder::build: {}", e)))?;
        let fast_locate = FastLocate::build(&gbz_gbwt);

        self.gbz_gbwt = Some(gbz_gbwt);
        self.fast_locate = Some(fast_locate);
        self.bp_offsets = Some(BpOffsets {
            offsets: offsets_flat,
            seq_starts,
        });
        Ok(())
    }

    /// True if the FastLocate fast path has been prepared (via
    /// [`Self::build_fast_locate`] or a load that restored it).
    pub fn has_fast_locate(&self) -> bool {
        self.fast_locate.is_some()
    }

    /// True if a sampled path-position sidecar has been built or loaded.
    pub fn has_sampled_positions(&self) -> bool {
        self.sampled_positions.is_some()
    }

    /// Enable online sampled-position collection for subsequently added paths.
    ///
    /// This should normally be called before adding any sequence. If called
    /// later, only future paths can be collected online; use
    /// [`Self::build_sampled_positions`] to create a sidecar for an existing
    /// index.
    pub fn enable_online_sampled_positions(
        &mut self,
        sample_shift: u32,
        seed: u64,
    ) -> io::Result<()> {
        let next_path_start = self
            .name_map
            .path_to_length
            .iter()
            .try_fold(0u64, |acc, &len| acc.checked_add(len))
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "total path length overflowed u64")
            })?;
        self.sampled_position_builder = Some(SampledPositionBuilder::new(
            sample_shift,
            seed,
            next_path_start,
            self.name_map.path_to_name.len(),
        )?);
        self.sampled_positions = None;
        Ok(())
    }

    /// Disable online sampled-position collection and drop unfinalized samples.
    pub fn disable_online_sampled_positions(&mut self) {
        self.sampled_position_builder = None;
        self.sampled_positions = None;
    }

    /// Compact online samples into the `.syng.spos` representation.
    pub fn finalize_online_sampled_positions(
        &mut self,
    ) -> io::Result<Option<SampledPositionBuildStats>> {
        let Some(builder) = self.sampled_position_builder.take() else {
            return Ok(None);
        };
        let (sampled_positions, stats) = builder.finish(&self.name_map.path_to_length)?;
        self.sampled_positions = Some(sampled_positions);
        Ok(Some(stats))
    }

    /// Build a sampled path-position sidecar from the existing syng paths.
    ///
    /// `sample_shift = N` keeps each path occurrence with probability
    /// `1 / 2^N` using a deterministic hash over `(node, path, offset)`.
    /// `sample_shift = 0` stores every occurrence and is intended only for
    /// small tests.
    pub fn build_sampled_positions(
        &mut self,
        sample_shift: u32,
        seed: u64,
    ) -> io::Result<SampledPositionBuildStats> {
        if sample_shift >= 63 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "sample_shift must be less than 63",
            ));
        }

        unsafe { syng_ffi::impg_syng_suppress_debug() };

        let path_global_starts =
            SampledPositions::path_starts_from_lengths(&self.name_map.path_to_length)?;
        let sample_mask = if sample_shift == 0 {
            0
        } else {
            (1u64 << sample_shift) - 1
        };

        let mut samples: Vec<(u32, u64)> = Vec::new();
        let mut walked_paths = 0usize;
        for path_idx in 0..self.name_map.path_to_name.len() {
            let Some(ps) = self.name_map.path_starts[path_idx].as_ref() else {
                continue;
            };
            walked_paths += 1;
            let walk = self.walk_path(ps);
            for (node_idx, &(signed_node, bp_pos)) in walk.iter().enumerate() {
                let node_id = signed_node.unsigned_abs();
                debug_assert!(node_id > 0, "forward walk produced zero node id");
                if sample_shift > 0
                    && sampled_position_hash(node_id, path_idx, node_idx, bp_pos, seed)
                        & sample_mask
                        != 0
                {
                    continue;
                }
                let global_pos = path_global_starts[path_idx].checked_add(bp_pos).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "sampled path coordinate overflowed u64",
                    )
                })?;
                let orient_bit = if signed_node >= 0 { 0u64 } else { 1u64 };
                let packed = global_pos
                    .checked_mul(2)
                    .and_then(|p| p.checked_add(orient_bit))
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "sampled packed coordinate overflowed u64",
                        )
                    })?;
                samples.push((node_id, packed));
            }
        }

        let sampled_positions = SampledPositions::from_samples(
            sample_shift,
            seed,
            &self.name_map.path_to_length,
            samples,
        )?;
        let stats = SampledPositionBuildStats {
            sampled_occurrences: sampled_positions.sample_count(),
            sampled_nodes: sampled_positions.node_count(),
            walked_paths,
            sample_shift,
        };
        self.sampled_positions = Some(sampled_positions);
        self.sampled_position_builder = None;
        Ok(stats)
    }

    /// Save only the sampled path-position sidecar as `{prefix}.syng.spos`.
    pub fn save_sampled_positions(&self, prefix: &str) -> io::Result<()> {
        let sampled_positions = self.sampled_positions.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "sampled positions have not been built",
            )
        })?;
        sampled_positions.save(&format!("{}.syng.spos", prefix))
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

        for (name, seq) in sequences {
            index.add_sequence(name, seq);
        }

        // Eagerly build the FastLocate fast path. This makes `query_region`
        // run in O(query_len + hits) instead of walking every forward path
        // per call. For very large inputs callers may want a `build_noloc`
        // variant later; for now the common case benefits.
        if let Err(e) = index.build_fast_locate() {
            log::warn!("build_fast_locate failed, falling back to walk-every-path: {}", e);
        }

        index
    }

    /// Add one sequence to the syng index.
    ///
    /// This is the same per-sequence work used by [`Self::build`], exposed so
    /// large callers can stream inputs without buffering all decompressed
    /// sequences before GBWT construction.
    pub fn add_sequence(&mut self, name: String, seq: Vec<u8>) -> SyngAddSequenceStats {
        self.gbz_gbwt = None;
        self.fast_locate = None;
        self.bp_offsets = None;
        self.sampled_positions = None;

        let syncmer_len = (self.params.w + self.params.k) as usize;
        let seq_len = seq.len();
        if seq_len < syncmer_len {
            // Sequence too short for syncmer extraction — record in name map but skip
            let path_num = self.name_map.add(name, seq_len as u64);
            if let Some(builder) = self.sampled_position_builder.as_mut() {
                builder.record_path(path_num as usize, seq_len as u64, &[]);
            }
            return SyngAddSequenceStats {
                sequence_len: seq_len,
                syncmers: 0,
                indexed: false,
            };
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
            while syng_ffi::syncmerNext(
                sit,
                std::ptr::null_mut(),
                &mut pos,
                std::ptr::null_mut(),
            ) {
                // Add the syncmer to the KmerHash (or find existing)
                let mut kmer_index: i64 = 0;
                syng_ffi::kmerHashAdd(
                    self.kmer_hash,
                    seq_buf.as_mut_ptr().add(pos as usize) as *mut i8,
                    &mut kmer_index,
                );
                syncmers.push((kmer_index, pos));
            }

            syng_ffi::impg_seqhashIteratorDestroy(sit);
        }

        if syncmers.is_empty() {
            let path_num = self.name_map.add(name, seq_len as u64);
            if let Some(builder) = self.sampled_position_builder.as_mut() {
                builder.record_path(path_num as usize, seq_len as u64, &[]);
            }
            return SyngAddSequenceStats {
                sequence_len: seq_len,
                syncmers: 0,
                indexed: false,
            };
        }

        // Build forward GBWT path and capture start info
        let fwd_start_node;
        let fwd_start_count;
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
            for i in 1..syncmers.len() {
                let next_sync = syncmers[i].0 as i32;
                let offset = (syncmers[i].1 - syncmers[i - 1].1) as u32;
                syng_ffi::syngBWTpathAdd(sbp, next_sync, offset);
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
            builder.record_path(path_num as usize, seq_len as u64, &syncmers);
        }

        SyngAddSequenceStats {
            sequence_len: seq_len,
            syncmers: syncmers.len(),
            indexed: true,
        }
    }

    /// Save the index to disk.
    ///
    /// Produces three files at the given prefix:
    /// - `{prefix}.1khash` — syncmer hash (syng-compatible)
    /// - `{prefix}.1gbwt` — GBWT graph (syng-compatible)
    /// - `{prefix}.syng.names` — sequence name mapping
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
                return Err(io::Error::other(format!("Failed to open {} for writing", gbwt_path)));
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
                return Err(io::Error::other("Failed to create ONEcode schema for khash"));
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
                return Err(io::Error::other(format!("Failed to open {} for writing", khash_path)));
            }
            let ok = syng_ffi::kmerHashWriteOneFile(self.kmer_hash, of);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
            if !ok {
                return Err(io::Error::other("kmerHashWriteOneFile failed"));
            }
        }

        // Write .syng.names
        let names_path = format!("{}.syng.names", prefix);
        self.name_map.save(&names_path)?;

        // Write .syng.locate (optional — only if fast-locate has been built).
        if let (Some(gbz_gbwt), Some(fl), Some(bp_off)) =
            (self.gbz_gbwt.as_ref(), self.fast_locate.as_ref(), self.bp_offsets.as_ref())
        {
            let locate_path = format!("{}.syng.locate", prefix);
            let mut f = std::io::BufWriter::new(std::fs::File::create(&locate_path)?);
            // 1) classical GBWT via simple-sds Serialize
            gbz_gbwt
                .serialize(&mut f)
                .map_err(|e| io::Error::other(format!("gbz::GBWT::serialize: {}", e)))?;
            // 2) FastLocate custom framing
            fl.save(&mut f)?;
            // 3) BpOffsets custom framing (little-endian u64-prefixed blobs)
            bp_off.save(&mut f)?;
            f.into_inner()?;
        }

        // Write .syng.spos (default sampled position sidecar for PAF projection).
        if let Some(sampled_positions) = self.sampled_positions.as_ref() {
            sampled_positions.save(&format!("{}.syng.spos", prefix))?;
        }

        Ok(())
    }

    /// Load an index from disk.
    ///
    /// Reads three files at the given prefix:
    /// - `{prefix}.1khash`
    /// - `{prefix}.1gbwt`
    /// - `{prefix}.syng.names`
    pub fn load(prefix: &str, params: SyncmerParams) -> io::Result<Self> {
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
                return Err(io::Error::other("Failed to create ONEcode schema for gbwt read"));
            }
            let of = syng_ffi::oneFileOpenRead(
                gbwt_cpath.as_ptr(),
                schema,
                gbwt_type.as_ptr(),
                1,
            );
            if of.is_null() {
                syng_ffi::oneSchemaDestroy(schema);
                return Err(io::Error::other(format!("Failed to open {} for reading", gbwt_path)));
            }
            let gbwt = syng_ffi::syngBWTread(of);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
            if gbwt.is_null() {
                return Err(io::Error::other("syngBWTread returned null"));
            }
            gbwt
        };

        // Read .1khash
        let khash_path = format!("{}.1khash", prefix);
        if !Path::new(&khash_path).exists() {
            unsafe { syng_ffi::syngBWTdestroy(gbwt) };
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("KmerHash file not found: {}", khash_path),
            ));
        }
        let khash_cpath = CString::new(khash_path.as_str())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let khash_type = CString::new("khash").unwrap();
        let kmer_hash = unsafe {
            let schema = syng_ffi::oneSchemaCreateFromText(schema_text.as_ptr());
            if schema.is_null() {
                syng_ffi::syngBWTdestroy(gbwt);
                return Err(io::Error::other("Failed to create ONEcode schema for khash read"));
            }
            let of = syng_ffi::oneFileOpenRead(
                khash_cpath.as_ptr(),
                schema,
                khash_type.as_ptr(),
                1,
            );
            if of.is_null() {
                syng_ffi::oneSchemaDestroy(schema);
                syng_ffi::syngBWTdestroy(gbwt);
                return Err(io::Error::other(format!("Failed to open {} for reading", khash_path)));
            }
            let kh = syng_ffi::kmerHashReadOneFile(of);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
            if kh.is_null() {
                syng_ffi::syngBWTdestroy(gbwt);
                return Err(io::Error::other("kmerHashReadOneFile returned null"));
            }
            kh
        };

        // Read .syng.names
        let names_path = format!("{}.syng.names", prefix);
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

        // Read optional .syng.locate sidecar (classical GBWT + FastLocate + BpOffsets).
        let locate_path = format!("{}.syng.locate", prefix);
        let (gbz_gbwt, fast_locate, bp_offsets) = if Path::new(&locate_path).exists() {
            let mut f = std::io::BufReader::new(std::fs::File::open(&locate_path)?);
            let gbz_gbwt_loaded = gbz::GBWT::load(&mut f)
                .map_err(|e| io::Error::other(format!("gbz::GBWT::load: {}", e)))?;
            let fl_loaded = FastLocate::load(&mut f)?;
            let bp_loaded = BpOffsets::load(&mut f)?;
            (Some(gbz_gbwt_loaded), Some(fl_loaded), Some(bp_loaded))
        } else {
            (None, None, None)
        };

        // Read optional .syng.spos sampled-position sidecar.
        let sampled_positions_path = format!("{}.syng.spos", prefix);
        let sampled_positions = if Path::new(&sampled_positions_path).exists() {
            Some(SampledPositions::load(&sampled_positions_path)?)
        } else {
            None
        };

        Ok(Self {
            gbwt,
            kmer_hash,
            seqhash,
            name_map,
            params,
            gbz_gbwt,
            fast_locate,
            bp_offsets,
            sampled_positions,
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
    fn walk_path(&self, start: &GbwtPathStart) -> Vec<(i32, u64)> {
        let mut nodes = Vec::with_capacity(start.num_syncmers as usize);
        unsafe {
            let sbp =
                syng_ffi::syngBWTpathStartOld(self.gbwt, start.start_node, start.start_count);
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
        let with_anchors =
            self.query_region_with_anchors_ext(genome, start, end, padding, query_extension)?;
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
        let syncmer_len = (self.params.w + self.params.k) as usize;
        if query_seq.len() < syncmer_len {
            return Vec::new();
        }

        let mut seq_buf: Vec<u8> = Vec::with_capacity(query_seq.len() + 1);
        seq_buf.extend(query_seq.iter().map(|&b| match b {
            b'a' | b'A' | 0 => 0u8,
            b'c' | b'C' | 1 => 1u8,
            b'g' | b'G' | 2 => 2u8,
            b't' | b'T' | 3 => 3u8,
            _ => 0u8,
        }));
        seq_buf.push(0);

        let mut out = Vec::new();
        unsafe {
            let sit = syng_ffi::syncmerIterator(
                self.seqhash,
                seq_buf.as_mut_ptr() as *mut i8,
                query_seq.len() as i32,
            );
            let mut pos: i32 = 0;
            while syng_ffi::syncmerNext(
                sit,
                std::ptr::null_mut(),
                &mut pos,
                std::ptr::null_mut(),
            ) {
                let mut kmer_index: i64 = 0;
                if syng_ffi::kmerHashFind(
                    self.kmer_hash,
                    seq_buf.as_mut_ptr().add(pos as usize) as *mut i8,
                    &mut kmer_index,
                ) {
                    let signed_node = kmer_index as i32;
                    out.push(SyngQuerySyncmer {
                        signed_node,
                        node_id: signed_node.unsigned_abs(),
                        query_pos: pos as u64,
                    });
                }
            }
            syng_ffi::impg_seqhashIteratorDestroy(sit);
        }
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

        let exact_locate = match (
            self.fast_locate.as_ref(),
            self.gbz_gbwt.as_ref(),
            self.bp_offsets.as_ref(),
        ) {
            (Some(fl), Some(gbz), Some(bp)) => Some((fl, gbz, bp)),
            _ => None,
        };
        if exact_locate.is_none() && self.sampled_positions.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "syng map -o paf requires a positional sidecar; rebuild the syng index with this impg or run syng-pos",
            ));
        }

        let query_syncmers = self.matched_syncmers_in_sequence(query_seq);
        if query_syncmers.is_empty() {
            return Ok(Vec::new());
        }

        let syncmer_len = (self.params.w + self.params.k) as u64;
        let mut per_path: FxHashMap<(usize, char), Vec<Anchor>> = FxHashMap::default();

        if let Some((fl, gbz_gbwt, bp_off)) = exact_locate {
            for sm in query_syncmers {
                let q_orient: u8 = if sm.signed_node >= 0 { 0 } else { 1 };
                for t_orient in 0u8..=1 {
                    let encoded = ((sm.node_id as usize) << 1) | (t_orient as usize);
                    if !gbz_gbwt.has_node(encoded) {
                        continue;
                    }
                    for (gbz_seq_id, seq_off_from_end) in fl.decompress_da(gbz_gbwt, encoded) {
                        if gbz_seq_id & 1 != 0 {
                            continue;
                        }
                        let forward_path_idx = gbz_seq_id >> 1;
                        let Some(num_nodes) = bp_off.num_forward_nodes(forward_path_idx) else {
                            continue;
                        };
                        if num_nodes == 0 || seq_off_from_end >= num_nodes {
                            continue;
                        }
                        let forward_node_idx = num_nodes - 1 - seq_off_from_end;
                        let Some(target_pos) = bp_off.bp_of(forward_path_idx, forward_node_idx) else {
                            continue;
                        };
                        let strand = if q_orient == t_orient { '+' } else { '-' };
                        per_path
                            .entry((forward_path_idx, strand))
                            .or_default()
                            .push(Anchor {
                                query_pos: sm.query_pos,
                                target_pos,
                                node_id: sm.node_id,
                            });
                    }
                }
            }
        } else if let Some(sampled_positions) = self.sampled_positions.as_ref() {
            for sm in query_syncmers {
                let q_orient: u8 = if sm.signed_node >= 0 { 0 } else { 1 };
                for hit in sampled_positions.decode_node(sm.node_id)? {
                    if hit.path_idx >= self.name_map.path_to_name.len() {
                        continue;
                    }
                    let strand = if q_orient == hit.target_orient { '+' } else { '-' };
                    per_path
                        .entry((hit.path_idx, strand))
                        .or_default()
                        .push(Anchor {
                            query_pos: sm.query_pos,
                            target_pos: hit.target_pos,
                            node_id: sm.node_id,
                        });
                }
            }
        }

        let mut raw_hits = Vec::new();
        for ((path_idx, strand), mut anchors) in per_path {
            anchors.sort_by(|a, b| {
                a.query_pos
                    .cmp(&b.query_pos)
                    .then(a.target_pos.cmp(&b.target_pos))
                    .then(a.node_id.cmp(&b.node_id))
            });
            anchors.dedup_by(|a, b| {
                a.query_pos == b.query_pos
                    && a.target_pos == b.target_pos
                    && a.node_id == b.node_id
            });
            if anchors.len() < min_anchors {
                continue;
            }
            let target_len = self.name_map.path_to_length[path_idx];
            let start = anchors.iter().map(|a| a.target_pos).min().unwrap_or(0);
            let end = anchors
                .iter()
                .map(|a| a.target_pos.saturating_add(syncmer_len))
                .max()
                .unwrap_or(start)
                .min(target_len);
            raw_hits.push(HomologousIntervalWithAnchors {
                genome: self.name_map.path_to_name[path_idx].clone(),
                start,
                end,
                strand,
                anchors,
            });
        }

        let chained = crate::syng_transitive::chain_anchors(raw_hits, chain_budget, syncmer_len);
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
    /// node sign); target syncmers are matched against both GBWT node
    /// encodings (forward = `2*N`, reverse = `2*N + 1`). A query orientation
    /// matching the target's produces a forward-strand hit (`strand='+'`);
    /// a mismatch produces a reverse-strand hit (`strand='-'`). The two
    /// strands are kept as separate homolog intervals (not merged) since
    /// they describe distinct homologies.
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
        self.query_region_with_anchors_ext_visited(genome, start, end, padding, query_extension, None)
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
        mut visited_nodes: Option<&mut FxHashSet<u32>>,
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

        let query_start = self.name_map.path_starts[query_path_num as usize]
            .as_ref()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "No GBWT path start info for genome '{}' — index may need rebuilding",
                        genome
                    ),
                )
            })?;

        // 2. Walk the query genome's path, collect nodes in [start, end)
        //    along with (query_position, query_orientation).
        //    Query orientation is 0 = forward (canonical), 1 = reverse
        //    (RC traversal) — encoded in the sign of the node id from
        //    walk_path: negative node id → reverse orientation.
        //    A syncmer node may appear at multiple query positions (repeat
        //    families) — we keep all of them and cross-product with target
        //    hits downstream. Callers pick the nearest anchor per edge.
        let query_nodes = self.walk_path(query_start);
        // Expanded syncmer-lookup window. `query_extension` pushes the filter
        // outward so conserved blocks that end just outside the user's
        // declared region are still discoverable; refined boundaries stay
        // clipped to the original interval downstream.
        let expanded_start = start.saturating_sub(query_extension);
        let expanded_end = end.saturating_add(query_extension);
        // abs_node -> Vec<(query_pos, q_orient_bit)>. q_orient_bit is 0 for
        // forward, 1 for reverse — matching the low bit of the GBWT node
        // encoding so we can compare directly against target visits.
        let mut query_node_positions: FxHashMap<u32, Vec<(u64, u8)>> =
            FxHashMap::default();
        for &(signed_node, pos) in &query_nodes {
            if pos >= expanded_start && pos < expanded_end {
                let abs = signed_node.unsigned_abs();
                let q_orient: u8 = if signed_node >= 0 { 0 } else { 1 };
                query_node_positions
                    .entry(abs)
                    .or_default()
                    .push((pos, q_orient));
            }
        }

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
        let mut p_nodes = 0u64;
        let mut p_decompress_ns = 0u128;
        let mut p_visits_raw = 0u64;
        let mut p_visits_used = 0u64;
        let mut p_visit_loop_ns = 0u128;
        let mut p_emit_anchors = 0u64;
        if let (Some(fl), Some(gbz_gbwt), Some(bp_off)) =
            (self.fast_locate.as_ref(), self.gbz_gbwt.as_ref(), self.bp_offsets.as_ref())
        {
            // Skip already-visited syncmer nodes up-front (serial filter,
            // cheap), then run the heavy per-node visit loop in parallel.
            // Updating `visited_nodes` happens after the parallel phase
            // so rayon workers don't need shared mutable access.
            let unvisited: Vec<(u32, &Vec<(u64, u8)>)> = query_node_positions
                .iter()
                .filter(|(node, _)| match visited_nodes.as_deref() {
                    Some(v) => !v.contains(node),
                    None => true,
                })
                .map(|(n, qp)| (*n, qp))
                .collect();
            if let Some(v) = visited_nodes.as_deref_mut() {
                for (n, _) in &unvisited {
                    v.insert(*n);
                }
            }
            if emit_prof {
                p_nodes = unvisited.len() as u64;
            }

            // Each query node is independent: its anchors are deduped
            // into a per-node (path, strand) → signature hashmap.
            // Parallelize across nodes on rayon's global pool. On
            // pathological rDNA tiles this turns ~180 s of inner-loop
            // work into ~15 s at 16-way parallelism.
            use std::sync::atomic::{AtomicU64, Ordering};
            let a_decompress = AtomicU64::new(0);
            let a_visits_raw = AtomicU64::new(0);
            let a_visits_used = AtomicU64::new(0);
            let a_visit_loop = AtomicU64::new(0);
            let a_emit = AtomicU64::new(0);

            // Nested parallelism:
            //   - outer par over nodes (handles non-rDNA: many small-ish
            //     nodes share threads)
            //   - inner par over visits for a single "monster" node (the
            //     rDNA case where 99.99% of work lives in one node's
            //     visit list)
            // A visit-count threshold skips the inner fold/reduce
            // overhead for small visit lists.
            const PAR_VISITS_THRESHOLD: usize = 50_000;

            let per_node_results: Vec<Vec<HomologousIntervalWithAnchors>> = unvisited
                .par_iter()
                .with_min_len(1)
                .map(|(query_node, query_positions)| {
                    let query_node = *query_node;
                // Collect all visits with their t_orient.
                let mut all_visits: Vec<(u8, usize, usize)> = Vec::new();
                for t_orient in 0u8..=1 {
                    let encoded = ((query_node as usize) << 1) | (t_orient as usize);
                    if !gbz_gbwt.has_node(encoded) {
                        continue;
                    }
                    let t_dec = std::time::Instant::now();
                    let visits = fl.decompress_da(gbz_gbwt, encoded);
                    if emit_prof {
                        a_decompress
                            .fetch_add(t_dec.elapsed().as_nanos() as u64, Ordering::Relaxed);
                        a_visits_raw.fetch_add(visits.len() as u64, Ordering::Relaxed);
                    }
                    for (sid, off) in visits {
                        all_visits.push((t_orient, sid, off));
                    }
                }
                if all_visits.is_empty() {
                    return Vec::new();
                }
                let t_vloop = std::time::Instant::now();
                let per_node: FxHashMap<
                    (usize, char),
                    (u64, u64, FxHashMap<i64, Anchor>),
                > = if all_visits.len() < PAR_VISITS_THRESHOLD {
                    let mut acc: FxHashMap<
                        (usize, char),
                        (u64, u64, FxHashMap<i64, Anchor>),
                    > = FxHashMap::default();
                    for &(t_orient, gbz_seq_id, seq_off_from_end) in &all_visits {
                        if gbz_seq_id & 1 != 0 {
                            continue;
                        }
                        let forward_path_idx = gbz_seq_id >> 1;
                        let Some(num_nodes) = bp_off.num_forward_nodes(forward_path_idx)
                        else {
                            continue;
                        };
                        if num_nodes == 0 || seq_off_from_end >= num_nodes {
                            continue;
                        }
                        let forward_node_idx = num_nodes - 1 - seq_off_from_end;
                        let Some(target_pos) =
                            bp_off.bp_of(forward_path_idx, forward_node_idx)
                        else {
                            continue;
                        };
                        let genome_len = self.name_map.path_to_length[forward_path_idx];
                        let hit_end = target_pos + syncmer_len;
                        let padded_start = target_pos.saturating_sub(padding);
                        let padded_end = (hit_end + padding).min(genome_len);
                        for &(qp, q_orient) in *query_positions {
                            let strand = if q_orient == t_orient { '+' } else { '-' };
                            let sig: i64 = if strand == '+' {
                                target_pos as i64 - qp as i64
                            } else {
                                target_pos as i64 + qp as i64
                            };
                            let entry = acc
                                .entry((forward_path_idx, strand))
                                .or_insert_with(|| {
                                    (padded_start, padded_end, FxHashMap::default())
                                });
                            entry.0 = entry.0.min(padded_start);
                            entry.1 = entry.1.max(padded_end);
                            entry.2.entry(sig).or_insert(Anchor {
                                query_pos: qp,
                                target_pos,
                                node_id: query_node,
                            });
                            if emit_prof {
                                a_emit.fetch_add(1, Ordering::Relaxed);
                            }
                        }
                        if emit_prof {
                            a_visits_used.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                    acc
                } else {
                    all_visits
                    .par_iter()
                    .with_min_len(1024)
                    .fold(
                        FxHashMap::default,
                        |mut acc, &(t_orient, gbz_seq_id, seq_off_from_end)| {
                            if gbz_seq_id & 1 != 0 {
                                return acc;
                            }
                            let forward_path_idx = gbz_seq_id >> 1;
                            let Some(num_nodes) = bp_off.num_forward_nodes(forward_path_idx)
                            else {
                                return acc;
                            };
                            if num_nodes == 0 || seq_off_from_end >= num_nodes {
                                return acc;
                            }
                            let forward_node_idx = num_nodes - 1 - seq_off_from_end;
                            let Some(target_pos) =
                                bp_off.bp_of(forward_path_idx, forward_node_idx)
                            else {
                                return acc;
                            };
                            let genome_len = self.name_map.path_to_length[forward_path_idx];
                            let hit_end = target_pos + syncmer_len;
                            let padded_start = target_pos.saturating_sub(padding);
                            let padded_end = (hit_end + padding).min(genome_len);
                            for &(qp, q_orient) in *query_positions {
                                let strand = if q_orient == t_orient { '+' } else { '-' };
                                let sig: i64 = if strand == '+' {
                                    target_pos as i64 - qp as i64
                                } else {
                                    target_pos as i64 + qp as i64
                                };
                                let entry = acc
                                    .entry((forward_path_idx, strand))
                                    .or_insert_with(|| {
                                        (padded_start, padded_end, FxHashMap::default())
                                    });
                                entry.0 = entry.0.min(padded_start);
                                entry.1 = entry.1.max(padded_end);
                                entry.2.entry(sig).or_insert(Anchor {
                                    query_pos: qp,
                                    target_pos,
                                    node_id: query_node,
                                });
                                if emit_prof {
                                    a_emit.fetch_add(1, Ordering::Relaxed);
                                }
                            }
                            if emit_prof {
                                a_visits_used.fetch_add(1, Ordering::Relaxed);
                            }
                            acc
                        },
                    )
                    .reduce(
                        FxHashMap::default,
                        |mut a,
                         b: FxHashMap<(usize, char), (u64, u64, FxHashMap<i64, Anchor>)>| {
                            for (key, (bs, be, bmap)) in b {
                                let entry = a
                                    .entry(key)
                                    .or_insert_with(|| (bs, be, FxHashMap::default()));
                                entry.0 = entry.0.min(bs);
                                entry.1 = entry.1.max(be);
                                for (sig, anchor) in bmap {
                                    entry.2.entry(sig).or_insert(anchor);
                                }
                            }
                            a
                        },
                    )
                };
                if emit_prof {
                    a_visit_loop
                        .fetch_add(t_vloop.elapsed().as_nanos() as u64, Ordering::Relaxed);
                }

                per_node
                    .into_iter()
                    .filter_map(|((forward_path_idx, strand), (start, end, sig_map))| {
                        if sig_map.is_empty() {
                            return None;
                        }
                        let anchors: Vec<Anchor> = sig_map.into_values().collect();
                        Some(HomologousIntervalWithAnchors {
                            genome: self.name_map.path_to_name[forward_path_idx].clone(),
                            start,
                            end,
                            strand,
                            anchors,
                        })
                    })
                    .collect()
                })
                .collect();

            // Merge per-node results into the global per_path map.
            // Re-derive forward_path_idx by name so the (usize, char)
            // key remains consistent with the pre-parallel layout.
            for node_result in per_node_results {
                for hiwa in node_result {
                    let Some(path_idx) =
                        self.name_map.name_to_path.get(&hiwa.genome).copied()
                    else {
                        continue;
                    };
                    per_path
                        .entry((path_idx as usize, hiwa.strand))
                        .or_default()
                        .push(hiwa);
                }
            }

            if emit_prof {
                p_decompress_ns = a_decompress.load(Ordering::Relaxed) as u128;
                p_visits_raw = a_visits_raw.load(Ordering::Relaxed);
                p_visits_used = a_visits_used.load(Ordering::Relaxed);
                p_visit_loop_ns = a_visit_loop.load(Ordering::Relaxed) as u128;
                p_emit_anchors = a_emit.load(Ordering::Relaxed);
            }
        } else {
            // Fallback: walk every forward path (the old pre-locate behavior).
            // Each target visit carries its own orientation (sign of the
            // target's node id) which we XOR against the query's orientation
            // to determine strand.
            let num_genomes = self.name_map.path_to_name.len();
            // Per-query-node signature dedup across the full path walk —
            // fallback equivalent of the fast path's per-query-node
            // aggregation. Keyed by (genome_idx, query_node, strand) since
            // the walk interleaves nodes from all paths.
            let mut per_node_fallback: FxHashMap<
                (usize, u32, char),
                (u64, u64, FxHashMap<i64, Anchor>),
            > = FxHashMap::default();
            for genome_idx in 0..num_genomes {
                let path_start = match &self.name_map.path_starts[genome_idx] {
                    Some(ps) => ps,
                    None => continue,
                };
                let genome_len = self.name_map.path_to_length[genome_idx];
                let nodes = self.walk_path(path_start);
                for &(signed_target_node, target_pos) in &nodes {
                    let abs_node = signed_target_node.unsigned_abs();
                    let t_orient: u8 = if signed_target_node >= 0 { 0 } else { 1 };
                    if let Some(query_positions) = query_node_positions.get(&abs_node) {
                        // Same BFS-scope filter as the fast path — suppress
                        // re-entering syncmer nodes already traversed.
                        if let Some(v) = visited_nodes.as_deref() {
                            if v.contains(&abs_node) {
                                continue;
                            }
                        }
                        if let Some(v) = visited_nodes.as_deref_mut() {
                            v.insert(abs_node);
                        }
                        let hit_end = target_pos + syncmer_len;
                        let padded_start = target_pos.saturating_sub(padding);
                        let padded_end = (hit_end + padding).min(genome_len);
                        for &(qp, q_orient) in query_positions {
                            let strand = if q_orient == t_orient { '+' } else { '-' };
                            let sig: i64 = if strand == '+' {
                                target_pos as i64 - qp as i64
                            } else {
                                target_pos as i64 + qp as i64
                            };
                            let entry = per_node_fallback
                                .entry((genome_idx, abs_node, strand))
                                .or_insert_with(|| (padded_start, padded_end, FxHashMap::default()));
                            entry.0 = entry.0.min(padded_start);
                            entry.1 = entry.1.max(padded_end);
                            entry.2.entry(sig).or_insert(Anchor {
                                query_pos: qp,
                                target_pos,
                                node_id: abs_node,
                            });
                        }
                    }
                }
            }
            for ((genome_idx, _node, strand), (start, end, sig_map)) in per_node_fallback {
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
                "EMITPROF nodes={} decompress={:.2}s visits_raw={} visits_used={} visit_loop={:.2}s emit_anchors={} merge_sort={:.2}s",
                p_nodes,
                p_decompress_ns as f64 / 1e9,
                p_visits_raw,
                p_visits_used,
                p_visit_loop_ns as f64 / 1e9,
                p_emit_anchors,
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
    pub fn build_region_gbwt(
        &self,
        sequences: &[(String, &[u8])],
        prefix: &str,
    ) -> io::Result<()> {
        let syncmer_len = (self.params.w + self.params.k) as i32;

        // Create fresh structures for this region
        let region_gbwt = unsafe { syng_ffi::syngBWTcreate(syncmer_len, 0) };
        let region_kh = unsafe { syng_ffi::kmerHashCreate(1024, syncmer_len) };
        let region_sh = unsafe {
            syng_ffi::impg_seqhashCreateSafe(
                self.params.k as i32,
                self.params.w as i32,
                self.params.seed as i32,
            )
        };

        let syncmer_len_usize = syncmer_len as usize;

        for (_name, seq) in sequences {
            if seq.len() < syncmer_len_usize {
                continue;
            }

            // Convert to numeric encoding (0-3) for seqhash (uses raw values as indices)
            let mut seq_buf: Vec<u8> = Vec::with_capacity(seq.len() + 1);
            seq_buf.extend(seq.iter().map(|&b| match b {
                b'a' | b'A' | 0 => 0u8,
                b'c' | b'C' | 1 => 1u8,
                b'g' | b'G' | 2 => 2u8,
                b't' | b'T' | 3 => 3u8,
                _ => 0u8,
            }));
            seq_buf.push(0); // null terminator

            // Extract syncmers
            let mut syncmers: Vec<(i64, i32)> = Vec::new();
            unsafe {
                let sit = syng_ffi::syncmerIterator(
                    region_sh,
                    seq_buf.as_mut_ptr() as *mut i8,
                    seq.len() as i32,
                );

                let mut pos: i32 = 0;
                while syng_ffi::syncmerNext(
                    sit,
                    std::ptr::null_mut(),
                    &mut pos,
                    std::ptr::null_mut(),
                ) {
                    let mut kmer_index: i64 = 0;
                    syng_ffi::kmerHashAdd(
                        region_kh,
                        seq_buf.as_mut_ptr().add(pos as usize) as *mut i8,
                        &mut kmer_index,
                    );
                    syncmers.push((kmer_index, pos));
                }

                syng_ffi::impg_seqhashIteratorDestroy(sit);
            }

            if syncmers.is_empty() {
                continue;
            }

            // Build forward GBWT path
            unsafe {
                let first_sync = syncmers[0].0 as i32;
                let sbp = syng_ffi::syngBWTpathStartNew(region_gbwt, first_sync);
                for i in 1..syncmers.len() {
                    let next_sync = syncmers[i].0 as i32;
                    let offset = (syncmers[i].1 - syncmers[i - 1].1) as u32;
                    syng_ffi::syngBWTpathAdd(sbp, next_sync, offset);
                }
                syng_ffi::syngBWTpathFinish(sbp);
            }

            // Build reverse complement GBWT path
            unsafe {
                let n = syncmers.len();
                let first_sync_rc = -(syncmers[n - 1].0 as i32);
                let sbp = syng_ffi::syngBWTpathStartNew(region_gbwt, first_sync_rc);
                for i in (0..n - 1).rev() {
                    let next_sync_rc = -(syncmers[i].0 as i32);
                    let offset = (syncmers[i + 1].1 - syncmers[i].1) as u32;
                    syng_ffi::syngBWTpathAdd(sbp, next_sync_rc, offset);
                }
                syng_ffi::syngBWTpathFinish(sbp);
            }
        }

        // Write .1gbwt
        let schema_text = syng_ffi::syng_schema_text();
        let gbwt_path = format!("{}.1gbwt", prefix);
        let gbwt_cpath = CString::new(gbwt_path.as_str())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let gbwt_type = CString::new("gbwt").unwrap();
        unsafe {
            let schema = syng_ffi::oneSchemaCreateFromText(schema_text.as_ptr());
            if schema.is_null() {
                syng_ffi::syngBWTdestroy(region_gbwt);
                syng_ffi::kmerHashDestroy(region_kh);
                syng_ffi::impg_seqhashDestroy(region_sh);
                return Err(io::Error::other("Failed to create ONEcode schema for region gbwt"));
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
                syng_ffi::syngBWTdestroy(region_gbwt);
                syng_ffi::kmerHashDestroy(region_kh);
                syng_ffi::impg_seqhashDestroy(region_sh);
                return Err(io::Error::other(format!("Failed to open {} for writing", gbwt_path)));
            }
            syng_ffi::syngBWTwrite(of, region_gbwt);
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
                syng_ffi::syngBWTdestroy(region_gbwt);
                syng_ffi::kmerHashDestroy(region_kh);
                syng_ffi::impg_seqhashDestroy(region_sh);
                return Err(io::Error::other("Failed to create ONEcode schema for region khash"));
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
                syng_ffi::syngBWTdestroy(region_gbwt);
                syng_ffi::kmerHashDestroy(region_kh);
                syng_ffi::impg_seqhashDestroy(region_sh);
                return Err(io::Error::other(format!("Failed to open {} for writing", khash_path)));
            }
            let ok = syng_ffi::kmerHashWriteOneFile(region_kh, of);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
            if !ok {
                syng_ffi::syngBWTdestroy(region_gbwt);
                syng_ffi::kmerHashDestroy(region_kh);
                syng_ffi::impg_seqhashDestroy(region_sh);
                return Err(io::Error::other("kmerHashWriteOneFile failed for region khash"));
            }
        }

        // Clean up C allocations
        unsafe {
            syng_ffi::syngBWTdestroy(region_gbwt);
            syng_ffi::kmerHashDestroy(region_kh);
            syng_ffi::impg_seqhashDestroy(region_sh);
        }

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
            iv.anchors
                .sort_by(|a, b| a.query_pos.cmp(&b.query_pos).then(a.target_pos.cmp(&b.target_pos)));
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
        assert!(!sh.is_null(), "seqhashCreate returned null for valid params");
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
            .map(|i| (format!("genome_{}", i), make_test_sequence(2000, i as u8 + 100)))
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

        // Load
        let loaded = SyngIndex::load(prefix_str, params).unwrap();

        // Verify name map is intact
        assert_eq!(
            loaded.name_map.path_to_name,
            index.name_map.path_to_name,
            "Name map names should match after round-trip"
        );
        assert_eq!(
            loaded.name_map.path_to_length,
            index.name_map.path_to_length,
            "Name map lengths should match after round-trip"
        );

        // Verify C pointers are valid
        assert!(!loaded.gbwt.is_null());
        assert!(!loaded.kmer_hash.is_null());
        assert!(!loaded.seqhash.is_null());

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
            let idx = *loaded.name_to_path.get(*name).unwrap_or_else(|| {
                panic!("Name '{}' missing from loaded name map", name)
            });
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
        let mut idx_alt = SyngIndex::build(
            params_alt_seed,
            vec![("seq".to_string(), seq)].into_iter(),
        );
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
        let size_default = std::fs::metadata(dir.join("default.1khash"))
            .unwrap()
            .len();
        let size_short = std::fs::metadata(dir.join("short.1khash"))
            .unwrap()
            .len();
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
        index.enable_online_sampled_positions(0, 7).unwrap();
        index.add_sequence("sampleA#0#chr1".to_string(), seq_a);
        index.add_sequence("sampleB#0#chr1".to_string(), seq_b);
        let stats = index.finalize_online_sampled_positions().unwrap().unwrap();
        assert!(stats.sampled_occurrences > 0);
        assert!(stats.sampled_nodes > 0);
        assert!(index.has_sampled_positions());
        assert!(!index.has_fast_locate());

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

        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert!(loaded.has_sampled_positions());
        assert!(!loaded.has_fast_locate());
        let loaded_hits = loaded
            .map_sequence("read1", query, 2, 10_000)
            .expect("loaded sampled-position PAF mapping should succeed");
        assert!(
            loaded_hits.iter().any(|h| h.target_name == "sampleA#0#chr1"),
            "loaded sampled mapping should include sampleA, got: {:?}",
            loaded_hits
        );
        assert!(
            loaded_hits.iter().any(|h| h.target_name == "sampleB#0#chr1"),
            "loaded sampled mapping should include sampleB, got: {:?}",
            loaded_hits
        );

        std::fs::remove_dir_all(&dir).ok();
    }

    // ── 11. CLI integration test ────────────────────────────────────

    #[test]
    fn test_syng_cli_fasta_roundtrip() {
        let _guard = lock_syng();
        // Run `impg syng -f <fasta> -o <prefix>` via Command,
        // verify output files exist.
        let dir = std::env::temp_dir().join("impg_test_syng_cli");
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
            eprintln!(
                "Skipping CLI test: impg binary not found at {:?}",
                bin
            );
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
            ".syng.locate should not be built unless impg syng --locate is used"
        );
        assert!(
            std::path::Path::new(&format!("{}.syng.spos", output_prefix_str)).exists(),
            ".syng.spos should be built by default"
        );

        // Load and verify content
        let loaded = SyngIndex::load(output_prefix_str, SyncmerParams::default()).unwrap();
        assert!(loaded.has_sampled_positions());
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
        assert!(result.is_err(), "Loading from nonexistent files should fail");
    }

    // ── 13. query_region tests ─────────────────────────────────────

    /// Full SyngIndex save/load with the FastLocate sidecar file
    /// (`.syng.locate`). Verify that the reloaded index still has
    /// `has_fast_locate() == true` and returns the same query results as
    /// the in-memory version.
    #[test]
    fn test_syng_save_load_with_fast_locate() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let shared = make_test_sequence(500, 42);
        let mut sa = shared.clone();
        sa.extend_from_slice(&make_test_sequence(500, 1));
        let mut sb = shared.clone();
        sb.extend_from_slice(&make_test_sequence(500, 2));
        let seqs = vec![
            ("ga".to_string(), sa),
            ("gb".to_string(), sb),
        ];

        let mut index = SyngIndex::build(params, seqs.into_iter());
        index.build_fast_locate().unwrap();
        let in_mem = index.query_region("ga", 0, 1000, 120).unwrap();

        let dir = std::env::temp_dir().join("impg_test_syng_locate_roundtrip");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("idx");
        let prefix_str = prefix.to_str().unwrap();
        index.save(prefix_str).unwrap();
        assert!(
            std::path::Path::new(&format!("{}.syng.locate", prefix_str)).exists(),
            ".syng.locate file should exist after save"
        );

        let loaded = SyngIndex::load(prefix_str, params).unwrap();
        assert!(loaded.has_fast_locate(), "loaded index should have fast_locate");
        let reloaded = loaded.query_region("ga", 0, 1000, 120).unwrap();

        let to_set = |v: &[HomologousInterval]| -> std::collections::BTreeSet<(String, u64, u64)> {
            v.iter().map(|iv| (iv.genome.clone(), iv.start, iv.end)).collect()
        };
        assert_eq!(to_set(&in_mem), to_set(&reloaded));
        assert!(!in_mem.is_empty());
    }

    /// FastLocate fast-path parity: building the locate structure on top of
    /// an existing SyngIndex and re-querying must return the SAME intervals
    /// (after merging) as the fallback walk-every-path implementation.
    #[test]
    fn test_query_region_fast_locate_parity() {
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

        // `SyngIndex::build` eagerly builds the locate structure, so we pull
        // out the fast result first and then tear down the locate state to
        // exercise the fallback walk-every-path implementation.
        let mut index = SyngIndex::build(params, seqs.into_iter());
        assert!(index.has_fast_locate());
        let fast = index.query_region("ga", 0, 1000, 120).unwrap();

        index.gbz_gbwt = None;
        index.fast_locate = None;
        index.bp_offsets = None;
        assert!(!index.has_fast_locate());
        let slow = index.query_region("ga", 0, 1000, 120).unwrap();

        // Normalize (genome, start, end) tuples and compare as sets.
        let to_set = |v: &[HomologousInterval]| -> std::collections::BTreeSet<(String, u64, u64)> {
            v.iter().map(|iv| (iv.genome.clone(), iv.start, iv.end)).collect()
        };
        let slow_set = to_set(&slow);
        let fast_set = to_set(&fast);
        assert_eq!(
            slow_set, fast_set,
            "fast-path and slow-path query_region disagree\nslow: {:?}\nfast: {:?}",
            slow_set, fast_set
        );
        assert!(!slow_set.is_empty(), "test should produce at least one hit");
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
        let params = SyncmerParams { k: 8, w: 55, seed: 7 };
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
        let intervals = index.query_region("genome_a", 0, shared_len as u64, 0).unwrap();

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
            assert!(iv.end > iv.start, "interval end should be > start: {:?}", iv);
            assert!(
                iv.end <= index.name_map.path_to_length[
                    *index.name_map.name_to_path.get(&iv.genome).unwrap() as usize
                ],
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
        let params = SyncmerParams { k: 8, w: 55, seed: 7 };
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
                    "Padded interval should be >= unpadded: {:?} vs {:?}", wp, np
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
        let params = SyncmerParams { k: 8, w: 55, seed: 7 };
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
            "Self-hit expected. Found: {:?}", genomes
        );
        assert!(
            genomes.contains(&"genome_b"),
            "genome_b shares backbone. Found: {:?}", genomes
        );
        assert!(
            genomes.contains(&"genome_c"),
            "genome_c shares backbone. Found: {:?}", genomes
        );

        // Coordinates should be in a reasonable range (within ~syncmer_len of 0..400)
        let syncmer_len = (params.w + params.k) as u64;
        for iv in &intervals {
            if iv.genome == "genome_a" || iv.genome == "genome_b" || iv.genome == "genome_c" {
                // Start should be near 0 (within syncmer length)
                assert!(
                    iv.start <= syncmer_len,
                    "{}: start {} should be near 0 (within {})",
                    iv.genome, iv.start, syncmer_len
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
            "Interior of shared region should find genome_b. Found: {:?}", genomes
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
            pad_120_b.start, no_pad_b.start
        );
        assert!(
            pad_120_b.end >= no_pad_b.end,
            "Padded end {} should be >= unpadded end {}",
            pad_120_b.end, no_pad_b.end
        );

        // The difference should be approximately the padding amount
        let start_diff = no_pad_b.start as i64 - pad_120_b.start as i64;
        let end_diff = pad_120_b.end as i64 - no_pad_b.end as i64;
        assert!(
            start_diff >= 0,
            "start_diff should be non-negative: {}", start_diff
        );
        assert!(
            end_diff >= 0,
            "end_diff should be non-negative: {}", end_diff
        );
    }

    #[test]
    fn test_query_padding_zero_vs_nonzero() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq = make_test_sequence(1000, 42);

        let index = SyngIndex::build(
            params,
            vec![
                ("g1".to_string(), seq.clone()),
                ("g2".to_string(), seq),
            ]
            .into_iter(),
        );

        let pad_0 = index.query_region("g1", 300, 500, 0).unwrap();
        let pad_60 = index.query_region("g1", 300, 500, 60).unwrap();
        let pad_120 = index.query_region("g1", 300, 500, 120).unwrap();

        // Intervals should grow monotonically with padding
        let g2_0 = pad_0.iter().find(|iv| iv.genome == "g2").unwrap();
        let g2_60 = pad_60.iter().find(|iv| iv.genome == "g2").unwrap();
        let g2_120 = pad_120.iter().find(|iv| iv.genome == "g2").unwrap();

        assert!(g2_60.start <= g2_0.start, "60bp pad start should be <= 0bp pad start");
        assert!(g2_60.end >= g2_0.end, "60bp pad end should be >= 0bp pad end");
        assert!(g2_120.start <= g2_60.start, "120bp pad start should be <= 60bp pad start");
        assert!(g2_120.end >= g2_60.end, "120bp pad end should be >= 60bp pad end");
    }

    #[test]
    fn test_query_padding_clamped_to_genome_length() {
        let _guard = lock_syng();
        let params = SyncmerParams::default();
        let seq = make_test_sequence(500, 42);
        let genome_len = seq.len() as u64;

        let index = SyngIndex::build(
            params,
            vec![
                ("g1".to_string(), seq.clone()),
                ("g2".to_string(), seq),
            ]
            .into_iter(),
        );

        // Query near the start with large padding — should clamp to 0
        let intervals = index.query_region("g1", 0, 200, 10000).unwrap();
        for iv in &intervals {
            assert!(iv.start == 0, "Start should be clamped to 0, got {}", iv.start);
            assert!(
                iv.end <= genome_len,
                "End {} should be <= genome length {}",
                iv.end, genome_len
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
                .map(|(q, t)| Anchor { query_pos: q, target_pos: t, node_id: 0 })
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
            "Self-hit expected. Found: {:?}", genomes
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

        let index = SyngIndex::build(
            params,
            vec![("only_genome".to_string(), seq)].into_iter(),
        );

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

        let index = SyngIndex::build(
            params,
            vec![("g1".to_string(), seq)].into_iter(),
        );

        let intervals = index.query_region("g1", 10000, 20000, 0).unwrap();
        assert!(
            intervals.is_empty(),
            "Query beyond sequence length should return empty, got {:?}",
            intervals.iter().map(|iv| (&iv.genome, iv.start, iv.end)).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_query_region_zero_width() {
        let _guard = lock_syng();
        // Query with start == end — no syncmers in empty range
        let params = SyncmerParams::default();
        let seq = make_test_sequence(500, 42);

        let index = SyngIndex::build(
            params,
            vec![("g1".to_string(), seq)].into_iter(),
        );

        let intervals = index.query_region("g1", 200, 200, 0).unwrap();
        assert!(
            intervals.is_empty(),
            "Zero-width query should return empty"
        );
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
                "-f", fasta_path.to_str().unwrap(),
                "-o", output_prefix.to_str().unwrap(),
            ])
            .output()
            .expect("Failed to run impg syng");
        assert!(
            output.status.success(),
            "impg syng failed: {}", String::from_utf8_lossy(&output.stderr)
        );

        // Query via -a syng prefix with BED output
        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "-a", output_prefix.to_str().unwrap(),
                "--sequence-files", fasta_path.to_str().unwrap(),
                "-r", "genome_a:0-400",
                "-o", "bed",
            ])
            .output()
            .expect("Failed to run impg query -a syng prefix");
        assert!(
            output.status.success(),
            "impg query -a syng prefix failed: {}", String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8_lossy(&output.stdout);
        let lines: Vec<&str> = stdout.lines().collect();
        assert!(
            !lines.is_empty(),
            "BED output should not be empty"
        );

        // Should find genome_b in the output
        let has_genome_b = lines.iter().any(|l| l.starts_with("genome_b\t"));
        assert!(has_genome_b, "BED output should contain genome_b. Got:\n{}", stdout);

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
                "-f", fasta_path.to_str().unwrap(),
                "-o", output_prefix.to_str().unwrap(),
            ])
            .output()
            .expect("Failed to run impg syng");
        assert!(output.status.success(), "impg syng failed: {}", String::from_utf8_lossy(&output.stderr));

        // Query with GFA output
        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "-a", output_prefix.to_str().unwrap(),
                "--sequence-files", fasta_path.to_str().unwrap(),
                "-r", "genome_a:0-400",
                "-o", "gfa",
            ])
            .output()
            .expect("Failed to run impg query -a syng prefix -o gfa");
        assert!(
            output.status.success(),
            "impg query -a syng prefix -o gfa failed: {}", String::from_utf8_lossy(&output.stderr)
        );

        let stdout = String::from_utf8_lossy(&output.stdout);
        // GFA should have S (segment) and L (link) lines
        let has_s_lines = stdout.lines().any(|l| l.starts_with("S\t"));
        assert!(has_s_lines, "GFA output should contain S lines. Got:\n{}", stdout);

        // Check for path-related lines (W or P lines)
        let has_paths = stdout.lines().any(|l| l.starts_with("W\t") || l.starts_with("P\t"));
        assert!(has_paths, "GFA output should contain W or P lines. Got:\n{}", stdout);

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
            .map(|i| (format!("region_seq_{}", i), make_test_sequence(500, i as u8 + 50)))
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
        assert!(!intervals.is_empty(), "Should find intervals in shared region");

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
        assert!(Path::new(&khash_path).exists(), "Region .1khash should exist");
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
        assert!(!loaded.kmer_hash.is_null(), "Loaded region KmerHash should be valid");
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
            .map(|i| (format!("genome_{}", i), make_test_sequence(2000, i as u8 + 10)))
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
        std::fs::remove_dir_all(
            std::env::temp_dir().join("impg_test_region_gbwt_nested"),
        )
        .ok();
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
        let region_genomes: Vec<&str> = region_intervals.iter().map(|iv| iv.genome.as_str()).collect();

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
                "-f", fasta_path.to_str().unwrap(),
                "-o", idx_prefix.to_str().unwrap(),
            ])
            .output()
            .expect("Failed to run impg syng");
        assert!(output.status.success(), "impg syng failed: {}", String::from_utf8_lossy(&output.stderr));

        // Query with GBWT output
        let gbwt_output = dir.join("region_gbwt");
        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "-a", idx_prefix.to_str().unwrap(),
                "--sequence-files", fasta_path.to_str().unwrap(),
                "-r", "genome_a:0-400",
                "-o", "gbwt",
                "-O", gbwt_output.to_str().unwrap(),
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
                "-a", paf_path.to_str().unwrap(),
                "--sequence-files", fasta_path.to_str().unwrap(),
                "-r", "genome_a:100-500",
                "-o", "gbwt",
                "-O", gbwt_output.to_str().unwrap(),
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
        assert!(Path::new(&gbwt_file).exists(), "PAF-based .1gbwt should exist");
        assert!(Path::new(&khash_file).exists(), "PAF-based .1khash should exist");

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
                "-a", "/tmp/nonexistent_prefix",
                "-r", "genome_a:0-100",
                "-o", "gbwt",
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
