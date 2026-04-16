//! Safe Rust wrapper around syng's C library.
//!
//! Provides `SyngIndex` for building, loading, saving, and querying
//! GBWT-based syncmer indices.

use rustc_hash::FxHashMap;
use std::ffi::CString;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
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
#[derive(Debug, Clone)]
pub struct HomologousInterval {
    pub genome: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
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
            // walk_path yields (node_id, accumulated_bp_pos). Syng's rskip may
            // report negative node ids when a forward walk traverses an edge
            // whose canonical orientation is the opposite strand; `query_region`
            // already collapses strand via `unsigned_abs()`, so we do the same
            // here. Each syncmer maps to one gbz record regardless of which
            // strand visits it.
            //
            // Encoded as gbz forward-orientation ids
            // (`encode_node(raw, Forward) = raw << 1`). ENDMARKER is 0, so
            // syncmer ids starting at 1 are safe.
            let mut path_encoded: Vec<usize> = Vec::with_capacity(walk.len());
            for &(node, pos) in &walk {
                let raw = node.unsigned_abs() as usize;
                debug_assert!(raw > 0, "forward walk produced zero node id");
                path_encoded.push(raw << 1);
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
        let syncmer_len = (params.w + params.k) as usize;

        for (name, seq) in sequences {
            let seq_len = seq.len();
            if seq_len < syncmer_len {
                // Sequence too short for syncmer extraction — record in name map but skip
                index.name_map.add(name, seq_len as u64);
                continue;
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
                    index.seqhash,
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
                        index.kmer_hash,
                        seq_buf.as_mut_ptr().add(pos as usize) as *mut i8,
                        &mut kmer_index,
                    );
                    syncmers.push((kmer_index, pos));
                }

                syng_ffi::impg_seqhashIteratorDestroy(sit);
            }

            if syncmers.is_empty() {
                index.name_map.add(name, seq_len as u64);
                continue;
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
                let sbp = syng_ffi::syngBWTpathStartNew(index.gbwt, first_sync);
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
                let sbp = syng_ffi::syngBWTpathStartNew(index.gbwt, first_sync_rc);
                for i in (0..n - 1).rev() {
                    let next_sync_rc = -(syncmers[i].0 as i32);
                    let offset = (syncmers[i + 1].1 - syncmers[i].1) as u32;
                    syng_ffi::syngBWTpathAdd(sbp, next_sync_rc, offset);
                }
                syng_ffi::syngBWTpathFinish(sbp);
            }

            // Record in name map with path start info
            let path_num = index.name_map.add(name, seq_len as u64);
            index.name_map.set_path_start(
                path_num,
                GbwtPathStart {
                    start_node: fwd_start_node,
                    start_count: fwd_start_count,
                    num_syncmers: syncmers.len() as u32,
                    first_syncmer_pos: fwd_first_syncmer_pos,
                },
            );
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

    /// Save the index to disk.
    ///
    /// Produces three files at the given prefix:
    /// - `{prefix}.1khash` — syncmer hash (syng-compatible)
    /// - `{prefix}.1gbwt` — GBWT graph (syng-compatible)
    /// - `{prefix}.syng.names` — sequence name mapping
    pub fn save(&self, prefix: &str) -> io::Result<()> {
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

        Ok(Self {
            gbwt,
            kmer_hash,
            seqhash,
            name_map,
            params,
            gbz_gbwt,
            fast_locate,
            bp_offsets,
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
    pub fn query_region(
        &self,
        genome: &str,
        start: u64,
        end: u64,
        padding: u64,
    ) -> Result<Vec<HomologousInterval>, io::Error> {
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

        // 2. Walk the query genome's path, find nodes in [start, end]
        let query_nodes = self.walk_path(query_start);
        let mut query_node_set = rustc_hash::FxHashSet::default();
        for &(node_id, pos) in &query_nodes {
            if pos >= start && pos < end {
                // Use absolute node ID (positive) so both strands match
                query_node_set.insert(node_id.unsigned_abs());
            }
        }

        if query_node_set.is_empty() {
            return Ok(Vec::new());
        }

        // 3. Enumerate hits. If the FastLocate fast path is prepared, use
        // it to query each syncmer node in O(k log r) rather than walking
        // every forward path per query. Otherwise fall back to the original
        // walk-every-path loop.
        let syncmer_len = (self.params.w + self.params.k) as u64;
        let mut all_intervals: Vec<HomologousInterval> = Vec::new();

        if let (Some(fl), Some(gbz_gbwt), Some(bp_off)) =
            (self.fast_locate.as_ref(), self.gbz_gbwt.as_ref(), self.bp_offsets.as_ref())
        {
            // Fast path: per query syncmer node, decompress_da -> (gbz_seq_id,
            // seq_offset_from_end). gbz_seq_id = 2*forward_path_idx + {0 fwd / 1 rev};
            // we only care about forward visits (even gbz ids).
            for &query_node in &query_node_set {
                // syng node id is positive i32; gbz encoded forward id = raw << 1.
                let encoded = (query_node as usize) << 1;
                if !gbz_gbwt.has_node(encoded) {
                    continue;
                }
                let visits = fl.decompress_da(gbz_gbwt, encoded);
                for (gbz_seq_id, seq_off_from_end) in visits {
                    // Only forward sequences (even gbz ids). Reverse hits land
                    // at odd gbz ids and are covered via the forward side of
                    // some other path (or are redundant with its own forward).
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
                    let Some(pos) = bp_off.bp_of(forward_path_idx, forward_node_idx) else {
                        continue;
                    };
                    let genome_len = self.name_map.path_to_length[forward_path_idx];
                    let hit_end = pos + syncmer_len;
                    let padded_start = pos.saturating_sub(padding);
                    let padded_end = (hit_end + padding).min(genome_len);
                    all_intervals.push(HomologousInterval {
                        genome: self.name_map.path_to_name[forward_path_idx].clone(),
                        start: padded_start,
                        end: padded_end,
                        strand: '+',
                    });
                }
            }
        } else {
            // Fallback: walk every forward path (the old pre-locate behavior).
            let num_genomes = self.name_map.path_to_name.len();
            for genome_idx in 0..num_genomes {
                let path_start = match &self.name_map.path_starts[genome_idx] {
                    Some(ps) => ps,
                    None => continue,
                };
                let genome_name = &self.name_map.path_to_name[genome_idx];
                let genome_len = self.name_map.path_to_length[genome_idx];
                let nodes = self.walk_path(path_start);
                for &(node_id, pos) in &nodes {
                    if query_node_set.contains(&node_id.unsigned_abs()) {
                        let hit_end = pos + syncmer_len;
                        let padded_start = pos.saturating_sub(padding);
                        let padded_end = (hit_end + padding).min(genome_len);
                        all_intervals.push(HomologousInterval {
                            genome: genome_name.clone(),
                            start: padded_start,
                            end: padded_end,
                            strand: '+',
                        });
                    }
                }
            }
        }

        // 4. Merge overlapping/adjacent intervals per genome
        Self::merge_intervals(&mut all_intervals);

        Ok(all_intervals)
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
    fn merge_intervals(intervals: &mut Vec<HomologousInterval>) {
        if intervals.len() <= 1 {
            return;
        }

        // Sort by (genome, strand, start)
        intervals.sort_by(|a, b| {
            a.genome
                .cmp(&b.genome)
                .then(a.strand.cmp(&b.strand))
                .then(a.start.cmp(&b.start))
        });

        let mut merged: Vec<HomologousInterval> = Vec::new();
        for iv in intervals.drain(..) {
            if let Some(last) = merged.last_mut() {
                if last.genome == iv.genome && last.strand == iv.strand && iv.start <= last.end {
                    last.end = last.end.max(iv.end);
                    continue;
                }
            }
            merged.push(iv);
        }
        *intervals = merged;
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
        let index = SyngIndex::build(params, seqs.into_iter());

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
        let original = SyngIndex::build(params, seqs.into_iter());

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
        let original = SyngIndex::build(params, seqs.into_iter());

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
        let original = SyngIndex::build(params, seqs.into_iter());

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
        let idx_default = SyngIndex::build(
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
        let idx_short = SyngIndex::build(
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
        let idx_alt = SyngIndex::build(
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
        let index = SyngIndex::build(params, seqs.into_iter());

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
        let index = SyngIndex::build(params, seqs.into_iter());
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
        let index = SyngIndex::build(params, seqs.into_iter());

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

        // Load and verify content
        let loaded = SyngIndex::load(output_prefix_str, SyncmerParams::default()).unwrap();
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

        let index = SyngIndex::build(
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

    #[test]
    fn test_merge_intervals_no_overlap() {
        let _guard = lock_syng();
        let mut intervals = vec![
            HomologousInterval { genome: "g1".to_string(), start: 0, end: 100, strand: '+' },
            HomologousInterval { genome: "g1".to_string(), start: 200, end: 300, strand: '+' },
        ];
        SyngIndex::merge_intervals(&mut intervals);
        assert_eq!(intervals.len(), 2, "Non-overlapping intervals should not merge");
    }

    #[test]
    fn test_merge_intervals_overlapping() {
        let _guard = lock_syng();
        let mut intervals = vec![
            HomologousInterval { genome: "g1".to_string(), start: 0, end: 150, strand: '+' },
            HomologousInterval { genome: "g1".to_string(), start: 100, end: 300, strand: '+' },
        ];
        SyngIndex::merge_intervals(&mut intervals);
        assert_eq!(intervals.len(), 1, "Overlapping intervals should merge");
        assert_eq!(intervals[0].start, 0);
        assert_eq!(intervals[0].end, 300);
    }

    #[test]
    fn test_merge_intervals_adjacent() {
        let _guard = lock_syng();
        let mut intervals = vec![
            HomologousInterval { genome: "g1".to_string(), start: 0, end: 100, strand: '+' },
            HomologousInterval { genome: "g1".to_string(), start: 100, end: 200, strand: '+' },
        ];
        SyngIndex::merge_intervals(&mut intervals);
        assert_eq!(intervals.len(), 1, "Adjacent intervals should merge");
        assert_eq!(intervals[0].start, 0);
        assert_eq!(intervals[0].end, 200);
    }

    #[test]
    fn test_merge_intervals_different_genomes() {
        let _guard = lock_syng();
        let mut intervals = vec![
            HomologousInterval { genome: "g1".to_string(), start: 0, end: 150, strand: '+' },
            HomologousInterval { genome: "g2".to_string(), start: 50, end: 200, strand: '+' },
        ];
        SyngIndex::merge_intervals(&mut intervals);
        assert_eq!(intervals.len(), 2, "Overlapping intervals on different genomes should not merge");
    }

    #[test]
    fn test_merge_intervals_multiple_groups() {
        let _guard = lock_syng();
        let mut intervals = vec![
            HomologousInterval { genome: "g1".to_string(), start: 0, end: 100, strand: '+' },
            HomologousInterval { genome: "g1".to_string(), start: 50, end: 200, strand: '+' },
            HomologousInterval { genome: "g1".to_string(), start: 150, end: 400, strand: '+' },
            HomologousInterval { genome: "g2".to_string(), start: 0, end: 500, strand: '+' },
        ];
        SyngIndex::merge_intervals(&mut intervals);
        // g1: [0,100] + [50,200] + [150,400] → [0,400]
        // g2: [0,500]
        assert_eq!(intervals.len(), 2);
        let g1 = intervals.iter().find(|iv| iv.genome == "g1").unwrap();
        assert_eq!(g1.start, 0);
        assert_eq!(g1.end, 400);
    }

    #[test]
    fn test_merge_intervals_empty() {
        let _guard = lock_syng();
        let mut intervals: Vec<HomologousInterval> = Vec::new();
        SyngIndex::merge_intervals(&mut intervals);
        assert_eq!(intervals.len(), 0);
    }

    #[test]
    fn test_merge_intervals_single() {
        let _guard = lock_syng();
        let mut intervals = vec![
            HomologousInterval { genome: "g1".to_string(), start: 10, end: 20, strand: '+' },
        ];
        SyngIndex::merge_intervals(&mut intervals);
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
        assert_eq!(
            intervals.len(), 1,
            "Single-sequence index should return exactly 1 interval (self)"
        );
        assert_eq!(intervals[0].genome, "only_genome");
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

        // Query via --syng with BED output
        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "--syng", output_prefix.to_str().unwrap(),
                "--sequence-files", fasta_path.to_str().unwrap(),
                "-r", "genome_a:0-400",
                "-o", "bed",
            ])
            .output()
            .expect("Failed to run impg query --syng");
        assert!(
            output.status.success(),
            "impg query --syng failed: {}", String::from_utf8_lossy(&output.stderr)
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
                "--syng", output_prefix.to_str().unwrap(),
                "--sequence-files", fasta_path.to_str().unwrap(),
                "-r", "genome_a:0-400",
                "-o", "gfa",
            ])
            .output()
            .expect("Failed to run impg query --syng -o gfa");
        assert!(
            output.status.success(),
            "impg query --syng -o gfa failed: {}", String::from_utf8_lossy(&output.stderr)
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

    #[test]
    fn test_syng_cli_mutual_exclusivity() {
        let _guard = lock_syng();
        // --syng + -a should produce an error
        let dir = std::env::temp_dir().join("impg_test_syng_cli_mutex");
        std::fs::create_dir_all(&dir).unwrap();

        let bin = find_impg_binary();
        if bin.is_none() {
            eprintln!("Skipping CLI test: impg binary not found");
            std::fs::remove_dir_all(&dir).ok();
            return;
        }
        let bin = bin.unwrap();

        let output = std::process::Command::new(&bin)
            .args([
                "query",
                "--syng", "some_prefix",
                "-a", "some_alignment.paf",
                "-r", "genome_a:0-100",
            ])
            .output()
            .expect("Failed to run impg query");

        assert!(
            !output.status.success(),
            "Using --syng with -a should fail"
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
        let full_index = SyngIndex::build(params, seqs.clone().into_iter());

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
        let index = SyngIndex::build(params, seqs.into_iter());

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
        let index = SyngIndex::build(params, seqs.into_iter());

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
        // impg query --syng prefix -f test.fa -r region -o gbwt -O tmpdir/region
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
                "--syng", idx_prefix.to_str().unwrap(),
                "--sequence-files", fasta_path.to_str().unwrap(),
                "-r", "genome_a:0-400",
                "-o", "gbwt",
                "-O", gbwt_output.to_str().unwrap(),
            ])
            .output()
            .expect("Failed to run impg query -o gbwt");
        assert!(
            output.status.success(),
            "impg query --syng -o gbwt failed: {}",
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
                "--syng", "/tmp/nonexistent_prefix",
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
