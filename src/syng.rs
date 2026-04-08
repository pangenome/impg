//! Safe Rust wrapper around syng's C library.
//!
//! Provides `SyngIndex` for building, loading, saving, and querying
//! GBWT-based syncmer indices.

use rustc_hash::FxHashMap;
use std::ffi::CString;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::syng_ffi;

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
    pub(crate) fn to_c(&self) -> syng_ffi::CSyncmerParams {
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
    /// `path_number\tname\tlength\tstart_node\tstart_count\tnum_syncmers`
    /// (last three columns are 0 if path start info is unavailable)
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
            let (sn, sc, ns) = match start {
                Some(info) => (info.start_node, info.start_count, info.num_syncmers),
                None => (0, 0, 0),
            };
            writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}", i, name, length, sn, sc, ns)?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Load from a `.syng.names` file.
    /// Supports both old 3-column format and new 6-column format.
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

            // Parse optional path start info (columns 3-5)
            if parts.len() >= 6 {
                let start_node: i32 = parts[3].parse().unwrap_or(0);
                let start_count: u32 = parts[4].parse().unwrap_or(0);
                let num_syncmers: u32 = parts[5].parse().unwrap_or(0);
                if num_syncmers > 0 {
                    name_map.set_path_start(
                        path_num,
                        GbwtPathStart {
                            start_node,
                            start_count,
                            num_syncmers,
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
}

// SAFETY: The C structures are only accessed through &self or &mut self,
// and the C library's read-only query functions are thread-safe.
unsafe impl Send for SyngIndex {}

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
            syng_ffi::seqhashCreate(params.k as i32, params.w as i32, params.seed as i32)
        };
        Self {
            gbwt,
            kmer_hash,
            seqhash,
            name_map: SyngNameMap::new(),
            params,
        }
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

            // Convert sequence to null-terminated C string for the C library.
            // The seqhash expects ASCII DNA characters (A/C/G/T).
            let mut seq_buf: Vec<u8> = seq
                .iter()
                .map(|&b| match b {
                    b'a' | b'A' => b'a',
                    b'c' | b'C' => b'c',
                    b'g' | b'G' => b'g',
                    b't' | b'T' => b't',
                    // Numeric encoding from AGC (0=A, 1=C, 2=G, 3=T)
                    0 => b'a',
                    1 => b'c',
                    2 => b'g',
                    3 => b't',
                    _ => b'a', // N -> a (syng convention)
                })
                .collect();
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
                },
            );
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
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "Failed to create ONEcode schema for gbwt",
                ));
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
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    format!("Failed to open {} for writing", gbwt_path),
                ));
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
                return Err(io::Error::new(
                    io::ErrorKind::Other,
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
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    format!("Failed to open {} for writing", khash_path),
                ));
            }
            let ok = syng_ffi::kmerHashWriteOneFile(self.kmer_hash, of);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
            if !ok {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "kmerHashWriteOneFile failed",
                ));
            }
        }

        // Write .syng.names
        let names_path = format!("{}.syng.names", prefix);
        self.name_map.save(&names_path)?;

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
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "Failed to create ONEcode schema for gbwt read",
                ));
            }
            let of = syng_ffi::oneFileOpenRead(
                gbwt_cpath.as_ptr(),
                schema,
                gbwt_type.as_ptr(),
                1,
            );
            if of.is_null() {
                syng_ffi::oneSchemaDestroy(schema);
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    format!("Failed to open {} for reading", gbwt_path),
                ));
            }
            let gbwt = syng_ffi::syngBWTread(of);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
            if gbwt.is_null() {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "syngBWTread returned null",
                ));
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
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "Failed to create ONEcode schema for khash read",
                ));
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
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    format!("Failed to open {} for reading", khash_path),
                ));
            }
            let kh = syng_ffi::kmerHashReadOneFile(of);
            syng_ffi::oneFileClose(of);
            syng_ffi::oneSchemaDestroy(schema);
            if kh.is_null() {
                syng_ffi::syngBWTdestroy(gbwt);
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    "kmerHashReadOneFile returned null",
                ));
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
            syng_ffi::seqhashCreate(params.k as i32, params.w as i32, params.seed as i32)
        };

        Ok(Self {
            gbwt,
            kmer_hash,
            seqhash,
            name_map,
            params,
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

    /// Walk a genome's forward GBWT path, returning (node_id, accumulated_position)
    /// for each syncmer in the path.
    fn walk_path(&self, start: &GbwtPathStart) -> Vec<(i32, u64)> {
        let mut nodes = Vec::with_capacity(start.num_syncmers as usize);
        unsafe {
            let sbp =
                syng_ffi::syngBWTpathStartOld(self.gbwt, start.start_node, start.start_count);
            // First node position is 0 (relative — actual genomic position comes from syncmer pos)
            // We track accumulated offset from path start
            let mut acc_pos: u64 = 0;
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

        // 3. Walk all genomes' forward paths, find overlapping nodes
        let num_genomes = self.name_map.path_to_name.len();
        let mut all_intervals: Vec<HomologousInterval> = Vec::new();

        for genome_idx in 0..num_genomes {
            let path_start = match &self.name_map.path_starts[genome_idx] {
                Some(ps) => ps,
                None => continue,
            };

            let genome_name = &self.name_map.path_to_name[genome_idx];
            let genome_len = self.name_map.path_to_length[genome_idx];
            let nodes = self.walk_path(path_start);

            // Find min/max positions of matching nodes on this genome
            let mut min_pos: Option<u64> = None;
            let mut max_pos: Option<u64> = None;

            for &(node_id, pos) in &nodes {
                if query_node_set.contains(&node_id.unsigned_abs()) {
                    min_pos = Some(min_pos.map_or(pos, |m: u64| m.min(pos)));
                    max_pos = Some(max_pos.map_or(pos, |m: u64| m.max(pos)));
                }
            }

            if let (Some(min_p), Some(max_p)) = (min_pos, max_pos) {
                // Add syncmer length to max_pos to get the end of the last syncmer
                let syncmer_len = (self.params.w + self.params.k) as u64;
                let interval_end = max_p + syncmer_len;

                // Apply padding
                let padded_start = min_p.saturating_sub(padding);
                let padded_end = (interval_end + padding).min(genome_len);

                all_intervals.push(HomologousInterval {
                    genome: genome_name.clone(),
                    start: padded_start,
                    end: padded_end,
                    strand: '+',
                });
            }
        }

        // 4. Merge overlapping intervals per genome
        Self::merge_intervals(&mut all_intervals);

        Ok(all_intervals)
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

    // ── 1. FFI smoke tests ──────────────────────────────────────────

    #[test]
    fn test_syng_ffi_seqhash_create_destroy() {
        // Create a Seqhash with default params (k=8, w=55, seed=7), verify non-null, free it.
        let sh = unsafe { syng_ffi::seqhashCreate(8, 55, 7) };
        assert!(!sh.is_null(), "seqhashCreate returned null for valid params");
        unsafe { syng_ffi::impg_seqhashDestroy(sh) };
    }

    #[test]
    fn test_syng_ffi_kmerhash_create_destroy() {
        // Create a KmerHash, verify non-null, destroy it.
        let syncmer_len = 55 + 8; // w + k = 63
        let kh = unsafe { syng_ffi::kmerHashCreate(1024, syncmer_len) };
        assert!(!kh.is_null(), "kmerHashCreate returned null");
        unsafe { syng_ffi::kmerHashDestroy(kh) };
    }

    #[test]
    fn test_syng_ffi_syngbwt_create_destroy() {
        // Create a SyngBWT, verify non-null, destroy it.
        let syncmer_len = 63;
        let sb = unsafe { syng_ffi::syngBWTcreate(syncmer_len, 0) };
        assert!(!sb.is_null(), "syngBWTcreate returned null");
        unsafe { syng_ffi::syngBWTdestroy(sb) };
    }

    // ── 2. SyncmerParams defaults ───────────────────────────────────

    #[test]
    fn test_syncmer_params_default() {
        let params = SyncmerParams::default();
        assert_eq!(params.k, 8);
        assert_eq!(params.w, 55);
        assert_eq!(params.seed, 7);
    }

    #[test]
    fn test_syncmer_params_to_c() {
        let params = SyncmerParams::default();
        let c_params = params.to_c();
        assert_eq!(c_params.k, 8);
        assert_eq!(c_params.w, 55);
        assert_eq!(c_params.seed, 7);
    }

    // ── 3. SyngIndex lifecycle ──────────────────────────────────────

    #[test]
    fn test_syng_index_create_drop() {
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
        let nm = SyngNameMap::new();
        assert!(nm.path_to_name.is_empty());
        assert!(nm.path_to_length.is_empty());
        assert!(nm.name_to_path.is_empty());
    }

    #[test]
    fn test_name_map_default() {
        // Default trait should produce the same as new()
        let nm = SyngNameMap::default();
        assert!(nm.path_to_name.is_empty());
        assert!(nm.path_to_length.is_empty());
        assert!(nm.name_to_path.is_empty());
    }

    #[test]
    fn test_name_map_add() {
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
        let index = SyngIndex::build(SyncmerParams::default(), std::iter::empty());
        assert!(index.name_map.path_to_name.is_empty());
    }

    #[test]
    fn test_build_short_sequence() {
        // Sequence shorter than syncmer length (63) — should be recorded but not panick
        let seqs = vec![("short".to_string(), b"ACGTACGT".to_vec())];
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 1);
        assert_eq!(index.name_map.path_to_name[0], "short");
    }

    #[test]
    fn test_build_single_sequence() {
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
        // A 0-length sequence should be handled gracefully (no panic)
        let seqs = vec![("empty".to_string(), Vec::new())];
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 1);
        assert_eq!(index.name_map.path_to_name[0], "empty");
        assert_eq!(index.name_map.path_to_length[0], 0);
    }

    #[test]
    fn test_syng_build_sequence_shorter_than_syncmer() {
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
        // Sequence exactly syncmer length (63bp) — borderline, should not panic
        let seq = make_test_sequence(63, 88);
        let seqs = vec![("exact".to_string(), seq)];
        let index = SyngIndex::build(SyncmerParams::default(), seqs.into_iter());
        assert_eq!(index.name_map.path_to_name.len(), 1);
        assert_eq!(index.name_map.path_to_length[0], 63);
    }

    #[test]
    fn test_syng_build_single_path_in_gbwt() {
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
        let result = SyngIndex::load("/tmp/nonexistent_prefix_xyz123", SyncmerParams::default());
        assert!(result.is_err(), "Loading from nonexistent files should fail");
    }

    // ── 13. query_region tests ─────────────────────────────────────

    /// Build an index from sequences that share some content, then query a region.
    #[test]
    fn test_query_region_basic() {
        // Create sequences that share a common prefix/suffix (to get shared syncmer nodes)
        // Use the same random seed for a shared backbone of ~500bp,
        // then diverge in the middle.
        let shared_len = 500;
        let total_len = 1000;
        let params = SyncmerParams { k: 8, w: 55, seed: 7 };
        let syncmer_len = (params.w + params.k) as usize;

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
        let params = SyncmerParams::default();
        let seq = make_test_sequence(500, 42);
        let index = SyngIndex::build(params, vec![("seq1".to_string(), seq)].into_iter());

        let result = index.query_region("nonexistent", 0, 100, 0);
        assert!(result.is_err(), "Should error for unknown genome");
    }

    #[test]
    fn test_query_region_padding() {
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
}
