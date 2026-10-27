//! THE EVIDENCE STORE — read evidence addressable by panel position.
//!
//! The projection emits reads in ARRIVAL order (`shard.<r>.bin.gz`), so
//! nothing can be seeked: the cut decompresses a whole region shard per
//! territory (~46 GB across all shards) to use ~0.03% of the reads. The
//! store is the missing layer, with BAM's own trick — POSITION ORDER
//! EQUALS FILE ORDER, so a sparse manifest is all the index needed.
//!
//! Two panel-side/sample-side artifacts, with different lifetimes:
//!
//! * `node -> (path, offset)` (`IMPG_PANEL_POSITIONS`, from
//!   `examples/panel_positions.rs`): ONE per panel, ~4 GB, reused by
//!   every sample. Keyed on `|node|` — a vertex's position is
//!   orientation-independent, and the panel walk emits most vertices
//!   only forward, so keying on the signed id silently dropped every
//!   reverse-strand read (measured: 45% placement before the fix,
//!   83% after, and 0 walk-bearing reads lost).
//!
//! * the store itself (`IMPG_EVIDENCE_STORE=<prefix>`): ONE per sample,
//!   22.4 GB with per-bucket zstd (168 GB raw; buckets are the seek unit
//!   and the median is 2.3 KB, so decompression is microseconds). A read
//!   is written into the bucket of EVERY walk node's position, so "does a
//!   walk touch node N" is answered by reading N's bucket — structurally
//!   exact, no window geometry and no rewind.
//!
//! GATE (measured): fetching a territory's node set and feeding it to the
//! cut's own scan produced `cut-census.jsonl` and `cut-cache.bin` BYTE-
//! IDENTICAL to the full-shard scan, at 307 s against 637-670 s — and the
//! scan's own counters matched (276,042 kept, 5,988 touching patterns).

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

const MAGIC_RAW: &[u8; 8] = b"EVDIDX01";
const MAGIC_ZSTD: &[u8; 8] = b"EVDIDX03";
const POS_MAGIC: &[u8; 8] = b"PNLPOS01";

/// `node -> (path, offset)`, sorted by `|node|`.
///
/// MMAPPED, not loaded: the table is 4.1 GB with 271M entries, and the
/// cut runs once per territory in its own process — reading it into a Vec
/// (and re-sorting it, though it is already sorted) would cost ~60s and
/// 8 GB on EVERY territory, eating the whole win. A binary search over
/// the mapping is a page-cache read.
pub struct PanelPositions {
    map: memmap2::Mmap,
    n: usize,
}

impl PanelPositions {
    pub fn load(path: &Path) -> io::Result<PanelPositions> {
        let file = File::open(path)?;
        let map = unsafe { memmap2::Mmap::map(&file)? };
        if map.len() < 40 || &map[..8] != POS_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{}: not a panel-positions table", path.display()),
            ));
        }
        let n = u64::from_le_bytes(map[8..16].try_into().unwrap()) as usize;
        let n = n.min(map.len().saturating_sub(40) / 16);
        Ok(PanelPositions { map, n })
    }

    pub fn len(&self) -> usize {
        self.n
    }

    #[inline]
    fn node_at(&self, index: usize) -> i32 {
        let at = 40 + index * 16;
        i32::from_le_bytes(self.map[at..at + 4].try_into().unwrap())
    }

    #[inline]
    fn entry_at(&self, index: usize) -> (u32, u64) {
        let at = 40 + index * 16;
        let path_idx = u32::from_le_bytes(self.map[at + 4..at + 8].try_into().unwrap());
        let offset = u64::from_le_bytes(self.map[at + 8..at + 16].try_into().unwrap());
        (path_idx, offset)
    }

    /// The panel position of a node, or `None` when the table does not
    /// cover it (the panel's own walk did not reach that vertex).
    pub fn position(&self, node: i32) -> Option<(u32, u64)> {
        // keyed on |node|: a vertex's position is orientation-independent
        let key = node.unsigned_abs() as i32;
        let (mut lo, mut hi) = (0usize, self.n);
        while lo < hi {
            let mid = (lo + hi) / 2;
            match self.node_at(mid).cmp(&key) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => return Some(self.entry_at(mid)),
            }
        }
        None
    }
}

/// The sample's position-ordered read evidence.
pub struct EvidenceStore {
    payload: File,
    payload_base: u64,
    width_bp: u64,
    compressed: bool,
    /// (bucket key, payload offset, payload length) - ascending by key.
    table: Vec<(u64, u64, u64)>,
    /// The read ordinal travels with each packet: a read whose walk spans
    /// two buckets is stored in BOTH, so a multi-bucket query must dedup.
    pub fetched_packets: u64,
    pub buckets_read: u64,
    pub bytes_read: u64,
}

/// The store's key: `(path, offset / width)` packed into a u64.
pub fn position_key(path_idx: u32, offset: u64, width_bp: u64) -> u64 {
    ((path_idx as u64) << 24) | (offset / width_bp.max(1))
}

impl EvidenceStore {
    /// Open `<prefix>.evd` + `<prefix>.evd.idx`.
    pub fn open(prefix: &str) -> io::Result<EvidenceStore> {
        let idx_path = format!("{prefix}.evd.idx");
        let mut head = Vec::new();
        File::open(&idx_path)?.read_to_end(&mut head)?;
        if head.len() < 40 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{idx_path}: truncated manifest"),
            ));
        }
        let compressed = if &head[..8] == MAGIC_ZSTD {
            true
        } else if &head[..8] == MAGIC_RAW {
            false
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{idx_path}: bad magic"),
            ));
        };
        let u64_at = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().unwrap());
        let width_bp = u64_at(8);
        let n_buckets = u64_at(16) as usize;
        let mut table = Vec::with_capacity(n_buckets);
        for i in 0..n_buckets {
            let base = 40 + i * 24; // magic(8) + width, n_buckets, n_reads, rewind (4 x u64)
            if base + 24 > head.len() {
                break;
            }
            table.push((u64_at(base), u64_at(base + 8), u64_at(base + 16)));
        }
        let payload = File::open(format!("{prefix}.evd"))?;
        Ok(EvidenceStore {
            payload,
            // THE PAYLOAD IS ITS OWN FILE: the header and the bucket table
            // live in `.evd.idx`, so `<prefix>.evd` starts at offset 0 and
            // every table offset is relative to it. (The single-file
            // prototype had header+table+payload together, which is where
            // the wrong 40 came from — the gate caught it as a zstd
            // "Unknown frame descriptor".)
            payload_base: 0,
            width_bp,
            compressed,
            table,
            fetched_packets: 0,
            buckets_read: 0,
            bytes_read: 0,
        })
    }

    pub fn width_bp(&self) -> u64 {
        self.width_bp
    }

    pub fn bucket_count(&self) -> usize {
        self.table.len()
    }

    /// The bucket keys covering a set of nodes, ascending and deduplicated.
    /// A node the position table does not cover contributes nothing (the
    /// panel-coverage gap, measured at 0.02%).
    pub fn buckets_for_nodes(
        &self,
        nodes: impl IntoIterator<Item = i32>,
        positions: &PanelPositions,
    ) -> (Vec<u64>, usize) {
        let mut keys = Vec::new();
        let mut unresolved = 0usize;
        for node in nodes {
            match positions.position(node) {
                Some((path_idx, offset)) => keys.push(position_key(path_idx, offset, self.width_bp)),
                None => unresolved += 1,
            }
        }
        keys.sort_unstable();
        keys.dedup();
        (keys, unresolved)
    }

    /// Fetch every packet in the given buckets, deduplicated by read
    /// ordinal and sorted by it — so the scan sees the reads in their
    /// ARRIVAL order, which is what keeps its output byte-identical.
    pub fn packets_for_buckets(&mut self, keys: &[u64]) -> io::Result<Vec<Vec<u8>>> {
        let mut seen: HashSet<u32> = HashSet::new();
        let mut packets: Vec<(u32, Vec<u8>)> = Vec::new();
        let mut buf = Vec::new();
        for &key in keys {
            let Ok(at) = self.table.binary_search_by_key(&key, |e| e.0) else {
                continue;
            };
            let (_, offset, len) = self.table[at];
            self.payload
                .seek(SeekFrom::Start(self.payload_base + offset))?;
            buf.resize(len as usize, 0);
            self.payload.read_exact(&mut buf)?;
            self.bytes_read += len;
            self.buckets_read += 1;
            let decoded;
            let raw: &[u8] = if self.compressed {
                decoded = zstd::stream::decode_all(std::io::Cursor::new(&buf[..]))?;
                &decoded
            } else {
                &buf[..]
            };
            let mut p = 0usize;
            while p + 8 <= raw.len() {
                let ordinal = u32::from_le_bytes(raw[p..p + 4].try_into().unwrap());
                let plen = u32::from_le_bytes(raw[p + 4..p + 8].try_into().unwrap()) as usize;
                let end = p + 8 + plen;
                if end > raw.len() {
                    break;
                }
                if seen.insert(ordinal) {
                    packets.push((ordinal, raw[p + 4..end].to_vec()));
                }
                p = end;
            }
        }
        packets.sort_unstable_by_key(|p| p.0);
        self.fetched_packets += packets.len() as u64;
        Ok(packets.into_iter().map(|(_, packet)| packet).collect())
    }
}

/// Stream packets from a raw shard (the committed source) — a thin reader
/// so the store path can hand the cut the SAME `Vec<Vec<u8>>` batches.
pub fn read_shard_packets(path: &Path) -> io::Result<Vec<Vec<u8>>> {
    let file = File::open(path)?;
    let gz = flate2::read::GzDecoder::new(BufReader::with_capacity(1 << 20, file));
    let mut reader = BufReader::with_capacity(1 << 20, gz);
    let mut packets = Vec::new();
    loop {
        let mut len_buf = [0u8; 4];
        match reader.read_exact(&mut len_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        }
        let len = u32::from_le_bytes(len_buf) as usize;
        let mut packet = vec![0u8; 4 + len];
        packet[..4].copy_from_slice(&len_buf);
        reader.read_exact(&mut packet[4..])?;
        packets.push(packet);
    }
    Ok(packets)
}
