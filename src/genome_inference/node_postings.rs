//! NODE -> READ POSTINGS (`.n2r`) — an ADDITIVE index beside the evidence store.
//!
//! What the store already answers: "which reads have a walk node whose panel
//! position lands in bucket B", by seeking bucket B. A bucket is a locus of
//! `width_bp` (~10 kb) panel positions, NOT a node, so the store cannot answer
//! "which reads touch node N" without reading N's whole bucket and re-testing
//! every record's walk. That difference is the whole reason this artifact
//! exists: a territory's node set resolves to a few thousand buckets holding
//! ~19M packets of which ~1.5% touch the territory (measured on T00048), and
//! the cut decodes and re-derives all of them.
//!
//! The store's per-bucket records DO carry the read's ordinal (a `u32` before
//! each packet, the store's own global arrival id), so read ids are
//! recoverable — but only by decompressing the bucket and reading the record
//! headers. The walk nodes needed to attribute a read to a NODE live in the
//! packet payload, so an exact node -> read index costs one full pass over the
//! store (22.4 GB read, 168 GB decompressed) and cannot be derived from
//! `<prefix>.evd.idx` alone.
//!
//! ```text
//! node_postings build <prefix>     # <prefix>.evd{,.idx} -> <prefix>.n2r
//! node_postings summary <n2r>      # header
//! node_postings postings <n2r> N   # one node's reads
//! ```
//!
//! LAYOUT (version 1), little-endian:
//!
//! ```text
//! magic "N2R1" (4) | version u32 | node_count u64 | postings_count u64
//! node_count x (node_key u64, offset u64, len u32)      # ascending node_key
//! postings blob: per node, reads as delta-varint u64s   # offsets are into the blob
//! ```
//!
//! `node_key` is `|node|`: a vertex's identity is orientation-independent, the
//! same fold the panel-positions table and the cut's territory set use.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

pub const MAGIC: &[u8; 4] = b"N2R1";
pub const VERSION: u32 = 1;

const HEADER_BYTES: usize = 24; // magic(4) + version(4) + node_count(8) + postings_count(8)
const DIRECTORY_ENTRY_BYTES: usize = 20; // node(8) + offset(8) + len(4)
const SPILL_RECORD_BYTES: usize = 12; // node(8) + read(4)
/// The store's postings are 22.4 GB compressed; one pass emits every
/// (node, read) incidence, far more than memory holds, so pass 1 partitions by
/// node id and pass 2 sorts one partition at a time. The partition count only
/// bounds peak memory and temporary disk.
const SPILLS: u64 = 256;

/// The evidence store's own manifest magic. Private in `evidence_store`, and
/// that file's format is not ours to touch, so the layout is mirrored here.
const EVD_IDX_ZSTD: &[u8; 8] = b"EVDIDX03";
const EVD_IDX_RAW: &[u8; 8] = b"EVDIDX01";

/// `<prefix>.evd.idx` — the fields the artifact build needs from the store.
///
/// magic(8) + width_bp u64 + n_buckets u64 + n_reads u64 + rewind u64, then
/// n_buckets x (bucket_key u64, payload_offset u64, payload_len u64) ascending,
/// where `bucket_key = (path << 24) | (offset / width_bp)`.
pub struct StoreManifest {
    pub width_bp: u64,
    pub n_reads: u64,
    pub compressed: bool,
    /// (bucket key, offset into `<prefix>.evd`, payload length).
    pub table: Vec<(u64, u64, u64)>,
}

impl StoreManifest {
    pub fn read(prefix: &str) -> io::Result<StoreManifest> {
        let idx_path = format!("{prefix}.evd.idx");
        let mut head = Vec::new();
        File::open(&idx_path)?.read_to_end(&mut head)?;
        if head.len() < 40 {
            return Err(bad(&format!("{idx_path}: truncated manifest")));
        }
        let compressed = if &head[..8] == EVD_IDX_ZSTD {
            true
        } else if &head[..8] == EVD_IDX_RAW {
            false
        } else {
            return Err(bad(&format!("{idx_path}: not an evidence-store manifest")));
        };
        let u64_at = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().unwrap());
        let width_bp = u64_at(8);
        let n_buckets = u64_at(16) as usize;
        let n_reads = u64_at(24);
        let mut table = Vec::with_capacity(n_buckets);
        for i in 0..n_buckets {
            let base = 40 + i * 24; // magic(8) + width, n_buckets, n_reads, rewind (4 x u64)
            if base + 24 > head.len() {
                break;
            }
            table.push((u64_at(base), u64_at(base + 8), u64_at(base + 16)));
        }
        Ok(StoreManifest {
            width_bp,
            n_reads,
            compressed,
            table,
        })
    }
}

fn bad(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_string())
}

/// The postings key of a walk node: `|node|`, matching the panel-positions
/// table (a vertex's position is the same in either orientation).
pub fn node_key(node: i32) -> u64 {
    node.unsigned_abs() as u64
}

/// The walk nodes of one streamed packet, folded to `node_key` and deduplicated
/// (a read that revisits a node is ONE posting for it).
///
/// Packet = u32 payload_len, then payload = u32 seq_len, seq, u32 n_records,
/// then per record (u32 n_tokens, tokens…, u32 n_walk, (i32 node, u32 pos)…).
fn walk_node_keys(packet: &[u8], out: &mut Vec<u64>) -> io::Result<()> {
    out.clear();
    if packet.len() < 8 {
        return Err(bad("short packet"));
    }
    let payload_len = u32::from_le_bytes(packet[..4].try_into().unwrap()) as usize;
    let bytes = packet
        .get(4..4 + payload_len)
        .ok_or_else(|| bad("truncated packet"))?;
    let mut off = 0usize;
    let take_u32 = |off: &mut usize| -> io::Result<u32> {
        let end = *off + 4;
        if end > bytes.len() {
            return Err(bad("truncated field"));
        }
        let v = u32::from_le_bytes(bytes[*off..end].try_into().unwrap());
        *off = end;
        Ok(v)
    };
    let seq_len = take_u32(&mut off)? as usize;
    off += seq_len;
    if off > bytes.len() {
        return Err(bad("truncated sequence"));
    }
    let n_records = take_u32(&mut off)? as usize;
    for _ in 0..n_records {
        let n_tokens = take_u32(&mut off)? as usize;
        off += 8 * n_tokens;
        let n_walk = take_u32(&mut off)? as usize;
        for _ in 0..n_walk {
            let end = off + 8;
            if end > bytes.len() {
                return Err(bad("truncated walk"));
            }
            let node = i32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
            out.push(node_key(node));
            off = end;
        }
    }
    out.sort_unstable();
    out.dedup();
    Ok(())
}

/// A bitmap over the store's read ordinals — the cut's candidate-read set.
pub struct ReadBitmap {
    words: Vec<u64>,
    ones: u64,
}

impl ReadBitmap {
    pub fn new() -> ReadBitmap {
        ReadBitmap {
            words: Vec::new(),
            ones: 0,
        }
    }

    pub fn set(&mut self, read: u32) {
        let word = read as usize / 64;
        let mask = 1u64 << (read % 64);
        if word >= self.words.len() {
            self.words.resize(word + 1, 0);
        }
        if self.words[word] & mask == 0 {
            self.words[word] |= mask;
            self.ones += 1;
        }
    }

    pub fn contains(&self, read: u32) -> bool {
        let word = read as usize / 64;
        word < self.words.len() && self.words[word] & (1u64 << (read % 64)) != 0
    }

    /// The number of distinct reads marked.
    pub fn ones(&self) -> u64 {
        self.ones
    }

    /// The marked reads, ascending.
    pub fn iter(&self) -> impl Iterator<Item = u32> + '_ {
        self.words.iter().enumerate().flat_map(|(index, &word)| {
            let mut word = word;
            std::iter::from_fn(move || {
                if word == 0 {
                    return None;
                }
                let bit = word.trailing_zeros();
                word &= word - 1;
                Some((index * 64 + bit as usize) as u32)
            })
        })
    }
}

impl Default for ReadBitmap {
    fn default() -> ReadBitmap {
        ReadBitmap::new()
    }
}

/// An open `<prefix>.n2r`, mmapped (the directory is binary-searched per node,
/// so the artifact is never loaded whole — the same reason the panel-positions
/// table is mmapped).
pub struct NodePostings {
    map: memmap2::Mmap,
    node_count: usize,
    postings_count: u64,
    blob_base: usize,
}

impl NodePostings {
    pub fn open(path: &Path) -> io::Result<NodePostings> {
        let file = File::open(path)?;
        let map = unsafe { memmap2::Mmap::map(&file)? };
        if map.len() < HEADER_BYTES || &map[..4] != MAGIC {
            return Err(bad(&format!(
                "{}: not a node-postings artifact",
                path.display()
            )));
        }
        let version = u32::from_le_bytes(map[4..8].try_into().unwrap());
        if version != VERSION {
            return Err(bad(&format!(
                "{}: node-postings version {version}, this build reads {VERSION}",
                path.display()
            )));
        }
        let node_count = u64::from_le_bytes(map[8..16].try_into().unwrap()) as usize;
        let postings_count = u64::from_le_bytes(map[16..24].try_into().unwrap());
        let blob_base = HEADER_BYTES + node_count.saturating_mul(DIRECTORY_ENTRY_BYTES);
        if blob_base > map.len() {
            return Err(bad(&format!("{}: truncated directory", path.display())));
        }
        Ok(NodePostings {
            map,
            node_count,
            postings_count,
            blob_base,
        })
    }

    pub fn node_count(&self) -> usize {
        self.node_count
    }

    pub fn postings_count(&self) -> u64 {
        self.postings_count
    }

    fn entry(&self, index: usize) -> (u64, u64, u32) {
        let at = HEADER_BYTES + index * DIRECTORY_ENTRY_BYTES;
        let u64_at = |at: usize| u64::from_le_bytes(self.map[at..at + 8].try_into().unwrap());
        (
            u64_at(at),
            u64_at(at + 8),
            u32::from_le_bytes(self.map[at + 16..at + 20].try_into().unwrap()),
        )
    }

    /// The read ordinals touching `node`, ascending; empty for a node the
    /// store never saw (or one the artifact does not cover).
    pub fn postings(&self, node: i32) -> Vec<u32> {
        let key = node_key(node);
        let (mut lo, mut hi) = (0usize, self.node_count);
        while lo < hi {
            let mid = (lo + hi) / 2;
            match self.entry(mid).0.cmp(&key) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    let (_, offset, len) = self.entry(mid);
                    let start = self.blob_base + offset as usize;
                    let bytes = self.map.get(start..start + len as usize).unwrap_or(&[]);
                    return decode_postings(bytes);
                }
            }
        }
        Vec::new()
    }

    /// Mark the postings of a whole node set into a read bitmap. Returns the
    /// marked incidence count (sum of the postings lengths, read-to-node
    /// incidences, which is not the same number as the distinct reads).
    pub fn mark(&self, nodes: impl IntoIterator<Item = i32>) -> (u64, ReadBitmap) {
        let mut bitmap = ReadBitmap::new();
        let mut incidences = 0u64;
        for node in nodes {
            for read in self.postings(node) {
                bitmap.set(read);
                incidences += 1;
            }
        }
        (incidences, bitmap)
    }
}

fn encode_postings(out: &mut Vec<u8>, reads: &[u32]) {
    let mut previous = 0u64;
    for (index, &read) in reads.iter().enumerate() {
        let value = read as u64;
        write_varint(out, if index == 0 { value } else { value - previous });
        previous = value;
    }
}

fn decode_postings(bytes: &[u8]) -> Vec<u32> {
    let mut reads = Vec::new();
    let mut at = 0usize;
    let mut previous = 0u64;
    while let Some(delta) = read_varint(bytes, &mut at) {
        previous += delta;
        reads.push(previous as u32);
    }
    reads
}

fn write_varint(out: &mut Vec<u8>, mut value: u64) {
    loop {
        let low = (value & 0x7f) as u8;
        value >>= 7;
        if value == 0 {
            out.push(low);
            return;
        }
        out.push(low | 0x80);
    }
}

fn read_varint(bytes: &[u8], at: &mut usize) -> Option<u64> {
    let mut value = 0u64;
    let mut shift = 0u32;
    loop {
        let byte = *bytes.get(*at)?;
        *at += 1;
        value |= ((byte & 0x7f) as u64) << shift;
        if byte & 0x80 == 0 {
            return Some(value);
        }
        shift += 7;
        if shift > 63 {
            return None;
        }
    }
}

pub struct BuildSummary {
    pub node_count: u64,
    pub postings_count: u64,
    pub artifact_bytes: u64,
    pub buckets: usize,
    pub records: u64,
    pub reads_with_walks: u64,
    pub incidences: u64,
    pub store_bytes_read: u64,
    pub store_bytes_decoded: u64,
}

/// Build `<prefix>.n2r` from `<prefix>.evd` + `<prefix>.evd.idx`.
///
/// The store is opened READ-ONLY and nothing about its format changes: this is
/// one sequential pass over every bucket, decompressing it and reading each
/// record's walk nodes. Two passes over the incidences (spill, then sort +
/// merge each partition) keep peak memory at one partition rather than the
/// whole sample's postings.
pub fn build(prefix: &str) -> io::Result<BuildSummary> {
    let t0 = Instant::now();
    let manifest = StoreManifest::read(prefix)?;
    let spill_dir = PathBuf::from(format!("{prefix}.n2r-spill"));
    std::fs::create_dir_all(&spill_dir)?;

    let mut spills: Vec<BufWriter<File>> = (0..SPILLS)
        .map(|i| File::create(spill_dir.join(format!("{i:03}.spill"))).map(BufWriter::new))
        .collect::<io::Result<_>>()?;
    let mut payload = File::open(format!("{prefix}.evd"))?;
    let mut raw: Vec<u8> = Vec::new();
    let mut keys: Vec<u64> = Vec::new();
    let (mut records, mut reads_with_walks, mut incidences) = (0u64, 0u64, 0u64);
    let mut store_bytes_read = 0u64;
    let mut store_bytes_decoded = 0u64;
    for &(_, offset, len) in &manifest.table {
        payload.seek(SeekFrom::Start(offset))?;
        raw.resize(len as usize, 0);
        payload.read_exact(&mut raw)?;
        store_bytes_read += len;
        let decoded;
        let bytes: &[u8] = if manifest.compressed {
            decoded = zstd::stream::decode_all(std::io::Cursor::new(&raw[..]))?;
            &decoded
        } else {
            &raw[..]
        };
        store_bytes_decoded += bytes.len() as u64;
        // a bucket's records: u32 read ordinal, then the packet (its own u32 length first)
        let mut at = 0usize;
        while at + 8 <= bytes.len() {
            let ordinal = u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap());
            let plen = u32::from_le_bytes(bytes[at + 4..at + 8].try_into().unwrap()) as usize;
            let end = at + 8 + plen;
            if end > bytes.len() {
                break;
            }
            walk_node_keys(&bytes[at + 4..end], &mut keys)?;
            if !keys.is_empty() {
                reads_with_walks += 1;
            }
            incidences += keys.len() as u64;
            for &key in &keys {
                let spill = &mut spills[(key % SPILLS) as usize];
                spill.write_all(&key.to_le_bytes())?;
                spill.write_all(&ordinal.to_le_bytes())?;
            }
            records += 1;
            at = end;
        }
    }
    drop(spills);

    // pass 2: sort each partition, then drop it (the sorted run replaces it)
    let mut runs: Vec<Option<BufReader<File>>> = Vec::with_capacity(SPILLS as usize);
    for i in 0..SPILLS {
        let spill = spill_dir.join(format!("{i:03}.spill"));
        let mut bytes = Vec::new();
        File::open(&spill)?.read_to_end(&mut bytes)?;
        let mut pairs: Vec<(u64, u32)> = Vec::with_capacity(bytes.len() / SPILL_RECORD_BYTES);
        let mut at = 0usize;
        while at + SPILL_RECORD_BYTES <= bytes.len() {
            pairs.push((
                u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap()),
                u32::from_le_bytes(bytes[at + 8..at + 12].try_into().unwrap()),
            ));
            at += SPILL_RECORD_BYTES;
        }
        // a read whose walk spans several walk nodes is stored in EVERY one's
        // bucket, so the same (node, read) incidence arrives once per placement
        pairs.sort_unstable();
        pairs.dedup();
        let run = spill_dir.join(format!("{i:03}.run"));
        {
            let mut out = BufWriter::new(File::create(&run)?);
            for (key, read) in &pairs {
                out.write_all(&key.to_le_bytes())?;
                out.write_all(&read.to_le_bytes())?;
            }
            out.flush()?;
        }
        drop(pairs);
        std::fs::remove_file(&spill)?;
        runs.push(Some(BufReader::with_capacity(1 << 20, File::open(&run)?)));
    }

    // merge the sorted runs: ascending node, ascending read within a node
    let dir_path = spill_dir.join("directory.tmp");
    let blob_path = spill_dir.join("postings.tmp");
    let mut directory = BufWriter::new(File::create(&dir_path)?);
    let mut blob = BufWriter::new(File::create(&blob_path)?);
    struct Run {
        reader: BufReader<File>,
        head: Option<(u64, u32)>,
    }
    fn advance(run: &mut Run) -> io::Result<()> {
        let mut record = [0u8; SPILL_RECORD_BYTES];
        match run.reader.read_exact(&mut record) {
            Ok(()) => {
                run.head = Some((
                    u64::from_le_bytes(record[..8].try_into().unwrap()),
                    u32::from_le_bytes(record[8..12].try_into().unwrap()),
                ));
                Ok(())
            }
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => {
                run.head = None;
                Ok(())
            }
            Err(e) => Err(e),
        }
    }
    let mut readers: Vec<Run> = runs
        .into_iter()
        .flatten()
        .map(|reader| Run { reader, head: None })
        .collect();
    let mut heap: BinaryHeap<Reverse<(u64, u32, usize)>> = BinaryHeap::new();
    for (index, run) in readers.iter_mut().enumerate() {
        advance(run)?;
        if let Some((key, read)) = run.head {
            heap.push(Reverse((key, read, index)));
        }
    }
    let (mut node_count, mut postings_count) = (0u64, 0u64);
    let mut blob_offset = 0u64;
    let mut encoded: Vec<u8> = Vec::new();
    let mut group: Vec<u32> = Vec::new();
    let mut group_key: Option<u64> = None;
    while let Some(Reverse((key, read, index))) = heap.pop() {
        if group_key != Some(key) {
            if let Some(previous) = group_key {
                encoded.clear();
                encode_postings(&mut encoded, &group);
                blob.write_all(&encoded)?;
                directory.write_all(&previous.to_le_bytes())?;
                directory.write_all(&blob_offset.to_le_bytes())?;
                directory.write_all(&(encoded.len() as u32).to_le_bytes())?;
                blob_offset += encoded.len() as u64;
                postings_count += group.len() as u64;
                node_count += 1;
                group.clear();
            }
            group_key = Some(key);
        }
        group.push(read);
        advance(&mut readers[index])?;
        if let Some((next_key, next_read)) = readers[index].head {
            heap.push(Reverse((next_key, next_read, index)));
        }
    }
    if let Some(previous) = group_key {
        encoded.clear();
        encode_postings(&mut encoded, &group);
        blob.write_all(&encoded)?;
        directory.write_all(&previous.to_le_bytes())?;
        directory.write_all(&blob_offset.to_le_bytes())?;
        directory.write_all(&(encoded.len() as u32).to_le_bytes())?;
        postings_count += group.len() as u64;
        node_count += 1;
    }
    directory.flush()?;
    blob.flush()?;
    drop(directory);
    drop(blob);
    drop(readers);

    // header, directory, blob — one artifact, written once
    let artifact = format!("{prefix}.n2r");
    {
        let mut out = BufWriter::new(File::create(&artifact)?);
        out.write_all(MAGIC)?;
        out.write_all(&VERSION.to_le_bytes())?;
        out.write_all(&node_count.to_le_bytes())?;
        out.write_all(&postings_count.to_le_bytes())?;
        io::copy(&mut File::open(&dir_path)?, &mut out)?;
        io::copy(&mut File::open(&blob_path)?, &mut out)?;
        out.flush()?;
    }
    std::fs::remove_dir_all(&spill_dir)?;

    let artifact_bytes = std::fs::metadata(&artifact)?.len();
    let _ = t0.elapsed();
    Ok(BuildSummary {
        node_count,
        postings_count,
        artifact_bytes,
        buckets: manifest.table.len(),
        records,
        reads_with_walks,
        incidences,
        store_bytes_read,
        store_bytes_decoded,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A synthetic evidence store: one bucket holding the given
    /// (read ordinal, walk nodes) records.
    fn write_store(dir: &Path, compressed: bool, records: &[(u32, &[i32])]) -> String {
        let prefix = dir.join("sample");
        let prefix = prefix.to_str().unwrap().to_string();
        let mut bucket = Vec::new();
        for &(ordinal, nodes) in records {
            let mut payload = Vec::new();
            payload.extend_from_slice(&4u32.to_le_bytes()); // seq_len
            payload.extend_from_slice(b"ACGT");
            payload.extend_from_slice(&1u32.to_le_bytes()); // n_records
            payload.extend_from_slice(&1u32.to_le_bytes()); // n_tokens
            payload.extend_from_slice(&7u64.to_le_bytes()); // one token
            payload.extend_from_slice(&(nodes.len() as u32).to_le_bytes());
            for &node in nodes {
                payload.extend_from_slice(&node.to_le_bytes());
                payload.extend_from_slice(&0u32.to_le_bytes());
            }
            bucket.extend_from_slice(&ordinal.to_le_bytes());
            bucket.extend_from_slice(&(payload.len() as u32).to_le_bytes());
            bucket.extend_from_slice(&payload);
        }
        let stored = if compressed {
            zstd::stream::encode_all(&bucket[..], 3).unwrap()
        } else {
            bucket
        };
        let mut evd = File::create(format!("{prefix}.evd")).unwrap();
        evd.write_all(&stored).unwrap();
        evd.flush().unwrap();
        let mut idx = File::create(format!("{prefix}.evd.idx")).unwrap();
        idx.write_all(if compressed {
            EVD_IDX_ZSTD
        } else {
            EVD_IDX_RAW
        })
        .unwrap();
        idx.write_all(&10_000u64.to_le_bytes()).unwrap(); // width_bp
        idx.write_all(&1u64.to_le_bytes()).unwrap(); // one bucket
        idx.write_all(&5u64.to_le_bytes()).unwrap(); // n_reads
        idx.write_all(&0u64.to_le_bytes()).unwrap(); // rewind
        idx.write_all(&(7u64).to_le_bytes()).unwrap(); // bucket key
        idx.write_all(&0u64.to_le_bytes()).unwrap(); // offset
        idx.write_all(&(stored.len() as u64).to_le_bytes()).unwrap();
        idx.flush().unwrap();
        prefix
    }

    /// 3 nodes, 5 reads: A={1,3}, B={3,5}, C={2,4}.
    const NODE_A: i32 = 100;
    const NODE_B: i32 = 200;
    const NODE_C: i32 = 300;
    const RECORDS: [(u32, &[i32]); 5] = [
        (1, &[NODE_A]),
        (2, &[NODE_C]),
        (3, &[NODE_A, NODE_B]),
        (4, &[NODE_C]),
        (5, &[NODE_B]),
    ];

    fn check(compressed: bool) {
        let dir = tempfile::tempdir().unwrap();
        let prefix = write_store(dir.path(), compressed, &RECORDS);
        let summary = build(&prefix).unwrap();
        assert_eq!(summary.node_count, 3);
        assert_eq!(summary.postings_count, 6);
        assert_eq!(summary.records, 5);
        let postings = NodePostings::open(&PathBuf::from(format!("{prefix}.n2r"))).unwrap();
        assert_eq!(postings.postings(NODE_A), vec![1, 3]);
        assert_eq!(postings.postings(NODE_B), vec![3, 5]);
        assert_eq!(postings.postings(NODE_C), vec![2, 4]);
        // the strand fold: the negative id is the same vertex
        assert_eq!(postings.postings(-NODE_A), vec![1, 3]);
        assert!(postings.postings(999).is_empty());
        let (incidences, bitmap) = postings.mark([NODE_A, NODE_B]);
        assert_eq!(incidences, 4);
        assert_eq!(bitmap.iter().collect::<Vec<u32>>(), vec![1, 3, 5]);
        assert_eq!(bitmap.ones(), 3);
    }

    #[test]
    fn node_postings_decode_the_synthetic_sets_from_a_compressed_store() {
        check(true);
    }

    #[test]
    fn node_postings_decode_the_synthetic_sets_from_a_raw_store() {
        check(false);
    }

    #[test]
    fn varint_round_trips_delta_coded_postings() {
        let reads = vec![0u32, 1, 3, 1 << 20, (1 << 31) - 1];
        let mut encoded = Vec::new();
        encode_postings(&mut encoded, &reads);
        assert_eq!(decode_postings(&encoded), reads);
    }

    #[test]
    fn read_bitmap_marks_and_iterates_ascending() {
        let mut bitmap = ReadBitmap::new();
        for read in [65u32, 3, 130, 3] {
            bitmap.set(read);
        }
        assert_eq!(bitmap.ones(), 3);
        assert!(bitmap.contains(130) && !bitmap.contains(131));
        assert_eq!(bitmap.iter().collect::<Vec<u32>>(), vec![3, 65, 130]);
    }

    /// The cut's candidate supply: the filter refuses records before they are
    /// copied out, so the survivors are a strict subset of the unfiltered
    /// fetch — same packets, same ordinal dedup, fewer of them.
    #[test]
    fn filtered_fetch_is_a_strict_subset_of_the_unfiltered_fetch() {
        use crate::genome_inference::evidence_store::EvidenceStore;
        let dir = tempfile::tempdir().unwrap();
        let prefix = write_store(dir.path(), true, &RECORDS);
        let (incidences, candidates) = {
            build(&prefix).unwrap();
            let postings = NodePostings::open(&PathBuf::from(format!("{prefix}.n2r"))).unwrap();
            postings.mark([NODE_A, NODE_B])
        };
        assert_eq!(incidences, 4);
        let mut store = EvidenceStore::open(&prefix).unwrap();
        let (unfiltered, all) = store.packets_for_buckets_filtered(&[7], |_| true).unwrap();
        assert_eq!(all.records, 5);
        assert_eq!(all.kept, 5);
        assert_eq!(all.skipped, 0);
        let (filtered, stats) = store
            .packets_for_buckets_filtered(&[7], |read| candidates.contains(read))
            .unwrap();
        assert_eq!(stats.records, 5);
        assert_eq!(stats.kept, 3);
        assert_eq!(stats.skipped, 2);
        assert_eq!(filtered.len(), 3);
        // 5 unfiltered + 3 filtered
        assert_eq!(store.fetched_packets, 8);
        // the two refused reads are the node-C ones (2, 4): their packets are
        // byte-identical to each other and absent from the survivors
        let node_c = unfiltered[1].clone();
        assert_eq!(unfiltered[3], node_c);
        assert!(!filtered.contains(&node_c));
        for packet in &filtered {
            assert!(unfiltered.contains(packet));
        }
    }
}
