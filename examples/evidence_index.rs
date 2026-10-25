//! THE EVIDENCE STORE — read evidence addressable by panel position.
//!
//! The projection emits reads in ARRIVAL order (`shard.<r>.bin.gz`: a gzip
//! stream of u32-length-prefixed packets, each a read with its MEM records
//! = token walks whose nodes ARE the panel coordinates). Because arrival
//! order is not position order, nothing can be seeked: the restricted cut
//! decompresses a whole region shard per territory (~46 GB per territory
//! across all shards) to use ~0.033% of the reads.
//!
//! This is the missing layer, with BAM's own trick: make POSITION ORDER
//! == FILE ORDER, then a sparse manifest is all the index you need.
//!
//! ```text
//! evidence_index build <dir> <shard-prefix>       # bucket the packets by panel position
//! evidence_index evidence <dir> <prefix> <lo> <hi> # the reads with walks in [lo,hi]
//! evidence_index summary  <dir> <prefix> <lo> <hi> # counts/mass, no packets touched
//! evidence_index scan     <dir> <prefix> <lo> <hi> # the full scan, for the identity gate
//! ```
//!
//! BUCKET == LOCUS: the bucket is BUCKET_NODES of panel position (~10 kb
//! at one node per base), so a locus query touches one or two buckets —
//! a seek, not a scan. `IMPG_EVIDENCE_BUCKET_NODES` overrides.
//!
//! WHAT IT DOES NOT DO YET: the restricted-basis filter (which reads
//! count for a territory) is a separate artifact (`basis(territory)`),
//! because it is a per-sample set, not a per-position property.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

const MAGIC: &[u8; 8] = b"EVDIDX01";
/// THE COMPRESSED STORE (per-bucket zstd). The store is 168 GB raw -
/// 2/3 of it the reads' own sequence - and compresses 7.4x at zstd-3
/// (measured at 0/40/80% of the file: 7.48/7.39/7.42), so the artifact
/// lands near 23 GB instead of 168. A bucket is the seek unit and the
/// median bucket is 2.3 KB, so the decompression is microseconds.
const MAGIC_Z: &[u8; 8] = b"EVDIDX03";
const POS_MAGIC: &[u8; 8] = b"PNLPOS01";

/// `node -> (path, offset)`, sorted by node. Panel-side, built by
/// `panel_positions`. This is the key that gives a walk LOCALITY: its 3
/// nodes sit at ~one position, so the read lands in ONE bucket, and the
/// bucket is a locus (the beds' own `path,start,end` coordinate).
struct Positions {
    entries: Vec<(i32, u32, u64)>,
}

impl Positions {
    fn load(path: &str) -> io::Result<Positions> {
        let mut bytes = Vec::new();
        File::open(path)?.read_to_end(&mut bytes)?;
        if bytes.len() < 32 || &bytes[..8] != POS_MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "bad positions table"));
        }
        let u64_at = |at: usize| u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
        let n = u64_at(8) as usize;
        let mut entries = Vec::with_capacity(n);
        let mut at = 40usize; // magic(8) + n,path_count,end_path,total_bp (4x8)
        for _ in 0..n {
            let node = i32::from_le_bytes(bytes[at..at + 4].try_into().unwrap());
            let path = u32::from_le_bytes(bytes[at + 4..at + 8].try_into().unwrap());
            let offset = u64::from_le_bytes(bytes[at + 8..at + 16].try_into().unwrap());
            entries.push((node, path, offset));
            at += 16;
        }
        entries.sort_unstable();
        Ok(Positions { entries })
    }

    /// The panel position of a node, or None when the panel has no walk
    /// for it (a node the table's path range does not cover).
    fn position(&self, node: i32) -> Option<(u32, u64)> {
        // THE KEY IS |node|: a vertex's position is the same in either
        // orientation, and the panel walk emits most vertices only once
        // (forward). Keying on the signed id silently dropped EVERY
        // reverse-strand read - measured: node 118354468 FOUND,
        // -118354468 ABSENT, so 20% of walk-bearing reads never placed.
        let node = node.unsigned_abs() as i32;
        match self.entries.binary_search_by_key(&node, |e| e.0) {
            Ok(at) => Some((self.entries[at].1, self.entries[at].2)),
            Err(_) => None,
        }
    }
}

/// The store's key: (path, offset/width) packed into a u64. Paths number
/// < 2^16 and offsets/width < 2^24 at any sane width, so this is exact.
fn position_key(path: u32, offset: u64, width: u64) -> u64 {
    ((path as u64) << 24) | (offset / width)
}


fn bucket_nodes() -> u64 {
    std::env::var("IMPG_EVIDENCE_BUCKET_NODES")
        .ok()
        .and_then(|v| v.trim().parse::<u64>().ok())
        .filter(|v| *v >= 256)
        .unwrap_or(10_000)
}

/// A panel position on EITHER strand lands in one bucket: the reverse
/// strand's negative node id is folded by magnitude, so a locus query
/// returns both orientations (which is what a locus query means).
fn bucket_of(node: i32, width: u64) -> u64 {
    (node.unsigned_abs() as u64) / width
}

/// One packet = u32 payload_len + payload; payload = u32 seq_len, seq,
/// u32 n_records, then per record (u32 n_tokens, tokens…) and
/// (u32 n_walk, (i32 node, u32 pos)…). The walk nodes are the panel
/// coordinates this store indexes.
fn walk_nodes(packet: &[u8]) -> io::Result<(Vec<i32>, u64, u64)> {
    let bad = |m: &str| io::Error::new(io::ErrorKind::InvalidData, m.to_string());
    if packet.len() < 8 {
        return Err(bad("short packet"));
    }
    let payload_len = u32::from_le_bytes(packet[..4].try_into().unwrap()) as usize;
    let bytes = packet
        .get(4..4 + payload_len)
        .ok_or_else(|| bad("truncated packet"))?;
    let mut off = 0usize;
    let mut u32_at = |off: &mut usize| -> io::Result<u32> {
        let end = *off + 4;
        if end > bytes.len() {
            return Err(bad("truncated field"));
        }
        let v = u32::from_le_bytes(bytes[*off..end].try_into().unwrap());
        *off = end;
        Ok(v)
    };
    let seq_len = u32_at(&mut off)? as usize;
    off += seq_len;
    if off > bytes.len() {
        return Err(bad("truncated sequence"));
    }
    let n_records = u32_at(&mut off)? as usize;
    let mut nodes = Vec::new();
    let mut tokens_total = 0u64;
    for _ in 0..n_records {
        let n_tokens = u32_at(&mut off)? as usize;
        tokens_total += n_tokens as u64;
        off += 8 * n_tokens;
        let n_walk = u32_at(&mut off)? as usize;
        for _ in 0..n_walk {
            let end = off + 8;
            if end > bytes.len() {
                return Err(bad("truncated walk"));
            }
            nodes.push(i32::from_le_bytes(bytes[off..off + 4].try_into().unwrap()));
            off = end;
        }
    }
    Ok((nodes, n_records as u64, tokens_total))
}

fn shard_gz(dir: &str, prefix: &str) -> PathBuf {
    Path::new(dir).join(format!("shard.{prefix}.bin.gz"))
}

fn read_packets(path: &Path) -> io::Result<Vec<Vec<u8>>> {
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

fn build_keyed(dir: &str, prefix: &str, positions_path: &str, width_bp: u64) -> io::Result<()> {
    let t0 = Instant::now();
    let positions = Positions::load(positions_path)?;
    let load_positions = t0.elapsed();
    let packets = read_packets(&shard_gz(dir, prefix))?;
    let read_time = t0.elapsed() - load_positions;

    // place each read in the bucket(s) of its walk nodes' PANEL POSITIONS
    let t1 = Instant::now();
    let mut buckets: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    let mut unplaced = 0usize;
    let mut multi = 0usize;
    for (index, packet) in packets.iter().enumerate() {
        let (nodes, _, _) = walk_nodes(packet)?;
        let mut keys: Vec<u64> = Vec::new();
        for &node in &nodes {
            if let Some((path, offset)) = positions.position(node) {
                keys.push(position_key(path, offset, width_bp));
            }
        }
        if keys.is_empty() {
            unplaced += 1;
            continue;
        }
        keys.sort_unstable();
        keys.dedup();
        if keys.len() > 1 {
            multi += 1;
        }
        for key in keys {
            buckets.entry(key).or_default().push(index);
        }
    }
    let place = t1.elapsed();

    let t2 = Instant::now();
    let evd = Path::new(dir).join(format!("shard.{prefix}.p.evd"));
    let idx = Path::new(dir).join(format!("shard.{prefix}.p.evd.idx"));
    let sum = Path::new(dir).join(format!("shard.{prefix}.p.evd.sum"));
    let mut table: Vec<(u64, u64, u64)> = Vec::with_capacity(buckets.len());
    {
        let mut out = BufWriter::new(File::create(&evd)?);
        let mut offset = 0u64;
        for (&bucket, members) in &buckets {
            let start = offset;
            for &index in members {
                // THE READ ORDINAL travels with the packet: a read whose
                // walk spans two buckets is stored in BOTH (the price of
                // locality), so a multi-bucket query must be able to
                // dedup - without an identity it double-counts (caught by
                // the gate: 19,270 vs the scan's 17,192).
                out.write_all(&(index as u32).to_le_bytes())?;
                out.write_all(&packets[index])?;
                offset += 4 + packets[index].len() as u64;
            }
            table.push((bucket, start, offset - start));
        }
        out.flush()?;
    }
    {
        let mut out = BufWriter::new(File::create(&idx)?);
        out.write_all(MAGIC)?;
        out.write_all(&width_bp.to_le_bytes())?;
        out.write_all(&(table.len() as u64).to_le_bytes())?;
        out.write_all(&(packets.len() as u64).to_le_bytes())?;
        out.write_all(&0u64.to_le_bytes())?; // rewind is 0: the key has locality
        for &(bucket, offset, len) in &table {
            out.write_all(&bucket.to_le_bytes())?;
            out.write_all(&offset.to_le_bytes())?;
            out.write_all(&len.to_le_bytes())?;
        }
    }
    {
        let mut out = BufWriter::new(File::create(&sum)?);
        for (&bucket, members) in &buckets {
            let (mut reads, mut records, mut tokens) = (0u64, 0u64, 0u64);
            for &index in members {
                let (nodes, recs, toks) = walk_nodes(&packets[index])?;
                if !nodes.is_empty() {
                    reads += 1;
                    records += recs;
                    tokens += toks;
                }
            }
            out.write_all(&bucket.to_le_bytes())?;
            out.write_all(&reads.to_le_bytes())?;
            out.write_all(&records.to_le_bytes())?;
            out.write_all(&tokens.to_le_bytes())?;
        }
    }
    let write = t2.elapsed();

    let evd_bytes = std::fs::metadata(&evd)?.len();
    let idx_bytes = std::fs::metadata(&idx)?.len();
    let sum_bytes = std::fs::metadata(&sum)?.len();
    println!("shard.{prefix} (POSITION-keyed, {width_bp} bp/bucket)");
    println!("  positions table    {} nodes ({:.1}s)", positions.entries.len(), load_positions.as_secs_f64());
    println!("  reads              {}", packets.len());
    println!("  placed             {} (unplaced {unplaced}, spanning >1 bucket {multi})", packets.len() - unplaced);
    println!("  buckets occupied   {}", table.len());
    println!("  read {:.1}s  place {:.1}s  write {:.1}s", read_time.as_secs_f64(), place.as_secs_f64(), write.as_secs_f64());
    println!(
        "  store {:.1} MB   index {:.2} MB   summary {:.2} MB   (overhead {:.3}%)",
        evd_bytes as f64 / 1e6,
        idx_bytes as f64 / 1e6,
        sum_bytes as f64 / 1e6,
        100.0 * (idx_bytes + sum_bytes) as f64 / evd_bytes.max(1) as f64
    );
    Ok(())
}

fn build(dir: &str, prefix: &str) -> io::Result<()> {
    let width = bucket_nodes();
    let t0 = Instant::now();
    let packets = read_packets(&shard_gz(dir, prefix))?;
    let read_time = t0.elapsed();

    let t1 = Instant::now();
    let mut buckets: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    let mut max_span_buckets = 0u64;
    for (index, packet) in packets.iter().enumerate() {
        let (nodes, _, _) = walk_nodes(packet)?;
        if nodes.is_empty() {
            continue;
        }
        let mut lo = u64::MAX;
        let mut hi = 0u64;
        for &node in &nodes {
            let b = bucket_of(node, width);
            lo = lo.min(b);
            hi = hi.max(b);
        }
        // keep the read in its LOWEST bucket; the query rewinds by the
        // widest span so a long walk is never missed
        buckets.entry(lo).or_default().push(index);
        max_span_buckets = max_span_buckets.max(hi - lo);
    }
    let place = t1.elapsed();

    // ---- write: packets verbatim, in bucket order ----
    let t2 = Instant::now();
    let evd = Path::new(dir).join(format!("shard.{prefix}.evd"));
    let idx = Path::new(dir).join(format!("shard.{prefix}.evd.idx"));
    let sum = Path::new(dir).join(format!("shard.{prefix}.evd.sum"));
    let mut table: Vec<(u64, u64, u64)> = Vec::with_capacity(buckets.len());
    {
        let mut out = BufWriter::new(File::create(&evd)?);
        let mut offset = 0u64;
        for (&bucket, members) in &buckets {
            let start = offset;
            for &index in members {
                // THE READ ORDINAL travels with the packet: a read whose
                // walk spans two buckets is stored in BOTH (the price of
                // locality), so a multi-bucket query must be able to
                // dedup - without an identity it double-counts (caught by
                // the gate: 19,270 vs the scan's 17,192).
                out.write_all(&(index as u32).to_le_bytes())?;
                out.write_all(&packets[index])?;
                offset += 4 + packets[index].len() as u64;
            }
            table.push((bucket, start, offset - start));
        }
        out.flush()?;
    }
    {
        let mut out = BufWriter::new(File::create(&idx)?);
        out.write_all(MAGIC)?;
        out.write_all(&width.to_le_bytes())?;
        out.write_all(&(table.len() as u64).to_le_bytes())?;
        out.write_all(&(packets.len() as u64).to_le_bytes())?;
        out.write_all(&max_span_buckets.to_le_bytes())?;
        for &(bucket, offset, len) in &table {
            out.write_all(&bucket.to_le_bytes())?;
            out.write_all(&offset.to_le_bytes())?;
            out.write_all(&len.to_le_bytes())?;
        }
    }
    {
        // the O(1) layer: per bucket, reads / records / tokens — enough to
        // answer summary() and plan() without touching a packet
        let mut out = BufWriter::new(File::create(&sum)?);
        for (&bucket, members) in &buckets {
            let (mut reads, mut records, mut tokens) = (0u64, 0u64, 0u64);
            for &index in members {
                let (nodes, recs, toks) = walk_nodes(&packets[index])?;
                if !nodes.is_empty() {
                    reads += 1;
                    records += recs;
                    tokens += toks;
                }
            }
            out.write_all(&bucket.to_le_bytes())?;
            out.write_all(&reads.to_le_bytes())?;
            out.write_all(&records.to_le_bytes())?;
            out.write_all(&tokens.to_le_bytes())?;
        }
    }
    let write = t2.elapsed();

    let evd_bytes = std::fs::metadata(&evd)?.len();
    let idx_bytes = std::fs::metadata(&idx)?.len();
    let sum_bytes = std::fs::metadata(&sum)?.len();
    println!("shard.{prefix}");
    println!("  bucket width       {width} nodes (~{} kb)", width / 1000);
    println!("  reads (with walks) {}", table.iter().map(|_| 0).count() * 0 + buckets.values().map(|v| v.len()).sum::<usize>());
    println!("  buckets            {}", table.len());
    println!("  max walk span      {max_span_buckets} buckets (the query's rewind)");
    println!("  read {:.1}s  place {:.1}s  write {:.1}s", read_time.as_secs_f64(), place.as_secs_f64(), write.as_secs_f64());
    println!(
        "  evd {:.1} MB   idx {:.2} MB   sum {:.2} MB   (index overhead {:.2}% of the store)",
        evd_bytes as f64 / 1e6,
        idx_bytes as f64 / 1e6,
        sum_bytes as f64 / 1e6,
        100.0 * (idx_bytes + sum_bytes) as f64 / evd_bytes.max(1) as f64
    );
    Ok(())
}

/// Whether this store's buckets are zstd frames (format EVDIDX03).
/// The store's file stem: the whole-sample builder writes `<prefix>.evd*`,
/// the per-shard prototype wrote `shard.<prefix>.p.evd*`.
fn store_stem(dir: &str, prefix: &str, ext: &str) -> PathBuf {
    let whole = Path::new(dir).join(format!("{prefix}.{ext}"));
    if whole.exists() {
        return whole;
    }
    Path::new(dir).join(format!("shard.{prefix}.p.{ext}"))
}

fn manifest_is_compressed(dir: &str, prefix: &str) -> io::Result<bool> {
    let mut head = [0u8; 8];
    File::open(store_stem(dir, prefix, "evd.idx"))?.read_exact(&mut head)?;
    Ok(&head == MAGIC_Z)
}

fn load_manifest(dir: &str, prefix: &str) -> io::Result<(u64, Vec<(u64, u64, u64)>, u64)> {
    let mut buf = Vec::new();
    File::open(Path::new(dir).join(format!("shard.{prefix}.evd.idx")))?.read_to_end(&mut buf)?;
    if buf.len() < 32 || (&buf[..8] != MAGIC && &buf[..8] != MAGIC_Z) {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "bad evd.idx"));
    }
    let u64_at = |at: usize| u64::from_le_bytes(buf[at..at + 8].try_into().unwrap());
    let width = u64_at(8);
    let n_buckets = u64_at(16) as usize;
    let rewind = u64_at(32);
    let mut table = Vec::with_capacity(n_buckets);
    for i in 0..n_buckets {
        let base = 40 + i * 24; // magic(8) + width,n_buckets,n_reads,rewind (4x8)
        table.push((u64_at(base), u64_at(base + 8), u64_at(base + 16)));
    }
    Ok((width, table, rewind))
}

/// The seek: read only the buckets covering [lo_bucket, hi_bucket].
fn read_buckets(
    dir: &str,
    prefix: &str,
    lo: u64,
    hi: u64,
) -> io::Result<(Vec<Vec<u8>>, usize, u64)> {
    let (_, table, rewind) = load_manifest(dir, prefix)?;
    let is_compressed = manifest_is_compressed(dir, prefix)?;
    let first_bucket = lo.saturating_sub(rewind);
    let spine: Vec<u64> = table.iter().map(|&(b, _, _)| b).collect();
    let start = spine.partition_point(|&b| b < first_bucket);
    let mut file = File::open(Path::new(dir).join(format!("shard.{prefix}.evd")))?;
    let mut packets = Vec::new();
    let mut buckets_read = 0usize;
    let mut bytes_read = 0u64;
    for &(bucket, offset, len) in &table[start..] {
        if bucket > hi {
            break;
        }
        use std::io::{Seek, SeekFrom};
        file.seek(SeekFrom::Start(offset))?;
        let mut buf = vec![0u8; len as usize];
        file.read_exact(&mut buf)?;
        bytes_read += len;
        buckets_read += 1;
        // split the bucket's bytes back into packets
        let mut at = 0usize;
        while at + 4 <= buf.len() {
            let plen = u32::from_le_bytes(buf[at..at + 4].try_into().unwrap()) as usize;
            let end = at + 4 + plen;
            if end > buf.len() {
                break;
            }
            packets.push(buf[at..end].to_vec());
            at = end;
        }
    }
    Ok((packets, buckets_read, bytes_read))
}

/// `evidence(locus)` — THE interface: a (path, bp range) in, the reads with
/// evidence there out, by seek.
fn evidence_pos(dir: &str, prefix: &str, positions_path: &str, path: u32, lo_bp: u64, hi_bp: u64) -> io::Result<()> {
    let positions = Positions::load(positions_path)?;
    let mut head = Vec::new();
    File::open(store_stem(dir, prefix, "evd.idx"))?.read_to_end(&mut head)?;
    let u64_at = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().unwrap());
    let is_compressed = &head[..8] == MAGIC_Z;
    let width = u64_at(8);
    let n_buckets = u64_at(16) as usize;
    let mut table = Vec::with_capacity(n_buckets);
    for i in 0..n_buckets {
        let base = 40 + i * 24;
        table.push((u64_at(base), u64_at(base + 8), u64_at(base + 16)));
    }
    let lo_key = position_key(path, lo_bp, width);
    let hi_key = position_key(path, hi_bp, width);
    let spine: Vec<u64> = table.iter().map(|&(b, _, _)| b).collect();
    let start = spine.partition_point(|&b| b < lo_key);
    let t0 = Instant::now();
    let mut file = File::open(store_stem(dir, prefix, "evd"))?;
    use std::io::{Seek, SeekFrom};
    let mut hits = 0usize;
    let mut examined = 0usize;
    let mut bytes_read = 0u64;
    let mut bucket_count = 0usize;
    // a read spanning two buckets is stored twice; count it once
    let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
    for &(key, offset, len) in &table[start..] {
        if key > hi_key {
            break;
        }
        file.seek(SeekFrom::Start(offset))?;
        let mut buf = vec![0u8; len as usize];
        file.read_exact(&mut buf)?;
        bytes_read += len;
        bucket_count += 1;
        let stored = if is_compressed {
            // one zstd frame per bucket: decode it, then read the records
            zstd::stream::decode_all(std::io::Cursor::new(&buf))?
        } else {
            buf
        };
        let buf = stored;
        let mut at = 0usize;
        while at + 8 <= buf.len() {
            let ordinal = u32::from_le_bytes(buf[at..at + 4].try_into().unwrap());
            let plen = u32::from_le_bytes(buf[at + 4..at + 8].try_into().unwrap()) as usize;
            let end = at + 8 + plen;
            if end > buf.len() { break; }
            let packet = &buf[at + 4..end];
            if seen.insert(ordinal) {
                examined += 1;
                let (nodes, _, _) = walk_nodes(packet)?;
                if nodes.iter().any(|&n| matches!(positions.position(n), Some((p, o)) if p == path && o >= lo_bp && o <= hi_bp)) {
                    hits += 1;
                }
            }
            at = end;
        }
    }
    println!("evidence(path {path}, {lo_bp}..{hi_bp}): {examined} reads examined, {hits} with evidence, {:.3} MB read, {bucket_count} buckets, {:.2} ms",
        bytes_read as f64/1e6, t0.elapsed().as_secs_f64()*1e3);
    Ok(())
}

/// The gate's other arm: the same predicate over every packet.
fn scan_pos(dir: &str, prefix: &str, positions_path: &str, path: u32, lo_bp: u64, hi_bp: u64) -> io::Result<()> {
    let t0 = Instant::now();
    let positions = Positions::load(positions_path)?;
    let packets = read_packets(&shard_gz(dir, prefix))?;
    let mut hits = 0usize;
    for packet in &packets {
        let (nodes, _, _) = walk_nodes(packet)?;
        if nodes.iter().any(|&n| matches!(positions.position(n), Some((p, o)) if p == path && o >= lo_bp && o <= hi_bp)) {
            hits += 1;
        }
    }
    println!("scan     (path {path}, {lo_bp}..{hi_bp}): {} reads examined, {hits} with evidence, {:.1} ms (whole shard decompressed)",
        packets.len(), t0.elapsed().as_secs_f64()*1e3);
    Ok(())
}

/// Stream one shard's packets, one at a time (the whole sample is 46 GB
/// gzipped; only the per-bucket spills are held on disk).
fn stream_packets(path: &Path) -> io::Result<Box<dyn FnMut() -> io::Result<Option<Vec<u8>>>>> {
    let file = File::open(path)?;
    let gz = flate2::read::GzDecoder::new(BufReader::with_capacity(1 << 20, file));
    let mut reader = BufReader::with_capacity(1 << 20, gz);
    Ok(Box::new(move || {
        let mut len_buf = [0u8; 4];
        match reader.read_exact(&mut len_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(e),
        }
        let len = u32::from_le_bytes(len_buf) as usize;
        let mut packet = vec![0u8; 4 + len];
        packet[..4].copy_from_slice(&len_buf);
        reader.read_exact(&mut packet[4..])?;
        Ok(Some(packet))
    }))
}

/// THE WHOLE-SAMPLE BUILDER: every region shard, one position-ordered
/// store. Two passes, because 46 GB of packets cannot be held: pass 1
/// spills each (key, read, packet) into one of NSPILL files chosen by
/// key range; pass 2 sorts each spill by key and appends it to the store
/// in order, recording the bucket boundaries. This is an external
/// bucketed sort, and it is the shape the projection emission should
/// eventually produce directly.
/// MATERIALISE ONE TERRITORY'S EVIDENCE as a single virtual shard.
///
/// This needs NO product change: the cut's scan already reads
/// `shard.*.bin.gz` from `--shards-dir`, so a directory holding one shard
/// built from the evidence store (the territory's reads, sorted by their
/// arrival ordinal so the scan sees them in the same relative order) makes
/// the existing pipeline fetch instead of scan. That is the measurement
/// the 66-hour question needs, and the gate is cut-census.jsonl +
/// cut-cache.bin byte-identity against the full-shard run.
fn materialize(
    store_prefix: &str,
    positions_path: &str,
    names_path: &str,
    beds_dir: &str,
    out_dir: &str,
    width_bp: u64,
) -> io::Result<()> {
    use std::collections::HashMap as Map;
    use std::io::{Seek, SeekFrom};
    let t0 = Instant::now();
    // name -> path index, from the panel's own name list (idx TAB name ...)
    let mut name_to_idx: Map<String, u32> = Map::new();
    for line in std::fs::read_to_string(names_path)?.lines() {
        let mut fields = line.split('\t');
        if let (Some(idx), Some(name)) = (fields.next(), fields.next()) {
            if let Ok(idx) = idx.parse::<u32>() {
                name_to_idx.insert(name.to_string(), idx);
            }
        }
    }
    // the store's manifest
    let mut head = Vec::new();
    File::open(format!("{store_prefix}.evd.idx"))?.read_to_end(&mut head)?;
    let is_compressed = &head[..8] == MAGIC_Z;
    let u64_at = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().unwrap());
    let width = u64_at(8);
    let n_buckets = u64_at(16) as usize;
    let mut table: Vec<(u64, u64, u64)> = Vec::with_capacity(n_buckets);
    for i in 0..n_buckets {
        let base = 40 + i * 24;
        table.push((u64_at(base), u64_at(base + 8), u64_at(base + 16)));
    }
    let spine: Vec<u64> = table.iter().map(|&(k, _, _)| k).collect();
    let _ = positions_path; // the windows are already in path coordinates

    // the territory's windows, from the partition beds (path_name TAB start TAB end)
    let mut windows: Vec<(u32, u64, u64)> = Vec::new();
    let mut unmatched = 0usize;
    for entry in std::fs::read_dir(beds_dir)? {
        let path = entry?.path();
        if path.extension().and_then(|e| e.to_str()) != Some("bed") {
            continue;
        }
        for line in std::fs::read_to_string(&path)?.lines() {
            let mut fields = line.split('\t');
            let (Some(name), Some(start), Some(end)) = (fields.next(), fields.next(), fields.next())
            else {
                continue;
            };
            // the bed names carry a ":start-end" suffix; the panel's do not
            let base = name.split(':').next().unwrap_or(name);
            let (Ok(start), Ok(end)) = (start.parse::<u64>(), end.parse::<u64>()) else {
                continue;
            };
            match name_to_idx.get(base) {
                Some(&idx) => windows.push((idx, start, end)),
                None => unmatched += 1,
            }
        }
    }
    println!(
        "  windows {} (unmatched path names {unmatched}) from {beds_dir}",
        windows.len()
    );

    // fetch: every bucket any window covers
    let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let mut packets: Vec<(u32, Vec<u8>)> = Vec::new();
    let mut file = File::open(format!("{store_prefix}.evd"))?;
    let mut bytes = 0u64;
    let mut bucket_visits = 0usize;
    for &(path_idx, lo_bp, hi_bp) in &windows {
        let lo_key = position_key(path_idx, lo_bp, width);
        let hi_key = position_key(path_idx, hi_bp, width);
        let start_at = spine.partition_point(|&k| k < lo_key);
        for &(key, offset, len) in &table[start_at..] {
            if key > hi_key {
                break;
            }
            file.seek(SeekFrom::Start(offset))?;
            let mut raw = vec![0u8; len as usize];
            file.read_exact(&mut raw)?;
            bytes += len as u64;
            bucket_visits += 1;
            let buf = if is_compressed {
                zstd::stream::decode_all(std::io::Cursor::new(&raw))?
            } else {
                raw
            };
            let mut at = 0usize;
            while at + 8 <= buf.len() {
                let ordinal = u32::from_le_bytes(buf[at..at + 4].try_into().unwrap());
                let plen = u32::from_le_bytes(buf[at + 4..at + 8].try_into().unwrap()) as usize;
                let end = at + 8 + plen;
                if end > buf.len() {
                    break;
                }
                if seen.insert(ordinal) {
                    packets.push((ordinal, buf[at + 4..end].to_vec()));
                }
                at = end;
            }
        }
    }
    println!(
        "  fetched {} packets ({} duplicated placements dropped), {:.1} MB read, {bucket_visits} bucket visits",
        packets.len(),
        seen.len() - packets.len(),
        bytes as f64 / 1e6
    );

    // the scan sees reads in ARRIVAL order: sort by ordinal, then write one shard
    packets.sort_unstable_by_key(|p| p.0);
    std::fs::create_dir_all(out_dir)?;
    {
        let out_path = Path::new(out_dir).join("shard.evidence.bin.gz");
        let mut encoder =
            flate2::write::GzEncoder::new(BufWriter::new(File::create(&out_path)?), flate2::Compression::fast());
        for (_, packet) in &packets {
            encoder.write_all(&packet[..4])?;
            encoder.write_all(&packet[4..])?;
        }
        encoder.finish()?;
        println!(
            "  wrote {} ({:.1} MB) in {:.1}s",
            out_path.display(),
            std::fs::metadata(&out_path)?.len() as f64 / 1e6,
            t0.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

/// MATERIALISE A TERRITORY BY ITS NODE SET - the exact query.
///
/// The bed windows are NOT the territory: the cut's predicate is "does a
/// walk node belong to the territory's node set", and that set extends
/// past the window edges (measured: bed windows fetched 488 of the 5,678
/// records the scan produces, a strict subset - zero extra). The node set
/// is what the cut itself computes as its expressibility map, BEFORE the
/// scan, so querying by it is both exact and cheap: each node resolves to
/// one (path, offset) via the panel position table, hence to one bucket,
/// and a read touching that node is stored in that bucket by construction.
fn materialize_nodes(
    store_prefix: &str,
    positions_path: &str,
    nodes_path: &str,
    out_dir: &str,
) -> io::Result<()> {
    use std::io::{Seek, SeekFrom};
    let t0 = Instant::now();
    let positions = Positions::load(positions_path)?;
    // the width comes from the store's own manifest
    let mut head = Vec::new();
    File::open(format!("{store_prefix}.evd.idx"))?.read_to_end(&mut head)?;
    let is_compressed = &head[..8] == MAGIC_Z;
    let u64_at = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().unwrap());
    let width = u64_at(8);
    let n_buckets = u64_at(16) as usize;
    let mut table: Vec<(u64, u64, u64)> = Vec::with_capacity(n_buckets);
    for i in 0..n_buckets {
        let base = 40 + i * 24;
        table.push((u64_at(base), u64_at(base + 8), u64_at(base + 16)));
    }
    let spine: Vec<u64> = table.iter().map(|&(k, _, _)| k).collect();

    // the territory's node ids -> the set of bucket keys they occupy
    let text = std::fs::read_to_string(nodes_path)?;
    let mut keys: Vec<u64> = Vec::new();
    let mut resolved = 0usize;
    let mut unresolved = 0usize;
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let Ok(node) = line.trim().parse::<i32>() else { continue };
        match positions.position(node) {
            Some((path, offset)) => {
                keys.push(position_key(path, offset, width));
                resolved += 1;
            }
            None => unresolved += 1,
        }
    }
    keys.sort_unstable();
    keys.dedup();
    println!(
        "  territory nodes resolved {resolved} (unresolved {unresolved}) -> {} distinct buckets",
        keys.len()
    );

    // fetch exactly those buckets
    let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let mut packets: Vec<(u32, Vec<u8>)> = Vec::new();
    let mut file = File::open(format!("{store_prefix}.evd"))?;
    let mut bytes = 0u64;
    for &key in &keys {
        let at = spine.partition_point(|&k| k < key);
        if at >= spine.len() || spine[at] != key {
            continue;
        }
        let (_, offset, len) = table[at];
        file.seek(SeekFrom::Start(offset))?;
        let mut raw = vec![0u8; len as usize];
        file.read_exact(&mut raw)?;
        bytes += len as u64;
        let buf = if is_compressed {
            zstd::stream::decode_all(std::io::Cursor::new(&raw))?
        } else {
            raw
        };
        let mut p = 0usize;
        while p + 8 <= buf.len() {
            let ordinal = u32::from_le_bytes(buf[p..p + 4].try_into().unwrap());
            let plen = u32::from_le_bytes(buf[p + 4..p + 8].try_into().unwrap()) as usize;
            let end = p + 8 + plen;
            if end > buf.len() {
                break;
            }
            if seen.insert(ordinal) {
                packets.push((ordinal, buf[p + 4..end].to_vec()));
            }
            p = end;
        }
    }
    println!(
        "  fetched {} packets, {:.1} MB read, {} buckets",
        packets.len(),
        bytes as f64 / 1e6,
        keys.len()
    );
    packets.sort_unstable_by_key(|p| p.0);
    std::fs::create_dir_all(out_dir)?;
    let out_path = Path::new(out_dir).join("shard.evidence.bin.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        BufWriter::new(File::create(&out_path)?),
        flate2::Compression::fast(),
    );
    for (_, packet) in &packets {
        encoder.write_all(packet)?;
    }
    encoder.finish()?;
    println!(
        "  wrote {} ({:.1} MB) in {:.1}s",
        out_path.display(),
        std::fs::metadata(&out_path)?.len() as f64 / 1e6,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn build_all(dir: &str, positions_path: &str, out_prefix: &str, width_bp: u64) -> io::Result<()> {
    const NSPILL: u64 = 256;
    let t0 = Instant::now();
    let positions = Positions::load(positions_path)?;
    println!("positions {} nodes ({:.1}s)", positions.entries.len(), t0.elapsed().as_secs_f64());

    let mut shards: Vec<PathBuf> = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name.starts_with("shard.") && name.ends_with(".bin.gz") {
            shards.push(path);
        }
    }
    shards.sort();
    println!("shards {}", shards.len());

    // ---- pass 1: spill by key range ----
    let t1 = Instant::now();
    let spill_dir = Path::new(dir).join("evd-spill");
    std::fs::create_dir_all(&spill_dir)?;
    let mut spills: Vec<BufWriter<File>> = (0..NSPILL)
        .map(|i| File::create(spill_dir.join(format!("{i:03}.spill"))).map(BufWriter::new))
        .collect::<io::Result<_>>()?;
    let mut read_id = 0u64;
    let mut placed = 0u64;
    let mut unplaced = 0u64;
    for shard in &shards {
        let mut next = stream_packets(shard)?;
        while let Some(packet) = next()? {
            let ordinal = read_id;
            read_id += 1;
            let (nodes, _, _) = walk_nodes(&packet)?;
            let mut keys: Vec<u64> = Vec::new();
            for &node in &nodes {
                if let Some((path, offset)) = positions.position(node) {
                    keys.push(position_key(path, offset, width_bp));
                }
            }
            if keys.is_empty() {
                unplaced += 1;
                continue;
            }
            keys.sort_unstable();
            keys.dedup();
            placed += 1;
            for key in keys {
                let spill = &mut spills[(key % NSPILL) as usize];
                spill.write_all(&key.to_le_bytes())?;
                spill.write_all(&ordinal.to_le_bytes())?;
                spill.write_all(&(packet.len() as u32).to_le_bytes())?;
                spill.write_all(&packet)?;
            }
        }
    }
    drop(spills);
    let spill_time = t1.elapsed();

    // ---- pass 2: merge each spill in key order ----
    let t2 = Instant::now();
    let evd = Path::new(dir).join(format!("{out_prefix}.evd"));
    let idx = Path::new(dir).join(format!("{out_prefix}.evd.idx"));
    let sum = Path::new(dir).join(format!("{out_prefix}.evd.sum"));
    let mut table: Vec<(u64, u64, u64)> = Vec::new();
    let mut offset = 0u64;
    {
        let mut out = BufWriter::new(File::create(&evd)?);
        for i in 0..NSPILL {
            let path = spill_dir.join(format!("{i:03}.spill"));
            let mut buf = Vec::new();
            File::open(&path)?.read_to_end(&mut buf)?;
            if buf.is_empty() {
                continue;
            }
            // decode (key, ordinal, packet) triples, then sort by key
            let mut entries: Vec<(u64, u64, Vec<u8>)> = Vec::new();
            let mut at = 0usize;
            while at + 16 <= buf.len() {
                let key = u64::from_le_bytes(buf[at..at + 8].try_into().unwrap());
                let ordinal = u64::from_le_bytes(buf[at + 8..at + 16].try_into().unwrap());
                let plen = u32::from_le_bytes(buf[at + 16..at + 20].try_into().unwrap()) as usize;
                let end = at + 20 + plen;
                if end > buf.len() {
                    break;
                }
                entries.push((key, ordinal, buf[at + 20..end].to_vec()));
                at = end;
            }
            entries.sort_unstable_by_key(|e| (e.0, e.1));
            let mut current: Option<(u64, u64)> = None; // (key, bucket start offset)
            let mut bucket: Vec<u8> = Vec::new();
            let mut flush = |bucket: &mut Vec<u8>,
                             key: u64,
                             start: u64,
                             out: &mut BufWriter<File>,
                             offset: &mut u64,
                             table: &mut Vec<(u64, u64, u64)>|
             -> io::Result<()> {
                let compressed = zstd::bulk::compress(bucket, 3)?;
                out.write_all(&compressed)?;
                *offset += compressed.len() as u64;
                table.push((key, start, compressed.len() as u64));
                bucket.clear();
                Ok(())
            };
            for (key, ordinal, packet) in entries {
                if current.map(|(k, _)| k) != Some(key) {
                    if let Some((k, start)) = current {
                        flush(&mut bucket, k, start, &mut out, &mut offset, &mut table)?;
                    }
                    current = Some((key, offset));
                }
                bucket.extend_from_slice(&(ordinal as u32).to_le_bytes());
                bucket.extend_from_slice(&packet);
            }
            if let Some((k, start)) = current {
                flush(&mut bucket, k, start, &mut out, &mut offset, &mut table)?;
            }
        }
        out.flush()?;
    }
    table.sort_unstable_by_key(|e| e.0);
    let merge_time = t2.elapsed();

    // ---- manifests ----
    let t3 = Instant::now();
    {
        let mut out = BufWriter::new(File::create(&idx)?);
        out.write_all(MAGIC_Z)?;
        out.write_all(&width_bp.to_le_bytes())?;
        out.write_all(&(table.len() as u64).to_le_bytes())?;
        out.write_all(&read_id.to_le_bytes())?;
        out.write_all(&0u64.to_le_bytes())?;
        for &(bucket, start, len) in &table {
            out.write_all(&bucket.to_le_bytes())?;
            out.write_all(&start.to_le_bytes())?;
            out.write_all(&len.to_le_bytes())?;
        }
    }
    {
        // the O(1) layer: per bucket, an approximate read count is enough
        // for plan(); exact counts come from the summary pass at emission
        let mut out = BufWriter::new(File::create(&sum)?);
        for &(bucket, start, len) in &table {
            out.write_all(&bucket.to_le_bytes())?;
            out.write_all(&0u64.to_le_bytes())?; // reads
            out.write_all(&len.to_le_bytes())?; // bytes (the usable number)
            out.write_all(&0u64.to_le_bytes())?; // tokens
        }
    }
    let write_time = t3.elapsed();
    let evd_bytes = std::fs::metadata(&evd)?.len();
    let idx_bytes = std::fs::metadata(&idx)?.len();
    let sum_bytes = std::fs::metadata(&sum)?.len();
    println!("  reads {read_id}  placed {placed}  unplaced {unplaced}");
    println!("  buckets {}", table.len());
    println!(
        "  pass1 spill {:.1}s  pass2 merge {:.1}s  manifests {:.1}s",
        spill_time.as_secs_f64(),
        merge_time.as_secs_f64(),
        write_time.as_secs_f64()
    );
    println!(
        "  store {:.1} MB   index {:.2} MB   summary {:.2} MB   (overhead {:.4}%)",
        evd_bytes as f64 / 1e6,
        idx_bytes as f64 / 1e6,
        sum_bytes as f64 / 1e6,
        100.0 * (idx_bytes + sum_bytes) as f64 / evd_bytes.max(1) as f64
    );
    Ok(())
}

fn evidence(dir: &str, prefix: &str, lo_node: u64, hi_node: u64) -> io::Result<()> {
    let width = {
        let (w, _, _) = load_manifest(dir, prefix)?;
        w
    };
    let (lo, hi) = (lo_node / width, hi_node / width);
    let t0 = Instant::now();
    let (packets, buckets_read, bytes_read) = read_buckets(dir, prefix, lo, hi)?;
    // filter: the read must have SOME walk node inside the requested range
    let mut hits = 0usize;
    for packet in &packets {
        let (nodes, _, _) = walk_nodes(packet)?;
        if nodes
            .iter()
            .any(|&n| (n.unsigned_abs() as u64) >= lo_node && (n.unsigned_abs() as u64) <= hi_node)
        {
            hits += 1;
        }
    }
    let elapsed = t0.elapsed();
    println!(
        "evidence [{}..{}] nodes (buckets {}..{}): {} reads examined, {} with evidence, {:.3} MB read, {} buckets, {:.2} ms",
        lo_node,
        hi_node,
        lo,
        hi,
        packets.len(),
        hits,
        bytes_read as f64 / 1e6,
        buckets_read,
        elapsed.as_secs_f64() * 1e3
    );
    Ok(())
}

fn summary(dir: &str, prefix: &str, lo_node: u64, hi_node: u64) -> io::Result<()> {
    let (width, table, rewind) = load_manifest(dir, prefix)?;
    let (lo, hi) = (lo_node / width, hi_node / width);
    let mut buf = Vec::new();
    File::open(Path::new(dir).join(format!("shard.{prefix}.evd.sum")))?.read_to_end(&mut buf)?;
    let u64_at = |at: usize| u64::from_le_bytes(buf[at..at + 8].try_into().unwrap());
    let first = lo.saturating_sub(rewind);
    let spine: Vec<u64> = table.iter().map(|&(b, _, _)| b).collect();
    let start = spine.partition_point(|&b| b < first);
    let (mut reads, mut records, mut tokens, mut buckets) = (0u64, 0u64, 0u64, 0usize);
    for (index, &(bucket, _, _)) in table.iter().enumerate().skip(start) {
        if bucket > hi {
            break;
        }
        let base = index * 32;
        if base + 32 <= buf.len() {
            reads += u64_at(base + 8);
            records += u64_at(base + 16);
            tokens += u64_at(base + 24);
            buckets += 1;
        }
    }
    let payload = idx_payload(dir, prefix)? as f64 / 1e6;
    println!(
        "summary [{}..{}]: {reads} reads, {records} records, {tokens} tokens, {buckets} buckets touched; query cost ~{} KB (no packets), full store {payload:.1} MB",
        lo_node, hi_node, buckets * 32 / 1024
    );
    Ok(())
}

fn idx_payload(dir: &str, prefix: &str) -> io::Result<u64> {
    Ok(std::fs::metadata(Path::new(dir).join(format!("shard.{prefix}.evd")))?.len())
}

fn scan(dir: &str, prefix: &str, lo_node: u64, hi_node: u64) -> io::Result<()> {
    // the identity gate's own arm: every packet, no index
    let t0 = Instant::now();
    let packets = read_packets(&shard_gz(dir, prefix))?;
    let mut hits = 0usize;
    for packet in &packets {
        let (nodes, _, _) = walk_nodes(packet)?;
        if nodes
            .iter()
            .any(|&n| (n.unsigned_abs() as u64) >= lo_node && (n.unsigned_abs() as u64) <= hi_node)
        {
            hits += 1;
        }
    }
    println!(
        "scan     [{}..{}]: {} reads examined, {} with evidence, in {:.1} ms (the whole shard, decompressed)",
        lo_node,
        hi_node,
        packets.len(),
        hits,
        t0.elapsed().as_secs_f64() * 1e3
    );
    Ok(())
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let arg = |i: usize| args.get(i).cloned().unwrap_or_default();
    match args.get(1).map(String::as_str) {
        Some("build") => build(&arg(2), &arg(3)),
        Some("build-keyed") => build_keyed(
            &arg(2),
            &arg(3),
            &arg(4),
            arg(5).parse().unwrap_or(10_000),
        ),
        Some("evidence") => evidence(&arg(2), &arg(3), arg(4).parse().unwrap_or(0), arg(5).parse().unwrap_or(0)),
        Some("build-all") => build_all(
            &arg(2),
            &arg(3),
            &arg(4),
            arg(5).parse().unwrap_or(10_000),
        ),
        Some("materialize-nodes") => materialize_nodes(&arg(2), &arg(3), &arg(4), &arg(5)),
        Some("materialize") => materialize(
            &arg(2),
            &arg(3),
            &arg(4),
            &arg(5),
            &arg(6),
            arg(7).parse().unwrap_or(10_000),
        ),
        Some("evidence-pos") => evidence_pos(&arg(2), &arg(3), &arg(4), arg(5).parse().unwrap_or(0), arg(6).parse().unwrap_or(0), arg(7).parse().unwrap_or(0)),
        Some("scan-pos") => scan_pos(&arg(2), &arg(3), &arg(4), arg(5).parse().unwrap_or(0), arg(6).parse().unwrap_or(0), arg(7).parse().unwrap_or(0)),
        Some("summary") => summary(&arg(2), &arg(3), arg(4).parse().unwrap_or(0), arg(5).parse().unwrap_or(0)),
        Some("scan") => scan(&arg(2), &arg(3), arg(4).parse().unwrap_or(0), arg(5).parse().unwrap_or(0)),
        _ => {
            eprintln!("usage: evidence_index build|evidence|summary|scan <dir> <shard-prefix> [lo hi]");
            Ok(())
        }
    }
}
