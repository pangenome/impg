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

fn load_manifest(dir: &str, prefix: &str) -> io::Result<(u64, Vec<(u64, u64, u64)>, u64)> {
    let mut buf = Vec::new();
    File::open(Path::new(dir).join(format!("shard.{prefix}.evd.idx")))?.read_to_end(&mut buf)?;
    if buf.len() < 32 || &buf[..8] != MAGIC {
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
    File::open(Path::new(dir).join(format!("shard.{prefix}.p.evd.idx")))?.read_to_end(&mut head)?;
    let u64_at = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().unwrap());
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
    let mut file = File::open(Path::new(dir).join(format!("shard.{prefix}.p.evd")))?;
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
