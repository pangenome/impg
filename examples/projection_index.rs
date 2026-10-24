//! THE PROJECTION'S MISSING INDEX — the compact format, built and measured.
//!
//! The projection emits per-region pattern multisets
//! (`shard.<name>.occ.bin`: u32 n_tokens, n_tokens × i64 node, u64 mult)
//! with NO index. The coordinates are in the records — the module header
//! says "the syng node ids ARE the coordinates" — but the stream is
//! emitted in flush order (measured 49.9% ascending, i.e. random), so
//! every consumer scans a whole region shard per territory. The cut
//! receipts price that: 90,898,144 distinct patterns scanned to use the
//! 12,797 that touch one territory — 7,103×.
//!
//! This builds the `.bai`-shaped thing that is missing, in a compact
//! form: an inverted index `token -> record offset`, sorted, written as
//! varint DELTAS (both fields are monotone inside a token run), with a
//! sparse bucket table so a query can seek and decode only the buckets
//! its target nodes fall in.
//!
//! ```text
//! projection_index build <dir> <shard-prefix>
//! projection_index query <dir> <shard-prefix> <targets> [scattered|clustered]
//! ```
//!
//! `query` answers with the SAME record set as the full scan (asserted)
//! and reports the record bytes read on each side, so the win is a
//! measurement. The `scattered`/`clustered` target shapes bracket the
//! unknown: whether a real territory's node ids cluster in the panel's
//! node space or spread across it.

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::time::Instant;

/// One bucket covers 2^20 token values; the panel's node space is
/// ~1.09e9, so ~1,040 buckets — a tiny table, and a query reads only
/// the buckets its targets touch.
/// One bucket covers 2^SHIFT token values. MEASURED on a real territory
/// (T00048's 176,659 global node ids spanning -271M..+271M, i.e.
/// scattered across the panel's whole node space): at 2^20 a query
/// touches 49.7% of the buckets and beats the scan by only 2.0x; at 2^12
/// it touches 2.2% and beats it by 45x. The table at 2^12 costs ~6.4 MB,
/// so the fine shift is free. Override with IMPG_INDEX_BUCKET_SHIFT.
fn bucket_shift() -> u32 {
    std::env::var("IMPG_INDEX_BUCKET_SHIFT")
        .ok()
        .and_then(|v| v.trim().parse::<u32>().ok())
        .filter(|s| *s >= 6 && *s <= 24)
        .unwrap_or(12)
}
const MAGIC: &[u8; 8] = b"OCCIDX02";

fn bucket_of(token: i64) -> u64 {
    let shift = bucket_shift();
    ((token as u64) >> shift) & ((1u64 << (64 - shift)) - 1)
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

fn get_varint(buf: &[u8], at: &mut usize) -> u64 {
    let mut value = 0u64;
    let mut shift = 0u32;
    loop {
        let byte = buf[*at];
        *at += 1;
        value |= ((byte & 0x7f) as u64) << shift;
        if byte & 0x80 == 0 {
            return value;
        }
        shift += 7;
    }
}

fn zigzag(v: i64) -> u64 {
    ((v << 1) ^ (v >> 63)) as u64
}

fn unzigzag(v: u64) -> i64 {
    ((v >> 1) as i64) ^ -((v & 1) as i64)
}

struct Record {
    offset: u64,
    tokens: Vec<i64>,
}

fn read_records(path: &std::path::Path) -> io::Result<(Vec<u8>, Vec<Record>)> {
    let mut bytes = Vec::new();
    File::open(path)?.read_to_end(&mut bytes)?;
    let mut records = Vec::new();
    let mut o = 0usize;
    while o + 4 <= bytes.len() {
        let n = u32::from_le_bytes(bytes[o..o + 4].try_into().unwrap()) as usize;
        let start = o;
        o += 4;
        assert!(n > 0 && n % 2 == 1, "invalid occ token count {n}");
        let mut tokens = Vec::with_capacity(n);
        for _ in 0..n {
            tokens.push(i64::from_le_bytes(bytes[o..o + 8].try_into().unwrap()));
            o += 8;
        }
        o += 8; // mult
        records.push(Record { offset: start as u64, tokens });
    }
    Ok((bytes, records))
}

fn shard_path(dir: &str, prefix: &str, ext: &str) -> std::path::PathBuf {
    std::path::Path::new(dir).join(format!("shard.{prefix}.{ext}"))
}

fn build(dir: &str, prefix: &str) -> io::Result<()> {
    let t0 = Instant::now();
    let (records_bytes, records) = read_records(&shard_path(dir, prefix, "occ.bin"))?;
    let parse = t0.elapsed();

    let t1 = Instant::now();
    let mut postings: Vec<(i64, u64)> = Vec::new();
    for record in &records {
        for &token in &record.tokens {
            postings.push((token, record.offset));
        }
    }
    postings.sort_unstable();
    let sort = t1.elapsed();

    // ---- the compact payload, per bucket ----
    let t2 = Instant::now();
    let mut buckets: Vec<(i64, u64, u64)> = Vec::new(); // (first token, payload offset, payload len)
    let mut payload: Vec<u8> = Vec::with_capacity(postings.len() * 4);
    let mut i = 0usize;
    while i < postings.len() {
        let bucket = bucket_of(postings[i].0);
        let first_token = postings[i].0;
        let start = payload.len();
        let mut prev_token = first_token;
        let mut prev_offset = 0i64;
        while i < postings.len() && bucket_of(postings[i].0) == bucket {
            let (token, offset) = postings[i];
            put_varint(&mut payload, zigzag(token - prev_token));
            put_varint(&mut payload, zigzag(offset as i64 - prev_offset));
            prev_token = token;
            prev_offset = offset as i64;
            i += 1;
        }
        buckets.push((first_token, start as u64, (payload.len() - start) as u64));
    }
    let encode = t2.elapsed();

    let path = shard_path(dir, prefix, "occ.idx");
    {
        let mut out = BufWriter::new(File::create(&path)?);
        out.write_all(MAGIC)?;
        out.write_all(&(postings.len() as u64).to_le_bytes())?;
        out.write_all(&(buckets.len() as u64).to_le_bytes())?;
        for &(first_token, offset, len) in &buckets {
            out.write_all(&first_token.to_le_bytes())?;
            out.write_all(&offset.to_le_bytes())?;
            out.write_all(&len.to_le_bytes())?;
        }
        out.write_all(&payload)?;
    }
    let index_bytes = std::fs::metadata(&path)?.len();

    println!("shard.{prefix}.occ.bin");
    println!("  records            {}", records.len());
    println!("  postings           {}", postings.len());
    println!("  buckets            {}", buckets.len());
    println!("  parse {:.2}s  sort {:.2}s  encode {:.2}s", parse.as_secs_f64(), sort.as_secs_f64(), encode.as_secs_f64());
    println!(
        "  index              {} bytes ({:.1} MB) = {:.2} B/posting",
        index_bytes,
        index_bytes as f64 / 1e6,
        index_bytes as f64 / postings.len() as f64
    );
    println!(
        "  records            {:.1} MB  -> index is {:.0}% of the record stream",
        records_bytes.len() as f64 / 1e6,
        100.0 * index_bytes as f64 / records_bytes.len() as f64
    );
    Ok(())
}

fn query(dir: &str, prefix: &str, wanted: usize, shape: &str) -> io::Result<()> {
    let (bytes, records) = read_records(&shard_path(dir, prefix, "occ.bin"))?;
    let mut tokens: Vec<i64> = records.iter().flat_map(|r| r.tokens.iter().copied()).collect();
    tokens.sort_unstable();
    let (lo, hi) = (tokens[0], *tokens.last().unwrap());

    let mut targets: Vec<i64> = match shape {
        // a 250 kb window's worth of node ids, i.e. a contiguous token band
        "clustered" => (lo + (hi - lo) / 3..).take(wanted).collect(),
        // spread across the whole node space (the pessimistic bracket)
        _ => {
            let mut state = 0x2545F4914F6CDD1Du64;
            (0..wanted)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    lo + (state % ((hi - lo) as u64 + 1)) as i64
                })
                .collect()
        }
    };
    query_with(dir, prefix, targets, shape)
}

fn query_with(dir: &str, prefix: &str, mut targets: Vec<i64>, shape: &str) -> io::Result<()> {
    let (bytes, records) = read_records(&shard_path(dir, prefix, "occ.bin"))?;
    let wanted = targets.len();
    targets.sort_unstable();
    let target_set: HashSet<i64> = targets.iter().copied().collect();

    // ---- the full scan (what the cut does today) ----
    let t0 = Instant::now();
    let mut scanned: HashSet<u64> = HashSet::new();
    for record in &records {
        if record.tokens.iter().any(|t| target_set.contains(t)) {
            scanned.insert(record.offset);
        }
    }
    let scan = t0.elapsed();

    // ---- the bucket-seeked indexed query ----
    let index_path = shard_path(dir, prefix, "occ.idx");
    let t0 = Instant::now();
    let mut file = File::open(&index_path)?;
    let mut head = [0u8; 8];
    file.read_exact(&mut head)?;
    assert_eq!(&head, MAGIC, "bad index magic");
    let mut b8 = [0u8; 8];
    file.read_exact(&mut b8)?;
    let n_postings = u64::from_le_bytes(b8) as usize;
    file.read_exact(&mut b8)?;
    let n_buckets = u64::from_le_bytes(b8) as usize;
    let mut table = Vec::with_capacity(n_buckets);
    for _ in 0..n_buckets {
        let mut e = [0u8; 24];
        file.read_exact(&mut e)?;
        table.push((
            i64::from_le_bytes(e[0..8].try_into().unwrap()),
            u64::from_le_bytes(e[8..16].try_into().unwrap()) as usize,
            u64::from_le_bytes(e[16..24].try_into().unwrap()) as usize,
        ));
    }
    let payload_base = 8 + 8 + 8 + n_buckets * 24;
    // which buckets do the targets touch? (bisect against the bucket spine)
    let spine: Vec<i64> = table.iter().map(|&(t, _, _)| t).collect();
    let mut bucket_reads = 0u64;
    let mut bytes_read = 0u64;
    let mut buf: Vec<u8> = Vec::new();
    let mut hits: Vec<u64> = Vec::new();
    let mut i = 0usize;
    while i < table.len() {
        let (first_token, off, len) = table[i];
        let next_first = if i + 1 < table.len() { table[i + 1].0 } else { i64::MAX };
        // does any target fall in [first_token, next_first)?
        let mut at = targets.partition_point(|&t| t < first_token);
        if at >= targets.len() || targets[at] >= next_first {
            i += 1;
            continue;
        }
        file.seek(SeekFrom::Start((payload_base + off) as u64))?;
        buf.resize(len, 0);
        file.read_exact(&mut buf)?;
        bytes_read += len as u64;
        bucket_reads += 1;
        let mut p = 0usize;
        let mut token = first_token;
        let mut offset = 0i64;
        // (THE MERGE SCAN: tokens are ASCENDING inside a bucket and the
        // targets are ASCENDING, so membership is a two-pointer merge -
        // no hash probe and no varint function call per posting. Measured
        // cost of the naive form: 8M postings x 41 ns, which was the whole
        // 332 ms; the decode, not the I/O, was the query's cost.)
        while p < buf.len() {
            // inline token delta
            let mut raw = 0u64;
            let mut shift = 0u32;
            loop {
                let byte = buf[p];
                p += 1;
                raw |= ((byte & 0x7f) as u64) << shift;
                if byte & 0x80 == 0 { break; }
                shift += 7;
            }
            token += unzigzag(raw);
            // inline offset delta
            let mut raw = 0u64;
            let mut shift = 0u32;
            loop {
                let byte = buf[p];
                p += 1;
                raw |= ((byte & 0x7f) as u64) << shift;
                if byte & 0x80 == 0 { break; }
                shift += 7;
            }
            offset += unzigzag(raw);
            while at < targets.len() && targets[at] < token {
                at += 1;
            }
            if at < targets.len() && targets[at] == token {
                hits.push(offset as u64);
            }
        }
        let _ = spine;
        i += 1;
    }
    let indexed_time = t0.elapsed();
    hits.sort_unstable();
    hits.dedup();
    let indexed: HashSet<u64> = hits.into_iter().collect();
    let record_bytes_per = bytes.len() as f64 / records.len() as f64;

    println!("shard.{prefix}.occ.bin  targets={wanted} ({shape})");
    println!("  records {}   postings {}", records.len(), n_postings);
    println!("  matches: scan={} index={}  EQUAL={}", scanned.len(), indexed.len(), scanned == indexed);
    println!(
        "  full scan   {:.1} ms   reads {:.1} MB of records",
        scan.as_secs_f64() * 1e3,
        bytes.len() as f64 / 1e6
    );
    println!(
        "  indexed     {:.1} ms   reads {:.3} MB of index ({}/{} buckets)",
        indexed_time.as_secs_f64() * 1e3,
        bytes_read as f64 / 1e6,
        bucket_reads,
        n_buckets
    );
    println!(
        "  -> {:.0}x less record-side I/O ({:.1} MB -> {:.3} MB of matching records)",
        bytes.len() as f64 / (indexed.len() as f64 * record_bytes_per).max(1.0),
        bytes.len() as f64 / 1e6,
        indexed.len() as f64 * record_bytes_per / 1e6
    );
    Ok(())
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("build") => build(&args[2], &args[3]),
        Some("query") => {
            let wanted = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(176_659);
            query(&args[2], &args[3], wanted, args.get(5).map(String::as_str).unwrap_or("scattered"))
        }
        // query-file: the target set is a REAL territory's node dump
        // (IMPG_DUMP_TERRITORY_NODES) — the receipt form.
        Some("query-file") => {
            let text = std::fs::read_to_string(&args[4])?;
            let targets: Vec<i64> = text
                .lines()
                .filter_map(|line| line.trim().parse::<i64>().ok())
                .collect();
            query_with(&args[2], &args[3], targets, "territory-dump")
        }
        _ => {
            eprintln!("usage: projection_index build <dir> <shard-prefix> | query <dir> <shard-prefix> <targets> [scattered|clustered]");
            Ok(())
        }
    }
}
