//! THE ONE WHOLE-SAMPLE PROJECTION (the owner's architecture ruling):
//! project the WHOLE sample onto the syng graph ONCE, then cut
//! locus-by-locus — the per-locality GRCh38-BAM fetch-and-reproject is
//! the retired crutch.
//!
//! Three receipt-side modes (public API only; no product/scoreboard
//! changes; the scoring machinery is untouched):
//!
//! `project` — THE ONE PROJECTION. Streams the sample's reads once
//! through the MEM machinery against the GLOBAL panel syng (the FM
//! index of walks in the GBWT): every read's complete maximal-MEM
//! records, each an (signed syng node, read position) walk — the syng
//! node ids ARE the coordinates. Emits per input shard:
//!   `shard.<name>.bin.gz` — per read, in stream order: sequence,
//!      quality, and every record's canonical tokens plus its
//!      own-frame anchor walk (the derivation
//!      `mem_records::tagged_mem_records` returns — the same accessor
//!      the committed census re-derivation and derive-cache builder
//!      use).
//!   `shard.<name>.occ.bin` — the shard's pattern multiset (distinct
//!      canonical MEM walks with read multiplicities), deduplicated
//!      in-shard with a bounded map (whole-map flushes; exact).
//!
//! `dedup` — merges the shard pattern multisets into the global
//! occurrence stream `occurrences.bin` (pattern + multiplicity, hash
//! partitioned so the resident map stays bounded) plus the merged
//! stats.
//!
//! `cut` — ONE locality's cut of the projected evidence. Loads the
//! locality substrate (panel syng, axis, partition beds, locus.fa,
//! refined bed) and the global syng, builds the committed component
//! territory index (the same construction the committed census used)
//! beside the per-path global<->locality node correspondence (both
//! frames' complete canonical extraction through both syngs over the
//! SAME verified sequence, paired by position), then scans the
//! shards: every record whose remapped node walk touches the
//! locality's territory contributes. Emits the cut census (the
//! committed JSONL schema, routed through the committed placement
//! machinery), the cut derive cache (the committed binary format, so
//! the unchanged CLI binds cut census records to cut cache keys
//! exactly), and the cut reads FASTA. No GRCh38 fetch, no MEM
//! machinery, no re-projection at cut time.
#![recursion_limit = "512"]

use clap::Parser;
use impg::agc_index::AgcIndex;
use impg::genome_inference::mem_records;
use impg::sample_mem_bwt::{canonical, encode_walk};
use impg::syng::{SyncmerParams, SyngIndex};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::Instant;

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

fn ensure(condition: bool, message: &str) -> io::Result<()> {
    if !condition {
        return Err(invalid(message));
    }
    Ok(())
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325;
    for &byte in bytes {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

fn rss_now_kb() -> u64 {
    if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
        for line in status.lines() {
            if let Some(rest) = line.strip_prefix("VmRSS:") {
                return rest
                    .trim()
                    .trim_end_matches("kB")
                    .trim()
                    .parse()
                    .unwrap_or(0);
            }
        }
    }
    0
}

/// A Send-safe stdin reader (`StdinLock` is not `Send`).
struct StdinReader {
    inner: io::Stdin,
}

impl StdinReader {
    fn new() -> Self {
        Self { inner: io::stdin() }
    }
}

impl Read for StdinReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}

#[derive(Parser)]
#[command(about = "the one whole-sample projection and its locality cuts")]
struct Options {
    #[command(subcommand)]
    mode: Mode,
}

#[derive(clap::Subcommand, Debug)]
enum Mode {
    /// The one projection: reads once through the MEM machinery
    /// against the global syng.
    Project {
        /// The GLOBAL panel syng prefix.
        #[arg(long)]
        panel: String,
        /// FASTQ inputs (paths; '-' = stdin). One shard file per input,
        /// named shard.<basename>.bin.gz.
        #[arg(long, required = true, num_args = 1..)]
        inputs: Vec<String>,
        /// Output directory (must not contain the shard files yet).
        #[arg(long)]
        out_dir: PathBuf,
        /// Per-shard pattern-map entry cap (bounded memory; whole-map
        /// flushes keep the multiset exact).
        #[arg(long, default_value_t = 4_000_000)]
        map_cap: usize,
        /// Progress line interval (reads).
        #[arg(long, default_value_t = 20_000_000)]
        progress: u64,
        /// Stop after N reads (smoke mode; 0 = all).
        #[arg(long, default_value_t = 0)]
        limit: u64,
    },
    /// Merge the shard pattern multisets into the occurrence stream.
    Dedup {
        /// The projection output directory (shard.*.occ.bin).
        #[arg(long)]
        out_dir: PathBuf,
        /// Occurrence stream output path.
        #[arg(long)]
        occurrences: PathBuf,
        /// Hash partitions for the bounded merge.
        #[arg(long, default_value_t = 32)]
        partitions: usize,
    },
    /// One locality's cut of the projected evidence.
    Cut {
        /// The GLOBAL panel syng prefix (the projection's panel).
        #[arg(long)]
        global_panel: String,
        /// The LOCALITY panel syng prefix (the committed substrate).
        #[arg(long)]
        panel: String,
        /// The locality axis JSON (the committed substrate).
        #[arg(long)]
        axis: PathBuf,
        /// The locality partition BED directory (the committed substrate).
        #[arg(long)]
        bed_directory: PathBuf,
        /// The component lane (the committed substrate).
        #[arg(long)]
        component: String,
        /// The locality locus.fa (the committed substrate; the territory
        /// source sequences).
        #[arg(long)]
        sources: PathBuf,
        /// The global panel AGC (the committed --sequence-files input).
        #[arg(long)]
        agc: PathBuf,
        /// The projection output directory (the shard files).
        #[arg(long)]
        shards_dir: PathBuf,
        /// Cut output directory (cut-census.jsonl, cut-cache.bin,
        /// cut-reads.fasta, cut-stats.json).
        #[arg(long)]
        out_dir: PathBuf,
    },
}

// ---------------------------------------------------------------------------
// The canonical record derivation (the committed census re-derivation
// and derive-cache builder's own accessor, verbatim semantics).
// ---------------------------------------------------------------------------

type OwnWalk = Vec<(i32, u32)>;

fn derive_read(panel: &SyngIndex, sequence: &[u8]) -> io::Result<Vec<(Vec<u64>, OwnWalk)>> {
    let tagged = mem_records::tagged_mem_records(panel, sequence)?;
    let mut out = Vec::with_capacity(tagged.len());
    for walk in &tagged {
        let encoded = encode_walk(walk)?;
        let tokens = canonical(&encoded);
        out.push((
            tokens,
            walk.iter().map(|&(node, pos)| (node, pos as u32)).collect(),
        ));
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// The shard record stream (little-endian). Each packet:
//   u32 payload_len, then the payload:
//   u32 seq_len, seq bytes, u32 n_records,
//   per record: u32 n_tokens, tokens, u32 n_walk, (i32 node, u32 pos)*
// ---------------------------------------------------------------------------

fn encode_read_packet(seq: &[u8], records: &[(Vec<u64>, OwnWalk)]) -> Vec<u8> {
    let mut payload = Vec::with_capacity(seq.len() + 64);
    payload.extend_from_slice(&(seq.len() as u32).to_le_bytes());
    payload.extend_from_slice(seq);
    payload.extend_from_slice(&(records.len() as u32).to_le_bytes());
    for (tokens, walk) in records {
        payload.extend_from_slice(&(tokens.len() as u32).to_le_bytes());
        for token in tokens {
            payload.extend_from_slice(&token.to_le_bytes());
        }
        payload.extend_from_slice(&(walk.len() as u32).to_le_bytes());
        for &(node, pos) in walk {
            payload.extend_from_slice(&node.to_le_bytes());
            payload.extend_from_slice(&pos.to_le_bytes());
        }
    }
    let mut packet = Vec::with_capacity(payload.len() + 4);
    packet.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    packet.extend_from_slice(&payload);
    packet
}

struct DecodedRead {
    seq: Vec<u8>,
    records: Vec<(Vec<u64>, OwnWalk)>,
}

fn decode_read_packet(packet: &[u8]) -> io::Result<DecodedRead> {
    ensure(packet.len() >= 4, "empty shard packet")?;
    let payload_len = u32::from_le_bytes(packet[..4].try_into().unwrap()) as usize;
    let bytes = packet
        .get(4..4 + payload_len)
        .ok_or_else(|| invalid("truncated shard packet"))?;
    let mut off = 0usize;
    let mut take = |n: usize| -> io::Result<&[u8]> {
        let end = off
            .checked_add(n)
            .filter(|&e| e <= bytes.len())
            .ok_or_else(|| invalid("truncated shard record"))?;
        let slice = &bytes[off..end];
        off = end;
        Ok(slice)
    };
    let seq_len = u32::from_le_bytes(take(4)?.try_into().unwrap()) as usize;
    let seq = take(seq_len)?.to_vec();
    let n_records = u32::from_le_bytes(take(4)?.try_into().unwrap()) as usize;
    let mut records = Vec::with_capacity(n_records);
    for _ in 0..n_records {
        let n_tokens = u32::from_le_bytes(take(4)?.try_into().unwrap()) as usize;
        let mut tokens = Vec::with_capacity(n_tokens);
        for _ in 0..n_tokens {
            tokens.push(u64::from_le_bytes(take(8)?.try_into().unwrap()));
        }
        let n_walk = u32::from_le_bytes(take(4)?.try_into().unwrap()) as usize;
        let mut walk = Vec::with_capacity(n_walk);
        for _ in 0..n_walk {
            let node = i32::from_le_bytes(take(4)?.try_into().unwrap());
            let pos = u32::from_le_bytes(take(4)?.try_into().unwrap());
            walk.push((node, pos));
        }
        records.push((tokens, walk));
    }
    ensure(off == bytes.len(), "trailing bytes in shard record")?;
    Ok(DecodedRead { seq, records })
}

/// Stream packets from one shard file in bounded batches.
fn shard_packet_stream(
    path: &Path,
    batch: usize,
) -> io::Result<impl Iterator<Item = io::Result<Vec<Vec<u8>>>>> {
    let (tx, rx) = crossbeam_channel::bounded::<io::Result<Vec<Vec<u8>>>>(4);
    let path = path.to_path_buf();
    std::thread::spawn(move || {
        let result = (|| -> io::Result<()> {
            let file = std::fs::File::open(&path)?;
            let gz = flate2::read::GzDecoder::new(BufReader::with_capacity(1 << 20, file));
            let mut reader = io::BufReader::with_capacity(1 << 20, gz);
            let mut current: Vec<Vec<u8>> = Vec::with_capacity(batch);
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
                current.push(packet);
                if current.len() == batch {
                    let out = std::mem::replace(&mut current, Vec::with_capacity(batch));
                    if tx.send(Ok(out)).is_err() {
                        return Ok(());
                    }
                }
            }
            if !current.is_empty() {
                let _ = tx.send(Ok(current));
            }
            Ok(())
        })();
        if let Err(e) = result {
            let _ = tx.send(Err(e));
        }
    });
    Ok(rx.into_iter())
}

/// One pattern-multiset entry: u32 n_tokens, tokens, u64 mult.
fn write_occ_entry(out: &mut impl Write, tokens: &[u64], mult: u64) -> io::Result<()> {
    out.write_all(&(tokens.len() as u32).to_le_bytes())?;
    for token in tokens {
        out.write_all(&token.to_le_bytes())?;
    }
    out.write_all(&mult.to_le_bytes())?;
    Ok(())
}

fn read_occ_entry(file: &mut impl Read) -> io::Result<Option<(Vec<u64>, u64)>> {
    let mut buf4 = [0u8; 4];
    if file.read_exact(&mut buf4).is_err() {
        return Ok(None);
    }
    let n = u32::from_le_bytes(buf4) as usize;
    ensure(n > 0 && n % 2 == 1, "invalid occ token count")?;
    let mut tokens = Vec::with_capacity(n);
    let mut buf8 = [0u8; 8];
    for _ in 0..n {
        file.read_exact(&mut buf8)?;
        tokens.push(u64::from_le_bytes(buf8));
    }
    file.read_exact(&mut buf8)?;
    Ok(Some((tokens, u64::from_le_bytes(buf8))))
}

// ---------------------------------------------------------------------------
// project
// ---------------------------------------------------------------------------

enum WriteMsg {
    Packet { bytes: Vec<u8>, tokens: Vec<Vec<u64>> },
    Finish,
}

#[derive(Serialize)]
struct ProjectStatsJson {
    reads: u64,
    bases: u64,
    reads_with_mems: u64,
    records: u64,
    read_lengths: BTreeMap<String, u64>,
    per_shard_reads: BTreeMap<String, u64>,
    pattern_flushes: u64,
    load_seconds: f64,
    wall_seconds: f64,
    peak_rss_kb: u64,
}

struct Totals {
    reads: u64,
    bases: u64,
    reads_with_mems: u64,
    records: u64,
    read_lengths: BTreeMap<usize, u64>,
    per_shard_reads: BTreeMap<String, u64>,
}

fn shard_base_name(input: &str) -> String {
    let trimmed = input.trim_start_matches('-');
    trimmed
        .rsplit('/')
        .next()
        .unwrap_or(trimmed)
        .replace(['.', '/'], "_")
}

#[allow(clippy::too_many_arguments)]
fn process_batch(
    panel: &SyngIndex,
    names: &[String],
    senders: &[crossbeam_channel::Sender<WriteMsg>],
    batch: &mut Vec<(usize, Vec<u8>, Vec<u8>)>,
    totals: &mut Totals,
) -> io::Result<()> {
    if batch.is_empty() {
        return Ok(());
    }
    // The whole per-read pipeline (MEM derivation AND packet encoding)
    // is pure per-read work; only the channel send and the stats are
    // sequential.
    struct Encoded {
        index: usize,
        seq_len: usize,
        reads_with_mems: bool,
        records: usize,
        tokens: Vec<Vec<u64>>,
        bytes: Vec<u8>,
    }
    let results: Vec<io::Result<Encoded>> = batch
        .par_iter()
        .map(|(index, seq, _qual)| {
            let records = derive_read(panel, seq)?;
            let tokens: Vec<Vec<u64>> = records.iter().map(|(t, _)| t.clone()).collect();
            let bytes = encode_read_packet(seq.as_slice(), &records);
            Ok(Encoded {
                index: *index,
                seq_len: seq.len(),
                reads_with_mems: !records.is_empty(),
                records: records.len(),
                tokens,
                bytes,
            })
        })
        .collect();
    batch.clear();
    for result in results {
        let encoded = result?;
        senders[encoded.index]
            .send(WriteMsg::Packet {
                bytes: encoded.bytes,
                tokens: encoded.tokens,
            })
            .map_err(|_| invalid("shard writer exited"))?;
        totals.reads += 1;
        totals.bases += encoded.seq_len as u64;
        if encoded.reads_with_mems {
            totals.reads_with_mems += 1;
        }
        totals.records += encoded.records as u64;
        *totals.read_lengths.entry(encoded.seq_len).or_insert(0) += 1;
        *totals
            .per_shard_reads
            .entry(names[encoded.index].clone())
            .or_insert(0) += 1;
    }
    Ok(())
}

fn run_project(
    panel_prefix: &str,
    inputs: &[String],
    out_dir: &Path,
    map_cap: usize,
    progress: u64,
    limit: u64,
) -> io::Result<()> {
    std::fs::create_dir_all(out_dir)?;
    let names: Vec<String> = inputs.iter().map(|i| shard_base_name(i)).collect();
    for name in &names {
        ensure(
            !out_dir.join(format!("shard.{name}.bin.gz")).exists(),
            "shard output already exists",
        )?;
    }
    let started = Instant::now();
    eprintln!("[project] loading the global syng {panel_prefix}");
    let load_started = Instant::now();
    let panel = SyngIndex::load(panel_prefix, SyncmerParams::default())?;
    let load_seconds = load_started.elapsed().as_secs_f64();
    eprintln!(
        "[project] syng loaded [{load_seconds:.1}s, rss {} kB, {} paths]",
        rss_now_kb(),
        panel.name_map.path_to_name.len()
    );

    // Writer threads: one per shard (gz compression + the pattern map).
    let mut senders: Vec<crossbeam_channel::Sender<WriteMsg>> = Vec::new();
    let mut writer_handles = Vec::new();
    for name in &names {
        let (tx, rx) = crossbeam_channel::bounded::<WriteMsg>(256);
        let shard_path = out_dir.join(format!("shard.{name}.bin.gz"));
        let occ_path = out_dir.join(format!("shard.{name}.occ.bin"));
        writer_handles.push(std::thread::spawn(move || -> io::Result<u64> {
            let mut gz = flate2::write::GzEncoder::new(
                io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&shard_path)?),
                flate2::Compression::fast(),
            );
            let mut occ =
                io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&occ_path)?);
            let mut map: HashMap<Vec<u64>, u64> = HashMap::new();
            let mut flushes = 0u64;
            while let Ok(message) = rx.recv() {
                match message {
                    WriteMsg::Packet { bytes, tokens } => {
                        gz.write_all(&bytes)?;
                        for tokens in tokens {
                            *map.entry(tokens).or_insert(0) += 1;
                        }
                        if map.len() >= map_cap {
                            for (tokens, mult) in map.drain() {
                                write_occ_entry(&mut occ, &tokens, mult)?;
                            }
                            flushes += 1;
                        }
                    }
                    WriteMsg::Finish => break,
                }
            }
            for (tokens, mult) in map.drain() {
                write_occ_entry(&mut occ, &tokens, mult)?;
            }
            occ.flush()?;
            drop(occ);
            gz.finish()?;
            Ok(flushes)
        }));
        senders.push(tx);
    }

    // Reader threads: one per input, streaming 4-line FASTQ records.
    let (read_tx, read_rx) = crossbeam_channel::bounded::<(usize, Vec<u8>, Vec<u8>)>(8192);
    let stop = std::sync::Arc::new(AtomicBool::new(false));
    let mut reader_handles = Vec::new();
    for (index, input) in inputs.iter().enumerate() {
        let tx = read_tx.clone();
        let input = input.clone();
        let stop = stop.clone();
        reader_handles.push(std::thread::spawn(move || -> io::Result<()> {
            let reader: Box<dyn Read + Send> = if input == "-" {
                Box::new(StdinReader::new())
            } else {
                Box::new(std::fs::File::open(&input)?)
            };
            let (mut decoder, _format) =
                niffler::get_reader(Box::new(reader)).map_err(io::Error::other)?;
            let mut lines = BufReader::new(decoder.as_mut()).lines();
            let mut current = match lines.next() {
                Some(line) => line?,
                None => return Ok(()),
            };
            loop {
                if stop.load(Ordering::Relaxed) {
                    break;
                }
                ensure(current.starts_with('@'), "reads input is not FASTQ")?;
                let seq = match lines.next() {
                    Some(line) => line?,
                    None => return Err(invalid("truncated FASTQ sequence")),
                };
                let plus = match lines.next() {
                    Some(line) => line?,
                    None => return Err(invalid("truncated FASTQ separator")),
                };
                let qual = match lines.next() {
                    Some(line) => line?,
                    None => return Err(invalid("truncated FASTQ quality")),
                };
                ensure(
                    plus.starts_with('+') && seq.len() == qual.len(),
                    "invalid FASTQ record",
                )?;
                tx.send((index, seq.into_bytes(), qual.into_bytes()))
                    .map_err(|_| invalid("reader channel closed"))?;
                match lines.next() {
                    Some(line) => current = line?,
                    None => break,
                }
            }
            Ok(())
        }));
    }
    drop(read_tx);

    // The main loop: batches -> parallel MEM derivation -> shard writers.
    let mut totals = Totals {
        reads: 0,
        bases: 0,
        reads_with_mems: 0,
        records: 0,
        read_lengths: BTreeMap::new(),
        per_shard_reads: BTreeMap::new(),
    };
    let peak_rss = AtomicU64::new(rss_now_kb());
    let mut batch: Vec<(usize, Vec<u8>, Vec<u8>)> = Vec::with_capacity(16384);
    let mut received = 0u64;
    while let Ok((index, seq, qual)) = read_rx.recv() {
        received += 1;
        batch.push((index, seq, qual));
        if limit > 0 && received >= limit {
            stop.store(true, Ordering::Relaxed);
        }
        if batch.len() == 16384 {
            process_batch(&panel, &names, &senders, &mut batch, &mut totals)?;
        }
        if received % progress == 0 {
            let rss = rss_now_kb();
            peak_rss.fetch_max(rss, Ordering::Relaxed);
            eprintln!(
                "[project] reads {} [{:.1}s, rss {rss} kB]",
                totals.reads,
                started.elapsed().as_secs_f64()
            );
        }
    }
    process_batch(&panel, &names, &senders, &mut batch, &mut totals)?;
    for sender in &senders {
        let _ = sender.send(WriteMsg::Finish);
    }
    drop(senders);
    let mut pattern_flushes = 0u64;
    for handle in writer_handles {
        pattern_flushes += handle.join().map_err(|_| invalid("writer panicked"))??;
    }
    for handle in reader_handles {
        handle
            .join()
            .map_err(|_| invalid("reader panicked"))??;
    }
    let wall_seconds = started.elapsed().as_secs_f64();
    peak_rss.fetch_max(rss_now_kb(), Ordering::Relaxed);
    let stats = ProjectStatsJson {
        reads: totals.reads,
        bases: totals.bases,
        reads_with_mems: totals.reads_with_mems,
        records: totals.records,
        read_lengths: totals
            .read_lengths
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect(),
        per_shard_reads: totals.per_shard_reads,
        pattern_flushes,
        load_seconds,
        wall_seconds,
        peak_rss_kb: peak_rss.load(Ordering::Relaxed),
    };
    let text = serde_json::to_string_pretty(&stats)?;
    std::fs::write(out_dir.join("projection-stats.json"), text + "\n")?;
    eprintln!(
        "[project] DONE reads {} bases {} records {} flushes {pattern_flushes} \
         [{wall_seconds:.1}s, peak rss {} kB]",
        stats.reads, stats.bases, stats.records, stats.peak_rss_kb
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// dedup
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct DedupStatsJson {
    distinct_patterns: u64,
    total_mult: u64,
    per_shard: BTreeMap<String, u64>,
    wall_seconds: f64,
    peak_rss_kb: u64,
}

fn run_dedup(out_dir: &Path, occurrences: &Path, partitions: usize) -> io::Result<()> {
    let started = Instant::now();
    ensure(partitions > 0, "partitions must be positive")?;
    let mut occ_paths: Vec<PathBuf> = Vec::new();
    for entry in std::fs::read_dir(out_dir)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name.starts_with("shard.") && name.ends_with(".occ.bin") {
            occ_paths.push(path);
        }
    }
    occ_paths.sort();
    ensure(!occ_paths.is_empty(), "no shard occ files found")?;
    let mut bins: Vec<io::BufWriter<std::fs::File>> = Vec::with_capacity(partitions);
    for p in 0..partitions {
        let path = occurrences.with_extension(format!("bin{p}"));
        std::fs::remove_file(&path).ok();
        bins.push(io::BufWriter::with_capacity(
            1 << 20,
            std::fs::File::create(&path)?,
        ));
    }
    let mut per_shard: BTreeMap<String, u64> = BTreeMap::new();
    for path in &occ_paths {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        let mut entries = 0u64;
        let file = std::fs::File::open(path)?;
        let mut reader = io::BufReader::with_capacity(1 << 20, file);
        while let Some((tokens, mult)) = read_occ_entry(&mut reader)? {
            let bin = (fnv1a64(
                &tokens
                    .iter()
                    .flat_map(|t| t.to_le_bytes())
                    .collect::<Vec<_>>(),
            ) as usize)
                % partitions;
            write_occ_entry(&mut bins[bin], &tokens, mult)?;
            entries += 1;
        }
        per_shard.insert(name, entries);
    }
    for mut bin in bins {
        bin.flush()?;
    }
    let mut out = io::BufWriter::with_capacity(1 << 20, std::fs::File::create(occurrences)?);
    let mut distinct = 0u64;
    let mut total_mult = 0u64;
    for p in 0..partitions {
        let path = occurrences.with_extension(format!("bin{p}"));
        let file = std::fs::File::open(&path)?;
        let mut reader = io::BufReader::with_capacity(1 << 20, file);
        let mut map: HashMap<Vec<u64>, u64> = HashMap::new();
        while let Some((tokens, mult)) = read_occ_entry(&mut reader)? {
            *map.entry(tokens).or_insert(0) += mult;
        }
        let sorted: BTreeMap<Vec<u64>, u64> = map.into_iter().collect();
        for (tokens, mult) in sorted {
            write_occ_entry(&mut out, &tokens, mult)?;
            distinct += 1;
            total_mult += mult;
        }
        std::fs::remove_file(&path)?;
    }
    out.flush()?;
    let stats = DedupStatsJson {
        distinct_patterns: distinct,
        total_mult,
        per_shard,
        wall_seconds: started.elapsed().as_secs_f64(),
        peak_rss_kb: rss_now_kb(),
    };
    let text = serde_json::to_string_pretty(&stats)?;
    std::fs::write(
        Path::new(occurrences).with_extension("stats.json"),
        text + "\n",
    )?;
    eprintln!(
        "[dedup] DONE distinct {distinct} total-mult {total_mult} [{:.1}s]",
        stats.wall_seconds
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// cut: the committed census machinery mirrors (verbatim semantics).
// ---------------------------------------------------------------------------

const SLACK: u64 = 150;

#[derive(serde::Deserialize)]
struct AxisFile {
    intervals: Vec<AxisInterval>,
}

#[derive(serde::Deserialize)]
struct AxisInterval {
    component: String,
    group: String,
}

struct UniverseRow {
    partition: u32,
    path_name: String,
    start: u64,
    end: u64,
}

struct Territory {
    partition: u32,
    path_idx: usize,
    start: u64,
    end: u64,
    steps: Vec<(u64, i32)>,
}

struct TerritoryIndex {
    territories: Vec<Territory>,
    entries: Vec<(i32, u32, u64)>,
    path_intervals: Vec<Vec<(u64, u64, u32)>>,
}

impl TerritoryIndex {
    fn node_range(&self, node: i32) -> (usize, usize) {
        let start = self.entries.partition_point(|&(n, _, _)| n < node);
        let end = self.entries.partition_point(|&(n, _, _)| n <= node);
        (start, end)
    }

    fn partitions_at(&self, path_idx: usize, bp: u64, k: u64) -> Vec<u32> {
        let mut out = Vec::new();
        for &(start, end, partition) in &self.path_intervals[path_idx] {
            if bp + k > start && bp < end {
                out.push(partition);
            }
        }
        out
    }
}

fn rc_frame_step(signed_hash: i32, q: u64, range_lo: u64, range_len: u64, k: u64) -> (u64, i32) {
    (range_lo + range_len - k - q, -signed_hash)
}

fn decode_tokens(tokens: &[u64]) -> io::Result<Vec<(i32, u64)>> {
    ensure(!tokens.is_empty() && tokens.len() % 2 == 1, "invalid record tokens")?;
    let mut anchors = Vec::with_capacity(tokens.len() / 2 + 1);
    let mut position = 0u64;
    for (i, &token) in tokens.iter().enumerate() {
        if i % 2 == 0 {
            let zigzag = token
                .checked_sub(2)
                .ok_or_else(|| invalid("record node token"))?
                / 2;
            let node = ((zigzag >> 1) as i64 ^ -(zigzag as i64 & 1)) as i32;
            anchors.push((node, position));
        } else {
            let gap = token
                .checked_sub(1)
                .ok_or_else(|| invalid("record gap token"))?
                / 2;
            position += gap;
        }
    }
    Ok(anchors)
}

fn reverse_complement_walk(anchors: &[(i32, u64)], k: u64) -> Vec<(i32, u64)> {
    let span = anchors.last().map(|&(_, p)| p + k).unwrap_or(k);
    let mut out = Vec::with_capacity(anchors.len());
    for &(node, pos) in anchors.iter().rev() {
        out.push((-node, span - k - pos));
    }
    out
}

fn verify_walk(walk: &[(i32, u64)], start: u64, territory: &Territory) -> bool {
    let steps = &territory.steps;
    for &(node, rel) in walk {
        let bp = match start.checked_add(rel) {
            Some(bp) => bp,
            None => return false,
        };
        let mut lo = 0usize;
        let mut hi = steps.len();
        let mut found = false;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            match steps[mid].0.cmp(&bp) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    found = steps[mid].1 == node;
                    break;
                }
            }
        }
        if !found {
            return false;
        }
    }
    true
}

fn occurrences_from_anchor(
    walk: &[(i32, u64)],
    anchor: usize,
    index: &TerritoryIndex,
    seen: &mut HashSet<(usize, u64)>,
) {
    let (start, end) = index.node_range(walk[anchor].0);
    let base = walk[anchor].1;
    for entry in &index.entries[start..end] {
        let territory_index = entry.1 as usize;
        let occurrence_start = entry.2.saturating_sub(base);
        let territory = &index.territories[territory_index];
        let key = (territory.path_idx, occurrence_start);
        if seen.contains(&key) {
            continue;
        }
        if verify_walk(walk, occurrence_start, territory) {
            seen.insert(key);
        }
    }
}

#[allow(clippy::type_complexity)]
fn route_record(
    anchors: &[(i32, u64)],
    index: &TerritoryIndex,
    k: u64,
) -> (
    BTreeMap<u32, u64>,
    u64,
    Vec<(usize, u64)>,
    Vec<(usize, u64)>,
) {
    let mut occurrences: BTreeMap<u32, u64> = BTreeMap::new();
    let mut total_occurrences = 0u64;
    let mut forward_positions: Vec<(usize, u64)> = Vec::new();
    let mut reverse_positions: Vec<(usize, u64)> = Vec::new();
    for (orientation_index, orientation) in [anchors, &reverse_complement_walk(anchors, k)]
        .into_iter()
        .enumerate()
    {
        let mut seen = HashSet::new();
        let mut best = 0usize;
        let mut best_entries = usize::MAX;
        for (anchor, &(node, _)) in orientation.iter().enumerate() {
            let (start, end) = index.node_range(node);
            let count = end - start;
            if count < best_entries {
                best_entries = count;
                best = anchor;
            }
        }
        occurrences_from_anchor(orientation, best, index, &mut seen);
        total_occurrences += seen.len() as u64;
        if orientation_index == 0 {
            forward_positions = seen.iter().copied().collect();
            forward_positions.sort_unstable();
        } else {
            reverse_positions = seen.iter().copied().collect();
            reverse_positions.sort_unstable();
        }
        for (path_idx, occurrence_start) in seen {
            let mut touched = BTreeSet::new();
            for &(_node, rel) in orientation {
                let bp = occurrence_start + rel;
                for partition in index.partitions_at(path_idx, bp, k) {
                    touched.insert(partition);
                }
            }
            for partition in touched {
                *occurrences.entry(partition).or_default() += 1;
            }
        }
    }
    (
        occurrences,
        total_occurrences,
        forward_positions,
        reverse_positions,
    )
}

const MULTI_CENSUS_PAIR_BINS: [&str; 8] = [
    "coalesced_same_window",
    "coalesced_adjacent_window",
    "coalesced_same_contig_far",
    "coalesced_other_contig",
    "same_window",
    "adjacent_window",
    "same_contig_far",
    "other_contig",
];

fn multi_census_pair_bin(
    nodes_a: &[u32],
    nodes_b: &[u32],
    partitions_a: &[u32],
    partitions_b: &[u32],
    same_component_contig: bool,
) -> usize {
    let coalesced = nodes_a.iter().any(|&node| nodes_b.contains(&node));
    let physical = if partitions_a
        .iter()
        .any(|&partition| partitions_b.contains(&partition))
    {
        0
    } else if partitions_a
        .iter()
        .any(|&left| partitions_b.iter().any(|&right| left.abs_diff(right) == 1))
    {
        1
    } else if same_component_contig {
        2
    } else {
        3
    };
    if coalesced {
        physical
    } else {
        physical + 4
    }
}

fn load_universe_component(
    axis: &AxisFile,
    bed_directory: &Path,
    component: &str,
) -> io::Result<Vec<UniverseRow>> {
    let mut rows = Vec::new();
    for (partition, interval) in axis
        .intervals
        .iter()
        .filter(|i| i.component == component)
        .enumerate()
    {
        let bed = bed_directory.join(format!("{}.bed", interval.group));
        for line in BufReader::new(std::fs::File::open(&bed)?).lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            ensure(fields.len() == 3, "invalid public BED3 group")?;
            let start: u64 = fields[1].parse().map_err(|_| invalid("BED start"))?;
            let end: u64 = fields[2].parse().map_err(|_| invalid("BED end"))?;
            ensure(start < end, "BED interval")?;
            rows.push(UniverseRow {
                partition: partition as u32,
                path_name: fields[0].to_string(),
                start,
                end,
            });
        }
    }
    ensure(!rows.is_empty(), "component universe is empty")?;
    Ok(rows)
}

/// The committed territory step construction (the canonical scheme),
/// verbatim semantics: forward steps from the panel's own stored path
/// walk (sampled sidecars exactly as the committed census saw them),
/// the rc frame's steps from the complete raw extraction, and the
/// per-position canonical-frame selection.
fn build_territory_index(
    panel: &SyngIndex,
    rows: &[UniverseRow],
    sources: &HashMap<String, Vec<u8>>,
    k: u64,
) -> io::Result<TerritoryIndex> {
    let path_of_name: HashMap<&str, usize> = panel
        .name_map
        .path_to_name
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();
    let path_count = panel.name_map.path_to_name.len();
    let syncmer_len = k;
    let territories: Vec<Territory> = rows
        .par_iter()
        .map(|row| {
            let path_idx = *path_of_name
                .get(row.path_name.as_str())
                .ok_or_else(|| invalid(&format!("universe path {} absent", row.path_name)))?;
            let lo = row.start.saturating_sub(SLACK);
            let hi = row.end + SLACK;
            let forward_steps: Vec<(u64, i32)> = panel
                .walk_path_range(path_idx, lo, hi)?
                .into_iter()
                .map(|(node, bp)| (bp, node))
                .collect();
            let seq_lo = lo.saturating_sub(syncmer_len);
            let seq_hi = hi + syncmer_len;
            let lane = sources
                .get(&row.path_name)
                .ok_or_else(|| invalid(&format!("source {} absent", row.path_name)))?;
            let hi_clip = seq_hi.min(lane.len() as u64);
            let seq: &[u8] = if seq_lo >= hi_clip {
                &[]
            } else {
                &lane[seq_lo as usize..hi_clip as usize]
            };
            let rc_seq = impg::graph::reverse_complement(seq);
            let reverse_steps: Vec<(u64, i32)> = mem_records::raw_matched_syncmers(panel, &rc_seq)?
                .into_iter()
                .map(|(signed_node, q)| {
                    rc_frame_step(signed_node, q, seq_lo, seq.len() as u64, syncmer_len)
                })
                .collect();
            let forward_map: BTreeMap<u64, i32> = forward_steps.into_iter().collect();
            let reverse_map: BTreeMap<u64, i32> = reverse_steps.into_iter().collect();
            let mut positions: BTreeSet<u64> = forward_map.keys().copied().collect();
            positions.extend(reverse_map.keys().copied());
            let mut steps = Vec::with_capacity(positions.len());
            for bp in positions {
                let start = (bp.saturating_sub(seq_lo)) as usize;
                let canonical_forward = seq
                    .get(start..start + syncmer_len as usize)
                    .map(mem_records::window_is_canonical_forward)
                    .unwrap_or(false);
                let chosen = if canonical_forward {
                    forward_map.get(&bp)
                } else {
                    reverse_map.get(&bp)
                };
                if let Some(&node) = chosen {
                    steps.push((bp, node));
                }
            }
            Ok::<_, io::Error>(Territory {
                partition: row.partition,
                path_idx,
                start: row.start,
                end: row.end,
                steps,
            })
        })
        .collect::<io::Result<_>>()?;
    let mut path_intervals: Vec<Vec<(u64, u64, u32)>> = vec![Vec::new(); path_count];
    for t in &territories {
        path_intervals[t.path_idx].push((t.start, t.end, t.partition));
    }
    for intervals in &mut path_intervals {
        intervals.sort_by_key(|(start, end, _)| (*start, *end));
    }
    let mut entries = Vec::with_capacity(territories.len() * 512);
    for (index, t) in territories.iter().enumerate() {
        for &(bp, node) in &t.steps {
            if bp + k > t.start && bp < t.end {
                entries.push((node, index as u32, bp));
            }
        }
    }
    entries.par_sort_unstable_by_key(|&(node, _, _)| node);
    Ok(TerritoryIndex {
        territories,
        entries,
        path_intervals,
    })
}

/// The derive cache writer (the committed binary format, verbatim).
fn write_derive_cache(
    path: &Path,
    key_tokens: &[Vec<u64>],
    key_reads: &[Vec<usize>],
    reads: &[Vec<u8>],
    read_multiplicity: &[u32],
    read_records: &[Vec<(u32, OwnWalk)>],
) -> io::Result<()> {
    let mut out = io::BufWriter::with_capacity(1 << 20, std::fs::File::create(path)?);
    out.write_all(&(key_tokens.len() as u64).to_le_bytes())?;
    for tokens in key_tokens {
        out.write_all(&(tokens.len() as u32).to_le_bytes())?;
        for token in tokens {
            out.write_all(&token.to_le_bytes())?;
        }
    }
    out.write_all(&(key_reads.len() as u64).to_le_bytes())?;
    for readers in key_reads {
        out.write_all(&(readers.len() as u32).to_le_bytes())?;
        for reader in readers {
            out.write_all(&(*reader as u32).to_le_bytes())?;
        }
    }
    out.write_all(&(reads.len() as u64).to_le_bytes())?;
    for seq in reads {
        out.write_all(&(seq.len() as u32).to_le_bytes())?;
        out.write_all(seq)?;
    }
    for &mult in read_multiplicity {
        out.write_all(&mult.to_le_bytes())?;
    }
    out.write_all(&(read_records.len() as u64).to_le_bytes())?;
    for records in read_records {
        out.write_all(&(records.len() as u32).to_le_bytes())?;
        for (key, walk) in records {
            out.write_all(&key.to_le_bytes())?;
            out.write_all(&(walk.len() as u32).to_le_bytes())?;
            for &(node, pos) in walk {
                out.write_all(&node.to_le_bytes())?;
                out.write_all(&pos.to_le_bytes())?;
            }
        }
    }
    out.flush()?;
    Ok(())
}

#[derive(Serialize)]
struct CutStatsJson {
    reads_scanned: u64,
    reads_kept: u64,
    records_scanned: u64,
    records_mappable: u64,
    records_unmappable_global_only: u64,
    reads_with_unmappable_records: u64,
    touching_patterns: u64,
    routed_patterns: u64,
    census_mult: u64,
    territory_rows: usize,
    territory_entries: usize,
    node_map_pairs: usize,
    offset_variants: BTreeMap<String, u64>,
    load_seconds: f64,
    scan_seconds: f64,
    census_seconds: f64,
    wall_seconds: f64,
    peak_rss_kb: u64,
}

/// One shard file's cut contribution.
struct FileCut {
    reads_scanned: u64,
    records_scanned: u64,
    records_mappable: u64,
    records_unmappable: u64,
    reads_with_unmappable: u64,
    reads: Vec<Vec<u8>>,
    pattern_mult: BTreeMap<Vec<u64>, u64>,
    key_tokens: Vec<Vec<u64>>,
    key_index: HashMap<Vec<u64>, u32>,
    key_reads: Vec<Vec<usize>>,
    /// per kept read (index): its (local key, remapped own walk) list
    read_records: Vec<Vec<(u32, OwnWalk)>>,
}

impl FileCut {
    fn intern_key(&mut self, tokens: &[u64]) -> u32 {
        if let Some(&index) = self.key_index.get(tokens) {
            return index;
        }
        let index = self.key_tokens.len() as u32;
        self.key_tokens.push(tokens.to_vec());
        self.key_reads.push(Vec::new());
        self.key_index.insert(tokens.to_vec(), index);
        index
    }
}

fn scan_shard_file(
    path: &Path,
    node_map: &HashMap<i32, i32>,
    territory_nodes_abs: &HashSet<u32>,
) -> io::Result<FileCut> {
    let mut cut = FileCut {
        reads_scanned: 0,
        records_scanned: 0,
        records_mappable: 0,
        records_unmappable: 0,
        reads_with_unmappable: 0,
        reads: Vec::new(),
        pattern_mult: BTreeMap::new(),
        key_tokens: Vec::new(),
        key_index: HashMap::new(),
        key_reads: Vec::new(),
        read_records: Vec::new(),
    };
    for batch in shard_packet_stream(path, 1024)? {
        let packets = batch?;
        cut.reads_scanned += packets.len() as u64;
        let per_read: Vec<io::Result<DecodedRead>> = packets
            .par_iter()
            .map(|packet| decode_read_packet(packet))
            .collect();
        let decoded: Vec<DecodedRead> = per_read
            .into_iter()
            .collect::<io::Result<_>>()?;
        let remapped: Vec<io::Result<Vec<Option<(Vec<u64>, OwnWalk, bool)>>>> = decoded
            .par_iter()
            .map(|read| {
                let mut out = Vec::with_capacity(read.records.len());
                cut_records(
                    &read.records,
                    node_map,
                    territory_nodes_abs,
                    &mut out,
                )?;
                Ok(out)
            })
            .collect();
        for (read, result) in decoded.into_iter().zip(remapped) {
            let records = result?;
            cut.records_scanned += records.len() as u64;
            let mut kept: Vec<(u32, OwnWalk)> = Vec::new();
            let mut any_touch = false;
            let mut had_unmappable = false;
            for record in &records {
                match record {
                    None => {
                        cut.records_unmappable += 1;
                        had_unmappable = true;
                    }
                    Some((tokens, walk, touch)) => {
                        cut.records_mappable += 1;
                        if *touch {
                            any_touch = true;
                            *cut.pattern_mult.entry(tokens.clone()).or_insert(0) += 1;
                        }
                        let key = cut.intern_key(tokens);
                        kept.push((key, walk.clone()));
                    }
                }
            }
            if had_unmappable {
                cut.reads_with_unmappable += 1;
            }
            if any_touch {
                let read_index = cut.reads.len();
                for (key, _walk) in &kept {
                    cut.key_reads[*key as usize].push(read_index);
                }
                cut.reads.push(read.seq);
                cut.read_records.push(std::mem::take(&mut kept));
            }
        }
    }
    Ok(cut)
}

/// Remap one read's records into locality node ids; None marks a record
/// with at least one node absent from the locality panel (the global
/// superset class — not expressible in the locality's node space).
fn cut_records(
    records: &[(Vec<u64>, OwnWalk)],
    node_map: &HashMap<i32, i32>,
    territory_nodes_abs: &HashSet<u32>,
    out: &mut Vec<Option<(Vec<u64>, OwnWalk, bool)>>,
) -> io::Result<()> {
    for (tokens, walk) in records {
        let anchors = decode_tokens(tokens)?;
        let mut mapped = Vec::with_capacity(anchors.len());
        let mut all_mapped = true;
        for &(node, pos) in &anchors {
            match node_map.get(&node) {
                Some(&local) => mapped.push((local, pos)),
                None => {
                    all_mapped = false;
                    break;
                }
            }
        }
        if !all_mapped {
            out.push(None);
            continue;
        }
        let new_tokens = canonical(&encode_walk(&mapped)?);
        let mut walk_mapped = Vec::with_capacity(walk.len());
        let mut walk_ok = true;
        for &(node, pos) in walk {
            match node_map.get(&node) {
                Some(&local) => walk_mapped.push((local, pos)),
                None => {
                    walk_ok = false;
                    break;
                }
            }
        }
        if !walk_ok {
            out.push(None);
            continue;
        }
        let touch = mapped
            .iter()
            .any(|&(node, _)| territory_nodes_abs.contains(&node.unsigned_abs()));
        out.push(Some((new_tokens, walk_mapped, touch)));
    }
    Ok(())
}

/// Parse a locus.fa record header of the form "path:start-end" (the
/// substrate's own provenance statement: the global panel path and the
/// interval the record was extracted from).
fn parse_locus_header(header: &str) -> io::Result<(String, u64, u64)> {
    let (name, range) = header
        .rsplit_once(':')
        .ok_or_else(|| invalid(&format!("locus header {header} carries no interval")))?;
    let (start, end) = range
        .split_once('-')
        .ok_or_else(|| invalid(&format!("locus header {header} carries no interval range")))?;
    let start: u64 = start
        .parse()
        .map_err(|_| invalid(&format!("locus header {header} start")))?;
    let end: u64 = end
        .parse()
        .map_err(|_| invalid(&format!("locus header {header} end")))?;
    Ok((name.to_string(), start, end))
}

#[allow(clippy::too_many_arguments)]
fn run_cut(
    global_panel_prefix: &str,
    panel_prefix: &str,
    axis_path: &Path,
    bed_directory: &Path,
    component: &str,
    sources_path: &Path,
    agc_path: &Path,
    shards_dir: &Path,
    out_dir: &Path,
) -> io::Result<()> {
    std::fs::create_dir_all(out_dir)?;
    for name in [
        "cut-census.jsonl",
        "cut-cache.bin",
        "cut-reads.fasta",
        "cut-stats.json",
    ] {
        ensure(
            !out_dir.join(name).exists(),
            &format!("cut output {name} already exists"),
        )?;
    }
    let started = Instant::now();
    eprintln!("[cut] loading panels");
    let load_started = Instant::now();
    let panel = SyngIndex::load(panel_prefix, SyncmerParams::default())?;
    let global = SyngIndex::load(global_panel_prefix, SyncmerParams::default())?;
    let agc = AgcIndex::build_from_files(&[agc_path.display().to_string()])?;
    let k = panel.syncmer_length_bp() as u64;
    ensure(
        global.syncmer_length_bp() as u64 == k,
        "global and locality syncmer lengths differ",
    )?;

    let axis: AxisFile =
        serde_json::from_str(&std::fs::read_to_string(axis_path)?).map_err(|e| invalid(format!("axis parse: {e}")))?;
    let universe = load_universe_component(&axis, bed_directory, component)?;
    eprintln!("[cut] universe rows: {}", universe.len());

    // The locus sources (the committed substrate's sequences).
    let mut sources: HashMap<String, Vec<u8>> = HashMap::new();
    {
        let text = std::fs::read_to_string(sources_path)?;
        let mut name: Option<String> = None;
        let mut seq: Vec<u8> = Vec::new();
        for line in text.lines() {
            if let Some(header) = line.strip_prefix('>') {
                if let Some(prev) = name.take() {
                    sources.insert(prev, std::mem::take(&mut seq));
                }
                name = Some(header.split_whitespace().next().unwrap_or("").to_string());
            } else if name.is_some() {
                seq.extend_from_slice(line.trim().as_bytes());
            }
        }
        if let Some(prev) = name {
            sources.insert(prev, seq);
        }
    }
    ensure(!sources.is_empty(), "no locus.fa records")?;

    let territory = build_territory_index(&panel, &universe, &sources, k)?;
    eprintln!(
        "[cut] territory: {} rows, {} entries",
        territory.territories.len(),
        territory.entries.len()
    );
    let territory_nodes_abs: HashSet<u32> = territory
        .entries
        .iter()
        .map(|&(node, _, _)| node.unsigned_abs())
        .collect();

    // The per-path node correspondence: the complete canonical walk of
    // the SAME AGC-verified sequence through both syngs, paired by
    // position. Each locus.fa header states its own provenance — the
    // global panel path and interval "path:start-end" (the refined-bed
    // rows and the axis-window records alike) — so the correspondence
    // is derived from the substrate's own statement, with the offset
    // convention (exact or one-base-earlier) verified by sequence
    // identity against the AGC.
    let mut node_map: HashMap<i32, i32> = HashMap::new();
    let mut node_map_pairs = 0usize;
    let mut global_only_anchor_total = 0usize;
    let mut offset_variants: BTreeMap<String, u64> = BTreeMap::new();
    for (header, record) in &sources {
        let (global_name, s, e) = parse_locus_header(header)?;
        ensure(
            (e - s) as usize == record.len(),
            &format!("locus record {header} length differs from its header interval"),
        )?;
        let mut delta: Option<i64> = None;
        for candidate in [0i64, -1] {
            let lo = (s as i64 + candidate).max(0) as usize;
            let hi = lo + record.len();
            if let Ok(fetched) = agc.fetch_sequence(&global_name, lo as i32, hi as i32) {
                if fetched == *record {
                    delta = Some(candidate);
                    break;
                }
            }
        }
        let delta = delta.ok_or_else(|| {
            invalid(format!(
                "locus record {header} is not identical to the AGC interval \
                 under either offset convention"
            ))
        })?;
        *offset_variants
            .entry(format!("delta={delta}"))
            .or_insert(0) += 1;
        let lo = (s as i64 + delta).max(0) as i32;
        let global_seq = agc
            .fetch_sequence(&global_name, lo, lo + record.len() as i32)
            .map_err(|e| invalid(format!("AGC fetch {global_name}: {e}")))?;
        ensure(global_seq == *record, "AGC refetch mismatch")?;
        let local_walk = mem_records::canonical_anchor_walk(&panel, record)?;
        let global_walk = mem_records::canonical_anchor_walk(&global, &global_seq)?;
        // The locality dictionary is a subset of the global one, so every
        // locality anchor position must appear in the global walk; the
        // global walk may carry EXTRA anchors (the superset class). Pair
        // on the locality positions (the true correspondence); the
        // global-only positions have no locality node and stay unmapped.
        let global_by_pos: HashMap<u64, i32> =
            global_walk.iter().map(|&(node, pos)| (pos, node)).collect();
        for &(lnode, lpos) in &local_walk {
            let gnode = *global_by_pos
                .get(&lpos)
                .ok_or_else(|| invalid("locality anchor position absent from the global walk"))?;
            match node_map.entry(gnode) {
                std::collections::hash_map::Entry::Occupied(o) => {
                    ensure(
                        *o.get() == lnode,
                        "conflicting node correspondence",
                    )?;
                }
                std::collections::hash_map::Entry::Vacant(v) => {
                    v.insert(lnode);
                    node_map_pairs += 1;
                }
            }
        }
        global_only_anchor_total += global_walk.len() - local_walk.len();
    }
    eprintln!(
        "[cut] node map: {node_map_pairs} pairs, {global_only_anchor_total} global-only anchors ({:?})",
        offset_variants
    );
    let load_seconds = load_started.elapsed().as_secs_f64();

    // The shard scan (parallel over files; deterministic merge order).
    let mut shard_paths: Vec<PathBuf> = Vec::new();
    for entry in std::fs::read_dir(shards_dir)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name.starts_with("shard.") && name.ends_with(".bin.gz") {
            shard_paths.push(path);
        }
    }
    shard_paths.sort();
    ensure(!shard_paths.is_empty(), "no shard files found")?;
    let node_map_ref = &node_map;
    let territory_nodes_ref = &territory_nodes_abs;
    let per_file: Vec<io::Result<FileCut>> = shard_paths
        .par_iter()
        .map(|path| scan_shard_file(path, node_map_ref, territory_nodes_ref))
        .collect();

    // Merge in file order (deterministic).
    struct Merged {
        reads: Vec<Vec<u8>>,
        pattern_mult: BTreeMap<Vec<u64>, u64>,
        key_tokens: Vec<Vec<u64>>,
        key_index: HashMap<Vec<u64>, u32>,
        key_reads: Vec<Vec<usize>>,
        read_records: Vec<Vec<(u32, OwnWalk)>>,
        reads_scanned: u64,
        records_scanned: u64,
        records_mappable: u64,
        records_unmappable: u64,
        reads_with_unmappable: u64,
    }
    let mut merged = Merged {
        reads: Vec::new(),
        pattern_mult: BTreeMap::new(),
        key_tokens: Vec::new(),
        key_index: HashMap::new(),
        key_reads: Vec::new(),
        read_records: Vec::new(),
        reads_scanned: 0,
        records_scanned: 0,
        records_mappable: 0,
        records_unmappable: 0,
        reads_with_unmappable: 0,
    };
    fn intern_merged(merged: &mut Merged, tokens: &[u64]) -> u32 {
        if let Some(&index) = merged.key_index.get(tokens) {
            return index;
        }
        let index = merged.key_tokens.len() as u32;
        merged.key_tokens.push(tokens.to_vec());
        merged.key_reads.push(Vec::new());
        merged.key_index.insert(tokens.to_vec(), index);
        index
    }
    for result in per_file {
        let cut = result?;
        merged.reads_scanned += cut.reads_scanned;
        merged.records_scanned += cut.records_scanned;
        merged.records_mappable += cut.records_mappable;
        merged.records_unmappable += cut.records_unmappable;
        merged.reads_with_unmappable += cut.reads_with_unmappable;
        let base = merged.reads.len();
        for (index, seq) in cut.reads.into_iter().enumerate() {
            merged.reads.push(seq);
            let read_index = base + index;
            let records = cut.read_records[index].clone();
            let mut rebased: Vec<(u32, OwnWalk)> = Vec::with_capacity(records.len());
            for (local_key, walk) in records {
                let tokens = cut.key_tokens[local_key as usize].clone();
                let new_key = intern_merged(&mut merged, &tokens);
                merged.key_reads[new_key as usize].push(read_index);
                rebased.push((new_key, walk));
            }
            merged.read_records.push(rebased);
        }
        for (tokens, mult) in cut.pattern_mult {
            *merged.pattern_mult.entry(tokens).or_insert(0) += mult;
        }
    }
    eprintln!(
        "[cut] scan done: {} reads scanned, {} kept, {} touching patterns",
        merged.reads_scanned,
        merged.reads.len(),
        merged.pattern_mult.len()
    );
    let scan_seconds = started.elapsed().as_secs_f64() - load_seconds;

    // The derive cache (the committed format): multiplicity by identical
    // sequence, first representative carries the count.
    let reads = merged.reads;
    let mut read_counts: HashMap<u64, u32> = HashMap::new();
    for seq in &reads {
        *read_counts.entry(fnv1a64(seq)).or_insert(0) += 1;
    }
    let mut read_multiplicity: Vec<u32> = Vec::with_capacity(reads.len());
    {
        let mut first_seen: HashMap<u64, usize> = HashMap::new();
        for seq in &reads {
            let hash = fnv1a64(seq);
            match first_seen.get(&hash) {
                Some(_) => read_multiplicity.push(0),
                None => {
                    first_seen.insert(hash, 0);
                    read_multiplicity.push(*read_counts.get(&hash).unwrap_or(&1));
                }
            }
        }
    }
    write_derive_cache(
        &out_dir.join("cut-cache.bin"),
        &merged.key_tokens,
        &merged.key_reads,
        &reads,
        &read_multiplicity,
        &merged.read_records,
    )?;
    {
        let mut out =
            io::BufWriter::new(std::fs::File::create(out_dir.join("cut-reads.fasta"))?);
        for (index, seq) in reads.iter().enumerate() {
            writeln!(out, ">cutread:{index}")?;
            out.write_all(seq)?;
            out.write_all(b"\n")?;
        }
        out.flush()?;
    }

    // The cut census: route every touching pattern through the committed
    // placement machinery; emit the committed JSONL schema.
    let census_started = Instant::now();
    let mut out = io::BufWriter::new(std::fs::File::create(out_dir.join("cut-census.jsonl"))?);
    let component_contig = component.rsplit('#').next().unwrap_or(component).to_string();
    let mut routed_patterns = 0u64;
    let mut census_mult_total = 0u64;
    let mut dropped_no_occurrence = 0u64;
    // The committed census record order: the pattern multiset's BTreeMap
    // order (sorted by tokens), dense record ids over the routed set.
    // The routing and per-occurrence covered-node construction are pure
    // per-pattern functions over shared immutable state, so they run in
    // parallel and the lines are written in the sorted order.
    let sorted: Vec<(Vec<u64>, u64)> = merged
        .pattern_mult
        .clone()
        .into_iter()
        .collect::<BTreeMap<_, _>>()
        .into_iter()
        .collect();
    let routed: Vec<io::Result<Option<(serde_json::Value, u64)>>> = sorted
        .par_iter()
        .map(|(tokens, multiplicity)| {
            let anchors = decode_tokens(tokens)?;
            let (occurrence_map, total_occurrences, forward_positions, reverse_positions) =
                route_record(&anchors, &territory, k);
            if occurrence_map.is_empty() {
                return Ok(None);
            }
            let anchor_count = anchors.len();
            struct CensusOcc {
                path: usize,
                start: u64,
                orientation: u8,
                intervals: Vec<(u64, u64)>,
                partitions: Vec<u32>,
            }
            let mut occurrences: Vec<CensusOcc> = Vec::new();
            let rc_anchors = reverse_complement_walk(&anchors, k);
            let orientations: [(&[(i32, u64)], &[(usize, u64)], u8); 2] = [
                (&anchors, &forward_positions, 0),
                (&rc_anchors, &reverse_positions, 1),
            ];
            for (walk, positions, orientation) in orientations {
                for &(path, start) in positions.iter() {
                    let mut windows: Vec<(u64, u64)> = walk
                        .iter()
                        .map(|&(_, rel)| {
                            (
                                start.saturating_add(rel),
                                start.saturating_add(rel).saturating_add(k),
                            )
                        })
                        .collect();
                    windows.sort_unstable();
                    let mut intervals: Vec<(u64, u64)> = Vec::with_capacity(windows.len());
                    for (lo, hi) in windows {
                        match intervals.last_mut() {
                            Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
                            _ => intervals.push((lo, hi)),
                        }
                    }
                    let mut partitions: Vec<u32> = Vec::new();
                    for &(_, rel) in walk {
                        for partition in territory.partitions_at(path, start + rel, k) {
                            if !partitions.contains(&partition) {
                                partitions.push(partition);
                            }
                        }
                    }
                    partitions.sort_unstable();
                    occurrences.push(CensusOcc {
                        path,
                        start,
                        orientation,
                        intervals,
                        partitions,
                    });
                }
            }
            let mut occurrence_nodes: Vec<Vec<Vec<u32>>> =
                vec![Vec::with_capacity(1); occurrences.len()];
            let mut by_path: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for (position, occ) in occurrences.iter().enumerate() {
                by_path.entry(occ.path).or_default().push(position);
            }
            for (path, members) in &by_path {
                let lo_min = members
                    .iter()
                    .flat_map(|&position| occurrences[position].intervals.iter().map(|&(lo, _)| lo))
                    .min()
                    .unwrap_or(0);
                let hi_max = members
                    .iter()
                    .flat_map(|&position| {
                        occurrences[position].intervals.iter().map(|&(_, hi)| hi)
                    })
                    .max()
                    .unwrap_or(0);
                if hi_max <= lo_min {
                    continue;
                }
                let mut steps: Vec<(u64, i32)> = panel
                    .walk_path_range(*path, lo_min, hi_max)?
                    .into_iter()
                    .map(|(node, bp)| (bp, node))
                    .collect();
                steps.sort_unstable_by_key(|&(bp, _)| bp);
                for &position in members {
                    for &(lo, hi) in &occurrences[position].intervals {
                        let start = steps.partition_point(|&(bp, _)| bp < lo);
                        let end = steps.partition_point(|&(bp, _)| bp.saturating_add(k) <= hi);
                        let contained = &steps[start..end.max(start)];
                        occurrence_nodes[position]
                            .push(contained.iter().map(|&(_, node)| node.unsigned_abs()).collect());
                    }
                }
            }
            let mut pair_bins = [0u64; MULTI_CENSUS_PAIR_BINS.len()];
            let contigs: Vec<bool> = occurrences
                .iter()
                .map(|occ| {
                    panel.name_map.path_to_name[occ.path]
                        .rsplit('#')
                        .next()
                        .map(|contig| contig == component_contig)
                        .unwrap_or(false)
                })
                .collect();
            let occurrence_node_sets: Vec<BTreeSet<u32>> = occurrence_nodes
                .iter()
                .map(|intervals| intervals.iter().flatten().copied().collect())
                .collect();
            for left in 0..occurrences.len() {
                for right in left + 1..occurrences.len() {
                    let bin = multi_census_pair_bin(
                        &occurrence_node_sets[left].iter().copied().collect::<Vec<_>>(),
                        &occurrence_node_sets[right].iter().copied().collect::<Vec<_>>(),
                        &occurrences[left].partitions,
                        &occurrences[right].partitions,
                        contigs[left] && contigs[right],
                    );
                    pair_bins[bin] += 1;
                }
            }
            let occurrences_json: Vec<serde_json::Value> = occurrences
                .iter()
                .enumerate()
                .map(|(position, occ)| {
                    serde_json::json!({
                        "path": occ.path,
                        "start": occ.start,
                        "orientation": occ.orientation,
                        "partitions": occ.partitions,
                        "intervals": occurrence_nodes[position]
                            .iter()
                            .map(|nodes| serde_json::json!(nodes))
                            .collect::<Vec<_>>(),
                    })
                })
                .collect();
            let pair_bins_json: serde_json::Value = MULTI_CENSUS_PAIR_BINS
                .iter()
                .enumerate()
                .map(|(bin, name)| (name.to_string(), serde_json::json!(pair_bins[bin])))
                .collect::<serde_json::Map<String, serde_json::Value>>()
                .into();
            let line = serde_json::json!({
                "multiplicity": multiplicity,
                "t_r": occurrence_map.len(),
                "share": *multiplicity as f64 / occurrence_map.len() as f64,
                "total_occurrences": total_occurrences,
                "anchors": anchor_count,
                "occurrences": occurrences_json,
                "pair_bins": pair_bins_json,
            });
            Ok(Some((line, *multiplicity)))
        })
        .collect();
    for result in routed {
        if let Some((mut line, multiplicity)) = result? {
            line["record"] = serde_json::json!(routed_patterns);
            writeln!(out, "{line}")?;
            routed_patterns += 1;
            census_mult_total += multiplicity;
        } else {
            dropped_no_occurrence += 1;
        }
    }
    out.flush()?;
    let census_seconds = census_started.elapsed().as_secs_f64();
    let wall_seconds = started.elapsed().as_secs_f64();
    let stats = CutStatsJson {
        reads_scanned: merged.reads_scanned,
        reads_kept: reads.len() as u64,
        records_scanned: merged.records_scanned,
        records_mappable: merged.records_mappable,
        records_unmappable_global_only: merged.records_unmappable,
        reads_with_unmappable_records: merged.reads_with_unmappable,
        touching_patterns: merged.pattern_mult.len() as u64,
        routed_patterns,
        census_mult: census_mult_total,
        territory_rows: territory.territories.len(),
        territory_entries: territory.entries.len(),
        node_map_pairs,
        offset_variants,
        load_seconds,
        scan_seconds,
        census_seconds,
        wall_seconds,
        peak_rss_kb: rss_now_kb(),
    };
    let text = serde_json::to_string_pretty(&stats)?;
    std::fs::write(out_dir.join("cut-stats.json"), text + "\n")?;
    eprintln!(
        "[cut] DONE: {routed_patterns} routed patterns ({} dropped, no occurrence), \
         mult {census_mult_total} [{wall_seconds:.1}s]",
        dropped_no_occurrence
    );
    Ok(())
}

fn main() -> io::Result<()> {
    let options = Options::parse();
    match options.mode {
        Mode::Project {
            panel,
            inputs,
            out_dir,
            map_cap,
            progress,
            limit,
        } => run_project(&panel, &inputs, &out_dir, map_cap, progress, limit),
        Mode::Dedup {
            out_dir,
            occurrences,
            partitions,
        } => run_dedup(&out_dir, &occurrences, partitions),
        Mode::Cut {
            global_panel,
            panel,
            axis,
            bed_directory,
            component,
            sources,
            agc,
            shards_dir,
            out_dir,
        } => run_cut(
            &global_panel,
            &panel,
            &axis,
            &bed_directory,
            &component,
            &sources,
            &agc,
            &shards_dir,
            &out_dir,
        ),
    }
}
