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
use impg::genome_inference::{mem_records, realign_score, territory_frame};
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
        /// THE LOCALITY-BATCHING LEVER (H1): reorder each 16,384-read
        /// arrival batch by input stream - the committed feed's own
        /// region grouping (one chromosome per stream, the cheapest
        /// possible pre-locator: zero pre-pass cost) - before the
        /// parallel dispatch, so every rayon chunk is a run of one
        /// stream's coordinate-adjacent reads (same-locus reads share
        /// their kmer-hash probe lines and GBWT node records: hot cache)
        /// while the packet sends stay spread across all shard channels
        /// exactly as in the committed scheduling. The default path
        /// (off) is the committed scheduling; the per-shard packet
        /// multiset is identical either way.
        #[arg(long, default_value_t = false)]
        region_batching: bool,
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
    /// One locality's cut of the projected evidence. The committed
    /// frame (default) is LOCAL: --panel names the locality's own
    /// syng and --sources its locus.fa, and the cut converts the
    /// projection's global records into the locality's node space
    /// (the node map). --frame global (the campaign's lever 1, the
    /// global-panel call mode) binds the cut to the ONE global panel
    /// alone: no --panel, no --sources — the universe rows translate
    /// to the global frame by their own record headers, the per-row
    /// slices come from the AGC by coordinate, and the projection's
    /// records are kept VERBATIM (no remap, no re-derivation); the
    /// census/cache are emitted in the global frame.
    Cut {
        /// The GLOBAL panel syng prefix (the projection's panel).
        #[arg(long)]
        global_panel: String,
        /// The cut's frame: local (the committed substrate) or global
        /// (the global-panel call mode; --panel/--sources omitted).
        #[arg(long, value_enum, default_value_t = CutFrame::Local)]
        frame: CutFrame,
        /// THE RECORD-ORIGIN CONVENTION (the global frame's
        /// content-origin parameter): the committed substrate's
        /// locus.fa content need not sit at its header's own interval
        /// — the L1-L4 substrates' content sits ONE BASE LATER (d=1,
        /// the uniform delta=1 AGC-identity receipts), the pilot's
        /// at the header interval exactly (d=0). The global frame,
        /// which deletes locus.fa by design, must be TOLD the
        /// convention: the universe rows translate and the
        /// expressibility content is fetched at `start + d`, so the
        /// cut's territory and census sit at the substrate's TRUE
        /// content (the emission-side back-translation cancels it
        /// exactly). Default 0 = the header's own interval. The
        /// LOCAL frame refuses it (it derives the convention per
        /// record against the AGC, fail-closed on neither
        /// convention matching — a stated origin is redundant and
        /// a conflicting one is a caller error).
        #[arg(long, default_value_t = 0)]
        record_origin: u64,
        /// THE CENSUS PROFILE (the global frame's product-default
        /// census, the owner's Ruling 2): restricted (the DEFAULT) =
        /// the committed routed read set — the scan keeps only the
        /// reads of the territory's committed kept-read set, stated
        /// by IMPG_CUT_READSET_FASTA=<fasta> (fail-closed without it:
        /// the product default demands the evidence-set statement,
        /// never a silent full-mass fallback); global-native = the
        /// full routed mass (the record form, byte-identical to the
        /// committed global-frame cut). The A/B receipt
        /// (docs/syng-gmem-bwt/census-profile-ab.md): an accuracy wash
        /// at 4 of 5 localities against a 1.77x-7.68x loci-phase cost,
        /// with ONE truth locus (L2 17) only the global-native census
        /// calls — traced to the conserved-elsewhere dictionary class
        /// of the restriction basis (§10.6.17), not a routing gap.
        /// The LOCAL frame refuses it (the committed local path is
        /// byte-for-byte; the read-set instrument still restricts
        /// there exactly as before).
        #[arg(long, value_enum, default_value_t = CensusProfile::Restricted)]
        census_profile: CensusProfile,
        /// The LOCALITY panel syng prefix (the committed substrate;
        /// required in the local frame).
        #[arg(long, required_unless_present = "frame")]
        panel: Option<String>,
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
        /// source sequences; required in the local frame).
        #[arg(long, required_unless_present = "frame")]
        sources: Option<PathBuf>,
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
    /// THE WHOLE-GENOME FLEET CUT (the amortized-load shape): ONE
    /// global syng + AGC load serving a BATCH of territories per
    /// process - every routed read's records derived by the committed
    /// local MEM accessor (IMPG_CUT_LOCAL_DERIVE=1 required), each
    /// territory's outputs byte-identical to the single cut's.
    CutBatch {
        /// The GLOBAL panel syng prefix (the projection's panel).
        #[arg(long)]
        global_panel: String,
        /// The global panel AGC (the committed --sequence-files input).
        #[arg(long)]
        agc: PathBuf,
        /// The projection output directory (the shard files).
        #[arg(long)]
        shards_dir: PathBuf,
        /// JSONL manifest: one territory per row {name, panel, axis,
        /// bed_directory, component, sources, out_dir} (the committed
        /// single-cut substrate inputs).
        #[arg(long)]
        manifest: PathBuf,
    },
    /// THE SHARED-INDEX CUT LANES (lever 2 — the global frame's
    /// batch-amortized-load form): ONE global syng + AGC load serving a
    /// BATCH of global-frame territories per process. The single global
    /// cut re-pays the ~325s global GBWT load pole every run; the batch
    /// amortizes it over the whole manifest. Each territory's outputs
    /// are byte-identical to the single global cut's (the gate).
    CutBatchGlobal {
        /// The GLOBAL panel syng prefix (the projection's panel).
        #[arg(long)]
        global_panel: String,
        /// The global panel AGC (the committed --sequence-files input).
        #[arg(long)]
        agc: PathBuf,
        /// The projection output directory (the shard files; ONE scan
        /// serves every manifest territory).
        #[arg(long)]
        shards_dir: PathBuf,
        /// JSONL manifest: one territory per row {name, axis,
        /// bed_directory, component, out_dir, record_origin} (the
        /// single global cut's own substrate inputs; record_origin
        /// defaults to 0 = the header's own interval).
        #[arg(long)]
        manifest: PathBuf,
    },
    /// THE SHARED-INDEX CUT+CALL SEAM (the combined process): ONE
    /// global panel + AGC load serving a BATCH of global-frame
    /// territories' cuts AND calls — the cut phase is cut-batch-global's
    /// own core over the loaded panel, then each manifest row carrying
    /// a "call" object is CALLED IN-PROCESS over the SAME loaded panel
    /// (the injected-index form of call-diplotypes: the call's
    /// per-process inputs pole — the panel's own load — vanishes into
    /// the shared load). The call rows' outputs are the committed CLI
    /// shape (calls.jsonl + the instrument receipts, the output-model
    /// discipline); the gate is byte-identity vs the committed-path
    /// receipts.
    CutCallBatchGlobal {
        /// The GLOBAL panel syng prefix (the projection's panel).
        #[arg(long)]
        global_panel: String,
        /// The global panel AGC (the committed --sequence-files input).
        #[arg(long)]
        agc: PathBuf,
        /// The projection output directory (the shard files; ONE scan
        /// serves every manifest territory).
        #[arg(long)]
        shards_dir: PathBuf,
        /// JSONL manifest: one territory per row {name, axis,
        /// bed_directory, component, out_dir, record_origin, call?}
        /// where call = {reads, partition_graphs, loci, out_dir}.
        #[arg(long)]
        manifest: PathBuf,
        /// The whole-panel routes inventory directory (graph.json +
        /// sources — the call's --routes input, shared by every call
        /// row).
        #[arg(long)]
        routes: PathBuf,
        /// Per-call resident-set guard in GiB (0 = unguarded — the
        /// combined process holds the shared panel by design, the
        /// cut batch's own rule; a per-call guard sized like the
        /// single-process default would kill the shared-load shape).
        #[arg(long, default_value_t = 0.0)]
        call_rss_budget_gib: f64,
    },
    /// THE FORK-AFTER-LOAD COW PANEL SHARE (SS10.6.9's form (a), the
    /// slot-ceiling lever): the parent loads the panel + AGC ONCE
    /// and stays the load owner; ONE CHILD PER MANIFEST ROW forks
    /// and runs its territory's cut + call over the inherited
    /// COW-shared panel, then exits. Receipts the per-child smaps
    /// COW arithmetic at open and close (the COW-hold gate).
    ForkCutCallBatchGlobal {
        /// The GLOBAL panel syng prefix (the projection's panel).
        #[arg(long)]
        global_panel: String,
        /// The global panel AGC (the committed --sequence-files input).
        #[arg(long)]
        agc: PathBuf,
        /// The projection output directory (the shard files; each
        /// child scans them for its own territory).
        #[arg(long)]
        shards_dir: PathBuf,
        /// JSONL manifest: one territory per row {name, axis,
        /// bed_directory, component, out_dir, record_origin, call?}
        /// where call = {reads, partition_graphs, loci, out_dir}.
        #[arg(long)]
        manifest: PathBuf,
        /// The whole-panel routes inventory directory (graph.json +
        /// sources — the call's --routes input, shared by every call
        /// row).
        #[arg(long)]
        routes: PathBuf,
        /// Per-call resident-set guard in GiB (0 = unguarded — the
        /// child holds the shared panel by design, the combined
        /// process's own rule).
        #[arg(long, default_value_t = 0.0)]
        call_rss_budget_gib: f64,
        /// Manifest rows per forked child (1 = one child per
        /// territory, the maximal-share form; the manifest's own
        /// batch width amortizes the cut scan's working set within a
        /// child — the fleet-shape parameter, never an algorithm
        /// constant).
        #[arg(long, default_value_t = 1)]
        rows_per_child: usize,
    },
    /// THE YIELD-ANOMALY AUTOPSY TRACE (receipt-side diagnostic; no
    /// product machinery): one scan of the projection shards that
    /// classifies, per locality and per read-class (sequence sets
    /// supplied as FASTA), every record's remap fate - touching,
    /// mappable-not-touching, unmappable - with the unmapped anchors'
    /// in-slice vs out-of-slice split (the dictionary-superset class
    /// vs the walk-exits-the-slices class) and the partial-remap
    /// recovery projection (would the mapped-subwalk touch?).
    Trace {
        /// The GLOBAL panel syng prefix (the projection's panel).
        #[arg(long)]
        global_panel: String,
        /// The global panel AGC.
        #[arg(long)]
        agc: PathBuf,
        /// The projection output directory (the shard files).
        #[arg(long)]
        shards_dir: PathBuf,
        /// JSON list of locality specs: [{name, panel, axis,
        /// bed_directory, component, sources}, ...].
        #[arg(long)]
        localities: PathBuf,
        /// JSON list of class specs: [{class: "<name>:<label>",
        /// fasta: path}, ...]. Every read whose sequence is in a class
        /// set is classified against that locality.
        #[arg(long)]
        classes: PathBuf,
        /// Optional FASTA of query reads: per-record detail dump.
        #[arg(long)]
        queries: Option<PathBuf>,
        /// Output directory (trace-agg.json, trace-detail.jsonl).
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
    region_batching: bool,
) -> io::Result<()> {
    if batch.is_empty() {
        return Ok(());
    }
    // THE LOCALITY-BATCHING SCHEDULING (the sorted form): reorder the
    // arrival batch by input stream (the committed feed's own region
    // grouping - one chromosome per stream) before the parallel
    // dispatch, so every rayon chunk is a run of ONE stream's
    // coordinate-adjacent reads (same-locus reads share their kmer-hash
    // probe lines and GBWT node records - hot cache) while the sends
    // stay spread across ALL shard channels exactly as in the committed
    // scheduling (no single-writer burst). Within-batch packet order is
    // already the committed send-in-par scheduling's own
    // nondeterminism; the per-shard packet multiset is unchanged.
    if region_batching {
        batch.sort_unstable_by_key(|(index, _, _)| *index);
    }
    // The whole per-read pipeline (MEM derivation AND packet encoding)
    // is pure per-read work; only the channel send and the stats are
    // sequential.
    struct Encoded {
        index: usize,
        seq_len: usize,
        reads_with_mems: bool,
        records: usize,
    }
    let results: Vec<io::Result<Encoded>> = batch
        .par_iter()
        .map(|(index, seq, _qual)| {
            let records = derive_read(panel, seq)?;
            let tokens: Vec<Vec<u64>> = records.iter().map(|(t, _)| t.clone()).collect();
            let bytes = encode_read_packet(seq.as_slice(), &records);
            // The channel send is thread-safe and moves the packet
            // entirely; the main thread never touches the payload.
            senders[*index]
                .send(WriteMsg::Packet { bytes, tokens })
                .map_err(|_| invalid("shard writer exited"))?;
            Ok(Encoded {
                index: *index,
                seq_len: seq.len(),
                reads_with_mems: !records.is_empty(),
                records: records.len(),
            })
        })
        .collect();
    batch.clear();
    for result in results {
        let encoded = result?;
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
    region_batching: bool,
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
        let (tx, rx) = crossbeam_channel::bounded::<WriteMsg>(8192);
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
            process_batch(&panel, &names, &senders, &mut batch, &mut totals, region_batching)?;
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
    process_batch(&panel, &names, &senders, &mut batch, &mut totals, region_batching)?;
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

/// The committed local frame's lane fetch over the locus.fa record
/// slices (clip at the record end, empty when the interval is empty) —
/// the exact semantics the inlined lookup had (byte-identical).
fn lane_fetch_over(sources: &HashMap<String, Vec<u8>>) -> impl Fn(&str, u64, u64) -> io::Result<Vec<u8>> + '_ {
    move |name: &str, lo: u64, hi: u64| {
        let lane = sources
            .get(name)
            .ok_or_else(|| invalid(&format!("source {name} absent")))?;
        let hi_clip = hi.min(lane.len() as u64);
        if lo >= hi_clip {
            Ok(Vec::new())
        } else {
            Ok(lane[lo as usize..hi_clip as usize].to_vec())
        }
    }
}

/// The committed territory step construction (the canonical scheme),
/// verbatim semantics: forward steps from the panel's own stored path
/// walk (sampled sidecars exactly as the committed census saw them),
/// the rc frame's steps from the complete raw extraction, and the
/// per-position canonical-frame selection.
/// THE LANE FETCH (the frame seam): `fetch_lane(name, lo, hi)` returns
/// the panel-frame bytes [lo, hi) clipped at the lane's end, empty
/// when the interval is empty — the local frame serves the locus.fa
/// record slices (the committed bytes), the global frame the AGC's
/// coordinate slices (the same bytes the locus.fa extraction produced,
/// by the cut's own AGC-identity receipt).
fn build_territory_index(
    panel: &SyngIndex,
    rows: &[UniverseRow],
    fetch_lane: &(dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>> + Sync),
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
            let seq: Vec<u8> = fetch_lane(&row.path_name, seq_lo, seq_hi)?;
            let rc_seq = impg::graph::reverse_complement(&seq);
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
    /// The cut's frame: "local" (the committed substrate path,
    /// node-map conversion) or "global" (the global-panel call mode,
    /// lever 1: the projection's records kept verbatim).
    frame: String,
    /// THE CENSUS PROFILE (the global frame's product-default census,
    /// the owner's Ruling 2): "restricted" (the committed routed read
    /// set) or "global-native" (the full routed mass, the record
    /// form). Emitted in the GLOBAL frame only — the local frame has
    /// no profile concept and keeps its committed stats bytes.
    #[serde(skip_serializing_if = "Option::is_none")]
    census_profile: Option<String>,
    /// The stated content-origin convention (the --record-origin d:
    /// the substrate's locus.fa content sits at `header_start + d`;
    /// the L1-L4 substrates: 1, the pilot: 0).
    record_origin: u64,
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
    /// The env-gated attribution mode receipts (false = the committed
    /// all-or-nothing remap, byte-identical behavior).
    partial_remap: bool,
    local_derive: bool,
    /// The global-frame receipt (lever 1: the projection's records
    /// kept verbatim; the census/cache are global-frame).
    global_frame: bool,
    reads_routed: u64,
    reads_unrouted: u64,
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
    /// local-derive mode: reads carrying the global routing signal
    /// (at least one globally-mappable record node).
    reads_routed: u64,
    /// local-derive mode: reads without the routing signal.
    reads_unrouted: u64,
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

/// The cut's attribution mode, all env-gated (default Strict = the
/// committed all-or-nothing remap, byte-identical).
/// - PartialSubwalk (IMPG_CUT_PARTIAL_REMAP=1): the locally-expressible
///   subwalk attribution (fragmentary patterns; binds on some
///   substrates, rejected by the unchanged CLI's binding on others).
/// - LocalDerive (IMPG_CUT_LOCAL_DERIVE=1): ROUTE GLOBALLY, ATTRIBUTE
///   LOCALLY - the global projection is the read ROUTER (a read is a
///   candidate iff any of its global records carries a
///   locality-mappable node), and every candidate read's records are
///   then derived by the committed local MEM accessor itself
///   (mem_records::tagged_mem_records against the locality panel,
///   verbatim the committed census re-derivation and derive-cache
///   builder), so the emitted census records are committed-structured
///   local MEM walks by construction.
/// THE CUT'S FRAME (the global-panel call mode's switch, lever 1:
/// docs/syng-gmem-bwt/global-panel-call-mode.md). Local (the default)
/// is the committed substrate path byte-for-byte: the locality's own
/// syng + locus.fa, and the projection's global records converted
/// into the locality's node space (the node map). Global binds the
/// cut to the ONE global panel alone — conversion-not-materialization:
/// the universe rows translate by their own record headers, the
/// per-row slices come from the AGC by coordinate, and the
/// projection's records are kept verbatim (no remap, no
/// re-derivation); the census/cache come out in the global frame.
#[derive(clap::ValueEnum, Clone, Copy, PartialEq, Eq, Debug)]
enum CutFrame {
    Local,
    Global,
}

/// THE CENSUS PROFILE (the owner's Ruling 2, executed: the
/// product-default census of the global-frame call mode). The
/// restricted profile is the committed routed read set — the scan
/// keeps only the reads of the territory's committed kept-read set
/// (the evidence mass the original local-frame receipts banked),
/// stated by the IMPG_CUT_READSET_FASTA instrument (the restriction
/// basis, untouched). The global-native profile is the full routed
/// mass (the receipts-of-record's superset census, kept for the
/// record). The A/B receipt that decided this
/// (docs/syng-gmem-bwt/census-profile-ab.md): a wash on accuracy at
/// 4 of 5 localities against a 1.77x-7.68x loci-phase cost, with
/// ONE truth locus (L2 17) that only the global-native census calls —
/// the autopsy (§10.6.17) traced that locus to the conserved-elsewhere
/// dictionary class, an evidence-set property of the restriction
/// basis, not a routing gap; no adjacency rule recovers it.
#[derive(clap::ValueEnum, Clone, Copy, PartialEq, Eq, Debug)]
enum CensusProfile {
    /// The committed routed read set (the product default; requires
    /// the IMPG_CUT_READSET_FASTA basis, fail-closed without it).
    Restricted,
    /// The full routed mass (the record form; byte-identical to the
    /// committed global-frame cut when no basis is stated).
    GlobalNative,
}

impl CensusProfile {
    fn label(self) -> &'static str {
        match self {
            CensusProfile::Restricted => "restricted",
            CensusProfile::GlobalNative => "global-native",
        }
    }

    /// THE PROFILE RESOLUTION (the three-step interface's own rule,
    /// fail-closed, no silent defaults): the frame owns the reach —
    /// the LOCAL frame keeps the committed behavior (the profile
    /// concept is the global frame's; the read-set instrument still
    /// restricts there exactly as before); the GLOBAL frame resolves
    /// the effective profile: a stated basis (IMPG_CUT_READSET_FASTA)
    /// overrides to the restricted scan in EITHER profile (the
    /// instrument is the override it always was), the restricted
    /// profile without a basis is a NAMED caller error (the product
    /// default demands the evidence-set statement — the committed
    /// kept-read set — never a silent full-mass fallback), and the
    /// global-native profile without a basis is the full-mass scan,
    /// byte-identical to the committed global-frame cut.
    fn resolve(
        frame: CutFrame,
        profile: CensusProfile,
        readset_basis: Option<&str>,
    ) -> io::Result<(CensusProfile, bool)> {
        match frame {
            CutFrame::Local => Ok((profile, false)),
            CutFrame::Global => match (profile, readset_basis) {
                (CensusProfile::Restricted, None) => Err(invalid(
                    "the restricted census profile (the default) requires the committed \
                     kept-read set: state it with IMPG_CUT_READSET_FASTA=<fasta>, or select \
                     --census-profile global-native (the record form, the full routed mass)",
                )),
                (CensusProfile::GlobalNative, Some(_)) => {
                    eprintln!(
                        "[cut] census profile: the IMPG_CUT_READSET_FASTA basis overrides \
                         --census-profile global-native — the scan runs RESTRICTED"
                    );
                    Ok((CensusProfile::Restricted, true))
                }
                (CensusProfile::GlobalNative, None) => Ok((CensusProfile::GlobalNative, false)),
                (CensusProfile::Restricted, Some(_)) => Ok((CensusProfile::Restricted, false)),
            },
        }
    }
}

impl CutFrame {
    /// The local-frame-only argument gate: a `--frame global` run
    /// refuses a locality substrate argument loudly (never ignores
    /// it); the local frame requires it.
    fn select<T>(self, value: Option<T>, flag: &str) -> io::Result<Option<T>> {
        match (self, value) {
            (CutFrame::Local, None) => Err(invalid(&format!(
                "the local frame requires {flag}"
            ))),
            (CutFrame::Global, Some(_)) => Err(invalid(&format!(
                "the global frame takes no {flag} (the ONE global panel \
                 is the frame; conversion-not-materialization)"
            ))),
            (_, value) => Ok(value),
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
enum RemapMode {
    Strict,
    PartialSubwalk,
    LocalDerive,
    /// The GLOBAL frame (CutFrame::Global): the projection's records
    /// are already panel-frame — kept verbatim; the territory touch is
    /// computed on their own nodes; nothing is remapped, nothing is
    /// re-derived. The census/cache are emitted in the global frame.
    GlobalFrame,
}

fn scan_shard_file(
    path: &Path,
    node_map: &HashMap<i32, i32>,
    territory_nodes_abs: &HashSet<u32>,
    mode: RemapMode,
    local_panel: Option<&SyngIndex>,
) -> io::Result<FileCut> {
    // THE CENSUS-SUPERSET SHADOW INSTRUMENT (receipt-side diagnostic,
    // default OFF — the committed behavior is byte-identical): with
    // IMPG_CUT_READSET_FASTA naming a FASTA of read sequences, the scan
    // keeps ONLY the reads whose sequence is in that set (by FNV hash
    // of the exact bytes), everything else is scanned-and-skipped. The
    // owner's census-superset question (§10.3: the global cut keeps
    // 2.8x the committed read mass) is an ACCURACY decision; this
    // instrument makes its PERFORMANCE shadow measurable without
    // touching the answer profile — the emitted census/cache are the
    // same machinery at the stated read mass, nothing more.
    let readset: Option<std::collections::HashSet<u64>> = std::env::var(
        "IMPG_CUT_READSET_FASTA",
    )
    .ok()
    .filter(|value| !value.is_empty())
    .map(|path| -> io::Result<_> {
        let text = std::fs::read_to_string(&path)?;
        let mut set = std::collections::HashSet::new();
        let mut seq: Vec<u8> = Vec::new();
        for line in text.lines() {
            if let Some(header) = line.strip_prefix('>') {
                let _ = header;
                if !seq.is_empty() {
                    set.insert(fnv1a64(&seq));
                }
                seq.clear();
            } else {
                seq.extend_from_slice(line.trim().as_bytes());
            }
        }
        if !seq.is_empty() {
            set.insert(fnv1a64(&seq));
        }
        eprintln!(
            "[cut] read-set shadow instrument: {} reads from {path}",
            set.len()
        );
        Ok(set)
    })
    .transpose()?;
    let readset_ref = readset.as_ref();
    let mut cut = FileCut {
        reads_scanned: 0,
        records_scanned: 0,
        records_mappable: 0,
        records_unmappable: 0,
        reads_with_unmappable: 0,
        reads_routed: 0,
        reads_unrouted: 0,
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
        let decoded: Vec<DecodedRead> = match readset_ref {
            Some(set) => decoded
                .into_iter()
                .filter(|read| set.contains(&fnv1a64(&read.seq)))
                .collect(),
            None => decoded,
        };
        let remapped: Vec<io::Result<(
            Vec<Option<(Vec<u64>, OwnWalk, bool)>>,
            bool,
        )>> = decoded
            .par_iter()
            .map(|read| {
                // (THE GLOBAL FRAME runs this same committed remap flow
                // with the IDENTITY node map over the territory's own
                // canonical anchors — the expressibility rule is the
                // node map's real role, preserved in the panel's own
                // frame; the tokens and walks come out verbatim.)
                if mode == RemapMode::LocalDerive {
                    // the global routing signal: any globally-mappable node
                    let routed = read.records.iter().any(|(_, walk)| {
                        walk.iter().any(|&(node, _)| node_map.contains_key(&node))
                    });
                    if !routed {
                        return Ok((Vec::new(), false));
                    }
                    // ROUTE GLOBALLY, ATTRIBUTE LOCALLY: the candidate's
                    // records are the committed local MEM derivation
                    // itself (verbatim committed structure).
                    let panel = local_panel
                        .ok_or_else(|| invalid("local-derive mode needs the locality panel"))?;
                    let local = derive_read(panel, &read.seq)?;
                    let out: Vec<Option<(Vec<u64>, OwnWalk, bool)>> = local
                        .into_iter()
                        .map(|(tokens, walk)| {
                            let touch = walk
                                .iter()
                                .any(|&(node, _)| {
                                    territory_nodes_abs.contains(&node.unsigned_abs())
                                });
                            Some((tokens, walk, touch))
                        })
                        .collect();
                    return Ok((out, true));
                }
                let mut out = Vec::with_capacity(read.records.len());
                cut_records(
                    &read.records,
                    node_map,
                    territory_nodes_abs,
                    mode == RemapMode::PartialSubwalk,
                    &mut out,
                )?;
                Ok((out, false))
            })
            .collect();
        for (read, result) in decoded.into_iter().zip(remapped) {
            let (records, routed_flag) = result?;
            if mode == RemapMode::LocalDerive {
                cut.records_scanned += read.records.len() as u64;
                if routed_flag {
                    cut.reads_routed += 1;
                    cut.records_mappable += records.len() as u64;
                } else {
                    cut.reads_unrouted += 1;
                    cut.reads_with_unmappable += 1;
                }
            } else {
                cut.records_scanned += records.len() as u64;
            }
            let mut kept: Vec<(u32, OwnWalk)> = Vec::new();
            let mut any_touch = false;
            for record in &records {
                match record {
                    None => {
                        cut.records_unmappable += 1;
                    }
                    Some((tokens, walk, touch)) => {
                        if mode != RemapMode::LocalDerive {
                            cut.records_mappable += 1;
                        }
                        if *touch {
                            any_touch = true;
                            *cut.pattern_mult.entry(tokens.clone()).or_insert(0) += 1;
                        }
                        let key = cut.intern_key(tokens);
                        kept.push((key, walk.clone()));
                    }
                }
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
/// With `partial` (the env-gated attribution correction
/// IMPG_CUT_PARTIAL_REMAP=1), a record whose walk carries global-only
/// anchors is attributed by its LOCALLY-EXPRESSIBLE SUBWALK instead of
/// being dropped: the OWN-FRAME walk is filtered to the steps whose
/// nodes the locality panel expresses (true own positions preserved),
/// and the canonical tokens are then DERIVED from the surviving own
/// steps by the CLI's own mirror formula - so the emitted
/// (tokens, own walk) pair satisfies the unchanged CLI's cache binding
/// by construction (verified inline for both mirror frames; a record
/// whose pair validates in neither frame is dropped). Records with no
/// expressible step stay dropped. No thresholds; the committed default
/// (unset) is the all-or-nothing remap, byte-identical to the committed
/// cuts (proven by the cut2-control-L4 byte-identity run).
fn cut_records(
    records: &[(Vec<u64>, OwnWalk)],
    node_map: &HashMap<i32, i32>,
    territory_nodes_abs: &HashSet<u32>,
    partial: bool,
    out: &mut Vec<Option<(Vec<u64>, OwnWalk, bool)>>,
) -> io::Result<()> {
    for (tokens, walk) in records {
        if partial {
            let pair = partial_remap(tokens, walk, node_map)?;
            match pair {
                Some((new_tokens, own_kept)) => {
                    let anchors = decode_tokens(&new_tokens)?;
                    let touch = anchors
                        .iter()
                        .any(|&(node, _)| territory_nodes_abs.contains(&node.unsigned_abs()));
                    out.push(Some((new_tokens, own_kept, touch)));
                }
                None => out.push(None),
            }
            continue;
        }
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

/// The partial-remap pair construction: filter the own-frame walk to
/// the locality-expressible steps, derive the canonical tokens from the
/// survivors by the mirror formula (both frames tried; the CLI's own
/// validation decides), and return the pair only when it binds.
fn partial_remap(
    tokens: &[u64],
    walk: &OwnWalk,
    node_map: &HashMap<i32, i32>,
) -> io::Result<Option<(Vec<u64>, OwnWalk)>> {
    let mut own_kept: OwnWalk = Vec::with_capacity(walk.len());
    for &(node, pos) in walk {
        if let Some(&local) = node_map.get(&node) {
            own_kept.push((local, pos));
        }
    }
    if own_kept.is_empty() {
        return Ok(None);
    }
    let e_own = match encode_walk(&own_kept.iter().map(|&(n, p)| (n, p as u64)).collect::<Vec<_>>())
    {
        Ok(encoded) => encoded,
        Err(_) => return Ok(None),
    };
    let mirrored = match encode_walk(
        &walk.iter().map(|&(n, p)| (n, p as u64)).collect::<Vec<_>>(),
    ) {
        Ok(original) => tokens != original.as_slice(),
        Err(_) => return Ok(None),
    };
    let m = own_kept.len();
    let w_lo = own_kept[0].1 as u64;
    let d: Vec<u64> = own_kept.iter().map(|&(_, p)| p as u64 - w_lo).collect();
    let d_last = *d.last().unwrap();
    for frame in [mirrored, !mirrored] {
        let mut canonical_list: Vec<(i32, u64)> = Vec::with_capacity(m);
        for j in 0..m {
            let index = if frame { m - 1 - j } else { j };
            let (node, _) = own_kept[index];
            canonical_list.push((
                if frame { -node } else { node },
                if frame { d_last - d[index] } else { d[index] },
            ));
        }
        let new_tokens = match encode_walk(&canonical_list) {
            Ok(encoded) => canonical(&encoded),
            Err(_) => continue,
        };
        if pair_binds(&new_tokens, &own_kept, &e_own)? {
            return Ok(Some((new_tokens, own_kept)));
        }
    }
    Ok(None)
}

/// The unchanged CLI's own cache-binding validation (realign_score's
/// key_shape), verbatim semantics: the own walk must mirror the
/// canonical tokens' walk (nodes and positions) in one of the two
/// frames. `e_own` is the own walk's encoded form (the CLI derives the
/// mirrored flag from exactly this comparison).
fn pair_binds(tokens: &[u64], own: &OwnWalk, e_own: &[u64]) -> io::Result<bool> {
    let mirrored = e_own != tokens;
    let canonical = decode_tokens(tokens)?;
    let m = canonical.len();
    if own.len() != m {
        return Ok(false);
    }
    let w_lo = own.first().map(|&(_, p)| p as u64).unwrap_or(0);
    let last = canonical.last().map(|&(_, p)| p).unwrap_or(0);
    for (t, &(own_node, own_pos)) in own.iter().enumerate() {
        let j = if mirrored { m - 1 - t } else { t };
        let expect_node = if mirrored { -canonical[j].0 } else { canonical[j].0 };
        let expect_pos = if mirrored {
            w_lo + (last - canonical[j].1)
        } else {
            w_lo + canonical[j].1
        };
        if own_node != expect_node || own_pos as u64 != expect_pos {
            return Ok(false);
        }
    }
    Ok(true)
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
    frame: CutFrame,
    record_origin: u64,
    census_profile: CensusProfile,
    panel_prefix: Option<&str>,
    axis_path: &Path,
    bed_directory: &Path,
    component: &str,
    sources_path: Option<&Path>,
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
    eprintln!(
        "[cut] loading panels (frame {:?})",
        frame
    );
    let load_started = Instant::now();
    let agc = AgcIndex::build_from_files(&[agc_path.display().to_string()])?;
    let global = SyngIndex::load(global_panel_prefix, SyncmerParams::default())?;
    // THE FRAME'S PANEL: the local frame loads the locality's own syng
    // beside the global one (the node-map conversion needs both); the
    // global frame binds to the ONE global panel alone — no locality
    // index is built, nothing is materialized per territory.
    let local_panel = match frame {
        CutFrame::Local => {
            let panel = SyngIndex::load(
                panel_prefix.ok_or_else(|| invalid("the local frame requires --panel"))?,
                SyncmerParams::default(),
            )?;
            ensure(
                global.syncmer_length_bp() as u64
                    == panel.syncmer_length_bp() as u64,
                "global and locality syncmer lengths differ",
            )?;
            Some(panel)
        }
        CutFrame::Global => None,
    };
    // The machinery's panel: whichever frame the cut runs in, the
    // census/cache below are emitted against THIS panel's frame.
    let panel: &SyngIndex = local_panel.as_ref().unwrap_or(&global);
    let k = panel.syncmer_length_bp() as u64;

    let axis: AxisFile =
        serde_json::from_str(&std::fs::read_to_string(axis_path)?).map_err(|e| invalid(format!("axis parse: {e}")))?;
    let universe = load_universe_component(&axis, bed_directory, component)?;
    eprintln!("[cut] universe rows: {}", universe.len());

    // THE FRAME SETUP: the local frame materializes the locus.fa
    // slices and derives the global→local node correspondence (the
    // committed path, byte-identical); the global frame TRANSLATES the
    // universe rows by their own record headers and fetches the
    // per-row slices from the AGC by coordinate — no locus.fa, no
    // locality syng, no node map (conversion-not-materialization).
    let mut node_map: HashMap<i32, i32> = HashMap::new();
    let mut node_map_pairs = 0usize;
    let mut global_only_anchor_total = 0usize;
    let mut offset_variants: BTreeMap<String, u64> = BTreeMap::new();
    let territory;
    let territory_nodes_abs: HashSet<u32>;
    let component_contig_source: String;
    match frame {
        CutFrame::Local => {
            // The locus sources (the committed substrate's sequences).
            let mut sources: HashMap<String, Vec<u8>> = HashMap::new();
            {
                let text = std::fs::read_to_string(
                    sources_path.ok_or_else(|| invalid("the local frame requires --sources"))?,
                )?;
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

            let sources_ref = &sources;
            let fetch_lane =
                move |name: &str, lo: u64, hi: u64| -> io::Result<Vec<u8>> {
                    let lane = sources_ref
                        .get(name)
                        .ok_or_else(|| invalid(&format!("source {name} absent")))?;
                    let hi_clip = hi.min(lane.len() as u64);
                    if lo >= hi_clip {
                        Ok(Vec::new())
                    } else {
                        Ok(lane[lo as usize..hi_clip as usize].to_vec())
                    }
                };
            territory = build_territory_index(panel, &universe, &fetch_lane, k)?;
            eprintln!(
                "[cut] territory: {} rows, {} entries",
                territory.territories.len(),
                territory.entries.len()
            );
            territory_nodes_abs = territory
                .entries
                .iter()
                .map(|&(node, _, _)| node.unsigned_abs())
                .collect();
            component_contig_source = component.to_string();

            // The per-path node correspondence: the complete canonical
            // walk of the SAME AGC-verified sequence through both syngs,
            // paired by position. Each locus.fa header states its own
            // provenance — the global panel path and interval
            // "path:start-end" (the refined-bed rows and the axis-window
            // records alike) — so the correspondence is derived from the
            // substrate's own statement, with the offset convention
            // (exact or one-base-earlier) verified by sequence identity
            // against the AGC.
            for (header, record) in &sources {
                let (global_name, s, e) = parse_locus_header(header)?;
                ensure(
                    (e - s) as usize == record.len(),
                    &format!("locus record {header} length differs from its header interval"),
                )?;
                // The offset candidates cover both committed extraction
                // conventions (the generic exact-bed intervals and the
                // pilot's one-base-earlier records) under the AGC range
                // fetch's own 1-based-inclusive coordinate system;
                // sequence identity decides.
                let mut delta: Option<i64> = None;
                for candidate in [0i64, -1, 1] {
                    let lo = (s as i64 + candidate).max(0) as usize;
                    let hi = lo + record.len();
                    if let Ok(fetched) =
                        agc.fetch_sequence(&global_name, lo as i32, hi as i32)
                    {
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
                let local_walk = mem_records::canonical_anchor_walk(panel, record)?;
                let global_walk = mem_records::canonical_anchor_walk(&global, &global_seq)?;
                // The locality dictionary is a subset of the global one, so
                // every locality anchor position must appear in the global
                // walk; the global walk may carry EXTRA anchors (the
                // superset class). Pair on the locality positions (the true
                // correspondence); the global-only positions have no
                // locality node and stay unmapped.
                let global_by_pos: HashMap<u64, i32> = global_walk
                    .iter()
                    .map(|&(node, pos)| (pos, node))
                    .collect();
                for &(lnode, lpos) in &local_walk {
                    let gnode = *global_by_pos.get(&lpos).ok_or_else(|| {
                        invalid("locality anchor position absent from the global walk")
                    })?;
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
        }
        CutFrame::Global => {
            // THE GLOBAL FRAME (lever 1): the universe rows translate to
            // the global panel's frame by their own record headers —
            // the re-anchor principle (k-independent panel
            // coordinates; conversion-not-materialization) — AT THE
            // SUBSTRATE'S STATED CONTENT ORIGIN `start + d` (the
            // --record-origin convention: d=1 = the committed L1-L4
            // substrates' one-base-later content, the uniform delta=1
            // AGC-identity receipts; d=0 = the pilot's exact headers).
            // The AGC is the sequence authority, fetched by coordinate
            // at the same origin — the true content the locus.fa
            // extraction produced.
            let universe_global: Vec<UniverseRow> = universe
                .iter()
                .map(|row| {
                    let translated = territory_frame::translate_row(
                        &row.path_name,
                        row.start,
                        row.end,
                        record_origin,
                    )?;
                    Ok::<_, io::Error>(UniverseRow {
                        partition: row.partition,
                        path_name: translated.path_name,
                        start: translated.start,
                        end: translated.end,
                    })
                })
                .collect::<io::Result<_>>()?;
            let path_lengths: HashMap<&str, u64> = global
                .name_map
                .path_to_name
                .iter()
                .zip(&global.name_map.path_to_length)
                .map(|(name, &length)| (name.as_str(), length))
                .collect();
            let agc_ref = &agc;
            let fetch_lane =
                move |name: &str, lo: u64, hi: u64| -> io::Result<Vec<u8>> {
                    let path_len = *path_lengths
                        .get(name)
                        .ok_or_else(|| invalid("universe path absent from the global panel"))?;
                    let hi = hi.min(path_len);
                    if lo >= hi {
                        return Ok(Vec::new());
                    }
                    // fetch(A, B) IS the 0-based [A, B) slice (the
                    // frame-offset probe's receipt: the canonical walk
                    // of fetch(s, s+len) aligns EXACTLY with the
                    // panel's own walk steps — 82/82 anchors — while
                    // the shifted candidates are off-by-one, 83/83).
                    Ok(agc_ref
                        .fetch_sequence(name, lo as i32, hi as i32)
                        .map_err(|e| invalid(format!("AGC fetch {name}: {e}")))?)
                };
            territory = build_territory_index(&global, &universe_global, &fetch_lane, k)?;
            eprintln!(
                "[cut] territory (global frame): {} rows, {} entries",
                territory.territories.len(),
                territory.entries.len()
            );
            territory_nodes_abs = territory
                .entries
                .iter()
                .map(|&(node, _, _)| node.unsigned_abs())
                .collect();
            // The pair-bins' contig comparison runs in the panel's own
            // frame: the translated component's contig suffix.
            component_contig_source =
                territory_frame::split_local_name(component)?.0.to_string();
            // THE EXPRESSIBILITY MAP (the node map's real role, in the
            // panel's own frame): the committed local cut's node map is
            // not merely an id conversion — it is the EVIDENCE FILTER
            // (a record survives iff EVERY anchor is one of the
            // substrate's own canonical anchors: the slice-context
            // canonical walk; the read-edge qualification class is
            // dropped by construction). The global frame derives the
            // SAME set from the global panel's own grid: the canonical
            // anchor walk of each record's content — fetched from the
            // AGC at the record's CONTENT origin `start + d` (the
            // stated --record-origin convention: the same bytes the
            // locus.fa extraction produced; a fetch that clips short
            // at the path end fails closed with a named error, never
            // a truncated walk) — and remaps every node to ITSELF.
            // The projection's records then ride the committed remap
            // flow unchanged: the expressible ones verbatim, the rest
            // dropped exactly as the strict local cut drops them.
            let mut expressibility_pairs = 0usize;
            let mut record_names: BTreeSet<&str> = BTreeSet::new();
            for row in &universe {
                record_names.insert(row.path_name.as_str());
            }
            for header in &record_names {
                let (global_name, start, end) = parse_locus_header(header)?;
                let length = (end - start) as usize;
                let content_start = start
                    .checked_add(record_origin)
                    .ok_or_else(|| invalid(format!(
                        "record {header}: the content origin start+{record_origin} overflows"
                    )))?;
                let content = agc
                    .fetch_sequence(
                        &global_name,
                        content_start as i32,
                        (content_start + length as u64) as i32,
                    )
                    .map_err(|e| invalid(format!("AGC fetch {global_name}: {e}")))?;
                ensure(
                    content.len() == length,
                    &format!(
                        "record {header}: the content fetch at origin +{record_origin} \
                         returned {} bytes, not the record's {length} (the content origin \
                         leaves the path; state the substrate's own convention)",
                        content.len()
                    ),
                )?;
                let walk = mem_records::canonical_anchor_walk(&global, &content)?;
                for &(node, _position) in &walk {
                    match node_map.entry(node) {
                        std::collections::hash_map::Entry::Occupied(o) => {
                            ensure(*o.get() == node, "conflicting identity entry")?;
                        }
                        std::collections::hash_map::Entry::Vacant(v) => {
                            v.insert(node);
                            expressibility_pairs += 1;
                        }
                    }
                }
            }
            node_map_pairs = expressibility_pairs;
            eprintln!(
                "[cut] global frame: the identity expressibility map over \
                 {} records: {} canonical anchors; component contig {}",
                record_names.len(),
                expressibility_pairs,
                component_contig_source
            );
        }
    }
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
    // The env-gated attribution modes (default OFF = the committed
    // all-or-nothing remap; the receipt fields in cut-stats.json name
    // the mode). The GLOBAL frame ignores them: the projection's
    // records are already panel-frame — nothing to remap or re-derive.
    let mode = if frame == CutFrame::Global {
        eprintln!(
            "[cut] global frame: the committed remap flow over the IDENTITY \
             expressibility map (the territory's own canonical anchors, \
             derived from the global panel; the census/cache are global-frame)"
        );
        RemapMode::GlobalFrame
    } else if std::env::var("IMPG_CUT_LOCAL_DERIVE").ok().as_deref() == Some("1") {
        eprintln!("[cut] local-derive remap: IMPG_CUT_LOCAL_DERIVE=1 (route globally, attribute locally - the committed local MEM derivation over the globally routed reads)");
        RemapMode::LocalDerive
    } else if std::env::var("IMPG_CUT_PARTIAL_REMAP").ok().as_deref() == Some("1") {
        eprintln!("[cut] partial remap: IMPG_CUT_PARTIAL_REMAP=1 (locally-expressible subwalk attribution)");
        RemapMode::PartialSubwalk
    } else {
        RemapMode::Strict
    };
    let per_file: Vec<io::Result<FileCut>> = shard_paths
        .par_iter()
        .map(|path| {
            scan_shard_file(
                path,
                node_map_ref,
                territory_nodes_ref,
                mode,
                if mode == RemapMode::LocalDerive {
                    Some(&panel)
                } else {
                    None
                },
            )
        })
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
        reads_routed: u64,
        reads_unrouted: u64,
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
        reads_routed: 0,
        reads_unrouted: 0,
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
        merged.reads_routed += cut.reads_routed;
        merged.reads_unrouted += cut.reads_unrouted;
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
    // (THE GLOBAL FRAME's pair-bins convention: the contig comparison
    // runs in the panel's own frame — the translated component's
    // contig suffix; the committed local frame keeps its own.)
    let component_contig = component_contig_source
        .rsplit('#')
        .next()
        .unwrap_or(&component_contig_source)
        .to_string();
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
        frame: match frame {
            CutFrame::Local => "local".into(),
            CutFrame::Global => "global".into(),
        },
        census_profile: match frame {
            CutFrame::Local => None,
            CutFrame::Global => Some(census_profile.label().into()),
        },
        record_origin,
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
        partial_remap: mode == RemapMode::PartialSubwalk,
        local_derive: mode == RemapMode::LocalDerive,
        global_frame: mode == RemapMode::GlobalFrame,
        reads_routed: merged.reads_routed,
        reads_unrouted: merged.reads_unrouted,
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

// ---------------------------------------------------------------------------
// cut-batch: THE WHOLE-GENOME FLEET CUT (the amortized-load shape - the
// owner's whole-genome 35x calling run). ONE global syng + AGC load
// serving a BATCH of territories per process: every routed read's
// records are derived by the committed local MEM accessor itself
// (IMPG_CUT_LOCAL_DERIVE semantics - route globally, attribute locally -
// REQUIRED for this mode), and each territory's census/cache/reads/
// stats outputs are byte-for-byte the single cut's (the shared code is
// the single cut's own, sliced verbatim; the runtime gate is
// byte-identity of a batch territory against the cut3 single-cut
// receipt). Receipt-side example only: no src/ changes, no product or
// scoreboard machinery.
// ---------------------------------------------------------------------------

/// One territory row of the batch manifest (JSONL: name, panel, axis,
/// bed_directory, component, sources, out_dir - the committed
/// substrate inputs of the single cut).
#[derive(serde::Deserialize)]
struct CutBatchSpec {
    name: String,
    panel: String,
    axis: PathBuf,
    bed_directory: PathBuf,
    component: String,
    sources: PathBuf,
    out_dir: PathBuf,
}

/// One GLOBAL-frame batch territory's manifest row (the shared-index
/// cut lanes): the single global cut's own substrate inputs — no
/// locality panel, no locus.fa (conversion-not-materialization) — plus
/// the substrate's stated content origin.
#[derive(serde::Deserialize, Clone)]
struct CutBatchGlobalSpec {
    name: String,
    axis: PathBuf,
    bed_directory: PathBuf,
    component: String,
    out_dir: PathBuf,
    /// The substrate's stated content origin (the single global
    /// cut's --record-origin d; 1 = the L1-L4 one-base-later
    /// convention, 0 = the pilot's exact headers).
    #[serde(default)]
    record_origin: u64,
    /// THE COMBINED CUT+CALL ROW (cut-call-batch-global only; ignored
    /// by cut-batch-global): when present, the SAME process calls the
    /// territory after its cut, over the ONE loaded panel — the
    /// shared-index cut+call seam (the call's per-process inputs pole
    /// — the panel's own load — vanishes into the shared load).
    #[serde(default)]
    call: Option<CutBatchGlobalCall>,
}

/// The per-territory call row of the combined cut+call manifest: the
/// call-diplotypes CLI's own inputs (the cut's census/cache are the
/// cut's own outputs, bound by construction).
#[derive(serde::Deserialize, Clone)]
struct CutBatchGlobalCall {
    /// The sample FASTQ (read identity + the quality-byte assertion).
    reads: PathBuf,
    /// The territory's locality map directory (window<N>.gfa.map.json
    /// — the committed per-territory substrate).
    partition_graphs: PathBuf,
    /// The scored loci (window ids of the component, comma-separated,
    /// distinct and sorted — the committed --loci form).
    loci: String,
    /// The call's output directory (calls.jsonl + the instrument
    /// receipts; refused if it exists — the CLI's own discipline).
    out_dir: PathBuf,
}

/// One loaded GLOBAL-frame batch territory: the shared panel is held
/// ONCE by the batch; per territory only the committed component
/// territory index and the identity expressibility map are loaded.
struct GlobalBatchTerritory {
    spec: CutBatchGlobalSpec,
    territory: TerritoryIndex,
    territory_nodes_abs: HashSet<u32>,
    /// The identity expressibility map (the single global cut's own
    /// node_map: the territory's own canonical anchors, every node
    /// remapped to itself — the evidence filter, in the panel's own
    /// frame).
    node_map: HashMap<i32, i32>,
    node_map_pairs: usize,
    component_contig_source: String,
}

/// One loaded batch territory: the locality panel, the committed
/// component territory index, and the global<->locality node map.
struct BatchTerritory {
    spec: CutBatchSpec,
    panel: SyngIndex,
    k: u64,
    territory: TerritoryIndex,
    territory_nodes_abs: HashSet<u32>,
    node_map: HashMap<i32, i32>,
    node_map_pairs: usize,
    global_only_anchor_total: usize,
    offset_variants: BTreeMap<String, u64>,
}

impl FileCut {
    fn empty() -> Self {
        FileCut {
            reads_scanned: 0,
            records_scanned: 0,
            records_mappable: 0,
            records_unmappable: 0,
            reads_with_unmappable: 0,
            reads_routed: 0,
            reads_unrouted: 0,
            reads: Vec::new(),
            pattern_mult: BTreeMap::new(),
            key_tokens: Vec::new(),
            key_index: HashMap::new(),
            key_reads: Vec::new(),
            read_records: Vec::new(),
        }
    }
}

    struct MergedBatch {
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
        reads_routed: u64,
        reads_unrouted: u64,
    }

impl MergedBatch {
    fn new() -> Self {
        MergedBatch {
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
            reads_routed: 0,
            reads_unrouted: 0,
        }
    }
}

    fn batch_intern_merged(merged: &mut MergedBatch, tokens: &[u64]) -> u32 {
        if let Some(&index) = merged.key_index.get(tokens) {
            return index;
        }
        let index = merged.key_tokens.len() as u32;
        merged.key_tokens.push(tokens.to_vec());
        merged.key_reads.push(Vec::new());
        merged.key_index.insert(tokens.to_vec(), index);
        index
    }

/// Merge one territory's per-file contributions in file order (the
/// single cut's merge loop, verbatim); returns the merged state.
fn batch_merge_file(mut merged: MergedBatch, files: Vec<FileCut>) -> MergedBatch {
    for cut in files {
        merged.reads_scanned += cut.reads_scanned;
        merged.records_scanned += cut.records_scanned;
        merged.records_mappable += cut.records_mappable;
        merged.records_unmappable += cut.records_unmappable;
        merged.reads_with_unmappable += cut.reads_with_unmappable;
        merged.reads_routed += cut.reads_routed;
        merged.reads_unrouted += cut.reads_unrouted;
        let base = merged.reads.len();
        for (index, seq) in cut.reads.into_iter().enumerate() {
            merged.reads.push(seq);
            let read_index = base + index;
            let records = cut.read_records[index].clone();
            let mut rebased: Vec<(u32, OwnWalk)> = Vec::with_capacity(records.len());
            for (local_key, walk) in records {
                let tokens = cut.key_tokens[local_key as usize].clone();
                let new_key = batch_intern_merged(&mut merged, &tokens);
                merged.key_reads[new_key as usize].push(read_index);
                rebased.push((new_key, walk));
            }
            merged.read_records.push(rebased);
        }
        for (tokens, mult) in cut.pattern_mult {
            *merged.pattern_mult.entry(tokens).or_insert(0) += mult;
        }
    }
    merged
}

/// The locus.fa records (the single cut's loader, verbatim).
fn batch_load_sources(sources_path: &Path) -> io::Result<HashMap<String, Vec<u8>>> {
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
    Ok(sources)
}

/// The per-path global<->locality node correspondence (the single
/// cut's construction, verbatim).
#[allow(clippy::type_complexity)]
fn batch_build_node_map(
    global: &SyngIndex,
    panel: &SyngIndex,
    agc: &AgcIndex,
    sources: &HashMap<String, Vec<u8>>,
) -> io::Result<(
    HashMap<i32, i32>,
    usize,
    usize,
    BTreeMap<String, u64>,
)> {
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
    for (header, record) in sources {
        let (global_name, s, e) = parse_locus_header(header)?;
        ensure(
            (e - s) as usize == record.len(),
            &format!("locus record {header} length differs from its header interval"),
        )?;
        // The offset candidates cover both committed extraction
        // conventions (the generic exact-bed intervals and the pilot's
        // one-base-earlier records) under the AGC range fetch's own
        // 1-based-inclusive coordinate system; sequence identity decides.
        let mut delta: Option<i64> = None;
        for candidate in [0i64, -1, 1] {
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
    Ok((
        node_map,
        node_map_pairs,
        global_only_anchor_total,
        offset_variants,
    ))
}

/// One territory's emission tail (the single cut's, verbatim): the
/// derive cache, the cut reads FASTA, the cut census, the stats.
#[allow(clippy::too_many_arguments)]
fn batch_emit_territory(
    t: BatchTerritory,
    files: Vec<FileCut>,
    load_seconds: f64,
    scan_seconds: f64,
    started: &Instant,
) -> io::Result<serde_json::Value> {
    let BatchTerritory {
        spec,
        panel,
        k,
        territory,
        territory_nodes_abs: _,
        node_map: _,
        node_map_pairs,
        global_only_anchor_total: _,
        offset_variants,
    } = t;
    let out_dir = spec.out_dir.as_path();
    let component = spec.component.as_str();
    let mode = RemapMode::LocalDerive;
    let merged = batch_merge_file(MergedBatch::new(), files);
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
        frame: "local-batch".into(),
        census_profile: None,
        record_origin: 0,
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
        partial_remap: mode == RemapMode::PartialSubwalk,
        local_derive: mode == RemapMode::LocalDerive,
        global_frame: mode == RemapMode::GlobalFrame,
        reads_routed: merged.reads_routed,
        reads_unrouted: merged.reads_unrouted,
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
    Ok(serde_json::json!({
        "name": spec.name,
        "reads_kept": stats.reads_kept,
        "routed_patterns": stats.routed_patterns,
        "census_mult": stats.census_mult,
        "load_seconds": load_seconds,
        "scan_seconds": scan_seconds,
        "census_seconds": stats.census_seconds,
        "wall_seconds": stats.wall_seconds,
    }))
}

/// THE BATCH CUT: one global syng + AGC load, one shard scan, many
/// territories' cut outputs (each byte-identical to the single cut's).
fn run_cut_batch(
    global_panel_prefix: &str,
    agc_path: &Path,
    shards_dir: &Path,
    manifest_path: &Path,
) -> io::Result<()> {
    ensure(
        std::env::var("IMPG_CUT_LOCAL_DERIVE").ok().as_deref() == Some("1"),
        "the batch cut is the local-derive fleet mode (route globally, attribute locally); set IMPG_CUT_LOCAL_DERIVE=1",
    )?;
    let mut specs: Vec<CutBatchSpec> = Vec::new();
    for line in BufReader::new(std::fs::File::open(manifest_path)?).lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        specs.push(
            serde_json::from_str(&line)
                .map_err(|e| invalid(format!("manifest row: {e}")))?,
        );
    }
    ensure(!specs.is_empty(), "empty batch manifest")?;
    ensure(specs.len() <= 65535, "batch wider than the u16 territory index")?;
    for spec in &specs {
        std::fs::create_dir_all(&spec.out_dir)?;
        for name in [
            "cut-census.jsonl",
            "cut-cache.bin",
            "cut-reads.fasta",
            "cut-stats.json",
        ] {
            ensure(
                !spec.out_dir.join(name).exists(),
                &format!(
                    "cut output {name} already exists in {}",
                    spec.out_dir.display()
                ),
            )?;
        }
    }
    let started = Instant::now();
    eprintln!(
        "[cut-batch] loading the global panel + AGC once for {} territories",
        specs.len()
    );
    let load_started = Instant::now();
    let global = SyngIndex::load(global_panel_prefix, SyncmerParams::default())?;
    let agc = AgcIndex::build_from_files(&[agc_path.display().to_string()])?;
    let global_k = global.syncmer_length_bp() as u64;

    // Per-territory setup (parallel; every locality panel is small).
    let setups: Vec<io::Result<BatchTerritory>> = specs
        .into_par_iter()
        .map(|spec| {
            let panel = SyngIndex::load(&spec.panel, SyncmerParams::default())?;
            let k = panel.syncmer_length_bp() as u64;
            ensure(
                global_k == k,
                "global and locality syncmer lengths differ",
            )?;
            let axis: AxisFile = serde_json::from_str(&std::fs::read_to_string(&spec.axis)?)
                .map_err(|e| invalid(format!("axis parse: {e}")))?;
            let universe = load_universe_component(&axis, &spec.bed_directory, &spec.component)?;
            eprintln!("[cut-batch] {} universe rows: {}", spec.name, universe.len());
            let sources = batch_load_sources(&spec.sources)?;
            let fetch_lane = lane_fetch_over(&sources);
            let territory = build_territory_index(&panel, &universe, &fetch_lane, k)?;
            eprintln!(
                "[cut-batch] {} territory: {} rows, {} entries",
                spec.name,
                territory.territories.len(),
                territory.entries.len()
            );
            let territory_nodes_abs: HashSet<u32> = territory
                .entries
                .iter()
                .map(|&(node, _, _)| node.unsigned_abs())
                .collect();
            let (node_map, node_map_pairs, global_only_anchor_total, offset_variants) =
                batch_build_node_map(&global, &panel, &agc, &sources)?;
            eprintln!(
                "[cut-batch] {} node map: {} pairs, {} global-only anchors",
                spec.name, node_map_pairs, global_only_anchor_total
            );
            Ok(BatchTerritory {
                spec,
                panel,
                k,
                territory,
                territory_nodes_abs,
                node_map,
                node_map_pairs,
                global_only_anchor_total,
                offset_variants,
            })
        })
        .collect();
    let territories: Vec<BatchTerritory> = setups.into_iter().collect::<io::Result<_>>()?;
    let n_terr = territories.len();
    // The per-territory fail-closed isolation bookkeeping: the emissions
    // loop uses this to name and ISOLATE a failing territory (the
    // 7479-class bind ambiguity fails ONE territory's cut, never the
    // batch) - the batch's other territories keep their cuts.
    let emit_specs: Vec<(String, std::path::PathBuf)> = territories
        .iter()
        .map(|t| (t.spec.name.clone(), t.spec.out_dir.clone()))
        .collect();

    // The inverted routing index: global node -> territory ids. The
    // routing signal is exactly the single cut's (a read is a candidate
    // for a territory iff any of its records' walk nodes is a key of
    // that territory's node map).
    let mut routing: HashMap<i32, Vec<u16>> = HashMap::new();
    for (index, t) in territories.iter().enumerate() {
        for &gnode in t.node_map.keys() {
            routing.entry(gnode).or_default().push(index as u16);
        }
    }
    eprintln!(
        "[cut-batch] routing index: {} global nodes over {} territories",
        routing.len(),
        n_terr
    );
    let load_seconds = load_started.elapsed().as_secs_f64();

    // The ONE shard scan (parallel over files; per read, every routed
    // territory's local derivation in parallel over reads; the
    // per-(file, territory) states accumulate in read order exactly as
    // the single cut's per-file state does).
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
    let routing_ref = &routing;
    let territories_ref = &territories;
    let scan_started = Instant::now();
    let per_file: Vec<io::Result<HashMap<u16, FileCut>>> = shard_paths
        .par_iter()
        .map(|path| {
            let mut per: HashMap<u16, FileCut> = HashMap::new();
            let mut file_reads = 0u64;
            let mut file_records = 0u64;
            for batch in shard_packet_stream(path, 1024)? {
                let packets = batch?;
                file_reads += packets.len() as u64;
                let decoded: Vec<DecodedRead> = packets
                    .par_iter()
                    .map(|packet| decode_read_packet(packet))
                    .collect::<io::Result<_>>()?;
                file_records += decoded
                    .iter()
                    .map(|read| read.records.len() as u64)
                    .sum::<u64>();
                // Per read: the routed territory set, then every routed
                // territory's LOCAL derivation (parallel over reads).
                let processed: Vec<
                    io::Result<
                        Option<(
                            Vec<u8>,
                            Vec<(u16, Vec<Option<(Vec<u64>, OwnWalk, bool)>>)>,
                        )>,
                    >,
                > = decoded
                    .into_par_iter()
                    .map(|read| {
                        let mut routed: Vec<u16> = Vec::new();
                        for (_, walk) in &read.records {
                            for &(node, _) in walk {
                                if let Some(list) = routing_ref.get(&node) {
                                    routed.extend_from_slice(list);
                                }
                            }
                        }
                        routed.sort_unstable();
                        routed.dedup();
                        if routed.is_empty() {
                            return Ok(None);
                        }
                        let mut out = Vec::with_capacity(routed.len());
                        for &ti in &routed {
                            let t = &territories_ref[ti as usize];
                            // ROUTE GLOBALLY, ATTRIBUTE LOCALLY (the
                            // committed local MEM derivation).
                            let local = derive_read(&t.panel, &read.seq)?;
                            let records: Vec<Option<(Vec<u64>, OwnWalk, bool)>> = local
                                .into_iter()
                                .map(|(tokens, walk)| {
                                    let touch = walk.iter().any(|&(node, _)| {
                                        t.territory_nodes_abs
                                            .contains(&node.unsigned_abs())
                                    });
                                    Some((tokens, walk, touch))
                                })
                                .collect();
                            out.push((ti, records));
                        }
                        Ok(Some((read.seq, out)))
                    })
                    .collect();
                // Fold in read order (the single cut's per-read
                // accumulation, per territory).
                for result in processed {
                    let Some((seq, per_terr)) = result? else {
                        continue;
                    };
                    for (ti, records) in per_terr {
                        let cut = per.entry(ti).or_insert_with(FileCut::empty);
                        cut.reads_routed += 1;
                        cut.records_mappable += records.len() as u64;
                        let mut kept: Vec<(u32, OwnWalk)> = Vec::new();
                        let mut any_touch = false;
                        for record in &records {
                            if let Some((tokens, walk, touch)) = record {
                                if *touch {
                                    any_touch = true;
                                    *cut.pattern_mult.entry(tokens.clone()).or_insert(0) += 1;
                                }
                                let key = cut.intern_key(tokens);
                                kept.push((key, walk.clone()));
                            }
                        }
                        if any_touch {
                            let read_index = cut.reads.len();
                            for (key, _walk) in &kept {
                                cut.key_reads[*key as usize].push(read_index);
                            }
                            cut.reads.push(seq.clone());
                            cut.read_records.push(std::mem::take(&mut kept));
                        }
                    }
                }
            }
            // The file's full-scan accounting, per territory (the single
            // cut counts every scanned read/record in every per-file
            // state; the batch reproduces that per territory).
            for index in 0..n_terr {
                let cut = per.entry(index as u16).or_insert_with(FileCut::empty);
                cut.reads_scanned = file_reads;
                cut.records_scanned = file_records;
                cut.reads_unrouted = file_reads - cut.reads_routed;
                cut.reads_with_unmappable = cut.reads_unrouted;
            }
            Ok(per)
        })
        .collect();
    let scan_seconds = scan_started.elapsed().as_secs_f64();
    let per_file: Vec<HashMap<u16, FileCut>> =
        per_file.into_iter().collect::<io::Result<_>>()?;
    eprintln!(
        "[cut-batch] scan done in {scan_seconds:.0}s; emitting {} territories",
        n_terr
    );

    // Per territory: the per-file states in file order, then the single
    // cut's emission tail (parallel over territories).
    let mut per_terr_files: Vec<Vec<FileCut>> = (0..n_terr).map(|_| Vec::new()).collect();
    for per in per_file {
        for (ti, cut) in per {
            per_terr_files[ti as usize].push(cut);
        }
    }
    let emissions: Vec<io::Result<serde_json::Value>> = territories
        .into_par_iter()
        .zip(per_terr_files)
        .map(|(t, files)| {
            batch_emit_territory(t, files, load_seconds, scan_seconds, &started)
        })
        .collect();
    let mut n_emit_failed = 0usize;
    for ((name, out_dir), result) in emit_specs.into_iter().zip(emissions) {
        match result {
            Ok(summary) => eprintln!("[cut-batch] {}: {}", summary["name"], summary),
            Err(e) => {
                // FAIL-CLOSED, ISOLATED: the territory is refused (its
                // partial outputs removed, the failure named on disk);
                // the batch's other territories proceed.
                n_emit_failed += 1;
                let _ = std::fs::write(out_dir.join("CUT-FAILED.txt"), format!("{e}\n"));
                for f in [
                    "cut-census.jsonl",
                    "cut-cache.bin",
                    "cut-reads.fasta",
                    "cut-stats.json",
                ] {
                    let _ = std::fs::remove_file(out_dir.join(f));
                }
                eprintln!(
                    "[cut-batch] {}: EMISSION FAILED (fail-closed, isolated): {e}",
                    name
                );
            }
        }
    }
    eprintln!(
        "[cut-batch] DONE: {} territories in {:.0}s (peak RSS {}kB, {} emission-failures)",
        n_terr,
        started.elapsed().as_secs_f64(),
        rss_now_kb(),
        n_emit_failed
    );
    Ok(())
}


// ---------------------------------------------------------------------------
// cut-batch-global: THE SHARED-INDEX CUT LANES (lever 2 — the 325s
// global GBWT load pole amortized over a batch of global-frame
// territories; one shard scan serves every territory). The per-file,
// per-read fold and the per-territory merge reproduce the single
// global cut's own accumulation order exactly (fold in read order per
// file; merge in file order per territory), so each territory's
// cut-census.jsonl / cut-cache.bin / cut-reads.fasta are byte-
// identical to the single cut's — the gate.
// ---------------------------------------------------------------------------
fn cut_batch_global_parse_manifest(
    manifest_path: &Path,
) -> io::Result<Vec<CutBatchGlobalSpec>> {
    let mut specs: Vec<CutBatchGlobalSpec> = Vec::new();
    for line in BufReader::new(std::fs::File::open(manifest_path)?).lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        specs.push(
            serde_json::from_str(&line)
                .map_err(|e| invalid(format!("manifest row: {e}")))?,
        );
    }
    ensure(!specs.is_empty(), "empty batch manifest")?;
    ensure(specs.len() <= 65535, "batch wider than the u16 territory index")?;
    for spec in &specs {
        std::fs::create_dir_all(&spec.out_dir)?;
        for name in [
            "cut-census.jsonl",
            "cut-cache.bin",
            "cut-reads.fasta",
            "cut-stats.json",
        ] {
            ensure(
                !spec.out_dir.join(name).exists(),
                &format!(
                    "cut output {name} already exists in {}",
                    spec.out_dir.display()
                ),
            )?;
        }
    }
    Ok(specs)
}

fn run_cut_batch_global(
    global_panel_prefix: &str,
    agc_path: &Path,
    shards_dir: &Path,
    manifest_path: &Path,
) -> io::Result<()> {
    let specs = cut_batch_global_parse_manifest(manifest_path)?;
    eprintln!(
        "[cut-batch-global] loading the global panel + AGC once for {} territories",
        specs.len()
    );
    let load_started = Instant::now();
    let global = SyngIndex::load(global_panel_prefix, SyncmerParams::default())?;
    let agc = AgcIndex::build_from_files(&[agc_path.display().to_string()])?;
    let load_seconds = load_started.elapsed().as_secs_f64();
    cut_batch_global_with(&global, &agc, shards_dir, specs, load_seconds)
}

/// THE SHARED-INDEX CORE: the batch cut over a CALLER-HELD loaded
/// panel + AGC (the load wall measured by the caller; the body is
/// run_cut_batch_global's own verbatim from the path-length map down)
/// — the combined cut+call process's cut phase, and the single-load
/// seam the call phase injects the same panel through.
fn cut_batch_global_with(
    global: &SyngIndex,
    agc: &AgcIndex,
    shards_dir: &Path,
    specs: Vec<CutBatchGlobalSpec>,
    load_seconds: f64,
) -> io::Result<()> {
    let started = Instant::now();
    let k = global.syncmer_length_bp() as u64;
    // The shared AGC lane fetch (the single global cut's own closure,
    // verbatim): fetch(A, B) is the 0-based [A, B) slice, clipped at
    // the path end, empty when the interval is empty.
    let path_lengths: HashMap<&str, u64> = global
        .name_map
        .path_to_name
        .iter()
        .zip(&global.name_map.path_to_length)
        .map(|(name, &length)| (name.as_str(), length))
        .collect();
    let agc_ref = agc;
    let fetch_lane =
        move |name: &str, lo: u64, hi: u64| -> io::Result<Vec<u8>> {
            let path_len = *path_lengths
                .get(name)
                .ok_or_else(|| invalid("universe path absent from the global panel"))?;
            let hi = hi.min(path_len);
            if lo >= hi {
                return Ok(Vec::new());
            }
            Ok(agc_ref
                .fetch_sequence(name, lo as i32, hi as i32)
                .map_err(|e| invalid(format!("AGC fetch {name}: {e}")))?)
        };

    // Per-territory setup (parallel; every territory's own universe,
    // territory index, and identity expressibility map — the single
    // global cut's Global-frame branch, verbatim).
    let global_ref = global;
    let setups: Vec<io::Result<GlobalBatchTerritory>> = specs
        .into_par_iter()
        .map(|spec| {
            let axis: AxisFile = serde_json::from_str(&std::fs::read_to_string(&spec.axis)?)
                .map_err(|e| invalid(format!("axis parse: {e}")))?;
            let universe =
                load_universe_component(&axis, &spec.bed_directory, &spec.component)?;
            eprintln!(
                "[cut-batch-global] {} universe rows: {}",
                spec.name,
                universe.len()
            );
            let record_origin = spec.record_origin;
            let universe_global: Vec<UniverseRow> = universe
                .iter()
                .map(|row| {
                    let translated = territory_frame::translate_row(
                        &row.path_name,
                        row.start,
                        row.end,
                        record_origin,
                    )?;
                    Ok::<_, io::Error>(UniverseRow {
                        partition: row.partition,
                        path_name: translated.path_name,
                        start: translated.start,
                        end: translated.end,
                    })
                })
                .collect::<io::Result<_>>()?;
            let territory = build_territory_index(global_ref, &universe_global, &fetch_lane, k)?;
            eprintln!(
                "[cut-batch-global] {} territory (global frame): {} rows, {} entries",
                spec.name,
                territory.territories.len(),
                territory.entries.len()
            );
            let territory_nodes_abs: HashSet<u32> = territory
                .entries
                .iter()
                .map(|&(node, _, _)| node.unsigned_abs())
                .collect();
            let component_contig_source =
                territory_frame::split_local_name(&spec.component)?.0.to_string();
            // THE IDENTITY EXPRESSIBILITY MAP (the single global cut's
            // own derivation, verbatim): the canonical anchor walk of
            // every record's content — fetched from the AGC at the
            // record's CONTENT origin start + d — every node remapped
            // to itself.
            let mut node_map: HashMap<i32, i32> = HashMap::new();
            let mut expressibility_pairs = 0usize;
            let mut record_names: BTreeSet<&str> = BTreeSet::new();
            for row in &universe {
                record_names.insert(row.path_name.as_str());
            }
            for header in &record_names {
                let (global_name, start, end) = parse_locus_header(header)?;
                let length = (end - start) as usize;
                let content_start = start
                    .checked_add(record_origin)
                    .ok_or_else(|| invalid(format!(
                        "record {header}: the content origin start+{record_origin} overflows"
                    )))?;
                let content = agc_ref
                    .fetch_sequence(
                        &global_name,
                        content_start as i32,
                        (content_start + length as u64) as i32,
                    )
                    .map_err(|e| invalid(format!("AGC fetch {global_name}: {e}")))?;
                ensure(
                    content.len() == length,
                    &format!(
                        "record {header}: the content fetch at origin +{record_origin} \
                         returned {} bytes, not the record's {length} (the content origin \
                         leaves the path; state the substrate's own convention)",
                        content.len()
                    ),
                )?;
                let walk = mem_records::canonical_anchor_walk(global_ref, &content)?;
                for &(node, _position) in &walk {
                    match node_map.entry(node) {
                        std::collections::hash_map::Entry::Occupied(o) => {
                            ensure(*o.get() == node, "conflicting identity entry")?;
                        }
                        std::collections::hash_map::Entry::Vacant(v) => {
                            v.insert(node);
                            expressibility_pairs += 1;
                        }
                    }
                }
            }
            eprintln!(
                "[cut-batch-global] {} identity expressibility map over {} records: \
                 {} canonical anchors",
                spec.name,
                record_names.len(),
                expressibility_pairs
            );
            Ok(GlobalBatchTerritory {
                spec,
                territory,
                territory_nodes_abs,
                node_map,
                node_map_pairs: expressibility_pairs,
                component_contig_source,
            })
        })
        .collect();
    let territories: Vec<GlobalBatchTerritory> =
        setups.into_iter().collect::<io::Result<_>>()?;
    let n_terr = territories.len();
    let emit_specs: Vec<(String, std::path::PathBuf)> = territories
        .iter()
        .map(|t| (t.spec.name.clone(), t.spec.out_dir.clone()))
        .collect();

    // The inverted routing index: global node -> territory ids. A read
    // is a candidate for a territory iff any of its records' walk
    // nodes is one of the territory's expressibility keys — the exact
    // completeness argument of the local batch's routing index (an
    // expressible record has EVERY node in the map, so a record that
    // touches has a key node; a read no territory routes contributes
    // only inexpressible records to every territory, exactly as the
    // single cut's all-or-nothing remap drops them).
    let mut routing: HashMap<i32, Vec<u16>> = HashMap::new();
    for (index, t) in territories.iter().enumerate() {
        for &gnode in t.node_map.keys() {
            routing.entry(gnode).or_default().push(index as u16);
        }
    }
    eprintln!(
        "[cut-batch-global] routing index: {} global nodes over {} territories",
        routing.len(),
        n_terr
    );

    // The ONE shard scan (parallel over files; per read, every routed
    // territory's committed identity remap; the per-(file, territory)
    // states accumulate in read order exactly as the single cut's
    // per-file state does).
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
    let routing_ref = &routing;
    let territories_ref = &territories;
    let scan_started = Instant::now();
    let per_file: Vec<io::Result<HashMap<u16, FileCut>>> = shard_paths
        .par_iter()
        .map(|path| {
            let mut per: HashMap<u16, FileCut> = HashMap::new();
            let mut file_reads = 0u64;
            let mut file_records = 0u64;
            // Per territory: the routed reads' record total (for the
            // file-end unmappable remainder — every record of every
            // UNROUTED read is inexpressible for the territory, exactly
            // the single cut's all-or-nothing accounting).
            let mut routed_records: HashMap<u16, u64> = HashMap::new();
            for batch in shard_packet_stream(path, 1024)? {
                let packets = batch?;
                file_reads += packets.len() as u64;
                let decoded: Vec<DecodedRead> = packets
                    .par_iter()
                    .map(|packet| decode_read_packet(packet))
                    .collect::<io::Result<_>>()?;
                file_records += decoded
                    .iter()
                    .map(|read| read.records.len() as u64)
                    .sum::<u64>();
                // Per read: the routed territory set, then every
                // routed territory's committed identity remap (parallel
                // over reads).
                let processed: Vec<
                    io::Result<
                        Option<(
                            Vec<u8>,
                            Vec<(u16, Vec<Option<(Vec<u64>, OwnWalk, bool)>>)>,
                        )>,
                    >,
                > = decoded
                    .into_par_iter()
                    .map(|read| {
                        let mut routed: Vec<u16> = Vec::new();
                        for (_, walk) in &read.records {
                            for &(node, _) in walk {
                                if let Some(list) = routing_ref.get(&node) {
                                    routed.extend_from_slice(list);
                                }
                            }
                        }
                        routed.sort_unstable();
                        routed.dedup();
                        if routed.is_empty() {
                            return Ok(None);
                        }
                        let mut out = Vec::with_capacity(routed.len());
                        for &ti in &routed {
                            let t = &territories_ref[ti as usize];
                            let mut records: Vec<Option<(Vec<u64>, OwnWalk, bool)>> =
                                Vec::with_capacity(read.records.len());
                            cut_records(
                                &read.records,
                                &t.node_map,
                                &t.territory_nodes_abs,
                                false,
                                &mut records,
                            )?;
                            out.push((ti, records));
                        }
                        Ok(Some((read.seq, out)))
                    })
                    .collect();
                // Fold in read order (the single cut's per-read
                // accumulation, per territory).
                for result in processed {
                    let Some((seq, per_terr)) = result? else {
                        continue;
                    };
                    for (ti, records) in per_terr {
                        let cut = per.entry(ti).or_insert_with(FileCut::empty);
                        *routed_records.entry(ti).or_insert(0) += records.len() as u64;
                        let mut kept: Vec<(u32, OwnWalk)> = Vec::new();
                        let mut any_touch = false;
                        for record in &records {
                            if let Some((tokens, walk, touch)) = record {
                                cut.records_mappable += 1;
                                if *touch {
                                    any_touch = true;
                                    *cut.pattern_mult.entry(tokens.clone()).or_insert(0) += 1;
                                }
                                let key = cut.intern_key(tokens);
                                kept.push((key, walk.clone()));
                            } else {
                                cut.records_unmappable += 1;
                            }
                        }
                        if any_touch {
                            let read_index = cut.reads.len();
                            for (key, _walk) in &kept {
                                cut.key_reads[*key as usize].push(read_index);
                            }
                            cut.reads.push(seq.clone());
                            cut.read_records.push(std::mem::take(&mut kept));
                        }
                    }
                }
            }
            // The file's full-scan accounting, per territory (the single
            // global cut counts every scanned read/record in its
            // per-file state; the batch reproduces that per territory,
            // with the unrouted reads' records counted unmappable —
            // every such record is inexpressible, the all-or-nothing
            // rule's own accounting).
            for index in 0..n_terr {
                let ti = index as u16;
                let cut = per.entry(ti).or_insert_with(FileCut::empty);
                cut.reads_scanned = file_reads;
                cut.records_scanned = file_records;
                let routed = routed_records.get(&ti).copied().unwrap_or(0);
                cut.records_unmappable += file_records - routed;
            }
            Ok(per)
        })
        .collect();
    let scan_seconds = scan_started.elapsed().as_secs_f64();
    let per_file: Vec<HashMap<u16, FileCut>> =
        per_file.into_iter().collect::<io::Result<_>>()?;
    eprintln!(
        "[cut-batch-global] scan done in {scan_seconds:.0}s; emitting {n_terr} territories"
    );

    // Per territory: the per-file states in file order, then the single
    // global cut's emission tail (parallel over territories).
    let mut per_terr_files: Vec<Vec<FileCut>> = (0..n_terr).map(|_| Vec::new()).collect();
    for per in per_file {
        for (ti, cut) in per {
            per_terr_files[ti as usize].push(cut);
        }
    }
    let emissions: Vec<io::Result<serde_json::Value>> = territories
        .into_par_iter()
        .zip(per_terr_files)
        .map(|(t, files)| {
            batch_global_emit_territory(&global, t, files, load_seconds, scan_seconds, &started)
        })
        .collect();
    let mut n_emit_failed = 0usize;
    for ((name, out_dir), result) in emit_specs.into_iter().zip(emissions) {
        match result {
            Ok(summary) => eprintln!("[cut-batch-global] {}: {}", summary["name"], summary),
            Err(e) => {
                // FAIL-CLOSED, ISOLATED (the local batch's own rule): the
                // territory is refused, its partial outputs removed, the
                // failure named on disk; the batch's other territories
                // proceed.
                n_emit_failed += 1;
                let _ = std::fs::write(out_dir.join("CUT-FAILED.txt"), format!("{e}\n"));
                for f in [
                    "cut-census.jsonl",
                    "cut-cache.bin",
                    "cut-reads.fasta",
                    "cut-stats.json",
                ] {
                    let _ = std::fs::remove_file(out_dir.join(f));
                }
                eprintln!(
                    "[cut-batch-global] {}: EMISSION FAILED (fail-closed, isolated): {e}",
                    name
                );
            }
        }
    }
    eprintln!(
        "[cut-batch-global] DONE: {} territories in {:.0}s (peak RSS {}kB, {} emission-failures)",
        n_terr,
        started.elapsed().as_secs_f64(),
        rss_now_kb(),
        n_emit_failed
    );
    Ok(())
}

/// THE COMBINED CUT+CALL PROCESS (the shared-index seam): parse the
/// manifest, load the ONE global panel + AGC, run the batch cut over
/// them (cut_batch_global_with), then call every row carrying a `call`
/// object IN-PROCESS over the SAME loaded panel — the injected-index
/// form of call-diplotypes (realign_score::run_injected): the call's
/// per-process inputs pole (the panel's own load) vanishes into the
/// shared load. The per-call options reproduce the committed CLI's
/// construction verbatim (the instrument receipts + the product
/// emission under the output-model discipline, --selection Locality,
/// --emit-sequences, --global-frame), with --panel naming the same
/// global panel the process holds (the identity check binds it to
/// the routes graph exactly as the loaded path does).
fn run_cut_call_batch_global(
    global_panel_prefix: &str,
    agc_path: &Path,
    shards_dir: &Path,
    manifest_path: &Path,
    routes_dir: &Path,
    call_rss_budget_gib: f64,
) -> io::Result<()> {
    let specs = cut_batch_global_parse_manifest(manifest_path)?;
    ensure(
        specs.iter().any(|spec| spec.call.is_some()),
        "no manifest row carries a call object (cut-batch-global is the cut-only form)",
    )?;
    // The call rows, extracted before the cut phase consumes the
    // specs: (name, component, record_origin, the cut's out_dir, the
    // call row itself).
    let call_rows: Vec<(String, String, u64, PathBuf, CutBatchGlobalCall)> = specs
        .iter()
        .filter_map(|spec| {
            spec.call.as_ref().map(|call| {
                (
                    spec.name.clone(),
                    spec.component.clone(),
                    spec.record_origin,
                    spec.out_dir.clone(),
                    call.clone(),
                )
            })
        })
        .collect();
    eprintln!(
        "[cut-call-batch-global] loading the global panel + AGC once for {} territories ({} to call in-process)",
        specs.len(),
        call_rows.len()
    );
    let load_started = Instant::now();
    let global = SyngIndex::load(global_panel_prefix, SyncmerParams::default())?;
    let agc = AgcIndex::build_from_files(&[agc_path.display().to_string()])?;
    let load_seconds = load_started.elapsed().as_secs_f64();
    let cut_started = Instant::now();
    cut_batch_global_with(&global, &agc, shards_dir, specs, load_seconds)?;
    eprintln!(
        "[cut-call-batch-global] cut phase done in {:.0}s; calling {} territories over the loaded panel",
        cut_started.elapsed().as_secs(),
        call_rows.len()
    );
    let calls_started = Instant::now();
    let shared_panel: &SyngIndex = &global;
    let call_results = run_combined_call_rows(
        shared_panel,
        global_panel_prefix,
        routes_dir,
        &call_rows,
        call_rss_budget_gib,
        load_seconds,
        "cut-call-batch-global",
    );
    eprintln!(
        "[cut-call-batch-global] call phase done in {:.0}s ({} concurrent)",
        calls_started.elapsed().as_secs(),
        call_rows.len()
    );
    for result in call_results {
        result?;
    }
    Ok(())
}

/// THE COMBINED FORM'S CALL PHASE, shared verbatim by the combined
/// process and the forked children (the fork-after-load driver runs
/// the SAME construction over the SAME loaded panel): the batch's K
/// calls run CONCURRENTLY — one scoped thread per call row (the
/// manifest's own width, never tuned). The shared state is read-only
/// by construction (the loaded panel, the routes dir, the per-call
/// options; every mutable thing a call touches is its own: its out
/// dir, its LaneCache, its Sources handles, its receipt buffers), and
/// the ONE mutable process-global the call path held (the read
/// length) is PER-CALL state — so each call sees exactly the
/// committed fresh-process trajectory at any concurrency; the
/// outputs are the committed machinery's own, byte-identical by the
/// thread-schedule-invariance every gate has proven. Error
/// discipline: the serial form failed at the first failing row; this
/// form returns the FIRST error in ROW ORDER (no failure is
/// swallowed). (THE SS10.6.4 cheapest-named step, landed with the
/// combined process; extracted verbatim for the fork form.)
fn run_combined_call_rows(
    shared_panel: &SyngIndex,
    global_panel_prefix: &str,
    routes_dir: &Path,
    call_rows: &[(String, String, u64, PathBuf, CutBatchGlobalCall)],
    call_rss_budget_gib: f64,
    load_seconds: f64,
    tag: &str,
) -> Vec<io::Result<()>> {
    std::thread::scope(|scope| {
        let handles: Vec<_> = call_rows
            .iter()
            .map(|(name, component, record_origin, cut_dir, call)| {
                let (name, component, cut_dir, call) =
                    (name.clone(), component.clone(), cut_dir.clone(), call.clone());
                let out_dir = call.out_dir.clone();
                scope.spawn(move || {
                    let call_started = Instant::now();
                    let result =
                        impg::genome_inference::with_output_model(&out_dir, realign_score::MODEL, || {
                            let options = realign_score::Options {
                                panel: global_panel_prefix.to_string(),
                                routes: routes_dir.to_path_buf(),
                                partition_graphs: call.partition_graphs.clone(),
                                component: component.clone(),
                                census: cut_dir.join("cut-census.jsonl"),
                                reads: call.reads.clone(),
                                derive_cache: cut_dir.join("cut-cache.bin"),
                                loci: call.loci.clone(),
                                out: out_dir.join("instrument-receipt.jsonl"),
                                exactness_out: out_dir.join("instrument-receipt.exactness.jsonl"),
                                exactness_sample: 100,
                                skeleton_out: Some(out_dir.join("instrument-receipt.skeleton.jsonl")),
                                anatomy_out: None,
                                frame_audit_only: false,
                                rss_budget_gib: call_rss_budget_gib,
                                selection: realign_score::Selection::Locality,
                                product_calls: Some(out_dir.join("calls.jsonl")),
                                truth_qv: Vec::new(),
                                product_sequences: true,
                                exclude_paths: Vec::new(),
                                flank_overlap: None,
                                region: None,
                                global_frame: true,
                                record_origin: *record_origin,
                                qv_biwfa: "genome/instrumented/qv-biwfa/target/release/qv-biwfa".into(),
                            };
                            realign_score::run_injected(options, Some(shared_panel))?;
                            let calls = out_dir.join("calls.jsonl");
                            let sequences = out_dir.join("calls.jsonl.sequences.jsonl");
                            Ok(serde_json::json!({
                                "calls": impg::genome_inference::reconstruction::fingerprint(&calls)?,
                                "calls_bytes": std::fs::metadata(&calls)?.len(),
                                "instrument_receipt": impg::genome_inference::reconstruction::fingerprint(
                                    &out_dir.join("instrument-receipt.jsonl"),
                                )?,
                                "sequences_emission": impg::genome_inference::reconstruction::fingerprint(
                                    &sequences,
                                )?,
                                "selection": "Locality",
                                "scope": "per-locus-diplotype-calls-with-model-internal-qual",
                                "shared_panel_load_seconds": load_seconds,
                            }))
                        });
                    eprintln!(
                        "[{tag}] {} called in {:.0}s ({})",
                        name,
                        call_started.elapsed().as_secs(),
                        if result.is_ok() { "ok" } else { "FAILED" },
                    );
                    result
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| match handle.join() {
                Ok(result) => result,
                Err(panic) => {
                    let message = panic
                        .downcast_ref::<&str>()
                        .map(|s| s.to_string())
                        .or_else(|| panic.downcast_ref::<String>().cloned())
                        .unwrap_or_else(|| "call thread panicked".to_string());
                    Err(io::Error::other(message))
                }
            })
            .collect()
    })
}

// ---------------------------------------------------------------------------
// THE FORK-AFTER-LOAD COW PANEL SHARE (SS10.6.9's lever): the smaps
// receipts and the fork driver. The parent is the ONE load owner (the
// receipted 241s load); every child inherits the loaded panel through
// copy-on-write and receipts its territory (cut + call) with the
// panel's pages staying SHARED while unwritten. The COW-hold gate is
// these smaps numbers: a child's PRIVATE growth over its whole
// cut+call is exactly the COW-break measure (SS10.6.9's named risk —
// post-load writes into the arenas break COW page by page).
// ---------------------------------------------------------------------------

/// One smaps_rollup reading (kB units, /proc's own fields).
#[derive(Default, Clone, Copy)]
struct ForkSmapsRollup {
    rss_kb: u64,
    pss_kb: u64,
    private_clean_kb: u64,
    private_dirty_kb: u64,
    shared_clean_kb: u64,
    shared_dirty_kb: u64,
}

impl ForkSmapsRollup {
    fn private_kb(&self) -> u64 {
        self.private_clean_kb + self.private_dirty_kb
    }
    fn shared_kb(&self) -> u64 {
        self.shared_clean_kb + self.shared_dirty_kb
    }
}

/// Parse the named kB field out of a /proc smaps* section.
fn fork_smaps_field(section: &str, name: &str) -> u64 {
    for line in section.lines() {
        let mut parts = line.split_whitespace();
        if parts.next() == Some(name) {
            if let Some(kb) = parts.next().and_then(|v| v.parse::<u64>().ok()) {
                return kb;
            }
        }
    }
    0
}

/// One big mapping of /proc/self/smaps (the C arenas and the mmap'd
/// sidecars), identified by its header line so a child's open and
/// close snapshots can be diffed mapping-by-mapping.
#[derive(Clone)]
struct ForkSmapsMapping {
    header: String,
    size_kb: u64,
    rss_kb: u64,
    pss_kb: u64,
    private_clean_kb: u64,
    private_dirty_kb: u64,
    shared_clean_kb: u64,
    shared_dirty_kb: u64,
    name: String,
}

impl ForkSmapsMapping {
    fn private_kb(&self) -> u64 {
        self.private_clean_kb + self.private_dirty_kb
    }
    fn shared_kb(&self) -> u64 {
        self.shared_clean_kb + self.shared_dirty_kb
    }
}

/// A smaps header line "lo-hi perm offset dev inode [path]": the
/// first token is two hex addresses joined by '-'.
fn fork_is_smaps_header(line: &str) -> bool {
    let Some((lo, hi)) = line.split_whitespace().next().unwrap_or("").split_once('-')
    else {
        return false;
    };
    u64::from_str_radix(lo, 16).is_ok() && u64::from_str_radix(hi, 16).is_ok()
}

/// Parse /proc/self/smaps into sections; keep every mapping whose
/// Size >= `min_size_kb` (the 12.13GB GBWT arena, the 4.31GB
/// KmerHash arena, the mmap'd sidecars — the shareable mass; smaller
/// mappings are summarized by the rollup alone).
fn fork_smaps_big_mappings(smaps: &str, min_size_kb: u64) -> Vec<ForkSmapsMapping> {
    let mut big = Vec::new();
    let mut header: Option<String> = None;
    let mut body: Vec<&str> = Vec::new();
    let mut flush = |header: &mut Option<String>, body: &mut Vec<&str>| {
        if let Some(header_line) = header.take() {
            let body_text = body.join("\n");
            let size_kb = fork_smaps_field(&body_text, "Size:");
            if size_kb >= min_size_kb {
                big.push(ForkSmapsMapping {
                    header: header_line.clone(),
                    size_kb,
                    rss_kb: fork_smaps_field(&body_text, "Rss:"),
                    pss_kb: fork_smaps_field(&body_text, "Pss:"),
                    private_clean_kb: fork_smaps_field(&body_text, "Private_Clean:"),
                    private_dirty_kb: fork_smaps_field(&body_text, "Private_Dirty:"),
                    shared_clean_kb: fork_smaps_field(&body_text, "Shared_Clean:"),
                    shared_dirty_kb: fork_smaps_field(&body_text, "Shared_Dirty:"),
                    name: header_line
                        .split_whitespace()
                        .nth(5)
                        .unwrap_or("[anon]")
                        .to_string(),
                });
            }
        }
        body.clear();
    };
    for line in smaps.lines() {
        if fork_is_smaps_header(line) {
            flush(&mut header, &mut body);
            header = Some(line.to_string());
        } else if header.is_some() {
            body.push(line);
        }
    }
    flush(&mut header, &mut body);
    big
}

/// THE COW SNAPSHOT: the rollup plus the big mappings, one
/// /proc/self read at a moment in a process's life.
struct ForkCowSnapshot {
    rollup: ForkSmapsRollup,
    big: Vec<ForkSmapsMapping>,
}

const FORK_COW_BIG_MAPPING_KB: u64 = 1 << 20; // >= 1 GiB

fn fork_cow_snapshot() -> io::Result<ForkCowSnapshot> {
    let rollup_text = std::fs::read_to_string("/proc/self/smaps_rollup")?;
    let rollup = ForkSmapsRollup {
        rss_kb: fork_smaps_field(&rollup_text, "Rss:"),
        pss_kb: fork_smaps_field(&rollup_text, "Pss:"),
        private_clean_kb: fork_smaps_field(&rollup_text, "Private_Clean:"),
        private_dirty_kb: fork_smaps_field(&rollup_text, "Private_Dirty:"),
        shared_clean_kb: fork_smaps_field(&rollup_text, "Shared_Clean:"),
        shared_dirty_kb: fork_smaps_field(&rollup_text, "Shared_Dirty:"),
    };
    let smaps = std::fs::read_to_string("/proc/self/smaps")?;
    Ok(ForkCowSnapshot {
        rollup,
        big: fork_smaps_big_mappings(&smaps, FORK_COW_BIG_MAPPING_KB),
    })
}

/// Print one snapshot as the [fork-cow] receipt lines.
fn fork_cow_print(label: &str, snap: &ForkCowSnapshot) {
    let r = &snap.rollup;
    eprintln!(
        "[fork-cow] {label} rollup: Rss={:.2}GB Pss={:.2}GB private={:.2}GB shared={:.2}GB",
        r.rss_kb as f64 / (1 << 20) as f64,
        r.pss_kb as f64 / (1 << 20) as f64,
        r.private_kb() as f64 / (1 << 20) as f64,
        r.shared_kb() as f64 / (1 << 20) as f64,
    );
    for m in &snap.big {
        eprintln!(
            "[fork-cow] {label} mapping {} size={:.2}GB rss={:.2}GB pss={:.2}GB private={:.2}GB shared={:.2}GB {}",
            m.header,
            m.size_kb as f64 / (1 << 20) as f64,
            m.rss_kb as f64 / (1 << 20) as f64,
            m.pss_kb as f64 / (1 << 20) as f64,
            m.private_kb() as f64 / (1 << 20) as f64,
            m.shared_kb() as f64 / (1 << 20) as f64,
            m.name,
        );
    }
}

/// The COW-hold receipt: the per-mapping and rollup PRIVATE growth of
/// a child's whole life (open snapshot -> close snapshot). The growth
/// in an INHERITED mapping's private pages IS the COW break for that
/// arena; the growth elsewhere is the child's own working set.
fn fork_cow_delta_print(label: &str, open: &ForkCowSnapshot, close: &ForkCowSnapshot) {
    // GB at two decimals for the big masses; exact MB below 1GB so a
    // page-class COW break is receipted as what it is, not rounded to
    // zero.
    let fmt_kb = |kb: u64| -> String {
        if kb >= (1 << 20) {
            format!("{:.2}GB", kb as f64 / (1 << 20) as f64)
        } else {
            format!("{}MB", kb / 1024)
        }
    };
    let fmt_delta = |kb_close: u64, kb_open: u64| -> String {
        if kb_close >= kb_open {
            fmt_kb(kb_close - kb_open)
        } else {
            format!("-{}", fmt_kb(kb_open - kb_close))
        }
    };
    let d_open = fmt_kb(open.rollup.private_kb());
    let d_close = fmt_kb(close.rollup.private_kb());
    eprintln!(
        "[fork-cow] {label} delta: private +{} (open {d_open} -> close {d_close}), pss {:+.2}GB, rss {:+.2}GB",
        fmt_delta(close.rollup.private_kb(), open.rollup.private_kb()),
        close.rollup.pss_kb as f64 / (1 << 20) as f64
            - open.rollup.pss_kb as f64 / (1 << 20) as f64,
        close.rollup.rss_kb as f64 / (1 << 20) as f64
            - open.rollup.rss_kb as f64 / (1 << 20) as f64,
    );
    for close_m in &close.big {
        if let Some(open_m) = open.big.iter().find(|m| m.header == close_m.header) {
            let grew = close_m.private_kb() > open_m.private_kb();
            eprintln!(
                "[fork-cow] {label} mapping {} private +{} (open {} -> close {}), shared {}{}",
                close_m.header,
                fmt_delta(close_m.private_kb(), open_m.private_kb()),
                fmt_kb(open_m.private_kb()),
                fmt_kb(close_m.private_kb()),
                fmt_delta(close_m.shared_kb(), open_m.shared_kb()),
                if grew { " [COW BREAK]" } else { "" },
            );
        } else {
            eprintln!(
                "[fork-cow] {label} mapping {} private +{} [CHILD-OWN]",
                close_m.header,
                fmt_kb(close_m.private_kb()),
            );
        }
    }
}

/// THE FORK-AFTER-LOAD DRIVER (SS10.6.9's designed form (a), the
/// slot-ceiling lever): the parent loads the panel ONCE and stays the
/// load owner; then ONE CHILD PER MANIFEST ROW forks and receipts its
/// territory (the 1-row batch cut — batch composition is the
/// receipted-invariant cut identity — plus the row's call over the
/// inherited panel), then exits. Byte-identity is by construction on
/// both sides: the cut's per-territory isolation is the batch gate's
/// own receipt (single == batch member), and the call is the combined
/// form's own construction (run_combined_call_rows) over the SAME
/// loaded panel.
///
/// PREFORK HYGIENE (the read-length-class globals lesson, applied to
/// fork): a child inherits the parent's whole address space at the
/// fork instant — every process-global the fresh-process trajectory
/// would RE-INITIALIZE must either be re-initialized in the child or
/// provably identical. The audit, every item:
/// - RAYON'S GLOBAL REGISTRY — THE FORK-DECISIVE ONE: a forked child
///   keeps only the calling thread, so a registry initialized in the
///   parent holds worker threads that DO NOT EXIST in the child
///   (rayon's global registry is set once per process and can never
///   be replaced); any par_iter that installs it would enqueue to
///   dead workers and hang — the first probe run receipted exactly
///   this hang (the call's seam threads' inner par_iters). THE FIX,
/// BY CONSTRUCTION: the parent NEVER touches rayon before the fork
/// (SyngIndex::load is rayon-free — the load path is the C reader's
/// own pthreads, mmap'd sidecars, and plain file reads — and the
/// parent builds NO AgcIndex, whose build would size its slots with
/// rayon::current_num_threads()); every child then initializes its
/// OWN global registry at its first rayon touch, exactly where the
/// committed fresh process does (the child's own AgcIndex build),
/// spawning the child's own 64 threads at the inherited
/// RAYON_NUM_THREADS. Every par_iter in the child — from the main
/// thread, the cut, the call's seam threads — joins the child's own
/// live registry: the committed fresh-process trajectory by
/// construction, no install() indirection.
/// - THE AGC INDEX: built per child (NOT inherited): the struct is
///   small (name maps + per-thread decompressor slots over the
///   page-cache-shared AGC file), the AGC's bulk is the page cache
///   (shared system-wide under every form), and the per-child build
///   is the same initialization-order point as the committed
///   process. The per-child build wall is receipted per child.
/// - THE READ LENGTH: per-call state (a local in run_injected,
///   initialized at the fresh-process 150 at every call entry) — the
///   forked child reproduces the committed fresh-process trajectory
///   BY CONSTRUCTION, no reset discipline.
/// - THE C DIAGNOSTIC COUNTERS (pathCount): set once at load
///   (impg_syng_suppress_debug) — the child inherits the
///   already-suppressed state, exactly the fresh process's own
///   post-load state.
/// - THE LOCATE CHECKPOINT OnceLock: never built (the panel's .spos
///   is present — SampledCheckpointIndex serves locate); if it were,
///   a fresh child would build its own from the same .pstep bytes.
/// - THREAD-LOCALS (aligners, AGC decompressor pools): per-OS-thread;
///   a child's threads are all fresh, so every thread-local
///   initializes in the child exactly as in a fresh process.
/// - MALLOC ARENAS: the child's small allocations are served from the
///   inherited arena chunks (free-list reuse dirties those pages);
///   the big C arenas are never written by the query paths (the COW
///   audit's code half — see SS10.6.12) and the smaps receipts measure
///   the total honestly.
fn run_fork_cut_call_batch_global(
    global_panel_prefix: &str,
    agc_path: &Path,
    shards_dir: &Path,
    manifest_path: &Path,
    routes_dir: &Path,
    call_rss_budget_gib: f64,
    rows_per_child: usize,
) -> io::Result<()> {
    let specs = cut_batch_global_parse_manifest(manifest_path)?;
    ensure(
        specs.iter().any(|spec| spec.call.is_some()),
        "no manifest row carries a call object (cut-batch-global is the cut-only form)",
    )?;
    ensure(rows_per_child >= 1, "--rows-per-child must be >= 1")?;
    eprintln!(
        "[fork-cut-call-batch-global] loading the global panel once (the parent is the load owner) for {} manifest rows in {} children (rows-per-child {})",
        specs.len(),
        specs.len().div_ceil(rows_per_child),
        rows_per_child
    );
    let load_started = Instant::now();
    let global = SyngIndex::load(global_panel_prefix, SyncmerParams::default())?;
    let load_seconds = load_started.elapsed().as_secs_f64();
    let parent_load_snap = fork_cow_snapshot()?;
    eprintln!(
        "[fork-cut-call-batch-global] panel loaded in {:.0}s; the children inherit it COW (each child builds its own AGC index, the fresh-process form)",
        load_seconds
    );
    fork_cow_print("parent-post-load", &parent_load_snap);

    // Prefork hygiene: no parent-buffered stdio byte re-emerges in a
    // child (the child _exits without flushing).
    let _ = std::io::stdout().flush();
    let _ = std::io::stderr().flush();

    let fork_started = Instant::now();
    let shared_panel: &SyngIndex = &global;
    let mut child_pids: Vec<(usize, String, libc::pid_t)> = Vec::new();
    for (chunk, chunk_specs) in specs.chunks(rows_per_child).enumerate() {
        let row_lo = chunk * rows_per_child;
        let row_hi = row_lo + chunk_specs.len() - 1;
        let pid = unsafe { libc::fork() };
        if pid < 0 {
            return Err(io::Error::last_os_error());
        }
        if pid == 0 {
            // ================= THE CHILD (never returns) ================
            let exit_code: i32 = 'child: {
                let name = chunk_specs[0].name.clone();
                // Per-child lane-trace file (receipt-side; default OFF):
                // when the parent's environment names a trace path,
                // each child appends ".{chunk}" so a multi-child wave's
                // traces stay attributable per child (each child is its
                // own process; the env rewrite is the child's own).
                if let Ok(trace) = std::env::var("IMPG_LANE_TRACE") {
                    if !trace.is_empty() {
                        std::env::set_var("IMPG_LANE_TRACE", format!("{trace}.{chunk}"));
                    }
                }
                eprintln!(
                    "[fork-cut-call-batch-global] child chunk={chunk} rows={row_lo}-{row_hi} ({name}) pid={} open",
                    std::process::id()
                );
                let open_snap = match fork_cow_snapshot() {
                    Ok(snap) => snap,
                    Err(e) => {
                        eprintln!("[fork-cow] child chunk={chunk} open snapshot failed: {e}");
                        break 'child 1;
                    }
                };
                fork_cow_print(
                    &format!("child-open chunk={chunk} rows={row_lo}-{row_hi} ({name})"),
                    &open_snap,
                );
                let child_started = Instant::now();
                let work = || -> io::Result<()> {
                    // The child's own AGC index: the AgcIndex struct is
                    // small (name maps + per-thread decompressor
                    // slots over the page-cache-shared AGC file), so
                    // inheriting it would buy nothing — and building
                    // it here means the child's rayon GLOBAL registry
                    // initializes at the same point in the sequence as
                    // the committed fresh process (the first rayon
                    // touch), spawning the child's OWN 64 threads at
                    // the inherited RAYON_NUM_THREADS.
                    let agc_started = Instant::now();
                    let agc = AgcIndex::build_from_files(&[agc_path.display().to_string()])?;
                    let agc_seconds = agc_started.elapsed().as_secs_f64();
                    eprintln!(
                        "[fork-cut-call-batch-global] child chunk={chunk} built its own AGC index in {:.0}s",
                        agc_seconds
                    );
                    // THE CUT: the child's own batch over its chunk of
                    // manifest rows (the batch's per-territory
                    // isolation gate receipt: the batch composition is
                    // byte-identity-invariant, single == batch member;
                    // rows-per-child = 1 is the committed single-cut
                    // form exactly, and the width-4 chunk is the
                    // committed batch shape).
                    let cut_started = Instant::now();
                    cut_batch_global_with(
                        shared_panel,
                        &agc,
                        shards_dir,
                        chunk_specs.to_vec(),
                        load_seconds,
                    )?;
                    eprintln!(
                        "[fork-cut-call-batch-global] child chunk={chunk} cut phase done in {:.0}s",
                        cut_started.elapsed().as_secs()
                    );
                    // MEMORY HYGIENE RECEIPT (the N=1 probe's finding:
                    // the cut scan's transient working set — the shard
                    // decode/remap arenas, ~35-40GB, width-insensitive
                    // — is retained by the allocator after the cut
                    // phase; this is freed memory, never live data).
                    // malloc_trim(0) asks glibc to return every
                    // arena's free tail to the OS: free-list state is
                    // not program semantics, and the receipt prints
                    // the before/after private so the effect is
                    // measured, never assumed.
                    let pre_trim = fork_cow_snapshot()?;
                    unsafe { libc::malloc_trim(0) };
                    let post_trim = fork_cow_snapshot()?;
                    let gb = |kb: u64| kb as f64 / (1 << 20) as f64;
                    eprintln!(
                        "[fork-cow] child chunk={chunk} malloc_trim: private {:.2}GB -> {:.2}GB (returned {:.2}GB)",
                        gb(pre_trim.rollup.private_kb()),
                        gb(post_trim.rollup.private_kb()),
                        gb(pre_trim
                            .rollup
                            .private_kb()
                            .saturating_sub(post_trim.rollup.private_kb())),
                    );
                    // THE CALL: the combined form's own construction
                    // over the inherited panel (the chunk's call rows,
                    // concurrent — the manifest's own width).
                    let call_rows: Vec<(String, String, u64, PathBuf, CutBatchGlobalCall)> =
                        chunk_specs
                            .iter()
                            .filter_map(|spec| {
                                spec.call.as_ref().map(|call| {
                                    (
                                        spec.name.clone(),
                                        spec.component.clone(),
                                        spec.record_origin,
                                        spec.out_dir.clone(),
                                        call.clone(),
                                    )
                                })
                            })
                            .collect();
                    if !call_rows.is_empty() {
                        let results = run_combined_call_rows(
                            shared_panel,
                            global_panel_prefix,
                            routes_dir,
                            &call_rows,
                            call_rss_budget_gib,
                            load_seconds,
                            "fork-cut-call-batch-global",
                        );
                        for result in results {
                            result?;
                        }
                    }
                    Ok(())
                };
                match work() {
                    Ok(()) => {
                        match fork_cow_snapshot() {
                            Ok(close_snap) => {
                                fork_cow_delta_print(
                                    &format!("child chunk={chunk} rows={row_lo}-{row_hi} ({name}) open->close"),
                                    &open_snap,
                                    &close_snap,
                                );
                                fork_cow_print(
                                    &format!("child-close chunk={chunk} rows={row_lo}-{row_hi} ({name})"),
                                    &close_snap,
                                );
                            }
                            Err(e) => {
                                eprintln!(
                                    "[fork-cow] child chunk={chunk} close snapshot failed: {e}"
                                )
                            }
                        }
                        eprintln!(
                            "[fork-cut-call-batch-global] child chunk={chunk} rows={row_lo}-{row_hi} ({name}) done in {:.0}s",
                            child_started.elapsed().as_secs()
                        );
                        0
                    }
                    Err(e) => {
                        eprintln!(
                            "[fork-cut-call-batch-global] child chunk={chunk} rows={row_lo}-{row_hi} ({name}) FAILED: {e}"
                        );
                        1
                    }
                }
            };
            unsafe { libc::_exit(exit_code) };
        }
        let child_name = chunk_specs[0].name.clone();
        child_pids.push((chunk, child_name.clone(), pid));
        eprintln!(
            "[fork-cut-call-batch-global] forked child chunk={chunk} rows={row_lo}-{row_hi} ({child_name}) pid={pid}"
        );
    }
    // THE PARENT: the load owner and the waiter. Children receipt
    // their own territories and exit; the parent fails at the FIRST
    // failing child in CHUNK ORDER (the serial form's first-failure
    // discipline).
    let mut first_failure: Option<(usize, String, i32)> = None;
    for (chunk, name, pid) in &child_pids {
        let mut status: libc::c_int = 0;
        if unsafe { libc::waitpid(*pid, &mut status, 0) } < 0 {
            return Err(io::Error::last_os_error());
        }
        let code = if libc::WIFEXITED(status) {
            libc::WEXITSTATUS(status)
        } else {
            -1
        };
        if code != 0 && first_failure.is_none() {
            first_failure = Some((*chunk, name.clone(), code));
        }
    }
    eprintln!(
        "[fork-cut-call-batch-global] all {} children exited in {:.0}s (load {:.0}s + children)",
        child_pids.len(),
        fork_started.elapsed().as_secs(),
        load_seconds
    );
    if let Some(fork_cow_close) = fork_cow_snapshot().ok() {
        fork_cow_print("parent-at-close", &fork_cow_close);
    }
    if let Some((chunk, name, code)) = first_failure {
        return Err(io::Error::other(format!(
            "child chunk={chunk} ({name}) failed with exit code {code}"
        )));
    }
    Ok(())
}

/// One global-frame batch territory's emission: the single global
/// cut's emission tail (derive cache, reads FASTA, census routing over
/// the territory index, stats), over the territory's merged per-file
/// state.
fn batch_global_emit_territory(
    panel: &SyngIndex,
    t: GlobalBatchTerritory,
    files: Vec<FileCut>,
    load_seconds: f64,
    scan_seconds: f64,
    started: &Instant,
) -> io::Result<serde_json::Value> {
    let GlobalBatchTerritory {
        spec,
        territory,
        territory_nodes_abs: _,
        node_map: _,
        node_map_pairs,
        component_contig_source,
    } = t;
    let out_dir = spec.out_dir.as_path();
    let record_origin = spec.record_origin;
    let merged = batch_merge_file(MergedBatch::new(), files);
    // The derive cache (the committed format): multiplicity by
    // identical sequence, first representative carries the count.
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

    // The cut census: route every touching pattern through the
    // committed placement machinery; emit the committed JSONL schema
    // (the single global cut's emission, verbatim — the component
    // contig comparison in the panel's own frame).
    let census_started = Instant::now();
    let k = panel.syncmer_length_bp() as u64;
    let mut out = io::BufWriter::new(std::fs::File::create(out_dir.join("cut-census.jsonl"))?);
    let component_contig = component_contig_source
        .rsplit('#')
        .next()
        .unwrap_or(&component_contig_source)
        .to_string();
    let mut routed_patterns = 0u64;
    let mut census_mult_total = 0u64;
    let mut dropped_no_occurrence = 0u64;
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
        frame: "global".into(),
        // The batch forms' own scan is the full routed mass (the
        // read-set instrument does not apply on the batch path; the
        // single cut's --census-profile is the three-step interface's
        // own switch — the named follow-up carries the profile into
        // the batch manifests the same way record_origin rides them).
        census_profile: Some("global-native".into()),
        record_origin,
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
        offset_variants: BTreeMap::new(),
        partial_remap: false,
        local_derive: false,
        global_frame: true,
        reads_routed: merged.reads_routed,
        reads_unrouted: merged.reads_unrouted,
        load_seconds,
        scan_seconds,
        census_seconds,
        wall_seconds,
        peak_rss_kb: rss_now_kb(),
    };
    let text = serde_json::to_string_pretty(&stats)?;
    std::fs::write(out_dir.join("cut-stats.json"), text + "\n")?;
    Ok(serde_json::json!({
        "name": spec.name,
        "reads_kept": stats.reads_kept,
        "routed_patterns": stats.routed_patterns,
        "census_mult": stats.census_mult,
    }))
}

// ---------------------------------------------------------------------------
// trace: THE YIELD-ANOMALY AUTOPSY DIAGNOSTIC (receipt-side only; the
// product/scoreboard machinery untouched). One scan of the projection
// shards; per locality and read-class (sequence sets given as FASTA),
// every record's remap fate is classified: touching, mappable without
// touch, or unmappable - with the unmapped anchors split into the
// in-slice class (the anchor position is a global-walk position of the
// locality's own slices: the dictionary-superset class) and the
// out-of-slice class (the record's walk anchors outside the locality's
// sequences entirely: the paralog/walk-exit class) - plus the
// partial-remap recovery projection (would the mapped subwalk touch?).
// Query reads get a per-record detail dump.
// ---------------------------------------------------------------------------

#[derive(serde::Deserialize)]
struct TraceLocalitySpec {
    name: String,
    panel: String,
    axis: PathBuf,
    bed_directory: PathBuf,
    component: String,
    sources: PathBuf,
}

#[derive(serde::Deserialize)]
struct TraceClassSpec {
    /// "<locality-name>:<label>".
    class: String,
    fasta: PathBuf,
}

struct TraceLocality {
    name: String,
    node_map: HashMap<i32, i32>,
    territory_nodes_abs: HashSet<u32>,
    /// Every global node id appearing in any locality slice's global
    /// canonical walk (the in-slice reference for the unmapped-anchor
    /// split).
    slice_global_nodes: HashSet<i32>,
}

#[derive(Default, Clone, Serialize)]
struct Agg {
    reads: u64,
    reads_no_records: u64,
    reads_touch: u64,
    touching_records: u64,
    notouch_records: u64,
    unmappable_records: u64,
    /// reads not kept whose failing records' unmapped anchors are all
    /// in-slice (the dictionary-superset loss class).
    reads_unm_inslice_only: u64,
    /// reads not kept with at least one out-of-slice unmapped anchor
    /// (the paralog/walk-exit loss class).
    reads_unm_any_outslice: u64,
    /// reads not kept, no unmappable record, but no touch either.
    reads_notouch_only: u64,
    unmapped_anchors: u64,
    unmapped_in_slice: u64,
    unmapped_out_slice: u64,
    /// unmappable records whose mapped subwalk would touch (the
    /// partial-remap recovery projection).
    recover_touch_records: u64,
    /// reads not kept that partial remap would recover.
    recovered_reads: u64,
    /// anchor-count histogram of touching records (16 bins, cap 15+).
    touch_anchor_hist: Vec<u64>,
    /// anchor-count histogram of recovered (partial-remap touch) records.
    recover_anchor_hist: Vec<u64>,
}

impl Agg {
    fn hist_bump(hist: &mut Vec<u64>, count: usize) {
        if hist.is_empty() {
            *hist = vec![0; 16];
        }
        hist[count.min(15)] += 1;
    }
}

struct FileTrace {
    reads_scanned: u64,
    agg: HashMap<(u32, u32), Agg>,
    detail: Vec<serde_json::Value>,
}

#[derive(Clone)]
enum RecStatus {
    Touch {
        anchors: usize,
    },
    NoTouch {
        anchors: usize,
    },
    Unmappable {
        mapped: Vec<(i32, u64)>,
        unmapped: Vec<(i32, bool)>,
    },
}

fn trace_classify_record(
    tokens: &[u64],
    walk: &OwnWalk,
    locality: &TraceLocality,
) -> io::Result<RecStatus> {
    let anchors = decode_tokens(tokens)?;
    let mut mapped = Vec::with_capacity(anchors.len());
    let mut unmapped = Vec::new();
    for &(node, pos) in &anchors {
        match locality.node_map.get(&node) {
            Some(&local) => mapped.push((local, pos)),
            None => unmapped.push((node, locality.slice_global_nodes.contains(&node))),
        }
    }
    let mut walk_mapped = 0usize;
    for &(node, _) in walk {
        if locality.node_map.contains_key(&node) {
            walk_mapped += 1;
        }
    }
    if unmapped.is_empty() && walk_mapped == walk.len() {
        let touch = mapped
            .iter()
            .any(|&(node, _)| locality.territory_nodes_abs.contains(&node.unsigned_abs()));
        if touch {
            Ok(RecStatus::Touch { anchors: mapped.len() })
        } else {
            Ok(RecStatus::NoTouch { anchors: mapped.len() })
        }
    } else {
        Ok(RecStatus::Unmappable { mapped, unmapped })
    }
}

fn scan_trace_file(
    path: &Path,
    localities: &[TraceLocality],
    memberships: &HashMap<Vec<u8>, Vec<(u32, u32)>>,
    query_set: &HashSet<Vec<u8>>,
) -> io::Result<FileTrace> {
    let mut trace = FileTrace {
        reads_scanned: 0,
        agg: HashMap::new(),
        detail: Vec::new(),
    };
    for batch in shard_packet_stream(path, 1024)? {
        let packets = batch?;
        trace.reads_scanned += packets.len() as u64;
        let per_read: Vec<io::Result<(
            Vec<u8>,
            Vec<Vec<RecStatus>>,
            Vec<(u32, u32)>,
            bool,
        )>> = packets
            .par_iter()
            .map(|packet| {
                let read = decode_read_packet(packet)?;
                let memberships = memberships.get(&read.seq);
                let is_query = query_set.contains(&read.seq);
                if memberships.is_none() && !is_query {
                    return Ok((read.seq, Vec::new(), Vec::new(), false));
                }
                let per_loc: Vec<Vec<RecStatus>> = localities
                    .iter()
                    .map(|locality| {
                        let mut out = Vec::with_capacity(read.records.len());
                        for (tokens, walk) in &read.records {
                            out.push(trace_classify_record(tokens, walk, locality)?);
                        }
                        Ok(out)
                    })
                    .collect::<io::Result<_>>()?;
                Ok((
                    read.seq,
                    per_loc,
                    memberships.cloned().unwrap_or_default(),
                    is_query,
                ))
            })
            .collect();
        for (seq, per_loc, memberships, is_query) in
            per_read.into_iter().collect::<io::Result<Vec<_>>>()? {
            if memberships.is_empty() && !is_query {
                continue;
            }
            let qhash = if is_query { Some(fnv1a64(&seq)) } else { None };
            let mut detail_loc: Vec<serde_json::Value> = Vec::new();
            for (loc_index, statuses) in per_loc.into_iter().enumerate() {
                let mut detail_recs: Option<Vec<serde_json::Value>> =
                    is_query.then(|| Some(Vec::new())).flatten();
                let mut any_touch = false;
                let mut any_unm = false;
                let mut any_notouch = false;
                let mut any_outslice = false;
                let mut recover_count = 0u64;
                let mut agg: Option<&mut Agg> = None;
                for &(loc_member, class_index) in &memberships {
                    if loc_member as usize == loc_index {
                        agg = Some(trace.agg.entry((loc_index as u32, class_index)).or_default());
                        let a = agg.as_mut().unwrap();
                        a.reads += 1;
                        if statuses.is_empty() {
                            a.reads_no_records += 1;
                        }
                        break;
                    }
                }
                for status in statuses {
                    match status {
                        RecStatus::Touch { anchors } => {
                            any_touch = true;
                            if let Some(a) = agg.as_deref_mut() {
                                a.touching_records += 1;
                                Agg::hist_bump(&mut a.touch_anchor_hist, anchors);
                            }
                            if let Some(recs) = detail_recs.as_mut() {
                                recs.push(serde_json::json!({"s": "touch", "a": anchors}));
                            }
                        }
                        RecStatus::NoTouch { anchors } => {
                            any_notouch = true;
                            if let Some(a) = agg.as_deref_mut() {
                                a.notouch_records += 1;
                            }
                            if let Some(recs) = detail_recs.as_mut() {
                                recs.push(serde_json::json!({"s": "notouch", "a": anchors}));
                            }
                        }
                        RecStatus::Unmappable { mapped, unmapped } => {
                            any_unm = true;
                            let recover = !mapped.is_empty()
                                && mapped.iter().any(|&(node, _)| {
                                    localities[loc_index]
                                        .territory_nodes_abs
                                        .contains(&node.unsigned_abs())
                                });
                            if recover {
                                recover_count += 1;
                            }
                            for &(node, in_slice) in &unmapped {
                                if !in_slice {
                                    any_outslice = true;
                                }
                                if let Some(a) = agg.as_deref_mut() {
                                    a.unmapped_anchors += 1;
                                    if in_slice {
                                        a.unmapped_in_slice += 1;
                                    } else {
                                        a.unmapped_out_slice += 1;
                                    }
                                }
                            }
                            if let Some(a) = agg.as_deref_mut() {
                                a.unmappable_records += 1;
                                if recover {
                                    a.recover_touch_records += 1;
                                    Agg::hist_bump(&mut a.recover_anchor_hist, mapped.len());
                                }
                            }
                            if let Some(recs) = detail_recs.as_mut() {
                                recs.push(serde_json::json!({
                                    "s": "unm",
                                    "r": recover,
                                    "m": mapped.iter().map(|&(n, p)| vec![n as i64, p as i64]).collect::<Vec<_>>(),
                                    "u": unmapped.iter().map(|&(n, s)| vec![n, if s {1} else {0}]).collect::<Vec<_>>(),
                                }));
                            }
                        }
                    }
                }
                if let Some(a) = agg.as_deref_mut() {
                    if any_touch {
                        a.reads_touch += 1;
                    } else {
                        if any_unm {
                            if any_outslice {
                                a.reads_unm_any_outslice += 1;
                            } else {
                                a.reads_unm_inslice_only += 1;
                            }
                        } else if any_notouch {
                            a.reads_notouch_only += 1;
                        }
                        if recover_count > 0 {
                            a.recovered_reads += 1;
                        }
                    }
                }
                if is_query {
                    if let Some(recs) = detail_recs {
                        detail_loc.push(serde_json::json!({
                            "loc": localities[loc_index].name,
                            "recs": recs,
                        }));
                    }
                }
            }
            if let Some(qhash) = qhash {
                trace.detail.push(serde_json::json!({
                    "q": qhash,
                    "locs": detail_loc,
                }));
            }
        }
    }
    Ok(trace)
}

fn load_fasta_sequences(path: &Path) -> io::Result<Vec<Vec<u8>>> {
    let mut seqs = Vec::new();
    let mut cur: Vec<u8> = Vec::new();
    for line in BufReader::new(std::fs::File::open(path)?).lines() {
        let line = line?;
        if let Some(_header) = line.strip_prefix('>') {
            if !cur.is_empty() {
                seqs.push(std::mem::take(&mut cur));
            }
        } else if !line.is_empty() {
            cur.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !cur.is_empty() {
        seqs.push(cur);
    }
    Ok(seqs)
}

/// The per-path node correspondence + territory for one trace locality
/// (the committed cut's own construction, with the slice-global-node set
/// additionally collected for the in-slice split).
fn trace_load_locality(
    spec: &TraceLocalitySpec,
    global: &SyngIndex,
    agc: &AgcIndex,
    k: u64,
) -> io::Result<TraceLocality> {
    let panel = SyngIndex::load(&spec.panel, SyncmerParams::default())?;
    ensure(
        global.syncmer_length_bp() as u64 == k,
        "global and locality syncmer lengths differ",
    )?;
    let axis: AxisFile = serde_json::from_str(&std::fs::read_to_string(&spec.axis)?)
        .map_err(|e| invalid(format!("axis parse: {e}")))?;
    let universe = load_universe_component(&axis, &spec.bed_directory, &spec.component)?;
    let mut sources: HashMap<String, Vec<u8>> = HashMap::new();
    {
        let text = std::fs::read_to_string(&spec.sources)?;
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
    let fetch_lane = lane_fetch_over(&sources);
    let territory = build_territory_index(&panel, &universe, &fetch_lane, k)?;
    let territory_nodes_abs: HashSet<u32> = territory
        .entries
        .iter()
        .map(|&(node, _, _)| node.unsigned_abs())
        .collect();
    let mut node_map: HashMap<i32, i32> = HashMap::new();
    let mut slice_global_nodes: HashSet<i32> = HashSet::new();
    for (header, record) in &sources {
        let (global_name, s, e) = parse_locus_header(header)?;
        ensure(
            (e - s) as usize == record.len(),
            &format!("locus record {header} length differs from its header interval"),
        )?;
        let mut delta: Option<i64> = None;
        for candidate in [0i64, -1, 1] {
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
        let lo = (s as i64 + delta).max(0) as i32;
        let global_seq = agc
            .fetch_sequence(&global_name, lo, lo + record.len() as i32)
            .map_err(|e| invalid(format!("AGC fetch {global_name}: {e}")))?;
        ensure(global_seq == *record, "AGC refetch mismatch")?;
        let local_walk = mem_records::canonical_anchor_walk(&panel, record)?;
        let global_walk = mem_records::canonical_anchor_walk(global, &global_seq)?;
        let global_by_pos: HashMap<u64, i32> =
            global_walk.iter().map(|&(node, pos)| (pos, node)).collect();
        for &(gnode, _pos) in &global_walk {
            slice_global_nodes.insert(gnode);
        }
        for &(lnode, lpos) in &local_walk {
            let gnode = *global_by_pos
                .get(&lpos)
                .ok_or_else(|| invalid("locality anchor position absent from the global walk"))?;
            match node_map.entry(gnode) {
                std::collections::hash_map::Entry::Occupied(o) => {
                    ensure(*o.get() == lnode, "conflicting node correspondence")?;
                }
                std::collections::hash_map::Entry::Vacant(v) => {
                    v.insert(lnode);
                }
            }
        }
    }
    eprintln!(
        "[trace] locality {}: {} node-map pairs, {} slice global nodes, {} territory rows",
        spec.name,
        node_map.len(),
        slice_global_nodes.len(),
        territory.territories.len()
    );
    Ok(TraceLocality {
        name: spec.name.clone(),
        node_map,
        territory_nodes_abs,
        slice_global_nodes,
    })
}

#[allow(clippy::too_many_arguments)]
fn run_trace(
    global_panel_prefix: &str,
    agc_path: &Path,
    shards_dir: &Path,
    localities_path: &Path,
    classes_path: &Path,
    queries_path: Option<&Path>,
    out_dir: &Path,
) -> io::Result<()> {
    std::fs::create_dir_all(out_dir)?;
    let started = Instant::now();
    eprintln!("[trace] loading global panel");
    let global = SyngIndex::load(global_panel_prefix, SyncmerParams::default())?;
    let agc = AgcIndex::build_from_files(&[agc_path.display().to_string()])?;
    let k = global.syncmer_length_bp() as u64;
    let locality_specs: Vec<TraceLocalitySpec> =
        serde_json::from_str(&std::fs::read_to_string(localities_path)?)
            .map_err(|e| invalid(format!("localities spec parse: {e}")))?;
    let class_specs: Vec<TraceClassSpec> =
        serde_json::from_str(&std::fs::read_to_string(classes_path)?)
            .map_err(|e| invalid(format!("classes spec parse: {e}")))?;
    let mut localities: Vec<TraceLocality> = Vec::with_capacity(locality_specs.len());
    for spec in &locality_specs {
        localities.push(trace_load_locality(spec, &global, &agc, k)?);
    }
    let locality_index: HashMap<&str, u32> = locality_specs
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i as u32))
        .collect();
    let mut class_labels: Vec<String> = Vec::new();
    let mut label_index: HashMap<String, u32> = HashMap::new();
    let mut memberships: HashMap<Vec<u8>, Vec<(u32, u32)>> = HashMap::new();
    let mut class_summary: Vec<serde_json::Value> = Vec::new();
    for spec in &class_specs {
        let (loc_name, label) = spec
            .class
            .split_once(':')
            .ok_or_else(|| invalid(&format!("class {} is not <locality>:<label>", spec.class)))?;
        let loc_index = *locality_index
            .get(loc_name)
            .ok_or_else(|| invalid(&format!("class {} names unknown locality", spec.class)))?;
        let class_index = match label_index.get(label) {
            Some(&index) => index,
            None => {
                let index = class_labels.len() as u32;
                class_labels.push(label.to_string());
                label_index.insert(label.to_string(), index);
                index
            }
        };
        let seqs = load_fasta_sequences(&spec.fasta)?;
        let mut members = 0u64;
        for seq in seqs {
            memberships
                .entry(seq)
                .or_default()
                .push((loc_index, class_index));
            members += 1;
        }
        class_summary.push(serde_json::json!({
            "class": spec.class,
            "sequences": members,
        }));
        eprintln!("[trace] class {}: {} sequences", spec.class, members);
    }
    let mut query_set: HashSet<Vec<u8>> = HashSet::new();
    if let Some(path) = queries_path {
        for seq in load_fasta_sequences(path)? {
            query_set.insert(seq);
        }
        eprintln!("[trace] {} query sequences", query_set.len());
    }
    let load_seconds = started.elapsed().as_secs_f64();

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
    let per_file: Vec<io::Result<FileTrace>> = shard_paths
        .par_iter()
        .map(|path| scan_trace_file(path, &localities, &memberships, &query_set))
        .collect();
    let mut merged = FileTrace {
        reads_scanned: 0,
        agg: HashMap::new(),
        detail: Vec::new(),
    };
    for result in per_file {
        let trace = result?;
        merged.reads_scanned += trace.reads_scanned;
        merged.detail.extend(trace.detail);
        for (key, agg) in trace.agg {
            let entry = merged.agg.entry(key).or_default();
            merge_agg(entry, agg);
        }
    }
    merged.detail.sort_by_key(|value| {
        value["q"].as_u64().unwrap_or(0)
    });
    {
        let mut out = io::BufWriter::new(std::fs::File::create(out_dir.join("trace-detail.jsonl"))?);
        for value in &merged.detail {
            serde_json::to_writer(&mut out, value)?;
            out.write_all(b"\n")?;
        }
        out.flush()?;
    }
    #[derive(Serialize)]
    struct TraceAggJson {
        reads_scanned: u64,
        classes: Vec<String>,
        class_summary: Vec<serde_json::Value>,
        localities: Vec<serde_json::Value>,
        load_seconds: f64,
        wall_seconds: f64,
        peak_rss_kb: u64,
    }
    let locality_json: Vec<serde_json::Value> = localities
        .iter()
        .enumerate()
        .map(|(loc_index, locality)| {
            let classes: Vec<serde_json::Value> = class_labels
                .iter()
                .enumerate()
                .map(|(class_index, label)| {
                    let agg = merged
                        .agg
                        .get(&(loc_index as u32, class_index as u32));
                    let value = match agg {
                        Some(a) => serde_json::to_value(a).unwrap_or(serde_json::Value::Null),
                        None => serde_json::Value::Null,
                    };
                    serde_json::json!({ "label": label, "agg": value })
                })
                .collect();
            serde_json::json!({ "name": locality.name, "classes": classes })
        })
        .collect();
    let stats = TraceAggJson {
        reads_scanned: merged.reads_scanned,
        classes: class_labels.clone(),
        class_summary,
        localities: locality_json,
        load_seconds,
        wall_seconds: started.elapsed().as_secs_f64(),
        peak_rss_kb: rss_now_kb(),
    };
    let text = serde_json::to_string_pretty(&stats)?;
    std::fs::write(out_dir.join("trace-agg.json"), text + "\n")?;
    eprintln!(
        "[trace] DONE: {} reads scanned, {} detail lines [{:.1}s]",
        merged.reads_scanned,
        merged.detail.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn merge_agg(dst: &mut Agg, src: Agg) {
    dst.reads += src.reads;
    dst.reads_no_records += src.reads_no_records;
    dst.reads_touch += src.reads_touch;
    dst.touching_records += src.touching_records;
    dst.notouch_records += src.notouch_records;
    dst.unmappable_records += src.unmappable_records;
    dst.reads_unm_inslice_only += src.reads_unm_inslice_only;
    dst.reads_unm_any_outslice += src.reads_unm_any_outslice;
    dst.reads_notouch_only += src.reads_notouch_only;
    dst.unmapped_anchors += src.unmapped_anchors;
    dst.unmapped_in_slice += src.unmapped_in_slice;
    dst.unmapped_out_slice += src.unmapped_out_slice;
    dst.recover_touch_records += src.recover_touch_records;
    dst.recovered_reads += src.recovered_reads;
    merge_hist(&mut dst.touch_anchor_hist, src.touch_anchor_hist);
    merge_hist(&mut dst.recover_anchor_hist, src.recover_anchor_hist);
}

fn merge_hist(dst: &mut Vec<u64>, src: Vec<u64>) {
    if dst.is_empty() {
        *dst = src;
    } else if !src.is_empty() {
        if dst.len() < src.len() {
            dst.resize(src.len(), 0);
        }
        for (i, v) in src.into_iter().enumerate() {
            dst[i] += v;
        }
    }
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
            region_batching,
        } => run_project(&panel, &inputs, &out_dir, map_cap, progress, limit, region_batching),
        Mode::Dedup {
            out_dir,
            occurrences,
            partitions,
        } => run_dedup(&out_dir, &occurrences, partitions),
        Mode::Cut {
            global_panel,
            frame,
            record_origin,
            census_profile,
            panel,
            axis,
            bed_directory,
            component,
            sources,
            agc,
            shards_dir,
            out_dir,
        } => {
            let panel = frame.select(panel, "--panel")?;
            let sources = frame.select(sources, "--sources")?;
            ensure(
                record_origin == 0 || frame == CutFrame::Global,
                "--record-origin belongs to the global frame (the committed local frame \
                 derives the content convention per record against the AGC, fail-closed)",
            )?;
            ensure(
                frame == CutFrame::Global || census_profile == CensusProfile::Restricted,
                "--census-profile belongs to the global frame (the committed local path \
                 is byte-for-byte; the read-set instrument still restricts there)",
            )?;
            // THE CENSUS PROFILE (the owner's Ruling 2, executed):
            // resolve BEFORE the panel load — the fail-closed form
            // must not pay the ~325s global load to learn it.
            let readset_basis = std::env::var("IMPG_CUT_READSET_FASTA")
                .ok()
                .filter(|value| !value.is_empty());
            let (census_profile, basis_override) =
                CensusProfile::resolve(frame, census_profile, readset_basis.as_deref())?;
            if frame == CutFrame::Global {
                eprintln!(
                    "[cut] census profile: {}{}",
                    census_profile.label(),
                    if basis_override { " (IMPG_CUT_READSET_FASTA override)" } else { "" }
                );
            }
            run_cut(
                &global_panel,
                frame,
                record_origin,
                census_profile,
                panel.as_deref(),
                &axis,
                &bed_directory,
                &component,
                sources.as_deref(),
                &agc,
                &shards_dir,
                &out_dir,
            )
        }
        Mode::CutBatch {
            global_panel,
            agc,
            shards_dir,
            manifest,
        } => run_cut_batch(&global_panel, &agc, &shards_dir, &manifest),
        Mode::CutBatchGlobal {
            global_panel,
            agc,
            shards_dir,
            manifest,
        } => run_cut_batch_global(&global_panel, &agc, &shards_dir, &manifest),
        Mode::CutCallBatchGlobal {
            global_panel,
            agc,
            shards_dir,
            manifest,
            routes,
            call_rss_budget_gib,
        } => run_cut_call_batch_global(
            &global_panel,
            &agc,
            &shards_dir,
            &manifest,
            &routes,
            call_rss_budget_gib,
        ),
        Mode::ForkCutCallBatchGlobal {
            global_panel,
            agc,
            shards_dir,
            manifest,
            routes,
            call_rss_budget_gib,
            rows_per_child,
        } => run_fork_cut_call_batch_global(
            &global_panel,
            &agc,
            &shards_dir,
            &manifest,
            &routes,
            call_rss_budget_gib,
            rows_per_child,
        ),
        Mode::Trace {
            global_panel,
            agc,
            shards_dir,
            localities,
            classes,
            queries,
            out_dir,
        } => run_trace(
            &global_panel,
            &agc,
            &shards_dir,
            &localities,
            &classes,
            queries.as_deref(),
            &out_dir,
        ),
    }
}

// ---------------------------------------------------------------------------
// THE GLOBAL-FRAME CONVERSION PROOF (the mini-panel test): the
// territory index against the ONE global panel over the TRANSLATED
// rows routes a read to exactly the occurrences the committed local
// substrate produces, modulo the coordinate rename — the
// conversion-not-materialization property at test scale (the re-anchor
// principle: node ids are opaque interned k-mer tokens; every semantic
// quantity is coordinate- and content-carried).
// ---------------------------------------------------------------------------
#[cfg(test)]
mod global_frame_tests {
    use super::*;
    use impg::syng::SyncmerParams;

    /// A deterministic ACGT generator (reproducible, no external deps).
    fn deterministic_sequence(seed: u64, len: usize) -> Vec<u8> {
        let mut state = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let mut out = Vec::with_capacity(len);
        while out.len() < len {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            out.push(b"ACGT"[(state >> 33) as usize % 4]);
        }
        out
    }

    #[test]
    fn the_global_frame_routes_what_the_local_substrate_routes() {
        // Small syncmers (k=4 inner, w=7 window -> 11bp syncmers) over
        // a mini "chromosome": the territory record is a mid-path slice
        // (the edge-band classes live at the record's head/tail; the
        // read is interior, the agreeing domain).
        let params = SyncmerParams { k: 4, w: 7, seed: 7 };
        let k = (params.k + params.w) as u64;
        let origin = 400u64;
        let record_len = 1000u64;
        let chromosome = deterministic_sequence(0x5eed, 2000);
        let record = chromosome[origin as usize..(origin + record_len) as usize].to_vec();
        let record_name = format!("A#0#chr1:{origin}-{}", origin + record_len);
        let global_name = "A#0#chr1";

        let global_panel = SyngIndex::build(
            params,
            vec![(global_name.to_string(), chromosome.clone())].into_iter(),
        );
        let local_panel =
            SyngIndex::build(params, vec![(record_name.clone(), record.clone())].into_iter());
        // the syng metadata-v2 discipline: both panels built by the
        // same fresh builder in-process.

        // the territory rows: the local frame's window grid (the
        // member row in the record's own coordinates) and its global
        // translation (the record header IS the frame statement).
        let window = (100u64, 900u64);
        let local_rows = vec![UniverseRow {
            partition: 0,
            path_name: record_name.clone(),
            start: window.0,
            end: window.1,
        }];
        let translated = territory_frame::translate_row(
            &record_name,
            window.0,
            window.1,
            0,
        )
        .unwrap();
        let global_rows = vec![UniverseRow {
            partition: 0,
            path_name: translated.path_name.clone(),
            start: translated.start,
            end: translated.end,
        }];
        assert_eq!(translated.start, origin + window.0);
        assert_eq!(translated.end, origin + window.1);

        // the lane fetches: the local frame serves the record slice
        // (the locus.fa bytes); the global frame serves the panel-path
        // slice by coordinate (the AGC's role, the same bytes).
        let mut local_sources: HashMap<String, Vec<u8>> = HashMap::new();
        local_sources.insert(record_name.clone(), record.clone());
        let local_fetch = lane_fetch_over(&local_sources);
        let mut path_sources: HashMap<String, Vec<u8>> = HashMap::new();
        path_sources.insert(global_name.to_string(), chromosome.clone());
        let global_fetch = lane_fetch_over(&path_sources);

        let local_territory =
            build_territory_index(&local_panel, &local_rows, &local_fetch, k).unwrap();
        let global_territory =
            build_territory_index(&global_panel, &global_rows, &global_fetch, k).unwrap();
        assert_eq!(local_territory.territories.len(), 1);
        assert_eq!(global_territory.territories.len(), 1);

        // an INTERIOR read (its every anchor is in the two grids'
        // agreeing domain: away from the record head/tail and the
        // window edges by more than a syncmer window).
        let read = &record[300 as usize..(300 + 120) as usize];
        let local_records = derive_read(&local_panel, read).unwrap();
        let global_records = derive_read(&global_panel, read).unwrap();
        assert!(!local_records.is_empty());

        let local_nodes: HashSet<u32> = local_territory
            .entries
            .iter()
            .map(|&(node, _, _)| node.unsigned_abs())
            .collect();
        let global_nodes: HashSet<u32> = global_territory
            .entries
            .iter()
            .map(|&(node, _, _)| node.unsigned_abs())
            .collect();

        for ((local_tokens, _), (global_tokens, _)) in
            local_records.iter().zip(global_records.iter())
        {
            let local_anchors = decode_tokens(local_tokens).unwrap();
            let global_anchors = decode_tokens(global_tokens).unwrap();
            // the same read content: the same anchor structure
            assert_eq!(local_anchors.len(), global_anchors.len());
            for (&(_, local_pos), &(_, global_pos)) in
                local_anchors.iter().zip(global_anchors.iter())
            {
                assert_eq!(local_pos, global_pos);
            }
            // the local touch (remapped onto the local grid) and the
            // global touch agree: the record touches the window in
            // both frames
            let local_touch = local_anchors
                .iter()
                .any(|&(node, _)| local_nodes.contains(&node.unsigned_abs()));
            let global_touch = global_anchors
                .iter()
                .any(|&(node, _)| global_nodes.contains(&node.unsigned_abs()));
            assert_eq!(local_touch, global_touch);
            if !global_touch {
                continue;
            }
            // THE ROUTING EQUIVALENCE: route the pattern against each
            // frame's territory index; the occurrences are the same
            // modulo the coordinate rename (path, start+origin) and
            // the window id.
            let (local_occ, local_total, _, _) =
                route_record(&local_anchors, &local_territory, k);
            let (global_occ, global_total, _, _) =
                route_record(&global_anchors, &global_territory, k);
            assert_eq!(local_total, global_total);
            assert_eq!(local_occ.len(), global_occ.len());
            for ((partition, local_mult), (gpartition, global_mult)) in
                local_occ.iter().zip(global_occ.iter())
            {
                assert_eq!(partition, gpartition);
                assert_eq!(local_mult, global_mult);
            }
        }

        // THE OCCURRENCE POSITIONS: route_record returns them in its
        // positional side channels; verify through the public routing
        // seam: the global frame's occurrences translated by -origin
        // equal the local frame's (the conversion is exactly the
        // coordinate rename).
        for ((local_tokens, _), (global_tokens, _)) in
            local_records.iter().zip(global_records.iter())
        {
            let local_anchors = decode_tokens(local_tokens).unwrap();
            let global_anchors = decode_tokens(global_tokens).unwrap();
            if !global_anchors.iter().any(|&(node, _)| {
                global_nodes.contains(&node.unsigned_abs())
            }) {
                continue;
            }
            let (_, _, local_forward, local_reverse) =
                route_record(&local_anchors, &local_territory, k);
            let (_, _, global_forward, global_reverse) =
                route_record(&global_anchors, &global_territory, k);
            let translate =
                |positions: &[(usize, u64)]| -> Vec<(usize, u64)> {
                    positions
                        .iter()
                        .map(|&(path, start)| (path, start + origin))
                        .collect()
                };
            assert_eq!(translate(&local_forward), global_forward);
            assert_eq!(translate(&local_reverse), global_reverse);
        }
    }

    /// THE ONE-BASE-LATER SUBSTRATE CLASS (the committed L1-L4
    /// convention, the uniform delta=1 AGC-identity receipts): the
    /// substrate's locus.fa content sits ONE BASE LATER than its
    /// record headers. The committed local frame is closed (its
    /// coordinates are slice-relative to the CONTENT, so it never
    /// sees the difference); the global frame must be TOLD the
    /// convention — with the stated `--record-origin 1` the
    /// conversion property holds exactly (the occurrences are the
    /// coordinate rename by `origin + 1`), and the emission-side
    /// back-translation cancels it (the round trip through the
    /// frame table is the identity).
    #[test]
    fn the_one_base_later_substrate_routes_at_the_stated_content_origin() {
        let params = SyncmerParams { k: 4, w: 7, seed: 7 };
        let k = (params.k + params.w) as u64;
        let origin = 400u64;
        let record_len = 1000u64;
        let chromosome = deterministic_sequence(0x5eed, 2000);
        // THE SUBSTRATE CONVENTION: the header states [origin, origin
        // + record_len) but the extraction carried the one-base-later
        // content [origin + 1, origin + 1 + record_len).
        let record = chromosome[(origin + 1) as usize..(origin + 1 + record_len) as usize].to_vec();
        let record_name = format!("A#0#chr1:{origin}-{}", origin + record_len);
        let global_name = "A#0#chr1";
        let record_origin = 1u64;

        let global_panel = SyngIndex::build(
            params,
            vec![(global_name.to_string(), chromosome.clone())].into_iter(),
        );
        // the committed local frame's own panel: built over the
        // CONTENT (the locus.fa bytes), slice-relative — closed.
        let local_panel =
            SyngIndex::build(params, vec![(record_name.clone(), record.clone())].into_iter());

        let window = (100u64, 900u64);
        let local_rows = vec![UniverseRow {
            partition: 0,
            path_name: record_name.clone(),
            start: window.0,
            end: window.1,
        }];
        // the global frame at the STATED content origin: the row
        // lands at origin + 1 + window (the true content position).
        let translated = territory_frame::translate_row(
            &record_name,
            window.0,
            window.1,
            record_origin,
        )
        .unwrap();
        assert_eq!(translated.start, origin + record_origin + window.0);
        assert_eq!(translated.end, origin + record_origin + window.1);
        let global_rows = vec![UniverseRow {
            partition: 0,
            path_name: translated.path_name.clone(),
            start: translated.start,
            end: translated.end,
        }];

        let mut local_sources: HashMap<String, Vec<u8>> = HashMap::new();
        local_sources.insert(record_name.clone(), record.clone());
        let local_fetch = lane_fetch_over(&local_sources);
        let mut path_sources: HashMap<String, Vec<u8>> = HashMap::new();
        path_sources.insert(global_name.to_string(), chromosome.clone());
        let global_fetch = lane_fetch_over(&path_sources);

        let local_territory =
            build_territory_index(&local_panel, &local_rows, &local_fetch, k).unwrap();
        let global_territory =
            build_territory_index(&global_panel, &global_rows, &global_fetch, k).unwrap();

        // an INTERIOR read of the substrate's own content.
        let read = &record[300usize..(300 + 120)];
        let local_records = derive_read(&local_panel, read).unwrap();
        let global_records = derive_read(&global_panel, read).unwrap();
        assert!(!local_records.is_empty());

        let global_nodes: HashSet<u32> = global_territory
            .entries
            .iter()
            .map(|&(node, _, _)| node.unsigned_abs())
            .collect();

        for ((local_tokens, _), (global_tokens, _)) in
            local_records.iter().zip(global_records.iter())
        {
            let local_anchors = decode_tokens(local_tokens).unwrap();
            let global_anchors = decode_tokens(global_tokens).unwrap();
            if !global_anchors
                .iter()
                .any(|&(node, _)| global_nodes.contains(&node.unsigned_abs()))
            {
                continue;
            }
            let (_, _, local_forward, local_reverse) =
                route_record(&local_anchors, &local_territory, k);
            let (_, _, global_forward, global_reverse) =
                route_record(&global_anchors, &global_territory, k);
            // THE ROUTING EQUIVALENCE AT THE STATED ORIGIN: the global
            // occurrences are the local ones renamed by origin +
            // record_origin (the true content position).
            let translate =
                |positions: &[(usize, u64)]| -> Vec<(usize, u64)> {
                    positions
                        .iter()
                        .map(|&(path, start)| (path, start + origin + record_origin))
                        .collect()
                };
            assert_eq!(translate(&local_forward), global_forward);
            assert_eq!(translate(&local_reverse), global_reverse);
        }

        // THE EMISSION CANCEL: a product row at a true content position
        // back-translates through the SAME stated convention to the
        // substrate's own slice-relative extent, and the record NAME
        // keeps the header's own interval (the committed output
        // shape) — the origin cancels exactly.
        let frame =
            impg::genome_inference::territory_frame::TerritoryFrame::from_local_names(
                [record_name.as_str()],
                record_origin,
            )
            .unwrap();
        let product = territory_frame::translate_row(&record_name, 250, 350, record_origin)
            .unwrap();
        let (name, start, end) = frame
            .local_row(&product.path_name, product.start, product.end)
            .unwrap();
        assert_eq!(name, record_name);
        assert_eq!((start, end), (250, 350));
    }
}

// ---------------------------------------------------------------------------
// THE FORK-COW RECEIPTS' OWN PARSER PROOF: the smaps parsing the
// COW-hold gate reads must be exact on the /proc format (header
// lines, field sections, the big-mapping filter) — a pure-function
// test, no /proc read.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod fork_cow_tests {
    use super::*;

    #[test]
    fn fork_smaps_parsing_keeps_the_big_mappings_and_their_fields() {
        let smaps = "\
7f2a00000000-7f2a30000000 rw-p 00000000 00:00 0                                  \n\
Size:           1048576 kB\n\
Rss:            1048576 kB\n\
Pss:             1048576 kB\n\
Shared_Clean:         0 kB\n\
Shared_Dirty:         0 kB\n\
Private_Clean:        0 kB\n\
Private_Dirty:  1048576 kB\n\
7f2b00000000-7f2b01000000 rw-p 00000000 00:00 0                                  \n\
Size:             65536 kB\n\
Rss:              65536 kB\n\
Pss:               32768 kB\n\
Shared_Clean:     32768 kB\n\
Shared_Dirty:         0 kB\n\
Private_Clean:        0 kB\n\
Private_Dirty:        0 kB\n\
7f2c00000000-7f2c30000000 r--p 00000000 fd:01 123456  /data/panel.syng.pstep\n\
Size:           1048577 kB\n\
Rss:              10000 kB\n\
Pss:               5000 kB\n\
Shared_Clean:     10000 kB\n\
Shared_Dirty:         0 kB\n\
Private_Clean:        0 kB\n\
Private_Dirty:        0 kB\n";
        // The header-line predicate: hex address ranges only.
        assert!(fork_is_smaps_header(
            "7f2a00000000-7f2a30000000 rw-p 00000000 00:00 0    "
        ));
        assert!(!fork_is_smaps_header("Size:           1048576 kB"));
        assert!(!fork_is_smaps_header("not-a-range rw-p 00000000 00:00 0"));

        let big = fork_smaps_big_mappings(smaps, 1 << 20);
        // The 1GiB anon arena and the just-over-1GiB file mapping are
        // kept; the 64MiB mapping is not.
        assert_eq!(big.len(), 2);
        assert_eq!(big[0].size_kb, 1 << 20);
        assert_eq!(big[0].private_kb(), 1 << 20);
        assert_eq!(big[0].shared_kb(), 0);
        assert_eq!(big[0].name, "[anon]");
        assert_eq!(big[1].size_kb, (1 << 20) + 1);
        assert_eq!(big[1].rss_kb, 10_000);
        assert_eq!(big[1].pss_kb, 5_000);
        assert_eq!(big[1].shared_kb(), 10_000);
        assert_eq!(big[1].name, "/data/panel.syng.pstep");

        // The field parser matches whole keys only (Shared_* vs
        // Shared_Clean, "Rss:" vs "Rss_Anon:").
        let section = "Rss_Anon:         1 kB\nRss:             2 kB\n";
        assert_eq!(fork_smaps_field(section, "Rss:"), 2);
        assert_eq!(fork_smaps_field(section, "Rss_Anon:"), 1);
        assert_eq!(fork_smaps_field(section, "Pss:"), 0);
    }
}

// ---------------------------------------------------------------------------
// THE CENSUS-PROFILE RESOLUTION PROOF (the owner's Ruling 2, executed:
// the restricted read set is the global-frame call's default census
// profile; the global-native form is the record opt-out). The
// resolution is a pure function of (frame, profile, basis) — every
// reachable state receipted: the local frame untouched, the global
// default fail-closed without a basis, the basis override, and the
// global-native record form.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod census_profile_tests {
    use super::*;

    #[test]
    fn the_local_frame_is_untouched_in_both_profiles() {
        // The local frame keeps the committed behavior: the profile is
        // not consulted (the read-set instrument still restricts there
        // exactly as before, by its own env).
        for profile in [CensusProfile::Restricted, CensusProfile::GlobalNative] {
            for basis in [None, Some("committed.fasta")] {
                let (effective, over) =
                    CensusProfile::resolve(CutFrame::Local, profile, basis).unwrap();
                assert_eq!(effective, profile);
                assert!(!over);
            }
        }
    }

    #[test]
    fn the_global_default_requires_the_basis_fail_closed() {
        let err = CensusProfile::resolve(
            CutFrame::Global,
            CensusProfile::Restricted,
            None,
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("IMPG_CUT_READSET_FASTA"),
            "the fail-closed error must name the basis instrument: {err}"
        );
        assert!(
            err.to_string().contains("global-native"),
            "the fail-closed error must name the record opt-out: {err}"
        );
    }

    #[test]
    fn the_global_restricted_profile_runs_restricted_with_a_basis() {
        let (effective, over) = CensusProfile::resolve(
            CutFrame::Global,
            CensusProfile::Restricted,
            Some("committed.fasta"),
        )
        .unwrap();
        assert_eq!(effective, CensusProfile::Restricted);
        assert!(!over);
    }

    #[test]
    fn the_global_native_record_form_is_the_full_mass() {
        let (effective, over) = CensusProfile::resolve(
            CutFrame::Global,
            CensusProfile::GlobalNative,
            None,
        )
        .unwrap();
        assert_eq!(effective, CensusProfile::GlobalNative);
        assert!(!over);
    }

    #[test]
    fn the_basis_overrides_the_global_native_profile() {
        // IMPG_CUT_READSET_FASTA is the restriction basis override it
        // always was: a stated basis restricts the scan even under the
        // global-native profile, and the resolution says so.
        let (effective, over) = CensusProfile::resolve(
            CutFrame::Global,
            CensusProfile::GlobalNative,
            Some("committed.fasta"),
        )
        .unwrap();
        assert_eq!(effective, CensusProfile::Restricted);
        assert!(over);
    }
}
