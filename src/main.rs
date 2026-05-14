#![allow(clippy::too_many_arguments)]
#![allow(clippy::type_complexity)]
use clap::Parser;
use coitrees::{Interval, IntervalTree};
use impg::alignment_record::{AlignmentFormat, AlignmentRecord, Strand};
use impg::commands::{graph, lace, partition, refine, similarity};
use impg::impg::{AdjustedInterval, CigarOp, Impg};
use impg::impg_index::{ImpgIndex, ImpgWrapper};
use impg::multi_impg::MultiImpg;
use impg::onealn::OneAlnParser;
use impg::seqidx::SequenceIndex;
use impg::sequence_index::{SequenceIndex as SeqIndexTrait, UnifiedSequenceIndex};
use impg::subset_filter::{apply_subset_filter, load_subset_filter, SubsetFilter};
use impg::tpa_parser::TpaParser;
use log::{debug, error, info, warn};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rustc_hash::FxHashMap;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::num::NonZeroUsize;
#[cfg(unix)]
use std::os::fd::{AsRawFd, RawFd};
use std::path::Path;
use std::time::{Duration, Instant};

/// Parse a size value with optional k/m/g suffix (case-insensitive)
/// Examples: "100" -> 100, "10k" -> 10000, "5M" -> 5000000, "1g" -> 1000000000
fn parse_size(s: &str) -> Result<u64, String> {
    let s = s.trim();
    if s.is_empty() {
        return Err("empty value".to_string());
    }

    let (num_str, multiplier) = match s.chars().last() {
        Some('k') | Some('K') => (&s[..s.len() - 1], 1_000u64),
        Some('m') | Some('M') => (&s[..s.len() - 1], 1_000_000u64),
        Some('g') | Some('G') => (&s[..s.len() - 1], 1_000_000_000u64),
        _ => (s, 1u64),
    };

    num_str
        .parse::<u64>()
        .map(|n| n * multiplier)
        .map_err(|e| format!("invalid number '{}': {}", num_str, e))
}

fn resolve_syng_syncmer_params(
    legacy_syncmer_k: Option<u32>,
    smer_length: Option<u32>,
    legacy_syncmer_w: Option<u32>,
    syncmer_length: Option<u32>,
    seed: u32,
) -> io::Result<impg::syng::SyncmerParams> {
    let s = smer_length.or(legacy_syncmer_k).unwrap_or(8);
    if s == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "syng smer length must be greater than 0",
        ));
    }

    let total_k = match (syncmer_length, legacy_syncmer_w) {
        (Some(total), None) => total,
        (None, Some(w)) => s.checked_add(w).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "syng syncmer length overflowed u32")
        })?,
        (None, None) => 63,
        (Some(_), Some(_)) => unreachable!("clap conflicts_with prevents this"),
    };

    if total_k == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "syng syncmer length must be greater than 0",
        ));
    }
    if total_k % 2 == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("syng syncmer length must be odd, got {}", total_k),
        ));
    }
    if total_k <= s {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "syng syncmer length ({}) must be greater than smer length ({})",
                total_k, s
            ),
        ));
    }

    Ok(impg::syng::SyncmerParams {
        k: s,
        w: total_k - s,
        seed,
    })
}

#[derive(Debug, Clone)]
struct SyngAgcRecord {
    sample: String,
    contig: String,
    name: String,
    length: usize,
}

fn syng_sequence_name(sample: &str, contig: &str) -> String {
    if contig.contains('#') {
        contig.to_string()
    } else {
        format!("{}@{}", contig, sample)
    }
}

fn ragc_numeric_to_ascii(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        })
        .collect()
}

fn read_syng_fasta_sequences(fasta_path: &str) -> io::Result<Vec<(String, Vec<u8>)>> {
    let file = File::open(fasta_path)?;
    let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
        io::Error::other(format!(
            "Failed to open reader for '{}': {}",
            fasta_path, e
        ))
    })?;
    let reader = BufReader::new(reader);

    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_name.is_empty() && !current_seq.is_empty() {
                sequences.push((
                    std::mem::take(&mut current_name),
                    std::mem::take(&mut current_seq),
                ));
            }
            current_name = line
                .strip_prefix('>')
                .unwrap_or("")
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_seq.clear();
        } else {
            current_seq.extend(line.trim().as_bytes());
        }
    }
    if !current_name.is_empty() && !current_seq.is_empty() {
        sequences.push((current_name, current_seq));
    }

    Ok(sequences)
}

fn collect_syng_agc_records(agc_path: &str) -> io::Result<Vec<SyngAgcRecord>> {
    let config = ragc_core::DecompressorConfig { verbosity: 0 };
    let mut decompressor = ragc_core::Decompressor::open(agc_path, config).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Failed to open AGC file '{}': {}", agc_path, e),
        )
    })?;

    let samples = decompressor.list_samples();
    let mut records = Vec::new();
    for sample in samples {
        let contigs = decompressor
            .list_contigs_names_only(&sample)
            .unwrap_or_default();
        for contig in contigs {
            let length = decompressor.get_contig_length(&sample, &contig).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Failed to get length for contig '{}@{}': {}", contig, sample, e),
                )
            })?;
            records.push(SyngAgcRecord {
                name: syng_sequence_name(&sample, &contig),
                sample: sample.clone(),
                contig,
                length,
            });
        }
    }
    Ok(records)
}

fn read_syng_agc_record(
    decompressor: &mut ragc_core::Decompressor,
    record: &SyngAgcRecord,
) -> io::Result<Vec<u8>> {
    let contig_data = decompressor
        .get_contig_range(&record.sample, &record.contig, 0, record.length)
        .map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Failed to decompress contig '{}@{}': {}",
                    record.contig, record.sample, e
                ),
            )
        })?;
    Ok(ragc_numeric_to_ascii(&contig_data))
}

fn format_duration(duration: Duration) -> String {
    if duration.as_secs() >= 60 {
        let mins = duration.as_secs() / 60;
        let secs = duration.as_secs_f64() - (mins * 60) as f64;
        format!("{}m{:.3}s", mins, secs)
    } else if duration.as_secs() > 0 {
        format!("{:.3}s", duration.as_secs_f64())
    } else if duration.as_millis() > 0 {
        format!("{:.3}ms", duration.as_secs_f64() * 1_000.0)
    } else {
        format!("{}us", duration.as_micros())
    }
}

fn parse_sequence_record_name(line: &str, marker: char) -> String {
    line.strip_prefix(marker)
        .unwrap_or("")
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

fn read_query_sequences(path: &str) -> io::Result<Vec<(String, Vec<u8>)>> {
    let file = File::open(path)?;
    let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
        io::Error::other(format!("Failed to open reader for '{}': {}", path, e))
    })?;
    let mut reader = BufReader::new(reader);

    let mut first_line = String::new();
    loop {
        first_line.clear();
        if reader.read_line(&mut first_line)? == 0 {
            return Ok(Vec::new());
        }
        if !first_line.trim().is_empty() {
            break;
        }
    }

    if first_line.starts_with('>') {
        read_fasta_records(reader, first_line)
    } else if first_line.starts_with('@') {
        read_fastq_records(reader, first_line)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("query file '{}' is not FASTA or FASTQ", path),
        ))
    }
}

fn read_fasta_records<R: BufRead>(
    reader: R,
    first_line: String,
) -> io::Result<Vec<(String, Vec<u8>)>> {
    let mut sequences = Vec::new();
    let mut current_name = parse_sequence_record_name(&first_line, '>');
    let mut current_seq = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            sequences.push((
                std::mem::take(&mut current_name),
                std::mem::take(&mut current_seq),
            ));
            current_name = parse_sequence_record_name(&line, '>');
            current_seq.clear();
        } else {
            current_seq.extend(line.trim().as_bytes());
        }
    }
    sequences.push((current_name, current_seq));
    Ok(sequences)
}

fn read_fastq_records<R: BufRead>(
    mut reader: R,
    mut header: String,
) -> io::Result<Vec<(String, Vec<u8>)>> {
    let mut sequences = Vec::new();
    loop {
        if !header.starts_with('@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ record header must start with '@'",
            ));
        }
        let name = parse_sequence_record_name(&header, '@');

        let mut seq = String::new();
        let mut plus = String::new();
        let mut qual = String::new();
        if reader.read_line(&mut seq)? == 0
            || reader.read_line(&mut plus)? == 0
            || reader.read_line(&mut qual)? == 0
        {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ record",
            ));
        }
        if !plus.starts_with('+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ separator line must start with '+'",
            ));
        }
        sequences.push((name, seq.trim().as_bytes().to_vec()));

        header.clear();
        loop {
            if reader.read_line(&mut header)? == 0 {
                return Ok(sequences);
            }
            if !header.trim().is_empty() {
                break;
            }
            header.clear();
        }
    }
}

fn write_syng_map_gaf<W: Write>(
    out: &mut W,
    query_name: &str,
    query_len: u64,
    syncmer_len: u64,
    syncmers: &[impg::syng::SyngQuerySyncmer],
) -> io::Result<()> {
    if syncmers.is_empty() {
        return Ok(());
    }
    let qstart = syncmers.iter().map(|s| s.query_pos).min().unwrap_or(0);
    let qend = syncmers
        .iter()
        .map(|s| s.query_pos.saturating_add(syncmer_len))
        .max()
        .unwrap_or(qstart)
        .min(query_len);
    let mut path = String::new();
    for sm in syncmers {
        let orient = if sm.signed_node >= 0 { '>' } else { '<' };
        path.push(orient);
        path.push_str(&sm.node_id.to_string());
    }
    let path_len = (syncmers.len() as u64).saturating_mul(syncmer_len);
    let matches = path_len.min(qend.saturating_sub(qstart));
    let block_len = qend.saturating_sub(qstart);
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t+\t{}\t{}\t0\t{}\t{}\t{}\t0\tan:i:{}\tsk:i:{}",
        query_name,
        query_len,
        qstart,
        qend,
        path,
        path_len,
        path_len,
        matches,
        block_len,
        syncmers.len(),
        syncmer_len
    )
}

fn write_syng_map_paf<W: Write>(
    out: &mut W,
    hit: &impg::syng::SyngMapHit,
    syncmer_len: u64,
) -> io::Result<()> {
    let matches = (hit.anchors as u64).saturating_mul(syncmer_len);
    let block_len = hit
        .query_end
        .saturating_sub(hit.query_start)
        .max(hit.target_end.saturating_sub(hit.target_start));
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\tan:i:{}\tsk:i:{}",
        hit.query_name,
        hit.query_len,
        hit.query_start,
        hit.query_end,
        hit.strand,
        hit.target_name,
        hit.target_len,
        hit.target_start,
        hit.target_end,
        matches,
        block_len,
        hit.anchors,
        syncmer_len
    )
}

fn emit_syng_map<W: Write>(
    out: &mut W,
    syng_index: &impg::syng::SyngIndex,
    queries: Vec<(String, Vec<u8>)>,
    output_format: &str,
    min_anchors: usize,
    chain_budget: u64,
    max_hits: usize,
) -> io::Result<()> {
    let syncmer_len = (syng_index.params.k + syng_index.params.w) as u64;
    match output_format {
        "gaf" => {
            for (query_name, query_seq) in queries {
                let mut syncmers = syng_index.matched_syncmers_in_sequence(&query_seq);
                if syncmers.len() < min_anchors {
                    continue;
                }
                syncmers.sort_by(|a, b| {
                    a.query_pos
                        .cmp(&b.query_pos)
                        .then(a.node_id.cmp(&b.node_id))
                });
                write_syng_map_gaf(
                    out,
                    &query_name,
                    query_seq.len() as u64,
                    syncmer_len,
                    &syncmers,
                )?;
            }
        }
        "paf" => {
            for (query_name, query_seq) in queries {
                let mut hits =
                    syng_index.map_sequence(&query_name, &query_seq, min_anchors, chain_budget)?;
                if max_hits > 0 && hits.len() > max_hits {
                    hits.truncate(max_hits);
                }
                for hit in hits {
                    write_syng_map_paf(out, &hit, syncmer_len)?;
                }
            }
        }
        other => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("unsupported map output format '{}'; expected 'gaf' or 'paf'", other),
            ));
        }
    }
    Ok(())
}

#[cfg(unix)]
fn silence_stdout_for_process() -> io::Result<RawFd> {
    io::stdout().flush()?;
    unsafe {
        libc::fflush(std::ptr::null_mut());
    }
    let saved_stdout = unsafe { libc::dup(libc::STDOUT_FILENO) };
    if saved_stdout < 0 {
        return Err(io::Error::last_os_error());
    }
    let devnull = std::fs::OpenOptions::new().write(true).open("/dev/null")?;
    if unsafe { libc::dup2(devnull.as_raw_fd(), libc::STDOUT_FILENO) } < 0 {
        let e = io::Error::last_os_error();
        unsafe {
            libc::close(saved_stdout);
        }
        return Err(e);
    }
    Ok(saved_stdout)
}

#[cfg(unix)]
fn write_all_to_fd(fd: RawFd, mut buf: &[u8]) -> io::Result<()> {
    while !buf.is_empty() {
        let written = unsafe { libc::write(fd, buf.as_ptr() as *const libc::c_void, buf.len()) };
        if written < 0 {
            let e = io::Error::last_os_error();
            if e.kind() == io::ErrorKind::Interrupted {
                continue;
            }
            return Err(e);
        }
        if written == 0 {
            return Err(io::Error::new(
                io::ErrorKind::WriteZero,
                "failed to write mapper output to stdout",
            ));
        }
        buf = &buf[written as usize..];
    }
    Ok(())
}

#[cfg(not(unix))]
fn silence_stdout_for_process() -> io::Result<()> {
    Ok(())
}

/// Resolve the `--temp-dir` value to an absolute path.
///
/// - `None` → current working directory (absolute).
/// - `Some("ramdisk")` → `/dev/shm` if available, else cwd.
/// - `Some(path)` → `path` (created if it doesn't exist), canonicalized to absolute.
///
/// The returned path is guaranteed to exist and to be absolute, so downstream
/// code can treat it as a stable, process-wide temp root.
fn resolve_temp_dir(temp_dir: Option<String>) -> io::Result<String> {
    let raw = match temp_dir.as_deref() {
        Some("ramdisk") => {
            let dev_shm = std::path::Path::new("/dev/shm");
            if dev_shm.is_dir() {
                "/dev/shm".to_string()
            } else {
                log::warn!(
                    "--temp-dir ramdisk requested but /dev/shm is not available, using current directory"
                );
                std::env::current_dir()?.to_string_lossy().into_owned()
            }
        }
        Some(other) => other.to_string(),
        None => std::env::current_dir()?.to_string_lossy().into_owned(),
    };

    // Ensure the directory exists before canonicalizing.
    let p = std::path::Path::new(&raw);
    if !p.exists() {
        std::fs::create_dir_all(p).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("failed to create --temp-dir '{}': {}", raw, e),
            )
        })?;
    } else if !p.is_dir() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("--temp-dir '{}' exists but is not a directory", raw),
        ));
    }

    let abs = std::fs::canonicalize(p).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("failed to canonicalize --temp-dir '{}': {}", raw, e),
        )
    })?;
    Ok(abs.to_string_lossy().into_owned())
}

/// Install `temp_dir` as the process-wide temp root for every tool that
/// reads `TMPDIR` (FastGA/GIXmake via fastga-rs, wfmash, gfaffix, Rust's
/// `tempfile` crate) **and** for seqwish's internal tempfile state.
fn setup_temp_dir(temp_dir: &str) -> io::Result<()> {
    // Sanity check
    let p = std::path::Path::new(temp_dir);
    if !p.is_dir() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "setup_temp_dir: '{}' does not exist or is not a directory",
                temp_dir
            ),
        ));
    }

    // 1. Environment: picked up by Rust's tempfile crate, FastGA-rs
    //    (which calls tempfile::Builder::new().tempdir() when no explicit
    //    temp_dir is passed), wfmash, and anything else that honors TMPDIR.
    std::env::set_var("TMPDIR", temp_dir);

    // 2. Seqwish has its own process-wide tempfile directory that must be
    //    configured separately — setting TMPDIR alone is not sufficient.
    seqwish::tempfile::set_dir(temp_dir);

    log::info!("[temp] Using temp directory: {}", temp_dir);
    Ok(())
}

use impg::{EngineOpts, GfaEngine};
use sweepga::knn_graph::SparsificationStrategy;

/// Index mode.
#[derive(Clone, Debug, Default, clap::ValueEnum)]
enum IndexMode {
    /// Auto-select based on file count
    #[default]
    Auto,
    /// Single combined index
    Single,
    /// One index per alignment file
    PerFile,
}

/// Basic common options shared between all commands
#[derive(Parser, Debug)]
#[command(next_help_heading = "General options")]
struct CommonOpts {
    /// Number of threads for parallel processing.
    #[clap(short = 't', long, value_parser, default_value_t = NonZeroUsize::new(4).unwrap())]
    threads: NonZeroUsize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[clap(short, long, default_value = "1")]
    verbose: u8,
}

/// Alignment options shared by all impg commands that run alignments.
///
/// Flattens sweepga's `AlnArgs` and adds one impg-specific filtering
/// knob, `--min-map-length`, which is a post-alignment drop threshold,
/// so it lives here rather than in sweepga. Field reads use `aln.sw.<field>`
/// to avoid collisions with impg-internal configs mirror these names.
#[derive(Parser, Debug, Clone)]
struct AlnOpts {
    #[clap(flatten)]
    sw: sweepga::AlnArgs,

    /// Drop mappings shorter than this before plane-sweep / scaffold filtering. Accepts k/m/g suffixes.
    #[clap(long = "min-map-length", value_parser = parse_size, default_value = "0",
           help_heading = "Basic filtering")]
    min_map_length: u64,
}

/// Seqwish graph induction options shared by the `graph`, `query`, and `partition` commands.
#[derive(Parser, Debug, Clone)]
#[command(next_help_heading = "Graph induction (seqwish)")]
struct SeqwishOpts {
    /// Maximum repeat count for transitive closure (0 = no limit)
    #[clap(long, value_parser, default_value_t = 0)]
    repeat_max: u64,

    /// Minimum distance between repeats
    #[clap(long, value_parser, default_value_t = 0)]
    min_repeat_dist: u64,

    /// Minimum match length filter for alignments
    #[clap(long, value_parser, default_value_t = 23)]
    min_match_len: u64,

    /// Sparse factor for input matches (0.0 = keep all)
    #[clap(long, value_parser, default_value_t = 0.0)]
    sparse_factor: f32,

    /// Batch size for transitive closure computation
    #[clap(long, value_parser, default_value_t = 10_000_000)]
    transclose_batch: u64,

    /// Use disk-backed interval trees (slower but lower memory)
    #[clap(long, action)]
    disk_backed: bool,
}

/// Smoothxg-style graph smoothing options shared by the `graph`, `query`, and `partition` commands.
#[derive(Parser, Debug, Clone)]
#[command(next_help_heading = "Graph smoothing (pggb engine)")]
struct SmoothOpts {
    /// Comma-separated target POA block lengths in bp, one per smoothing pass
    /// (default: "700,1100" matching pggb's `-G 700,1100`)
    #[clap(long, value_parser, default_value = "700,1100")]
    target_poa_length: String,

    /// Maximum node length before chopping (default: 100)
    #[clap(long, value_parser, default_value_t = 100)]
    max_node_length: usize,

    /// Fraction of avg block length to use as POA padding (default: 0.001)
    #[clap(long, value_parser, default_value_t = 0.001)]
    poa_padding_fraction: f64,
}

impl SmoothOpts {
    /// Parse `--target-poa-length` into a `Vec<usize>`.
    fn parse_target_poa_lengths(&self) -> io::Result<Vec<usize>> {
        self.target_poa_length
            .split(',')
            .map(|s| {
                s.trim().parse::<usize>().map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "--target-poa-length: '{}' is not a valid positive integer",
                            s.trim()
                        ),
                    )
                })
            })
            .collect()
    }
}

/// Engine + graph-building options shared by query, partition, similarity, and graph.
#[derive(Parser, Debug)]
struct EngineCliOpts {
    /// GFA engine: 'pggb' (default), 'seqwish', or 'poa'.
    /// Append ':WINDOW' to enable partitioned mode, e.g. 'pggb:10000'
    /// splits into 10kb windows, builds per-window, laces, and normalizes.
    #[arg(help_heading = "Output options")]
    #[clap(
        long = "gfa-engine",
        default_value = "pggb",
        value_name = "ENGINE[:WINDOW]"
    )]
    engine_raw: String,

    /// POA alignment scores as match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2
    #[arg(help_heading = "Alignment options")]
    #[clap(long, value_parser, default_value = "5,4,6,2,24,1")]
    poa_scoring: String,

    #[clap(flatten)]
    aln: AlnOpts,

    #[clap(flatten)]
    seqwish: SeqwishOpts,

    #[clap(flatten)]
    smooth: SmoothOpts,

    /// Save intermediate debug files (PAFs, FASTAs, sub-GFAs, chunk info) to this directory.
    #[clap(long, value_parser)]
    debug_dir: Option<String>,
}

impl EngineCliOpts {
    /// Parse `--gfa-engine` value into (GfaEngine, Option<partition_size>).
    ///
    /// Accepted forms: `pggb`, `seqwish`, `poa`, `pggb:10000`, `seqwish:5000`, etc.
    /// A bare colon (`pggb:`) is an error.
    fn parse_engine(&self) -> io::Result<(GfaEngine, Option<usize>)> {
        let raw = self.engine_raw.trim();
        let (name, partition_size) = if let Some(idx) = raw.find(':') {
            let engine_str = &raw[..idx];
            let ps_str = &raw[idx + 1..];
            if ps_str.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Invalid --gfa-engine '{}': expected a window size after ':', e.g. '{}:10000'",
                        raw, engine_str
                    ),
                ));
            }
            let ps: usize = ps_str.parse().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Invalid --gfa-engine '{}': '{}' is not a valid window size (expected integer)",
                        raw, ps_str
                    ),
                )
            })?;
            if ps < 1_000 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Invalid --gfa-engine '{}': window size must be at least 1000 bp",
                        raw
                    ),
                ));
            }
            (engine_str, Some(ps))
        } else {
            (raw, None)
        };

        let engine = match name {
            "pggb" => GfaEngine::Pggb,
            "seqwish" => GfaEngine::Seqwish,
            "poa" => GfaEngine::Poa,
            "syng-native" | "syng_native" | "syng" => GfaEngine::SyngNative,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Unknown GFA engine '{}'. Valid engines: pggb, seqwish, poa, syng-native",
                        name
                    ),
                ));
            }
        };

        Ok((engine, partition_size))
    }

    /// Parse POA scoring parameters.
    fn parse_poa_scoring(&self) -> io::Result<(u8, u8, u8, u8, u8, u8)> {
        parse_poa_scoring_string(&self.poa_scoring)
    }

    /// Validate that CLI options are compatible with the selected engine.
    fn validate_engine_params(&self, engine: GfaEngine) -> io::Result<()> {
        match engine {
            GfaEngine::Poa => {
                // The POA engine builds the graph directly from input
                // sequences withoyut calling an external pairwise aligner
                // (FastGA/wfmash), so it has no plane-sweep/scaffold filter stage.

                if self.aln.sw.sparsify != SparsificationStrategy::None {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--sparsify controls external-aligner pair selection; \
                         --gfa-engine poa has no pair-selection step",
                    ));
                }
                if self.aln.sw.no_filter {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--no-filter disables the post-alignment PAF filter; \
                         --gfa-engine poa has no filtering step",
                    ));
                }
            }
            GfaEngine::Seqwish | GfaEngine::Pggb => {}
            GfaEngine::SyngNative => {
                // Syng-native will generate its own alignments internally
                // from syncmer anchors + BiWFA; the external-aligner knobs
                // don't apply. Sparsification is driven by syng anchor
                // density and sweepga::knn_graph, not --sparsify.
                if self.aln.sw.sparsify != SparsificationStrategy::None {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--sparsify controls external-aligner pair selection; \
                         --gfa-engine syng-native selects pairs from syng anchor counts",
                    ));
                }
                if self.aln.sw.no_filter {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--no-filter disables the post-alignment PAF filter; \
                         --gfa-engine syng-native has no external filter step",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Resolve and build an `EngineOpts`.
    fn build(&self, num_threads: usize) -> io::Result<EngineOpts> {
        let (engine, partition_size) = self.parse_engine()?;
        self.validate_engine_params(engine)?;

        let sparsify = self.aln.sw.sparsify.clone();
        let mash_params = sweepga::knn_graph::MashParams {
            kmer_size: self.aln.sw.mash_kmer_size,
            sketch_size: self.aln.sw.mash_sketch_size,
        };
        let temp_dir = resolve_temp_dir(self.aln.sw.tempdir.clone())?;
        setup_temp_dir(&temp_dir)?;
        build_engine_opts(
            engine,
            num_threads,
            &self.aln,
            sparsify,
            mash_params,
            &self.seqwish,
            &self.smooth,
            self.debug_dir.clone(),
            Some(temp_dir),
            partition_size,
        )
    }
}

/// Alignment file and index options for commands that work with alignments
#[derive(Parser, Debug)]
struct AlignmentOpts {
    /// Path to the alignment files (PAF, 1ALN, or TPA format).
    #[arg(help_heading = "Alignment input")]
    #[clap(short = 'a', long, value_parser, required = false, num_args = 1.., conflicts_with = "alignment_list")]
    alignment_files: Vec<String>,

    /// Path to a text file containing paths to alignment files (one per line, PAF, 1ALN, or TPA format).
    #[arg(help_heading = "Alignment input")]
    #[clap(
        long,
        value_parser,
        required = false,
        conflicts_with = "alignment_files"
    )]
    alignment_list: Option<String>,

    /// Path to the IMPG index file (implies single-index mode).
    #[arg(help_heading = "Index options")]
    #[clap(short = 'i', long, value_parser)]
    index: Option<String>,

    /// Force the regeneration of the index, even if it already exists.
    #[arg(help_heading = "Index options")]
    #[clap(short = 'f', long, action)]
    force_reindex: bool,

    /// Index mode: single, per-file or auto (single, but per-file if >= 100 files).
    #[arg(help_heading = "Index options")]
    #[clap(long, value_enum, default_value_t = IndexMode::Auto)]
    index_mode: IndexMode,

    /// Trace spacing for .1aln alignment files (used when converting tracepoints to CIGAR)
    #[arg(help_heading = "Alignment options")]
    #[clap(long, value_parser, default_value = "100")]
    trace_spacing: u32,

    /// Disable bidirectional alignment interpretation (default: bidirectional enabled).
    /// By default, every alignment A->B creates implicit reverse mapping B->A,
    /// guaranteeing all genomes are bridges for transitive queries.
    #[arg(help_heading = "Index options")]
    #[clap(long)]
    unidirectional: bool,
}

/// Sequence file options for commands that need FASTA/AGC files
#[derive(Parser, Debug)]
#[command(next_help_heading = "Sequence input")]
struct SequenceOpts {
    /// List of sequence file paths (FASTA or AGC) (required for 'gfa', 'maf', and 'fasta')
    #[clap(long, value_parser, num_args = 1.., value_delimiter = ' ', conflicts_with_all = &["sequence_list"])]
    sequence_files: Option<Vec<String>>,

    /// Path to a text file containing paths to sequence files (FASTA or AGC) (required for 'gfa', 'maf', and 'fasta')
    #[clap(long, value_parser, conflicts_with_all = &["sequence_files"])]
    sequence_list: Option<String>,
}

impl SequenceOpts {
    /// Resolve sequence files from either --sequence-files or --sequence-list
    fn resolve_sequence_files(&self) -> io::Result<Vec<String>> {
        match (&self.sequence_files, &self.sequence_list) {
            // Handle --sequence-files option
            (Some(files), None) => Ok(files.clone()),
            // Handle --sequence-list option
            (None, Some(list_file)) => {
                let content = std::fs::read_to_string(list_file).map_err(|e| {
                    io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("Failed to read sequence list file '{list_file}': {e}"),
                    )
                })?;

                Ok(content
                    .lines()
                    .filter(|line| !line.trim().is_empty() && !line.trim().starts_with('#'))
                    .map(|line| line.trim().to_string())
                    .collect())
            }
            (None, None) => Ok(Vec::new()),
            (Some(_), Some(_)) => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Cannot specify both --sequence-files and --sequence-list",
            )),
        }
    }

    /// Build sequence index if files are provided
    fn build_sequence_index(&self) -> io::Result<Option<UnifiedSequenceIndex>> {
        let seq_files = self.resolve_sequence_files()?;

        if seq_files.is_empty() {
            Ok(None)
        } else {
            let file_type = if seq_files.iter().any(|f| f.ends_with(".agc")) {
                "AGC"
            } else {
                "FASTA"
            };
            let num_files = seq_files.len();
            info!(
                "Building {file_type} index for {num_files} file{}",
                if num_files == 1 { "" } else { "s" }
            );

            match UnifiedSequenceIndex::from_files(&seq_files) {
                Ok(index) => Ok(Some(index)),
                Err(e) => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Failed to build sequence index: {e}"),
                )),
            }
        }
    }
}

/// Common sequence and output options for GFA/MAF/FASTA output formats
#[derive(Parser, Debug)]
struct GfaMafFastaOpts {
    #[clap(flatten)]
    sequence: SequenceOpts,

    /// Reverse complement reverse strand sequences (for 'fasta' output)
    #[arg(help_heading = "Output options")]
    #[clap(long, action)]
    reverse_complement: bool,

    /// Force processing of large regions (>10kbp) with maf/gfa output formats
    #[arg(help_heading = "Output options")]
    #[clap(long, action)]
    force_large_region: bool,
}

impl GfaMafFastaOpts {
    /// Helper to validate and setup POA/sequence resources for a given output format, including sequence index for PAF with original coordinates or for 1aln files
    fn setup_output_resources(
        self,
        output_format: &str,
        original_sequence_coordinates: bool,
        alignment_files: &[String],
        approximate_mode: bool,
        poa_scoring: &str,
    ) -> io::Result<(
        Option<UnifiedSequenceIndex>,
        Option<(u8, u8, u8, u8, u8, u8)>,
    )> {
        // Check if any of the alignment files are tracepoint-based (.1aln/.tpa, which require sequence data for tracepoint conversion)
        let has_tracepoint_files = alignment_files
            .iter()
            .any(|f| f.ends_with(".1aln") || f.ends_with(".tpa"));

        // In approximate mode with bed/bedpe output, tracepoint files don't need sequences
        let tracepoint_needs_sequences = has_tracepoint_files
            && !(approximate_mode && (output_format == "bed" || output_format == "bedpe"));

        let needs_sequence_mandatory = matches!(
            output_format,
            "gfa" | "maf" | "fasta" | "fasta-aln" | "fasta+paf" | "gbwt"
        ) || tracepoint_needs_sequences;
        let needs_sequence_optional = output_format == "paf" && original_sequence_coordinates;
        // POA scoring is needed for gfa (all engines may use it), maf, and fasta-aln
        let needs_poa = matches!(output_format, "gfa" | "maf" | "fasta-aln");

        let scoring_params = if needs_poa {
            Some(parse_poa_scoring_string(poa_scoring)?)
        } else {
            None
        };

        let sequence_index = if needs_sequence_mandatory || needs_sequence_optional {
            let index = self.sequence.build_sequence_index()?;
            if index.is_none() && needs_sequence_mandatory {
                let file_types = "FASTA/AGC";

                let msg = if has_tracepoint_files {
                    format!("Sequence files ({file_types}) are required for tracepoint alignment files (.1aln/.tpa) to convert tracepoints to CIGAR strings. Use --sequence-files or --sequence-list")
                } else {
                    format!("Sequence files ({file_types}) are required for '{output_format}' output format. Use --sequence-files or --sequence-list")
                };

                return Err(io::Error::new(io::ErrorKind::InvalidInput, msg));
            }
            index
        } else {
            None
        };

        Ok((sequence_index, scoring_params))
    }
}

/// Transitive query options
#[derive(Parser, Debug, Clone)]
struct TransitiveOpts {
    /// Use Depth-First Search instead of BFS for transitive queries (slower, but returns fewer overlapping results)
    #[arg(help_heading = "Transitive query options")]
    #[clap(long, action)]
    transitive_dfs: bool,

    /// Maximum recursion depth for transitive overlaps (0 for no limit)
    #[arg(help_heading = "Transitive query options")]
    #[clap(short = 'm', long, value_parser, default_value_t = 2)]
    max_depth: u16,

    /// Minimum region size to consider for transitive queries (default: 101; for fastga --approximate, set > trace_spacing)
    #[arg(help_heading = "Transitive query options")]
    #[clap(long, value_parser)]
    min_transitive_len: Option<i32>,

    /// Minimum distance between transitive ranges to consider on the same sequence
    #[arg(help_heading = "Transitive query options")]
    #[clap(long, value_parser, default_value_t = 10)]
    min_distance_between_ranges: i32,
}

impl TransitiveOpts {
    /// Get the effective min_transitive_len value (default: 101)
    fn effective_min_transitive_len(&self) -> i32 {
        self.min_transitive_len.unwrap_or(101)
    }
}

/// Common query and filtering options
#[derive(Parser, Debug, Clone)]
struct QueryOpts {
    /// Target range in the format `seq_name:start-end`
    #[arg(help_heading = "Query region")]
    #[clap(short = 'r', long, value_parser, conflicts_with = "target_bed")]
    target_range: Option<String>,

    /// Path to the BED file containing target regions
    #[arg(help_heading = "Query region")]
    #[clap(short = 'b', long, value_parser, conflicts_with = "target_range")]
    target_bed: Option<String>,

    /// Maximum distance between regions to merge
    #[arg(help_heading = "Filtering and merging")]
    #[clap(
        short = 'd',
        long,
        value_parser,
        conflicts_with = "no_merge",
        default_value_t = 0
    )]
    merge_distance: i32,

    /// Disable merging for all output formats
    #[arg(help_heading = "Filtering and merging")]
    #[clap(long, action, conflicts_with = "merge_distance")]
    no_merge: bool,

    /// Minimum gap-compressed identity threshold (0.0-1.0)
    #[arg(help_heading = "Filtering and merging")]
    #[clap(long = "min-result-identity", value_parser)]
    min_result_identity: Option<f64>,

    /// Minimum output length: filter results shorter than this (bp)
    #[arg(help_heading = "Filtering and merging")]
    #[clap(short = 'l', long, value_parser)]
    min_output_length: Option<i32>,

    /// Path to a file listing sequence names to include (one per line)
    #[arg(help_heading = "Filtering and merging")]
    #[clap(long, value_parser)]
    subset_sequence_list: Option<String>,

    /// Enable transitive queries (uses BFS by default, or DFS with --transitive-dfs)
    #[arg(help_heading = "Transitive queries")]
    #[clap(short = 'x', long, action)]
    transitive: bool,

    #[clap(flatten)]
    transitive_opts: TransitiveOpts,

    /// Update coordinates to original sequences when input sequences are subsequences (seq_name:start-end) for 'bed', 'bedpe', and 'paf'
    #[arg(help_heading = "Coordinate options")]
    #[clap(long, action)]
    original_sequence_coordinates: bool,

    /// Use approximate mode for faster queries with tracepoint files (.1aln/.tpa, only bed/bedpe output)
    #[arg(help_heading = "Performance")]
    #[clap(long, action)]
    approximate: bool,

    /// Consider strandness when merging output intervals (defaults: merge strands for bed/gfa/maf; keep separate for fasta/fasta-aln)
    #[arg(help_heading = "Filtering and merging")]
    #[clap(long, action)]
    consider_strandness: bool,
}

impl QueryOpts {
    /// Get effective merge distance (-1 if merging is disabled)
    fn effective_merge_distance(&self) -> i32 {
        if self.no_merge {
            -1
        } else {
            self.merge_distance
        }
    }

    /// Whether merged intervals should collapse opposite strands for a given output format
    fn merge_strands_for_output(&self, output_format: &str) -> bool {
        // Default behavior per output format
        let default = match output_format {
            "fasta" | "fasta-aln" => false,
            "maf" | "gfa" | "bed" => true,
            _ => true,
        };

        if self.consider_strandness {
            // When considering strandness, keep strands separate
            false
        } else {
            default
        }
    }
}

/// Refinement-specific options
#[derive(Parser, Debug)]
#[command(next_help_heading = "Refinement options")]
struct RefineOpts {
    #[clap(flatten)]
    query: QueryOpts,

    #[clap(flatten)]
    sequence: SequenceOpts,

    /// Minimum number of bases that supporting samples must span at each region boundary
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_parser, default_value_t = 1000)]
    span_bp: i32,

    /// Maximum per-side extension explored when maximizing boundary support.
    /// Values <= 1 are treated as fractions of the locus length; values > 1 as absolute bp.
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_parser, default_value_t = 0.5)]
    max_extension: f64,

    /// PanSN aggregation mode when counting support (sample/haplotype)
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_enum)]
    pansn_mode: Option<sweepga::pansn::PanSnLevel>,

    /// Step size for expanding flanks (bp)
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_parser, default_value_t = 1000)]
    extension_step: i32,

    /// Optional BED file capturing the entities that span the refined region
    #[arg(help_heading = "Output options")]
    #[clap(long, value_parser)]
    support_output: Option<String>,

    /// BED file with ranges to blacklist when counting maximum possible entities
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_parser)]
    blacklist_bed: Option<String>,
}

impl RefineOpts {
    fn validate(&self) -> io::Result<()> {
        if self.span_bp < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--span-bp must be >= 0",
            ));
        }
        if self.max_extension < 0.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--max-extension must be >= 0",
            ));
        }
        if self.extension_step <= 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--extension-step must be > 0",
            ));
        }

        if self.query.approximate {
            validate_approximate_mode_min_length(
                self.query.transitive_opts.min_transitive_len,
                self.query.transitive,
            )?;
        }

        Ok(())
    }
}

/// Parse sequence name to extract subsequence coordinates and original sequence name
/// Format: `seq_name:start-end` -> (original_seq_name, start_offset)
/// Detect whether a path points to a syng index. Returns the syng
/// prefix (base name before `.1khash` / `.1gbwt`) if one is found.
///
/// Accepts:
///   - Explicit file path ending in `.1khash`, `.1gbwt`, `.spos`, `.pstep`, `.names`, or `.meta` → strip suffix
///   - Bare prefix with sibling `.1khash` on disk → return prefix as-is
fn detect_syng_prefix(path: &str) -> Option<String> {
    if let Some(stem) = path.strip_suffix(".1khash") {
        return Some(stem.to_string());
    }
    if let Some(stem) = path.strip_suffix(".1gbwt") {
        return Some(stem.to_string());
    }
    for suffix in [".spos", ".pstep", ".names", ".meta"] {
        if let Some(stem) = path.strip_suffix(suffix) {
            if std::path::Path::new(&format!("{stem}.1khash")).exists() {
                return Some(stem.to_string());
            }
        }
    }
    for suffix in [".syng.spos", ".syng.pstep", ".syng.names", ".syng.meta"] {
        if let Some(stem) = path.strip_suffix(suffix) {
            if std::path::Path::new(&format!("{stem}.1khash")).exists() {
                return Some(stem.to_string());
            }
            return Some(stem.to_string());
        }
    }
    if std::path::Path::new(&format!("{path}.1khash")).exists() {
        return Some(path.to_string());
    }
    None
}

/// Resolve the effective syng prefix for a command: if the single
/// alignment argument looks like a syng index path, use that.
fn resolve_syng_prefix(alignment: &AlignmentOpts) -> Option<String> {
    if alignment.alignment_files.len() == 1 {
        detect_syng_prefix(&alignment.alignment_files[0])
    } else {
        None
    }
}

fn parse_subsequence_coordinates(seq_name: &str) -> Option<(String, i32)> {
    // Find the last colon to handle formats like "sample#hap#chr:start-end"
    if let Some(colon_pos) = seq_name.rfind(':') {
        let (base_name, range_part) = seq_name.split_at(colon_pos);
        let range_part = &range_part[1..]; // Remove the colon

        // Check if the range part contains a dash
        if let Some(dash_pos) = range_part.find('-') {
            let (start_str, _end_str) = range_part.split_at(dash_pos);

            // Parse the start coordinate
            if let Ok(start_offset) = start_str.parse::<i32>() {
                return Some((base_name.to_string(), start_offset));
            }
        }
    }
    None
}

/// Transform coordinates from subsequence space to original sequence space
fn transform_coordinates_to_original(
    seq_name: &str,
    start: u32,
    end: u32,
    original_coordinates: bool,
) -> (String, u32, u32) {
    if !original_coordinates {
        return (seq_name.to_string(), start, end);
    }

    if let Some((original_name, offset)) = parse_subsequence_coordinates(seq_name) {
        let offset = offset as u32;
        (original_name, start + offset, end + offset)
    } else {
        (seq_name.to_string(), start, end)
    }
}

/// Get the original sequence length when using original_sequence_coordinates
fn get_original_sequence_length(
    original_seq_name: &str,
    external_seq_index: Option<&UnifiedSequenceIndex>,
) -> usize {
    // If we have an external sequence index, try to get the length from it
    if let Some(ext_index) = external_seq_index {
        match ext_index.get_sequence_length(original_seq_name) {
            Ok(length) => return length,
            Err(_) => {
                // Emit warning when sequence not found in index
                warn!(
                    "Sequence '{original_seq_name}' not found in sequence index, using 0 as length"
                );
            }
        }
    } else {
        // Emit warning when no index is provided
        warn!("No sequence index provided, using 0 as length for sequence '{original_seq_name}'");
    }

    0 // Return 0 if the sequence is not found or no index is provided
}

/// Command-line tool for querying pangenome alignment
#[derive(Parser, Debug)]
#[command(author, version, about, disable_help_subcommand = true)]
enum Args {
    /// Create an IMPG index
    Index {
        // --- Input ---
        #[clap(flatten)]
        alignment: AlignmentOpts,

        #[clap(flatten)]
        sequence: SequenceOpts,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Lace files together
    Lace {
        /// List of input files (space-separated)
        #[clap(short = 'f', long, value_parser, num_args = 1.., value_delimiter = ' ', conflicts_with = "file_list")]
        files: Option<Vec<String>>,

        /// Text file containing input file paths (one per line)
        #[clap(short = 'l', long, value_parser, conflicts_with = "files")]
        file_list: Option<String>,

        /// Input file format (gfa, vcf, auto)
        #[clap(long, value_parser, default_value = "auto")]
        format: String,

        #[clap(flatten)]
        sequence: SequenceOpts,

        /// Reference (FASTA or AGC) file for validating contig lengths in VCF files
        #[clap(long, value_parser)]
        reference: Option<String>,

        /// Output file path
        #[clap(short, long, value_parser)]
        output: String,

        /// Output compression format (none, gzip, bgzip, zstd, auto)
        #[clap(long, value_parser, default_value = "auto")]
        compress: String,

        /// Gap filling mode: 0=none, 1=middle gaps only, 2=all gaps (requires --sequence-files or --sequence-list for end gaps, GFA mode only)
        #[clap(long, default_value = "0")]
        fill_gaps: u8,

        /// Skip path range length validation (faster but may miss data integrity issues)
        #[clap(long, default_value = "false")]
        skip_validation: bool,

        /// Directory for temporary files [default: $TMPDIR or cwd; use "ramdisk" for /dev/shm]
        #[clap(long, value_parser)]
        temp_dir: Option<String>,

        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Partition the alignment
    Partition {
        // --- Input ---
        #[clap(flatten)]
        alignment: AlignmentOpts,

        /// Boundary padding in bp for syng queries (default: 120 = 2× syncmer length)
        #[arg(help_heading = "Syng input")]
        #[clap(long, value_parser, default_value_t = 120)]
        syng_padding: u64,

        /// Minimum anchor count required for a syng chain to become a
        /// partition edge. Raising this drops weak/paralog-noise chains
        /// before partition's union-find closes over them, preventing the
        /// transitive-collapse catch-all on TE/subtelomere-rich regions.
        /// 0 disables chain filtering (raw per-syncmer intervals).
        #[arg(help_heading = "Syng input")]
        #[clap(long, value_parser, default_value_t = 0)]
        syng_min_chain_anchors: usize,

        /// Minimum fraction of the query window a syng chain must cover to
        /// become a partition edge. 0.0 disables the extent filter.
        #[arg(help_heading = "Syng input")]
        #[clap(long, value_parser, default_value_t = 0.0)]
        syng_min_chain_fraction: f64,

        /// Disable post-partition sliver rehoming. By default, singleton
        /// partitions (slivers created by greedy masking) are iteratively
        /// reassigned to their flank partitions — the partition that owns
        /// the biologically contiguous context. This is source-agnostic and
        /// runs for both syng- and alignment-backed partitioning.
        #[arg(help_heading = "Partition options")]
        #[clap(long, action)]
        no_rehome_singletons: bool,

        // --- Partition-specific ---
        /// Window size for partitioning
        #[arg(help_heading = "Partition options")]
        #[clap(short = 'w', long, value_parser)]
        window_size: usize,

        /// Path to the file with sequence names to start with (one per line)
        #[arg(help_heading = "Partition options")]
        #[clap(long, value_parser)]
        starting_sequences_file: Option<String>,

        #[arg(help_heading = "Partition options")]
        #[clap(
            long,
            value_parser,
            default_value = "longest",
            help = "Selection mode for next sequence:\n\
                - \"longest\": sequence with longest single missing region\n\
                - \"total\": sequence with highest total missing regions\n\
                - \"sample[,separator]\": sample with highest total missing regions\n\
                - \"haplotype[,separator]\": haplotype highest total missing regions\n\
                The sample/haplotype modes assume PanSN naming; '#' is the default separator."
        )]
        selection_mode: String,

        /// Minimum region size for missing regions
        #[arg(help_heading = "Partition options")]
        #[clap(long, value_parser, default_value_t = 3000)]
        min_missing_size: i32,

        /// Minimum distance from sequence start/end - closer regions will be extended to the boundaries
        #[arg(help_heading = "Partition options")]
        #[clap(long, value_parser, default_value_t = 3000)]
        min_boundary_distance: i32,

        /// Output separate files for each partition when 'bed'
        #[arg(help_heading = "Partition options")]
        #[clap(long, action)]
        separate_files: bool,

        // --- Filtering and merging ---
        /// Maximum distance between regions to merge
        #[arg(help_heading = "Filtering and merging")]
        #[clap(short = 'd', long, value_parser, default_value_t = 100000)]
        merge_distance: i32,

        /// Minimum gap-compressed identity threshold (0.0-1.0)
        #[arg(help_heading = "Filtering and merging")]
        #[clap(long = "min-result-identity", value_parser)]
        min_result_identity: Option<f64>,

        #[clap(flatten)]
        transitive_opts: TransitiveOpts,

        // --- Performance ---
        /// Use approximate mode for faster queries with tracepoint files (.1aln/.tpa, only bed/bedpe output)
        #[arg(help_heading = "Performance")]
        #[clap(long, action)]
        approximate: bool,

        // --- Output ---
        /// Output format: 'bed', 'gfa' (v1.0), 'maf', or 'fasta' ('gfa', 'maf', and 'fasta' require --sequence-files or --sequence-list)
        #[arg(help_heading = "Output options")]
        #[clap(short = 'o', long, value_parser, default_value = "bed")]
        output_format: String,

        /// Output folder for partition files (default: current directory)
        #[arg(help_heading = "Output options")]
        #[clap(long, value_parser)]
        output_folder: Option<String>,

        #[clap(flatten)]
        gfa_maf_fasta: GfaMafFastaOpts,

        #[clap(flatten)]
        engine_cli: EngineCliOpts,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Query overlaps in the alignment
    Query {
        // --- Input ---
        #[clap(flatten)]
        alignment: AlignmentOpts,

        /// Boundary padding in bp for syng queries (default: 120 = 2× syncmer length)
        #[arg(help_heading = "Syng input")]
        #[clap(long, value_parser, default_value_t = 120)]
        syng_padding: u64,

        /// Source-side window extension (bp) for syncmer lookup. Widens the
        /// query interval during syncmer discovery only; boundaries remain
        /// clipped by BiWFA refinement. Helps when the query lands just past
        /// the end of a conserved syncmer block (default: 0)
        #[arg(help_heading = "Syng input")]
        #[clap(long, value_parser, default_value_t = 0)]
        syng_extension: u64,

        /// Boundary-extension budget (bp) per cluster in syng queries.
        /// `refine_homolog_by_alignment` expands the anchor-supported
        /// sub-query and the refinement target window outward by this
        /// amount on each side before BiWFA. Target extension clips to
        /// neighbor-cluster bounds on the same `(genome, strand)`;
        /// source extension clips to the user's query region and to
        /// half the clipped target span (keeps BiWFA's two inputs
        /// similarly-sized, avoiding the O(s²) pathology of a wide
        /// query forced into a narrow target). The default 1 kb
        /// recovers most indel-bounded block boundaries (typical
        /// intra-homology indels are &lt; 500 bp) while staying
        /// tractable on repeat-dense queries with tens of thousands
        /// of clusters. Raise for sparser pangenomes with longer
        /// conserved blocks.
        #[arg(help_heading = "Syng input")]
        #[clap(long, value_parser, default_value_t = 1_000)]
        syng_extend_budget: u64,

        /// Minimum anchor count for a syng chain to be emitted. Chains
        /// with fewer anchors than this are dropped at emission. Default
        /// 2 filters singletons (one-anchor chains with no colinear
        /// mutual-best partner — weakest possible evidence, mostly
        /// noise in repeat-dense regions). Set to 1 to keep singletons
        /// (useful when you want every syncmer match reported).
        #[arg(help_heading = "Syng input")]
        #[clap(long, value_parser, default_value_t = 2)]
        syng_min_chain_anchors: usize,

        /// Minimum chain query-extent as a fraction of the queried
        /// range (0.0 to 1.0). Chains whose anchor span on the query
        /// axis covers less than `fraction × query_range_len` are
        /// dropped. Default 0.0 keeps every chain (maximum paralog
        /// discovery). Set 0.5 to keep only chains covering at least
        /// half the query (good for partitioning — filters out
        /// fragment chains in repeat regions and speeds up hard tiles
        /// dramatically).
        #[arg(help_heading = "Syng input")]
        #[clap(long, value_parser, default_value_t = 0.0)]
        syng_min_chain_fraction: f64,

        /// Debug-only: skip boundary realignment and emit raw syncmer-resolution
        /// intervals from syng's query_region. The default syng path runs
        /// BiWFA boundary realignment for base-pair-precise edges (and iterates
        /// under --transitive). Use --syng-raw to inspect the underlying syng
        /// hits without refinement.
        #[arg(help_heading = "Syng input")]
        #[clap(long, default_value_t = false)]
        syng_raw: bool,

        // --- Query-specific ---
        #[clap(flatten)]
        query: QueryOpts,

        // --- Output ---
        /// Output format: 'auto' ('bed' for -r, 'bedpe' for -b), 'bed', 'bedpe', 'paf', 'gfa', 'maf', 'fasta', 'fasta+paf', or 'gbwt' (region-specific GBWT, requires --sequence-files)
        #[arg(help_heading = "Output options")]
        #[clap(short = 'o', long, value_parser, default_value = "auto")]
        output_format: String,

        /// Prefix for output file (automatically appends the extension based on format; required for 'gbwt' output)
        #[clap(short = 'O', long, value_parser, default_value = None)]
        output_prefix: Option<String>,

        #[clap(flatten)]
        gfa_maf_fasta: GfaMafFastaOpts,

        #[clap(flatten)]
        engine_cli: EngineCliOpts,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Refine loci to maximize the number of samples that span both ends of the region
    Refine {
        // --- Input ---
        #[clap(flatten)]
        alignment: AlignmentOpts,

        // --- Refine-specific ---
        #[clap(flatten)]
        refine: RefineOpts,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Compute pairwise similarity between sequences in a region
    Similarity {
        // --- Input ---
        #[clap(flatten)]
        alignment: AlignmentOpts,

        // --- Query-specific ---
        #[clap(flatten)]
        query: QueryOpts,

        // --- Output ---
        /// Output distances instead of similarities
        #[arg(help_heading = "Output options")]
        #[clap(long, action)]
        distances: bool,

        /// Emit entries for all pairs of groups, including those with zero intersection
        #[arg(help_heading = "Output options")]
        #[clap(long, action, default_value_t = false)]
        all: bool,

        /// The part of each path name before this delimiter is a group identifier
        #[arg(help_heading = "Output options")]
        #[clap(long, value_parser)]
        delim: Option<char>,

        /// Consider the N-th occurrence of the delimiter (1-indexed, default: 1)
        #[arg(help_heading = "Output options")]
        #[clap(long, value_parser, default_value_t = 1)]
        delim_pos: u16,

        #[clap(flatten)]
        gfa_maf_fasta: GfaMafFastaOpts,

        #[clap(flatten)]
        engine_cli: EngineCliOpts,

        // --- PCA ---
        /// Perform PCA/MDS dimensionality reduction on the distance matrix
        #[arg(help_heading = "PCA options")]
        #[clap(long, action)]
        pca: bool,

        /// Number of PCA components to output (default: 2)
        #[arg(help_heading = "PCA options")]
        #[clap(long, value_parser, requires = "pca", default_value_t = 2)]
        pca_components: usize,

        /// Number of previous regions to use for adaptive polarization (0 to disable)
        #[arg(help_heading = "PCA options")]
        #[clap(
            long,
            value_parser,
            requires = "pca",
            conflicts_with = "polarize_guide_samples",
            default_value_t = 3
        )]
        polarize_n_prev: usize,

        /// Comma-separated names of the samples to use for adaptive polarization
        #[arg(help_heading = "PCA options")]
        #[clap(
            long,
            value_parser,
            conflicts_with = "polarize_n_prev",
            value_delimiter = ','
        )]
        polarize_guide_samples: Option<Vec<String>>,

        /// Similarity measure to use for PCA distance matrix ("jaccard", "cosine", or "dice")
        #[arg(help_heading = "PCA options")]
        #[clap(long, value_parser, requires = "pca", default_value = "jaccard")]
        pca_measure: String,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },

    /// Print alignment statistics
    Stats {
        // --- Input ---
        #[clap(flatten)]
        alignment: AlignmentOpts,

        #[clap(flatten)]
        sequence: SequenceOpts,

        // --- Stats-specific ---
        /// List sequence names and lengths (skip overlap statistics)
        #[clap(long)]
        list_sequences: bool,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },

    /// Build a pangenome graph from FASTA sequences using sweepga+seqwish
    Graph {
        // --- Input ---
        #[clap(flatten)]
        fasta_input: SequenceOpts,

        /// Input PAF file (skip alignment step if provided)
        #[clap(short = 'a', long, value_parser)]
        paf_file: Option<String>,

        // --- Output ---
        /// Output GFA file path (use "-" for stdout)
        #[clap(short = 'g', long, value_parser, default_value = "-")]
        output: String,

        #[clap(flatten)]
        engine_cli: EngineCliOpts,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },

    /// Generate alignment pairs with sparsification strategies
    Align {
        // --- Input ---
        #[clap(flatten)]
        fasta_input: SequenceOpts,

        // --- Output ---
        /// Output directory for alignments
        #[clap(short = 'o', long, value_parser, default_value = "alignments")]
        output_dir: String,

        /// Output format: paf, 1aln, or joblist
        #[clap(long, value_parser, default_value = "joblist")]
        format: String,

        /// Execute an existing joblist file (one shell command per line)
        /// in parallel with progress/ETA logging. When set, sparsification
        /// and sequence inputs are ignored.
        #[clap(long, value_parser)]
        run_joblist: Option<String>,

        /// Parallel slots for `--run-joblist`. Each slot runs one command
        /// at a time; per-command thread count is whatever the joblist line
        /// already specifies. Defaults to `--threads`.
        #[clap(long, value_parser)]
        jobs: Option<usize>,

        // --- Alignment ---
        #[clap(flatten)]
        aln: AlnOpts,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },

    /// Map sequences to a syng index using exact shared syncmers
    Map {
        /// Syng index prefix or .1khash/.1gbwt/.spos/.pstep/.names/.meta path
        #[clap(short = 'a', long, value_parser)]
        index: String,

        /// Query FASTA or FASTQ file
        #[clap(short = 'q', long, value_parser)]
        query: String,

        /// Output format: gaf (syncmer-node support) or paf (projected genome coordinates)
        #[clap(short = 'o', long, value_parser, default_value = "paf")]
        output_format: String,

        /// Output file path (default: stdout)
        #[clap(short = 'O', long, value_parser)]
        output: Option<String>,

        /// Minimum shared syncmer anchors required to emit a mapping
        #[clap(long, value_parser, default_value_t = 2)]
        min_anchors: usize,

        /// Anchor chaining budget for PAF projection
        #[clap(long, value_parser, default_value_t = 10000)]
        chain_budget: u64,

        /// Maximum PAF hits to emit per query (0 = no limit)
        #[clap(long, value_parser, default_value_t = 0)]
        max_hits: usize,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },

    /// Build a GBWT syncmer index from sequences
    Syng {
        // --- Input (one required) ---
        /// AGC archive input
        #[clap(long, value_parser, conflicts_with = "fasta")]
        agc: Option<String>,

        /// FASTA input file
        #[clap(short = 'f', long, value_parser, conflicts_with = "agc")]
        fasta: Option<String>,

        // --- Output ---
        /// Output file prefix (produces .1khash, .1gbwt, and syng sidecars; avoids .syng.syng.* if prefix ends in .syng)
        #[clap(short = 'o', long, value_parser)]
        output: String,

        /// Sample every 2^N syncmer steps per path for positional sidecars.
        #[arg(help_heading = "Position sampling")]
        #[clap(long, value_parser, default_value_t = impg::syng::DEFAULT_POSITION_SAMPLE_SHIFT)]
        position_sample_shift: u32,

        /// Compatibility seed stored in sidecars; regular-grid sampling ignores it.
        #[arg(help_heading = "Position sampling")]
        #[clap(long, value_parser, default_value_t = impg::syng::DEFAULT_POSITION_SAMPLE_SEED)]
        position_sample_seed: u64,

        /// Build a deterministic global syncmer dictionary in a parallel prepass, then replay paths through it.
        #[arg(help_heading = "Construction")]
        #[clap(long, action)]
        parallel_dictionary: bool,

        // --- Syncmer parameters ---
        /// Inner k-mer length for syncmer extraction
        #[arg(help_heading = "Syncmer parameters")]
        #[clap(long, value_parser, conflicts_with = "smer_length")]
        syncmer_k: Option<u32>,

        /// Inner smer length for syncmer extraction (`s` in the syng paper; default 8)
        #[arg(help_heading = "Syncmer parameters")]
        #[clap(long, value_parser, conflicts_with = "syncmer_k")]
        smer_length: Option<u32>,

        /// Window length for syncmer extraction (total syncmer length = w + k)
        #[arg(help_heading = "Syncmer parameters")]
        #[clap(long, value_parser, conflicts_with = "syncmer_length")]
        syncmer_w: Option<u32>,

        /// Total syncmer length (`k` in the syng paper; must be odd; default 63)
        #[arg(help_heading = "Syncmer parameters")]
        #[clap(long, value_parser, conflicts_with = "syncmer_w")]
        syncmer_length: Option<u32>,

        /// Hash function seed for syncmer extraction
        #[arg(help_heading = "Syncmer parameters")]
        #[clap(long, value_parser, default_value_t = 7)]
        syncmer_seed: u32,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },

    /// Dump a syng index to GFA (S = syncmer, L = adjacency, P/W = path).
    ///
    /// Loads `<prefix>.1khash`, `<prefix>.1gbwt`, `<prefix>.syng.names`,
    /// `<prefix>.syng.meta` (and `.syng.spos`) and writes a GFA with one
    /// segment per syncmer plus one segment per inter-syncmer gap. Gaps
    /// are filled with real DNA when `--sequence-files` is provided
    /// (sequence names must match the syng path names); otherwise they
    /// are filled with `N`s and a warning is emitted.
    Syng2gfa {
        /// Syng index prefix (the same prefix passed to `impg syng -o`).
        #[clap(short = 'a', long = "syng-prefix", value_parser)]
        syng_prefix: String,

        /// Output GFA path (default: stdout; use "-" for stdout)
        #[clap(short = 'o', long, value_parser, default_value = "-")]
        output: String,

        /// GFA spec version: `1.0` (P lines, default) or `1.1` (W lines, PanSN-parsed).
        #[clap(long, value_parser, default_value = "1.0")]
        gfa_version: String,

        /// Sequence input for gap filling. Without these, gaps are filled with `N`s.
        #[clap(flatten)]
        sequence: SequenceOpts,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },

    /// Repair or create sampled positional sidecars for an existing syng index
    SyngRepair {
        /// Syng index prefix or .1khash/.1gbwt/.spos/.pstep/.names/.meta path
        #[clap(short = 'a', long, value_parser)]
        index: String,

        /// Resample and overwrite existing position sidecars.
        #[clap(long, action)]
        force: bool,

        /// Sample every 2^N syncmer steps per path for positional sidecars.
        #[arg(help_heading = "Position sampling")]
        #[clap(long, value_parser, default_value_t = impg::syng::DEFAULT_POSITION_SAMPLE_SHIFT)]
        position_sample_shift: u32,

        /// Compatibility seed stored in sidecars; regular-grid sampling ignores it.
        #[arg(help_heading = "Position sampling")]
        #[clap(long, value_parser, default_value_t = impg::syng::DEFAULT_POSITION_SAMPLE_SEED)]
        position_sample_seed: u64,

        /// Rebuild sidecars serially instead of walking paths in parallel.
        #[arg(help_heading = "Position sampling")]
        #[clap(long, action)]
        serial_position_sampling: bool,

        /// Log positional repair progress after this many completed paths (0 disables).
        #[arg(help_heading = "Position sampling")]
        #[clap(long, value_parser, default_value_t = 1)]
        position_progress_interval: usize,

        // --- General ---
        #[clap(flatten)]
        common: CommonOpts,
    },
}

fn main() {
    if let Err(e) = run() {
        error!("{}", e);
        std::process::exit(1);
    }
}

fn run() -> io::Result<()> {
    let args = Args::parse();

    match args {
        Args::Index {
            common,
            alignment,
            sequence,
        } => {
            initialize_threads_and_log(&common);
            let alignment_files = resolve_alignment_files(&alignment)?;
            let sequence_files = sequence.resolve_sequence_files()?;
            let _ = initialize_index(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files,
            )?;

            info!("Index created successfully");
        }
        Args::Lace {
            common,
            sequence,
            files,
            file_list,
            format,
            output,
            compress,
            fill_gaps,
            skip_validation,
            temp_dir,
            reference,
        } => {
            initialize_threads_and_log(&common);
            let temp_dir = resolve_temp_dir(temp_dir)?;
            setup_temp_dir(&temp_dir)?;

            // Check that at least one input is provided
            if files.is_none() && file_list.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --files or --file-list must be provided",
                ));
            }

            // Validate gap filling mode
            if fill_gaps > 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "fill_gaps must be 0, 1, or 2",
                ));
            }

            // Validate format
            let valid_formats = ["gfa", "vcf", "auto"];
            if !valid_formats.contains(&format.as_str()) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Invalid format '{}'. Must be one of: {}",
                        format,
                        valid_formats.join(", ")
                    ),
                ));
            }

            // Validate compression format
            let valid_compress = ["none", "gzip", "bgzip", "zstd", "auto"];
            if !valid_compress.contains(&compress.as_str()) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Invalid compression format '{}'. Must be one of: {}",
                        compress,
                        valid_compress.join(", ")
                    ),
                ));
            }

            // Determine the actual format (auto-detect if needed)
            let actual_format = determine_file_format(&format, &files, &file_list)?;

            if actual_format == "vcf" {
                // Build reference sequence index if provided
                let reference_index = if let Some(ref_file) = reference {
                    Some(UnifiedSequenceIndex::from_files(&[ref_file])?)
                } else {
                    None
                };

                // VCF lacing mode
                lace::run_vcf_lace(
                    files,
                    file_list,
                    &output,
                    &compress,
                    common.threads.get(),
                    common.verbose,
                    reference_index.as_ref(),
                )?;
            } else {
                // GFA lacing mode (existing functionality)
                // Build sequence index for sequence fetching (always build if sequence files provided)
                let sequence_index = sequence.build_sequence_index()?;

                lace::run_gfa_lace(
                    files,
                    file_list,
                    &output,
                    &compress,
                    fill_gaps,
                    skip_validation,
                    Some(temp_dir),
                    sequence_index.as_ref(),
                    common.verbose,
                )?;
            }
        }
        Args::Partition {
            common,
            alignment,
            syng_padding,
            syng_min_chain_anchors,
            syng_min_chain_fraction,
            no_rehome_singletons,
            window_size,
            output_format,
            output_folder,
            gfa_maf_fasta,
            engine_cli,
            merge_distance,
            min_result_identity,
            transitive_opts,
            starting_sequences_file,
            selection_mode,
            min_missing_size,
            min_boundary_distance,
            separate_files,
            approximate,
        } => {
            initialize_threads_and_log(&common);

            validate_selection_mode(&selection_mode)?;
            validate_output_format(&output_format, &["bed", "gfa", "maf", "fasta"])?;

            // Validate --approximate mode compatibility
            if approximate {
                if output_format != "bed" {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "--approximate mode is only compatible with 'bed' output format for partition, not '{}'",
                            output_format
                        ),
                    ));
                }
                validate_approximate_mode_min_length(transitive_opts.min_transitive_len, true)?;
                // Partition always uses transitive queries
            }

            // Parse engine spec early (engine + optional partition size)
            let (parsed_engine, parsed_partition_size) = engine_cli.parse_engine()?;

            // Validate partitioned mode + --separate-files are mutually exclusive
            if parsed_partition_size.is_some() && separate_files {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Partitioned --gfa-engine (e.g. 'pggb:10000') and --separate-files are mutually exclusive",
                ));
            }

            // For size validation, flat POA on "gfa" needs the same limit as "gfa-poa"
            let size_check_format = if output_format == "gfa" && parsed_engine == GfaEngine::Poa {
                "gfa-poa"
            } else {
                &output_format
            };
            // Skip region size validation when partitioned mode is active
            if parsed_partition_size.is_none() {
                validate_region_size(
                    0,
                    window_size as i32,
                    size_check_format,
                    merge_distance,
                    gfa_maf_fasta.force_large_region,
                )?;
            }

            // Validate single-file output compatibility
            // When partitioned mode is active, single-file GFA output is allowed (laced together)
            if !separate_files && output_format != "bed" && parsed_partition_size.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Single-file output is only supported for BED format. Use --separate-files for {} format.",
                        output_format.to_uppercase()
                    ),
                ));
            }

            // Allow `-a <prefix>` or `-a <prefix>.1khash` to route to syng.
            let effective_syng = resolve_syng_prefix(&alignment);

            // ─── Syng-based partition path ──────────────────────────────────
            if let Some(ref syng_prefix) = effective_syng {
                if approximate {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "--approximate mode is not supported with syng index input",
                    ));
                }

                info!("Loading syng index from prefix: {}", syng_prefix);
                let syng_index = impg::syng::SyngIndex::load(syng_prefix, impg::syng::SyncmerParams::default())?;
                let seq_index = syng_index.build_seq_index();

                let reverse_complement = gfa_maf_fasta.reverse_complement;
                let needs_sequences = matches!(output_format.as_str(), "gfa" | "maf" | "fasta");
                let sequence_index = if needs_sequences {
                    let si = gfa_maf_fasta.sequence.build_sequence_index()?;
                    if si.is_none() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("Sequence files are required for '{}' output with syng index input. Use --sequence-files or --sequence-list", output_format),
                        ));
                    }
                    si
                } else {
                    None
                };

                let scoring_params = if matches!(output_format.as_str(), "gfa" | "maf") {
                    Some(parse_poa_scoring_string(&engine_cli.poa_scoring)?)
                } else {
                    None
                };

                let engine_config = engine_cli.build(common.threads.get())?;
                let wrapper = {
                    let base = impg::SyngImpgWrapper::new(syng_index, seq_index, syng_padding);
                    if syng_min_chain_anchors > 0 {
                        base.with_chain_filter(syng_min_chain_anchors, syng_min_chain_fraction)
                    } else {
                        base
                    }
                };

                partition::partition_alignments(
                    &wrapper,
                    window_size,
                    starting_sequences_file.as_deref(),
                    &selection_mode,
                    merge_distance,
                    min_result_identity,
                    min_missing_size,
                    min_boundary_distance,
                    transitive_opts.transitive_dfs,
                    transitive_opts.max_depth,
                    transitive_opts.effective_min_transitive_len(),
                    transitive_opts.min_distance_between_ranges,
                    &output_format,
                    output_folder.as_deref(),
                    sequence_index.as_ref(),
                    scoring_params,
                    reverse_complement,
                    common.verbose > 1,
                    separate_files,
                    false, // approximate always false for syng
                    &engine_config,
                    !no_rehome_singletons,
                )?;
            } else {
            // ─── Normal (alignment-based) partition path ────────────────────

            // Resolve alignment files once (supports process substitution inputs)
            let alignment_files = resolve_alignment_files(&alignment)?;

            // Extract fields and resolve sequence files before moving gfa_maf_fasta
            let reverse_complement = gfa_maf_fasta.reverse_complement;
            let sequence_files_for_impg = gfa_maf_fasta.sequence.resolve_sequence_files()?;

            // Setup POA/sequence resources
            let (sequence_index, scoring_params) = gfa_maf_fasta.setup_output_resources(
                &output_format,
                false,
                alignment_files.as_slice(),
                approximate,
                &engine_cli.poa_scoring,
            )?;

            // Build engine config (resolves temp_dir and sparsify internally)
            let engine_config = engine_cli.build(common.threads.get())?;

            // Initialize impg after validation
            let impg = initialize_index(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files_for_impg,
            )?;

            partition::partition_alignments(
                &impg,
                window_size,
                starting_sequences_file.as_deref(),
                &selection_mode,
                merge_distance,
                min_result_identity,
                min_missing_size,
                min_boundary_distance,
                transitive_opts.transitive_dfs,
                transitive_opts.max_depth,
                transitive_opts.effective_min_transitive_len(),
                transitive_opts.min_distance_between_ranges,
                &output_format,
                output_folder.as_deref(),
                sequence_index.as_ref(),
                scoring_params,
                reverse_complement,
                common.verbose > 1,
                separate_files,
                approximate,
                &engine_config,
                !no_rehome_singletons,
            )?;
            }
        }
        Args::Query {
            common,
            alignment,
            syng_padding,
            syng_extension,
            syng_extend_budget,
            syng_min_chain_anchors,
            syng_min_chain_fraction,
            syng_raw,
            query,
            output_format,
            output_prefix,
            gfa_maf_fasta,
            engine_cli,
        } => {
            initialize_threads_and_log(&common);

            // Migration errors for removed format strings
            if output_format == "gfa-poa" {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Format 'gfa-poa' has been removed. Use '-o gfa --gfa-engine poa' instead.",
                ));
            }
            if output_format == "gfa-seqwish" {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Format 'gfa-seqwish' has been removed. Use '-o gfa --gfa-engine seqwish' instead.",
                ));
            }

            validate_output_format(
                &output_format,
                &[
                    "auto",
                    "bed",
                    "bedpe",
                    "paf",
                    "gfa",
                    "maf",
                    "fasta",
                    "fasta+paf",
                    "fasta-aln",
                    "gbwt",
                ],
            )?;

            // Check that either --target-range or --target-bed is provided (cheap validation)
            if query.target_range.is_none() && query.target_bed.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --target-range or --target-bed must be provided",
                ));
            }

            // Allow `-a <prefix>` or `-a <prefix>.1khash` to auto-route
            // to the syng backend.
            let effective_syng = resolve_syng_prefix(&alignment);

            // ─── Syng query path ──────────────────────────────────────────
            if let Some(ref syng_prefix) = effective_syng {
                // Validate that alignment files are NOT also provided
                // (conflicts_with_all handles this at clap level, but be explicit)

                // Only bed, bedpe, gfa, fasta, and gbwt output are supported for syng queries
                let resolved_format = if output_format == "auto" { "bed" } else { output_format.as_str() };
                if !matches!(resolved_format, "bed" | "bedpe" | "gfa" | "fasta" | "gbwt") {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "syng index queries currently support 'bed', 'bedpe', 'gfa', 'fasta', and 'gbwt' output formats, not '{}'",
                            resolved_format
                        ),
                    ));
                }

                info!("Loading syng index from prefix: {}", syng_prefix);
                let syng_index = impg::syng::SyngIndex::load(syng_prefix, impg::syng::SyncmerParams::default())?;
                let seq_index_built = syng_index.build_seq_index();
                // Wrap once; wrapper owns syng_index and seq_index for the rest of the path.
                let wrapper = impg::SyngImpgWrapper::new(syng_index, seq_index_built, syng_padding);

                // Parse target ranges
                let target_ranges: Vec<(String, (i32, i32), String)> = if let Some(ref target_range_str) = query.target_range {
                    if target_range_str.contains(':') {
                        vec![partition::parse_target_range(target_range_str)?]
                    } else {
                        // Whole sequence
                        let seq_name = target_range_str.as_str();
                        let seq_len = wrapper
                            .seq_index()
                            .get_id(seq_name)
                            .and_then(|id| wrapper.seq_index().get_len_from_id(id))
                            .ok_or_else(|| {
                                io::Error::new(
                                    io::ErrorKind::NotFound,
                                    format!("Sequence '{}' not found in syng index", seq_name),
                                )
                            })? as i32;
                        let name = format!("{}:{}-{}", seq_name, 0, seq_len);
                        vec![(seq_name.to_string(), (0, seq_len), name)]
                    }
                } else if let Some(ref target_bed) = query.target_bed {
                    partition::parse_bed_file(target_bed)?
                } else {
                    unreachable!();
                };

                // Boundary-realignment default path: needs a UnifiedSequenceIndex to
                // fetch small windows around each fuzzy edge. --syng-raw opts out and
                // uses the current syncmer-resolution pass-through. All output
                // formats (bed, fasta, gbwt, gfa) route through realignment by
                // default — fragmented raw intervals would otherwise feed straight
                // into the GFA partitioning pipeline and produce a fragmented graph.
                let use_boundary_realign = !syng_raw;
                let syng_max_depth = if query.transitive {
                    query.transitive_opts.max_depth.max(1)
                } else {
                    1
                };
                // Flows straight through to `cluster_by_signature`
                // on the hot syng path. There it's interpreted as the
                // maximum signature-space gap tolerated within one
                // homology block: within-block jitter is tens of bp
                // (local indels), paralog copies and long insertions
                // sit at structurally different signatures (kb-scale
                // on yeast). A typical `-d` up to a few kb cleanly
                // Setup output resources for GFA/FASTA/GBWT (need sequence files).
                // Boundary realignment also needs them for edge-window fetches.
                let sequence_files_supplied =
                    !gfa_maf_fasta.sequence.resolve_sequence_files()?.is_empty();
                let needs_sequences =
                    matches!(resolved_format, "gfa" | "fasta" | "gbwt")
                        || use_boundary_realign
                        || sequence_files_supplied;
                let sequence_index = if needs_sequences {
                    let si = gfa_maf_fasta.sequence.build_sequence_index()?;
                    if si.is_none() {
                        let reason = if use_boundary_realign && resolved_format == "bed" {
                            "boundary realignment (default syng behavior). Use --syng-raw for the syncmer-resolution pass-through, or".to_string()
                        } else {
                            format!("'{}' output with syng index input.", resolved_format)
                        };
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("Sequence files are required for {} Use --sequence-files or --sequence-list", reason),
                        ));
                    }
                    si
                } else {
                    None
                };

                // Validate that -O is provided for gbwt output
                if resolved_format == "gbwt" && output_prefix.is_none() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Output prefix (-O) is required for 'gbwt' output format",
                    ));
                }

                let scoring_params = if resolved_format == "gfa" {
                    Some(parse_poa_scoring_string(&engine_cli.poa_scoring)?)
                } else {
                    None
                };

                let query_raw_syng = |
                    target_name: &str,
                    range_start: i32,
                    range_end: i32,
                    padding: u64,
                    extension: u64,
                | -> io::Result<Vec<impg::syng::HomologousInterval>> {
                    let start = range_start.max(0) as u64;
                    let end = range_end.max(range_start).max(0) as u64;
                    wrapper.syng_index().query_region_ext(
                        target_name,
                        start,
                        end,
                        padding,
                        extension,
                    )
                };

                // Process each target range. Serial outer — the
                // per-tile query already uses the rayon pool
                // internally (anchor emission + chain projection);
                // an outer par_iter over tiles oversubscribes and
                // actually slowed things down in measurement.
                for (target_name, (range_start, range_end), name) in &target_ranges {
                    info!("Syng query: {} ({}:{}-{})", name, target_name, range_start, range_end);

                    match resolved_format {
                        "bed" | "bedpe" => {
                            let intervals = if use_boundary_realign {
                                impg::syng_transitive::query_transitive_ext(
                                    wrapper.syng_index(),
                                    target_name,
                                    *range_start as u64,
                                    *range_end as u64,
                                    syng_padding,
                                    syng_max_depth,
                                    syng_extension,
                                    syng_extend_budget,
                                    syng_min_chain_anchors,
                                    syng_min_chain_fraction,
                                    sequence_index.as_ref().unwrap(),
                                )?
                            } else {
                                query_raw_syng(
                                    target_name,
                                    *range_start,
                                    *range_end,
                                    syng_padding,
                                    syng_extension,
                                )?
                            };
                            // Unified path: convert syng output to AdjustedInterval,
                            // flow through the shared merge + emit pipeline so `-d`
                            // honors gap-tolerant 2D chaining.
                            let mut results = syng_intervals_to_adjusted(
                                &intervals,
                                target_name,
                                *range_start,
                                *range_end,
                                wrapper.seq_index(),
                            );
                            let ext = if resolved_format == "bedpe" { "bedpe" } else { "bed" };
                            let mut out = find_output_stream(&output_prefix, ext)?;
                            if resolved_format == "bedpe" {
                                output_results_bedpe(
                                    &wrapper,
                                    &mut results,
                                    &mut out,
                                    &name,
                                    query.effective_merge_distance(),
                                    query.original_sequence_coordinates,
                                )?;
                            } else {
                                output_results_bed(
                                    &wrapper,
                                    &mut results,
                                    &mut out,
                                    &name,
                                    query.effective_merge_distance(),
                                    query.merge_strands_for_output("bed"),
                                    query.original_sequence_coordinates,
                                )?;
                            }
                        }
                        "gfa" => {
                            let engine_opts = engine_cli.build(common.threads.get())?;

                            if let Some(partition_size) = engine_opts.partition_size {
                                // ─── Sub-windowed path: syng-at-the-outer-level ───
                                //
                                // Split the query range into `partition_size`-bp
                                // sub-windows, call `query_region` per sub-window
                                // (which returns tight, small-scale intervals),
                                // collect per-window query intervals, then run the
                                // partitioned GFA pipeline which does a fresh local
                                // alignment + graph induction per partition and
                                // laces them together with a single final
                                // gfaffix. Structurally mirrors the alignment
                                // path's `output_results_gfa_partitioned`.
                                let ps = partition_size as i32;
                                let mut partitions: Vec<Vec<Interval<u32>>> = Vec::new();
                                let mut window_start = *range_start;
                                let mut window_idx = 0usize;
                                while window_start < *range_end {
                                    let window_end = (window_start + ps).min(*range_end);
                                    let intervals = if use_boundary_realign {
                                        impg::syng_transitive::query_transitive_ext(
                                            wrapper.syng_index(),
                                            target_name,
                                            window_start as u64,
                                            window_end as u64,
                                            syng_padding,
                                            syng_max_depth,
                                            syng_extension,
                                            syng_extend_budget,
                                            syng_min_chain_anchors,
                                            syng_min_chain_fraction,
                                            sequence_index.as_ref().unwrap(),
                                        )?
                                    } else {
                                        query_raw_syng(target_name, window_start, window_end, syng_padding, syng_extension)?
                                    };
                                    let window_intervals: Vec<Interval<u32>> = intervals
                                        .iter()
                                        .filter_map(|iv| {
                                            let id = wrapper.seq_index().get_id(&iv.genome)?;
                                            Some(Interval::new(iv.start as i32, iv.end as i32, id))
                                        })
                                        .collect();
                                    info!(
                                        "  [syng sub-window {}] {}:{}-{} → {} intervals",
                                        window_idx, target_name, window_start, window_end,
                                        window_intervals.len()
                                    );
                                    if !window_intervals.is_empty() {
                                        partitions.push(window_intervals);
                                    }
                                    window_start = window_end;
                                    window_idx += 1;
                                }

                                if partitions.is_empty() {
                                    let mut out = find_output_stream(&output_prefix, "gfa")?;
                                    writeln!(out, "H\tVN:Z:1.0")?;
                                    continue;
                                }

                                let gfa_output = impg::partitioned_gfa_pipeline(
                                    &partitions,
                                    &wrapper,
                                    sequence_index.as_ref().unwrap(),
                                    scoring_params,
                                    &engine_opts,
                                )?;
                                let mut out = find_output_stream(&output_prefix, "gfa")?;
                                writeln!(out, "{gfa_output}")?;
                            } else {
                                // ─── Flat path: one query_region, one engine run ───
                                let intervals = if use_boundary_realign {
                                    impg::syng_transitive::query_transitive_ext(
                                        wrapper.syng_index(),
                                        target_name,
                                        *range_start as u64,
                                        *range_end as u64,
                                        syng_padding,
                                        syng_max_depth,
                                        syng_extension,
                                        syng_extend_budget,
                                        syng_min_chain_anchors,
                                        syng_min_chain_fraction,
                                        sequence_index.as_ref().unwrap(),
                                    )?
                                } else {
                                    query_raw_syng(target_name, *range_start, *range_end, syng_padding, syng_extension)?
                                };
                                let query_intervals: Vec<coitrees::Interval<u32>> = intervals
                                    .iter()
                                    .filter_map(|iv| {
                                        let id = wrapper.seq_index().get_id(&iv.genome)?;
                                        Some(Interval::new(iv.start as i32, iv.end as i32, id))
                                    })
                                    .collect();

                                if query_intervals.is_empty() {
                                    let mut out = find_output_stream(&output_prefix, "gfa")?;
                                    writeln!(out, "H\tVN:Z:1.0")?;
                                    continue;
                                }

                                let gfa_output = impg::dispatch_gfa_engine(
                                    &wrapper,
                                    &query_intervals,
                                    sequence_index.as_ref().unwrap(),
                                    scoring_params,
                                    &engine_opts,
                                )?;
                                let mut out = find_output_stream(&output_prefix, "gfa")?;
                                writeln!(out, "{gfa_output}")?;
                            }
                        }
                        "fasta" => {
                            let seq_idx = sequence_index.as_ref().unwrap();
                            let intervals = if use_boundary_realign {
                                impg::syng_transitive::query_transitive_ext(
                                    wrapper.syng_index(),
                                    target_name,
                                    *range_start as u64,
                                    *range_end as u64,
                                    syng_padding,
                                    syng_max_depth,
                                    syng_extension,
                                    syng_extend_budget,
                                    syng_min_chain_anchors,
                                    syng_min_chain_fraction,
                                    seq_idx,
                                )?
                            } else {
                                query_raw_syng(target_name, *range_start, *range_end, syng_padding, syng_extension)?
                            };
                            let mut out = find_output_stream(&output_prefix, "fa")?;
                            for iv in &intervals {
                                let sequence = seq_idx.fetch_sequence(&iv.genome, iv.start as i32, iv.end as i32)?;
                                writeln!(out, ">{}:{}-{}({})", iv.genome, iv.start, iv.end, iv.strand)?;
                                out.write_all(&sequence)?;
                                writeln!(out)?;
                            }
                        }
                        "gbwt" => {
                            let seq_idx = sequence_index.as_ref().unwrap();
                            let intervals = if use_boundary_realign {
                                impg::syng_transitive::query_transitive_ext(
                                    wrapper.syng_index(),
                                    target_name,
                                    *range_start as u64,
                                    *range_end as u64,
                                    syng_padding,
                                    syng_max_depth,
                                    syng_extension,
                                    syng_extend_budget,
                                    syng_min_chain_anchors,
                                    syng_min_chain_fraction,
                                    seq_idx,
                                )?
                            } else {
                                query_raw_syng(target_name, *range_start, *range_end, syng_padding, syng_extension)?
                            };
                            let gbwt_prefix = output_prefix.as_ref().unwrap();

                            // Fetch sequences for all intervals
                            let fetched: Vec<(String, Vec<u8>)> = intervals
                                .iter()
                                .map(|iv| {
                                    let sequence = seq_idx.fetch_sequence(&iv.genome, iv.start as i32, iv.end as i32)?;
                                    let seq_name = format!("{}:{}-{}({})", iv.genome, iv.start, iv.end, iv.strand);
                                    Ok((seq_name, sequence))
                                })
                                .collect::<io::Result<Vec<_>>>()?;

                            let seq_refs: Vec<(String, &[u8])> = fetched
                                .iter()
                                .map(|(name, seq)| (name.clone(), seq.as_slice()))
                                .collect();

                            wrapper.syng_index().build_region_gbwt(&seq_refs, gbwt_prefix)?;
                            info!("Wrote region GBWT: {}.1gbwt + {}.1khash", gbwt_prefix, gbwt_prefix);
                        }
                        _ => unreachable!(),
                    }
                }
                // Skip the rest of the normal query path
            } else {
            // ─── Normal (alignment-based) query path ──────────────────────

            // Resolve alignment files once (handles process substitution inputs)
            let alignment_files = resolve_alignment_files(&alignment)?;

            // Early validation for approximate mode (before expensive operations)
            if query.approximate {
                // Check that all alignment files are tracepoint-based (.1aln or .tpa)
                let non_tracepoint_files: Vec<&String> = alignment_files
                    .iter()
                    .filter(|f| !f.ends_with(".1aln") && !f.ends_with(".tpa"))
                    .collect();

                if !non_tracepoint_files.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "--approximate mode only works with tracepoint alignment files (.1aln or .tpa), found {} incompatible files.",
                            non_tracepoint_files.len()
                        ),
                    ));
                }

                validate_approximate_mode_min_length(
                    query.transitive_opts.min_transitive_len,
                    query.transitive,
                )?;
            }

            // Extract sequence files before consuming gfa_maf_fasta
            let sequence_files_for_impg = gfa_maf_fasta.sequence.resolve_sequence_files()?;

            // Load subset filter if provided
            let subset_filter = load_subset_filter_if_provided(&query.subset_sequence_list)?;

            // Initialize impg after validation but before target range validation (which needs seq_index)
            let impg = initialize_index(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files_for_impg,
            )?;

            // Parse engine spec early (engine + optional partition size)
            let (parsed_engine, parsed_partition_size) = engine_cli.parse_engine()?;

            // For size validation, flat POA on "gfa" needs the same limit as the old "gfa-poa"
            let size_check_format = if output_format == "gfa" && parsed_engine == GfaEngine::Poa {
                "gfa-poa"
            } else {
                &output_format
            };

            // Parse and validate all target ranges, tracking which parameter was used
            let (target_ranges, from_range_param) = if let Some(target_range_str) =
                &query.target_range
            {
                let (target_name, target_range, name) = if target_range_str.contains(':') {
                    partition::parse_target_range(target_range_str)?
                } else {
                    // No interval specified: use the whole sequence [0, len)
                    let seq_name = target_range_str;
                    let seq_id = impg.seq_index().get_id(seq_name).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::NotFound,
                            format!("Sequence '{seq_name}' not found in index"),
                        )
                    })?;
                    let seq_len = impg.seq_index().get_len_from_id(seq_id).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Could not get length for sequence '{seq_name}'"),
                        )
                    })? as i32;
                    let name = format!("{}:{}-{}", seq_name, 0, seq_len);
                    (seq_name.to_string(), (0, seq_len), name)
                };
                // Validate sequence exists and range is within bounds
                validate_sequence_range(
                    &target_name,
                    target_range.0,
                    target_range.1,
                    impg.seq_index(),
                )?;
                validate_range_min_length(
                    target_range.0,
                    target_range.1,
                    &name,
                    query.transitive_opts.effective_min_transitive_len(),
                )?;
                // Skip region size validation when --partition-size is set (each sub-window is within limits)
                if parsed_partition_size.is_none() {
                    validate_region_size(
                        target_range.0,
                        target_range.1,
                        size_check_format,
                        query.effective_merge_distance(),
                        gfa_maf_fasta.force_large_region,
                    )?;
                }
                (vec![(target_name, target_range, name)], true)
            } else if let Some(target_bed) = &query.target_bed {
                let targets = partition::parse_bed_file(target_bed)?;
                // Validate all entries in the BED file
                for (seq_name, (start, end), name) in &targets {
                    validate_sequence_range(seq_name, *start, *end, impg.seq_index())?;
                    validate_range_min_length(
                        *start,
                        *end,
                        name,
                        query.transitive_opts.effective_min_transitive_len(),
                    )?;
                    // Skip region size validation when --partition-size is set
                    if parsed_partition_size.is_none() {
                        validate_region_size(
                            *start,
                            *end,
                            size_check_format,
                            query.effective_merge_distance(),
                            gfa_maf_fasta.force_large_region,
                        )?;
                    }
                }
                (targets, false)
            } else {
                unreachable!("Already validated that either target_range or target_bed is present");
            };

            // Resolve output format based on 'auto' and parameter used
            let resolved_output_format = if output_format == "auto" {
                if from_range_param {
                    "bed"
                } else {
                    "bedpe"
                }
            } else {
                output_format.as_str()
            };

            // Validate --approximate mode output format compatibility
            if query.approximate
                && resolved_output_format != "bed"
                && resolved_output_format != "bedpe"
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "--approximate mode is only compatible with 'bed' and 'bedpe' output formats, not '{}'",
                        resolved_output_format
                    ),
                ));
            }

            // Validate gbwt output requires -O prefix
            if resolved_output_format == "gbwt" && output_prefix.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Output prefix (-O) is required for 'gbwt' output format",
                ));
            }

            // Extract reverse_complement before moving gfa_maf_fasta
            let reverse_complement = gfa_maf_fasta.reverse_complement;

            // Setup POA/sequence resources
            let (sequence_index, scoring_params) = gfa_maf_fasta.setup_output_resources(
                resolved_output_format,
                query.original_sequence_coordinates,
                alignment_files.as_slice(),
                query.approximate,
                &engine_cli.poa_scoring,
            )?;

            // Process all target ranges in a unified loop
            info!("Querying target ranges");
            for (target_name, target_range, name) in target_ranges {
                let mut results = perform_query(
                    &impg,
                    &target_name,
                    target_range,
                    resolved_output_format == "paf" || resolved_output_format == "bedpe", // Store CIGAR for PAF/BEDPE output
                    query.min_result_identity,
                    query.min_output_length,
                    query.transitive,
                    query.transitive_opts.transitive_dfs,
                    &query.transitive_opts,
                    sequence_index.as_ref(),
                    query.approximate,
                    subset_filter.as_ref(),
                )?;

                // Output results based on the resolved format
                match resolved_output_format {
                    "bed" => {
                        // BED format - include the first element
                        output_results_bed(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "bed")?,
                            &name,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("bed"),
                            query.original_sequence_coordinates,
                        )?;
                    }
                    "bedpe" => {
                        // Skip the first element (the input range) for BEDPE output
                        results.remove(0);
                        output_results_bedpe(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "bed")?,
                            &name,
                            query.effective_merge_distance(),
                            query.original_sequence_coordinates,
                        )?;
                    }
                    "paf" => {
                        // Skip the first element (the input range) for PAF output
                        results.remove(0);
                        output_results_paf(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "paf")?,
                            &name,
                            query.effective_merge_distance(),
                            query.original_sequence_coordinates,
                            sequence_index.as_ref(),
                        )?;
                    }
                    "gfa" => {
                        let engine_opts = engine_cli.build(common.threads.get())?;
                        if let Some(ps) = engine_opts.partition_size {
                            // Partitioned mode: split query region into sub-windows
                            output_results_gfa_partitioned(
                                &impg,
                                &mut find_output_stream(&output_prefix, "gfa")?,
                                sequence_index.as_ref().unwrap(),
                                scoring_params,
                                &engine_opts,
                                &target_name,
                                target_range,
                                ps,
                                &query,
                                subset_filter.as_ref(),
                            )?;
                        } else {
                            output_results_gfa(
                                &impg,
                                &mut results,
                                &mut find_output_stream(&output_prefix, "gfa")?,
                                sequence_index.as_ref().unwrap(),
                                &name,
                                query.effective_merge_distance(),
                                query.merge_strands_for_output("gfa"),
                                scoring_params,
                                &engine_opts,
                            )?;
                        }
                    }
                    "maf" => {
                        output_results_maf(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "maf")?,
                            sequence_index.as_ref().unwrap(),
                            &name,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("maf"),
                            scoring_params.unwrap(),
                        )?;
                    }
                    "fasta" => {
                        output_results_fasta(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "fa")?,
                            sequence_index.as_ref().unwrap(),
                            &name,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("fasta"),
                            reverse_complement,
                        )?;
                    }
                    "fasta+paf" => {
                        output_results_fasta(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "fa")?,
                            sequence_index.as_ref().unwrap(),
                            &name,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("fasta"),
                            reverse_complement,
                        )?;
                        // Skip the first element (the input range) for PAF output
                        results.remove(0);
                        output_results_paf(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "paf")?,
                            &name,
                            query.effective_merge_distance(),
                            query.original_sequence_coordinates,
                            sequence_index.as_ref(),
                        )?;
                    }
                    "fasta-aln" => {
                        output_results_fasta_aln(
                            &impg,
                            &mut results,
                            sequence_index.as_ref().unwrap(),
                            name.clone(),
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("fasta-aln"),
                            scoring_params.unwrap(),
                        )?;
                    }
                    "gbwt" => {
                        let seq_idx = sequence_index.as_ref().unwrap();
                        let gbwt_prefix = output_prefix.as_ref().ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidInput,
                                "Output prefix (-O) is required for 'gbwt' output format",
                            )
                        })?;

                        // Merge intervals before fetching
                        merge_query_adjusted_intervals(
                            &mut results,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("gbwt"),
                        );

                        // Fetch sequences for all result intervals
                        let fetched: Vec<(String, Vec<u8>)> = results
                            .iter()
                            .map(|(qi, _, _)| {
                                let qname = impg.seq_index().get_name(qi.metadata).unwrap();
                                let (start, end, strand) = if qi.first <= qi.last {
                                    (qi.first, qi.last, '+')
                                } else {
                                    (qi.last, qi.first, '-')
                                };
                                let sequence = seq_idx.fetch_sequence(qname, start, end)?;
                                let seq_name = format!("{}:{}-{}({})", qname, start, end, strand);
                                Ok((seq_name, sequence))
                            })
                            .collect::<io::Result<Vec<_>>>()?;

                        let seq_refs: Vec<(String, &[u8])> = fetched
                            .iter()
                            .map(|(name, seq)| (name.clone(), seq.as_slice()))
                            .collect();

                        // Build region GBWT using default syncmer params
                        let region_index = impg::syng::SyngIndex::new(impg::syng::SyncmerParams::default());
                        region_index.build_region_gbwt(&seq_refs, gbwt_prefix)?;
                        info!("Wrote region GBWT: {}.1gbwt + {}.1khash", gbwt_prefix, gbwt_prefix);
                    }
                    _ => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("Invalid output format: {resolved_output_format}"),
                        ));
                    }
                }
            }
            } // end of else (normal query path)
        }
        Args::Refine {
            alignment,
            refine,
            common,
        } => {
            initialize_threads_and_log(&common);
            refine.validate()?;

            let alignment_files = resolve_alignment_files(&alignment)?;
            let sequence_files = refine.sequence.resolve_sequence_files()?;

            // Check if we have .1aln files and validate sequence files are provided
            let has_onealn_files = alignment_files.iter().any(|f| f.ends_with(".1aln"));
            if has_onealn_files && sequence_files.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Sequence files (FASTA or AGC) are required for the 'refine' command with .1aln alignment files to convert tracepoints to CIGAR strings. Use --sequence-files or --sequence-list".to_string(),
                ));
            }

            let impg = initialize_index(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files,
            )?;
            let subset_filter = load_subset_filter_if_provided(&refine.query.subset_sequence_list)?;

            let target_ranges = if let Some(target_range_str) = &refine.query.target_range {
                let (target_name, target_range, name) =
                    partition::parse_target_range(target_range_str)?;
                validate_sequence_range(
                    &target_name,
                    target_range.0,
                    target_range.1,
                    impg.seq_index(),
                )?;
                validate_range_min_length(
                    target_range.0,
                    target_range.1,
                    &name,
                    refine.query.transitive_opts.effective_min_transitive_len(),
                )?;
                vec![(target_name, target_range, name)]
            } else if let Some(target_bed) = &refine.query.target_bed {
                let targets = partition::parse_bed_file(target_bed)?;
                let mut validated = Vec::with_capacity(targets.len());
                for (seq_name, (start, end), name) in targets {
                    validate_sequence_range(&seq_name, start, end, impg.seq_index())?;
                    validate_range_min_length(
                        start,
                        end,
                        &name,
                        refine.query.transitive_opts.effective_min_transitive_len(),
                    )?;
                    validated.push((seq_name, (start, end), name));
                }
                validated
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --target-range or --target-bed must be provided",
                ));
            };

            // Parse blacklist BED file if provided
            let blacklist = if let Some(ref blacklist_path) = refine.blacklist_bed {
                Some(refine::parse_blacklist_bed(blacklist_path)?)
            } else {
                None
            };

            let support_level = refine
                .pansn_mode
                .unwrap_or(sweepga::pansn::PanSnLevel::Sequence);

            let config = refine::RefineConfig {
                span_bp: refine.span_bp,
                max_extension: refine.max_extension,
                extension_step: refine.extension_step,
                support_level,
                merge_distance: refine.query.effective_merge_distance(),
                min_identity: refine.query.min_result_identity,
                use_transitive_bfs: refine.query.transitive,
                use_transitive_dfs: refine.query.transitive_opts.transitive_dfs,
                max_transitive_depth: refine.query.transitive_opts.max_depth,
                min_transitive_len: refine.query.transitive_opts.effective_min_transitive_len(),
                min_distance_between_ranges: refine
                    .query
                    .transitive_opts
                    .min_distance_between_ranges,
                subset_filter: subset_filter.as_ref(),
                blacklist: blacklist.as_ref(),
                approximate_mode: refine.query.approximate,
            };

            let mut records = refine::run_refine(&impg, &target_ranges, config)?;
            info!(
                "Refining {} targets with max_extension={} (level: {:?})",
                target_ranges.len(),
                refine.max_extension,
                support_level
            );
            let mut writer = BufWriter::new(io::stdout());
            let mut support_writer = if let Some(path) = &refine.support_output {
                let support_path = Path::new(path);
                if let Some(parent) = support_path.parent() {
                    if !parent.as_os_str().is_empty() {
                        std::fs::create_dir_all(parent)?;
                    }
                }
                Some(BufWriter::new(File::create(support_path)?))
            } else {
                None
            };

            // Write header
            writeln!(
                writer,
                "#chrom\tstart\tend\tname\toriginal.support\tnew.support\tleft.extension.bp\tright.extension.bp"
            )?;

            for record in records.drain(..) {
                let original_range = format!(
                    "{}:{}-{}",
                    record.chrom, record.original_start, record.original_end
                );
                let mut name_field = record.label.clone();
                if name_field.trim().is_empty() || name_field == "." {
                    name_field = original_range.clone();
                }

                // Emit an informative BED-like row: chrom start end name new_support old_support left_extension right_extension.
                // Maximizing the sample count while minimizing the expansion helps avoid loci that start or end inside SVs.
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    record.chrom,
                    record.refined_start,
                    record.refined_end,
                    name_field,
                    record.original_support_count,
                    record.support_count,
                    record.applied_left_extension,
                    record.applied_right_extension
                )?;

                if let Some(ref mut support_out) = support_writer {
                    for entity in &record.support_entities {
                        writeln!(
                            support_out,
                            "{}\t{}\t{}\t{}",
                            entity.sequence, entity.start, entity.end, name_field
                        )?;
                    }
                }
            }
            writer.flush()?;
            if let Some(mut support_out) = support_writer {
                support_out.flush()?;
            }
        }
        Args::Similarity {
            common,
            alignment,
            query,
            gfa_maf_fasta,
            distances,
            all,
            delim,
            delim_pos,
            pca,
            pca_components,
            pca_measure,
            polarize_n_prev,
            polarize_guide_samples,
            engine_cli,
        } => {
            initialize_threads_and_log(&common);

            // Validate delim_pos
            if delim_pos < 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--delim-pos must be greater than 0",
                ));
            }

            if pca {
                // Validate components
                if pca_components == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Number of components must be greater than 0",
                    ));
                }

                // Validate pca_measure
                if !["jaccard", "cosine", "dice"].contains(&pca_measure.as_str()) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Invalid PCA similarity measure '{pca_measure}'. Must be one of: jaccard, cosine, dice"
                        ),
                    ));
                }
            }

            // Check that either --target-range or --target-bed is provided
            if query.target_range.is_none() && query.target_bed.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --target-range or --target-bed must be provided",
                ));
            }

            // Parse engine spec early
            let (parsed_engine, parsed_partition_size) = engine_cli.parse_engine()?;

            if parsed_partition_size.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Partitioned mode is not yet supported for the similarity command",
                ));
            }

            if parsed_engine == GfaEngine::Seqwish || parsed_engine == GfaEngine::Pggb {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "The 'seqwish' and 'pggb' engines are not yet implemented for the similarity command. Use 'poa'.",
                ));
            }
            engine_cli.validate_engine_params(parsed_engine)?;

            // Resolve alignment files once (supports process substitution inputs)
            let alignment_files = resolve_alignment_files(&alignment)?;

            // Extract fields and resolve sequence files before moving gfa_maf_fasta
            let force_large_region = gfa_maf_fasta.force_large_region;
            let sequence_files_for_impg = gfa_maf_fasta.sequence.resolve_sequence_files()?;

            // Setup POA/sequence resources (always required for similarity)
            let (sequence_index, scoring_params) = gfa_maf_fasta.setup_output_resources(
                "gfa",
                false,
                alignment_files.as_slice(),
                false,
                &engine_cli.poa_scoring,
            )?;
            let sequence_index = sequence_index.unwrap(); // Safe since "gfa" always requires sequence files
            let scoring_params = scoring_params.unwrap(); // Safe since "gfa" always requires POA scoring

            let subset_filter = load_subset_filter_if_provided(&query.subset_sequence_list)?;

            // Initialize impg after validation but before target range validation (which needs seq_index)
            let impg = initialize_index(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files_for_impg,
            )?;

            // Validate target_range and target_bed before any expensive operations
            let target_ranges = {
                let mut targets = Vec::new();

                if let Some(target_range_str) = &query.target_range {
                    let (target_name, target_range, name) = if target_range_str.contains(':') {
                        partition::parse_target_range(target_range_str)?
                    } else {
                        // No interval specified: use the whole sequence [0, len)
                        let seq_name = target_range_str;
                        let seq_id = impg.seq_index().get_id(seq_name).ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::NotFound,
                                format!("Sequence '{seq_name}' not found in index"),
                            )
                        })?;
                        let seq_len = impg.seq_index().get_len_from_id(seq_id).ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("Could not get length for sequence '{seq_name}'"),
                            )
                        })? as i32;
                        let name = format!("{}:{}-{}", seq_name, 0, seq_len);
                        (seq_name.to_string(), (0, seq_len), name)
                    };
                    validate_sequence_range(
                        &target_name,
                        target_range.0,
                        target_range.1,
                        impg.seq_index(),
                    )?;
                    validate_region_size(
                        target_range.0,
                        target_range.1,
                        "gfa",
                        query.effective_merge_distance(),
                        force_large_region,
                    )?;
                    targets.push((target_name, target_range, name));
                }

                if let Some(target_bed) = &query.target_bed {
                    let bed_targets = partition::parse_bed_file(target_bed)?;
                    for (target_name, target_range, name) in bed_targets {
                        validate_sequence_range(
                            &target_name,
                            target_range.0,
                            target_range.1,
                            impg.seq_index(),
                        )?;
                        validate_region_size(
                            target_range.0,
                            target_range.1,
                            "gfa",
                            query.effective_merge_distance(),
                            force_large_region,
                        )?;
                        targets.push((target_name, target_range, name));
                    }
                }

                targets
            };
            // Note: Already validated that either target_range or target_bed is present

            info!("Parsed {} target ranges from BED file", target_ranges.len());

            // Query all regions serially (already parallelized internally)
            let mut all_query_data = Vec::new();
            for (target_name, target_range, _name) in target_ranges.into_iter() {
                let mut results = perform_query(
                    &impg,
                    &target_name,
                    target_range,
                    false, // Don't need CIGAR for similarity
                    query.min_result_identity,
                    query.min_output_length,
                    query.transitive,
                    query.transitive_opts.transitive_dfs,
                    &query.transitive_opts,
                    Some(&sequence_index),
                    query.approximate,
                    subset_filter.as_ref(),
                )?;

                let region_label = format!("{}:{}-{}", target_name, target_range.0, target_range.1);

                // Merge intervals if needed
                merge_query_adjusted_intervals(
                    &mut results,
                    query.effective_merge_distance(),
                    query.merge_strands_for_output("similarity"),
                );

                // Extract query intervals
                let query_intervals: Vec<Interval<u32>> = results
                    .iter()
                    .map(|(query_interval, _, _)| *query_interval)
                    .collect();

                all_query_data.push((query_intervals, region_label));
            }

            // Process all regions in parallel
            similarity::compute_and_output_similarities(
                &impg,
                all_query_data,
                &sequence_index,
                scoring_params,
                distances,
                all,
                delim,
                delim_pos,
                pca,
                pca_components,
                &pca_measure,
                polarize_n_prev,
                polarize_guide_samples.as_deref(),
            )?;
        }
        Args::Stats {
            common,
            alignment,
            sequence,
            list_sequences,
        } => {
            initialize_threads_and_log(&common);
            let alignment_files = resolve_alignment_files(&alignment)?;
            let sequence_files = sequence.resolve_sequence_files()?;
            let impg = initialize_index(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files,
            )?;

            if list_sequences {
                print_list_sequences(&impg);
            } else {
                print_stats(&impg);
            }
        }
        Args::Graph {
            fasta_input,
            paf_file,
            output,
            engine_cli,
            common,
        } => {
            initialize_threads_and_log(&common);
            let (parsed_engine, parsed_partition_size) = engine_cli.parse_engine()?;
            engine_cli.validate_engine_params(parsed_engine)?;
            let temp_dir = resolve_temp_dir(engine_cli.aln.sw.tempdir.clone())?;
            setup_temp_dir(&temp_dir)?;

            let fasta_files = fasta_input.resolve_sequence_files()?;
            if fasta_files.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "No sequence files specified. Use --sequence-files or --sequence-list",
                ));
            }

            // Wfmash-density sparsification only makes sense with the wfmash backend
            if matches!(
                engine_cli.aln.sw.sparsify,
                SparsificationStrategy::WfmashDensity(_)
            ) && engine_cli.aln.sw.aligner != "wfmash"
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Wfmash density sparsification ({}) requires --aligner wfmash, but '{}' was specified",
                        engine_cli.aln.sw.sparsify, engine_cli.aln.sw.aligner
                    ),
                ));
            }
            if engine_cli.aln.sw.aligner == "wfmash" && engine_cli.aln.sw.frequency.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--fastga-frequency is only supported with --aligner fastga; wfmash uses its own frequency estimation",
                ));
            }

            let poa_scoring = engine_cli.parse_poa_scoring()?;

            let show_progress = common.verbose > 0;
            let num_threads = common.threads.get();

            if let Some(ps) = parsed_partition_size {
                // Partitioned mode: align → IMPG → partition → per-partition engine → lace → gfaffix
                let engine_opts = engine_cli.build(num_threads)?;
                let graph_config = build_graph_config(
                    &engine_cli,
                    &engine_cli.aln.sw.sparsify,
                    paf_file,
                    &temp_dir,
                    num_threads,
                    show_progress,
                )?;

                let scoring = Some(poa_scoring);
                graph::run_graph_build_partitioned(
                    fasta_files,
                    &output,
                    &graph_config,
                    &engine_opts,
                    scoring,
                    ps,
                )?;
            } else {
                let graph_config = build_graph_config(
                    &engine_cli,
                    &engine_cli.aln.sw.sparsify,
                    paf_file,
                    &temp_dir,
                    num_threads,
                    show_progress,
                )?;

                match parsed_engine {
                    GfaEngine::Poa => {
                        let scoring = poa_scoring;

                        if output == "-" {
                            let stdout = io::stdout();
                            let mut out = BufWriter::with_capacity(1024 * 1024, stdout.lock());
                            graph::run_graph_build_poa(
                                fasta_files,
                                &mut out,
                                scoring,
                                &graph_config,
                            )?;
                        } else {
                            let mut out =
                                BufWriter::with_capacity(1024 * 1024, File::create(&output)?);
                            graph::run_graph_build_poa(
                                fasta_files,
                                &mut out,
                                scoring,
                                &graph_config,
                            )?;
                        }
                    }
                    GfaEngine::Seqwish => {
                        graph::run_graph_build(fasta_files, &output, graph_config)?;
                    }
                    GfaEngine::SyngNative => {
                        // `impg graph` is the flat whole-FASTA entry point and
                        // doesn't produce syng partitions. Syng-native only has
                        // a meaning inside `partition --gfa-engine syng-native`.
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "--gfa-engine syng-native is only available under `partition` (requires a syng index); use `seqwish` or `pggb` here",
                        ));
                    }
                    GfaEngine::Pggb => {
                        let target_poa_lengths = engine_cli.smooth.parse_target_poa_lengths()?;
                        if output == "-" {
                            let stdout = io::stdout();
                            let mut out = BufWriter::with_capacity(1024 * 1024, stdout.lock());
                            graph::run_graph_build_pggb(
                                fasta_files,
                                &mut out,
                                &graph_config,
                                target_poa_lengths,
                                engine_cli.smooth.max_node_length,
                                engine_cli.smooth.poa_padding_fraction,
                            )?;
                        } else {
                            let mut out =
                                BufWriter::with_capacity(1024 * 1024, File::create(&output)?);
                            graph::run_graph_build_pggb(
                                fasta_files,
                                &mut out,
                                &graph_config,
                                target_poa_lengths,
                                engine_cli.smooth.max_node_length,
                                engine_cli.smooth.poa_padding_fraction,
                            )?;
                        }
                    }
                }
            }
        }
        Args::Align {
            fasta_input,
            output_dir,
            format,
            run_joblist,
            jobs,
            aln,
            common,
        } => {
            initialize_threads_and_log(&common);

            // Short-circuit: execute a pre-generated joblist and exit.
            if let Some(joblist_path) = run_joblist {
                let jobs = jobs.unwrap_or_else(|| common.threads.get());
                impg::commands::align::run_joblist(
                    std::path::Path::new(&joblist_path),
                    jobs,
                )?;
                return Ok(());
            }

            let temp_dir = resolve_temp_dir(aln.sw.tempdir.clone())?;
            setup_temp_dir(&temp_dir)?;

            let fasta_files = fasta_input.resolve_sequence_files()?;
            if fasta_files.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "No sequence files specified. Use --sequence-files or --sequence-list",
                ));
            }

            if aln.sw.aligner == "wfmash" && aln.sw.frequency.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--fastga-frequency is only supported with --aligner fastga; wfmash uses its own frequency estimation",
                ));
            }

            // Parse output format
            let output_format = match format.to_lowercase().as_str() {
                "paf" => impg::commands::align::AlignOutputFormat::Paf,
                "1aln" | "onealn" => impg::commands::align::AlignOutputFormat::OneAln,
                "joblist" | "jobs" => impg::commands::align::AlignOutputFormat::JobList,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Unknown output format: {}. Valid: paf, 1aln, joblist",
                            format
                        ),
                    ));
                }
            };

            // `aln.sw.min_identity` is a String (preset or "ani50-2" etc); for
            // the `impg align` config we want the numeric form — default 0.0
            // means "let sweepga decide".
            let min_identity_f64 = aln.sw.min_identity.parse::<f64>().unwrap_or(0.0);

            let config = impg::commands::align::AlignConfig {
                num_threads: common.threads.get(),
                sparsify: aln.sw.sparsify.clone(),
                mash_params: sweepga::knn_graph::MashParams {
                    kmer_size: aln.sw.mash_kmer_size,
                    sketch_size: aln.sw.mash_sketch_size,
                },
                frequency_multiplier: aln.sw.fastga_frequency_multiplier,
                frequency: aln.sw.frequency,
                min_aln_length: aln.sw.block_length.unwrap_or(0),
                output_format,
                show_progress: common.verbose > 0,
                aligner: aln.sw.aligner.clone(),
                temp_dir: Some(temp_dir),
                batch_bytes: aln.sw.batch_bytes.clone(),
                no_filter: aln.sw.no_filter,
                num_mappings: aln.sw.num_mappings.clone(),
                scaffold_jump: aln.sw.scaffold_jump,
                scaffold_mass: aln.sw.scaffold_mass,
                scaffold_filter: aln.sw.scaffold_filter.clone(),
                overlap: aln.sw.overlap,
                min_identity: min_identity_f64,
                scaffold_dist: aln.sw.scaffold_dist,
                min_map_length: aln.min_map_length,
            };

            impg::commands::align::run_align(fasta_files, &output_dir, config)?;
        }
        Args::Map {
            index,
            query,
            output_format,
            output,
            min_anchors,
            chain_budget,
            max_hits,
            common,
        } => {
            initialize_threads_and_log(&common);

            let syng_prefix = detect_syng_prefix(&index).unwrap_or(index);
            info!("Loading syng index from prefix: {}", syng_prefix);
            #[cfg(unix)]
            let saved_stdout = silence_stdout_for_process()?;
            #[cfg(not(unix))]
            silence_stdout_for_process()?;

            let syng_index =
                impg::syng::SyngIndex::load(&syng_prefix, impg::syng::SyncmerParams::default())?;
            let queries = read_query_sequences(&query)?;

            if let Some(path) = output {
                let mut out = BufWriter::new(File::create(path)?);
                emit_syng_map(
                    &mut out,
                    &syng_index,
                    queries,
                    &output_format,
                    min_anchors,
                    chain_budget,
                    max_hits,
                )?;
                out.flush()?;
                #[cfg(unix)]
                unsafe {
                    libc::close(saved_stdout);
                }
            } else {
                let mut out = Vec::new();
                emit_syng_map(
                    &mut out,
                    &syng_index,
                    queries,
                    &output_format,
                    min_anchors,
                    chain_budget,
                    max_hits,
                )?;
                #[cfg(unix)]
                {
                    write_all_to_fd(saved_stdout, &out)?;
                    unsafe {
                        libc::close(saved_stdout);
                    }
                }
                #[cfg(not(unix))]
                {
                    let mut stdout = io::stdout();
                    stdout.write_all(&out)?;
                    stdout.flush()?;
                }
            }
        }
        Args::Syng {
            agc,
            fasta,
            output,
            syncmer_k,
            smer_length,
            syncmer_w,
            syncmer_length,
            syncmer_seed,
            position_sample_shift,
            position_sample_seed,
            parallel_dictionary,
            common,
        } => {
            initialize_threads_and_log(&common);

            if agc.is_none() && fasta.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --agc or --fasta must be provided",
                ));
            }

            let params = resolve_syng_syncmer_params(
                syncmer_k,
                smer_length,
                syncmer_w,
                syncmer_length,
                syncmer_seed,
            )?;
            let paper_k = params.k + params.w;
            let paper_s = params.k;

            info!(
                "Building syng index with syncmer params: paper k={}, paper s={}, internal k={}, w={}, seed={}, threads={}",
                paper_k,
                paper_s,
                params.k,
                params.w,
                syncmer_seed,
                common.threads.get(),
            );
            let configure_position_sampling =
                |index: &mut impg::syng::SyngIndex| -> io::Result<()> {
                    index.enable_online_sampled_positions(
                        position_sample_shift,
                        position_sample_seed,
                    )?;
                    info!(
                        "Online regular sampled position sidecars enabled: sample_shift={} (every {} syncmer steps per path)",
                        position_sample_shift,
                        1u64 << position_sample_shift,
                    );
                    Ok(())
                };
            let total_start = Instant::now();

            let mut index = if let Some(agc_path) = agc {
                if parallel_dictionary {
                    let input_start = Instant::now();
                    info!(
                        "Reading AGC sequence catalog for parallel dictionary build: {} at +{}",
                        agc_path,
                        format_duration(total_start.elapsed())
                    );
                    let records = collect_syng_agc_records(&agc_path)?;
                    let total_bp: u64 = records.iter().map(|record| record.length as u64).sum();
                    info!(
                        "Found {} AGC sequences ({} bp) for dictionary prepass",
                        records.len(),
                        total_bp
                    );

                    let dictionary_start = Instant::now();
                    let agc_path_for_threads = agc_path.clone();
                    let packed_syncmers = records
                        .par_iter()
                        .map_init(
                            move || {
                                let config = ragc_core::DecompressorConfig { verbosity: 0 };
                                ragc_core::Decompressor::open(&agc_path_for_threads, config)
                                    .map_err(|e| {
                                        io::Error::new(
                                            io::ErrorKind::InvalidData,
                                            format!(
                                                "Failed to open AGC file '{}': {}",
                                                agc_path_for_threads, e
                                            ),
                                        )
                                    })
                            },
                            |decompressor, record| {
                                let decompressor = match decompressor {
                                    Ok(decompressor) => decompressor,
                                    Err(e) => {
                                        return Err(io::Error::new(e.kind(), e.to_string()));
                                    }
                                };
                                let seq = read_syng_agc_record(decompressor, record)?;
                                impg::syng_parallel::extract_packed_syncmers(params, &seq)
                            },
                        )
                        .try_reduce(Vec::new, |mut acc, mut part| {
                            acc.append(&mut part);
                            Ok(acc)
                        })?;
                    let raw_syncmers = packed_syncmers.len();
                    let dictionary =
                        impg::syng_parallel::sort_dedup_packed_syncmers(packed_syncmers);
                    info!(
                        "Built packed syncmer dictionary in {}: {} raw syncmers, {} unique syncmer nodes at +{}",
                        format_duration(dictionary_start.elapsed()),
                        raw_syncmers,
                        dictionary.len(),
                        format_duration(total_start.elapsed())
                    );

                    let preload_start = Instant::now();
                    let mut index = impg::syng::SyngIndex::new_with_packed_syncmer_dictionary(
                        params,
                        &dictionary,
                    )?;
                    configure_position_sampling(&mut index)?;
                    info!(
                        "Preloaded {} syncmer nodes into KmerHash in {}",
                        dictionary.len(),
                        format_duration(preload_start.elapsed())
                    );

                    let config = ragc_core::DecompressorConfig { verbosity: 0 };
                    let mut decompressor =
                        ragc_core::Decompressor::open(&agc_path, config).map_err(|e| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("Failed to open AGC file '{}': {}", agc_path, e),
                            )
                        })?;
                    let mut sequence_count = 0usize;
                    let mut indexed_bp = 0u64;
                    let mut total_syncmers = 0usize;
                    for (record_idx, record) in records.iter().enumerate() {
                        let contig_start = Instant::now();
                        let seq = read_syng_agc_record(&mut decompressor, record)?;
                        info!("  Processing {} ({} bp)", record.name, seq.len());
                        let stats = index
                            .add_sequence_with_existing_syncmers(record.name.clone(), seq)?;
                        sequence_count += 1;
                        indexed_bp += stats.sequence_len as u64;
                        total_syncmers += stats.syncmers;
                        info!(
                            "  Indexed contig {}/{} {}: {} bp, {} syncmers, indexed={}, elapsed {}, total {} sequences / {} bp / {} syncmers at +{}",
                            record_idx + 1,
                            records.len(),
                            record.name,
                            stats.sequence_len,
                            stats.syncmers,
                            stats.indexed,
                            format_duration(contig_start.elapsed()),
                            sequence_count,
                            indexed_bp,
                            total_syncmers,
                            format_duration(total_start.elapsed())
                        );
                    }

                    info!(
                        "Built syng paths from {} sequences ({} bp, {} syncmers) in {}",
                        sequence_count,
                        indexed_bp,
                        total_syncmers,
                        format_duration(input_start.elapsed())
                    );
                    index
                } else {
                // Stream sequences from AGC
                let input_start = Instant::now();
                info!(
                    "Reading sequences from AGC: {} at +{}",
                    agc_path,
                    format_duration(total_start.elapsed())
                );
                let config = ragc_core::DecompressorConfig { verbosity: 0 };
                let mut decompressor =
                    ragc_core::Decompressor::open(&agc_path, config).map_err(|e| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Failed to open AGC file '{}': {}", agc_path, e),
                        )
                    })?;

                let samples = decompressor.list_samples();
                info!("Found {} AGC samples", samples.len());
                let mut index = impg::syng::SyngIndex::new(params);
                configure_position_sampling(&mut index)?;
                let mut sequence_count = 0usize;
                let mut total_bp = 0u64;
                let mut total_syncmers = 0usize;
                for (sample_idx, sample) in samples.iter().enumerate() {
                    let sample_start = Instant::now();
                    let contigs = decompressor
                        .list_contigs_names_only(sample)
                        .unwrap_or_default();
                    info!(
                        "Including sample {}/{}: {} ({} contigs) at +{}",
                        sample_idx + 1,
                        samples.len(),
                        sample,
                        contigs.len(),
                        format_duration(total_start.elapsed())
                    );
                    let mut sample_sequences = 0usize;
                    let mut sample_bp = 0u64;
                    let mut sample_syncmers = 0usize;
                    for (contig_idx, contig) in contigs.iter().enumerate() {
                        let contig_start = Instant::now();
                        // NOTE: Do NOT use `decompressor.get_contig()` here — it has a
                        // bug in ragc-core where it skips the detail-reload after a
                        // `list_contigs_names_only()` call (it checks contig count, which
                        // is populated, instead of `are_details_loaded()`), and silently
                        // returns an empty Vec. Use get_contig_length + get_contig_range
                        // instead — both correctly reload details. Same pattern as
                        // src/agc_index.rs.
                        let length = decompressor
                            .get_contig_length(sample, contig)
                            .map_err(|e| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!(
                                        "Failed to get length for contig '{}@{}': {}",
                                        contig, sample, e
                                    ),
                                )
                            })?;
                        let contig_data = decompressor
                            .get_contig_range(sample, contig, 0, length)
                            .map_err(|e| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!(
                                        "Failed to decompress contig '{}@{}': {}",
                                        contig, sample, e
                                    ),
                                )
                            })?;

                        let seq: Vec<u8> = contig_data
                            .iter()
                            .map(|&b| match b {
                                0 => b'A',
                                1 => b'C',
                                2 => b'G',
                                3 => b'T',
                                _ => b'N',
                            })
                            .collect();

                        // If the contig name already follows PanSN convention
                        // (sample#haplotype#contig), use it as-is — it already
                        // embeds the sample, so appending @sample would create
                        // a redundant suffix that makes query -r lookups fail.
                        // Otherwise (raw contig name like `chr1`), prepend the
                        // sample via `contig@sample` so names stay unique across
                        // samples that share a contig name.
                        let name = if contig.contains('#') {
                            contig.clone()
                        } else {
                            format!("{}@{}", contig, sample)
                        };
                        let seq_len = seq.len();
                        info!("  Processing {} ({} bp)", name, seq_len);
                        let stats = index.add_sequence(name.clone(), seq);
                        sequence_count += 1;
                        total_bp += stats.sequence_len as u64;
                        total_syncmers += stats.syncmers;
                        sample_sequences += 1;
                        sample_bp += stats.sequence_len as u64;
                        sample_syncmers += stats.syncmers;
                        info!(
                            "  Indexed contig {}/{} {}: {} bp, {} syncmers, indexed={}, elapsed {}, total {} sequences / {} bp / {} syncmers at +{}",
                            contig_idx + 1,
                            contigs.len(),
                            name,
                            stats.sequence_len,
                            stats.syncmers,
                            stats.indexed,
                            format_duration(contig_start.elapsed()),
                            sequence_count,
                            total_bp,
                            total_syncmers,
                            format_duration(total_start.elapsed())
                        );
                    }
                    info!(
                        "Finished sample {}/{}: {} ({} sequences, {} bp, {} syncmers) in {} at +{}",
                        sample_idx + 1,
                        samples.len(),
                        sample,
                        sample_sequences,
                        sample_bp,
                        sample_syncmers,
                        format_duration(sample_start.elapsed()),
                        format_duration(total_start.elapsed())
                    );
                }

                info!(
                    "Built syng paths from {} sequences ({} bp, {} syncmers) in {}",
                    sequence_count,
                    total_bp,
                    total_syncmers,
                    format_duration(input_start.elapsed())
                );
                index
                }
            } else if let Some(fasta_path) = fasta {
                if parallel_dictionary {
                    let input_start = Instant::now();
                    info!(
                        "Reading sequences from FASTA for parallel dictionary build: {} at +{}",
                        fasta_path,
                        format_duration(total_start.elapsed())
                    );
                    let sequences = read_syng_fasta_sequences(&fasta_path)?;
                    let total_sequences = sequences.len();
                    let total_bp: u64 =
                        sequences.iter().map(|(_, seq)| seq.len() as u64).sum();
                    info!(
                        "Read {} FASTA sequences ({} bp) for dictionary prepass",
                        sequences.len(),
                        total_bp
                    );

                    let dictionary_start = Instant::now();
                    let dictionary =
                        impg::syng_parallel::build_packed_syncmer_dictionary(params, &sequences)?;
                    info!(
                        "Built packed syncmer dictionary in {}: {} unique syncmer nodes at +{}",
                        format_duration(dictionary_start.elapsed()),
                        dictionary.len(),
                        format_duration(total_start.elapsed())
                    );

                    let preload_start = Instant::now();
                    let mut index = impg::syng::SyngIndex::new_with_packed_syncmer_dictionary(
                        params,
                        &dictionary,
                    )?;
                    configure_position_sampling(&mut index)?;
                    info!(
                        "Preloaded {} syncmer nodes into KmerHash in {}",
                        dictionary.len(),
                        format_duration(preload_start.elapsed())
                    );

                    let mut sequence_count = 0usize;
                    let mut total_syncmers = 0usize;
                    for (seq_idx, (name, seq)) in sequences.into_iter().enumerate() {
                        let seq_start = Instant::now();
                        let seq_len = seq.len();
                        info!("  Processing {} ({} bp)", name, seq_len);
                        let stats = index.add_sequence_with_existing_syncmers(name.clone(), seq)?;
                        sequence_count += 1;
                        total_syncmers += stats.syncmers;
                        info!(
                            "  Indexed sequence {}/{} {}: {} bp, {} syncmers, indexed={}, elapsed {}, total {} sequences / {} bp / {} syncmers at +{}",
                            seq_idx + 1,
                            total_sequences,
                            name,
                            stats.sequence_len,
                            stats.syncmers,
                            stats.indexed,
                            format_duration(seq_start.elapsed()),
                            sequence_count,
                            total_bp,
                            total_syncmers,
                            format_duration(total_start.elapsed())
                        );
                    }

                    info!(
                        "Built syng paths from {} sequences ({} bp, {} syncmers) in {}",
                        sequence_count,
                        total_bp,
                        total_syncmers,
                        format_duration(input_start.elapsed())
                    );
                    index
                } else {
                // Read sequences from FASTA
                let input_start = Instant::now();
                info!(
                    "Reading sequences from FASTA: {} at +{}",
                    fasta_path,
                    format_duration(total_start.elapsed())
                );

                let file = File::open(&fasta_path)?;
                let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
                    io::Error::other(format!(
                        "Failed to open reader for '{}': {}",
                        fasta_path, e
                    ))
                })?;
                let reader = BufReader::new(reader);

                let mut index = impg::syng::SyngIndex::new(params);
                configure_position_sampling(&mut index)?;
                let mut current_name = String::new();
                let mut current_seq = Vec::new();
                let mut sequence_count = 0usize;
                let mut total_bp = 0u64;
                let mut total_syncmers = 0usize;

                for line in reader.lines() {
                    let line = line?;
                    if line.starts_with('>') {
                        if !current_name.is_empty() && !current_seq.is_empty() {
                            let seq_start = Instant::now();
                            let seq_len = current_seq.len();
                            info!(
                                "  Processing {} ({} bp)",
                                current_name,
                                seq_len
                            );
                            let name = std::mem::take(&mut current_name);
                            let stats = index.add_sequence(
                                name.clone(),
                                std::mem::take(&mut current_seq),
                            );
                            sequence_count += 1;
                            total_bp += stats.sequence_len as u64;
                            total_syncmers += stats.syncmers;
                            info!(
                                "  Indexed {}: {} bp, {} syncmers, indexed={}, elapsed {}, total {} sequences / {} bp / {} syncmers at +{}",
                                name,
                                stats.sequence_len,
                                stats.syncmers,
                                stats.indexed,
                                format_duration(seq_start.elapsed()),
                                sequence_count,
                                total_bp,
                                total_syncmers,
                                format_duration(total_start.elapsed())
                            );
                        }
                        current_name = line
                            .strip_prefix('>')
                            .unwrap_or("")
                            .split_whitespace()
                            .next()
                            .unwrap_or("")
                            .to_string();
                        current_seq.clear();
                    } else {
                        current_seq.extend(line.trim().as_bytes());
                    }
                }
                if !current_name.is_empty() && !current_seq.is_empty() {
                    let seq_start = Instant::now();
                    let seq_len = current_seq.len();
                    info!(
                        "  Processing {} ({} bp)",
                        current_name,
                        seq_len
                    );
                    let stats = index.add_sequence(current_name.clone(), current_seq);
                    sequence_count += 1;
                    total_bp += stats.sequence_len as u64;
                    total_syncmers += stats.syncmers;
                    info!(
                        "  Indexed {}: {} bp, {} syncmers, indexed={}, elapsed {}, total {} sequences / {} bp / {} syncmers at +{}",
                        current_name,
                        stats.sequence_len,
                        stats.syncmers,
                        stats.indexed,
                        format_duration(seq_start.elapsed()),
                        sequence_count,
                        total_bp,
                        total_syncmers,
                        format_duration(total_start.elapsed())
                    );
                }

                info!(
                    "Built syng paths from {} sequences ({} bp, {} syncmers) in {}",
                    sequence_count,
                    total_bp,
                    total_syncmers,
                    format_duration(input_start.elapsed())
                );
                index
                }
            } else {
                unreachable!()
            };

            let position_start = Instant::now();
            if let Some(stats) = index.finalize_online_sampled_positions()? {
                info!(
                    "Regular sampled position sidecars compacted in {}: {} node-position samples across {} nodes; {} path-step checkpoints across {} paths from {} walked paths (sample_shift={}, period={} syncmer steps) at +{}",
                    format_duration(position_start.elapsed()),
                    stats.sampled_occurrences,
                    stats.sampled_nodes,
                    stats.sampled_path_steps,
                    stats.sampled_step_paths,
                    stats.walked_paths,
                    stats.sample_shift,
                    1u64 << stats.sample_shift,
                    format_duration(total_start.elapsed()),
                );
            }

            let save_start = Instant::now();
            info!("Saving index to prefix: {}", output);
            index.save(&output)?;
            info!(
                "Index saved in {}: {}.1khash, {}.1gbwt, {}, {}, {}, {}",
                format_duration(save_start.elapsed()),
                output,
                output,
                impg::syng::syng_names_path(&output),
                impg::syng::syng_spos_path(&output),
                impg::syng::syng_pstep_path(&output),
                impg::syng::syng_meta_path(&output),
            );
            info!(
                "Name map contains {} sequences",
                index.name_map.path_to_name.len()
            );
            info!("Total syng build time: {}", format_duration(total_start.elapsed()));
        }
        Args::Syng2gfa {
            syng_prefix,
            output,
            gfa_version,
            sequence,
            common,
        } => {
            initialize_threads_and_log(&common);
            let version = impg::commands::syng2gfa::GfaVersion::parse(&gfa_version)?;
            let sequence_files = sequence.resolve_sequence_files()?;
            impg::commands::syng2gfa::run(&syng_prefix, &output, version, &sequence_files)?;
        }
        Args::SyngRepair {
            index,
            force,
            position_sample_shift,
            position_sample_seed,
            serial_position_sampling,
            position_progress_interval,
            common,
        } => {
            initialize_threads_and_log(&common);
            let prefix = detect_syng_prefix(&index).unwrap_or(index);
            let spos_path = impg::syng::syng_spos_path(&prefix);
            let pstep_path = impg::syng::syng_pstep_path(&prefix);
            let have_spos = std::path::Path::new(&spos_path).exists();
            let have_pstep = std::path::Path::new(&pstep_path).exists();
            if !force && have_spos && have_pstep {
                info!(
                    "Syng position sidecars already exist: {}, {}. Use --force to rebuild regular-grid sidecars.",
                    spos_path, pstep_path
                );
                return Ok(());
            }

            let total_start = Instant::now();
            info!(
                "Loading syng index for positional sidecar repair from prefix: {}",
                prefix
            );
            #[cfg(unix)]
            let saved_stdout = silence_stdout_for_process()?;
            #[cfg(not(unix))]
            silence_stdout_for_process()?;
            let mut syng_index = impg::syng::SyngIndex::load_for_repair(
                &prefix,
                impg::syng::SyncmerParams::default(),
            )?;
            #[cfg(unix)]
            {
                if unsafe { libc::dup2(saved_stdout, libc::STDOUT_FILENO) } < 0 {
                    let e = io::Error::last_os_error();
                    unsafe {
                        libc::close(saved_stdout);
                    }
                    return Err(e);
                }
                unsafe {
                    libc::close(saved_stdout);
                }
            }

            let rebuild_start = Instant::now();
            if position_sample_shift >= 63 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--position-sample-shift must be less than 63",
                ));
            }
            let sample_period = 1u64 << position_sample_shift;
            if have_spos && !force {
                info!(
                    "Rebuilding regular sampled path-step sidecar: sample_shift={} (every {} syncmer steps per path)",
                    position_sample_shift,
                    sample_period,
                );
                let stats = if serial_position_sampling {
                    syng_index.rebuild_sampled_path_steps_from_gbwt(
                        position_sample_shift,
                        position_sample_seed,
                    )?
                } else {
                    syng_index.rebuild_sampled_path_steps_from_gbwt_parallel(
                        position_sample_shift,
                        position_sample_seed,
                        position_progress_interval,
                    )?
                };
                info!(
                    "Rebuilt sampled path-step sidecar in {}: {} checkpoints across {} paths from {} walked paths",
                    format_duration(rebuild_start.elapsed()),
                    stats.sampled_path_steps,
                    stats.sampled_step_paths,
                    stats.walked_paths
                );
                let save_start = Instant::now();
                syng_index.save_path_step_sidecar(&prefix)?;
                info!(
                    "Saved path-step sidecar in {}: {}",
                    format_duration(save_start.elapsed()),
                    impg::syng::syng_pstep_path(&prefix)
                );
            } else {
                info!(
                    "Rebuilding regular sampled position sidecars: sample_shift={} (every {} syncmer steps per path)",
                    position_sample_shift,
                    sample_period,
                );
                let stats = if serial_position_sampling {
                    syng_index.rebuild_sampled_position_indexes_from_gbwt(
                        position_sample_shift,
                        position_sample_seed,
                    )?
                } else {
                    syng_index.rebuild_sampled_position_indexes_from_gbwt_parallel(
                        position_sample_shift,
                        position_sample_seed,
                        position_progress_interval,
                    )?
                };
                info!(
                    "Rebuilt sampled position sidecars in {}: {} node-position samples across {} nodes; {} path-step checkpoints across {} paths from {} walked paths",
                    format_duration(rebuild_start.elapsed()),
                    stats.sampled_occurrences,
                    stats.sampled_nodes,
                    stats.sampled_path_steps,
                    stats.sampled_step_paths,
                    stats.walked_paths
                );
                let save_start = Instant::now();
                syng_index.save_position_sidecars(&prefix)?;
                info!(
                    "Saved positional sidecars in {}: {}, {}",
                    format_duration(save_start.elapsed()),
                    impg::syng::syng_spos_path(&prefix),
                    impg::syng::syng_pstep_path(&prefix)
                );
            }
            info!(
                "Total syng sidecar repair time: {}",
                format_duration(total_start.elapsed())
            );
        }
    }

    Ok(())
}

fn build_graph_config(
    engine_cli: &EngineCliOpts,
    sparsify: &SparsificationStrategy,
    paf_file: Option<String>,
    temp_dir: &str,
    num_threads: usize,
    show_progress: bool,
) -> io::Result<graph::GraphBuildConfig> {
    Ok(graph::GraphBuildConfig {
        num_threads,
        frequency_multiplier: engine_cli.aln.sw.fastga_frequency_multiplier,
        frequency: engine_cli.aln.sw.frequency,
        min_aln_length: engine_cli.aln.sw.block_length.unwrap_or(0),
        repeat_max: engine_cli.seqwish.repeat_max,
        min_repeat_dist: engine_cli.seqwish.min_repeat_dist,
        min_match_len: engine_cli.seqwish.min_match_len,
        sparse_factor: engine_cli.seqwish.sparse_factor,
        transclose_batch: engine_cli.seqwish.transclose_batch,
        disk_backed: engine_cli.seqwish.disk_backed,
        show_progress,
        temp_dir: Some(temp_dir.to_string()),
        input_paf: paf_file,
        aligner: engine_cli.aln.sw.aligner.clone(),
        no_filter: engine_cli.aln.sw.no_filter,
        num_mappings: engine_cli.aln.sw.num_mappings.clone(),
        scaffold_jump: engine_cli.aln.sw.scaffold_jump,
        scaffold_mass: engine_cli.aln.sw.scaffold_mass,
        scaffold_filter: engine_cli.aln.sw.scaffold_filter.clone(),
        overlap: engine_cli.aln.sw.overlap,
        min_identity: sweepga::parse_identity_value(&engine_cli.aln.sw.min_identity, None)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?,
        scaffold_dist: engine_cli.aln.sw.scaffold_dist,
        min_map_length: engine_cli.aln.min_map_length,
        debug_dir: None,
        sparsify: sparsify.clone(),
        mash_params: sweepga::knn_graph::MashParams {
            kmer_size: engine_cli.aln.sw.mash_kmer_size,
            sketch_size: engine_cli.aln.sw.mash_sketch_size,
        },
        batch_bytes: engine_cli.aln.sw.batch_bytes.clone(),
        map_pct_identity: engine_cli
            .aln
            .sw
            .map_pct_identity
            .clone()
            .or_else(|| Some("90".to_string())),
    })
}

fn build_engine_opts(
    engine: GfaEngine,
    num_threads: usize,
    aln: &AlnOpts,
    sparsify: SparsificationStrategy,
    mash_params: sweepga::knn_graph::MashParams,
    seqwish: &SeqwishOpts,
    smooth: &SmoothOpts,
    debug_dir: Option<String>,
    temp_dir: Option<String>,
    partition_size: Option<usize>,
) -> io::Result<EngineOpts> {
    let pipeline = graph::GraphBuildConfig {
        num_threads,
        no_filter: aln.sw.no_filter,
        debug_dir,
        sparsify,
        mash_params,
        aligner: aln.sw.aligner.clone(),
        num_mappings: aln.sw.num_mappings.clone(),
        scaffold_jump: aln.sw.scaffold_jump,
        scaffold_mass: aln.sw.scaffold_mass,
        scaffold_filter: aln.sw.scaffold_filter.clone(),
        overlap: aln.sw.overlap,
        min_identity: sweepga::parse_identity_value(&aln.sw.min_identity, None)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?,
        scaffold_dist: aln.sw.scaffold_dist,
        min_map_length: aln.min_map_length,
        min_aln_length: aln.sw.block_length.unwrap_or(0),
        frequency_multiplier: aln.sw.fastga_frequency_multiplier,
        repeat_max: seqwish.repeat_max,
        min_repeat_dist: seqwish.min_repeat_dist,
        min_match_len: seqwish.min_match_len,
        sparse_factor: seqwish.sparse_factor,
        transclose_batch: seqwish.transclose_batch,
        disk_backed: seqwish.disk_backed,
        temp_dir,
        batch_bytes: aln.sw.batch_bytes.clone(),
        // The partitioned pipeline already logs its own progress
        show_progress: false,
        ..graph::GraphBuildConfig::default()
    };

    Ok(EngineOpts {
        engine,
        pipeline,
        partition_size,
        target_poa_lengths: smooth.parse_target_poa_lengths()?,
        max_node_length: smooth.max_node_length,
        poa_padding_fraction: smooth.poa_padding_fraction,
    })
}

fn validate_approximate_mode_min_length(
    min_transitive_len_opt: Option<i32>,
    is_transitive: bool,
) -> io::Result<()> {
    // For transitive queries in approximate mode, warn if min-transitive-len is not set
    // (coordinates are coarse-grained from tracepoint boundaries)
    if is_transitive && min_transitive_len_opt.is_none() {
        log::warn!(
            "--approximate mode with transitive queries uses default --min-transitive-len=101. \
            For fastga tracepoints, consider setting it greater than your trace_spacing."
        );
    }

    Ok(())
}

/// Validate that a range meets the minimum length requirement
fn validate_range_min_length(
    start: i32,
    end: i32,
    range_name: &str,
    min_transitive_len: i32,
) -> io::Result<()> {
    let length = end - start;
    if length < min_transitive_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Range '{range_name}' ({length} bp) is below minimum of {min_transitive_len} bp. Lower --min-transitive-len or use a longer range"
            ),
        ));
    }
    Ok(())
}

fn validate_selection_mode(mode: &str) -> io::Result<()> {
    match mode {
        "longest" | "total" => Ok(()),
        mode if mode == "sample" || mode == "haplotype"
            || mode.starts_with("sample,") || mode.starts_with("haplotype,") => Ok(()),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Invalid selection mode. Must be 'longest', 'total', 'sample[,sep]', or 'haplotype[,sep]'."
        ))
    }
}

/// Parse POA scoring string "match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2"
fn parse_poa_scoring_string(s: &str) -> io::Result<(u8, u8, u8, u8, u8, u8)> {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() != 6 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "POA scoring format should be 'match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2'",
        ));
    }
    let parse_u8 = |s: &str, name: &str| {
        s.parse::<u8>().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidInput, format!("Invalid {name} value"))
        })
    };
    Ok((
        parse_u8(parts[0], "match score")?,
        parse_u8(parts[1], "mismatch cost")?,
        parse_u8(parts[2], "gap opening 1 cost")?,
        parse_u8(parts[3], "gap extension 1 cost")?,
        parse_u8(parts[4], "gap opening 2 cost")?,
        parse_u8(parts[5], "gap extension 2 cost")?,
    ))
}

fn validate_output_format(format: &str, valid_formats: &[&str]) -> io::Result<()> {
    if valid_formats.contains(&format) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Invalid output format '{}'. Must be one of: {}",
                format,
                valid_formats.join(", ")
            ),
        ))
    }
}

/// Validates that a sequence name exists in the index and the range is within bounds
fn validate_sequence_range(
    seq_name: &str,
    start: i32,
    end: i32,
    seq_index: &SequenceIndex,
) -> io::Result<()> {
    // Check if sequence exists in index
    let seq_id = seq_index.get_id(seq_name).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!("Sequence '{seq_name}' not found in index"),
        )
    })?;

    // Get sequence length
    let seq_len = seq_index.get_len_from_id(seq_id).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Could not get length for sequence '{seq_name}'"),
        )
    })? as i32;

    // Validate range bounds
    if start < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Start position {start} cannot be negative"),
        ));
    }

    if end < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("End position {end} cannot be negative"),
        ));
    }

    if start >= end {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Start position {start} must be less than end position {end}"),
        ));
    }

    if end > seq_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "End position {end} exceeds sequence length {seq_len} for sequence '{seq_name}'"
            ),
        ));
    }

    Ok(())
}

/// Validate region size
fn validate_region_size(
    start: i32,
    end: i32,
    output_format: &str,
    merge_distance: i32,
    force_large_region: bool,
) -> io::Result<()> {
    let region_size = (end - start).unsigned_abs() as u64;
    const SIZE_LIMIT: u64 = 10_000; // 10kbp limit
    const MERGE_DISTANCE_LIMIT: i32 = 1000; // 1k limit

    // Check if this is a format that uses flat single-pass SPOA (may be slow/OOM on large regions)
    let uses_spoa = matches!(output_format, "maf" | "gfa-poa" | "fasta-aln");

    if uses_spoa && !force_large_region {
        if region_size > SIZE_LIMIT {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Region size ({region_size} bp) exceeds 10kbp for '{output_format}' output format, which may require large time and memory. Use --force-large-region to proceed anyway."
                ),
            ));
        }

        if merge_distance > MERGE_DISTANCE_LIMIT {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Merge distance ({merge_distance} bp) exceeds 1kbp for '{output_format}' output format, which may require large time and memory. Use --force-large-region to proceed anyway."
                ),
            ));
        }
    }

    Ok(())
}

/// Initialize thread pool and logger with the specified number of threads and verbosity
fn initialize_threads_and_log(common: &CommonOpts) {
    // Initialize logger based on verbosity
    env_logger::Builder::new()
        .filter_level(match common.verbose {
            0 => log::LevelFilter::Error,
            1 => log::LevelFilter::Info,
            _ => log::LevelFilter::Debug,
        })
        .init();

    // Configure thread pool
    ThreadPoolBuilder::new()
        .num_threads(common.threads.into())
        .build_global()
        .unwrap();
}

/// Determine file format from explicit format or auto-detection
fn determine_file_format(
    format: &str,
    files: &Option<Vec<String>>,
    file_list: &Option<String>,
) -> io::Result<String> {
    if format != "auto" {
        return Ok(format.to_string());
    }

    // Auto-detect from first file
    let first_file = match (files, file_list) {
        (Some(files), _) if !files.is_empty() => files[0].clone(),
        (_, Some(list_file)) => {
            // Read first line from list file
            let content = std::fs::read_to_string(list_file)?;
            let first_line = content
                .lines()
                .find(|line| !line.trim().is_empty() && !line.trim().starts_with('#'))
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "No valid files found in list file",
                    )
                })?;
            first_line.trim().to_string()
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "No input files provided",
            ))
        }
    };

    // Auto-detect based on file extension
    if first_file.ends_with(".vcf")
        || first_file.ends_with(".vcf.gz")
        || first_file.ends_with(".vcf.zst")
        || first_file.ends_with(".vcf.bgz")
    {
        Ok("vcf".to_string())
    } else if first_file.ends_with(".gfa")
        || first_file.ends_with(".gfa.gz")
        || first_file.ends_with(".gfa.zst")
        || first_file.ends_with(".gfa.bgz")
    {
        Ok("gfa".to_string())
    } else {
        // Try to detect by reading first few lines through decompression
        let reader = get_auto_reader(&first_file)?;

        for line in reader.lines().take(10) {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            // Check for VCF header
            if line.starts_with("#CHROM") {
                return Ok("vcf".to_string());
            }

            // Check for GFA segments, links, or paths
            if line.starts_with('S')
                || line.starts_with('L')
                || line.starts_with('P')
                || line.starts_with('H')
            {
                return Ok("gfa".to_string());
            }
        }

        // Default to GFA if we can't determine
        warn!("Could not auto-detect file format for '{first_file}', defaulting to GFA");
        Ok("gfa".to_string())
    }
}

// Helper function for auto-detection with compression support
fn get_auto_reader(path: &str) -> io::Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(path)?;
    let (reader, _format) = niffler::get_reader(Box::new(file))
        .map_err(|e| io::Error::other(format!("Failed to open reader: {e}")))?;
    Ok(Box::new(BufReader::new(reader)))
}

/// Helper function to return a Write implementer that is either standard output or a file with the
/// appropriate basename and extension. When no basename is provided, uses standard output.
fn find_output_stream(basename: &Option<String>, extension: &str) -> io::Result<Box<dyn Write>> {
    match basename {
        Some(name) => {
            let filename = format!("{}.{}", name, extension);
            let file = File::create(filename)?;
            Ok(Box::new(BufWriter::new(file)))
        }
        None => Ok(Box::new(BufWriter::new(io::stdout()))),
    }
}

/// Load/generate index based on common and alignment options
fn initialize_index(
    common: &CommonOpts,
    alignment: &AlignmentOpts,
    alignment_files: &[String],
    sequence_files: Vec<String>,
) -> io::Result<ImpgWrapper> {
    // The list of alignment files (PAF, 1ALN, or TPA) is pre-resolved
    info!("Found {} alignment file(s)", alignment_files.len());

    let seq_files_opt = if sequence_files.is_empty() {
        None
    } else {
        Some(sequence_files.as_slice())
    };

    // Determine index mode: explicit mode takes priority, then auto threshold
    let use_per_file = match alignment.index_mode {
        IndexMode::PerFile => {
            if alignment.index.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Cannot use -i/--index with --index-mode per-file",
                ));
            }
            true
        }
        IndexMode::Single => false,
        IndexMode::Auto => {
            if alignment.index.is_some() {
                // Explicit -i path implies single index
                false
            } else if alignment_files.len() >= 100 {
                warn!(
                    "Auto-switching to per-file indexing ({} files >= 100). \
                     Override with --index-mode single or --index-mode per-file.",
                    alignment_files.len(),
                );
                true
            } else {
                false
            }
        }
    };

    if use_per_file {
        load_or_build_per_file_index(
            alignment_files,
            common.threads,
            alignment.force_reindex,
            seq_files_opt,
            alignment.alignment_list.as_deref(),
            !alignment.unidirectional,
        )
    } else {
        let impg = load_or_build_single_index(
            alignment_files,
            common.threads,
            alignment.index.as_deref(),
            seq_files_opt,
            alignment.force_reindex,
            !alignment.unidirectional,
        )?;
        Ok(ImpgWrapper::from_single(impg))
    }
}

/// Initialize MultiImpg from per-file indices
fn load_or_build_per_file_index(
    alignment_files: &[String],
    threads: NonZeroUsize,
    force_reindex: bool,
    sequence_files: Option<&[String]>,
    alignment_list: Option<&str>,
    bidirectional: bool,
) -> io::Result<ImpgWrapper> {
    use indicatif::{ProgressBar, ProgressStyle};
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicUsize, Ordering};

    // Generate per-file index paths
    let index_paths: Vec<PathBuf> = alignment_files
        .iter()
        .map(|f| PathBuf::from(format!("{}.impg", f)))
        .collect();

    // Check which indices need to be built
    let indices_to_build: Vec<usize> = if force_reindex {
        (0..alignment_files.len()).collect()
    } else {
        index_paths
            .iter()
            .enumerate()
            .filter(|(_, path)| !path.exists())
            .map(|(i, _)| i)
            .collect()
    };

    let direction = if bidirectional {
        "bidirectional"
    } else {
        "unidirectional"
    };
    let total = alignment_files.len();
    let to_build = indices_to_build.len();
    if to_build == 0 {
        info!(
            "Using per-file indexing mode ({} {} index file(s) up to date)",
            total, direction
        );
    } else if force_reindex {
        info!("Using per-file indexing mode (force rebuild)");
        info!(
            "Building {} of {} {} index file(s)...",
            to_build, total, direction
        );
    } else {
        info!("Using per-file indexing mode");
        info!(
            "Building {} of {} {} index file(s) ({} already up to date)...",
            to_build,
            total,
            direction,
            total - to_build
        );
    }

    if !indices_to_build.is_empty() {
        // Create progress bar at info level (not at error-only or debug level)
        let pb = if log::log_enabled!(log::Level::Info) && !log::log_enabled!(log::Level::Debug) {
            let progress_bar = ProgressBar::new(indices_to_build.len() as u64);
            progress_bar.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
                    .unwrap()
                    .progress_chars("#>-")
            );
            Some(progress_bar)
        } else {
            None
        };

        let completed = AtomicUsize::new(0);

        // Build indices in parallel
        let build_results: Vec<io::Result<()>> = indices_to_build
            .par_iter()
            .map(|&i| {
                let aln_file = &alignment_files[i];
                let index_path = &index_paths[i];

                debug!("Building index for {}", aln_file);

                // Build index for this single alignment file (no progress bar, tracked at per-file level)
                let single_file = vec![aln_file.clone()];
                let impg = build_single_index(
                    &single_file,
                    threads,
                    Some(index_path.to_string_lossy().as_ref()),
                    sequence_files,
                    false,
                    bidirectional,
                )?;

                // Save the index
                let file = File::create(index_path)?;
                let mut writer = BufWriter::new(file);
                impg.serialize_with_forest_map(&mut writer)?;

                debug!("Built index: {:?}", index_path);

                // Update progress
                completed.fetch_add(1, Ordering::Relaxed);
                if let Some(ref progress_bar) = pb {
                    progress_bar.inc(1);
                }

                Ok(())
            })
            .collect();

        if let Some(progress_bar) = pb {
            progress_bar.finish_with_message("Completed building indices");
        }

        // Check for errors
        for result in build_results {
            result?;
        }
    }

    // Load MultiImpg from all per-file indices
    let multi = if let Some(list_path) = alignment_list {
        // Use cache when alignment list file is available
        MultiImpg::load_with_cache(
            &index_paths,
            alignment_files,
            sequence_files,
            Path::new(list_path),
        )?
    } else {
        // No list file, can't use cache (e.g., --alignment-files mode)
        MultiImpg::load_from_files(&index_paths, alignment_files, sequence_files)?
    };

    Ok(ImpgWrapper::from_multi(multi))
}

/// Resolve the list of alignment files (PAF, 1ALN, or TPA) from either --alignment-files or --alignment-list
fn resolve_alignment_files(alignment: &AlignmentOpts) -> io::Result<Vec<String>> {
    let alignment_files = if !alignment.alignment_files.is_empty() {
        alignment.alignment_files.clone()
    } else if let Some(alignment_list_file) = &alignment.alignment_list {
        // Read alignment files from the list file
        let file = File::open(alignment_list_file)?;
        let reader = BufReader::new(file);
        let mut files = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                files.push(trimmed.to_string());
            }
        }

        if files.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("No valid alignment files found in list file: {alignment_list_file}"),
            ));
        }

        files
    } else {
        // Neither alignment_files nor alignment_list provided
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Either --alignment-files or --alignment-list must be provided",
        ));
    };

    // Check if the number of PAF files exceeds u32::MAX
    if alignment_files.len() > u32::MAX as usize {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Too many PAF files specified: {} (maximum allowed: {})",
                alignment_files.len(),
                u32::MAX
            ),
        ));
    }

    Ok(alignment_files)
}

fn load_or_build_single_index(
    alignment_files: &[String],
    threads: NonZeroUsize,
    custom_index: Option<&str>,
    sequence_files: Option<&[String]>,
    force_reindex: bool,
    bidirectional: bool,
) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(alignment_files, custom_index);

    let direction = if bidirectional {
        "bidirectional"
    } else {
        "unidirectional"
    };
    let per_file_hint = if custom_index.is_some() && alignment_files.len() >= 100 {
        format!(
            " (drop -i or use --index-mode per-file for per-file indexing with {} files)",
            alignment_files.len()
        )
    } else {
        String::new()
    };

    if force_reindex {
        info!(
            "Using single indexing mode (force rebuild){}",
            per_file_hint
        );
        info!(
            "Building 1 {} index file processing {} alignment file(s)...",
            direction,
            alignment_files.len()
        );
        return build_single_index(
            alignment_files,
            threads,
            custom_index,
            sequence_files,
            true,
            bidirectional,
        );
    }

    if !std::path::Path::new(&index_file).exists() {
        info!("Using single indexing mode{}", per_file_hint);
        info!(
            "Building 1 {} index file processing {} alignment file(s)...",
            direction,
            alignment_files.len()
        );
        return build_single_index(
            alignment_files,
            threads,
            custom_index,
            sequence_files,
            true,
            bidirectional,
        );
    }

    // Load existing index
    info!("Reading IMPG index from {index_file}");

    // Check if all alignment files are newer than the index
    let index_file_metadata = std::fs::metadata(&index_file)?;
    let index_file_ts = index_file_metadata.modified().ok();

    if let Some(index_ts) = index_file_ts {
        alignment_files.par_iter().for_each(|alignment_file| {
            if let Ok(alignment_file_metadata) = std::fs::metadata(alignment_file) {
                if let Ok(alignment_file_ts) = alignment_file_metadata.modified() {
                    if alignment_file_ts > index_ts {
                        warn!(
                            "WARNING:\tAlignment file {alignment_file} has been modified since impg index creation."
                        );
                    }
                }
            }
        });
    }

    // Load from embedded format
    let file = File::open(&index_file)?;
    let reader = BufReader::new(file);

    Impg::load_from_file(reader, alignment_files, index_file, sequence_files)
}

fn build_single_index(
    alignment_files: &[String],
    threads: NonZeroUsize,
    custom_index: Option<&str>,
    sequence_files: Option<&[String]>,
    show_progress: bool,
    bidirectional: bool,
) -> io::Result<Impg> {
    use indicatif::{ProgressBar, ProgressStyle};

    let index_file = get_combined_index_filename(alignment_files, custom_index);
    debug!("Creating index: {index_file}");

    if bidirectional {
        debug!("Building bidirectional index");
    } else {
        debug!("Building unidirectional index");
    }

    let num_alignment_files = alignment_files.len();

    // Create progress bar at info level (not at error-only or debug level)
    let pb = if show_progress
        && log::log_enabled!(log::Level::Info)
        && !log::log_enabled!(log::Level::Debug)
    {
        let progress_bar = ProgressBar::new(num_alignment_files as u64);
        progress_bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap()
                .progress_chars("#>-")
        );
        Some(progress_bar)
    } else {
        None
    };

    // Process alignment files in parallel, using a per-file local SequenceIndex to avoid global locking
    // We will merge names after parsing and remap IDs in a parallel pass.
    let mut records_by_file_with_local_index: Vec<(Vec<AlignmentRecord>, String, SequenceIndex)> =
        (0..alignment_files.len())
            .into_par_iter()
            .map(
                |file_index| -> io::Result<(Vec<AlignmentRecord>, String, SequenceIndex)> {
                    let aln_file = &alignment_files[file_index];

                    debug!(
                        "Processing alignment file ({}/{num_alignment_files}): {aln_file}",
                        file_index + 1
                    );

                    // Local sequence index for this file only
                    let mut local_seq_index = SequenceIndex::new();

                    // Detect file format and parse accordingly
                    let format = AlignmentFormat::from_path(aln_file).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("Unsupported alignment format: {aln_file}"),
                        )
                    })?;

                    let records = match format {
                        AlignmentFormat::Paf => {
                            let file = File::open(aln_file).map_err(|e| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("Failed to open PAF file: {}", e),
                                )
                            })?;
                            impg::paf::parse_paf_file(aln_file, file, threads, &mut local_seq_index)
                                .map_err(|e| {
                                    io::Error::new(
                                        io::ErrorKind::InvalidData,
                                        format!("Failed to parse PAF records: {}", e),
                                    )
                                })?
                        }
                        AlignmentFormat::OneAln => {
                            let parser = OneAlnParser::new(aln_file.clone(), sequence_files)
                                .map_err(|e| {
                                    io::Error::new(
                                        io::ErrorKind::InvalidData,
                                        format!("Failed to create 1aln parser: {}", e),
                                    )
                                })?;
                            parser.parse_alignments(&mut local_seq_index).map_err(|e| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("Failed to parse 1aln records: {}", e),
                                )
                            })?
                        }
                        AlignmentFormat::Tpa => {
                            let parser = TpaParser::new(aln_file.clone()).map_err(|e| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("Failed to create TPA parser: {}", e),
                                )
                            })?;
                            parser
                                .parse_alignments(&mut local_seq_index, threads.get())
                                .map_err(|e| {
                                    io::Error::new(
                                        io::ErrorKind::InvalidData,
                                        format!("Failed to parse TPA records: {}", e),
                                    )
                                })?
                        }
                    };
                    debug!(
                        "Parsed {} alignment records from file: {aln_file}",
                        records.len()
                    );

                    // Update progress bar
                    if let Some(ref progress_bar) = pb {
                        progress_bar.inc(1);
                    }

                    Ok((records, aln_file.clone(), local_seq_index))
                },
            )
            .collect::<Result<Vec<_>, _>>()?; // Propagate any errors

    if let Some(progress_bar) = pb {
        progress_bar.finish_with_message("Completed parsing alignment files");
    }

    // Build a temporary union of all sequence names and lengths from local indices
    let mut tmp_union_index = SequenceIndex::new();
    for (_, _, local_idx) in &records_by_file_with_local_index {
        for (name, id) in &local_idx.name_to_id {
            let length = local_idx.get_len_from_id(*id).unwrap_or(0);
            tmp_union_index.get_or_insert_id(name, Some(length));
        }
    }

    // Sort sequence names to ensure deterministic order
    let mut sequence_names = tmp_union_index
        .name_to_id
        .keys()
        .cloned()
        .collect::<Vec<String>>();
    sequence_names.par_sort_unstable(); // Order of identical sequence names is irrelevant

    // Create a deterministic global SequenceIndex
    let mut seq_index = SequenceIndex::new();
    for (name, id) in &tmp_union_index.name_to_id {
        let length = tmp_union_index.get_len_from_id(*id).unwrap_or(0);
        seq_index.get_or_insert_id(name, Some(length));
    }

    // Update query and target IDs with the new SequenceIndex
    records_by_file_with_local_index
        .par_iter_mut()
        .for_each(|(records, _path, local_idx)| {
            for record in records.iter_mut() {
                let query_name = local_idx.get_name(record.query_id).unwrap();
                record.query_id = seq_index.get_id(query_name).unwrap();

                let target_name = local_idx.get_name(record.target_id).unwrap();
                record.target_id = seq_index.get_id(target_name).unwrap();
            }
        });

    // Drop local indices and keep only (records, path)
    let records_by_file: Vec<(Vec<AlignmentRecord>, String)> = records_by_file_with_local_index
        .into_iter()
        .map(|(records, path, _)| (records, path))
        .collect();

    let impg = Impg::from_multi_alignment_records(
        &records_by_file,
        seq_index,
        sequence_files,
        bidirectional,
    )
    .map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Failed to create index: {e}"),
        )
    })?;

    // Serialize the index with embedded forest map
    let index_file_path = index_file.clone();
    let file = File::create(&index_file_path)?;
    let mut writer = BufWriter::new(file);
    impg.serialize_with_forest_map(&mut writer)?;

    Ok(impg)
}

fn get_combined_index_filename(alignment_files: &[String], custom_index: Option<&str>) -> String {
    if let Some(index) = custom_index {
        return index.to_string();
    }

    if alignment_files.len() == 1 {
        format!("{}.impg", alignment_files[0])
    } else {
        // For multiple files, create a hash of the sorted filenames

        let mut file_refs: Vec<&str> = alignment_files.iter().map(|s| s.as_str()).collect();
        file_refs.sort();

        let mut hasher = DefaultHasher::new();
        for file in &file_refs {
            file.hash(&mut hasher);
        }

        format!("combined_{:016x}.impg", hasher.finish())
    }
}

fn perform_query(
    impg: &impl ImpgIndex,
    target_name: &str,
    target_range: (i32, i32),
    store_cigar: bool,
    min_result_identity: Option<f64>,
    min_output_length: Option<i32>,
    transitive: bool,
    transitive_dfs: bool,
    transitive_opts: &TransitiveOpts,
    sequence_index: Option<&UnifiedSequenceIndex>,
    approximate_mode: bool,
    subset_filter: Option<&SubsetFilter>,
) -> io::Result<Vec<AdjustedInterval>> {
    let (target_start, target_end) = target_range;
    let target_id = impg.seq_index().get_id(target_name).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Target sequence '{target_name}' not found in index"),
        )
    })?;
    let target_length = impg.seq_index().get_len_from_id(target_id).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Target sequence '{target_name}' length not found in index"),
        )
    })?;
    if target_end > target_length as i32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Target range end ({target_end}) exceeds the target sequence length ({target_length})"
            ),
        ));
    }

    let results = if transitive {
        if transitive_dfs {
            impg.query_transitive_dfs(
                target_id,
                target_start,
                target_end,
                None,
                transitive_opts.max_depth,
                transitive_opts.effective_min_transitive_len(),
                transitive_opts.min_distance_between_ranges,
                min_output_length,
                store_cigar,
                min_result_identity,
                sequence_index,
                approximate_mode,
                subset_filter,
            )
        } else {
            impg.query_transitive_bfs(
                target_id,
                target_start,
                target_end,
                None,
                transitive_opts.max_depth,
                transitive_opts.effective_min_transitive_len(),
                transitive_opts.min_distance_between_ranges,
                min_output_length,
                store_cigar,
                min_result_identity,
                sequence_index,
                approximate_mode,
                subset_filter,
            )
        }
    } else {
        let mut res = impg.query(
            target_id,
            target_start,
            target_end,
            store_cigar,
            min_result_identity,
            sequence_index,
            approximate_mode,
        );
        // Filter by minimum output length for regular queries
        if let Some(min_len) = min_output_length {
            res.retain(|(query_interval, _, _)| {
                let length = (query_interval.last - query_interval.first).abs();
                length >= min_len
            });
        }

        // Apply subset filter for non-transitive queries (transitive queries filter during exploration)
        if let Some(filter) = subset_filter {
            apply_subset_filter(impg, target_id, &mut res, Some(filter));
        }

        res
    };

    info!(
        "Collected {} results (excluding input range)",
        results.len() - 1
    ); // Exclude the first element (the input range itself)

    Ok(results)
}

/// Load subset filter if path is provided
fn load_subset_filter_if_provided(path: &Option<String>) -> io::Result<Option<SubsetFilter>> {
    if let Some(ref list_path) = path {
        let filter = load_subset_filter(list_path)?;
        info!(
            "Loaded {} sequence names from subset list {}",
            filter.entry_count(),
            list_path
        );
        Ok(Some(filter))
    } else {
        Ok(None)
    }
}

/// Convert syng `HomologousInterval`s into `AdjustedInterval`s so syng output
/// can flow through the common merge + emit pipeline (output_results_bed,
/// output_results_bedpe). CIGAR is always empty — syng has no per-base
/// alignment. The first interval is the homologous result interval, matching
/// the alignment-backed path's BED/BEDPE emitters.
fn syng_intervals_to_adjusted(
    intervals: &[impg::syng::HomologousInterval],
    query_name: &str,
    range_start: i32,
    range_end: i32,
    seq_index: &SequenceIndex,
) -> Vec<AdjustedInterval> {
    let Some(query_id) = seq_index.get_id(query_name) else {
        return Vec::new();
    };
    intervals
        .iter()
        .filter_map(|iv| {
            let homolog_id = seq_index.get_id(&iv.genome)?;
            let homolog = if iv.strand == '-' {
                Interval::new(iv.end as i32, iv.start as i32, homolog_id)
            } else {
                Interval::new(iv.start as i32, iv.end as i32, homolog_id)
            };
            let target = Interval::new(range_start, range_end, query_id);
            Some((homolog, Vec::<CigarOp>::new(), target))
        })
        .collect()
}

fn output_results_bed(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    name: &str,
    merge_distance: i32,
    merge_strands: bool,
    original_coordinates: bool,
) -> io::Result<()> {
    let any_empty_cigar = results.iter().any(|(_, c, _)| c.is_empty());
    if any_empty_cigar {
        // 2D merge first — collapses fragmented chains on each target, then
        // falls through to the query-axis merge so downstream BED dedupe
        // across targets still works.
        merge_adjusted_intervals_gap_2d(results, merge_distance);
    }
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    for (query_interval, _, _) in results {
        let query_name = impg.seq_index().get_name(query_interval.metadata).unwrap();
        let (first, last, strand) = if query_interval.first <= query_interval.last {
            (query_interval.first, query_interval.last, '+')
        } else {
            (query_interval.last, query_interval.first, '-')
        };

        // Transform coordinates to original sequence space if requested
        let (transformed_name, transformed_first, transformed_last) =
            transform_coordinates_to_original(
                query_name,
                first as u32,
                last as u32,
                original_coordinates,
            );

        writeln!(
            out,
            "{}\t{}\t{}\t{}\t.\t{}",
            transformed_name, transformed_first, transformed_last, name, strand
        )?;
    }

    Ok(())
}

fn output_results_bedpe(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    name: &str,
    merge_distance: i32,
    original_coordinates: bool,
) -> io::Result<()> {
    // If any row lacks a CIGAR (syng output), use the gap-tolerant 2D merge
    // so fragmented chains on the same (query, target, strand) collapse into
    // one alignment. Otherwise keep the existing CIGAR-faithful merge.
    let any_empty_cigar = results.iter().any(|(_, c, _)| c.is_empty());
    if any_empty_cigar {
        merge_adjusted_intervals_gap_2d(results, merge_distance);
    } else {
        merge_adjusted_intervals(results, merge_distance);
    }

    for (overlap_query, cigar, overlap_target) in results {
        let query_name = impg.seq_index().get_name(overlap_query.metadata).unwrap();
        let target_name = impg.seq_index().get_name(overlap_target.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };

        // Transform coordinates to original sequence space if requested
        let (transformed_query_name, transformed_first, transformed_last) =
            transform_coordinates_to_original(
                query_name,
                first as u32,
                last as u32,
                original_coordinates,
            );
        let (transformed_target_name, transformed_target_first, transformed_target_last) =
            transform_coordinates_to_original(
                target_name,
                overlap_target.first as u32,
                overlap_target.last as u32,
                original_coordinates,
            );

        // Calculate gap-compressed-identity and block-identity from CIGAR
        let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, _block_len) =
            cigar.iter().fold(
                (0, 0, 0, 0, 0, 0, 0),
                |(m, mm, i, i_bp, d, d_bp, bl), op| {
                    let len = op.len();
                    match op.op() {
                        'M' => (m + len, mm, i, i_bp, d, d_bp, bl + len), // We overestimate num. of matches by assuming 'M' represents matches for simplicity
                        '=' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
                        'X' => (m, mm + len, i, i_bp, d, d_bp, bl + len),
                        'I' => (m, mm, i + 1, i_bp + len, d, d_bp, bl + len),
                        'D' => (m, mm, i, i_bp, d + 1, d_bp + len, bl + len),
                        _ => (m, mm, i, i_bp, d, d_bp, bl),
                    }
                },
            );
        let gap_compressed_identity =
            (matches as f32) / (matches + mismatches + insertions + deletions) as f32;

        let edit_distance = mismatches + inserted_bp + deleted_bp;
        let block_identity = (matches as f32) / (matches + edit_distance) as f32;

        // Format gi and bi fields without trailing zeros
        let gi_str = format!("{gap_compressed_identity:.6}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();
        let bi_str = format!("{block_identity:.6}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();

        // BEDPE supports any number of additional, fields after the standard 10 fields.
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t+\tgi:f:{}\tbi:f:{}",
            transformed_query_name,
            transformed_first,
            transformed_last,
            transformed_target_name,
            transformed_target_first,
            transformed_target_last,
            name,
            strand,
            gi_str,
            bi_str
        )?;
    }

    Ok(())
}

fn output_results_paf(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    name: &str,
    merge_distance: i32,
    original_coordinates: bool,
    sequence_index: Option<&UnifiedSequenceIndex>,
) -> io::Result<()> {
    merge_adjusted_intervals(results, merge_distance);

    for (overlap_query, cigar, overlap_target) in results {
        let query_name = impg.seq_index().get_name(overlap_query.metadata).unwrap();
        let target_name = impg.seq_index().get_name(overlap_target.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };

        // Transform coordinates to original sequence space if requested
        let (transformed_query_name, transformed_first, transformed_last) =
            transform_coordinates_to_original(
                query_name,
                first as u32,
                last as u32,
                original_coordinates,
            );
        let (transformed_target_name, transformed_target_first, transformed_target_last) =
            transform_coordinates_to_original(
                target_name,
                overlap_target.first as u32,
                overlap_target.last as u32,
                original_coordinates,
            );

        // Get original sequence lengths when original_sequence_coordinates is enabled
        let query_length = if original_coordinates {
            get_original_sequence_length(&transformed_query_name, sequence_index)
        } else {
            impg.seq_index()
                .get_len_from_id(overlap_query.metadata)
                .unwrap()
        };

        let target_length = if original_coordinates {
            get_original_sequence_length(&transformed_target_name, sequence_index)
        } else {
            impg.seq_index()
                .get_len_from_id(overlap_target.metadata)
                .unwrap()
        };

        let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len) =
            cigar.iter().fold(
                (0, 0, 0, 0, 0, 0, 0),
                |(m, mm, i, i_bp, d, d_bp, bl), op| {
                    let len = op.len();
                    match op.op() {
                        'M' => (m + len, mm, i, i_bp, d, d_bp, bl + len), // We overestimate num. of matches by assuming 'M' represents matches for simplicity
                        '=' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
                        'X' => (m, mm + len, i, i_bp, d, d_bp, bl + len),
                        'I' => (m, mm, i + 1, i_bp + len, d, d_bp, bl + len),
                        'D' => (m, mm, i, i_bp, d + 1, d_bp + len, bl + len),
                        _ => (m, mm, i, i_bp, d, d_bp, bl),
                    }
                },
            );
        let gap_compressed_identity =
            (matches as f32) / (matches + mismatches + insertions + deletions) as f32;

        let edit_distance = mismatches + inserted_bp + deleted_bp;
        let block_identity = (matches as f32) / (matches + edit_distance) as f32;

        // Format bi and gi fields without trailing zeros
        let gi_str = format!("{gap_compressed_identity:.6}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();
        let bi_str = format!("{block_identity:.6}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();

        let cigar_str: String = cigar
            .iter()
            .map(|op| format!("{}{}", op.len(), op.op()))
            .collect();

        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tgi:f:{}\tbi:f:{}\tcg:Z:{}\tan:Z:{}",
            transformed_query_name,
            query_length,
            transformed_first,
            transformed_last,
            strand,
            transformed_target_name,
            target_length,
            transformed_target_first,
            transformed_target_last,
            matches,
            block_len,
            255,
            gi_str,
            bi_str,
            cigar_str,
            name
        )?;
    }

    Ok(())
}

fn output_results_gfa(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    sequence_index: &UnifiedSequenceIndex,
    _name: &str,
    merge_distance: i32,
    merge_strands: bool,
    scoring_params: Option<(u8, u8, u8, u8, u8, u8)>,
    engine_opts: &EngineOpts,
) -> io::Result<()> {
    // Merge intervals if needed
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    // Extract query intervals by consuming results - no cloning
    let query_intervals: Vec<Interval<u32>> = results
        .drain(..)
        .map(|(query_interval, _, _)| query_interval)
        .collect();

    let gfa_output = impg::dispatch_gfa_engine(
        impg,
        &query_intervals,
        sequence_index,
        scoring_params,
        engine_opts,
    )?;
    out.write_all(gfa_output.as_bytes())?;

    Ok(())
}

/// Partitioned GFA output for the query command: splits the query region into
/// sub-windows of `partition_size` bp, queries each, then runs the partitioned
/// GFA pipeline (per-partition engine → lace → gfaffix).
fn output_results_gfa_partitioned(
    impg: &impl ImpgIndex,
    out: &mut dyn Write,
    sequence_index: &UnifiedSequenceIndex,
    scoring_params: Option<(u8, u8, u8, u8, u8, u8)>,
    engine_opts: &EngineOpts,
    target_name: &str,
    target_range: (i32, i32),
    partition_size: usize,
    query: &QueryOpts,
    subset_filter: Option<&SubsetFilter>,
) -> io::Result<()> {
    let (start, end) = target_range;
    let ps = partition_size as i32;

    // Split into sub-windows
    let mut partitions: Vec<Vec<Interval<u32>>> = Vec::new();
    let mut window_start = start;

    while window_start < end {
        let window_end = (window_start + ps).min(end);

        // Query this sub-window
        let mut results = perform_query(
            impg,
            target_name,
            (window_start, window_end),
            false, // no CIGAR needed for GFA
            query.min_result_identity,
            query.min_output_length,
            query.transitive,
            query.transitive_opts.transitive_dfs,
            &query.transitive_opts,
            Some(sequence_index),
            query.approximate,
            subset_filter,
        )?;

        if !results.is_empty() {
            // Merge intervals
            merge_query_adjusted_intervals(
                &mut results,
                query.effective_merge_distance(),
                query.merge_strands_for_output("gfa"),
            );

            // Extract query intervals
            let query_intervals: Vec<Interval<u32>> =
                results.drain(..).map(|(qi, _, _)| qi).collect();

            partitions.push(query_intervals);
        }

        window_start = window_end;
    }

    let gfa_output = impg::partitioned_gfa_pipeline(
        &partitions,
        impg,
        sequence_index,
        scoring_params,
        engine_opts,
    )?;
    out.write_all(gfa_output.as_bytes())?;

    Ok(())
}

fn output_results_fasta(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    sequence_index: &UnifiedSequenceIndex,
    _name: &str,
    merge_distance: i32,
    merge_strands: bool,
    reverse_complement: bool,
) -> io::Result<()> {
    // Merge intervals if needed
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    // Parallelize sequence fetching and processing
    let sequence_data: Vec<(String, String)> = results
        .par_iter()
        .map(|(query_interval, _, _)| -> io::Result<(String, String)> {
            let query_name = impg.seq_index().get_name(query_interval.metadata).unwrap();

            // Determine actual start and end based on orientation
            let (start, end, strand) = if query_interval.first <= query_interval.last {
                (query_interval.first, query_interval.last, '+')
            } else {
                (query_interval.last, query_interval.first, '-')
            };

            // Fetch the sequence
            let sequence = sequence_index.fetch_sequence(query_name, start, end)?;

            // If reverse strand and reverse complementing, reverse complement the sequence
            let sequence = if strand == '-' && reverse_complement {
                impg::graph::reverse_complement(&sequence)
            } else {
                sequence
            };

            // Create header
            let header_suffix = if strand == '-' && reverse_complement {
                "/rc"
            } else {
                ""
            };
            let header = format!(">{query_name}:{start}-{end}{header_suffix}");

            // Convert sequence to string with line breaks every 80 characters
            let sequence_str = String::from_utf8_lossy(&sequence);
            let formatted_sequence = sequence_str
                .as_bytes()
                .chunks(80)
                .map(|chunk| String::from_utf8_lossy(chunk).to_string())
                .collect::<Vec<_>>()
                .join("\n");

            Ok((header, formatted_sequence))
        })
        .collect::<Result<Vec<_>, _>>()?;

    // Output sequences sequentially to maintain order
    for (header, sequence) in sequence_data {
        writeln!(out, "{header}")?;
        writeln!(out, "{sequence}")?;
    }

    Ok(())
}

fn output_results_maf(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    sequence_index: &UnifiedSequenceIndex,
    _name: &str,
    merge_distance: i32,
    merge_strands: bool,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<()> {
    // Merge intervals if needed
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    let query_intervals: Vec<Interval<u32>> = results
        .drain(..)
        .map(|(query_interval, _, _)| query_interval)
        .collect();

    let maf_output = impg::graph::generate_maf_from_intervals(
        impg,
        &query_intervals,
        sequence_index,
        scoring_params,
    );
    writeln!(out, "{maf_output}")?;

    Ok(())
}

fn output_results_fasta_aln(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    sequence_index: &UnifiedSequenceIndex,
    _name: String,
    merge_distance: i32,
    merge_strands: bool,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<()> {
    // Merge intervals as for MAF/GFA (collapse per-query coords, merge strands)
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    // Drain query intervals
    let query_intervals: Vec<coitrees::Interval<u32>> =
        results.drain(..).map(|(q, _, _)| q).collect();

    // Ask graph layer to generate aligned FASTA from SPOA's MSA
    let fasta_aln = impg::graph::generate_fasta_alignment_from_intervals(
        impg,
        &query_intervals,
        sequence_index,
        scoring_params,
    );
    print!("{fasta_aln}");
    Ok(())
}

// Merge adjusted intervals by ignoring the target intervals (optimized for simple genomic interval merging in BED and GFA formats)
fn merge_query_adjusted_intervals(
    results: &mut Vec<AdjustedInterval>,
    merge_distance: i32,
    merge_strands: bool,
) {
    if results.len() > 1 && (merge_distance >= 0 || merge_strands) {
        // Sort by sequence ID, start position, and strand (forward first)
        results.par_sort_by_key(|(query_interval, _, _)| {
            let is_forward = query_interval.first <= query_interval.last;
            let start = if is_forward {
                query_interval.first
            } else {
                query_interval.last
            };

            (
                query_interval.metadata, // First sort by sequence ID
                start,                   // Then by actual start position
                !is_forward,             // Finally by strand orientation (forward first)
            )
        });

        let mut write_idx = 0;
        for read_idx in 1..results.len() {
            let (curr_interval, _, _) = &results[write_idx];
            let (next_interval, _, _) = &results[read_idx];

            // Orientation flags
            let curr_is_forward = curr_interval.first <= curr_interval.last;
            let next_is_forward = next_interval.first <= next_interval.last;

            // Extract actual start/end positions based on orientation
            let (curr_start, curr_end) = if curr_is_forward {
                (curr_interval.first, curr_interval.last)
            } else {
                (curr_interval.last, curr_interval.first)
            };

            let (next_start, next_end) = if next_is_forward {
                (next_interval.first, next_interval.last)
            } else {
                (next_interval.last, next_interval.first)
            };

            // Only merge if same sequence and within merge distance (if merge_distance >= 0),
            // and orientation is compatible with merge_strands choice.
            if merge_distance < 0
                || curr_interval.metadata != next_interval.metadata
                || (!merge_strands && curr_is_forward != next_is_forward)
                || next_start > curr_end + merge_distance
            {
                write_idx += 1;
                if write_idx != read_idx {
                    results.swap(write_idx, read_idx);
                }
            } else {
                // Merge intervals, possibly across strands.
                let merged_start = curr_start.min(next_start);
                let merged_end = curr_end.max(next_end);

                // When merging across strands, choose the orientation with the larger span; on ties, keep existing.
                let merged_is_forward = if merge_strands && curr_is_forward != next_is_forward {
                    let curr_len = curr_end.saturating_sub(curr_start);
                    let next_len = next_end.saturating_sub(next_start);
                    if next_len > curr_len {
                        next_is_forward
                    } else {
                        curr_is_forward
                    }
                } else {
                    curr_is_forward
                };

                if merged_is_forward {
                    results[write_idx].0.first = merged_start;
                    results[write_idx].0.last = merged_end;
                } else {
                    results[write_idx].0.first = merged_end;
                    results[write_idx].0.last = merged_start;
                }
            }
        }
        results.truncate(write_idx + 1);

        info!("Collected {} merged intervals", results.len());
    }
}

// Merge adjusted intervals by considering both query and target intervals and the corresponding CIGAR operations
fn merge_adjusted_intervals(results: &mut Vec<AdjustedInterval>, merge_distance: i32) {
    if results.len() > 1 && merge_distance >= 0 {
        // Sort by query ID, query position, target ID, target position
        results.par_sort_by_key(|(query_interval, _, target_interval)| {
            let query_forward = query_interval.first < query_interval.last;

            (
                query_interval.metadata, // Group by query sequence ID
                query_forward,           // Group by orientation (keep same orientations together)
                if query_forward {
                    // Use appropriate position based on orientation
                    query_interval.first // Forward: use start position
                } else {
                    query_interval.last // Reverse: use end position
                },
                target_interval.metadata, // Group by target sequence ID
                target_interval.first,    // Target always forward
            )
        });

        let num_results = results.len();

        // Drain off all entries, taking ownership
        let mut results_iter = results.drain(..);

        // Take the first as the current "in-progress" interval
        if let Some((mut current_query, mut current_cigar, mut current_target)) =
            results_iter.next()
        {
            // Create a new vector to store merged results
            let mut merged_results = Vec::with_capacity(num_results);

            // Iterate through remaining elements
            for (next_query, next_cigar, next_target) in results_iter {
                // Determine orientations
                let query_forward = current_query.first <= current_query.last;
                let next_query_forward = next_query.first <= next_query.last;

                let target_forward = current_target.first <= current_target.last;
                let next_target_forward = next_target.first <= next_target.last;
                if !target_forward || !next_target_forward {
                    panic!("Target intervals should always be in forward!");
                }

                // Check if sequences match and orientations are the same
                if current_query.metadata != next_query.metadata
                    || current_target.metadata != next_target.metadata
                    || query_forward != next_query_forward
                {
                    // Store current interval
                    merged_results.push((current_query, current_cigar, current_target));
                    // Clone the next as the new current
                    (current_query, current_cigar, current_target) =
                        (next_query, next_cigar, next_target);
                    continue;
                }

                // Check contiguity or overlap
                let (query_contiguous, target_contiguous, query_overlap, target_overlap) =
                    if query_forward {
                        let q_contig = current_query.last == next_query.first;
                        let t_contig = current_target.last == next_target.first;
                        let q_overlap = current_query.last > next_query.first;
                        let t_overlap = current_target.last > next_target.first;
                        (q_contig, t_contig, q_overlap, t_overlap)
                    } else {
                        // Reverse orientation (remember that first > last in reverse, so we swap first/last)
                        let q_contig = current_query.first == next_query.last;
                        let t_contig = current_target.first == next_target.last;
                        let q_overlap = current_query.first > next_query.last;
                        let t_overlap = current_target.first < next_target.last;
                        (q_contig, t_contig, q_overlap, t_overlap)
                    };

                // Handle perfect contiguity (existing logic)
                if query_contiguous && target_contiguous {
                    // debug!(
                    //     "Merge contiguous! Query: current {}:{}-{}({}), next {}:{}-{}({}); Target: current {}:{}-{}({}), next {}:{}-{}({})",
                    //     current_query.metadata,
                    //     current_query.first,
                    //     current_query.last,
                    //     if query_forward { "+" } else { "-" },
                    //     next_query.metadata,
                    //     next_query.first,
                    //     next_query.last,
                    //     if next_query_forward { "+" } else { "-" },
                    //     current_target.metadata,
                    //     current_target.first,
                    //     current_target.last,
                    //     if target_forward { "+" } else { "-" },
                    //     next_target.metadata,
                    //     next_target.first,
                    //     next_target.last,
                    //     if next_target_forward { "+" } else { "-" },
                    // );

                    // Merge intervals and CIGAR operations
                    if query_forward {
                        current_query.last = next_query.last;
                        current_target.last = next_target.last;
                        current_cigar.extend_from_slice(&next_cigar);
                    } else {
                        current_query.first = next_query.first;
                        current_target.first = next_target.first;

                        let mut new_cigar =
                            Vec::with_capacity(current_cigar.len() + next_cigar.len());
                        new_cigar.extend_from_slice(&next_cigar);
                        new_cigar.extend_from_slice(&current_cigar);
                        current_cigar = new_cigar;
                    }
                    merge_consecutive_cigar_ops(&mut current_cigar);
                    continue;
                }

                // Handle overlap case
                if query_overlap && target_overlap {
                    // Calculate overlap lengths
                    let (query_overlap_len, target_overlap_len) = if query_forward {
                        (
                            next_query.first - current_query.last,
                            next_target.first - current_target.last,
                        )
                    } else {
                        // Reverse orientation (remember that first > last in reverse, so we swap first/last)
                        (
                            next_query.last - current_query.first,
                            current_target.first - next_target.last,
                        )
                    };

                    // Check if overlaps are proportional (same alignment)
                    if query_overlap_len > 0 && target_overlap_len > 0 {
                        // Check if CIGAR strings are identical in the overlap region
                        let overlap_matches = check_cigar_overlap_match(
                            &current_cigar,
                            &next_cigar,
                            query_overlap_len,
                            query_forward,
                        );

                        if overlap_matches {
                            // debug!(
                            //     "Merge overlapping! Overlap: query={}, target={}, Query: current {}:{}-{}({}), next {}:{}-{}({}); Target: current {}:{}-{}({}), next {}:{}-{}({})",
                            //     query_overlap_len,
                            //     target_overlap_len,
                            //     current_query.metadata,
                            //     current_query.first,
                            //     current_query.last,
                            //     if query_forward { "+" } else { "-" },
                            //     next_query.metadata,
                            //     next_query.first,
                            //     next_query.last,
                            //     if next_query_forward { "+" } else { "-" },
                            //     current_target.metadata,
                            //     current_target.first,
                            //     current_target.last,
                            //     if target_forward { "+" } else { "-" },
                            //     next_target.metadata,
                            //     next_target.first,
                            //     next_target.last,
                            //     if next_target_forward { "+" } else { "-" },
                            // );

                            // Trim the overlap from the next interval and merge
                            let trimmed_next_cigar = trim_cigar_prefix(
                                &next_cigar,
                                query_overlap_len,
                                target_overlap_len,
                            );

                            if query_forward {
                                current_query.last = next_query.last;
                                current_target.last = next_target.last;
                                current_cigar.extend(trimmed_next_cigar);
                            } else {
                                current_query.first = next_query.first;
                                current_target.first = next_target.first;

                                let mut new_cigar = Vec::with_capacity(
                                    trimmed_next_cigar.len() + current_cigar.len(),
                                );
                                new_cigar.extend(trimmed_next_cigar);
                                new_cigar.extend_from_slice(&current_cigar);
                                current_cigar = new_cigar;
                            }
                            continue;
                        }
                    }
                }

                // Handle gaps within merge distance
                if !query_overlap && !target_overlap {
                    let (query_gap, target_gap) = if query_forward {
                        (
                            next_query.first - current_query.last,
                            next_target.first - current_target.last,
                        )
                    } else {
                        (
                            current_query.first - next_query.last,
                            current_target.first - next_target.last,
                        )
                    };

                    // Check if gaps are within merge distance and at least one gap exists
                    if query_gap >= 0
                        && target_gap >= 0
                        && (query_gap > 0 || target_gap > 0)
                        && query_gap <= merge_distance
                        && target_gap <= merge_distance
                    {
                        // debug!(
                        //     "Merge gaps! Query gap: {}, Target gap: {}, Query: current {}:{}-{}({}), next {}:{}-{}({}); Target: current {}:{}-{}({}), next {}:{}-{}({})",
                        //     query_gap,
                        //     target_gap,
                        //     current_query.metadata,
                        //     current_query.first,
                        //     current_query.last,
                        //     if query_forward { "+" } else { "-" },
                        //     next_query.metadata,
                        //     next_query.first,
                        //     next_query.last,
                        //     if next_query_forward { "+" } else { "-" },
                        //     current_target.metadata,
                        //     current_target.first,
                        //     current_target.last,
                        //     if target_forward { "+" } else { "-" },
                        //     next_target.metadata,
                        //     next_target.first,
                        //     next_target.last,
                        //     if next_target_forward { "+" } else { "-" },
                        // );

                        // Create gap-filling CIGAR operations
                        let mut gap_cigar = Vec::new();

                        if query_gap > 0 {
                            gap_cigar.push(CigarOp::new(query_gap, 'I'));
                        }
                        if target_gap > 0 {
                            gap_cigar.push(CigarOp::new(target_gap, 'D'));
                        }

                        // Merge intervals and CIGAR
                        if query_forward {
                            current_query.last = next_query.last;
                            current_target.last = next_target.last;
                            current_cigar.extend(gap_cigar);
                            current_cigar.extend_from_slice(&next_cigar);
                        } else {
                            current_query.first = next_query.first;
                            current_target.first = next_target.first;

                            let mut new_cigar = Vec::with_capacity(
                                current_cigar.len() + gap_cigar.len() + next_cigar.len(),
                            );
                            new_cigar.extend_from_slice(&next_cigar);
                            new_cigar.extend(gap_cigar);
                            new_cigar.extend_from_slice(&current_cigar);
                            current_cigar = new_cigar;
                        }
                        merge_consecutive_cigar_ops(&mut current_cigar);
                        continue;
                    }
                }

                // No merge possible - store current and move to next
                merged_results.push((current_query, current_cigar, current_target));
                (current_query, current_cigar, current_target) =
                    (next_query, next_cigar, next_target);
            }

            // Don't forget to add the last current element
            merged_results.push((current_query, current_cigar, current_target));

            // Replace original results with merged results
            *results = merged_results;

            info!("Collected {} merged intervals", results.len());
        }
    }
}

/// Gap-tolerant 2D merge for adjusted intervals.
///
/// Groups by `(query_id, target_id, query_strand)` and runs a forward-progress
/// union-find: two intervals i,j with `i.q_start < j.q_start` merge iff
/// `q_gap ≤ d AND t_gap ≤ d` on the strand-appropriate axes, with overlap
/// allowed (negative gaps pass). Union-find makes merging transitive.
///
/// Safe on empty CIGARs (syng) — merged intervals keep an empty CIGAR.
/// Non-empty CIGARs are concatenated in query order; the concatenation is
/// only semantically meaningful when the merged chunks were contiguous, so
/// PAF output should keep using `merge_adjusted_intervals`.
fn merge_adjusted_intervals_gap_2d(
    results: &mut Vec<AdjustedInterval>,
    merge_distance: i32,
) {
    if results.len() <= 1 || merge_distance < 0 {
        return;
    }
    let d = merge_distance as i64;
    use rustc_hash::FxHashMap as Map;

    let mut groups: Map<(u32, u32, bool), Vec<usize>> = Map::default();
    for (i, (q, _, t)) in results.iter().enumerate() {
        let strand_fwd = q.first <= q.last;
        groups
            .entry((q.metadata, t.metadata, strand_fwd))
            .or_default()
            .push(i);
    }

    let n = results.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn uf_find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    for ((_, _, strand_fwd), mut indices) in groups {
        indices.sort_by_key(|&i| {
            let q = &results[i].0;
            if strand_fwd { q.first } else { -q.first }
        });
        for a_pos in 0..indices.len() {
            let ia = indices[a_pos];
            let qa = results[ia].0;
            let ta = results[ia].2;
            let (qa_start, qa_end) = if strand_fwd {
                (qa.first as i64, qa.last as i64)
            } else {
                (qa.last as i64, qa.first as i64)
            };
            let (ta_start, ta_end) = (ta.first as i64, ta.last as i64);
            for b_pos in (a_pos + 1)..indices.len() {
                let ib = indices[b_pos];
                let qb = results[ib].0;
                let tb = results[ib].2;
                let qb_start = if strand_fwd { qb.first as i64 } else { qb.last as i64 };
                // Reject strictly-backward starts only. Equal starts are fine
                // (syng: every row shares the full query tile). Forward progress
                // on the target axis below does the real work of ordering.
                if qb_start < qa_start {
                    continue;
                }
                let q_gap = qb_start - qa_end;
                if q_gap > d {
                    break;
                }
                let (tb_start, tb_end) = (tb.first as i64, tb.last as i64);
                let (t_gap, t_forward) = if strand_fwd {
                    (tb_start - ta_end, tb_start > ta_start)
                } else {
                    (ta_start - tb_end, tb_end < ta_end)
                };
                if !t_forward || t_gap > d {
                    continue;
                }
                let ra = uf_find(&mut parent, ia);
                let rb = uf_find(&mut parent, ib);
                if ra != rb {
                    parent[ra] = rb;
                }
            }
        }
    }

    let mut buckets: Map<usize, Vec<usize>> = Map::default();
    for i in 0..n {
        let r = uf_find(&mut parent, i);
        buckets.entry(r).or_default().push(i);
    }

    let mut merged: Vec<AdjustedInterval> = Vec::with_capacity(buckets.len());
    let mut taken = vec![false; n];
    for i in 0..n {
        if taken[i] {
            continue;
        }
        let r = uf_find(&mut parent, i);
        let Some(members) = buckets.remove(&r) else { continue };
        for &m in &members {
            taken[m] = true;
        }
        let strand_fwd = {
            let q = &results[members[0]].0;
            q.first <= q.last
        };
        let mut ordered = members.clone();
        ordered.sort_by_key(|&idx| {
            let q = &results[idx].0;
            if strand_fwd { q.first } else { -q.first }
        });
        let first = &results[ordered[0]];
        let (mut q_lo, mut q_hi) = (first.0.first, first.0.last);
        let (mut t_lo, mut t_hi) = (first.2.first, first.2.last);
        let q_meta = first.0.metadata;
        let t_meta = first.2.metadata;
        let mut cigar: Vec<CigarOp> = Vec::new();
        for &idx in &ordered {
            let (q, c, t) = &results[idx];
            if strand_fwd {
                q_lo = q_lo.min(q.first);
                q_hi = q_hi.max(q.last);
            } else {
                q_lo = q_lo.max(q.first);
                q_hi = q_hi.min(q.last);
            }
            t_lo = t_lo.min(t.first);
            t_hi = t_hi.max(t.last);
            cigar.extend_from_slice(c);
        }
        merge_consecutive_cigar_ops(&mut cigar);
        let q_iv = coitrees::Interval { first: q_lo, last: q_hi, metadata: q_meta };
        let t_iv = coitrees::Interval { first: t_lo, last: t_hi, metadata: t_meta };
        merged.push((q_iv, cigar, t_iv));
    }

    let count = merged.len();
    *results = merged;
    info!("Collected {} merged intervals (gap-2d, d={})", count, merge_distance);
}

// Merge consecutive operations of the same type
fn merge_consecutive_cigar_ops(cigar: &mut Vec<CigarOp>) {
    if cigar.len() <= 1 {
        return;
    }

    let mut write_idx = 0;
    for read_idx in 1..cigar.len() {
        if cigar[write_idx].op() == cigar[read_idx].op() {
            // Same operation type - merge by adding lengths
            let combined_len = cigar[write_idx].len() + cigar[read_idx].len();
            cigar[write_idx] = CigarOp::new(combined_len, cigar[write_idx].op());
        } else {
            // Different operation types - keep separate
            write_idx += 1;
            if write_idx != read_idx {
                cigar[write_idx] = cigar[read_idx].clone();
            }
        }
    }
    cigar.truncate(write_idx + 1);
}

// Check if CIGAR strings match in the overlap region
fn check_cigar_overlap_match(
    current_cigar: &[CigarOp],
    next_cigar: &[CigarOp],
    query_overlap_len: i32,
    query_forward: bool,
) -> bool {
    // Extract the suffix of current CIGAR that corresponds to the overlap
    let current_suffix = extract_cigar_suffix(current_cigar, query_overlap_len, query_forward);

    // Extract the prefix of next CIGAR that corresponds to the overlap
    let next_prefix = extract_cigar_prefix(next_cigar, query_overlap_len, query_forward);

    // Compare if they're identical
    current_suffix == next_prefix
}

// Extract the last part of CIGAR that covers query_len bases
fn extract_cigar_suffix(cigar: &[CigarOp], query_len: i32, forward: bool) -> Vec<CigarOp> {
    let mut result = Vec::new();
    let mut remaining_query = query_len;

    // Traverse CIGAR from end to beginning
    for op in cigar.iter().rev() {
        if remaining_query <= 0 {
            break;
        }

        let query_delta = op
            .query_delta(if forward {
                Strand::Forward
            } else {
                Strand::Reverse
            })
            .abs();

        if query_delta <= remaining_query {
            // Include entire operation
            result.push(op.clone());
            remaining_query -= query_delta;
        } else if query_delta > 0 {
            // Include partial operation
            let scale = remaining_query as f32 / query_delta as f32;
            let new_len = (op.len() as f32 * scale) as i32;
            let partial_op = CigarOp::new(new_len, op.op());
            result.push(partial_op);
            remaining_query = 0;
        }
    }

    // Reverse to get correct order
    result.reverse();
    result
}

// Extract the first part of CIGAR that covers query_len bases
fn extract_cigar_prefix(cigar: &[CigarOp], query_len: i32, forward: bool) -> Vec<CigarOp> {
    let mut result = Vec::new();
    let mut remaining_query = query_len;

    for op in cigar.iter() {
        if remaining_query <= 0 {
            break;
        }

        let query_delta = op
            .query_delta(if forward {
                Strand::Forward
            } else {
                Strand::Reverse
            })
            .abs();

        if query_delta <= remaining_query {
            // Include entire operation
            result.push(op.clone());
            remaining_query -= query_delta;
        } else if query_delta > 0 {
            // Include partial operation
            let scale = remaining_query as f32 / query_delta as f32;
            let new_len = (op.len() as f32 * scale) as i32;
            let partial_op = CigarOp::new(new_len, op.op());
            result.push(partial_op);
            remaining_query = 0;
        }
    }

    result
}

// Trim the prefix of CIGAR by removing operations that cover the first query_len/target_len bases
fn trim_cigar_prefix(cigar: &[CigarOp], query_len: i32, target_len: i32) -> Vec<CigarOp> {
    let mut result = Vec::new();
    let mut query_consumed = 0;
    let mut target_consumed = 0;
    let mut start_idx = 0;

    // Find where to start after trimming
    for (idx, op) in cigar.iter().enumerate() {
        let q_delta = op.query_delta(Strand::Forward).abs();
        let t_delta = op.target_delta();

        if query_consumed + q_delta > query_len || target_consumed + t_delta > target_len {
            // This operation partially overlaps - need to trim it
            let query_remaining = query_len - query_consumed;
            let target_remaining = target_len - target_consumed;

            // Calculate how much of this operation to skip
            let skip_ratio = if q_delta > 0 && t_delta > 0 {
                (query_remaining as f32 / q_delta as f32)
                    .min(target_remaining as f32 / t_delta as f32)
            } else if q_delta > 0 {
                query_remaining as f32 / q_delta as f32
            } else if t_delta > 0 {
                target_remaining as f32 / t_delta as f32
            } else {
                0.0
            };

            let skip_len = (op.len() as f32 * skip_ratio) as i32;

            if skip_len < op.len() {
                // Create partial operation with remaining length
                let partial_op = CigarOp::new(op.len() - skip_len, op.op());
                result.push(partial_op);
            }

            // Add all remaining operations
            start_idx = idx + 1;
            break;
        }

        query_consumed += q_delta;
        target_consumed += t_delta;

        if query_consumed >= query_len && target_consumed >= target_len {
            start_idx = idx + 1;
            break;
        }
    }

    // Add all remaining operations
    result.extend_from_slice(&cigar[start_idx..]);
    result
}

fn print_list_sequences(impg: &impl ImpgIndex) {
    let seq_index = impg.seq_index();
    let num_sequences = seq_index.len();

    println!("Sequence\tLength");
    for i in 0..num_sequences as u32 {
        if let (Some(name), Some(len)) = (seq_index.get_name(i), seq_index.get_len_from_id(i)) {
            println!("{}\t{}", name, len);
        }
    }
}

fn print_stats(impg: &impl ImpgIndex) {
    // Basic stats
    let num_sequences = impg.seq_index().len();
    let total_sequence_length: usize = (0..num_sequences as u32)
        .into_par_iter()
        .filter_map(|id| impg.seq_index().get_len_from_id(id))
        .sum();

    // Compute overlap stats
    let num_targets = impg.num_targets();

    info!("Computing statistics for {} trees...", num_targets);

    // Collect target IDs to process
    let target_ids: Vec<u32> = impg.target_ids();

    // Process trees in parallel, computing stats without keeping them in memory
    let results: Vec<(u32, usize)> = target_ids
        .par_iter()
        .map(|&target_id| {
            // Get tree and count, removing it from memory for efficiency
            let count = if let Some(tree) = impg.get_or_load_tree(target_id) {
                let len = tree.len();
                impg.remove_cached_tree(target_id);
                len
            } else {
                0
            };

            (target_id, count)
        })
        .collect();

    // Collect results
    let mut num_overlaps = 0;
    let mut overlaps_per_seq = FxHashMap::default();

    for (target_id, count) in results {
        num_overlaps += count;
        overlaps_per_seq.insert(target_id, count);
    }

    info!(
        "Processed {} trees, total overlaps: {}",
        target_ids.len(),
        num_overlaps
    );

    println!("Number of query+target sequences: {num_sequences}");
    println!("Total query+target sequence length: {total_sequence_length} bp");
    println!("Number of overlaps: {num_overlaps}");

    let mut entries: Vec<(u32, usize)> = overlaps_per_seq.into_iter().collect();
    if !entries.is_empty() {
        entries.par_sort_by(|a, b| b.1.cmp(&a.1));

        // Calculate mean and median overlaps
        let sum: usize = entries.par_iter().map(|(_, count)| count).sum();
        let mean = sum as f64 / entries.len() as f64;

        let median = if entries.is_empty() {
            0.0
        } else if entries.len().is_multiple_of(2) {
            let mid = entries.len() / 2;
            (entries[mid - 1].1 + entries[mid].1) as f64 / 2.0
        } else {
            entries[entries.len() / 2].1 as f64
        };
        println!("\nMean overlaps per sequence: {mean:.2}");
        println!("Median overlaps per sequence: {median:.2}");

        println!("\nTop target sequences by number of overlaps:");
        for (idx, (seq_id, count)) in entries.iter().take(5).enumerate() {
            if let Some(name) = impg.seq_index().get_name(*seq_id) {
                println!("{}. {}: {} overlaps", idx + 1, name, count);
            }
        }
    }

    // Bridge genome statistics
    info!("Computing bridge genome statistics...");

    use std::collections::HashSet;
    let target_ids_set: HashSet<u32> = impg.target_ids().into_iter().collect();

    // Collect query IDs from all trees
    let query_ids_set: HashSet<u32> = impg
        .target_ids()
        .par_iter()
        .flat_map(|&target_id| {
            if let Some(tree) = impg.get_or_load_tree(target_id) {
                let ids: Vec<u32> = tree
                    .iter()
                    .map(|interval| interval.metadata.query_id())
                    .collect();
                impg.remove_cached_tree(target_id);
                ids
            } else {
                Vec::new()
            }
        })
        .collect();

    let bridges: HashSet<u32> = target_ids_set
        .intersection(&query_ids_set)
        .copied()
        .collect();
    let target_only = target_ids_set.len() - bridges.len();
    let query_only = query_ids_set.len() - bridges.len();

    println!("\nBridge genome coverage:");
    println!("  Total sequences: {}", num_sequences);
    println!("  Target sequences: {}", target_ids_set.len());
    println!("  Query sequences: {}", query_ids_set.len());
    println!(
        "  Bridge sequences (both): {} ({:.1}%)",
        bridges.len(),
        100.0 * bridges.len() as f64 / num_sequences as f64
    );
    println!("  Target-only: {}", target_only);
    println!("  Query-only: {}", query_only);

    if bridges.len() < target_ids_set.len() / 2 {
        warn!(
            "Low bridge coverage ({:.1}%) may limit transitive query reach. \
             Rebuild without --unidirectional for better coverage (bidirectional mode is default).",
            100.0 * bridges.len() as f64 / target_ids_set.len() as f64
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_subsequence_coordinates() {
        // Test parsing subsequence coordinates
        let result = parse_subsequence_coordinates("HG002#1#chr1:5116130-6116563");
        assert_eq!(result, Some(("HG002#1#chr1".to_string(), 5116130)));

        let result = parse_subsequence_coordinates("GRCh38#0#chr1:5477602-6474357");
        assert_eq!(result, Some(("GRCh38#0#chr1".to_string(), 5477602)));

        // Test with no subsequence coordinates
        let result = parse_subsequence_coordinates("chr1");
        assert_eq!(result, None);

        // Test with invalid format
        let result = parse_subsequence_coordinates("chr1:invalid");
        assert_eq!(result, None);
    }

    #[test]
    fn test_resolve_syng_syncmer_params_defaults_match_paper() {
        let params = resolve_syng_syncmer_params(None, None, None, None, 7).unwrap();
        assert_eq!(params.k, 8);
        assert_eq!(params.w, 55);
        assert_eq!(params.k + params.w, 63);
        assert_eq!(params.seed, 7);
    }

    #[test]
    fn test_resolve_syng_syncmer_params_paper_names() {
        let params = resolve_syng_syncmer_params(None, Some(8), None, Some(63), 7).unwrap();
        assert_eq!(params.k, 8);
        assert_eq!(params.w, 55);
        assert_eq!(params.k + params.w, 63);
    }

    #[test]
    fn test_resolve_syng_syncmer_params_legacy_names() {
        let params = resolve_syng_syncmer_params(Some(8), None, Some(55), None, 7).unwrap();
        assert_eq!(params.k, 8);
        assert_eq!(params.w, 55);
        assert_eq!(params.k + params.w, 63);
    }

    #[test]
    fn test_resolve_syng_syncmer_params_rejects_even_total_length() {
        let err = resolve_syng_syncmer_params(None, Some(8), None, Some(64), 7).unwrap_err();
        assert!(err.to_string().contains("must be odd"));
    }

    #[test]
    fn test_transform_coordinates_to_original() {
        // Test with original coordinates enabled
        let (name, start, end) =
            transform_coordinates_to_original("HG002#1#chr1:5116130-6116563", 45803, 45861, true);
        assert_eq!(name, "HG002#1#chr1");
        assert_eq!(start, 5116130 + 45803);
        assert_eq!(end, 5116130 + 45861);

        // Test with original coordinates disabled
        let (name, start, end) =
            transform_coordinates_to_original("HG002#1#chr1:5116130-6116563", 45803, 45861, false);
        assert_eq!(name, "HG002#1#chr1:5116130-6116563");
        assert_eq!(start, 45803);
        assert_eq!(end, 45861);

        // Test with sequence name that doesn't contain subsequence coordinates
        let (name, start, end) = transform_coordinates_to_original("chr1", 45803, 45861, true);
        assert_eq!(name, "chr1");
        assert_eq!(start, 45803);
        assert_eq!(end, 45861);
    }
}
