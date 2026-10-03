//! THE ANCHOR-PROJECTION LAYER (slice A of the local-realignment
//! evidence model; owner architecture ruling 2026-11-06;
//! assessment-side only — no product file is touched, no instrument of
//! record is modified, no threshold enters anything).
//!
//! THE RULING: the local-realignment evidence model EXTENDS the
//! existing syng GBWT MEM-to-read projection — it does NOT replace it.
//! The MEMs are the anchor layer: they place each read exactly on
//! interned nodes (that machinery stands). The realignment projects
//! THE BASES THE MEMs SKIP — the read positions outside the verified
//! anchor windows — through the graph along the anchored placement,
//! against each candidate path's spelled sequence at the projected
//! position. NO new aligner, NO parallel machinery — the existing MEM
//! mapping is projected through.
//!
//! THIS SLICE builds and receipts ONLY the anchor-projection layer:
//!   * PER-READ ANCHORED PLACEMENTS: every committed multi-census
//!     occurrence (the read-matched stage's verified occurrence set)
//!     is a placement of the record's canonical anchor chain at
//!     (path, start, orientation); each read carrying the record
//!     inherits the placement with its own PHYSICAL strand
//!     (the census orientation refers to the canonical token walk,
//!     which may be that read's own walk or its lexicographic mirror,
//!     so the physical strand is (orientation == 0) != mirrored) and
//!     its own skipped bases (the read positions outside the merged
//!     anchor windows — the flanks; the committed census measured
//!     zero interior gaps over all occurrences, and the merge is
//!     computed honestly here, gaps included).
//!   * THE PROJECTION EQUATIONS (own read position i, pattern span S,
//!     hull start w_lo, window length k): physically forward places
//!     read position i at path position start + (i - w_lo); physically
//!     reverse places it at start + S - 1 - (i - w_lo). Both are
//!     asserted at every anchor window against the panel for every
//!     occurrence, and the read hull is asserted equal to the placed
//!     segment (flanks are NOT asserted — they are the skipped
//!     evidence the scoring slice votes).
//!   * PER-PATH GRAPH CONTEXT: for every candidate row (the member
//!     rows of the axis partition of every window the occurrence
//!     touches — the receipts' own locality) whose walk shares at
//!     least one anchor node with the placement's required sign, the
//!     anchor correspondence (the row-side window-start position of
//!     each anchor; repeated nodes resolved by the unique monotone
//!     assignment rule, otherwise named ambiguous) and, per skipped
//!     base, the row-side coordinate, the covering segment(s) (the
//!     signed node ids whose windows cover it; none = an inter-window
//!     gap pocket) or the out-of-extent flag, and the row's base —
//!     the lookups the scoring slice consumes. The scoring itself is
//!     NOT built here.
//!   * RECEIPTS: `--out` carries the full placement layer (per
//!     record: the canonical anchor chain, per read variant the
//!     mirror flag + flank bases + instance count, per occurrence the
//!     census-verbatim placement + the re-derived merged read-matched
//!     interval count, cross-validated against the committed census
//!     covered-node lists for EVERY occurrence); `--context-out`
//!     carries the per-path graph context for an evenly-spaced stated
//!     sample of occurrences (the full-fidelity per-path context over
//!     every occurrence x every candidate row measures in the tens of
//!     GB on chrI and is NOT emitted; the scoring slice computes
//!     contexts in-process with this same committed machinery).
//!
//! Inputs: the panel syng, the route graph's sources (panel sequence
//! fetch), the partition-graph maps (the committed partition-embedded
//! build), the committed multi-matching census receipt, and the sample
//! FASTQ (read identity). The record derivation cache is a pure
//! function of the FASTQ and the panel (component-independent).
//!
//! Usage (assessment-side):
//!   partition_anchor_projection --panel <syng-prefix> --routes <dir> \
//!     --partition-graphs <dir> --component S288C#0#chrMT \
//!     --census <cosine-multi-census.jsonl> \
//!     --reads <reads.fastq.gz> --derive-cache <cache.bin> \
//!     --out <anchor-projection.jsonl> --context-out <context.jsonl>

#![recursion_limit = "512"]

use clap::Parser;
use impg::genome_inference::{panel_routes as routes, read_json, PanelIdentity};
use impg::syng::{SyncmerParams, SyngIndex};
use rayon::prelude::*;
use serde::Deserialize;
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::time::Instant;

const READ_LENGTH: usize = 150;

// ---------------------------------------------------------------------------
// Options and small utilities.
// ---------------------------------------------------------------------------

#[derive(Parser)]
struct Options {
    /// The panel syng prefix.
    #[arg(long)]
    panel: String,
    /// The route graph directory (graph.json + sources).
    #[arg(long)]
    routes: PathBuf,
    /// The partition-embedded graph directory (partition<N>.gfa.map.json).
    #[arg(long)]
    partition_graphs: PathBuf,
    /// The component lane (e.g. S288C#0#chrMT) — the axis path whose
    /// member rows define the window-to-partition mapping.
    #[arg(long)]
    component: String,
    /// The committed multi-matching census receipt (records, placements).
    #[arg(long)]
    census: PathBuf,
    /// The sample FASTQ (read identity).
    #[arg(long)]
    reads: PathBuf,
    /// The record-derivation cache (written when absent, loaded when
    /// present) — a pure function of the FASTQ and the panel.
    #[arg(long)]
    derive_cache: PathBuf,
    /// The anchored-projection receipt path.
    #[arg(long)]
    out: PathBuf,
    /// The per-path graph-context receipt path (the stated sample).
    #[arg(long)]
    context_out: PathBuf,
    /// The evenly-spaced occurrence sample count for the per-path
    /// context receipt (a stated measurement sample; not a threshold —
    /// nothing downstream of the receipts depends on which occurrences
    /// are sampled).
    #[arg(long, default_value_t = 1000)]
    sample_count: u64,
    /// Resident-set guard in GiB (the 64 GiB discipline; 0 = no guard).
    #[arg(long, default_value_t = 64.0)]
    rss_budget_gib: f64,
}

fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, s.to_string())
}

fn ensure(ok: bool, s: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(s))
    }
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

struct RssGuard {
    budget_bytes: Option<u64>,
}

impl RssGuard {
    fn new(budget_gib: f64) -> Self {
        RssGuard {
            budget_bytes: (budget_gib > 0.0).then(|| (budget_gib * (1u64 << 30) as f64) as u64),
        }
    }
    fn probe(&self, stage: &str) -> io::Result<u64> {
        let text = std::fs::read_to_string("/proc/self/status")?;
        let rss_kb = text
            .lines()
            .find(|l| l.starts_with("VmRSS:"))
            .and_then(|l| l.split_whitespace().nth(1).and_then(|n| n.parse::<u64>().ok()))
            .unwrap_or(0);
        if let Some(budget) = self.budget_bytes {
            ensure(
                rss_kb * 1024 <= budget,
                &format!("RSS guard exceeded at {stage}: {rss_kb} kB"),
            )?;
        }
        eprintln!("[anchor] rss {stage}: {rss_kb} kB");
        Ok(rss_kb)
    }
}

// ---------------------------------------------------------------------------
// Receipt inputs.
// ---------------------------------------------------------------------------

/// One verified occurrence in the committed multi-matching census.
#[derive(Deserialize)]
struct CensusOccurrence {
    path: usize,
    start: u64,
    orientation: u8,
    /// The touched component windows (the routing's own territory-touch
    /// rule; these are WINDOW ids, not partition ids).
    #[serde(default)]
    partitions: Vec<u32>,
    /// Per read-matched interval, the contained covered node abs ids in
    /// path order (the committed containment rule).
    #[serde(default)]
    intervals: Vec<Vec<u32>>,
}

#[derive(Deserialize)]
struct CensusRecordLine {
    record: u64,
    multiplicity: u64,
    #[serde(default)]
    anchors: u64,
    #[serde(default)]
    occurrences: Vec<CensusOccurrence>,
}

/// One member row of a partition graph map (the alignment-induced
/// membership; the partition BEDs are source-forward).
#[derive(Clone, Deserialize)]
struct MemberRow {
    path_name: String,
    start: u64,
    end: u64,
}

#[derive(Clone, Deserialize)]
struct PartitionMap {
    partition: u32,
    members: Vec<MemberRow>,
}

// ---------------------------------------------------------------------------
// The projection core (pure functions; unit-tested).
// ---------------------------------------------------------------------------

/// The canonical decoded walk of a record's token pattern: anchors
/// (signed node, canonical position), first anchor at position 0.
fn decode_record_tokens(tokens: &[u64]) -> io::Result<Vec<(i32, u64)>> {
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

/// The record's span: the last anchor's window end in canonical
/// coordinates (the pattern's full hull length).
fn walk_span(walk: &[(i32, u64)], k: u64) -> u64 {
    walk.last().map(|&(_, p)| p + k).unwrap_or(k)
}

/// The read's OWN anchor positions for one mirror state, per canonical
/// anchor index: an unmirrored read spells the canonical pattern (own
/// anchor j at w_lo + c_j); a mirrored read spells its reverse
/// complement (canonical anchor j at w_lo + (c_last - c_j)).
fn own_positions_of(canonical: &[(i32, u64)], mirror: bool, w_lo: u64) -> Vec<u64> {
    let last = canonical.last().map(|&(_, p)| p).unwrap_or(0);
    canonical
        .iter()
        .map(|&(_, pos)| {
            if mirror {
                w_lo + (last - pos)
            } else {
                w_lo + pos
            }
        })
        .collect()
}

/// The read's physical strand at one occurrence: the census orientation
/// refers to the canonical token walk, which is the read's own walk
/// exactly when the read is not mirrored.
fn physical_forward(orientation: u8, mirror: bool) -> bool {
    (orientation == 0) != mirror
}

/// The panel-side path coordinate of read own position i under one
/// placement (start = the pattern's position-0 window start on the
/// path; span = the pattern hull length; w_lo = the hull start in own
/// read coordinates):
///   physically forward:  start + (i - w_lo)
///   physically reverse:  start + span - 1 - (i - w_lo)
#[allow(dead_code)] // the tested core equation; the run paths use the
                    // anchor-relative forms derived from it.
fn project_own_position(i: u64, start: u64, span: u64, w_lo: u64, forward: bool) -> i64 {
    if forward {
        start as i64 + (i as i64 - w_lo as i64)
    } else {
        start as i64 + span as i64 - 1 - (i as i64 - w_lo as i64)
    }
}

/// The panel-side window start of canonical anchor j under one
/// placement (the projection equation at the anchor windows):
///   orientation 0:  start + c_j
///   orientation 1:  start + span - k - c_j
fn anchor_window_start(canonical_pos: u64, start: u64, span: u64, k: u64, orientation: u8) -> u64 {
    if orientation == 0 {
        start + canonical_pos
    } else {
        start + span - k - canonical_pos
    }
}

/// Merge windows (lo, hi): overlapping and exactly-abutting windows
/// merge; a positive gap splits. Returns the merged intervals.
fn merge_windows(mut windows: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    windows.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(windows.len());
    for (lo, hi) in windows {
        match merged.last_mut() {
            Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
            _ => merged.push((lo, hi)),
        }
    }
    merged
}

/// The complement of the merged own-coordinate anchor windows within
/// the read: the SKIPPED base positions (the bases the MEMs do not
/// cover).
fn skipped_positions(own_positions: &[u64], k: u64, read_len: u64) -> Vec<(u64, u64)> {
    let windows: Vec<(u64, u64)> = own_positions.iter().map(|&p| (p, p + k)).collect();
    let merged = merge_windows(windows);
    let mut skipped: Vec<(u64, u64)> = Vec::new();
    let mut cursor = 0u64;
    for (lo, hi) in merged {
        if lo > cursor {
            skipped.push((cursor, lo));
        }
        cursor = cursor.max(hi);
    }
    if cursor < read_len {
        skipped.push((cursor, read_len));
    }
    skipped
}

/// The panel-side read-matched intervals of one occurrence: the merged
/// anchor windows on the path (the census's own recipe).
fn occurrence_intervals(canonical: &[(i32, u64)], k: u64, start: u64, orientation: u8) -> Vec<(u64, u64)> {
    let span = walk_span(canonical, k);
    let windows: Vec<(u64, u64)> = canonical
        .iter()
        .map(|&(_, pos)| {
            let lo = anchor_window_start(pos, start, span, k, orientation);
            (lo, lo + k)
        })
        .collect();
    merge_windows(windows)
}

/// The required panel-path sign of canonical anchor j at one
/// occurrence: orientation 0 spells the pattern (sign as decoded);
/// orientation 1 spells its reverse complement (flipped).
fn required_sign(node: i32, orientation: u8) -> i32 {
    if orientation == 0 {
        node
    } else {
        -node
    }
}

/// The anchor correspondence of one placement against one candidate
/// row's node positions (node abs id -> (window-start bp, sign)):
/// per canonical anchor, the row-side window-start position with the
/// required sign. A node occurring several times admits several
/// monotone assignments (strictly increasing for orientation 0,
/// decreasing for orientation 1, in canonical anchor order); a
/// position is PINNED iff it is identical across ALL monotone
/// assignments, otherwise it is AMBIGUOUS (None, named). Anchors
/// absent from the row are None but not ambiguous.
fn anchor_correspondence(
    canonical: &[(i32, u64)],
    orientation: u8,
    row_positions: &dyn Fn(u32) -> Vec<(u64, i32)>,
) -> (Vec<Option<u64>>, Vec<usize>) {
    let m = canonical.len();
    let mut candidates: Vec<Vec<u64>> = Vec::with_capacity(m);
    for &(node, _) in canonical {
        let required = required_sign(node, orientation);
        let positions = row_positions(required.unsigned_abs());
        candidates.push(
            positions
                .iter()
                .filter(|&&(_, sign)| sign == required.signum())
                .map(|&(bp, _)| bp)
                .collect(),
        );
    }
    let mut pinned: Vec<Option<u64>> = vec![None; m];
    let mut pinned_seen: Vec<bool> = vec![false; m];
    let mut assignments = 0u64;
    // The chain of (anchor index, row-side position) for the anchors
    // present in the row; absent anchors are skipped (they constrain
    // nothing, and the monotone rule applies to the surviving chain).
    let mut chain: Vec<(usize, u64)> = Vec::with_capacity(m);
    fn enumerate(
        index: usize,
        candidates: &[Vec<u64>],
        increasing: bool,
        chain: &mut Vec<(usize, u64)>,
        pinned: &mut Vec<Option<u64>>,
        pinned_seen: &mut Vec<bool>,
        assignments: &mut u64,
    ) {
        if index == candidates.len() {
            *assignments += 1;
            for &(j, bp) in chain.iter() {
                if !pinned_seen[j] {
                    pinned[j] = Some(bp);
                    pinned_seen[j] = true;
                } else if pinned[j] != Some(bp) {
                    pinned[j] = None;
                }
            }
            return;
        }
        if candidates[index].is_empty() {
            enumerate(index + 1, candidates, increasing, chain, pinned, pinned_seen, assignments);
            return;
        }
        for &bp in &candidates[index] {
            let ok = match chain.last() {
                Some(&(_, last)) => {
                    if increasing {
                        bp > last
                    } else {
                        bp < last
                    }
                }
                None => true,
            };
            if ok {
                chain.push((index, bp));
                enumerate(index + 1, candidates, increasing, chain, pinned, pinned_seen, assignments);
                chain.pop();
            }
        }
    }
    enumerate(
        0,
        &candidates,
        orientation == 0,
        &mut chain,
        &mut pinned,
        &mut pinned_seen,
        &mut assignments,
    );
    let ambiguous: Vec<usize> = (0..m)
        .filter(|&j| !candidates[j].is_empty() && pinned[j].is_none())
        .collect();
    let _ = assignments;
    (pinned, ambiguous)
}

// ---------------------------------------------------------------------------
// FASTQ streaming and the record-derivation cache (the cache is a pure
// function of the FASTQ and the panel; format identical to the
// prior realignment attempt's cache).
// ---------------------------------------------------------------------------

fn stream_fastq(
    path: &PathBuf,
    mut visit: impl FnMut(&[u8], &[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let (reader, _) = niffler::get_reader(Box::new(File::open(path)?)).map_err(io::Error::other)?;
    let mut lines = BufReader::new(reader).lines();
    let mut current = lines
        .next()
        .transpose()?
        .ok_or_else(|| invalid("empty reads file"))?;
    loop {
        ensure(current.starts_with('@'), "reads file is not FASTQ")?;
        let seq = lines
            .next()
            .transpose()?
            .ok_or_else(|| invalid("truncated FASTQ sequence"))?;
        let plus = lines
            .next()
            .transpose()?
            .ok_or_else(|| invalid("truncated FASTQ separator"))?;
        let qual = lines
            .next()
            .transpose()?
            .ok_or_else(|| invalid("truncated FASTQ quality"))?;
        ensure(
            plus.starts_with('+') && seq.len() == qual.len(),
            "invalid FASTQ record",
        )?;
        visit(seq.as_bytes(), qual.as_bytes())?;
        match lines.next().transpose()? {
            Some(line) => current = line,
            None => break,
        }
    }
    Ok(())
}

#[allow(clippy::type_complexity)]
fn write_derive_cache(
    path: &std::path::Path,
    key_tokens: &[Vec<u64>],
    key_reads: &[Vec<usize>],
    reads: &[Vec<u8>],
    read_multiplicity: &[u32],
    read_records: &[Vec<(u32, Vec<(i32, u32)>)>],
) -> io::Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    out.write_all(&(key_tokens.len() as u64).to_le_bytes())?;
    for tokens in key_tokens {
        out.write_all(&(tokens.len() as u32).to_le_bytes())?;
        for &token in tokens {
            out.write_all(&token.to_le_bytes())?;
        }
    }
    out.write_all(&(key_reads.len() as u64).to_le_bytes())?;
    for readers in key_reads {
        out.write_all(&(readers.len() as u32).to_le_bytes())?;
        for &index in readers {
            out.write_all(&(index as u32).to_le_bytes())?;
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
        for &(key, ref walk) in records {
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

#[allow(clippy::type_complexity)]
fn read_derive_cache(
    path: &std::path::Path,
) -> io::Result<(Vec<Vec<u64>>, Vec<Vec<usize>>, Vec<Vec<u8>>, Vec<u32>, Vec<Vec<(u32, Vec<(i32, u32)>)>>)> {
    let mut file = BufReader::new(File::open(path)?);
    let mut u64buf = [0u8; 8];
    let mut u32buf = [0u8; 4];
    let mut take_u64 = |file: &mut BufReader<File>| -> io::Result<u64> {
        file.read_exact(&mut u64buf)?;
        Ok(u64::from_le_bytes(u64buf))
    };
    let mut take_u32 = |file: &mut BufReader<File>| -> io::Result<u32> {
        file.read_exact(&mut u32buf)?;
        Ok(u32::from_le_bytes(u32buf))
    };
    let n_keys = take_u64(&mut file)?;
    let mut key_tokens = Vec::with_capacity(n_keys as usize);
    for _ in 0..n_keys {
        let len = take_u32(&mut file)? as usize;
        let mut tokens = Vec::with_capacity(len);
        for _ in 0..len {
            tokens.push(take_u64(&mut file)?);
        }
        key_tokens.push(tokens);
    }
    let n_lists = take_u64(&mut file)?;
    let mut key_reads = Vec::with_capacity(n_lists as usize);
    for _ in 0..n_lists {
        let len = take_u32(&mut file)? as usize;
        let mut readers = Vec::with_capacity(len);
        for _ in 0..len {
            readers.push(take_u32(&mut file)? as usize);
        }
        key_reads.push(readers);
    }
    let n_reads = take_u64(&mut file)?;
    let mut reads = Vec::with_capacity(n_reads as usize);
    for _ in 0..n_reads {
        let len = take_u32(&mut file)? as usize;
        let mut seq = vec![0u8; len];
        file.read_exact(&mut seq)?;
        reads.push(seq);
    }
    let mut read_multiplicity = Vec::with_capacity(n_reads as usize);
    for _ in 0..n_reads {
        read_multiplicity.push(take_u32(&mut file)?);
    }
    let n_record_lists = take_u64(&mut file)?;
    let mut read_records = Vec::with_capacity(n_record_lists as usize);
    for _ in 0..n_record_lists {
        let len = take_u32(&mut file)? as usize;
        let mut records = Vec::with_capacity(len);
        for _ in 0..len {
            let key = take_u32(&mut file)?;
            let walk_len = take_u32(&mut file)? as usize;
            let mut walk = Vec::with_capacity(walk_len);
            for _ in 0..walk_len {
                let mut local = [0u8; 4];
                file.read_exact(&mut local)?;
                let node = i32::from_le_bytes(local);
                let pos = take_u32(&mut file)?;
                walk.push((node, pos));
            }
            records.push((key, walk));
        }
        read_records.push(records);
    }
    Ok((key_tokens, key_reads, reads, read_multiplicity, read_records))
}

// ---------------------------------------------------------------------------
// Candidate rows (the partition graphs' member rows, deduplicated
// across partitions by (path, start, end)).
// ---------------------------------------------------------------------------

struct RowFold {
    path: usize,
    start: u64,
    end: u64,
    seq: Vec<u8>,
    /// All walk steps overlapping the row range (absolute bp, signed
    /// node), sorted by bp — the same walk the committed partition
    /// GFAs' P lines spell (the committed checker's spell-equality
    /// proof).
    walk: Vec<(u64, i32)>,
    node_positions: HashMap<u32, Vec<(u64, i32)>>,
    /// (partition, member index) pairs of every partition holding this
    /// row (rows are co-members of many partitions).
    members: Vec<(u32, usize)>,
}

impl RowFold {
    fn positions_of(&self, node_abs: u32) -> Vec<(u64, i32)> {
        self.node_positions
            .get(&node_abs)
            .cloned()
            .unwrap_or_default()
    }

    /// The covering signed node ids at one row-side path coordinate
    /// (windows [bp, bp + k) containing it; empty = an inter-window
    /// gap pocket).
    /// The index range [lo, hi) of the row's walk steps whose windows
    /// cover the coordinate (contiguous in the bp-sorted walk; empty =
    /// an inter-window gap pocket). The walk is the committed partition
    /// GFA's own step sequence, so the indices are exact segment
    /// references.
    fn covering_index_range(&self, coord: u64, k: u64) -> (usize, usize) {
        let lo = self.walk.partition_point(|&(bp, _)| bp + k <= coord);
        let hi = lo + self.walk[lo..].iter().take_while(|&&(bp, _)| bp <= coord).count();
        (lo, hi)
    }

    /// The row's own (path-forward) base at one row-side path
    /// coordinate; None when the coordinate is outside the row extent
    /// (out-of-extent: the partition graph does not spell it).
    fn base_at(&self, coord: i64) -> Option<u8> {
        if coord < self.start as i64 || coord >= self.end as i64 {
            return None;
        }
        self.seq.get((coord - self.start as i64) as usize).copied()
    }
}

/// The axis-partition row store: rows loaded lazily per partition,
/// deduplicated by (path, start, end) across partitions.
struct RowStore<'a> {
    panel: &'a SyngIndex,
    path_of_name: &'a HashMap<String, usize>,
    fetch: &'a dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>>,
    folds: Vec<RowFold>,
    dedup: HashMap<(usize, u64, u64), usize>,
    by_partition: HashMap<u32, Vec<usize>>,
    by_partition_path: HashMap<(u32, usize), Vec<usize>>,
}

impl<'a> RowStore<'a> {
    fn new(
        panel: &'a SyngIndex,
        path_of_name: &'a HashMap<String, usize>,
        fetch: &'a dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>>,
    ) -> Self {
        RowStore {
            panel,
            path_of_name,
            fetch,
            folds: Vec::new(),
            dedup: HashMap::new(),
            by_partition: HashMap::new(),
            by_partition_path: HashMap::new(),
        }
    }

    fn ensure_partition(&mut self, members: &[MemberRow], partition: u32) -> io::Result<()> {
        if self.by_partition.contains_key(&partition) {
            return Ok(());
        }
        let mut indices = Vec::with_capacity(members.len());
        for (member_index, member) in members.iter().enumerate() {
            let path_idx = *self
                .path_of_name
                .get(&member.path_name)
                .ok_or_else(|| invalid("member path absent from the syng"))?;
            let end = member.end.min(self.panel.name_map.path_to_length[path_idx]);
            ensure(member.start < end, "degenerate member row")?;
            let dedup_key = (path_idx, member.start, end);
            let fold_index = match self.dedup.get(&dedup_key) {
                Some(&index) => index,
                None => {
                    let seq = (self.fetch)(&member.path_name, member.start, end)?;
                    let mut walk: Vec<(u64, i32)> = self
                        .panel
                        .walk_path_range(path_idx, member.start, end)?
                        .into_iter()
                        .map(|(node, bp)| (bp, node))
                        .collect();
                    walk.sort_unstable_by_key(|&(bp, _)| bp);
                    let mut node_positions: HashMap<u32, Vec<(u64, i32)>> = HashMap::new();
                    for &(bp, node) in &walk {
                        node_positions
                            .entry(node.unsigned_abs())
                            .or_default()
                            .push((bp, node.signum()));
                    }
                    let index = self.folds.len();
                    self.folds.push(RowFold {
                        path: path_idx,
                        start: member.start,
                        end,
                        seq,
                        walk,
                        node_positions,
                        members: Vec::new(),
                    });
                    self.dedup.insert(dedup_key, index);
                    index
                }
            };
            self.folds[fold_index].members.push((partition, member_index));
            self.by_partition_path
                .entry((partition, path_idx))
                .or_default()
                .push(fold_index);
            indices.push(fold_index);
        }
        self.by_partition.insert(partition, indices);
        Ok(())
    }

    fn rows_of_partition_path(&self, partition: u32, path: usize) -> &[usize] {
        self.by_partition_path
            .get(&(partition, path))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

fn main() -> io::Result<()> {
    let started = Instant::now();
    let options = Options::parse();
    let rss = RssGuard::new(options.rss_budget_gib);

    // ------------------------------------------------------------- inputs
    let identity = PanelIdentity::read(&options.panel)?;
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    let k = panel.syncmer_length_bp() as u64;
    let graph: routes::Graph = read_json(&options.routes.join("graph.json"))?;
    if graph.panel != identity || !graph.generation_complete {
        return Err(invalid("incompatible route inventory"));
    }
    let lanes = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.clone(), lane.length))
        .collect::<Vec<_>>();
    let sources = routes::Sources::open(&graph.source_paths, lanes)?;
    let path_of_name: HashMap<String, usize> = panel
        .name_map
        .name_to_path
        .iter()
        .map(|(name, &path)| (name.clone(), path as usize))
        .collect();
    let source_of_name: HashMap<&str, usize> = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.as_str(), lane.id))
        .collect();
    let fetch_seq = |name: &str, lo: u64, hi: u64| -> io::Result<Vec<u8>> {
        let path_idx = *path_of_name
            .get(name)
            .ok_or_else(|| invalid("panel path absent from the syng name map"))?;
        let path_len = panel.name_map.path_to_length[path_idx];
        let hi = hi.min(path_len);
        let source = *source_of_name
            .get(name)
            .ok_or_else(|| invalid("panel path absent from the route sources"))?;
        if lo >= hi {
            return Ok(Vec::new());
        }
        sources.fetch(source, lo, hi)
    };

    // The census receipt.
    let mut census_records: Vec<CensusRecordLine> = Vec::new();
    {
        let (reader, _) = niffler::get_reader(Box::new(File::open(&options.census)?))
            .map_err(io::Error::other)?;
        for line in BufReader::new(reader).lines() {
            census_records.push(serde_json::from_str(&line?)?);
        }
    }
    ensure(
        census_records
            .windows(2)
            .all(|w| w[0].record + 1 == w[1].record),
        "census record ids are not dense",
    )?;

    // The partition graph maps; the AXIS partitions (holding a member
    // row on the component path) ranked by axis-row start map to the
    // component's window ids in order.
    let mut maps: BTreeMap<u32, PartitionMap> = BTreeMap::new();
    let mut axis_partitions: Vec<(u32, u64)> = Vec::new();
    for entry in std::fs::read_dir(&options.partition_graphs)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if !name.ends_with(".gfa.map.json") {
            continue;
        }
        let text = std::fs::read_to_string(&path)?;
        let map: PartitionMap = serde_json::from_str(&text)?;
        if let Some(axis) = map
            .members
            .iter()
            .find(|m| m.path_name == options.component)
        {
            axis_partitions.push((map.partition, axis.start));
        }
        maps.insert(map.partition, map);
    }
    axis_partitions.sort_by_key(|&(partition, start)| (start, partition));
    let n_windows = census_records
        .iter()
        .flat_map(|line| line.occurrences.iter())
        .flat_map(|occ| occ.partitions.iter().copied())
        .max()
        .map(|w| w as usize + 1)
        .unwrap_or(0);
    ensure(
        axis_partitions.len() == n_windows,
        "axis partition count does not match the census window span",
    )?;
    let window_partition: Vec<u32> =
        axis_partitions.iter().map(|&(partition, _)| partition).collect();
    eprintln!(
        "[anchor] inputs: {} census records, k {k}, {} windows, axis partitions {:?}",
        census_records.len(),
        n_windows,
        window_partition,
    );

    // ------------------------------------- the read-record re-derivation
    let derive_started = Instant::now();
    let (key_tokens, key_reads, reads, read_multiplicity, raw_read_records): (
        Vec<Vec<u64>>,
        Vec<Vec<usize>>,
        Vec<Vec<u8>>,
        Vec<u32>,
        Vec<Vec<(u32, Vec<(i32, u32)>)>>,
    ) = if options.derive_cache.exists() {
        eprintln!(
            "[anchor] loading the derivation cache {}",
            options.derive_cache.display()
        );
        let (tokens, lists, seqs, mults, records) = read_derive_cache(&options.derive_cache)?;
        ensure(
            seqs.iter().all(|seq| seq.len() == READ_LENGTH),
            "cached read is not L150",
        )?;
        (tokens, lists, seqs, mults, records)
    } else {
        let scheme = impg::genome_inference::mem_records::anchor_scheme();
        let mut key_index: BTreeMap<Vec<u64>, u32> = BTreeMap::new();
        let mut key_tokens: Vec<Vec<u64>> = Vec::new();
        let mut key_reads: Vec<Vec<usize>> = Vec::new();
        let mut reads: Vec<Vec<u8>> = Vec::new();
        let mut read_multiplicity: Vec<u32> = Vec::new();
        let mut read_records: Vec<Vec<(u32, Vec<(i32, u32)>)>> = Vec::new();
        {
            let mut read_counts: HashMap<u64, u32> = HashMap::new();
            let mut pending: Vec<Vec<u8>> = Vec::new();
            let mut drain = |pending: &mut Vec<Vec<u8>>| -> io::Result<()> {
                let per_read: Vec<Vec<(Vec<u64>, Vec<(i32, u32)>)>> = pending
                    .par_iter()
                    .map(|read| {
                        ensure(read.len() == READ_LENGTH, "sample read is not L150")?;
                        let tagged = impg::genome_inference::mem_records::tagged_mem_records_scheme(
                            &panel, read, scheme,
                        )?;
                        tagged
                            .iter()
                            .map(|walk| {
                                let encoded = impg::sample_mem_bwt::encode_walk(walk)?;
                                let tokens = impg::sample_mem_bwt::canonical(&encoded);
                                Ok((
                                    tokens,
                                    walk.iter()
                                        .map(|&(node, pos)| (node, pos as u32))
                                        .collect::<Vec<_>>(),
                                ))
                            })
                            .collect::<io::Result<Vec<_>>>()
                    })
                    .collect::<io::Result<_>>()?;
                for (index, seq) in pending.drain(..).enumerate() {
                    let base = reads.len();
                    reads.push(seq);
                    *read_counts.entry(fnv1a64(&reads[base])).or_default() += 1;
                    let mut own_records: Vec<(u32, Vec<(i32, u32)>)> =
                        Vec::with_capacity(per_read[index].len());
                    for (tokens, walk) in &per_read[index] {
                        let key = *key_index.entry(tokens.clone()).or_insert_with(|| {
                            key_tokens.push(tokens.clone());
                            key_reads.push(Vec::new());
                            key_tokens.len() as u32 - 1
                        });
                        key_reads[key as usize].push(base);
                        own_records.push((key, walk.clone()));
                    }
                    read_records.push(own_records);
                }
                Ok(())
            };
            stream_fastq(&options.reads, |seq, _qual| {
                pending.push(seq.to_vec());
                if pending.len() == 4096 {
                    drain(&mut pending)?;
                }
                Ok(())
            })?;
            if !pending.is_empty() {
                drain(&mut pending)?;
            }
            let mut first_seen: HashMap<u64, usize> = HashMap::new();
            for (index, seq) in reads.iter().enumerate() {
                let hash = fnv1a64(seq);
                match first_seen.get(&hash) {
                    Some(_) => read_multiplicity.push(0),
                    None => {
                        first_seen.insert(hash, index);
                        read_multiplicity.push(*read_counts.get(&hash).unwrap_or(&1));
                    }
                }
            }
        }
        write_derive_cache(
            &options.derive_cache,
            &key_tokens,
            &key_reads,
            &reads,
            &read_multiplicity,
            &read_records,
        )?;
        eprintln!(
            "[anchor] derivation cache written: {} ({:.1}s)",
            options.derive_cache.display(),
            derive_started.elapsed().as_secs_f64()
        );
        (key_tokens, key_reads, reads, read_multiplicity, read_records)
    };
    rss.probe("derivation")?;
    eprintln!(
        "[anchor] reads re-derived: {} reads, {} keys ({:.1}s)",
        reads.len(),
        key_tokens.len(),
        derive_started.elapsed().as_secs_f64(),
    );
    // The per-read record instances (own walks with absolute read
    // positions).
    struct ReadRecord {
        key: u32,
        walk: Vec<(i32, u32)>,
    }
    let read_records: Vec<Vec<ReadRecord>> = raw_read_records
        .into_iter()
        .map(|records| {
            records
                .into_iter()
                .map(|(key, walk)| ReadRecord { key, walk })
                .collect()
        })
        .collect();

    /// The routing's own canonical-scheme territory steps around one
    /// panel position (ported from build_territory_index_rows: the
    /// forward frame's steps are the syng's stored path walk; the rc
    /// frame's qualifying k-mers come from the raw extraction on the
    /// fetched range's reverse complement; per position, the frame that
    /// spells the k-mer's canonical form forward decides). The census
    /// occurrence starts live in THIS convention.
    let canonical_steps_near = |path: usize, around: u64| -> io::Result<Vec<(u64, i32)>> {
        let lo = around.saturating_sub(150);
        let hi = around + 150;
        let forward: Vec<(u64, i32)> = panel
            .walk_path_range(path, lo, hi)?
            .into_iter()
            .map(|(node, bp)| (bp, node))
            .collect();
        let seq_lo = lo.saturating_sub(k);
        let seq_hi = hi + k;
        let seq = fetch_seq(&panel.name_map.path_to_name[path], seq_lo, seq_hi)?;
        let rc_seq = revcomp(&seq);
        let reverse: Vec<(u64, i32)> =
            impg::genome_inference::mem_records::raw_matched_syncmers(&panel, &rc_seq)?
                .into_iter()
                .map(|(signed, q)| {
                    (
                        seq_lo + seq.len() as u64 - k - q,
                        -signed,
                    )
                })
                .collect();
        let forward_map: BTreeMap<u64, i32> = forward.into_iter().collect();
        let reverse_map: BTreeMap<u64, i32> = reverse.into_iter().collect();
        let mut positions: BTreeSet<u64> = forward_map.keys().copied().collect();
        positions.extend(reverse_map.keys().copied());
        let mut steps = Vec::with_capacity(positions.len());
        for bp in positions {
            let rel = (bp - seq_lo) as usize;
            let canonical_forward = seq
                .get(rel..rel + k as usize)
                .map(|window| {
                    impg::genome_inference::mem_records::window_is_canonical_forward(window)
                })
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
        Ok(steps)
    };

    // --------------------- bind every census record to its derived key
    // The pattern's position-0 window starts exactly at the occurrence
    // start: at an orientation-0 occurrence the panel step there is the
    // pattern's FIRST anchor node; at an orientation-1 occurrence the
    // path spells the reverse complement, so the step there is the
    // FLIPPED LAST anchor node. Index the derived keys by both endpoint
    // nodes and the anchor count, then verify by sequence.
    let decode_node = |token: u64| -> io::Result<i32> {
        let zigzag = token
            .checked_sub(2)
            .ok_or_else(|| invalid("record node token"))?
            / 2;
        Ok(((zigzag >> 1) as i64 ^ -(zigzag as i64 & 1)) as i32)
    };
    let mut key_by_first: HashMap<(i32, u64), Vec<u32>> = HashMap::new();
    let mut key_by_last: HashMap<(i32, u64), Vec<u32>> = HashMap::new();
    for (key, tokens) in key_tokens.iter().enumerate() {
        let first_node = decode_node(tokens[0])?;
        let last_node = decode_node(*tokens.last().unwrap())?;
        let anchor_count = (tokens.len() / 2 + 1) as u64;
        key_by_first
            .entry((first_node, anchor_count))
            .or_default()
            .push(key as u32);
        key_by_last
            .entry((last_node, anchor_count))
            .or_default()
            .push(key as u32);
    }
    // One key's shape, computed on demand: the canonical walk, span,
    // and one representative read's own record instance.
    struct KeyShape {
        canonical: Vec<(i32, u64)>,
        span: u64,
        read: usize,
        mirrored: bool,
        /// The own anchor positions per canonical anchor index.
        own_positions: Vec<u64>,
    }
    fn key_shape(
        key: u32,
        key_tokens: &[Vec<u64>],
        key_reads: &[Vec<usize>],
        read_records: &[Vec<ReadRecord>],
        k: u64,
    ) -> io::Result<KeyShape> {
        let tokens = &key_tokens[key as usize];
        let canonical = decode_record_tokens(tokens)?;
        ensure(!canonical.is_empty(), "empty record walk")?;
        let span = walk_span(&canonical, k);
        let read = *key_reads[key as usize]
            .first()
            .ok_or_else(|| invalid("record key has no reads"))?;
        let own = read_records[read]
            .iter()
            .find(|record| record.key == key)
            .map(|record| record.walk.clone())
            .ok_or_else(|| invalid("read record mismatch"))?;
        let encoded = impg::sample_mem_bwt::encode_walk(
            &own
                .iter()
                .map(|&(node, pos)| (node, pos as u64))
                .collect::<Vec<_>>(),
        )?;
        let mirrored = encoded != tokens.as_slice();
        let w_lo = own.first().map(|&(_, p)| p as u64).unwrap_or(0);
        let own_positions = own_positions_of(&canonical, mirrored, w_lo);
        // The mirror formula is validated against the representative
        // read's own walk: own anchor t is canonical anchor m-1-t with
        // flipped sign at the computed own position.
        let m = canonical.len();
        ensure(own.len() == m, "own walk length differs from the canonical walk")?;
        for (t, &(own_node, own_pos)) in own.iter().enumerate() {
            let j = if mirrored { m - 1 - t } else { t };
            let expect_node = if mirrored { -canonical[j].0 } else { canonical[j].0 };
            ensure(
                own_node == expect_node && own_pos as u64 == own_positions[j],
                "the own walk disagrees with the mirror formula",
            )?;
        }
        Ok(KeyShape {
            canonical,
            span,
            read,
            mirrored,
            own_positions,
        })
    }
    // Sequence verification of one key at one occurrence: the anchor
    // windows must sit at their projected panel positions and the
    // representative read's own k-mers must equal them (directly when
    // the read is physically forward, as reverse complements when
    // physically reverse).
    /// Verify one key at one occurrence with the pattern origin at
    /// occ.start + shift: every anchor window (the panel-side recipe)
    /// must equal the representative read's own k-mer (directly when
    /// the read is physically forward, as its reverse complement when
    /// physically reverse). The panel segment is fetched from the AGC
    /// — the sequence source of truth — with one base of slack on each
    /// side so the shifted windows stay inside.
    fn verify_key_at(
        shape: &KeyShape,
        occ_path_name: &str,
        occ: &CensusOccurrence,
        reads: &[Vec<u8>],
        fetch: &dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>>,
        k: u64,
        shift: i64,
    ) -> io::Result<bool> {
        let lo = (occ.start as i64 + shift - 1).max(0) as u64;
        let hi = occ.start + shape.span + 1;
        let seq = fetch(occ_path_name, lo, hi)?;
        let origin = occ.start as i64 + shift;
        let read = &reads[shape.read];
        let forward = physical_forward(occ.orientation, shape.mirrored);
        for (j, &(node, canonical_pos)) in shape.canonical.iter().enumerate() {
            let _ = node;
            let window_lo =
                anchor_window_start(canonical_pos, origin as u64, shape.span, k, occ.orientation);
            let rel_lo = window_lo as i64 - lo as i64;
            if rel_lo < 0 || rel_lo + k as i64 > seq.len() as i64 {
                return Ok(false);
            }
            let window = &seq[rel_lo as usize..(rel_lo + k as i64) as usize];
            let own_pos = shape.own_positions[j] as usize;
            let kmer = &read[own_pos..own_pos + k as usize];
            let ok = if forward {
                window == kmer
            } else {
                revcomp(window) == kmer
            };
            if !ok {
                return Ok(false);
            }
        }
        Ok(true)
    }
    let mut census_key: Vec<Option<u32>> = vec![None; census_records.len()];
    // The verified origin shift per occurrence (census start + shift =
    // the sequence-verified pattern origin).
    let mut record_shifts: Vec<Vec<i64>> = Vec::with_capacity(census_records.len());
    for (record, line) in census_records.iter().enumerate() {
        let occ = line
            .occurrences
            .first()
            .ok_or_else(|| invalid("census record has no occurrences"))?;
        // The pattern's first/last anchor node, as spelled by the
        // routing's own canonical-scheme steps AT the census start (the
        // occurrence start is the pattern's position-0 window start in
        // the territory convention; the sequence verification then
        // decides the true origin shift).
        let steps = canonical_steps_near(occ.path, occ.start)?
            .into_iter()
            .filter(|&(bp, _)| bp == occ.start)
            .map(|(_, node)| node)
            .collect::<Vec<i32>>();
        if steps.is_empty() {
            eprintln!(
                "[anchor] NO STEP record {record}: occ path {} start {} orientation {} \
                 anchors {} multiplicity {} raw steps {:?}",
                occ.path,
                occ.start,
                occ.orientation,
                line.anchors,
                line.multiplicity,
                panel
                    .walk_path_range(occ.path, occ.start.saturating_sub(80), occ.start + 80)?
                    .into_iter()
                    .collect::<Vec<(i32, u64)>>(),
            );
        }
        ensure(!steps.is_empty(), "census occurrence start has no panel step")?;
        let mut candidates: BTreeSet<u32> = BTreeSet::new();
        for &step in &steps {
            // Orientation 0: the step is the pattern's FIRST anchor;
            // orientation 1: the step is the flipped LAST anchor. The
            // other endpoint is unconstrained by this lookup — the
            // sequence verification decides.
            let index: &HashMap<(i32, u64), Vec<u32>> = if occ.orientation == 0 {
                &key_by_first
            } else {
                &key_by_last
            };
            let node = if occ.orientation == 0 { step } else { -step };
            if let Some(keys) = index.get(&(node, line.anchors)) {
                for &key in keys {
                    candidates.insert(key);
                }
            }
        }
        if candidates.is_empty() {
            eprintln!(
                "[anchor] BIND FAIL record {record}: occ path {} start {} orientation {} \
                 steps {:?} anchors {} multiplicity {}",
                occ.path,
                occ.start,
                occ.orientation,
                steps,
                line.anchors,
                line.multiplicity,
            );
            // One-off diagnosis: every derived key whose walk touches the
            // occurrence's step nodes, with anchor counts and instances.
            for &step in &steps {
                for (key, tokens) in key_tokens.iter().enumerate() {
                    let walk = decode_record_tokens(tokens)?;
                    if walk.iter().any(|&(node, _)| node == step || node == -step) {
                        eprintln!(
                            "  key {key}: anchors {} instances {} walk {:?}",
                            walk.len(),
                            key_reads[key].len(),
                            walk,
                        );
                    }
                }
            }
        }
        // The verified origin shift per occurrence: the unique shift in
        // {-1, 0, +1} at which every anchor window matches the panel
        // sequence (the AGC fetch is the source of truth).
        let shifts = [0i64, 1, -1];
        let verify_all = |key: u32| -> io::Result<Option<Vec<i64>>> {
            let shape = key_shape(key, &key_tokens, &key_reads, &read_records, k)?;
            let mut per_occurrence = Vec::with_capacity(line.occurrences.len());
            for occ in &line.occurrences {
                let mut found: Vec<i64> = Vec::new();
                for &shift in &shifts {
                    if verify_key_at(
                        &shape,
                        &panel.name_map.path_to_name[occ.path],
                        occ,
                        &reads,
                        &fetch_seq,
                        k,
                        shift,
                    )? {
                        found.push(shift);
                    }
                }
                match found.len() {
                    0 => return Ok(None),
                    1 => per_occurrence.push(found[0]),
                    _ => return Err(invalid("key verifies at two origin shifts")),
                }
            }
            Ok(Some(per_occurrence))
        };
        let n_candidates = candidates.len();
        let mut surviving: Vec<(u32, Vec<i64>)> = Vec::new();
        for key in candidates {
            if let Some(per_occurrence) = verify_all(key)? {
                surviving.push((key, per_occurrence));
            }
        }
        if surviving.len() != 1 {
            eprintln!(
                "[anchor] BIND {} record {record}: {} surviving keys of {} candidates \
                 (occ path {} start {} orientation {} anchors {})",
                if surviving.is_empty() { "FAIL" } else { "AMBIGUOUS" },
                surviving.len(),
                n_candidates,
                occ.path,
                occ.start,
                occ.orientation,
                line.anchors,
            );
        }
        ensure(
            surviving.len() == 1,
            "census record has no unique verifying key",
        )?;
        census_key[record] = Some(surviving[0].0);
        record_shifts.push(surviving.swap_remove(0).1);
    }
    eprintln!(
        "[anchor] census records bound: {} ({:.1}s)",
        census_records.len(),
        started.elapsed().as_secs_f64(),
    );
    rss.probe("binding")?;

    // ----------------------- per-record placements and read variants
    // One read variant: the physical state of one instance (mirror
    // flag, hull start, full read sequence identity).
    struct ReadVariant {
        mirror: bool,
        read: usize,
        w_lo: u64,
        count: u64, // instance count (duplicate reads weighted)
    }
    let mut record_variants: Vec<Vec<ReadVariant>> = Vec::with_capacity(census_records.len());
    for (record, line) in census_records.iter().enumerate() {
        let key = census_key[record].unwrap();
        // Group every instance (key_reads carries one entry per record
        // instance) by (mirror, hull start, full read sequence).
        let mut variants: BTreeMap<(bool, u64, u64), ReadVariant> = BTreeMap::new();
        // Per read: which matching record instance we are consuming
        // (a read can carry one record at several positions).
        let mut instance_cursor: HashMap<usize, usize> = HashMap::new();
        for &read in &key_reads[key as usize] {
            let cursor = instance_cursor.entry(read).or_insert(0);
            let mut matching = read_records[read]
                .iter()
                .filter(|record| record.key == key)
                .skip(*cursor);
            let entry = matching
                .next()
                .ok_or_else(|| invalid("key_reads entry without a matching record"))?;
            *cursor += 1;
            let own = &entry.walk;
            let encoded = impg::sample_mem_bwt::encode_walk(
                &own
                    .iter()
                    .map(|&(node, pos)| (node, pos as u64))
                    .collect::<Vec<_>>(),
            )?;
            let mirrored = encoded != key_tokens[key as usize].as_slice();
            let w_lo = own.first().map(|&(_, p)| p as u64).unwrap_or(0);
            let vkey = (mirrored, w_lo, fnv1a64(&reads[read]));
            let weight = read_multiplicity[read] as u64;
            match variants.get_mut(&vkey) {
                Some(variant) => variant.count += weight,
                None => {
                    variants.insert(
                        vkey,
                        ReadVariant {
                            mirror: mirrored,
                            read,
                            w_lo,
                            count: weight,
                        },
                    );
                }
            }
        }
        let variants: Vec<ReadVariant> = variants.into_values().collect();
        let total: u64 = variants.iter().map(|v| v.count).sum();
        ensure(
            total == line.multiplicity,
            "variant instance counts do not sum to the census multiplicity",
        )?;
        record_variants.push(variants);
    }
    rss.probe("variants")?;

    // ------------------ the full-domain placement cross-validation
    let mut store = RowStore::new(&panel, &path_of_name, &fetch_seq);
    let mut occurrence_total = 0u64;
    let mut interval_match = 0u64;
    let mut node_list_match = 0u64;
    let mut own_row_checked = 0u64;
    let mut own_row_missing = 0u64;
    let mut hull_checked = 0u64;
    let mut hull_fetches = 0u64;
    {
        let mut report = BufWriter::new(File::create(&options.out)?);
        for (record, line) in census_records.iter().enumerate() {
            let key = census_key[record].unwrap();
            let canonical = decode_record_tokens(&key_tokens[key as usize])?;
            let span = walk_span(&canonical, k);
            let variants = &record_variants[record];
            let shifts = &record_shifts[record];
            ensure(span <= READ_LENGTH as u64, "record span exceeds the read length")?;
            // Per occurrence: the re-derived intervals + the census
            // cross-validation (one path walk per path over the
            // record's occurrence hulls on that path — the census's
            // own recipe).
            let mut by_path: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for (position, occ) in line.occurrences.iter().enumerate() {
                by_path.entry(occ.path).or_default().push(position);
            }
            for (path, positions) in &by_path {
                let mut hulls: Vec<(u64, u64)> = Vec::new();
                for &position in positions {
                    let occ = &line.occurrences[position];
                    for interval in
                        occurrence_intervals(&canonical, k, occ.start, occ.orientation)
                    {
                        hulls.push(interval);
                    }
                }
                let lo_min = hulls.iter().map(|&(lo, _)| lo).min().unwrap_or(0);
                let hi_max = hulls.iter().map(|&(_, hi)| hi).max().unwrap_or(0);
                if hi_max <= lo_min {
                    continue;
                }
                let mut steps: Vec<(u64, i32)> = panel
                    .walk_path_range(*path, lo_min, hi_max)?
                    .into_iter()
                    .map(|(node, bp)| (bp, node))
                    .collect();
                steps.sort_unstable_by_key(|&(bp, _)| bp);
                for &position in positions {
                    let occ = &line.occurrences[position];
                    let intervals =
                        occurrence_intervals(&canonical, k, occ.start, occ.orientation);
                    ensure(
                        intervals.len() == occ.intervals.len(),
                        "re-derived interval count differs from the census",
                    )?;
                    interval_match += 1;
                    let mut per_interval: Vec<Vec<u32>> = Vec::with_capacity(intervals.len());
                    for &(lo, hi) in &intervals {
                        let start = steps.partition_point(|&(bp, _)| bp < lo);
                        let end = steps.partition_point(|&(bp, _)| bp.saturating_add(k) <= hi);
                        per_interval.push(
                            steps[start..end.max(start)]
                                .iter()
                                .map(|&(_, n)| n.unsigned_abs())
                                .collect(),
                        );
                    }
                    ensure(
                        per_interval == occ.intervals,
                        "re-derived covered-node lists differ from the census",
                    )?;
                    node_list_match += 1;
                }
            }
            // The own-row, anchor-window, and hull assertions per
            // occurrence, at the SEQUENCE-VERIFIED origin (the census
            // start plus the verified shift).
            for (position, occ) in line.occurrences.iter().enumerate() {
                occurrence_total += 1;
                let origin = (occ.start as i64 + shifts[position]).max(0) as u64;
                let mut seg: Option<Vec<u8>> = None;
                for &window in &occ.partitions {
                    let partition = window_partition[window as usize];
                    let map = maps
                        .get(&partition)
                        .ok_or_else(|| invalid("axis partition map missing"))?;
                    store.ensure_partition(&map.members, partition)?;
                    for fold_index in store.rows_of_partition_path(partition, occ.path).to_vec()
                    {
                        let row = &store.folds[fold_index];
                        if row.start <= origin && origin + span <= row.end {
                            // The containing row's SEQUENCE at the
                            // projected span is the hull check below
                            // (the row walks carry only one frame's
                            // step selection, so a walk-position
                            // assertion does not hold in general; the
                            // sequence is the truth).
                            seg = Some(
                                row.seq[(origin - row.start) as usize
                                    ..(origin - row.start + span) as usize]
                                    .to_vec(),
                            );
                            own_row_checked += 1;
                            break;
                        }
                    }
                    if seg.is_some() {
                        break;
                    }
                }
                let seg = match seg {
                    Some(seg) => seg,
                    None => {
                        own_row_missing += 1;
                        hull_fetches += 1;
                        fetch_seq(
                            &panel.name_map.path_to_name[occ.path],
                            origin,
                            origin + span,
                        )?
                    }
                };
                // The hull equality per variant (flanks are NOT
                // asserted — they are the skipped evidence).
                for variant in variants {
                    let read = &reads[variant.read];
                    let hull = &read[variant.w_lo as usize..(variant.w_lo + span) as usize];
                    let ok = if physical_forward(occ.orientation, variant.mirror) {
                        hull == seg.as_slice()
                    } else {
                        hull == revcomp(&seg).as_slice()
                    };
                    ensure(ok, "the read hull differs from the placed segment")?;
                    hull_checked += 1;
                }
            }
            // The receipt line.
            let occurrences: Vec<serde_json::Value> = line
                .occurrences
                .iter()
                .enumerate()
                .map(|(position, occ)| {
                    json!({
                        "path": occ.path,
                        "start": occ.start,
                        "shift": shifts[position],
                        "orientation": occ.orientation,
                        "partitions": occ.partitions,
                        "intervals": occ.intervals.len(),
                        "pf": variants.iter()
                            .map(|v| physical_forward(occ.orientation, v.mirror))
                            .collect::<Vec<bool>>(),
                    })
                })
                .collect();
            let variant_values: Vec<serde_json::Value> = variants
                .iter()
                .map(|v| {
                    let read = &reads[v.read];
                    json!([
                        if v.mirror { 1 } else { 0 },
                        String::from_utf8_lossy(&read[..v.w_lo as usize]),
                        String::from_utf8_lossy(&read[(v.w_lo + span) as usize..]),
                        v.count,
                    ])
                })
                .collect();
            // The skipped positions are the read flanks plus any
            // interior hull gaps (measured zero on this sample; emitted
            // honestly per mirror state when present).
            let mut hull_gaps: Vec<serde_json::Value> = Vec::new();
            for mirror in [false, true] {
                let own_positions = own_positions_of(&canonical, mirror, 0);
                let skipped = skipped_positions(&own_positions, k, span);
                let interior: Vec<(u64, u64)> = skipped
                    .into_iter()
                    .filter(|&(lo, hi)| lo > 0 && hi < span)
                    .collect();
                if !interior.is_empty() {
                    hull_gaps.push(json!([
                        if mirror { 1 } else { 0 },
                        interior.iter().map(|&(lo, hi)| json!([lo, hi])).collect::<Vec<_>>(),
                    ]));
                }
            }
            serde_json::to_writer(
                &mut report,
                &json!({
                    "record": record as u64,
                    "multiplicity": line.multiplicity,
                    "k": k,
                    "span": span,
                    "anchors": canonical.iter()
                        .map(|&(node, pos)| json!([node, pos]))
                        .collect::<Vec<_>>(),
                    "reads": variant_values,
                    "hull_gaps": hull_gaps,
                    "occurrences": occurrences,
                }),
            )?;
            writeln!(report)?;
        }
        report.flush()?;
    }
    eprintln!(
        "[anchor] placement layer: {occurrence_total} occurrences, \
         interval-count matches {interval_match}, node-list matches {node_list_match}, \
         own-row checked {own_row_checked} (missing {own_row_missing}, fetches {hull_fetches}), \
         hull checks {hull_checked}",
    );
    rss.probe("placement_layer")?;

    // ------------------------ the per-path graph-context sample
    // The evenly-spaced stated sample of occurrences (census order,
    // stride = ceil(total / sample_count)).
    let total_occ: u64 = census_records
        .iter()
        .map(|l| l.occurrences.len() as u64)
        .sum();
    let stride = ((total_occ + options.sample_count - 1) / options.sample_count).max(1);
    let mut context = BufWriter::new(File::create(&options.context_out)?);
    let mut context_lines = 0u64;
    let mut context_rows = 0u64;
    let mut context_lookups = 0u64;
    let mut flat = 0u64;
    for (record, line) in census_records.iter().enumerate() {
        let key = census_key[record].unwrap();
        let canonical = decode_record_tokens(&key_tokens[key as usize])?;
        let span = walk_span(&canonical, k);
        let variants = &record_variants[record];
        for (position, occ) in line.occurrences.iter().enumerate() {
            if flat % stride != 0 {
                flat += 1;
                continue;
            }
            flat += 1;
            // The candidate rows: the member rows of the axis
            // partitions of every touched window, deduplicated.
            let mut candidate_folds: BTreeSet<usize> = BTreeSet::new();
            for &window in &occ.partitions {
                let partition = window_partition[window as usize];
                let map = maps
                    .get(&partition)
                    .ok_or_else(|| invalid("axis partition map missing"))?;
                store.ensure_partition(&map.members, partition)?;
                for &fold_index in &store.by_partition[&partition] {
                    candidate_folds.insert(fold_index);
                }
            }
            let mut row_values: Vec<serde_json::Value> = Vec::new();
            for &fold_index in &candidate_folds {
                let row = &store.folds[fold_index];
                // Traversal: the row's walk shares at least one
                // anchor node with the required sign.
                let traverses = canonical.iter().any(|&(node, _)| {
                    let required = required_sign(node, occ.orientation);
                    row.positions_of(required.unsigned_abs())
                        .iter()
                        .any(|&(_, sign)| sign == required.signum())
                });
                if !traverses {
                    continue;
                }
                // The correspondence is orientation-based (the
                // pattern's required signs) — identical for every
                // variant; the skipped lookups are mirror-based and
                // hull-relative, so all variants of one mirror state
                // share them.
                let positions_of = |node: u32| row.positions_of(node);
                let (pinned, ambiguous) =
                    anchor_correspondence(&canonical, occ.orientation, &positions_of);
                // Sequence-verify every pinned anchor against the row's
                // own sequence: the pinned position is the true window
                // start (the row walks are one frame's step selection,
                // so a pinned anchor not visible in the walk is simply
                // unshared; one that IS visible must agree with the
                // sequence).
                let mut per_mirror: Vec<serde_json::Value> = Vec::new();
                for mirror in [false, true] {
                    let forward = physical_forward(occ.orientation, mirror);
                    let own_positions = own_positions_of(&canonical, mirror, 0);
                    // The hull-relative skipped d-ranges shared by every
                    // variant of this mirror state: the flank unions over
                    // the variants' hull starts (d = own position minus
                    // w_lo, so a variant's left flank is d in
                    // [-w_lo, 0) and its right flank d in
                    // [span, 150 - w_lo)) plus the interior hull gaps.
                    let wlos: Vec<u64> = variants
                        .iter()
                        .filter(|v| v.mirror == mirror)
                        .map(|v| v.w_lo)
                        .collect();
                    if wlos.is_empty() {
                        continue;
                    }
                    let (min_w, max_w) = (
                        wlos.iter().copied().min().unwrap(),
                        wlos.iter().copied().max().unwrap(),
                    );
                    let mut skipped_hull: Vec<(i64, i64)> = Vec::new();
                    skipped_hull.push((-(max_w as i64), 0));
                    skipped_hull.extend(
                        skipped_positions(&own_positions, k, span)
                            .into_iter()
                            .filter(|&(lo, hi)| lo > 0 && hi < span)
                            .map(|(lo, hi)| (lo as i64, hi as i64)),
                    );
                    let right_lo = span as i64;
                    let right_hi = READ_LENGTH as i64 - min_w as i64;
                    if right_hi > right_lo {
                        skipped_hull.push((right_lo, right_hi));
                    }
                    skipped_hull.sort_unstable();
                    if skipped_hull.is_empty() {
                        continue;
                    }
                    let mut lookups: Vec<serde_json::Value> = Vec::new();
                    for (lo, hi) in &skipped_hull {
                        for d in *lo..*hi {
                            // The serving anchor: the pinned anchor with
                            // the greatest own position <= the base (the
                            // left bracket; the outermost pinned anchors
                            // serve the flanks), else the least. d is
                            // hull-relative (own position minus w_lo), so
                            // the left flank carries negative d.
                            let serving = {
                                let mut best: Option<(usize, i64)> = None;
                                for (index, pinned_bp) in pinned.iter().enumerate() {
                                    if pinned_bp.is_none() {
                                        continue;
                                    }
                                    let own_pos = own_positions[index] as i64;
                                    let better = match best {
                                        None => true,
                                        Some((_, best_pos)) => {
                                            let (own_left, best_left) =
                                                (own_pos <= d, best_pos <= d);
                                            if own_left && !best_left {
                                                true
                                            } else if own_left == best_left {
                                                if own_left {
                                                    own_pos > best_pos
                                                } else {
                                                    own_pos < best_pos
                                                }
                                            } else {
                                                false
                                            }
                                        }
                                    };
                                    if better {
                                        best = Some((index, own_pos));
                                    }
                                }
                                best
                            };
                            let (coord, range, base) = match serving {
                                Some((index, own_pos)) => {
                                    let r = pinned[index].unwrap();
                                    let coord = if forward {
                                        r as i64 + (d - own_pos)
                                    } else {
                                        r as i64 + (own_pos + k as i64 - 1 - d)
                                    };
                                    let (range, base) = if coord < 0 {
                                        ((0usize, 0usize), None)
                                    } else {
                                        (
                                            row.covering_index_range(coord as u64, k),
                                            row.base_at(coord),
                                        )
                                    };
                                    (coord, range, base)
                                }
                                None => (-1i64, (0usize, 0usize), None),
                            };
                            context_lookups += 1;
                            lookups.push(json!([
                                d,
                                coord,
                                range.0,
                                range.1,
                                base.map(|b| b as char),
                            ]));
                        }
                    }
                    per_mirror.push(json!({
                        "mirror": if mirror { 1 } else { 0 },
                        "anchors": pinned,
                        "ambiguous": ambiguous,
                        "lookups": lookups,
                    }));
                }
                context_rows += 1;
                row_values.push(json!({
                    "path": panel.name_map.path_to_name[row.path],
                    "start": row.start,
                    "end": row.end,
                    "members": row.members,
                    "per_mirror": per_mirror,
                }));
            }
            serde_json::to_writer(
                &mut context,
                &json!({
                    "record": record as u64,
                    "occ": position,
                    "path": occ.path,
                    "start": occ.start,
                    "shift": record_shifts[record][position],
                    "orientation": occ.orientation,
                    "span": span,
                    "rows": row_values,
                }),
            )?;
            writeln!(context)?;
            context_lines += 1;
        }
    }
    context.flush()?;
    eprintln!(
        "[anchor] context sample: stride {stride}, {context_lines} sampled occurrences, \
         {context_rows} traversing rows, {context_lookups} per-base lookups",
    );
    rss.probe("context_sample")?;
    eprintln!("[anchor] done in {:.1}s", started.elapsed().as_secs_f64(),);
    Ok(())
}

// ---------------------------------------------------------------------------
// Unit tests: the strand logic, the projection equations, the skipped
// extraction, and the correspondence rule (the previous attempt died
// fixing exactly this — it is proven first, on synthetic walks).
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const K: u64 = 8;

    fn canonical() -> Vec<(i32, u64)> {
        // Three anchors at canonical positions 0, 10, 30; span = 38.
        vec![(5, 0), (7, 10), (-9, 30)]
    }

    #[test]
    fn own_positions_mirror_the_canonical() {
        let walk = canonical();
        let span = walk_span(&walk, K);
        assert_eq!(span, 38);
        // Unmirrored: own position = w_lo + c_j.
        let own = own_positions_of(&walk, false, 4);
        assert_eq!(own, vec![4, 14, 34]);
        // Mirrored: canonical anchor j sits at w_lo + (c_last - c_j),
        // and the own ORDER reverses (the mirrored read spells the
        // pattern's reverse complement).
        let own = own_positions_of(&walk, true, 4);
        assert_eq!(own, vec![4 + 30, 4 + 20, 4 + 0]);
    }

    #[test]
    fn physical_strand_is_orientation_xor_mirror() {
        // The census orientation refers to the canonical token walk.
        assert!(physical_forward(0, false));
        assert!(!physical_forward(1, false));
        assert!(physical_forward(1, true));
        assert!(!physical_forward(0, true));
    }

    #[test]
    fn projection_equations_hold_at_anchor_windows() {
        // For every (orientation, mirror) state: the own window of
        // canonical anchor j maps to exactly the panel window whose
        // start the anchor equation gives — the projection and the
        // anchor equations are the same map.
        let walk = canonical();
        let span = walk_span(&walk, K);
        let (start, w_lo) = (1000u64, 6u64);
        for orientation in [0u8, 1] {
            for mirror in [false, true] {
                let forward = physical_forward(orientation, mirror);
                let own = own_positions_of(&walk, mirror, w_lo);
                for (j, &p) in own.iter().enumerate() {
                    let c_j = walk[j].1;
                    let window_start = if forward {
                        project_own_position(p, start, span, w_lo, true)
                    } else {
                        project_own_position(p + K - 1, start, span, w_lo, false)
                    };
                    let expect =
                        anchor_window_start(c_j, start, span, K, orientation) as i64;
                    assert_eq!(window_start, expect, "orientation {orientation} mirror {mirror} j {j}");
                    // The other window boundary projects K-1 away on
                    // the same side.
                    let window_end = if forward {
                        project_own_position(p + K - 1, start, span, w_lo, true)
                    } else {
                        project_own_position(p, start, span, w_lo, false)
                    };
                    assert_eq!((window_end - window_start).abs(), (K - 1) as i64);
                }
            }
        }
    }

    #[test]
    fn reverse_projection_is_the_mirror_equation() {
        // start + span - 1 - (i - w_lo): hand-checked cases.
        assert_eq!(project_own_position(0, 100, 38, 6, false), 143);
        assert_eq!(project_own_position(6, 100, 38, 6, false), 137);
        assert_eq!(project_own_position(37, 100, 38, 6, false), 106);
        // The hull end projects to the pattern's position-0 base.
        assert_eq!(project_own_position(43, 100, 38, 6, false), 100);
        assert_eq!(project_own_position(5, 100, 38, 6, true), 99);
    }

    #[test]
    fn skipped_positions_are_the_window_complement() {
        // Windows [4,12), [14,22) with a gap: the complement inside
        // the read carries the gap and both flanks.
        let skipped = skipped_positions(&[4, 14], K, 30);
        assert_eq!(skipped, vec![(0, 4), (12, 14), (22, 30)]);
        // Contiguous windows leave only the flanks.
        let skipped = skipped_positions(&[0, K], K, 30);
        assert_eq!(skipped, vec![(2 * K, 30)]);
        // A full-coverage walk skips nothing.
        let skipped = skipped_positions(&[0, 5, 10, 15], K, 23);
        assert_eq!(skipped, vec![]);
    }

    #[test]
    fn occurrence_intervals_merge_panel_windows() {
        let walk = canonical();
        let span = walk_span(&walk, K);
        // Orientation 0: windows at start + {0,10,30}: the first two
        // abut (10 == 0+8+2? no: gap 2) — actually [1000,1008),
        // [1010,1018), [1030,1038): the first two are 2 apart, so
        // three separate intervals.
        let intervals = occurrence_intervals(&walk, K, 1000, 0);
        assert_eq!(
            intervals,
            vec![(1000, 1008), (1010, 1018), (1030, 1038)]
        );
        // Orientation 1: windows at start + span - K - c_j
        // = {1030, 1020, 1000} — the mirrored window set.
        let intervals = occurrence_intervals(&walk, K, 1000, 1);
        assert_eq!(
            intervals,
            vec![(1000, 1008), (1020, 1028), (1030, 1038)]
        );
        let _ = span;
    }

    #[test]
    fn correspondence_pins_unique_monotone_assignments() {
        let walk = canonical();
        // A row walk containing each anchor node exactly once, in
        // order, at shifted coordinates (an indel between the second
        // and third anchor).
        let row_positions = |node: u32| -> Vec<(u64, i32)> {
            match node {
                5 => vec![(500, 1)],
                7 => vec![(510, 1)],
                9 => vec![(545, -1)],
                _ => vec![],
            }
        };
        let (pinned, ambiguous) = anchor_correspondence(&walk, 0, &row_positions);
        assert_eq!(pinned, vec![Some(500), Some(510), Some(545)]);
        assert!(ambiguous.is_empty());
        // Orientation 1 requires the flipped signs and holds the
        // anchors in REVERSED row order (the rc pattern's first anchor
        // sits at the highest row coordinate).
        let row_positions = |node: u32| -> Vec<(u64, i32)> {
            match node {
                5 => vec![(530, -1)],
                7 => vec![(520, -1)],
                9 => vec![(500, 1)],
                _ => vec![],
            }
        };
        let (pinned, ambiguous) = anchor_correspondence(&walk, 1, &row_positions);
        assert_eq!(pinned, vec![Some(530), Some(520), Some(500)]);
        assert!(ambiguous.is_empty());
        // A repeated middle node with TWO monotone options: both
        // options keep the assignment monotone, so the middle anchor
        // is AMBIGUOUS while the outer anchors stay pinned.
        let row_positions = |node: u32| -> Vec<(u64, i32)> {
            match node {
                5 => vec![(500, 1)],
                7 => vec![(510, 1), (520, 1)],
                9 => vec![(545, -1)],
                _ => vec![],
            }
        };
        let (pinned, ambiguous) = anchor_correspondence(&walk, 0, &row_positions);
        assert_eq!(pinned, vec![Some(500), None, Some(545)]);
        assert_eq!(ambiguous, vec![1]);
        // An anchor absent from the row is None but not ambiguous.
        let row_positions = |node: u32| -> Vec<(u64, i32)> {
            match node {
                5 => vec![(500, 1)],
                9 => vec![(545, -1)],
                _ => vec![],
            }
        };
        let (pinned, ambiguous) = anchor_correspondence(&walk, 0, &row_positions);
        assert_eq!(pinned, vec![Some(500), None, Some(545)]);
        assert!(ambiguous.is_empty());
    }

    #[test]
    fn token_roundtrip_and_span() {
        // encode -> canonical -> decode reproduces the walk positions
        // (rebased to 0) for both mirror states.
        let walk = canonical();
        let encoded = impg::sample_mem_bwt::encode_walk(&walk).unwrap();
        let decoded = decode_record_tokens(&encoded).unwrap();
        assert_eq!(decoded, walk);
        let rc = impg::sample_mem_bwt::reverse_complement(&encoded);
        assert_eq!(decode_record_tokens(&rc).unwrap().len(), walk.len());
    }
}
