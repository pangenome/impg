//! THE SCORING LAYER (slice B of the local-realignment evidence model;
//! owner direction 2026-11-06, slice A's anchor-projection layer stands
//! as the placement machinery of record). Assessment-side only.
//!
//! THE MODEL, derived and stated (no thresholds, no tuning constants):
//!
//! 1. THE SCORING RATE comes from the reads' own base qualities: the
//!    sample's FASTQ carries one uniform quality byte ('I' = Phred 40,
//!    asserted by streaming the whole file), so the per-base error
//!    probability is eps = 10^(-40/10) = 1e-4 and the substitution
//!    model gives P(match) = 1 - eps, P(mismatch) = eps/3 (uniform over
//!    the three non-read bases). The per-base log weights are
//!    A = ln(1 - eps) and B = ln(eps/3). The measured error profile is
//!    stated beside them: slice A measured zero unmatched placement bp
//!    over 1,063,522 verified occurrences (every base of every chained
//!    placement span lies inside an exactly-matched anchor window), so
//!    the quality byte is the sample's own statement of its error rate.
//!
//! 2. PER (RECORD, CANDIDATE ROW) SKIPPED-BASE VOTES, computed
//!    in-process from slice A's committed machinery (the anchor
//!    correspondence, the serving-anchor projection, the mirror
//!    equations — the same functions the slice A receipts and checker
//!    validated; the tens-of-GB full-fidelity per-path context is
//!    computed here, never materialized). A candidate row's
//!    log-likelihood for one read instance is
//!
//!        LL(read | row) = -ln(2 * (len_row - 149))
//!                         + logsumexp over anchored placements of
//!                           [ backbone term + skipped-base evidence ]
//!
//!    * THE BACKBONE TERM is the MEM-anchored shared material: the
//!      merged bp of the PINNED anchor windows (the anchors the row's
//!      own walk shares with the required sign, placed by the
//!      unique-monotone correspondence), each a match by node identity
//!      and orientation sign. It is identical across near-identical
//!      paths (rows sharing the anchor chain share the backbone) and
//!      cancels in their comparisons; it is computed once per
//!      (record, placement), never per base.
//!    * THE SKIPPED-BASE EVIDENCE is the read's bases outside its
//!      merged anchor windows (the flanks; interior hull gaps measured
//!      zero on this sample, emitted honestly per mirror state when
//!      present): each skipped base is projected onto the row through
//!      its serving anchor (the bracketing pinned anchor, slice A's
//!      committed rule) and votes A on a match, B on a mismatch
//!      against the row's spelled sequence (the read base compared as
//!      its complement at physically-reverse placements). Bases
//!      projecting outside the row's extent carry no vote (the row
//!      spells nothing there) — and any placement whose voted bases do
//!      not all project inside the extent is not a placement of the
//!      model (the read must lie fully inside the row, the uniform
//!      placement prior's own support).
//!    * THE PLACEMENT PRIOR is uniform over the row's full-read
//!      placements, both strands: 2 * (len_row - 149) of them — the
//!      derived count, stated; the anchored placements are the
//!      dominant terms of the sum (the exactness receipt measures
//!      the dominance directly).
//!    * Bases inside anchor windows the row does NOT pin (absent or
//!      ambiguous anchors) abstain — the named remainder from slice A
//!      ("the within-anchor variant positions"), counted and reported
//!      per locus, not silently scored.
//!
//! 3. NO SHARE-SPLITTING: each candidate row is an independent
//!    hypothesis; a record votes ONCE per row (its instances'
//!    likelihoods under that row are computed independently of every
//!    other row and of how many occurrences discovered the placement —
//!    multiple occurrences of one record pin the SAME placement on a
//!    row and count once). The record's instances are its evidence:
//!    every distinct read variant (mirror state, hull start, sequence)
//!    votes with its own instance count.
//!
//! 4. THE LOCALITY RULE (the ruling for occurrences outside built
//!    axis-partition rows, measured before implementation): a locus's
//!    evidence set is the records with at least one occurrence
//!    TOUCHING the locus's window (the routing's own territory-touch
//!    rule — the same convention that defines the current instrument's
//!    per-locus record sets), and the occurrences that vote at the
//!    locus are exactly those touching occurrences, whether or not
//!    the occurrence's own placement is contained in a built
//!    axis-partition row. The measured outside class is row-edge
//!    overhangs (the read's anchor hull crosses its own path's
//!    partition-membership row boundary) plus territory-extension
//!    material — both window-seam artifacts of the alignment-induced
//!    membership, and the likelihood a read supports a row is a
//!    property of (read, row sequence), not of a row boundary; the
//!    placement-validity rule already guarantees such a read testifies
//!    only where the row actually spells its bases.
//!
//! 5. CANDIDATES fold by identical sequence AND walk (coalesced copies
//!    are ONE hypothesis — the identical-through-graph pair
//!    BTE#3/#4 block28_contig1 folds by construction; the fold is
//!    verified against the committed GFAs by the checker). Classes are
//!    all unordered fold pairs i <= j; class LL = sum over the locus's
//!    evidence units (record, variant) of count * ln(0.5 e^{LL_i} +
//!    0.5 e^{LL_j}) — the equal 1/2 mixture is the balanced 15x/15x
//!    diploid prior. QUAL is the existing cluster machinery verbatim
//!    (spectrum_knee / cluster_form_qual / qual_p / qual_from_p and
//!    the material-distance spectrum, mirrored pure functions from
//!    panel_route_spine/cosine_probe.rs) over the likelihood ratios.
//!
//! 6. THE EXACTNESS PROOF: for EVERY scored (unit, fold, placement) in
//!    the pilot domain the factorized score is asserted equal to the
//!    direct per-base comparison (the backbone verified base-by-base
//!    against the row's sequence, not trusted from node identity) —
//!    integer (match, mismatch) equality, loudly failing. A measured
//!    evenly-spaced sample of (read, candidate row) pairs — plus the
//!    named strata: the truth folds, the identical-through-graph fold,
//!    and pairs with mismatch votes (the variant-pocket class) — is
//!    emitted with vote-level detail for the independent checker,
//!    which re-derives the scores from the committed GFAs, runs the
//!    bounded-offset dominance scan and the bounded edit-distance
//!    alignment on the sample.

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
// Options and small utilities (the slice A runner conventions).
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
    /// The component lane (e.g. S288C#0#chrI) — the axis path whose
    /// member rows define the window-to-partition mapping.
    #[arg(long)]
    component: String,
    /// The committed multi-matching census receipt (records, placements).
    #[arg(long)]
    census: PathBuf,
    /// The sample FASTQ (read identity + the quality-byte assertion).
    #[arg(long)]
    reads: PathBuf,
    /// The record-derivation cache (the slice A binary; a pure
    /// function of the FASTQ and the panel).
    #[arg(long)]
    derive_cache: PathBuf,
    /// The scored loci (window ids of the component, comma-separated,
    /// e.g. 2,4,7) — the pilot set.
    #[arg(long)]
    loci: String,
    /// The likelihood receipt path (sidecars are derived from it).
    #[arg(long)]
    out: PathBuf,
    /// The exactness-sample receipt path (vote-level detail for the
    /// independent checker).
    #[arg(long)]
    exactness_out: PathBuf,
    /// The evenly-spaced exactness sample count per locus (a stated
    /// measurement sample; the in-process assertion covers the whole
    /// pilot domain regardless).
    #[arg(long, default_value_t = 100)]
    exactness_sample: u64,
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

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
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
        eprintln!("[score] rss {stage}: {rss_kb} kB");
        Ok(rss_kb)
    }
}

// ---------------------------------------------------------------------------
// Receipt inputs (slice A structures, verbatim).
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
// The projection core (slice A's committed machinery, verbatim).
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
// FASTQ streaming and the record-derivation cache (slice A, verbatim;
// the cache is a pure function of the FASTQ and the panel).
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
// Candidate rows (slice A's RowStore, verbatim semantics).
// ---------------------------------------------------------------------------

struct RowFold {
    #[allow(dead_code)]
    path: usize,
    start: u64,
    end: u64,
    seq: Vec<u8>,
    /// All walk steps overlapping the row range (absolute bp, signed
    /// node), sorted by bp — the committed partition GFAs' own P-line
    /// walks.
    walk: Vec<(u64, i32)>,
    /// (partition, member index) pairs of every partition holding this
    /// row.
    members: Vec<(u32, usize)>,
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
                    let index = self.folds.len();
                    self.folds.push(RowFold {
                        path: path_idx,
                        start: member.start,
                        end,
                        seq,
                        walk,
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
// The cluster machinery (mirrored pure functions from
// panel_route_spine/cosine_probe.rs — the existing instrument of
// record's own functions; verbatim semantics, no threshold altered).
// ---------------------------------------------------------------------------

type GraphNode = u32;
type PackedEdge = u64;

fn pack_edge(left: GraphNode, right: GraphNode) -> PackedEdge {
    (left as u64) << 32 | right as u64
}

fn multiset_overlap(left: &[(u64, u32)], right: &[(u64, u32)]) -> f64 {
    let (mut i, mut j, mut total) = (0usize, 0usize, 0u64);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                total += left[i].1 as u64 * right[j].1 as u64;
                i += 1;
                j += 1;
            }
        }
    }
    total as f64
}

fn merged_multiset(left: &[(u64, u32)], right: &[(u64, u32)]) -> Vec<(u64, u32)> {
    let (mut i, mut j) = (0usize, 0usize);
    let mut merged: Vec<(u64, u32)> = Vec::with_capacity(left.len() + right.len());
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => {
                merged.push(left[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                merged.push(right[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                merged.push((left[i].0, left[i].1 + right[j].1));
                i += 1;
                j += 1;
            }
        }
    }
    merged.extend_from_slice(&left[i..]);
    merged.extend_from_slice(&right[j..]);
    merged
}

fn differing_observed_mass(
    left: &[(u64, u32)],
    right: &[(u64, u32)],
    observed: &dyn Fn(u64) -> f64,
) -> f64 {
    let (mut i, mut j, mut total) = (0usize, 0usize, 0.0f64);
    while i < left.len() && j < right.len() {
        match left[i].0.cmp(&right[j].0) {
            std::cmp::Ordering::Less => {
                total += observed(left[i].0) * left[i].1 as f64;
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                total += observed(right[j].0) * right[j].1 as f64;
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                let delta = (left[i].1 as i64 - right[j].1 as i64).unsigned_abs();
                if delta > 0 {
                    total += observed(left[i].0) * delta as f64;
                }
                i += 1;
                j += 1;
            }
        }
    }
    while i < left.len() {
        total += observed(left[i].0) * left[i].1 as f64;
        i += 1;
    }
    while j < right.len() {
        total += observed(right[j].0) * right[j].1 as f64;
        j += 1;
    }
    total
}

fn signature_cosine_distance(
    nodes_a: &[(u64, u32)],
    edges_a: &[(u64, u32)],
    nodes_b: &[(u64, u32)],
    edges_b: &[(u64, u32)],
) -> f64 {
    let norm = |nodes: &[(u64, u32)], edges: &[(u64, u32)]| {
        nodes.iter().map(|&(_, mult)| mult as f64 * mult as f64).sum::<f64>()
            + edges.iter().map(|&(_, mult)| mult as f64 * mult as f64).sum::<f64>()
    };
    let (norm_a, norm_b) = (norm(nodes_a, edges_a), norm(nodes_b, edges_b));
    if norm_a == 0.0 && norm_b == 0.0 {
        return 0.0;
    }
    if norm_a == 0.0 || norm_b == 0.0 {
        return 1.0;
    }
    let overlap = multiset_overlap(nodes_a, nodes_b) + multiset_overlap(edges_a, edges_b);
    1.0 - overlap / (norm_a * norm_b).sqrt()
}

fn class_physical_pairs(same_row: bool, members_a: usize, members_b: usize) -> u64 {
    if same_row {
        members_a as u64 * (members_a as u64 + 1) / 2
    } else {
        members_a as u64 * members_b as u64
    }
}

fn median_of(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).expect("finite relative jumps"));
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        values[middle]
    } else {
        (values[middle - 1] + values[middle]) / 2.0
    }
}

struct KneeCut {
    has_knee: bool,
    cut: f64,
}

fn spectrum_knee(sorted: &[f64]) -> KneeCut {
    let positives: Vec<f64> = sorted.iter().copied().filter(|&d| d > 0.0).collect();
    if positives.len() < 2 {
        return KneeCut { has_knee: false, cut: 0.0 };
    }
    let jumps: Vec<f64> = positives
        .windows(2)
        .map(|pair| (pair[1] - pair[0]) / pair[1])
        .collect();
    let mut center = jumps.clone();
    let typical = median_of(&mut center);
    let (mut best, mut best_index) = (f64::NEG_INFINITY, 0usize);
    for (index, &jump) in jumps.iter().enumerate() {
        if jump >= best {
            best = jump;
            best_index = index;
        }
    }
    if best > typical {
        KneeCut { has_knee: true, cut: positives[best_index] }
    } else {
        KneeCut { has_knee: false, cut: 0.0 }
    }
}

fn qual_p(s_win: f64, called_classes: usize, alternative: Option<f64>) -> Option<f64> {
    if called_classes == 0 {
        return None;
    }
    let alternative = alternative.unwrap_or(0.0);
    let denominator = called_classes as f64 * s_win + alternative;
    (denominator > 0.0).then(|| s_win / denominator)
}

fn qual_from_p(p: f64) -> Option<f64> {
    (p < 1.0).then(|| -10.0 * (1.0 - p).log10())
}

struct ClusterQual {
    knee: Option<f64>,
    cut: f64,
    shape: &'static str,
    cluster_size: usize,
    k: usize,
    excluded: Vec<bool>,
    alternative: Option<f64>,
    p: Option<f64>,
}

fn cluster_form_qual(s_win: f64, spectrum: &[(f64, f64)], tied: &[(usize, &[f64])]) -> ClusterQual {
    let mut sorted: Vec<f64> = spectrum.iter().map(|&(distance, _)| distance).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite distances"));
    let knee_cut = spectrum_knee(&sorted);
    let cut = knee_cut.cut;
    let mut excluded = vec![false; spectrum.len()];
    let mut band = 0usize;
    for (index, &(distance, _)) in spectrum.iter().enumerate() {
        if distance <= cut {
            excluded[index] = true;
            band += 1;
        }
    }
    for (_, band_distances) in tied {
        for (index, &distance) in band_distances.iter().enumerate() {
            if distance <= cut {
                excluded[index] = true;
            }
        }
    }
    let mut parent: Vec<usize> = (0..=tied.len()).collect();
    let root = |parent: &mut Vec<usize>, node: usize| {
        let mut current = node;
        while parent[current] != current {
            current = parent[current];
        }
        let mut walked = node;
        while parent[walked] != current {
            let next = parent[walked];
            parent[walked] = current;
            walked = next;
        }
        current
    };
    for (left, (index_left, _)) in tied.iter().enumerate() {
        if spectrum[*index_left].0 <= cut {
            let (a, b) = (root(&mut parent, 0), root(&mut parent, left + 1));
            if a != b {
                parent[a] = b;
            }
        }
        for (right, (index_right, _)) in tied.iter().enumerate().skip(left + 1) {
            if tied[left].1[*index_right] <= cut {
                let (a, b) = (root(&mut parent, left + 1), root(&mut parent, right + 1));
                if a != b {
                    parent[a] = b;
                }
            }
        }
    }
    let mut roots: Vec<usize> = (0..=tied.len()).map(|node| root(&mut parent, node)).collect();
    roots.sort_unstable();
    roots.dedup();
    let k = roots.len();
    let mut alternative: Option<f64> = None;
    for (index, &(_, score)) in spectrum.iter().enumerate() {
        if !excluded[index] {
            alternative = Some(alternative.map_or(score, |value| value.max(score)));
        }
    }
    let shape = if spectrum.is_empty() {
        "single_class"
    } else if knee_cut.has_knee {
        "knee"
    } else {
        "no_knee"
    };
    let p = qual_p(s_win, k, alternative);
    ClusterQual {
        knee: knee_cut.has_knee.then_some(cut),
        cut,
        shape,
        cluster_size: 1 + band,
        k,
        excluded,
        alternative,
        p,
    }
}

/// The canonical diploid pair mixture: P(read | generated from the pair)
/// with the equal 1/2 homolog prior (the balanced 15x/15x sample).
#[inline]
fn mix_logsumexp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY && b == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let m = a.max(b);
    m + (0.5 * (a - m).exp() + 0.5 * (b - m).exp()).ln()
}

/// Inverse of flat = j*(j+1)/2 + i (i <= j).
fn unrank_pair(flat: usize) -> (usize, usize) {
    let mut j = (((8.0 * flat as f64 + 1.0).sqrt() - 1.0) / 2.0).floor() as usize;
    while (j + 1) * (j + 2) / 2 <= flat {
        j += 1;
    }
    while j * (j + 1) / 2 > flat {
        j -= 1;
    }
    let i = flat - j * (j + 1) / 2;
    (i, j)
}

// ---------------------------------------------------------------------------
// The scoring core (the new likelihood layer).
// ---------------------------------------------------------------------------

/// The uniform substitution scoring constants derived from the reads'
/// own quality bytes (epsilon and phred are carried for the receipts).
#[derive(Clone, Copy)]
#[allow(dead_code)]
struct Scoring {
    a: f64,
    b: f64,
    epsilon: f64,
    phred: u8,
}

impl Scoring {
    /// The canonical per-placement closed form: m matching bases and
    /// c mismatching bases. Both the factorized path and the direct
    /// per-base walk produce the integers (m, c); this one expression
    /// turns them into the placement score, so integer equality
    /// implies bit-exact score equality.
    #[inline]
    fn placement(&self, m: u32, c: u32) -> f64 {
        m as f64 * self.a + c as f64 * self.b
    }

    /// The model's derived per-read MINIMUM likelihood: every base of
    /// the read mismatching (READ_LENGTH * ln(eps/3)) — the rigorous
    /// lower bound of the read's likelihood under ANY placement, hence
    /// the floor substituted for a read the row cannot anchor (the
    /// anchored evidence model has no placement there; the exactness
    /// receipt measures the true unanchored floor on the sample). As
    /// a logsumexp term it underflows to +0.0 exactly wherever an
    /// anchored placement exists (150 * ln(eps/3) ~ -1546 against
    /// placement scores ~ -1), so it is a bit-exact no-op there.
    #[inline]
    fn floor(&self) -> f64 {
        READ_LENGTH as f64 * self.b
    }
}

/// One anchored placement of a record's pattern onto one candidate
/// fold: the census orientation that discovered it and the
/// FOLD-RELATIVE pinned window starts per canonical anchor (None =
/// absent from the fold's walk or ambiguous under the monotone rule).
struct Pin {
    orientation: u8,
    pinned: Vec<Option<i64>>,
    /// The merged bp of the pinned anchor windows in hull coordinates
    /// (mirror-invariant: the mirrored window set is the reflection of
    /// the forward set, so the merged coverage length is identical).
    backbone_bp: u64,
}

/// One folded candidate: identical sequence AND walk (relative) —
/// coalesced copies are ONE hypothesis.
struct Fold {
    seq: Vec<u8>,
    len: u64,
    /// The representative row's walk, relative to its start (the same
    /// walk the committed partition GFAs' P lines spell).
    walk: Vec<(u64, i32)>,
    /// node abs id -> relative (bp, sign) occurrences.
    node_positions: HashMap<u32, Vec<(u64, i32)>>,
    members: Vec<MemberRow>,
    /// Contained-step node usage multiset (the QUAL spectrum's
    /// material signature), sorted by key.
    nodes: Vec<(u64, u32)>,
    edges: Vec<(u64, u32)>,
}

impl Fold {
    fn positions_of(&self, node_abs: u32) -> Vec<(u64, i32)> {
        self.node_positions
            .get(&node_abs)
            .cloned()
            .unwrap_or_default()
    }
}

/// One evidence unit: a distinct read variant of one record with its
/// instance count (the record-once discipline: the record's instances
/// are its evidence, each instance voting once per row).
struct Unit {
    record: usize,
    variant: usize,
    mirror: bool,
    read: usize,
    w_lo: u64,
    count: u64,
}

/// The bracketing pinned anchor for one hull-relative skipped
/// position d (slice A's committed serving rule: the pinned anchor
/// with the greatest own position <= d, else the least).
fn serving_anchor(pinned: &[Option<i64>], own_hull: &[u64], d: i64) -> Option<(usize, i64)> {
    let mut best: Option<(usize, i64)> = None;
    for (index, &r) in pinned.iter().enumerate() {
        let Some(_r) = r else { continue };
        let own_pos = own_hull[index] as i64;
        let better = match best {
            None => true,
            Some((_, best_pos)) => {
                let (own_left, best_left) = (own_pos <= d, best_pos <= d);
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
}

/// The FACTORIZED per-placement score for one read variant: the
/// backbone (merged pinned-window bp, matches by node identity and
/// sign — asserted per base by the direct check) plus the skipped-base
/// votes against the fold's spelled sequence. Returns None when the
/// placement is invalid for this variant (a voted base projects
/// outside the fold's extent — the uniform placement prior's support
/// requires the read fully inside).
#[allow(clippy::too_many_arguments)]
fn placement_score(
    unit: &Unit,
    canonical: &[(i32, u64)],
    pin: &Pin,
    fold_seq: &[u8],
    fold_len: u64,
    reads: &[Vec<u8>],
    k: u64,
) -> Option<(u32, u32)> {
    let read = &reads[unit.read];
    let forward = physical_forward(pin.orientation, unit.mirror);
    let own_positions = own_positions_of(canonical, unit.mirror, unit.w_lo);
    let own_hull = own_positions_of(canonical, unit.mirror, 0);
    let skipped = skipped_positions(&own_positions, k, read.len() as u64);
    let mut matches = pin.backbone_bp as u32;
    let mut mismatches = 0u32;
    for &(lo, hi) in &skipped {
        for i in lo..hi {
            let d = i as i64 - unit.w_lo as i64;
            let (index, own_pos) = serving_anchor(&pin.pinned, &own_hull, d)?;
            let r = pin.pinned[index].unwrap();
            let coord = if forward {
                r + (d - own_pos)
            } else {
                r + (own_pos + k as i64 - 1 - d)
            };
            if coord < 0 || coord >= fold_len as i64 {
                return None;
            }
            let base = fold_seq[coord as usize];
            let read_base = if forward {
                read[i as usize]
            } else {
                complement(read[i as usize])
            };
            if read_base == base {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }
    Some((matches, mismatches))
}

/// The DIRECT per-base recomputation of one placement's score for one
/// variant (the brute-force side of the exactness proof, computed
/// in-process over the whole pilot domain): the backbone is verified
/// base-by-base against the fold's sequence instead of trusted from
/// node identity, and the skipped votes are recomputed. Integer
/// equality with the factorized (m, c) is the exactness assertion.
#[allow(clippy::too_many_arguments)]
fn direct_score(
    unit: &Unit,
    canonical: &[(i32, u64)],
    pin: &Pin,
    fold_seq: &[u8],
    fold_len: u64,
    reads: &[Vec<u8>],
    k: u64,
) -> io::Result<(u32, u32)> {
    let read = &reads[unit.read];
    let forward = physical_forward(pin.orientation, unit.mirror);
    let own_positions = own_positions_of(canonical, unit.mirror, unit.w_lo);
    let own_hull = own_positions_of(canonical, unit.mirror, 0);
    let mut matches = 0u32;
    let mut mismatches = 0u32;
    // The backbone, per base over the MERGED pinned-window coverage:
    // syncmer windows overlap heavily (k=63, spacing <= w), so the
    // per-anchor loops would count the overlap bases several times.
    // Each base is verified through ONE covering pinned anchor (the
    // greatest own position at or before the base whose window covers
    // it) against the fold's sequence by direct comparison.
    let backbone_windows: Vec<(u64, u64)> = pin
        .pinned
        .iter()
        .enumerate()
        .filter_map(|(j, &r)| r.map(|_| (own_hull[j], own_hull[j] + k)))
        .collect();
    for (lo, hi) in merge_windows(backbone_windows) {
        for d in lo..hi {
            let mut cover: Option<usize> = None;
            for (j, &r) in pin.pinned.iter().enumerate() {
                if r.is_none() {
                    continue;
                }
                let p = own_hull[j];
                if p <= d && d < p + k {
                    match cover {
                        None => cover = Some(j),
                        Some(current) if p > own_hull[current] => cover = Some(j),
                        _ => {}
                    }
                }
            }
            let j = cover.expect("merged coverage base without a covering anchor");
            let r = pin.pinned[j].unwrap();
            let t = d - own_hull[j];
            let read_i = unit.w_lo as usize + d as usize;
            let read_base = if forward {
                read[read_i]
            } else {
                complement(read[read_i])
            };
            let coord = if forward {
                r + t as i64
            } else {
                r + k as i64 - 1 - t as i64
            };
            ensure(
                coord >= 0 && (coord as usize) < fold_seq.len(),
                "direct backbone outside the fold",
            )?;
            if read_base == fold_seq[coord as usize] {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }
    // The skipped votes, per base.
    let skipped = skipped_positions(&own_positions, k, read.len() as u64);
    for &(lo, hi) in &skipped {
        for i in lo..hi {
            let d = i as i64 - unit.w_lo as i64;
            let (index, own_pos) = serving_anchor(&pin.pinned, &own_hull, d)
                .ok_or_else(|| invalid("direct vote without a serving anchor"))?;
            let r = pin.pinned[index].unwrap();
            let coord = if forward {
                r + (d - own_pos)
            } else {
                r + (own_pos + k as i64 - 1 - d)
            };
            ensure(
                coord >= 0 && coord < fold_len as i64,
                "direct vote outside the fold",
            )?;
            let base = fold_seq[coord as usize];
            let read_base = if forward {
                read[i as usize]
            } else {
                complement(read[i as usize])
            };
            if read_base == base {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }
    let _ = own_positions;
    Ok((matches, mismatches))
}

/// The per-unit per-fold log-likelihood: the uniform placement prior
/// (both strands, 2 * (len - 149) full-read placements) plus the
/// logsumexp over the valid anchored placements. NEG_INFINITY when no
/// placement is valid (the fold shares no anchor, or every pin leaves
/// the read partially outside the fold — the anchored evidence model
/// has no placement there; the exactness receipt measures what the
/// unanchored floor would be on the sample).
fn unit_log_likelihood(
    scores: &[f64],
    fold_len: u64,
) -> f64 {
    if scores.is_empty() || fold_len < READ_LENGTH as u64 {
        return f64::NEG_INFINITY;
    }
    let prior = -(2.0 * (fold_len - READ_LENGTH as u64 + 1) as f64).ln();
    let best = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = scores.iter().map(|&s| (s - best).exp()).sum();
    prior + best + sum.ln()
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

fn main() -> io::Result<()> {
    let started = Instant::now();
    let options = Options::parse();
    let rss = RssGuard::new(options.rss_budget_gib);
    let loci: Vec<u32> = options
        .loci
        .split(',')
        .map(|s| s.trim().parse::<u32>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| invalid("invalid loci list"))?;
    ensure(!loci.is_empty(), "no loci requested")?;
    {
        let mut sorted = loci.clone();
        sorted.sort_unstable();
        sorted.dedup();
        ensure(sorted == loci, "loci must be distinct and sorted")?;
    }

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

    // The census receipt (slice A's input, verbatim).
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
    // component's window ids in order (slice A's derivation).
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
    for &locus in &loci {
        ensure((locus as usize) < n_windows, "locus outside the window span")?;
    }
    eprintln!(
        "[score] inputs: {} census records, k {k}, {} windows, loci {loci:?}",
        census_records.len(),
        n_windows,
    );

    // ------------------------------------------- the scoring rate derivation
    let quality_started = Instant::now();
    let mut quality_bytes: BTreeSet<u8> = BTreeSet::new();
    let mut total_reads = 0u64;
    stream_fastq(&options.reads, |_seq, qual| {
        total_reads += 1;
        for &q in qual {
            quality_bytes.insert(q);
        }
        Ok(())
    })?;
    ensure(
        quality_bytes.len() == 1,
        "the sample does not carry uniform read qualities",
    )?;
    let phred = *quality_bytes.iter().next().unwrap() - b'!';
    let epsilon = 10f64.powf(-(phred as f64) / 10.0);
    let scoring = Scoring {
        a: (1.0 - epsilon).ln(),
        b: (epsilon / 3.0).ln(),
        epsilon,
        phred,
    };
    eprintln!(
        "[score] scoring rate: uniform Phred {phred} over {total_reads} reads \
         (epsilon {epsilon}, A {}, B {}) [{:.1}s]",
        scoring.a,
        scoring.b,
        quality_started.elapsed().as_secs_f64(),
    );
    rss.probe("quality")?;

    // ------------------------------------- the read-record re-derivation
    let derive_started = Instant::now();
    ensure(options.derive_cache.exists(), "derive cache absent")?;
    let (key_tokens, key_reads, reads, read_multiplicity, raw_read_records): (
        Vec<Vec<u64>>,
        Vec<Vec<usize>>,
        Vec<Vec<u8>>,
        Vec<u32>,
        Vec<Vec<(u32, Vec<(i32, u32)>)>>,
    ) = {
        let (tokens, lists, seqs, mults, records) = read_derive_cache(&options.derive_cache)?;
        ensure(
            seqs.iter().all(|seq| seq.len() == READ_LENGTH),
            "cached read is not L150",
        )?;
        (tokens, lists, seqs, mults, records)
    };
    eprintln!(
        "[score] derive cache loaded: {} reads, {} keys ({:.1}s)",
        reads.len(),
        key_tokens.len(),
        derive_started.elapsed().as_secs_f64(),
    );
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
    rss.probe("derivation")?;

    // The routing's own canonical-scheme territory steps around one
    // panel position (ported from build_territory_index_rows, slice
    // A's binding lesson: the census occurrence starts live in THIS
    // convention, not the stored path walk's single frame).
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

    // ------------------------- the pilot record set and the binding
    // The locus's evidence set: records with at least one occurrence
    // TOUCHING a scored window (the routing's own territory-touch rule).
    let mut pilot_records: BTreeSet<usize> = BTreeSet::new();
    for (record, line) in census_records.iter().enumerate() {
        for occ in &line.occurrences {
            if occ.partitions.iter().any(|w| loci.contains(w)) {
                pilot_records.insert(record);
                break;
            }
        }
    }
    eprintln!("[score] pilot records: {}", pilot_records.len());

    // Bind every pilot record to its derived key (slice A's binding,
    // verbatim, restricted to the pilot set).
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
    struct KeyShape {
        canonical: Vec<(i32, u64)>,
        span: u64,
        read: usize,
        mirrored: bool,
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
    // Sequence verification of one key at one occurrence (slice A,
    // verbatim): every anchor window must sit at its projected panel
    // position and equal the representative read's own k-mer.
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
    // The anchor window start under one placement (slice A, verbatim).
    fn anchor_window_start(canonical_pos: u64, start: u64, span: u64, k: u64, orientation: u8) -> u64 {
        if orientation == 0 {
            start + canonical_pos
        } else {
            start + span - k - canonical_pos
        }
    }
    let mut census_key: Vec<Option<u32>> = vec![None; census_records.len()];
    let mut record_shifts: HashMap<usize, Vec<i64>> = HashMap::new();
    for &record in &pilot_records {
        let line = &census_records[record];
        let occ = line
            .occurrences
            .first()
            .ok_or_else(|| invalid("census record has no occurrences"))?;
        let steps = canonical_steps_near(occ.path, occ.start)?
            .into_iter()
            .filter(|&(bp, _)| bp == occ.start)
            .map(|(_, node)| node)
            .collect::<Vec<i32>>();
        ensure(!steps.is_empty(), "census occurrence start has no panel step")?;
        let mut candidates: BTreeSet<u32> = BTreeSet::new();
        for &step in &steps {
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
        let mut surviving: Vec<(u32, Vec<i64>)> = Vec::new();
        for key in candidates {
            if let Some(per_occurrence) = verify_all(key)? {
                surviving.push((key, per_occurrence));
            }
        }
        ensure(
            surviving.len() == 1,
            "census record has no unique verifying key",
        )?;
        census_key[record] = Some(surviving[0].0);
        record_shifts.insert(record, surviving.swap_remove(0).1);
    }
    eprintln!(
        "[score] pilot records bound: {} ({:.1}s)",
        pilot_records.len(),
        started.elapsed().as_secs_f64(),
    );
    rss.probe("binding")?;

    // ----------------------- per-record read variants (slice A, verbatim)
    struct ReadVariant {
        mirror: bool,
        read: usize,
        w_lo: u64,
        count: u64,
    }
    let mut record_variants: HashMap<usize, Vec<ReadVariant>> = HashMap::new();
    for &record in &pilot_records {
        let line = &census_records[record];
        let key = census_key[record].unwrap();
        let mut variants: BTreeMap<(bool, u64, u64), ReadVariant> = BTreeMap::new();
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
        record_variants.insert(record, variants);
    }
    rss.probe("variants")?;

    // ------------------------------------------------------------------ receipts
    let mut report = BufWriter::new(File::create(&options.out)?);
    let mut exactness = BufWriter::new(File::create(&options.exactness_out)?);
    let ingredients_path = format!("{}.ingredients.jsonl", options.out.display());
    let mut ingredients = BufWriter::new(File::create(&ingredients_path)?);
    let records_path = format!("{}.records.jsonl", options.out.display());
    let mut records_file = BufWriter::new(File::create(&records_path)?);

    let mut store = RowStore::new(&panel, &path_of_name, &fetch_seq);
    let contig = options
        .component
        .splitn(3, '#')
        .nth(2)
        .ok_or_else(|| invalid("component lacks a contig suffix"))?
        .to_string();
    // The axis partitions of every window the pilot records touch
    // (the locality classification needs their rows).
    for &partition in &window_partition {
        let map = maps
            .get(&partition)
            .ok_or_else(|| invalid("axis partition map missing"))?;
        store.ensure_partition(&map.members, partition)?;
    }
    rss.probe("rows")?;

    let mut total_factorized_checked = 0u64;
    let mut total_factorized_equal = 0u64;

    for &locus in &loci {
        let locus_started = Instant::now();
        let window = locus as usize;
        let partition = window_partition[window];
        let map = &maps[&partition];

        // ------------------------------------------------------- the folds
        // Fold the partition's member rows by identical sequence AND
        // walk (coalesced copies are ONE hypothesis).
        let fold_started = Instant::now();
        let mut fold_map: BTreeMap<(Vec<u8>, Vec<(u64, i32)>), usize> = BTreeMap::new();
        let mut folds: Vec<Fold> = Vec::new();
        for member in &map.members {
            let row_index = *store
                .dedup
                .get(&{
                    let path_idx = path_of_name[&member.path_name];
                    (path_idx, member.start, member.end.min(panel.name_map.path_to_length[path_idx]))
                })
                .ok_or_else(|| invalid("member row missing from the store"))?;
            let row = &store.folds[row_index];
            let relative: Vec<(u64, i32)> = row
                .walk
                .iter()
                .map(|&(bp, node)| (bp - row.start, node))
                .collect();
            let key = (row.seq.clone(), relative);
            let index = *fold_map.entry(key).or_insert_with(|| {
                let mut node_positions: HashMap<u32, Vec<(u64, i32)>> = HashMap::new();
                for &(bp, node) in &row.walk {
                    node_positions
                        .entry(node.unsigned_abs())
                        .or_default()
                        .push((bp - row.start, node.signum()));
                }
                let contained: Vec<(u64, i32)> = row
                    .walk
                    .iter()
                    .copied()
                    .filter(|&(bp, _)| bp >= row.start && bp + k <= row.end)
                    .map(|(bp, node)| (bp - row.start, node))
                    .collect();
                let mut nodes: Vec<(u64, u32)> = contained
                    .iter()
                    .map(|&(_, node)| (node.unsigned_abs() as u64, 1))
                    .collect();
                nodes.sort_unstable();
                nodes.dedup_by(|a, b| {
                    if a.0 == b.0 {
                        b.1 += a.1;
                        true
                    } else {
                        false
                    }
                });
                let mut edges: Vec<(u64, u32)> = contained
                    .windows(2)
                    .map(|w| (pack_edge(w[0].1.unsigned_abs(), w[1].1.unsigned_abs()), 1))
                    .collect();
                edges.sort_unstable();
                edges.dedup_by(|a, b| {
                    if a.0 == b.0 {
                        b.1 += a.1;
                        true
                    } else {
                        false
                    }
                });
                folds.push(Fold {
                    seq: row.seq.clone(),
                    len: (row.end - row.start) as u64,
                    walk: row
                        .walk
                        .iter()
                        .map(|&(bp, node)| (bp - row.start, node))
                        .collect(),
                    node_positions,
                    members: Vec::new(),
                    nodes,
                    edges,
                });
                folds.len() - 1
            });
            folds[index].members.push(member.clone());
        }
        let n_folds = folds.len();
        eprintln!(
            "[score] locus {locus}: partition {partition}, {} members -> {n_folds} folds [{:.1}s]",
            map.members.len(),
            fold_started.elapsed().as_secs_f64(),
        );

        // ---------------------------------- the locus's records and units
        let mut locus_records: Vec<usize> = Vec::new();
        for &record in &pilot_records {
            if census_records[record]
                .occurrences
                .iter()
                .any(|occ| occ.partitions.contains(&locus))
            {
                locus_records.push(record);
            }
        }
        // Per record: the touching occurrences' orientations, covered
        // keys, and the locality classification.
        struct TouchingData {
            orientations: BTreeSet<u8>,
            covered_nodes: BTreeSet<u64>,
            covered_edges: BTreeSet<u64>,
            in_axis: u64,
            overhang: u64,
            extension: u64,
        }
        let universe_nodes: BTreeSet<u64> =
            folds.iter().flat_map(|f| f.nodes.iter().map(|n| n.0)).collect();
        let universe_edges: BTreeSet<u64> =
            folds.iter().flat_map(|f| f.edges.iter().map(|e| e.0)).collect();
        let mut touching: HashMap<usize, TouchingData> = HashMap::new();
        let mut locality_examples: Vec<serde_json::Value> = Vec::new();
        for &record in &locus_records {
            let line = &census_records[record];
            let span = walk_span(
                &key_shape(
                    census_key[record].unwrap(),
                    &key_tokens,
                    &key_reads,
                    &read_records,
                    k,
                )?
                .canonical,
                k,
            );
            let shifts = &record_shifts[&record];
            let mut data = TouchingData {
                orientations: BTreeSet::new(),
                covered_nodes: BTreeSet::new(),
                covered_edges: BTreeSet::new(),
                in_axis: 0,
                overhang: 0,
                extension: 0,
            };
            for (position, occ) in line.occurrences.iter().enumerate() {
                if !occ.partitions.contains(&locus) {
                    continue;
                }
                data.orientations.insert(occ.orientation);
                for interval in &occ.intervals {
                    for &node in interval {
                        if universe_nodes.contains(&(node as u64)) {
                            data.covered_nodes.insert(node as u64);
                        }
                    }
                    for pair in interval.windows(2) {
                        let edge = pack_edge(pair[0], pair[1]);
                        if universe_edges.contains(&edge) {
                            data.covered_edges.insert(edge);
                        }
                    }
                }
                // The locality classification (slice A's own-row rule):
                // containment in a member row of any touched window's
                // axis partition, on the occurrence's own path.
                let origin = (occ.start as i64 + shifts[position]).max(0) as u64;
                let mut contained = false;
                let mut overlapped = false;
                for &w in &occ.partitions {
                    if (w as usize) >= n_windows {
                        continue;
                    }
                    for &row_index in
                        store.rows_of_partition_path(window_partition[w as usize], occ.path)
                    {
                        let row = &store.folds[row_index];
                        if row.start <= origin && origin + span <= row.end {
                            contained = true;
                        }
                        if row.start < origin + span && origin < row.end {
                            overlapped = true;
                        }
                    }
                }
                if contained {
                    data.in_axis += 1;
                } else if overlapped {
                    data.overhang += 1;
                    if locality_examples.len() < 8 {
                        locality_examples.push(json!({
                            "record": record, "occ": position,
                            "path": occ.path, "origin": origin, "span": span,
                            "class": "overhang",
                        }));
                    }
                } else {
                    data.extension += 1;
                    if locality_examples.len() < 8 {
                        locality_examples.push(json!({
                            "record": record, "occ": position,
                            "path": occ.path, "origin": origin, "span": span,
                            "class": "extension",
                        }));
                    }
                }
            }
            touching.insert(record, data);
        }
        // The evidence units: (record, variant).
        let mut units: Vec<Unit> = Vec::new();
        for &record in &locus_records {
            for (variant, v) in record_variants[&record].iter().enumerate() {
                units.push(Unit {
                    record,
                    variant,
                    mirror: v.mirror,
                    read: v.read,
                    w_lo: v.w_lo,
                    count: v.count,
                });
            }
        }
        let n_units = units.len();
        eprintln!(
            "[score] locus {locus}: {} records, {n_units} units",
            locus_records.len(),
        );

        // ------------------------------- the scoring matrix (parallel folds)
        let scoring_started = Instant::now();
        struct FoldScore {
            ll: Vec<f64>,
            pinned: Vec<bool>,
            mismatches: Vec<u32>,
            factorized_checked: u64,
            factorized_equal: u64,
            edge_discarded: u64,
            no_pin_records: u64,
            abstained_hull_bp: u64,
        }
        let scored: Vec<FoldScore> = (0..n_folds)
            .into_par_iter()
            .map(|fold_index| {
                let fold = &folds[fold_index];
                let mut ll = Vec::with_capacity(n_units);
                let mut pinned: Vec<bool> = Vec::with_capacity(n_units);
                let mut mismatches = Vec::with_capacity(n_units);
                let mut factorized_checked = 0u64;
                let mut factorized_equal = 0u64;
                let mut edge_discarded = 0u64;
                let mut no_pin_records = 0u64;
                let mut abstained_hull_bp = 0u64;
                let positions_of = |node: u32| fold.positions_of(node);
                for unit in &units {
                    let record = unit.record;
                    let canonical = decode_record_tokens(&key_tokens[census_key[record].unwrap() as usize])
                        .expect("bound record");
                    // The pins per orientation present among the
                    // touching occurrences (fold-independent — the
                    // convention cannot bias fold comparisons).
                    let mut pins: Vec<Option<Pin>> = vec![None, None];
                    let mut any_pin = false;
                    for &orientation in &touching[&record].orientations {
                        let (pinned, _ambiguous) =
                            anchor_correspondence(&canonical, orientation, &positions_of);
                        if pinned.iter().all(|p| p.is_none()) {
                            continue;
                        }
                        // The merged pinned-window backbone (hull
                        // coordinates, mirror-invariant).
                        let own_flat = own_positions_of(&canonical, false, 0);
                        let windows: Vec<(u64, u64)> = pinned
                            .iter()
                            .enumerate()
                            .filter_map(|(j, &p)| p.map(|r| (own_flat[j], r as u64)))
                            .map(|(p, _)| (p, p + k))
                            .collect();
                        let backbone_bp: u64 =
                            merge_windows(windows).iter().map(|&(lo, hi)| hi - lo).sum();
                        let windows_inside = pinned.iter().flatten().all(|&r| {
                            r + k <= fold.len
                        });
                        // The abstained hull bp (the named remainder:
                        // read hull bases inside unpinned anchor
                        // windows).
                        let hull_windows: Vec<(u64, u64)> =
                            own_flat.iter().map(|&p| (p, p + k)).collect();
                        let hull_covered: u64 = merge_windows(hull_windows)
                            .iter()
                            .map(|&(lo, hi)| hi - lo)
                            .sum();
                        abstained_hull_bp += hull_covered - backbone_bp;
                        if !windows_inside {
                            continue;
                        }
                        any_pin = true;
                        let pinned_i: Vec<Option<i64>> =
                            pinned.iter().map(|p| p.map(|v| v as i64)).collect();
                        pins[orientation as usize] = Some(Pin {
                            orientation,
                            pinned: pinned_i,
                            backbone_bp,
                        });
                    }
                    if !any_pin {
                        // No anchored placement: the model's derived
                        // per-read minimum (the all-mismatch floor).
                        no_pin_records += 1;
                        ll.push(unit_log_likelihood(&[scoring.floor()], fold.len));
                        pinned.push(false);
                        mismatches.push(0);
                        continue;
                    }
                    let mut scores: Vec<f64> = vec![scoring.floor()];
                    let mut unit_mismatches = 0u32;
                    let mut pinned_here = false;
                    for pin in pins.iter().flatten() {
                        let Some((m, c)) = placement_score(
                            unit, &canonical, pin, &fold.seq, fold.len, &reads, k,
                        ) else {
                            edge_discarded += 1;
                            continue;
                        };
                        let (dm, dc) = direct_score(
                            unit, &canonical, pin, &fold.seq, fold.len, &reads, k,
                        )
                        .expect("direct score");
                        factorized_checked += 1;
                        if (m, c) == (dm, dc) {
                            factorized_equal += 1;
                        }
                        unit_mismatches = unit_mismatches.max(c);
                        pinned_here = true;
                        scores.push(scoring.placement(m, c));
                    }
                    pinned.push(pinned_here);
                    ll.push(unit_log_likelihood(&scores, fold.len));
                    mismatches.push(unit_mismatches);
                }
                FoldScore {
                    ll,
                    pinned,
                    mismatches,
                    factorized_checked,
                    factorized_equal,
                    edge_discarded,
                    no_pin_records,
                    abstained_hull_bp,
                }
            })
            .collect();
        let scoring_seconds = scoring_started.elapsed().as_secs_f64();
        let matrix: Vec<Vec<f64>> = scored.iter().map(|s| s.ll.clone()).collect();
        let pinned: Vec<Vec<bool>> = scored.iter().map(|s| s.pinned.clone()).collect();
        let mismatches: Vec<Vec<u32>> = scored.iter().map(|s| s.mismatches.clone()).collect();
        let factorized_checked: u64 = scored.iter().map(|s| s.factorized_checked).sum();
        let factorized_equal: u64 = scored.iter().map(|s| s.factorized_equal).sum();
        let edge_discarded: u64 = scored.iter().map(|s| s.edge_discarded).sum();
        let no_pin: u64 = scored.iter().map(|s| s.no_pin_records).sum();
        let abstained_hull_bp: u64 = scored.iter().map(|s| s.abstained_hull_bp).sum();
        total_factorized_checked += factorized_checked;
        total_factorized_equal += factorized_equal;
        ensure(
            factorized_checked == factorized_equal,
            &format!(
                "FACTORIZATION GATE FAILED at locus {locus}: \
                 {factorized_equal} of {factorized_checked} equal"
            ),
        )?;
        eprintln!(
            "[score] locus {locus}: factorized {factorized_equal}/{factorized_checked} exact, \
             edge-discarded placements {edge_discarded}, no-pin (unit,fold) {no_pin}, \
             abstained hull bp {abstained_hull_bp} [{scoring_seconds:.1}s]",
        );
        rss.probe(&format!("locus_{locus}_scored"))?;

        // ------------------------------------------- the exhaustive classes
        let class_started = Instant::now();
        let pair_count = n_folds * (n_folds + 1) / 2;
        let unit_counts: Vec<f64> = units.iter().map(|u| u.count as f64).collect();
        let class_lls: Vec<f64> = (0..pair_count)
            .into_par_iter()
            .map(|flat| {
                let (first, second) = unrank_pair(flat);
                let ll_a = &matrix[first];
                let ll_b = &matrix[second];
                let mut total = 0.0f64;
                for (u, &count) in unit_counts.iter().enumerate() {
                    total += count * mix_logsumexp(ll_a[u], ll_b[u]);
                }
                total
            })
            .collect();
        let class_seconds = class_started.elapsed().as_secs_f64();
        let best_ll = class_lls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        ensure(best_ll.is_finite(), "no class has anchored evidence at this locus")?;
        let called: Vec<(usize, usize)> = (0..pair_count)
            .filter(|&flat| class_lls[flat] == best_ll)
            .map(unrank_pair)
            .collect();
        let winner_pair = called[0];

        // -------------------------------------------------- the truth folds
        let fold_of_path = |name: &str| -> io::Result<Option<usize>> {
            let hits: Vec<usize> = folds
                .iter()
                .enumerate()
                .filter(|(_, fold)| fold.members.iter().any(|m| m.path_name == name))
                .map(|(index, _)| index)
                .collect();
            ensure(hits.len() <= 1, "multiple folds carry the truth path")?;
            Ok(hits.first().copied())
        };
        let truth_first = fold_of_path(&options.component)?;
        let truth_second = fold_of_path(&format!("SK1#0#{contig}"))?;
        let truth_pair = match (truth_first, truth_second) {
            (Some(a), Some(b)) => Some((a.min(b), a.max(b))),
            _ => None,
        };
        let truth_flat = truth_pair.map(|(a, b)| b * (b + 1) / 2 + a);
        let truth_log_likelihood = truth_flat.map(|flat| class_lls[flat]);
        let truth_rank = truth_log_likelihood.map(|value| {
            1 + class_lls
                .iter()
                .filter(|&&ll| ll > value)
                .count()
        });
        let truth_tied_classes = truth_log_likelihood.map(|value| {
            class_lls
                .iter()
                .enumerate()
                .filter(|&(flat, &ll)| ll == value && Some(flat) != truth_flat)
                .count()
        });
        let log_gap = match (best_ll, truth_log_likelihood) {
            (best, Some(truth)) => Some(best - truth),
            _ => None,
        };

        // ------------------------------------------- observed mass (no shares)
        let observed_nodes: HashMap<u64, f64> = {
            let mut map: HashMap<u64, f64> = HashMap::new();
            for &record in &locus_records {
                let multiplicity = census_records[record].multiplicity as f64;
                for &node in &touching[&record].covered_nodes {
                    *map.entry(node).or_default() += multiplicity;
                }
            }
            map
        };
        let observed_edges: HashMap<u64, f64> = {
            let mut map: HashMap<u64, f64> = HashMap::new();
            for &record in &locus_records {
                let multiplicity = census_records[record].multiplicity as f64;
                for &edge in &touching[&record].covered_edges {
                    *map.entry(edge).or_default() += multiplicity;
                }
            }
            map
        };
        let observed_node_fn = |key: u64| observed_nodes.get(&key).copied().unwrap_or(0.0);
        let observed_edge_fn = |key: u64| observed_edges.get(&key).copied().unwrap_or(0.0);
        let observed_node_mass: f64 = observed_nodes.values().sum();
        let observed_edge_mass: f64 = observed_edges.values().sum();

        // ------------------------------------------------- the QUAL cluster
        let qual_started = Instant::now();
        let winner_nodes = {
            let (first, second) = winner_pair;
            merged_multiset(&folds[first].nodes, &folds[second].nodes)
        };
        let winner_edges = {
            let (first, second) = winner_pair;
            merged_multiset(&folds[first].edges, &folds[second].edges)
        };
        let called_score = Some(best_ll);
        let spectrum: Vec<(f64, f64, f64, usize, usize, usize)> =
            match (winner_pair, called_score) {
                ((winner_first, winner_second), Some(score)) => {
                    let mut spectrum: Vec<(f64, f64, f64, usize, usize, usize)> = (0..pair_count)
                        .into_par_iter()
                        .filter_map(|flat| {
                            let (first, second) = unrank_pair(flat);
                            if first == winner_first && second == winner_second {
                                return None;
                            }
                            let ll = class_lls[flat];
                            let nodes =
                                merged_multiset(&folds[first].nodes, &folds[second].nodes);
                            let edges =
                                merged_multiset(&folds[first].edges, &folds[second].edges);
                            let distance = differing_observed_mass(&nodes, &winner_nodes, &observed_node_fn)
                                + differing_observed_mass(&edges, &winner_edges, &observed_edge_fn);
                            let signature =
                                signature_cosine_distance(&nodes, &edges, &winner_nodes, &winner_edges);
                            Some((distance, (ll - score).exp(), signature, first, second, flat))
                        })
                        .collect();
                    spectrum
                        .sort_by(|a, b| a.partial_cmp(b).expect("finite spectrum entries"));
                    spectrum
                }
                _ => Vec::new(),
            };
        let mut tied_bands: Vec<(usize, Vec<f64>)> = Vec::new();
        for &(first, second) in called.iter().skip(1) {
            let index = spectrum
                .iter()
                .position(|&(_, _, _, a, b, _)| a == first && b == second)
                .expect("called class missing from the winner's spectrum");
            let nodes = merged_multiset(&folds[first].nodes, &folds[second].nodes);
            let edges = merged_multiset(&folds[first].edges, &folds[second].edges);
            let band: Vec<f64> = spectrum
                .par_iter()
                .map(|&(_, _, _, other_first, other_second, _)| {
                    let other_nodes =
                        merged_multiset(&folds[other_first].nodes, &folds[other_second].nodes);
                    let other_edges =
                        merged_multiset(&folds[other_first].edges, &folds[other_second].edges);
                    differing_observed_mass(&nodes, &other_nodes, &observed_node_fn)
                        + differing_observed_mass(&edges, &other_edges, &observed_edge_fn)
                })
                .collect();
            tied_bands.push((index, band));
        }
        let tied_refs: Vec<(usize, &[f64])> =
            tied_bands.iter().map(|(index, band)| (*index, band.as_slice())).collect();
        let cluster = called_score
            .filter(|_| !spectrum.is_empty() || !called.is_empty())
            .map(|_| {
                let pairs: Vec<(f64, f64)> = spectrum
                    .iter()
                    .map(|&(distance, score, _, _, _, _)| (distance, score))
                    .collect();
                cluster_form_qual(1.0, &pairs, &tied_refs)
            });
        let cluster_knee = cluster.as_ref().and_then(|state| state.knee);
        let cluster_shape = cluster.as_ref().map(|state| state.shape);
        let cluster_size = cluster.as_ref().map(|state| state.cluster_size);
        let cluster_k = cluster.as_ref().map(|state| state.k);
        let cluster_alternative = cluster.as_ref().and_then(|state| state.alternative);
        let cluster_excluded = cluster
            .as_ref()
            .map(|state| state.excluded.iter().filter(|e| **e).count());
        let confidence = cluster.as_ref().and_then(|state| state.p);
        let qual = confidence.and_then(qual_from_p);
        let qual_unbounded = confidence.is_some_and(|value| value >= 1.0);
        let posterior_logsumexp = called_score.map(|score| {
            let total: f64 = class_lls.iter().map(|&ll| (ll - score).exp()).sum();
            score + total.ln()
        });
        let top_cluster_posterior = match (&cluster, posterior_logsumexp) {
            (Some(state), Some(logsumexp)) => {
                let mut mass = 0.0f64;
                for &(distance, score, _, _, _, _) in &spectrum {
                    if distance <= state.cut {
                        mass += score;
                    }
                }
                Some((1.0 + mass) / (logsumexp - called_score.unwrap_or(0.0)).exp())
            }
            _ => None,
        };
        let truth_in_called_set = truth_flat.map(|flat| {
            let (first, second) = unrank_pair(flat);
            called.contains(&(first, second))
        });
        // The best divergent rival outside the excluded band.
        let best_outside_log_likelihood = cluster
            .as_ref()
            .map(|state| {
                spectrum
                    .iter()
                    .zip(&state.excluded)
                    .filter(|(_, excluded)| !**excluded)
                    .map(|(entry, _)| class_lls[entry.5])
                    .fold(f64::NEG_INFINITY, f64::max)
            })
            .filter(|value| value.is_finite());
        let mut qual_alternative_classes: Vec<serde_json::Value> = Vec::new();
        let alternative_log_gap =
            best_outside_log_likelihood.map(|outside| best_ll - outside);
        if let (Some(state), Some(best_outside)) = (&cluster, best_outside_log_likelihood) {
            for (index, entry) in spectrum.iter().enumerate() {
                let &(distance, score, signature_distance, first, second, flat) = entry;
                if !state.excluded[index] && class_lls[flat] == best_outside {
                    qual_alternative_classes.push(json!({
                        "spectrum_index": index,
                        "distance": distance,
                        "signature_distance": signature_distance,
                        "relative_likelihood": score,
                        "log_likelihood": class_lls[flat],
                        "fold_indices": [first, second],
                    }));
                }
            }
        }
        let nearest_rival_distance = spectrum
            .iter()
            .find(|&&(distance, ..)| distance > 0.0)
            .map(|&(distance, ..)| distance);
        let spectrum_zero_distance = called_score
            .map(|_| spectrum.iter().filter(|&&(distance, ..)| distance == 0.0).count());
        let qual_seconds = qual_started.elapsed().as_secs_f64();

        // ---------------------------------------- the exactness sample
        // Strata: the evenly-spaced valid (unit, fold) pairs; the
        // evenly-spaced pairs WITH mismatch votes (the variant-pocket
        // class); the truth folds; the identical-through-graph fold
        // (the named BTE#3/#4 block28_contig1 pair, when it is a member
        // of this locus's partition). The in-process factorization
        // assertion already covers the WHOLE pilot domain; the sample
        // carries vote-level detail for the independent checker.
        let exactness_started = Instant::now();
        let valid_pairs: Vec<(usize, usize)> = (0..n_folds)
            .flat_map(|f| (0..n_units).map(move |u| (f, u)))
            .filter(|&(f, u)| pinned[f][u])
            .collect();
        let mismatch_pairs: Vec<(usize, usize)> = valid_pairs
            .iter()
            .copied()
            .filter(|&(f, u)| mismatches[f][u] > 0)
            .collect();
        let stride_of = |len: usize, want: u64| -> usize {
            if len == 0 || want == 0 {
                usize::MAX
            } else {
                ((len + want as usize - 1) / want as usize).max(1)
            }
        };
        let mut sampled: BTreeSet<(usize, usize)> = BTreeSet::new();
        let mut stratum_counts: BTreeMap<String, u64> = BTreeMap::new();
        {
            let stride = stride_of(valid_pairs.len(), options.exactness_sample);
            let mut taken = 0u64;
            for (index, &pair) in valid_pairs.iter().enumerate() {
                if index % stride == 0 {
                    sampled.insert(pair);
                    taken += 1;
                }
            }
            stratum_counts.insert("evenly_spaced_valid".to_string(), taken);
        }
        {
            let stride = stride_of(mismatch_pairs.len(), options.exactness_sample);
            let mut taken = 0u64;
            for (index, &pair) in mismatch_pairs.iter().enumerate() {
                if index % stride == 0 {
                    if sampled.insert(pair) {
                        taken += 1;
                    }
                }
            }
            stratum_counts.insert("variant_pocket_mismatch".to_string(), taken);
        }
        let mut named_folds: Vec<(String, usize)> = Vec::new();
        if let Some(first) = truth_first {
            named_folds.push(("truth_first".to_string(), first));
        }
        if let Some(second) = truth_second {
            named_folds.push(("truth_second".to_string(), second));
        }
        if let Some(index) = folds.iter().position(|fold| {
            fold.members
                .iter()
                .any(|m| m.path_name == "BTE#3#block28_contig1")
        }) {
            named_folds.push(("identical_pair".to_string(), index));
        }
        for (name, fold_index) in &named_folds {
            let pairs: Vec<(usize, usize)> = (0..n_units)
                .filter(|&u| pinned[*fold_index][u])
                .map(|u| (*fold_index, u))
                .collect();
            let stride = stride_of(pairs.len(), 25);
            let mut taken = 0u64;
            for (index, pair) in pairs.iter().enumerate() {
                if index % stride == 0 && sampled.insert(*pair) {
                    taken += 1;
                }
            }
            stratum_counts.insert(format!("named_{name}"), taken);
        }
        for &(fold_index, unit_index) in &sampled {
            let unit = &units[unit_index];
            let record = unit.record;
            let canonical =
                decode_record_tokens(&key_tokens[census_key[record].unwrap() as usize])?;
            let fold = &folds[fold_index];
            let positions_of = |node: u32| fold.positions_of(node);
            let mut placements: Vec<serde_json::Value> = Vec::new();
            for &orientation in &touching[&record].orientations {
                let (pinned, _ambiguous) =
                    anchor_correspondence(&canonical, orientation, &positions_of);
                if pinned.iter().all(|p| p.is_none()) {
                    continue;
                }
                let windows_inside =
                    pinned.iter().flatten().all(|&r| r + k <= fold.len);
                if !windows_inside {
                    continue;
                }
                let own_flat = own_positions_of(&canonical, false, 0);
                let windows: Vec<(u64, u64)> = pinned
                    .iter()
                    .enumerate()
                    .filter_map(|(j, &p)| p.map(|r| (own_flat[j], r as u64)))
                    .map(|(p, _)| (p, p + k))
                    .collect();
                let backbone_bp: u64 =
                    merge_windows(windows).iter().map(|&(lo, hi)| hi - lo).sum();
                let pinned_i: Vec<Option<i64>> =
                    pinned.iter().map(|p| p.map(|v| v as i64)).collect();
                let pin = Pin {
                    orientation,
                    pinned: pinned_i,
                    backbone_bp,
                };
                let Some((m, c)) = placement_score(
                    unit, &canonical, &pin, &fold.seq, fold.len, &reads, k,
                ) else {
                    continue;
                };
                // The vote-level detail (the skipped bases, their
                // projected coordinates and their votes).
                let forward = physical_forward(orientation, unit.mirror);
                let own_positions = own_positions_of(&canonical, unit.mirror, unit.w_lo);
                let own_hull = own_positions_of(&canonical, unit.mirror, 0);
                let skipped = skipped_positions(&own_positions, k, READ_LENGTH as u64);
                let mut votes: Vec<serde_json::Value> = Vec::new();
                for &(lo, hi) in &skipped {
                    for i in lo..hi {
                        let d = i as i64 - unit.w_lo as i64;
                        let (index, own_pos) =
                            serving_anchor(&pin.pinned, &own_hull, d).expect("serving anchor");
                        let r = pin.pinned[index].unwrap();
                        let coord = if forward {
                            r + (d - own_pos)
                        } else {
                            r + (own_pos + k as i64 - 1 - d)
                        };
                        let base = fold.seq[coord as usize];
                        let read_base = if forward {
                            reads[unit.read][i as usize]
                        } else {
                            complement(reads[unit.read][i as usize])
                        };
                        votes.push(json!([i, coord, if read_base == base { 1 } else { 0 }]));
                    }
                }
                placements.push(json!({
                    "orientation": orientation,
                    "pinned": pin.pinned,
                    "backbone_bp": backbone_bp,
                    "m": m,
                    "c": c,
                    "score": scoring.placement(m, c),
                    "votes": votes,
                }));
            }
            let read = &reads[unit.read];
            serde_json::to_writer(
                &mut exactness,
                &json!({
                    "locus": locus,
                    "unit": unit_index,
                    "record": record,
                    "variant": unit.variant,
                    "mirror": if unit.mirror { 1 } else { 0 },
                    "w_lo": unit.w_lo,
                    "count": unit.count,
                    "read_fnv": fnv1a64(read),
                    "read": String::from_utf8(read.clone()).unwrap(),
                    "fold": fold_index,
                    "ll": matrix[fold_index][unit_index],
                    "placements": placements,
                }),
            )?;
            writeln!(exactness)?;
        }
        let exactness_seconds = exactness_started.elapsed().as_secs_f64();
        eprintln!(
            "[score] locus {locus}: exactness sample {} pairs of {} valid \
             ({} with mismatch votes) [{exactness_seconds:.1}s]",
            sampled.len(),
            valid_pairs.len(),
            mismatch_pairs.len(),
        );

        // ------------------------------------------------- the named classes
        let mut named: Vec<(usize, usize)> = Vec::new();
        named.push(winner_pair);
        if let Some(flat) = truth_flat {
            let pair = unrank_pair(flat);
            if !named.contains(&pair) {
                named.push(pair);
            }
        }
        if let (Some(state), Some(_)) = (&cluster, best_outside_log_likelihood) {
            for (index, entry) in spectrum.iter().enumerate() {
                if !state.excluded[index] && class_lls[entry.5] == best_outside_log_likelihood.unwrap() {
                    let pair = (entry.3, entry.4);
                    if !named.contains(&pair) {
                        named.push(pair);
                    }
                    break;
                }
            }
        }
        let mut classes_evidence: Vec<serde_json::Value> = Vec::new();
        for &(first, second) in &named {
            let per_unit: Vec<serde_json::Value> = units
                .iter()
                .enumerate()
                .map(|(u, unit)| {
                    json!({
                        "record": unit.record,
                        "variant": unit.variant,
                        "count": unit.count,
                        "fold_ll": [matrix[first][u], matrix[second][u]],
                        "pair_term": mix_logsumexp(matrix[first][u], matrix[second][u]),
                    })
                })
                .collect();
            classes_evidence.push(json!({
                "fold_indices": [first, second],
                "log_likelihood": class_lls[second * (second + 1) / 2 + first],
                "per_unit": per_unit,
            }));
        }
        serde_json::to_writer(&mut records_file, &json!({
            "locus": locus,
            "classes": classes_evidence,
        }))?;
        writeln!(records_file)?;

        // ----------------------------------------------------- the ingredients
        serde_json::to_writer(&mut ingredients, &json!({
            "locus": locus,
            "partition": partition,
            "folds": folds.iter().map(|fold| json!({
                "members": fold.members.iter().map(|m| json!({
                    "path_name": m.path_name, "start": m.start, "end": m.end,
                })).collect::<Vec<_>>(),
                "length": fold.len,
                "sequence": String::from_utf8(fold.seq.clone()).unwrap(),
                "walk": fold.walk.iter().map(|&(bp, node)| json!([bp, node]))
                    .collect::<Vec<_>>(),
                "nodes": fold.nodes,
                "edges": fold.edges,
            })).collect::<Vec<_>>(),
            "units": units.iter().map(|unit| json!({
                "record": unit.record,
                "variant": unit.variant,
                "mirror": if unit.mirror { 1 } else { 0 },
                "w_lo": unit.w_lo,
                "count": unit.count,
                "read": String::from_utf8(reads[unit.read].clone()).unwrap(),
            })).collect::<Vec<_>>(),
            // The per-unit per-fold log-likelihood matrix (folds x
            // units, row-major; null = no valid anchored placement).
            "ll_matrix": matrix,
        }))?;
        writeln!(ingredients)?;

        // ---------------------------------------------------- the main receipt
        let called_classes_json: Vec<serde_json::Value> = called
            .iter()
            .map(|&(first, second)| {
                json!({
                    "fold_indices": [first, second],
                    "physical_pairs": class_physical_pairs(
                        first == second,
                        folds[first].members.len(),
                        folds[second].members.len(),
                    ),
                })
            })
            .collect();
        let fold_identities: Vec<serde_json::Value> = folds
            .iter()
            .map(|fold| {
                json!({
                    "members": fold.members.iter().map(|m| json!({
                        "path_name": m.path_name, "start": m.start, "end": m.end,
                    })).collect::<Vec<_>>(),
                    "length": fold.len,
                })
            })
            .collect();
        let class_fold_pairs: Vec<[usize; 2]> = (0..pair_count)
            .map(|flat| {
                let (first, second) = unrank_pair(flat);
                [first, second]
            })
            .collect();
        let identical_pair_fold = folds.iter().position(|fold| {
            fold.members.iter().any(|m| m.path_name == "BTE#3#block28_contig1")
                && fold.members
                    .iter()
                    .any(|m| m.path_name == "BTE#4#block28_contig1")
        });
        // The log-gap floor attribution: how much of the winner-truth
        // gap is the all-mismatch floor (reads the pair can place on
        // neither homolog) versus per-base evidence on the reads both
        // place. The floor's magnitude is the model's derived minimum,
        // not the true unanchored likelihood (named).
        let pair_places = |i: usize, j: usize, u: usize| pinned[i][u] || pinned[j][u];
        let mut truth_unplaced_units = 0u64;
        let mut truth_unplaced_mass = 0.0f64;
        let mut winner_unplaced_units = 0u64;
        let mut winner_unplaced_mass = 0.0f64;
        let mut both_place_gap = 0.0f64;
        if let (Some((ti, tj)), Some((wi, wj))) = (truth_pair, Some(winner_pair)) {
            for (u, unit) in units.iter().enumerate() {
                let truth_places = pair_places(ti, tj, u);
                let winner_places = pair_places(wi, wj, u);
                if !truth_places {
                    truth_unplaced_units += 1;
                    truth_unplaced_mass += unit.count as f64;
                }
                if !winner_places {
                    winner_unplaced_units += 1;
                    winner_unplaced_mass += unit.count as f64;
                }
                if truth_places && winner_places {
                    both_place_gap += unit.count as f64
                        * (mix_logsumexp(matrix[wi][u], matrix[wj][u])
                            - mix_logsumexp(matrix[ti][u], matrix[tj][u]));
                }
            }
        }
        let locality_totals = {
            let mut in_axis = 0u64;
            let mut overhang = 0u64;
            let mut extension = 0u64;
            for &record in &locus_records {
                let data = &touching[&record];
                in_axis += data.in_axis;
                overhang += data.overhang;
                extension += data.extension;
            }
            (in_axis, overhang, extension)
        };
        let locus_seconds = locus_started.elapsed().as_secs_f64();
        let rss_kb = rss.probe(&format!("locus_{locus}"))?;
        serde_json::to_writer(&mut report, &json!({
            "locus": locus,
            "partition": partition,
            "model": "anchor-realign-v1",
            "scoring": {
                "phred": phred,
                "epsilon": epsilon,
                "match_log_prob": scoring.a,
                "mismatch_log_prob": scoring.b,
                "read_length": READ_LENGTH,
                "uniform_qualities": true,
            },
            "record_count": locus_records.len(),
            "unit_count": n_units,
            "folds": n_folds,
            "physical_members": map.members.len(),
            "class_count": pair_count,
            "class_fold_pairs": class_fold_pairs,
            // null entries = no valid anchored placement for any unit
            // of the class (the anchored evidence model's floor).
            "class_log_likelihoods": class_lls,
            "best_log_likelihood": best_ll,
            "best_fold_indices": [winner_pair.0, winner_pair.1],
            "called_class_count": called.len(),
            "called_classes": called_classes_json,
            "fold_identities": fold_identities,
            "truth_folds": truth_pair.map(|(a, b)| [a, b]),
            "truth_pair_expressible": truth_pair.is_some(),
            "truth_log_likelihood": truth_log_likelihood,
            "truth_rank": truth_rank,
            "truth_tied_classes": truth_tied_classes,
            "truth_in_called_set": truth_in_called_set,
            "log_gap": log_gap,
            "identical_pair_fold": identical_pair_fold,
            "qual": qual,
            "qual_unbounded": qual_unbounded,
            "qual_spectrum_shape": cluster_shape,
            "qual_knee_distance": cluster_knee,
            "qual_cluster_size": cluster_size,
            "qual_cluster_k": cluster_k,
            "qual_spectrum_classes": called_score.map(|_| spectrum.len()),
            "qual_spectrum_zero_distance_classes": spectrum_zero_distance,
            "qual_nearest_rival_distance": nearest_rival_distance,
            "qual_excluded_class_count": cluster_excluded,
            "qual_distance_spectrum": spectrum.iter().map(|s| s.0).collect::<Vec<_>>(),
            "qual_distance_spectrum_scores": spectrum.iter().map(|s| s.1).collect::<Vec<_>>(),
            "qual_signature_distance_spectrum": spectrum.iter().map(|s| s.2).collect::<Vec<_>>(),
            "qual_alternative_similarity": cluster_alternative,
            "qual_alternative_log_gap": alternative_log_gap,
            "qual_alternative_classes": qual_alternative_classes,
            "qual_posterior_logsumexp": posterior_logsumexp,
            "qual_posterior_top_cluster": top_cluster_posterior,
            "observed_node_mass": observed_node_mass,
            "observed_edge_mass": observed_edge_mass,
            "observed_node_keys": observed_nodes.len(),
            "observed_edge_keys": observed_edges.len(),
            "log_gap_decomposition": {
                "truth_unplaced_units": truth_unplaced_units,
                "truth_unplaced_mass": truth_unplaced_mass,
                "winner_unplaced_units": winner_unplaced_units,
                "winner_unplaced_mass": winner_unplaced_mass,
                "both_placed_evidence_gap": both_place_gap,
                "note": "unplaced reads score at the model's derived \
                         all-mismatch minimum (150*ln(eps/3)), the rigorous \
                         lower bound of their unanchored likelihood; the \
                         truth rank is floor-magnitude-robust, the gap \
                         magnitudes are floor-dominated where unplaced \
                         counts differ",
            },
            "locality": {
                "touching_in_axis": locality_totals.0,
                "touching_overhang": locality_totals.1,
                "touching_extension": locality_totals.2,
                "examples": locality_examples,
            },
            "validation": {
                "factorized_checked": factorized_checked,
                "factorized_equal": factorized_equal,
                "edge_discarded_placements": edge_discarded,
                "floor_only_unit_fold": no_pin,
                "abstained_hull_bp": abstained_hull_bp,
                "exactness_sample_pairs": sampled.len(),
                "exactness_valid_pairs": valid_pairs.len(),
                "exactness_mismatch_pairs": mismatch_pairs.len(),
                "exactness_strata": stratum_counts,
            },
            "walls": {
                "folds_seconds": fold_started.elapsed().as_secs_f64() - scoring_seconds - class_seconds - qual_seconds - exactness_seconds,
                "scoring_seconds": scoring_seconds,
                "classes_seconds": class_seconds,
                "qual_seconds": qual_seconds,
                "exactness_seconds": exactness_seconds,
                "locus_seconds": locus_seconds,
            },
            "rss_kb": rss_kb,
        }))?;
        writeln!(report)?;
        eprintln!(
            "[score] locus {locus} done: {n_folds} folds, {n_units} units, {pair_count} classes, \
             truth_rank {:?}, log_gap {log_gap:?}, qual {qual:?} [{locus_seconds:.1}s]",
            truth_rank,
        );
    }

    report.flush()?;
    exactness.flush()?;
    ingredients.flush()?;
    records_file.flush()?;
    ensure(
        total_factorized_checked == total_factorized_equal,
        "the factorization gate failed somewhere",
    )?;
    eprintln!(
        "[score] complete: loci {loci:?}, factorized placements {total_factorized_equal}/\
         {total_factorized_checked} exact, total {:.1}s",
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// Unit tests: the factorized exact form against the direct per-base walk,
// the strand symmetry, the projection equations, and the machinery.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const K: u64 = 8;

    fn fold_seq() -> Vec<u8> {
        // 60 bp, deterministic and non-repetitive enough for the tests.
        (0..60)
            .map(|i| match i % 5 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b"ACGT"[(i / 5) % 4],
            })
            .collect()
    }

    /// The canonical pattern: anchors (node 5 @ 0), (node 7 @ 10);
    /// span 18; the fold carries node 5 forward at rel 4 and node 7
    /// forward at rel 14. The forward read F (30 bp, w_lo 2) spells the
    /// fold from rel 2 to 32 with ONE deliberate mismatch in the
    /// 2 bp interior hull gap.
    fn forward_read(seq: &[u8]) -> (Vec<u8>, usize) {
        let mut read = vec![0u8; 30];
        for i in 0..30usize {
            read[i] = seq[i + 2];
        }
        read[11] = match read[11] {
            b'A' => b'T',
            other => b'A' + (other == b'A') as u8,
        };
        (read, 13) // the mismatching fold coordinate (rel 13)
    }

    fn pin() -> Pin {
        Pin {
            orientation: 0,
            pinned: vec![Some(4), Some(14)],
            backbone_bp: 16,
        }
    }

    fn canonical() -> Vec<(i32, u64)> {
        vec![(5, 0), (7, 10)]
    }


    #[test]
    fn factorized_equals_direct_on_the_forward_placement() {
        let seq = fold_seq();
        let (read, mismatch_coord) = forward_read(&seq);
        let reads = vec![read];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: false,
            read: 0,
            w_lo: 2,
            count: 1,
        };
        let pin = pin();
        let (m, c) = placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K)
            .expect("valid placement");
        // Backbone 16 matches + 14 skipped votes, one of them the
        // deliberate mismatch at fold rel 13.
        assert_eq!((m, c), (29, 1));
        let (dm, dc) = direct_score(&unit, &canonical(), &pin, &seq, 60, &reads, K).unwrap();
        assert_eq!((m, c), (dm, dc));
        // The mismatch is where we placed it.
        let canonical = canonical();
        let own_positions = own_positions_of(&canonical, false, 2);
        let skipped = skipped_positions(&own_positions, K, 30);
        assert_eq!(skipped, vec![(0, 2), (10, 12), (20, 30)]);
        let fold_base = seq[mismatch_coord];
        let read_base = reads[0][11];
        assert_ne!(read_base, fold_base);
    }

    #[test]
    fn the_reverse_strand_placement_is_strand_symmetric() {
        // M = revcomp(F): the same physical read; its own orientation
        // extraction sees the pattern mirrored (w_lo 10), and at the
        // same orientation-0 pin the placement is physically reverse.
        // The per-base evidence must be IDENTICAL: (m, c) equal to the
        // forward case, every voted base the same physical comparison.
        let seq = fold_seq();
        let (forward, _) = forward_read(&seq);
        let mirrored = revcomp(&forward);
        let reads = vec![forward, mirrored];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: true,
            read: 1,
            w_lo: 10,
            count: 1,
        };
        // Mirror own positions: canonical anchor 0 at hull 10, anchor 1
        // at hull 0 — the correspondence is canonical-order-based and
        // mirror-blind, so the pin is the same.
        let own = own_positions_of(&canonical(), true, 10);
        assert_eq!(own, vec![20, 10]);
        let pin = pin();
        let (m, c) = placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K)
            .expect("valid placement");
        assert_eq!((m, c), (29, 1));
        let (dm, dc) = direct_score(&unit, &canonical(), &pin, &seq, 60, &reads, K).unwrap();
        assert_eq!((m, c), (dm, dc));
    }

    #[test]
    fn placements_with_out_of_extent_votes_are_invalid() {
        // A fold whose extent ends inside the read's right flank: the
        // voted bases do not all project inside, so the placement is
        // not a placement of the model (the read must lie fully
        // inside).
        let seq = fold_seq();
        let (read, _) = forward_read(&seq);
        let reads = vec![read];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: false,
            read: 0,
            w_lo: 2,
            count: 1,
        };
        let pin = Pin {
            orientation: 0,
            pinned: vec![Some(4), Some(14)],
            backbone_bp: 16,
        };
        let placed = placement_score(&unit, &canonical(), &pin, &seq, 24, &reads, K);
        assert!(placed.is_none());
        // Backbone windows partly outside the extent are not
        // placements either (the windows_inside check at pin
        // construction discards them).
        let pin_outside = [Some(196i64), Some(206)];
        assert!(pin_outside.iter().flatten().any(|&r| r + K as i64 > 24));
    }

    #[test]
    fn the_backbone_is_mirror_invariant() {
        // The merged pinned-window coverage in hull coordinates has
        // the same total length under both mirror states (the mirrored
        // window set is the reflection of the forward set).
        let canonical = vec![(5, 0), (7, 10), (-9, 14)];
        for mirror in [false, true] {
            let own = own_positions_of(&canonical, mirror, 0);
            let windows: Vec<(u64, u64)> = own.iter().map(|&p| (p, p + K)).collect();
            let total: u64 = merge_windows(windows).iter().map(|&(lo, hi)| hi - lo).sum();
            assert_eq!(total, 20); // windows [0,8),[10,18),[14,22) -> merged [0,8)+[10,22)
        }
    }

    #[test]
    fn unit_likelihood_is_prior_plus_logsumexp() {
        // No placements: the anchored floor.
        assert_eq!(unit_log_likelihood(&[], 200), f64::NEG_INFINITY);
        assert_eq!(unit_log_likelihood(&[0.0], 100), f64::NEG_INFINITY);
        let a = -3.0f64;
        let b = -5.0f64;
        let ll = unit_log_likelihood(&[a, b], 200);
        let expect = -(2.0 * (200 - 149) as f64).ln() + (a.exp() + b.exp()).ln();
        assert!((ll - expect).abs() < 1e-12);
        // Bit-exact when the placements are equal (the identical-fold
        // case): the sum is the placement plus ln 2.
        let ll = unit_log_likelihood(&[a, a], 200);
        let expect = -(2.0 * (200 - 149) as f64).ln() + (2.0 * a.exp()).ln();
        assert!((ll - expect).abs() < 1e-12);
    }

    #[test]
    fn the_serving_anchor_brackets_the_skipped_position() {
        // Pinned anchors at hull positions 0 and 10: the left flank is
        // served by the LEAST, the interior gap and the right flank by
        // the GREATEST pinned anchor at or before the base.
        let pinned = vec![Some(4i64), Some(14)];
        let own_hull = vec![0u64, 10];
        // No pinned anchor at or before d: the LEAST serves.
        assert_eq!(serving_anchor(&pinned, &own_hull, -2), Some((0, 0)));
        assert_eq!(serving_anchor(&pinned, &own_hull, 0), Some((0, 0)));
        assert_eq!(serving_anchor(&pinned, &own_hull, 8), Some((0, 0)));
        assert_eq!(serving_anchor(&pinned, &own_hull, 10), Some((1, 10)));
        assert_eq!(serving_anchor(&pinned, &own_hull, 25), Some((1, 10)));
    }

    #[test]
    fn unrank_pair_is_the_inverse_of_flat() {
        for i in 0..40usize {
            for j in i..40usize {
                let flat = j * (j + 1) / 2 + i;
                assert_eq!(unrank_pair(flat), (i, j));
            }
        }
    }

    #[test]
    fn the_diploid_mixture_handles_the_floor() {
        assert_eq!(mix_logsumexp(f64::NEG_INFINITY, f64::NEG_INFINITY), f64::NEG_INFINITY);
        let a = -2.0f64;
        // The homozygous class: the mixture of a fold with itself is
        // the fold's own likelihood (0.5 P + 0.5 P = P).
        assert!((mix_logsumexp(a, a) - a).abs() < 1e-12);
        assert!((mix_logsumexp(a, f64::NEG_INFINITY) - (a - (2f64).ln())).abs() < 1e-12);
        let mixed = mix_logsumexp(-2.0, -3.0);
        assert!((mixed - ((-2.0f64).exp() * 0.5 + (-3.0f64).exp() * 0.5).ln()).abs() < 1e-12);
    }

    #[test]
    fn the_knee_cut_is_the_max_relative_jump_above_the_typical() {
        // A spectrum with a clean knee: distances 1,1,2,2,3,3,90.
        let sorted = vec![1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 90.0];
        let knee = spectrum_knee(&sorted);
        assert!(knee.has_knee);
        assert_eq!(knee.cut, 3.0);
        // No knee in a geometric spectrum (every relative jump equal
        // to the typical): the derived no-knee case.
        let geometric = vec![1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0];
        assert!(!spectrum_knee(&geometric).has_knee);
    }

    #[test]
    fn qual_is_the_cluster_form_over_the_excluded_band() {
        // One divergent rival far outside a tight band: the knee cuts
        // the band, k = 2, p = 1/(2 + a).
        let spectrum: Vec<(f64, f64)> = vec![
            (0.0, 1.0),
            (1.0, 0.5),
            (2.0, 0.1),
            (100.0, 1e-40),
        ];
        let state = cluster_form_qual(1.0, &spectrum, &[]);
        assert!(state.knee.is_some());
        // A single called class: k = 1, and p = 1/(1 + a) is capped by
        // the alternative's likelihood mass.
        assert_eq!(state.k, 1);
        assert_eq!(state.alternative, Some(1e-40));
        let p = state.p.expect("finite p");
        assert!((p - 1.0 / (1.0 + 1e-40)).abs() < 1e-12);
    }
}
