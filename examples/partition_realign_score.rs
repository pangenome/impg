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
//!
//! 7. THE GENERATED-HERE/GENERATED-ELSEWHERE MARGINALIZATION (slice C,
//!    the owner's rho-squared principle: conserved material cancels;
//!    locally-variable material discriminates). A locus's record set
//!    is defined by the routing's territory touch, so conserved
//!    pockets bring foreign-row-generated and cross-window reads into
//!    the locale's evidence; a candidate that spells such material
//!    would otherwise collect it at full match score while a candidate
//!    that does not pays the all-mismatch floor — even though the
//!    read's own genome generated it somewhere the candidate models
//!    only as "outside this locale". Each read's likelihood under a
//!    candidate therefore MARGINALIZES the two generation branches:
//!
//!        LL(read | candidate) = logsumexp( local-spell branch ,
//!                                          E(read) )
//!
//!    where the local-spell branch is the per-(record, row) likelihood
//!    of rule 2 (the uniform placement prior + the anchored
//!    placements' backbone + skipped-base votes) and E(read) is the
//!    CANDIDATE-INDEPENDENT elsewhere branch — the read explained by
//!    the genome outside this locale at the derived background rate
//!    (rule: the read matches its origin at the MEASURED per-base rate
//!    — the census's zero-unmatched-bp measurement over 1,063,522
//!    verified occurrences — under the max-entropy uniform origin
//!    over the donor's diploid genome, the same max-entropy spirit as
//!    the old Poisson model's beta = M/|U|: the unit of match
//!    likelihood spread uniformly over the universe of placements;
//!    the genome size measured from the panel, the instrument's own
//!    truth-free genome universe). THE FLOOR for unexplained reads is
//!    E(read), not the all-mismatch 150*B (which survives only as the
//!    per-placement lower bound inside the local branch's logsumexp,
//!    a bit-exact no-op wherever an anchor pins). A read whose local
//!    placements score no better than E has ~zero discriminating
//!    power between candidates — it does not vote beyond E; the
//!    logsumexp is MONOTONE in the local branch, so a candidate that
//!    spells a read strictly better than every rival under the local
//!    branch still wins that read after the marginalization (the
//!    L4-control preservation property, unit-proven).
//!
//!    THE POCKET-READ ANATOMY (measured first, at chrI L7 — the
//!    form-deciding gate): the winner's pocket placements are 100%
//!    collinear (one continuous path through the pocket row's spell)
//!    with m mean 148.7/150 — the CONTINUOUS world, so no separate
//!    placement-continuity rule is imposed; the marginalization
//!    carries the whole weight. The anatomy's second finding is
//!    recorded in the receipts and the docs note: the pocket reads'
//!    FULL-READ donor occurrences sit on the sample's own S288C/SK1
//!    chrI paths — 66% INSIDE the window's own truth rows, 27% in the
//!    neighboring seam — i.e. the reads are (mostly) locally
//!    generated, and the truth's failure to place them is the pin
//!    skeleton's stored-walk frame blindness (slice A's binding
//!    lesson: rc-frame-qualified anchors are absent from the stored
//!    path walk), NOT an inability of the truth's rows to spell them.
//!    The marginalization bounds that artifact's cost at E; the frame
//!    repair itself is named as the next lever.

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
use std::sync::atomic::{AtomicU64, Ordering};
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
    /// The pocket-read anatomy sidecar (slice C: per unit, classified
    /// by whether the winner and truth pairs place it, the placement
    /// anatomy on the winner and truth folds plus the read's donor-
    /// path full-read verification). Pure emission; no scoring change.
    #[arg(long)]
    anatomy_out: Option<PathBuf>,
    /// The pin-skeleton sidecar (slice D: per fold, the canonical-
    /// scheme steps with per-step frame tags beside the stored path
    /// walk's steps and the position diff — the frame-repair audit).
    #[arg(long)]
    skeleton_out: Option<PathBuf>,
    /// The frame-audit-only mode (slice D's diagnosis): emit the
    /// per-row skeleton diff (stored walk vs canonical scheme) for
    /// every axis-partition row of the given loci and exit — no
    /// reads, no census, no scoring. The output path is --out.
    #[arg(long, default_value_t = false)]
    frame_audit_only: bool,
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

/// The anchor window start under one placement (slice A, verbatim).
fn anchor_window_start(canonical_pos: u64, start: u64, span: u64, k: u64, orientation: u8) -> u64 {
    if orientation == 0 {
        start + canonical_pos
    } else {
        start + span - k - canonical_pos
    }
}

/// One candidate key's canonical anchor chain with the projection
/// geometry the binding verifies against (slice A).
struct KeyShape {
    canonical: Vec<(i32, u64)>,
    span: u64,
    read: usize,
    mirrored: bool,
    own_positions: Vec<u64>,
}

/// The per-shift crop of a hoisted (union-fetched) occurrence range
/// (phase 1, lever 1): exactly the bytes the pre-rebuild per-shift
/// fetch returned. All three origin shifts share the high bound
/// `start + span + 1` (clamped to the path length by `fetch_seq`
/// before the call) and differ only in the low bound
/// `max(0, start + shift - 1)`, so from the union buffer — fetched
/// once per (occurrence, candidate) starting at panel coordinate
/// `seq_lo` — the per-shift crop is the buffer's tail from the
/// shift's low bound (empty when that low bound reaches past the
/// clamped high end, the empty-crop case of the old fetch).
fn verify_shift_crop<'a>(seq: &'a [u8], seq_lo: u64, occ_start: u64, shift: i64) -> &'a [u8] {
    let lo = (occ_start as i64 + shift - 1).max(0) as u64;
    let from = lo.saturating_sub(seq_lo) as usize;
    seq.get(from..).unwrap_or(&[])
}

/// Sequence verification of one key at one occurrence (slice A,
/// verbatim semantics; phase 1 hoisted the fetch): every anchor
/// window must sit at its projected panel position and equal the
/// representative read's own k-mer. `seq` is the once-per-
/// (occurrence, candidate) union range fetch covering all three
/// origin shifts, starting at panel coordinate `seq_lo`; each shift
/// verifies against the exact slice its own per-shift fetch
/// returned before the rebuild (unit-proven; the identity gate
/// proves the receipts unchanged end-to-end).
fn verify_key_at(
    shape: &KeyShape,
    occ: &CensusOccurrence,
    reads: &[Vec<u8>],
    seq: &[u8],
    seq_lo: u64,
    k: u64,
    shift: i64,
) -> bool {
    let lo = (occ.start as i64 + shift - 1).max(0) as u64;
    let seq = verify_shift_crop(seq, seq_lo, occ.start, shift);
    let origin = occ.start as i64 + shift;
    let read = &reads[shape.read];
    let forward = physical_forward(occ.orientation, shape.mirrored);
    for (j, &(node, canonical_pos)) in shape.canonical.iter().enumerate() {
        let _ = node;
        let window_lo =
            anchor_window_start(canonical_pos, origin as u64, shape.span, k, occ.orientation);
        let rel_lo = window_lo as i64 - lo as i64;
        if rel_lo < 0 || rel_lo + k as i64 > seq.len() as i64 {
            return false;
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
            return false;
        }
    }
    true
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
// The canonical-scheme pin skeleton (slice D's frame repair).
//
// Slice A's binding lesson, applied to the slice-B pin machinery: the
// syng's stored path walk (`walk_path_range`, what the partition GFAs'
// P lines spell) carries only ONE frame's syncmer selection, so anchor
// k-mers that are rc-frame-qualified on a path are ABSENT from that
// walk — a fold pin skeleton built from it is frame-blind and cannot
// pin records whose anchors qualify in the rc frame (the slice-C
// anatomy's named cause: 2,059 in-window full-match units at chrI L7
// with no valid truth placement). The repair reproduces the routing's
// own canonical-scheme territory step extraction (ported from
// build_territory_index_rows): per position, the frame that spells the
// k-mer's canonical form min(K, rc(K)) forward decides; the chosen
// frame's step node is kept (node identity is frame-independent — the
// same interned syncmer id, the sign carries the frame's orientation)
// and a position whose canonical frame did not qualify carries NO
// step. The census records are derived under this same scheme, so
// every anchor of every record is present at every true occurrence of
// either orientation.
// ---------------------------------------------------------------------------

/// The per-position canonical-scheme selection over the two frames'
/// step maps (positions absolute bp, step nodes signed). Frame tag:
/// 0 = the stored forward frame's selection, 1 = the rc frame's
/// qualification.
fn canonical_scheme_steps(
    forward_map: &BTreeMap<u64, i32>,
    reverse_map: &BTreeMap<u64, i32>,
    seq: &[u8],
    seq_lo: u64,
    k: u64,
) -> Vec<(u64, i32, u8)> {
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
            let frame = if canonical_forward { 0u8 } else { 1u8 };
            steps.push((bp, node, frame));
        }
    }
    steps
}

/// The canonical-scheme steps overlapping one row range [start, end)
/// — the pin skeleton's step source. Every emitted step is verified BY
/// SEQUENCE against the AGC-fetched panel sequence (slice A's
/// exactness discipline): the window at the claimed position must
/// equal the claimed node's interned syncmer sequence, orientation-
/// aware per the sign (`syncmer_seq(negative)` is the reverse
/// complement), and the frame decision must agree with the window's
/// canonical form. Returns the steps plus the census of the stored
/// walk's diff (the diagnosis): `added` = positions the stored walk
/// lacks (the rc-frame-only anchors — the blinded class),
/// `replaced` = positions where the canonical scheme keeps a different
/// node than the stored walk, `dropped` = stored positions the
/// canonical scheme excludes (their canonical frame did not qualify —
/// no census record can anchor there).
#[allow(clippy::too_many_arguments)]
fn canonical_steps_overlapping(
    panel: &SyngIndex,
    fetch: &(dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>> + Sync),
    path_name: &str,
    start: u64,
    end: u64,
    k: u64,
    stored: &[(u64, i32)],
) -> io::Result<(Vec<(u64, i32, u8)>, u64, u64, u64)> {
    // The rc-symmetric augmentation (the raw single-frame extraction on
    // the fetched range's reverse complement, mapped back to the
    // path's own coordinates with the frame's orientation on the
    // sign). Fetch with one syncmer length of context on each side so
    // every overlapping step's full k-mer is inside the fetched
    // sequence for the canonical-form check.
    let seq_lo = start.saturating_sub(k);
    let seq_hi = end + k;
    let seq = fetch(path_name, seq_lo, seq_hi)?;
    let rc_seq = revcomp(&seq);
    let reverse: Vec<(u64, i32)> =
        impg::genome_inference::mem_records::raw_matched_syncmers(panel, &rc_seq)?
            .into_iter()
            .map(|(signed, q)| (seq_lo + seq.len() as u64 - k - q, -signed))
            .collect();
    // Both frames' selections restricted to the stored walk's own
    // domain: steps whose windows overlap [start, end).
    let keep = |bp: u64| bp + k > start && bp < end;
    let forward_map: BTreeMap<u64, i32> = stored
        .iter()
        .copied()
        .filter(|&(bp, _)| keep(bp))
        .collect();
    let reverse_map: BTreeMap<u64, i32> =
        reverse.into_iter().filter(|&(bp, _)| keep(bp)).collect();
    let steps = canonical_scheme_steps(&forward_map, &reverse_map, &seq, seq_lo, k);
    let mut verified = 0u64;
    let mut added = 0u64;
    let mut replaced = 0u64;
    for &(bp, node, frame) in &steps {
        let rel = (bp - seq_lo) as usize;
        let window = &seq[rel..rel + k as usize];
        let mut interned = panel.syncmer_seq(node);
        // The syng's kmerHashSeq may return lowercase bases; the GFA
        // writer's own convention uppercases them (write_segment's
        // make_ascii_uppercase) — match it.
        interned.make_ascii_uppercase();
        ensure(
            window == interned.as_slice(),
            &format!(
                "canonical-scheme step fails sequence verification against the fetched \
                 sequence: {path_name} bp {bp} node {node}"
            ),
        )?;
        let canonical_forward =
            impg::genome_inference::mem_records::window_is_canonical_forward(window);
        ensure(
            canonical_forward == (frame == 0),
            "canonical-scheme frame decision disagrees with the window's canonical form",
        )?;
        verified += 1;
        match forward_map.get(&bp) {
            Some(&stored_node) if stored_node == node => {}
            Some(_) => replaced += 1,
            None => added += 1,
        }
    }
    // Completeness of the forward part: every stored step at a
    // canonical-forward position must survive with the stored node.
    for (&bp, &node) in forward_map.iter() {
        let rel = (bp - seq_lo) as usize;
        let canonical_forward = seq
            .get(rel..rel + k as usize)
            .map(|window| {
                impg::genome_inference::mem_records::window_is_canonical_forward(window)
            })
            .unwrap_or(false);
        if canonical_forward {
            ensure(
                steps.iter().any(|&(s_bp, s_node, _)| s_bp == bp && s_node == node),
                "canonical-scheme skeleton dropped a canonical-forward stored step",
            )?;
        }
    }
    let chosen_positions: BTreeSet<u64> = steps.iter().map(|&(bp, _, _)| bp).collect();
    let dropped = forward_map
        .keys()
        .filter(|&&bp| !chosen_positions.contains(&bp))
        .count() as u64;
    Ok((steps, verified + added, added + replaced, dropped))
}

// ---------------------------------------------------------------------------
// Candidate rows (slice A's RowStore, verbatim semantics; slice D: the
// pin skeleton's step source is the canonical scheme by default).
// ---------------------------------------------------------------------------

/// The per-row skeleton diff over the CONTAINED steps (windows fully
/// inside the row extent — the pin-skeleton domain): the canonical-
/// scheme steps with frame tags beside the stored walk's steps, plus
/// the position census (added = the rc-frame-only anchors the stored
/// walk lacks; replaced = positions where the canonical scheme keeps
/// a different node; dropped = stored positions the canonical scheme
/// excludes).
struct SkeletonDiff {
    canonical: Vec<(u64, i32)>,
    frames: Vec<(u64, u8)>,
    stored: Vec<(u64, i32)>,
    added: u64,
    replaced: u64,
    dropped: u64,
    kept_forward: u64,
    kept_reverse: u64,
}

fn contained_skeleton_diff(row: &RowFold, k: u64) -> SkeletonDiff {
    let contained_of = |walk: &[(u64, i32)]| {
        walk.iter()
            .filter(|&&(bp, _)| bp >= row.start && bp + k <= row.end)
            .map(|&(bp, node)| (bp - row.start, node))
            .collect::<Vec<(u64, i32)>>()
    };
    let stored = contained_of(&row.stored_walk);
    let canonical = contained_of(&row.walk);
    let frames = row
        .walk
        .iter()
        .zip(row.frames.iter())
        .filter(|(&(bp, _), _)| bp >= row.start && bp + k <= row.end)
        .map(|(&(bp, _), &frame)| (bp - row.start, frame))
        .collect::<Vec<(u64, u8)>>();
    let stored_positions: BTreeSet<u64> = stored.iter().map(|&(bp, _)| bp).collect();
    let canonical_positions: BTreeSet<u64> = canonical.iter().map(|&(bp, _)| bp).collect();
    let stored_nodes: BTreeMap<u64, i32> = stored.iter().copied().collect();
    let mut added = 0u64;
    let mut replaced = 0u64;
    let mut dropped = 0u64;
    let mut kept_forward = 0u64;
    let mut kept_reverse = 0u64;
    for &(bp, node) in &canonical {
        match stored_nodes.get(&bp) {
            None => added += 1,
            Some(&stored_node) if stored_node == node => {}
            Some(_) => replaced += 1,
        }
    }
    for &bp in &stored_positions {
        if !canonical_positions.contains(&bp) {
            dropped += 1;
        }
    }
    for &(_, frame) in &frames {
        if frame == 0 {
            kept_forward += 1;
        } else {
            kept_reverse += 1;
        }
    }
    SkeletonDiff {
        canonical,
        frames,
        stored,
        added,
        replaced,
        dropped,
        kept_forward,
        kept_reverse,
    }
}

struct RowFold {
    #[allow(dead_code)]
    path: usize,
    start: u64,
    end: u64,
    seq: Vec<u8>,
    /// The pin skeleton's steps overlapping the row range (absolute
    /// bp, signed node), sorted by bp — the CANONICAL-SCHEME selection
    /// by default (slice D's frame repair), or the stored path walk
    /// under the identity gate (the slice-B/C before-record).
    walk: Vec<(u64, i32)>,
    /// The frame tag per walk step (0 = the stored forward frame's
    /// selection, 1 = the rc frame's qualification) — the audit view;
    /// empty under the identity gate.
    frames: Vec<u8>,
    /// The stored path walk (the committed partition GFAs' own P-line
    /// steps) — the audit's comparison baseline.
    stored_walk: Vec<(u64, i32)>,
    /// (partition, member index) pairs of every partition holding this
    /// row.
    members: Vec<(u32, usize)>,
}

/// The axis-partition row store: rows loaded lazily per partition,
/// deduplicated by (path, start, end) across partitions.
struct RowStore<'a> {
    panel: &'a SyngIndex,
    path_of_name: &'a HashMap<String, usize>,
    fetch: &'a (dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>> + Sync),
    k: u64,
    /// The identity gate: keep the stored path walk as the pin
    /// skeleton (env IMPG_REALIGN_STORED_WALK_SKELETON; assessment-side
    /// diagnostic for the before/after pairing — the repair is the
    /// default).
    stored_frame: bool,
    folds: Vec<RowFold>,
    dedup: HashMap<(usize, u64, u64), usize>,
    by_partition: HashMap<u32, Vec<usize>>,
    by_partition_path: HashMap<(u32, usize), Vec<usize>>,
}

impl<'a> RowStore<'a> {
    fn new(
        panel: &'a SyngIndex,
        path_of_name: &'a HashMap<String, usize>,
        fetch: &'a (dyn Fn(&str, u64, u64) -> io::Result<Vec<u8>> + Sync),
        k: u64,
        stored_frame: bool,
    ) -> Self {
        RowStore {
            panel,
            path_of_name,
            fetch,
            k,
            stored_frame,
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
                    let mut stored_walk: Vec<(u64, i32)> = self
                        .panel
                        .walk_path_range(path_idx, member.start, end)?
                        .into_iter()
                        .map(|(node, bp)| (bp, node))
                        .collect();
                    stored_walk.sort_unstable_by_key(|&(bp, _)| bp);
                    let (walk, frames) = if self.stored_frame {
                        (stored_walk.clone(), Vec::new())
                    } else {
                        let (steps, _verified, _changed, _dropped) = canonical_steps_overlapping(
                            self.panel,
                            self.fetch,
                            &member.path_name,
                            member.start,
                            end,
                            self.k,
                            &stored_walk,
                        )?;
                        (
                            steps.iter().map(|&(bp, node, _)| (bp, node)).collect(),
                            steps.iter().map(|&(_, _, frame)| frame).collect(),
                        )
                    };
                    let index = self.folds.len();
                    self.folds.push(RowFold {
                        path: path_idx,
                        start: member.start,
                        end,
                        seq,
                        walk,
                        frames,
                        stored_walk,
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

    /// The model's derived per-placement MINIMUM likelihood: every
    /// base of the read mismatching (READ_LENGTH * ln(eps/3)) — the
    /// rigorous lower bound of the read's likelihood under ANY
    /// placement. It remains the stated lower bound of the local
    /// branch's per-placement terms; as a logsumexp term it
    /// underflows to +0.0 exactly wherever an anchored placement
    /// exists (150 * ln(eps/3) ~ -1546 against placement scores ~ -1).
    #[inline]
    fn floor(&self) -> f64 {
        READ_LENGTH as f64 * self.b
    }
}

/// The plain two-term logsumexp.
#[inline]
fn logsumexp2(a: f64, b: f64) -> f64 {
    let max = a.max(b);
    if max == f64::NEG_INFINITY {
        return max;
    }
    max + ((a - max).exp() + (b - max).exp()).ln()
}

/// THE ELSEWHERE BRANCH E(read) (slice C, derived — no tuning
/// constants): each read was generated at ONE origin, uniform
/// (max-entropy) over the donor's diploid genome; a candidate models
/// the locale only, so the read's likelihood under a candidate
/// marginalizes generated-here (the candidate's local spell) against
/// generated-elsewhere (the genome outside this locale).
///
/// The elsewhere branch: the read matches its origin at the MEASURED
/// per-base rate (the census measured zero unmatched placement bp
/// over 1,063,522 verified occurrences — reads are exact segments of
/// their origins, at the reads' own Phred-derived rate A), and the
/// origin is uniform over the genome's 2*(G-149) both-strand
/// full-read placements — the same max-entropy spirit as the old
/// Poisson model's beta = M/|U|: the unit of match likelihood spread
/// uniformly over the universe of placements, the universe measured
/// from the panel (the instrument's genome universe, truth-free: the
/// mean per-strain total path length, doubled for the diploid).
///
///   E(read) = len(read) * A - ln(2 * (G_diploid - len(read) + 1))
///
/// With uniform L150 reads and uniform quality bytes this is ONE
/// derived constant per run. THE FLOOR for a read a candidate cannot
/// place is E — not the all-mismatch 150*B (which survives only as
/// the per-placement lower bound inside the local branch, a bit-exact
/// no-op wherever an anchor pins).
#[inline]
fn elsewhere_log_prob(read_len: usize, match_log_prob: f64, genome_diploid_bp: f64) -> f64 {
    read_len as f64 * match_log_prob
        - (2.0 * (genome_diploid_bp - read_len as f64 + 1.0)).ln()
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
    /// The representative row's CONTAINED walk steps (windows fully
    /// inside the row extent), relative to its start — the pinnable
    /// anchor skeleton; the committed partition GFAs' P lines spell
    /// these steps among their edge-overlapping and gap steps.
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

/// The COLLINEARITY verdict of one placement's pins (the slice C
/// anatomy): the pinned anchors must sit at ONE offset for the
/// placement to be ONE CONTINUOUS path through the fold's spell (a
/// substitution-only generation spells the read along a single
/// offset; the piecewise serving of slice B can collect flanks at
/// mutually incompatible offsets — the flank-indel class). For
/// physically-forward placements the invariant is r_j - own_hull_j;
/// for physically-reverse placements the reflection r_j + own_hull_j
/// (the reverse geometry). A TRUE placement is always collinear: the
/// read is a contiguous segment of the row's source genome, the true
/// assignment is among the monotone assignments, so every pinned
/// anchor sits at its true position.
fn placement_offset(
    pinned: &[Option<i64>],
    own_hull: &[u64],
    forward: bool,
) -> (bool, Option<i64>) {
    let mut offsets: Vec<i64> = Vec::new();
    for (j, &r) in pinned.iter().enumerate() {
        if let Some(r) = r {
            let own = own_hull[j] as i64;
            offsets.push(if forward { r - own } else { r + own });
        }
    }
    match offsets.first() {
        None => (false, None),
        Some(&first) => (offsets.iter().all(|&o| o == first), Some(first)),
    }
}

/// The ONE-CONTINUOUS-PATH score: every base of the read projected
/// at the placement's single offset — the whole-read form of the
/// backbone + serving-anchor votes, INCLUDING the bases slice B
/// abstained inside unpinned anchor windows (at a collinear offset
/// the whole read projects, pinned or not). Requires collinear pins.
/// The third return is the out-of-extent base count (a measurement
/// for the anatomy; the placement-validity rule itself is
/// unchanged).
#[allow(clippy::too_many_arguments)]
fn single_offset_score(
    unit: &Unit,
    offset_const: i64,
    forward: bool,
    fold_seq: &[u8],
    fold_len: u64,
    reads: &[Vec<u8>],
    k: u64,
) -> (u32, u32, u32) {
    let read = &reads[unit.read];
    let mut matches = 0u32;
    let mut mismatches = 0u32;
    let mut out_of_extent = 0u32;
    for i in 0..read.len() as i64 {
        let d = i - unit.w_lo as i64;
        let coord = if forward {
            offset_const + d
        } else {
            offset_const + k as i64 - 1 - d
        };
        if coord < 0 || coord >= fold_len as i64 {
            out_of_extent += 1;
            continue;
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
    (matches, mismatches, out_of_extent)
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

// (Phase 3 — the parallelism pilot: THE SEAM MAP. The two measured-
// dominant seams (the per-record binding loop and the per-locus loop)
// run on DEDICATED std::threads at a configurable width, NOT on the
// rayon pool. MEASURED REASON, both effects caught by the first 2-wide
// pilot run and its control: (1) the committed serial runs were never
// single-threaded in the loci loop — the per-locus scoring/classes/
// spectrum par_iters ride rayon's default pool (the whole box), and
// bounding RAYON_NUM_THREADS to the rung width silently serialized
// THOSE phases too (chrMT loci 1.8s serial-default -> 8.2s at a
// 2-wide pool); (2) nesting the seam par_iter in the SAME pool adds
// contention on top (18.5s: two concurrent loci each enqueueing
// thousands of inner jobs onto the same 2 workers). The seam map
// instead spawns `width` plain scoped threads, stride-assigned over
// the items; from a non-rayon thread every inner par_iter installs
// rayon's GLOBAL pool — the committed inner behavior, unchanged. The
// seam bodies are pure per-item functions over immutable shared state
// (compiler-enforced: the map takes `&(dyn Fn + Sync)`); each result
// lands in its own pre-sized slot and is read out in item order, so
// the outputs are interleaving-independent by construction — the
// byte-identity gate against the serial receipts is the race
// detector.)
fn seam_map<T: Sync, R: Send>(
    items: &[T],
    width: usize,
    f: &(dyn Fn(&T) -> io::Result<R> + Sync),
) -> io::Result<Vec<R>> {
    ensure(width >= 1, "the seam width must be at least 1")?;
    let n = items.len();
    if n == 0 || width == 1 {
        return items.iter().map(f).collect();
    }
    let out: std::sync::Mutex<Vec<Option<R>>> =
        std::sync::Mutex::new((0..n).map(|_| None).collect());
    let errors: std::sync::Mutex<Vec<(usize, io::Error)>> =
        std::sync::Mutex::new(Vec::new());
    let out_ref = &out;
    let errors_ref = &errors;
    std::thread::scope(|scope| {
        for w in 0..width {
            scope.spawn(move || {
                let mut i = w;
                while i < n {
                    match f(&items[i]) {
                        Ok(value) => out_ref.lock().unwrap()[i] = Some(value),
                        Err(error) => errors_ref.lock().unwrap().push((i, error)),
                    }
                    i += width;
                }
            });
        }
    });
    let errors = errors.into_inner().unwrap();
    if let Some((_, error)) = errors.into_iter().min_by_key(|(i, _)| *i) {
        return Err(error);
    }
    let out = out.into_inner().unwrap();
    Ok(out
        .into_iter()
        .map(|slot| slot.expect("the seam map left an unset slot"))
        .collect())
}

fn main() -> io::Result<()> {
    let started = Instant::now();
    let options = Options::parse();
    let rss = RssGuard::new(options.rss_budget_gib);
    // (Phase 3 — the parallelism pilot: THE SEAM SWITCH. The committed
    // default is the serial code path — byte-for-byte the phase-1
    // behavior. With IMPG_REALIGN_PARALLEL_SEAMS set (and the width in
    // IMPG_REALIGN_SEAM_WIDTH, default 2), the two measured-dominant
    // seams (the per-record binding loop and the per-locus loop) run
    // on dedicated seam threads — see `seam_map` for the measured
    // reasons the seams do NOT share the rayon pool with the per-locus
    // scoring/classes par_iters. The seam bodies are pure per-record /
    // per-locus computations whose results land in per-item slots read
    // out in item order, so the receipts are thread-interleaving-
    // independent BY CONSTRUCTION; the byte-identity gate against the
    // serial receipts is the race detector.)
    let seam_parallel = std::env::var("IMPG_REALIGN_PARALLEL_SEAMS").is_ok();
    let seam_width: usize = std::env::var("IMPG_REALIGN_SEAM_WIDTH")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(2);
    ensure(seam_width >= 1, "the seam width must be at least 1")?;
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
    // (Phase instrumentation for the dominance measurement: timers and
    // counters only, all emitted to stderr; no receipt field changes,
    // no behavior change.)
    let inputs_started = Instant::now();
    let fetch_calls = AtomicU64::new(0);
    let fetch_bytes = AtomicU64::new(0);
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
    // (Phase 1, lever 2 — the call-bound fix: a per-lane in-memory
    // sequence cache over the route sources. Phase 0 measured the
    // binding verify loop at ~90%/93% of the whole instrument wall,
    // CALL-bound at ~129-190us per ~72bp crop — the per-call
    // random-access seek/decompress of Sources::fetch into the AGC —
    // and measured the touched subset: one 93.3MB AGC file (3.34GB
    // uncompressed over 9,901 lanes), of which the exhaustive chrI
    // run touches 840 occurrence lanes (384.5MB) and the partition-
    // map member rows 1,286 lanes (514.7MB total). DERIVED CHOICE:
    // load each touched LANE once, whole, through the same validated
    // Sources::fetch path (one sequential decompression per lane —
    // the per-locality batch-prefetch option taken at its natural
    // locality, the contig; a window/block cache would still pay a
    // decompression call per 64KB window), then serve every crop as
    // a memory slice of the cached buffer. Byte-identity: the lane
    // buffer IS the 0..len crop Sources::fetch returns (length- and
    // alphabet-validated, uppercased); a crop of it is a crop of the
    // answer, so every fetch_seq call returns the same bytes as the
    // per-call AGC fetch it replaces — the identity gate proves it
    // end-to-end. The cache holds only touched lanes (<= ~515MB at
    // chrI), never the 3.34GB whole panel.)
    struct LaneCache {
        sources: routes::Sources,
        lane_lengths: Vec<u64>,
        seqs: Vec<std::sync::Mutex<Option<std::sync::Arc<Vec<u8>>>>>,
        loads: AtomicU64,
        bytes: AtomicU64,
        nanos: AtomicU64,
    }
    impl LaneCache {
        fn lane(&self, source: usize) -> io::Result<std::sync::Arc<Vec<u8>>> {
            let started = Instant::now();
            let mut slot = self.seqs[source].lock().unwrap();
            if let Some(seq) = slot.as_ref() {
                return Ok(std::sync::Arc::clone(seq));
            }
            let len = self.lane_lengths[source];
            let seq = std::sync::Arc::new(self.sources.fetch(source, 0, len)?);
            self.loads.fetch_add(1, Ordering::Relaxed);
            self.bytes.fetch_add(len, Ordering::Relaxed);
            self.nanos
                .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            *slot = Some(std::sync::Arc::clone(&seq));
            Ok(seq)
        }
        fn fetch(&self, source: usize, start: u64, end: u64) -> io::Result<Vec<u8>> {
            // The Sources::fetch contract, mirrored exactly.
            ensure(
                source < self.lane_lengths.len()
                    && start < end
                    && end <= self.lane_lengths[source],
                "invalid source crop",
            )?;
            let lane = self.lane(source)?;
            Ok(lane[start as usize..end as usize].to_vec())
        }
    }
    let lane_cache = LaneCache {
        lane_lengths: sources.lanes.iter().map(|&(_, len)| len).collect(),
        sources,
        seqs: (0..graph.lanes.len())
            .map(|_| std::sync::Mutex::new(None))
            .collect(),
        loads: AtomicU64::new(0),
        bytes: AtomicU64::new(0),
        nanos: AtomicU64::new(0),
    };
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
        fetch_calls.fetch_add(1, Ordering::Relaxed);
        let seq = lane_cache.fetch(source, lo, hi)?;
        fetch_bytes.fetch_add(seq.len() as u64, Ordering::Relaxed);
        Ok(seq)
    };

    // The census receipt (slice A's input, verbatim; not needed in the
    // frame-audit-only mode — the skeleton audit touches no reads).
    let mut census_records: Vec<CensusRecordLine> = Vec::new();
    if !options.frame_audit_only {
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
    }

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
    let n_windows = if options.frame_audit_only {
        axis_partitions.len()
    } else {
        census_records
            .iter()
            .flat_map(|line| line.occurrences.iter())
            .flat_map(|occ| occ.partitions.iter().copied())
            .max()
            .map(|w| w as usize + 1)
            .unwrap_or(0)
    };
    ensure(
        axis_partitions.len() == n_windows,
        "axis partition count does not match the census window span",
    )?;
    let window_partition: Vec<u32> =
        axis_partitions.iter().map(|&(partition, _)| partition).collect();
    for &locus in &loci {
        ensure((locus as usize) < n_windows, "locus outside the window span")?;
    }
    let inputs_seconds = inputs_started.elapsed().as_secs_f64();
    eprintln!(
        "[score] phase inputs: panel+routes+census+partition maps \
         ({} census records, {} windows) [{inputs_seconds:.1}s]",
        census_records.len(),
        n_windows,
    );
    rss.probe("inputs")?;

    // ------------------------------------- slice D gate 1: the frame audit
    // (the frame-blindness diagnosis, pure emission: per axis-partition
    // row of the given loci, the canonical-scheme pin skeleton beside
    // the stored path walk, with the per-position diff — the class
    // census of the frame-blinded anchors. Every canonical step is
    // sequence-verified in-process against the fetched panel sequence.)
    if options.frame_audit_only {
        let audit_started = Instant::now();
        let mut store = RowStore::new(&panel, &path_of_name, &fetch_seq, k, false);
        for &locus in &loci {
            let partition = window_partition[locus as usize];
            let map = maps
                .get(&partition)
                .ok_or_else(|| invalid("axis partition map missing"))?;
            store.ensure_partition(&map.members, partition)?;
        }
        let mut audit = BufWriter::new(File::create(&options.out)?);
        let mut rows_total = 0u64;
        let mut rows_affected = 0u64;
        let mut steps_forward = 0u64;
        let mut steps_reverse = 0u64;
        let mut added_total = 0u64;
        let mut replaced_total = 0u64;
        let mut dropped_total = 0u64;
        for row in &store.folds {
            let diff = contained_skeleton_diff(row, k);
            rows_total += 1;
            steps_forward += diff.kept_forward;
            steps_reverse += diff.kept_reverse;
            added_total += diff.added;
            replaced_total += diff.replaced;
            dropped_total += diff.dropped;
            if diff.added > 0 || diff.replaced > 0 || diff.dropped > 0 {
                rows_affected += 1;
            }
            serde_json::to_writer(
                &mut audit,
                &json!({
                    "partitions": row.members.iter().map(|&(p, _)| p).collect::<Vec<_>>(),
                    "path_name": panel.name_map.path_to_name[row.path],
                    "start": row.start,
                    "end": row.end,
                    "length": row.end - row.start,
                    "steps": diff.canonical.iter().zip(diff.frames.iter())
                        .map(|(&(bp, node), &(_, frame))| json!([bp, node, frame]))
                        .collect::<Vec<_>>(),
                    "stored": diff.stored.iter()
                        .map(|&(bp, node)| json!([bp, node]))
                        .collect::<Vec<_>>(),
                    "added": diff.added,
                    "replaced": diff.replaced,
                    "dropped": diff.dropped,
                    "kept_forward": diff.kept_forward,
                    "kept_reverse": diff.kept_reverse,
                }),
            )?;
            writeln!(audit)?;
        }
        audit.flush()?;
        let partitions: BTreeSet<u32> = loci
            .iter()
            .map(|&l| window_partition[l as usize])
            .collect();
        eprintln!(
            "[score] frame audit: {} partitions, {} rows, {} affected \
             (added {added_total}, replaced {replaced_total}, dropped {dropped_total}; \
             kept forward {steps_forward}, reverse {steps_reverse}) [{:.1}s]",
            partitions.len(),
            rows_total,
            rows_affected,
            audit_started.elapsed().as_secs_f64(),
        );
        rss.probe("frame_audit")?;
        return Ok(());
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
    let quality_seconds = quality_started.elapsed().as_secs_f64();
    eprintln!(
        "[score] scoring rate: uniform Phred {phred} over {total_reads} reads \
         (epsilon {epsilon}, A {}, B {}) [{quality_seconds:.1}s]",
        scoring.a,
        scoring.b,
    );
    rss.probe("quality")?;

    // ------------------------- the ELSEWHERE branch derivation (slice C)
    // The donor's genome size, measured from the panel (the
    // instrument's own genome universe — truth-free): the panel's
    // paths are grouped per (strain, haplotype) — "S288C#0#chrI",
    // "ATV#3#block43_contig1" — so the measured unit is the mean
    // total path length per HAPLOTYPE path set, doubled for the
    // diploid generation model (the balanced diploid sample's
    // documented construction: two haploid homologs).
    let mut strain_totals: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    let elsewhere_started = Instant::now();
    for (name, &len) in panel
        .name_map
        .path_to_name
        .iter()
        .zip(panel.name_map.path_to_length.iter())
    {
        let haplotype: String = name
            .splitn(3, '#')
        .take(2)
            .collect::<Vec<_>>()
            .join("#");
        let entry = strain_totals.entry(haplotype).or_insert((0, 0));
        entry.0 += 1;
        entry.1 += len;
    }
    ensure(!strain_totals.is_empty(), "the panel has no paths")?;
    let panel_mean_haploid_bp: f64 = strain_totals
        .values()
        .map(|&(_, bp)| bp as f64)
        .sum::<f64>()
        / strain_totals.len() as f64;
    let genome_diploid_bp = 2.0 * panel_mean_haploid_bp;
    let elsewhere = elsewhere_log_prob(READ_LENGTH, scoring.a, genome_diploid_bp);
    let strain_list: Vec<serde_json::Value> = strain_totals
        .iter()
        .map(|(strain, &(paths, bp))| json!({"strain": strain, "paths": paths, "bp": bp}))
        .collect();
    eprintln!(
        "[score] elsewhere branch: E = {elsewhere:.6} \
         (panel strains {}, mean haploid {panel_mean_haploid_bp:.1} bp, \
         diploid {genome_diploid_bp:.1} bp)",
        strain_totals.len(),
    );
    eprintln!(
        "[score] phase elsewhere: E derivation over {} haplotype sets \
         [{:.3}s]",
        strain_totals.len(),
        elsewhere_started.elapsed().as_secs_f64(),
    );
    let elsewhere_seconds = elsewhere_started.elapsed().as_secs_f64();

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
    let derive_seconds = derive_started.elapsed().as_secs_f64();
    eprintln!(
        "[score] derive cache loaded: {} reads, {} keys ({derive_seconds:.1}s)",
        reads.len(),
        key_tokens.len(),
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
    // (Phase instrumentation: the binding sub-timers split the phase
    // into its natural cost drivers — the per-record canonical-scheme
    // step extraction, the per-candidate key-shape re-derivation, and
    // the per-occurrence AGC fetch verification.)
    let binding_started = Instant::now();
    let binding_fetch_calls = fetch_calls.load(Ordering::Relaxed);
    let binding_fetch_bytes = fetch_bytes.load(Ordering::Relaxed);
    let binding_lane_loads = lane_cache.loads.load(Ordering::Relaxed);
    let binding_lane_bytes = lane_cache.bytes.load(Ordering::Relaxed);
    struct BindTimers {
        shape_seconds: f64,
        verify_seconds: f64,
        shape_calls: u64,
        verify_calls: u64,
        range_fetches: u64,
    }
    let mut bind_timers = BindTimers {
        shape_seconds: 0.0,
        verify_seconds: 0.0,
        shape_calls: 0,
        verify_calls: 0,
        range_fetches: 0,
    };
    let mut bind_steps_seconds = 0.0f64;
    let mut bind_steps_calls = 0u64;
    let mut bind_candidates = 0u64;
    let mut bind_occurrences = 0u64;
    let mut bind_shaped_keys: BTreeSet<u32> = BTreeSet::new();
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
    // Sequence verification of one key at one occurrence: the hoisted
    // (phase 1, lever 1) union-range verify — see `verify_key_at` at
    // the top level. The fetch is performed once per (occurrence,
    // candidate) below, covering all three origin shifts.
    let mut census_key: Vec<Option<u32>> = vec![None; census_records.len()];
    let mut record_shifts: Vec<Vec<i64>> = vec![Vec::new(); census_records.len()];
    // (Phase 3 — the parallelism pilot, SEAM 1: the per-record binding
    // loop. The body below is a PURE function of the record: it reads
    // only immutable shared state (the census lines, the key indexes,
    // the derive cache, the panel, the reads) plus the thread-safe
    // fetch path (the atomic fetch counters and the mutex-per-lane
    // LaneCache, whose loads are serialized per lane by construction
    // and idempotent), and writes only its own result slot. Results
    // are collected in pilot_records order and merged serially, so
    // the bound keys and shifts are thread-interleaving-independent
    // by construction; the serial path is the committed phase-1 code
    // path with identical semantics, and the byte-identity gate
    // against the serial receipts is the race detector.)
    struct BindOne {
        key: Option<u32>,
        shifts: Vec<i64>,
        occurrences: u64,
        candidates: u64,
        shaped_keys: Vec<u32>,
        timers: BindTimers,
        steps_seconds: f64,
    }
    let bind_one = |record: usize| -> io::Result<BindOne> {
        let line = &census_records[record];
        let occ = line
            .occurrences
            .first()
            .ok_or_else(|| invalid("census record has no occurrences"))?;
        let occurrences = line.occurrences.len() as u64;
        let steps_started = Instant::now();
        let steps = canonical_steps_near(occ.path, occ.start)?
            .into_iter()
            .filter(|&(bp, _)| bp == occ.start)
            .map(|(_, node)| node)
            .collect::<Vec<i32>>();
        let steps_seconds = steps_started.elapsed().as_secs_f64();
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
        let candidates_count = candidates.len() as u64;
        let shifts = [0i64, 1, -1];
        let mut timers = BindTimers {
            shape_seconds: 0.0,
            verify_seconds: 0.0,
            shape_calls: 0,
            verify_calls: 0,
            range_fetches: 0,
        };
        let verify_all = |key: u32, t: &mut BindTimers| -> io::Result<Option<Vec<i64>>> {
            let shape_started = Instant::now();
            let shape = key_shape(key, &key_tokens, &key_reads, &read_records, k)?;
            t.shape_seconds += shape_started.elapsed().as_secs_f64();
            t.shape_calls += 1;
            let mut per_occurrence = Vec::with_capacity(line.occurrences.len());
            for occ in &line.occurrences {
                let mut found: Vec<i64> = Vec::new();
                // (Phase 1, lever 1 — hoist the re-fetch: ONE range
                // fetch per (occurrence, candidate). Phase 0 measured
                // the three-shift loop as a 3.01x pure re-fetch
                // redundancy — the three per-shift ranges differ only
                // in the low bound, all inside
                // [max(0, start-2), start + span + 1). The union fetch
                // replaces the three per-shift fetches; each shift
                // then verifies against the exact slice its own
                // fetch returned before the rebuild.)
                let verify_started = Instant::now();
                let union_lo = (occ.start as i64 - 2).max(0) as u64;
                let union_hi = occ.start + shape.span + 1;
                let seq = fetch_seq(
                    &panel.name_map.path_to_name[occ.path],
                    union_lo,
                    union_hi,
                )?;
                for &shift in &shifts {
                    let verified =
                        verify_key_at(&shape, occ, &reads, &seq, union_lo, k, shift);
                    t.verify_calls += 1;
                    if verified {
                        found.push(shift);
                    }
                }
                t.verify_seconds += verify_started.elapsed().as_secs_f64();
                t.range_fetches += 1;
                match found.len() {
                    0 => return Ok(None),
                    1 => per_occurrence.push(found[0]),
                    _ => return Err(invalid("key verifies at two origin shifts")),
                }
            }
            Ok(Some(per_occurrence))
        };
        let mut surviving: Vec<(u32, Vec<i64>)> = Vec::new();
        // Every CANDIDATE key is shaped before its verify (the serial
        // phase-1 code inserted each candidate into bind_shaped_keys
        // before the verify, surviving or not — the SET is the
        // semantic, collected here and unioned serially in the merge).
        let shaped_keys: Vec<u32> = candidates.iter().copied().collect();
        for key in candidates {
            if let Some(per_occurrence) = verify_all(key, &mut timers)? {
                surviving.push((key, per_occurrence));
            }
        }
        ensure(
            surviving.len() == 1,
            "census record has no unique verifying key",
        )?;
        Ok(BindOne {
            key: Some(surviving[0].0),
            shifts: surviving.swap_remove(0).1,
            occurrences,
            candidates: candidates_count,
            shaped_keys,
            timers,
            steps_seconds,
        })
    };
    let pilot_list: Vec<usize> = pilot_records.iter().copied().collect();
    let bound: Vec<BindOne> = if seam_parallel {
        seam_map(&pilot_list, seam_width, &|&record| bind_one(record))?
    } else {
        pilot_list
            .iter()
            .map(|&record| bind_one(record))
            .collect::<io::Result<Vec<_>>>()?
    };
    for (&record, one) in pilot_records.iter().zip(bound) {
        census_key[record] = one.key;
        record_shifts[record] = one.shifts;
        bind_occurrences += one.occurrences;
        bind_candidates += one.candidates;
        bind_steps_seconds += one.steps_seconds;
        bind_steps_calls += 1;
        bind_shaped_keys.extend(one.shaped_keys);
        bind_timers.shape_seconds += one.timers.shape_seconds;
        bind_timers.verify_seconds += one.timers.verify_seconds;
        bind_timers.shape_calls += one.timers.shape_calls;
        bind_timers.verify_calls += one.timers.verify_calls;
        bind_timers.range_fetches += one.timers.range_fetches;
    }
    eprintln!(
        "[score] pilot records bound: {} ({:.1}s)",
        pilot_records.len(),
        started.elapsed().as_secs_f64(),
    );
    eprintln!(
        "[score] phase binding: {} records / {} occurrences / {} candidate keys \
         ({} distinct shaped) / {} shape derivations / {} verify probes over \
         {} range fetches; steps {:.1}s ({} extractions), shapes {:.1}s, verify \
         {:.1}s, fetches {} ({:.1} MB), lane loads {} ({:.1} MB), phase {:.1}s",
        pilot_records.len(),
        bind_occurrences,
        bind_candidates,
        bind_shaped_keys.len(),
        bind_timers.shape_calls,
        bind_timers.verify_calls,
        bind_timers.range_fetches,
        bind_steps_seconds,
        bind_steps_calls,
        bind_timers.shape_seconds,
        bind_timers.verify_seconds,
        fetch_calls.load(Ordering::Relaxed) - binding_fetch_calls,
        (fetch_bytes.load(Ordering::Relaxed) - binding_fetch_bytes) as f64 / (1024.0 * 1024.0),
        lane_cache.loads.load(Ordering::Relaxed) - binding_lane_loads,
        (lane_cache.bytes.load(Ordering::Relaxed) - binding_lane_bytes) as f64
            / (1024.0 * 1024.0),
        binding_started.elapsed().as_secs_f64(),
    );
    let binding_seconds = binding_started.elapsed().as_secs_f64();
    rss.probe("binding")?;

    // ----------------------- per-record read variants (slice A, verbatim)
    struct ReadVariant {
        mirror: bool,
        read: usize,
        w_lo: u64,
        count: u64,
    }
    let mut record_variants: HashMap<usize, Vec<ReadVariant>> = HashMap::new();
    let variants_started = Instant::now();
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
    eprintln!(
        "[score] phase variants: {} records, per-record read variants \
         [{:.1}s]",
        record_variants.len(),
        variants_started.elapsed().as_secs_f64(),
    );
    let variants_seconds = variants_started.elapsed().as_secs_f64();
    rss.probe("variants")?;

    // ------------------------------------------------------------------ receipts
    let mut report = BufWriter::new(File::create(&options.out)?);
    let mut exactness = BufWriter::new(File::create(&options.exactness_out)?);
    let mut anatomy = match &options.anatomy_out {
        Some(path) => Some(BufWriter::new(File::create(path)?)),
        None => None,
    };
    let ingredients_path = format!("{}.ingredients.jsonl", options.out.display());
    let mut ingredients = BufWriter::new(File::create(&ingredients_path)?);
    let records_path = format!("{}.records.jsonl", options.out.display());
    let mut records_file = BufWriter::new(File::create(&records_path)?);

    // The identity gate (slice D): under IMPG_REALIGN_STORED_WALK_SKELETON
    // the pin skeleton stays the stored path walk (the slice-B/C
    // before-record; assessment-side diagnostic for the before/after
    // pairing). The repair — the canonical-scheme skeleton — is the
    // default.
    let stored_frame = std::env::var("IMPG_REALIGN_STORED_WALK_SKELETON").is_ok();
    let mut skeleton = match &options.skeleton_out {
        Some(path) => Some(BufWriter::new(File::create(path)?)),
        None => None,
    };
    let mut store = RowStore::new(&panel, &path_of_name, &fetch_seq, k, stored_frame);
    let contig = options
        .component
        .splitn(3, '#')
        .nth(2)
        .ok_or_else(|| invalid("component lacks a contig suffix"))?
        .to_string();
    // The axis partitions of every window the pilot records touch
    // (the locality classification needs their rows).
    // (Phase instrumentation: this is the CONTEXT ASSEMBLY phase —
    // per-row sequence fetch, stored-walk extraction, and the
    // canonical-scheme skeleton extraction with in-process sequence
    // verification.)
    let rows_started = Instant::now();
    let rows_fetch_calls = fetch_calls.load(Ordering::Relaxed);
    let rows_fetch_bytes = fetch_bytes.load(Ordering::Relaxed);
    let rows_lane_loads = lane_cache.loads.load(Ordering::Relaxed);
    let rows_lane_bytes = lane_cache.bytes.load(Ordering::Relaxed);
    for &partition in &window_partition {
        let map = maps
            .get(&partition)
            .ok_or_else(|| invalid("axis partition map missing"))?;
        store.ensure_partition(&map.members, partition)?;
    }
    eprintln!(
        "[score] phase rows: {} distinct rows over {} partitions, context \
         assembly (canonical-scheme skeletons; fetches {} / {:.1} MB, \
         lane loads {} / {:.1} MB) [{:.1}s]",
        store.folds.len(),
        window_partition.len(),
        fetch_calls.load(Ordering::Relaxed) - rows_fetch_calls,
        (fetch_bytes.load(Ordering::Relaxed) - rows_fetch_bytes) as f64 / (1024.0 * 1024.0),
        lane_cache.loads.load(Ordering::Relaxed) - rows_lane_loads,
        (lane_cache.bytes.load(Ordering::Relaxed) - rows_lane_bytes) as f64
            / (1024.0 * 1024.0),
        rows_started.elapsed().as_secs_f64(),
    );
    let rows_seconds = rows_started.elapsed().as_secs_f64();
    rss.probe("rows")?;

    let mut total_factorized_checked = 0u64;
    let mut total_factorized_equal = 0u64;
    let mut loci_seconds_total = 0.0f64;
    let skeleton_enabled = skeleton.is_some();
    let anatomy_enabled = anatomy.is_some();

    // (Phase 3 — the parallelism pilot, SEAM 2: the per-locus loop.
    // The body below is a PURE function of the locus: it reads only
    // immutable shared state (the partition maps, the assembled row
    // store, the census lines, the bound keys/shifts/variants, the
    // reads, the derive cache, the panel) plus the thread-safe fetch
    // path (atomic counters; the mutex-per-lane LaneCache) and the
    // read-only RSS guard — and it writes only its own per-locus
    // buffers. Every sidecar line is buffered per locus and written
    // in LOCUS ORDER after the loop, so parallel execution cannot
    // interleave lines and the receipts are thread-interleaving-
    // independent by construction. The serial path is the committed
    // phase-1 code path; the byte-identity gate against the serial
    // receipts is the race detector.)
    struct LocusOutput {
        skeleton: Vec<u8>,
        exactness: Vec<u8>,
        records: Vec<u8>,
        ingredients: Vec<u8>,
        report: Vec<u8>,
        anatomy: Vec<u8>,
        factorized_checked: u64,
        factorized_equal: u64,
        locus_seconds: f64,
    }
    let compute_locus = |locus: u32| -> io::Result<LocusOutput> {
        let locus_started = Instant::now();
        let mut skeleton_buf: Vec<u8> = Vec::new();
        let mut exactness_buf: Vec<u8> = Vec::new();
        let mut records_buf: Vec<u8> = Vec::new();
        let mut ingredients_buf: Vec<u8> = Vec::new();
        let mut report_buf: Vec<u8> = Vec::new();
        let mut anatomy_buf: Vec<u8> = Vec::new();
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
            // The CONTAINED steps only (windows fully inside the row
            // extent): an anchor whose window overlaps the row edge
            // cannot be a pinned backbone anchor (the read must lie
            // fully inside), and keeping it would place its window
            // start at a NEGATIVE relative position — the u64
            // subtraction would wrap and poison the monotone
            // enumeration. Edge-overlapping anchors abstain (the named
            // remainder class).
            let relative: Vec<(u64, i32)> = row
                .walk
                .iter()
                .filter(|&&(bp, _)| bp >= row.start && bp + k <= row.end)
                .map(|&(bp, node)| (bp - row.start, node))
                .collect();
            let key = (row.seq.clone(), relative.clone());
            let index = *fold_map.entry(key).or_insert_with(|| {
                let mut node_positions: HashMap<u32, Vec<(u64, i32)>> = HashMap::new();
                for &(bp, node) in &row.walk {
                    if bp >= row.start && bp + k <= row.end {
                        node_positions
                            .entry(node.unsigned_abs())
                            .or_default()
                            .push((bp - row.start, node.signum()));
                    }
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
                    walk: relative,
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
        // The pin-skeleton sidecar (slice D's frame-repair audit): per
        // fold, the canonical-scheme contained steps with frame tags
        // beside the stored walk's, with the position diff — and the
        // per-locus totals for the receipt. Under the identity gate
        // (stored_frame) the canonical rows were not extracted for
        // scoring, so the sidecar is not emitted there.
        let mut skeleton_rows = 0u64;
        let mut skeleton_affected = 0u64;
        let mut skeleton_added = 0u64;
        let mut skeleton_replaced = 0u64;
        let mut skeleton_dropped = 0u64;
        let mut skeleton_forward = 0u64;
        let mut skeleton_reverse = 0u64;
        if skeleton_enabled {
            for (fold_index, fold) in folds.iter().enumerate() {
                let first = &fold.members[0];
                let path_idx = path_of_name[&first.path_name];
                let clipped_end =
                    first.end.min(panel.name_map.path_to_length[path_idx]);
                let row_index = *store
                    .dedup
                    .get(&(path_idx, first.start, clipped_end))
                    .ok_or_else(|| invalid("fold row missing from the store"))?;
                let row = &store.folds[row_index];
                let diff = contained_skeleton_diff(row, k);
                skeleton_rows += 1;
                skeleton_added += diff.added;
                skeleton_replaced += diff.replaced;
                skeleton_dropped += diff.dropped;
                skeleton_forward += diff.kept_forward;
                skeleton_reverse += diff.kept_reverse;
                if diff.added > 0 || diff.replaced > 0 || diff.dropped > 0 {
                    skeleton_affected += 1;
                }
                serde_json::to_writer(
                    &mut skeleton_buf,
                    &json!({
                        "locus": locus,
                        "partition": partition,
                        "fold": fold_index,
                        "members": fold.members.iter().map(|m| json!({
                            "path_name": m.path_name, "start": m.start, "end": m.end,
                        })).collect::<Vec<_>>(),
                        "length": fold.len,
                        "steps": diff.canonical.iter().zip(diff.frames.iter())
                            .map(|(&(bp, node), &(_, frame))| json!([bp, node, frame]))
                            .collect::<Vec<_>>(),
                        "stored": diff.stored.iter()
                            .map(|&(bp, node)| json!([bp, node]))
                            .collect::<Vec<_>>(),
                        "added": diff.added,
                        "replaced": diff.replaced,
                        "dropped": diff.dropped,
                        "kept_forward": diff.kept_forward,
                        "kept_reverse": diff.kept_reverse,
                    }),
                )?;
                skeleton_buf.push(b'\n');
            }
        }
        eprintln!(
            "[score] locus {locus}: partition {partition}, {} members -> {n_folds} folds [{:.1}s]",
            map.members.len(),
            fold_started.elapsed().as_secs_f64(),
        );

        // ---------------------------------- the locus's records and units
        let units_started = Instant::now();
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
            let shifts = &record_shifts[record];
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
            "[score] locus {locus}: {} records, {n_units} units \
             (touching/locality classification) [{:.1}s]",
            locus_records.len(),
            units_started.elapsed().as_secs_f64(),
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
                        // No anchored placement: the ELSEWHERE branch
                        // (the generated-elsewhere marginalization's
                        // floor; a read a candidate cannot place is
                        // explained by the genome outside the locale at
                        // the derived background rate, not by the
                        // all-mismatch minimum).
                        no_pin_records += 1;
                        ll.push(elsewhere);
                        pinned.push(false);
                        mismatches.push(0);
                        continue;
                    }
                    let mut scores: Vec<f64> = Vec::new();
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
                    // THE MARGINALIZATION: LL(read | fold) =
                    // logsumexp( the local-spell branch (the uniform
                    // placement prior + the anchored placements),
                    // E(read) — the candidate-independent
                    // generated-elsewhere branch ). A local branch
                    // that scores no better than E is absorbed by it
                    // (the read does not vote beyond E); the
                    // logsumexp2 is monotone in the local branch, so
                    // per-unit orderings the local branch preserves
                    // (a full-match truth placement vs a rival's
                    // worse placement) survive the marginalization.
                    ll.push(logsumexp2(
                        unit_log_likelihood(&scores, fold.len),
                        elsewhere,
                    ));
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
            let mut local_scores: Vec<f64> = Vec::new();
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
                local_scores.push(scoring.placement(m, c));
            }
            let read = &reads[unit.read];
            serde_json::to_writer(
                &mut exactness_buf,
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
                    "ll_local": unit_log_likelihood(&local_scores, fold.len),
                    "elsewhere": elsewhere,
                    "placements": placements,
                }),
            )?;
            exactness_buf.push(b'\n');
        }
        let exactness_seconds = exactness_started.elapsed().as_secs_f64();
        eprintln!(
            "[score] locus {locus}: exactness sample {} pairs of {} valid \
             ({} with mismatch votes) [{exactness_seconds:.1}s]",
            sampled.len(),
            valid_pairs.len(),
            mismatch_pairs.len(),
        );

        // -------------------------------------- the pocket-read anatomy
        // (slice C gate 1: THE ANATOMY THAT DECIDES THE FORM of the
        // generated-here/generated-elsewhere marginalization. For
        // every unit, classified by whether the winner and truth
        // pairs place it, the placement anatomy on the winner and
        // truth folds: the anchored/skipped/abstained decomposition,
        // the collinearity verdict, the slice-B piecewise (m, c)
        // BESIDE the one-continuous-path whole-read score, and the
        // flank votes with their projected coordinates — plus the
        // read's full verified-occurrence census and, on the
        // sample's own donor paths (S288C/SK1, the balanced 15x/15x
        // construction), the FULL-READ match at each donor
        // occurrence: the generated-elsewhere evidence. Pure
        // emission; the scoring is unchanged.)
        let anatomy_seconds = if anatomy_enabled {
            let anatomy_started = Instant::now();
            let (winner_first, winner_second) = winner_pair;
            let mut anatomy_folds: Vec<usize> = vec![winner_first, winner_second];
            if let Some((a, b)) = truth_pair {
                for fold in [a, b] {
                    if !anatomy_folds.contains(&fold) {
                        anatomy_folds.push(fold);
                    }
                }
            }
            let truth_places_unit = |u: usize| {
                truth_pair
                    .map(|(a, b)| pinned[a][u] || pinned[b][u])
                    .unwrap_or(false)
            };
            let winner_places_unit = |u: usize| pinned[winner_first][u] || pinned[winner_second][u];
            let is_donor_path = |name: &str| {
                name.starts_with("S288C#0#") || name.starts_with("SK1#0#")
            };
            let mut anatomy_units = 0u64;
            for (u, unit) in units.iter().enumerate() {
                let record = unit.record;
                let canonical =
                    decode_record_tokens(&key_tokens[census_key[record].unwrap() as usize])?;
                let line = &census_records[record];
                let shifts = &record_shifts[record];
                let span = walk_span(&canonical, k);
                let class = match (winner_places_unit(u), truth_places_unit(u)) {
                    (true, false) => "winner_only",
                    (true, true) => "both",
                    (false, true) => "truth_only",
                    (false, false) => "neither",
                };
                let pinning_folds: Vec<usize> = (0..n_folds)
                    .filter(|&f| pinned[f][u])
                    .collect();
                // The donor-path full-read verification: at each
                // verified occurrence of this record on a donor
                // path, the whole read (all 150 bases, strand-aware,
                // origin-shift-corrected) against the donor path's
                // own sequence.
                let mut donor_checks: Vec<serde_json::Value> = Vec::new();
                for (position, occ) in line.occurrences.iter().enumerate() {
                    let name = &panel.name_map.path_to_name[occ.path];
                    if !is_donor_path(name) {
                        continue;
                    }
                    let origin = (occ.start as i64 + shifts[position]).max(0) as u64;
                    let path_len = panel.name_map.path_to_length[occ.path];
                    let forward = physical_forward(occ.orientation, unit.mirror);
                    let lo = origin.saturating_sub(300);
                    let hi = (origin + span + 300).min(path_len);
                    let seq = fetch_seq(name, lo, hi)?;
                    let read = &reads[unit.read];
                    let mut full_m = 0u32;
                    let mut full_c = 0u32;
                    let mut full_oob = 0u32;
                    for i in 0..read.len() as i64 {
                        let d = i - unit.w_lo as i64;
                        let coord = if forward {
                            origin as i64 + d
                        } else {
                            origin as i64 + span as i64 - 1 - d
                        };
                        if coord < 0
                            || coord >= path_len as i64
                            || (coord as u64) < lo
                            || coord as u64 >= hi
                        {
                            full_oob += 1;
                            continue;
                        }
                        let base = seq[(coord - lo as i64) as usize];
                        let read_base = if forward {
                            read[i as usize]
                        } else {
                            complement(read[i as usize])
                        };
                        if read_base == base {
                            full_m += 1;
                        } else {
                            full_c += 1;
                        }
                    }
                    donor_checks.push(json!({
                        "path": name, "origin": origin,
                        "orientation": occ.orientation,
                        "m": full_m, "c": full_c, "oob": full_oob,
                    }));
                }
                // The placement anatomy on the winner and truth
                // folds (the pocket rows and the rows the truth
                // spells).
                let emit_votes = class == "winner_only" || u % 16 == 0;
                let mut placements: Vec<serde_json::Value> = Vec::new();
                for &fold_index in &anatomy_folds {
                    let fold = &folds[fold_index];
                    let positions_of = |node: u32| fold.positions_of(node);
                    for &orientation in &touching[&record].orientations {
                        let (pinned_vec, _ambiguous) =
                            anchor_correspondence(&canonical, orientation, &positions_of);
                        // The per-anchor candidate census (the diagnosis of
                        // all-None correspondences: which anchors the
                        // fold's walk carries, with which signs).
                        let anchor_candidates: Vec<serde_json::Value> = canonical
                            .iter()
                            .map(|&(node, _)| {
                                let required = required_sign(node, orientation);
                                let positions = positions_of(required.unsigned_abs());
                                let with_sign = positions
                                    .iter()
                                    .filter(|&&(_, sign)| sign == required.signum())
                                    .count();
                                json!([
                                    required,
                                    with_sign,
                                    positions.len(),
                                ])
                            })
                            .collect();
                        if pinned_vec.iter().all(|p| p.is_none()) {
                            placements.push(json!({
                                "fold": fold_index,
                                "orientation": orientation,
                                "valid": false,
                                "windows_inside": false,
                                "pinned": pinned_vec,
                                "anchor_candidates": anchor_candidates,
                                "backbone_bp": 0,
                                "abstained_bp": 0,
                                "skipped_bp": 0,
                                "collinear": false,
                                "offset": null,
                                "m": 0, "c": 0,
                                "single": null,
                                "votes": [],
                            }));
                            continue;
                        }
                        let windows_inside =
                            pinned_vec.iter().flatten().all(|&r| r + k <= fold.len);
                        let own_flat = own_positions_of(&canonical, false, 0);
                        let backbone_windows: Vec<(u64, u64)> = pinned_vec
                            .iter()
                            .enumerate()
                            .filter_map(|(j, &r)| r.map(|r2| (own_flat[j], r2)))
                            .map(|(p, _)| (p, p + k))
                            .collect();
                        let backbone_bp: u64 =
                            merge_windows(backbone_windows).iter().map(|&(lo, hi)| hi - lo).sum();
                        let hull_windows: Vec<(u64, u64)> =
                            own_flat.iter().map(|&p| (p, p + k)).collect();
                        let hull_covered: u64 = merge_windows(hull_windows)
                            .iter()
                            .map(|&(lo, hi)| hi - lo)
                            .sum();
                        let pinned_i: Vec<Option<i64>> =
                            pinned_vec.iter().map(|p| p.map(|v| v as i64)).collect();
                        let pin = Pin {
                            orientation,
                            pinned: pinned_i.clone(),
                            backbone_bp,
                        };
                        let forward = physical_forward(orientation, unit.mirror);
                        let own_hull = own_positions_of(&canonical, unit.mirror, 0);
                        let (collinear, offset_const) =
                            placement_offset(&pin.pinned, &own_hull, forward);
                        let own_positions =
                            own_positions_of(&canonical, unit.mirror, unit.w_lo);
                        let skipped = skipped_positions(&own_positions, k, READ_LENGTH as u64);
                        let skipped_bp: u64 =
                            skipped.iter().map(|&(lo, hi)| hi - lo).sum();
                        let piecewise = placement_score(
                            unit, &canonical, &pin, &fold.seq, fold.len, &reads, k,
                        );
                        let (m, c, valid) = match piecewise {
                            Some((m, c)) => (m, c, true),
                            None => (0, 0, false),
                        };
                        let single = if collinear {
                            Some(single_offset_score(
                                unit,
                                offset_const.unwrap(),
                                forward,
                                &fold.seq,
                                fold.len,
                                &reads,
                                k,
                            ))
                        } else {
                            None
                        };
                        let mut votes: Vec<serde_json::Value> = Vec::new();
                        if emit_votes {
                            if let Some((m, c)) = piecewise {
                                let _ = (m, c);
                                for &(lo, hi) in &skipped {
                                    for i in lo..hi {
                                        let d = i as i64 - unit.w_lo as i64;
                                        let (index, own_pos) = match serving_anchor(
                                            &pin.pinned, &own_hull, d,
                                        ) {
                                            Some(found) => found,
                                            None => continue,
                                        };
                                        let r = match pin.pinned[index] {
                                            Some(r) => r,
                                            None => continue,
                                        };
                                        let coord = if forward {
                                            r + (d - own_pos)
                                        } else {
                                            r + (own_pos + k as i64 - 1 - d)
                                        };
                                        if coord < 0 || coord >= fold.len as i64 {
                                            votes.push(json!([i, -1, 0]));
                                            continue;
                                        }
                                        let base = fold.seq[coord as usize];
                                        let read_base = if forward {
                                            reads[unit.read][i as usize]
                                        } else {
                                            complement(reads[unit.read][i as usize])
                                        };
                                        votes.push(json!([
                                            i,
                                            coord,
                                            if read_base == base { 1 } else { 0 },
                                        ]));
                                    }
                                }
                            }
                        }
                        placements.push(json!({
                            "fold": fold_index,
                            "orientation": orientation,
                            "valid": valid,
                            "windows_inside": windows_inside,
                            "pinned": pin.pinned,
                            "anchor_candidates": anchor_candidates,
                            "backbone_bp": backbone_bp,
                            "abstained_bp": hull_covered - backbone_bp,
                            "skipped_bp": skipped_bp,
                            "collinear": collinear,
                            "offset": offset_const,
                            "m": m, "c": c,
                            "single": single,
                            "votes": votes,
                        }));
                    }
                }
                let read = &reads[unit.read];
                serde_json::to_writer(
                    &mut anatomy_buf,
                    &json!({
                        "type": "unit",
                        "locus": locus,
                        "unit": u,
                        "record": record,
                        "variant": unit.variant,
                        "mirror": if unit.mirror { 1 } else { 0 },
                        "w_lo": unit.w_lo,
                        "count": unit.count,
                        "read_fnv": fnv1a64(read),
                        "read": String::from_utf8(read.clone()).unwrap(),
                        "class": class,
                        "pinning_folds": pinning_folds,
                        "donor_checks": donor_checks,
                        "placements": placements,
                    }),
                )?;
                anatomy_buf.push(b'\n');
                anatomy_units += 1;
            }
            eprintln!(
                "[score] locus {locus}: anatomy {anatomy_units} units \
                 [{:.1}s]",
                anatomy_started.elapsed().as_secs_f64(),
            );
            anatomy_started.elapsed().as_secs_f64()
        } else {
            0.0
        };

        // ------------------------------------------------- the named classes
        // (Phase instrumentation: the receipts/IO phase — the named
        // classes, the ingredients (with the full ll_matrix), and the
        // main receipt write.)
        let receipts_started = Instant::now();
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
        serde_json::to_writer(&mut records_buf, &json!({
            "locus": locus,
            "classes": classes_evidence,
        }))?;
        records_buf.push(b'\n');

        // ----------------------------------------------------- the ingredients
        // (Phase 1, lever 3 — the streaming serializer: the ll_matrix is
        // serialized DIRECTLY from the typed matrix, not materialized
        // as a serde_json value tree first. Phase 0 measured the L16
        // tree at +1.15GB live RSS over the scoring baseline (34.5M
        // entries); the streamed bytes are identical — the json! macro
        // emits object keys in serde_json's default (non-preserved)
        // map order, i.e. alphabetical, and every field here
        // serializes through the same writer in that same order. The
        // identity gate proves the sidecar byte-identical.)
        let folds_json: Vec<serde_json::Value> = folds
            .iter()
            .map(|fold| {
                json!({
                    "members": fold.members.iter().map(|m| json!({
                        "path_name": m.path_name, "start": m.start, "end": m.end,
                    })).collect::<Vec<_>>(),
                    "length": fold.len,
                    "sequence": String::from_utf8(fold.seq.clone()).unwrap(),
                    "walk": fold.walk.iter().map(|&(bp, node)| json!([bp, node]))
                        .collect::<Vec<_>>(),
                    "nodes": fold.nodes,
                    "edges": fold.edges,
                })
            })
            .collect();
        let units_json: Vec<serde_json::Value> = units
            .iter()
            .map(|unit| {
                json!({
                    "record": unit.record,
                    "variant": unit.variant,
                    "mirror": if unit.mirror { 1 } else { 0 },
                    "w_lo": unit.w_lo,
                    "count": unit.count,
                    "read": String::from_utf8(reads[unit.read].clone()).unwrap(),
                })
            })
            .collect();
        ingredients_buf.write_all(b"{\"folds\":")?;
        serde_json::to_writer(&mut ingredients_buf, &folds_json)?;
        // The per-unit per-fold log-likelihood matrix (folds x units,
        // row-major; entries at E = no valid anchored placement).
        ingredients_buf.write_all(b",\"ll_matrix\":")?;
        serde_json::to_writer(&mut ingredients_buf, &matrix)?;
        ingredients_buf.write_all(b",\"locus\":")?;
        serde_json::to_writer(&mut ingredients_buf, &locus)?;
        ingredients_buf.write_all(b",\"partition\":")?;
        serde_json::to_writer(&mut ingredients_buf, &partition)?;
        ingredients_buf.write_all(b",\"units\":")?;
        serde_json::to_writer(&mut ingredients_buf, &units_json)?;
        ingredients_buf.write_all(b"}")?;
        ingredients_buf.push(b'\n');
        rss.probe(&format!("locus_{locus}_receipts"))?;

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
        serde_json::to_writer(&mut report_buf, &json!({
            "locus": locus,
            "partition": partition,
            "model": if stored_frame {
                "anchor-realign-v2-marginal"
            } else {
                "anchor-realign-v2-marginal-frame"
            },
            "skeleton": {
                "scheme": if stored_frame { "stored_walk" } else { "canonical_scheme" },
                "rows": skeleton_rows,
                "rows_affected": skeleton_affected,
                "added": skeleton_added,
                "replaced": skeleton_replaced,
                "dropped": skeleton_dropped,
                "kept_forward": skeleton_forward,
                "kept_reverse": skeleton_reverse,
                "note": "the canonical-scheme pin skeleton (slice D's frame \
                         repair: per position the frame that spells the \
                         k-mer's canonical form forward decides; every step \
                         sequence-verified against the fetched panel sequence); \
                         added = rc-frame-only anchors the stored walk lacks",
            },
            "scoring": {
                "phred": phred,
                "epsilon": epsilon,
                "match_log_prob": scoring.a,
                "mismatch_log_prob": scoring.b,
                "read_length": READ_LENGTH,
                "uniform_qualities": true,
                "elsewhere_log_prob": elsewhere,
                "genome_diploid_bp": genome_diploid_bp,
                "panel_mean_haploid_bp": panel_mean_haploid_bp,
                "panel_strains": strain_totals.len(),
                "panel_strain_list": strain_list,
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
                "note": "unplaced reads score at the derived ELSEWHERE \
                         branch E (the generated-elsewhere marginalization; \
                         the genome outside this locale at the measured \
                         per-base rate under the max-entropy uniform origin), \
                         not the all-mismatch 150*B; ranks are bounded by E, \
                         gap magnitudes are E-dominated where unplaced \
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
                "folds_seconds": fold_started.elapsed().as_secs_f64() - scoring_seconds - class_seconds - qual_seconds - exactness_seconds - anatomy_seconds,
                "scoring_seconds": scoring_seconds,
                "classes_seconds": class_seconds,
                "qual_seconds": qual_seconds,
                "exactness_seconds": exactness_seconds,
                "anatomy_seconds": anatomy_seconds,
                "locus_seconds": locus_seconds,
            },
            "rss_kb": rss_kb,
        }))?;
        report_buf.push(b'\n');
        eprintln!(
            "[score] locus {locus}: receipts/IO (named classes + ingredients + \
             the main receipt) [{:.1}s]",
            receipts_started.elapsed().as_secs_f64(),
        );
        eprintln!(
            "[score] locus {locus} done: {n_folds} folds, {n_units} units, {pair_count} classes, \
             truth_rank {:?}, log_gap {log_gap:?}, qual {qual:?} [{locus_seconds:.1}s]",
            truth_rank,
        );
        Ok(LocusOutput {
            skeleton: skeleton_buf,
            exactness: exactness_buf,
            records: records_buf,
            ingredients: ingredients_buf,
            report: report_buf,
            anatomy: anatomy_buf,
            factorized_checked,
            factorized_equal,
            locus_seconds,
        })
    };
    let outputs: Vec<LocusOutput> = if seam_parallel {
        seam_map(&loci, seam_width, &|&locus| compute_locus(locus))?
    } else {
        loci
            .iter()
            .map(|&locus| compute_locus(locus))
            .collect::<io::Result<Vec<_>>>()?
    };
    // The ordered emission: the per-locus buffers are written in LOCUS
    // ORDER — the receipt line order of the committed serial runs.
    for output in outputs {
        if let Some(writer) = skeleton.as_mut() {
            writer.write_all(&output.skeleton)?;
        }
        exactness.write_all(&output.exactness)?;
        records_file.write_all(&output.records)?;
        ingredients.write_all(&output.ingredients)?;
        report.write_all(&output.report)?;
        if let Some(writer) = anatomy.as_mut() {
            writer.write_all(&output.anatomy)?;
        }
        total_factorized_checked += output.factorized_checked;
        total_factorized_equal += output.factorized_equal;
        loci_seconds_total += output.locus_seconds;
    }

    report.flush()?;
    exactness.flush()?;
    ingredients.flush()?;
    records_file.flush()?;
    if let Some(writer) = anatomy.as_mut() {
        writer.flush()?;
    }
    ensure(
        total_factorized_checked == total_factorized_equal,
        "the factorization gate failed somewhere",
    )?;
    eprintln!(
        "[score] complete: loci {loci:?}, factorized placements {total_factorized_equal}/\
         {total_factorized_checked} exact, total {:.1}s",
        started.elapsed().as_secs_f64(),
    );
    eprintln!(
        "[score] lane cache: {} lanes loaded ({:.1} MB, {:.1}s of load time)",
        lane_cache.loads.load(Ordering::Relaxed),
        lane_cache.bytes.load(Ordering::Relaxed) as f64 / (1024.0 * 1024.0),
        lane_cache.nanos.load(Ordering::Relaxed) as f64 / 1e9,
    );
    eprintln!(
        "[score] phase summary: inputs {inputs_seconds:.1}s, quality (FASTQ scan) \
         {quality_seconds:.1}s, elsewhere {elsewhere_seconds:.3}s, derive cache \
         {derive_seconds:.1}s, binding {binding_seconds:.1}s, variants \
         {variants_seconds:.1}s, rows/context assembly {rows_seconds:.1}s, \
         loci (folds+units+scoring+classes+qual+exactness+receipts) \
         {loci_seconds_total:.1}s; AGC fetches {} total ({:.1} MB)",
        fetch_calls.load(Ordering::Relaxed),
        fetch_bytes.load(Ordering::Relaxed) as f64 / (1024.0 * 1024.0),
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
    fn collinearity_holds_for_true_placements_both_strands() {
        // Forward: r_j - own_hull_j invariant.
        let pinned = vec![Some(4i64), Some(14)];
        let own_forward = vec![0u64, 10];
        let (collinear, offset) = placement_offset(&pinned, &own_forward, true);
        assert!(collinear);
        assert_eq!(offset, Some(4));
        // Physically reverse: the reflection r_j + own_hull_j invariant
        // (the mirrored read's own hull positions are reflected).
        let own_mirror = vec![10u64, 0];
        let (collinear, offset) = placement_offset(&pinned, &own_mirror, false);
        assert!(collinear);
        assert_eq!(offset, Some(14));
        // An indel between read and row pins the anchors at mutually
        // incompatible offsets: NOT one continuous path.
        let pinned_indel = vec![Some(4i64), Some(20)];
        let (collinear, _) = placement_offset(&pinned_indel, &own_forward, true);
        assert!(!collinear);
    }

    #[test]
    fn the_one_continuous_path_score_equals_the_placement_score_when_collinear() {
        // The slice-B piecewise (m, c) and the whole-read single-offset
        // score agree exactly on a collinear placement where every
        // anchor pins: the piecewise serving projects the flanks at
        // the same offset, and the previously-abstaining within-window
        // bases agree because the read is a contiguous segment.
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
        let pin = pin();
        let (m, c) = placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K)
            .expect("valid placement");
        assert_eq!((m, c), (29, 1));
        let (sm, sc, oob) = single_offset_score(&unit, 4, true, &seq, 60, &reads, K);
        assert_eq!((sm, sc, oob), (29, 1, 0));
        // The mirrored read at the same pin: the reflected offset
        // const 14 gives the identical per-base evidence.
        let (forward_read_bytes, _) = forward_read(&seq);
        let mirrored = revcomp(&forward_read_bytes);
        let reads = vec![forward_read_bytes, mirrored];
        let unit = Unit {
            record: 0,
            variant: 0,
            mirror: true,
            read: 1,
            w_lo: 10,
            count: 1,
        };
        let (m, c) = placement_score(&unit, &canonical(), &pin, &seq, 60, &reads, K)
            .expect("valid placement");
        assert_eq!((m, c), (29, 1));
        let (sm, sc, oob) = single_offset_score(&unit, 14, false, &seq, 60, &reads, K);
        assert_eq!((sm, sc, oob), (29, 1, 0));
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
    fn the_elsewhere_branch_is_derived_and_bounds_the_floor() {
        // E(read) = len*A - ln(2*(G - len + 1)): the measured per-base
        // match rate under the max-entropy uniform origin over the
        // genome's both-strand full-read placements.
        let a = (1.0 - 1e-4f64).ln();
        let e = elsewhere_log_prob(150, a, 24_400_000.0);
        let expect = 150.0 * a - (2.0f64 * (24_400_000.0 - 150.0 + 1.0)).ln();
        assert!((e - expect).abs() < 1e-9);
        // The floor for an unplaceable read is E — far above the
        // all-mismatch minimum, and the mixture with a missing local
        // branch (NEG_INFINITY) is exactly E.
        assert_eq!(logsumexp2(f64::NEG_INFINITY, e), e);
        // A local branch scoring below E is absorbed by it: the read
        // does not vote beyond E (a full-mismatch or heavily-mismatched
        // placement carries no credit).
        let bad_local = -(2.0f64 * (10_000.0 - 149.0)).ln() + 140.0 * a + 10.0 * (1e-4f64 / 3.0).ln();
        assert!(logsumexp2(bad_local, e) < e + 1e-6);
        // A full-match local placement on a 10kb row survives and
        // beats E by the prior ratio (the locale's concentration).
        let good_local = -(2.0f64 * (10_000.0 - 149.0)).ln() + 150.0 * a;
        assert!(good_local > e);
        assert!(logsumexp2(good_local, e) > e);
    }

    #[test]
    fn the_marginalization_is_monotone_and_preserves_ties() {
        // The L4-control preservation property: logsumexp2 is monotone
        // in the local branch, so a candidate whose local branch is
        // at least every rival's on every unit keeps its per-unit
        // ordering after the marginalization; identical local branches
        // stay identical (the identical-fold non-votes stay non-votes).
        let e = elsewhere_log_prob(150, (1.0 - 1e-4f64).ln(), 24_400_000.0);
        let a_local = -9.9f64;
        let b_local = -25.0f64;
        assert!(logsumexp2(a_local, e) > logsumexp2(b_local, e));
        assert_eq!(logsumexp2(a_local, e), logsumexp2(a_local, e));
        // Deep below E the branches are indistinguishable: the read
        // does not vote.
        let x = logsumexp2(-1546.0f64, e);
        let y = logsumexp2(-2000.0f64, e);
        assert!((x - y).abs() < 1e-12);
        assert!((x - e).abs() < 1e-12);
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
    #[test]
    fn canonical_scheme_selects_per_window_form() {
        // The canonical-scheme selection (slice D's frame repair): per
        // position the frame that spells the k-mer's canonical form
        // forward decides; a position whose canonical frame did not
        // qualify carries NO step; the chosen frame's node is kept.
        let k = 4u64;
        // A 16 bp sequence: window at 0 is canonical-forward (ACGT),
        // the windows at 6 and 12 are rc-canonical (T leads and the
        // complement of the trailing base loses the comparison).
        let seq = b"ACGTGTTAACGTTTGAC".to_vec();
        let seq_lo = 100u64;
        assert!(impg::genome_inference::mem_records::window_is_canonical_forward(
            &seq[0..4]
        ));
        assert!(!impg::genome_inference::mem_records::window_is_canonical_forward(
            &seq[6..10]
        ));
        let mut forward = BTreeMap::new();
        forward.insert(seq_lo + 0, 11);
        forward.insert(seq_lo + 6, 13); // forward-qualified, rc-canonical
        let mut reverse = BTreeMap::new();
        reverse.insert(seq_lo + 6, -13); // the rc frame's node at the same position
        reverse.insert(seq_lo + 12, -17); // rc-only position, rc-canonical window
        assert!(!impg::genome_inference::mem_records::window_is_canonical_forward(
            &seq[12..16]
        ));
        let steps = canonical_scheme_steps(&forward, &reverse, &seq, seq_lo, k);
        // The canonical-forward position keeps the forward node; the
        // rc-canonical positions take the reverse frame's node (or
        // drop when only the forward frame qualified); the forward-only
        // rc-canonical position is REPLACED by the reverse node here.
        assert_eq!(
            steps,
            vec![(seq_lo + 0, 11, 0u8), (seq_lo + 6, -13, 1u8), (seq_lo + 12, -17, 1u8)]
        );
    }

    #[test]
    fn canonical_scheme_drops_forward_only_rc_canonical_positions() {
        // A position where ONLY the forward frame qualified and the
        // window is rc-canonical carries NO step: no census record
        // (derived under the same scheme) can anchor there.
        let k = 4u64;
        let seq = b"ACGTGTTAACGTTTGAC".to_vec();
        let seq_lo = 0u64;
        let mut forward = BTreeMap::new();
        forward.insert(6, 13);
        let reverse: BTreeMap<u64, i32> = BTreeMap::new();
        let steps = canonical_scheme_steps(&forward, &reverse, &seq, seq_lo, k);
        assert!(steps.is_empty());
        // With the rc frame qualifying there too, the position keeps
        // the rc frame's node.
        let mut reverse = BTreeMap::new();
        reverse.insert(6, -13);
        let steps = canonical_scheme_steps(&forward, &reverse, &seq, seq_lo, k);
        assert_eq!(steps, vec![(6, -13, 1u8)]);
    }

    /// The pre-rebuild per-shift verification (phase 0's code, the
    /// reference for the hoisted union-range verify): one fetch per
    /// (occurrence x shift), with the fetch closure carrying the
    /// fetch_seq/Sources clamping semantics (high bound clamped to the
    /// path length; empty crop when the low bound reaches it).
    fn old_verify_key_at(
        shape: &KeyShape,
        occ: &CensusOccurrence,
        reads: &[Vec<u8>],
        fetch: &dyn Fn(u64, u64) -> Vec<u8>,
        k: u64,
        shift: i64,
    ) -> bool {
        let lo = (occ.start as i64 + shift - 1).max(0) as u64;
        let hi = occ.start + shape.span + 1;
        let seq = fetch(lo, hi);
        let origin = occ.start as i64 + shift;
        let read = &reads[shape.read];
        let forward = physical_forward(occ.orientation, shape.mirrored);
        for (j, &(node, canonical_pos)) in shape.canonical.iter().enumerate() {
            let _ = node;
            let window_lo =
                anchor_window_start(canonical_pos, origin as u64, shape.span, k, occ.orientation);
            let rel_lo = window_lo as i64 - lo as i64;
            if rel_lo < 0 || rel_lo + k as i64 > seq.len() as i64 {
                return false;
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
                return false;
            }
        }
        true
    }

    #[test]
    fn hoisted_verify_reproduces_the_per_shift_fetch_verification() {
        // Phase 1, lever 1: the once-per-(occurrence, candidate) union
        // range fetch must reproduce, for every origin shift, exactly
        // the verification the three per-shift fetches performed —
        // including the edge cases where the clamped crop is short or
        // empty (occurrence start near 0; range overhanging the path
        // end), both orientations, both mirror states, and verdicts of
        // both true and false.
        let k = 4u64;
        let path: Vec<u8> = b"ACGTACGTACGTTTGACACGTGCATATCGGATTCGCA".to_vec();
        let path_len = path.len() as u64;
        let fetch = |lo: u64, hi: u64| -> Vec<u8> {
            let hi = hi.min(path_len);
            if lo >= hi {
                Vec::new()
            } else {
                path[lo as usize..hi as usize].to_vec()
            }
        };
        // anchors at canonical positions 0 and 8, own k-mers at read
        // positions 2 and 10; span 8 - 0 + 4 = 12.
        let make_shape = |mirrored: bool| KeyShape {
            canonical: vec![(7, 0), (9, 8)],
            span: 12,
            read: 0,
            mirrored,
            own_positions: vec![2, 10],
        };
        let mut read = path[2..16].to_vec();
        let reads_match = vec![read.clone()];
        read[3] = match read[3] {
            b'A' => b'T',
            other => b'A' + (other == b'A') as u8,
        };
        let reads_mismatch = vec![read];
        // Interior starts verify true (orientation 0, unmirrored);
        // near-zero starts exercise the low-bound clamp, path-end
        // starts exercise the high-bound clamp and the empty crop.
        for start in [0u64, 1, 2, 3, 5, 20, 30, 33, 34, 36, 38, 39, 41] {
            for orientation in [0u8, 1] {
                for mirrored in [false, true] {
                    let occ = CensusOccurrence {
                        path: 0,
                        start,
                        orientation,
                        partitions: vec![],
                        intervals: vec![],
                    };
                    let shape = make_shape(mirrored);
                    for reads in [&reads_match, &reads_mismatch] {
                        // The union range fetch (the hoist): one fetch
                        // covering all three shifts' ranges.
                        let union_lo = (occ.start as i64 - 2).max(0) as u64;
                        let union_hi = occ.start + shape.span + 1;
                        let seq = fetch(union_lo, union_hi);
                        for shift in [-1i64, 0, 1] {
                            assert_eq!(
                                verify_key_at(&shape, &occ, reads, &seq, union_lo, k, shift),
                                old_verify_key_at(&shape, &occ, reads, &fetch, k, shift),
                                "hoisted verify disagrees with the per-shift fetch \
                                 (start {start}, orientation {orientation}, mirrored \
                                 {mirrored}, shift {shift})",
                            );
                        }
                    }
                }
            }
        }
        // And a true verdict exists (the equivalence is not vacuous):
        // the read spells path[2..16], its first anchor's own k-mer at
        // read position 2, so the occurrence start is 4.
        let occ = CensusOccurrence {
            path: 0,
            start: 4,
            orientation: 0,
            partitions: vec![],
            intervals: vec![],
        };
        let shape = make_shape(false);
        let union_lo = 0u64;
        let seq = fetch(union_lo, occ.start + shape.span + 1);
        assert!(verify_key_at(&shape, &occ, &reads_match, &seq, union_lo, k, 0));
        assert!(!verify_key_at(&shape, &occ, &reads_mismatch, &seq, union_lo, k, 0));
    }
}
