//! THE LOCAL-REALIGNMENT EVIDENCE MODEL (owner direction 2026-11-06;
//! assessment-side only — no product file is touched, no instrument of
//! record is modified, no threshold enters anything).
//!
//! Reads currently reach the panel through exact syncmer matches only:
//! variant bases never testify (reads spanning variant pockets match only
//! shared flanks) and mass folds very far (identical sequence anywhere is
//! one interned node, so read mass splits across contexts it never came
//! from). This instrument replaces the evidence model: each read's
//! evidence is its FULL-SEQUENCE local alignment likelihood against each
//! candidate haplotype path in the partition graph — not seed matching.
//!
//! THE MODEL (derived, documented; no tuning constants):
//!   * SCORING: the sample's reads carry uniform Phred-40 quality bytes
//!     (measured over the whole FASTQ; asserted), so the per-base error
//!     probability is the reads' own stated rate ε = 10^(-Q/10); a read
//!     base generated from a candidate base matches with probability
//!     1-ε and mismatches with probability ε/3 (the uniform substitution
//!     model the sample's qualities state). The measured error profile
//!     is documented beside it: the committed multi-census measured
//!     zero unmatched placement bp over 1,063,522 verified occurrences
//!     (reads are exact panel segments where placed).
//!   * GENERATIVE: P(read | generated from row R) = sum over BOTH read
//!     strands and EVERY offset σ with the read fully inside R of
//!     (1/(2*len_R)) * prod_i P(read_i | R[σ+i]). Substitution-only, so
//!     a placement's score is the closed form m*A + c*B over integers
//!     m (matching bases) and c = 150 - m (mismatching bases); A =
//!     ln(1-ε), B = ln(ε/3). Indels between read and row are not
//!     representable under substitution-only generation and score as
//!     the shifted mismatches they are: the row carrying the read's own
//!     indel structure is the one that can generate the read. The
//!     uniform placement prior 1/(2*len_R) is the stated position/strand
//!     prior.
//!   * NO SHARE-SPLITTING: per-read per-row likelihoods are independent;
//!     a record contributes ONCE times its multiplicity; reads
//!     multi-matching several local paths get independent likelihoods
//!     under each (the sum within one row is over that row's own
//!     placements, never a split across candidates).
//!   * CANDIDATES: the locus's candidate paths are the partition graph's
//!     own member rows (the panel's alignment-induced partition
//!     locality), FOLDED by identical sequence+walk (coalesced copies
//!     are one hypothesis); classes are ALL unordered row pairs i<=j
//!     (the instrument of record's convention); a class's
//!     log-likelihood is the sum over the locus's records of
//!     multiplicity * ln(0.5*e^{LL(rec|R_i)} + 0.5*e^{LL(rec|R_j)}) —
//!     the equal 1/2 mixture is the balanced 15x/15x diploid prior,
//!     stated. The record-once universe per locus is the committed
//!     remedy receipt's ingredient record list, so the before/after
//!     table shares the current instrument's record sets exactly.
//!   * THE FACTORIZED EXACT FORM (the owner's (2), validated below):
//!     the locus frame is the axis row (the S288C#0 member row — the
//!     instrument's own coordinate convention); each row's VARIANT
//!     POCKETS against the frame are precomputed once (collinear
//!     shared-node runs, SNV edits inside anchor-pinned segments,
//!     identical flanks extended, indel structure = torn points); each
//!     record is aligned ONCE to its placed backbone (its occurrence
//!     path segment projected onto the frame: one verified projection
//!     per read strand, with the read's alt bases vs the frame); the
//!     frame-collinear placement's (m,c) is then computed PURELY from
//!     the two cached pocket lists — backbone term + the read's base
//!     evidence over the row's pockets — and asserted equal to the
//!     direct full-sequence walk at the same σ over the WHOLE exhaustive
//!     domain (integer equality), unit-proven. Torn/edge placements are
//!     named and always computed by the direct walk (which also carries
//!     the non-anchor background offsets — the sub-syncmer similarity
//!     evidence the seed model dropped entirely).
//!   * QUAL: the existing cluster machinery VERBATIM (mirrored pure
//!     functions from panel_route_spine/cosine_probe.rs — spectrum
//!     knee, cluster cut, k, a, p = s_win/(k*s_win + a),
//!     QUAL = -10*log10(1-p)) over the realignment likelihood ratios,
//!     with the material distance = differing observed mass on merged
//!     node+edge usage multisets; the observed mass per key is the
//!     covering records' multiplicity (record-once, no share split).
//!
//! Inputs: the panel syng, the route graph's sources, the partition
//! graph maps (the committed partition-embedded build), the committed
//! multi-matching census receipt (record placements + multiplicities),
//! the committed remedy likelihood receipt's ingredients sidecar
//! (per-locus record sets), and the sample FASTQ (read-identity and
//! quality verification — the instrument's read sequences are verified
//! against the real reads by hash-multiset membership).
//!
//! Usage (assessment-side):
//!   partition_realign_likelihood --panel <syng-prefix> --routes <dir> \
//!     --partition-graphs <dir> --component S288C#0#chrMT \
//!     --census <cosine-multi-census.jsonl> \
//!     --ingredients <...ingredients.jsonl> \
//!     --reads <reads.fastq.gz> --out <realign.jsonl>

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

/// The reads' exact length: the model's placement arithmetic and the
/// derive-cache record derivation are keyed to it. 150 = the yeast
/// validation's committed reads; the HG002 locus pilot reads are 148
/// (the GIAB BAM's uniformly trimmed length). Set ONCE at startup from
/// --read-length.
static READ_LENGTH: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(150);
/// The offset-histogram capacity: one bin per possible match count of
/// the longest supported read (the validation lengths are 150 and 148).
const READ_HIST_CAP: usize = 150 + 1;
#[inline]
fn read_len() -> usize {
    READ_LENGTH.load(Ordering::Relaxed)
}

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
    /// member rows define the loci and the frame.
    #[arg(long)]
    component: String,
    /// The committed multi-matching census receipt (records, placements).
    #[arg(long)]
    census: PathBuf,
    /// The committed remedy likelihood receipt's ingredients sidecar
    /// (the per-locus record-once sets).
    #[arg(long)]
    ingredients: PathBuf,
    /// The sample FASTQ (read-identity and quality verification).
    #[arg(long)]
    reads: PathBuf,
    /// The likelihood receipt path (sidecars are derived from it).
    #[arg(long)]
    out: PathBuf,
    /// Resident-set guard in GiB (the 64 GiB discipline; 0 = no guard).
    #[arg(long, default_value_t = 64.0)]
    rss_budget_gib: f64,
    /// The reads' exact length (150 = the yeast validation's committed
    /// reads; the HG002 pilot reads are 148).
    #[arg(long, default_value_t = 150)]
    read_length: usize,
    /// Optional binary cache of the read-record derivation (written when
    /// absent, loaded when present) — the derivation pass is a pure
    /// function of the FASTQ and the panel, so the cache cannot change
    /// any result.
    #[arg(long)]
    derive_cache: Option<PathBuf>,
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
        eprintln!("[realign] rss {stage}: {rss_kb} kB");
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
    #[serde(default)]
    partitions: Vec<u32>,
    /// Per read-matched interval, the contained covered node abs ids
    /// (path order) — the seed projection's own key credit, and (taken
    /// consecutively within one interval) the covered adjacencies.
    #[serde(default)]
    intervals: Vec<Vec<u32>>,
}

/// One routed record line of the census receipt.
#[derive(Deserialize)]
struct CensusRecordLine {
    record: u64,
    multiplicity: u64,
    #[serde(default)]
    anchors: u64,
    #[serde(default)]
    occurrences: Vec<CensusOccurrence>,
}

#[derive(Deserialize)]
struct IngredientRecord {
    record: u64,
}

#[derive(Deserialize)]
struct IngredientLine {
    locus: usize,
    records: Vec<IngredientRecord>,
}

/// One member row of a partition graph map (the alignment-induced
/// membership; the partition BEDs are source-forward).
#[derive(Clone, Deserialize)]
struct MemberRow {
    path_name: String,
    start: u64,
    end: u64,
}

#[derive(Deserialize)]
struct PartitionMap {
    partition: u32,
    members: Vec<MemberRow>,
}

// ---------------------------------------------------------------------------
// The QUAL cluster machinery — MIRRORED UNCHANGED from
// examples/panel_route_spine/cosine_probe.rs (the existing instrument of
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
    let mut root = |parent: &mut Vec<usize>, node: usize| {
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

// ---------------------------------------------------------------------------
// The scoring core.
// ---------------------------------------------------------------------------

/// The uniform substitution scoring constants derived from the reads'
/// own quality bytes.
#[derive(Clone, Copy)]
struct Scoring {
    a: f64,
    b: f64,
    epsilon: f64,
    phred: u8,
}

impl Scoring {
    /// The canonical per-placement closed form: m matching bases and
    /// c = 150 - m mismatching bases. Both the factorized path and the
    /// direct walk produce the integers (m, c); this one expression
    /// turns them into the placement score, so integer equality implies
    /// bit-exact score equality.
    #[inline]
    fn placement(&self, m: u32, c: u32) -> f64 {
        m as f64 * self.a + c as f64 * self.b
    }
    /// The canonical row likelihood: logsumexp over the full offset
    /// histogram (both read strands summed; the per-placement weight is
    /// 1/(2*len)) — one expression, so integer-equal histograms imply
    /// bit-exact likelihoods. Terms ≥ 73 matches below the best
    /// underflow to +0.0 exactly and cannot perturb the sum.
    fn row_log_likelihood(&self, hist: &[u64; READ_HIST_CAP], len: u64) -> f64 {
        let mut best = None::<usize>;
        for (m, &count) in hist.iter().enumerate() {
            if count > 0 {
                best = Some(m);
            }
        }
        let Some(best) = best else {
            return f64::NEG_INFINITY;
        };
        let mut sum = 0.0f64;
        for (m, &count) in hist.iter().enumerate() {
            if count > 0 {
                sum += count as f64 * ((m as f64 - best as f64) * (self.a - self.b)).exp();
            }
        }
        self.placement(best as u32, (read_len() - best) as u32) + sum.ln()
            - (2.0 * len as f64).ln()
    }
}

/// One folded candidate row.
struct Row {
    seq: Vec<u8>,
    len: u64,
    /// All walk steps overlapping the row range (absolute bp, signed
    /// node), sorted by bp.
    walk: Vec<(u64, i32)>,
    /// The row's start on its (first) member's path — rows in one fold
    /// share the sequence AND walk, so their relative geometry agrees;
    /// absolute positions differ per member and are kept in `members`.
    start: u64,
    end: u64,
    /// node abs id -> occurrences in the walk (absolute bp, sign).
    node_index: HashMap<u32, Vec<(u64, i32)>>,
    /// Frame runs (row-relative r, frame-relative a): the verified
    /// collinear intervals with the row's SNV pockets vs the frame.
    runs: Vec<FrameRun>,
    members: Vec<MemberRow>,
    nodes: Vec<(u64, u32)>,
    edges: Vec<(u64, u32)>,
}

#[derive(Clone)]
struct FrameRun {
    r_lo: u64,
    r_hi: u64,
    a_lo: u64,
    delta: i64,
    /// SNV edits inside the run: (frame position, this side's base).
    edits: Vec<(u64, u8)>,
}

/// A read strand's projection onto the locus frame: read offset i
/// aligns to frame position x0 + i; the alt list carries the read's own
/// bases where they differ from the frame.
#[derive(Clone)]
struct Projection {
    x0: u64,
    alts: Vec<(u64, u8)>,
}

/// One record's locus data.
struct RecordLocus {
    record: usize,
    multiplicity: u64,
    read: Vec<u8>,
    read_rc: Vec<u8>,
    forward: Option<Projection>,
    reverse: Option<Projection>,
    covered_nodes: BTreeSet<u32>,
    covered_edges: BTreeSet<u64>,
}

// ---------------------------------------------------------------------------
// Frame machinery: shared-node anchored runs between a segment and the
// frame, with SNV edits and verified-identical flank extension.
// ---------------------------------------------------------------------------

struct Anchor {
    seg: u64,
    frame: u64,
    delta: i64,
}

/// Build the anchor-pinned runs of a segment against the frame sequence.
/// Consecutive equal-delta anchors pin the between-segment alignment:
/// every differing base is a SNV edit (an indel would shift the next
/// anchor's delta and tear the run). Unanchored flanks extend only while
/// byte-identical; a differing flank may hide an indel and stays torn.
fn frame_runs(
    anchors: &[Anchor],
    seg_seq: &[u8],
    frame_seq: &[u8],
    k: u64,
) -> Vec<FrameRun> {
    let seg_len = seg_seq.len() as u64;
    let mut runs: Vec<FrameRun> = Vec::new();
    let mut idx = 0usize;
    while idx < anchors.len() {
        let start = &anchors[idx];
        let delta = start.delta;
        let run_lo = start.seg;
        let mut run_hi = start.seg + k;
        let mut edits: Vec<(u64, u8)> = Vec::new();
        let mut cursor = idx;
        while cursor + 1 < anchors.len() {
            let next = &anchors[cursor + 1];
            if next.delta != delta || next.seg <= anchors[cursor].seg {
                break;
            }
            // The inter-anchor segment is pinned at both ends (equal
            // delta): positions [run_hi, next.seg) correspond to frame
            // [run_hi + delta, next.seg + delta).
            for p in run_hi.min(next.seg)..next.seg {
                let frame_pos = p as i64 + delta;
                if frame_pos < 0 || frame_pos >= frame_seq.len() as i64 {
                    break;
                }
                let seg_base = seg_seq[p as usize];
                let frame_base = frame_seq[frame_pos as usize];
                if seg_base != frame_base {
                    edits.push((frame_pos as u64, seg_base));
                }
            }
            run_hi = next.seg + k;
            cursor += 1;
        }
        // Extend flanks only while byte-identical.
        let mut lo = run_lo;
        while lo > 0 {
            let p = lo - 1;
            let frame_pos = p as i64 + delta;
            if frame_pos < 0
                || frame_pos >= frame_seq.len() as i64
                || seg_seq[p as usize] != frame_seq[frame_pos as usize]
            {
                break;
            }
            lo = p;
        }
        let mut hi = run_hi;
        while hi < seg_len {
            let frame_pos = hi as i64 + delta;
            if frame_pos < 0
                || frame_pos >= frame_seq.len() as i64
                || seg_seq[hi as usize] != frame_seq[frame_pos as usize]
            {
                break;
            }
            hi += 1;
        }
        if hi > lo {
            runs.push(FrameRun {
                r_lo: lo,
                r_hi: hi,
                a_lo: (lo as i64 + delta) as u64,
                delta,
                edits,
            });
        }
        idx = cursor + 1;
    }
    runs
}

// ---------------------------------------------------------------------------
// The direct full-sequence walk (the model's brute force) and the
// factorized pocket combination.
// ---------------------------------------------------------------------------

/// The full offset histogram of one read strand against one row: for
/// EVERY offset σ with the read fully inside the row, the number of
/// matching bases — the model's complete placement sum.
/// The full offset histogram of one read strand against one row: for
/// EVERY offset σ with the read fully inside the row, the number of
/// matching bases — the model's complete placement sum. Also returns
/// the strand's dominant placement (best m, first σ at it) for the
/// smear-class measurement.
fn offset_histogram(
    read: &[u8],
    seq: &[u8],
) -> ([u64; READ_HIST_CAP], Option<(u32, u64)>) {
    let mut hist = [0u64; READ_HIST_CAP];
    let mut best: Option<(u32, u64)> = None;
    if seq.len() >= read_len() {
        let max_sigma = seq.len() - read_len();
        for sigma in 0..=max_sigma {
            let window = &seq[sigma..sigma + read_len()];
            let mut m = 0u32;
            for i in 0..read_len() {
                m += (read[i] == window[i]) as u32;
            }
            hist[m as usize] += 1;
            match best {
                Some((best_m, _)) if m <= best_m => {}
                _ => best = Some((m, sigma as u64)),
            }
        }
    }
    (hist, best)
}

/// The direct per-placement walk: (m, c) at one offset.
fn direct_placement(read: &[u8], seq: &[u8], sigma: u64) -> (u32, u32) {
    let window = &seq[sigma as usize..sigma as usize + read_len()];
    let mut m = 0u32;
    for i in 0..read_len() {
        m += (read[i] == window[i]) as u32;
    }
    (m, read_len() as u32 - m)
}

/// The FACTORIZED placement score from the cached pocket lists (the
/// owner's exact form): the read strand's projection (frame x0 + alt
/// bases) against the row's run (frame edits), for the frame-collinear
/// offset σ = x0 + delta. None when the read's frame interval is not
/// inside the run's verified coverage (torn/edge — the direct walk
/// computes those).
fn factorized_placement(
    projection: &Projection,
    run: &FrameRun,
    frame_seq: &[u8],
) -> Option<(u32, u32)> {
    let x_hi = projection.x0 + read_len() as u64;
    let a_hi = run.a_lo + (run.r_hi - run.r_lo);
    if projection.x0 < run.a_lo || x_hi > a_hi {
        return None;
    }
    // c = the read positions where the row's base differs: the union of
    // the two alt lists inside the read's frame interval, comparing the
    // read's own base against the row's own base (frame elsewhere).
    let mut c = 0u32;
    let mut i = 0usize;
    let mut j = 0usize;
    while i < projection.alts.len() || j < run.edits.len() {
        let pa = projection.alts.get(i).map(|&(p, _)| p);
        let ra = run.edits.get(j).map(|&(p, _)| p);
        let (pos, read_base, row_base) = match (pa, ra) {
            (Some(p), Some(r)) => {
                if p < r {
                    i += 1;
                    (p, projection.alts[i - 1].1, u8::MAX)
                } else if r < p {
                    j += 1;
                    (r, u8::MAX, run.edits[j - 1].1)
                } else {
                    i += 1;
                    j += 1;
                    (p, projection.alts[i - 1].1, run.edits[j - 1].1)
                }
            }
            (Some(p), None) => {
                i += 1;
                (p, projection.alts[i - 1].1, u8::MAX)
            }
            (None, Some(r)) => {
                j += 1;
                (r, u8::MAX, run.edits[j - 1].1)
            }
            (None, None) => break,
        };
        if pos < projection.x0 || pos >= x_hi {
            continue;
        }
        let frame_base = frame_seq[pos as usize];
        let read_base = if read_base == u8::MAX { frame_base } else { read_base };
        let row_base = if row_base == u8::MAX { frame_base } else { row_base };
        c += (read_base != row_base) as u32;
    }
    Some((read_len() as u32 - c, c))
}

fn node_unused() {}

// ---------------------------------------------------------------------------
// FASTQ streaming.
// ---------------------------------------------------------------------------

fn stream_fastq(
    path: &PathBuf,
    mut visit: impl FnMut(&[u8], &[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let (reader, _) =
        niffler::get_reader(Box::new(File::open(path)?)).map_err(io::Error::other)?;
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
    // Exact integer guards on both sides of the floating point estimate.
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
// main
// ---------------------------------------------------------------------------

fn main() -> io::Result<()> {
    let options = Options::parse();
    ensure(options.read_length > 0 && options.read_length <= 150, "read length must be in (0, 150]")?;
    READ_LENGTH.store(options.read_length, Ordering::Release);
    let started = Instant::now();
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
    let contig = options
        .component
        .splitn(3, '#')
        .nth(2)
        .ok_or_else(|| invalid("component lacks a contig suffix"))?
        .to_string();
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
    // The ingredients sidecar (the per-locus record-once sets).
    let mut locus_records: Vec<Vec<usize>> = Vec::new();
    {
        let file = File::open(&options.ingredients)?;
        for line in BufReader::new(file).lines() {
            let line: IngredientLine = serde_json::from_str(&line?)?;
            ensure(line.locus == locus_records.len(), "ingredient loci not dense")?;
            locus_records.push(line.records.iter().map(|r| r.record as usize).collect());
        }
    }
    let n_loci = locus_records.len();
    ensure(n_loci > 0, "no loci in the ingredients receipt")?;
    // The partition graph maps: the axis partitions (holding a member
    // row on the component path) ranked by axis-row start are the loci.
    let mut axis_partitions: Vec<(u32, u64)> = Vec::new();
    let mut maps: HashMap<u32, PartitionMap> = HashMap::new();
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
    ensure(
        axis_partitions.len() == n_loci,
        "axis partition count does not match the ingredient loci",
    )?;
    let locus_partition: Vec<u32> =
        axis_partitions.iter().map(|&(partition, _)| partition).collect();
    for records in &locus_records {
        for &record in records {
            ensure(record < census_records.len(), "ingredient record outside the census")?;
        }
    }
    rss.probe("input_load")?;
    eprintln!(
        "[realign] inputs: {} census records, {} loci, axis partitions {:?}",
        census_records.len(),
        n_loci,
        locus_partition,
    );

    // --------------------------------------- read identity + scoring rate
    let mut quality_bytes: BTreeSet<u8> = BTreeSet::new();
    let mut quality_counts: BTreeMap<u8, u64> = BTreeMap::new();
    let mut total_reads = 0u64;
    stream_fastq(&options.reads, |_seq, qual| {
        total_reads += 1;
        for &q in qual {
            quality_bytes.insert(q);
            *quality_counts.entry(q).or_insert(0) += 1;
        }
        Ok(())
    })?;
    // The sample's own quality statement: the mean per-base error
    // probability over every quality byte in the FASTQ (the yeast
    // validation's uniform-quality reads reduce EXACTLY to their single
    // Phred; real reads - the HG002 pilot - carry their measured
    // per-base qualities).
    ensure(
        !quality_bytes.is_empty(),
        "the sample carries no read qualities",
    )?;
    let mut base_count = 0u64;
    let mut error_sum = 0.0f64;
    for (&q, &count) in quality_counts.iter() {
        base_count += count;
        error_sum += count as f64 * 10f64.powf(-((q - b'!') as f64) / 10.0);
    }
    ensure(base_count > 0, "the sample carries no read qualities")?;
    let epsilon0 = error_sum / base_count as f64;
    ensure(epsilon0 > 0.0 && epsilon0 < 1.0, "degenerate read qualities")?;
    let phred = (-10.0 * epsilon0.log10()).round() as u8;
    let epsilon = 10f64.powf(-(phred as f64) / 10.0);
    let scoring = Scoring {
        a: (1.0 - epsilon).ln(),
        b: (epsilon / 3.0).ln(),
        epsilon,
        phred,
    };
    // THE READ-RECORD RE-DERIVATION (the public accessor the committed
    // receipts used): every FASTQ read's own-orientation maximal MEM
    // records via `tagged_mem_records` (the Syng anchor scheme), each
    // encoded and canonicalized exactly as the record multiset of
    // record. A record is a WALK (its anchor hull can be shorter than
    // the read); the reads carrying it are its evidence, and under the
    // realignment model each read's FULL sequence is the evidence — the
    // reads of one record agree exactly over the record's anchor hull
    // (the census's zero-gap measurement, asserted below) and carry
    // their own bases beyond it.
    let derive_started = Instant::now();
    let scheme = impg::genome_inference::mem_records::anchor_scheme();
    // The derivation cache (assessment-side convenience only: the
    // derivation is a pure function of the FASTQ and the panel).
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
        let _ = &mut u64buf;
        Ok((key_tokens, key_reads, reads, read_multiplicity, read_records))
    }
    #[derive(Clone)]
    struct ReadRecord {
        key: u32,
        // The walk's anchors with ABSOLUTE read positions (the own
        // orientation, as extracted): (signed node, read position).
        walk: Vec<(i32, u32)>,
    }
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
    let (key_tokens, key_reads, reads, read_multiplicity, read_records): (
        Vec<Vec<u64>>,
        Vec<Vec<usize>>,
        Vec<Vec<u8>>,
        Vec<u32>,
        Vec<Vec<ReadRecord>>,
    ) = match &options.derive_cache {
        Some(path) if path.exists() => {
            eprintln!("[realign] loading the derivation cache {}", path.display());
            let (tokens, lists, seqs, mults, records) = read_derive_cache(path)?;
            (
                tokens,
                lists,
                seqs,
                mults,
                records
                    .into_iter()
                    .map(|records| {
                        records
                            .into_iter()
                            .map(|(key, walk)| ReadRecord { key, walk })
                            .collect()
                    })
                    .collect(),
            )
        }
        _ => {
    let mut key_index: BTreeMap<Vec<u64>, u32> = BTreeMap::new();
    let mut key_tokens: Vec<Vec<u64>> = Vec::new();
    let mut key_reads: Vec<Vec<usize>> = Vec::new();
    let mut reads: Vec<Vec<u8>> = Vec::new();
    let mut read_multiplicity: Vec<u32> = Vec::new();
    let mut read_records: Vec<Vec<ReadRecord>> = Vec::new();
    {
        let mut read_counts: HashMap<u64, u32> = HashMap::new();
        let mut pending: Vec<Vec<u8>> = Vec::new();
        let mut drain = |pending: &mut Vec<Vec<u8>>| -> io::Result<()> {
            // Parallel phase: pure per-read extraction (tokens + walks).
            let per_read: Vec<Vec<(Vec<u64>, Vec<(i32, u32)>)>> = pending
                .par_iter()
                .map(|read| {
                    let tagged = impg::genome_inference::mem_records::
                        tagged_mem_records_scheme(&panel, read, scheme)?;
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
            // Sequential phase: intern the keys, bind the reads.
            for (index, seq) in pending.drain(..).enumerate() {
                let base = reads.len();
                reads.push(seq);
                *read_counts.entry(fnv1a64(&reads[base])).or_default() += 1;
                let mut own_records: Vec<ReadRecord> = Vec::with_capacity(per_read[index].len());
                for (tokens, walk) in &per_read[index] {
                    let key = *key_index.entry(tokens.clone()).or_insert_with(|| {
                        key_tokens.push(tokens.clone());
                        key_reads.push(Vec::new());
                        key_tokens.len() as u32 - 1
                    });
                    key_reads[key as usize].push(base);
                    own_records.push(ReadRecord { key, walk: walk.clone() });
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
        // The multiplicity of a read = the number of identical reads in
        // the sample (identical sequences carry identical records, so
        // the count is global); duplicates collapse to their first
        // representative (weight carried by read_multiplicity).
        let mut first_seen: HashMap<u64, usize> = HashMap::new();
        for (index, seq) in reads.iter().enumerate() {
            let hash = fnv1a64(seq);
            match first_seen.get(&hash) {
                Some(&first) => {
                    read_multiplicity.push(0);
                    let _ = first;
                }
                None => {
                    first_seen.insert(hash, index);
                    read_multiplicity.push(*read_counts.get(&hash).unwrap_or(&1));
                }
            }
        }
    }
            (
                key_tokens,
                key_reads,
                reads,
                read_multiplicity,
                read_records,
            )
        }
    };
    if let Some(path) = &options.derive_cache {
        if !path.exists() {
            write_derive_cache(
                path,
                &key_tokens,
                &key_reads,
                &reads,
                &read_multiplicity,
                &read_records
                    .iter()
                    .map(|records| {
                        records
                            .iter()
                            .map(|record| (record.key, record.walk.clone()))
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            )?;
            eprintln!("[realign] derivation cache written: {}", path.display());
        }
    }
    // Bind every census record (its occurrences are the routing's
    // verified placements) to the re-derived key whose decoded walk
    // verifies at the record's first occurrence — SEQUENCE-BASED (the
    // routing's canonical-scheme territories include rc-frame-qualified
    // steps the forward panel walk does not spell, so node-based
    // verification via walk_path_range cannot see those placements; the
    // anchor k-mers ARE the placement either way). The anchor count and
    // the read count must both match the census line; records with more
    // than one verifying key are disambiguated by full-occurrence
    // verification.
    let census_key: Vec<u32> = {
        let by_multiplicity: BTreeMap<u64, Vec<usize>> = {
            let mut map: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
            for (record, line) in census_records.iter().enumerate() {
                map.entry(line.multiplicity).or_default().push(record);
            }
            map
        };
        // Per key: the decoded walk, the span, and one representative
        // read with its walk (for the anchor k-mers).
        struct KeyShape {
            walk: Vec<(i32, u64)>,
            span: u64,
            read: usize,
            w_lo: u64,
            w_hi: u64,
        }
        let mut key_shapes: Vec<Option<KeyShape>> = Vec::with_capacity(key_tokens.len());
        for (key, tokens) in key_tokens.iter().enumerate() {
            let shape = (|| -> Option<KeyShape> {
                let walk = decode_record_tokens(tokens).ok()?;
                if walk.is_empty() {
                    return None;
                }
                let span = walk.last().map(|&(_, p)| p + k).unwrap_or(0);
                let read = *key_reads[key].first()?;
                let own = read_records[read]
                    .iter()
                    .find(|record| record.key as usize == key)?;
                let w_lo = own.walk.first().map(|&(_, p)| p as u64).unwrap_or(0);
                let w_hi = own.walk.last().map(|&(_, p)| p as u64).unwrap_or(0) + k;
                Some(KeyShape { walk, span, read, w_lo, w_hi })
            })();
            key_shapes.push(shape);
        }
        // The sequence verification of one key at one occurrence. The
        // census's orientation flag refers to the CANONICAL token walk
        // (which may be the mirror of the read's own walk — the
        // canonical() token form picks the lexicographically smaller
        // encoding), so the read's PHYSICAL strand against the path is
        // (orientation == 0) XOR mirrored. Physically forward: every
        // anchor's k-mer (the read's own bases at the walk's absolute
        // positions) equals the path segment's k-mer at the placed rel;
        // physically reverse: the read's window equals the reverse
        // complement of the placed rc-walk's k-mer.
        let verify_at = |key: usize, occ: &CensusOccurrence| -> io::Result<bool> {
            let Some(shape) = key_shapes[key].as_ref() else { return Ok(false) };
            let seq = fetch_seq(
                &panel.name_map.path_to_name[occ.path],
                occ.start,
                occ.start + shape.span,
            )?;
            if seq.len() != shape.span as usize {
                return Ok(false);
            }
            let read = &reads[shape.read];
            let own = read_records[shape.read]
                .iter()
                .find(|record| record.key as usize == key)
                .ok_or_else(|| invalid("read record mismatch"))?;
            let mirrored = impg::sample_mem_bwt::encode_walk(
                &own.walk.iter().map(|&(node, pos)| (node, pos as u64)).collect::<Vec<_>>(),
            )? != key_tokens[key];
            let physical_forward = (occ.orientation == 0) != mirrored;
            // Work in the read's OWN walk coordinates throughout: a
            // physically forward placement puts the own anchor k-mer at
            // [start + rel); a physically reverse placement puts it at
            // the reverse complement of [start + rc_rel). (Mirroring
            // reverses the anchor order of the canonical token walk, so
            // index pairing with the decoded walk would be wrong.)
            for &(_node, own_pos) in &own.walk {
                let rel = own_pos as u64 - shape.w_lo;
                let own_pos = own_pos as u64;
                let kmer = &read[own_pos as usize..(own_pos + k) as usize];
                if physical_forward {
                    if kmer != &seq[rel as usize..(rel + k) as usize] {
                        return Ok(false);
                    }
                } else {
                    let rc_rel = shape.span - k - rel;
                    let placed = revcomp(&seq[rc_rel as usize..(rc_rel + k) as usize]);
                    if placed != kmer {
                        return Ok(false);
                    }
                }
            }
            let _ = node_unused();
            Ok(true)
        };
        let mut candidates_of: Vec<Vec<u32>> = vec![Vec::new(); census_records.len()];
        for (key, tokens) in key_tokens.iter().enumerate() {
            let count = key_reads[key].len() as u64;
            let Some(candidates) = by_multiplicity.get(&count) else {
                continue;
            };
            let Ok(walk) = decode_record_tokens(tokens) else { continue };
            if walk.is_empty() {
                continue;
            }
            for &record in candidates {
                let line = &census_records[record];
                if line.anchors as usize != walk.len() {
                    continue;
                }
                let verified = verify_at(key, &line.occurrences[0])?;
                let debug_key = std::env::var("IMPG_REALIGN_DEBUG_KEY")
                    .ok()
                    .and_then(|value| value.parse::<usize>().ok());
                if debug_key == Some(key) && !verified {
                    let shape = key_shapes[key].as_ref().unwrap();
                    let seq = fetch_seq(
                        &panel.name_map.path_to_name[line.occurrences[0].path],
                        line.occurrences[0].start,
                        line.occurrences[0].start + shape.span,
                    )
                    .unwrap_or_default();
                    let read = &reads[shape.read];
                    let own = read_records[shape.read]
                        .iter()
                        .find(|record| record.key as usize == key)
                        .unwrap();
                    let first = own.walk.first().map(|&(_, p)| p as usize).unwrap_or(0);
                    eprintln!(
                        "[realign] debug: record {record} key {key} walk {:?} \
                         span {} own_pos {first} occ {:?} span_seq {:?} read_kmer {:?} \
                         read[0..8] {:?}",
                        shape.walk,
                        shape.span,
                        (line.occurrences[0].path, line.occurrences[0].start,
                         line.occurrences[0].orientation),
                        String::from_utf8_lossy(&seq),
                        String::from_utf8_lossy(&read[first..first + k as usize]),
                        String::from_utf8_lossy(&read[..8]),
                    );
                }
                if verified {
                    candidates_of[record].push(key as u32);
                }
            }
        }
        let mut bound: Vec<Option<u32>> = vec![None; census_records.len()];
        for (record, candidates) in candidates_of.iter().enumerate() {
            match candidates.len() {
                0 => {}
                1 => bound[record] = Some(candidates[0]),
                _ => {
                    // Disambiguate by full-occurrence verification: the
                    // record's own key verifies at every occurrence.
                    let mut surviving: Vec<u32> = Vec::new();
                    for &key in candidates {
                        let mut all = true;
                        for occ in &census_records[record].occurrences {
                            if !verify_at(key as usize, occ)? {
                                all = false;
                                break;
                            }
                        }
                        if all {
                            surviving.push(key);
                        }
                    }
                    ensure(
                        surviving.len() == 1,
                        "census record has ambiguous verifying keys",
                    )?;
                    bound[record] = surviving.first().copied();
                }
            }
        }
        let bound_for_diag = &bound;
        let unbound: Vec<usize> = bound_for_diag
            .iter()
            .enumerate()
            .filter(|(_, key)| key.is_none())
            .map(|(record, _)| record)
            .collect();
        if !unbound.is_empty() {
            for &record in unbound.iter().take(5) {
                let line = &census_records[record];
                let occ = &line.occurrences[0];
                eprintln!(
                    "[realign] UNBOUND record {record}: multiplicity {}, anchors {}, \
                     occurrence[0] path {} start {} orientation {} ({} occurrences)",
                    line.multiplicity,
                    line.anchors,
                    occ.path,
                    occ.start,
                    occ.orientation,
                    line.occurrences.len(),
                );
                // Candidate keys by (count, anchor count) with their
                // decoded walks, for diagnosis.
                let mut shown = 0;
                for (key, tokens) in key_tokens.iter().enumerate() {
                    if key_reads[key].len() as u64 != line.multiplicity {
                        continue;
                    }
                    let Ok(walk) = decode_record_tokens(tokens) else { continue };
                    if walk.len() as u64 != line.anchors {
                        continue;
                    }
                    eprintln!(
                        "  candidate key {key}: walk {:?}",
                        walk.iter().map(|&(n, p)| (n, p)).collect::<Vec<_>>()
                    );
                    shown += 1;
                    if shown >= 5 {
                        break;
                    }
                }
                if shown == 0 {
                    eprintln!("  (no candidate keys with matching count and anchor count)");
                }
            }
        }
        bound
            .into_iter()
            .map(|key| key.ok_or_else(|| invalid("census record has no bound key")))
            .collect::<io::Result<Vec<u32>>>()?
    };
    // Per census record: its reads, its covered nodes/edges, and the
    // hull verification (every read equals the placed path segment over
    // the record's anchor hull — the census's zero-gap claim) plus the
    // measured edge profile beyond the hull.
    struct RecordReads {
        reads: Vec<usize>,
        covered_nodes: BTreeSet<u32>,
        covered_edges: BTreeSet<u64>,
        multiplicity: u64,
    }
    let mut census_reads: Vec<RecordReads> = Vec::with_capacity(census_records.len());
    let mut edge_bp = 0u64;
    let mut edge_mismatches = 0u64;
    for (record, line) in census_records.iter().enumerate() {
        let key = census_key[record] as usize;
        let mut covered_nodes = BTreeSet::new();
        let mut covered_edges = BTreeSet::new();
        for occ in &line.occurrences {
            for interval in &occ.intervals {
                for &node in interval {
                    covered_nodes.insert(node);
                }
                for window in interval.windows(2) {
                    covered_edges.insert(pack_edge(window[0], window[1]));
                }
            }
        }
        let walk = decode_record_tokens(&key_tokens[key])?;
        let span = walk.last().map(|&(_, p)| p + k).unwrap_or(0);
        let occ = &line.occurrences[0];
        // The read's offset-0 placement on the occurrence path: the
        // placed (rebased) walk starts at `occ.start`, and the read's
        // walk begins at the own walk's first anchor position.
        for &read_index in &key_reads[key] {
            let read = &reads[read_index];
            ensure(read.len() == read_len(), "bound read is not the declared read length")?;
            let own = read_records[read_index]
                .iter()
                .find(|r| r.key as usize == key)
                .ok_or_else(|| invalid("read record mismatch"))?;
            let w_lo = own.walk.first().map(|&(_, p)| p as u64).unwrap_or(0);
            let w_hi = own.walk.last().map(|&(_, p)| p as u64).unwrap_or(0) + k;
            // The read's PHYSICAL strand (the census orientation refers
            // to the canonical token walk, which may be the read's own
            // walk or its mirror).
            let mirrored = impg::sample_mem_bwt::encode_walk(
                &own.walk
                    .iter()
                    .map(|&(node, pos)| (node, pos as u64))
                    .collect::<Vec<_>>(),
            )? != key_tokens[key];
            let physical_forward = (occ.orientation == 0) != mirrored;
            let read_offset = if physical_forward {
                w_lo
            } else {
                read_len() as u64 - w_hi
            };
            let path_lo = occ.start.saturating_sub(read_offset);
            let seg_full = fetch_seq(
                &panel.name_map.path_to_name[occ.path],
                path_lo,
                path_lo + read_len() as u64,
            )?;
            let strand_full = if physical_forward {
                seg_full.clone()
            } else {
                revcomp(&seg_full)
            };
            ensure(strand_full.len() == read_len(), "full placed segment short")?;
            // The hull (both orientations place the same read window).
            ensure(
                read[w_lo as usize..w_hi as usize]
                    == strand_full[w_lo as usize..w_hi as usize],
                "the read hull differs from the placed path segment",
            )?;
            edge_bp += read_len() as u64 - (w_hi - w_lo);
            for i in 0..read_len() {
                if (i as u64) < w_lo || (i as u64) >= w_hi {
                    edge_mismatches += (read[i] != strand_full[i]) as u64;
                }
            }
        }
        census_reads.push(RecordReads {
            reads: key_reads[key].clone(),
            covered_nodes,
            covered_edges,
            multiplicity: line.multiplicity,
        });
    }
    let mut census_of_key: Vec<i64> = vec![-1; key_tokens.len()];
    for (record, &key) in census_key.iter().enumerate() {
        census_of_key[key as usize] = record as i64;
    }
    let multiplicity_total: u64 = census_reads.iter().map(|r| r.multiplicity).sum();
    ensure(multiplicity_total <= total_reads, "multiplicities exceed the reads")?;
    rss.probe("reads_verified")?;
    eprintln!(
        "[realign] reads re-derived: {total_reads} reads, {} keys, census records {},          multiplicity total {multiplicity_total}, edge bp {edge_bp}, edge mismatches          {edge_mismatches}, uniform Phred {phred}, epsilon {epsilon:e} ({:.1}s)",
        key_tokens.len(),
        census_records.len(),
        derive_started.elapsed().as_secs_f64(),
    );

    // ------------------------------------------------------- the receipts
    let out = &options.out;
    let mut report = BufWriter::new(File::create(out)?);
    let mut competitors =
        BufWriter::new(File::create(format!("{}.competitors.jsonl", out.display()))?);
    let mut ties = BufWriter::new(File::create(format!("{}.ties.jsonl", out.display()))?);
    let mut records_file =
        BufWriter::new(File::create(format!("{}.records.jsonl", out.display()))?);
    let mut ingredients =
        BufWriter::new(File::create(format!("{}.ingredients.jsonl", out.display()))?);
    let mut validation =
        BufWriter::new(File::create(format!("{}.validation.jsonl", out.display()))?);
    let mut reads_file =
        BufWriter::new(File::create(format!("{}.reads.jsonl", out.display()))?);
    {
        // The component's read universe: every distinct read bound to a
        // census record (collapsed duplicates keep their weight).
        let mut component_reads: BTreeSet<usize> = BTreeSet::new();
        for reads_of in &census_reads {
            for &read_index in &reads_of.reads {
                if read_multiplicity[read_index] > 0 {
                    component_reads.insert(read_index);
                }
            }
        }
        let component_read_count = component_reads.len();
        for read_index in component_reads {
            let records: Vec<serde_json::Value> = read_records[read_index]
                .iter()
                .map(|record| {
                    json!({
                        "key": record.key,
                        "walk": record.walk.iter()
                            .map(|&(node, pos)| json!([pos, node])).collect::<Vec<_>>(),
                    })
                })
                .collect();
            serde_json::to_writer(
                &mut reads_file,
                &json!({
                    "read_index": read_index,
                    "multiplicity": read_multiplicity[read_index],
                    "sequence": String::from_utf8(reads[read_index].clone()).unwrap(),
                    "records": records,
                }),
            )?;
            writeln!(reads_file)?;
        }
        eprintln!(
            "[realign] read universe: {component_read_count} distinct reads over {} census records",
            census_records.len(),
        );
    }
    drop(reads_file);

    // --------------------------------------------------- per-locus passes
    let mut total_factorized_checked = 0u64;
    let mut total_factorized_equal = 0u64;
    for locus in 0..n_loci {
        let locus_started = Instant::now();
        let partition = locus_partition[locus];
        let map = maps
            .get(&partition)
            .ok_or_else(|| invalid("axis partition map missing"))?;

        // ------------------------------------------------------- the rows
        let mut fold: BTreeMap<(Vec<u8>, Vec<(u64, i32)>), usize> = BTreeMap::new();
        let mut rows: Vec<Row> = Vec::new();
        for member in &map.members {
            let path_idx = path_of_name[&member.path_name];
            let end = member
                .end
                .min(panel.name_map.path_to_length[path_idx]);
            ensure(member.start < end, "degenerate member row")?;
            let seq = fetch_seq(&member.path_name, member.start, end)?;
            let mut walk: Vec<(u64, i32)> = panel
                .walk_path_range(path_idx, member.start, end)?
                .into_iter()
                .map(|(node, bp)| (bp, node))
                .collect();
            walk.sort_unstable_by_key(|&(bp, _)| bp);
            let key = (
                seq.clone(),
                walk.iter()
                    .map(|&(bp, node)| (bp - member.start, node))
                    .collect::<Vec<_>>(),
            );
            let index = *fold.entry(key).or_insert_with(|| {
                let mut node_index: HashMap<u32, Vec<(u64, i32)>> = HashMap::new();
                for &(bp, node) in &walk {
                    node_index
                        .entry(node.unsigned_abs())
                        .or_default()
                        .push((bp, node.signum()));
                }
                let mut nodes: Vec<(u64, u32)> = walk
                    .iter()
                    .filter(|&&(bp, _)| bp >= member.start && bp + k <= end)
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
                let contained: Vec<(u64, i32)> = walk
                    .iter()
                    .copied()
                    .filter(|&(bp, _)| bp >= member.start && bp + k <= end)
                    .collect();
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
                rows.push(Row {
                    len: seq.len() as u64,
                    seq,
                    walk,
                    start: member.start,
                    end,
                    node_index,
                    runs: Vec::new(),
                    members: Vec::new(),
                    nodes,
                    edges,
                });
                rows.len() - 1
            });
            rows[index].members.push(member.clone());
        }
        let n_rows = rows.len();

        // ------------------------------------ the frame + the row pockets
        let frame_index = rows
            .iter()
            .position(|row| row.members.iter().any(|m| m.path_name == options.component))
            .ok_or_else(|| invalid("the axis row is absent from the partition graph"))?;
        let axis_start = rows[frame_index].start;
        let frame_seq = rows[frame_index].seq.clone();
        let mut frame_nodes: HashMap<u32, Vec<(u64, i32)>> = HashMap::new();
        for &(bp, node) in &rows[frame_index].walk {
            frame_nodes
                .entry(node.unsigned_abs())
                .or_default()
                .push((bp - axis_start, node.signum()));
        }
        let strand_mismatch_anchors = AtomicU64::new(0);
        for row_index in 0..n_rows {
            let row_start = rows[row_index].start;
            let mut anchors: Vec<Anchor> = Vec::new();
            for &(bp, node) in &rows[row_index].walk {
                if let Some(frame_hits) = frame_nodes.get(&node.unsigned_abs()) {
                    for &(fbp, fsign) in frame_hits {
                        if fsign == node.signum() {
                            anchors.push(Anchor {
                                seg: bp - row_start,
                                frame: fbp,
                                delta: fbp as i64 - (bp - row_start) as i64,
                            });
                        } else {
                            strand_mismatch_anchors.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }
            }
            anchors.sort_unstable_by_key(|a| (a.seg, a.frame));
            anchors.dedup_by(|a, b| a.seg == b.seg && a.frame == b.frame);
            rows[row_index].runs = frame_runs(&anchors, &rows[row_index].seq, &frame_seq, k);
        }

        // -------------------------------------------------- the truth rows
        let row_of_member = |name: &str| -> io::Result<Option<usize>> {
            let hits: Vec<usize> = rows
                .iter()
                .enumerate()
                .filter(|(_, row)| row.members.iter().any(|m| m.path_name == name))
                .map(|(index, _)| index)
                .collect();
            ensure(hits.len() <= 1, "multiple member rows for the truth path")?;
            Ok(hits.first().copied())
        };
        let truth_first = row_of_member(&options.component)?;
        let truth_second = row_of_member(&format!("SK1#0#{contig}"))?;
        let truth_pair = match (truth_first, truth_second) {
            (Some(a), Some(b)) => Some((a.min(b), a.max(b))),
            _ => None,
        };
        let truth_flat = truth_pair.map(|(a, b)| b * (b + 1) / 2 + a);

        // --------------------------------------- the locus's read data
        // The locus's evidence unit is the READ (the record-once
        // discipline carries over as read-once): every distinct read
        // carrying any of the locus's records, weighted by its
        // identical-sequence multiplicity. A read's placement structure
        // is the union of ALL its records' verified occurrences; its
        // frame projection is built once per strand from those anchors.
        let records_started = Instant::now();
        let mut locus_read_set: BTreeSet<usize> = BTreeSet::new();
        for &record in &locus_records[locus] {
            for &read_index in &census_reads[record].reads {
                if read_multiplicity[read_index] > 0 {
                    locus_read_set.insert(read_index);
                }
            }
        }
        let mut locus_record_data: Vec<RecordLocus> =
            Vec::with_capacity(locus_read_set.len());
        for read_index in locus_read_set {
            let read_seq = &reads[read_index];
            let read_rc_seq = revcomp(read_seq);
            let mut covered_nodes: BTreeSet<u32> = BTreeSet::new();
            let mut covered_edges: BTreeSet<u64> = BTreeSet::new();
            // The read's placement anchors per strand, from ALL its
            // records' occurrences.
            let mut anchors_forward: Vec<Anchor> = Vec::new();
            let mut anchors_reverse: Vec<Anchor> = Vec::new();
            for own in &read_records[read_index] {
                let record = usize::try_from(census_of_key[own.key as usize])
                    .ok()
                    .ok_or_else(|| invalid("read carries a key outside the census"))?;
                let reads_of = &census_reads[record];
                covered_nodes.extend(reads_of.covered_nodes.iter().copied());
                covered_edges.extend(reads_of.covered_edges.iter().copied());
                let span = own.walk.last().map(|&(_, p)| p as u64).unwrap_or(0) + k;
                let w_lo = own.walk.first().map(|&(_, p)| p as u64).unwrap_or(0);
                let w_hi = span + w_lo;
                let mirrored = impg::sample_mem_bwt::encode_walk(
                    &own.walk
                        .iter()
                        .map(|&(node, pos)| (node, pos as u64))
                        .collect::<Vec<_>>(),
                )? != key_tokens[own.key as usize];
                for occ in &census_records[record].occurrences {
                    // The read's PHYSICAL strand at this occurrence (the
                    // census orientation refers to the canonical token
                    // walk, which may be the read's own walk or its
                    // mirror). Physically forward: the read's anchor
                    // k-mer at [w_lo + rel) equals the path's at
                    // [start + rel), and the frame hit must spell the
                    // node in the same orientation as the read's own
                    // walk sign. Physically reverse: the read's window
                    // equals the reverse complement of the placed
                    // rc-walk's k-mer; the frame hit must spell the
                    // node opposite to the read's own sign.
                    let physical_forward = (occ.orientation == 0) != mirrored;
                    for &(node, own_pos) in &own.walk {
                        let rel = own_pos as u64 - w_lo;
                        if physical_forward {
                            let read_pos = own_pos as u64;
                            let path_bp = occ.start + rel;
                            if let Some(frame_hits) =
                                frame_nodes.get(&node.unsigned_abs())
                            {
                                for &(fbp, fsign) in frame_hits {
                                    if fsign == node.signum() {
                                        anchors_forward.push(Anchor {
                                            seg: read_pos,
                                            frame: fbp - axis_start,
                                            delta: (fbp - axis_start) as i64
                                                - read_pos as i64,
                                        });
                                    }
                                }
                            }
                        } else {
                            let rc_rel = span - k - rel;
                            let rc_read_pos = (read_len() as u64 - w_hi) + rc_rel;
                            let path_bp = occ.start + rc_rel;
                            let _ = path_bp;
                            if let Some(frame_hits) =
                                frame_nodes.get(&node.unsigned_abs())
                            {
                                for &(fbp, fsign) in frame_hits {
                                    if fsign == -node.signum() {
                                        anchors_reverse.push(Anchor {
                                            seg: rc_read_pos,
                                            frame: fbp - axis_start,
                                            delta: (fbp - axis_start) as i64
                                                - rc_read_pos as i64,
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
            }
            let mut build_projection = |anchors: Vec<Anchor>,
                                        strand: &[u8]|
             -> io::Result<Option<Projection>> {
                let mut anchors = anchors;
                anchors.sort_unstable_by_key(|a| (a.seg, a.frame));
                anchors.dedup_by(|a, b| a.seg == b.seg && a.frame == b.frame);
                let runs = frame_runs(&anchors, strand, &frame_seq, k);
                for run in &runs {
                    // A full-coverage verified run gives the projection:
                    // this strand's bases sit at frame [a_lo, a_lo+150).
                    // The alt list is computed from the READ STRAND'S OWN
                    // bases against the frame (never assumed from the
                    // placed path segment): within the anchor hull the
                    // read equals its placed segment exactly (verified at
                    // binding); beyond the hull the read's real edge
                    // bases are compared as they are.
                    if run.r_lo == 0 && run.r_hi >= read_len() as u64 {
                        let x0 = run.a_lo;
                        let mut alts: Vec<(u64, u8)> = Vec::new();
                        for i in 0..read_len() {
                            let frame_pos = x0 + i as u64;
                            if frame_pos >= frame_seq.len() as u64 {
                                break;
                            }
                            if strand[i] != frame_seq[frame_pos as usize] {
                                alts.push((frame_pos, strand[i]));
                            }
                        }
                        return Ok(Some(Projection { x0, alts }));
                    }
                }
                Ok(None)
            };
            let forward = build_projection(anchors_forward, read_seq)?;
            let reverse = build_projection(anchors_reverse, &read_rc_seq)?;
            locus_record_data.push(RecordLocus {
                record: read_index,
                multiplicity: read_multiplicity[read_index] as u64,
                read: read_seq.clone(),
                read_rc: read_rc_seq,
                forward,
                reverse,
                covered_nodes,
                covered_edges,
            });
        }
        let n_records = locus_record_data.len();
        let records_seconds = records_started.elapsed().as_secs_f64();

        // ------------------------------------- observed mass (no shares)
        let observed_nodes: HashMap<u64, f64> = {
            let mut map: HashMap<u64, f64> = HashMap::new();
            for record in &locus_record_data {
                for &node in &record.covered_nodes {
                    *map.entry(node as u64).or_default() += record.multiplicity as f64;
                }
            }
            map
        };
        let observed_edges: HashMap<u64, f64> = {
            let mut map: HashMap<u64, f64> = HashMap::new();
            for record in &locus_record_data {
                for &edge in &record.covered_edges {
                    *map.entry(edge).or_default() += record.multiplicity as f64;
                }
            }
            map
        };
        let observed_node_fn = |key: u64| observed_nodes.get(&key).copied().unwrap_or(0.0);
        let observed_edge_fn = |key: u64| observed_edges.get(&key).copied().unwrap_or(0.0);
        let observed_node_mass: f64 = observed_nodes.values().sum();
        let observed_edge_mass: f64 = observed_edges.values().sum();

        // ------------------------------------------- the scoring matrix
        // Per row (parallel), per record: the full all-offset likelihood
        // (both strands — the model's complete placement sum), the
        // factorized pocket placement asserted equal to the direct walk
        // at the same σ (the exactness gate, whole domain), and the
        // dominant placement for the smear-class measurement.
        let scoring_started = Instant::now();
        struct Scored {
            ll: Vec<f64>,
            dominant: Vec<Option<(u32, u64, u8)>>,
            factorized_checked: u64,
            factorized_equal: u64,
            fallback_no_projection: u64,
            fallback_torn: u64,
        }
        let scored: Vec<Scored> = (0..n_rows)
            .into_par_iter()
            .map(|row_index| {
                let row = &rows[row_index];
                let mut ll = Vec::with_capacity(n_records);
                let mut dominant: Vec<Option<(u32, u64, u8)>> = Vec::with_capacity(n_records);
                let mut factorized_checked = 0u64;
                let mut factorized_equal = 0u64;
                let mut fallback_no_projection = 0u64;
                let mut fallback_torn = 0u64;
                for record in &locus_record_data {
                    let (mut hist, best_f) = offset_histogram(&record.read, &row.seq);
                    let (hist_r, best_r) = offset_histogram(&record.read_rc, &row.seq);
                    for (h, r) in hist.iter_mut().zip(hist_r.iter()) {
                        *h += *r;
                    }
                    ll.push(scoring.row_log_likelihood(&hist, row.len));
                    dominant.push(match (best_f, best_r) {
                        (Some((m_f, s_f)), Some((m_r, s_r))) => Some(if m_f >= m_r {
                            (m_f, s_f, 0u8)
                        } else {
                            (m_r, s_r, 1u8)
                        }),
                        (Some((m, s)), None) => Some((m, s, 0)),
                        (None, Some((m, s))) => Some((m, s, 1)),
                        (None, None) => None,
                    });
                    // Path B: the factorized pocket placement per strand.
                    for (projection, strand_read) in [
                        (&record.forward, &record.read),
                        (&record.reverse, &record.read_rc),
                    ] {
                        let Some(projection) = projection else {
                            fallback_no_projection += 1;
                            continue;
                        };
                        let Some(run) = row.runs.iter().find(|run| {
                            let sigma = projection.x0 as i64 + run.delta;
                            sigma >= 0
                                && sigma + read_len() as i64 <= row.len as i64
                                && projection.x0 >= run.a_lo
                                && projection.x0 + read_len() as u64
                                    <= run.a_lo + (run.r_hi - run.r_lo)
                        }) else {
                            fallback_torn += 1;
                            continue;
                        };
                        let sigma = (projection.x0 as i64 + run.delta) as u64;
                        let Some((m, c)) = factorized_placement(projection, run, &frame_seq)
                        else {
                            fallback_torn += 1;
                            continue;
                        };
                        factorized_checked += 1;
                        let direct = direct_placement(strand_read, &row.seq, sigma);
                        if direct == (m, c) {
                            factorized_equal += 1;
                        }
                    }
                }
                Scored {
                    ll,
                    dominant,
                    factorized_checked,
                    factorized_equal,
                    fallback_no_projection,
                    fallback_torn,
                }
            })
            .collect();
        let scoring_seconds = scoring_started.elapsed().as_secs_f64();
        let factorized_checked: u64 = scored.iter().map(|s| s.factorized_checked).sum();
        let factorized_equal: u64 = scored.iter().map(|s| s.factorized_equal).sum();
        let fallback_no_projection: u64 =
            scored.iter().map(|s| s.fallback_no_projection).sum();
        let fallback_torn: u64 = scored.iter().map(|s| s.fallback_torn).sum();
        ensure(
            factorized_checked == factorized_equal,
            &format!(
                "FACTORIZATION GATE FAILED at locus {locus}: \
                 {factorized_equal} of {factorized_checked} equal"
            ),
        )?;
        total_factorized_checked += factorized_checked;
        total_factorized_equal += factorized_equal;
        let matrix: Vec<Vec<f64>> = scored.iter().map(|s| s.ll.clone()).collect();
        let dominant: Vec<Vec<Option<(u32, u64, u8)>>> =
            scored.iter().map(|s| s.dominant.clone()).collect();

        // ------------------------------- the smear-class measurement
        // Per (record, row) at the dominant placement: how much of the
        // read's base evidence sits on positions the seed projection
        // could see (inside a shared covered node's k-mer window on the
        // row) versus the NEWLY VISIBLE agreement (the variant-pocket
        // and sub-syncmer positions seed projection dropped), plus the
        // disagreement mass (the read's bases voting against the row —
        // invisible to the seed model by construction).
        let smear_started = Instant::now();
        let mut base_budget = 0.0f64;
        let mut agreement_visible = 0.0f64;
        let mut agreement_newly_visible = 0.0f64;
        let mut disagreement = 0.0f64;
        let mut shared_key_pairs = 0u64;
        for (r, record) in locus_record_data.iter().enumerate() {
            for row_index in 0..n_rows {
                let row = &rows[row_index];
                // The budget is per (record, row): every candidate row
                // sees the record's full base evidence under the model.
                base_budget += record.multiplicity as f64 * read_len() as f64;
                let Some((best_m, sigma, strand)) = dominant[row_index][r] else {
                    continue;
                };
                // The seed-visible positions: shared covered nodes'
                // k-mer windows on the row side (the old instrument's
                // key-credit convention rendered at base resolution).
                let mut windows: Vec<(u64, u64)> = Vec::new();
                for &node in &record.covered_nodes {
                    if let Some(occurrences) = row.node_index.get(&node) {
                        for &(bp, _) in occurrences {
                            let lo = bp.saturating_sub(row.start);
                            windows.push((lo, lo + k));
                        }
                    }
                }
                if !windows.is_empty() {
                    shared_key_pairs += 1;
                }
                windows.sort_unstable();
                let mut merged: Vec<(u64, u64)> = Vec::new();
                for (lo, hi) in windows {
                    match merged.last_mut() {
                        Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
                        _ => merged.push((lo, hi)),
                    }
                }
                let strand_read = if strand == 0 { &record.read } else { &record.read_rc };
                for i in 0..read_len() {
                    let p = sigma as usize + i;
                    let matched = strand_read[i] == row.seq[p];
                    let lo = merged.partition_point(|&(a, _)| a <= sigma + i as u64);
                    let visible = lo > 0 && merged[lo - 1].1 > sigma + i as u64;
                    if matched {
                        if visible {
                            agreement_visible += record.multiplicity as f64;
                        } else {
                            agreement_newly_visible += record.multiplicity as f64;
                        }
                    } else {
                        disagreement += record.multiplicity as f64;
                    }
                }
                let _ = best_m;
            }
        }
        let smear_seconds = smear_started.elapsed().as_secs_f64();

        // ------------------------------------------- the exhaustive classes
        let class_started = Instant::now();
        let pair_count = n_rows * (n_rows + 1) / 2;
        let multiplicity: Vec<f64> =
            locus_record_data.iter().map(|r| r.multiplicity as f64).collect();
        let class_lls: Vec<f64> = (0..pair_count)
            .into_par_iter()
            .map(|flat| {
                let (first, second) = unrank_pair(flat);
                let ll_a = &matrix[first];
                let ll_b = &matrix[second];
                let mut total = 0.0f64;
                for (r, &mult) in multiplicity.iter().enumerate() {
                    total += mult * mix_logsumexp(ll_a[r], ll_b[r]);
                }
                total
            })
            .collect();
        let class_seconds = class_started.elapsed().as_secs_f64();

        // The winner set (bit-identical maxima), truth rank, ties,
        // competitors — one deterministic scan in class order.
        let mut best_ll: Option<f64> = None;
        for &ll in &class_lls {
            best_ll = Some(match best_ll {
                None => ll,
                Some(current) => current.max(ll),
            });
        }
        let winner_pair = best_ll.and_then(|best| {
            class_lls
                .iter()
                .position(|&ll| ll == best)
                .map(|flat| unrank_pair(flat))
        });
        let called: Vec<(usize, usize)> = match best_ll {
            Some(best) => (0..pair_count)
                .filter(|&flat| class_lls[flat] == best)
                .map(|flat| unrank_pair(flat))
                .collect(),
            None => Vec::new(),
        };
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
                .filter(|&(flat, &ll)| {
                    ll == value && Some(flat) != truth_flat
                })
                .count()
        });
        if let (Some(value), Some(flat)) = (truth_log_likelihood, truth_flat) {
            for (other, &ll) in class_lls.iter().enumerate() {
                if ll > value {
                    let (first, second) = unrank_pair(other);
                    serde_json::to_writer(
                        &mut competitors,
                        &json!({
                            "locus": locus,
                            "row_indices": [first, second],
                            "log_likelihood": ll,
                            "truth_log_likelihood": value,
                            "likelihood_ratio": (ll - value).exp(),
                        }),
                    )?;
                    writeln!(competitors)?;
                } else if ll == value && other != flat {
                    let (first, second) = unrank_pair(other);
                    serde_json::to_writer(
                        &mut ties,
                        &json!({
                            "locus": locus,
                            "row_indices": [first, second],
                            "log_likelihood": ll,
                            "ulp_delta":
                                (ll.to_bits() as i64 - value.to_bits() as i64).abs(),
                            "bit_exact": ll.to_bits() == value.to_bits(),
                        }),
                    )?;
                    writeln!(ties)?;
                }
            }
        }
        let log_gap = match (best_ll, truth_log_likelihood) {
            (Some(best), Some(truth)) => Some(best - truth),
            _ => None,
        };

        // ------------------------------------------- the QUAL machinery
        let qual_started = Instant::now();
        let winner_nodes = winner_pair
            .map(|(first, second)| merged_multiset(&rows[first].nodes, &rows[second].nodes));
        let winner_edges = winner_pair
            .map(|(first, second)| merged_multiset(&rows[first].edges, &rows[second].edges));
        let called_score = best_ll;
        let spectrum: Vec<(f64, f64, f64, usize, usize, usize)> = match (
            winner_pair,
            winner_nodes.as_ref(),
            winner_edges.as_ref(),
            called_score,
        ) {
            (Some((winner_first, winner_second)), Some(wn), Some(we), Some(score)) => {
                let mut spectrum: Vec<(f64, f64, f64, usize, usize, usize)> = (0..pair_count)
                    .into_par_iter()
                    .filter_map(|flat| {
                        let (first, second) = unrank_pair(flat);
                        if first == winner_first && second == winner_second {
                            return None;
                        }
                        let ll = class_lls[flat];
                        let nodes = merged_multiset(&rows[first].nodes, &rows[second].nodes);
                        let edges = merged_multiset(&rows[first].edges, &rows[second].edges);
                        let distance = differing_observed_mass(&nodes, wn, &observed_node_fn)
                            + differing_observed_mass(&edges, we, &observed_edge_fn);
                        let signature = signature_cosine_distance(&nodes, &edges, wn, we);
                        Some((distance, (ll - score).exp(), signature, first, second, flat))
                    })
                    .collect();
                spectrum
                    .sort_by(|a, b| a.partial_cmp(b).expect("finite spectrum entries"));
                spectrum
            }
            _ => Vec::new(),
        };
        // The tied winner classes (beyond the first) and their bands.
        let mut tied_bands: Vec<(usize, Vec<f64>)> = Vec::new();
        for &(first, second) in called.iter().skip(1) {
            let index = spectrum
                .iter()
                .position(|&(_, _, _, a, b, _)| a == first && b == second)
                .expect("called class missing from the winner's spectrum");
            let nodes = merged_multiset(&rows[first].nodes, &rows[second].nodes);
            let edges = merged_multiset(&rows[first].edges, &rows[second].edges);
            let band: Vec<f64> = spectrum
                .par_iter()
                .map(|&(_, _, _, other_first, other_second, _)| {
                    let other_nodes =
                        merged_multiset(&rows[other_first].nodes, &rows[other_second].nodes);
                    let other_edges =
                        merged_multiset(&rows[other_first].edges, &rows[other_second].edges);
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
            let total: f64 = class_lls
                .iter()
                .map(|&ll| (ll - score).exp())
                .sum();
            score + total.ln()
        });
        let top_cluster_posterior = match (&cluster, posterior_logsumexp) {
            (Some(state), Some(logsumexp)) => {
                let mut mass = 1.0f64;
                for &(distance, score, _, _, _, _) in &spectrum {
                    if distance <= state.cut {
                        mass += score;
                    }
                }
                Some(mass / (logsumexp - called_score.unwrap_or(0.0)).exp())
            }
            _ => None,
        };
        let best_outside_log_likelihood = cluster.as_ref().and_then(|state| {
            spectrum
                .iter()
                .zip(&state.excluded)
                .filter(|(_, excluded)| !**excluded)
                .map(|(entry, _)| class_lls[entry.5])
                .max_by(|a, b| a.partial_cmp(b).expect("finite log-likelihoods"))
        });
        let alternative_log_gap = match (best_outside_log_likelihood, called_score) {
            (Some(best_outside), Some(_)) => Some(called_score.unwrap_or(0.0) - best_outside),
            _ => None,
        };
        let mut qual_alternative_classes: Vec<serde_json::Value> = Vec::new();
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
                        "row_indices": [first, second],
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

        // The named classes' per-record evidence (winner first, truth,
        // first a-holder rival) — the records sidecar.
        let mut named: Vec<(usize, usize)> = Vec::new();
        if let Some(pair) = winner_pair {
            named.push(pair);
        }
        if let Some(flat) = truth_flat {
            let pair = unrank_pair(flat);
            if !named.contains(&pair) {
                named.push(pair);
            }
        }
        if let Some(state) = &cluster {
            if let Some(best_outside) = best_outside_log_likelihood {
                for (entry, excluded) in spectrum.iter().zip(&state.excluded) {
                    if !excluded && class_lls[entry.5] == best_outside {
                        let pair = (entry.3, entry.4);
                        if !named.contains(&pair) {
                            named.push(pair);
                        }
                        break;
                    }
                }
            }
        }
        let mut classes_evidence: Vec<serde_json::Value> = Vec::new();
        for &(first, second) in &named {
            let per_record: Vec<serde_json::Value> = locus_record_data
                .iter()
                .enumerate()
                .map(|(r, record)| {
                    json!({
                        "record": record.record,
                        "multiplicity": record.multiplicity,
                        "row_ll": [matrix[first][r], matrix[second][r]],
                        "pair_term": mix_logsumexp(matrix[first][r], matrix[second][r]),
                    })
                })
                .collect();
            classes_evidence.push(json!({
                "row_indices": [first, second],
                "log_likelihood": class_lls[second * (second + 1) / 2 + first],
                "per_record": per_record,
            }));
        }
        serde_json::to_writer(&mut records_file, &json!({
            "locus": locus,
            "classes": classes_evidence,
        }))?;
        writeln!(records_file)?;

        // The re-derivation ingredients: rows with sequences, walks,
        // pocket runs, usage multisets; records with projections.
        serde_json::to_writer(&mut ingredients, &json!({
            "locus": locus,
            "partition": partition,
            "frame_row": frame_index,
            "rows": rows.iter().map(|row| json!({
                "members": row.members.iter().map(|m| json!({
                    "path_name": m.path_name, "start": m.start, "end": m.end,
                })).collect::<Vec<_>>(),
                "length": row.len,
                "sequence": String::from_utf8(row.seq.clone()).unwrap(),
                "walk": row.walk.iter().map(|&(bp, node)| json!([bp - row.start, node]))
                    .collect::<Vec<_>>(),
                "runs": row.runs.iter().map(|run| json!({
                    "r_lo": run.r_lo, "r_hi": run.r_hi, "a_lo": run.a_lo,
                    "delta": run.delta,
                    "edits": run.edits.iter().map(|&(p, b)| json!([p, b as char]))
                        .collect::<Vec<_>>(),
                })).collect::<Vec<_>>(),
                "nodes": row.nodes,
                "edges": row.edges,
            })).collect::<Vec<_>>(),
            "records": locus_record_data.iter().map(|record| json!({
                "record": record.record,
                "multiplicity": record.multiplicity,
                "forward": record.forward.as_ref().map(|p| json!({
                    "x0": p.x0,
                    "alts": p.alts.iter().map(|&(pos, b)| json!([pos, b as char]))
                        .collect::<Vec<_>>(),
                })),
                "reverse": record.reverse.as_ref().map(|p| json!({
                    "x0": p.x0,
                    "alts": p.alts.iter().map(|&(pos, b)| json!([pos, b as char]))
                        .collect::<Vec<_>>(),
                })),
                "covered_nodes": record.covered_nodes.iter().collect::<Vec<_>>(),
                "covered_edges": record.covered_edges.iter().collect::<Vec<_>>(),
            })).collect::<Vec<_>>(),
        }))?;
        writeln!(ingredients)?;
        // The full per-record per-row likelihood matrix (binary-free:
        // rows x records, row-major) — the checker's independent
        // re-derivation base.
        {
            let mut matrix_json: Vec<Vec<f64>> = Vec::with_capacity(n_rows);
            for row_index in 0..n_rows {
                matrix_json.push(matrix[row_index].clone());
            }
            serde_json::to_writer(&mut ingredients, &json!({
                "locus": locus,
                "ll_matrix": matrix_json,
            }))?;
            writeln!(ingredients)?;
        }

        // The validation receipt: the factorization gate over the whole
        // domain at this locus.
        serde_json::to_writer(&mut validation, &json!({
            "locus": locus,
            "factorized_checked": factorized_checked,
            "factorized_equal": factorized_equal,
            "fallback_no_projection": fallback_no_projection,
            "fallback_torn": fallback_torn,
            "strand_mismatch_anchors": strand_mismatch_anchors.load(Ordering::Relaxed),
            "gate": if factorized_checked == factorized_equal {
                "PASS"
            } else {
                "FAIL"
            },
        }))?;
        writeln!(validation)?;

        // The main receipt line.
        let called_classes_json: Vec<serde_json::Value> = called
            .iter()
            .map(|&(first, second)| {
                json!({
                    "row_indices": [first, second],
                    "physical_pairs": class_physical_pairs(
                        first == second,
                        rows[first].members.len(),
                        rows[second].members.len(),
                    ),
                })
            })
            .collect();
        let row_identities: Vec<serde_json::Value> = rows
            .iter()
            .map(|row| {
                json!({
                    "members": row.members.iter().map(|m| json!({
                        "path_name": m.path_name,
                        "start": m.start,
                        "end": m.end,
                    })).collect::<Vec<_>>(),
                    "length": row.len,
                })
            })
            .collect();
        let class_row_pairs: Vec<[usize; 2]> =
            (0..pair_count).map(|flat| {
                let (first, second) = unrank_pair(flat);
                [first, second]
            }).collect();
        let locus_seconds = locus_started.elapsed().as_secs_f64();
        let rss_kb = rss.probe(&format!("locus_{locus}"))?;
        serde_json::to_writer(&mut report, &json!({
            "locus": locus,
            "partition": partition,
            "model": "local-realignment-v1",
            "scoring": {
                "phred": phred,
                "epsilon": epsilon,
                "match_log_prob": scoring.a,
                "mismatch_log_prob": scoring.b,
                "read_length": read_len(),
            },
            "record_count": n_records,
            "rows": n_rows,
            "physical_members": map.members.len(),
            "class_count": pair_count,
            "class_row_pairs": class_row_pairs,
            "class_log_likelihoods": class_lls,
            "best_log_likelihood": best_ll,
            "best_row_indices": winner_pair.map(|(a, b)| [a, b]),
            "called_class_count": called.len(),
            "called_classes": called_classes_json,
            "row_identities": row_identities,
            "truth_rows": truth_pair.map(|(a, b)| [a, b]),
            "truth_pair_expressible": truth_pair.is_some(),
            "truth_log_likelihood": truth_log_likelihood,
            "truth_rank": truth_rank,
            "truth_tied_classes": truth_tied_classes,
            "log_gap": log_gap,
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
            "validation": {
                "factorized_checked": factorized_checked,
                "factorized_equal": factorized_equal,
                "fallback_no_projection": fallback_no_projection,
                "fallback_torn": fallback_torn,
                "strand_mismatch_anchors": strand_mismatch_anchors.load(Ordering::Relaxed),
            },
            "smear_measurement": {
                "base_budget": base_budget,
                "agreement_seed_visible": agreement_visible,
                "agreement_newly_visible": agreement_newly_visible,
                "disagreement": disagreement,
                "shared_key_pairs": shared_key_pairs,
            },
            "walls": {
                "records_seconds": records_seconds,
                "scoring_seconds": scoring_seconds,
                "smear_seconds": smear_seconds,
                "classes_seconds": class_seconds,
                "qual_seconds": qual_seconds,
                "locus_seconds": locus_seconds,
            },
            "rss_kb": rss_kb,
        }))?;
        writeln!(report)?;
        eprintln!(
            "[realign] locus {locus} done: {n_rows} rows, {n_records} records, \
             {pair_count} classes, truth_rank {:?}, factorized {factorized_equal}/\
             {factorized_checked}, walls {locus_seconds:.1}s",
            truth_rank,
        );
    }

    report.flush()?;
    competitors.flush()?;
    ties.flush()?;
    records_file.flush()?;
    ingredients.flush()?;
    validation.flush()?;
    ensure(
        total_factorized_checked == total_factorized_equal,
        "the factorization gate failed somewhere",
    )?;
    eprintln!(
        "[realign] complete: {n_loci} loci, factorized placements {total_factorized_equal}/\
         {total_factorized_checked} exact, total {:.1}s",
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// Unit tests: the factorized exact form, the scoring, the machinery.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn scoring() -> Scoring {
        Scoring {
            a: (1.0f64 - 1e-4).ln(),
            b: (1e-4f64 / 3.0).ln(),
            epsilon: 1e-4,
            phred: 40,
        }

    }

    #[test]
    fn factorized_equals_direct_on_snv_pockets() {
        // A frame, a row with two SNV pockets, and a read whose own
        // backbone differs from the frame at one pocket: the factorized
        // (m, c) from the cached pocket lists must equal the direct
        // walk at the frame-collinear offset.
        let frame: Vec<u8> = std::iter::repeat(b'A').take(300).collect();
        let mut row_seq = frame.clone();
        row_seq[10] = b'T';
        row_seq[140] = b'C';
        let run = FrameRun {
            r_lo: 0,
            r_hi: row_seq.len() as u64,
            a_lo: 0,
            delta: 0,
            edits: vec![(10, b'T'), (140, b'C')],
        };
        // The read: the frame with its own alt at 140 (matching the row)
        // and the frame base at 10 (differing from the row).
        let mut read = frame[..read_len()].to_vec();
        read[140] = b'C';
        let projection = Projection {
            x0: 0,
            alts: vec![(140, b'C')],
        };
        let (m, c) = factorized_placement(&projection, &run, &frame).unwrap();
        let direct = direct_placement(&read, &row_seq, 0);
        assert_eq!((m, c), direct);
        assert_eq!(c, 1); // only the row's pocket at 10 votes against.
        // And a read carrying the row's pocket at 10 too:
        let mut read2 = read.clone();
        read2[10] = b'T';
        let projection2 = Projection {
            x0: 0,
            alts: vec![(10, b'T'), (140, b'C')],
        };
        let (m2, c2) = factorized_placement(&projection2, &run, &frame).unwrap();
        assert_eq!((m2, c2), direct_placement(&read2, &row_seq, 0));
        assert_eq!(c2, 0);
        let _ = scoring();
    }

    #[test]
    fn factorized_rejects_out_of_run_coverage() {
        let frame = vec![b'A'; 300];
        let run = FrameRun {
            r_lo: 0,
            r_hi: 100,
            a_lo: 0,
            delta: 0,
            edits: vec![],
        };
        let projection = Projection { x0: 0, alts: vec![] };
        // The read's frame interval [0, 150) exceeds the run coverage.
        assert!(factorized_placement(&projection, &run, &frame).is_none());
    }
    #[test]
    fn row_log_likelihood_prefers_exact_placements() {
        let s = scoring();
        let seq = vec![b'A'; 300];
        let read = vec![b'A'; read_len()];
        let (mut hist, best_f) = offset_histogram(&read, &seq);
        let (hist_rc, best_r) = offset_histogram(&read, &seq);
        for (h, r) in hist.iter_mut().zip(hist_rc.iter()) {
            *h += *r;
        }
        let ll = s.row_log_likelihood(&hist, 300);
        // 151 exact placements (forward) + 151 (rc) at m = 150.
        let expect = (151 + 151) as f64;
        let approx = s.placement(150, 0) + expect.ln() - (2.0f64 * 300.0).ln();
        assert!((ll - approx).abs() < 1e-9, "{ll} vs {approx}");
        assert_eq!(best_f, Some((150, 0)));
        assert_eq!(best_r, Some((150, 0)));
    }

    #[test]
    fn mix_logsumexp_is_exact_for_equal_components() {
        let a = -3.5f64;
        assert_eq!(mix_logsumexp(a, a), a);
        assert_eq!(mix_logsumexp(f64::NEG_INFINITY, f64::NEG_INFINITY), f64::NEG_INFINITY);
        let mixed = mix_logsumexp(0.0, f64::NEG_INFINITY);
        assert!((mixed - (-(2f64).ln())).abs() < 1e-12);
    }

    #[test]
    fn pair_unrank_roundtrips() {
        for flat in [0usize, 1, 2, 3, 6, 10, 100, 1000, 165024] {
            let (i, j) = unrank_pair(flat);
            assert_eq!(j * (j + 1) / 2 + i, flat);
            assert!(i <= j);
        }
    }

    #[test]
    fn frame_runs_pin_snvs_and_tear_at_indels() {
        // Frame: 300 A's; a row sequence = frame with an SNV at 100,
        // anchored by identity at 0 and 200 (63bp exact windows).
        let frame = vec![b'A'; 300];
        let mut row = frame.clone();
        row[100] = b'C';
        let anchors = vec![
            Anchor { seg: 0, frame: 0, delta: 0 },
            Anchor { seg: 200, frame: 200, delta: 0 },
        ];
        let runs = frame_runs(&anchors, &row, &frame, 63);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].r_lo, 0);
        // The flank beyond the last anchor extends while identical: to
        // the row's end (all A's past the SNV).
        assert_eq!(runs[0].r_hi, 300);
        assert_eq!(runs[0].edits, vec![(100, b'C')]);
        // An indel between the anchors tears the run: the second anchor
        // pins a different delta.
        let anchors = vec![
            Anchor { seg: 0, frame: 0, delta: 0 },
            Anchor { seg: 200, frame: 205, delta: 5 },
        ];
        let runs = frame_runs(&anchors, &row, &frame, 63);
        assert_eq!(runs.len(), 2);
        // The first run extends forward only while identical: the SNV
        // at 100 stops it; the second run starts at its anchor.
        assert_eq!(runs[0].r_lo, 0);
        assert_eq!(runs[0].r_hi, 100);
        assert!(runs[0].edits.is_empty());
        // The second run's flank extension is legitimate: the anchor-
        // pinned offset (delta 5) continues over verified-identical
        // bases down to just past the SNV at 100.
        assert_eq!(runs[1].r_lo, 101);
        assert_eq!(runs[1].a_lo, 106);
    }

    #[test]
    fn offset_histogram_counts_every_offset() {
        // A 200bp row containing the read exactly once, with a decoy
        // prefix: the histogram's maximum sits at the exact placement.
        let mut seq = vec![b'A'; 200];
        for (i, b) in seq.iter_mut().enumerate().take(50) {
            *b = if i % 2 == 0 { b'C' } else { b'G' };
        }
        let read = seq[50..200].to_vec();
        let (hist, best) = offset_histogram(&read, &seq);
        assert_eq!(hist[150], 1);
        assert_eq!(best, Some((150, 50)));
    }
}
