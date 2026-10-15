//! THE CHAIN/PHASING LAYER — the thin DP that links the CLI's
//! per-locus diplotype calls into TWO WHOLE-CHROMOSOME MOLECULES per
//! component (the owner's product ruling from the campaign's start:
//! the diploid product is two molecules; per-locus genotyping is the
//! product, the chain is a thin phasing layer over it).
//!
//! THE MODEL, derived and stated (no thresholds, no tuning constants):
//!
//! 1. INPUT: the CLOSED per-locus calls (`call-diplotypes`'s
//!    calls.jsonl — the two called paths per locus, the class
//!    identity, the model-internal QUAL, E, the alternative log-gap).
//!    The chain NEVER writes the calls; it only reads them, so the
//!    per-locus product is bit-identical pre/post chain by
//!    construction (the gate proves it by fingerprint). The
//!    boundary evidence comes from two truth-free channels:
//!    (a) THE ADJACENCY LAYER — the called folds' own member rows:
//!    shared panel paths whose rows continue across the window
//!    junction (row_f.end == row_g.start, either direction) or
//!    overlap in shared material;
//!    (b) THE EDGE LAYER — the junction-spanning reads' votes (the
//!    layer COSIGT never had): census occurrences touching BOTH
//!    consecutive windows are verified placements of a read pattern
//!    spanning the junction; each is attributed to the called folds
//!    its covered panel material (its contained segment-node
//!    windows on its own path) overlaps, and votes on which homolog
//!    of one locus connects to which homolog of the next.
//!
//! 2. THE CHAIN: a two-state orientation DP. The called pair per
//!    locus is UNORDERED; a chain state s_j in {0,1} says which
//!    emitted fold of locus j is carried on molecule 1. Each boundary
//!    between consecutive called loci carries two hypotheses —
//!    "same" (s_{j+1} = s_j) and "swap" (s_{j+1} = 1 - s_j) — with
//!    measured evidence for both. The DECISION RULE is the stated
//!    lexicographic order over the measured channels:
//!    (i) the edge votes (record-once: a census record's read
//!    instances are its mass, split uniformly over its crossing
//!    occurrences at the boundary and over the attributed fold
//!    pairs; both-side attribution is required — one-sided material
//!    is ORPHAN crossing mass, counted and reported, never a vote,
//!    because the calls do not express that molecule's continuity
//!    there);
//!    (ii) if the votes tie exactly, the adjacency continuation
//!    (exact junction touches dominate, then shared-overlap bp);
//!    (iii) if both tie, the boundary is UNIDENTIFIABLE and the
//!    stated deterministic rule carries the chain through ("same"),
//!    with the tie flagged — the chain may not invent certainty.
//!    A tie is exact equality of the accumulated f64 masses:
//!    structurally symmetric boundaries (homozygous or identical
//!    called pairs) produce bitwise-identical accumulations by
//!    construction, and every boundary's raw margin is reported to
//!    full precision so a near-tie is visible as what it is.
//!    With only relative orientation evidence the 2-state Viterbi
//!    maximizes each boundary independently; the accumulation
//!    s_{j+1} = s_j XOR (boundary = swap) builds the two molecules.
//!    No locus call is changed: the molecules carry exactly the
//!    emitted pairs, reoriented only.
//!
//! 3. AMBIGUOUS LOCI, carried honestly: a locus whose called pair
//!    carries tied classes (called_class_count > 1) enters the
//!    linkage through its EMITTED winner pair; every boundary at it
//!    is flagged through_ambiguous_locus (an alternative tied pair
//!    could reorient it). A homozygous emitted pair (equal fold
//!    indices) makes both adjacent boundaries structurally
//!    unidentifiable (the evidence is symmetric by construction).
//!    The link confidence reports each endpoint's model-internal
//!    QUAL — the per-locus call certainty weighting the uncertain
//!    links — and never a synthesized posterior.
//!
//! 4. EMISSION: molecules.jsonl, one record per component — the two
//!    molecules (the per-locus orientation bits, the per-locus fold
//!    material with junction accounting, the spelled sequences), and
//!    every boundary's evidence, orientation, decision channel, and
//!    link confidence.
//!
//! 5. THE TEST MODE (assessment-side ONLY, the owner's ruling): the
//!    --truth-qv-file input consumes the product run's OWN
//!    calls.jsonl.truth-qv.jsonl (never the default output). Per
//!    boundary, the truth orientation comes from the committed
//!    assignment yardstick (the emitted assignment; recomputed
//!    rate-and-edit ties bracket the boundary — the truth side's own
//!    unidentifiability, never silently resolved). The SWITCH count
//!    (chain orientation != truth orientation over assessable
//!    boundaries) and the PER-LOCUS accuracy under best-case
//!    assignment (carried UNCHANGED from the artifact into separate
//!    columns) are reported side by side, NEVER conflated — the
//!    owner's evaluation ruling since the haploid era.

use crate::genome_inference::{panel_routes as routes, read_json, PanelIdentity};
use crate::syng::{SyncmerParams, SyngIndex};
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

pub const MODEL: &str = "diplotype-chain-orientation-v1";

/// The sample's uniform read length (asserted L150 by the product
/// run's own derive-cache load; a census occurrence's covered hull is
/// bounded by it — the walk span never exceeds the read).
const READ_LENGTH: u64 = 150;

pub struct Options {
    /// The CLI product's closed per-locus calls (calls.jsonl)
    pub calls: PathBuf,
    /// The panel syng prefix (path name/id resolution, node walks)
    pub panel: String,
    /// The route graph directory (graph.json + sources; sequence spelling)
    pub routes: PathBuf,
    /// The locality map directory (window<N>.gfa.map.json; axis extents)
    pub partition_graphs: PathBuf,
    /// The multi-matching census receipt (the same truth-free input
    /// the product run consumed; the crossing records)
    pub census: PathBuf,
    /// THE TEST MODE (assessment-side only): the product run's
    /// calls.jsonl.truth-qv.jsonl — the truth-referenced per-locus
    /// artifact, consumed read-only into molecules.jsonl.truth-qv.jsonl
    pub truth_qv: Option<PathBuf>,
    /// Resident-set guard in GiB (the 64 GiB discipline; 0 = no guard)
    pub rss_budget_gib: f64,
    /// New output directory (refused if it exists)
    pub out_dir: PathBuf,
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
        eprintln!("[chain] rss {stage}: {rss_kb} kB");
        Ok(rss_kb)
    }
}

// ---------------------------------------------------------------------------
// The closed calls (the product record, the fields the chain consumes).
// ---------------------------------------------------------------------------

#[derive(Clone, Deserialize)]
struct MemberRow {
    path_name: String,
    start: u64,
    end: u64,
}

#[derive(Deserialize)]
struct CallFold {
    #[allow(dead_code)]
    fold_index: u64,
    length: u64,
    members: Vec<MemberRow>,
    #[serde(default)]
    territory_image: Option<Vec<u64>>,
}

#[derive(Deserialize)]
struct CallDiplotype {
    fold_indices: Vec<u64>,
    paths: Vec<CallFold>,
}

#[derive(Deserialize)]
struct CallRecord {
    locus: u32,
    partition: u32,
    component: String,
    diplotype: CallDiplotype,
    called_class_count: u64,
    #[serde(default)]
    qual: Option<f64>,
}

/// One called locus, chain-ready: the unordered pair of folds with
/// their member rows keyed by panel path id, plus the per-locus
/// honesty fields (the emitted pair is NEVER altered).
struct LocusCall {
    locus: u32,
    partition: u32,
    fold_indices: [u64; 2],
    /// Per emitted fold: the member rows (path id, start, end capped
    /// at the path length), grouped by path id for attribution.
    fold_rows: [Vec<(usize, u64, u64)>; 2],
    /// The fold's spelled-row representative (the first member; the
    /// fold's rows are identical-sequence coalesced copies by the
    /// product's construction), per emitted fold.
    fold_row_head: [(String, u64, u64); 2],
    fold_length: [u64; 2],
    homozygous: bool,
    tied_classes: bool,
    qual: Option<f64>,
}

// ---------------------------------------------------------------------------
// The census (the committed multi-matching census structures,
// mirrored from the scoring layer verbatim).
// ---------------------------------------------------------------------------

/// One verified occurrence in the committed multi-matching census.
#[derive(Deserialize)]
struct CensusOccurrence {
    path: usize,
    start: u64,
    /// The touched component windows (the routing's own
    /// territory-touch rule; WINDOW ids, not partition ids).
    #[serde(default)]
    partitions: Vec<u32>,
    /// Per read-matched interval, the contained covered node abs ids
    /// in path order (the committed containment rule) — the
    /// occurrence's exact covered panel material.
    #[serde(default)]
    intervals: Vec<Vec<u32>>,
}

#[derive(Deserialize)]
struct CensusRecordLine {
    #[allow(dead_code)]
    record: u64,
    multiplicity: u64,
    #[serde(default)]
    occurrences: Vec<CensusOccurrence>,
}

// ---------------------------------------------------------------------------
// The pure chain core (unit-tested).
// ---------------------------------------------------------------------------

/// The adjacency continuation between two folds' member rows:
/// (exact junction touches, shared-overlap bp) over every pair of
/// rows on the SAME panel path.
fn fold_pair_adjacency(a: &[(usize, u64, u64)], b: &[(usize, u64, u64)]) -> (u64, u64) {
    let mut touches = 0u64;
    let mut overlap_bp = 0u64;
    for &(path_a, start_a, end_a) in a {
        for &(path_b, start_b, end_b) in b {
            if path_a != path_b {
                continue;
            }
            if end_a == start_b || end_b == start_a {
                touches += 1;
            }
            let lo = start_a.max(start_b);
            let hi = end_a.min(end_b);
            if hi > lo {
                overlap_bp += hi - lo;
            }
        }
    }
    (touches, overlap_bp)
}

/// Does any row contain (overlap) any material window?
fn rows_touch_material(rows: &[(u64, u64)], windows: &[(u64, u64)]) -> bool {
    for &(row_start, row_end) in rows {
        for &(window_lo, window_hi) in windows {
            if row_start < window_hi && row_end > window_lo {
                return true;
            }
        }
    }
    false
}

/// THE BOUNDARY DECISION (the stated lexicographic rule): returns
/// (swap?, decided_by, unidentifiable). Votes first (the edge layer);
/// exact vote ties fall to the adjacency continuation (touches
/// dominate, then shared-overlap bp); both tied = unidentifiable,
/// carried through as "same" (the stated deterministic rule).
fn decide_boundary(
    votes_same: f64,
    votes_swap: f64,
    adjacency_same: (u64, u64),
    adjacency_swap: (u64, u64),
) -> (bool, &'static str, bool) {
    if votes_same > votes_swap {
        (false, "crossing", false)
    } else if votes_swap > votes_same {
        (true, "crossing", false)
    } else if adjacency_same > adjacency_swap {
        (false, "adjacency", false)
    } else if adjacency_swap > adjacency_same {
        (true, "adjacency", false)
    } else {
        (false, "tie", true)
    }
}

/// The orientation accumulation: state 0 = the first emitted fold on
/// molecule 1; each swap boundary flips the state.
fn accumulate_orientations(swaps: &[bool]) -> Vec<u8> {
    let mut states = Vec::with_capacity(swaps.len() + 1);
    let mut state = 0u8;
    states.push(state);
    for &swap in swaps {
        state ^= swap as u8;
        states.push(state);
    }
    states
}

// ---------------------------------------------------------------------------
// The test mode (assessment-side only): the truth-QV artifact of the
// product run, consumed read-only.
// ---------------------------------------------------------------------------

#[derive(Clone, Default, Deserialize)]
struct TruthPairCounts {
    #[serde(default)]
    edits: u64,
    #[serde(default)]
    columns: u64,
}

#[derive(Clone, Default, Deserialize)]
struct TruthPairs {
    #[serde(default)]
    called0_truth0: TruthPairCounts,
    #[serde(default)]
    called1_truth1: TruthPairCounts,
    #[serde(default)]
    called0_truth1: TruthPairCounts,
    #[serde(default)]
    called1_truth0: TruthPairCounts,
}

#[derive(Deserialize)]
struct TruthQvRecord {
    locus: u32,
    truth_pair_expressible: bool,
    #[serde(default)]
    assignment: Option<String>,
    #[serde(default)]
    pairs: Option<TruthPairs>,
    #[serde(default)]
    rank1: Option<bool>,
    #[serde(default)]
    perfect: Option<bool>,
    #[serde(default)]
    truth_rank: Option<u64>,
    #[serde(default)]
    error: Option<f64>,
    #[serde(default)]
    qv: Option<f64>,
}

/// The truth-assignment state of one locus for the switch table:
/// either the emitted assignment's orientation sign (true = the
/// emitted fold 0 sits on truth haplotype 1 — a "crossed"
/// assignment), or the bracket reason (the truth side's own
/// unidentifiability, never silently resolved).
enum TruthAssignment {
    Swapped(bool),
    Bracket(&'static str),
}

fn truth_assignment(record: &TruthQvRecord) -> TruthAssignment {
    if !record.truth_pair_expressible {
        return TruthAssignment::Bracket("truth_pair_not_expressible");
    }
    let Some(assignment) = record.assignment.as_deref() else {
        return TruthAssignment::Bracket("truth_assignment_absent");
    };
    // The committed yardstick's own tie, recomputed exactly from the
    // emitted pair counts: equal rates AND equal total edits — the
    // two orientations indistinguishable at this locus (a homozygous
    // called pair, or identical assignments both ways).
    if let Some(pairs) = &record.pairs {
        let direct_edits = pairs.called0_truth0.edits + pairs.called1_truth1.edits;
        let direct_columns = pairs.called0_truth0.columns + pairs.called1_truth1.columns;
        let crossed_edits = pairs.called0_truth1.edits + pairs.called1_truth0.edits;
        let crossed_columns = pairs.called0_truth1.columns + pairs.called1_truth0.columns;
        let direct_rate = if direct_columns > 0 {
            direct_edits as f64 / direct_columns as f64
        } else {
            0.0
        };
        let crossed_rate = if crossed_columns > 0 {
            crossed_edits as f64 / crossed_columns as f64
        } else {
            0.0
        };
        if direct_rate == crossed_rate && direct_edits == crossed_edits {
            return TruthAssignment::Bracket("truth_assignment_tied");
        }
    }
    match assignment {
        "crossed" => TruthAssignment::Swapped(true),
        _ => TruthAssignment::Swapped(false),
    }
}

// ---------------------------------------------------------------------------
// The chain run.
// ---------------------------------------------------------------------------

pub fn run(options: Options) -> io::Result<()> {
    let started = Instant::now();
    let rss = RssGuard::new(options.rss_budget_gib);

    // ------------------------------------------------------- the closed calls
    let calls_started = Instant::now();
    let mut raw_calls: Vec<CallRecord> = Vec::new();
    for line in BufReader::new(File::open(&options.calls)?).lines() {
        raw_calls.push(serde_json::from_str(&line?)?);
    }
    ensure(!raw_calls.is_empty(), "the calls file is empty")?;
    let component = raw_calls[0].component.clone();
    ensure(
        raw_calls.iter().all(|record| record.component == component),
        "the calls file mixes components",
    )?;
    raw_calls.sort_by_key(|record| record.locus);
    ensure(
        raw_calls.windows(2).all(|pair| pair[0].locus < pair[1].locus),
        "the calls file repeats a locus",
    )?;

    // The panel: path name/id resolution and the node walks; the route
    // sources: the sequence spelling. (Both are the product run's own
    // inputs, unchanged.)
    let identity = PanelIdentity::read(&options.panel)?;
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    let k = panel.syncmer_length_bp() as u64;
    let path_of_name: HashMap<String, usize> = panel
        .name_map
        .name_to_path
        .iter()
        .map(|(name, &path)| (name.clone(), path as usize))
        .collect();
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
    let source_of_name: HashMap<&str, usize> = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.as_str(), lane.id))
        .collect();
    let lane_lengths: Vec<u64> = sources.lanes.iter().map(|&(_, len)| len).collect();
    let fetch_seq = |name: &str, start: u64, end: u64| -> io::Result<Vec<u8>> {
        let path_idx = *path_of_name
            .get(name)
            .ok_or_else(|| invalid("panel path absent from the syng name map"))?;
        let end = end.min(panel.name_map.path_to_length[path_idx]);
        let source = *source_of_name
            .get(name)
            .ok_or_else(|| invalid("panel path absent from the route sources"))?;
        if start >= end || source >= lane_lengths.len() {
            return Ok(Vec::new());
        }
        sources.fetch(source, start, end)
    };
    rss.probe("inputs")?;

    // The window axis extents (the committed derivation, mirrored:
    // the component's own member rows ranked by start map to the
    // window ids in order — slice A's derivation, one window per axis
    // row).
    let mut axis_rows: Vec<(u32, u64, u64)> = Vec::new();
    for entry in std::fs::read_dir(&options.partition_graphs)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if !name.ends_with(".gfa.map.json") {
            continue;
        }
        let text = std::fs::read_to_string(&path)?;
        #[derive(Deserialize)]
        struct PartitionMap {
            #[allow(dead_code)]
            partition: u32,
            members: Vec<MemberRow>,
        }
        let map: PartitionMap = serde_json::from_str(&text)?;
        for axis in map
            .members
            .iter()
            .filter(|m| m.path_name == component)
        {
            axis_rows.push((map.partition, axis.start, axis.end));
        }
    }
    axis_rows.sort_by_key(|&(partition, start, end)| (start, end, partition));
    let window_axis: Vec<(u64, u64)> = axis_rows
        .iter()
        .map(|&(_, start, end)| (start, end))
        .collect();
    ensure(
        raw_calls.iter().all(|record| (record.locus as usize) < window_axis.len()),
        "a called locus is outside the axis window span",
    )?;

    // Chain-ready loci: the emitted pair per locus, member rows keyed
    // by panel path id (ends capped at the path length, the scoring
    // layer's own convention).
    let mut loci: Vec<LocusCall> = Vec::with_capacity(raw_calls.len());
    for record in &raw_calls {
        ensure(
            record.diplotype.paths.len() == 2
                && record.diplotype.fold_indices.len() == 2
                && record.diplotype.paths[0].members.len() > 0
                && record.diplotype.paths[1].members.len() > 0,
            "a called locus lacks its two-path diplotype",
        )?;
        let mut fold_rows: [Vec<(usize, u64, u64)>; 2] = [Vec::new(), Vec::new()];
        let mut fold_row_head: [(String, u64, u64); 2] = Default::default();
        let mut fold_length = [0u64; 2];
        for side in 0..2 {
            let fold = &record.diplotype.paths[side];
            let head = &fold.members[0];
            let path_idx = *path_of_name
                .get(&head.path_name)
                .ok_or_else(|| invalid("a fold row's path is absent from the syng"))?;
            let capped_end = head.end.min(panel.name_map.path_to_length[path_idx]);
            fold_row_head[side] = (head.path_name.clone(), head.start, capped_end);
            fold_length[side] = capped_end.saturating_sub(head.start);
            for member in &fold.members {
                let path_idx = *path_of_name
                    .get(&member.path_name)
                    .ok_or_else(|| invalid("a fold row's path is absent from the syng"))?;
                let end = member.end.min(panel.name_map.path_to_length[path_idx]);
                if member.start < end {
                    fold_rows[side].push((path_idx, member.start, end));
                }
            }
        }
        loci.push(LocusCall {
            locus: record.locus,
            partition: record.partition,
            fold_indices: [
                record.diplotype.fold_indices[0],
                record.diplotype.fold_indices[1],
            ],
            fold_rows,
            fold_row_head,
            fold_length,
            homozygous: record.diplotype.fold_indices[0] == record.diplotype.fold_indices[1],
            tied_classes: record.called_class_count > 1,
            qual: record.qual,
        });
    }
    eprintln!(
        "[chain] calls loaded: {} loci [{:.1}s]",
        loci.len(),
        calls_started.elapsed().as_secs_f64(),
    );
    rss.probe("calls")?;

    // ------------------------------------------------- the boundary evidence
    // Consecutive called loci define the boundaries; the census's
    // window-touch field identifies the crossing occurrences.
    let n_boundaries = loci.len().saturating_sub(1);
    // window id -> its index among the called loci; the boundary
    // between consecutive called loci j and j+1 is the index j.
    let locus_index: HashMap<u32, usize> = loci
        .iter()
        .enumerate()
        .map(|(index, locus)| (locus.locus, index))
        .collect();
    let successor: HashMap<u32, u32> = (0..n_boundaries)
        .map(|index| (loci[index].locus, loci[index + 1].locus))
        .collect();

    struct BoundaryEvidence {
        /// votes[fold of locus j][fold of locus j+1]: the record-once
        /// crossing mass attributed to that link.
        votes: [[f64; 2]; 2],
        crossing_records: u64,
        crossing_occurrences: u64,
        /// One-sided crossing occurrences (material escapes the
        /// called folds on one side): the molecule's continuity is
        /// not expressed by the calls there — counted, never voted.
        orphan_occurrences: u64,
        /// Crossing occurrences whose material escapes the called
        /// folds on BOTH sides.
        unattributed_occurrences: u64,
    }
    let mut evidence: Vec<BoundaryEvidence> = (0..n_boundaries)
        .map(|_| BoundaryEvidence {
            votes: [[0.0; 2]; 2],
            crossing_records: 0,
            crossing_occurrences: 0,
            orphan_occurrences: 0,
            unattributed_occurrences: 0,
        })
        .collect();

    let census_started = Instant::now();
    let (reader, _) = niffler::get_reader(Box::new(File::open(&options.census)?))
        .map_err(io::Error::other)?;
    let mut census_records = 0u64;
    for line in BufReader::new(reader).lines() {
        let line = line?;
        let record: CensusRecordLine = serde_json::from_str(&line)?;
        census_records += 1;
        // The record's crossing occurrences, grouped by boundary
        // (record-once: the record's mass splits over its crossing
        // occurrences AT ONE BOUNDARY).
        let mut by_boundary: BTreeMap<usize, Vec<&CensusOccurrence>> = BTreeMap::new();
        for occurrence in &record.occurrences {
            for &window in &occurrence.partitions {
                if let Some(&next) = successor.get(&window) {
                    if occurrence.partitions.contains(&next) {
                        let boundary = locus_index[&window];
                        by_boundary.entry(boundary).or_default().push(occurrence);
                    }
                }
            }
        }
        for (boundary, occurrences) in by_boundary {
            let slot = &mut evidence[boundary];
            slot.crossing_records += 1;
            slot.crossing_occurrences += occurrences.len() as u64;
            let share = record.multiplicity as f64 / occurrences.len() as f64;
            for occurrence in occurrences {
                // The occurrence's covered panel material: its
                // contained node windows on its own path (the walk
                // over the read-length hull; the census's own
                // containment rule guarantees the windows lie inside).
                let mut material: Vec<(u64, u64)> = Vec::new();
                if !occurrence.intervals.is_empty() {
                    let nodes: HashSet<u32> =
                        occurrence.intervals.iter().flatten().copied().collect();
                    if !nodes.is_empty() {
                        for (signed_node, bp) in panel.walk_path_range(
                            occurrence.path,
                            occurrence.start,
                            occurrence.start + READ_LENGTH,
                        )? {
                            if nodes.contains(&signed_node.unsigned_abs()) {
                                material.push((bp, bp + k));
                            }
                        }
                        material.sort_unstable();
                        material.dedup();
                    }
                }
                if material.is_empty() {
                    slot.unattributed_occurrences += 1;
                    continue;
                }
                // Fold attribution on each side: the called folds whose
                // member rows on the occurrence's own path contain
                // any material window.
                let attributed = |locus: &LocusCall| -> usize {
                    let mut mask = 0usize;
                    for side in 0..2 {
                        let rows: Vec<(u64, u64)> = locus.fold_rows[side]
                            .iter()
                            .filter(|&&(path, _, _)| path == occurrence.path)
                            .map(|&(_, start, end)| (start, end))
                            .collect();
                        if rows_touch_material(&rows, &material) {
                            mask |= 1 << side;
                        }
                    }
                    mask
                };
                let mask_left = attributed(&loci[boundary]);
                let mask_right = attributed(&loci[boundary + 1]);
                match (mask_left, mask_right) {
                    (0, 0) => slot.unattributed_occurrences += 1,
                    (0, _) | (_, 0) => slot.orphan_occurrences += 1,
                    _ => {
                        let sides_left = mask_left.count_ones() as f64;
                        let sides_right = mask_right.count_ones() as f64;
                        for left in 0..2 {
                            for right in 0..2 {
                                if mask_left & (1 << left) != 0 && mask_right & (1 << right) != 0 {
                                    slot.votes[left][right] +=
                                        share / (sides_left * sides_right);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    eprintln!(
        "[chain] census streamed: {census_records} records [{:.1}s]",
        census_started.elapsed().as_secs_f64(),
    );
    rss.probe("census")?;

    // ------------------------------------------ the decisions and the chain
    // Per boundary: both evidence channels, the decision, and the
    // accumulated orientation bits.
    struct BoundaryDecision {
        locus_left: u32,
        locus_right: u32,
        window_gap_bp: i64,
        adjacency_same: (u64, u64),
        adjacency_swap: (u64, u64),
        votes_same: f64,
        votes_swap: f64,
        swap: bool,
        decided_by: &'static str,
        unidentifiable: bool,
    }
    let mut decisions: Vec<BoundaryDecision> = Vec::with_capacity(n_boundaries);
    for index in 0..n_boundaries {
        let (left, right) = (&loci[index], &loci[index + 1]);
        let adjacency = |left_side: usize, right_side: usize| {
            fold_pair_adjacency(&left.fold_rows[left_side], &right.fold_rows[right_side])
        };
        let adjacency_same = (
            adjacency(0, 0).0 + adjacency(1, 1).0,
            adjacency(0, 0).1 + adjacency(1, 1).1,
        );
        let adjacency_swap = (
            adjacency(0, 1).0 + adjacency(1, 0).0,
            adjacency(0, 1).1 + adjacency(1, 0).1,
        );
        let slot = &evidence[index];
        let votes_same = slot.votes[0][0] + slot.votes[1][1];
        let votes_swap = slot.votes[0][1] + slot.votes[1][0];
        let (swap, decided_by, unidentifiable) =
            decide_boundary(votes_same, votes_swap, adjacency_same, adjacency_swap);
        let window_gap_bp = window_axis
            .get(right.locus as usize)
            .zip(window_axis.get(left.locus as usize))
            .map(|(right_axis, left_axis)| right_axis.0 as i64 - left_axis.1 as i64);
        decisions.push(BoundaryDecision {
            locus_left: left.locus,
            locus_right: right.locus,
            window_gap_bp: window_gap_bp.unwrap_or(0),
            adjacency_same,
            adjacency_swap,
            votes_same,
            votes_swap,
            swap,
            decided_by,
            unidentifiable,
        });
    }
    let swaps: Vec<bool> = decisions.iter().map(|d| d.swap).collect();
    let states = accumulate_orientations(&swaps);

    // ------------------------------------------------- the two molecules
    // Per locus: molecule 1 carries the emitted fold
    // fold_indices[state ^ 0], molecule 2 the other. The spelled
    // sequences concatenate the per-locus fold rows (the fold's own
    // identical-sequence member rows), with the per-locus junction
    // accounting against the axis window extent.
    let spell_started = Instant::now();
    let mut sequence_cache: BTreeMap<(String, u64, u64), Vec<u8>> = BTreeMap::new();
    let mut molecules_json: Vec<serde_json::Value> = Vec::with_capacity(2);
    for molecule in 0..2u8 {
        let mut sequence: Vec<u8> = Vec::new();
        let mut locus_entries: Vec<serde_json::Value> = Vec::with_capacity(loci.len());
        for (index, locus) in loci.iter().enumerate() {
            let side = ((states[index] ^ molecule) & 1) as usize;
            let (name, start, end) = &locus.fold_row_head[side];
            let key = (name.clone(), *start, *end);
            let seq = match sequence_cache.get(&key) {
                Some(seq) => seq.clone(),
                None => {
                    let seq = fetch_seq(name, *start, *end)?;
                    sequence_cache.insert(key, seq.clone());
                    seq
                }
            };
            let axis_len = window_axis
                .get(locus.locus as usize)
                .map(|&(lo, hi)| hi.saturating_sub(lo))
                .unwrap_or(0);
            let fold_len = seq.len() as u64;
            locus_entries.push(serde_json::json!({
                "locus": locus.locus,
                "partition": locus.partition,
                "fold_index": locus.fold_indices[side],
                "members": locus.fold_rows[side].iter().map(|&(path, row_start, row_end)| {
                    serde_json::json!({
                        "path_name": panel.name_map.path_to_name[path],
                        "start": row_start,
                        "end": row_end,
                    })
                }).collect::<Vec<_>>(),
                "length": fold_len,
                "axis_window_length": axis_len,
                "junction_gap_bp": axis_len as i64 - fold_len as i64,
            }));
            sequence.extend_from_slice(&seq);
        }
        molecules_json.push(serde_json::json!({
            "molecule": molecule + 1,
            "length": sequence.len(),
            "loci": locus_entries,
            "sequence": String::from_utf8_lossy(&sequence),
        }));
    }
    eprintln!(
        "[chain] molecules spelled [{:.1}s]",
        spell_started.elapsed().as_secs_f64(),
    );
    rss.probe("spelling")?;

    // -------------------------------------------------------- the emission
    let out = options.out_dir.join("molecules.jsonl");
    let mut writer = BufWriter::new(File::create(&out)?);
    let boundaries_json: Vec<serde_json::Value> = decisions
        .iter()
        .enumerate()
        .map(|(index, decision)| {
            let slot = &evidence[index];
            let left = &loci[index];
            let right = &loci[index + 1];
            let per_link = serde_json::json!({
                "0,0": slot.votes[0][0], "0,1": slot.votes[0][1],
                "1,0": slot.votes[1][0], "1,1": slot.votes[1][1],
            });
            serde_json::json!({
                "loci": [decision.locus_left, decision.locus_right],
                "window_gap_bp": decision.window_gap_bp,
                "adjacency": {
                    "same": {"touches": decision.adjacency_same.0,
                             "overlap_bp": decision.adjacency_same.1},
                    "swap": {"touches": decision.adjacency_swap.0,
                             "overlap_bp": decision.adjacency_swap.1},
                },
                "crossing": {
                    "records": slot.crossing_records,
                    "occurrences": slot.crossing_occurrences,
                    "votes_same": decision.votes_same,
                    "votes_swap": decision.votes_swap,
                    "per_link": per_link,
                    "orphan_occurrences": slot.orphan_occurrences,
                    "unattributed_occurrences": slot.unattributed_occurrences,
                },
                "orientation": if decision.swap { "swap" } else { "same" },
                "decided_by": decision.decided_by,
                "unidentifiable": decision.unidentifiable,
                "link_confidence": {
                    "votes_same": decision.votes_same,
                    "votes_swap": decision.votes_swap,
                    "vote_margin": decision.votes_same - decision.votes_swap,
                    "adjacency_same": [decision.adjacency_same.0,
                                       decision.adjacency_same.1],
                    "adjacency_swap": [decision.adjacency_swap.0,
                                       decision.adjacency_swap.1],
                    "endpoint_qual": [left.qual, right.qual],
                    "endpoint_tied_class": [left.tied_classes, right.tied_classes],
                    "endpoint_homozygous": [left.homozygous, right.homozygous],
                    "orphan_occurrences": slot.orphan_occurrences,
                    "window_gap_bp": decision.window_gap_bp,
                },
            })
        })
        .collect();
    let loci_json: Vec<serde_json::Value> = loci
        .iter()
        .enumerate()
        .map(|(index, locus)| {
            serde_json::json!({
                "locus": locus.locus,
                "partition": locus.partition,
                "fold_indices": [locus.fold_indices[0], locus.fold_indices[1]],
                "orientation_bit": states[index],
                "molecule1_fold": locus.fold_indices[(states[index] & 1) as usize],
                "homozygous": locus.homozygous,
                "tied_classes": locus.tied_classes,
                "qual": locus.qual,
            })
        })
        .collect();
    let calls_fingerprint = crate::genome_inference::reconstruction::fingerprint(&options.calls)?;
    let record = serde_json::json!({
        "component": component,
        "model": MODEL,
        "calls": calls_fingerprint,
        "loci": loci_json,
        "boundaries": boundaries_json,
        "molecules": molecules_json,
    });
    serde_json::to_writer(&mut writer, &record)?;
    writer.write_all(b"\n")?;
    writer.flush()?;
    rss.probe("emission")?;

    // ------------------------------------ THE TEST MODE (assessment-side)
    // The product run's own truth-QV artifact, consumed read-only:
    // the switch table (chain orientation vs the committed
    // assignment yardstick, recomputed ties bracketing the boundary)
    // and the per-locus accuracy under best-case assignment CARRIED
    // UNCHANGED — separate columns, never conflated.
    if let Some(truth_path) = &options.truth_qv {
        let mut truth_records: BTreeMap<u32, TruthQvRecord> = BTreeMap::new();
        for line in BufReader::new(File::open(truth_path)?).lines() {
            let record: TruthQvRecord = serde_json::from_str(&line?)?;
            truth_records.insert(record.locus, record);
        }
        let assignments: BTreeMap<u32, TruthAssignment> = truth_records
            .iter()
            .map(|(&locus, record)| (locus, truth_assignment(record)))
            .collect();
        let mut assessable = 0u64;
        let mut switches = 0u64;
        let mut brackets: BTreeMap<&'static str, u64> = BTreeMap::new();
        let mut truth_lines: Vec<serde_json::Value> = Vec::new();
        for decision in &decisions {
            let left = &assignments[&decision.locus_left];
            let right = &assignments[&decision.locus_right];
            let chain_swap = decision.swap;
            let (truth_swap, bracket) = match (left, right) {
                (TruthAssignment::Swapped(a), TruthAssignment::Swapped(b)) => {
                    (a != b, None)
                }
                (TruthAssignment::Bracket(reason), TruthAssignment::Swapped(_))
                | (TruthAssignment::Swapped(_), TruthAssignment::Bracket(reason)) => {
                    (false, Some(*reason))
                }
                (TruthAssignment::Bracket(a), TruthAssignment::Bracket(b)) => {
                    (false, Some(if a == b { *a } else { "truth_pair_not_expressible" }))
                }
            };
            let switch_error = bracket.is_none() && chain_swap != truth_swap;
            if let Some(reason) = bracket {
                *brackets.entry(reason).or_insert(0) += 1;
            } else {
                assessable += 1;
            }
            if switch_error {
                switches += 1;
            }
            truth_lines.push(serde_json::json!({
                "loci": [decision.locus_left, decision.locus_right],
                "chain_orientation": if chain_swap { "swap" } else { "same" },
                "chain_decided_by": decision.decided_by,
                "chain_unidentifiable": decision.unidentifiable,
                "truth_orientation": match bracket {
                    None => serde_json::json!(if truth_swap { "swap" } else { "same" }),
                    Some(reason) => serde_json::json!({
                        "bracket": reason,
                    }),
                },
                "switch_error": switch_error,
                "assessable": bracket.is_none(),
                // The per-locus accuracy under best-case assignment,
                // carried UNCHANGED from the artifact (separate
                // columns, never conflated with the switch count).
                "endpoint_accuracy": [
                    accuracy_columns(truth_records.get(&decision.locus_left)),
                    accuracy_columns(truth_records.get(&decision.locus_right)),
                ],
            }));
        }
        let expressible = truth_records
            .values()
            .filter(|record| record.truth_pair_expressible)
            .count();
        let rank1 = truth_records
            .values()
            .filter(|record| record.rank1 == Some(true))
            .count();
        let perfect = truth_records
            .values()
            .filter(|record| record.perfect == Some(true) && record.truth_pair_expressible)
            .count();
        let summary = serde_json::json!({
            "component": component,
            "model": MODEL,
            "summary": true,
            "boundaries": decisions.len(),
            "assessable_boundaries": assessable,
            "switches": switches,
            "bracketed_boundaries": brackets,
            // The per-locus product accuracy, carried unchanged and
            // reported SEPARATELY (the owner's evaluation ruling:
            // per-partition accuracy under best-case assignment AND
            // switch errors, never conflated).
            "loci": truth_records.len(),
            "truth_pair_expressible": expressible,
            "truth_rank1": rank1,
            "perfect": perfect,
        });
        let truth_out = options.out_dir.join("molecules.jsonl.truth-qv.jsonl");
        let mut writer = BufWriter::new(File::create(&truth_out)?);
        for line in &truth_lines {
            serde_json::to_writer(&mut writer, line)?;
            writer.write_all(b"\n")?;
        }
        serde_json::to_writer(&mut writer, &summary)?;
        writer.write_all(b"\n")?;
        writer.flush()?;
        eprintln!(
            "[chain] test mode: {assessable} assessable boundaries, \
             {switches} switches",
        );
    }

    eprintln!(
        "[chain] complete: {} loci, {} boundaries [{:.1}s]",
        loci.len(),
        decisions.len(),
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}

/// The per-locus accuracy columns carried unchanged from the product
/// run's truth-QV artifact (the best-case-assignment yardstick's own
/// numbers; null where the artifact itself does not state them).
fn accuracy_columns(record: Option<&TruthQvRecord>) -> serde_json::Value {
    match record {
        None => serde_json::Value::Null,
        Some(record) => serde_json::json!({
            "locus": record.locus,
            "truth_pair_expressible": record.truth_pair_expressible,
            "truth_rank": record.truth_rank,
            "rank1": record.rank1,
            "error": record.error,
            "qv": record.qv,
            "perfect": record.perfect,
        }),
    }
}

// ---------------------------------------------------------------------------
// Unit tests (source-named, the pure core).
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adjacency_counts_touches_and_overlap_on_shared_paths() {
        let a = [(0usize, 10u64, 20u64), (1, 5, 15)];
        // Same path, exact junction touch (a ends where b starts).
        let b = [(0, 20, 30)];
        assert_eq!(fold_pair_adjacency(&a, &b), (1, 0));
        // Overlapping material on the same path.
        let c = [(0, 15, 25)];
        assert_eq!(fold_pair_adjacency(&a, &c), (0, 5));
        // Touch in the other direction (b ends where a starts).
        let d = [(0, 0, 10)];
        assert_eq!(fold_pair_adjacency(&a, &d), (1, 0));
        // Different paths: no continuation evidence.
        let e = [(2, 10, 20)];
        assert_eq!(fold_pair_adjacency(&a, &e), (0, 0));
    }

    #[test]
    fn rows_touch_material_overlaps_any_window() {
        let rows = [(100u64, 200u64)];
        assert!(rows_touch_material(&rows, &[(150, 213)]));
        assert!(!rows_touch_material(&rows, &[(200, 260)]));
        assert!(rows_touch_material(&rows, &[(199, 262)]));
    }

    #[test]
    fn boundary_decision_is_lexicographic_votes_then_adjacency() {
        // Votes dominate.
        assert_eq!(
            decide_boundary(3.0, 1.0, (0, 0), (9, 9)),
            (false, "crossing", false)
        );
        assert_eq!(
            decide_boundary(1.0, 3.0, (9, 9), (0, 0)),
            (true, "crossing", false)
        );
        // Exact vote ties fall to adjacency (touches before overlap).
        assert_eq!(
            decide_boundary(0.0, 0.0, (2, 0), (1, 50)),
            (false, "adjacency", false)
        );
        assert_eq!(
            decide_boundary(0.0, 0.0, (1, 0), (1, 50)),
            (true, "adjacency", false)
        );
        // Both channels tied: unidentifiable, carried as "same".
        assert_eq!(
            decide_boundary(0.0, 0.0, (1, 5), (1, 5)),
            (false, "tie", true)
        );
    }

    #[test]
    fn orientations_accumulate_through_swaps() {
        assert_eq!(accumulate_orientations(&[]), vec![0u8]);
        assert_eq!(
            accumulate_orientations(&[false, true, false, true, true]),
            vec![0, 0, 1, 1, 0, 1]
        );
    }

    #[test]
    fn truth_assignment_brackets_ties_and_inexpressible() {
        let base = |expressible, assignment, pairs, edits_direct, edits_crossed| TruthQvRecord {
            locus: 0,
            truth_pair_expressible: expressible,
            assignment,
            pairs,
            rank1: None,
            perfect: None,
            truth_rank: None,
            error: None,
            qv: None,
        };
        // Non-expressible: bracketed.
        assert!(matches!(
            truth_assignment(&base(false, None, None, 0, 0)),
            TruthAssignment::Bracket(_)
        ));
        // Direct, untied.
        let pairs = Some(TruthPairs {
            called0_truth0: TruthPairCounts { edits: 1, columns: 100 },
            called1_truth1: TruthPairCounts { edits: 1, columns: 100 },
            called0_truth1: TruthPairCounts { edits: 9, columns: 100 },
            called1_truth0: TruthPairCounts { edits: 9, columns: 100 },
        });
        assert!(matches!(
            truth_assignment(&base(true, Some("direct".into()), pairs.clone(), 2, 18)),
            TruthAssignment::Swapped(false)
        ));
        assert!(matches!(
            truth_assignment(&base(true, Some("crossed".into()), pairs, 2, 18)),
            TruthAssignment::Swapped(true)
        ));
        // A rate-and-edit tie: bracketed (the yardstick's own tie).
        let tied = Some(TruthPairs {
            called0_truth0: TruthPairCounts { edits: 5, columns: 100 },
            called1_truth1: TruthPairCounts { edits: 5, columns: 100 },
            called0_truth1: TruthPairCounts { edits: 5, columns: 100 },
            called1_truth0: TruthPairCounts { edits: 5, columns: 100 },
        });
        assert!(matches!(
            truth_assignment(&base(true, Some("direct".into()), tied, 10, 10)),
            TruthAssignment::Bracket("truth_assignment_tied")
        ));
    }
}
