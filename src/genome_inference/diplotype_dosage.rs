//! THE DOSAGE SURFACE — the copy-number layer of the two-molecule
//! emission (the owner's product chain: per-locus calls -> two
//! molecules -> the per-segment copy number the two molecules
//! imply). The dosage is DERIVED from the chain layer's
//! molecules.jsonl and from nothing else on the product side: the
//! same material, counted — never re-selected, never re-inferred.
//!
//! THE MODEL, derived and stated (no thresholds, no tuning
//! constants):
//!
//! 1. THE DOSAGE EMISSION (truth-free, product-side): per panel
//!    path, the two molecules' member rows decompose into
//!    ELEMENTARY SEGMENTS (the breakpoints are every row's start
//!    and end, so every row contains every elementary segment it
//!    touches — no partial overlaps exist by construction). A
//!    segment's COPY COUNT is the number of traversals of it by
//!    molecule 1 plus the number by molecule 2, one traversal per
//!    (molecule, locus, member row) containing the segment: a
//!    segment both molecules traverse is 2, one molecule 1, and
//!    REPEAT-VISIT MULTIPLICITY is carried — a molecule traversing
//!    a segment twice (two loci, or two rows of one fold) counts
//!    its true multiplicity. The copy-bp identity is exact: the
//!    segments tile the rows, so SUM(dosage * length) over the
//!    emission equals the sum of member-row lengths. Copy count 0
//!    is the complement by construction: the dosage domain is the
//!    molecules' own material (stated per path as covered_bp
//!    against the path length), and panel material outside it
//!    carries 0 copies. Per segment the territory/axis extents are
//!    stated: every visit carries its locus's axis window (the
//!    partition-graphs' own derivation, the same window spine the
//!    chain consumes), and the segment carries the merged axis
//!    extent of its visits.
//!
//! 2. THE DEPTH-CONSISTENCY QC (truth-free, a derived field — the
//!    ratio stated, never binned, no thresholds): the observed
//!    per-segment census mass vs the expected under the emitted
//!    dosage. The likelihood's own depth model is the balanced
//!    diploid prior — the class likelihood's EQUAL 1/2 MIXTURE
//!    WEIGHTS over the two homologs are the statement that each
//!    molecule's covered material carries the same expected read
//!    mass per covered bp, per covered copy. The per-covered-copy
//!    expectation is therefore uniform per bp: mu = attributed
//!    mass / emitted copy-bp, DERIVED per component from the
//!    census and the emission itself (never a tuned constant),
//!    and expected_mass(s) = mu * dosage(s) * length(s). The
//!    observed mass is the census's own record-once material
//!    attribution: a record's multiplicity-weighted covered bp —
//!    the union over its occurrences of the contained-node
//!    material windows on the occurrence's path (the same
//!    (bp, bp+k) windows the chain's edge layer walks) —
//!    attributed to a segment by overlap bp. The QC field is the
//!    ratio observed/expected per segment, stated to full
//!    precision; a segment emitted at 2 copies that carries only
//!    one copy's mass states its own ratio, and the aggregate
//!    states the totals (the census's multi-matching spread
//!    included as unattributed material mass — never hidden).
//!
//! 3. THE TEST MODE (assessment-side ONLY, the --truth-qv-file
//!    convention): the product run's own calls.jsonl.truth-qv.jsonl
//!    names the truth pair per locus (fold indices into the
//!    product run's own instrument receipt, whose fold identities
//!    carry the truth folds' member rows — the receipt is a
//!    truth-free product artifact; the truth REFERENCE is the test
//!    mode's own input). The truth pair's per-segment dosage is
//!    derived by the same rule as the emission, and the agreement
//!    table is per component: per locus the emitted pair's dosage
//!    vs the truth pair's over their union universe, with the
//!    honest classes the data name —
//!      * correct_dosage (the copy vectors agree),
//!      * dosage_inherited_from_non_rank1_call (the call itself is
//!        not rank-1; the dosage error is the call's, not the
//!        emission's; the mechanism stated: territory_split vs
//!        multiplicity_mismatch),
//!      * dosage_specific_emission_defect (the locus IS rank-1 and
//!        the dosage still differs — necessarily a winner-pair
//!        mismatch at a tie-degenerate locus, mechanistically a
//!        territory split (the segment sets differ) or a
//!        multiplicity mismatch (the same segments, different
//!        counts)),
//!      * non_rank1_call_dosage_agrees (a wrong call whose
//!        traversal structure coincides with the truth's — stated
//!        honestly, never counted as correct),
//!      * bracketed (the truth pair is not expressible at the
//!        locus; no truth dosage exists to compare);
//!    and per segment the component-wide agreement with the same
//!    decomposition (differences at rank-1 loci vs inherited from
//!    non-rank-1 calls vs segments a bracketed locus touches,
//!    which no truth total can judge). The per-locus product
//!    accuracy columns are CARRIED UNCHANGED from the artifact in
//!    separate columns, never conflated — the owner's evaluation
//!    ruling since the haploid era.
//!
//! 4. THE THIN-LAYER PROPERTY: the dosage NEVER writes the
//!    molecules (or the calls); it consumes molecules.jsonl and
//!    emits its own record carrying the consumed molecules'
//!    fingerprint, and every segment, visit, copy count, and
//!    per-locus dosage view is reproducible BIT-FOR-BIT from
//!    molecules.jsonl alone (the gate re-derives them in
//!    independence and compares). The default output carries ZERO
//!    truth keys.

use crate::genome_inference::{PanelIdentity};
use crate::syng::{SyncmerParams, SyngIndex};
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

pub const MODEL: &str = "diplotype-dosage-emission-v1";

/// The sample's uniform read length (the census occurrence's
/// covered hull is bounded by it — the chain layer's own constant).
const READ_LENGTH: u64 = 150;

pub struct Options {
    /// The chain layer's two-molecule product (molecules.jsonl) —
    /// the dosage is derived from it and nothing else product-side.
    pub molecules: PathBuf,
    /// The panel syng prefix (path name/id resolution, node walks)
    pub panel: String,
    /// The locality map directory (window<N>.gfa.map.json; the axis
    /// window extents — the chain's own derivation)
    pub partition_graphs: PathBuf,
    /// The multi-matching census receipt (the same truth-free input
    /// the product and chain runs consumed; the observed mass)
    pub census: PathBuf,
    /// THE TEST MODE (assessment-side only): the product run's
    /// calls.jsonl.truth-qv.jsonl — the truth pair per locus,
    /// consumed read-only into dosage.jsonl.truth-qv.jsonl
    pub truth_qv_file: Option<PathBuf>,
    /// The product run's instrument receipt (fold_identities: the
    /// per-locus fold member rows — the truth folds' material; a
    /// truth-free product artifact, required by the test mode)
    pub folds: Option<PathBuf>,
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
        eprintln!("[dosage] rss {stage}: {rss_kb} kB");
        Ok(rss_kb)
    }
}

// ---------------------------------------------------------------------------
// The consumed molecules (the chain layer's product record).
// ---------------------------------------------------------------------------

#[derive(Clone, Deserialize)]
struct MemberRow {
    path_name: String,
    start: u64,
    end: u64,
}

#[derive(Deserialize)]
struct MoleculeLocus {
    locus: u32,
    partition: u32,
    fold_index: u64,
    members: Vec<MemberRow>,
}

#[derive(Deserialize)]
struct Molecule {
    #[allow(dead_code)]
    molecule: u64,
    loci: Vec<MoleculeLocus>,
}

#[derive(Deserialize)]
struct LocusMeta {
    locus: u32,
    partition: u32,
    fold_indices: [u64; 2],
    #[serde(default)]
    homozygous: bool,
    #[serde(default)]
    tied_classes: bool,
    #[serde(default)]
    qual: Option<f64>,
}

#[derive(Deserialize)]
struct MoleculesRecord {
    component: String,
    #[allow(dead_code)]
    model: String,
    loci: Vec<LocusMeta>,
    molecules: Vec<Molecule>,
}

// ---------------------------------------------------------------------------
// The pure dosage core (unit-tested).
// ---------------------------------------------------------------------------

/// Row material on one panel path: (start, end).
type Rows = Vec<(u64, u64)>;
/// Rows keyed by panel path name.
type RowsByPath = BTreeMap<String, Rows>;

/// The elementary segment decomposition of row material per path:
/// the breakpoints are every row's start and end, and the
/// segments are the COVERED consecutive breakpoint pairs — a
/// row's own bounds bound it, so every row CONTAINS every
/// elementary segment it touches (no partial overlaps exist by
/// construction, and the gaps between disjoint rows are not
/// material: they carry 0 copies by construction, not segments).
fn elementary_segments(per_path: &RowsByPath) -> Vec<(String, u64, u64)> {
    let mut segments = Vec::new();
    for (path, rows) in per_path {
        let mut bounds: Vec<u64> = Vec::with_capacity(rows.len() * 2);
        for &(start, end) in rows {
            bounds.push(start);
            bounds.push(end);
        }
        bounds.sort_unstable();
        bounds.dedup();
        // The union of the rows: only the covered consecutive
        // breakpoint pairs are segments.
        let mut covered = rows.clone();
        covered.sort_unstable();
        merge_intervals(&mut covered);
        let mut bound_index = 0usize;
        for &(lo, hi) in &covered {
            while bound_index < bounds.len() && bounds[bound_index] < lo {
                bound_index += 1;
            }
            let first = bound_index;
            while bound_index < bounds.len() && bounds[bound_index] <= hi {
                bound_index += 1;
            }
            for pair in bounds[first..bound_index].windows(2) {
                segments.push((path.clone(), pair[0], pair[1]));
            }
        }
    }
    segments
}

/// The number of rows containing one elementary segment (one
/// traversal per containing row — the repeat-visit multiplicity).
fn containing_rows(rows: &Rows, start: u64, end: u64) -> u64 {
    rows.iter().filter(|&&(s, e)| s <= start && e >= end).count() as u64
}

/// The dosage view of one locus's pair of folds over the pair's own
/// restricted universe: per elementary segment the copy count
/// (fold-1 traversals + fold-2 traversals), the segment count, the
/// covered bp, and the flat-class-of-2 flag (every covered segment
/// carried by BOTH folds — one coverage class at 2 copies).
#[derive(Clone, Default)]
struct LocusDosageView {
    segment_count: u64,
    covered_bp: u64,
    flat_class_of_2: bool,
    /// copy-count histogram: copies -> segments at that count.
    histogram: BTreeMap<u64, u64>,
}

fn locus_dosage_view(fold_rows: [&RowsByPath; 2]) -> LocusDosageView {
    let mut per_path: RowsByPath = RowsByPath::new();
    for fold in &fold_rows {
        for (path, rows) in fold.iter() {
            per_path.entry(path.clone()).or_default().extend(rows);
        }
    }
    let segments = elementary_segments(&per_path);
    let mut view = LocusDosageView {
        segment_count: segments.len() as u64,
        flat_class_of_2: !segments.is_empty(),
        ..Default::default()
    };
    for (path, start, end) in &segments {
        let copies = fold_rows
            .iter()
            .map(|fold| {
                fold.get(path)
                    .map(|rows| containing_rows(rows, *start, *end))
                    .unwrap_or(0)
            })
            .sum::<u64>();
        view.covered_bp += end - start;
        *view.histogram.entry(copies).or_default() += 1;
        if copies != 2 {
            view.flat_class_of_2 = false;
        }
    }
    view
}

/// The expected mass under the emitted dosage: the per-covered-copy
/// expectation (uniform per bp — the balanced diploid prior's own
/// statement) scaled by the segment's copy count and length. The
/// normalizer mu is DERIVED per component (attributed mass / emitted
/// copy-bp), never tuned.
fn expected_mass(mu: f64, dosage: u64, length: u64) -> f64 {
    mu * (dosage * length) as f64
}

/// The number of entries before the first for which `predicate` is
/// false (the standard partition_point over a sorted vec).
fn partition_point<T>(sorted: &[T], predicate: impl Fn(&T) -> bool) -> usize {
    let mut lo = 0usize;
    let mut hi = sorted.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if predicate(&sorted[mid]) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

/// Merge the overlapping/touching intervals of a sorted list in
/// place (the union of a record's material windows).
fn merge_intervals(spans: &mut Vec<(u64, u64)>) {
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(spans.len());
    for &(lo, hi) in spans.iter() {
        match merged.last_mut() {
            Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
            _ => merged.push((lo, hi)),
        }
    }
    *spans = merged;
}

// ---------------------------------------------------------------------------
// The test mode (assessment-side only): the truth-QV artifact and
// the product receipt's fold identities.
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct TruthFoldRef {
    index: u64,
}

#[derive(Deserialize)]
struct TruthQvRecord {
    locus: u32,
    #[serde(default)]
    truth_pair_expressible: bool,
    #[serde(default)]
    rank1: Option<bool>,
    #[serde(default)]
    perfect: Option<bool>,
    #[serde(default)]
    truth_folds: Option<Vec<TruthFoldRef>>,
}

#[derive(Deserialize)]
struct FoldIdentity {
    #[allow(dead_code)]
    length: u64,
    members: Vec<MemberRow>,
}

#[derive(Deserialize)]
struct ReceiptLocus {
    locus: u32,
    fold_identities: Vec<FoldIdentity>,
}

// ---------------------------------------------------------------------------
// The dosage run.
// ---------------------------------------------------------------------------

pub fn run(options: Options) -> io::Result<()> {
    let started = Instant::now();
    let rss = RssGuard::new(options.rss_budget_gib);
    ensure(options.molecules.is_file(), "the molecules file is absent")?;
    let text = std::fs::read_to_string(&options.molecules)?;
    let record: MoleculesRecord = serde_json::from_str(&text)?;
    ensure(record.molecules.len() == 2, "the chain emits two molecules")?;
    let component = record.component.clone();
    ensure(
        record.molecules[0].loci.len() == record.molecules[1].loci.len()
            && !record.molecules[0].loci.is_empty(),
        "the two molecules disagree on the locus set",
    )?;
    for (index, meta) in record.loci.iter().enumerate() {
        ensure(
            record.molecules[0].loci[index].locus == meta.locus
                && record.molecules[1].loci[index].locus == meta.locus,
            "the molecules and the locus table disagree on a locus",
        )?;
    }
    ensure(
        record
            .molecules
            .iter()
            .flat_map(|m| m.loci.iter())
            .all(|l| l.members.iter().all(|m| m.start < m.end)),
        "a member row is empty",
    )?;
    ensure(
        options.truth_qv_file.is_some() == options.folds.is_some(),
        "the test mode requires BOTH --truth-qv-file and the fold identities receipt",
    )?;
    rss.probe("molecules")?;

    // The panel: path name/id resolution and lengths (the census
    // walks; the test mode's row capping — the chain's own
    // convention).
    let _identity = PanelIdentity::read(&options.panel)?;
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    let k = panel.syncmer_length_bp() as u64;
    let path_of_name: HashMap<String, usize> = panel
        .name_map
        .name_to_path
        .iter()
        .map(|(name, &path)| (name.clone(), path as usize))
        .collect();
    rss.probe("panel")?;

    // The window axis extents (the chain's own derivation: the
    // component's member rows ranked by start map to the window ids
    // in order — one window per axis row).
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
        for axis in map.members.iter().filter(|m| m.path_name == component) {
            axis_rows.push((map.partition, axis.start, axis.end));
        }
    }
    axis_rows.sort_by_key(|&(partition, start, end)| (start, end, partition));
    let window_axis: Vec<(u64, u64)> =
        axis_rows.iter().map(|&(_, start, end)| (start, end)).collect();
    ensure(
        record
            .molecules
            .iter()
            .flat_map(|m| m.loci.iter())
            .all(|l| (l.locus as usize) < window_axis.len()),
        "a molecule locus is outside the axis window span",
    )?;
    rss.probe("axis")?;

    // ---------------------------------------------------- the derivation
    // Per molecule, per locus: the carried fold's rows keyed by
    // panel path. The PAIR at a locus is molecule 1's entry plus
    // molecule 2's entry.
    let fold_rows: [Vec<RowsByPath>; 2] = [
        record.molecules[0]
            .loci
            .iter()
            .map(|locus| rows_by_path(&locus.members))
            .collect(),
        record.molecules[1]
            .loci
            .iter()
            .map(|locus| rows_by_path(&locus.members))
            .collect(),
    ];
    // The component universe: every row of both molecules.
    let mut universe_paths: RowsByPath = RowsByPath::new();
    for molecule in &fold_rows {
        for locus_rows in molecule {
            for (path, rows) in locus_rows {
                universe_paths
                    .entry(path.clone())
                    .or_default()
                    .extend(rows);
            }
        }
    }
    let universe = elementary_segments(&universe_paths);
    // Per path: the sorted segment bounds (the census attribution's
    // binary search) and the per-segment observed-mass accumulator.
    struct PathSegments {
        starts: Vec<u64>,
        ends: Vec<u64>,
        mass: Vec<u64>,
    }
    let mut per_path: Vec<(String, PathSegments)> = Vec::new();
    let mut slot_of_name: HashMap<&str, usize> = HashMap::new();
    // The global universe index of each slot's first segment (the
    // segments within a path are contiguous in the universe vector).
    let mut slot_base: Vec<usize> = Vec::new();
    // Per universe segment: its slot (path) and position within it.
    let mut segment_slot: Vec<usize> = Vec::with_capacity(universe.len());
    let mut segment_pos: Vec<usize> = Vec::with_capacity(universe.len());
    for (index, (path, start, end)) in universe.iter().enumerate() {
        let slot = match slot_of_name.get(path.as_str()) {
            Some(&slot) => slot,
            None => {
                per_path.push((
                    path.clone(),
                    PathSegments {
                        starts: Vec::new(),
                        ends: Vec::new(),
                        mass: Vec::new(),
                    },
                ));
                slot_base.push(index);
                slot_of_name.insert(path.as_str(), per_path.len() - 1);
                per_path.len() - 1
            }
        };
        let segments = &mut per_path[slot].1;
        segment_slot.push(slot);
        segment_pos.push(segments.starts.len());
        segments.starts.push(*start);
        segments.ends.push(*end);
        segments.mass.push(0);
    }
    rss.probe("derivation")?;

    // The visits: one per (molecule, locus, row) containing the
    // segment — the repeat-visit multiplicity is the traversal count.
    struct Visit {
        molecule: u8,
        locus: u32,
        partition: u32,
        fold_index: u64,
        axis_window: (u64, u64),
    }
    let mut visits: Vec<Vec<Visit>> =
        (0..universe.len()).map(|_| Vec::new()).collect();
    for (molecule_index, molecule) in fold_rows.iter().enumerate() {
        for (locus_index, locus_rows) in molecule.iter().enumerate() {
            let entry = &record.molecules[molecule_index].loci[locus_index];
            let axis_window = window_axis
                .get(entry.locus as usize)
                .copied()
                .unwrap_or((0, 0));
            for (path, rows) in locus_rows {
                let Some(&slot) = slot_of_name.get(path.as_str()) else {
                    continue;
                };
                let starts = &per_path[slot].1.starts;
                let ends = &per_path[slot].1.ends;
                let base = slot_base[slot];
                for &(row_start, row_end) in rows {
                    // The segments a row contains are a contiguous
                    // sorted index range within the path.
                    let lo = partition_point(starts, |&s| s < row_start);
                    let hi = partition_point(starts, |&s| s < row_end);
                    for position in lo..hi {
                        if ends[position] <= row_end {
                            visits[base + position].push(Visit {
                                molecule: molecule_index as u8 + 1,
                                locus: entry.locus,
                                partition: entry.partition,
                                fold_index: entry.fold_index,
                                axis_window,
                            });
                        }
                    }
                }
            }
        }
    }
    let dosage: Vec<u64> = visits
        .iter()
        .map(|v| v.len() as u64)
        .collect();
    // THE COPY-BP IDENTITY: the segments tile the rows, so the
    // emitted copy-bp equals the summed member-row lengths exactly.
    let mut row_bp = 0u64;
    for molecule in &fold_rows {
        for locus_rows in molecule {
            for rows in locus_rows.values() {
                for &(start, end) in rows {
                    row_bp += end - start;
                }
            }
        }
    }
    let copy_bp: u64 = universe
        .iter()
        .zip(&dosage)
        .map(|((_, start, end), &copies)| copies * (end - start))
        .sum();
    ensure(copy_bp == row_bp, "the copy-bp identity failed")?;
    rss.probe("visits")?;

    // ---------------------------------------- the depth-consistency QC
    // The observed per-segment mass: record-once — a record's
    // multiplicity-weighted covered bp, the union over its
    // occurrences of the contained-node material windows on the
    // occurrence's path, attributed per segment by overlap bp.
    let census_started = Instant::now();
    #[derive(Deserialize)]
    struct Occurrence {
        path: usize,
        start: u64,
        #[serde(default)]
        intervals: Vec<Vec<u32>>,
    }
    #[derive(Deserialize)]
    struct CensusLine {
        multiplicity: u64,
        #[serde(default)]
        occurrences: Vec<Occurrence>,
    }
    let (reader, _) = niffler::get_reader(Box::new(File::open(&options.census)?))
        .map_err(io::Error::other)?;
    let mut census_records = 0u64;
    let mut census_mass = 0u64;
    let mut material_mass = 0u64;
    let mut attributed_mass = 0u64;
    for line in BufReader::new(reader).lines() {
        let record: CensusLine = serde_json::from_str(&line?)?;
        census_records += 1;
        census_mass += record.multiplicity;
        // The union of the record's material windows, per path,
        // merged across occurrences (record-once).
        let mut windows: BTreeMap<usize, Vec<(u64, u64)>> = BTreeMap::new();
        for occurrence in &record.occurrences {
            let nodes: HashSet<u32> = occurrence.intervals.iter().flatten().copied().collect();
            if nodes.is_empty() {
                continue;
            }
            let mut spans: Vec<(u64, u64)> = panel
                .walk_path_range(
                    occurrence.path,
                    occurrence.start,
                    occurrence.start + READ_LENGTH,
                )?
                .into_iter()
                .filter(|&(node, _)| nodes.contains(&node.unsigned_abs()))
                .map(|(_, bp)| (bp, bp + k))
                .collect();
            if spans.is_empty() {
                continue;
            }
            spans.sort_unstable();
            windows.entry(occurrence.path).or_default().extend(spans);
        }
        for spans in windows.values_mut() {
            spans.sort_unstable();
            spans.dedup();
            merge_intervals(spans);
        }
        for (path, spans) in &windows {
            let Some(name) = panel
                .name_map
                .path_to_name
                .get(*path)
                .map(|n| n.as_str())
            else {
                continue;
            };
            let Some(&slot) = slot_of_name.get(name) else {
                continue;
            };
            let segments = &mut per_path[slot].1;
            for &(lo, hi) in spans {
                material_mass += record.multiplicity * (hi - lo);
                // The segments overlapping [lo, hi): a contiguous
                // sorted index range from the first candidate.
                let mut index =
                    partition_point(&segments.starts, |&s| s < lo).saturating_sub(1);
                while index < segments.starts.len() && segments.starts[index] < hi {
                    let overlap_lo = segments.starts[index].max(lo);
                    let overlap_hi = segments.ends[index].min(hi);
                    if overlap_hi > overlap_lo {
                        let mass = record.multiplicity * (overlap_hi - overlap_lo);
                        segments.mass[index] += mass;
                        attributed_mass += mass;
                    }
                    index += 1;
                }
            }
        }
    }
    // mu: the DERIVED per-covered-copy per-bp expectation.
    let mu = if copy_bp > 0 {
        attributed_mass as f64 / copy_bp as f64
    } else {
        0.0
    };
    eprintln!(
        "[dosage] census streamed: {census_records} records, attributed {attributed_mass} of {material_mass} material mass [{:.1}s]",
        census_started.elapsed().as_secs_f64(),
    );
    rss.probe("census")?;

    // -------------------------------------------- the per-locus views
    let mut locus_views: Vec<LocusDosageView> = Vec::with_capacity(record.loci.len());
    for index in 0..record.loci.len() {
        let pair = [
            rows_by_path(&record.molecules[0].loci[index].members),
            rows_by_path(&record.molecules[1].loci[index].members),
        ];
        locus_views.push(locus_dosage_view([&pair[0], &pair[1]]));
    }

    // -------------------------------------------------- the emission
    let out = options.out_dir.join("dosage.jsonl");
    let mut writer = BufWriter::new(File::create(&out)?);
    let segments_json: Vec<serde_json::Value> = universe
        .iter()
        .enumerate()
        .map(|(index, (path, start, end))| {
            let observed = per_path[segment_slot[index]].1.mass[segment_pos[index]];
            let expected = expected_mass(mu, dosage[index], end - start);
            serde_json::json!({
                "path_name": path,
                "start": start,
                "end": end,
                "length": end - start,
                "dosage": dosage[index],
                "visits_molecule1":
                    visits[index].iter().filter(|v| v.molecule == 1).count(),
                "visits_molecule2":
                    visits[index].iter().filter(|v| v.molecule == 2).count(),
                "visits": visits[index].iter().map(|v| serde_json::json!({
                    "molecule": v.molecule,
                    "locus": v.locus,
                    "partition": v.partition,
                    "fold_index": v.fold_index,
                    "axis_window": [v.axis_window.0, v.axis_window.1],
                })).collect::<Vec<_>>(),
                "axis_extent": [
                    visits[index].iter().map(|v| v.axis_window.0).min()
                        .unwrap_or(*start),
                    visits[index].iter().map(|v| v.axis_window.1).max()
                        .unwrap_or(*end),
                ],
                "observed_mass": observed,
                "expected_mass": expected,
                "mass_ratio": if expected > 0.0 {
                    serde_json::json!(observed as f64 / expected)
                } else {
                    serde_json::Value::Null
                },
            })
        })
        .collect();
    let paths_json: Vec<serde_json::Value> = per_path
        .iter()
        .map(|(path, segments)| {
            let covered_bp: u64 = segments
                .starts
                .iter()
                .zip(&segments.ends)
                .map(|(&s, &e)| e - s)
                .sum();
            let copy_bp: u64 = universe
                .iter()
                .zip(&dosage)
                .filter(|((p, _, _), _)| p == path)
                .map(|((_, s, e), &c)| c * (e - s))
                .sum();
            let path_length = path_of_name
                .get(path)
                .map(|&idx| panel.name_map.path_to_length[idx])
                .unwrap_or(0);
            serde_json::json!({
                "path_name": path,
                "path_length": path_length,
                "segments": segments.starts.len(),
                "covered_bp": covered_bp,
                "copy_bp": copy_bp,
                "zero_copy_bp": path_length.saturating_sub(covered_bp),
            })
        })
        .collect();
    let loci_json: Vec<serde_json::Value> = record
        .loci
        .iter()
        .zip(&locus_views)
        .map(|(meta, view)| {
            serde_json::json!({
                "locus": meta.locus,
                "partition": meta.partition,
                "fold_indices": [meta.fold_indices[0], meta.fold_indices[1]],
                "homozygous_pair": meta.fold_indices[0] == meta.fold_indices[1],
                "homozygous_flag": meta.homozygous,
                "tied_classes": meta.tied_classes,
                "qual": meta.qual,
                "segment_count": view.segment_count,
                "covered_bp": view.covered_bp,
                "flat_class_of_2": view.flat_class_of_2,
                "dosage_histogram": view.histogram,
            })
        })
        .collect();
    let molecules_fingerprint =
        crate::genome_inference::reconstruction::fingerprint(&options.molecules)?;
    let emission = serde_json::json!({
        "component": component,
        "model": MODEL,
        "molecules": molecules_fingerprint,
        "census_material_window_bp": k,
        "loci": loci_json,
        "paths": paths_json,
        "segments": segments_json,
        "depth_qc": {
            "per_covered_copy_bp_expectation": mu,
            "census_records": census_records,
            "census_mass": census_mass,
            "material_mass": material_mass,
            "attributed_mass": attributed_mass,
            "unattributed_material_mass": material_mass.saturating_sub(attributed_mass),
            "emitted_copy_bp": copy_bp,
            "row_bp": row_bp,
        },
    });
    serde_json::to_writer(&mut writer, &emission)?;
    writer.write_all(b"\n")?;
    writer.flush()?;
    rss.probe("emission")?;

    // ------------------------------------ THE TEST MODE (assessment-side)
    if let (Some(truth_path), Some(folds_path)) = (&options.truth_qv_file, &options.folds) {
        let mut truth: BTreeMap<u32, TruthQvRecord> = BTreeMap::new();
        for line in BufReader::new(File::open(truth_path)?).lines() {
            let record: TruthQvRecord = serde_json::from_str(&line?)?;
            truth.insert(record.locus, record);
        }
        let mut receipts: BTreeMap<u32, ReceiptLocus> = BTreeMap::new();
        for line in BufReader::new(File::open(folds_path)?).lines() {
            let record: ReceiptLocus = serde_json::from_str(&line?)?;
            receipts.insert(record.locus, record);
        }
        ensure(
            truth.len() == record.loci.len()
                && truth.len() == receipts.len()
                && truth.keys().eq(record.loci.iter().map(|m| &m.locus)),
            "the test-mode artifacts do not cover the locus set",
        )?;

        // The truth pair's rows per locus (the receipt's fold
        // identities at the truth indices; the chain's own capping
        // convention at the path length).
        let cap_rows = |members: &[MemberRow]| -> RowsByPath {
            let mut per_path: RowsByPath = RowsByPath::new();
            for member in members {
                if let Some(&path_idx) = path_of_name.get(&member.path_name) {
                    let end = member
                        .end
                        .min(panel.name_map.path_to_length[path_idx]);
                    if member.start < end {
                        per_path
                            .entry(member.path_name.clone())
                            .or_default()
                            .push((member.start, end));
                    }
                }
            }
            per_path
        };
        // The component-wide union universe: the emitted material
        // plus the truth pairs' material at the expressible loci.
        let mut union_paths: RowsByPath = universe_paths.clone();
        let mut truth_pair_rows: Vec<[RowsByPath; 2]> =
            vec![Default::default(); record.loci.len()];
        let mut truth_homozygous: Vec<bool> = vec![false; record.loci.len()];
        for (index, meta) in record.loci.iter().enumerate() {
            let entry = &truth[&meta.locus];
            if !entry.truth_pair_expressible {
                continue;
            }
            let indices: Vec<u64> = entry
                .truth_folds
                .as_ref()
                .map(|folds| folds.iter().map(|f| f.index).collect())
                .unwrap_or_default();
            ensure(indices.len() == 2, "the truth pair is not a pair")?;
            truth_homozygous[index] = indices[0] == indices[1];
            let receipt = &receipts[&meta.locus];
            for side in 0..2 {
                let fold = indices[side] as usize;
                ensure(
                    fold < receipt.fold_identities.len(),
                    "a truth fold index is outside the receipt's fold list",
                )?;
                let rows = cap_rows(&receipt.fold_identities[fold].members);
                for (path, material) in &rows {
                    union_paths
                        .entry(path.clone())
                        .or_default()
                        .extend(material);
                }
                truth_pair_rows[index][side] = rows;
            }
        }
        let union = elementary_segments(&union_paths);
        // Per union segment: its path slot and position (the row
        // enumeration's binary search; the segments within a path
        // are contiguous in the union vector, so the global index is
        // the slot's base plus the position).
        let mut union_slot_of_name: HashMap<&str, usize> = HashMap::new();
        let mut union_bounds: Vec<(Vec<u64>, Vec<u64>)> = Vec::new();
        let mut union_slot_base: Vec<usize> = Vec::new();
        for (index, (path, start, end)) in union.iter().enumerate() {
            let slot = match union_slot_of_name.get(path.as_str()) {
                Some(&slot) => slot,
                None => {
                    union_bounds.push((Vec::new(), Vec::new()));
                    union_slot_base.push(index);
                    union_slot_of_name.insert(path.as_str(), union_bounds.len() - 1);
                    union_bounds.len() - 1
                }
            };
            union_bounds[slot].0.push(*start);
            union_bounds[slot].1.push(*end);
        }
        // The union segments a row contains (one index per
        // traversal — the same enumeration the emission uses; the
        // segments within a path are contiguous in the union vector,
        // so the global index is the slot's base plus the position).
        let union_indices_of_row = |path: &str, row_start: u64, row_end: u64| -> Vec<usize> {
            let Some(&slot) = union_slot_of_name.get(path) else {
                return Vec::new();
            };
            let (starts, ends) = (&union_bounds[slot].0, &union_bounds[slot].1);
            let lo = partition_point(starts, |&s| s < row_start);
            let hi = partition_point(starts, |&s| s < row_end);
            (lo..hi)
                .filter(|&pos| ends[pos] <= row_end)
                .map(|pos| union_slot_base[slot] + pos)
                .collect()
        };

        // Per locus: the emitted pair's dosage vs the truth pair's
        // over the pair union universe; the per-(locus, segment)
        // copies feed the component-wide table.
        struct LocusComparison {
            class: &'static str,
            mechanism: &'static str,
            truth_flat2: bool,
            truth_homozygous: bool,
            equal: bool,
            only_emitted: u64,
            only_truth: u64,
            count_mismatch: u64,
        }
        let mut comparisons: Vec<LocusComparison> =
            Vec::with_capacity(record.loci.len());
        let mut emitted_locus_copies: HashMap<(usize, usize), u64> = HashMap::new();
        let mut truth_locus_copies: HashMap<(usize, usize), u64> = HashMap::new();
        for (index, meta) in record.loci.iter().enumerate() {
            let emitted = [&fold_rows[0][index], &fold_rows[1][index]];
            let entry = &truth[&meta.locus];
            let truth_folds = &truth_pair_rows[index];
            let mut pair_paths: RowsByPath = RowsByPath::new();
            for fold in emitted.iter().copied().chain(truth_folds.iter()) {
                for (path, rows) in fold.iter() {
                    pair_paths.entry(path.clone()).or_default().extend(rows);
                }
            }
            let pair_universe = elementary_segments(&pair_paths);
            let mut equal = true;
            let mut only_emitted = 0u64;
            let mut only_truth = 0u64;
            let mut count_mismatch = 0u64;
            for (path, start, end) in &pair_universe {
                let emitted_copies: u64 = emitted
                    .iter()
                    .map(|fold| {
                        fold.get(path)
                            .map(|rows| containing_rows(rows, *start, *end))
                            .unwrap_or(0)
                    })
                    .sum();
                let truth_copies: u64 = if entry.truth_pair_expressible {
                    truth_folds
                        .iter()
                        .map(|fold| {
                            fold.get(path)
                                .map(|rows| containing_rows(rows, *start, *end))
                                .unwrap_or(0)
                        })
                        .sum()
                } else {
                    0
                };
                if !entry.truth_pair_expressible {
                    continue;
                }
                if emitted_copies > 0 && truth_copies == 0 {
                    only_emitted += 1;
                    equal = false;
                } else if emitted_copies == 0 && truth_copies > 0 {
                    only_truth += 1;
                    equal = false;
                } else if emitted_copies != truth_copies {
                    count_mismatch += 1;
                    equal = false;
                }
            }
            // The per-(locus, union-segment) copies for the
            // component-wide table (rows enumerate their contained
            // union segments; one count per containing row).
            for side in 0..2 {
                for (path, rows) in emitted[side].iter() {
                    for &(row_start, row_end) in rows {
                        for union_index in union_indices_of_row(path, row_start, row_end)
                        {
                            *emitted_locus_copies
                                .entry((index, union_index))
                                .or_default() += 1;
                        }
                    }
                }
            }
            if entry.truth_pair_expressible {
                for side in 0..2 {
                    for (path, rows) in &truth_folds[side] {
                        for &(row_start, row_end) in rows {
                            for union_index in
                                union_indices_of_row(path, row_start, row_end)
                            {
                                *truth_locus_copies
                                    .entry((index, union_index))
                                    .or_default() += 1;
                            }
                        }
                    }
                }
            }
            let truth_view = if entry.truth_pair_expressible {
                locus_dosage_view([
                    &truth_pair_rows[index][0],
                    &truth_pair_rows[index][1],
                ])
            } else {
                LocusDosageView::default()
            };            let rank1 = entry.rank1 == Some(true);
            let mechanism = if only_emitted > 0 || only_truth > 0 {
                "territory_split"
            } else {
                "multiplicity_mismatch"
            };
            let (class, mechanism) = if !entry.truth_pair_expressible {
                ("bracketed_truth_pair_not_expressible", "none")
            } else if rank1 && equal {
                ("correct_dosage", "none")
            } else if rank1 {
                ("dosage_specific_emission_defect", mechanism)
            } else if equal {
                ("non_rank1_call_dosage_agrees", "none")
            } else {
                ("dosage_inherited_from_non_rank1_call", mechanism)
            };
            comparisons.push(LocusComparison {
                class,
                mechanism,
                truth_flat2: truth_view.flat_class_of_2,
                truth_homozygous: truth_homozygous[index],
                equal,
                only_emitted,
                only_truth,
                count_mismatch,
            });
        }

        // The component-wide per-segment agreement over the union
        // universe, with the honest decomposition: assessable iff
        // every locus touching the segment is expressible. The
        // touching loci per segment are collected once from the
        // per-(locus, segment) copy maps.
        let mut touching: Vec<Vec<usize>> = vec![Vec::new(); union.len()];
        for (&(locus_index, union_index), _) in emitted_locus_copies.iter() {
            touching[union_index].push(locus_index);
        }
        for (&(locus_index, union_index), _) in truth_locus_copies.iter() {
            if emitted_locus_copies
                .get(&(locus_index, union_index))
                .is_none()
            {
                touching[union_index].push(locus_index);
            }
        }
        for list in touching.iter_mut() {
            list.sort_unstable();
            list.dedup();
        }
        let mut segment_classes: BTreeMap<&'static str, u64> = BTreeMap::new();
        let mut differing_segments: Vec<serde_json::Value> = Vec::new();
        for (union_index, (path, start, end)) in union.iter().enumerate() {
            let touching = &touching[union_index];
            if touching
                .iter()
                .any(|&index| !truth[&record.loci[index].locus].truth_pair_expressible)
            {
                *segment_classes
                    .entry("bracketed_locus_involved")
                    .or_default() += 1;
                continue;
            }
            let emitted_total: u64 = touching
                .iter()
                .map(|&index| {
                    emitted_locus_copies
                        .get(&(index, union_index))
                        .copied()
                        .unwrap_or(0)
                })
                .sum();
            let truth_total: u64 = touching
                .iter()
                .map(|&index| {
                    truth_locus_copies
                        .get(&(index, union_index))
                        .copied()
                        .unwrap_or(0)
                })
                .sum();
            if emitted_total == truth_total {
                *segment_classes.entry("agree").or_default() += 1;
                continue;
            }
            let rank1_divergence = touching.iter().any(|&index| {
                let locus = &record.loci[index];
                truth[&locus.locus].rank1 == Some(true)
                    && emitted_locus_copies.get(&(index, union_index)).copied().unwrap_or(0)
                        != truth_locus_copies
                            .get(&(index, union_index))
                            .copied()
                            .unwrap_or(0)
            });
            let class = if rank1_divergence {
                "dosage_specific_at_rank1_locus"
            } else {
                "inherited_from_non_rank1_calls"
            };
            *segment_classes.entry(class).or_default() += 1;
            if differing_segments.len() < 10_000 {
                differing_segments.push(serde_json::json!({
                    "path_name": path,
                    "start": start,
                    "end": end,
                    "emitted_copies": emitted_total,
                    "truth_copies": truth_total,
                    "class": class,
                    "loci": touching.iter()
                        .map(|&index| record.loci[index].locus).collect::<Vec<_>>(),
                }));
            }
        }

        // The per-locus table lines and the component summary.
        let mut class_counts: BTreeMap<&'static str, u64> = BTreeMap::new();
        for comparison in &comparisons {
            *class_counts.entry(comparison.class).or_default() += 1;
        }
        let expressible = truth
            .values()
            .filter(|entry| entry.truth_pair_expressible)
            .count() as u64;
        let rank1 = truth
            .values()
            .filter(|entry| entry.rank1 == Some(true))
            .count() as u64;
        let correct = *class_counts.get("correct_dosage").unwrap_or(&0);
        let rank1_agree = record
            .loci
            .iter()
            .zip(&comparisons)
            .filter(|(meta, comparison)| {
                truth[&meta.locus].rank1 == Some(true) && comparison.equal
            })
            .count() as u64;
        let emitted_flat2 = locus_views
            .iter()
            .filter(|view| view.flat_class_of_2)
            .count() as u64;
        let truth_flat2 = comparisons
            .iter()
            .filter(|c| c.truth_flat2)
            .count() as u64;
        let truth_out = options.out_dir.join("dosage.jsonl.truth-qv.jsonl");
        let mut writer = BufWriter::new(File::create(&truth_out)?);
        for ((meta, view), comparison) in
            record.loci.iter().zip(&locus_views).zip(&comparisons)
        {
            let entry = &truth[&meta.locus];
            serde_json::to_writer(
                &mut writer,
                &serde_json::json!({
                    "locus": meta.locus,
                    "partition": meta.partition,
                    "truth_pair_expressible": entry.truth_pair_expressible,
                    "rank1": entry.rank1,
                    "perfect": entry.perfect,
                    "emitted": {
                        "homozygous_pair":
                            meta.fold_indices[0] == meta.fold_indices[1],
                        "flat_class_of_2": view.flat_class_of_2,
                        "segment_count": view.segment_count,
                        "dosage_histogram": view.histogram,
                    },
                    "truth": if entry.truth_pair_expressible {
                        serde_json::json!({
                            "homozygous_pair": comparison.truth_homozygous,
                            "flat_class_of_2": comparison.truth_flat2,
                        })
                    } else {
                        serde_json::Value::Null
                    },
                    "dosage_equal": comparison.equal,
                    "class": comparison.class,
                    "mechanism": comparison.mechanism,
                    "segments_only_emitted": comparison.only_emitted,
                    "segments_only_truth": comparison.only_truth,
                    "segments_count_mismatch": comparison.count_mismatch,
                }),
            )?;
            writer.write_all(b"\n")?;
        }
        serde_json::to_writer(
            &mut writer,
            &serde_json::json!({
                "component": component,
                "model": MODEL,
                "summary": true,
                "loci": record.loci.len(),
                "truth_pair_expressible": expressible,
                "truth_rank1": rank1,
                "correct_dosage": correct,
                "correct_dosage_fraction": if expressible > 0 {
                    serde_json::json!(correct as f64 / expressible as f64)
                } else {
                    serde_json::Value::Null
                },
                "rank1_dosage_agree": rank1_agree,
                "class_counts": class_counts,
                "emitted_flat_class_of_2_loci": emitted_flat2,
                "truth_flat_class_of_2_loci": truth_flat2,
                "segment_classes": segment_classes,
                "union_segments": union.len(),
                "differing_segment_count": differing_segments.len(),
                "differing_segments": differing_segments,
            }),
        )?;
        writer.write_all(b"\n")?;
        writer.flush()?;
        eprintln!(
            "[dosage] test mode: {correct}/{expressible} expressible agree, \
             {rank1_agree}/{rank1} rank-1 agree",
        );
    }

    eprintln!(
        "[dosage] complete: {} segments, copy-bp {copy_bp} [{:.1}s]",
        universe.len(),
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}


fn rows_by_path(members: &[MemberRow]) -> RowsByPath {
    let mut per_path: RowsByPath = RowsByPath::new();
    for member in members {
        if member.start < member.end {
            per_path
                .entry(member.path_name.clone())
                .or_default()
                .push((member.start, member.end));
        }
    }
    per_path
}

// ---------------------------------------------------------------------------
// Unit tests (source-named, the pure core).
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(pairs: &[(u64, u64)]) -> Rows {
        pairs.to_vec()
    }

    #[test]
    fn elementary_segments_split_at_every_row_bound() {
        let mut per_path = RowsByPath::new();
        per_path.insert("p".to_string(), rows(&[(0, 10), (5, 20)]));
        let segments = elementary_segments(&per_path);
        assert_eq!(
            segments,
            vec![
                ("p".to_string(), 0, 5),
                ("p".to_string(), 5, 10),
                ("p".to_string(), 10, 20),
            ]
        );
        // Every row contains every segment it touches.
        for &(_, start, end) in &segments {
            let in_first = end <= 10;
            let in_second = start >= 5;
            assert!(in_first || in_second);
        }
    }

    #[test]
    fn elementary_segments_exclude_disjoint_row_gaps() {
        // The gap between two disjoint rows is not material: no
        // segment is emitted there (0 copies by construction).
        let mut per_path = RowsByPath::new();
        per_path.insert("p".to_string(), rows(&[(0, 10), (20, 30)]));
        let segments = elementary_segments(&per_path);
        assert_eq!(
            segments,
            vec![("p".to_string(), 0, 10), ("p".to_string(), 20, 30)]
        );
        // Nested rows split at every bound, covered pairs only.
        let mut per_path = RowsByPath::new();
        per_path.insert("p".to_string(), rows(&[(0, 10), (20, 30), (5, 25)]));
        let segments = elementary_segments(&per_path);
        assert_eq!(
            segments,
            vec![
                ("p".to_string(), 0, 5),
                ("p".to_string(), 5, 10),
                ("p".to_string(), 10, 20),
                ("p".to_string(), 20, 25),
                ("p".to_string(), 25, 30),
            ]
        );
    }

    #[test]
    fn repeat_visits_count_their_true_multiplicity() {
        // Two rows of one fold over the same extent: a molecule
        // traversing the segment twice carries multiplicity 2.
        let repeat = rows(&[(0, 10), (0, 10)]);
        assert_eq!(containing_rows(&repeat, 0, 5), 2);
        let mut first = RowsByPath::new();
        first.insert("p".to_string(), rows(&[(0, 10), (0, 10)]));
        let mut second = RowsByPath::new();
        second.insert("p".to_string(), rows(&[(0, 10)]));
        let view = locus_dosage_view([&first, &second]);
        assert_eq!(view.segment_count, 1);
        assert_eq!(view.histogram.get(&3), Some(&1));
        assert!(!view.flat_class_of_2);
    }

    #[test]
    fn flat_class_of_2_is_full_double_coverage() {
        let mut same = RowsByPath::new();
        same.insert("p".to_string(), rows(&[(0, 10)]));
        let mut other = RowsByPath::new();
        other.insert("q".to_string(), rows(&[(0, 10)]));
        // Both folds covering every segment: the flat class of 2.
        let view = locus_dosage_view([&same, &same]);
        assert!(view.flat_class_of_2);
        assert_eq!(view.histogram.get(&2), Some(&1));
        // Distinct material: mixed dosage (the heterozygous shape).
        let view = locus_dosage_view([&same, &other]);
        assert!(!view.flat_class_of_2);
        assert_eq!(view.histogram.get(&1), Some(&2));
    }

    #[test]
    fn expected_mass_scales_with_copies_and_length() {
        let mu = 2.0;
        assert_eq!(expected_mass(mu, 1, 100), 200.0);
        assert_eq!(expected_mass(mu, 2, 100), 400.0);
        assert_eq!(expected_mass(mu, 3, 50), 300.0);
    }

    #[test]
    fn merged_intervals_union_touching_windows() {
        let mut spans = vec![(0, 5), (5, 10), (20, 25), (22, 30)];
        merge_intervals(&mut spans);
        assert_eq!(spans, vec![(0, 10), (20, 30)]);
    }

    #[test]
    fn partition_point_finds_the_sorted_range() {
        let sorted = vec![0u64, 10, 20, 30];
        assert_eq!(partition_point(&sorted, |&s| s < 0), 0);
        assert_eq!(partition_point(&sorted, |&s| s < 10), 1);
        assert_eq!(partition_point(&sorted, |&s| s < 25), 3);
        assert_eq!(partition_point(&sorted, |&s| s < 100), 4);
    }
}
