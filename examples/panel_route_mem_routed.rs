//! MEM-routed evidence prototype (experimental, slice-only).
//!
//! Owner-decided redesign of the evidence layer: route read-MEMs through the
//! pangenome instead of simulating reads over candidates.
//!
//! - Observed side: one streaming pass over the sample's read-MEM record
//!   multiset routes each record to the partitions whose pangenome territory
//!   it touches; a record touching k partitions contributes its read
//!   multiplicity as an equal share w/k to each (marginalization prior),
//!   replacing pooled orbit counts. (The stored WeightedBwt is count-only and
//!   exposes no record enumeration without changing source-bound compiler
//!   identities, so the record multiset is re-derived from the base reads
//!   with the public accessor and verified against the index's stats.)
//! - Predicted side: geometric q — a candidate's features are its contained
//!   panel-path anchor subwalks with window-incidence counts from interval
//!   arithmetic against the L150 window grid; no per-window MEM simulation,
//!   no profile cache, no disk artifacts beyond this run's output.
//!
//! The slow path (`panel_route_genome_smoke`) is untouched and remains the
//! comparison oracle; this example reuses its module tree read-only.
#![recursion_limit = "512"]
#[path = "panel_route_diploid_search/mod.rs"]
mod search;
#[path = "panel_route_spine/mod.rs"]
mod spine;

use clap::Parser;
use impg::{
    genome_inference::{
        mem_records, panel_routes as routes, read_json, sample, PanelIdentity, COUNT_POLICY,
    },
    sample_mem_bwt::{canonical, encode_walk, reverse_complement},
    syng::{SyncmerParams, SyngIndex},
};
use rayon::prelude::*;
use search::genome_wide as genome;
use search::partition::{FeatureKey, Profile, ScoreModel, SourceRange};
use serde::Serialize;
use std::sync::atomic::{AtomicU64, Ordering};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fs::File,
    io::{self, BufRead, BufReader},
    path::PathBuf,
    time::Instant,
};

const READ_LENGTH: usize = 150;
const SLACK: u64 = 150;
const MAX_FEATURES: usize = genome::MAX_FEATURES;

#[derive(Parser)]
struct Options {
    #[arg(long)]
    panel: String,
    #[arg(long)]
    routes: PathBuf,
    #[arg(long)]
    sample: PathBuf,
    /// Base reads the sample index was built from (record re-derivation).
    #[arg(long)]
    reads: PathBuf,
    #[arg(long)]
    axis: PathBuf,
    #[arg(long)]
    bed_directory: PathBuf,
    #[arg(long, default_value = "S288C#0#chrIII")]
    component: String,
    /// Loci [START, END) of the loaded component axis (full-component
    /// numbering), matching the oracle slice. Omitted: the full component
    /// (every partition) — resolved after the component partition load.
    #[arg(long, num_args = 2, value_names = ["START", "END"])]
    locus_range: Option<Vec<usize>>,
    #[arg(long, default_value_t = 10.0)]
    depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
    /// Partition universe for support distribution: `component` (the
    /// inferred component's axis partitions) or `genome` (every axis
    /// partition; the honest marginalization, costs a full-territory index).
    #[arg(long, default_value = "component")]
    routing_universe: String,
    /// Window-domain extension (owner-ruled Policy A + minimal universe
    /// extension): each window's candidate rows add every component-family
    /// row from any other group whose interval numerically overlaps the
    /// window's coordinate range, charged against the row's OWNING group's
    /// routed share; pure-new groups enter the COMPONENT routing universe as
    /// appended partitions (territory = exactly their window-overlapping
    /// component-family rows). The genome universe is not extended.
    #[arg(long)]
    window_domain_extension: bool,
    /// Run only the routing stage (inputs -> domain -> records -> territory
    /// index -> routing -> share accumulation -> genome-universe placement
    /// pass), print the routing summary JSON (the reconciliation gate's
    /// share-diff and t_r_genome instrumentation) and exit.
    #[arg(long)]
    routing_diagnostics_only: bool,
    /// Full-component loci (comma-separated) whose forward candidates enter
    /// the geometric-vs-oracle equivalence sample.
    #[arg(long, default_value = "12,19,26")]
    equivalence_loci: String,
    /// Forward candidates per equivalence locus (plus reverse and exit
    /// segment samples).
    #[arg(long, default_value_t = 12)]
    equivalence_per_locus: usize,
    /// Skip the exit-cut test and the equivalence sample (validation stages
    /// already measured); run only the routed-evidence DP.
    #[arg(long)]
    dp_only: bool,
    /// LOCAL-FIRST SPINE (owner redesign): exhaustive local sweep + exact
    /// chain DP + targeted splits + authoritative rescore. The beam path
    /// (below) is quarantined and stays untouched.
    #[arg(long)]
    spine: bool,
    /// With --spine: stop after Stage 1 (the exhaustive local sweep) and
    /// report the per-locus winner table plus the model-vs-search verdict.
    #[arg(long)]
    spine_sweep_only: bool,
    /// With --spine: replace the port-word-viability chain DP with the
    /// correlation-phasing chain DP (panel co-occurrence + boundary MEM
    /// evidence as the transitions; no port-word legality; both homolog
    /// matchings carried as ordered states). Option A partial rows are
    /// generated at structurally detected continuity failures.
    #[arg(long, requires = "spine")]
    spine_phasing: bool,
    /// Standalone falsification-census mode: a JSON file of junction
    /// descriptors {"left": {source,start,end,reverse}, "right": {...},
    /// optional "partitions": [p,...]} whose spanning-read counts (and, when
    /// partitions are given, restricted charges) are reported under the
    /// within-read adjacency index, then exit. No DP, no scoring changes.
    #[arg(long)]
    spine_census_junctions: Option<PathBuf>,
    /// With --spine: a previously selected chain's route JSON (the rescore
    /// `selected_route` format: per-locus [slot-0 {identity, segments},
    /// "empty-slot2"]) whose slot-0 rows are decomposed per-feature against
    /// the truth reference's copy-0 (mosaic) pieces in the HAPLOID
    /// single-allele local table — the Step-A measurement (no model change;
    /// the decomposition reuses the run's own routed shares, measured
    /// backgrounds and oracle profiles, the rescore convention).
    #[arg(long, requires = "spine")]
    decompose_selected: Option<PathBuf>,
    /// With --decompose-selected: output path for the per-feature
    /// decomposition table (compact JSON).
    #[arg(long, requires = "decompose_selected")]
    decompose_out: Option<PathBuf>,
    /// Truth-free bounded split refinement evidence (previous pass JSON),
    /// matching the oracle slow path's admission inputs. When present the
    /// routed-evidence DP runs.
    #[arg(long)]
    refinement_evidence: Option<PathBuf>,
    /// Split prefilter ranks checkpoint (previous pass JSON); restricts split
    /// admission exactly like the oracle slow path.
    #[arg(long)]
    refinement_ranks: Option<PathBuf>,
    #[arg(long, default_value_t = 100.0)]
    refinement_margin: f64,
    #[arg(long, default_value_t = 8)]
    refinement_top_k: usize,
    #[arg(long, default_value_t = 4)]
    refinement_cut_top_n: usize,
    #[arg(long, default_value_t = 256)]
    beam_width: usize,
    #[arg(long, default_value_t = 1)]
    tie_rotations: usize,
    /// Distinct finalists retained for external rescoring per tie rotation.
    #[arg(long, default_value_t = 16)]
    finalists: usize,
    /// Assessment-side reference pair to rescore under the new model:
    /// LABEL:ROUTE_A.json:ROUTE_B.json (repeatable; truth used only for
    /// validation).
    #[arg(long = "reference-pair")]
    reference_pairs: Vec<String>,
    /// With --spine --spine-phasing: a previously measured chain's route JSON
    /// (the rescore `selected_route` format: per-locus [slot-0 {segments},
    /// "empty-slot2"]) mapped into the finalist tables. Its surrogate/L totals
    /// and its MODEL-OF-RECORD total (the re-rank evaluator) are reported:
    /// the k-best finalist re-ranking's bit-identity anchors and the
    /// surrogate-rank measurement.
    #[arg(long = "probe-route")]
    probe_routes: Vec<String>,
    /// In-process RSS guard budget in GiB: the run fails closed with a
    /// stage diagnostic instead of a kernel kill if resident memory exceeds
    /// it. 0 disables the guard.
    #[arg(long, default_value_t = 64.0)]
    rss_budget_gib: f64,
}

fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s)
}

/// Sample the guard at a stage boundary and echo the observed resident peak
/// to stderr, so every run's log carries the per-stage memory curve. Fails
/// closed (Err) when the budget is exceeded.
fn rss_probe(guard: &mut genome::PeriodicRssGuard, stage: &str) -> io::Result<()> {
    guard.sample(stage)?;
    eprintln!("[rss] {stage}: resident-peak {} MiB", guard.peak_bytes() >> 20);
    Ok(())
}

fn ensure(ok: bool, s: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(s))
    }
}

/// Fractional-observed variant of `ScoreModel::loss`, mirroring its exact
/// f64 arithmetic (bit-identical when `observed` is integral).
fn loss_fractional(model: &ScoreModel, q: u64, observed: f64) -> io::Result<f64> {
    let product = q
        .checked_mul(model.histogram)
        .ok_or_else(|| invalid("partition exposure product overflow"))?;
    ensure(product <= 1 << 53, "partition exposure conversion precision limit")?;
    let signal = product as f64 * model.depth / model.denominator;
    let loss = signal - observed * (signal / model.background).ln_1p();
    ensure(loss.is_finite(), "nonfinite fractional-observed loss")?;
    Ok(loss)
}

// ---------------------------------------------------------------------------
// Per-feature multiplicity backgrounds/attribution from the panel (the
// (a)/(b)/(d) family; the earlier panel-scale beta subtraction (a) was
// REJECTED on measurement — see genome/multiplicity-background.md).
// The PRIMARY variant (owner-promoted 2026-09-23, option (b), also the
// env-unset code default since that promotion) raises the feature's
// background to beta_f = base + (m_f - 1) * per_count with the credit
// unchanged. The DIAGNOSTIC variant (owner option (d) + domain hygiene
// (c)) instead attributes the Poisson credit of a feature's observed
// support across its panel contexts — an observation of a feature present
// in m_f equally-likely contexts is attributable to THIS context with
// probability 1/m_f, so the loss's credit term carries the attribution
// weight w_f = 1/m_f (the signal term and the flat background stay
// unchanged). A unique feature (m = 1) and a panel-absent feature (m = 0)
// keep full credit; a feature in 232 contexts earns 1/232 per context.
// No tuning constants: m_f is measured by one panel scan (the independent
// census cross-check validated it 16/16), w_f = 1/max(m_f, 1) is fully
// derived.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MultiplicityVariant {
    /// PRIMARY (owner-promoted 2026-09-23, option (b)): the cross-support
    /// subtraction at the corrected per-count unit —
    /// beta_f = base + (m_f - 1) * per_count, credit unchanged.
    GentleBeta,
    /// Diagnostic (owner option (d)): credit attribution weight w_f =
    /// 1/m_f; measured to leave real pairs positive and flip
    /// truth-vs-native at two tract loci.
    Attribution,
}

/// The code default follows the owner's promotion of (b) to primary
/// (2026-09-23): an unset IMPG_MULTIPLICITY_VARIANT now selects GentleBeta
/// — the env-unset default contradicting the promotion is what made the p1o
/// run silently execute under (d). (d) Attribution stays explicitly
/// selectable with IMPG_MULTIPLICITY_VARIANT=attribution.
pub(crate) fn multiplicity_variant() -> MultiplicityVariant {
    static VARIANT: std::sync::OnceLock<MultiplicityVariant> = std::sync::OnceLock::new();
    VARIANT
        .get_or_init(|| {
            match std::env::var("IMPG_MULTIPLICITY_VARIANT").as_deref() {
                Ok("attribution") => MultiplicityVariant::Attribution,
                _ => MultiplicityVariant::GentleBeta,
            }
        })
        .clone()
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct FeatureBackgroundEntry {
    /// Distinct panel paths whose forward walk contains the feature.
    pub(crate) paths: u32,
    /// Sum over those paths of the feature's per-path L150-window incidence.
    pub(crate) incidence: u64,
    /// Attribution weight w_f = 1/max(paths, 1) (owner option (d), diagnostic).
    pub(crate) weight: f64,
    /// The (a)/(b) family's per-feature background
    /// beta_f = base + (paths - 1) * per_count (PRIMARY under (b)).
    pub(crate) beta: f64,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct BackgroundScanStats {
    /// Distinct features requested (already-known features are not rescanned).
    pub(crate) features: u64,
    /// Features with at least one panel-path occurrence.
    pub(crate) matched: u64,
    /// Positional pattern occurrences across all panel paths (both
    /// orientations), the census quantity of `panel_route_feature_occurrences`.
    pub(crate) occurrences: u64,
    /// Paths whose encoded walk was matched (panel size in paths).
    pub(crate) paths: u64,
    pub(crate) wall_seconds: f64,
}

pub(crate) struct FeatureBackgrounds {
    map: HashMap<FeatureKey, FeatureBackgroundEntry>,
    base: f64,
    /// Encoded forward walk tape + anchor bp positions per panel path,
    /// built on first scan and reused by every later extension scan.
    tapes: Vec<Option<(Vec<u64>, Vec<u64>)>>,
    pub(crate) scans: Vec<BackgroundScanStats>,
}

/// The PRIMARY per-feature background (owner-promoted option (b)): the
/// cross-support subtraction at the corrected per-count unit —
/// beta_f = base + (m_f - 1) * per_count (per window-incidence unit).
pub(crate) fn multiplicity_background_beta(base: f64, paths: u32, per_count: f64) -> f64 {
    if paths == 0 {
        return base;
    }
    base + (paths as f64 - 1.0) * per_count
}

/// The (d) diagnostic's attribution weight: an observation of a feature
/// present in m_f equally-likely panel contexts is attributable to this
/// context with probability 1/m_f; unique (m = 1) and panel-absent
/// (m = 0) features keep full credit.
pub(crate) fn multiplicity_attribution_weight(paths: u32) -> f64 {
    1.0 / paths.max(1) as f64
}

// ---------------------------------------------------------------------------
// The SAMPLE-SIDE cross-support background for the HAPLOID single-allele
// table — the GENOME-WIDE form (owner decision 2026-09-24, candidate family
// 1, "genome-wide touched set"; measured first — see
// genome/genome-wide-cross-support.md): beta_f = base + Sum_r m_r*(1 -
// 1/t_r_genome), summed over every multiset record r realizing feature f,
// where m_r is the record's genome-wide multiplicity and t_r_genome its
// touched count across ALL components' partitions. The OBSERVED share C
// stays component-scoped exactly as before. Why the wide form: the read
// multiset merges the sample's identical repeat reads genome-wide; under
// the component-scoped complementary share the flooding records' touched
// sets truncate to the component (t_r = 1), so their whole multiplicity
// floods one partition with NO background and a long raw-BED row over the
// shared segment harvests ~-1,150 per repeat feature (measured at locus 21:
// 798/837 observed repeat features had complement 0). The genome-wide
// touched count makes the flood's routed-outside mass the background and
// the repeat credit collapses toward 0 (C ~ Sum m_r vs beta ~ base +
// Sum m_r (1 - 1/t_r_genome)). beta_f is a FEATURE-LEVEL quantity — the
// stated form subtracts no share at the scored partition — so the same
// beta applies at every charge site (locus charge and boundary seam).
// No new tuning constants: base, m_r and t_r are all pre-existing or
// measured per run (zero-cache doctrine). The DIPLOID pair table keeps the
// panel background (its (b) acceptance criteria were established on pairs);
// restricted junction charges stay flat (background-invariant).
// ---------------------------------------------------------------------------

/// The sample-side cross-support backgrounds for the haploid track: the
/// per-feature genome-wide cross-support map (built in the genome-universe
/// placement pass) plus the model's base.
pub(crate) struct SampleSideBackgrounds<'a> {
    cross: &'a HashMap<FeatureKey, f64>,
    base: f64,
}

impl<'a> SampleSideBackgrounds<'a> {
    pub(crate) fn new(cross: &'a HashMap<FeatureKey, f64>, base: f64) -> Self {
        Self { cross, base }
    }

    /// The haploid background of a feature — the same feature-level value at
    /// every charge site (locus or boundary seam): base + the genome-wide
    /// cross support Sum_r m_r*(1 - 1/t_r_genome) over the records
    /// realizing it.
    pub(crate) fn beta(&self, feature: &FeatureKey) -> f64 {
        self.base + self.cross.get(feature).copied().unwrap_or(0.0)
    }
}

/// One record's genome-wide cross-support contribution to each of its
/// features: m_r * (1 - 1/t_r_genome) — the record's multiplicity routed
/// outside each one of its own placements. t = 0 (no genome placement) and
/// t = 1 (the whole multiplicity floods one partition, fully observed
/// there) contribute exactly 0.
pub(crate) fn cross_support_contribution(multiplicity: u64, touched: usize) -> f64 {
    if touched == 0 {
        return 0.0;
    }
    multiplicity as f64 * (1.0 - 1.0 / touched as f64)
}

/// Accumulate one record's genome-wide cross support into `map`: every
/// feature of the record's subwalk list gains the record's contribution.
/// No within-record dedup — the same convention as the pooled/share
/// accumulation (a duplicated subwalk key contributes per occurrence in
/// the list, both there and here).
fn accumulate_cross_support(
    map: &mut HashMap<FeatureKey, f64>,
    tokens: &[u64],
    multiplicity: u64,
    touched: usize,
) {
    let contribution = cross_support_contribution(multiplicity, touched);
    if contribution == 0.0 || tokens.is_empty() {
        return;
    }
    if let Ok(subwalks) = enumerate_subwalks(tokens) {
        for feature in &subwalks {
            *map.entry(feature.clone()).or_default() += contribution;
        }
    }
}

/// The per-feature loss under an explicit background (the (b) form with the
/// sample-side beta): signal - observed * ln(1 + signal / beta).
fn loss_fractional_sample(
    model: &ScoreModel,
    q: u64,
    observed: f64,
    beta: f64,
) -> io::Result<f64> {
    let product = q
        .checked_mul(model.histogram)
        .ok_or_else(|| invalid("partition exposure product overflow"))?;
    ensure(product <= 1 << 53, "partition exposure conversion precision limit")?;
    let signal = product as f64 * model.depth / model.denominator;
    ensure(
        beta.is_finite() && beta > 0.0,
        "sample-side background must be positive",
    )?;
    let loss = signal - observed * (signal / beta).ln_1p();
    ensure(loss.is_finite(), "nonfinite sample-side loss")?;
    Ok(loss)
}

/// The HAPLOID single-allele local charge of one profile against one
/// partition's routed shares under the sample-side backgrounds (the
/// merged-with-empty convention: the empty second slot contributes
/// nothing, so the charge is the profile's own terms in key order).
fn merged_single_loss_sample(
    profile: &Profile,
    obs: &HashMap<FeatureKey, f64>,
    sample: &SampleSideBackgrounds,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = 0.0;
    for (feature, &q) in profile {
        let beta = sample.beta(feature);
        loss += loss_fractional_sample(
            model,
            q,
            obs.get(feature).copied().unwrap_or(0.0),
            beta,
        )?;
    }
    Ok(loss)
}

/// The boundary-seam charge under the sample-side backgrounds (the
/// haploid analogue of `profile_loss_boundary_multiplicity`).
fn profile_loss_boundary_sample(
    profile: &Profile,
    obs_left: &HashMap<FeatureKey, f64>,
    obs_right: &HashMap<FeatureKey, f64>,
    sample: &SampleSideBackgrounds,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = 0.0;
    for (feature, &q) in profile {
        let beta = sample.beta(feature);
        let observed = obs_left.get(feature).copied().unwrap_or(0.0)
            + obs_right.get(feature).copied().unwrap_or(0.0);
        loss += loss_fractional_sample(model, q, observed, beta)?;
    }
    Ok(loss)
}

/// Per-feature multiplicity loss: the flat loss with the feature's measured
/// panel entry applied — the PRIMARY (b) per-feature background
/// beta_f = base + (m_f - 1) * per_count, or the (d) diagnostic's
/// attribution weight w_f on the credit term (signal and flat background
/// unchanged). `None` (never scanned / panel-absent) is the
/// flat loss bit-identically.
fn loss_fractional_entry(
    model: &ScoreModel,
    q: u64,
    observed: f64,
    entry: Option<FeatureBackgroundEntry>,
) -> io::Result<f64> {
    let product = q
        .checked_mul(model.histogram)
        .ok_or_else(|| invalid("partition exposure product overflow"))?;
    ensure(product <= 1 << 53, "partition exposure conversion precision limit")?;
    let signal = product as f64 * model.depth / model.denominator;
    let loss = match entry {
        None => signal - observed * (signal / model.background).ln_1p(),
        Some(entry) => match multiplicity_variant() {
            MultiplicityVariant::Attribution => {
                signal - (observed * entry.weight) * (signal / model.background).ln_1p()
            }
            MultiplicityVariant::GentleBeta => {
                let beta = if entry.beta.is_finite() && entry.beta > 0.0 {
                    entry.beta
                } else {
                    model.background
                };
                signal - observed * (signal / beta).ln_1p()
            }
        },
    };
    ensure(loss.is_finite(), "nonfinite fractional-observed loss")?;
    Ok(loss)
}

/// L150-window incidence of one subwalk occurrence at anchor positions
/// [p_first, p_last] on a `length`-bp path: window starts s contain the
/// subwalk's anchors iff s <= p_first and s >= p_last + k - L, clipped to
/// [0, length - L] — the model's own containment arithmetic
/// (`anchor_contain_tail`).
pub(crate) fn window_incidence(p_first: u64, p_last: u64, k: u64, length: u64) -> u64 {
    if length < READ_LENGTH as u64 {
        return 0;
    }
    let lo = p_last.saturating_add(k).saturating_sub(READ_LENGTH as u64);
    let hi = p_first.min(length - READ_LENGTH as u64);
    if hi >= lo {
        hi - lo + 1
    } else {
        0
    }
}

impl FeatureBackgrounds {
    pub(crate) fn new(base: f64) -> Self {
        Self {
            map: HashMap::new(),
            base,
            tapes: Vec::new(),
            scans: Vec::new(),
        }
    }

    /// The feature's measured multiplicity entry; features never scanned
    /// (and features with no panel occurrence) are absent — the loss then
    /// applies full credit (w = 1) and the flat background.
    pub(crate) fn entry(&self, feature: &FeatureKey) -> Option<FeatureBackgroundEntry> {
        self.map.get(feature).copied()
    }

    pub(crate) fn contains(&self, feature: &FeatureKey) -> bool {
        self.map.contains_key(feature)
    }

    pub(crate) fn len(&self) -> usize {
        self.map.len()
    }

    pub(crate) fn iter_entries(
        &self,
    ) -> impl Iterator<Item = (&FeatureKey, FeatureBackgroundEntry)> + '_ {
        self.map.iter().map(|(key, &entry)| (key, entry))
    }

    /// Measure backgrounds for every feature not already known, by matching
    /// the features (both orientations, canonical keys) against the encoded
    /// forward walks of all panel paths. Per match the occurrence's L150
    /// window incidence is accumulated with the model's own containment
    /// arithmetic (a window start s contains the subwalk's anchors iff
    /// s <= p_first and s >= p_last + k - L, clipped to [0, len - L]).
    pub(crate) fn scan_extend(
        &mut self,
        panel: &SyngIndex,
        features: Vec<FeatureKey>,
        k: u64,
        model: &ScoreModel,
    ) -> io::Result<BackgroundScanStats> {
        let targets: Vec<FeatureKey> = features
            .into_iter()
            .filter(|feature| !self.map.contains_key(feature))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let stats = self.scan(panel, &targets, k, model)?;
        self.scans.push(stats);
        Ok(stats)
    }

    fn scan(
        &mut self,
        panel: &SyngIndex,
        targets: &[FeatureKey],
        k: u64,
        model: &ScoreModel,
    ) -> io::Result<BackgroundScanStats> {
        let started = Instant::now();
        let path_count = panel.name_map.path_to_name.len();
        if self.tapes.len() != path_count {
            self.tapes = vec![None; path_count];
        }
        // Build the encoded walk tape (+ anchor bp positions) for every panel
        // path once; later extension scans reuse the cache.
        let missing: Vec<usize> = (0..path_count)
            .filter(|&path| self.tapes[path].is_none())
            .collect();
        let built: Vec<(usize, (Vec<u64>, Vec<u64>))> = missing
            .par_iter()
            .map(|&path| {
                let start = panel.name_map.path_starts[path]
                    .as_ref()
                    .ok_or_else(|| invalid("panel path lacks a forward start"))?;
                let walk = panel.walk_forward_path(start);
                let tape = encode_walk(&walk)?;
                let positions: Vec<u64> = walk.iter().map(|&(_, bp)| bp).collect();
                Ok((path, (tape, positions)))
            })
            .collect::<io::Result<_>>()?;
        for (path, entry) in built {
            self.tapes[path] = Some(entry);
        }
        // Pattern index: singleton features by node token, longer features by
        // their (node, gap) leading pair; a feature's reverse complement is
        // indexed too unless identical (canonical keys may match a path only
        // in that form).
        let mut patterns: Vec<Vec<u64>> = Vec::with_capacity(2 * targets.len());
        let mut owners: Vec<u32> = Vec::with_capacity(2 * targets.len());
        for (id, feature) in targets.iter().enumerate() {
            let reverse = reverse_complement(feature);
            patterns.push(feature.clone());
            owners.push(id as u32);
            if reverse != *feature {
                patterns.push(reverse);
                owners.push(id as u32);
            }
        }
        let mut singleton_index: HashMap<u64, Vec<u32>> = HashMap::new();
        let mut pair_index: HashMap<(u64, u64), Vec<u32>> = HashMap::new();
        for (id, pattern) in patterns.iter().enumerate() {
            match pattern.len() {
                0 => {}
                1 => singleton_index.entry(pattern[0]).or_default().push(id as u32),
                _ => pair_index
                    .entry((pattern[0], pattern[1]))
                    .or_default()
                    .push(id as u32),
            }
        }
        let occurrences = std::sync::atomic::AtomicU64::new(0);
        // Per feature: (paths containing it, summed per-path window incidence).
        let mut paths_seen = vec![0u32; targets.len()];
        let mut incidence = vec![0u64; targets.len()];
        let mut present = vec![false; targets.len()];
        let per_path: Vec<Vec<(u32, u64, u64)>> = (0..path_count)
            .into_par_iter()
            .map(|path| {
                let (tape, positions) = self.tapes[path].as_ref().expect("tape cached");
                let anchors = positions.len();
                let length = panel.name_map.path_to_length[path];
                let mut local: HashMap<u32, (bool, u64)> = HashMap::new();
                let mut local_occurrences = 0u64;
                for anchor in 0..anchors {
                    let base = 2 * anchor;
                    let node = tape[base];
                    let mut matched: Vec<u32> = Vec::new();
                    if let Some(ids) = singleton_index.get(&node) {
                        matched.extend_from_slice(ids);
                    }
                    if base + 1 < tape.len() {
                        if let Some(ids) = pair_index.get(&(node, tape[base + 1])) {
                            matched.extend_from_slice(ids);
                        }
                    }
                    for pattern in matched {
                        let tokens = &patterns[pattern as usize];
                        if base + tokens.len() > tape.len() {
                            continue;
                        }
                        if tape[base..base + tokens.len()] != tokens[..] {
                            continue;
                        }
                        local_occurrences += 1;
                        let id = owners[pattern as usize];
                        let anchor_count = (tokens.len() + 1) / 2;
                        let p_first = positions[anchor];
                        let p_last = positions[anchor + anchor_count - 1];
                        let windows = window_incidence(p_first, p_last, k, length);
                        let entry = local.entry(id).or_insert((false, 0));
                        entry.0 = true;
                        entry.1 = entry
                            .1
                            .checked_add(windows)
                            .ok_or_else(|| invalid("window incidence overflow"))?;
                    }
                }
                occurrences.fetch_add(local_occurrences, Ordering::Relaxed);
                Ok::<_, io::Error>(
                    local
                        .into_iter()
                        .map(|(id, (seen, inc))| (id, u64::from(seen), inc))
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<io::Result<Vec<_>>>()?;
        for local in &per_path {
            for &(id, seen, inc) in local {
                paths_seen[id as usize] += seen as u32;
                incidence[id as usize] = incidence[id as usize]
                    .checked_add(inc)
                    .ok_or_else(|| invalid("panel incidence overflow"))?;
                present[id as usize] = true;
            }
        }
        let per_count = model.histogram as f64 * model.depth / model.denominator;
        let mut matched = 0u64;
        for (id, feature) in targets.iter().enumerate() {
            let paths = if present[id] {
                matched += 1;
                let paths = paths_seen[id];
                ensure(paths >= 1, "presence implies at least one path")?;
                paths
            } else {
                0
            };
            let entry = FeatureBackgroundEntry {
                paths,
                incidence: if paths > 0 { incidence[id] } else { 0 },
                weight: multiplicity_attribution_weight(paths),
                beta: multiplicity_background_beta(self.base, paths, per_count),
            };
            self.map.insert(feature.clone(), entry);
        }
        Ok(BackgroundScanStats {
            features: targets.len() as u64,
            matched,
            occurrences: occurrences.load(Ordering::Relaxed),
            paths: path_count as u64,
            wall_seconds: started.elapsed().as_secs_f64(),
        })
    }
}

// ---------------------------------------------------------------------------
// Territory: the partition universe's pangenome territory and its step index.
// ---------------------------------------------------------------------------

struct UniverseRow {
    /// Universe partition id (axis-order index within the universe scope).
    partition: u32,
    path_name: String,
    start: u64,
    end: u64,
}

struct Territory {
    partition: u32,
    path_idx: usize,
    #[allow(dead_code)]
    start: u64,
    #[allow(dead_code)]
    end: u64,
    /// Padded syncmer steps (bp, node), sorted by bp, covering
    /// [start - SLACK, end + SLACK) intersected with the path.
    steps: Vec<(u64, i32)>,
}

const NODE_KEY_OFFSET: i64 = 1 << 30;

/// Bound-prune tolerance, mirroring the oracle DP's BOUND_PRUNE_EPSILON: the
/// beam score is an f64 fold while the bound is a direct sum, so near-tie
/// states are kept to avoid pruning an optimum by summation noise.
const BOUND_PRUNE_EPSILON: f64 = 1e-6;

struct TerritoryIndex {
    territories: Vec<Territory>,
    /// (node, territory index, bp) sorted by node; entries only for anchors
    /// whose syncmer span overlaps the territory interval.
    entries: Vec<(i32, u32, u64)>,
    /// node_ranges[key] = number of entries with node key < key, so the
    /// range of node with key `k` is [node_ranges[k], node_ranges[k+1]).
    node_ranges: Vec<u64>,
    /// Per path: sorted (bp_start, bp_end, partition) interval lists.
    path_intervals: Vec<Vec<(u64, u64, u32)>>,
}

fn node_key(node: i32) -> usize {
    (node as i64 + NODE_KEY_OFFSET) as usize
}

impl TerritoryIndex {
    fn node_range(&self, node: i32) -> (usize, usize) {
        let key = node_key(node);
        if key + 1 >= self.node_ranges.len() {
            return (0, 0);
        }
        let start = self.node_ranges[key] as usize;
        let end = self.node_ranges[key + 1] as usize;
        (start, end)
    }

    /// Partitions whose territory interval on `path` is overlapped by the
    /// anchor span [bp, bp + k).
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

/// Universe rows plus the component-locus -> universe-partition mapping
/// (identity for the component universe; global axis index for genome).
struct Universe {
    rows: Vec<UniverseRow>,
    component_locus_to_partition: Vec<u32>,
    partitions: usize,
}

fn load_universe(
    axis: &genome::AxisFile,
    bed_directory: &PathBuf,
    lanes: &[(String, u64)],
    component: &str,
    universe: &str,
    component_loci: usize,
    // The window-domain extension (owner-ruled minimal universe extension):
    // when present, the COMPONENT universe appends one partition per
    // pure-new contributing group (sorted-name order, ids component_loci +
    // ordinal) whose territory is EXACTLY the group's window-overlapping
    // component-family rows. The GENOME universe is never extended — the
    // genome-universe placement pass of the genome-wide background stays
    // exactly as is (gap records keep t_r_genome = 0, the ruled measured
    // boundary).
    extension: Option<&genome::WindowDomainExtension>,
) -> io::Result<Universe> {
    let scoped: Vec<(usize, &genome::AxisInterval)> = match universe {
        "component" => axis
            .intervals
            .iter()
            .enumerate()
            .filter(|(_, i)| i.component == component)
            .collect(),
        "genome" => axis.intervals.iter().enumerate().collect(),
        _ => return Err(invalid("routing universe must be component or genome")),
    };
    let mut component_locus_to_partition = vec![u32::MAX; component_loci];
    let mut rows = Vec::new();
    for (partition, (_axis_index, interval)) in scoped.iter().enumerate() {
        if interval.component == component {
            let within = component_locus_to_partition
                .iter()
                .filter(|&&p| p != u32::MAX)
                .count();
            if within < component_loci {
                component_locus_to_partition[within] = partition as u32;
            }
        }
        let bed = bed_directory.join(format!("{}.bed", interval.group));
        for line in BufReader::new(File::open(&bed)?).lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            ensure(fields.len() == 3, "invalid public BED3 group")?;
            let start: u64 = fields[1].parse().map_err(|_| invalid("BED start"))?;
            let end: u64 = fields[2].parse().map_err(|_| invalid("BED end"))?;
            ensure(start < end, "invalid BED interval")?;
            rows.push(UniverseRow {
                partition: partition as u32,
                path_name: fields[0].to_string(),
                start,
                end,
            });
        }
    }
    let mut partitions = scoped.len();
    if let Some(extension) = extension {
        ensure(
            universe == "component",
            "the window-domain extension requires the component routing universe",
        )?;
        ensure(
            extension.component_loci == component_loci,
            "window-domain extension locus cardinality mismatch",
        )?;
        for (ordinal, (_group, extension_rows)) in extension.pure_new.iter().enumerate() {
            let partition = (partitions + ordinal) as u32;
            for &(source, start, end) in extension_rows {
                rows.push(UniverseRow {
                    partition,
                    path_name: lanes[source].0.clone(),
                    start,
                    end,
                });
            }
        }
        partitions += extension.pure_new.len();
    }
    ensure(!rows.is_empty(), "empty routing universe")?;
    ensure(
        component_locus_to_partition.iter().all(|&p| p != u32::MAX),
        "component loci missing from the routing universe",
    )?;
    let partitions = scoped.len();
    Ok(Universe {
        rows,
        component_locus_to_partition,
        partitions,
    })
}

/// The rc frame's step for one matched rc-view syncmer (the rc-symmetric
/// index's unit): the rc k-mer at rc position `q` covers the path's forward
/// [b, b+k) with b = range_lo + range_len - k - q, and its step node is
/// -signed_hash = the path's own forward-frame node id for the k-mer it
/// carries at b (the node identity is frame-independent; only the
/// qualification is frame-dependent).
fn rc_frame_step(signed_hash: i32, q: u64, range_lo: u64, range_len: u64, k: u64) -> (u64, i32) {
    (range_lo + range_len - k - q, -signed_hash)
}

fn build_territory_index(
    panel: &SyngIndex,
    universe: &Universe,
    path_of_name: &HashMap<String, usize>,
    syncmer_len: u64,
    fetch_path_seq: &(dyn Fn(usize, u64, u64) -> io::Result<Vec<u8>> + Sync),
    both_frames: bool,
) -> io::Result<TerritoryIndex> {
    let rows = &universe.rows;
    let path_count = panel.name_map.path_to_name.len();
    let territories: Vec<Territory> = rows
        .par_iter()
        .map(|row| {
            let path_idx = *path_of_name.get(&row.path_name).ok_or_else(|| {
                invalid(&format!("universe path {} absent from panel", row.path_name))
            })?;
            let lo = row.start.saturating_sub(SLACK);
            let hi = row.end + SLACK;
            let mut steps: Vec<(u64, i32)> = panel
                .walk_path_range(path_idx, lo, hi)?
                .into_iter()
                .map(|(node, bp)| (bp, node))
                .collect();
            // The rc-symmetric index (the strand-asymmetry fix): the path's
            // step list above carries only the forward frame's qualifying
            // k-mers — the syncmer scheme is not rc-symmetric, so a locus
            // where the path stores a k-mer's rc (the rc strand's k-mer
            // qualifies, the forward strand's does not) has NO step at all,
            // and every record whose occurrence there is an rc occurrence is
            // unplaceable (the route_record rc search's anchors are exactly
            // the forward hash's negation, which the verify can never meet).
            // The fix emits the RC frame's qualifying k-mers as steps too:
            // run the same matched-syncmer extraction on the range's reverse
            // complement; its k-mer at rc position q covers the path's
            // forward [b, b+k) with b = len - k - q, its signed hash m =
            // -hash(path[b..b+k)), and the step node -m = hash(path[b..b+k))
            // — the SAME node id the forward frame would assign (the node
            // identity is frame-independent; only the qualification is
            // frame-dependent). With both frames' steps, a record's rc-
            // orientation search verifies at rc-stored loci (its anchors'
            // k-mers qualify on the rc strand there by construction — they
            // are the record's own selected anchors), and the rc-frame
            // records' forward search verifies symmetrically. No query-side
            // change, no record change, no remap change.
            if both_frames {
                let seq = fetch_path_seq(path_idx, lo, hi)?;
                let seq_len = seq.len() as u64;
                let rc_seq = impg::graph::reverse_complement(&seq);
                // The RAW single-frame extraction (no best-orientation
                // selection — the selector discards the losing frame, and
                // the losing frame is exactly where the rc-frame-only
                // qualifying k-mers live).
                for (signed_node, q) in mem_records::raw_matched_syncmers(panel, &rc_seq)? {
                    steps.push(rc_frame_step(signed_node, q, lo, seq_len, syncmer_len));
                }
            }
            steps.sort_by_key(|&(bp, _)| bp);
            steps.dedup_by(|a, b| a.0 == b.0 && a.1 == b.1);
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
            if bp + syncmer_len > t.start && bp < t.end {
                entries.push((node, index as u32, bp));
            }
        }
    }
    entries.par_sort_unstable_by_key(|&(node, _, _)| node);
    let max_node = panel.num_syncmer_nodes() as i64;
    ensure(
        NODE_KEY_OFFSET > max_node,
        "node key offset collides with node space",
    )?;
    let table_len = node_key(max_node as i32) + 2;
    let mut counts = vec![0u64; table_len];
    for &(node, _, _) in &entries {
        counts[node_key(node)] += 1;
    }
    let mut node_ranges = vec![0u64; table_len];
    let mut acc = 0u64;
    for key in 0..table_len {
        node_ranges[key] = acc;
        acc += counts[key];
    }
    Ok(TerritoryIndex {
        territories,
        entries,
        node_ranges,
        path_intervals,
    })
}

// ---------------------------------------------------------------------------
// Base reads streaming (sample record re-derivation and the within-read
// adjacency chains both consume it).
// ---------------------------------------------------------------------------

fn stream_fastq(path: &PathBuf, mut visit: impl FnMut(&[u8]) -> io::Result<()>) -> io::Result<()> {
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
        visit(seq.as_bytes())?;
        match lines.next().transpose()? {
            Some(line) => current = line,
            None => break,
        }
    }
    Ok(())
}


// ---------------------------------------------------------------------------
// Routing: record occurrences -> touched partitions -> shares.
// ---------------------------------------------------------------------------

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
    // The rc walk's anchor for the k-mer [rel, rel+k) is the mirrored
    // k-mer's START in the rc frame: span - k - rel (the rc of the walk's
    // sequence reverses the k-mer order AND mirrors each k-mer's start;
    // the previous span - rel form was the mirrored k-mer's END — off by
    // k — which made every multi-anchor rc verification fail by geometry).
    for &(node, pos) in anchors.iter().rev() {
        out.push((-node, span - k - pos));
    }
    out
}

/// Verified occurrences of one orientation of a record, deduplicated by
/// (path, occurrence start). Search anchors: either the least-indexed anchor
/// (default; can miss occurrences whose search anchor falls outside every
/// territory interval) or every anchor (exact diagnostic mode).
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

/// partition -> deduplicated partition-touch count per verified occurrence,
/// plus the total number of verified occurrences (both orientations), for one
/// record. Touch semantics: an occurrence touches partition p when any of its
/// anchor spans overlaps one of p's territory intervals on the occurrence's
/// path; one occurrence counts at most once per partition.
fn route_record(
    anchors: &[(i32, u64)],
    index: &TerritoryIndex,
    k: u64,
    every_anchor: bool,
) -> (
    BTreeMap<u32, u64>,
    u64,
    Vec<(usize, u64)>,
    Vec<(usize, u64)>,
) {
    let mut occurrences: BTreeMap<u32, u64> = BTreeMap::new();
    let mut total_occurrences = 0u64;
    // The FORWARD orientation's verified occurrence positions (path index,
    // occurrence start), sorted — the placement positions the in-window
    // feature attribution (the omission charge's universes) attributes by.
    let mut forward_positions: Vec<(usize, u64)> = Vec::new();
    // The RC orientation's verified occurrence positions (the same form as
    // the forward positions; the strand-asymmetry fix: the placement sets
    // and every downstream consumer must see both orientations' real
    // occurrences).
    let mut reverse_positions: Vec<(usize, u64)> = Vec::new();
    for (orientation_index, orientation) in [anchors, &reverse_complement_walk(anchors, k)]
        .into_iter()
        .enumerate()
    {
        let mut seen = HashSet::new();
        if every_anchor {
            for anchor in 0..orientation.len() {
                occurrences_from_anchor(orientation, anchor, index, &mut seen);
            }
        } else {
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
        }
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
    (occurrences, total_occurrences, forward_positions, reverse_positions)
}

// ---------------------------------------------------------------------------
// ROUTER-PLACEMENT DIAGNOSIS (env-gated: IMPG_ROUTER_PLACEMENT_DIAG=spec.json).
//
// The instance-exoneration's measured premise failure (the flood records'
// global placement sets omit the spelled true home) is diagnosed at the
// placement construction itself. route_record searches verified full-walk
// occurrences from ONE anchor (the least-populated node in the territory
// index) and keeps only full-walk-verified occurrences inside routed
// territory intervals. This diagnostic (a) reproduces each targeted record's
// production placements, (b) re-searches with EVERY anchor (the exact mode),
// and (c) force-verifies the record's walk at the spelled loci's territories
// — separating a SEARCH-SPACE artifact (the walk verifies at the locus but
// the single-anchor search never looked) from a VERIFICATION failure (the
// read genuinely does not occur there) from a TERRITORY gap (the locus lies
// outside every routed territory interval). Writes its report and exits
// before any downstream stage; unset env = zero cost.
// ---------------------------------------------------------------------------

#[derive(serde::Deserialize)]
struct RouterDiagSpec {
    /// The targeted records' known feature spans [path, lo, hi] (from the
    /// probe report's instance_dump) identifying the flood records.
    target_spans: Vec<[u64; 3]>,
    /// The spelled true-home loci [path, lo, hi] to force-verify against.
    target_loci: Vec<[u64; 3]>,
    #[serde(default = "router_diag_default_records")]
    max_target_records: usize,
    #[serde(default = "router_diag_default_stride")]
    census_stride: usize,
    #[serde(default = "router_diag_default_cap")]
    census_cap: usize,
    /// Anchors whose node's territory-index entry count exceeds this are
    /// skipped in the census (their full scan is the production search's own
    /// cost multiplied by the entry count); the skips are reported.
    #[serde(default = "router_diag_default_entry_cap")]
    census_anchor_entry_cap: usize,
    out: String,
}

fn router_diag_default_records() -> usize {
    8
}
fn router_diag_default_stride() -> usize {
    8
}
fn router_diag_default_cap() -> usize {
    2048
}
fn router_diag_default_entry_cap() -> usize {
    30_000
}

/// (matched, mismatched, missing) walk steps at `start` inside `territory` —
/// verify_walk's own test, instrumented.
fn walk_match_detail(walk: &[(i32, u64)], start: u64, territory: &Territory) -> (usize, usize, usize) {
    let steps = &territory.steps;
    let mut matched = 0usize;
    let mut mismatched = 0usize;
    let mut missing = 0usize;
    for &(node, rel) in walk {
        let bp = match start.checked_add(rel) {
            Some(bp) => bp,
            None => {
                missing += 1;
                continue;
            }
        };
        let mut lo = 0usize;
        let mut hi = steps.len();
        let mut found = None;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            match steps[mid].0.cmp(&bp) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    found = Some(steps[mid].1);
                    break;
                }
            }
        }
        match found {
            None => missing += 1,
            Some(n) if n == node => matched += 1,
            Some(_) => mismatched += 1,
        }
    }
    (matched, mismatched, missing)
}

fn router_placement_diagnosis(
    spec_path: &str,
    routed_records: &[RoutedRecord],
    territory: &TerritoryIndex,
    panel: &SyngIndex,
    k: u64,
) -> io::Result<()> {
    let spec: RouterDiagSpec = read_json(std::path::Path::new(spec_path))?;
    let targets: HashSet<(u64, u64, u64)> = spec
        .target_spans
        .iter()
        .map(|&span| (span[0], span[1], span[2]))
        .collect();

    // -- census: production (single-anchor) vs every-anchor FORWARD placement
    // sets on a deterministic stride sample; position-level, not partitions.
    let mut census_records = 0usize;
    let mut census_records_with_extra = 0usize;
    let mut census_total_extra = 0usize;
    let mut census_max_extra = 0usize;
    let mut census_heavy_anchor_skips = 0usize;
    let mut census_records_skipped_heavy = 0usize;
    let mut census_extra_on_target_locus_paths = 0usize;
    let target_locus_paths: HashSet<u64> = spec.target_loci.iter().map(|l| l[0]).collect();
    for (index, record) in routed_records.iter().enumerate() {
        if index % spec.census_stride != 0 || census_records >= spec.census_cap {
            continue;
        }
        census_records += 1;
        let Ok(anchors) = decode_tokens(&record.tokens) else {
            continue;
        };
        let mut seen = HashSet::new();
        let mut heavy_anchors = 0usize;
        for anchor in 0..anchors.len() {
            let (start, end) = territory.node_range(anchors[anchor].0);
            if end - start > spec.census_anchor_entry_cap {
                heavy_anchors += 1;
                continue;
            }
            occurrences_from_anchor(&anchors, anchor, territory, &mut seen);
        }
        census_heavy_anchor_skips += heavy_anchors;
        let production: HashSet<(usize, u64)> = record.forward_positions.iter().copied().collect();
        let extra: Vec<&(usize, u64)> = seen
            .iter()
            .filter(|position| !production.contains(*position))
            .collect();
        if !extra.is_empty() {
            census_records_with_extra += 1;
            census_total_extra += extra.len();
            census_max_extra = census_max_extra.max(extra.len());
            if extra.iter().any(|&&(p, _)| target_locus_paths.contains(&(p as u64))) {
                census_extra_on_target_locus_paths += 1;
            }
        }
        if heavy_anchors > 0 {
            census_records_skipped_heavy += 1;
        }
    }

    // -- targeted trace: the flood records' placements and the forced
    // verification at the spelled loci.
    let mut target_reports = Vec::new();
    let mut diag_feature_nodes: HashSet<i32> = HashSet::new();
    for record in routed_records.iter() {
        if target_reports.len() >= spec.max_target_records {
            break;
        }
        let Ok(anchors) = decode_tokens(&record.tokens) else {
            continue;
        };
        // Production placement spans of the record's single anchors (the
        // i==j subwalks) — the shapes the dumped feature spans take.
        let mut is_target = false;
        for &(p, occ) in &record.forward_positions {
            for &(node, rel) in &anchors {
                if targets.contains(&((p as u64), occ + rel, occ + rel + k)) {
                    is_target = true;
                    diag_feature_nodes.insert(node);
                }
            }
            if is_target {
                break;
            }
        }
        if !is_target {
            continue;
        }
        // The exact mode's full forward search.
        let (every_occurrences, every_total, every_forward, _) = route_record(&anchors, territory, k, true);
        let production_set: HashSet<(usize, u64)> = record.forward_positions.iter().copied().collect();
        let mut extra_positions: Vec<[u64; 2]> = every_forward
            .iter()
            .filter(|position| !production_set.contains(*position))
            .take(40)
            .map(|&(p, occ)| [p as u64, occ])
            .collect();
        extra_positions.sort_unstable();

        // The production search's own anchor choice (least-populated node).
        let mut best = 0usize;
        let mut best_entries = usize::MAX;
        for (anchor, &(node, _)) in anchors.iter().enumerate() {
            let (start, end) = territory.node_range(node);
            let count = end - start;
            if count < best_entries {
                best_entries = count;
                best = anchor;
            }
        }
        let chosen_node = anchors[best].0;

        // Forced verification at the spelled loci.
        let mut locus_reports = Vec::new();
        for locus in &spec.target_loci {
            let (locus_path, locus_lo, locus_hi) = (locus[0], locus[1], locus[2]);
            let territories: Vec<&Territory> = territory
                .territories
                .iter()
                .filter(|t| t.path_idx as u64 == locus_path)
                .collect();
            let mut territory_reports = Vec::new();
            let mut any_full_verify = false;
            let mut best_overall = (0usize, 0usize, 0usize);
            // The panel's own steps at the locus (the ground truth of what
            // sequence/nodes live there).
            let panel_steps_near: Vec<[i64; 2]> = panel
                .walk_path_range(locus_path as usize, locus_lo.saturating_sub(100), locus_hi + 100)
                .map(|steps| {
                    steps
                        .iter()
                        .take(64)
                        .map(|&(node, bp)| [node as i64, bp as i64])
                        .collect()
                })
                .unwrap_or_else(|_| Vec::new());
            // Where the flood features' own k-mer nodes actually occur on
            // this path (the spelled span's claimed locus vs the node's real
            // positions — the spelled-coordinate defect's direct measure).
            let feature_nodes_here: Vec<i32> = diag_feature_nodes
                .iter()
                .copied()
                .collect();
            let feature_node_occurrences: Vec<[i64; 2]> = panel
                .walk_path_range(locus_path as usize, 0, u64::MAX)
                .map(|steps| {
                    steps
                        .iter()
                        .filter(|&(node, _)| {
                            diag_feature_nodes.contains(node)
                                || diag_feature_nodes.contains(&-*node)
                        })
                        .map(|&(node, bp)| [node as i64, bp as i64])
                        .take(64)
                        .collect()
                })
                .unwrap_or_else(|_| Vec::new());
            // The routed universe's own homes for these nodes (the territory
            // index's entries — where the k-mer actually verifies inside
            // routed territories), with the panel paths' NAMES — the
            // lane-name -> panel-path binding test.
            let mut index_homes: Vec<[String; 3]> = Vec::new();
            for &node in &feature_nodes_here {
                for sign in [node, -node] {
                    let (start, end) = territory.node_range(sign);
                    for entry in &territory.entries[start..end.min(start + 4096)] {
                        let t = &territory.territories[entry.1 as usize];
                        index_homes.push([
                            panel.name_map.path_to_name[t.path_idx].clone(),
                            entry.2.to_string(),
                            sign.to_string(),
                        ]);
                        if index_homes.len() >= 256 {
                            break;
                        }
                    }
                }
            }
            index_homes.sort_unstable();
            index_homes.dedup();
            let _ = feature_nodes_here;
            let watched_path_names: Vec<String> = spec
                .target_loci
                .iter()
                .map(|locus| panel.name_map.path_to_name[locus[0] as usize].clone())
                .collect();
            for t in &territories {
                let chosen_entries = t
                    .steps
                    .iter()
                    .filter(|&&(bp, node)| node == chosen_node && bp + k > t.start && bp < t.end)
                    .count();
                // Candidate occurrence starts: the record's single-anchor
                // feature span [occ + rel, occ + rel + k) overlapping the
                // locus, for every anchor of the walk.
                let mut best_here = (0usize, 0usize, 0usize);
                let mut full_verify_at: Option<u64> = None;
                for &(_, rel_i) in anchors.iter() {
                    // occ in [lo + 1 - k - rel, hi - 1 - rel] (i64 math).
                    let occ_lo_i = locus_lo as i64 + 1 - k as i64 - rel_i as i64;
                    let occ_hi_i = locus_hi as i64 - 1 - rel_i as i64;
                    if occ_lo_i > occ_hi_i {
                        continue;
                    }
                    let occ_lo = occ_lo_i.max(0) as u64;
                    let occ_hi = occ_hi_i.max(0) as u64;
                    for occ in occ_lo..=occ_hi.min(occ_lo + 511) {
                        let detail = walk_match_detail(&anchors, occ, t);
                        if detail.0 > best_here.0 {
                            best_here = detail;
                        }
                        if detail.0 == anchors.len() && detail.1 == 0 && detail.2 == 0 {
                            full_verify_at = Some(occ);
                            any_full_verify = true;
                            break;
                        }
                    }
                    if full_verify_at.is_some() {
                        break;
                    }
                }
                // Exhaustive phase-aligned verification: every syncmer step
                // bp of the territory is a candidate occurrence start for
                // every anchor offset (the complete candidate set — no
                // search-space assumption at all), in BOTH orientations (the
                // production router keeps only the forward orientation's
                // occurrences in forward_positions).
                let rc_walk = reverse_complement_walk(&anchors, k);
                let mut exhaustive_full_verify_at: Option<u64> = None;
                let mut exhaustive_best = (0usize, 0usize, 0usize);
                let mut exhaustive_orientation = 0usize;
                if t.start <= locus_lo && locus_hi <= t.end {
                    'outer: for (orientation_index, walk) in [&anchors, &rc_walk].into_iter().enumerate() {
                        for &(bp, _) in t.steps.iter() {
                            for &(_, rel_i) in walk.iter() {
                                let Some(occ) = bp.checked_sub(rel_i) else {
                                    continue;
                                };
                                let detail = walk_match_detail(walk, occ, t);
                                if (detail.0, -(detail.1 as i64))
                                    > (exhaustive_best.0, -(exhaustive_best.1 as i64))
                                {
                                    exhaustive_best = detail;
                                    exhaustive_orientation = orientation_index;
                                }
                                if detail.0 == walk.len() && detail.1 == 0 && detail.2 == 0 {
                                    exhaustive_full_verify_at = Some(occ);
                                    any_full_verify = true;
                                    exhaustive_orientation = orientation_index;
                                    break 'outer;
                                }
                            }
                        }
                    }
                }
                if (best_here.0, -(best_here.1 as i64)) > (best_overall.0, -(best_overall.1 as i64)) {
                    best_overall = best_here;
                }
                territory_reports.push(serde_json::json!({
                    "territory_start": t.start,
                    "territory_end": t.end,
                    "covers_locus": t.start <= locus_lo && locus_hi <= t.end,
                    "steps": t.steps.len(),
                    "chosen_anchor_entries_in_interval": chosen_entries,
                    "best_match": best_here,
                    "full_verify_at": full_verify_at,
                    "exhaustive_best_match": exhaustive_best,
                    "exhaustive_full_verify_at": exhaustive_full_verify_at,
                    "exhaustive_orientation": exhaustive_orientation,
                }));
            }
            locus_reports.push(serde_json::json!({
                "locus": [locus_path, locus_lo, locus_hi],
                "territories_on_path": territories.len(),
                "any_full_verify": any_full_verify,
                "best_match_any_territory": best_overall,
                "walk_len": anchors.len(),
                "panel_steps_near_locus": panel_steps_near,
                "feature_node_occurrences_on_locus_path": feature_node_occurrences,
                "feature_node_index_homes": index_homes,
                "watched_path_names": watched_path_names,
                "territories": territory_reports,
            }));
        }
        target_reports.push(serde_json::json!({
            "multiplicity": record.multiplicity,
            "anchors": anchors.len(),
            "walk": anchors.iter().take(24).map(|&(node, rel)| [node as i64, rel as i64]).collect::<Vec<_>>(),
            "anchors": anchors.len(),
            "total_occurrences": record.total_occurrences,
            "production_forward_positions": record.forward_positions.len(),
            "production_reverse_positions": record.reverse_positions.len(),
            "production_reverse_positions_list": record
                .reverse_positions
                .iter()
                .take(40)
                .map(|&(p, occ)| [p as u64, occ])
                .collect::<Vec<_>>(),
            "production_forward_positions_list": record
                .forward_positions
                .iter()
                .take(40)
                .map(|&(p, occ)| [p as u64, occ])
                .collect::<Vec<_>>(),
            "every_anchor_forward_positions": every_forward.len(),
            "every_anchor_total_occurrences": every_total,
            "production_positions_on_locus_paths": spec
                .target_loci
                .iter()
                .map(|locus| {
                    let p = locus[0] as usize;
                    let lo = locus[1];
                    let hi = locus[2];
                    serde_json::json!({
                        "path": p,
                        "positions": record
                            .forward_positions
                            .iter()
                            .filter(|&&(pp, occ)| {
                                pp == p
                                    && occ
                                        + anchors
                                            .iter()
                                            .map(|&(_, rel)| rel)
                                            .min()
                                            .unwrap_or(0)
                                    <= hi
                                    && occ
                                        + anchors
                                            .iter()
                                            .map(|&(_, rel)| rel)
                                            .max()
                                            .unwrap_or(0)
                                            + k
                                        >= lo
                            })
                            .map(|&(pp, occ)| occ as u64)
                            .take(16)
                            .collect::<Vec<_>>(),
                        "reverse_positions": record
                            .reverse_positions
                            .iter()
                            .filter(|&&(pp, occ)| {
                                pp == p
                                    && occ
                                        + anchors
                                            .iter()
                                            .map(|&(_, rel)| rel)
                                            .min()
                                            .unwrap_or(0)
                                    <= hi
                                    && occ
                                        + anchors
                                            .iter()
                                            .map(|&(_, rel)| rel)
                                            .max()
                                            .unwrap_or(0)
                                            + k
                                        >= lo
                            })
                            .map(|&(pp, occ)| occ as u64)
                            .take(16)
                            .collect::<Vec<_>>(),
                    })
                })
                .collect::<Vec<_>>(),
            "every_anchor_extra_positions": extra_positions,
            "every_anchor_partitions": every_occurrences.len(),
            "chosen_anchor": {"index": best, "node": chosen_node, "entries": best_entries},
            "loci": locus_reports,
        }));
    }

    let report = serde_json::json!({
        "census": {
            "stride": spec.census_stride,
            "records": census_records,
            "records_with_extra_positions": census_records_with_extra,
            "total_extra_positions": census_total_extra,
            "max_extra_positions_one_record": census_max_extra,
            "records_with_extra_on_target_locus_paths": census_extra_on_target_locus_paths,
            "heavy_anchor_skips": census_heavy_anchor_skips,
            "records_skipped_heavy": census_records_skipped_heavy,
            "anchor_entry_cap": spec.census_anchor_entry_cap,
        },
        "targets": target_reports,
    });
    let out_file = std::fs::File::create(&spec.out)?;
    serde_json::to_writer_pretty(&out_file, &report)?;
    eprintln!("[router-diag] wrote {} ({} target records)", spec.out, target_reports.len());
    Ok(())
}

/// Per-(partition, feature) routed observed support (fractional shares).
type RoutedObs = HashMap<u32, HashMap<FeatureKey, f64>>;

fn merge_routed(target: &mut RoutedObs, source: &RoutedObs) {
    for (partition, features) in source {
        let into = target.entry(*partition).or_default();
        for (feature, share) in features {
            *into.entry(feature.clone()).or_default() += *share;
        }
    }
}

// ---------------------------------------------------------------------------
// RECORD-ONCE SITE OBSERVED MAPS (owner-approved Fix 1, 2026-09-24 — the
// red-line enforcement of the window-domain extension's mixed-owner
// charging). The extension's union-sum convention (sum over the DISTINCT
// owning partitions' routed-share maps) re-pools record mass: a record
// touching k charged owners contributes m_r/t_r to EACH owner's map, so a
// site that sums the maps attributes the record k times. The fixed form: a
// charge site's observed support for a feature = the SUM over the records
// realizing it of the record's equal-share contribution m_r/t_r, each
// record attributed ONCE PER SITE — a record counts at a site iff it has a
// placement inside the site's charged owners' territory. Fully derived from
// the routing's placement structure (multiplicity, touched count, touched
// partitions); no constants. For a single-owner site the form reduces to
// exactly the existing per-partition routed-share map (same records, same
// order, same values — bit-identical), so single-owner sites keep their
// precomputed maps and the builder is invoked only for owner sets with two
// or more owners (the mixed-owner sites the extension introduced).
// ---------------------------------------------------------------------------

/// One routed sample record: its token walk, genome-wide multiplicity, and
/// the per-universe-partition verified occurrence counts (both orientations,
/// one occurrence counted at most once per partition), plus the total
/// verified occurrence count.
pub(crate) struct RoutedRecord {
    tokens: Vec<u64>,
    multiplicity: u64,
    occurrences: BTreeMap<u32, u64>,
    total_occurrences: u64,
    /// The FORWARD orientation's verified occurrence positions (path index,
    /// occurrence start), sorted (the in-window feature attribution's
    /// placement positions).
    forward_positions: Vec<(usize, u64)>,
    /// The RC orientation's verified occurrence positions (path index,
    /// occurrence start of the RC walk's anchors), sorted — the strand
    /// asymmetry fix: the record's placement sets must see the read's
    /// rc-orientated occurrences too.
    reverse_positions: Vec<(usize, u64)>,
}

/// One sample record's placement structure for the record-once site maps.
pub(crate) struct SiteRecord {
    /// The record's genome-wide multiplicity (the read multiset count).
    pub(crate) multiplicity: u64,
    /// The record's touched-partition count — the established equal-share
    /// denominator t_r of the routed-share accumulation.
    pub(crate) touched: u64,
    /// The record's touched universe partitions (sorted, deduplicated).
    pub(crate) partitions: Vec<u32>,
    /// The record's token walk (its feature subwalks enumerate from it).
    pub(crate) tokens: Vec<u64>,
}

/// Structural accessor trait so `SiteObserved::build` accepts the routing's
/// local record type without hoisting it to module scope.
pub(crate) trait SiteRecordSource {
    fn site_multiplicity(&self) -> u64;
    fn site_touched(&self) -> u64;
    fn site_partitions(&self) -> impl Iterator<Item = &u32>;
    fn site_tokens(&self) -> &Vec<u64>;
}

/// The record-once site-map builder over the slice-touching records (in the
/// share accumulation's order, so a single-owner site's map is the
/// per-partition map bit-identically).
pub(crate) struct SiteObserved {
    records: Vec<SiteRecord>,
    /// Universe partition -> indices into `records` of the records touching
    /// it (accumulation order preserved).
    by_partition: HashMap<u32, Vec<u32>>,
}

impl SiteObserved {
    /// Build from the routed records (accumulation order) restricted to the
    /// slice's partitions.
    pub(crate) fn build(routed: &[RoutedRecord], slice_partitions: &BTreeSet<u32>) -> Self {
        let mut records = Vec::new();
        let mut by_partition: HashMap<u32, Vec<u32>> = HashMap::new();
        for record in routed {
            let partitions: Vec<u32> = record
                .occurrences
                .keys()
                .filter(|&partition| slice_partitions.contains(partition))
                .copied()
                .collect();
            if partitions.is_empty() {
                continue;
            }
            let index = records.len() as u32;
            for &partition in &partitions {
                by_partition.entry(partition).or_default().push(index);
            }
            records.push(SiteRecord {
                multiplicity: record.multiplicity,
                touched: record.occurrences.len() as u64,
                partitions,
                tokens: record.tokens.clone(),
            });
        }
        Self {
            records,
            by_partition,
        }
    }

    /// The record-once observed map for a set of owning UNIVERSE partitions:
    /// every record with a touched partition in the set contributes
    /// m_r/t_r ONCE, per feature of its walk. Deterministic: owners ascend,
    /// records keep accumulation order.
    pub(crate) fn site_map<I: IntoIterator<Item = u32>>(
        &self,
        owners: I,
    ) -> HashMap<FeatureKey, f64> {
        let owners: BTreeSet<u32> = owners.into_iter().collect();
        let mut map: HashMap<FeatureKey, f64> = HashMap::new();
        if owners.is_empty() {
            return map;
        }
        let mut seen = vec![false; self.records.len()];
        for &owner in &owners {
            let indices = match self.by_partition.get(&owner) {
                Some(indices) => indices,
                None => continue,
            };
            for &index in indices {
                if seen[index as usize] {
                    continue;
                }
                seen[index as usize] = true;
                let record = &self.records[index as usize];
                let share = record.multiplicity as f64 / record.touched as f64;
                if let Ok(subwalks) = enumerate_subwalks(&record.tokens) {
                    for feature in subwalks {
                        *map.entry(feature).or_default() += share;
                    }
                }
            }
        }
        map
    }
}

/// The OMISSION CHARGE addend (owner-approved Fix 2, 2026-09-24): the
/// haploid-track candidate's cost for the window's observed DNA it does not
/// predict. For each feature of the WINDOW'S FULL OBSERVED PROFILE (the
/// record-once routed shares over the window's owning partitions) that the
/// row's profile omits, the row pays the pure-background Poisson cost
/// beta_f - C_f*ln(beta_f) — the absolute cost of explaining the feature's
/// observed reads as pure background, growing with the unexplained observed
/// mass. This closes the loophole that made ignoring observed DNA free (the
/// measured pathology: 162bp donor fragments winning tract windows while
/// leaving ~9.8k share of observed mass unexplained). The row's EXPLAINED
/// features keep the current charge bit-exactly (the established relative
/// form over its own owner-resolved observed side), so every existing
/// discrimination is preserved and the omitted-observed-mass cost is new;
/// between two rows omitting the same feature the term cancels, and between
/// an omittor and an explainer the difference gains exactly the feature's
/// omission cost. No constants: beta_f is the charged genome-wide
/// sample-side background and C_f the record-once routed share.
/// ln(Gamma(x + 1)) for x > -1 via the Lanczos approximation (g = 7, the
/// standard nine-coefficient form; numerical infrastructure, not a model
/// constant). The omission charge's Poisson count term needs the count's
/// log-factorial at fractional observed shares (the equal-share routed
/// supports are fractional), which the stable standard library does not
/// expose (`float_gamma` is unstable).
fn ln_gamma_observation(observed: f64) -> f64 {
    const LANCZOS_G: f64 = 7.0;
    const LANCZOS: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_13,
        -176.615_029_162_140_59,
        12.507_343_278_68_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_576e-6,
        1.505_632_735_149_311_6e-7,
    ];
    // ln Gamma(z) = 0.5*ln(2*pi) + (z - 0.5)*ln(t) - t + ln(series),
    // t = z + g - 0.5, for z > 0.5 (here z = observed + 1 > 1). The Lanczos
    // correction series is a SUM OF RECIPROCAL TERMS:
    // series = c0 + sum_{k=1..8} c_k / (z - 1 + k).
    let z = observed + 1.0;
    let mut series = LANCZOS[0];
    for k in 1..LANCZOS.len() {
        series += LANCZOS[k] / (z - 1.0 + k as f64);
    }
    let t = z + LANCZOS_G - 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z - 0.5) * t.ln() - t + series.ln()
}

pub(crate) fn omission_charge_addend(
    profile: &Profile,
    window_obs: &HashMap<FeatureKey, f64>,
    sample: &SampleSideBackgrounds,
) -> f64 {
    let omitted: Vec<&FeatureKey> = window_obs
        .keys()
        .filter(|feature| !profile.contains_key(*feature))
        .collect();
    omission_term_sum(&omitted, window_obs, sample)
}

/// The CHAIN-LEVEL once-per-material form of the omission addend (supervisor
/// ruling, genome/stitching-omission-alignment 2026-09-24): the window's
/// omission covers the in-window observed features the row does not predict
/// AND the chain explains NOWHERE along the whole chain (the union of the
/// charged per-window profiles' features). A feature the chain spells at any
/// window is exonerated here — conserved material charged at its first window
/// no longer pays the full window omission at its later windows (the measured
/// artifact: the group-shared-once filter emptied the truth's charged profile
/// at loci 6/9/11/13 and the truth paid 201,598 of full-window omission for
/// material it explains at its first window). An EMPTY explained set
/// reproduces the per-window form bit-exactly (the DP's surrogate layer keeps
/// that form; the model of record — the references, the selected chain's
/// chain-level M1 and the decompose — uses this one).
pub(crate) fn omission_charge_addend_except(
    profile: &Profile,
    explained: &std::collections::HashSet<FeatureKey>,
    window_obs: &HashMap<FeatureKey, f64>,
    sample: &SampleSideBackgrounds,
) -> f64 {
    if explained.is_empty() {
        return omission_charge_addend(profile, window_obs, sample);
    }
    let omitted: Vec<&FeatureKey> = window_obs
        .keys()
        .filter(|feature| !profile.contains_key(*feature) && !explained.contains(*feature))
        .collect();
    omission_term_sum(&omitted, window_obs, sample)
}

/// ONE omission term's value (the term arithmetic's single spelling; see
/// `omission_term_sum` for the model semantics). Shared by the per-window and
/// the chain-level addends AND by the pairwise-coupling search objective's
/// per-window term tables, so every exonation channel prices a feature with
/// the SAME value.
pub(crate) fn omission_term_value(
    feature: &FeatureKey,
    observed: f64,
    sample: &SampleSideBackgrounds,
) -> f64 {
    let beta = sample.beta(feature);
    // The pure-background Poisson charge beta_f - C_f*ln(beta_f) plus
    // the count's log-factorial ln(Gamma(C_f+1)): the negative
    // log-likelihood of observing C_f reads under pure background must
    // be bounded below and grow with C_f (the Poisson deviance shape
    // C*(ln(C/beta) - 1) + ...). Without the log-factorial the bare
    // form goes NEGATIVE once C_f exceeds beta_f/ln(beta_f) — it would
    // REWARD omitting a high-background (flood) feature with large
    // observed support (measured: -371M haploid totals on try 3). The
    // log-factorial is C-dependent and fully derived, not a tuning
    // constant; the EXPLAINED features keep the established baselined
    // form bit-exactly (its credit term never sums C unsandwiched, so
    // it never needed the correction).
    beta - observed * beta.ln() + ln_gamma_observation(observed)
}

/// The omission terms over an explicit feature list, ascending key order (the
/// observed maps are hash maps, whose iteration order must never enter an f64
/// sum). Shared by the per-window and the chain-level addends so the term
/// arithmetic has exactly one spelling.
fn omission_term_sum(
    omitted: &[&FeatureKey],
    window_obs: &HashMap<FeatureKey, f64>,
    sample: &SampleSideBackgrounds,
) -> f64 {
    let mut omitted: Vec<&FeatureKey> = omitted.to_vec();
    // Ascending key order: the observed maps are hash maps, whose iteration
    // order must never enter an f64 sum.
    omitted.sort();
    let mut total = 0.0;
    for feature in omitted {
        let observed = window_obs[feature];
        total += omission_term_value(feature, observed, sample);
    }
    total
}

/// A window's observed-feature index: the features sorted ascending, their
/// dense indices, and each feature's omission-term value at this window. The
/// pairwise-coupling search objective's exonation channel works on these
/// indices (the exoneration lists are sorted index vectors; the weighted
/// sums recover the terms by index), so the per-pair work is a merge join
/// over integers instead of a hash sweep.
pub(crate) struct WindowObsIndex {
    pub(crate) features: Vec<FeatureKey>,
    pub(crate) index: HashMap<FeatureKey, u32>,
    pub(crate) terms: Vec<f64>,
}

pub(crate) fn build_window_obs_index(
    window_obs: &HashMap<FeatureKey, f64>,
    sample: &SampleSideBackgrounds,
) -> WindowObsIndex {
    let mut features: Vec<FeatureKey> = window_obs.keys().cloned().collect();
    features.sort();
    let index: HashMap<FeatureKey, u32> = features
        .iter()
        .enumerate()
        .map(|(i, key)| (key.clone(), i as u32))
        .collect();
    let terms: Vec<f64> = features
        .iter()
        .map(|feature| omission_term_value(feature, window_obs[feature], sample))
        .collect();
    WindowObsIndex {
        features,
        index,
        terms,
    }
}

/// The IN-WINDOW observed profiles (the omission charge's universes — owner
/// ruling 2026-09-24, reading (B)): per slice locus, per feature, the
/// observed support PHYSICALLY INSIDE the window's territory. A record
/// contributes its equal share m_r/t_r for each FORWARD subwalk instance of
/// the feature whose anchors' k-spans touch one of the window's owning
/// partitions' territory intervals (the placement-fraction derivation,
/// applied per feature: m_r * |feature placements of r inside the site| /
/// t_r). The established record-touching maps charge a touching record's
/// WHOLE walk features — observed mass that physically lives elsewhere —
/// which made every candidate's omission cost ~50-65k/locus (uniform, no
/// discrimination; measured run try-2). The in-window form is the honest
/// in-window observed DNA. No constants: share, multiplicity, touched count
/// and the territory intervals are all existing measured quantities.
pub(crate) fn build_in_window_obs(
    routed: &[RoutedRecord],
    index: &TerritoryIndex,
    k: u64,
    window_owner_sets: &[BTreeSet<u32>],
    // OUT (the S2 instance-level exoneration): the records' GLOBAL
    // placement spans per feature (indexed by record) and the windows'
    // contributing record lists. Built in the SAME pass with the SAME
    // enumeration and the SAME map accumulation order, so the observed
    // maps are bit-identical with and without it.
    spans_out: Option<&mut RecordPlacementSpans>,
    records_out: Option<&mut WindowRecordLists>,
) -> io::Result<Vec<HashMap<FeatureKey, f64>>> {
    // partition -> the loci whose owner set contains it (a partition owns
    // rows in few windows; dual-role groups appear in their anchored window
    // plus every window whose added rows they own).
    let mut partition_loci: HashMap<u32, Vec<usize>> = HashMap::new();
    for (locus, owners) in window_owner_sets.iter().enumerate() {
        for &owner in owners {
            partition_loci.entry(owner).or_default().push(locus);
        }
    }
    let mut maps: Vec<HashMap<FeatureKey, f64>> =
        (0..window_owner_sets.len()).map(|_| HashMap::new()).collect();
    let collect_instances = spans_out.is_some();
    let mut record_spans: RecordPlacementSpans = Vec::with_capacity(routed.len());
    let mut window_records: WindowRecordLists =
        (0..window_owner_sets.len()).map(|_| HashMap::new()).collect();
    for (record_index, record) in routed.iter().enumerate() {
        // The record's placement-span slot is pushed FIRST so the index
        // alignment with the routed slice survives the skip paths.
        record_spans.push(HashMap::new());
        let anchors = match decode_tokens(&record.tokens) {
            Ok(anchors) => anchors,
            Err(_) => continue,
        };
        let n = anchors.len();
        if n == 0 || record.forward_positions.is_empty() {
            continue;
        }
        let share = record.multiplicity as f64 / record.occurrences.len() as f64;
        // Per (record, feature): the touched-partition UNION across the
        // record's forward occurrences. The site attribution is the
        // ESTABLISHED share semantics restricted to in-window features: a
        // record contributes m_r/t_r ONCE per site for a feature it
        // realizes in-window (the established maps count each record once
        // per touched partition, regardless of its occurrence count there;
        // counting per occurrence would multiply the flood records'
        // features by their in-window placement count — measured: 29M/locus
        // omission totals on try 4).
        let mut per_feature: HashMap<FeatureKey, Vec<u32>> = HashMap::new();
        // The per-instance placement structure (collected only when the
        // instance-level exoneration consumes it): per (occurrence,
        // subwalk) the subwalk anchors' k-mer extent on its path plus the
        // instance's own touched partitions.
        // The record's GLOBAL placement spans per feature (the S2
        // coverage side: the record's own forward-placement set for the
        // feature, across all paths and occurrences).
        let mut per_feature_spans: HashMap<FeatureKey, Vec<(usize, u64, u64)>> =
            HashMap::new();
        // BOTH frames' placements (the strand-asymmetry fix): the record's
        // forward occurrences with its own anchors, and the rc orientation's
        // occurrences with the rc walk's anchors (the rc-frame mem
        // extraction's own qualifying k-mers). An rc occurrence is a real
        // occurrence of the same features: the in-window attribution and the
        // S2 record spans must see it.
        let rc_anchors = reverse_complement_walk(&anchors, k);
        let orientations: [(&[(i32, u64)], &[(usize, u64)]); 2] = [
            (&anchors[..], &record.forward_positions[..]),
            (&rc_anchors[..], &record.reverse_positions[..]),
        ];
        for (walk, positions) in orientations {
            let n_walk = walk.len();
            for &(path_idx, occurrence_start) in positions {
                // Per-anchor partition lists, computed once per occurrence.
                let anchor_parts: Vec<Vec<u32>> = walk
                    .iter()
                    .map(|&(_node, rel)| index.partitions_at(path_idx, occurrence_start + rel, k))
                    .collect();
                for i in 0..n_walk {
                    for j in i..n_walk {
                        let encoded = match encode_walk(&walk[i..=j]) {
                            Ok(encoded) => encoded,
                            Err(_) => continue,
                        };
                        let feature = canonical(&encoded);
                        let entry = per_feature.entry(feature.clone()).or_default();
                        for list in &anchor_parts[i..=j] {
                            for &partition in list {
                                if !entry.contains(&partition) {
                                    entry.push(partition);
                                }
                            }
                        }
                        if collect_instances {
                            per_feature_spans.entry(feature).or_default().push((
                                path_idx,
                                occurrence_start + walk[i].1,
                                occurrence_start + walk[j].1 + k,
                            ));
                        }
                    }
                }
            }
        }
        if collect_instances {
            record_spans[record_index] = per_feature_spans;
        }
        for (feature, touched) in per_feature {
            // The loci whose owner sets the record's in-window feature
            // instances touch; the record contributes its share ONCE per
            // site (the established per-(record, partition) dedup).
            let mut loci_hit: Vec<usize> = Vec::new();
            for &partition in &touched {
                if let Some(loci) = partition_loci.get(&partition) {
                    for &locus in loci {
                        if !loci_hit.contains(&locus) {
                            loci_hit.push(locus);
                        }
                    }
                }
            }
            for locus in loci_hit {
                *maps[locus].entry(feature.clone()).or_default() += share;
                if collect_instances {
                    window_records[locus]
                        .entry(feature.clone())
                        .or_default()
                        .push(record_index as u32);
                }
            }
        }
    }
    if let Some(out) = spans_out {
        *out = record_spans;
    }
    if let Some(out) = records_out {
        *out = window_records;
    }
    Ok(maps)
}

// ---------------------------------------------------------------------------
// INSTANCE-LEVEL EXONERATION — S2, RECORD-LEVEL PLACEMENT COVERAGE (the
// corrected Ruling 2, supervisor-approved on the measured placement
// structure 2026-09-24): a record's share of a feature at window W is
// covered iff the chain spells an anchor OF THAT FEATURE overlapping one of
// THE RECORD'S OWN placement spans (its global forward-placement set for
// that feature) — one read = one observation, already share-split by the
// routing. The cases the ruling distinguishes, honored by the same rule:
// (a) the same DNA piece serving several overlapping windows: the windows'
// records place AT the piece's coordinates, one spelling covers them all —
// charged once, as today; (b) a repeat feature at DISTINCT loci with
// DIFFERENT READS: a record that never places at the spelled copy stays
// charged — the patchwork's free-riding dies at READ granularity; the
// flood records (one read placing across many panel paths/loci — measured:
// the pooled's top truth-advantage features): covered by the one spelling
// at the read's true home (the literal-positional form S1 charged the same
// read at every window — measured +1M inflation on every chain and the
// truth's margin vs the patchwork INVERTING; refuted by the pooled
// external's own read-level semantics). Spelled side: the chain's pieces'
// read windows' realized records (the same event-run machinery the
// profiles charge with, carrying each subwalk's anchor positions).
// Coverage test = interval overlap of the anchor extents on the same panel
// path — the same overlap test the in-window attribution applies to
// territory.
// ---------------------------------------------------------------------------

/// One spelled feature position: the chain spells the feature's anchors
/// across the path interval [lo, hi) on the panel path `path`.
#[derive(Clone, Debug)]
pub(crate) struct SpelledSpan {
    pub(crate) path: usize,
    pub(crate) lo: u64,
    pub(crate) hi: u64,
}

/// The records' GLOBAL placement spans per feature (the S2 coverage side):
/// per routed record, per feature the record realizes, its forward-
/// placement subwalk anchors' k-mer extents across all paths. Indexed by
/// the record's position in the routed slice (the window record lists
/// reference these indices).
pub(crate) type RecordPlacementSpans = Vec<HashMap<FeatureKey, Vec<(usize, u64, u64)>>>;

/// Per window, per feature: the contributing RECORDS' indices (one entry
/// per (record, locus, feature) — the established record-once mass
/// semantics; the record's equal share rides its index).
pub(crate) type WindowRecordLists = Vec<HashMap<FeatureKey, Vec<u32>>>;

/// The S2 instance structure as one bundle (the threaded parameter): the
/// records' global placement spans, their equal shares, and the windows'
/// contributing record lists.
#[derive(Default)]
pub(crate) struct InstanceStructure {
    pub(crate) record_spans: RecordPlacementSpans,
    pub(crate) record_shares: Vec<f64>,
    pub(crate) window_records: WindowRecordLists,
}

/// The coverage test (ONE spelling for the S2 rule): the record's own
/// placement spans overlap a spelled span on the same panel path.
pub(crate) fn record_covered(
    record_spans: &[(usize, u64, u64)],
    spelled: &[SpelledSpan],
) -> bool {
    record_spans
        .iter()
        .any(|&(path, lo, hi)| {
            spelled
                .iter()
                .any(|s| s.path == path && s.lo < hi && lo < s.hi)
        })
}

/// The spelled-piece span memo's type (keyed WITH the piece's orientation:
/// the span-to-path mapping is orientation dependent, unlike the canonical
/// profiles' memo key).
pub(crate) type SpelledMemo = HashMap<
    (usize, u64, u64, bool),
    std::rc::Rc<BTreeMap<FeatureKey, Vec<(u64, u64)>>>,
>;

/// The sources -> panel-path map (the spelled spans' path coordinates and
/// the observed instances' path indices must name the same path space for
/// the coverage test; derived from the lanes' names, the same derivation
/// the run's own `path_of_source` uses).
fn sources_path_map(panel: &SyngIndex, sources: &routes::Sources) -> io::Result<Vec<usize>> {
    let path_of_name: HashMap<String, usize> = panel
        .name_map
        .name_to_path
        .iter()
        .map(|(name, &path)| (name.clone(), path as usize))
        .collect();
    sources
        .lanes
        .iter()
        .map(|(name, _)| {
            path_of_name
                .get(name)
                .copied()
                .ok_or_else(|| invalid("route lane absent from panel"))
        })
        .collect()
}

/// Per-feature PIECE-RELATIVE anchor spans of one spelled segment sequence:
/// where the sequence's read windows realize the feature's subwalks. The
/// event-run invariance (every read start in a run realizes the same
/// records with the same within-read anchor positions) makes the run's
/// spelled extent for one (record, subwalk) the run's start range crossed
/// with the anchors' read extent. Same enumeration the profiles charge
/// with (`profile_event_runs`); the anchor positions come from the tagged
/// records (the input-forward positions `canonical_mem_records`
/// discards). No model arithmetic here — positions only.
pub(crate) fn oracle_segment_spans(
    panel: &SyngIndex,
    sequence: &[u8],
) -> io::Result<BTreeMap<FeatureKey, Vec<(u64, u64)>>> {
    let read_length = READ_LENGTH;
    let mut spans: BTreeMap<FeatureKey, Vec<(u64, u64)>> = BTreeMap::new();
    if sequence.len() < read_length {
        return Ok(spans);
    }
    let events = genome::event_boundaries(
        panel,
        sequence,
        read_length,
        0,
        sequence.len() - read_length + 1,
    )?;
    let k = panel.syncmer_length_bp() as u64;
    for run in events.windows(2) {
        let start = run[0] as usize;
        let run_end = run[1] as usize; // exclusive; read starts start..run_end-1
        for record in mem_records::tagged_mem_records(panel, &sequence[start..start + read_length])? {
            // Every entry is an anchor (signed node, input-forward k-mer
            // start); subwalks are the contiguous slices.
            for i in 0..record.len() {
                for j in i..record.len() {
                    let encoded = encode_walk(&record[i..=j])?;
                    let feature = canonical(&encoded);
                    // The anchor's PATH position is CONSTANT across the
                    // run: each read start of the run re-observes the SAME
                    // genomic occurrence (the within-read position shifts
                    // with the window; piece_start + run_start + read_pos
                    // is invariant). The spelled extent = the subwalk's
                    // anchor k-mer extent at that occurrence.
                    let lo = start as u64 + record[i].1;
                    let hi = start as u64 + record[j].1 + k;
                    spans.entry(feature).or_default().push((lo, hi));
                }
            }
        }
    }
    Ok(spans)
}

/// The seam context's per-feature spans, SPLIT so no span crosses the
/// left/right boundary (a junction-straddling subwalk's anchors distribute
/// between the two pieces' sources). Spans are (is_left, lo, hi) in
/// CONTEXT-relative coordinates; the caller maps each side through its
/// piece. Same window restriction `profile_event_seam` charges with (only
/// read starts crossing the boundary).
fn profile_event_seam_spans(
    panel: &SyngIndex,
    left: &[u8],
    right: &[u8],
    read_length: usize,
) -> io::Result<BTreeMap<FeatureKey, Vec<(bool, u64, u64)>>> {
    let mut spans: BTreeMap<FeatureKey, Vec<(bool, u64, u64)>> = BTreeMap::new();
    if read_length <= 1 || left.is_empty() || right.is_empty() {
        return Ok(spans);
    }
    let left_take = left.len().min(read_length - 1);
    let right_take = right.len().min(read_length - 1);
    let mut context = left[left.len() - left_take..].to_vec();
    let boundary = context.len();
    context.extend_from_slice(&right[..right_take]);
    if context.len() < read_length {
        return Ok(spans);
    }
    let lo = boundary.saturating_add(1).saturating_sub(read_length);
    let hi = boundary.min(context.len() - read_length + 1);
    if lo >= hi {
        return Ok(spans);
    }
    let events = genome::event_boundaries(panel, &context, read_length, lo, hi)?;
    let k = panel.syncmer_length_bp() as u64;
    for run in events.windows(2) {
        let start = run[0] as usize;
        let run_end = run[1] as usize;
        for record in mem_records::tagged_mem_records(panel, &context[start..start + read_length])? {
            for i in 0..record.len() {
                for j in i..record.len() {
                    let encoded = encode_walk(&record[i..=j])?;
                    let feature = canonical(&encoded);
                    // Split the subwalk's anchors at the boundary: maximal
                    // same-side anchor runs, one span per run.
                    let mut a = i;
                    while a <= j {
                        let is_left = record[a].1 < boundary as u64;
                        let mut b = a;
                        while b < j
                            && (record[b + 1].1 < boundary as u64) == is_left
                        {
                            b += 1;
                        }
                        let span_lo = start as u64 + record[a].1;
                        let span_hi = start as u64 + record[b].1 + k;
                        spans
                            .entry(feature.clone())
                            .or_default()
                            .push((is_left, span_lo, span_hi));
                        a = b + 1;
                    }
                }
            }
        }
    }
    Ok(spans)
}

/// A piece-relative span mapped to PATH coordinates (piece tuple: source,
/// start, end, reverse). The reverse piece's rc span [o_lo, o_hi) maps to
/// the source interval [end - o_hi, end - o_lo): rc(S)[o..o+k] =
/// rc(S[L-o-k..L-o]), so the rc k-mer at offset o sits at source
/// [end-o-k, end-o).
fn piece_span_coords(piece: (usize, u64, u64, bool), lo: u64, hi: u64) -> (u64, u64) {
    if piece.3 {
        (piece.2 - hi, piece.2 - lo)
    } else {
        (piece.1 + lo, piece.1 + hi)
    }
}

/// The spelled-span derivation diagnosis (env-gated,
/// IMPG_SPELLED_SPAN_DIAG=spec.json): for the configured features, when a
/// spelled piece produces spans for one of them, dump the FULL derivation
/// chain — the piece tuple, the piece-relative spans, the piece sequence's
/// own anchored syncmers around the span (what the tagged records saw), and
/// the panel path's steps at the mapped coordinates (what the placements
/// and the coverage test see) — then exit. The phantom-span question (a
/// spelled span on a path whose steps lack the feature's nodes) is decided
/// by comparing the two node lists at the same coordinates.
fn spelled_span_derivation_diagnosis(
    piece_spans: &BTreeMap<FeatureKey, Vec<(u64, u64)>>,
    piece: (usize, u64, u64, bool),
    sequence: &[u8],
    panel: &SyngIndex,
    path_of_source: &[usize],
) -> io::Result<()> {
    let Ok(spec_path) = std::env::var("IMPG_SPELLED_SPAN_DIAG") else {
        return Ok(());
    };
    #[derive(serde::Deserialize)]
    struct Spec {
        features: Vec<Vec<u64>>,
        /// When set, only pieces mapped to these panel paths are dumped (the
        /// phantom-span pieces); empty = any target-feature piece.
        #[serde(default)]
        watch_paths: Vec<u64>,
        #[serde(default = "spelled_diag_default_cap")]
        max_rows: usize,
        out: String,
    }
    fn spelled_diag_default_cap() -> usize {
        24
    }
    let spec: Spec = read_json(std::path::Path::new(&spec_path))?;
    let path = path_of_source[piece.0];
    if !spec.watch_paths.is_empty() && !spec.watch_paths.contains(&(path as u64)) {
        return Ok(());
    }
    let mut rows = Vec::new();
    for (feature, spans) in piece_spans {
        if !spec.features.contains(feature) {
            continue;
        }
        let Some(&(lo, hi)) = spans.first() else {
            continue;
        };
        let (plo, phi) = piece_span_coords(piece, lo, hi);
        // The piece sequence's own anchored syncmers around the span (the
        // spelled side's view: what the tagged records anchored against).
        let seq_lo = lo.saturating_sub(2 * k_window(panel) as u64) as usize;
        let seq_hi = (hi + 2 * k_window(panel) as u64).min(sequence.len() as u64) as usize;
        let slice = &sequence[seq_lo..seq_hi];
        let piece_side_syncmers: Vec<[i64; 2]> = panel
            .matched_syncmers_in_sequence(slice)
            .iter()
            .take(32)
            .map(|s| [s.signed_node as i64, (seq_lo as u64 + s.query_pos) as i64])
            .collect();
        // The panel path's steps at the mapped coordinates (the coverage
        // test's view).
        let path_side_steps: Vec<[i64; 2]> = panel
            .walk_path_range(path, plo.saturating_sub(2 * k_window(panel) as u64), phi + 2 * k_window(panel) as u64)
            .map(|steps| {
                steps
                    .iter()
                    .take(32)
                    .map(|&(node, bp)| [node as i64, bp as i64])
                    .collect()
            })
            .unwrap_or_else(|_| Vec::new());
        rows.push(serde_json::json!({
            "feature": feature,
            "piece": {"source": piece.0, "start": piece.1, "end": piece.2, "reverse": piece.3},
            "piece_relative_span": [lo, hi],
            "path": path,
            "path_span": [plo, phi],
            "piece_relative_spans_all": spans.iter().take(8).cloned().collect::<Vec<_>>(),
            "piece_sequence_window": String::from_utf8_lossy(slice).to_string(),
            "piece_sequence_window_offset": seq_lo,
            "piece_side_syncmers": piece_side_syncmers,
            "path_side_steps": path_side_steps,
        }));
        if rows.len() >= 8 {
            break;
        }
    }
    if rows.is_empty() {
        return Ok(());
    }
    let mut existing: Vec<serde_json::Value> = std::fs::read(&spec.out)
        .ok()
        .and_then(|bytes| serde_json::from_slice(&bytes).ok())
        .unwrap_or_default();
    let seen_pieces: std::collections::HashSet<[i64; 4]> = existing
        .iter()
        .filter_map(|row| row.get("piece").and_then(|p| p.as_object()))
        .filter_map(|p| {
            Some([
                p.get("source")?.as_i64()?,
                p.get("start")?.as_i64()?,
                p.get("end")?.as_i64()?,
                p.get("reverse")?.as_i64()?,
            ])
        })
        .collect();
    if seen_pieces.contains(&[piece.0 as i64, piece.1 as i64, piece.2 as i64, piece.3 as i64]) {
        return Ok(());
    }
    existing.extend(rows);
    std::fs::write(&spec.out, serde_json::to_vec_pretty(&existing)?)?;
    eprintln!("[spelled-diag] appended ({} distinct pieces)", existing.len());
    if existing.len() >= spec.max_rows {
        eprintln!("[spelled-diag] complete, exiting");
        std::process::exit(0);
    }
    Ok(())
}

fn k_window(panel: &SyngIndex) -> u64 {
    panel.syncmer_length_bp() as u64
}

/// A spelled piece's piece-relative spans mapped to PATH coordinates and
/// merged into the chain's spelled map.
fn piece_spans_to_path(
    path_of_source: &[usize],
    piece: (usize, u64, u64, bool),
    spans: &BTreeMap<FeatureKey, Vec<(u64, u64)>>,
    out: &mut HashMap<FeatureKey, Vec<SpelledSpan>>,
) -> io::Result<()> {
    let path = *path_of_source
        .get(piece.0)
        .ok_or_else(|| invalid("piece source absent from the path map"))?;
    for (feature, list) in spans {
        let entry = out.entry(feature.clone()).or_default();
        for &(lo, hi) in list {
            let (plo, phi) = piece_span_coords(piece, lo, hi);
            entry.push(SpelledSpan { path, lo: plo, hi: phi });
        }
    }
    Ok(())
}

/// The instance-level exoneration's PER-FEATURE probe report (the ladder's
/// arithmetic-first validation): for the configured features, per window,
/// the observed mass split into covered/uncovered under the evaluated
/// chain's spelled anchors, plus the feature-level (old Ruling 2) charge
/// for comparison — the old form exonerates a feature at EVERY window once
/// the chain spells it anywhere; the corrected form charges the uncovered
/// instances' mass. The old charge derives from the same spelled map (a
/// feature with any spelled span is old-exonerated) — the sets coincide by
/// construction (both enumerate the same records/subwalks).
pub(crate) fn instance_probe_report(
    window_obs: &[HashMap<FeatureKey, f64>],
    instances: &InstanceStructure,
    spelled: &HashMap<FeatureKey, Vec<SpelledSpan>>,
    sample: &SampleSideBackgrounds,
    features: &[FeatureKey],
) -> serde_json::Value {
    let rows: Vec<serde_json::Value> = features
        .iter()
        .map(|feature| {
            let spans = spelled.get(feature);
            let old_charge = if spans.is_some() {
                0.0
            } else {
                window_obs
                    .iter()
                    .filter_map(|w| w.get(feature))
                    .map(|&full| omission_term_value(feature, full, sample))
                    .sum::<f64>()
            };
            let mut new_charge = 0.0f64;
            let windows: Vec<serde_json::Value> = window_obs
                .iter()
                .enumerate()
                .filter_map(|(locus, w)| {
                    let full = *w.get(feature)?;
                    let covered = match (instances.window_records.get(locus), spans) {
                        (Some(list), Some(spans)) if !spans.is_empty() => list
                            .get(feature)
                            .map(|records| {
                                records
                                    .iter()
                                    .filter(|&record_index| {
                                        record_covered(
                                            instances
                                                .record_spans
                                                .get(*record_index as usize)
                                                .and_then(|m| m.get(feature))
                                                .map(|v| v.as_slice())
                                                .unwrap_or(&[]),
                                            spans,
                                        )
                                    })
                                    .map(|&record_index| {
                                        instances
                                            .record_shares
                                            .get(record_index as usize)
                                            .copied()
                                            .unwrap_or(0.0)
                                    })
                                    .sum::<f64>()
                            })
                            .unwrap_or(0.0),
                        _ => 0.0,
                    };
                    let uncovered = (full - covered).max(0.0);
                    if uncovered > 0.0 {
                        new_charge += omission_term_value(feature, uncovered, sample);
                    }
                    Some(serde_json::json!({
                        "locus": locus,
                        "full": full,
                        "covered": covered,
                        "uncovered": uncovered,
                    }))
                })
                .collect();
            // The record-atomic coverage choice's exposure: records with
            // MULTIPLE placement spans for the feature (any-covered ->
            // covered). Measured, reported — not assumed.
            let mut multi_span_entries = 0usize;
            let mut total_entries = 0usize;
            for spans_map in &instances.record_spans {
                if let Some(spans) = spans_map.get(feature) {
                    total_entries += 1;
                    multi_span_entries += usize::from(spans.len() > 1);
                }
            }
            let spelled_dump: Vec<serde_json::Value> = spans
                .map(|list| {
                    list.iter()
                        .take(6)
                        .map(|s| serde_json::json!({"path": s.path, "lo": s.lo, "hi": s.hi}))
                        .collect()
                })
                .unwrap_or_default();
            let instance_dump: Vec<serde_json::Value> = instances
                .record_spans
                .iter()
                .filter_map(|map| map.get(feature))
                .take(6)
                .map(|spans| {
                    serde_json::json!({
                        "spans": spans.iter().take(6).map(|(p, lo, hi)|
                            serde_json::json!({"path": p, "lo": lo, "hi": hi}))
                            .collect::<Vec<_>>(),
                    })
                })
                .collect();
            serde_json::json!({
                "feature": feature,
                "spelled_span_count": spans.map(|s| s.len()).unwrap_or(0),
                "spelled_dump": spelled_dump,
                "instance_dump": instance_dump,
                "old_charge": old_charge,
                "new_charge": new_charge,
                "new_minus_old": new_charge - old_charge,
                "instance_entries": total_entries,
                "multi_span_entries": multi_span_entries,
                "windows": windows,
            })
        })
        .collect();
    serde_json::json!({ "features": rows })
}

/// The HAPLOID single-allele local charge under the omission charge: the
/// row's current charge (its own features against its own owner-resolved
/// observed side — bit-identical to the established form) plus the omission
/// addend over the window's full observed profile.
pub(crate) fn merged_single_loss_sample_with_omission(
    profile: &Profile,
    obs_row: &HashMap<FeatureKey, f64>,
    window_obs: &HashMap<FeatureKey, f64>,
    sample: &SampleSideBackgrounds,
    model: &ScoreModel,
) -> io::Result<f64> {
    let empty = std::collections::HashSet::new();
    merged_single_loss_sample_with_omission_except(
        profile,
        obs_row,
        &empty,
        window_obs,
        sample,
        model,
    )
}

/// The chain-level form (see `omission_charge_addend_except`): the row keeps
/// its current charge bit-exactly and pays the omission cost of the window's
/// observed features it neither predicts nor is exonerated from by the
/// chain's explained set.
pub(crate) fn merged_single_loss_sample_with_omission_except(
    profile: &Profile,
    obs_row: &HashMap<FeatureKey, f64>,
    explained: &std::collections::HashSet<FeatureKey>,
    window_obs: &HashMap<FeatureKey, f64>,
    sample: &SampleSideBackgrounds,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = merged_single_loss_sample(profile, obs_row, sample, model)?;
    loss += omission_charge_addend_except(profile, explained, window_obs, sample);
    Ok(loss)
}

/// The INSTANCE-LEVEL chain-level omission addend (the corrected Ruling 2):
/// the row keeps its current charge bit-exactly; the window's omission
/// covers the in-window observed features the row does not predict, at the
/// observed mass the chain's spelled anchors do NOT cover — a feature's
/// in-window observed mass is covered iff the chain spells an anchor of
/// that feature overlapping the observed instance's placement span on the
/// same path. The term is charged at the UNCOVERED mass (the Poisson
/// count's own argument); the covered/uncovered split conserves the
/// window's observed mass exactly (uncovered = full - covered, so the
/// nothing-covered case charges the established full-mass term
/// bit-exactly). No instance structure (or no spelled positions) falls
/// back to the established feature-level form: full mass, exonerated iff
/// the explained set contains the feature.
pub(crate) fn merged_single_loss_sample_with_omission_instances(
    profile: &Profile,
    obs_row: &HashMap<FeatureKey, f64>,
    window_obs: &HashMap<FeatureKey, f64>,
    // The S2 coverage inputs: this window's contributing records (their
    // indices) and the records' global placement spans + equal shares.
    window_records: Option<&HashMap<FeatureKey, Vec<u32>>>,
    instances: Option<&InstanceStructure>,
    spelled: &HashMap<FeatureKey, Vec<SpelledSpan>>,
    sample: &SampleSideBackgrounds,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = merged_single_loss_sample(profile, obs_row, sample, model)?;
    let mut omitted: Vec<&FeatureKey> = window_obs
        .keys()
        .filter(|feature| !profile.contains_key(*feature))
        .collect();
    // Ascending key order: the observed maps are hash maps, whose
    // iteration order must never enter an f64 sum.
    omitted.sort();
    for feature in omitted {
        let full = window_obs[feature];
        let covered = match (window_records, instances, spelled.get(feature)) {
            (Some(list), Some(instances), Some(spans)) if !spans.is_empty() => {
                match list.get(feature) {
                    Some(records) => records
                        .iter()
                        .filter(|&record_index| {
                            record_covered(
                                instances.record_spans
                                    .get(*record_index as usize)
                                    .and_then(|m| m.get(feature))
                                    .map(|v| v.as_slice())
                                    .unwrap_or(&[]),
                                spans,
                            )
                        })
                        .map(|&record_index| {
                            instances
                                .record_shares
                                .get(record_index as usize)
                                .copied()
                                .unwrap_or(0.0)
                        })
                        .sum::<f64>(),
                    // Observed mass without a record list (a decode-skipped
                    // record cannot contribute mass, so this is not reached
                    // by the builder; the safe fallback charges full).
                    None => 0.0,
                }
            }
            _ => 0.0,
        };
        let uncovered = (full - covered).max(0.0);
        if uncovered > 0.0 {
            loss += omission_term_value(feature, uncovered, sample);
        }
    }
    Ok(loss)
}

/// All subwalk features of a record's token walk (the per-feature feature
/// list the share accumulation and the site maps charge).
fn enumerate_subwalks(tokens: &[u64]) -> io::Result<Vec<FeatureKey>> {
    let anchors = decode_tokens(tokens)?;
    let n = anchors.len();
    let mut out = Vec::with_capacity(n * (n + 1) / 2);
    for i in 0..n {
        for j in i..n {
            let encoded = encode_walk(&anchors[i..=j])?;
            out.push(canonical(&encoded));
        }
    }
    Ok(out)
}
// ---------------------------------------------------------------------------
// Geometric predicted profiles.
// ---------------------------------------------------------------------------

fn contained_path_anchors(
    panel: &SyngIndex,
    path_idx: usize,
    start: u64,
    end: u64,
    k: u64,
) -> io::Result<Vec<(i32, u64)>> {
    let steps = panel.walk_path_range(path_idx, start, end + k.saturating_sub(1))?;
    Ok(steps
        .into_iter()
        .filter(|&(_node, bp)| bp >= start && bp + k <= end)
        .map(|(node, bp)| (node, bp - start))
        .collect())
}

fn merge_profiles(profiles: &[&Profile]) -> io::Result<Profile> {
    let mut merged = Profile::new();
    for profile in profiles {
        for (key, &count) in profile.iter() {
            let value = merged.entry(key.clone()).or_default();
            *value = value
                .checked_add(count)
                .ok_or_else(|| invalid("profile merge overflow"))?;
        }
    }
    Ok(merged)
}

// ---------------------------------------------------------------------------
// Geometric inverted-repeat part (the RC view): one padded reverse-complement
// MEM extraction per candidate territory, then a joint window-event sweep with
// the self view and per-window maximality re-filtering.
// ---------------------------------------------------------------------------

/// Number of window starts s that contain the anchor-pair span, given an
/// anchor at candidate position `pos`: contained iff s <= pos and
/// pos + k <= s + READ_LENGTH, i.e. s in [pos + k - READ_LENGTH, pos].
const fn anchor_contain_tail(k: u64) -> u64 {
    READ_LENGTH as u64 - k
}

/// RC-view maximal MEM records of one padded candidate extraction, in
/// candidate input coordinates (positions relative to the candidate start).
/// The extraction is OWN-ORIENTATION on the padded reverse complement (the
/// raw syncmer phase of the RC query, matched against the panel and walked
/// through the GBWT): this surfaces the inverted-duplicate loci's records,
/// which the public best-orientation matcher hides on long queries. Records
/// are then reversed, node-negated, and position-mapped to input-forward
/// coordinates exactly like `collect_tagged_read`'s reverse branch; anchors
/// outside every candidate window's containable range are trimmed from the
/// record ends.
fn geometric_ir_records(
    panel: &SyngIndex,
    sources: &routes::Sources,
    source: usize,
    start: u64,
    end: u64,
    k: u64,
) -> io::Result<Vec<Vec<(i32, u64)>>> {
    let core_len = end - start;
    if core_len < READ_LENGTH as u64 {
        return Ok(Vec::new());
    }
    let lane_len = sources.lanes[source].1;
    let pad_lo = start.saturating_sub(SLACK);
    let pad_hi = (end + SLACK).min(lane_len);
    let padded = sources.fetch(source, pad_lo, pad_hi)?;
    let padded_len = padded.len() as u64;
    ensure(padded_len >= core_len, "padded fetch shorter than candidate")?;
    let core_lo = start - pad_lo;
    // Anchors a candidate window can contain: pos in
    // [core_lo, core_lo + core_len - k] (the last window starts at
    // core_len - READ_LENGTH and ends at core_len).
    let anchor_hi = core_lo + core_len - k;
    let reverse = impg::graph::reverse_complement(&padded);
    let raw = mem_records::own_orientation_mem_records(panel, &reverse)?;
    let mut records = Vec::new();
    for record in raw {
        let record: Vec<(i32, u64)> = record
            .iter()
            .rev()
            .map(|&(node, pos)| (-node, padded_len - k - pos))
            .collect();
        // Positions increase, so out-of-range anchors sit only at the
        // record ends; trim to the maximal contiguous in-range slice (a
        // candidate window can never contain an anchor outside it, and the
        // sweep clips per window anyway).
        let first = record.iter().position(|&(_, pos)| pos >= core_lo);
        let last = record.iter().rposition(|&(_, pos)| pos <= anchor_hi);
        if let (Some(first), Some(last)) = (first, last) {
            if last >= first {
                records.push(
                    record[first..=last]
                        .iter()
                        .map(|&(node, pos)| (node, pos - core_lo))
                        .collect(),
                );
            }
        }
    }
    records.sort();
    records.dedup();
    Ok(records)
}

/// Joint window-event sweep over the self view (contained path anchors) and
/// the RC-view records: per window the record set is the union of the
/// per-window clips, filtered by the `maximal_content_records` rule (a
/// record whose (node, position) sequence is contained in a longer present
/// record is dropped), and every surviving record contributes all its
/// node-to-node subwalks with the run's window multiplicity. Window contain
/// events (anchor enters at pos + k - READ_LENGTH, leaves at pos + 1) match
/// the oracle's event-compression boundaries exactly.
fn geometric_sweep_profile(
    self_anchors: &[(i32, u64)],
    ir_records: &[Vec<(i32, u64)>],
    len: u64,
    k: u64,
    label: &str,
) -> io::Result<Profile> {
    let mut profile = Profile::new();
    if len < READ_LENGTH as u64 {
        return Ok(profile);
    }
    let max_start = len - READ_LENGTH as u64;
    let tail = anchor_contain_tail(k);
    let mut events = vec![0u64, max_start + 1];
    let mut push_anchor_events = |pos: u64, events: &mut Vec<u64>| {
        let enter = pos.saturating_sub(tail);
        if enter > 0 && enter <= max_start {
            events.push(enter);
        }
        let leave = pos.saturating_add(1);
        if leave <= max_start {
            events.push(leave);
        }
    };
    for &(_, pos) in self_anchors {
        push_anchor_events(pos, &mut events);
    }
    for record in ir_records {
        for &(_, pos) in record {
            push_anchor_events(pos, &mut events);
        }
    }
    events.sort_unstable();
    events.dedup();
    for pair in events.windows(2) {
        let start = pair[0];
        let multiplicity = pair[1] - start;
        let lo = start;
        let hi = start + tail;
        let mut clips: Vec<Vec<(i32, u64)>> = Vec::new();
        let self_lo = self_anchors.partition_point(|&(_, pos)| pos < lo);
        let self_hi = self_anchors.partition_point(|&(_, pos)| pos <= hi);
        if self_hi > self_lo {
            clips.push(self_anchors[self_lo..self_hi].to_vec());
        }
        for record in ir_records {
            let rlo = record.partition_point(|&(_, pos)| pos < lo);
            let rhi = record.partition_point(|&(_, pos)| pos <= hi);
            if rhi > rlo {
                clips.push(record[rlo..rhi].to_vec());
            }
        }
        if clips.is_empty() {
            continue;
        }
        clips.sort();
        clips.dedup();
        let kept: Vec<&[(i32, u64)]> = clips
            .iter()
            .filter(|clip| {
                !clips
                    .iter()
                    .any(|other| {
                        other.len() > clip.len()
                            && other.windows(clip.len()).any(|window| window == clip.as_slice())
                    })
            })
            .map(|clip| clip.as_slice())
            .collect();
        for clip in kept {
            for i in 0..clip.len() {
                for j in i..clip.len() {
                    let encoded = encode_walk(&clip[i..=j])?;
                    let key = canonical(&encoded);
                    let value = profile.entry(key).or_default();
                    *value = value.saturating_add(multiplicity);
                }
            }
        }
    }
    ensure(
        profile.len() <= MAX_FEATURES,
        &format!(
            "geometric profile over budget: {label}: features={} len={len}",
            profile.len()
        ),
    )?;
    Ok(profile)
}

/// Full geometric candidate profile: contained path anchors (self view) plus
/// the RC-view records (inverted-repeat part), jointly swept. Returns the
/// profile, the model's anchor-node universe (for diff classification), and
/// the number of RC-view records the extraction produced.
fn geometric_candidate_profile(
    panel: &SyngIndex,
    sources: &routes::Sources,
    source: usize,
    start: u64,
    end: u64,
    path_idx: usize,
    k: u64,
    label: &str,
) -> io::Result<(Profile, BTreeSet<i32>, usize)> {
    let anchors = contained_path_anchors(panel, path_idx, start, end, k)?;
    let ir = geometric_ir_records(panel, sources, source, start, end, k)?;
    let ir_count = ir.len();
    let mut nodes = BTreeSet::new();
    for &(node, _) in &anchors {
        nodes.insert(node);
    }
    for record in &ir {
        for &(node, _) in record {
            nodes.insert(node);
        }
    }
    let profile = geometric_sweep_profile(&anchors, &ir, end - start, k, label)?;
    Ok((profile, nodes, ir_count))
}

// ---------------------------------------------------------------------------
// Exit-cut test helpers.
// ---------------------------------------------------------------------------

fn fetch_segment(
    sources: &routes::Sources,
    source: usize,
    start: u64,
    end: u64,
) -> io::Result<Vec<u8>> {
    sources.fetch(source, start, end)
}

fn oracle_segment_profile(panel: &SyngIndex, sequence: &[u8]) -> io::Result<Profile> {
    let (profile, _) = genome::profile_event_interior(panel, sequence, READ_LENGTH, MAX_FEATURES)?;
    Ok(profile)
}

fn oracle_split_profile(
    panel: &SyngIndex,
    sources: &routes::Sources,
    segments: &[(usize, u64, u64)],
) -> io::Result<Profile> {
    ensure(segments.len() == 2, "exit split must be two segments")?;
    let left = fetch_segment(sources, segments[0].0, segments[0].1, segments[0].2)?;
    let right = fetch_segment(sources, segments[1].0, segments[1].1, segments[1].2)?;
    let interior_left = oracle_segment_profile(panel, &left)?;
    let interior_right = oracle_segment_profile(panel, &right)?;
    let (seam, _) =
        genome::profile_event_seam(panel, &left, &right, READ_LENGTH, MAX_FEATURES)?;
    merge_profiles(&[&interior_left, &interior_right, &seam])
}

fn composed_pair_score(
    candidate: &Profile,
    native: &Profile,
    observed: &dyn Fn(&FeatureKey) -> io::Result<f64>,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut keys: Vec<&FeatureKey> = candidate.keys().chain(native.keys()).collect();
    keys.sort_unstable();
    keys.dedup();
    let mut loss = 0.0;
    for key in keys {
        let q = candidate.get(key).copied().unwrap_or(0)
            + native.get(key).copied().unwrap_or(0);
        loss += loss_fractional(model, q, observed(key)?)?;
    }
    ensure(loss.is_finite(), "nonfinite composed score")?;
    Ok(loss)
}

#[derive(Serialize)]
struct FeatureRow {
    feature: FeatureKey,
    q_plateau: u64,
    q_truth: u64,
    pooled_observed: u64,
    routed_observed: f64,
    pooled_delta: f64,
    routed_delta: f64,
    delta_shift: f64,
}

fn feature_deltas(
    plateau: &Profile,
    truth: &Profile,
    native: &Profile,
    pooled_observed: &dyn Fn(&FeatureKey) -> io::Result<u64>,
    routed_observed: &dyn Fn(&FeatureKey) -> f64,
    model: &ScoreModel,
) -> io::Result<Vec<FeatureRow>> {
    let mut keys: Vec<&FeatureKey> = plateau
        .keys()
        .chain(truth.keys())
        .chain(native.keys())
        .collect();
    keys.sort_unstable();
    keys.dedup();
    let mut rows = Vec::with_capacity(keys.len());
    for key in keys {
        let q_plateau = plateau.get(key).copied().unwrap_or(0)
            + native.get(key).copied().unwrap_or(0);
        let q_truth = truth.get(key).copied().unwrap_or(0)
            + native.get(key).copied().unwrap_or(0);
        let pooled = pooled_observed(key)?;
        let routed = routed_observed(key);
        let pooled_delta = loss_fractional(model, q_plateau, pooled as f64)?
            - loss_fractional(model, q_truth, pooled as f64)?;
        let routed_delta = loss_fractional(model, q_plateau, routed)?
            - loss_fractional(model, q_truth, routed)?;
        rows.push(FeatureRow {
            feature: key.clone(),
            q_plateau,
            q_truth,
            pooled_observed: pooled,
            routed_observed: routed,
            pooled_delta,
            routed_delta,
            delta_shift: routed_delta - pooled_delta,
        });
    }
    rows.sort_by(|a, b| {
        b.delta_shift
            .abs()
            .total_cmp(&a.delta_shift.abs())
            .then_with(|| a.feature.cmp(&b.feature))
    });
    Ok(rows)
}

// ---------------------------------------------------------------------------
// Equivalence (geometric q vs oracle q).
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct ProfileDiff {
    identity: String,
    segments: Vec<(usize, u64, u64)>,
    reverse: bool,
    oracle_features: usize,
    geometric_features: usize,
    exact_features: usize,
    oracle_only_features: usize,
    geometric_only_features: usize,
    sum_abs_delta: f64,
    max_abs_delta: u64,
    /// Per-feature diff classification: every differing feature is assigned
    /// exactly one cause class from its node membership and presence pattern
    /// (see `classify_feature`).
    class_counts: BTreeMap<String, usize>,
    class_sum_abs_delta: BTreeMap<String, f64>,
    class_examples: BTreeMap<String, Vec<(FeatureKey, u64, u64)>>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    examples: Vec<(FeatureKey, u64, u64)>,
}

/// Classify one differing feature against the geometric model's anchor-node
/// universe. `model_nodes` collects every node of the self view (contained
/// path anchors) and every node of the RC-view records, so:
/// - `*_model_anchors`: every feature node is in the model universe (up to
///   sign — features are RC-canonicalized): the feature IS derived by the
///   geometric model, and any difference is a per-window multiplicity
///   effect — the oracle's read-local syncmer selection loses or shifts
///   anchors near window edges, and the public matcher's best-orientation
///   pick can drop the self record for a window entirely, both of which the
///   full-context geometric sweep counts through;
/// - `*_read_local_anchors`: some node is outside the model universe — the
///   feature exists only because the oracle's read-local syncmer selection
///   matched a panel node the self/RC derivations never select (a sanity
///   flag: after the RC-view fix this class is expected to be empty).
fn classify_feature(
    key: &FeatureKey,
    oracle: u64,
    geometric: u64,
    model_nodes: &BTreeSet<i32>,
) -> String {
    let mut nodes = Vec::with_capacity((key.len() + 1) / 2);
    if let Ok(decoded) = decode_tokens(key) {
        for (node, _) in decoded {
            nodes.push(node);
        }
    }
    // Features are RC-canonicalized, so a subwalk may be stored in its
    // reversed-negated orientation; membership is up to sign.
    let in_model = nodes
        .iter()
        .all(|node| model_nodes.contains(node) || model_nodes.contains(&-*node));
    let side = if oracle > 0 && geometric > 0 {
        "count_delta"
    } else if oracle > 0 {
        "oracle_only"
    } else {
        "geometric_only"
    };
    format!(
        "{}_{}",
        side,
        if in_model {
            "model_anchors"
        } else {
            "read_local_anchors"
        }
    )
}

#[derive(Serialize)]
struct WindowDiag {
    identity: String,
    window_start: u64,
    window_len: usize,
    oracle_records: Vec<Vec<u64>>,
    geometric_records: Vec<Vec<u64>>,
    matched_anchor_nodes: Vec<i32>,
    path_anchor_nodes: Vec<i32>,
    matched_anchor_positions: Vec<u64>,
    path_anchor_positions: Vec<u64>,
    oracle_matches_geometric: bool,
}

fn profile_equivalence(
    identity: &str,
    segments: &[(usize, u64, u64)],
    reverse: bool,
    oracle: &Profile,
    geometric: &Profile,
    model_nodes: &BTreeSet<i32>,
) -> ProfileDiff {
    let mut exact = 0;
    let mut oracle_only = 0;
    let mut geometric_only = 0;
    let mut sum_abs = 0.0;
    let mut max_abs = 0u64;
    let mut examples = Vec::new();
    let mut class_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut class_sum_abs: BTreeMap<String, f64> = BTreeMap::new();
    let mut class_examples: BTreeMap<String, Vec<(FeatureKey, u64, u64)>> = BTreeMap::new();
    let mut keys: BTreeSet<&FeatureKey> = oracle.keys().chain(geometric.keys()).collect();
    for key in keys.into_iter() {
        let (a, b) = (
            oracle.get(key).copied().unwrap_or(0),
            geometric.get(key).copied().unwrap_or(0),
        );
        if a == b {
            exact += 1;
        } else {
            if a == 0 {
                geometric_only += 1;
            }
            if b == 0 {
                oracle_only += 1;
            }
            let delta = a.max(b) - a.min(b);
            sum_abs += delta as f64;
            max_abs = max_abs.max(delta);
            if examples.len() < 24 {
                examples.push((key.clone(), a, b));
            }
            let class = classify_feature(key, a, b, model_nodes);
            *class_counts.entry(class.clone()).or_default() += 1;
            *class_sum_abs.entry(class.clone()).or_default() += delta as f64;
            let bucket = class_examples.entry(class).or_default();
            if bucket.len() < 12 {
                bucket.push((key.clone(), a, b));
            }
        }
    }
    ProfileDiff {
        identity: identity.to_string(),
        segments: segments.to_vec(),
        reverse,
        oracle_features: oracle.len(),
        geometric_features: geometric.len(),
        exact_features: exact,
        oracle_only_features: oracle_only,
        geometric_only_features: geometric_only,
        sum_abs_delta: sum_abs,
        max_abs_delta: max_abs,
        class_counts,
        class_sum_abs_delta: class_sum_abs,
        class_examples,
        examples,
    }
}

// ---------------------------------------------------------------------------
// Routed-evidence DP (the M1 additive model of Milestone 0.5).
// ---------------------------------------------------------------------------

fn spans_feasible_local(segments: &[SourceRange]) -> bool {
    for i in 0..segments.len() {
        for j in i + 1..segments.len() {
            if segments[i].source == segments[j].source
                && segments[i].start < segments[j].end
                && segments[j].start < segments[i].end
            {
                return false;
            }
        }
    }
    true
}

fn allele_sequence(
    sources: &routes::Sources,
    segments: &[SourceRange],
) -> io::Result<Vec<u8>> {
    let mut sequence = Vec::new();
    for segment in segments {
        if segment.start == segment.end {
            continue;
        }
        let mut part = sources.fetch(segment.source, segment.start, segment.end)?;
        if segment.reverse {
            part = impg::graph::reverse_complement(&part);
        }
        sequence.extend(part);
    }
    Ok(sequence)
}

fn allele_endpoints(
    sources: &routes::Sources,
    memo: &FlankMemo,
    segments: &[SourceRange],
    flank: usize,
) -> io::Result<(Vec<u8>, Vec<u8>)> {
    if segments.len() == 1 {
        let deletion = &segments[0];
        if deletion.start == deletion.end {
            let end = deletion
                .start
                .saturating_add(flank as u64)
                .min(sources.lanes[deletion.source].1);
            let head = fetch_oriented(
                sources,
                memo,
                deletion.source,
                deletion.reverse,
                deletion.start,
                end,
            )?;
            return Ok((head, Vec::new()));
        }
    }
    let mut head = Vec::with_capacity(flank);
    for segment in segments {
        if head.len() == flank {
            break;
        }
        let take = (flank - head.len()).min((segment.end - segment.start) as usize) as u64;
        head.extend(fetch_oriented(
            sources,
            memo,
            segment.source,
            segment.reverse,
            if segment.reverse {
                segment.end - take
            } else {
                segment.start
            },
            if segment.reverse {
                segment.end
            } else {
                segment.start + take
            },
        )?);
    }
    let mut tail_parts = Vec::new();
    let mut tail_len = 0usize;
    for segment in segments.iter().rev() {
        if tail_len == flank {
            break;
        }
        let take = (flank - tail_len).min((segment.end - segment.start) as usize) as u64;
        let part = fetch_oriented(
            sources,
            memo,
            segment.source,
            segment.reverse,
            if segment.reverse {
                segment.start
            } else {
                segment.end - take
            },
            if segment.reverse {
                segment.start + take
            } else {
                segment.end
            },
        )?;
        tail_len += part.len();
        tail_parts.push(part);
    }
    tail_parts.reverse();
    Ok((head, tail_parts.concat()))
}

/// Routed-observed loss of one profile against one partition's equal shares.
fn profile_loss_routed(
    profile: &Profile,
    obs: &HashMap<FeatureKey, f64>,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = 0.0;
    for (key, &q) in profile {
        loss += loss_fractional(model, q, obs.get(key).copied().unwrap_or(0.0))?;
    }
    ensure(loss.is_finite(), "nonfinite routed profile loss")?;
    Ok(loss)
}

/// Loss of a boundary seam profile: each feature is charged against the SUM
/// of the two adjacent partitions' routed shares (a junction-window record's
/// anchors straddle the boundary, so its read multiplicity distributes
/// between exactly those two partitions; see Milestone 0.5).
fn profile_loss_boundary(
    profile: &Profile,
    obs_left: &HashMap<FeatureKey, f64>,
    obs_right: &HashMap<FeatureKey, f64>,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = 0.0;
    for (key, &q) in profile {
        let observed = obs_left.get(key).copied().unwrap_or(0.0)
            + obs_right.get(key).copied().unwrap_or(0.0);
        loss += loss_fractional(model, q, observed)?;
    }
    ensure(loss.is_finite(), "nonfinite boundary seam loss")?;
    Ok(loss)
}

/// Per-feature multiplicity variant of `profile_loss_boundary`: each seam
/// feature is charged against its own measured panel multiplicity entry.
fn profile_loss_boundary_multiplicity(
    profile: &Profile,
    obs_left: &HashMap<FeatureKey, f64>,
    obs_right: &HashMap<FeatureKey, f64>,
    model: &ScoreModel,
    backgrounds: &FeatureBackgrounds,
) -> io::Result<f64> {
    let mut loss = 0.0;
    for (key, &q) in profile {
        let observed = obs_left.get(key).copied().unwrap_or(0.0)
            + obs_right.get(key).copied().unwrap_or(0.0);
        loss += loss_fractional_entry(model, q, observed, backgrounds.entry(key))?;
    }
    ensure(loss.is_finite(), "nonfinite boundary seam loss")?;
    Ok(loss)
}

/// Loss of a diploid pair of class profiles against one partition's shares.
fn merged_pair_loss(
    first: &Profile,
    second: &Profile,
    obs: &HashMap<FeatureKey, f64>,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = 0.0;
    let mut left = first.iter().peekable();
    let mut right = second.iter().peekable();
    while left.peek().is_some() || right.peek().is_some() {
        let (q, key) = match (left.peek(), right.peek()) {
            (Some(&(a_key, &a)), Some(&(b_key, &b))) if a_key == b_key => {
                left.next();
                right.next();
                (a + b, a_key)
            }
            (Some(&(a_key, &a)), Some(&(b_key, _))) if a_key < b_key => {
                left.next();
                (a, a_key)
            }
            (Some(_), Some(_)) => {
                let (b_key, &b) = right.next().unwrap();
                (b, b_key)
            }
            (Some(&(a_key, &a)), None) => {
                left.next();
                (a, a_key)
            }
            (None, Some(&(b_key, &b))) => {
                right.next();
                (b, b_key)
            }
            (None, None) => unreachable!(),
        };
        loss += loss_fractional(model, q, obs.get(key).copied().unwrap_or(0.0))?;
    }
    ensure(loss.is_finite(), "nonfinite pair loss")?;
    Ok(loss)
}

/// Per-feature multiplicity variant of `merged_pair_loss`: each pair-merged
/// feature is charged against its own measured panel multiplicity entry.
fn merged_pair_loss_multiplicity(
    first: &Profile,
    second: &Profile,
    obs: &HashMap<FeatureKey, f64>,
    model: &ScoreModel,
    backgrounds: &FeatureBackgrounds,
) -> io::Result<f64> {
    let mut loss = 0.0;
    let mut left = first.iter().peekable();
    let mut right = second.iter().peekable();
    while left.peek().is_some() || right.peek().is_some() {
        let (q, key) = match (left.peek(), right.peek()) {
            (Some(&(a_key, &a)), Some(&(b_key, &b))) if a_key == b_key => {
                left.next();
                right.next();
                (a + b, a_key)
            }
            (Some(&(a_key, &a)), Some(&(b_key, _))) if a_key < b_key => {
                left.next();
                (a, a_key)
            }
            (Some(_), Some(_)) => {
                let (b_key, &b) = right.next().unwrap();
                (b, b_key)
            }
            (Some(&(a_key, &a)), None) => {
                left.next();
                (a, a_key)
            }
            (None, Some(&(b_key, &b))) => {
                right.next();
                (b, b_key)
            }
            (None, None) => unreachable!(),
        };
        loss += loss_fractional_entry(
            model,
            q,
            obs.get(key).copied().unwrap_or(0.0),
            backgrounds.entry(key),
        )?;
    }
    ensure(loss.is_finite(), "nonfinite pair loss")?;
    Ok(loss)
}

/// The window-domain extension's owner-resolved charging (the
/// no-duplication red line, owner-ruled): a charge site's observed support
/// for a feature is the sum over the DISTINCT owning universe partitions of
/// the charged material's candidate rows, of those partitions' routed
/// shares. Every row identity has exactly one owning group (the panel
/// tiling is disjoint per sequence), so no row is charged twice and no
/// row's material is charged against a foreign group's share. When the
/// charged material's rows share one owner (every site of the
/// pre-extension model), the map is exactly that partition's routed share
/// and every value is bit-identical to the unextended model.
///
/// Owner -> universe-partition resolution: an owner below the component
/// locus count is a component-local axis slot (resolve through the
/// component locus -> universe partition map); an owner at or above it is
/// an appended extension partition (already a universe id).
pub(crate) fn owner_universe_partition(owner: u32, component_locus_to_partition: &[u32]) -> u32 {
    if (owner as usize) < component_locus_to_partition.len() {
        component_locus_to_partition[owner as usize]
    } else {
        owner
    }
}

static EMPTY_ROUTED_OBS: std::sync::LazyLock<HashMap<FeatureKey, f64>> =
    std::sync::LazyLock::new(HashMap::new);

/// A traversal's owning universe partition: every segment shares one owner
/// (mixed-owner split admission is excluded at admission).
pub(crate) fn traversal_owner(traversal: &genome::SpanningTraversal) -> u32 {
    traversal.segments[0].partition as u32
}

pub(crate) fn owner_routed_obs<'a>(
    routed_equal: &'a RoutedObs,
    owner: u32,
    component_locus_to_partition: &[u32],
) -> &'a HashMap<FeatureKey, f64> {
    routed_equal
        .get(&owner_universe_partition(owner, component_locus_to_partition))
        .unwrap_or(&EMPTY_ROUTED_OBS)
}

/// Diploid pair loss with owner-resolved observed sides: when both classes'
/// rows share one owning partition the merged single-map charge is
/// bit-identical to `merged_pair_loss`; across distinct owners the merged
/// feature count is charged against the SUM of the distinct owners' shares
/// (the boundary convention's form: the junction-straddling record mass
/// distributes across exactly the touched partitions).
fn merged_pair_loss_owned(
    first: &Profile,
    second: &Profile,
    obs_first: &HashMap<FeatureKey, f64>,
    obs_second: &HashMap<FeatureKey, f64>,
    owners_equal: bool,
    model: &ScoreModel,
) -> io::Result<f64> {
    if owners_equal {
        return merged_pair_loss(first, second, obs_first, model);
    }
    let mut loss = 0.0;
    let mut left = first.iter().peekable();
    let mut right = second.iter().peekable();
    while left.peek().is_some() || right.peek().is_some() {
        let (q, key) = match (left.peek(), right.peek()) {
            (Some(&(a_key, &a)), Some(&(b_key, &b))) if a_key == b_key => {
                left.next();
                right.next();
                (a + b, a_key)
            }
            (Some(&(a_key, &a)), Some(&(b_key, _))) if a_key < b_key => {
                left.next();
                (a, a_key)
            }
            (Some(_), Some(_)) => {
                let (b_key, &b) = right.next().unwrap();
                (b, b_key)
            }
            (Some(&(a_key, &a)), None) => {
                left.next();
                (a, a_key)
            }
            (None, Some(&(b_key, &b))) => {
                right.next();
                (b, b_key)
            }
            (None, None) => unreachable!(),
        };
        let observed = obs_first.get(key).copied().unwrap_or(0.0)
            + obs_second.get(key).copied().unwrap_or(0.0);
        loss += loss_fractional(model, q, observed)?;
    }
    ensure(loss.is_finite(), "nonfinite owned pair loss")?;
    Ok(loss)
}

/// Per-feature multiplicity variant of `merged_pair_loss_owned`.
fn merged_pair_loss_multiplicity_owned(
    first: &Profile,
    second: &Profile,
    obs_first: &HashMap<FeatureKey, f64>,
    obs_second: &HashMap<FeatureKey, f64>,
    owners_equal: bool,
    model: &ScoreModel,
    backgrounds: &FeatureBackgrounds,
) -> io::Result<f64> {
    if owners_equal {
        return merged_pair_loss_multiplicity(first, second, obs_first, model, backgrounds);
    }
    let mut loss = 0.0;
    let mut left = first.iter().peekable();
    let mut right = second.iter().peekable();
    while left.peek().is_some() || right.peek().is_some() {
        let (q, key) = match (left.peek(), right.peek()) {
            (Some(&(a_key, &a)), Some(&(b_key, &b))) if a_key == b_key => {
                left.next();
                right.next();
                (a + b, a_key)
            }
            (Some(&(a_key, &a)), Some(&(b_key, _))) if a_key < b_key => {
                left.next();
                (a, a_key)
            }
            (Some(_), Some(_)) => {
                let (b_key, &b) = right.next().unwrap();
                (b, b_key)
            }
            (Some(&(a_key, &a)), None) => {
                left.next();
                (a, a_key)
            }
            (None, Some(&(b_key, &b))) => {
                right.next();
                (b, b_key)
            }
            (None, None) => unreachable!(),
        };
        let observed = obs_first.get(key).copied().unwrap_or(0.0)
            + obs_second.get(key).copied().unwrap_or(0.0);
        loss += loss_fractional_entry(model, q, observed, backgrounds.entry(key))?;
    }
    ensure(loss.is_finite(), "nonfinite owned pair loss")?;
    Ok(loss)
}

/// Boundary seam loss with owner-resolved observed sides: the sum of the
/// two adjacent candidates' OWNING partitions' shares (same form as the
/// pre-extension boundary); when both sides share one owner the record
/// mass touches that partition once, so the single map is charged alone.
fn profile_loss_boundary_owned(
    profile: &Profile,
    obs_left: &HashMap<FeatureKey, f64>,
    obs_right: &HashMap<FeatureKey, f64>,
    owners_equal: bool,
    model: &ScoreModel,
) -> io::Result<f64> {
    if owners_equal {
        return profile_loss_boundary(profile, obs_left, obs_left, model);
    }
    profile_loss_boundary(profile, obs_left, obs_right, model)
}

/// Per-feature multiplicity variant of `profile_loss_boundary_owned`.
fn profile_loss_boundary_multiplicity_owned(
    profile: &Profile,
    obs_left: &HashMap<FeatureKey, f64>,
    obs_right: &HashMap<FeatureKey, f64>,
    owners_equal: bool,
    model: &ScoreModel,
    backgrounds: &FeatureBackgrounds,
) -> io::Result<f64> {
    if owners_equal {
        return profile_loss_boundary_multiplicity(profile, obs_left, obs_left, model, backgrounds);
    }
    profile_loss_boundary_multiplicity(profile, obs_left, obs_right, model, backgrounds)
}

#[derive(Serialize)]
struct DpFinalist {
    proposal_loss: f64,
    m1_oracle_rescore: f64,
    pooled_external_rescore: f64,
    route: Vec<[Vec<serde_json::Value>; 2]>,
    selected: bool,
}

struct DpState {
    score: f64,
    pair: [usize; 2],
    pred: Option<u32>,
    /// Conflict-window history INCLUDING the current pair.
    history: Vec<[usize; 2]>,
    /// Fingerprint of `history[1..]` (the next state's dedup tail).
    suffix_hash: (u64, u64),
}

fn history_fingerprint(history: &[[usize; 2]]) -> (u64, u64) {
    let mut first = 0xcbf29ce484222325u64;
    let mut second = 0x9e3779b97f4a7c15u64;
    for pair in history {
        for &value in pair {
            first ^= value as u64;
            first = first.wrapping_mul(0x100000001b3);
            second = (second ^ (value as u64)).wrapping_mul(0x2545f4914f6cdd1d);
        }
    }
    (first, second)
}

/// Total-order wrapper for f64 scores under `total_cmp` (port of the
/// class-graph DP's `ScoreKey` in `genome.rs`).
#[derive(Clone, Copy)]
struct ScoreKey(f64);

impl PartialEq for ScoreKey {
    fn eq(&self, other: &Self) -> bool {
        self.0.total_cmp(&other.0).is_eq()
    }
}
impl Eq for ScoreKey {}
impl PartialOrd for ScoreKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for ScoreKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

/// Routed beam state key: the diploid allele pair plus the predecessor
/// history-tail fingerprint (the old materialize-then-sort prune's dedup
/// key, carried verbatim into the ported bounded layer).
type LayerKey = ([usize; 2], (u64, u64));

struct LayerCandidate {
    score: f64,
    pred: u32,
}

/// Bounded beam layer, ported from the class-graph DP's `BoundedSuccessors`
/// (`genome.rs`): an incrementally maintained top-`limit` set keyed by the
/// routed state key, ranked by (score, key) — exactly the kept-set order of
/// the materialize-then-sort prune it replaces (score, then natural key;
/// full-key ties always fell to the first-generated candidate, which the
/// existing-key Equal branch reproduces). Bulk drops and the cached worst
/// entry make the per-group cutoff test O(1); drop accounting follows the
/// slow path's incremental semantics (cutoff at drop time, ties
/// reclassified when the cutoff improves) instead of the old prune's
/// post-hoc sorted-stream classification — an accounting-only difference,
/// documented in the run report.
struct BoundedLayer {
    limit: usize,
    candidates: HashMap<LayerKey, LayerCandidate>,
    ranking: BTreeSet<(ScoreKey, LayerKey)>,
    worst: Option<(ScoreKey, LayerKey)>,
    dropped_ties: u64,
    dropped_non_ties: u64,
    insert_attempts: u64,
}

impl BoundedLayer {
    fn new(limit: usize) -> Self {
        Self {
            limit,
            candidates: HashMap::with_capacity(limit * 2),
            ranking: BTreeSet::new(),
            worst: None,
            dropped_ties: 0,
            dropped_non_ties: 0,
            insert_attempts: 0,
        }
    }

    fn is_full(&self) -> bool {
        self.candidates.len() >= self.limit
    }

    /// Current worst kept entry (score, key) — defined exactly when the
    /// layer is full.
    fn worst_entry(&self) -> (ScoreKey, LayerKey) {
        self.worst.expect("full bounded layer")
    }

    fn refresh_cutoff(&mut self, previous: Option<ScoreKey>) {
        if let (Some(previous), Some(current)) = (previous, self.worst) {
            if current.0.cmp(&previous) == std::cmp::Ordering::Less {
                self.dropped_non_ties = self.dropped_non_ties.saturating_add(self.dropped_ties);
                self.dropped_ties = 0;
            }
        }
    }

    fn record_bulk_drop(&mut self, score: f64, count: u64) {
        if self
            .worst
            .is_some_and(|worst| score.total_cmp(&worst.0 .0).is_eq())
        {
            self.dropped_ties = self.dropped_ties.saturating_add(count);
        } else {
            self.dropped_non_ties = self.dropped_non_ties.saturating_add(count);
        }
    }

    fn insert(&mut self, key: LayerKey, score: f64, pred: u32) {
        self.insert_attempts += 1;
        let previous_worst = self.worst.map(|(score, _)| score);
        if let Some(existing) = self.candidates.get_mut(&key) {
            if score.total_cmp(&existing.score) == std::cmp::Ordering::Less {
                self.ranking.remove(&(ScoreKey(existing.score), key));
                existing.score = score;
                existing.pred = pred;
                self.ranking.insert((ScoreKey(score), key));
                self.worst = self.ranking.last().copied();
                self.refresh_cutoff(previous_worst);
            }
            // Equal keeps the first-generated candidate (the old prune's
            // stable-sort dedup); a later worse copy is ignored.
            return;
        }
        if self.candidates.len() < self.limit {
            self.ranking.insert((ScoreKey(score), key));
            self.candidates.insert(key, LayerCandidate { score, pred });
            self.worst = self.ranking.last().copied();
            return;
        }
        let worst = self.worst.expect("full bounded layer");
        if (ScoreKey(score), key).cmp(&worst) == std::cmp::Ordering::Less {
            self.ranking.remove(&worst);
            self.candidates.remove(&worst.1);
            self.ranking.insert((ScoreKey(score), key));
            self.candidates.insert(key, LayerCandidate { score, pred });
            self.worst = self.ranking.last().copied();
            self.refresh_cutoff(Some(worst.0));
            self.record_bulk_drop(worst.0 .0, 1);
        } else {
            self.record_bulk_drop(score, 1);
        }
    }
}

struct BeamOutcome {
    layer: Vec<DpState>,
    dropped_ties: u64,
    dropped_non_ties: u64,
    bound_pruned_states: u64,
    insert_attempts: u64,
}

impl BoundedLayer {
    /// Finish the layer: emit the kept states in (score, key) order — the
    /// old prune's layer order — with the same per-state admissible-bound
    /// rejection before each push.
    fn finish(
        mut self,
        layer: &[DpState],
        conflict_window: usize,
        future_bound: f64,
        incumbent: f64,
    ) -> io::Result<BeamOutcome> {
        let mut next: Vec<DpState> = Vec::with_capacity(self.candidates.len());
        let mut bound_pruned = 0u64;
        for (rank, key) in self.ranking {
            let Some(candidate) = self.candidates.get(&key) else {
                continue;
            };
            let score = candidate.score;
            let pred = candidate.pred;
            if score + future_bound > incumbent + BOUND_PRUNE_EPSILON {
                bound_pruned += 1;
                continue;
            }
            let mut history = if pred == u32::MAX {
                Vec::new()
            } else {
                layer[pred as usize].history.clone()
            };
            history.push(key.0);
            if history.len() > conflict_window {
                let excess = history.len() - conflict_window;
                history.drain(0..excess);
            }
            let suffix_hash = history_fingerprint(&history[1..]);
            next.push(DpState {
                score,
                pair: key.0,
                pred: (pred != u32::MAX).then_some(pred),
                history,
                suffix_hash,
            });
        }
        Ok(BeamOutcome {
            layer: next,
            dropped_ties: self.dropped_ties,
            dropped_non_ties: self.dropped_non_ties,
            bound_pruned_states: bound_pruned,
            insert_attempts: self.insert_attempts,
        })
    }
}

/// Per-locus folded class machinery: the locus's feature universe interned
/// to dense fids (fid = rank in feature-sorted order, so a fid-ordered merge
/// reproduces the profile merge's summation order BIT-EXACTLY), each class's
/// prepared (fid, q, observed) array, and the per-fid maxima/observations
/// the admissible bound family consumes. This is the routed adaptation of
/// the pooled DP's converted-profile/feature-universe family (`genome.rs`):
/// under M1 the per-class-pair loss constant is the WHOLE pair loss (every
/// feature is fresh at every state — the pooled model's telescoping
/// reduces to per-pair constants, the documented Milestone 0.5 consequence),
/// so the interner's job here is to make each constant's streamed merge
/// dense (u32 fid compares) instead of BTreeMap/Vec-key comparisons, and to
/// memoize constants per (locus, class pair) so only reached pairs are ever
/// computed.
struct LocusFolded {
    /// (fid, q, observed) per class, fid-ascending. `observed` is the
    /// class's OWNING partition's routed share for the feature (the
    /// window-domain extension's owner-resolved charging; every class of the
    /// pre-extension model shares the locus's axis partition, so the values
    /// are bit-identical there).
    prepared: Vec<Vec<(u32, u64, f64)>>,
    /// Per-class owning universe partition (component-universe encoding;
    /// see `owner_universe_partition`).
    owners: Vec<u32>,
    /// True iff every class's rows share one owning partition — the
    /// pre-extension shape, where the per-fid loss tables and the merged
    /// single-map pair charge stay bit-identical.
    single_owner: bool,
    /// Per-fid maximum class observed share (admissible-bound input; the
    /// largest credit any class's charge can see for the feature).
    obs_by_fid: Vec<f64>,
    /// Per-fid maximum class count at this locus (admissible-bound input).
    max_add: Vec<u64>,
    /// Per-fid measured panel multiplicity entries (the attribution model);
    /// EMPTY under the quarantined flat path, where every fid falls back to
    /// the flat loss bit-identically.
    entries_by_fid: Vec<Option<FeatureBackgroundEntry>>,
}

impl LocusFolded {
    fn build(
        classes: &[Profile],
        owners: &[u32],
        routed_equal: &RoutedObs,
        component_locus_to_partition: &[u32],
    ) -> io::Result<Self> {
        Self::build_inner(
            classes,
            owners,
            routed_equal,
            component_locus_to_partition,
            None,
        )
    }

    /// The multiplicity build: per-fid measured panel entries (the interior
    /// charging path's local tables).
    fn build_backgrounds(
        classes: &[Profile],
        owners: &[u32],
        routed_equal: &RoutedObs,
        component_locus_to_partition: &[u32],
        backgrounds: &FeatureBackgrounds,
    ) -> io::Result<Self> {
        Self::build_inner(
            classes,
            owners,
            routed_equal,
            component_locus_to_partition,
            Some(backgrounds),
        )
    }

    fn build_inner(
        classes: &[Profile],
        owners: &[u32],
        routed_equal: &RoutedObs,
        component_locus_to_partition: &[u32],
        backgrounds: Option<&FeatureBackgrounds>,
    ) -> io::Result<Self> {
        ensure(
            classes.len() == owners.len(),
            "folded class/owner cardinality mismatch",
        )?;
        let single_owner = owners.iter().all(|&owner| owner == owners.first().copied().unwrap_or(0));
        let class_obs = owners
            .iter()
            .map(|&owner| owner_routed_obs(routed_equal, owner, component_locus_to_partition))
            .collect::<Vec<_>>();
        let mut interned: HashMap<&FeatureKey, u32> = HashMap::new();
        let mut order: Vec<&FeatureKey> = Vec::new();
        for profile in classes {
            for (feature, _) in profile.iter() {
                if let std::collections::hash_map::Entry::Vacant(entry) = interned.entry(feature)
                {
                    entry.insert(u32::MAX);
                    order.push(feature);
                }
            }
        }
        // fid = rank in feature-sorted order: the merged constant's f64
        // accumulation then visits features in exactly the profile order the
        // old BTreeMap merge used, so every lazily computed value is
        // bit-identical to the old fully precomputed table.
        order.sort();
        for (rank, feature) in order.iter().enumerate() {
            interned.insert(*feature, rank as u32);
        }
        let mut prepared = Vec::with_capacity(classes.len());
        for (profile, obs) in classes.iter().zip(class_obs.iter()) {
            let mut entries = Vec::with_capacity(profile.len());
            for (feature, &q) in profile.iter() {
                let fid = interned
                    .get(feature)
                    .copied()
                    .expect("interned locus feature");
                entries.push((fid, q, obs.get(feature).copied().unwrap_or(0.0)));
            }
            prepared.push(entries);
        }
        // Per-fid maximum observed share across the locus's classes (the
        // single-owner case reproduces the unextended per-fid map value
        // bit-for-bit).
        let mut obs_by_fid = vec![0.0f64; order.len()];
        for entries in &prepared {
            for &(fid, _, observed) in entries {
                if observed > obs_by_fid[fid as usize] {
                    obs_by_fid[fid as usize] = observed;
                }
            }
        }
        let entries_by_fid = match backgrounds {
            Some(backgrounds) => order
                .iter()
                .map(|feature| backgrounds.entry(feature))
                .collect::<Vec<_>>(),
            None => Vec::new(),
        };
        let mut max_add = vec![0u64; order.len()];
        for entries in &prepared {
            for &(fid, q, _) in entries {
                max_add[fid as usize] = max_add[fid as usize].max(q);
            }
        }
        Ok(Self {
            prepared,
            owners: owners.to_vec(),
            single_owner,
            max_add,
            obs_by_fid,
            entries_by_fid,
        })
    }

    #[inline]
    fn entry(&self, fid: u32) -> Option<FeatureBackgroundEntry> {
        self.entries_by_fid.get(fid as usize).copied().flatten()
    }

    /// Admissible per-class weight bound: b(c) = Σ_f min over the achievable
    /// merged range [q_c(f), q_c(f) + max_f] of the fractional loss — a lower
    /// bound on the pair loss of (c, ANY partner class at this locus), so
    /// max(b(c1), b(c2)) lower-bounds every pair (c1, c2). Used only to
    /// order and pre-skip class pairs (the pooled DP's weight-table family);
    /// every skip comparison carries BOUND_PRUNE_EPSILON so f64 regrouping in
    /// the bound can never drop a pair the exact score would have kept.
    fn class_weight(&self, class: usize, model: &ScoreModel) -> io::Result<f64> {
        let mut total = 0.0f64;
        for &(fid, q, observed) in &self.prepared[class] {
            total += loss_fractional_entry(model, q, observed, self.entry(fid))?
                + feature_potential_entry(q, self.max_add[fid as usize], observed, model, self.entry(fid));
        }
        Ok(total)
    }

    /// Per-fid memoized loss terms, the routed port of the pooled DP's
    /// `component_terms` family: loss(q, observed) is a pure function of
    /// (fid, q), so every distinct count reachable by a single-side or merged
    /// entry (0..=2·max_f) is evaluated ONCE with the identical
    /// `loss_fractional` call the inline loop made — the stored f64 is
    /// bit-identical, and the merge's summation order is unchanged, so every
    /// pair constant keeps its exact bits. Tables are built for the
    /// narrow-count fids first under a fixed entry budget (deterministic
    /// ascending max order); wider fids keep the direct evaluation (empty
    /// table) so memory stays bounded.
    fn loss_tables(&self, model: &ScoreModel) -> io::Result<(Vec<Vec<f64>>, u64)> {
        const TABLE_ENTRY_BUDGET: u64 = 4 << 20;
        // The per-fid tables key a count at ONE fixed observed share; with
        // mixed class owners the observed varies per class, so the tables
        // stay empty and every term takes the direct (bit-exact) evaluation.
        if !self.single_owner {
            return Ok((vec![Vec::<f64>::new(); self.max_add.len()], 0));
        }
        let mut order: Vec<u32> = (0..self.max_add.len() as u32).collect();
        order.sort_by_key(|&fid| self.max_add[fid as usize]);
        let mut tables = vec![Vec::<f64>::new(); self.max_add.len()];
        let mut budget = TABLE_ENTRY_BUDGET;
        let mut entries = 0u64;
        for fid in order {
            let max_add = self.max_add[fid as usize];
            let size = 2 * max_add + 1;
            if size > budget {
                break;
            }
            budget -= size;
            entries += size;
            let observed = self.obs_by_fid[fid as usize];
            let entry = self.entry(fid);
            let mut table = Vec::with_capacity(size as usize);
            for count in 0..=2 * max_add {
                table.push(loss_fractional_entry(model, count, observed, entry)?);
            }
            tables[fid as usize] = table;
        }
        Ok((tables, entries))
    }
}

/// Memoized per-fid loss term: the identical f64 the inline
/// `loss_fractional` loop produced for (fid, q), read from the dense table
/// when one was budgeted for that fid, recomputed directly otherwise. The
/// per-fid multiplicity entry (None = the quarantined flat path) applies only
/// to the direct evaluation; a budgeted table already bakes it in.
#[inline]
fn folded_term(
    fid: u32,
    q: u64,
    observed: f64,
    entry: Option<FeatureBackgroundEntry>,
    loss_tables: &[Vec<f64>],
    model: &ScoreModel,
) -> io::Result<f64> {
    if let Some(&value) = loss_tables[fid as usize].get(q as usize) {
        Ok(value)
    } else {
        loss_fractional_entry(model, q, observed, entry)
    }
}

/// Streamed merge of two prepared fid arrays: the identical accumulation of
/// `prepared_pair_loss` (same feature order, same observed side per branch),
/// so a lazily computed constant is bit-identical to the old precomputed
/// table entry. `entries_by_fid` carries the per-fid measured panel
/// multiplicity entries (empty = the quarantined flat path).
fn prepared_pair_loss_folded(
    first: &[(u32, u64, f64)],
    second: &[(u32, u64, f64)],
    owners_equal: bool,
    entries_by_fid: &[Option<FeatureBackgroundEntry>],
    loss_tables: &[Vec<f64>],
    model: &ScoreModel,
) -> io::Result<f64> {
    let entry = |fid: u32| entries_by_fid.get(fid as usize).copied().flatten();
    let mut loss = 0.0;
    let (mut i, mut j) = (0usize, 0usize);
    while i < first.len() || j < second.len() {
        match (first.get(i), second.get(j)) {
            (Some(&(fa, aq, ao)), Some(&(fb, bq, bo))) => {
                if fa == fb {
                    if owners_equal {
                        loss += folded_term(fa, aq + bq, ao, entry(fa), loss_tables, model)?;
                    } else {
                        // Distinct owners: the two copies' observed sides
                        // differ, so the merged count charges the SUM of the
                        // distinct owners' shares (the boundary convention's
                        // form) by direct evaluation.
                        loss += loss_fractional_entry(
                            model,
                            aq + bq,
                            ao + bo,
                            entry(fa),
                        )?;
                    }
                    i += 1;
                    j += 1;
                } else if fa < fb {
                    loss += folded_term(fa, aq, ao, entry(fa), loss_tables, model)?;
                    i += 1;
                } else {
                    loss += folded_term(fb, bq, bo, entry(fb), loss_tables, model)?;
                    j += 1;
                }
            }
            (Some(&(_, aq, ao)), None) => {
                loss += folded_term(first[i].0, aq, ao, entry(first[i].0), loss_tables, model)?;
                i += 1;
            }
            (None, Some(&(_, bq, bo))) => {
                loss += folded_term(second[j].0, bq, bo, entry(second[j].0), loss_tables, model)?;
                j += 1;
            }
            (None, None) => unreachable!(),
        }
    }
    ensure(loss.is_finite(), "nonfinite folded pair loss")?;
    Ok(loss)
}

/// Memoized per-locus class-pair loss constants: entries are computed on
/// first query by the SAME streamed merge the old full precompute used, so
/// every value is bit-identical; only the compute schedule changes (queried
/// pairs only, replacing the every-pair precompute that dominated the old
/// pair_tables stage).
struct LazyPairTable<'a> {
    classes: usize,
    prepared: &'a [Vec<(u32, u64, f64)>],
    owners: &'a [u32],
    entries_by_fid: &'a [Option<FeatureBackgroundEntry>],
    loss_tables: &'a [Vec<f64>],
    values: Vec<f64>,
    computed: Vec<bool>,
    computed_pairs: u64,
    model: &'a ScoreModel,
}

impl<'a> LazyPairTable<'a> {
    fn new(
        classes: usize,
        prepared: &'a [Vec<(u32, u64, f64)>],
        owners: &'a [u32],
        entries_by_fid: &'a [Option<FeatureBackgroundEntry>],
        loss_tables: &'a [Vec<f64>],
        model: &'a ScoreModel,
    ) -> Self {
        let size = classes * (classes + 1) / 2;
        Self {
            classes,
            prepared,
            owners,
            entries_by_fid,
            loss_tables,
            values: vec![0.0; size],
            computed: vec![false; size],
            computed_pairs: 0,
            model,
        }
    }

    fn from_full(
        classes: usize,
        prepared: &'a [Vec<(u32, u64, f64)>],
        owners: &'a [u32],
        entries_by_fid: &'a [Option<FeatureBackgroundEntry>],
        loss_tables: &'a [Vec<f64>],
        values: Vec<f64>,
        model: &'a ScoreModel,
    ) -> io::Result<Self> {
        let size = classes * (classes + 1) / 2;
        ensure(
            values.len() == size,
            "precomputed pair table cardinality mismatch",
        )?;
        Ok(Self {
            classes,
            prepared,
            owners,
            entries_by_fid,
            loss_tables,
            values,
            computed: vec![true; size],
            computed_pairs: size as u64,
            model,
        })
    }

    fn entry(&mut self, first: usize, second: usize) -> io::Result<f64> {
        let (lo, hi) = (first.min(second), first.max(second));
        let index = hi * (hi + 1) / 2 + lo;
        if !self.computed[index] {
            self.values[index] = prepared_pair_loss_folded(
                &self.prepared[lo],
                &self.prepared[hi],
                self.owners[lo] == self.owners[hi],
                self.entries_by_fid,
                self.loss_tables,
                self.model,
            )?;
            self.computed[index] = true;
            self.computed_pairs += 1;
        }
        Ok(self.values[index])
    }
}

/// Minimum of the fractional-observed loss over the achievable count
/// interval [q, q + max_add], relative to loss(q). Port of the pooled DP's
/// `feature_potential` (`genome.rs`) with the routed model's fractional
/// observed side; used only for admissible bounds.
fn feature_potential_entry(
    q: u64,
    max_add: u64,
    observed: f64,
    model: &ScoreModel,
    entry: Option<FeatureBackgroundEntry>,
) -> f64 {
    let loss = |count: u64| -> f64 {
        loss_fractional_entry(model, count, observed, entry).expect("admissible-bound loss")
    };
    let hi = q.saturating_add(max_add);
    let mut best = loss(q);
    if hi != q {
        best = best.min(loss(hi));
    }
    let per_count = model.histogram as f64 * model.depth / model.denominator;
    // The interior optimum of the per-variant loss: attribution
    // s* = w_f·C - background; gentle-beta s* = C - beta_f; flat
    // s* = C - background.
    let optimum = match entry {
        None => (observed - model.background) / per_count,
        Some(entry) => match multiplicity_variant() {
            MultiplicityVariant::Attribution => {
                (observed * entry.weight - model.background) / per_count
            }
            MultiplicityVariant::GentleBeta => (observed - entry.beta) / per_count,
        },
    };
    if optimum.is_finite() && optimum > 0.0 {
        let floored = (optimum as u64).max(q).min(hi);
        if floored != q {
            best = best.min(loss(floored));
        }
        let ceiled = (optimum.ceil() as u64).max(q).min(hi);
        if ceiled != floored {
            best = best.min(loss(ceiled));
        }
    }
    best - loss(q)
}

/// Admissible lower bound on the per-locus minimum pair loss (the pooled
/// DP's fresh-potential family, adapted to the additive M1 objective) on the
/// folded per-fid arrays: every pair's loss is at least the per-feature
/// minimum over the achievable merged range [0, 2·max_f], summed over the
/// locus's full feature universe — a superset of any pair's features, and
/// each c_f(0) <= 0, so the sum lower-bounds every pair. The old
/// implementation used the exact minimum of the fully precomputed pair
/// table; this port replaces that exact minimum (which required the
/// every-pair precompute) with the admissible bound — bound_pruned counters
/// change accordingly, kept-state survivors can only grow, and the
/// selection comparison against the validated c12/c14 run is the acceptance
/// test.
fn locus_min_pair_lower_bound_folded(folded: &LocusFolded, model: &ScoreModel) -> f64 {
    let mut total = 0.0f64;
    for fid in 0..folded.max_add.len() {
        total += feature_potential_entry(
            0,
            2 * folded.max_add[fid],
            folded.obs_by_fid[fid],
            model,
            folded.entry(fid as u32),
        );
    }
    total
}

/// Commutative content hash of one (feature, count) profile entry, ported
/// from the pooled DP's `profile_value_fingerprint` (`genome.rs`); XOR
/// accumulation over entries makes the class-signature fold
/// order-independent.
fn profile_entry_hash(feature: &FeatureKey, value: u64) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    let mut fnv_update = |value: &mut u64, bytes: &[u8]| {
        for &byte in bytes {
            *value ^= byte as u64;
            *value = value.wrapping_mul(0x100000001b3);
        }
    };
    let len = feature.len() as u64;
    fnv_update(&mut hash, &len.to_le_bytes());
    for token in feature {
        fnv_update(&mut hash, &token.to_le_bytes());
    }
    fnv_update(&mut hash, &value.to_le_bytes());
    hash
}

/// Hash-first class signature: the commutative XOR fold of the profile's
/// entry hashes. Classing verifies exactly on hash hit, so the partition
/// is exact (hash collisions only cost a comparison).
fn profile_signature_hash(profile: &Profile) -> u64 {
    let mut acc = 0u64;
    for (feature, &count) in profile.iter() {
        acc ^= profile_entry_hash(feature, count);
    }
    acc
}


/// Same-source overlap conflict window, mirroring `genome::audit_chain`'s
/// derivation without the per-seam profile-signature sweep.
fn conflict_window_of(ranges: &[Vec<genome::SpanningTraversal>]) -> usize {
    let mut by_source: BTreeMap<usize, Vec<(usize, u64, u64)>> = BTreeMap::new();
    for (locus, locus_ranges) in ranges.iter().enumerate() {
        for traversal in locus_ranges {
            for segment in &traversal.segments {
                if segment.start < segment.end {
                    by_source
                        .entry(segment.source)
                        .or_default()
                        .push((locus, segment.start, segment.end));
                }
            }
        }
    }
    let mut window = 0usize;
    for spans in by_source.values_mut() {
        spans.sort_by_key(|&(locus, start, end)| (start, end, locus));
        for left in 0..spans.len() {
            for right in left + 1..spans.len() {
                if spans[right].1 >= spans[left].2 {
                    break;
                }
                if spans[left].1 < spans[right].2 {
                    window = window.max(spans[left].0.abs_diff(spans[right].0));
                }
            }
        }
    }
    window
}

/// Geometric profile of one split partial (sub-interval of a locus parent
/// interval), swept from the parent's self anchors and RC-view records.
/// Window-relevant anchor decisions are context-identical within the parent's
/// padded extraction (an anchor a partial window can contain lies at least
/// READ_LENGTH - k inside the parent's core, so its syncmer decision context
/// is complete in the parent's extraction), and clipping a maximal record per
/// window preserves the per-window record set, so the partial's profile is
/// exactly the parent-data sweep over the partial's window grid.
fn partial_profile_from_parents(
    parent_data: &HashMap<(usize, u64, u64), (Vec<(i32, u64)>, Vec<Vec<(i32, u64)>>)>,
    key: (usize, u64, u64),
    k: u64,
) -> io::Result<Profile> {
    let (source, pstart, pend) = key;
    let mut chosen: Option<(&(usize, u64, u64), u64)> = None;
    for parent in parent_data.keys() {
        if parent.0 == source && parent.1 <= pstart && pend <= parent.2 {
            let slack = (parent.2 - parent.1) - (pend - pstart);
            if chosen.is_none_or(|(_, best)| slack < best) {
                chosen = Some((parent, slack));
            }
        }
    }
    let (parent, _) = chosen.ok_or_else(|| invalid("split partial lacks a containing parent"))?;
    let origin = pstart - parent.1;
    let (anchors, records) = &parent_data[parent];
    let shift = |pos: u64| pos.checked_sub(origin);
    let self_anchors: Vec<(i32, u64)> = anchors
        .iter()
        .filter_map(|&(node, pos)| shift(pos).map(|q| (node, q)))
        .collect();
    let ir: Vec<Vec<(i32, u64)>> = records
        .iter()
        .map(|record| {
            record
                .iter()
                .filter_map(|&(node, pos)| shift(pos).map(|q| (node, q)))
                .collect()
        })
        .collect();
    geometric_sweep_profile(
        &self_anchors,
        &ir,
        pend - pstart,
        k,
        &format!("partial:{source}:{pstart}-{pend}"),
    )
}

/// Oracle (slow-path) profile of one allele, for the M1 external rescore,
/// the reference scorer, and (owner ruling, dp-evaluator-consistency) the
/// haploid DP's per-window losses: `profile_event_interior` over singles,
/// partial interiors plus the surviving L149 internal-seam machinery over
/// splits.
/// The per-feature spelled PATH spans of one DP candidate allele (the
/// instance-level exoneration's per-allele spelled side for the coupling
/// channel): the union of the segments' piece spans and each consecutive
/// pair's seam spans — the SAME shape `oracle_allele_profile` charges with
/// (segments' oracle profiles plus every interior pair's seam), positions
/// only. Memoized per traversal identity; the per-segment spans share the
/// orientation-keyed piece memo.
pub(crate) fn oracle_allele_spans(
    panel: &SyngIndex,
    sources: &routes::Sources,
    path_of_source: &[usize],
    traversal: &genome::SpanningTraversal,
    memo: &mut HashMap<String, HashMap<FeatureKey, Vec<SpelledSpan>>>,
    piece_spans_memo: &mut SpelledMemo,
) -> io::Result<HashMap<FeatureKey, Vec<SpelledSpan>>> {
    if let Some(spans) = memo.get(&traversal.identity) {
        return Ok(spans.clone());
    }
    let mut out: HashMap<FeatureKey, Vec<SpelledSpan>> = HashMap::new();
    let mut piece_spans: Vec<(usize, u64, u64, bool, std::rc::Rc<BTreeMap<FeatureKey, Vec<(u64, u64)>>>)> = Vec::new();
    for segment in &traversal.segments {
        if segment.start >= segment.end {
            continue;
        }
        let piece = (segment.source, segment.start, segment.end, segment.reverse);
        let spans = match piece_spans_memo.get(&piece) {
            Some(spans) => spans.clone(),
            None => {
                let mut sequence = sources.fetch(segment.source, segment.start, segment.end)?;
                if segment.reverse {
                    sequence = impg::graph::reverse_complement(&sequence);
                }
                let spans = std::rc::Rc::new(oracle_segment_spans(panel, &sequence)?);
                piece_spans_memo.insert(piece, spans.clone());
                spans
            }
        };
        piece_spans.push((segment.source, segment.start, segment.end, segment.reverse, spans));
    }
    for (source, start, end, reverse, spans) in &piece_spans {
        piece_spans_to_path_from(path_of_source, (*source, *start, *end, *reverse), spans, &mut out)?;
    }
    // The interior seams (every consecutive pair — the profile's own
    // convention): the seam spans split at the boundary; each side maps
    // through its own piece.
    for pair in piece_spans.windows(2) {
        let left = (pair[0].0, pair[0].1, pair[0].2, pair[0].3);
        let right = (pair[1].0, pair[1].1, pair[1].2, pair[1].3);
        let left_seq = sources.fetch(left.0, left.1, left.2)?;
        let left_seq = if left.3 { impg::graph::reverse_complement(&left_seq) } else { left_seq };
        let right_seq = sources.fetch(right.0, right.1, right.2)?;
        let right_seq = if right.3 { impg::graph::reverse_complement(&right_seq) } else { right_seq };
        let seam_spans = profile_event_seam_spans(panel, &left_seq, &right_seq, READ_LENGTH)?;
        let left_len = left_seq.len() as u64;
        let left_take = left_len.min(READ_LENGTH as u64 - 1);
        for (feature, list) in &seam_spans {
            let entry = out.entry(feature.clone()).or_default();
            for &(is_left, lo, hi) in list {
                let (piece, off_lo, off_hi) = if is_left {
                    let base = left_len - left_take;
                    (left, base + lo, base + hi)
                } else {
                    (right, lo - left_take, hi - left_take)
                };
                let (plo, phi) = piece_span_coords(piece, off_lo, off_hi);
                let path = *path_of_source
                    .get(piece.0)
                    .ok_or_else(|| invalid("piece source absent from the path map"))?;
                entry.push(SpelledSpan { path, lo: plo, hi: phi });
            }
        }
    }
    memo.insert(traversal.identity.clone(), out.clone());
    Ok(out)
}

/// The piece spans -> path mapping for externally built piece-span tables
/// (the allele spans' own pieces; same mapping as the chain-level helper).
fn piece_spans_to_path_from(
    path_of_source: &[usize],
    piece: (usize, u64, u64, bool),
    spans: &BTreeMap<FeatureKey, Vec<(u64, u64)>>,
    out: &mut HashMap<FeatureKey, Vec<SpelledSpan>>,
) -> io::Result<()> {
    let path = *path_of_source
        .get(piece.0)
        .ok_or_else(|| invalid("piece source absent from the path map"))?;
    for (feature, list) in spans {
        let entry = out.entry(feature.clone()).or_default();
        for &(lo, hi) in list {
            let (plo, phi) = piece_span_coords(piece, lo, hi);
            entry.push(SpelledSpan { path, lo: plo, hi: phi });
        }
    }
    Ok(())
}

pub(crate) fn oracle_allele_profile(
    panel: &SyngIndex,
    sources: &routes::Sources,
    traversal: &genome::SpanningTraversal,
    memo: &mut HashMap<String, Profile>,
) -> io::Result<Profile> {
    if let Some(profile) = memo.get(&traversal.identity) {
        return Ok(profile.clone());
    }
    let profile = if traversal.segments.len() == 1 {
        let sequence = allele_sequence(sources, &traversal.segments)?;
        oracle_segment_profile(panel, &sequence)?
    } else {
        // The established same-owner seam convention at every interior seam:
        // the segments' oracle profiles plus each consecutive pair's
        // profile_event_seam juxtaposition (the co-occurring/pooled form —
        // admitted stitched chains are adjacent same-source rows, so every
        // interior seam is the panel-joined real junction; the split
        // candidates' novel-restricted handling lives in the classing, not
        // here — the oracle profile is the rescore's own charged shape).
        let mut parts: Vec<Profile> = Vec::with_capacity(traversal.segments.len() + 1);
        let mut sequences: Vec<Vec<u8>> = Vec::with_capacity(traversal.segments.len());
        for segment in &traversal.segments {
            let sequence = allele_sequence(sources, std::slice::from_ref(segment))?;
            parts.push(oracle_segment_profile(panel, &sequence)?);
            sequences.push(sequence);
        }
        for pair in sequences.windows(2) {
            let (seam, _) =
                genome::profile_event_seam(panel, &pair[0], &pair[1], READ_LENGTH, MAX_FEATURES)?;
            parts.push(seam);
        }
        merge_profiles(&parts.iter().collect::<Vec<&Profile>>())?
    };
    memo.insert(traversal.identity.clone(), profile.clone());
    Ok(profile)
}

/// Per-run memo for orientation-resolved L149 flanks, keyed by
/// (source, strand, crop start, crop end) — the pooled boundary-fetch fix
/// family: endpoint flanks repeat heavily across candidates, loci and the
/// split-junction/rescore consumers, and the base `Sources` memo keys only
/// forward-strand crops, so the reverse-complement and the clone otherwise
/// repeat on every request.
type FlankMemo = std::sync::Mutex<HashMap<(usize, bool, u64, u64), Vec<u8>>>;

/// Orientation-resolved source fetch through the flank memo: the stored
/// value is the post-reverse-complement crop, so repeated requests for the
/// same (source, strand, position) flank are a single lookup.
fn fetch_oriented(
    sources: &routes::Sources,
    memo: &FlankMemo,
    source: usize,
    reverse: bool,
    start: u64,
    end: u64,
) -> io::Result<Vec<u8>> {
    let key = (source, reverse, start, end);
    if let Ok(memo) = memo.lock() {
        if let Some(bytes) = memo.get(&key) {
            return Ok(bytes.clone());
        }
    }
    let mut part = sources.fetch(source, start, end)?;
    if reverse {
        part = impg::graph::reverse_complement(&part);
    }
    if let Ok(mut memo) = memo.lock() {
        memo.insert(key, part.clone());
    }
    Ok(part)
}

/// L149 flank of one segment (head or tail), fetching only the needed crop:
/// `profile_event_seam` consumes at most the first/last READ_LENGTH-1 bases
/// of each side, so split internal seams never need the full segments.
fn segment_flank(
    sources: &routes::Sources,
    memo: &FlankMemo,
    segment: &SourceRange,
    head: bool,
    flank: usize,
) -> io::Result<Vec<u8>> {
    if segment.start == segment.end {
        return Ok(Vec::new());
    }
    let length = segment.end - segment.start;
    let take = (flank as u64).min(length);
    let (start, end) = if segment.reverse == head {
        (segment.end - take, segment.end)
    } else {
        (segment.start, segment.start + take)
    };
    fetch_oriented(sources, memo, segment.source, segment.reverse, start, end)
}

/// Loss of a diploid pair of prepared class vectors (feature-sorted, with
/// each feature's routed share resolved once per class), merged without
/// per-pair hash lookups.
fn prepared_pair_loss(
    first: &[(&FeatureKey, u64, f64)],
    second: &[(&FeatureKey, u64, f64)],
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = 0.0;
    let (mut i, mut j) = (0usize, 0usize);
    while i < first.len() || j < second.len() {
        let take_first = match (first.get(i), second.get(j)) {
            (Some(&(_, aq, ao)), Some(&(_, _, _))) if j >= second.len() => true,
            (Some(&(_, aq, ao)), Some(&(_, bq, bo))) => {
                if first[i].0 == second[j].0 {
                    i += 1;
                    j += 1;
                    loss += loss_fractional(model, aq + bq, ao)?;
                    continue;
                }
                first[i].0 < second[j].0
            }
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => unreachable!(),
        };
        if take_first {
            let (_, aq, ao) = first[i];
            i += 1;
            loss += loss_fractional(model, aq, ao)?;
        } else {
            let (_, bq, bo) = second[j];
            j += 1;
            loss += loss_fractional(model, bq, bo)?;
        }
    }
    ensure(loss.is_finite(), "nonfinite prepared pair loss")?;
    Ok(loss)
}

#[allow(clippy::too_many_arguments)]
fn run_routed_dp(
    panel: &SyngIndex,
    sources: &routes::Sources,
    graph: &routes::Graph,
    routes_dir: &PathBuf,
    path_of_source: &[usize],
    axis_slice: &[genome::AxisInterval],
    mut ranges: Vec<Vec<genome::SpanningTraversal>>,
    territory: &[Vec<SourceRange>],
    k: u64,
    routed_equal: &RoutedObs,
    component_locus_to_partition: &[u32],
    locus_offset: usize,
    component_loci: usize,
    target_source: usize,
    target_length: u64,
    model: &ScoreModel,
    sample_counts: &impg::sample_mem_bwt::WeightedBwt,
    evidence_path: &PathBuf,
    ranks_path: Option<&std::path::Path>,
    refinement_margin: f64,
    refinement_top_k: usize,
    beam_width: usize,
    tie_rotations: usize,
    finalist_count: usize,
    reference_specs: &[(String, PathBuf, PathBuf)],
    // The record-once site-map builder (Fix 1's mixed-owner charging) and
    // the windows' full observed profiles (Fix 2's omission universes; the
    // DP's reference ladder charges the haploid references through them).
    site: &SiteObserved,
    window_obs: &[HashMap<FeatureKey, f64>],
    rescore_directory: &std::path::Path,
    rss: &mut genome::PeriodicRssGuard,
) -> io::Result<serde_json::Value> {
    let started = Instant::now();
    let locus_count = ranges.len();
    ensure(locus_count > 0, "empty DP domain")?;
    ensure(
        axis_slice.len() == locus_count,
        "axis/DP domain cardinality mismatch",
    )?;

    // ------------------------------------------------ split admission (oracle parity)
    let admission_started = Instant::now();
    let evidence_text = std::fs::read_to_string(evidence_path)?;
    let evidence: serde_json::Value = serde_json::from_str(
        &evidence_text[evidence_text
            .find('{')
            .ok_or_else(|| invalid("refinement evidence has no JSON"))?..],
    )
    .map_err(io::Error::other)?;
    let local = evidence["accounting"]["local_pair_evidence"]
        .as_array()
        .ok_or_else(|| invalid("refinement evidence lacks local pair diagnostics"))?;
    let rank_checkpoint = match ranks_path {
        Some(path) => {
            let text = std::fs::read_to_string(path)?;
            let value: serde_json::Value = serde_json::from_str(
                &text[text
                    .find('{')
                    .ok_or_else(|| invalid("rank checkpoint has no JSON"))?..],
            )
            .map_err(io::Error::other)?;
            Some(
                serde_json::from_value::<Vec<genome::SplitPrefilterStats>>(
                    value["split_prefilter_stats"].clone(),
                )
                .map_err(io::Error::other)?,
            )
        }
        None => None,
    };
    let checkpoint = rank_checkpoint
        .as_ref()
        .ok_or_else(|| invalid("the routed DP requires the ranks checkpoint (oracle parity)"))?;
    let mut expand = vec![false; locus_count];
    let mut selected_indices = vec![Vec::<usize>::new(); locus_count];
    for row in local {
        let full_locus = row["locus"]
            .as_u64()
            .ok_or_else(|| invalid("invalid refinement locus"))?
            as usize;
        if let Some(locus) = full_locus
            .checked_sub(locus_offset)
            .filter(|&l| l < locus_count)
        {
            if row["native_minus_best"].as_f64().unwrap_or(0.0) >= refinement_margin {
                expand[locus] = true;
            }
        }
    }
    if let Some(copies) = evidence["final_result"]["selected_physical_alleles"].as_array() {
        for copy in copies {
            let alleles = copy
                .as_array()
                .ok_or_else(|| invalid("invalid selected refinement route"))?;
            for (full_locus, allele) in alleles.iter().enumerate() {
                if let Some(locus) = full_locus
                    .checked_sub(locus_offset)
                    .filter(|&l| l < locus_count)
                {
                    if let Some(identity) = allele["identity"].as_str() {
                        if let Some(index) = ranges[locus]
                            .iter()
                            .position(|candidate| candidate.identity == identity)
                        {
                            if !selected_indices[locus].contains(&index) {
                                selected_indices[locus].push(index);
                            }
                        }
                    }
                }
            }
            for boundary in 0..alleles.len().saturating_sub(1) {
                let left = alleles[boundary]["segments"][0]["source"].as_u64();
                let right = alleles[boundary + 1]["segments"][0]["source"].as_u64();
                if left != right {
                    for full_locus in boundary.saturating_sub(1)
                        ..=(boundary + 2).min(locus_count + locus_offset - 1)
                    {
                        if let Some(locus) = full_locus
                            .checked_sub(locus_offset)
                            .filter(|&l| l < locus_count)
                        {
                            expand[locus] = true;
                        }
                    }
                }
            }
        }
    }
    let mut ports = routes::Ports::open_without_global_verification(routes_dir, graph)?;
    let mut split_candidate_counts = vec![0usize; locus_count];
    for locus in 0..locus_count {
        if !expand[locus] {
            continue;
        }
        let full_locus = locus_offset + locus;
        let mut admitted: Vec<usize> = local
            .get(full_locus)
            .ok_or_else(|| invalid("refinement evidence lacks a slice locus"))?
            ["top_class_representatives"]
            .as_array()
            .ok_or_else(|| invalid("refinement evidence lacks top classes"))?
            .iter()
            .take(refinement_top_k)
            .filter_map(|value| value.as_u64().map(|value| value as usize))
            .filter(|&index| {
                index < ranges[locus].len()
                    && ranges[locus][index].segments.len() == 1
                    && ranges[locus][index].segments[0].start < ranges[locus][index].segments[0].end
            })
            .collect();
        for &selected in &selected_indices[locus] {
            let candidate = &ranges[locus][selected];
            if candidate.segments.len() == 1
                && candidate.segments[0].start < candidate.segments[0].end
                && !admitted.contains(&selected)
            {
                if admitted.len() == refinement_top_k {
                    admitted.pop();
                }
                admitted.push(selected);
            }
        }
        // Split admission stays within the locus's own axis partition: the
        // split machinery composes an allele's two pieces from the anchor
        // group's own candidate rows (public graph geometry). The
        // window-domain extension's added rows from other owning partitions
        // enter the domain as standalone single-segment candidates
        // (documented semantic choice; mixed-owner split charging is not
        // defined in this run).
        admitted.retain(|&index| {
            ranges[locus][index]
                .segments
                .iter()
                .all(|segment| segment.partition == full_locus)
        });
        let stats = checkpoint
            .get(full_locus)
            .ok_or_else(|| invalid("rank checkpoint locus mismatch"))?;
        let retained = stats
            .pairs
            .iter()
            .flat_map(|pair| {
                pair.retained
                    .iter()
                    .map(|(_, identity, _)| identity.as_str())
            })
            .collect::<BTreeSet<_>>();
        let split: Vec<genome::SpanningTraversal> = genome::split_candidates(
            full_locus,
            &ranges[locus],
            &admitted,
            graph,
            &mut ports,
        )?
        .into_iter()
        .filter(|candidate| retained.contains(candidate.identity.as_str()))
        .collect();
        ensure(
            split.len() <= stats.retained_cuts,
            "rank checkpoint regeneration produced extra candidates"
        )?;
        split_candidate_counts[locus] = split.len();
        ranges[locus].extend(split);
    }
    genome::retain_native_endpoint_candidates_in_component_range(
        &mut ranges,
        target_source,
        target_length,
        locus_offset,
        component_loci - 1,
    )?;
    let admission_seconds = admission_started.elapsed().as_secs_f64();
    rss_probe(rss, "dp_admission")?;
    eprintln!(
        "[dp] admission done: alleles {:?}, splits {:?}",
        ranges.iter().map(Vec::len).collect::<Vec<_>>(),
        split_candidate_counts
    );

    // Orientation-resolved L149 flank memo (P3): shared by the classing
    // junction keys and the boundary endpoint extraction, keyed by
    // (source, strand, crop), so repeated endpoint-flank fetches cost one
    // lookup instead of a base fetch + reverse-complement.
    let flank_memo: FlankMemo = std::sync::Mutex::new(HashMap::new());
    // ------------------------------------------------ per-locus geometric classing
    let profiling_started = Instant::now();
    let partition_obs: Vec<HashMap<FeatureKey, f64>> = (0..locus_count)
        .map(|locus| {
            let partition = component_locus_to_partition[locus_offset + locus];
            ensure(partition != u32::MAX, "DP locus outside the routing universe")?;
            Ok(routed_equal.get(&partition).cloned().unwrap_or_default())
        })
        .collect::<io::Result<_>>()?;
    struct LocusClasses {
        profiles: Vec<Profile>,
        membership: Vec<usize>,
        /// Per-class owning universe partition (the window-domain
        /// extension's owner-resolved charging; the owning partition enters
        /// the class key, so classes stay owner-homogeneous).
        class_owners: Vec<u32>,
    }
    let mut locus_classes: Vec<LocusClasses> = Vec::with_capacity(locus_count);
    let mut allele_counts = Vec::with_capacity(locus_count);
    let mut split_seam_queries = 0u64;
    let mut ir_extractions = 0u64;
    // Per-locus classing sub-stage accounting (P4 instrumentation): parents
    // (anchor + IR extraction), sweeps (parent profile construction),
    // partials+seams (split partials and L149 junction seams), and classing
    // itself, so the next measurement attributes the classing stage
    // precisely instead of reporting one undifferentiated number.
    let mut substage_parents = 0.0f64;
    let mut substage_sweeps = 0.0f64;
    let mut substage_partials_seams = 0.0f64;
    let mut substage_classing = 0.0f64;
    for (locus, locus_ranges) in ranges.iter().enumerate() {
        let parents_started = Instant::now();
        let mut parent_intervals: BTreeSet<(usize, u64, u64)> = BTreeSet::new();
        for traversal in locus_ranges {
            for segment in &traversal.segments {
                if segment.start < segment.end {
                    parent_intervals.insert((segment.source, segment.start, segment.end));
                }
            }
        }
        let parent_list: Vec<(usize, u64, u64)> = parent_intervals.into_iter().collect();
        let parent_data: HashMap<(usize, u64, u64), (Vec<(i32, u64)>, Vec<Vec<(i32, u64)>>)> =
            parent_list
                .par_iter()
                .map(|&(source, start, end)| {
                    let anchors =
                        contained_path_anchors(panel, path_of_source[source], start, end, k)?;
                    let records = geometric_ir_records(panel, sources, source, start, end, k)?;
                    Ok(((source, start, end), (anchors, records)))
                })
                .collect::<io::Result<_>>()?;
        ir_extractions += parent_data.len() as u64;
        let parents_seconds = parents_started.elapsed().as_secs_f64();
        let sweeps_started = Instant::now();
        let parent_profiles: HashMap<(usize, u64, u64), Profile> = parent_data
            .par_iter()
            .map(|(key, (anchors, records))| {
                let profile = geometric_sweep_profile(
                    anchors,
                    records,
                    key.2 - key.1,
                    k,
                    &format!("locus-{locus}-interval-{}-{}-{}", key.0, key.1, key.2),
                )?;
                Ok((*key, profile))
            })
            .collect::<io::Result<_>>()?;
        let sweeps_seconds = sweeps_started.elapsed().as_secs_f64();
        // Distinct split partials and junction L149 flank pairs, computed in
        // parallel once per locus (splits share partials and junctions).
        let partials_started = Instant::now();
        let mut partial_keys: BTreeSet<(usize, u64, u64)> = BTreeSet::new();
        let mut junctions: BTreeSet<(Vec<u8>, Vec<u8>)> = BTreeSet::new();
        for traversal in locus_ranges {
            if traversal.segments.len() == 2 {
                for segment in &traversal.segments {
                    if segment.start < segment.end {
                        partial_keys.insert((segment.source, segment.start, segment.end));
                    }
                }
                junctions.insert((
                    segment_flank(
                        sources,
                        &flank_memo,
                        &traversal.segments[0],
                        false,
                        READ_LENGTH - 1,
                    )?,
                    segment_flank(
                        sources,
                        &flank_memo,
                        &traversal.segments[1],
                        true,
                        READ_LENGTH - 1,
                    )?,
                ));
            }
        }
        let partials: HashMap<(usize, u64, u64), std::sync::Arc<Profile>> = partial_keys
            .par_iter()
            .map(|&key| {
                let profile = if let Some(full) = parent_profiles.get(&key) {
                    full.clone()
                } else {
                    partial_profile_from_parents(&parent_data, key, k)?
                };
                Ok((key, std::sync::Arc::new(profile)))
            })
            .collect::<io::Result<_>>()?;
        let seam_profiles_by_junction: HashMap<(Vec<u8>, Vec<u8>), std::sync::Arc<Profile>> =
            junctions
                .par_iter()
                .map(|junction| {
                    let (seam, _) = genome::profile_event_seam(
                        panel,
                        &junction.0,
                        &junction.1,
                        READ_LENGTH,
                        MAX_FEATURES,
                    )?;
                    Ok((junction.clone(), std::sync::Arc::new(seam)))
                })
                .collect::<io::Result<_>>()?;
        split_seam_queries += seam_profiles_by_junction.len() as u64;
        let partials_seconds = partials_started.elapsed().as_secs_f64();
        // Per-allele profile + commutative-hash class signature,
        // chunk-parallel; classed sequentially so only one chunk's profiles
        // are resident. Hash-first classing (P4): the signature is the 64-bit
        // XOR fold of per-entry fingerprints (`profile_signature_hash`), and
        // a hash hit is verified by exact profile comparison, so the class
        // partition stays exact — identical to the old full-signature
        // comparison (first-seen order wins in both) — while the per-allele
        // key-clone signature construction disappears.
        let classing_started = Instant::now();
        let mut class_profiles: Vec<Profile> = Vec::new();
        let mut class_owners: Vec<u32> = Vec::new();
        let mut class_buckets: HashMap<u64, Vec<usize>> = HashMap::new();
        let mut membership = vec![0usize; locus_ranges.len()];
        for (chunk_base, chunk) in locus_ranges.chunks(2048).enumerate() {
            let built: Vec<(Profile, u64, u32)> = chunk
                .par_iter()
                .map(|traversal| {
                    // The candidate's OWNING universe partition: every
                    // segment of a well-formed traversal shares one owner
                    // (single-segment domain rows and the anchor-owned
                    // split candidates alike; mixed-owner split admission
                    // is excluded at admission).
                    let owner = {
                        let first = traversal.segments[0].partition as u32;
                        ensure(
                            traversal
                                .segments
                                .iter()
                                .all(|segment| segment.partition as u32 == first),
                            "DP allele mixes candidate rows of distinct owning partitions",
                        )?;
                        first
                    };
                    let profile = if traversal.segments.len() == 1 {
                        let segment = &traversal.segments[0];
                        if segment.start == segment.end {
                            Profile::new()
                        } else {
                            parent_profiles
                                .get(&(segment.source, segment.start, segment.end))
                                .cloned()
                                .ok_or_else(|| invalid("allele interval missing parent data"))?
                        }
                    } else {
                        ensure(
                            traversal.segments.len() == 2,
                            "DP allele has unsupported segment arity"
                        )?;
                        let left_segment = &traversal.segments[0];
                        let right_segment = &traversal.segments[1];
                        let mut parts: Vec<&Profile> = Vec::with_capacity(3);
                        if left_segment.start < left_segment.end {
                            let key =
                                (left_segment.source, left_segment.start, left_segment.end);
                            parts.push(
                                partials
                                    .get(&key)
                                    .ok_or_else(|| invalid("split partial missing"))?,
                            );
                        }
                        if right_segment.start < right_segment.end {
                            let key =
                                (right_segment.source, right_segment.start, right_segment.end);
                            parts.push(
                                partials
                                    .get(&key)
                                    .ok_or_else(|| invalid("split partial missing"))?,
                            );
                        }
                        let junction = (
                            segment_flank(
                                sources,
                                &flank_memo,
                                left_segment,
                                false,
                                READ_LENGTH - 1,
                            )?,
                            segment_flank(
                                sources,
                                &flank_memo,
                                right_segment,
                                true,
                                READ_LENGTH - 1,
                            )?,
                        );
                        parts.push(
                            seam_profiles_by_junction
                                .get(&junction)
                                .ok_or_else(|| invalid("split junction seam missing"))?,
                        );
                        merge_profiles(&parts)?
                    };
                    let signature = profile_signature_hash(&profile);
                    Ok((profile, signature, owner))
                })
                .collect::<io::Result<_>>()?;
            for (offset, (profile, signature, owner)) in built.into_iter().enumerate() {
                // The owning partition enters the class key: identical
                // profiles from rows of different owning partitions stay
                // distinct classes, each charged against its own owner's
                // routed share (the red line).
                let class = if let Some(class) = class_buckets
                    .get(&signature)
                    .and_then(|bucket| {
                        bucket
                            .iter()
                            .find(|&&class| {
                                class_profiles[class] == profile && class_owners[class] == owner
                            })
                            .copied()
                    }) {
                    class
                } else {
                    let class = class_profiles.len();
                    class_profiles.push(profile);
                    class_owners.push(owner);
                    class_buckets.entry(signature).or_default().push(class);
                    class
                };
                membership[chunk_base * 2048 + offset] = class;
            }
        }
        let classing_seconds = classing_started.elapsed().as_secs_f64();
        substage_parents += parents_seconds;
        substage_sweeps += sweeps_seconds;
        substage_partials_seams += partials_seconds;
        substage_classing += classing_seconds;
        allele_counts.push(locus_ranges.len());
        locus_classes.push(LocusClasses {
            profiles: class_profiles,
            membership,
            class_owners,
        });
        eprintln!(
            "[dp] locus {locus} classed: alleles {} classes {} splits-in-locus {} \
             | parents {parents_seconds:.2}s sweeps {sweeps_seconds:.2}s \
             partials+seams {partials_seconds:.2}s classing {classing_seconds:.2}s",
            locus_ranges.len(),
            locus_classes[locus].profiles.len(),
            locus_ranges.iter().filter(|t| t.segments.len() == 2).count()
        );
        rss_probe(rss, &format!("dp_classing_locus_{locus}"))?;
    }
    let class_counts: Vec<usize> = locus_classes.iter().map(|l| l.profiles.len()).collect();
    if std::env::var("IMPG_DP_DEBUG").is_ok() {
        for locus in 0..locus_count {
            let native = ranges[locus]
                .iter()
                .position(|traversal| {
                    traversal.segments.len() == 1
                        && traversal.segments[0].source == target_source
                        && !traversal.segments[0].reverse
                        && traversal.segments[0].start < traversal.segments[0].end
                });
            if let Some(native) = native {
                let class = locus_classes[locus].membership[native];
                let profile = &locus_classes[locus].profiles[class];
                let direct = merged_pair_loss(profile, profile, &partition_obs[locus], model)?;
                let obs = &partition_obs[locus];
                let prepared: Vec<(&FeatureKey, u64, f64)> = profile
                    .iter()
                    .map(|(key, &q)| (key, q, obs.get(key).copied().unwrap_or(0.0)))
                    .collect();
                let via_prepared = prepared_pair_loss(&prepared, &prepared, model)?;
                eprintln!(
                    "[dp-dbg] locus {locus} native class {class} class-features {} direct-pair-loss {direct:.2} prepared-pair-loss {via_prepared:.2}",
                    profile.len()
                );
            } else {
                eprintln!("[dp-dbg] locus {locus} has no native allele");
            }
        }
    }
    ensure(
        class_counts.iter().all(|&c| c <= 6_000),
        "locus class count exceeds the pair-table budget"
    )?;
    let profiling_seconds = profiling_started.elapsed().as_secs_f64();

    // ------------------------------------------------ boundary seams (L149 endpoint classes)
    let boundary_started = Instant::now();
    let endpoints_started = Instant::now();
    let endpoints: Vec<Vec<(Vec<u8>, Vec<u8>)>> = ranges
        .iter()
        .map(|locus| {
            locus
                .par_iter()
                .map(|traversal| {
                    allele_endpoints(sources, &flank_memo, &traversal.segments, 149)
                })
                .collect::<io::Result<_>>()
        })
        .collect::<io::Result<_>>()?;
    let boundary_endpoints_seconds = endpoints_started.elapsed().as_secs_f64();
    let links_started = Instant::now();
    let links = genome::port_word_seams(axis_slice, &ranges, graph, &mut ports)?;
    ensure(
        links.len() + 1 == locus_count,
        "boundary link cardinality mismatch"
    )?;
    let boundary_links_seconds = links_started.elapsed().as_secs_f64();
    let compose_started = Instant::now();
    let mut seam_profiles: Vec<Vec<std::sync::Arc<Profile>>> = Vec::with_capacity(links.len());
    let mut successors: Vec<Vec<Vec<(usize, f64)>>> = Vec::with_capacity(links.len());
    let mut seam_by_pair: Vec<HashMap<(usize, usize), std::sync::Arc<Profile>>> =
        Vec::with_capacity(links.len());
    let mut boundary_compositions = 0usize;
    for (boundary, boundary_links) in links.iter().enumerate() {
        let left_ends: Vec<Vec<u8>> = endpoints[boundary]
            .iter()
            .map(|(_, tail)| tail.clone())
            .collect();
        let right_starts: Vec<Vec<u8>> = endpoints[boundary + 1]
            .iter()
            .map(|(head, _)| head.clone())
            .collect();
        let left_owners: Vec<u32> = ranges[boundary]
            .iter()
            .map(|traversal| traversal.segments[0].partition as u32)
            .collect();
        let right_owners: Vec<u32> = ranges[boundary + 1]
            .iter()
            .map(|traversal| traversal.segments[0].partition as u32)
            .collect();
        let (profiles, stats) = genome::deduplicated_boundary_profiles(
            boundary_links,
            &left_ends,
            &right_starts,
            |left_end, right_start| {
                let (profile, _) = genome::profile_event_seam(
                    panel,
                    left_end,
                    right_start,
                    READ_LENGTH,
                    MAX_FEATURES,
                )?;
                Ok(profile)
            },
        )?;
        boundary_compositions += stats.class_compositions;
        let mut per_link_scores: Vec<f64> = Vec::with_capacity(boundary_links.len());
        let mut succ: Vec<Vec<(usize, f64)>> = vec![Vec::new(); ranges[boundary].len()];
        let mut by_pair = HashMap::with_capacity(boundary_links.len());
        for (index, (&(left, right), profile)) in
            boundary_links.iter().zip(profiles.iter()).enumerate()
        {
            // The seam's observed sides are the two adjacent candidates'
            // OWNING partitions' shares (owner-resolved boundary form).
            let obs_left = owner_routed_obs(
                routed_equal,
                left_owners[left],
                component_locus_to_partition,
            );
            let obs_right = owner_routed_obs(
                routed_equal,
                right_owners[right],
                component_locus_to_partition,
            );
            per_link_scores.push(profile_loss_boundary_owned(
                profile,
                obs_left,
                obs_right,
                left_owners[left] == right_owners[right],
                model,
            )?);
            succ[left].push((right, per_link_scores[index]));
            by_pair.insert((left, right), std::sync::Arc::clone(profile));
        }
        seam_profiles.push(profiles);
        successors.push(succ);
        seam_by_pair.push(by_pair);
        rss_probe(rss, &format!("dp_boundary_{boundary}"))?;
    }
    let physical_links: u64 = links.iter().map(Vec::len).sum::<usize>() as u64;
    let boundary_compose_seconds = compose_started.elapsed().as_secs_f64();
    drop(links);
    drop(endpoints);
    let boundary_seconds = boundary_started.elapsed().as_secs_f64();
    rss_probe(rss, "dp_boundary_links")?;
    eprintln!(
        "[dp] boundaries done: links {physical_links} \
         | endpoints {boundary_endpoints_seconds:.2}s links {boundary_links_seconds:.2}s \
         compose+score {boundary_compose_seconds:.2}s"
    );

    // ------------------------------------------------ viability + conflict window
    let conflict_window = conflict_window_of(&ranges).max(1);
    let viable: Vec<Vec<bool>> = {
        let mut forward: Vec<Vec<bool>> = ranges.iter().map(|l| vec![true; l.len()]).collect();
        for locus in (0..locus_count - 1).rev() {
            for (allele, list) in successors[locus].iter().enumerate() {
                forward[locus][allele] = list.iter().any(|&(next, _)| forward[locus + 1][next]);
            }
        }
        let mut backward: Vec<Vec<bool>> = ranges.iter().map(|l| vec![false; l.len()]).collect();
        backward[0].iter_mut().for_each(|value| *value = true);
        for locus in 1..locus_count {
            for (allele, list) in successors[locus - 1].iter().enumerate() {
                if backward[locus - 1][allele] {
                    for &(next, _) in list {
                        backward[locus][next] = true;
                    }
                }
            }
        }
        (0..locus_count)
            .map(|locus| {
                ranges[locus]
                    .iter()
                    .enumerate()
                    .map(|(allele, traversal)| {
                        forward[locus][allele]
                            && backward[locus][allele]
                            && spans_feasible_local(&traversal.segments)
                    })
                    .collect()
            })
            .collect()
    };

    // ------------------------------------------------ beam search (per tie rotation)
    let beam_started = Instant::now();
    // Per-locus folded class machinery (P2): fid-interned prepared classes
    // shared by every tie rotation; pair constants are computed lazily per
    // (locus, class pair) on first query, replacing the every-pair
    // precompute that dominated the old pair_tables stage. IMPG_DP_FULL_TABLES
    // restores the old full precompute (parallel, bit-identical values) for
    // A/B measurement.
    let pair_table_started = Instant::now();
    let folded: Vec<LocusFolded> = (0..locus_count)
        .into_par_iter()
        .map(|locus| {
            LocusFolded::build(
                &locus_classes[locus].profiles,
                &locus_classes[locus].class_owners,
                routed_equal,
                component_locus_to_partition,
            )
        })
        .collect::<io::Result<_>>()?;
    let loss_tables: Vec<Vec<Vec<f64>>> = (0..locus_count)
        .into_par_iter()
        .map(|locus| folded[locus].loss_tables(model).map(|(tables, _)| tables))
        .collect::<io::Result<_>>()?;
    let loss_table_entries: u64 = (0..locus_count)
        .map(|locus| {
            loss_tables[locus]
                .iter()
                .map(|table| table.len() as u64)
                .sum::<u64>()
        })
        .sum::<u64>();
    let mut pair_tables: Vec<LazyPairTable> = if std::env::var("IMPG_DP_FULL_TABLES").is_ok() {
        (0..locus_count)
            .into_par_iter()
            .map(|locus| -> io::Result<LazyPairTable> {
                let class_count = folded[locus].prepared.len();
                let rows: Vec<Vec<(usize, usize, f64)>> = (0..class_count)
                    .into_par_iter()
                    .map(|first| {
                        let mut row = Vec::with_capacity(class_count - first);
                        for second in first..class_count {
                            row.push((
                                first,
                                second,
                                prepared_pair_loss_folded(
                                    &folded[locus].prepared[first],
                                    &folded[locus].prepared[second],
                                    folded[locus].owners[first] == folded[locus].owners[second],
                                    &folded[locus].entries_by_fid,
                                    &loss_tables[locus],
                                    model,
                                )?,
                            ));
                        }
                        Ok(row)
                    })
                    .collect::<io::Result<_>>()?;
                let mut table = vec![0.0f64; class_count * (class_count + 1) / 2];
                for row in rows {
                    for (first, second, loss) in row {
                        table[second * (second + 1) / 2 + first] = loss;
                    }
                }
                LazyPairTable::from_full(
                    class_count,
                    &folded[locus].prepared,
                    &folded[locus].owners,
                    &folded[locus].entries_by_fid,
                    &loss_tables[locus],
                    table,
                    model,
                )
            })
            .collect::<io::Result<_>>()?
    } else {
        folded
            .iter()
            .zip(loss_tables.iter())
            .map(|(tables, losses)| {
                LazyPairTable::new(
                    tables.prepared.len(),
                    &tables.prepared,
                    &tables.owners,
                    &tables.entries_by_fid,
                    losses,
                    model,
                )
            })
            .collect()
    };
    let pair_table_seconds = pair_table_started.elapsed().as_secs_f64();
    rss_probe(rss, "dp_pair_tables")?;
    // Admissible suffix bounds: the least loss any completion of a state at
    // locus L can still add (per-locus minimum pair loss, per-boundary
    // minimum seam score for both copies). Bound pruning then removes every
    // state that cannot beat the incumbent backbone — the A*-bound piece of
    // the oracle machinery, adapted to the additive M1 objective. The
    // per-locus minimum pair loss uses the admissible folded bound
    // (locus_min_pair_lower_bound_folded) instead of the old exact minimum
    // of the fully precomputed table — the pre-registered accounting-only
    // delta (bound values shift down, bound_pruned counters shrink, kept
    // survivors can only grow; the c12/c14 selection comparison decides).
    let min_pair_of: Vec<f64> = folded
        .iter()
        .map(|tables| locus_min_pair_lower_bound_folded(tables, model))
        .collect();
    let suffix_bound: Vec<f64> = {
        let mut bound = vec![0.0f64; locus_count + 1];
        for locus in (0..locus_count).rev() {
            let min_pair = min_pair_of[locus];
            let min_seam = if locus + 1 < locus_count {
                successors[locus]
                    .iter()
                    .flat_map(|list| list.iter().map(|&(_, score)| score))
                    .reduce(f64::min)
                    .unwrap_or(0.0)
            } else {
                0.0
            };
            bound[locus] = bound[locus + 1] + min_pair + 2.0 * min_seam;
        }
        bound
    };
    // Native backbone (seed incumbent): the same-source-contiguous native
    // chain, scored through the class profiles and seam links. Its total is
    // the inherited incumbent every state must beat; its own states are
    // seeded into every layer so the backbone survives as a finalist.
    let backbone: Option<(Vec<usize>, Vec<f64>)> = {
        let first = (0..locus_count)
            .map(|locus| {
                ranges[locus].iter().position(|traversal| {
                    traversal.segments.len() == 1
                        && traversal.segments[0].source == target_source
                        && !traversal.segments[0].reverse
                        && traversal.segments[0].start < traversal.segments[0].end
                })
            })
            .collect::<Vec<Option<usize>>>();
        let mut chain: Vec<usize> = Vec::new();
        let mut cumulative: Vec<f64> = Vec::new();
        let mut total = 0.0f64;
        let mut broken = false;
        for locus in 0..locus_count {
            let native = if locus == 0 {
                first[0]
            } else {
                chain
                    .last()
                    .and_then(|&previous| {
                        successors[locus - 1][previous]
                            .iter()
                            .filter(|&&(next, _)| first[locus] == Some(next))
                            .map(|&(next, score)| (next, score))
                            .next()
                            .map(|(next, _)| next)
                    })
                    .or_else(|| {
                        chain.last().and_then(|&previous| {
                            successors[locus - 1][previous]
                                .iter()
                                .find(|&&(next, _)| {
                                    ranges[locus][next].segments.len() == 1
                                        && ranges[locus][next].segments[0].source == target_source
                                        && !ranges[locus][next].segments[0].reverse
                                        && ranges[locus][next].segments[0].start
                                            < ranges[locus][next].segments[0].end
                                })
                                .map(|&(next, _)| next)
                        })
                    })
            };
            let Some(native) = native else {
                broken = true;
                break;
            };
            let profile = &locus_classes[locus].profiles[locus_classes[locus].membership[native]];
            total += merged_pair_loss(profile, profile, &partition_obs[locus], model)?;
            if locus > 0 {
                let previous = *chain.last().expect("chain");
                let seam = successors[locus - 1][previous]
                    .iter()
                    .find(|&&(next, _)| next == native)
                    .map(|&(_, score)| score)
                    .unwrap_or(f64::NAN);
                total += 2.0 * seam;
            }
            chain.push(native);
            cumulative.push(total);
        }
        (!broken).then_some((chain, cumulative))
    };
    let incumbent = backbone
        .as_ref()
        .map(|(_, cumulative)| *cumulative.last().expect("backbone"))
        .unwrap_or(f64::INFINITY);
    // future_bound[L]: the least loss any completion of a state AT locus L
    // can still add (suffix_bound[L] minus the local term already charged).
    let future_bound: Vec<f64> = (0..locus_count)
        .map(|locus| suffix_bound[locus] - min_pair_of[locus])
        .collect();
    let mut bound_pruned_states = 0u64;
    let mut bound_pruned_edges = 0u64;
    let mut occupancy: Vec<Vec<usize>> = Vec::new();
    let mut dropped_ties = 0u64;
    let mut dropped_non_ties = 0u64;
    let mut candidates_total = 0u64;
    let mut beam_candidate_counts: Vec<u64> = Vec::new();
    let mut beam_initial_seconds = 0.0f64;
    let mut beam_eligible_seconds = 0.0f64;
    let mut beam_expansion_seconds = 0.0f64;
    let mut insert_attempts_total = 0u64;
    let mut rotation_best: Vec<Option<f64>> = Vec::new();
    let mut finalists: Vec<(f64, Vec<[usize; 2]>)> = Vec::new();
    if std::env::var("IMPG_DP_DEBUG").is_ok() {
        if let Some((chain, cumulative)) = &backbone {
            eprintln!(
                "[dp-dbg] native backbone: total {:.2} chain {:?}",
                cumulative.last().copied().unwrap_or(f64::NAN),
                chain
            );
        } else {
            eprintln!("[dp-dbg] native backbone: BROKEN");
        }
    }
    // Admissible per-class weight bounds for the initial layer (locus 0):
    // max(b(c1), b(c2)) lower-bounds every physical pair of the class pair
    // (LocusFolded::class_weight), so the weight-ordered initial enumeration
    // can skip the pair-constant computation for class pairs that provably
    // sit beyond the beam cutoff — kept-set neutral, because a skipped pair's
    // score exceeds the cutoff by more than BOUND_PRUNE_EPSILON (ties
    // included) and the kept set is an exact top-K by (score, key).
    let initial_class_weights: Vec<f64> = (0..locus_classes[0].profiles.len())
        .map(|class| folded[0].class_weight(class, model))
        .collect::<io::Result<_>>()?;
    let initial_viable_by_class: Vec<Vec<usize>> = {
        let membership = &locus_classes[0].membership;
        let mut members = vec![Vec::<usize>::new(); locus_classes[0].profiles.len()];
        for allele in 0..ranges[0].len() {
            if viable[0][allele] {
                members[membership[allele]].push(allele);
            }
        }
        members
    };
    let initial_class_pairs: Vec<(f64, usize, usize)> = {
        let count = initial_class_weights.len();
        let mut pairs = Vec::with_capacity(count * (count + 1) / 2);
        for first in 0..count {
            if initial_viable_by_class[first].is_empty() {
                continue;
            }
            for second in first..count {
                if initial_viable_by_class[second].is_empty() {
                    continue;
                }
                pairs.push((
                    initial_class_weights[first].max(initial_class_weights[second]),
                    first,
                    second,
                ));
            }
        }
        // Cheapest admissible bound first, so the layer fills with the best
        // class pairs and the weight skip fires as early as possible.
        pairs.sort_by(|a, b| {
            a.0.total_cmp(&b.0)
                .then_with(|| (a.1, a.2).cmp(&(b.1, b.2)))
        });
        pairs
    };
    for rotation in 0..tie_rotations.max(1) {
        // Beam pruning is the bounded layer (above): exact kept-set
        // reproduction of the old materialize-then-sort prune, fed folded
        // class-pair edges only (P1).
        // Initial layer: every viable diploid pair at locus 0, enumerated as
        // weight-ordered CLASS PAIRS whose physical members all share one
        // lazily computed pair constant — the grouped shape the transitions
        // use below, applied where the old code pushed all ~3.2M pairs one
        // by one.
        let initial_started = Instant::now();
        let table = &mut pair_tables[0];
        let mut collector = BoundedLayer::new(beam_width);
        let mut layer_candidates = 0u64;
        for &(weight, class_first, class_second) in &initial_class_pairs {
            let first_members = &initial_viable_by_class[class_first];
            let second_members = &initial_viable_by_class[class_second];
            let feasible = if class_first == class_second {
                let count = first_members.len() as u64;
                count * (count + 1) / 2
            } else {
                first_members.len() as u64 * second_members.len() as u64
            };
            // Candidate accounting follows the old push count exactly:
            // every viable physical pair of the class pair was pushed, and
            // the incumbent bound never rejected an initial state in the
            // validated runs (bound_pruned_states 0).
            layer_candidates += feasible;
            if collector.is_full() {
                let worst = collector.worst_entry();
                // Admissible weight skip, tie-safe by the epsilon margin.
                if weight > worst.0 .0 + BOUND_PRUNE_EPSILON {
                    collector.record_bulk_drop(f64::INFINITY, feasible);
                    continue;
                }
            }
            let score = table.entry(class_first, class_second)?;
            if collector.is_full() {
                let worst = collector.worst_entry();
                // Minimal physical key of the class pair: members are
                // ascending, so the smallest pair is the sorted combination
                // of the two smallest members ((a0, a0) within one class).
                let minimum_pair = if class_first == class_second {
                    [first_members[0], first_members[0]]
                } else {
                    [first_members[0].min(second_members[0]), first_members[0].max(second_members[0])]
                };
                if (ScoreKey(score), (minimum_pair, (0, 0))).cmp(&worst)
                    != std::cmp::Ordering::Less
                {
                    // Every member pair sits at or beyond the beam cutoff
                    // exactly as per-pair rejection would decide.
                    collector.record_bulk_drop(score, feasible);
                    continue;
                }
            }
            // Ascending member-pair enumeration with the monotone early
            // break: keys ascend within the class pair, so the first
            // beyond-cutoff pair rejects the whole remaining tail.
            let mut pairs: Vec<[usize; 2]> = Vec::new();
            if class_first == class_second {
                for (offset, &first) in first_members.iter().enumerate() {
                    for &second in &first_members[offset..] {
                        pairs.push([first, second]);
                    }
                }
            } else {
                for &first in first_members {
                    for &second in second_members {
                        pairs.push([first.min(second), first.max(second)]);
                    }
                }
                pairs.sort_unstable();
            }
            let mut enumerated = 0u64;
            'members: for [first, second] in pairs {
                enumerated += 1;
                if collector.is_full() {
                    let worst = collector.worst_entry();
                    if (ScoreKey(score), ([first, second], (0, 0))).cmp(&worst)
                        == std::cmp::Ordering::Greater
                    {
                        let remainder = feasible.saturating_sub(enumerated - 1);
                        collector.record_bulk_drop(score, remainder);
                        break 'members;
                    }
                }
                collector.insert(([first, second], (0, 0)), score, u32::MAX);
            }
        }
        let pruned = collector.finish(&[], conflict_window, future_bound[0], incumbent)?;
        candidates_total += layer_candidates;
        beam_candidate_counts.push(layer_candidates);
        dropped_ties += pruned.dropped_ties;
        dropped_non_ties += pruned.dropped_non_ties;
        bound_pruned_states += pruned.bound_pruned_states;
        insert_attempts_total += pruned.insert_attempts;
        beam_initial_seconds += initial_started.elapsed().as_secs_f64();
        let mut layer = pruned.layer;
        let mut backbone_indices: Vec<Option<u32>> = vec![None; locus_count];
        if let Some((chain, cumulative)) = &backbone {
            let pair = [chain[0], chain[0]];
            let history = vec![pair];
            match layer
                .iter()
                .position(|state| state.pair == pair && state.history.len() == 1)
            {
                Some(index) => backbone_indices[0] = Some(index as u32),
                None => {
                    backbone_indices[0] = Some(layer.len() as u32);
                    layer.push(DpState {
                        score: cumulative[0],
                        pair,
                        pred: None,
                        history: history.clone(),
                        suffix_hash: history_fingerprint(&[]),
                    });
                }
            }
        }
        occupancy.push(vec![layer.len()]);
        let initial_seconds = initial_started.elapsed().as_secs_f64();
        // Transitions.
        let mut transition_seconds = 0.0f64;
        let mut layer_pairs: Vec<Vec<([usize; 2], Option<u32>)>> = Vec::with_capacity(locus_count);
        layer_pairs.push(layer.iter().map(|s| (s.pair, s.pred)).collect());
        for locus in 1..locus_count {
            eprintln!("[dp] rotation {rotation} transition locus {locus}");
            let locus_started = Instant::now();
            let membership = &locus_classes[locus].membership;
            let next_count = ranges[locus].len();
            let mut spans_by_source: BTreeMap<usize, Vec<(u64, u64, usize)>> = BTreeMap::new();
            for (member, traversal) in ranges[locus].iter().enumerate() {
                for segment in &traversal.segments {
                    if segment.start < segment.end {
                        spans_by_source
                            .entry(segment.source)
                            .or_default()
                            .push((segment.start, segment.end, member));
                    }
                }
            }
            for spans in spans_by_source.values_mut() {
                spans.sort();
            }
            // Folded successor groups (P1): per PREVIOUS-locus left allele,
            // the viable next-locus members grouped by (next class, seam
            // score) — the class-graph DP's folded edges. Every member of a
            // group pair contributes the SAME lazily memoized class-pair
            // constant and seam terms, so the state computes ONE constant per
            // group pair and enumerates physical members only in ascending
            // order with the bulk-drop early break; the old unfolded
            // expansion's ~492M per-pair push attempts collapse to
            // O(states × group pairs) constant work.
            let grouped: Vec<Vec<(usize, f64, Vec<usize>)>> = successors[locus - 1]
                .iter()
                .map(|list| {
                    let mut groups: Vec<(usize, f64, Vec<usize>)> = Vec::new();
                    for &(next, seam) in list {
                        if !viable[locus][next] {
                            continue;
                        }
                        let class = membership[next];
                        match groups.iter_mut().find(|group| {
                            group.0 == class && group.1.to_bits() == seam.to_bits()
                        }) {
                            Some(group) => group.2.push(next),
                            None => groups.push((class, seam, vec![next])),
                        }
                    }
                    groups.sort_by(|a, b| {
                        a.0.cmp(&b.0).then_with(|| a.1.to_bits().cmp(&b.1.to_bits()))
                    });
                    for group in &mut groups {
                        group.2.sort_unstable();
                    }
                    groups
                })
                .collect();
            let mut collector = BoundedLayer::new(beam_width);
            let table = &mut pair_tables[locus];
            let mut layer_candidates = 0u64;
            let expansion_started = Instant::now();
            let mut locus_eligible = 0.0f64;
            for (state_index, state) in layer.iter().enumerate() {
                let eligible_started = Instant::now();
                let first_locus = locus - state.history.len();
                let mut eligible = [vec![true; next_count], vec![true; next_count]];
                for (offset, previous) in state.history.iter().enumerate() {
                    let previous_locus = first_locus + offset;
                    for copy in 0..2 {
                        for segment in &ranges[previous_locus][previous[copy]].segments {
                            let Some(spans) = spans_by_source.get(&segment.source) else {
                                continue;
                            };
                            let prefix =
                                spans.partition_point(|&(start, _, _)| start < segment.end);
                            for &(start, end, member) in &spans[..prefix] {
                                if !eligible[copy][member] || end <= segment.start {
                                    continue;
                                }
                                if ranges[locus][member].segments.iter().any(|right| {
                                    segment.source == right.source
                                        && segment.start < right.end
                                        && right.start < segment.end
                                }) {
                                    eligible[copy][member] = false;
                                }
                            }
                        }
                    }
                }
                locus_eligible += eligible_started.elapsed().as_secs_f64();
                rss.checkpoint("dp_transition_hot_loop")?;
                let member_stats = |copy: usize, members: &[usize]| -> (u64, usize) {
                    let mut count = 0u64;
                    let mut first = usize::MAX;
                    for &member in members {
                        if eligible[copy][member] {
                            count += 1;
                            if first == usize::MAX {
                                first = member;
                            }
                        }
                    }
                    (count, first)
                };
                let first_groups = &grouped[state.pair[0]];
                let second_groups = &grouped[state.pair[1]];
                for (class_first, first_seam, first_members) in first_groups {
                    let (first_count, first_start) = member_stats(0, first_members);
                    if first_count == 0 {
                        continue;
                    }
                    for (class_second, second_seam, second_members) in second_groups {
                        let (second_count, second_start) = member_stats(1, second_members);
                        let feasible = first_count.saturating_mul(second_count);
                        if feasible == 0 {
                            continue;
                        }
                        let score = state.score
                            + first_seam
                            + second_seam
                            + table.entry(*class_first, *class_second)?;
                        if score + future_bound[locus] > incumbent + BOUND_PRUNE_EPSILON {
                            // The old expansion counted one bound-pruned
                            // edge per eligible pair; the group carries the
                            // same count in one step.
                            bound_pruned_edges += feasible;
                            continue;
                        }
                        // Candidate accounting: the old expansion pushed
                        // every eligible, viable, bound-passing pair.
                        layer_candidates += feasible;
                        if collector.is_full() {
                            let worst = collector.worst_entry();
                            let minimum_key =
                                ([first_start, second_start], state.suffix_hash);
                            if (ScoreKey(score), minimum_key).cmp(&worst)
                                != std::cmp::Ordering::Less
                            {
                                // Every member pair of this group pair sits
                                // at or beyond the beam cutoff exactly as
                                // per-pair rejection would decide.
                                collector.record_bulk_drop(score, feasible);
                                continue;
                            }
                        }
                        // Ascending member-pair enumeration with the
                        // monotone early break: keys ascend within the group
                        // pair, so the first beyond-cutoff pair rejects the
                        // whole remaining tail in one bulk drop.
                        let mut enumerated = 0u64;
                        'members: for &first in first_members {
                            if !eligible[0][first] {
                                continue;
                            }
                            for &second in second_members {
                                if !eligible[1][second] {
                                    continue;
                                }
                                enumerated += 1;
                                if collector.is_full() {
                                    let worst = collector.worst_entry();
                                    if (ScoreKey(score), ([first, second], state.suffix_hash))
                                        .cmp(&worst)
                                        == std::cmp::Ordering::Greater
                                    {
                                        let remainder =
                                            feasible.saturating_sub(enumerated - 1);
                                        collector.record_bulk_drop(score, remainder);
                                        break 'members;
                                    }
                                }
                                collector.insert(
                                    ([first, second], state.suffix_hash),
                                    score,
                                    state_index as u32,
                                );
                            }
                        }
                    }
                }
            }
            let transition_elapsed = expansion_started.elapsed().as_secs_f64();
            beam_eligible_seconds += locus_eligible;
            beam_expansion_seconds += transition_elapsed - locus_eligible;
            let pruned = collector.finish(&layer, conflict_window, future_bound[locus], incumbent)?;
            candidates_total += layer_candidates;
            beam_candidate_counts.push(layer_candidates);
            dropped_ties += pruned.dropped_ties;
            dropped_non_ties += pruned.dropped_non_ties;
            bound_pruned_states += pruned.bound_pruned_states;
            insert_attempts_total += pruned.insert_attempts;
            let mut next_layer = pruned.layer;
            std::mem::swap(&mut layer, &mut next_layer);
            if let Some((chain, cumulative)) = &backbone {
                let pair = [chain[locus], chain[locus]];
                let history_start = locus.saturating_sub(conflict_window - 1);
                let history: Vec<[usize; 2]> = (history_start..=locus)
                    .map(|l| [chain[l], chain[l]])
                    .collect();
                match layer.iter().position(|state| {
                    state.pair == pair && state.history == history
                }) {
                    Some(index) => backbone_indices[locus] = Some(index as u32),
                    None => {
                        backbone_indices[locus] = Some(layer.len() as u32);
                        layer.push(DpState {
                            score: cumulative[locus],
                            pair,
                            pred: backbone_indices[locus - 1],
                            suffix_hash: history_fingerprint(&history[1..]),
                            history,
                        });
                    }
                }
            }
            layer_pairs.push(layer.iter().map(|s| (s.pair, s.pred)).collect());
            occupancy.last_mut().unwrap().push(layer.len());
            if std::env::var("IMPG_DP_DEBUG").is_ok() {
                let best = layer.first().map(|state| state.score);
                let native_in_layer = backbone
                    .as_ref()
                    .and_then(|(chain, _)| {
                        let native = chain[locus];
                        layer
                            .iter()
                            .position(|state| state.pair[0] == native && state.pair[1] == native)
                    });
                eprintln!(
                    "[dp-dbg] layer {locus}: states {} best {:?} native-pair-state {native_in_layer:?}",
                    layer.len(),
                    best.map(|score| format!("{score:.2}")),
                );
            }
            transition_seconds += locus_started.elapsed().as_secs_f64();
            drop(next_layer);
            rss_probe(rss, &format!("dp_transition_locus_{locus}"))?;
        }
        let _ = initial_seconds;
        rotation_best.push(layer.first().map(|state| state.score));
        // Finalists: best terminal states, backtracked through layer_pairs.
        let mut ranked: Vec<(usize, f64)> = (0..layer.len())
            .map(|index| (index, layer[index].score))
            .collect();
        ranked.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
        for (index, score) in ranked.iter().take(finalist_count) {
            let mut route = Vec::with_capacity(locus_count);
            let mut cursor = (*index, locus_count - 1);
            loop {
                let (state_index, locus) = cursor;
                let (pair, pred) = layer_pairs[locus][state_index];
                route.push(pair);
                match pred {
                    Some(previous) if locus > 0 => cursor = (previous as usize, locus - 1),
                    _ => break,
                }
            }
            route.reverse();
            finalists.push((*score, route));
        }
        let _ = transition_seconds;
    }
    let beam_seconds = beam_started.elapsed().as_secs_f64();
    // Deduplicate finalists across rotations (identical routes).
    finalists.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    finalists.dedup_by(|a, b| a.1 == b.1);

    // ------------------------------------------------ finalist external rescoring
    let rescore_started = Instant::now();
    let mut oracle_memo: HashMap<String, Profile> = HashMap::new();
    let mut m1_scores: Vec<f64> = Vec::with_capacity(finalists.len());
    for (_, route) in &finalists {
        ensure(route.len() == locus_count, "finalist route length mismatch")?;
        let mut total = 0.0;
        for locus in 0..locus_count {
            let first = oracle_allele_profile(
                panel,
                sources,
                &ranges[locus][route[locus][0]],
                &mut oracle_memo,
            )?;
            let second = oracle_allele_profile(
                panel,
                sources,
                &ranges[locus][route[locus][1]],
                &mut oracle_memo,
            )?;
            let owner_first = traversal_owner(&ranges[locus][route[locus][0]]);
            let owner_second = traversal_owner(&ranges[locus][route[locus][1]]);
            if owner_first == owner_second {
                total += merged_pair_loss(
                    &first,
                    &second,
                    owner_routed_obs(routed_equal, owner_first, component_locus_to_partition),
                    model,
                )?;
            } else {
                // Fix 1 (record-once mixed-owner charging): the merged pair
                // charge's observed side attributes each record ONCE — the
                // summed two-map form counted a record touching both owners
                // twice.
                let obs_site = site.site_map([
                    owner_universe_partition(owner_first, component_locus_to_partition),
                    owner_universe_partition(owner_second, component_locus_to_partition),
                ]);
                total += merged_pair_loss(&first, &second, &obs_site, model)?;
            }
        }
        for boundary in 0..locus_count - 1 {
            for copy in 0..2 {
                let key = (route[boundary][copy], route[boundary + 1][copy]);
                let profile = seam_by_pair[boundary]
                    .get(&key)
                    .ok_or_else(|| invalid("finalist boundary pair lacks a seam profile"))?;
                let owner_left = traversal_owner(&ranges[boundary][route[boundary][copy]]);
                let owner_right = traversal_owner(&ranges[boundary + 1][route[boundary + 1][copy]]);
                // The two-side junction seam keeps the boundary convention's
                // summed form (a junction-straddling record's anchors
                // distribute between exactly the two adjacent partitions —
                // the established, pre-extension bit-reconciled junction
                // semantics; the record-once form is for the per-window
                // LOCAL charges and the pooled composition maps).
                total += profile_loss_boundary_owned(
                    profile,
                    owner_routed_obs(routed_equal, owner_left, component_locus_to_partition),
                    owner_routed_obs(routed_equal, owner_right, component_locus_to_partition),
                    owner_left == owner_right,
                    model,
                )?;
            }
        }
        m1_scores.push(total);
        rss_probe(rss, &format!("dp_rescore_{}", m1_scores.len()))?;
    }
    let mut selected_index = 0usize;
    for (index, score) in m1_scores.iter().enumerate() {
        if score.total_cmp(&m1_scores[selected_index]).is_lt() {
            selected_index = index;
        }
    }
    // Unchanged pooled external rescore of the selected route (continuity
    // with the old references; run outputs only).
    let pooled_external: f64 = {
        let route = &finalists[selected_index].1;
        let mut sequences = Vec::with_capacity(2);
        for copy in 0..2 {
            let mut sequence = Vec::new();
            for locus in 0..locus_count {
                sequence.extend(allele_sequence(
                    sources,
                    &ranges[locus][route[locus][copy]].segments,
                )?);
            }
            sequences.push(sequence);
        }
        let directory = rescore_directory.join("selected-pooled");
        let (loss, cost, runs) = genome::external_rescore(
            panel,
            [&sequences[0], &sequences[1]],
            sample_counts,
            model,
            &directory,
        )?;
        let _ = (cost, runs);
        loss
    };
    let rescore_seconds = rescore_started.elapsed().as_secs_f64();

    // ------------------------------------------------ reference scoring (new model)
    let references_started = Instant::now();
    let mut reference_scores = Vec::new();
    for (index, (label, path_a, path_b)) in reference_specs.iter().enumerate() {
        let route_a: routes::Route = read_json(path_a)?;
        let route_b: routes::Route = read_json(path_b)?;
        let (loss, _) = score_reference_m1(
            panel,
            sources,
            territory,
            [&route_a, &route_b],
            routed_equal,
            site,
            window_obs,
            // Quarantined beam path: no instance structure (the
            // feature-level chain-level form, unchanged).
            None,
            component_locus_to_partition,
            locus_count,
            locus_offset,
            model,
            // Quarantined beam path: pooled junction scoring unchanged.
            None,
            // Quarantined beam path: flat backgrounds unchanged.
            None,
            // Quarantined beam path: panel convention everywhere.
            None,
            // Quarantined beam path: legacy group-shared behavior.
            false,
            None,
            None,
            None,
        )?;
        // Continuity number: the UNCHANGED pooled external rescore of the
        // same reference pair, assembling each copy from the route's OWN
        // segments in file order (the assessment convention; no per-locus
        // decomposition).
        let mut sequences = Vec::with_capacity(2);
        for route in [&route_a, &route_b] {
            let mut sequence = Vec::new();
            for segment in &route.segments {
                if segment.start >= segment.end {
                    continue;
                }
                let mut part = sources.fetch(segment.source, segment.start, segment.end)?;
                if segment.reverse {
                    part = impg::graph::reverse_complement(&part);
                }
                sequence.extend(part);
            }
            sequences.push(sequence);
        }
        let directory = rescore_directory.join(format!("reference-{index}-{label}"));
        let (pooled, pooled_cost, pooled_runs) = genome::external_rescore(
            panel,
            [&sequences[0], &sequences[1]],
            sample_counts,
            model,
            &directory,
        )?;
        reference_scores.push(serde_json::json!({
            "label": label,
            "m1_oracle_rescore": loss,
            "pooled_external_rescore": pooled,
            "pooled_mem_queries": pooled_cost.mem_queries,
            "pooled_initial_runs": pooled_runs,
        }));
        rss_probe(rss, &format!("dp_reference_{label}"))?;
    }
    let references_seconds = references_started.elapsed().as_secs_f64();

    let finalist_rows: Vec<serde_json::Value> = finalists
        .iter()
        .zip(m1_scores.iter())
        .enumerate()
        .map(|(index, ((proposal, route), m1))| {
            let route_json: Vec<[serde_json::Value; 2]> = (0..locus_count)
                .map(|locus| {
                    [
                        traversal_json(&ranges[locus][route[locus][0]]),
                        traversal_json(&ranges[locus][route[locus][1]]),
                    ]
                })
                .collect();
            serde_json::json!({
                "index": index,
                "proposal_loss": proposal,
                "m1_oracle_rescore": m1,
                "selected": index == selected_index,
                "route": route_json,
            })
        })
        .collect();
    Ok(serde_json::json!({
        "model": "routed-evidence-dp-m1-v1",
        "loci": locus_count,
        "locus_offset": locus_offset,
        "beam_width": beam_width,
        "tie_rotations": tie_rotations.max(1),
        "conflict_window": conflict_window,
        "allele_counts": allele_counts,
        "class_counts": class_counts,
        "split_candidate_counts": split_candidate_counts,
        "physical_links": physical_links,
        "boundary_compositions": boundary_compositions,
        "geometric_ir_extractions": ir_extractions,
        "split_seam_queries": split_seam_queries,
        "classing_substage_seconds": {
            "parents": substage_parents,
            "sweeps": substage_sweeps,
            "partials_seams": substage_partials_seams,
            "classing": substage_classing,
        },
        "boundary_substage_seconds": {
            "endpoints": boundary_endpoints_seconds,
            "links": boundary_links_seconds,
            "compose_score": boundary_compose_seconds,
        },
        "pair_entries_computed": pair_tables.iter().map(|t| t.computed_pairs).sum::<u64>(),
        "pair_entries_computed_per_locus": pair_tables
            .iter()
            .map(|t| t.computed_pairs)
            .collect::<Vec<_>>(),
        "beam_substage_seconds": {
            "initial": beam_initial_seconds,
            "transition_eligible": beam_eligible_seconds,
            "transition_expansion": beam_expansion_seconds,
        },
        "insert_attempts": insert_attempts_total,
        "loss_table_entries": loss_table_entries,
        "dropped_ties": dropped_ties,
        "dropped_non_ties": dropped_non_ties,
        "candidates_total": candidates_total,
        "beam_candidate_counts": beam_candidate_counts,
        "bound_pruned_states": bound_pruned_states,
        "bound_pruned_edges": bound_pruned_edges,
        "backbone_incumbent": incumbent,
        "backbone_present": backbone.is_some(),
        "occupancy": occupancy,
        "rotation_best_scores": rotation_best,
        "finalists": finalist_rows,
        "selected_index": selected_index,
        "selected_m1_oracle_rescore": m1_scores[selected_index],
        "selected_pooled_external_rescore": pooled_external,
        "rss_peak_bytes": rss.peak_bytes(),
        "references": reference_scores,
        "stage_wall_seconds": {
            "admission": admission_seconds,
            "classing": profiling_seconds,
            "boundaries": boundary_seconds,
            "beam": beam_seconds,
            "pair_tables": pair_table_seconds,
            "rescore": rescore_seconds,
            "references": references_seconds,
            "total": started.elapsed().as_secs_f64(),
        },
    }))
}

fn traversal_json(traversal: &genome::SpanningTraversal) -> serde_json::Value {
    serde_json::json!({
        "identity": traversal.identity,
        "segments": traversal.segments.iter().map(|segment| serde_json::json!({
            "source": segment.source,
            "start": segment.start,
            "end": segment.end,
            "reverse": segment.reverse,
        })).collect::<Vec<_>>(),
    })
}

/// The haploid chain's EXPLAINED-feature set (supervisor ruling,
/// genome/stitching-omission-alignment 2026-09-24): the union of the CHARGED
/// per-window run profiles' features across the whole route — built with the
/// SAME construction the charging pass uses (the piece merge rule, the
/// group-shared-once filter in route order with its OWN filter state, the
/// interior seams merged exactly where the charging pass merges them:
/// co-occurring seams always, novel seams only when no junction context is
/// available; novel seams under junction context are EXCLUDED — their
/// restricted charge is the row's own payment, so their features are not
/// spelled by the charged profile). The omission addend of every window
/// covers the window's observed features the chain explains NOWHERE
/// (Ruling 2's once-per-material alignment: material the chain spells at its
/// first window is no longer re-charged as omission at its later windows).
/// `piece_memo` is shared with the charging pass (the pre-pass warms it).
fn haploid_chain_explained_features(
    panel: &SyngIndex,
    sources: &routes::Sources,
    local_pieces: &[[Vec<(usize, u64, u64, bool, u32)>; 2]],
    piece_memo: &mut HashMap<(usize, u64, u64), Profile>,
    junction_context: bool,
    // The instance-level correction's switch: when false, no positional
    // work happens (the established feature-level form, bit-identically).
    collect_spans: bool,
    path_of_source: &[usize],
    // The spelled pieces' piece-relative anchor spans, memoized per piece
    // WITH its orientation (the span-to-path mapping is orientation
    // dependent, unlike the canonical profiles' memo key).
    spelled_memo: &mut SpelledMemo,
    // OUT: the chain's spelled anchor positions per feature (the
    // instance-level exoneration's spelled side). Union over the FILTERED
    // pieces (the shared-once filter's own convention: a piece charged at
    // its first window spells at its own coordinates; a later window of the
    // same group observes the SAME coordinates, so one spelling covers it —
    // the group-shared one-instance case) and over the merged seams.
    spans_out: &mut HashMap<FeatureKey, Vec<SpelledSpan>>,
) -> io::Result<std::collections::HashSet<FeatureKey>> {
    let mut explained = std::collections::HashSet::new();
    let mut seen_pieces: std::collections::HashSet<(usize, u64, u64, bool)> =
        std::collections::HashSet::new();
    for locus_pieces in local_pieces {
        // The merge rule, exactly the charging pass's (same source and
        // orientation with no gap and the same owning partition).
        let mut merged_pieces: Vec<(usize, u64, u64, bool, u32)> = Vec::new();
        for piece in &locus_pieces[0] {
            match merged_pieces.last_mut() {
                Some(last)
                    if last.0 == piece.0
                        && last.3 == piece.3
                        && last.4 == piece.4
                        && last.2 >= piece.1
                        && last.2 <= piece.2 =>
                {
                    last.2 = last.2.max(piece.2);
                }
                _ => merged_pieces.push(*piece),
            }
        }
        // The shared-once filter in route order (own state — the charging
        // pass keeps its own; the first-window charge decision is identical).
        let pieces: Vec<(usize, u64, u64, bool, u32)> = merged_pieces
            .iter()
            .filter(|piece| seen_pieces.insert((piece.0, piece.1, piece.2, piece.3)))
            .copied()
            .collect();
        let mut run = Profile::new();
        let mut previous: Option<(Vec<u8>, (usize, u64, u64, bool, u32))> = None;
        for &(source, start, end, reverse, owner) in &pieces {
            let mut sequence = sources.fetch(source, start, end)?;
            if reverse {
                sequence = impg::graph::reverse_complement(&sequence);
            }
            let piece_spans: std::rc::Rc<BTreeMap<FeatureKey, Vec<(u64, u64)>>> = if collect_spans {
                match spelled_memo.get(&(source, start, end, reverse)) {
                    Some(spans) => spans.clone(),
                    None => {
                        let spans = std::rc::Rc::new(oracle_segment_spans(panel, &sequence)?);
                        spelled_memo.insert((source, start, end, reverse), spans.clone());
                        spans
                    }
                }
            } else {
                std::rc::Rc::new(BTreeMap::new())
            };
            if collect_spans {
                spelled_span_derivation_diagnosis(
                    &piece_spans,
                    (source, start, end, reverse),
                    &sequence,
                    panel,
                    path_of_source,
                )?;
                piece_spans_to_path(path_of_source, (source, start, end, reverse), &piece_spans, spans_out)?;
            }
            let piece_profile = match piece_memo.get(&(source, start, end)) {
                Some(profile) => profile.clone(),
                None => {
                    let profile = oracle_segment_profile(panel, &sequence)?;
                    piece_memo.insert((source, start, end), profile.clone());
                    profile
                }
            };
            if let Some((previous_sequence, previous_piece)) = previous {
                let (seam, _) = genome::profile_event_seam(
                    panel,
                    &previous_sequence,
                    &sequence,
                    READ_LENGTH,
                    MAX_FEATURES,
                )?;
                let cooccurring = spine::junction::segment_pair_gap(
                    &spine::junction::junction_range(
                        previous_piece.0,
                        previous_piece.1,
                        previous_piece.2,
                        previous_piece.3,
                        u32::MAX,
                    ),
                    &spine::junction::junction_range(source, start, end, reverse, u32::MAX),
                )
                .is_some_and(|gap| gap >= 0);
                let merge_seam = cooccurring || !junction_context;
                if merge_seam && collect_spans {
                    // The merged seam's features are spelled at the junction:
                    // their anchor positions split between the two pieces'
                    // sources (a junction-straddling subwalk's anchors
                    // distribute on both sides; each side's anchor run is
                    // one span on its own path).
                    let seam_spans =
                        profile_event_seam_spans(panel, &previous_sequence, &sequence, READ_LENGTH)?;
                    let left_len = previous_sequence.len() as u64;
                    let left_take = left_len.min(READ_LENGTH as u64 - 1);
                    for (feature, list) in &seam_spans {
                        let entry = spans_out.entry(feature.clone()).or_default();
                        for &(is_left, lo, hi) in list {
                            let (piece, off_lo, off_hi) = if is_left {
                                let base = left_len - left_take;
                                (
                                    (
                                        previous_piece.0,
                                        previous_piece.1,
                                        previous_piece.2,
                                        previous_piece.3,
                                    ),
                                    base + lo,
                                    base + hi,
                                )
                            } else {
                                ((source, start, end, reverse), lo - left_take, hi - left_take)
                            };
                            let (plo, phi) = piece_span_coords(piece, off_lo, off_hi);
                            let path = *path_of_source
                                .get(piece.0)
                                .ok_or_else(|| invalid("piece source absent from the path map"))?;
                            entry.push(SpelledSpan { path, lo: plo, hi: phi });
                        }
                    }
                }
                if previous_piece.4 == owner {
                    run = merge_profiles(&[&run, &piece_profile])?;
                    if merge_seam {
                        run = merge_profiles(&[&run, &seam])?;
                    }
                } else {
                    // Owner change, the haploid merged shape: a co-occurring
                    // seam merges into the run across the owner change.
                    if merge_seam {
                        run = merge_profiles(&[&run, &seam])?;
                    }
                    run = merge_profiles(&[&run, &piece_profile])?;
                }
            } else {
                run = merge_profiles(&[&run, &piece_profile])?;
            }
            previous = Some((sequence, (source, start, end, reverse, owner)));
        }
        for key in run.keys() {
            explained.insert(key.clone());
        }
    }
    Ok(explained)
}

/// Assessment-side reference scorer under the new model: decomposes each
/// reference route per locus via the component's public BED territory
/// intervals (geographic intersection on matching sources, route order
/// preserved), profiles each copy's local pieces with the ORACLE scorer
/// (interiors plus the surviving L149 seam machinery between consecutive
/// pieces and across locus boundaries), and charges per-locus merged pair
/// profiles against the partition's routed shares plus per-boundary seam
/// profiles against the adjacent-share sum — the M1 loss with oracle q.
fn score_reference_m1(
    panel: &SyngIndex,
    sources: &routes::Sources,
    territory: &[Vec<SourceRange>],
    route_pair: [&routes::Route; 2],
    routed_equal: &RoutedObs,
    // The record-once site-map builder (Fix 1's mixed-owner charging) and
    // the windows' full observed profiles (Fix 2's omission universes,
    // indexed by the slice-local locus).
    site: &SiteObserved,
    window_obs: &[HashMap<FeatureKey, f64>],
    // The windows' observed INSTANCE structure (the S2 record-level
    // coverage's observed side). `None` keeps the established feature-level
    // chain-level form bit-identically (the quarantined beam path and the
    // diploid branch).
    window_instances: Option<&InstanceStructure>,
    component_locus_to_partition: &[u32],
    // The scored slice's locus count (the owner-resolved pieces and charges
    // cover exactly the slice's component loci).
    locus_count: usize,
    locus_offset: usize,
    model: &ScoreModel,
    // The within-read adjacency index: when present, NOVEL junctions (the
    // reference's real recombination points) are paid by the reads that
    // actually cross them; co-occurring real panel junctions keep the
    // pooled charge bit-identically. `None` keeps the quarantined beam
    // path's pooled scoring unchanged.
    junction: Option<(&spine::junction::JunctionSpanIndex, &[usize], &[u32])>,
    // Per-feature Poisson backgrounds (the multiplicity background model).
    // `None` keeps the quarantined beam path's flat-background scoring
    // bit-identically; `Some` charges every interior feature against its
    // own panel-multiplicity background, extending the scan first for any
    // charged feature not yet measured (one bounded recursion: the second
    // pass finds every feature present and does not recurse further).
    mut backgrounds: Option<&mut FeatureBackgrounds>,
    // The sample-side cross-support backgrounds (haploid single-allele
    // convention): a reference whose SECOND route is empty is a HAPLOID
    // statement, and its copy-0 charges use the sample-side backgrounds
    // (the same convention as the phasing haploid track — Step B,
    // owner-approved 2026-09-23). `None` keeps the panel convention
    // everywhere (the quarantined beam path).
    sample: Option<&SampleSideBackgrounds>,
    // GROUP-SHARED MATERIAL CHARGED ONCE (supervisor decision 2026-09-23,
    // assessment-harness bookkeeping — not model semantics): a partition
    // group can serve several axis windows, so the same route material
    // can be extracted as a piece at more than one locus (measured:
    // 9564:[80018,90183) at loci 8 and 17, both partition110). Under this
    // convention each distinct piece is charged at its FIRST locus window
    // only; later windows of the same group spell nothing for that copy
    // (the measured double charge was 2 x +8,621.4 on the haploid truth).
    // This changes reported margins only — never the selection. The legacy
    // behavior (charge at every window) is selected with `false` and
    // reported alongside during the transition.
    charge_group_shared_once: bool,
    // SHARED PIECE-PROFILE CACHE across repeated evaluations of the same
    // slice (the finalist re-ranking evaluates many candidate chains; the
    // chains' territory pieces recur). The memo maps (source, start, end)
    // to the SAME oracle profile either way, so values are bit-identical
    // with and without sharing; only the profile re-computation is saved.
    // `None` keeps a call-local memo (the ladder references' behavior).
    mut shared_piece_memo: Option<&mut HashMap<(usize, u64, u64), Profile>>,
    // SHARED SPELLED-SPAN cache across repeated evaluations (the same
    // sharing argument as the piece-profile memo: the chains' territory
    // pieces recur; the spelled spans are positions only, values
    // bit-identical either way). `None` keeps a call-local memo.
    mut shared_spelled_memo: Option<&mut SpelledMemo>,
    // OUT (the instance-probe diagnostic): the evaluated chain's spelled
    // anchor positions per feature. `None` (the default) shares nothing.
    mut spelled_out: Option<&mut HashMap<FeatureKey, Vec<SpelledSpan>>>,
) -> io::Result<(f64, Vec<serde_json::Value>)> {
    // A reference whose second route is empty is a HAPLOID statement (the
    // ploidy-track convention): its copy-0 charges use the sample-side
    // backgrounds when provided.
    let haploid_reference = sample.is_some() && route_pair[1].segments.is_empty();
    // Distinct pieces already charged, per copy (the group-shared-once
    // convention's bookkeeping).
    let mut charged_pieces: [HashSet<(usize, u64, u64, bool)>; 2] =
        [HashSet::new(), HashSet::new()];
    let mut local_piece_memo: HashMap<(usize, u64, u64), Profile> = HashMap::new();
    let piece_memo: &mut HashMap<(usize, u64, u64), Profile> = match shared_piece_memo {
        Some(map) => map,
        None => &mut local_piece_memo,
    };
    let mut local_spelled_memo: SpelledMemo = HashMap::new();
    let spelled_memo: &mut SpelledMemo = match shared_spelled_memo {
        Some(map) => map,
        None => &mut local_spelled_memo,
    };
    let mut total = 0.0;
    let mut census: Vec<serde_json::Value> = Vec::new();
    // Features charged by this reference that the background scan has not
    // measured yet (collected only under `Some(backgrounds)`); when any
    // exist the scan is extended and the scoring recurses ONCE on the
    // complete backgrounds (the second pass finds every feature present
    // and does not collect or recurse further).
    let mut missing: Vec<FeatureKey> = Vec::new();
    // Per (locus, copy): local pieces in route order, each with the OWNING
    // universe partition of the territory row it intersects (the
    // window-domain extension's owner-resolved charging; every territory
    // row of the pre-extension model is owned by the locus's axis
    // partition, so the pieces are bit-identical there).
    let mut local_pieces: Vec<[Vec<(usize, u64, u64, bool, u32)>; 2]> =
        vec![[Vec::new(), Vec::new()]; locus_count];
    for (copy, route) in route_pair.iter().enumerate() {
        for segment in &route.segments {
            for locus in 0..locus_count {
                for interval in &territory[locus_offset + locus] {
                    if interval.source != segment.source {
                        continue;
                    }
                    let lo = segment.start.max(interval.start);
                    let hi = segment.end.min(interval.end);
                    if lo < hi {
                        local_pieces[locus][copy].push((
                            segment.source,
                            lo,
                            hi,
                            segment.reverse,
                            interval.partition as u32,
                        ));
                    }
                }
            }
        }
    }
    let mut previous_tails: [Vec<u8>; 2] = [Vec::new(), Vec::new()];
    let mut previous_pieces: [Option<(usize, u64, u64, bool, u32)>; 2] = [None, None];
    let mut chain_spelled: HashMap<FeatureKey, Vec<SpelledSpan>> = HashMap::new();
    // Ruling 2 (the chain-level once-per-material omission): the chain's
    // explained-feature set, built BEFORE the charging loop with the same
    // construction (a window's omission addend covers the window's observed
    // features the chain explains NOWHERE along the whole route). The
    // diploid branch's post-hoc addend keeps its per-window form (the
    // frozen-track rule) and uses no explained set.
    let chain_explained = if haploid_reference {
        // The instance-level correction activates only with the windows'
        // instance structure present; without it the established
        // feature-level form runs bit-identically (no positional work).
        let collect_spans = window_instances.is_some();
        let source_paths = sources_path_map(panel, sources)?;
        let mut spans_out: HashMap<FeatureKey, Vec<SpelledSpan>> = HashMap::new();
        let explained = haploid_chain_explained_features(
            panel,
            sources,
            &local_pieces,
            &mut *piece_memo,
            junction.is_some(),
            collect_spans,
            &source_paths,
            spelled_memo,
            &mut spans_out,
        )?;
        chain_spelled = spans_out;
        if let Some(out) = spelled_out.as_deref_mut() {
            *out = chain_spelled.clone();
        }
        explained
    } else {
        std::collections::HashSet::new()
    };
    for locus in 0..locus_count {
        // Per copy: the piece stream decomposed into OWNER RUNS — maximal
        // stretches of same-owner pieces whose interior co-occurring seams
        // merge into the run's profile exactly as the pre-extension copy
        // profile did. Seams at owner changes and every novel junction are
        // charged separately (the boundary convention across the distinct
        // owners). With a single owner per locus (the pre-extension shape)
        // there is exactly one run per copy and every charge is
        // bit-identical to the unextended model.
        let mut copy_runs: [Vec<(u32, Profile)>; 2] = [Vec::new(), Vec::new()];
        let mut tails: [Vec<u8>; 2] = [Vec::new(), Vec::new()];
        let mut copy_novel_charges = [0.0f64; 2];
        for copy in 0..2 {
            // Merge route-contiguous pieces: same source and orientation with
            // no gap AND the same owning partition — a contiguous route
            // segment crossing several territory rows of one owner is ONE
            // piece (no internal seam; seams only at real route junctions or
            // owner changes).
            let raw_pieces: Vec<(usize, u64, u64, bool, u32)> = local_pieces[locus][copy].clone();
            let mut merged_pieces: Vec<(usize, u64, u64, bool, u32)> = Vec::new();
            for piece in raw_pieces {
                match merged_pieces.last_mut() {
                    Some(last)
                        if last.0 == piece.0
                            && last.3 == piece.3
                            && last.4 == piece.4
                            && last.2 >= piece.1
                            && last.2 <= piece.2 =>
                    {
                        last.2 = last.2.max(piece.2);
                    }
                    _ => merged_pieces.push(piece),
                }
            }
            // Group-shared material charged ONCE: a piece extracted at an
            // earlier window of the same partition group is not spelled
            // again here (the supervisor-approved harness convention; the
            // piece is charged at its first window only).
            let pieces: Vec<(usize, u64, u64, bool, u32)> = if charge_group_shared_once {
                merged_pieces
                    .iter()
                    .filter(|piece| {
                        charged_pieces[copy]
                            .insert((piece.0, piece.1, piece.2, piece.3))
                    })
                    .copied()
                    .collect()
            } else {
                merged_pieces
            };
            let pieces = &pieces;
            let mut runs: Vec<(u32, Profile)> = Vec::new();
            let mut previous_sequence: Option<Vec<u8>> = None;
            let mut previous_piece: Option<(usize, u64, u64, bool, u32)> = None;
            let mut novel_charge = 0.0f64;
            let mut head: Option<Vec<u8>> = None;
            let mut head_piece: Option<(usize, u64, u64, bool, u32)> = None;
            let mut tail: Vec<u8> = Vec::new();
            for &(source, start, end, reverse, owner) in pieces {
                let mut sequence = sources.fetch(source, start, end)?;
                if reverse {
                    sequence = impg::graph::reverse_complement(&sequence);
                }
                let piece_profile = match piece_memo.get(&(source, start, end)) {
                    Some(profile) => profile.clone(),
                    None => {
                        let profile = oracle_segment_profile(panel, &sequence)?;
                        piece_memo.insert((source, start, end), profile.clone());
                        profile
                    }
                };
                if let Some(previous) = previous_sequence {
                    let (seam, _) = genome::profile_event_seam(
                        panel,
                        &previous,
                        &sequence,
                        READ_LENGTH,
                        MAX_FEATURES,
                    )?;
                    let prev = previous_piece.expect("piece tracked with sequence");
                    let piece = (source, start, end, reverse, owner);
                    let cooccurring = spine::junction::segment_pair_gap(
                        &spine::junction::junction_range(prev.0, prev.1, prev.2, prev.3, u32::MAX),
                        &spine::junction::junction_range(source, start, end, reverse, u32::MAX),
                    )
                    .is_some_and(|gap| gap >= 0);
                    // Same-owner interior junction: exactly the
                    // pre-extension handling (pooled seams merge into the
                    // run profile; novel junctions are restricted-charged
                    // against the owner).
                    if prev.4 == owner {
                        let run = runs
                            .last_mut()
                            .expect("same-owner piece continues a run");
                        match (cooccurring, junction) {
                            (true, _) => {
                                run.1 = merge_profiles(&[&run.1, &piece_profile])?;
                                run.1 = merge_profiles(&[&run.1, &seam])?;
                            }
                            (false, Some((span, path_of_source, _partitions))) => {
                                run.1 = merge_profiles(&[&run.1, &piece_profile])?;
                                let owner_partition =
                                    owner_universe_partition(owner, component_locus_to_partition);
                                let outcome = span.restricted_charge(
                                    &seam,
                                    &[(
                                        spine::junction::junction_range(
                                            prev.0, prev.1, prev.2, prev.3, owner_partition,
                                        ),
                                        spine::junction::junction_range(
                                            source, start, end, reverse, owner_partition,
                                        ),
                                    )],
                                    sources,
                                    path_of_source,
                                    &[owner_partition],
                                    model,
                                )?;
                                novel_charge += outcome.charge;
                                census.push(serde_json::json!({
                                    "copy": copy,
                                    "locus": locus,
                                    "kind": "interior",
                                    "left": [prev.0, prev.1, prev.2, prev.3],
                                    "right": [source, start, end, reverse],
                                    "spanning_reads": outcome.spanning_reads,
                                    "events": outcome.events,
                                    "span_features": outcome.span_feature_count,
                                    "restricted_charge": outcome.charge,
                                }));
                            }
                            (false, None) => {
                                run.1 = merge_profiles(&[&run.1, &piece_profile])?;
                                run.1 = merge_profiles(&[&run.1, &seam])?;
                            }
                        }
                    } else if haploid_reference {
                        // Fix 1 + Fix 2 restructure (owner-approved
                        // 2026-09-24): the haploid local charge is ONE site
                        // per window — the merged copy profile against the
                        // window's record-once observed map — so a
                        // co-occurring seam merges into the merged profile
                        // regardless of the owner change (the decompose's
                        // shape; the per-owner boundary charge here would
                        // pool the two owners' maps — the union-sum form the
                        // red line forbids). Novel junctions stay restricted
                        // (background-invariant, unchanged).
                        if let Some(bg) = &backgrounds {
                            for key in seam.keys() {
                                if !bg.contains(key) {
                                    missing.push(key.clone());
                                }
                            }
                        }
                        match (cooccurring, junction) {
                            (false, Some((span, path_of_source, _partitions))) => {
                                let prev_partition =
                                    owner_universe_partition(prev.4, component_locus_to_partition);
                                let this_partition =
                                    owner_universe_partition(owner, component_locus_to_partition);
                                let outcome = span.restricted_charge(
                                    &seam,
                                    &[(
                                        spine::junction::junction_range(
                                            prev.0, prev.1, prev.2, prev.3, prev_partition,
                                        ),
                                        spine::junction::junction_range(
                                            source, start, end, reverse, this_partition,
                                        ),
                                    )],
                                    sources,
                                    path_of_source,
                                    &[prev_partition, this_partition],
                                    model,
                                )?;
                                novel_charge += outcome.charge;
                                census.push(serde_json::json!({
                                    "copy": copy,
                                    "locus": locus,
                                    "kind": "interior",
                                    "left": [prev.0, prev.1, prev.2, prev.3],
                                    "right": [source, start, end, reverse],
                                    "spanning_reads": outcome.spanning_reads,
                                    "events": outcome.events,
                                    "span_features": outcome.span_feature_count,
                                    "restricted_charge": outcome.charge,
                                }));
                            }
                            _ => {
                                if let Some(run) = runs.last_mut() {
                                    run.1 = merge_profiles(&[&run.1, &seam])?;
                                }
                            }
                        }
                        runs.push((owner, piece_profile));
                    } else {
                        // Owner change inside the locus: the seam between
                        // the two owners' material charges as a boundary
                        // (pooled forms) or a restricted novel charge
                        // across both owners; the piece starts a new run.
                        if let Some(bg) = &backgrounds {
                            for key in seam.keys() {
                                if !bg.contains(key) {
                                    missing.push(key.clone());
                                }
                            }
                        }
                        let obs_prev = owner_routed_obs(
                            routed_equal,
                            prev.4,
                            component_locus_to_partition,
                        );
                        let obs_this =
                            owner_routed_obs(routed_equal, owner, component_locus_to_partition);
                        let owners_equal = prev.4 == owner;
                        match (cooccurring, junction) {
                            (false, Some((span, path_of_source, _partitions))) => {
                                let prev_partition =
                                    owner_universe_partition(prev.4, component_locus_to_partition);
                                let this_partition =
                                    owner_universe_partition(owner, component_locus_to_partition);
                                let outcome = span.restricted_charge(
                                    &seam,
                                    &[(
                                        spine::junction::junction_range(
                                            prev.0, prev.1, prev.2, prev.3, prev_partition,
                                        ),
                                        spine::junction::junction_range(
                                            source, start, end, reverse, this_partition,
                                        ),
                                    )],
                                    sources,
                                    path_of_source,
                                    &[prev_partition, this_partition],
                                    model,
                                )?;
                                novel_charge += outcome.charge;
                                census.push(serde_json::json!({
                                    "copy": copy,
                                    "locus": locus,
                                    "kind": "interior",
                                    "left": [prev.0, prev.1, prev.2, prev.3],
                                    "right": [source, start, end, reverse],
                                    "spanning_reads": outcome.spanning_reads,
                                    "events": outcome.events,
                                    "span_features": outcome.span_feature_count,
                                    "restricted_charge": outcome.charge,
                                }));
                            }
                            _ => {
                                total += match (&backgrounds, haploid_reference) {
                                    (Some(_), true) => profile_loss_boundary_sample(
                                        &seam,
                                        obs_prev,
                                        obs_this,
                                        sample.expect("haploid reference implies sample"),
                                        model,
                                    )?,
                                    (Some(bg), false) => {
                                        profile_loss_boundary_multiplicity_owned(
                                            &seam,
                                            obs_prev,
                                            obs_this,
                                            owners_equal,
                                            model,
                                            bg,
                                        )?
                                    }
                                    (None, _) => profile_loss_boundary_owned(
                                        &seam,
                                        obs_prev,
                                        obs_this,
                                        owners_equal,
                                        model,
                                    )?,
                                };
                            }
                        }
                        runs.push((owner, piece_profile));
                    }
                } else {
                    runs.push((owner, piece_profile));
                }
                if head.is_none() {
                    head = Some(sequence[..sequence.len().min(149)].to_vec());
                    head_piece = Some((source, start, end, reverse, owner));
                }
                tail = sequence[sequence.len().saturating_sub(149)..].to_vec();
                previous_sequence = Some(sequence);
                previous_piece = Some((source, start, end, reverse, owner));
            }
            copy_runs[copy] = runs;
            copy_novel_charges[copy] = novel_charge;
            tails[copy] = tail;
            let head = head.unwrap_or_default();
            if locus > 0 && !previous_tails[copy].is_empty() && !head.is_empty() {
                let (seam, _) = genome::profile_event_seam(
                    panel,
                    &previous_tails[copy],
                    &head,
                    READ_LENGTH,
                    MAX_FEATURES,
                )?;
                let prev = previous_pieces[copy].expect("piece tracked with tail");
                let next = head_piece.expect("piece tracked with head");
                let cooccurring = spine::junction::segment_pair_gap(
                    &spine::junction::junction_range(prev.0, prev.1, prev.2, prev.3, u32::MAX),
                    &spine::junction::junction_range(next.0, next.1, next.2, next.3, u32::MAX),
                )
                .is_some_and(|gap| gap >= 0);
                match (cooccurring, junction) {
                    // Novel junction with the span index: paid by the reads
                    // that actually cross it (restricted charge, flat
                    // background — the junction-evidence machinery), across
                    // the two pieces' OWNING partitions.
                    (false, Some((span, path_of_source, _partitions))) => {
                        let prev_partition =
                            owner_universe_partition(prev.4, component_locus_to_partition);
                        let next_partition =
                            owner_universe_partition(next.4, component_locus_to_partition);
                        let outcome = span.restricted_charge(
                            &seam,
                            &[(
                                spine::junction::junction_range(
                                    prev.0, prev.1, prev.2, prev.3, prev_partition,
                                ),
                                spine::junction::junction_range(
                                    next.0, next.1, next.2, next.3, next_partition,
                                ),
                            )],
                            sources,
                            path_of_source,
                            &[prev_partition, next_partition],
                            model,
                        )?;
                        total += outcome.charge;
                        census.push(serde_json::json!({
                            "copy": copy,
                            "locus": locus,
                            "kind": "boundary",
                            "left": [prev.0, prev.1, prev.2, prev.3],
                            "right": [next.0, next.1, next.2, next.3],
                            "spanning_reads": outcome.spanning_reads,
                            "events": outcome.events,
                            "span_features": outcome.span_feature_count,
                            "restricted_charge": outcome.charge,
                        }));
                    }
                    // Co-occurring real panel junction (pooled, merged), and
                    // novel junctions under the quarantined pooled model:
                    // the interior charge, per-feature backgrounds when
                    // present, flat bit-identically when absent. A haploid
                    // reference charges the seam sample-side. The observed
                    // sides are the two pieces' OWNING partitions' shares.
                    _ => {
                        if let Some(bg) = &backgrounds {
                            for key in seam.keys() {
                                if !bg.contains(key) {
                                    missing.push(key.clone());
                                }
                            }
                        }
                        let obs_left = owner_routed_obs(
                            routed_equal,
                            prev.4,
                            component_locus_to_partition,
                        );
                        let obs_right = owner_routed_obs(
                            routed_equal,
                            next.4,
                            component_locus_to_partition,
                        );
                        let owners_equal = prev.4 == next.4;
                        // The boundary layer keeps the established
                        // convention UNCHANGED (Fix 2's omission charge
                        // applies at the per-window candidate charges only):
                        // the two-side junction seam's observed side is the
                        // boundary convention's summed form (a
                        // junction-straddling record's anchors distribute
                        // between exactly the two adjacent partitions), and
                        // the charge keeps the baselined sample-side form.
                        total += match (&backgrounds, haploid_reference) {
                            (Some(_), true) => {
                                if owners_equal {
                                    profile_loss_boundary_sample(
                                        &seam,
                                        obs_left,
                                        &EMPTY_ROUTED_OBS,
                                        sample.expect("haploid reference implies sample"),
                                        model,
                                    )?
                                } else {
                                    profile_loss_boundary_sample(
                                        &seam,
                                        obs_left,
                                        obs_right,
                                        sample.expect("haploid reference implies sample"),
                                        model,
                                    )?
                                }
                            }
                            (Some(bg), false) => profile_loss_boundary_multiplicity_owned(
                                &seam,
                                obs_left,
                                obs_right,
                                owners_equal,
                                model,
                                bg,
                            )?,
                            (None, _) => profile_loss_boundary_owned(
                                &seam,
                                obs_left,
                                obs_right,
                                owners_equal,
                                model,
                            )?,
                        };
                    }
                }
            }
            previous_tails[copy] = tails[copy].clone();
            previous_pieces[copy] = pieces.last().copied();
        }
        // Per-locus local charge. DIPLOID references: the copies' owner runs
        // grouped by owning partition; each owner's merged copy profiles
        // charge that owner's routed shares (single owner = the
        // pre-extension merged charge, bit-identically). HAPLOID references
        // (Fix 1 + Fix 2 restructure): ONE site per window — the merged
        // copy-0 profile (all owner runs) charged against the window's full
        // record-once observed map, unbaselined (the omission charge). The
        // per-owner-run splitting would pay the window's omission constant
        // once PER owner — the union-sum pathology re-entered.
        let mut local = 0.0f64;
        if haploid_reference {
            // Fix 1: ONE charge site per window for the haploid reference —
            // the merged copy-0 profile (all owner runs; the interior
            // co-occurring seams merged across owner changes above, the
            // pre-extension shape) charged against the RECORD-ONCE observed
            // map over the profile's owning partitions (each record once —
            // the union-sum form re-pooled record mass).
            let mut merged = Profile::new();
            let mut site_owners: BTreeSet<u32> = BTreeSet::new();
            for (run_owner, run_profile) in &copy_runs[0] {
                merged = merge_profiles(&[&merged, run_profile])?;
                site_owners.insert(*run_owner);
            }
            if let Some(bg) = &backgrounds {
                for key in merged.keys() {
                    if !bg.contains(key) {
                        missing.push(key.clone());
                    }
                }
            }
            let obs_row = site.site_map(
                site_owners
                    .iter()
                    .map(|&owner| owner_universe_partition(owner, component_locus_to_partition)),
            );
            // Fix 2: the omission charge — the window's observed DNA the
            // reference does not spell, at its pure-background Poisson cost;
            // exonerated where the chain explains the feature SOMEWHERE
            // (Ruling 2's chain-level once-per-material form).
            let obs_window = window_obs
                .get(locus)
                .ok_or_else(|| invalid("window observed profile missing for haploid locus"))?;
            local += if window_instances.is_some() && haploid_reference {
                // The corrected Ruling 2: instance-level exoneration — the
                // window's omission covers the observed mass the chain's
                // spelled anchors do not cover at the instances' own
                // positions.
                merged_single_loss_sample_with_omission_instances(
                    &merged,
                    &obs_row,
                    obs_window,
                    window_instances.map(|w| &w.window_records[locus]),
                    window_instances,
                    &chain_spelled,
                    sample.expect("haploid reference implies sample"),
                    model,
                )?
            } else {
                merged_single_loss_sample_with_omission_except(
                    &merged,
                    &obs_row,
                    &chain_explained,
                    obs_window,
                    sample.expect("haploid reference implies sample"),
                    model,
                )?
            };
        } else {
        let mut owners: BTreeSet<u32> = BTreeSet::new();
        for copy in 0..2 {
            for (owner, _) in &copy_runs[copy] {
                owners.insert(*owner);
            }
        }
        for owner in owners {
            let obs = owner_routed_obs(routed_equal, owner, component_locus_to_partition);
            let merged = |copy: usize| -> io::Result<Profile> {
                let mut profile = Profile::new();
                for (run_owner, run_profile) in &copy_runs[copy] {
                    if *run_owner == owner {
                        profile = merge_profiles(&[&profile, run_profile])?;
                    }
                }
                Ok(profile)
            };
            let first = merged(0)?;
            let second = merged(1)?;
            match (&backgrounds, haploid_reference) {
                (Some(bg), true) => {
                    for profile in [&first, &second] {
                        for key in profile.keys() {
                            if !bg.contains(key) {
                                missing.push(key.clone());
                            }
                        }
                    }
                    local += merged_single_loss_sample(
                        &first,
                        obs,
                        sample.expect("haploid reference implies sample"),
                        model,
                    )?;
                }
                (Some(bg), false) => {
                    for profile in [&first, &second] {
                        for key in profile.keys() {
                            if !bg.contains(key) {
                                missing.push(key.clone());
                            }
                        }
                    }
                    local += merged_pair_loss_multiplicity(&first, &second, obs, model, bg)?;
                }
                (None, _) => local += merged_pair_loss(&first, &second, obs, model)?,
            }
        }
        }
        // Fix 2 (the omission charge, reference scope): a DIPLOID
        // reference's union pair profile pays the in-window observed
        // features neither copy explains (the reported score shifts by
        // exactly this addend — the owner-ratified decomposable shift). The
        // haploid branch's addend is inside its merged charge above.
        if !haploid_reference {
            let mut pair_union = Profile::new();
            for copy in 0..2 {
                for (_, run_profile) in &copy_runs[copy] {
                    pair_union = merge_profiles(&[&pair_union, run_profile])?;
                }
            }
            if let Some(sample) = sample {
                local += omission_charge_addend(&pair_union, &window_obs[locus], sample);
            }
        }
        local += copy_novel_charges[0] + copy_novel_charges[1];
        if std::env::var("IMPG_DP_DEBUG").is_ok() {
            eprintln!(
                "[ref-dbg] locus {locus} pieces {} / {} local {local:.2}",
                local_pieces[locus][0].len(),
                local_pieces[locus][1].len()
            );
        }
        total += local;
    }
    // Complete the backgrounds for any charged feature the scan has not
    // measured (piece/seam features beyond the candidate universe), then
    // recurse once: the second pass charges every feature against a
    // measured background and collects nothing.
    if let Some(bg) = &mut backgrounds {
        if !missing.is_empty() {
            let k = panel.syncmer_length_bp() as u64;
            bg.scan_extend(panel, missing, k, model)?;
            return score_reference_m1(
                panel,
                sources,
                territory,
                route_pair,
                routed_equal,
                site,
                window_obs,
                window_instances,
                component_locus_to_partition,
                locus_count,
                locus_offset,
                model,
                junction,
                Some(&mut **bg),
                sample,
                charge_group_shared_once,
                None,
                None,
                spelled_out,
            );
        }
    }
    ensure(total.is_finite(), "nonfinite reference score")?;
    Ok((total, census))
}

fn main() -> io::Result<()> {
    let started = Instant::now();
    let options = Options::parse();
    let (locus_lo, mut locus_hi) = match &options.locus_range {
        Some(range) => {
            ensure(
                range.len() == 2 && range[0] < range[1],
                "locus range needs start < end",
            )?;
            (range[0], range[1])
        }
        None => (0, usize::MAX),
    };
    let rss_budget_bytes: Option<u64> = if options.rss_budget_gib > 0.0 {
        Some((options.rss_budget_gib * (1u64 << 30) as f64) as u64)
    } else {
        None
    };
    let mut rss = genome::PeriodicRssGuard::new(rss_budget_bytes, 1);

    // ---------------------------------------------------------------- inputs
    let input_started = Instant::now();
    let identity = PanelIdentity::read(&options.panel)?;
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    let k = panel.syncmer_length_bp() as u64;
    let graph: routes::Graph = read_json(&options.routes.join("graph.json"))?;
    if graph.panel != identity
        || !graph.generation_complete
        || graph.count_policy != COUNT_POLICY
        || graph.lanes.len() != panel.name_map.path_to_name.len()
    {
        return Err(invalid("incompatible route inventory"));
    }
    let sample_index = sample::SampleIndex::load(&options.sample, &identity)?;
    let histogram = *sample_index
        .stats
        .read_lengths
        .get(&READ_LENGTH)
        .ok_or_else(|| invalid("sample lacks L150 histogram"))?;
    let model = ScoreModel {
        read_length: READ_LENGTH as u64,
        histogram,
        denominator: (READ_LENGTH as u64 * histogram) as f64,
        depth: options.depth,
        background: options.background,
    };
    model.validate()?;
    let lanes = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.clone(), lane.length))
        .collect::<Vec<_>>();
    let sources = routes::Sources::open(&graph.source_paths, lanes.clone())?;
    let axis = genome::load_axis(&options.axis)?;
    let path_of_name: HashMap<String, usize> = panel
        .name_map
        .name_to_path
        .iter()
        .map(|(name, &path)| (name.clone(), path as usize))
        .collect();
    let path_of_source: Vec<usize> = graph
        .lanes
        .iter()
        .map(|lane| {
            path_of_name
                .get(&lane.name)
                .copied()
                .ok_or_else(|| invalid("route lane absent from panel"))
        })
        .collect::<io::Result<_>>()?;
    let source_of_name: HashMap<&str, usize> = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.as_str(), lane.id))
        .collect();
    // The territory index's rc-symmetric augmentation fetches each universe
    // row's path sequence through the route graph's sources (the same
    // sequence binding the spelled side spells with).
    let fetch_path_seq = |path_idx: usize, lo: u64, hi: u64| -> io::Result<Vec<u8>> {
        let name = &panel.name_map.path_to_name[path_idx];
        let source = *source_of_name.get(name.as_str()).ok_or_else(|| {
            io::Error::other(format!(
                "panel path {name} absent from the route graph's sources"
            ))
        })?;
        let hi = hi.min(sources.lanes[source].1);
        if lo >= hi {
            return Ok(Vec::new());
        }
        sources.fetch(source, lo, hi)
    };
    let input_seconds = input_started.elapsed().as_secs_f64();
    rss_probe(&mut rss, "input_load")?;
    let mut dp_summary_global: Option<serde_json::Value> = None;

    // ------------------------------------------------- slice candidate domain
    let domain_started = Instant::now();
    // The window-domain extension (owner-ruled Policy A + minimal universe
    // extension). Flag-gated: without it every domain, universe and share
    // is bit-identical to the pre-extension model (the regression guard's
    // baseline).
    ensure(
        !options.window_domain_extension || options.routing_universe == "component",
        "the window-domain extension requires the component routing universe",
    )?;
    let (component_axis, mut component_ranges) = genome::load_bed_axis_partitions(
        &axis,
        &options.bed_directory,
        &lanes,
        &options.component,
        None,
        None,
    )?;
    let target = graph
        .lanes
        .iter()
        .find(|lane| lane.name == options.component)
        .ok_or_else(|| invalid("component is not a route lane"))?;
    let component_suffix = target
        .name
        .splitn(3, '#')
        .nth(2)
        .ok_or_else(|| invalid("route lane lacks component suffix"))?;
    let component_sources = graph
        .lanes
        .iter()
        .filter(|lane| lane.name.splitn(3, '#').nth(2) == Some(component_suffix))
        .map(|lane| lane.id)
        .collect::<BTreeSet<_>>();
    let extension = if options.window_domain_extension {
        Some(genome::build_window_domain_extension(
            &axis,
            &options.bed_directory,
            &lanes,
            &options.component,
            &component_sources,
        )?)
    } else {
        None
    };
    if let Some(extension) = &extension {
        eprintln!(
            "[window-domain] extension: {} pure-new groups, {} dual-role groups, \
             {} added row pairs, duplicate row identities {}, overlapping row pairs {}",
            extension.pure_new.len(),
            extension.dual_role.len(),
            extension.added.iter().map(Vec::len).sum::<usize>() / 2,
            extension.duplicate_row_identities,
            extension.overlapping_row_pairs,
        );
    }
    let mut bed_territory: Vec<Vec<SourceRange>> = component_ranges.clone();
    // Forward source-path completion runs on the ANCHOR rows only (semantic
    // choice: the completed domain stays the axis rows' own path geometry;
    // the extension's added rows enter verbatim AFTER completion, so the
    // candidate domain is exactly Policy A's anchor BED + overlapping
    // component-family rows).
    genome::complete_forward_source_paths(&mut component_ranges, &component_sources)?;
    if let Some(extension) = &extension {
        for (locus, added) in extension.added.iter().enumerate() {
            component_ranges[locus].extend(added.iter().cloned());
            bed_territory[locus].extend(added.iter().cloned());
        }
    }
    let mut traversals = component_ranges
        .into_iter()
        .map(|ranges| {
            ranges
                .into_iter()
                .map(genome::SpanningTraversal::single)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    genome::retain_native_endpoint_candidates(&mut traversals, target.id, target.length)?;
    if locus_hi == usize::MAX {
        locus_hi = traversals.len();
    }
    ensure(
        locus_lo < locus_hi && locus_hi <= traversals.len(),
        "locus range outside component partitions",
    )?;
    // Forward single-segment candidates of the slice (deduplicated by
    // interval; reverse candidates are profiled by the mirrored geometry).
    let mut forward_candidates: Vec<(usize, Vec<(usize, u64, u64)>, String)> = Vec::new();
    let mut reverse_sample: Vec<(usize, Vec<(usize, u64, u64)>, String)> = Vec::new();
    let mut seen_intervals: BTreeSet<(usize, u64, u64)> = BTreeSet::new();
    for (locus, locus_traversals) in traversals[locus_lo..locus_hi].iter().enumerate() {
        for traversal in locus_traversals {
            if traversal.segments.len() != 1 {
                continue;
            }
            let range = &traversal.segments[0];
            let segment = (range.source, range.start, range.end);
            if range.reverse {
                if reverse_sample.len() < 4 {
                    reverse_sample.push((locus_lo + locus, vec![segment], traversal.identity.clone()));
                }
            } else if seen_intervals.insert(segment) {
                forward_candidates.push((locus_lo + locus, vec![segment], traversal.identity.clone()));
            }
        }
    }
    ensure(!forward_candidates.is_empty(), "slice has no forward candidates")?;
    let slice_loci = locus_hi - locus_lo;
    let domain_seconds = domain_started.elapsed().as_secs_f64();
    rss_probe(&mut rss, "candidate_domain")?;

    // ------------------------------------------------- sample record pass
    let records_started = Instant::now();
    // Sample record re-derivation WITH the within-read adjacency chains
    // (one streaming pass over the base reads; the chains feed the junction
    // span index — the junction-evidence restriction's per-run derivation).
    let mut reads_batch: Vec<Vec<u8>> = Vec::new();
    stream_fastq(&options.reads, |seq| {
        reads_batch.push(seq.to_vec());
        Ok(())
    })?;
    ensure(
        reads_batch.len() as u64 == sample_index.stats.reads,
        "reads count mismatch",
    )?;
    rss_probe(&mut rss, "sample_read_buffer")?;
    let read_chains = spine::junction::derive_read_chains(
        &panel,
        &reads_batch,
        &sample_index.stats,
        &mut rss,
    )?;
    let records: BTreeMap<Vec<u64>, u64> = read_chains
        .record_tokens
        .iter()
        .zip(read_chains.record_counts.iter())
        .map(|(tokens, &count)| (tokens.clone(), count))
        .collect();
    let mut spot_checked = 0usize;
    for (record, &multiplicity) in records.iter() {
        if spot_checked >= 512 {
            break;
        }
        spot_checked += 1;
        ensure(
            sample_index.counts.count(record)? >= multiplicity,
            "derived record missing from stored index",
        )?;
    }
    let records_seconds = records_started.elapsed().as_secs_f64();
    rss_probe(&mut rss, "sample_record_rederivation")?;

    // ------------------------------------------------- territory index
    let territory_started = Instant::now();
    let universe = load_universe(
        &axis,
        &options.bed_directory,
        &lanes,
        &options.component,
        &options.routing_universe,
        traversals.len(),
        extension.as_ref(),
    )?;
    // The strand-asymmetry fix's measured-yield gate: when
    // IMPG_STRANDFIX_YIELD names an output path, a SECOND territory index is
    // built WITHOUT the rc-frame augmentation and every routed record is
    // routed against both — the completeness delta (placements gained/lost)
    // is the fix's measured yield over the routed records.
    let strand_yield_path = std::env::var("IMPG_STRANDFIX_YIELD").ok();
    let yield_index = strand_yield_path
        .as_ref()
        .map(|_| {
            build_territory_index(&panel, &universe, &path_of_name, k, &fetch_path_seq, false)
        })
        .transpose()?;
    let territory = build_territory_index(&panel, &universe, &path_of_name, k, &fetch_path_seq, true)?;
    let territory_seconds = territory_started.elapsed().as_secs_f64();
    rss_probe(&mut rss, "territory_index")?;

    // ------------------------------------------------- routing pass
    let routing_started = Instant::now();
    let mut slice_partitions: BTreeSet<u32> = (locus_lo..locus_hi)
        .filter_map(|locus| {
            let partition = universe.component_locus_to_partition[locus];
            (partition != u32::MAX).then_some(partition)
        })
        .collect();
    // The window-domain extension's owner partitions: every owning universe
    // partition of an added candidate row in the slice's domains charges its
    // own routed share, so it joins the slice's accumulation set.
    if let Some(extension) = &extension {
        for locus_traversals in &traversals[locus_lo..locus_hi] {
            for traversal in locus_traversals {
                for segment in &traversal.segments {
                    if segment.start < segment.end {
                        slice_partitions.insert(segment.partition as u32);
                    }
                }
            }
        }
        // The added rows' owning partitions even where no traversal uses
        // them: their territory is inside the slice's windows, so their
        // routed shares belong to the windows' observed profiles (the
        // omission charge's universes must see the windows' full observed
        // mass, gap material included).
        for locus in locus_lo..locus_hi {
            for row in &extension.added[locus] {
                slice_partitions.insert(owner_universe_partition(
                    row.partition as u32,
                    &universe.component_locus_to_partition,
                ));
            }
        }
    }
    let record_vec: Vec<(&Vec<u64>, u64)> = records.iter().map(|(r, &w)| (r, w)).collect();
    let yield_records_total = AtomicU64::new(0);
    let yield_slice_records_total = AtomicU64::new(0);
    let yield_records_gained = AtomicU64::new(0);
    let yield_slice_records_gained = AtomicU64::new(0);
    let yield_placements_gained = AtomicU64::new(0);
    let yield_placements_lost = AtomicU64::new(0);
    let routed_records: Vec<RoutedRecord> = record_vec
        .par_iter()
        .filter(|(tokens, _)| !tokens.is_empty())
        .filter_map(|(tokens, multiplicity)| {
            let multiplicity = *multiplicity;
            let anchors = decode_tokens(tokens).ok()?;
            let (occurrences, total_occurrences, forward_positions, reverse_positions) =
                route_record(&anchors, &territory, k, false);
            if let Some(yield_index) = &yield_index {
                let (_, _, old_forward, old_reverse) =
                    route_record(&anchors, yield_index, k, false);
                let new_set: HashSet<(usize, u64)> = forward_positions
                    .iter()
                    .chain(reverse_positions.iter())
                    .copied()
                    .collect();
                let old_set: HashSet<(usize, u64)> = old_forward
                    .iter()
                    .chain(old_reverse.iter())
                    .copied()
                    .collect();
                let gained = new_set.difference(&old_set).count() as u64;
                let lost = old_set.difference(&new_set).count() as u64;
                yield_placements_gained.fetch_add(gained, Ordering::Relaxed);
                yield_placements_lost.fetch_add(lost, Ordering::Relaxed);
                let touches_slice = occurrences.keys().any(|p| slice_partitions.contains(p));
                yield_records_total.fetch_add(1, Ordering::Relaxed);
                if touches_slice {
                    yield_slice_records_total.fetch_add(1, Ordering::Relaxed);
                }
                if gained > 0 {
                    yield_records_gained.fetch_add(1, Ordering::Relaxed);
                    if touches_slice {
                        yield_slice_records_gained.fetch_add(1, Ordering::Relaxed);
                    }
                }
            }
            (!occurrences.is_empty()).then(|| {
                RoutedRecord {
                    tokens: (*tokens).clone(),
                    multiplicity,
                    occurrences,
                    total_occurrences,
                    forward_positions,
                    reverse_positions,
                }
            })
        })
        .collect();
    let records_touching_universe = routed_records.len() as u64;
    // The strand-fix yield summary (env-gated): the completeness delta of the
    // rc-symmetric index over the routed records — the fraction of records
    // (and of slice-touching records) that GAINED placements they lacked,
    // plus the gross gained/lost placement counts.
    if let Some(yield_path) = &strand_yield_path {
        let total = yield_records_total.load(Ordering::Relaxed);
        let slice_total = yield_slice_records_total.load(Ordering::Relaxed);
        let summary = serde_json::json!({
            "records_routed": total,
            "slice_records_routed": slice_total,
            "records_gained_placements": yield_records_gained.load(Ordering::Relaxed),
            "slice_records_gained_placements": yield_slice_records_gained.load(Ordering::Relaxed),
            "placements_gained": yield_placements_gained.load(Ordering::Relaxed),
            "placements_lost": yield_placements_lost.load(Ordering::Relaxed),
            "record_gain_fraction": if total > 0 {
                yield_records_gained.load(Ordering::Relaxed) as f64 / total as f64
            } else {
                0.0
            },
            "slice_record_gain_fraction": if slice_total > 0 {
                yield_slice_records_gained.load(Ordering::Relaxed) as f64 / slice_total as f64
            } else {
                0.0
            },
        });
        let text = serde_json::to_string_pretty(&summary)?;
        std::fs::write(yield_path, text.clone() + "\n")?;
        eprintln!("[strand-yield] {text}");
    }
    // The twin index is dead weight after the routing pass.
    drop(yield_index);
    // The router-placement diagnosis (env-gated): traces the targeted
    // records' placements through route_record's search and force-verifies
    // the walks at the spelled loci — then exits before any downstream stage.
    if let Ok(spec_path) = std::env::var("IMPG_ROUTER_PLACEMENT_DIAG") {
        router_placement_diagnosis(&spec_path, &routed_records, &territory, &panel, k)?;
        std::process::exit(0);
    }
    let route_search_seconds = routing_started.elapsed().as_secs_f64();
    rss_probe(&mut rss, "route_search")?;

    // Equal-share and occurrence-weighted accumulation for slice partitions.
    let accumulation_started = Instant::now();
    let mut routed_equal: RoutedObs = HashMap::new();
    let mut routed_occurrence: RoutedObs = HashMap::new();
    // The LEGACY sample-side pooled support per feature (decompose-table
    // diagnostic columns only since the genome-wide decision 2026-09-24;
    // no charging path reads it): the sample's TOTAL realized multiplicity
    // of each feature across the routing universe over slice-touching
    // records — the pooled C whose component-scoped complementary share was
    // the Step-B background. Only built when a decomposition (or a run that
    // consumes the diagnostics) is requested.
    let mut pooled_equal: HashMap<FeatureKey, f64> = HashMap::new();
    let build_pooled = options.spine_phasing
        || options.decompose_selected.is_some()
        || options.routing_diagnostics_only;
    let mut records_routed_to_slice = 0u64;
    // Reconciliation-gate instrumentation (owner ruling): per-slice-partition
    // total routed-share mass (the before/after diff baseline of the axis
    // partitions' routed shares under the extended universe), and the
    // t_r_genome buckets of the records realizing the extension's new
    // territory (the gap-realizing records; whether gap material carries
    // genome-wide repeat mass the genome pass misses).
    let mut routed_equal_partition_mass: BTreeMap<u32, f64> = BTreeMap::new();
    let extension_axis_partitions = universe
        .partitions
        .saturating_sub(extension.as_ref().map_or(0, |ext| ext.pure_new.len()));
    // The gap-realizing records (those touching an appended extension
    // partition), identified by token multiset for the genome pass's
    // t_r_genome bucketing (built only when the genome pass runs).
    let extension_record_tokens: std::sync::Arc<std::collections::HashSet<Vec<u64>>> =
        if extension.is_some() {
            std::sync::Arc::new(
                routed_records
                    .iter()
                    .filter(|record| {
                        record
                            .occurrences
                            .keys()
                            .any(|&partition| partition as usize >= extension_axis_partitions)
                    })
                    .map(|record| record.tokens.clone())
                    .collect(),
            )
        } else {
            std::sync::Arc::new(std::collections::HashSet::new())
        };
    let mut extension_records = 0u64;
    for record in &routed_records {
        let touched: BTreeSet<u32> = record.occurrences.keys().copied().collect();
        let total_occurrences = record.total_occurrences;
        let mut touches_slice = false;
        for partition in record.occurrences.keys() {
            if slice_partitions.contains(partition) {
                touches_slice = true;
            }
        }
        if !touches_slice {
            continue;
        }
        records_routed_to_slice += 1;
        let subwalks = enumerate_subwalks(&record.tokens)?;
        let denominator = touched.len() as f64;
        if record
            .occurrences
            .keys()
            .any(|&partition| partition as usize >= extension_axis_partitions)
        {
            extension_records += 1;
        }
        if build_pooled {
            let multiplicity = record.multiplicity as f64;
            for feature in &subwalks {
                *pooled_equal.entry(feature.clone()).or_default() += multiplicity;
            }
        }
        for (partition, count) in &record.occurrences {
            if !slice_partitions.contains(partition) {
                continue;
            }
            let equal_share = record.multiplicity as f64 / denominator;
            let occurrence_share = if total_occurrences > 0 {
                record.multiplicity as f64 * (*count as f64) / total_occurrences as f64
            } else {
                0.0
            };
            let equal_map = routed_equal.entry(*partition).or_default();
            let occ_map = routed_occurrence.entry(*partition).or_default();
            let mut partition_mass = 0.0f64;
            for feature in &subwalks {
                *equal_map.entry(feature.clone()).or_default() += equal_share;
                *occ_map.entry(feature.clone()).or_default() += occurrence_share;
                partition_mass += equal_share;
            }
            *routed_equal_partition_mass.entry(*partition).or_default() += partition_mass;
        }
    }
    let accumulation_seconds = accumulation_started.elapsed().as_secs_f64();
    rss_probe(&mut rss, "share_accumulation")?;

    // The record-once site-map builder (Fix 1) and the per-window owner sets
    // (the axis partition + every added row's owning partition) driving the
    // in-window observed profiles (Fix 2's omission universes).
    let site_observed = SiteObserved::build(&routed_records, &slice_partitions);
    let window_owner_sets: Vec<BTreeSet<u32>> = (locus_lo..locus_hi)
        .map(|component_locus| {
            let mut owners: BTreeSet<u32> = BTreeSet::new();
            let axis = universe.component_locus_to_partition[component_locus];
            if axis != u32::MAX {
                owners.insert(axis);
            }
            if let Some(extension) = &extension {
                for row in &extension.added[component_locus] {
                    owners.insert(owner_universe_partition(
                        row.partition as u32,
                        &universe.component_locus_to_partition,
                    ));
                }
            }
            owners
        })
        .collect();
    // The windows' IN-WINDOW observed profiles (the omission charge's
    // universes, owner ruling 2026-09-24): the per-feature positional
    // attribution of the observed mass physically inside each window's
    // territory. The SAME pass also emits the windows' observed INSTANCE
    // structure (the corrected Ruling 2's observed side: per feature, the
    // contributing records' in-window placement spans + equal shares).
    // The instance-probe diagnostic's feature list (the pooled
    // decomposition's top truth-advantage features; JSON array of arrays).
    // Operational diagnostic input, not a model constant.
    let instance_probe_features: Vec<FeatureKey> = std::env::var("IMPG_INSTANCE_PROBE_FEATURES")
        .ok()
        .map(|raw| serde_json::from_str::<Vec<FeatureKey>>(&raw))
        .transpose()?
        .unwrap_or_default();
    let mut instance_spans: RecordPlacementSpans = Vec::new();
    let mut window_records: WindowRecordLists = Vec::new();
    let window_obs: Vec<HashMap<FeatureKey, f64>> = build_in_window_obs(
        &routed_records,
        &territory,
        k,
        &window_owner_sets,
        Some(&mut instance_spans),
        Some(&mut window_records),
    )?;
    // The records' equal shares (the mass each record's index carries; the
    // same share the observed maps accumulated).
    let record_shares: Vec<f64> = routed_records
        .iter()
        .map(|record| record.multiplicity as f64 / record.occurrences.len() as f64)
        .collect();
    let window_instances = InstanceStructure {
        record_spans: instance_spans,
        record_shares,
        window_records,
    };
    rss_probe(&mut rss, "site_observed_maps")?;

    // --------------------------------------------- genome-universe placement pass
    // (the haploid GENOME-WIDE cross-support background — the owner-decided
    // candidate family 1, "genome-wide touched set"): per-record touched
    // counts across ALL components' partitions and the per-feature
    // background map Sum_r m_r*(1 - 1/t_r_genome) over every record
    // realizing each feature. Computed per run (zero-cache doctrine); built
    // only when a consumer is present (phasing or decomposition — the same
    // condition as the pooled diagnostic map). The OBSERVED shares stay
    // component-scoped exactly as before: this pass changes only the
    // haploid backgrounds. When the run's routing universe is already the
    // genome, the main placement pass IS the genome pass (no second
    // territory index is built).
    const GENOME_PASS_CHUNK: usize = 8192;
    let genome_started = Instant::now();
    let mut genome_background: HashMap<FeatureKey, f64> = HashMap::new();
    let mut genome_partitions = 0usize;
    let mut genome_records_placed = 0u64;
    let mut genome_territory_seconds = 0.0f64;
    let mut genome_seconds = 0.0f64;
    // t_r_genome buckets of the gap-realizing records (the reconciliation
    // gate's report; component-universe runs only, genome pass present).
    let mut extension_t_r_genome_buckets = [0u64; 8];
    let mut extension_records_placed = 0u64;
    if build_pooled {
        let chunk_maps: Vec<(HashMap<FeatureKey, f64>, u64, [u64; 8], u64)> =
            if options.routing_universe == "genome" {
                genome_partitions = universe.partitions;
                routed_records
                    .par_chunks(GENOME_PASS_CHUNK)
                    .map(|chunk| {
                        let mut map = HashMap::new();
                        let mut placed = 0u64;
                        for record in chunk {
                            if record.occurrences.is_empty() {
                                continue;
                            }
                            placed += 1;
                            accumulate_cross_support(
                                &mut map,
                                &record.tokens,
                                record.multiplicity,
                                record.occurrences.len(),
                            );
                        }
                        (map, placed, [0u64; 8], 0u64)
                    })
                    .collect()
            } else {
                let genome_territory_started = Instant::now();
                let genome_universe = load_universe(
                    &axis,
                    &options.bed_directory,
                    &lanes,
                    &options.component,
                    "genome",
                    traversals.len(),
                    // The genome-universe placement pass stays EXACTLY as is
                    // (owner ruling): the extension never enters it, so gap
                    // records keep t_r_genome = 0 -> beta = base.
                    None,
                )?;
                genome_partitions = genome_universe.partitions;
                let genome_index = build_territory_index(
                    &panel,
                    &genome_universe,
                    &path_of_name,
                    k,
                    &fetch_path_seq,
                    true,
                )?;
                genome_territory_seconds =
                    genome_territory_started.elapsed().as_secs_f64();
                record_vec
                    .par_chunks(GENOME_PASS_CHUNK)
                    .map(|chunk| {
                        let mut map = HashMap::new();
                        let mut placed = 0u64;
                        let mut gap_buckets = [0u64; 8];
                        let mut gap_placed = 0u64;
                        for &(tokens, multiplicity) in chunk {
                            if tokens.is_empty() {
                                continue;
                            }
                            let anchors = match decode_tokens(tokens) {
                                Ok(anchors) => anchors,
                                Err(_) => continue,
                            };
                            let (occurrences, _, _, _) =
                                route_record(&anchors, &genome_index, k, false);
                            if occurrences.is_empty() {
                                continue;
                            }
                            placed += 1;
                            if extension_record_tokens.contains(tokens) {
                                gap_placed += 1;
                                let bucket = match occurrences.len() {
                                    1 => 1,
                                    2 => 2,
                                    3..=4 => 3,
                                    5..=8 => 4,
                                    9..=16 => 5,
                                    17..=32 => 6,
                                    _ => 7,
                                };
                                gap_buckets[bucket] += 1;
                            }
                            accumulate_cross_support(
                                &mut map,
                                tokens,
                                multiplicity,
                                occurrences.len(),
                            );
                        }
                        (map, placed, gap_buckets, gap_placed)
                    })
                    .collect()
            };
        for (map, placed, gap_buckets, gap_placed) in chunk_maps {
            genome_records_placed += placed;
            extension_records_placed += gap_placed;
            for (bucket, count) in gap_buckets.iter().enumerate() {
                extension_t_r_genome_buckets[bucket] += count;
            }
            for (feature, value) in map {
                *genome_background.entry(feature).or_default() += value;
            }
        }
        genome_seconds = genome_started.elapsed().as_secs_f64();
        rss_probe(&mut rss, "genome_background_pass")?;
    }

    // Boundary-crossing miss diagnostic: min-anchor search vs every-anchor
    // search on a deterministic sample of records.
    let mut multi_anchor_checked = 0usize;
    let mut multi_anchor_touched_set_mismatches = 0usize;
    for (index, record) in routed_records.iter().enumerate() {
        if index % 1024 != 0 || multi_anchor_checked >= 512 {
            continue;
        }
        multi_anchor_checked += 1;
        let anchors = decode_tokens(&record.tokens)?;
        let (exact, _, _, _) = route_record(&anchors, &territory, k, true);
        let keys_exact: BTreeSet<u32> = exact.keys().copied().collect();
        let keys_fast: BTreeSet<u32> = record.occurrences.keys().copied().collect();
        if keys_exact != keys_fast {
            multi_anchor_touched_set_mismatches += 1;
        }
    }
    let routing_seconds = routing_started.elapsed().as_secs_f64();
    let routing_summary = serde_json::json!({
        "universe": options.routing_universe,
        "universe_partitions": universe.partitions,
        "universe_territory_intervals": universe.rows.len(),
        "territory_entries": territory.entries.len(),
        "records_total": records.len() as u64,
        "records_touching_universe": records_touching_universe,
        "records_routed_to_slice": records_routed_to_slice,
        "routed_equal_feature_partitions": routed_equal.len(),
        "routed_occurrence_feature_partitions": routed_occurrence.len(),
        "multi_anchor_checked": multi_anchor_checked,
        "multi_anchor_touched_set_mismatches": multi_anchor_touched_set_mismatches,
        // Reconciliation-gate instrumentation: the extension's shape and the
        // per-slice-partition routed-share mass (the before/after axis-share
        // diff baseline; the axis partitions' ids are the component-local
        // locus indices).
        "window_domain": extension.as_ref().map(|ext| serde_json::json!({
            "pure_new_groups": ext.pure_new.len(),
            "dual_role_groups": ext.dual_role.len(),
            "added_row_pairs_total": ext.added.iter().map(Vec::len).sum::<usize>() / 2,
            "added_rows_per_locus": ext.added.iter().map(|rows| rows.len() / 2).collect::<Vec<_>>(),
            "duplicate_row_identities": ext.duplicate_row_identities,
            "overlapping_row_pairs": ext.overlapping_row_pairs,
            "extension_records": extension_records,
            "extension_t_r_genome_buckets": if build_pooled {
                serde_json::json!({
                    "placed": extension_records_placed,
                    "t1": extension_t_r_genome_buckets[1],
                    "t2": extension_t_r_genome_buckets[2],
                    "t3_4": extension_t_r_genome_buckets[3],
                    "t5_8": extension_t_r_genome_buckets[4],
                    "t9_16": extension_t_r_genome_buckets[5],
                    "t17_32": extension_t_r_genome_buckets[6],
                    "t33plus": extension_t_r_genome_buckets[7],
                })
            } else {
                serde_json::Value::Null
            },
        })),
        "routed_equal_partition_mass": routed_equal_partition_mass,
    });
    if options.routing_diagnostics_only {
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "model": "routing-diagnostics-v1",
                "component": options.component,
                "locus_range": [locus_lo, locus_hi],
                "routing_universe": options.routing_universe,
                "rss_peak_bytes": rss.peak_bytes(),
                "routing": routing_summary,
                "genome_background": {
                    "partitions": genome_partitions,
                    "records_placed": genome_records_placed,
                    "features": genome_background.len(),
                },
                "stage_wall_seconds": {
                    "input_load": input_seconds,
                    "candidate_domain": domain_seconds,
                    "sample_record_rederivation": records_seconds,
                    "territory_index": territory_seconds,
                    "route_search": route_search_seconds,
                    "share_accumulation": accumulation_seconds,
                    "genome_background_pass": genome_seconds,
                },
            }))?
        );
        return Ok(());
    }

    // ------------------------------------------------- within-read adjacency index
    // (the junction-evidence restriction): built for the spine paths and the
    // standalone census mode; the quarantined beam path never builds it.
    let span_index: Option<spine::junction::JunctionSpanIndex> =
        if options.spine || options.spine_census_junctions.is_some() {
            let junction_started = Instant::now();
            let touched_by_tokens: HashMap<Vec<u64>, BTreeSet<u32>> = routed_records
                .iter()
                .map(|record| {
                    (
                        record.tokens.clone(),
                        record.occurrences.keys().copied().collect(),
                    )
                })
                .collect();
            let source_of_path: Vec<usize> = {
                let mut map = vec![usize::MAX; panel.name_map.path_to_name.len()];
                for (source, &path) in path_of_source.iter().enumerate() {
                    ensure(map[path] == usize::MAX, "lane paths are not unique")?;
                    map[path] = source;
                }
                map
            };
            let index = spine::junction::build_junction_span_index(
                &territory,
                &sources,
                &source_of_path,
                &reads_batch,
                &read_chains,
                &touched_by_tokens,
                k,
                &mut rss,
            )?;
            eprintln!(
                "[junction] span index built in {:.1}s",
                junction_started.elapsed().as_secs_f64()
            );
            Some(index)
        } else {
            None
        };

    // ------------------------------------------------- standalone census mode
    if let Some(census_path) = &options.spine_census_junctions {
        let span = span_index
            .as_ref()
            .ok_or_else(|| invalid("census mode requires the span index"))?;
        let census_input: serde_json::Value = read_json(census_path)?;
        let rows: Vec<serde_json::Value> = census_input["junctions"]
            .as_array()
            .ok_or_else(|| invalid("census input needs a junctions array"))?
            .iter()
            .map(|row| {
                let field = |name: &str| -> io::Result<(usize, u64, u64, bool)> {
                    let node = &row[name];
                    Ok((
                        node["source"].as_u64().ok_or_else(|| invalid("source"))? as usize,
                        node["start"].as_u64().ok_or_else(|| invalid("start"))?,
                        node["end"].as_u64().ok_or_else(|| invalid("end"))?,
                        node["reverse"].as_bool().ok_or_else(|| invalid("reverse"))?,
                    ))
                };
                let left = field("left")?;
                let right = field("right")?;
                let partitions: Vec<u32> = row["partitions"]
                    .as_array()
                    .map(|list| {
                        list.iter()
                            .filter_map(|value| value.as_u64().map(|value| value as u32))
                            .collect()
                    })
                    .unwrap_or_default();
                let label = row["label"].as_str().unwrap_or("").to_string();
                let left_range = spine::junction::junction_range(
                    left.0, left.1, left.2, left.3, partitions.first().copied().unwrap_or(u32::MAX),
                );
                let right_range = spine::junction::junction_range(
                    right.0,
                    right.1,
                    right.2,
                    right.3,
                    partitions.last().copied().unwrap_or(u32::MAX),
                );
                let crossing =
                    span.crossing_reads(&left_range, &right_range, &sources, &path_of_source)?;
                let trace = if row.get("trace").and_then(|t| t.as_bool()).unwrap_or(false) {
                    Some(span.trace_junction(&left_range, &right_range, &path_of_source))
                } else {
                    None
                };
                // The junction's seam profile (the same L149 juxtaposition
                // prediction the seam machinery scores) and both charges:
                // restricted (crossing reads only) and pooled (the OLD
                // model's charge, for the before/after comparison).
                let flank_memo: FlankMemo = std::sync::Mutex::new(HashMap::new());
                let tail = segment_flank(&sources, &flank_memo, &left_range, false, READ_LENGTH - 1)?;
                let head = segment_flank(&sources, &flank_memo, &right_range, true, READ_LENGTH - 1)?;
                let (profile, _) =
                    genome::profile_event_seam(&panel, &tail, &head, READ_LENGTH, MAX_FEATURES)?;
                let (restricted, pooled) = if partitions.is_empty() {
                    (None, None)
                } else {
                    let restricted = span
                        .restricted_charge(
                            &profile,
                            &[(left_range.clone(), right_range.clone())],
                            &sources,
                            &path_of_source,
                            &partitions,
                            &model,
                        )?
                        .charge;
                    let mut pooled_obs: HashMap<FeatureKey, f64> = HashMap::new();
                    for partition in &partitions {
                        let empty_partition: HashMap<FeatureKey, f64> = HashMap::new();
                        for (feature, share) in routed_equal
                            .get(partition)
                            .unwrap_or(&empty_partition)
                        {
                            *pooled_obs.entry(feature.clone()).or_default() += share;
                        }
                    }
                    let empty: HashMap<FeatureKey, f64> = HashMap::new();
                    let pooled =
                        profile_loss_boundary(&profile, &pooled_obs, &empty, &model)?;
                    (Some(restricted), Some(pooled))
                };
                Ok::<_, io::Error>(serde_json::json!({
                    "label": label,
                    "left": [left.0, left.1, left.2, left.3],
                    "right": [right.0, right.1, right.2, right.3],
                    "partitions": partitions,
                    "spanning_reads": crossing.distinct_reads,
                    "events": crossing.events,
                    "read_gaps": crossing.read_gaps,
                    "restricted_charge": restricted,
                    "pooled_charge_comparison": pooled,
                    "trace": trace,
                }))
            })
            .collect::<io::Result<_>>()?;
        let trace_reads: Option<Vec<serde_json::Value>> =
            census_input["trace_reads"]
                .as_array()
                .map(|list| {
                    list.iter()
                        .filter_map(|value| value.as_u64().map(|value| value as u32))
                        .collect::<Vec<u32>>()
                })
                .map(|reads| span.trace_reads(&reads));
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "model": "junction-span-census-v1",
                "junction_span_index": span.stats,
                "trace_reads": trace_reads,
                "census": rows,
            }))?
        );
        return Ok(());
    }

    // ------------------------------------------------- local-first spine
    if options.spine {
        let reference_specs: Vec<(String, PathBuf, PathBuf)> = options
            .reference_pairs
            .iter()
            .filter_map(|spec| {
                let mut parts = spec.splitn(3, ':');
                let label = parts.next()?.to_string();
                let path_a = PathBuf::from(parts.next()?);
                let path_b = PathBuf::from(parts.next()?);
                Some((label, path_a, path_b))
            })
            .collect();
        ensure(
            options.reference_pairs.len() == reference_specs.len(),
            "--reference-pair needs LABEL:A.json:B.json"
        )?;
        let probe_specs: Vec<(String, PathBuf)> = options
            .probe_routes
            .iter()
            .filter_map(|spec| {
                let mut parts = spec.splitn(2, ':');
                let label = parts.next()?.to_string();
                let path = PathBuf::from(parts.next()?);
                Some((label, path))
            })
            .collect();
        ensure(
            options.probe_routes.len() == probe_specs.len(),
            "--probe-route needs LABEL:ROUTE.json"
        )?;
        let decompose_request = match (&options.decompose_selected, &options.decompose_out) {
            (Some(selected), Some(out)) => Some((selected.clone(), out.clone())),
            (None, None) => None,
            _ => {
                return Err(invalid(
                    "--decompose-selected and --decompose-out must be given together"
                ))
            }
        };
        let spine_summary = spine::run_local_first_spine(
            &panel,
            &sources,
            &graph,
            &options.routes,
            &path_of_source,
            &component_axis[locus_lo..locus_hi],
            traversals[locus_lo..locus_hi].to_vec(),
            &bed_territory,
            k,
            &routed_equal,
            &site_observed,
            &window_obs,
            Some(&window_instances),
            &instance_probe_features,
            &universe.component_locus_to_partition,
            locus_lo,
            target.id,
            target.length,
            &model,
            &sample_index.counts,
            &reference_specs,
            options.spine_sweep_only,
            options.spine_phasing,
            decompose_request.as_ref(),
            // The LEGACY component-scoped pooled support — decompose-table
            // diagnostic columns only (the old-vs-new delta distribution);
            // no charging path reads it.
            if build_pooled {
                Some(&pooled_equal)
            } else {
                None
            },
            // The haploid track's CHARGED background: the genome-wide
            // cross-support map (beta_f = base + Sum_r m_r*(1 - 1/t_r_genome)).
            if build_pooled {
                Some(&genome_background)
            } else {
                None
            },
            span_index
                .as_ref()
                .ok_or_else(|| invalid("spine requires the junction span index"))?,
            std::path::Path::new("spine-rescore"),
            &mut rss,
            &probe_specs,
        )?;
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "model": "local-first-spine-v1",
                "spine_phasing": options.spine_phasing,
                "component": options.component,
                "locus_range": [locus_lo, locus_hi],
                "depth": options.depth,
                "background": options.background,
                "routing_universe": options.routing_universe,
                "rss_peak_bytes": rss.peak_bytes(),
                "top_level_wall_seconds": {
                    "input_load": input_seconds,
                    "candidate_domain": domain_seconds,
                    "sample_record_rederivation": records_seconds,
                    "territory_index": territory_seconds,
                    "route_search": route_search_seconds,
                    "share_accumulation": accumulation_seconds,
                    "genome_background_pass": genome_seconds,
                },
                "genome_background": {
                    "partitions": genome_partitions,
                    "records_placed": genome_records_placed,
                    "territory_seconds": genome_territory_seconds,
                    "features": genome_background.len(),
                },
                "spine": spine_summary,
                "routing": routing_summary,
            }))?
        );
        return Ok(());
    }

    // ------------------------------------------------- routed-evidence DP stage
    if let Some(evidence_path) = &options.refinement_evidence {
        let reference_specs: Vec<(String, PathBuf, PathBuf)> = options
            .reference_pairs
            .iter()
            .filter_map(|spec| {
                let mut parts = spec.splitn(3, ':');
                let label = parts.next()?.to_string();
                let path_a = PathBuf::from(parts.next()?);
                let path_b = PathBuf::from(parts.next()?);
                Some((label, path_a, path_b))
            })
            .collect();
        ensure(
            options.reference_pairs.len() == reference_specs.len(),
            "--reference-pair needs LABEL:A.json:B.json"
        )?;
        let dp_summary = run_routed_dp(
            &panel,
            &sources,
            &graph,
            &options.routes,
            &path_of_source,
            &component_axis[locus_lo..locus_hi],
            traversals[locus_lo..locus_hi].to_vec(),
            &bed_territory,
            k,
            &routed_equal,
            &universe.component_locus_to_partition,
            locus_lo,
            traversals.len(),
            target.id,
            target.length,
            &model,
            &sample_index.counts,
            evidence_path,
            options.refinement_ranks.as_deref(),
            options.refinement_margin,
            options.refinement_top_k,
            options.beam_width,
            options.tie_rotations,
            options.finalists,
            &reference_specs,
            &site_observed,
            &window_obs,
            std::path::Path::new("dp-rescore"),
            &mut rss,
        )?;
        if options.dp_only {
            println!(
                "{}",
                serde_json::to_string_pretty(&serde_json::json!({
                    "model": "experimental-mem-routed-evidence-v1",
                    "component": options.component,
                    "locus_range": [locus_lo, locus_hi],
                    "depth": options.depth,
                    "background": options.background,
                    "routing_universe": options.routing_universe,
                    "rss_peak_bytes": rss.peak_bytes(),
                    "stage_wall_seconds": {
                        "input_load": input_seconds,
                        "candidate_domain": domain_seconds,
                        "sample_record_rederivation": records_seconds,
                        "territory_index": territory_seconds,
                        "route_search": route_search_seconds,
                        "share_accumulation": accumulation_seconds,
                        "total": started.elapsed().as_secs_f64(),
                    },
                    "dp": dp_summary,
                    "routing": routing_summary,
                }))?
            );
            return Ok(());
        }
        dp_summary_global = Some(dp_summary);
    }

    // ------------------------------------------------- geometric profiles
    let geometric_started = Instant::now();
    let geometric: Vec<(Profile, usize)> = forward_candidates
        .par_iter()
        .map(|(_, segments, identity)| {
            let &(source, start, end) = &segments[0];
            let (profile, _nodes, ir_count) = geometric_candidate_profile(
                &panel,
                &sources,
                source,
                start,
                end,
                path_of_source[source],
                k,
                identity,
            )?;
            Ok::<_, io::Error>((profile, ir_count))
        })
        .collect::<io::Result<_>>()?;
    let geometric_seconds = geometric_started.elapsed().as_secs_f64();
    let geometric_feature_total: u64 = geometric.iter().map(|(p, _)| p.len() as u64).sum();
    let geometric_ir_records_total: u64 =
        geometric.iter().map(|(_, records)| *records as u64).sum();

    // ------------------------------------------------- exit-cut test (locus 26)
    let exit_started = Instant::now();
    let exit_locus = 26usize;
    ensure(
        exit_locus >= locus_lo && exit_locus < locus_hi,
        "exit locus 26 outside slice",
    )?;
    let s288c = *source_of_name
        .get("S288C#0#chrIII")
        .ok_or_else(|| invalid("S288C#0#chrIII absent"))?;
    let sk1 = *source_of_name
        .get("SK1#0#chrIII")
        .ok_or_else(|| invalid("SK1#0#chrIII absent"))?;
    // Historical identities (chrIII-full-prefilter-ranks-v1, locus 26):
    // plateau ...@204179>...@224054 loss -34276.450882011006
    // truth    ...@207743>...@227721 loss -34150.53126098728
    let historical_plateau_loss = -34276.450882011006;
    let historical_truth_loss = -34150.53126098728;
    let native_segment = (s288c, 220062u64, 230049u64);
    let truth_exit: Vec<(usize, u64, u64)> = vec![(sk1, 200187, 207743), (s288c, 227721, 230049)];
    let plateau_exit: Vec<(usize, u64, u64)> =
        vec![(sk1, 200187, 204179), (s288c, 224054, 230049)];
    ensure(
        forward_candidates.iter().any(|(locus, segments, identity)| {
            *locus == exit_locus
                && segments.len() == 1
                && segments[0] == native_segment
                && identity.starts_with("single:260000212:")
        }),
        "locus-26 native allele missing from the slice domain",
    )?;

    let native_sequence = fetch_segment(&sources, native_segment.0, native_segment.1, native_segment.2)?;
    let native_oracle = oracle_segment_profile(&panel, &native_sequence)?;
    let (native_geometric, _native_nodes, _native_ir) = geometric_candidate_profile(
        &panel,
        &sources,
        s288c,
        native_segment.1,
        native_segment.2,
        path_of_source[s288c],
        k,
        "native26",
    )?;

    let truth_oracle = oracle_split_profile(&panel, &sources, &truth_exit)?;
    let plateau_oracle = oracle_split_profile(&panel, &sources, &plateau_exit)?;
    let seam_truth = {
        let left = fetch_segment(&sources, truth_exit[0].0, truth_exit[0].1, truth_exit[0].2)?;
        let right = fetch_segment(&sources, truth_exit[1].0, truth_exit[1].1, truth_exit[1].2)?;
        let (seam, _) = genome::profile_event_seam(&panel, &left, &right, READ_LENGTH, MAX_FEATURES)?;
        seam
    };
    let seam_plateau = {
        let left = fetch_segment(&sources, plateau_exit[0].0, plateau_exit[0].1, plateau_exit[0].2)?;
        let right = fetch_segment(&sources, plateau_exit[1].0, plateau_exit[1].1, plateau_exit[1].2)?;
        let (seam, _) =
            genome::profile_event_seam(&panel, &left, &right, READ_LENGTH, MAX_FEATURES)?;
        seam
    };
    let truth_geometric = {
        let (left, _nodes_l, _ir_l) = geometric_candidate_profile(
            &panel,
            &sources,
            truth_exit[0].0,
            truth_exit[0].1,
            truth_exit[0].2,
            path_of_source[truth_exit[0].0],
            k,
            "truth-left",
        )?;
        let (right, _nodes_r, _ir_r) = geometric_candidate_profile(
            &panel,
            &sources,
            truth_exit[1].0,
            truth_exit[1].1,
            truth_exit[1].2,
            path_of_source[truth_exit[1].0],
            k,
            "truth-right",
        )?;
        merge_profiles(&[&left, &right, &seam_truth])
    }?;
    let plateau_geometric = {
        let (left, _nodes_l, _ir_l) = geometric_candidate_profile(
            &panel,
            &sources,
            plateau_exit[0].0,
            plateau_exit[0].1,
            plateau_exit[0].2,
            path_of_source[plateau_exit[0].0],
            k,
            "plateau-left",
        )?;
        let (right, _nodes_r, _ir_r) = geometric_candidate_profile(
            &panel,
            &sources,
            plateau_exit[1].0,
            plateau_exit[1].1,
            plateau_exit[1].2,
            path_of_source[plateau_exit[1].0],
            k,
            "plateau-right",
        )?;
        merge_profiles(&[&left, &right, &seam_plateau])
    }?;

    let exit_partition = universe.component_locus_to_partition[exit_locus];
    let pooled = |feature: &FeatureKey| -> io::Result<u64> { sample_index.counts.count(feature) };
    let routed_equal_26 = routed_equal
        .get(&exit_partition)
        .cloned()
        .unwrap_or_default();
    let routed_occ_26 = routed_occurrence
        .get(&exit_partition)
        .cloned()
        .unwrap_or_default();
    let pooled_f = |feature: &FeatureKey| -> io::Result<f64> {
        Ok(sample_index.counts.count(feature)? as f64)
    };
    let routed_f = |feature: &FeatureKey| -> io::Result<f64> {
        Ok(routed_equal_26.get(feature).copied().unwrap_or(0.0))
    };
    let routed_occ_f = |feature: &FeatureKey| -> io::Result<f64> {
        Ok(routed_occ_26.get(feature).copied().unwrap_or(0.0))
    };

    let pooled_truth = composed_pair_score(&truth_oracle, &native_oracle, &pooled_f, &model)?;
    let pooled_plateau = composed_pair_score(&plateau_oracle, &native_oracle, &pooled_f, &model)?;
    let routed_truth = composed_pair_score(&truth_oracle, &native_oracle, &routed_f, &model)?;
    let routed_plateau = composed_pair_score(&plateau_oracle, &native_oracle, &routed_f, &model)?;
    let routed_occ_truth =
        composed_pair_score(&truth_oracle, &native_oracle, &routed_occ_f, &model)?;
    let routed_occ_plateau =
        composed_pair_score(&plateau_oracle, &native_oracle, &routed_occ_f, &model)?;
    let geo_pooled_truth =
        composed_pair_score(&truth_geometric, &native_geometric, &pooled_f, &model)?;
    let geo_pooled_plateau =
        composed_pair_score(&plateau_geometric, &native_geometric, &pooled_f, &model)?;
    let geo_routed_truth =
        composed_pair_score(&truth_geometric, &native_geometric, &routed_f, &model)?;
    let geo_routed_plateau =
        composed_pair_score(&plateau_geometric, &native_geometric, &routed_f, &model)?;

    let pooled_reproduces = (pooled_plateau - historical_plateau_loss).abs() < 1e-6
        && (pooled_truth - historical_truth_loss).abs() < 1e-6;
    let delta_rows = feature_deltas(
        &plateau_oracle,
        &truth_oracle,
        &native_oracle,
        &pooled,
        &|feature| routed_f(feature).unwrap_or(0.0),
        &model,
    )?;
    let top_shift: Vec<&FeatureRow> = delta_rows.iter().take(20).collect();
    let pooled_support_total: u64 = {
        let mut keys: Vec<&FeatureKey> = plateau_oracle
            .keys()
            .chain(truth_oracle.keys())
            .chain(native_oracle.keys())
            .collect();
        keys.sort_unstable();
        keys.dedup();
        let mut total = 0u64;
        for key in keys {
            total += pooled(key)?;
        }
        total
    };
    let routed_support_total: f64 = {
        let mut keys: Vec<&FeatureKey> = plateau_oracle
            .keys()
            .chain(truth_oracle.keys())
            .chain(native_oracle.keys())
            .collect();
        keys.sort_unstable();
        keys.dedup();
        keys.iter()
            .map(|key| routed_f(key).unwrap_or(0.0))
            .sum::<f64>()
    };
    let exit_seconds = exit_started.elapsed().as_secs_f64();

    // ------------------------------------------------- equivalence sample
    let equivalence_started = Instant::now();
    let mut equivalence_loci = Vec::new();
    for item in options.equivalence_loci.split(',') {
        let locus: usize = item
            .trim()
            .parse()
            .map_err(|_| invalid("equivalence locus"))?;
        ensure(
            locus >= locus_lo && locus < locus_hi,
            "equivalence locus outside slice",
        )?;
        equivalence_loci.push(locus);
    }
    let mut sample_candidates: Vec<(usize, Vec<(usize, u64, u64)>, String, bool)> = Vec::new();
    for locus in &equivalence_loci {
        let per_locus: Vec<_> = forward_candidates
            .iter()
            .filter(|(c_locus, _, _)| c_locus == locus)
            .take(options.equivalence_per_locus)
            .cloned()
            .collect();
        for (c_locus, segments, identity) in per_locus {
            sample_candidates.push((c_locus, segments, identity, false));
        }
    }
    for (locus, segments, identity) in reverse_sample.iter() {
        sample_candidates.push((*locus, segments.clone(), identity.clone(), true));
    }
    let mut exit_segments: Vec<(usize, Vec<(usize, u64, u64)>, String, bool)> = Vec::new();
    for (name, segment) in [
        ("native26", native_segment),
        ("truth-left", truth_exit[0]),
        ("truth-right", truth_exit[1]),
        ("plateau-left", plateau_exit[0]),
        ("plateau-right", plateau_exit[1]),
    ] {
        exit_segments.push((
            exit_locus,
            vec![segment],
            format!("exit-segment:{name}"),
            false,
        ));
    }
    sample_candidates.extend(exit_segments);

    let mut profile_diffs = Vec::new();
    let mut window_diags = Vec::new();
    let mut oracle_profile_count = 0usize;
    let mut geometric_exact_candidates = 0usize;
    for (index, (locus, segments, identity, reverse)) in sample_candidates.iter().enumerate() {
        let &(source, start, end) = &segments[0];
        let sequence = fetch_segment(&sources, source, start, end)?;
        let sequence = if *reverse {
            impg::graph::reverse_complement(&sequence)
        } else {
            sequence
        };
        let oracle = oracle_segment_profile(&panel, &sequence)?;
        oracle_profile_count += 1;
        // Mirrored geometry for reverse candidates: the RC candidate's
        // per-window record set equals the forward candidate's (both views
        // of the same sequence are extracted by `collect_tagged_read`), so
        // the forward interval's geometric profile is the reverse
        // candidate's profile.
        let (geometric, model_nodes, _ir) = geometric_candidate_profile(
            &panel,
            &sources,
            source,
            start,
            end,
            path_of_source[source],
            k,
            identity,
        )?;
        let diff = profile_equivalence(identity, segments, *reverse, &oracle, &geometric, &model_nodes);
        if diff.oracle_features == diff.exact_features {
            geometric_exact_candidates += 1;
        }
        profile_diffs.push(diff);
        if index < 4 && !*reverse && end - start >= READ_LENGTH as u64 {
            let anchors =
                contained_path_anchors(&panel, path_of_source[source], start, end, k)?;
            let len = end - start;
            let stride = ((len as usize - READ_LENGTH) / 16).max(1);
            for offset in (0..=(len as usize - READ_LENGTH)).step_by(stride).take(16) {
                let window = &sequence[offset..offset + READ_LENGTH];
                let oracle_records =
                    mem_records::canonical_mem_records(&panel, window)?;
                let lo = start + offset as u64;
                let hi = lo + READ_LENGTH as u64;
                let window_anchors: Vec<(i32, u64)> = anchors
                    .iter()
                    .filter(|&&(_, bp)| bp >= offset as u64 && bp + k <= offset as u64 + READ_LENGTH as u64)
                    .map(|&(node, bp)| (node, bp - offset as u64))
                    .collect();
                let geometric_records: Vec<Vec<u64>> = if window_anchors.is_empty() {
                    Vec::new()
                } else {
                    vec![canonical(&encode_walk(&window_anchors)?)]
                };
                let matched: Vec<(i32, u64)> = panel
                    .matched_syncmers_in_sequence(window)
                    .iter()
                    .map(|m| (m.signed_node, m.query_pos))
                    .collect();
                window_diags.push(WindowDiag {
                    identity: identity.clone(),
                    window_start: lo,
                    window_len: READ_LENGTH,
                    oracle_records,
                    geometric_records,
                    matched_anchor_nodes: matched.iter().map(|&(node, _)| node).collect(),
                    path_anchor_nodes: window_anchors.iter().map(|&(node, _)| node).collect(),
                    matched_anchor_positions: matched.iter().map(|&(_, pos)| pos).collect(),
                    path_anchor_positions: window_anchors
                        .iter()
                        .map(|&(_, bp)| bp)
                        .collect(),
                    oracle_matches_geometric: false,
                });
                let last = window_diags.last_mut().expect("just pushed");
                let mut oracle_sorted = last.oracle_records.clone();
                oracle_sorted.sort();
                let mut geometric_sorted = last.geometric_records.clone();
                geometric_sorted.sort();
                last.oracle_matches_geometric = oracle_sorted == geometric_sorted;
            }
        }
    }
    let equivalence_seconds = equivalence_started.elapsed().as_secs_f64();

    // ------------------------------------------------- report
    let new_path_seconds = input_seconds
        + domain_seconds
        + records_seconds
        + territory_seconds
        + geometric_seconds
        + (route_search_seconds + accumulation_seconds);
    let exit_summary = serde_json::json!({
        "historical_plateau_loss": historical_plateau_loss,
        "historical_truth_loss": historical_truth_loss,
        "pooled_reproduces_historical": pooled_reproduces,
        "pooled": {
            "plateau_loss": pooled_plateau,
            "truth_loss": pooled_truth,
            "plateau_minus_truth": pooled_plateau - pooled_truth,
        },
        "routed_equal_share": {
            "plateau_loss": routed_plateau,
            "truth_loss": routed_truth,
            "plateau_minus_truth": routed_plateau - routed_truth,
        },
        "routed_occurrence_weighted": {
            "plateau_loss": routed_occ_plateau,
            "truth_loss": routed_occ_truth,
            "plateau_minus_truth": routed_occ_plateau - routed_occ_truth,
        },
        "geometric_predicted_pooled": {
            "plateau_loss": geo_pooled_plateau,
            "truth_loss": geo_pooled_truth,
            "plateau_minus_truth": geo_pooled_plateau - geo_pooled_truth,
        },
        "geometric_predicted_routed_equal_share": {
            "plateau_loss": geo_routed_plateau,
            "truth_loss": geo_routed_truth,
            "plateau_minus_truth": geo_routed_plateau - geo_routed_truth,
        },
        "composed_feature_pooled_support_total": pooled_support_total,
        "composed_feature_routed_support_total": routed_support_total,
        "top_feature_shifts": top_shift,
        "inversion_fixed": (pooled_plateau - pooled_truth) > 0.0
            && (routed_plateau - routed_truth) <= 0.0,
    });
    let equivalence_summary = serde_json::json!({
        "sample_candidates": profile_diffs.len(),
        "oracle_profiles_computed": oracle_profile_count,
        "geometric_exact_candidates": geometric_exact_candidates,
        "profile_diffs": profile_diffs,
        "window_diagnostics": window_diags,
    });
    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::json!({
            "model": "experimental-mem-routed-evidence-v1",
            "component": options.component,
            "locus_range": [locus_lo, locus_hi],
            "depth": options.depth,
            "background": options.background,
            "routing_universe": options.routing_universe,
            "rss_peak_bytes": rss.peak_bytes(),
            "slice_loci": slice_loci,
            "forward_candidates": forward_candidates.len(),
            "geometric_feature_total": geometric_feature_total,
            "geometric_ir_extractions": forward_candidates.len() as u64,
            "geometric_ir_records": geometric_ir_records_total,
            "sample_stats": {
                "reads": sample_index.stats.reads,
                "mem_records": sample_index.stats.mem_records,
                "distinct_mems": sample_index.stats.distinct_mems,
                "spot_checked_against_index": spot_checked,
            },
            "routing": routing_summary,
            "dp": dp_summary_global,
            "exit_cut_test": exit_summary,
            "equivalence": equivalence_summary,
            "stage_wall_seconds": {
                "input_load": input_seconds,
                "candidate_domain": domain_seconds,
                "sample_record_rederivation": records_seconds,
                "territory_index": territory_seconds,
                "route_search": route_search_seconds,
                "share_accumulation": accumulation_seconds,
                "geometric_profiles": geometric_seconds,
                "exit_oracle_scoring": exit_seconds,
                "equivalence_sample": equivalence_seconds,
                "new_path_total": new_path_seconds,
                "total": started.elapsed().as_secs_f64(),
            },
        }))?
    );
    Ok(())
}

#[cfg(test)]
mod multiplicity_background_tests {
    use super::*;

    /// The genome-wide cross-support contribution (the owner-decided
    /// candidate family 1): a record with no placement or a single
    /// placement contributes exactly 0 (the latter's whole multiplicity
    /// floods one partition and is fully observed there); otherwise the
    /// record contributes its multiplicity routed outside each one of its
    /// own placements, m_r*(1 - 1/t_r_genome).
    #[test]
    fn cross_support_contribution_formula() {
        assert_eq!(cross_support_contribution(200, 0), 0.0);
        assert_eq!(cross_support_contribution(200, 1), 0.0);
        assert!((cross_support_contribution(4, 2) - 2.0).abs() < 1e-12);
        assert!((cross_support_contribution(9, 3) - 6.0).abs() < 1e-12);
        // The locus-21 flooding shape: multiplicity ~200 spread over ~200
        // genome placements leaves ~199 routed outside each one.
        assert!((cross_support_contribution(200, 200) - 199.0).abs() < 1e-9);
    }

    /// The haploid sample-side background is the feature-level genome-wide
    /// cross support plus the base — the same value at every charge site
    /// (locus and boundary); features no record realizes keep the base.
    #[test]
    fn sample_side_background_beta_genome_wide() {
        let mut cross: HashMap<FeatureKey, f64> = HashMap::new();
        cross.insert(vec![7u64], 58.7);
        let sample = SampleSideBackgrounds::new(&cross, 0.1);
        assert!((sample.beta(&vec![7u64]) - 58.8).abs() < 1e-9);
        assert!((sample.beta(&vec![9u64]) - 0.1).abs() < 1e-12);
    }

    fn model(background: f64) -> ScoreModel {
        ScoreModel {
            read_length: 150,
            histogram: 100,
            denominator: 15_000.0,
            depth: 10.0,
            background,
        }
    }

    fn entry(paths: u32) -> FeatureBackgroundEntry {
        let per_count = 100.0 * 10.0 / 15_000.0;
        FeatureBackgroundEntry {
            paths,
            incidence: 0,
            weight: multiplicity_attribution_weight(paths),
            beta: multiplicity_background_beta(0.1, paths, per_count),
        }
    }

    /// The attribution weight (owner decision (d)): unique and panel-absent
    /// features keep full credit; a feature in m contexts is attributable to
    /// this context with probability 1/m.
    #[test]
    fn attribution_weight_is_one_over_multiplicity() {
        assert_eq!(multiplicity_attribution_weight(0), 1.0);
        assert_eq!(multiplicity_attribution_weight(1), 1.0);
        assert!((multiplicity_attribution_weight(232) - 1.0 / 232.0).abs() < 1e-12);
        assert!((multiplicity_attribution_weight(3375) - 1.0 / 3375.0).abs() < 1e-12);
    }

    /// The (b) diagnostic background: beta_f = base + (m-1)*per_count.
    #[test]
    fn gentle_beta_formula() {
        let base = 0.1;
        let per_count = 100.0 * 10.0 / 15_000.0;
        assert_eq!(multiplicity_background_beta(base, 0, per_count), base);
        assert_eq!(multiplicity_background_beta(base, 1, per_count), base);
        let beta = multiplicity_background_beta(base, 232, per_count);
        let expected = base + 231.0 * per_count;
        assert!((beta - expected).abs() < 1e-9, "{beta} vs {expected}");
    }

    /// The attribution loss keeps unique features bit-identical to the flat
    /// loss and scales the credit term of pooled features by 1/m.
    #[test]
    fn attribution_loss_scales_only_pooled_credit() {
        let m = model(0.1);
        for q in [0u64, 1, 5, 86, 1000] {
            for observed in [0.0f64, 0.5, 3.25, 137.0] {
                let flat = loss_fractional(&m, q, observed).unwrap();
                // m = 0 (never scanned) and m = 1 (unique): bit-identical.
                assert_eq!(
                    flat.to_bits(),
                    loss_fractional_entry(&m, q, observed, None).unwrap().to_bits()
                );
                assert_eq!(
                    flat.to_bits(),
                    loss_fractional_entry(&m, q, observed, Some(entry(1)))
                        .unwrap()
                        .to_bits()
                );
                // A pooled feature earns strictly less credit (higher loss)
                // for positive observed support, and its signal term is
                // unchanged.
                let pooled = loss_fractional_entry(&m, q, observed, Some(entry(232)))
                    .unwrap();
                if observed > 0.0 && q > 0 {
                    // The credit term is scaled by 1/232: strictly less
                    // credit, same signal.
                    assert!(pooled > flat);
                } else {
                    assert_eq!(pooled.to_bits(), flat.to_bits());
                }
            }
        }
    }

    /// Window containment arithmetic: a subwalk spanning [p_first, p_last]
    /// with k-mers fits the L150 windows starting in
    /// [p_last + k - 150, p_first] clipped to [0, len - 150].
    #[test]
    fn window_incidence_matches_containment_arithmetic() {
        let k = 63u64;
        let len = 10_000u64;
        // Interior occurrence spanning 50bp: 151 - k - (p_last - p_first).
        assert_eq!(window_incidence(300, 350, k, len), 38);
        // A single anchor (p_first == p_last): 151 - k = 88 windows.
        assert_eq!(window_incidence(300, 300, k, len), 88);
        // Clipped at the path start.
        assert_eq!(window_incidence(10, 30, k, len), 11);
        // Clipped at the path end: an occurrence reaching the final 5bp of
        // the path cannot fit any window (its last anchor needs a window
        // start >= p_last + k - 150, beyond the final start len - 150).
        assert_eq!(window_incidence(9_990, 9_995, k, len), 0);
        // An occurrence whose final anchor sits exactly at the last
        // window's containment edge: exactly one containing window start.
        assert_eq!(window_incidence(9_900, 9_937, k, len), 1);
        // An occurrence wider than a read window admits none.
        assert_eq!(window_incidence(100, 300, k, len), 0);
        // A path shorter than a read has no windows at all.
        assert_eq!(window_incidence(0, 0, k, 100), 0);
    }

    /// Fix 1 (record-once mixed-owner charging): a single-owner site map is
    /// EXACTLY the per-partition routed-share accumulation (same records,
    /// same order, same values — bit-identical), and a mixed-owner site map
    /// attributes a record touching k charged owners ONCE — not k times as
    /// the union-sum form did.
    #[test]
    fn record_once_site_map_single_owner_bit_identity_and_cross_owner_dedup() {
        let walk = vec![(5i32, 0u64), (7, 10), (-3, 20)];
        let tokens = impg::sample_mem_bwt::encode_walk(&walk).unwrap();
        // Record A: multiplicity 10, touched {1, 2}. Record B: multiplicity
        // 4, touched {2} only.
        let records = vec![
            RoutedRecord {
                tokens: tokens.clone(),
                multiplicity: 10,
                occurrences: BTreeMap::from([(1u32, 1u64), (2u32, 1u64)]),
                total_occurrences: 2,
                forward_positions: Vec::new(),
                reverse_positions: Vec::new(),
            },
            RoutedRecord {
                tokens: tokens.clone(),
                multiplicity: 4,
                occurrences: BTreeMap::from([(2u32, 1u64)]),
                total_occurrences: 1,
                forward_positions: Vec::new(),
                reverse_positions: Vec::new(),
            },
        ];
        let slice: BTreeSet<u32> = BTreeSet::from([1u32, 2u32]);
        let site = SiteObserved::build(&records, &slice);
        let features = enumerate_subwalks(&tokens).unwrap();
        assert!(!features.is_empty());
        // Single-owner maps: exactly the accumulation form's values.
        let map1 = site.site_map([1u32]);
        let map2 = site.site_map([2u32]);
        for feature in &features {
            assert_eq!(map1[feature], 5.0f64); // A: 10 / 2
            assert_eq!(map2[feature], 9.0f64); // A: 5 + B: 4
        }
        // The mixed-owner site: each record ONCE (A: 5, B: 4) — the
        // union-sum form would have charged A twice (10) here.
        let mixed = site.site_map([1u32, 2u32]);
        for feature in &features {
            assert_eq!(mixed[feature], 9.0f64);
        }
        // Determinism: the same owner set rebuilds the same values.
        let mixed_again = site.site_map([2u32, 1u32]);
        for feature in &features {
            assert_eq!(mixed[feature], mixed_again[feature]);
        }
    }

    /// Fix 2 (the omission charge): an omitted observed feature contributes
    /// exactly beta_f - C_f*ln(beta_f); an explained feature reduces to the
    /// current relative form plus that same per-feature constant; the
    /// omittor-vs-explainer difference equals the current discriminative
    /// term; and two rows omitting the same feature keep their differences
    /// unchanged.
    #[test]
    fn omission_charge_unbaselined_term_and_preserved_discriminations() {
        let model = model(0.1);
        let mut cross: HashMap<FeatureKey, f64> = HashMap::new();
        cross.insert(vec![7u64], 0.0);
        cross.insert(vec![9u64], 0.0);
        let sample = SampleSideBackgrounds::new(&cross, 0.1);
        let mut obs: HashMap<FeatureKey, f64> = HashMap::new();
        obs.insert(vec![7u64], 5.0);
        // The pure omittor: beta_f - C_f*ln(beta_f).
        let omittor = merged_single_loss_sample_with_omission(
            &Profile::new(),
            &obs,
            &obs,
            &sample,
            &model,
        )
        .unwrap();
        let expected = 0.1 - 5.0 * 0.1f64.ln() + ln_gamma_observation(5.0);
        assert!((omittor - expected).abs() < 1e-9);
        // The explainer: its current charge BIT-EXACTLY (nothing omitted).
        let mut profile = Profile::new();
        profile.insert(vec![7u64], 2);
        let explainer = merged_single_loss_sample_with_omission(
            &profile,
            &obs,
            &obs,
            &sample,
            &model,
        )
        .unwrap();
        let old_explainer =
            merged_single_loss_sample(&profile, &obs, &sample, &model).unwrap();
        assert_eq!(explainer, old_explainer);
        // omittor - explainer == the current discriminative term PLUS the
        // feature's omission cost (the new discrimination).
        assert!(((omittor - explainer) + old_explainer - expected).abs() < 1e-9);
        // Two rows omitting f1 but differing only on an UNOBSERVED feature:
        // the common omission term cancels, the signal difference is
        // preserved bit-exactly.
        let mut profile_g = Profile::new();
        profile_g.insert(vec![9u64], 1);
        let row_a = merged_single_loss_sample_with_omission(
            &profile_g,
            &obs,
            &obs,
            &sample,
            &model,
        )
        .unwrap();
        let row_b = merged_single_loss_sample_with_omission(
            &Profile::new(),
            &obs,
            &obs,
            &sample,
            &model,
        )
        .unwrap();
        let s_g = 1.0 * (model.histogram as f64) * model.depth / model.denominator;
        assert!(((row_a - row_b) - s_g).abs() < 1e-9);
    }

    /// Ruling 2 (the chain-level once-per-material omission,
    /// genome/stitching-omission-alignment): a window's omission addend
    /// covers the window's observed features the CHAIN explains NOWHERE —
    /// a feature in the chain's explained set (spelled at any window) is
    /// exonerated here; an EMPTY explained set reproduces the per-window
    /// form bit-exactly (the DP's surrogate layer); the row's own charge is
    /// untouched in both forms.
    #[test]
    fn chain_level_omission_exonerates_chain_explained_features() {
        let model = model(0.1);
        let mut cross: HashMap<FeatureKey, f64> = HashMap::new();
        cross.insert(vec![7u64], 0.0);
        cross.insert(vec![9u64], 0.0);
        cross.insert(vec![11u64], 0.0);
        let sample = SampleSideBackgrounds::new(&cross, 0.1);
        let mut window_obs: HashMap<FeatureKey, f64> = HashMap::new();
        window_obs.insert(vec![7u64], 5.0);
        window_obs.insert(vec![9u64], 5.0);
        // The shared-once-emptied row shape: the chain's charged profile at
        // this window is EMPTY (its material was charged at its first
        // window), yet the chain explains feature 7 somewhere.
        let profile = Profile::new();
        let per_window = omission_charge_addend(&profile, &window_obs, &sample);
        // The per-window form charges BOTH observed features.
        let only_7 = {
            let mut map: HashMap<FeatureKey, f64> = HashMap::new();
            map.insert(vec![7u64], 5.0);
            omission_charge_addend(&profile, &map, &sample)
        };
        let only_9 = {
            let mut map: HashMap<FeatureKey, f64> = HashMap::new();
            map.insert(vec![9u64], 5.0);
            omission_charge_addend(&profile, &map, &sample)
        };
        // Empty explained set: bit-identical to the per-window form.
        let empty: std::collections::HashSet<FeatureKey> = std::collections::HashSet::new();
        assert_eq!(
            omission_charge_addend_except(&profile, &empty, &window_obs, &sample).to_bits(),
            per_window.to_bits()
        );
        // The chain explains feature 7 at its first window: exonerated
        // here — exactly feature 9's omission term remains.
        let explained: std::collections::HashSet<FeatureKey> =
            [vec![7u64]].into_iter().collect();
        assert_eq!(
            omission_charge_addend_except(&profile, &explained, &window_obs, &sample).to_bits(),
            only_9.to_bits()
        );
        assert!((per_window - (only_7 + only_9)).abs() < 1e-9);
        // An explained feature the window does not observe changes nothing.
        let explained_elsewhere: std::collections::HashSet<FeatureKey> =
            [vec![11u64]].into_iter().collect();
        assert_eq!(
            omission_charge_addend_except(&profile, &explained_elsewhere, &window_obs, &sample)
                .to_bits(),
            per_window.to_bits()
        );
        // The row's own charge is untouched: with a NON-empty profile the
        // The row's own charge under the chain-level form (below).
        let mut predicted: HashMap<FeatureKey, f64> = HashMap::new();
        predicted.insert(vec![7u64], 5.0);
        let mut row = Profile::new();
        row.insert(vec![7u64], 2);
        let plain = merged_single_loss_sample_with_omission(
            &row,
            &predicted,
            &window_obs,
            &sample,
            &model,
        )
        .unwrap();
        let row_charge = merged_single_loss_sample(&row, &predicted, &sample, &model).unwrap();
        // With BOTH window features in the chain's explained set (7 predicted
        // here, 9 spelled at another window) the chain-level form is exactly
        // the row's own charge; the per-window form charges exactly feature
        // 9's omission term on top.
        let explained_both: std::collections::HashSet<FeatureKey> =
            [vec![7u64], vec![9u64]].into_iter().collect();
        let chain_level = merged_single_loss_sample_with_omission_except(
            &row,
            &predicted,
            &explained_both,
            &window_obs,
            &sample,
            &model,
        )
        .unwrap();

        assert_eq!(chain_level.to_bits(), row_charge.to_bits());
        assert_eq!(plain.to_bits(), (row_charge + only_9).to_bits());
    }

    /// The finalist objective's admissibility shape (supervisor ruling,
    /// genome/finalist-reranking): the model-of-record per-window charge (the
    /// chain-level except-form) is AT LEAST the allele's own base charge —
    /// the surrogate loss minus the allele's own in-window omission addend —
    /// because the omission channel's refund charges only window-observed
    /// features the allele does not explain in-window. The L objective
    /// therefore lower-bounds the model-of-record total up to the measured
    /// profile-form residual, which the run reports per enumerated finalist.
    #[test]
    fn finalist_l_form_lower_bounds_chain_level_charge() {
        let model = model(0.1);
        let mut cross: HashMap<FeatureKey, f64> = HashMap::new();
        cross.insert(vec![7u64], 0.0);
        cross.insert(vec![9u64], 0.0);
        let sample = SampleSideBackgrounds::new(&cross, 0.1);
        let mut window_obs: HashMap<FeatureKey, f64> = HashMap::new();
        window_obs.insert(vec![7u64], 5.0);
        window_obs.insert(vec![9u64], 5.0);
        let mut predicted: HashMap<FeatureKey, f64> = HashMap::new();
        predicted.insert(vec![7u64], 5.0);
        let mut row = Profile::new();
        row.insert(vec![7u64], 2);
        // The allele's own base charge (the L form's per-window value) and
        // its surrogate loss (base + the in-window omission addend).
        let base = merged_single_loss_sample(&row, &predicted, &sample, &model).unwrap();
        let surrogate_loss = merged_single_loss_sample_with_omission(
            &row,
            &predicted,
            &window_obs,
            &sample,
            &model,
        )
        .unwrap();
        let addend = omission_charge_addend(&row, &window_obs, &sample);
        assert_eq!(
            surrogate_loss.to_bits(),
            (base + addend).to_bits(),
            "the surrogate loss is exactly base + the allele's own addend"
        );
        assert!(addend >= 0.0, "the omission terms are nonnegative");
        // The model-of-record charge with the chain explaining EVERYTHING the
        // window observes beyond the allele (feature 9 spelled elsewhere):
        // exactly the base charge — the addend was fully refundable.
        let explained: std::collections::HashSet<FeatureKey> =
            [vec![7u64], vec![9u64]].into_iter().collect();
        let model_charge = merged_single_loss_sample_with_omission_except(
            &row,
            &predicted,
            &explained,
            &window_obs,
            &sample,
            &model,
        )
        .unwrap();
        assert_eq!(model_charge.to_bits(), base.to_bits());
        // With a PARTIAL refund (the chain spells nothing beyond the allele):
        // the full per-window form — still at least the base charge.
        let empty: std::collections::HashSet<FeatureKey> = std::collections::HashSet::new();
        let unexonerated = merged_single_loss_sample_with_omission_except(
            &row,
            &predicted,
            &empty,
            &window_obs,
            &sample,
            &model,
        )
        .unwrap();
        assert_eq!(unexonerated.to_bits(), surrogate_loss.to_bits());
        assert!(unexonerated >= base);
    }

    /// The CORRECTED Ruling 2 under S2 (record-level placement coverage,
    /// supervisor-approved on the measured placement structure): a record's
    /// share of a feature is covered iff the chain spells the feature
    /// overlapping one of THE RECORD'S OWN placement spans. The
    /// repeat-at-distinct-loci case (different reads): a record that never
    /// places at the spelled copy stays charged — the unspelled instance's
    /// omission RETURNS at its own mass. The flood-record case (one read
    /// placing at several loci): the one spelling at the read's own
    /// placement covers its share at EVERY window it touches (one read =
    /// one observation). MASS CONSERVATION asserted throughout.
    #[test]
    fn instance_exoneration_charges_distinct_loci_instances() {
        let model = model(0.1);
        let mut cross: HashMap<FeatureKey, f64> = HashMap::new();
        cross.insert(vec![7u64], 0.0);
        let sample = SampleSideBackgrounds::new(&cross, 0.1);
        // Two windows on DISTINCT paths (distinct loci), each observing the
        // same walk's features at its own DNA position.
        let index = TerritoryIndex {
            territories: Vec::new(),
            entries: Vec::new(),
            node_ranges: Vec::new(),
            path_intervals: vec![
                vec![(100u64, 200u64, 1u32)],
                vec![(100u64, 200u64, 2u32)],
            ],
        };
        let walk = vec![(5i32, 0u64), (7, 10), (-3, 20)];
        let tokens = impg::sample_mem_bwt::encode_walk(&walk).unwrap();
        let features = enumerate_subwalks(&tokens).unwrap();
        // Record A: one copy on path 0 (owner partition 1); record B: one
        // copy on path 1 (owner partition 2). SAME feature, DIFFERENT reads.
        let record_a = RoutedRecord {
            tokens: tokens.clone(),
            multiplicity: 10,
            occurrences: BTreeMap::from([(1u32, 1u64)]),
            total_occurrences: 1,
            forward_positions: vec![(0usize, 90u64)],
            reverse_positions: Vec::new(),
        };
        let record_b = RoutedRecord {
            tokens,
            multiplicity: 8,
            occurrences: BTreeMap::from([(2u32, 1u64)]),
            total_occurrences: 1,
            forward_positions: vec![(1usize, 90u64)],
            reverse_positions: Vec::new(),
        };
        let owner_sets = vec![BTreeSet::from([1u32]), BTreeSet::from([2u32])];
        let mut instance_spans: RecordPlacementSpans = Vec::new();
        let mut window_records: WindowRecordLists = Vec::new();
        let maps = build_in_window_obs(
            &[record_a, record_b],
            &index,
            63,
            &owner_sets,
            Some(&mut instance_spans),
            Some(&mut window_records),
        )
        .unwrap();
        let record_shares = vec![10.0f64, 8.0];
        let instances = InstanceStructure {
            record_spans: instance_spans,
            record_shares,
            window_records,
        };
        // MASS CONSERVATION: the contributing records' shares sum to the
        // observed map's mass, exactly.
        for (locus, map) in maps.iter().enumerate() {
            for (feature, mass) in map {
                let sum: f64 = instances.window_records[locus][feature]
                    .iter()
                    .map(|&r| instances.record_shares[r as usize])
                    .sum();
                assert_eq!(sum.to_bits(), mass.to_bits());
                assert!(!instances.window_records[locus][feature].is_empty());
            }
        }
        let profile = Profile::new();
        let feature = &features[features.len() - 1];
        // The chain spells the feature at record A's position ONLY (path 0,
        // the occurrence's anchor extent).
        let mut spelled: HashMap<FeatureKey, Vec<SpelledSpan>> = HashMap::new();
        spelled.insert(
            feature.clone(),
            vec![SpelledSpan {
                path: 0,
                lo: 90,
                hi: 90 + 20 + 63,
            }],
        );
        // Window 0 (record A's copy): A's own placements include the
        // spelled span — covered, no omission for it.
        let covered_window = merged_single_loss_sample_with_omission_instances(
            &profile,
            &maps[0],
            &maps[0],
            Some(&instances.window_records[0]),
            Some(&instances),
            &spelled,
            &sample,
            &model,
        )
        .unwrap();
        let explained_a: std::collections::HashSet<FeatureKey> =
            [feature.clone()].into_iter().collect();
        let old_covered = merged_single_loss_sample_with_omission_except(
            &profile,
            &maps[0],
            &explained_a,
            &maps[0],
            &sample,
            &model,
        )
        .unwrap();
        assert_eq!(
            covered_window.to_bits(),
            old_covered.to_bits(),
            "the covered record's omission is wiped exactly as the old form's"
        );
        // Window 1 (record B's copy at a DIFFERENT DNA position): B never
        // places at the spelled copy — the omission RETURNS at B's own mass.
        let other_window = merged_single_loss_sample_with_omission_instances(
            &profile,
            &maps[1],
            &maps[1],
            Some(&instances.window_records[1]),
            Some(&instances),
            &spelled,
            &sample,
            &model,
        )
        .unwrap();
        let share_b = 8.0f64;
        let base_b = merged_single_loss_sample(&profile, &maps[1], &sample, &model).unwrap();
        let mut keys: Vec<&FeatureKey> = maps[1].keys().collect();
        keys.sort();
        let mut expected = base_b;
        for f in keys {
            let mass = maps[1][f];
            let uncovered = if *f == *feature { share_b } else { mass };
            expected += omission_term_value(f, uncovered, &sample);
        }
        assert_eq!(
            other_window.to_bits(),
            expected.to_bits(),
            "the unspelled record's omission returns at its own mass"
        );
        // The OLD form (feature-level exoneration) wiped it entirely.
        let old = merged_single_loss_sample_with_omission_except(
            &profile,
            &maps[1],
            &explained_a,
            &maps[1],
            &sample,
            &model,
        )
        .unwrap();
        assert!(
            other_window > old,
            "the corrected form charges what the feature-level form wiped"
        );
        // THE FLOOD-RECORD CASE: one read placing at BOTH loci; the chain's
        // one spelling (at the read's own placement on path 0) covers its
        // share at EVERY window it touches — one read = one observation.
        let record_flood = RoutedRecord {
            tokens: impg::sample_mem_bwt::encode_walk(&walk).unwrap(),
            multiplicity: 12,
            occurrences: BTreeMap::from([(1u32, 1u64), (2u32, 1u64)]),
            total_occurrences: 2,
            forward_positions: vec![(0usize, 90u64), (1usize, 90u64)],
            reverse_positions: Vec::new(),
        };
        let mut instance_spans2: RecordPlacementSpans = Vec::new();
        let mut window_records2: WindowRecordLists = Vec::new();
        let maps2 = build_in_window_obs(
            &[record_flood],
            &index,
            63,
            &owner_sets,
            Some(&mut instance_spans2),
            Some(&mut window_records2),
        )
        .unwrap();
        let instances2 = InstanceStructure {
            record_spans: instance_spans2,
            record_shares: vec![6.0f64],
            window_records: window_records2,
        };
        for locus in 0..2 {
            let flood_window = merged_single_loss_sample_with_omission_instances(
                &profile,
                &maps2[locus],
                &maps2[locus],
                Some(&instances2.window_records[locus]),
                Some(&instances2),
                &spelled,
                &sample,
                &model,
            )
            .unwrap();
            let old_flood = merged_single_loss_sample_with_omission_except(
                &profile,
                &maps2[locus],
                &explained_a,
                &maps2[locus],
                &sample,
                &model,
            )
            .unwrap();
            assert_eq!(
                flood_window.to_bits(),
                old_flood.to_bits(),
                "the flood record's share is covered at every window by the one spelling at its own placement"
            );
        }
    }

    /// The group-shared ONE-INSTANCE case: the same DNA piece serving
    /// several overlapping windows observes the feature at the SAME
    /// coordinates in every window — one spelling covers them all (the
    /// records' own placements include the spelled piece), and the
    /// corrected form charges the omission exactly as the established
    /// once-per-material form does.
    #[test]
    fn instance_exoneration_group_shared_one_instance_stays_once() {
        let model = model(0.1);
        let mut cross: HashMap<FeatureKey, f64> = HashMap::new();
        cross.insert(vec![7u64], 0.0);
        let sample = SampleSideBackgrounds::new(&cross, 0.1);
        // TWO windows, SAME owner partition, SAME territory interval on the
        // same path: the piece's DNA is observed by both windows' maps at
        // the same coordinates.
        let index = TerritoryIndex {
            territories: Vec::new(),
            entries: Vec::new(),
            node_ranges: Vec::new(),
            path_intervals: vec![vec![(100u64, 200u64, 1u32)]],
        };
        let walk = vec![(5i32, 0u64), (7, 10), (-3, 20)];
        let tokens = impg::sample_mem_bwt::encode_walk(&walk).unwrap();
        let features = enumerate_subwalks(&tokens).unwrap();
        let record = RoutedRecord {
            tokens,
            multiplicity: 10,
            occurrences: BTreeMap::from([(1u32, 1u64)]),
            total_occurrences: 1,
            forward_positions: vec![(0usize, 90u64)],
            reverse_positions: Vec::new(),
        };
        let owner_sets = vec![BTreeSet::from([1u32]), BTreeSet::from([1u32])];
        let mut instance_spans: RecordPlacementSpans = Vec::new();
        let mut window_records: WindowRecordLists = Vec::new();
        let maps = build_in_window_obs(
            &[record],
            &index,
            63,
            &owner_sets,
            Some(&mut instance_spans),
            Some(&mut window_records),
        )
        .unwrap();
        let instances = InstanceStructure {
            record_spans: instance_spans,
            record_shares: vec![10.0f64],
            window_records,
        };
        let profile = Profile::new();
        let feature = &features[features.len() - 1];
        let mut spelled: HashMap<FeatureKey, Vec<SpelledSpan>> = HashMap::new();
        spelled.insert(
            feature.clone(),
            vec![SpelledSpan {
                path: 0,
                lo: 90,
                hi: 90 + 20 + 63,
            }],
        );
        for locus in 0..2 {
            let corrected = merged_single_loss_sample_with_omission_instances(
                &profile,
                &maps[locus],
                &maps[locus],
                Some(&instances.window_records[locus]),
                Some(&instances),
                &spelled,
                &sample,
                &model,
            )
            .unwrap();
            let explained: std::collections::HashSet<FeatureKey> =
                [feature.clone()].into_iter().collect();
            let old = merged_single_loss_sample_with_omission_except(
                &profile,
                &maps[locus],
                &explained,
                &maps[locus],
                &sample,
                &model,
            )
            .unwrap();
            assert_eq!(
                corrected.to_bits(),
                old.to_bits(),
                "the one-instance case keeps the established once-per-material charge"
            );
        }
    }

    /// The reverse piece's rc-span -> source-coordinate mapping: the rc
    /// span [o_lo, o_hi) maps to [end - o_hi, end - o_lo) (rc(S)[o..o+k] =
    /// rc(S[L-o-k..L-o])), and the forward mapping is the plain offset.
    #[test]
    fn piece_span_coords_maps_reverse_pieces() {
        let piece = (0usize, 1000u64, 2000u64, true);
        assert_eq!(piece_span_coords(piece, 100, 163), (2000 - 163, 2000 - 100));
        let forward = (0usize, 1000u64, 2000u64, false);
        assert_eq!(piece_span_coords(forward, 100, 163), (1100, 1163));
    }

    /// Fix 2's universe (owner ruling 2026-09-24, reading (B)): the
    /// in-window observed profile attributes a feature's share only for
    /// subwalk instances whose anchors lie inside the window's territory —
    /// a record whose verified occurrence is elsewhere contributes nothing,
    /// even though the established record-touching maps charge it.
    #[test]
    fn in_window_observed_attribution_is_positional() {
        // Window 0: owner partition 1, territory [100, 200) on path 0.
        // Window 1: owner partition 2, no territory on path 0.
        let index = TerritoryIndex {
            territories: Vec::new(),
            entries: Vec::new(),
            node_ranges: Vec::new(),
            path_intervals: vec![vec![(100u64, 200u64, 1u32)]],
        };
        let walk = vec![(5i32, 0u64), (7, 10), (-3, 20)];
        let tokens = impg::sample_mem_bwt::encode_walk(&walk).unwrap();
        let features = enumerate_subwalks(&tokens).unwrap();
        assert_eq!(features.len(), 6);
        // Record A: an occurrence at bp 90 on path 0 — its anchors (90, 100,
        // 110) all touch partition 1's [100, 200) with k = 63. The record
        // also touches partition 2 elsewhere (t = 2), so its equal share is
        // 10 / 2 = 5 per instance.
        let record_a = RoutedRecord {
            tokens: tokens.clone(),
            multiplicity: 10,
            occurrences: BTreeMap::from([(1u32, 1u64), (2u32, 1u64)]),
            total_occurrences: 1,
            forward_positions: vec![(0usize, 90u64)],
            reverse_positions: Vec::new(),
        };
        // Record B: its occurrence sits outside the window (bp 500); the
        // record-touching form would charge it (its walk's features), the
        // in-window form must not.
        let record_b = RoutedRecord {
            tokens,
            multiplicity: 4,
            occurrences: BTreeMap::from([(1u32, 1u64)]),
            total_occurrences: 1,
            forward_positions: vec![(0usize, 500u64)],
            reverse_positions: Vec::new(),
        };
        let owner_sets = vec![BTreeSet::from([1u32]), BTreeSet::from([2u32])];
        let maps =
            build_in_window_obs(&[record_a, record_b], &index, 63, &owner_sets, None, None).unwrap();
        assert_eq!(maps.len(), 2);
        for feature in &features {
            assert_eq!(maps[0][feature], 5.0f64); // A: 10 / 1 touched partition
        }
        assert!(maps[0].len() == features.len());
        assert!(maps[1].is_empty());
    }

    /// The in-window attribution counts each feature ONCE PER READ
    /// OCCURRENCE: a walk whose repeated region yields several anchor spans
    /// with the same canonical feature contributes its share once, not once
    /// per span (the (i, j) subwalk enumeration would multiply repeat
    /// features by their span count).
    #[test]
    fn in_window_observed_deduplicates_repeat_features_per_occurrence() {
        let index = TerritoryIndex {
            territories: Vec::new(),
            entries: Vec::new(),
            node_ranges: Vec::new(),
            path_intervals: vec![vec![(100u64, 300u64, 1u32)]],
        };
        // Node 5 recurs: the singleton subwalks (0,0) and (2,2) encode the
        // same canonical feature.
        let walk = vec![(5i32, 0u64), (9, 10), (5, 20)];
        let tokens = impg::sample_mem_bwt::encode_walk(&walk).unwrap();
        let features = enumerate_subwalks(&tokens).unwrap();
        let distinct: BTreeSet<&FeatureKey> = features.iter().collect();
        let record = RoutedRecord {
            tokens,
            multiplicity: 8,
            occurrences: BTreeMap::from([(1u32, 1u64)]),
            total_occurrences: 1,
            forward_positions: vec![(0usize, 150u64)],
            reverse_positions: Vec::new(),
        };
        let owner_sets = vec![BTreeSet::from([1u32])];
        let maps = build_in_window_obs(&[record], &index, 63, &owner_sets, None, None).unwrap();
        assert_eq!(maps[0].len(), distinct.len());
        for feature in &distinct {
            assert_eq!(maps[0][*feature], 8.0f64);
        }
    }

    /// The omission term is a valid negative log-likelihood: positive, and
    /// growing in the observed share, including for high-background (flood)
    /// features where the bare beta - C*ln(beta) form goes negative.
    #[test]
    fn omission_term_is_poisson_complete() {
        let mut cross: HashMap<FeatureKey, f64> = HashMap::new();
        let mut cross2: HashMap<FeatureKey, f64> = HashMap::new();
        cross.insert(vec![7u64], 467.9); // beta = 468.0 (a flood feature)
        let sample = SampleSideBackgrounds::new(&cross, 0.1);
        let mut window_obs: HashMap<FeatureKey, f64> = HashMap::new();
        window_obs.insert(vec![7u64], 50.0);
        let mut profile = Profile::new();
        let at_50 = omission_charge_addend(&profile, &window_obs, &sample);
        *window_obs.get_mut(&vec![7u64]).unwrap() = 500.0;
        let at_500 = omission_charge_addend(&profile, &window_obs, &sample);
        assert!(at_50 > 0.0, "omission term must be positive, got {at_50}");
        assert!(at_500 > 0.0, "omission term must stay positive, got {at_500}");
        // The bare form at beta=468, C=500 is 468 - 500*ln(468) < 0: the
        // log-factorial correction is what keeps the charge a valid
        // (positive, lower-bounded) negative log-likelihood. Observing C
        // near the background mean beta is genuinely cheap to omit — the
        // background EXPECTS that mass — so the term is not monotone in C
        // around beta; the discriminating cost sits at the ordinary
        // (low-background) features whose observed support the background
        // does not expect.
        assert!(468.0 - 500.0 * 468.0f64.ln() < 0.0);
        // A low-background feature's term: ~beta + |C*ln(beta)|, smaller
        // than the high-background feature's at the same share scale.
        cross2.insert(vec![9u64], 0.0);
        let sample2 = SampleSideBackgrounds::new(&cross2, 0.1);
        let mut window_obs2: HashMap<FeatureKey, f64> = HashMap::new();
        window_obs2.insert(vec![9u64], 5.0);
        let at_low = omission_charge_addend(&profile, &window_obs2, &sample2);
        // beta=0.1, C=5: 0.1 + |5*ln(0.1)| + lnGamma(6) ~= 16.5.
        assert!(at_low > 0.0 && at_low < 20.0);
    }

    /// Test-only TerritoryIndex assembled from raw steps (the routing's own
    /// entry filter and node table construction, for hand-built step lists).
    fn index_with_steps(steps: Vec<(u64, i32)>, interval: (u64, u64, u32)) -> TerritoryIndex {
        let k = 63u64;
        let (start, end, partition) = interval;
        let mut entries: Vec<(i32, u32, u64)> = steps
            .iter()
            .filter(|&&(bp, _)| bp + k > start && bp < end)
            .map(|&(bp, node)| (node, 0u32, bp))
            .collect();
        entries.sort_by_key(|&(node, _, _)| node_key(node));
        let max_node = entries
            .iter()
            .map(|&(node, _, _)| node)
            .max()
            .unwrap_or(i32::MAX);
        let table_len = node_key(max_node) + 2;
        let mut counts = vec![0u64; table_len];
        for &(node, _, _) in &entries {
            counts[node_key(node)] += 1;
        }
        let mut node_ranges = vec![0u64; table_len];
        let mut acc = 0u64;
        for key in 0..table_len {
            node_ranges[key] = acc;
            acc += counts[key];
        }
        TerritoryIndex {
            territories: vec![Territory {
                partition,
                path_idx: 0,
                start,
                end,
                steps,
            }],
            entries,
            node_ranges,
            path_intervals: vec![vec![interval]],
        }
    }

    /// THE RC-SYMMETRY TEST (the strand-asymmetry fix's contract): a
    /// canonical feature's k-mer and its reverse complement BOTH index —
    /// a path storing the k-mer forward yields the forward frame's step,
    /// a path storing the same feature's rc yields the rc frame's step with
    /// the SAME node identity the forward frame would assign — so a record
    /// places at BOTH loci, the forward one through its forward orientation
    /// (forward_positions) and the rc-stored one through its rc orientation
    /// (reverse_positions). Under the forward-only index the rc-stored locus
    /// is unplaceable (the pre-fix behavior the flood channel measured).
    #[test]
    fn rc_symmetric_indexing_places_kmer_and_rc_stored_locus() {
        // Node identity sanity of the mapping itself: the rc-view matched
        // syncmer (signed hash m, rc position q) maps to the path's forward
        // bp = lo + len - k - q with the node -m (the forward frame's own id
        // for the k-mer the path carries there).
        assert_eq!(rc_frame_step(-5, 37, 0, 300, 63), (200u64, 5i32));
        assert_eq!(rc_frame_step(7, 0, 1000, 463, 63), (1400u64, -7i32));

        // One canonical feature (node 5's k-mer mod rc): stored forward at
        // bp 100 (its k-mer qualifies, the forward frame's step), and stored
        // rc at bp 163 (only the rc frame qualifies; the rc-symmetric step
        // carries the rc-stored string's own forward-frame node id -5).
        let single_frame = index_with_steps(vec![(100u64, 5i32)], (0, 1000, 1));
        let both_frames = index_with_steps(vec![(100u64, 5i32), (163u64, -5i32)], (0, 1000, 1));
        let anchors = vec![(5i32, 0u64)];

        // Forward-only index: only the forward locus is placeable.
        let (occ, total, forward, reverse) = route_record(&anchors, &single_frame, 63, false);
        assert_eq!(total, 1);
        assert_eq!(forward, vec![(0usize, 100u64)]);
        assert!(reverse.is_empty());
        assert_eq!(occ, BTreeMap::from([(1u32, 1u64)]));

        // Both frames: the rc-stored locus becomes placeable via the
        // record's rc orientation, and the S2 spans' reverse_positions
        // carries it.
        let (occ, total, forward, reverse) = route_record(&anchors, &both_frames, 63, false);
        assert_eq!(total, 2);
        assert_eq!(forward, vec![(0usize, 100u64)]);
        assert_eq!(reverse, vec![(0usize, 163u64)]);
        assert_eq!(occ, BTreeMap::from([(1u32, 2u64)]));
    }

    /// The rc orientation's mirrored geometry (the alpha fix) plus the
    /// rc-symmetric steps place a MULTI-anchor record at its rc-stored
    /// occurrence: the rc walk's anchors (-node, span-k-rel) meet the
    /// augmented steps exactly at the mirrored positions.
    #[test]
    fn route_record_rc_orientation_mirrored_geometry_places_multi_anchor_rc_stored_walk() {
        // Forward record [(5,0), (9,10)], span 73. The path stores the
        // walk's rc with occurrence alignment X = 500: the rc walk is
        // [(-9, 0), (-5, 10)], so the steps are (500, -9) and (510, -5);
        // the forward locus sits at [100, 110) with the plain steps.
        let index = index_with_steps(
            vec![(100u64, 5i32), (110u64, 9i32), (500u64, -9i32), (510u64, -5i32)],
            (0, 1000, 1),
        );
        let anchors = vec![(5i32, 0u64), (9i32, 10u64)];
        let (occ, total, forward, reverse) = route_record(&anchors, &index, 63, false);
        assert_eq!(total, 2);
        assert_eq!(forward, vec![(0usize, 100u64)]);
        assert_eq!(reverse, vec![(0usize, 500u64)]);
        assert_eq!(occ, BTreeMap::from([(1u32, 2u64)]));
    }
}
