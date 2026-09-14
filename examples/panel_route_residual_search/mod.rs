//! B1 defaults plus opt-in B2 geometry; no global pricing/support certificate.
pub mod chain;
pub mod geometry;
pub mod oriented;
pub mod score;
use clap::Parser;
use geometry::Region;
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::SyngIndex,
};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fs::{self, File},
    io::{self, Write},
    path::PathBuf,
    time::Instant,
};
fn invalid(s: &str) -> io::Error {
    io::Error::other(s)
}
fn ensure(ok: bool, s: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(s))
    }
}
fn close(a: f64, b: f64) -> io::Result<()> {
    ensure(
        a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-7 * (1.0 + a.abs().max(b.abs())),
        "exact/public objective parity failed",
    )
}
#[derive(Parser, Clone, Debug, Serialize)]
pub struct Options {
    /// Experimental oriented exchanges and return-hub two-donor construction.
    #[arg(long)]
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    pub geometry_b2: bool,
    #[arg(long)]
    pub panel: String,
    #[arg(long)]
    pub routes: PathBuf,
    #[arg(long)]
    pub sample: PathBuf,
    #[arg(long)]
    pub out_dir: PathBuf,
    #[arg(long, default_value_t = 2_000_000)]
    pub max_work: u64,
    #[arg(long, default_value_t = 2_000_000)]
    pub max_profile_work: u64,
    #[arg(long, default_value_t = 2048)]
    pub max_scores: u64,
    #[arg(long, default_value_t = 1024)]
    pub max_validations: u64,
    #[arg(long, default_value_t = 128)]
    pub max_tasks: usize,
    #[arg(long, default_value_t = 64)]
    pub max_ties: usize,
    #[arg(long, default_value_t = 64)]
    pub max_segments: usize,
    #[arg(long, default_value_t = 50_000)]
    pub max_features: usize,
    #[arg(long, default_value_t = 134_217_728)]
    pub max_state_bytes: u64,
    #[arg(long, default_value_t = 4)]
    pub max_epochs: u8,
    #[arg(long, default_value_t = 8)]
    pub max_level: u8,
}
#[derive(Serialize)]
pub struct Meter {
    #[serde(skip)]
    pub limits: Options,
    pub work_used: u64,
    pub profile_work: u64,
    /// Charged admissions, including calls that fail without a complete objective.
    pub score_attempts: u64,
    pub public_score_attempts: u64,
    pub native_score_attempts: u64,
    pub exact_scores: u64,
    pub exact_delta_scores: u64,
    pub public_scores: u64,
    pub native_scores: u64,
    pub validations: u64,
    pub explicit_profile_attempts: u64,
    pub public_evaluation_route_reservations: u64,
    /// Successful route_counts returns plus routes in successful evaluate returns.
    /// Internal profile progress in a failed evaluate is not exposed by the API.
    pub completed_route_profiles_observed: u64,
    pub route_cache_hits: u64,
    pub corrected_starts: u64,
    pub port_reads: u64,
    pub source_index_opens: u64,
    pub state_bytes: u64,
    pub peak_state_bytes: u64,
    pub cleanup_work_prepaid: u64,
    #[serde(skip)]
    pub accounting: Option<File>,
}
impl Meter {
    pub fn new(o: &Options) -> Self {
        Self {
            limits: o.clone(),
            work_used: 0,
            profile_work: 0,
            score_attempts: 0,
            public_score_attempts: 0,
            native_score_attempts: 0,
            exact_scores: 0,
            exact_delta_scores: 0,
            public_scores: 0,
            native_scores: 0,
            validations: 0,
            explicit_profile_attempts: 0,
            public_evaluation_route_reservations: 0,
            completed_route_profiles_observed: 0,
            route_cache_hits: 0,
            corrected_starts: 0,
            port_reads: 0,
            source_index_opens: 0,
            state_bytes: 0,
            peak_state_bytes: 0,
            cleanup_work_prepaid: 0,
            accounting: None,
        }
    }
    pub fn work(&mut self, n: u64) -> io::Result<()> {
        ensure(
            n <= self.limits.max_work.saturating_sub(self.work_used),
            "budget: work",
        )?;
        self.work_used += n;
        Ok(())
    }
    pub fn reserve(&mut self, n: u64) -> io::Result<()> {
        ensure(
            n <= self.limits.max_state_bytes.saturating_sub(self.state_bytes),
            "budget: state_bytes",
        )?;
        let cleanup = n.div_ceil(64);
        self.work(cleanup)?;
        self.cleanup_work_prepaid += cleanup;
        self.state_bytes += n;
        self.peak_state_bytes = self.peak_state_bytes.max(self.state_bytes);
        if let Some(f) = &mut self.accounting {
            line(
                f,
                &json!({"kind":"state_reservation","bytes":n,"total":self.state_bytes,"cleanup_work_prepaid":cleanup}),
            )?;
        }
        Ok(())
    }
    pub fn profile(&mut self, length: u64, lengths: usize, label: &str) -> io::Result<()> {
        let n = length
            .checked_mul(lengths as u64)
            .ok_or_else(|| invalid("profile work overflow"))?;
        ensure(
            n <= self
                .limits
                .max_profile_work
                .saturating_sub(self.profile_work),
            "budget: profile_work",
        )?;
        self.profile_work += n;
        if let Some(f) = &mut self.accounting {
            line(
                f,
                &json!({"kind":"profile_charge","label":label,"complete_bp":length,"read_lengths":lengths,"units":n,"cumulative":self.profile_work}),
            )?;
        }
        Ok(())
    }
    pub fn score_room(&self, n: u64) -> io::Result<()> {
        ensure(
            n <= self.limits.max_scores.saturating_sub(self.score_attempts),
            "budget: exact_scores",
        )
    }
    pub fn admit_score(&mut self) -> io::Result<()> {
        self.score_room(1)?;
        self.score_attempts += 1;
        Ok(())
    }
    pub fn validation(&mut self, a: &routes::Assignment, k: u64, ports: u64) -> io::Result<()> {
        ensure(
            self.validations < self.limits.max_validations,
            "budget: validations",
        )?;
        ensure(
            a.routes
                .iter()
                .all(|r| r.segments.len() <= self.limits.max_segments),
            "budget: route_segments",
        )?;
        let n = a
            .routes
            .iter()
            .map(|r| r.segments.len() as u64)
            .sum::<u64>();
        self.work(n * n + n * (2 * (64 - ports.max(1).leading_zeros() as u64) + k + 8))?;
        self.validations += 1;
        Ok(())
    }
}
#[derive(Clone, Debug, Serialize)]
enum Task {
    Coarse {
        base: usize,
        slot: usize,
        level: u8,
        cell: u64,
    },
    Region {
        base: usize,
        region: Region,
        guided: bool,
    },
    Pairs {
        base: usize,
        region: Region,
        guided: bool,
        lo: u64,
        hi: u64,
        left: (u64, u64),
        right: (u64, u64),
        i: u64,
        j: u64,
    },
    Candidate {
        base: usize,
        region: Region,
        guided: bool,
        lo: u64,
        hi: u64,
        donor: routes::Segment,
        other: usize,
    },
}
#[derive(Clone, Serialize)]
struct Advice {
    score: f64,
    region: Region,
}
/// Match only the unchanged public profile/union cap diagnostics at this adapter.
/// Other public validation, provenance and numerical errors remain fatal.
pub fn public_feature_limit(error: io::Error) -> io::Error {
    if error.kind() == io::ErrorKind::InvalidData
        && matches!(
            error.to_string().as_str(),
            "native feature resource cap exhausted (no truncation)"
                | "realized feature union exceeds resource cap; no partial objective"
        )
    {
        invalid(&format!("budget: public_feature_terms: {error}"))
    } else {
        error
    }
}
fn line(w: &mut File, v: &impl Serialize) -> io::Result<()> {
    serde_json::to_writer(&mut *w, v)?;
    w.write_all(b"\n")
}
struct Engine<'a, 'b> {
    e: &'a mut routes::Evaluator<'b>,
    meter: Meter,
    observed: BTreeMap<[u64; 3], u64>,
    histogram: Vec<u64>,
    denominator: u64,
    bases: Vec<Box<score::Baseline>>,
    tasks: VecDeque<Task>,
    advice: Vec<Vec<Advice>>,
    batch: Vec<usize>,
    family_best: Vec<f64>,
    epochs: Vec<u8>,
    ties: Vec<routes::Assignment>,
    best: f64,
    ledger: File,
    events: File,
    active: Option<Task>,
    native_done: usize,
    confirmed: usize,
    refinements: usize,
    conflicts: usize,
    compound_confirmed: usize,
    unsupported_opposite_strand: usize,
    chains: VecDeque<chain::Cursor>,
    active_chain: Option<chain::Cursor>,
    next_chain_key: (usize, usize),
    chain_primitives: u64,
    chain_primitives_completed: u64,
    chain_confirmed: usize,
    oriented_rejections: usize,
    stop: String,
}
impl<'a, 'b> Engine<'a, 'b> {
    fn push(&mut self, t: Task) -> io::Result<()> {
        self.meter.work(1)?;
        ensure(
            self.tasks.len() < self.meter.limits.max_tasks,
            "budget: pending_tasks",
        )?;
        self.tasks.push_back(t);
        Ok(())
    }
    fn public(&mut self, a: &routes::Assignment, native: bool) -> io::Result<routes::Evaluation> {
        self.meter.score_room(1)?;
        self.meter
            .validation(a, self.e.graph.k, self.e.graph.port_count)?;
        let length = a.routes.iter().try_fold(0u64, |n, r| {
            n.checked_add(r.length()?)
                .ok_or_else(|| invalid("genome length overflow"))
        })?;
        self.meter.profile(
            length,
            self.histogram.len(),
            if native {
                "public_native"
            } else {
                "public_confirmation"
            },
        )?;
        self.meter.admit_score()?;
        self.meter.public_score_attempts += 1;
        if native {
            self.meter.native_score_attempts += 1;
        }
        self.meter.public_evaluation_route_reservations += a.routes.len() as u64;
        let hits = self.e.cache_hits;
        let attempted = self.e.evaluate(a);
        self.meter.route_cache_hits += self.e.cache_hits - hits;
        let actual = attempted.map_err(public_feature_limit)?;
        self.meter.exact_scores += 1;
        self.meter.public_scores += 1;
        if native {
            self.meter.native_scores += 1;
        }
        self.meter.completed_route_profiles_observed += a.routes.len() as u64;
        self.meter.corrected_starts += actual.corrected_starts;
        Ok(actual)
    }
    fn remember(&mut self, e: &routes::Evaluation) -> io::Result<()> {
        if e.relative_objective < self.best - 1e-9 {
            self.best = e.relative_objective;
            self.ties.clear();
        }
        if (e.relative_objective - self.best).abs() <= 1e-9 && !self.ties.contains(&e.assignment) {
            ensure(
                self.ties.len() < self.meter.limits.max_ties,
                "budget: retained_ties (full scored ledger retained)",
            )?;
            self.ties.push(e.assignment.clone());
        }
        Ok(())
    }
    fn add_base(&mut self, e: routes::Evaluation, epoch: u8) -> io::Result<usize> {
        let q = score::baseline(e, epoch, self.histogram.len());
        let cells = q.counts.iter().map(|q| q.len()).sum::<usize>();
        let seg = q
            .assignment
            .routes
            .iter()
            .map(|r| r.segments.capacity())
            .sum::<usize>();
        self.meter.reserve(
            (cells * 160 + seg * 64 + q.assignment.routes.capacity() * 128 + 4096) as u64,
        )?;
        let id = self.bases.len();
        line(
            &mut self.events,
            &json!({"event":"baseline","baseline":id,"epoch":epoch,"assignment":q.assignment,"objective":q.score}),
        )?;
        self.bases.push(Box::new(q));
        self.advice.push(Vec::with_capacity(2));
        self.batch.push(0);
        self.push(Task::Coarse {
            base: id,
            slot: 0,
            level: 0,
            cell: 0,
        })?;
        Ok(id)
    }
    fn initialize(&mut self, family: usize) -> io::Result<()> {
        let a = self.e.graph.native_assignment(family)?;
        let actual = self.public(&a, true)?;
        line(&mut self.ledger, &json!({"kind":"native","public":actual}))?;
        self.remember(&actual)?;
        self.family_best[family] = actual.relative_objective;
        self.native_done += 1;
        self.add_base(actual, 0)?;
        Ok(())
    }
    fn feedback(&mut self, base: usize, region: &Region, objective: f64) -> io::Result<()> {
        self.batch[base] += 1;
        let advice = Advice {
            score: objective,
            region: region.clone(),
        };
        let rows = &mut self.advice[base];
        rows.push(advice);
        rows.sort_by(|a, b| {
            a.score
                .total_cmp(&b.score)
                .then(a.region.left.cmp(&b.region.left))
                .then(a.region.right.cmp(&b.region.right))
        });
        if rows.len() > 2 {
            rows.pop();
        } // Only bounded refinement summaries; every score stays in ledger.
        if self.batch[base] % 8 == 0 {
            let selected = std::mem::replace(&mut self.advice[base], Vec::with_capacity(2));
            for a in selected {
                let children = a.region.refinements(
                    self.bases[base].assignment.routes[a.region.slot].length()?,
                    self.meter.limits.max_level,
                );
                line(
                    &mut self.events,
                    &json!({"event":"refinement_choice","baseline":base,"batch":self.batch[base],"selected":a,"children":children,"non_improving_eligible":true}),
                )?;
                for region in children {
                    self.push(Task::Region {
                        base,
                        region,
                        guided: true,
                    })?;
                    self.refinements += 1;
                }
            }
        }
        Ok(())
    }
    fn step(&mut self, t: Task) -> io::Result<()> {
        self.meter.work(1)?;
        match t {
            Task::Coarse {
                base,
                slot,
                level,
                cell,
            } => {
                let lengths = self.bases[base].assignment.routes.len();
                if slot >= lengths {
                    return Ok(());
                }
                let n = 1u64 << level;
                let length = self.bases[base].assignment.routes[slot].length()?;
                let region = Region {
                    slot,
                    left: length * cell / (2 * n),
                    right: length * (cell + 2) / (2 * n),
                    level,
                };
                let (next_slot, next_level, next_cell) = if cell + 1 < 2 * n - 1 {
                    (slot, level, cell + 1)
                } else if level < self.meter.limits.max_level {
                    (slot, level + 1, 0)
                } else {
                    (slot + 1, 0, 0)
                };
                if region.left < region.right {
                    self.push(Task::Region {
                        base,
                        region,
                        guided: false,
                    })?;
                }
                self.push(Task::Coarse {
                    base,
                    slot: next_slot,
                    level: next_level,
                    cell: next_cell,
                })?;
            }
            Task::Region {
                base,
                region,
                guided,
            } => {
                let route = &self.bases[base].assignment.routes[region.slot];
                let left = geometry::nearest(self.e, route, region.left, &mut self.meter)?;
                let right = geometry::nearest(self.e, route, region.right, &mut self.meter)?;
                line(
                    &mut self.events,
                    &json!({"event":"geometry_attempt","baseline":base,"region":region,"guided":guided,"left":left,"right":right}),
                )?;
                if let (Some((lo, p)), Some((hi, q))) = (left, right) {
                    if lo < hi {
                        let left = geometry::bucket(self.e, &p.word, &mut self.meter)?;
                        let right = geometry::bucket(self.e, &q.word, &mut self.meter)?;
                        self.push(Task::Pairs {
                            base,
                            region,
                            guided,
                            lo,
                            hi,
                            left,
                            right,
                            i: left.0,
                            j: right.0,
                        })?;
                    }
                }
            }
            Task::Pairs {
                base,
                region,
                guided,
                lo,
                hi,
                left,
                right,
                i,
                j,
            } => {
                if i >= left.1 || j >= right.1 {
                    return Ok(());
                }
                let a = geometry::read_port(
                    &mut self.e.ports.global,
                    self.e.graph.k,
                    i,
                    &mut self.meter,
                )?;
                let b = geometry::read_port(
                    &mut self.e.ports.global,
                    self.e.graph.k,
                    j,
                    &mut self.meter,
                )?;
                let (ni, nj) = if j + 1 < right.1 {
                    (i, j + 1)
                } else {
                    (i + 1, right.0)
                };
                self.push(Task::Pairs {
                    base,
                    region: region.clone(),
                    guided,
                    lo,
                    hi,
                    left,
                    right,
                    i: ni,
                    j: nj,
                })?;
                if a.source == b.source
                    && a.reverse == b.reverse
                    && if a.reverse {
                        a.cut > b.cut
                    } else {
                        a.cut < b.cut
                    }
                {
                    let donor = routes::Segment {
                        source: a.source,
                        start: a.cut.min(b.cut),
                        end: a.cut.max(b.cut),
                        reverse: a.reverse,
                    };
                    self.push(Task::Candidate {
                        base,
                        region,
                        guided,
                        lo,
                        hi,
                        donor,
                        other: 0,
                    })?;
                }
            }
            Task::Candidate {
                base,
                region,
                guided,
                lo,
                hi,
                donor,
                other,
            } => {
                let assignment = &self.bases[base].assignment;
                let slots = assignment.routes.len();
                self.meter.work(
                    4 * assignment
                        .routes
                        .iter()
                        .map(|r| r.segments.len() as u64)
                        .sum::<u64>()
                        + 8,
                )?;
                let candidate = if other == 0 {
                    let mut a = assignment.clone();
                    a.routes[region.slot] =
                        geometry::splice(&a.routes[region.slot], lo, hi, vec![donor.clone()])?;
                    Some(a)
                } else if self.meter.limits.geometry_b2 {
                    oriented::reciprocal_proposal(
                        assignment,
                        region.slot,
                        lo,
                        hi,
                        &donor,
                        other - 1,
                    )?
                } else {
                    match geometry::reciprocal(assignment, region.slot, lo, hi, &donor, other - 1)?
                    {
                        geometry::Reciprocal::Supported(a) => Some(a),
                        geometry::Reciprocal::NotContained => None,
                        geometry::Reciprocal::UnsupportedOppositeStrand => {
                            self.unsupported_opposite_strand += 1;
                            line(
                                &mut self.events,
                                &json!({"event":"unsupported_opposite_strand_reciprocal","baseline":base,"region":region,"donor":donor,"other_slot":other-1,"biological_infeasibility":false}),
                            )?;
                            None
                        }
                    }
                };
                // The reciprocal successor is retained regardless of hypothetical single validity.
                if other < slots {
                    self.push(Task::Candidate {
                        base,
                        region: region.clone(),
                        guided,
                        lo,
                        hi,
                        donor: donor.clone(),
                        other: other + 1,
                    })?;
                }
                let Some(candidate) = candidate else {
                    return Ok(());
                };
                self.confirm_candidate(base, region, guided, lo, hi, other, candidate)?;
            }
        }
        Ok(())
    }
    fn confirm_candidate(
        &mut self,
        base: usize,
        region: Region,
        guided: bool,
        lo: u64,
        hi: u64,
        other: usize,
        candidate: routes::Assignment,
    ) -> io::Result<()> {
        if candidate == self.bases[base].assignment {
            return Ok(());
        }
        if self.meter.limits.geometry_b2 {
            if let Some(rejection) = oriented::seam_rejection(self.e, &candidate, &mut self.meter)?
            {
                self.oriented_rejections += 1;
                line(
                    &mut self.ledger,
                    &json!({"kind":"oriented_construction_rejection","baseline":base,"region":region,"assignment":candidate,"rejection":rejection,"global_infeasibility_claim":false}),
                )?;
                return Ok(());
            }
        }
        self.meter.score_room(2)?;
        self.meter
            .validation(&candidate, self.e.graph.k, self.e.graph.port_count)?;
        let candidate = match self.e.validate_assignment(&candidate) {
            Ok(a) => a,
            Err(e) if e.to_string().contains("conflicting canonical") => {
                self.conflicts += 1;
                line(
                    &mut self.ledger,
                    &json!({"kind":"capacity_rejection","baseline":base,"region":region,"compound":other>0,"assignment":candidate,"error":e.to_string()}),
                )?;
                return Ok(());
            }
            Err(e) => {
                line(
                    &mut self.ledger,
                    &json!({"kind":"unexpected_validation_error","task":self.active,"baseline":base,"candidate":candidate,"error":e.to_string()}),
                )?;
                return Err(e);
            }
        };
        let change = score::delta(
            self.e,
            &self.bases[base],
            &candidate,
            &self.observed,
            &self.histogram,
            self.denominator,
            &mut self.meter,
        )?;
        line(
            &mut self.ledger,
            &json!({"kind":"unconfirmed_exact_delta","baseline":base,"assignment":candidate,"delta":change,"region":region,"compound":other>0}),
        )?;
        let actual = self.public(&candidate, false)?;
        line(
            &mut self.ledger,
            &json!({"kind":"public_confirmation","baseline":base,"public":actual}),
        )?;
        score::parity(
            &self.bases[base],
            &change,
            &actual,
            &self.observed,
            &mut self.meter,
        )?;
        self.confirmed += 1;
        if other > 0 {
            self.compound_confirmed += 1;
        }
        line(
            &mut self.ledger,
            &json!({"kind":"candidate","baseline":base,"baseline_objective":self.bases[base].score,"region":region,"snapped":[lo,hi],"guided":guided,"compound":other>0,"delta":change,"public_minus_sparse":actual.relative_objective-change.objective,"public":actual}),
        )?;
        self.remember(&actual)?;
        self.feedback(base, &region, actual.relative_objective)?;
        let family = actual.assignment.family;
        if actual.relative_objective < self.family_best[family] - 1e-9 {
            self.family_best[family] = actual.relative_objective;
            if self.epochs[family] < self.meter.limits.max_epochs {
                self.epochs[family] += 1;
                self.add_base(actual, self.epochs[family])?;
            }
        }
        Ok(())
    }
    fn admit_chains(&mut self) -> io::Result<()> {
        if !self.meter.limits.geometry_b2 {
            return Ok(());
        }
        while self.chains.len() < 32 && self.next_chain_key.0 < self.bases.len() {
            self.meter.work(1)?;
            let (base, slot) = self.next_chain_key;
            if slot == self.bases[base].assignment.routes.len() {
                self.next_chain_key = (base + 1, 0);
            } else {
                self.chains.push_back(chain::Cursor::new(base, slot));
                self.next_chain_key.1 += 1;
            }
        }
        Ok(())
    }
    fn chain_step(&mut self) -> io::Result<()> {
        let mut cursor = self.chains.pop_front().unwrap();
        self.active_chain = Some(cursor.clone());
        ensure(self.chain_primitives < 100_000, "budget: chain_primitives")?;
        self.meter.work(1)?;
        self.chain_primitives += 1;
        let (done, candidate) =
            cursor.advance(self.e, &self.bases[cursor.base].assignment, &mut self.meter)?;
        self.chain_primitives_completed += 1;
        let base = cursor.base;
        if !done {
            self.chains.push_back(cursor);
        }
        if let Some((region, candidate)) = candidate {
            line(
                &mut self.events,
                &json!({"event":"chain_complete_geometry","baseline":base,"cursor":self.active_chain,"region":region,"assignment":candidate,"right_discovered_from_return_hub":true}),
            )?;
            let before = self.confirmed;
            let (lo, hi) = (region.left, region.right);
            let confirmation = self.confirm_candidate(base, region, false, lo, hi, 0, candidate);
            // Confirmation precedes fallible retention/feedback. Classify any
            // completed public/parity result even when that later work stops.
            self.chain_confirmed += self.confirmed - before;
            confirmation?;
        }
        self.active_chain = None;
        Ok(())
    }
}
/// mechanism_family is TEST-ONLY injected fixed baseline, never a CLI parameter.
pub fn run(options: Options, mechanism_family: Option<usize>) -> io::Result<serde_json::Value> {
    fs::create_dir(&options.out_dir)?;
    genome::write_json(&options.out_dir.join("options.json"), &options)?;
    genome::write_json(
        &options.out_dir.join("status.json"),
        &json!({"status":"running","global_bound":null}),
    )?;
    let start = Instant::now();
    let result = (|| {
        ensure(
            options.max_level <= 8
                && options.max_epochs <= 4
                && options.max_features > 0
                && options.max_features <= 50_000
                && options.max_tasks > 0
                && options.max_tasks <= 128
                && options.max_ties > 0
                && options.max_ties <= 64
                && options.max_segments <= 64
                && options.max_work <= 2_000_000
                && options.max_profile_work <= 2_000_000
                && options.max_scores <= 2048
                && options.max_validations <= 1024
                && options.max_state_bytes <= 134_217_728,
            "invalid experimental limits",
        )?;
        let metadata = fs::metadata(options.routes.join("graph.json"))?.len();
        ensure(
            metadata.saturating_mul(8) < options.max_state_bytes,
            "budget: graph_metadata",
        )?;
        let identity = genome::PanelIdentity::read(&options.panel)?;
        let panel = SyngIndex::load(&options.panel, Default::default())?;
        let graph = routes::Graph::load(&options.routes, &identity)?;
        let sample = sample::SampleIndex::load(&options.sample, &identity)?;
        let observed = sample.counts.observed_pairs()?;
        ensure(
            observed.len() <= options.max_features,
            "budget: observed_terms",
        )?;
        ensure(
            observed.values().all(|&n| n <= 1 << 53),
            "observed integer conversion precision limit",
        )?;
        let histogram: Vec<_> = graph
            .read_lengths
            .iter()
            .map(|l| {
                sample
                    .stats
                    .read_lengths
                    .get(&(*l as usize))
                    .copied()
                    .unwrap_or(0)
            })
            .collect();
        let denominator =
            graph
                .read_lengths
                .iter()
                .zip(&histogram)
                .try_fold(0u64, |n, (&l, &h)| {
                    n.checked_add(
                        l.checked_mul(h)
                            .ok_or_else(|| invalid("histogram overflow"))?,
                    )
                    .ok_or_else(|| invalid("histogram overflow"))
                })?;
        let mut meter = Meter::new(&options);
        meter.accounting = Some(File::create(options.out_dir.join("accounting.jsonl"))?);
        // Prepay public opaque profile/cache temporaries, sparse map clones, candidate
        // transfers, retained queue capacities and support slots. Not a native RSS bound.
        meter.reserve(
            metadata * 8
                + (options.max_features as u64 + 1)
                    * 160
                    * (graph.read_lengths.len() as u64 * 6 + 4)
                + (options.max_tasks as u64) * 4096
                + (options.max_ties as u64)
                    * graph.lanes.len() as u64
                    * options.max_segments as u64
                    * 64
                + graph.families.len() as u64 * 4096
                + graph.lanes.len() as u64 * 1024,
        )?;
        if options.geometry_b2 {
            // 32 queued cursors plus active and transfer copies. Each cursor
            // owns at most three k-byte words; reserve eight word buffers per
            // slot for clones/read/RC temporaries, plus fixed fields/allocator
            // capacity. Candidate transfers stay in B1's assignment workspace.
            meter.reserve(
                34 * (8192
                    + graph
                        .k
                        .checked_mul(8)
                        .ok_or_else(|| invalid("cursor word size overflow"))?),
            )?;
        }
        let mut evaluator = routes::Evaluator::new(
            &options.routes,
            &graph,
            &panel,
            &sample,
            10.0,
            0.1,
            options.max_features,
            options.max_features,
        )?;
        genome::write_json(
            &options.out_dir.join("provenance.json"),
            &json!({"graph":graph.digest()?,"compiler":graph.compiler_identity,"count_policy":graph.count_policy,"sample":genome::reconstruction::fingerprint(&options.sample)?,"histogram":histogram,"read_lengths":graph.read_lengths,"verified_source_bytes":graph.source_files.iter().map(|f|f["bytes"].as_u64().unwrap_or(0)).sum::<u64>(),"global_port_bytes":graph.ports.bytes,"mechanism_only_fixed_family":mechanism_family,"fractional_genomes":false,"new_two_donor_chain_generation":options.geometry_b2}),
        )?;
        let nf = graph.families.len();
        let max_bases = nf * (options.max_epochs as usize + 1);
        let mut engine = Engine {
            e: &mut evaluator,
            meter,
            observed,
            histogram,
            denominator,
            bases: Vec::with_capacity(max_bases),
            tasks: VecDeque::with_capacity(options.max_tasks),
            advice: Vec::with_capacity(max_bases),
            batch: Vec::with_capacity(max_bases),
            family_best: vec![f64::INFINITY; nf],
            epochs: vec![0; nf],
            ties: Vec::with_capacity(options.max_ties),
            best: f64::INFINITY,
            ledger: File::create(options.out_dir.join("scores.jsonl"))?,
            events: File::create(options.out_dir.join("events.jsonl"))?,
            active: None,
            native_done: 0,
            confirmed: 0,
            refinements: 0,
            conflicts: 0,
            compound_confirmed: 0,
            unsupported_opposite_strand: 0,
            chains: VecDeque::with_capacity(if options.geometry_b2 { 32 } else { 0 }),
            active_chain: None,
            next_chain_key: (0, 0),
            chain_primitives: 0,
            chain_primitives_completed: 0,
            chain_confirmed: 0,
            oriented_rejections: 0,
            stop: "declared_neighborhood_schedule_finished_not_global".into(),
        };
        let families: Vec<_> = mechanism_family
            .map(|f| vec![f])
            .unwrap_or_else(|| (0..nf).collect());
        if mechanism_family.is_some() {
            ensure(options.max_epochs == 0, "mechanism must fix baseline epoch")?;
        }
        let mut next = 0;
        let mut init_turn = true;
        let mut chain_turn = false;
        let mut active_native_family = None;
        loop {
            if let Err(e) = engine.admit_chains() {
                if e.to_string().starts_with("budget:") {
                    engine.stop = e.to_string();
                    break;
                }
                return Err(e);
            }
            let servicing_chain = !engine.chains.is_empty()
                && (chain_turn || (next == families.len() && engine.tasks.is_empty()));
            let attempt = if servicing_chain {
                engine.chain_step()
            } else if next < families.len() && (init_turn || engine.tasks.is_empty()) {
                let family = families[next];
                active_native_family = Some(family);
                engine.initialize(family).map(|()| {
                    next += 1;
                    active_native_family = None;
                })
            } else if let Some(task) = engine.tasks.pop_front() {
                engine.active = Some(task.clone());
                engine.step(task)
            } else {
                break;
            };
            // Chain turns must not consume the ordinary channel's independent
            // native/task alternation (B1 has no chain turns).
            if !servicing_chain {
                init_turn = !init_turn;
            }
            chain_turn = !chain_turn;
            if let Err(e) = attempt {
                if e.to_string().starts_with("budget:") {
                    engine.stop = e.to_string();
                    break;
                }
                return Err(e);
            }
            engine.active = None;
        }
        engine.ledger.flush()?;
        engine.events.flush()?;
        engine.meter.accounting.as_mut().unwrap().flush()?;
        let mut pending = json!({"active_at_stop":engine.active,"active_native_family":active_native_family,"pending_native_families":&families[next..],"queued":engine.tasks,"immutable_baselines":engine.bases.iter().enumerate().map(|(id,b)|json!({"id":id,"epoch":b.epoch,"assignment":b.assignment,"objective":b.score})).collect::<Vec<_>>(),"resume_supported":false});
        if options.geometry_b2 {
            pending["b2"] = json!({"active_chain":engine.active_chain,"queued_chains":engine.chains,"next_chain_admission_key":engine.next_chain_key,"primitive_admissions":engine.chain_primitives,"live_cursor_cap":32,"key_admission_is_lazy":true});
        }
        genome::write_json(&options.out_dir.join("pending.json"), &pending)?;
        let mut result = json!({"status":engine.stop,"native_initializations":engine.native_done,"native_initialization_complete":next==families.len(),"confirmed_candidates":engine.confirmed,"compound_confirmations":engine.compound_confirmed,"unsupported_opposite_strand_reciprocals":engine.unsupported_opposite_strand,"capacity_rejections":engine.conflicts,"refinement_tasks_generated":engine.refinements,"best_public_objective":engine.best.is_finite().then_some(engine.best),"retained_correlated_best_found":engine.ties,"meter":engine.meter,"elapsed_seconds":start.elapsed().as_secs_f64(),"pending_tasks":engine.tasks.len(),"global_bound":null,"global_gap":null,"correlated_support_complete":false,"sequence_emission_authorized":false,"remaining_gate_b":"oriented coupled construction, two-donor/no-simple-return generation and reviewed genome-scale development experiments","remaining_finer_work":true});
        if options.geometry_b2 {
            result["b2"] = json!({"chain_primitive_admissions":engine.chain_primitives,"chain_primitives_completed":engine.chain_primitives_completed,"chain_confirmations":engine.chain_confirmed,"oriented_construction_rejections":engine.oriented_rejections,"live_chains":engine.chains.len(),"transported_interior_positive_witness_established":false,"diploid_dosage_or_phasing_validated":false,"full_b2_complete":false});
            result["remaining_gate_b"]=json!("positive transported-interior witness, full B2 acceptance and reviewed genome-scale development experiments");
        }
        genome::write_json(&options.out_dir.join("result.json"), &result)?;
        Ok(result)
    })();
    genome::write_json(
        &options.out_dir.join("status.json"),
        &match &result {
            Ok(v) => {
                json!({"status":"succeeded_bounded_experiment","stop":v["status"],"global_bound":null})
            }
            Err(e) => {
                json!({"status":"failed","error":e.to_string(),"global_bound":null,"support_complete":false})
            }
        },
    )?;
    result
}
