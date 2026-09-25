//! Fixed two ordered whole-copy baselines; shared geometry, complete paired feedback.
#[path = "genome.rs"]
pub mod genome_wide;
#[allow(dead_code)]
#[path = "../panel_route_residual_search/mod.rs"]
pub mod haploid;
pub mod partition;
pub mod score;
use clap::Parser;
use haploid::{
    chain,
    continuation::{self, Task, TaskQueue},
    geometry::Region,
    oriented,
};
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::SyngIndex,
};
use score::{Pair, Scored};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, VecDeque},
    fs::{self, File},
    io::{self, Write},
    path::PathBuf,
    time::Instant,
};
fn invalid(s: &str) -> io::Error {
    io::Error::other(s)
}
const JOINT_PROJECTION_WINDOW: u64 = 16;
fn ensure(ok: bool, s: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(s))
    }
}
fn line(f: &mut File, v: &impl Serialize) -> io::Result<()> {
    serde_json::to_writer(&mut *f, v)?;
    f.write_all(b"\n")
}
#[derive(Clone, Debug, Parser, Serialize)]
pub struct Options {
    #[arg(long)]
    pub panel: String,
    #[arg(long)]
    pub routes: PathBuf,
    #[arg(long)]
    pub sample: PathBuf,
    #[arg(long)]
    pub out_dir: PathBuf,
    #[arg(long, default_value_t = 4_000_000)]
    pub max_work: u64,
    #[arg(long, default_value_t = 8_000_000)]
    pub max_profile_work: u64,
    #[arg(long, default_value_t = 2048)]
    pub max_scores: u64,
    #[arg(long, default_value_t = 4096)]
    pub max_haploid_scores: u64,
    #[arg(long, default_value_t = 8192)]
    pub max_validations: u64,
    #[arg(long, default_value_t = 128)]
    pub max_tasks: usize,
    #[arg(long, default_value_t = 64)]
    pub max_ties: usize,
    #[arg(long, default_value_t = 32)]
    pub max_baselines: usize,
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
    #[arg(long, default_value_t = 32)]
    pub max_cursors: usize,
    #[arg(long, default_value_t = 100_000)]
    pub max_chain_primitives: u64,
}
impl Options {
    pub fn validate(&self) -> io::Result<()> {
        ensure(
            self.max_work <= 4_000_000
                && self.max_profile_work <= 8_000_000
                && self.max_scores <= 2048
                && self.max_haploid_scores <= 4096
                && self.max_validations <= 8192
                && self.max_tasks > 0
                && self.max_tasks <= 128
                && self.max_ties > 0
                && self.max_ties <= 64
                && self.max_baselines > 0
                && self.max_baselines <= 32
                && self.max_segments <= 64
                && self.max_features > 0
                && self.max_features <= 50_000
                && self.max_state_bytes <= 134_217_728
                && self.max_epochs <= 4
                && self.max_level <= 8
                && self.max_cursors > 0
                && self.max_cursors <= 32
                && self.max_chain_primitives <= 100_000,
            "invalid fixed paired experimental limits",
        )
    }
    fn meter_options(&self) -> haploid::Options {
        haploid::Options {
            geometry_b2: true,
            panel: self.panel.clone(),
            routes: self.routes.clone(),
            sample: self.sample.clone(),
            out_dir: self.out_dir.clone(),
            max_work: self.max_work,
            max_profile_work: self.max_profile_work,
            max_scores: self.max_haploid_scores,
            max_validations: self.max_validations,
            max_tasks: self.max_tasks,
            max_ties: self.max_ties,
            max_segments: self.max_segments,
            max_features: self.max_features,
            max_state_bytes: self.max_state_bytes,
            max_epochs: self.max_epochs,
            max_level: self.max_level,
        }
    }
}
#[derive(Clone, Serialize)]
struct PairedTask {
    copy: usize,
    task: Task,
}
/// Only metadata lives in a deferred producer. Complete scores remain in baselines.
#[derive(Clone, Default, Serialize)]
struct Emission {
    receipt: u64,
    source_base: Option<usize>,
    batch: usize,
    advice: [Option<Advice>; 2],
    advice_index: usize,
    child_index: usize,
    seeds: [Option<usize>; 2],
    seed_index: usize,
    seed_copy: usize,
    seed_turn: bool,
}
impl Emission {
    fn has_intents(&self) -> bool {
        self.advice.iter().any(Option::is_some) || self.seeds.iter().any(Option::is_some)
    }
}
#[derive(Clone, Serialize)]
struct NativeHead {
    families: [usize; 2],
    next: Option<[usize; 2]>,
    seed: Option<Emission>,
}
#[derive(Clone, Serialize)]
enum Work {
    Geometry(PairedTask),
    Emission(Emission),
    Native(NativeHead),
}
impl Work {
    fn class(&self) -> usize {
        match self {
            Self::Native(n) if n.seed.is_none() => 0,
            Self::Geometry(_) => 1,
            _ => 2,
        }
    }
}
/// The shared sink sees only its reserved push quota, not unreserved free space.
/// Its original input is an unfinished remainder, never an extra runnable task.
struct BoundQueue<'a> {
    copy: usize,
    tasks: &'a mut VecDeque<Work>,
    credits: &'a mut usize,
    emitted: &'a mut usize,
    quota: usize,
    cap: usize,
}
impl TaskQueue for BoundQueue<'_> {
    fn len(&self) -> usize {
        self.cap - (self.quota - *self.emitted)
    }
    fn push_back(&mut self, task: Task) {
        assert!(*self.emitted < self.quota && *self.credits > 0);
        self.tasks.push_back(Work::Geometry(PairedTask {
            copy: self.copy,
            task,
        }));
        *self.credits -= 1;
        *self.emitted += 1;
    }
}
#[derive(Clone, Serialize)]
struct Advice {
    score: f64,
    copy: usize,
    region: Region,
}
#[derive(Serialize)]
struct Baseline {
    id: usize,
    scored: Scored,
    depth: u8,
    batch: usize,
    advice: Vec<Advice>,
    exploratory: Vec<Scored>,
    seed_owner: u64,
    seed_emitted: [bool; 2],
    native_seen: bool,
    geometry_services: u64,
}
#[derive(Clone, Serialize)]
struct PairedChain {
    copy: usize,
    cursor: chain::Cursor,
}
#[derive(Serialize)]
struct Receipt {
    scored: Scored,
    stage: &'static str,
}
struct EmissionPlan {
    next: Emission,
    action: Option<(usize, usize, Region, bool)>, // base, copy, region, seed
    remaining: usize,
}
fn next_native(pair: [usize; 2], families: usize, fixed: bool) -> Option<[usize; 2]> {
    if fixed {
        None
    } else if pair[1] + 1 < families {
        Some([pair[0], pair[1] + 1])
    } else if pair[0] + 1 < families {
        Some([pair[0] + 1, pair[0] + 1])
    } else {
        None
    }
}
/// Physical equality only. Count-equivalent assignments remain expansion-eligible.
pub fn expansion_eligible(pair: &Pair, admitted: &[Pair]) -> bool {
    let key = pair.key();
    !admitted.iter().any(|p| p.key() == key)
}
fn route_slice(route: &routes::Route, from: u64, to: u64) -> io::Result<Vec<routes::Segment>> {
    ensure(
        from <= to && to <= route.length()?,
        "invalid route projection slice",
    )?;
    let mut offset = 0;
    let mut out = Vec::new();
    for segment in &route.segments {
        let length = segment.end - segment.start;
        let left = from.saturating_sub(offset).min(length);
        let right = to.saturating_sub(offset).min(length);
        if left < right {
            out.push(if segment.reverse {
                routes::Segment {
                    source: segment.source,
                    start: segment.end - right,
                    end: segment.end - left,
                    reverse: true,
                }
            } else {
                routes::Segment {
                    source: segment.source,
                    start: segment.start + left,
                    end: segment.start + right,
                    reverse: false,
                }
            });
        }
        offset += length;
    }
    ensure(
        offset == route.length()?,
        "incomplete route projection slice",
    )?;
    Ok(out)
}
fn segment_slice(segment: &routes::Segment, from: u64, to: u64) -> io::Result<routes::Segment> {
    let length = segment.end - segment.start;
    ensure(
        from < to && to <= length,
        "invalid physical segment projection slice",
    )?;
    Ok(if segment.reverse {
        routes::Segment {
            source: segment.source,
            start: segment.end - to,
            end: segment.end - from,
            reverse: true,
        }
    } else {
        routes::Segment {
            source: segment.source,
            start: segment.start + from,
            end: segment.start + to,
            reverse: false,
        }
    })
}
/// Project one physically generated segment at its traversal coordinates into another
/// complete assignment. Public validation remains authoritative for the new joins.
fn project_segment(
    target: &routes::Assignment,
    slot: usize,
    segment: &routes::Segment,
    left: u64,
    right: u64,
) -> io::Result<Option<routes::Assignment>> {
    if slot >= target.routes.len()
        || segment.end - segment.start != right - left
        || right > target.routes[slot].length()?
    {
        return Ok(None);
    }
    let mut projected = target.clone();
    let mut pieces = route_slice(&target.routes[slot], 0, left)?;
    pieces.push(segment.clone());
    pieces.extend(route_slice(
        &target.routes[slot],
        right,
        target.routes[slot].length()?,
    )?);
    let mut segments: Vec<routes::Segment> = Vec::with_capacity(pieces.len());
    for piece in pieces {
        if let Some(previous) = segments.last_mut() {
            let contiguous = if piece.reverse {
                previous.reverse && previous.source == piece.source && previous.start == piece.end
            } else {
                !previous.reverse && previous.source == piece.source && previous.end == piece.start
            };
            if contiguous {
                if piece.reverse {
                    previous.start = piece.start;
                } else {
                    previous.end = piece.end;
                }
                continue;
            }
        }
        segments.push(piece);
    }
    projected.routes[slot] = routes::Route { segments };
    Ok((projected != *target).then_some(projected))
}
struct Engine<'a, 'b> {
    e: &'a mut routes::Evaluator<'b>,
    options: Options,
    meter: haploid::Meter,
    histogram: u64,
    bases: Vec<Baseline>,
    tasks: VecDeque<Work>,
    chains: VecDeque<PairedChain>,
    next_chain_key: (usize, usize, usize),
    chain_attempts: u64,
    chain_completed: u64,
    chain_confirmed: u64,
    pair_attempts: u64,
    pair_completed: u64,
    native_completed: u64,
    confirmed: u64,
    conflicts: u64,
    compound: u64,
    best: f64,
    ties: Vec<Pair>,
    ledger: File,
    events: File,
    active: Option<Work>,
    active_chain: Option<PairedChain>,
    active_pair: Option<Pair>,
    chain_result_binding: Option<(usize, usize, Region)>,
    expansion: Option<Scored>,
    unsupported: usize,
    chain_head: bool,
    credits: usize,
    peak_owned: usize,
    shared_emitted: usize,
    receipt: Option<Receipt>,
    emission: Option<Emission>,
    ordinary_turn: usize,
    chain_turn: bool,
    native_cursor: [usize; 2],
    native_done: bool,
    family_count: usize,
    fixed_native: bool,
    blocked: Vec<Value>,
}
impl Engine<'_, '_> {
    /// One production dispatch, or finite retirement; errors preserve the current remainder.
    fn scheduler_step(&mut self) -> io::Result<bool> {
        self.admit_chains()?;
        let chain_ready = !self.chains.is_empty() && self.owned() < self.options.max_tasks;
        let ordinary = if chain_ready && self.chain_turn {
            None
        } else {
            self.ready()?
        };
        if chain_ready && (self.chain_turn || ordinary.is_none()) {
            self.meter.work(1)?;
            self.credits = 1;
            self.shared_emitted = 0;
            self.account("chain_result_reserved_before_advance")?;
            self.chain()?;
            if let Some(p) = self.emission.take().filter(Emission::has_intents) {
                self.transfer(Work::Emission(p));
            }
            self.chain_turn = false;
        } else if let Some((index, grant, class)) = ordinary {
            let work = self.begin(index, grant)?;
            line(
                &mut self.events,
                &json!({"event":"frontier_dispatch","class":class,"input":work,"grant":grant}),
            )?;
            self.dispatch(work)?;
            self.ordinary_turn = (class + 1) % 3;
            self.chain_turn = true;
        } else if self.tasks.is_empty()
            && self.chains.is_empty()
            && self.next_chain_key.0 >= self.bases.len()
        {
            self.chain_head = false;
            self.account("finite_schedule_retired")?;
            return Ok(false);
        } else {
            let owned = self.owned();
            line(
                &mut self.events,
                &json!({"event":"frontier_no_progress","owned":owned,"blocked":self.blocked,"chain_result_requires":usize::from(!self.chains.is_empty())}),
            )?;
            return Err(invalid("budget: frontier_no_progress"));
        }
        self.release_current()?;
        Ok(true)
    }
    fn write_pending(&self) -> io::Result<()> {
        genome::write_json(
            &self.options.out_dir.join("pending.json"),
            &json!({"next_native_pair_cursor":self.native_cursor,"native_complete":self.native_done,"family_count":self.family_count,"active_work_snapshot":self.active,"active_task":self.active.as_ref().and_then(|w| if let Work::Geometry(t)=w {Some(t)} else {None}),"active_snapshot_is_not_additional_work":true,"shared_outputs_accepted":self.shared_emitted,"completion_receipt":self.receipt,"chain_result_binding":self.chain_result_binding,"active_emission":self.emission,"owned128":self.owned(),"reserved_or_remainder_entitlements":self.credits,"chain_admission_head":self.chain_head,"blocked":self.blocked,"queued":self.tasks,"active_chain":self.active_chain,"queued_chains":self.chains,"next_chain_key":self.next_chain_key,"active_pair":self.active_pair,"pending_expansion":self.expansion,"immutable_ordered_baselines":self.bases,"resume_supported":false}),
        )
    }
    fn owned(&self) -> usize {
        usize::from(self.chain_head) + self.tasks.len() + self.credits
    }
    fn account(&mut self, stage: &str) -> io::Result<()> {
        let owned = self.owned();
        ensure(
            owned <= self.options.max_tasks,
            "frontier ownership invariant",
        )?;
        self.peak_owned = self.peak_owned.max(owned);
        line(
            &mut self.events,
            &json!({"event":"frontier_ownership","stage":stage,"owned":owned,"queued":self.tasks.len(),"active_or_reserved":self.credits,"chain_head":self.chain_head,"chain_cursors":self.chains.len()+usize::from(self.active_chain.is_some()),"peak":self.peak_owned}),
        )
    }
    fn transfer(&mut self, work: Work) {
        assert!(self.credits > 0);
        self.credits -= 1;
        self.tasks.push_back(work);
    }
    fn begin(&mut self, index: usize, grant: usize) -> io::Result<Work> {
        let grant = grant.max(1);
        ensure(
            self.owned() + grant - 1 <= self.options.max_tasks,
            "frontier pre-entry reservation invariant",
        )?;
        // VecDeque::remove shifts the shorter side, preserving all waiting ages.
        self.meter
            .work(1 + index.min(self.tasks.len() - index - 1) as u64)?;
        let w = self.tasks.remove(index).unwrap();
        self.credits = grant;
        self.active = Some(w.clone());
        self.shared_emitted = 0;
        self.account("reserved_before_entry")?;
        Ok(w)
    }
    fn release_current(&mut self) -> io::Result<()> {
        self.credits = 0;
        self.shared_emitted = 0;
        self.active = None;
        self.active_pair = None;
        self.chain_result_binding = None;
        self.receipt = None;
        self.emission = None;
        self.account("completed_release")
    }
    fn geometry_grant(&self, t: &PairedTask) -> io::Result<usize> {
        let base = t.task.base();
        ensure(
            base < self.bases.len() && t.copy < 2,
            "immutable paired task binding",
        )?;
        let slots = self.bases[base].scored.pair.copies[t.copy].routes.len();
        Ok(match &t.task {
            Task::Coarse { slot, .. } => {
                if *slot >= slots {
                    1
                } else {
                    2
                }
            }
            Task::Pairs {
                i, j, left, right, ..
            } => {
                if *i >= left.1 || *j >= right.1 {
                    1
                } else {
                    2
                }
            }
            Task::Candidate { other, .. } => {
                if *other < slots {
                    2
                } else {
                    1
                }
            }
            Task::Region { .. } => 1,
        })
    }
    fn emission_plan(&mut self, p: &Emission) -> io::Result<EmissionPlan> {
        let mut count = 0;
        let mut first = None;
        for i in p.advice_index..2 {
            if let Some(a) = &p.advice[i] {
                self.meter.work(4)?;
                let base = p
                    .source_base
                    .ok_or_else(|| invalid("refinement without source baseline"))?;
                let children = a.region.refinements(
                    self.bases[base].scored.pair.copies[a.copy].routes[a.region.slot].length()?,
                    self.options.max_level,
                );
                let start = if i == p.advice_index {
                    p.child_index
                } else {
                    0
                };
                count += children.len().saturating_sub(start);
                if first.is_none() && start < children.len() {
                    first = Some((i, start, base, a.copy, children[start].clone()));
                }
            }
        }
        let mut seed_count = 0;
        let mut seed = None;
        for i in p.seed_index..2 {
            if let Some(base) = p.seeds[i] {
                let copy = if i == p.seed_index { p.seed_copy } else { 0 };
                seed_count += 2 - copy;
                if copy < 2 && seed.is_none() {
                    seed = Some((i, base, copy));
                }
            }
        }
        let remaining = count + seed_count;
        let mut next = p.clone();
        let action = if seed.is_some() && (p.seed_turn || first.is_none()) {
            let (i, base, copy) = seed.unwrap();
            next.seed_index = i + usize::from(copy == 1);
            next.seed_copy = if copy == 1 { 0 } else { 1 };
            next.seed_turn = false;
            Some((
                base,
                copy,
                Region {
                    slot: 0,
                    left: 0,
                    right: 0,
                    level: 0,
                },
                true,
            ))
        } else if let Some((i, j, base, copy, region)) = first {
            next.advice_index = i;
            next.child_index = j + 1;
            next.seed_turn = true;
            Some((base, copy, region, false))
        } else {
            None
        };
        Ok(EmissionPlan {
            next,
            action,
            remaining,
        })
    }
    fn grant(&mut self, work: &Work) -> io::Result<usize> {
        match work {
            Work::Geometry(t) => self.geometry_grant(t),
            Work::Native(n) if n.seed.is_none() => Ok(1),
            Work::Native(n) => {
                let p = self.emission_plan(n.seed.as_ref().unwrap())?;
                Ok(
                    if p.remaining > 1 || (p.remaining == 1 && n.next.is_some()) {
                        2
                    } else {
                        1
                    },
                )
            }
            Work::Emission(e) => Ok(if self.emission_plan(e)?.remaining > 1 {
                2
            } else {
                1
            }),
        }
    }
    fn ready(&mut self) -> io::Result<Option<(usize, usize, usize)>> {
        self.blocked.clear();
        for offset in 0..3 {
            let class = (self.ordinary_turn + offset) % 3;
            for index in 0..self.tasks.len() {
                self.meter.work(1)?;
                if self.tasks[index].class() != class {
                    continue;
                }
                let grant = self.grant(&self.tasks[index].clone())?;
                if self.owned() + grant - 1 <= self.options.max_tasks {
                    return Ok(Some((index, grant, class)));
                }
                // At most128 records: one per work item, only in its own class.
                self.blocked
                    .push(json!({"index":index,"class":class,"additional_required":grant-1}));
            }
        }
        Ok(None)
    }
    fn service_emission(&mut self, p: Emission, native: Option<NativeHead>) -> io::Result<()> {
        let plan = self.emission_plan(&p)?;
        if plan.action.is_some() {
            self.meter.work(1)?;
        }
        // No fallible operation between ordinal/seed-bit mutation and transfers.
        if let Some((base, copy, region, seed)) = &plan.action {
            let task = if *seed {
                assert_eq!(self.bases[*base].seed_owner, p.receipt);
                assert!(!self.bases[*base].seed_emitted[*copy]);
                self.bases[*base].seed_emitted[*copy] = true;
                Task::Coarse {
                    base: *base,
                    slot: 0,
                    level: 0,
                    cell: 0,
                }
            } else {
                Task::Region {
                    base: *base,
                    region: region.clone(),
                    guided: true,
                }
            };
            self.transfer(Work::Geometry(PairedTask { copy: *copy, task }));
        }
        if plan.remaining > 1 {
            self.transfer(match native {
                Some(mut n) => {
                    n.seed = Some(plan.next.clone());
                    Work::Native(n)
                }
                None => Work::Emission(plan.next.clone()),
            });
        } else if let Some(n) = native {
            if let Some(pair) = n.next {
                self.native_cursor = pair;
                self.transfer(Work::Native(NativeHead {
                    families: pair,
                    next: next_native(pair, self.family_count, self.fixed_native),
                    seed: None,
                }));
            } else {
                self.native_done = true;
            }
        }
        line(
            &mut self.events,
            &json!({"event":"frontier_emission","producer":p,"action":plan.action,"remaining_after":plan.remaining.saturating_sub(1),"next":plan.next}),
        )?;
        self.account("emission_transferred")
    }
    fn dispatch(&mut self, work: Work) -> io::Result<()> {
        match work {
            Work::Geometry(t) => {
                self.ordinary(t)?;
                if let Some(p) = self.emission.take().filter(Emission::has_intents) {
                    self.transfer(Work::Emission(p));
                }
            }
            Work::Emission(p) => self.service_emission(p, None)?,
            Work::Native(n) => {
                if let Some(p) = n.seed.clone() {
                    self.service_emission(p, Some(n))?;
                } else {
                    self.native(n.families)?;
                    if let Some(p) = self.emission.take().filter(Emission::has_intents) {
                        let mut n = n;
                        n.seed = Some(p);
                        self.transfer(Work::Native(n));
                    } else if let Some(pair) = n.next {
                        self.native_cursor = pair;
                        self.transfer(Work::Native(NativeHead {
                            families: pair,
                            next: next_native(pair, self.family_count, self.fixed_native),
                            seed: None,
                        }));
                    } else {
                        self.native_done = true;
                    }
                }
            }
        }
        Ok(())
    }
    fn score(&mut self, pair: &Pair, native: bool) -> io::Result<Scored> {
        self.active_pair = Some(pair.clone());
        score::Scorer {
            e: self.e,
            meter: &mut self.meter,
            ledger: &mut self.ledger,
            pair_attempts: &mut self.pair_attempts,
            pair_completed: &mut self.pair_completed,
            pair_cap: self.options.max_scores,
            histogram: self.histogram,
        }
        .evaluate(pair, native)
    }
    fn eligible(&self, pair: &Pair) -> bool {
        !self.bases.iter().any(|b| b.scored.pair.key() == pair.key())
    }
    fn add_base(&mut self, scored: Scored, depth: u8, reason: &str) -> io::Result<()> {
        if let Some(b) = self
            .bases
            .iter_mut()
            .find(|b| b.scored.pair.key() == scored.pair.key())
        {
            b.native_seen |= reason == "native_pair";
            return Ok(());
        }
        self.expansion = Some(scored.clone());
        ensure(
            self.bases.len() < self.options.max_baselines,
            "budget: paired_baselines",
        )?;
        ensure(
            self.credits > 0,
            "baseline requires owned producer entitlement",
        )?;
        let p = self
            .emission
            .as_mut()
            .ok_or_else(|| invalid("baseline without emission owner"))?;
        let slot = p
            .seeds
            .iter()
            .position(Option::is_none)
            .ok_or_else(|| invalid("receipt seed bound"))?;
        let id = self.bases.len();
        p.seeds[slot] = Some(id);
        self.bases.push(Baseline {
            id,
            scored,
            depth,
            batch: 0,
            advice: Vec::with_capacity(3),
            exploratory: Vec::with_capacity(16),
            seed_owner: p.receipt,
            seed_emitted: [false; 2],
            native_seen: reason == "native_pair",
            geometry_services: 0,
        });
        self.expansion = None;
        line(
            &mut self.events,
            &json!({"event":"paired_baseline","id":id,"depth":depth,"reason":reason,"scored":self.bases[id].scored,"slots_are_immutable":true,"seed_owner":self.bases[id].seed_owner,"committed_awaiting_emission":true}),
        )
    }
    fn completed_receipt(&mut self, s: &Scored) {
        self.receipt = Some(Receipt {
            scored: s.clone(),
            stage: "score_returned",
        });
        self.emission = Some(Emission {
            receipt: self.pair_completed,
            ..Emission::default()
        });
    }
    fn remember(&mut self, s: &Scored) -> io::Result<bool> {
        let improved = s.objective < self.best - 1e-9;
        if improved {
            self.best = s.objective;
            self.ties.clear();
        }
        let key = s.pair.key();
        if (s.objective - self.best).abs() <= 1e-9 && !self.ties.contains(&key) {
            ensure(
                self.ties.len() < self.options.max_ties,
                "budget: paired_retained_ties",
            )?;
            self.ties.push(key);
        }
        Ok(improved)
    }
    fn native(&mut self, families: [usize; 2]) -> io::Result<()> {
        let pair = Pair {
            copies: [
                self.e.graph.native_assignment(families[0])?,
                self.e.graph.native_assignment(families[1])?,
            ],
        };
        let s = self.score(&pair, true)?;
        self.native_completed += 1;
        self.completed_receipt(&s);
        self.remember(&s)?;
        self.receipt.as_mut().unwrap().stage = "retention_applied";
        self.add_base(s, 0, "native_pair")?;
        self.receipt.as_mut().unwrap().stage = "commitments_done";
        self.active_pair = None;
        Ok(())
    }
    fn feedback(
        &mut self,
        base: usize,
        copy: usize,
        region: Region,
        s: &Scored,
        improved: bool,
    ) -> io::Result<()> {
        let eligible = self.eligible(&s.pair);
        let b = &mut self.bases[base];
        b.batch += 1;
        b.advice.push(Advice {
            score: s.objective,
            copy,
            region,
        });
        b.advice.sort_by(|a, b| {
            a.score
                .total_cmp(&b.score)
                .then(a.copy.cmp(&b.copy))
                .then(a.region.slot.cmp(&b.region.slot))
                .then(a.region.left.cmp(&b.region.left))
                .then(a.region.right.cmp(&b.region.right))
        });
        b.advice.truncate(2);
        if eligible {
            b.exploratory.push(s.clone());
        }
        let (depth, batch) = (b.depth, b.batch);
        let producer = self
            .emission
            .as_mut()
            .ok_or_else(|| invalid("feedback without owned receipt"))?;
        producer.source_base = Some(base);
        producer.batch = batch;
        if batch % 8 == 0 {
            for (i, a) in std::mem::replace(&mut b.advice, Vec::with_capacity(3))
                .into_iter()
                .enumerate()
            {
                producer.advice[i] = Some(a);
            }
        }
        self.receipt.as_mut().unwrap().stage = "window_folded";
        line(
            &mut self.events,
            &json!({"event":"paired_feedback_frozen","baseline":base,"batch":batch,"producer":self.emission,"window":self.bases[base].exploratory.iter().map(|s| s.objective).collect::<Vec<_>>(),"complete_paired_objective":true}),
        )?;
        if depth < self.options.max_epochs {
            if improved {
                self.add_base(s.clone(), depth + 1, "joint_improvement")?;
            }
            if batch % 16 == 0 {
                self.bases[base].exploratory.sort_by(|a, b| {
                    a.objective
                        .total_cmp(&b.objective)
                        .then(a.pair.key().cmp(&b.pair.key()))
                });
                let selected = self.bases[base]
                    .exploratory
                    .iter()
                    .find(|s| self.eligible(&s.pair))
                    .cloned();
                if let Some(s) = selected {
                    self.add_base(s, depth + 1, "non_improving_physical_exploration")?;
                }
                self.bases[base].exploratory.clear();
            }
        } else if batch % 16 == 0 {
            self.bases[base].exploratory.clear();
        }
        self.receipt.as_mut().unwrap().stage = "commitments_done";
        Ok(())
    }
    fn confirm(&mut self, base: usize, copy: usize, p: continuation::Proposal) -> io::Result<()> {
        ensure(p.base == base && copy < 2, "immutable paired task binding")?;
        let mut pair = self.bases[base].scored.pair.clone();
        pair.copies[copy] = p.candidate;
        if pair == self.bases[base].scored.pair {
            return Ok(());
        }
        self.active_pair = Some(pair.clone());
        if let Some(rejection) =
            oriented::seam_rejection(self.e, &pair.copies[copy], &mut self.meter)?
        {
            line(
                &mut self.ledger,
                &json!({"kind":"oriented_rejection","baseline":base,"copy":copy,"pair":pair,"rejection":rejection}),
            )?;
            self.active_pair = None;
            return Ok(());
        }
        for slot in 0..2 {
            self.meter
                .validation(&pair.copies[slot], self.e.graph.k, self.e.graph.port_count)?;
            match self.e.validate_assignment(&pair.copies[slot]) {
                Ok(a) => pair.copies[slot] = a,
                Err(e)
                    if e.to_string() == "conflicting canonical source-span reuse across genome" =>
                {
                    self.conflicts += 1;
                    line(
                        &mut self.ledger,
                        &json!({"kind":"within_copy_capacity_rejection","baseline":base,"copy":slot,"pair":pair,"compound":p.other>0,"error":e.to_string(),"partner_cannot_repair":true}),
                    )?;
                    self.active_pair = None;
                    return Ok(());
                }
                Err(e) => {
                    line(
                        &mut self.ledger,
                        &json!({"kind":"unexpected_validation_error","pair":pair,"error":e.to_string()}),
                    )?;
                    return Err(e);
                }
            }
        }
        let mut selected = self.score(&pair, false)?;
        self.confirmed += 1;
        if p.other > 0 {
            self.compound += 1;
        }
        let mut joint_scored = 0u64;
        let alternative_count = self.bases[base].exploratory.len();
        if alternative_count > 0 {
            // A bounded physical batch admits at most one fixed baseline window of
            // public-valid cross-copy projections per returned complete candidate.
            // Admission does not require the triggering single-copy score to improve;
            // the existing joint objective chooses among the completed alternatives.
            let other = 1 - copy;
            // Generated complete alternatives contribute only physical breakpoints and
            // source segments. Adjacent breakpoint windows can expose a legal local
            // allele that no single coarse proposal isolated.
            let mut boundaries = pair.copies[other]
                .routes
                .iter()
                .map(|route| Ok(vec![0, route.length()?]))
                .collect::<io::Result<Vec<_>>>()?;
            for alternative_index in 0..alternative_count {
                let alternative = &self.bases[base].exploratory[alternative_index].pair;
                for assignment in &alternative.copies {
                    for (slot, route) in assignment.routes.iter().enumerate() {
                        if slot >= boundaries.len()
                            || route.length()? != pair.copies[other].routes[slot].length()?
                        {
                            continue;
                        }
                        let mut offset = 0;
                        for segment in &route.segments {
                            offset += segment.end - segment.start;
                            boundaries[slot].push(offset);
                        }
                    }
                }
            }
            for cuts in &mut boundaries {
                cuts.sort_unstable();
                cuts.dedup();
            }
            let mut seen = BTreeMap::new();
            seen.insert(selected.pair.key(), ());
            ensure(
                seen.len() == 1 + joint_scored as usize,
                "joint projection key bound",
            )?;
            'alternatives: for alternative_index in 0..alternative_count {
                let alternative = self.bases[base].exploratory[alternative_index].pair.clone();
                for source_copy in 0..2 {
                    for (slot, route) in alternative.copies[source_copy].routes.iter().enumerate() {
                        if slot >= pair.copies[other].routes.len()
                            || route.length()? != pair.copies[other].routes[slot].length()?
                        {
                            continue;
                        }
                        let mut segment_left = 0;
                        for segment in &route.segments {
                            let segment_right = segment_left + segment.end - segment.start;
                            let cuts = boundaries[slot]
                                .iter()
                                .copied()
                                .filter(|cut| *cut >= segment_left && *cut <= segment_right)
                                .collect::<Vec<_>>();
                            for window in cuts.windows(2) {
                                self.meter.work(1)?;
                                let (left, right) = (window[0], window[1]);
                                let physical = segment_slice(
                                    segment,
                                    left - segment_left,
                                    right - segment_left,
                                )?;
                                let Some(projected) = project_segment(
                                    &pair.copies[other],
                                    slot,
                                    &physical,
                                    left,
                                    right,
                                )?
                                else {
                                    continue;
                                };
                                let mut joint = pair.clone();
                                joint.copies[other] = projected;
                                if let Some(rejection) = oriented::seam_rejection(
                                    self.e,
                                    &joint.copies[other],
                                    &mut self.meter,
                                )? {
                                    line(
                                        &mut self.ledger,
                                        &json!({"kind":"joint_projection_oriented_rejection","baseline":base,"copy":other,"pair":joint,"rejection":rejection}),
                                    )?;
                                    continue;
                                }
                                self.meter.validation(
                                    &joint.copies[other],
                                    self.e.graph.k,
                                    self.e.graph.port_count,
                                )?;
                                match self.e.validate_assignment(&joint.copies[other]) {
                                    Ok(validated) => joint.copies[other] = validated,
                                    Err(error) => {
                                        let capacity = error.to_string()
                                            == "conflicting canonical source-span reuse across genome";
                                        self.conflicts += u64::from(capacity);
                                        line(
                                            &mut self.ledger,
                                            &json!({"kind":"joint_projection_validation_rejection","baseline":base,"copy":other,"pair":joint,"capacity":capacity,"error":error.to_string()}),
                                        )?;
                                        continue;
                                    }
                                }
                                let key = joint.key();
                                if !self.eligible(&joint) || seen.contains_key(&key) {
                                    continue;
                                }
                                if self.pair_attempts >= self.options.max_scores {
                                    break 'alternatives;
                                }
                                ensure(
                                    seen.len() == 1 + joint_scored as usize
                                        && joint_scored < JOINT_PROJECTION_WINDOW,
                                    "joint projection key bound",
                                )?;
                                let scored = self.score(&joint, false)?;
                                ensure(
                                    scored.pair.key() == key,
                                    "joint projection normalization changed",
                                )?;
                                self.confirmed += 1;
                                self.compound += 1;
                                joint_scored += 1;
                                ensure(
                                    seen.insert(key, ()).is_none()
                                        && seen.len() == 1 + joint_scored as usize
                                        && joint_scored <= JOINT_PROJECTION_WINDOW,
                                    "joint projection key bound",
                                )?;
                                line(
                                    &mut self.events,
                                    &json!({"event":"paired_joint_projection","baseline":base,"trigger_copy":copy,"projected_copy":other,"source_alternative":alternative_index,"source_copy":source_copy,"slot":slot,"range":[left,right],"scored":scored}),
                                )?;
                                if scored.objective < selected.objective
                                    || (scored.objective - selected.objective).abs() <= 1e-9
                                        && scored.pair.key() < selected.pair.key()
                                {
                                    selected = scored;
                                }
                                if joint_scored == JOINT_PROJECTION_WINDOW {
                                    break 'alternatives;
                                }
                            }
                            segment_left = segment_right;
                        }
                    }
                }
            }
        }
        self.completed_receipt(&selected);
        line(
            &mut self.events,
            &json!({"event":"paired_confirmation","baseline":base,"copy":copy,"region":p.region,"guided":p.guided,"snapped":[p.lo,p.hi],"compound":p.other>0||joint_scored>0,"joint_projection_scores":joint_scored,"scored":selected,"baseline_objective":self.bases[base].scored.objective}),
        )?;
        let improved = self.remember(&selected)?;
        self.receipt.as_mut().unwrap().stage = "retention_applied";
        self.feedback(base, copy, p.region, &selected, improved)?;
        self.active_pair = None;
        Ok(())
    }
    fn ordinary(&mut self, t: PairedTask) -> io::Result<()> {
        let base = t.task.base();
        ensure(
            base < self.bases.len() && t.copy < 2,
            "immutable paired task binding",
        )?;
        line(
            &mut self.events,
            &json!({"event":"paired_geometry_task","baseline":base,"copy":t.copy,"task":t.task}),
        )?;
        let quota = match &t.task {
            Task::Candidate { other, .. } => {
                usize::from(*other < self.bases[base].scored.pair.copies[t.copy].routes.len())
            }
            _ => self.credits,
        };
        self.bases[base].geometry_services += 1;
        let result = continuation::Context {
            assignment: &self.bases[base].scored.pair.copies[t.copy],
            e: self.e,
            meter: &mut self.meter,
            tasks: &mut BoundQueue {
                copy: t.copy,
                tasks: &mut self.tasks,
                credits: &mut self.credits,
                emitted: &mut self.shared_emitted,
                quota,
                cap: self.options.max_tasks,
            },
            events: &mut self.events,
            unsupported_opposite_strand: &mut self.unsupported,
        }
        .step(t.task);
        self.account("shared_return_or_failed_remainder")?;
        let p = result?;
        if let Some(p) = p {
            self.confirm(base, t.copy, p)?;
        }
        Ok(())
    }
    fn admit_chains(&mut self) -> io::Result<()> {
        while self.chains.len() < self.options.max_cursors
            && self.next_chain_key.0 < self.bases.len()
        {
            self.meter.work(1)?;
            let (base, copy, slot) = self.next_chain_key;
            if slot == self.bases[base].scored.pair.copies[copy].routes.len() {
                self.next_chain_key = if copy == 0 {
                    (base, 1, 0)
                } else {
                    (base + 1, 0, 0)
                };
            } else {
                self.chains.push_back(PairedChain {
                    copy,
                    cursor: chain::Cursor::new(base, slot),
                });
                self.next_chain_key.2 += 1;
            }
        }
        Ok(())
    }
    fn chain(&mut self) -> io::Result<()> {
        ensure(
            self.credits == 1,
            "chain requires reserved result entitlement",
        )?;
        let mut c = self.chains.pop_front().unwrap();
        self.active_chain = Some(c.clone());
        ensure(
            self.chain_attempts < self.options.max_chain_primitives,
            "budget: paired_chain_primitives",
        )?;
        self.meter.work(1)?;
        self.chain_attempts += 1;
        let base = c.cursor.base;
        let copy = c.copy;
        let result = c.cursor.advance(
            self.e,
            &self.bases[base].scored.pair.copies[copy],
            &mut self.meter,
        );
        self.active_chain = Some(c.clone());
        let (done, p) = result?;
        self.chain_completed += 1;
        self.active_chain = None;
        if !done {
            self.chains.push_back(c);
        }
        if let Some((region, candidate)) = p {
            self.chain_result_binding = Some((base, copy, region.clone()));
            line(
                &mut self.events,
                &json!({"event":"paired_chain_geometry","baseline":base,"copy":copy,"region":region,"candidate":candidate}),
            )?;
            let (lo, hi) = (region.left, region.right);
            let before = self.confirmed;
            let result = self.confirm(
                base,
                copy,
                continuation::Proposal {
                    base,
                    region,
                    guided: false,
                    lo,
                    hi,
                    other: 0,
                    candidate,
                },
            );
            self.chain_confirmed += self.confirmed - before;
            result?;
        }
        self.active_chain = None;
        Ok(())
    }
}
/// Fixed native family is only a causal-test hook, never a CLI input.
pub fn run(options: Options, fixed_family: Option<usize>) -> io::Result<Value> {
    fs::create_dir(&options.out_dir)?;
    genome::write_json(&options.out_dir.join("options.json"), &options)?;
    let started = Instant::now();
    let result = (|| {
        options.validate()?;
        ensure(
            fixed_family.is_none() || options.max_epochs == 0,
            "causal baseline must disable expansion",
        )?;
        let metadata = fs::metadata(options.routes.join("graph.json"))?.len();
        ensure(
            metadata.saturating_mul(8) < options.max_state_bytes,
            "budget: graph_metadata",
        )?;
        let identity = genome::PanelIdentity::read(&options.panel)?;
        let panel = SyngIndex::load(&options.panel, Default::default())?;
        let graph = routes::Graph::load(&options.routes, &identity)?;
        ensure(
            graph.read_lengths == [150],
            "paired experiment requires L150",
        )?;
        ensure(
            graph.families.iter().all(|f| f.paths.len() <= 4)
                && graph.lanes.iter().all(|l| l.length <= 4096),
            "paired molecule inventory/length cap",
        )?;
        let sample = sample::SampleIndex::load(&options.sample, &identity)?;
        let histogram = *sample
            .stats
            .read_lengths
            .get(&150)
            .ok_or_else(|| invalid("missing L150 histogram"))?;
        let mut meter = haploid::Meter::new(&options.meter_options());
        meter.accounting = Some(File::create(options.out_dir.join("accounting.jsonl"))?);
        // Bounded numerical/public/cache scratch plus all assignment/advice/cursor slots.
        // Logical reservations, not native allocation/RSS guarantees.
        let pair_bytes = 2 * 4 * options.max_segments as u64 * 64 + 4096;
        let joint_projection_key_bytes = pair_bytes * (JOINT_PROJECTION_WINDOW + 1);
        meter.reserve(
            metadata * 8
                + (options.max_features as u64 + 1) * 160 * 12
                + options.max_tasks as u64 * 8192
                + options.max_ties as u64 * pair_bytes
                + options.max_baselines as u64 * (pair_bytes * 18 + 8192)
                + (options.max_cursors as u64 + 2) * (8192 + graph.k * 8)
                + joint_projection_key_bytes,
        )?;
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
        let mut engine = Engine {
            e: &mut evaluator,
            options: options.clone(),
            meter,
            histogram,
            bases: Vec::with_capacity(options.max_baselines),
            tasks: VecDeque::with_capacity(options.max_tasks),
            chains: VecDeque::with_capacity(options.max_cursors),
            next_chain_key: (0, 0, 0),
            chain_attempts: 0,
            chain_completed: 0,
            chain_confirmed: 0,
            pair_attempts: 0,
            pair_completed: 0,
            native_completed: 0,
            confirmed: 0,
            conflicts: 0,
            compound: 0,
            best: f64::INFINITY,
            ties: Vec::with_capacity(options.max_ties),
            ledger: File::create(options.out_dir.join("scores.jsonl"))?,
            events: File::create(options.out_dir.join("events.jsonl"))?,
            active: None,
            active_chain: None,
            active_pair: None,
            chain_result_binding: None,
            expansion: None,
            unsupported: 0,
            chain_head: true,
            credits: 0,
            peak_owned: 2,
            shared_emitted: 0,
            receipt: None,
            emission: None,
            ordinary_turn: 0,
            chain_turn: false,
            native_cursor: fixed_family.map(|f| [f, f]).unwrap_or([0, 0]),
            native_done: false,
            family_count: graph.families.len(),
            fixed_native: fixed_family.is_some(),
            blocked: Vec::with_capacity(options.max_tasks),
        };
        genome::write_json(
            &options.out_dir.join("provenance.json"),
            &json!({"graph":graph.digest()?,"panel":identity,"sample":genome::reconstruction::fingerprint(&options.sample)?,"count_policy":graph.count_policy,"histogram":histogram,"nominal_per_copy_depth":10,"background_once":0.1,"read_length":150,"actual_reads":sample.stats.reads,"actual_bases":sample.stats.bases,"fixed_family_test_only":fixed_family}),
        )?;
        ensure(options.max_tasks >= 2, "budget: frontier_initial_heads")?;
        ensure(
            std::mem::size_of::<Work>() + std::mem::size_of::<Value>() < 8192
                && 4 * std::mem::size_of::<Region>() <= 256,
            "compact frontier layout exceeds reservation",
        )?;
        engine.tasks.push_back(Work::Native(NativeHead {
            families: engine.native_cursor,
            next: next_native(
                engine.native_cursor,
                engine.family_count,
                engine.fixed_native,
            ),
            seed: None,
        }));
        engine.account("initial_heads")?;
        let mut stop = "declared_schedule_finished_not_global".to_string();
        let mut failure = None;
        let execution = (|| -> io::Result<()> {
            while engine.scheduler_step()? {}
            Ok(())
        })();
        if let Err(e) = execution {
            failure = Some(e);
        }
        if let Some(e) = &failure {
            stop = e.to_string();
        }
        engine.ledger.flush()?;
        engine.events.flush()?;
        engine.meter.accounting.as_mut().unwrap().flush()?;
        engine.write_pending()?;
        let result = json!({"status":stop,"native_initialization_complete":engine.native_done && engine.bases.iter().filter(|b| b.native_seen).all(|b| b.seed_emitted == [true;2]),"native_pairs_completed":engine.native_completed,"paired_score_admissions":engine.pair_attempts,"paired_objectives_completed":engine.pair_completed,"public_haploid_losses_completed":engine.meter.public_scores,"confirmed_candidates":engine.confirmed,"compound_confirmations":engine.compound,"within_copy_capacity_rejections":engine.conflicts,"chain_admissions":engine.chain_attempts,"chain_completions":engine.chain_completed,"chain_confirmations":engine.chain_confirmed,"best_paired_objective":engine.best.is_finite().then_some(engine.best),"retained_correlated_best_found":engine.ties,"baseline_count":engine.bases.len(),"frontier_peak_owned":engine.peak_owned,"frontier_owned":engine.owned(),"frontier_layout_bytes":std::mem::size_of::<Work>(),"meter":engine.meter,"elapsed_seconds":started.elapsed().as_secs_f64(),"support_complete":false,"global_bound":null,"sequence_emission_authorized":false,"full_gate_b_complete":false});
        genome::write_json(&options.out_dir.join("result.json"), &result)?;
        if let Some(e) = failure {
            if !e.to_string().starts_with("budget:") {
                return Err(e);
            }
        }
        Ok(result)
    })();
    genome::write_json(
        &options.out_dir.join("status.json"),
        &match &result {
            Ok(v) => json!({"status":"bounded_stop","stop":v["status"]}),
            Err(e) => json!({"status":"failed","error":e.to_string(),"support_complete":false}),
        },
    )?;
    result
}

/// Test-only score-controlled mechanism: real public-valid physical pairs, no CLI input.
#[cfg(test)]
pub fn test_intermediate_improvement_window(
    e: &mut routes::Evaluator<'_>,
    options: Options,
) -> io::Result<bool> {
    let first = e.validate_assignment(&e.graph.native_assignment(0)?)?;
    let second = e.validate_assignment(&e.graph.native_assignment(1)?)?;
    let initial = Pair {
        copies: [first.clone(), first.clone()],
    };
    let candidate = Pair {
        copies: [first, second],
    };
    ensure(
        initial.key() != candidate.key(),
        "test requires distinct physical pair keys",
    )?;
    let mut engine = test_engine(e, options, initial, 15)?;
    let original = engine.bases[0].scored.pair.clone();
    for (base, copy) in [(99, 0), (0, 2)] {
        let error = engine
            .ordinary(PairedTask {
                copy,
                task: Task::Coarse {
                    base,
                    slot: 0,
                    level: 0,
                    cell: 0,
                },
            })
            .err()
            .unwrap();
        ensure(
            error.to_string() == "immutable paired task binding",
            "stale task guard failed",
        )?;
    }
    // 5 improves its baseline10, but not the global best0; global-admission is false.
    let scored = Scored {
        pair: candidate.clone(),
        objective: 5.0,
    };
    engine.completed_receipt(&scored);
    engine.feedback(
        0,
        1,
        Region {
            slot: 0,
            left: 200,
            right: 400,
            level: 8,
        },
        &scored,
        false,
    )?;
    ensure(
        engine.bases[0].id == 0 && engine.bases[0].scored.pair == original,
        "baseline mutated by expansion",
    )?;
    ensure(
        engine.tasks.is_empty() && engine.emission.as_ref().unwrap().seeds[0] == Some(1),
        "new tasks retargeted old baseline",
    )?;
    Ok(engine
        .bases
        .iter()
        .any(|b| b.scored.pair.key() == candidate.key()))
}

#[cfg(test)]
fn test_engine<'a, 'b>(
    e: &'a mut routes::Evaluator<'b>,
    options: Options,
    initial: Pair,
    batch: usize,
) -> io::Result<Engine<'a, 'b>> {
    Ok(Engine {
        e,
        options: options.clone(),
        meter: haploid::Meter::new(&options.meter_options()),
        histogram: 1,
        bases: vec![Baseline {
            id: 0,
            scored: Scored {
                pair: initial,
                objective: 10.0,
            },
            depth: 0,
            batch,
            seed_owner: 0,
            seed_emitted: [true; 2],
            native_seen: true,
            geometry_services: 0,
            advice: Vec::with_capacity(3),
            exploratory: Vec::with_capacity(16),
        }],
        tasks: VecDeque::new(),
        chains: VecDeque::new(),
        next_chain_key: (0, 0, 0),
        chain_attempts: 0,
        chain_completed: 0,
        chain_confirmed: 0,
        pair_attempts: 0,
        pair_completed: 0,
        native_completed: 0,
        confirmed: 0,
        conflicts: 0,
        compound: 0,
        best: 0.0,
        ties: Vec::new(),
        ledger: File::create(options.out_dir.join("mechanism-scores.jsonl"))?,
        events: File::create(options.out_dir.join("mechanism-events.jsonl"))?,
        active: None,
        active_chain: None,
        active_pair: None,
        chain_result_binding: None,
        expansion: None,
        unsupported: 0,
        chain_head: true,
        credits: 1,
        peak_owned: 2,
        shared_emitted: 0,
        receipt: None,
        emission: None,
        ordinary_turn: 0,
        chain_turn: false,
        native_cursor: [0, 0],
        native_done: true,
        family_count: 2,
        fixed_native: true,
        blocked: Vec::with_capacity(options.max_tasks),
    })
}

/// Controlled transition fixtures only: never an input path to run/CLI.
#[cfg(test)]
pub fn test_frontier_transitions(
    e: &mut routes::Evaluator<'_>,
    options: Options,
) -> io::Result<()> {
    let a = e.validate_assignment(&e.graph.native_assignment(0)?)?;
    let b = e.validate_assignment(&e.graph.native_assignment(1)?)?;
    let pair = Pair {
        copies: [a.clone(), a.clone()],
    };
    let candidate = Pair { copies: [a, b] };
    for label in [
        "full-drain",
        "partial-shared",
        "all-growing",
        "windows",
        "two-commit-stop",
        "chain-block",
        "fairness",
    ] {
        let mut o = options.clone();
        o.out_dir = options.out_dir.join(label);
        fs::create_dir(&o.out_dir)?;
        if label == "two-commit-stop" {
            o.max_baselines = 2;
        }
        let mut h = test_engine(e, o, pair.clone(), 15)?;
        h.credits = 0;
        let coarse = Work::Geometry(PairedTask {
            copy: 0,
            task: Task::Coarse {
                base: 0,
                slot: 0,
                level: 0,
                cell: 0,
            },
        });
        match label {
            "full-drain" | "partial-shared" | "all-growing" => {
                for _ in 0..127 {
                    h.tasks.push_back(coarse.clone());
                }
                assert_eq!(h.owned(), 128);
                assert!(h.ready()?.is_none());
                assert_eq!(h.pair_attempts, 0);
                if label != "all-growing" {
                    h.tasks[126] = Work::Geometry(PairedTask {
                        copy: 1,
                        task: Task::Pairs {
                            base: 0,
                            region: Region {
                                slot: 0,
                                left: 0,
                                right: 10,
                                level: 0,
                            },
                            guided: false,
                            lo: 0,
                            hi: 10,
                            left: (0, 1),
                            right: (0, 1),
                            i: 1,
                            j: 0,
                        },
                    });
                    let (i, g, _) = h.ready()?.unwrap();
                    assert_eq!((i, g), (126, 1));
                    let work = h.begin(i, g)?;
                    h.dispatch(work)?;
                    h.release_current()?;
                    assert_eq!(h.owned(), 127);
                    assert_eq!(h.pair_completed, 0);
                    if label == "partial-shared" {
                        let work = h.begin(0, 2)?;
                        h.meter.limits.max_work = h.meter.work_used + 2; // step + first push only
                        let err = h.dispatch(work).unwrap_err();
                        assert_eq!(err.to_string(), "budget: work");
                        assert_eq!(h.shared_emitted, 1);
                        assert_eq!(h.credits, 1);
                        assert_eq!(h.owned(), 128);
                        assert!(matches!(
                            h.tasks.back(),
                            Some(Work::Geometry(PairedTask {
                                task: Task::Region { .. },
                                ..
                            }))
                        ));
                    }
                }
            }
            "windows" | "two-commit-stop" => {
                h.credits = 1;
                let first = Scored {
                    pair: candidate.clone(),
                    objective: 5.0,
                };
                if label == "two-commit-stop" {
                    h.bases[0].exploratory.push(first.clone());
                    let second = Scored {
                        pair: Pair {
                            copies: [candidate.copies[1].clone(), candidate.copies[1].clone()],
                        },
                        objective: 1.0,
                    };
                    h.completed_receipt(&second);
                    let err = h
                        .feedback(
                            0,
                            1,
                            Region {
                                slot: 0,
                                left: 0,
                                right: 384,
                                level: 3,
                            },
                            &second,
                            true,
                        )
                        .unwrap_err();
                    assert_eq!(err.to_string(), "budget: paired_baselines");
                    assert_eq!(h.bases.len(), 2);
                    assert_eq!(h.emission.as_ref().unwrap().seeds, [Some(1), None]);
                    assert!(h.expansion.is_some());
                    assert_eq!(h.bases[1].seed_emitted, [false; 2]);
                } else {
                    h.completed_receipt(&first);
                    h.feedback(
                        0,
                        1,
                        Region {
                            slot: 0,
                            left: 0,
                            right: 384,
                            level: 3,
                        },
                        &first,
                        false,
                    )?;
                    let frozen = h.emission.clone().unwrap();
                    assert_eq!(h.bases[1].seed_owner, frozen.receipt);
                    let p = h.emission_plan(&frozen)?;
                    assert_eq!(p.remaining, 6); // four children + two roots
                    assert_eq!(p.action.as_ref().unwrap().1, 1);
                    let second = h.emission_plan(&p.next)?;
                    assert!(second.action.as_ref().unwrap().3);
                    assert_eq!(second.action.as_ref().unwrap().1, 0);
                    h.bases[0].batch = 23;
                    h.completed_receipt(&first);
                    h.feedback(
                        0,
                        0,
                        Region {
                            slot: 1,
                            left: 10,
                            right: 30,
                            level: 2,
                        },
                        &first,
                        false,
                    )?;
                    assert_eq!(frozen.batch, 16);
                    assert_eq!(frozen.advice[0].as_ref().unwrap().copy, 1);
                    assert_eq!(h.emission.as_ref().unwrap().batch, 24);
                    assert_eq!(h.bases.len(), 2); // duplicate physical key cannot mint another seed owner
                    h.emission = None;
                    let mut cursor = frozen;
                    let mut emissions = 0;
                    loop {
                        let remaining = h.emission_plan(&cursor)?.remaining;
                        h.credits = if remaining > 1 { 2 } else { 1 };
                        h.service_emission(cursor, None)?;
                        emissions += 1;
                        if remaining <= 1 {
                            break;
                        }
                        cursor = match h.tasks.pop_back().unwrap() {
                            Work::Emission(p) => p,
                            _ => panic!("missing producer"),
                        };
                    }
                    assert_eq!(emissions, 6);
                    assert_eq!(h.bases[1].seed_emitted, [true; 2]);
                    assert_eq!(h.tasks.len(), 6);
                    assert_eq!(h.pair_completed, 0); // controlled scores, no evaluator/objective invocation
                }
            }
            "chain-block" => {
                for _ in 0..127 {
                    h.tasks.push_back(coarse.clone());
                }
                h.admit_chains()?;
                let before = serde_json::to_value(&h.chains)?;
                assert_eq!(h.owned(), 128);
                assert!(!h.chains.is_empty());
                assert!(!(h.owned() < h.options.max_tasks));
                assert_eq!(before, serde_json::to_value(&h.chains)?);
                assert_eq!((h.chain_attempts, h.chain_completed), (0, 0));
            }
            "fairness" => {
                h.tasks.push_back(Work::Native(NativeHead {
                    families: [0, 0],
                    next: None,
                    seed: None,
                }));
                h.tasks.push_back(coarse);
                h.tasks.push_back(Work::Emission(Emission::default()));
                for class in 0..3 {
                    h.ordinary_turn = class;
                    assert_eq!(h.ready()?.unwrap().2, class);
                }
                // Final seed needs one slot only if there is no future native pair.
                let p = Emission {
                    seeds: [Some(0), None],
                    seed_copy: 1,
                    ..Emission::default()
                };
                assert_eq!(
                    h.grant(&Work::Native(NativeHead {
                        families: [0, 0],
                        next: None,
                        seed: Some(p.clone())
                    }))?,
                    1
                );
                assert_eq!(
                    h.grant(&Work::Native(NativeHead {
                        families: [0, 0],
                        next: Some([0, 1]),
                        seed: Some(p)
                    }))?,
                    2
                );
            }
            _ => unreachable!(),
        }
        genome::write_json(
            &h.options.out_dir.join("mechanical-result.json"),
            &json!({"case":label,"owned":h.owned(),"credits":h.credits,"queued":h.tasks,"active_snapshot":h.active,"emission":h.emission,"receipt":h.receipt,"baselines":h.bases,"expansion":h.expansion,"meter":h.meter,"controlled_transition_fixture_not_scientific_run":true}),
        )?;
    }
    Ok(())
}

/// Exercises the same dispatch and pending writer as run; controlled states are not CLI inputs.
#[cfg(test)]
pub fn test_scheduler_transitions(
    e: &mut routes::Evaluator<'_>,
    options: Options,
    diagnostics: bool,
) -> io::Result<()> {
    let a = e.validate_assignment(&e.graph.native_assignment(0)?)?;
    let b = e.validate_assignment(&e.graph.native_assignment(1)?)?;
    let pair = Pair {
        copies: [a.clone(), a.clone()],
    };
    let coarse = Work::Geometry(PairedTask {
        copy: 0,
        task: Task::Coarse {
            base: 0,
            slot: 0,
            level: 0,
            cell: 0,
        },
    });
    let terminal = Work::Geometry(PairedTask {
        copy: 1,
        task: Task::Pairs {
            base: 0,
            region: Region {
                slot: 0,
                left: 0,
                right: 10,
                level: 0,
            },
            guided: false,
            lo: 0,
            hi: 10,
            left: (0, 1),
            right: (0, 1),
            i: 1,
            j: 0,
        },
    });
    let labels: &[&str] = if diagnostics {
        &[
            "inactive",
            "successful-retirement",
            "partial-shared-failure",
            "chain-failure",
        ]
    } else {
        &[
            "recurring-turns",
            "full-drain-seed",
            "blocked-producer",
            "all-growing-stop",
            "dormant-retirement",
            "chain-32",
        ]
    };
    let mut diagnostic_failures = Vec::new();
    for &label in labels {
        let mut o = options.clone();
        o.out_dir = options.out_dir.join(label);
        fs::create_dir(&o.out_dir)?;
        let mut h = test_engine(e, o, pair.clone(), 15)?;
        h.credits = 0;
        h.next_chain_key = (1, 0, 0);
        let mut expected_outputs = 0;
        match label {
            "inactive" => {}
            "successful-retirement" => {
                h.tasks.push_back(coarse.clone());
                assert!(h.scheduler_step()?);
                assert_eq!(h.tasks.len(), 2); // both actual shared pushes accepted
                assert!(h.active.is_none());
                assert_eq!(h.credits, 0);
            }
            "partial-shared-failure" => {
                h.tasks.push_back(coarse.clone());
                // Two readiness visits, removal, shared entry and first push; second push fails.
                h.meter.limits.max_work = 5;
                assert_eq!(h.scheduler_step().unwrap_err().to_string(), "budget: work");
                assert_eq!(h.tasks.len(), 1);
                assert_eq!(h.credits, 1);
                assert!(matches!(h.active, Some(Work::Geometry(_))));
                assert!(matches!(
                    h.tasks.back(),
                    Some(Work::Geometry(PairedTask {
                        task: Task::Region { .. },
                        ..
                    }))
                ));
                expected_outputs = 1;
            }
            "chain-failure" => {
                // Independently test initialization, not just the successful-release reset.
                h.shared_emitted = 9;
                h.chains.push_back(PairedChain {
                    copy: 1,
                    cursor: chain::Cursor::new(0, 0),
                });
                h.options.max_chain_primitives = 0;
                assert_eq!(
                    h.scheduler_step().unwrap_err().to_string(),
                    "budget: paired_chain_primitives"
                );
                assert!(h.active.is_none());
                assert!(h.active_chain.is_some());
                assert_eq!(h.credits, 1);
            }
            "recurring-turns" => {
                h.next_chain_key = (0, 0, 0);
                h.native_done = false;
                h.fixed_native = false;
                // Native scoring heads are class 0; a native head owing roots is class 2.
                h.tasks.push_back(Work::Native(NativeHead {
                    families: [0, 0],
                    next: Some([0, 1]),
                    seed: None,
                }));
                for _ in 0..3 {
                    h.tasks.push_back(terminal.clone());
                }
                for _ in 0..2 {
                    h.tasks.push_back(Work::Emission(Emission::default()));
                }
                let mut turns = Vec::new();
                for step in 0..12 {
                    let before = h.chain_attempts;
                    assert!(h.scheduler_step()?);
                    let chain = h.chain_attempts != before;
                    assert_eq!(chain, step % 2 == 1, "outer turn {step}");
                    if !chain {
                        turns.push((h.ordinary_turn + 2) % 3);
                    }
                    assert_eq!(h.chain_turn, !chain);
                }
                assert_eq!(turns, [0, 1, 2, 0, 1, 2]);
                assert_eq!(h.bases[0].seed_emitted, [true; 2]);
                assert!(!h.native_done); // geometry and chains ran without a native barrier
                assert_eq!(h.pair_completed, 2);
                assert_eq!(h.meter.public_scores, 4);
            }
            "full-drain-seed" | "blocked-producer" => {
                let scored = Scored {
                    pair: Pair {
                        copies: [a.clone(), b.clone()],
                    },
                    objective: 5.0,
                };
                h.credits = 1;
                h.completed_receipt(&scored);
                h.feedback(
                    0,
                    1,
                    Region {
                        slot: 0,
                        left: 0,
                        right: 384,
                        level: 3,
                    },
                    &scored,
                    false,
                )?;
                let frozen = h.emission.take().unwrap();
                let owner = frozen.receipt;
                assert_eq!(frozen.batch, 16);
                assert_eq!(h.bases[1].seed_owner, owner);
                assert_eq!(h.bases[1].seed_emitted, [false; 2]);
                h.transfer(Work::Emission(frozen.clone()));
                h.release_current()?;
                for _ in 0..124 {
                    h.tasks.push_back(coarse.clone());
                }
                h.tasks.push_back(terminal.clone());
                h.ordinary_turn = 2; // initial state only; all subsequent turns are production updates
                h.next_chain_key = (0, 0, 0);
                assert_eq!(h.owned(), 127);
                assert!(h.scheduler_step()?); // first refinement fills the pool
                assert_eq!(h.owned(), 128);
                assert_eq!(h.bases[1].seed_emitted, [false; 2]);
                // The extra fixture starts here: select a capacity-blocked producer at128.
                // After this initial-state choice, only production dispatch updates the turns.
                let blocked_producer = if label == "blocked-producer" {
                    h.ordinary_turn = 2;
                    let index = h
                        .tasks
                        .iter()
                        .position(|w| matches!(w, Work::Emission(p) if p.receipt == owner))
                        .unwrap();
                    Some((index, serde_json::to_value(&h.tasks[index])?))
                } else {
                    None
                };
                let before = h.chain_attempts;
                assert!(h.scheduler_step()?); // skips blocked work; terminal geometry drains
                assert_eq!(h.owned(), 127);
                assert_eq!(h.chain_attempts, before);
                if let Some((index, producer)) = blocked_producer {
                    assert_eq!(
                        h.blocked.first(),
                        Some(&json!({
                            "index":index,"class":2,"additional_required":1
                        }))
                    );
                    let after = h
                        .tasks
                        .iter()
                        .find(|w| matches!(w, Work::Emission(p) if p.receipt == owner))
                        .unwrap();
                    assert_eq!(serde_json::to_value(after)?, producer);
                    assert_eq!(h.ordinary_turn, 2); // terminal class1 selected by production
                    assert_eq!(h.bases[1].seed_emitted, [false; 2]);
                    assert_eq!(h.bases[1].seed_owner, owner);
                    genome::write_json(
                        &h.options.out_dir.join("blocked-producer-proof.json"),
                        &json!({"producer_before":producer,"producer_after":after,
                            "blocked":h.blocked,"owned_before":128,"owned_after":h.owned(),
                            "ordinary_turn_after":h.ordinary_turn,"seed_owner":owner,
                            "seed_emitted":h.bases[1].seed_emitted}),
                    )?;
                }
                assert_eq!(h.bases[0].geometry_services, 1);
                assert!(h.scheduler_step()?); // outer chain share once a result credit is available
                assert_eq!(h.chain_attempts, before + 1);
                assert!(h.scheduler_step()?); // seed arm, before remaining refinement siblings
                assert_eq!(h.owned(), 128);
                assert_eq!(h.bases[1].seed_emitted, [true, false]);
                assert_eq!(h.bases[1].seed_owner, owner);
                let p = match h.tasks.back().unwrap() {
                    Work::Emission(p) => p,
                    _ => panic!("producer missing"),
                };
                assert_eq!(p.batch, 16);
                assert_eq!(
                    serde_json::to_value(&p.advice)?,
                    serde_json::to_value(&frozen.advice)?
                );
                assert_eq!(h.pair_completed, 0);
            }
            "all-growing-stop" | "chain-32" => {
                for _ in 0..127 {
                    h.tasks.push_back(coarse.clone());
                }
                if label == "chain-32" {
                    h.next_chain_key = (0, 0, 0);
                    for i in 0..32 {
                        h.chains.push_back(PairedChain {
                            copy: i % 2,
                            cursor: chain::Cursor::new(0, 0),
                        });
                    }
                }
                let before = serde_json::to_value(&h.chains)?;
                let key = h.next_chain_key;
                assert_eq!(
                    h.scheduler_step().unwrap_err().to_string(),
                    "budget: frontier_no_progress"
                );
                assert_eq!(h.owned(), 128);
                assert_eq!(h.blocked.len(), 127);
                assert_eq!(h.next_chain_key, key);
                assert_eq!(serde_json::to_value(&h.chains)?, before);
                assert_eq!(
                    (
                        h.chain_attempts,
                        h.chain_completed,
                        h.pair_completed,
                        h.meter.public_scores
                    ),
                    (0, 0, 0, 0)
                );
                assert!(h.active.is_none());
                assert_eq!(h.credits, 0);
                if label == "chain-32" {
                    h.write_pending()?;
                    fs::copy(
                        h.options.out_dir.join("pending.json"),
                        h.options.out_dir.join("denied-pending.json"),
                    )?;
                    h.tasks[126] = terminal.clone();
                    assert!(h.scheduler_step()?); // actual drain, not snapshot-only backpressure
                    assert_eq!(h.owned(), 127);
                    h.options.max_chain_primitives = 0;
                    assert_eq!(
                        h.scheduler_step().unwrap_err().to_string(),
                        "budget: paired_chain_primitives"
                    );
                    assert_eq!(h.chains.len(), 31);
                    assert!(h.active_chain.is_some());
                    assert_eq!(
                        h.chains.len() + usize::from(h.active_chain.is_some()),
                        h.options.max_cursors
                    );
                    assert_eq!(h.next_chain_key, key);
                }
            }
            "dormant-retirement" => {
                assert!(!h.scheduler_step()?);
                assert!(!h.chain_head);
                assert_eq!(h.owned(), 0);
                assert_eq!(h.pair_completed, 0);
            }
            _ => unreachable!(),
        }
        h.write_pending()?;
        let pending: Value =
            serde_json::from_reader(File::open(h.options.out_dir.join("pending.json"))?)?;
        assert_eq!(pending["owned128"], h.owned());
        assert_eq!(pending["reserved_or_remainder_entitlements"], h.credits);
        assert_eq!(pending["resume_supported"], false);
        if diagnostics && pending["shared_outputs_accepted"] != expected_outputs {
            diagnostic_failures.push(format!(
                "{label}: expected {expected_outputs}, got {}",
                pending["shared_outputs_accepted"]
            ));
        }
        genome::write_json(
            &h.options.out_dir.join("scheduler-result.json"),
            &json!({"case":label,
            "owned":h.owned(),"chain_attempts":h.chain_attempts,"chain_completed":h.chain_completed,
            "paired_objectives":h.pair_completed,"public_haploid_returns":h.meter.public_scores,
            "controlled_fixture_not_scientific_run":true}),
        )?;
    }
    assert!(
        diagnostic_failures.is_empty(),
        "pending diagnostics: {}",
        diagnostic_failures.join("; ")
    );
    Ok(())
}
