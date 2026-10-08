//! Command-owned fair sample-ranked lazy coupled search. No backend identity inputs.
mod adapter;
mod checkpoint;
mod continuation;
mod machine;
mod source_order;
mod transition;
mod v1_schema;
mod validation;
mod validation_containers;
use crate::genome_inference::{self as genome, panel_routes as routes};
use crate::sample_mem_bwt::invalid;
use adapter::{Permutation, Port, SourcePermutation};
use routes::{Assignment, Evaluator, Route, Segment};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};
use std::path::{Path, PathBuf};

pub const POLICY_VERSION: &str = "fair-sample-ranked-source-pair-order-v3";
fn ensure(ok: bool, message: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(message))
    }
}
fn plus(n: u64, x: u64) -> io::Result<u64> {
    n.checked_add(x)
        .ok_or_else(|| invalid("policy counter overflow"))
}
fn bytes(value: &impl Serialize) -> io::Result<u64> {
    Ok(serde_json::to_vec(value).map_err(io::Error::other)?.len() as u64)
}
pub fn identity() -> BTreeMap<String, String> {
    [
        ("version", POLICY_VERSION),
        ("policy", include_str!("mod.rs")),
        ("machine", include_str!("machine.rs")),
        ("adapter", include_str!("adapter.rs")),
        ("checkpoint", include_str!("checkpoint.rs")),
        ("cli", include_str!("../../genome_infer.rs")),
        ("continuation", include_str!("continuation.rs")),
        ("transition", include_str!("transition.rs")),
        ("source_order", include_str!("source_order.rs")),
        ("v1_schema", include_str!("v1_schema.rs")),
        ("validation", include_str!("validation.rs")),
        (
            "validation_containers",
            include_str!("validation_containers.rs"),
        ),
        ("supported_v1_policy", include_str!("supported_v1/mod.rs")),
        (
            "supported_v1_machine",
            include_str!("supported_v1/machine.rs"),
        ),
        (
            "supported_v1_adapter",
            include_str!("supported_v1/adapter.rs"),
        ),
        (
            "supported_v1_checkpoint",
            include_str!("supported_v1/checkpoint.rs"),
        ),
        (
            "supported_v1_cli",
            include_str!("supported_v1/genome_infer.rs"),
        ),
    ]
    .into_iter()
    .map(|(name, source)| {
        (
            name.into(),
            if name == "version" {
                source.into()
            } else {
                format!(
                    "fnv1a64-source-{:016x}",
                    genome::checksum(source.as_bytes())
                )
            },
        )
    })
    .collect()
}
#[derive(Clone, Debug, clap::Args)]
pub struct Options {
    #[arg(long, default_value_t = 100000)]
    pub max_work: u64,
    #[arg(long, default_value_t = 10000)]
    pub max_evaluations: u64,
    /// Retained logical occupancy, NOT RSS; see policy documentation
    #[arg(long, default_value_t = 268435456)]
    pub max_state_bytes: u64,
    #[arg(long, default_value_t = 1000)]
    pub max_optima: usize,
    #[arg(long, default_value_t = 1e-9)]
    pub tie_epsilon: f64,
    /// Immutable previous output directory; current --out-dir must be NEW
    #[arg(long)]
    pub resume_from: Option<PathBuf>,
    /// Explicit conversion only; no search or rescoring. Requires an exact supported v1 parent.
    #[arg(long, requires_all = ["resume_from", "v1_ancestry_manifest"])]
    pub convert_exact_v1_to_v2: bool,
    /// Independently frozen exact old checkpoint/ledger pins, including every ancestor.
    #[arg(long, requires = "convert_exact_v1_to_v2")]
    pub v1_ancestry_manifest: Option<PathBuf>,
    /// Update a compatible v2 checkpoint to paired source order only; no search or rescoring.
    #[arg(
        long,
        requires = "resume_from",
        conflicts_with = "convert_exact_v1_to_v2"
    )]
    pub update_source_pair_order: bool,
    /// Explicit authorization to increase cumulative work/evaluation/storage/support caps
    #[arg(long, requires = "resume_from")]
    pub extend_budgets: bool,
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
struct Limits {
    max_work: u64,
    max_evaluations: u64,
    max_state_bytes: u64,
    max_optima: usize,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Scored {
    assignment: Assignment,
    objective_bits: u64,
}
impl Scored {
    fn objective(&self) -> f64 {
        f64::from_bits(self.objective_bits)
    }
    fn json(&self) -> serde_json::Value {
        serde_json::json!({"assignment":self.assignment,"relative_objective":self.objective()})
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Context {
    family: usize,
    slot: usize,
    completed: Vec<Route>,
    segments: Vec<Segment>,
    source: usize,
    cut: u64,
    reverse: bool,
}
impl Context {
    fn spans(&self) -> impl Iterator<Item = &Segment> {
        self.completed
            .iter()
            .flat_map(|r| &r.segments)
            .chain(&self.segments)
    }
    fn piece(&self, to: u64) -> Segment {
        Segment {
            source: self.source,
            start: self.cut.min(to),
            end: self.cut.max(to),
            reverse: self.reverse,
        }
    }
    fn depth(&self) -> usize {
        self.segments.len()
            + self
                .completed
                .iter()
                .map(|r| r.segments.len().saturating_sub(1))
                .sum::<usize>()
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum AfterCheck {
    Close,
    Hub(Port),
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum Op {
    Start,
    SourceBound {
        lo: u64,
        hi: u64,
    },
    SourceScan {
        base: u64,
        permutation: SourcePermutation,
    },
    Check {
        piece: Segment,
        index: usize,
        next: AfterCheck,
    },
    Close {
        piece: Segment,
    },
    HubBound {
        port: Port,
        lower: Option<u64>,
        lo: u64,
        hi: u64,
    },
    HubScan {
        port: Port,
        base: u64,
        permutation: Permutation,
    },
    Child {
        piece: Segment,
        donor: Port,
    },
    Probe {
        routes: Vec<Route>,
        next_slot: usize,
    },
    ProbeCheck {
        assignment: Assignment,
        i: usize,
        j: usize,
    },
    Evaluate {
        assignment: Assignment,
    },
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Task {
    id: u64,
    ready: u64,
    depth: usize,
    prev: Option<u64>,
    next: Option<u64>,
    context: Context,
    op: Op,
}
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct FamilyStats {
    work: u64,
    visits: u64,
    fresh: u64,
    mixed: u64,
    modes: [u64; 3],
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct HistoricalAccounting {
    newest_ready_modes: Vec<[u64; 3]>,
    task_bytes_v1: u64,
    peak_state_bytes_v1: u64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct State {
    native: Vec<Scored>,
    ranking: Vec<usize>,
    family_cursor: usize,
    quantum: u8,
    #[serde(deserialize_with = "validation::unique_map")]
    tasks: BTreeMap<u64, Task>,
    #[serde(deserialize_with = "validation::unique_set")]
    fifo: BTreeSet<(usize, u64, u64)>,
    #[serde(deserialize_with = "validation::unique_set")]
    shallow: BTreeSet<(usize, usize, u64)>,
    focus: Vec<Option<u64>>,
    historical_accounting: Option<HistoricalAccounting>,
    next_id: u64,
    next_ready: u64,
    task_bytes: u64,
    score_bytes: u64,
    work: u64,
    evaluations: u64,
    visits: u64,
    native_reuses: u64,
    donors: u64,
    probes_started: u64,
    probe_conflicts: u64,
    probes_completed: u64,
    mixed: u64,
    mixed_identity: u64,
    non_native: u64,
    switches: BTreeMap<usize, u64>,
    visited_sources: BTreeSet<usize>,
    mixed_donor_sources: BTreeSet<usize>,
    families: Vec<FamilyStats>,
    incumbent: Option<Scored>,
    support: Vec<Scored>,
    lost_support: bool,
    peak_state_bytes: u64,
    peak_tasks: usize,
}
impl State {
    fn new(families: usize) -> Self {
        Self {
            native: vec![],
            ranking: vec![],
            family_cursor: 0,
            quantum: 32,
            tasks: BTreeMap::new(),
            fifo: BTreeSet::new(),
            shallow: BTreeSet::new(),
            focus: vec![None; families],
            historical_accounting: None,
            next_id: 0,
            next_ready: 0,
            task_bytes: 0,
            score_bytes: 0,
            work: 0,
            evaluations: 0,
            visits: 0,
            native_reuses: 0,
            donors: 0,
            probes_started: 0,
            probe_conflicts: 0,
            probes_completed: 0,
            mixed: 0,
            mixed_identity: 0,
            non_native: 0,
            switches: BTreeMap::new(),
            visited_sources: BTreeSet::new(),
            mixed_donor_sources: BTreeSet::new(),
            families: vec![FamilyStats::default(); families],
            incumbent: None,
            support: vec![],
            lost_support: false,
            peak_state_bytes: 0,
            peak_tasks: 0,
        }
    }
    // Logical occupancy: compact serialized task payload + 256 bytes for ownership,
    // two ordered indices and live links and allocator/node allowances, with no tombstones.
    fn task_weight(task: &Task) -> io::Result<u64> {
        plus(bytes(task)?, 256)
    }
    fn score_weight(score: &Scored) -> io::Result<u64> {
        plus(bytes(score)?, 64)
    }
    fn occupancy(&self) -> io::Result<u64> {
        plus(
            plus(self.task_bytes, self.score_bytes)?,
            bytes(&self.historical_accounting)?
                + 4096
                + self.families.len() as u64 * 512
                + (self.switches.len()
                    + self.visited_sources.len()
                    + self.mixed_donor_sources.len()) as u64
                    * 64,
        )
    }
    fn insert(&mut self, mut task: Task) -> io::Result<()> {
        task.ready = self.next_ready;
        self.next_ready = plus(self.next_ready, 1)?;
        self.task_bytes = plus(self.task_bytes, Self::task_weight(&task)?)?;
        let f = task.context.family;
        self.fifo.insert((f, task.ready, task.id));
        self.shallow.insert((f, task.depth, task.id));
        ensure(
            self.tasks.insert(task.id, task).is_none(),
            "duplicate task ownership",
        )
    }
    fn spawn(&mut self, context: Context, op: Op) -> io::Result<()> {
        let id = self.next_id;
        self.next_id = plus(id, 1)?;
        self.insert(Task {
            id,
            ready: 0,
            depth: context.depth(),
            prev: None,
            next: None,
            context,
            op,
        })
    }
    fn remove(&mut self, id: u64) -> io::Result<Task> {
        let task = self
            .tasks
            .remove(&id)
            .ok_or_else(|| invalid("missing indexed task"))?;
        let f = task.context.family;
        ensure(
            self.fifo.remove(&(f, task.ready, id)) && self.shallow.remove(&(f, task.depth, id)),
            "inconsistent task indices",
        )?;
        self.task_bytes = self
            .task_bytes
            .checked_sub(Self::task_weight(&task)?)
            .ok_or_else(|| invalid("occupancy underflow"))?;
        Ok(task)
    }
    fn select(&self) -> Option<(usize, u8, u64, usize)> {
        for skip in 0..self.ranking.len() {
            let cursor = (self.family_cursor + skip) % self.ranking.len();
            let f = self.ranking[cursor];
            let fifo = self.fifo.range((f, 0, 0)..=(f, u64::MAX, u64::MAX)).next();
            if let Some(&(_, _, first)) = fifo {
                let mode = (self.families[f].work % 3) as usize;
                let id = match mode {
                    0 => first,
                    1 => {
                        self.shallow
                            .range((f, 0, 0)..=(f, usize::MAX, u64::MAX))
                            .next()
                            .unwrap()
                            .2
                    }
                    _ => self.focus[f].unwrap_or(first),
                };
                return Some((cursor, if skip == 0 { self.quantum } else { 32 }, id, mode));
            }
        }
        None
    }
    fn record(
        &mut self,
        scored: Scored,
        reuse: bool,
        initialization: bool,
        epsilon: f64,
        limits: &Limits,
        e: &Evaluator<'_>,
        ledger: &mut std::fs::File,
        task: Option<u64>,
    ) -> io::Result<()> {
        self.visits = plus(self.visits, 1)?;
        let family = scored.assignment.family;
        self.families[family].visits += 1;
        if reuse {
            self.native_reuses = plus(self.native_reuses, 1)?;
        } else {
            self.evaluations = plus(self.evaluations, 1)?;
            self.families[family].fresh += 1;
        }
        let switches = scored
            .assignment
            .routes
            .iter()
            .map(|r| r.segments.len() - 1)
            .sum::<usize>();
        let mixed = scored
            .assignment
            .routes
            .iter()
            .any(|r| r.segments.iter().any(|s| s.source != r.segments[0].source));
        let identity = scored.assignment.routes.iter().any(|r| {
            r.segments.iter().any(|s| {
                e.graph.lanes[s.source].family != e.graph.lanes[r.segments[0].source].family
            })
        });
        if !reuse {
            *self.switches.entry(switches).or_default() += 1;
            if switches > 0 {
                self.non_native += 1;
            }
            if mixed {
                self.mixed += 1;
                self.families[family].mixed += 1;
            }
            if identity {
                self.mixed_identity += 1;
            }
            for r in &scored.assignment.routes {
                for s in &r.segments {
                    self.visited_sources.insert(s.source);
                    if s.source != r.segments[0].source {
                        self.mixed_donor_sources.insert(s.source);
                    }
                }
            }
        }
        if self
            .incumbent
            .as_ref()
            .is_none_or(|i| scored.objective() < i.objective())
        {
            self.incumbent = Some(scored.clone());
            self.support
                .retain(|s| s.objective() - scored.objective() <= epsilon);
        }
        if scored.objective() - self.incumbent.as_ref().unwrap().objective() <= epsilon
            && !self
                .support
                .iter()
                .any(|s| s.assignment == scored.assignment)
        {
            if self.support.len() < limits.max_optima {
                self.support.push(scored.clone());
            } else {
                self.lost_support = true;
            }
        }
        if initialization {
            self.native.push(scored.clone());
        }
        self.score_bytes = 0;
        for s in self
            .native
            .iter()
            .chain(self.incumbent.iter())
            .chain(&self.support)
        {
            self.score_bytes = plus(self.score_bytes, Self::score_weight(s)?)?;
        }
        let row = serde_json::json!({"visit":self.visits,"evaluation":self.evaluations,"work":self.work,
            "task":task,"kind":if initialization {"native-initialization"} else if reuse {"native-score-reuse"} else {"fresh-evaluation"},
            "switches":switches,"cross_source_route":mixed,"cross_identity_route":identity,"scored":scored.json(),
            "assignment_checksum":format!("{:016x}",genome::checksum(&serde_json::to_vec(&scored.assignment).map_err(io::Error::other)?))});
        serde_json::to_writer(&mut *ledger, &row).map_err(io::Error::other)?;
        ledger.write_all(b"\n")
    }
}

pub fn run(
    e: &mut Evaluator<'_>,
    sample_checksum: String,
    cache_terms: usize,
    options: Options,
    out: &Path,
) -> io::Result<serde_json::Value> {
    ensure(
        options.max_work > 0
            && options.max_evaluations > 0
            && options.max_state_bytes > 0
            && options.max_optima > 0
            && options.tie_epsilon.is_finite()
            && options.tie_epsilon >= 0.0,
        "invalid guided search budgets",
    )?;
    let limits = Limits {
        max_work: options.max_work,
        max_evaluations: options.max_evaluations,
        max_state_bytes: options.max_state_bytes,
        max_optima: options.max_optima,
    };
    let bindings = checkpoint::Bindings {
        policy: identity(),
        backend: routes::compiler_identity(),
        graph: e.graph.digest()?,
        sample: sample_checksum,
        depth_bits: e.depth.to_bits(),
        background_bits: e.background.to_bits(),
        tie_bits: options.tie_epsilon.to_bits(),
        max_feature_terms: e.max_terms,
        cache_terms,
        count_policy: genome::COUNT_POLICY.into(),
    };
    ensure(!options.convert_exact_v1_to_v2,
        "direct v1 conversion is not supported by v3; use the preserved v2 binary first, then --update-source-pair-order")?;
    if options.update_source_pair_order {
        return source_order::update(
            e,
            options.resume_from.as_ref().unwrap(),
            out,
            bindings,
            limits,
            options.extend_budgets,
        );
    }
    let (mut state, budgets, parent) = if let Some(path) = &options.resume_from {
        let (s, b, p) = checkpoint::load(path, &bindings, &limits, options.extend_budgets)?;
        (s, b, Some(p))
    } else {
        (
            State::new(e.graph.families.len()),
            vec![limits.clone()],
            None,
        )
    };
    validation::graph(&state, e.graph)?;
    let metadata_bytes = bytes(&(&bindings, &budgets, &parent))?;
    ensure(
        plus(state.occupancy()?, metadata_bytes)? <= limits.max_state_bytes,
        "state budget cannot hold restored state and checkpoint metadata",
    )?;
    // A budget stop may occur before any operation updates the high-water mark.
    state.peak_state_bytes = state.peak_state_bytes.max(state.occupancy()?);
    let mut ledger = checkpoint::ledger(out, parent.as_ref())?;
    let mut stop = "domain-exhausted";
    loop {
        let initializing = state.native.len() < e.graph.families.len();
        let selected = if initializing { None } else { state.select() };
        if !initializing && selected.is_none() {
            break;
        }
        if state.work >= limits.max_work {
            stop = "work-budget-exhausted";
            break;
        }
        let (family, task_size) = if initializing {
            (state.native.len(), 0)
        } else {
            let t = &state.tasks[&selected.unwrap().2];
            (t.context.family, State::task_weight(t)?)
        };
        let native = if initializing {
            e.graph.native_assignment(family)?
        } else {
            state.native[family].assignment.clone()
        };
        // Preflight reserves the maximum retained growth of ONE operation: at most
        // three task/score copies, one native tail, indices and bounded counters.
        // Eight copies plus slack deliberately overcharge; never a claimed RSS cap.
        let reserve = reservation(task_size, &native, e.graph.k)?;
        if plus(plus(state.occupancy()?, metadata_bytes)?, reserve)? > limits.max_state_bytes {
            stop = "state-budget-exhausted";
            break;
        }
        let needs_evaluation = initializing
            || selected.is_some_and(|(_, _, id, _)| match &state.tasks[&id].op {
                Op::Evaluate { assignment } => {
                    !state.native.iter().any(|s| s.assignment == *assignment)
                }
                _ => false,
            });
        if needs_evaluation && state.evaluations >= limits.max_evaluations {
            stop = "evaluation-budget-exhausted";
            break;
        }
        // Reserve clock/counter widths before ownership, selectors or work mutate.
        validation::counter_room(&state)?;
        state.work = plus(state.work, 1)?;
        if initializing {
            let result = e.evaluate(&native)?;
            state.record(
                Scored {
                    assignment: result.assignment,
                    objective_bits: result.relative_objective.to_bits(),
                },
                false,
                true,
                options.tie_epsilon,
                &limits,
                e,
                &mut ledger,
                None,
            )?;
            state.spawn(machine::initial(e.graph, family, 0, vec![]), Op::Start)?;
            if state.native.len() == e.graph.families.len() {
                state.ranking = (0..e.graph.families.len()).collect();
                state.ranking.sort_by(|&a, &b| {
                    state.native[a]
                        .objective()
                        .total_cmp(&state.native[b].objective())
                        .then_with(|| {
                            e.graph.families[a]
                                .identity
                                .cmp(&e.graph.families[b].identity)
                        })
                });
            }
        } else {
            let (cursor, quantum, id, mode) = selected.unwrap();
            state.family_cursor = cursor;
            state.quantum = quantum;
            state.families[family].work += 1;
            state.families[family].modes[mode] += 1;
            if mode == 2 {
                state.seed_focus(id);
            }
            let task = state.remove(id)?;
            machine::advance(
                &mut state,
                task,
                e,
                &limits,
                options.tie_epsilon,
                &mut ledger,
            )?;
            state.quantum -= 1;
            if state.quantum == 0 {
                state.quantum = 32;
                state.family_cursor = (cursor + 1) % state.ranking.len();
            }
        }
        ensure(
            plus(state.occupancy()?, metadata_bytes)? <= limits.max_state_bytes,
            "internal reservation bound exceeded",
        )?;
        state.peak_state_bytes = state.peak_state_bytes.max(state.occupancy()?);
        state.peak_tasks = state.peak_tasks.max(state.tasks.len());
    }
    let initialized = state.native.len() == e.graph.families.len();
    let exhausted = initialized && state.tasks.is_empty();
    if let Some(s) = &state.incumbent {
        genome::write_json(&out.join("incumbent-assignment.json"), &s.assignment)?;
    }
    let mut result = serde_json::json!({"method":POLICY_VERSION,"policy_identity":bindings.policy,"status":stop,
        "budget":limits,"work":state.work,"evaluations":state.evaluations,"complete_assignment_visits":state.visits,
        "native_score_reuses":state.native_reuses,"native_families_evaluated":state.native.len(),
        "native_initializations_complete":initialized,"native_ranked_family_order":state.ranking,
        "native_ranked_identities":state.ranking.iter().map(|&f| &e.graph.families[f].identity).collect::<Vec<_>>(),
        "native_scores":state.native.iter().map(Scored::json).collect::<Vec<_>>(),"family_statistics":state.families,
        "donor_transitions_examined":state.donors,"mixed_assignments_evaluated":state.mixed,
        "completion_probes_started":state.probes_started,"completion_probe_conflicts":state.probe_conflicts,
        "completion_probes_completed":state.probes_completed,
        "mixed_identity_assignments_evaluated":state.mixed_identity,"non_native_assignments_evaluated":state.non_native,
        "maximum_switches_evaluated":state.switches.keys().next_back().copied().unwrap_or(0),
        "switch_distribution":state.switches,"visited_source_ids":state.visited_sources,
        "frontier_tasks":state.tasks.len(),"peak_frontier_tasks":state.peak_tasks,
        "retained_state_bytes":plus(state.occupancy()?,metadata_bytes)?,"checkpoint_metadata_bytes":metadata_bytes,
        "peak_state_bytes":plus(state.peak_state_bytes,metadata_bytes)?,
        "occupancy_kind":"logical-serialized-payload-plus-live-links-v2;not-RSS;atomic-evaluator-excluded",
        "search_exhausted":exhausted,"global_optimum_certified":exhausted,
        "correlated_optima_complete":exhausted && !state.lost_support,"sticky_support_loss":state.lost_support,
        "incumbent":state.incumbent.as_ref().map(Scored::json),"correlated_optima":state.support.iter().map(Scored::json).collect::<Vec<_>>(),
        "lower_bound":e.lower_bound(),"generation_rule_complete":true,"native_evaluation_complete":true,
        "biological_topology_complete":false,"sequence_emission_authorized":false});
    result["historical_accounting"] =
        serde_json::to_value(&state.historical_accounting).map_err(io::Error::other)?;
    result["focused_services_by_family"] = serde_json::to_value(
        state
            .families
            .iter()
            .enumerate()
            .map(|(i, f)| {
                f.modes[2]
                    - state
                        .historical_accounting
                        .as_ref()
                        .map_or(0, |h| h.newest_ready_modes[i][2])
            })
            .collect::<Vec<_>>(),
    )
    .map_err(io::Error::other)?;
    result["mixed_donor_source_ids"] =
        serde_json::to_value(&state.mixed_donor_sources).map_err(io::Error::other)?;
    checkpoint::save(out, bindings, parent, budgets, state, &mut ledger, None)?;
    Ok(result)
}

fn reservation(task_size: u64, native: &Assignment, k: u64) -> io::Result<u64> {
    plus(
        plus(task_size, bytes(native)?)?,
        k.checked_mul(16)
            .ok_or_else(|| invalid("reservation overflow"))?,
    )?
    .checked_mul(8)
    .and_then(|n| n.checked_add(8192))
    .ok_or_else(|| invalid("reservation overflow"))
}
