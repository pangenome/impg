//! Native-seeded, budgeted successive-switch DFS. A path state traverses positive
//! source spans and branches through indexed DNA hubs, so mixed donor chains do
//! not require a preexisting single-donor bridge or chromosome alternative list.
use super::*;
use storage::Port;
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SearchBudget {
    pub max_work: u64,
    pub max_evaluations: u64,
    pub max_frontier: usize,
    pub max_optima: usize,
    pub tie_epsilon: f64,
}
#[derive(Clone, Debug, Serialize)]
pub struct Scored {
    pub assignment: Assignment,
    pub relative_objective: f64,
}
#[derive(Debug, Serialize)]
pub struct SearchResult {
    pub version: u32,
    pub model: String,
    pub graph_checksum: String,
    pub method: String,
    pub objective_kind: String,
    pub universe: String,
    pub budget: SearchBudget,
    pub status: String,
    pub work: u64,
    pub evaluations: u64,
    pub native_families_evaluated: usize,
    pub native_initializations_complete: bool,
    /// Hub-successor work, including same-source jumps; not evidence of mixing.
    pub donor_transitions_examined: u64,
    /// At least one complete route traverses more than one source path ID.
    pub mixed_assignments_evaluated: u64,
    /// At least one complete route traverses more than one sample#haplotype identity.
    pub mixed_identity_assignments_evaluated: u64,
    pub non_native_assignments_evaluated: u64,
    pub maximum_switches_evaluated: usize,
    pub family_searches_exhausted: usize,
    pub search_exhausted: bool,
    pub global_optimum_certified: bool,
    pub correlated_optima_complete: bool,
    pub incumbent: Option<Scored>,
    pub correlated_optima: Vec<Scored>,
    pub lower_bound: Option<f64>,
    pub lower_bound_kind: String,
    pub generation_rule_complete: bool,
    pub native_evaluation_complete: bool,
    pub biological_topology_complete: bool,
    pub sequence_emission_authorized: bool,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
struct State {
    family: usize,
    slot: usize,
    completed: Vec<Route>,
    segments: Vec<Segment>,
    source: usize,
    cut: u64,
    reverse: bool,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
struct Frame {
    state: State,
    terminal_pending: bool,
    next_port: i64,
    step: i64,
    pending: Option<(Port, u64, u64)>,
}
enum Action {
    Continue,
    Child(State),
    Evaluate(Assignment),
    Done,
}
impl Frame {
    fn new(state: State, e: &mut Evaluator<'_>) -> io::Result<Self> {
        let g = e.graph;
        let shift = if state.reverse {
            g.k - g.cut_offset
        } else {
            g.cut_offset
        };
        let threshold = if state.reverse {
            state.cut.saturating_sub(shift)
        } else {
            state
                .cut
                .saturating_sub(shift)
                .saturating_add(u64::from(state.cut >= shift))
        };
        let count = g.lanes[state.source].port_count;
        let mut file = e.ports.source_file(g, state.source)?;
        let (mut lo, mut hi) = (0, count);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if storage::port_at(&mut file, g.k as usize, mid)?.anchor < threshold {
                lo = mid + 1
            } else {
                hi = mid
            }
        }
        let next_port = if state.reverse {
            lo as i64 - 1
        } else {
            lo as i64
        };
        let step = if state.reverse { -1 } else { 1 };
        Ok(Self {
            state,
            terminal_pending: true,
            next_port,
            step,
            pending: None,
        })
    }
    fn feasible(&self, piece: &Segment) -> bool {
        evaluate::resources_feasible(
            self.state
                .completed
                .iter()
                .flat_map(|r| r.segments.clone())
                .chain(self.state.segments.clone())
                .chain(std::iter::once(piece.clone())),
        )
    }
    fn traversal(&self, to: u64) -> Segment {
        Segment {
            source: self.state.source,
            start: self.state.cut.min(to),
            end: self.state.cut.max(to),
            reverse: self.state.reverse,
        }
    }
    /// One unit scans at most one native port or hub member. Skipped repeat and
    /// infeasible alternatives still consume work; none are removed from the rule.
    fn advance(&mut self, e: &mut Evaluator<'_>, donors: &mut u64) -> io::Result<Action> {
        let g = e.graph;
        let target = g.families[self.state.family].paths[self.state.slot];
        if self.terminal_pending {
            self.terminal_pending = false;
            if self.state.source == target
                && !self.state.reverse
                && self.state.cut < g.lanes[target].length
            {
                let piece = self.traversal(g.lanes[target].length);
                if self.feasible(&piece) {
                    let mut segments = self.state.segments.clone();
                    segments.push(piece);
                    let mut completed = self.state.completed.clone();
                    completed.push(Route { segments }.normalize());
                    if self.state.slot + 1 == g.families[self.state.family].paths.len() {
                        return Ok(Action::Evaluate(Assignment {
                            version: VERSION,
                            model: MODEL.into(),
                            graph_checksum: e.graph_checksum.clone(),
                            family: self.state.family,
                            routes: completed,
                        }));
                    }
                    let slot = self.state.slot + 1;
                    return Ok(Action::Child(State {
                        family: self.state.family,
                        slot,
                        completed,
                        segments: Vec::new(),
                        source: g.families[self.state.family].paths[slot],
                        cut: 0,
                        reverse: false,
                    }));
                }
            }
            return Ok(Action::Continue);
        }
        if let Some((port, next, end)) = &mut self.pending {
            if *next < *end {
                let donor = storage::port_at(&mut e.ports.global, g.k as usize, *next)?;
                *next += 1;
                add(donors, 1)?;
                let outgoing = port.clone();
                require(
                    donor.word == outgoing.word
                        && donor.source < g.lanes.len()
                        && donor.anchor + g.k <= g.lanes[donor.source].length,
                    "corrupt hub traversal",
                )?;
                let cut = outgoing.cut(g.k);
                if donor.source == outgoing.source
                    && donor.cut(g.k) == cut
                    && donor.reverse == outgoing.reverse
                {
                    return Ok(Action::Continue);
                }
                let piece = self.traversal(cut);
                if !self.feasible(&piece) {
                    return Ok(Action::Continue);
                }
                let mut segments = self.state.segments.clone();
                segments.push(piece);
                return Ok(Action::Child(State {
                    family: self.state.family,
                    slot: self.state.slot,
                    completed: self.state.completed.clone(),
                    segments,
                    source: donor.source,
                    cut: donor.cut(g.k),
                    reverse: donor.reverse,
                }));
            }
            self.pending = None;
            return Ok(Action::Continue);
        }
        let count = g.lanes[self.state.source].port_count as i64;
        if self.next_port < 0 || self.next_port >= count {
            return Ok(Action::Done);
        }
        let mut file = e.ports.source_file(g, self.state.source)?;
        let p = storage::port_at(&mut file, g.k as usize, self.next_port as u64)?;
        self.next_port += self.step;
        require(
            p.source == self.state.source && p.anchor + g.k <= g.lanes[p.source].length,
            "corrupt source port index",
        )?;
        let cut = p.cut(g.k);
        if p.reverse != self.state.reverse
            || (!p.reverse && cut <= self.state.cut)
            || (p.reverse && cut >= self.state.cut)
        {
            return Ok(Action::Continue);
        }
        let (lo, hi) = storage::bucket(&mut e.ports.global, g.k as usize, g.port_count, &p.word)?;
        require(lo < hi, "source switch missing hub membership")?;
        self.pending = Some((p, lo, hi));
        Ok(Action::Continue)
    }
}
pub fn search(e: &mut Evaluator<'_>, budget: SearchBudget, out: &Path) -> io::Result<SearchResult> {
    require(
        budget.max_work > 0
            && budget.max_evaluations > 0
            && budget.max_frontier > 0
            && budget.max_optima > 0
            && budget.tie_epsilon.is_finite()
            && budget.tie_epsilon >= 0.0,
        "invalid route exploration budget",
    )?;
    let mut ledger = BufWriter::new(File::create(out.join("evaluations.jsonl"))?);
    let mut incumbent: Option<Scored> = None;
    let mut optima: Vec<Scored> = Vec::new();
    let mut discarded = false;
    let (mut work, mut evaluations, mut native, mut donors, mut mixed, mut max_switches) =
        (0u64, 0u64, 0usize, 0u64, 0u64, 0usize);
    let (mut non_native, mut mixed_identity) = (0u64, 0u64);
    let mut record = |a: Assignment, e: &mut Evaluator<'_>| -> io::Result<()> {
        let result = e.evaluate(&a)?;
        let switches: usize = result
            .assignment
            .routes
            .iter()
            .map(|r| r.segments.len() - 1)
            .sum();
        if switches > 0 {
            non_native += 1;
        }
        let cross_source = result
            .assignment
            .routes
            .iter()
            .any(|r| r.segments.iter().any(|s| s.source != r.segments[0].source));
        let cross_identity = result.assignment.routes.iter().any(|r| {
            let first = e.graph.lanes[r.segments[0].source].family;
            r.segments
                .iter()
                .any(|s| e.graph.lanes[s.source].family != first)
        });
        if cross_source {
            mixed += 1;
        }
        if cross_identity {
            mixed_identity += 1;
        }
        max_switches = max_switches.max(switches);
        evaluations += 1;
        let scored = Scored {
            assignment: result.assignment,
            relative_objective: result.relative_objective,
        };
        storage::line(
            &mut ledger,
            &serde_json::json!({"evaluation":evaluations,"native_evaluation_complete":true,"switches":switches,"cross_source_route":cross_source,"cross_identity_route":cross_identity,"scored":scored}),
        )?;
        if incumbent
            .as_ref()
            .is_none_or(|i| scored.relative_objective < i.relative_objective)
        {
            incumbent = Some(scored.clone());
            optima
                .retain(|o| o.relative_objective - scored.relative_objective <= budget.tie_epsilon);
        }
        if scored.relative_objective - incumbent.as_ref().unwrap().relative_objective
            <= budget.tie_epsilon
            && !optima.iter().any(|o| o.assignment == scored.assignment)
        {
            if optima.len() < budget.max_optima {
                optima.push(scored);
            } else {
                discarded = true;
            }
        }
        Ok(())
    };
    // Every family is an initialization, including fragmented and zero-anchor
    // inventories. A too-small budget reports native initialization incompleteness.
    for f in 0..e.graph.families.len() {
        if work >= budget.max_work || native as u64 >= budget.max_evaluations {
            break;
        }
        work += 1;
        record(e.graph.native_assignment(f)?, e)?;
        native += 1;
    }
    let mut family_done = 0;
    let mut stack: Vec<Frame> = Vec::new();
    let mut pending_child = None;
    let mut stop = "work-budget-exhausted";
    // record owns evaluation count while mutably borrowed; track leaves separately.
    let mut leaves = native as u64;
    if native == e.graph.families.len() {
        'families: for family in 0..e.graph.families.len() {
            let source = e.graph.families[family].paths[0];
            stack.push(Frame::new(
                State {
                    family,
                    slot: 0,
                    completed: Vec::new(),
                    segments: Vec::new(),
                    source,
                    cut: 0,
                    reverse: false,
                },
                e,
            )?);
            while let Some(frame) = stack.last_mut() {
                if work >= budget.max_work {
                    stop = "work-budget-exhausted";
                    break 'families;
                }
                if leaves >= budget.max_evaluations {
                    stop = "evaluation-budget-exhausted";
                    break 'families;
                }
                work += 1;
                match frame.advance(e, &mut donors)? {
                    Action::Continue => (),
                    Action::Done => {
                        stack.pop();
                    }
                    Action::Child(state) => {
                        if stack.len() >= budget.max_frontier {
                            pending_child = Some(state);
                            stop = "frontier-budget-exhausted";
                            break 'families;
                        }
                        stack.push(Frame::new(state, e)?);
                    }
                    Action::Evaluate(a) => {
                        record(a, e)?;
                        leaves += 1;
                    }
                }
            }
            family_done += 1;
        }
    } else if native as u64 >= budget.max_evaluations {
        stop = "evaluation-budget-exhausted";
    }
    drop(record);
    ledger.flush()?;
    let exhausted = family_done == e.graph.families.len();
    let complete = exhausted && !discarded;
    // Exact traversal cursors and the not-yet-pushed child are preserved. Resume
    // is not exposed yet; this is an auditable checkpoint, never claimed success.
    write_json(
        &out.join("frontier.json"),
        &serde_json::json!({"version":VERSION,"model":MODEL,"graph_checksum":e.graph.digest()?,"family_searches_exhausted":family_done,"native_families_evaluated":native,"frames":stack,"pending_child":pending_child,"resume_supported":false,"budget":budget}),
    )?;
    Ok(SearchResult{version:VERSION,model:MODEL.into(),graph_checksum:e.graph.digest()?,method:"native-seeded-indexed-successive-donor-switch-DFS-v1".into(),objective_kind:"fixed-universe-background-relative-NLL-only".into(),universe:UNIVERSE.into(),status:if exhausted{if complete{"domain-exhausted-optimum-and-correlated-support-certified"}else{"domain-exhausted-optimum-certified-support-incomplete"}}else{stop}.into(),budget,work,evaluations,native_families_evaluated:native,native_initializations_complete:native==e.graph.families.len(),donor_transitions_examined:donors,mixed_assignments_evaluated:mixed,mixed_identity_assignments_evaluated:mixed_identity,non_native_assignments_evaluated:non_native,maximum_switches_evaluated:max_switches,family_searches_exhausted:family_done,search_exhausted:exhausted,global_optimum_certified:exhausted,correlated_optima_complete:complete,incumbent,correlated_optima:optima,lower_bound:e.lower_bound(),lower_bound_kind:"independent nonnegative-signal relaxation; ln(r)<=r-1; outward-rounded -(c-beta)^2/beta; null means negative-infinity overflow".into(),generation_rule_complete:true,native_evaluation_complete:true,biological_topology_complete:false,sequence_emission_authorized:false})
}
