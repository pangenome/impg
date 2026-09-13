//! Alternating discovery/active service and producer/closure service, with backpressure.
use super::*;
#[derive(Serialize)]
enum SourceStage {
    Recover,
    Bound {
        lo: u64,
        hi: u64,
    },
    Scan {
        base: u64,
        permutation: SourcePermutation,
    },
}
#[derive(Serialize)]
enum HubStage {
    Bound {
        lower: Option<u64>,
        lo: u64,
        hi: u64,
    },
    Members {
        next: u64,
        end: u64,
    },
}
#[derive(Serialize)]
struct Hub {
    port: Port,
    stage: HubStage,
}
#[derive(Serialize)]
struct Producer {
    task: u64,
    original: Context,
    next_entry: usize,
    current: Option<Context>,
    source: SourceStage,
    hub: Option<Hub>,
    direct: Option<Port>,
    weight: u64,
}
#[derive(Serialize)]
enum ChainStage {
    Tail { slot: usize },
    Capacity { i: usize, j: usize },
    Evaluate,
}
#[derive(Serialize)]
struct Chain {
    task: u64,
    entry: usize,
    assignment: Assignment,
    stage: ChainStage,
    weight: u64,
}
struct Pending {
    task: Task,
    index: usize,
}
pub struct Engine<'a> {
    options: &'a Options,
    pub snapshot: snapshot::Snapshot,
    mem: Memory,
    counts: Counts,
    producers: VecDeque<Producer>,
    chains: VecDeque<Chain>,
    pending: Option<Pending>,
    cache: BTreeMap<Assignment, u64>,
    scores: BTreeSet<u64>,
    input_first: bool,
    input_end: bool,
    producer_turn: bool,
    discovery_turn: bool,
    stop: String,
    ledger: File,
}
impl<'a> Engine<'a> {
    pub fn new(
        options: &'a Options,
        snapshot: snapshot::Snapshot,
        mut mem: Memory,
    ) -> io::Result<Self> {
        let slots = options
            .producer_slots
            .checked_add(options.chain_slots)
            .ok_or_else(|| invalid("queue size overflow"))?;
        mem.reserve(
            (slots as u64)
                .checked_mul(4096)
                .ok_or_else(|| invalid("queue reservation overflow"))?,
        )?;
        let ledger = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(options.out_dir.join("proposals.jsonl"))?;
        let mut scores = BTreeSet::new();
        for n in &snapshot.native {
            mem.reserve(128)?;
            scores.insert(n.objective_bits);
        }
        Ok(Self {
            options,
            snapshot,
            mem,
            counts: Counts::default(),
            producers: VecDeque::with_capacity(options.producer_slots),
            chains: VecDeque::with_capacity(options.chain_slots),
            pending: None,
            cache: BTreeMap::new(),
            scores,
            input_first: true,
            input_end: false,
            producer_turn: true,
            discovery_turn: true,
            stop: String::new(),
            ledger,
        })
    }
    fn discovery(&mut self, input: &mut stream::Json, e: &Evaluator<'_>) -> io::Result<()> {
        if self.pending.is_none() {
            if !input.next(&mut self.input_first, b'}')? {
                self.input_end = true;
                return Ok(());
            }
            let key = input.key()?;
            let task: Task = input.value()?;
            ensure(key == task.id.to_string(), "discovery task ID mismatch")?;
            self.counts.inspected += 1;
            // Scored assignments and native tails are not proposal seeds.
            if matches!(
                task.op,
                Op::Evaluate { .. } | Op::Probe { .. } | Op::ProbeCheck { .. }
            ) {
                self.counts.skipped_scored += 1;
                return Ok(());
            }
            self.pending = Some(Pending { task, index: 0 });
            return Ok(());
        }
        let p = self.pending.as_mut().unwrap();
        let c = &p.task.context;
        // ONE span or route-boundary geometry check per charged service. The
        // existing index owns this cursor; no additional retained tracking grows.
        let changed = if let Some(changed) = eligibility(c, p.index, |slot| {
            let source = e.graph.families[c.family].paths[slot];
            (source, e.graph.lanes[source].length)
        }) {
            p.index += 1;
            changed
        } else {
            self.counts.skipped_native += 1;
            self.pending = None;
            return Ok(());
        };
        if !changed {
            return Ok(());
        }
        let w = weight(plus(
            bytes(c)?,
            e.graph
                .k
                .checked_mul(2)
                .ok_or_else(|| invalid("port reservation overflow"))?,
        )?)?;
        self.mem.reserve(w)?; // original + current prefix + cursor/port slack, before retaining
        let p = self.pending.take().unwrap();
        let direct = match p.task.op {
            Op::HubBound { port, .. } | Op::HubScan { port, .. } => Some(port),
            _ => None,
        };
        if direct.is_some() {
            self.counts.deferred += 1;
        }
        self.producers.push_back(Producer {
            task: p.task.id,
            original: p.task.context,
            next_entry: 0,
            current: None,
            source: SourceStage::Recover,
            hub: None,
            direct,
            weight: w,
        });
        self.counts.admitted += 1;
        Ok(())
    }
    pub fn drive(&mut self, input: &mut stream::Json, e: &mut Evaluator<'_>) -> io::Result<()> {
        loop {
            let discovery = !self.input_end && self.producers.len() < self.options.producer_slots;
            let active = !self.producers.is_empty() || !self.chains.is_empty();
            if !discovery && !active {
                self.stop = "proposal-subset-exhausted".into();
                break;
            }
            if self.counts.work >= self.options.max_work {
                self.stop = "work-budget-exhausted".into();
                break;
            }
            self.counts.work += 1;
            let result = if discovery && (self.discovery_turn || !active) {
                self.discovery_turn = false;
                self.discovery(input, e)
            } else {
                self.discovery_turn = true;
                self.active(e)
            };
            if let Err(err) = result {
                if [
                    "state-budget-exhausted",
                    "evaluation-budget-exhausted",
                    "distinct-budget-exhausted",
                    "input-record-byte-budget-exhausted",
                ]
                .contains(&err.to_string().as_str())
                {
                    self.stop = err.to_string();
                    break;
                }
                return Err(err);
            }
        }
        self.ledger.flush()?;
        self.ledger.sync_all()
    }
    fn active(&mut self, e: &mut Evaluator<'_>) -> io::Result<()> {
        if !self.producers.is_empty() && (self.producer_turn || self.chains.is_empty()) {
            self.producer_turn = false;
            let mut p = self.producers.pop_front().unwrap();
            let result = self.produce(&mut p, e);
            if matches!(result, Ok(true)) {
                self.mem.release(p.weight);
            } else {
                self.producers.push_back(p);
            }
            result.map(|_| ())
        } else {
            self.producer_turn = true;
            let mut c = self.chains.pop_front().unwrap();
            let result = self.close(&mut c, e);
            if matches!(result, Ok(true)) {
                self.mem.release(c.weight);
            } else {
                self.chains.push_back(c);
            }
            result.map(|_| ())
        }
    }
    fn produce(&mut self, p: &mut Producer, e: &mut Evaluator<'_>) -> io::Result<bool> {
        let g = e.graph;
        if p.hub.is_some() {
            return self.member(p, e).map(|_| false);
        }
        match &mut p.source {
            SourceStage::Recover => {
                // Deferred positive exits receive service without first exhausting reconstructed entries.
                if let Some(port) = p.direct.take() {
                    p.current = Some(p.original.clone());
                    if positive(p.current.as_ref().unwrap(), &port, g.k)? {
                        p.hub = Some(Hub {
                            port,
                            stage: HubStage::Bound {
                                lower: None,
                                lo: 0,
                                hi: g.port_count,
                            },
                        });
                    }
                    return Ok(false);
                }
                let i = p.next_entry;
                if i > p.original.segments.len() {
                    return Ok(true);
                }
                p.next_entry += 1;
                self.counts.prefix_positions += 1;
                // Position zero has no donor entry; it cannot manufacture a changed prefix.
                if i == 0 && p.original.completed.is_empty() {
                    return Ok(false);
                }
                let c = entry_context(&p.original, i)?;
                if i < p.original.segments.len() {
                    self.counts.recovered_earlier += 1;
                }
                p.source = SourceStage::Bound {
                    lo: 0,
                    hi: g.lanes[c.source].port_count,
                };
                p.current = Some(c);
            }
            SourceStage::Bound { lo, hi } => {
                self.counts.source_bounds += 1;
                let c = p.current.as_ref().unwrap();
                if *lo < *hi {
                    let mid = *lo + (*hi - *lo) / 2;
                    let mut f = e.ports.source_file(g, c.source)?;
                    let port = adapter::read(&mut f, g.k, g.lanes[c.source].port_count, mid)?;
                    adapter::verified(&port, g)?;
                    ensure(port.source == c.source, "source index mismatch")?;
                    let shift = if c.reverse {
                        g.k - g.cut_offset
                    } else {
                        g.cut_offset
                    };
                    let threshold = if c.reverse {
                        c.cut.saturating_sub(shift)
                    } else {
                        plus(c.cut.saturating_sub(shift), u64::from(c.cut >= shift))?
                    };
                    if port.anchor < threshold {
                        *lo = mid + 1
                    } else {
                        *hi = mid
                    }
                } else {
                    let (base, len) = if c.reverse {
                        (0, *lo)
                    } else {
                        (*lo, g.lanes[c.source].port_count - *lo)
                    };
                    p.source = SourceStage::Scan {
                        base,
                        permutation: SourcePermutation::new(len),
                    };
                }
            }
            SourceStage::Scan { base, permutation } => {
                self.counts.source_ordinals += 1;
                let c = p.current.as_ref().unwrap();
                match permutation.next() {
                    None => {
                        p.source = SourceStage::Recover;
                        p.current = None;
                    }
                    Some(None) => (),
                    Some(Some(index)) => {
                        let mut f = e.ports.source_file(g, c.source)?;
                        let port = adapter::read(
                            &mut f,
                            g.k,
                            g.lanes[c.source].port_count,
                            plus(*base, index)?,
                        )?;
                        adapter::verified(&port, g)?;
                        ensure(port.source == c.source, "source index mismatch")?;
                        if positive(c, &port, g.k)? {
                            p.hub = Some(Hub {
                                port,
                                stage: HubStage::Bound {
                                    lower: None,
                                    lo: 0,
                                    hi: g.port_count,
                                },
                            });
                        }
                    }
                }
            }
        }
        Ok(false)
    }
    fn member(&mut self, p: &mut Producer, e: &mut Evaluator<'_>) -> io::Result<()> {
        let g = e.graph;
        let c = p.current.as_ref().unwrap();
        let target = g.families[c.family].paths[c.slot];
        let hub = p.hub.as_mut().unwrap();
        match &mut hub.stage {
            HubStage::Bound { lower, lo, hi } => {
                self.counts.target_bounds += 1;
                if *lo < *hi {
                    let mid = *lo + (*hi - *lo) / 2;
                    let q = adapter::read(&mut e.ports.global, g.k, g.port_count, mid)?;
                    adapter::verified(&q, g)?;
                    // Port::Ord is word THEN DECODED NUMERIC source, never raw LE bytes.
                    if target_before(&q, &hub.port.word, target, lower.is_some()) {
                        *lo = mid + 1
                    } else {
                        *hi = mid
                    }
                } else if let Some(base) = *lower {
                    hub.stage = HubStage::Members {
                        next: base,
                        end: *lo,
                    };
                } else {
                    *lower = Some(*lo);
                    *hi = g.port_count;
                }
            }
            HubStage::Members { next, end } => {
                if *next == *end {
                    p.hub = None;
                    return Ok(());
                }
                if self.chains.len() == self.options.chain_slots {
                    return Ok(());
                } // charged backpressure, cursor retained
                let q = adapter::read(&mut e.ports.global, g.k, g.port_count, *next)?;
                adapter::verified(&q, g)?;
                ensure(
                    q.word == hub.port.word && q.source == target,
                    "numeric target interval mismatch",
                )?;
                let cut = q.cut(g.k)?;
                if q.reverse || cut >= g.lanes[target].length {
                    *next += 1;
                    self.counts.target_members += 1;
                    return Ok(());
                }
                let w = weight(plus(
                    bytes(c)?,
                    bytes(&self.snapshot.native[c.family].assignment)?,
                )?)?;
                self.mem.reserve(w)?; // entire native-tail chain reserved BEFORE clone/push
                let mut segments = c.segments.clone();
                segments.push(piece(c, hub.port.cut(g.k)?));
                segments.push(Segment {
                    source: target,
                    start: cut,
                    end: g.lanes[target].length,
                    reverse: false,
                });
                let mut all = c.completed.clone();
                all.push(Route { segments }.normalize());
                let assignment = Assignment {
                    version: routes::VERSION,
                    model: routes::MODEL.into(),
                    graph_checksum: self.snapshot.bindings.graph.clone(),
                    family: c.family,
                    routes: all,
                };
                self.chains.push_back(Chain {
                    task: p.task,
                    entry: p.next_entry.saturating_sub(1),
                    assignment,
                    stage: ChainStage::Tail { slot: c.slot + 1 },
                    weight: w,
                });
                *next += 1;
                self.counts.target_members += 1;
                self.counts.proposals += 1;
            }
        }
        Ok(())
    }
    fn close(&mut self, c: &mut Chain, e: &mut Evaluator<'_>) -> io::Result<bool> {
        match &mut c.stage {
            ChainStage::Tail { slot } => {
                let family = &e.graph.families[c.assignment.family];
                if *slot < family.paths.len() {
                    let source = family.paths[*slot];
                    c.assignment.routes.push(Route {
                        segments: vec![Segment {
                            source,
                            start: 0,
                            end: e.graph.lanes[source].length,
                            reverse: false,
                        }],
                    });
                    *slot += 1;
                } else {
                    c.stage = ChainStage::Capacity { i: 0, j: 1 };
                }
            }
            ChainStage::Capacity { i, j } => {
                let a = c.assignment.routes.iter().flat_map(|r| &r.segments).nth(*i);
                if let Some(a) = a {
                    if let Some(b) = c.assignment.routes.iter().flat_map(|r| &r.segments).nth(*j) {
                        self.counts.capacity_comparisons += 1;
                        if !disjoint(a, b) {
                            self.counts.capacity_rejections += 1;
                            return Ok(true);
                        }
                        *j += 1;
                    } else {
                        *i += 1;
                        *j = *i + 1;
                    }
                } else {
                    c.stage = ChainStage::Evaluate;
                }
            }
            ChainStage::Evaluate => {
                let native = &self.snapshot.native[c.assignment.family];
                let saved = if native.assignment == c.assignment {
                    Some((native.objective_bits, "native"))
                } else {
                    self.cache.get(&c.assignment).map(|b| (*b, "cache"))
                };
                if saved.is_none() {
                    ensure(
                        self.counts.fresh_evaluations < self.options.max_evaluations,
                        "evaluation-budget-exhausted",
                    )?;
                }
                let distinct = !self.cache.contains_key(&c.assignment);
                if distinct {
                    ensure(
                        self.cache.len() < self.options.max_distinct,
                        "distinct-budget-exhausted",
                    )?;
                    self.mem.reserve(weight(bytes(&c.assignment)?)?)?;
                }
                self.counts.evaluation_visits += 1;
                let (bits, kind) = if let Some(saved) = saved {
                    if saved.1 == "native" {
                        self.counts.native_reuses += 1
                    } else {
                        self.counts.cached_reuses += 1
                    }
                    saved
                } else {
                    // Atomic unchanged evaluator: validates whole topology/DNA/capacity and exact counts.
                    let result = e.evaluate(&c.assignment)?;
                    ensure(
                        result.assignment == c.assignment,
                        "proposal normalization mismatch",
                    )?;
                    self.counts.fresh_evaluations += 1;
                    (result.relative_objective.to_bits(), "fresh")
                };
                let score_tie = distinct && !self.scores.insert(bits);
                if distinct {
                    self.cache.insert(c.assignment.clone(), bits);
                    self.counts.distinct += 1;
                }
                if score_tie {
                    self.counts.score_ties += 1;
                }
                if !self.input_end {
                    self.counts.closures_before_input_end += 1;
                }
                serde_json::to_writer(&mut self.ledger,&json!({"visit":self.counts.evaluation_visits,"work":self.counts.work,"task":c.task,"entry":c.entry,"kind":kind,"distinct_normalized_physical_assignment":distinct,"score_tie":score_tie,"assignment":c.assignment,"objective_bits":bits,"relative_objective":f64::from_bits(bits),"native_delta_bits":self.snapshot.native.iter().map(|n|(f64::from_bits(bits)-n.objective()).to_bits()).collect::<Vec<_>>(),"input_inspected":self.counts.inspected,"discovery_unfinished":!self.input_end})).map_err(io::Error::other)?;
                self.ledger.write_all(b"\n")?;
                return Ok(true);
            }
        }
        Ok(false)
    }
    pub fn finish(mut self, position: u64) -> io::Result<()> {
        self.ledger.flush()?;
        create_json(
            &self.options.out_dir.join("result.json"),
            &json!({"method":"experimental-snapshot-completion-probe-v1","status":self.stop,"proposal_subset_only":true,"candidate_domain_complete":false,"original_frontier_immutable":true,"one_shot_no_sidecar_resume":true,"sequence_emission_authorized":false,"counts":self.counts,"native_baselines_reused":self.snapshot.native.len(),"native_initialization_evaluations":0,"input_position":position,"input_end":self.input_end,"input_tasks_preflight":self.snapshot.tasks,"pending_input":self.pending.as_ref().map(|p|json!({"task":p.task.id,"eligibility_position":p.index})),"pending_producers":self.producers.iter().map(|p|json!({"task":p.task,"next_entry":p.next_entry,"source_stage":p.source,"hub_stage":p.hub.as_ref().map(|h|&h.stage)})).collect::<Vec<_>>(),"pending_chains":self.chains.iter().map(|c|json!({"task":c.task,"stage":c.stage})).collect::<Vec<_>>(),"logical_reserved_bytes":self.mem.used,"peak_logical_reserved_bytes":self.mem.peak,"logical_limit_bytes":self.mem.limit,"output_fnv1a64":stream::hash_file(&self.options.out_dir.join("proposals.jsonl"))?}),
        )
    }
}
fn positive(c: &Context, p: &Port, k: u64) -> io::Result<bool> {
    ensure(c.source == p.source, "deferred source mismatch")?;
    let cut = p.cut(k)?;
    Ok(c.reverse == p.reverse && if c.reverse { cut < c.cut } else { cut > c.cut })
}
fn piece(c: &Context, to: u64) -> Segment {
    Segment {
        source: c.source,
        start: c.cut.min(to),
        end: c.cut.max(to),
        reverse: c.reverse,
    }
}
fn disjoint(a: &Segment, b: &Segment) -> bool {
    a.source != b.source || a.end <= b.start || b.end <= a.start
}
fn entry_context(original: &Context, i: usize) -> io::Result<Context> {
    ensure(i <= original.segments.len(), "entry cursor outside prefix")?;
    if i == original.segments.len() {
        return Ok(original.clone());
    }
    let s = &original.segments[i];
    Ok(Context {
        family: original.family,
        slot: original.slot,
        completed: original.completed.clone(),
        segments: original.segments[..i].to_vec(),
        source: s.source,
        cut: s.entry(),
        reverse: s.reverse,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn earlier_entry_discards_detour_and_keeps_arbitrary_multidonor_prefix() {
        let segments = (0..5)
            .map(|source| Segment {
                source,
                start: source as u64 * 10,
                end: source as u64 * 10 + 5,
                reverse: source == 3,
            })
            .collect::<Vec<_>>();
        let c = Context {
            family: 0,
            slot: 1,
            completed: vec![Route {
                segments: vec![segments[0].clone()],
            }],
            segments: segments.clone(),
            source: 5,
            cut: 50,
            reverse: false,
        };
        for i in 0..5 {
            let p = entry_context(&c, i).unwrap();
            assert_eq!(p.segments, segments[..i]);
            assert_eq!(p.completed, c.completed);
            assert_eq!(p.cut, segments[i].entry());
            assert_eq!(p.source, i);
        }
        assert_eq!(entry_context(&c, 3).unwrap().segments.len(), 3);
        let mut mem = Memory {
            used: 10,
            peak: 10,
            limit: 20,
        };
        assert!(mem.reserve(11).is_err());
        assert_eq!((mem.used, mem.peak), (10, 10));
        assert!(!disjoint(
            &segments[0],
            &Segment {
                reverse: true,
                ..segments[0].clone()
            }
        ));
    }
}
fn target_before(q: &Port, word: &Vec<u8>, target: usize, upper: bool) -> bool {
    let ord = (&q.word, q.source).cmp(&(word, target));
    ord.is_lt() || (upper && ord.is_eq())
}

/// Mechanism-only injected contexts, never used by the end-to-end gate/CLI.
#[cfg(test)]
pub(crate) fn mechanism_checks(mut o: Options, e: &mut Evaluator<'_>) -> io::Result<Value> {
    fs::create_dir(&o.out_dir)?;
    // Frozen mechanism budget independent of the meaningful mosaic gate.
    o.max_work = 20000;
    o.max_evaluations = 1000;
    o.producer_slots = 2;
    o.chain_slots = 2;
    create_json(
        &o.out_dir.join("mechanism-predeclaration.json"),
        &json!({"options":o,"injected_contexts_only":true,"requirements":["dead current tip","two earlier distinct donors retained","deferred positive prefix","closure before 10000 native distractors exhausted","full capacity rejection before evaluator"]}),
    )?;
    let g = e.graph;
    ensure(
        g.lanes.len() == 4 && g.families.len() == 4,
        "mechanism fixture shape",
    )?;
    let target = g.families[0].paths[0];
    let sources = g.families.iter().map(|f| f.paths[0]).collect::<Vec<_>>();
    let mut ports = vec![];
    for i in 0..g.port_count {
        ports.push(adapter::read(&mut e.ports.global, g.k, g.port_count, i)?);
    }
    let mut common = ports
        .iter()
        .filter(|p| {
            p.source == sources[3]
                && !p.reverse
                && sources.iter().all(|s| {
                    ports.iter().any(|q| {
                        q.source == *s && !q.reverse && q.word == p.word && q.anchor == p.anchor
                    })
                })
        })
        .cloned()
        .collect::<Vec<_>>();
    common.sort_by_key(|p| p.anchor);
    ensure(common.len() >= 4, "mechanism lacks shared ports")?;
    let cuts = [
        common[0].cut(g.k)?,
        common[common.len() / 2].cut(g.k)?,
        common.last().unwrap().cut(g.k)?,
    ];
    ensure(
        cuts.windows(2).all(|w| w[0] < w[1]),
        "mechanism needs separated cuts",
    )?;
    let dead = Context {
        family: 0,
        slot: 0,
        completed: vec![],
        segments: (0..3)
            .map(|i| Segment {
                source: sources[i],
                start: if i == 0 { 0 } else { cuts[i - 1] },
                end: cuts[i],
                reverse: false,
            })
            .collect(),
        source: sources[3],
        cut: cuts[2],
        reverse: false,
    };
    let dead_returns = ports
        .iter()
        .filter(|p| p.source == dead.source && !p.reverse && p.cut(g.k).unwrap() > dead.cut)
        .filter(|p| {
            ports
                .iter()
                .any(|q| q.source == target && !q.reverse && q.word == p.word)
        })
        .count();
    ensure(dead_returns == 0, "current tip unexpectedly has a return")?;
    // All saved prefix joins are actual full-word ports; earlier complete routes
    // are not fabricated, and the last source-entry port exists normally.
    for (i, s) in dead.segments.iter().enumerate() {
        ensure(
            ports
                .iter()
                .any(|p| p.source == s.source && !p.reverse && p.cut(g.k).unwrap() == cuts[i]),
            "undeclared injected prefix exit",
        )?;
    }
    let deferred = entry_context(&dead, 2)?;
    let port = ports
        .iter()
        .find(|p| p.source == deferred.source && !p.reverse && p.cut(g.k).unwrap() == cuts[2])
        .unwrap()
        .clone();
    let mut native = vec![];
    for f in 0..g.families.len() {
        let r = e.evaluate(&g.native_assignment(f)?)?;
        native.push(Scored {
            assignment: r.assignment,
            objective_bits: r.relative_objective.to_bits(),
        });
    }
    let bindings = snapshot::Bindings {
        policy: identity::policy_identity(),
        backend: routes::compiler_identity(),
        graph: g.digest()?,
        sample: "mechanism-only".into(),
        depth_bits: 10.0f64.to_bits(),
        background_bits: 0.1f64.to_bits(),
        tie_bits: 1e-9f64.to_bits(),
        max_feature_terms: e.max_terms,
        cache_terms: 1000000,
        count_policy: genome::COUNT_POLICY.into(),
    };
    let snapshot = snapshot::Snapshot {
        bindings,
        seal: Value::Null,
        parent: Value::Null,
        budgets: Value::Null,
        transition: Value::Null,
        order_update: Value::Null,
        native,
        tasks_position: 0,
        tasks: 10000,
        preflight_bytes: 0,
    };
    let mut mem = Memory {
        used: 0,
        peak: 0,
        limit: o.max_state_bytes,
    };
    mem.reserve(weight(o.record_bytes as u64)? + 131072)?;
    mem.reserve(weight(bytes(&snapshot.native)?)?)?;
    let mut engine = Engine::new(&o, snapshot, mem)?;
    for (id, context, direct) in [(1, dead, None), (2, deferred, Some(port))] {
        let w = weight(bytes(&context)? + 2 * g.k)?;
        engine.mem.reserve(w)?;
        engine.producers.push_back(Producer {
            task: id,
            original: context,
            next_entry: 0,
            current: None,
            source: SourceStage::Recover,
            hub: None,
            direct,
            weight: w,
        });
        engine.counts.admitted += 1;
    }
    engine.counts.deferred = 1;
    // Check an opposite-orientation canonical conflict across WHOLE routes, before evaluation.
    let mut conflict = engine.snapshot.native[0].assignment.clone();
    let mut other = conflict.routes[0].clone();
    other.segments[0].reverse = true;
    conflict.routes.push(other);
    let mut chain = Chain {
        task: 0,
        entry: 0,
        assignment: conflict,
        stage: ChainStage::Capacity { i: 0, j: 1 },
        weight: 0,
    };
    ensure(
        engine.close(&mut chain, e)? && engine.counts.fresh_evaluations == 0,
        "capacity conflict reached evaluator",
    )?;
    let input_path = o.out_dir.join("native-distractors.json");
    let mut f = File::create(&input_path)?;
    f.write_all(b"{")?;
    let native_context = Context {
        family: 0,
        slot: 0,
        completed: vec![],
        segments: vec![],
        source: target,
        cut: 0,
        reverse: false,
    };
    for id in 100..10100 {
        if id != 100 {
            f.write_all(b",")?;
        }
        serde_json::to_writer(&mut f, &id.to_string()).map_err(io::Error::other)?;
        f.write_all(b":")?;
        serde_json::to_writer(
            &mut f,
            &Task {
                id,
                ready: id,
                depth: 0,
                prev: None,
                next: None,
                context: native_context.clone(),
                op: Op::Start,
            },
        )
        .map_err(io::Error::other)?;
    }
    f.write_all(b"}")?;
    drop(f);
    let mut input = stream::Json::open(&input_path, o.record_bytes)?;
    input.expect(b'{')?;
    engine.drive(&mut input, e)?;
    ensure(
        engine.counts.closures_before_input_end > 0
            && !engine.input_end
            && engine.counts.skipped_native > 0
            && engine.counts.admitted == 2,
        "distractors blocked closure service or filled queue",
    )?;
    let recovered = engine.counts.recovered_earlier;
    engine.finish(input.position)?;
    let mut proposals = stream::Json::open(&o.out_dir.join("proposals.jsonl"), o.record_bytes)?;
    let len = fs::metadata(o.out_dir.join("proposals.jsonl"))?.len();
    let mut earlier = 0;
    let mut multidonor = 0;
    let mut deferred_visits = 0;
    loop {
        proposals.ws()?;
        if proposals.position == len {
            break;
        }
        let row: Value = proposals.value()?;
        let a: Assignment =
            serde_json::from_value(row["assignment"].clone()).map_err(io::Error::other)?;
        if row["task"] == 1 {
            earlier += 1;
        }
        if row["task"] == 2 {
            deferred_visits += 1;
        }
        let donors = a
            .routes
            .iter()
            .flat_map(|r| &r.segments)
            .filter(|s| s.source != target)
            .map(|s| s.source)
            .collect::<BTreeSet<_>>();
        if donors.len() >= 2 {
            multidonor += 1;
        }
    }
    ensure(
        earlier > 0 && multidonor > 0 && deferred_visits > 0 && recovered > 0,
        "missing earlier/deferred/multi-donor closure",
    )?;
    let evidence = json!({"injected_mechanism_not_end_to_end":true,"dead_tip_target_returns":dead_returns,"earlier_prefix_closures":earlier,"multidonor_closures":multidonor,"deferred_closures":deferred_visits,"capacity_conflict_pre_evaluator":true,"unfinished_discovery":true});
    create_json(&o.out_dir.join("mechanism-evidence.json"), &evidence)?;
    Ok(evidence)
}

#[cfg(test)]
mod numeric_tests {
    use super::*;
    #[test]
    fn target_bounds_compare_numeric_ids_across_little_endian_boundary() {
        let mut f = tempfile::tempfile().unwrap();
        for source in [1u64, 255, 256, 256, 257, 65536] {
            f.write_all(b"ACG").unwrap();
            f.write_all(&source.to_le_bytes()).unwrap();
            f.write_all(&0u64.to_le_bytes()).unwrap();
            f.write_all(&[0]).unwrap();
        }
        let mut bounds = vec![];
        for upper in [false, true] {
            let (mut lo, mut hi) = (0, 6);
            while lo < hi {
                let mid = lo + (hi - lo) / 2;
                let q = adapter::read(&mut f, 3, 6, mid).unwrap();
                if target_before(&q, &b"ACG".to_vec(), 256, upper) {
                    lo = mid + 1
                } else {
                    hi = mid
                }
            }
            bounds.push(lo);
        }
        assert_eq!(bounds, vec![2, 4]);
        assert!(
            256u64.to_le_bytes() < 255u64.to_le_bytes(),
            "test must distinguish numeric from byte-string order"
        );
    }
}

/// Some(false) inspected one native span/boundary; Some(true) found changed
/// geometry; None completed a native-only inspection. Storage splits are not
/// biological changes. Complete routes each reset at zero with their own extent.
fn eligibility(c: &Context, index: usize, native: impl Fn(usize) -> (usize, u64)) -> Option<bool> {
    let (target, _) = native(c.slot);
    if index == 0 && (c.source != target || c.reverse) {
        return Some(true);
    }
    let position = c
        .completed
        .iter()
        .enumerate()
        .flat_map(|(slot, r)| (0..=r.segments.len()).map(move |i| (slot, &r.segments, i, true)))
        .chain((0..c.segments.len()).map(|i| (c.slot, &c.segments, i, false)))
        .nth(index);
    if let Some((slot, segments, i, complete)) = position {
        let (source, extent) = native(slot);
        if let Some(s) = segments.get(i) {
            let start = if i == 0 { 0 } else { segments[i - 1].end };
            Some(s.source != source || s.reverse || s.start != start)
        } else {
            debug_assert!(complete);
            Some(segments.last().is_none_or(|s| s.end != extent))
        }
    } else {
        let end = c.segments.last().map_or(0, |s| s.end);
        (c.cut != end).then_some(true)
    }
}

#[cfg(test)]
mod eligibility_tests {
    use super::*;
    fn native(slot: usize) -> (usize, u64) {
        [(3, 100), (7, 60), (9, 90)][slot]
    }
    fn segment(source: usize, start: u64, end: u64) -> Segment {
        Segment {
            source,
            start,
            end,
            reverse: false,
        }
    }
    fn context() -> Context {
        Context {
            family: 0,
            slot: 2,
            completed: vec![
                Route {
                    segments: vec![segment(3, 0, 40), segment(3, 40, 100)],
                },
                Route {
                    segments: vec![segment(7, 0, 20), segment(7, 20, 60)],
                },
            ],
            segments: vec![segment(9, 0, 10), segment(9, 10, 30)],
            source: 9,
            cut: 30,
            reverse: false,
        }
    }
    fn first_change(c: &Context) -> Option<usize> {
        for i in 0..20 {
            match eligibility(c, i, native) {
                Some(true) => return Some(i),
                Some(false) => (),
                None => return None,
            }
        }
        panic!("eligibility did not terminate")
    }
    #[test]
    fn contiguous_same_source_storage_splits_remain_native_per_molecule() {
        let c = context();
        for i in 0..8 {
            assert_eq!(eligibility(&c, i, native), Some(false));
        }
        assert_eq!(eligibility(&c, 8, native), None);
        assert_eq!(first_change(&c), None);
    }
    #[test]
    fn same_source_gaps_jumps_live_cut_and_completed_extents_are_changed() {
        let mut c = context();
        c.segments[1].start = 11;
        assert_eq!(first_change(&c), Some(7));
        c = context();
        c.segments[1].start = 9;
        assert_eq!(first_change(&c), Some(7));
        c = context();
        c.cut = 31;
        assert_eq!(first_change(&c), Some(8));
        c = context();
        c.completed[1].segments[1].end = 59;
        assert_eq!(first_change(&c), Some(5));
        c = context();
        c.completed[1].segments[0].start = 1;
        assert_eq!(first_change(&c), Some(3));
        c = context();
        c.completed[1].segments[0].source = 3;
        assert_eq!(first_change(&c), Some(3));
        c = context();
        c.completed[1].segments.clear();
        assert_eq!(first_change(&c), Some(3));
    }
}
