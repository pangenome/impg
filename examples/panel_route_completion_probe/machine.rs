//! Alternating discovery/active service and producer/closure service, with backpressure.
use super::*;
use std::collections::LinkedList;
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
// Membership is owned here, not copied/indexed by repeated scans of the pool.
struct ReadyFamily {
    producers: LinkedList<Producer>,
    grants: u64,
    native_order: u64,
    in_flight: bool,
}
struct ReadyFamilies {
    families: Vec<ReadyFamily>,
    eligible: BTreeSet<(u64, u64, usize)>,
    producers: usize,
}
impl ReadyFamilies {
    fn new(native: &[Scored], slots: usize, mem: &mut Memory) -> io::Result<Self> {
        // One family row/key plus prepaid node/transfer slack per live producer.
        // LinkedList owns one node per queued element and pop_front deallocates
        // that node: empty/partly drained families retain no historical capacity.
        // 8*Producer+256 covers two links/alignment, a local moved Producer and
        // even one extra transient node under the total live producer width.
        // Producer payload allocations are separately reserved, never cloned here.
        let rows = (native.len() as u64)
            .checked_mul(256)
            .ok_or_else(|| invalid("ready metadata overflow"))?;
        let queues = (slots as u64)
            .checked_mul(8 * std::mem::size_of::<Producer>() as u64 + 256)
            .ok_or_else(|| invalid("ready queue overflow"))?;
        mem.reserve(plus(rows, queues)?)?;
        Ok(Self {
            families: native
                .iter()
                .map(|n| ReadyFamily {
                    producers: LinkedList::new(),
                    grants: 0,
                    in_flight: false,
                    // Monotone total order on finite saved f64 scores, including -0.
                    native_order: if n.objective_bits >> 63 == 1 {
                        !n.objective_bits
                    } else {
                        n.objective_bits ^ (1 << 63)
                    },
                })
                .collect(),
            eligible: BTreeSet::new(),
            producers: 0,
        })
    }
    fn key(&self, family: usize) -> (u64, u64, usize) {
        let f = &self.families[family];
        (f.grants, f.native_order, family)
    }
    fn refresh(&mut self, family: usize) {
        let f = &self.families[family];
        if !f.in_flight && !f.producers.is_empty() {
            self.eligible.insert(self.key(family));
        }
    }
    fn retain(&mut self, p: Producer) {
        let family = p.original.family;
        self.families[family].producers.push_back(p);
        self.producers += 1;
        self.refresh(family);
    }
    fn take(&mut self) -> io::Result<Producer> {
        let key = *self.eligible.first().expect("ready family available");
        let next = plus(key.0, 1)?;
        self.eligible.remove(&key);
        let f = &mut self.families[key.2];
        f.grants = next; // every granted member attempt, regardless of outcome
        self.producers -= 1;
        Ok(f.producers
            .pop_front()
            .expect("eligible family owns producer"))
    }
    fn launched(&mut self, family: usize) {
        assert!(!self.families[family].in_flight);
        self.families[family].in_flight = true;
    }
    fn retired(&mut self, family: usize) {
        assert!(self.families[family].in_flight);
        self.families[family].in_flight = false;
        self.refresh(family);
    }
}
fn member_ready(p: &Producer) -> bool {
    matches!(p.hub.as_ref().map(|h| &h.stage), Some(HubStage::Members { next, end }) if next < end)
}
pub struct Engine<'a> {
    options: &'a Options,
    pub snapshot: snapshot::Snapshot,
    mem: Memory,
    counts: Counts,
    producers: VecDeque<Producer>,
    ready: ReadyFamilies,
    member_turn: bool,
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
        let ready = ReadyFamilies::new(&snapshot.native, options.producer_slots, &mut mem)?;
        Ok(Self {
            ready,
            member_turn: true,
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
            let discovery = !self.input_end
                && self.producers.len() + self.ready.producers < self.options.producer_slots;
            let active =
                !self.producers.is_empty() || self.ready.producers > 0 || !self.chains.is_empty();
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
    fn retain_producer(&mut self, p: Producer) {
        if member_ready(&p) {
            self.ready.retain(p);
        } else {
            self.producers.push_back(p);
        }
    }
    fn active(&mut self, e: &mut Evaluator<'_>) -> io::Result<()> {
        let ready = self.chains.len() < self.options.chain_slots && !self.ready.eligible.is_empty();
        let producing = !self.producers.is_empty() || ready;
        if producing && (self.producer_turn || self.chains.is_empty()) {
            self.producer_turn = false;
            if ready && (self.member_turn || self.producers.is_empty()) {
                self.member_turn = false;
                let mut p = self.ready.take()?;
                let family = p.original.family;
                let before = self.chains.len();
                let result = self.member(&mut p, e);
                if self.chains.len() > before {
                    self.ready.launched(family);
                }
                self.retain_producer(p);
                self.ready.refresh(family);
                result
            } else {
                self.member_turn = true;
                let mut p = self.producers.pop_front().unwrap();
                let result = self.produce(&mut p, e);
                if matches!(result, Ok(true)) {
                    self.mem.release(p.weight);
                } else {
                    self.retain_producer(p);
                }
                result.map(|_| ())
            }
        } else {
            self.producer_turn = true;
            let mut c = self.chains.pop_front().unwrap();
            let result = self.close(&mut c, e);
            if matches!(result, Ok(true)) {
                self.mem.release(c.weight);
                self.ready.retired(c.assignment.family);
            } else {
                // In particular, a resource-stopped chain still owns its family.
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
            &json!({"method":METHOD,"ready_family_service":self.ready.families.iter().enumerate().map(|(family,f)|json!({"family":family,"member_attempt_grants":f.grants,"in_flight":f.in_flight,"ready_producers":f.producers.len()})).collect::<Vec<_>>(),"status":self.stop,"proposal_subset_only":true,"candidate_domain_complete":false,"original_frontier_immutable":true,"one_shot_no_sidecar_resume":true,"sequence_emission_authorized":false,"counts":self.counts,"native_baselines_reused":self.snapshot.native.len(),"native_initialization_evaluations":0,"input_position":position,"input_end":self.input_end,"input_tasks_preflight":self.snapshot.tasks,"pending_input":self.pending.as_ref().map(|p|json!({"task":p.task.id,"eligibility_position":p.index})),"pending_producers":self.producers.iter().chain(self.ready.families.iter().flat_map(|f|f.producers.iter())).map(|p|json!({"task":p.task,"family":p.original.family,"member_ready":member_ready(p),"next_entry":p.next_entry,"source_stage":p.source,"hub_stage":p.hub.as_ref().map(|h|&h.stage)})).collect::<Vec<_>>(),"pending_chains":self.chains.iter().map(|c|json!({"task":c.task,"family":c.assignment.family,"stage":c.stage})).collect::<Vec<_>>(),"logical_reserved_bytes":self.mem.used,"peak_logical_reserved_bytes":self.mem.peak,"logical_limit_bytes":self.mem.limit,"output_fnv1a64":stream::hash_file(&self.options.out_dir.join("proposals.jsonl"))?}),
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

#[cfg(test)]
mod ready_tests {
    use super::*;
    fn native(bits: &[f64]) -> Vec<Scored> {
        bits.iter()
            .enumerate()
            .map(|(family, n)| Scored {
                assignment: Assignment {
                    version: 1,
                    model: routes::MODEL.into(),
                    graph_checksum: String::new(),
                    family,
                    routes: vec![],
                },
                objective_bits: n.to_bits(),
            })
            .collect()
    }
    fn producer(family: usize, task: u64) -> Producer {
        Producer {
            task,
            original: Context {
                family,
                slot: 0,
                completed: vec![],
                segments: vec![],
                source: 0,
                cut: 0,
                reverse: false,
            },
            next_entry: 0,
            current: None,
            source: SourceStage::Recover,
            hub: Some(Hub {
                port: Port {
                    word: vec![],
                    source: 0,
                    anchor: 0,
                    reverse: false,
                },
                stage: HubStage::Members { next: 0, end: 2 },
            }),
            direct: None,
            weight: 0,
        }
    }
    fn check(r: &ReadyFamilies) {
        let expected = r
            .families
            .iter()
            .enumerate()
            .filter(|(_, f)| !f.in_flight && !f.producers.is_empty())
            .map(|(i, _)| r.key(i))
            .collect::<BTreeSet<_>>();
        assert_eq!(r.eligible, expected);
        assert_eq!(
            r.producers,
            r.families.iter().map(|f| f.producers.len()).sum::<usize>()
        );
    }
    #[test]
    fn ready_family_local_keys_fifo_unready_and_all_outcome_grants() {
        let mut mem = Memory {
            used: 0,
            peak: 0,
            limit: 1 << 20,
        };
        let mut r = ReadyFamilies::new(&native(&[-100.0, -10.0, -1.0]), 16, &mut mem).unwrap();
        // Best family 0 is absent and cannot impose an all-family barrier.
        r.retain(producer(1, 11));
        r.retain(producer(1, 12));
        r.retain(producer(2, 21));
        check(&r);
        let a = r.take().unwrap();
        assert_eq!(a.task, 11);
        r.launched(1);
        r.refresh(1);
        check(&r);
        let mut b = r.take().unwrap();
        assert_eq!(b.task, 21); // A is better-ranked but already in flight.
                                // A rejected member consumes the same grant, advances only its own cursor.
        if let Some(Hub {
            stage: HubStage::Members { next, .. },
            ..
        }) = &mut b.hub
        {
            *next += 1;
        }
        r.retain(b);
        check(&r);
        assert_eq!(r.families[2].grants, 1);
        let b = r.take().unwrap();
        assert_eq!(b.task, 21);
        r.launched(2);
        r.refresh(2);
        check(&r);
        // Cache/native/capacity terminal outcomes all retire identically; no bonus
        // or reset for a non-fresh evaluation or a nonproductive closure.
        r.retired(1);
        check(&r);
        let a = r.take().unwrap();
        assert_eq!(a.task, 12);
        r.launched(1);
        r.refresh(1);
        assert_eq!(r.families[1].grants, 2);
        r.retired(1);
        r.retired(2);
        check(&r);
        assert!(r.eligible.is_empty());
        r.retain(producer(1, 13));
        check(&r);
        assert_eq!(r.take().unwrap().task, 13);
        r.refresh(1);
        check(&r);
        assert_eq!(r.families[1].grants, 3);
    }
    #[test]
    fn ready_family_storage_churn_releases_nodes_not_historical_capacity() {
        const FAMILIES: usize = 32;
        const WIDTH: usize = 16;
        const PASSES: usize = 4;
        const ROTATIONS: usize = 128;
        let p_bytes = std::mem::size_of::<Producer>();
        let alignment = std::mem::align_of::<Producer>().max(std::mem::align_of::<usize>());
        // std LinkedList owns one separately allocated node per element, with
        // two links; pop_front deallocates the removed node, not a reusable
        // per-family backing buffer. Allow two alignment pads conservatively.
        // This is owned collection storage, NOT allocator caching/RSS.
        let node_upper = p_bytes + 2 * std::mem::size_of::<usize>() + 2 * alignment;
        let per_slot = 8 * p_bytes + 256;
        assert!(2 * node_upper + p_bytes <= per_slot);
        let mut mem = Memory {
            used: 0,
            peak: 0,
            limit: 1 << 20,
        };
        let mut ready = ReadyFamilies::new(&native(&[-1.0; FAMILIES]), WIDTH, &mut mem).unwrap();
        let reservation = mem.used as usize;
        assert_eq!(reservation, 256 * FAMILIES + WIDTH * per_slot);
        let mut old: Vec<VecDeque<Producer>> = (0..FAMILIES).map(|_| VecDeque::new()).collect();
        let mut task = 0;
        let mut peak_nodes = 0;
        let mut peak_bound = 0;
        let mut partial = false;
        let mut most_families = 0;
        let mut observe = |r: &ReadyFamilies, old: &[VecDeque<Producer>]| {
            check(r);
            let nodes = r.families.iter().map(|f| f.producers.len()).sum::<usize>();
            assert!(nodes <= WIDTH);
            assert_eq!(nodes, old.iter().map(VecDeque::len).sum::<usize>());
            for (f, q) in r.families.iter().zip(old) {
                // Same FIFO contents throughout the adversarial storage history.
                assert_eq!(
                    f.producers.iter().map(|p| p.task).collect::<Vec<_>>(),
                    q.iter().map(|p| p.task).collect::<Vec<_>>()
                );
            }
            // LinkedList length equals owned nodes specifically because removed
            // nodes are deallocated. VecDeque length would NOT justify this bound.
            let bound = (nodes + 1) * node_upper + p_bytes;
            assert!(bound <= WIDTH * per_slot);
            peak_nodes = peak_nodes.max(nodes);
            peak_bound = peak_bound.max(bound);
            partial |= r
                .families
                .iter()
                .any(|f| (2..WIDTH).contains(&f.producers.len()));
            most_families = most_families.max(
                r.families
                    .iter()
                    .filter(|f| !f.producers.is_empty())
                    .count(),
            );
        };
        for _ in 0..PASSES {
            for family in 0..FAMILIES {
                for _ in 0..WIDTH {
                    task += 1;
                    ready.retain(producer(family, task));
                    old[family].push_back(producer(family, task));
                    observe(&ready, &old);
                }
                for _ in 0..WIDTH {
                    let p = ready.take().unwrap();
                    assert_eq!(p.original.family, family);
                    assert_eq!(old[family].pop_front().unwrap().task, p.task);
                    ready.refresh(family);
                    observe(&ready, &old);
                }
            }
        }
        let old_capacity_slots = old.iter().map(VecDeque::capacity).sum::<usize>();
        let old_retained_bytes = old_capacity_slots * p_bytes;
        assert!(old.iter().all(VecDeque::is_empty));
        assert!(old_capacity_slots >= FAMILIES * WIDTH);
        // Reproduce the P1 against ACTUAL old capacities, not live lengths.
        assert!(old_retained_bytes > reservation);
        assert_eq!(ready.producers, 0);
        // Keep a partially drained large family alongside a moving small family.
        for _ in 0..WIDTH {
            task += 1;
            ready.retain(producer(0, task));
            old[0].push_back(producer(0, task));
        }
        for _ in 0..FAMILIES {
            let mut p = ready.take().unwrap();
            let from = p.original.family;
            let mut q = old[from].pop_front().unwrap();
            ready.refresh(from);
            observe(&ready, &old);
            p.original.family = (from + 1) % FAMILIES;
            q.original.family = p.original.family;
            old[p.original.family].push_back(q);
            ready.retain(p);
            observe(&ready, &old);
        }
        while ready.producers > 0 {
            let p = ready.take().unwrap();
            old[p.original.family].pop_front().unwrap();
            ready.refresh(p.original.family);
            observe(&ready, &old);
        }
        // Sixteen nonempty families coexist; move one element at a time without
        // ever increasing the aggregate width, including between earlier families.
        for family in 0..WIDTH {
            task += 1;
            ready.retain(producer(family, task));
            old[family].push_back(producer(family, task));
        }
        observe(&ready, &old);
        for _ in FAMILIES..ROTATIONS {
            let mut p = ready.take().unwrap();
            let from = p.original.family;
            let mut q = old[from].pop_front().unwrap();
            ready.refresh(from);
            observe(&ready, &old);
            p.original.family = (from + WIDTH) % FAMILIES;
            q.original.family = p.original.family;
            old[p.original.family].push_back(q);
            ready.retain(p);
            observe(&ready, &old);
        }
        while ready.producers > 0 {
            let p = ready.take().unwrap();
            old[p.original.family].pop_front().unwrap();
            ready.refresh(p.original.family);
            observe(&ready, &old);
        }
        drop(observe);
        assert!(partial);
        assert_eq!(most_families, WIDTH);
        assert_eq!(peak_nodes, WIDTH);
        assert_eq!(mem.used as usize, reservation);
        assert!(ready.families.iter().all(|f| f.producers.is_empty()));
        let evidence = json!({"families":FAMILIES,"width":WIDTH,"passes":PASSES,"rotations":ROTATIONS,
            "producer_bytes":p_bytes,"node_upper_bytes":node_upper,"prepaid_ready_bytes":reservation,
            "old_retained_capacity_slots_after_full_drain":old_capacity_slots,"old_backing_bytes_after_full_drain":old_retained_bytes,
            "old_reservation_violation":old_retained_bytes>reservation,"new_peak_owned_nodes":peak_nodes,
            "new_peak_nodes_and_transfer_bound":peak_bound,"node_allowance":WIDTH*per_slot,
            "new_owned_nodes_after_full_drain":0,"partially_drained_families_checked":partial,
            "max_simultaneous_nonempty_families":most_families,"fifo_and_key_conservation":true,
            "ownership_basis":"std LinkedList one owned node per element; pop_front deallocates nodes; no historical per-family backing capacity; excludes allocator caching/RSS"});
        if let Some(path) = std::env::var_os("IMPG_TEST_READY_STORAGE_OUTPUT") {
            create_json(Path::new(&path), &evidence).unwrap();
        }
        println!("{evidence}");
    }
    #[test]
    fn ready_metadata_reservation_stops_before_growth() {
        let mut mem = Memory {
            used: 7,
            peak: 7,
            limit: 8,
        };
        assert!(ReadyFamilies::new(&native(&[-1.0]), 16, &mut mem).is_err());
        assert_eq!((mem.used, mem.peak), (7, 7));
    }
}

#[cfg(test)]
fn fixture_engine<'a>(
    o: &'a Options,
    native: &[Scored],
    g: &routes::Graph,
) -> io::Result<Engine<'a>> {
    fs::create_dir(&o.out_dir)?;
    let snapshot = snapshot::Snapshot {
        bindings: snapshot::Bindings {
            policy: identity::policy_identity(),
            backend: routes::compiler_identity(),
            graph: g.digest()?,
            sample: "ready-family-mechanism-only".into(),
            depth_bits: 10.0f64.to_bits(),
            background_bits: 0.1f64.to_bits(),
            tie_bits: 1e-9f64.to_bits(),
            max_feature_terms: 50000000,
            cache_terms: 1000000,
            count_policy: genome::COUNT_POLICY.into(),
        },
        seal: Value::Null,
        parent: Value::Null,
        budgets: Value::Null,
        transition: Value::Null,
        order_update: Value::Null,
        native: native.to_vec(),
        tasks_position: 0,
        tasks: 10000,
        preflight_bytes: 0,
    };
    let mut mem = Memory {
        used: 0,
        peak: 0,
        limit: o.max_state_bytes,
    };
    mem.reserve(weight(o.record_bytes as u64)? + 131072 + weight(bytes(&snapshot.native)?)?)?;
    Engine::new(o, snapshot, mem)
}
#[cfg(test)]
fn conservation(e: &Engine<'_>) {
    let mut ids = BTreeSet::new();
    for p in &e.producers {
        assert!(!member_ready(p));
        assert!(ids.insert(p.task));
    }
    for (family, f) in e.ready.families.iter().enumerate() {
        for p in &f.producers {
            assert!(member_ready(p));
            assert_eq!(p.original.family, family);
            assert!(ids.insert(p.task));
        }
        assert_eq!(
            e.chains
                .iter()
                .filter(|c| c.assignment.family == family)
                .count(),
            usize::from(f.in_flight)
        );
        assert_eq!(
            e.ready.eligible.contains(&e.ready.key(family)),
            !f.in_flight && !f.producers.is_empty()
        );
    }
    assert_eq!(ids.len(), e.producers.len() + e.ready.producers);
    assert!(ids.len() <= e.options.producer_slots);
    assert!(e.chains.len() <= e.options.chain_slots);
}

/// Isolated bb981c8 selector only. Uses the identical member/tail/capacity/evaluator
/// primitives, with the original producer RR and shared chain pool; no fair index.
#[cfg(test)]
fn legacy_drive(
    engine: &mut Engine<'_>,
    input: &mut stream::Json,
    e: &mut Evaluator<'_>,
) -> io::Result<()> {
    loop {
        let discovery = !engine.input_end && engine.producers.len() < engine.options.producer_slots;
        let active = !engine.producers.is_empty() || !engine.chains.is_empty();
        if !discovery && !active {
            engine.stop = "proposal-subset-exhausted".into();
            break;
        }
        if engine.counts.work >= engine.options.max_work {
            engine.stop = "work-budget-exhausted".into();
            break;
        }
        engine.counts.work += 1;
        let result = if discovery && (engine.discovery_turn || !active) {
            engine.discovery_turn = false;
            engine.discovery(input, e)
        } else {
            engine.discovery_turn = true;
            if !engine.producers.is_empty() && (engine.producer_turn || engine.chains.is_empty()) {
                engine.producer_turn = false;
                let mut p = engine.producers.pop_front().unwrap();
                let result = engine.produce(&mut p, e);
                if matches!(result, Ok(true)) {
                    engine.mem.release(p.weight);
                } else {
                    engine.producers.push_back(p);
                }
                result.map(|_| ())
            } else {
                engine.producer_turn = true;
                let mut c = engine.chains.pop_front().unwrap();
                let result = engine.close(&mut c, e);
                if matches!(result, Ok(true)) {
                    engine.mem.release(c.weight);
                } else {
                    engine.chains.push_back(c);
                }
                result.map(|_| ())
            }
        };
        if let Err(err) = result {
            if err.to_string() == "evaluation-budget-exhausted" {
                engine.stop = err.to_string();
                break;
            }
            return Err(err);
        }
    }
    engine.ledger.flush()?;
    engine.ledger.sync_all()
}

#[cfg(test)]
pub(crate) fn ready_family_checks(mut o: Options, e: &mut Evaluator<'_>) -> io::Result<Value> {
    fs::create_dir(&o.out_dir)?;
    o.max_work = 20000;
    o.max_evaluations = 2; // frozen mechanism budget; widths/defaults unchanged
    let root = o.out_dir.clone();
    create_json(
        &root.join("predeclaration.json"),
        &json!({"options":o,"selector_fixture":"already-ready A,A,D; A native strictly better than D; live B source and10000 undiscovered native contexts","outcome_fixture":"one attempt at a time; fresh cap1,work2000: fresh/cache/native/capacity-reject/member-reject","before":"isolated bb981c8 producer-RR selection, unchanged primitives"}),
    )?;
    let g = e.graph;
    ensure(g.families.len() == 4, "four-family mechanism required")?;
    let sources = g.families.iter().map(|f| f.paths[0]).collect::<Vec<_>>();
    let mut ports = vec![];
    for i in 0..g.port_count {
        ports.push(adapter::read(&mut e.ports.global, g.k, g.port_count, i)?);
    }
    let mut common = ports
        .iter()
        .filter(|p| {
            p.source == sources[2]
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
    ensure(common.len() >= 4, "separated common ports required")?;
    let entry = common[0].cut(g.k)?;
    let exits = [
        common[common.len() / 2].clone(),
        common.last().unwrap().clone(),
    ];
    let mut native = vec![];
    for f in 0..4 {
        let r = e.evaluate(&g.native_assignment(f)?)?;
        native.push(Scored {
            assignment: r.assignment,
            objective_bits: r.relative_objective.to_bits(),
        });
    }
    ensure(
        native[0].objective() < native[3].objective(),
        "waiting family must have worse native rank",
    )?;
    let make = |family: usize, task: u64, exit: usize| -> io::Result<Producer> {
        let source = sources[2];
        let target = sources[family];
        let port = exits[exit].clone();
        let member = ports
            .iter()
            .position(|q| q.source == target && !q.reverse && q.word == port.word)
            .ok_or_else(|| invalid("missing test return"))? as u64;
        let c = Context {
            family,
            slot: 0,
            completed: vec![],
            segments: vec![Segment {
                source: target,
                start: 0,
                end: entry,
                reverse: false,
            }],
            source,
            cut: entry,
            reverse: false,
        };
        Ok(Producer {
            task,
            weight: weight(bytes(&c)? + 2 * g.k)?,
            original: c.clone(),
            next_entry: 2,
            current: Some(c),
            source: SourceStage::Recover,
            hub: Some(Hub {
                port,
                stage: HubStage::Members {
                    next: member,
                    end: member + 1,
                },
            }),
            direct: None,
        })
    };
    let input_path = root.join("distractors.json");
    let mut input_file = File::create(&input_path)?;
    input_file.write_all(b"{")?;
    for id in 100..10100 {
        if id != 100 {
            input_file.write_all(b",")?;
        }
        serde_json::to_writer(&mut input_file, &id.to_string()).map_err(io::Error::other)?;
        input_file.write_all(b":")?;
        serde_json::to_writer(
            &mut input_file,
            &Task {
                id,
                ready: id,
                depth: 0,
                prev: None,
                next: None,
                context: Context {
                    family: 0,
                    slot: 0,
                    completed: vec![],
                    segments: vec![],
                    source: sources[0],
                    cut: 0,
                    reverse: false,
                },
                op: Op::Start,
            },
        )
        .map_err(io::Error::other)?;
    }
    input_file.write_all(b"}")?;
    drop(input_file);
    let mut summaries = vec![];
    for legacy in [true, false] {
        o.out_dir = root.join(if legacy {
            "before-old-selector"
        } else {
            "after-ready-family"
        });
        let mut engine = fixture_engine(&o, &native, g)?;
        for p in [make(0, 1, 0)?, make(0, 2, 1)?, make(3, 3, 0)?] {
            engine.mem.reserve(p.weight)?;
            if legacy {
                engine.producers.push_back(p)
            } else {
                engine.retain_producer(p)
            }
            engine.counts.admitted += 1;
        }
        let mut live = make(1, 4, 0)?;
        live.hub = None;
        live.source = SourceStage::Bound {
            lo: 0,
            hi: g.lanes[live.original.source].port_count,
        };
        engine.mem.reserve(live.weight)?;
        engine.producers.push_back(live);
        engine.counts.admitted += 1;
        let mut input = stream::Json::open(&input_path, o.record_bytes)?;
        input.expect(b'{')?;
        if legacy {
            legacy_drive(&mut engine, &mut input, e)?
        } else {
            engine.drive(&mut input, e)?;
            conservation(&engine);
        }
        let ledger = std::fs::read_to_string(o.out_dir.join("proposals.jsonl"))?;
        let families = ledger
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).unwrap())
            .filter(|r| r["kind"] == "fresh")
            .map(|r| r["assignment"]["family"].as_u64().unwrap())
            .collect::<Vec<_>>();
        let live = engine
            .producers
            .iter()
            .chain(engine.ready.families.iter().flat_map(|f| &f.producers))
            .any(|p| p.task == 4);
        let row = json!({"selector":if legacy{"isolated-bb981c8-RR"}else{METHOD},"fresh_families":families,"counts":engine.counts,"status":engine.stop,"source_producer_live":live,"discovery_unfinished":!engine.input_end,"input_position":input.position});
        create_json(&o.out_dir.join("causal-result.json"), &row)?;
        ensure(
            engine.stop == "evaluation-budget-exhausted"
                && families.len() == 2
                && live
                && !engine.input_end
                && engine.counts.inspected > 0,
            "causal fixture did not retain live work",
        )?;
        if legacy {
            ensure(
                families == [0, 0],
                "old selector did not reproduce waiting-family miss",
            )?;
        } else {
            ensure(
                families.contains(&3),
                "ready-family selector missed waiting family",
            )?;
            ensure(
                engine
                    .chains
                    .iter()
                    .all(|c| engine.ready.families[c.assignment.family].in_flight),
                "stopped chain lost family ownership",
            )?;
            engine.finish(input.position)?;
        }
        summaries.push(row);
    }
    // A failed chain reservation retains the exact member and ready ownership.
    o.out_dir = root.join("ready-storage-stop");
    let mut stopped = fixture_engine(&o, &native, g)?;
    let p = make(0, 9, 0)?;
    let original_hub = bytes(&p.hub)?;
    let cursor = serde_json::to_value(&p.hub).map_err(io::Error::other)?;
    stopped.mem.reserve(p.weight)?;
    stopped.retain_producer(p);
    stopped.mem.limit = stopped.mem.used;
    let occupied = stopped.mem.used;
    stopped.counts.work += 1;
    let err = stopped.active(e).unwrap_err();
    ensure(
        err.to_string() == "state-budget-exhausted",
        "wrong ready reservation stop",
    )?;
    conservation(&stopped);
    ensure(
        stopped.mem.used == occupied
            && stopped.ready.families[0].grants == 1
            && stopped.counts.target_members == 0
            && stopped.chains.is_empty(),
        "ready stop changed occupancy/cursor service",
    )?;
    let retained = stopped.ready.families[0].producers.front().unwrap();
    ensure(
        bytes(&retained.hub)? == original_hub
            && serde_json::to_value(&retained.hub).map_err(io::Error::other)? == cursor,
        "ready storage stop advanced member",
    )?;
    stopped.stop = err.to_string();
    stopped.finish(0)?;
    // Same service counter for rejected members and rejected/native/cached chains.
    o.out_dir = root.join("outcome-accounting");
    o.max_evaluations = 1;
    o.max_work = 2000;
    let mut engine = fixture_engine(&o, &native, g)?;
    for (i, kind) in ["fresh", "cache", "native", "capacity", "member"]
        .iter()
        .enumerate()
    {
        let mut p = make(0, 10 + i as u64, 0)?;
        if *kind == "native" {
            p.original.source = sources[0];
            p.current.as_mut().unwrap().source = sources[0];
            p.hub.as_mut().unwrap().port.source = sources[0];
        }
        if *kind == "capacity" {
            p.current.as_mut().unwrap().segments.push(Segment {
                source: sources[2],
                start: entry,
                end: exits[0].cut(g.k)?,
                reverse: false,
            });
            p.original = p.current.as_ref().unwrap().clone();
            // This explicitly one-shot outcome fixture has no unvisited entry.
            p.next_entry = p.original.segments.len() + 1;
            p.weight = weight(bytes(&p.original)? + 2 * g.k)?;
        }
        if *kind == "member" {
            let (index, q) = ports
                .iter()
                .enumerate()
                .find(|(_, q)| q.source == sources[0] && q.reverse)
                .unwrap();
            let h = p.hub.as_mut().unwrap();
            h.port.word = q.word.clone();
            h.stage = HubStage::Members {
                next: index as u64,
                end: index as u64 + 1,
            };
        }
        engine.mem.reserve(p.weight)?;
        engine.retain_producer(p);
        for _ in 0..400 {
            if engine.producers.is_empty()
                && engine.ready.producers == 0
                && engine.chains.is_empty()
            {
                break;
            }
            engine.counts.work += 1;
            engine.active(e)?;
            conservation(&engine);
        }
        ensure(
            engine.producers.is_empty() && engine.ready.producers == 0 && engine.chains.is_empty(),
            "outcome service did not drain",
        )?;
        create_json(
            &o.out_dir.join(format!("{i}-{kind}.json")),
            &json!({"counts":engine.counts,"grants":engine.ready.families[0].grants,"expected_grants":i+1}),
        )?;
        ensure(
            engine.ready.families[0].grants == i as u64 + 1,
            "outcome changed grant weight",
        )?;
    }
    ensure(
        engine.counts.fresh_evaluations == 1
            && engine.counts.cached_reuses == 1
            && engine.counts.native_reuses == 1
            && engine.counts.capacity_rejections == 1
            && engine.counts.target_members == 5,
        "outcome accounting mismatch",
    )?;
    create_json(
        &o.out_dir.join("outcome-evidence.json"),
        &json!({"counts":engine.counts,"grants":engine.ready.families[0].grants,"conservation":true}),
    )?;
    let evidence = json!({"injected_mechanism_only":true,"native_A_better_than_waiting_D":true,"before_after":summaries,"all_five_outcomes_charged_equally":true});
    create_json(&root.join("causal-evidence.json"), &evidence)?;
    Ok(evidence)
}
