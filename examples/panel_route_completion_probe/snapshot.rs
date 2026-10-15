//! Trusted-local immutable v3 input, with bounded syntax and native-ledger checks.
use super::*;
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Bindings {
    pub policy: BTreeMap<String, String>,
    pub backend: String,
    pub graph: String,
    pub sample: String,
    pub depth_bits: u64,
    pub background_bits: u64,
    pub tie_bits: u64,
    pub max_feature_terms: usize,
    pub cache_terms: usize,
    pub count_policy: String,
}
pub struct Snapshot {
    pub bindings: Bindings,
    pub seal: Value,
    pub parent: Value,
    pub budgets: Value,
    pub transition: Value,
    pub order_update: Value,
    pub native: Vec<Scored>,
    pub tasks_position: u64,
    pub tasks: u64,
    pub preflight_bytes: u64,
}
fn unique(seen: &mut BTreeSet<String>, key: &str, allowed: &[&str]) -> io::Result<()> {
    ensure(
        allowed.contains(&key) && seen.insert(key.into()),
        "unknown/duplicate snapshot field",
    )
}
pub fn task_valid(t: &Task, g: &routes::Graph) -> io::Result<()> {
    let c = &t.context;
    ensure(
        c.family < g.families.len()
            && c.slot < g.families[c.family].paths.len()
            && c.completed.len() == c.slot
            && c.source < g.lanes.len()
            && c.cut <= g.lanes[c.source].length,
        "invalid snapshot context",
    )?;
    for p in c
        .completed
        .iter()
        .flat_map(|r| &r.segments)
        .chain(&c.segments)
    {
        ensure(
            p.source < g.lanes.len() && p.start < p.end && p.end <= g.lanes[p.source].length,
            "invalid committed segment",
        )?;
    }
    match &t.op {
        Op::SourceScan { base, permutation } => {
            ensure(permutation.paired(), "unconverted v2 source cursor")?;
            permutation.validate(*base, g.lanes[c.source].port_count)?;
        }
        Op::SourceBound { lo, hi } => ensure(
            lo <= hi && *hi <= g.lanes[c.source].port_count,
            "invalid source bound",
        )?,
        Op::HubBound {
            port,
            lower,
            lo,
            hi,
        } => {
            adapter::verified(port, g)?;
            ensure(
                lo <= hi && *hi <= g.port_count && lower.is_none_or(|b| b <= *lo),
                "invalid hub bound",
            )?;
        }
        Op::HubScan {
            port,
            base,
            permutation,
        } => {
            adapter::verified(port, g)?;
            permutation.validate(*base, g.port_count)?;
        }
        _ => (),
    }
    Ok(())
}
pub fn preflight(
    o: &Options,
    g: &routes::Graph,
    sample: &str,
    mem: &mut Memory,
) -> io::Result<Snapshot> {
    mem.reserve(16384)?; // metadata indices/counters, before any root/state key growth
    let path = o.snapshot.join("checkpoint.json");
    let mut seal_input =
        stream::Json::open(&o.snapshot.join("checkpoint-seal.json"), o.record_bytes)?;
    let seal: Value = seal_input.value()?;
    seal_input.end()?;
    ensure(
        seal["version"] == 1 && seal.as_object().is_some_and(|m| m.len() == 3),
        "invalid checkpoint seal",
    )?;
    ensure(
        seal["checkpoint"] == stream::hash_file(&path)?
            && seal["ledger"] == stream::hash_file(&o.snapshot.join("evaluations.jsonl"))?,
        "checkpoint/ledger seal mismatch",
    )?;
    let mut j = stream::Json::open(&path, o.record_bytes)?;
    j.expect(b'{')?;
    let mut first = true;
    let mut seen = BTreeSet::new();
    let mut bindings = None;
    let mut version = 0;
    let mut native = Vec::new();
    let mut tasks_position = None;
    let mut tasks = 0;
    let mut parent = Value::Null;
    let mut transition = Value::Null;
    let mut order_update = Value::Null;
    let mut budgets = Value::Null;
    let mut state_visits = None;
    let mut state_evals = None;
    while j.next(&mut first, b'}')? {
        let key = j.key()?;
        unique(
            &mut seen,
            &key,
            &[
                "version",
                "bindings",
                "parent",
                "budgets",
                "state",
                "transition",
                "order_update",
            ],
        )?;
        match key.as_str() {
            "version" => version = j.value::<u32>()?,
            "bindings" => {
                let raw = j.raw()?;
                mem.reserve(weight(raw.len() as u64)?)?;
                bindings =
                    Some(serde_json::from_slice::<Bindings>(&raw).map_err(io::Error::other)?);
            }
            "parent" | "transition" | "order_update" | "budgets" => {
                let raw = j.raw()?;
                mem.reserve(weight(raw.len() as u64)?)?;
                let v = serde_json::from_slice(&raw).map_err(io::Error::other)?;
                match key.as_str() {
                    "parent" => parent = v,
                    "transition" => transition = v,
                    "budgets" => budgets = v,
                    _ => order_update = v,
                }
            }
            "state" => {
                j.expect(b'{')?;
                let mut first = true;
                let mut fields = BTreeSet::new();
                while j.next(&mut first, b'}')? {
                    let key = j.key()?;
                    unique(
                        &mut fields,
                        &key,
                        &[
                            "native",
                            "ranking",
                            "family_cursor",
                            "quantum",
                            "tasks",
                            "fifo",
                            "shallow",
                            "focus",
                            "historical_accounting",
                            "next_id",
                            "next_ready",
                            "task_bytes",
                            "score_bytes",
                            "work",
                            "evaluations",
                            "visits",
                            "native_reuses",
                            "donors",
                            "probes_started",
                            "probe_conflicts",
                            "probes_completed",
                            "mixed",
                            "mixed_identity",
                            "non_native",
                            "switches",
                            "visited_sources",
                            "mixed_donor_sources",
                            "families",
                            "incumbent",
                            "support",
                            "lost_support",
                            "peak_state_bytes",
                            "peak_tasks",
                        ],
                    )?;
                    match key.as_str() {
                        "native" => {
                            j.expect(b'[')?;
                            let mut first = true;
                            while j.next(&mut first, b']')? {
                                let raw = j.raw()?;
                                mem.reserve(weight(raw.len() as u64)?)?;
                                let n: Scored =
                                    serde_json::from_slice(&raw).map_err(io::Error::other)?;
                                ensure(
                                    native.len() < g.families.len()
                                        && n.assignment == g.native_assignment(native.len())?
                                        && n.objective().is_finite(),
                                    "native baseline/graph mismatch",
                                )?;
                                native.push(n);
                            }
                        }
                        "tasks" => {
                            j.expect(b'{')?;
                            tasks_position = Some(j.position);
                            let mut first = true;
                            let mut last = None;
                            while j.next(&mut first, b'}')? {
                                let key = j.key()?;
                                let id = key.parse::<u64>().map_err(io::Error::other)?;
                                ensure(
                                    key == id.to_string() && last.is_none_or(|n| id > n),
                                    "duplicate/unordered task ID",
                                )?;
                                let t: Task = j.value()?;
                                ensure(t.id == id, "task key/ID mismatch")?;
                                task_valid(&t, g)?;
                                last = Some(id);
                                tasks = plus(tasks, 1)?;
                            }
                        }
                        "visits" => state_visits = Some(j.value::<u64>()?),
                        "evaluations" => state_evals = Some(j.value::<u64>()?),
                        _ => j.skip()?,
                    }
                }
                ensure(fields.len() == 33, "missing v3 state field")?;
            }
            _ => j.skip()?,
        }
    }
    j.end()?;
    ensure(
        version == 3
            && transition.is_null()
            && [
                "version",
                "bindings",
                "parent",
                "budgets",
                "state",
                "transition",
            ]
            .iter()
            .all(|k| seen.contains(*k)),
        "unsupported/incomplete v3 snapshot",
    )?;
    ensure(
        budgets.as_array().is_some_and(|a| !a.is_empty()),
        "missing cumulative snapshot budgets",
    )?;
    let bindings = bindings.ok_or_else(|| invalid("missing bindings"))?;
    ensure(
        bindings.policy == identity::policy_identity()
            && bindings.backend == routes::compiler_identity()
            && bindings.graph == g.digest()?
            && bindings.sample == sample
            && bindings.count_policy == genome::COUNT_POLICY
            && bindings.depth_bits == 10.0f64.to_bits()
            && bindings.background_bits == 0.1f64.to_bits()
            && f64::from_bits(bindings.tie_bits).is_finite()
            && f64::from_bits(bindings.tie_bits) >= 0.0
            && bindings.max_feature_terms > 0,
        "snapshot scientific/policy binding mismatch",
    )?;
    ensure(
        native.len() == g.families.len(),
        "native initialization incomplete",
    )?;
    // Stream every ledger row for EOF and visit continuity. Only native rows survive.
    let mut ledger = stream::Json::open(&o.snapshot.join("evaluations.jsonl"), o.record_bytes)?;
    let ledger_len = fs::metadata(o.snapshot.join("evaluations.jsonl"))?.len();
    let mut natives = 0;
    let mut visits = 0;
    let mut evaluations = 0;
    loop {
        ledger.ws()?;
        if ledger.position == ledger_len {
            break;
        }
        let raw = ledger.raw()?;
        ledger.line_end()?;
        let row: Value = serde_json::from_slice(&raw).map_err(io::Error::other)?;
        visits += 1;
        let kind = row["kind"]
            .as_str()
            .ok_or_else(|| invalid("missing ledger kind"))?;
        ensure(
            [
                "native-initialization",
                "native-score-reuse",
                "fresh-evaluation",
            ]
            .contains(&kind)
                && row["visit"] == visits,
            "ledger visit/kind mismatch",
        )?;
        if kind != "native-score-reuse" {
            evaluations += 1;
        }
        ensure(
            row["evaluation"] == evaluations,
            "ledger evaluation mismatch",
        )?;
        if kind == "native-initialization" {
            let n = native
                .get(natives)
                .ok_or_else(|| invalid("extra native ledger row"))?;
            let assignment: Assignment =
                serde_json::from_value(row["scored"]["assignment"].clone())
                    .map_err(io::Error::other)?;
            let spelling = format!(
                "\"relative_objective\":{}}}",
                serde_json::to_string(&n.objective()).map_err(io::Error::other)?
            );
            ensure(
                assignment == n.assignment
                    && row["task"].is_null()
                    && std::str::from_utf8(&raw)
                        .map_err(io::Error::other)?
                        .contains(&spelling),
                "native objective bits/ledger mismatch",
            )?;
            natives += 1;
        }
    }
    ledger.end()?;
    ensure(
        natives == native.len() && Some(visits) == state_visits && Some(evaluations) == state_evals,
        "ledger/native/state correspondence mismatch",
    )?;
    Ok(Snapshot {
        bindings,
        seal,
        parent,
        budgets,
        transition,
        order_update,
        native,
        tasks_position: tasks_position.ok_or_else(|| invalid("missing tasks"))?,
        tasks,
        preflight_bytes: j.position
            + ledger.position
            + seal_input.position
            + fs::metadata(&path)?.len()
            + ledger_len,
    })
}
