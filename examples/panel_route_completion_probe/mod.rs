//! Standalone bounded proposal subset. Nothing here is production search state.
mod identity;
mod machine;
mod schema;
mod snapshot;
mod stream;
// Reuse the current read-only verified-handle adapter, not a scheduler/archive.
#[path = "../../src/commands/genome_infer/panel_route_search_policy/adapter.rs"]
mod adapter;
use adapter::{Permutation, Port, SourcePermutation};
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::{SyncmerParams, SyngIndex},
};
use routes::{Assignment, Evaluator, Route, Segment};
use schema::{Context, Op, Task};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fs::{self, File, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
};
fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s)
}
fn ensure(ok: bool, s: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(s))
    }
}
fn plus(n: u64, x: u64) -> io::Result<u64> {
    n.checked_add(x).ok_or_else(|| invalid("counter overflow"))
}
fn bytes(v: &impl Serialize) -> io::Result<u64> {
    Ok(serde_json::to_vec(v).map_err(io::Error::other)?.len() as u64)
}
#[derive(Clone, Debug, clap::Parser, Serialize)]
pub struct Options {
    #[arg(long)]
    pub panel: String,
    #[arg(long)]
    pub routes: PathBuf,
    #[arg(long)]
    pub sample: PathBuf,
    #[arg(long)]
    pub snapshot: PathBuf,
    #[arg(long)]
    pub out_dir: PathBuf,
    #[arg(long, default_value_t = 5_000_000)]
    pub max_work: u64,
    #[arg(long, default_value_t = 100_000)]
    pub max_evaluations: u64,
    #[arg(long, default_value_t = 134_217_728)]
    pub max_state_bytes: u64,
    #[arg(long, default_value_t = 16)]
    pub producer_slots: usize,
    #[arg(long, default_value_t = 16)]
    pub chain_slots: usize,
    #[arg(long, default_value_t = 1_048_576)]
    pub record_bytes: usize,
    #[arg(long, default_value_t = 200_000)]
    pub max_distinct: usize,
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
}
#[derive(Default, Serialize)]
struct Counts {
    work: u64,
    inspected: u64,
    skipped_native: u64,
    skipped_scored: u64,
    admitted: u64,
    prefix_positions: u64,
    recovered_earlier: u64,
    deferred: u64,
    source_bounds: u64,
    source_ordinals: u64,
    target_bounds: u64,
    target_members: u64,
    proposals: u64,
    capacity_comparisons: u64,
    capacity_rejections: u64,
    evaluation_visits: u64,
    fresh_evaluations: u64,
    native_reuses: u64,
    cached_reuses: u64,
    distinct: u64,
    score_ties: u64,
    closures_before_input_end: u64,
}
/// Logical reservations, not allocator/RSS: wire payload *16 + fixed node slack.
/// Transient record parsing is separately reserved before opening the input.
struct Memory {
    used: u64,
    peak: u64,
    limit: u64,
}
impl Memory {
    fn reserve(&mut self, n: u64) -> io::Result<()> {
        let used = plus(self.used, n)?;
        ensure(used <= self.limit, "state-budget-exhausted")?;
        self.used = used;
        self.peak = self.peak.max(used);
        Ok(())
    }
    fn release(&mut self, n: u64) {
        self.used = self.used.checked_sub(n).expect("reservation ownership");
    }
}
fn weight(wire: u64) -> io::Result<u64> {
    wire.checked_mul(16)
        .and_then(|n| n.checked_add(4096))
        .ok_or_else(|| invalid("payload weight overflow"))
}
fn create_json(path: &Path, value: &impl Serialize) -> io::Result<()> {
    let mut f = OpenOptions::new().write(true).create_new(true).open(path)?;
    serde_json::to_writer_pretty(&mut f, value).map_err(io::Error::other)?;
    f.write_all(b"\n")?;
    f.sync_all()
}
pub fn run(o: Options) -> io::Result<()> {
    ensure(
        o.producer_slots > 0 && o.chain_slots > 0 && o.record_bytes > 0 && o.max_distinct > 0,
        "zero storage/queue cap",
    )?;
    // A new empty output is mandatory. Failure evidence stays here; never reuse it.
    fs::create_dir(&o.out_dir)?;
    create_json(&o.out_dir.join("options.json"), &o)?;
    let result = run_inner(&o);
    if let Err(e) = &result {
        create_json(
            &o.out_dir.join("failure.json"),
            &json!({"error":e.to_string(),"one_shot":true,"not_exhaustion":true}),
        )?;
    }
    result
}
fn run_inner(o: &Options) -> io::Result<()> {
    let mut mem = Memory {
        used: 0,
        peak: 0,
        limit: o.max_state_bytes,
    };
    let transient = weight(o.record_bytes as u64)?
        .checked_add(131072)
        .ok_or_else(|| invalid("transient overflow"))?;
    mem.reserve(transient)?;
    let began = std::time::Instant::now();
    let identity = genome::PanelIdentity::read(&o.panel)?;
    let graph = routes::Graph::load(&o.routes, &identity)?;
    ensure(
        graph.read_lengths == [150],
        "prototype requires exactly L150",
    )?;
    let panel = SyngIndex::load(&o.panel, SyncmerParams::default())?;
    let (sample, checksum) = sample::SampleIndex::load_with_checksum(&o.sample, &identity)?;
    ensure(
        sample.stats.read_lengths.len() == 1 && sample.stats.read_lengths.contains_key(&150),
        "prototype requires L150 sample",
    )?;
    let snap = snapshot::preflight(o, &graph, &checksum, &mut mem)?;
    let mut e = Evaluator::new(
        &o.routes,
        &graph,
        &panel,
        &sample,
        10.0,
        0.1,
        snap.bindings.max_feature_terms,
        snap.bindings.cache_terms,
    )?;
    let preflight_seconds = began.elapsed().as_secs_f64();
    create_json(
        &o.out_dir.join("provenance.json"),
        &json!({"snapshot":fs::canonicalize(&o.snapshot)?,"checkpoint_seal":snap.seal,"bindings":snap.bindings,"parent":snap.parent,"budgets":snap.budgets,"transition":snap.transition,"order_update":snap.order_update,"native":snap.native,"preflight_bytes":snap.preflight_bytes,"preflight_seconds":preflight_seconds,"preflight_is_not_discovery":true,"transient_reserved_bytes":transient,"memory_excludes":"graph/panel/sample/evaluator route cache and atomic evaluator scratch; RSS measured externally"}),
    )?;
    let mut input = stream::Json::open(&o.snapshot.join("checkpoint.json"), o.record_bytes)?;
    input.seek(snap.tasks_position)?;
    let mut engine = machine::Engine::new(o, snap, mem)?;
    engine.drive(&mut input, &mut e)?;
    // Detect mutation during use as well as before admission. No artifact authorizes resumption.
    ensure(
        stream::hash_file(&o.snapshot.join("checkpoint.json"))?
            == engine.snapshot.seal["checkpoint"].as_str().unwrap(),
        "snapshot changed during probe",
    )?;
    ensure(
        stream::hash_file(&o.snapshot.join("evaluations.jsonl"))?
            == engine.snapshot.seal["ledger"].as_str().unwrap(),
        "ledger changed during probe",
    )?;
    engine.finish(input.position)
}
#[cfg(test)]
pub(crate) use machine::mechanism_checks;
