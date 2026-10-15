//! Versioned, immutable-parent checkpoints; checksums detect corruption, not malicious forgery.
use super::*;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub(super) struct Bindings {
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
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Parent {
    directory: PathBuf,
    checkpoint: String,
    ledger: String,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Seal {
    version: u32,
    checkpoint: String,
    ledger: String,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Checkpoint {
    version: u32,
    bindings: Bindings,
    parent: Option<Parent>,
    budgets: Vec<Limits>,
    state: State,
}
pub(super) fn hash_file(path: &Path) -> io::Result<String> {
    let mut f = File::open(path)?;
    let mut hash = 0xcbf29ce484222325u64;
    let mut buf = [0; 65536];
    let mut size = 0u64;
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        size = size
            .checked_add(n as u64)
            .ok_or_else(|| invalid("file size overflow"))?;
        for b in &buf[..n] {
            hash = (hash ^ *b as u64).wrapping_mul(0x100000001b3);
        }
    }
    Ok(format!("fnv1a64:{size}:{hash:016x}"))
}
fn verify(directory: &Path) -> io::Result<(Checkpoint, Seal)> {
    let seal: Seal = genome::read_json(&directory.join("checkpoint-seal.json"))?;
    ensure(
        seal.version == 1
            && seal.checkpoint == hash_file(&directory.join("checkpoint.json"))?
            && seal.ledger == hash_file(&directory.join("evaluations.jsonl"))?,
        "checkpoint/ledger seal mismatch",
    )?;
    let checkpoint: Checkpoint = genome::read_json(&directory.join("checkpoint.json"))?;
    ensure(checkpoint.version == 1, "unsupported checkpoint version")?;
    Ok((checkpoint, seal))
}
pub(super) fn load(
    directory: &Path,
    bindings: &Bindings,
    limits: &Limits,
    extend: bool,
) -> io::Result<(State, Vec<Limits>, Parent)> {
    let directory = fs::canonicalize(directory)?;
    let (checkpoint, seal) = verify(&directory)?;
    ensure(
        &checkpoint.bindings == bindings,
        "resume policy/backend/graph/sample/parameter mismatch",
    )?;
    let old = checkpoint
        .budgets
        .last()
        .ok_or_else(|| invalid("missing cumulative budgets"))?;
    ensure(
        limits.max_work >= old.max_work
            && limits.max_evaluations >= old.max_evaluations
            && limits.max_state_bytes >= old.max_state_bytes
            && limits.max_optima >= old.max_optima,
        "resume budgets cannot decrease",
    )?;
    ensure(
        limits == old || extend,
        "budget extension requires --extend-budgets",
    )?;
    let mut ancestor = checkpoint.parent.clone();
    let mut seen = BTreeSet::from([directory.clone()]);
    while let Some(parent) = ancestor {
        ensure(
            seen.insert(fs::canonicalize(&parent.directory)?),
            "checkpoint ancestry cycle",
        )?;
        let (cp, s) = verify(&parent.directory)?;
        ensure(
            s.checkpoint == parent.checkpoint
                && s.ledger == parent.ledger
                && cp.bindings == *bindings,
            "immutable parent checkpoint/ledger mismatch",
        )?;
        ancestor = cp.parent;
    }
    let extended = limits != old;
    let mut budgets = checkpoint.budgets;
    if extended {
        budgets.push(limits.clone());
    }
    Ok((
        checkpoint.state,
        budgets,
        Parent {
            directory,
            checkpoint: seal.checkpoint,
            ledger: seal.ledger,
        },
    ))
}
pub(super) fn ledger(out: &Path, parent: Option<&Parent>) -> io::Result<File> {
    let mut ledger = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(out.join("evaluations.jsonl"))?;
    if let Some(p) = parent {
        io::copy(
            &mut File::open(p.directory.join("evaluations.jsonl"))?,
            &mut ledger,
        )?;
        // Detect changed parent bytes during the copy before any new operation.
        ledger.flush()?;
        ensure(
            hash_file(&out.join("evaluations.jsonl"))? == p.ledger,
            "parent ledger changed during copy",
        )?;
    }
    Ok(ledger)
}
pub(super) fn save(
    out: &Path,
    bindings: Bindings,
    parent: Option<Parent>,
    budgets: Vec<Limits>,
    state: State,
    ledger: &mut File,
) -> io::Result<()> {
    ledger.flush()?;
    ledger.sync_all()?;
    genome::write_json(
        &out.join("checkpoint.json"),
        &Checkpoint {
            version: 1,
            bindings,
            parent,
            budgets,
            state,
        },
    )?;
    genome::write_json(
        &out.join("checkpoint-seal.json"),
        &Seal {
            version: 1,
            checkpoint: hash_file(&out.join("checkpoint.json"))?,
            ledger: hash_file(&out.join("evaluations.jsonl"))?,
        },
    )?;
    #[cfg(unix)]
    File::open(out)?.sync_all()?;
    Ok(())
}
