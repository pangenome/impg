//! Explicit exact-a520 conversion, never a policy-agnostic resume or search restart.
use super::*;
use checkpoint::{Bindings, Parent};
use std::fs::File;

pub(super) fn supported_identity() -> BTreeMap<String, String> {
    let version_line = include_str!("supported_v1/mod.rs")
        .lines()
        .find(|l| l.starts_with("pub const POLICY_VERSION:"))
        .expect("pinned v1 version declaration");
    let version = version_line
        .split('"')
        .nth(1)
        .expect("pinned v1 version literal");
    [
        ("version", version),
        ("policy", include_str!("supported_v1/mod.rs")),
        ("machine", include_str!("supported_v1/machine.rs")),
        ("adapter", include_str!("supported_v1/adapter.rs")),
        ("checkpoint", include_str!("supported_v1/checkpoint.rs")),
        ("cli", include_str!("supported_v1/genome_infer.rs")),
    ]
    .into_iter()
    .map(|(name, s)| {
        (
            name.into(),
            if name == "version" {
                s.into()
            } else {
                format!("fnv1a64-source-{:016x}", genome::checksum(s.as_bytes()))
            },
        )
    })
    .collect()
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct OldCheckpoint {
    version: u32,
    bindings: Bindings,
    parent: Option<Parent>,
    budgets: Vec<Limits>,
    state: v1_schema::OldState,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Receipt {
    version: u32,
    declaration_path: PathBuf,
    declaration_hash: String,
    conversion_only: bool,
    old_schema: u32,
    new_schema: u32,
    old_bindings: Bindings,
    new_bindings: Bindings,
    transform_source: String,
    ancestor: Parent,
    ledger_prefix: String,
    old_task_bytes: u64,
    old_peak_state_bytes: u64,
    old_accounting: String,
    new_accounting: String,
    work_baseline: u64,
    evaluations_baseline: u64,
    visits_baseline: u64,
    historical_modes: Vec<[u64; 3]>,
    budget_history: Vec<Limits>,
    authorized_limits: Limits,
    converted_state: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Declaration {
    version: u32,
    checkpoints: Vec<Parent>,
}
fn declaration(path: &Path, expected_hash: &str) -> io::Result<Declaration> {
    ensure(
        checkpoint::hash_file(path)? == expected_hash,
        "old ancestry declaration changed",
    )?;
    let d: Declaration = serde_json::from_reader(std::io::BufReader::new(File::open(path)?))
        .map_err(io::Error::other)?;
    ensure(
        d.version == 1 && !d.checkpoints.is_empty(),
        "unsupported or empty old ancestry declaration",
    )?;
    let mut seen = BTreeSet::new();
    for p in &d.checkpoints {
        ensure(
            p.directory.is_absolute()
                && std::fs::canonicalize(&p.directory)? == p.directory
                && seen.insert(p.directory.clone()),
            "noncanonical or duplicate old ancestry pin",
        )?;
    }
    ensure(
        checkpoint::hash_file(path)? == expected_hash,
        "old ancestry declaration changed during read",
    )?;
    Ok(d)
}
fn pinned(d: &Declaration, path: &Path, seal: &checkpoint::Seal) -> io::Result<()> {
    let path = std::fs::canonicalize(path)?;
    let pin = d
        .checkpoints
        .iter()
        .find(|p| p.directory == path)
        .ok_or_else(|| invalid("missing independently pinned old ancestor"))?;
    ensure(
        pin.checkpoint == seal.checkpoint && pin.ledger == seal.ledger,
        "independent old ancestry pin mismatch",
    )
}
fn old_verify(
    directory: &Path,
    expected: &Bindings,
) -> io::Result<(OldCheckpoint, checkpoint::Seal)> {
    let seal: checkpoint::Seal = genome::read_json(&directory.join("checkpoint-seal.json"))?;
    ensure(
        seal.version == 1
            && seal.checkpoint == checkpoint::hash_file(&directory.join("checkpoint.json"))?
            && seal.ledger == checkpoint::hash_file(&directory.join("evaluations.jsonl"))?,
        "old checkpoint/ledger seal mismatch",
    )?;
    let cp: OldCheckpoint = serde_json::from_reader(std::io::BufReader::new(File::open(
        directory.join("checkpoint.json"),
    )?))
    .map_err(io::Error::other)?;
    ensure(
        cp.version == 1 && cp.bindings == *expected,
        "unsupported old schema/policy/backend/graph/sample/parameters",
    )?;
    validate_budgets(&cp.budgets)?;
    Ok((cp, seal))
}
fn validate_budgets(history: &[Limits]) -> io::Result<()> {
    ensure(!history.is_empty(), "missing cumulative budget history")?;
    for b in history {
        ensure(
            b.max_work > 0 && b.max_evaluations > 0 && b.max_state_bytes > 0 && b.max_optima > 0,
            "invalid historical budgets",
        )?;
    }
    for w in history.windows(2) {
        authorize(&w[0], &w[1], true)?;
    }
    Ok(())
}
fn authorize(old: &Limits, new: &Limits, extend: bool) -> io::Result<()> {
    ensure(
        new.max_work >= old.max_work
            && new.max_evaluations >= old.max_evaluations
            && new.max_state_bytes >= old.max_state_bytes
            && new.max_optima >= old.max_optima,
        "conversion budgets cannot decrease",
    )?;
    ensure(
        new == old || extend,
        "budget extension requires --extend-budgets",
    )
}
fn load_old(
    directory: &Path,
    expected: &Bindings,
    declaration_path: &Path,
    declaration_hash: &str,
) -> io::Result<(State, Receipt)> {
    ensure(
        expected.policy == supported_identity(),
        "unsupported v1 source identity",
    )?;
    let declaration = declaration(declaration_path, declaration_hash)?;
    let directory = std::fs::canonicalize(directory)?;
    let (cp, seal) = old_verify(&directory, expected)?;
    pinned(&declaration, &directory, &seal)?;
    let old_task_bytes = cp.state.task_bytes;
    let old_peak_state_bytes = cp.state.peak_state_bytes;
    let mut ancestor = cp.parent.clone();
    let mut seen = BTreeSet::from([directory.clone()]);
    while let Some(parent) = ancestor {
        ensure(
            seen.insert(std::fs::canonicalize(&parent.directory)?),
            "old checkpoint ancestry cycle",
        )?;
        let (a, s) = old_verify(&parent.directory, expected)?;
        pinned(&declaration, &parent.directory, &s)?;
        ensure(
            s.checkpoint == parent.checkpoint && s.ledger == parent.ledger,
            "immutable old ancestor mismatch",
        )?;
        let state = a.state.convert()?;
        validation::state(&state)?;
        validation::ledger(&state, &parent.directory.join("evaluations.jsonl"))?;
        ancestor = a.parent;
    }
    ensure(
        seen.len() == declaration.checkpoints.len(),
        "incomplete or extra old ancestry pins",
    )?;
    let state = cp.state.convert()?;
    validation::state(&state)?;
    validation::ledger(&state, &directory.join("evaluations.jsonl"))?;
    let limits = cp.budgets.last().unwrap().clone();
    ensure(
        state.work <= limits.max_work
            && state.evaluations <= limits.max_evaluations
            && state.support.len() <= limits.max_optima,
        "old counters exceed budgets",
    )?;
    let receipt = Receipt {
        version: 1,
        declaration_path: std::fs::canonicalize(declaration_path)?,
        declaration_hash: declaration_hash.into(),
        conversion_only: true,
        old_schema: 1,
        new_schema: 2,
        old_bindings: cp.bindings.clone(),
        new_bindings: cp.bindings,
        transform_source: identity()["transition"].clone(),
        ancestor: Parent {
            directory,
            checkpoint: seal.checkpoint,
            ledger: seal.ledger.clone(),
        },
        ledger_prefix: seal.ledger,
        old_task_bytes,
        old_peak_state_bytes,
        old_accounting: "logical-serialized-payload-plus-index-allowances-v1".into(),
        new_accounting: "logical-serialized-payload-plus-live-links-v2".into(),
        work_baseline: state.work,
        evaluations_baseline: state.evaluations,
        visits_baseline: state.visits,
        historical_modes: state.families.iter().map(|f| f.modes).collect(),
        budget_history: cp.budgets,
        authorized_limits: limits,
        converted_state: validation::fingerprint(&state)?,
    };
    Ok((state, receipt))
}
pub(super) fn verify_boundary(cp: &checkpoint::Checkpoint, ledger: &str) -> io::Result<()> {
    let receipt = cp
        .transition
        .as_ref()
        .ok_or_else(|| invalid("missing transition receipt"))?;
    ensure(
        cp.version == 2
            && source_order::is_v2(&cp.bindings.policy)
            && receipt.new_bindings == cp.bindings,
        "historical transition target policy mismatch",
    )?;
    let mut old = cp.bindings.clone();
    old.policy = supported_identity();
    let (state, mut expected) = load_old(
        &receipt.ancestor.directory,
        &old,
        &receipt.declaration_path,
        &receipt.declaration_hash,
    )?;
    authorize(
        expected.budget_history.last().unwrap(),
        &receipt.authorized_limits,
        true,
    )?;
    expected.new_bindings = cp.bindings.clone();
    // Historical policy source is provenance, not the current implementation.
    expected.transform_source = cp
        .bindings
        .policy
        .get("transition")
        .ok_or_else(|| invalid("missing historical transition source"))?
        .clone();
    expected.authorized_limits = receipt.authorized_limits.clone();
    ensure(
        validation::fingerprint(&expected)? == validation::fingerprint(receipt)?,
        "transition receipt mismatch",
    )?;
    ensure(
        ledger == receipt.ledger_prefix,
        "conversion ledger differs from saved history",
    )?;
    ensure(
        validation::fingerprint(&state)? == validation::fingerprint(&cp.state)?,
        "conversion changed search payload/accounting",
    )?;
    ensure(
        validation::fingerprint(&cp.parent)? == validation::fingerprint(&Some(&receipt.ancestor))?,
        "transition parent mismatch",
    )?;
    let mut budgets = receipt.budget_history.clone();
    if budgets.last() != Some(&receipt.authorized_limits) {
        budgets.push(receipt.authorized_limits.clone());
    }
    ensure(budgets == cp.budgets, "transition budget history mismatch")
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn genuine_supported_v1_source_identity() {
        let actual = supported_identity();
        for (k, v) in [
            ("policy", "4a4ebebca0b9b545"),
            ("machine", "43d1f15c7ab65a37"),
            ("adapter", "82d407db2690469d"),
            ("checkpoint", "a2250f6b87cc6393"),
            ("cli", "b49f2d6f418a441c"),
        ] {
            assert_eq!(actual[k], format!("fnv1a64-source-{v}"));
        }
        assert_eq!(actual["version"], "fair-sample-ranked-lazy-coupled-v1");
    }
}
