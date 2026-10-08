//! Explicit trusted-local v2 -> v3 source-order update; no search or native evaluation.
use super::*;
use checkpoint::{Bindings, Parent};

pub(super) fn is_v2(policy: &BTreeMap<String, String>) -> bool {
    policy
        .get("version")
        .is_some_and(|v| v == "fair-sample-ranked-live-continuation-v2")
        && [
            "policy",
            "machine",
            "adapter",
            "checkpoint",
            "cli",
            "continuation",
            "transition",
            "v1_schema",
            "validation",
            "validation_containers",
        ]
        .iter()
        .all(|k| {
            policy
                .get(*k)
                .is_some_and(|s| s.starts_with("fnv1a64-source-"))
        })
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Receipt {
    version: u32,
    old_bindings: Bindings,
    new_bindings: Bindings,
    ancestor: Parent,
    work: u64,
    evaluations: u64,
    visits: u64,
    old_task_bytes: u64,
    old_peak_state_bytes: u64,
}
pub(super) fn validate_scans(state: &State, paired: bool) -> io::Result<()> {
    for task in state.tasks.values() {
        if let Op::SourceScan { permutation, .. } = &task.op {
            ensure(
                permutation.paired() == paired,
                "source scan order/schema mismatch",
            )?;
        }
    }
    Ok(())
}
fn transform(state: &mut State) -> io::Result<()> {
    validate_scans(state, false)?;
    for task in state.tasks.values_mut() {
        if let Op::SourceScan { permutation, .. } = &mut task.op {
            permutation.update()?;
        }
    }
    state.task_bytes = state
        .tasks
        .values()
        .try_fold(0, |n, t| plus(n, State::task_weight(t)?))?;
    state.peak_state_bytes = state.peak_state_bytes.max(state.occupancy()?);
    Ok(())
}
fn old_bindings(old: &Bindings, new: &Bindings) -> io::Result<()> {
    ensure(
        is_v2(&old.policy),
        "source order update requires compatible live-continuation-v2 policy",
    )?;
    let mut expected = new.clone();
    expected.policy = old.policy.clone();
    ensure(
        *old == expected,
        "source order update backend/graph/sample/parameter mismatch",
    )
}
pub(super) fn verify_boundary(cp: &checkpoint::Checkpoint, ledger: &str) -> io::Result<()> {
    let receipt = cp.order_update.as_ref().unwrap();
    ensure(
        receipt.version == 1 && receipt.new_bindings == cp.bindings,
        "source order update receipt mismatch",
    )?;
    old_bindings(&receipt.old_bindings, &cp.bindings)?;
    ensure(
        validation::fingerprint(&cp.parent)? == validation::fingerprint(&Some(&receipt.ancestor))?,
        "source order update parent mismatch",
    )?;
    ensure(
        ledger == receipt.ancestor.ledger,
        "update-only ledger changed",
    )?;
    let limits = cp
        .budgets
        .last()
        .ok_or_else(|| invalid("missing update budgets"))?;
    let (mut state, budgets, parent) = checkpoint::load(
        &receipt.ancestor.directory,
        &receipt.old_bindings,
        limits,
        true,
    )?;
    ensure(
        validation::fingerprint(&parent)? == validation::fingerprint(&receipt.ancestor)?
            && budgets == cp.budgets,
        "source order update history mismatch",
    )?;
    ensure(
        (
            state.work,
            state.evaluations,
            state.visits,
            state.task_bytes,
            state.peak_state_bytes,
        ) == (
            receipt.work,
            receipt.evaluations,
            receipt.visits,
            receipt.old_task_bytes,
            receipt.old_peak_state_bytes,
        ),
        "source order update accounting mismatch",
    )?;
    transform(&mut state)?;
    ensure(
        validation::fingerprint(&state)? == validation::fingerprint(&cp.state)?,
        "source order update changed preserved search state",
    )
}
pub(super) fn update(
    e: &Evaluator<'_>,
    from: &Path,
    out: &Path,
    bindings: Bindings,
    limits: Limits,
    extend: bool,
) -> io::Result<serde_json::Value> {
    let (cp, _) = checkpoint::verify(from)?;
    ensure(
        cp.version == 2,
        "source order update requires an unupdated v2 checkpoint",
    )?;
    old_bindings(&cp.bindings, &bindings)?;
    let old = cp.bindings;
    drop(cp.state);
    let (mut state, budgets, parent) = checkpoint::load(from, &old, &limits, extend)?;
    validation::graph(&state, e.graph)?;
    let receipt = Receipt {
        version: 1,
        old_bindings: old,
        new_bindings: bindings.clone(),
        ancestor: parent.clone(),
        work: state.work,
        evaluations: state.evaluations,
        visits: state.visits,
        old_task_bytes: state.task_bytes,
        old_peak_state_bytes: state.peak_state_bytes,
    };
    transform(&mut state)?;
    validation::state(&state)?;
    let parent = Some(parent);
    let metadata = bytes(&(&bindings, &budgets, &parent, &Some(&receipt)))?;
    ensure(
        plus(state.occupancy()?, metadata)? <= limits.max_state_bytes,
        "updated state and metadata exceed authorized state cap",
    )?;
    let result = serde_json::json!({"method":POLICY_VERSION,"status":"updated-source-pair-order-no-search",
        "update_only":true,"work":state.work,"evaluations":state.evaluations,"complete_assignment_visits":state.visits,
        "retained_state_bytes":plus(state.occupancy()?,metadata)?,"receipt":receipt,
        "sequence_emission_authorized":false});
    let mut ledger = checkpoint::ledger(out, parent.as_ref())?;
    checkpoint::save(
        out,
        bindings,
        parent,
        budgets,
        state,
        &mut ledger,
        Some(receipt),
    )?;
    Ok(result)
}
