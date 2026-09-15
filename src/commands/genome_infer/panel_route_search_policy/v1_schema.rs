//! Exact a520 schema, separated from the v2 scheduler; no old search is compiled.
use super::*;
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct OldTask {
    id: u64,
    ready: u64,
    depth: usize,
    context: Context,
    op: Op,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct OldState {
    native: Vec<Scored>,
    ranking: Vec<usize>,
    family_cursor: usize,
    quantum: u8,
    #[serde(deserialize_with = "validation::unique_map")]
    tasks: BTreeMap<u64, OldTask>,
    #[serde(deserialize_with = "validation::unique_set")]
    fifo: BTreeSet<(usize, u64, u64)>,
    #[serde(deserialize_with = "validation::unique_set")]
    shallow: BTreeSet<(usize, usize, u64)>,
    #[serde(deserialize_with = "validation::unique_set")]
    progress: BTreeSet<(usize, u64, u64)>,
    next_id: u64,
    next_ready: u64,
    pub(super) task_bytes: u64,
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
    pub(super) peak_state_bytes: u64,
    peak_tasks: usize,
}

impl OldState {
    pub(super) fn convert(self) -> io::Result<State> {
        ensure(self.progress == self.fifo, "old progress index mismatch")?;
        let old_bytes = self
            .tasks
            .values()
            .try_fold(0u64, |n, t| plus(n, plus(bytes(t)?, 256)?))?;
        ensure(old_bytes == self.task_bytes, "old task occupancy mismatch")?;
        let old_occupancy = plus(
            plus(self.task_bytes, self.score_bytes)?,
            4096 + self.families.len() as u64 * 512
                + (self.switches.len()
                    + self.visited_sources.len()
                    + self.mixed_donor_sources.len()) as u64
                    * 64,
        )?;
        ensure(
            self.peak_state_bytes >= old_occupancy && self.peak_tasks >= self.tasks.len(),
            "invalid historical v1 occupancy high water",
        )?;
        let historical_accounting = Some(HistoricalAccounting {
            newest_ready_modes: self.families.iter().map(|f| f.modes).collect(),
            task_bytes_v1: self.task_bytes,
            peak_state_bytes_v1: self.peak_state_bytes,
        });
        let mut state = State {
            historical_accounting,
            native: self.native,
            ranking: self.ranking,
            family_cursor: self.family_cursor,
            quantum: self.quantum,
            tasks: self
                .tasks
                .into_iter()
                .map(|(id, t)| {
                    (
                        id,
                        Task {
                            id: t.id,
                            ready: t.ready,
                            depth: t.depth,
                            context: t.context,
                            op: t.op,
                            prev: None,
                            next: None,
                        },
                    )
                })
                .collect(),
            fifo: self.fifo,
            shallow: self.shallow,
            focus: vec![None; self.families.len()],
            next_id: self.next_id,
            next_ready: self.next_ready,
            task_bytes: self.task_bytes,
            score_bytes: self.score_bytes,
            work: self.work,
            evaluations: self.evaluations,
            visits: self.visits,
            native_reuses: self.native_reuses,
            donors: self.donors,
            probes_started: self.probes_started,
            probe_conflicts: self.probe_conflicts,
            probes_completed: self.probes_completed,
            mixed: self.mixed,
            mixed_identity: self.mixed_identity,
            non_native: self.non_native,
            switches: self.switches,
            visited_sources: self.visited_sources,
            mixed_donor_sources: self.mixed_donor_sources,
            families: self.families,
            incumbent: self.incumbent,
            support: self.support,
            lost_support: self.lost_support,
            peak_state_bytes: self.peak_state_bytes,
            peak_tasks: self.peak_tasks,
        };
        state.task_bytes = state
            .tasks
            .values()
            .try_fold(0, |n, t| plus(n, State::task_weight(t)?))?;
        state.peak_state_bytes = state.occupancy()?;
        Ok(state)
    }
}
