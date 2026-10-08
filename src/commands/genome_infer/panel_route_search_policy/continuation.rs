//! Live task-owned continuation/return links; no consumed-task history or extra index.
use super::*;

impl State {
    fn linked(&self, task: &Task) -> bool {
        self.focus[task.context.family] == Some(task.id) || task.prev.is_some()
    }
    fn links(&mut self, id: u64, prev: Option<u64>, next: Option<u64>) -> io::Result<()> {
        let t = self
            .tasks
            .get_mut(&id)
            .ok_or_else(|| invalid("missing continuation neighbor"))?;
        let old = Self::task_weight(t)?;
        t.prev = prev;
        t.next = next;
        self.task_bytes = plus(
            self.task_bytes
                .checked_sub(old)
                .ok_or_else(|| invalid("link occupancy underflow"))?,
            Self::task_weight(t)?,
        )?;
        Ok(())
    }
    pub(super) fn seed_focus(&mut self, id: u64) {
        let f = self.tasks[&id].context.family;
        if self.focus[f].is_none() {
            self.focus[f] = Some(id);
        }
    }
    // Called after exactly one primitive. The old task has been removed and any
    // retained task/emission is back in ordinary indices before the local splice.
    pub(super) fn handoff(&mut self, old: &Task, emitted: Option<u64>) -> io::Result<()> {
        if !self.linked(old) {
            return Ok(());
        }
        let active = self.focus[old.context.family] == Some(old.id);
        let retained = self.tasks.contains_key(&old.id);
        let producer = matches!(old.op, Op::SourceScan { .. } | Op::HubScan { .. });
        let changed = old.context.depth() > 0;
        let mut chain = Vec::with_capacity(2);
        if producer {
            if active {
                chain.extend(emitted);
            }
            if retained {
                chain.push(old.id);
            }
        } else {
            match old.op {
                Op::Start => {
                    if changed && retained {
                        chain.push(old.id);
                    }
                    chain.extend(emitted);
                }
                Op::Close { .. } => {
                    if changed && retained {
                        chain.push(old.id);
                    }
                    chain.extend(emitted);
                }
                Op::Check {
                    next: AfterCheck::Close,
                    ..
                }
                | Op::Probe { .. }
                | Op::ProbeCheck { .. }
                | Op::Evaluate { .. }
                    if !changed => {}
                _ => {
                    if retained {
                        chain.push(old.id);
                    }
                    chain.extend(emitted);
                }
            }
        }
        if retained {
            self.links(old.id, None, None)?;
        }
        let first = chain.first().copied().or(old.next);
        let last = chain.last().copied().or(old.prev);
        if let Some(prev) = old.prev {
            let before = self.tasks[&prev].prev;
            self.links(prev, before, first)?;
        } else {
            self.focus[old.context.family] = first;
        }
        if let Some(next) = old.next {
            let after = self.tasks[&next].next;
            self.links(next, last, after)?;
        }
        for (i, &id) in chain.iter().enumerate() {
            self.links(
                id,
                if i == 0 { old.prev } else { Some(chain[i - 1]) },
                chain.get(i + 1).copied().or(old.next),
            )?;
        }
        Ok(())
    }
    pub(super) fn validate_links(&self) -> io::Result<()> {
        ensure(
            self.focus.len() == self.families.len(),
            "focus family count mismatch",
        )?;
        let mut members = BTreeSet::new();
        for (family, &head) in self.focus.iter().enumerate() {
            let mut prev = None;
            let mut cursor = head;
            while let Some(id) = cursor {
                ensure(
                    members.insert(id),
                    "continuation cycle or repeated membership",
                )?;
                let t = self
                    .tasks
                    .get(&id)
                    .ok_or_else(|| invalid("dangling continuation"))?;
                ensure(
                    t.context.family == family && t.prev == prev,
                    "nonreciprocal or cross-family continuation",
                )?;
                prev = cursor;
                cursor = t.next;
            }
        }
        for (&id, t) in &self.tasks {
            ensure(
                members.contains(&id) || (t.prev.is_none() && t.next.is_none()),
                "orphan continuation links",
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "continuation_tests.rs"]
mod tests;
