//! Non-scoring checkpoint consistency checks, stronger than the historical v1 loader.
use super::*;
use std::io::BufRead;

pub(super) fn state(s: &State) -> io::Result<()> {
    s.validate_links()?;
    ensure(s.quantum > 0 && s.quantum <= 32, "invalid family quantum")?;
    ensure(
        s.native.len() <= s.families.len(),
        "invalid initialization cursor",
    )?;
    if s.native.len() == s.families.len() {
        ensure(
            s.ranking.iter().copied().collect::<BTreeSet<_>>() == (0..s.families.len()).collect()
                && s.ranking.len() == s.families.len()
                && s.family_cursor < s.ranking.len(),
            "invalid family ranking/cursor",
        )?;
    } else {
        ensure(s.ranking.is_empty(), "ranking before initialization")?;
    }
    ensure(
        s.fifo.len() == s.tasks.len() && s.shallow.len() == s.tasks.len(),
        "task/index cardinality mismatch",
    )?;
    let mut task_bytes = 0;
    for (&id, t) in &s.tasks {
        ensure(
            id == t.id
                && id < s.next_id
                && t.ready < s.next_ready
                && t.context.family < s.native.len()
                && t.depth == t.context.depth(),
            "invalid task identity/context/clock",
        )?;
        ensure(
            s.fifo.contains(&(t.context.family, t.ready, id))
                && s.shallow.contains(&(t.context.family, t.depth, id)),
            "task/index mismatch",
        )?;
        task_bytes = plus(task_bytes, State::task_weight(t)?)?;
    }
    ensure(task_bytes == s.task_bytes, "task occupancy mismatch")?;
    let score_bytes = s
        .native
        .iter()
        .chain(s.incumbent.iter())
        .chain(&s.support)
        .try_fold(0, |n, x| {
            ensure(x.objective().is_finite(), "nonfinite stored objective")?;
            plus(n, State::score_weight(x)?)
        })?;
    ensure(score_bytes == s.score_bytes, "score occupancy mismatch")?;
    if let Some(h) = &s.historical_accounting {
        ensure(
            h.newest_ready_modes.len() == s.families.len(),
            "historical family count mismatch",
        )?;
        for (old, f) in h.newest_ready_modes.iter().zip(&s.families) {
            ensure(
                old.iter().zip(f.modes).all(|(&a, b)| a <= b),
                "historical modes exceed cumulative service",
            )?;
        }
    }
    let mut work = s.native.len() as u64;
    let mut visits = 0;
    let mut evaluations = 0;
    let mut mixed = 0;
    for f in &s.families {
        ensure(
            f.modes
                == [
                    f.work / 3 + u64::from(f.work % 3 > 0),
                    f.work / 3 + u64::from(f.work % 3 > 1),
                    f.work / 3,
                ],
            "service mode counter mismatch",
        )?;
        work = plus(work, f.work)?;
        visits = plus(visits, f.visits)?;
        evaluations = plus(evaluations, f.fresh)?;
        mixed = plus(mixed, f.mixed)?;
    }
    ensure(
        work == s.work
            && visits == s.visits
            && evaluations == s.evaluations
            && mixed == s.mixed
            && plus(s.evaluations, s.native_reuses)? == s.visits,
        "cumulative counter mismatch",
    )?;
    ensure(
        s.peak_tasks >= s.tasks.len() && s.peak_state_bytes >= s.occupancy()?,
        "invalid occupancy high water",
    )
}

pub(super) fn ledger(s: &State, path: &Path) -> io::Result<()> {
    let mut visits = 0u64;
    let mut evaluations = 0u64;
    let mut natives = 0usize;
    let mut work = 0;
    let mut family_visits = vec![0u64; s.families.len()];
    let mut family_fresh = family_visits.clone();
    for line in std::io::BufReader::new(std::fs::File::open(path)?).lines() {
        let line = line?;
        let row: serde_json::Value = serde_json::from_str(&line).map_err(io::Error::other)?;
        visits = plus(visits, 1)?;
        let kind = row["kind"]
            .as_str()
            .ok_or_else(|| invalid("missing ledger kind"))?;
        ensure(
            matches!(
                kind,
                "native-initialization" | "native-score-reuse" | "fresh-evaluation"
            ),
            "unknown ledger kind",
        )?;
        if kind != "native-score-reuse" {
            evaluations = plus(evaluations, 1)?;
        }
        let a: Assignment = serde_json::from_value(row["scored"]["assignment"].clone())
            .map_err(io::Error::other)?;
        ensure(a.family < s.families.len(), "ledger family out of range")?;
        family_visits[a.family] += 1;
        if kind != "native-score-reuse" {
            family_fresh[a.family] += 1;
        }
        ensure(
            row["visit"] == visits && row["evaluation"] == evaluations,
            "ledger counters mismatch",
        )?;
        let w = row["work"]
            .as_u64()
            .ok_or_else(|| invalid("missing ledger work"))?;
        ensure(w > work && w <= s.work, "ledger work mismatch")?;
        work = w;
        ensure(
            row["assignment_checksum"]
                == format!(
                    "{:016x}",
                    genome::checksum(&serde_json::to_vec(&a).map_err(io::Error::other)?)
                ),
            "ledger assignment checksum mismatch",
        )?;
        if kind == "native-initialization" {
            let n = s
                .native
                .get(natives)
                .ok_or_else(|| invalid("extra native ledger row"))?;
            ensure(
                a == n.assignment && a.family == natives && row["task"].is_null(),
                "native ledger assignment mismatch",
            )?;
            // Compare the original numeric token, avoiding a second floating parse.
            ensure(
                line.contains(&format!(
                    "\"relative_objective\":{}}}",
                    serde_json::to_string(&n.objective()).map_err(io::Error::other)?
                )),
                "native objective bits/ledger mismatch",
            )?;
            natives += 1;
        } else {
            ensure(
                row["task"].as_u64().is_some_and(|id| id < s.next_id),
                "ledger task clock mismatch",
            )?;
        }
    }
    ensure(
        visits == s.visits && evaluations == s.evaluations && natives == s.native.len(),
        "ledger/state totals mismatch",
    )?;
    for (i, f) in s.families.iter().enumerate() {
        ensure(
            f.visits == family_visits[i] && f.fresh == family_fresh[i],
            "ledger family totals mismatch",
        )?;
    }
    Ok(())
}

// Streaming semantic fingerprint: no second all-state JSON string/Value allocation.
pub(super) fn fingerprint(value: &impl Serialize) -> io::Result<String> {
    struct Sink(u64, u64);
    impl Write for Sink {
        fn write(&mut self, b: &[u8]) -> io::Result<usize> {
            for x in b {
                self.0 = (self.0 ^ *x as u64).wrapping_mul(0x100000001b3);
            }
            self.1 += b.len() as u64;
            Ok(b.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let mut sink = Sink(0xcbf29ce484222325, 0);
    serde_json::to_writer(&mut sink, value).map_err(io::Error::other)?;
    Ok(format!("fnv1a64:{}:{:016x}", sink.1, sink.0))
}

pub(super) use super::validation_containers::{unique_map, unique_set};
pub(super) fn counter_room(s: &State) -> io::Result<()> {
    for n in [
        s.work,
        s.next_id,
        s.next_ready,
        s.evaluations,
        s.visits,
        s.native_reuses,
        s.donors,
        s.probes_started,
        s.probe_conflicts,
        s.probes_completed,
        s.mixed,
        s.mixed_identity,
        s.non_native,
    ] {
        plus(n, 3)?;
    }
    for f in &s.families {
        for n in [f.work, f.visits, f.fresh, f.mixed]
            .into_iter()
            .chain(f.modes)
        {
            plus(n, 3)?;
        }
    }
    for &n in s.switches.values() {
        plus(n, 1)?;
    }
    Ok(())
}
pub(super) fn graph(s: &State, g: &routes::Graph) -> io::Result<()> {
    ensure(
        s.families.len() == g.families.len(),
        "checkpoint graph family count mismatch",
    )?;
    for (i, n) in s.native.iter().enumerate() {
        ensure(
            n.assignment == g.native_assignment(i)?,
            "native assignment differs from graph",
        )?;
    }
    if !s.ranking.is_empty() {
        let mut rank = (0..g.families.len()).collect::<Vec<_>>();
        rank.sort_by(|&a, &b| {
            s.native[a]
                .objective()
                .total_cmp(&s.native[b].objective())
                .then_with(|| g.families[a].identity.cmp(&g.families[b].identity))
        });
        ensure(s.ranking == rank, "native-ranked order mismatch")?;
    }
    let segment = |p: &Segment| {
        ensure(
            p.source < g.lanes.len() && p.start < p.end && p.end <= g.lanes[p.source].length,
            "checkpoint segment outside graph",
        )
    };
    let assignment = |a: &Assignment| -> io::Result<()> {
        ensure(
            a.version == routes::VERSION
                && a.model == routes::MODEL
                && a.graph_checksum == g.digest()?
                && a.family < g.families.len()
                && a.routes.len() == g.families[a.family].paths.len(),
            "pending assignment graph mismatch",
        )?;
        for r in &a.routes {
            ensure(!r.segments.is_empty(), "empty pending route")?;
            for p in &r.segments {
                segment(p)?;
            }
        }
        Ok(())
    };
    for t in s.tasks.values() {
        let c = &t.context;
        ensure(
            c.family < g.families.len()
                && c.slot < g.families[c.family].paths.len()
                && c.source < g.lanes.len()
                && c.cut <= g.lanes[c.source].length
                && c.completed.len() == c.slot,
            "pending context outside graph",
        )?;
        for p in c.spans() {
            segment(p)?;
        }
        match &t.op {
            Op::Start => {}
            Op::SourceBound { lo, hi } => ensure(
                lo <= hi && *hi <= g.lanes[c.source].port_count,
                "invalid source bounds",
            )?,
            Op::SourceScan { base, permutation } => {
                permutation.validate(*base, g.lanes[c.source].port_count)?
            }
            Op::Check { piece, index, next } => {
                segment(piece)?;
                ensure(*index <= c.spans().count(), "invalid feasibility cursor")?;
                if let AfterCheck::Hub(p) = next {
                    adapter::verified(p, g)?;
                }
            }
            Op::Close { piece } => segment(piece)?,
            Op::HubBound {
                port,
                lower,
                lo,
                hi,
            } => {
                adapter::verified(port, g)?;
                ensure(
                    lo <= hi && *hi <= g.port_count && lower.is_none_or(|x| x <= *lo),
                    "invalid hub bounds",
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
            Op::Child { piece, donor } => {
                segment(piece)?;
                adapter::verified(donor, g)?;
            }
            Op::Probe { routes, next_slot } => {
                ensure(
                    routes.len() == *next_slot && *next_slot <= g.families[c.family].paths.len(),
                    "invalid native-tail cursor",
                )?;
                for r in routes {
                    for p in &r.segments {
                        segment(p)?;
                    }
                }
            }
            Op::ProbeCheck {
                assignment: a,
                i,
                j,
            } => {
                assignment(a)?;
                let n = a.routes.iter().map(|r| r.segments.len()).sum::<usize>();
                ensure(
                    *i <= n && *j <= n + 1 && *j > *i,
                    "invalid probe comparison cursor",
                )?;
            }
            Op::Evaluate { assignment: a } => assignment(a)?,
        }
    }
    for x in s.incumbent.iter().chain(&s.support) {
        assignment(&x.assignment)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strict_duplicate_containers_and_resource_overflow() {
        #[derive(Deserialize)]
        struct Map {
            #[serde(deserialize_with = "unique_map")]
            tasks: BTreeMap<u64, u64>,
        }
        #[derive(Deserialize)]
        struct Set {
            #[serde(deserialize_with = "unique_set")]
            fifo: BTreeSet<u64>,
        }
        assert!(serde_json::from_str::<Map>(r#"{"tasks":{"1":0,"1":1}}"#).is_err());
        assert!(serde_json::from_str::<Set>(r#"{"fifo":[1,1]}"#).is_err());
        assert_eq!(
            serde_json::from_str::<Map>(r#"{"tasks":{"1":0}}"#)
                .unwrap()
                .tasks
                .len(),
            1
        );
        assert_eq!(
            serde_json::from_str::<Set>(r#"{"fifo":[1]}"#)
                .unwrap()
                .fifo
                .len(),
            1
        );
        let mut s = State::new(1);
        s.next_ready = u64::MAX;
        let before = fingerprint(&s).unwrap();
        assert!(counter_room(&s).is_err());
        assert_eq!(fingerprint(&s).unwrap(), before);
        let a = Assignment {
            version: 1,
            model: routes::MODEL.into(),
            graph_checksum: "test".into(),
            family: 0,
            routes: vec![],
        };
        assert!(reservation(u64::MAX, &a, 63).is_err());
        assert!(reservation(0, &a, u64::MAX).is_err());
    }
    #[test]
    fn internal_counter_occupancy_and_link_checks_are_independent_of_seals() {
        let mut s = State::new(1);
        s.peak_state_bytes = s.occupancy().unwrap();
        state(&s).unwrap();
        s.work = 1;
        assert!(state(&s).is_err());
        s.work = 0;
        s.task_bytes = 1;
        assert!(state(&s).is_err());
        s.task_bytes = 0;
        s.focus[0] = Some(9);
        assert!(state(&s).is_err());
    }
}
