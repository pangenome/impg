//! One charged resumable primitive per service; completion probes never own ordinary work.
use super::*;
pub(super) fn initial(
    g: &routes::Graph,
    family: usize,
    slot: usize,
    completed: Vec<Route>,
) -> Context {
    Context {
        family,
        slot,
        completed,
        segments: vec![],
        source: g.families[family].paths[slot],
        cut: 0,
        reverse: false,
    }
}
fn assignment(checksum: &str, family: usize, routes: Vec<Route>) -> io::Result<Assignment> {
    Ok(Assignment {
        version: routes::VERSION,
        model: routes::MODEL.into(),
        graph_checksum: checksum.into(),
        family,
        routes,
    })
}
fn disjoint(a: &Segment, b: &Segment) -> bool {
    a.source != b.source || a.end <= b.start || b.end <= a.start
}

pub(super) fn advance(
    s: &mut State,
    mut task: Task,
    e: &mut Evaluator<'_>,
    limits: &Limits,
    epsilon: f64,
    ledger: &mut std::fs::File,
) -> io::Result<()> {
    let old = task.clone();
    let emitted_id = s.next_id;
    let g = e.graph;
    let c = &task.context;
    let target = g.families[c.family].paths[c.slot];
    let mut again = false;
    match task.op.clone() {
        Op::Start => {
            s.spawn(
                c.clone(),
                Op::SourceBound {
                    lo: 0,
                    hi: g.lanes[c.source].port_count,
                },
            )?;
            if c.source == target && !c.reverse && c.cut < g.lanes[target].length {
                task.op = Op::Check {
                    piece: c.piece(g.lanes[target].length),
                    index: 0,
                    next: AfterCheck::Close,
                };
                again = true;
            }
        }
        Op::SourceBound { mut lo, mut hi } => {
            if lo < hi {
                let mid = lo + (hi - lo) / 2;
                let mut file = e.ports.source_file(g, c.source)?;
                let p = adapter::read(&mut file, g.k, g.lanes[c.source].port_count, mid)?;
                adapter::verified(&p, g)?;
                ensure(p.source == c.source, "source index lane mismatch")?;
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
                if p.anchor < threshold {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
                task.op = Op::SourceBound { lo, hi };
            } else {
                let (base, len) = if c.reverse {
                    (0, lo)
                } else {
                    (lo, g.lanes[c.source].port_count - lo)
                };
                task.op = Op::SourceScan {
                    base,
                    permutation: SourcePermutation::new(len),
                };
            }
            again = true;
        }
        Op::SourceScan {
            base,
            mut permutation,
        } => {
            if let Some(index) = permutation.next() {
                if let Some(index) = index {
                    let mut file = e.ports.source_file(g, c.source)?;
                    let p = adapter::read(
                        &mut file,
                        g.k,
                        g.lanes[c.source].port_count,
                        plus(base, index)?,
                    )?;
                    adapter::verified(&p, g)?;
                    ensure(p.source == c.source, "source index lane mismatch")?;
                    let cut = p.cut(g.k)?;
                    if p.reverse == c.reverse && if c.reverse { cut < c.cut } else { cut > c.cut } {
                        s.spawn(
                            c.clone(),
                            Op::Check {
                                piece: c.piece(cut),
                                index: 0,
                                next: AfterCheck::Hub(p),
                            },
                        )?;
                    }
                }
                task.op = Op::SourceScan { base, permutation };
                again = true;
            }
        }
        Op::Check { piece, index, next } => {
            ensure(piece.start < piece.end, "nonpositive proposed traversal")?;
            if let Some(prior) = c.spans().nth(index) {
                if disjoint(prior, &piece) {
                    task.op = Op::Check {
                        piece,
                        index: index + 1,
                        next,
                    };
                    again = true;
                }
            } else {
                task.op = match next {
                    AfterCheck::Close => Op::Close { piece },
                    AfterCheck::Hub(port) => Op::HubBound {
                        port,
                        lower: None,
                        lo: 0,
                        hi: g.port_count,
                    },
                };
                again = true;
            }
        }
        Op::Close { piece } => {
            let mut segments = c.segments.clone();
            segments.push(piece);
            let mut completed = c.completed.clone();
            completed.push(Route { segments }.normalize());
            let next_slot = c.slot + 1;
            if next_slot == g.families[c.family].paths.len() {
                task.op = Op::Evaluate {
                    assignment: assignment(
                        &s.native[c.family].assignment.graph_checksum,
                        c.family,
                        completed,
                    )?,
                };
                again = true;
            } else {
                // Ordinary coupled continuation is unconditional after committed-span
                // feasibility. A failed hypothetical native tail destroys ONLY its probe.
                s.probes_started = plus(s.probes_started, 1)?;
                s.spawn(
                    initial(g, c.family, next_slot, completed.clone()),
                    Op::Start,
                )?;
                task.op = Op::Probe {
                    routes: completed,
                    next_slot,
                };
                again = true;
            }
        }
        Op::HubBound {
            port,
            lower,
            mut lo,
            mut hi,
        } => {
            if lo < hi {
                let mid = lo + (hi - lo) / 2;
                let p = adapter::read(&mut e.ports.global, g.k, g.port_count, mid)?;
                adapter::verified(&p, g)?;
                if p.word < port.word || (lower.is_some() && p.word == port.word) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
                task.op = Op::HubBound {
                    port,
                    lower,
                    lo,
                    hi,
                };
            } else if let Some(base) = lower {
                ensure(base < lo, "source port missing hub membership")?;
                task.op = Op::HubScan {
                    port,
                    base,
                    permutation: Permutation::new(lo - base),
                };
            } else {
                task.op = Op::HubBound {
                    port,
                    lower: Some(lo),
                    lo,
                    hi: g.port_count,
                };
            }
            again = true;
        }
        Op::HubScan {
            port,
            base,
            mut permutation,
        } => {
            if let Some(index) = permutation.next() {
                if let Some(index) = index {
                    let donor =
                        adapter::read(&mut e.ports.global, g.k, g.port_count, plus(base, index)?)?;
                    adapter::verified(&donor, g)?;
                    ensure(donor.word == port.word, "hub word mismatch")?;
                    s.donors = plus(s.donors, 1)?;
                    if donor.source != port.source
                        || donor.cut(g.k)? != port.cut(g.k)?
                        || donor.reverse != port.reverse
                    {
                        s.spawn(
                            c.clone(),
                            Op::Child {
                                piece: c.piece(port.cut(g.k)?),
                                donor,
                            },
                        )?;
                    }
                }
                task.op = Op::HubScan {
                    port,
                    base,
                    permutation,
                };
                again = true;
            }
        }
        Op::Child { piece, donor } => {
            let mut child = c.clone();
            child.segments.push(piece);
            child.source = donor.source;
            child.cut = donor.cut(g.k)?;
            child.reverse = donor.reverse;
            s.spawn(child, Op::Start)?;
        }
        Op::Probe {
            mut routes,
            next_slot,
        } => {
            if next_slot < g.families[c.family].paths.len() {
                let source = g.families[c.family].paths[next_slot];
                routes.push(Route {
                    segments: vec![Segment {
                        source,
                        start: 0,
                        end: g.lanes[source].length,
                        reverse: false,
                    }],
                });
                task.op = Op::Probe {
                    routes,
                    next_slot: next_slot + 1,
                };
            } else {
                task.op = Op::ProbeCheck {
                    assignment: assignment(
                        &s.native[c.family].assignment.graph_checksum,
                        c.family,
                        routes,
                    )?,
                    i: 0,
                    j: 1,
                };
            }
            again = true;
        }
        Op::ProbeCheck { assignment, i, j } => {
            let mut spans = assignment.routes.iter().flat_map(|r| &r.segments);
            if let Some(a) = spans.nth(i) {
                if let Some(b) = assignment.routes.iter().flat_map(|r| &r.segments).nth(j) {
                    if disjoint(a, b) {
                        task.op = Op::ProbeCheck {
                            assignment,
                            i,
                            j: j + 1,
                        };
                        again = true;
                    } else {
                        s.probe_conflicts = plus(s.probe_conflicts, 1)?;
                    }
                } else {
                    task.op = Op::ProbeCheck {
                        assignment,
                        i: i + 1,
                        j: i + 2,
                    };
                    again = true;
                }
            } else {
                s.probes_completed = plus(s.probes_completed, 1)?;
                task.op = Op::Evaluate { assignment };
                again = true;
            }
        }
        Op::Evaluate { assignment } => {
            // All completed routes were backend-normalized at closure. Exact equality
            // with a previously validated native assignment needs no new evaluator call.
            let recorded = s
                .native
                .iter()
                .find(|n| n.assignment == assignment)
                .cloned();
            let reuse = recorded.is_some();
            let scored = if let Some(score) = recorded {
                score
            } else {
                let result = e.evaluate(&assignment)?;
                Scored {
                    assignment: result.assignment,
                    objective_bits: result.relative_objective.to_bits(),
                }
            };
            s.record(
                scored,
                reuse,
                false,
                epsilon,
                limits,
                e,
                ledger,
                Some(task.id),
            )?;
        }
    }
    if again {
        s.insert(task)?;
    }
    s.handoff(&old, (s.next_id > emitted_id).then_some(emitted_id))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn canonical_opposite_orientations_conflict_only_on_committed_spans() {
        let a = Segment {
            source: 0,
            start: 1,
            end: 5,
            reverse: false,
        };
        assert!(!disjoint(
            &a,
            &Segment {
                reverse: true,
                ..a.clone()
            }
        ));
        assert!(disjoint(
            &a,
            &Segment {
                start: 5,
                end: 9,
                reverse: true,
                ..a
            }
        ));
    }
    #[test]
    fn fifo_share_and_preference_indices_have_single_ownership() {
        let mut s = State::new(2);
        s.ranking = vec![1, 0];
        for family in 0..2 {
            for depth in 0..5 {
                let c = Context {
                    family,
                    slot: 0,
                    completed: vec![],
                    segments: vec![
                        Segment {
                            source: 0,
                            start: 0,
                            end: 1,
                            reverse: false
                        };
                        depth
                    ],
                    source: 0,
                    cut: 0,
                    reverse: false,
                };
                s.spawn(c, Op::Start).unwrap();
            }
        }
        for _ in 0..320 {
            let (cursor, q, id, mode) = s.select().unwrap();
            let f = s.tasks[&id].context.family;
            s.family_cursor = cursor;
            s.quantum = q;
            s.families[f].work += 1;
            s.families[f].modes[mode] += 1;
            let t = s.remove(id).unwrap();
            s.insert(t).unwrap();
            s.quantum -= 1;
            if s.quantum == 0 {
                s.quantum = 32;
                s.family_cursor = (cursor + 1) % 2;
            }
        }
        for f in &s.families {
            assert_eq!(f.work, 160);
            assert_eq!(f.modes, [54, 53, 53]);
        }
        assert_eq!(s.tasks.len(), 10);
        assert_eq!(s.fifo.len(), 10);
        assert_eq!(s.shallow.len(), 10);
        s.validate_links().unwrap();
    }
}
