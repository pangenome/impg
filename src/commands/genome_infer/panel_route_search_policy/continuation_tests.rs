//! Causal scheduling harness and local handoff adversaries, not biological evaluations.
use super::*;
fn context(family: usize) -> Context {
    Context {
        family,
        slot: 0,
        completed: vec![],
        segments: vec![],
        source: 0,
        cut: 0,
        reverse: false,
    }
}
fn port() -> Port {
    Port {
        word: vec![b'A'; 63],
        source: 1,
        anchor: 64,
        reverse: false,
    }
}
fn scanner(s: &mut State, f: usize) -> u64 {
    let id = s.next_id;
    s.spawn(
        context(f),
        Op::SourceScan {
            base: 0,
            permutation: SourcePermutation::new(1 << 20),
        },
    )
    .unwrap();
    id
}
fn service(s: &mut State, id: u64, emit: Option<Op>, retain: bool) {
    let old = s.remove(id).unwrap();
    let emitted = emit.map(|op| {
        let id = s.next_id;
        s.spawn(old.context.clone(), op).unwrap();
        id
    });
    if retain {
        s.insert(old.clone()).unwrap();
    }
    s.handoff(&old, emitted).unwrap();
    s.validate_links().unwrap();
}
#[test]
fn head_middle_tail_peer_and_waiting_producer_handoffs() {
    let mut s = State::new(2);
    let a = scanner(&mut s, 0);
    let peer = scanner(&mut s, 0);
    let other = scanner(&mut s, 1);
    s.seed_focus(a);
    let b = s.next_id;
    service(
        &mut s,
        a,
        Some(Op::HubBound {
            port: port(),
            lower: None,
            lo: 0,
            hi: 8,
        }),
        true,
    );
    assert_eq!(s.focus[0], Some(b));
    assert_eq!(s.tasks[&b].next, Some(a));
    let before = s.focus.clone();
    service(&mut s, peer, Some(Op::Start), true);
    service(&mut s, other, Some(Op::Start), true);
    service(&mut s, a, Some(Op::Start), true);
    assert_eq!(s.focus, before);
    // Nonproducer head -> replacement, waiting tail consumption repairs its neighbor.
    let replacement = s.next_id;
    service(
        &mut s,
        b,
        Some(Op::Child {
            piece: context(0).piece(64),
            donor: port(),
        }),
        false,
    );
    assert_eq!(s.focus[0], Some(replacement));
    service(&mut s, a, None, false);
    assert_eq!(s.tasks[&replacement].next, None);
    service(&mut s, replacement, None, false);
    assert_eq!(s.focus[0], None);
    s.validate_links().unwrap();
}
#[test]
fn waiting_nonproducer_replacement_preserves_position_and_native_completion_unpins() {
    let mut s = State::new(1);
    let a = scanner(&mut s, 0);
    s.seed_focus(a);
    let b = s.next_id;
    service(
        &mut s,
        a,
        Some(Op::Child {
            piece: context(0).piece(64),
            donor: port(),
        }),
        true,
    );
    // Put a third live task before the current head; then replace the middle task.
    let c = s.next_id;
    s.spawn(context(0), Op::SourceBound { lo: 0, hi: 8 })
        .unwrap();
    s.links(c, None, Some(b)).unwrap();
    s.links(b, Some(c), Some(a)).unwrap();
    s.focus[0] = Some(c);
    let replacement = s.next_id;
    service(&mut s, b, Some(Op::Start), false);
    assert_eq!(s.focus[0], Some(c));
    assert_eq!(s.tasks[&c].next, Some(replacement));
    assert_eq!(s.tasks[&a].prev, Some(replacement));
    service(&mut s, replacement, None, false);
    service(&mut s, c, None, false);
    assert_eq!(s.focus[0], Some(a));
    let old = s.remove(a).unwrap();
    let mut native = old.clone();
    native.op = Op::Probe {
        routes: vec![],
        next_slot: 0,
    };
    s.insert(native).unwrap();
    let before = s.remove(a).unwrap();
    s.insert(before.clone()).unwrap();
    s.handoff(&before, None).unwrap();
    assert!(s.focus[0].is_none());
    assert!(s.tasks.contains_key(&a));
    s.validate_links().unwrap();
}
#[test]
fn changed_geometry_includes_same_source_orientation_and_completed_routes() {
    let mut c = context(0);
    assert_eq!(c.depth(), 0);
    c.segments.push(c.piece(10));
    c.cut = 20;
    assert!(c.depth() > 0);
    c.reverse = true;
    assert!(c.depth() > 0);
    c.segments.clear();
    c.completed.push(Route {
        segments: vec![
            Segment {
                source: 0,
                start: 0,
                end: 10,
                reverse: false,
            },
            Segment {
                source: 0,
                start: 20,
                end: 30,
                reverse: false,
            },
        ],
    });
    assert!(c.depth() > 0);
}
#[test]
fn malformed_links_rejected_without_stale_cleanup() {
    let mut s = State::new(2);
    let a = scanner(&mut s, 0);
    let b = scanner(&mut s, 1);
    s.seed_focus(a);
    s.links(a, None, Some(b)).unwrap();
    s.links(b, Some(a), None).unwrap();
    assert!(s.validate_links().is_err());
    s.links(b, None, None).unwrap();
    s.links(a, None, Some(a)).unwrap();
    assert!(s.validate_links().is_err());
}
fn causal(policy: usize) -> [u64; 2] {
    let mut s = State::new(2);
    s.ranking = vec![1, 0];
    for f in 0..2 {
        scanner(&mut s, f);
        scanner(&mut s, f);
    }
    let mut complete = [0; 2];
    let mut stages = [0u64; 8];
    for _ in 0..8192 {
        let (cursor, q, mut id, mode) = s.select().unwrap();
        let f = s.tasks[&id].context.family;
        if mode == 2 && policy != 2 {
            id = if policy == 0 {
                s.fifo
                    .range((f, 0, 0)..=(f, u64::MAX, u64::MAX))
                    .next_back()
                    .unwrap()
                    .2
            } else {
                *s.tasks
                    .iter()
                    .filter(|(_, t)| t.context.family == f)
                    .next_back()
                    .unwrap()
                    .0
            };
        }
        if policy == 2 && mode == 2 {
            s.seed_focus(id);
        }
        s.family_cursor = cursor;
        s.quantum = q;
        s.families[f].work += 1;
        s.families[f].modes[mode] += 1;
        let old = s.remove(id).unwrap();
        let mut t = old.clone();
        let emitted = s.next_id;
        let mut retain = true;
        match t.op.clone() {
            Op::SourceScan {
                base,
                mut permutation,
            } => {
                stages[0] += 1;
                assert!(permutation.next().is_some());
                s.spawn(
                    t.context.clone(),
                    Op::Check {
                        piece: t.context.piece(64),
                        index: 0,
                        next: AfterCheck::Hub(port()),
                    },
                )
                .unwrap();
                t.op = Op::SourceScan { base, permutation };
            }
            Op::Check { piece, index, next } => {
                stages[1] += 1;
                t.op = if index < 8 {
                    Op::Check {
                        piece,
                        index: index + 1,
                        next,
                    }
                } else {
                    Op::HubBound {
                        port: port(),
                        lower: None,
                        lo: 0,
                        hi: 20,
                    }
                };
            }
            Op::HubBound {
                port,
                lower,
                lo,
                hi,
            } => {
                stages[2 + usize::from(lower.is_some())] += 1;
                t.op = if lo < hi {
                    Op::HubBound {
                        port,
                        lower,
                        lo: lo + 1,
                        hi,
                    }
                } else if lower.is_none() {
                    Op::HubBound {
                        port,
                        lower: Some(0),
                        lo: 0,
                        hi: 20,
                    }
                } else {
                    Op::HubScan {
                        port,
                        base: 0,
                        permutation: Permutation::new(1 << 20),
                    }
                };
            }
            Op::HubScan {
                port,
                base,
                mut permutation,
            } => {
                stages[4] += 1;
                assert!(permutation.next().is_some());
                s.spawn(
                    t.context.clone(),
                    Op::Child {
                        piece: t.context.piece(64),
                        donor: port.clone(),
                    },
                )
                .unwrap();
                t.op = Op::HubScan {
                    port,
                    base,
                    permutation,
                };
            }
            Op::Child { piece, .. } => {
                stages[5] += 1;
                let mut c = t.context.clone();
                c.segments.push(piece);
                s.spawn(
                    c,
                    Op::Close {
                        piece: t.context.piece(128),
                    },
                )
                .unwrap();
                retain = false;
            }
            Op::Close { .. } => {
                stages[6] += 1;
                t.op = Op::Evaluate {
                    assignment: Assignment {
                        version: 1,
                        model: routes::MODEL.into(),
                        graph_checksum: "harness-not-graph".into(),
                        family: f,
                        routes: vec![],
                    },
                };
            }
            Op::Evaluate { .. } => {
                stages[7] += 1;
                complete[f] += 1;
                retain = false;
            }
            _ => panic!("unexpected harness stage"),
        }
        if retain {
            s.insert(t).unwrap();
        }
        if policy == 2 {
            s.handoff(&old, (s.next_id > emitted).then_some(emitted))
                .unwrap();
        }
        s.quantum -= 1;
        if s.quantum == 0 {
            s.quantum = 32;
            s.family_cursor = (cursor + 1) % 2;
        }
        s.validate_links().unwrap();
    }
    for f in &s.families {
        assert_eq!(f.work, 4096);
        assert_eq!(f.modes, [1366, 1365, 1365]);
    }
    eprintln!("CAUSAL policy={policy} source/check/lower/upper/member/child/close/harness-complete={stages:?} completions={complete:?}; these are NOT biological assignments");
    complete
}
#[test]
fn long_producer_causal_red_red_green() {
    assert_eq!(causal(0), [0, 0], "newest-ready regression must remain RED");
    assert_eq!(
        causal(1),
        [0, 0],
        "greatest-creation regression must remain RED"
    );
    assert!(
        causal(2).iter().all(|&n| n > 0),
        "live continuation must be GREEN"
    );
}

#[test]
fn start_closure_probe_success_and_conflict_keep_ordinary_returns() {
    for changed in [false, true] {
        let mut s = State::new(1);
        let mut c = context(0);
        if changed {
            c.segments.push(c.piece(10));
            c.cut = 20;
        }
        s.spawn(c.clone(), Op::Start).unwrap();
        s.seed_focus(0);
        service(&mut s, 0, Some(Op::SourceBound { lo: 0, hi: 32 }), true);
        assert_eq!(s.focus[0], Some(if changed { 0 } else { 1 }));
        if changed {
            assert_eq!(s.tasks[&0].next, Some(1));
        }
    }
    for conflict in [false, true] {
        let mut s = State::new(1);
        let mut c = context(0);
        c.segments.push(c.piece(10));
        s.spawn(c.clone(), Op::Close { piece: c.piece(20) })
            .unwrap();
        s.seed_focus(0);
        service(&mut s, 0, Some(Op::Start), true);
        assert_eq!(s.tasks[&0].next, Some(1));
        let mut probe = s.remove(0).unwrap();
        probe.op = Op::ProbeCheck {
            assignment: Assignment {
                version: 1,
                model: routes::MODEL.into(),
                graph_checksum: "harness".into(),
                family: 0,
                routes: vec![],
            },
            i: 0,
            j: 1,
        };
        s.insert(probe).unwrap();
        if !conflict {
            service(&mut s, 0, None, true);
            assert_eq!(s.focus[0], Some(0));
            let mut eval = s.remove(0).unwrap();
            eval.op = Op::Evaluate {
                assignment: Assignment {
                    version: 1,
                    model: routes::MODEL.into(),
                    graph_checksum: "harness".into(),
                    family: 0,
                    routes: vec![],
                },
            };
            s.insert(eval).unwrap();
        }
        service(&mut s, 0, None, false);
        assert_eq!(s.focus[0], Some(1));
        assert!(s.tasks.contains_key(&1));
    }
}
#[test]
fn long_dead_chain_and_exhaustion_have_no_history_or_depth_cap() {
    let mut s = State::new(1);
    let root = scanner(&mut s, 0);
    s.seed_focus(root);
    for _ in 0..2048 {
        let id = s.focus[0].unwrap();
        service(
            &mut s,
            id,
            Some(Op::SourceScan {
                base: 0,
                permutation: SourcePermutation::new(1 << 20),
            }),
            true,
        );
    }
    assert_eq!(s.tasks.len(), 2049);
    while let Some(id) = s.focus[0] {
        service(&mut s, id, None, false);
    }
    assert!(s.tasks.is_empty() && s.fifo.is_empty() && s.shallow.is_empty());
    assert_eq!(s.task_bytes, 0);
}
#[test]
#[ignore = "intentional RED proof: current newest-ready cannot complete the frozen causal harness"]
fn newest_ready_intentional_red() {
    assert!(causal(0).iter().all(|&n| n > 0));
}
#[test]
#[ignore = "intentional RED proof: naive greatest-creation-ID cannot complete the frozen causal harness"]
fn greatest_creation_intentional_red() {
    assert!(causal(1).iter().all(|&n| n > 0));
}
