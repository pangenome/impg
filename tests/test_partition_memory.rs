// SIMULATED unit-test intervals only; not biological results.
use coitrees::Interval;
use impg::commands::partition::rehome_singleton_slivers;
use log::info;
fn original_rehome(collected_partitions: &mut Vec<(usize, Vec<Interval<u32>>)>) {
    if collected_partitions.is_empty() {
        return;
    }
    type Row = (u32, i32, i32, usize, Interval<u32>);
    let mut rows: Vec<Row> = Vec::new();
    for (pidx, (_, ivs)) in collected_partitions.iter().enumerate() {
        for iv in ivs {
            let (s, e) = if iv.first <= iv.last {
                (iv.first, iv.last)
            } else {
                (iv.last, iv.first)
            };
            rows.push((iv.metadata, s, e, pidx, *iv));
        }
    }
    rows.sort_by_key(|r| (r.0, r.1, r.2));

    let mut counts = vec![0usize; collected_partitions.len()];
    for r in &rows {
        counts[r.3] += 1;
    }
    let initial_singletons = counts.iter().filter(|&&c| c == 1).count();
    if initial_singletons == 0 {
        return;
    }

    let mut pass = 0;
    loop {
        pass += 1;
        let singletons: std::collections::HashSet<usize> = counts
            .iter()
            .enumerate()
            .filter(|(_, &c)| c == 1)
            .map(|(i, _)| i)
            .collect();
        let mut pending: Vec<(usize, usize)> = Vec::new();
        for i in 0..rows.len() {
            let (c, s, e, pidx, _) = rows[i];
            if !singletons.contains(&pidx) {
                continue;
            }
            let left = if i > 0 && rows[i - 1].0 == c && rows[i - 1].2 == s {
                Some(rows[i - 1].3)
            } else {
                None
            };
            let right = if i + 1 < rows.len() && rows[i + 1].0 == c && rows[i + 1].1 == e {
                Some(rows[i + 1].3)
            } else {
                None
            };
            let ls = left.map(|p| !singletons.contains(&p)).unwrap_or(false);
            let rs = right.map(|p| !singletons.contains(&p)).unwrap_or(false);
            let target = if ls && rs {
                let lp = left.unwrap();
                let rp = right.unwrap();
                if counts[lp] >= counts[rp] {
                    lp
                } else {
                    rp
                }
            } else if ls {
                left.unwrap()
            } else if rs {
                right.unwrap()
            } else {
                continue;
            };
            if target != pidx {
                pending.push((i, target));
            }
        }
        if pending.is_empty() || pass > 100 {
            break;
        }
        for (row_idx, new_pidx) in &pending {
            let old_pidx = rows[*row_idx].3;
            counts[old_pidx] -= 1;
            counts[*new_pidx] += 1;
            rows[*row_idx].3 = *new_pidx;
        }
    }
    let remaining = counts.iter().filter(|&&c| c == 1).count();
    info!(
        "Rehomed {} singleton slivers into flank partitions ({} → {} singletons, {} passes)",
        initial_singletons - remaining,
        initial_singletons,
        remaining,
        pass
    );

    let mut new_intervals: Vec<Vec<Interval<u32>>> = (0..collected_partitions.len())
        .map(|_| Vec::new())
        .collect();
    for (_, _, _, pidx, iv) in rows.drain(..) {
        new_intervals[pidx].push(iv);
    }
    let rebuilt: Vec<(usize, Vec<Interval<u32>>)> = collected_partitions
        .iter()
        .enumerate()
        .filter_map(|(i, (pnum, _))| {
            let ivs = std::mem::take(&mut new_intervals[i]);
            if ivs.is_empty() {
                None
            } else {
                Some((*pnum, ivs))
            }
        })
        .collect();
    *collected_partitions = rebuilt;
}

fn signature(p: &[(usize, Vec<Interval<u32>>)]) -> Vec<(usize, Vec<(u32, i32, i32)>)> {
    p.iter()
        .map(|(n, v)| {
            (
                *n,
                v.iter().map(|x| (x.metadata, x.first, x.last)).collect(),
            )
        })
        .collect()
}
#[test]
fn simulated_rehome_matches_original() {
    // Exercise both flanks, ties, reverse orientation, empty partitions,
    // equal coordinate sort keys, cascades past the 100-application cutoff.
    for seed in 0..80usize {
        let mut p = vec![(0, vec![]); 150];
        for i in 0..150 {
            p[i].0 = i * 3;
        }
        for i in 0..450usize {
            let part = if seed == 0 {
                (i / 3).min(149)
            } else {
                (i * (seed + 1) + i / 7) % 150
            };
            let first = (i / 3) as i32 * 10;
            let (a, b) = if (i + seed) % 2 == 0 {
                (first, first + 10)
            } else {
                (first + 10, first)
            };
            p[part].1.push(Interval {
                first: a,
                last: b,
                metadata: (i % 3) as u32,
            });
        }
        let mut expected = p.clone();
        original_rehome(&mut expected);
        rehome_singleton_slivers(&mut p);
        assert_eq!(signature(&p), signature(&expected), "seed {seed}");
    }
    // SIMULATED mixed occupancy, ties, reverse intervals and multiple sequences.
    for seed in 1..41u64 {
        let mut state = seed;
        let mut p: Vec<_> = (0..400).map(|i| (i * 7, Vec::new())).collect();
        for i in 0..800 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let part = (state >> 32) as usize % p.len();
            let start = (i / 4) * 10;
            let (first, last) = if state & 1 == 0 {
                (start, start + 10)
            } else {
                (start + 10, start)
            };
            p[part].1.push(Interval {
                first,
                last,
                metadata: (i % 4) as u32,
            });
        }
        let mut expected = p.clone();
        original_rehome(&mut expected);
        rehome_singleton_slivers(&mut p);
        assert_eq!(signature(&p), signature(&expected), "mixed seed {seed}");
    }
    let mut tie = vec![
        (
            10,
            vec![
                Interval {
                    first: 0,
                    last: 1,
                    metadata: 0,
                },
                Interval {
                    first: 0,
                    last: 1,
                    metadata: 1,
                },
            ],
        ),
        (
            20,
            vec![Interval {
                first: 2,
                last: 1,
                metadata: 0,
            }],
        ),
        (
            30,
            vec![
                Interval {
                    first: 2,
                    last: 3,
                    metadata: 0,
                },
                Interval {
                    first: 0,
                    last: 1,
                    metadata: 2,
                },
            ],
        ),
    ];
    let mut expected = tie.clone();
    original_rehome(&mut expected);
    rehome_singleton_slivers(&mut tie);
    assert_eq!(signature(&tie), signature(&expected));
    assert_eq!(tie[0].1.len(), 3); // equal flanks choose left
    let mut p: Vec<_> = (0..250)
        .map(|i| {
            (
                i,
                vec![Interval {
                    first: i as i32,
                    last: i as i32 + 1,
                    metadata: 0,
                }],
            )
        })
        .collect();
    p[0].1.push(Interval {
        first: 0,
        last: 1,
        metadata: 1,
    });
    let mut expected = p.clone();
    original_rehome(&mut expected);
    rehome_singleton_slivers(&mut p);
    assert_eq!(signature(&p), signature(&expected));
}

#[test]
#[ignore = "requires IMPG_MEMORY_REAL_LIST pointing at indexed real PAFs"]
fn real_multi_cache_cold_reload() {
    use impg::impg_index::ImpgIndex;
    use impg::multi_impg::MultiImpg;
    use std::path::PathBuf;
    let list = std::env::var("IMPG_MEMORY_REAL_LIST").expect("real indexed PAF list required");
    let files: Vec<String> = std::fs::read_to_string(list)
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect();
    let paths: Vec<PathBuf> = files
        .iter()
        .map(|f| PathBuf::from(format!("{f}.impg")))
        .collect();
    let index = MultiImpg::load_from_files(&paths, &files, None).unwrap();
    assert_eq!(index.cache_stats(), (0, 0, 0));
    let mut targets = index.target_ids();
    targets.sort();
    let target = targets[0];
    let query = || {
        index
            .query(target, 0, 5000, true, None, None, false)
            .unwrap()
    };
    let first = format!("{:?}", query());
    let loaded = index.cache_stats();
    assert!(
        loaded.0 > 1 && loaded.1 > 1 && loaded.2 > 0,
        "must exercise all sub-indices: {loaded:?}"
    );
    assert_eq!(first, format!("{:?}", query()));
    assert_eq!(loaded, index.cache_stats());
    index.remove_cached_tree(target);
    assert_eq!(index.cache_stats(), (loaded.0, 0, 0));
    // Visit a different sequence, then return with a cold tree cache.
    index
        .query(targets[1], 0, 5000, false, None, None, false)
        .unwrap();
    index.clear_cached_trees();
    assert_eq!(first, format!("{:?}", query()));
    assert_eq!(index.cache_stats().1, loaded.1);
    index.clear_cached_trees();
    assert_eq!(index.cache_stats().1, 0);
}

#[test]
fn simulated_in_memory_index_is_not_evicted() {
    use impg::alignment_record::AlignmentRecord;
    use impg::impg::Impg;
    use impg::impg_index::ImpgIndex;
    use impg::seqidx::SequenceIndex;
    // SIMULATED metadata-only fixture: no real alignment or result file.
    let mut sequences = SequenceIndex::new();
    let q = sequences.get_or_insert_id("SIMULATED-query", Some(100));
    let t = sequences.get_or_insert_id("SIMULATED-target", Some(100));
    let records = vec![(
        vec![AlignmentRecord {
            query_id: q,
            query_start: 0,
            query_end: 100,
            target_id: t,
            target_start: 0,
            target_end: 100,
            strand_and_data_offset: 0,
            data_bytes: 0,
        }],
        String::from("SIMULATED.paf"),
    )];
    let index = Impg::from_multi_alignment_records(&records, sequences, None, false).unwrap();
    let before = index.cache_stats();
    assert_eq!(before.1, 1);
    index.remove_cached_tree(t);
    index.clear_cached_trees();
    assert_eq!(index.cache_stats(), before);
    assert!(index.get_or_load_tree(t).is_some());
}
