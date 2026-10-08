//! Deterministic work-count fixtures for the real partition scheduler.
use super::*;
use crate::impg::{AdjustedInterval, QueryMetadata};
use crate::seqidx::SequenceIndex as NameIndex;
use crate::subset_filter::SubsetFilter;
use coitrees::BasicCOITree;
use std::sync::{Arc, Mutex};

type Window = (u32, i32, i32);
struct CountingIndex {
    names: NameIndex,
    hits: FxHashMap<Window, Vec<Window>>,
    calls: Mutex<Vec<(Window, i32)>>,
    fail: Option<Window>,
}
impl CountingIndex {
    fn new(lengths: &[usize]) -> Self {
        let mut names = NameIndex::new();
        for (id, &len) in lengths.iter().enumerate() {
            names.get_or_insert_id(&format!("sample{id}#0#chr1"), Some(len));
        }
        Self {
            names,
            hits: FxHashMap::default(),
            calls: Mutex::new(Vec::new()),
            fail: None,
        }
    }
    fn dispatch(
        &self,
        id: u32,
        start: i32,
        end: i32,
        masks: Option<&FxHashMap<u32, SortedRanges>>,
    ) -> io::Result<Vec<AdjustedInterval>> {
        let covered: i32 = masks
            .and_then(|m| m.get(&id))
            .map(|r| {
                r.iter()
                    .map(|&(s, e)| (end.min(e) - start.max(s)).max(0))
                    .sum()
            })
            .unwrap_or(0);
        self.calls
            .lock()
            .unwrap()
            .push(((id, start, end), end - start - covered));
        if self.fail == Some((id, start, end)) {
            return Err(io::Error::other("counting backend failure"));
        }
        Ok(self
            .hits
            .get(&(id, start, end))
            .into_iter()
            .flatten()
            .map(|&(id, s, e)| {
                let iv = Interval {
                    first: s,
                    last: e,
                    metadata: id,
                };
                (iv, Vec::new(), iv)
            })
            .collect())
    }
}
#[allow(unused_variables)]
impl ImpgIndex for CountingIndex {
    fn seq_index(&self) -> &NameIndex {
        &self.names
    }
    fn query(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
    ) -> io::Result<Vec<AdjustedInterval>> {
        unreachable!("partition must use transitive dispatch")
    }
    fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) -> io::Result<Vec<AdjustedInterval>> {
        unreachable!("partition must use transitive dispatch")
    }
    fn populate_cigar_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cache: &mut FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) {
    }
    fn query_transitive_dfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        min_output_length: Option<i32>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> io::Result<Vec<AdjustedInterval>> {
        self.dispatch(target_id, range_start, range_end, masked_regions)
    }
    fn query_transitive_bfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        min_output_length: Option<i32>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> io::Result<Vec<AdjustedInterval>> {
        self.dispatch(target_id, range_start, range_end, masked_regions)
    }
    fn get_or_load_tree(&self, target_id: u32) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        None
    }
    fn target_ids(&self) -> Vec<u32> {
        (0..self.names.len() as u32).collect()
    }
    fn remove_cached_tree(&self, target_id: u32) {}
    fn sequence_files(&self) -> &[String] {
        &[]
    }
}

fn run(index: &CountingIndex, starting_ids: &[u32], mode: &str, dfs: bool) -> io::Result<String> {
    run_with_min_missing(index, starting_ids, mode, dfs, 0)
}

fn run_with_min_missing(
    index: &CountingIndex,
    starting_ids: &[u32],
    mode: &str,
    dfs: bool,
    min_missing_size: i32,
) -> io::Result<String> {
    let dir = tempfile::tempdir()?;
    let seeds = dir.path().join("seeds.txt");
    std::fs::write(
        &seeds,
        starting_ids
            .iter()
            .map(|&id| format!("{}\n", index.names.get_name(id).unwrap()))
            .collect::<String>(),
    )?;
    let engine = EngineOpts {
        engine: crate::GfaEngine::Pggb,
        syng_gfa_mode: None,
        syng_params: None,
        syng_gfa_frequency_mask: crate::commands::syng2gfa::SyngGfaFrequencyMask::disabled(),
        pipeline: crate::commands::graph::GraphBuildConfig::default(),
        target_poa_lengths: vec![700],
        max_node_length: 100,
        poa_padding_fraction: 0.001,
        partition_size: None,
        crush_config: None,
        smooth_after_crush: None,
        graph_sort_pipeline: None,
    };
    partition_alignments(
        index,
        100,
        Some(seeds.to_str().unwrap()),
        mode,
        0,
        None,
        min_missing_size,
        0,
        dfs,
        1,
        0,
        0,
        "bed",
        Some(dir.path().to_str().unwrap()),
        None,
        None,
        false,
        false,
        false,
        false,
        &engine,
        false,
    )?;
    std::fs::read_to_string(dir.path().join("partitions.bed"))
}

#[test]
fn stale_queued_window_dispatches_even_when_it_emits_nothing() {
    for dfs in [false, true] {
        let mut index = CountingIndex::new(&[100, 100]);
        index.hits.insert((0, 0, 100), vec![(1, 0, 100)]);
        let bed = run(&index, &[0, 1], "haplotype", dfs).unwrap();
        assert_eq!(
            *index.calls.lock().unwrap(),
            vec![((0, 0, 100), 100), ((1, 0, 100), 0)]
        );
        assert_eq!(bed.lines().count(), 2);
        assert!(bed.lines().all(|l| l.ends_with("\t0")));
    }
}

#[test]
fn mostly_covered_selected_group_keeps_partial_residual_query_context() {
    for mode in ["sample", "haplotype", "total"] {
        let mut index = CountingIndex::new(&[100, 1000, 100]);
        index
            .hits
            .insert((0, 0, 100), vec![(1, 0, 990), (2, 0, 95)]);
        // Ten residual bp alone do not carry this context-dependent homolog.
        index.hits.insert((1, 900, 1000), vec![(2, 0, 100)]);
        let bed = run(&index, &[0], mode, true).unwrap();
        let calls = index.calls.lock().unwrap();
        assert_eq!(calls.len(), 11);
        assert_eq!(calls.iter().filter(|(_, core)| *core == 0).count(), 9);
        assert_eq!(calls.last(), Some(&((1, 900, 1000), 10)));
        assert!(bed.contains("sample1#0#chr1\t990\t1000\t1"));
        assert!(bed.contains("sample2#0#chr1\t95\t100\t1"));
        drop(calls);
        assert_eq!(index.dispatch(1, 990, 1000, None).unwrap().len(), 0);
    }
}

#[test]
fn covered_b_discovers_c_that_a_and_c_do_not_discover() {
    for dfs in [false, true] {
        let mut index = CountingIndex::new(&[100, 100, 100, 100]);
        index.hits.insert((0, 0, 100), vec![(1, 0, 100)]);
        index
            .hits
            .insert((1, 0, 100), vec![(2, 0, 100), (3, 0, 100)]);
        let bed = run(&index, &[0, 1], "haplotype", dfs).unwrap();
        assert_eq!(
            *index.calls.lock().unwrap(),
            vec![((0, 0, 100), 100), ((1, 0, 100), 0)]
        );
        assert!(bed.contains("sample2#0#chr1\t0\t100\t1"));
        assert!(bed.contains("sample3#0#chr1\t0\t100\t1"));
        // Model dropping B's discovery: source union is unchanged, grouping/recall is not.
        index.hits.remove(&(1, 0, 100));
        let without_discovery = run(&index, &[0], "haplotype", dfs).unwrap();
        let source_rows = |text: &str| {
            let mut rows: Vec<String> = text
                .lines()
                .map(|l| l.rsplit_once('\t').unwrap().0.to_owned())
                .collect();
            rows.sort();
            rows
        };
        assert_eq!(source_rows(&bed), source_rows(&without_discovery));
        let partition_for = |text: &str, name: &str| -> String {
            text.lines()
                .find(|l| l.starts_with(name))
                .unwrap()
                .rsplit_once('\t')
                .unwrap()
                .1
                .to_owned()
        };
        assert_eq!(
            partition_for(&bed, "sample2#"),
            partition_for(&bed, "sample3#")
        );
        assert_ne!(
            partition_for(&without_discovery, "sample2#"),
            partition_for(&without_discovery, "sample3#")
        );
    }
}

#[test]
fn repeated_distinct_occurrences_survive_and_duplicate_windows_repeat_work() {
    let mut index = CountingIndex::new(&[100, 300]);
    index
        .hits
        .insert((0, 0, 100), vec![(1, 0, 100), (1, 300, 200)]);
    let bed = run(&index, &[0, 0], "haplotype", true).unwrap();
    let calls = index.calls.lock().unwrap();
    assert_eq!(calls.len(), 5);
    assert_eq!(calls[1], ((0, 0, 100), 0));
    assert_eq!(calls.iter().filter(|(_, core)| *core == 0).count(), 3);
    assert!(bed.contains("sample1#0#chr1\t0\t100\t0"));
    assert!(bed.contains("sample1#0#chr1\t200\t300\t0"));
    assert!(bed.contains("sample1#0#chr1\t100\t200\t1"));
}

#[test]
fn covered_window_backend_errors_remain_errors() {
    let mut index = CountingIndex::new(&[100, 100]);
    index.hits.insert((0, 0, 100), vec![(1, 0, 100)]);
    index.fail = Some((1, 0, 100));
    assert_eq!(
        run(&index, &[0, 1], "haplotype", true)
            .unwrap_err()
            .to_string(),
        "counting backend failure"
    );
    assert_eq!(index.calls.lock().unwrap().len(), 2);
}

#[test]
fn minimum_missing_size_absorbs_small_residual_without_omitting_source() {
    let mut index = CountingIndex::new(&[100, 100]);
    index.hits.insert((0, 0, 100), vec![(1, 0, 95)]);
    let bed = run_with_min_missing(&index, &[0], "haplotype", true, 6).unwrap();
    assert_eq!(*index.calls.lock().unwrap(), vec![((0, 0, 100), 100)]);
    assert_eq!(bed.lines().count(), 2);
    assert!(bed.contains("sample1#0#chr1\t0\t100\t0"));
}
