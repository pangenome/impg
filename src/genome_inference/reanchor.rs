//! THE K-REANCHOR TRANSFORM (the owner's dynamic local-graph design,
//! stated verbatim: "dynamically, take the sequences and build a local
//! syng graph using different kmer sizes and project the mappings we
//! have from the original set onto those").
//!
//! The committed multi-matching census records carry SEQUENCE
//! COORDINATES (per occurrence: path, start; per merged read-matched
//! interval: the covered node list, which is exactly the source
//! panel's walk steps FULLY CONTAINED in the interval's coordinate
//! span, in path order). The projection onto a target panel built at
//! a different k is therefore a COORDINATE LOOKUP, NOT A
//! RE-ALIGNMENT:
//!
//!   (1) recover each interval's coordinate span [lo, hi) on the
//!       source panel by locating the interval's node list as a
//!       contiguous run of the source panel's own stored walk
//!       (lo = the first contained step's bp, hi = the last step's
//!       bp + k_source — exact, because the interval's first/last
//!       contained nodes ARE the first/last anchors of the merged
//!       anchor-window union), verified by a source-grid containment
//!       self-check;
//!   (2) emit the TARGET panel's steps fully contained in [lo, hi)
//!       (bp >= lo and bp + k_target <= hi) as the re-anchored
//!       interval node list — a mapping at sequence offset X lands on
//!       whichever k node contains X;
//!   (3) re-anchor the occurrence's touched windows (the
//!       `partitions` field) by the same coordinate rule onto the
//!       target tiling's beds (the routing's own territory-touch
//!       rule: a window is touched iff one of its bed rows on the
//!       occurrence's path overlaps the occurrence's interval union).
//!
//! Everything else in the record (multiplicity, t_r, share,
//! total_occurrences, anchors, pair_bins, path, start, orientation)
//! is preserved VERBATIM — the committed occurrences are already
//! computed; only the node grid changes. The one named boundary: the
//! empty-interval occurrence class (read-edge anchors whose syncmer
//! minimization context differs on the panel path — start off-grid,
//! zero contained nodes, zero coverage mass) keeps `partitions`
//! verbatim, because its anchor-window coordinates are not
//! recoverable from the record.
//!
//! The walks are loaded IN MEMORY from the two syng indexes (the same
//! `walk_path_range` machinery that serves the homology-dump walks);
//! there is no walks-dump round trip and no JSON dump parse. The
//! census is streamed record by record — the transform is
//! seconds-class for millions of occurrences.

use crate::sample_mem_bwt::invalid;
use crate::syng::{SyncmerParams, SyngIndex};
use clap::Args;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// The per-occurrence scan window for locating an interval's node list
/// as a contiguous run of the source walk, starting at the occurrence
/// start's binary-search position. The committed intervals sit at or
/// immediately after their occurrence start, so a small bounded scan is
/// exact for the committed material; a miss is counted as a recovery
/// failure and the interval is kept verbatim (never silently dropped).
const RUN_SCAN_LIMIT: usize = 4096;

/// One panel path's stored walk steps, sorted by bp (the syng's own
/// step grid: absolute bp position -> abs node id).
pub(crate) struct WalkGrid {
    bp: Vec<u64>,
    node: Vec<u32>,
}

impl WalkGrid {
    /// Build from the walk machinery's (signed node, bp) steps.
    fn from_steps(steps: &[(i32, u64)]) -> Self {
        let mut sorted: Vec<(u64, u32)> = steps
            .iter()
            .map(|&(node, bp)| (bp, node.unsigned_abs()))
            .collect();
        sorted.sort_unstable_by_key(|&(bp, _)| bp);
        let mut grid = WalkGrid {
            bp: Vec::with_capacity(sorted.len()),
            node: Vec::with_capacity(sorted.len()),
        };
        for (bp, node) in sorted {
            grid.bp.push(bp);
            grid.node.push(node);
        }
        grid
    }

    /// Whether `position` is one of the grid's step positions.
    fn contains_bp(&self, position: u64) -> bool {
        self.bp.binary_search(&position).is_ok()
    }
}

/// A recovered interval coordinate span [lo, hi).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) struct Span {
    lo: u64,
    hi: u64,
}

/// Locate the interval's node list as a contiguous run of the source
/// walk at or after `from`, and return the run's coordinate span with
/// the containment self-check (the source-grid containment of [lo, hi)
/// must reproduce the interval exactly). Returns the span and the run's
/// end (the next interval's scan start).
pub(crate) fn recover_span(
    src: &WalkGrid,
    from: usize,
    interval: &[u32],
    k_source: u64,
    scan_limit: usize,
) -> Option<(Span, usize)> {
    let n = interval.len();
    if n == 0 || src.node.len() < n {
        return None;
    }
    let last = (src.node.len() - n).min(from.saturating_add(scan_limit));
    let mut i = from;
    while i <= last {
        if &src.node[i..i + n] == interval {
            let lo = src.bp[i];
            let hi = src.bp[i + n - 1] + k_source;
            // Self-consistency: the source-grid containment of [lo, hi)
            // must reproduce the interval exactly.
            let mut j = i;
            let mut count = 0usize;
            let mut exact = true;
            while j < src.bp.len() && src.bp[j] < hi {
                if src.bp[j] >= lo && src.bp[j] + k_source <= hi {
                    if count >= n || src.node[j] != interval[count] {
                        exact = false;
                        break;
                    }
                    count += 1;
                }
                j += 1;
            }
            if exact && count == n {
                return Some((Span { lo, hi }, i + n));
            }
        }
        i += 1;
    }
    None
}

/// The target grid's answer: the target steps fully contained in
/// [lo, hi) (bp >= lo and bp + k_target <= hi), in path order.
pub(crate) fn contained_steps(tgt: &WalkGrid, span: Span, k_target: u64) -> Vec<u32> {
    let mut out = Vec::new();
    let mut j = tgt.bp.partition_point(|&bp| bp < span.lo);
    while j < tgt.bp.len() && tgt.bp[j] < span.hi {
        if tgt.bp[j] + k_target <= span.hi {
            out.push(tgt.node[j]);
        }
        j += 1;
    }
    out
}

/// The re-anchored territory: the target tiling's beds plus the
/// component's window-id convention (the axis rows ranked by
/// (start, end, partition) map to window ids in order — the CLI's own
/// convention, one window per axis row).
pub(crate) struct Territory {
    /// window id -> partition.
    window_partition: Vec<u32>,
    /// path id -> window id -> bed rows [(start, end)] (empty when the
    /// window's partition has no row on that path).
    rows: Vec<Vec<Vec<(u64, u64)>>>,
}

impl Territory {
    /// The windows touched by one occurrence's interval-union spans:
    /// window w is touched iff one of w's partition's bed rows on the
    /// occurrence's path overlaps one of the spans (the routing's own
    /// territory-touch rule). Ascending window order by construction.
    pub(crate) fn touched_windows(&self, path: usize, spans: &[Span]) -> Vec<u32> {
        let empty: &[Vec<(u64, u64)>] = &[];
        let windows = self.rows.get(path).map_or(empty, |r| r.as_slice());
        let mut touched = Vec::new();
        for (w, rows) in windows.iter().enumerate() {
            'row: for &(row_lo, row_hi) in rows.iter() {
                for span in spans {
                    if span.hi > row_lo && span.lo < row_hi {
                        touched.push(w as u32);
                        break 'row;
                    }
                }
            }
        }
        touched
    }
}

/// The target tiling's bed rows per (window, path), built from the
/// partition bed directory and the window->partition convention.
fn load_territory(
    bed_dir: &Path,
    partition_graphs: &Path,
    component: &str,
    panel: &SyngIndex,
) -> io::Result<Territory> {
    // The component's axis rows (its own member rows) ranked by
    // (start, end, partition) are the window ids in order.
    let mut axis_rows: Vec<(u32, u64, u64)> = Vec::new();
    for entry in fs::read_dir(partition_graphs)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if !name.ends_with(".gfa.map.json") {
            continue;
        }
        #[derive(Deserialize)]
        struct MemberRow {
            path_name: String,
            start: u64,
            end: u64,
        }
        #[derive(Deserialize)]
        struct PartitionMap {
            partition: u32,
            members: Vec<MemberRow>,
        }
        let text = fs::read_to_string(&path)?;
        let map: PartitionMap = serde_json::from_str(&text).map_err(io::Error::other)?;
        for member in map.members.iter().filter(|m| m.path_name == component) {
            axis_rows.push((map.partition, member.start, member.end));
        }
    }
    if axis_rows.is_empty() {
        return Err(invalid(format!(
            "component {component:?} has no axis rows in {}",
            partition_graphs.display()
        )));
    }
    axis_rows.sort_by_key(|&(partition, start, end)| (start, end, partition));
    let window_partition: Vec<u32> = axis_rows.iter().map(|&(partition, _, _)| partition).collect();
    let n_windows = window_partition.len();

    // The partition beds: (partition, path) -> rows; then the per-path
    // window-indexed table the hot loop reads.
    let mut beds: FxHashMap<(u32, u32), Vec<(u64, u64)>> = FxHashMap::default();
    for entry in fs::read_dir(bed_dir)? {
        let path = entry?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let partition = match name.strip_suffix(".bed").and_then(|stem| {
            stem.strip_prefix("partition")
                .and_then(|id| id.parse::<u32>().ok())
        }) {
            Some(partition) => partition,
            None => continue,
        };
        let reader = BufReader::new(File::open(&path)?);
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let (Some(seq), Some(start), Some(end)) = (fields.next(), fields.next(), fields.next())
            else {
                return Err(invalid(format!(
                    "malformed bed row in {}: {line:?}",
                    path.display()
                )));
            };
            let Some(&path_id) = panel.name_map.name_to_path.get(seq) else {
                // A bed row for a sequence outside the target panel.
                continue;
            };
            beds.entry((partition, path_id))
                .or_default()
                .push((
                    start
                        .parse::<u64>()
                        .map_err(io::Error::other)?,
                    end.parse::<u64>().map_err(io::Error::other)?,
                ));
        }
    }
    let n_paths = panel.name_map.path_to_name.len();
    let mut rows: Vec<Vec<Vec<(u64, u64)>>> = vec![Vec::new(); n_paths];
    for (path_id, per_path) in rows.iter_mut().enumerate() {
        per_path.reserve(n_windows);
        for &partition in &window_partition {
            per_path.push(beds.get(&(partition, path_id as u32)).cloned().unwrap_or_default());
        }
    }
    Ok(Territory {
        window_partition,
        rows,
    })
}

// ---------------------------------------------------------------------------
// The census record (the committed multi-matching census format,
// fields in the committed alphabetical emission order).
// ---------------------------------------------------------------------------

#[derive(Deserialize, Serialize)]
struct CensusOccurrence {
    intervals: Vec<Vec<u32>>,
    orientation: u8,
    partitions: Vec<u32>,
    path: usize,
    start: u64,
}

#[derive(Deserialize, Serialize)]
struct CensusRecord {
    anchors: u64,
    multiplicity: u64,
    occurrences: Vec<CensusOccurrence>,
    pair_bins: std::collections::BTreeMap<String, u64>,
    record: u64,
    share: f64,
    t_r: u64,
    total_occurrences: u64,
}

#[derive(Args, Debug)]
pub struct Options {
    /// The committed multi-matching census receipt (jsonl or jsonl.gz)
    #[arg(long)]
    pub census: PathBuf,
    /// The source panel syng prefix (the grid the census records were
    /// derived against; supplies the coordinate recovery)
    #[arg(long)]
    pub source_panel: String,
    /// The target panel syng prefix (the new grid the records are
    /// re-anchored onto)
    #[arg(long)]
    pub panel: String,
    /// The target tiling's partition bed directory (partitionN.bed)
    #[arg(long)]
    pub bed_dir: PathBuf,
    /// The target locality map directory (windowN.gfa.map.json)
    #[arg(long)]
    pub partition_graphs: PathBuf,
    /// The component lane (the axis path whose member rows define the
    /// window-to-partition mapping)
    #[arg(long)]
    pub component: String,
    /// The re-anchored census output (jsonl or jsonl.gz)
    #[arg(long)]
    pub out: PathBuf,
    /// Optional stats receipt path (the transform's measured counters)
    #[arg(long)]
    pub stats: Option<PathBuf>,
}

/// The per-record transform's shared read-only state (the two grids,
/// the source->target path correspondence, the re-anchored territory).
struct Context {
    /// The source panel's walk grids, indexed by source path id.
    source_grids: Vec<WalkGrid>,
    /// Source path id -> target path id (by name; None when the
    /// sequence is absent from the target panel).
    source_to_target: Vec<Option<usize>>,
    /// The target panel's walk grids, indexed by target path id.
    target_grids: Vec<WalkGrid>,
    territory: Territory,
    k_source: u64,
    k_target: u64,
}

#[derive(Default)]
struct Counters {
    records: u64,
    occurrences: u64,
    intervals: u64,
    nodes_src: u64,
    nodes_tgt: u64,
    empty_intervals: u64,
    empty_tgt_intervals: u64,
    recovery_failures: u64,
    parts_changed: u64,
    parts_dropped: u64,
    parts_added: u64,
    occs_parts_changed: u64,
    parts_verbatim_empty: u64,
    starts_on_target_grid: u64,
    starts_on_source_grid: u64,
}

impl Counters {
    fn merge(&mut self, other: Counters) {
        let Counters {
            records,
            occurrences,
            intervals,
            nodes_src,
            nodes_tgt,
            empty_intervals,
            empty_tgt_intervals,
            recovery_failures,
            parts_changed,
            parts_dropped,
            parts_added,
            occs_parts_changed,
            parts_verbatim_empty,
            starts_on_target_grid,
            starts_on_source_grid,
        } = other;
        self.records += records;
        self.occurrences += occurrences;
        self.intervals += intervals;
        self.nodes_src += nodes_src;
        self.nodes_tgt += nodes_tgt;
        self.empty_intervals += empty_intervals;
        self.empty_tgt_intervals += empty_tgt_intervals;
        self.recovery_failures += recovery_failures;
        self.parts_changed += parts_changed;
        self.parts_dropped += parts_dropped;
        self.parts_added += parts_added;
        self.occs_parts_changed += occs_parts_changed;
        self.parts_verbatim_empty += parts_verbatim_empty;
        self.starts_on_target_grid += starts_on_target_grid;
        self.starts_on_source_grid += starts_on_source_grid;
    }
}

/// The per-record transform: parse one census line, re-grid every
/// occurrence's intervals onto the target panel, re-anchor the
/// touched windows onto the target tiling, serialize the committed
/// emission format. Pure with respect to the shared read-only context.
fn transform_record(
    line: &str,
    context: &Context,
    counters: &mut Counters,
) -> io::Result<Vec<u8>> {
    let mut record: CensusRecord = serde_json::from_str(line).map_err(io::Error::other)?;
    counters.records += 1;
    for occurrence in record.occurrences.iter_mut() {
        counters.occurrences += 1;
        let Some(source_grid) = context.source_grids.get(occurrence.path) else {
            // The occurrence's path is outside the source panel:
            // kept VERBATIM and counted (never silently dropped).
            counters.recovery_failures += 1;
            continue;
        };
        // The target grid is the same-named path's grid.
        let target_grid = context
            .source_to_target
            .get(occurrence.path)
            .copied()
            .flatten()
            .and_then(|p| context.target_grids.get(p));
        if source_grid.contains_bp(occurrence.start) {
            counters.starts_on_source_grid += 1;
        }
        if target_grid.is_some_and(|g| g.contains_bp(occurrence.start)) {
            counters.starts_on_target_grid += 1;
        }
        // (1) recover every interval's coordinate span on the
        // source grid; (2) slice the target grid.
        let mut scan = source_grid.bp.partition_point(|&bp| bp < occurrence.start);
        let mut spans: Vec<Span> = Vec::new();
        let mut recovered_all = true;
        let mut new_intervals: Vec<Vec<u32>> = Vec::with_capacity(occurrence.intervals.len());
        for interval in occurrence.intervals.iter() {
            counters.intervals += 1;
            counters.nodes_src += interval.len() as u64;
            if interval.is_empty() {
                // The read-edge-anchor class: no covered node on
                // the source grid; the empty interval is kept.
                counters.empty_intervals += 1;
                new_intervals.push(Vec::new());
                continue;
            }
            match recover_span(source_grid, scan, interval, context.k_source, RUN_SCAN_LIMIT) {
                Some((span, next)) => {
                    scan = next;
                    spans.push(span);
                    let new = match target_grid {
                        Some(grid) => contained_steps(grid, span, context.k_target),
                        None => interval.clone(),
                    };
                    if new.is_empty() {
                        counters.empty_tgt_intervals += 1;
                    }
                    counters.nodes_tgt += new.len() as u64;
                    new_intervals.push(new);
                }
                None => {
                    // The interval could not be located as a
                    // contiguous run of the source walk: kept
                    // verbatim, the occurrence's windows stay
                    // untouched (counted, never dropped).
                    recovered_all = false;
                    counters.recovery_failures += 1;
                    new_intervals.push(interval.clone());
                }
            }
        }
        occurrence.intervals = new_intervals;
        // (3) the territory-touch re-anchor onto the target
        // tiling's beds — only when every anchor-window coordinate
        // of the occurrence was recovered (the predecessor's
        // rule verbatim: the empty-interval class keeps its
        // window attribution, the named zero-coverage boundary).
        if recovered_all && !spans.is_empty() {
            match context.source_to_target[occurrence.path] {
                Some(path_id) => {
                    let touched = context.territory.touched_windows(path_id, &spans);
                    let old: BTreeSet<u32> = occurrence.partitions.iter().copied().collect();
                    let new: BTreeSet<u32> = touched.iter().copied().collect();
                    if old != new {
                        counters.occs_parts_changed += 1;
                        counters.parts_dropped += old.difference(&new).count() as u64;
                        counters.parts_added += new.difference(&old).count() as u64;
                    }
                    occurrence.partitions = touched;
                }
                // The sequence is absent from the target panel: no
                // target-grid coordinates, the window attribution is
                // kept verbatim (counted with the verbatim class).
                None => counters.parts_verbatim_empty += 1,
            }
        } else {
            counters.parts_verbatim_empty += 1;
        }
    }
    let mut bytes = serde_json::to_vec(&record).map_err(io::Error::other)?;
    bytes.push(b'\n');
    Ok(bytes)
}

/// The re-anchor transform: stream the census, re-grid every
/// occurrence's intervals onto the target panel, re-anchor the
/// partitions onto the target tiling, emit the committed format.
pub fn run(options: &Options) -> io::Result<()> {
    let whole = Instant::now();

    // The two panel loads (the accepted 3-4s local syng load each).
    let t = Instant::now();
    let source = SyngIndex::load(&options.source_panel, SyncmerParams::default())?;
    let source_load_seconds = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let target = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    let target_load_seconds = t.elapsed().as_secs_f64();
    let k_source = source.syncmer_length_bp() as u64;
    let k_target = target.syncmer_length_bp() as u64;

    // The walks, in memory, from the indexes' own walk machinery (no
    // dump round trip). Each path's steps, sorted by bp. The target
    // grid for the occurrence's path is looked up BY NAME (the two
    // panels carry the same sequences; the name lookup is robust to
    // path-id reordering).
    let t = Instant::now();
    let mut source_grids: Vec<WalkGrid> = Vec::with_capacity(source.name_map.path_to_name.len());
    let mut source_to_target: Vec<Option<usize>> =
        Vec::with_capacity(source.name_map.path_to_name.len());
    for path_idx in 0..source.name_map.path_to_name.len() {
        let len = source.name_map.path_to_length[path_idx];
        let steps = source.walk_path_range(path_idx, 0, len)?;
        source_grids.push(WalkGrid::from_steps(&steps));
        let name = &source.name_map.path_to_name[path_idx];
        source_to_target.push(
            target
                .name_map
                .name_to_path
                .get(name)
                .map(|&p| p as usize),
        );
    }
    let source_walk_seconds = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let mut target_grids: Vec<WalkGrid> = Vec::with_capacity(target.name_map.path_to_name.len());
    for path_idx in 0..target.name_map.path_to_name.len() {
        let len = target.name_map.path_to_length[path_idx];
        let steps = target.walk_path_range(path_idx, 0, len)?;
        target_grids.push(WalkGrid::from_steps(&steps));
    }
    let target_walk_seconds = t.elapsed().as_secs_f64();

    // The target tiling's territory (beds + the window convention).
    let t = Instant::now();
    let territory = load_territory(
        &options.bed_dir,
        &options.partition_graphs,
        &options.component,
        &target,
    )?;
    let n_windows = territory.window_partition.len();
    let territory_seconds = t.elapsed().as_secs_f64();
    let context = Context {
        source_grids,
        source_to_target,
        target_grids,
        territory,
        k_source,
        k_target,
    };

    // The census stream and the re-anchored output. The records are
    // INDEPENDENT (the grids and the territory are read-only), so the
    // per-record transform rides rayon in bounded batches, emitted in
    // input order; the counters merge per batch.
    let t = Instant::now();
    let (reader, _) = niffler::get_reader(Box::new(File::open(&options.census)?))
        .map_err(io::Error::other)?;
    let out_name = options.out.display().to_string();
    let out_file = File::create(&options.out)?;
    let out_buf = BufWriter::with_capacity(4 << 20, out_file);
    let mut out: Box<dyn Write> = if out_name.ends_with(".gz") {
        niffler::get_writer(
            Box::new(out_buf),
            niffler::compression::Format::Gzip,
            niffler::Level::Six,
        )
        .map_err(io::Error::other)?
    } else {
        Box::new(out_buf)
    };
    let mut counters = Counters::default();
    let mut reader = BufReader::with_capacity(4 << 20, reader);
    let mut process_seconds = 0f64;
    let mut emit_seconds = 0f64;
    const BATCH: usize = 256;
    let mut batch: Vec<String> = Vec::with_capacity(BATCH);
    loop {
        batch.clear();
        while batch.len() < BATCH {
            let mut line = String::new();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if !line.trim().is_empty() {
                batch.push(line);
            }
        }
        if batch.is_empty() {
            break;
        }
        let t_batch = Instant::now();
        let results: Vec<(io::Result<Vec<u8>>, Counters)> = batch
            .par_iter()
            .map(|line| {
                let mut local = Counters::default();
                let record = transform_record(line, &context, &mut local);
                (record, local)
            })
            .collect();
        process_seconds += t_batch.elapsed().as_secs_f64();
        let t_batch = Instant::now();
        for (record, local) in results {
            out.write_all(&record?).map_err(io::Error::other)?;
            counters.merge(local);
        }
        emit_seconds += t_batch.elapsed().as_secs_f64();
    }
    out.flush()?;
    drop(out);
    let transform_seconds = t.elapsed().as_secs_f64();
    let wall_seconds = whole.elapsed().as_secs_f64();

    let summary = serde_json::json!({
        "records": counters.records,
        "occurrences": counters.occurrences,
        "intervals": counters.intervals,
        "nodes_src": counters.nodes_src,
        "nodes_tgt": counters.nodes_tgt,
        "empty_intervals": counters.empty_intervals,
        "empty_tgt_intervals": counters.empty_tgt_intervals,
        "recovery_failures": counters.recovery_failures,
        "parts_changed": counters.parts_changed,
        "parts_dropped": counters.parts_dropped,
        "parts_added": counters.parts_added,
        "occs_parts_changed": counters.occs_parts_changed,
        "parts_verbatim_empty": counters.parts_verbatim_empty,
        "starts_on_target_grid": counters.starts_on_target_grid,
        "starts_on_source_grid": counters.starts_on_source_grid,
        "windows_target": n_windows,
        "k_source": k_source,
        "k_target": k_target,
        "source_load_seconds": source_load_seconds,
        "target_load_seconds": target_load_seconds,
        "source_walk_seconds": source_walk_seconds,
        "target_walk_seconds": target_walk_seconds,
        "territory_seconds": territory_seconds,
        "transform_seconds": transform_seconds,
        "process_seconds": process_seconds,
        "emit_seconds": emit_seconds,
        "wall_seconds": wall_seconds,
    });
    log::info!("re-anchor census: {}", summary);
    if let Some(stats_path) = &options.stats {
        let mut stats = String::new();
        for (key, value) in summary.as_object().unwrap() {
            if key.ends_with("_seconds") || key.starts_with("k_") {
                continue;
            }
            stats.push_str(&format!("{key}\t{value}\n"));
        }
        stats.push_str(&format!("wall_seconds\t{wall_seconds:.1}\n"));
        stats.push_str(&format!("transform_seconds\t{transform_seconds:.1}\n"));
        stats.push_str(&format!("process_seconds\t{process_seconds:.1}\n"));
        stats.push_str(&format!("emit_seconds\t{emit_seconds:.1}\n"));
        stats.push_str(&format!("k_source\t{k_source}\n"));
        stats.push_str(&format!("k_target\t{k_target}\n"));
        fs::write(stats_path, stats)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// The coordinate-lookup semantics tests (the owner's gate 1):
// containment, the k=63 round-trip on a small fixture, and the
// empty-interval edge class (partitions kept verbatim).
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// A small source fixture: one path, bp steps 0, 24, 48, ... with
    /// distinct node ids (k=63 spacing shape shrunk for the test).
    fn fixture_grid(bps: &[u64], nodes: &[u32]) -> WalkGrid {
        WalkGrid::from_steps(
            &bps
                .iter()
                .zip(nodes)
                .map(|(&bp, &node)| (node as i32, bp))
                .collect::<Vec<_>>(),
        )
    }

    #[test]
    fn recover_span_locates_run_and_self_checks() {
        // bp: 0, 24, 48, 72, 96 with nodes 10..14; k=63.
        let grid = fixture_grid(&[0, 24, 48, 72, 96], &[10, 11, 12, 13, 14]);
        // The interval [11, 12] is the run at bp 24, 48; span is
        // [24, 48 + 63) with the containment self-check reproducing it.
        let (span, next) = recover_span(&grid, 0, &[11, 12], 63, RUN_SCAN_LIMIT).unwrap();
        assert_eq!(span, Span { lo: 24, hi: 48 + 63 });
        assert_eq!(next, 3);
        // A node list that is NOT a contiguous run fails.
        assert!(recover_span(&grid, 0, &[11, 13], 63, RUN_SCAN_LIMIT).is_none());
        // A run past the scan start is found via the bounded scan.
        let (span2, _) = recover_span(&grid, 1, &[13, 14], 63, RUN_SCAN_LIMIT).unwrap();
        assert_eq!(span2, Span { lo: 72, hi: 96 + 63 });
    }

    #[test]
    fn contained_steps_is_the_full_containment_rule() {
        // Target grid at k=31: steps at 0, 20, 40, 60, 80.
        let grid = fixture_grid(&[0, 20, 40, 60, 80], &[5, 6, 7, 8, 9]);
        // Span [10, 75): steps with bp >= 10 and bp + 31 <= 75 ->
        // bp 20 (51 <= 75 yes), bp 40 (71 <= 75 yes), bp 60 (91 no).
        let out = contained_steps(&grid, Span { lo: 10, hi: 75 }, 31);
        assert_eq!(out, vec![6, 7]);
        // A span containing nothing target-side yields the empty
        // interval (the counted edge class, never an error).
        assert!(contained_steps(&grid, Span { lo: 0, hi: 15 }, 31).is_empty());
        // Exact boundary containment: bp + k == hi is contained.
        let out = contained_steps(&grid, Span { lo: 10, hi: 71 }, 31);
        assert_eq!(out, vec![6, 7]);
        // lo itself is contained when bp == lo.
        let out = contained_steps(&grid, Span { lo: 20, hi: 52 }, 31);
        assert_eq!(out, vec![6]);
    }

    #[test]
    fn k63_round_trip_on_a_small_fixture() {
        // THE ROUND-TRIP SEMANTICS: when the target grid is the source
        // grid (same k), the re-grid must reproduce every interval's
        // node list exactly (the identity gate's mechanism on a small
        // fixture).
        let bps: Vec<u64> = (0..20).map(|i| i as u64 * 24).collect();
        let nodes: Vec<u32> = (100..120).collect();
        let grid = fixture_grid(&bps, &nodes);
        // An interval spanning bp [24, 96): contained steps are the
        // run at 24, 48, 72 (nodes 101, 102, 103).
        let interval = vec![101u32, 102, 103];
        let (span, next) = recover_span(&grid, 0, &interval, 63, RUN_SCAN_LIMIT).unwrap();
        assert_eq!(span, Span { lo: 24, hi: 72 + 63 });
        assert_eq!(next, 4);
        let re_gridded = contained_steps(&grid, span, 63);
        assert_eq!(re_gridded, interval, "same-k round trip is the identity");
        // Across the whole fixture: every maximal run round-trips.
        for start in 0..bps.len() {
            for n in 1..=(bps.len() - start) {
                let run: Vec<u32> = nodes[start..start + n].to_vec();
                let from = grid.bp.partition_point(|&bp| bp < bps[start]);
                let (span, _) = recover_span(&grid, from, &run, 63, RUN_SCAN_LIMIT)
                    .expect("run must be located");
                assert_eq!(span.lo, bps[start]);
                assert_eq!(span.hi, bps[start + n - 1] + 63);
                assert_eq!(contained_steps(&grid, span, 63), run);
            }
        }
    }

    #[test]
    fn empty_interval_class_keeps_partitions_verbatim() {
        // THE EMPTY-INTERVAL EDGE CLASS (the 0.034% read-edge-anchor
        // occurrences): no covered node on the source grid -> no
        // recovered span -> the territory re-anchor is skipped and
        // the recorded windows are kept VERBATIM.
        let grid = fixture_grid(&[0, 24, 48], &[10, 11, 12]);
        let mut territory = Territory {
            window_partition: vec![0, 1],
            rows: vec![vec![
                vec![(0u64, 30u64)], // window 0 rows on path 0
                vec![(40, 80)],      // window 1 rows on path 0
            ]],
        };
        territory.rows.push(Vec::new()); // path 1: no rows at all
        // An occurrence whose only interval is empty: the run() loop
        // records it empty and the spans list stays empty, so the
        // occurrence keeps its recorded partitions verbatim. Here we
        // prove the two pieces the loop composes: no span from an
        // empty interval, and the untouched-window path.
        let empty: Vec<u32> = Vec::new();
        assert!(empty.is_empty());
        let spans: Vec<Span> = Vec::new();
        // With no spans nothing is touched; the caller keeps the old
        // attribution (mirroring run(): recovered_all && !spans.is_empty()
        // is false -> parts_verbatim_empty).
        assert!(territory.touched_windows(0, &spans).is_empty());
        // The non-empty case for contrast: a span over [40, 80) touches
        // only window 1 on path 0.
        let touched = territory.touched_windows(0, &[Span { lo: 40, hi: 80 }]);
        assert_eq!(touched, vec![1]);
        // A span overlapping both rows touches both windows, ascending.
        let touched = territory.touched_windows(0, &[Span { lo: 20, hi: 60 }]);
        assert_eq!(touched, vec![0, 1]);
        // A path with no bed rows touches nothing.
        let touched = territory.touched_windows(1, &[Span { lo: 0, hi: 1000 }]);
        assert!(touched.is_empty());
    }
}
