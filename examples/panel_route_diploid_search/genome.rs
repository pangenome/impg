//! Experimental genome-scale partition chaining.
//!
//! Public axis and source geometry define candidates and seams; sample data only
//! supplies counts. Chaining uses exact min-plus DP when the physical exit pair
//! is sufficient and an explicitly incomplete bounded conflict-window beam when
//! a larger exact legality ledger does not fit the approved experimental slice.
use super::partition::{
    add_record_subwalks_weighted, FeatureKey, Profile, ScoreModel, SourceRange,
};
use super::{ensure, invalid};
use impg::{
    genome_inference::{mem_records, panel_routes as routes},
    sample_mem_bwt::WeightedBwt,
    syng::SyngIndex,
};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fs::{self, File, OpenOptions},
    io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    sync::Arc,
    time::Instant,
};

/// Experimental genome-path ceiling: actual event-run MEM queries. The finite
/// prototype retains its unchanged 8M ceiling.
pub const MAX_PROFILE_WORK: u64 = 600_000_000;
pub const MAX_STATE_BYTES: u64 = 134_217_728;
pub const MAX_FEATURES: usize = 50_000;
pub const MAX_COMPLETE_PAIRS: usize = 2_048;
const CONFLICT_BEAM_WIDTH: usize = 128;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct AxisFile {
    pub version: u32,
    pub coordinate_system: String,
    pub intervals: Vec<AxisInterval>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct AxisInterval {
    pub component: String,
    pub start: u64,
    pub end: u64,
    pub group: String,
    pub reference_occurrence: usize,
    pub reference_strand: String,
    pub orientations: BTreeMap<String, String>,
}

pub fn load_axis(path: &Path) -> io::Result<AxisFile> {
    let axis: AxisFile =
        serde_json::from_reader(BufReader::new(File::open(path)?)).map_err(io::Error::other)?;
    ensure(
        axis.version == 1
            && axis.coordinate_system
                == "declared-S288C#0-assembly-source-axis-v1; no genotype prior or truth-conditioned candidate selection"
                || axis.version == 1
                    && axis.coordinate_system == "declared-S288C#0-assembly-source-axis-v1",
        "incompatible public reference axis",
    )?;
    ensure(!axis.intervals.is_empty(), "empty public reference axis")?;
    for (i, interval) in axis.intervals.iter().enumerate() {
        ensure(
            interval.start < interval.end
                && !interval.group.is_empty()
                && interval.reference_strand == "+",
            "invalid public reference-axis interval",
        )?;
        if i > 0 && axis.intervals[i - 1].component == interval.component {
            ensure(
                axis.intervals[i - 1].end == interval.start,
                "non-contiguous public reference axis",
            )?;
        }
    }
    Ok(axis)
}

/// Load only axis-selected public BED groups. BED3 is the catalog's public input;
/// null strand therefore retains both orientations. `max_occurrences` is a public
/// smoke-test subsample bound, not a truth-conditioned candidate selector. The
/// reference source is always retained when present.
pub fn load_bed_axis_partitions(
    axis: &AxisFile,
    bed_directory: &Path,
    lanes: &[(String, u64)],
    component: &str,
    max_partitions: Option<usize>,
    max_occurrences: Option<usize>,
) -> io::Result<(Vec<AxisInterval>, Vec<Vec<SourceRange>>)> {
    let source_by_name = lanes
        .iter()
        .enumerate()
        .map(|(source, (name, _))| (name.as_str(), source))
        .collect::<BTreeMap<_, _>>();
    let selected = axis
        .intervals
        .iter()
        .filter(|interval| interval.component == component)
        .take(max_partitions.unwrap_or(usize::MAX))
        .cloned()
        .collect::<Vec<_>>();
    ensure(!selected.is_empty(), "axis component has no partitions")?;
    let mut partitions = Vec::with_capacity(selected.len());
    for (partition, interval) in selected.iter().enumerate() {
        let bed = bed_directory.join(format!("{}.bed", interval.group));
        let mut physical = Vec::new();
        for (line_number, line) in BufReader::new(File::open(&bed)?).lines().enumerate() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            ensure(fields.len() == 3, "invalid public BED3 group")?;
            let source = *source_by_name
                .get(fields[0])
                .ok_or_else(|| invalid("BED source absent from route graph"))?;
            let start: u64 = fields[1].parse().map_err(io::Error::other)?;
            let end: u64 = fields[2].parse().map_err(io::Error::other)?;
            ensure(
                start < end && end <= lanes[source].1,
                "BED interval outside public source",
            )?;
            // BED line order is stable public catalog input identity. It is not
            // the catalog occurrence id, so use a collision-free scoped id.
            let occurrence = partition
                .checked_mul(10_000_000)
                .and_then(|x| x.checked_add(line_number))
                .ok_or_else(|| invalid("scoped occurrence id overflow"))?;
            for reverse in [false, true] {
                physical.push(SourceRange {
                    partition,
                    occurrence,
                    source,
                    start,
                    end,
                    reverse,
                });
            }
        }
        ensure(!physical.is_empty(), "empty public BED group")?;
        if let Some(limit) = max_occurrences {
            ensure(limit > 0, "zero occurrence subsample")?;
            // Keep complete forward/reverse pairs and deterministic source order.
            physical.sort_by_key(|r| (r.source, r.start, r.end, r.reverse));
            physical.truncate(limit.saturating_mul(2));
        }
        partitions.push(physical);
    }
    Ok((selected, partitions))
}

pub const MAX_COMPLETION_BRIDGE_BASES: u64 = 30_000;

/// The window-domain extension (owner-ruled: genome/window-domain-completeness.md,
/// Ruling 1 minimal universe extension + Ruling 2 Policy A row granularity).
///
/// A component run's per-window candidate rows are the anchor group's full BED
/// (unchanged bit-for-bit) PLUS every component-family row from any OTHER group
/// whose interval numerically overlaps the window's coordinate range
/// (the component's axis [start,end), the same family definition
/// `component_sources` uses), deduplicated by row identity (source,start,end).
///
/// Every contributing group keeps a single OWNING universe partition:
/// - a group that anchors a component window (dual-role) keeps that EXISTING
///   axis slot — the FIRST (lowest component-local locus) window it anchors
///   when several; and
/// - every other contributing group (pure-new) enters the COMPONENT routing
///   universe as an appended partition (id component_loci + ordinal, in
///   sorted group-name order) whose territory is EXACTLY its
///   window-overlapping component-family rows — derived from the axis window
///   ranges, no constants.
///
/// The genome routing universe is NOT extended (the genome-universe placement
/// pass of the genome-wide background stays exactly as is; gap records keep
/// t_r_genome = 0 — the ruled measured boundary).
pub struct WindowDomainExtension {
    /// Component-local axis locus count (the existing axis slots).
    pub component_loci: usize,
    /// Pure-new groups in sorted-name order with the DISTINCT
    /// component-family rows (source, start, end) overlapping some window of
    /// this component — exactly the appended universe partitions' territory.
    pub pure_new: Vec<(String, Vec<(usize, u64, u64)>)>,
    /// Dual-role group name -> its existing component-local axis slot.
    pub dual_role: BTreeMap<String, usize>,
    /// Per component locus: the window's ADDED rows as forward+reverse
    /// SourceRange pairs (the anchor loader's shape). `partition` carries the
    /// owner in the component-universe encoding: a component-local locus
    /// index for dual-role owners, `component_loci + ordinal` for pure-new.
    pub added: Vec<Vec<SourceRange>>,
    /// Row identities (source,start,end) found in MORE THAN ONE group BED on
    /// a component-family sequence (the panel tiling's disjointness check;
    /// the owner's no-duplication red line gate expects 0).
    pub duplicate_row_identities: u64,
    /// Component-family row pairs on one sequence that overlap but are not
    /// identical (the same disjointness check, adjacent-interval form after
    /// sorting; expected 0).
    pub overlapping_row_pairs: u64,
}

/// Scoped occurrence-id base for extension-added rows: far above the anchor
/// loader's `partition * 10_000_000 + line` block (axis partitions number in
/// the low thousands, so anchor ids stay below 10^11) and above the
/// completion machinery's per-locus 9_000_000+offset block.
const WINDOW_DOMAIN_ADDED_OCCURRENCE_BASE: usize = 100_000_000_000;

/// Scan the BED directory once and build the window-domain extension for one
/// component: the window-overlapping component-family rows per group, the
/// dual-role/pure-new ownership split, and the per-locus added-row lists.
/// The scan is bounded by the component's family sources (exact suffix match
/// already computed by the caller as `component_sources`) and the component's
/// axis window ranges; rows of other families are skipped without parsing.
pub fn build_window_domain_extension(
    axis: &AxisFile,
    bed_directory: &Path,
    lanes: &[(String, u64)],
    component: &str,
    component_sources: &BTreeSet<usize>,
) -> io::Result<WindowDomainExtension> {
    let source_by_name = lanes
        .iter()
        .enumerate()
        .map(|(source, (name, _))| (name.as_str(), source))
        .collect::<BTreeMap<_, _>>();
    let selected = axis
        .intervals
        .iter()
        .filter(|interval| interval.component == component)
        .cloned()
        .collect::<Vec<_>>();
    ensure(!selected.is_empty(), "axis component has no partitions")?;
    let component_loci = selected.len();
    let anchor_group_of_window = selected
        .iter()
        .map(|interval| interval.group.clone())
        .collect::<Vec<_>>();
    let anchor_groups = anchor_group_of_window
        .iter()
        .cloned()
        .collect::<BTreeSet<_>>();
    // Group -> rows (source, start, end), deduplicated by row identity; a
    // row identity belongs to exactly one group (the disjoint tiling) so the
    // per-group lists are disjoint by construction when the census is clean.
    let mut group_rows: BTreeMap<String, Vec<(usize, u64, u64)>> = BTreeMap::new();
    let mut duplicate_row_identities = 0u64;
    // Row identities seen across ALL groups (sorted scan order): the panel
    // tiling is disjoint per sequence, so a second group claiming the same
    // row identity is a tiling violation — counted for the no-duplication
    // gate and dropped from the later group (first sorted group wins
    // ownership, deterministically).
    let mut seen_rows: BTreeSet<(usize, u64, u64)> = BTreeSet::new();
    let mut bed_names = Vec::new();
    for entry in std::fs::read_dir(bed_directory)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|ext| ext.to_str()) == Some("bed") {
            bed_names.push(path.file_name().expect("file name").to_string_lossy().into_owned());
        }
    }
    bed_names.sort();
    for name in &bed_names {
        let group = name.strip_suffix(".bed").expect("bed suffix").to_string();
        for line in BufReader::new(File::open(bed_directory.join(name))?).lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() != 3 {
                continue;
            }
            let Some(&source) = source_by_name.get(fields[0]) else {
                continue;
            };
            if !component_sources.contains(&source) {
                continue;
            }
            let start: u64 = fields[1].parse().map_err(io::Error::other)?;
            let end: u64 = fields[2].parse().map_err(io::Error::other)?;
            ensure(
                start < end && end <= lanes[source].1,
                "BED interval outside public source",
            )?;
            // Only rows overlapping some window of this component contribute.
            if !selected
                .iter()
                .any(|interval| start < interval.end && end > interval.start)
            {
                continue;
            }
            if !seen_rows.insert((source, start, end)) {
                duplicate_row_identities += 1;
                continue;
            }
            let rows = group_rows.entry(group.clone()).or_default();
            rows.push((source, start, end));
            rows.sort_unstable();
        }
    }
    // Disjointness census: per component-family source, sorted adjacent
    // intervals must be disjoint (zero overlap) — the owner's red-line
    // structural verification, measured at extension-build time.
    let mut overlapping_row_pairs = 0u64;
    {
        let mut per_source: BTreeMap<usize, Vec<(u64, u64)>> = BTreeMap::new();
        for rows in group_rows.values() {
            for &(source, start, end) in rows {
                per_source.entry(source).or_default().push((start, end));
            }
        }
        for (_, mut rows) in per_source {
            rows.sort_unstable();
            for window in rows.windows(2) {
                if window[0].1 > window[1].0 {
                    overlapping_row_pairs += 1;
                }
            }
        }
    }
    // Ownership: dual-role groups keep their existing axis slot (the FIRST
    // window they anchor); pure-new groups get appended slots in sorted-name
    // order (BTreeMap iteration order).
    let mut dual_role = BTreeMap::new();
    let mut pure_new = Vec::new();
    for (group, rows) in &group_rows {
        if anchor_groups.contains(group) {
            let slot = anchor_group_of_window
                .iter()
                .position(|anchor| anchor == group)
                .expect("anchor group present");
            dual_role.insert(group.clone(), slot);
        } else {
            pure_new.push((group.clone(), rows.clone()));
        }
    }
    let pure_new_index = pure_new
        .iter()
        .enumerate()
        .map(|(ordinal, (group, _))| (group.clone(), ordinal))
        .collect::<BTreeMap<_, _>>();
    // Per-locus added rows: every contributing group's rows overlapping the
    // window, EXCLUDING the window's own anchor rows (already present),
    // owner-encoded, forward+reverse pairs, deterministic occurrence ids.
    let mut added = vec![Vec::<SourceRange>::new(); component_loci];
    for (group, rows) in &group_rows {
        let owner = if let Some(&slot) = dual_role.get(group) {
            slot
        } else {
            let ordinal = pure_new_index
                .get(group)
                .copied()
                .expect("contributing group classified");
            component_loci + ordinal
        };
        let ordinal = if let Some(&slot) = dual_role.get(group) {
            slot
        } else {
            pure_new_index.get(group).copied().expect("classified")
        };
        for (line_number, &(source, start, end)) in rows.iter().enumerate() {
            let occurrence = WINDOW_DOMAIN_ADDED_OCCURRENCE_BASE
                .checked_add(ordinal.checked_mul(10_000_000).expect("ordinal block"))
                .and_then(|value| value.checked_add(line_number))
                .ok_or_else(|| invalid("window-domain occurrence id overflow"))?;
            for interval in selected.iter().enumerate() {
                let (locus, window) = interval;
                if window.group == *group {
                    continue;
                }
                if start < window.end && end > window.start {
                    for reverse in [false, true] {
                        added[locus].push(SourceRange {
                            partition: owner,
                            occurrence,
                            source,
                            start,
                            end,
                            reverse,
                        });
                    }
                }
            }
        }
    }
    Ok(WindowDomainExtension {
        component_loci,
        pure_new,
        dual_role,
        added,
        duplicate_row_identities,
        overlapping_row_pairs,
    })
}

#[derive(Clone, Debug, Serialize)]
pub struct DomainCompletionException {
    pub source: usize,
    pub left_locus: usize,
    pub right_locus: usize,
    pub bridge_bases: u64,
    pub reason: String,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct DomainCompletionStats {
    pub bridge_candidates: usize,
    pub deletion_candidates: usize,
    pub additions_per_locus: Vec<usize>,
    pub exceptions: Vec<DomainCompletionException>,
}

/// Complete a component-family candidate chain from public source-path order.
///
/// BED groups can omit a panel source where that source uses alternate groups or
/// has a deletion relative to the reference-axis group. Between two observed,
/// forward occurrences of the same public source, attach the exact intervening
/// source span to the left observed locus and add zero-length pass-through
/// candidates at skipped axis loci. No cross-source coordinate projection is
/// performed: every added coordinate is an endpoint from that source's BED path.
pub fn complete_forward_source_paths(
    partitions: &mut [Vec<SourceRange>],
    family_sources: &BTreeSet<usize>,
) -> io::Result<DomainCompletionStats> {
    let mut stats = DomainCompletionStats {
        additions_per_locus: vec![0; partitions.len()],
        ..DomainCompletionStats::default()
    };
    for &source in family_sources {
        let observed = partitions
            .iter()
            .enumerate()
            .filter_map(|(locus, ranges)| {
                let mut forward = ranges
                    .iter()
                    .filter(|range| {
                        range.source == source && !range.reverse && range.start < range.end
                    })
                    .cloned()
                    .collect::<Vec<_>>();
                forward.sort_by_key(|range| (range.start, range.end, range.occurrence));
                (!forward.is_empty()).then_some((locus, forward))
            })
            .collect::<Vec<_>>();
        for pair in observed.windows(2) {
            let (left_locus, left_ranges) = &pair[0];
            let (right_locus, right_ranges) = &pair[1];
            if *right_locus <= *left_locus {
                continue;
            }
            let compatible = left_ranges
                .iter()
                .flat_map(|left| {
                    right_ranges
                        .iter()
                        .filter(move |right| left.end <= right.start)
                        .map(move |right| (left, right, right.start - left.end))
                })
                .min_by_key(|(left, right, gap)| (*gap, left.end, right.start));
            let Some((left, right, gap)) = compatible else {
                continue;
            };
            if gap > MAX_COMPLETION_BRIDGE_BASES {
                stats.exceptions.push(DomainCompletionException {
                    source,
                    left_locus: *left_locus,
                    right_locus: *right_locus,
                    bridge_bases: gap,
                    reason: "public source bridge exceeds profileable completion bound".into(),
                });
                continue;
            }
            if gap > 0 {
                let occurrence = left_locus
                    .checked_mul(10_000_000)
                    .and_then(|value| value.checked_add(9_000_000))
                    .and_then(|value| value.checked_add(source))
                    .ok_or_else(|| invalid("completed occurrence id overflow"))?;
                let bridge = SourceRange {
                    partition: *left_locus,
                    occurrence,
                    source,
                    start: left.start,
                    end: right.start,
                    reverse: false,
                };
                if !partitions[*left_locus].iter().any(|range| range == &bridge) {
                    partitions[*left_locus].push(bridge);
                    stats.bridge_candidates += 1;
                    stats.additions_per_locus[*left_locus] += 1;
                }
            }
            for locus in *left_locus + 1..*right_locus {
                let occurrence = locus
                    .checked_mul(10_000_000)
                    .and_then(|value| value.checked_add(9_000_000))
                    .and_then(|value| value.checked_add(source))
                    .ok_or_else(|| invalid("deletion occurrence id overflow"))?;
                let deletion = SourceRange {
                    partition: locus,
                    occurrence,
                    source,
                    start: right.start,
                    end: right.start,
                    reverse: false,
                };
                if !partitions[locus].iter().any(|range| range == &deletion) {
                    partitions[locus].push(deletion);
                    stats.deletion_candidates += 1;
                    stats.additions_per_locus[locus] += 1;
                }
            }
        }
    }
    Ok(stats)
}

/// Enforce production native endpoint pairing before profiling or scoring. This
/// is public graph geometry, never a sample-guided repair.
#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct SpanningTraversal {
    pub partition: usize,
    pub identity: String,
    pub segments: Vec<SourceRange>,
}

impl SpanningTraversal {
    pub fn single(range: SourceRange) -> Self {
        Self {
            partition: range.partition,
            identity: format!(
                "single:{}:{}:{}-{}:{}",
                range.occurrence, range.source, range.start, range.end, range.reverse
            ),
            segments: vec![range],
        }
    }
    fn entry(&self) -> &SourceRange {
        &self.segments[0]
    }
    fn exit(&self) -> &SourceRange {
        self.segments.last().unwrap()
    }
}

pub fn retain_native_endpoint_candidates(
    partitions: &mut [Vec<SpanningTraversal>],
    target: usize,
    target_length: u64,
) -> io::Result<()> {
    let last = partitions.len().saturating_sub(1);
    retain_native_endpoint_candidates_in_component_range(
        partitions,
        target,
        target_length,
        0,
        last,
    )
}

/// Endpoint legality measured against full-component locus positions, so a
/// locus-range slice does not demand native endpoints at its artificial edges.
pub fn retain_native_endpoint_candidates_in_component_range(
    partitions: &mut [Vec<SpanningTraversal>],
    target: usize,
    target_length: u64,
    first_component_locus: usize,
    last_component_locus: usize,
) -> io::Result<()> {
    ensure(!partitions.is_empty(), "empty endpoint candidate chain")?;
    for (offset, alleles) in partitions.iter_mut().enumerate() {
        let component_locus = first_component_locus + offset;
        alleles.retain(|allele| {
            let entry = allele.entry();
            let exit = allele.exit();
            let starts_legally = component_locus != 0
                || (entry.source == target && !entry.reverse && entry.start == 0);
            let ends_legally = component_locus != last_component_locus
                || (exit.source == target && !exit.reverse && exit.end == target_length);
            starts_legally && ends_legally
        });
        ensure(
            !alleles.is_empty(),
            "no production-endpoint candidate at axis locus",
        )?;
    }
    Ok(())
}

fn port_word_seams_with<F>(
    axis: &[AxisInterval],
    partitions: &[Vec<SpanningTraversal>],
    mut port_word: F,
) -> io::Result<Vec<Vec<(usize, usize)>>>
where
    F: FnMut(&SourceRange, bool) -> io::Result<Option<Vec<u8>>>,
{
    ensure(
        !partitions.is_empty() && partitions.len() == axis.len(),
        "axis/partition cardinality mismatch",
    )?;
    let mut result = Vec::with_capacity(partitions.len().saturating_sub(1));
    for boundary in 0..partitions.len().saturating_sub(1) {
        ensure(
            axis[boundary].component == axis[boundary + 1].component
                && axis[boundary].end == axis[boundary + 1].start,
            "seam requested across molecule boundary",
        )?;
        let mut right_by_word = BTreeMap::<Vec<u8>, Vec<usize>>::new();
        for (right_index, right) in partitions[boundary + 1].iter().enumerate() {
            if let Some(word) = port_word(right.entry(), false)? {
                right_by_word.entry(word).or_default().push(right_index);
            }
        }
        let mut links = Vec::new();
        for (left_index, left) in partitions[boundary].iter().enumerate() {
            // A deletion marker profiles its incoming seam using downstream
            // context. It may therefore leave only by exact same-source
            // continuity; a cross-source exit here would omit that exit seam.
            if left.exit().start < left.exit().end {
                if let Some(word) = port_word(left.exit(), true)? {
                    if let Some(rights) = right_by_word.get(&word) {
                        links.extend(rights.iter().map(|&right_index| (left_index, right_index)));
                    }
                }
            }
            // A contiguous same-source traversal normalizes to one segment and
            // therefore requires no production switch port.
            links.extend(
                partitions[boundary + 1]
                    .iter()
                    .enumerate()
                    .filter(|(_, right)| {
                        let left = left.exit();
                        let right = right.entry();
                        left.source == right.source
                            && left.reverse == right.reverse
                            && ((!left.reverse && left.end == right.start)
                                || (left.reverse && left.start == right.end))
                    })
                    .map(|(right_index, _)| (left_index, right_index)),
            );
        }
        links.sort_unstable();
        links.dedup();
        ensure(
            !links.is_empty(),
            "partition boundary has no public legal seam",
        )?;
        result.push(links);
    }
    Ok(result)
}

/// Public production seam legality: a contiguous same-source traversal, or an
/// oriented left exit port word exactly equal to the right entry port word. The
/// route artifact's verified per-lane port indexes supply both cuts and words.
pub fn port_word_seams(
    axis: &[AxisInterval],
    partitions: &[Vec<SpanningTraversal>],
    graph: &routes::Graph,
    ports: &mut routes::Ports,
) -> io::Result<Vec<Vec<(usize, usize)>>> {
    port_word_seams_with(axis, partitions, |range, exit| {
        let cut = match (range.reverse, exit) {
            (false, false) => range.start,
            (false, true) => range.end,
            (true, false) => range.end,
            (true, true) => range.start,
        };
        Ok(ports
            .at_cut(graph, range.source, cut, range.reverse)?
            .map(|port| port.word))
    })
}

/// Deterministic public-only smoke subsample: retain up to `max_paths`
/// same-source/orientation paths that cross every selected axis boundary.
fn legal_split_cut_order(
    left: &SourceRange,
    right: &SourceRange,
    left_cut: u64,
    right_cut: u64,
) -> bool {
    left.source != right.source || left_cut < right_cut
}

pub fn split_candidates(
    locus: usize,
    coarse: &[SpanningTraversal],
    admitted: &[usize],
    graph: &routes::Graph,
    ports: &mut routes::Ports,
) -> io::Result<Vec<SpanningTraversal>> {
    let mut port_sets = BTreeMap::<usize, Vec<routes::Port>>::new();
    for &index in admitted {
        ensure(index < coarse.len(), "split refinement class outside locus")?;
        let traversal = &coarse[index];
        ensure(
            traversal.segments.len() == 1 && !traversal.segments[0].reverse,
            "split refinement currently requires forward single-segment occurrences",
        )?;
        let range = &traversal.segments[0];
        port_sets.insert(
            index,
            ports.forward_ports_inside(graph, range.source, range.start, range.end)?,
        );
    }
    let mut result = BTreeMap::<String, SpanningTraversal>::new();
    for &left_index in admitted {
        let left = &coarse[left_index].segments[0];
        for &right_index in admitted {
            let right = &coarse[right_index].segments[0];
            let mut right_by_word = BTreeMap::<&[u8], Vec<u64>>::new();
            for port in &port_sets[&right_index] {
                right_by_word
                    .entry(port.word.as_slice())
                    .or_default()
                    .push(port.cut(graph.k));
            }
            for left_port in &port_sets[&left_index] {
                let Some(right_cuts) = right_by_word.get(left_port.word.as_slice()) else {
                    continue;
                };
                let left_cut = left_port.cut(graph.k);
                for &right_cut in right_cuts {
                    if !legal_split_cut_order(left, right, left_cut, right_cut) {
                        continue;
                    }
                    let identity = format!(
                        "split:{}:{}:{}-{}@{}>{}:{}-{}@{}",
                        locus,
                        left.occurrence,
                        left.start,
                        left.end,
                        left_cut,
                        right.occurrence,
                        right.start,
                        right.end,
                        right_cut
                    );
                    result
                        .entry(identity.clone())
                        .or_insert_with(|| SpanningTraversal {
                            partition: locus,
                            identity,
                            segments: vec![
                                SourceRange {
                                    end: left_cut,
                                    ..left.clone()
                                },
                                SourceRange {
                                    start: right_cut,
                                    ..right.clone()
                                },
                            ],
                        });
                }
            }
        }
    }
    Ok(result.into_values().collect())
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct SplitPairPrefilterStats {
    pub left_occurrence: usize,
    pub right_occurrence: usize,
    pub total_cuts: usize,
    pub retained: Vec<(usize, String, f64)>,
    pub all_ranked: Vec<(usize, String, f64)>,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct SplitPrefilterStats {
    pub occurrence_pairs: usize,
    pub total_cuts: usize,
    pub retained_cuts: usize,
    pub dropped_cuts: usize,
    pub retained_best_ties: usize,
    pub pairs: Vec<SplitPairPrefilterStats>,
}

/// Same-owner stitched spanning candidates (supervisor ruling,
/// genome/stitching-omission-alignment 2026-09-24): chains of ADJACENT
/// same-source, same-orientation, forward candidate rows whose union contains
/// the window's axis interval. The donor's route-boundary-fragmented pieces
/// (e.g. the SK1 chrIII rows [118626,128754) + [128754,138845) at the locus
/// window [120043,130060)) offer no single-row allele spanning the window
/// domain, while full-domain single rows of other sources win the tract
/// windows under the omission charge — the stitched form gives multi-row
/// sources the same spanning capability. Admission is derived, no constants:
/// adjacency is the panel tiling's own same-source continuity (the interior
/// seam between adjacent rows is the source's own contiguous DNA — the
/// co-occurring/pooled convention; gapped same-source pairs would be novel
/// junctions and are NOT admitted here), and the spanning test is the
/// window's own axis interval. Multi-SOURCE chains stay out of scope (the
/// owner-flagged narrowed capability); reverse orientations stay single-row
/// (the classing's mirror convention).
pub fn stitched_candidates(
    locus: usize,
    rows: &[SpanningTraversal],
    window_start: u64,
    window_end: u64,
) -> Vec<SpanningTraversal> {
    // Per source: the distinct forward non-empty rows keyed by interval, in
    // coordinate order (the BTreeSet dedups the anchor/added row duality);
    // the stitched segments CLONE the rows (owner partition and occurrence
    // identity preserved — the charging resolves owners per segment).
    let mut by_source: BTreeMap<usize, BTreeMap<(u64, u64), SourceRange>> = BTreeMap::new();
    for traversal in rows {
        if traversal.segments.len() != 1 {
            continue;
        }
        let segment = &traversal.segments[0];
        if segment.reverse || segment.start >= segment.end {
            continue;
        }
        by_source
            .entry(segment.source)
            .or_default()
            .entry((segment.start, segment.end))
            .or_insert_with(|| segment.clone());
    }
    let mut result = BTreeMap::<String, SpanningTraversal>::new();
    for (source, intervals) in by_source {
        let ordered: Vec<(u64, u64)> = intervals.keys().copied().collect();
        // Maximal runs of adjacent rows (left end == right start).
        let mut run_start = 0usize;
        for position in 1..=ordered.len() {
            let boundary = position == ordered.len()
                || ordered[position].0 != ordered[position - 1].1;
            if !boundary {
                continue;
            }
            // Every contiguous subchain [a..=b] of the run whose union
            // [ordered[a].0, ordered[b].1) contains the window.
            for first in run_start..position {
                if ordered[first].0 > window_start {
                    break;
                }
                let mut covering = None;
                for last in first..position {
                    if ordered[last].1 >= window_end {
                        covering = Some(last);
                        break;
                    }
                }
                let Some(last) = covering else { continue };
                for end in last..position {
                    // A single row is already a candidate — only chains of
                    // two or more segments are stitched forms.
                    if end == first {
                        continue;
                    }
                    let identity = format!(
                        "stitch:{}:{}:{}-{}:{}",
                        locus,
                        source,
                        ordered[first].0,
                        ordered[end].1,
                        end - first + 1,
                    );
                    result.entry(identity.clone()).or_insert_with(|| {
                        let segments: Vec<SourceRange> = ordered[first..=end]
                            .iter()
                            .map(|key| intervals[key].clone())
                            .collect();
                        let partition = segments[0].partition;
                        SpanningTraversal {
                            partition,
                            identity,
                            segments,
                        }
                    });
                }
            }
            run_start = position;
        }
    }
    result.into_values().collect()
}

pub fn split_candidates_prefiltered<F, S>(
    locus: usize,
    coarse: &[SpanningTraversal],
    admitted: &[usize],
    graph: &routes::Graph,
    ports: &mut routes::Ports,
    top_n: usize,
    report_all_ranks: bool,
    mut score: F,
) -> io::Result<(Vec<SpanningTraversal>, SplitPrefilterStats)>
where
    F: FnMut(&SpanningTraversal) -> io::Result<(f64, S)>,
    S: Clone + Eq,
{
    ensure(top_n > 0, "split seam prefilter top-N must be positive")?;
    let mut port_sets = BTreeMap::<usize, Vec<routes::Port>>::new();
    for &index in admitted {
        ensure(index < coarse.len(), "split refinement class outside locus")?;
        let traversal = &coarse[index];
        ensure(
            traversal.segments.len() == 1 && !traversal.segments[0].reverse,
            "split refinement currently requires forward single-segment occurrences",
        )?;
        let range = &traversal.segments[0];
        port_sets.insert(
            index,
            ports.forward_ports_inside(graph, range.source, range.start, range.end)?,
        );
    }
    let mut result = BTreeMap::<String, SpanningTraversal>::new();
    let mut stats = SplitPrefilterStats::default();
    for &left_index in admitted {
        let left = &coarse[left_index].segments[0];
        for &right_index in admitted {
            let right = &coarse[right_index].segments[0];
            let mut right_by_word = BTreeMap::<&[u8], Vec<u64>>::new();
            for port in &port_sets[&right_index] {
                right_by_word
                    .entry(port.word.as_slice())
                    .or_default()
                    .push(port.cut(graph.k));
            }
            let mut ranked = Vec::new();
            for left_port in &port_sets[&left_index] {
                let Some(right_cuts) = right_by_word.get(left_port.word.as_slice()) else {
                    continue;
                };
                let left_cut = left_port.cut(graph.k);
                for &right_cut in right_cuts {
                    if !legal_split_cut_order(left, right, left_cut, right_cut) {
                        continue;
                    }
                    let identity = format!(
                        "split:{}:{}:{}-{}@{}>{}:{}-{}@{}",
                        locus,
                        left.occurrence,
                        left.start,
                        left.end,
                        left_cut,
                        right.occurrence,
                        right.start,
                        right.end,
                        right_cut
                    );
                    let traversal = SpanningTraversal {
                        partition: locus,
                        identity,
                        segments: vec![
                            SourceRange {
                                end: left_cut,
                                ..left.clone()
                            },
                            SourceRange {
                                start: right_cut,
                                ..right.clone()
                            },
                        ],
                    };
                    let (loss, signature) = score(&traversal)?;
                    ranked.push((loss, signature, traversal));
                }
            }
            if ranked.is_empty() {
                continue;
            }
            stats.occurrence_pairs += 1;
            let total_here = ranked.len();
            stats.total_cuts += total_here;
            ranked.sort_by(|a, b| {
                a.0.total_cmp(&b.0)
                    .then_with(|| a.2.identity.cmp(&b.2.identity))
            });
            let best_signature = ranked[0].1.clone();
            let all_ranked = if report_all_ranks {
                ranked
                    .iter()
                    .enumerate()
                    .map(|(rank, (loss, _, traversal))| (rank, traversal.identity.clone(), *loss))
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            };
            let mut retained_here = 0usize;
            let mut retained = Vec::new();
            for (rank, (loss, signature, traversal)) in ranked.into_iter().enumerate() {
                let best_tie = signature == best_signature;
                if rank < top_n || best_tie {
                    if rank >= top_n && best_tie {
                        stats.retained_best_ties += 1;
                    }
                    retained_here += 1;
                    retained.push((rank, traversal.identity.clone(), loss));
                    result
                        .entry(traversal.identity.clone())
                        .or_insert(traversal);
                }
            }
            stats.retained_cuts += retained_here;
            stats.dropped_cuts += total_here - retained_here;
            stats.pairs.push(SplitPairPrefilterStats {
                left_occurrence: left.occurrence,
                right_occurrence: right.occurrence,
                total_cuts: total_here,
                retained,
                all_ranked,
            });
        }
    }
    Ok((result.into_values().collect(), stats))
}

pub fn connected_path_subsample(
    _axis: &[AxisInterval],
    partitions: &[Vec<SpanningTraversal>],
    max_paths: usize,
) -> io::Result<Vec<Vec<SpanningTraversal>>> {
    ensure(
        max_paths > 0 && !partitions.is_empty(),
        "invalid connected subsample",
    )?;
    let mut paths = Vec::new();
    for first in &partitions[0] {
        let mut path = vec![first.clone()];
        for partition in partitions.iter().skip(1) {
            let left = path.last().unwrap().exit();
            let next = partition.iter().find(|right| {
                let right = right.entry();
                left.source == right.source
                    && left.reverse == right.reverse
                    && ((!left.reverse && left.end == right.start)
                        || (left.reverse && left.start == right.end))
            });
            let Some(next) = next else {
                path.clear();
                break;
            };
            path.push(next.clone());
        }
        if path.len() == partitions.len() {
            paths.push(path);
            if paths.len() == max_paths {
                break;
            }
        }
    }
    ensure(!paths.is_empty(), "no public path crosses selected axis")?;
    let mut reduced = vec![Vec::new(); partitions.len()];
    for path in paths {
        for (partition, range) in path.into_iter().enumerate() {
            if !reduced[partition].contains(&range) {
                reduced[partition].push(range);
            }
        }
    }
    Ok(reduced)
}

#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProfileCost {
    /// Exact calls to the MEM extractor: one per invariant event run.
    pub mem_queries: u64,
    /// Complete read starts represented by those calls; diagnostic only.
    pub integrated_windows: u64,
}

fn add_cost(a: &mut ProfileCost, b: ProfileCost) -> io::Result<()> {
    a.mem_queries = a
        .mem_queries
        .checked_add(b.mem_queries)
        .ok_or_else(|| invalid("profile query count overflow"))?;
    a.integrated_windows = a
        .integrated_windows
        .checked_add(b.integrated_windows)
        .ok_or_else(|| invalid("integrated window count overflow"))?;
    Ok(())
}

pub(crate) fn event_boundaries(
    panel: &SyngIndex,
    sequence: &[u8],
    read_length: usize,
    start_lo: usize,
    start_hi: usize,
) -> io::Result<Vec<u64>> {
    ensure(
        read_length > 0
            && sequence.len() >= read_length
            && start_lo < start_hi
            && start_hi <= sequence.len() - read_length + 1,
        "invalid event profile start range",
    )?;
    let views = mem_records::raw_mem_anchor_positions(panel, sequence)?;
    let length = read_length as u64;
    let k = panel.syncmer_length_bp() as u64;
    let lo = start_lo as u64;
    let hi = start_hi as u64;
    let mut events = vec![lo, hi];
    for positions in &views {
        for &position in positions {
            for event in [
                (position + k).saturating_sub(length),
                position.saturating_add(1),
            ] {
                if event > lo && event < hi {
                    events.push(event);
                }
            }
        }
    }
    events.sort_unstable();
    events.dedup();
    Ok(events)
}

/// Count exact event-run MEM operations without executing MEM queries. Used for
/// fail-closed chromosome preflight against the unchanged 8M ceiling.
pub fn event_profile_cost(
    panel: &SyngIndex,
    sequence: &[u8],
    read_length: usize,
    start_lo: usize,
    start_hi: usize,
) -> io::Result<ProfileCost> {
    let events = event_boundaries(panel, sequence, read_length, start_lo, start_hi)?;
    Ok(ProfileCost {
        mem_queries: (events.len() - 1) as u64,
        integrated_windows: (start_hi - start_lo) as u64,
    })
}

/// Exact event-compressed full-subwalk profile over `[start_lo,start_hi)` read
/// starts. Enter/leave events hold both raw anchor views invariant, so one MEM
/// query represents every start in a run and is multiplied by that run length.
pub fn profile_event_runs(
    panel: &SyngIndex,
    sequence: &[u8],
    read_length: usize,
    start_lo: usize,
    start_hi: usize,
    max_features: usize,
) -> io::Result<(Profile, ProfileCost)> {
    ensure(
        max_features > 0 && max_features <= MAX_FEATURES,
        "invalid event profile limits",
    )?;
    let events = event_boundaries(panel, sequence, read_length, start_lo, start_hi)?;
    let mut profile = Profile::new();
    for run in events.windows(2) {
        let start = run[0] as usize;
        let multiplicity = run[1] - run[0];
        for record in
            mem_records::canonical_mem_records(panel, &sequence[start..start + read_length])?
        {
            add_record_subwalks_weighted(&mut profile, &record, multiplicity, max_features)?;
        }
    }
    Ok((
        profile,
        ProfileCost {
            mem_queries: (events.len() - 1) as u64,
            integrated_windows: (start_hi - start_lo) as u64,
        },
    ))
}

pub fn profile_event_interior(
    panel: &SyngIndex,
    sequence: &[u8],
    read_length: usize,
    max_features: usize,
) -> io::Result<(Profile, ProfileCost)> {
    if sequence.len() < read_length {
        return Ok((Profile::new(), ProfileCost::default()));
    }
    profile_event_runs(
        panel,
        sequence,
        read_length,
        0,
        sequence.len() - read_length + 1,
        max_features,
    )
}

pub fn profile_event_seam(
    panel: &SyngIndex,
    left: &[u8],
    right: &[u8],
    read_length: usize,
    max_features: usize,
) -> io::Result<(Profile, ProfileCost)> {
    if read_length <= 1 || left.is_empty() || right.is_empty() {
        return Ok((Profile::new(), ProfileCost::default()));
    }
    let left_take = left.len().min(read_length - 1);
    let right_take = right.len().min(read_length - 1);
    let mut context = left[left.len() - left_take..].to_vec();
    let boundary = context.len();
    context.extend_from_slice(&right[..right_take]);
    if context.len() < read_length {
        return Ok((Profile::new(), ProfileCost::default()));
    }
    let lo = boundary.saturating_add(1).saturating_sub(read_length);
    let hi = boundary.min(context.len() - read_length + 1);
    if lo >= hi {
        return Ok((Profile::new(), ProfileCost::default()));
    }
    profile_event_runs(panel, &context, read_length, lo, hi, max_features)
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct CacheRecord {
    key: String,
    entries: Vec<(FeatureKey, u64)>,
    cost: ProfileCost,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct CacheSeal {
    bytes: u64,
    fnv1a64: String,
    records: usize,
}

/// Dev-only fast-slice flags: set by the named variable or by the
/// IMPG_DEV_SLICE_FAST umbrella. These change reported accounting
/// (greedy_incumbent, work, local_pair_evidence) and, for
/// IMPG_DEV_TELESCOPED_EXTEND, the extend arithmetic itself; they are
/// authorized for the DEV slice loop ONLY. They must be explicitly off in
/// production and in any run whose accounting is quoted; a quoted dev run
/// must carry the flag name in its label.
fn dev_flag(name: &str) -> bool {
    std::env::var_os("IMPG_DEV_SLICE_FAST").is_some() || std::env::var_os(name).is_some()
}

fn fnv_update(mut value: u64, bytes: &[u8]) -> u64 {
    for &byte in bytes {
        value ^= byte as u64;
        value = value.wrapping_mul(0x100000001b3);
    }
    value
}

/// Append-only JSONL profile cache with a rewritten FNV seal. A key is either
/// present and reused exactly or absent and profiled once; corrupt/truncated
/// caches fail closed.
#[derive(Clone, Copy, Debug)]
struct CacheIndex {
    offset: u64,
    length: usize,
}

#[derive(Deserialize)]
struct CacheHeader {
    key: String,
}

pub struct ProfileCache {
    path: PathBuf,
    entries: BTreeMap<String, CacheIndex>,
    /// Total archive record count (unfiltered): the seal and sidecar describe
    /// the whole archive even when this open retained only a filtered subset.
    records: usize,
    bytes: u64,
    hash: u64,
    /// Sidecar offset index is loaded and maintained; future opens skip the
    /// whole-archive scan.
    indexed: bool,
    /// Binary records sidecar (key -> record payload) covering the archive
    /// up to `records_covered` bytes. None when absent or unrecoverably
    /// stale; lookups then take the JSONL path.
    records_sidecar: Option<RecordsSidecar>,
    /// Whole-archive JSONL bytes the binary records sidecar covers (0 when
    /// absent): reported so appends since the last sync stay visible.
    records_covered: u64,
    /// One-line diagnostics for the smoke's stage report.
    records_sidecar_status: String,
    records_sidecar_seconds: f64,
}

/// Binary records sidecar: a parse-free mirror of the JSONL archive's
/// records under `profiles.jsonl.rec`. Layout:
///
/// ```text
/// [record payloads...] [index block v1] [u64 ptr]
/// [record payloads...] [index block v2] [u64 ptr]   (after tail syncs)
/// ```
///
/// Each payload is a little-endian self-describing record (key, feature
/// entries, profile cost). An index block is a cumulative map from the
/// record's JSONL offset (unique per archive record) to its payload
/// location, prefixed by the archive byte/record counts the block covers;
/// the file's last u64 points at the newest complete block, so a crash
/// mid-sync leaves the previous block authoritative and the next open
/// re-syncs the tail. The JSONL archive remains the sole source of truth:
/// the sidecar is written only at open time, appends just make it stale,
/// and any structural defect falls back to the verified JSONL path.
struct RecordsSidecar {
    path: PathBuf,
    /// Sorted by JSONL record offset (records are appended in offset order).
    index: Vec<RecordSlot>,
    covered_bytes: u64,
    covered_records: u64,
}

#[derive(Clone, Copy, Debug)]
struct RecordSlot {
    jsonl_offset: u64,
    record_offset: u64,
    record_len: u32,
}

fn records_sidecar_path(path: &Path) -> PathBuf {
    path.with_extension("jsonl.rec")
}

fn put_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_u64(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn encode_record_binary(record: &CacheRecord) -> Vec<u8> {
    let mut out = Vec::with_capacity(16 + record.key.len() + 16 * record.entries.len());
    put_u32(&mut out, record.key.len() as u32);
    out.extend_from_slice(record.key.as_bytes());
    put_u32(&mut out, record.entries.len() as u32);
    for (key, count) in &record.entries {
        put_u32(&mut out, key.len() as u32);
        for word in key {
            put_u64(&mut out, *word);
        }
        put_u64(&mut out, *count);
    }
    put_u64(&mut out, record.cost.mem_queries);
    put_u64(&mut out, record.cost.integrated_windows);
    out
}

fn take_u32(payload: &[u8], cursor: &mut usize, what: &str) -> io::Result<u32> {
    let end = cursor.checked_add(4).filter(|end| end <= &payload.len());
    let Some(end) = end else {
        return Err(invalid(&format!("truncated records sidecar {what}")));
    };
    let value = u32::from_le_bytes(payload[*cursor..end].try_into().expect("four bytes"));
    *cursor = end;
    Ok(value)
}

fn take_u64(payload: &[u8], cursor: &mut usize, what: &str) -> io::Result<u64> {
    let end = cursor.checked_add(8).filter(|end| end <= &payload.len());
    let Some(end) = end else {
        return Err(invalid(&format!("truncated records sidecar {what}")));
    };
    let value = u64::from_le_bytes(payload[*cursor..end].try_into().expect("eight bytes"));
    *cursor = end;
    Ok(value)
}

fn take_bytes<'a>(payload: &'a [u8], cursor: &mut usize, len: usize, what: &str) -> io::Result<&'a [u8]> {
    let end = cursor.checked_add(len).filter(|end| end <= &payload.len());
    let Some(end) = end else {
        return Err(invalid(&format!("truncated records sidecar {what}")));
    };
    let value = &payload[*cursor..end];
    *cursor = end;
    Ok(value)
}

fn decode_record_binary(payload: &[u8], key: &str) -> io::Result<CacheRecord> {
    let mut cursor = 0usize;
    let key_len = take_u32(payload, &mut cursor, "key length")? as usize;
    let stored = take_bytes(payload, &mut cursor, key_len, "key")?;
    ensure(
        stored == key.as_bytes(),
        "records sidecar key mismatch",
    )?;
    let entry_count = take_u32(payload, &mut cursor, "entry count")? as usize;
    let mut entries = Vec::with_capacity(entry_count.min(1 << 20));
    for _ in 0..entry_count {
        let word_count = take_u32(payload, &mut cursor, "feature word count")? as usize;
        let mut feature = Vec::with_capacity(word_count.min(1 << 16));
        for _ in 0..word_count {
            feature.push(take_u64(payload, &mut cursor, "feature word")?);
        }
        let count = take_u64(payload, &mut cursor, "feature count")?;
        entries.push((feature, count));
    }
    let mem_queries = take_u64(payload, &mut cursor, "profile cost queries")?;
    let integrated_windows = take_u64(payload, &mut cursor, "profile cost windows")?;
    ensure(cursor == payload.len(), "records sidecar trailing bytes")?;
    Ok(CacheRecord {
        key: key.to_owned(),
        entries,
        cost: ProfileCost {
            mem_queries,
            integrated_windows,
        },
    })
}

/// Read and parse one JSONL archive record at a known offset.
fn read_jsonl_record(path: &Path, index: CacheIndex, key: &str) -> io::Result<CacheRecord> {
    let mut reader = File::open(path)?;
    reader.seek(SeekFrom::Start(index.offset))?;
    let mut line = vec![0; index.length];
    reader.read_exact(&mut line)?;
    let record: CacheRecord = serde_json::from_slice(&line).map_err(io::Error::other)?;
    ensure(record.key == key, "profile cache index mismatch")?;
    Ok(record)
}

/// Append one cumulative index block (coverage seal + slots by JSONL offset)
/// and the 8-byte pointer that makes it the newest complete block.
fn write_index_block(
    writer: &mut dyn Write,
    covered_bytes: u64,
    covered_records: u64,
    slots: &[RecordSlot],
    index_block_offset: u64,
) -> io::Result<()> {
    let mut block = Vec::with_capacity(24 + 20 * slots.len());
    put_u64(&mut block, covered_bytes);
    put_u64(&mut block, covered_records);
    put_u64(&mut block, slots.len() as u64);
    for slot in slots {
        put_u64(&mut block, slot.jsonl_offset);
        put_u64(&mut block, slot.record_offset);
        put_u32(&mut block, slot.record_len);
    }
    writer.write_all(&block)?;
    writer.write_all(&index_block_offset.to_le_bytes())?;
    Ok(())
}

impl RecordsSidecar {
    /// Locate a record slot by its unique JSONL record offset.
    fn find(&self, jsonl_offset: u64) -> Option<RecordSlot> {
        self.index
            .binary_search_by(|slot| slot.jsonl_offset.cmp(&jsonl_offset))
            .ok()
            .map(|position| self.index[position])
    }

    /// Read the newest complete index block. Any structural defect returns
    /// None so the caller falls back to a verified rebuild.
    fn load(path: &Path) -> Option<(Vec<RecordSlot>, u64, u64)> {
        let metadata = fs::metadata(path).ok()?;
        let mut file = File::open(path).ok()?;
        if metadata.len() < 8 {
            return None;
        }
        file.seek(SeekFrom::End(-8)).ok()?;
        let mut pointer = [0u8; 8];
        file.read_exact(&mut pointer).ok()?;
        let index_offset = u64::from_le_bytes(pointer);
        if index_offset >= metadata.len() - 8 {
            return None;
        }
        let index_len = (metadata.len() - 8 - index_offset) as usize;
        file.seek(SeekFrom::Start(index_offset)).ok()?;
        let mut block = vec![0; index_len];
        file.read_exact(&mut block).ok()?;
        let mut cursor = 0usize;
        let covered_bytes = take_u64(&block, &mut cursor, "index coverage bytes").ok()?;
        let covered_records = take_u64(&block, &mut cursor, "index coverage records").ok()?;
        let count = take_u64(&block, &mut cursor, "index slot count").ok()?;
        if count != covered_records || count > (1 << 24) {
            return None;
        }
        let mut slots = Vec::with_capacity(count as usize);
        for _ in 0..count {
            let jsonl_offset = take_u64(&block, &mut cursor, "slot jsonl offset").ok()?;
            let record_offset = take_u64(&block, &mut cursor, "slot record offset").ok()?;
            let record_len = take_u32(&block, &mut cursor, "slot record length").ok()?;
            if record_offset >= index_offset {
                return None;
            }
            slots.push(RecordSlot {
                jsonl_offset,
                record_offset,
                record_len,
            });
        }
        if cursor != block.len() || !slots.windows(2).all(|pair| pair[0].jsonl_offset < pair[1].jsonl_offset)
        {
            return None;
        }
        Some((slots, covered_bytes, covered_records))
    }
}

/// Sidecar offset-index seal: a snapshot of the profile-cache seal taken when
/// the index was last known to cover the archive exactly.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct CacheSidecarSeal {
    bytes: u64,
    fnv1a64: String,
    records: usize,
}

fn sidecar_path(path: &Path) -> PathBuf {
    path.with_extension("jsonl.idx")
}

fn sidecar_seal_path(path: &Path) -> PathBuf {
    path.with_extension("jsonl.idx.fnv")
}

impl ProfileCache {
    pub fn open(path: &Path) -> io::Result<Self> {
        Self::open_inner(path, &|_| true, false)
    }

    /// Open retaining only the keys accepted by `retain`. When a valid sidecar
    /// offset index exists, the archive JSONL is not scanned at all: a
    /// component-local slice run loads only its own subset of entries and
    /// seeks to those records on demand. Without a sidecar, one full scan
    /// builds it (atomically) and future opens are sidecar-fast. The retain
    /// predicate must accept every key this run can insert: a filtered-out key
    /// that is later inserted would append a duplicate archive record, which
    /// the next open detects and rejects (fail closed).
    pub fn open_filtered(path: &Path, retain: &dyn Fn(&str) -> bool) -> io::Result<Self> {
        Self::open_inner(path, retain, true)
    }

    pub fn retained_entries(&self) -> usize {
        self.entries.len()
    }

    pub fn loaded_from_sidecar(&self) -> bool {
        self.indexed
    }

    fn open_inner(
        path: &Path,
        retain: &dyn Fn(&str) -> bool,
        write_sidecar: bool,
    ) -> io::Result<Self> {
        Self::open_guarded(path, retain, write_sidecar, None)
    }

    /// Like [`Self::open_filtered`], with a hot-loop-sampled RSS budget for
    /// the one-time whole-archive scan (the scan allocates the full record
    /// index and previously had no checkpoint between stage boundaries).
    pub fn open_filtered_guarded(
        path: &Path,
        retain: &dyn Fn(&str) -> bool,
        rss_budget_bytes: Option<u64>,
    ) -> io::Result<Self> {
        Self::open_guarded(path, retain, true, rss_budget_bytes)
    }

    fn open_guarded(
        path: &Path,
        retain: &dyn Fn(&str) -> bool,
        write_sidecar: bool,
        rss_budget_bytes: Option<u64>,
    ) -> io::Result<Self> {
        let seal_path = path.with_extension("jsonl.fnv");
        let mut entries = BTreeMap::new();
        let mut byte_count = 0u64;
        let mut hash = 0xcbf29ce484222325u64;
        let mut indexed = false;
        let mut records = 0usize;
        let mut records_sidecar = None;
        let mut records_covered = 0u64;
        let mut records_sidecar_status = "absent".to_string();
        let mut records_sidecar_seconds = 0.0f64;
        if path.exists() || seal_path.exists() {
            ensure(
                path.exists() && seal_path.exists(),
                "incomplete profile cache",
            )?;
            let seal: CacheSeal =
                serde_json::from_reader(BufReader::new(File::open(&seal_path)?))
                    .map_err(io::Error::other)?;
            records = seal.records;
            if Self::load_sidecar(path, &seal, &mut entries)? {
                // Sidecar entries cover the archive exactly; the running hash
                // and byte count come from the matched seal snapshot.
                byte_count = seal.bytes;
                hash = u64::from_str_radix(&seal.fnv1a64, 16)
                    .map_err(io::Error::other)?;
                indexed = true;
            } else {
                let mut scan_guard = PeriodicRssGuard::new(rss_budget_bytes, 262_144);
                let mut reader = BufReader::new(File::open(path)?);
                let mut line = Vec::new();
                loop {
                    let offset = byte_count;
                    line.clear();
                    let length = reader.read_until(b'\n', &mut line)?;
                    if length == 0 {
                        break;
                    }
                    scan_guard.checkpoint("profile_cache_full_scan")?;
                    ensure(
                        line.last() == Some(&b'\n'),
                        "truncated profile cache record",
                    )?;
                    hash = fnv_update(hash, &line);
                    byte_count = byte_count
                        .checked_add(length as u64)
                        .ok_or_else(|| invalid("profile cache byte count overflow"))?;
                    let header: CacheHeader =
                        serde_json::from_slice(&line).map_err(io::Error::other)?;
                    ensure(
                        entries
                            .insert(header.key, CacheIndex { offset, length })
                            .is_none(),
                        "duplicate profile cache key",
                    )?;
                }
                ensure(
                    seal.bytes == byte_count && seal.fnv1a64 == format!("{hash:016x}"),
                    "profile cache seal mismatch",
                )?;
                ensure(
                    entries.len() == seal.records,
                    "profile cache record mismatch",
                )?;
                if write_sidecar {
                    Self::rebuild_sidecar(path, &seal, &entries)?;
                    // Appends in this process must keep the fresh sidecar in
                    // step with the archive.
                    indexed = true;
                }
            }
            // Parse-free records sidecar: brought in step with the sealed
            // archive before the filter is applied, so it covers the whole
            // archive for every future open regardless of filter.
            if seal.records > 0 {
                (
                    records_sidecar,
                    records_covered,
                    records_sidecar_status,
                    records_sidecar_seconds,
                ) = Self::sync_records_sidecar(path, &seal, &entries);
            }
            entries.retain(|key, _| retain(key));
        }
        Ok(Self {
            path: path.into(),
            entries,
            records,
            bytes: byte_count,
            hash,
            indexed,
            records_sidecar,
            records_covered,
            records_sidecar_status,
            records_sidecar_seconds,
        })
    }

    /// Binary records sidecar state for the smoke's stage report: (used,
    /// covered JSONL bytes, status, sync seconds).
    pub fn records_sidecar_report(&self) -> (bool, u64, &str, f64) {
        (
            self.records_sidecar.is_some(),
            self.records_covered,
            &self.records_sidecar_status,
            self.records_sidecar_seconds,
        )
    }

    /// Load the sidecar index when it provably covers the current archive:
    /// its seal snapshot must equal the archive seal and its entries must end
    /// exactly at the archive byte count. Returns false (and leaves `entries`
    /// untouched) whenever any check fails, so the caller falls back to the
    /// full verified scan.
    fn load_sidecar(
        path: &Path,
        seal: &CacheSeal,
        entries: &mut BTreeMap<String, CacheIndex>,
    ) -> io::Result<bool> {
        let index_path = sidecar_path(path);
        let seal_snapshot_path = sidecar_seal_path(path);
        if !index_path.exists() || !seal_snapshot_path.exists() {
            return Ok(false);
        }
        let snapshot: CacheSidecarSeal =
            serde_json::from_reader(BufReader::new(File::open(&seal_snapshot_path)?))
                .map_err(io::Error::other)?;
        if snapshot.bytes != seal.bytes
            || snapshot.records != seal.records
            || snapshot.fnv1a64 != seal.fnv1a64
        {
            return Ok(false);
        }
        // Any structural defect in the sidecar falls back to the verified
        // full archive scan; the sidecar is a pure accelerator.
        let mut loaded = BTreeMap::new();
        let mut end = 0u64;
        let mut count = 0usize;
        for line in BufReader::new(File::open(&index_path)?).lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let (Some(key), Some(offset), Some(length), None) =
                (fields.next(), fields.next(), fields.next(), fields.next())
            else {
                return Ok(false);
            };
            let (Ok(offset), Ok(length)) = (offset.parse::<u64>(), length.parse::<u64>()) else {
                return Ok(false);
            };
            if offset != end {
                return Ok(false);
            }
            end = offset
                .checked_add(length)
                .ok_or_else(|| invalid("profile cache sidecar byte count overflow"))?;
            let Ok(length) = usize::try_from(length) else {
                return Ok(false);
            };
            if loaded
                .insert(
                    key.to_owned(),
                    CacheIndex { offset, length },
                )
                .is_some()
            {
                return Ok(false);
            }
            count = count.saturating_add(1);
        }
        if count != seal.records || end != seal.bytes {
            return Ok(false);
        }
        *entries = loaded;
        Ok(true)
    }

    /// Bring the binary records sidecar in step with the sealed archive and
    /// return it for parse-free record lookups. The JSONL archive stays the
    /// source of truth: the sidecar is written only here (at open, before any
    /// filter is applied, so it covers the whole archive), a fresh tail is
    /// appended in place when the archive grew since the last sync, and a
    /// stale/corrupt sidecar is rebuilt atomically. Never fails the open:
    /// errors degrade to the JSONL path with a diagnostic status.
    fn sync_records_sidecar(
        path: &Path,
        seal: &CacheSeal,
        entries: &BTreeMap<String, CacheIndex>,
    ) -> (Option<RecordsSidecar>, u64, String, f64) {
        let started = Instant::now();
        let sidecar_path = records_sidecar_path(path);
        let mut slots: Vec<RecordSlot> = Vec::new();
        let mut covered_bytes = 0u64;
        let mut covered_records = 0u64;
        let mut mode;
        match RecordsSidecar::load(&sidecar_path) {
            Some((loaded, bytes, records)) => {
                if bytes == seal.bytes && records == seal.records as u64 {
                    return (
                        Some(RecordsSidecar {
                            path: sidecar_path,
                            index: loaded,
                            covered_bytes: bytes,
                            covered_records: records,
                        }),
                        bytes,
                        format!("loaded_{}_records", records),
                        started.elapsed().as_secs_f64(),
                    );
                } else if bytes < seal.bytes && records <= seal.records as u64 {
                    slots = loaded;
                    covered_bytes = bytes;
                    covered_records = records;
                    mode = "tail_synced";
                } else {
                    mode = "rebuilt_stale";
                }
            }
            None => {
                mode = if sidecar_path.exists() {
                    "rebuilt_corrupt"
                } else {
                    "rebuilt_absent"
                };
            }
        }
        // Rebuild from scratch discards any previously loaded slots.
        if mode.starts_with("rebuilt") {
            slots.clear();
            covered_bytes = 0;
            covered_records = 0;
        }
        // Records appended since the covered prefix, in JSONL offset order.
        let mut tail = entries
            .iter()
            .filter(|(_, index)| index.offset >= covered_bytes)
            .map(|(key, index)| (key.as_str(), *index))
            .collect::<Vec<_>>();
        tail.sort_by_key(|(_, index)| index.offset);
        if mode.starts_with("rebuilt")
            && (tail.len() as u64) != seal.records as u64
        {
            return (
                None,
                0,
                "rejected_index_seal_mismatch".into(),
                started.elapsed().as_secs_f64(),
            );
        }
        if !mode.starts_with("rebuilt")
            && (tail.len() as u64) != seal.records as u64 - covered_records
        {
            return (
                None,
                0,
                "rejected_tail_seal_mismatch".into(),
                started.elapsed().as_secs_f64(),
            );
        }
        if tail.is_empty() {
            // No archive growth to mirror: keep whatever loaded state is
            // coherent (a loaded sidecar returned earlier; here only the
            // empty-archive case remains).
            return (
                (mode.starts_with("rebuilt")).then(|| RecordsSidecar {
                    path: sidecar_path,
                    index: Vec::new(),
                    covered_bytes: 0,
                    covered_records: 0,
                }),
                0,
                format!("{mode}_empty"),
                started.elapsed().as_secs_f64(),
            );
        }
        // Write: rebuilds go through an atomic tmp+rename; tail syncs append
        // in place (payloads, then a new cumulative index block, then the
        // pointer — a crash before the pointer leaves the previous block
        // authoritative and the next open re-syncs the same tail). Payloads
        // are encoded and streamed one record at a time: a component-scale
        // rebuild must never hold the sidecar in memory.
        let write = || -> io::Result<(u64, Vec<RecordSlot>)> {
            let rebuild = mode.starts_with("rebuilt");
            let mut slots = if rebuild {
                Vec::with_capacity(tail.len())
            } else {
                slots
            };
            let mut record_offset = if rebuild {
                0
            } else {
                fs::metadata(&sidecar_path)?.len()
            };
            let writer: Box<dyn Write> = if rebuild {
                let temporary = sidecar_path.with_extension("rec.tmp");
                Box::new(BufWriter::new(File::create(&temporary)?))
            } else {
                Box::new(BufWriter::new(
                    OpenOptions::new().create(true).append(true).open(&sidecar_path)?,
                ))
            };
            let mut writer = writer;
            let mut written = 0u64;
            for (key, index) in &tail {
                let record = read_jsonl_record(path, *index, key)?;
                let encoded = encode_record_binary(&record);
                slots.push(RecordSlot {
                    jsonl_offset: index.offset,
                    record_offset,
                    record_len: encoded.len() as u32,
                });
                record_offset = record_offset.saturating_add(encoded.len() as u64);
                written = written.saturating_add(encoded.len() as u64);
                writer.write_all(&encoded)?;
            }
            let index_block_offset = if rebuild {
                written
            } else {
                record_offset
            };
            write_index_block(
                &mut writer,
                seal.bytes,
                seal.records as u64,
                &slots,
                index_block_offset,
            )?;
            writer.flush()?;
            drop(writer);
            if rebuild {
                // fsync before rename: the renamed file must be durable.
                let temporary = sidecar_path.with_extension("rec.tmp");
                File::open(&temporary)?.sync_all()?;
                fs::rename(&temporary, &sidecar_path)?;
            } else {
                File::open(&sidecar_path)?.sync_all()?;
            }
            Ok((record_offset, slots))
        };
        let slots = match write() {
            Ok((_, slots)) => slots,
            Err(_) => {
                return (
                    None,
                    0,
                    format!("{mode}_write_failed"),
                    started.elapsed().as_secs_f64(),
                )
            }
        };
        (
            Some(RecordsSidecar {
                path: sidecar_path,
                index: slots,
                covered_bytes: seal.bytes,
                covered_records: seal.records as u64,
            }),
            seal.bytes,
            format!("{}_{}_records", mode, tail.len()),
            started.elapsed().as_secs_f64(),
        )
    }

    fn rebuild_sidecar(
        path: &Path,
        seal: &CacheSeal,
        entries: &BTreeMap<String, CacheIndex>,
    ) -> io::Result<()> {
        let index_path = sidecar_path(path);
        let temporary = index_path.with_extension("jsonl.idx.tmp");
        {
            let mut writer = BufWriter::new(File::create(&temporary)?);
            let mut by_offset = entries
                .iter()
                .map(|(key, index)| (index.offset, index.length, key.as_str()))
                .collect::<Vec<_>>();
            by_offset.sort_unstable();
            for (offset, length, key) in by_offset {
                writeln!(writer, "{key}\t{offset}\t{length}")?;
            }
            writer.flush()?;
            writer.get_ref().sync_all()?;
        }
        fs::rename(&temporary, &index_path)?;
        let snapshot = CacheSidecarSeal {
            bytes: seal.bytes,
            fnv1a64: seal.fnv1a64.clone(),
            records: seal.records,
        };
        let snapshot_path = sidecar_seal_path(path);
        let snapshot_temporary = snapshot_path.with_extension("fnv.tmp");
        serde_json::to_writer(BufWriter::new(File::create(&snapshot_temporary)?), &snapshot)
            .map_err(io::Error::other)?;
        fs::rename(&snapshot_temporary, &snapshot_path)?;
        Ok(())
    }

    /// Read-only cache hit: parse the record if the key is indexed, None
    /// otherwise. The parallel interior loader uses this; misses fall back to
    /// the sequential mutable get_or_insert_with path (which may derive and
    /// append). The returned u64 is the record's parsed byte length (for
    /// parse-rate accounting): the binary payload length when the records
    /// sidecar served the record, the JSONL record length otherwise.
    pub fn get_if_cached(&self, key: &str) -> io::Result<Option<(Profile, ProfileCost, u64)>> {
        let Some(&index) = self.entries.get(key) else {
            return Ok(None);
        };
        let (record, bytes) = self
            .read_indexed_record(key, index)?
            .expect("indexed profile cache entry lacks a readable record");
        Ok(Some((record.entries.into_iter().collect(), record.cost, bytes)))
    }

    /// Fetch one indexed record through the binary records sidecar when it
    /// covers the record's JSONL offset, else through the JSONL itself.
    /// Both paths verify the record's key against the requested key.
    fn read_indexed_record(
        &self,
        key: &str,
        index: CacheIndex,
    ) -> io::Result<Option<(CacheRecord, u64)>> {
        if let Some(sidecar) = &self.records_sidecar {
            if let Some(slot) = sidecar.find(index.offset) {
                let mut reader = File::open(&sidecar.path)?;
                reader.seek(SeekFrom::Start(slot.record_offset))?;
                let mut payload = vec![0; slot.record_len as usize];
                reader.read_exact(&mut payload)?;
                let record = decode_record_binary(&payload, key)?;
                return Ok(Some((record, payload.len() as u64)));
            }
        }
        let mut reader = File::open(&self.path)?;
        reader.seek(SeekFrom::Start(index.offset))?;
        let mut line = vec![0; index.length];
        reader.read_exact(&mut line)?;
        let record: CacheRecord = serde_json::from_slice(&line).map_err(io::Error::other)?;
        ensure(record.key == key, "profile cache index mismatch")?;
        Ok(Some((record, index.length as u64)))
    }

    pub fn get_or_insert_with(
        &mut self,
        key: &str,
        profile: impl FnOnce() -> io::Result<(Profile, ProfileCost)>,
    ) -> io::Result<(Profile, ProfileCost, bool)> {
        if let Some(index) = self.entries.get(key).copied() {
            let (record, _) = self
                .read_indexed_record(key, index)?
                .expect("indexed profile cache entry lacks a readable record");
            return Ok((record.entries.into_iter().collect(), record.cost, true));
        }
        let (value, cost) = profile()?;
        ensure(
            value.len() <= MAX_FEATURES,
            "budget: maximal_mem_subwalk_features",
        )?;
        if let Some(parent) = self.path.parent() {
            fs::create_dir_all(parent)?;
        }
        let record = CacheRecord {
            key: key.into(),
            entries: value.iter().map(|(k, &v)| (k.clone(), v)).collect(),
            cost,
        };
        let mut line = serde_json::to_vec(&record).map_err(io::Error::other)?;
        line.push(b'\n');
        let mut writer = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.path)?;
        writer.write_all(&line)?;
        writer.sync_all()?;
        let index = CacheIndex {
            offset: self.bytes,
            length: line.len(),
        };
        self.entries.insert(key.to_owned(), index);
        self.records = self.records.saturating_add(1);
        self.bytes = self
            .bytes
            .checked_add(line.len() as u64)
            .ok_or_else(|| invalid("profile cache byte count overflow"))?;
        self.hash = fnv_update(self.hash, &line);
        if self.indexed {
            // Keep the sidecar offset index exactly in step with the archive:
            // append the new entry, then rewrite its seal snapshot to match the
            // rewritten archive seal. A failure here leaves the sidecar stale,
            // and the next open falls back to the verified full scan.
            let index_path = sidecar_path(&self.path);
            let mut writer = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&index_path)?;
            writeln!(writer, "{key}\t{}\t{}", index.offset, index.length)?;
            writer.sync_all()?;
            let snapshot = CacheSidecarSeal {
                bytes: self.bytes,
                fnv1a64: format!("{:016x}", self.hash),
                records: self.records,
            };
            let snapshot_path = sidecar_seal_path(&self.path);
            let snapshot_temporary = snapshot_path.with_extension("fnv.tmp");
            serde_json::to_writer(
                BufWriter::new(File::create(&snapshot_temporary)?),
                &snapshot,
            )
            .map_err(io::Error::other)?;
            fs::rename(&snapshot_temporary, &snapshot_path)?;
        }
        let seal = CacheSeal {
            bytes: self.bytes,
            fnv1a64: format!("{:016x}", self.hash),
            records: self.records,
        };
        let seal_path = self.path.with_extension("jsonl.fnv");
        let temporary = seal_path.with_extension("fnv.tmp");
        serde_json::to_writer(BufWriter::new(File::create(&temporary)?), &seal)
            .map_err(io::Error::other)?;
        fs::rename(temporary, seal_path)?;
        Ok((value, cost, false))
    }
}

#[derive(Clone, Debug)]
pub struct GenomeAllele {
    pub traversal: SpanningTraversal,
    pub profile: Profile,
}

#[derive(Clone, Debug)]
pub struct GenomeSeam {
    pub left: usize,
    pub right: usize,
    pub profile: Arc<Profile>,
}

pub fn intern_profile(
    profile: Profile,
    buckets: &mut HashMap<u64, Vec<Arc<Profile>>>,
) -> Arc<Profile> {
    let mut hash = 0xcbf29ce484222325u64;
    for (feature, count) in &profile {
        hash = fnv_update(hash, &(feature.len() as u64).to_le_bytes());
        for token in feature {
            hash = fnv_update(hash, &token.to_le_bytes());
        }
        hash = fnv_update(hash, &count.to_le_bytes());
    }
    if let Some(existing) = buckets
        .get(&hash)
        .and_then(|bucket| bucket.iter().find(|existing| existing.as_ref() == &profile))
    {
        return Arc::clone(existing);
    }
    let profile = Arc::new(profile);
    buckets.entry(hash).or_default().push(Arc::clone(&profile));
    profile
}

#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct BoundaryCompositionStats {
    pub physical_links: usize,
    pub class_compositions: usize,
    pub reused_physical_links: usize,
    pub per_physical_pair_compositions: usize,
}

pub fn deduplicated_boundary_profiles<F>(
    links: &[(usize, usize)],
    left_ends: &[Vec<u8>],
    right_starts: &[Vec<u8>],
    mut compose: F,
) -> io::Result<(Vec<Arc<Profile>>, BoundaryCompositionStats)>
where
    F: FnMut(&[u8], &[u8]) -> io::Result<Profile>,
{
    let mut classes = HashMap::<(u64, u64), Vec<(Vec<u8>, Vec<u8>, Arc<Profile>)>>::new();
    let mut output = Vec::with_capacity(links.len());
    let mut stats = BoundaryCompositionStats {
        physical_links: links.len(),
        ..BoundaryCompositionStats::default()
    };
    for &(left, right) in links {
        let left_end = &left_ends[left];
        let right_start = &right_starts[right];
        let key = (
            fnv_update(0xcbf29ce484222325, left_end),
            fnv_update(0xcbf29ce484222325, right_start),
        );
        if let Some(profile) = classes.get(&key).and_then(|bucket| {
            bucket
                .iter()
                .find(|(known_left, known_right, _)| {
                    known_left == left_end && known_right == right_start
                })
                .map(|(_, _, profile)| Arc::clone(profile))
        }) {
            stats.reused_physical_links = stats.reused_physical_links.saturating_add(1);
            output.push(profile);
            continue;
        }
        let profile = Arc::new(compose(left_end, right_start)?);
        stats.class_compositions = stats.class_compositions.saturating_add(1);
        classes.entry(key).or_default().push((
            left_end.clone(),
            right_start.clone(),
            Arc::clone(&profile),
        ));
        output.push(profile);
    }
    Ok((output, stats))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SeenAt {
    None,
    One(usize),
    Multiple,
}

impl SeenAt {
    fn observe(&mut self, at: usize) {
        *self = match *self {
            Self::None => Self::One(at),
            Self::One(previous) if previous == at => Self::One(previous),
            Self::One(_) | Self::Multiple => Self::Multiple,
        };
    }
}

#[derive(Clone, Copy, Debug)]
struct Incidence {
    witness: u64,
    interior: SeenAt,
    boundary: SeenAt,
    collision: bool,
}

#[derive(Clone, Debug)]
pub struct IncidenceTable {
    entries: HashMap<u64, Incidence>,
    pub collisions: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum LocalClass {
    Owned(usize),
    Boundary(usize),
    Deferred,
}

#[derive(Clone, Debug, Serialize)]
pub struct LocalPairEvidence {
    pub locus: usize,
    pub best_sources: [usize; 2],
    pub best_ranges: [(u64, u64, bool); 2],
    pub best_loss: f64,
    pub native_loss: Option<f64>,
    pub native_minus_best: Option<f64>,
    pub best_owned_observed: u64,
    pub best_owned_predicted: u64,
    pub native_owned_observed: Option<u64>,
    pub native_owned_predicted: Option<u64>,
    pub top_class_representatives: Vec<usize>,
    pub top_homozygous_losses: Vec<f64>,
    pub selected_loss: Option<f64>,
    pub selected_minus_best: Option<f64>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ChainAccounting {
    pub work: u64,
    pub profile_work: u64,
    pub integrated_windows: u64,
    pub state_bytes: u64,
    pub occupancy: Vec<usize>,
    pub tie_backpointers: u64,
    pub dropped_ties: u64,
    pub dropped_non_ties: u64,
    pub global_early_stop_class_pairs: u64,
    pub beam_width: usize,
    pub phase_cuts: Vec<usize>,
    pub exit_classes: Vec<usize>,
    pub conflict_window: usize,
    pub exact_blocks: usize,
    pub bounded_blocks: usize,
    pub incidence_hashes: usize,
    pub incidence_collisions: u64,
    pub sample_count_queries: u64,
    pub sample_count_cache_resets: u64,
    /// A*-bound pruning: states removed at a locus start because
    /// score + admissible future lower bound cannot beat the incumbent.
    #[serde(default)]
    pub bound_pruned_states: u64,
    /// A*-bound pruning: (state, successor-class) edges skipped for the same
    /// reason before any physical member materialization.
    #[serde(default)]
    pub bound_pruned_edges: u64,
    /// Greedy incumbent's completed score, when the greedy class chain
    /// reached the terminal locus; bound pruning is active only then.
    #[serde(default)]
    pub greedy_incumbent: Option<f64>,
    pub local_pair_evidence: Vec<LocalPairEvidence>,
}

#[derive(Clone, Debug, Serialize)]
pub struct GenomeChain {
    pub complete: bool,
    pub stop_reason: String,
    pub choices: Vec<[Vec<usize>; 2]>,
    pub proposal_losses: Vec<f64>,
    pub accounting: ChainAccounting,
    /// Best internal score per tie-rotation pass (main pass and each rotated
    /// pass), reported for tie-diversity diagnostics. Empty when the chain
    /// stopped before completing any pass (budget stops) or for the pure
    /// exact-DP path.
    pub rotation_best_scores: Vec<Option<f64>>,
}

#[derive(Clone, Debug, Serialize)]
pub struct GenomeAudit {
    pub exit_classes: Vec<usize>,
    pub conflict_window: usize,
    pub boundary_state_products: Vec<u128>,
}

fn profile_signature(profile: &Profile) -> Vec<(FeatureKey, u64)> {
    profile
        .iter()
        .map(|(key, &count)| (key.clone(), count))
        .collect()
}

pub fn profile_classes(alleles: &[GenomeAllele]) -> Vec<Vec<usize>> {
    let mut classes: BTreeMap<Vec<(FeatureKey, u64)>, Vec<usize>> = BTreeMap::new();
    for (index, allele) in alleles.iter().enumerate() {
        classes
            .entry(profile_signature(&allele.profile))
            .or_default()
            .push(index);
    }
    classes.into_values().collect()
}

fn class_for_alleles(alleles: &[GenomeAllele]) -> (Vec<Profile>, Vec<usize>) {
    class_for_profiles(
        &alleles
            .iter()
            .map(|allele| allele.profile.clone())
            .collect::<Vec<_>>(),
    )
}

fn class_for_profiles(physical: &[Profile]) -> (Vec<Profile>, Vec<usize>) {
    let mut classes: BTreeMap<Vec<(FeatureKey, u64)>, Vec<usize>> = BTreeMap::new();
    for (index, profile) in physical.iter().enumerate() {
        classes
            .entry(profile_signature(profile))
            .or_default()
            .push(index);
    }
    let mut profiles = Vec::with_capacity(classes.len());
    let mut membership = vec![0; physical.len()];
    for (class, members) in classes.into_values().enumerate() {
        profiles.push(physical[members[0]].clone());
        for member in members {
            membership[member] = class;
        }
    }
    (profiles, membership)
}

fn spans_feasible(ranges: &[SourceRange]) -> bool {
    for i in 0..ranges.len() {
        for j in i + 1..ranges.len() {
            if ranges[i].source == ranges[j].source
                && ranges[i].start < ranges[j].end
                && ranges[j].start < ranges[i].end
            {
                return false;
            }
        }
    }
    true
}

fn key_hash(key: &FeatureKey, seed: u64) -> u64 {
    let mut hash = seed ^ (key.len() as u64);
    for &token in key {
        for byte in token.to_le_bytes() {
            hash ^= byte as u64;
            hash = hash.wrapping_mul(0x100000001b3);
        }
    }
    hash
}

fn observe_incidence(
    table: &mut HashMap<u64, Incidence>,
    collisions: &mut u64,
    key: &FeatureKey,
    at: usize,
    interior: bool,
) {
    let primary = key_hash(key, 0xcbf29ce484222325);
    let witness = key_hash(key, 0x9e3779b97f4a7c15);
    let entry = table.entry(primary).or_insert(Incidence {
        witness,
        interior: SeenAt::None,
        boundary: SeenAt::None,
        collision: false,
    });
    if entry.witness != witness {
        if !entry.collision {
            *collisions += 1;
        }
        entry.collision = true;
        return;
    }
    if interior {
        entry.interior.observe(at);
    } else {
        entry.boundary.observe(at);
    }
}

pub fn incidence_table(
    partitions: &[Vec<GenomeAllele>],
    seams: &[Vec<GenomeSeam>],
) -> IncidenceTable {
    let mut entries = HashMap::new();
    let mut collisions = 0;
    for (locus, alleles) in partitions.iter().enumerate() {
        for allele in alleles {
            for key in allele.profile.keys() {
                observe_incidence(&mut entries, &mut collisions, key, locus, true);
            }
        }
    }
    for (boundary, joins) in seams.iter().enumerate() {
        for seam in joins {
            for key in seam.profile.keys() {
                observe_incidence(&mut entries, &mut collisions, key, boundary, false);
            }
        }
    }
    IncidenceTable {
        entries,
        collisions,
    }
}

impl IncidenceTable {
    pub fn empty() -> Self {
        Self {
            entries: HashMap::new(),
            collisions: 0,
        }
    }

    pub fn observe_interior(&mut self, locus: usize, profile: &Profile) {
        for key in profile.keys() {
            observe_incidence(&mut self.entries, &mut self.collisions, key, locus, true);
        }
    }

    pub fn observe_boundary(&mut self, boundary: usize, profile: &Profile) {
        for key in profile.keys() {
            observe_incidence(
                &mut self.entries,
                &mut self.collisions,
                key,
                boundary,
                false,
            );
        }
    }

    fn classify(&self, key: &FeatureKey) -> LocalClass {
        let primary = key_hash(key, 0xcbf29ce484222325);
        let witness = key_hash(key, 0x9e3779b97f4a7c15);
        let Some(entry) = self.entries.get(&primary) else {
            return LocalClass::Deferred;
        };
        if entry.collision || entry.witness != witness {
            return LocalClass::Deferred;
        }
        match (entry.interior, entry.boundary) {
            (SeenAt::One(locus), SeenAt::None) => LocalClass::Owned(locus),
            (SeenAt::None, SeenAt::One(boundary)) => LocalClass::Boundary(boundary),
            _ => LocalClass::Deferred,
        }
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn classification_label(&self, key: &FeatureKey) -> String {
        match self.classify(key) {
            LocalClass::Owned(locus) => format!("owned:{locus}"),
            LocalClass::Boundary(boundary) => format!("boundary:{boundary}"),
            LocalClass::Deferred => "deferred".into(),
        }
    }

    pub fn deferred(&self) -> usize {
        self.entries
            .values()
            .filter(|entry| {
                entry.collision
                    || !matches!(
                        (entry.interior, entry.boundary),
                        (SeenAt::One(_), SeenAt::None) | (SeenAt::None, SeenAt::One(_))
                    )
            })
            .count()
    }
}

struct CountCache<'a> {
    sample: &'a WeightedBwt,
    /// Compact hash shards keep each materialized map below 50k entries.
    values: Vec<HashMap<u64, (u64, u64)>>,
    entries: usize,
    queries: u64,
    resets: u64,
}
impl CountCache<'_> {
    fn get(&mut self, key: &FeatureKey) -> io::Result<u64> {
        let primary = key_hash(key, 0xcbf29ce484222325);
        let witness = key_hash(key, 0x9e3779b97f4a7c15);
        let shard = (primary as usize) & 255;
        if let Some(&(stored_witness, count)) = self.values[shard].get(&primary) {
            if stored_witness == witness {
                return Ok(count);
            }
            // A detected primary collision is queried directly and never cached.
            self.queries += 1;
            return self.sample.count(key);
        }
        if self.entries == 4_000_000 {
            for values in &mut self.values {
                values.clear();
            }
            self.entries = 0;
            self.resets += 1;
        }
        let count = self.sample.count(key)?;
        self.values[shard].insert(primary, (witness, count));
        self.entries += 1;
        self.queries += 1;
        Ok(count)
    }
}

type PreparedProfile = Vec<(FeatureKey, u64, u64)>;

fn prepare_profile(
    profile: &Profile,
    wanted: LocalClass,
    incidence: &IncidenceTable,
    counts: &mut CountCache<'_>,
) -> io::Result<PreparedProfile> {
    profile
        .iter()
        .filter(|(key, _)| incidence.classify(key) == wanted)
        .map(|(key, &predicted)| Ok((key.clone(), predicted, counts.get(key)?)))
        .collect()
}

fn score_prepared_pair(
    first: &PreparedProfile,
    second: &PreparedProfile,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = 0.0;
    let mut left = first.iter().peekable();
    let mut right = second.iter().peekable();
    while left.peek().is_some() || right.peek().is_some() {
        let (q, observed) = match (left.peek(), right.peek()) {
            (Some(a), Some(b)) if a.0 == b.0 => {
                let a = left.next().unwrap();
                let b = right.next().unwrap();
                ensure(a.2 == b.2, "inconsistent cached observation")?;
                (
                    a.1.checked_add(b.1)
                        .ok_or_else(|| invalid("diploid profile overflow"))?,
                    a.2,
                )
            }
            (Some(a), Some(b)) if a.0 < b.0 => {
                let a = left.next().unwrap();
                (a.1, a.2)
            }
            (Some(_), Some(_)) => {
                let b = right.next().unwrap();
                (b.1, b.2)
            }
            (Some(_), None) => {
                let a = left.next().unwrap();
                (a.1, a.2)
            }
            (None, Some(_)) => {
                let b = right.next().unwrap();
                (b.1, b.2)
            }
            (None, None) => unreachable!(),
        };
        loss += model.loss(q, observed)?;
    }
    ensure(loss.is_finite(), "nonfinite prepared local score")?;
    Ok(loss)
}

fn prepared_pair_totals(
    first: &PreparedProfile,
    second: &PreparedProfile,
) -> io::Result<(u64, u64)> {
    let mut predicted = 0u64;
    let mut observed = 0u64;
    let mut left = first.iter().peekable();
    let mut right = second.iter().peekable();
    while left.peek().is_some() || right.peek().is_some() {
        let (q, count) = match (left.peek(), right.peek()) {
            (Some(a), Some(b)) if a.0 == b.0 => {
                let a = left.next().unwrap();
                let b = right.next().unwrap();
                ensure(a.2 == b.2, "inconsistent cached observation")?;
                (
                    a.1.checked_add(b.1)
                        .ok_or_else(|| invalid("diploid profile overflow"))?,
                    a.2,
                )
            }
            (Some(a), Some(b)) if a.0 < b.0 => {
                let a = left.next().unwrap();
                (a.1, a.2)
            }
            (Some(_), Some(_)) => {
                let b = right.next().unwrap();
                (b.1, b.2)
            }
            (Some(_), None) => {
                let a = left.next().unwrap();
                (a.1, a.2)
            }
            (None, Some(_)) => {
                let b = right.next().unwrap();
                (b.1, b.2)
            }
            (None, None) => unreachable!(),
        };
        predicted = predicted
            .checked_add(q)
            .ok_or_else(|| invalid("owned prediction overflow"))?;
        observed = observed
            .checked_add(count)
            .ok_or_else(|| invalid("owned observation overflow"))?;
    }
    Ok((observed, predicted))
}

fn score_local_pair(
    first: &Profile,
    second: &Profile,
    wanted: LocalClass,
    incidence: &IncidenceTable,
    counts: &mut CountCache<'_>,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = 0.0;
    let mut keys = first.keys().peekable();
    let mut other = second.keys().peekable();
    while keys.peek().is_some() || other.peek().is_some() {
        let key = match (keys.peek(), other.peek()) {
            (Some(a), Some(b)) if a == b => {
                other.next();
                keys.next().unwrap()
            }
            (Some(a), Some(b)) if a < b => keys.next().unwrap(),
            (Some(_), Some(_)) => other.next().unwrap(),
            (Some(_), None) => keys.next().unwrap(),
            (None, Some(_)) => other.next().unwrap(),
            (None, None) => unreachable!(),
        };
        if incidence.classify(key) != wanted {
            continue;
        }
        let q = first
            .get(key)
            .copied()
            .unwrap_or(0)
            .checked_add(second.get(key).copied().unwrap_or(0))
            .ok_or_else(|| invalid("diploid profile overflow"))?;
        loss += model.loss(q, counts.get(key)?)?;
    }
    ensure(loss.is_finite(), "nonfinite local score")?;
    Ok(loss)
}

fn merged_profile(parts: &[&Profile]) -> io::Result<Profile> {
    let mut merged = Profile::new();
    for part in parts {
        for (feature, &value) in *part {
            let total = merged.entry(feature.clone()).or_default();
            *total = total
                .checked_add(value)
                .ok_or_else(|| invalid("cumulative profile overflow"))?;
        }
    }
    Ok(merged)
}

fn score_complete_profile(
    profile: &Profile,
    counts: &mut CountCache<'_>,
    model: &ScoreModel,
) -> io::Result<f64> {
    profile.iter().try_fold(0.0, |loss, (feature, &predicted)| {
        Ok(loss + model.loss(predicted, counts.get(feature)?)?)
    })
}

/// Persistent fid -> count map (HAMT: 32-way branching, 5-bit chunks of
/// the u32 fid, max depth 7 — distinct fids can never collide past bit 31).
/// Structural sharing: `apply` copies only the root-to-leaf paths of the
/// delta's features, so a child state's map shares every other node with
/// its parent. This replaces the per-state flat `FxHashMap` clones that
/// were both the layer-commit's O(profile) cost and the full-component RSS
/// breach: extend and commit are O(delta), lookups are O(depth <= 7), and
/// memory is O(delta) per state instead of O(profile).
#[derive(Clone, Debug, Default)]
struct PersistentCountMap {
    root: Option<Arc<CountNode>>,
}

#[derive(Debug, Clone)]
struct CountNode {
    bitmap: u32,
    slots: Vec<CountSlot>,
}

#[derive(Debug, Clone)]
enum CountSlot {
    Leaf(u32, u64),
    Child(Arc<CountNode>),
}

const COUNT_CHUNK_BITS: usize = 5;
const COUNT_FANOUT: usize = 32;

fn count_chunk(fid: u32, depth: usize) -> u32 {
    (fid >> (COUNT_CHUNK_BITS * depth)) & ((1 << COUNT_CHUNK_BITS) - 1)
}

impl PersistentCountMap {
    fn get(&self, fid: u32) -> u64 {
        let mut node = self.root.as_deref();
        let mut depth = 0usize;
        while let Some(current) = node {
            let chunk = count_chunk(fid, depth);
            if current.bitmap & (1 << chunk) == 0 {
                return 0;
            }
            let index = (current.bitmap & ((1 << chunk) - 1)).count_ones() as usize;
            match &current.slots[index] {
                CountSlot::Leaf(key, count) => return if *key == fid { *count } else { 0 },
                CountSlot::Child(child) => {
                    node = Some(child);
                    depth += 1;
                }
            }
        }
        0
    }

    /// Bulk-build from a fid-sorted, fid-distinct entry list (a distribution
    /// sort by chunk; entries within a bucket keep no order requirement).
    fn from_sorted(entries: &[(u32, u64)]) -> Self {
        Self {
            root: (!entries.is_empty()).then(|| Self::build_node(entries, 0)),
        }
    }

    fn build_node(entries: &[(u32, u64)], depth: usize) -> Arc<CountNode> {
        if entries.len() == 1 {
            return Arc::new(CountNode {
                bitmap: 1 << count_chunk(entries[0].0, depth),
                slots: vec![CountSlot::Leaf(entries[0].0, entries[0].1)],
            });
        }
        let mut buckets = vec![Vec::new(); COUNT_FANOUT];
        for &(fid, count) in entries {
            buckets[count_chunk(fid, depth) as usize].push((fid, count));
        }
        let mut bitmap = 0u32;
        let mut slots = Vec::new();
        for (chunk, bucket) in buckets.into_iter().enumerate() {
            if bucket.is_empty() {
                continue;
            }
            bitmap |= 1 << chunk;
            if bucket.len() == 1 {
                slots.push(CountSlot::Leaf(bucket[0].0, bucket[0].1));
            } else {
                slots.push(CountSlot::Child(Self::build_node(&bucket, depth + 1)));
            }
        }
        Arc::new(CountNode { bitmap, slots })
    }

    /// Apply a fid-sorted delta (fid, added) with structural sharing: every
    /// touched root-to-leaf path is copied once, all other nodes are shared
    /// with the parent map. O(delta x depth).
    fn apply_delta(&self, delta: &[(u32, u64)]) -> io::Result<Self> {
        Ok(Self {
            root: Self::apply_root(self.root.as_ref(), delta, 0)?,
        })
    }

    fn apply_root(
        node: Option<&Arc<CountNode>>,
        delta: &[(u32, u64)],
        depth: usize,
    ) -> io::Result<Option<Arc<CountNode>>> {
        if delta.is_empty() {
            // Untouched subtree: share it unchanged (Arc clone).
            return Ok(node.cloned());
        }
        let Some(current) = node else {
            // No node yet: the delta's values ARE the cumulative counts.
            return Ok(Some(Self::build_node(delta, depth)));
        };
        let mut child_deltas = vec![Vec::new(); COUNT_FANOUT];
        for &(fid, added) in delta {
            let chunk = count_chunk(fid, depth) as usize;
            if child_deltas[chunk]
                .last()
                .is_some_and(|(last, _)| *last == fid)
            {
                return Err(invalid("delta contains a duplicate feature"));
            }
            child_deltas[chunk].push((fid, added));
        }
        let mut bitmap = current.bitmap;
        let mut slots = current.slots.clone();
        for chunk in 0..COUNT_FANOUT {
            let child_delta = &child_deltas[chunk];
            if child_delta.is_empty() {
                continue;
            }
            let bit = 1 << chunk;
            // Index from the UPDATING bitmap: a fresh insert at a lower
            // chunk shifts this chunk's slot, and the updated bitmap already
            // counts every inserted bit below `chunk`.
            let index = (bitmap & (bit - 1)).count_ones() as usize;
            bitmap |= bit;
            if current.bitmap & bit == 0 {
                // Fresh slot: no parent entry exists under this chunk, so
                // the delta's values are the cumulative counts.
                slots.insert(index, CountSlot::Child(Self::build_node(child_delta, depth + 1)));
                continue;
            }
            match &slots[index] {
                CountSlot::Leaf(key, count) => {
                    // One existing entry plus (possibly several) new ones:
                    // fold the leaf into the sub-delta and rebuild the slot.
                    let mut combined = vec![(*key, *count)];
                    for &(fid, added) in child_delta {
                        if let Some(entry) = combined.iter_mut().find(|(f, _)| *f == fid) {
                            entry.1 = entry
                                .1
                                .checked_add(added)
                                .ok_or_else(|| invalid("cumulative profile overflow"))?;
                        } else {
                            combined.push((fid, added));
                        }
                    }
                    if combined.len() == 1 {
                        slots[index] = CountSlot::Leaf(combined[0].0, combined[0].1);
                    } else {
                        slots[index] =
                            CountSlot::Child(Self::build_node(&combined, depth + 1));
                    }
                }
                CountSlot::Child(child) => {
                    let updated = Self::apply_root(Some(child), child_delta, depth + 1)?
                        .expect("apply never removes a node");
                    slots[index] = CountSlot::Child(updated);
                }
            }
        }
        Ok(Some(Arc::new(CountNode { bitmap, slots })))
    }

    /// Iterate every (fid, count) entry (depth-first; order unspecified).
    fn iter(&self) -> CountMapIter<'_> {
        CountMapIter {
            stack: self
                .root
                .as_ref()
                .map(|root| vec![(root.as_ref(), 0usize)])
                .unwrap_or_default(),
            pending: Vec::new(),
        }
    }

    fn len(&self) -> usize {
        fn walk(node: &CountNode) -> usize {
            node.slots
                .iter()
                .map(|slot| match slot {
                    CountSlot::Leaf(_, _) => 1,
                    CountSlot::Child(child) => walk(child),
                })
                .sum()
        }
        self.root.as_deref().map(walk).unwrap_or(0)
    }
}

struct CountMapIter<'a> {
    stack: Vec<(&'a CountNode, usize)>,
    pending: Vec<(u32, u64)>,
}

impl<'a> Iterator for CountMapIter<'a> {
    type Item = (u32, u64);
    fn next(&mut self) -> Option<(u32, u64)> {
        loop {
            if let Some(entry) = self.pending.pop() {
                return Some(entry);
            }
            let (node, slot) = self.stack.last_mut()?;
            if *slot >= node.slots.len() {
                self.stack.pop();
                continue;
            }
            let index = *slot;
            *slot += 1;
            match &node.slots[index] {
                CountSlot::Leaf(fid, count) => self.pending.push((*fid, *count)),
                CountSlot::Child(child) => self.stack.push((child, 0)),
            }
        }
    }
}

#[derive(Debug)]
struct PersistentProfile {
    id: u64,
    parent: Option<Arc<PersistentProfile>>,
    /// The transition's delta as its (up to four) component profiles, kept
    /// fid-sorted and UNMERGED. Merged delta BTreeMaps materialized per class
    /// pair peaked at tens of GB (the RSS guard breach at a wide-reachable
    /// state); the merge is streamed wherever the delta is consumed.
    components: [Arc<ConvertedProfile>; 4],
    fingerprint: u64,
    entries: usize,
    score: f64,
}

impl PersistentProfile {

    /// Stream this chain's cumulative (fid, count) pairs in fid order
    /// WITHOUT materializing a map: a k-way merge over every component of
    /// every ancestor node, allocation-free. Exact interner dedup equality
    /// compares two such streams element-wise, so a fingerprint-bucket
    /// collision no longer pays a full-profile materialization (the O(profile)
    /// spike inside the extend hot path).
    fn cumulative_pairs(&self) -> CumulativePairs<'_> {
        let mut cursors = Vec::new();
        let mut heap = std::collections::BinaryHeap::new();
        let mut node = Some(self);
        while let Some(profile) = node {
            for component in &profile.components {
                if let Some(&(fid, _)) = component.first() {
                    cursors.push((component.as_slice(), 0usize));
                    heap.push(std::cmp::Reverse((fid, cursors.len() - 1)));
                }
            }
            node = profile.parent.as_deref();
        }
        CumulativePairs { cursors, heap }
    }

    /// Exact equality of two chains' cumulative profiles, streamed with early
    /// exit on the first divergence. u64 addition is order-independent, so
    /// the streamed sums equal the materialized map's entries bit for bit.
    fn cumulative_equal(&self, other: &PersistentProfile) -> bool {
        let mut left = self.cumulative_pairs();
        let mut right = other.cumulative_pairs();
        loop {
            match (left.next(), right.next()) {
                (Some(a), Some(b)) if a == b => continue,
                (None, None) => return true,
                _ => return false,
            }
        }
    }

}

/// K-way merge iterator over a PersistentProfile chain's components,
/// yielding (fid, summed count) in fid order.
struct CumulativePairs<'a> {
    cursors: Vec<(&'a [(u32, u64)], usize)>,
    heap: std::collections::BinaryHeap<std::cmp::Reverse<(u32, usize)>>,
}

impl<'a> Iterator for CumulativePairs<'a> {
    type Item = (u32, u64);
    fn next(&mut self) -> Option<(u32, u64)> {
        let std::cmp::Reverse((fid, _)) = *self.heap.peek()?;
        let mut total = 0u64;
        // Pop and sum every cursor parked at this fid, advancing each; a
        // cursor re-pushes its successor entry. Saturating addition can never
        // engage here: counts are checked for overflow at every extend.
        while let Some(std::cmp::Reverse((top, index))) = self.heap.peek().copied() {
            if top != fid {
                break;
            }
            self.heap.pop();
            let (stream, position) = self.cursors[index];
            total = total.saturating_add(stream[position].1);
            let next = position + 1;
            self.cursors[index] = (stream, next);
            if let Some(&(next_fid, _)) = stream.get(next) {
                self.heap.push(std::cmp::Reverse((next_fid, index)));
            }
        }
        Some((fid, total))
    }
}

// Temporary hot-path probe: relaxed atomic counters, printed per locus under
// IMPG_DB. Used only to calibrate the per-extend cost breakdown; zero cost
// when the probe is off (one boolean check per extend).
/// DP sub-stage wall accumulators (nanoseconds), summed across every
/// `bounded_conflict_chain` pass (seed ladder, main pass, tie rotations) and
/// the per-locus evidence stages of `exact_chain_streaming`. Probe-only: no
/// search semantics depends on these counters; the smoke reports them under
/// stage_probe to itemize the DP stage.
static DP_STAGE_NANOS: [std::sync::atomic::AtomicU64; 13] = [
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
];
static DP_STAGE_PASSES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

pub struct DpStageProbe;

impl DpStageProbe {
    fn add(slot: usize, started: Instant) {
        DP_STAGE_NANOS[slot].fetch_add(
            started.elapsed().as_nanos().min(u64::MAX as u128) as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }

    fn note_pass() {
        DP_STAGE_PASSES.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }

    /// (per-slot seconds, pass count). Slot order: successor registries,
    /// feature universe, initial layer, seed ladder, per-locus weight tables,
    /// suffix max-add tables, expansion loop, layer commit, suffix rebase,
    /// terminal collection, per-locus class evidence, streamed profile load,
    /// audit_chain re-run at chain entry.
    pub fn snapshot() -> (Vec<f64>, u64) {
        let seconds = DP_STAGE_NANOS
            .iter()
            .map(|counter| {
                counter.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9
            })
            .collect();
        (
            seconds,
            DP_STAGE_PASSES.load(std::sync::atomic::Ordering::Relaxed),
        )
    }
}

struct DpProbe;

impl DpProbe {
    fn on() -> bool {
        static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
        *ON.get_or_init(|| std::env::var("IMPG_DB").is_ok())
    }

    fn counter() -> &'static [std::sync::atomic::AtomicU64; 9] {
        use std::sync::atomic::AtomicU64;
        static V: [AtomicU64; 9] = [
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
            AtomicU64::new(0),
        ];
        &V
    }

    fn snapshot() -> [u64; 9] {
        use std::sync::atomic::Ordering::Relaxed;
        let counters = Self::counter();
        [
            counters[0].load(Relaxed),
            counters[1].load(Relaxed),
            counters[2].load(Relaxed),
            counters[3].load(Relaxed),
            counters[4].load(Relaxed),
            counters[5].load(Relaxed),
            counters[6].load(Relaxed),
            counters[7].load(Relaxed),
            counters[8].load(Relaxed),
        ]
    }

    fn print_delta(stage: &str, before: &[u64; 9]) {
        let after = Self::snapshot();
        let d = |i: usize| after[i] - before[i];
        eprintln!(
            "PROBE {stage} extends {} features {} loop_ms {:.1} depth_avg {:.1} dedup_ms {:.1} mat {} buckets {} tight_ms {:.1} tight_calls {}",
            d(0),
            d(1),
            d(2) as f64 / 1e6,
            d(3) as f64 / d(0).max(1) as f64,
            d(4) as f64 / 1e6,
            d(5),
            d(6),
            d(7) as f64 / 1e6,
            d(8),
        );
    }
}

/// Feature identity for the DP hot path. The universe is built once per
/// chain from every profile-class and seam-signature feature, sorted in
/// FeatureKey order, so `fid` order EQUALS FeatureKey order and every
/// fid-ordered merge iterates features exactly in the merged BTreeMap order
/// the previous interner used — the telescoped score summation order (and
/// therefore every loss bit) is preserved.
struct FeatureUniverse {
    /// fid -> feature, sorted by feature.
    keys: Vec<FeatureKey>,
    /// feature -> fid.
    ids: FxHashMap<FeatureKey, u32>,
    /// Per-fid FNV base (feature tokens only); a per-(feature, value)
    /// fingerprint is one extra FNV step from this base, bit-identical to
    /// profile_value_fingerprint.
    base_hash: Vec<u64>,
    /// Per-fid observed sample count, resolved once at construction.
    observed: Vec<u64>,
}

impl FeatureUniverse {
    fn build(
        profile_classes: &[Vec<Profile>],
        locus_registries: &[(
            Vec<SuccessorClass>,
            BTreeMap<usize, BTreeMap<usize, Vec<usize>>>,
        )],
        counts: &mut CountCache<'_>,
    ) -> io::Result<Self> {
        let mut keys = Vec::new();
        for classes in profile_classes {
            for profile in classes {
                keys.extend(profile.keys().cloned());
            }
        }
        // The seam signatures enter the DP only through the per-locus
        // registries, whose distinct (profile class, signature) classes
        // already deduplicate the per-seam profiles; iterating the raw seams
        // (1.7M on the slice) cloned every feature key once per seam and
        // dominated the DP stage wall clock.
        for (classes, _) in locus_registries.iter() {
            for (_, signature) in classes.iter() {
                for (feature, _) in signature.iter() {
                    keys.push(feature.clone());
                }
            }
        }
        keys.sort_unstable();
        keys.dedup();
        let total = keys.len();
        let mut ids = FxHashMap::default();
        ids.reserve(total);
        let mut base_hash = Vec::with_capacity(total);
        let mut observed = Vec::with_capacity(total);
        for (fid, key) in keys.iter().enumerate() {
            ids.insert(key.clone(), fid as u32);
            let mut hash = 0xcbf29ce484222325u64;
            hash = fnv_update(hash, &(key.len() as u64).to_le_bytes());
            for token in key {
                hash = fnv_update(hash, &token.to_le_bytes());
            }
            base_hash.push(hash);
            observed.push(counts.get(key)?);
        }
        Ok(Self {
            keys,
            ids,
            base_hash,
            observed,
        })
    }

    fn id(&self, feature: &FeatureKey) -> io::Result<u32> {
        self.ids
            .get(feature)
            .copied()
            .ok_or_else(|| invalid("dp feature outside the chain feature universe"))
    }

    /// Fingerprint contribution of (fid, value): one FNV step from the
    /// feature's base hash — bit-identical to hashing the feature tokens and
    /// the value sequentially.
    fn fingerprint_value(&self, fid: u32, value: u64) -> u64 {
        fnv_update(self.base_hash[fid as usize], &value.to_le_bytes())
    }
}

/// A profile converted to fid-sorted (fid, count) pairs. Conversion preserves
/// the source BTreeMap's FeatureKey order, which equals fid order.
type ConvertedProfile = Vec<(u32, u64)>;

fn convert_profile(
    profile: &Profile,
    universe: &FeatureUniverse,
) -> io::Result<Arc<ConvertedProfile>> {
    let mut converted = Vec::with_capacity(profile.len());
    for (feature, &count) in profile {
        converted.push((universe.id(feature)?, count));
    }
    Ok(Arc::new(converted))
}

/// A seam signature is already a (feature, count) vector in feature order.
fn convert_signature(
    signature: &[(FeatureKey, u64)],
    universe: &FeatureUniverse,
) -> io::Result<Arc<ConvertedProfile>> {
    let mut converted = Vec::with_capacity(signature.len());
    for (feature, count) in signature {
        converted.push((universe.id(feature)?, *count));
    }
    Ok(Arc::new(converted))
}

/// Empty shared component (the initial delta's absent seam slots).
fn empty_component() -> Arc<ConvertedProfile> {
    Arc::new(Vec::new())
}

/// Merge up to four fid-sorted components into one fid-ordered stream of
/// (fid, total count). This replaces materializing the merged delta profile:
/// the merge is streamed, allocation-free, and iterates in exactly the order
/// the merged BTreeMap would.
struct ComponentMerge<'a> {
    streams: [&'a [(u32, u64)]; 4],
    heads: [usize; 4],
}

impl<'a> ComponentMerge<'a> {
    fn new(components: &'a [Arc<ConvertedProfile>; 4]) -> Self {
        Self {
            streams: [
                components[0].as_slice(),
                components[1].as_slice(),
                components[2].as_slice(),
                components[3].as_slice(),
            ],
            heads: [0; 4],
        }
    }

    fn next(&mut self) -> io::Result<Option<(u32, u64)>> {
        let mut min_fid = u32::MAX;
        for stream_index in 0..4 {
            let head = self.heads[stream_index];
            if let Some(&(fid, _)) = self.streams[stream_index].get(head) {
                min_fid = min_fid.min(fid);
            }
        }
        if min_fid == u32::MAX {
            return Ok(None);
        }
        let mut total = 0u64;
        for stream_index in 0..4 {
            let head = self.heads[stream_index];
            if let Some(&(fid, count)) = self.streams[stream_index].get(head) {
                if fid == min_fid {
                    total = total
                        .checked_add(count)
                        .ok_or_else(|| invalid("cumulative profile overflow"))?;
                    self.heads[stream_index] = head + 1;
                }
            }
        }
        Ok(Some((min_fid, total)))
    }

    /// Like `next`, but also reports the single contributing stream and its
    /// head position when exactly one component contributed the feature.
    /// The extend hot path uses the tag to add the contributing component's
    /// PRECOMPUTED fresh loss term (terms[s][h] == model.loss(count,
    /// observed[fid]) — the identical value the inline computation produced)
    /// instead of re-evaluating the loss per feature. Only multi-component
    /// features fall back to an inline loss evaluation of the merged total.
    fn next_tagged(&mut self) -> io::Result<Option<(u32, u64, Option<(usize, usize)>)>> {
        let mut min_fid = u32::MAX;
        for stream_index in 0..4 {
            let head = self.heads[stream_index];
            if let Some(&(fid, _)) = self.streams[stream_index].get(head) {
                min_fid = min_fid.min(fid);
            }
        }
        if min_fid == u32::MAX {
            return Ok(None);
        }
        let mut total = 0u64;
        let mut contributors = 0usize;
        let mut single: Option<(usize, usize)> = None;
        for stream_index in 0..4 {
            let head = self.heads[stream_index];
            if let Some(&(fid, count)) = self.streams[stream_index].get(head) {
                if fid == min_fid {
                    total = total
                        .checked_add(count)
                        .ok_or_else(|| invalid("cumulative profile overflow"))?;
                    self.heads[stream_index] = head + 1;
                    contributors += 1;
                    if contributors == 1 {
                        single = Some((stream_index, head));
                    } else {
                        single = None;
                    }
                }
            }
        }
        Ok(Some((min_fid, total, single)))
    }
}

/// Remaining-loss potential of one feature at cumulative count `q` when the
/// future can add at most `max_add`: c(q) = min_{x in [q, q+max_add]} loss(x)
/// - loss(q). The loss is convex in q, so the minimum sits at an endpoint or
/// at the clamped continuous optimum; every candidate is evaluated (dropping
/// the clamped-to-endpoint points made an earlier version of this bound
/// inadmissible for under-covered features). c(q) = 0 exactly when the state
/// has already saturated the feature.
fn feature_potential(
    q: u64,
    max_add: u64,
    observed: u64,
    histogram: u64,
    depth: f64,
    denominator: f64,
    background: f64,
) -> f64 {
    let loss = |count: u64| -> f64 {
        let signal = (count * histogram) as f64 * depth / denominator;
        signal - observed as f64 * (signal / background).ln_1p()
    };
    let hi = q.saturating_add(max_add);
    let mut best = loss(q);
    if hi != q {
        best = best.min(loss(hi));
    }
    let per_count = histogram as f64 * depth / denominator;
    let optimum = (observed as f64 - background) / per_count;
    if optimum.is_finite() && optimum > 0.0 {
        let floored = (optimum as u64).max(q).min(hi);
        if floored != q {
            best = best.min(loss(floored));
        }
        let ceiled = (optimum.ceil() as u64).max(q).min(hi);
        if ceiled != floored {
            best = best.min(loss(ceiled));
        }
    }
    best - loss(q)
}

fn profile_value_fingerprint(feature: &FeatureKey, value: u64) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    hash = fnv_update(hash, &(feature.len() as u64).to_le_bytes());
    for token in feature {
        hash = fnv_update(hash, &token.to_le_bytes());
    }
    fnv_update(hash, &value.to_le_bytes())
}

/// Per-(locus, class-pair) fresh-feature constants, computed once per pair
/// per run and shared by every DP pass: for a state that has touched NONE of
/// the pair's delta features, every merged delta feature f contributes
/// exactly loss(d_f) to the score (loss(0) == 0.0 exactly), fp(fid, d_f) to
/// the commutative XOR fingerprint, and +1 entry — all independent of the
/// state. k_hash/fresh_entries regroupings are bit-identical to the
/// per-feature fold (XOR and integer addition are associative); the k_score
/// fold is algebraically identical to the per-feature telescoping but
/// regroups the f64 additions, so only the dev telescoped extend consumes
/// it (production keeps the per-feature fold and stays bit-identical).
#[derive(Clone, Copy)]
struct PairFreshConstants {
    k_score: f64,
    k_hash: u64,
    fresh_entries: usize,
}

#[derive(Default)]
struct PersistentProfileInterner {
    next_id: u64,
    extensions: HashMap<(u64, [usize; 4]), std::sync::Weak<PersistentProfile>>,
    buckets: HashMap<(u64, usize), Vec<std::sync::Weak<PersistentProfile>>>,
    /// Lazily computed per-component fresh loss terms keyed by component
    /// Arc address: terms[i] == model.loss(component[i].1, observed[fid]) —
    /// the identical f64 the per-feature loop evaluated inline, so adding
    /// the indexed term preserves the telescoped score's bits exactly.
    component_terms: FxHashMap<usize, Arc<Vec<f64>>>,
    /// IMPG_DEV_TELESCOPED_EXTEND (dev-only): extends consume the per-pair
    /// constants plus O(overlap) corrections instead of the per-feature
    /// fold. Algebraically identical; NOT bit-identical (the f64 additions
    /// regroup), so it must never be on in production or in any run whose
    /// accounting is quoted.
    telescoped: bool,
}

impl PersistentProfileInterner {
    fn clear_search_indexes(&mut self) {
        self.extensions.clear();
        self.buckets.clear();
        self.component_terms.clear();
    }

    fn component_terms(
        &mut self,
        component: &Arc<ConvertedProfile>,
        universe: &FeatureUniverse,
        model: &ScoreModel,
    ) -> io::Result<Arc<Vec<f64>>> {
        let key = Arc::as_ptr(component) as usize;
        if let Some(terms) = self.component_terms.get(&key) {
            return Ok(Arc::clone(terms));
        }
        let mut terms = Vec::with_capacity(component.len());
        for &(fid, count) in component.iter() {
            terms.push(model.loss(count, universe.observed[fid as usize])?);
        }
        let terms = Arc::new(terms);
        self.component_terms.insert(key, Arc::clone(&terms));
        Ok(terms)
    }

    #[allow(clippy::too_many_arguments)]
    fn extend(
        &mut self,
        parent: Option<&Arc<PersistentProfile>>,
        // The state's touched counts restricted to this boundary's delta
        // feature universe: {fid: q > 0} ∩ (every reachable pair's delta
        // features). A delta feature absent here has q == 0 and is fresh by
        // construction. None exactly when parent is None (initial layer).
        small: Option<&FxHashMap<u32, u64>>,
        components: [Arc<ConvertedProfile>; 4],
        base_score: f64,
        constants: PairFreshConstants,
        universe: &FeatureUniverse,
        model: &ScoreModel,
    ) -> io::Result<(f64, Arc<PersistentProfile>)> {
        let extension_key = (
            parent.map_or(u64::MAX, |profile| profile.id),
            [
                Arc::as_ptr(&components[0]) as usize,
                Arc::as_ptr(&components[1]) as usize,
                Arc::as_ptr(&components[2]) as usize,
                Arc::as_ptr(&components[3]) as usize,
            ],
        );
        if let Some(existing) = self
            .extensions
            .get(&extension_key)
            .and_then(|p| p.upgrade())
        {
            return Ok((existing.score, existing));
        }
        let probe = DpProbe::on();
        let probe_start = std::time::Instant::now();
        let mut score = base_score;
        let mut fingerprint = parent.map_or(0, |profile| profile.fingerprint);
        // The pair's fresh constants carry every fresh feature's entry gain
        // and fingerprint contribution; overlap features correct both back.
        let mut entries = parent.map_or(0, |profile| profile.entries) + constants.fresh_entries;
        let mut depth = 0usize;
        if let Some(parent) = parent {
            let mut node = Some(parent.as_ref());
            while let Some(profile) = node {
                depth += 1;
                node = profile.parent.as_deref();
            }
        }
        if self.telescoped {
            // Dev-only telescoped extend: child = base + k_score plus an
            // O(|touched ∩ delta|) correction per feature the state already
            // touched (loss/fingerprint/entries each). Algebraically
            // identical to the fold below; regroups the f64 additions, so
            // dev runs only.
            score += constants.k_score;
            if let Some(small) = small {
                for (&fid, &q) in small.iter() {
                    let mut total = 0u64;
                    let mut present = false;
                    for component in components.iter() {
                        if let Ok(index) = component.binary_search_by_key(&fid, |&(f, _)| f) {
                            present = true;
                            total = total
                                .checked_add(component[index].1)
                                .ok_or_else(|| invalid("cumulative profile overflow"))?;
                        }
                    }
                    if !present {
                        continue;
                    }
                    let observed = universe.observed[fid as usize];
                    let next = q
                        .checked_add(total)
                        .ok_or_else(|| invalid("cumulative profile overflow"))?;
                    score += model.loss(next, observed)?
                        - model.loss(q, observed)?
                        - model.loss(total, observed)?;
                    fingerprint ^= universe.fingerprint_value(fid, q)
                        ^ universe.fingerprint_value(fid, next)
                        ^ universe.fingerprint_value(fid, total);
                    entries -= 1;
                }
            }
        } else {
            // Production fold: every feature's score increment is added in
            // merged-fid order exactly as the historical per-feature loop
            // did, with bit-identical values (the indexed term IS the loss
            // the inline evaluation produced; loss(0) == 0.0 exactly), so
            // the telescoped score keeps its bits. Fresh features skip the
            // inline loss evaluation and the observed[] load entirely;
            // their entry gain and fingerprint contribution come from the
            // pair constants. Overlap features (the state already touched
            // them) keep the full inline correction.
            let terms = [
                self.component_terms(&components[0], universe, model)?,
                self.component_terms(&components[1], universe, model)?,
                self.component_terms(&components[2], universe, model)?,
                self.component_terms(&components[3], universe, model)?,
            ];
            let empty_small_map;
            let small = match small {
                Some(map) => map,
                None => {
                    empty_small_map = FxHashMap::default();
                    &empty_small_map
                }
            };
            let mut merge = ComponentMerge::new(&components);
            while let Some((fid, added, single)) = merge.next_tagged()? {
                let previous = small.get(&fid).copied().unwrap_or(0);
                let next = previous
                    .checked_add(added)
                    .ok_or_else(|| invalid("cumulative profile overflow"))?;
                if previous == 0 {
                    score += match single {
                        Some((slot, head)) => terms[slot][head],
                        None => model.loss(added, universe.observed[fid as usize])?,
                    };
                } else {
                    let observed = universe.observed[fid as usize];
                    score += model.loss(next, observed)? - model.loss(previous, observed)?;
                    // The pair constant already folded fp(fid, added) for
                    // this feature (the fresh assumption); XOR it back out
                    // and fold the true contribution fp(previous) ^
                    // fp(next).
                    fingerprint ^= universe.fingerprint_value(fid, previous)
                        ^ universe.fingerprint_value(fid, next)
                        ^ universe.fingerprint_value(fid, added);
                    entries -= 1;
                }
            }
        }
        fingerprint ^= constants.k_hash;
        if probe {
            use std::sync::atomic::Ordering::Relaxed;
            let counters = DpProbe::counter();
            let features: u64 = components
                .iter()
                .map(|component| component.len() as u64)
                .sum();
            counters[0].fetch_add(1, Relaxed);
            counters[1].fetch_add(features, Relaxed);
            counters[2].fetch_add(probe_start.elapsed().as_nanos() as u64, Relaxed);
            counters[3].fetch_add(depth as u64, Relaxed);
        }
        let candidate = Arc::new(PersistentProfile {
            id: self.next_id,
            parent: parent.cloned(),
            components,
            fingerprint,
            entries,
            score,
        });
        let dedup_start = std::time::Instant::now();
        let bucket_key = (fingerprint, entries);
        if let Some(bucket) = self.buckets.get_mut(&bucket_key) {
            bucket.retain(|profile| profile.strong_count() > 0);
            let live: Vec<_> = bucket.iter().filter_map(|w| w.upgrade()).collect();
            if !live.is_empty() {
                if probe {
                    use std::sync::atomic::Ordering::Relaxed;
                    DpProbe::counter()[6].fetch_add(1, Relaxed);
                }
                for existing in live {
                    if probe {
                        use std::sync::atomic::Ordering::Relaxed;
                        DpProbe::counter()[5].fetch_add(1, Relaxed);
                    }
                    if candidate.cumulative_equal(&existing) {
                        self.extensions
                            .insert(extension_key, Arc::downgrade(&existing));
                        if probe {
                            use std::sync::atomic::Ordering::Relaxed;
                            DpProbe::counter()[4]
                                .fetch_add(dedup_start.elapsed().as_nanos() as u64, Relaxed);
                        }
                        return Ok((existing.score, existing));
                    }
                }
            }
            if probe {
                use std::sync::atomic::Ordering::Relaxed;
                DpProbe::counter()[4].fetch_add(dedup_start.elapsed().as_nanos() as u64, Relaxed);
            }
        }
        self.next_id = self.next_id.saturating_add(1);
        self.buckets
            .entry(bucket_key)
            .or_default()
            .push(Arc::downgrade(&candidate));
        self.extensions
            .insert(extension_key, Arc::downgrade(&candidate));
        Ok((score, candidate))
    }
}

fn exit_memberships(
    partitions: &[Vec<GenomeAllele>],
    seams: &[Vec<GenomeSeam>],
) -> (Vec<Vec<usize>>, Vec<usize>) {
    let mut memberships = vec![Vec::new(); partitions.len()];
    let mut counts = vec![0; partitions.len()];
    for locus in (0..partitions.len()).rev() {
        let mut signatures: BTreeMap<Vec<(usize, Vec<(FeatureKey, u64)>)>, usize> = BTreeMap::new();
        let mut membership = vec![0; partitions[locus].len()];
        for allele in 0..partitions[locus].len() {
            let mut exits = if locus + 1 < partitions.len() {
                seams[locus]
                    .iter()
                    .filter(|seam| seam.left == allele)
                    .map(|seam| {
                        (
                            memberships[locus + 1][seam.right],
                            profile_signature(&seam.profile),
                        )
                    })
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            };
            exits.sort();
            exits.dedup();
            let next = signatures.len();
            membership[allele] = *signatures.entry(exits).or_insert(next);
        }
        counts[locus] = signatures.len();
        memberships[locus] = membership;
    }
    (memberships, counts)
}

pub fn audit_chain(
    partitions: &[Vec<GenomeAllele>],
    seams: &[Vec<GenomeSeam>],
) -> io::Result<GenomeAudit> {
    ensure(
        !partitions.is_empty() && seams.len() + 1 == partitions.len(),
        "invalid genome partition chain",
    )?;
    let (_, exit_classes) = exit_memberships(partitions, seams);
    let mut by_source: BTreeMap<usize, Vec<(usize, &SourceRange)>> = BTreeMap::new();
    for (locus, alleles) in partitions.iter().enumerate() {
        for allele in alleles {
            for segment in &allele.traversal.segments {
                by_source
                    .entry(segment.source)
                    .or_default()
                    .push((locus, segment));
            }
        }
    }
    let mut conflict_window = 0;
    for ranges in by_source.values_mut() {
        ranges.sort_by_key(|(locus, range)| (range.start, range.end, *locus));
        for left in 0..ranges.len() {
            for right in left + 1..ranges.len() {
                if ranges[right].1.start >= ranges[left].1.end {
                    break;
                }
                if ranges[left].1.start < ranges[right].1.end {
                    conflict_window = conflict_window.max(ranges[left].0.abs_diff(ranges[right].0));
                }
            }
        }
    }
    let boundary_state_products = exit_classes
        .windows(2)
        .map(|e| {
            let a = e[0] as u128;
            let b = e[1] as u128;
            a * a * b * b
        })
        .collect();
    Ok(GenomeAudit {
        exit_classes,
        conflict_window,
        boundary_state_products,
    })
}

#[derive(Clone, Debug)]
struct PhysicalNode {
    score: f64,
    choices: [usize; 2],
    predecessors: Vec<usize>,
}

fn cumulative_profile_bytes(profiles: &BTreeMap<usize, Arc<PersistentProfile>>) -> u64 {
    let mut seen_nodes = std::collections::HashSet::new();
    let mut seen_components = std::collections::HashSet::new();
    let mut bytes = 0u64;
    for profile in profiles.values() {
        let mut node = Some(profile.as_ref());
        while let Some(profile) = node {
            if !seen_nodes.insert(profile.id) {
                break;
            }
            bytes = bytes.saturating_add(std::mem::size_of::<PersistentProfile>() as u64);
            for component in &profile.components {
                if seen_components.insert(Arc::as_ptr(component) as usize) {
                    bytes = bytes.saturating_add(
                        component
                            .len()
                            .saturating_mul(24)
                            .saturating_add(64) as u64,
                    );
                }
            }
            node = profile.parent.as_deref();
        }
    }
    bytes
}

fn physical_state_bytes(states: usize, nodes: &[PhysicalNode]) -> u64 {
    let ties: usize = nodes.iter().map(|node| node.predecessors.len()).sum();
    states
        .saturating_mul(128)
        .saturating_add(
            nodes
                .len()
                .saturating_mul(std::mem::size_of::<PhysicalNode>()),
        )
        .saturating_add(ties.saturating_mul(std::mem::size_of::<usize>())) as u64
}

fn backtrack_physical(
    nodes: &[PhysicalNode],
    node: usize,
    reverse: &mut Vec<[usize; 2]>,
    output: &mut Vec<[Vec<usize>; 2]>,
    truncated: &mut bool,
) {
    if output.len() >= MAX_COMPLETE_PAIRS {
        *truncated = true;
        return;
    }
    reverse.push(nodes[node].choices);
    if nodes[node].predecessors.is_empty() {
        let history = reverse.iter().rev().collect::<Vec<_>>();
        output.push([
            history.iter().map(|pair| pair[0]).collect(),
            history.iter().map(|pair| pair[1]).collect(),
        ]);
    } else {
        for &previous in &nodes[node].predecessors {
            backtrack_physical(nodes, previous, reverse, output, truncated);
        }
    }
    reverse.pop();
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct BeamKey {
    history: Vec<[usize; 2]>,
    profile_id: u64,
}

#[derive(Clone, Debug)]
struct BeamCandidate {
    score: f64,
    choices: [usize; 2],
    predecessors: Vec<usize>,
    cumulative_profile: Arc<PersistentProfile>,
}

#[derive(Clone, Copy, Debug)]
struct ScoreKey(f64);

impl PartialEq for ScoreKey {
    fn eq(&self, other: &Self) -> bool {
        self.0.total_cmp(&other.0).is_eq()
    }
}
impl Eq for ScoreKey {}
impl PartialOrd for ScoreKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for ScoreKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

/// Seeded 64-bit fingerprint of a beam key, used only to reorder equal-score
/// tie classes across tie-rotation passes.
fn tie_fingerprint(rotation: u64, key: &BeamKey) -> u64 {
    let mut hash = 0xcbf29ce484222325u64 ^ rotation.wrapping_mul(0x9e3779b97f4a7c15);
    for pair in &key.history {
        for &value in pair {
            hash ^= value as u64;
            hash = hash.wrapping_mul(0x100000001b3);
        }
    }
    hash ^= key.profile_id;
    hash.wrapping_mul(0x100000001b3)
}

/// Order two beam keys under a tie rotation. Rotation 0 is the historical
/// BeamKey order (bit-identical retention); rotation r > 0 orders tied states
/// by a seeded fingerprint so each rotation retains a different — still
/// deterministic — subset of an over-width tie class. Score always dominates;
/// only equal-score ties rotate.
fn cmp_ranked(rotation: u64, left: &BeamKey, right: &BeamKey) -> std::cmp::Ordering {
    if rotation == 0 {
        left.cmp(right)
    } else {
        tie_fingerprint(rotation, left)
            .cmp(&tie_fingerprint(rotation, right))
            .then_with(|| left.cmp(right))
    }
}

/// Ranked beam key: a BeamKey plus the tie rotation its layer runs under.
/// Both sides of any comparison share the rotation, so the ordering is
/// `cmp_ranked` by construction.
#[derive(Clone, Debug)]
struct RankKey {
    rotation: u64,
    key: BeamKey,
}

impl PartialEq for RankKey {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}
impl Eq for RankKey {}
impl PartialOrd for RankKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for RankKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        cmp_ranked(self.rotation, &self.key, &other.key)
    }
}

struct BoundedSuccessors {
    limit: usize,
    rotation: u64,
    candidates: BTreeMap<BeamKey, BeamCandidate>,
    ranking: BTreeSet<(ScoreKey, RankKey)>,
    dropped_ties: u64,
    dropped_non_ties: u64,
}

impl BoundedSuccessors {
    fn new(limit: usize, rotation: u64) -> Self {
        Self {
            limit,
            rotation,
            candidates: BTreeMap::new(),
            ranking: BTreeSet::new(),
            dropped_ties: 0,
            dropped_non_ties: 0,
        }
    }

    fn cutoff(&self) -> Option<f64> {
        self.ranking.last().map(|entry| entry.0 .0)
    }

    fn refresh_cutoff(&mut self, previous: Option<f64>) {
        if let (Some(previous), Some(current)) = (previous, self.cutoff()) {
            if current.total_cmp(&previous).is_lt() {
                self.dropped_non_ties = self.dropped_non_ties.saturating_add(self.dropped_ties);
                self.dropped_ties = 0;
            }
        }
    }

    fn record_drop(&mut self, score: f64) {
        if self
            .cutoff()
            .is_some_and(|cutoff| score.total_cmp(&cutoff).is_eq())
        {
            self.dropped_ties = self.dropped_ties.saturating_add(1);
        } else {
            self.dropped_non_ties = self.dropped_non_ties.saturating_add(1);
        }
    }

    fn record_bulk_drop(&mut self, score: f64, count: u64) {
        if self
            .cutoff()
            .is_some_and(|cutoff| score.total_cmp(&cutoff).is_eq())
        {
            self.dropped_ties = self.dropped_ties.saturating_add(count);
        } else {
            self.dropped_non_ties = self.dropped_non_ties.saturating_add(count);
        }
    }

    fn insert(
        &mut self,
        key: BeamKey,
        score: f64,
        choices: [usize; 2],
        predecessor: usize,
        cumulative_profile: Arc<PersistentProfile>,
        accounting: &mut ChainAccounting,
    ) {
        let previous_cutoff = self.cutoff();
        if let Some(existing) = self.candidates.get_mut(&key) {
            match score.total_cmp(&existing.score) {
                std::cmp::Ordering::Less => {
                    self.ranking
                        .remove(&(ScoreKey(existing.score), RankKey {
                            rotation: self.rotation,
                            key: key.clone(),
                        }));
                    existing.score = score;
                    existing.choices = choices;
                    existing.predecessors.clear();
                    existing.predecessors.push(predecessor);
                    existing.cumulative_profile = cumulative_profile;
                    self.ranking.insert((
                        ScoreKey(score),
                        RankKey {
                            rotation: self.rotation,
                            key,
                        },
                    ));
                    self.refresh_cutoff(previous_cutoff);
                }
                std::cmp::Ordering::Equal => {
                    existing.predecessors.push(predecessor);
                    accounting.tie_backpointers = accounting.tie_backpointers.saturating_add(1);
                }
                std::cmp::Ordering::Greater => {}
            }
            return;
        }
        if self.candidates.len() < self.limit {
            self.ranking.insert((
                ScoreKey(score),
                RankKey {
                    rotation: self.rotation,
                    key: key.clone(),
                },
            ));
            self.candidates.insert(
                key,
                BeamCandidate {
                    score,
                    choices,
                    predecessors: vec![predecessor],
                    cumulative_profile,
                },
            );
            return;
        }
        let worst = self
            .ranking
            .last()
            .cloned()
            .expect("nonempty bounded layer");
        let rank = RankKey {
            rotation: self.rotation,
            key: key.clone(),
        };
        if (ScoreKey(score), &rank) < (worst.0, &worst.1) {
            let previous_cutoff = Some(worst.0 .0);
            let worst_key = worst.1.key.clone();
            self.ranking.remove(&(worst.0, worst.1));
            let evicted = self
                .candidates
                .remove(&worst_key)
                .expect("ranked bounded candidate");
            self.ranking.insert((ScoreKey(score), rank));
            self.candidates.insert(
                key,
                BeamCandidate {
                    score,
                    choices,
                    predecessors: vec![predecessor],
                    cumulative_profile,
                },
            );
            self.refresh_cutoff(previous_cutoff);
            self.record_drop(evicted.score);
        } else {
            self.record_drop(score);
        }
    }
}

fn conflict_copy_eligible(
    partitions: &[Vec<GenomeAllele>],
    locus: usize,
    history: &[[usize; 2]],
    copy: usize,
    allele: usize,
) -> bool {
    if !spans_feasible(&partitions[locus][allele].traversal.segments) {
        return false;
    }
    let first_locus = locus - history.len();
    for (offset, previous) in history.iter().enumerate() {
        let previous_locus = first_locus + offset;
        let left = &partitions[previous_locus][previous[copy]]
            .traversal
            .segments;
        let right = &partitions[locus][allele].traversal.segments;
        if left.iter().any(|left| {
            right.iter().any(|right| {
                left.source == right.source && left.start < right.end && right.start < left.end
            })
        }) {
            return false;
        }
    }
    true
}

fn conflict_eligible(
    partitions: &[Vec<GenomeAllele>],
    locus: usize,
    history: &[[usize; 2]],
    pair: [usize; 2],
) -> bool {
    (0..2).all(|copy| conflict_copy_eligible(partitions, locus, history, copy, pair[copy]))
}

/// Exact per-copy conflict eligibility for every member allele at `locus`
/// against one state's conflict-window history, matching
/// `conflict_copy_eligible` member-by-member. Candidates are narrowed by the
/// per-source span index: an overlap requires the member's span to start
/// below the past segment's end, and spanning axis intervals are disjoint, so
/// only split traversals reaching into the past segment's partition qualify.
fn conflict_eligible_members(
    partitions: &[Vec<GenomeAllele>],
    locus: usize,
    history: &[[usize; 2]],
    copy: usize,
    spans_by_source: &BTreeMap<usize, Vec<(u64, u64, usize)>>,
) -> Vec<bool> {
    let mut eligible = partitions[locus]
        .iter()
        .map(|allele| spans_feasible(&allele.traversal.segments))
        .collect::<Vec<_>>();
    let first_locus = locus - history.len();
    for (offset, previous) in history.iter().enumerate() {
        let previous_locus = first_locus + offset;
        for segment in &partitions[previous_locus][previous[copy]].traversal.segments {
            let Some(spans) = spans_by_source.get(&segment.source) else {
                continue;
            };
            let prefix = spans.partition_point(|&(start, _, _)| start < segment.end);
            for &(_start, end, member) in spans[..prefix].iter() {
                if !eligible[member] || end <= segment.start {
                    continue;
                }
                if partitions[locus][member]
                    .traversal
                    .segments
                    .iter()
                    .any(|right| {
                        segment.source == right.source
                            && segment.start < right.end
                            && right.start < segment.end
                    })
                {
                    eligible[member] = false;
                }
            }
        }
    }
    eligible
}

fn prune_beam(
    candidates: BTreeMap<BeamKey, BeamCandidate>,
    accounting: &mut ChainAccounting,
    beam_width: usize,
    tie_rotation: u64,
) -> Vec<(BeamKey, BeamCandidate)> {
    let mut ranked = candidates.into_iter().collect::<Vec<_>>();
    ranked.sort_by(|left, right| {
        left.1
            .score
            .total_cmp(&right.1.score)
            .then_with(|| cmp_ranked(tie_rotation, &left.0, &right.0))
    });
    if ranked.len() > beam_width {
        let cutoff = ranked[beam_width - 1].1.score;
        for (_, candidate) in &ranked[beam_width..] {
            if candidate.score.total_cmp(&cutoff).is_eq() {
                accounting.dropped_ties = accounting.dropped_ties.saturating_add(1);
            } else {
                accounting.dropped_non_ties = accounting.dropped_non_ties.saturating_add(1);
            }
        }
        ranked.truncate(beam_width);
    }
    ranked
}

fn completion_viability(
    partitions: &[Vec<GenomeAllele>],
    seams: &[Vec<GenomeSeam>],
) -> Vec<Vec<bool>> {
    let mut viable = partitions
        .iter()
        .map(|alleles| vec![false; alleles.len()])
        .collect::<Vec<_>>();
    viable.last_mut().unwrap().fill(true);
    for boundary in (0..seams.len()).rev() {
        for seam in &seams[boundary] {
            if viable[boundary + 1][seam.right] {
                viable[boundary][seam.left] = true;
            }
        }
    }
    viable
}

fn chain_resident_bytes() -> io::Result<u64> {
    let status = fs::read_to_string("/proc/self/status")?;
    let kib = status
        .lines()
        .find(|line| line.starts_with("VmRSS:"))
        .and_then(|line| line.split_whitespace().nth(1))
        .ok_or_else(|| invalid("VmRSS absent from /proc/self/status"))?
        .parse::<u64>()
        .map_err(io::Error::other)?;
    kib.checked_mul(1024)
        .ok_or_else(|| invalid("VmRSS overflow"))
}

/// statm resident-page count (second field) times the 4 KiB base page size:
/// a single short read and parse, cheap enough for hot-loop sampling.
fn statm_resident_bytes() -> io::Result<u64> {
    const BASE_PAGE_BYTES: u64 = 4096;
    let statm = fs::read_to_string("/proc/self/statm")?;
    let pages = statm
        .split_whitespace()
        .nth(1)
        .and_then(|pages| pages.parse::<u64>().ok())
        .ok_or_else(|| invalid("resident pages absent from /proc/self/statm"))?;
    pages
        .checked_mul(BASE_PAGE_BYTES)
        .ok_or_else(|| invalid("statm resident bytes overflow"))
}

/// Hot-loop RSS guard: samples resident bytes only once every `every`
/// iterations, so a stage that grows 30 GB inside one checkpoint interval is
/// still caught by the loop itself instead of waiting for the next stage
/// boundary. Tracks the observed peak and fails closed over the budget.
#[derive(Debug)]
pub struct PeriodicRssGuard {
    budget_bytes: Option<u64>,
    every: u64,
    iterations: u64,
    peak_bytes: u64,
}

impl PeriodicRssGuard {
    pub fn new(budget_bytes: Option<u64>, every: u64) -> Self {
        Self {
            budget_bytes,
            // A zero interval would never sample; clamp to unconditional.
            every: every.max(1),
            iterations: 0,
            peak_bytes: 0,
        }
    }

    pub fn peak_bytes(&self) -> u64 {
        self.peak_bytes
    }

    /// One hot-loop iteration. Samples /proc/self/statm only when the
    /// iteration count reaches the sampling interval.
    pub fn checkpoint(&mut self, stage: &str) -> io::Result<()> {
        self.iterations = self.iterations.saturating_add(1);
        if self.iterations % self.every != 0 {
            return Ok(());
        }
        self.sample(stage)
    }

    /// Unconditional sample (stage boundaries use this; the peak stays shared
    /// with the hot-loop samples).
    pub fn sample(&mut self, stage: &str) -> io::Result<()> {
        let Some(budget) = self.budget_bytes else {
            return Ok(());
        };
        let rss = statm_resident_bytes()?;
        self.peak_bytes = self.peak_bytes.max(rss);
        if rss > budget {
            return Err(invalid(&format!(
                "RSS guard exceeded during {stage}: {rss} > {budget}"
            )));
        }
        Ok(())
    }
}

fn enforce_chain_rss(budget: Option<u64>, stage: &str) -> io::Result<()> {
    let Some(budget) = budget else {
        return Ok(());
    };
    let rss = chain_resident_bytes()?;
    if rss > budget {
        return Err(invalid(&format!(
            "RSS guard exceeded during {stage}: {rss} > {budget}"
        )));
    }
    Ok(())
}

type SuccessorClass = (usize, Vec<(FeatureKey, u64)>);
type LocusSuccessors = (
    Vec<SuccessorClass>,
    BTreeMap<usize, BTreeMap<usize, Vec<usize>>>,
);

/// Per-locus machinery that is identical for every bounded pass over the
/// same partitions/seams (ladder seeds, main pass, tie rotations): the
/// successor registries, the registry weight tables, and the suffix
/// max-addition tables. Building them once per exact_chain_streaming call
/// and sharing them down the ladder removed the per-pass rebuilds that
/// dominated the slice DP wall (the seed ladder alone was >half of it).
struct SharedLocusMachinery {
    registries: Vec<LocusSuccessors>,
    weights: Vec<(Vec<f64>, Vec<f64>)>,
    max_additions: Vec<FxHashMap<u32, u64>>,
    /// Aggregate max-addition over all suffix loci (loci 1..=last): each
    /// pass starts from this total and subtracts departed loci as it goes.
    suffix_max_add_total: FxHashMap<u32, u64>,
    /// Lazily computed per-(locus, class-pair) fresh constants (see
    /// `PairFreshConstants`), shared by every DP pass over this machinery
    /// (seed ladder, main pass, rotations). The DP is single-threaded.
    pair_constants:
        std::cell::RefCell<FxHashMap<(usize, usize, usize), PairFreshConstants>>,
}

/// Compute the per-pair fresh constants with one streamed merge pass over
/// the delta: every merged feature contributes the state-independent
/// constants loss(d_f) / fp(fid, d_f) / +1 entry (loss(0) == 0.0 exactly).
fn compute_pair_fresh_constants(
    components: &[Arc<ConvertedProfile>; 4],
    universe: &FeatureUniverse,
    model: &ScoreModel,
) -> io::Result<PairFreshConstants> {
    let mut k_score = 0.0f64;
    let mut k_hash = 0u64;
    let mut fresh_entries = 0usize;
    let mut merge = ComponentMerge::new(components);
    while let Some((fid, added)) = merge.next()? {
        k_score += model.loss(added, universe.observed[fid as usize])?;
        k_hash ^= universe.fingerprint_value(fid, added);
        fresh_entries += 1;
    }
    Ok(PairFreshConstants {
        k_score,
        k_hash,
        fresh_entries,
    })
}

/// Fetch (or compute once) the per-(locus, class-pair) fresh constants.
/// With shared machinery the cache lives for the whole run (the components
/// of a given (locus, class1, class2) triple are content-identical in every
/// pass); without it a per-call map still dedups across the states of one
/// pass.
#[allow(clippy::too_many_arguments)]
fn pair_fresh_constants(
    shared: Option<&SharedLocusMachinery>,
    owned: &mut FxHashMap<(usize, usize, usize), PairFreshConstants>,
    locus: usize,
    first_class: usize,
    second_class: usize,
    components: &[Arc<ConvertedProfile>; 4],
    universe: &FeatureUniverse,
    model: &ScoreModel,
) -> io::Result<PairFreshConstants> {
    let key = (locus, first_class, second_class);
    if let Some(shared) = shared {
        if let Some(&constants) = shared.pair_constants.borrow().get(&key) {
            return Ok(constants);
        }
    } else if let Some(&constants) = owned.get(&key) {
        return Ok(constants);
    }
    let constants = compute_pair_fresh_constants(components, universe, model)?;
    if let Some(shared) = shared {
        shared.pair_constants.borrow_mut().insert(key, constants);
    } else {
        owned.insert(key, constants);
    }
    Ok(constants)
}

/// Build the locus successor structure: physical left allele -> registry
/// class id -> ascending right members. Registry classes are distinct
/// (profile class, seam signature) pairs across the boundary.
fn build_locus_successors(
    seams: &[GenomeSeam],
    memberships: &[usize],
    viable_right: &[bool],
) -> io::Result<(Vec<SuccessorClass>, BTreeMap<usize, BTreeMap<usize, Vec<usize>>>)> {
    let mut class_ids = HashMap::<SuccessorClass, usize>::new();
    let mut classes: Vec<SuccessorClass> = Vec::new();
    let mut successors = BTreeMap::<usize, BTreeMap<usize, Vec<usize>>>::new();
    for seam in seams {
        if viable_right[seam.right] {
            let key = (memberships[seam.right], profile_signature(&seam.profile));
            let next_id = classes.len();
            let class = *class_ids.entry(key.clone()).or_insert(next_id);
            if class == next_id {
                classes.push(key);
            }
            successors
                .entry(seam.left)
                .or_default()
                .entry(class)
                .or_default()
                .push(seam.right);
        }
    }
    for members in successors
        .values_mut()
        .flat_map(|classes| classes.values_mut())
    {
        // Ascending member order makes the physical tie-set enumeration
        // deterministic in BeamKey order, so the beam cutoff tail can be
        // bulk-dropped exactly.
        members.sort_unstable();
        members.dedup();
    }
    Ok((classes, successors))
}

/// Admissible future-cost ingredients for one locus transition: the linear
/// Poisson lower bound weights w_f = d * (1 - observed_f / background) make
/// every loss difference loss(q + delta) - loss(q) >= sum_f delta_f * w_f
/// (the loss derivative in q is increasing), so class and signature weight
/// sums lower-bound the delta and boundary contributions independently of
/// the unknown cumulative profile.
fn locus_weight_tables(
    profile_classes: &[Profile],
    classes: &[SuccessorClass],
    counts: &mut CountCache<'_>,
    model: &ScoreModel,
) -> io::Result<(Vec<f64>, Vec<f64>)> {
    let per_count = model.depth / 150.0;
    let mut observed = HashMap::<FeatureKey, f64>::new();
    let mut weight = |counts: &mut CountCache<'_>, feature: &FeatureKey| -> io::Result<f64> {
        let obs = *observed
            .entry(feature.clone())
            .or_insert(counts.get(feature)? as f64);
        Ok(per_count * (1.0 - obs / model.background))
    };
    let mut class_weights = vec![0.0; profile_classes.len()];
    for (class, profile) in profile_classes.iter().enumerate() {
        let mut total = 0.0;
        for (feature, count) in profile.iter() {
            total += weight(counts, feature)? * *count as f64;
        }
        class_weights[class] = total;
    }
    let mut signature_weights = vec![0.0; classes.len()];
    for (k, (_, signature)) in classes.iter().enumerate() {
        let mut total = 0.0;
        for (feature, count) in signature.iter() {
            total += weight(counts, feature)? * *count as f64;
        }
        signature_weights[k] = total;
    }
    Ok((class_weights, signature_weights))
}

/// Per-locus fid table for the tight A* bound: the maximum count any single
/// transition can add for the feature. A transition's delta is
/// class1 + class2 + seam1 + seam2, so the maximum addition is
/// 2 * max profile-class count + 2 * max seam-signature count.
fn locus_max_additions(
    profile_classes: &[Profile],
    classes: &[SuccessorClass],
    universe: &FeatureUniverse,
) -> io::Result<FxHashMap<u32, u64>> {
    let mut max_class = FxHashMap::<u32, u64>::default();
    for profile in profile_classes {
        for (feature, &count) in profile.iter() {
            let fid = universe.id(feature)?;
            let entry = max_class.entry(fid).or_default();
            *entry = (*entry).max(count);
        }
    }
    let mut max_signature = FxHashMap::<u32, u64>::default();
    for (_, signature) in classes {
        for (feature, count) in signature.iter() {
            let fid = universe.id(feature)?;
            let entry = max_signature.entry(fid).or_default();
            *entry = (*entry).max(*count);
        }
    }
    let mut table = FxHashMap::<u32, u64>::default();
    for (fid, count) in max_class {
        let signature_add = 2 * max_signature.get(&fid).copied().unwrap_or(0);
        table.insert(fid, 2 * count + signature_add);
    }
    for (fid, count) in max_signature {
        *table.entry(fid).or_default() += 2 * count;
    }
    Ok(table)
}

/// Fresh-feature potentials for one suffix state: c_f(0) =
/// min_{x in [0, S_f]} loss(x) per suffix feature, plus the constant
/// C = sum_f c_f(0). The state-aware remaining-cost bound is
/// tight(state) = C + correction(state), where correction accumulates
/// sum_f [c_f(q_f) - c_f(0)] over the state's touched suffix features along
/// the chain (admissible: the per-feature remaining cost
/// loss(final) - loss(q_f) is at least c_f(q_f) because counts only grow and
/// future additions are bounded by the suffix max-addition table).
fn fresh_potentials(
    suffix_max_add: &FxHashMap<u32, u64>,
    universe: &FeatureUniverse,
    model: &ScoreModel,
) -> (FxHashMap<u32, f64>, f64) {
    let mut c0 = FxHashMap::default();
    c0.reserve(suffix_max_add.len());
    let mut constant = 0.0f64;
    for (&fid, &max_add) in suffix_max_add {
        let value = feature_potential(
            0,
            max_add,
            universe.observed[fid as usize],
            model.histogram,
            model.depth,
            model.denominator,
            model.background,
        );
        c0.insert(fid, value);
        constant += value;
    }
    (c0, constant)
}

/// Bound-prune tolerance. The DP score is a telescoped sum of per-feature
/// loss increments while the bound recomputes losses directly from the
/// cumulative counts, so the two differ by f64 summation-order noise
/// (measured 5.3e-9 absolute on the slice's selected route). A near-optimal
/// realizable incumbent can sit within that noise of the true optimum, so a
/// zero-tolerance comparison prunes the optimal path. With the tolerance, a
/// pruned state's true completion satisfies
/// completion > incumbent + EPS - eta, where eta is the bound's own
/// overestimate — itself only summation noise bounded by EPS — so completion
/// > incumbent: no state that can beat the incumbent is ever pruned, and an
/// exhausted search is exact (complete: true is honest).
const BOUND_PRUNE_EPSILON: f64 = 1e-6;

fn bounded_conflict_chain(
    partitions: &[Vec<GenomeAllele>],
    seams: &[Vec<GenomeSeam>],
    viable: &[Vec<bool>],
    _incidence: &IncidenceTable,
    model: &ScoreModel,
    profile_memberships: &[Vec<usize>],
    profile_classes: &[Vec<Profile>],
    _local_losses: &[BTreeMap<(usize, usize), f64>],
    counts: &mut CountCache<'_>,
    mut accounting: ChainAccounting,
    force_initial_beam: bool,
    beam_width: usize,
    rss_budget_bytes: Option<u64>,
    inherited_incumbent: Option<f64>,
    seed_mode: bool,
    tie_rotation: u64,
    shared_universe: Option<&FeatureUniverse>,
    shared_machinery: Option<&SharedLocusMachinery>,
) -> io::Result<GenomeChain> {
    ensure(beam_width > 0, "zero conflict beam width")?;
    enforce_chain_rss(rss_budget_bytes, "bounded_chain_start")?;
    DpStageProbe::note_pass();
    accounting.beam_width = beam_width;
    accounting.bounded_blocks = 1;
    let state_bounded = force_initial_beam || accounting.conflict_window == 0;
    // The cumulative profile identity is part of the search key, so only the
    // suffix needed for physical overlap conflicts must remain in the key.
    let history_limit = accounting.conflict_window.max(1);
    let mut profile_interner = PersistentProfileInterner {
        telescoped: dev_flag("IMPG_DEV_TELESCOPED_EXTEND"),
        ..PersistentProfileInterner::default()
    };
    // Per-(locus, class-pair) fresh constants for callers without shared
    // machinery (tests); with shared machinery the run-wide RefCell cache
    // in SharedLocusMachinery is used instead.
    let mut owned_pair_constants =
        FxHashMap::<(usize, usize, usize), PairFreshConstants>::default();
    // Per-locus successor registries, built once per pass: the feature
    // universe, the suffix max-addition tables, and the expansion loop all
    // share them.
    let registries_started = Instant::now();
    let mut owned_registries;
    let locus_registries: &[LocusSuccessors] = match shared_machinery {
        Some(shared) => &shared.registries,
        None => {
            owned_registries = Vec::with_capacity(partitions.len());
            owned_registries.push((Vec::new(), BTreeMap::new()));
            for locus in 1..partitions.len() {
                owned_registries.push(build_locus_successors(
                    &seams[locus - 1],
                    &profile_memberships[locus],
                    &viable[locus],
                )?);
            }
            &owned_registries
        }
    };
    DpStageProbe::add(0, registries_started);
    // The universe resolves every feature's observed count once; rebuilding
    // it per seed pass re-queried the full feature set through the count
    // cache each time (and thrashed its capacity cap on component-scale
    // feature universes), so the top-level build is shared down the ladder.
    let universe_started = Instant::now();
    let owned_universe;
    let universe = match shared_universe {
        Some(shared) => shared,
        None => {
            owned_universe = FeatureUniverse::build(profile_classes, &locus_registries, counts)?;
            &owned_universe
        }
    };
    DpStageProbe::add(1, universe_started);
    let empty_component = empty_component();
    let initial_started = Instant::now();
    let mut converted_locus0 = vec![None; profile_classes[0].len()];
    let mut initial = BTreeMap::new();
    for first in 0..partitions[0].len() {
        if !viable[0][first] || !spans_feasible(&partitions[0][first].traversal.segments) {
            continue;
        }
        for second in first..partitions[0].len() {
            if !viable[0][second] || !spans_feasible(&partitions[0][second].traversal.segments) {
                continue;
            }
            let pair = [first, second];
            let first_class = profile_memberships[0][first];
            let second_class = profile_memberships[0][second];
            // The initial delta is the pair of locus-0 class profiles; its
            // component quadruple needs no merged materialization.
            let first_component = converted_locus0[first_class]
                .get_or_insert_with(|| {
                    convert_profile(&profile_classes[0][first_class], &universe)
                        .expect("converted locus-0 class")
                })
                .clone();
            let second_component = converted_locus0[second_class]
                .get_or_insert_with(|| {
                    convert_profile(&profile_classes[0][second_class], &universe)
                        .expect("converted locus-0 class")
                })
                .clone();
            let components = [
                first_component,
                second_component,
                Arc::clone(&empty_component),
                Arc::clone(&empty_component),
            ];
            let constants = pair_fresh_constants(
                shared_machinery,
                &mut owned_pair_constants,
                0,
                first_class,
                second_class,
                &components,
                &universe,
                model,
            )?;
            let (score, cumulative) = profile_interner.extend(
                None,
                None,
                components,
                0.0,
                constants,
                &universe,
                model,
            )?;
            initial.insert(
                BeamKey {
                    history: vec![pair],
                    profile_id: cumulative.id,
                },
                BeamCandidate {
                    score,
                    choices: pair,
                    predecessors: Vec::new(),
                    cumulative_profile: cumulative,
                },
            );
        }
    }
    let mut nodes = Vec::<PhysicalNode>::new();
    let mut layer = BTreeMap::<BeamKey, usize>::new();
    let initial = if state_bounded {
        prune_beam(initial, &mut accounting, beam_width, tie_rotation)
    } else {
        initial.into_iter().collect()
    };
    let mut layer_profiles = BTreeMap::<usize, Arc<PersistentProfile>>::new();
    // Per-state cumulative count maps (fid -> count) are built AFTER the
    // logical-state budget check: an unbounded initial layer must fail the
    // budget honestly before any per-state map allocates.
    for (key, candidate) in initial {
        let node = nodes.len();
        layer_profiles.insert(node, Arc::clone(&candidate.cumulative_profile));
        nodes.push(PhysicalNode {
            score: candidate.score,
            choices: candidate.choices,
            predecessors: candidate.predecessors,
        });
        layer.insert(key, node);
    }
    accounting.occupancy.push(layer.len());
    accounting.state_bytes = physical_state_bytes(layer.len(), &nodes)
        .saturating_add(cumulative_profile_bytes(&layer_profiles));
    if accounting.state_bytes > MAX_STATE_BYTES {
        return Ok(GenomeChain {
            complete: false,
            stop_reason: "budget: exact_conflict_prefix_state".into(),
            choices: Vec::new(),
            proposal_losses: Vec::new(),
            rotation_best_scores: Vec::new(),
            accounting,
        });
    }
    let mut layer_maps = BTreeMap::<usize, PersistentCountMap>::new();
    for &node in layer.values() {
        let profile = layer_profiles
            .get(&node)
            .expect("bounded state cumulative profile");
        let mut entries = Vec::new();
        let mut merge = ComponentMerge::new(&profile.components);
        loop {
            match merge.next().expect("initial state map merge") {
                Some((fid, count)) => entries.push((fid, count)),
                None => break,
            }
        }
        layer_maps.insert(node, PersistentCountMap::from_sorted(&entries));
    }
    profile_interner.clear_search_indexes();
    DpStageProbe::add(2, initial_started);
    // A* bounds: a realizable incumbent plus per-state tight admissible
    // the minimum of the (convex) Poisson loss over the achievable count
    // interval [q, q + suffix max-addition]; see tight_profile_future. States
    // and edges whose bounded completion cannot beat the incumbent are
    // skipped and counted; pruning is by admissible bound with a summation
    // noise tolerance (BOUND_PRUNE_EPSILON), so every state that can still
    // beat (or tie) the incumbent survives and an exhausted search is exact.
    // IMPG_NO_ASTAR_BOUND disables the incumbent and every bound check: the
    // pure exhaustive beam DP used to adversarially verify completeness.
    let bound_disabled = std::env::var("IMPG_NO_ASTAR_BOUND").is_ok();
    let mut incumbent: Option<f64> = if bound_disabled {
        None
    } else {
        // The inherited incumbent is a realizable completed chain's score
        // (an externally rescored reference route for seeded runs, or the
        // previous ladder seed's / rotation pass's best), so the admissible
        // bound can prune against it honestly.
        accounting.greedy_incumbent = inherited_incumbent;
        inherited_incumbent
    };
    // Seed ladder: progressive narrow-beam probes give the main pass an
    // incumbent the admissible state bound can actually engage against. Every
    // seed is a real beam DP over physical conflict-checked states, so its
    // completed chain is realizable and its score is a valid incumbent. Each
    // seed inherits the previous seed's incumbent, so with the bound engaged
    // the wider seeds are themselves pruned passes. A seed's drops and bound
    // counters belong to the probe, not the reported search, and are
    // discarded with its accounting.
    // IMPG_DEV_SKIP_SEED_LADDER (dev-only): skip the ladder entirely. At
    // slice scale unseeded, the ladder's incumbents never prune a main-pass
    // state (bound_pruned_states = 0 in every unseeded slice run so far), so
    // it is pure overhead in the dev loop; at full-component scale it is the
    // mechanism that makes the A* bound engage, so production keeps it and
    // reported greedy_incumbent stays honest.
    if !seed_mode && !bound_disabled && !dev_flag("IMPG_DEV_SKIP_SEED_LADDER") {
        let ladder_started = Instant::now();
        for &seed_width in &[1usize, 8, 64, 256] {
            if seed_width >= beam_width {
                break;
            }
            let seed_accounting = accounting.clone();
            let seed = bounded_conflict_chain(
                partitions,
                seams,
                viable,
                _incidence,
                model,
                profile_memberships,
                profile_classes,
                _local_losses,
                counts,
                seed_accounting,
                force_initial_beam,
                seed_width,
                rss_budget_bytes,
                incumbent,
                true,
                tie_rotation,
                Some(universe),
                shared_machinery,
            )?;
            if let Some(score) = seed.proposal_losses.first().copied() {
                let improved = incumbent
                    .map_or(true, |current| score.total_cmp(&current).is_lt());
                if improved {
                    incumbent = Some(score);
                    accounting.greedy_incumbent = Some(score);
                }
            } else {
                // The seed dead-ended at this width; wider seeds may still
                // complete, so continue up the ladder.
                continue;
            }
        }
        DpStageProbe::add(3, ladder_started);
    }
    // Running fid-keyed suffix max-addition table: while expanding locus l it
    // holds, per feature, the maximum count any transition from locus l
    // onward can add. It starts as the full sum over loci 1..=last and each
    // departed locus is subtracted at the end of its iteration.
    let suffix_tables_started = Instant::now();
    let mut owned_max_add;
    let mut suffix_max_add = FxHashMap::<u32, u64>::default();
    let locus_max_add: &[FxHashMap<u32, u64>] = match shared_machinery {
        Some(shared) => {
            suffix_max_add = shared.suffix_max_add_total.clone();
            &shared.max_additions
        }
        None => {
            owned_max_add = Vec::with_capacity(partitions.len());
            owned_max_add.push(FxHashMap::default());
            for locus in 1..partitions.len() {
                let (classes, _) = &locus_registries[locus];
                let table =
                    locus_max_additions(&profile_classes[locus], classes, &universe)?;
                for (fid, max_add) in &table {
                    *suffix_max_add.entry(*fid).or_default() += max_add;
                }
                owned_max_add.push(table);
            }
            &owned_max_add
        }
    };
    // State-aware remaining-cost bound. tight(state) = C + correction(state):
    // C is the suffix's fresh-feature potential constant, correction(state)
    // accumulates sum_f [c_f(q_f) - c_f(0)] over the state's touched suffix
    // features incrementally along the chain, and both are re-based when the
    // suffix shrinks at a locus boundary. Admissibility: the remaining cost
    // loss(final) - loss(q) is at least sum_f c_f(q_f) because counts only
    // grow and future additions per feature are bounded by the suffix
    // max-addition table; the correction's per-feature increments re-base
    // exactly at every suffix shrink, so the invariant holds for the current
    // suffix. Prune comparisons carry BOUND_PRUNE_EPSILON (f64 summation
    // noise between the telescoped score and the recomputed bound).
    let (mut c0_current, mut suffix_constant) = fresh_potentials(&suffix_max_add, &universe, model);
    DpStageProbe::add(5, suffix_tables_started);
    let mut corrections = BTreeMap::<usize, f64>::new();
    for (_, &node) in &layer {
        let map = &layer_maps[&node];
        let mut correction = 0.0f64;
        for (fid, q) in map.iter() {
            if q == 0 {
                continue;
            }
            if let Some(&max_add) = suffix_max_add.get(&fid) {
                let base = c0_current
                    .get(&fid)
                    .copied()
                    .expect("suffix feature has a fresh potential");
                correction += feature_potential(
                    q,
                    max_add,
                    universe.observed[fid as usize],
                    model.histogram,
                    model.depth,
                    model.denominator,
                    model.background,
                ) - base;
            }
        }
        corrections.insert(node, correction);
    }
    for locus in 1..partitions.len() {
        enforce_chain_rss(rss_budget_bytes, "bounded_chain_locus_start")?;
        let probe_logging = DpProbe::on();
        let probe_snapshot = DpProbe::snapshot();
        // Class-graph registry: distinct (profile class, seam signature)
        // successor classes for this boundary. Expansion iterates
        // (state-profile-identity x successor-class) pairs and memoizes each
        // distinct transition; physical members materialize only as ordered
        // tie-set enumerations under the beam cutoff. There is no global C^2
        // class-pair plan: each state expands exactly the class pairs its own
        // physical exit pair can reach.
        let (classes, successors) = &locus_registries[locus];
        let weight_tables_started = Instant::now();
        let owned_weights;
        let (class_weights, signature_weights): (&[f64], &[f64]) = match shared_machinery {
            Some(shared) => (&shared.weights[locus].0, &shared.weights[locus].1),
            None => {
                owned_weights = locus_weight_tables(&profile_classes[locus], classes, counts, model)?;
                (&owned_weights.0, &owned_weights.1)
            }
        };
        DpStageProbe::add(4, weight_tables_started);
        // Registry-class weight w_k = class weight + signature weight: the
        // linear per-class-pair delta lower bound is w_k1 + w_k2 (used only
        // to order each state's class pairs; the state-level bound below is
        // the tight one).
        let mut registry_weights = vec![0.0f64; classes.len()];
        for (k, (profile_class, _)) in classes.iter().enumerate() {
            registry_weights[k] = class_weights[*profile_class] + signature_weights[k];
        }
        // Lazily converted profile classes and seam signatures for this
        // locus: the transition delta is the component quadruple
        // (class1, class2, signature1, signature2), never a merged map.
        let mut converted_profile_classes = vec![None; profile_classes[locus].len()];
        let mut converted_signatures = vec![None; classes.len()];
        // Per-source member spans sorted by span start: a past allele segment
        // can only overlap members whose span starts below the segment end.
        let mut spans_by_source = BTreeMap::<usize, Vec<(u64, u64, usize)>>::new();
        if accounting.conflict_window > 0 {
            for (member, allele) in partitions[locus].iter().enumerate() {
                let mut start = u64::MAX;
                let mut end = 0u64;
                for segment in &allele.traversal.segments {
                    start = start.min(segment.start);
                    end = end.max(segment.end);
                }
                spans_by_source
                    .entry(allele.traversal.entry().source)
                    .or_default()
                    .push((start, end, member));
            }
            for spans in spans_by_source.values_mut() {
                spans.sort_unstable();
            }
        }
        let mut candidates = BoundedSuccessors::new(beam_width, tie_rotation);
        // Hot-loop RSS guard: sampled every 1,024 (state, class) expansions via
        // /proc/self/statm, so a single allocation spike inside one locus is
        // caught by the loop rather than at the next stage boundary.
        let mut expansion_rss = PeriodicRssGuard::new(rss_budget_bytes, 1_024);
        // Per-state conflict eligibility, hoisted out of the physical member
        // enumeration: eligibility of a member allele depends only on the
        // state's conflict-window history and the member's own segments.
        let mut eligibility_cache = HashMap::<usize, [Vec<bool>; 2]>::new();
        let mut evaluated_state_classes = 0usize;
        let mut transition_plan_len_dbg = 0usize;
        // Probe-only: per-state size of the touched-recurring set
        // {fid: q > 0} intersected with the suffix feature table — the set
        // an O(overlap) telescoped extend would have to iterate per (state,
        // class pair).
        let mut recur_states = 0usize;
        let mut recur_sum = 0usize;
        let mut recur_max = 0usize;
        let mut map_len_sum = 0usize;
        let mut map_len_max = 0usize;
        // Probe-only: per-state size of {fid: q > 0} intersected with THIS
        // boundary's delta feature universe (the union of every successor
        // class's profile-class and seam-signature features) — the true
        // per-(state, locus) overlap set an O(overlap) telescoped extend
        // would iterate per class pair, after one per-state filter pass.
        let mut overlap_states = 0usize;
        let mut overlap_sum = 0usize;
        let mut overlap_max = 0usize;
        let mut locus_delta_universe: Option<FxHashSet<u32>> = None;
        // Reusable per-state reachable class-pair buffer, sorted by the
        // admissible weight bound so the early break is monotone.
        let mut reachable = Vec::<(f64, usize, usize)>::new();
        let expansion_started = Instant::now();
        for (state_key, &previous_node) in &layer {
            let history = &state_key.history;
            let previous_pair = history.last().expect("nonempty conflict history");
            let Some(first_adjacent) = successors.get(&previous_pair[0]) else {
                continue;
            };
            let Some(second_adjacent) = successors.get(&previous_pair[1]) else {
                continue;
            };
            let previous_profile = layer_profiles
                .get(&previous_node)
                .expect("bounded state cumulative profile");
            // This boundary's delta feature universe: the union of every
            // registry class's profile-class and seam-signature features,
            // built from the converted (fid-keyed) components. Any delta
            // feature of any reachable pair is in this set by construction,
            // so a state's touched count that is absent from the small map
            // below is exactly q == 0 (fresh) — the extend hot path probes
            // the small map instead of the persistent radix map.
            let locus_universe = locus_delta_universe.get_or_insert_with(|| {
                let mut set = FxHashSet::default();
                for (profile_class, _) in classes.iter() {
                    let component = converted_profile_classes[*profile_class]
                        .get_or_insert_with(|| {
                            convert_profile(&profile_classes[locus][*profile_class], &universe)
                                .expect("converted profile class")
                        })
                        .clone();
                    for &(fid, _) in component.iter() {
                        set.insert(fid);
                    }
                }
                for (class_index, (_, signature)) in classes.iter().enumerate() {
                    let component = converted_signatures[class_index]
                        .get_or_insert_with(|| {
                            convert_signature(signature, &universe)
                                .expect("converted seam signature")
                        })
                        .clone();
                    for &(fid, _) in component.iter() {
                        set.insert(fid);
                    }
                }
                set
            });
            // Per-state touched counts restricted to this boundary's delta
            // universe: the O(overlap) correction set. Measured tiny exactly
            // where the extend volume concentrates (loci 1-2 of the slice:
            // 73 and 68 per state vs ~6.8k-feature deltas).
            let mut small = FxHashMap::<u32, u64>::default();
            {
                let map = layer_maps
                    .get(&previous_node)
                    .expect("bounded state cumulative map");
                for (fid, q) in map.iter() {
                    if q > 0 && locus_universe.contains(&fid) {
                        small.insert(fid, q);
                    }
                }
            }
            if probe_logging {
                let map = layer_maps
                    .get(&previous_node)
                    .expect("bounded state cumulative map");
                let map_len = map.len();
                let mut recur = 0usize;
                for (fid, q) in map.iter() {
                    if q > 0 && suffix_max_add.contains_key(&fid) {
                        recur += 1;
                    }
                }
                recur_states += 1;
                recur_sum += recur;
                recur_max = recur_max.max(recur);
                overlap_states += 1;
                overlap_sum += small.len();
                overlap_max = overlap_max.max(small.len());
                map_len_sum += map_len;
                map_len_max = map_len_max.max(map_len);
            }
            if let Some(incumbent) = incumbent {
                // State-aware admissible bound: tight = suffix fresh constant
                // plus the state's accumulated touched-feature correction.
                // A state whose score plus that bound cannot beat the
                // incumbent is skipped. Exact-best ties survive: the
                // comparison carries the summation-noise tolerance
                // BOUND_PRUNE_EPSILON, and the bound's own overestimate is
                // within that tolerance, so a pruned state's true completion
                // is strictly worse than the incumbent.
                let tight = suffix_constant
                    + corrections
                        .get(&previous_node)
                        .copied()
                        .expect("layer state carries a bound correction");
                if nodes[previous_node].score + tight > incumbent + BOUND_PRUNE_EPSILON {
                    if probe_logging {
                        eprintln!(
                            "prune locus {} score {:.4} tight {:.4} incumbent {:.4}",
                            locus,
                            nodes[previous_node].score,
                            tight,
                            incumbent
                        );
                    }
                    accounting.bound_pruned_states += 1;
                    // Every class pair this state would have expanded is a
                    // pruned edge; the adjacency product is exactly the
                    // reachable list the expansion would have built.
                    accounting.bound_pruned_edges = accounting
                        .bound_pruned_edges
                        .saturating_add((first_adjacent.len() * second_adjacent.len()) as u64);
                    continue;
                }
            }
            reachable.clear();
            for (&first_class, _) in first_adjacent {
                let first_weight = registry_weights[first_class];
                for (&second_class, _) in second_adjacent {
                    let weight_bound = first_weight + registry_weights[second_class];
                    reachable.push((weight_bound, first_class, second_class));
                }
            }
            // The inner break only prunes suffixes of the second adjacency for
            // a fixed first class; sort so the cheapest pairs come first and
            // the monotone early break applies across the whole state.
            reachable.sort_by(|a, b| {
                a.0.total_cmp(&b.0)
                    .then_with(|| (a.1, a.2).cmp(&(b.1, b.2)))
            });
            transition_plan_len_dbg += reachable.len();
            for &(_weight_bound, first_class, second_class) in &reachable {
                expansion_rss.checkpoint("bounded_chain_expansion_hot_loop")?;
                evaluated_state_classes += 1;
                let first_members = &first_adjacent[&first_class];
                let second_members = &second_adjacent[&second_class];
                // The transition delta is the component quadruple (class1,
                // class2, signature1, signature2) — no merged profile is
                // materialized; the interner streams the four-way merge and
                // memoizes each distinct (parent profile, component
                // quadruple) computation.
                let (first_profile_class, first_signature) = &classes[first_class];
                let (second_profile_class, second_signature) = &classes[second_class];
                let first_component = converted_profile_classes[*first_profile_class]
                    .get_or_insert_with(|| {
                        convert_profile(&profile_classes[locus][*first_profile_class], &universe)
                            .expect("converted profile class")
                    })
                    .clone();
                let second_component = converted_profile_classes[*second_profile_class]
                    .get_or_insert_with(|| {
                        convert_profile(&profile_classes[locus][*second_profile_class], &universe)
                            .expect("converted profile class")
                    })
                    .clone();
                let first_boundary = converted_signatures[first_class]
                    .get_or_insert_with(|| {
                        convert_signature(first_signature, &universe)
                            .expect("converted seam signature")
                    })
                    .clone();
                let second_boundary = converted_signatures[second_class]
                    .get_or_insert_with(|| {
                        convert_signature(second_signature, &universe)
                            .expect("converted seam signature")
                    })
                    .clone();
                let components = [
                    first_component,
                    second_component,
                    first_boundary,
                    second_boundary,
                ];
                let constants = pair_fresh_constants(
                    shared_machinery,
                    &mut owned_pair_constants,
                    locus,
                    first_class,
                    second_class,
                    &components,
                    &universe,
                    model,
                )?;
                let (total, cumulative_profile) = profile_interner.extend(
                    Some(previous_profile),
                    Some(&small),
                    components,
                    nodes[previous_node].score,
                    constants,
                    &universe,
                    model,
                )?;
                // Exact per-copy member eligibility for this state's history.
                let eligible = if accounting.conflict_window > 0 {
                    Some(
                        eligibility_cache
                            .entry(previous_node)
                            .or_insert_with(|| {
                                [
                                    conflict_eligible_members(
                                        partitions,
                                        locus,
                                        history,
                                        0,
                                        &spans_by_source,
                                    ),
                                    conflict_eligible_members(
                                        partitions,
                                        locus,
                                        history,
                                        1,
                                        &spans_by_source,
                                    ),
                                ]
                            }),
                    )
                } else {
                    None
                };
                let eligible = eligible.as_deref();
                let member_stats = |copy: usize, members: &[usize]| -> (u64, usize) {
                    match eligible {
                        None => (
                            members.len() as u64,
                            members.first().copied().unwrap_or(usize::MAX),
                        ),
                        Some(eligible) => {
                            let mut count = 0u64;
                            let mut first = usize::MAX;
                            for &member in members {
                                if eligible[copy][member] {
                                    count += 1;
                                    if first == usize::MAX {
                                        first = member;
                                    }
                                }
                            }
                            (count, first)
                        }
                    }
                };
                let (first_count, first_start) = member_stats(0, first_members);
                let (second_count, second_start) = member_stats(1, second_members);
                let feasible = first_count.saturating_mul(second_count);
                if candidates.candidates.len() == beam_width {
                    let worst = candidates
                        .ranking
                        .last()
                        .expect("full bounded successor layer");
                    let below_cutoff = first_start != usize::MAX
                        && second_start != usize::MAX
                        && {
                            let mut minimum_history = history.clone();
                            minimum_history.push([first_start, second_start]);
                            if minimum_history.len() > history_limit {
                                minimum_history.remove(0);
                            }
                            let minimum_key = BeamKey {
                                history: minimum_history,
                                profile_id: cumulative_profile.id,
                            };
                            ScoreKey(total).cmp(&worst.0).is_lt()
                                || {
                                    // Same score: the tie order (rotated)
                                    // decides, exactly as the historical
                                    // lexicographic tuple compare did at
                                    // rotation 0.
                                    ScoreKey(total) == worst.0
                                        && cmp_ranked(
                                            tie_rotation,
                                            &minimum_key,
                                            &worst.1.key,
                                        )
                                        .is_lt()
                                }
                        };
                    if !below_cutoff {
                        // Every physical member pair is at or beyond the beam
                        // cutoff exactly as per-pair rejection would decide.
                        candidates.record_bulk_drop(total, feasible);
                        accounting.work = accounting.work.saturating_add(feasible);
                        continue;
                    }
                }
                accounting.work = accounting.work.saturating_add(feasible);
                let mut enumerated = 0u64;
                'members: for &first in first_members {
                    if eligible.is_some_and(|eligible| !eligible[0][first]) {
                        continue;
                    }
                    for &second in second_members {
                        if eligible.is_some_and(|eligible| !eligible[1][second]) {
                            continue;
                        }
                        enumerated += 1;
                        let pair = [first, second];
                        let mut next_history = history.clone();
                        next_history.push(pair);
                        if next_history.len() > history_limit {
                            next_history.remove(0);
                        }
                        let key = BeamKey {
                            history: next_history,
                            profile_id: cumulative_profile.id,
                        };
                        if candidates.candidates.len() == beam_width {
                            let worst = candidates
                                .ranking
                                .last()
                                .expect("full bounded successor layer");
                            let beyond_cutoff = ScoreKey(total).cmp(&worst.0).is_gt()
                                || {
                                    // Same score: the tie order (rotated)
                                    // decides.
                                    ScoreKey(total) == worst.0
                                        && cmp_ranked(tie_rotation, &key, &worst.1.key).is_gt()
                                };
                            if beyond_cutoff {
                                // Ascending key order: this pair and every
                                // remaining pair are strictly beyond the beam
                                // cutoff, so reject the tail in one step.
                                let remainder = feasible.saturating_sub(enumerated - 1);
                                candidates.record_bulk_drop(total, remainder);
                                break 'members;
                            }
                        }
                        candidates.insert(
                            key,
                            total,
                            pair,
                            previous_node,
                            Arc::clone(&cumulative_profile),
                            &mut accounting,
                        );
                    }
                }
            }
        }
        DpStageProbe::add(6, expansion_started);
        if candidates.candidates.is_empty() {
            accounting.sample_count_queries = counts.queries;
            accounting.sample_count_cache_resets = counts.resets;
            return Ok(GenomeChain {
                complete: false,
                stop_reason: "bounded_conflict_beam_no_retained_legal_state".into(),
                choices: Vec::new(),
                proposal_losses: Vec::new(),
                rotation_best_scores: Vec::new(),
                accounting,
            });
        }
        if std::env::var("IMPG_DB").is_ok() {
            eprintln!(
                "PASS{} locus {} states {} plan_pairs {} kept {} drops_t {} drops_nt {} suffix_feats {} recur_states {} recur_sum {} recur_max {} map_len_sum {} map_len_max {} ov_states {} ov_sum {} ov_max {}",
                beam_width,
                locus,
                layer.len(),
                transition_plan_len_dbg,
                candidates.candidates.len(),
                candidates.dropped_ties,
                candidates.dropped_non_ties,
                suffix_max_add.len(),
                recur_states,
                recur_sum,
                recur_max,
                map_len_sum,
                map_len_max,
                overlap_states,
                overlap_sum,
                overlap_max,
            );
            DpProbe::print_delta(&format!("locus {locus}"), &probe_snapshot);
        }
        accounting.dropped_ties = accounting
            .dropped_ties
            .saturating_add(candidates.dropped_ties);
        accounting.dropped_non_ties = accounting
            .dropped_non_ties
            .saturating_add(candidates.dropped_non_ties);
        let mut next = BTreeMap::new();
        let commit_started = Instant::now();
        let mut next_profiles = BTreeMap::new();
        let mut next_maps = BTreeMap::<usize, PersistentCountMap>::new();
        let mut next_corrections = BTreeMap::<usize, f64>::new();
        for (_, key) in candidates.ranking {
            let candidate = candidates
                .candidates
                .remove(&key.key)
                .expect("ranked bounded successor");
            let node = nodes.len();
            // Child count map and bound correction: both derive from the
            // first-inserted predecessor's state (the stored cumulative
            // profile is that predecessor extended by the transition's
            // components). One streamed merge pass applies the delta to the
            // parent map and accumulates the touched-feature correction
            // increments [c(q+d) - c(q)] over the current suffix.
            let predecessor_node = candidate
                .predecessors
                .first()
                .copied()
                .expect("bounded successor carries a predecessor");
            let parent_map = layer_maps
                .get(&predecessor_node)
                .expect("predecessor state map");
            let parent_correction = corrections
                .get(&predecessor_node)
                .copied()
                .expect("predecessor bound correction");
            let mut correction = parent_correction;
            let mut merge = ComponentMerge::new(&candidate.cumulative_profile.components);
            let mut delta = Vec::new();
            loop {
                let next_entry = merge.next()?;
                let Some((fid, added)) = next_entry else {
                    break;
                };
                let previous = parent_map.get(fid);
                let updated = previous
                    .checked_add(added)
                    .ok_or_else(|| invalid("cumulative profile overflow"))?;
                // The persistent map applies the ADDITION itself (its leaf
                // fold adds to the parent's stored count), so the delta
                // carries `added`, not the cumulative value.
                delta.push((fid, added));
                if previous == 0 && updated == 0 {
                    continue;
                }
                if let Some(&max_add) = suffix_max_add.get(&fid) {
                    if max_add > 0 {
                        let observed = universe.observed[fid as usize];
                        correction += feature_potential(
                            updated,
                            max_add,
                            observed,
                            model.histogram,
                            model.depth,
                            model.denominator,
                            model.background,
                        ) - feature_potential(
                            previous,
                            max_add,
                            observed,
                            model.histogram,
                            model.depth,
                            model.denominator,
                            model.background,
                        );
                    }
                }
            }
            // Structural sharing: the child map copies only the delta's
            // root-to-leaf paths and shares every other node with the
            // parent map. O(delta), not O(profile).
            next_maps.insert(node, parent_map.apply_delta(&delta)?);
            next_corrections.insert(node, correction);
            next_profiles.insert(node, Arc::clone(&candidate.cumulative_profile));
            nodes.push(PhysicalNode {
                score: candidate.score,
                choices: candidate.choices,
                predecessors: candidate.predecessors,
            });
            next.insert(key.key, node);
        }
        accounting.occupancy.push(next.len());
        accounting.state_bytes = accounting.state_bytes.max(
            physical_state_bytes(next.len(), &nodes)
                .saturating_add(cumulative_profile_bytes(&next_profiles)),
        );
        if accounting.state_bytes > MAX_STATE_BYTES {
            accounting.sample_count_queries = counts.queries;
            accounting.sample_count_cache_resets = counts.resets;
            return Ok(GenomeChain {
                complete: false,
                stop_reason: "budget: bounded_conflict_beam_state".into(),
                choices: Vec::new(),
                proposal_losses: Vec::new(),
                rotation_best_scores: Vec::new(),
                accounting,
            });
        }
        layer = next;
        layer_profiles = next_profiles;
        layer_maps = next_maps;
        corrections = next_corrections;
        profile_interner.clear_search_indexes();
        DpStageProbe::add(7, commit_started);
        // Shrink the suffix by the departed locus and re-base the layer's
        // bound corrections to the new suffix: a departed touched feature
        // loses its whole term, a retained one re-bases both its endpoint
        // potential and its fresh baseline.
        let rebase_started = Instant::now();
        let mut shrunk = Vec::<(u32, u64, u64, f64, f64, u64)>::new();
        for (fid, &locus_add) in &locus_max_add[locus] {
            let old_add = suffix_max_add
                .get(fid)
                .copied()
                .expect("departed locus feature is in the suffix");
            let new_add = old_add.saturating_sub(locus_add);
            if new_add == old_add {
                continue;
            }
            suffix_max_add.insert(*fid, new_add);
            let old_base = c0_current
                .get(fid)
                .copied()
                .expect("suffix feature has a fresh potential");
            let observed = universe.observed[*fid as usize];
            let new_base = feature_potential(
                0,
                new_add,
                observed,
                model.histogram,
                model.depth,
                model.denominator,
                model.background,
            );
            shrunk.push((*fid, old_add, new_add, old_base, new_base, observed));
        }
        for (_, &node) in &layer {
            let map = layer_maps.get(&node).expect("layer state map");
            let Some(correction) = corrections.get_mut(&node) else {
                continue;
            };
            for &(fid, old_add, new_add, old_base, new_base, observed) in &shrunk {
                let q = map.get(fid);
                if q == 0 {
                    continue;
                }
                if new_add == 0 {
                    // Departed: the feature can no longer change, so its
                    // whole term leaves the suffix sum.
                    *correction -= feature_potential(
                        q,
                        old_add,
                        observed,
                        model.histogram,
                        model.depth,
                        model.denominator,
                        model.background,
                    ) - old_base;
                } else {
                    *correction += feature_potential(
                        q,
                        new_add,
                        observed,
                        model.histogram,
                        model.depth,
                        model.denominator,
                        model.background,
                    ) - feature_potential(
                        q,
                        old_add,
                        observed,
                        model.histogram,
                        model.depth,
                        model.denominator,
                        model.background,
                    ) - (new_base - old_base);
                }
            }
        }
        let (next_c0, next_constant) = fresh_potentials(&suffix_max_add, &universe, model);
        c0_current = next_c0;
        suffix_constant = next_constant;
        DpStageProbe::add(8, rebase_started);
    }
    let terminal_started = Instant::now();
    let best = layer
        .values()
        .map(|&node| nodes[node].score)
        .min_by(f64::total_cmp)
        .ok_or_else(|| invalid("empty bounded DP terminal layer"))?;
    let mut choices = Vec::new();
    let mut truncated = false;
    for node in layer
        .values()
        .copied()
        .filter(|&node| nodes[node].score.total_cmp(&best).is_eq())
    {
        backtrack_physical(&nodes, node, &mut Vec::new(), &mut choices, &mut truncated);
    }
    if let Some(selected) = choices.first() {
        // The per-locus selected-pair reporting needs the all-pairs local
        // loss table; dev-lazy-evidence runs carry no evidence rows and skip
        // the reporting entirely.
        if !accounting.local_pair_evidence.is_empty() {
            for locus in 0..selected[0].len() {
                let first = profile_memberships[locus][selected[0][locus]];
                let second = profile_memberships[locus][selected[1][locus]];
                let key = if first <= second {
                    (first, second)
                } else {
                    (second, first)
                };
                let loss = _local_losses[locus][&key];
                accounting.local_pair_evidence[locus].selected_loss = Some(loss);
                accounting.local_pair_evidence[locus].selected_minus_best =
                    Some(loss - accounting.local_pair_evidence[locus].best_loss);
            }
        }
    }
    accounting.sample_count_queries = counts.queries;
    accounting.sample_count_cache_resets = counts.resets;
    // Honest completeness: with NO beam drops anywhere (the initial layer's
    // forced prune and every locus's successor beam only record drops when
    // they actually truncate), no tie-enumeration truncation, and no budget
    // stop, the only pruning was admissible-bound pruning under
    // BOUND_PRUNE_EPSILON — the returned optimum is the true optimum of the
    // class-graph objective and complete is honest. Any beam drop, tie
    // truncation, or budget stop stays incomplete.
    let exhaustive = !truncated
        && accounting.dropped_ties == 0
        && accounting.dropped_non_ties == 0;
    DpStageProbe::add(9, terminal_started);
    Ok(GenomeChain {
        complete: exhaustive,
        stop_reason: if truncated {
            "budget: bounded_physical_beam_tie_enumeration".into()
        } else if exhaustive {
            "bounded_conflict_beam_exhausted_within_admissible_bound".into()
        } else if state_bounded {
            "bounded_physical_state_beam".into()
        } else {
            "bounded_conflict_beam".into()
        },
        proposal_losses: vec![best; choices.len()],
        choices,
        accounting,
        rotation_best_scores: vec![Some(best)],
    })
}

/// Run the bounded conflict chain under `tie_rotations` deterministic tie
/// orders, starting from a realizable `inherited_incumbent`. Rotation 0 is
/// the historical BeamKey order (bit-identical single-pass behavior when
/// tie_rotations <= 1); each further rotation re-runs the main-width pass
/// with a different — still deterministic — retained subset of over-width
/// tie classes, skipping the seed ladder (the incumbent from the previous
/// pass is already realizable and at least as good as any ladder seed).
/// All rotations' finalists are pooled into one candidate set so the smoke's
/// count-distinct external rescore selection sees every tie-diverse route;
/// the 2048-finalist cap applies downstream. Counters from later rotations
/// are added onto the rotation-0 accounting so reported work/drops/prunes
/// stay honest, and `complete` requires every rotation to be exhaustive.
#[allow(clippy::too_many_arguments)]
fn bounded_chain_rotations(
    partitions: &[Vec<GenomeAllele>],
    seams: &[Vec<GenomeSeam>],
    viable: &[Vec<bool>],
    incidence: &IncidenceTable,
    model: &ScoreModel,
    profile_memberships: &[Vec<usize>],
    profile_classes: &[Vec<Profile>],
    local_losses: &[BTreeMap<(usize, usize), f64>],
    counts: &mut CountCache<'_>,
    accounting: ChainAccounting,
    force_initial_beam: bool,
    beam_width: usize,
    rss_budget_bytes: Option<u64>,
    inherited_incumbent: Option<f64>,
    tie_rotations: usize,
) -> io::Result<GenomeChain> {
    let original = accounting;
    let rotations = tie_rotations.max(1) as u64;
    let mut incumbent = inherited_incumbent;
    let mut rotation_best_scores = Vec::new();
    let mut merged: Option<GenomeChain> = None;
    // Shared per-locus machinery: registries, weight tables, and suffix
    // max-addition tables are identical for every pass over the same
    // partitions/seams (seed ladder, main pass, rotations), so they are
    // built once here together with the feature universe and shared down.
    let machinery_started = Instant::now();
    let mut registries = Vec::with_capacity(partitions.len());
    registries.push((Vec::new(), BTreeMap::new()));
    for locus in 1..partitions.len() {
        registries.push(build_locus_successors(
            &seams[locus - 1],
            &profile_memberships[locus],
            &viable[locus],
        )?);
    }
    let universe = FeatureUniverse::build(profile_classes, &registries, counts)?;
    let mut weights = Vec::with_capacity(partitions.len());
    weights.push((Vec::new(), Vec::new()));
    for locus in 1..partitions.len() {
        let (classes, _) = &registries[locus];
        weights.push(locus_weight_tables(&profile_classes[locus], classes, counts, model)?);
    }
    let mut max_additions = Vec::with_capacity(partitions.len());
    max_additions.push(FxHashMap::default());
    for locus in 1..partitions.len() {
        let (classes, _) = &registries[locus];
        max_additions
            .push(locus_max_additions(&profile_classes[locus], classes, &universe)?);
    }
    let mut suffix_max_add_total = FxHashMap::<u32, u64>::default();
    for table in &max_additions {
        for (fid, max_add) in table {
            *suffix_max_add_total.entry(*fid).or_default() += max_add;
        }
    }
    let machinery = SharedLocusMachinery {
        registries,
        weights,
        max_additions,
        suffix_max_add_total,
        pair_constants: std::cell::RefCell::new(FxHashMap::default()),
    };
    DpStageProbe::add(0, machinery_started);
    for rotation in 0..rotations {
        let chain = bounded_conflict_chain(
            partitions,
            seams,
            viable,
            incidence,
            model,
            profile_memberships,
            profile_classes,
            local_losses,
            counts,
            original.clone(),
            force_initial_beam,
            beam_width,
            rss_budget_bytes,
            incumbent,
            // seed_mode skips the seed ladder: rotation 0 is the main pass and
            // runs it exactly as the single-pass configuration did; later
            // rotations inherit the realized incumbent instead.
            rotation != 0,
            rotation,
            Some(&universe),
            Some(&machinery),
        )?;
        let best = chain.proposal_losses.first().copied();
        if let Some(score) = best {
            incumbent = Some(match incumbent {
                Some(current) if current.total_cmp(&score).is_lt() => current,
                _ => score,
            });
        }
        rotation_best_scores.push(best);
        if merged.is_none() {
            merged = Some(chain);
            continue;
        }
        let merged = merged.as_mut().expect("rotation 0 completed first");
        {
            let accounting = &mut merged.accounting;
            let pass = &chain.accounting;
            accounting.work = accounting
                .work
                .saturating_add(pass.work.saturating_sub(original.work));
            accounting.dropped_ties = accounting
                .dropped_ties
                .saturating_add(pass.dropped_ties.saturating_sub(original.dropped_ties));
            accounting.dropped_non_ties = accounting.dropped_non_ties.saturating_add(
                pass.dropped_non_ties
                    .saturating_sub(original.dropped_non_ties),
            );
            accounting.tie_backpointers = accounting.tie_backpointers.saturating_add(
                pass.tie_backpointers
                    .saturating_sub(original.tie_backpointers),
            );
            accounting.bound_pruned_states = accounting.bound_pruned_states.saturating_add(
                pass.bound_pruned_states
                    .saturating_sub(original.bound_pruned_states),
            );
            accounting.bound_pruned_edges = accounting.bound_pruned_edges.saturating_add(
                pass.bound_pruned_edges
                    .saturating_sub(original.bound_pruned_edges),
            );
            accounting.state_bytes = accounting.state_bytes.max(pass.state_bytes);
            accounting.sample_count_queries = accounting
                .sample_count_queries
                .max(pass.sample_count_queries);
            accounting.sample_count_cache_resets = accounting
                .sample_count_cache_resets
                .max(pass.sample_count_cache_resets);
            if pass
                .greedy_incumbent
                .is_some_and(|score| {
                    accounting
                        .greedy_incumbent
                        .is_none_or(|current| score.total_cmp(&current).is_lt())
                })
            {
                accounting.greedy_incumbent = pass.greedy_incumbent;
            }
        }
        merged.complete &= chain.complete;
        merged.choices.extend(chain.choices);
        merged.proposal_losses.extend(chain.proposal_losses);
    }
    let mut merged = merged.expect("rotation 0 always runs");
    merged.rotation_best_scores = rotation_best_scores;
    Ok(merged)
}

/// Exact min-plus DP over ordered physical exit pairs when `W=0`; for `W>0`,
/// the approved bounded beam retains exact in-chain conflict eligibility and
/// reports incompleteness and every dropped cutoff tie.
pub fn exact_chain_streaming<F>(
    partitions: &[Vec<GenomeAllele>],
    seams: &[Vec<GenomeSeam>],
    incidence: &IncidenceTable,
    sample: &WeightedBwt,
    model: &ScoreModel,
    profile_cost: ProfileCost,
    force_initial_beam: bool,
    beam_width: usize,
    rss_budget_bytes: Option<u64>,
    inherited_incumbent: Option<f64>,
    tie_rotations: usize,
    mut load_profiles: F,
) -> io::Result<GenomeChain>
where
    F: FnMut(usize) -> io::Result<Vec<Profile>>,
{
    let audit = audit_chain(partitions, seams)?;
    exact_chain_streaming_with_audit(
        partitions,
        seams,
        &audit,
        incidence,
        sample,
        model,
        profile_cost,
        force_initial_beam,
        beam_width,
        rss_budget_bytes,
        inherited_incumbent,
        tie_rotations,
        load_profiles,
    )
}

/// `exact_chain_streaming` with the chain audit supplied by the caller.
/// `audit_chain` walks every seam once per call (26 s on the tract slice:
/// per-seam profile signatures under `exit_memberships`); callers that
/// already computed the audit for their own reporting pass it here instead
/// of paying the sweep a second time inside the chain stage.
#[allow(clippy::too_many_arguments)]
pub fn exact_chain_streaming_with_audit<F>(
    partitions: &[Vec<GenomeAllele>],
    seams: &[Vec<GenomeSeam>],
    audit: &GenomeAudit,
    incidence: &IncidenceTable,
    sample: &WeightedBwt,
    model: &ScoreModel,
    profile_cost: ProfileCost,
    force_initial_beam: bool,
    beam_width: usize,
    rss_budget_bytes: Option<u64>,
    inherited_incumbent: Option<f64>,
    tie_rotations: usize,
    mut load_profiles: F,
) -> io::Result<GenomeChain>
where
    F: FnMut(usize) -> io::Result<Vec<Profile>>,
{
    ensure(
        !partitions.is_empty() && seams.len() + 1 == partitions.len(),
        "invalid genome partition chain",
    )?;
    let mut accounting = ChainAccounting {
        profile_work: profile_cost.mem_queries,
        integrated_windows: profile_cost.integrated_windows,
        exit_classes: audit.exit_classes.clone(),
        conflict_window: audit.conflict_window,
        incidence_hashes: incidence.len(),
        incidence_collisions: incidence.collisions,
        ..ChainAccounting::default()
    };
    if accounting.profile_work > MAX_PROFILE_WORK {
        return Ok(GenomeChain {
            complete: false,
            stop_reason: "budget: experimental_genome_profile_work".into(),
            choices: Vec::new(),
            proposal_losses: Vec::new(),
            rotation_best_scores: Vec::new(),
            accounting,
        });
    }
    let mut counts = CountCache {
        sample,
        values: (0..256).map(|_| HashMap::new()).collect(),
        entries: 0,
        queries: 0,
        resets: 0,
    };
    let target_source = partitions[0][0].traversal.entry().source;
    let mut profile_memberships = Vec::with_capacity(partitions.len());
    let mut profile_classes = Vec::with_capacity(partitions.len());
    let mut local_losses = Vec::with_capacity(partitions.len());
    // Hoisted path decision (the same expressions the post-loop branches
    // used, computed once): whether the bounded rotation machinery will run
    // this chain. IMPG_DEV_LAZY_CLASS_EVIDENCE (dev-only) skips the per-locus
    // all-pairs class evidence — prepare_profile plus every class-pair loss
    // — which is reporting-only data on the bounded path (the DP search
    // consumes profile classes, not local losses). It is legal only when the
    // bounded path is guaranteed, because the exact path searches over
    // local_losses. Dev runs report an empty local_pair_evidence array.
    let viable = completion_viability(partitions, seams);
    let first_alleles = viable[0].iter().filter(|&&is_viable| is_viable).count();
    let first_states = first_alleles.saturating_mul(first_alleles.saturating_add(1)) / 2;
    let first_state_bytes =
        first_states.saturating_mul(128usize.saturating_add(std::mem::size_of::<PhysicalNode>()))
            as u64;
    let bounded_path =
        audit.conflict_window > 0 || force_initial_beam || first_state_bytes > MAX_STATE_BYTES;
    let dev_lazy_evidence =
        bounded_path && dev_flag("IMPG_DEV_LAZY_CLASS_EVIDENCE");
    for (locus, alleles) in partitions.iter().enumerate() {
        let locus_setup = std::env::var("IMPG_DB").is_ok()
            .then(std::time::Instant::now);
        let load_started = Instant::now();
        let physical_profiles = load_profiles(locus)?;
        DpStageProbe::add(11, load_started);
        let evidence_started = Instant::now();
        if let Some(started) = locus_setup {
            eprintln!(
                "SETUP locus {} load_ms {:.1}",
                locus,
                started.elapsed().as_secs_f64() * 1e3
            );
        }
        ensure(
            physical_profiles.len() == alleles.len(),
            "streamed locus profile count mismatch",
        )?;
        let (profiles, membership) = class_for_profiles(&physical_profiles);
        if dev_lazy_evidence {
            profile_memberships.push(membership);
            profile_classes.push(profiles);
            local_losses.push(BTreeMap::new());
            DpStageProbe::add(10, evidence_started);
            continue;
        }
        let prepared = profiles
            .iter()
            .map(|profile| {
                prepare_profile(profile, LocalClass::Owned(locus), incidence, &mut counts)
            })
            .collect::<io::Result<Vec<_>>>()?;
        let pairs = (0..profiles.len())
            .flat_map(|first| (first..profiles.len()).map(move |second| (first, second)))
            .collect::<Vec<_>>();
        accounting.work = accounting.work.saturating_add(pairs.len() as u64);
        let losses = pairs
            .into_par_iter()
            .map(|pair| {
                score_prepared_pair(&prepared[pair.0], &prepared[pair.1], model)
                    .map(|loss| (pair, loss))
            })
            .collect::<io::Result<BTreeMap<_, _>>>()?;
        if let Some(started) = locus_setup {
            eprintln!(
                "SETUP locus {} classes {} pairs {} score_ms {:.1}",
                locus,
                profiles.len(),
                losses.len(),
                started.elapsed().as_secs_f64() * 1e3
            );
        }
        let (&best_pair, &best_loss) = losses
            .iter()
            .min_by(|a, b| a.1.total_cmp(b.1))
            .ok_or_else(|| invalid("empty local pair loss table"))?;
        let mut representatives = vec![usize::MAX; profiles.len()];
        for (physical, &class) in membership.iter().enumerate() {
            representatives[class] = representatives[class].min(physical);
        }
        let native_class = alleles
            .iter()
            .enumerate()
            .find(|(_, allele)| {
                allele.traversal.segments.len() == 1
                    && allele.traversal.entry().source == target_source
                    && !allele.traversal.entry().reverse
            })
            .map(|(physical, _)| membership[physical]);
        let native_pair = native_class.map(|class| (class, class));
        let native_loss = native_pair.map(|pair| losses[&pair]);
        let mut ranked_classes = (0..profiles.len())
            .map(|class| (class, losses[&(class, class)]))
            .collect::<Vec<_>>();
        ranked_classes.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
        ranked_classes.truncate(8);
        let (best_owned_observed, best_owned_predicted) =
            prepared_pair_totals(&prepared[best_pair.0], &prepared[best_pair.1])?;
        let native_totals = native_pair
            .map(|pair| prepared_pair_totals(&prepared[pair.0], &prepared[pair.1]))
            .transpose()?;
        let best_ranges = [
            alleles[representatives[best_pair.0]].traversal.entry(),
            alleles[representatives[best_pair.1]].traversal.entry(),
        ];
        accounting.local_pair_evidence.push(LocalPairEvidence {
            locus,
            best_sources: [best_ranges[0].source, best_ranges[1].source],
            best_ranges: [
                (
                    best_ranges[0].start,
                    best_ranges[0].end,
                    best_ranges[0].reverse,
                ),
                (
                    best_ranges[1].start,
                    best_ranges[1].end,
                    best_ranges[1].reverse,
                ),
            ],
            best_loss,
            native_loss,
            native_minus_best: native_loss.map(|loss| loss - best_loss),
            best_owned_observed,
            best_owned_predicted,
            native_owned_observed: native_totals.map(|totals| totals.0),
            native_owned_predicted: native_totals.map(|totals| totals.1),
            top_class_representatives: ranked_classes
                .iter()
                .map(|(class, _)| representatives[*class])
                .collect(),
            top_homozygous_losses: ranked_classes.iter().map(|(_, loss)| *loss).collect(),
            selected_loss: None,
            selected_minus_best: None,
        });
        profile_memberships.push(membership);
        profile_classes.push(profiles);
        local_losses.push(losses);
        DpStageProbe::add(10, evidence_started);
    }
    if bounded_path {
        return bounded_chain_rotations(
            partitions,
            seams,
            &viable,
            incidence,
            model,
            &profile_memberships,
            &profile_classes,
            &local_losses,
            &mut counts,
            accounting,
            force_initial_beam,
            beam_width,
            rss_budget_bytes,
            inherited_incumbent,
            tie_rotations,
        );
    }
    let local_loss = |locus: usize, pair: [usize; 2]| {
        let first = profile_memberships[locus][pair[0]];
        let second = profile_memberships[locus][pair[1]];
        let key = if first <= second {
            (first, second)
        } else {
            (second, first)
        };
        local_losses[locus][&key]
    };
    let mut nodes = Vec::<PhysicalNode>::new();
    let mut layer = BTreeMap::<[usize; 2], usize>::new();
    for first in 0..partitions[0].len() {
        if !viable[0][first] {
            continue;
        }
        for second in first..partitions[0].len() {
            if !viable[0][second] {
                continue;
            }
            let pair = [first, second];
            let node = nodes.len();
            nodes.push(PhysicalNode {
                score: local_loss(0, pair),
                choices: pair,
                predecessors: Vec::new(),
            });
            layer.insert(pair, node);
        }
    }
    accounting.occupancy.push(layer.len());
    accounting.state_bytes = physical_state_bytes(layer.len(), &nodes);
    if accounting.state_bytes > MAX_STATE_BYTES {
        return Ok(GenomeChain {
            complete: false,
            stop_reason: "budget: exact_physical_dp_state".into(),
            choices: Vec::new(),
            proposal_losses: Vec::new(),
            rotation_best_scores: Vec::new(),
            accounting,
        });
    }
    for locus in 1..partitions.len() {
        enforce_chain_rss(rss_budget_bytes, "exact_chain_locus_start")?;
        let mut successors: BTreeMap<usize, Vec<(usize, Vec<(FeatureKey, u64)>)>> = BTreeMap::new();
        for seam in &seams[locus - 1] {
            if viable[locus][seam.right] {
                successors
                    .entry(seam.left)
                    .or_default()
                    .push((seam.right, profile_signature(&seam.profile)));
            }
        }
        let mut boundary_cache = BTreeMap::new();
        let mut next = BTreeMap::<[usize; 2], usize>::new();
        for (previous_pair, &previous_node) in &layer {
            let Some(first_successors) = successors.get(&previous_pair[0]) else {
                continue;
            };
            let Some(second_successors) = successors.get(&previous_pair[1]) else {
                continue;
            };
            for (first, first_signature) in first_successors {
                for (second, second_signature) in second_successors {
                    accounting.work = accounting.work.saturating_add(1);
                    let pair = [*first, *second];
                    let mut cache_key = (first_signature.clone(), second_signature.clone());
                    if cache_key.1 < cache_key.0 {
                        cache_key = (cache_key.1, cache_key.0);
                    }
                    let boundary_loss = if let Some(loss) = boundary_cache.get(&cache_key) {
                        *loss
                    } else {
                        let first_profile: Profile = first_signature.iter().cloned().collect();
                        let second_profile: Profile = second_signature.iter().cloned().collect();
                        let loss = score_local_pair(
                            &first_profile,
                            &second_profile,
                            LocalClass::Boundary(locus - 1),
                            incidence,
                            &mut counts,
                            model,
                        )?;
                        boundary_cache.insert(cache_key, loss);
                        loss
                    };
                    let total =
                        nodes[previous_node].score + local_loss(locus, pair) + boundary_loss;
                    if let Some(&node) = next.get(&pair) {
                        match total.total_cmp(&nodes[node].score) {
                            std::cmp::Ordering::Less => {
                                nodes[node].score = total;
                                nodes[node].predecessors.clear();
                                nodes[node].predecessors.push(previous_node);
                            }
                            std::cmp::Ordering::Equal => {
                                nodes[node].predecessors.push(previous_node);
                                accounting.tie_backpointers += 1;
                            }
                            std::cmp::Ordering::Greater => {}
                        }
                    } else {
                        let node = nodes.len();
                        nodes.push(PhysicalNode {
                            score: total,
                            choices: pair,
                            predecessors: vec![previous_node],
                        });
                        next.insert(pair, node);
                    }
                }
            }
        }
        accounting.occupancy.push(next.len());
        accounting.state_bytes = accounting
            .state_bytes
            .max(physical_state_bytes(next.len(), &nodes));
        if accounting.state_bytes > MAX_STATE_BYTES {
            accounting.sample_count_queries = counts.queries;
            accounting.sample_count_cache_resets = counts.resets;
            accounting.bounded_blocks = 1;
            return Ok(GenomeChain {
                complete: false,
                stop_reason: "budget: exact_physical_dp_state".into(),
                choices: Vec::new(),
                proposal_losses: Vec::new(),
                rotation_best_scores: Vec::new(),
                accounting,
            });
        }
        if next.is_empty() {
            accounting.sample_count_queries = counts.queries;
            accounting.sample_count_cache_resets = counts.resets;
            return Ok(GenomeChain {
                complete: true,
                stop_reason: "no public legal chained state".into(),
                choices: Vec::new(),
                proposal_losses: Vec::new(),
                rotation_best_scores: Vec::new(),
                accounting,
            });
        }
        layer = next;
    }
    let best = layer
        .values()
        .map(|&node| nodes[node].score)
        .min_by(f64::total_cmp)
        .ok_or_else(|| invalid("empty exact DP terminal layer"))?;
    let mut choices = Vec::new();
    let mut truncated = false;
    for node in layer
        .values()
        .copied()
        .filter(|&node| nodes[node].score.total_cmp(&best).is_eq())
    {
        backtrack_physical(&nodes, node, &mut Vec::new(), &mut choices, &mut truncated);
    }
    accounting.sample_count_queries = counts.queries;
    accounting.sample_count_cache_resets = counts.resets;
    accounting.exact_blocks = 1;
    if truncated {
        accounting.bounded_blocks = 1;
    }
    Ok(GenomeChain {
        complete: !truncated,
        stop_reason: if truncated {
            "budget: exact_physical_tie_enumeration".into()
        } else {
            "exact_min_plus_dp".into()
        },
        proposal_losses: vec![best; choices.len()],
        choices,
        accounting,
        rotation_best_scores: vec![Some(best)],
    })
}

pub fn exact_chain(
    partitions: &[Vec<GenomeAllele>],
    seams: &[Vec<GenomeSeam>],
    incidence: &IncidenceTable,
    sample: &WeightedBwt,
    model: &ScoreModel,
    profile_cost: ProfileCost,
) -> io::Result<GenomeChain> {
    exact_chain_streaming(
        partitions,
        seams,
        incidence,
        sample,
        model,
        profile_cost,
        false,
        CONFLICT_BEAM_WIDTH,
        None,
        None,
        1,
        |locus| {
            Ok(partitions[locus]
                .iter()
                .map(|allele| allele.profile.clone())
                .collect())
        },
    )
}

fn write_run(path: &Path, profile: Profile) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for entry in profile {
        serde_json::to_writer(&mut writer, &entry).map_err(io::Error::other)?;
        writer.write_all(b"\n")?;
    }
    writer.flush()
}

fn next_run(reader: &mut BufReader<File>) -> io::Result<Option<(FeatureKey, u64)>> {
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        return Ok(None);
    }
    serde_json::from_str(&line)
        .map(Some)
        .map_err(io::Error::other)
}

fn merge_runs(inputs: &[PathBuf], output: &Path) -> io::Result<()> {
    let mut readers = inputs
        .iter()
        .map(|path| File::open(path).map(BufReader::new))
        .collect::<io::Result<Vec<_>>>()?;
    let mut heads = readers
        .iter_mut()
        .map(next_run)
        .collect::<io::Result<Vec<_>>>()?;
    let mut writer = BufWriter::new(File::create(output)?);
    while let Some(key) = heads
        .iter()
        .filter_map(|entry| entry.as_ref().map(|(key, _)| key))
        .min()
        .cloned()
    {
        let mut total = 0u64;
        for index in 0..heads.len() {
            if heads[index]
                .as_ref()
                .is_some_and(|(candidate, _)| *candidate == key)
            {
                total = total
                    .checked_add(heads[index].as_ref().unwrap().1)
                    .ok_or_else(|| invalid("external profile count overflow"))?;
                heads[index] = next_run(&mut readers[index])?;
            }
        }
        serde_json::to_writer(&mut writer, &(key, total)).map_err(io::Error::other)?;
        writer.write_all(b"\n")?;
    }
    writer.flush()
}

/// Exact complete-copy full-subwalk rescore using bounded profile chunks and a
/// bounded-fan-in external merge. Candidate-zero terms are omitted because their
/// background-relative loss is exactly zero.
pub fn external_rescore(
    panel: &SyngIndex,
    sequences: [&[u8]; 2],
    sample: &WeightedBwt,
    model: &ScoreModel,
    directory: &Path,
) -> io::Result<(f64, ProfileCost, usize)> {
    external_rescore_components(panel, &[sequences], sample, model, directory)
}

/// Jointly rescore chromosome-reset copies. Every chromosome emits independent
/// bounded sorted runs, preventing fictitious cross-chromosome read starts; all
/// runs are then merged before one diploid Poisson loss.
pub fn external_rescore_components(
    panel: &SyngIndex,
    components: &[[&[u8]; 2]],
    sample: &WeightedBwt,
    model: &ScoreModel,
    directory: &Path,
) -> io::Result<(f64, ProfileCost, usize)> {
    fs::create_dir_all(directory)?;
    let read_length = usize::try_from(model.read_length)
        .map_err(|_| invalid("read length exceeds address space"))?;
    let mut runs = Vec::new();
    let mut cost = ProfileCost::default();
    for (component, sequences) in components.iter().enumerate() {
        for (copy, &sequence) in sequences.iter().enumerate() {
            if sequence.len() < read_length {
                continue;
            }
            let starts = sequence.len() - read_length + 1;
            let mut lo = 0;
            while lo < starts {
                let mut hi = (lo + 512).min(starts);
                let (profile, next) = loop {
                    match profile_event_runs(panel, sequence, read_length, lo, hi, MAX_FEATURES) {
                        Ok(value) => break value,
                        Err(error)
                            if error.to_string() == "budget: maximal_mem_subwalk_features"
                                && hi > lo + 1 =>
                        {
                            hi = lo + (hi - lo) / 2;
                        }
                        Err(error) => return Err(error),
                    }
                };
                add_cost(&mut cost, next)?;
                ensure(
                    cost.mem_queries <= MAX_PROFILE_WORK,
                    "budget: experimental_genome_profile_work",
                )?;
                let path = directory.join(format!(
                    "component-{component}-copy-{copy}-{:08}.jsonl",
                    runs.len()
                ));
                write_run(&path, profile)?;
                runs.push(path);
                lo = hi;
            }
        }
    }
    let initial_runs = runs.len();
    let mut round = 0usize;
    while runs.len() > 1 {
        let mut next = Vec::new();
        for (group, chunk) in runs.chunks(16).enumerate() {
            let output = directory.join(format!("merge-{round}-{group:08}.jsonl"));
            merge_runs(chunk, &output)?;
            next.push(output);
        }
        for path in runs {
            fs::remove_file(path)?;
        }
        runs = next;
        round += 1;
    }
    let mut loss = 0.0;
    if let Some(path) = runs.pop() {
        let mut reader = BufReader::new(File::open(&path)?);
        while let Some((key, predicted)) = next_run(&mut reader)? {
            loss += model.loss(predicted, sample.count(&key)?)?;
        }
        fs::remove_file(path)?;
    }
    ensure(loss.is_finite(), "nonfinite external complete rescore")?;
    Ok((loss, cost, initial_runs))
}

pub fn coalesced_route(ranges: &[SourceRange]) -> io::Result<routes::Route> {
    ensure(!ranges.is_empty(), "empty genome route")?;
    ensure(
        spans_feasible(ranges),
        "overlapping genome route source spans",
    )?;
    let mut segments: Vec<routes::Segment> = Vec::new();
    for range in ranges {
        if range.start == range.end {
            continue;
        }
        if let Some(previous) = segments.last_mut() {
            let contiguous = previous.source == range.source
                && previous.reverse == range.reverse
                && ((!range.reverse && previous.end == range.start)
                    || (range.reverse && previous.start == range.end));
            if contiguous {
                previous.start = previous.start.min(range.start);
                previous.end = previous.end.max(range.end);
                continue;
            }
        }
        segments.push(routes::Segment {
            source: range.source,
            start: range.start,
            end: range.end,
            reverse: range.reverse,
        });
    }
    ensure(segments.len() <= 64, "budget: partition_segments")?;
    Ok(routes::Route { segments })
}

#[cfg(test)]
mod tests {
    use super::*;
    use impg::{sample_mem_bwt::WeightedBwt, syng::SyncmerParams};
    use tempfile::tempdir;

    fn dna(length: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..length)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x as usize) & 3]
            })
            .collect()
    }

    #[test]
    fn binary_records_sidecar_serves_rebuilds_tail_syncs_and_corruption_fallback() {
        let directory = tempdir().unwrap();
        let path = directory.path().join("profiles.jsonl");
        let cost = ProfileCost {
            mem_queries: 3,
            integrated_windows: 5,
        };
        let profile_a: Profile = [(vec![1u64], 2u64), (vec![3], 4)].into_iter().collect();
        let profile_b: Profile = [(vec![7u64], 1u64)].into_iter().collect();
        let profile_c: Profile = [(vec![9u64], 8u64), (vec![11], 2)].into_iter().collect();
        // First run: create the archive with two records.
        {
            let mut cache = ProfileCache::open(&path).unwrap();
            let (stored, stored_cost, reused) = cache
                .get_or_insert_with("allele:a", || Ok((profile_a.clone(), cost)))
                .unwrap();
            assert!(!reused);
            assert_eq!(stored, profile_a);
            assert_eq!(stored_cost, cost);
            cache
                .get_or_insert_with("allele:b", || Ok((profile_b.clone(), cost)))
                .unwrap();
        }
        // Second open: the sidecar is built from scratch and serves reads.
        {
            let cache = ProfileCache::open(&path).unwrap();
            let (used, covered, status, _) = cache.records_sidecar_report();
            assert!(used, "sidecar must serve after a clean rebuild");
            assert_eq!(status, "rebuilt_absent_2_records");
            assert!(covered > 0);
            let (stored, stored_cost, _) = cache.get_if_cached("allele:a").unwrap().unwrap();
            assert_eq!(stored, profile_a);
            assert_eq!(stored_cost, cost);
            assert_eq!(cache.get_if_cached("allele:b").unwrap().unwrap().0, profile_b);
            assert!(cache.get_if_cached("allele:missing").unwrap().is_none());
        }
        // Third run: an append makes the sidecar stale; the next open tail-syncs.
        {
            let mut cache = ProfileCache::open(&path).unwrap();
            cache
                .get_or_insert_with("allele:c", || Ok((profile_c.clone(), cost)))
                .unwrap();
            // In-process, the appended record is beyond the sidecar's
            // coverage and is served through the JSONL path.
            let (stored, _, reused) = cache
                .get_or_insert_with("allele:c", || panic!("appended record requeried"))
                .unwrap();
            assert!(reused);
            assert_eq!(stored, profile_c);
        }
        {
            let cache = ProfileCache::open(&path).unwrap();
            let (used, _, status, _) = cache.records_sidecar_report();
            assert!(used);
            assert_eq!(status, "tail_synced_1_records");
            for (key, profile) in [
                ("allele:a", &profile_a),
                ("allele:b", &profile_b),
                ("allele:c", &profile_c),
            ] {
                assert_eq!(cache.get_if_cached(key).unwrap().unwrap().0, *profile);
            }
        }
        // Fourth open with no changes: the sidecar loads without writes.
        {
            let cache = ProfileCache::open(&path).unwrap();
            let (_, _, status, _) = cache.records_sidecar_report();
            assert_eq!(status, "loaded_3_records");
        }
        // A corrupt sidecar falls back to a verified rebuild, and reads stay
        // correct.
        {
            let sidecar = records_sidecar_path(&path);
            let corrupted = b"garbage-not-a-sidecar";
            fs::write(&sidecar, corrupted).unwrap();
            let cache = ProfileCache::open(&path).unwrap();
            let (used, _, status, _) = cache.records_sidecar_report();
            assert!(used, "corruption must rebuild, not disable");
            assert_eq!(status, "rebuilt_corrupt_3_records");
            assert_eq!(cache.get_if_cached("allele:a").unwrap().unwrap().0, profile_a);
            assert_eq!(cache.get_if_cached("allele:c").unwrap().unwrap().0, profile_c);
        }
        // A key mismatch inside the sidecar fails closed on read.
        {
            let cache = ProfileCache::open(&path).unwrap();
            let sidecar = records_sidecar_path(&path);
            let bytes = fs::read(&sidecar).unwrap();
            // Flip the first record's stored key bytes: the record payloads
            // start at offset 0, key length first.
            let mut corrupted = bytes;
            let key_len = u32::from_le_bytes(corrupted[0..4].try_into().unwrap()) as usize;
            for byte in corrupted[4..4 + key_len].iter_mut() {
                *byte ^= 0xff;
            }
            fs::write(&sidecar, corrupted).unwrap();
            let error = cache.get_if_cached("allele:a").unwrap_err();
            assert!(error.to_string().contains("records sidecar key mismatch"));
        }
    }

    #[test]
    fn event_compressed_profiles_equal_every_window_full_subwalks_and_cache_once() {
        let sequence = dna(2200, 731);
        let panel = SyngIndex::build(
            SyncmerParams::default(),
            [("s".into(), sequence.clone())].into_iter(),
        );
        let (event, cost) = profile_event_interior(&panel, &sequence, 150, MAX_FEATURES).unwrap();
        let exhaustive =
            super::super::partition::profile_interior(&panel, &sequence, 150, MAX_FEATURES)
                .unwrap();
        assert_eq!(event, exhaustive);
        assert_eq!(cost.integrated_windows, 2051);
        assert!(cost.mem_queries < cost.integrated_windows);

        let directory = tempdir().unwrap();
        let path = directory.path().join("profiles.jsonl");
        let mut cache = ProfileCache::open(&path).unwrap();
        let (_, first_cost, reused) = cache
            .get_or_insert_with("allele:0", || Ok((event.clone(), cost)))
            .unwrap();
        assert!(!reused);
        let (_, second_cost, reused) = cache
            .get_or_insert_with("allele:0", || panic!("cached profile was requeried"))
            .unwrap();
        assert!(reused);
        assert_eq!(first_cost, second_cost);
        drop(cache);
        assert_eq!(ProfileCache::open(&path).unwrap().entries.len(), 1);

        let sample = WeightedBwt::build(&event).unwrap();
        let model = ScoreModel {
            read_length: 150,
            histogram: 1,
            denominator: 150.0,
            depth: 150.0,
            background: 0.1,
        };
        let expected = event
            .iter()
            .try_fold(0.0, |loss, (key, &count)| {
                Ok::<_, io::Error>(loss + model.loss(count * 2, sample.count(key)?)?)
            })
            .unwrap();
        let (external, _, runs) = external_rescore(
            &panel,
            [&sequence, &sequence],
            &sample,
            &model,
            &directory.path().join("runs"),
        )
        .unwrap();
        assert!(runs > 2);
        assert!((external - expected).abs() < 1e-9);

        let expected_components = event
            .iter()
            .try_fold(0.0, |loss, (key, &count)| {
                Ok::<_, io::Error>(loss + model.loss(count * 4, sample.count(key)?)?)
            })
            .unwrap();
        let (component_external, _, component_runs) = external_rescore_components(
            &panel,
            &[[&sequence, &sequence], [&sequence, &sequence]],
            &sample,
            &model,
            &directory.path().join("component-runs"),
        )
        .unwrap();
        assert!(component_runs > runs);
        assert!((component_external - expected_components).abs() < 1e-9);
    }

    #[test]
    fn deletion_marker_exits_only_by_same_source_continuity() {
        let range = |partition, occurrence, source, start, end| SourceRange {
            partition,
            occurrence,
            source,
            start,
            end,
            reverse: false,
        };
        let axis = vec![
            AxisInterval {
                component: "target".into(),
                start: 0,
                end: 10,
                group: "a".into(),
                reference_occurrence: 0,
                reference_strand: "+".into(),
                orientations: BTreeMap::new(),
            },
            AxisInterval {
                component: "target".into(),
                start: 10,
                end: 20,
                group: "b".into(),
                reference_occurrence: 1,
                reference_strand: "+".into(),
                orientations: BTreeMap::new(),
            },
        ];
        let partitions = vec![
            vec![SpanningTraversal::single(range(0, 0, 7, 10, 10))],
            vec![
                SpanningTraversal::single(range(1, 1, 7, 10, 20)),
                SpanningTraversal::single(range(1, 2, 8, 100, 110)),
            ],
        ];
        let links =
            port_word_seams_with(&axis, &partitions, |_, _| Ok(Some(b"MATCH".to_vec()))).unwrap();
        assert_eq!(links, vec![vec![(0, 0)]]);
    }

    #[test]
    fn same_source_reverse_order_split_cut_is_rejected_before_admission() {
        let range = |occurrence, source| SourceRange {
            partition: 0,
            occurrence,
            source,
            start: 100,
            end: 200,
            reverse: false,
        };
        let left = range(1, 7);
        let right = range(2, 7);
        assert!(legal_split_cut_order(&left, &right, 140, 160));
        assert!(!legal_split_cut_order(&left, &right, 160, 140));
        assert!(!legal_split_cut_order(&left, &right, 150, 150));
        assert!(legal_split_cut_order(&left, &range(3, 8), 160, 140));
    }

    #[test]
    fn source_path_completion_fills_missing_axis_loci_without_coordinate_projection() {
        let range = |partition, occurrence, source, start, end| SourceRange {
            partition,
            occurrence,
            source,
            start,
            end,
            reverse: false,
        };
        let mut partitions = vec![
            vec![range(0, 0, 7, 0, 10), range(0, 1, 8, 100, 110)],
            vec![range(1, 2, 7, 10, 20)],
            vec![range(2, 3, 7, 20, 30)],
            vec![range(3, 4, 7, 30, 40), range(3, 5, 8, 120, 130)],
        ];
        let stats =
            complete_forward_source_paths(&mut partitions, &BTreeSet::from([7usize, 8usize]))
                .unwrap();
        assert_eq!(stats.bridge_candidates, 1);
        assert_eq!(stats.deletion_candidates, 2);
        assert!(partitions[0]
            .iter()
            .any(|range| range.source == 8 && (range.start, range.end) == (100, 120)));
        for locus in [1, 2] {
            assert!(partitions[locus]
                .iter()
                .any(|range| range.source == 8 && (range.start, range.end) == (120, 120)));
        }
    }

    #[test]
    fn chr_iii_public_sk1_gaps_receive_bridge_and_deletion_candidates() {
        let range = |partition, occurrence, start, end| SourceRange {
            partition,
            occurrence,
            source: 9602,
            start,
            end,
            reverse: false,
        };
        let mut partitions = vec![Vec::new(); 38];
        partitions[15].push(range(15, 1, 128754, 138845));
        partitions[20].push(range(20, 2, 156413, 166355));
        partitions[21].push(range(21, 3, 172876, 178652));
        partitions[24].push(range(24, 4, 178652, 190177));
        let stats =
            complete_forward_source_paths(&mut partitions, &BTreeSet::from([9602usize])).unwrap();
        for locus in [16, 18, 19, 22, 23] {
            assert!(partitions[locus].iter().any(|range| range.source == 9602));
        }
        assert!(partitions[15]
            .iter()
            .any(|range| (range.start, range.end) == (128754, 156413)));
        assert!(partitions[20]
            .iter()
            .any(|range| (range.start, range.end) == (156413, 172876)));
        assert_eq!(stats.bridge_candidates, 2);
        assert_eq!(stats.deletion_candidates, 6);
    }

    #[test]
    fn native_endpoint_admission_rejects_full_donor_routes_before_scoring() {
        let range = |partition, occurrence, source, start, end, reverse| SourceRange {
            partition,
            occurrence,
            source,
            start,
            end,
            reverse,
        };
        let mut partitions = vec![
            vec![
                range(0, 0, 7, 0, 10, false),
                range(0, 1, 8, 0, 10, false),
                range(0, 2, 7, 0, 10, true),
            ]
            .into_iter()
            .map(SpanningTraversal::single)
            .collect(),
            vec![range(1, 3, 7, 10, 20, false), range(1, 4, 8, 10, 19, false)]
                .into_iter()
                .map(SpanningTraversal::single)
                .collect(),
        ];
        retain_native_endpoint_candidates(&mut partitions, 7, 20).unwrap();
        assert_eq!(partitions[0].len(), 1);
        assert_eq!(partitions[0][0].entry().source, 7);
        assert_eq!(partitions[1].len(), 1);
        assert_eq!(partitions[1][0].entry().source, 7);
    }

    #[test]
    fn port_word_equal_cross_source_seams_form_native_anchored_chain() {
        let range = |partition, occurrence, source, start, end| SourceRange {
            partition,
            occurrence,
            source,
            start,
            end,
            reverse: false,
        };
        let axis = vec![
            AxisInterval {
                component: "target".into(),
                start: 0,
                end: 10,
                group: "a".into(),
                reference_occurrence: 0,
                reference_strand: "+".into(),
                orientations: BTreeMap::new(),
            },
            AxisInterval {
                component: "target".into(),
                start: 10,
                end: 20,
                group: "b".into(),
                reference_occurrence: 1,
                reference_strand: "+".into(),
                orientations: BTreeMap::new(),
            },
            AxisInterval {
                component: "target".into(),
                start: 20,
                end: 30,
                group: "c".into(),
                reference_occurrence: 2,
                reference_strand: "+".into(),
                orientations: BTreeMap::new(),
            },
        ];
        let mut partitions = vec![
            vec![range(0, 0, 7, 0, 10)]
                .into_iter()
                .map(SpanningTraversal::single)
                .collect(),
            vec![range(1, 1, 7, 10, 20), range(1, 2, 8, 100, 110)]
                .into_iter()
                .map(SpanningTraversal::single)
                .collect(),
            vec![range(2, 3, 7, 20, 30)]
                .into_iter()
                .map(SpanningTraversal::single)
                .collect(),
        ];
        retain_native_endpoint_candidates(&mut partitions, 7, 30).unwrap();
        let links = port_word_seams_with(&axis, &partitions, |range, exit| {
            let word = match (range.source, range.start, exit) {
                (7, 0, true) | (8, 100, false) => b"LEFT".to_vec(),
                (8, 100, true) | (7, 20, false) => b"RIGHT".to_vec(),
                (7, 10, false) => b"OTHER-ENTRY".to_vec(),
                (7, 10, true) => b"OTHER-EXIT".to_vec(),
                _ => return Ok(None),
            };
            Ok(Some(word))
        })
        .unwrap();
        assert_eq!(links, vec![vec![(0, 0), (0, 1)], vec![(0, 0), (1, 0)]]);
        assert_eq!(partitions[1][1].entry().source, 8);
    }

    #[test]
    fn split_refinement_recovers_supported_two_segment_middle_allele() {
        let allele =
            |partition, occurrence, segments: Vec<SourceRange>, profile: Profile| GenomeAllele {
                traversal: SpanningTraversal {
                    partition,
                    identity: format!("fixture-{occurrence}"),
                    segments,
                },
                profile,
            };
        let segment = |partition, occurrence, source, start, end| SourceRange {
            partition,
            occurrence,
            source,
            start,
            end,
            reverse: false,
        };
        let partitions = vec![
            vec![allele(0, 0, vec![segment(0, 0, 7, 0, 10)], Profile::new())],
            vec![
                allele(
                    1,
                    1,
                    vec![segment(1, 1, 7, 10, 20)],
                    BTreeMap::from([(vec![6], 1)]),
                ),
                allele(
                    1,
                    2,
                    vec![segment(1, 2, 7, 10, 15), segment(1, 3, 8, 100, 105)],
                    BTreeMap::from([(vec![4], 1)]),
                ),
            ],
            vec![allele(2, 4, vec![segment(2, 4, 7, 20, 30)], Profile::new())],
        ];
        let seams = vec![
            vec![
                GenomeSeam {
                    left: 0,
                    right: 0,
                    profile: Arc::new(Profile::new()),
                },
                GenomeSeam {
                    left: 0,
                    right: 1,
                    profile: Arc::new(Profile::new()),
                },
            ],
            vec![
                GenomeSeam {
                    left: 0,
                    right: 0,
                    profile: Arc::new(Profile::new()),
                },
                GenomeSeam {
                    left: 1,
                    right: 0,
                    profile: Arc::new(Profile::new()),
                },
            ],
        ];
        let sample = WeightedBwt::build(&BTreeMap::from([(vec![4], 2)])).unwrap();
        let model = ScoreModel {
            read_length: 150,
            histogram: 1,
            denominator: 150.0,
            depth: 150.0,
            background: 0.1,
        };
        let incidence = incidence_table(&partitions, &seams);
        let result = exact_chain(
            &partitions,
            &seams,
            &incidence,
            &sample,
            &model,
            ProfileCost::default(),
        )
        .unwrap();
        assert!(result.choices[0][0][1] == 1 || result.choices[0][1][1] == 1);
        assert_eq!(partitions[1][1].traversal.segments.len(), 2);
    }

    #[test]
    fn bounded_conflict_beam_reports_dropped_ties_and_rejects_overlapping_history() {
        let sample = WeightedBwt::build(&BTreeMap::from([(vec![4], 1)])).unwrap();
        let first = (0..17)
            .map(|occurrence| GenomeAllele {
                traversal: SpanningTraversal::single(SourceRange {
                    partition: 0,
                    occurrence,
                    source: 0,
                    start: 0,
                    end: 10,
                    reverse: false,
                }),
                profile: Profile::new(),
            })
            .collect::<Vec<_>>();
        let second = vec![GenomeAllele {
            traversal: SpanningTraversal::single(SourceRange {
                partition: 1,
                occurrence: 17,
                source: 0,
                start: 5,
                end: 15,
                reverse: false,
            }),
            profile: Profile::new(),
        }];
        let seams = vec![(0..17)
            .map(|left| GenomeSeam {
                left,
                right: 0,
                profile: Arc::new(Profile::new()),
            })
            .collect::<Vec<_>>()];
        let partitions = vec![first, second];
        let incidence = incidence_table(&partitions, &seams);
        let model = ScoreModel {
            read_length: 150,
            histogram: 1,
            denominator: 150.0,
            depth: 150.0,
            background: 0.1,
        };
        let result = exact_chain(
            &partitions,
            &seams,
            &incidence,
            &sample,
            &model,
            ProfileCost::default(),
        )
        .unwrap();
        assert!(!result.complete);
        assert_eq!(result.accounting.conflict_window, 1);
        assert_eq!(result.accounting.beam_width, CONFLICT_BEAM_WIDTH);
        assert_eq!(result.accounting.occupancy, vec![153]);
        assert_eq!(result.accounting.dropped_ties, 0);
        assert_eq!(result.accounting.dropped_non_ties, 0);
        assert!(result.choices.is_empty());
        assert_eq!(
            result.stop_reason,
            "bounded_conflict_beam_no_retained_legal_state"
        );
    }

    #[test]
    fn beam_pruning_counts_cutoff_ties() {
        let mut accounting = ChainAccounting::default();
        let candidates = (0..153)
            .map(|index| {
                (
                    BeamKey {
                        history: vec![[index, index]],
                        profile_id: index as u64,
                    },
                    BeamCandidate {
                        score: 0.0,
                        choices: [index, index],
                        predecessors: Vec::new(),
                        cumulative_profile: Arc::new(PersistentProfile {
                            id: index as u64,
                            parent: None,
                            components: [
                                Arc::new(Vec::new()),
                                Arc::new(Vec::new()),
                                Arc::new(Vec::new()),
                                Arc::new(Vec::new()),
                            ],
                            fingerprint: 0,
                            entries: 0,
                            score: 0.0,
                        }),
                    },
                )
            })
            .collect();
        assert_eq!(
            prune_beam(candidates, &mut accounting, CONFLICT_BEAM_WIDTH, 0).len(),
            128
        );
        assert_eq!(accounting.dropped_ties, 25);
        assert_eq!(accounting.dropped_non_ties, 0);
    }

    #[test]
    fn consecutive_split_boundaries_never_profile_physical_pairs() {
        let links = (0..3)
            .flat_map(|left| (0..3).map(move |right| (left, right)))
            .collect::<Vec<_>>();
        let ends = vec![b"ACGTACGT".to_vec(); 3];
        let starts = vec![b"TGCATGCA".to_vec(); 3];
        let mut profiling_calls = 0usize;
        for _boundary in 0..2 {
            let (profiles, stats) =
                deduplicated_boundary_profiles(&links, &ends, &starts, |_, _| {
                    profiling_calls += 1;
                    Ok(BTreeMap::from([(vec![4], 1)]))
                })
                .unwrap();
            assert_eq!(profiles.len(), links.len());
            assert_eq!(stats.physical_links, 9);
            assert_eq!(stats.class_compositions, 1);
            assert_eq!(stats.reused_physical_links, 8);
            assert_eq!(stats.per_physical_pair_compositions, 0);
        }
        assert_eq!(profiling_calls, 2);
    }

    #[test]
    fn large_beam_class_guided_path_is_byte_identical_to_exact_physical_path() {
        let records = (0..4)
            .flat_map(|locus| (0..3).map(move |allele| (vec![100 + locus * 10 + allele * 2], 1)))
            .collect::<BTreeMap<_, _>>();
        let sample = WeightedBwt::build(&records).unwrap();
        let partitions = (0..4)
            .map(|locus| {
                (0..3)
                    .map(|allele| GenomeAllele {
                        traversal: SpanningTraversal::single(SourceRange {
                            partition: locus as usize,
                            occurrence: locus as usize * 3 + allele as usize,
                            source: allele as usize,
                            start: locus * 10,
                            end: (locus + 1) * 10,
                            reverse: false,
                        }),
                        profile: BTreeMap::from([(vec![100 + locus * 10 + allele * 2], 1)]),
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let seams = (0..3)
            .map(|_| {
                (0..3)
                    .flat_map(|left| {
                        (0..3).map(move |right| GenomeSeam {
                            left,
                            right,
                            profile: Arc::new(Profile::new()),
                        })
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let incidence = incidence_table(&partitions, &seams);
        let model = ScoreModel {
            read_length: 150,
            histogram: 1,
            denominator: 150.0,
            depth: 150.0,
            background: 0.1,
        };
        let exact = exact_chain(
            &partitions,
            &seams,
            &incidence,
            &sample,
            &model,
            ProfileCost::default(),
        )
        .unwrap();
        let bounded = exact_chain_streaming(
            &partitions,
            &seams,
            &incidence,
            &sample,
            &model,
            ProfileCost::default(),
            true,
            8192,
            None,
            None,
            1,
            |locus| {
                Ok(partitions[locus]
                    .iter()
                    .map(|allele| allele.profile.clone())
                    .collect())
            },
        )
        .unwrap();
        assert_eq!(bounded.accounting.dropped_ties, 0);
        assert_eq!(bounded.accounting.dropped_non_ties, 0);
        assert_eq!(bounded.accounting.global_early_stop_class_pairs, 0);
        let mut bounded_choices = bounded.choices.clone();
        let mut exact_choices = exact.choices.clone();
        bounded_choices.sort();
        exact_choices.sort();
        assert_eq!(
            serde_json::to_vec(&bounded_choices).unwrap(),
            serde_json::to_vec(&exact_choices).unwrap()
        );
        assert_eq!(
            bounded
                .proposal_losses
                .iter()
                .map(|loss| loss.to_bits())
                .collect::<Vec<_>>(),
            exact
                .proposal_losses
                .iter()
                .map(|loss| loss.to_bits())
                .collect::<Vec<_>>()
        );
        let narrow = exact_chain_streaming(
            &partitions,
            &seams,
            &incidence,
            &sample,
            &model,
            ProfileCost::default(),
            true,
            1,
            None,
            None,
            1,
            |locus| {
                Ok(partitions[locus]
                    .iter()
                    .map(|allele| allele.profile.clone())
                    .collect())
            },
        )
        .unwrap();
        assert!(narrow.accounting.dropped_ties + narrow.accounting.dropped_non_ties > 0);
    }

    #[test]
    fn persistent_count_map_matches_flat_map_under_random_deltas() {
        let xorshift = |x: &mut u64| {
            *x ^= *x << 13;
            *x ^= *x >> 7;
            *x ^= *x << 17;
            *x
        };
        let mut state = 0x9e3779b97f4a7c15u64;
        // Random insertion/merge sequences vs a flat reference map.
        let mut flat: FxHashMap<u32, u64> = FxHashMap::default();
        let mut persistent = PersistentCountMap::default();
        for round in 0..200 {
            let mut delta = Vec::new();
            let count = 1 + xorshift(&mut state) % 50;
            for _ in 0..count {
                let fid = (xorshift(&mut state) % 100_000) as u32;
                let added = 1 + xorshift(&mut state) % 7;
                delta.push((fid, added));
            }
            delta.sort_unstable_by_key(|(fid, _)| *fid);
            delta.dedup_by_key(|(fid, _)| *fid);
            for &(fid, added) in &delta {
                *flat.entry(fid).or_default() += added;
            }
            persistent = persistent.apply_delta(&delta).unwrap();
            for &(fid, _) in &delta {
                assert_eq!(
                    persistent.get(fid),
                    flat.get(&fid).copied().unwrap_or(0),
                    "round {round} fid {fid}"
                );
            }
            if round % 37 == 0 {
                let mut streamed: Vec<(u32, u64)> = persistent.iter().collect();
                streamed.sort_unstable();
                let mut flat_entries: Vec<(u32, u64)> =
                    flat.iter().map(|(&f, &c)| (f, c)).collect();
                flat_entries.sort_unstable();
                assert_eq!(streamed, flat_entries, "round {round} full compare");
            }
        }
    }

    #[test]
    fn incremental_extend_fingerprint_and_score_match_from_scratch() {
        // The extend path's incremental invariants, verified at every node
        // of randomized chains against from-scratch recomputation:
        // (i) the commutative XOR content hash maintained per-delta-feature
        //     equals the XOR of per-(feature, cumulative count) hashes over
        //     the full streamed cumulative profile, independent of the
        //     component-slot assignment and the delta decomposition;
        // (ii) the telescoped score equals the summed loss over the full
        //     profile (to f64 summation-order tolerance);
        // (iii) extending the same parent with permuted component slots
        //     interns to the SAME node;
        // (iv) a from-None extension carrying the identical full content
        //     interns to that same node via the bucket verification path.
        let xorshift = |x: &mut u64| {
            *x ^= *x << 13;
            *x ^= *x >> 7;
            *x ^= *x << 17;
            *x
        };
        let mut state = 0x2545f4914f6cdd1du64;
        const FEATURES: usize = 48;
        let keys: Vec<FeatureKey> = (0..FEATURES).map(|f| vec![f as u64]).collect();
        let mut ids = FxHashMap::default();
        for (fid, key) in keys.iter().enumerate() {
            ids.insert(key.clone(), fid as u32);
        }
        let mut base_hash = Vec::with_capacity(FEATURES);
        for key in &keys {
            let mut hash = 0xcbf29ce484222325u64;
            hash = fnv_update(hash, &(key.len() as u64).to_le_bytes());
            for token in key {
                hash = fnv_update(hash, &token.to_le_bytes());
            }
            base_hash.push(hash);
        }
        let observed: Vec<u64> = (0..FEATURES)
            .map(|_| 1 + xorshift(&mut state) % 50)
            .collect();
        let universe = FeatureUniverse {
            keys,
            ids,
            base_hash,
            observed: observed.clone(),
        };
        let model = ScoreModel {
            read_length: 150,
            histogram: 1,
            denominator: 150.0,
            depth: 150.0,
            background: 4.0,
        };
        let mut interner = PersistentProfileInterner::default();
        let random_component = |state: &mut u64| -> Arc<ConvertedProfile> {
            let count = xorshift(state) % 6;
            let mut fids = (0..count)
                .map(|_| xorshift(state) % FEATURES as u64)
                .collect::<Vec<_>>();
            fids.sort_unstable();
            fids.dedup();
            Arc::new(
                fids.into_iter()
                    .map(|fid| (fid as u32, 1 + xorshift(state) % 5))
                    .collect::<Vec<_>>(),
            )
        };
        // Chain of random extends; the parent map is rebuilt from the
        // child's full streamed cumulative profile each round, which also
        // cross-checks from_sorted against the incremental content.
        let mut parent: Option<Arc<PersistentProfile>> = None;
        let mut parent_map = PersistentCountMap::default();
        let mut round_fingerprints = Vec::new();
        let mut owned_pair_constants = FxHashMap::default();
        for _round in 0..40 {
            let components = [
                random_component(&mut state),
                random_component(&mut state),
                random_component(&mut state),
                random_component(&mut state),
            ];
            let base_score = parent.as_ref().map_or(0.0, |profile| profile.score);
            let constants = pair_fresh_constants(
                None,
                &mut owned_pair_constants,
                0,
                _round,
                0,
                &components,
                &universe,
                &model,
            )
            .unwrap();
            // The production call path's small map: the parent's touched
            // counts (the test universe covers every delta feature).
            let small: FxHashMap<u32, u64> = {
                let mut entries = parent_map.iter().collect::<Vec<_>>();
                entries.sort_unstable();
                entries
                    .into_iter()
                    .filter(|&(_, count)| count > 0)
                    .collect()
            };
            let (score, cumulative) = interner
                .extend(
                    parent.as_ref(),
                    Some(&small),
                    components.clone(),
                    base_score,
                    constants,
                    &universe,
                    &model,
                )
                .unwrap();
            assert_eq!(score, cumulative.score);
            // From-scratch recompute over the full streamed profile.
            let mut full = Vec::new();
            {
                let mut pairs = cumulative.cumulative_pairs();
                while let Some((fid, count)) = pairs.next() {
                    full.push((fid, count));
                }
            }
            let mut scratch_hash = 0u64;
            let mut scratch_score = 0.0f64;
            for &(fid, count) in &full {
                scratch_hash ^= universe.fingerprint_value(fid, count);
                scratch_score += model.loss(count, observed[fid as usize]).unwrap();
            }
            assert_eq!(
                cumulative.fingerprint, scratch_hash,
                "round {_round}: incremental fingerprint != from-scratch XOR"
            );
            assert_eq!(cumulative.entries, full.len());
            assert!(
                (score - scratch_score).abs() <= 1e-6 * scratch_score.abs().max(1.0),
                "round {_round}: telescoped score {score} != from-scratch {scratch_score}"
            );
            round_fingerprints.push(scratch_hash);
            // (iii) commutativity across component slots: the same parent
            // extended with the quadruple in a rotated slot order must
            // intern to the same node.
            let rotated = [
                Arc::clone(&components[2]),
                Arc::clone(&components[0]),
                Arc::clone(&components[3]),
                Arc::clone(&components[1]),
            ];
            let rotated_constants =
                compute_pair_fresh_constants(&rotated, &universe, &model).unwrap();
            let (_, rotated_node) = interner
                .extend(
                    parent.as_ref(),
                    Some(&small),
                    rotated,
                    base_score,
                    rotated_constants,
                    &universe,
                    &model,
                )
                .unwrap();
            assert!(
                Arc::ptr_eq(&cumulative, &rotated_node),
                "round {_round}: permuted slots interned to a different node"
            );
            parent = Some(Arc::clone(&cumulative));
            parent_map = PersistentCountMap::from_sorted(&full);
        }
        let chain_tip = parent.expect("nonempty chain");
        // (iv) identical full content extended from None must dedup to the
        // chain tip through the (fingerprint, entries) bucket + streamed
        // byte verification path.
        let mut full = Vec::new();
        {
            let mut pairs = chain_tip.cumulative_pairs();
            while let Some((fid, count)) = pairs.next() {
                full.push((fid, count));
            }
        }
        let full_components = [
            Arc::new(full.clone()),
            empty_component(),
            empty_component(),
            empty_component(),
        ];
        let from_scratch = interner
            .extend(
                None,
                None,
                full_components.clone(),
                0.0,
                compute_pair_fresh_constants(&full_components, &universe, &model).unwrap(),
                &universe,
                &model,
            )
            .unwrap()
            .1;
        assert!(
            Arc::ptr_eq(&chain_tip, &from_scratch),
            "from-None identical content did not intern to the chain tip"
        );
        let mut distinct = Vec::new();
        for &value in &round_fingerprints {
            if !distinct.contains(&value) {
                distinct.push(value);
            }
        }
        // A degenerate constant fingerprint would make the dedup checks
        // vacuous; the randomized chain must actually move through content
        // space.
        assert!(distinct.len() >= 32, "degenerate fingerprint sequence");
    }

    #[test]
    fn extend_overlap_and_fresh_deltas_match_naive_path() {
        // Targeted equivalence for the telescoped extend paths against the
        // naive per-feature reference loop, over randomized parents and
        // deltas with both overlapping (the state already touched the
        // feature) and non-overlapping (fresh) features:
        // (i) the production fold path matches the naive loop BIT-EXACTLY
        //     (score, fingerprint, entries);
        // (ii) the dev telescoped path matches fingerprint and entries
        //     EXACTLY and the score to f64 summation-order tolerance
        //     (identical k_pair + identical corrections in identical order,
        //     so count-identical children stay exactly tied);
        // (iii) the dev telescoped path is deterministic: the same extend
        //     recomputed bit-identically (tie-set structure survives).
        let xorshift = |x: &mut u64| {
            *x ^= *x << 13;
            *x ^= *x >> 7;
            *x ^= *x << 17;
            *x
        };
        let mut state = 0x853c49e6748fea9bu64;
        const FEATURES: usize = 32;
        let keys: Vec<FeatureKey> = (0..FEATURES).map(|f| vec![f as u64]).collect();
        let mut ids = FxHashMap::default();
        for (fid, key) in keys.iter().enumerate() {
            ids.insert(key.clone(), fid as u32);
        }
        let mut base_hash = Vec::with_capacity(FEATURES);
        for key in &keys {
            let mut hash = 0xcbf29ce484222325u64;
            hash = fnv_update(hash, &(key.len() as u64).to_le_bytes());
            for token in key {
                hash = fnv_update(hash, &token.to_le_bytes());
            }
            base_hash.push(hash);
        }
        let observed: Vec<u64> = (0..FEATURES)
            .map(|_| 1 + xorshift(&mut state) % 30)
            .collect();
        let universe = FeatureUniverse {
            keys,
            ids,
            base_hash,
            observed: observed.clone(),
        };
        let model = ScoreModel {
            read_length: 150,
            histogram: 1,
            denominator: 150.0,
            depth: 150.0,
            background: 4.0,
        };
        let random_component = |state: &mut u64| -> Arc<ConvertedProfile> {
            let count = xorshift(state) % 8;
            let mut fids = (0..count)
                .map(|_| xorshift(state) % FEATURES as u64)
                .collect::<Vec<_>>();
            fids.sort_unstable();
            fids.dedup();
            Arc::new(
                fids.into_iter()
                    .map(|fid| (fid as u32, 1 + xorshift(state) % 5))
                    .collect::<Vec<_>>(),
            )
        };
        // Random cumulative parents (mix of touched and untouched features)
        // with random component quadruples.
        let mut overlap_rounds = 0usize;
        for round in 0..200 {
            let mut parent_map = PersistentCountMap::default();
            let mut parent_counts: Vec<(u32, u64)> = Vec::new();
            let touched = xorshift(&mut state) % (FEATURES as u64 / 2);
            let mut fids = (0..touched)
                .map(|_| xorshift(&mut state) % FEATURES as u64)
                .collect::<Vec<_>>();
            fids.sort_unstable();
            fids.dedup();
            for fid in fids {
                parent_counts.push((fid as u32, 1 + xorshift(&mut state) % 9));
            }
            if !parent_counts.is_empty() {
                parent_map = PersistentCountMap::from_sorted(&parent_counts);
            }
            let small: FxHashMap<u32, u64> = parent_counts.iter().copied().collect();
            let base_score = (1 + xorshift(&mut state) % 97) as f64 / 7.0;
            let components = [
                random_component(&mut state),
                random_component(&mut state),
                random_component(&mut state),
                random_component(&mut state),
            ];
            let constants = compute_pair_fresh_constants(&components, &universe, &model).unwrap();
            // Naive reference: the historical per-feature loop over the
            // merged stream.
            let naive = {
                let mut score = base_score;
                let mut fingerprint = 0u64;
                let mut entries = 0usize;
                let mut merge = ComponentMerge::new(&components);
                let mut touched = FxHashMap::<u32, u64>::default();
                for &(fid, count) in &parent_counts {
                    touched.insert(fid, count);
                }
                while let Some((fid, added)) = merge.next().unwrap() {
                    let previous = touched.get(&fid).copied().unwrap_or(0);
                    let next = previous + added;
                    score += model.loss(next, observed[fid as usize]).unwrap()
                        - model.loss(previous, observed[fid as usize]).unwrap();
                    if previous == 0 {
                        entries += 1;
                    } else {
                        fingerprint ^= universe.fingerprint_value(fid, previous);
                    }
                    fingerprint ^= universe.fingerprint_value(fid, next);
                }
                (score, fingerprint, entries)
            };
            // Production interner (fold path): bit-exact vs naive.
            let mut production = PersistentProfileInterner::default();
            let (score, node) = production
                .extend(
                    None,
                    Some(&small),
                    components.clone(),
                    base_score,
                    constants,
                    &universe,
                    &model,
                )
                .unwrap();
            assert_eq!(score, naive.0, "round {round}: production score bits");
            assert_eq!(
                node.fingerprint, naive.1,
                "round {round}: production fingerprint"
            );
            assert_eq!(node.entries, naive.2, "round {round}: production entries");
            // Dev telescoped interner: fingerprint/entries exact, score to
            // f64 tolerance.
            let mut telescoped = PersistentProfileInterner {
                telescoped: true,
                ..PersistentProfileInterner::default()
            };
            let (tscore, tnode) = telescoped
                .extend(
                    None,
                    Some(&small),
                    components.clone(),
                    base_score,
                    constants,
                    &universe,
                    &model,
                )
                .unwrap();
            assert_eq!(
                tnode.fingerprint, naive.1,
                "round {round}: telescoped fingerprint"
            );
            assert_eq!(tnode.entries, naive.2, "round {round}: telescoped entries");
            assert!(
                (tscore - naive.0).abs() <= 1e-9 * naive.0.abs().max(1.0),
                "round {round}: telescoped score {tscore} vs naive {}",
                naive.0
            );
            // Determinism: the same telescoped extend recomputates
            // bit-identically (count-identical states stay exactly tied).
            let mut telescoped_repeat = PersistentProfileInterner {
                telescoped: true,
                ..PersistentProfileInterner::default()
            };
            let (rscore, rnode) = telescoped_repeat
                .extend(
                    None,
                    Some(&small),
                    components.clone(),
                    base_score,
                    constants,
                    &universe,
                    &model,
                )
                .unwrap();
            assert_eq!(tscore, rscore, "round {round}: telescoped determinism");
            assert_eq!(tnode.fingerprint, rnode.fingerprint);
            assert_eq!(tnode.entries, rnode.entries);
            // The correction actually exercised overlap in some rounds;
            // otherwise the randomized test would be vacuous for the
            // overlap path.
            let overlap_engaged = components.iter().any(|component| component
                .iter()
                .any(|&(fid, _)| small.contains_key(&fid)));
            if overlap_engaged {
                overlap_rounds += 1;
            }
        }
        assert!(
            overlap_rounds >= 100,
            "randomized overlap engaged in only {overlap_rounds}/200 rounds"
        );
    }

    #[test]
    fn tie_rotation_passes_pool_finalists_and_merge_accounting_honestly() {
        // Seventeen identical-profile alleles per locus over disjoint spans
        // with all-legal seams: every state ties, so a width-4 beam truncates
        // a 153-member tie class at the initial layer. Rotation 0 must keep
        // the historical tie order; later rotations retain different — still
        // deterministic — subsets, and the pooled finalist set must contain
        // every rotation's routes.
        let sample = WeightedBwt::build(&BTreeMap::from([(vec![4], 1)])).unwrap();
        let allele = |partition: usize, occurrence: usize, start: u64| GenomeAllele {
            traversal: SpanningTraversal::single(SourceRange {
                partition,
                occurrence,
                source: 0,
                start,
                end: start + 10,
                reverse: false,
            }),
            profile: Profile::new(),
        };
        let partitions = vec![
            (0..17)
                .map(|occurrence| allele(0, occurrence, 0))
                .collect::<Vec<_>>(),
            (0..17)
                .map(|occurrence| allele(1, 17 + occurrence, 10))
                .collect::<Vec<_>>(),
        ];
        let seams = vec![(0..17)
            .flat_map(|left| {
                (0..17).map(move |right| GenomeSeam {
                    left,
                    right,
                    profile: Arc::new(Profile::new()),
                })
            })
            .collect::<Vec<_>>()];
        let incidence = incidence_table(&partitions, &seams);
        let model = ScoreModel {
            read_length: 150,
            histogram: 1,
            denominator: 150.0,
            depth: 150.0,
            background: 0.1,
        };
        let run = |rotations: usize| {
            exact_chain_streaming(
                &partitions,
                &seams,
                &incidence,
                &sample,
                &model,
                ProfileCost::default(),
                true,
                4,
                None,
                None,
                rotations,
                |locus: usize| {
                    Ok(partitions[locus]
                        .iter()
                        .map(|allele| allele.profile.clone())
                        .collect())
                },
            )
            .unwrap()
        };
        let single = run(1);
        assert!(!single.complete);
        assert!(single.accounting.dropped_ties > 0);
        assert_eq!(single.rotation_best_scores.len(), 1);
        let rotated = run(3);
        assert_eq!(rotated.rotation_best_scores.len(), 3);
        // Identical profiles make every tied state's future identical, so
        // each rotation's best internal score is bit-identical.
        for score in &rotated.rotation_best_scores {
            assert_eq!(
                score.map(f64::to_bits),
                single.rotation_best_scores[0].map(f64::to_bits)
            );
        }
        // Honest accounting: every rotation pass drops its own tie truncation
        // set on top of rotation 0's counters (rotation 0 additionally runs
        // the seed ladder, so the growth is strictly monotone, not exactly 3x).
        assert!(rotated.accounting.dropped_ties > single.accounting.dropped_ties);
        assert!(rotated.accounting.work > single.accounting.work);
        // The pooled finalist set contains rotation 0's routes and grows with
        // the tie-diverse rotations.
        let contains = |haystack: &Vec<[Vec<usize>; 2]>, needle: &[Vec<usize>; 2]| {
            haystack.iter().any(|choice| choice == needle)
        };
        for choice in &single.choices {
            assert!(contains(&rotated.choices, choice));
        }
        assert!(rotated.choices.len() > single.choices.len());
        // Determinism: the same rotations reproduce the same pooled set.
        let replay = run(3);
        assert_eq!(
            serde_json::to_vec(&rotated.choices).unwrap(),
            serde_json::to_vec(&replay.choices).unwrap()
        );
    }

    #[test]
    fn midscale_eighty_kilobase_exact_dp_recovers_count_supported_switch() {
        const PARTITIONS: usize = 8;
        const BP: u64 = 10_000;
        let records = (0..PARTITIONS)
            .flat_map(|partition| {
                [
                    (vec![4 + partition as u64 * 4], 1),
                    (vec![6 + partition as u64 * 4], 1),
                ]
            })
            .collect::<BTreeMap<_, _>>();
        let sample = WeightedBwt::build(&records).unwrap();
        let mut partitions = Vec::new();
        for partition in 0..PARTITIONS {
            let mut alleles = Vec::new();
            for source in 0..2 {
                let source_key = if (partition < 4) == (source == 0) {
                    0
                } else {
                    1
                };
                let key = vec![4 + partition as u64 * 4 + source_key * 2];
                alleles.push(GenomeAllele {
                    traversal: SpanningTraversal::single(SourceRange {
                        partition,
                        occurrence: partition * 2 + source,
                        source,
                        start: partition as u64 * BP,
                        end: (partition as u64 + 1) * BP,
                        reverse: false,
                    }),
                    profile: BTreeMap::from([(key, 1)]),
                });
            }
            partitions.push(alleles);
        }
        let seams = (0..PARTITIONS - 1)
            .map(|_| {
                vec![
                    GenomeSeam {
                        left: 0,
                        right: 0,
                        profile: Arc::new(Profile::new()),
                    },
                    GenomeSeam {
                        left: 1,
                        right: 1,
                        profile: Arc::new(Profile::new()),
                    },
                ]
            })
            .collect::<Vec<_>>();
        let incidence = incidence_table(&partitions, &seams);
        let model = ScoreModel {
            read_length: 150,
            histogram: 1,
            denominator: 150.0,
            depth: 150.0,
            background: 0.1,
        };
        let result = exact_chain(
            &partitions,
            &seams,
            &incidence,
            &sample,
            &model,
            ProfileCost {
                mem_queries: 16,
                integrated_windows: 2 * (80_000 - 149),
            },
        )
        .unwrap();
        assert!(result.complete);
        assert_eq!(result.stop_reason, "exact_min_plus_dp");
        assert_eq!(result.accounting.exact_blocks, 1);
        assert_eq!(result.accounting.occupancy.len(), PARTITIONS);
        assert!(!result.choices.is_empty());
        // Physical source continuity is exact; the best pair contains both
        // public 80 kb molecules and their profiles jointly explain both counts.
        let best = &result.choices[0];
        assert!(best[0].iter().all(|&choice| choice == best[0][0]));
        assert!(best[1].iter().all(|&choice| choice == best[1][0]));
        assert_ne!(best[0][0], best[1][0]);
        assert_eq!(result.accounting.bounded_blocks, 0);
    }

    /// Same-owner stitched spanning candidates (supervisor ruling,
    /// genome/stitching-omission-alignment 2026-09-24): adjacent same-source
    /// forward rows chain into a spanning candidate covering the window's
    /// axis interval; the stitched segments carry the ROWS' own owner
    /// partitions (the charging resolves owners per segment); multi-source
    /// adjacencies, gapped same-source pairs, reverse rows and single-row
    /// "chains" (already candidates) are never stitched.
    #[test]
    fn stitched_candidates_admit_only_adjacent_same_source_spanning_chains() {
        let row = |partition: usize, source: usize, start: u64, end: u64, reverse: bool| {
            SpanningTraversal::single(SourceRange {
                partition,
                occurrence: partition * 10_000_000 + start as usize,
                source,
                start,
                end,
                reverse,
            })
        };
        // Source 5: two ADJACENT rows with DIFFERENT owner partitions (the
        // measured donor shape: SK1 rows under dual-role groups). Source 6:
        // one row already spanning. Source 7: a GAPPED pair. Plus a reverse
        // row of source 5.
        let rows = vec![
            row(2, 5, 100, 200, false),
            row(3, 5, 200, 300, false),
            row(4, 6, 100, 300, false),
            row(5, 7, 100, 200, false),
            row(5, 7, 350, 450, false),
            row(2, 5, 100, 200, true),
        ];
        // The window the two source-5 rows jointly span: exactly one chain,
        // two segments, each keeping its own owner partition and occurrence.
        let stitched = stitched_candidates(9, &rows, 150, 250);
        assert_eq!(stitched.len(), 1);
        let chain = &stitched[0];
        assert_eq!(chain.segments.len(), 2);
        assert_eq!(chain.segments[0].source, 5);
        assert_eq!(chain.segments[0].start, 100);
        assert_eq!(chain.segments[0].end, 200);
        assert_eq!(chain.segments[0].partition, 2);
        assert_eq!(chain.segments[1].source, 5);
        assert_eq!(chain.segments[1].start, 200);
        assert_eq!(chain.segments[1].end, 300);
        assert_eq!(chain.segments[1].partition, 3);
        assert!(!chain.segments[0].reverse && !chain.segments[1].reverse);
        // A window fully inside one row still admits the run's covering
        // chains (they jointly span it; the single row is already a
        // candidate and competes on merits).
        assert_eq!(stitched_candidates(9, &rows, 210, 260).len(), 1);
        assert_eq!(stitched_candidates(9, &rows, 120, 180).len(), 1);
        // The full source-5 span: the same single chain.
        assert_eq!(stitched_candidates(9, &rows, 100, 300).len(), 1);
        // Multi-source adjacency is never stitched (the mixed-source
        // exclusion stays): dropping source 5's second row leaves nothing
        // even though source 6's row is adjacent-compatible by coordinates.
        let single = vec![row(2, 5, 100, 200, false), row(4, 6, 200, 300, false)];
        assert!(stitched_candidates(9, &single, 150, 250).is_empty());
        // Gapped same-source rows are not stitched here (a gap would be a
        // novel junction, outside this step's admission).
        assert!(stitched_candidates(9, &rows, 150, 400).is_empty());
        // Reverse rows never stitch.
        assert!(stitched_candidates(9, &[row(2, 5, 100, 200, true), row(3, 5, 200, 300, true)], 150, 250).is_empty());
    }
}
