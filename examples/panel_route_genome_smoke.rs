//! Public-input chromosome-subset smoke for experimental genome-wide chaining.
#![recursion_limit = "512"]
#[path = "panel_route_diploid_search/mod.rs"]
mod search;

use clap::Parser;
use impg::{
    genome_inference::{panel_routes as routes, read_json, sample, PanelIdentity, COUNT_POLICY},
    graph::reverse_complement,
    syng::{SyncmerParams, SyngIndex},
};
use search::genome_wide::{self as genome, GenomeAllele, GenomeSeam, ProfileCost};
use std::{io, path::PathBuf, sync::atomic::{AtomicU64, Ordering}, time::Instant};
use search::partition::Profile;
use rayon::prelude::*;

#[derive(Parser)]
struct Options {
    #[arg(long)]
    panel: String,
    #[arg(long)]
    routes: PathBuf,
    #[arg(long)]
    sample: PathBuf,
    #[arg(long)]
    axis: PathBuf,
    #[arg(long)]
    bed_directory: PathBuf,
    #[arg(long)]
    cache: PathBuf,
    #[arg(long, default_value = "S288C#0#chrI")]
    component: String,
    #[arg(long, default_value_t = 2)]
    partitions: usize,
    #[arg(long, default_value_t = 2)]
    paths: usize,
    /// Use every public physical occurrence instead of the smoke subsample.
    #[arg(long)]
    all_paths: bool,
    #[arg(long, default_value_t = 10.0)]
    depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
    /// Count exact interior event runs for the complete component and stop.
    #[arg(long)]
    preflight_only: bool,
    /// Profile candidates and report exit-class/conflict audit without chaining.
    #[arg(long)]
    audit_only: bool,
    /// Previous pass JSON used for truth-free bounded split refinement.
    #[arg(long)]
    refinement_evidence: Option<PathBuf>,
    #[arg(long, default_value_t = 100.0)]
    refinement_margin: f64,
    #[arg(long, default_value_t = 8)]
    refinement_top_k: usize,
    #[arg(long, default_value_t = 4)]
    refinement_cut_top_n: usize,
    #[arg(long, default_value_t = 1024)]
    beam_width: usize,
    #[arg(long, default_value_t = 16)]
    rss_budget_gib: u64,
    #[arg(long)]
    refinement_candidates_only: bool,
    /// Reuse a completed full-profile cut-ranking checkpoint with the same inputs.
    #[arg(long)]
    refinement_ranks: Option<PathBuf>,
    /// Assessment-only exact objective audit: four route files define A0/A1/B0/B1.
    #[arg(long, num_args = 4)]
    objective_audit_routes: Vec<PathBuf>,
    /// Assessment-only candidate-domain reachability trace for one public route.
    #[arg(long)]
    trace_route: Option<PathBuf>,
    /// Run the component machinery over loci [START, END) of the loaded
    /// component axis (tract-slice development loop; also enables
    /// component-local loading).
    #[arg(long, num_args = 2, value_names = ["START", "END"])]
    locus_range: Option<Vec<usize>>,
    /// Component-local loading: open the profile cache filtered to the
    /// selected loci (building its sidecar offset index once) and skip the
    /// whole-artifact global ports verification that this smoke's per-cut
    /// per-source port lookups never use.
    #[arg(long)]
    component_local_load: bool,
    /// Force the bounded initial beam without refinement evidence: the initial
    /// physical pair layer is beam-pruned (state_bounded, reported honestly)
    /// instead of requiring the whole unexpanded domain to fit the logical
    /// state budget.
    #[arg(long)]
    force_initial_beam: bool,
    /// Reference route pair (two JSON segment files): rescored externally
    /// before the DP, used as the ladder's inherited realizable incumbent (the
    /// admissible bound then prunes every state that cannot beat it), and
    /// retained as a finalist in the external count-distinct selection so the
    /// reported result can only improve on it.
    #[arg(long, num_args = 2)]
    seed_incumbent_routes: Vec<PathBuf>,
    /// Deterministic tie-rotation passes over the main-width beam: rotation 0
    /// is the historical tie order; each further rotation retains a different
    /// deterministic subset of over-width tie classes, pooling tie-diverse
    /// finalists for the external selection. Each extra rotation re-runs the
    /// main DP pass.
    #[arg(long, default_value_t = 1)]
    tie_rotations: usize,
}

fn resident_bytes() -> io::Result<u64> {
    let status = std::fs::read_to_string("/proc/self/status")?;
    let line = status
        .lines()
        .find(|line| line.starts_with("VmRSS:"))
        .ok_or_else(|| io::Error::other("VmRSS absent from /proc/self/status"))?;
    let kib = line
        .split_whitespace()
        .nth(1)
        .ok_or_else(|| io::Error::other("invalid VmRSS line"))?
        .parse::<u64>()
        .map_err(io::Error::other)?;
    kib.checked_mul(1024)
        .ok_or_else(|| io::Error::other("VmRSS overflow"))
}

fn guard_rss(
    peak: &mut u64,
    budget: u64,
    stage: &str,
    allocations: serde_json::Value,
) -> io::Result<()> {
    let rss = resident_bytes()?;
    *peak = (*peak).max(rss);
    if rss > budget {
        let diagnostic = serde_json::json!({
            "model": "genome-smoke-rss-guard-v1",
            "stage": stage,
            "rss_bytes": rss,
            "budget_bytes": budget,
            "largest_allocation_categories": allocations,
        });
        eprintln!(
            "{}",
            serde_json::to_string_pretty(&diagnostic).map_err(io::Error::other)?
        );
        return Err(io::Error::other(format!(
            "RSS guard exceeded during {stage}: {rss} > {budget}"
        )));
    }
    Ok(())
}

fn traversal_endpoints(
    sources: &routes::Sources,
    traversal: &genome::SpanningTraversal,
    flank: usize,
) -> io::Result<(Vec<u8>, Vec<u8>)> {
    if traversal.segments.len() == 1 {
        let deletion = &traversal.segments[0];
        if deletion.start == deletion.end {
            let end = deletion
                .start
                .saturating_add(flank as u64)
                .min(sources.lanes[deletion.source].1);
            let mut head = sources.fetch(deletion.source, deletion.start, end)?;
            if deletion.reverse {
                head = reverse_complement(&head);
            }
            return Ok((head, Vec::new()));
        }
    }
    let mut head = Vec::with_capacity(flank);
    for segment in &traversal.segments {
        if head.len() == flank {
            break;
        }
        let take = (flank - head.len()).min((segment.end - segment.start) as usize) as u64;
        let mut part = if segment.reverse {
            sources.fetch(segment.source, segment.end - take, segment.end)?
        } else {
            sources.fetch(segment.source, segment.start, segment.start + take)?
        };
        if segment.reverse {
            part = reverse_complement(&part);
        }
        head.extend(part);
    }
    let mut tail_parts = Vec::new();
    let mut tail_len = 0usize;
    for segment in traversal.segments.iter().rev() {
        if tail_len == flank {
            break;
        }
        let take = (flank - tail_len).min((segment.end - segment.start) as usize) as u64;
        let mut part = if segment.reverse {
            sources.fetch(segment.source, segment.start, segment.start + take)?
        } else {
            sources.fetch(segment.source, segment.end - take, segment.end)?
        };
        if segment.reverse {
            part = reverse_complement(&part);
        }
        tail_len += part.len();
        tail_parts.push(part);
    }
    tail_parts.reverse();
    let tail = tail_parts.concat();
    Ok((head, tail))
}

fn traversal_sequence(
    sources: &routes::Sources,
    traversal: &genome::SpanningTraversal,
) -> io::Result<Vec<u8>> {
    let mut sequence = Vec::new();
    for segment in &traversal.segments {
        if segment.start == segment.end {
            continue;
        }
        let mut part = sources.fetch(segment.source, segment.start, segment.end)?;
        if segment.reverse {
            part = reverse_complement(&part);
        }
        sequence.extend(part);
    }
    Ok(sequence)
}

#[derive(serde::Deserialize)]
struct AuditRoute {
    segments: Vec<routes::Segment>,
}

fn assemble_audit_route(sources: &routes::Sources, path: &std::path::Path) -> io::Result<Vec<u8>> {
    let route: AuditRoute = read_json(path)?;
    let mut sequence = Vec::new();
    for segment in route.segments {
        let mut part = sources.fetch(segment.source, segment.start, segment.end)?;
        if segment.reverse {
            part = reverse_complement(&part);
        }
        sequence.extend(part);
    }
    Ok(sequence)
}

fn full_pair_profile(
    panel: &SyngIndex,
    sequences: [&[u8]; 2],
) -> io::Result<search::partition::Profile> {
    let mut total = search::partition::Profile::new();
    for sequence in sequences {
        let starts = sequence.len().saturating_sub(150).saturating_add(1);
        let mut lo = 0usize;
        while lo < starts {
            let hi = (lo + 512).min(starts);
            let (profile, _) =
                genome::profile_event_runs(panel, sequence, 150, lo, hi, genome::MAX_FEATURES)?;
            add_profile(&mut total, &profile)?;
            lo = hi;
        }
    }
    Ok(total)
}

fn add_profile(
    into: &mut search::partition::Profile,
    profile: &search::partition::Profile,
) -> io::Result<()> {
    for (key, &count) in profile {
        let value = into.entry(key.clone()).or_default();
        *value = value
            .checked_add(count)
            .ok_or_else(|| io::Error::other("composed traversal profile overflow"))?;
    }
    Ok(())
}

fn traversal_profile(
    panel: &SyngIndex,
    sources: &routes::Sources,
    cache: &mut genome::ProfileCache,
    traversal: &genome::SpanningTraversal,
) -> io::Result<(search::partition::Profile, ProfileCost)> {
    if traversal.segments.len() == 1 {
        let sequence = traversal_sequence(sources, traversal)?;
        let key = format!("allele:{}:L150", traversal.identity);
        let (profile, cost, reused) = cache.get_or_insert_with(&key, || {
            genome::profile_event_interior(panel, &sequence, 150, genome::MAX_FEATURES)
        })?;
        return Ok((profile, if reused { ProfileCost::default() } else { cost }));
    }
    // A multi-segment (split) traversal's merged profile is itself cached
    // under the traversal's own allele key: the merge of two partial records
    // plus the internal seam record costs three JSON parses and two BTreeMap
    // merges per candidate per run, which dominated the interior stage on
    // refined domains (split candidates outnumber physical alleles there).
    // On a hit the partial/seam derivation below is skipped entirely — the
    // merged record fully determines the profile.
    let full_key = format!("allele:{}:L150", traversal.identity);
    if let Some((profile, _, _)) = cache.get_if_cached(&full_key)? {
        return Ok((profile, ProfileCost::default()));
    }
    let mut combined = search::partition::Profile::new();
    let mut charged = ProfileCost::default();
    let mut sequences = Vec::new();
    for (side, segment) in traversal.segments.iter().enumerate() {
        let part = genome::SpanningTraversal {
            partition: traversal.partition,
            identity: format!(
                "partial:{side}:{}:{}-{}:{}",
                segment.occurrence, segment.start, segment.end, segment.reverse
            ),
            segments: vec![segment.clone()],
        };
        let sequence = traversal_sequence(sources, &part)?;
        let key = format!("allele:{}:L150", part.identity);
        let (profile, cost, reused) = cache.get_or_insert_with(&key, || {
            genome::profile_event_interior(panel, &sequence, 150, genome::MAX_FEATURES)
        })?;
        if !reused {
            charged.mem_queries = charged.mem_queries.saturating_add(cost.mem_queries);
            charged.integrated_windows = charged
                .integrated_windows
                .saturating_add(cost.integrated_windows);
        }
        add_profile(&mut combined, &profile)?;
        sequences.push(sequence);
    }
    let left = &traversal.segments[0];
    let right = &traversal.segments[1];
    let seam_key = format!(
        "internal:{}:{}>{}:{}:L150",
        left.occurrence, left.end, right.occurrence, right.start
    );
    let (profile, cost, reused) = cache.get_or_insert_with(&seam_key, || {
        genome::profile_event_seam(
            panel,
            &sequences[0],
            &sequences[1],
            150,
            genome::MAX_FEATURES,
        )
    })?;
    if !reused {
        charged.mem_queries = charged.mem_queries.saturating_add(cost.mem_queries);
        charged.integrated_windows = charged
            .integrated_windows
            .saturating_add(cost.integrated_windows);
    }
    add_profile(&mut combined, &profile)?;
    // Cache the merged split profile under its own allele key so later runs
    // parse one record instead of re-deriving from partials + seam. The
    // cache's per-record feature cap (MAX_FEATURES) applies to stored
    // records only; an oversized merge stays uncached and is re-derived,
    // preserving today's semantics for domains the cap would reject.
    if combined.len() <= genome::MAX_FEATURES {
        let charged_snapshot = charged;
        let (stored, _, _) = cache.get_or_insert_with(&full_key, || Ok((combined, charged)))?;
        return Ok((stored, charged_snapshot));
    }
    Ok((combined, charged))
}

fn profile_fingerprint_pair(profile: &search::partition::Profile) -> (u64, u64) {
    // 128-bit content fingerprint (two independent hashes) for per-locus
    // incidence dedup: observing the same profile twice at one locus is a
    // no-op (SeenAt::One(locus) is idempotent), so only distinct profiles
    // need observation. Two hashes make a silent collision astronomically
    // unlikely, and any error would surface as an accounting.feature_hashes
    // change in the equivalence check.
    let mut first = 0xcbf29ce484222325u64;
    let mut second = 0x9e3779b97f4a7c15u64;
    for (key, &count) in profile {
        first ^= key.len() as u64;
        second = (second ^ (key.len() as u64)).wrapping_mul(0x2545f4914f6cdd1d);
        for &word in key {
            first = (first ^ word).wrapping_mul(0x100000001b3);
            second = (second ^ word).wrapping_mul(0x100000001b3);
        }
        first = (first ^ count).wrapping_mul(0x100000001b3);
        second = (second ^ count).wrapping_mul(0x100000001b3);
    }
    (first, second)
}

fn charge(total: &mut ProfileCost, cost: ProfileCost, reused: bool) -> io::Result<()> {
    if !reused {
        total.mem_queries = total
            .mem_queries
            .checked_add(cost.mem_queries)
            .ok_or_else(|| io::Error::other("profile query count overflow"))?;
    }
    total.integrated_windows = total
        .integrated_windows
        .checked_add(cost.integrated_windows)
        .ok_or_else(|| io::Error::other("integrated window count overflow"))?;
    Ok(())
}

fn main() -> io::Result<()> {
    let started = Instant::now();
    let options = Options::parse();
    let rss_budget_bytes = options
        .rss_budget_gib
        .checked_mul(1024 * 1024 * 1024)
        .ok_or_else(|| io::Error::other("RSS budget overflow"))?;
    let mut peak_rss_bytes = 0u64;
    let identity = PanelIdentity::read(&options.panel)?;
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    // Development route artifacts are bound to their build compiler identity;
    // this read-only smoke needs only their public lane/source inventory.
    let graph: routes::Graph = read_json(&options.routes.join("graph.json"))?;
    if graph.panel != identity
        || !graph.generation_complete
        || graph.count_policy != COUNT_POLICY
        || graph.lanes.len() != panel.name_map.path_to_name.len()
    {
        return Err(io::Error::other("incompatible smoke route inventory"));
    }
    let sample = sample::SampleIndex::load(&options.sample, &identity)?;
    let input_seconds = started.elapsed().as_secs_f64();
    let sources = routes::Sources::open(
        &graph.source_paths,
        graph
            .lanes
            .iter()
            .map(|lane| (lane.name.clone(), lane.length))
            .collect(),
    )?;
    let axis = genome::load_axis(&options.axis)?;
    let (axis_intervals, all_ranges) = genome::load_bed_axis_partitions(
        &axis,
        &options.bed_directory,
        &graph
            .lanes
            .iter()
            .map(|lane| (lane.name.clone(), lane.length))
            .collect::<Vec<_>>(),
        &options.component,
        // A locus-range slice always loads the whole component axis so native
        // endpoint legality is evaluated against the true component ends;
        // component locality applies to profile loading and the DP, not the
        // public axis/BED inventory.
        (!options.preflight_only && !options.all_paths && options.locus_range.is_none())
            .then_some(options.partitions),
        None,
    )?;
    let target = graph
        .lanes
        .iter()
        .find(|lane| lane.name == options.component)
        .ok_or_else(|| io::Error::other("axis component is not a route lane"))?;
    let component_suffix = target
        .name
        .splitn(3, '#')
        .nth(2)
        .ok_or_else(|| io::Error::other("route lane lacks public component suffix"))?;
    let component_sources = graph
        .lanes
        .iter()
        .filter(|lane| lane.name.splitn(3, '#').nth(2) == Some(component_suffix))
        .map(|lane| lane.id)
        .collect::<std::collections::BTreeSet<_>>();
    let mut all_ranges = all_ranges;
    let domain_completion =
        genome::complete_forward_source_paths(&mut all_ranges, &component_sources)?;
    let mut all_ranges = all_ranges
        .into_iter()
        .map(|ranges| {
            ranges
                .into_iter()
                .map(genome::SpanningTraversal::single)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    genome::retain_native_endpoint_candidates(&mut all_ranges, target.id, target.length)?;
    if options.preflight_only {
        let mut cost = ProfileCost::default();
        let mut processed = 0usize;
        'partitions: for ranges in &all_ranges {
            for range in ranges {
                let sequence = traversal_sequence(&sources, range)?;
                if sequence.len() >= 150 {
                    let next = genome::event_profile_cost(
                        &panel,
                        &sequence,
                        150,
                        0,
                        sequence.len() - 149,
                    )?;
                    cost.mem_queries += next.mem_queries;
                    cost.integrated_windows += next.integrated_windows;
                }
                processed += 1;
                if cost.mem_queries > genome::MAX_PROFILE_WORK {
                    break 'partitions;
                }
            }
        }
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "model": "experimental-partition-genome-profile-preflight-v1",
                "component": options.component,
                "partitions": all_ranges.len(),
                "physical_alleles_total": all_ranges.iter().map(Vec::len).sum::<usize>(),
                "physical_alleles_processed": processed,
                "profile_work": cost.mem_queries,
                "integrated_windows": cost.integrated_windows,
                "ceiling": genome::MAX_PROFILE_WORK,
                "complete": cost.mem_queries <= genome::MAX_PROFILE_WORK,
                "stop_reason": (cost.mem_queries > genome::MAX_PROFILE_WORK)
                    .then_some("budget: partition_profile_work"),
            }))?
        );
        return Ok(());
    }
    let mut ranges = if options.all_paths {
        all_ranges
    } else {
        genome::connected_path_subsample(&axis_intervals, &all_ranges, options.paths)?
    };
    // Tract-slice mode: run the component machinery over [start, end) loci of
    // the loaded component axis. Candidate domains, native-endpoint retention,
    // and completion are applied to the whole component first, so sliced loci
    // hold exactly the candidates the full run would hold. Identity strings
    // (split candidates, cache keys, rank checkpoints) keep full-component
    // locus numbering via `slice_offset`; only array indexing is renumbered.
    let component_partition_count = axis_intervals.len();
    let (mut axis_intervals, mut ranges, slice_offset) = match options.locus_range.as_deref() {
        Some(&[start, end]) => {
            if start >= end {
                return Err(io::Error::other("locus range start must precede end"));
            }
            if end > ranges.len() || end > axis_intervals.len() {
                return Err(io::Error::other("locus range outside component partitions"));
            }
            (
                axis_intervals[start..end].to_vec(),
                ranges[start..end].to_vec(),
                start,
            )
        }
        Some(other) => {
            return Err(io::Error::other(format!(
                "--locus-range needs exactly two values, got {other:?}"
            )))
        }
        None => (axis_intervals, ranges, 0),
    };
    let slice_len = ranges.len();
    let slice_locus = move |full_locus: usize| {
        full_locus
            .checked_sub(slice_offset)
            .filter(|&locus| locus < slice_len)
    };
    let histogram = *sample
        .stats
        .read_lengths
        .get(&150)
        .ok_or_else(|| io::Error::other("sample lacks L150 histogram"))?;
    let model = search::partition::ScoreModel {
        read_length: 150,
        histogram,
        denominator: (150 * histogram) as f64,
        depth: options.depth,
        background: options.background,
    };
    let cache_started = Instant::now();
    let mut cache = if options.component_local_load {
        // Component-local profile loading: retain seam-class and internal-seam
        // classes (tiny) plus exactly the allele identities the selected loci
        // request. Allele cache keys carry full-component locus numbering
        // (single/partial carry the partition-scaled occurrence, split
        // candidates the full locus), so the slice filter is a superset of
        // every key this run can look up or insert.
        let ranges_slice_len = ranges.len();
        let slice_retains = move |key: &str| {
            if !key.starts_with("allele:") || !key.ends_with(":L150") {
                // seam-class:*, internal:*, and any future class families are
                // tiny shared classes.
                return true;
            }
            let identity = &key["allele:".len()..key.len() - ":L150".len()];
            let mut fields = identity.split(':');
            let partition = match fields.next() {
                Some("single") => fields
                    .next()
                    .and_then(|occ| occ.parse::<u64>().ok())
                    .map(|occ| (occ / 10_000_000) as usize),
                Some("partial") => fields
                    .nth(1)
                    .and_then(|occ| occ.parse::<u64>().ok())
                    .map(|occ| (occ / 10_000_000) as usize),
                Some("split") => fields
                    .next()
                    .and_then(|locus| locus.parse::<usize>().ok()),
                // Unattributed identities (fixtures, future forms) stay.
                _ => None,
            };
            match partition {
                Some(partition) => {
                    partition >= slice_offset && partition < slice_offset + ranges_slice_len
                }
                None => true,
            }
        };
        genome::ProfileCache::open_filtered_guarded(&options.cache, &slice_retains, Some(rss_budget_bytes))?
    } else {
        genome::ProfileCache::open(&options.cache)?
    };
    let cache_open_seconds = cache_started.elapsed().as_secs_f64();
    guard_rss(
        &mut peak_rss_bytes,
        rss_budget_bytes,
        "cache_open",
        serde_json::json!({"cache": options.cache}),
    )?;
    let ports_started = Instant::now();
    let mut ports = if options.component_local_load {
        routes::Ports::open_without_global_verification(&options.routes, &graph)?
    } else {
        routes::Ports::open(&options.routes, &graph)?
    };
    let ports_open_seconds = ports_started.elapsed().as_secs_f64();
    let mut expanded_blocks = Vec::new();
    let mut split_candidate_counts = vec![0usize; ranges.len()];
    let mut split_prefilter_stats = vec![genome::SplitPrefilterStats::default(); ranges.len()];
    let mut prefilter_cost = ProfileCost::default();
    let mut checkpoint_policy_dropped = 0usize;
    // Temporary probe: sub-stage wall-clock breakdown of the stages that
    // "interior_load_or_profile" lumps together (refinement evidence parse,
    // rank checkpoint parse, split-candidate machinery, interior profile
    // loop). Printed into the final JSON under stage_probe; remove once the
    // slice budget is calibrated.
    let stage_probe_started = Instant::now();
    let mut probe_rank_parse_seconds = 0.0f64;
    let mut probe_split_seconds = 0.0f64;
    let mut probe_interior_loop_seconds = 0.0f64;
    let mut probe_interior_candidates = 0usize;
    // Sub-probes inside the interior profile loop: the parallel read-only
    // cache hits (JSON record parse), the sequential miss derivations, and
    // the incidence observation. Determines whether the loop is parse-bound
    // or bookkeeping-bound before any cache-format work.
    let mut probe_interior_cached_seconds = 0.0f64;
    let mut probe_interior_miss_seconds = 0.0f64;
    let mut probe_interior_miss_count = 0usize;
    let mut probe_interior_incidence_seconds = 0.0f64;
    let mut probe_interior_bytes = 0u64;
    let mut probe_evidence_parse_seconds = 0.0f64;
    let mut probe_retain_seconds = 0.0f64;
    if let Some(evidence_path) = &options.refinement_evidence {
        if options.refinement_top_k == 0 {
            return Err(io::Error::other("refinement top-K must be positive"));
        }
        let evidence_parse_started = Instant::now();
        let text = std::fs::read_to_string(evidence_path)?;
        let evidence: serde_json::Value = serde_json::from_str(
            &text[text
                .find('{')
                .ok_or_else(|| io::Error::other("refinement evidence has no JSON"))?..],
        )
        .map_err(io::Error::other)?;
        let local = evidence["accounting"]["local_pair_evidence"]
            .as_array()
            .ok_or_else(|| io::Error::other("refinement evidence lacks local pair diagnostics"))?;
        probe_evidence_parse_seconds += evidence_parse_started.elapsed().as_secs_f64();
        let rank_checkpoint = if let Some(path) = &options.refinement_ranks {
            let rank_started = Instant::now();
            let text = std::fs::read_to_string(path)?;
            let value: serde_json::Value = serde_json::from_str(
                &text[text
                    .find('{')
                    .ok_or_else(|| io::Error::other("rank checkpoint has no JSON"))?..],
            )
            .map_err(io::Error::other)?;
            let parsed = Some(
                serde_json::from_value::<Vec<genome::SplitPrefilterStats>>(
                    value["split_prefilter_stats"].clone(),
                )
                .map_err(io::Error::other)?,
            );
            probe_rank_parse_seconds += rank_started.elapsed().as_secs_f64();
            parsed
        } else {
            None
        };
        let mut expand = vec![false; ranges.len()];
        let mut selected_indices = vec![Vec::<usize>::new(); ranges.len()];
        for row in local {
            let full_locus = row["locus"]
                .as_u64()
                .ok_or_else(|| io::Error::other("invalid refinement locus"))?
                as usize;
            let Some(locus) = slice_locus(full_locus) else {
                continue;
            };
            if row["native_minus_best"].as_f64().unwrap_or(0.0) >= options.refinement_margin {
                expand[locus] = true;
            }
        }
        if let Some(copies) = evidence["final_result"]["selected_physical_alleles"].as_array() {
            for copy in copies {
                let alleles = copy
                    .as_array()
                    .ok_or_else(|| io::Error::other("invalid selected refinement route"))?;
                for (full_locus, allele) in alleles.iter().enumerate() {
                    let Some(locus) = slice_locus(full_locus) else {
                        continue;
                    };
                    if let Some(identity) = allele["identity"].as_str() {
                        if let Some(index) = ranges[locus]
                            .iter()
                            .position(|candidate| candidate.identity == identity)
                        {
                            if !selected_indices[locus].contains(&index) {
                                selected_indices[locus].push(index);
                            }
                        }
                    }
                }
                for boundary in 0..alleles.len().saturating_sub(1) {
                    let left = alleles[boundary]["segments"][0]["source"].as_u64();
                    let right = alleles[boundary + 1]["segments"][0]["source"].as_u64();
                    if left != right {
                        for full_locus in
                            boundary.saturating_sub(1)..=(boundary + 2).min(ranges.len() + slice_offset - 1)
                        {
                            if let Some(locus) = slice_locus(full_locus) {
                                expand[locus] = true;
                            }
                        }
                    }
                }
            }
        }
        for locus in 0..ranges.len() {
            let full_locus = slice_offset + locus;
            if !expand[locus] {
                continue;
            }
            let split_started = Instant::now();
            let mut admitted = local[full_locus]["top_class_representatives"]
                .as_array()
                .ok_or_else(|| io::Error::other("refinement evidence lacks top classes"))?
                .iter()
                .take(options.refinement_top_k)
                .filter_map(|value| value.as_u64().map(|value| value as usize))
                .filter(|&index| {
                    index < ranges[locus].len()
                        && ranges[locus][index].segments.len() == 1
                        && ranges[locus][index].segments[0].start
                            < ranges[locus][index].segments[0].end
                })
                .collect::<Vec<_>>();
            for &selected in &selected_indices[locus] {
                if ranges[locus][selected].segments.len() == 1
                    && ranges[locus][selected].segments[0].start
                        < ranges[locus][selected].segments[0].end
                    && !admitted.contains(&selected)
                {
                    if admitted.len() == options.refinement_top_k {
                        admitted.pop();
                    }
                    admitted.push(selected);
                }
            }
            let native_profile = ranges[locus]
                .iter()
                .find(|candidate| {
                    candidate.segments.len() == 1
                        && candidate.segments[0].source == target.id
                        && !candidate.segments[0].reverse
                })
                .map(|candidate| traversal_profile(&panel, &sources, &mut cache, candidate))
                .transpose()?
                .map(|(profile, _)| profile)
                .ok_or_else(|| io::Error::other("refinement locus lacks native baseline"))?;
            let (split, stats) = if let Some(checkpoint) = &rank_checkpoint {
                let stats = checkpoint
                    .get(full_locus)
                    .ok_or_else(|| io::Error::other("rank checkpoint locus mismatch"))?
                    .clone();
                let retained = stats
                    .pairs
                    .iter()
                    .flat_map(|pair| {
                        pair.retained
                            .iter()
                            .map(|(_, identity, _)| identity.as_str())
                    })
                    .collect::<std::collections::BTreeSet<_>>();
                let split =
                    genome::split_candidates(full_locus, &ranges[locus], &admitted, &graph, &mut ports)?
                        .into_iter()
                        .filter(|candidate| retained.contains(candidate.identity.as_str()))
                        .collect::<Vec<_>>();
                if split.len() != stats.retained_cuts {
                    // The checkpoint must be reproducible under the current
                    // binary: every regenerated candidate must have been
                    // retained by the checkpoint. A strict subset is accepted
                    // when the missing candidates are exactly the ones current
                    // legality rules reject (e.g. same-source reverse-order
                    // cuts with internally overlapping spans, excluded after
                    // the checkpoint was frozen); they would fail DP admission
                    // anyway. Extra candidates mean real domain drift and fail
                    // closed.
                    let regenerated = split
                        .iter()
                        .map(|candidate| candidate.identity.as_str())
                        .collect::<std::collections::BTreeSet<_>>();
                    let extras = regenerated
                        .difference(&retained)
                        .take(8)
                        .collect::<Vec<_>>();
                    if !extras.is_empty() {
                        return Err(io::Error::other(format!(
                            "rank checkpoint regeneration produced candidates outside the checkpoint at locus {full_locus}; extras: {extras:?}"
                        )));
                    }
                    checkpoint_policy_dropped += stats.retained_cuts - split.len();
                }
                (split, stats)
            } else {
                genome::split_candidates_prefiltered(
                    full_locus,
                    &ranges[locus],
                    &admitted,
                    &graph,
                    &mut ports,
                    options.refinement_cut_top_n,
                    options.refinement_candidates_only,
                    |traversal| {
                        let (profile, profile_cost) =
                            traversal_profile(&panel, &sources, &mut cache, traversal)?;
                        charge(&mut prefilter_cost, profile_cost, false)?;
                        let mut loss = 0.0;
                        let mut keys = profile
                            .keys()
                            .chain(native_profile.keys())
                            .collect::<Vec<_>>();
                        keys.sort();
                        keys.dedup();
                        for feature in keys {
                            let predicted = profile
                                .get(feature)
                                .copied()
                                .unwrap_or(0)
                                .checked_add(native_profile.get(feature).copied().unwrap_or(0))
                                .ok_or_else(|| io::Error::other("diploid prefilter overflow"))?;
                            loss += model.loss(predicted, sample.counts.count(feature)?)?;
                        }
                        Ok((loss, profile))
                    },
                )?
            };
            split_candidate_counts[locus] = split.len();
            split_prefilter_stats[locus] = stats;
            ranges[locus].extend(split);
            expanded_blocks.push(locus);
            probe_split_seconds += split_started.elapsed().as_secs_f64();
        }
        let retain_started = Instant::now();
        if options.locus_range.is_some() {
            genome::retain_native_endpoint_candidates_in_component_range(
                &mut ranges,
                target.id,
                target.length,
                slice_offset,
                component_partition_count - 1,
            )?;
        } else {
            genome::retain_native_endpoint_candidates(&mut ranges, target.id, target.length)?;
        }
        probe_retain_seconds += retain_started.elapsed().as_secs_f64();
    }
    if options.refinement_candidates_only {
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "model": "experimental-bounded-split-refinement-candidates-v1",
                "component": options.component,
                "refinement_margin": options.refinement_margin,
                "refinement_top_k": options.refinement_top_k,
                "refinement_cut_top_n": options.refinement_cut_top_n,
                "refinement_cut_ranking": "full_composed_profile_plus_native",
                "expanded_blocks": expanded_blocks,
                "domain_completion": domain_completion,
                "split_candidate_counts": split_candidate_counts,
                "split_prefilter_stats": split_prefilter_stats,
                "physical_alleles_per_block": ranges.iter().map(Vec::len).collect::<Vec<_>>(),
                "bounded_refinement": true,
            }))?
        );
        return Ok(());
    }
    guard_rss(
        &mut peak_rss_bytes,
        rss_budget_bytes,
        "candidate_construction",
        serde_json::json!({
            "physical_alleles": ranges.iter().map(Vec::len).sum::<usize>(),
            "split_candidates": split_candidate_counts.iter().sum::<usize>(),
        }),
    )?;
    let seam_derivation_started = Instant::now();
    let links = genome::port_word_seams(&axis_intervals, &ranges, &graph, &mut ports)?;
    let seam_derivation_seconds = seam_derivation_started.elapsed().as_secs_f64();
    guard_rss(
        &mut peak_rss_bytes,
        rss_budget_bytes,
        "port_word_seams",
        serde_json::json!({
            "physical_links": links.iter().map(Vec::len).sum::<usize>(),
        }),
    )?;
    if let Some(path) = &options.trace_route {
        let desired: routes::Route = read_json(path)?;
        let allowed = ranges
            .iter()
            .map(|locus| {
                locus
                    .iter()
                    .map(|candidate| {
                        candidate.segments.iter().all(|segment| {
                            desired.segments.iter().any(|target| {
                                segment.source == target.source
                                    && segment.reverse == target.reverse
                                    && target.start <= segment.start
                                    && segment.end <= target.end
                            })
                        })
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut reachable = allowed[0].clone();
        let mut predecessors = vec![vec![None; ranges[0].len()]];
        let mut first_unreachable_locus = None;
        for (boundary, boundary_links) in links.iter().enumerate() {
            let mut next = vec![false; ranges[boundary + 1].len()];
            let mut previous = vec![None; ranges[boundary + 1].len()];
            for &(left, right) in boundary_links {
                if reachable[left] && allowed[boundary + 1][right] {
                    next[right] = true;
                    previous[right].get_or_insert(left);
                }
            }
            if !next.iter().any(|&value| value) && first_unreachable_locus.is_none() {
                first_unreachable_locus = Some(boundary + 1);
            }
            reachable = next;
            predecessors.push(previous);
        }
        let terminal = reachable.iter().position(|&value| value);
        let mut selected = Vec::new();
        if let Some(mut allele) = terminal {
            selected.push(allele);
            for locus in (1..ranges.len()).rev() {
                allele = predecessors[locus][allele]
                    .ok_or_else(|| io::Error::other("trace predecessor missing"))?;
                selected.push(allele);
            }
            selected.reverse();
        }
        let selected_traversals = selected
            .iter()
            .enumerate()
            .map(|(locus, &allele)| &ranges[locus][allele])
            .collect::<Vec<_>>();
        let selected_segments = selected_traversals
            .iter()
            .flat_map(|candidate| candidate.segments.iter())
            .collect::<Vec<_>>();
        let route_cut_needles = desired
            .segments
            .windows(2)
            .map(|pair| (format!("@{}>", pair[0].end), format!("@{}", pair[1].start)))
            .collect::<Vec<_>>();
        let retained_route_cuts = split_prefilter_stats
            .iter()
            .flat_map(|stats| stats.pairs.iter())
            .flat_map(|pair| pair.retained.iter())
            .filter(|(_, identity, _)| {
                route_cut_needles.iter().any(|(left, right)| {
                    identity.contains(left.as_str()) && identity.contains(right.as_str())
                })
            })
            .collect::<Vec<_>>();
        let mut overlaps = Vec::new();
        for first in 0..selected_segments.len() {
            for second in first + 1..selected_segments.len() {
                let a = selected_segments[first];
                let b = selected_segments[second];
                if a.source == b.source && a.start < b.end && b.start < a.end {
                    overlaps.push((a.clone(), b.clone()));
                }
            }
        }
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "model": "assessment-route-domain-trace-v1",
                "component": options.component,
                "all_paths": options.all_paths,
                "refinement_ranks": options.refinement_ranks,
                "route": path,
                "domain_completion": domain_completion,
                "allowed_candidates_per_locus": allowed.iter().map(|row| row.iter().filter(|&&value| value).count()).collect::<Vec<_>>(),
                "allowed_candidates": allowed.iter().enumerate().map(|(locus, row)| serde_json::json!({
                    "locus": locus,
                    "candidates": ranges[locus].iter().zip(row).filter(|(_, &value)| value).map(|(candidate, _)| serde_json::json!({"identity": candidate.identity, "segments": candidate.segments})).collect::<Vec<_>>()
                })).collect::<Vec<_>>(),
                "unmatched_loci": allowed.iter().enumerate().filter(|(_, row)| !row.iter().any(|&value| value)).map(|(locus, _)| serde_json::json!({
                    "locus": locus,
                    "candidates": ranges[locus].iter().map(|candidate| serde_json::json!({"identity": candidate.identity, "segments": candidate.segments})).collect::<Vec<_>>()
                })).collect::<Vec<_>>(),
                "reachable": terminal.is_some(),
                "first_unreachable_locus": first_unreachable_locus,
                "selected_alleles": selected_traversals.iter().map(|candidate| serde_json::json!({"identity": candidate.identity, "segments": candidate.segments})).collect::<Vec<_>>(),
                "same_source_overlaps": overlaps,
                "retained_route_cuts": retained_route_cuts,
            }))?
        );
        return Ok(());
    }
    let mut cost = prefilter_cost;
    let mut partitions = Vec::new();
    let mut incidence = genome::IncidenceTable::empty();
    let mut locus_profile_work = Vec::new();
    // Parsed interior profiles, kept per locus so the DP's streaming loader
    // reuses this parse instead of re-reading and re-parsing every cache
    // record (the parse dominated the interior stage). The pool is bounded:
    // component-scale runs (169k alleles on full chrIII) would hold ~20GB of
    // parsed profiles, so beyond POOL_ALLELE_LIMIT the interior pass drops
    // the profiles and the DP loader re-parses from the cache in parallel.
    // Consumed locus by locus by exact_chain_streaming.
    const POOL_ALLELE_LIMIT: usize = 80_000;
    let pool_profiles = ranges.iter().map(Vec::len).sum::<usize>() <= POOL_ALLELE_LIMIT;
    let mut interior_profiles = Vec::<Option<Vec<Profile>>>::with_capacity(ranges.len());
    for ranges in &ranges {
        let interior_locus_started = Instant::now();
        let mut alleles = Vec::new();
        let mut locus_cost = ProfileCost::default();
        // Distinct-profile set for this locus: identical profiles repeat the
        // same (feature, locus) observations, which are idempotent, so only
        // the first occurrence of each distinct profile needs observation.
        let mut observed = std::collections::HashSet::new();
        // Parallel read-only cache hits for single-segment candidates and for
        // multi-segment (split) candidates whose merged profile was cached
        // under the traversal's own allele key; uncached splits and any cache
        // miss fall back to the sequential mutable derivation path below.
        let cached_started = Instant::now();
        let locus_cached_bytes = AtomicU64::new(0);
        let mut cached: Vec<Option<Profile>> = ranges
            .par_iter()
            .map(|range| {
                let key = format!("allele:{}:L150", range.identity);
                Ok(cache
                    .get_if_cached(&key)?
                    .map(|(profile, _, bytes)| {
                        locus_cached_bytes.fetch_add(bytes, Ordering::Relaxed);
                        profile
                    }))
            })
            .collect::<io::Result<Vec<_>>>()?;
        probe_interior_cached_seconds += cached_started.elapsed().as_secs_f64();
        probe_interior_bytes += locus_cached_bytes.load(Ordering::Relaxed);
        let mut locus_profiles = Vec::<Profile>::with_capacity(ranges.len());
        for (candidate_index, range) in ranges.iter().enumerate() {
            if candidate_index % 32 == 0 {
                guard_rss(
                    &mut peak_rss_bytes,
                    rss_budget_bytes,
                    "interior_profile_hot_loop",
                    serde_json::json!({
                        "locus": partitions.len(),
                        "candidate": candidate_index,
                        "locus_candidates": ranges.len(),
                    }),
                )?;
            }
            let (profile, profile_cost) = match cached[candidate_index].take() {
                Some(profile) => (profile, ProfileCost::default()),
                None => {
                    probe_interior_miss_count += 1;
                    let miss_started = Instant::now();
                    let derived =
                        traversal_profile(&panel, &sources, &mut cache, range).map_err(|error| {
                            io::Error::other(format!("profile {}: {error}", range.identity))
                        })?;
                    probe_interior_miss_seconds += miss_started.elapsed().as_secs_f64();
                    derived
                }
            };
            charge(&mut cost, profile_cost, false)?;
            charge(&mut locus_cost, profile_cost, false)?;
            if observed.insert(profile_fingerprint_pair(&profile)) {
                let incidence_started = Instant::now();
                incidence.observe_interior(partitions.len(), &profile);
                probe_interior_incidence_seconds += incidence_started.elapsed().as_secs_f64();
            }
            if pool_profiles {
                locus_profiles.push(profile);
            }
            alleles.push(GenomeAllele {
                traversal: range.clone(),
                profile: Default::default(),
            });
        }
        interior_profiles.push(pool_profiles.then_some(locus_profiles));
        partitions.push(alleles);
        locus_profile_work.push(locus_cost.mem_queries);
        probe_interior_loop_seconds += interior_locus_started.elapsed().as_secs_f64();
        probe_interior_candidates += ranges.len();
        guard_rss(
            &mut peak_rss_bytes,
            rss_budget_bytes,
            "interior_profiles",
            serde_json::json!({
                "completed_loci": partitions.len(),
                "physical_alleles": partitions.iter().map(Vec::len).sum::<usize>(),
            }),
        )?;
    }
    let interior_seconds = cache_started.elapsed().as_secs_f64() - cache_open_seconds;
    let boundary_started = Instant::now();
    let boundary_endpoints_started = Instant::now();
    let endpoints = ranges
        .iter()
        .map(|locus| {
            locus
                .iter()
                .map(|traversal| traversal_endpoints(&sources, traversal, 149))
                .collect::<io::Result<Vec<_>>>()
        })
        .collect::<io::Result<Vec<_>>>()?;
    let mut probe_boundary_endpoints_seconds = boundary_endpoints_started.elapsed().as_secs_f64();
    if std::env::var("IMPG_DB").is_ok() {
        let (calls, hits) = routes::sources_fetch_probe_stats();
        eprintln!("FETCH memo_miss {calls} memo_hits {hits}");
    }
    let mut probe_boundary_composition_seconds = 0.0f64;
    let mut probe_boundary_materialize_seconds = 0.0f64;
    let mut boundaries = Vec::new();
    let mut boundary_profile_work = Vec::new();
    let mut boundary_composition_stats = Vec::new();
    for (boundary, boundary_links) in links.iter().enumerate() {
        let left_ends = endpoints[boundary]
            .iter()
            .map(|(_, tail)| tail.clone())
            .collect::<Vec<_>>();
        let right_starts = endpoints[boundary + 1]
            .iter()
            .map(|(head, _)| head.clone())
            .collect::<Vec<_>>();
        let mut boundary_cost = ProfileCost::default();
        let mut composed_classes = 0usize;
        let composition_started = Instant::now();
        let (profiles, stats) = genome::deduplicated_boundary_profiles(
            boundary_links,
            &left_ends,
            &right_starts,
            |left_end, right_start| {
                if composed_classes % 32 == 0 {
                    guard_rss(
                        &mut peak_rss_bytes,
                        rss_budget_bytes,
                        "boundary_profile_hot_loop",
                        serde_json::json!({
                            "boundary": boundary,
                            "composed_classes": composed_classes,
                            "physical_links": boundary_links.len(),
                        }),
                    )?;
                }
                composed_classes += 1;
                let key = format!(
                    "seam-class:{}>{}:L150",
                    String::from_utf8_lossy(left_end),
                    String::from_utf8_lossy(right_start)
                );
                let (profile, profile_cost, reused) = cache.get_or_insert_with(&key, || {
                    genome::profile_event_seam(
                        &panel,
                        left_end,
                        right_start,
                        150,
                        genome::MAX_FEATURES,
                    )
                })?;
                charge(&mut cost, profile_cost, reused)?;
                charge(&mut boundary_cost, profile_cost, reused)?;
                Ok(profile)
            },
        )?;
        probe_boundary_composition_seconds += composition_started.elapsed().as_secs_f64();
        let materialize_started = Instant::now();
        let mut seams = Vec::with_capacity(boundary_links.len());
        for (physical_index, (&(left, right), profile)) in
            boundary_links.iter().zip(profiles).enumerate()
        {
            if physical_index % 65_536 == 0 {
                guard_rss(
                    &mut peak_rss_bytes,
                    rss_budget_bytes,
                    "boundary_materialization_hot_loop",
                    serde_json::json!({
                        "boundary": boundary,
                        "physical_link": physical_index,
                        "physical_links": boundary_links.len(),
                    }),
                )?;
            }
            incidence.observe_boundary(boundary, &profile);
            seams.push(GenomeSeam {
                left,
                right,
                profile,
            });
        }
        boundaries.push(seams);
        probe_boundary_materialize_seconds += materialize_started.elapsed().as_secs_f64();
        boundary_profile_work.push(boundary_cost.mem_queries);
        boundary_composition_stats.push(stats);
        guard_rss(
            &mut peak_rss_bytes,
            rss_budget_bytes,
            "boundary_profiles",
            serde_json::json!({
                "completed_boundaries": boundaries.len(),
                "completed_physical_links": boundaries.iter().map(Vec::len).sum::<usize>(),
                "class_compositions": boundary_composition_stats.iter().map(|stats| stats.class_compositions).sum::<usize>(),
            }),
        )?;
    }
    let boundary_seconds = boundary_started.elapsed().as_secs_f64();
    // Physical link pairs and endpoint words are construction-only. Keeping
    // them alongside class-interned seams duplicated the retained domain at
    // chain time and was enough to trip the RSS guard.
    drop(links);
    drop(endpoints);
    let audit_pre_started = Instant::now();
    let audit = genome::audit_chain(&partitions, &boundaries)?;
    let audit_chain_pre_seconds = audit_pre_started.elapsed().as_secs_f64();
    if !options.objective_audit_routes.is_empty() {
        if options.objective_audit_routes.len() != 4 {
            return Err(io::Error::other("objective audit requires A0 A1 B0 B1"));
        }
        let assembled = options
            .objective_audit_routes
            .iter()
            .map(|path| assemble_audit_route(&sources, path))
            .collect::<io::Result<Vec<_>>>()?;
        let profile_a = full_pair_profile(&panel, [&assembled[0], &assembled[1]])?;
        let profile_b = full_pair_profile(&panel, [&assembled[2], &assembled[3]])?;
        let mut keys = profile_a
            .keys()
            .chain(profile_b.keys())
            .cloned()
            .collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        let mut rows = Vec::new();
        let mut by_class = std::collections::BTreeMap::<String, (usize, f64, f64)>::new();
        let mut evaluator_delta = 0.0;
        let mut local_delta = 0.0;
        for feature in keys {
            let predicted_a = profile_a.get(&feature).copied().unwrap_or(0);
            let predicted_b = profile_b.get(&feature).copied().unwrap_or(0);
            if predicted_a == predicted_b {
                continue;
            }
            let observed = sample.counts.count(&feature)?;
            let loss_a = model.loss(predicted_a, observed)?;
            let loss_b = model.loss(predicted_b, observed)?;
            let delta = loss_a - loss_b;
            let class = incidence.classification_label(&feature);
            let attributed = if class == "deferred" { 0.0 } else { delta };
            evaluator_delta += delta;
            local_delta += attributed;
            let entry = by_class.entry(class.clone()).or_default();
            entry.0 += 1;
            entry.1 += delta;
            entry.2 += attributed;
            if delta.to_bits() != attributed.to_bits() {
                rows.push(serde_json::json!({
                    "feature": feature,
                    "class": class,
                    "predicted_a": predicted_a,
                    "predicted_b": predicted_b,
                    "observed": observed,
                    "loss_a": loss_a,
                    "loss_b": loss_b,
                    "evaluator_delta": delta,
                    "local_delta": attributed,
                    "missing_delta": delta - attributed,
                }));
            }
        }
        rows.sort_by(|a, b| {
            b["missing_delta"]
                .as_f64()
                .unwrap()
                .abs()
                .total_cmp(&a["missing_delta"].as_f64().unwrap().abs())
        });
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "model": "assessment-local-objective-decomposition-v1",
                "routes": options.objective_audit_routes,
                "evaluator_delta": evaluator_delta,
                "current_local_attributed_delta": local_delta,
                "missing_delta": evaluator_delta - local_delta,
                "by_class": by_class,
                "mismatching_features": rows,
                "feature_hashes": incidence.len(),
                "deferred_hashes": incidence.deferred(),
                "boundary_composition_stats": boundary_composition_stats,
            }))?
        );
        return Ok(());
    }
    if options.audit_only {
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "model": "experimental-partition-genome-exact-dp-audit-v1",
                "component": options.component,
                "partitions": partitions.len(),
                "physical_alleles": partitions.iter().map(Vec::len).sum::<usize>(),
                "legal_seams": boundaries.iter().map(Vec::len).sum::<usize>(),
                "seams_per_boundary": boundaries.iter().map(Vec::len).collect::<Vec<_>>(),
                "refinement_margin": options.refinement_margin,
                "refinement_top_k": options.refinement_top_k,
                "refinement_cut_top_n": options.refinement_cut_top_n,
                "refinement_cut_ranking": "full_composed_profile_plus_native",
                "expanded_blocks": expanded_blocks,
                "domain_completion": domain_completion,
                "split_candidate_counts": split_candidate_counts,
                "split_prefilter_stats": split_prefilter_stats,
                "locus_profile_work": locus_profile_work,
                "boundary_profile_work": boundary_profile_work,
                "profile_work": cost.mem_queries,
                "profile_ceiling": genome::MAX_PROFILE_WORK,
                "audit": audit,
                "complete": cost.mem_queries <= genome::MAX_PROFILE_WORK,
            }))?
        );
        return Ok(());
    }
    let incidence_seconds = 0.0;
    // Reference-route seeding: rescore the seed pair with the authoritative
    // external evaluator before the DP, validate its production port-word
    // seam legality, and hand the realized loss to the chain as the ladder's
    // inherited incumbent — the admissible bound then prunes every state
    // that cannot beat a realizable route. The pair is also retained as a
    // finalist so the external count-distinct selection can only improve on
    // it. An illegal seed is reported but neither seeds the bound nor competes
    // as a finalist (the current pipeline cannot realize it).
    struct SeedIncumbent {
        loss: Option<f64>,
        profile_work: u64,
        integrated_windows: u64,
        initial_runs: usize,
        routes: [routes::Route; 2],
        all_legal: bool,
        seam_reports: Vec<serde_json::Value>,
        native_endpoints_ok: bool,
        seconds: f64,
    }
    let seed_incumbent = if options.seed_incumbent_routes.is_empty() {
        None
    } else {
        if options.seed_incumbent_routes.len() != 2 {
            return Err(io::Error::other(
                "--seed-incumbent-routes needs exactly two route files",
            ));
        }
        let seed_started = Instant::now();
        let first: routes::Route = read_json(&options.seed_incumbent_routes[0])?;
        let second: routes::Route = read_json(&options.seed_incumbent_routes[1])?;
        let seed_routes = [first, second];
        let mut legality = Vec::new();
        let mut assembled = [Vec::new(), Vec::new()];
        let mut all_legal = true;
        let mut native_endpoints_ok = true;
        for (copy, route) in seed_routes.iter().enumerate() {
            let segments = &route.segments;
            if segments.is_empty() {
                return Err(io::Error::other("seed route is empty"));
            }
            let first = &segments[0];
            let last = segments.last().expect("nonempty");
            // A slice-local reference route (locus-range run) covers only
            // the slice's loci, so component-wide native endpoints do not
            // apply; its port-word seams are still validated below.
            let endpoints_ok = options.locus_range.is_some()
                || (first.start == 0
                    && !first.reverse
                    && !last.reverse
                    && last.end == graph.lanes[last.source].length);
            native_endpoints_ok &= endpoints_ok;
            let mut seam_reports = Vec::new();
            let mut copy_legal = true;
            for pair in segments.windows(2) {
                let left = &pair[0];
                let right = &pair[1];
                let left_port =
                    ports.at_cut(&graph, left.source, left.exit(), left.reverse)?;
                let right_port =
                    ports.at_cut(&graph, right.source, right.entry(), right.reverse)?;
                let legal = matches!((left_port, right_port), (Some(l), Some(r)) if l.word == r.word);
                copy_legal &= legal;
                seam_reports.push(serde_json::json!({
                    "left": {"source": left.source, "exit": left.exit()},
                    "right": {"source": right.source, "entry": right.entry()},
                    "legal": legal,
                }));
            }
            all_legal &= copy_legal & endpoints_ok;
            legality.push(serde_json::json!({
                "copy": copy,
                "segments": segments.len(),
                "native_endpoints_ok": endpoints_ok,
                "all_port_word_seams_legal": copy_legal,
                "seams": seam_reports,
            }));
            for segment in segments {
                let mut sequence = sources.fetch(segment.source, segment.start, segment.end)?;
                if segment.reverse {
                    sequence = reverse_complement(&sequence);
                }
                assembled[copy].extend(sequence);
            }
        }
        let (loss, rescore_cost, initial_runs) = if all_legal {
            let run_directory = options
                .cache
                .parent()
                .unwrap_or_else(|| std::path::Path::new("."))
                .join("seed-incumbent-rescore-runs");
            if run_directory.exists() {
                std::fs::remove_dir_all(&run_directory)?;
            }
            let (loss, cost, runs) = genome::external_rescore(
                &panel,
                [&assembled[0], &assembled[1]],
                &sample.counts,
                &model,
                &run_directory,
            )?;
            (Some(loss), cost, runs)
        } else {
            (None, ProfileCost::default(), 0)
        };
        Some(SeedIncumbent {
            loss,
            profile_work: rescore_cost.mem_queries,
            integrated_windows: rescore_cost.integrated_windows,
            initial_runs,
            routes: seed_routes,
            all_legal,
            seam_reports: legality,
            native_endpoints_ok,
            seconds: seed_started.elapsed().as_secs_f64(),
        })
    };
    let chain_started = Instant::now();
    // The DP consumes the interior pass's parsed profiles locus by locus
    // when they were pooled; otherwise it re-parses the locus's cache records
    // in parallel (the previous closure re-read and re-parsed sequentially).
    let stored_loci = std::cell::RefCell::new(interior_profiles);
    let chain_loader_cache = &mut cache;
    let chain_loader_ranges = &ranges;
    let chain_loader_panel = &panel;
    let chain_loader_sources = &sources;
    let chain = genome::exact_chain_streaming_with_audit(
        &partitions,
        &boundaries,
        &audit,
        &incidence,
        &sample.counts,
        &model,
        cost,
        // The refinement-evidence lineage keeps the forced initial beam: the
        // all-paths initial layer exceeds the logical-state budget, and the
        // forced beam is the documented honest bounded state for refined
        // runs (an unforced huge initial layer would also build a per-state
        // count map per initial state before the budget check fires). The
        // chain still reports its own completeness honestly.
        options.refinement_evidence.is_some() || options.force_initial_beam,
        options.beam_width,
        Some(rss_budget_bytes),
        seed_incumbent.as_ref().and_then(|seed| seed.loss),
        options.tie_rotations,
        move |locus| {
            if let Some(stored) = stored_loci.borrow_mut()[locus].take() {
                return Ok(stored);
            }
            let mut parsed: Vec<Option<Profile>> = chain_loader_ranges[locus]
                .par_iter()
                .map(|range| {
                    if range.segments.len() != 1 {
                        return Ok(None);
                    }
                    let key = format!("allele:{}:L150", range.identity);
                    Ok(chain_loader_cache
                        .get_if_cached(&key)?
                        .map(|(profile, _, _)| profile))
                })
                .collect::<io::Result<Vec<_>>>()?;
            chain_loader_ranges[locus]
                .iter()
                .enumerate()
                .map(|(index, range)| match parsed[index].take() {
                    Some(profile) => Ok(profile),
                    None => {
                        traversal_profile(chain_loader_panel, chain_loader_sources, chain_loader_cache, range)
                            .map(|(profile, _)| profile)
                            .map_err(|error| {
                                io::Error::other(format!("profile {}: {error}", range.identity))
                            })
                    }
                })
                .collect()
        },
    )?;
    let chain_seconds = chain_started.elapsed().as_secs_f64();
    guard_rss(
        &mut peak_rss_bytes,
        rss_budget_bytes,
        "chain_complete",
        serde_json::json!({
            "occupancy": &chain.accounting.occupancy,
            "logical_state_bytes": chain.accounting.state_bytes,
            "physical_nodes": chain.accounting.occupancy.iter().sum::<usize>(),
        }),
    )?;
    let final_started = Instant::now();
    let final_result = if chain.choices.is_empty() {
        // The DP found no completable route that can beat the inherited
        // incumbent (seeded runs may dead-end honestly under the bound). A
        // legal seed pair is still the best KNOWN realizable selection and
        // is reported as such, flagged selected_seed_reference.
        seed_incumbent.as_ref().and_then(|seed| {
            let loss = seed.loss?;
            seed.all_legal.then_some(serde_json::json!({
                    "exact_full_subwalk_loss": loss,
                    "profile_work": seed.profile_work,
                    "integrated_windows": seed.integrated_windows,
                    "external_initial_runs": seed.initial_runs,
                    "public_physical_routes": seed.routes,
                    "selected_physical_alleles": serde_json::Value::Null,
                    "selected_pair_index": 0,
                    "selected_seed_reference": true,
                    "physical_proposal_pairs": 0,
                    "count_distinct_pairs": 0,
                    "rescored_count_distinct_pairs": 1,
                    "rescored_candidate_pairs": [{
                        "candidate": serde_json::Value::Null,
                        "seed_reference": true,
                        "physical_members": [],
                        "proposal_loss": serde_json::Value::Null,
                        "seed_port_word_seams_legal": seed.all_legal,
                        "exact_full_subwalk_loss": loss,
                        "profile_work": seed.profile_work,
                        "integrated_windows": seed.integrated_windows,
                        "external_initial_runs": seed.initial_runs,
                        "public_physical_routes": seed.routes,
                        "selected_physical_alleles": serde_json::Value::Null,
                    }],
                    "validation": "seed reference pair (DP retained no state that can beat the inherited incumbent); authoritative Evaluator validation occurs after whole-axis assembly",
                }))
        })
    } else {
        let mut distinct = Vec::<(search::partition::Profile, Vec<usize>, f64)>::new();
        let mut ranked_candidates = (0..chain.choices.len()).collect::<Vec<_>>();
        ranked_candidates.sort_by(|&a, &b| {
            chain.proposal_losses[a]
                .total_cmp(&chain.proposal_losses[b])
                .then_with(|| a.cmp(&b))
        });
        for candidate in ranked_candidates {
            let choices = &chain.choices[candidate];
            let mut signature = search::partition::Profile::new();
            for copy in 0..2 {
                for (locus, &allele) in choices[copy].iter().enumerate() {
                    let (profile, _) =
                        traversal_profile(&panel, &sources, &mut cache, &ranges[locus][allele])?;
                    add_profile(&mut signature, &profile)?;
                }
                for boundary in 0..boundaries.len() {
                    let seam = boundaries[boundary]
                        .iter()
                        .find(|seam| {
                            seam.left == choices[copy][boundary]
                                && seam.right == choices[copy][boundary + 1]
                        })
                        .ok_or_else(|| io::Error::other("selected route lacks boundary seam"))?;
                    add_profile(&mut signature, &seam.profile)?;
                }
            }
            if let Some((_, members, _)) = distinct
                .iter_mut()
                .find(|(existing, _, _)| *existing == signature)
            {
                members.push(candidate);
            } else {
                distinct.push((signature, vec![candidate], chain.proposal_losses[candidate]));
            }
        }
        distinct.sort_by(|a, b| a.2.total_cmp(&b.2).then_with(|| a.1[0].cmp(&b.1[0])));
        let physical_proposal_pairs = chain.choices.len();
        let count_distinct_pairs = distinct.len();
        let retained_distinct = distinct
            .into_iter()
            .take(genome::MAX_COMPLETE_PAIRS)
            .collect::<Vec<_>>();
        let mut rescored = Vec::with_capacity(retained_distinct.len());
        let mut best_index = 0usize;
        let mut best_loss = f64::INFINITY;
        let mut selected_physical_classes = std::collections::BTreeSet::new();
        for (_signature, physical_members, proposal_loss) in retained_distinct {
            let physical_class = |candidate: usize| {
                (0..2)
                    .map(|copy| {
                        chain.choices[candidate][copy]
                            .iter()
                            .enumerate()
                            .map(|(locus, &allele)| {
                                ranges[locus][allele]
                                    .segments
                                    .iter()
                                    .map(|segment| segment.source)
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            };
            let candidate = physical_members
                .iter()
                .copied()
                .find(|&candidate| !selected_physical_classes.contains(&physical_class(candidate)))
                .unwrap_or(physical_members[0]);
            selected_physical_classes.insert(physical_class(candidate));
            let choices = &chain.choices[candidate];
            let selected = [
                choices[0]
                    .iter()
                    .enumerate()
                    .map(|(locus, &allele)| ranges[locus][allele].clone())
                    .collect::<Vec<_>>(),
                choices[1]
                    .iter()
                    .enumerate()
                    .map(|(locus, &allele)| ranges[locus][allele].clone())
                    .collect::<Vec<_>>(),
            ];
            let selected_segments = [
                selected[0]
                    .iter()
                    .flat_map(|allele| allele.segments.iter().cloned())
                    .collect::<Vec<_>>(),
                selected[1]
                    .iter()
                    .flat_map(|allele| allele.segments.iter().cloned())
                    .collect::<Vec<_>>(),
            ];
            let routes = [
                genome::coalesced_route(&selected_segments[0])?,
                genome::coalesced_route(&selected_segments[1])?,
            ];
            let mut assembled = [Vec::new(), Vec::new()];
            for copy in 0..2 {
                for (locus, &allele) in choices[copy].iter().enumerate() {
                    assembled[copy].extend(traversal_sequence(&sources, &ranges[locus][allele])?);
                }
            }
            let run_directory = options
                .cache
                .parent()
                .unwrap_or_else(|| std::path::Path::new("."))
                .join(format!("final-rescore-runs-{candidate}"));
            if run_directory.exists() {
                std::fs::remove_dir_all(&run_directory)?;
            }
            let (loss, rescore_cost, initial_runs) = genome::external_rescore(
                &panel,
                [&assembled[0], &assembled[1]],
                &sample.counts,
                &model,
                &run_directory,
            )?;
            if loss.total_cmp(&best_loss).is_lt() {
                best_loss = loss;
                best_index = rescored.len();
            }
            rescored.push(serde_json::json!({
                "candidate": candidate,
                "physical_members": physical_members,
                "proposal_loss": proposal_loss,
                "exact_full_subwalk_loss": loss,
                "profile_work": rescore_cost.mem_queries,
                "integrated_windows": rescore_cost.integrated_windows,
                "external_initial_runs": initial_runs,
                "public_physical_routes": routes,
                "selected_physical_alleles": selected,
            }));
        }
        // The seeded reference pair competes as a finalist only when its
        // production port-word seams are legal (loss is None otherwise), so
        // the selected external loss can only improve on the seed.
        if let Some(seed) = &seed_incumbent {
            if seed
                .loss
                .is_some_and(|loss| loss.total_cmp(&best_loss).is_lt())
            {
                best_loss = seed.loss.expect("checked");
                best_index = rescored.len();
            }
            rescored.push(serde_json::json!({
                "candidate": serde_json::Value::Null,
                "seed_reference": true,
                "physical_members": [],
                "proposal_loss": serde_json::Value::Null,
                "seed_port_word_seams_legal": seed.all_legal,
                "exact_full_subwalk_loss": seed.loss,
                "profile_work": seed.profile_work,
                "integrated_windows": seed.integrated_windows,
                "external_initial_runs": seed.initial_runs,
                "public_physical_routes": seed.routes,
                "selected_physical_alleles": serde_json::Value::Null,
            }));
        }
        let best = &rescored[best_index];
        Some(serde_json::json!({
            "exact_full_subwalk_loss": best["exact_full_subwalk_loss"],
            "profile_work": best["profile_work"],
            "integrated_windows": best["integrated_windows"],
            "external_initial_runs": best["external_initial_runs"],
            "public_physical_routes": best["public_physical_routes"],
            "selected_physical_alleles": best["selected_physical_alleles"],
            "selected_pair_index": best_index,
            "selected_seed_reference": best["seed_reference"].as_bool() == Some(true),
            "physical_proposal_pairs": physical_proposal_pairs,
            "count_distinct_pairs": count_distinct_pairs,
            "rescored_count_distinct_pairs": rescored.len(),
            "rescored_candidate_pairs": rescored,
            "validation": "public coordinate seams, native endpoint pairing, and within-copy span capacity validated; authoritative Evaluator validation occurs after whole-axis assembly",
        }))
    };
    let final_seconds = final_started.elapsed().as_secs_f64();
    let dp_probe = genome::DpStageProbe::snapshot();
    let records_sidecar = cache.records_sidecar_report();
    let seed_incumbent_report = seed_incumbent.as_ref().map(|seed| {
        serde_json::json!({
            "loss": seed.loss,
            "all_legal": seed.all_legal,
            "native_endpoints_ok": seed.native_endpoints_ok,
            "profile_work": seed.profile_work,
            "integrated_windows": seed.integrated_windows,
            "external_initial_runs": seed.initial_runs,
            "rescore_seconds": seed.seconds,
            "legality": seed.seam_reports,
        })
    });
    guard_rss(
        &mut peak_rss_bytes,
        rss_budget_bytes,
        "final_rescore_complete",
        serde_json::json!({
            "physical_proposals": chain.choices.len(),
            "logical_state_bytes": chain.accounting.state_bytes,
        }),
    )?;
    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::json!({
            "model": "experimental-partition-genome-exact-dp-v1",
            "component": options.component,
            "locus_range": options.locus_range,
            "component_local_load": options.component_local_load,
            "cache_entries_retained": cache.retained_entries(),
            "cache_index_sidecar": cache.loaded_from_sidecar(),
            "records_sidecar": {
                "used": records_sidecar.0,
                "covered_bytes": records_sidecar.1,
                "status": records_sidecar.2,
                "sync_seconds": records_sidecar.3,
            },
            "partitions": partitions.len(),
            "physical_alleles": partitions.iter().map(Vec::len).sum::<usize>(),
            "legal_seams": boundaries.iter().map(Vec::len).sum::<usize>(),
            "seams_per_boundary": boundaries.iter().map(Vec::len).collect::<Vec<_>>(),
            "refinement_margin": options.refinement_margin,
            "refinement_top_k": options.refinement_top_k,
            "refinement_cut_top_n": options.refinement_cut_top_n,
            "refinement_cut_ranking": "full_composed_profile_plus_native",
            "beam_width": options.beam_width,
            "tie_rotations": options.tie_rotations,
            "tie_rotation_best_scores": chain.rotation_best_scores,
            "seed_incumbent_routes": options.seed_incumbent_routes,
            "seed_incumbent": seed_incumbent_report,
            "rss_budget_bytes": rss_budget_bytes,
            "peak_rss_bytes": peak_rss_bytes,
            "expanded_blocks": expanded_blocks,
            "domain_completion": domain_completion,
            "checkpoint_policy_dropped": checkpoint_policy_dropped,
            "split_candidate_counts": split_candidate_counts,
            "split_prefilter_stats": split_prefilter_stats,
            "bounded_refinement": options.refinement_evidence.is_some(),
            "forced_initial_beam": options.force_initial_beam,
            "feature_hashes": incidence.len(),
            "deferred_hashes": incidence.deferred(),
            "hash_collisions": incidence.collisions,
            "locus_profile_work": locus_profile_work,
            "boundary_profile_work": boundary_profile_work,
            "boundary_composition_stats": boundary_composition_stats,
            "profile_ceiling": genome::MAX_PROFILE_WORK,
            "audit": audit,
            "stage_wall_seconds": {
                "input_load": input_seconds,
                "cache_open": cache_open_seconds,
                "ports_open": ports_open_seconds,
                "port_word_seam_derivation": seam_derivation_seconds,
                "interior_load_or_profile": interior_seconds,
                "boundary_load_or_profile": boundary_seconds,
                "incidence": incidence_seconds,
                "exact_dp_and_physical_expansion": chain_seconds,
                "external_final_rescore": final_seconds,
                "total": started.elapsed().as_secs_f64(),
            },
            "stage_probe": {
                "rank_checkpoint_parse_seconds": probe_rank_parse_seconds,
                "split_candidates_seconds": probe_split_seconds,
                "interior_loop_seconds": probe_interior_loop_seconds,
                "interior_cached_seconds": probe_interior_cached_seconds,
                "interior_cached_bytes": probe_interior_bytes,
                "interior_miss_seconds": probe_interior_miss_seconds,
                "interior_miss_count": probe_interior_miss_count,
                "interior_incidence_seconds": probe_interior_incidence_seconds,
                "evidence_parse_seconds": probe_evidence_parse_seconds,
                "retain_seconds": probe_retain_seconds,
                "audit_chain_pre_seconds": audit_chain_pre_seconds,
                "interior_candidates": probe_interior_candidates,
                "boundary_endpoints_seconds": probe_boundary_endpoints_seconds,
                "boundary_composition_seconds": probe_boundary_composition_seconds,
                "boundary_materialize_seconds": probe_boundary_materialize_seconds,
                "pre_dp_total_seconds": stage_probe_started.elapsed().as_secs_f64(),
                "dp_probe": {
                    "passes": dp_probe.1,
                    "successor_registries_seconds": dp_probe.0[0],
                    "feature_universe_seconds": dp_probe.0[1],
                    "initial_layer_seconds": dp_probe.0[2],
                    "seed_ladder_seconds": dp_probe.0[3],
                    "locus_weight_tables_seconds": dp_probe.0[4],
                    "suffix_max_add_tables_seconds": dp_probe.0[5],
                    "expansion_loop_seconds": dp_probe.0[6],
                    "layer_commit_seconds": dp_probe.0[7],
                    "suffix_rebase_seconds": dp_probe.0[8],
                    "terminal_collection_seconds": dp_probe.0[9],
                    "locus_class_evidence_seconds": dp_probe.0[10],
                    "streamed_profile_load_seconds": dp_probe.0[11],
                    "audit_chain_seconds": dp_probe.0[12],
                },
            },
            "candidate_pairs": chain.choices.len(),
            "final_result": final_result,
            "complete": chain.complete
                && options.refinement_evidence.is_none()
                && !options.force_initial_beam,
            "stop_reason": if options.refinement_evidence.is_some() {
                format!("bounded_refined:{}", chain.stop_reason)
            } else if options.force_initial_beam {
                format!("forced_initial_beam:{}", chain.stop_reason)
            } else {
                chain.stop_reason
            },
            "accounting": chain.accounting,
        }))?
    );
    Ok(())
}

#[cfg(test)]
mod rss_guard_tests {
    use super::*;

    #[test]
    fn rss_guard_fails_closed_at_a_hot_loop_checkpoint() {
        let mut peak = 0;
        let error = guard_rss(
            &mut peak,
            0,
            "test_hot_loop",
            serde_json::json!({"allocation": "provably resident process"}),
        )
        .unwrap_err();
        assert!(peak > 0);
        assert!(error
            .to_string()
            .contains("RSS guard exceeded during test_hot_loop"));
    }

    /// The periodic hot-loop guard must fire from inside the loop it samples,
    /// not only at stage boundaries: with a zero budget and interval 8, the
    /// first seven checkpoints of a resident process pass untouched and the
    /// eighth (the first actual sample) breaches.
    #[test]
    fn periodic_hot_loop_guard_breaches_at_its_sampling_interval() {
        let mut guard = search::genome_wide::PeriodicRssGuard::new(Some(0), 8);
        for iteration in 1..8 {
            guard
                .checkpoint("periodic_hot_loop")
                .expect("unsampled iterations never fire");
            assert_eq!(guard.peak_bytes(), 0, "iteration {iteration} must not sample");
        }
        let error = guard
            .checkpoint("periodic_hot_loop")
            .expect_err("sampled hot-loop checkpoint must fire");
        assert!(error
            .to_string()
            .contains("RSS guard exceeded during periodic_hot_loop"));
        assert!(guard.peak_bytes() > 0);
    }

    /// Unconditional sampling with a generous budget records the peak without
    /// firing, so stage boundaries and hot loops share one honest peak.
    #[test]
    fn periodic_guard_tracks_peak_without_false_breach() {
        let mut guard = search::genome_wide::PeriodicRssGuard::new(Some(u64::MAX), 1);
        for _ in 0..4 {
            guard.checkpoint("peak_tracking").expect("no breach under u64::MAX");
        }
        assert!(guard.peak_bytes() > 0);
    }
}
