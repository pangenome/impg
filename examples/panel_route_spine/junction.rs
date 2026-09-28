//! WITHIN-READ MEM ADJACENCY — the junction-evidence restriction (owner
//! model fix, 2026-09): "follow the actual junctions". A candidate
//! junction's seam evidence comes ONLY from reads that actually cross it.
//!
//! The L149 seam machinery predicts the features of read windows spanning a
//! candidate juxtaposition; the OLD model charged those predictions against
//! the POOLED routed shares of the adjacent partitions, so a novel junction
//! whose window is conserved yeast sequence was rewarded for observations
//! produced by unrelated reads elsewhere (the measured M1 chain defect:
//! novel junction charges predominantly negative, median about -2.2k).
//! Under the fix, a novel junction's predicted seam profile is charged
//! against observed evidence contributed ONLY by the reads whose
//! within-read MEM adjacency chain realizes the juxtaposition: a record
//! placed on the left row's source lane reaching the exit flank window,
//! immediately followed (in read order) by a record placed on the right
//! row's source lane reaching the entry flank window, each record's
//! occurrence strand matching the corresponding row's orientation — or,
//! for a sequence-continuous junction (same source, same orientation,
//! gap 0), a single record whose one placement spans the junction point
//! (overlaps both flank windows): a read over a continuous boundary
//! matches the lane as one maximal MEM crossing it, with no within-read
//! adjacency at the boundary, so the pair form alone can never see it
//! (the census defect this corrects).
//!
//! The evidence is derived PER RUN from the base reads (`reads.fastq.gz`)
//! in the same streaming pass that re-derives the sample records (zero
//! cache discipline: nothing new is persisted; the compact index is rebuilt
//! each run and its size is reported). Co-occurring (real panel) junctions
//! keep the pooled charge bit-identically; only novel junctions are
//! restricted (see the callers).
//!
//! INDEX DESIGN (measured on the first build attempt: the naive
//! placement-pair event materialization reached 605M events / ~40GB from
//! just 1.12M placements, because multi-copy records cross-multiply):
//! placements are stored ONCE per canonical record (one-sided), with a
//! per-path start-sorted index and per-record read instances; crossing
//! queries walk one left row's window once and then filter per right row —
//! the cross product is never materialized.
//!
//! Crossing semantics (documented choices, no tuning constants):
//! - "compatible offset" is the flank-window reach: the left record's
//!   verified occurrence interval must overlap the exit flank window
//!   ([end-take, end) forward / [start, start+take) reverse, take = the
//!   segment length capped at L149, mirroring `segment_flank`'s crop), and
//!   the right record's the entry flank window ([start, start+take) forward
//!   / [end-take, end) reverse). The read's inter-record gap is NOT required
//!   to equal the candidate juxtaposition gap (every real junction has an
//!   unmatchable seam of a few bases; equality would reject the truth's own
//!   reads) — the gap distribution of crossing events is REPORTED so the
//!   choice is measurable.
//! - A read pays at most once per junction composition, with the
//!   equal-share routed mass of its crossing records at the junction's
//!   adjacent partitions (the same routing machinery as the pooled shares,
//!   restricted to the crossing reads).
//! - Crossing evidence pools by seam COMPOSITION (the byte-level junction
//!   identity the seam machinery already interns), matching the prediction
//!   pooling; per-pair crossing counts are reported alongside.

use crate::*;

// ---------------------------------------------------------------------------
// Per-read chains and the record table.
// ---------------------------------------------------------------------------

/// One record's placement inside one read (the within-read MEM chain).
#[derive(Clone, Copy)]
pub struct ChainEntry {
    /// Canonical record id (the orientation-canonical token vector's id).
    pub record: u32,
    /// The instance's walk FORM relative to the record's canonical walk:
    /// true when the instance's walk is the mirrored (reverse-complement)
    /// encoding. Node signs in a walk are the panel's per-k-mer storage
    /// orientations (a path's own steps carry mixed signs), so the form is
    /// NOT a match strand; the match strand of a (chain entry, placement)
    /// pair is resolved against the placement's orientation at query time
    /// (see `JunctionSpanIndex`).
    pub mirrored: bool,
    /// Read coordinates of the record's extent [start, end).
    pub read_start: u32,
    pub read_end: u32,
}

/// One CROSSING-RECORD-ATTESTED composition (ITEM 2's expressibility
/// bound, owner ruling 2026-09-28): the junction between `left_source`'s
/// material exiting at `left_exit` and `right_source`'s entering at
/// `right_entry`, attested by `attesting_reads` distinct observed reads
/// (each read's consecutive record pair placing anchor mass on both
/// sides). The orientation flags are the placement forwardness per side
/// (a false-forward side reads its source right-to-left).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct AttestedComposition {
    pub left_source: usize,
    pub right_source: usize,
    pub left_forward: bool,
    pub right_forward: bool,
    /// The reads attesting this source pair's junction (existence, not
    /// magnitude — the placement-ambiguity cross-product is AGGREGATED
    /// into this one entry, never materialized: the junction is the
    /// source pair, and the exits/entries are the read-side anchor-point
    /// ranges over the attesting placement pairs).
    pub attesting_reads: u64,
    /// The left-side exit-point range over the attesting pairs (the left
    /// record's anchored extent ends; the true left cut sits in
    /// [left_exit_min, left_exit_max + gap_max]).
    pub left_exit_min: u64,
    pub left_exit_max: u64,
    /// The right-side entry-point range (the right record's anchored
    /// extents start; the true right cut sits in
    /// [right_entry_min - gap_max, right_entry_max]).
    pub right_entry_min: u64,
    pub right_entry_max: u64,
    /// The WIDEST inter-record read gap over the attesting pairs (the
    /// read's own gap between the left record's end and the right's
    /// start — the bracket's constant-free enclosure bound).
    pub gap_max: i64,
    /// The window row's locus covering the side's role points when the
    /// side is window-anchored (None = the side's points fall in no
    /// window row of the run's loci). At least one side is always
    /// anchored (the census's own bound); the materialization's locus
    /// rule reads the pair.
    pub left_window_locus: Option<usize>,
    pub right_window_locus: Option<usize>,
}

/// The per-run read-adjacency derivation: per read, the position-ordered
/// record chain, over a canonical-record table. Verifies against the sample
/// index's stats exactly like the record re-derivation it replaces (reads,
/// MEM records, distinct MEMs).
pub struct ReadChains {
    pub chains: Vec<Vec<ChainEntry>>,
    /// Canonical record tokens, by id.
    pub record_tokens: Vec<Vec<u64>>,
    /// Multiplicity per record id (read incidences).
    pub record_counts: Vec<u64>,
    /// One representative tagged walk per record id, ALWAYS in the record's
    /// canonical form (mirrored instances are reverse-complemented at
    /// interning): positions are those of the interning read, and the
    /// placement geometry is shared because the walk's relative gaps are.
    pub record_walks: Vec<Vec<(i32, u64)>>,
}

/// Derive the per-read chains and canonical record multiset in one pass over
/// the base reads (the process-wide anchor scheme; see
/// `derive_read_chains_scheme` for the scheme-conditioned guards).
pub fn derive_read_chains(
    panel: &SyngIndex,
    reads: &[Vec<u8>],
    expected: &sample::SampleStats,
    rss: &mut genome::PeriodicRssGuard,
) -> io::Result<ReadChains> {
    derive_read_chains_scheme(panel, reads, expected, rss, mem_records::anchor_scheme())
}

/// `derive_read_chains` with an explicit anchor scheme. Under `Syng` the
/// re-derived record multiset is pinned to the stored sample index's stats
/// (the re-derivation consistency guards). Under `Canonical` the record
/// set legitimately differs (the canonical qualification changes anchor
/// positions, hence MEM boundaries and multiplicities — the scheme of
/// record for the collapsed index), so the stored-scheme totals are not a
/// valid guard; only the reads-count identity and internal self-consistency
/// hold (the sample index remains base data for the read-length histogram,
/// and the stored-count pooled diagnostics are quarantined as non-gate).
pub fn derive_read_chains_scheme(
    panel: &SyngIndex,
    reads: &[Vec<u8>],
    expected: &sample::SampleStats,
    rss: &mut genome::PeriodicRssGuard,
    scheme: mem_records::AnchorScheme,
) -> io::Result<ReadChains> {
    ensure(reads.len() as u64 == expected.reads, "reads count mismatch")?;
    let mut record_index: BTreeMap<Vec<u64>, u32> = BTreeMap::new();
    let mut record_tokens: Vec<Vec<u64>> = Vec::new();
    let mut record_counts: Vec<u64> = Vec::new();
    let mut record_walks: Vec<Vec<(i32, u64)>> = Vec::new();
    let mut chains: Vec<Vec<ChainEntry>> = Vec::with_capacity(reads.len());
    let k = panel.syncmer_length_bp() as u64;
    for chunk in reads.chunks(4096) {
        // Parallel per-read extraction (positions absolute in the read).
        let per_read: Vec<(Vec<(Vec<u64>, bool, u32, u32)>, Vec<Vec<(i32, u64)>>)> = chunk
            .par_iter()
            .map(|read| {
                let tagged = mem_records::tagged_mem_records_scheme(panel, read, scheme)?;
                let mut entries: Vec<(Vec<u64>, bool, u32, u32)> =
                    Vec::with_capacity(tagged.len());
                let mut walks: Vec<Vec<(i32, u64)>> = Vec::with_capacity(tagged.len());
                for walk in &tagged {
                    let encoded = encode_walk(walk)?;
                    let tokens = canonical(&encoded);
                    // The instance's form: mirrored iff its own encoding is
                    // the non-canonical orientation of the record.
                    let mirrored = encoded != tokens;
                    let start = walk.first().expect("record").1 as u32;
                    let end = (walk.last().expect("record").1 + k) as u32;
                    entries.push((tokens, mirrored, start, end));
                    walks.push(walk.clone());
                }
                // The within-read chain: position order (the adjacency
                // order the junction test uses), ties broken by extent and
                // content for determinism.
                let mut order: Vec<usize> = (0..entries.len()).collect();
                order.sort_by(|&a, &b| {
                    (entries[a].2, entries[a].3, &entries[a].0).cmp(&(entries[b].2, entries[b].3, &entries[b].0))
                });
                Ok((
                    order.iter().map(|&i| entries[i].clone()).collect(),
                    order.iter().map(|&i| walks[i].clone()).collect(),
                ))
            })
            .collect::<io::Result<_>>()?;
        for (entries, walks) in per_read {
            let mut chain: Vec<ChainEntry> = Vec::with_capacity(entries.len());
            for ((tokens, mirrored, read_start, read_end), walk) in
                entries.into_iter().zip(walks.iter())
            {
                let id = match record_index.get(&tokens) {
                    Some(&id) => id,
                    None => {
                        let id = record_tokens.len() as u32;
                        record_index.insert(tokens.clone(), id);
                        record_tokens.push(tokens);
                        record_counts.push(0);
                        // Store the canonical-form walk: mirrored instances
                        // are reverse-complemented so every record's walk
                        // has one fixed orientation.
                        record_walks.push(if mirrored {
                            reverse_complement_walk(walk, k)
                        } else {
                            walk.clone()
                        });
                        id
                    }
                };
                record_counts[id as usize] += 1;
                chain.push(ChainEntry {
                    record: id,
                    mirrored,
                    read_start,
                    read_end,
                });
            }
            chains.push(chain);
        }
        if chains.len() % (4096 * 8) < 4096 {
            rss_probe(rss, "junction_read_chains")?;
        }
    }
    let total: u64 = record_counts.iter().sum();
    if scheme == mem_records::AnchorScheme::Syng {
        ensure(total == expected.mem_records, "re-derived MEM record total mismatch")?;
        ensure(
            record_tokens.len() == expected.distinct_mems,
            "re-derived distinct MEM mismatch",
        )?;
    }
    Ok(ReadChains {
        chains,
        record_tokens,
        record_counts,
        record_walks,
    })
}

// ---------------------------------------------------------------------------
// The span index: one-sided record placements + crossing queries.
// ---------------------------------------------------------------------------

/// Per canonical record: the partitions its routed occurrences touch and
/// its feature multiset (subwalks), for the crossing reads' payment.
struct SpanRecord {
    touched: Vec<u32>,
    features: Vec<FeatureKey>,
}

pub struct JunctionSpanIndex {
    /// The CROSSING-RECORD-ATTESTED composition set (ITEM 2's expressibility
    /// bound, owner ruling 2026-09-28): every (left source, left exit
    /// point, right source, right entry point, forwardness) composition
    /// attested by at least one observed record — a read's CONSECUTIVE
    /// record pair whose placements put the left record's anchor mass on
    /// the left source and the right record's on the right source, the
    /// pair straddling the junction (the record-once adjacency the crossing
    /// queries price). Existence, not magnitude — constant-free and tight
    /// by construction (the bound is crossing records x per-side placement
    /// ambiguity). The `left_forward`/`right_forward` flags are the
    /// placement orientations (the record's walk form matching its
    /// placement): a false-forward side reads its source right-to-left.
    pub(in crate) attested: Vec<AttestedComposition>,
    /// Per record: verified in-window placements `(path, start, mirrored)`,
    /// sorted by (path, start, mirrored). `mirrored` is the orientation the
    /// placement was verified in: the record's canonical walk (false) or
    /// its reverse complement (true). The placement's bp interval is
    /// `[start + extents[record][mirrored].0, start + extents[record][mirrored].1)`.
    placements: Vec<Vec<(u32, u32, bool)>>,
    /// Per record: the anchor bp extent offsets of its canonical walk
    /// (index 0) and of its mirrored walk (index 1).
    extents: Vec<[(u32, u32); 2]>,
    /// Per path: `(start, record, mirrored)` sorted by (start, record).
    index: HashMap<u32, Vec<(u32, u32, bool)>>,
    /// Per record: read instances (read, chain position).
    record_reads: Vec<Vec<(u32, u32)>>,
    /// The within-read chains (right-side adjacency lookups and gap stats).
    chains: Vec<Vec<ChainEntry>>,
    /// Payment data per placed record (eager: only the placed-record set,
    /// which the first build measured at ~9k of 776k records).
    records: Vec<SpanRecord>,
    pub stats: serde_json::Value,
}

/// The bp extent offsets (lo, hi) of a walk: the interval `[lo, hi)` spans
/// the walk's anchors plus the syncmer length.
fn walk_extent(walk: &[(i32, u64)], k: u64) -> (u32, u32) {
    let mut lo = u64::MAX;
    let mut hi = 0u64;
    for &(_, rel) in walk {
        lo = lo.min(rel);
        hi = hi.max(rel);
    }
    (lo as u32, (hi + k) as u32)
}

/// Candidate placements of one record: every anchor's territory entries,
/// deduplicated before verification (each distinct placement verified once
/// — the same placement is derivable through several of its anchors).
fn candidate_placements(
    path_walk: &[(i32, u64)],
    territory: &TerritoryIndex,
) -> HashSet<(usize, u64)> {
    let mut candidates: HashSet<(usize, u64)> = HashSet::new();
    for anchor in 0..path_walk.len() {
        let mut seen: HashSet<(usize, u64)> = HashSet::new();
        occurrences_from_anchor(path_walk, anchor, territory, &mut seen);
        candidates.extend(seen);
    }
    candidates
}

/// Build the per-run junction span index from the read chains.
#[allow(clippy::too_many_arguments)]
pub fn build_junction_span_index(
    territory: &TerritoryIndex,
    sources: &routes::Sources,
    source_of_path: &[usize],
    reads: &[Vec<u8>],
    chains: &ReadChains,
    // The run's per-locus window rows (the extension's own anchor bound:
    // a cross-source composition matters only when a side's anchor point
    // lands in a window row of one of the run's loci — the per-locus
    // domain being extended; constant-free, the window structure IS the
    // bound).
    window_rows: &[Vec<SourceRange>],
    // Per canonical record tokens: the partitions its routed occurrences
    // touch (from the run's routing pass, min-anchor convention, exactly
    // the shares the pooled observations use).
    touched_by_tokens: &HashMap<Vec<u64>, BTreeSet<u32>>,
    k: u64,
    rss: &mut genome::PeriodicRssGuard,
) -> io::Result<JunctionSpanIndex> {
    let started = Instant::now();
    // Placements + extents per canonical record, parallel over records.
    // A record's content can occur on the lanes in EITHER orientation —
    // node signs are per-k-mer storage orientations, so BOTH the canonical
    // walk and its reverse complement are searched (mirroring
    // `route_record`'s two-orientation occurrence search): a placement
    // verified in the canonical form is a content-forward occurrence (its
    // read instances with canonical-form walks matched the lane forward);
    // a placement verified in the mirrored form is a reverse-strand
    // occurrence. Only records with at least one in-window anchor (in
    // either orientation) are searched.
    let per_record: Vec<(Vec<(u32, u32, bool)>, [(u32, u32); 2])> = (0..chains
        .record_tokens.len())
        .into_par_iter()
        .map(|id| {
            let walk = &chains.record_walks[id];
            let mirrored_walk = reverse_complement_walk(walk, k);
            let extent = [walk_extent(walk, k), walk_extent(&mirrored_walk, k)];
            let has_entries = walk
                .iter()
                .chain(mirrored_walk.iter())
                .any(|&(node, _)| territory.node_range(node).1 > 0);
            let mut placements: Vec<(u32, u32, bool)> = Vec::new();
            if has_entries {
                for (path_walk, mirrored) in [(walk, false), (&mirrored_walk, true)] {
                    for (path, start) in candidate_placements(path_walk, territory) {
                        placements.push((path as u32, start as u32, mirrored));
                    }
                }
                placements.sort_unstable();
                placements.dedup();
            }
            Ok::<_, io::Error>((placements, extent))
        })
        .collect::<io::Result<_>>()?;
    let mut placements: Vec<Vec<(u32, u32, bool)>> = Vec::with_capacity(per_record.len());
    let mut extents: Vec<[(u32, u32); 2]> = Vec::with_capacity(per_record.len());
    for (list, extent) in per_record {
        placements.push(list);
        extents.push(extent);
    }
    let records_with_placements = placements.iter().filter(|p| !p.is_empty()).count();
    let total_placements: u64 = placements.iter().map(|p| p.len() as u64).sum();
    rss_probe(rss, "junction_span_placements")?;

    // Per-record read instances (chain positions).
    let mut record_reads: Vec<Vec<(u32, u32)>> = vec![Vec::new(); chains.record_tokens.len()];
    for (read, chain) in chains.chains.iter().enumerate() {
        for (position, entry) in chain.iter().enumerate() {
            record_reads[entry.record as usize].push((read as u32, position as u32));
        }
    }

    // Per-path start index over all placements.
    let mut index: HashMap<u32, Vec<(u32, u32, bool)>> = HashMap::new();
    for (record, list) in placements.iter().enumerate() {
        for &(path, start, mirrored) in list {
            index
                .entry(path)
                .or_default()
                .push((start, record as u32, mirrored));
        }
    }
    let index_entries: u64 = index.values().map(|list| list.len() as u64).sum();
    for list in index.values_mut() {
        list.sort_unstable();
    }
    rss_probe(rss, "junction_span_index")?;

    // Payment data for every PLACED record (the crossing reads pay with
    // these records' features and routed shares).
    let mut records: Vec<SpanRecord> = Vec::with_capacity(chains.record_tokens.len());
    let mut feature_memo: HashMap<Vec<u64>, Vec<FeatureKey>> = HashMap::new();
    let mut paid_feature_incidences = 0u64;
    for id in 0..chains.record_tokens.len() {
        if placements[id].is_empty() {
            records.push(SpanRecord {
                touched: Vec::new(),
                features: Vec::new(),
            });
            continue;
        }
        let tokens = &chains.record_tokens[id];
        let touched: Vec<u32> = touched_by_tokens
            .get(tokens)
            .map(|set| set.iter().copied().collect())
            .unwrap_or_default();
        let features = match feature_memo.get(tokens) {
            Some(features) => features.clone(),
            None => {
                let features = enumerate_subwalks(tokens)?;
                feature_memo.insert(tokens.clone(), features.clone());
                features
            }
        };
        paid_feature_incidences += features.len() as u64;
        records.push(SpanRecord {
            touched,
            features,
        });
    }
    drop(feature_memo);

    // Spot-checks vs the raw reads (validation rung a): for sampled
    // placements, the read's bytes at the record's read extent must equal
    // the source-lane bytes at the verified placement — directly when the
    // instance's walk form matches the placement's orientation (a
    // forward lane match), reverse-complemented otherwise.
    let mut spot_checks: Vec<serde_json::Value> = Vec::new();
    let mut spot_source = 0usize;
    'outer: for (record, list) in placements.iter().enumerate() {
        for &(path, start, mirrored) in list {
            if spot_checks.len() >= 8 {
                break 'outer;
            }
            spot_source += 1;
            if spot_source % 97 != 0 {
                continue;
            }
            let Some(&(read, position)) = record_reads[record].first() else {
                continue;
            };
            let entry = &chains.chains[read as usize][position as usize];
            let read_seq = &reads[read as usize];
            let source = source_of_path[path as usize];
            let lane_length = sources.lanes[source].1;
            let (extent_lo, extent_hi) = extents[record][usize::from(mirrored)];
            let interval_lo = start as u64 + extent_lo as u64;
            let interval_hi = start as u64 + extent_hi as u64;
            let forward = entry.mirrored == mirrored;
            let in_bounds = interval_hi <= lane_length
                && (entry.read_end as usize) <= read_seq.len()
                && entry.read_start < entry.read_end;
            let verified = if in_bounds {
                let expected = sources.fetch(source, interval_lo, interval_hi)?;
                let slice = &read_seq[entry.read_start as usize..entry.read_end as usize];
                if forward {
                    slice == expected.as_slice()
                } else {
                    impg::graph::reverse_complement(slice) == expected
                }
            } else {
                false
            };
            spot_checks.push(serde_json::json!({
                "read": read,
                "record": record,
                "path": path,
                "interval": [interval_lo, interval_hi],
                "forward_lane_match": forward,
                "read_interval": [entry.read_start, entry.read_end],
                "verified": verified,
            }));
        }
    }
    let spot_checked = spot_checks.len() as u64;
    let spot_failures = spot_checks
        .iter()
        .filter(|row| !row["verified"].as_bool().unwrap_or(false))
        .count() as u64;

    let index_bytes = (placements.len() * 8) as u64
        + total_placements * 12
        + index_entries * 12
        + record_reads.iter().map(|l| l.len() as u64 * 8).sum::<u64>()
        + paid_feature_incidences * 24;
    let mut stats = serde_json::json!({
        "reads": chains.chains.len(),
        "distinct_records": chains.record_tokens.len(),
        "records_with_placements": records_with_placements,
        "placements": total_placements,
        "index_entries": index_entries,
        "index_bytes": index_bytes,
        "paid_feature_incidences": paid_feature_incidences,
        "spot_checks": {
            "checked": spot_checked,
            "failures": spot_failures,
            "rows": spot_checks,
        },
        "wall_seconds": started.elapsed().as_secs_f64(),
    });
    eprintln!(
        "[junction] span index: {} reads, {} records ({} placed, {} placements, \
         {} index entries) ({:.1}s)",
        chains.chains.len(),
        chains.record_tokens.len(),
        records_with_placements,
        total_placements,
        index_entries,
        started.elapsed().as_secs_f64()
    );
    // The CROSSING-RECORD-ATTESTED composition set (the expressibility
    // bound's census): every cross-source junction attested by a read's
    // consecutive record pair — the left record placing on one source, the
    // right on another (the same-source novel junctions are the rung-(i)
    // co-occurrence case, not attested here). Existence, not magnitude;
    // the per-composition read count is kept for the report (the bound is
    // crossing records x per-side placement ambiguity).
    let attested_started = Instant::now();
    let mut attested_map: std::collections::BTreeMap<
        (usize, usize, bool, bool),
        (
            u64,
            i64,
            u64,
            u64,
            u64,
            u64,
            Option<usize>,
            Option<usize>,
        ),
    > = std::collections::BTreeMap::new();
    // The routed-attestation filter's per-record touched sets (the
    // routing pass's own min-anchor disambiguation — the same convention
    // the observed side charges with; owner ruling 2026-09-28): a
    // placement is ROUTED-ATTRIBUTED iff the partition covering its
    // anchor extent is in its record's routed touched set. Precomputed
    // per record id once.
    let routed: Vec<Option<&std::collections::BTreeSet<u32>>> = (0..chains
        .record_tokens
        .len())
        .map(|id| touched_by_tokens.get(&chains.record_tokens[id]))
        .collect();
    // The partition covering a placement's anchor extent on its path
    // (the interval containing the extent's low end; the territory's own
    // per-path interval lists).
    let covering_partition = |path: u32, start: u32, lo: u32| -> Option<u32> {
        let point = start as u64 + lo as u64;
        territory
            .path_intervals
            .get(path as usize)?
            .iter()
            .find(|&&(bs, be, _)| bs <= point && point < be)
            .map(|&(_, _, partition)| partition)
    };
    // The distinct consecutive (record, form) instances with their
    // per-gap read counts: the placement enumeration does not depend on
    // the read, so it runs once per distinct form pair; the per-read
    // geometry (the gap, hence the co-linear offset) is what the
    // continuity test needs per instance.
    let mut pair_instances: std::collections::BTreeMap<
        (u32, u32, bool, bool),
        std::collections::BTreeMap<i64, u64>,
    > = std::collections::BTreeMap::new();
    for chain in chains.chains.iter() {
        for pair in chain.windows(2) {
            let (left_entry, right_entry) = (&pair[0], &pair[1]);
            if placements[left_entry.record as usize].is_empty()
                || placements[right_entry.record as usize].is_empty()
            {
                continue;
            }
            if routed[left_entry.record as usize].is_none()
                || routed[right_entry.record as usize].is_none()
            {
                continue;
            }
            let gap = right_entry.read_start as i64 - left_entry.read_end as i64;
            pair_instances
                .entry((
                    left_entry.record,
                    right_entry.record,
                    left_entry.mirrored,
                    right_entry.mirrored,
                ))
                .or_default()
                .entry(gap)
                .and_modify(|count| *count += 1)
                .or_insert(1);
        }
    }
    rss_probe(rss, "junction_attested_pairs")?;
    // Per record: the ROUTED placements (the covering partition is in
    // the record's routed touched set) with the WINDOW-ANCHOR bits for
    // both role points. A placement with no window anchor can still
    // pair with an anchored partner on the other side, so both are kept
    // and the OR test applies at pair time.
    let mut window_by_source: std::collections::HashMap<usize, Vec<(u64, u64, usize)>> =
        std::collections::HashMap::new();
    for (locus, locus_rows) in window_rows.iter().enumerate() {
        for row in locus_rows {
            if row.reverse || row.start >= row.end {
                continue;
            }
            window_by_source
                .entry(row.source)
                .or_default()
                .push((row.start, row.end, locus));
        }
    }
    for rows in window_by_source.values_mut() {
        rows.sort_unstable();
    }
    let window_locus_at = |source: usize, point: u64| -> Option<usize> {
        window_by_source
            .get(&source)?
            .iter()
            .rev()
            .find(|&&(start, end, _)| start <= point && point < end)
            .map(|&(_, _, locus)| locus)
    };
    struct RoutedPlacement {
        start: u32,
        mirrored: bool,
        source: usize,
        lo_locus: Option<usize>,
        hi_locus: Option<usize>,
    }
    let routed_placements: Vec<Vec<RoutedPlacement>> = (0..chains.record_tokens.len())
        .into_par_iter()
        .map(|id| {
            let touched = match routed[id] {
                Some(touched) => touched,
                None => return Vec::new(),
            };
            let mut out: Vec<RoutedPlacement> = Vec::new();
            for &(path, start, mirrored) in &placements[id] {
                let source = source_of_path[path as usize];
                if source == usize::MAX {
                    continue;
                }
                let (lo, hi) = extents[id][usize::from(mirrored)];
                let Some(partition) = covering_partition(path, start, lo) else {
                    continue;
                };
                if !touched.contains(&partition) {
                    continue;
                }
                let lo_locus = window_locus_at(source, start as u64 + lo as u64);
                let hi_locus = window_locus_at(source, start as u64 + hi as u64);
                out.push(RoutedPlacement {
                    start,
                    mirrored,
                    source,
                    lo_locus,
                    hi_locus,
                });
            }
            out
        })
        .collect();
    rss_probe(rss, "junction_routed_placements")?;
    let mut continuity_killed_reads = 0u64;
    let mut novel_pairs_enumerated = 0u64;
    let mut unanchored_pair_skips = 0u64;
    for ((left_record, right_record, left_mirrored, right_mirrored), gaps) in pair_instances {
        let left_list = &placements[left_record as usize];
        let right_list = &placements[right_record as usize];
        let left_routed = &routed_placements[left_record as usize];
        let right_routed = &routed_placements[right_record as usize];
        // The SAME-SOURCE CONTINUITY test per observed gap (the read's
        // own co-linear offset): a same-source placement pair carrying
        // both records at that offset is the homologous continuation —
        // those reads attest no junction. The right placements are
        // looked up per path (binary search over the sorted list), so
        // the scan is linear in the left placements, not their product.
        let left_extent = extents[left_record as usize][usize::from(left_mirrored)];
        let mut surviving: Vec<(i64, u64)> = Vec::new();
        let mut all_continuation = true;
        for (&gap, &reads) in gaps.iter() {
            let offset = gap + (left_extent.1 - left_extent.0) as i64;
            let mut continuation = false;
            'pair_scan: for &(left_path, left_start, _) in left_list {
                let lo = right_list.partition_point(|&(p, _, _)| p < left_path);
                let hi = right_list.partition_point(|&(p, _, _)| p <= left_path);
                for &(right_path, right_start, _) in &right_list[lo..hi] {
                    let _ = right_path;
                    if right_start as i64 - left_start as i64 == offset {
                        continuation = true;
                        break 'pair_scan;
                    }
                }
            }
            if continuation {
                continuity_killed_reads += reads;
            } else {
                all_continuation = false;
                surviving.push((gap, reads));
            }
        }
        if all_continuation {
            continue;
        }
        novel_pairs_enumerated += 1;
        let reads_total: u64 = surviving.iter().map(|&(_, r)| r).sum();
        let gap_max = surviving
            .iter()
            .map(|&(g, _)| g)
            .max()
            .unwrap_or(i64::MIN);
        // The per-side, per-(source, forwardness) aggregation of the
        // routed placements: the placement-ambiguity cross-product is
        // AGGREGATED (min/max of the role points), never materialized.
        // The anchor locus is the window row covering the side's ROLE
        // point (the exit for the left role, the entry for the right).
        let mut lefts: std::collections::BTreeMap<
            (usize, bool),
            (u64, u64, Option<usize>),
        > = std::collections::BTreeMap::new();
        for left in left_routed {
            let left_forward = left_mirrored == left.mirrored;
            let (lo, hi) = extents[left_record as usize][usize::from(left.mirrored)];
            let exit = if left_forward {
                left.start as u64 + hi as u64
            } else {
                left.start as u64 + lo as u64
            };
            let anchor_locus = if left_forward {
                left.hi_locus
            } else {
                left.lo_locus
            };
            let entry = lefts
                .entry((left.source, left_forward))
                .or_insert((exit, exit, anchor_locus));
            entry.0 = entry.0.min(exit);
            entry.1 = entry.1.max(exit);
        }
        let mut rights: std::collections::BTreeMap<
            (usize, bool),
            (u64, u64, Option<usize>),
        > = std::collections::BTreeMap::new();
        for right in right_routed {
            let right_forward = right_mirrored == right.mirrored;
            let (rlo, rhi) = extents[right_record as usize][usize::from(right.mirrored)];
            let entry_point = if right_forward {
                right.start as u64 + rlo as u64
            } else {
                right.start as u64 + rhi as u64
            };
            let anchor_locus = if right_forward {
                right.lo_locus
            } else {
                right.hi_locus
            };
            let entry = rights
                .entry((right.source, right_forward))
                .or_insert((entry_point, entry_point, anchor_locus));
            entry.0 = entry.0.min(entry_point);
            entry.1 = entry.1.max(entry_point);
        }
        // The cross-source source-pair groups with the window-anchor OR:
        // at least one side's role points land in a window row of the
        // run's loci (the per-locus domain being extended — the
        // constant-free bound; the unanchored pairs carry no locus and
        // cannot be materialized anywhere).
        for (&(left_source, left_forward), &(lexit_min, lexit_max, llocus)) in &lefts {
            for (&(right_source, right_forward), &(rentry_min, rentry_max, rlocus)) in &rights {
                if right_source == left_source {
                    continue;
                }
                if llocus.is_none() && rlocus.is_none() {
                    unanchored_pair_skips += 1;
                    continue;
                }
                let key = (left_source, right_source, left_forward, right_forward);
                let entry = attested_map
                    .entry(key)
                    .or_insert((0, i64::MIN, lexit_min, lexit_max, rentry_min, rentry_max, None, None));
                entry.0 += reads_total;
                entry.1 = entry.1.max(gap_max);
                entry.2 = entry.2.min(lexit_min);
                entry.3 = entry.3.max(lexit_max);
                entry.4 = entry.4.min(rentry_min);
                entry.5 = entry.5.max(rentry_max);
                if llocus.is_some() {
                    entry.6 = llocus;
                }
                if rlocus.is_some() {
                    entry.7 = rlocus;
                }
            }
        }
    }
    rss_probe(rss, "junction_attested_compositions")?;
    let attested: Vec<AttestedComposition> = attested_map
        .into_iter()
        .map(
            |(
                (left_source, right_source, left_forward, right_forward),
                (
                    attesting_reads,
                    gap_max,
                    left_exit_min,
                    left_exit_max,
                    right_entry_min,
                    right_entry_max,
                    left_window_locus,
                    right_window_locus,
                ),
            )| AttestedComposition {
                left_source,
                right_source,
                left_forward,
                right_forward,
                attesting_reads,
                left_exit_min,
                left_exit_max,
                right_entry_min,
                right_entry_max,
                gap_max,
                left_window_locus,
                right_window_locus,
            },
        )
        .collect();
    let attested_seconds = attested_started.elapsed().as_secs_f64();
    eprintln!(
        "[junction] crossing-attested compositions: {} (cross-source, read-attested; \
         {attested_seconds:.1}s; novel pairs {}, continuity-killed reads {}, \
         unanchored skips {})",
        attested.len(),
        novel_pairs_enumerated,
        continuity_killed_reads,
        unanchored_pair_skips
    );
    stats["attested_novel_pairs"] = novel_pairs_enumerated.into();
    stats["attested_continuity_killed_reads"] = continuity_killed_reads.into();
    stats["attested_unanchored_skips"] = unanchored_pair_skips.into();
    Ok(JunctionSpanIndex {
        placements,
        extents,
        index,
        record_reads,
        chains: chains.chains.clone(),
        records,
        attested,
        stats,
    })
}

// ---------------------------------------------------------------------------
// Crossing queries and the restricted junction charge.
// ---------------------------------------------------------------------------

/// The public same-source order-compatible gap between two segments (the
/// co-occurrence relation at segment granularity): `Some(gap >= 0)` when
/// the segments share source and orientation and are order-compatible in
/// that source's own coordinate frame (a real panel junction), `Some(negative)`
/// when they overlap, `None` when they cannot co-occur on any panel
/// chromosome. Mirrors `phasing::cooccurrence_gap` with the exit/entry
/// segments given directly.
pub fn segment_pair_gap(left: &SourceRange, right: &SourceRange) -> Option<i64> {
    if left.source != right.source || left.reverse != right.reverse {
        return None;
    }
    if left.reverse {
        Some(left.start as i64 - right.end as i64)
    } else {
        Some(right.start as i64 - left.end as i64)
    }
}

/// A bare (source, start, end, reverse) piece as a `SourceRange` for the
/// crossing query (partition/occurrence fields are not junction evidence).
pub fn junction_range(
    source: usize,
    start: u64,
    end: u64,
    reverse: bool,
    partition: u32,
) -> SourceRange {
    SourceRange {
        partition: partition as usize,
        occurrence: 0,
        source,
        start,
        end,
        reverse,
    }
}

/// The reads crossing one candidate junction: the matched (read, left
/// record, right record) adjacency realizations, the distinct-read count,
/// and the per-match inter-record read gaps (the "compatible offset"
/// measurement).
pub struct CrossingReads {
    /// (read, left record, right record) — one entry per matched adjacency
    /// realization (a read can realize a composition through more than one
    /// of its record pairs; a continuous single-record crossing has its
    /// left and right record equal).
    pub matches: Vec<(u32, u32, u32)>,
    pub distinct_reads: u64,
    pub events: u64,
    /// Inter-record read gaps of the PAIR realizations only (continuous
    /// single-record realizations have no inter-record gap).
    pub read_gaps: Vec<i64>,
}

/// The restricted charge of one junction composition.
pub struct RestrictedCharge {
    pub charge: f64,
    /// Distinct reads crossing any pair of the composition.
    pub spanning_reads: u64,
    pub events: u64,
    /// Predicted features with nonzero crossing-read evidence.
    pub span_feature_count: usize,
}

/// The maximal bp reach of a record's placement interval past its start
/// (records live in L150 reads; the margin covers the k-mer overhang and
/// the per-record extent offsets).
const SPAN_SCAN_SLACK: u64 = READ_LENGTH as u64 + 128;

/// Bound on the measured seam-shared context between two rows' lanes (the
/// single-record crossing forms only need context comparable to a record's
/// reach, which is bounded by the read length; larger sharing is
/// indistinguishable for L150 evidence).
const SHARED_CONTEXT_CAP: u64 = READ_LENGTH as u64 + 128;

impl JunctionSpanIndex {
    /// Exit flank window of a left segment, in its path's coordinates: the
    /// last `take` bases of the segment's material in molecule orientation,
    /// take = min(L149, segment length) — mirroring `segment_flank`'s crop
    /// (empty for zero-length segments, whose seam flank is empty too).
    fn exit_window(left: &SourceRange) -> (u64, u64) {
        let take = (READ_LENGTH as u64 - 1).min(left.end.saturating_sub(left.start));
        if left.reverse {
            (left.start, left.start + take)
        } else {
            (left.end.saturating_sub(take), left.end)
        }
    }

    /// Entry flank window of a right segment (its first `take` bases in
    /// molecule orientation).
    fn entry_window(right: &SourceRange) -> (u64, u64) {
        let take = (READ_LENGTH as u64 - 1).min(right.end.saturating_sub(right.start));
        if right.reverse {
            (right.end.saturating_sub(take), right.end)
        } else {
            (right.start, right.start + take)
        }
    }

    /// The read-adjacency candidates whose left side reaches `left`'s exit
    /// flank window on `left`'s own source lane with the row's orientation:
    /// `(read, chain position)` of the left-reaching record instances.
    /// Grouped by left row so a boundary matrix walks each left row's window
    /// once and filters per right row (the placement cross product is never
    /// materialized). Deduplicated per (read, position): a read pays at
    /// most once per junction composition, even when its record has
    /// several placements in the window.
    pub fn crossings_from(
        &self,
        left: &SourceRange,
        path_of_source: &[usize],
    ) -> Vec<(u32, u32)> {
        let window = Self::exit_window(left);
        let mut out: Vec<(u32, u32)> = Vec::new();
        if window.0 >= window.1 {
            return out;
        }
        let path = path_of_source[left.source] as u32;
        let Some(list) = self.index.get(&path) else {
            return out;
        };
        let lower = list.partition_point(|&(start, _, _)| {
            (start as u64) < window.0.saturating_sub(SPAN_SCAN_SLACK)
        });
        for &(start, record, mirrored) in &list[lower..] {
            if start as u64 >= window.1 {
                break;
            }
            let (_, hi) = self.extents[record as usize][usize::from(mirrored)];
            let interval_hi = start as u64 + hi as u64;
            if interval_hi <= window.0 {
                continue;
            }
            for &(read, position) in &self.record_reads[record as usize] {
                let entry = &self.chains[read as usize][position as usize];
                // The instance matches the lane forward iff its walk form
                // matches the placement's orientation; keep it iff that
                // match direction realizes the row's orientation.
                let forward = entry.mirrored == mirrored;
                if forward == left.reverse {
                    continue;
                }
                out.push((read, position));
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// The traversal point where `row`'s material exits into the junction
    /// (its last base in traversal order is the base before this point).
    fn traversal_exit_point(row: &SourceRange) -> u64 {
        if row.reverse {
            row.start
        } else {
            row.end
        }
    }

    /// The traversal point where `row`'s material enters from the junction.
    fn traversal_entry_point(row: &SourceRange) -> u64 {
        if row.reverse {
            row.end
        } else {
            row.start
        }
    }

    /// `row`'s lane bytes of traversal length `len` ENDING at the row's
    /// traversal boundary point `at`, in traversal orientation (a reverse
    /// row reads its lane right-to-left).
    fn traversal_bytes_before(
        sources: &routes::Sources,
        row: &SourceRange,
        at: u64,
        len: u64,
    ) -> io::Result<Vec<u8>> {
        let lane_length = sources.lanes[row.source].1;
        let (lo, hi) = if row.reverse {
            (at.min(lane_length), (at + len).min(lane_length))
        } else {
            (at.saturating_sub(len), at.min(lane_length))
        };
        if hi <= lo {
            return Ok(Vec::new());
        }
        let mut dna = sources.fetch(row.source, lo, hi)?;
        if row.reverse {
            dna = impg::graph::reverse_complement(&dna);
        }
        Ok(dna)
    }

    /// `row`'s lane bytes of traversal length `len` STARTING at the row's
    /// traversal boundary point `at`, in traversal orientation.
    fn traversal_bytes_after(
        sources: &routes::Sources,
        row: &SourceRange,
        at: u64,
        len: u64,
    ) -> io::Result<Vec<u8>> {
        let lane_length = sources.lanes[row.source].1;
        let (lo, hi) = if row.reverse {
            (at.saturating_sub(len), at.min(lane_length))
        } else {
            (at.min(lane_length), (at + len).min(lane_length))
        };
        if hi <= lo {
            return Ok(Vec::new());
        }
        let mut dna = sources.fetch(row.source, lo, hi)?;
        if row.reverse {
            dna = impg::graph::reverse_complement(&dna);
        }
        Ok(dna)
    }

    /// The maximal context length `L` (capped) such that the `L` traversal
    /// bytes ending at the RIGHT row's entry point on its own lane equal the
    /// `L` traversal bytes ending at the LEFT row's exit point on its lane:
    /// the length over which the lanes are sequence-continuous across the
    /// candidate junction's seam. A read whose single through-record is
    /// placed on the right lane realizes the juxtaposition exactly when its
    /// pre-seam reach fits within this shared context. Clamped to the left
    /// row's own material and both lanes' bounds.
    fn shared_entry_context(
        sources: &routes::Sources,
        left: &SourceRange,
        right: &SourceRange,
        x_l: u64,
        x_r: u64,
    ) -> io::Result<u64> {
        let mut cap = SHARED_CONTEXT_CAP.min(if left.reverse {
            left.end.saturating_sub(x_l)
        } else {
            x_l.saturating_sub(left.start)
        });
        cap = cap.min(sources.lanes[left.source].1).min(sources.lanes[right.source].1);
        while cap > 0 {
            let left_bytes = Self::traversal_bytes_before(sources, left, x_l, cap)?;
            let right_bytes = Self::traversal_bytes_before(sources, right, x_r, cap)?;
            if left_bytes == right_bytes {
                return Ok(cap);
            }
            // fall back to the longest common suffix of the fetched windows
            let mut shared = 0u64;
            while (shared as usize) < left_bytes.len().min(right_bytes.len())
                && left_bytes[left_bytes.len() - 1 - shared as usize]
                    == right_bytes[right_bytes.len() - 1 - shared as usize]
            {
                shared += 1;
            }
            return Ok(shared);
        }
        Ok(0)
    }

    /// The maximal context length `R` (capped) such that the `R` traversal
    /// bytes starting at the LEFT row's exit point on its lane equal the
    /// `R` traversal bytes starting at the RIGHT row's entry point on its
    /// lane (the mirror of `shared_entry_context`). Clamped to the right
    /// row's own material and both lanes' bounds.
    fn shared_exit_context(
        sources: &routes::Sources,
        left: &SourceRange,
        right: &SourceRange,
        x_l: u64,
        x_r: u64,
    ) -> io::Result<u64> {
        let mut cap = SHARED_CONTEXT_CAP.min(if right.reverse {
            x_r.saturating_sub(right.start)
        } else {
            right.end.saturating_sub(x_r)
        });
        cap = cap.min(sources.lanes[left.source].1).min(sources.lanes[right.source].1);
        while cap > 0 {
            let left_bytes = Self::traversal_bytes_after(sources, left, x_l, cap)?;
            let right_bytes = Self::traversal_bytes_after(sources, right, x_r, cap)?;
            if left_bytes == right_bytes {
                return Ok(cap);
            }
            let mut shared = 0u64;
            while (shared as usize) < left_bytes.len().min(right_bytes.len())
                && left_bytes[shared as usize] == right_bytes[shared as usize]
            {
                shared += 1;
            }
            return Ok(shared);
        }
        Ok(0)
    }

    /// Filter one left row's crossing candidates down to the reads that
    /// realize the junction into `right`, in three forms:
    ///
    /// - PAIR: the candidate's following chain record is placed in
    ///   `right`'s entry flank window on `right`'s own source lane, with
    ///   both records' occurrence orientations realizing the rows'
    ///   orientations. The generic form when the seam's lanes diverge
    ///   within the read (MEMs break at the seam).
    /// - RIGHT SINGLE: a record whose one placement on `right`'s lane
    ///   strictly covers the traversal entry point, with its pre-seam
    ///   traversal reach no longer than the lanes' measured shared entry
    ///   context. The generic form for a seam embedded in shared context
    ///   (the truth's mosaic junctions are cut at shared boundary words,
    ///   so the seam read matches one lane CONTINUOUSLY through the seam
    ///   and the other side's MEM is pruned as a subwalk of the longer
    ///   through-record): the read's pre-seam bytes equal the left row's
    ///   own flank bytes for exactly that shared context, so the read
    ///   realizes the juxtaposition; a read reaching past the shared
    ///   context matches the right lane's OWN continuation instead and
    ///   must not pay.
    /// - LEFT SINGLE: the mirror through-record on `left`'s lane strictly
    ///   covering the traversal exit point, with its post-seam traversal
    ///   reach no longer than the measured shared exit context. For a
    ///   same-source gap-0 junction both contexts are maximal (the lane is
    ///   continuous), so the single forms reproduce the
    ///   sequence-continuous crossing exactly; for a juxtaposition with
    ///   no shared context only the PAIR form can fire.
    pub fn filter_crossings(
        &self,
        crossings: &[(u32, u32)],
        left: &SourceRange,
        right: &SourceRange,
        sources: &routes::Sources,
        path_of_source: &[usize],
    ) -> io::Result<CrossingReads> {
        let window = Self::entry_window(right);
        let path = path_of_source[right.source] as u32;
        let x_l = Self::traversal_exit_point(left);
        let x_r = Self::traversal_entry_point(right);
        let entry_ctx = if window.0 < window.1 {
            Self::shared_entry_context(sources, left, right, x_l, x_r)?
        } else {
            0
        };
        let exit_ctx = Self::shared_exit_context(sources, left, right, x_l, x_r)?;
        let mut matches: Vec<(u32, u32, u32)> = Vec::new();
        let mut gaps: Vec<i64> = Vec::new();
        let mut events = 0u64;
        if window.0 < window.1 {
            for &(read, position) in crossings {
                let chain = &self.chains[read as usize];
                let entry = &chain[position as usize];
                // Pair case: the following chain record reaches the entry
                // flank window on the right row's path realizing the row's
                // orientation.
                if let Some(next) = chain.get(position as usize + 1) {
                    let list = &self.placements[next.record as usize];
                    let lower = list.partition_point(|&(p, start, _)| {
                        p < path
                            || (p == path
                                && (start as u64) < window.0.saturating_sub(SPAN_SCAN_SLACK))
                    });
                    for &(p, start, next_mirrored) in &list[lower..] {
                        if p > path || start as u64 >= window.1 {
                            break;
                        }
                        if p < path {
                            continue;
                        }
                        let next_forward = next.mirrored == next_mirrored;
                        if next_forward == right.reverse {
                            continue;
                        }
                        let (_, hi) = self.extents[next.record as usize]
                            [usize::from(next_mirrored)];
                        if start as u64 + hi as u64 <= window.0 {
                            continue;
                        }
                        events += 1;
                        gaps.push(next.read_start as i64 - entry.read_end as i64);
                        matches.push((read, entry.record, next.record));
                        break;
                    }
                }
                // Left single case: this record's own placement on the LEFT
                // row's lane strictly covers the traversal exit point, and
                // its post-seam traversal reach fits the shared exit
                // context (the read's post-seam bytes then equal the right
                // row's own entry flank bytes).
                if exit_ctx > 0 {
                    let list = &self.placements[entry.record as usize];
                    let lower = list.partition_point(|&(p, start, _)| {
                        p < path_of_source[left.source] as u32
                            || (p == path_of_source[left.source] as u32
                                && (start as u64)
                                    < x_l.saturating_sub(SPAN_SCAN_SLACK))
                    });
                    for &(p, start, self_mirrored) in &list[lower..] {
                        if p > path_of_source[left.source] as u32
                            || start as u64 >= x_l + SPAN_SCAN_SLACK
                        {
                            break;
                        }
                        if p != path_of_source[left.source] as u32 {
                            continue;
                        }
                        let (lo, hi) = self.extents[entry.record as usize]
                            [usize::from(self_mirrored)];
                        let interval_lo = start as u64 + lo as u64;
                        let interval_hi = start as u64 + hi as u64;
                        if !(interval_lo < x_l && interval_hi > x_l) {
                            continue;
                        }
                        let forward = entry.mirrored == self_mirrored;
                        if forward == left.reverse {
                            continue;
                        }
                        let reach = if left.reverse {
                            x_l - interval_lo
                        } else {
                            interval_hi - x_l
                        };
                        if reach == 0 || reach > exit_ctx {
                            continue;
                        }
                        events += 1;
                        matches.push((read, entry.record, entry.record));
                        break;
                    }
                }
            }
            // Right single case: records whose one placement on the RIGHT
            // row's lane strictly covers the traversal entry point, with the
            // pre-seam traversal reach within the shared entry context.
            if entry_ctx > 0 {
                let right_path = path_of_source[right.source] as u32;
                let empty: &[(u32, u32, bool)] = &[];
                let list = self
                    .index
                    .get(&right_path)
                    .map(|list| list.as_slice())
                    .unwrap_or(empty);
                let lower = list.partition_point(|&(start, _, _)| {
                    (start as u64) < x_r.saturating_sub(SPAN_SCAN_SLACK)
                });
                let mut seen: HashSet<(u32, u32)> = HashSet::new();
                for &(start, record, mirrored) in &list[lower..] {
                    if start as u64 >= x_r + SPAN_SCAN_SLACK {
                        break;
                    }
                    let (lo, hi) = self.extents[record as usize][usize::from(mirrored)];
                    let interval_lo = start as u64 + lo as u64;
                    let interval_hi = start as u64 + hi as u64;
                    if !(interval_lo < x_r && interval_hi > x_r) {
                        continue;
                    }
                    let reach = if right.reverse {
                        interval_hi - x_r
                    } else {
                        x_r - interval_lo
                    };
                    if reach == 0 || reach > entry_ctx {
                        continue;
                    }
                    for &(read, position) in &self.record_reads[record as usize] {
                        if !seen.insert((read, position)) {
                            continue;
                        }
                        let entry = &self.chains[read as usize][position as usize];
                        let forward = entry.mirrored == mirrored;
                        if forward == right.reverse {
                            continue;
                        }
                        events += 1;
                        matches.push((read, record, record));
                    }
                }
            }
        }
        let mut distinct: HashSet<u32> = HashSet::new();
        for &(read, _, _) in &matches {
            distinct.insert(read);
        }
        Ok(CrossingReads {
            matches,
            distinct_reads: distinct.len() as u64,
            events,
            read_gaps: gaps,
        })
    }

    /// The reads crossing the junction between `left` (exit segment) and
    /// `right` (entry segment): consecutive chain records placed on the
    /// rows' own source lanes, reaching the junction's flank windows, with
    /// occurrence strands matching the rows' orientations.
    pub fn crossing_reads(
        &self,
        left: &SourceRange,
        right: &SourceRange,
        sources: &routes::Sources,
        path_of_source: &[usize],
    ) -> io::Result<CrossingReads> {
        let candidates = self.crossings_from(left, path_of_source);
        self.filter_crossings(&candidates, left, right, sources, path_of_source)
    }

    /// Census diagnostic for specific reads: the full within-read chain
    /// with every record's placements (path, start, interval), exposing
    /// exactly what the crossing queries can see per read.
    pub fn trace_reads(&self, reads: &[u32]) -> Vec<serde_json::Value> {
        let mut out: Vec<serde_json::Value> = Vec::new();
        for &read in reads {
            let Some(chain) = self.chains.get(read as usize) else {
                out.push(serde_json::json!({"read": read, "missing": true}));
                continue;
            };
            let entries: Vec<serde_json::Value> = chain
                .iter()
                .map(|entry| {
                    let placements: Vec<serde_json::Value> = self.placements
                        [entry.record as usize]
                        .iter()
                        .map(|&(path, start, mirrored)| {
                            let (lo, hi) = self.extents[entry.record as usize]
                                [usize::from(mirrored)];
                            serde_json::json!([
                                path,
                                start,
                                mirrored,
                                start as u64 + lo as u64,
                                start as u64 + hi as u64,
                            ])
                        })
                        .collect();
                    serde_json::json!({
                        "record": entry.record,
                        "mirrored": entry.mirrored,
                        "read_interval": [entry.read_start, entry.read_end],
                        "placements": placements,
                    })
                })
                .collect();
            out.push(serde_json::json!({
                "read": read,
                "records": entries,
            }));
        }
        out
    }

    /// Census diagnostic for one candidate junction: every record with a
    /// placement on a row's own path reaching that row's flank window, with
    /// its placements on the path, its extent, and its read instances'
    /// chain contexts (previous/next records with read intervals). Bounded:
    /// one window scan per side, no cross products.
    pub fn trace_junction(
        &self,
        left: &SourceRange,
        right: &SourceRange,
        path_of_source: &[usize],
    ) -> serde_json::Value {
        let exit = Self::exit_window(left);
        let entry = Self::entry_window(right);
        serde_json::json!({
            "left": [left.source, left.start, left.end, left.reverse],
            "right": [right.source, right.start, right.end, right.reverse],
            "exit_window": [exit.0, exit.1],
            "entry_window": [entry.0, entry.1],
            "exit_side": self.trace_reaching(left, exit, path_of_source),
            "entry_side": self.trace_reaching(right, entry, path_of_source),
        })
    }

    /// Records whose placement on `row`'s own path overlaps `window`:
    /// `(record, start, interval, extent, all placements on the path, read
    /// instances with chain context)`.
    fn trace_reaching(
        &self,
        row: &SourceRange,
        window: (u64, u64),
        path_of_source: &[usize],
    ) -> Vec<serde_json::Value> {
        let mut out: Vec<serde_json::Value> = Vec::new();
        if window.0 >= window.1 {
            return out;
        }
        let path = path_of_source[row.source] as u32;
        let Some(list) = self.index.get(&path) else {
            return out;
        };
        let lower = list.partition_point(|&(start, _, _)| {
            (start as u64) < window.0.saturating_sub(SPAN_SCAN_SLACK)
        });
        for &(start, record, mirrored) in &list[lower..] {
            if start as u64 >= window.1 {
                break;
            }
            let (lo, hi) = self.extents[record as usize][usize::from(mirrored)];
            let interval = (start as u64 + lo as u64, start as u64 + hi as u64);
            if interval.1 <= window.0 {
                continue;
            }
            let on_path: Vec<(u32, u32, bool)> = self.placements[record as usize]
                .iter()
                .copied()
                .filter(|&(p, _, _)| p == path)
                .collect();
            let mut instances: Vec<serde_json::Value> = Vec::new();
            for &(read, position) in &self.record_reads[record as usize] {
                let chain = &self.chains[read as usize];
                let entry = &chain[position as usize];
                let next = chain
                    .get(position as usize + 1)
                    .map(|n| {
                        serde_json::json!([n.record, n.mirrored, n.read_start, n.read_end])
                    })
                    .unwrap_or(serde_json::Value::Null);
                instances.push(serde_json::json!({
                    "read": read,
                    "position": position,
                    "mirrored": entry.mirrored,
                    "forward_lane_match": entry.mirrored == mirrored,
                    "read_interval": [entry.read_start, entry.read_end],
                    "next": next,
                }));
            }
            out.push(serde_json::json!({
                "record": record,
                "start": start,
                "mirrored": mirrored,
                "interval": [interval.0, interval.1],
                "extent": [lo, hi],
                "placements_on_path": on_path,
                "instances": instances,
            }));
        }
        out
    }

    /// The restricted junction charge: the seam profile's predicted features
    /// charged against observed evidence contributed ONLY by the crossing
    /// reads of the junction's seam composition (the composition's pairs
    /// unioned, grouped by left row; a read pays at most once per
    /// composition, with its distinct crossing records' equal-share routed
    /// mass at `partitions`).
    pub fn restricted_charge(
        &self,
        profile: &Profile,
        pairs: &[(SourceRange, SourceRange)],
        sources: &routes::Sources,
        path_of_source: &[usize],
        partitions: &[u32],
        model: &ScoreModel,
    ) -> io::Result<RestrictedCharge> {
        let mut charge = 0.0f64;
        let mut events = 0u64;
        if profile.is_empty() {
            return Ok(RestrictedCharge {
                charge,
                spanning_reads: 0,
                events,
                span_feature_count: 0,
            });
        }
        // Group the composition's pairs by left segment: each left row's
        // window is walked once.
        let mut by_left: BTreeMap<(usize, u64, u64, bool), Vec<&SourceRange>> = BTreeMap::new();
        for (left, right) in pairs {
            by_left
                .entry((left.source, left.start, left.end, left.reverse))
                .or_default()
                .push(right);
        }
        let mut by_read: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
        for ((source, start, end, reverse), rights) in by_left {
            let left_seg = SourceRange {
                partition: 0,
                occurrence: 0,
                source,
                start,
                end,
                reverse,
            };
            let crossings = self.crossings_from(&left_seg, path_of_source);
            for right in rights {
                let matched =
                    self.filter_crossings(&crossings, &left_seg, right, sources, path_of_source)?;
                events += matched.events;
                for (read, left_record, right_record) in matched.matches {
                    by_read.entry(read).or_default().insert(left_record);
                    by_read.entry(read).or_default().insert(right_record);
                }
            }
        }
        let mut span: HashMap<&FeatureKey, f64> = HashMap::new();
        for record_set in by_read.values() {
            for &record in record_set {
                let data = &self.records[record as usize];
                if data.touched.is_empty() {
                    continue;
                }
                let share = data
                    .touched
                    .iter()
                    .filter(|partition| partitions.contains(partition))
                    .count() as f64
                    / data.touched.len() as f64;
                if share == 0.0 {
                    continue;
                }
                for feature in &data.features {
                    *span.entry(feature).or_default() += share;
                }
            }
        }
        let mut span_feature_count = 0usize;
        for (feature, &q) in profile.iter() {
            let observed = span.get(feature).copied().unwrap_or(0.0);
            if observed > 0.0 {
                span_feature_count += 1;
            }
            charge += loss_fractional(model, q, observed)?;
        }
        ensure(charge.is_finite(), "nonfinite restricted junction charge")?;
        Ok(RestrictedCharge {
            charge,
            spanning_reads: by_read.len() as u64,
            events,
            span_feature_count,
        })
    }
}
