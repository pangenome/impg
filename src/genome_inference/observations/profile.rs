//! Exact error-free source replay. Enter/leave events keep BOTH raw anchor sets
//! fixed. Between events every relative position translates by the same offset;
//! node/gap MEM queries, selection ties and strict coordinate subwalk pruning are
//! invariant. We query once and multiply its integer, coordinate-attributed count.
use super::super::sample::{collect_raw_views, normalize_reverse, TaggedWalk};
use super::*;

pub fn raw_views(panel: &SyngIndex, dna: &[u8]) -> io::Result<[TaggedWalk; 2]> {
    let f = panel.raw_syncmers_in_sequence(dna)?;
    let rc = crate::graph::reverse_complement(dna);
    let r = normalize_reverse(
        &panel.raw_syncmers_in_sequence(&rc)?,
        dna.len() as u64,
        panel.syncmer_length_bp() as u64,
    );
    Ok([f, r])
}
fn window(w: &[(i32, u64)], start: u64, length: u64, k: u64) -> &[(i32, u64)] {
    let lo = w.partition_point(|&(_, p)| p < start);
    let hi = w.partition_point(|&(_, p)| p + k <= start + length);
    if hi <= lo {
        &[]
    } else {
        &w[lo..hi]
    }
}
pub fn restrict(w: &[(i32, u64)], start: u64, length: u64, k: u64) -> TaggedWalk {
    window(w, start, length, k)
        .iter()
        .map(|&(n, p)| (n, p - start))
        .collect()
}

/// One extraction for a physical core plus read-only source-context halos.
/// Coordinates in both raw streams are relative to crop_start, never interleaved.
pub struct SourceContext {
    pub source: usize,
    pub crop_start: u64,
    pub crop_end: u64,
    pub views: [TaggedWalk; 2],
}

#[derive(Default, Debug, Serialize)]
pub struct CompilationStats {
    pub reused_context_profiles: u64,
    pub fallback_context_fetches: u64,
    pub fallback_context_bp_fetched: u64,
    pub no_start_profiles: u64,
    pub max_profile_context_bp: u64,
    pub max_profile_context_anchors: usize,
}
pub fn add(a: &mut u64, b: u64) -> io::Result<()> {
    *a = a
        .checked_add(b)
        .ok_or_else(|| invalid("integer observation contribution overflow"))?;
    Ok(())
}
pub fn starts(start: u64, end: u64, source_length: u64, length: u64) -> Option<(u64, u64)> {
    if length > source_length || end - start > length {
        return None;
    }
    let lo = end.saturating_sub(length);
    let hi = start.min(source_length - length);
    (lo <= hi).then_some((lo, hi))
}

pub fn total(
    panel: &SyngIndex,
    views: &[TaggedWalk; 2],
    crop_start: u64,
    start: u64,
    end: u64,
    source_length: u64,
    length: u64,
    tokens: &[u64; 3],
) -> io::Result<(u64, u64)> {
    let Some((lo, hi)) = starts(start, end, source_length, length) else {
        return Ok((0, 0));
    };
    let k = panel.syncmer_length_bp() as u64;
    require(
        crop_start <= lo,
        "profile context starts after required read envelope",
    )?;
    let mut events = vec![lo, hi + 1];
    // Even a caller supplying long core traces only visits anchors inside this
    // occurrence's spanning-read envelope, not every anchor in the core.
    for view in views {
        for &(_, p) in window(view, lo - crop_start, hi + length - lo, k) {
            let p = p + crop_start;
            for event in [(p + k).saturating_sub(length), p + 1] {
                if event > lo && event <= hi {
                    events.push(event);
                }
            }
        }
    }
    events.sort_unstable();
    events.dedup();
    let mut sum = 0;
    for range in events.windows(2) {
        let s = range[0];
        let f = restrict(&views[0], s - crop_start, length, k);
        let r = restrict(&views[1], s - crop_start, length, k);
        let records = collect_raw_views(panel, &f, &r, length)?;
        let mut q = 0;
        for record in records {
            for pair in record.windows(2) {
                if pair[0].1 + s == start
                    && pair[1].1 + s + k == end
                    && canonical(&encode_walk(pair)?) == tokens
                {
                    add(&mut q, 1)?;
                }
            }
        }
        add(
            &mut sum,
            q.checked_mul(range[1] - s)
                .ok_or_else(|| invalid("event contribution overflow"))?,
        )?;
    }
    Ok((sum, (events.len() - 1) as u64))
}

/// Bounded fallback/reference path, used for endpoint hits without a live core.
#[cfg(test)]
pub fn compile(
    panel: &SyngIndex,
    sources: &input::Sources,
    catalog: &catalog::Catalog,
    cores: &[Vec<input::Core>],
    lengths: &[u64],
    definition: &input::Definition,
    hit: &mut Incidence,
) -> io::Result<u64> {
    compile_with_context(
        panel,
        sources,
        catalog,
        cores,
        lengths,
        definition,
        hit,
        None,
        &mut CompilationStats::default(),
    )
}

pub fn compile_with_context(
    panel: &SyngIndex,
    sources: &input::Sources,
    catalog: &catalog::Catalog,
    cores: &[Vec<input::Core>],
    lengths: &[u64],
    definition: &input::Definition,
    hit: &mut Incidence,
    shared: Option<&SourceContext>,
    stats: &mut CompilationStats,
) -> io::Result<u64> {
    let source_length = catalog.sources[hit.source].length;
    let mut envelope: Option<(u64, u64)> = None;
    for &length in lengths {
        if let Some((lo, hi)) = starts(hit.start, hit.end, source_length, length) {
            envelope =
                Some(envelope.map_or((lo, hi + length), |(a, b)| (a.min(lo), b.max(hi + length))));
        }
    }
    let source_cores = &cores[hit.source];
    let owner = &source_cores[source_cores.partition_point(|c| c.end <= hit.start)];
    hit.group = (hit.end <= owner.end).then_some(owner.group);
    // Local validity must also cover contexts that would become available at
    // an unvalidated mosaic junction, including beyond original source ends.
    // Keep conditional profiles clipped below, but never saturate away a
    // negative required support bound or treat L > source length as local zero.
    for &length in lengths {
        if length >= hit.end - hit.start {
            let left = i128::from(hit.end) - i128::from(length);
            let right = hit
                .start
                .checked_add(length)
                .ok_or_else(|| invalid("context envelope overflow"))?;
            hit.context_nonlocal |= left < i128::from(owner.start) || right > owner.end;
            hit.source_terminal_context |= left < 0 || right > source_length;
        }
    }
    let Some((lo, end)) = envelope else {
        hit.contributions = vec![[0, 0]; lengths.len()];
        add(&mut stats.no_start_profiles, 1)?;
        return Ok(0);
    };
    let views = if let Some(context) = shared {
        require(
            context.source == hit.source && context.crop_start <= lo && context.crop_end >= end,
            "shared core context does not cover required source/read envelope",
        )?;
        add(&mut stats.reused_context_profiles, 1)?;
        // Binary-slice BOTH traces before event construction. Rebase only this
        // bounded envelope; cloning/scanning the entire core per hit is forbidden.
        let k = panel.syncmer_length_bp() as u64;
        [
            restrict(&context.views[0], lo - context.crop_start, end - lo, k),
            restrict(&context.views[1], lo - context.crop_start, end - lo, k),
        ]
    } else {
        add(&mut stats.fallback_context_fetches, 1)?;
        add(&mut stats.fallback_context_bp_fetched, end - lo)?;
        let dna = sources.fetch(catalog, hit.source, lo, end)?;
        raw_views(panel, &dna)?
    };
    stats.max_profile_context_bp = stats.max_profile_context_bp.max(end - lo);
    stats.max_profile_context_anchors = stats
        .max_profile_context_anchors
        .max(views[0].len() + views[1].len());
    let mut events = 0;
    for &length in lengths {
        let (q, n) = total(
            panel,
            &views,
            lo,
            hit.start,
            hit.end,
            source_length,
            length,
            &definition.tokens,
        )?;
        // Collection explicitly queries both input orientations and canonicalizes
        // their coordinate union: reversing the input swaps the two passes.
        hit.contributions.push([q, q]);
        add(&mut events, n)?;
    }
    Ok(events)
}
