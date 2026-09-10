//! Bounded, opt-in partition diagnostics; never a scheduling or discovery filter.
use crate::impg::SortedRanges;
use log::{info, warn};
use rustc_hash::FxHashMap;
use std::cell::RefCell;
use std::time::{Duration, Instant};

#[derive(Default)]
pub(crate) struct BackendMetrics {
    pub raw_hits: usize,
    pub raw_anchors: Option<usize>,
    pub chains: Option<usize>,
    pub passing_chains: Option<usize>,
    pub lookup: Duration,
    pub chaining: Option<Duration>,
}

thread_local! {
    // The wrapper runs synchronously on the partition caller, including after
    // its parallel chaining joins. No per-read or per-anchor state is retained.
    static BACKEND: RefCell<Option<Option<BackendMetrics>>> = const { RefCell::new(None) };
}

pub(crate) fn backend_enabled() -> bool {
    BACKEND.with(|b| b.borrow().is_some())
}

pub(crate) fn record_backend(metrics: BackendMetrics) {
    BACKEND.with(|b| {
        if b.borrow().is_some() {
            *b.borrow_mut() = Some(Some(metrics));
        }
    });
}

/// SortedRanges and partition BED coordinates are half-open, despite the
/// coitrees Interval container used elsewhere. Touching endpoints add no bp.
pub(crate) fn missing_core_bp(ranges: Option<&SortedRanges>, start: i32, end: i32) -> u64 {
    let Some(ranges) = ranges else { return 0 };
    let first = ranges.ranges.partition_point(|&(_, e)| e <= start);
    ranges.ranges[first..]
        .iter()
        .take_while(|&&(s, _)| s < end)
        .map(|&(s, e)| (e.min(end) - s.max(start)).max(0) as u64)
        .sum()
}

pub(crate) struct PartitionProfile {
    skip: u64,
    limit: u64,
    queued: u64,
    dispatched: u64,
    sampled: u64,
    visits: FxHashMap<(u32, i32, i32), u64>,
}

impl PartitionProfile {
    pub fn from_env() -> Option<Self> {
        let limit = std::env::var_os("IMPG_PARTITION_PROFILE_MAX_QUERIES")?;
        let skip =
            std::env::var("IMPG_PARTITION_PROFILE_SKIP_QUERIES").unwrap_or_else(|_| "0".into());
        match (limit.to_string_lossy().parse::<u64>(), skip.parse::<u64>()) {
            (Ok(limit @ 1..=10_000), Ok(skip)) => Some(Self::new(skip, limit)),
            _ => {
                warn!("Partition profiling disabled: MAX_QUERIES must be 1..=10000 and SKIP_QUERIES a nonnegative integer");
                None
            }
        }
    }

    fn new(skip: u64, limit: u64) -> Self {
        Self {
            skip,
            limit,
            queued: 0,
            dispatched: 0,
            sampled: 0,
            visits: FxHashMap::default(),
        }
    }

    pub fn queued(&mut self, count: usize) {
        self.queued += count as u64;
    }

    pub fn dispatch(
        &mut self,
        key: (u32, i32, i32),
        missing: Option<&SortedRanges>,
    ) -> Option<QueryProfile> {
        self.dispatched += 1;
        if self.dispatched <= self.skip || self.sampled >= self.limit {
            return None;
        }
        self.sampled += 1;
        let visit = self.visits.entry(key).or_default();
        *visit += 1;
        let previous = BACKEND.with(|b| b.replace(Some(None)));
        Some(QueryProfile {
            ordinal: self.dispatched,
            queued: self.queued,
            key,
            missing: missing_core_bp(missing, key.1, key.2),
            visit: *visit,
            start: Instant::now(),
            previous,
        })
    }
}

impl Drop for PartitionProfile {
    fn drop(&mut self) {
        info!("partition_profile_summary queued={} dispatched={} sampled={} unique_query_tuples_in_sample={} skip={} limit={}",
            self.queued, self.dispatched, self.sampled, self.visits.len(), self.skip, self.limit);
    }
}

pub(crate) struct QueryProfile {
    ordinal: u64,
    queued: u64,
    key: (u32, i32, i32),
    missing: u64,
    visit: u64,
    start: Instant,
    previous: Option<Option<BackendMetrics>>,
}

impl QueryProfile {
    pub fn started(&self, name: &str) {
        info!("partition_profile_dispatch dispatched={} queued={} seq_id={} name={} start={} end={} query_bp={} missing_core_bp={} visit_in_sample={}",
            self.ordinal, self.queued, self.key.0, name, self.key.1, self.key.2,
            self.key.2 - self.key.1, self.missing, self.visit);
    }

    pub fn elapsed(&self) -> Duration {
        self.start.elapsed()
    }

    pub fn finish(
        &self,
        query: Duration,
        postprocess: Duration,
        mask: Duration,
        backend_overlaps: usize,
        emitted_bp: u64,
        error: bool,
    ) {
        info!("partition_profile_query dispatched={} queued={} seq_id={} start={} end={} query_bp={} missing_core_bp={} visit_in_sample={} query_us={} postprocess_us={} mask_us={} output_us={} backend_overlaps={} emitted_new_source_bp={} error={}",
            self.ordinal, self.queued, self.key.0, self.key.1, self.key.2, self.key.2 - self.key.1,
            self.missing, self.visit, query.as_micros(), postprocess.as_micros(), mask.as_micros(),
            self.elapsed().saturating_sub(query + postprocess).as_micros(), backend_overlaps, emitted_bp, error);
        BACKEND.with(|b| {
            // Non-syng backends and failed lookups do not publish metrics.
            if let Some(Some(m)) = b.borrow().as_ref() {
                info!("partition_profile_syng dispatched={} raw_hits={} raw_anchors={:?} chains={:?} passing_chains={:?} lookup_us={} chaining_us={:?}",
                    self.ordinal, m.raw_hits, m.raw_anchors, m.chains, m.passing_chains,
                    m.lookup.as_micros(), m.chaining.map(|d| d.as_micros()));
            }
        });
    }
}

impl Drop for QueryProfile {
    fn drop(&mut self) {
        BACKEND.with(|b| {
            b.replace(self.previous.take());
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_core_uses_half_open_intersections() {
        let mut r = SortedRanges::new(100, 0);
        r.insert((10, 20));
        r.insert((30, 31));
        r.insert((40, 50));
        for (s, e, bp) in [
            (0, 10, 0),
            (20, 30, 0),
            (20, 31, 1),
            (19, 41, 3),
            (0, 100, 21),
            (50, 100, 0),
            (10, 10, 0),
        ] {
            assert_eq!(missing_core_bp(Some(&r), s, e), bp);
        }
        assert_eq!(missing_core_bp(None, 0, 100), 0);
    }

    #[test]
    fn diagnostics_bound_storage_and_distinguish_exact_tuple_revisits() {
        let mut p = PartitionProfile::new(1, 3);
        p.queued(10);
        assert!(p.dispatch((0, 0, 10), None).is_none());
        assert!(!backend_enabled());
        for (key, visit) in [((0, 0, 10), 1), ((0, 10, 20), 1), ((0, 0, 10), 2)] {
            let q = p.dispatch(key, None).unwrap();
            assert_eq!(q.visit, visit);
            assert!(backend_enabled());
            drop(q);
            assert!(!backend_enabled());
        }
        for _ in 0..100 {
            assert!(p.dispatch((1, 0, 10), None).is_none());
        }
        assert_eq!(p.visits.len(), 2);
        assert_eq!(p.sampled, 3);
        assert_eq!(p.dispatched, 104);
    }
}
