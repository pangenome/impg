//! Complete canonical maximal-MEM records for one read, for experimental
//! variable-length subwalk scorers.
//!
//! This accessor lives in its own module on purpose: `panel_routes` and
//! `observations` `compiler_identity` values are source-bound over
//! `sample.rs`, so adding public API there would invalidate already-built
//! route artifacts. This module is outside that hashed source set.
use crate::sample_mem_bwt::{canonical, encode_walk};
use crate::syng::SyngIndex;
use std::io;

/// The record-anchor qualification scheme (the collapsed machinery's
/// scheme of record). `Syng` is the historical extraction: per-read
/// best-orientation selection of syng's own (rc-ASYNCHRONOUS, measured
/// one-base-shifted) closed-syncmer frames. `Canonical` is the collapse:
/// a k-mer qualifies iff its canonical form min(K, rc(K)) qualifies — per
/// position the syng selection of the frame that spells the canonical
/// form forward — which is frame-free (the same physical k-mer gets the
/// same decision in both frames), so one index serves both strands and
/// the read side needs no orientation selection at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnchorScheme {
    Syng,
    Canonical,
}

impl AnchorScheme {
    pub fn parse(value: &str) -> io::Result<Self> {
        match value {
            "syng" => Ok(AnchorScheme::Syng),
            "canonical" => Ok(AnchorScheme::Canonical),
            _ => Err(io::Error::other(
                "IMPG_ANCHOR_SCHEME must be 'syng' or 'canonical'",
            )),
        }
    }
}

/// The process-wide anchor scheme. Default: `canonical` (the collapsed
/// production path). Set IMPG_ANCHOR_SCHEME=syng to run the historical
/// extraction (the transition oracle; the strandfix-era baseline).
pub fn anchor_scheme() -> AnchorScheme {
    static SCHEME: std::sync::OnceLock<AnchorScheme> = std::sync::OnceLock::new();
    *SCHEME.get_or_init(|| {
        std::env::var("IMPG_ANCHOR_SCHEME")
            .ok()
            .and_then(|value| AnchorScheme::parse(&value).ok())
            .unwrap_or(AnchorScheme::Canonical)
    })
}

/// Is the window's forward spelling its own canonical form (the same
/// A<C<G<T order syng's own `isCanonical` uses; equality is impossible at
/// odd window length, and non-ACGT windows never qualify)?
pub fn window_is_canonical_forward(window: &[u8]) -> bool {
    let len = window.len();
    let mut lo = 0usize;
    let mut hi = len - 1;
    loop {
        let (a, b) = (window[lo], window[hi]);
        let comp_b = match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => return false,
        };
        match a {
            b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => {}
            _ => return false,
        }
        if a != comp_b {
            return a < comp_b;
        }
        if lo + 1 >= hi {
            return true;
        }
        lo += 1;
        hi -= 1;
    }
}

/// The canonical scheme's matched anchor walk of one sequence in its OWN
/// frame: every position whose 63-mer qualifies under the canonical rule
/// (its canonical form qualifies — the syng selection of the frame that
/// spells the canonical form forward) AND is present in the panel's syncmer
/// dictionary. Node signs are the dictionary's own (spelled in THIS frame);
/// positions increase. Frame-free by construction: the reverse complement
/// view yields exactly the mirrored, negated walk.
pub fn canonical_anchor_walk(panel: &SyngIndex, sequence: &[u8]) -> io::Result<Vec<(i32, u64)>> {
    let k = panel.syncmer_length_bp() as usize;
    if sequence.len() < k {
        return Ok(Vec::new());
    }
    let length = sequence.len() as u64;
    let forward: std::collections::BTreeMap<u64, i32> = raw_matched_syncmers(panel, sequence)?
        .into_iter()
        .map(|(node, pos)| (pos, node))
        .collect();
    let reverse_sequence = crate::graph::reverse_complement(sequence);
    let reverse: std::collections::BTreeMap<u64, i32> =
        raw_matched_syncmers(panel, &reverse_sequence)?
            .into_iter()
            .map(|(node, q)| (length - k as u64 - q, -node))
            .collect();
    let mut positions: std::collections::BTreeSet<u64> = forward.keys().copied().collect();
    positions.extend(reverse.keys().copied());
    let mut out = Vec::new();
    for position in positions {
        let start = position as usize;
        let Some(window) = sequence.get(start..start + k) else {
            continue;
        };
        let chosen = if window_is_canonical_forward(window) {
            forward.get(&position)
        } else {
            reverse.get(&position)
        };
        if let Some(&node) = chosen {
            out.push((node, position));
        }
    }
    Ok(out)
}

/// The canonical scheme's maximal-MEM records of one sequence in its own
/// frame: GBWT walks over the canonical anchor walk (the same chaining
/// `collect_tagged_read` applies, without any orientation selection — the
/// canonical anchors are frame-free, so both orientation views yield the
/// same records up to mirroring, measured zero additions on real reads).
pub fn canonical_tagged_mem_records(
    panel: &SyngIndex,
    sequence: &[u8],
) -> io::Result<Vec<super::sample::TaggedWalk>> {
    let walk: Vec<crate::syng::SyngWalkStep> = canonical_anchor_walk(panel, sequence)?
        .into_iter()
        .map(|(signed_node, bp_pos)| crate::syng::SyngWalkStep {
            signed_node,
            bp_pos,
        })
        .collect();
    let mut records = Vec::new();
    for mem in panel.gbwt_mems_for_walk(&walk)? {
        records.push(
            walk[mem.step_start..mem.step_end]
                .iter()
                .map(|step| (step.signed_node, step.bp_pos))
                .collect(),
        );
    }
    Ok(records)
}

/// Return the complete canonical maximal-MEM records retained for one read.
/// Coordinates and read identity are intentionally absent, matching
/// `WeightedBwt` construction exactly. Experimental variable-length scorers
/// should derive node-to-node subwalk keys from these records rather than
/// `observed_pairs`.
pub fn canonical_mem_records(panel: &SyngIndex, sequence: &[u8]) -> io::Result<Vec<Vec<u64>>> {
    let tagged = match anchor_scheme() {
        AnchorScheme::Syng => super::sample::collect_tagged_read(panel, sequence)?,
        AnchorScheme::Canonical => canonical_tagged_mem_records(panel, sequence)?,
    };
    tagged
        .iter()
        .map(|r| encode_walk(r).map(|t| canonical(&t)))
        .collect()
}

/// The complete maximal-MEM records of one read WITH read coordinates and
/// strand (the within-read adjacency-chain view): each record is a walk of
/// (signed node, input-forward read position), positions increasing along
/// the read; negative nodes mean the read matches the reverse complement of
/// the panel path region. This is the same extraction `canonical_mem_records`
/// wraps, before canonicalization discards coordinates.
pub fn tagged_mem_records(
    panel: &SyngIndex,
    sequence: &[u8],
) -> io::Result<Vec<Vec<(i32, u64)>>> {
    tagged_mem_records_scheme(panel, sequence, anchor_scheme())
}

/// `tagged_mem_records` with an explicit scheme: every production record
/// extraction (the routing pass's per-run re-derivation, the candidate
/// predicted profiles' read-window simulations, the spell diagnostics)
/// goes through this dispatcher so the anchor qualification is coherent
/// across the whole machinery under either scheme.
pub fn tagged_mem_records_scheme(
    panel: &SyngIndex,
    sequence: &[u8],
    scheme: AnchorScheme,
) -> io::Result<Vec<Vec<(i32, u64)>>> {
    match scheme {
        AnchorScheme::Syng => super::sample::collect_tagged_read(panel, sequence),
        AnchorScheme::Canonical => canonical_tagged_mem_records(panel, sequence),
    }
}

/// The RAW single-frame matched-syncmer extraction of one sequence AS GIVEN:
/// every position whose k-mer qualifies in the sequence's own frame AND
/// resolves in the panel's node space, with its signed node and its own-frame
/// position. Unlike `SyngIndex::matched_syncmers_in_sequence` this performs
/// NO best-orientation selection — the orientation selector discards the
/// losing frame, and the losing frame is exactly where a sequence's
/// rc-frame-only qualifying k-mers live (the strand-asymmetry defect's
/// extraction face). The router's rc-symmetric territory augmentation needs
/// the losing frame's raw selection on a path range's reverse complement.
/// This wrapper lives in this module on purpose: mem_records.rs is outside
/// every source-bound compiler_identity hash set (see the module header), so
/// widening the raw extractor's reach cannot invalidate built artifacts.
pub fn raw_matched_syncmers(panel: &SyngIndex, sequence: &[u8]) -> io::Result<Vec<(i32, u64)>> {
    panel.raw_syncmers_in_sequence(sequence)
}

/// Input-forward anchor positions for both raw sketch orientations. These are
/// sufficient to derive exact enter/leave boundaries for event-compressed replay.
pub fn raw_mem_anchor_positions(panel: &SyngIndex, sequence: &[u8]) -> io::Result<[Vec<u64>; 2]> {
    let reverse = crate::graph::reverse_complement(sequence);
    let k = panel.syncmer_length_bp() as u64;
    let length = sequence.len() as u64;
    let forward = panel
        .raw_syncmers_in_sequence(sequence)?
        .into_iter()
        .map(|(_, position)| position)
        .collect();
    let mut reverse_positions = panel
        .raw_syncmers_in_sequence(&reverse)?
        .into_iter()
        .map(|(_, position)| length - k - position)
        .collect::<Vec<_>>();
    reverse_positions.sort_unstable();
    Ok([forward, reverse_positions])
}

/// Own-orientation maximal MEM records of one query sequence, WITHOUT the
/// public matcher's best-orientation selection: the raw syncmers syng
/// extracts from the sequence exactly as given, kept when their k-mer
/// matches a panel syncmer node, walked through the GBWT. Records are
/// (signed node, position) walks in query coordinates, positions increasing.
///
/// Experimental accessor for the geometric predicted-profile derivation: the
/// reverse-complement view of a candidate must see the panel matches of the
/// query's own syncmer phase (the inverted-duplicate loci), which the public
/// best-orientation matcher hides whenever the forward-matching orientation
/// has more matched syncmers over the whole query.
pub fn own_orientation_mem_records(
    panel: &crate::syng::SyngIndex,
    sequence: &[u8],
) -> io::Result<Vec<Vec<(i32, u64)>>> {
    own_orientation_mem_records_scheme(panel, sequence, anchor_scheme())
}

/// `own_orientation_mem_records` with an explicit scheme: under
/// `Canonical` the anchors are the canonical scheme's own-frame walk (the
/// records the candidate predicted-profile machinery must simulate to stay
/// coherent with the per-run re-derived evidence records).
pub fn own_orientation_mem_records_scheme(
    panel: &crate::syng::SyngIndex,
    sequence: &[u8],
    scheme: AnchorScheme,
) -> io::Result<Vec<Vec<(i32, u64)>>> {
    let raw = match scheme {
        AnchorScheme::Syng => panel.raw_syncmers_in_sequence(sequence)?,
        AnchorScheme::Canonical => canonical_anchor_walk(panel, sequence)?,
    };
    let walk: Vec<crate::syng::SyngWalkStep> = raw
        .iter()
        .map(|&(signed_node, bp_pos)| crate::syng::SyngWalkStep {
            signed_node,
            bp_pos,
        })
        .collect();
    Ok(panel
        .gbwt_mems_for_walk(&walk)?
        .into_iter()
        .map(|mem| {
            walk[mem.step_start..mem.step_end]
                .iter()
                .map(|step| (step.signed_node, step.bp_pos))
                .collect()
        })
        .collect())
}
