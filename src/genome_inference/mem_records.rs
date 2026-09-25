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

/// Return the complete canonical maximal-MEM records retained for one read.
/// Coordinates and read identity are intentionally absent, matching
/// `WeightedBwt` construction exactly. Experimental variable-length scorers
/// should derive node-to-node subwalk keys from these records rather than
/// `observed_pairs`.
pub fn canonical_mem_records(panel: &SyngIndex, sequence: &[u8]) -> io::Result<Vec<Vec<u64>>> {
    super::sample::collect_tagged_read(panel, sequence)?
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
    super::sample::collect_tagged_read(panel, sequence)
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
    let raw = panel.raw_syncmers_in_sequence(sequence)?;
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
