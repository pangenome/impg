//! The canonical anchor-scheme tests (the collapsed index's scheme of
//! record): frame-freeness, single-frame density, and placement
//! preservation of both orientation views against a canonical step index.
//!
//! The scheme: a k-mer qualifies iff its canonical form min(K, rc(K))
//! qualifies — per position, the syng closed-syncmer selection of the frame
//! that spells the canonical form forward. Consequences under test:
//! 1. frame-freeness: the canonical anchor walk of rc(S) is exactly the
//!    mirrored, negated walk of S (the rule is content-determined);
//! 2. the canonical anchor set is a subset of the two frames' raw union
//!    with single-frame-scale density (the twin index's collapse);
//! 3. both orientation views of a record verify against a canonical step
//!    index (a step set built by the canonical rule), where the syng-era
//!    best-orientation record anchors could not (positions differ between
//!    frames by the measured one-base shift).

use impg::genome_inference::mem_records::{
    anchor_scheme, canonical_anchor_walk, raw_matched_syncmers, window_is_canonical_forward,
    AnchorScheme,
};
use impg::graph::reverse_complement;
use impg::syng::{SyncmerParams, SyngIndex};

fn dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut state = seed;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        out.push(bases[((state >> 33) & 3) as usize]);
    }
    out
}

fn test_panel(seeds: &[u64]) -> SyngIndex {
    let sequences: Vec<(String, Vec<u8>)> = seeds
        .iter()
        .enumerate()
        .map(|(i, &seed)| (format!("S#{i}"), dna(3000, seed)))
        .collect();
    SyngIndex::build(SyncmerParams::default(), sequences.into_iter())
}

#[test]
fn window_canonicality_matches_lexicographic_convention() {
    // rc("AC") = "GT": forward spelling is canonical (A<C<G<T order).
    assert!(window_is_canonical_forward(b"AC"));
    // rc("TG") = "CA": not canonical.
    assert!(!window_is_canonical_forward(b"TG"));
    // Non-ACGT windows never qualify.
    assert!(!window_is_canonical_forward(b"ANN"));
    // Even-length palindromic windows tie (never used at the odd syncmer
    // length 63, where a window equal to its own rc is impossible); our
    // tie resolves to forward-canonical by construction.
    assert!(window_is_canonical_forward(b"ACGT"));
    // A long asymmetric window: first differing inner pair decides
    // (A vs comp(A)=T at the outer pair -> forward spelling is canonical).
    let mut w = dna(63, 7);
    w[0] = b'A';
    w[62] = b'A';
    assert!(window_is_canonical_forward(&w));
    w[0] = b'T';
    w[62] = b'T'; // T vs comp(T)=A -> forward spelling is NOT canonical
    assert!(!window_is_canonical_forward(&w));
}

#[test]
fn canonical_anchor_walk_is_frame_free() {
    let panel = test_panel(&[11, 17, 29]);
    let k = panel.syncmer_length_bp() as u64;
    for seed in [101u64, 202, 303] {
        let seq = dna(2000, seed);
        let rc = reverse_complement(&seq);
        let length = seq.len() as u64;
        let forward = canonical_anchor_walk(&panel, &seq).unwrap();
        let reverse = canonical_anchor_walk(&panel, &rc).unwrap();
        // Mirror: rc position q maps to forward length - k - q, node negates.
        let mirrored: Vec<(i32, u64)> = reverse
            .iter()
            .map(|&(node, q)| (-node, length - k - q))
            .collect();
        assert_eq!(
            forward, mirrored,
            "the canonical anchor walk must be frame-free (mirrored rc walk identical)"
        );
    }
}

#[test]
fn canonical_anchors_are_single_frame_subset_of_twin_union() {
    let panel = test_panel(&[11, 17, 29]);
    let k = panel.syncmer_length_bp() as usize;
    // Queries are panel-path substrings, so their 63-mers are dictionary
    // members and both frames' selections surface matched positions.
    for window in [100usize, 400] {
        let seq = dna(3000, 11)[window..window + 2500].to_vec();
        let rc = reverse_complement(&seq);
        let length = seq.len() as u64;
        let fwd: std::collections::BTreeSet<(i32, u64)> =
            raw_matched_syncmers(&panel, &seq).unwrap().into_iter().collect();
        let rev: std::collections::BTreeSet<(i32, u64)> =
            raw_matched_syncmers(&panel, &rc)
                .unwrap()
                .into_iter()
                .map(|(node, q)| (-node, length - k as u64 - q))
                .collect();
        let mut union: std::collections::BTreeSet<(u64, i32)> =
            fwd.iter().map(|&(n, p)| (p, n)).collect();
        union.extend(rev.iter().map(|&(n, p)| (p, n)));
        let canonical = canonical_anchor_walk(&panel, &seq).unwrap();
        for (node, position) in &canonical {
            assert!(
                union.contains(&(*position, *node)),
                "canonical anchor ({node},{position}) must come from a frame's own selection"
            );
        }
        // Single-frame density: the canonical set is no larger than either
        // raw frame's matched set (up to small edge variation) and strictly
        // smaller than the union on random sequence.
        assert!(
            canonical.len() <= fwd.len().max(rev.len()) + 8,
            "canonical anchors must have single-frame scale: {} vs {}/{}",
            canonical.len(),
            fwd.len(),
            rev.len()
        );
        assert!(
            canonical.len() * 2 <= union.len() + 16,
            "the canonical selection must collapse the twin union: {} vs union {}",
            canonical.len(),
            union.len()
        );
    }
}

#[test]
fn canonical_records_verify_against_canonical_steps_both_orientations() {
    // The placement-preservation core: a canonical step index (built by the
    // canonical rule over panel paths) verifies EVERY canonical anchor of a
    // read occurring there — in the read's own orientation on a forward
    // copy, and in the rc orientation on an inverted copy of the same
    // sequence (the twin index's raison d'être, now served by ONE
    // single-frame selection).
    let seeds = [11u64, 23];
    let forward_source = dna(3000, seeds[0]);
    let inverted_source = reverse_complement(&forward_source);
    let sequences: Vec<(String, Vec<u8>)> = vec![
        ("S#0#fwd".to_string(), forward_source.clone()),
        ("S#0#inv".to_string(), inverted_source),
    ];
    let panel = SyngIndex::build(SyncmerParams::default(), sequences.into_iter());
    let k = panel.syncmer_length_bp() as u64;

    // Canonical steps per path (the build_territory_index row logic,
    // miniature): forward GBWT steps + raw rc extraction, per-position
    // canonical-frame choice.
    let canonical_steps = |path_idx: usize, source: &[u8]| -> Vec<(u64, i32)> {
        let path_len = source.len() as u64;
        let forward_steps: Vec<(u64, i32)> = panel
            .walk_path_range(path_idx, 0, path_len)
            .unwrap()
            .into_iter()
            .map(|(node, bp)| (bp, node))
            .collect();
        let rc_source = reverse_complement(source);
        let reverse_steps: Vec<(u64, i32)> = raw_matched_syncmers(&panel, &rc_source)
            .unwrap()
            .into_iter()
            .map(|(node, q)| (source.len() as u64 - k - q, -node))
            .collect();
        let fwd_map: std::collections::BTreeMap<u64, i32> =
            forward_steps.into_iter().collect();
        let rev_map: std::collections::BTreeMap<u64, i32> =
            reverse_steps.into_iter().collect();
        let mut positions: std::collections::BTreeSet<u64> =
            fwd_map.keys().copied().collect();
        positions.extend(rev_map.keys().copied());
        let mut steps = Vec::new();
        for bp in positions {
            let start = bp as usize;
            let canonical_forward =
                window_is_canonical_forward(&source[start..start + k as usize]);
            let chosen = if canonical_forward {
                fwd_map.get(&bp)
            } else {
                rev_map.get(&bp)
            };
            if let Some(&node) = chosen {
                steps.push((bp, node));
            }
        }
        steps
    };
    let forward_steps = canonical_steps(0, &forward_source);
    let inverted_steps = canonical_steps(1, &reverse_complement(&forward_source));
    let fwd_set: std::collections::BTreeSet<(u64, i32)> =
        forward_steps.iter().copied().collect();
    let inv_set: std::collections::BTreeSet<(u64, i32)> =
        inverted_steps.iter().copied().collect();

    let reverse_complement_walk = |anchors: &[(i32, u64)]| -> Vec<(i32, u64)> {
        let span = anchors.last().map(|&(_, p)| p + k).unwrap_or(k);
        let mut out = Vec::with_capacity(anchors.len());
        for &(node, pos) in anchors.iter().rev() {
            out.push((-node, span - k - pos));
        }
        out
    };
    let verify = |anchors: &[(i32, u64)], set: &std::collections::BTreeSet<(u64, i32)>| -> bool {
        anchors
            .iter()
            .all(|&(node, rel)| set.contains(&(rel, node)))
    };

    for (lo, hi) in [(137usize, 789usize), (555, 1301), (911, 1511)] {
        let read = forward_source[lo..hi].to_vec();
        let anchors = canonical_anchor_walk(&panel, &read).unwrap();
        assert!(!anchors.is_empty(), "read must have canonical anchors");
        // FORWARD occurrence on the forward copy: the read's anchor at
        // read position rel sits at locus position lo + rel.
        let shifted: Vec<(i32, u64)> = anchors
            .iter()
            .map(|&(n, r)| (n, lo as u64 + r))
            .collect();
        assert!(
            verify(&shifted, &fwd_set),
            "forward-view anchors must verify against the canonical steps"
        );
        // RC occurrence on the inverted copy: the inverted path spells
        // rc(read) at [len - hi, len - lo). The read's anchor [rel, rel+k)
        // corresponds to rc-read position L - k - rel, i.e. locus position
        // (len - hi) + (L - k - rel) = len - hi - rel + (L - k) ... using the
        // mirrored geometry: locus bp = (len - hi) + ((hi - lo) - k - rel),
        // with the negated node (the dictionary's sign for the rc spelling).
        let inv_len = forward_source.len() as u64;
        let read_len = (hi - lo) as u64;
        let rc_shifted: Vec<(i32, u64)> = anchors
            .iter()
            .map(|&(n, r)| (-n, inv_len - hi as u64 + read_len - k - r))
            .collect();
        assert!(
            verify(&rc_shifted, &inv_set),
            "rc-view anchors must verify against the inverted copy's canonical steps"
        );
    }
}

#[test]
fn anchor_scheme_dispatch_reads_env_once() {
    // The default is the collapsed production scheme; the env override is
    // parsed once and stable in-process.
    let scheme = anchor_scheme();
    assert!(
        scheme == AnchorScheme::Canonical || scheme == AnchorScheme::Syng,
        "scheme must be one of the two"
    );
    assert_eq!(anchor_scheme(), scheme, "scheme is process-stable");
    assert!(AnchorScheme::parse("syng").is_ok());
    assert!(AnchorScheme::parse("canonical").is_ok());
    assert!(AnchorScheme::parse("bogus").is_err());
}
