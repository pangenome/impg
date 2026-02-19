//! End-to-end integration tests for the realize engine.
//!
//! These tests exercise `realize_from_sequences` with various inputs
//! and configurations, validating GFA output structure and statistics.


use impg::graph::SequenceMetadata;
use impg::realize::{realize_from_sequences, RealizeConfig, RealizeResult};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build a (sequence, metadata) pair for testing.
fn make_seq(name: &str, seq: &str, start: i32, total_length: usize) -> (String, SequenceMetadata) {
    (
        seq.to_string(),
        SequenceMetadata {
            name: name.to_string(),
            start,
            size: seq.len() as i32,
            strand: '+',
            total_length,
        },
    )
}

/// Build a (sequence, metadata) pair on the reverse strand.
fn make_seq_rev(
    name: &str,
    seq: &str,
    start: i32,
    total_length: usize,
) -> (String, SequenceMetadata) {
    (
        seq.to_string(),
        SequenceMetadata {
            name: name.to_string(),
            start,
            size: seq.len() as i32,
            strand: '-',
            total_length,
        },
    )
}

/// Generate a random-ish DNA sequence of the given length (deterministic from seed).
fn make_dna(len: usize, seed: u64) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut out = Vec::with_capacity(len);
    let mut state = seed;
    for _ in 0..len {
        // Simple LCG
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        out.push(bases[((state >> 33) % 4) as usize]);
    }
    String::from_utf8(out).unwrap()
}

/// Mutate a DNA string: introduce a SNP every `interval` bases.
fn mutate_snps(seq: &str, interval: usize) -> String {
    let mut out: Vec<u8> = seq.bytes().collect();
    for i in (0..out.len()).step_by(interval) {
        out[i] = match out[i] {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => out[i],
        };
    }
    String::from_utf8(out).unwrap()
}

/// Default config for small test sequences: no sorting (avoids gfasort dependency),
/// small thresholds.
fn test_config() -> RealizeConfig {
    RealizeConfig {
        poa_threshold: 1000,
        chunk_size: 500,
        padding: 0,
        max_depth: 10,
        num_threads: 1,
        scoring_params: (5, 4, 6, 2, 24, 1),
        temp_dir: None,
        sort_output: false,
    }
}

/// Count lines starting with the given prefix in a GFA string.
fn count_gfa_lines(gfa: &str, prefix: &str) -> usize {
    gfa.lines().filter(|l| l.starts_with(prefix)).count()
}

/// Extract path names from a GFA string.
fn path_names(gfa: &str) -> Vec<String> {
    gfa.lines()
        .filter(|l| l.starts_with("P\t"))
        .filter_map(|l| l.split('\t').nth(1).map(String::from))
        .collect()
}

/// Validate basic GFA structural properties.
fn validate_gfa(result: &RealizeResult) {
    let gfa = &result.gfa;
    assert!(
        gfa.contains("H\t") || gfa.starts_with("H\t"),
        "GFA should have a header line"
    );

    let segments = count_gfa_lines(gfa, "S\t");
    let paths = count_gfa_lines(gfa, "P\t");
    let links = count_gfa_lines(gfa, "L\t");

    assert!(segments > 0, "GFA should have at least one segment");
    assert!(paths > 0, "GFA should have at least one path");

    // Each path should reference segments that exist.
    let segment_ids: std::collections::HashSet<String> = gfa
        .lines()
        .filter(|l| l.starts_with("S\t"))
        .filter_map(|l| l.split('\t').nth(1).map(String::from))
        .collect();

    for line in gfa.lines().filter(|l| l.starts_with("P\t")) {
        let steps_field = line.split('\t').nth(2).unwrap_or("");
        for step in steps_field.split(',') {
            let node_id = step.trim_end_matches('+').trim_end_matches('-');
            assert!(
                segment_ids.contains(node_id),
                "Path step references non-existent segment {}: full GFA:\n{}",
                node_id,
                gfa
            );
        }
    }

    // Links should reference existing segments.
    for line in gfa.lines().filter(|l| l.starts_with("L\t")) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 4 {
            assert!(
                segment_ids.contains(fields[1]),
                "Link references non-existent from-segment {}: full GFA:\n{}",
                fields[1],
                gfa
            );
            assert!(
                segment_ids.contains(fields[3]),
                "Link references non-existent to-segment {}: full GFA:\n{}",
                fields[3],
                gfa
            );
        }
    }

    // If there are multiple segments, there should be at least one link.
    if segments > 1 {
        assert!(
            links > 0,
            "GFA with {} segments should have links, got 0. GFA:\n{}",
            segments,
            gfa
        );
    }
}

/// Concatenate the segment sequences traversed by a path (following orientations).
/// This reconstructs the "spelled" sequence from the GFA.
fn spell_path(gfa: &str, path_name: &str) -> Option<String> {
    // Build segment map: id -> sequence
    let segments: std::collections::HashMap<String, String> = gfa
        .lines()
        .filter(|l| l.starts_with("S\t"))
        .filter_map(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 3 {
                Some((fields[1].to_string(), fields[2].to_string()))
            } else {
                None
            }
        })
        .collect();

    // Find the path line
    let path_line = gfa
        .lines()
        .filter(|l| l.starts_with("P\t"))
        .find(|l| l.split('\t').nth(1) == Some(path_name))?;

    let steps_field = path_line.split('\t').nth(2)?;

    let mut result = String::new();
    for step in steps_field.split(',') {
        let (node_id, is_reverse) = if let Some(s) = step.strip_suffix('-') {
            (s, true)
        } else if let Some(s) = step.strip_suffix('+') {
            (s, false)
        } else {
            (step, false)
        };

        let seg_seq = segments.get(node_id)?;
        if is_reverse {
            // Reverse complement
            let rc: String = seg_seq
                .bytes()
                .rev()
                .map(|b| match b {
                    b'A' | b'a' => 'T',
                    b'T' | b't' => 'A',
                    b'C' | b'c' => 'G',
                    b'G' | b'g' => 'C',
                    _ => b as char,
                })
                .collect();
            result.push_str(&rc);
        } else {
            result.push_str(seg_seq);
        }
    }

    Some(result)
}

// ===========================================================================
// Tests: empty and trivial inputs
// ===========================================================================

#[test]
fn test_realize_empty_input() {
    let config = test_config();
    let result = realize_from_sequences(&[], &config).unwrap();
    assert!(result.gfa.contains("H\tVN:Z:1.0"));
    assert_eq!(result.stats.num_sequences, 0);
    assert_eq!(result.stats.poa_calls, 0);
    assert_eq!(result.stats.sweepga_calls, 0);
}

#[test]
fn test_realize_single_sequence() {
    let config = test_config();
    let seqs = vec![make_seq("chr1", "ACGTACGTACGT", 0, 100)];
    let result = realize_from_sequences(&seqs, &config).unwrap();

    validate_gfa(&result);
    assert_eq!(result.stats.num_sequences, 1);
    assert_eq!(result.stats.poa_calls, 1);
    assert_eq!(result.stats.sweepga_calls, 0);

    let names = path_names(&result.gfa);
    assert_eq!(names.len(), 1);
    assert!(names[0].contains("chr1"), "Path should reference chr1");

    // The spelled sequence should match the input.
    let spelled = spell_path(&result.gfa, &names[0]).unwrap();
    assert_eq!(spelled, "ACGTACGTACGT");
}

// ===========================================================================
// Tests: POA base case (below threshold)
// ===========================================================================

#[test]
fn test_realize_two_identical_sequences() {
    let config = test_config();
    let seqs = vec![
        make_seq("s1", "ACGTACGTACGTACGTACGT", 0, 100),
        make_seq("s2", "ACGTACGTACGTACGTACGT", 10, 100),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.num_sequences, 2);
    assert_eq!(result.stats.poa_calls, 1);
    assert_eq!(result.stats.sweepga_calls, 0);

    let names = path_names(&result.gfa);
    assert_eq!(names.len(), 2);

    // Both paths spell the same sequence.
    let s1 = spell_path(&result.gfa, &names[0]).unwrap();
    let s2 = spell_path(&result.gfa, &names[1]).unwrap();
    assert_eq!(s1, "ACGTACGTACGTACGTACGT");
    assert_eq!(s2, "ACGTACGTACGTACGTACGT");
}

#[test]
fn test_realize_two_sequences_with_snp() {
    let config = test_config();
    // Two sequences differing by a single SNP in the middle.
    let seq1 = "AAAAACCCCCGGGGGAAAAACCCCC";
    let mut seq2_bytes: Vec<u8> = seq1.bytes().collect();
    seq2_bytes[12] = b'T'; // G→T at position 12
    let seq2 = String::from_utf8(seq2_bytes).unwrap();

    let seqs = vec![
        make_seq("s1", seq1, 0, 100),
        make_seq("s2", &seq2, 0, 100),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.num_sequences, 2);

    // A SNP should produce a bubble: at least 2 segments.
    let segments = count_gfa_lines(&result.gfa, "S\t");
    assert!(
        segments >= 2,
        "SNP should create multiple segments, got {}",
        segments
    );

    // Paths should spell back the original sequences.
    let names = path_names(&result.gfa);
    let s1 = spell_path(&result.gfa, &names[0]).unwrap();
    let s2 = spell_path(&result.gfa, &names[1]).unwrap();
    assert_eq!(s1, seq1);
    assert_eq!(s2, seq2);
}

#[test]
fn test_realize_three_sequences_with_indel() {
    let config = test_config();
    // Base:    AAAAACCCCCGGGGG
    // Insert:  AAAAACCCTTCCCGGGGG  (TT inserted after position 8)
    // Delete:  AAAAAGGGGG          (CCCCC deleted)
    let seqs = vec![
        make_seq("base", "AAAAACCCCCGGGGG", 0, 100),
        make_seq("ins", "AAAAACCCTTCCCGGGGG", 0, 100),
        make_seq("del", "AAAAAGGGGG", 0, 100),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.num_sequences, 3);
    let names = path_names(&result.gfa);
    assert_eq!(names.len(), 3);

    // Each path should spell its original sequence.
    for (name, expected_seq) in names.iter().zip(["AAAAACCCCCGGGGG", "AAAAACCCTTCCCGGGGG", "AAAAAGGGGG"].iter()) {
        let spelled = spell_path(&result.gfa, name).unwrap();
        assert_eq!(
            &spelled, expected_seq,
            "Path {} spelled {} but expected {}",
            name, spelled, expected_seq
        );
    }
}

#[test]
fn test_realize_reverse_strand_metadata() {
    let config = test_config();
    let seqs = vec![
        make_seq("s1", "ACGTACGTACGTACGTACGT", 0, 100),
        make_seq_rev("s2", "ACGTACGTACGTACGTACGT", 80, 100),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.num_sequences, 2);
    let names = path_names(&result.gfa);
    assert_eq!(names.len(), 2);
}

// ===========================================================================
// Tests: config variations
// ===========================================================================

#[test]
fn test_realize_with_padding() {
    // With padding > 0, padded POA is used.
    let mut config = test_config();
    config.padding = 5;

    let seqs = vec![
        make_seq("s1", "AAAAACCCCCGGGGG", 5, 100),
        make_seq("s2", "AAAAACCCCCGGGGG", 10, 100),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);
    assert_eq!(result.stats.num_sequences, 2);
}

#[test]
fn test_realize_max_depth_zero_forces_poa() {
    // max_depth=0 means the very first call is at depth 0 >= max_depth, so POA.
    let mut config = test_config();
    config.max_depth = 0;

    // Even a long sequence should not trigger sweepga.
    let long_seq = make_dna(2000, 42);
    let long_seq2 = mutate_snps(&long_seq, 50);

    let seqs = vec![
        make_seq("s1", &long_seq, 0, 3000),
        make_seq("s2", &long_seq2, 0, 3000),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.sweepga_calls, 0, "max_depth=0 should skip sweepga");
    assert!(result.stats.poa_calls >= 1, "Should have at least 1 POA call");
}

#[test]
fn test_realize_high_poa_threshold_skips_recursion() {
    // With poa_threshold very high, everything goes to POA.
    let mut config = test_config();
    config.poa_threshold = 100_000;
    config.padding = 0;

    let seq = make_dna(500, 123);
    let seq2 = mutate_snps(&seq, 25);

    let seqs = vec![
        make_seq("s1", &seq, 0, 1000),
        make_seq("s2", &seq2, 0, 1000),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.sweepga_calls, 0);
    assert_eq!(result.stats.poa_calls, 1);
}

// ===========================================================================
// Tests: large regions triggering sweepga + recursion
// ===========================================================================

#[test]
fn test_realize_large_region_triggers_sweepga() {
    // Two similar sequences longer than poa_threshold → sweepga + partition.
    let mut config = test_config();
    config.poa_threshold = 500;
    config.chunk_size = 300;
    config.padding = 50;
    config.max_depth = 5;


    // Generate a 2kb sequence and a variant with ~2% divergence.
    let base_seq = make_dna(2000, 77);
    let variant_seq = mutate_snps(&base_seq, 50);

    let seqs = vec![
        make_seq("anchor", &base_seq, 0, 3000),
        make_seq("variant", &variant_seq, 0, 3000),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.num_sequences, 2);
    assert!(
        result.stats.sweepga_calls >= 1,
        "Expected sweepga to be called for a 2kb region (threshold=500), got {} sweepga calls",
        result.stats.sweepga_calls
    );
    assert!(
        result.stats.poa_calls >= 1,
        "Recursion should bottom out in POA"
    );
    // sweepga + POA means the recursive engine did its job: aligned, partitioned,
    // and resolved chunks via POA (possibly as fallbacks).
    assert!(
        result.stats.sweepga_calls + result.stats.poa_calls >= 2,
        "Expected at least one sweepga call and one POA call, got sweepga={} poa={}",
        result.stats.sweepga_calls,
        result.stats.poa_calls
    );

    let names = path_names(&result.gfa);
    assert_eq!(names.len(), 2, "Should have paths for both sequences");
}

#[test]
fn test_realize_three_sequences_large_region() {
    // Three sequences above the threshold → sweepga + partition + lace.
    let mut config = test_config();
    config.poa_threshold = 400;
    config.chunk_size = 250;
    config.padding = 30;
    config.max_depth = 5;


    let base = make_dna(1500, 99);
    let var1 = mutate_snps(&base, 40);
    let var2 = mutate_snps(&base, 60);

    let seqs = vec![
        make_seq("s1", &base, 0, 2000),
        make_seq("s2", &var1, 0, 2000),
        make_seq("s3", &var2, 0, 2000),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.num_sequences, 3);
    let names = path_names(&result.gfa);
    assert_eq!(names.len(), 3);

    // Segments should exist and paths should reference them.
    let segments = count_gfa_lines(&result.gfa, "S\t");
    assert!(segments >= 2, "Expected multiple segments, got {}", segments);
}

// ===========================================================================
// Tests: GFA structural validation
// ===========================================================================

#[test]
fn test_realize_gfa_has_header() {
    let config = test_config();
    let seqs = vec![make_seq("s1", "ACGTACGT", 0, 100)];
    let result = realize_from_sequences(&seqs, &config).unwrap();
    assert!(
        result.gfa.lines().any(|l| l.starts_with("H\t")),
        "GFA must have header"
    );
}

#[test]
fn test_realize_gfa_path_names_match_input() {
    let config = test_config();
    let seqs = vec![
        make_seq("chr1", "ACGTACGTACGTACGTACGT", 100, 1000),
        make_seq("chr2", "ACGTACGTACGTACGTACGT", 200, 1000),
    ];
    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    let names = path_names(&result.gfa);
    assert_eq!(names.len(), 2);

    // Path names should contain the chromosome names with coordinate ranges.
    assert!(
        names.iter().any(|n| n.contains("chr1")),
        "Expected a path for chr1, got: {:?}",
        names
    );
    assert!(
        names.iter().any(|n| n.contains("chr2")),
        "Expected a path for chr2, got: {:?}",
        names
    );
}

#[test]
fn test_realize_stats_consistency() {
    let mut config = test_config();
    config.poa_threshold = 300;
    config.chunk_size = 200;
    config.padding = 20;
    config.max_depth = 3;


    let base = make_dna(1000, 55);
    let var = mutate_snps(&base, 30);

    let seqs = vec![
        make_seq("s1", &base, 0, 2000),
        make_seq("s2", &var, 0, 2000),
    ];

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.num_sequences, 2);
    assert!(result.stats.total_ms < 60_000, "Should complete in <60s");
    assert!(
        result.stats.poa_calls + result.stats.sweepga_calls >= 1,
        "Should have at least one call"
    );
    assert!(
        result.stats.max_depth_reached <= config.max_depth,
        "Depth {} should not exceed max {}",
        result.stats.max_depth_reached,
        config.max_depth
    );
}

// ===========================================================================
// Tests: many sequences
// ===========================================================================

#[test]
fn test_realize_five_sequences_poa() {
    let config = test_config();

    let base = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32bp
    let seqs: Vec<(String, SequenceMetadata)> = (0..5)
        .map(|i| {
            let mut s: Vec<u8> = base.bytes().collect();
            // Introduce a unique SNP per sequence.
            s[10 + i] = b'N';
            make_seq(
                &format!("s{}", i),
                &String::from_utf8(s).unwrap(),
                i as i32 * 10,
                200,
            )
        })
        .collect();

    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.num_sequences, 5);
    let names = path_names(&result.gfa);
    assert_eq!(names.len(), 5);
}

// ===========================================================================
// Tests: edge cases
// ===========================================================================

#[test]
fn test_realize_very_short_sequences() {
    // Sequences at the minimum viable length for POA.
    let config = test_config();
    let seqs = vec![
        make_seq("s1", "ACGT", 0, 10),
        make_seq("s2", "ACGT", 0, 10),
    ];
    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);
    assert_eq!(result.stats.num_sequences, 2);
}

#[test]
fn test_realize_sequences_of_different_lengths() {
    let config = test_config();
    // One long and one short sequence.
    let seqs = vec![
        make_seq("long", "ACGTACGTACGTACGTACGTACGTACGT", 0, 100),
        make_seq("short", "ACGTACGT", 0, 50),
    ];
    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);
    assert_eq!(result.stats.num_sequences, 2);
}

#[test]
fn test_realize_repeated_sequence_content() {
    // All sequences identical — graph should be minimal.
    let config = test_config();
    let seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 28bp of A's
    let seqs = vec![
        make_seq("s1", seq, 0, 100),
        make_seq("s2", seq, 10, 100),
        make_seq("s3", seq, 20, 100),
    ];
    let result = realize_from_sequences(&seqs, &config).unwrap();
    validate_gfa(&result);

    assert_eq!(result.stats.num_sequences, 3);
    let names = path_names(&result.gfa);
    assert_eq!(names.len(), 3);

    // All sequences identical → each path should spell the same sequence.
    for name in &names {
        let spelled = spell_path(&result.gfa, name).unwrap();
        assert_eq!(
            spelled, seq,
            "Path {} should spell the identical input sequence",
            name
        );
    }
}

