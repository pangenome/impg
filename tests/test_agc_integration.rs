#![cfg(feature = "agc")]

use impg::agc_index::AgcIndex;
use impg::faidx::FastaIndex;
use impg::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use std::io;

#[test]
fn test_agc_vs_fasta_same_content() -> io::Result<()> {
    // Test data paths
    let test_data_dir = "tests/test_data";
    let agc_file = format!("{test_data_dir}/test.agc");

    // Build AGC index
    let agc_index = AgcIndex::build_from_files(&[agc_file])?;

    // Test sequences with unambiguous names (using sample@contig format)
    // to avoid depending on library-specific ordering when contig names are duplicated
    let test_cases = vec![
        // (agc_query, fasta_file, fasta_contig, start, end)
        ("chr1@ref", "ref.fa", "chr1", 0, 10),      // First 10 bp of chr1 from ref
        ("chr1@ref", "ref.fa", "chr1", 5, 15),      // Middle section of chr1 from ref
        ("chr1@b", "b.fa", "chr1", 0, 9),           // chr1 from sample b (9 bp)
        ("chr1a", "a.fa", "chr1a", 0, 5),           // chr1a from sample a (unique)
        ("1", "c.fa", "1", 0, 10),                  // sequence "1" from sample c (unique)
    ];

    for (agc_query, fasta_file, fasta_contig, start, end) in test_cases {
        // Build a FASTA index for just this file
        let fasta_path = format!("{}/{}", test_data_dir, fasta_file);
        let fasta_index = FastaIndex::build_from_files(&[fasta_path])?;

        // Fetch from FASTA
        let fasta_seq = fasta_index.fetch_sequence(fasta_contig, start, end)?;

        // Fetch from AGC
        let agc_seq = agc_index.fetch_sequence(agc_query, start, end)?;

        // Compare sequences
        assert_eq!(
            fasta_seq, agc_seq,
            "Sequences differ for AGC:{agc_query} vs FASTA:{fasta_file}/{fasta_contig}:{start}-{end}"
        );

        // Also check that content is uppercase ASCII
        assert!(fasta_seq
            .iter()
            .all(|&b| b.is_ascii_uppercase() || b == b'N'));
        assert!(agc_seq.iter().all(|&b| b.is_ascii_uppercase() || b == b'N'));
    }

    Ok(())
}

#[test]
fn test_unified_sequence_index_fasta() -> io::Result<()> {
    let test_data_dir = "tests/test_data";
    let fasta_files = vec![
        format!("{}/ref.fa", test_data_dir),
        format!("{}/a.fa", test_data_dir),
    ];

    let index = UnifiedSequenceIndex::from_files(&fasta_files)?;

    // Test fetching
    let seq = index.fetch_sequence("chr1", 0, 10)?;
    assert_eq!(seq.len(), 10);
    assert!(seq.iter().all(|&b| b.is_ascii_uppercase() || b == b'N'));

    Ok(())
}

#[test]
fn test_unified_sequence_index_agc() -> io::Result<()> {
    let test_data_dir = "tests/test_data";
    let agc_files = vec![format!("{}/test.agc", test_data_dir)];

    let index = UnifiedSequenceIndex::from_files(&agc_files)?;

    // Test fetching with sample@contig format
    let seq = index.fetch_sequence("chr1@ref", 0, 10)?;
    assert_eq!(seq.len(), 10);
    assert!(seq.iter().all(|&b| b.is_ascii_uppercase() || b == b'N'));

    // Test fetching with just contig name (chr1a is unique)
    let seq2 = index.fetch_sequence("chr1a", 0, 5)?;
    assert_eq!(seq2.len(), 5);

    Ok(())
}

#[test]
fn test_mixed_file_types_error() {
    let test_data_dir = "tests/test_data";
    let mixed_files = vec![
        format!("{}/ref.fa", test_data_dir),
        format!("{}/test.agc", test_data_dir),
    ];

    let result = UnifiedSequenceIndex::from_files(&mixed_files);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Mixed file types"));
}
