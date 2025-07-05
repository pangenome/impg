use impg::agc_index::AgcIndex;
use impg::faidx::FastaIndex;
use impg::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use std::io;

#[test]
fn test_agc_vs_fasta_same_content() -> io::Result<()> {
    // Test data paths
    let test_data_dir = "test_data";
    let fasta_files = vec![
        format!("{}/ref.fa", test_data_dir),
        format!("{}/a.fa", test_data_dir),
        format!("{}/b.fa", test_data_dir),
        format!("{}/c.fa", test_data_dir),
    ];
    let agc_file = format!("{}/test.agc", test_data_dir);

    // Build both indices
    let fasta_index = FastaIndex::build_from_files(&fasta_files)?;
    let agc_index = AgcIndex::build_from_files(&[agc_file])?;

    // Test sequences to check
    let test_cases = vec![
        ("chr1", 0, 10),      // First 10 bp of chr1
        ("chr1", 5, 15),      // Middle section  
        ("chr1a", 0, 5),      // chr1a from sample a
        ("1", 0, 10),         // sequence "1" from sample c
    ];

    for (seq_name, start, end) in test_cases {
        // Fetch from FASTA
        let fasta_seq = fasta_index.fetch_sequence(seq_name, start, end)?;
        
        // Fetch from AGC - try with just contig name first
        let agc_seq = match agc_index.fetch_sequence(seq_name, start, end) {
            Ok(seq) => seq,
            Err(_) => {
                // If contig name alone doesn't work, try different sample@contig combinations
                // This might happen if contig names are not unique across samples
                let samples = ["ref", "a", "b", "c"];
                let mut found = false;
                let mut result = Vec::new();
                
                for sample in &samples {
                    let query = format!("{}@{}", seq_name, sample);
                    if let Ok(seq) = agc_index.fetch_sequence(&query, start, end) {
                        result = seq;
                        found = true;
                        break;
                    }
                }
                
                if !found {
                    return Err(io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("Sequence {} not found in AGC", seq_name),
                    ));
                }
                result
            }
        };

        // Compare sequences
        assert_eq!(
            fasta_seq, agc_seq,
            "Sequences differ for {}:{}-{}", seq_name, start, end
        );
        
        // Also check that content is uppercase ASCII
        assert!(fasta_seq.iter().all(|&b| b.is_ascii_uppercase() || b == b'N'));
        assert!(agc_seq.iter().all(|&b| b.is_ascii_uppercase() || b == b'N'));
    }

    Ok(())
}

#[test]
fn test_unified_sequence_index_fasta() -> io::Result<()> {
    let test_data_dir = "test_data";
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
    let test_data_dir = "test_data";
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
    let test_data_dir = "test_data";
    let mixed_files = vec![
        format!("{}/ref.fa", test_data_dir),
        format!("{}/test.agc", test_data_dir),
    ];

    let result = UnifiedSequenceIndex::from_files(&mixed_files);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Mixed file types"));
}