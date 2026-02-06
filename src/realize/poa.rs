use crate::graph::{reverse_complement, SequenceMetadata};
use crate::impg_index::ImpgIndex;
use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use coitrees::Interval;
use spoa_rs::{AlignmentEngine, AlignmentType as SpoaAlignmentType, Graph as SpoaGraph};
use std::io;

/// Result of padded POA: a GFA string with padding trimmed away.
pub struct PaddedPoaResult {
    /// GFA graph as a string (padding regions removed)
    pub gfa: String,
    /// Metadata for each sequence (original, unpadded coordinates)
    pub sequence_metadata: Vec<SequenceMetadata>,
}

/// Information about how a sequence was padded.
struct PaddedSequence {
    /// The padded sequence string (including flanking context)
    sequence: String,
    /// Metadata for the original (unpadded) region
    metadata: SequenceMetadata,
    /// Number of bases of left padding actually fetched
    left_pad: usize,
    /// Number of bases of right padding actually fetched
    right_pad: usize,
}

/// Run SPOA partial order alignment with boundary padding, then trim padding
/// from the resulting GFA.
///
/// # Algorithm
///
/// 1. For each interval, extend coordinates by `padding` bp on each side
///    (clamped to contig boundaries).
/// 2. Fetch extended sequences and track padding amounts.
/// 3. Build SPOA graph from padded sequences (longest first).
/// 4. Generate MSA from the SPOA graph to identify padding columns.
/// 5. Determine which MSA columns are "core" (non-padding) vs "padding-only".
/// 6. Regenerate GFA from only the core portion of the alignment.
pub fn padded_poa(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
    padding: usize,
) -> io::Result<PaddedPoaResult> {
    if results.is_empty() {
        return Ok(PaddedPoaResult {
            gfa: String::new(),
            sequence_metadata: Vec::new(),
        });
    }

    // Step 1-2: Fetch padded sequences
    let mut padded_sequences = fetch_padded_sequences(impg, results, sequence_index, padding)?;

    // Sort by total (padded) sequence length, longest first — best for SPOA quality
    padded_sequences.sort_by(|a, b| b.sequence.len().cmp(&a.sequence.len()));

    // Step 3: Build SPOA graph
    let (match_score, mismatch, gap_open1, gap_extend1, gap_open2, gap_extend2) = scoring_params;
    let mut graph = SpoaGraph::new();
    let mut engine = AlignmentEngine::new_convex(
        SpoaAlignmentType::kSW,
        match_score as i8,
        -(mismatch as i8),
        -(gap_open1 as i8),
        -(gap_extend1 as i8),
        -(gap_open2 as i8),
        -(gap_extend2 as i8),
    );

    for ps in &padded_sequences {
        let weights = vec![1u32; ps.sequence.len()];
        let (_, alignment) = engine.align(&ps.sequence, &graph);
        graph.add_alignment_with_weights(alignment, &ps.sequence, &weights);
    }

    // If no padding was applied (all zeros), skip trimming
    let any_padding = padded_sequences
        .iter()
        .any(|ps| ps.left_pad > 0 || ps.right_pad > 0);

    if !any_padding {
        // No padding — generate GFA directly (same as existing generate_gfa_from_intervals)
        let headers = make_headers(&padded_sequences);
        let gfa = graph.generate_gfa(&headers, false);
        let metadata: Vec<SequenceMetadata> =
            padded_sequences.into_iter().map(|ps| ps.metadata).collect();
        let gfa = post_process_gfa_for_strands(gfa, &metadata);
        return Ok(PaddedPoaResult {
            gfa,
            sequence_metadata: metadata,
        });
    }

    // Step 4: Generate MSA to find padding column boundaries
    let msa = graph.generate_msa();

    // Step 5: Find the core (non-padding) column range
    let (core_start, core_end) = find_core_column_range(&msa, &padded_sequences);

    // Step 6: Build a new SPOA graph from only the trimmed (core) sequences
    let trimmed_gfa =
        build_trimmed_gfa(&msa, &padded_sequences, core_start, core_end, scoring_params)?;

    let metadata: Vec<SequenceMetadata> =
        padded_sequences.into_iter().map(|ps| ps.metadata).collect();

    Ok(PaddedPoaResult {
        gfa: trimmed_gfa,
        sequence_metadata: metadata,
    })
}

/// Run padded POA directly on pre-prepared sequences (no impg/index needed).
/// This is the interface used by the recursive realize engine, which already
/// has extracted sequences.
pub fn padded_poa_from_sequences(
    sequences: &[(String, SequenceMetadata)],
    scoring_params: (u8, u8, u8, u8, u8, u8),
    padding: usize,
) -> io::Result<PaddedPoaResult> {
    if sequences.is_empty() {
        return Ok(PaddedPoaResult {
            gfa: String::new(),
            sequence_metadata: Vec::new(),
        });
    }

    // When called from the realize engine, sequences are already extracted at
    // the correct coordinates (possibly with padding applied upstream).
    // The `padding` parameter here is informational — if the caller already
    // padded, they pass padding=0 here. If they want us to do column-level
    // trimming on already-padded data, they provide the padding amount and
    // we use the MSA approach.
    //
    // For the common case in recursive realize: sequences come pre-padded
    // and we need to trim `padding` bases from each end in MSA space.

    let padded_sequences: Vec<PaddedSequence> = sequences
        .iter()
        .map(|(seq, meta)| {
            // Assume symmetric padding was applied upstream
            let left_pad = padding.min(seq.len() / 2);
            let right_pad = padding.min(seq.len() - left_pad);
            PaddedSequence {
                sequence: seq.clone(),
                metadata: meta.clone(),
                left_pad,
                right_pad,
            }
        })
        .collect();

    // Build SPOA graph
    let (match_score, mismatch, gap_open1, gap_extend1, gap_open2, gap_extend2) = scoring_params;
    let mut graph = SpoaGraph::new();
    let mut engine = AlignmentEngine::new_convex(
        SpoaAlignmentType::kSW,
        match_score as i8,
        -(mismatch as i8),
        -(gap_open1 as i8),
        -(gap_extend1 as i8),
        -(gap_open2 as i8),
        -(gap_extend2 as i8),
    );

    for ps in &padded_sequences {
        let weights = vec![1u32; ps.sequence.len()];
        let (_, alignment) = engine.align(&ps.sequence, &graph);
        graph.add_alignment_with_weights(alignment, &ps.sequence, &weights);
    }

    let any_padding = padded_sequences
        .iter()
        .any(|ps| ps.left_pad > 0 || ps.right_pad > 0);

    if !any_padding {
        let headers = make_headers(&padded_sequences);
        let gfa = graph.generate_gfa(&headers, false);
        let metadata: Vec<SequenceMetadata> =
            padded_sequences.into_iter().map(|ps| ps.metadata).collect();
        let gfa = post_process_gfa_for_strands(gfa, &metadata);
        return Ok(PaddedPoaResult {
            gfa,
            sequence_metadata: metadata,
        });
    }

    let msa = graph.generate_msa();
    let (core_start, core_end) = find_core_column_range(&msa, &padded_sequences);

    let trimmed_gfa =
        build_trimmed_gfa(&msa, &padded_sequences, core_start, core_end, scoring_params)?;

    let metadata: Vec<SequenceMetadata> =
        padded_sequences.into_iter().map(|ps| ps.metadata).collect();

    Ok(PaddedPoaResult {
        gfa: trimmed_gfa,
        sequence_metadata: metadata,
    })
}

/// Fetch sequences with padding extended on each side, clamped to contig boundaries.
fn fetch_padded_sequences(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    padding: usize,
) -> io::Result<Vec<PaddedSequence>> {
    use rayon::prelude::*;

    let padding_i32 = padding as i32;

    results
        .par_iter()
        .map(|interval| -> io::Result<PaddedSequence> {
            let seq_name = impg
                .seq_index()
                .get_name(interval.metadata)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("Sequence name not found for ID {}", interval.metadata),
                    )
                })?;

            let total_length = impg
                .seq_index()
                .get_len_from_id(interval.metadata)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("Sequence length not found for ID {}", interval.metadata),
                    )
                })?;

            // Determine orientation
            let (start, end, strand) = if interval.first <= interval.last {
                (interval.first, interval.last, '+')
            } else {
                (interval.last, interval.first, '-')
            };

            // Extend coordinates by padding, clamped to [0, total_length]
            let padded_start = (start - padding_i32).max(0);
            let padded_end = (end + padding_i32).min(total_length as i32);

            let actual_left_pad = (start - padded_start) as usize;
            let actual_right_pad = (padded_end - end) as usize;

            // Fetch padded sequence
            let sequence = sequence_index.fetch_sequence(seq_name, padded_start, padded_end)?;

            // Reverse complement if needed
            // Note: for RC, left/right padding swap because the sequence is reversed
            let (sequence, left_pad, right_pad) = if strand == '-' {
                (
                    reverse_complement(&sequence),
                    actual_right_pad, // after RC, the right pad becomes left
                    actual_left_pad,  // and left becomes right
                )
            } else {
                (sequence, actual_left_pad, actual_right_pad)
            };

            let sequence_str = String::from_utf8_lossy(&sequence).to_string();
            let seq_size = end - start;

            // MAF-style start for the ORIGINAL (unpadded) region
            let maf_start = if strand == '-' {
                (total_length as i32) - end
            } else {
                start
            };

            let metadata = SequenceMetadata {
                name: seq_name.to_string(),
                start: maf_start,
                size: seq_size,
                strand,
                total_length,
            };

            Ok(PaddedSequence {
                sequence: sequence_str,
                metadata,
                left_pad,
                right_pad,
            })
        })
        .collect()
}

/// Determine the MSA column range that corresponds to the core (non-padding) region.
///
/// For each sequence in the MSA, we know how many bases of left/right padding
/// it has. We walk the MSA columns, counting non-gap characters per sequence,
/// to find:
/// - `core_start`: the first column where ALL sequences have consumed their left padding
/// - `core_end`: the last column before ANY sequence enters its right padding
///
/// This ensures the core region contains only "real" (non-padding) aligned bases.
fn find_core_column_range(msa: &[String], padded_sequences: &[PaddedSequence]) -> (usize, usize) {
    let ncols = msa.first().map(|s| s.len()).unwrap_or(0);
    if ncols == 0 {
        return (0, 0);
    }

    let msa_bytes: Vec<&[u8]> = msa.iter().map(|s| s.as_bytes()).collect();
    let nseqs = msa_bytes.len().min(padded_sequences.len());

    // For each sequence, count how many non-gap bases appear in each column.
    // We need to find the column where sequence i has consumed left_pad[i] non-gap chars
    // (i.e., the left padding is "used up") and the column where it would start consuming
    // right padding chars.

    // Phase 1: find core_start — the maximum over all sequences of the column
    // where left padding ends.
    let mut core_start = 0;
    for i in 0..nseqs {
        let left_pad = padded_sequences[i].left_pad;
        if left_pad == 0 {
            continue;
        }
        let mut non_gap_count = 0;
        for col in 0..ncols {
            if msa_bytes[i][col] != b'-' {
                non_gap_count += 1;
                if non_gap_count == left_pad {
                    // The next column is where core starts for this sequence
                    core_start = core_start.max(col + 1);
                    break;
                }
            }
        }
    }

    // Phase 2: find core_end — the minimum over all sequences of the column
    // where right padding begins.
    let mut core_end = ncols;
    for i in 0..nseqs {
        let right_pad = padded_sequences[i].right_pad;
        if right_pad == 0 {
            continue;
        }
        // Total non-gap bases in this MSA row
        let total_bases: usize = msa_bytes[i].iter().filter(|&&b| b != b'-').count();
        // The core bases end at position (total_bases - right_pad) in non-gap counting
        let core_base_count = total_bases.saturating_sub(right_pad);
        if core_base_count == 0 {
            core_end = 0;
            break;
        }

        let mut non_gap_count = 0;
        for col in 0..ncols {
            if msa_bytes[i][col] != b'-' {
                non_gap_count += 1;
                if non_gap_count == core_base_count {
                    // This column is the last core base; core_end is the next column
                    core_end = core_end.min(col + 1);
                    break;
                }
            }
        }
    }

    // Ensure valid range
    if core_start >= core_end {
        // Degenerate case: padding consumed everything. Return the middle.
        let mid = ncols / 2;
        return (mid, mid);
    }

    (core_start, core_end)
}

/// Build a trimmed GFA by extracting core (non-padding) subsequences from
/// the MSA and re-running POA on just those.
///
/// This approach is clean: we extract only the bases each sequence contributes
/// to the core alignment columns, then build a fresh SPOA graph. This avoids
/// the complexity of trying to surgically remove nodes from a GFA.
fn build_trimmed_gfa(
    msa: &[String],
    padded_sequences: &[PaddedSequence],
    core_start: usize,
    core_end: usize,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<String> {
    if core_start >= core_end {
        return Ok(String::new());
    }

    let msa_bytes: Vec<&[u8]> = msa.iter().map(|s| s.as_bytes()).collect();
    let nseqs = msa_bytes.len().min(padded_sequences.len());

    // Extract core subsequences: for each MSA row, take only non-gap characters
    // from columns [core_start, core_end)
    let mut core_sequences: Vec<String> = Vec::with_capacity(nseqs);
    for i in 0..nseqs {
        let core_seq: String = msa_bytes[i][core_start..core_end]
            .iter()
            .filter(|&&b| b != b'-')
            .map(|&b| b as char)
            .collect();
        core_sequences.push(core_seq);
    }

    // Build fresh SPOA graph from core sequences
    let (match_score, mismatch, gap_open1, gap_extend1, gap_open2, gap_extend2) = scoring_params;
    let mut graph = SpoaGraph::new();
    let mut engine = AlignmentEngine::new_convex(
        SpoaAlignmentType::kSW,
        match_score as i8,
        -(mismatch as i8),
        -(gap_open1 as i8),
        -(gap_extend1 as i8),
        -(gap_open2 as i8),
        -(gap_extend2 as i8),
    );

    // Keep original order (already sorted longest-first)
    for seq in &core_sequences {
        if seq.is_empty() {
            continue;
        }
        let weights = vec![1u32; seq.len()];
        let (_, alignment) = engine.align(seq, &graph);
        graph.add_alignment_with_weights(alignment, seq, &weights);
    }

    // Generate GFA with headers matching the original (unpadded) metadata
    let headers = make_headers(padded_sequences);
    let gfa = graph.generate_gfa(&headers, false);

    // Post-process for strand information
    let metadata: Vec<SequenceMetadata> =
        padded_sequences.iter().map(|ps| ps.metadata.clone()).collect();
    Ok(post_process_gfa_for_strands(gfa, &metadata))
}

/// Build GFA path headers from sequence metadata.
fn make_headers(padded_sequences: &[PaddedSequence]) -> Vec<String> {
    padded_sequences
        .iter()
        .map(|ps| {
            let meta = &ps.metadata;
            format!(
                "{}:{}-{}",
                meta.name,
                if meta.strand == '+' {
                    meta.start
                } else {
                    (meta.total_length as i32) - meta.start - meta.size
                },
                if meta.strand == '+' {
                    meta.start + meta.size
                } else {
                    (meta.total_length as i32) - meta.start
                }
            )
        })
        .collect()
}

/// Post-process GFA to handle strand information.
/// Reverse strand paths need their segments reversed and orientations flipped.
fn post_process_gfa_for_strands(gfa: String, metadata: &[SequenceMetadata]) -> String {
    let mut output = String::new();

    // Build a map of header -> strand
    let strand_map: std::collections::HashMap<String, char> = metadata
        .iter()
        .map(|meta| {
            let header = format!(
                "{}:{}-{}",
                meta.name,
                if meta.strand == '+' {
                    meta.start
                } else {
                    (meta.total_length as i32) - meta.start - meta.size
                },
                if meta.strand == '+' {
                    meta.start + meta.size
                } else {
                    (meta.total_length as i32) - meta.start
                }
            );
            (header, meta.strand)
        })
        .collect();

    for line in gfa.lines() {
        if line.starts_with("P\t") {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let path_name = parts[1];
                let path_segments = parts[2];

                if let Some(&strand) = strand_map.get(path_name) {
                    if strand == '-' {
                        let segments: Vec<&str> = path_segments.split(',').collect();
                        let reversed_path: Vec<String> = segments
                            .iter()
                            .rev()
                            .map(|seg| {
                                if let Some(seg_stripped) = seg.strip_suffix('+') {
                                    format!("{seg_stripped}-")
                                } else if let Some(seg_stripped) = seg.strip_suffix('-') {
                                    format!("{seg_stripped}+")
                                } else {
                                    seg.to_string()
                                }
                            })
                            .collect();

                        output.push_str(&format!(
                            "P\t{}\t{}",
                            path_name,
                            reversed_path.join(",")
                        ));
                        if parts.len() > 3 {
                            output.push('\t');
                            output.push_str(&parts[3..].join("\t"));
                        }
                        output.push('\n');
                        continue;
                    }
                }
            }
        }

        output.push_str(line);
        output.push('\n');
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_core_column_range_no_padding() {
        // All padding is 0 → core is entire MSA
        let msa = vec!["ACGTACGT".to_string(), "ACGTACGT".to_string()];
        let padded = vec![
            PaddedSequence {
                sequence: "ACGTACGT".to_string(),
                metadata: SequenceMetadata {
                    name: "s1".to_string(),
                    start: 0,
                    size: 8,
                    strand: '+',
                    total_length: 100,
                },
                left_pad: 0,
                right_pad: 0,
            },
            PaddedSequence {
                sequence: "ACGTACGT".to_string(),
                metadata: SequenceMetadata {
                    name: "s2".to_string(),
                    start: 0,
                    size: 8,
                    strand: '+',
                    total_length: 100,
                },
                left_pad: 0,
                right_pad: 0,
            },
        ];

        let (start, end) = find_core_column_range(&msa, &padded);
        assert_eq!(start, 0);
        assert_eq!(end, 8);
    }

    #[test]
    fn test_find_core_column_range_symmetric_padding() {
        // Each sequence has 2bp left pad and 2bp right pad
        // MSA: "PPCCCCPP" (P=padding, C=core) — no gaps for simplicity
        let msa = vec!["AACCCCTT".to_string(), "AACCCCTT".to_string()];
        let padded = vec![
            PaddedSequence {
                sequence: "AACCCCTT".to_string(),
                metadata: SequenceMetadata {
                    name: "s1".to_string(),
                    start: 2,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
                left_pad: 2,
                right_pad: 2,
            },
            PaddedSequence {
                sequence: "AACCCCTT".to_string(),
                metadata: SequenceMetadata {
                    name: "s2".to_string(),
                    start: 2,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
                left_pad: 2,
                right_pad: 2,
            },
        ];

        let (start, end) = find_core_column_range(&msa, &padded);
        assert_eq!(start, 2); // after 2 non-gap chars of padding
        assert_eq!(end, 6); // before last 2 non-gap chars of padding
    }

    #[test]
    fn test_find_core_column_range_with_gaps() {
        // Seq 1: "AA--CCCCTT" (left_pad=2, right_pad=2)
        // Seq 2: "AAGGCCCC--" (left_pad=2, right_pad=0)
        // For seq1: 2 non-gap left padding bases consumed at col 0,1 → core starts at col 2
        // For seq2: 2 non-gap left padding bases consumed at col 0,1 → core starts at col 2
        // But seq2's gap at col 2,3 means we need to look at non-gap columns carefully
        let msa = vec!["AA--CCCCTT".to_string(), "AAGGCCCC--".to_string()];
        let padded = vec![
            PaddedSequence {
                sequence: "AACCCCTT".to_string(),
                metadata: SequenceMetadata {
                    name: "s1".to_string(),
                    start: 2,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
                left_pad: 2,
                right_pad: 2,
            },
            PaddedSequence {
                sequence: "AAGGCCCC".to_string(),
                metadata: SequenceMetadata {
                    name: "s2".to_string(),
                    start: 2,
                    size: 6,
                    strand: '+',
                    total_length: 100,
                },
                left_pad: 2,
                right_pad: 0,
            },
        ];

        let (start, end) = find_core_column_range(&msa, &padded);
        // Seq1 left pad: 2 non-gap chars → cols 0,1 → core starts at 2
        // Seq2 left pad: 2 non-gap chars → cols 0,1 → core starts at 2
        // core_start = max(2, 2) = 2
        assert_eq!(start, 2);

        // Seq1 right pad: 2 bases; total non-gap = 8; core_base_count = 6
        //   non-gap counts: col0=1, col1=2, col4=3, col5=4, col6=5, col7=6 → core_end from seq1 = 8
        // Seq2 right pad: 0 → core_end not constrained
        // core_end = min(8, 10) = 8
        assert_eq!(end, 8);
    }

    #[test]
    fn test_find_core_column_range_asymmetric_padding() {
        // Seq1 has 3bp left pad, 1bp right pad
        // Seq2 has 1bp left pad, 3bp right pad
        let msa = vec!["AAACCCCT".to_string(), "ACCCCGGG".to_string()];
        let padded = vec![
            PaddedSequence {
                sequence: "AAACCCCT".to_string(),
                metadata: SequenceMetadata {
                    name: "s1".to_string(),
                    start: 3,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
                left_pad: 3,
                right_pad: 1,
            },
            PaddedSequence {
                sequence: "ACCCCGGG".to_string(),
                metadata: SequenceMetadata {
                    name: "s2".to_string(),
                    start: 1,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
                left_pad: 1,
                right_pad: 3,
            },
        ];

        let (start, end) = find_core_column_range(&msa, &padded);
        // Seq1 left: 3 non-gap → core starts at col 3
        // Seq2 left: 1 non-gap → core starts at col 1
        // core_start = max(3, 1) = 3
        assert_eq!(start, 3);

        // Seq1 right: 1 base → total=8, core_count=7 → 7th non-gap at col 6 → core_end=7
        // Seq2 right: 3 bases → total=8, core_count=5 → 5th non-gap at col 4 → core_end=5
        // core_end = min(7, 5) = 5
        assert_eq!(end, 5);
    }

    #[test]
    fn test_find_core_column_range_empty_msa() {
        let msa: Vec<String> = vec![];
        let padded: Vec<PaddedSequence> = vec![];
        let (start, end) = find_core_column_range(&msa, &padded);
        assert_eq!(start, 0);
        assert_eq!(end, 0);
    }

    #[test]
    fn test_padded_poa_from_sequences_no_padding() {
        // Two identical sequences, no padding → should produce valid GFA
        let sequences = vec![
            (
                "ACGTACGTACGT".to_string(),
                SequenceMetadata {
                    name: "seq1".to_string(),
                    start: 0,
                    size: 12,
                    strand: '+',
                    total_length: 100,
                },
            ),
            (
                "ACGTACGTACGT".to_string(),
                SequenceMetadata {
                    name: "seq2".to_string(),
                    start: 10,
                    size: 12,
                    strand: '+',
                    total_length: 100,
                },
            ),
        ];

        let scoring = (1, 4, 6, 2, 26, 1);
        let result = padded_poa_from_sequences(&sequences, scoring, 0).unwrap();

        // GFA should contain S (segment) and P (path) lines
        assert!(result.gfa.contains("S\t"));
        assert!(result.gfa.contains("P\t"));
        // Should have paths for both sequences
        assert!(result.gfa.contains("seq1:0-12"));
        assert!(result.gfa.contains("seq2:10-22"));
        assert_eq!(result.sequence_metadata.len(), 2);
    }

    #[test]
    fn test_padded_poa_from_sequences_with_padding() {
        // Core region is "CCCC" with 2bp padding on each side
        // Seq1: "AA" + "CCCC" + "TT" = "AACCCCT T"
        // Seq2: "GG" + "CCCC" + "AA" = "GGCCCCAA"
        // After POA with padding, the trimmed result should only reflect the core "CCCC"
        let sequences = vec![
            (
                "AACCCCTT".to_string(),
                SequenceMetadata {
                    name: "seq1".to_string(),
                    start: 2,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
            ),
            (
                "GGCCCCAA".to_string(),
                SequenceMetadata {
                    name: "seq2".to_string(),
                    start: 2,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
            ),
        ];

        let scoring = (1, 4, 6, 2, 26, 1);
        let result = padded_poa_from_sequences(&sequences, scoring, 2).unwrap();

        // Should produce valid GFA
        assert!(result.gfa.contains("S\t"));
        assert!(result.gfa.contains("P\t"));
        assert_eq!(result.sequence_metadata.len(), 2);
    }

    #[test]
    fn test_padded_poa_from_sequences_with_variation() {
        // Sequences have a SNP in the core region, with padding context
        // Seq1: "AA" + "CCGCC" + "TT" = "AACCGCCTT"
        // Seq2: "AA" + "CCACC" + "TT" = "AACCACCTT"
        // The padding helps SPOA align the flanking context correctly
        let sequences = vec![
            (
                "AACCGCCTT".to_string(),
                SequenceMetadata {
                    name: "seq1".to_string(),
                    start: 2,
                    size: 5,
                    strand: '+',
                    total_length: 100,
                },
            ),
            (
                "AACCACCTT".to_string(),
                SequenceMetadata {
                    name: "seq2".to_string(),
                    start: 2,
                    size: 5,
                    strand: '+',
                    total_length: 100,
                },
            ),
        ];

        let scoring = (1, 4, 6, 2, 26, 1);
        let result = padded_poa_from_sequences(&sequences, scoring, 2).unwrap();

        // Should produce a valid bubble graph reflecting the SNP
        assert!(result.gfa.contains("S\t"));
        assert!(result.gfa.contains("P\t"));
        // With a SNP, we expect more than one segment
        let segment_count = result.gfa.lines().filter(|l| l.starts_with("S\t")).count();
        assert!(
            segment_count >= 2,
            "Expected multiple segments for SNP, got {segment_count}"
        );
    }

    #[test]
    fn test_padded_poa_empty_input() {
        let sequences: Vec<(String, SequenceMetadata)> = vec![];
        let scoring = (1, 4, 6, 2, 26, 1);
        let result = padded_poa_from_sequences(&sequences, scoring, 10).unwrap();
        assert!(result.gfa.is_empty());
        assert!(result.sequence_metadata.is_empty());
    }

    #[test]
    fn test_padded_poa_single_sequence() {
        // Single sequence → trivial graph
        let sequences = vec![(
            "AACCCCTT".to_string(),
            SequenceMetadata {
                name: "seq1".to_string(),
                start: 2,
                size: 4,
                strand: '+',
                total_length: 100,
            },
        )];

        let scoring = (1, 4, 6, 2, 26, 1);
        let result = padded_poa_from_sequences(&sequences, scoring, 2).unwrap();
        assert!(result.gfa.contains("S\t"));
        assert!(result.gfa.contains("P\t"));
        assert_eq!(result.sequence_metadata.len(), 1);
    }

    #[test]
    fn test_padded_poa_reverse_strand_metadata() {
        // Reverse strand sequence — check that metadata is preserved correctly
        let sequences = vec![
            (
                "AACCCCTT".to_string(),
                SequenceMetadata {
                    name: "seq1".to_string(),
                    start: 10,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
            ),
            (
                "AACCCCTT".to_string(),
                SequenceMetadata {
                    name: "seq2".to_string(),
                    start: 20,
                    size: 4,
                    strand: '-',
                    total_length: 100,
                },
            ),
        ];

        let scoring = (1, 4, 6, 2, 26, 1);
        let result = padded_poa_from_sequences(&sequences, scoring, 2).unwrap();
        assert!(result.gfa.contains("S\t"));
        assert!(result.gfa.contains("P\t"));
        assert_eq!(result.sequence_metadata.len(), 2);
        assert_eq!(result.sequence_metadata[0].strand, '+');
        assert_eq!(result.sequence_metadata[1].strand, '-');
    }

    #[test]
    fn test_padded_poa_large_padding_clamped() {
        // Padding larger than sequence → clamped to half the sequence length
        let sequences = vec![
            (
                "ACGT".to_string(),
                SequenceMetadata {
                    name: "seq1".to_string(),
                    start: 0,
                    size: 4,
                    strand: '+',
                    total_length: 4,
                },
            ),
            (
                "ACGT".to_string(),
                SequenceMetadata {
                    name: "seq2".to_string(),
                    start: 0,
                    size: 4,
                    strand: '+',
                    total_length: 4,
                },
            ),
        ];

        let scoring = (1, 4, 6, 2, 26, 1);
        // padding=100 but sequences are only 4bp — should clamp gracefully
        let result = padded_poa_from_sequences(&sequences, scoring, 100).unwrap();
        // With padding clamped to half of 4 = 2, and right_pad = min(100, 4-2) = 2,
        // the entire sequence is considered padding. The core range will be degenerate.
        // The function should handle this gracefully without panicking.
        assert!(result.sequence_metadata.len() == 2);
    }

    #[test]
    fn test_find_core_column_range_all_padding() {
        // Edge case: padding consumes everything → degenerate range at midpoint
        let msa = vec!["AAAA".to_string(), "AAAA".to_string()];
        let padded = vec![
            PaddedSequence {
                sequence: "AAAA".to_string(),
                metadata: SequenceMetadata {
                    name: "s1".to_string(),
                    start: 0,
                    size: 0,
                    strand: '+',
                    total_length: 4,
                },
                left_pad: 2,
                right_pad: 2,
            },
            PaddedSequence {
                sequence: "AAAA".to_string(),
                metadata: SequenceMetadata {
                    name: "s2".to_string(),
                    start: 0,
                    size: 0,
                    strand: '+',
                    total_length: 4,
                },
                left_pad: 2,
                right_pad: 2,
            },
        ];

        let (start, end) = find_core_column_range(&msa, &padded);
        // Degenerate: core_start >= core_end, returns midpoint
        assert_eq!(start, end);
    }

    #[test]
    fn test_padded_poa_three_sequences_with_indel() {
        // Three sequences with an indel in the core region
        // Seq1: pad("AA") + core("CCGCC") + pad("TT")
        // Seq2: pad("AA") + core("CC-CC") + pad("TT") → "AACCCCTT" (deletion)
        // Seq3: pad("AA") + core("CCAGCC") + pad("TT") → "AACCAGCCTT" (insertion)
        let sequences = vec![
            (
                "AACCGCCTT".to_string(),
                SequenceMetadata {
                    name: "seq1".to_string(),
                    start: 2,
                    size: 5,
                    strand: '+',
                    total_length: 100,
                },
            ),
            (
                "AACCCCTT".to_string(),
                SequenceMetadata {
                    name: "seq2".to_string(),
                    start: 2,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
            ),
            (
                "AACCAGCCTT".to_string(),
                SequenceMetadata {
                    name: "seq3".to_string(),
                    start: 2,
                    size: 6,
                    strand: '+',
                    total_length: 100,
                },
            ),
        ];

        let scoring = (1, 4, 6, 2, 26, 1);
        let result = padded_poa_from_sequences(&sequences, scoring, 2).unwrap();
        assert!(result.gfa.contains("S\t"));
        assert!(result.gfa.contains("P\t"));
        assert_eq!(result.sequence_metadata.len(), 3);
        // With indels, we expect multiple segments forming a bubble
        let segment_count = result.gfa.lines().filter(|l| l.starts_with("S\t")).count();
        assert!(
            segment_count >= 2,
            "Expected multiple segments for indel variation, got {segment_count}"
        );
    }

    #[test]
    fn test_post_process_gfa_for_strands_forward() {
        // Forward strand paths should remain unchanged
        let gfa = "S\t1\tACGT\nP\tseq1:0-4\t1+\t*\n".to_string();
        let metadata = vec![SequenceMetadata {
            name: "seq1".to_string(),
            start: 0,
            size: 4,
            strand: '+',
            total_length: 100,
        }];
        let result = post_process_gfa_for_strands(gfa.clone(), &metadata);
        assert_eq!(result, gfa);
    }

    #[test]
    fn test_post_process_gfa_for_strands_reverse() {
        // Reverse strand paths should have segments reversed and orientations flipped
        // For metadata: name="seq1", start=86, size=14, strand='-', total_length=100
        // The header is: seq1:(100-86-14)-(100-86) = seq1:0-14
        let gfa = "S\t1\tACGT\nS\t2\tTTTT\nP\tseq1:0-14\t1+,2+\t*\n".to_string();
        let metadata = vec![SequenceMetadata {
            name: "seq1".to_string(),
            start: 86,
            size: 14,
            strand: '-',
            total_length: 100,
        }];
        let result = post_process_gfa_for_strands(gfa, &metadata);
        // Path should be reversed: 2-,1-
        assert!(result.contains("2-,1-"), "Expected reversed path, got: {result}");
    }

    #[test]
    fn test_make_headers_forward_and_reverse() {
        let padded = vec![
            PaddedSequence {
                sequence: "ACGT".to_string(),
                metadata: SequenceMetadata {
                    name: "chr1".to_string(),
                    start: 10,
                    size: 4,
                    strand: '+',
                    total_length: 100,
                },
                left_pad: 0,
                right_pad: 0,
            },
            PaddedSequence {
                sequence: "ACGT".to_string(),
                metadata: SequenceMetadata {
                    name: "chr2".to_string(),
                    start: 20,
                    size: 4,
                    strand: '-',
                    total_length: 100,
                },
                left_pad: 0,
                right_pad: 0,
            },
        ];

        let headers = make_headers(&padded);
        assert_eq!(headers[0], "chr1:10-14");
        // For reverse strand: start=20, size=4, total=100
        // header = name:(total - start - size)-(total - start) = chr2:76-80
        assert_eq!(headers[1], "chr2:76-80");
    }
}
