use std::io::{self, Write};
use spoa_rs::{Graph as SpoaGraph, AlignmentEngine, AlignmentType as SpoaAlignmentType};
use crate::faidx::FastaIndex;
use crate::impg::{AdjustedInterval, Impg};

struct SequenceMetadata {
    name: String,
    start: i32,
    size: i32,
    strand: char,
    total_length: usize
}

pub fn generate_gfa_from_intervals(
    impg: &Impg,
    results: &[AdjustedInterval],
    fasta_index: &FastaIndex,
    scoring_params: (u8, u8, u8),
) -> String {
    // Prepare POA graph and sequences
    let (graph, sequence_metadata) = prepare_poa_graph_and_sequences(
        impg, results, fasta_index, scoring_params
    ).unwrap();

    // Generate headers for GFA
    let headers: Vec<String> = sequence_metadata.iter()
        .map(|meta| format!("{}:{}-{}", meta.name, 
            if meta.strand == '+' { meta.start } else { (meta.total_length as i32) - meta.start - meta.size },
            if meta.strand == '+' { meta.start + meta.size } else { (meta.total_length as i32) - meta.start }
        ))
        .collect();
    
    // Generate GFA directly from the graph
    graph.generate_gfa(&headers, false)
}

pub fn generate_maf_from_intervals(
    impg: &Impg,
    results: &[AdjustedInterval],
    fasta_index: &FastaIndex,
    scoring_params: (u8, u8, u8),
) -> String {
    // Prepare POA graph and sequences
    let (graph, sequence_metadata) = prepare_poa_graph_and_sequences(
        impg, results, fasta_index, scoring_params
    ).unwrap();

    // Generate MSA from the SPOA graph
    let msa = graph.generate_msa();
    
    // Print MAF format to stdout
    format_maf_from_msa(&msa, &sequence_metadata, None)
}

fn prepare_poa_graph_and_sequences(
    impg: &Impg,
    results: &[AdjustedInterval],
    fasta_index: &FastaIndex,
    scoring_params: (u8, u8, u8),
) -> io::Result<(SpoaGraph, Vec<SequenceMetadata>)> {
    // Create a SPOA graph
    let mut graph = SpoaGraph::new();
    
    // Create scoring parameters for alignment
    let (mismatch, gap_open, gap_extend) = scoring_params;
    
    // Create an alignment engine with affine gap penalties
    let mut engine = AlignmentEngine::new_affine(
        SpoaAlignmentType::kNW,  // Global alignment (Needleman-Wunsch)
        5,                        // match score (positive)
        -(mismatch as i8),        // mismatch penalty (negative)
        -(gap_open as i8),        // gap open penalty (negative)
        -(gap_extend as i8),      // gap extend penalty (negative)
    );

    // Collect sequences and metadata for each interval
    let mut sequence_metadata = Vec::new();
    
    for (interval, _, _) in results.iter() {
        let seq_name = impg.seq_index.get_name(interval.metadata)
            .ok_or_else(|| io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence name not found for ID {}", interval.metadata)
            ))?;
            
        // Get total sequence length
        let total_length = impg.seq_index.get_len_from_id(interval.metadata)
            .ok_or_else(|| io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence length not found for ID {}", interval.metadata)
            ))? as usize;
            
        // Determine actual start and end based on orientation
        let (start, end, strand) = if interval.first <= interval.last {
            (interval.first, interval.last, '+')
        } else {
            (interval.last, interval.first, '-')
        };
        
        // Fetch the sequence
        let sequence = fasta_index.fetch_sequence(seq_name, start, end)?;
        
        // If reverse strand, reverse complement the sequence
        let sequence = if strand == '-' {
            reverse_complement(&sequence)
        } else {
            sequence
        };
        
        let sequence_str = String::from_utf8_lossy(&sequence).to_string();
        let seq_size = (end - start) as i32;
        
        // For MAF format, if strand is "-", start is relative to reverse-complemented sequence
        let maf_start = if strand == '-' {
            (total_length as i32) - end
        } else {
            start
        };
        
        sequence_metadata.push(SequenceMetadata {
            name: seq_name.to_string(),
            start: maf_start,
            size: seq_size,
            strand,
            total_length,
        });
        
        // Add to SPOA graph
        let weights = vec![1u32; sequence_str.len()];
        let (_, alignment) = engine.align(&sequence_str, &graph);
        graph.add_alignment_with_weights(alignment, &sequence_str, &weights);
    }
    
    Ok((graph, sequence_metadata))
}

fn format_maf_from_msa(
    msa: &[String],
    sequence_metadata: &[SequenceMetadata],
    block_name: Option<String>,
) -> String {
    let mut output = String::new();
    
    // Write MAF header
    output.push_str("##maf version=1 scoring=spoa\n");
    if let Some(ref name) = block_name {
        output.push_str(&format!("# {}\n", name));
    }
    
    // Find trimming positions (remove all-gap columns at start and end)
    let msa_len = msa.first().map(|s| s.len()).unwrap_or(0);
    let mut start_trim = 0;
    let mut end_trim = msa_len;
    
    // Find first non-gap column
    'outer: for pos in 0..msa_len {
        for seq in msa {
            if seq.chars().nth(pos).unwrap_or('-') != '-' {
                start_trim = pos;
                break 'outer;
            }
        }
    }
    
    // Find last non-gap column
    'outer2: for pos in (0..msa_len).rev() {
        for seq in msa {
            if seq.chars().nth(pos).unwrap_or('-') != '-' {
                end_trim = pos + 1;
                break 'outer2;
            }
        }
    }
    
    // Write alignment block
    output.push('\n'); // blank line before block
    output.push_str("a score=0.0\n"); // We don't have a meaningful score from SPOA
    
    // Write sequence lines
    for (msa_seq, meta) in msa.iter().zip(sequence_metadata.iter()) {
        let trimmed_seq = &msa_seq[start_trim..end_trim];
        
        // Count non-gap characters to get the actual aligned size
        let aligned_size = trimmed_seq.chars().filter(|&c| c != '-').count() as i32;
        
        output.push_str(&format!("s {} {} {} {} {} {}\n",
            meta.name,
            meta.start,
            aligned_size,
            meta.strand,
            meta.total_length,
            trimmed_seq
        ));
    }
    
    output.push('\n'); // blank line after block
    output
}

/// Given a raw GFAv1.1 string and a `Write` target, convert it to GFAv1.0:
///  - Rewrite `H` lines to `H\tVN:Z:1.0`
///  - Strip leading `s` from every `S` and bump the numeric ID by 1
///  - Strip leading `s` from every `L` endpoint and bump both numeric IDs by 1
///  - Convert every `W` line into a `P` line
fn _convert_and_write_gfa<R: AsRef<str>, W: Write>(
    raw_gfa: R,
    writer: &mut W,
) -> io::Result<()> {
    for line in raw_gfa.as_ref().lines() {
        // 1) Header → force "H\tVN:Z:1.0"
        if line.starts_with("H\t") {
            writeln!(writer, "H\tVN:Z:1.0")?;
            continue;
        }

        // 2) Segment lines: S\t<old>\t<seq>[\t<TAGs>...]
        if line.starts_with("S\t") {
            // Split out “S”, old-segName, sequence, and any trailing tags
            let fields: Vec<&str> = line.splitn(4, '\t').collect();
            if fields.len() < 3 {
                // Malformed S line: echo raw
                writeln!(writer, "{}", line)?;
                continue;
            }

            let old_seg = fields[1];
            // Strip leading 's' and bump by 1 if it parses as an integer
            let new_id = if let Some(stripped) = old_seg.strip_prefix('s') {
                if let Ok(x) = stripped.parse::<i32>() {
                    (x + 1).to_string()
                } else {
                    old_seg.to_string()
                }
            } else {
                old_seg.to_string()
            };

            // Reconstruct:
            if fields.len() == 3 {
                // No extra tags
                writeln!(writer, "S\t{}\t{}", new_id, fields[2])?;
            } else {
                // fields[3] holds “all trailing content”
                writeln!(writer, "S\t{}\t{}\t{}", new_id, fields[2], fields[3])?;
            }
            continue;
        }

        // 3) Link lines: L\t<oldFrom>\t<fromOri>\t<oldTo>\t<toOri>\t<overlap>[\t<TAGs>...]
        if line.starts_with("L\t") {
            let fields: Vec<&str> = line.splitn(7, '\t').collect();
            if fields.len() < 6 {
                // Malformed L line: echo raw
                writeln!(writer, "{}", line)?;
                continue;
            }

            let old_from = fields[1];
            let old_to   = fields[3];

            let new_from = if let Some(stripped) = old_from.strip_prefix('s') {
                if let Ok(x) = stripped.parse::<i32>() {
                    (x + 1).to_string()
                } else {
                    old_from.to_string()
                }
            } else {
                old_from.to_string()
            };

            let new_to = if let Some(stripped) = old_to.strip_prefix('s') {
                if let Ok(x) = stripped.parse::<i32>() {
                    (x + 1).to_string()
                } else {
                    old_to.to_string()
                }
            } else {
                old_to.to_string()
            };

            // Reassemble:
            if fields.len() == 6 {
                // No extra tags
                writeln!(
                    writer,
                    "L\t{}\t{}\t{}\t{}\t{}",
                    new_from, fields[2], new_to, fields[4], fields[5]
                )?;
            } else {
                // fields[6] holds trailing tags
                writeln!(
                    writer,
                    "L\t{}\t{}\t{}\t{}\t{}\t{}",
                    new_from, fields[2], new_to, fields[4], fields[5], fields[6]
                )?;
            }
            continue;
        }

        // 4) Walk lines: W\t<sampleId>\t<hapIndex>\t<seqId>\t<seqStart>\t<seqEnd>\t<walkStr>[\t<TAGs>...]
        if line.starts_with("W\t") {
            let fields: Vec<&str> = line.splitn(8, '\t').collect();
            if fields.len() < 7 {
                // Malformed W line: echo raw
                writeln!(writer, "{}", line)?;
                continue;
            }

            // Use the “seqId” field (fields[3]) as the PathName
            let path_name = fields[3];
            let walk_str  = fields[6]; // e.g. “>s0>s2>s3...”

            // Parse walk_str into pairs of (orientation, segmentName),
            // strip leading 's', bump numeric part by 1, then reattach '+' or '-'.
            let mut segs_oriented = Vec::new();
            let mut chars = walk_str.chars().peekable();
            while let Some(c) = chars.next() {
                if c == '>' || c == '<' {
                    let orient = if c == '>' { '+' } else { '-' };
                    let mut segname = String::new();
                    while let Some(&next) = chars.peek() {
                        if next == '>' || next == '<' {
                            break;
                        }
                        segname.push(next);
                        chars.next();
                    }
                    // segname is e.g. “s0”. Strip 's' and bump if possible:
                    let new_seg = if let Some(stripped) = segname.strip_prefix('s') {
                        if let Ok(x) = stripped.parse::<i32>() {
                            (x + 1).to_string()
                        } else {
                            segname.clone()
                        }
                    } else {
                        segname.clone()
                    };
                    segs_oriented.push(format!("{}{}", new_seg, orient));
                }
            }

            let joined = segs_oriented.join(",");
            // Emit as a `P` line, with Overlaps = “*”
            writeln!(writer, "P\t{}\t{}\t*", path_name, joined)?;
            continue;
        }

        // 5) All other lines (including C/J/P/L/W tags not shown above, or comment lines): pass through:
        writeln!(writer, "{}", line)?;
    }

    Ok(())
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => base,
        })
        .collect()
}
