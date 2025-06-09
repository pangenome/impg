use crate::faidx::FastaIndex;
use crate::impg::Impg;
use coitrees::Interval;
use spoa_rs::{AlignmentEngine, AlignmentType as SpoaAlignmentType, Graph as SpoaGraph};
use std::io::{self, Write};

struct SequenceMetadata {
    name: String,
    start: i32,
    size: i32,
    strand: char,
    total_length: usize,
}

pub fn generate_gfa_from_intervals(
    impg: &Impg,
    results: &[Interval<u32>],
    fasta_index: &FastaIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> String {
    // Prepare POA graph and sequences
    let (graph, sequence_metadata) =
        prepare_poa_graph_and_sequences(impg, results, fasta_index, scoring_params).unwrap();

    // Generate headers for GFA
    let headers: Vec<String> = sequence_metadata
        .iter()
        .map(|meta| {
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
        .collect();

    // Generate GFA directly from the graph
    let gfa_output = graph.generate_gfa(&headers, false);

    // Post-process GFA to handle strand information
    post_process_gfa_for_strands(gfa_output, &sequence_metadata)
}

fn post_process_gfa_for_strands(gfa: String, sequence_metadata: &[SequenceMetadata]) -> String {
    let mut output = String::new();

    // Create a map of sequence names to their strand information
    let strand_map: std::collections::HashMap<String, char> = sequence_metadata
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

    // Process each line
    for line in gfa.lines() {
        if line.starts_with("P\t") {
            // Path line - need to adjust orientations for reverse strand sequences
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let path_name = parts[1];
                let path_segments = parts[2];

                // Check if this path corresponds to a reverse strand sequence
                if let Some(&strand) = strand_map.get(path_name) {
                    if strand == '-' {
                        // Reverse the path and flip all orientations
                        let segments: Vec<&str> = path_segments.split(',').collect();
                        let reversed_path: Vec<String> = segments
                            .iter()
                            .rev()
                            .map(|seg| {
                                if let Some(seg_stripped) = seg.strip_suffix('+') {
                                    format!("{}-", seg_stripped)
                                } else if let Some(seg_stripped) = seg.strip_suffix('-') {
                                    format!("{}+", seg_stripped)
                                } else {
                                    panic!("Missing segment orientation in path: {}", path_name);
                                }
                            })
                            .collect();

                        output.push_str(&format!("P\t{}\t{}", path_name, reversed_path.join(",")));
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

        // For all other lines (or unmodified path lines), output as-is
        output.push_str(line);
        output.push('\n');
    }

    output
}

pub fn generate_maf_from_intervals(
    impg: &Impg,
    results: &[Interval<u32>],
    fasta_index: &FastaIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> String {
    // Prepare POA graph and sequences
    let (graph, sequence_metadata) =
        prepare_poa_graph_and_sequences(impg, results, fasta_index, scoring_params).unwrap();

    // Generate MSA from the SPOA graph
    let msa = graph.generate_msa();

    // Print MAF format to stdout
    format_maf_from_msa(&msa, &sequence_metadata, None)
}

fn prepare_poa_graph_and_sequences(
    impg: &Impg,
    results: &[Interval<u32>],
    fasta_index: &FastaIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<(SpoaGraph, Vec<SequenceMetadata>)> {
    // Create a SPOA graph
    let mut graph = SpoaGraph::new();

    // Create scoring parameters for alignment
    let (match_score, mismatch, gap_open1, gap_extend1, gap_open2, gap_extend2) = scoring_params;

    // Create an alignment engine with affine gap penalties
    let mut engine = AlignmentEngine::new_convex(
        SpoaAlignmentType::kSW, // Local alignment (Smith-Waterman)
        match_score as i8,      // match score (positive)
        -(mismatch as i8),      // mismatch penalty (negative)
        -(gap_open1 as i8),     // gap open penalty (negative)
        -(gap_extend1 as i8),   // gap extend penalty (negative)
        -(gap_open2 as i8),     // gap open penalty (negative)
        -(gap_extend2 as i8),   // gap extend penalty (negative)
    );

    // Collect sequences and metadata for each interval
    let mut sequence_metadata = Vec::new();

    for interval in results.iter() {
        let seq_name = impg.seq_index.get_name(interval.metadata).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence name not found for ID {}", interval.metadata),
            )
        })?;

        // Get total sequence length
        let total_length = impg
            .seq_index
            .get_len_from_id(interval.metadata)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("Sequence length not found for ID {}", interval.metadata),
                )
            })?;

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
        let seq_size = end - start;

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

        output.push_str(&format!(
            "s {} {} {} {} {} {}\n",
            meta.name, meta.start, aligned_size, meta.strand, meta.total_length, trimmed_seq
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
fn _convert_and_write_gfa<R: AsRef<str>, W: Write>(raw_gfa: R, writer: &mut W) -> io::Result<()> {
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
            let old_to = fields[3];

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
            let walk_str = fields[6]; // e.g. “>s0>s2>s3...”

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

use std::collections::BTreeMap;

#[derive(Clone)]
struct GroupInfo {
    name: String,
    sequence_indices: Vec<usize>,
    total_length: usize,
}

pub fn compute_and_output_similarities(
    impg: &Impg,
    results: &[Interval<u32>],
    fasta_index: &FastaIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
    emit_distances: bool,
    emit_all_pairs: bool,
    delim: Option<char>,
    delim_pos: u16,
) -> io::Result<()> {
    // Generate POA graph and get sequences
    let (graph, sequence_metadata) = prepare_poa_graph_and_sequences(
        impg,
        results,
        fasta_index,
        scoring_params,
    )?;

    let msa = graph.generate_msa();

    // Convert MSA strings to Vec<char> once for efficient indexing
    let msa_chars: Vec<Vec<char>> = msa.iter().map(|s| s.chars().collect()).collect();

    // Create groups and calculate their lengths
    let groups = create_groups(&sequence_metadata, delim, delim_pos)?;

    // Print header
    println!("group.a\tgroup.b\tgroup.a.length\tgroup.b.length\tintersection\t{}", 
        if emit_distances {
            "jaccard.distance\tcosine.distance\tdice.distance\testimated.difference.rate"
        } else {
            "jaccard.similarity\tcosine.similarity\tdice.similarity\testimated.identity"
        }
    );

    // Compute only upper triangle (including diagonal) and output both (i,j) and (j,i)
    for i in 0..groups.len() {
        for j in i..groups.len() {
            let group_a = &groups[i];
            let group_b = &groups[j];
            
            let intersection = if delim.is_none() && group_a.sequence_indices.len() == 1 && group_b.sequence_indices.len() == 1 {
                // Fast path for individual sequences (no grouping)
                calculate_pairwise_intersection(&msa_chars[group_a.sequence_indices[0]], &msa_chars[group_b.sequence_indices[0]])
            } else {
                // General case for groups
                calculate_group_intersection(&msa_chars, &group_a.sequence_indices, &group_b.sequence_indices)
            };
            
            if intersection == 0 && !emit_all_pairs {
                // Skip both (i,j) and (j,i) if no intersection
                continue;
            }

            // Compute metrics once
            let union = (group_a.total_length + group_b.total_length).saturating_sub(intersection);
            let jaccard = if union > 0 { intersection as f64 / union as f64 } else { 0.0 };
            let cosine = if group_a.total_length > 0 && group_b.total_length > 0 { 
                (intersection as f64) / ((group_a.total_length as f64).sqrt() * (group_b.total_length as f64).sqrt())
            } else { 0.0 };
            let dice = if (group_a.total_length + group_b.total_length) > 0 { 
                2.0 * (intersection as f64) / (group_a.total_length + group_b.total_length) as f64 
            } else { 0.0 };
            let estimated_identity = if jaccard > 0.0 { 2.0 * jaccard / (1.0 + jaccard) } else { 0.0 };

            // Output (i,j)
            output_similarity_line(
                &group_a.name, &group_b.name, 
                group_a.total_length, group_b.total_length, 
                intersection, jaccard, cosine, dice, estimated_identity, 
                emit_distances
            );

            // Output (j,i) if different from (i,j)
            if i != j {
                output_similarity_line(
                    &group_b.name, &group_a.name, 
                    group_b.total_length, group_a.total_length, 
                    intersection, jaccard, cosine, dice, estimated_identity, 
                    emit_distances
                );
            }
        }
    }

    Ok(())
}

fn output_similarity_line(
    name_a: &str, name_b: &str,
    len_a: usize, len_b: usize,
    intersection: usize,
    jaccard: f64, cosine: f64, dice: f64, estimated_identity: f64,
    emit_distances: bool
) {
    print!("{}\t{}\t{}\t{}\t{}\t", name_a, name_b, len_a, len_b, intersection);

    if emit_distances {
        println!("{}\t{}\t{}\t{}",
            format_similarity_value(1.0 - jaccard),
            format_similarity_value(1.0 - cosine),
            format_similarity_value(1.0 - dice),
            format_similarity_value(1.0 - estimated_identity)
        );
    } else {
        println!("{}\t{}\t{}\t{}",
            format_similarity_value(jaccard),
            format_similarity_value(cosine),
            format_similarity_value(dice),
            format_similarity_value(estimated_identity)
        );
    }
}

fn create_groups(
    sequence_metadata: &[SequenceMetadata],
    delim: Option<char>,
    delim_pos: u16,
) -> io::Result<Vec<GroupInfo>> {
    if let Some(delim_char) = delim {
        // Group by delimiter
        let mut group_map: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        
        for (seq_idx, meta) in sequence_metadata.iter().enumerate() {
            let group_name = extract_group_name(&meta.name, delim_char, delim_pos)?;
            group_map.entry(group_name).or_default().push(seq_idx);
        }
        
        let mut groups = Vec::new();
        for (group_name, indices) in group_map {
            // Simply sum all metadata.size values in the group
            let total_length: usize = indices.iter()
                .map(|&idx| sequence_metadata[idx].size as usize)
                .sum();
            
            groups.push(GroupInfo {
                name: group_name,
                sequence_indices: indices,
                total_length,
            });
        }
        
        Ok(groups)
    } else {
        // Each sequence is its own group
        let groups = sequence_metadata
            .iter()
            .enumerate()
            .map(|(i, meta)| {
                let name = format!("{}:{}-{}", meta.name, meta.start, meta.start + meta.size);
                GroupInfo {
                    name,
                    sequence_indices: vec![i],
                    total_length: meta.size as usize,
                }
            })
            .collect();
        
        Ok(groups)
    }
}

fn extract_group_name(path_name: &str, delim: char, delim_pos: u16) -> io::Result<String> {
    let mut delimiter_positions = Vec::new();
    let mut start = 0;
    
    // Find all delimiter positions
    while let Some(pos) = path_name[start..].find(delim) {
        delimiter_positions.push(start + pos);
        start += pos + 1;
    }
    
    let target_index = (delim_pos as usize).saturating_sub(1);
    
    if target_index < delimiter_positions.len() {
        Ok(path_name[..delimiter_positions[target_index]].to_string())
    } else {
        // Not enough delimiters, use the whole name
        Ok(path_name.to_string())
    }
}

// Fast path for comparing two individual sequences
fn calculate_pairwise_intersection(seq_a: &[char], seq_b: &[char]) -> usize {
    seq_a.iter()
        .zip(seq_b.iter())
        .filter(|&(&char_a, &char_b)| char_a != '-' && char_b != '-' && char_a == char_b)
        .count()
}

// General case for comparing groups of sequences
fn calculate_group_intersection(msa_chars: &[Vec<char>], group_a_indices: &[usize], group_b_indices: &[usize]) -> usize {    
    let msa_len = msa_chars[0].len();
    let mut intersection = 0;
    
    for pos in 0..msa_len {
        // For each position, check all pairwise comparisons between groups and count matches
        let mut position_matches = 0;
        
        for &idx_a in group_a_indices {
            let char_a = msa_chars[idx_a][pos];
            if char_a != '-' {
                for &idx_b in group_b_indices {
                    let char_b = msa_chars[idx_b][pos];
                    if char_b != '-' && char_a == char_b {
                        position_matches += 1;
                    }
                }
            }
        }
        
        // For this position, add the minimum of:
        // - number of non-gap characters in group A
        // - number of non-gap characters in group B  
        // - number of actual matches found
        let group_a_count = group_a_indices.iter()
            .filter(|&&idx| msa_chars[idx][pos] != '-')
            .count();
        let group_b_count = group_b_indices.iter()
            .filter(|&&idx| msa_chars[idx][pos] != '-')
            .count();
        
        intersection += position_matches.min(group_a_count).min(group_b_count);
    }
    
    intersection
}

fn format_similarity_value(value: f64) -> String {
    let formatted = format!("{:.7}", value);
    // Trim trailing zeros
    let trimmed = formatted.trim_end_matches('0');
    if trimmed.ends_with('.') {
        // If we removed all decimal places, remove the decimal point too
        trimmed.trim_end_matches('.').to_string()
    } else {
        trimmed.to_string()
    }
}