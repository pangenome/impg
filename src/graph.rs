use crate::impg_index::ImpgIndex;
use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use coitrees::Interval;
use spoa_rs::{AlignmentEngine, AlignmentType as SpoaAlignmentType, Graph as SpoaGraph};
use std::io::{self, BufWriter, Write};

// Gfasort imports for graph sorting
use gfasort::gfa_parser::load_gfa;
use gfasort::ygs::{unchop_only, ygs_sort, YgsParams};

#[derive(Clone)]
pub struct SequenceMetadata {
    pub name: String,
    pub start: i32,
    pub size: i32,
    pub strand: char,
    pub total_length: usize,
}

impl SequenceMetadata {
    /// Canonical GFA path name in forward-strand coordinates: `name:fwd_start-fwd_end`.
    ///
    /// For `+` strand this is `name:start-(start+size)`.
    /// For `-` strand this converts from MAF-style RC-frame coordinates back to
    /// forward-strand: `name:(total-start-size)-(total-start)`.
    pub fn path_name(&self) -> String {
        let (fwd_start, fwd_end) = if self.strand == '+' {
            (self.start, self.start + self.size)
        } else {
            (
                (self.total_length as i32) - self.start - self.size,
                (self.total_length as i32) - self.start,
            )
        };
        format!("{}:{}-{}", self.name, fwd_start, fwd_end)
    }

    /// Name used in FASTA headers for alignment tools (sweepga/FastGA).
    /// Includes strand to distinguish orientations: `name:start-end(strand)`.
    /// Uses raw metadata coordinates (MAF-style for `-` strand).
    pub fn alignment_name(&self) -> String {
        format!(
            "{}:{}-{}({})",
            self.name,
            self.start,
            self.start + self.size,
            self.strand
        )
    }
}

pub fn generate_gfa_from_intervals(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<String> {
    let (graph, sequence_metadata) =
        prepare_poa_graph_and_sequences(impg, results, sequence_index, scoring_params).unwrap();
    spoa_graph_to_unchoped_gfa(graph, &sequence_metadata)
}

/// Convert an SPOA graph + metadata into an unchopped GFA string.
///
/// Pipeline: generate GFA → post-process strands → unchop.
/// Callers apply gfaffix + sort via normalize_and_sort.
pub(crate) fn spoa_graph_to_unchoped_gfa(
    graph: SpoaGraph,
    sequence_metadata: &[SequenceMetadata],
) -> io::Result<String> {
    let headers: Vec<String> = sequence_metadata
        .iter()
        .map(|meta| meta.path_name())
        .collect();
    let gfa_output = graph.generate_gfa(&headers, false);
    let gfa_output = post_process_gfa_for_strands(gfa_output, sequence_metadata);
    unchop_gfa(&gfa_output)
}

pub fn post_process_gfa_for_strands(gfa: String, sequence_metadata: &[SequenceMetadata]) -> String {
    let mut output = String::new();

    // Create a map of sequence names to their strand information
    let strand_map: std::collections::HashMap<String, char> = sequence_metadata
        .iter()
        .map(|meta| (meta.path_name(), meta.strand))
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
                                    format!("{seg_stripped}-")
                                } else if let Some(seg_stripped) = seg.strip_suffix('-') {
                                    format!("{seg_stripped}+")
                                } else {
                                    seg.to_string()
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
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> String {
    // Prepare POA graph and sequences
    let (graph, sequence_metadata) =
        prepare_poa_graph_and_sequences(impg, results, sequence_index, scoring_params).unwrap();

    // Generate MSA from the SPOA graph
    let msa = graph.generate_msa();

    // Print MAF format to stdout
    format_maf_from_msa(&msa, &sequence_metadata, None)
}

/// Find the first and last non-gap columns in an MSA, returning `(left, right)` where
/// `right` is exclusive. If the MSA is empty or all-gap, returns `(0, msa_len)`.
fn find_msa_trim_bounds(msa: &[String]) -> (usize, usize) {
    let msa_len = msa.first().map(|s| s.len()).unwrap_or(0);
    if msa_len == 0 {
        return (0, 0);
    }
    let mut left = 0usize;
    let mut right = msa_len;

    'find_left: for i in 0..msa_len {
        for s in msa {
            if s.as_bytes()[i] != b'-' {
                left = i;
                break 'find_left;
            }
        }
    }
    'find_right: for i in (0..msa_len).rev() {
        for s in msa {
            if s.as_bytes()[i] != b'-' {
                right = i + 1;
                break 'find_right;
            }
        }
    }
    if right < left {
        (0, msa_len)
    } else {
        (left, right)
    }
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
        output.push_str(&format!("# {name}\n"));
    }

    let (start_trim, end_trim) = find_msa_trim_bounds(msa);

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

pub fn generate_fasta_alignment_from_intervals(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> String {
    // Prepare POA graph and sequences (same as MAF)
    let (graph, sequence_metadata) =
        prepare_poa_graph_and_sequences(impg, results, sequence_index, scoring_params)
            .expect("failed to prepare SPOA graph/sequences");

    // Generate MSA from the SPOA graph
    let msa = graph.generate_msa();

    // Format as aligned FASTA
    // - trim all-gap columns on both ends (like MAF)
    // - wrap lines at 80 columns (standard FASTA convention)
    format_fasta_alignment_from_msa(&msa, &sequence_metadata, 80, true)
}

fn format_fasta_alignment_from_msa(
    msa: &[String],
    sequence_metadata: &[SequenceMetadata],
    line_width: usize,  // e.g., 80
    trim_all_gap: bool, // true to trim leading/trailing all-gap columns
) -> String {
    let mut out = String::new();

    if msa.is_empty() || sequence_metadata.is_empty() {
        return out;
    }

    // Determine trimming window (remove all-gap columns at the ends)
    let aln_len = msa.first().map(|s| s.len()).unwrap_or(0);
    let (left, right) = if trim_all_gap && aln_len > 0 {
        find_msa_trim_bounds(msa)
    } else {
        (0, aln_len)
    };

    // Emit aligned FASTA:
    // Header encodes original label + coordinates consistent with MAF logic:
    // meta.start is already "MAF-style" start (reverse-complemented for '-' strand),
    // so end = start + size is consistent for both '+' and '-'.
    for (seq, meta) in msa.iter().zip(sequence_metadata.iter()) {
        let start = meta.start;
        let end = meta.start + meta.size;

        // Example header: >name:start-end(strand)
        // You can adjust to omit strand or change style if you prefer.
        out.push('>');
        out.push_str(&meta.name);
        out.push(':');
        out.push_str(&start.to_string());
        out.push('-');
        out.push_str(&end.to_string());
        out.push('(');
        out.push(meta.strand);
        out.push(')');
        out.push('\n');

        let slice = &seq[left..right];
        if line_width == 0 {
            // no wrapping
            out.push_str(slice);
            out.push('\n');
        } else {
            // wrap at line_width (80 by default)
            let bytes = slice.as_bytes();
            let mut i = 0usize;
            while i < bytes.len() {
                let j = (i + line_width).min(bytes.len());
                // safe because MSA is ASCII-only (A,C,G,T,N,'-')
                out.push_str(&String::from_utf8_lossy(&bytes[i..j]));
                out.push('\n');
                i = j;
            }
        }
    }

    out
}

/// Build a new SPOA graph and alignment engine with the given scoring parameters.
pub(crate) fn build_spoa_engine(
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> (SpoaGraph, AlignmentEngine) {
    let (match_score, mismatch, gap_open1, gap_extend1, gap_open2, gap_extend2) = scoring_params;
    let graph = SpoaGraph::new();
    let engine = AlignmentEngine::new_convex(
        SpoaAlignmentType::kSW,
        match_score as i8,
        -(mismatch as i8),
        -(gap_open1 as i8),
        -(gap_extend1 as i8),
        -(gap_open2 as i8),
        -(gap_extend2 as i8),
    );
    (graph, engine)
}

/// Add sequences to a SPOA graph via `engine`. Skips empty sequences.
pub(crate) fn feed_sequences_to_graph(
    engine: &mut AlignmentEngine,
    graph: &mut SpoaGraph,
    sequences: impl Iterator<Item = impl AsRef<str>>,
) {
    for seq in sequences {
        let seq = seq.as_ref();
        if seq.is_empty() {
            continue;
        }
        let weights = vec![1u32; seq.len()];
        let (_, alignment) = engine.align(seq, graph);
        graph.add_alignment_with_weights(alignment, seq, &weights);
    }
}

pub fn prepare_poa_graph_and_sequences(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<(SpoaGraph, Vec<SequenceMetadata>)> {
    let mut processed_sequences = prepare_sequences(impg, results, sequence_index)?;
    // SPOA benefits from longest-first feeding order
    processed_sequences.sort_by(|a, b| b.0.len().cmp(&a.0.len()));

    let (mut graph, mut engine) = build_spoa_engine(scoring_params);

    // Collect metadata and feed sequences to SPOA graph (SPOA is not thread-safe)
    let sequence_metadata: Vec<SequenceMetadata> =
        processed_sequences.iter().map(|(_, m)| m.clone()).collect();
    feed_sequences_to_graph(
        &mut engine,
        &mut graph,
        processed_sequences.iter().map(|(s, _)| s.as_str()),
    );

    Ok((graph, sequence_metadata))
}

pub fn prepare_sequences(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
) -> std::io::Result<Vec<(String, SequenceMetadata)>> {
    use rayon::prelude::*;

    // Fetch, strand-normalize, and annotate each interval
    let pairs: Vec<(String, SequenceMetadata)> = results
        .par_iter()
        .map(|interval| -> std::io::Result<(String, SequenceMetadata)> {
            // Resolve sequence name
            let seq_name = impg
                .seq_index()
                .get_name(interval.metadata)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("Sequence name not found for ID {}", interval.metadata),
                    )
                })?;

            // Resolve total contig length
            let total_length = impg
                .seq_index()
                .get_len_from_id(interval.metadata)
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::NotFound,
                        format!("Sequence length not found for ID {}", interval.metadata),
                    )
                })?;

            // Determine orientation from the interval
            let (start, end, strand) = if interval.first <= interval.last {
                (interval.first, interval.last, '+')
            } else {
                (interval.last, interval.first, '-')
            };

            // Fetch subsequence on the forward strand coordinates [start, end)
            let mut seq_bytes = sequence_index.fetch_sequence(seq_name, start, end)?;

            // Reverse-complement if the interval is on the reverse strand
            if strand == '-' {
                seq_bytes = reverse_complement(&seq_bytes);
            }

            // Convert to UTF-8 String for downstream code that expects Strings
            let sequence_str = String::from_utf8_lossy(&seq_bytes).to_string();

            // Size of the fetched slice in original coordinates
            let seq_size = end - start;

            // MAF-style start: if reverse, origin is from the RC frame
            let maf_start = if strand == '-' {
                (total_length as i32) - end
            } else {
                start
            };

            let meta = SequenceMetadata {
                name: seq_name.to_string(),
                start: maf_start,
                size: seq_size,
                strand,
                total_length,
            };

            Ok((sequence_str, meta))
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(pairs)
}

/// Extract a multiple sequence alignment (MSA) from a GFA graph string.
///
/// Parses segment (S) and path (P) lines from the GFA. Uses the sorted segment
/// order as the column order. For each path, outputs the segment sequence at
/// positions the path visits, or a gap string of equal length where it doesn't.
///
/// Returns a vector of (path_name, msa_string) pairs, in the order paths appear
/// in the GFA.
pub fn gfa_to_msa(gfa: &str) -> Vec<(String, String)> {
    use std::collections::HashMap;

    // Parse segments: id -> sequence
    let mut segments: HashMap<String, String> = HashMap::new();
    // Ordered list of segment IDs (preserving GFA order for column layout)
    let mut segment_order: Vec<String> = Vec::new();
    // Parse paths: name -> vec of (segment_id, orientation)
    let mut paths: Vec<(String, Vec<(String, char)>)> = Vec::new();

    for line in gfa.lines() {
        if line.starts_with("S\t") {
            let fields: Vec<&str> = line.splitn(4, '\t').collect();
            if fields.len() >= 3 {
                let id = fields[1].to_string();
                let seq = fields[2].to_string();
                if !segments.contains_key(&id) {
                    segment_order.push(id.clone());
                }
                segments.insert(id, seq);
            }
        } else if line.starts_with("P\t") {
            let fields: Vec<&str> = line.splitn(4, '\t').collect();
            if fields.len() >= 3 {
                let path_name = fields[1].to_string();
                let steps: Vec<(String, char)> = fields[2]
                    .split(',')
                    .map(|step| {
                        if let Some(stripped) = step.strip_suffix('+') {
                            (stripped.to_string(), '+')
                        } else if let Some(stripped) = step.strip_suffix('-') {
                            (stripped.to_string(), '-')
                        } else {
                            (step.to_string(), '+')
                        }
                    })
                    .collect();
                paths.push((path_name, steps));
            }
        }
    }

    if segments.is_empty() || paths.is_empty() {
        return Vec::new();
    }

    // For each path, build a set of visited segments
    // Then construct the MSA row
    let mut result = Vec::with_capacity(paths.len());

    for (path_name, steps) in &paths {
        // Track which segments this path visits and in what order
        let mut visited: HashMap<&str, char> = HashMap::new();
        for (seg_id, orient) in steps {
            visited.insert(seg_id.as_str(), *orient);
        }

        let mut msa_row = String::new();
        for seg_id in &segment_order {
            let seq = &segments[seg_id];
            if let Some(&orient) = visited.get(seg_id.as_str()) {
                if orient == '-' {
                    // Reverse complement
                    let rc: String = seq
                        .bytes()
                        .rev()
                        .map(|b| match b {
                            b'A' | b'a' => 'T',
                            b'T' | b't' => 'A',
                            b'C' | b'c' => 'G',
                            b'G' | b'g' => 'C',
                            _ => 'N',
                        })
                        .collect();
                    msa_row.push_str(&rc);
                } else {
                    msa_row.push_str(seq);
                }
            } else {
                // Gap: same length as segment
                for _ in 0..seq.len() {
                    msa_row.push('-');
                }
            }
        }

        result.push((path_name.clone(), msa_row));
    }

    result
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
                writeln!(writer, "{line}")?;
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
                writeln!(writer, "{line}")?;
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
                writeln!(writer, "{line}")?;
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
                    segs_oriented.push(format!("{new_seg}{orient}"));
                }
            }

            let joined = segs_oriented.join(",");
            // Emit as a `P` line, with Overlaps = “*”
            writeln!(writer, "P\t{path_name}\t{joined}\t*")?;
            continue;
        }

        // 5) All other lines (including C/J/P/L/W tags not shown above, or comment lines): pass through:
        writeln!(writer, "{line}")?;
    }

    Ok(())
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let len = seq.len();
    let mut result = Vec::with_capacity(len);
    for &base in seq.iter().rev() {
        result.push(match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => base,
        });
    }
    result
}

/// Sort a GFA string using gfasort's Ygs pipeline (path-guided SGD + grooming + topological sort)
///
/// This is equivalent to running `gfasort -p Ygs` on the command line.
/// It produces a well-ordered graph with nodes arranged to minimize path distances.
/// Compact consecutive nodes (unchop) without sorting.
/// Reduces 1-bp SPOA nodes into longer segments early, shrinking the graph for downstream steps.
pub(crate) fn unchop_gfa(gfa_content: &str) -> io::Result<String> {
    let temp_gfa = tempfile::Builder::new()
        .suffix(".gfa")
        .tempfile()
        .map_err(|e| io::Error::other(format!("unchop: failed to create temp file: {}", e)))?;
    std::fs::write(temp_gfa.path(), gfa_content)
        .map_err(|e| io::Error::other(format!("unchop: failed to write temp file: {}", e)))?;

    let mut graph = load_gfa(temp_gfa.path())
        .map_err(|e| io::Error::other(format!("unchop: failed to load GFA: {}", e)))?;

    if graph.nodes.iter().filter(|n| n.is_some()).count() <= 1 {
        return Ok(gfa_content.to_string());
    }

    unchop_only(&mut graph, 0);

    let mut output = Vec::new();
    graph
        .write_gfa(&mut output)
        .map_err(|e| io::Error::other(format!("unchop: failed to write GFA: {}", e)))?;
    String::from_utf8(output).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unchop: invalid UTF-8 in output: {}", e),
        )
    })
}

/// Apply gfaffix normalization then SGD sort.
///
/// Shared final step for all engines. Fails with an error if gfaffix is not
/// found or exits non-zero — callers must ensure the binary is available.
pub fn normalize_and_sort(gfa: String, num_threads: usize) -> io::Result<String> {
    let normalized = run_gfaffix(&gfa, num_threads)?;
    sort_gfa(&normalized, num_threads)
}

pub fn sort_gfa(gfa_content: &str, num_threads: usize) -> io::Result<String> {
    // Write GFA to temp file (gfasort's load_gfa requires a file path)
    let temp_gfa = tempfile::Builder::new()
        .suffix(".gfa")
        .tempfile()
        .map_err(|e| io::Error::other(format!("Failed to create temp GFA file: {}", e)))?;

    std::fs::write(temp_gfa.path(), gfa_content)
        .map_err(|e| io::Error::other(format!("Failed to write temp GFA: {}", e)))?;

    // Load GFA into gfasort's graph structure
    let mut graph = load_gfa(temp_gfa.path())
        .map_err(|e| io::Error::other(format!("Failed to load GFA for sorting: {}", e)))?;

    // Skip sorting for trivial graphs (0 or 1 node)
    let node_count = graph.nodes.iter().filter(|n| n.is_some()).count();
    if node_count <= 1 {
        return Ok(gfa_content.to_string());
    }

    // Create Ygs parameters from the graph
    let params = YgsParams::from_graph(&graph, 0, num_threads);

    // Sort the graph using Ygs pipeline (path-guided SGD + grooming + topological sort)
    ygs_sort(&mut graph, &params);

    // Write sorted GFA to string
    let mut sorted_output = Vec::new();
    graph
        .write_gfa(&mut sorted_output)
        .map_err(|e| io::Error::other(format!("Failed to write sorted GFA: {}", e)))?;

    String::from_utf8(sorted_output).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid UTF-8 in sorted GFA: {}", e),
        )
    })
}

/// Run gfaffix graph normalization on a GFA string.
///
/// Searches for the `gfaffix` binary in PATH first, then next to the current executable
/// (so it works when both are built together in `target/release/`).
/// Returns an error if the binary is not found, allowing callers to fall back gracefully.
pub fn run_gfaffix(gfa_content: &str, _num_threads: usize) -> io::Result<String> {
    use std::process::Command;

    // Resolve binary: sibling of current exe first (e.g. target/release/gfaffix),
    // then fall back to PATH lookup.
    let gfaffix_bin = std::env::current_exe()
        .ok()
        .map(|exe| exe.with_file_name("gfaffix"))
        .filter(|p| p.exists())
        .or_else(|| {
            std::env::var_os("PATH").and_then(|paths| {
                std::env::split_paths(&paths)
                    .map(|dir| dir.join("gfaffix"))
                    .find(|p| p.exists())
            })
        })
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                "gfaffix binary not found in PATH or next to impg executable",
            )
        })?;

    // Write input GFA to temp file
    let input_file = tempfile::Builder::new()
        .suffix(".gfa")
        .tempfile()
        .map_err(|e| io::Error::other(format!("gfaffix: failed to create temp input: {}", e)))?;
    std::fs::write(input_file.path(), gfa_content.as_bytes())
        .map_err(|e| io::Error::other(format!("gfaffix: failed to write temp input: {}", e)))?;

    // Output GFA temp file
    let output_file = tempfile::Builder::new()
        .suffix(".gfa")
        .tempfile()
        .map_err(|e| io::Error::other(format!("gfaffix: failed to create temp output: {}", e)))?;

    // Run gfaffix
    let status = Command::new(&gfaffix_bin)
        .arg(input_file.path())
        .arg("-o")
        .arg(output_file.path())
        .status()
        .map_err(|e| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("gfaffix failed to start: {}", e),
            )
        })?;

    if !status.success() {
        return Err(io::Error::other(format!(
            "gfaffix exited with non-zero status: {:?}",
            status.code()
        )));
    }

    std::fs::read_to_string(output_file.path())
        .map_err(|e| io::Error::other(format!("gfaffix: failed to read output: {}", e)))
}

/// Configuration for seqwish-based GFA generation
pub struct SeqwishConfig {
    /// Number of threads for parallel processing
    pub num_threads: usize,
    /// K-mer frequency multiplier (frequency = num_sequences * multiplier)
    pub frequency_multiplier: usize,
    /// Minimum alignment length for the aligner
    pub min_aln_length: u64,
    /// Optional temp directory for intermediate files
    pub temp_dir: Option<String>,
    /// Skip PAF filtering (faster but may produce broken graphs)
    pub no_filter: bool,
    /// Aligner backend: "wfmash" or "fastga"
    pub aligner: String,
    /// n:m-best mappings kept in query:target dimensions (e.g., "1:1", "many:many")
    pub num_mappings: String,
    /// Scaffold jump/gap distance (0 = disable scaffolding)
    pub scaffold_jump: u64,
    /// Minimum scaffold chain length
    pub scaffold_mass: u64,
    /// Scaffold filter mode (e.g., "1:1", "many:many")
    pub scaffold_filter: String,
    /// Maximum overlap ratio for plane sweep filtering
    pub overlap: f64,
    /// Minimum identity threshold (0.0-1.0)
    pub min_identity: f64,
    /// Maximum scaffold deviation distance (0 = no limit)
    pub scaffold_dist: u64,
    /// Minimum mapping length for post-alignment filtering
    pub min_map_length: u64,
    /// Optional directory to save intermediate debug files (FASTA, raw PAF, filtered PAF).
    pub debug_dir: Option<String>,
    /// Unified sparsification strategy.
    pub sparsify: sweepga::knn_graph::SparsificationStrategy,
    /// Mash distance parameters for sparsification sketching.
    pub mash_params: sweepga::knn_graph::MashParams,
    /// Maximum repeat count for transitive closure (0 = no limit)
    pub repeat_max: u64,
    /// Minimum distance between repeats
    pub min_repeat_dist: u64,
    /// Minimum match length filter for alignments
    pub min_match_len: u64,
    /// Sparse factor for input matches (0.0 = keep all)
    pub sparse_factor: f32,
    /// Batch size for transitive closure computation
    pub transclose_batch: u64,
    /// Use in-memory interval trees (false = disk-backed, slower but lower memory)
    pub use_in_memory: bool,
}

impl Default for SeqwishConfig {
    fn default() -> Self {
        SeqwishConfig {
            num_threads: 4,
            frequency_multiplier: 10,
            min_aln_length: 0,
            temp_dir: None,
            no_filter: false,
            aligner: "wfmash".to_string(),
            num_mappings: "many:many".to_string(),
            scaffold_jump: 50_000,
            scaffold_mass: 10_000,
            scaffold_filter: "many:many".to_string(),
            overlap: 0.95,
            min_identity: 0.0,
            scaffold_dist: 0,
            min_map_length: 0,
            debug_dir: None,
            sparsify: sweepga::knn_graph::SparsificationStrategy::None,
            mash_params: sweepga::knn_graph::MashParams::default(),
            repeat_max: 0,
            min_repeat_dist: 0,
            min_match_len: 23,
            sparse_factor: 0.0,
            transclose_batch: 10_000_000,
            use_in_memory: true,
        }
    }
}

/// Generate a variation graph from intervals using sweepga+seqwish (graph induction
/// via transitive closure). The resulting graph is raw (unsmoothed).
///
/// Delegates to `commands::graph::build_graph()` so that query and graph
/// subcommands use exactly the same alignment, filtering, and graph-building
/// pipeline.
pub fn generate_gfa_seqwish_from_intervals(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    config: &SeqwishConfig,
) -> io::Result<String> {
    if results.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }

    // 1) Extract sequences from intervals
    let sequences = prepare_sequences(impg, results, sequence_index)?;

    if sequences.is_empty() {
        return Ok(String::from("H\tVN:Z:1.0\n"));
    }

    // 2) Write sequences to a temp FASTA file
    let combined_fasta = tempfile::Builder::new()
        .suffix(".fa")
        .tempfile()
        .map_err(|e| io::Error::other(format!("Failed to create temp file: {}", e)))?;

    {
        let mut writer = BufWriter::new(&combined_fasta);
        for (seq, meta) in &sequences {
            let header = format!(">{}", meta.path_name());
            writeln!(writer, "{}", header)?;
            writeln!(writer, "{}", seq)?;
        }
        writer.flush()?;
    }

    // 3) Build graph using the same pipeline as `graph --engine seqwish`
    let graph_config = crate::commands::graph::GraphBuildConfig {
        num_threads: config.num_threads,
        frequency_multiplier: config.frequency_multiplier,
        min_aln_length: config.min_aln_length,
        temp_dir: config.temp_dir.clone(),
        no_filter: config.no_filter,
        aligner: config.aligner.clone(),
        num_mappings: config.num_mappings.clone(),
        scaffold_jump: config.scaffold_jump,
        scaffold_mass: config.scaffold_mass,
        scaffold_filter: config.scaffold_filter.clone(),
        overlap: config.overlap,
        min_identity: config.min_identity,
        scaffold_dist: config.scaffold_dist,
        min_map_length: config.min_map_length,
        debug_dir: config.debug_dir.clone(),
        sparsify: config.sparsify.clone(),
        mash_params: config.mash_params.clone(),
        repeat_max: config.repeat_max,
        min_repeat_dist: config.min_repeat_dist,
        min_match_len: config.min_match_len,
        sparse_factor: config.sparse_factor,
        transclose_batch: config.transclose_batch,
        use_in_memory: config.use_in_memory,
        show_progress: false,
        ..crate::commands::graph::GraphBuildConfig::default()
    };

    let fasta_path = combined_fasta
        .path()
        .to_str()
        .ok_or_else(|| io::Error::other("Temp FASTA path is not valid UTF-8"))?
        .to_string();

    let mut gfa_output = Vec::new();
    crate::commands::graph::build_graph(&[fasta_path], &mut gfa_output, &graph_config)?;

    String::from_utf8(gfa_output).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid UTF-8 in GFA: {}", e),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gfa_to_msa_basic() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
S\t2\tTG
S\t3\tAATT
P\tseq1:0-10\t1+,2+,3+\t*
P\tseq2:0-8\t1+,3+\t*
";
        let msa = gfa_to_msa(gfa);
        assert_eq!(msa.len(), 2);
        assert_eq!(msa[0].0, "seq1:0-10");
        assert_eq!(msa[0].1, "ACGTTGAATT");
        assert_eq!(msa[1].0, "seq2:0-8");
        assert_eq!(msa[1].1, "ACGT--AATT");
    }

    #[test]
    fn test_gfa_to_msa_empty() {
        let gfa = "H\tVN:Z:1.0\n";
        let msa = gfa_to_msa(gfa);
        assert!(msa.is_empty());
    }

    #[test]
    fn test_gfa_to_msa_single_path() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
P\tseq1:0-4\t1+\t*
";
        let msa = gfa_to_msa(gfa);
        assert_eq!(msa.len(), 1);
        assert_eq!(msa[0].1, "ACGT");
    }

    #[test]
    fn test_gfa_to_msa_reverse_orient() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
P\tseq1:0-4\t1-\t*
";
        let msa = gfa_to_msa(gfa);
        assert_eq!(msa.len(), 1);
        // Reverse complement of ACGT is ACGT
        assert_eq!(msa[0].1, "ACGT");
    }

    #[test]
    fn test_gfa_to_msa_reverse_orient_asymmetric() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tACG
P\tseq1:0-3\t1+\t*
P\tseq2:0-3\t1-\t*
";
        let msa = gfa_to_msa(gfa);
        assert_eq!(msa.len(), 2);
        assert_eq!(msa[0].1, "ACG"); // forward
        assert_eq!(msa[1].1, "CGT"); // reverse complement of ACG
    }
}
