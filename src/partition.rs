use crate::impg::CigarOp;
use crate::impg::Impg;
use crate::impg::SortedRanges;
use coitrees::Interval;
use log::{debug, info};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
//use std::time::Instant;

pub fn partition_alignments(
    impg: &Impg,
    window_size: usize,
    sequence_prefix: &str,
    merge_distance: i32,
    min_missing_size: i32,
    min_boundary_distance: i32,
    max_depth: u16,
    min_transitive_len: i32,
    min_distance_between_ranges: i32,
    selection_mode: Option<&str>,
    debug: bool,
) -> io::Result<()> {
    // Get all sequences with the given prefix
    let mut sample_regions = Vec::<(u32, i32, i32)>::new();
    let mut total_sequence_length = 0;
    for seq_id in 0..impg.seq_index.len() as u32 {
        let seq_name = impg.seq_index.get_name(seq_id).unwrap();
        let seq_length = impg.seq_index.get_len_from_id(seq_id).unwrap();
        if seq_name.starts_with(sequence_prefix) {
            sample_regions.push((seq_id, 0, seq_length as i32));
        }
        total_sequence_length += seq_length;
    }
    if sample_regions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "No sequences with prefix {} found in sequence index",
                sequence_prefix
            ),
        ));
    }
    // Natural sort by sequence name
    sample_regions.sort_by(|a, b| {
        let chrom_a = impg.seq_index.get_name(a.0).unwrap();
        let chrom_b = impg.seq_index.get_name(b.0).unwrap();
        natord::compare(chrom_a, chrom_b)
    });
    info!("Total sequence length: {} bp", total_sequence_length);

    if debug {
        debug!(
            "Found {} sequences with prefix {}",
            sample_regions.len(),
            sequence_prefix
        );
        for (seq_id, start, end) in &sample_regions {
            let chrom = impg.seq_index.get_name(*seq_id).unwrap();
            debug!(
                "  Sequence: {}:{}-{}, len: {}",
                chrom,
                start,
                end,
                end - start
            );
        }
    }

    // Create windows from sample regions
    let mut windows = Vec::<(u32, i32, i32)>::new();
    for (seq_id, start, end) in sample_regions {
        let mut pos = start;
        while pos < end as i32 {
            let window_end = std::cmp::min(pos + window_size as i32, end as i32);

            // If this tail‐window is too small, merge it into the last one
            if window_end - pos < window_size as i32
                && !windows.is_empty()
                && windows.last().unwrap().0 == seq_id
            {
                let last = windows.last_mut().unwrap();
                last.2 = end as i32;
                break;
            }

            windows.push((seq_id, pos, window_end));
            pos = window_end;
        }
    }

    if debug {
        debug!("Starting with {} windows:", windows.len());
        for (seq_id, start, end) in &windows {
            let chrom = impg.seq_index.get_name(*seq_id).unwrap();
            debug!(
                "  Window: {}:{}-{}, len: {}",
                chrom,
                start,
                end,
                end - start
            );
        }
    }

    // Initialize masked regions
    let mut masked_regions: FxHashMap<u32, SortedRanges> = (0..impg.seq_index.len() as u32)
        .map(|id| {
            let len = impg.seq_index.get_len_from_id(id).unwrap();
            (id, SortedRanges::new(len as i32, 0))
        })
        .collect();

    // Initialize missing regions from sequence index
    let mut missing_regions: FxHashMap<u32, SortedRanges> = (0..impg.seq_index.len() as u32)
        .map(|id| {
            let len = impg.seq_index.get_len_from_id(id).unwrap();
            let mut ranges = SortedRanges::new(len as i32, 0);
            ranges.insert((0, len as i32));
            (id, ranges)
        })
        .collect();

    let mut partition_num = 0;
    let mut total_partitioned_length = 0;

    info!("Partitioning");

    while !windows.is_empty() {
        debug!("Processing new window set");

        for (seq_id, start, end) in windows.iter() {
            let chrom = impg.seq_index.get_name(*seq_id).unwrap();

            if debug {
                debug!(
                    "  Missing {} regions in {} sequences",
                    missing_regions
                        .values()
                        .map(|ranges| ranges.len())
                        .sum::<usize>(),
                    missing_regions.len()
                );
                for (chrom, ranges) in &missing_regions {
                    let chrom = impg.seq_index.get_name(*chrom).unwrap();
                    for &(start, end) in ranges.iter() {
                        debug!(
                            "    Region: {}:{}-{}, len: {}",
                            chrom,
                            start,
                            end,
                            end - start
                        );
                    }
                }

                debug!(
                    "  Masked {} regions in {} sequences",
                    masked_regions
                        .values()
                        .map(|ranges| ranges.len())
                        .sum::<usize>(),
                    masked_regions.len()
                );
                for (chrom, ranges) in &masked_regions {
                    let chrom = impg.seq_index.get_name(*chrom).unwrap();
                    for &(start, end) in ranges.iter() {
                        debug!(
                            "    Region: {}:{}-{}, len: {}",
                            chrom,
                            start,
                            end,
                            end - start
                        );
                    }
                }
            }

            // Query overlaps for current window
            //let query_start = Instant::now();
            let mut overlaps = impg.query_transitive(
                *seq_id,
                *start,
                *end,
                Some(&masked_regions),
                max_depth,
                min_transitive_len,
                min_distance_between_ranges,
                false, // Don't store CIGAR strings during partitioning
            );
            //let query_time = query_start.elapsed();
            debug!("  Collected {} query overlaps", overlaps.len());

            // Ignore CIGAR strings and target intervals.
            //debug!("  Merging overlaps closer than {}bp", merge_distance); // bedtools sort | bedtools merge -d merge_distance
            //let merge_start = Instant::now();
            merge_overlaps(&mut overlaps, merge_distance);
            //let merge_time = merge_start.elapsed();
            debug!(
                "  Collected {} query overlaps after merging those closer than {}bp",
                overlaps.len(),
                merge_distance
            );

            if min_boundary_distance > 0 {
                //debug!("  Extending short intervals");
                //let extend_start = Instant::now();
                extend_to_close_boundaries(&mut overlaps, impg, min_boundary_distance);
                //let extend_time = extend_start.elapsed();
                debug!(
                    "  Collected {} query overlaps after extending those close to boundaries",
                    overlaps.len()
                );
            }

            //debug!("  Excluding masked regions"); // bedtools subtract -a "partition$num.tmp.bed" -b "$MASK_BED"
            //let mask_start = Instant::now();
            overlaps = mask_and_update_regions(
                &mut overlaps,
                &mut masked_regions,
                &mut missing_regions,
                min_missing_size,
            );
            //let mask_time = mask_start.elapsed();

            if !overlaps.is_empty() {
                debug!(
                    "  Collected {} query overlaps after masking",
                    overlaps.len()
                );
                merge_overlaps(&mut overlaps, 0); // Final merge to ensure no overlaps remain
                debug!(
                    "  Collected {} query overlaps after re-merging",
                    overlaps.len()
                );

                // Calculate current partition length
                let current_partition_length: u64 = overlaps
                    .iter()
                    .map(|(interval, _, _)| (interval.last - interval.first).abs() as u64)
                    .sum();
                total_partitioned_length += current_partition_length;

                // Calculate percentages
                let current_percentage =
                    (current_partition_length as f64 / total_sequence_length as f64) * 100.0;
                let total_percentage =
                    (total_partitioned_length as f64 / total_sequence_length as f64) * 100.0;

                // Create formatted percentage strings with conditional scientific notation
                let current_percentage_str = if current_percentage < 0.0001 {
                    format!("{:.4e}%", current_percentage)
                } else {
                    format!("{:.4}%", current_percentage)
                };
                let total_percentage_str = if total_percentage < 0.0001 {
                    format!("{:.4e}%", total_percentage)
                } else {
                    format!("{:.4}%", total_percentage)
                };
                info!("  Writing partition {} with {} regions (query {}:{}-{}, len: {}) - {} of total sequence ({} so far)", 
                    partition_num,
                    overlaps.len(),
                    chrom,
                    start,
                    end,
                    end - start,
                    current_percentage_str,
                    total_percentage_str
                );

                //let write_start = Instant::now();
                write_partition(partition_num, &overlaps, impg)?;
                //let write_time = write_start.elapsed();

                partition_num += 1;

                // info!("  Partition {} timings: query={:?}, merge={:?}, mask={:?}, update={:?}, extend={:?}, write={:?}",
                //     partition_num, query_time, merge_time, mask_time, update_time, extend_time, write_time);
            } else {
                debug!(
                    "  No overlaps found for region {}:{}-{}, len: {}",
                    chrom,
                    start,
                    end,
                    end - start
                );
            }
        }

        // Clear existing windows but keep the allocation
        windows.clear();

        // Vector to collect ranges for windowing
        let mut ranges_to_window: Vec<(u32, i32, i32)> = Vec::new();
                
        // Determine which selection mode to use
        let use_longest_missing = selection_mode.map_or(false, |mode| mode == "none");
        let pan_sn_config = selection_mode.filter(|&mode| mode != "none");

        // Find regions to process based on the selection strategy
        if use_longest_missing {
            // Select longest single missing region
            let mut max_length = 0;
            let mut longest_region: Option<(u32, i32, i32)> = None;
            
            // Scan through missing regions
            for (seq_id, ranges) in missing_regions.iter() {
                for &(start, end) in ranges.iter() {
                    let length = end - start;
                    if length > max_length
                        || (length == max_length
                            && longest_region.map_or(true, |(curr_seq_id, curr_start, _)| {
                                (*seq_id < curr_seq_id)
                                    || (*seq_id == curr_seq_id && start < curr_start)
                            }))
                    {
                        max_length = length;
                        longest_region = Some((*seq_id, start, end));
                    }
                }
            }
            
            // Add the single longest region to windowing list
            if let Some(region) = longest_region {
                ranges_to_window.push(region);
            }
        } else if let Some(config) = pan_sn_config {
            // PanSN-based selection: find sample/haplotype with most missing sequence
            
            // Parse config: field[,separator]
            let mut parts = config.split(',');
            let field_type = parts.next().unwrap_or("sample");
            let separator = parts.next().unwrap_or("#");
            
            // Determine field count based on field_type
            let field_count = match field_type {
                "sample" => 1,
                "haplotype" => 2,
                _ => 1, // Default to sample if unrecognized
            };
            
            // Group sequences by sample/haplotype prefix
            let mut prefix_to_seqs: FxHashMap<String, Vec<u32>> = FxHashMap::default();
            
            for seq_id in 0..impg.seq_index.len() as u32 {
                if let Some(seq_name) = impg.seq_index.get_name(seq_id) {
                    let parts: Vec<&str> = seq_name.split(separator).collect();
                    if parts.len() >= field_count {
                        let prefix = parts[..field_count].join(separator);
                        prefix_to_seqs.entry(prefix).or_default().push(seq_id);
                    }
                }
            }
            
            // Calculate missing sequence per sample/haplotype
            let mut max_missing = 0;
            let mut prefix_with_most_missing: Option<String> = None;
            
            for (prefix, seqs) in &prefix_to_seqs {
                let mut total_missing = 0;
                
                for &seq_id in seqs {
                    if let Some(ranges) = missing_regions.get(&seq_id) {
                        total_missing += ranges.iter().map(|&(start, end)| end - start).sum::<i32>();
                    }
                }
                
                if total_missing > max_missing {
                    max_missing = total_missing;
                    prefix_with_most_missing = Some(prefix.clone());
                }
            }
            
            // If we found a prefix, add all its sequences to the windowing list
            if let Some(prefix) = prefix_with_most_missing {
                if let Some(seqs) = prefix_to_seqs.get(&prefix) {
                    // Get all sequences for this prefix and sort by length
                    let mut seqs_with_length: Vec<(u32, usize)> = seqs
                        .iter()
                        .filter_map(|&seq_id| {
                            impg.seq_index.get_len_from_id(seq_id).map(|len| (seq_id, len))
                        })
                        .collect();
                    
                    // Sort by length (descending)
                    seqs_with_length.sort_by(|a, b| b.1.cmp(&a.1));
                    
                    debug!("Selected {} sequences from {} with total missing {}bp",
                        seqs_with_length.len(), prefix, max_missing);
                    
                    // Add full range for each sequence
                    for (seq_id, _) in seqs_with_length {
                        let seq_length = impg.seq_index.get_len_from_id(seq_id).unwrap() as i32;
                        ranges_to_window.push((seq_id, 0, seq_length));
                    }
                }
            }
        } else {
            // Default strategy: select sequence with highest total missing
            let mut max_total_missing = 0;
            let mut seq_with_most_missing: Option<u32> = None;
            
            // Calculate total missing sequence for each seq_id
            for (seq_id, ranges) in missing_regions.iter() {
                let total_missing: i32 = ranges.iter().map(|&(start, end)| end - start).sum();
                
                if total_missing > max_total_missing {
                    max_total_missing = total_missing;
                    seq_with_most_missing = Some(*seq_id);
                }
            }
            
            // Add the sequence with most missing to windowing list
            if let Some(seq_id) = seq_with_most_missing {
                let seq_length = impg.seq_index.get_len_from_id(seq_id).unwrap() as i32;
                
                debug!("Selected sequence with ID {} for having most missing sequence ({}bp)",
                    seq_id, max_total_missing);
                
                ranges_to_window.push((seq_id, 0, seq_length));
            }
        }

        // Create windows from all collected ranges
        for (seq_id, start, end) in ranges_to_window {
            let mut pos = start;
            while pos < end {
                let window_end = std::cmp::min(pos + window_size as i32, end);
                
                // If this tail-window is too small, merge it into the last one
                if window_end - pos < window_size as i32
                    && !windows.is_empty()
                    && windows.last().unwrap().0 == seq_id
                {
                    let last = windows.last_mut().unwrap();
                    last.2 = end;
                    break;
                }
                
                windows.push((seq_id, pos, window_end));
                pos = window_end;
            }
        }

    }

    info!("Partitioned into {} regions", partition_num);

    Ok(())
}

fn merge_overlaps(overlaps: &mut Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)>, max_gap: i32) {
    if overlaps.len() > 1 && max_gap >= 0 {
        // Sort by sequence ID and start position
        overlaps.sort_unstable_by_key(|(query_interval, _, _)| {
            (
                query_interval.metadata,
                std::cmp::min(query_interval.first, query_interval.last),
            )
        });

        let mut write_idx = 0;
        for read_idx in 1..overlaps.len() {
            let (curr_interval, _, _) = &overlaps[write_idx];
            let (next_interval, _, _) = &overlaps[read_idx];

            let curr_min = std::cmp::min(curr_interval.first, curr_interval.last);
            let curr_max = std::cmp::max(curr_interval.first, curr_interval.last);
            let next_min = std::cmp::min(next_interval.first, next_interval.last);
            let next_max = std::cmp::max(next_interval.first, next_interval.last);

            if curr_interval.metadata != next_interval.metadata || next_min > curr_max + max_gap {
                write_idx += 1;
                if write_idx != read_idx {
                    overlaps.swap(write_idx, read_idx);
                }
            } else {
                overlaps[write_idx].0.first = curr_min.min(next_min);
                overlaps[write_idx].0.last = curr_max.max(next_max);
            }
        }
        overlaps.truncate(write_idx + 1);
    }
}

fn mask_and_update_regions(
    overlaps: &mut Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)>,
    masked_regions: &mut FxHashMap<u32, SortedRanges>,
    missing_regions: &mut FxHashMap<u32, SortedRanges>,
    min_fragment_size: i32,
) -> Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)> {
    let mut result = Vec::new();

    // First, identify and store extensions separately, indexed by sequence ID
    let mut extensions_by_seq_id: FxHashMap<u32, Vec<(i32, i32)>> = FxHashMap::default();

    // Collect extension info before draining overlaps
    for (query_interval, _, _) in overlaps.iter() {
        let seq_id = query_interval.metadata;
        let (mask_start, mask_end) = if query_interval.first <= query_interval.last {
            (query_interval.first, query_interval.last)
        } else {
            (query_interval.last, query_interval.first)
        };

        if let Some(missing) = missing_regions.get(&seq_id) {
            for &(miss_start, miss_end) in missing.iter() {
                // Check if mask would create a small fragment at start
                if mask_start > miss_start
                    && mask_start < miss_end
                    && mask_start - miss_start < min_fragment_size
                    && mask_start - miss_start > 0
                {
                    extensions_by_seq_id
                        .entry(seq_id)
                        .or_default()
                        .push((miss_start, mask_start));
                }
                // Check if mask would create a small fragment at end
                if mask_end > miss_start
                    && mask_end < miss_end
                    && miss_end - mask_end < min_fragment_size
                    && miss_end - mask_end > 0
                {
                    extensions_by_seq_id
                        .entry(seq_id)
                        .or_default()
                        .push((mask_end, miss_end));
                }
            }
        }
    }

    // Process each overlap
    for (query_interval, cigar, target_interval) in overlaps.drain(..) {
        let seq_id = query_interval.metadata;
        let (mut start, mut end) = if query_interval.first <= query_interval.last {
            (query_interval.first, query_interval.last)
        } else {
            (query_interval.last, query_interval.first)
        };

        // Apply extensions if any
        if let Some(extensions) = extensions_by_seq_id.get(&seq_id) {
            for &(ext_start, ext_end) in extensions {
                // Only apply extension if it's relevant to this interval
                if (ext_end == start) || (ext_start == end) {
                    if ext_start < start {
                        start = ext_start;
                    }
                    if ext_end > end {
                        end = ext_end;
                    }
                }
            }
        }

        // Adjust query and target intervals based on extensions
        let adjusted_query = if query_interval.first <= query_interval.last {
            Interval {
                first: start,
                last: end,
                metadata: query_interval.metadata,
            }
        } else {
            Interval {
                first: end,
                last: start,
                metadata: query_interval.metadata,
            }
        };

        // Proportionally adjust target interval
        let original_span = if query_interval.first <= query_interval.last {
            query_interval.last - query_interval.first
        } else {
            query_interval.first - query_interval.last
        } as f64;

        let new_span = if adjusted_query.first <= adjusted_query.last {
            adjusted_query.last - adjusted_query.first
        } else {
            adjusted_query.first - adjusted_query.last
        } as f64;

        let scale = new_span / original_span;
        let target_span = (target_interval.last - target_interval.first) as f64;

        let adjusted_target = Interval {
            first: target_interval.first,
            last: target_interval.first + (target_span * scale) as i32,
            metadata: target_interval.metadata,
        };

        // Filter against existing masked regions
        if let Some(masks) = masked_regions.get(&seq_id) {
            let mut curr_pos = start;
            let mut unmasked_segments = Vec::new();

            for &(mask_start, mask_end) in masks.iter() {
                if mask_start >= end {
                    break;
                }
                if mask_end <= curr_pos {
                    continue;
                }

                if curr_pos < mask_start {
                    unmasked_segments.push((curr_pos, mask_start));
                }
                curr_pos = mask_end;
            }

            if curr_pos < end {
                unmasked_segments.push((curr_pos, end));
            }

            // Create new intervals for unmasked segments
            for (seg_start, seg_end) in unmasked_segments {
                // Create adjusted intervals
                let segment_ratio = (seg_end - seg_start) as f64 / (end - start) as f64;
                let segment_target_span = target_span * segment_ratio;
                let segment_offset =
                    (seg_start - start) as f64 / (end - start) as f64 * target_span;

                let new_target = Interval {
                    first: target_interval.first + segment_offset as i32,
                    last: target_interval.first + (segment_offset + segment_target_span) as i32,
                    metadata: target_interval.metadata,
                };

                let new_query = if query_interval.first <= query_interval.last {
                    Interval {
                        first: seg_start,
                        last: seg_end,
                        metadata: query_interval.metadata,
                    }
                } else {
                    Interval {
                        first: seg_end,
                        last: seg_start,
                        metadata: query_interval.metadata,
                    }
                };

                result.push((new_query, Vec::new(), new_target));

                // Update masked regions
                let masked = masked_regions.entry(seq_id).or_default();
                masked.insert((seg_start, seg_end));
            }
        } else {
            // No masks for this sequence - keep original interval with extensions
            result.push((adjusted_query, cigar, adjusted_target));

            // Update masked regions
            let masked = masked_regions.entry(seq_id).or_default();
            masked.insert((start, end));
        }
    }

    // Update missing regions based on new masked regions
    for (seq_id, masked) in masked_regions.iter() {
        if let Some(missing) = missing_regions.get_mut(seq_id) {
            let mut new_ranges = Vec::new();

            for &(miss_start, miss_end) in missing.iter() {
                let mut current_ranges = vec![(miss_start, miss_end)];

                for &(mask_start, mask_end) in masked.iter() {
                    let mut next_ranges = Vec::new();

                    for (curr_start, curr_end) in current_ranges {
                        if curr_end <= mask_start || curr_start >= mask_end {
                            next_ranges.push((curr_start, curr_end));
                        } else {
                            if curr_start < mask_start {
                                next_ranges.push((curr_start, mask_start));
                            }
                            if curr_end > mask_end {
                                next_ranges.push((mask_end, curr_end));
                            }
                        }
                    }
                    current_ranges = next_ranges;
                }

                new_ranges.extend(current_ranges);
            }

            missing.ranges.clear();
            for range in new_ranges {
                missing.insert(range);
            }

            if missing.is_empty() {
                missing_regions.remove(seq_id);
            }
        }
    }

    result
}

// Function to extend intervals that are close to sequence boundaries
fn extend_to_close_boundaries(
    overlaps: &mut [(Interval<u32>, Vec<CigarOp>, Interval<u32>)],
    impg: &Impg,
    min_boundary_distance: i32,
) {
    for (query_interval, _, target_interval) in overlaps.iter_mut() {
        let seq_len = impg
            .seq_index
            .get_len_from_id(query_interval.metadata)
            .unwrap() as i32;
        let is_forward = query_interval.first <= query_interval.last;

        if is_forward {
            // Extend to start if close to beginning
            if query_interval.first < min_boundary_distance {
                let shift = query_interval.first;
                query_interval.first = 0;
                target_interval.first -= shift;
            }
            // Extend to end if close to sequence end
            if seq_len - query_interval.last < min_boundary_distance {
                let shift = seq_len - query_interval.last;
                query_interval.last = seq_len;
                target_interval.last += shift;
            }
        } else {
            // Reverse strand - same logic with first/last swapped
            if query_interval.last < min_boundary_distance {
                let shift = query_interval.last;
                query_interval.last = 0;
                target_interval.first -= shift;
            }
            if seq_len - query_interval.first < min_boundary_distance {
                let shift = seq_len - query_interval.first;
                query_interval.first = seq_len;
                target_interval.last += shift;
            }
        }
    }
}

fn write_partition(
    partition_num: usize,
    overlaps: &[(Interval<u32>, Vec<CigarOp>, Interval<u32>)],
    impg: &Impg,
) -> io::Result<()> {
    // Only performing an actual file system write (a system call) when:
    // - The buffer is full
    // - flush() is called
    // - The BufWriter is dropped

    // Create file with buffer
    let file = File::create(format!("partition{}.bed", partition_num))?;
    let mut writer = BufWriter::new(file);

    // Write all overlaps to buffered writer
    for (query_interval, _, _) in overlaps {
        let name = impg.seq_index.get_name(query_interval.metadata).unwrap();
        let (start, end) = if query_interval.first <= query_interval.last {
            (query_interval.first, query_interval.last)
        } else {
            (query_interval.last, query_interval.first)
        };

        writeln!(writer, "{}\t{}\t{}", name, start, end)?;
    }

    // Ensure all data is written to file
    writer.flush()?;
    Ok(())
}
