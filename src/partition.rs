use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use crate::impg::Impg;
use coitrees::Interval;
use crate::impg::CigarOp;
use crate::impg::SortedRanges;
use log::{debug, info};
//use std::time::Instant;

pub fn partition_alignments(
    impg: &Impg,
    window_size: usize,
    sequence_prefix: &str,
    min_length: usize,
    merge_distance: usize,
    max_depth: Option<u32>,
    debug: bool,
) -> io::Result<()> {
    // Get all sequences with the given prefix
    let mut sample_regions = Vec::new();
    for seq_id in 0..impg.seq_index.len() as u32 {
        let seq_name = impg.seq_index.get_name(seq_id).unwrap();
        if seq_name.starts_with(sequence_prefix) {
            let seq_length = impg.seq_index.get_len_from_id(seq_id).unwrap();
            sample_regions.push((seq_id, 0, seq_length));
        }
    }
    if sample_regions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("No sequences with prefix {} found in sequence index", sequence_prefix),
        ));
    }
    // Natural sort by sequence name
    sample_regions.sort_by(|a, b| {
        let chrom_a = impg.seq_index.get_name(a.0).unwrap();
        let chrom_b = impg.seq_index.get_name(b.0).unwrap();
        natord::compare(chrom_a, chrom_b)
    });

    if debug {
        debug!("Found {} sequences with prefix {}", sample_regions.len(), sequence_prefix);
        for (seq_id, start, end) in &sample_regions {
            let chrom = impg.seq_index.get_name(*seq_id).unwrap();
            debug!("  Sequence: {}:{}-{}", chrom, start, end);
        }
    }

    // Create windows from sample regions
    let mut windows = Vec::<(u32, i32, i32)>::new();
    for (seq_id, start, end) in sample_regions {
        let mut pos = start;
        while pos < end as i32 {
            let window_end = std::cmp::min(pos + window_size as i32, end as i32);
            windows.push((seq_id, pos, window_end));
            pos = window_end;
        }
    }

    if debug {
        debug!("Starting with {} windows:", windows.len());
        for (chrom, start, end) in &windows {
            debug!("  Window: {}:{}-{}", chrom, start, end);
        }
    }

    // Initialize masked regions
    let mut masked_regions: FxHashMap<u32, SortedRanges> = FxHashMap::default();
    
    // Initialize missing regions from sequence index
    let mut missing_regions: FxHashMap<u32, SortedRanges> = (0..impg.seq_index.len() as u32)
        .map(|id| {
            let len = impg.seq_index.get_len_from_id(id).unwrap();
            let mut ranges = SortedRanges::new();
            ranges.insert((0, len as i32));
            (id, ranges)
        })
        .collect();

    let mut partition_num = 0;

    info!("Partitioning");

    while !windows.is_empty() {
        for (seq_id, start, end) in windows.iter() {     
            let chrom = impg.seq_index.get_name(*seq_id).unwrap();

            if debug {
                debug!("Processing new window set");

                debug!("  Querying region {}:{}-{}", chrom, start, end);

                debug!("  Missing {} regions in {} sequences", 
                    missing_regions.values().map(|ranges| ranges.len()).sum::<usize>(),
                    missing_regions.len()
                );
                for (chrom, ranges) in &missing_regions {
                    for &(start, end) in ranges.iter() {
                        debug!("    Region: {}:{}-{}", chrom, start, end);
                    }
                }
                
                debug!("  Masked {} regions in {} sequences", 
                    masked_regions.values().map(|ranges| ranges.len()).sum::<usize>(),
                    masked_regions.len()
                );
                for (chrom, ranges) in &masked_regions {
                    for &(start, end) in ranges.iter() {
                        debug!("    Region: {}:{}-{}", chrom, start, end);
                    }
                }
            }

            // Query overlaps for current window
            //let query_start = Instant::now();
            let mut overlaps = impg.query_transitive(*seq_id, *start, *end, Some(&masked_regions), max_depth);
            //let query_time = query_start.elapsed();
            debug!("  Collected {} query overlaps", overlaps.len());

            // Ignore CIGAR strings and target intervals.
            debug!("  Merging overlaps closer than {}bp", merge_distance); // bedtools sort | bedtools merge -d merge_distance
            //let merge_start = Instant::now();
            merge_overlaps(&mut overlaps, merge_distance as i32);
            //let merge_time = merge_start.elapsed();
            debug!("  Collected {} query overlaps after merging", overlaps.len());

            debug!("  Excluding masked regions"); // bedtools subtract -a "partition$num.tmp.bed" -b "$MASK_BED"
            //let mask_start = Instant::now();
            overlaps = subtract_masked_regions(&mut overlaps, &masked_regions);
            //let mask_time = mask_start.elapsed();

            if !overlaps.is_empty() {
                debug!("  Collected {} query overlaps in partition {}", overlaps.len(), partition_num);

                debug!("  Updating mask and missing regions");
                //let update_start = Instant::now();
                update_masked_and_missing_regions(&mut masked_regions, &mut missing_regions, &overlaps);            
                //let update_time = update_start.elapsed();

                debug!("  Extending short intervals");
                //let extend_start = Instant::now();
                extend_short_intervals(&mut overlaps, impg, min_length as i32);
                //let extend_time = extend_start.elapsed();

                info!("  Writing partition {} with {} regions (query {}:{}-{})", partition_num, overlaps.len(), chrom, start, end);
                //let write_start = Instant::now();
                write_partition(partition_num, &overlaps, impg)?;
                //let write_time = write_start.elapsed();

                partition_num += 1;

                // info!("  Partition {} timings: query={:?}, merge={:?}, mask={:?}, update={:?}, extend={:?}, write={:?}",
                //     partition_num, query_time, merge_time, mask_time, update_time, extend_time, write_time);
            } else {
                debug!("  No overlaps found for region {}:{}-{}", chrom, start, end);
            }
        }

        // If no missing regions remain, we're done
        if missing_regions.is_empty() {
            break;
        }

        // Find longest remaining region with smallest seq_id and start position
        let mut longest_region: Option<(u32, i32, i32)> = None;
        let mut max_length = 0;
        // Scan through missing regions
        for (seq_id, ranges) in missing_regions.iter() {
            for &(start, end) in ranges.iter() {
                let length = end - start;
                if length > max_length || (length == max_length && longest_region.map_or(true, |(curr_seq_id, curr_start, _)| 
                    (*seq_id < curr_seq_id) || (*seq_id == curr_seq_id && start < curr_start))) 
                {
                    max_length = length;
                    longest_region = Some((*seq_id, start, end));
                }
            }
        }

        // Clear existing windows but keep the allocation
        windows.clear();

        // Create new windows from the selected region
        if let Some((seq_id, start, end)) = longest_region {
            let mut pos = start;
            while pos < end {
                let window_end = std::cmp::min(pos + window_size as i32, end);
                windows.push((seq_id, pos, window_end));
                pos = window_end;
            }
        }
    }

    info!("Partitioned into {} regions", partition_num);

    Ok(())
}

fn merge_overlaps(
    overlaps: &mut Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)>,
    max_gap: i32
) {
    if overlaps.len() > 1 {
        // Sort by sequence ID and start position
        overlaps.sort_unstable_by_key(|(query_interval, _, _)| {
            (query_interval.metadata, std::cmp::min(query_interval.first, query_interval.last))
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

fn subtract_masked_regions(
    overlaps: &mut Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)>,
    masked_regions: &FxHashMap<u32, SortedRanges>
) -> Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)> {
    let mut result = Vec::new();

    for (query_interval, cigar, target_interval) in overlaps.drain(..) {       
        // Get masked regions for this sequence
        if let Some(masks) = masked_regions.get(&query_interval.metadata) {
            let (start, end) = if query_interval.first <= query_interval.last {
                (query_interval.first, query_interval.last)
            } else {
                (query_interval.last, query_interval.first)
            };

            // Track unmasked segments
            let mut curr_pos = start;
            let mut unmasked_segments = Vec::new();

            for &(mask_start, mask_end) in masks.iter() {
                // If mask starts after current segment ends, we're done
                if mask_start >= end {
                    break;
                }

                // If mask ends before current position, skip it
                if mask_end <= curr_pos {
                    continue;
                }

                // Add unmasked segment before mask if it exists
                if curr_pos < mask_start {
                    unmasked_segments.push((curr_pos, mask_start));
                }

                // Move current position to end of mask
                curr_pos = mask_end;
            }

            // Add final unmasked segment if needed
            if curr_pos < end {
                unmasked_segments.push((curr_pos, end));
            }

            // Create new intervals for unmasked segments
            for (seg_start, seg_end) in unmasked_segments {
                let new_query = Interval {
                    first: seg_start,
                    last: seg_end,
                    metadata: query_interval.metadata,
                };
                
                // Adjust target interval proportionally
                let query_len = (end - start) as f64;
                let seg_frac_start = (seg_start - start) as f64 / query_len;
                let seg_frac_end = (seg_end - start) as f64 / query_len;
                
                let target_span = (target_interval.last - target_interval.first) as f64;
                let new_target = Interval {
                    first: target_interval.first + (seg_frac_start * target_span) as i32,
                    last: target_interval.first + (seg_frac_end * target_span) as i32,
                    metadata: target_interval.metadata,
                };

                result.push((new_query, cigar.clone(), new_target));
            }
        } else {
            // No masks for this sequence - keep original interval
            result.push((query_interval, cigar, target_interval));
        }
    }

    result
}

fn update_masked_and_missing_regions(
    masked_regions: &mut FxHashMap<u32, SortedRanges>,
    missing_regions: &mut FxHashMap<u32, SortedRanges>,
    overlaps: &Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)>
) {
    // First, collect all new regions to be masked by sequence
    let mut new_masks: FxHashMap<u32, Vec<(i32, i32)>> = FxHashMap::default();
    for (query_interval, _, _) in overlaps {
        let (start, end) = if query_interval.first <= query_interval.last {
            (query_interval.first, query_interval.last)
        } else {
            (query_interval.last, query_interval.first)
        };
        new_masks.entry(query_interval.metadata).or_default().push((start, end));
    }

    // Update masked regions with new masks
    for (seq_id, ranges) in new_masks {
        let masked = masked_regions.entry(seq_id).or_default();
        for range in ranges {
            masked.insert((range.0, range.1));
        }

        // Update missing regions for this sequence
        if let Some(missing) = missing_regions.get_mut(&seq_id) {
            // Create temporary vector to store new ranges
            let mut new_ranges = Vec::new();

            // For each missing range, subtract all masked ranges
            for &(miss_start, miss_end) in missing.iter() {
                let mut current_ranges = vec![(miss_start, miss_end)];
                
                for &(mask_start, mask_end) in masked.iter() {
                    let mut next_ranges = Vec::new();
                    
                    for (curr_start, curr_end) in current_ranges {
                        // If current range doesn't overlap with mask
                        if curr_end <= mask_start || curr_start >= mask_end {
                            next_ranges.push((curr_start, curr_end));
                        }
                        // If ranges overlap
                        else {
                            // Add portion before mask if it exists
                            if curr_start < mask_start {
                                next_ranges.push((curr_start, mask_start));
                            }
                            // Add portion after mask if it exists
                            if curr_end > mask_end {
                                next_ranges.push((mask_end, curr_end));
                            }
                        }
                    }
                    current_ranges = next_ranges;
                }
                
                new_ranges.extend(current_ranges);
            }
            
            // Clear existing ranges and insert new ones in-place
            missing.ranges.clear();
            for range in new_ranges {
                missing.insert(range);
            }
            
            // Remove sequence from missing_regions if no ranges remain
            if missing.is_empty() {
                missing_regions.remove(&seq_id);
            }
        }
    }
}

fn extend_short_intervals(
    overlaps: &mut Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)>,
    impg: &Impg,
    min_length: i32,
) {   
    for (query_interval, _, target_interval) in overlaps.iter_mut() {
        let len = (query_interval.last - query_interval.first).abs();
        
        if len < min_length {
            let extension_needed = min_length - len;
            // Add 1 to the first side if extension_needed is odd
            let extension_per_side = extension_needed / 2;
            let first_side_extension = extension_per_side + (extension_needed % 2);
            let second_side_extension = extension_per_side;
            
            let seq_len = impg.seq_index.get_len_from_id(query_interval.metadata).unwrap() as i32;
            if query_interval.first <= query_interval.last {
                // Forward strand
                let mut start_extension = first_side_extension;
                let mut end_extension = second_side_extension;
                
                // If we can't extend fully on the start side, add the remainder to the end
                if query_interval.first < start_extension {
                    let remaining = start_extension - query_interval.first;
                    start_extension = query_interval.first;
                    end_extension += remaining;
                }
                
                // If we can't extend fully on the end side, add the remainder to the start
                if query_interval.last + end_extension > seq_len {
                    let remaining = (query_interval.last + end_extension) - seq_len;
                    end_extension = seq_len - query_interval.last;
                    // Only use remaining if we can extend more at the start
                    if query_interval.first >= remaining {
                        start_extension += remaining;
                    }
                }
                
                query_interval.first -= start_extension;
                query_interval.last += end_extension;
                target_interval.first -= start_extension;
                target_interval.last += end_extension;
            } else {
                // Reverse strand
                let mut start_extension = first_side_extension;
                let mut end_extension = second_side_extension;
                
                // For reverse strand, query_interval.last is the start position
                if query_interval.last < start_extension {
                    let remaining = start_extension - query_interval.last;
                    start_extension = query_interval.last;
                    end_extension += remaining;
                }
                
                if query_interval.first + end_extension > seq_len {
                    let remaining = (query_interval.first + end_extension) - seq_len;
                    end_extension = seq_len - query_interval.first;
                    // Only use remaining if we can extend more at the start
                    if query_interval.last >= remaining {
                        start_extension += remaining;
                    }
                }
                
                query_interval.last -= start_extension;
                query_interval.first += end_extension;
                target_interval.first -= start_extension;
                target_interval.last += end_extension;
            }
        }
    }

    overlaps.sort_unstable_by_key(|(query_interval, _, _)| {
        (query_interval.metadata, std::cmp::min(query_interval.first, query_interval.last))
    });
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
