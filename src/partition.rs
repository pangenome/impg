use crate::impg::CigarOp;
use crate::impg::Impg;
use crate::impg::SortedRanges;
use coitrees::Interval;
use log::{debug, info};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
//use std::time::Instant;

pub fn partition_alignments(
    impg: &Impg,
    window_size: usize,
    starting_sequences_file: Option<&str>,
    selection_mode: &str,
    merge_distance: i32,
    min_missing_size: i32,
    min_boundary_distance: i32,
    max_depth: u16,
    min_transitive_len: i32,
    min_distance_between_ranges: i32,
    debug: bool,
) -> io::Result<()> {
    // Initialize windows from starting sequences if provided
    let mut windows = Vec::<(u32, i32, i32)>::new();
    if let Some(path) = starting_sequences_file {
        // Read sequences from starting-sequences file
        let file = File::open(path).map_err(|e| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Could not open starting sequences file {}: {}", path, e),
            )
        })?;

        let reader = BufReader::new(file);
        let mut starting_sequences = Vec::new();

        // Read all sequences from the file
        for line_result in reader.lines() {
            let line = line_result?;
            if let Some(seq_name) = line.split('\t').next() {
                let trimmed_name = seq_name.trim();
                if trimmed_name.is_empty() || trimmed_name.starts_with('#') {
                    continue; // Skip empty lines and comments
                }

                if let Some(seq_id) = impg.seq_index.get_id(trimmed_name) {
                    let seq_length = impg.seq_index.get_len_from_id(seq_id).unwrap() as i32;
                    starting_sequences.push((seq_id, 0, seq_length));
                } else if debug {
                    debug!(
                        "Sequence {} from starting file not found in index",
                        trimmed_name
                    );
                }
            }
        }

        info!(
            "Read {} sequences from starting file {}",
            starting_sequences.len(),
            path
        );

        // Create windows from starting sequences
        for (seq_id, start, end) in starting_sequences {
            let seq_name = impg.seq_index.get_name(seq_id).unwrap();
            debug!(
                "Creating windows for sequence {} ({}bp)",
                seq_name,
                end - start
            );

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

    // Initialize masked regions
    let mut masked_regions: FxHashMap<u32, SortedRanges> = (0..impg.seq_index.len() as u32)
        .into_par_iter()
        .map(|id| {
            let len = impg.seq_index.get_len_from_id(id).unwrap();
            (id, SortedRanges::new(len as i32, 0))
        })
        .collect();

    // Initialize missing regions from sequence index
    let mut missing_regions: FxHashMap<u32, SortedRanges> = (0..impg.seq_index.len() as u32)
        .into_par_iter()
        .map(|id| {
            let len = impg.seq_index.get_len_from_id(id).unwrap();
            let mut ranges = SortedRanges::new(len as i32, 0);
            ranges.insert((0, len as i32));
            (id, ranges)
        })
        .collect();

    let mut partition_num = 0;
    let mut total_partitioned_length = 0;
    let total_sequence_length: u64 = (0..impg.seq_index.len() as u32)
        .into_par_iter()
        .filter_map(|id| impg.seq_index.get_len_from_id(id))
        .sum::<usize>() as u64;

    info!("Partitioning");

    // If no windows are yet defined, select initial windows based on selection_mode
    if windows.is_empty() {
        select_and_window_sequences(
            &mut windows,
            impg,
            &missing_regions,
            selection_mode,
            window_size,
        )?;
    }

    while !windows.is_empty() {
        if debug {
            debug!("Processing new set of {} windows", windows.len());
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

            //let extend_start = Instant::now();
            if min_boundary_distance > 0 {
                //debug!("  Extending short intervals");
                extend_to_close_boundaries(&mut overlaps, impg, min_boundary_distance);
                debug!(
                    "  Collected {} query overlaps after extending those close to boundaries",
                    overlaps.len()
                );
            }
            //let extend_time = extend_start.elapsed();

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
                //let merge2_start = Instant::now();
                merge_overlaps(&mut overlaps, 0); // Final merge to ensure no overlaps remain
                debug!(
                    "  Collected {} query overlaps after re-merging",
                    overlaps.len()
                );
                //let merge2_time = merge2_start.elapsed();

                //let write_start = Instant::now();
                write_partition(partition_num, &overlaps, impg)?;
                //let write_time = write_start.elapsed();

                partition_num += 1;

                //let calc_start = Instant::now();
                // Calculate current partition length
                let num_regions = overlaps.len();
                let current_partition_length: u64 = overlaps
                    .into_par_iter()
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
                //let calc_time = calc_start.elapsed();

                info!("  Wrote partition {} with {} regions: {} bp this partition ({}), {} bp total ({}) - (query {}:{}-{}, len: {})", 
                    partition_num,
                    num_regions,
                    current_partition_length,   // Current partition size in bp
                    current_percentage_str,     // Current percentage of total sequence
                    total_partitioned_length,   // Total bp written so far
                    total_percentage_str,       // Total percentage of pangenome
                    chrom,
                    start,
                    end,
                    end - start,
                );


                //info!("  Partition {} timings: query={:?}, merge={:?}, merge2={:?}, extend={:?}, mask={:?}, calc={:?}, write={:?}",
                //    partition_num, query_time, merge_time, merge2_time, extend_time, mask_time, calc_time, write_time);
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

        //let window_start = Instant::now();
        select_and_window_sequences(
            &mut windows,
            impg,
            &missing_regions,
            selection_mode,
            window_size,
        )?;
        //let window_time = window_start.elapsed();
        //info!("  select_and_window_sequences={:?}", window_time);
    }

    // Calculate final percentage
    let final_percentage = (total_partitioned_length as f64 / total_sequence_length as f64) * 100.0;
    // Create formatted percentage string with conditional scientific notation
    let final_percentage_str = if final_percentage < 0.0001 {
        format!("{:.4e}%", final_percentage)
    } else {
        format!("{:.4}%", final_percentage)
    };

    info!(
        "Partitioned into {} regions: {} bp total written / {} bp total sequence ({})",
        partition_num, total_partitioned_length, total_sequence_length, final_percentage_str
    );

    Ok(())
}

// Helper function to select and window sequences based on selection_mode
fn select_and_window_sequences(
    windows: &mut Vec<(u32, i32, i32)>,
    impg: &Impg,
    missing_regions: &FxHashMap<u32, SortedRanges>,
    selection_mode: &str,
    window_size: usize,
) -> io::Result<()> {
    let mut ranges_to_window = Vec::new();

    match selection_mode {
        "longest" => {
            // Select longest single missing region
            let longest_region = missing_regions
                .par_iter()
                .flat_map(|(seq_id, ranges)| {
                    // For each sequence, convert its ranges to (seq_id, start, end, length) tuples
                    ranges.iter().map(|&(start, end)| {
                        let length = end - start;
                        (*seq_id, start, end, length)
                    }).collect::<Vec<_>>()
                })
                .max_by(|&(id1, _, _, len1), &(id2, _, _, len2)| {
                    // First compare by length
                    len1.cmp(&len2)
                        // If lengths are equal, compare by sequence ID for deterministic results
                        .then_with(|| id1.cmp(&id2))
                })
                .map(|(seq_id, start, end, _)| (seq_id, start, end));

            if let Some((seq_id, start, end)) = longest_region {
                let seq_name = impg.seq_index.get_name(seq_id).unwrap();

                debug!(
                    "Selected longest missing region {}:{}-{} ({}bp)",
                    seq_name,
                    start,
                    end,
                    end - start
                );

                ranges_to_window.push((seq_id, start, end));
            }
        }
        "total" => {
            // Select sequence with highest total missing
            let seq_with_most_missing = missing_regions
                .par_iter()
                .map(|(seq_id, ranges)| {
                    // i64 to avoid overflow for large ranges
                    let total_missing: i64 = ranges
                        .iter()
                        .map(|&(start, end)| (end - start) as i64)
                        .sum();
                    (*seq_id, total_missing)
                })
                .max_by(|&(id1, missing1), &(id2, missing2)| {
                    // First compare by total missing
                    missing1.cmp(&missing2)
                        // If total missing is equal, compare by sequence ID for deterministic results
                        .then_with(|| id1.cmp(&id2))
                });

            if let Some((seq_id, max_total_missing)) = seq_with_most_missing {
                let seq_length = impg.seq_index.get_len_from_id(seq_id).unwrap() as i32;
                let seq_name = impg.seq_index.get_name(seq_id).unwrap();

                debug!(
                    "Selected sequence {} with most missing sequence ({}bp)",
                    seq_name, max_total_missing
                );

                ranges_to_window.push((seq_id, 0, seq_length));
            }
        }
        mode if mode == "sample"
            || mode == "haplotype"
            || mode.starts_with("sample,")
            || mode.starts_with("haplotype,") =>
        {
            // Parse <field_type>[,<separator>]
            let mut parts = mode.splitn(2, ',');
            let field_type = parts.next().unwrap_or("sample");
            let separator = parts.next().unwrap_or("#");
            let field_count = match field_type {
                "sample" => 1,
                "haplotype" => 2,
                _ => 1, // Default to sample if unrecognized
            };

            // Group only sequences that still have missing regions.
            let prefix_to_seqs: FxHashMap<String, Vec<u32>> = missing_regions.keys()
                .par_bridge()
                .fold(
                || FxHashMap::with_capacity_and_hasher(0, Default::default()),
                |mut map: FxHashMap<String, Vec<u32>>, seq_id: &u32| {
                    let seq_id = *seq_id; // Dereference to get u32
                    if let Some(name) = impg.seq_index.get_name(seq_id) {
                        let mut split = name.split(separator);
                        let prefix = match field_count {
                            1 => split.next().unwrap_or(name).to_string(),
                            2 => {
                                let p1 = split.next().unwrap_or(name);
                                let p2 = split.next().unwrap_or("");
                                format!("{p1}{separator}{p2}")
                            }
                            _ => name.to_string(),
                        };
                        map.entry(prefix).or_default().push(seq_id);
                    }
                    map
                },
            ).reduce(
                || FxHashMap::with_capacity_and_hasher(0, Default::default()),
                |mut map1, map2| {
                    for (key, vec2) in map2 {
                        map1.entry(key).or_default().extend(vec2);
                    }
                    map1
                },
            );
                
            // Find the prefix with the most missing bases in total.
            if !prefix_to_seqs.is_empty() {
                let prefix_with_missing: Vec<(String, i64)> = prefix_to_seqs
                    .par_iter()
                    .map(|(prefix, ids)| {
                        let missing: i64 = ids
                            .par_iter()
                            .filter_map(|id| missing_regions.get(id))
                            .map(|ranges| ranges.iter().map(|&(s, e)| (e - s) as i64).sum::<i64>())
                            .sum();
                        (prefix.clone(), missing)
                    })
                    .collect();

                // Find the prefix with maximum missing (sequential)
                if let Some((best_prefix, best_missing)) = prefix_with_missing
                    .par_iter()
                    .max_by(|&(prefix1, missing1), &(prefix2, missing2)| {
                        // First compare by missing amount
                        missing1.cmp(&missing2)
                            // If missing amounts are equal, compare by prefix name for deterministic results
                            .then_with(|| prefix1.cmp(prefix2))
                    })
                    .map(|(p, m)| (p.clone(), *m))
                {
                    if let Some(best_seqs) = prefix_to_seqs.get(&best_prefix) {
                        if !best_seqs.is_empty() {
                            debug!(
                                "Selected {} sequences from group {} ({} bp missing)",
                                best_seqs.len(),
                                best_prefix,
                                best_missing
                            );

                            // Parallelize fetching lengths
                            let mut seqs_with_len: Vec<(u32, usize)> = best_seqs
                                .par_iter()
                                .filter_map(|&id| {
                                    impg.seq_index.get_len_from_id(id).map(|l| (id, l))
                                })
                                .collect();

                            // Sort by length descending (this remains sequential)
                            seqs_with_len.sort_unstable_by(|a, b| b.1.cmp(&a.1));

                            ranges_to_window.extend(
                                seqs_with_len
                                    .into_iter()
                                    .map(|(id, len)| (id, 0, len as i32)),
                            );
                        }
                    }
                }
            }
        }
        _ => {
            // Invalid selection mode
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid selection mode. Must be 'longest', 'total', 'sample[,sep]', or 'haplotype[,sep]'."
            ));
        }
    }

    // If no more missing regions, we're done
    if ranges_to_window.is_empty() {
        return Ok(());
    }

    // Create windows from ranges
    let new_windows: Vec<(u32, i32, i32)> = ranges_to_window.par_iter()
        .flat_map(|(seq_id, start, end)| {
            let mut range_windows: Vec<(u32, i32, i32)> = Vec::new();
            let mut pos = *start;
            while pos < *end {
                let window_end = std::cmp::min(pos + window_size as i32, *end);
                
                // If this tail-window is too small, merge it into the last one
                if window_end - pos < window_size as i32 && !range_windows.is_empty() {
                    let last = range_windows.last_mut().unwrap();
                    last.2 = *end;
                } else {
                    range_windows.push((*seq_id, pos, window_end));
                }
                
                pos = window_end;
            }
            range_windows
        })
        .collect();
        
    windows.extend(new_windows);
    
    Ok(())
}

fn merge_overlaps(overlaps: &mut Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)>, max_gap: i32) {
    if overlaps.len() > 1 && max_gap >= 0 {
        // Sort by sequence ID and start position
        overlaps.par_sort_by_key(|(query_interval, _, _)| {
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

    // Group overlaps by sequence ID for batch processing
    let mut overlaps_by_seq_id: FxHashMap<u32, Vec<(Interval<u32>, Vec<CigarOp>, Interval<u32>)>> =
        FxHashMap::default();
    for (interval, cigar, target) in overlaps.drain(..) {
        overlaps_by_seq_id
            .entry(interval.metadata)
            .or_default()
            .push((interval, cigar, target));
    }

    // Process each sequence's overlaps as a batch
    for (seq_id, seq_overlaps) in overlaps_by_seq_id {
        // Step 1: Collect all potential extensions in one pass
        let mut extensions = Vec::new();

        if let Some(missing) = missing_regions.get(&seq_id) {
            for (query_interval, _, _) in &seq_overlaps {
                let (mask_start, mask_end) = if query_interval.first <= query_interval.last {
                    (query_interval.first, query_interval.last)
                } else {
                    (query_interval.last, query_interval.first)
                };

                // Use binary search to find relevant missing regions
                let pos = match missing
                    .ranges
                    .binary_search_by_key(&mask_start, |&(s, _)| s)
                {
                    Ok(pos) => pos,
                    Err(pos) => {
                        if pos > 0 {
                            // Check if previous range might overlap
                            let (_prev_start, prev_end) = missing.ranges[pos - 1];
                            if prev_end > mask_start {
                                pos - 1
                            } else {
                                pos
                            }
                        } else {
                            pos
                        }
                    }
                };

                // Check only relevant missing regions
                for i in pos..missing.ranges.len() {
                    let (miss_start, miss_end) = missing.ranges[i];

                    // Stop if we've gone past potential overlaps
                    if miss_start > mask_end {
                        break;
                    }

                    // Check for small fragment at start
                    if mask_start > miss_start
                        && mask_start < miss_end
                        && mask_start - miss_start < min_fragment_size
                        && mask_start - miss_start > 0
                    {
                        extensions.push((miss_start, mask_start));
                    }

                    // Check for small fragment at end
                    if mask_end > miss_start
                        && mask_end < miss_end
                        && miss_end - mask_end < min_fragment_size
                        && miss_end - mask_end > 0
                    {
                        extensions.push((mask_end, miss_end));
                    }
                }
            }
        }

        // Step 2: Sort and merge extensions
        if !extensions.is_empty() {
            extensions.sort_by_key(|&(start, _)| start);

            let mut merged_extensions = Vec::with_capacity(extensions.len());
            let mut current = extensions[0];

            for &(start, end) in &extensions[1..] {
                if start <= current.1 {
                    // Merge overlapping extensions
                    current.1 = current.1.max(end);
                } else {
                    merged_extensions.push(current);
                    current = (start, end);
                }
            }
            merged_extensions.push(current);
            extensions = merged_extensions;
        }

        // Step 3: Process each overlap with precomputed extensions
        let mut new_masked_ranges = Vec::new();

        for (query_interval, cigar, target_interval) in seq_overlaps {
            let (mut start, mut end) = if query_interval.first <= query_interval.last {
                (query_interval.first, query_interval.last)
            } else {
                (query_interval.last, query_interval.first)
            };

            // Apply relevant extensions
            for &(ext_start, ext_end) in &extensions {
                // Check if extension applies to this interval's boundaries
                if (ext_end >= start && ext_start <= start) || (ext_start <= end && ext_end >= end)
                {
                    if ext_start < start {
                        start = ext_start;
                    }
                    if ext_end > end {
                        end = ext_end;
                    }
                }
            }

            // Add to masked ranges for this sequence
            new_masked_ranges.push((start, end));

            // Adjust query interval based on extensions
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

            // Find unmasked segments efficiently using binary search
            if let Some(masks) = masked_regions.get(&seq_id) {
                let mut curr_pos = start;
                let mut unmasked_segments = Vec::new();

                // Binary search for starting position in masks
                let mut idx = match masks.ranges.binary_search_by_key(&curr_pos, |&(s, _)| s) {
                    Ok(pos) => pos,
                    Err(pos) => {
                        if pos > 0 {
                            // Check if previous range might overlap
                            let (_, prev_end) = masks.ranges[pos - 1];
                            if prev_end > curr_pos {
                                pos - 1
                            } else {
                                pos
                            }
                        } else {
                            pos
                        }
                    }
                };

                // Process only relevant masks
                while idx < masks.ranges.len() {
                    let (mask_start, mask_end) = masks.ranges[idx];

                    if mask_start > end {
                        break;
                    } // No more relevant masks
                    if mask_end <= curr_pos {
                        idx += 1;
                        continue;
                    }

                    if curr_pos < mask_start {
                        unmasked_segments.push((curr_pos, mask_start));
                    }

                    curr_pos = curr_pos.max(mask_end);
                    idx += 1;

                    if curr_pos >= end {
                        break;
                    }
                }

                // Check for final unmasked segment
                if curr_pos < end {
                    unmasked_segments.push((curr_pos, end));
                }

                // Create adjusted intervals for each unmasked segment
                for (seg_start, seg_end) in unmasked_segments {
                    // Calculate proportional target interval
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
                }
            } else {
                // No existing masks - keep the entire interval
                result.push((adjusted_query, cigar, adjusted_target));
            }
        }

        // Step 4: Update masked regions for this sequence
        let masked = masked_regions.entry(seq_id).or_default();
        for (start, end) in new_masked_ranges {
            masked.insert((start, end));
        }

        // Step 5: Efficiently update missing regions for this sequence
        if let Some(missing) = missing_regions.get_mut(&seq_id) {
            if let Some(masked) = masked_regions.get(&seq_id) {
                // Clone missing ranges to avoid borrow issues
                let original_missing = missing.ranges.clone();

                // Clear current missing ranges and rebuild
                missing.ranges.clear();

                for (miss_start, miss_end) in original_missing {
                    let mut current = miss_start;

                    // Find masks that could affect this missing range with binary search
                    let mut idx = match masked.ranges.binary_search_by_key(&miss_start, |&(s, _)| s)
                    {
                        Ok(pos) => pos,
                        Err(pos) => {
                            if pos > 0 {
                                let (_, prev_end) = masked.ranges[pos - 1];
                                if prev_end > miss_start {
                                    pos - 1
                                } else {
                                    pos
                                }
                            } else {
                                pos
                            }
                        }
                    };

                    // Process all relevant masks
                    while idx < masked.ranges.len() && current < miss_end {
                        let (mask_start, mask_end) = masked.ranges[idx];

                        if mask_start > miss_end {
                            break;
                        }
                        if mask_end <= current {
                            idx += 1;
                            continue;
                        }

                        if current < mask_start {
                            // Found a gap - this part is still missing
                            missing.insert((current, mask_start));
                        }

                        current = current.max(mask_end);
                        idx += 1;
                    }

                    // Check for remaining gap at the end
                    if current < miss_end {
                        missing.insert((current, miss_end));
                    }
                }

                // Remove sequence if no more missing regions
                if missing.is_empty() {
                    missing_regions.remove(&seq_id);
                }
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
