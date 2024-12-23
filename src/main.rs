use clap::Parser;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::num::NonZeroUsize;
use noodles::bgzf;
use impg::impg::{Impg, SerializableImpg, AdjustedInterval};
use coitrees::IntervalTree;
use impg::paf;
use rayon::ThreadPoolBuilder;
use std::io::BufRead;
use log::{debug, info, warn, error};
use std::collections::HashSet;

/// Common options shared between all commands
#[derive(Parser, Debug)]
struct CommonOpts {
    /// Path to the PAF file. If specified without an index, the tool will look for or generate an associated index file.
    #[clap(short = 'p', long, value_parser)]
    paf_file: String,

    /// Force the regeneration of the index, even if it already exists.
    #[clap(short = 'I', long, action)]
    force_reindex: bool,

    /// Number of threads for parallel processing.
    #[clap(short = 't', long, value_parser, default_value_t = NonZeroUsize::new(1).unwrap())]
    num_threads: NonZeroUsize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[clap(short, long, default_value = "0")]
    verbose: u8,
}

/// Command-line tool for querying overlaps in PAF files.
#[derive(Parser, Debug)]
#[command(author, version, about, disable_help_subcommand = true)]
enum Args {
    /// Partition the alignment
    Partition {
        #[clap(flatten)]
        common: CommonOpts,

        /// Window size for partitioning.
        #[clap(short = 'w', long, value_parser)]
        window_size: usize,
        
        /// Sequence name prefix to start - all sequences starting with this prefix will be included
        #[clap(short = 's', long, value_parser)]
        sequence_prefix: String,

        /// Minimum length for intervals (default: 5000).
        #[clap(long, value_parser, default_value_t = 5000)]
        min_length: usize,
    },
    /// Query overlaps in the alignment.
    Query {
        #[clap(flatten)]
        common: CommonOpts,
        
        /// Target range in the format `seq_name:start-end`.
        #[clap(short = 'r', long, value_parser)]
        target_range: Option<String>,

        /// Path to the BED file containing target regions.
        #[clap(short = 'b', long, value_parser)]
        target_bed: Option<String>,

        /// Enable transitive overlap requests.
        #[clap(short = 'x', long, action)]
        transitive: bool,

        /// Output results in PAF format.
        #[clap(short = 'P', long, action)]
        output_paf: bool,

        /// Check the projected intervals, reporting the wrong ones (slow, useful for debugging).
        #[clap(short = 'c', long, action)]
        check_intervals: bool,
    },
    /// Print alignment statistics.
    Stats {
        #[clap(flatten)]
        common: CommonOpts,
    },
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    match args {
        Args::Partition {
            common,
            window_size,
            sequence_prefix,
            min_length,
        } => {
            let impg = initialize_impg(&common)?;
            partition_alignments(&impg, window_size, &sequence_prefix, min_length, common.verbose > 1)?;
        },
        Args::Query {
            common,
            target_range,
            target_bed,
            transitive,
            output_paf,
            check_intervals,
        } => {
            let impg = initialize_impg(&common)?;

            if let Some(target_range) = target_range {
                let (target_name, target_range) = parse_target_range(&target_range)?;
                let results = perform_query(&impg, &target_name, target_range, transitive);
                if check_intervals {
                    let invalid_cigars = impg::impg::check_intervals(&impg, &results);
                    if !invalid_cigars.is_empty() {
                        for (row, error_reason) in invalid_cigars {
                            error!("{}; {}", error_reason, row);
                        }
                        panic!("Invalid intervals encountered.");
                    }
                }
                
                if output_paf {
                    // Skip the first element (the input range) for PAF
                    output_results_paf(&impg, results.into_iter().skip(1), &target_name, None);
                } else {
                    output_results_bed(&impg, results);
                }
            } else if let Some(target_bed) = target_bed {
                let targets = parse_bed_file(&target_bed)?;
                for (target_name, target_range, name) in targets {
                    let results = perform_query(&impg, &target_name, target_range, transitive);
                    if check_intervals {
                        let invalid_cigars = impg::impg::check_intervals(&impg, &results);
                        if !invalid_cigars.is_empty() {
                            for (row, error_reason) in invalid_cigars {
                                error!("{}; {}", error_reason, row);
                            }
                            panic!("Invalid intervals encountered.");
                        }
                    }

                    // Skip the first element (the input range) for both PAF and BEDPE
                    let results_iter = results.into_iter().skip(1);
                    if output_paf {
                        output_results_paf(&impg, results_iter, &target_name, name);
                    } else {
                        output_results_bedpe(&impg, results_iter, &target_name, name);
                    }
                }
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --target-range or --target-bed must be provided for query subcommand",
                ));
            }
        },
        Args::Stats { common } => {
            let impg = initialize_impg(&common)?;

            print_stats(&impg);
        },
    }

    Ok(())
}

/// Initialize thread pool and load/generate index based on common options
fn initialize_impg(common: &CommonOpts) -> io::Result<Impg> {
    // Initialize logger based on verbosity
    env_logger::Builder::new()
    .filter_level(match common.verbose {
        0 => log::LevelFilter::Error,
        1 => log::LevelFilter::Info,
        _ => log::LevelFilter::Debug,
    })
    .init();

    // Configure thread pool
    ThreadPoolBuilder::new()
        .num_threads(common.num_threads.into())
        .build_global()
        .unwrap();

    // Load or generate index
    if common.force_reindex {
        generate_index(&common.paf_file, common.num_threads)
    } else {
        load_or_generate_index(&common.paf_file, common.num_threads)
    }
}

fn load_or_generate_index(paf_file: &str, num_threads: NonZeroUsize) -> io::Result<Impg> {
    let index_file = format!("{}.impg", paf_file);
    if std::path::Path::new(&index_file).exists() {
        load_index(paf_file)
    } else {
        generate_index(paf_file, num_threads)
    }
}

fn load_index(paf_file: &str) -> io::Result<Impg> {
    let index_file = format!("{}.impg", paf_file);
    
    let paf_file_metadata = std::fs::metadata(paf_file)?;
    let index_file_metadata = std::fs::metadata(index_file.clone())?;
    if let (Ok(paf_file_ts), Ok(index_file_ts)) = (paf_file_metadata.modified(), index_file_metadata.modified()) {
        if paf_file_ts > index_file_ts
        {
            warn!("WARNING:\tPAF file has been modified since impg index creation.");
        }
    } else {
        warn!("WARNING:\tUnable to compare timestamps of PAF file and impg index file. PAF file may have been modified since impg index creation.");
    }

    let file = File::open(index_file)?;
    let reader = BufReader::new(file);
    let serializable: SerializableImpg = bincode::deserialize_from(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to deserialize index: {:?}", e)))?;
    Ok(Impg::from_paf_and_serializable(paf_file, serializable))
}

fn generate_index(paf_file: &str, num_threads: NonZeroUsize) -> io::Result<Impg> {
    let file = File::open(paf_file)?;
    let reader: Box<dyn io::Read> = if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
        Box::new(bgzf::MultithreadedReader::with_worker_count(num_threads, file))
    } else {
        Box::new(file)
    };
    let reader = BufReader::new(reader);
    let records = paf::parse_paf(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to parse PAF records: {:?}", e)))?;
    let impg = Impg::from_paf_records(&records, paf_file).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to create index: {:?}", e)))?;

    let index_file = format!("{}.impg", paf_file);
    let serializable = impg.to_serializable();
    let file = File::create(index_file)?;
    let writer = BufWriter::new(file);
    bincode::serialize_into(writer, &serializable).map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Failed to serialize index: {:?}", e)))?;

    Ok(impg)
}

fn partition_alignments(
    impg: &Impg,
    window_size: usize,
    sequence_prefix: &str,
    min_length: usize,
    debug: bool,
) -> io::Result<()> {
    // Get all sequences with the given prefix
    let mut sample_regions = Vec::new();
    for seq_id in 0..impg.seq_index.len() as u32 {
        let seq_name = impg.seq_index.get_name(seq_id).unwrap();
        if seq_name.starts_with(sequence_prefix) {
            let seq_length = impg.seq_index.get_len_from_id(seq_id).unwrap();
            sample_regions.push((seq_name.to_string(), 0, seq_length));
        }
    }
    if sample_regions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("No sequences with prefix {} found in sequence index", sequence_prefix),
        ));
    }

    if debug {
        debug!("Starting sequences:");
        for (chrom, start, end) in &sample_regions {
            debug!("  Region: {}:{}-{}", chrom, start, end);
        }
    }

    // Create windows from sample regions
    let mut windows = Vec::new();
    for (chrom, start, end) in sample_regions {
        let mut pos = start;
        while pos < end {
            let window_end = std::cmp::min(pos + window_size, end);
            windows.push((chrom.clone(), pos, window_end));
            pos = window_end;
        }
    }

    if debug {
        debug!("Starting windows:");
        for (chrom, start, end) in &windows {
            debug!("  Window: {}:{}-{}", chrom, start, end);
        }
    }

    // Initialize masked regions
    let mut masked_regions = HashSet::new();
    
    // Initialize missing regions from sequence index
    let mut missing_regions: HashSet<_> = (0..impg.seq_index.len() as u32)
        .map(|id| {
            let name = impg.seq_index.get_name(id).unwrap();
            let len = impg.seq_index.get_len_from_id(id).unwrap();
            (name.to_string(), 0, len)
        })
        .collect();

    let mut partition_num = 0;
    
    while !windows.is_empty() {
        debug!("Processing new window set");
        
        if debug {
            debug!("  Missing regions:");
            for (chrom, start, end) in &missing_regions {
                debug!("    Region: {}:{}-{}", chrom, start, end);
            }
        }

        let mut new_windows = Vec::new();
        
        for (chrom, start, end) in windows.iter() {
            let region = format!("{}:{}-{}", chrom, start, end);
            debug!("  Querying region {}", region);
            
            // Query overlaps for current window
            let overlaps = impg.query_transitive(
                impg.seq_index.get_id(chrom).unwrap(), 
                *start as i32, 
                *end as i32
            );

            // TODO: Merge overlaps that are close to each other
            // bedtools sort | bedtools merge -d 10000

            // Apply mask
            // bedtools subtract -a "partition$num.tmp.bed" -b "$MASK_BED"
            // Convert overlaps to regions, excluding masked areas
            let mut partition_regions = Vec::new();
            for (interval, _, _) in overlaps.into_iter().skip(1) { // Skip input interval
                let name = impg.seq_index.get_name(interval.metadata).unwrap().to_string();
                let region = (name, interval.first as usize, interval.last as usize);
                if !masked_regions.contains(&region) {
                    partition_regions.push(region);
                }
            }

            if debug {
                debug!("  Collected {} overlaps for region {}:", partition_regions.len(), region);
                // for (chrom, start, end) in &partition_regions {
                //     debug!("    Region: {}:{}-{}", chrom, start, end);
                // }
            }

            if !partition_regions.is_empty() {
                info!("  Processing partition {}", partition_num);
                
                // Write partition to file
                let mut partition_file = File::create(format!("partition{}.bed", partition_num))?;
                for (name, start, end) in &partition_regions {
                    writeln!(partition_file, "{}\t{}\t{}", name, start, end)?;
                }

                // Update masked regions and missing regions
                for region in partition_regions {
                    masked_regions.insert(region.clone());
                    missing_regions.remove(&region);
                }

                // Handle short intervals
                let mut short_regions = Vec::new();
                let mut long_regions = Vec::new();
                
                for region in &masked_regions {
                    let len = region.2 - region.1;
                    if len < min_length {
                        short_regions.push(region.clone());
                    } else {
                        long_regions.push(region.clone());
                    }
                }

                // Extend short intervals
                for region in short_regions {
                    if let Some(seq_id) = impg.seq_index.get_id(&region.0) {
                        if let Some(seq_len) = impg.seq_index.get_len_from_id(seq_id) {
                            let new_start = if region.1 >= min_length { region.1 - min_length } else { 0 };
                            let new_end = std::cmp::min(region.2 + min_length, seq_len);
                            long_regions.push((region.0, new_start, new_end));
                        }
                    }
                }

                partition_num += 1;
            } else {
                panic!("No overlaps found for region {}", region);
            }
        }

        // If no missing regions remain, we're done
        if missing_regions.is_empty() {
            break;
        }

        // Find longest remaining region for new windows
        if let Some(longest_region) = missing_regions.iter()
            .max_by_key(|r| r.2 - r.1) {
                let (chrom, start, end) = longest_region;
                let mut pos = *start;
                while pos < *end {
                    let window_end = std::cmp::min(pos + window_size, *end);
                    new_windows.push((chrom.clone(), pos, window_end));
                    pos = window_end;
                }
        }

        windows = new_windows;
    }

    Ok(())
}

fn parse_target_range(target_range: &str) -> io::Result<(String, (i32, i32))> {
    let parts: Vec<&str> = target_range.rsplitn(2, ':').collect();
    if parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Target range format should be `seq_name:start-end`"));
    }

    let (start, end) = parse_range(&parts[0].split('-').collect::<Vec<_>>())?;
    Ok((parts[1].to_string(), (start, end)))
}

fn parse_bed_file(bed_file: &str) -> io::Result<Vec<(String, (i32, i32), Option<String>)>> {
    let file = File::open(bed_file)?;
    let reader = BufReader::new(file);
    let mut ranges = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid BED file format"));
        }

        let (start, end) = parse_range(&parts[1..=2])?;
        let name = parts.get(3).map(|s| s.to_string());
        ranges.push((parts[0].to_string(), (start, end), name));
    }

    Ok(ranges)
}

fn parse_range(range_parts: &[&str]) -> io::Result<(i32, i32)> {
    if range_parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Range format should be `start-end`"));
    }

    let start = range_parts[0].parse::<i32>().map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid start value"))?;
    let end = range_parts[1].parse::<i32>().map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid end value"))?;

    if start >= end {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Start value must be less than end value"));
    }

    Ok((start, end))
}

fn perform_query(impg: &Impg, target_name: &str, target_range: (i32, i32), transitive: bool) -> Vec<AdjustedInterval> {
    let (target_start, target_end) = target_range;
    let target_id = impg.seq_index.get_id(target_name).expect("Target name not found in index");
    let target_length = impg.seq_index.get_len_from_id(target_id).expect("Target length not found in index");
    if target_end > target_length as i32 {
        panic!("Target range end ({}) exceeds the target sequence length ({})", target_end, target_length);
    }
    if transitive {
        impg.query_transitive(target_id, target_start, target_end)
    } else {
        impg.query(target_id, target_start, target_end)
    }
}

fn output_results_bed(impg: &Impg, results: Vec<AdjustedInterval>) {
    for (overlap, _, _) in results {
        let overlap_name = impg.seq_index.get_name(overlap.metadata).unwrap();
        let (first, last, strand) = if overlap.first <= overlap.last {
            (overlap.first, overlap.last, '+')
        } else {
            (overlap.last, overlap.first, '-')
        };
        println!("{}\t{}\t{}\t.\t.\t{}", overlap_name, first, last, strand);
    }
}

fn output_results_bedpe<I>(impg: &Impg, results: I, target_name: &str, name: Option<String>)
where
    I: Iterator<Item = AdjustedInterval>
{
    for (overlap_query, _, overlap_target) in results {
        let overlap_name = impg.seq_index.get_name(overlap_query.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t+",
                 overlap_name, first, last,
                 target_name, overlap_target.first, overlap_target.last,
                 name.as_deref().unwrap_or("."), strand);
    }
}

fn output_results_paf<I>(impg: &Impg, results: I, target_name: &str, name: Option<String>)
where
    I: Iterator<Item = AdjustedInterval>
{
    let target_length = impg.seq_index.get_len_from_id(impg.seq_index.get_id(target_name).unwrap()).unwrap();  
    for (overlap_query, cigar, overlap_target) in results {
        let overlap_name = impg.seq_index.get_name(overlap_query.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };

        let query_length = impg.seq_index.get_len_from_id(overlap_query.metadata).unwrap();  

        let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len) = cigar.iter().fold((0, 0, 0, 0, 0, 0, 0), |(m, mm, i, i_bp, d, d_bp, bl), op| {
            let len = op.len();
            match op.op() {
                'M' => (m + len, mm, i, i_bp, d, d_bp, bl + len), // We overestimate num. of matches by assuming 'M' represents matches for simplicity
                '=' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
                'X' => (m, mm + len, i, i_bp, d, d_bp, bl + len),
                'I' => (m, mm, i + 1, i_bp + len, d, d_bp, bl + len),
                'D' => (m, mm, i, i_bp, d + 1, d_bp + len, bl + len),
                _ => (m, mm, i, i_bp, d, d_bp, bl),
            }
        });
        let gap_compressed_identity = (matches as f64) / (matches + mismatches + insertions + deletions) as f64;
        
        let edit_distance = mismatches + inserted_bp + deleted_bp;
        let block_identity = (matches as f64) / (matches + edit_distance) as f64;

        // Format bi and gi fields without trailing zeros
        let gi_str = format!("{:.6}", gap_compressed_identity).trim_end_matches('0').trim_end_matches('.').to_string();
        let bi_str = format!("{:.6}", block_identity).trim_end_matches('0').trim_end_matches('.').to_string();

        let cigar_str : String = cigar.iter().map(|op| format!("{}{}", op.len(), op.op())).collect();

        match name {
            Some(ref name) => println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tgi:f:{}\tbi:f:{}\tcg:Z:{}\tan:Z:{}",
                                    overlap_name, query_length, first, last, strand,
                                    target_name, target_length, overlap_target.first, overlap_target.last,
                                    matches, block_len, 255, gi_str, bi_str, cigar_str, name),
            None => println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tgi:f:{}\tbi:f:{}\tcg:Z:{}",
                                overlap_name, query_length, first, last, strand,
                                target_name, target_length, overlap_target.first, overlap_target.last,
                                matches, block_len, 255, gi_str, bi_str, cigar_str),
        }
    }
}

fn print_stats(impg: &Impg) {
    println!("Number of sequences: {}", impg.seq_index.len());
    println!("Number of overlaps: {}", impg.trees.values().map(|tree| tree.len()).sum::<usize>());
}
