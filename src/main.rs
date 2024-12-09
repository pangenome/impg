use clap::{Parser};
use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use std::num::NonZeroUsize;
use std::collections::HashMap;
use noodles::bgzf;
use impg::impg::{Impg, SerializableImpg, AdjustedInterval, check_intervals};
use coitrees::IntervalTree;
use impg::paf;
use rayon::ThreadPoolBuilder;
use std::io::BufRead;

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
}

/// Command-line tool for querying overlaps in PAF files.
#[derive(Parser, Debug)]
#[command(author, version, about, disable_help_subcommand = true)]
enum Args {
    /// Query overlaps in the PAF file.
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
        check_intervals_bool: bool,
    },

    /// Perform partitioning on the PAF file.
    Partition {
        #[clap(flatten)]
        common: CommonOpts,

        /// Path to the FASTA index file (.fai).
        #[clap(short = 'f', long, value_parser)]
        fasta_index: String,

        /// Average window size for partitioning.
        #[clap(short = 'w', long, value_parser)]
        window_size: usize,

        /// Sample name to partition.
        #[clap(short = 'n', long, value_parser)]
        sample_name: String,
    },

    /// Print stats about the index.
    Stats {
        #[clap(flatten)]
        common: CommonOpts,
    },
}

/// Initialize thread pool and load/generate index based on common options
fn initialize_impg(common: &CommonOpts) -> io::Result<Impg> {
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

fn main() -> io::Result<()> {
    let args = Args::parse();

    match args {
        Args::Query {
            common,
            target_range,
            target_bed,
            transitive,
            output_paf,
            check_intervals_bool,
        } => {
            let impg = initialize_impg(&common)?;

            // Handle query command
            if let Some(target_range) = target_range {
                let (target_name, target_range) = parse_target_range(&target_range)?;
                let results = perform_query(&impg, &target_name, target_range, transitive);
                if check_intervals_bool {
                    let invalid_cigars = check_intervals(&impg, &results);
                    if !invalid_cigars.is_empty() {
                        for (row, error_reason) in invalid_cigars {
                            eprintln!("{}; {}", error_reason, row);
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
                    if check_intervals_bool {
                        let invalid_cigars = check_intervals(&impg, &results);
                        if !invalid_cigars.is_empty() {
                            for (row, error_reason) in invalid_cigars {
                                eprintln!("{}; {}", error_reason, row);
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
        Args::Partition {
            common,
            fasta_index,
            window_size,
            sample_name,
        } => {
            let impg = initialize_impg(&common)?;

            // Handle partition command
            run_partitioning(&impg, &fasta_index, window_size, &sample_name)?;
        },
        Args::Stats { common } => {
            let impg = initialize_impg(&common)?;

            // Print stats
            print_stats(&impg);
        }
    }

    Ok(())
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

fn parse_target_range(target_range: &str) -> io::Result<(String, (i32, i32))> {
    let parts: Vec<&str> = target_range.rsplitn(2, ':').collect();
    if parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Target range format should be `seq_name:start-end`"));
    }

    let (start, end) = parse_range(&parts[0].split('-').collect::<Vec<_>>())?;
    Ok((parts[1].to_string(), (start, end)))
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

fn load_or_generate_index(paf_file: &str, num_threads: NonZeroUsize) -> io::Result<Impg> {
    let index_file = format!("{}.impg", paf_file);
    if std::path::Path::new(&index_file).exists() {
        load_index(paf_file)
    } else {
        generate_index(paf_file, num_threads)
    }
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

fn load_index(paf_file: &str) -> io::Result<Impg> {
    let index_file = format!("{}.impg", paf_file);
    
    let paf_file_metadata = std::fs::metadata(paf_file)?;
    let index_file_metadata = std::fs::metadata(index_file.clone())?;
    if let (Ok(paf_file_ts), Ok(index_file_ts)) = (paf_file_metadata.modified(), index_file_metadata.modified()) {
        if paf_file_ts > index_file_ts
        {
            eprintln!("WARNING:\tPAF file has been modified since impg index creation.");
        }
    } else {
        eprintln!("WARNING:\tUnable to compare timestamps of PAF file and impg index file. PAF file may have been modified since impg index creation.");
    }

    let file = File::open(index_file)?;
    let reader = BufReader::new(file);
    let serializable: SerializableImpg = bincode::deserialize_from(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to deserialize index: {:?}", e)))?;
    Ok(Impg::from_paf_and_serializable(paf_file, serializable))
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

fn parse_fasta_index(fasta_index: &str) -> io::Result<HashMap<String, usize>> {
    let file = File::open(fasta_index)?;
    let reader = BufReader::new(file);
    let mut lengths = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let seq_name = parts[0].to_string();
            let seq_length = parts[1]
                .parse::<usize>()
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Invalid sequence length"))?;
            lengths.insert(seq_name, seq_length);
        }
    }
    Ok(lengths)
}

fn subtract_intervals(
    missing: &[(usize, usize)],
    partitions: &[(usize, usize)],
) -> Vec<(usize, usize)> {
    let mut result = missing.to_vec();
    for &(p_start, p_end) in partitions {
        result = result
            .into_iter()
            .flat_map(|(m_start, m_end)| {
                if p_end <= m_start || p_start >= m_end {
                    vec![(m_start, m_end)]
                } else {
                    let mut intervals = Vec::new();
                    if p_start > m_start {
                        intervals.push((m_start, p_start));
                    }
                    if p_end < m_end {
                        intervals.push((p_end, m_end));
                    }
                    intervals
                }
            })
            .collect();
    }
    result
}

fn run_partitioning(
    impg: &Impg,
    fasta_index: &str,
    window_size: usize,
    sample_name: &str,
) -> io::Result<()> {
    // Parse fasta index and get sample length
    let sample_length = parse_fasta_index(fasta_index)?
        .get(sample_name)
        .copied()
        .ok_or_else(|| io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Sample '{}' not found in FASTA index", sample_name),
        ))?;

    // Generate initial windows
    let mut windows = Vec::new();
    let mut start = 0;
    while start < sample_length {
        let end = usize::min(start + window_size, sample_length);
        windows.push((sample_name.to_string(), start, end));
        start = end;
    }

    // Initialize mask and missing regions using Vec for interval tracking
    let mut mask_intervals: Vec<(usize, usize)> = Vec::new();
    let mut missing_regions = vec![(0, sample_length)];
    let mut partition_num = 0;

    while !windows.is_empty() {
        println!("Processing new window set");

        for (chrom, start, end) in &windows {
            let region = format!("{}:{}-{}", chrom, start, end);
            println!("-- Querying region {}", region);

            // Perform transitive query
            let results = perform_query(impg, chrom, (*start as i32, *end as i32), true);

            // Extract intervals from results
            let partitions: Vec<(usize, usize)> = results
                .into_iter()
                .skip(1) // Skip the input range
                .map(|(interval, _, _)| {
                    (
                        usize::min(interval.first as usize, interval.last as usize),
                        usize::max(interval.first as usize, interval.last as usize),
                    )
                })
                .collect();

            // Apply mask to partitions
            let filtered_partitions: Vec<_> = partitions
                .into_iter()
                .filter(|&(p_start, p_end)| {
                    !mask_intervals.iter().any(|&(m_start, m_end)| {
                        p_start < m_end && p_end > m_start
                    })
                })
                .collect();

            if !filtered_partitions.is_empty() {
                println!("-- Processing partition {}", partition_num);

                // Update mask with new intervals
                mask_intervals.extend(filtered_partitions.iter().copied());
                // Merge overlapping intervals for efficiency
                mask_intervals.sort_by_key(|k| k.0);
                let mut merged: Vec<(usize, usize)> = Vec::new();
                for &interval in &mask_intervals {
                    if let Some(&mut (ref mut last_start, ref mut last_end)) = merged.last_mut() {
                        if interval.0 <= *last_end {
                            *last_end = (*last_end).max(interval.1);
                        } else {
                            merged.push(interval);
                        }
                    } else {
                        merged.push(interval);
                    }
                }
                mask_intervals = merged;

                // Update missing regions
                missing_regions = subtract_intervals(&missing_regions, &filtered_partitions);

                // Output partitions
                // for &(p_start, p_end) in &filtered_partitions {
                //     println!("{}\t{}\t{}", sample_name, p_start, p_end);
                // }

                partition_num += 1;
            }
        }

        // Check for remaining missing regions
        if missing_regions.is_empty() {
            break;
        }

        // Find longest remaining region
        let &(longest_start, longest_end) = missing_regions
            .iter()
            .max_by_key(|&&(s, e)| e - s)
            .unwrap();

        // Generate new windows from longest remaining region
        windows = Vec::new();
        let mut start = longest_start;
        while start < longest_end {
            let end = usize::min(start + window_size, longest_end);
            windows.push((sample_name.to_string(), start, end));
            start = end;
        }
    }

    Ok(())
}

fn print_stats(impg: &Impg) {
    println!("Number of sequences: {}", impg.seq_index.len());
    println!("Number of overlaps: {}", impg.trees.values().map(|tree| tree.len()).sum::<usize>());
}
