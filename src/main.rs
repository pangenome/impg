use clap::Parser;
use coitrees::IntervalTree;
use impg::impg::{AdjustedInterval, Impg, SerializableImpg};
use impg::paf;
use impg::partition::partition_alignments;
use log::{info, debug, warn};
use noodles::bgzf;
use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::BufRead;
use std::io::{self, BufReader, BufWriter};
use std::num::NonZeroUsize;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use crate::paf::PafRecord;
use rayon::prelude::*;


/// Common options shared between all commands
#[derive(Parser, Debug)]
struct CommonOpts {
    /// Path to the PAF files.
    #[clap(short = 'p', long, value_parser, required = false, num_args = 1.., conflicts_with = "paf_list")]
    paf_files: Vec<String>,

    /// Path to a text file containing paths to PAF files (one per line).
    #[clap(long, value_parser, required = false, conflicts_with = "paf_files")]
    paf_list: Option<String>,

    /// Path to the IMPG index file.
    #[clap(short = 'i', long, value_parser)]
    index: Option<String>,

    /// Force the regeneration of the index, even if it already exists.
    #[clap(short = 'f', long, action)]
    force_reindex: bool,

    /// Number of threads for parallel processing.
    #[clap(short = 't', long, value_parser, default_value_t = NonZeroUsize::new(4).unwrap())]
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

        /// Window size for partitioning
        #[clap(short = 'w', long, value_parser)]
        window_size: usize,

        /// Maximum distance between regions to merge
        #[clap(short = 'd', long, value_parser, default_value_t = 100000)]
        merge_distance: i32,

        /// Maximum recursion depth for transitive overlaps (0 for no limit)
        #[clap(short = 'm', long, value_parser, default_value_t = 2)]
        max_depth: u16,

        /// Minimum region size to consider for transitive queries
        #[clap(short = 'l', long, value_parser, default_value_t = 10)]
        min_transitive_len: i32,

        /// Minimum distance between transitive ranges to consider on the same sequence
        #[clap(long, value_parser, default_value_t = 10)]
        min_distance_between_ranges: i32,

        /// Path to the file with sequence names to start with (one per line)
        #[clap(long, value_parser)]
        starting_sequences_file: Option<String>,

        #[clap(
            long,
            value_parser,
            default_value = "longest",
            help = "Selection mode for next sequence:\n\
                - \"longest\": sequence with longest single missing region\n\
                - \"total\": sequence with highest total missing regions\n\
                - \"sample[,separator]\": sample with highest total missing regions\n\
                - \"haplotype[,separator]\": haplotype highest total missing regions\n\
                The sample/haplotype modes assume PanSN naming; '#' is the default separator."
        )]
        selection_mode: String,

        /// Minimum region size for missing regions
        #[clap(long, value_parser, default_value_t = 3000)]
        min_missing_size: i32,

        /// Minimum distance from sequence start/end - closer regions will be extended to the boundaries
        #[clap(long, value_parser, default_value_t = 3000)]
        min_boundary_distance: i32,
    },
    /// Query overlaps in the alignment
    Query {
        #[clap(flatten)]
        common: CommonOpts,

        /// Target range in the format `seq_name:start-end`
        #[clap(short = 'r', long, value_parser)]
        target_range: Option<String>,

        /// Path to the BED file containing target regions
        #[clap(short = 'b', long, value_parser)]
        target_bed: Option<String>,

        /// Output format: 'auto' (BED for -r, BEDPE for -b), 'bed', 'bedpe', or 'paf'
        #[clap(short = 'o', long, value_parser, default_value = "auto")]
        output_format: String,
        
        /// Enable transitive queries
        #[clap(short = 'x', long, action)]
        transitive: bool,

        /// Enable transitive queries with breadth-first search (faster, but returns more overlapping results)
        #[clap(long, action)]
        transitive_bfs: bool,

        /// Maximum recursion depth for transitive overlaps (0 for no limit)
        #[clap(short = 'm', long, value_parser, default_value_t = 0)]
        max_depth: u16,

        /// Minimum region size to consider for transitive queries
        #[clap(short = 'l', long, value_parser, default_value_t = 0)]
        min_transitive_len: i32,

        /// Minimum distance between transitive ranges to consider on the same sequence
        #[clap(long, value_parser, default_value_t = 0)]
        min_distance_between_ranges: i32,
    },
    /// Print alignment statistics
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
            merge_distance,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            starting_sequences_file,
            selection_mode,
            min_missing_size,
            min_boundary_distance,
        } => {
            validate_selection_mode(&selection_mode)?;

            let impg = initialize_impg(&common)?;
            partition_alignments(
                &impg,
                window_size,
                starting_sequences_file.as_deref(),
                &selection_mode,
                merge_distance,
                min_missing_size,
                min_boundary_distance,
                max_depth,
                min_transitive_len,
                min_distance_between_ranges,
                common.verbose > 1,
            )?;
        }
        Args::Query {
            common,
            target_range,
            target_bed,
            output_format,
            transitive,
            transitive_bfs,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
        } => {
            validate_output_format(&output_format)?;

            let impg = initialize_impg(&common)?;

            if let Some(target_range) = target_range {
                let (target_name, target_range) = parse_target_range(&target_range)?;
                let results = perform_query(
                    &impg,
                    &target_name,
                    target_range,
                    transitive,
                    transitive_bfs,
                    max_depth,
                    min_transitive_len,
                    min_distance_between_ranges,
                );

                // Output results based on the format
                match output_format.as_str() {
                    "bedpe" => {
                        // Skip the first element (the input range) for BEDPE output
                        output_results_bedpe(&impg, results.into_iter().skip(1), None);
                    },
                    "paf" => {
                        // Skip the first element (the input range) for PAF output
                        output_results_paf(&impg, results.into_iter().skip(1), None);
                    },
                    _ => { // 'auto' or 'bed'
                        // BED format - include the first element
                        output_results_bed(&impg, results);
                    }
                }
            } else if let Some(target_bed) = target_bed {
                let targets = parse_bed_file(&target_bed)?;
                info!("Parsed {} target ranges from BED file", targets.len());
                for (target_name, target_range, name) in targets {
                    let results = perform_query(
                        &impg,
                        &target_name,
                        target_range,
                        transitive,
                        transitive_bfs,
                        max_depth,
                        min_transitive_len,
                        min_distance_between_ranges,
                    );

                    // Output results based on the format
                    match output_format.as_str() {
                        "bed" => {
                            // BED format - include the first element
                            output_results_bed(&impg, results);
                        },
                        "paf" => {
                            // Skip the first element (the input range) for PAF output
                            output_results_paf(&impg, results.into_iter().skip(1), name);
                        },
                        _ => { // 'auto' or 'bedpe'
                            // Skip the first element (the input range) for BEDPE output
                            output_results_bedpe(&impg, results.into_iter().skip(1), name);                            
                        }
                    }
                }
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --target-range or --target-bed must be provided for query subcommand",
                ));
            }
        }
        Args::Stats { common } => {
            let impg = initialize_impg(&common)?;

            print_stats(&impg);
        }
    }

    Ok(())
}

fn validate_selection_mode(mode: &str) -> io::Result<()> {
    match mode {
        "longest" | "total" => Ok(()),
        mode if mode == "sample" || mode == "haplotype" 
            || mode.starts_with("sample,") || mode.starts_with("haplotype,") => Ok(()),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Invalid selection mode. Must be 'longest', 'total', 'sample[,sep]', or 'haplotype[,sep]'."
        ))
    }
}

fn validate_output_format(mode: &str) -> io::Result<()> {
    match mode {
        "auto" | "bed" | "bedpe" | "paf" => Ok(()),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Invalid output format. Must be 'auto', 'bed', 'bedpe', or 'paf'."
        ))
    }
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

    // Resolve the list of PAF files
    let paf_files = resolve_paf_files(common)?;
    info!("Found {} PAF files", paf_files.len());

    // Load or generate index
    if common.force_reindex {
        generate_multi_index(&paf_files, common.num_threads, common.index.as_deref())
    } else {
        load_or_generate_multi_index(&paf_files, common.num_threads, common.index.as_deref())
    }
}

/// Resolve the list of PAF files from either --paf-files or --paf-list
fn resolve_paf_files(common: &CommonOpts) -> io::Result<Vec<String>> {
    if !common.paf_files.is_empty() {
        return Ok(common.paf_files.clone());
    }
    
    if let Some(paf_list_file) = &common.paf_list {
        // Read PAF files from the list file
        let file = File::open(paf_list_file)?;
        let reader = BufReader::new(file);
        let mut paf_files = Vec::new();
        
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                paf_files.push(trimmed.to_string());
            }
        }
        
        if paf_files.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("No valid PAF files found in list file: {}", paf_list_file),
            ));
        }
        
        return Ok(paf_files);
    }
    
    // Neither paf_files nor paf_list provided
    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        "Either --paf-files or --paf-list must be provided",
    ))
}

fn load_or_generate_multi_index(paf_files: &[String], num_threads: NonZeroUsize, custom_index: Option<&str>) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(paf_files, custom_index);
    if std::path::Path::new(&index_file).exists() {
        load_multi_index(paf_files, custom_index)
    } else {
        generate_multi_index(paf_files, num_threads, custom_index)
    }
}

fn load_multi_index(paf_files: &[String], custom_index: Option<&str>) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(paf_files, custom_index);
    info!("Reading IMPG index from {}", index_file);

    // Check if all PAF files are newer than the index
    let index_file_metadata = std::fs::metadata(&index_file)?;
    let index_file_ts = index_file_metadata.modified().ok();
    
    if let Some(index_ts) = index_file_ts {
        for paf_file in paf_files {
            let paf_file_metadata = std::fs::metadata(paf_file)?;
            if let Ok(paf_file_ts) = paf_file_metadata.modified() {
                if paf_file_ts > index_ts {
                    warn!("WARNING:\tPAF file {} has been modified since impg index creation.", paf_file);
                }
            }
        }
    }

    let file = File::open(index_file)?;
    let reader = BufReader::new(file);
    let serializable: SerializableImpg = bincode::deserialize_from(reader).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Failed to deserialize index: {:?}", e),
        )
    })?;
    
    Ok(Impg::from_multi_paf_and_serializable(paf_files, serializable))
}

fn generate_multi_index(paf_files: &[String], num_threads: NonZeroUsize, custom_index: Option<&str>) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(paf_files, custom_index);
    info!("No index found at {}. Creating it now.", index_file);

    let num_paf_files = paf_files.len();

    // Process PAF files in parallel using Rayon
    let records_by_file: Vec<_> = (0..paf_files.len())
        .into_par_iter()
        .map(|file_index| -> io::Result<(Vec<PafRecord>, String, u32)> {
            let paf_file = &paf_files[file_index];
            // Print file path and num
            debug!("Processing PAF file ({}/{}): {}", file_index + 1, num_paf_files, paf_file);

            let file = File::open(paf_file)?;
            let reader: Box<dyn io::Read> = if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
                Box::new(bgzf::MultithreadedReader::with_worker_count(
                    num_threads,
                    file,
                ))
            } else {
                Box::new(file)
            };
            let reader = BufReader::new(reader);
            let records = paf::parse_paf(reader).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Failed to parse PAF records from {}: {:?}", paf_file, e),
                )
            })?;
            
            Ok((records, paf_file.clone(), file_index as u32))
        })
        .collect::<Result<Vec<_>, _>>()?;  // Propagate any errors
    
    let impg = Impg::from_multi_paf_records(&records_by_file).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Failed to create index: {:?}", e),
        )
    })?;

    let serializable = impg.to_serializable();
    let file = File::create(index_file)?;
    let writer = BufWriter::new(file);
    bincode::serialize_into(writer, &serializable).map_err(|e| {
        io::Error::new(
            io::ErrorKind::Other,
            format!("Failed to serialize index: {:?}", e),
        )
    })?;

    Ok(impg)
}

fn parse_target_range(target_range: &str) -> io::Result<(String, (i32, i32))> {
    let parts: Vec<&str> = target_range.rsplitn(2, ':').collect();
    if parts.len() != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Target range format should be `seq_name:start-end`",
        ));
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
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid BED file format",
            ));
        }

        let (start, end) = parse_range(&parts[1..=2])?;
        let name = parts.get(3).map(|s| s.to_string());
        ranges.push((parts[0].to_string(), (start, end), name));
    }

    Ok(ranges)
}

fn parse_range(range_parts: &[&str]) -> io::Result<(i32, i32)> {
    if range_parts.len() != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Range format should be `start-end`",
        ));
    }

    let start = range_parts[0]
        .parse::<i32>()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid start value"))?;
    let end = range_parts[1]
        .parse::<i32>()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid end value"))?;

    if start >= end {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Start value must be less than end value",
        ));
    }

    Ok((start, end))
}

fn get_combined_index_filename(paf_files: &[String], custom_index: Option<&str>) -> String {
    if let Some(index) = custom_index {
        return index.to_string();
    }
    
    if paf_files.len() == 1 {
        format!("{}.impg", paf_files[0])
    } else {
        // For multiple files, create a hash of the sorted filenames
        
        let mut file_refs: Vec<&str> = paf_files.iter().map(|s| s.as_str()).collect();
        file_refs.sort();
        
        let mut hasher = DefaultHasher::new();
        for file in &file_refs {
            file.hash(&mut hasher);
        }
        
        format!("combined_{:016x}.impg", hasher.finish())
    }
}

fn perform_query(
    impg: &Impg,
    target_name: &str,
    target_range: (i32, i32),
    transitive: bool,
    transitive_bfs: bool,
    max_depth: u16,
    min_transitive_len: i32,
    min_distance_between_ranges: i32,
) -> Vec<AdjustedInterval> {
    let (target_start, target_end) = target_range;
    let target_id = impg
        .seq_index
        .get_id(target_name)
        .expect("Target name not found in index");
    let target_length = impg
        .seq_index
        .get_len_from_id(target_id)
        .expect("Target length not found in index");
    if target_end > target_length as i32 {
        panic!(
            "Target range end ({}) exceeds the target sequence length ({})",
            target_end, target_length
        );
    }
    if transitive {
        impg.query_transitive(
            target_id,
            target_start,
            target_end,
            None,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            true,
        )
    } else if transitive_bfs {
        impg.query_transitive_bfs(
            target_id,
            target_start,
            target_end,
            None,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            false,
        )
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

fn output_results_bedpe<I>(impg: &Impg, results: I, name: Option<String>)
where
    I: Iterator<Item = AdjustedInterval>,
{
    for (overlap_query, _, overlap_target) in results {
        let query_name = impg.seq_index.get_name(overlap_query.metadata).unwrap();
        let target_name = impg.seq_index.get_name(overlap_target.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t+",
            query_name,
            first,
            last,
            target_name,
            overlap_target.first,
            overlap_target.last,
            name.as_deref().unwrap_or("."),
            strand
        );
    }
}

fn output_results_paf<I>(impg: &Impg, results: I, name: Option<String>)
where
    I: Iterator<Item = AdjustedInterval>,
{
    for (overlap_query, cigar, overlap_target) in results {
        let query_name = impg.seq_index.get_name(overlap_query.metadata).unwrap();
        let target_name = impg.seq_index.get_name(overlap_target.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };

        let query_length = impg
            .seq_index
            .get_len_from_id(overlap_query.metadata)
            .unwrap();
        let target_length = impg
            .seq_index
            .get_len_from_id(overlap_target.metadata)
            .unwrap();

        let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len) =
            cigar.iter().fold(
                (0, 0, 0, 0, 0, 0, 0),
                |(m, mm, i, i_bp, d, d_bp, bl), op| {
                    let len = op.len();
                    match op.op() {
                        'M' => (m + len, mm, i, i_bp, d, d_bp, bl + len), // We overestimate num. of matches by assuming 'M' represents matches for simplicity
                        '=' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
                        'X' => (m, mm + len, i, i_bp, d, d_bp, bl + len),
                        'I' => (m, mm, i + 1, i_bp + len, d, d_bp, bl + len),
                        'D' => (m, mm, i, i_bp, d + 1, d_bp + len, bl + len),
                        _ => (m, mm, i, i_bp, d, d_bp, bl),
                    }
                },
            );
        let gap_compressed_identity =
            (matches as f64) / (matches + mismatches + insertions + deletions) as f64;

        let edit_distance = mismatches + inserted_bp + deleted_bp;
        let block_identity = (matches as f64) / (matches + edit_distance) as f64;

        // Format bi and gi fields without trailing zeros
        let gi_str = format!("{:.6}", gap_compressed_identity)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();
        let bi_str = format!("{:.6}", block_identity)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();

        let cigar_str: String = cigar
            .iter()
            .map(|op| format!("{}{}", op.len(), op.op()))
            .collect();

        match name {
            Some(ref name) => println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tgi:f:{}\tbi:f:{}\tcg:Z:{}\tan:Z:{}",
                                query_name, query_length, first, last, strand,
                                target_name, target_length, overlap_target.first, overlap_target.last,
                                matches, block_len, 255, gi_str, bi_str, cigar_str, name),
            None => println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tgi:f:{}\tbi:f:{}\tcg:Z:{}",
                                query_name, query_length, first, last, strand,
                                target_name, target_length, overlap_target.first, overlap_target.last,
                                matches, block_len, 255, gi_str, bi_str, cigar_str),
        }
    }
}

use rustc_hash::FxHashMap;
fn print_stats(impg: &Impg) {
    // Basic stats
    let num_sequences = impg.seq_index.len();
    let total_sequence_length: usize = (0..num_sequences as u32)
        .filter_map(|id| impg.seq_index.get_len_from_id(id))
        .sum();
    let num_overlaps = impg.trees.values().map(|tree| tree.len()).sum::<usize>();
    println!("Number of sequences: {}", num_sequences);
    println!("Total sequence length: {} bp", total_sequence_length);
    println!("Number of overlaps: {}", num_overlaps);

    // Overlap distribution stats
    let mut overlaps_per_seq: FxHashMap<u32, usize> = FxHashMap::default();
    for (&target_id, tree) in &impg.trees {
        overlaps_per_seq.insert(target_id, tree.len());
    }

    let mut entries: Vec<(u32, usize)> = overlaps_per_seq.into_iter().collect();
    entries.sort_by(|a, b| b.1.cmp(&a.1));

    if !entries.is_empty() {
        // Calculate mean and median overlaps
        let sum: usize = entries.iter().map(|(_, count)| count).sum();
        let mean = sum as f64 / entries.len() as f64;

        let median = if entries.is_empty() {
            0.0
        } else if entries.len() % 2 == 0 {
            let mid = entries.len() / 2;
            (entries[mid - 1].1 + entries[mid].1) as f64 / 2.0
        } else {
            entries[entries.len() / 2].1 as f64
        };
        println!("\nMean overlaps per sequence: {:.2}", mean);
        println!("Median overlaps per sequence: {:.2}", median);

        println!("\nTop sequences by number of overlaps:");
        for (idx, (seq_id, count)) in entries.iter().take(5).enumerate() {
            if let Some(name) = impg.seq_index.get_name(*seq_id) {
                println!("{}. {}: {} overlaps", idx + 1, name, count);
            }
        }
    }
}
