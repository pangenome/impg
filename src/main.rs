use clap::Parser;
use coitrees::{Interval, IntervalTree};
use impg::faidx::FastaIndex;
use impg::impg::{AdjustedInterval, CigarOp, Impg, SerializableImpg};
use impg::paf::{PartialPafRecord, Strand};
use impg::partition::partition_alignments;
use impg::seqidx::SequenceIndex;
use log::{debug, error, info, warn};
use noodles::bgzf;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rustc_hash::FxHashMap;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter};

use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

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
    threads: NonZeroUsize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[clap(short, long, default_value = "0")]
    verbose: u8,
}

/// Common FASTA and POA scoring options
#[derive(Parser, Debug)]
struct GfaMafFastaOpts {
    /// List of FASTA file paths (required for 'gfa', 'maf', and 'fasta')
    #[clap(long, value_parser, num_args = 1.., value_delimiter = ' ', conflicts_with_all = &["fasta_list"])]
    fasta_files: Option<Vec<String>>,

    /// Path to a text file containing paths to FASTA files (required for 'gfa', 'maf', and 'fasta')
    #[clap(long, value_parser, conflicts_with_all = &["fasta_files"])]
    fasta_list: Option<String>,

    /// POA alignment scores as match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2 (for 'gfa' and 'maf')
    #[clap(long, value_parser, default_value = "1,4,6,2,26,1")]
    poa_scoring: String,

    /// Reverse complement reverse strand sequences (for 'fasta' output)
    #[clap(long, action)]
    reverse_complement: bool,
}

impl GfaMafFastaOpts {
    /// Resolve FASTA files from either --fasta-files or --fasta-list
    fn resolve_fasta_files(self) -> io::Result<Vec<String>> {
        match (self.fasta_files, self.fasta_list) {
            // Handle --fasta-files option - no clone needed!
            (Some(files), None) => Ok(files),
            // Handle --fasta-list option
            (None, Some(list_file)) => {
                let content = std::fs::read_to_string(&list_file).map_err(|e| {
                    io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("Failed to read FASTA list file '{}': {}", list_file, e),
                    )
                })?;

                Ok(content
                    .lines()
                    .filter(|line| !line.trim().is_empty() && !line.trim().starts_with('#'))
                    .map(|line| line.trim().to_string())
                    .collect())
            }
            (None, None) => Ok(Vec::new()),
            (Some(_), Some(_)) => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Cannot specify both --fasta-files and --fasta-list",
            )),
        }
    }

    /// Parse POA scoring parameters
    fn parse_poa_scoring(&self) -> io::Result<(u8, u8, u8, u8, u8, u8)> {
        let parts: Vec<&str> = self.poa_scoring.split(',').collect();
        if parts.len() != 6 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "POA scoring format should be 'match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2'",
            ));
        }

        let parse_u8 = |s: &str, name: &str| {
            s.parse::<u8>().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Invalid {} value", name),
                )
            })
        };

        Ok((
            parse_u8(parts[0], "match score")?,
            parse_u8(parts[1], "mismatch cost")?,
            parse_u8(parts[2], "gap opening 1 cost")?,
            parse_u8(parts[3], "gap extension 1 cost")?,
            parse_u8(parts[4], "gap opening 2 cost")?,
            parse_u8(parts[5], "gap extension 2 cost")?,
        ))
    }

    /// Build FASTA index if files are provided
    fn build_fasta_index(self) -> io::Result<Option<FastaIndex>> {
        let fasta_files = self.resolve_fasta_files()?;

        if fasta_files.is_empty() {
            Ok(None)
        } else {
            match FastaIndex::build_from_files(&fasta_files) {
                Ok(index) => {
                    info!(
                        "Built FASTA index for {} files with {} sequences",
                        index.fasta_paths.len(),
                        index.path_key_to_fasta.len()
                    );
                    Ok(Some(index))
                }
                Err(e) => {
                    error!("Failed to build FASTA index: {}", e);
                    Err(e)
                }
            }
        }
    }
}

/// Common query and filtering options
#[derive(Parser, Debug, Clone)]
struct QueryOpts {
    /// Target range in the format `seq_name:start-end`
    #[clap(short = 'r', long, value_parser)]
    target_range: Option<String>,

    /// Path to the BED file containing target regions
    #[clap(short = 'b', long, value_parser)]
    target_bed: Option<String>,

    /// Maximum distance between regions to merge
    #[clap(
        short = 'd',
        long,
        value_parser,
        conflicts_with = "no_merge",
        default_value_t = 0
    )]
    merge_distance: i32,

    /// Disable merging for all output formats
    #[clap(long, action, conflicts_with = "merge_distance")]
    no_merge: bool,

    /// Minimum gap-compressed identity threshold (0.0-1.0)
    #[clap(long, value_parser)]
    min_identity: Option<f64>,

    /// Enable transitive queries (with Breadth-First Search)
    #[clap(short = 'x', long, action, conflicts_with = "transitive_dfs")]
    transitive: bool,

    /// Enable transitive queries with Depth-First Search (slower, but returns fewer overlapping results)
    #[clap(long, action, conflicts_with = "transitive")]
    transitive_dfs: bool,

    /// Maximum recursion depth for transitive overlaps (0 for no limit)
    #[clap(short = 'm', long, value_parser, default_value_t = 0)]
    max_depth: u16,

    /// Minimum region size to consider for transitive queries
    #[clap(short = 'l', long, value_parser, default_value_t = 0)]
    min_transitive_len: i32,

    /// Minimum distance between transitive ranges to consider on the same sequence
    #[clap(long, value_parser, default_value_t = 0)]
    min_distance_between_ranges: i32,
}

impl QueryOpts {
    /// Get effective merge distance (-1 if merging is disabled)
    fn effective_merge_distance(&self) -> i32 {
        if self.no_merge {
            -1
        } else {
            self.merge_distance
        }
    }
}

/// Command-line tool for querying overlaps in PAF files.
#[derive(Parser, Debug)]
#[command(author, version, about, disable_help_subcommand = true)]
enum Args {
    /// Create an IMPG index
    Index {
        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Partition the alignment
    Partition {
        #[clap(flatten)]
        common: CommonOpts,

        /// Window size for partitioning
        #[clap(short = 'w', long, value_parser)]
        window_size: usize,

        /// Output format: 'bed', 'gfa' (v1.0), 'maf', or 'fasta' ('gfa', 'maf', and 'fasta' require --fasta-files or --fasta-list)
        #[clap(short = 'o', long, value_parser, default_value = "bed")]
        output_format: String,

        #[clap(flatten)]
        gfa_maf_fasta: GfaMafFastaOpts,

        /// Maximum distance between regions to merge
        #[clap(short = 'd', long, value_parser, default_value_t = 100000)]
        merge_distance: i32,

        /// Minimum gap-compressed identity threshold (0.0-1.0)
        #[clap(long, value_parser)]
        min_identity: Option<f64>,

        /// Enable transitive queries with Depth-First Search (slower, but returns fewer overlapping results)
        #[clap(long, action)]
        transitive_dfs: bool,

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

        /// Force processing of large regions (>10kbp) with maf/gfa output formats
        #[clap(long, action)]
        force_large_region: bool,
    },
    /// Query overlaps in the alignment
    Query {
        #[clap(flatten)]
        common: CommonOpts,

        #[clap(flatten)]
        query: QueryOpts,

        /// Output format: 'auto' ('bed' for -r, 'bedpe' for -b), 'bed', 'bedpe', 'paf', 'gfa' (v1.0), 'maf', or 'fasta' ('gfa', 'maf', and 'fasta' require --fasta-files or --fasta-list)
        #[clap(short = 'o', long, value_parser, default_value = "auto")]
        output_format: String,

        #[clap(flatten)]
        gfa_maf_fasta: GfaMafFastaOpts,

        /// Force processing of large regions (>10kbp) with maf/gfa output formats
        #[clap(long, action)]
        force_large_region: bool,
    },
    /// Compute pairwise similarity between sequences in a region
    Similarity {
        #[clap(flatten)]
        common: CommonOpts,

        #[clap(flatten)]
        query: QueryOpts,

        #[clap(flatten)]
        gfa_maf_fasta: GfaMafFastaOpts,

        /// Output distances instead of similarities
        #[clap(long, action)]
        distances: bool,

        /// Emit entries for all pairs of groups, including those with zero intersection
        #[clap(short = 'a', long, action, default_value_t = false)]
        all: bool,

        /// The part of each path name before this delimiter is a group identifier
        #[clap(long, value_parser)]
        delim: Option<char>,

        /// Consider the N-th occurrence of the delimiter (1-indexed, default: 1)
        #[clap(long, value_parser, default_value_t = 1)]
        delim_pos: u16,

        /// Perform PCA/MDS dimensionality reduction on the distance matrix
        #[clap(long, action)]
        pca: bool,

        /// Number of PCA components to output (default: 2)
        #[clap(long, value_parser, requires = "pca", default_value_t = 2)]
        pca_components: usize,

        /// Number of previous regions to use for adaptive polarization (0 to disable)
        #[clap(
            long,
            value_parser,
            requires = "pca",
            conflicts_with = "polarize_guide_samples",
            default_value_t = 3
        )]
        polarize_n_prev: usize,

        /// Comma-separated names of the samples to use for adaptive polarization
        #[clap(
            long,
            value_parser,
            conflicts_with = "polarize_n_prev",
            value_delimiter = ','
        )]
        polarize_guide_samples: Option<Vec<String>>,

        /// Similarity measure to use for PCA distance matrix ("jaccard", "cosine", or "dice")
        #[clap(long, value_parser, requires = "pca", default_value = "jaccard")]
        pca_measure: String,

        /// Force processing of large regions (>10kbp)
        #[clap(long, action)]
        force_large_region: bool,
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
        Args::Index { common } => {
            let _ = initialize_impg(&common)?;

            info!("Index created successfully");
        }
        Args::Partition {
            common,
            window_size,
            output_format,
            gfa_maf_fasta,
            merge_distance,
            min_identity,
            transitive_dfs,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            starting_sequences_file,
            selection_mode,
            min_missing_size,
            min_boundary_distance,
            force_large_region,
        } => {
            validate_selection_mode(&selection_mode)?;
            validate_output_format(&output_format, &["bed", "gfa", "maf", "fasta"])?;

            validate_region_size((0, window_size as i32), &output_format, force_large_region)?;

            // Extract reverse_complement before moving gfa_maf_fasta
            let reverse_complement = gfa_maf_fasta.reverse_complement;

            // Parse POA scoring parameters if GFA/MAF output is requested
            let scoring_params = if output_format == "gfa" || output_format == "maf" {
                Some(gfa_maf_fasta.parse_poa_scoring()?)
            } else {
                None
            };

            // Build FASTA index if GFA/MAF/FASTA output is requested
            let fasta_index = if output_format == "gfa"
                || output_format == "maf"
                || output_format == "fasta"
            {
                let index = gfa_maf_fasta.build_fasta_index()?;
                if index.is_none() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("FASTA files are required for '{}' output format. Use --fasta-files or --fasta-list", output_format),
                    ));
                }
                index
            } else {
                None
            };

            let impg = initialize_impg(&common)?;

            partition_alignments(
                &impg,
                window_size,
                starting_sequences_file.as_deref(),
                &selection_mode,
                merge_distance,
                min_identity,
                min_missing_size,
                min_boundary_distance,
                transitive_dfs,
                max_depth,
                min_transitive_len,
                min_distance_between_ranges,
                &output_format,
                fasta_index.as_ref(),
                scoring_params,
                reverse_complement,
                common.verbose > 1,
            )?;
        }
        Args::Query {
            common,
            query,
            output_format,
            gfa_maf_fasta,
            force_large_region,
        } => {
            validate_output_format(
                &output_format,
                &["auto", "bed", "bedpe", "paf", "gfa", "maf", "fasta"],
            )?;

            // Early validation for target_range before ANY expensive operations
            if let Some(target_range_str) = &query.target_range {
                let (_, target_range) = parse_target_range(target_range_str)?;
                validate_region_size(target_range, &output_format, force_large_region)?;
            }

            // Early validation for target_bed and cache parsed data to avoid re-parsing
            let cached_targets = if let Some(target_bed) = &query.target_bed {
                let targets = impg::partition::parse_bed_file(target_bed)?;
                for (_, target_range, _) in &targets {
                    validate_region_size(*target_range, &output_format, force_large_region)?;
                }
                Some(targets)
            } else {
                None
            };

            // Extract reverse_complement before moving gfa_maf_fasta
            let reverse_complement = gfa_maf_fasta.reverse_complement;

            // Parse POA scoring parameters if GFA/MAF output is requested
            let scoring_params = if output_format == "gfa" || output_format == "maf" {
                Some(gfa_maf_fasta.parse_poa_scoring()?)
            } else {
                None
            };

            // Build FASTA index if GFA/MAF/FASTA output is requested
            let fasta_index = if output_format == "gfa"
                || output_format == "maf"
                || output_format == "fasta"
            {
                let index = gfa_maf_fasta.build_fasta_index()?;
                if index.is_none() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("FASTA files are required for '{}' output format. Use --fasta-files or --fasta-list", output_format),
                    ));
                }
                index
            } else {
                None
            };

            let impg = initialize_impg(&common)?;

            if let Some(target_range) = &query.target_range {
                let (target_name, target_range) = parse_target_range(target_range)?;

                let mut results = perform_query(
                    &impg,
                    &target_name,
                    target_range,
                    output_format == "paf" || output_format == "bedpe", // Store CIGAR for PAF/BEDPE output
                    query.min_identity,
                    query.transitive,
                    query.transitive_dfs,
                    query.max_depth,
                    query.min_transitive_len,
                    query.min_distance_between_ranges,
                )?;

                // Output results based on the format
                match output_format.as_str() {
                    "bedpe" => {
                        // Skip the first element (the input range) for BEDPE output
                        results.remove(0);
                        output_results_bedpe(
                            &impg,
                            &mut results,
                            None,
                            query.effective_merge_distance(),
                        );
                    }
                    "paf" => {
                        // Skip the first element (the input range) for PAF output
                        results.remove(0);
                        output_results_paf(
                            &impg,
                            &mut results,
                            None,
                            query.effective_merge_distance(),
                        );
                    }
                    "gfa" => {
                        output_results_gfa(
                            &impg,
                            &mut results,
                            &fasta_index.unwrap(),
                            None,
                            query.effective_merge_distance(),
                            scoring_params.unwrap(),
                        )?;
                    }
                    "maf" => {
                        output_results_maf(
                            &impg,
                            &mut results,
                            &fasta_index.unwrap(),
                            None,
                            query.effective_merge_distance(),
                            scoring_params.unwrap(),
                        )?;
                    }
                    "fasta" => {
                        output_results_fasta(
                            &impg,
                            &mut results,
                            &fasta_index.unwrap(),
                            None,
                            query.effective_merge_distance(),
                            reverse_complement,
                        )?;
                    }
                    _ => {
                        // 'auto' or 'bed'
                        // BED format - include the first element
                        output_results_bed(&impg, &mut results, query.effective_merge_distance());
                    }
                }
            } else if let Some(targets) = cached_targets {
                info!("Parsed {} target ranges from BED file", targets.len());
                for (target_name, target_range, name) in targets {
                    let mut results = perform_query(
                        &impg,
                        &target_name,
                        target_range,
                        output_format == "paf" || output_format == "bedpe", // Store CIGAR for PAF/BEDPE output
                        query.min_identity,
                        query.transitive,
                        query.transitive_dfs,
                        query.max_depth,
                        query.min_transitive_len,
                        query.min_distance_between_ranges,
                    )?;

                    // Output results based on the format
                    match output_format.as_str() {
                        "bed" => {
                            // BED format - include the first element
                            output_results_bed(
                                &impg,
                                &mut results,
                                query.effective_merge_distance(),
                            );
                        }
                        "paf" => {
                            // Skip the first element (the input range) for PAF output
                            results.remove(0);
                            output_results_paf(
                                &impg,
                                &mut results,
                                name,
                                query.effective_merge_distance(),
                            );
                        }
                        "gfa" => {
                            output_results_gfa(
                                &impg,
                                &mut results,
                                fasta_index.as_ref().unwrap(),
                                name,
                                query.effective_merge_distance(),
                                scoring_params.unwrap(),
                            )?;
                        }
                        "maf" => {
                            output_results_maf(
                                &impg,
                                &mut results,
                                fasta_index.as_ref().unwrap(),
                                name,
                                query.effective_merge_distance(),
                                scoring_params.unwrap(),
                            )?;
                        }
                        "fasta" => {
                            output_results_fasta(
                                &impg,
                                &mut results,
                                fasta_index.as_ref().unwrap(),
                                name,
                                query.effective_merge_distance(),
                                reverse_complement,
                            )?;
                        }
                        _ => {
                            // 'auto' or 'bedpe'
                            // Skip the first element (the input range) for BEDPE output
                            results.remove(0);
                            output_results_bedpe(
                                &impg,
                                &mut results,
                                name,
                                query.effective_merge_distance(),
                            );
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
        Args::Similarity {
            common,
            query,
            gfa_maf_fasta,
            distances,
            all,
            delim,
            delim_pos,
            pca,
            pca_components,
            pca_measure,
            polarize_n_prev,
            polarize_guide_samples,
            force_large_region,
        } => {
            // Validate delim_pos
            if delim_pos < 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "delim-pos must be greater than 0",
                ));
            }

            if pca {
                // Validate components
                if pca_components == 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Number of components must be greater than 0",
                    ));
                }

                // Validate pca_measure
                if !["jaccard", "cosine", "dice"].contains(&pca_measure.as_str()) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Invalid PCA similarity measure '{}'. Must be one of: jaccard, cosine, dice",
                            pca_measure
                        ),
                    ));
                }
            }

            // Early validation for target_range before expensive operations
            if let Some(target_range_str) = &query.target_range {
                let (_, target_range) = parse_target_range(target_range_str)?;
                validate_region_size(target_range, "gfa", force_large_region)?; // Similarity always uses SPOA
            }

            // Early validation for target_bed and cache parsed data to avoid re-parsing
            let cached_targets = if let Some(target_bed) = &query.target_bed {
                let targets = impg::partition::parse_bed_file(target_bed)?;
                for (_, target_range, _) in &targets {
                    // Similarity always uses SPOA
                    validate_region_size(*target_range, "gfa", force_large_region)?;
                }
                Some(targets)
            } else {
                None
            };

            // Parse POA scoring parameters
            let scoring_params = gfa_maf_fasta.parse_poa_scoring()?;

            // Build FASTA index (always required for similarity)
            let fasta_index = gfa_maf_fasta.build_fasta_index()?.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "FASTA files are required for similarity computation",
                )
            })?;

            let impg = initialize_impg(&common)?;

            if let Some(target_range) = &query.target_range {
                let (target_name, target_range) = parse_target_range(target_range)?;

                let mut results = perform_query(
                    &impg,
                    &target_name,
                    target_range,
                    false, // Don't need CIGAR for similarity
                    query.min_identity,
                    query.transitive,
                    query.transitive_dfs,
                    query.max_depth,
                    query.min_transitive_len,
                    query.min_distance_between_ranges,
                )?;

                // Merge intervals if needed
                merge_query_adjusted_intervals(
                    &mut results,
                    query.effective_merge_distance(),
                    true,
                );

                // Extract query intervals
                let query_intervals: Vec<Interval<u32>> = results
                    .iter()
                    .map(|(query_interval, _, _)| *query_interval)
                    .collect();
                let region = format!("{}:{}-{}", target_name, target_range.0, target_range.1);
                let query_data = vec![(query_intervals, region)];

                // Compute and output similarities
                impg::similarity::compute_and_output_similarities(
                    &impg,
                    query_data,
                    &fasta_index,
                    scoring_params,
                    distances,
                    all,
                    delim,
                    delim_pos,
                    pca,
                    pca_components,
                    &pca_measure,
                    0,    // No polarization for single query
                    None, // No polarization for single query
                )?;
            } else if let Some(targets) = cached_targets {
                info!("Parsed {} target ranges from BED file", targets.len());

                // Query all regions serially (already parallelized internally)
                let mut all_query_data = Vec::new();
                for (target_name, target_range, _name) in targets {
                    let mut results = perform_query(
                        &impg,
                        &target_name,
                        target_range,
                        false, // Don't need CIGAR for similarity
                        query.min_identity,
                        query.transitive,
                        query.transitive_dfs,
                        query.max_depth,
                        query.min_transitive_len,
                        query.min_distance_between_ranges,
                    )?;

                    // Merge intervals if needed
                    merge_query_adjusted_intervals(
                        &mut results,
                        query.effective_merge_distance(),
                        true,
                    );

                    // Extract query intervals
                    let query_intervals: Vec<Interval<u32>> = results
                        .iter()
                        .map(|(query_interval, _, _)| *query_interval)
                        .collect();

                    let region = format!("{}:{}-{}", target_name, target_range.0, target_range.1);
                    all_query_data.push((query_intervals, region));
                }

                // Process all regions in parallel
                impg::similarity::compute_and_output_similarities(
                    &impg,
                    all_query_data,
                    &fasta_index,
                    scoring_params,
                    distances,
                    all,
                    delim,
                    delim_pos,
                    pca,
                    pca_components,
                    &pca_measure,
                    polarize_n_prev,
                    polarize_guide_samples.as_deref(),
                )?;
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

fn validate_output_format(format: &str, valid_formats: &[&str]) -> io::Result<()> {
    if valid_formats.contains(&format) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Invalid output format '{}'. Must be one of: {}",
                format,
                valid_formats.join(", ")
            ),
        ))
    }
}

/// Validate region size
fn validate_region_size(
    target_range: (i32, i32),
    output_format: &str,
    force_large_region: bool,
) -> io::Result<()> {
    let (start, end) = target_range;

    // Check that start and end are non-negative
    if start < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Target range start ({}) cannot be negative", start),
        ));
    }

    if end < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Target range end ({}) cannot be negative", end),
        ));
    }

    // Check that start <= end
    if start > end {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Target range start ({}) cannot be greater than end ({})",
                start, end
            ),
        ));
    }

    let region_size = (target_range.1 - target_range.0).unsigned_abs() as u64;
    const SIZE_LIMIT: u64 = 10_000; // 10kbp limit

    // Check if this is a maf/gfa output format that uses SPOA
    let uses_spoa = matches!(output_format, "maf" | "gfa");

    if uses_spoa && region_size > SIZE_LIMIT && !force_large_region {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Region size ({} bp) exceeds 10kbp for '{}' output format, which may require large time and memory. Use --force-large-region to proceed anyway.",
                region_size, output_format
            ),
        ));
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
        .num_threads(common.threads.into())
        .build_global()
        .unwrap();

    // Resolve the list of PAF files
    let paf_files = resolve_paf_files(common)?;
    info!("Found {} PAF files", paf_files.len());

    // Load or generate index
    if common.force_reindex {
        generate_multi_index(&paf_files, common.threads, common.index.as_deref())
    } else {
        load_or_generate_multi_index(&paf_files, common.threads, common.index.as_deref())
    }
}

/// Resolve the list of PAF files from either --paf-files or --paf-list
fn resolve_paf_files(common: &CommonOpts) -> io::Result<Vec<String>> {
    let paf_files = if !common.paf_files.is_empty() {
        common.paf_files.clone()
    } else if let Some(paf_list_file) = &common.paf_list {
        // Read PAF files from the list file
        let file = File::open(paf_list_file)?;
        let reader = BufReader::new(file);
        let mut files = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                files.push(trimmed.to_string());
            }
        }

        if files.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("No valid PAF files found in list file: {}", paf_list_file),
            ));
        }

        files
    } else {
        // Neither paf_files nor paf_list provided
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Either --paf-files or --paf-list must be provided",
        ));
    };

    // Check if the number of PAF files exceeds u16::MAX
    if paf_files.len() > u16::MAX as usize {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Too many PAF files specified: {} (maximum allowed: {})",
                paf_files.len(),
                u16::MAX
            ),
        ));
    }

    Ok(paf_files)
}

fn load_or_generate_multi_index(
    paf_files: &[String],
    threads: NonZeroUsize,
    custom_index: Option<&str>,
) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(paf_files, custom_index);
    if std::path::Path::new(&index_file).exists() {
        load_multi_index(paf_files, custom_index)
    } else {
        generate_multi_index(paf_files, threads, custom_index)
    }
}

fn load_multi_index(paf_files: &[String], custom_index: Option<&str>) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(paf_files, custom_index);
    info!("Reading IMPG index from {}", index_file);

    // Check if all PAF files are newer than the index
    let index_file_metadata = std::fs::metadata(&index_file)?;
    let index_file_ts = index_file_metadata.modified().ok();

    if let Some(index_ts) = index_file_ts {
        paf_files.par_iter().for_each(|paf_file| {
            if let Ok(paf_file_metadata) = std::fs::metadata(paf_file) {
                if let Ok(paf_file_ts) = paf_file_metadata.modified() {
                    if paf_file_ts > index_ts {
                        warn!(
                            "WARNING:\tPAF file {} has been modified since impg index creation.",
                            paf_file
                        );
                    }
                }
            }
        });
    }

    let file = File::open(index_file)?;
    let mut reader = BufReader::new(file);
    let serializable: SerializableImpg =
        bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard()).map_err(
            |e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Failed to deserialize index: {:?}", e),
                )
            },
        )?;

    Ok(Impg::from_multi_paf_and_serializable(
        paf_files,
        serializable,
    ))
}

fn generate_multi_index(
    paf_files: &[String],
    threads: NonZeroUsize,
    custom_index: Option<&str>,
) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(paf_files, custom_index);
    info!("No index found at {}. Creating it now.", index_file);

    let num_paf_files = paf_files.len();
    // Thread-safe counter for tracking progress
    let files_processed = AtomicUsize::new(0);

    // Create a shared, thread-safe index
    let tmp_seq_index = Arc::new(Mutex::new(SequenceIndex::new()));

    // Process PAF files in parallel using Rayon
    let mut records_by_file: Vec<(Vec<PartialPafRecord>, String)> = (0..paf_files.len())
        .into_par_iter()
        .map(
            |file_index| -> io::Result<(Vec<PartialPafRecord>, String)> {
                let paf_file = &paf_files[file_index];

                // Increment the counter and get the new value atomically
                let current_count = files_processed.fetch_add(1, Ordering::SeqCst) + 1;
                // Print progress with sequential counter
                debug!(
                    "Processing PAF file ({}/{}): {}",
                    current_count, num_paf_files, paf_file
                );

                let file = File::open(paf_file)?;
                let reader: Box<dyn io::Read> =
                    if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
                        Box::new(bgzf::io::MultithreadedReader::with_worker_count(
                            threads, file,
                        ))
                    } else {
                        Box::new(file)
                    };
                let reader = BufReader::new(reader);

                // Lock, get IDs, build records
                let mut seq_index_guard = tmp_seq_index.lock().unwrap();
                let records = impg::paf::parse_paf(reader, &mut seq_index_guard).map_err(|e| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("Failed to parse PAF records from {}: {:?}", paf_file, e),
                    )
                })?;

                Ok((records, paf_file.clone()))
            },
        )
        .collect::<Result<Vec<_>, _>>()?; // Propagate any errors

    // Take back ownership of the SequenceIndex
    let tmp_seq_index = Arc::try_unwrap(tmp_seq_index)
        .unwrap_or_else(|_| panic!("Failed to unwrap SequenceIndex"))
        .into_inner()
        .unwrap_or_else(|_| panic!("Failed to get inner SequenceIndex"));

    // Sort sequence names to ensure deterministic order
    let mut sequence_names = tmp_seq_index
        .name_to_id
        .keys()
        .cloned()
        .collect::<Vec<String>>();
    sequence_names.par_sort_unstable(); // Order of identical sequence names is irrelevant

    // Create a deterministic SequenceIndex
    let mut seq_index = SequenceIndex::new();
    for (name, id) in &tmp_seq_index.name_to_id {
        let length = tmp_seq_index.get_len_from_id(*id).unwrap();
        seq_index.get_or_insert_id(name, Some(length));
    }

    // Update query and target IDs with the new SequenceIndex
    records_by_file.par_iter_mut().for_each(|(records, _)| {
        for record in records.iter_mut() {
            let query_name = tmp_seq_index.get_name(record.query_id).unwrap();
            record.query_id = seq_index.get_id(query_name).unwrap();

            let target_name = tmp_seq_index.get_name(record.target_id).unwrap();
            record.target_id = seq_index.get_id(target_name).unwrap();
        }
    });

    let impg = Impg::from_multi_paf_records(&records_by_file, seq_index).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Failed to create index: {:?}", e),
        )
    })?;

    let serializable = impg.to_serializable();
    let file = File::create(index_file)?;
    let mut writer = BufWriter::new(file);
    bincode::serde::encode_into_std_write(&serializable, &mut writer, bincode::config::standard())
        .map_err(|e| io::Error::other(format!("Failed to serialize index: {:?}", e)))?;

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

    let (start, end) = impg::partition::parse_range(&parts[0].split('-').collect::<Vec<_>>())?;
    Ok((parts[1].to_string(), (start, end)))
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
    store_cigar: bool,
    min_identity: Option<f64>,
    transitive: bool,
    transitive_dfs: bool,
    max_depth: u16,
    min_transitive_len: i32,
    min_distance_between_ranges: i32,
) -> io::Result<Vec<AdjustedInterval>> {
    let (target_start, target_end) = target_range;
    let target_id = impg.seq_index.get_id(target_name).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Target sequence '{}' not found in index", target_name),
        )
    })?;
    let target_length = impg.seq_index.get_len_from_id(target_id).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Target sequence '{}' length not found in index",
                target_name
            ),
        )
    })?;
    if target_end > target_length as i32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Target range end ({}) exceeds the target sequence length ({})",
                target_end, target_length
            ),
        ));
    }

    let results = if transitive {
        impg.query_transitive_bfs(
            target_id,
            target_start,
            target_end,
            None,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            store_cigar,
            min_identity,
        )
    } else if transitive_dfs {
        impg.query_transitive_dfs(
            target_id,
            target_start,
            target_end,
            None,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            store_cigar,
            min_identity,
        )
    } else {
        impg.query(
            target_id,
            target_start,
            target_end,
            store_cigar,
            min_identity,
        )
    };

    info!("Collected {} results", results.len());

    Ok(results)
}

fn output_results_bed(impg: &Impg, results: &mut Vec<AdjustedInterval>, merge_distance: i32) {
    merge_query_adjusted_intervals(results, merge_distance, false);

    for (query_interval, _, _) in results {
        let query_name = impg.seq_index.get_name(query_interval.metadata).unwrap();
        let (first, last, strand) = if query_interval.first <= query_interval.last {
            (query_interval.first, query_interval.last, '+')
        } else {
            (query_interval.last, query_interval.first, '-')
        };
        println!("{}\t{}\t{}\t.\t.\t{}", query_name, first, last, strand);
    }
}

fn output_results_bedpe(
    impg: &Impg,
    results: &mut Vec<AdjustedInterval>,
    name: Option<String>,
    merge_distance: i32,
) {
    merge_adjusted_intervals(results, merge_distance);

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

fn output_results_paf(
    impg: &Impg,
    results: &mut Vec<AdjustedInterval>,
    name: Option<String>,
    merge_distance: i32,
) {
    merge_adjusted_intervals(results, merge_distance);

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

pub fn output_results_gfa(
    impg: &Impg,
    results: &mut Vec<AdjustedInterval>,
    fasta_index: &FastaIndex,
    _name: Option<String>,
    merge_distance: i32,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<()> {
    // Merge intervals if needed
    merge_query_adjusted_intervals(results, merge_distance, true);

    // Extract query intervals by consuming results - no cloning
    let query_intervals: Vec<Interval<u32>> = results
        .drain(..)
        .map(|(query_interval, _, _)| query_interval)
        .collect();
    let gfa_output = impg::graph::generate_gfa_from_intervals(
        impg,
        &query_intervals,
        fasta_index,
        scoring_params,
    );
    print!("{}", gfa_output);

    Ok(())
}

pub fn output_results_fasta(
    impg: &Impg,
    results: &mut Vec<AdjustedInterval>,
    fasta_index: &FastaIndex,
    _name: Option<String>,
    merge_distance: i32,
    reverse_complement: bool,
) -> io::Result<()> {
    // Merge intervals if needed
    merge_query_adjusted_intervals(results, merge_distance, false);

    // Output each sequence as FASTA
    for (query_interval, _, _) in results {
        let query_name = impg.seq_index.get_name(query_interval.metadata).unwrap();

        // Determine actual start and end based on orientation
        let (start, end, strand) = if query_interval.first <= query_interval.last {
            (query_interval.first, query_interval.last, '+')
        } else {
            (query_interval.last, query_interval.first, '-')
        };

        // Fetch the sequence
        let sequence = fasta_index.fetch_sequence(query_name, start, end)?;

        // If reverse strand and reverse complementing, reverse complement the sequence
        let sequence = if strand == '-' && reverse_complement {
            impg::graph::reverse_complement(&sequence)
        } else {
            sequence
        };

        // Output FASTA format
        let header_suffix = if strand == '-' && reverse_complement {
            "/rc"
        } else {
            ""
        };
        println!(">{}:{}-{}{}", query_name, start, end, header_suffix);

        // Print sequence in lines of 80 characters
        let sequence_str = String::from_utf8_lossy(&sequence);
        for line in sequence_str.as_bytes().chunks(80) {
            println!("{}", String::from_utf8_lossy(line));
        }
    }

    Ok(())
}

pub fn output_results_maf(
    impg: &Impg,
    results: &mut Vec<AdjustedInterval>,
    fasta_index: &FastaIndex,
    _name: Option<String>,
    merge_distance: i32,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<()> {
    // Merge intervals if needed
    merge_query_adjusted_intervals(results, merge_distance, true);

    let query_intervals: Vec<Interval<u32>> = results
        .drain(..)
        .map(|(query_interval, _, _)| query_interval)
        .collect();

    let maf_output = impg::graph::generate_maf_from_intervals(
        impg,
        &query_intervals,
        fasta_index,
        scoring_params,
    );
    print!("{}", maf_output);

    Ok(())
}

// Merge adjusted intervals by ignoring the target intervals (optimized for simple genomic interval merging in BED and GFA formats)
fn merge_query_adjusted_intervals(
    results: &mut Vec<AdjustedInterval>,
    merge_distance: i32,
    merge_strands: bool,
) {
    if results.len() > 1 && (merge_distance >= 0 || merge_strands) {
        // Sort by sequence ID, strand orientation, and start position
        results.par_sort_by_key(|(query_interval, _, _)| {
            let is_forward = query_interval.first <= query_interval.last;
            let start = if is_forward {
                query_interval.first
            } else {
                query_interval.last
            };

            (
                query_interval.metadata, // First sort by sequence ID
                start,                   // Then by actual start position
                !is_forward,             // Finally by strand orientation (forward first)
            )
        });

        let mut write_idx = 0;
        for read_idx in 1..results.len() {
            let (curr_interval, _, _) = &results[write_idx];
            let (next_interval, _, _) = &results[read_idx];

            // Check if both intervals are on the same sequence and have same orientation
            let curr_is_forward = curr_interval.first <= curr_interval.last;
            let next_is_forward = next_interval.first <= next_interval.last;

            // Extract actual start/end positions based on orientation
            let (curr_start, curr_end) = if curr_is_forward {
                (curr_interval.first, curr_interval.last)
            } else {
                (curr_interval.last, curr_interval.first)
            };

            let (next_start, next_end) = if next_is_forward {
                (next_interval.first, next_interval.last)
            } else {
                (next_interval.last, next_interval.first)
            };

            // Check if they represent the same region (different strands)
            if merge_strands
                && curr_interval.metadata == next_interval.metadata
                && curr_start == next_start
                && curr_end == next_end
            {
                // Keep the forward strand version by skipping the reversed one (don't increment write_idx)
                continue;
            }

            // Only merge if same sequence, same orientation, and within merge distance (if merge_distance >= 0)
            if merge_distance < 0
                || curr_interval.metadata != next_interval.metadata
                || curr_is_forward != next_is_forward
                || next_start > curr_end + merge_distance
            {
                write_idx += 1;
                if write_idx != read_idx {
                    results.swap(write_idx, read_idx);
                }
            } else {
                // Merge while preserving orientation
                if curr_is_forward {
                    // Forward orientation
                    results[write_idx].0.first = curr_start.min(next_start);
                    results[write_idx].0.last = curr_end.max(next_end);
                } else {
                    // Reverse orientation
                    results[write_idx].0.first = curr_end.max(next_end);
                    results[write_idx].0.last = curr_start.min(next_start);
                }
            }
        }
        results.truncate(write_idx + 1);

        info!("Collected {} merged intervals", results.len());
    }
}

// Merge adjusted intervals by considering both query and target intervals and the corresponding CIGAR operations
fn merge_adjusted_intervals(results: &mut Vec<AdjustedInterval>, merge_distance: i32) {
    if results.len() > 1 && merge_distance >= 0 {
        // Sort by query ID, query position, target ID, target position
        results.par_sort_by_key(|(query_interval, _, target_interval)| {
            let query_forward = query_interval.first < query_interval.last;

            (
                query_interval.metadata, // Group by query sequence ID
                query_forward,           // Group by orientation (keep same orientations together)
                if query_forward {
                    // Use appropriate position based on orientation
                    query_interval.first // Forward: use start position
                } else {
                    query_interval.last // Reverse: use end position
                },
                target_interval.metadata, // Group by target sequence ID
                target_interval.first,    // Target always forward
            )
        });

        let num_results = results.len();

        // Drain off all entries, taking ownership
        let mut results_iter = results.drain(..);

        // Take the first as the current "in-progress" interval
        if let Some((mut current_query, mut current_cigar, mut current_target)) =
            results_iter.next()
        {
            // Create a new vector to store merged results
            let mut merged_results = Vec::with_capacity(num_results);

            // Iterate through remaining elements
            for (next_query, next_cigar, next_target) in results_iter {
                // Determine orientations
                let query_forward = current_query.first <= current_query.last;
                let next_query_forward = next_query.first <= next_query.last;

                let target_forward = current_target.first <= current_target.last;
                let next_target_forward = next_target.first <= next_target.last;
                if !target_forward || !next_target_forward {
                    panic!("Target intervals should always be in forward!");
                }

                // Check if sequences match and orientations are the same
                if current_query.metadata != next_query.metadata
                    || current_target.metadata != next_target.metadata
                    || query_forward != next_query_forward
                {
                    // Store current interval
                    merged_results.push((current_query, current_cigar, current_target));
                    // Clone the next as the new current
                    (current_query, current_cigar, current_target) =
                        (next_query, next_cigar, next_target);
                    continue;
                }

                // Check contiguity or overlap
                let (query_contiguous, target_contiguous, query_overlap, target_overlap) =
                    if query_forward {
                        let q_contig = current_query.last == next_query.first;
                        let t_contig = current_target.last == next_target.first;
                        let q_overlap = current_query.last > next_query.first;
                        let t_overlap = current_target.last > next_target.first;
                        (q_contig, t_contig, q_overlap, t_overlap)
                    } else {
                        // Reverse orientation (remember that first > last in reverse, so we swap first/last)
                        let q_contig = current_query.first == next_query.last;
                        let t_contig = current_target.first == next_target.last;
                        let q_overlap = current_query.first > next_query.last;
                        let t_overlap = current_target.first < next_target.last;
                        (q_contig, t_contig, q_overlap, t_overlap)
                    };

                // Handle perfect contiguity (existing logic)
                if query_contiguous && target_contiguous {
                    debug!(
                        "Merge contiguous! Query: current {}:{}-{}({}), next {}:{}-{}({}); Target: current {}:{}-{}({}), next {}:{}-{}({})",
                        current_query.metadata,
                        current_query.first,
                        current_query.last,
                        if query_forward { "+" } else { "-" },
                        next_query.metadata,
                        next_query.first,
                        next_query.last,
                        if next_query_forward { "+" } else { "-" },
                        current_target.metadata,
                        current_target.first,
                        current_target.last,
                        if target_forward { "+" } else { "-" },
                        next_target.metadata,
                        next_target.first,
                        next_target.last,
                        if next_target_forward { "+" } else { "-" },
                    );

                    // Merge intervals and CIGAR operations
                    if query_forward {
                        current_query.last = next_query.last;
                        current_target.last = next_target.last;
                        current_cigar.extend_from_slice(&next_cigar);
                    } else {
                        current_query.first = next_query.first;
                        current_target.first = next_target.first;

                        let mut new_cigar =
                            Vec::with_capacity(current_cigar.len() + next_cigar.len());
                        new_cigar.extend_from_slice(&next_cigar);
                        new_cigar.extend_from_slice(&current_cigar);
                        current_cigar = new_cigar;
                    }
                    merge_consecutive_cigar_ops(&mut current_cigar);
                    continue;
                }

                // Handle overlap case
                if query_overlap && target_overlap {
                    // Calculate overlap lengths
                    let (query_overlap_len, target_overlap_len) = if query_forward {
                        (
                            next_query.first - current_query.last,
                            next_target.first - current_target.last,
                        )
                    } else {
                        // Reverse orientation (remember that first > last in reverse, so we swap first/last)
                        (
                            next_query.last - current_query.first,
                            current_target.first - next_target.last,
                        )
                    };

                    // Check if overlaps are proportional (same alignment)
                    if query_overlap_len > 0 && target_overlap_len > 0 {
                        // Check if CIGAR strings are identical in the overlap region
                        let overlap_matches = check_cigar_overlap_match(
                            &current_cigar,
                            &next_cigar,
                            query_overlap_len,
                            query_forward,
                        );

                        if overlap_matches {
                            debug!(
                                "Merge overlapping! Overlap: query={}, target={}, Query: current {}:{}-{}({}), next {}:{}-{}({}); Target: current {}:{}-{}({}), next {}:{}-{}({})",
                                query_overlap_len,
                                target_overlap_len,
                                current_query.metadata,
                                current_query.first,
                                current_query.last,
                                if query_forward { "+" } else { "-" },
                                next_query.metadata,
                                next_query.first,
                                next_query.last,
                                if next_query_forward { "+" } else { "-" },
                                current_target.metadata,
                                current_target.first,
                                current_target.last,
                                if target_forward { "+" } else { "-" },
                                next_target.metadata,
                                next_target.first,
                                next_target.last,
                                if next_target_forward { "+" } else { "-" },
                            );

                            // Trim the overlap from the next interval and merge
                            let trimmed_next_cigar = trim_cigar_prefix(
                                &next_cigar,
                                query_overlap_len,
                                target_overlap_len,
                            );

                            if query_forward {
                                current_query.last = next_query.last;
                                current_target.last = next_target.last;
                                current_cigar.extend(trimmed_next_cigar);
                            } else {
                                current_query.first = next_query.first;
                                current_target.first = next_target.first;

                                let mut new_cigar = Vec::with_capacity(
                                    trimmed_next_cigar.len() + current_cigar.len(),
                                );
                                new_cigar.extend(trimmed_next_cigar);
                                new_cigar.extend_from_slice(&current_cigar);
                                current_cigar = new_cigar;
                            }
                            continue;
                        }
                    }
                }

                // Handle gaps within merge distance
                if !query_overlap && !target_overlap {
                    let (query_gap, target_gap) = if query_forward {
                        (
                            next_query.first - current_query.last,
                            next_target.first - current_target.last,
                        )
                    } else {
                        (
                            current_query.first - next_query.last,
                            current_target.first - next_target.last,
                        )
                    };

                    // Check if gaps are within merge distance and at least one gap exists
                    if query_gap >= 0
                        && target_gap >= 0
                        && (query_gap > 0 || target_gap > 0)
                        && query_gap <= merge_distance
                        && target_gap <= merge_distance
                    {
                        debug!(
                            "Merge gaps! Query gap: {}, Target gap: {}, Query: current {}:{}-{}({}), next {}:{}-{}({}); Target: current {}:{}-{}({}), next {}:{}-{}({})",
                            query_gap,
                            target_gap,
                            current_query.metadata,
                            current_query.first,
                            current_query.last,
                            if query_forward { "+" } else { "-" },
                            next_query.metadata,
                            next_query.first,
                            next_query.last,
                            if next_query_forward { "+" } else { "-" },
                            current_target.metadata,
                            current_target.first,
                            current_target.last,
                            if target_forward { "+" } else { "-" },
                            next_target.metadata,
                            next_target.first,
                            next_target.last,
                            if next_target_forward { "+" } else { "-" },
                        );

                        // Create gap-filling CIGAR operations
                        let mut gap_cigar = Vec::new();

                        if query_gap > 0 {
                            gap_cigar.push(CigarOp::new(query_gap, 'I'));
                        }
                        if target_gap > 0 {
                            gap_cigar.push(CigarOp::new(target_gap, 'D'));
                        }

                        // Merge intervals and CIGAR
                        if query_forward {
                            current_query.last = next_query.last;
                            current_target.last = next_target.last;
                            current_cigar.extend(gap_cigar);
                            current_cigar.extend_from_slice(&next_cigar);
                        } else {
                            current_query.first = next_query.first;
                            current_target.first = next_target.first;

                            let mut new_cigar = Vec::with_capacity(
                                current_cigar.len() + gap_cigar.len() + next_cigar.len(),
                            );
                            new_cigar.extend_from_slice(&next_cigar);
                            new_cigar.extend(gap_cigar);
                            new_cigar.extend_from_slice(&current_cigar);
                            current_cigar = new_cigar;
                        }
                        merge_consecutive_cigar_ops(&mut current_cigar);
                        continue;
                    }
                }

                // No merge possible - store current and move to next
                merged_results.push((current_query, current_cigar, current_target));
                (current_query, current_cigar, current_target) =
                    (next_query, next_cigar, next_target);
            }

            // Don't forget to add the last current element
            merged_results.push((current_query, current_cigar, current_target));

            // Replace original results with merged results
            *results = merged_results;

            info!("Collected {} merged intervals", results.len());
        }
    }
}

// Merge consecutive operations of the same type
fn merge_consecutive_cigar_ops(cigar: &mut Vec<CigarOp>) {
    if cigar.len() <= 1 {
        return;
    }

    let mut write_idx = 0;
    for read_idx in 1..cigar.len() {
        if cigar[write_idx].op() == cigar[read_idx].op() {
            // Same operation type - merge by adding lengths
            let combined_len = cigar[write_idx].len() + cigar[read_idx].len();
            cigar[write_idx] = CigarOp::new(combined_len, cigar[write_idx].op());
        } else {
            // Different operation types - keep separate
            write_idx += 1;
            if write_idx != read_idx {
                cigar[write_idx] = cigar[read_idx].clone();
            }
        }
    }
    cigar.truncate(write_idx + 1);
}

// Check if CIGAR strings match in the overlap region
fn check_cigar_overlap_match(
    current_cigar: &[CigarOp],
    next_cigar: &[CigarOp],
    query_overlap_len: i32,
    query_forward: bool,
) -> bool {
    // Extract the suffix of current CIGAR that corresponds to the overlap
    let current_suffix = extract_cigar_suffix(current_cigar, query_overlap_len, query_forward);

    // Extract the prefix of next CIGAR that corresponds to the overlap
    let next_prefix = extract_cigar_prefix(next_cigar, query_overlap_len, query_forward);

    // Compare if they're identical
    current_suffix == next_prefix
}

// Extract the last part of CIGAR that covers query_len bases
fn extract_cigar_suffix(cigar: &[CigarOp], query_len: i32, forward: bool) -> Vec<CigarOp> {
    let mut result = Vec::new();
    let mut remaining_query = query_len;

    // Traverse CIGAR from end to beginning
    for op in cigar.iter().rev() {
        if remaining_query <= 0 {
            break;
        }

        let query_delta = op
            .query_delta(if forward {
                Strand::Forward
            } else {
                Strand::Reverse
            })
            .abs();

        if query_delta <= remaining_query {
            // Include entire operation
            result.push(op.clone());
            remaining_query -= query_delta;
        } else if query_delta > 0 {
            // Include partial operation
            let scale = remaining_query as f32 / query_delta as f32;
            let new_len = (op.len() as f32 * scale) as i32;
            let partial_op = CigarOp::new(new_len, op.op());
            result.push(partial_op);
            remaining_query = 0;
        }
    }

    // Reverse to get correct order
    result.reverse();
    result
}

// Extract the first part of CIGAR that covers query_len bases
fn extract_cigar_prefix(cigar: &[CigarOp], query_len: i32, forward: bool) -> Vec<CigarOp> {
    let mut result = Vec::new();
    let mut remaining_query = query_len;

    for op in cigar.iter() {
        if remaining_query <= 0 {
            break;
        }

        let query_delta = op
            .query_delta(if forward {
                Strand::Forward
            } else {
                Strand::Reverse
            })
            .abs();

        if query_delta <= remaining_query {
            // Include entire operation
            result.push(op.clone());
            remaining_query -= query_delta;
        } else if query_delta > 0 {
            // Include partial operation
            let scale = remaining_query as f32 / query_delta as f32;
            let new_len = (op.len() as f32 * scale) as i32;
            let partial_op = CigarOp::new(new_len, op.op());
            result.push(partial_op);
            remaining_query = 0;
        }
    }

    result
}

// Trim the prefix of CIGAR by removing operations that cover the first query_len/target_len bases
fn trim_cigar_prefix(cigar: &[CigarOp], query_len: i32, target_len: i32) -> Vec<CigarOp> {
    let mut result = Vec::new();
    let mut query_consumed = 0;
    let mut target_consumed = 0;
    let mut start_idx = 0;

    // Find where to start after trimming
    for (idx, op) in cigar.iter().enumerate() {
        let q_delta = op.query_delta(Strand::Forward).abs();
        let t_delta = op.target_delta();

        if query_consumed + q_delta > query_len || target_consumed + t_delta > target_len {
            // This operation partially overlaps - need to trim it
            let query_remaining = query_len - query_consumed;
            let target_remaining = target_len - target_consumed;

            // Calculate how much of this operation to skip
            let skip_ratio = if q_delta > 0 && t_delta > 0 {
                (query_remaining as f32 / q_delta as f32)
                    .min(target_remaining as f32 / t_delta as f32)
            } else if q_delta > 0 {
                query_remaining as f32 / q_delta as f32
            } else if t_delta > 0 {
                target_remaining as f32 / t_delta as f32
            } else {
                0.0
            };

            let skip_len = (op.len() as f32 * skip_ratio) as i32;

            if skip_len < op.len() {
                // Create partial operation with remaining length
                let partial_op = CigarOp::new(op.len() - skip_len, op.op());
                result.push(partial_op);
            }

            // Add all remaining operations
            start_idx = idx + 1;
            break;
        }

        query_consumed += q_delta;
        target_consumed += t_delta;

        if query_consumed >= query_len && target_consumed >= target_len {
            start_idx = idx + 1;
            break;
        }
    }

    // Add all remaining operations
    result.extend_from_slice(&cigar[start_idx..]);
    result
}

fn print_stats(impg: &Impg) {
    // Basic stats
    let num_sequences = impg.seq_index.len();
    let total_sequence_length: usize = (0..num_sequences as u32)
        .into_par_iter()
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
        let sum: usize = entries.par_iter().map(|(_, count)| count).sum();
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
