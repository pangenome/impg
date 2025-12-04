use clap::Parser;
use coitrees::{Interval, IntervalTree};
use impg::alignment_record::{AlignmentFormat, AlignmentRecord, Strand};
use impg::commands::{lace, partition, refine, similarity};
use impg::impg::{AdjustedInterval, CigarOp, Impg};
use impg::impg_index::{ImpgIndex, ImpgWrapper};
use impg::multi_impg::MultiImpg;
use impg::onealn::OneAlnParser;
use impg::seqidx::SequenceIndex;
use impg::sequence_index::{SequenceIndex as SeqIndexTrait, UnifiedSequenceIndex};
use impg::subset_filter::{apply_subset_filter, load_subset_filter, SubsetFilter};
use log::{debug, error, info, warn};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rustc_hash::FxHashMap;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Basic common options shared between all commands
#[derive(Parser, Debug)]
#[command(next_help_heading = "General options")]
struct CommonOpts {
    /// Number of threads for parallel processing.
    #[clap(short = 't', long, value_parser, default_value_t = NonZeroUsize::new(4).unwrap())]
    threads: NonZeroUsize,

    /// Verbosity level (0 = error, 1 = info, 2 = debug)
    #[clap(short, long, default_value = "1")]
    verbose: u8,
}

/// Alignment file and index options for commands that work with alignments
#[derive(Parser, Debug)]
struct AlignmentOpts {
    /// Path to the alignment files (PAF or .1aln format).
    #[arg(help_heading = "Alignment input")]
    #[clap(short = 'a', long, value_parser, required = false, num_args = 1.., conflicts_with = "alignment_list")]
    alignment_files: Vec<String>,

    /// Path to a text file containing paths to alignment files (one per line, PAF or .1aln format).
    #[arg(help_heading = "Alignment input")]
    #[clap(
        long,
        value_parser,
        required = false,
        conflicts_with = "alignment_files"
    )]
    alignment_list: Option<String>,

    /// Path to the IMPG index file.
    #[arg(help_heading = "Index options")]
    #[clap(short = 'i', long, value_parser)]
    index: Option<String>,

    /// Force the regeneration of the index, even if it already exists.
    #[arg(help_heading = "Index options")]
    #[clap(short = 'f', long, action)]
    force_reindex: bool,

    /// Use per-file indexing (one .impg file per alignment file).
    /// Enables incremental updates and parallel index building.
    #[arg(help_heading = "Index options")]
    #[clap(long, action)]
    per_file_index: bool,

    /// Trace spacing for .1aln alignment files (used when converting tracepoints to CIGAR)
    #[arg(help_heading = "Alignment options")]
    #[clap(long, value_parser, default_value = "100")]
    trace_spacing: u32,
}

/// Sequence file options for commands that need FASTA/AGC files
#[derive(Parser, Debug)]
#[command(next_help_heading = "Sequence input")]
struct SequenceOpts {
    /// List of sequence file paths (FASTA or AGC) (required for 'gfa', 'maf', and 'fasta')
    #[clap(long, value_parser, num_args = 1.., value_delimiter = ' ', conflicts_with_all = &["sequence_list"])]
    sequence_files: Option<Vec<String>>,

    /// Path to a text file containing paths to sequence files (FASTA or AGC) (required for 'gfa', 'maf', and 'fasta')
    #[clap(long, value_parser, conflicts_with_all = &["sequence_files"])]
    sequence_list: Option<String>,
}

impl SequenceOpts {
    /// Resolve sequence files from either --sequence-files or --sequence-list
    fn resolve_sequence_files(&self) -> io::Result<Vec<String>> {
        match (&self.sequence_files, &self.sequence_list) {
            // Handle --sequence-files option
            (Some(files), None) => Ok(files.clone()),
            // Handle --sequence-list option
            (None, Some(list_file)) => {
                let content = std::fs::read_to_string(list_file).map_err(|e| {
                    io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("Failed to read sequence list file '{list_file}': {e}"),
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
                "Cannot specify both --sequence-files and --sequence-list",
            )),
        }
    }

    /// Build sequence index if files are provided
    fn build_sequence_index(&self) -> io::Result<Option<UnifiedSequenceIndex>> {
        let seq_files = self.resolve_sequence_files()?;

        if seq_files.is_empty() {
            Ok(None)
        } else {
            let file_type = if seq_files.iter().any(|f| f.ends_with(".agc")) {
                "AGC"
            } else {
                "FASTA"
            };
            let num_files = seq_files.len();
            info!(
                "Building {file_type} index for {num_files} file{}",
                if num_files == 1 { "" } else { "s" }
            );

            match UnifiedSequenceIndex::from_files(&seq_files) {
                Ok(index) => Ok(Some(index)),
                Err(e) => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Failed to build sequence index: {e}"),
                )),
            }
        }
    }
}

/// Common sequence and POA scoring options
#[derive(Parser, Debug)]
struct GfaMafFastaOpts {
    #[clap(flatten)]
    sequence: SequenceOpts,

    /// POA alignment scores as match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2 (for 'gfa' and 'maf')
    #[arg(help_heading = "Alignment options")]
    #[clap(long, value_parser, default_value = "1,4,6,2,26,1")]
    poa_scoring: String,

    /// Reverse complement reverse strand sequences (for 'fasta' output)
    #[arg(help_heading = "Output options")]
    #[clap(long, action)]
    reverse_complement: bool,

    /// Force processing of large regions (>10kbp) with maf/gfa output formats
    #[arg(help_heading = "Output options")]
    #[clap(long, action)]
    force_large_region: bool,
}

impl GfaMafFastaOpts {
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
                io::Error::new(io::ErrorKind::InvalidInput, format!("Invalid {name} value"))
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

    /// Helper to validate and setup POA/sequence resources for a given output format, including sequence index for PAF with original coordinates or for 1aln files
    fn setup_output_resources(
        self,
        output_format: &str,
        original_sequence_coordinates: bool,
        alignment_files: &[String],
        approximate_mode: bool,
    ) -> io::Result<(
        Option<UnifiedSequenceIndex>,
        Option<(u8, u8, u8, u8, u8, u8)>,
    )> {
        // Check if any of the alignment files are .1aln files (which require sequence data for tracepoint conversion)
        let has_onealn_files = alignment_files.iter().any(|f| f.ends_with(".1aln"));

        // In approximate mode with bed/bedpe output, .1aln files don't need sequences
        let onealn_needs_sequences = has_onealn_files
            && !(approximate_mode && (output_format == "bed" || output_format == "bedpe"));

        let needs_sequence_mandatory = matches!(
            output_format,
            "gfa" | "maf" | "fasta" | "fasta-aln" | "fasta+paf"
        ) || onealn_needs_sequences;
        let needs_sequence_optional = output_format == "paf" && original_sequence_coordinates;
        let needs_poa = matches!(output_format, "gfa" | "maf" | "fasta-aln");

        let scoring_params = if needs_poa {
            Some(self.parse_poa_scoring()?)
        } else {
            None
        };

        let sequence_index = if needs_sequence_mandatory || needs_sequence_optional {
            let index = self.sequence.build_sequence_index()?;
            if index.is_none() && needs_sequence_mandatory {
                let file_types = "FASTA/AGC";

                let msg = if has_onealn_files {
                    format!("Sequence files ({file_types}) are required for .1aln alignment files to convert tracepoints to CIGAR strings. Use --sequence-files or --sequence-list")
                } else {
                    format!("Sequence files ({file_types}) are required for '{output_format}' output format. Use --sequence-files or --sequence-list")
                };

                return Err(io::Error::new(io::ErrorKind::InvalidInput, msg));
            }
            index
        } else {
            None
        };

        Ok((sequence_index, scoring_params))
    }
}

/// Transitive query options
#[derive(Parser, Debug, Clone)]
struct TransitiveOpts {
    /// Enable transitive queries with Depth-First Search (slower, but returns fewer overlapping results)
    #[arg(help_heading = "Transitive query options")]
    #[clap(long, action, conflicts_with = "transitive")]
    transitive_dfs: bool,

    /// Maximum recursion depth for transitive overlaps (0 for no limit)
    #[arg(help_heading = "Transitive query options")]
    #[clap(short = 'm', long, value_parser, default_value_t = 2)]
    max_depth: u16,

    /// Minimum region size to consider for transitive queries (required > trace_spacing when using --approximate )
    #[arg(help_heading = "Transitive query options")]
    #[clap(long, value_parser)]
    min_transitive_len: Option<i32>,

    /// Minimum distance between transitive ranges to consider on the same sequence
    #[arg(help_heading = "Transitive query options")]
    #[clap(long, value_parser, default_value_t = 10)]
    min_distance_between_ranges: i32,
}

impl TransitiveOpts {
    /// Get the effective min_transitive_len value (default: 100 for non-approximate mode)
    fn effective_min_transitive_len(&self) -> i32 {
        self.min_transitive_len.unwrap_or(100)
    }
}

/// Common query and filtering options
#[derive(Parser, Debug, Clone)]
struct QueryOpts {
    /// Target range in the format `seq_name:start-end`
    #[arg(help_heading = "Query region")]
    #[clap(short = 'r', long, value_parser, conflicts_with = "target_bed")]
    target_range: Option<String>,

    /// Path to the BED file containing target regions
    #[arg(help_heading = "Query region")]
    #[clap(short = 'b', long, value_parser, conflicts_with = "target_range")]
    target_bed: Option<String>,

    /// Maximum distance between regions to merge
    #[arg(help_heading = "Filtering and merging")]
    #[clap(
        short = 'd',
        long,
        value_parser,
        conflicts_with = "no_merge",
        default_value_t = 0
    )]
    merge_distance: i32,

    /// Disable merging for all output formats
    #[arg(help_heading = "Filtering and merging")]
    #[clap(long, action, conflicts_with = "merge_distance")]
    no_merge: bool,

    /// Minimum gap-compressed identity threshold (0.0-1.0)
    #[arg(help_heading = "Filtering and merging")]
    #[clap(long, value_parser)]
    min_identity: Option<f64>,

    /// Minimum output length: filter results shorter than this (bp)
    #[arg(help_heading = "Filtering and merging")]
    #[clap(short = 'l', long, value_parser)]
    min_output_length: Option<i32>,

    /// Path to a file listing sequence names to include (one per line)
    #[arg(help_heading = "Filtering and merging")]
    #[clap(long, value_parser)]
    subset_sequence_list: Option<String>,

    /// Enable transitive queries (with Breadth-First Search)
    #[arg(help_heading = "Transitive queries")]
    #[clap(short = 'x', long, action, conflicts_with = "transitive_dfs")]
    transitive: bool,

    #[clap(flatten)]
    transitive_opts: TransitiveOpts,

    /// Update coordinates to original sequences when input sequences are subsequences (seq_name:start-end) for 'bed', 'bedpe', and 'paf'
    #[arg(help_heading = "Coordinate options")]
    #[clap(long, action)]
    original_sequence_coordinates: bool,

    /// Use approximate mode for faster queries with 1aln files (only bed/bedpe output)
    #[arg(help_heading = "Performance")]
    #[clap(long, action)]
    approximate: bool,

    /// Consider strandness when merging output intervals (defaults: merge strands for bed/gfa/maf; keep separate for fasta/fasta-aln)
    #[arg(help_heading = "Filtering and merging")]
    #[clap(long, action)]
    consider_strandness: bool,
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

    /// Whether merged intervals should collapse opposite strands for a given output format
    fn merge_strands_for_output(&self, output_format: &str) -> bool {
        // Default behavior per output format
        let default = match output_format {
            "fasta" | "fasta-aln" => false,
            "maf" | "gfa" | "bed" => true,
            _ => true,
        };

        if self.consider_strandness {
            // When considering strandness, keep strands separate
            false
        } else {
            default
        }
    }
}

/// Refinement-specific options
#[derive(Parser, Debug)]
#[command(next_help_heading = "Refinement options")]
struct RefineOpts {
    #[clap(flatten)]
    query: QueryOpts,

    #[clap(flatten)]
    sequence: SequenceOpts,

    /// Minimum number of bases that supporting samples must span at each region boundary
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_parser, default_value_t = 1000)]
    span_bp: i32,

    /// Maximum per-side extension explored when maximizing boundary support.
    /// Values <= 1 are treated as fractions of the locus length; values > 1 as absolute bp.
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_parser, default_value_t = 0.5)]
    max_extension: f64,

    /// PanSN aggregation mode when counting support (sample/haplotype)
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_enum)]
    pansn_mode: Option<refine::RefineSupportArg>,

    /// Step size for expanding flanks (bp)
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_parser, default_value_t = 1000)]
    extension_step: i32,

    /// Optional BED file capturing the entities that span the refined region
    #[arg(help_heading = "Output options")]
    #[clap(long, value_parser)]
    support_output: Option<String>,

    /// BED file with ranges to blacklist when counting maximum possible entities
    #[arg(help_heading = "Refinement options")]
    #[clap(long, value_parser)]
    blacklist_bed: Option<String>,
}

impl RefineOpts {
    fn validate(&self) -> io::Result<()> {
        if self.span_bp < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--span-bp must be >= 0",
            ));
        }
        if self.max_extension < 0.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--max-extension must be >= 0",
            ));
        }
        if self.extension_step <= 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--extension-step must be > 0",
            ));
        }

        if self.query.approximate {
            validate_approximate_mode_min_length(
                self.query.transitive_opts.min_transitive_len,
                self.query.transitive || self.query.transitive_opts.transitive_dfs,
            )?;
        }

        Ok(())
    }
}

/// Parse sequence name to extract subsequence coordinates and original sequence name
/// Format: `seq_name:start-end` -> (original_seq_name, start_offset)
fn parse_subsequence_coordinates(seq_name: &str) -> Option<(String, i32)> {
    // Find the last colon to handle formats like "sample#hap#chr:start-end"
    if let Some(colon_pos) = seq_name.rfind(':') {
        let (base_name, range_part) = seq_name.split_at(colon_pos);
        let range_part = &range_part[1..]; // Remove the colon

        // Check if the range part contains a dash
        if let Some(dash_pos) = range_part.find('-') {
            let (start_str, _end_str) = range_part.split_at(dash_pos);

            // Parse the start coordinate
            if let Ok(start_offset) = start_str.parse::<i32>() {
                return Some((base_name.to_string(), start_offset));
            }
        }
    }
    None
}

/// Transform coordinates from subsequence space to original sequence space
fn transform_coordinates_to_original(
    seq_name: &str,
    start: u32,
    end: u32,
    original_coordinates: bool,
) -> (String, u32, u32) {
    if !original_coordinates {
        return (seq_name.to_string(), start, end);
    }

    if let Some((original_name, offset)) = parse_subsequence_coordinates(seq_name) {
        let offset = offset as u32;
        (original_name, start + offset, end + offset)
    } else {
        (seq_name.to_string(), start, end)
    }
}

/// Get the original sequence length when using original_sequence_coordinates
fn get_original_sequence_length(
    original_seq_name: &str,
    external_seq_index: Option<&UnifiedSequenceIndex>,
) -> usize {
    // If we have an external sequence index, try to get the length from it
    if let Some(ext_index) = external_seq_index {
        match ext_index.get_sequence_length(original_seq_name) {
            Ok(length) => return length,
            Err(_) => {
                // Emit warning when sequence not found in index
                warn!(
                    "Sequence '{original_seq_name}' not found in sequence index, using 0 as length"
                );
            }
        }
    } else {
        // Emit warning when no index is provided
        warn!("No sequence index provided, using 0 as length for sequence '{original_seq_name}'");
    }

    0 // Return 0 if the sequence is not found or no index is provided
}

/// Command-line tool for querying pangenome alignment
#[derive(Parser, Debug)]
#[command(author, version, about, disable_help_subcommand = true)]
enum Args {
    /// Create an IMPG index
    Index {
        #[clap(flatten)]
        alignment: AlignmentOpts,

        #[clap(flatten)]
        sequence: SequenceOpts,

        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Lace files together
    Lace {
        /// List of input files (space-separated)
        #[clap(short = 'f', long, value_parser, num_args = 1.., value_delimiter = ' ', conflicts_with = "file_list")]
        files: Option<Vec<String>>,

        /// Text file containing input file paths (one per line)
        #[clap(short = 'l', long, value_parser, conflicts_with = "files")]
        file_list: Option<String>,

        /// Input file format (gfa, vcf, auto)
        #[clap(long, value_parser, default_value = "auto")]
        format: String,

        /// Output file path
        #[clap(short, long, value_parser)]
        output: String,

        /// Output compression format (none, gzip, bgzip, zstd, auto)
        #[clap(long, value_parser, default_value = "auto")]
        compress: String,

        /// Gap filling mode: 0=none, 1=middle gaps only, 2=all gaps (requires --sequence-files or --sequence-list for end gaps, GFA mode only)
        #[clap(long, default_value = "0")]
        fill_gaps: u8,

        /// Skip path range length validation (faster but may miss data integrity issues)
        #[clap(long, default_value = "false")]
        skip_validation: bool,

        #[clap(flatten)]
        sequence: SequenceOpts,

        /// Directory for temporary files
        #[clap(long, value_parser)]
        temp_dir: Option<String>,

        /// Reference (FASTA or AGC) file for validating contig lengths in VCF files
        #[clap(long, value_parser)]
        reference: Option<String>,

        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Partition the alignment
    Partition {
        #[clap(flatten)]
        alignment: AlignmentOpts,

        /// Window size for partitioning
        #[arg(help_heading = "Partition options")]
        #[clap(short = 'w', long, value_parser)]
        window_size: usize,

        /// Output format: 'bed', 'gfa' (v1.0), 'maf', or 'fasta' ('gfa', 'maf', and 'fasta' require --sequence-files or --sequence-list)
        #[arg(help_heading = "Output options")]
        #[clap(short = 'o', long, value_parser, default_value = "bed")]
        output_format: String,

        /// Output folder for partition files (default: current directory)
        #[arg(help_heading = "Output options")]
        #[clap(long, value_parser)]
        output_folder: Option<String>,

        #[clap(flatten)]
        gfa_maf_fasta: GfaMafFastaOpts,

        /// Maximum distance between regions to merge
        #[arg(help_heading = "Filtering and merging")]
        #[clap(short = 'd', long, value_parser, default_value_t = 100000)]
        merge_distance: i32,

        /// Minimum gap-compressed identity threshold (0.0-1.0)
        #[arg(help_heading = "Filtering and merging")]
        #[clap(long, value_parser)]
        min_identity: Option<f64>,

        #[clap(flatten)]
        transitive_opts: TransitiveOpts,

        /// Path to the file with sequence names to start with (one per line)
        #[arg(help_heading = "Partition options")]
        #[clap(long, value_parser)]
        starting_sequences_file: Option<String>,

        #[arg(help_heading = "Partition options")]
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
        #[arg(help_heading = "Partition options")]
        #[clap(long, value_parser, default_value_t = 3000)]
        min_missing_size: i32,

        /// Minimum distance from sequence start/end - closer regions will be extended to the boundaries
        #[arg(help_heading = "Partition options")]
        #[clap(long, value_parser, default_value_t = 3000)]
        min_boundary_distance: i32,

        /// Output separate files for each partition when 'bed'
        #[arg(help_heading = "Partition options")]
        #[clap(long, action)]
        separate_files: bool,

        /// Use approximate mode for faster queries with 1aln files (only bed/bedpe output)
        #[arg(help_heading = "Performance")]
        #[clap(long, action)]
        approximate: bool,

        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Query overlaps in the alignment
    Query {
        #[clap(flatten)]
        alignment: AlignmentOpts,

        #[clap(flatten)]
        query: QueryOpts,

        /// Output format: 'auto' ('bed' for -r, 'bedpe' for -b), 'bed', 'bedpe', 'paf', 'gfa' (v1.0), 'maf', 'fasta', or 'fasta+paf' ('gfa', 'maf', 'fasta', and 'fasta+paf' require --sequence-files or --sequence-list)
        #[arg(help_heading = "Output options")]
        #[clap(short = 'o', long, value_parser, default_value = "auto")]
        output_format: String,

        /// Prefix for output file (automatically appends the extension based on format)
        #[clap(short = 'O', long, value_parser, default_value = None)]
        output_prefix: Option<String>,

        #[clap(flatten)]
        gfa_maf_fasta: GfaMafFastaOpts,

        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Refine loci to maximize the number of samples that span both ends of the region
    Refine {
        #[clap(flatten)]
        alignment: AlignmentOpts,

        #[clap(flatten)]
        refine: RefineOpts,

        #[clap(flatten)]
        common: CommonOpts,
    },
    /// Compute pairwise similarity between sequences in a region
    Similarity {
        #[clap(flatten)]
        alignment: AlignmentOpts,

        #[clap(flatten)]
        query: QueryOpts,

        #[clap(flatten)]
        gfa_maf_fasta: GfaMafFastaOpts,

        /// Show the progress bar
        #[arg(help_heading = "Output options")]
        #[clap(long, action)]
        progress_bar: bool,

        /// Output distances instead of similarities
        #[arg(help_heading = "Output options")]
        #[clap(long, action)]
        distances: bool,

        /// Emit entries for all pairs of groups, including those with zero intersection
        #[arg(help_heading = "Output options")]
        #[clap(long, action, default_value_t = false)]
        all: bool,

        /// The part of each path name before this delimiter is a group identifier
        #[arg(help_heading = "Output options")]
        #[clap(long, value_parser)]
        delim: Option<char>,

        /// Consider the N-th occurrence of the delimiter (1-indexed, default: 1)
        #[arg(help_heading = "Output options")]
        #[clap(long, value_parser, default_value_t = 1)]
        delim_pos: u16,

        /// Perform PCA/MDS dimensionality reduction on the distance matrix
        #[arg(help_heading = "PCA options")]
        #[clap(long, action)]
        pca: bool,

        /// Number of PCA components to output (default: 2)
        #[arg(help_heading = "PCA options")]
        #[clap(long, value_parser, requires = "pca", default_value_t = 2)]
        pca_components: usize,

        /// Number of previous regions to use for adaptive polarization (0 to disable)
        #[arg(help_heading = "PCA options")]
        #[clap(
            long,
            value_parser,
            requires = "pca",
            conflicts_with = "polarize_guide_samples",
            default_value_t = 3
        )]
        polarize_n_prev: usize,

        /// Comma-separated names of the samples to use for adaptive polarization
        #[arg(help_heading = "PCA options")]
        #[clap(
            long,
            value_parser,
            conflicts_with = "polarize_n_prev",
            value_delimiter = ','
        )]
        polarize_guide_samples: Option<Vec<String>>,

        /// Similarity measure to use for PCA distance matrix ("jaccard", "cosine", or "dice")
        #[arg(help_heading = "PCA options")]
        #[clap(long, value_parser, requires = "pca", default_value = "jaccard")]
        pca_measure: String,

        #[clap(flatten)]
        common: CommonOpts,
    },

    /// Print alignment statistics
    Stats {
        #[clap(flatten)]
        alignment: AlignmentOpts,

        #[clap(flatten)]
        sequence: SequenceOpts,

        #[clap(flatten)]
        common: CommonOpts,
    },
}

fn main() {
    if let Err(e) = run() {
        error!("{}", e);
        std::process::exit(1);
    }
}

fn run() -> io::Result<()> {
    let args = Args::parse();

    match args {
        Args::Index {
            common,
            alignment,
            sequence,
        } => {
            initialize_threads_and_log(&common);
            let alignment_files = resolve_alignment_files(&alignment)?;
            let sequence_files = sequence.resolve_sequence_files()?;
            let _ = initialize_impg(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files,
            )?;

            info!("Index created successfully");
        }
        Args::Lace {
            common,
            sequence,
            files,
            file_list,
            format,
            output,
            compress,
            fill_gaps,
            skip_validation,
            temp_dir,
            reference,
        } => {
            initialize_threads_and_log(&common);

            // Check that at least one input is provided
            if files.is_none() && file_list.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --files or --file-list must be provided",
                ));
            }

            // Validate gap filling mode
            if fill_gaps > 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "fill_gaps must be 0, 1, or 2",
                ));
            }

            // Validate format
            let valid_formats = ["gfa", "vcf", "auto"];
            if !valid_formats.contains(&format.as_str()) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Invalid format '{}'. Must be one of: {}",
                        format,
                        valid_formats.join(", ")
                    ),
                ));
            }

            // Validate compression format
            let valid_compress = ["none", "gzip", "bgzip", "zstd", "auto"];
            if !valid_compress.contains(&compress.as_str()) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Invalid compression format '{}'. Must be one of: {}",
                        compress,
                        valid_compress.join(", ")
                    ),
                ));
            }

            // Determine the actual format (auto-detect if needed)
            let actual_format = determine_file_format(&format, &files, &file_list)?;

            if actual_format == "vcf" {
                // Build reference sequence index if provided
                let reference_index = if let Some(ref_file) = reference {
                    Some(UnifiedSequenceIndex::from_files(&[ref_file])?)
                } else {
                    None
                };

                // VCF lacing mode
                lace::run_vcf_lace(
                    files,
                    file_list,
                    &output,
                    &compress,
                    common.threads.get(),
                    common.verbose,
                    reference_index.as_ref(),
                )?;
            } else {
                // GFA lacing mode (existing functionality)
                // Build sequence index for sequence fetching (always build if sequence files provided)
                let sequence_index = sequence.build_sequence_index()?;

                lace::run_gfa_lace(
                    files,
                    file_list,
                    &output,
                    &compress,
                    fill_gaps,
                    skip_validation,
                    temp_dir,
                    sequence_index.as_ref(),
                    common.verbose,
                )?;
            }
        }
        Args::Partition {
            common,
            alignment,
            window_size,
            output_format,
            output_folder,
            gfa_maf_fasta,
            merge_distance,
            min_identity,
            transitive_opts,
            starting_sequences_file,
            selection_mode,
            min_missing_size,
            min_boundary_distance,
            separate_files,
            approximate,
        } => {
            initialize_threads_and_log(&common);

            validate_selection_mode(&selection_mode)?;
            validate_output_format(&output_format, &["bed", "gfa", "maf", "fasta", "fasta-aln"])?;

            // Validate --approximate mode compatibility
            if approximate {
                if output_format != "bed" {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "--approximate mode is only compatible with 'bed' output format for partition, not '{}'",
                            output_format
                        ),
                    ));
                }
                validate_approximate_mode_min_length(transitive_opts.min_transitive_len, true)?;
                // Partition always uses transitive queries
            }

            validate_region_size(
                0,
                window_size as i32,
                &output_format,
                merge_distance,
                gfa_maf_fasta.force_large_region,
            )?;

            // Validate single-file output compatibility
            if !separate_files && output_format != "bed" {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Single-file output is only supported for BED format. Use --separate-files for {} format.",
                        output_format.to_uppercase()
                    ),
                ));
            }

            // Resolve alignment files once (supports process substitution inputs)
            let alignment_files = resolve_alignment_files(&alignment)?;

            // Extract fields and resolve sequence files before moving gfa_maf_fasta
            let reverse_complement = gfa_maf_fasta.reverse_complement;
            let sequence_files_for_impg = gfa_maf_fasta.sequence.resolve_sequence_files()?;

            // Setup POA/sequence resources
            let (sequence_index, scoring_params) = gfa_maf_fasta.setup_output_resources(
                &output_format,
                false,
                alignment_files.as_slice(),
                approximate,
            )?;

            // Initialize impg after validation
            let impg = initialize_impg(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files_for_impg,
            )?;

            partition::partition_alignments(
                &impg,
                window_size,
                starting_sequences_file.as_deref(),
                &selection_mode,
                merge_distance,
                min_identity,
                min_missing_size,
                min_boundary_distance,
                transitive_opts.transitive_dfs,
                transitive_opts.max_depth,
                transitive_opts.effective_min_transitive_len(),
                transitive_opts.min_distance_between_ranges,
                &output_format,
                output_folder.as_deref(),
                sequence_index.as_ref(),
                scoring_params,
                reverse_complement,
                common.verbose > 1,
                separate_files,
                approximate,
            )?;
        }
        Args::Query {
            common,
            alignment,
            query,
            output_format,
            output_prefix,
            gfa_maf_fasta,
        } => {
            initialize_threads_and_log(&common);

            validate_output_format(
                &output_format,
                &[
                    "auto",
                    "bed",
                    "bedpe",
                    "paf",
                    "gfa",
                    "maf",
                    "fasta",
                    "fasta+paf",
                    "fasta-aln",
                ],
            )?;

            // Check that either --target-range or --target-bed is provided (cheap validation)
            if query.target_range.is_none() && query.target_bed.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --target-range or --target-bed must be provided",
                ));
            }

            // Resolve alignment files once (handles process substitution inputs)
            let alignment_files = resolve_alignment_files(&alignment)?;

            // Early validation for approximate mode (before expensive operations)
            if query.approximate {
                // Check that all alignment files are .1aln files
                let non_onealn_files: Vec<&String> = alignment_files
                    .iter()
                    .filter(|f| !f.ends_with(".1aln"))
                    .collect();

                if !non_onealn_files.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "--approximate mode only works with .1aln alignment files (found {} non-.1aln files).",
                            non_onealn_files.len()
                        ),
                    ));
                }

                validate_approximate_mode_min_length(
                    query.transitive_opts.min_transitive_len,
                    query.transitive || query.transitive_opts.transitive_dfs,
                )?;
            }

            // Extract sequence files before consuming gfa_maf_fasta
            let sequence_files_for_impg = gfa_maf_fasta.sequence.resolve_sequence_files()?;

            // Load subset filter if provided
            let subset_filter = load_subset_filter_if_provided(&query.subset_sequence_list)?;

            // Initialize impg after validation but before target range validation (which needs seq_index)
            let impg = initialize_impg(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files_for_impg,
            )?;

            // Parse and validate all target ranges, tracking which parameter was used
            let (target_ranges, from_range_param) = if let Some(target_range_str) =
                &query.target_range
            {
                let (target_name, target_range, name) = if target_range_str.contains(':') {
                    partition::parse_target_range(target_range_str)?
                } else {
                    // No interval specified: use the whole sequence [0, len)
                    let seq_name = target_range_str;
                    let seq_id = impg.seq_index().get_id(seq_name).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::NotFound,
                            format!("Sequence '{seq_name}' not found in index"),
                        )
                    })?;
                    let seq_len = impg.seq_index().get_len_from_id(seq_id).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Could not get length for sequence '{seq_name}'"),
                        )
                    })? as i32;
                    let name = format!("{}:{}-{}", seq_name, 0, seq_len);
                    (seq_name.to_string(), (0, seq_len), name)
                };
                // Validate sequence exists and range is within bounds
                validate_sequence_range(
                    &target_name,
                    target_range.0,
                    target_range.1,
                    impg.seq_index(),
                )?;
                validate_range_min_length(
                    target_range.0,
                    target_range.1,
                    &name,
                    query.transitive_opts.effective_min_transitive_len(),
                )?;
                validate_region_size(
                    target_range.0,
                    target_range.1,
                    &output_format,
                    query.effective_merge_distance(),
                    gfa_maf_fasta.force_large_region,
                )?;
                (vec![(target_name, target_range, name)], true)
            } else if let Some(target_bed) = &query.target_bed {
                let targets = partition::parse_bed_file(target_bed)?;
                // Validate all entries in the BED file
                for (seq_name, (start, end), name) in &targets {
                    validate_sequence_range(seq_name, *start, *end, impg.seq_index())?;
                    validate_range_min_length(
                        *start,
                        *end,
                        name,
                        query.transitive_opts.effective_min_transitive_len(),
                    )?;
                    validate_region_size(
                        *start,
                        *end,
                        &output_format,
                        query.effective_merge_distance(),
                        gfa_maf_fasta.force_large_region,
                    )?;
                }
                (targets, false)
            } else {
                unreachable!("Already validated that either target_range or target_bed is present");
            };

            // Resolve output format based on 'auto' and parameter used
            let resolved_output_format = if output_format == "auto" {
                if from_range_param {
                    "bed"
                } else {
                    "bedpe"
                }
            } else {
                output_format.as_str()
            };

            // Validate --approximate mode output format compatibility
            if query.approximate
                && resolved_output_format != "bed"
                && resolved_output_format != "bedpe"
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "--approximate mode is only compatible with 'bed' and 'bedpe' output formats, not '{}'",
                        resolved_output_format
                    ),
                ));
            }

            // Extract reverse_complement before moving gfa_maf_fasta
            let reverse_complement = gfa_maf_fasta.reverse_complement;

            // Setup POA/sequence resources
            let (sequence_index, scoring_params) = gfa_maf_fasta.setup_output_resources(
                resolved_output_format,
                query.original_sequence_coordinates,
                alignment_files.as_slice(),
                query.approximate,
            )?;

            // Process all target ranges in a unified loop
            info!("Querying target ranges");
            for (target_name, target_range, name) in target_ranges {
                let mut results = perform_query(
                    &impg,
                    &target_name,
                    target_range,
                    resolved_output_format == "paf" || resolved_output_format == "bedpe", // Store CIGAR for PAF/BEDPE output
                    query.min_identity,
                    query.min_output_length,
                    query.transitive,
                    query.transitive_opts.transitive_dfs,
                    &query.transitive_opts,
                    sequence_index.as_ref(),
                    query.approximate,
                    subset_filter.as_ref(),
                )?;

                // Output results based on the resolved format
                match resolved_output_format {
                    "bed" => {
                        // BED format - include the first element
                        output_results_bed(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "bed")?,
                            &name,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("bed"),
                            query.original_sequence_coordinates,
                        )?;
                    }
                    "bedpe" => {
                        // Skip the first element (the input range) for BEDPE output
                        results.remove(0);
                        output_results_bedpe(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "bed")?,
                            &name,
                            query.effective_merge_distance(),
                            query.original_sequence_coordinates,
                        )?;
                    }
                    "paf" => {
                        // Skip the first element (the input range) for PAF output
                        results.remove(0);
                        output_results_paf(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "paf")?,
                            &name,
                            query.effective_merge_distance(),
                            query.original_sequence_coordinates,
                            sequence_index.as_ref(),
                        )?;
                    }
                    "gfa" => {
                        output_results_gfa(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "gfa")?,
                            sequence_index.as_ref().unwrap(),
                            &name,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("gfa"),
                            scoring_params.unwrap(),
                        )?;
                    }
                    "maf" => {
                        output_results_maf(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "maf")?,
                            sequence_index.as_ref().unwrap(),
                            &name,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("maf"),
                            scoring_params.unwrap(),
                        )?;
                    }
                    "fasta" => {
                        output_results_fasta(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "fa")?,
                            sequence_index.as_ref().unwrap(),
                            &name,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("fasta"),
                            reverse_complement,
                        )?;
                    }
                    "fasta+paf" => {
                        output_results_fasta(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "fa")?,
                            sequence_index.as_ref().unwrap(),
                            &name,
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("fasta"),
                            reverse_complement,
                        )?;
                        // Skip the first element (the input range) for PAF output
                        results.remove(0);
                        output_results_paf(
                            &impg,
                            &mut results,
                            &mut find_output_stream(&output_prefix, "paf")?,
                            &name,
                            query.effective_merge_distance(),
                            query.original_sequence_coordinates,
                            sequence_index.as_ref(),
                        )?;
                    }
                    "fasta-aln" => {
                        output_results_fasta_aln(
                            &impg,
                            &mut results,
                            sequence_index.as_ref().unwrap(),
                            name.clone(),
                            query.effective_merge_distance(),
                            query.merge_strands_for_output("fasta-aln"),
                            scoring_params.unwrap(),
                        )?;
                    }
                    _ => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("Invalid output format: {resolved_output_format}"),
                        ));
                    }
                }
            }
        }
        Args::Refine {
            alignment,
            refine,
            common,
        } => {
            initialize_threads_and_log(&common);
            refine.validate()?;

            let alignment_files = resolve_alignment_files(&alignment)?;
            let sequence_files = refine.sequence.resolve_sequence_files()?;

            // Check if we have .1aln files and validate sequence files are provided
            let has_onealn_files = alignment_files.iter().any(|f| f.ends_with(".1aln"));
            if has_onealn_files && sequence_files.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Sequence files (FASTA or AGC) are required for the 'refine' command with .1aln alignment files to convert tracepoints to CIGAR strings. Use --sequence-files or --sequence-list".to_string(),
                ));
            }

            let impg = initialize_impg(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files,
            )?;
            let subset_filter = load_subset_filter_if_provided(&refine.query.subset_sequence_list)?;

            let target_ranges = if let Some(target_range_str) = &refine.query.target_range {
                let (target_name, target_range, name) =
                    partition::parse_target_range(target_range_str)?;
                validate_sequence_range(
                    &target_name,
                    target_range.0,
                    target_range.1,
                    impg.seq_index(),
                )?;
                validate_range_min_length(
                    target_range.0,
                    target_range.1,
                    &name,
                    refine.query.transitive_opts.effective_min_transitive_len(),
                )?;
                vec![(target_name, target_range, name)]
            } else if let Some(target_bed) = &refine.query.target_bed {
                let targets = partition::parse_bed_file(target_bed)?;
                let mut validated = Vec::with_capacity(targets.len());
                for (seq_name, (start, end), name) in targets {
                    validate_sequence_range(&seq_name, start, end, impg.seq_index())?;
                    validate_range_min_length(
                        start,
                        end,
                        &name,
                        refine.query.transitive_opts.effective_min_transitive_len(),
                    )?;
                    validated.push((seq_name, (start, end), name));
                }
                validated
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --target-range or --target-bed must be provided",
                ));
            };

            // Parse blacklist BED file if provided
            let blacklist = if let Some(ref blacklist_path) = refine.blacklist_bed {
                Some(refine::parse_blacklist_bed(blacklist_path)?)
            } else {
                None
            };

            let config = refine::RefineConfig {
                span_bp: refine.span_bp,
                max_extension: refine.max_extension,
                extension_step: refine.extension_step,
                support_mode: refine
                    .pansn_mode
                    .map(Into::into)
                    .unwrap_or(refine::SupportMode::Sequence),
                merge_distance: refine.query.effective_merge_distance(),
                min_identity: refine.query.min_identity,
                use_transitive_bfs: refine.query.transitive,
                use_transitive_dfs: refine.query.transitive_opts.transitive_dfs,
                max_transitive_depth: refine.query.transitive_opts.max_depth,
                min_transitive_len: refine.query.transitive_opts.effective_min_transitive_len(),
                min_distance_between_ranges: refine
                    .query
                    .transitive_opts
                    .min_distance_between_ranges,
                subset_filter: subset_filter.as_ref(),
                blacklist: blacklist.as_ref(),
                approximate_mode: refine.query.approximate,
            };

            let mut records = refine::run_refine(&impg, &target_ranges, config)?;
            info!(
                "Refining {} targets with max_extension={} (mode: {:?})",
                target_ranges.len(),
                refine.max_extension,
                refine
                    .pansn_mode
                    .map(Into::into)
                    .unwrap_or(refine::SupportMode::Sequence)
            );
            let mut writer = BufWriter::new(io::stdout());
            let mut support_writer = if let Some(path) = &refine.support_output {
                let support_path = Path::new(path);
                if let Some(parent) = support_path.parent() {
                    if !parent.as_os_str().is_empty() {
                        std::fs::create_dir_all(parent)?;
                    }
                }
                Some(BufWriter::new(File::create(support_path)?))
            } else {
                None
            };

            // Write header
            writeln!(
                writer,
                "#chrom\tstart\tend\tname\toriginal.support\tnew.support\tleft.extension.bp\tright.extension.bp"
            )?;

            for record in records.drain(..) {
                let original_range = format!(
                    "{}:{}-{}",
                    record.chrom, record.original_start, record.original_end
                );
                let mut name_field = record.label.clone();
                if name_field.trim().is_empty() || name_field == "." {
                    name_field = original_range.clone();
                }

                // Emit an informative BED-like row: chrom start end name new_support old_support left_extension right_extension.
                // Maximizing the sample count while minimizing the expansion helps avoid loci that start or end inside SVs.
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    record.chrom,
                    record.refined_start,
                    record.refined_end,
                    name_field,
                    record.original_support_count,
                    record.support_count,
                    record.applied_left_extension,
                    record.applied_right_extension
                )?;

                if let Some(ref mut support_out) = support_writer {
                    for entity in &record.support_entities {
                        writeln!(
                            support_out,
                            "{}\t{}\t{}\t{}",
                            entity.sequence, entity.start, entity.end, name_field
                        )?;
                    }
                }
            }
            writer.flush()?;
            if let Some(mut support_out) = support_writer {
                support_out.flush()?;
            }
        }
        Args::Similarity {
            common,
            alignment,
            query,
            gfa_maf_fasta,
            progress_bar,
            distances,
            all,
            delim,
            delim_pos,
            pca,
            pca_components,
            pca_measure,
            polarize_n_prev,
            polarize_guide_samples,
        } => {
            initialize_threads_and_log(&common);

            // Validate delim_pos
            if delim_pos < 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--delim-pos must be greater than 0",
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
                            "Invalid PCA similarity measure '{pca_measure}'. Must be one of: jaccard, cosine, dice"
                        ),
                    ));
                }
            }

            // Check that either --target-range or --target-bed is provided
            if query.target_range.is_none() && query.target_bed.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Either --target-range or --target-bed must be provided",
                ));
            }

            // Resolve alignment files once (supports process substitution inputs)
            let alignment_files = resolve_alignment_files(&alignment)?;

            // Extract fields and resolve sequence files before moving gfa_maf_fasta
            let force_large_region = gfa_maf_fasta.force_large_region;
            let sequence_files_for_impg = gfa_maf_fasta.sequence.resolve_sequence_files()?;

            // Setup POA/sequence resources (always required for similarity)
            let (sequence_index, scoring_params) = gfa_maf_fasta.setup_output_resources(
                "gfa",
                false,
                alignment_files.as_slice(),
                false,
            )?;
            let sequence_index = sequence_index.unwrap(); // Safe since "gfa" always requires sequence files
            let scoring_params = scoring_params.unwrap(); // Safe since "gfa" always requires POA

            let subset_filter = load_subset_filter_if_provided(&query.subset_sequence_list)?;

            // Initialize impg after validation but before target range validation (which needs seq_index)
            let impg = initialize_impg(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files_for_impg,
            )?;

            // Validate target_range and target_bed before any expensive operations
            let target_ranges = {
                let mut targets = Vec::new();

                if let Some(target_range_str) = &query.target_range {
                    let (target_name, target_range, name) = if target_range_str.contains(':') {
                        partition::parse_target_range(target_range_str)?
                    } else {
                        // No interval specified: use the whole sequence [0, len)
                        let seq_name = target_range_str;
                        let seq_id = impg.seq_index().get_id(seq_name).ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::NotFound,
                                format!("Sequence '{seq_name}' not found in index"),
                            )
                        })?;
                        let seq_len = impg.seq_index().get_len_from_id(seq_id).ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("Could not get length for sequence '{seq_name}'"),
                            )
                        })? as i32;
                        let name = format!("{}:{}-{}", seq_name, 0, seq_len);
                        (seq_name.to_string(), (0, seq_len), name)
                    };
                    validate_sequence_range(
                        &target_name,
                        target_range.0,
                        target_range.1,
                        impg.seq_index(),
                    )?;
                    validate_region_size(
                        target_range.0,
                        target_range.1,
                        "gfa",
                        query.effective_merge_distance(),
                        force_large_region,
                    )?;
                    targets.push((target_name, target_range, name));
                }

                if let Some(target_bed) = &query.target_bed {
                    let bed_targets = partition::parse_bed_file(target_bed)?;
                    for (target_name, target_range, name) in bed_targets {
                        validate_sequence_range(
                            &target_name,
                            target_range.0,
                            target_range.1,
                            impg.seq_index(),
                        )?;
                        validate_region_size(
                            target_range.0,
                            target_range.1,
                            "gfa",
                            query.effective_merge_distance(),
                            force_large_region,
                        )?;
                        targets.push((target_name, target_range, name));
                    }
                }

                targets
            };
            // Note: Already validated that either target_range or target_bed is present

            info!("Parsed {} target ranges from BED file", target_ranges.len());

            // Query all regions serially (already parallelized internally)
            let mut all_query_data = Vec::new();
            for (target_name, target_range, _name) in target_ranges.into_iter() {
                let mut results = perform_query(
                    &impg,
                    &target_name,
                    target_range,
                    false, // Don't need CIGAR for similarity
                    query.min_identity,
                    query.min_output_length,
                    query.transitive,
                    query.transitive_opts.transitive_dfs,
                    &query.transitive_opts,
                    Some(&sequence_index),
                    query.approximate,
                    subset_filter.as_ref(),
                )?;

                let region_label = format!("{}:{}-{}", target_name, target_range.0, target_range.1);

                // Merge intervals if needed
                merge_query_adjusted_intervals(
                    &mut results,
                    query.effective_merge_distance(),
                    query.merge_strands_for_output("similarity"),
                );

                // Extract query intervals
                let query_intervals: Vec<Interval<u32>> = results
                    .iter()
                    .map(|(query_interval, _, _)| *query_interval)
                    .collect();

                all_query_data.push((query_intervals, region_label));
            }

            // Process all regions in parallel
            similarity::compute_and_output_similarities(
                &impg,
                all_query_data,
                &sequence_index,
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
                progress_bar,
            )?;
        }
        Args::Stats {
            common,
            alignment,
            sequence,
        } => {
            initialize_threads_and_log(&common);
            let alignment_files = resolve_alignment_files(&alignment)?;
            let sequence_files = sequence.resolve_sequence_files()?;
            let impg = initialize_impg(
                &common,
                &alignment,
                alignment_files.as_slice(),
                sequence_files,
            )?;

            print_stats(&impg);
        }
    }

    Ok(())
}

fn validate_approximate_mode_min_length(
    min_transitive_len_opt: Option<i32>,
    is_transitive: bool,
) -> io::Result<()> {
    // Only require explicit --min-transitive-len for transitive queries in approximate mode
    if is_transitive && min_transitive_len_opt.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--approximate mode with transitive queries requires explicitly setting --min-transitive-len. \
            Set it to greater than your trace_spacing (e.g., 101 for trace_spacing=100).",
        ));
    }

    Ok(())
}

/// Validate that a range meets the minimum length requirement
fn validate_range_min_length(
    start: i32,
    end: i32,
    range_name: &str,
    min_transitive_len: i32,
) -> io::Result<()> {
    let length = end - start;
    if length < min_transitive_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Range '{range_name}' ({length} bp) is below minimum of {min_transitive_len} bp. Lower --min-transitive-len or use a longer range"
            ),
        ));
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

/// Validates that a sequence name exists in the index and the range is within bounds
fn validate_sequence_range(
    seq_name: &str,
    start: i32,
    end: i32,
    seq_index: &SequenceIndex,
) -> io::Result<()> {
    // Check if sequence exists in index
    let seq_id = seq_index.get_id(seq_name).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!("Sequence '{seq_name}' not found in index"),
        )
    })?;

    // Get sequence length
    let seq_len = seq_index.get_len_from_id(seq_id).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Could not get length for sequence '{seq_name}'"),
        )
    })? as i32;

    // Validate range bounds
    if start < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Start position {start} cannot be negative"),
        ));
    }

    if end < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("End position {end} cannot be negative"),
        ));
    }

    if start >= end {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Start position {start} must be less than end position {end}"),
        ));
    }

    if end > seq_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "End position {end} exceeds sequence length {seq_len} for sequence '{seq_name}'"
            ),
        ));
    }

    Ok(())
}

/// Validate region size
fn validate_region_size(
    start: i32,
    end: i32,
    output_format: &str,
    merge_distance: i32,
    force_large_region: bool,
) -> io::Result<()> {
    let region_size = (end - start).unsigned_abs() as u64;
    const SIZE_LIMIT: u64 = 10_000; // 10kbp limit
    const MERGE_DISTANCE_LIMIT: i32 = 1000; // 1k limit

    // Check if this is a maf/gfa output format that uses SPOA
    let uses_spoa = matches!(output_format, "maf" | "gfa" | "fasta-aln");

    if uses_spoa && !force_large_region {
        if region_size > SIZE_LIMIT {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Region size ({region_size} bp) exceeds 10kbp for '{output_format}' output format, which may require large time and memory. Use --force-large-region to proceed anyway."
                ),
            ));
        }

        if merge_distance > MERGE_DISTANCE_LIMIT {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Merge distance ({merge_distance} bp) exceeds 1kbp for '{output_format}' output format, which may require large time and memory. Use --force-large-region to proceed anyway."
                ),
            ));
        }
    }

    Ok(())
}

/// Initialize thread pool and logger with the specified number of threads and verbosity
fn initialize_threads_and_log(common: &CommonOpts) {
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
}

/// Determine file format from explicit format or auto-detection
fn determine_file_format(
    format: &str,
    files: &Option<Vec<String>>,
    file_list: &Option<String>,
) -> io::Result<String> {
    if format != "auto" {
        return Ok(format.to_string());
    }

    // Auto-detect from first file
    let first_file = match (files, file_list) {
        (Some(files), _) if !files.is_empty() => files[0].clone(),
        (_, Some(list_file)) => {
            // Read first line from list file
            let content = std::fs::read_to_string(list_file)?;
            let first_line = content
                .lines()
                .find(|line| !line.trim().is_empty() && !line.trim().starts_with('#'))
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "No valid files found in list file",
                    )
                })?;
            first_line.trim().to_string()
        }
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "No input files provided",
            ))
        }
    };

    // Auto-detect based on file extension
    if first_file.ends_with(".vcf")
        || first_file.ends_with(".vcf.gz")
        || first_file.ends_with(".vcf.zst")
        || first_file.ends_with(".vcf.bgz")
    {
        Ok("vcf".to_string())
    } else if first_file.ends_with(".gfa")
        || first_file.ends_with(".gfa.gz")
        || first_file.ends_with(".gfa.zst")
        || first_file.ends_with(".gfa.bgz")
    {
        Ok("gfa".to_string())
    } else {
        // Try to detect by reading first few lines through decompression
        let reader = get_auto_reader(&first_file)?;

        for line in reader.lines().take(10) {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            // Check for VCF header
            if line.starts_with("#CHROM") {
                return Ok("vcf".to_string());
            }

            // Check for GFA segments, links, or paths
            if line.starts_with('S')
                || line.starts_with('L')
                || line.starts_with('P')
                || line.starts_with('H')
            {
                return Ok("gfa".to_string());
            }
        }

        // Default to GFA if we can't determine
        warn!("Could not auto-detect file format for '{first_file}', defaulting to GFA");
        Ok("gfa".to_string())
    }
}

// Helper function for auto-detection with compression support
fn get_auto_reader(path: &str) -> io::Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(path)?;
    let (reader, _format) = niffler::get_reader(Box::new(file))
        .map_err(|e| io::Error::other(format!("Failed to open reader: {e}")))?;
    Ok(Box::new(BufReader::new(reader)))
}

/// Helper function to return a Write implementer that is either standard output or a file with the
/// appropriate basename and extension. When no basename is provided, uses standard output.
fn find_output_stream(basename: &Option<String>, extension: &str) -> io::Result<Box<dyn Write>> {
    match basename {
        Some(name) => {
            let filename = format!("{}.{}", name, extension);
            let file = File::create(filename)?;
            Ok(Box::new(BufWriter::new(file)))
        }
        None => Ok(Box::new(BufWriter::new(io::stdout()))),
    }
}

/// Load/generate index based on common and alignment options
fn initialize_impg(
    common: &CommonOpts,
    alignment: &AlignmentOpts,
    alignment_files: &[String],
    sequence_files: Vec<String>,
) -> io::Result<ImpgWrapper> {
    // The list of alignment files (PAF or .1aln) is pre-resolved
    info!("Found {} alignment file(s)", alignment_files.len());

    let seq_files_opt = if sequence_files.is_empty() {
        None
    } else {
        Some(sequence_files.as_slice())
    };

    // Check if per-file indexing is requested
    if alignment.per_file_index {
        // Per-file indexing mode: each alignment file gets its own .impg index
        return initialize_multi_impg(
            alignment_files,
            common.threads,
            alignment.force_reindex,
            seq_files_opt,
        );
    }

    // Load or generate combined index (original behavior)
    let impg = if alignment.force_reindex {
        generate_multi_index(
            alignment_files,
            common.threads,
            alignment.index.as_deref(),
            seq_files_opt,
        )
    } else {
        load_or_generate_multi_index(
            alignment_files,
            common.threads,
            alignment.index.as_deref(),
            seq_files_opt,
        )
    }?;

    Ok(ImpgWrapper::from_single(impg))
}

/// Initialize MultiImpg from per-file indices
fn initialize_multi_impg(
    alignment_files: &[String],
    threads: NonZeroUsize,
    force_reindex: bool,
    sequence_files: Option<&[String]>,
) -> io::Result<ImpgWrapper> {
    use std::path::PathBuf;

    info!("Using per-file indexing mode with {} alignment files", alignment_files.len());

    // Generate per-file index paths
    let index_paths: Vec<PathBuf> = alignment_files
        .iter()
        .map(|f| PathBuf::from(format!("{}.impg", f)))
        .collect();

    // Check which indices need to be built
    let indices_to_build: Vec<usize> = if force_reindex {
        (0..alignment_files.len()).collect()
    } else {
        index_paths
            .iter()
            .enumerate()
            .filter(|(_, path)| !path.exists())
            .map(|(i, _)| i)
            .collect()
    };

    if !indices_to_build.is_empty() {
        info!("Building {} per-file indices...", indices_to_build.len());

        // Build indices in parallel
        let build_results: Vec<io::Result<()>> = indices_to_build
            .par_iter()
            .map(|&i| {
                let aln_file = &alignment_files[i];
                let index_path = &index_paths[i];

                debug!("Building index for {}", aln_file);

                // Build single-file index
                let single_file = vec![aln_file.clone()];
                let impg = generate_multi_index(
                    &single_file,
                    threads,
                    Some(index_path.to_string_lossy().as_ref()),
                    sequence_files,
                )?;

                // Save the index
                let file = File::create(index_path)?;
                let mut writer = BufWriter::new(file);
                impg.serialize_with_forest_map(&mut writer)?;

                debug!("Built index: {:?}", index_path);
                Ok(())
            })
            .collect();

        // Check for errors
        for result in build_results {
            result?;
        }

        info!("Finished building {} per-file indices", indices_to_build.len());
    }

    // Load MultiImpg from all per-file indices
    info!("Loading per-file indices into MultiImpg...");
    let multi = MultiImpg::load_from_files(&index_paths, alignment_files, sequence_files)?;

    Ok(ImpgWrapper::from_multi(multi))
}

/// Resolve the list of alignment files (PAF or .1aln) from either --alignment-files or --alignment-list
fn resolve_alignment_files(alignment: &AlignmentOpts) -> io::Result<Vec<String>> {
    let alignment_files = if !alignment.alignment_files.is_empty() {
        alignment.alignment_files.clone()
    } else if let Some(alignment_list_file) = &alignment.alignment_list {
        // Read alignment files from the list file
        let file = File::open(alignment_list_file)?;
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
                format!("No valid alignment files found in list file: {alignment_list_file}"),
            ));
        }

        files
    } else {
        // Neither alignment_files nor alignment_list provided
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Either --alignment-files or --alignment-list must be provided",
        ));
    };

    // Check if the number of PAF files exceeds u32::MAX
    if alignment_files.len() > u32::MAX as usize {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Too many PAF files specified: {} (maximum allowed: {})",
                alignment_files.len(),
                u32::MAX
            ),
        ));
    }

    Ok(alignment_files)
}

fn load_or_generate_multi_index(
    alignment_files: &[String],
    threads: NonZeroUsize,
    custom_index: Option<&str>,
    sequence_files: Option<&[String]>,
) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(alignment_files, custom_index);
    if std::path::Path::new(&index_file).exists() {
        load_multi_index(alignment_files, custom_index, sequence_files)
    } else {
        generate_multi_index(alignment_files, threads, custom_index, sequence_files)
    }
}

fn load_multi_index(
    alignment_files: &[String],
    custom_index: Option<&str>,
    sequence_files: Option<&[String]>,
) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(alignment_files, custom_index);
    info!("Reading IMPG index from {index_file}");

    // Check if all alignment files are newer than the index
    let index_file_metadata = std::fs::metadata(&index_file)?;
    let index_file_ts = index_file_metadata.modified().ok();

    if let Some(index_ts) = index_file_ts {
        alignment_files.par_iter().for_each(|alignment_file| {
            if let Ok(alignment_file_metadata) = std::fs::metadata(alignment_file) {
                if let Ok(alignment_file_ts) = alignment_file_metadata.modified() {
                    if alignment_file_ts > index_ts {
                        warn!(
                            "WARNING:\tAlignment file {alignment_file} has been modified since impg index creation."
                        );
                    }
                }
            }
        });
    }

    // Load from embedded format
    let file = File::open(&index_file)?;
    let reader = BufReader::new(file);

    Impg::load_from_file(reader, alignment_files, index_file, sequence_files)
}

fn generate_multi_index(
    alignment_files: &[String],
    threads: NonZeroUsize,
    custom_index: Option<&str>,
    sequence_files: Option<&[String]>,
) -> io::Result<Impg> {
    let index_file = get_combined_index_filename(alignment_files, custom_index);
    info!("No index found at {index_file}. Creating it now.");

    let num_alignment_files = alignment_files.len();
    // Thread-safe counter for tracking progress
    let files_processed = AtomicUsize::new(0);

    // Process alignment files in parallel, using a per-file local SequenceIndex to avoid global locking
    // We will merge names after parsing and remap IDs in a parallel pass.
    let mut records_by_file_with_local_index: Vec<(Vec<AlignmentRecord>, String, SequenceIndex)> =
        (0..alignment_files.len())
            .into_par_iter()
            .map(|file_index| -> io::Result<(Vec<AlignmentRecord>, String, SequenceIndex)> {
                let aln_file = &alignment_files[file_index];

                // Increment the counter and get the new value atomically
                let current_count = files_processed.fetch_add(1, Ordering::SeqCst) + 1;
                debug!("Processing alignment file ({current_count}/{num_alignment_files}): {aln_file}");

                // Local sequence index for this file only
                let mut local_seq_index = SequenceIndex::new();

            // Detect file format and parse accordingly
            let format = AlignmentFormat::from_path(aln_file).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Unsupported alignment format: {aln_file}"),
                )
            })?;

                let records = match format {
                    AlignmentFormat::Paf => {
                        let file = File::open(aln_file).map_err(|e| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("Failed to open PAF file: {}", e),
                            )
                        })?;
                    impg::paf::parse_paf_file(aln_file, file, threads, &mut local_seq_index)
                        .map_err(|e| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("Failed to parse PAF records: {}", e),
                            )
                        })?
                    }
                    AlignmentFormat::OneAln => {
                        let parser =
                            OneAlnParser::new(aln_file.clone(), sequence_files).map_err(|e| {
                                io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("Failed to create 1aln parser: {}", e),
                                )
                            })?;
                        parser.parse_alignments(&mut local_seq_index).map_err(|e| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("Failed to parse 1aln records: {}", e),
                            )
                        })?
                    }
                };
                debug!(
                    "Parsed {} alignment records from file: {aln_file}",
                    records.len()
                );
                Ok((records, aln_file.clone(), local_seq_index))
            })
        .collect::<Result<Vec<_>, _>>()?; // Propagate any errors

    // Build a temporary union of all sequence names and lengths from local indices
    let mut tmp_union_index = SequenceIndex::new();
    for (_, _, local_idx) in &records_by_file_with_local_index {
        for (name, id) in &local_idx.name_to_id {
            let length = local_idx.get_len_from_id(*id).unwrap_or(0);
            tmp_union_index.get_or_insert_id(name, Some(length));
        }
    }

    // Sort sequence names to ensure deterministic order
    let mut sequence_names = tmp_union_index
        .name_to_id
        .keys()
        .cloned()
        .collect::<Vec<String>>();
    sequence_names.par_sort_unstable(); // Order of identical sequence names is irrelevant

    // Create a deterministic global SequenceIndex
    let mut seq_index = SequenceIndex::new();
    for (name, id) in &tmp_union_index.name_to_id {
        let length = tmp_union_index.get_len_from_id(*id).unwrap_or(0);
        seq_index.get_or_insert_id(name, Some(length));
    }

    // Update query and target IDs with the new SequenceIndex
    records_by_file_with_local_index
        .par_iter_mut()
        .for_each(|(records, _path, local_idx)| {
            for record in records.iter_mut() {
                let query_name = local_idx.get_name(record.query_id).unwrap();
                record.query_id = seq_index.get_id(query_name).unwrap();

                let target_name = local_idx.get_name(record.target_id).unwrap();
                record.target_id = seq_index.get_id(target_name).unwrap();
            }
        });

    // Drop local indices and keep only (records, path)
    let records_by_file: Vec<(Vec<AlignmentRecord>, String)> = records_by_file_with_local_index
        .into_iter()
        .map(|(records, path, _)| (records, path))
        .collect();

    let impg = Impg::from_multi_alignment_records(&records_by_file, seq_index, sequence_files)
        .map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Failed to create index: {e}"),
            )
        })?;

    // Serialize the index with embedded forest map
    let index_file_path = index_file.clone();
    let file = File::create(&index_file_path)?;
    let mut writer = BufWriter::new(file);
    impg.serialize_with_forest_map(&mut writer)?;

    Ok(impg)
}

fn get_combined_index_filename(alignment_files: &[String], custom_index: Option<&str>) -> String {
    if let Some(index) = custom_index {
        return index.to_string();
    }

    if alignment_files.len() == 1 {
        format!("{}.impg", alignment_files[0])
    } else {
        // For multiple files, create a hash of the sorted filenames

        let mut file_refs: Vec<&str> = alignment_files.iter().map(|s| s.as_str()).collect();
        file_refs.sort();

        let mut hasher = DefaultHasher::new();
        for file in &file_refs {
            file.hash(&mut hasher);
        }

        format!("combined_{:016x}.impg", hasher.finish())
    }
}

fn perform_query(
    impg: &impl ImpgIndex,
    target_name: &str,
    target_range: (i32, i32),
    store_cigar: bool,
    min_identity: Option<f64>,
    min_output_length: Option<i32>,
    transitive: bool,
    transitive_dfs: bool,
    transitive_opts: &TransitiveOpts,
    sequence_index: Option<&UnifiedSequenceIndex>,
    approximate_mode: bool,
    subset_filter: Option<&SubsetFilter>,
) -> io::Result<Vec<AdjustedInterval>> {
    let (target_start, target_end) = target_range;
    let target_id = impg.seq_index().get_id(target_name).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Target sequence '{target_name}' not found in index"),
        )
    })?;
    let target_length = impg.seq_index().get_len_from_id(target_id).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Target sequence '{target_name}' length not found in index"),
        )
    })?;
    if target_end > target_length as i32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Target range end ({target_end}) exceeds the target sequence length ({target_length})"
            ),
        ));
    }

    let results = if transitive {
        impg.query_transitive_bfs(
            target_id,
            target_start,
            target_end,
            None,
            transitive_opts.max_depth,
            transitive_opts.effective_min_transitive_len(),
            transitive_opts.min_distance_between_ranges,
            min_output_length,
            store_cigar,
            min_identity,
            sequence_index,
            approximate_mode,
            subset_filter,
        )
    } else if transitive_dfs {
        impg.query_transitive_dfs(
            target_id,
            target_start,
            target_end,
            None,
            transitive_opts.max_depth,
            transitive_opts.effective_min_transitive_len(),
            transitive_opts.min_distance_between_ranges,
            min_output_length,
            store_cigar,
            min_identity,
            sequence_index,
            approximate_mode,
            subset_filter,
        )
    } else {
        let mut res = impg.query(
            target_id,
            target_start,
            target_end,
            store_cigar,
            min_identity,
            sequence_index,
            approximate_mode,
        );
        // Filter by minimum output length for regular queries
        if let Some(min_len) = min_output_length {
            res.retain(|(query_interval, _, _)| {
                let length = (query_interval.last - query_interval.first).abs();
                length >= min_len
            });
        }

        // Apply subset filter for non-transitive queries (transitive queries filter during exploration)
        if let Some(filter) = subset_filter {
            apply_subset_filter(impg, target_id, &mut res, Some(filter));
        }

        res
    };

    info!(
        "Collected {} results (excluding input range)",
        results.len() - 1
    ); // Exclude the first element (the input range itself)

    Ok(results)
}

/// Load subset filter if path is provided
fn load_subset_filter_if_provided(path: &Option<String>) -> io::Result<Option<SubsetFilter>> {
    if let Some(ref list_path) = path {
        let filter = load_subset_filter(list_path)?;
        info!(
            "Loaded {} sequence names from subset list {}",
            filter.entry_count(),
            list_path
        );
        Ok(Some(filter))
    } else {
        Ok(None)
    }
}

fn output_results_bed(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    name: &str,
    merge_distance: i32,
    merge_strands: bool,
    original_coordinates: bool,
) -> io::Result<()> {
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    for (query_interval, _, _) in results {
        let query_name = impg.seq_index().get_name(query_interval.metadata).unwrap();
        let (first, last, strand) = if query_interval.first <= query_interval.last {
            (query_interval.first, query_interval.last, '+')
        } else {
            (query_interval.last, query_interval.first, '-')
        };

        // Transform coordinates to original sequence space if requested
        let (transformed_name, transformed_first, transformed_last) =
            transform_coordinates_to_original(
                query_name,
                first as u32,
                last as u32,
                original_coordinates,
            );

        writeln!(
            out,
            "{}\t{}\t{}\t{}\t.\t{}",
            transformed_name, transformed_first, transformed_last, name, strand
        )?;
    }

    Ok(())
}

fn output_results_bedpe(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    name: &str,
    merge_distance: i32,
    original_coordinates: bool,
) -> io::Result<()> {
    merge_adjusted_intervals(results, merge_distance);

    for (overlap_query, cigar, overlap_target) in results {
        let query_name = impg.seq_index().get_name(overlap_query.metadata).unwrap();
        let target_name = impg.seq_index().get_name(overlap_target.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };

        // Transform coordinates to original sequence space if requested
        let (transformed_query_name, transformed_first, transformed_last) =
            transform_coordinates_to_original(
                query_name,
                first as u32,
                last as u32,
                original_coordinates,
            );
        let (transformed_target_name, transformed_target_first, transformed_target_last) =
            transform_coordinates_to_original(
                target_name,
                overlap_target.first as u32,
                overlap_target.last as u32,
                original_coordinates,
            );

        // Calculate gap-compressed-identity and block-identity from CIGAR
        let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, _block_len) =
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
            (matches as f32) / (matches + mismatches + insertions + deletions) as f32;

        let edit_distance = mismatches + inserted_bp + deleted_bp;
        let block_identity = (matches as f32) / (matches + edit_distance) as f32;

        // Format gi and bi fields without trailing zeros
        let gi_str = format!("{gap_compressed_identity:.6}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();
        let bi_str = format!("{block_identity:.6}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();

        // BEDPE supports any number of additional, fields after the standard 10 fields.
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t+\tgi:f:{}\tbi:f:{}",
            transformed_query_name,
            transformed_first,
            transformed_last,
            transformed_target_name,
            transformed_target_first,
            transformed_target_last,
            name,
            strand,
            gi_str,
            bi_str
        )?;
    }

    Ok(())
}

fn output_results_paf(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    name: &str,
    merge_distance: i32,
    original_coordinates: bool,
    sequence_index: Option<&UnifiedSequenceIndex>,
) -> io::Result<()> {
    merge_adjusted_intervals(results, merge_distance);

    for (overlap_query, cigar, overlap_target) in results {
        let query_name = impg.seq_index().get_name(overlap_query.metadata).unwrap();
        let target_name = impg.seq_index().get_name(overlap_target.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };

        // Transform coordinates to original sequence space if requested
        let (transformed_query_name, transformed_first, transformed_last) =
            transform_coordinates_to_original(
                query_name,
                first as u32,
                last as u32,
                original_coordinates,
            );
        let (transformed_target_name, transformed_target_first, transformed_target_last) =
            transform_coordinates_to_original(
                target_name,
                overlap_target.first as u32,
                overlap_target.last as u32,
                original_coordinates,
            );

        // Get original sequence lengths when original_sequence_coordinates is enabled
        let query_length = if original_coordinates {
            get_original_sequence_length(&transformed_query_name, sequence_index)
        } else {
            impg.seq_index()
                .get_len_from_id(overlap_query.metadata)
                .unwrap()
        };

        let target_length = if original_coordinates {
            get_original_sequence_length(&transformed_target_name, sequence_index)
        } else {
            impg.seq_index()
                .get_len_from_id(overlap_target.metadata)
                .unwrap()
        };

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
            (matches as f32) / (matches + mismatches + insertions + deletions) as f32;

        let edit_distance = mismatches + inserted_bp + deleted_bp;
        let block_identity = (matches as f32) / (matches + edit_distance) as f32;

        // Format bi and gi fields without trailing zeros
        let gi_str = format!("{gap_compressed_identity:.6}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();
        let bi_str = format!("{block_identity:.6}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string();

        let cigar_str: String = cigar
            .iter()
            .map(|op| format!("{}{}", op.len(), op.op()))
            .collect();

        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tgi:f:{}\tbi:f:{}\tcg:Z:{}\tan:Z:{}",
            transformed_query_name,
            query_length,
            transformed_first,
            transformed_last,
            strand,
            transformed_target_name,
            target_length,
            transformed_target_first,
            transformed_target_last,
            matches,
            block_len,
            255,
            gi_str,
            bi_str,
            cigar_str,
            name
        )?;
    }

    Ok(())
}

fn output_results_gfa(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    sequence_index: &UnifiedSequenceIndex,
    _name: &str,
    merge_distance: i32,
    merge_strands: bool,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<()> {
    // Merge intervals if needed
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    // Extract query intervals by consuming results - no cloning
    let query_intervals: Vec<Interval<u32>> = results
        .drain(..)
        .map(|(query_interval, _, _)| query_interval)
        .collect();
    let gfa_output = impg::graph::generate_gfa_from_intervals(
        impg,
        &query_intervals,
        sequence_index,
        scoring_params,
    );
    writeln!(out, "{gfa_output}")?;

    Ok(())
}

fn output_results_fasta(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    sequence_index: &UnifiedSequenceIndex,
    _name: &str,
    merge_distance: i32,
    merge_strands: bool,
    reverse_complement: bool,
) -> io::Result<()> {
    // Merge intervals if needed
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    // Parallelize sequence fetching and processing
    let sequence_data: Vec<(String, String)> = results
        .par_iter()
        .map(|(query_interval, _, _)| -> io::Result<(String, String)> {
            let query_name = impg.seq_index().get_name(query_interval.metadata).unwrap();

            // Determine actual start and end based on orientation
            let (start, end, strand) = if query_interval.first <= query_interval.last {
                (query_interval.first, query_interval.last, '+')
            } else {
                (query_interval.last, query_interval.first, '-')
            };

            // Fetch the sequence
            let sequence = sequence_index.fetch_sequence(query_name, start, end)?;

            // If reverse strand and reverse complementing, reverse complement the sequence
            let sequence = if strand == '-' && reverse_complement {
                impg::graph::reverse_complement(&sequence)
            } else {
                sequence
            };

            // Create header
            let header_suffix = if strand == '-' && reverse_complement {
                "/rc"
            } else {
                ""
            };
            let header = format!(">{query_name}:{start}-{end}{header_suffix}");

            // Convert sequence to string with line breaks every 80 characters
            let sequence_str = String::from_utf8_lossy(&sequence);
            let formatted_sequence = sequence_str
                .as_bytes()
                .chunks(80)
                .map(|chunk| String::from_utf8_lossy(chunk).to_string())
                .collect::<Vec<_>>()
                .join("\n");

            Ok((header, formatted_sequence))
        })
        .collect::<Result<Vec<_>, _>>()?;

    // Output sequences sequentially to maintain order
    for (header, sequence) in sequence_data {
        writeln!(out, "{header}")?;
        writeln!(out, "{sequence}")?;
    }

    Ok(())
}

fn output_results_maf(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    out: &mut dyn Write,
    sequence_index: &UnifiedSequenceIndex,
    _name: &str,
    merge_distance: i32,
    merge_strands: bool,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<()> {
    // Merge intervals if needed
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    let query_intervals: Vec<Interval<u32>> = results
        .drain(..)
        .map(|(query_interval, _, _)| query_interval)
        .collect();

    let maf_output = impg::graph::generate_maf_from_intervals(
        impg,
        &query_intervals,
        sequence_index,
        scoring_params,
    );
    writeln!(out, "{maf_output}")?;

    Ok(())
}

fn output_results_fasta_aln(
    impg: &impl ImpgIndex,
    results: &mut Vec<AdjustedInterval>,
    sequence_index: &UnifiedSequenceIndex,
    _name: String,
    merge_distance: i32,
    merge_strands: bool,
    scoring_params: (u8, u8, u8, u8, u8, u8),
) -> io::Result<()> {
    // Merge intervals as for MAF/GFA (collapse per-query coords, merge strands)
    merge_query_adjusted_intervals(results, merge_distance, merge_strands);

    // Drain query intervals
    let query_intervals: Vec<coitrees::Interval<u32>> =
        results.drain(..).map(|(q, _, _)| q).collect();

    // Ask graph layer to generate aligned FASTA from SPOA's MSA
    let fasta_aln = impg::graph::generate_fasta_alignment_from_intervals(
        impg,
        &query_intervals,
        sequence_index,
        scoring_params,
    );
    print!("{fasta_aln}");
    Ok(())
}

// Merge adjusted intervals by ignoring the target intervals (optimized for simple genomic interval merging in BED and GFA formats)
fn merge_query_adjusted_intervals(
    results: &mut Vec<AdjustedInterval>,
    merge_distance: i32,
    merge_strands: bool,
) {
    if results.len() > 1 && (merge_distance >= 0 || merge_strands) {
        // Sort by sequence ID, start position, and strand (forward first)
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

            // Orientation flags
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

            // Only merge if same sequence and within merge distance (if merge_distance >= 0),
            // and orientation is compatible with merge_strands choice.
            if merge_distance < 0
                || curr_interval.metadata != next_interval.metadata
                || (!merge_strands && curr_is_forward != next_is_forward)
                || next_start > curr_end + merge_distance
            {
                write_idx += 1;
                if write_idx != read_idx {
                    results.swap(write_idx, read_idx);
                }
            } else {
                // Merge intervals, possibly across strands.
                let merged_start = curr_start.min(next_start);
                let merged_end = curr_end.max(next_end);

                // When merging across strands, choose the orientation with the larger span; on ties, keep existing.
                let merged_is_forward = if merge_strands && curr_is_forward != next_is_forward {
                    let curr_len = curr_end.saturating_sub(curr_start);
                    let next_len = next_end.saturating_sub(next_start);
                    if curr_len == next_len {
                        curr_is_forward
                    } else if curr_len > next_len {
                        curr_is_forward
                    } else {
                        next_is_forward
                    }
                } else {
                    curr_is_forward
                };

                if merged_is_forward {
                    results[write_idx].0.first = merged_start;
                    results[write_idx].0.last = merged_end;
                } else {
                    results[write_idx].0.first = merged_end;
                    results[write_idx].0.last = merged_start;
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
                    // debug!(
                    //     "Merge contiguous! Query: current {}:{}-{}({}), next {}:{}-{}({}); Target: current {}:{}-{}({}), next {}:{}-{}({})",
                    //     current_query.metadata,
                    //     current_query.first,
                    //     current_query.last,
                    //     if query_forward { "+" } else { "-" },
                    //     next_query.metadata,
                    //     next_query.first,
                    //     next_query.last,
                    //     if next_query_forward { "+" } else { "-" },
                    //     current_target.metadata,
                    //     current_target.first,
                    //     current_target.last,
                    //     if target_forward { "+" } else { "-" },
                    //     next_target.metadata,
                    //     next_target.first,
                    //     next_target.last,
                    //     if next_target_forward { "+" } else { "-" },
                    // );

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
                            // debug!(
                            //     "Merge overlapping! Overlap: query={}, target={}, Query: current {}:{}-{}({}), next {}:{}-{}({}); Target: current {}:{}-{}({}), next {}:{}-{}({})",
                            //     query_overlap_len,
                            //     target_overlap_len,
                            //     current_query.metadata,
                            //     current_query.first,
                            //     current_query.last,
                            //     if query_forward { "+" } else { "-" },
                            //     next_query.metadata,
                            //     next_query.first,
                            //     next_query.last,
                            //     if next_query_forward { "+" } else { "-" },
                            //     current_target.metadata,
                            //     current_target.first,
                            //     current_target.last,
                            //     if target_forward { "+" } else { "-" },
                            //     next_target.metadata,
                            //     next_target.first,
                            //     next_target.last,
                            //     if next_target_forward { "+" } else { "-" },
                            // );

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
                        // debug!(
                        //     "Merge gaps! Query gap: {}, Target gap: {}, Query: current {}:{}-{}({}), next {}:{}-{}({}); Target: current {}:{}-{}({}), next {}:{}-{}({})",
                        //     query_gap,
                        //     target_gap,
                        //     current_query.metadata,
                        //     current_query.first,
                        //     current_query.last,
                        //     if query_forward { "+" } else { "-" },
                        //     next_query.metadata,
                        //     next_query.first,
                        //     next_query.last,
                        //     if next_query_forward { "+" } else { "-" },
                        //     current_target.metadata,
                        //     current_target.first,
                        //     current_target.last,
                        //     if target_forward { "+" } else { "-" },
                        //     next_target.metadata,
                        //     next_target.first,
                        //     next_target.last,
                        //     if next_target_forward { "+" } else { "-" },
                        // );

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

fn print_stats(impg: &impl ImpgIndex) {
    // Basic stats
    let num_sequences = impg.seq_index().len();
    let total_sequence_length: usize = (0..num_sequences as u32)
        .into_par_iter()
        .filter_map(|id| impg.seq_index().get_len_from_id(id))
        .sum();

    // Compute overlap stats
    let num_targets = impg.num_targets();

    info!("Computing statistics for {} trees...", num_targets);

    // Collect target IDs to process
    let target_ids: Vec<u32> = impg.target_ids();

    // Process trees in parallel, computing stats without keeping them in memory
    let results: Vec<(u32, usize)> = target_ids
        .par_iter()
        .map(|&target_id| {
            // Get tree and count, removing it from memory for efficiency
            let count = if let Some(tree) = impg.get_or_load_tree(target_id) {
                let len = tree.len();
                impg.remove_cached_tree(target_id);
                len
            } else {
                0
            };

            (target_id, count)
        })
        .collect();

    // Collect results
    let mut num_overlaps = 0;
    let mut overlaps_per_seq = FxHashMap::default();

    for (target_id, count) in results {
        num_overlaps += count;
        overlaps_per_seq.insert(target_id, count);
    }

    info!(
        "Processed {} trees, total overlaps: {}",
        target_ids.len(),
        num_overlaps
    );

    println!("Number of query+target sequences: {num_sequences}");
    println!("Total query+target sequence length: {total_sequence_length} bp");
    println!("Number of overlaps: {num_overlaps}");

    let mut entries: Vec<(u32, usize)> = overlaps_per_seq.into_iter().collect();
    if !entries.is_empty() {
        entries.par_sort_by(|a, b| b.1.cmp(&a.1));

        // Calculate mean and median overlaps
        let sum: usize = entries.par_iter().map(|(_, count)| count).sum();
        let mean = sum as f64 / entries.len() as f64;

        let median = if entries.is_empty() {
            0.0
        } else if entries.len().is_multiple_of(2) {
            let mid = entries.len() / 2;
            (entries[mid - 1].1 + entries[mid].1) as f64 / 2.0
        } else {
            entries[entries.len() / 2].1 as f64
        };
        println!("\nMean overlaps per sequence: {mean:.2}");
        println!("Median overlaps per sequence: {median:.2}");

        println!("\nTop target sequences by number of overlaps:");
        for (idx, (seq_id, count)) in entries.iter().take(5).enumerate() {
            if let Some(name) = impg.seq_index().get_name(*seq_id) {
                println!("{}. {}: {} overlaps", idx + 1, name, count);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_subsequence_coordinates() {
        // Test parsing subsequence coordinates
        let result = parse_subsequence_coordinates("HG002#1#chr1:5116130-6116563");
        assert_eq!(result, Some(("HG002#1#chr1".to_string(), 5116130)));

        let result = parse_subsequence_coordinates("GRCh38#0#chr1:5477602-6474357");
        assert_eq!(result, Some(("GRCh38#0#chr1".to_string(), 5477602)));

        // Test with no subsequence coordinates
        let result = parse_subsequence_coordinates("chr1");
        assert_eq!(result, None);

        // Test with invalid format
        let result = parse_subsequence_coordinates("chr1:invalid");
        assert_eq!(result, None);
    }

    #[test]
    fn test_transform_coordinates_to_original() {
        // Test with original coordinates enabled
        let (name, start, end) =
            transform_coordinates_to_original("HG002#1#chr1:5116130-6116563", 45803, 45861, true);
        assert_eq!(name, "HG002#1#chr1");
        assert_eq!(start, 5116130 + 45803);
        assert_eq!(end, 5116130 + 45861);

        // Test with original coordinates disabled
        let (name, start, end) =
            transform_coordinates_to_original("HG002#1#chr1:5116130-6116563", 45803, 45861, false);
        assert_eq!(name, "HG002#1#chr1:5116130-6116563");
        assert_eq!(start, 45803);
        assert_eq!(end, 45861);

        // Test with sequence name that doesn't contain subsequence coordinates
        let (name, start, end) = transform_coordinates_to_original("chr1", 45803, 45861, true);
        assert_eq!(name, "chr1");
        assert_eq!(start, 45803);
        assert_eq!(end, 45861);
    }
}
