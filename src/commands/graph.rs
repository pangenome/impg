// Graph building command using sweepga + seqwish integration
//
// This module provides functionality to build pangenome graphs from FASTA sequences
// by combining sweepga (for alignment) and seqwish (for graph induction).

use log::info;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::{Arc, Mutex, RwLock};
use std::time::Instant;

use bitvec::prelude::*;

// Import gfasort for graph sorting
use crate::graph::sort_gfa;

// Import from sweepga
use sweepga::fastga_integration::FastGAIntegration;
use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};

// Import from seqwish
use seqwish::alignments::unpack_paf_alignments;
use seqwish::compact::compact_nodes;
use seqwish::gfa::emit_gfa;
use seqwish::intervaltree::{AdaptiveTree, IntervalTree};
use seqwish::links::{derive_links, RankSelectBitVector};
use seqwish::seqindex::SeqIndex;
use seqwish::transclosure::compute_transitive_closures;

/// Configuration for graph building
pub struct GraphBuildConfig {
    /// Number of threads for parallel processing
    pub num_threads: usize,
    /// K-mer frequency multiplier (genomes * this value)
    pub frequency_multiplier: usize,
    /// Explicit frequency override (if set, overrides multiplier)
    pub frequency: Option<usize>,
    /// Minimum alignment length for FastGA
    pub min_alignment_length: u64,
    /// Maximum repeat count for transitive closure
    pub repeat_max: u64,
    /// Minimum repeat distance for transitive closure
    pub min_repeat_dist: u64,
    /// Minimum match length filter for alignments
    pub min_match_len: u64,
    /// Sparse factor for input matches
    pub sparse_factor: f32,
    /// Batch size for transitive closure computation
    pub transclose_batch: u64,
    /// Whether to use in-memory interval trees
    pub use_in_memory: bool,
    /// Show progress during graph building
    pub show_progress: bool,
    /// Directory for temporary files
    pub temp_dir: Option<String>,

    // Sweepga filtering options
    /// Disable all filtering (default: filtering enabled)
    pub no_filter: bool,
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
}

impl Default for GraphBuildConfig {
    fn default() -> Self {
        GraphBuildConfig {
            num_threads: 4,
            frequency_multiplier: 10,
            frequency: None,
            min_alignment_length: 100,
            repeat_max: 0,
            min_repeat_dist: 0,
            min_match_len: 0,
            sparse_factor: 0.0,
            transclose_batch: 1_000_000,
            use_in_memory: true,
            show_progress: true,
            temp_dir: None,
            // Filtering options with sensible defaults
            no_filter: false,
            num_mappings: "1:1".to_string(),
            scaffold_jump: 50_000,       // 50kb default scaffold gap
            scaffold_mass: 10_000,        // 10kb minimum scaffold length
            scaffold_filter: "1:1".to_string(),
            overlap: 0.95,
            min_identity: 0.0,
        }
    }
}

/// Count sequences and genomes in FASTA files
/// Returns (num_sequences, num_genomes)
/// For PanSN-spec names like SAMPLE#HAPLOTYPE#CONTIG, genomes are unique SAMPLE#HAPLOTYPE prefixes
fn count_sequences_and_genomes_in_fasta(fasta_paths: &[String]) -> io::Result<(usize, usize)> {
    use std::collections::HashSet;
    let mut seq_count = 0;
    let mut genome_prefixes: HashSet<String> = HashSet::new();

    for path in fasta_paths {
        let file = File::open(path)?;
        // Use niffler to auto-detect compression
        let (reader, _format) = niffler::get_reader(Box::new(file))
            .map_err(|e| io::Error::other(format!("Failed to open reader for '{}': {}", path, e)))?;
        let reader = BufReader::new(reader);

        for line in reader.lines() {
            let line: String = line?;
            if line.starts_with('>') {
                seq_count += 1;
                // Extract genome prefix: everything before the last # for PanSN names
                // For ">SAMPLE#HAPLOTYPE#CONTIG", prefix is "SAMPLE#HAPLOTYPE"
                let name = line[1..].split_whitespace().next().unwrap_or("");
                let parts: Vec<&str> = name.split('#').collect();
                let prefix = if parts.len() >= 2 {
                    // PanSN format: use SAMPLE#HAPLOTYPE (first two parts)
                    format!("{}#{}", parts[0], parts[1])
                } else {
                    // Not PanSN: use whole name as prefix
                    name.to_string()
                };
                genome_prefixes.insert(prefix);
            }
        }
    }

    let genome_count = genome_prefixes.len().max(1); // At least 1 genome
    Ok((seq_count, genome_count))
}

/// Parse filter mode string (e.g., "1:1", "many:many", "5:3") into FilterMode and limits
fn parse_filter_mode(s: &str) -> (FilterMode, Option<usize>, Option<usize>) {
    let s_lower = s.to_lowercase();
    if s_lower == "many:many" || s_lower == "n:n" {
        return (FilterMode::ManyToMany, None, None);
    }

    // Parse M:N format
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 2 {
        // Default to 1:1 if invalid format
        return (FilterMode::OneToOne, Some(1), Some(1));
    }

    let query_max = if parts[0] == "many" || parts[0] == "n" {
        None
    } else {
        parts[0].parse().ok()
    };

    let target_max = if parts[1] == "many" || parts[1] == "n" {
        None
    } else {
        parts[1].parse().ok()
    };

    // FilterMode determines the filtering strategy:
    // - OneToOne: enforces strict 1:1 mapping
    // - OneToMany: 1:N (one per query, N per target)
    // - ManyToMany: N:N (no strict filtering, uses overlap-based plane sweep)
    // The actual limits are controlled by max_per_query/max_per_target
    match (query_max, target_max) {
        (Some(1), Some(1)) => (FilterMode::OneToOne, Some(1), Some(1)),
        (Some(1), _) => (FilterMode::OneToMany, Some(1), target_max),
        (_, Some(1)) => (FilterMode::OneToMany, query_max, Some(1)),
        (Some(q), Some(t)) => (FilterMode::ManyToMany, Some(q), Some(t)),
        (Some(q), None) => (FilterMode::ManyToMany, Some(q), None),
        (None, Some(t)) => (FilterMode::ManyToMany, None, Some(t)),
        (None, None) => (FilterMode::ManyToMany, None, None),
    }
}

/// Build a pangenome graph from FASTA sequences
///
/// # Arguments
/// * `fasta_files` - List of FASTA file paths (can be gzip compressed)
/// * `output_gfa` - Output GFA file path (or "-" for stdout)
/// * `config` - Configuration options for graph building
///
/// # Returns
/// The number of nodes in the resulting graph
pub fn build_graph<W: Write>(
    fasta_files: &[String],
    output: &mut W,
    config: &GraphBuildConfig,
) -> io::Result<usize> {
    let start_time = Instant::now();

    // Validate inputs
    if fasta_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No FASTA files provided",
        ));
    }

    for path in fasta_files {
        if !Path::new(path).exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("FASTA file not found: {}", path),
            ));
        }
    }

    // Set up temp directory for all temp file operations (seqwish, sweepga, tempfile crate)
    // Uses system default (/tmp) unless --temp-dir is specified
    // For better I/O performance on systems with RAM-backed tmpfs, use --temp-dir /dev/shm
    if let Some(ref temp_dir) = config.temp_dir {
        // Set TMPDIR environment variable so Rust's tempfile crate uses it
        std::env::set_var("TMPDIR", temp_dir);
        // Also configure seqwish's internal temp file handling
        seqwish::tempfile::set_dir(temp_dir);
        if config.show_progress {
            info!(
                "[graph::temp] Using temp directory: {}",
                temp_dir
            );
        }
    }

    // 1) Count sequences and genomes, determine k-mer frequency
    if config.show_progress {
        info!(
            "[graph::count] {:.3}s Counting sequences in {} FASTA file(s)",
            start_time.elapsed().as_secs_f64(),
            fasta_files.len()
        );
    }

    let (num_sequences, num_genomes) = count_sequences_and_genomes_in_fasta(fasta_files)?;
    // Use num_genomes (not num_sequences) for frequency - matches sweepga behavior
    let kmer_frequency = config
        .frequency
        .unwrap_or(num_genomes * config.frequency_multiplier);

    if config.show_progress {
        info!(
            "[graph::count] {:.3}s Found {} sequences in {} genomes, using k-mer frequency {}",
            start_time.elapsed().as_secs_f64(),
            num_sequences,
            num_genomes,
            kmer_frequency
        );
    }

    // 2) Create combined FASTA for alignment
    // sweepga's FastGA needs a single combined FASTA file for all-vs-all alignment
    // FastGA requires .fa extension to recognize the file format
    let combined_fasta = tempfile::Builder::new()
        .suffix(".fa")
        .tempfile()?;
    {
        let mut writer = BufWriter::new(&combined_fasta);
        for path in fasta_files {
            let file = File::open(path)?;
            // Use niffler to auto-detect compression
            let (reader, _format) = niffler::get_reader(Box::new(file))
                .map_err(|e| io::Error::other(format!("Failed to open reader for '{}': {}", path, e)))?;
            let reader = BufReader::new(reader);
            for line in reader.lines() {
                let line: String = line?;
                writeln!(writer, "{}", line)?;
            }
        }
        writer.flush()?;
    }

    if config.show_progress {
        info!(
            "[graph::combine] {:.3}s Combined FASTA files for alignment",
            start_time.elapsed().as_secs_f64()
        );
    }

    // 3) Run sweepga/FastGA alignment
    if config.show_progress {
        info!(
            "[graph::align] {:.3}s Running FastGA alignment with -f {}",
            start_time.elapsed().as_secs_f64(),
            kmer_frequency
        );
    }

    let fastga = FastGAIntegration::new(
        Some(kmer_frequency),
        config.num_threads,
        config.min_alignment_length,
        config.temp_dir.clone(),
    );

    // Run all-vs-all alignment (query = target = combined FASTA)
    let paf_temp = fastga
        .align_to_temp_paf(combined_fasta.path(), combined_fasta.path())
        .map_err(|e| io::Error::other(format!("FastGA alignment failed: {}", e)))?;

    if config.show_progress {
        info!(
            "[graph::align] {:.3}s Alignment complete",
            start_time.elapsed().as_secs_f64()
        );
    }

    // 3.5) Apply sweepga filtering to alignments (unless --no-filter)
    let filtered_paf = if config.no_filter {
        if config.show_progress {
            info!(
                "[graph::filter] {:.3}s Filtering disabled (--no-filter)",
                start_time.elapsed().as_secs_f64()
            );
        }
        paf_temp
    } else {
        if config.show_progress {
            info!(
                "[graph::filter] {:.3}s Filtering alignments with sweepga (n={}, scaffold_jump={}, scaffold_mass={}, scaffold_filter={})",
                start_time.elapsed().as_secs_f64(),
                config.num_mappings,
                config.scaffold_jump,
                config.scaffold_mass,
                config.scaffold_filter
            );
        }

        // Parse filter modes
        let (mapping_mode, mapping_per_query, mapping_per_target) =
            parse_filter_mode(&config.num_mappings);
        let (scaffold_mode, scaffold_per_query, scaffold_per_target) =
            parse_filter_mode(&config.scaffold_filter);

        // Create filter configuration
        let filter_config = FilterConfig {
            chain_gap: 0,
            min_block_length: 0,
            mapping_filter_mode: mapping_mode,
            mapping_max_per_query: mapping_per_query,
            mapping_max_per_target: mapping_per_target,
            plane_sweep_secondaries: 0,
            scaffold_filter_mode: scaffold_mode,
            scaffold_max_per_query: scaffold_per_query,
            scaffold_max_per_target: scaffold_per_target,
            overlap_threshold: config.overlap,
            sparsity: 1.0,
            no_merge: true,
            scaffold_gap: config.scaffold_jump,
            min_scaffold_length: config.scaffold_mass,
            scaffold_overlap_threshold: 0.5,
            scaffold_max_deviation: 0,
            prefix_delimiter: '#',
            skip_prefix: false,
            scoring_function: ScoringFunction::LogLengthIdentity,
            min_identity: config.min_identity,
            min_scaffold_identity: config.min_identity,
        };

        // Create filtered PAF temp file
        let filtered_paf_file = tempfile::Builder::new()
            .suffix(".filtered.paf")
            .tempfile()?;

        // Apply filtering
        let filter = PafFilter::new(filter_config).with_keep_self(false);
        filter
            .filter_paf(paf_temp.path(), filtered_paf_file.path())
            .map_err(|e| io::Error::other(format!("Filtering failed: {}", e)))?;

        if config.show_progress {
            info!(
                "[graph::filter] {:.3}s Filtering complete",
                start_time.elapsed().as_secs_f64()
            );
        }

        filtered_paf_file
    };

    // 4) Build sequence index for seqwish
    if config.show_progress {
        info!(
            "[graph::seqindex] {:.3}s Building sequence index",
            start_time.elapsed().as_secs_f64()
        );
    }

    let mut seqidx = SeqIndex::new();
    seqidx
        .build_index(combined_fasta.path().to_str().unwrap())
        .map_err(io::Error::other)?;
    let seqidx = Arc::new(seqidx);

    if config.show_progress {
        info!(
            "[graph::seqindex] {:.3}s Indexed {} sequences",
            start_time.elapsed().as_secs_f64(),
            seqidx.n_seqs()
        );
    }

    // 5) Index alignments into interval tree
    if config.show_progress {
        info!(
            "[graph::alignments] {:.3}s Loading and indexing alignments",
            start_time.elapsed().as_secs_f64()
        );
    }

    let aln_iitree_idx = seqwish::tempfile::create("seqwish-", ".sqa")?;
    let mut aln_iitree_obj = if config.use_in_memory {
        AdaptiveTree::new_memory()?
    } else {
        AdaptiveTree::new_disk(&aln_iitree_idx)?
    };
    aln_iitree_obj.open_writer()?;
    let aln_iitree = Arc::new(Mutex::new(aln_iitree_obj));

    unpack_paf_alignments(
        filtered_paf.path().to_str().unwrap(),
        Arc::clone(&aln_iitree),
        Arc::clone(&seqidx),
        config.min_match_len,
        config.sparse_factor,
        config.num_threads,
    )?;

    aln_iitree.lock().unwrap().index()?;

    if config.show_progress {
        info!(
            "[graph::alignments] {:.3}s Alignments indexed",
            start_time.elapsed().as_secs_f64()
        );
    }

    // Unwrap Mutex - alignment tree is read-only during graph construction
    let aln_iitree_readonly = Arc::new(
        Arc::try_unwrap(aln_iitree)
            .map_err(|_| io::Error::other("Failed to unwrap alignment tree Arc"))?
            .into_inner()
            .unwrap(),
    );

    // 6) Compute transitive closures
    if config.show_progress {
        info!(
            "[graph::transclosure] {:.3}s Computing transitive closures",
            start_time.elapsed().as_secs_f64()
        );
    }

    let seq_v_file = seqwish::tempfile::create("seqwish-", ".sqs")?;
    let node_iitree_idx = seqwish::tempfile::create("seqwish-", ".sqn")?;
    let path_iitree_idx = seqwish::tempfile::create("seqwish-", ".sqp")?;

    let mut node_iitree_obj = if config.use_in_memory {
        AdaptiveTree::new_memory()?
    } else {
        AdaptiveTree::new_disk(&node_iitree_idx)?
    };
    node_iitree_obj.open_writer()?;
    let node_iitree = Arc::new(RwLock::new(node_iitree_obj));

    let mut path_iitree_obj = if config.use_in_memory {
        AdaptiveTree::new_memory()?
    } else {
        AdaptiveTree::new_disk(&path_iitree_idx)?
    };
    path_iitree_obj.open_writer()?;
    let path_iitree = Arc::new(RwLock::new(path_iitree_obj));

    let graph_length = compute_transitive_closures(
        Arc::clone(&seqidx),
        Arc::clone(&aln_iitree_readonly),
        seq_v_file.to_str().unwrap(),
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        config.repeat_max,
        config.min_repeat_dist,
        config.transclose_batch,
        config.show_progress,
        config.num_threads,
    )?;

    if config.show_progress {
        info!(
            "[graph::transclosure] {:.3}s Transitive closure complete (graph length: {})",
            start_time.elapsed().as_secs_f64(),
            graph_length
        );
    }

    // 7) Compact nodes
    if config.show_progress {
        info!(
            "[graph::compact] {:.3}s Compacting nodes",
            start_time.elapsed().as_secs_f64()
        );
    }

    let mut seq_id_bv = BitVec::<u64, Lsb0>::repeat(false, graph_length + 1);

    compact_nodes(
        Arc::clone(&seqidx),
        graph_length,
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        &mut seq_id_bv,
        config.num_threads,
    )?;

    // Build rank/select structure
    let seq_id_cbv =
        RankSelectBitVector::from_bitvec(&seq_id_bv.iter().by_vals().collect::<Vec<bool>>());
    let node_count = seq_id_cbv.rank(graph_length);
    drop(seq_id_bv); // Free memory

    if config.show_progress {
        info!(
            "[graph::compact] {:.3}s Compaction complete ({} nodes)",
            start_time.elapsed().as_secs_f64(),
            node_count
        );
    }

    // 8) Derive links between nodes
    if config.show_progress {
        info!(
            "[graph::links] {:.3}s Finding graph links",
            start_time.elapsed().as_secs_f64()
        );
    }

    let link_set = derive_links(
        Arc::clone(&seqidx),
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        &seq_id_cbv,
        config.num_threads,
    )?;

    if config.show_progress {
        info!(
            "[graph::links] {:.3}s Links derived ({} links)",
            start_time.elapsed().as_secs_f64(),
            link_set.len()
        );
    }

    // 9) Emit GFA output to buffer
    if config.show_progress {
        info!(
            "[graph::gfa] {:.3}s Generating GFA output",
            start_time.elapsed().as_secs_f64()
        );
    }

    let mut gfa_buffer = Vec::new();
    {
        let mut buffer_writer = BufWriter::new(&mut gfa_buffer);
        emit_gfa(
            &mut buffer_writer,
            graph_length,
            seq_v_file.to_str().unwrap(),
            Arc::clone(&node_iitree),
            Arc::clone(&path_iitree),
            &seq_id_cbv,
            Arc::clone(&seqidx),
            link_set.links(),
            config.num_threads,
        )?;
    }

    let gfa_string = String::from_utf8(gfa_buffer)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Invalid UTF-8 in GFA: {}", e)))?;

    // 10) Sort GFA using gfasort's Ygs pipeline (path-guided SGD + grooming + topological sort)
    if config.show_progress {
        info!(
            "[graph::sort] {:.3}s Sorting graph with gfasort (Ygs pipeline)",
            start_time.elapsed().as_secs_f64()
        );
    }

    let sorted_gfa = sort_gfa(&gfa_string, config.num_threads)?;

    // Write sorted GFA to output
    output.write_all(sorted_gfa.as_bytes())?;

    if config.show_progress {
        info!(
            "[graph::gfa] {:.3}s Graph building complete",
            start_time.elapsed().as_secs_f64()
        );
    }

    // Clean up seqwish temp files before returning
    // This is critical when building multiple graphs in the same process
    seqwish::tempfile::cleanup();

    Ok(node_count)
}

/// Run the graph building command
pub fn run_graph_build(
    fasta_files: Vec<String>,
    fasta_list: Option<String>,
    output: &str,
    config: GraphBuildConfig,
) -> io::Result<()> {
    // Resolve FASTA files
    let fasta_files = resolve_fasta_files(fasta_files, fasta_list)?;

    if fasta_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No FASTA files specified. Use --fasta-files or --fasta-list",
        ));
    }

    info!("Building graph from {} FASTA file(s)", fasta_files.len());

    // Build graph with output
    if output == "-" {
        let stdout = io::stdout();
        let mut out = BufWriter::with_capacity(1024 * 1024, stdout.lock());
        build_graph(&fasta_files, &mut out, &config)?;
    } else {
        let mut out = BufWriter::with_capacity(1024 * 1024, File::create(output)?);
        build_graph(&fasta_files, &mut out, &config)?;
    }

    Ok(())
}

/// Resolve FASTA files from either direct list or file list
fn resolve_fasta_files(
    fasta_files: Vec<String>,
    fasta_list: Option<String>,
) -> io::Result<Vec<String>> {
    match (fasta_files.is_empty(), fasta_list) {
        (false, None) => {
            // Validate all files exist
            for file in &fasta_files {
                if !Path::new(file).exists() {
                    return Err(io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("FASTA file '{}' not found", file),
                    ));
                }
            }
            Ok(fasta_files)
        }
        (true, Some(list_file)) => {
            let file = File::open(&list_file)?;
            let reader = std::io::BufReader::new(file);
            let mut files = Vec::new();

            use std::io::BufRead;
            for line in reader.lines() {
                let line = line?;
                let trimmed = line.trim();
                if !trimmed.is_empty() && !trimmed.starts_with('#') {
                    if !Path::new(trimmed).exists() {
                        return Err(io::Error::new(
                            io::ErrorKind::NotFound,
                            format!("FASTA file '{}' not found", trimmed),
                        ));
                    }
                    files.push(trimmed.to_string());
                }
            }

            if files.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("No valid FASTA files found in list file: {}", list_file),
                ));
            }

            Ok(files)
        }
        (true, None) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Either --fasta-files or --fasta-list must be provided",
        )),
        (false, Some(_)) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Cannot specify both --fasta-files and --fasta-list",
        )),
    }
}
