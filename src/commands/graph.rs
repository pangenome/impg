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

/// PGGB-style auto-sparsification heuristic: ln(n)/n * 10.
pub(crate) fn auto_sparsify(n_haps: usize) -> Option<f64> {
    if n_haps <= 1 {
        return None;
    }
    let n = n_haps as f64;
    let frac = n.ln() / n * 10.0;
    if frac >= 1.0 { None } else { Some(frac) }
}

/// Round a value to a "nice" multiple based on its magnitude:
/// ≤500 → multiple of 50, ≤1000 → 100, ≤3000 → 200, >3000 → 500.
fn round_nice(v: u64) -> u64 {
    if v == 0 {
        return 0;
    }
    let step = if v <= 500 {
        50
    } else if v <= 1000 {
        100
    } else if v <= 3000 {
        200
    } else {
        500
    };
    ((v + step / 2) / step * step).max(step)
}

// Import gfasort for graph sorting
use crate::graph::sort_gfa;

// Import from sweepga
use sweepga::aligner::Aligner;
use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};

use crate::commands::create_aligner_adaptive;

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
    /// Input PAF file (skip alignment if provided)
    pub input_paf: Option<String>,
    /// Aligner backend: "wfmash" or "fastga"
    pub aligner: String,

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
    /// Maximum scaffold deviation distance (0 = no limit)
    pub scaffold_dist: u64,
    /// Minimum mapping length to include in filtering
    pub min_mapping_length: u64,
    /// Optional directory to save intermediate files (FASTA, raw PAF, filtered PAF).
    pub debug_dir: Option<String>,
    /// Wfmash mapping sparsification fraction (0.0-1.0). None = keep all.
    pub sparsify: Option<String>,
}

impl Default for GraphBuildConfig {
    fn default() -> Self {
        GraphBuildConfig {
            num_threads: 4,
            frequency_multiplier: 10,
            frequency: None,
            min_alignment_length: 0,
            repeat_max: 0,
            min_repeat_dist: 0,
            min_match_len: 23,
            sparse_factor: 0.0,
            transclose_batch: 10_000_000,
            use_in_memory: true,
            show_progress: true,
            temp_dir: None,
            input_paf: None,
            aligner: "wfmash".to_string(),
            // Filtering options
            no_filter: false,
            num_mappings: "many:many".to_string(),
            scaffold_jump: 50_000,              // 50kb default scaffold gap
            scaffold_mass: 10_000,              // 10kb minimum scaffold length
            scaffold_filter: "many:many".to_string(),
            overlap: 0.95,
            min_identity: 0.0,
            scaffold_dist: 0,      // No deviation limit by default
            min_mapping_length: 0, // No minimum mapping length by default
            debug_dir: None,
            sparsify: None,
        }
    }
}

/// Count sequences and genomes in FASTA files.
/// Returns `(num_sequences, num_genomes)` using PanSN naming convention.
fn count_sequences_and_genomes_in_fasta(fasta_paths: &[String]) -> io::Result<(usize, usize)> {
    crate::commands::count_sequences_and_genomes(fasta_paths)
}

/// Parse filter mode string (e.g., "1:1", "many:many", "5:3") into FilterMode and limits
pub fn parse_filter_mode(s: &str) -> (FilterMode, Option<usize>, Option<usize>) {
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
            info!("[graph::temp] Using temp directory: {}", temp_dir);
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
    let combined_fasta = tempfile::Builder::new().suffix(".fa").tempfile()?;
    let mut total_bases: u64 = 0;
    {
        let mut writer = BufWriter::new(&combined_fasta);
        for path in fasta_files {
            let file = File::open(path)?;
            // Use niffler to auto-detect compression
            let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
                io::Error::other(format!("Failed to open reader for '{}': {}", path, e))
            })?;
            let reader = BufReader::new(reader);
            for line in reader.lines() {
                let line: String = line?;
                if !line.starts_with('>') {
                    total_bases += line.trim().len() as u64;
                }
                writeln!(writer, "{}", line)?;
            }
        }
        writer.flush()?;
    }

    let avg_seq_len = if num_sequences > 0 {
        total_bases / num_sequences as u64
    } else {
        0
    };

    if config.show_progress {
        info!(
            "[graph::combine] {:.3}s Combined FASTA files for alignment (avg seq len: {} bp)",
            start_time.elapsed().as_secs_f64(),
            avg_seq_len
        );
    }

    // Debug: save combined FASTA
    if let Some(ref debug_dir) = config.debug_dir {
        let dst = format!("{}/combined.fa", debug_dir);
        let _ = std::fs::copy(combined_fasta.path(), &dst);
        info!("[graph::debug] Saved combined FASTA to {}", dst);
    }

    // Create FASTA index (.fai) — required by wfmash
    if config.aligner == "wfmash" {
        rust_htslib::faidx::Reader::from_path(combined_fasta.path()).map_err(|e| {
            io::Error::other(format!("Failed to create FASTA index: {e}"))
        })?;
    }

    // 3) Get PAF alignments - either from input file or run aligner
    let paf_temp: tempfile::NamedTempFile = if let Some(ref input_paf) = config.input_paf {
        // Use provided PAF file - copy to temp file to match expected type
        if config.show_progress {
            info!(
                "[graph::align] {:.3}s Using provided PAF file: {}",
                start_time.elapsed().as_secs_f64(),
                input_paf
            );
        }
        let paf_temp = tempfile::Builder::new().suffix(".paf").tempfile()?;
        std::fs::copy(input_paf, paf_temp.path())?;
        paf_temp
    } else {
        // Run alignment
        if config.show_progress {
            info!(
                "[graph::align] {:.3}s Running {} alignment (f={})",
                start_time.elapsed().as_secs_f64(),
                config.aligner,
                kmer_frequency
            );
        }

        // Segment length is adapted automatically by sweepga from avg_seq_len.
        let segment_length = None;

        // Resolve sparsification: "auto" → pggb heuristic, float → use directly
        let sparsify = match config.sparsify.as_deref() {
            Some("auto") => {
                let frac = auto_sparsify(num_genomes);
                if config.show_progress {
                    if let Some(f) = frac {
                        info!(
                            "[graph::align] {:.3}s Auto-sparsification: keeping {:.1}% of mappings ({} genomes)",
                            start_time.elapsed().as_secs_f64(), f * 100.0, num_genomes
                        );
                    }
                }
                frac
            }
            Some(val) => Some(val.parse::<f64>().map_err(|_| {
                io::Error::other(format!("Invalid --sparsify value '{}': expected 'auto' or a float 0.0-1.0", val))
            })?),
            None => None,
        };

        let aligner: Box<dyn Aligner> = create_aligner_adaptive(
            &config.aligner,
            kmer_frequency,
            config.num_threads,
            config.min_alignment_length,
            Some("90".to_string()),
            config.temp_dir.clone(),
            segment_length,
            Some(avg_seq_len),
            sparsify,
            None, // num_mappings: use wfmash default (-n 1)
        )?;

        // Run all-vs-all alignment (query = target = combined FASTA)
        let paf_temp = aligner
            .align_to_temp_paf(combined_fasta.path(), combined_fasta.path())
            .map_err(|e| io::Error::other(format!("{} alignment failed: {}", config.aligner, e)))?;

        if config.show_progress {
            info!(
                "[graph::align] {:.3}s Alignment complete",
                start_time.elapsed().as_secs_f64()
            );
        }
        paf_temp
    };

    // Debug: save raw PAF
    if let Some(ref debug_dir) = config.debug_dir {
        let dst = format!("{}/raw.paf", debug_dir);
        let _ = std::fs::copy(paf_temp.path(), &dst);
        info!("[graph::debug] Saved raw PAF to {}", dst);
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
        // Parse filter modes
        let (mapping_mode, mapping_per_query, mapping_per_target) =
            parse_filter_mode(&config.num_mappings);
        let (scaffold_mode, scaffold_per_query, scaffold_per_target) =
            parse_filter_mode(&config.scaffold_filter);

        // Auto-adapt scaffold parameters to input sequence sizes.
        // Defaults (scaffold_mass=10kb, scaffold_jump=50kb) are tuned for
        // whole-genome alignments. For short sequences (e.g. 1kb excerpts
        // from `impg query -o fasta`), these thresholds would filter out
        // every alignment. Clamp to 80% of the average sequence length.
        let scaffold_mass = if avg_seq_len > 0 {
            round_nice(config.scaffold_mass.min(avg_seq_len * 3 / 5))
        } else {
            config.scaffold_mass
        };
        let scaffold_jump = if avg_seq_len > 0 {
            config.scaffold_jump.min(avg_seq_len * 10)
        } else {
            config.scaffold_jump
        };

        if config.show_progress {
            info!(
                "[graph::filter] {:.3}s Filtering alignments with sweepga (n={}, scaffold_jump={}, scaffold_mass={}, scaffold_filter={})",
                start_time.elapsed().as_secs_f64(),
                config.num_mappings,
                scaffold_jump,
                scaffold_mass,
                config.scaffold_filter
            );
        }

        // Create filter configuration
        let filter_config = FilterConfig {
            chain_gap: 0,
            min_block_length: config.min_mapping_length,
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
            scaffold_gap: scaffold_jump,
            min_scaffold_length: scaffold_mass,
            scaffold_overlap_threshold: 0.5,
            scaffold_max_deviation: config.scaffold_dist,
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

    // Debug: save filtered PAF
    if let Some(ref debug_dir) = config.debug_dir {
        let dst = format!("{}/filtered.paf", debug_dir);
        let _ = std::fs::copy(filtered_paf.path(), &dst);
        info!("[graph::debug] Saved filtered PAF to {}", dst);
    }

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

    let gfa_string = String::from_utf8(gfa_buffer).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid UTF-8 in GFA: {}", e),
        )
    })?;

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

/// Build a pangenome graph from FASTA sequences using the recursive realize engine
///
/// Instead of the seqwish graph induction pipeline, this uses the recursive
/// realize engine (sweepga + POA with lacing) to build the variation graph.
pub fn run_graph_build_realize<W: Write>(
    fasta_files: Vec<String>,
    fasta_list: Option<String>,
    output: &mut W,
    config: &crate::realize::RealizeConfig,
) -> io::Result<()> {
    let start_time = Instant::now();

    // Resolve FASTA files
    let fasta_files = resolve_fasta_files(fasta_files, fasta_list)?;

    if fasta_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No FASTA files specified. Use --fasta-files or --fasta-list",
        ));
    }

    info!(
        "Building graph from {} FASTA file(s) using realize engine",
        fasta_files.len()
    );

    // Read all sequences from FASTA files into memory
    let sequences = read_fasta_sequences(&fasta_files)?;

    info!(
        "[realize] {:.3}s Loaded {} sequences from FASTA files",
        start_time.elapsed().as_secs_f64(),
        sequences.len()
    );

    if sequences.is_empty() {
        output.write_all(b"H\tVN:Z:1.0\n")?;
        return Ok(());
    }

    // Run the realize engine
    let result = crate::realize::realize_from_sequences(&sequences, config)?;

    info!(
        "[realize] {:.3}s Realize complete: {} sequences, max_depth={}, poa_calls={}, sweepga_calls={}, seqwish_calls={}, {}ms",
        start_time.elapsed().as_secs_f64(),
        result.stats.num_sequences,
        result.stats.max_depth_reached,
        result.stats.poa_calls,
        result.stats.sweepga_calls,
        result.stats.seqwish_calls,
        result.stats.total_ms,
    );

    output.write_all(result.gfa.as_bytes())?;

    Ok(())
}

/// Build a POA graph from FASTA sequences.
///
/// Reads all sequences, runs single-pass SPOA, emits GFA, then sorts.
pub fn run_graph_build_poa<W: Write>(
    fasta_files: Vec<String>,
    fasta_list: Option<String>,
    output: &mut W,
    scoring_params: (u8, u8, u8, u8, u8, u8),
    num_threads: usize,
) -> io::Result<()> {
    // Resolve FASTA files
    let fasta_files = resolve_fasta_files(fasta_files, fasta_list)?;

    if fasta_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No FASTA files specified. Use --fasta-files or --fasta-list",
        ));
    }

    info!(
        "Building graph from {} FASTA file(s) using poa engine",
        fasta_files.len()
    );

    let sequences = read_fasta_sequences(&fasta_files)?;

    info!("Loaded {} sequences from FASTA files", sequences.len());

    if sequences.is_empty() {
        output.write_all(b"H\tVN:Z:1.0\n")?;
        return Ok(());
    }

    // Build SPOA graph
    let (mut graph, mut engine) = crate::graph::build_spoa_engine(scoring_params);
    let metadata: Vec<crate::graph::SequenceMetadata> =
        sequences.iter().map(|(_, m)| m.clone()).collect();
    crate::graph::feed_sequences_to_graph(
        &mut engine,
        &mut graph,
        sequences.iter().map(|(s, _)| s.as_str()),
    );

    // Generate GFA, post-process strands, sort
    let sorted = crate::graph::spoa_graph_to_sorted_gfa(graph, &metadata, num_threads)?;
    output.write_all(sorted.as_bytes())?;

    Ok(())
}

/// Build a pangenome graph using the PGGB pipeline: sweepga + seqwish + smoothxg + gfaffix.
///
/// Equivalent to running `--engine seqwish` followed by graph smoothing (smoothxg-style
/// block decomposition + per-block POA) and optional gfaffix normalization.
pub fn run_graph_build_pggb<W: Write>(
    fasta_files: Vec<String>,
    fasta_list: Option<String>,
    output: &mut W,
    config: &GraphBuildConfig,
    target_poa_length: usize,
    max_node_length: usize,
    poa_padding_fraction: f64,
) -> io::Result<()> {
    let start_time = Instant::now();

    let fasta_files = resolve_fasta_files(fasta_files, fasta_list)?;

    if fasta_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No FASTA files specified. Use --fasta-files or --fasta-list",
        ));
    }

    info!(
        "[pggb] Building graph from {} FASTA file(s)",
        fasta_files.len()
    );

    // Count genomes to set n_haps for smoothing
    let (_num_sequences, num_genomes) = count_sequences_and_genomes_in_fasta(&fasta_files)?;
    let n_haps = num_genomes.max(1);

    // Step 1: Run the seqwish pipeline, capturing the raw GFA
    let mut gfa_buffer: Vec<u8> = Vec::new();
    build_graph(&fasta_files, &mut gfa_buffer, config)?;
    let raw_gfa = String::from_utf8(gfa_buffer)
        .map_err(|e| io::Error::other(format!("seqwish GFA is not valid UTF-8: {}", e)))?;

    info!(
        "[pggb] {:.3}s Seqwish done, smoothing (n_haps={}, target_poa_length={})",
        start_time.elapsed().as_secs_f64(),
        n_haps,
        target_poa_length
    );

    // Step 2: Smooth
    let smooth_config = crate::smooth::SmoothConfig {
        n_haps,
        target_poa_length,
        max_block_weight: target_poa_length * n_haps,
        max_poa_length: 2 * target_poa_length,
        max_node_length,
        poa_padding_fraction,
        num_threads: config.num_threads,
        temp_dir: config.temp_dir.clone(),
        ..crate::smooth::SmoothConfig::new(n_haps)
    };
    let smoothed = crate::smooth::smooth_gfa(&raw_gfa, &smooth_config)?;

    info!(
        "[pggb] {:.3}s Smoothing done",
        start_time.elapsed().as_secs_f64()
    );

    // Step 3: gfaffix normalization (optional, graceful fallback if binary not found)
    let normalized = crate::graph::run_gfaffix(&smoothed, config.num_threads)
        .unwrap_or_else(|e| {
            log::warn!("[pggb] gfaffix not available ({}), skipping normalization", e);
            smoothed
        });

    // Step 4: Final sort
    let sorted = sort_gfa(&normalized, config.num_threads)?;

    info!(
        "[pggb] {:.3}s Done",
        start_time.elapsed().as_secs_f64()
    );

    output.write_all(sorted.as_bytes())?;
    Ok(())
}

/// Read sequences from FASTA files into (sequence, metadata) pairs for the realize engine.
fn read_fasta_sequences(
    fasta_paths: &[String],
) -> io::Result<Vec<(String, crate::graph::SequenceMetadata)>> {
    let mut sequences = Vec::new();

    for path in fasta_paths {
        let file = File::open(path)?;
        // Use niffler to auto-detect compression
        let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
            io::Error::other(format!("Failed to open reader for '{}': {}", path, e))
        })?;
        let reader = BufReader::new(reader);

        let mut current_name = String::new();
        let mut current_seq = String::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                // Save previous sequence if any
                if !current_name.is_empty() && !current_seq.is_empty() {
                    let seq_len = current_seq.len();
                    sequences.push((
                        std::mem::take(&mut current_seq),
                        metadata_from_fasta_header(&current_name, seq_len),
                    ));
                }
                // Start new sequence
                current_name = line[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                current_seq.clear();
            } else {
                current_seq.push_str(line.trim());
            }
        }

        // Don't forget the last sequence
        if !current_name.is_empty() && !current_seq.is_empty() {
            let seq_len = current_seq.len();
            sequences.push((current_seq, metadata_from_fasta_header(&current_name, seq_len)));
        }
    }

    // Sort by length descending (expected by realize engine)
    sequences.sort_by(|a, b| b.0.len().cmp(&a.0.len()));

    Ok(sequences)
}

/// Parse a FASTA header into SequenceMetadata.
///
/// If the header contains a `name:start-end` suffix (e.g. from a prior
/// `impg query -o fasta`), parse the coordinates so that `path_name()`
/// reproduces the original name without double-encoding.
/// Otherwise, treat the whole header as the name with start=0.
fn metadata_from_fasta_header(header: &str, seq_len: usize) -> crate::graph::SequenceMetadata {
    // Try to parse trailing :start-end
    if let Some(last_colon) = header.rfind(':') {
        let (key, range_str) = header.split_at(last_colon);
        let range_str = &range_str[1..]; // skip ':'
        if let Some((start_str, end_str)) = range_str.split_once('-') {
            if let (Ok(start), Ok(end)) = (start_str.parse::<i32>(), end_str.parse::<i32>()) {
                if end > start {
                    return crate::graph::SequenceMetadata {
                        name: key.to_string(),
                        start,
                        size: seq_len as i32,
                        strand: '+',
                        total_length: (start + seq_len as i32) as usize,
                    };
                }
            }
        }
    }
    // Fallback: treat whole header as name
    crate::graph::SequenceMetadata {
        name: header.to_string(),
        start: 0,
        size: seq_len as i32,
        strand: '+',
        total_length: seq_len,
    }
}

/// Resolve FASTA files from either direct list or file list
fn resolve_fasta_files(
    fasta_files: Vec<String>,
    fasta_list: Option<String>,
) -> io::Result<Vec<String>> {
    crate::commands::resolve_file_list(fasta_files, fasta_list, "FASTA")
}
