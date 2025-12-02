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

// Import from sweepga
use sweepga::fastga_integration::FastGAIntegration;

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
        }
    }
}

/// Count the number of unique sequences (genomes/haplotypes) in FASTA files
fn count_sequences_in_fasta(fasta_paths: &[String]) -> io::Result<usize> {
    let mut count = 0;
    for path in fasta_paths {
        let file = File::open(path)?;
        // Use niffler to auto-detect compression
        let (reader, _format) = niffler::get_reader(Box::new(file))
            .map_err(|e| io::Error::other(format!("Failed to open reader for '{}': {}", path, e)))?;
        let reader = BufReader::new(reader);

        for line in reader.lines() {
            let line: String = line?;
            if line.starts_with('>') {
                count += 1;
            }
        }
    }
    Ok(count)
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

    // Set up temp directory for seqwish
    if let Some(ref temp_dir) = config.temp_dir {
        seqwish::tempfile::set_dir(temp_dir);
    }

    // 1) Count sequences and determine k-mer frequency
    if config.show_progress {
        info!(
            "[graph::count] {:.3}s Counting sequences in {} FASTA file(s)",
            start_time.elapsed().as_secs_f64(),
            fasta_files.len()
        );
    }

    let num_sequences = count_sequences_in_fasta(fasta_files)?;
    let kmer_frequency = config
        .frequency
        .unwrap_or(num_sequences * config.frequency_multiplier);

    if config.show_progress {
        info!(
            "[graph::count] {:.3}s Found {} sequences, using k-mer frequency {}",
            start_time.elapsed().as_secs_f64(),
            num_sequences,
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
    );

    // Run all-vs-all alignment (query = target = combined FASTA)
    let paf_temp = fastga
        .align_to_temp_paf(combined_fasta.path(), combined_fasta.path())
        .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("FastGA alignment failed: {}", e)))?;

    if config.show_progress {
        info!(
            "[graph::align] {:.3}s Alignment complete",
            start_time.elapsed().as_secs_f64()
        );
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
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
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
        paf_temp.path().to_str().unwrap(),
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
            .map_err(|_| io::Error::new(io::ErrorKind::Other, "Failed to unwrap alignment tree Arc"))?
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

    // 9) Emit GFA output
    if config.show_progress {
        info!(
            "[graph::gfa] {:.3}s Writing GFA output",
            start_time.elapsed().as_secs_f64()
        );
    }

    emit_gfa(
        output,
        graph_length,
        seq_v_file.to_str().unwrap(),
        Arc::clone(&node_iitree),
        Arc::clone(&path_iitree),
        &seq_id_cbv,
        Arc::clone(&seqidx),
        link_set.links(),
        config.num_threads,
    )?;

    if config.show_progress {
        info!(
            "[graph::gfa] {:.3}s Graph building complete",
            start_time.elapsed().as_secs_f64()
        );
    }

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
