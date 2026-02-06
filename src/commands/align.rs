// Alignment command for setting up pangenome alignments with sweepga
//
// This module provides functionality to generate alignment pairs from input sequences,
// with various sparsification strategies to reduce the number of alignments needed.

use log::info;
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

// Import sweepga for direct alignment execution
use sweepga::fastga_integration::FastGAIntegration;
use sweepga::paf_filter::{FilterConfig, PafFilter, ScoringFunction};

/// Sparsification strategy for pair selection
#[derive(Clone, Debug)]
pub enum SparsificationStrategy {
    /// Align all pairs - O(n²) complexity
    None,
    /// Random subsampling with given fraction (0.0-1.0)
    Random(f64),
    /// Random subsampling with giant component guarantee (Erdős-Rényi model)
    /// Parameter is the probability of connectivity (e.g., 0.99)
    Connectivity(f64),
    /// Tree-based sampling combining k-nearest, k-farthest (stranger-joining), and random
    /// Parameters: (k_nearest, k_farthest, random_fraction, kmer_size)
    TreeSampling {
        k_nearest: usize,
        k_farthest: usize,
        random_fraction: f64,
        kmer_size: usize,
    },
}

impl SparsificationStrategy {
    /// Parse sparsification strategy from string
    /// Formats:
    /// - "none" or "all" - all pairs
    /// - "random:0.5" - random 50%
    /// - "giant:0.99" or "connectivity:0.99" - giant component with 99% probability
    /// - "tree:2:1:0.1" or "tree:2:1:0.1:15" - k-nearest:k-farthest:random_frac[:kmer_size]
    pub fn parse(s: &str) -> io::Result<Self> {
        let s = s.to_lowercase();
        let parts: Vec<&str> = s.split(':').collect();

        match parts[0] {
            "none" | "all" => Ok(SparsificationStrategy::None),
            "random" => {
                if parts.len() != 2 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Random format: random:<fraction> (e.g., random:0.5)",
                    ));
                }
                let fraction: f64 = parts[1].parse().map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidInput, "Invalid fraction")
                })?;
                if !(0.0..=1.0).contains(&fraction) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Fraction must be between 0.0 and 1.0",
                    ));
                }
                Ok(SparsificationStrategy::Random(fraction))
            }
            "giant" | "connectivity" => {
                if parts.len() != 2 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Connectivity format: giant:<probability> (e.g., giant:0.99)",
                    ));
                }
                let prob: f64 = parts[1].parse().map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidInput, "Invalid probability")
                })?;
                if !(0.0..1.0).contains(&prob) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Probability must be between 0.0 and 1.0 (exclusive)",
                    ));
                }
                Ok(SparsificationStrategy::Connectivity(prob))
            }
            "tree" | "knn" => {
                if parts.len() < 4 || parts.len() > 5 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Tree format: tree:<k_nearest>:<k_farthest>:<random_frac>[:<kmer_size>]",
                    ));
                }
                let k_nearest: usize = parts[1].parse().map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidInput, "Invalid k_nearest")
                })?;
                let k_farthest: usize = parts[2].parse().map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidInput, "Invalid k_farthest")
                })?;
                let random_fraction: f64 = parts[3].parse().map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidInput, "Invalid random_fraction")
                })?;
                let kmer_size: usize = if parts.len() == 5 {
                    parts[4].parse().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidInput, "Invalid kmer_size")
                    })?
                } else {
                    15 // default
                };

                if !(0.0..=1.0).contains(&random_fraction) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Random fraction must be between 0.0 and 1.0",
                    ));
                }
                if !(3..=31).contains(&kmer_size) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "K-mer size must be between 3 and 31",
                    ));
                }

                Ok(SparsificationStrategy::TreeSampling {
                    k_nearest,
                    k_farthest,
                    random_fraction,
                    kmer_size,
                })
            }
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Unknown sparsification strategy: '{}'. Valid: none, random:<frac>, giant:<prob>, tree:<k_near>:<k_far>:<rand_frac>",
                    parts[0]
                ),
            )),
        }
    }

    /// Get a description of the strategy
    pub fn description(&self) -> String {
        match self {
            SparsificationStrategy::None => "all pairs (no sparsification)".to_string(),
            SparsificationStrategy::Random(f) => format!("random {:.1}%", f * 100.0),
            SparsificationStrategy::Connectivity(p) => {
                format!("giant component (p={:.2})", p)
            }
            SparsificationStrategy::TreeSampling {
                k_nearest,
                k_farthest,
                random_fraction,
                kmer_size,
            } => {
                format!(
                    "tree sampling (k_near={}, k_far={}, rand={:.1}%, k={})",
                    k_nearest,
                    k_farthest,
                    random_fraction * 100.0,
                    kmer_size
                )
            }
        }
    }
}

/// Configuration for alignment command
pub struct AlignConfig {
    pub num_threads: usize,
    pub sparsification: SparsificationStrategy,
    pub frequency_multiplier: usize,
    pub frequency: Option<usize>,
    pub min_alignment_length: u64,
    pub output_format: AlignOutputFormat,
    pub show_progress: bool,
    // Sweepga filtering options
    pub no_filter: bool,
    pub num_mappings: String,
    pub scaffold_jump: u64,
    pub scaffold_mass: u64,
    pub scaffold_filter: String,
    pub overlap: f64,
    pub min_identity: f64,
    pub scaffold_dist: u64,
    pub min_mapping_length: u64,
}

impl Default for AlignConfig {
    fn default() -> Self {
        AlignConfig {
            num_threads: 4,
            sparsification: SparsificationStrategy::Connectivity(0.99),
            frequency_multiplier: 10,
            frequency: None,
            min_alignment_length: 100,
            output_format: AlignOutputFormat::Paf,
            show_progress: true,
            // Filtering defaults - match graph command
            no_filter: false,
            num_mappings: "1:1".to_string(),
            scaffold_jump: 50_000,
            scaffold_mass: 10_000,
            scaffold_filter: "1:1".to_string(),
            overlap: 0.95,
            min_identity: 0.0,
            scaffold_dist: 0,
            min_mapping_length: 0,
        }
    }
}

/// Output format for alignments
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AlignOutputFormat {
    /// PAF format (standard)
    Paf,
    /// 1aln format (compact, for impg)
    OneAln,
    /// Job list (commands only, for cluster execution)
    JobList,
}

/// Sequence info for distance computation
struct SequenceInfo {
    name: String,
    path: String,
    sketch: Vec<u64>,
}

/// Compute MinHash sketch for a sequence
fn sketch_sequence(sequence: &[u8], k: usize, sketch_size: usize) -> Vec<u64> {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let mut hashes = Vec::new();

    if sequence.len() < k {
        return vec![];
    }

    for i in 0..=(sequence.len() - k) {
        let kmer = &sequence[i..i + k];

        // Skip k-mers with non-ACGT characters
        if kmer
            .iter()
            .any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
        {
            continue;
        }

        // Compute canonical k-mer (lexicographically smaller of forward/reverse)
        let rc_kmer: Vec<u8> = kmer
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                _ => b'N',
            })
            .collect();

        let canonical = if kmer < rc_kmer.as_slice() {
            kmer
        } else {
            &rc_kmer
        };

        let mut hasher = DefaultHasher::new();
        canonical.hash(&mut hasher);
        hashes.push(hasher.finish());
    }

    // Sort and take smallest (MinHash)
    hashes.sort_unstable();
    hashes.truncate(sketch_size);
    hashes
}

/// Compute Mash distance between two sketches
fn mash_distance(sketch1: &[u64], sketch2: &[u64], k: usize) -> f64 {
    if sketch1.is_empty() || sketch2.is_empty() {
        return 1.0;
    }

    // Compute Jaccard index using sorted merge
    let mut i = 0;
    let mut j = 0;
    let mut intersection = 0;
    let mut union_size = 0;

    while i < sketch1.len() && j < sketch2.len() {
        if sketch1[i] == sketch2[j] {
            intersection += 1;
            union_size += 1;
            i += 1;
            j += 1;
        } else if sketch1[i] < sketch2[j] {
            union_size += 1;
            i += 1;
        } else {
            union_size += 1;
            j += 1;
        }
    }
    union_size += sketch1.len() - i;
    union_size += sketch2.len() - j;

    if union_size == 0 {
        return 1.0;
    }

    let jaccard = intersection as f64 / union_size as f64;

    if jaccard == 0.0 {
        return 1.0;
    }

    // Mash distance formula
    let mash = -1.0 / (k as f64) * (2.0 * jaccard / (1.0 + jaccard)).ln();
    mash.clamp(0.0, 1.0)
}

/// Compute probability for giant component (Erdős-Rényi model)
fn compute_connectivity_probability(n: usize, connectivity_prob: f64) -> f64 {
    if n <= 1 {
        return 1.0;
    }

    // Special cases for small n
    match n {
        2 => return 1.0,
        3 => return 0.8,
        4 => return 0.7,
        5 => return 0.6,
        6..=10 => return 0.5,
        _ => {}
    }

    // p = (log n + c) / n where c = -log(-log(connectivity_prob))
    let c = -(-connectivity_prob.ln()).ln();
    let p = ((n as f64).ln() + c) / (n as f64);

    p.clamp(0.001, 1.0)
}

/// Generate alignment pairs from named in-memory sequences.
///
/// This is the version used by the realize engine's `sweepga_align`, where
/// sequences are provided as `(name, bytes)` pairs rather than `SequenceInfo`.
/// For `TreeSampling`, mash sketches are computed on the fly.
pub fn generate_pairs_for_sequences(
    sequences: &[(String, &[u8])],
    strategy: &SparsificationStrategy,
) -> Vec<(usize, usize)> {
    let n = sequences.len();
    if n <= 1 {
        return vec![];
    }

    match strategy {
        SparsificationStrategy::None => {
            (0..n)
                .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
                .collect()
        }

        SparsificationStrategy::Random(fraction) => {
            use std::collections::hash_map::DefaultHasher;
            use std::hash::{Hash, Hasher};

            let mut pairs = Vec::new();
            for i in 0..n {
                for j in (i + 1)..n {
                    let mut hasher = DefaultHasher::new();
                    format!("{}:{}", sequences[i].0, sequences[j].0).hash(&mut hasher);
                    let hash = hasher.finish();
                    if (hash as f64 / u64::MAX as f64) < *fraction {
                        pairs.push((i, j));
                    }
                }
            }
            pairs
        }

        SparsificationStrategy::Connectivity(prob) => {
            let keep_fraction = compute_connectivity_probability(n, *prob);
            generate_pairs_for_sequences(
                sequences,
                &SparsificationStrategy::Random(keep_fraction),
            )
        }

        SparsificationStrategy::TreeSampling {
            k_nearest,
            k_farthest,
            random_fraction,
            kmer_size,
        } => {
            // Compute mash sketches on the fly
            let sketches: Vec<Vec<u64>> = sequences
                .iter()
                .map(|(_, seq)| sketch_sequence(seq, *kmer_size, 1000))
                .collect();

            let distances: Vec<Vec<f64>> = (0..n)
                .into_par_iter()
                .map(|i| {
                    (0..n)
                        .map(|j| {
                            if i == j {
                                0.0
                            } else {
                                mash_distance(&sketches[i], &sketches[j], *kmer_size)
                            }
                        })
                        .collect()
                })
                .collect();

            let mut pairs = HashSet::new();

            if *k_nearest > 0 {
                for (i, dist_row) in distances.iter().enumerate() {
                    let mut neighbors: Vec<(usize, f64)> = (0..n)
                        .filter(|&j| i != j)
                        .map(|j| (j, dist_row[j]))
                        .collect();
                    neighbors.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
                    for (j, _) in neighbors.iter().take(*k_nearest) {
                        let pair = if i < *j { (i, *j) } else { (*j, i) };
                        pairs.insert(pair);
                    }
                }
            }

            if *k_farthest > 0 {
                for (i, dist_row) in distances.iter().enumerate() {
                    let mut neighbors: Vec<(usize, f64)> = (0..n)
                        .filter(|&j| i != j)
                        .map(|j| (j, dist_row[j]))
                        .collect();
                    neighbors.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
                    for (j, _) in neighbors.iter().take(*k_farthest) {
                        let pair = if i < *j { (i, *j) } else { (*j, i) };
                        pairs.insert(pair);
                    }
                }
            }

            if *random_fraction > 0.0 {
                let random_pairs = generate_pairs_for_sequences(
                    sequences,
                    &SparsificationStrategy::Random(*random_fraction),
                );
                for pair in random_pairs {
                    pairs.insert(pair);
                }
            }

            let mut result: Vec<_> = pairs.into_iter().collect();
            result.sort_unstable();
            result
        }
    }
}

/// Generate pairs using specified sparsification strategy
fn generate_pairs(
    sequences: &[SequenceInfo],
    strategy: &SparsificationStrategy,
) -> Vec<(usize, usize)> {
    let n = sequences.len();
    if n <= 1 {
        return vec![];
    }

    match strategy {
        SparsificationStrategy::None => {
            // All pairs
            (0..n)
                .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
                .collect()
        }

        SparsificationStrategy::Random(fraction) => {
            use std::collections::hash_map::DefaultHasher;
            use std::hash::{Hash, Hasher};

            let mut pairs = Vec::new();
            for i in 0..n {
                for j in (i + 1)..n {
                    // Deterministic hash based on sequence names
                    let mut hasher = DefaultHasher::new();
                    format!("{}:{}", sequences[i].name, sequences[j].name).hash(&mut hasher);
                    let hash = hasher.finish();

                    if (hash as f64 / u64::MAX as f64) < *fraction {
                        pairs.push((i, j));
                    }
                }
            }
            pairs
        }

        SparsificationStrategy::Connectivity(prob) => {
            let keep_fraction = compute_connectivity_probability(n, *prob);
            generate_pairs(sequences, &SparsificationStrategy::Random(keep_fraction))
        }

        SparsificationStrategy::TreeSampling {
            k_nearest,
            k_farthest,
            random_fraction,
            kmer_size: _,
        } => {
            // Compute distance matrix
            let distances: Vec<Vec<f64>> = (0..n)
                .into_par_iter()
                .map(|i| {
                    (0..n)
                        .map(|j| {
                            if i == j {
                                0.0
                            } else {
                                mash_distance(&sequences[i].sketch, &sequences[j].sketch, 15)
                            }
                        })
                        .collect()
                })
                .collect();

            let mut pairs = HashSet::new();

            // K-nearest neighbors for each sequence
            if *k_nearest > 0 {
                for (i, dist_row) in distances.iter().enumerate() {
                    let mut neighbors: Vec<(usize, f64)> = (0..n)
                        .filter(|&j| i != j)
                        .map(|j| (j, dist_row[j]))
                        .collect();
                    neighbors.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

                    for (j, _) in neighbors.iter().take(*k_nearest) {
                        let pair = if i < *j { (i, *j) } else { (*j, i) };
                        pairs.insert(pair);
                    }
                }
            }

            // K-farthest neighbors (stranger-joining)
            if *k_farthest > 0 {
                for (i, dist_row) in distances.iter().enumerate() {
                    let mut neighbors: Vec<(usize, f64)> = (0..n)
                        .filter(|&j| i != j)
                        .map(|j| (j, dist_row[j]))
                        .collect();
                    neighbors.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap()); // Descending

                    for (j, _) in neighbors.iter().take(*k_farthest) {
                        let pair = if i < *j { (i, *j) } else { (*j, i) };
                        pairs.insert(pair);
                    }
                }
            }

            // Random pairs
            if *random_fraction > 0.0 {
                let random_pairs =
                    generate_pairs(sequences, &SparsificationStrategy::Random(*random_fraction));
                for pair in random_pairs {
                    pairs.insert(pair);
                }
            }

            let mut result: Vec<_> = pairs.into_iter().collect();
            result.sort_unstable();
            result
        }
    }
}

/// Load sequences from FASTA files and compute sketches
fn load_sequences(
    fasta_files: &[String],
    kmer_size: usize,
    sketch_size: usize,
    show_progress: bool,
) -> io::Result<Vec<SequenceInfo>> {
    let sequences = Mutex::new(Vec::new());

    fasta_files
        .par_iter()
        .try_for_each(|path| -> io::Result<()> {
            let file = File::open(path)?;
            let (reader, _format) = niffler::get_reader(Box::new(file))
                .map_err(|e| io::Error::other(format!("Failed to open {}: {}", path, e)))?;
            let reader = BufReader::new(reader);

            let mut current_name: Option<String> = None;
            let mut current_seq = Vec::new();

            for line in reader.lines() {
                let line = line?;
                if let Some(header) = line.strip_prefix('>') {
                    // Process previous sequence
                    if let Some(name) = current_name.take() {
                        let sketch = sketch_sequence(&current_seq, kmer_size, sketch_size);
                        sequences.lock().unwrap().push(SequenceInfo {
                            name,
                            path: path.clone(),
                            sketch,
                        });
                    }
                    current_name = Some(header.split_whitespace().next().unwrap_or("").to_string());
                    current_seq.clear();
                } else {
                    current_seq.extend(line.trim().as_bytes());
                }
            }

            // Don't forget last sequence
            if let Some(name) = current_name {
                let sketch = sketch_sequence(&current_seq, kmer_size, sketch_size);
                sequences.lock().unwrap().push(SequenceInfo {
                    name,
                    path: path.clone(),
                    sketch,
                });
            }

            Ok(())
        })?;

    let mut result = sequences.into_inner().unwrap();

    if show_progress {
        info!(
            "Loaded {} sequences from {} files",
            result.len(),
            fasta_files.len()
        );
    }

    // Sort by name for reproducibility
    result.sort_by(|a, b| a.name.cmp(&b.name));

    Ok(result)
}

/// Estimate pair count for a given strategy and sequence count
#[allow(dead_code)]
fn estimate_pair_count(n: usize, strategy: &SparsificationStrategy) -> usize {
    if n <= 1 {
        return 0;
    }

    let total_pairs = n * (n - 1) / 2;

    match strategy {
        SparsificationStrategy::None => total_pairs,
        SparsificationStrategy::Random(fraction) => ((total_pairs as f64) * fraction) as usize,
        SparsificationStrategy::Connectivity(prob) => {
            let keep_fraction = compute_connectivity_probability(n, *prob);
            ((total_pairs as f64) * keep_fraction) as usize
        }
        SparsificationStrategy::TreeSampling {
            k_nearest,
            k_farthest,
            random_fraction,
            ..
        } => {
            // Each sequence contributes k_nearest + k_farthest edges (with dedup)
            let tree_pairs = n * (k_nearest + k_farthest);
            let random_pairs = ((total_pairs as f64) * random_fraction) as usize;
            // Rough estimate (edges can overlap)
            (tree_pairs / 2 + random_pairs).min(total_pairs)
        }
    }
}

/// Write job list for cluster execution
fn write_job_list<W: Write>(
    pairs: &[(usize, usize)],
    sequences: &[SequenceInfo],
    output_dir: &str,
    writer: &mut W,
    config: &AlignConfig,
) -> io::Result<()> {
    let freq_arg = if let Some(f) = config.frequency {
        format!("-f{}", f)
    } else {
        format!("-f{}", sequences.len() * config.frequency_multiplier)
    };

    for (i, j) in pairs {
        let seq_i = &sequences[*i];
        let seq_j = &sequences[*j];

        // Create output filename
        let output_name = format!("{}_vs_{}.paf", seq_i.name, seq_j.name);
        let output_path = format!("{}/{}", output_dir, output_name);

        // Write command
        writeln!(
            writer,
            "FastGA {} -T{} -l{} {} {} > {}",
            freq_arg,
            config.num_threads,
            config.min_alignment_length,
            seq_i.path,
            seq_j.path,
            output_path
        )?;
    }

    Ok(())
}

/// Collect unique file pairs from sequence pairs.
/// Returns deduplicated (path_i, path_j) pairs where path_i < path_j lexicographically,
/// plus self-pairs where path_i == path_j (for multi-sequence FASTA files).
fn collect_file_pairs(
    pairs: &[(usize, usize)],
    sequences: &[SequenceInfo],
) -> Vec<(String, String)> {
    let mut file_pairs = HashSet::new();
    for &(i, j) in pairs {
        let path_i = &sequences[i].path;
        let path_j = &sequences[j].path;
        // Normalize ordering for deduplication
        let pair = if path_i <= path_j {
            (path_i.clone(), path_j.clone())
        } else {
            (path_j.clone(), path_i.clone())
        };
        file_pairs.insert(pair);
    }
    let mut result: Vec<_> = file_pairs.into_iter().collect();
    result.sort();
    result
}

/// Count sequences and genomes across FASTA files (for k-mer frequency calculation).
/// Returns (num_sequences, num_genomes) using PanSN naming convention.
fn count_sequences_and_genomes(fasta_files: &[String]) -> io::Result<(usize, usize)> {
    let mut seq_count = 0;
    let mut genome_prefixes: HashSet<String> = HashSet::new();

    for path in fasta_files {
        let file = File::open(path)?;
        let (reader, _format) = niffler::get_reader(Box::new(file)).map_err(|e| {
            io::Error::other(format!("Failed to open reader for '{}': {}", path, e))
        })?;
        let reader = BufReader::new(reader);

        for line in reader.lines() {
            let line: String = line?;
            if line.starts_with('>') {
                seq_count += 1;
                let name = line[1..].split_whitespace().next().unwrap_or("");
                let parts: Vec<&str> = name.split('#').collect();
                let prefix = if parts.len() >= 2 {
                    format!("{}#{}", parts[0], parts[1])
                } else {
                    name.to_string()
                };
                genome_prefixes.insert(prefix);
            }
        }
    }

    let genome_count = genome_prefixes.len().max(1);
    Ok((seq_count, genome_count))
}

/// Run sweepga alignments for selected pairs and write output.
///
/// For each unique file pair, runs FastGA alignment, applies filtering (unless --no-filter),
/// and writes results to the output directory. For PAF format, produces a combined
/// output.paf file. For 1aln format, produces individual .1aln files per pair.
fn run_alignments(
    pairs: &[(usize, usize)],
    sequences: &[SequenceInfo],
    output_dir: &str,
    config: &AlignConfig,
) -> io::Result<()> {
    let start_time = Instant::now();

    // Collect unique FASTA file paths for genome counting
    let unique_files: Vec<String> = {
        let mut files: HashSet<String> = HashSet::new();
        for seq in sequences {
            files.insert(seq.path.clone());
        }
        files.into_iter().collect()
    };

    // Determine k-mer frequency
    let (_num_sequences, num_genomes) = count_sequences_and_genomes(&unique_files)?;
    let kmer_frequency = config
        .frequency
        .unwrap_or(num_genomes * config.frequency_multiplier);

    if config.show_progress {
        info!(
            "[align] {:.3}s {} genomes, using k-mer frequency {}",
            start_time.elapsed().as_secs_f64(),
            num_genomes,
            kmer_frequency
        );
    }

    // Collect unique file pairs
    let file_pairs = collect_file_pairs(pairs, sequences);

    if config.show_progress {
        info!(
            "[align] {:.3}s Running {} pairwise alignments ({} unique file pairs)",
            start_time.elapsed().as_secs_f64(),
            pairs.len(),
            file_pairs.len()
        );
    }

    // Create FastGA integration
    let fastga = FastGAIntegration::new(
        Some(kmer_frequency),
        config.num_threads,
        config.min_alignment_length,
        None, // temp_dir uses system default
    );

    // Run alignments for each unique file pair
    let paf_output_path = format!("{}/alignments.paf", output_dir);
    let mut combined_writer = BufWriter::new(File::create(&paf_output_path)?);

    for (pair_idx, (query_path, target_path)) in file_pairs.iter().enumerate() {
        if config.show_progress {
            info!(
                "[align] {:.3}s Aligning pair {}/{}: {} vs {}",
                start_time.elapsed().as_secs_f64(),
                pair_idx + 1,
                file_pairs.len(),
                Path::new(query_path)
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy(),
                Path::new(target_path)
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy(),
            );
        }

        match config.output_format {
            AlignOutputFormat::OneAln => {
                // For 1aln format, use align_to_temp_1aln and copy to output
                let temp_aln = fastga
                    .align_to_temp_1aln(Path::new(query_path), Path::new(target_path))
                    .map_err(|e| {
                        io::Error::other(format!(
                            "FastGA alignment failed for {} vs {}: {}",
                            query_path, target_path, e
                        ))
                    })?;

                let query_stem = Path::new(query_path)
                    .file_stem()
                    .unwrap_or_default()
                    .to_string_lossy();
                let target_stem = Path::new(target_path)
                    .file_stem()
                    .unwrap_or_default()
                    .to_string_lossy();
                let output_name = format!("{}_vs_{}.1aln", query_stem, target_stem);
                let output_path = format!("{}/{}", output_dir, output_name);
                std::fs::copy(temp_aln.path(), &output_path)?;
            }
            AlignOutputFormat::Paf => {
                // Run alignment to temp PAF
                let paf_temp = fastga
                    .align_to_temp_paf(Path::new(query_path), Path::new(target_path))
                    .map_err(|e| {
                        io::Error::other(format!(
                            "FastGA alignment failed for {} vs {}: {}",
                            query_path, target_path, e
                        ))
                    })?;

                // Apply filtering unless --no-filter
                let result_paf = if config.no_filter {
                    paf_temp
                } else {
                    apply_paf_filter(paf_temp, config)?
                };

                // Append to combined output
                let paf_data = std::fs::read(result_paf.path())?;
                combined_writer.write_all(&paf_data)?;
            }
            AlignOutputFormat::JobList => unreachable!(),
        }
    }

    combined_writer.flush()?;

    if config.show_progress {
        match config.output_format {
            AlignOutputFormat::Paf => {
                info!(
                    "[align] {:.3}s Alignments complete, wrote {}",
                    start_time.elapsed().as_secs_f64(),
                    paf_output_path,
                );
            }
            AlignOutputFormat::OneAln => {
                info!(
                    "[align] {:.3}s Alignments complete, wrote {} .1aln files to {}",
                    start_time.elapsed().as_secs_f64(),
                    file_pairs.len(),
                    output_dir,
                );
            }
            AlignOutputFormat::JobList => unreachable!(),
        }
    }

    Ok(())
}

/// Apply sweepga filtering to a PAF temp file, returning a new filtered temp file.
fn apply_paf_filter(
    paf_temp: tempfile::NamedTempFile,
    config: &AlignConfig,
) -> io::Result<tempfile::NamedTempFile> {
    use super::graph::parse_filter_mode;

    let (mapping_mode, mapping_per_query, mapping_per_target) =
        parse_filter_mode(&config.num_mappings);
    let (scaffold_mode, scaffold_per_query, scaffold_per_target) =
        parse_filter_mode(&config.scaffold_filter);

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
        scaffold_gap: config.scaffold_jump,
        min_scaffold_length: config.scaffold_mass,
        scaffold_overlap_threshold: 0.5,
        scaffold_max_deviation: config.scaffold_dist,
        prefix_delimiter: '#',
        skip_prefix: false,
        scoring_function: ScoringFunction::LogLengthIdentity,
        min_identity: config.min_identity,
        min_scaffold_identity: config.min_identity,
    };

    let filtered_paf_file = tempfile::Builder::new()
        .suffix(".filtered.paf")
        .tempfile()?;

    let filter = PafFilter::new(filter_config).with_keep_self(false);
    filter
        .filter_paf(paf_temp.path(), filtered_paf_file.path())
        .map_err(|e| io::Error::other(format!("PAF filtering failed: {}", e)))?;

    Ok(filtered_paf_file)
}

/// Configuration for in-memory sweepga alignment (used by the realize engine).
pub struct SweepgaAlignConfig {
    /// Number of threads for alignment
    pub num_threads: usize,
    /// K-mer frequency for FastGA
    pub kmer_frequency: usize,
    /// Minimum alignment length
    pub min_alignment_length: u64,
    /// Whether to skip filtering
    pub no_filter: bool,
    /// Filter: n:m-best mappings
    pub num_mappings: String,
    /// Filter: scaffold jump distance
    pub scaffold_jump: u64,
    /// Filter: minimum scaffold chain length
    pub scaffold_mass: u64,
    /// Filter: scaffold filter mode
    pub scaffold_filter: String,
    /// Filter: max overlap ratio
    pub overlap: f64,
    /// Filter: minimum identity
    pub min_identity: f64,
    /// Filter: max scaffold deviation
    pub scaffold_dist: u64,
    /// Filter: minimum mapping length
    pub min_mapping_length: u64,
    /// Optional temp directory
    pub temp_dir: Option<String>,
    /// Sparsification strategy for pair selection.
    /// When not None, only alignments between selected pairs are retained.
    pub sparsification: SparsificationStrategy,
}

impl Default for SweepgaAlignConfig {
    fn default() -> Self {
        SweepgaAlignConfig {
            num_threads: 4,
            kmer_frequency: 10,
            min_alignment_length: 100,
            no_filter: false,
            num_mappings: "1:1".to_string(),
            scaffold_jump: 50_000,
            scaffold_mass: 10_000,
            scaffold_filter: "1:1".to_string(),
            overlap: 0.95,
            min_identity: 0.0,
            scaffold_dist: 0,
            min_mapping_length: 0,
            temp_dir: None,
            sparsification: SparsificationStrategy::None,
        }
    }
}

/// Run sweepga alignment on in-memory sequences with optional sparsification.
///
/// When sparsification is `None`, runs all-vs-all alignment. Otherwise,
/// generates selected pairs using the sparsification strategy and runs
/// pairwise alignments only for those pairs.
///
/// # Arguments
/// * `sequences` - Named sequences: `(name, sequence_bytes)` pairs
/// * `config` - Alignment, filtering, and sparsification configuration
///
/// # Returns
/// A `NamedTempFile` containing the (optionally filtered) PAF alignments.
pub fn sweepga_align(
    sequences: &[(String, &[u8])],
    config: &SweepgaAlignConfig,
) -> io::Result<tempfile::NamedTempFile> {
    if sequences.len() < 2 {
        // Nothing to align — return empty PAF
        let temp = tempfile::Builder::new().suffix(".paf").tempfile()?;
        return Ok(temp);
    }

    // Set up temp directory
    if let Some(ref temp_dir) = config.temp_dir {
        std::env::set_var("TMPDIR", temp_dir);
    }

    // Generate pairs based on sparsification strategy
    let pairs = generate_pairs_for_sequences(sequences, &config.sparsification);
    let total_possible = sequences.len() * (sequences.len() - 1) / 2;

    if pairs.is_empty() {
        let temp = tempfile::Builder::new().suffix(".paf").tempfile()?;
        return Ok(temp);
    }

    // If all pairs are selected (or None strategy), use fast all-vs-all
    let paf_temp = if pairs.len() == total_possible {
        sweepga_align_all_vs_all(sequences, config)?
    } else {
        log::info!(
            "sweepga: sparsification selected {} of {} pairs ({:.1}%)",
            pairs.len(),
            total_possible,
            (pairs.len() as f64 / total_possible as f64) * 100.0
        );
        sweepga_align_pairwise(sequences, &pairs, config)?
    };

    // Apply filtering unless disabled
    if config.no_filter {
        Ok(paf_temp)
    } else {
        let align_config = AlignConfig {
            num_threads: config.num_threads,
            sparsification: SparsificationStrategy::None,
            frequency_multiplier: 10,
            frequency: Some(config.kmer_frequency),
            min_alignment_length: config.min_alignment_length,
            output_format: AlignOutputFormat::Paf,
            show_progress: false,
            no_filter: config.no_filter,
            num_mappings: config.num_mappings.clone(),
            scaffold_jump: config.scaffold_jump,
            scaffold_mass: config.scaffold_mass,
            scaffold_filter: config.scaffold_filter.clone(),
            overlap: config.overlap,
            min_identity: config.min_identity,
            scaffold_dist: config.scaffold_dist,
            min_mapping_length: config.min_mapping_length,
        };
        apply_paf_filter(paf_temp, &align_config)
    }
}

/// All-vs-all alignment: write all sequences to one FASTA, align against itself.
fn sweepga_align_all_vs_all(
    sequences: &[(String, &[u8])],
    config: &SweepgaAlignConfig,
) -> io::Result<tempfile::NamedTempFile> {
    let mut combined_fasta = tempfile::Builder::new().suffix(".fa").tempfile()?;
    {
        let mut writer = BufWriter::new(&mut combined_fasta);
        for (name, seq) in sequences {
            writeln!(writer, ">{}", name)?;
            writer.write_all(seq)?;
            writeln!(writer)?;
        }
        writer.flush()?;
    }

    let fastga = FastGAIntegration::new(
        Some(config.kmer_frequency),
        config.num_threads,
        config.min_alignment_length,
        config.temp_dir.clone(),
    );

    fastga
        .align_to_temp_paf(combined_fasta.path(), combined_fasta.path())
        .map_err(|e| io::Error::other(format!("FastGA alignment failed: {}", e)))
}

/// Pairwise alignment: write individual FASTA files per pair, align each pair,
/// and combine results into a single PAF.
fn sweepga_align_pairwise(
    sequences: &[(String, &[u8])],
    pairs: &[(usize, usize)],
    config: &SweepgaAlignConfig,
) -> io::Result<tempfile::NamedTempFile> {
    let fastga = FastGAIntegration::new(
        Some(config.kmer_frequency),
        config.num_threads,
        config.min_alignment_length,
        config.temp_dir.clone(),
    );

    let mut combined_paf = tempfile::Builder::new().suffix(".paf").tempfile()?;

    for &(i, j) in pairs {
        // Write query FASTA
        let mut query_fasta = tempfile::Builder::new().suffix(".fa").tempfile()?;
        {
            let mut writer = BufWriter::new(&mut query_fasta);
            writeln!(writer, ">{}", sequences[i].0)?;
            writer.write_all(sequences[i].1)?;
            writeln!(writer)?;
            writer.flush()?;
        }

        // Write target FASTA
        let mut target_fasta = tempfile::Builder::new().suffix(".fa").tempfile()?;
        {
            let mut writer = BufWriter::new(&mut target_fasta);
            writeln!(writer, ">{}", sequences[j].0)?;
            writer.write_all(sequences[j].1)?;
            writeln!(writer)?;
            writer.flush()?;
        }

        // Align this pair
        match fastga.align_to_temp_paf(query_fasta.path(), target_fasta.path()) {
            Ok(pair_paf) => {
                // Append this pair's PAF to the combined output
                let contents = std::fs::read(pair_paf.path())?;
                if !contents.is_empty() {
                    combined_paf.write_all(&contents)?;
                }
            }
            Err(e) => {
                log::warn!(
                    "sweepga: pairwise alignment failed for {} vs {}: {}",
                    sequences[i].0,
                    sequences[j].0,
                    e
                );
            }
        }
    }

    combined_paf.flush()?;
    Ok(combined_paf)
}

/// Run alignment command
pub fn run_align(
    fasta_files: Vec<String>,
    fasta_list: Option<String>,
    output_dir: &str,
    config: AlignConfig,
) -> io::Result<()> {
    // Resolve FASTA files
    let fasta_files = resolve_fasta_files(fasta_files, fasta_list)?;

    if fasta_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No FASTA files specified",
        ));
    }

    // Create output directory
    std::fs::create_dir_all(output_dir)?;

    // Determine k-mer size for sketching
    let kmer_size = match &config.sparsification {
        SparsificationStrategy::TreeSampling { kmer_size, .. } => *kmer_size,
        _ => 15, // default
    };

    // Load sequences and compute sketches
    if config.show_progress {
        info!(
            "Loading sequences and computing sketches (k={})...",
            kmer_size
        );
    }
    let sequences = load_sequences(&fasta_files, kmer_size, 1000, config.show_progress)?;

    let n = sequences.len();
    let total_pairs = n * (n - 1) / 2;

    // Generate pairs
    if config.show_progress {
        info!(
            "Generating alignment pairs using {}...",
            config.sparsification.description()
        );
    }
    let pairs = generate_pairs(&sequences, &config.sparsification);

    if config.show_progress {
        info!(
            "Selected {} pairs out of {} total ({:.1}%)",
            pairs.len(),
            total_pairs,
            (pairs.len() as f64 / total_pairs as f64) * 100.0
        );
    }

    // Write output based on format
    match config.output_format {
        AlignOutputFormat::JobList => {
            let job_file = format!("{}/align_jobs.txt", output_dir);
            let mut writer = BufWriter::new(File::create(&job_file)?);
            write_job_list(&pairs, &sequences, output_dir, &mut writer, &config)?;
            info!("Wrote {} alignment jobs to {}", pairs.len(), job_file);
        }
        AlignOutputFormat::Paf | AlignOutputFormat::OneAln => {
            run_alignments(&pairs, &sequences, output_dir, &config)?;
        }
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
            let reader = BufReader::new(file);
            let mut files = Vec::new();

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sparsification() {
        // None
        assert!(matches!(
            SparsificationStrategy::parse("none").unwrap(),
            SparsificationStrategy::None
        ));

        // Random
        if let SparsificationStrategy::Random(f) =
            SparsificationStrategy::parse("random:0.5").unwrap()
        {
            assert!((f - 0.5).abs() < 0.001);
        } else {
            panic!("Expected Random");
        }

        // Connectivity
        if let SparsificationStrategy::Connectivity(p) =
            SparsificationStrategy::parse("giant:0.99").unwrap()
        {
            assert!((p - 0.99).abs() < 0.001);
        } else {
            panic!("Expected Connectivity");
        }

        // Tree
        if let SparsificationStrategy::TreeSampling {
            k_nearest,
            k_farthest,
            random_fraction,
            kmer_size,
        } = SparsificationStrategy::parse("tree:2:1:0.1:15").unwrap()
        {
            assert_eq!(k_nearest, 2);
            assert_eq!(k_farthest, 1);
            assert!((random_fraction - 0.1).abs() < 0.001);
            assert_eq!(kmer_size, 15);
        } else {
            panic!("Expected TreeSampling");
        }
    }

    #[test]
    fn test_connectivity_probability() {
        // Small n should give high probability
        assert!(compute_connectivity_probability(2, 0.99) >= 0.9);
        assert!(compute_connectivity_probability(5, 0.99) >= 0.5);

        // Large n should give smaller probability
        let p_100 = compute_connectivity_probability(100, 0.99);
        let p_1000 = compute_connectivity_probability(1000, 0.99);
        assert!(p_1000 < p_100);
    }

    #[test]
    fn test_sketch_and_distance() {
        let seq1 = b"ACGTACGTACGTACGT";
        let seq2 = b"ACGTACGTACGTACGT";
        let seq3 = b"GGGGGGGGGGGGGGGG";

        let sketch1 = sketch_sequence(seq1, 5, 100);
        let sketch2 = sketch_sequence(seq2, 5, 100);
        let sketch3 = sketch_sequence(seq3, 5, 100);

        // Identical sequences should have distance 0
        assert!(mash_distance(&sketch1, &sketch2, 5) < 0.01);

        // Different sequences should have higher distance
        assert!(mash_distance(&sketch1, &sketch3, 5) > 0.5);
    }
}
