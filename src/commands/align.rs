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
use sweepga::aligner::Aligner;
use sweepga::paf_filter::PafFilter;

use crate::commands::create_aligner_adaptive;

pub use sweepga::knn_graph::SparsificationStrategy;

/// Configuration for alignment command
pub struct AlignConfig {
    pub num_threads: usize,
    pub sparsify: SparsificationStrategy,
    pub mash_params: sweepga::knn_graph::MashParams,
    pub frequency_multiplier: usize,
    pub frequency: Option<usize>,
    pub min_aln_length: u64,
    pub output_format: AlignOutputFormat,
    pub show_progress: bool,
    /// Aligner backend: "wfmash" or "fastga"
    pub aligner: String,
    /// Directory for temporary files
    pub temp_dir: Option<String>,
    /// Batch alignment: max resource usage per batch (e.g., "2G", "500M")
    pub batch_bytes: Option<String>,
    // Sweepga filtering options
    pub no_filter: bool,
    pub num_mappings: String,
    pub scaffold_jump: u64,
    pub scaffold_mass: u64,
    pub scaffold_filter: String,
    pub overlap: f64,
    pub min_identity: f64,
    pub scaffold_dist: u64,
    pub min_map_length: u64,
}

impl Default for AlignConfig {
    fn default() -> Self {
        AlignConfig {
            num_threads: 4,
            sparsify: SparsificationStrategy::Connectivity(0.99),
            mash_params: sweepga::knn_graph::MashParams::default(),
            frequency_multiplier: 10,
            frequency: None,
            min_aln_length: 0,
            output_format: AlignOutputFormat::Paf,
            show_progress: true,
            aligner: "wfmash".to_string(),
            temp_dir: None,
            batch_bytes: None,
            // Filtering defaults - must match GraphBuildConfig and CLI defaults
            no_filter: false,
            num_mappings: "many:many".to_string(),
            scaffold_jump: 50_000,
            scaffold_mass: 10_000,
            scaffold_filter: "many:many".to_string(),
            overlap: 0.95,
            min_identity: 0.0,
            scaffold_dist: 0,
            min_map_length: 0,
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

/// Sequence info for distance computation (used by align command for file-based pair selection)
struct SequenceInfo {
    name: String,
    path: String,
    sketch: sweepga::mash::KmerSketch,
}

/// Group sequence indices by PanSN haplotype key (`SAMPLE#HAPLOTYPE`).
/// For non-PanSN names the whole name acts as its own key, so every contig
/// becomes its own "haplotype" and the caller falls back to contig-level
/// behavior.
fn group_indices_by_haplotype(names: &[&str]) -> Vec<Vec<usize>> {
    use std::collections::BTreeMap;
    use sweepga::pansn::{extract_pansn_key, PanSnLevel};
    let mut map: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, name) in names.iter().enumerate() {
        let key =
            extract_pansn_key(name, PanSnLevel::Haplotype).unwrap_or_else(|| (*name).to_string());
        map.entry(key).or_default().push(i);
    }
    map.into_values().collect()
}

/// Merge MinHash sketches by unioning their minimizer vectors and re-truncating
/// to `sketch_size`. MinHash bottom-k has the mergeability property that the
/// bottom-k of the union equals the bottom-k of the merged bottom-k sketches.
fn merge_sketches(
    parts: &[&sweepga::mash::KmerSketch],
    sketch_size: usize,
) -> sweepga::mash::KmerSketch {
    let k = parts
        .first()
        .map(|s| s.k)
        .unwrap_or(sweepga::mash::DEFAULT_KMER_SIZE);
    let mut minimizers: Vec<u64> = parts
        .iter()
        .flat_map(|s| s.minimizers.iter().copied())
        .collect();
    let length: usize = parts.iter().map(|s| s.length).sum();
    minimizers.sort_unstable();
    minimizers.dedup();
    minimizers.truncate(sketch_size);
    sweepga::mash::KmerSketch {
        minimizers,
        k,
        length,
    }
}

/// Expand haplotype-level pairs into contig-level pairs: cross-product across
/// selected haplotype pairs plus all intra-haplotype contig pairs.
fn expand_haplotype_pairs(
    hap_pairs: &[(usize, usize)],
    hap_groups: &[Vec<usize>],
) -> Vec<(usize, usize)> {
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    for &(hi, hj) in hap_pairs {
        for &ci in &hap_groups[hi] {
            for &cj in &hap_groups[hj] {
                if ci == cj {
                    continue;
                }
                let pair = if ci < cj { (ci, cj) } else { (cj, ci) };
                seen.insert(pair);
            }
        }
    }
    for contigs in hap_groups {
        for a in 0..contigs.len() {
            for b in a + 1..contigs.len() {
                let (ci, cj) = (contigs[a], contigs[b]);
                let pair = if ci < cj { (ci, cj) } else { (cj, ci) };
                seen.insert(pair);
            }
        }
    }
    let mut out: Vec<(usize, usize)> = seen.into_iter().collect();
    out.sort_unstable();
    out
}

/// Run `select_pairs_from_sketches` at PanSN haplotype granularity.
///
/// Groups input contigs by their `SAMPLE#HAPLOTYPE` prefix, merges per-contig
/// MinHash sketches into one sketch per haplotype, runs the sparsification
/// strategy on haplotype sketches, then expands the selected haplotype pairs
/// to contig pairs (cross-product + all intra-haplotype pairs).
///
/// Falls back to plain contig-level selection when every name is already its
/// own haplotype (non-PanSN inputs).
fn select_pairs_haplotype_aware(
    names: &[&str],
    contig_sketches: &[sweepga::mash::KmerSketch],
    strategy: &SparsificationStrategy,
    mash_params: &sweepga::knn_graph::MashParams,
) -> Vec<(usize, usize)> {
    let hap_groups = group_indices_by_haplotype(names);
    if hap_groups.len() == contig_sketches.len() {
        return sweepga::knn_graph::select_pairs_from_sketches(contig_sketches, strategy);
    }
    let hap_sketches: Vec<sweepga::mash::KmerSketch> = hap_groups
        .iter()
        .map(|idxs| {
            let parts: Vec<&sweepga::mash::KmerSketch> =
                idxs.iter().map(|&i| &contig_sketches[i]).collect();
            merge_sketches(&parts, mash_params.sketch_size)
        })
        .collect();
    let hap_pairs = sweepga::knn_graph::select_pairs_from_sketches(&hap_sketches, strategy);
    expand_haplotype_pairs(&hap_pairs, &hap_groups)
}

/// Run sketch-free pair selection (None/Random/WfmashDensity) at haplotype
/// granularity, then expand to contig pairs. `Random(f)` samples `f` of the
/// haplotype pairs, which — once expanded — gives `f` of the haplotype-pair
/// fraction, not `f` of contig pairs. That is the intended semantics: the
/// user's sparsification fraction operates on biologically meaningful units.
fn select_pairs_haplotype_aware_no_sketch(
    names: &[&str],
    strategy: &SparsificationStrategy,
    mash_params: &sweepga::knn_graph::MashParams,
) -> Vec<(usize, usize)> {
    let n = names.len();
    let hap_groups = group_indices_by_haplotype(names);
    if hap_groups.len() == n {
        return sweepga::knn_graph::select_pairs(n, None, strategy, mash_params);
    }
    let hap_pairs =
        sweepga::knn_graph::select_pairs(hap_groups.len(), None, strategy, mash_params);
    expand_haplotype_pairs(&hap_pairs, &hap_groups)
}

/// Generate alignment pairs from named in-memory sequences.
pub fn generate_pairs_for_sequences(
    sequences: &[(String, &[u8])],
    strategy: &SparsificationStrategy,
    mash_params: &sweepga::knn_graph::MashParams,
) -> Vec<(usize, usize)> {
    let n = sequences.len();
    if n <= 1 {
        return vec![];
    }
    let names: Vec<&str> = sequences.iter().map(|(name, _)| name.as_str()).collect();

    match strategy {
        SparsificationStrategy::None
        | SparsificationStrategy::Random(_)
        | SparsificationStrategy::WfmashDensity(_) => {
            select_pairs_haplotype_aware_no_sketch(&names, strategy, mash_params)
        }
        _ => {
            let raw_seqs: Vec<Vec<u8>> = sequences.iter().map(|(_, s)| s.to_vec()).collect();
            let sketches = sweepga::mash::compute_sketches_parallel(
                &raw_seqs,
                mash_params.kmer_size,
                mash_params.sketch_size,
            );
            select_pairs_haplotype_aware(&names, &sketches, strategy, mash_params)
        }
    }
}

/// Generate pairs from pre-loaded SequenceInfo (with pre-computed sketches).
fn generate_pairs(
    sequences: &[SequenceInfo],
    strategy: &SparsificationStrategy,
    mash_params: &sweepga::knn_graph::MashParams,
) -> Vec<(usize, usize)> {
    let n = sequences.len();
    if n <= 1 {
        return vec![];
    }
    let names: Vec<&str> = sequences.iter().map(|s| s.name.as_str()).collect();

    match strategy {
        SparsificationStrategy::None
        | SparsificationStrategy::Random(_)
        | SparsificationStrategy::WfmashDensity(_) => {
            select_pairs_haplotype_aware_no_sketch(&names, strategy, mash_params)
        }
        _ => {
            let sketches: Vec<sweepga::mash::KmerSketch> =
                sequences.iter().map(|s| s.sketch.clone()).collect();
            select_pairs_haplotype_aware(&names, &sketches, strategy, mash_params)
        }
    }
}

/// Load sequences from FASTA files and compute sketches using sweepga's mash module
fn load_sequences(
    fasta_files: &[String],
    show_progress: bool,
    mash_params: &sweepga::knn_graph::MashParams,
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
                        let sketch = sweepga::mash::KmerSketch::from_sequence(
                            &current_seq,
                            mash_params.kmer_size,
                            mash_params.sketch_size,
                        );
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
                let sketch = sweepga::mash::KmerSketch::from_sequence(
                    &current_seq,
                    mash_params.kmer_size,
                    mash_params.sketch_size,
                );
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

/// Replace filesystem-hostile characters (including PanSN `#`) in a string
/// so it's safe to use as a filename component.
fn sanitize_for_filename(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '/' | '\\' | '#' | ':' | ' ' | '\t' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            other => other,
        })
        .collect()
}

/// Write job list for cluster execution.
///
/// For `wfmash`: groups the selected pairs by PanSN haplotype
/// (`SAMPLE#HAPLOTYPE`) and emits one `wfmash` command per haplotype pair,
/// using `-T` (target prefix) and `-Q` (query prefix) to filter the input.
/// This matches wfmash's "one genome pair per invocation" model and lets the
/// user parallelize at haplotype granularity with `xargs -P` / `parallel`.
///
/// For `fastga`: FastGA has no prefix filter, so we emit one command per
/// unique `(target_file, query_file)` pair — not per contig pair — which
/// keeps single-file PanSN inputs from exploding into millions of identical
/// lines.
fn write_job_list<W: Write>(
    pairs: &[(usize, usize)],
    sequences: &[SequenceInfo],
    output_dir: &str,
    writer: &mut W,
    config: &AlignConfig,
) -> io::Result<usize> {
    match config.aligner.as_str() {
        "wfmash" => write_wfmash_joblist(pairs, sequences, output_dir, writer, config),
        _ => write_fastga_joblist(pairs, sequences, output_dir, writer, config),
    }
}

/// Emit one `wfmash` command per unique haplotype pair using `-T`/`-Q`.
fn write_wfmash_joblist<W: Write>(
    pairs: &[(usize, usize)],
    sequences: &[SequenceInfo],
    output_dir: &str,
    writer: &mut W,
    config: &AlignConfig,
) -> io::Result<usize> {
    use std::collections::BTreeMap;
    use sweepga::pansn::{extract_pansn_key, PanSnLevel};

    // Map each contig index to its haplotype key.
    let hap_of: Vec<String> = sequences
        .iter()
        .map(|s| {
            extract_pansn_key(&s.name, PanSnLevel::Haplotype).unwrap_or_else(|| s.name.clone())
        })
        .collect();

    // For each hap, remember one representative FASTA path (assume each hap
    // lives in a single file — the common pangenome case).
    let mut hap_to_file: BTreeMap<String, String> = BTreeMap::new();
    for (i, s) in sequences.iter().enumerate() {
        hap_to_file
            .entry(hap_of[i].clone())
            .or_insert_with(|| s.path.clone());
    }

    // Collect unique (target_hap, query_hap) pairs. Normalize ordering so
    // (H_a, H_b) and (H_b, H_a) map to the same command.
    let mut hap_pairs: std::collections::BTreeSet<(String, String)> =
        std::collections::BTreeSet::new();
    for &(i, j) in pairs {
        let (a, b) = (&hap_of[i], &hap_of[j]);
        let pair = if a <= b {
            (a.clone(), b.clone())
        } else {
            (b.clone(), a.clone())
        };
        hap_pairs.insert(pair);
    }

    let mut count = 0usize;
    for (target_hap, query_hap) in hap_pairs {
        let target_file = hap_to_file
            .get(&target_hap)
            .cloned()
            .unwrap_or_else(|| sequences[0].path.clone());
        let query_file = hap_to_file
            .get(&query_hap)
            .cloned()
            .unwrap_or_else(|| target_file.clone());

        let output_name = format!(
            "{}_vs_{}.paf",
            sanitize_for_filename(&target_hap),
            sanitize_for_filename(&query_hap),
        );
        let output_path = format!("{}/{}", output_dir, output_name);

        let mut cmd = format!("wfmash -t {}", config.num_threads);
        if config.min_aln_length > 0 {
            cmd.push_str(&format!(" -l {}", config.min_aln_length));
        }
        cmd.push_str(&format!(" -T {} -Q {}", target_hap, query_hap));
        cmd.push_str(&format!(" {}", target_file));
        // Only pass query.fa when it's a different file — otherwise let wfmash self-map.
        if query_file != target_file {
            cmd.push_str(&format!(" {}", query_file));
        }
        cmd.push_str(&format!(" > {}", output_path));

        writeln!(writer, "{}", cmd)?;
        count += 1;
    }

    Ok(count)
}

/// Emit one `FastGA` command per unique `(target_file, query_file)` pair.
fn write_fastga_joblist<W: Write>(
    pairs: &[(usize, usize)],
    sequences: &[SequenceInfo],
    output_dir: &str,
    writer: &mut W,
    config: &AlignConfig,
) -> io::Result<usize> {
    let freq_arg = if let Some(f) = config.frequency {
        format!("-f{}", f)
    } else {
        format!("-f{}", sequences.len() * config.frequency_multiplier)
    };

    let mut count = 0usize;
    for (file_i, file_j) in collect_file_pairs(pairs, sequences) {
        let stem_i = Path::new(&file_i)
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned();
        let stem_j = Path::new(&file_j)
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned();
        let output_name = format!(
            "{}_vs_{}.paf",
            sanitize_for_filename(&stem_i),
            sanitize_for_filename(&stem_j),
        );
        let output_path = format!("{}/{}", output_dir, output_name);

        writeln!(
            writer,
            "FastGA {} -T{} -l{} {} {} > {}",
            freq_arg, config.num_threads, config.min_aln_length, file_i, file_j, output_path,
        )?;
        count += 1;
    }

    Ok(count)
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
    crate::commands::count_sequences_and_genomes(fasta_files)
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

    // Check if pairs are a proper subset (sparsified) vs all-vs-all
    let total_possible = sequences.len() * (sequences.len() - 1) / 2;
    let is_sparsified = pairs.len() < total_possible;

    // Collect unique file pairs
    let file_pairs = collect_file_pairs(pairs, sequences);

    if config.show_progress {
        info!(
            "[align] {:.3}s Running {} pairwise alignments ({} unique file pairs{})",
            start_time.elapsed().as_secs_f64(),
            pairs.len(),
            file_pairs.len(),
            if is_sparsified {
                ", using --pairs-file"
            } else {
                ""
            },
        );
    }

    // Compute average sequence length for adaptive wfmash parameters
    let avg_len: Option<u64> = {
        let mut total_len: u64 = 0;
        let mut count: u64 = 0;
        for path in unique_files.iter() {
            if let Ok(fai) = rust_htslib::faidx::Reader::from_path(path) {
                let n = fai.n_seqs();
                for i in 0..n {
                    if let Ok(name) = fai.seq_name(i as i32) {
                        total_len += fai.fetch_seq_len(&name);
                        count += 1;
                    }
                }
            }
        }
        if count > 0 {
            Some(total_len / count)
        } else {
            None
        }
    };

    // Resolve wfmash mapping density from the strategy
    let wfmash_density =
        sweepga::orchestrator::resolve_wfmash_density(&config.sparsify, sequences.len());

    // When sparsified and using wfmash, write a pairs TSV so wfmash filters at L1
    // instead of running unfiltered all-vs-all per file pair.
    let pairs_file: Option<tempfile::NamedTempFile> = if is_sparsified && config.aligner == "wfmash"
    {
        let mut pairs_tsv = tempfile::Builder::new().suffix(".pairs.tsv").tempfile()?;
        {
            let mut writer = BufWriter::new(&mut pairs_tsv);
            writeln!(writer, "# query_name\ttarget_name")?;
            for &(i, j) in pairs {
                writeln!(writer, "{}\t{}", sequences[i].name, sequences[j].name)?;
                writeln!(writer, "{}\t{}", sequences[j].name, sequences[i].name)?;
            }
            writer.flush()?;
        }
        Some(pairs_tsv)
    } else {
        None
    };

    let pairs_file_path = pairs_file.as_ref().map(|f| f.path().to_path_buf());

    // Create aligner with adaptive parameters
    let aligner: Box<dyn Aligner> = create_aligner_adaptive(
        &config.aligner,
        kmer_frequency,
        config.num_threads,
        config.min_aln_length,
        Some("90".to_string()),
        config.temp_dir.clone(),
        None, // segment_length: let sweepga adapt from avg_len
        avg_len,
        wfmash_density,
        None, // num_mappings: use wfmash default
        pairs_file_path.clone(),
    )?;

    // Run alignments for each unique file pair
    let paf_output_path = format!("{}/alignments.paf", output_dir);
    let mut combined_writer = BufWriter::new(File::create(&paf_output_path)?);

    // Self-alignment (single file pair, query==target): delegate entirely to
    // sweepga::align_self_paf() which handles both batched and unbatched paths.
    let is_self_alignment = file_pairs.len() == 1 && file_pairs[0].0 == file_pairs[0].1;

    if is_self_alignment {
        let fasta_path = Path::new(&file_pairs[0].0);
        if config.show_progress {
            info!(
                "[align] {:.3}s Self-aligning {}",
                start_time.elapsed().as_secs_f64(),
                fasta_path.file_name().unwrap_or_default().to_string_lossy(),
            );
        }

        let paf_temp = match config.batch_bytes.as_deref() {
            None => sweepga::align_self_paf_direct(aligner.as_ref(), fasta_path),
            Some(batch_bytes) => sweepga::align_self_paf_batched(
                fasta_path,
                &config.aligner,
                kmer_frequency,
                config.num_threads,
                config.min_aln_length,
                Some("90".to_string()),
                config.temp_dir.clone(),
                wfmash_density,
                pairs_file_path.clone(),
                batch_bytes,
                !config.show_progress,
            ),
        }
        .map_err(|e| io::Error::other(format!("{} alignment failed: {}", config.aligner, e)))?;

        // Apply filtering unless --no-filter
        let result_paf = if config.no_filter {
            paf_temp
        } else {
            apply_paf_filter(paf_temp, config, avg_len.unwrap_or(0))?
        };

        let paf_data = std::fs::read(result_paf.path())?;
        combined_writer.write_all(&paf_data)?;
    } else {
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
                    let temp_aln = aligner
                        .align_to_temp_1aln(Path::new(query_path), Path::new(target_path))
                        .map_err(|e| {
                            io::Error::other(format!(
                                "{} alignment failed for {} vs {}: {}",
                                config.aligner, query_path, target_path, e
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
                    let paf_temp = aligner
                        .align_to_temp_paf(Path::new(query_path), Path::new(target_path))
                        .map_err(|e| {
                            io::Error::other(format!(
                                "{} alignment failed for {} vs {}: {}",
                                config.aligner, query_path, target_path, e
                            ))
                        })?;

                    // Apply filtering unless --no-filter
                    let result_paf = if config.no_filter {
                        paf_temp
                    } else {
                        apply_paf_filter(paf_temp, config, avg_len.unwrap_or(0))?
                    };

                    // Append to combined output
                    let paf_data = std::fs::read(result_paf.path())?;
                    combined_writer.write_all(&paf_data)?;
                }
                AlignOutputFormat::JobList => unreachable!(),
            }
        }

        // Keep pairs_file alive until all alignments are done
        drop(pairs_file);
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
    avg_seq_len: u64,
) -> io::Result<tempfile::NamedTempFile> {
    let filter_config = super::build_filter_config(
        &super::FilterParams {
            num_mappings: config.num_mappings.clone(),
            scaffold_jump: config.scaffold_jump,
            scaffold_mass: config.scaffold_mass,
            scaffold_filter: config.scaffold_filter.clone(),
            overlap: config.overlap,
            min_identity: config.min_identity,
            scaffold_dist: config.scaffold_dist,
            min_map_length: config.min_map_length,
        },
        avg_seq_len,
    );

    let filtered_paf_file = tempfile::Builder::new()
        .suffix(".filtered.paf")
        .tempfile()?;

    let filter = PafFilter::new(filter_config).with_keep_self(false);
    filter
        .filter_paf(paf_temp.path(), filtered_paf_file.path())
        .map_err(|e| io::Error::other(format!("PAF filtering failed: {}", e)))?;

    Ok(filtered_paf_file)
}

/// Configuration for in-memory sweepga alignment.
pub struct SweepgaAlignConfig {
    /// Number of threads for alignment
    pub num_threads: usize,
    /// K-mer frequency for FastGA
    pub kmer_frequency: usize,
    /// Minimum alignment length
    pub min_aln_length: u64,
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
    pub min_map_length: u64,
    /// Optional temp directory
    pub temp_dir: Option<String>,
    /// Unified sparsification strategy (pair selection + mapping density).
    pub sparsify: SparsificationStrategy,
    /// Mash distance parameters for sparsification sketching.
    pub mash_params: sweepga::knn_graph::MashParams,
    /// Aligner backend: "wfmash" or "fastga"
    pub aligner: String,
    /// Minimum mapping identity for wfmash (e.g. "70"). None = wfmash auto-estimates.
    pub map_pct_identity: Option<String>,
    /// Batch alignment: max resource usage per batch (e.g., "2G", "500M")
    pub batch_bytes: Option<String>,
}

impl Default for SweepgaAlignConfig {
    fn default() -> Self {
        SweepgaAlignConfig {
            num_threads: 4,
            kmer_frequency: 10,
            min_aln_length: 0,
            no_filter: false,
            num_mappings: "many:many".to_string(),
            scaffold_jump: 50_000,
            scaffold_mass: 10_000,
            scaffold_filter: "many:many".to_string(),
            overlap: 0.95,
            min_identity: 0.0,
            scaffold_dist: 0,
            min_map_length: 0,
            temp_dir: None,
            sparsify: SparsificationStrategy::None,
            mash_params: sweepga::knn_graph::MashParams::default(),
            aligner: "wfmash".to_string(),
            map_pct_identity: None,
            batch_bytes: None,
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
    // Generate pairs based on sparsification strategy
    let pairs = generate_pairs_for_sequences(sequences, &config.sparsify, &config.mash_params);
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
        let avg_seq_len = if !sequences.is_empty() {
            sequences.iter().map(|(_, s)| s.len() as u64).sum::<u64>() / sequences.len() as u64
        } else {
            0
        };
        let align_config = AlignConfig {
            num_threads: config.num_threads,
            sparsify: SparsificationStrategy::None,
            mash_params: config.mash_params.clone(),
            frequency_multiplier: 10,
            frequency: Some(config.kmer_frequency),
            min_aln_length: config.min_aln_length,
            output_format: AlignOutputFormat::Paf,
            show_progress: false,
            aligner: config.aligner.clone(),
            temp_dir: config.temp_dir.clone(),
            batch_bytes: None,
            no_filter: config.no_filter,
            num_mappings: config.num_mappings.clone(),
            scaffold_jump: config.scaffold_jump,
            scaffold_mass: config.scaffold_mass,
            scaffold_filter: config.scaffold_filter.clone(),
            overlap: config.overlap,
            min_identity: config.min_identity,
            scaffold_dist: config.scaffold_dist,
            min_map_length: config.min_map_length,
        };
        apply_paf_filter(paf_temp, &align_config, avg_seq_len)
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

    // Create FASTA index (.fai) — required by wfmash
    if config.aligner == "wfmash" {
        rust_htslib::faidx::Reader::from_path(combined_fasta.path())
            .map_err(|e| io::Error::other(format!("Failed to create FASTA index: {e}")))?;
    }

    // Compute avg sequence length for adaptive wfmash parameters
    let avg_len = if !sequences.is_empty() {
        Some(sequences.iter().map(|(_, s)| s.len() as u64).sum::<u64>() / sequences.len() as u64)
    } else {
        None
    };

    // Resolve wfmash mapping density from the strategy
    let wfmash_density =
        sweepga::orchestrator::resolve_wfmash_density(&config.sparsify, sequences.len());

    let aligner: Box<dyn Aligner> = create_aligner_adaptive(
        &config.aligner,
        config.kmer_frequency,
        config.num_threads,
        config.min_aln_length,
        config.map_pct_identity.clone(),
        config.temp_dir.clone(),
        None, // segment_length: let sweepga adapt from avg_len
        avg_len,
        wfmash_density,
        None, // num_mappings: use wfmash default (-n 1)
        None, // pairs_file
    )?;

    // Run all-vs-all alignment, with optional batching
    match config.batch_bytes.as_deref() {
        None => sweepga::align_self_paf_direct(aligner.as_ref(), combined_fasta.path()),
        Some(batch_bytes) => sweepga::align_self_paf_batched(
            combined_fasta.path(),
            &config.aligner,
            config.kmer_frequency,
            config.num_threads,
            config.min_aln_length,
            config.map_pct_identity.clone(),
            config.temp_dir.clone(),
            wfmash_density,
            None, // pairs_file: sweepga_align_all_vs_all is unsparsified
            batch_bytes,
            true, // quiet
        ),
    }
    .map_err(|e| io::Error::other(format!("{} alignment failed: {}", config.aligner, e)))
}

/// Pairwise alignment: for wfmash, write all sequences to a single FASTA and
/// use --pairs-file to restrict which pairs are aligned (single invocation).
/// For other aligners, fall back to per-pair alignment.
fn sweepga_align_pairwise(
    sequences: &[(String, &[u8])],
    pairs: &[(usize, usize)],
    config: &SweepgaAlignConfig,
) -> io::Result<tempfile::NamedTempFile> {
    // Adapt wfmash parameters to input sequence sizes.
    let avg_len = if !sequences.is_empty() {
        Some(sequences.iter().map(|(_, s)| s.len() as u64).sum::<u64>() / sequences.len() as u64)
    } else {
        None
    };

    // Resolve wfmash mapping density from the strategy
    let wfmash_density =
        sweepga::orchestrator::resolve_wfmash_density(&config.sparsify, sequences.len());

    if config.aligner == "wfmash" {
        // Optimized path: single wfmash invocation with --pairs-file
        sweepga_align_pairwise_wfmash(sequences, pairs, config, avg_len, wfmash_density)
    } else {
        // Fallback: per-pair alignment for aligners without pairs-file support
        sweepga_align_pairwise_generic(sequences, pairs, config, avg_len, wfmash_density)
    }
}

/// Optimized pairwise alignment using wfmash --pairs-file:
/// write all sequences to one FASTA, write pairs to a TSV, single wfmash call.
fn sweepga_align_pairwise_wfmash(
    sequences: &[(String, &[u8])],
    pairs: &[(usize, usize)],
    config: &SweepgaAlignConfig,
    avg_len: Option<u64>,
    wfmash_density: Option<f64>,
) -> io::Result<tempfile::NamedTempFile> {
    // 1. Write all sequences to a single combined FASTA
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

    // Create FASTA index (.fai) — required by wfmash
    rust_htslib::faidx::Reader::from_path(combined_fasta.path())
        .map_err(|e| io::Error::other(format!("Failed to create FASTA index: {e}")))?;

    // 2. Write pairs to a temp TSV file (both directions for bidirectional alignment)
    let mut pairs_tsv = tempfile::Builder::new().suffix(".pairs.tsv").tempfile()?;
    {
        let mut writer = BufWriter::new(&mut pairs_tsv);
        writeln!(writer, "# query_name\ttarget_name")?;
        for &(i, j) in pairs {
            // Write both directions so wfmash aligns A→B and B→A
            writeln!(writer, "{}\t{}", sequences[i].0, sequences[j].0)?;
            writeln!(writer, "{}\t{}", sequences[j].0, sequences[i].0)?;
        }
        writer.flush()?;
    }

    log::info!(
        "sweepga: wfmash batch alignment: {} sequences, {} pairs, pairs file: {}",
        sequences.len(),
        pairs.len(),
        pairs_tsv.path().display()
    );

    // 3. Create aligner with pairs_file set
    let aligner: Box<dyn Aligner> = create_aligner_adaptive(
        &config.aligner,
        config.kmer_frequency,
        config.num_threads,
        config.min_aln_length,
        config.map_pct_identity.clone(),
        config.temp_dir.clone(),
        None, // segment_length: let sweepga adapt from avg_len
        avg_len,
        wfmash_density,
        None, // num_mappings: use wfmash default
        Some(pairs_tsv.path().to_path_buf()),
    )?;

    // 4. Single alignment call: combined FASTA against itself, filtered by pairs
    aligner
        .align_to_temp_paf(combined_fasta.path(), combined_fasta.path())
        .map_err(|e| io::Error::other(format!("wfmash batch alignment failed: {}", e)))
}

/// Fallback per-pair alignment for aligners without --pairs-file support.
fn sweepga_align_pairwise_generic(
    sequences: &[(String, &[u8])],
    pairs: &[(usize, usize)],
    config: &SweepgaAlignConfig,
    avg_len: Option<u64>,
    wfmash_density: Option<f64>,
) -> io::Result<tempfile::NamedTempFile> {
    let aligner: Box<dyn Aligner> = create_aligner_adaptive(
        &config.aligner,
        config.kmer_frequency,
        config.num_threads,
        config.min_aln_length,
        config.map_pct_identity.clone(),
        config.temp_dir.clone(),
        None,
        avg_len,
        wfmash_density,
        None,
        None, // no pairs_file
    )?;

    let mut combined_paf = tempfile::Builder::new().suffix(".paf").tempfile()?;

    // Pre-write one FASTA file per unique sequence index to avoid O(pairs) temp file creation.
    let unique_indices: std::collections::BTreeSet<usize> =
        pairs.iter().flat_map(|&(i, j)| [i, j]).collect();
    let mut fasta_files: std::collections::HashMap<usize, tempfile::NamedTempFile> =
        std::collections::HashMap::new();
    for idx in unique_indices {
        let mut fasta = tempfile::Builder::new().suffix(".fa").tempfile()?;
        {
            let mut writer = BufWriter::new(&mut fasta);
            writeln!(writer, ">{}", sequences[idx].0)?;
            writer.write_all(sequences[idx].1)?;
            writeln!(writer)?;
            writer.flush()?;
        }
        // Create FASTA index (.fai) — required by wfmash
        if config.aligner == "wfmash" {
            rust_htslib::faidx::Reader::from_path(fasta.path())
                .map_err(|e| io::Error::other(format!("Failed to create FASTA index: {e}")))?;
        }
        fasta_files.insert(idx, fasta);
    }

    for &(i, j) in pairs {
        // Align this pair using pre-written FASTA files
        match aligner.align_to_temp_paf(fasta_files[&i].path(), fasta_files[&j].path()) {
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
    output_dir: &str,
    config: AlignConfig,
) -> io::Result<()> {
    if fasta_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No sequence files specified. Use --sequence-files or --sequence-list",
        ));
    }

    // Create output directory
    std::fs::create_dir_all(output_dir)?;

    // Load sequences and compute sketches
    if config.show_progress {
        info!("Loading sequences and computing sketches...");
    }
    let sequences = load_sequences(&fasta_files, config.show_progress, &config.mash_params)?;

    let n = sequences.len();
    let total_pairs = n * (n - 1) / 2;

    // Generate pairs
    if config.show_progress {
        info!(
            "Generating alignment pairs using {}...",
            config.sparsify.description()
        );
    }
    let pairs = generate_pairs(&sequences, &config.sparsify, &config.mash_params);

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
            let written = write_job_list(&pairs, &sequences, output_dir, &mut writer, &config)?;
            info!("Wrote {} alignment jobs to {}", written, job_file);
        }
        AlignOutputFormat::Paf | AlignOutputFormat::OneAln => {
            run_alignments(&pairs, &sequences, output_dir, &config)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sparsification() {
        // Parsing is delegated to sweepga's FromStr — verify it works through the re-export
        let s: SparsificationStrategy = "none".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::None);

        let s: SparsificationStrategy = "giant:0.99".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::Connectivity(0.99));

        let s: SparsificationStrategy = "tree:2:1:0.1".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::TreeSampling(2, 1, 0.1));

        let s: SparsificationStrategy = "wfmash:auto".parse().unwrap();
        assert_eq!(s, SparsificationStrategy::WfmashDensity(None));
    }

    #[test]
    fn test_haplotype_grouping_pansn() {
        let names = [
            "HG01#1#chr1",
            "HG01#1#chr2",
            "HG01#2#chr1",
            "HG02#1#chr1",
            "chm13#0#chrX",
        ];
        let groups = group_indices_by_haplotype(&names);
        // 4 haplotypes: HG01#1, HG01#2, HG02#1, chm13#0
        assert_eq!(groups.len(), 4);
        // HG01#1 should group indices 0 and 1
        let hg01_1 = groups.iter().find(|g| g.len() == 2).unwrap();
        assert_eq!(hg01_1, &vec![0, 1]);
    }

    #[test]
    fn test_haplotype_grouping_non_pansn() {
        // Non-PanSN names: each contig is its own haplotype
        let names = ["chr1", "chr2", "chr3"];
        let groups = group_indices_by_haplotype(&names);
        assert_eq!(groups.len(), 3);
    }

    #[test]
    fn test_expand_haplotype_pairs_cross_and_intra() {
        // Two haplotypes: hap0 = [0, 1], hap1 = [2, 3, 4]
        let groups = vec![vec![0, 1], vec![2, 3, 4]];
        // Select the single hap pair (0, 1)
        let pairs = expand_haplotype_pairs(&[(0, 1)], &groups);
        // Expect: 2*3 cross pairs + 1 intra-hap0 + 3 intra-hap1 = 10 pairs
        assert_eq!(pairs.len(), 10);
        assert!(pairs.contains(&(0, 1))); // intra hap0
        assert!(pairs.contains(&(2, 3))); // intra hap1
        assert!(pairs.contains(&(0, 2))); // cross
        assert!(pairs.contains(&(1, 4))); // cross
    }

    #[test]
    fn test_select_pairs_haplotype_aware_none_strategy() {
        // With None strategy and PanSN input, haplotype expansion must yield
        // the full set of contig pairs (equivalent to contig-level "all pairs").
        let names = ["S1#1#a", "S1#1#b", "S2#1#a", "S2#1#b"];
        let mash_params = sweepga::knn_graph::MashParams::default();
        let pairs = select_pairs_haplotype_aware_no_sketch(
            &names,
            &SparsificationStrategy::None,
            &mash_params,
        );
        // 4 choose 2 = 6 pairs total
        assert_eq!(pairs.len(), 6);
    }

    #[test]
    fn test_merge_sketches_union_semantics() {
        let a = sweepga::mash::KmerSketch::from_sequence(b"ACGTACGTACGTACGT", 5, 100);
        let b = sweepga::mash::KmerSketch::from_sequence(b"TTTTGGGGAAAACCCC", 5, 100);
        let merged = merge_sketches(&[&a, &b], 100);
        assert_eq!(merged.k, 5);
        assert_eq!(merged.length, a.length + b.length);
        // Every minimizer in merged must come from a or b
        let union: std::collections::HashSet<u64> =
            a.minimizers.iter().chain(b.minimizers.iter()).copied().collect();
        for m in &merged.minimizers {
            assert!(union.contains(m));
        }
    }

    #[test]
    fn test_mash_sketch_and_distance() {
        let seq1 = b"ACGTACGTACGTACGT";
        let seq2 = b"ACGTACGTACGTACGT";
        let seq3 = b"GGGGGGGGGGGGGGGG";

        let sketch1 = sweepga::mash::KmerSketch::from_sequence(seq1, 5, 100);
        let sketch2 = sweepga::mash::KmerSketch::from_sequence(seq2, 5, 100);
        let sketch3 = sweepga::mash::KmerSketch::from_sequence(seq3, 5, 100);

        // Identical sequences should have distance 0
        assert!(sketch1.mash_distance(&sketch2) < 0.01);

        // Different sequences should have higher distance
        assert!(sketch1.mash_distance(&sketch3) > 0.5);
    }
}
