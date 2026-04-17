use crate::graph::{prepare_poa_graph_and_sequences, prepare_sequences};
use crate::impg_index::ImpgIndex;
use crate::sequence_index::UnifiedSequenceIndex;
use coitrees::Interval;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::io;

use std::io::Write;
use std::sync::{Arc, Mutex};

use sweepga::mash::KmerSketch;

/// Hard floor on the adaptive sketch size so that Jaccard estimator variance stays bounded
/// even when the median group length is vanishingly small. At sketch size 100 the std-dev of
/// Jaccard at J=0.99 is ~1%, which is already noisy; going lower would dominate any real
/// signal in the pairwise identity column.
const FAST_SKETCH_FLOOR: usize = 100;

/// Which back-end the `similarity` command uses to produce pairwise signatures.
#[derive(Clone, Copy, Debug)]
pub enum SimilarityEngine {
    /// SPOA-based multiple sequence alignment; intersection = matching aligned columns.
    Poa {
        scoring: (u8, u8, u8, u8, u8, u8),
    },
    /// MinHash k-mer sketching; intersection = shared minimizers between group sketches.
    Mash {
        kmer_k: usize,
        sketch_cap: usize,
    },
}

/// Per-group signatures produced by the chosen engine. Exactly one variant is used per region.
enum GroupSignatures {
    /// Aligned rows (one per input sequence). Indexed by `GroupInfo::sequence_indices`.
    Msa { chars: Vec<Vec<char>> },
    /// One deduplicated, sorted minimizer set per group. Indexed by `GroupInfo::sketch_index`.
    Sketches {
        /// Sorted, deduplicated minimizers per group.
        sketches: Vec<Vec<u64>>,
    },
}

#[derive(Clone)]
struct GroupInfo {
    name: String,
    /// Indices into the per-sequence MSA rows. Used by the POA path.
    sequence_indices: Vec<usize>,
    /// POA path: sum of input sequence lengths. Mash path: size of the group's deduplicated sketch.
    total_length: usize,
    /// Mash path only: index into `GroupSignatures::Sketches::sketches`.
    sketch_index: Option<usize>,
}

/// Counts extracted from a pairwise comparison between two groups' signatures.
///
/// `matches` is reported directly in the per-pair output as `intersection` (POA) or
/// `shared.minimizers` (Mash). `mismatches` and `gap_events` are only populated in POA
/// mode. Together they drive **gap-compressed identity**:
/// `matches / (matches + mismatches + gap_events)`.
///
/// A `gap_event` is a *maximal contiguous run* of MSA columns where one group has at
/// least one non-gap base and the other has none — a single 40 bp deletion counts once,
/// not forty times. This is the standard alignment-field convention; per-base gap
/// penalties would harshly over-weight long indels, per-alignment penalties would
/// ignore them entirely. Gap-compressed sits between.
///
/// Mash mode leaves `mismatches = gap_events = 0`; its identity estimator
/// `(2J/(1+J))^(1/k)` draws its indel sensitivity from k-mer loss instead.
#[derive(Clone, Copy, Debug, Default)]
struct AlignmentCounts {
    matches: usize,
    mismatches: usize,
    gap_events: usize,
}

#[derive(Clone)]
struct SimilarityMetrics {
    jaccard: f32,
    cosine: f32,
    dice: f32,
    estimated_identity: f32,
}

impl SimilarityMetrics {
    /// Build metrics from a pair of alignment counts and two "sizes".
    ///
    /// `counts.matches` is the set intersection count used for jaccard/cosine/dice.
    ///
    /// `estimated_identity` semantics:
    /// * Mash mode (`identity_kmer_k = Some(k)`): apply the Mash-paper correction
    ///   `(2J/(1+J))^(1/k)` to recover per-base nucleotide identity from k-mer Jaccard.
    /// * POA mode (`identity_kmer_k = None`): gap-compressed alignment identity —
    ///   `matches / (matches + mismatches + gap_events)`. One indel of any length
    ///   counts as a single event.
    fn new(
        counts: AlignmentCounts,
        len_a: usize,
        len_b: usize,
        identity_kmer_k: Option<usize>,
    ) -> Self {
        let intersection = counts.matches;
        // Perfect-match shortcut: both sides the same length, every position matches,
        // no mismatches and no indel events were observed.
        let is_perfect_match = len_a == len_b
            && intersection == len_a
            && counts.mismatches == 0
            && counts.gap_events == 0;

        let union = (len_a + len_b).saturating_sub(intersection);
        let jaccard = if is_perfect_match {
            1.0
        } else if union > 0 {
            intersection as f32 / union as f32
        } else {
            0.0
        };
        let cosine = if is_perfect_match {
            1.0
        } else if len_a > 0 && len_b > 0 {
            (intersection as f32) / ((len_a as f32).sqrt() * (len_b as f32).sqrt())
        } else {
            0.0
        };
        let dice = if is_perfect_match {
            1.0
        } else if (len_a + len_b) > 0 {
            2.0 * (intersection as f32) / (len_a + len_b) as f32
        } else {
            0.0
        };
        let estimated_identity = if is_perfect_match {
            1.0
        } else if let Some(k) = identity_kmer_k {
            // Mash mode: per-base identity from k-mer Jaccard via the Mash formula.
            if jaccard > 0.0 {
                let base = 2.0 * jaccard / (1.0 + jaccard);
                if k > 1 {
                    base.powf(1.0 / k as f32)
                } else {
                    base
                }
            } else {
                0.0
            }
        } else {
            // POA mode: gap-compressed identity.
            let denom = counts.matches + counts.mismatches + counts.gap_events;
            if denom > 0 {
                counts.matches as f32 / denom as f32
            } else {
                0.0
            }
        };

        Self {
            jaccard,
            cosine,
            dice,
            estimated_identity,
        }
    }

    fn get_by_name(&self, name: &str) -> f32 {
        match name {
            "jaccard" => self.jaccard,
            "cosine" => self.cosine,
            "dice" => self.dice,
            _ => self.jaccard,
        }
    }
}

pub fn compute_and_output_similarities(
    impg: &impl ImpgIndex,
    query_data: Vec<(Vec<Interval<u32>>, String)>,
    sequence_index: &UnifiedSequenceIndex,
    engine: SimilarityEngine,
    emit_distances: bool,
    emit_all_pairs: bool,
    delim: Option<char>,
    delim_pos: u16,
    perform_pca: bool,
    n_components: usize,
    pca_similarity: &str,
    polarize_n_prev: usize,
    guide_samples: Option<&[String]>,
) -> io::Result<()> {
    // Column names for the three per-pair size/intersection columns depend on the engine:
    // POA counts aligned bases, Mash counts shared minimizers.
    let (len_a_col, len_b_col, intersection_col) = match engine {
        SimilarityEngine::Poa { .. } => ("group.a.length", "group.b.length", "intersection"),
        SimilarityEngine::Mash { .. } => (
            "group.a.sketch.size",
            "group.b.sketch.size",
            "shared.minimizers",
        ),
    };

    if !perform_pca {
        info!("Computing similarities for {} regions", query_data.len());

        // Case 1: No PCA - compute similarities in parallel and write directly

        // Write header once before parallel processing
        println!(
            "chrom\tstart\tend\tgroup.a\tgroup.b\t{len_a_col}\t{len_b_col}\t{intersection_col}\t{}",
            if emit_distances {
                "jaccard.distance\tcosine.distance\tdice.distance\testimated.difference.rate"
            } else {
                "jaccard.similarity\tcosine.similarity\tdice.similarity\testimated.identity"
            }
        );

        // Create progress bar at info level (not at error-only or debug level)
        let pb = if log::log_enabled!(log::Level::Info) && !log::log_enabled!(log::Level::Debug) {
            let progress_bar = ProgressBar::new(query_data.len() as u64);
            progress_bar.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
                    .unwrap()
                    .progress_chars("#>-")
            );
            Some(Arc::new(progress_bar))
        } else {
            None
        };

        // Create a mutex to protect stdout
        let stdout_mutex = Arc::new(Mutex::new(io::stdout()));

        // Process in parallel and write results directly
        query_data
            .par_iter()
            .try_for_each(|(query_intervals, region)| -> io::Result<()> {
                // Compute similarities for this region
                let similarity_output = compute_similarities_for_region(
                    impg,
                    query_intervals,
                    sequence_index,
                    engine,
                    emit_distances,
                    emit_all_pairs,
                    delim,
                    delim_pos,
                    region,
                )?;

                // Lock stdout and write both header and results
                let stdout = stdout_mutex.lock().unwrap();
                let mut handle = stdout.lock();
                write!(handle, "{similarity_output}")?;
                drop(handle);
                drop(stdout);

                // Update progress
                if let Some(ref progress_bar) = pb {
                    progress_bar.inc(1);
                }

                Ok(())
            })?;

        if let Some(progress_bar) = pb {
            progress_bar.finish_with_message("Completed computing similarities");
        }
    } else {
        info!("Performing PCA for {} regions", query_data.len());

        // Create progress bar at info level (not at error-only or debug level)
        let pb = if log::log_enabled!(log::Level::Info) && !log::log_enabled!(log::Level::Debug) {
            let progress_bar = ProgressBar::new(query_data.len() as u64);
            progress_bar.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta}) PCA")
                    .unwrap()
                    .progress_chars("#>-")
            );
            Some(Arc::new(progress_bar))
        } else {
            None
        };

        // Case 2 & 3: PCA with or without polarization
        let mut pca_results: Vec<_> = query_data
            .par_iter()
            .map(|(query_intervals, _)| {
                let result = compute_pca_for_region(
                    impg,
                    query_intervals,
                    sequence_index,
                    engine,
                    emit_all_pairs,
                    delim,
                    delim_pos,
                    n_components,
                    pca_similarity,
                );
                if let Some(ref progress_bar) = pb {
                    progress_bar.inc(1);
                }
                result
            })
            .collect::<Result<Vec<_>, _>>()?;

        if let Some(progress_bar) = pb {
            progress_bar.finish_with_message("Completed PCA computation");
        }

        // Apply polarization if needed
        if let Some(guide_samples) = guide_samples {
            // Guide sample polarization
            polarize_pca_result_with_guides(&mut pca_results, guide_samples)?;
        } else if polarize_n_prev > 0 {
            // Adaptive polarization
            let mut polarization_window: Vec<PolarizationData> = Vec::new();

            for pca_result in pca_results.iter_mut() {
                let polarization_data = polarize_pca_result(pca_result, &polarization_window);
                polarization_window.push(polarization_data);
                if polarization_window.len() > polarize_n_prev {
                    polarization_window.remove(0);
                }
            }
        }

        // Output results
        for (idx, pca_result) in pca_results.iter().enumerate() {
            print!("{}", pca_result.to_tsv(Some(&query_data[idx].1), idx == 0));
        }
    }

    Ok(())
}

// Helper function for computing similarities for a single region
fn compute_similarities_for_region(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    engine: SimilarityEngine,
    emit_distances: bool,
    emit_all_pairs: bool,
    delim: Option<char>,
    delim_pos: u16,
    region: &str,
) -> io::Result<String> {
    debug!("Computing similarities for region {region:?}");

    let (groups, signatures) = prepare_groups_and_signatures(
        impg,
        results,
        sequence_index,
        engine,
        delim,
        delim_pos,
    )?;

    // Identity correction: POA mode uses the k=1 special case, Mash mode uses the k-th root.
    let identity_k: Option<usize> = match engine {
        SimilarityEngine::Poa { .. } => None,
        SimilarityEngine::Mash { kmer_k, .. } => Some(kmer_k),
    };

    // Parse region once
    let (chrom, start, end) = parse_region_string(region);

    let mut output = String::new();

    // Compute only upper triangle (including diagonal) and output both (i,j) and (j,i)
    for i in 0..groups.len() {
        for j in i..groups.len() {
            let group_a = &groups[i];
            let group_b = &groups[j];

            let counts =
                calculate_intersection(&signatures, group_a, group_b, delim.is_none());

            if counts.matches == 0 && !emit_all_pairs {
                continue;
            }

            let metrics = SimilarityMetrics::new(
                counts,
                group_a.total_length,
                group_b.total_length,
                identity_k,
            );

            // Output (i,j)
            format_similarity_line(
                &mut output,
                &chrom,
                &start,
                &end,
                &group_a.name,
                &group_b.name,
                group_a.total_length,
                group_b.total_length,
                counts.matches,
                &metrics,
                emit_distances,
            );

            // Output (j,i) if different
            if i != j {
                format_similarity_line(
                    &mut output,
                    &chrom,
                    &start,
                    &end,
                    &group_b.name,
                    &group_a.name,
                    group_b.total_length,
                    group_a.total_length,
                    counts.matches,
                    &metrics,
                    emit_distances,
                );
            }
        }
    }

    Ok(output)
}

// Helper function for computing PCA for a single region
fn compute_pca_for_region(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    engine: SimilarityEngine,
    emit_all_pairs: bool,
    delim: Option<char>,
    delim_pos: u16,
    n_components: usize,
    pca_similarity: &str,
) -> io::Result<PcaResult> {
    debug!("Computing PCA for region with {} intervals", results.len());

    let (groups, signatures) = prepare_groups_and_signatures(
        impg,
        results,
        sequence_index,
        engine,
        delim,
        delim_pos,
    )?;

    let identity_k: Option<usize> = match engine {
        SimilarityEngine::Poa { .. } => None,
        SimilarityEngine::Mash { kmer_k, .. } => Some(kmer_k),
    };

    // Collect all similarities for PCA
    let mut similarities = Vec::new();

    for i in 0..groups.len() {
        for j in i..groups.len() {
            let group_a = &groups[i];
            let group_b = &groups[j];

            let counts =
                calculate_intersection(&signatures, group_a, group_b, delim.is_none());

            if counts.matches > 0 || emit_all_pairs {
                let metrics = SimilarityMetrics::new(
                    counts,
                    group_a.total_length,
                    group_b.total_length,
                    identity_k,
                );
                let similarity = metrics.get_by_name(pca_similarity);

                similarities.push((group_a.name.clone(), group_b.name.clone(), similarity));
                if i != j {
                    similarities.push((group_b.name.clone(), group_a.name.clone(), similarity));
                }
            }
        }
    }

    // Build distance matrix and perform MDS
    let (distance_matrix, labels) = build_distance_matrix(&similarities, pca_similarity)?;
    let mds = ClassicalMDS::new(n_components);

    match mds.fit_transform(&distance_matrix, labels) {
        Ok(pca_result) => Ok(pca_result),
        Err(e) => {
            if e.contains("fewer than 2 samples") || e.contains("theoretical maximum") {
                let labels: Vec<String> = groups.iter().map(|g| g.name.clone()).collect();
                Ok(PcaResult {
                    coordinates: vec![vec![0.0; n_components]; labels.len()],
                    eigenvalues: vec![0.0; n_components],
                    explained_variance_ratio: vec![0.0; n_components],
                    sample_labels: labels,
                    n_components: 0,
                })
            } else {
                Err(io::Error::other(format!("PCA/MDS failed: {e}")))
            }
        }
    }
}

// Helper to prepare groups and engine-specific signatures (MSA rows or sketches).
fn prepare_groups_and_signatures(
    impg: &impl ImpgIndex,
    results: &[Interval<u32>],
    sequence_index: &UnifiedSequenceIndex,
    engine: SimilarityEngine,
    delim: Option<char>,
    delim_pos: u16,
) -> io::Result<(Vec<GroupInfo>, GroupSignatures)> {
    match engine {
        SimilarityEngine::Poa { scoring } => {
            let (graph, sequence_metadata) =
                prepare_poa_graph_and_sequences(impg, results, sequence_index, scoring)?;

            let msa = graph.generate_msa();
            let msa_chars: Vec<Vec<char>> = msa.iter().map(|s| s.chars().collect()).collect();
            let groups = create_groups(&sequence_metadata, delim, delim_pos)?;

            Ok((groups, GroupSignatures::Msa { chars: msa_chars }))
        }
        SimilarityEngine::Mash {
            kmer_k,
            sketch_cap,
        } => {
            // 1. Fetch strand-normalized sequences + metadata via the shared helper.
            let pairs = prepare_sequences(impg, results, sequence_index)?;
            let sequences: Vec<String> = pairs.iter().map(|(s, _)| s.clone()).collect();
            let sequence_metadata: Vec<_> = pairs.into_iter().map(|(_, m)| m).collect();

            // 2. Build groups over sequence metadata (same logic as POA path).
            let mut groups = create_groups(&sequence_metadata, delim, delim_pos)?;

            // 3. Derive adaptive sketch size from the median group's k-mer capacity.
            //    Rule: `effective_sketch = min(user_cap, max(FLOOR, median_num_kmers))`.
            //    The intent is the opposite of "compress aggressively": if a group has fewer
            //    k-mers than the cap we keep them ALL (exact set Jaccard, no sampling noise);
            //    only when groups are long enough that we'd actually be throwing k-mers away
            //    does the cap kick in. Sampling too aggressively on small regions was hiding
            //    real pair-to-pair variation under Jaccard estimator noise — observed on chr6
            //    1-5 kb runs where POA disagreed with a highly-quantized sketch identity.
            let mut group_num_kmers: Vec<usize> = groups
                .iter()
                .map(|g| {
                    // Count only k-mers that can't be destroyed by the `N*k` separator.
                    // For a single sequence of length L this is max(0, L - k + 1); for
                    // concatenated group members we sum per-sequence capacities (the
                    // separator blocks any cross-boundary k-mers).
                    g.sequence_indices
                        .iter()
                        .map(|&idx| sequences[idx].len().saturating_sub(kmer_k - 1))
                        .sum()
                })
                .collect();
            group_num_kmers.sort_unstable();
            let median_num_kmers = if group_num_kmers.is_empty() {
                0
            } else {
                group_num_kmers[group_num_kmers.len() / 2]
            };
            let effective_sketch = median_num_kmers
                .max(FAST_SKETCH_FLOOR)
                .min(sketch_cap.max(FAST_SKETCH_FLOOR));
            debug!(
                "Mash sketching: {} groups, median k-mer count {}, effective sketch size {} (cap {}, k={})",
                groups.len(),
                median_num_kmers,
                effective_sketch,
                sketch_cap,
                kmer_k
            );

            // 4. Warn about sequences too short to produce any k-mers.
            for (idx, meta) in sequence_metadata.iter().enumerate() {
                if sequences[idx].len() < kmer_k {
                    debug!(
                        "Sequence '{}' is shorter than k={} ({} bp); contributes no minimizers",
                        meta.path_name(),
                        kmer_k,
                        sequences[idx].len()
                    );
                }
            }

            // 5. Sketch each group in parallel: concat sequences with an N*k separator so no
            //    k-mer can span a boundary, then build a bottom-k MinHash.
            //    Store sketches deduplicated and sorted for O(n+m) merge-intersection.
            let sketches: Vec<Vec<u64>> = groups
                .par_iter()
                .map(|g| {
                    let mut concat: Vec<u8> = Vec::new();
                    let sep = vec![b'N'; kmer_k];
                    for (i, &seq_idx) in g.sequence_indices.iter().enumerate() {
                        if i > 0 {
                            concat.extend_from_slice(&sep);
                        }
                        concat.extend_from_slice(sequences[seq_idx].as_bytes());
                    }
                    let sketch =
                        KmerSketch::from_sequence(&concat, kmer_k, effective_sketch);
                    let mut mins = sketch.minimizers;
                    // Already sorted by bottom-k; dedupe so set-Jaccard is unambiguous.
                    mins.dedup();
                    mins
                })
                .collect();

            // 6. Update GroupInfo: total_length now reports the (deduplicated) sketch size
            //    and sketch_index points at the per-group sketch.
            for (i, group) in groups.iter_mut().enumerate() {
                group.total_length = sketches[i].len();
                group.sketch_index = Some(i);
            }

            Ok((groups, GroupSignatures::Sketches { sketches }))
        }
    }
}

// Calculate intersection between two groups using the engine-specific signatures.
fn calculate_intersection(
    signatures: &GroupSignatures,
    group_a: &GroupInfo,
    group_b: &GroupInfo,
    is_individual: bool,
) -> AlignmentCounts {
    match signatures {
        GroupSignatures::Msa { chars } => {
            if is_individual
                && group_a.sequence_indices.len() == 1
                && group_b.sequence_indices.len() == 1
            {
                calculate_pairwise_intersection(
                    &chars[group_a.sequence_indices[0]],
                    &chars[group_b.sequence_indices[0]],
                )
            } else {
                calculate_group_intersection(
                    chars,
                    &group_a.sequence_indices,
                    &group_b.sequence_indices,
                )
            }
        }
        GroupSignatures::Sketches { sketches } => {
            let idx_a = group_a
                .sketch_index
                .expect("Mash signatures require sketch_index on every group");
            let idx_b = group_b
                .sketch_index
                .expect("Mash signatures require sketch_index on every group");
            // Mismatches and gap events are not recoverable from k-mer sketches; identity
            // is computed via the Mash k-th-root formula in SimilarityMetrics::new.
            AlignmentCounts {
                matches: count_shared_sorted(&sketches[idx_a], &sketches[idx_b]),
                mismatches: 0,
                gap_events: 0,
            }
        }
    }
}

/// Count the number of shared hashes between two sorted, deduplicated minimizer vectors
/// using a merge walk. O(|a| + |b|), zero allocation.
fn count_shared_sorted(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j, mut count) = (0usize, 0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Equal => {
                count += 1;
                i += 1;
                j += 1;
            }
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
        }
    }
    count
}

// Unified function for formatting similarity output
fn format_similarity_line(
    output: &mut String,
    chrom: &str,
    start: &str,
    end: &str,
    name_a: &str,
    name_b: &str,
    len_a: usize,
    len_b: usize,
    intersection: usize,
    metrics: &SimilarityMetrics,
    emit_distances: bool,
) {
    output.push_str(&format!(
        "{chrom}\t{start}\t{end}\t{name_a}\t{name_b}\t{len_a}\t{len_b}\t{intersection}\t"
    ));

    if emit_distances {
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            format_value(1.0 - metrics.jaccard),
            format_value(1.0 - metrics.cosine),
            format_value(1.0 - metrics.dice),
            format_value(1.0 - metrics.estimated_identity)
        ));
    } else {
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            format_value(metrics.jaccard),
            format_value(metrics.cosine),
            format_value(metrics.dice),
            format_value(metrics.estimated_identity)
        ));
    }
}

fn format_value(value: f32) -> String {
    let formatted = format!("{value:.7}");
    let trimmed = formatted.trim_end_matches('0');
    if trimmed.ends_with('.') {
        trimmed.trim_end_matches('.').to_string()
    } else {
        trimmed.to_string()
    }
}

fn create_groups(
    sequence_metadata: &[crate::graph::SequenceMetadata],
    delim: Option<char>,
    delim_pos: u16,
) -> io::Result<Vec<GroupInfo>> {
    if let Some(delim_char) = delim {
        // Group by delimiter
        let mut group_map: BTreeMap<String, Vec<usize>> = BTreeMap::new();

        for (seq_idx, meta) in sequence_metadata.iter().enumerate() {
            let group_name = extract_group_name(&meta.name, delim_char, delim_pos)?;
            group_map.entry(group_name).or_default().push(seq_idx);
        }

        let mut groups = Vec::new();
        for (group_name, indices) in group_map {
            // Simply sum all metadata.size values in the group
            let total_length: usize = indices
                .iter()
                .map(|&idx| sequence_metadata[idx].size as usize)
                .sum();

            groups.push(GroupInfo {
                name: group_name,
                sequence_indices: indices,
                total_length,
                sketch_index: None,
            });
        }

        Ok(groups)
    } else {
        // Each sequence is its own group
        let groups = sequence_metadata
            .iter()
            .enumerate()
            .map(|(i, meta)| {
                let name = meta.path_name();
                GroupInfo {
                    name,
                    sequence_indices: vec![i],
                    total_length: meta.size as usize,
                    sketch_index: None,
                }
            })
            .collect();

        Ok(groups)
    }
}

fn extract_group_name(path_name: &str, delim: char, delim_pos: u16) -> io::Result<String> {
    let mut delimiter_positions = Vec::new();
    let mut start = 0;

    // Find all delimiter positions
    while let Some(pos) = path_name[start..].find(delim) {
        delimiter_positions.push(start + pos);
        start += pos + 1;
    }

    let target_index = (delim_pos as usize).saturating_sub(1);

    if target_index < delimiter_positions.len() {
        Ok(path_name[..delimiter_positions[target_index]].to_string())
    } else {
        // Not enough delimiters, use the whole name
        Ok(path_name.to_string())
    }
}

// Fast path for comparing two individual sequences. Walks the two aligned rows column by
// column and classifies each column:
//   * both non-gap, equal       → matches++
//   * both non-gap, different   → mismatches++
//   * exactly one is a gap      → part of an indel run; one maximal run = one gap_event
//   * both gapped               → ignored (gap-vs-gap MSA columns carry no pair signal)
// Gap-compressed identity: indel *events* are penalized, not indel *bases*.
fn calculate_pairwise_intersection(seq_a: &[char], seq_b: &[char]) -> AlignmentCounts {
    #[derive(PartialEq, Eq, Clone, Copy)]
    enum ColState {
        Aligned,
        GapInA,
        GapInB,
        BothGap,
    }

    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut gap_events = 0usize;
    let mut prev_state: Option<ColState> = None;

    for (&char_a, &char_b) in seq_a.iter().zip(seq_b.iter()) {
        let state = match (char_a, char_b) {
            ('-', '-') => ColState::BothGap,
            ('-', _) => ColState::GapInA,
            (_, '-') => ColState::GapInB,
            (a, b) => {
                if a == b {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
                ColState::Aligned
            }
        };

        // Count gap events on entry into a GapInA / GapInB run. Consecutive gapped
        // columns on the same side collapse to a single event; a transition from
        // GapInA directly to GapInB (unusual but possible in MSAs) counts as two.
        let is_new_indel_run = matches!(state, ColState::GapInA | ColState::GapInB)
            && prev_state.as_ref() != Some(&state);
        if is_new_indel_run {
            gap_events += 1;
        }
        prev_state = Some(state);
    }

    AlignmentCounts {
        matches,
        mismatches,
        gap_events,
    }
}

// General case for comparing groups of sequences. Classifies each MSA column into one of
// four states based on which groups have at least one non-gap base:
//   * Both groups present → scale matches/mismatches by min(g_a_count, g_b_count) so the
//     per-column contribution equals the effective coverage of the smaller group.
//   * Exactly one group present → part of an indel run at the group level.
//   * Both fully gapped → ignored (no alignment signal for this pair).
// Gap events are counted at transitions *into* a group-level indel state.
fn calculate_group_intersection(
    msa_chars: &[Vec<char>],
    group_a_indices: &[usize],
    group_b_indices: &[usize],
) -> AlignmentCounts {
    #[derive(PartialEq, Eq, Clone, Copy)]
    enum ColState {
        Aligned,
        GapInA,
        GapInB,
        BothGap,
    }

    let msa_len = msa_chars[0].len();
    let mut total_matches = 0usize;
    let mut total_mismatches = 0usize;
    let mut gap_events = 0usize;
    let mut prev_state: Option<ColState> = None;

    for pos in 0..msa_len {
        let mut position_matches = 0usize;
        let mut group_a_count = 0usize;
        let mut group_b_count = 0usize;

        for &idx_a in group_a_indices {
            let char_a = msa_chars[idx_a][pos];
            if char_a == '-' {
                continue;
            }
            group_a_count += 1;
            for &idx_b in group_b_indices {
                let char_b = msa_chars[idx_b][pos];
                if char_b == '-' {
                    continue;
                }
                if char_a == char_b {
                    position_matches += 1;
                }
            }
        }
        for &idx_b in group_b_indices {
            if msa_chars[idx_b][pos] != '-' {
                group_b_count += 1;
            }
        }

        let state = match (group_a_count, group_b_count) {
            (0, 0) => ColState::BothGap,
            (0, _) => ColState::GapInA,
            (_, 0) => ColState::GapInB,
            (ga, gb) => {
                let coverage = ga.min(gb);
                let scaled_matches = position_matches.min(coverage);
                let scaled_mismatches = coverage - scaled_matches;
                total_matches += scaled_matches;
                total_mismatches += scaled_mismatches;
                ColState::Aligned
            }
        };

        let is_new_indel_run = matches!(state, ColState::GapInA | ColState::GapInB)
            && prev_state.as_ref() != Some(&state);
        if is_new_indel_run {
            gap_events += 1;
        }
        prev_state = Some(state);
    }

    AlignmentCounts {
        matches: total_matches,
        mismatches: total_mismatches,
        gap_events,
    }
}

// PCA/MDS related structures and implementations
use nalgebra::{DMatrix, SymmetricEigen};

#[derive(Debug, Clone)]
pub struct PcaResult {
    pub coordinates: Vec<Vec<f32>>,
    pub eigenvalues: Vec<f32>,
    pub explained_variance_ratio: Vec<f32>,
    pub sample_labels: Vec<String>,
    pub n_components: usize,
}

impl PcaResult {
    pub fn to_tsv(&self, region: Option<&str>, include_header: bool) -> String {
        let mut output = String::new();

        // Parse region to extract chrom, start, end
        let (chrom, start, end) = if let Some(region_str) = region {
            parse_region_string(region_str)
        } else {
            ("unknown".to_string(), "0".to_string(), "0".to_string())
        };

        // Header
        if include_header {
            output.push_str("chrom\tstart\tend\tgroup\tPC.rank\tPC.value\n");
        }

        // Data rows - one row per group per PC component
        for (group_idx, group_name) in self.sample_labels.iter().enumerate() {
            for pc_idx in 0..self.n_components {
                let pc_rank = pc_idx + 1;
                let pc_value = self.coordinates[group_idx][pc_idx];

                output.push_str(&format!(
                    "{chrom}\t{start}\t{end}\t{group_name}\t{pc_rank}\t{pc_value:.7}\n"
                ));
            }
        }

        output
    }
}

#[derive(Clone)]
pub struct PolarizationData {
    pub polarizer_indices: Vec<usize>, // One per PC component
    pub polarizer_signs: Vec<bool>,    // One per PC component
}

fn polarize_pca_result(
    pca_result: &mut PcaResult,
    polarization_window: &[PolarizationData],
) -> PolarizationData {
    let mut polarizer_indices = Vec::new();
    let mut polarizer_signs = Vec::new();

    for pc_idx in 0..pca_result.n_components {
        debug!("Processing PC{}", pc_idx + 1);

        let pc_values: Vec<f32> = pca_result
            .coordinates
            .iter()
            .map(|coords| coords.get(pc_idx).copied().unwrap_or(0.0))
            .collect();

        // Find current polarizer
        let (current_polarizer_idx, _) = pc_values
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .unwrap();

        debug!(
            "Current polarizer: [{}] {} (value: {:.3})",
            current_polarizer_idx,
            pca_result.sample_labels.get(current_polarizer_idx).unwrap(),
            pc_values[current_polarizer_idx]
        );

        if polarization_window.is_empty()
            || pc_idx >= polarization_window[0].polarizer_indices.len()
        {
            polarizer_indices.push(current_polarizer_idx);
            polarizer_signs.push(pc_values[current_polarizer_idx] > 0.0);
        } else {
            // Find most frequent polarizer from previous windows
            let mut polarizer_counts: std::collections::HashMap<usize, usize> =
                std::collections::HashMap::new();
            for prev_data in polarization_window {
                if let Some(&prev_idx) = prev_data.polarizer_indices.get(pc_idx) {
                    *polarizer_counts.entry(prev_idx).or_insert(0) += 1;
                }
            }

            let most_frequent_polarizer = polarizer_counts
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(&idx, _)| idx)
                .unwrap();

            // Determine flipping
            let mut flip_votes = 0;
            let mut total_votes = 0;

            if current_polarizer_idx == most_frequent_polarizer {
                let current_sign = pc_values[current_polarizer_idx] > 0.0;

                for prev_data in polarization_window.iter() {
                    if let Some(&prev_idx) = prev_data.polarizer_indices.get(pc_idx) {
                        if prev_idx == current_polarizer_idx {
                            if let Some(&prev_sign) = prev_data.polarizer_signs.get(pc_idx) {
                                if current_sign != prev_sign {
                                    flip_votes += 1;
                                }
                                total_votes += 1;
                            }
                        }
                    }
                }
            } else if let Some(&value) = pc_values.get(most_frequent_polarizer) {
                let current_sign_at_prev = value > 0.0;

                for prev_data in polarization_window.iter() {
                    if let Some(&prev_idx) = prev_data.polarizer_indices.get(pc_idx) {
                        if prev_idx == most_frequent_polarizer {
                            if let Some(&prev_sign) = prev_data.polarizer_signs.get(pc_idx) {
                                if current_sign_at_prev != prev_sign {
                                    flip_votes += 1;
                                }
                                total_votes += 1;
                            }
                        }
                    }
                }
            }

            let should_flip = total_votes > 0 && flip_votes > total_votes / 2;
            if should_flip {
                debug!("Flipping all PC{} values", pc_idx + 1);
                for coords in pca_result.coordinates.iter_mut() {
                    if let Some(val) = coords.get_mut(pc_idx) {
                        *val *= -1.0;
                    }
                }
                polarizer_indices.push(current_polarizer_idx);
                polarizer_signs.push(pc_values[current_polarizer_idx] <= 0.0);
            } else {
                polarizer_indices.push(current_polarizer_idx);
                polarizer_signs.push(pc_values[current_polarizer_idx] > 0.0);
            }
        }
    }

    PolarizationData {
        polarizer_indices,
        polarizer_signs,
    }
}

fn polarize_pca_result_with_guides(
    pca_results: &mut [PcaResult],
    guide_samples: &[String],
) -> io::Result<()> {
    // Find guide sample indices
    let mut guide_indices: Vec<Vec<Option<usize>>> = Vec::new();

    for guide_name in guide_samples {
        let mut indices = Vec::new();
        for result in pca_results.iter() {
            let idx = result
                .sample_labels
                .iter()
                .position(|label| label == guide_name);
            indices.push(idx);
        }
        guide_indices.push(indices);
    }

    // Check if all guide samples exist in at least one window
    for (i, guide_name) in guide_samples.iter().enumerate() {
        if guide_indices[i].iter().all(|idx| idx.is_none()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Guide sample '{guide_name}' not found in any window"),
            ));
        }
    }

    // Process each PC component
    let n_components = pca_results.first().map(|r| r.n_components).unwrap();

    for pc_idx in 0..n_components {
        // For each guide sample, track flip decisions
        let mut guide_flip_decisions: Vec<Vec<i32>> =
            vec![vec![0; pca_results.len()]; guide_samples.len()];

        for (guide_idx, guide_sample_indices) in guide_indices.iter().enumerate() {
            let mut prev_value: Option<f32> = None;

            for (window_idx, result) in pca_results.iter().enumerate() {
                if window_idx == 0 {
                    // First window: no flip
                    guide_flip_decisions[guide_idx][window_idx] = 0;

                    // Set prev_value if guide sample exists in this window
                    if let Some(sample_idx) = guide_sample_indices[window_idx] {
                        prev_value = result.coordinates[sample_idx].get(pc_idx).copied();
                    }
                } else if let Some(sample_idx) = guide_sample_indices[window_idx] {
                    if let Some(current_value) = result.coordinates[sample_idx].get(pc_idx).copied()
                    {
                        if let Some(prev) = prev_value {
                            // Check if closer to flipped or unflipped previous value
                            let dist_to_prev = (current_value - prev).abs();
                            let dist_to_flipped = (current_value - (-prev)).abs();

                            if dist_to_flipped < dist_to_prev {
                                guide_flip_decisions[guide_idx][window_idx] = 1; // Flip
                                prev_value = Some(-current_value); // Update with flipped value
                            } else {
                                guide_flip_decisions[guide_idx][window_idx] = -1; // Don't flip
                                prev_value = Some(current_value);
                            }
                        } else {
                            guide_flip_decisions[guide_idx][window_idx] = 0;
                        }
                    } else {
                        guide_flip_decisions[guide_idx][window_idx] = 0;
                    }
                } else {
                    // Guide sample not present in this window
                    guide_flip_decisions[guide_idx][window_idx] = 0;
                }
            }
        }

        // Get consensus flip decision for each window
        for window_idx in 0..pca_results.len() {
            let consensus: i32 = guide_flip_decisions
                .iter()
                .map(|decisions| decisions[window_idx])
                .sum();

            if consensus > 0 {
                // Flip this window's PC values
                for coords in pca_results[window_idx].coordinates.iter_mut() {
                    if let Some(val) = coords.get_mut(pc_idx) {
                        *val *= -1.0;
                    }
                }
            }
        }
    }

    Ok(())
}

pub struct ClassicalMDS {
    n_components: usize,
}

impl ClassicalMDS {
    pub fn new(n_components: usize) -> Self {
        Self { n_components }
    }

    pub fn fit_transform(
        &self,
        distance_matrix: &DMatrix<f32>,
        labels: Vec<String>,
    ) -> Result<PcaResult, String> {
        let n = distance_matrix.nrows();

        if n != distance_matrix.ncols() {
            return Err("Distance matrix must be square".to_string());
        }

        if labels.len() != n {
            return Err("Number of labels must match matrix dimensions".to_string());
        }

        // Theoretical maximum number of components for N points is N-1
        let max_components = n.saturating_sub(1);
        let actual_components = self.n_components.min(max_components);

        if self.n_components > max_components {
            warn!(
                "Warning: Requested {} components, but theoretical maximum for {} samples is {}. Using {} components.",
                self.n_components, n, max_components, actual_components
            );
        }

        if actual_components == 0 {
            return Err("Cannot perform MDS with fewer than 2 samples".to_string());
        }

        // Step 1: Square the distance matrix element-wise
        let d_squared = distance_matrix.map(|x| x * x);

        // Step 2: Double centering: B = -0.5 * C * D^2 * C
        let ones_matrix = DMatrix::from_element(n, n, 1.0 / n as f32);
        let identity = DMatrix::identity(n, n);
        let centering_matrix = identity - ones_matrix;
        let b_matrix = -0.5 * &centering_matrix * d_squared * &centering_matrix;

        // Step 3: Eigenvalue decomposition
        let eigen_decomp = SymmetricEigen::new(b_matrix);
        let eigenvalues = eigen_decomp.eigenvalues;
        let eigenvectors = eigen_decomp.eigenvectors;

        // Step 4: Sort eigenvalues and eigenvectors in descending order
        let mut eigen_pairs: Vec<(f32, usize)> = eigenvalues
            .iter()
            .enumerate()
            .map(|(i, &val)| (val, i))
            .collect();
        eigen_pairs.sort_by(|a, b| {
            match b.0.partial_cmp(&a.0) {
                Some(ordering) => ordering,
                None => {
                    // Handle NaN values - treat NaN as smaller than any finite value
                    if b.0.is_nan() && a.0.is_nan() {
                        std::cmp::Ordering::Equal
                    } else if b.0.is_nan() {
                        std::cmp::Ordering::Greater
                    } else {
                        std::cmp::Ordering::Less
                    }
                }
            }
        });

        // Step 5: Compute coordinates using positive eigenvalues only
        let mut coordinates = vec![vec![0.0; actual_components]; n];
        let mut result_eigenvalues = Vec::new();
        let mut component_idx = 0;

        for (eigenval, original_idx) in eigen_pairs.iter() {
            if *eigenval > 1e-10 && component_idx < actual_components {
                let sqrt_eigenval = eigenval.sqrt();
                let eigenvec = eigenvectors.column(*original_idx);

                for i in 0..n {
                    coordinates[i][component_idx] = eigenvec[i] * sqrt_eigenval;
                }

                result_eigenvalues.push(*eigenval);
                component_idx += 1;
            } else if component_idx >= actual_components {
                break; // Stop once we have enough components
            }
        }

        // If we have fewer positive eigenvalues than requested components,
        // pad with small positive values (not zeros) to maintain dimensions
        while result_eigenvalues.len() < actual_components {
            result_eigenvalues.push(1e-15); // Very small positive value instead of zero
        }

        // Calculate explained variance ratios using all positive eigenvalues
        let positive_eigenvalues: Vec<f32> = eigen_pairs
            .iter()
            .filter(|(val, _)| *val > 1e-10)
            .map(|(val, _)| *val)
            .collect();

        let total_variance: f32 = positive_eigenvalues.iter().sum();

        let explained_variance_ratio: Vec<f32> = result_eigenvalues
            .iter()
            .map(|val| {
                if total_variance > 0.0 {
                    val / total_variance
                } else {
                    0.0
                }
            })
            .collect();

        Ok(PcaResult {
            coordinates,
            eigenvalues: result_eigenvalues,
            explained_variance_ratio,
            sample_labels: labels,
            n_components: actual_components,
        })
    }
}

pub fn build_distance_matrix(
    similarities: &[(String, String, f32)],
    similarity_type: &str,
) -> Result<(DMatrix<f32>, Vec<String>), io::Error> {
    // Collect unique labels
    let mut labels = std::collections::BTreeSet::new();
    for (a, b, _) in similarities {
        labels.insert(a.clone());
        labels.insert(b.clone());
    }
    let labels: Vec<String> = labels.into_iter().collect();
    let n = labels.len();

    // Create label to index mapping
    let label_to_idx: std::collections::HashMap<String, usize> = labels
        .iter()
        .enumerate()
        .map(|(i, label)| (label.clone(), i))
        .collect();

    // Initialize distance matrix
    let mut distance_matrix = DMatrix::zeros(n, n);

    // Fill the distance matrix
    for (label_a, label_b, similarity) in similarities {
        let i = label_to_idx[label_a];
        let j = label_to_idx[label_b];

        // Convert similarity to distance based on similarity type
        let distance = match similarity_type {
            "jaccard" | "dice" | "cosine" => 1.0 - similarity,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Unknown similarity type: {similarity_type}"),
                ))
            }
        };

        distance_matrix[(i, j)] = distance;
    }

    Ok((distance_matrix, labels))
}

// Add this function to the similarity module
fn parse_region_string(region_str: &str) -> (String, String, String) {
    if let Some(colon_pos) = region_str.rfind(':') {
        let chrom = &region_str[..colon_pos];
        let range_part = &region_str[colon_pos + 1..];
        if let Some(dash_pos) = range_part.find('-') {
            let start_str = &range_part[..dash_pos];
            let end_str = &range_part[dash_pos + 1..];
            (
                chrom.to_string(),
                start_str.to_string(),
                end_str.to_string(),
            )
        } else {
            (region_str.to_string(), "0".to_string(), "0".to_string())
        }
    } else {
        (region_str.to_string(), "0".to_string(), "0".to_string())
    }
}
