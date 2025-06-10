use crate::faidx::FastaIndex;
use crate::graph::prepare_poa_graph_and_sequences;
use crate::impg::Impg;
use coitrees::Interval;
use log::warn;
use std::collections::BTreeMap;
use std::io;

#[derive(Clone)]
struct GroupInfo {
    name: String,
    sequence_indices: Vec<usize>,
    total_length: usize,
}

pub fn compute_and_output_similarities(
    impg: &Impg,
    results: &[Interval<u32>],
    fasta_index: &FastaIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
    emit_distances: bool,
    emit_all_pairs: bool,
    delim: Option<char>,
    delim_pos: u16,
    perform_pca: bool,
    n_components: usize,
    pca_similarity: &str,
    region: Option<&str>,
    include_header: bool,
) -> io::Result<()> {
    // Generate POA graph and get sequences
    let (graph, sequence_metadata) =
        prepare_poa_graph_and_sequences(impg, results, fasta_index, scoring_params)?;

    let msa = graph.generate_msa();

    // Convert MSA strings to Vec<char> once for efficient indexing
    let msa_chars: Vec<Vec<char>> = msa.iter().map(|s| s.chars().collect()).collect();

    // Create groups and calculate their lengths
    let groups = create_groups(&sequence_metadata, delim, delim_pos)?;

    if perform_pca {
        // Collect all similarities for PCA
        let mut similarities = Vec::new();

        for i in 0..groups.len() {
            for j in i..groups.len() {
                let group_a = &groups[i];
                let group_b = &groups[j];

                let intersection = if delim.is_none()
                    && group_a.sequence_indices.len() == 1
                    && group_b.sequence_indices.len() == 1
                {
                    // Fast path for individual sequences (no grouping)
                    calculate_pairwise_intersection(
                        &msa_chars[group_a.sequence_indices[0]],
                        &msa_chars[group_b.sequence_indices[0]],
                    )
                } else {
                    // General case for groups
                    calculate_group_intersection(
                        &msa_chars,
                        &group_a.sequence_indices,
                        &group_b.sequence_indices,
                    )
                };

                if intersection > 0 || emit_all_pairs {
                    // Compute the requested similarity measure
                    let similarity = match pca_similarity {
                        "jaccard" => {
                            let union = (group_a.total_length + group_b.total_length)
                                .saturating_sub(intersection);
                            if union > 0 {
                                intersection as f64 / union as f64
                            } else {
                                0.0
                            }
                        }
                        "cosine" => {
                            if group_a.total_length > 0 && group_b.total_length > 0 {
                                (intersection as f64)
                                    / ((group_a.total_length as f64).sqrt()
                                        * (group_b.total_length as f64).sqrt())
                            } else {
                                0.0
                            }
                        }
                        "dice" => {
                            if (group_a.total_length + group_b.total_length) > 0 {
                                2.0 * (intersection as f64)
                                    / (group_a.total_length + group_b.total_length) as f64
                            } else {
                                0.0
                            }
                        }
                        _ => {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidInput,
                                format!("Unknown similarity measure: {}", pca_similarity),
                            ))
                        }
                    };

                    // Add both (i,j) and (j,i) if different
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
            Ok(pca_result) => {
                print!("{}", pca_result.to_tsv(region, include_header));
            }
            Err(e) => {
                // Check if it's an insufficient samples error
                if e.contains("fewer than 2 samples") || e.contains("theoretical maximum") {
                    // Create empty result with appropriate labels
                    let labels: Vec<String> = groups.iter().map(|g| g.name.clone()).collect();
                    let empty_pca_result = PcaResult {
                        coordinates: vec![vec![0.0; n_components]; labels.len()],
                        eigenvalues: vec![0.0; n_components],
                        explained_variance_ratio: vec![0.0; n_components],
                        sample_labels: labels,
                        n_components: 0,
                    };
                    print!("{}", empty_pca_result.to_tsv(region, include_header));
                } else {
                    return Err(io::Error::other(
                        format!("PCA/MDS failed: {}", e),
                    ));
                }
            }
        }
    } else {
        // Output all similarities
        println!(
            "group.a\tgroup.b\tgroup.a.length\tgroup.b.length\tintersection\t{}",
            if emit_distances {
                "jaccard.distance\tcosine.distance\tdice.distance\testimated.difference.rate"
            } else {
                "jaccard.similarity\tcosine.similarity\tdice.similarity\testimated.identity"
            }
        );

        // Compute only upper triangle (including diagonal) and output both (i,j) and (j,i)
        for i in 0..groups.len() {
            for j in i..groups.len() {
                let group_a = &groups[i];
                let group_b = &groups[j];

                let intersection = if delim.is_none()
                    && group_a.sequence_indices.len() == 1
                    && group_b.sequence_indices.len() == 1
                {
                    // Fast path for individual sequences (no grouping)
                    calculate_pairwise_intersection(
                        &msa_chars[group_a.sequence_indices[0]],
                        &msa_chars[group_b.sequence_indices[0]],
                    )
                } else {
                    // General case for groups
                    calculate_group_intersection(
                        &msa_chars,
                        &group_a.sequence_indices,
                        &group_b.sequence_indices,
                    )
                };

                if intersection == 0 && !emit_all_pairs {
                    // Skip both (i,j) and (j,i) if no intersection
                    continue;
                }

                // Compute metrics once
                let union =
                    (group_a.total_length + group_b.total_length).saturating_sub(intersection);
                let jaccard = if union > 0 {
                    intersection as f64 / union as f64
                } else {
                    0.0
                };
                let cosine = if group_a.total_length > 0 && group_b.total_length > 0 {
                    (intersection as f64)
                        / ((group_a.total_length as f64).sqrt()
                            * (group_b.total_length as f64).sqrt())
                } else {
                    0.0
                };
                let dice = if (group_a.total_length + group_b.total_length) > 0 {
                    2.0 * (intersection as f64)
                        / (group_a.total_length + group_b.total_length) as f64
                } else {
                    0.0
                };
                let estimated_identity = if jaccard > 0.0 {
                    2.0 * jaccard / (1.0 + jaccard)
                } else {
                    0.0
                };

                // Output (i,j)
                output_similarity_line(
                    &group_a.name,
                    &group_b.name,
                    group_a.total_length,
                    group_b.total_length,
                    intersection,
                    jaccard,
                    cosine,
                    dice,
                    estimated_identity,
                    emit_distances,
                );

                // Output (j,i) if different from (i,j)
                if i != j {
                    output_similarity_line(
                        &group_b.name,
                        &group_a.name,
                        group_b.total_length,
                        group_a.total_length,
                        intersection,
                        jaccard,
                        cosine,
                        dice,
                        estimated_identity,
                        emit_distances,
                    );
                }
            }
        }
    }

    Ok(())
}

fn output_similarity_line(
    name_a: &str,
    name_b: &str,
    len_a: usize,
    len_b: usize,
    intersection: usize,
    jaccard: f64,
    cosine: f64,
    dice: f64,
    estimated_identity: f64,
    emit_distances: bool,
) {
    print!(
        "{}\t{}\t{}\t{}\t{}\t",
        name_a, name_b, len_a, len_b, intersection
    );

    if emit_distances {
        println!(
            "{}\t{}\t{}\t{}",
            format_similarity_value(1.0 - jaccard),
            format_similarity_value(1.0 - cosine),
            format_similarity_value(1.0 - dice),
            format_similarity_value(1.0 - estimated_identity)
        );
    } else {
        println!(
            "{}\t{}\t{}\t{}",
            format_similarity_value(jaccard),
            format_similarity_value(cosine),
            format_similarity_value(dice),
            format_similarity_value(estimated_identity)
        );
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
            });
        }

        Ok(groups)
    } else {
        // Each sequence is its own group
        let groups = sequence_metadata
            .iter()
            .enumerate()
            .map(|(i, meta)| {
                let name = format!("{}:{}-{}", meta.name, meta.start, meta.start + meta.size);
                GroupInfo {
                    name,
                    sequence_indices: vec![i],
                    total_length: meta.size as usize,
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

// Fast path for comparing two individual sequences
fn calculate_pairwise_intersection(seq_a: &[char], seq_b: &[char]) -> usize {
    seq_a
        .iter()
        .zip(seq_b.iter())
        .filter(|&(&char_a, &char_b)| char_a != '-' && char_b != '-' && char_a == char_b)
        .count()
}

// General case for comparing groups of sequences
fn calculate_group_intersection(
    msa_chars: &[Vec<char>],
    group_a_indices: &[usize],
    group_b_indices: &[usize],
) -> usize {
    let msa_len = msa_chars[0].len();
    let mut intersection = 0;

    for pos in 0..msa_len {
        // For each position, check all pairwise comparisons between groups and count matches
        let mut position_matches = 0;

        for &idx_a in group_a_indices {
            let char_a = msa_chars[idx_a][pos];
            if char_a != '-' {
                for &idx_b in group_b_indices {
                    let char_b = msa_chars[idx_b][pos];
                    if char_b != '-' && char_a == char_b {
                        position_matches += 1;
                    }
                }
            }
        }

        // For this position, add the minimum of:
        // - number of non-gap characters in group A
        // - number of non-gap characters in group B
        // - number of actual matches found
        let group_a_count = group_a_indices
            .iter()
            .filter(|&&idx| msa_chars[idx][pos] != '-')
            .count();
        let group_b_count = group_b_indices
            .iter()
            .filter(|&&idx| msa_chars[idx][pos] != '-')
            .count();

        intersection += position_matches.min(group_a_count).min(group_b_count);
    }

    intersection
}

fn format_similarity_value(value: f64) -> String {
    let formatted = format!("{:.7}", value);
    // Trim trailing zeros
    let trimmed = formatted.trim_end_matches('0');
    if trimmed.ends_with('.') {
        // If we removed all decimal places, remove the decimal point too
        trimmed.trim_end_matches('.').to_string()
    } else {
        trimmed.to_string()
    }
}

use nalgebra::{DMatrix, SymmetricEigen};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcaResult {
    pub coordinates: Vec<Vec<f64>>,
    pub eigenvalues: Vec<f64>,
    pub explained_variance_ratio: Vec<f64>,
    pub sample_labels: Vec<String>,
    pub n_components: usize,
}

impl PcaResult {
    pub fn to_tsv(&self, region: Option<&str>, include_header: bool) -> String {
        let mut output = String::new();

        // Parse region to extract chrom, start, end
        let (chrom, start, end) = if let Some(region_str) = region {
            if let Some(colon_pos) = region_str.rfind(':') {
                let chrom = &region_str[..colon_pos];
                let range_part = &region_str[colon_pos + 1..];
                if let Some(dash_pos) = range_part.find('-') {
                    let start_str = &range_part[..dash_pos];
                    let end_str = &range_part[dash_pos + 1..];
                    (chrom.to_string(), start_str.to_string(), end_str.to_string())
                } else {
                    (region_str.to_string(), "0".to_string(), "0".to_string())
                }
            } else {
                (region_str.to_string(), "0".to_string(), "0".to_string())
            }
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
                    "{}\t{}\t{}\t{}\t{}\t{:.7}\n",
                    chrom, start, end, group_name, pc_rank, pc_value
                ));
            }
        }

        output
    }
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
        distance_matrix: &DMatrix<f64>,
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
        let ones_matrix = DMatrix::from_element(n, n, 1.0 / n as f64);
        let identity = DMatrix::identity(n, n);
        let centering_matrix = identity - ones_matrix;
        let b_matrix = -0.5 * &centering_matrix * d_squared * &centering_matrix;

        // Step 3: Eigenvalue decomposition
        let eigen_decomp = SymmetricEigen::new(b_matrix);
        let eigenvalues = eigen_decomp.eigenvalues;
        let eigenvectors = eigen_decomp.eigenvectors;

        // Step 4: Sort eigenvalues and eigenvectors in descending order
        let mut eigen_pairs: Vec<(f64, usize)> = eigenvalues
            .iter()
            .enumerate()
            .map(|(i, &val)| (val, i))
            .collect();
        eigen_pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

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
        let positive_eigenvalues: Vec<f64> = eigen_pairs
            .iter()
            .filter(|(val, _)| *val > 1e-10)
            .map(|(val, _)| *val)
            .collect();

        let total_variance: f64 = positive_eigenvalues.iter().sum();

        let explained_variance_ratio: Vec<f64> = result_eigenvalues
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
    similarities: &[(String, String, f64)],
    similarity_type: &str,
) -> Result<(DMatrix<f64>, Vec<String>), io::Error> {
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
                    format!("Unknown similarity type: {}", similarity_type),
                ))
            }
        };

        distance_matrix[(i, j)] = distance;
        distance_matrix[(j, i)] = distance; // Ensure symmetry
    }

    Ok((distance_matrix, labels))
}
