use crate::graph::prepare_poa_graph_and_sequences;
use crate::impg::Impg;
use crate::sequence_index::UnifiedSequenceIndex;
use coitrees::Interval;
use log::{debug, info, warn};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

#[derive(Clone, Debug)]
struct Variant {
    position: usize,
    ref_allele: char,
    alt_allele: char,
    #[allow(dead_code)]
    ref_count: usize,
    #[allow(dead_code)]
    alt_count: usize,
}

#[derive(Clone)]
struct GroupInfo {
    #[allow(dead_code)]
    name: String,
    sequence_indices: Vec<usize>,
    total_sequences: usize,
}

#[derive(Clone, Debug)]
struct SelectionMetrics {
    fst: f64,
    pi_ratio: f64,
    xp_clr: f64,
    pi_group1: f64,
    pi_group2: f64,
}

#[derive(Clone)]
struct WindowMetrics {
    start: usize,
    end: usize,
    metrics: SelectionMetrics,
    variant_count: usize,
}

fn load_groups_from_file(groups_file_path: &str) -> io::Result<HashMap<String, String>> {
    let file = File::open(groups_file_path).map_err(|e| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!("Failed to open groups file '{}': {}", groups_file_path, e),
        )
    })?;
    
    let reader = BufReader::new(file);
    let mut groups_map = HashMap::new();
    let mut line_num = 0;
    
    for line_result in reader.lines() {
        line_num += 1;
        let line = line_result?;
        
        // Skip empty lines and comments
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Invalid format in groups file at line {}: expected 2 tab-separated columns (sequenceprefix, group), found {}",
                    line_num, parts.len()
                ),
            ));
        }
        
        let sequence_prefix = parts[0].trim().to_string();
        let group_name = parts[1].trim().to_string();
        
        if sequence_prefix.is_empty() || group_name.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Empty sequence prefix or group name at line {}", line_num),
            ));
        }
        
        groups_map.insert(sequence_prefix, group_name);
    }
    
    if groups_map.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Groups file is empty or contains no valid entries",
        ));
    }
    
    debug!("Loaded {} group mappings from file", groups_map.len());
    for (prefix, group) in &groups_map {
        debug!("  '{}' -> '{}'", prefix, group);
    }
    
    Ok(groups_map)
}

fn find_group_for_sequence(sequence_name: &str, groups_map: &HashMap<String, String>) -> Option<String> {
    // Try exact match first
    if let Some(group) = groups_map.get(sequence_name) {
        return Some(group.clone());
    }
    
    // Try prefix matching - find the longest matching prefix
    let mut best_match = None;
    let mut best_length = 0;
    
    for (prefix, group) in groups_map {
        if sequence_name.starts_with(prefix) && prefix.len() > best_length {
            best_match = Some(group.clone());
            best_length = prefix.len();
        }
    }
    
    best_match
}

fn create_groups(
    sequence_names: &[String],
    groups_map: &HashMap<String, String>,
) -> io::Result<Vec<GroupInfo>> {
    let mut group_map: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    let mut unassigned_sequences = Vec::new();
    
    debug!("Creating groups from {} sequences using groups mapping", sequence_names.len());
    
    for (idx, name) in sequence_names.iter().enumerate() {
        let group_name = match find_group_for_sequence(name, groups_map) {
            Some(group) => group,
            None => {
                unassigned_sequences.push((idx, name.clone()));
                continue;
            }
        };
        
        if idx < 5 {
            debug!("  Sequence {}: '{}' -> group '{}'", idx, name, group_name);
        }
        
        group_map
            .entry(group_name.clone())
            .or_insert_with(Vec::new)
            .push(idx);
    }
    
    if !unassigned_sequences.is_empty() {
        warn!("Found {} sequences without group assignment:", unassigned_sequences.len());
        for (idx, name) in &unassigned_sequences {
            warn!("  Sequence {}: '{}'", idx, name);
        }
    }
    
    debug!("Group mapping created: {} unique groups", group_map.len());
    
    let groups: Vec<GroupInfo> = group_map
        .into_iter()
        .map(|(name, indices)| {
            debug!("  Group '{}': {} sequences (indices: {:?})", 
                   name, indices.len(), 
                   if indices.len() <= 5 { &indices[..] } else { &indices[..5] });
            GroupInfo {
                total_sequences: indices.len(),
                name: name.clone(),
                sequence_indices: indices,
            }
        })
        .collect();
    
    if groups.len() < 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Need at least 2 groups for selection scan, found {}. Check your groups file and sequence naming", groups.len()),
        ));
    }
    
    Ok(groups)
}

fn detect_variants_from_msa(
    msa_chars: &[Vec<char>],
    min_allele_freq: f64,
) -> Vec<Variant> {
    if msa_chars.is_empty() {
        debug!("MSA is empty, returning no variants");
        return Vec::new();
    }
    
    let msa_len = msa_chars[0].len();
    let total_seqs = msa_chars.len();
    let min_count = (total_seqs as f64 * min_allele_freq).ceil() as usize;
    
    debug!("Detecting variants from MSA: {} sequences, {} positions, min_count={}", 
           total_seqs, msa_len, min_count);
    
    let variants: Vec<Variant> = (0..msa_len)
        .into_par_iter()
        .filter_map(|pos| {
            let mut allele_counts: HashMap<char, usize> = HashMap::new();
            let mut gap_count = 0;
            let mut n_count = 0;
            
            for seq in msa_chars {
                if pos < seq.len() {
                    let allele = seq[pos].to_ascii_uppercase();
                    if allele == '-' {
                        gap_count += 1;
                    } else if allele == 'N' {
                        n_count += 1;
                    } else {
                        *allele_counts.entry(allele).or_insert(0) += 1;
                    }
                }
            }
            
            let mut alleles: Vec<(char, usize)> = allele_counts.into_iter().collect();
            alleles.sort_by(|a, b| b.1.cmp(&a.1));
            
            // Debug logging for first few positions
            if pos < 5 {
                let allele_str: Vec<String> = alleles.iter()
                    .map(|(a, c)| format!("{}:{}", a, c))
                    .collect();
                debug!("Position {}: alleles=[{}], gaps={}, Ns={}", 
                       pos, allele_str.join(","), gap_count, n_count);
            }
            
            if alleles.len() == 2 {
                let (ref_allele, ref_count) = alleles[0];
                let (alt_allele, alt_count) = alleles[1];
                
                if ref_count >= min_count && alt_count >= min_count {
                    if pos < 10 || pos % 100 == 0 {
                        debug!("Variant at position {}: {}({}) vs {}({})",
                               pos, ref_allele, ref_count, alt_allele, alt_count);
                    }
                    return Some(Variant {
                        position: pos,
                        ref_allele,
                        alt_allele,
                        ref_count,
                        alt_count,
                    });
                }
            }
            
            None
        })
        .collect();
    
    debug!("Found {} biallelic variants", variants.len());
    variants
}

fn calculate_allele_frequencies(
    msa_chars: &[Vec<char>],
    variants: &[Variant],
    group: &GroupInfo,
) -> Vec<(f64, f64)> {
    debug!("Calculating allele frequencies for group with {} sequences, {} variants",
           group.sequence_indices.len(), variants.len());
    
    let freqs: Vec<(f64, f64)> = variants
        .iter()
        .enumerate()
        .map(|(i, var)| {
            let mut ref_count = 0;
            let mut alt_count = 0;
            let mut missing_count = 0;
            
            for &seq_idx in &group.sequence_indices {
                if seq_idx < msa_chars.len() && var.position < msa_chars[seq_idx].len() {
                    let allele = msa_chars[seq_idx][var.position].to_ascii_uppercase();
                    if allele == var.ref_allele {
                        ref_count += 1;
                    } else if allele == var.alt_allele {
                        alt_count += 1;
                    } else {
                        missing_count += 1;
                    }
                }
            }
            
            let total = (ref_count + alt_count) as f64;
            let freq = if total > 0.0 {
                (ref_count as f64 / total, alt_count as f64 / total)
            } else {
                (0.0, 0.0)
            };
            
            // Debug first few variants
            if i < 5 {
                debug!("  Variant {} at pos {}: ref={}({}/{}), alt={}({}/{}), missing={}, freq=({:.3},{:.3})",
                       i, var.position, var.ref_allele, ref_count, total as usize,
                       var.alt_allele, alt_count, total as usize, missing_count,
                       freq.0, freq.1);
            }
            
            freq
        })
        .collect();
    
    debug!("Frequency calculation complete");
    freqs
}

fn calculate_nucleotide_diversity(allele_freqs: &[(f64, f64)]) -> f64 {
    if allele_freqs.is_empty() {
        debug!("No allele frequencies provided for pi calculation");
        return 0.0;
    }
    
    let pi_per_site: f64 = allele_freqs
        .iter()
        .map(|(p, q)| 2.0 * p * q)
        .sum();
    
    let pi = pi_per_site / allele_freqs.len() as f64;
    
    debug!("Nucleotide diversity (pi): sum={:.6}, n_sites={}, pi={:.6}",
           pi_per_site, allele_freqs.len(), pi);
    
    pi
}

fn calculate_fst_hudson(
    freq1: &[(f64, f64)],
    freq2: &[(f64, f64)],
    n1: usize,
    n2: usize,
) -> f64 {
    if freq1.len() != freq2.len() || freq1.is_empty() {
        debug!("FST calculation: incompatible frequency arrays or empty");
        return 0.0;
    }
    
    debug!("Calculating FST: {} variants, n1={}, n2={}", freq1.len(), n1, n2);
    
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    let mut valid_sites = 0;
    
    for i in 0..freq1.len() {
        let p1 = freq1[i].0;
        let p2 = freq2[i].0;
        
        // Skip monomorphic sites
        if (p1 == 0.0 || p1 == 1.0) && (p2 == 0.0 || p2 == 1.0) {
            continue;
        }
        
        let p_total = (n1 as f64 * p1 + n2 as f64 * p2) / (n1 + n2) as f64;
        
        let msp = (p1 - p2).powi(2) - (p1 * (1.0 - p1)) / (n1 - 1).max(1) as f64 
                  - (p2 * (1.0 - p2)) / (n2 - 1).max(1) as f64;
        let msg = p_total * (1.0 - p_total);
        
        if i < 5 {
            debug!("  Site {}: p1={:.3}, p2={:.3}, p_total={:.3}, msp={:.6}, msg={:.6}",
                   i, p1, p2, p_total, msp, msg);
        }
        
        numerator += msp;
        denominator += msp + msg;
        valid_sites += 1;
    }
    
    let fst = if denominator > 0.0 {
        numerator / denominator
    } else {
        0.0
    };
    
    debug!("FST result: numerator={:.6}, denominator={:.6}, valid_sites={}, FST={:.6}",
           numerator, denominator, valid_sites, fst);
    
    fst
}

fn calculate_xp_clr(
    freq1: &[(f64, f64)],
    freq2: &[(f64, f64)],
    _window_size: usize,
) -> f64 {
    if freq1.len() != freq2.len() || freq1.is_empty() {
        debug!("XP-CLR: incompatible frequency arrays or empty");
        return 0.0;
    }
    
    debug!("Calculating XP-CLR for {} sites", freq1.len());
    
    let mut clr_sum = 0.0;
    let mut valid_sites = 0;
    
    for i in 0..freq1.len() {
        let p1 = freq1[i].0.max(0.01).min(0.99);
        let p2 = freq2[i].0.max(0.01).min(0.99);
        
        let lr = (p1 * (1.0 - p2)) / ((1.0 - p1) * p2);
        if lr > 0.0 {
            let log_lr = lr.ln();
            clr_sum += log_lr;
            valid_sites += 1;
            
            if i < 5 {
                debug!("  Site {}: p1={:.3}, p2={:.3}, LR={:.6}, log(LR)={:.6}",
                       i, p1, p2, lr, log_lr);
            }
        }
    }
    
    let xp_clr = if valid_sites > 0 {
        clr_sum / valid_sites as f64
    } else {
        0.0
    };
    
    debug!("XP-CLR result: sum={:.6}, valid_sites={}, XP-CLR={:.6}",
           clr_sum, valid_sites, xp_clr);
    
    xp_clr
}

fn calculate_window_metrics(
    msa_chars: &[Vec<char>],
    variants: &[Variant],
    group1: &GroupInfo,
    group2: &GroupInfo,
    window_start: usize,
    window_end: usize,
) -> SelectionMetrics {
    let window_variants: Vec<Variant> = variants
        .iter()
        .filter(|v| v.position >= window_start && v.position < window_end)
        .cloned()
        .collect();
    
    if window_variants.is_empty() {
        return SelectionMetrics {
            fst: 0.0,
            pi_ratio: 0.0,
            xp_clr: 0.0,
            pi_group1: 0.0,
            pi_group2: 0.0,
        };
    }
    
    if window_start == 0 {
        debug!("Window {}-{}: calculating metrics for {} variants",
               window_start, window_end, window_variants.len());
    }
    
    let freq1 = calculate_allele_frequencies(msa_chars, &window_variants, group1);
    let freq2 = calculate_allele_frequencies(msa_chars, &window_variants, group2);
    
    let pi1 = calculate_nucleotide_diversity(&freq1);
    let pi2 = calculate_nucleotide_diversity(&freq2);
    
    let pi_ratio = if pi2 > 0.0 { pi1 / pi2 } else { 0.0 };
    
    let fst = calculate_fst_hudson(
        &freq1,
        &freq2,
        group1.total_sequences,
        group2.total_sequences,
    );
    
    let xp_clr = calculate_xp_clr(&freq1, &freq2, window_variants.len());
    
    SelectionMetrics {
        fst,
        pi_ratio,
        xp_clr,
        pi_group1: pi1,
        pi_group2: pi2,
    }
}

pub fn compute_selection_scan(
    impg: &Impg,
    query_data: Vec<(Vec<Interval<u32>>, String)>,
    sequence_index: &UnifiedSequenceIndex,
    scoring_params: (u8, u8, u8, u8, u8, u8),
    groups_file_path: &str,
    window_size: usize,
    step_size: usize,
    min_allele_freq: f64,
    _output_format: &str,
) -> io::Result<()> {
    info!("Computing selection scans for {} regions", query_data.len());
    
    // Load group mappings from file
    let groups_map = load_groups_from_file(groups_file_path)?;
    info!("Loaded {} group mappings from {}", groups_map.len(), groups_file_path);
    
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    
    writeln!(
        handle,
        "region\twindow_start\twindow_end\tvariant_count\tFST\tpi_ratio\tXP-CLR\tpi_group1\tpi_group2"
    )?;
    
    for (intervals, name) in &query_data {
        debug!("Processing region: {}", name);
        
        let (graph, sequence_metadata) =
            prepare_poa_graph_and_sequences(impg, intervals, sequence_index, scoring_params)?;
        
        if sequence_metadata.is_empty() {
            warn!("No sequences found for region {}", name);
            continue;
        }
        
        let msa = graph.generate_msa();
        let msa_chars: Vec<Vec<char>> = msa
            .iter()
            .map(|s| s.chars().collect())
            .collect();
        
        let sequence_names: Vec<String> = sequence_metadata
            .iter()
            .map(|m| m.name.clone())
            .collect();
        
        debug!("Sequence names extracted: {:?}", 
               if sequence_names.len() <= 10 { &sequence_names[..] } else { &sequence_names[..10] });
        
        let groups = create_groups(&sequence_names, &groups_map)?;
        
        debug!("Created {} groups for region {}", groups.len(), name);
        for (i, group) in groups.iter().enumerate() {
            debug!("  Group {}: {} sequences", i, group.total_sequences);
        }
        
        if groups.len() < 2 {
            warn!("Insufficient groups for region {} (found {}), skipping", name, groups.len());
            continue;
        }
        
        debug!("MSA dimensions: {} sequences x {} positions", msa_chars.len(), 
               if !msa_chars.is_empty() { msa_chars[0].len() } else { 0 });
        
        let variants = detect_variants_from_msa(&msa_chars, min_allele_freq);
        info!("Found {} biallelic variants in region {} (min_allele_freq={:.3})", 
              variants.len(), name, min_allele_freq);
        
        if variants.is_empty() {
            warn!("No variants found in region {}, skipping", name);
            continue;
        }
        
        let msa_length = msa_chars[0].len();
        let group1 = &groups[0];
        let group2 = &groups[1];
        
        let windows: Vec<WindowMetrics> = (0..msa_length)
            .step_by(step_size)
            .par_bridge()
            .map(|start| {
                let end = (start + window_size).min(msa_length);
                let window_variants: Vec<Variant> = variants
                    .iter()
                    .filter(|v| v.position >= start && v.position < end)
                    .cloned()
                    .collect();
                
                let window_variants_count = variants
                    .iter()
                    .filter(|v| v.position >= start && v.position < end)
                    .count();
                
                if window_variants_count > 0 && start == 0 {
                    debug!("Window {}-{}: {} variants", start, end, window_variants_count);
                }
                
                let metrics = calculate_window_metrics(
                    &msa_chars,
                    &variants,
                    group1,
                    group2,
                    start,
                    end,
                );
                
                WindowMetrics {
                    start,
                    end,
                    metrics,
                    variant_count: window_variants.len(),
                }
            })
            .collect();
        
        for window in windows {
            if window.variant_count > 0 {
                writeln!(
                    handle,
                    "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                    name,
                    window.start,
                    window.end,
                    window.variant_count,
                    window.metrics.fst,
                    window.metrics.pi_ratio,
                    window.metrics.xp_clr,
                    window.metrics.pi_group1,
                    window.metrics.pi_group2,
                )?;
            }
        }
    }
    
    Ok(())
}