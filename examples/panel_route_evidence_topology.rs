//! Axis-free evidence-topology scout (assessment-side diagnostic).
//!
//! Classifies where and how a sample's occurrence-count evidence lives and
//! dies against the panel, without the axis-organized DP:
//!   1. per-homolog observation depth per catalog group (dosage signature),
//!   2. total depth trajectory along the axis,
//!   3. pair-tie spectrum: per-locus dosage hypotheses {2,0}/{1,1}/{1,0}
//!      scored with the exact production Poisson loss,
//!   4. observed-vs-panel and observed-vs-catalog feature universes,
//!   5. break taxonomy and binned boundary localization.
//!
//! The sample representation stores canonical subwalks of maximal MEM records;
//! MEMs match panel paths by construction, so a stored 3-token pair outside the
//! panel universe is structurally impossible. Novel-junction adjacency is not
//! representable in this index; junctions appear as depth steps instead. This
//! example verifies that claim empirically and reports it.
//!
//! No production behavior depends on this example. No truth information is
//! accepted on any input path.
use clap::Parser;
use impg::{
    genome_inference::{sample, PanelIdentity},
    sample_mem_bwt::{encode_walk, reverse_complement},
    syng::{SyncmerParams, SyngIndex},
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs, io,
    path::PathBuf,
};

#[derive(Parser)]
struct Options {
    #[arg(long)]
    panel: String,
    #[arg(long)]
    sample: PathBuf,
    #[arg(long)]
    axis: PathBuf,
    /// Directory holding the extracted catalog artifacts
    /// (sources.tsv, groups.tsv, occurrences.tsv, occ-features.bin,
    ///  occ-feat-offsets.bin, pairs.bin, loc-bins.bin).
    #[arg(long)]
    artifacts: PathBuf,
    #[arg(long)]
    out_dir: PathBuf,
    #[arg(long)]
    label: String,
    #[arg(long, default_value_t = 10.0)]
    depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
    /// Bin size in source-path bases for boundary localization.
    #[arg(long, default_value_t = 500)]
    bin_size: u32,
    /// How many alternative sources to evaluate per locus.
    #[arg(long, default_value_t = 8)]
    top_alts: usize,
    /// Abort when /proc VmRSS exceeds this many GiB.
    #[arg(long, default_value_t = 16.0)]
    rss_budget_gib: f64,
}

#[derive(Deserialize)]
struct Axis {
    intervals: Vec<AxisInterval>,
}
#[derive(Deserialize)]
struct AxisInterval {
    component: String,
    #[allow(dead_code)]
    start: u64,
    #[allow(dead_code)]
    end: u64,
    group: String,
    reference_occurrence: usize,
}

#[derive(Serialize, Clone)]
struct HypothesisReport {
    label: String,
    loss: f64,
}

#[derive(Serialize)]
struct LocusReport {
    index: usize,
    axis_start: u64,
    axis_end: u64,
    component: String,
    group: String,
    reference_source: String,
    ref_rho: f64,
    occurrences: usize,
    group_union_features: u64,
    group_union_mass: f64,
    unexplained_mass: f64,
    best: HypothesisReport,
    top_hypotheses: Vec<HypothesisReport>,
    top_alternatives: Vec<AltReport>,
}

#[derive(Serialize)]
struct AltReport {
    source: String,
    occurrence: usize,
    observed_mass: f64,
    feature_incidences: u64,
    rho: f64,
}

#[derive(Serialize)]
struct BinTransition {
    source: String,
    component: String,
    kind: String,
    profile: String,
    junction_interval_start: u64,
    junction_interval_end: u64,
    level_left: f64,
    level_right: f64,
    bins_left: u32,
    bins_right: u32,
}

#[derive(Serialize)]
struct BinProfile {
    source: String,
    component: String,
    bin_size: u32,
    /// [bin_start, mass, count, unique_mass, unique_count]
    bins: Vec<[u64; 5]>,
}

#[derive(Serialize)]
struct UniverseReport {
    observed_distinct_pairs: u64,
    observed_total_mass: u64,
    catalog_features: u64,
    observed_covered_by_catalog: u64,
    observed_in_panel_not_catalog: u64,
    observed_not_in_panel: u64,
    not_in_panel_examples: Vec<[u64; 3]>,
    panel_walk_paths: usize,
    structural_novelty_note: String,
}

#[derive(Serialize)]
struct Report {
    model: String,
    label: String,
    panel: String,
    sample: String,
    read_length: u64,
    histogram: u64,
    universe: UniverseReport,
    median_ref_rho: f64,
    loci: Vec<LocusReport>,
    bin_transitions: Vec<BinTransition>,
    bin_profiles: Vec<BinProfile>,
    query_accounting: serde_json::Value,
}

fn canonical3(mut pattern: [u64; 3]) -> [u64; 3] {
    let rc = reverse_complement(&pattern);
    let rc: [u64; 3] = [rc[0], rc[1], rc[2]];
    if rc < pattern {
        pattern = rc;
    }
    pattern
}

fn guard_rss(budget_gib: f64, stage: &str) -> io::Result<()> {
    let status = fs::read_to_string("/proc/self/status")?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmRSS:") {
            let kib: f64 = rest
                .trim()
                .trim_end_matches("kB")
                .trim()
                .parse()
                .map_err(io::Error::other)?;
            let gib = kib / (1024.0 * 1024.0);
            if gib > budget_gib {
                return Err(io::Error::other(format!(
                    "RSS guard exceeded at stage {stage}: {gib:.2} GiB > {budget_gib} GiB"
                )));
            }
            return Ok(());
        }
    }
    Err(io::Error::other("VmRSS missing"))
}

struct Artifacts {
    sources: Vec<String>,
    groups: Vec<(String, Vec<u32>)>,
    occurrences: Vec<OccRow>,
    offsets: Vec<u64>,
    features_file: fs::File,
    features_len: u64,
    catalog_pairs: Vec<[u64; 3]>,
    loc_rows: Vec<[u32; 4]>,
}

#[derive(Clone)]
struct OccRow {
    group: u32,
    source: u32,
    path: String,
    start: u64,
    end: u64,
}

fn load_artifacts(dir: &PathBuf) -> io::Result<Artifacts> {
    let mut sources = Vec::new();
    for line in fs::read_to_string(dir.join("sources.tsv"))?.lines() {
        let mut it = line.split('\t');
        let _id: usize = it.next().unwrap().parse().map_err(io::Error::other)?;
        let path = it.next().unwrap().to_string();
        sources.push(path);
    }
    let mut groups = Vec::new();
    for line in fs::read_to_string(dir.join("groups.tsv"))?.lines() {
        let mut it = line.split('\t');
        let _idx: usize = it.next().unwrap().parse().map_err(io::Error::other)?;
        let name = it.next().unwrap().to_string();
        let occs: Vec<u32> = match it.next() {
            Some("") | None => Vec::new(),
            Some(rest) => rest
                .split(' ')
                .map(|t| t.parse().map_err(io::Error::other))
                .collect::<io::Result<_>>()?,
        };
        groups.push((name, occs));
    }
    let mut occurrences = Vec::new();
    for line in fs::read_to_string(dir.join("occurrences.tsv"))?.lines() {
        let mut it = line.split('\t');
        let _id: u32 = it.next().unwrap().parse().map_err(io::Error::other)?;
        let group: u32 = it.next().unwrap().parse().map_err(io::Error::other)?;
        let source: u32 = it.next().unwrap().parse().map_err(io::Error::other)?;
        let path = it.next().unwrap().to_string();
        let start: u64 = it.next().unwrap().parse().map_err(io::Error::other)?;
        let end: u64 = it.next().unwrap().parse().map_err(io::Error::other)?;
        occurrences.push(OccRow {
            group,
            source,
            path,
            start,
            end,
        });
    }
    let offsets_bytes = fs::read(dir.join("occ-feat-offsets.bin"))?;
    if offsets_bytes.len() % 8 != 0 || offsets_bytes.len() / 8 != occurrences.len() {
        return Err(io::Error::other("occurrence offset table mismatch"));
    }
    let offsets: Vec<u64> = offsets_bytes
        .chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    let features_file = fs::File::open(dir.join("occ-features.bin"))?;
    let features_len = features_file.metadata()?.len();
    let pairs_bytes = fs::read(dir.join("pairs.bin"))?;
    if pairs_bytes.len() % 24 != 0 {
        return Err(io::Error::other("pairs.bin size mismatch"));
    }
    let catalog_pairs: Vec<[u64; 3]> = pairs_bytes
        .chunks_exact(24)
        .map(|c| {
            [
                u64::from_le_bytes(c[0..8].try_into().unwrap()),
                u64::from_le_bytes(c[8..16].try_into().unwrap()),
                u64::from_le_bytes(c[16..24].try_into().unwrap()),
            ]
        })
        .collect();
    let loc_bytes = fs::read(dir.join("loc-bins.bin"))?;
    let loc_rows: Vec<[u32; 4]> = loc_bytes
        .chunks_exact(16)
        .map(|c| {
            [
                u32::from_le_bytes(c[0..4].try_into().unwrap()),
                u32::from_le_bytes(c[4..8].try_into().unwrap()),
                u32::from_le_bytes(c[8..12].try_into().unwrap()),
                u32::from_le_bytes(c[12..16].try_into().unwrap()),
            ]
        })
        .collect();
    Ok(Artifacts {
        sources,
        groups,
        occurrences,
        offsets,
        features_file,
        features_len,
        catalog_pairs,
        loc_rows,
    })
}

impl Artifacts {
    /// (feature_id, multiplicity) run of one occurrence, read by offset.
    fn occurrence_features(&self, occ: usize, scratch: &mut Vec<[u32; 2]>) -> io::Result<()> {
        let start = self.offsets[occ] as usize;
        let end = if occ + 1 < self.offsets.len() {
            self.offsets[occ + 1] as usize
        } else {
            self.features_len as usize
        };
        if end < start || (end - start) % 8 != 0 {
            return Err(io::Error::other("corrupt occurrence run"));
        }
        let count = (end - start) / 8;
        scratch.clear();
        scratch.reserve(count);
        use std::io::Read;
        let mut file = &self.features_file;
        file.seek(io::SeekFrom::Start(start as u64))?;
        let mut buf = vec![0u8; end - start];
        file.read_exact(&mut buf)?;
        for c in buf.chunks_exact(8) {
            scratch.push([
                u32::from_le_bytes(c[0..4].try_into().unwrap()),
                u32::from_le_bytes(c[4..8].try_into().unwrap()),
            ]);
        }
        Ok(())
    }
}

use std::io::Seek;

/// Exact production Poisson factor (mirrors search::partition::ScoreModel::loss):
/// signal = q*histogram*depth/denominator, loss = signal - observed*ln1p(signal/background).
fn loss_factor(q: u64, observed: u64, histogram: u64, depth: f64, background: f64) -> f64 {
    let denominator = (150.0 * histogram as f64) as f64;
    let signal = (q * histogram) as f64 * depth / denominator;
    signal - observed as f64 * (signal / background).ln_1p()
}

fn main() -> io::Result<()> {
    let options = Options::parse();
    fs::create_dir_all(&options.out_dir)?;
    let rss_budget = options.rss_budget_gib;

    let identity = PanelIdentity::read(&options.panel)?;
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    guard_rss(rss_budget, "panel loaded")?;
    let sample = sample::SampleIndex::load(&options.sample, &identity)?;
    guard_rss(rss_budget, "sample loaded")?;
    let histogram = *sample
        .stats
        .read_lengths
        .get(&150)
        .ok_or_else(|| io::Error::other("sample lacks L150 histogram"))?;

    let observed: BTreeMap<[u64; 3], u64> = sample.counts.observed_pairs()?;
    guard_rss(rss_budget, "observed pairs enumerated")?;
    let observed_total_mass: u64 = observed.values().sum();
    eprintln!(
        "observed universe: {} distinct pairs, total mass {}",
        observed.len(),
        observed_total_mass
    );

    let art = load_artifacts(&options.artifacts)?;
    guard_rss(rss_budget, "artifacts loaded")?;
    eprintln!(
        "catalog: {} features, {} occurrences, {} groups, {} sources, {} loc rows",
        art.catalog_pairs.len(),
        art.occurrences.len(),
        art.groups.len(),
        art.sources.len(),
        art.loc_rows.len()
    );

    // ---- universe join ------------------------------------------------
    // Canonical catalog pair multiset for membership.
    let mut catalog_canonical: Vec<[u64; 3]> = art
        .catalog_pairs
        .iter()
        .map(|&p| canonical3(p))
        .collect();
    catalog_canonical.sort_unstable();
    catalog_canonical.dedup();
    let mut covered_by_catalog = 0u64;
    let mut in_panel_not_catalog: Vec<[u64; 3]> = Vec::new();
    let mut observed_keys: HashSet<[u64; 3]> = HashSet::new();
    for (&pair, _) in observed.iter() {
        observed_keys.insert(pair);
        if catalog_canonical.binary_search(&pair).is_ok() {
            covered_by_catalog += 1;
        } else {
            in_panel_not_catalog.push(pair);
        }
    }
    // Panel-walk check: every observed pair must occur in some panel path.
    // MEM subwalks match panel paths by construction, so misses should be
    // structurally impossible; verify empirically on the full observed set.
    let mut panel_hits: HashSet<[u64; 3]> = HashSet::new();
    let mut walked_paths = 0usize;
    'paths: for (_path, start) in panel.name_map.path_starts.iter().enumerate() {
        let Some(start) = start else { continue };
        let tape = encode_walk(&panel.walk_forward_path(start))?;
        walked_paths += 1;
        let mut i = 0usize;
        while i + 2 < tape.len() {
            let triple = [tape[i], tape[i + 1], tape[i + 2]];
            if observed_keys.contains(&triple) || observed_keys.contains(&canonical3(triple)) {
                panel_hits.insert(canonical3(triple));
                if panel_hits.len() == observed_keys.len() {
                    break 'paths;
                }
            }
            i += 2;
        }
    }
    guard_rss(rss_budget, "panel walk")?;
    let not_in_panel: Vec<[u64; 3]> = observed
        .keys()
        .filter(|p| !panel_hits.contains(*p))
        .copied()
        .collect();
    eprintln!(
        "universe join: covered_by_catalog={}; in_panel_not_catalog={}; not_in_panel={}",
        covered_by_catalog,
        in_panel_not_catalog.len(),
        not_in_panel.len()
    );

    // ---- C_f cache over catalog features -------------------------------
    let n_features = art.catalog_pairs.len();
    let mut c_cache: Vec<u64> = vec![u64::MAX; n_features];
    let count_of = |fid: u32,
                    c_cache: &mut Vec<u64>,
                    catalog_pairs: &[[u64; 3]]|
     -> Option<u64> {
        let fid = fid as usize;
        if fid >= c_cache.len() {
            return None;
        }
        if c_cache[fid] == u64::MAX {
            let key = canonical3(catalog_pairs[fid]);
            let c = observed.get(&key).copied().unwrap_or(0);
            c_cache[fid] = c;
        }
        Some(c_cache[fid])
    };

    // Feature panel positional multiplicity: one sequential pass over the
    // occurrence feature runs. Repeat-family features pool observations
    // across the panel; unique-feature profiles carry the clean boundary signal.
    let mut panel_multiplicity: Vec<u32> = vec![0; n_features];
    {
        use std::io::{Read, Seek, SeekFrom};
        let mut file = &art.features_file;
        file.seek(SeekFrom::Start(0))?;
        let mut remaining = art.features_len as usize;
        let mut buf = vec![0u8; 1 << 23];
        while remaining > 0 {
            let want = buf.len().min(remaining);
            file.read_exact(&mut buf[..want])?;
            for c in buf[..want].chunks_exact(8) {
                let fid = u32::from_le_bytes(c[0..4].try_into().unwrap()) as usize;
                let mult = u32::from_le_bytes(c[4..8].try_into().unwrap());
                panel_multiplicity[fid] = panel_multiplicity[fid].saturating_add(mult);
            }
            remaining -= want;
        }
    }
    let unique_threshold: u32 = 4;
    eprintln!("panel multiplicity pass complete");

    let axis: Axis = serde_json::from_reader(fs::File::open(&options.axis)?)
        .map_err(io::Error::other)?;
    let group_index: HashMap<&str, usize> = art
        .groups
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.as_str(), i))
        .collect();

    // ---- per-locus scan -------------------------------------------------
    let mut loci: Vec<LocusReport> = Vec::new();
    let mut scratch: Vec<[u32; 2]> = Vec::new();
    let mut feature_gen: Vec<u32> = vec![0; n_features];
    let mut gen: u32 = 0;
    let mut median_ref_rho_basis: Vec<(usize, f64)> = Vec::new();

    for (index, interval) in axis.intervals.iter().enumerate() {
        gen += 1;
        let gidx = *group_index
            .get(interval.group.as_str())
            .ok_or_else(|| io::Error::other(format!("axis group missing: {}", interval.group)))?;
        let ref_occ = interval.reference_occurrence;
        if ref_occ >= art.occurrences.len() {
            return Err(io::Error::other("axis reference occurrence out of range"));
        }
        let group_occs = &art.groups[gidx].1;
        // Per-occurrence masses.
        let mut occ_mass: Vec<(usize, f64, u64)> = Vec::with_capacity(group_occs.len());
        let mut union_mass = 0f64;
        let mut union_features = 0u64;
        for &occ in group_occs {
            let occ = occ as usize;
            art.occurrence_features(occ, &mut scratch)?;
            let mut mass = 0f64;
            let mut incidences = 0u64;
            for &[fid, mult] in scratch.iter() {
                let Some(c) = count_of(fid, &mut c_cache, &art.catalog_pairs) else {
                    return Err(io::Error::other("feature id out of range"));
                };
                mass += c as f64 * mult as f64;
                incidences += mult as u64;
                if feature_gen[fid as usize] != gen {
                    feature_gen[fid as usize] = gen;
                    union_features += 1;
                    union_mass += c as f64;
                }
            }
            occ_mass.push((occ, mass, incidences));
        }
        let (ref_mass, ref_incidences) = occ_mass
            .iter()
            .find(|(occ, _, _)| *occ == ref_occ)
            .map(|(_, m, f)| (*m, *f))
            .unwrap_or((0.0, 0));
        let ref_rho = if ref_incidences > 0 {
            ref_mass / ref_incidences as f64
        } else {
            0.0
        };
        median_ref_rho_basis.push((index, ref_rho));

        // Alternative sources by observed mass, distinct source paths.
        let mut by_source: HashMap<u32, (usize, f64, u64)> = HashMap::new();
        for (occ, mass, incidences) in occ_mass.iter() {
            if *occ == ref_occ {
                continue;
            }
            let src = art.occurrences[*occ].source;
            let slot = by_source.entry(src).or_insert((*occ, 0.0, 0));
            slot.1 += *mass;
            slot.2 += *incidences;
        }
        let mut alts: Vec<(u32, usize, f64, u64)> = by_source
            .into_iter()
            .map(|(src, (occ, mass, inc))| (src, occ, mass, inc))
            .collect();
        alts.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
        alts.truncate(options.top_alts);

        // Dosage hypotheses over each candidate pair's own features
        // (production candidate-conditional convention).
        let mut hypotheses: Vec<HypothesisReport> = Vec::new();
        let mut evaluate = |_label: String,
                            copies: &[(usize, u64)]|
         -> io::Result<f64> {
            // union features with per-occurrence multiplicities
            let mut union: HashMap<u32, (u64, u64)> = HashMap::new();
            for (occ, _mult) in copies.iter() {
                art.occurrence_features(*occ, &mut scratch)?;
                for &[fid, m] in scratch.iter() {
                    let slot = union.entry(fid).or_insert((0, 0));
                    if *occ == copies[0].0 {
                        slot.0 += m as u64 * copies[0].1;
                    } else {
                        slot.1 += m as u64 * copies[1].1;
                    }
                }
            }
            let mut total = 0f64;
            for (fid, (q0, q1)) in union.iter() {
                let Some(c) = count_of(*fid, &mut c_cache, &art.catalog_pairs) else {
                    return Err(io::Error::other("feature id out of range"));
                };
                let q = q0 + q1;
                total += loss_factor(q, c, histogram, options.depth, options.background);
            }
            Ok(total)
        };
        hypotheses.push(HypothesisReport {
            label: format!("{}x2", art.occurrences[ref_occ].path),
            loss: evaluate(
                format!("{}x2", art.occurrences[ref_occ].path),
                &[(ref_occ, 2)],
            )?,
        });
        hypotheses.push(HypothesisReport {
            label: format!("{}x1", art.occurrences[ref_occ].path),
            loss: evaluate(
                format!("{}x1", art.occurrences[ref_occ].path),
                &[(ref_occ, 1)],
            )?,
        });
        let mut pair_losses: Vec<(usize, HypothesisReport)> = Vec::new();
        let mut homo_alt: Option<HypothesisReport> = None;
        for (rank, (src, occ, _mass, _inc)) in alts.iter().enumerate() {
            let label = format!(
                "{}+{}",
                art.occurrences[ref_occ].path, art.sources[*src as usize]
            );
            let loss = evaluate(label.clone(), &[(ref_occ, 1), (*occ, 1)])?;
            pair_losses.push((
                rank,
                HypothesisReport {
                    label,
                    loss,
                },
            ));
            if rank == 0 {
                let label2 = format!("{}x2", art.sources[*src as usize]);
                let loss2 = evaluate(label2.clone(), &[(*occ, 2)])?;
                homo_alt = Some(HypothesisReport {
                    label: label2,
                    loss: loss2,
                });
            }
        }
        if let Some(h) = homo_alt {
            hypotheses.push(h);
        }
        hypotheses.extend(pair_losses.into_iter().map(|(_, h)| h));
        hypotheses.sort_by(|a, b| a.loss.partial_cmp(&b.loss).unwrap());
        let best = hypotheses.first().cloned().unwrap();
        let top_hypotheses: Vec<HypothesisReport> =
            hypotheses.iter().take(4).cloned().collect();

        // Unexplained mass: group union mass minus the best pair's explained
        // mass is not directly subtractable across different feature sets;
        // report the union mass over features of the two best candidates.
        let mut explained = 0f64;
        let mut explained_union: HashSet<u32> = HashSet::new();
        for occ in [ref_occ, alts.first().map(|a| a.1).unwrap_or(ref_occ)] {
            art.occurrence_features(occ, &mut scratch)?;
            for &[fid, _m] in scratch.iter() {
                if explained_union.insert(fid) {
                    if let Some(c) = count_of(fid, &mut c_cache, &art.catalog_pairs) {
                        explained += c as f64;
                    }
                }
            }
        }
        loci.push(LocusReport {
            index,
            axis_start: interval.start,
            axis_end: interval.end,
            component: interval.component.clone(),
            group: interval.group.clone(),
            reference_source: art.occurrences[ref_occ].path.clone(),
            ref_rho: ref_rho,
            occurrences: group_occs.len(),
            group_union_features: union_features,
            group_union_mass: union_mass,
            unexplained_mass: (union_mass - explained).max(0.0),
            best,
            top_hypotheses,
            top_alternatives: alts
                .iter()
                .map(|(src, occ, mass, inc)| AltReport {
                    source: art.sources[*src as usize].clone(),
                    occurrence: *occ,
                    observed_mass: *mass,
                    feature_incidences: *inc,
                    rho: if *inc > 0 { *mass / *inc as f64 } else { 0.0 },
                })
                .collect(),
        });
        if index % 200 == 0 {
            eprintln!("locus {index}/{} done", axis.intervals.len());
        }
    }
    let mut rhos: Vec<f64> = median_ref_rho_basis.iter().map(|(_, r)| *r).collect();
    rhos.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_ref_rho = rhos[rhos.len() / 2];
    eprintln!("median ref rho = {median_ref_rho}");

    // ---- binned boundary localization ------------------------------------
    // loc_rows: (feature_id, source_id, start, end) for S288C#0#*/SK1#0#*.
    let mut bin_transitions: Vec<BinTransition> = Vec::new();
    let mut bin_profiles: Vec<BinProfile> = Vec::new();
    {
        // group rows by (source, chromosome-from-path)
        let mut by_src: BTreeMap<u32, Vec<[u32; 4]>> = BTreeMap::new();
        for row in art.loc_rows.iter() {
            by_src.entry(row[1]).or_default().push(*row);
        }
        for (src, rows) in by_src {
            let path = art
                .sources
                .get(src as usize)
                .cloned()
                .unwrap_or_else(|| format!("src{src}"));
            let component = path.clone();
            let max_start = rows.iter().map(|r| r[2]).max().unwrap_or(0) as usize;
            let nbins = max_start / options.bin_size as usize + 2;
            let mut mass = vec![0f64; nbins];
            let mut count = vec![0u64; nbins];
            let mut umass = vec![0f64; nbins];
            let mut ucount = vec![0u64; nbins];
            for row in rows {
                let Some(c) = count_of(row[0], &mut c_cache, &art.catalog_pairs) else {
                    continue;
                };
                let b = row[2] as usize / options.bin_size as usize;
                mass[b] += c as f64;
                count[b] += 1;
                if (panel_multiplicity[row[0] as usize] as u64)
                    <= unique_threshold as u64
                {
                    umass[b] += c as f64;
                    ucount[b] += 1;
                }
            }
            let w = 3usize; // window in bins on each side
            let min_count = 20u64;
            let mut profile_bins: Vec<[u64; 5]> = (0..nbins)
                .map(|b| {
                    [
                        b as u64 * options.bin_size as u64,
                        mass[b] as u64,
                        count[b],
                        umass[b] as u64,
                        ucount[b],
                    ]
                })
                .collect();
            bin_profiles.push(BinProfile {
                source: path.clone(),
                component: component.clone(),
                bin_size: options.bin_size,
                bins: profile_bins,
            });
            let mut detect = |profile: &str, mass: &[f64], count: &[u64], min_count: u64, w: usize| {
                for b in w..nbins.saturating_sub(w) {
                    let (mut lm, mut lc) = (0f64, 0u64);
                    let (mut rm, mut rc) = (0f64, 0u64);
                    for k in 0..w {
                        lm += mass[b - 1 - k];
                        lc += count[b - 1 - k];
                        rm += mass[b + k];
                        rc += count[b + k];
                    }
                    if lc < min_count || rc < min_count {
                        continue;
                    }
                    let left = lm / lc as f64;
                    let right = rm / rc as f64;
                    let kind = if left <= 0.02 && right > 0.1 {
                        "rise"
                    } else if right <= 0.02 && left > 0.1 {
                        "fall"
                    } else if right < left * 0.6 {
                        "drop"
                    } else if right > left * 1.6 && left > 0.02 {
                        "rise"
                    } else {
                        continue;
                    };
                    bin_transitions.push(BinTransition {
                        source: path.clone(),
                        component: component.clone(),
                        kind: kind.to_string(),
                        profile: profile.to_string(),
                        junction_interval_start: (b.saturating_sub(1)
                            * options.bin_size as usize) as u64,
                        junction_interval_end: ((b + w) * options.bin_size as usize) as u64,
                        level_left: left,
                        level_right: right,
                        bins_left: lc as u32,
                        bins_right: rc as u32,
                    });
                }
            };
            detect("unique", &umass, &ucount, 8, 4);
            detect("all", &mass, &count, 20, 3);
        }
    }
    guard_rss(rss_budget, "bin profiles")?;
    eprintln!("bin transitions: {}", bin_transitions.len());

    let report = Report {
        model: "evidence-topology-scout-v1".into(),
        label: options.label.clone(),
        panel: options.panel.clone(),
        sample: options
            .sample
            .to_string_lossy()
            .into_owned(),
        read_length: 150,
        histogram,
        universe: UniverseReport {
            observed_distinct_pairs: observed.len() as u64,
            observed_total_mass,
            catalog_features: art.catalog_pairs.len() as u64,
            observed_covered_by_catalog: covered_by_catalog,
            observed_in_panel_not_catalog: in_panel_not_catalog.len() as u64,
            observed_not_in_panel: not_in_panel.len() as u64,
            not_in_panel_examples: not_in_panel.iter().take(32).copied().collect(),
            panel_walk_paths: walked_paths,
            structural_novelty_note: "The sample index stores canonical subwalks of maximal MEM records; MEMs match panel paths by construction, so 3-token observations outside the panel universe are structurally impossible. Novel-junction adjacency is not representable in this index: junctions appear as MEM termination (depth steps), not as novel graph pairs.".into(),
        },
        median_ref_rho,
        loci,
        bin_transitions,
        bin_profiles,
        query_accounting: serde_json::json!({
            "count_queries": 0,
            "observed_pair_enumerations": 1,
            "note": "C_f values come from one observed_pairs() enumeration joined against the catalog; no per-feature BWT queries.",
        }),
    };
    let bytes = serde_json::to_vec_pretty(&report).map_err(io::Error::other)?;
    fs::write(options.out_dir.join("report.json"), &bytes)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&serde_json::json!({
            "model": report.model,
            "label": report.label,
            "observed_pairs": report.universe.observed_distinct_pairs,
            "not_in_panel": report.universe.observed_not_in_panel,
            "median_ref_rho": report.median_ref_rho,
            "loci": report.loci.len(),
            "bin_transitions": report.bin_transitions.len(),
            "report": options.out_dir.join("report.json").to_string_lossy(),
        }))
        .map_err(io::Error::other)?
    );
    Ok(())
}
