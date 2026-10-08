//! POOLED-DIVERGENCE DECOMPOSITION AT SEQUENCE LEVEL (owner direction,
//! genome/finalist-reranking close-out): map where the pooled external's
//! truth preference lives. The pooled scorer is SEPARABLE over features
//! (total = Sum_f model.loss(predicted_f, observed_f), predicted_f = the
//! feature's count summed over every L150 window of the spelled sequence),
//! so a region's pooled contribution is measured EXACTLY by truncating the
//! spelled sequence at that region's boundary and differencing the per-
//! feature losses. No model change: a measurement over the run's own
//! conventions (the assessment concat: the truth route's segments in file
//! order; the selected chain's per-locus alleles in locus order).
#![recursion_limit = "512"]
#[path = "panel_route_diploid_search/mod.rs"]
mod search;

use clap::Parser;
use impg::{
    genome_inference::{panel_routes as routes, read_json, sample, PanelIdentity},
    graph::reverse_complement,
    syng::{SyncmerParams, SyngIndex},
};
use search::genome_wide as genome;
use search::partition::{FeatureKey, ScoreModel};
use std::collections::BTreeMap;
use std::io::{self};
use std::path::PathBuf;

const READ_LENGTH: usize = 150;
const MAX_FEATURES: usize = genome::MAX_FEATURES;

#[derive(Parser)]
struct Options {
    #[arg(long)]
    panel: String,
    #[arg(long)]
    routes: PathBuf,
    #[arg(long)]
    sample: PathBuf,
    #[arg(long, default_value_t = 10.0)]
    depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
    /// The SELECTED chain's route (the rescore `selected_route` format:
    /// per-locus [slot-0 {segments}, "empty-slot2"]); spelled in LOCUS order.
    #[arg(long)]
    selected_route: PathBuf,
    /// The TRUTH route (a `routes::Route` JSON); spelled in FILE order (the
    /// assessment convention its pooled ladder number used).
    #[arg(long)]
    truth_route: PathBuf,
    /// Truncate the selected spelled sequence at SOURCE:COORD (its spelled
    /// material beyond the coordinate is dropped; later segments are not
    /// spelled). Repeatable; each truncation is scored separately.
    #[arg(long = "truncate-selected")]
    truncate_selected: Vec<String>,
    /// Truncate the truth spelled sequence at SOURCE:COORD.
    #[arg(long = "truncate-truth")]
    truncate_truth: Vec<String>,
    #[arg(long)]
    out: PathBuf,
}

fn invalid(s: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s)
}

/// A spelled sequence as an ordered list of pieces with the concat's byte
/// offsets, so a truncation at a piece coordinate maps to a byte offset.
struct SpelledSequence {
    pieces: Vec<(usize, u64, u64, bool, usize, usize)>, // source, start, end, rev, byte_lo, byte_hi
    bases: Vec<u8>,
}

impl SpelledSequence {
    /// The bases up to `coordinate` on `source` (end-exclusive): everything
    /// before the coordinate's piece keeps its bytes; the coordinate's piece
    /// is cut; later pieces are dropped.
    fn truncate_at(&self, source: usize, coordinate: u64) -> io::Result<Vec<u8>> {
        for &(piece_source, start, end, _, byte_lo, byte_hi) in &self.pieces {
            if piece_source == source && coordinate > start && coordinate <= end {
                let within = (coordinate - start) as usize;
                let cut = byte_lo + within;
                if cut > byte_hi {
                    return Err(invalid("truncation cut outside its piece"));
                }
                return Ok(self.bases[..cut].to_vec());
            }
        }
        Err(invalid("truncation coordinate inside no spelled piece"))
    }
}

/// The per-feature predicted counts over the whole spelled sequence: the
/// L150 windows in batches (the pooled machinery's own batching), the
/// features accumulated IN MEMORY (the external merge is a plain sum).
fn accumulated_profile(
    panel: &SyngIndex,
    sequence: &[u8],
) -> io::Result<BTreeMap<FeatureKey, u64>> {
    let read_length = READ_LENGTH;
    let starts = sequence.len().saturating_sub(read_length) + 1;
    let mut total: BTreeMap<FeatureKey, u64> = BTreeMap::new();
    let mut lo = 0usize;
    while lo < starts {
        let hi = (lo + 512).min(starts);
        let (profile, _) =
            genome::profile_event_runs(panel, sequence, read_length, lo, hi, MAX_FEATURES)?;
        for (key, count) in profile {
            *total.entry(key).or_default() += count;
        }
        lo = hi;
    }
    Ok(total)
}

/// The pooled loss over an accumulated profile, ascending key order (the
/// external scorer's own order — bit-identical totals).
fn pooled_loss(
    profile: &BTreeMap<FeatureKey, u64>,
    sample_counts: &impg::sample_mem_bwt::WeightedBwt,
    model: &ScoreModel,
) -> io::Result<f64> {
    let mut loss = 0.0f64;
    for (key, predicted) in profile {
        loss += model.loss(*predicted, sample_counts.count(key)?)?;
    }
    Ok(loss)
}

fn parse_truncation(spec: &str) -> io::Result<(usize, u64)> {
    let mut parts = spec.splitn(2, ':');
    let source = parts
        .next()
        .and_then(|value| value.parse::<usize>().ok())
        .ok_or_else(|| invalid("--truncate needs SOURCE:COORD"))?;
    let coordinate = parts
        .next()
        .and_then(|value| value.parse::<u64>().ok())
        .ok_or_else(|| invalid("--truncate needs SOURCE:COORD"))?;
    Ok((source, coordinate))
}

fn main() -> io::Result<()> {
    let started = std::time::Instant::now();
    let options = Options::parse();
    let identity = PanelIdentity::read(&options.panel)?;
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    let graph: routes::Graph = read_json(&options.routes.join("graph.json"))?;
    let lanes = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.clone(), lane.length))
        .collect::<Vec<_>>();
    let sources = routes::Sources::open(&graph.source_paths, lanes)?;
    let sample_index = sample::SampleIndex::load(&options.sample, &identity)?;
    let histogram = *sample_index
        .stats
        .read_lengths
        .get(&READ_LENGTH)
        .ok_or_else(|| invalid("sample lacks L150 histogram"))?;
    let model = ScoreModel {
        read_length: READ_LENGTH as u64,
        histogram,
        denominator: (READ_LENGTH as u64 * histogram) as f64,
        depth: options.depth,
        background: options.background,
    };
    model.validate()?;

    // The selected chain: per-locus slot-0 alleles, spelled in LOCUS order.
    let selected_json: serde_json::Value = read_json(&options.selected_route)?;
    let selected_rows = selected_json
        .as_array()
        .ok_or_else(|| invalid("selected route must be an array"))?;
    let mut selected_pieces: Vec<(usize, u64, u64, bool)> = Vec::new();
    for row in selected_rows {
        let slots = row
            .as_array()
            .ok_or_else(|| invalid("selected route row must be an array"))?;
        if slots.len() != 2 || slots[1].as_str() != Some("empty-slot2") {
            return Err(invalid(
                "selected route row must be a haploid slot-0/empty-slot2 pair",
            ));
        }
        for segment in slots[0]
            .get("segments")
            .and_then(serde_json::Value::as_array)
            .ok_or_else(|| invalid("selected route row lacks segments"))?
        {
            selected_pieces.push((
                segment
                    .get("source")
                    .and_then(serde_json::Value::as_u64)
                    .ok_or_else(|| invalid("segment lacks source"))? as usize,
                segment
                    .get("start")
                    .and_then(serde_json::Value::as_u64)
                    .ok_or_else(|| invalid("segment lacks start"))?,
                segment
                    .get("end")
                    .and_then(serde_json::Value::as_u64)
                    .ok_or_else(|| invalid("segment lacks end"))?,
                segment
                    .get("reverse")
                    .and_then(serde_json::Value::as_bool)
                    .ok_or_else(|| invalid("segment lacks reverse"))?,
            ));
        }
    }
    // The truth: the route's own segments in FILE order.
    let truth_route: routes::Route = read_json(&options.truth_route)?;
    let truth_pieces: Vec<(usize, u64, u64, bool)> = truth_route
        .segments
        .iter()
        .filter(|segment| segment.start < segment.end)
        .map(|segment| (segment.source, segment.start, segment.end, segment.reverse))
        .collect();

    let spell = |pieces: &[(usize, u64, u64, bool)]| -> io::Result<SpelledSequence> {
        let mut bases = Vec::new();
        let mut layout = Vec::new();
        for &(source, start, end, reverse) in pieces {
            let mut part = sources.fetch(source, start, end)?;
            if reverse {
                part = reverse_complement(&part);
            }
            let byte_lo = bases.len();
            bases.extend_from_slice(&part);
            layout.push((source, start, end, reverse, byte_lo, bases.len()));
        }
        Ok(SpelledSequence {
            pieces: layout,
            bases,
        })
    };
    let selected_sequence = spell(&selected_pieces)?;
    let truth_sequence = spell(&truth_pieces)?;
    eprintln!(
        "[pooled-decompose] spelled: selected {} bases ({} pieces), truth {} bases ({} pieces) ({:.1}s)",
        selected_sequence.bases.len(),
        selected_pieces.len(),
        truth_sequence.bases.len(),
        truth_pieces.len(),
        started.elapsed().as_secs_f64()
    );

    // The two full-length anchors (must reproduce the run's pooled numbers
    // bit-exactly: the same accumulation, the same ascending-key loss order).
    let selected_profile = accumulated_profile(&panel, &selected_sequence.bases)?;
    let truth_profile = accumulated_profile(&panel, &truth_sequence.bases)?;
    let selected_full = pooled_loss(&selected_profile, &sample_index.counts, &model)?;
    let truth_full = pooled_loss(&truth_profile, &sample_index.counts, &model)?;
    eprintln!(
        "[pooled-decompose] anchors: selected {selected_full:.6} truth {truth_full:.6} ({:.1}s)",
        started.elapsed().as_secs_f64()
    );

    // The truncation experiments: each region's pooled contribution = the
    // full-minus-truncated difference, with the per-feature anatomy (the
    // features whose predicted counts changed).
    let experiment = |name: &str,
                      sequence: &SpelledSequence,
                      full_profile: &BTreeMap<FeatureKey, u64>,
                      full_loss: f64,
                      spec: &str|
     -> io::Result<serde_json::Value> {
        let (source, coordinate) = parse_truncation(spec)?;
        let truncated_bases = sequence.truncate_at(source, coordinate)?;
        let truncated_profile = accumulated_profile(&panel, &truncated_bases)?;
        let truncated_loss = pooled_loss(&truncated_profile, &sample_index.counts, &model)?;
        // The region's features: the keys whose predicted counts dropped.
        let mut rows: Vec<(FeatureKey, u64, u64, f64, f64)> = Vec::new();
        let mut region_penalty = 0.0f64;
        let mut dropped_features = 0usize;
        let mut weakened_features = 0usize;
        for (key, full_count) in full_profile {
            let truncated_count = truncated_profile.get(key).copied().unwrap_or(0);
            if truncated_count == *full_count {
                continue;
            }
            let observed = sample_index.counts.count(key)?;
            let full_term = model.loss(*full_count, observed)?;
            let truncated_term = model.loss(truncated_count, observed)?;
            let delta = truncated_term - full_term;
            region_penalty += delta;
            if truncated_count == 0 {
                dropped_features += 1;
            } else {
                weakened_features += 1;
            }
            rows.push((key.clone(), *full_count, truncated_count, full_term, delta));
        }
        rows.sort_by(|left, right| right.4.total_cmp(&left.4));
        let shown: Vec<serde_json::Value> = rows
            .iter()
            .take(24)
            .map(|(key, full_count, truncated_count, full_term, delta)| {
                serde_json::json!({
                    "feature": key,
                    "full_count": full_count,
                    "truncated_count": truncated_count,
                    "observed": sample_index.counts.count(key).unwrap_or(0),
                    "full_term": full_term,
                    "region_penalty_contribution": delta,
                })
            })
            .collect();
        eprintln!(
            "[pooled-decompose] {name}: full {full_loss:.4} truncated {truncated_loss:.4} \
             region penalty {region_penalty:.4} (dropped {dropped_features}, weakened \
             {weakened_features})",
        );
        Ok(serde_json::json!({
            "name": name,
            "truncation": spec,
            "full_loss": full_loss,
            "truncated_loss": truncated_loss,
            "region_pooled_penalty": region_penalty,
            "dropped_features": dropped_features,
            "weakened_features": weakened_features,
            "top_contributions": shown,
        }))
    };

    let mut experiments = Vec::new();
    for spec in &options.truncate_selected {
        experiments.push(experiment(
            "selected_truncated",
            &selected_sequence,
            &selected_profile,
            selected_full,
            spec,
        )?);
    }
    for spec in &options.truncate_truth {
        experiments.push(experiment(
            "truth_truncated",
            &truth_sequence,
            &truth_profile,
            truth_full,
            spec,
        )?);
    }

    // The FULL truth-vs-selected per-feature difference table (where the
    // pooled's preference lives, feature by feature).
    let mut keys: Vec<FeatureKey> = selected_profile
        .keys()
        .chain(truth_profile.keys())
        .cloned()
        .collect();
    keys.sort();
    keys.dedup();
    let mut difference_rows: Vec<(FeatureKey, u64, u64, u64, f64, f64)> = Vec::new();
    for key in keys {
        let selected_count = selected_profile.get(&key).copied().unwrap_or(0);
        let truth_count = truth_profile.get(&key).copied().unwrap_or(0);
        if selected_count == truth_count {
            continue;
        }
        let observed = sample_index.counts.count(&key)?;
        let selected_term = model.loss(selected_count, observed)?;
        let truth_term = model.loss(truth_count, observed)?;
        difference_rows.push((
            key,
            selected_count,
            truth_count,
            observed,
            selected_term,
            truth_term,
        ));
    }
    // The truth's preference per feature: truth_term - selected_term
    // (NEGATIVE = the truth is cheaper there = the truth's pooled advantage).
    difference_rows.sort_by(|left, right| {
        (left.5 - left.4).total_cmp(&(right.5 - right.4))
    });
    let truth_advantage_total: f64 = difference_rows
        .iter()
        .map(|&(_, _, _, _, selected_term, truth_term)| truth_term - selected_term)
        .sum();
    let truth_better_count = difference_rows
        .iter()
        .filter(|&&(_, _, _, _, selected_term, truth_term)| truth_term < selected_term)
        .count();
    let selected_better_count = difference_rows
        .iter()
        .filter(|&&(_, _, _, _, selected_term, truth_term)| selected_term < truth_term)
        .count();
    let top_truth_advantage: Vec<serde_json::Value> = difference_rows
        .iter()
        .take(32)
        .map(
            |&(ref key, selected_count, truth_count, observed, selected_term, truth_term)| {
                serde_json::json!({
                    "feature": key,
                    "selected_count": selected_count,
                    "truth_count": truth_count,
                    "observed": observed,
                    "selected_term": selected_term,
                    "truth_term": truth_term,
                    "truth_advantage": truth_term - selected_term,
                })
            },
        )
        .collect();
    eprintln!(
        "[pooled-decompose] truth-vs-selected: {} differing features (truth better on {}, \
         selected better on {}), net truth advantage {truth_advantage_total:.4}",
        difference_rows.len(),
        truth_better_count,
        selected_better_count
    );

    let table = serde_json::json!({
        "what": "pooled-divergence decomposition at sequence level: the pooled external scorer \
            is separable over features, so each region's pooled contribution is the \
            full-minus-truncated difference of the per-feature losses (owner direction, \
            genome/finalist-reranking close-out)",
        "anchors": {
            "selected_full": selected_full,
            "truth_full": truth_full,
            "truth_minus_selected": truth_full - selected_full,
        },
        "experiments": experiments,
        "truth_vs_selected": {
            "differing_features": difference_rows.len(),
            "truth_better_features": truth_better_count,
            "selected_better_features": selected_better_count,
            "net_truth_advantage": truth_advantage_total,
            "top_truth_advantages": top_truth_advantage,
        },
        "wall_seconds": started.elapsed().as_secs_f64(),
    });
    let mut writer = std::fs::File::create(&options.out)?;
    serde_json::to_writer_pretty(&mut writer, &table)?;
    writer.sync_all()?;
    println!(
        "{}",
        serde_json::to_string(&serde_json::json!({
            "anchors": table["anchors"],
            "experiments": table["experiments"],
            "truth_vs_selected_net": truth_advantage_total,
        }))?
    );
    Ok(())
}
