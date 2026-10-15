//! Assessment-side feature decomposition for two public port-matched seam cuts.
#[path = "panel_route_diploid_search/mod.rs"]
mod search;

use clap::Parser;
use impg::{
    genome_inference::{panel_routes as routes, read_json, sample, PanelIdentity},
    syng::{SyncmerParams, SyngIndex},
};
use search::{genome_wide, partition::ScoreModel};
use serde::Serialize;
use std::{collections::BTreeMap, io};

#[derive(Parser)]
struct Options {
    #[arg(long)]
    panel: String,
    #[arg(long)]
    sample: std::path::PathBuf,
    #[arg(long)]
    routes: std::path::PathBuf,
    #[arg(long)]
    left_source: usize,
    #[arg(long)]
    left_start: u64,
    #[arg(long)]
    left_cut_a: u64,
    #[arg(long)]
    right_source: usize,
    #[arg(long)]
    right_cut_a: u64,
    #[arg(long)]
    right_end: u64,
    #[arg(long)]
    left_cut_b: u64,
    #[arg(long)]
    right_cut_b: u64,
    #[arg(long, default_value_t = 10.0)]
    depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
    #[arg(long)]
    full: bool,
}

#[derive(Serialize)]
struct Row {
    feature: Vec<u64>,
    predicted_a: u64,
    predicted_b: u64,
    observed: u64,
    signal_a: f64,
    signal_b: f64,
    loss_a: f64,
    loss_b: f64,
    a_minus_b: f64,
}

fn main() -> io::Result<()> {
    let o = Options::parse();
    let identity = PanelIdentity::read(&o.panel)?;
    let panel = SyngIndex::load(&o.panel, SyncmerParams::default())?;
    let sample = sample::SampleIndex::load(&o.sample, &identity)?;
    let histogram = *sample
        .stats
        .read_lengths
        .get(&150)
        .ok_or_else(|| io::Error::other("sample lacks L150"))?;
    let model = ScoreModel {
        read_length: 150,
        histogram,
        denominator: (150 * histogram) as f64,
        depth: o.depth,
        background: o.background,
    };
    let graph: routes::Graph = read_json(&o.routes.join("graph.json"))?;
    let sources = routes::Sources::open(
        &graph.source_paths,
        graph
            .lanes
            .iter()
            .map(|lane| (lane.name.clone(), lane.length))
            .collect(),
    )?;
    let sequence = |source, start, end| sources.fetch(source, start, end);
    let left_a = sequence(o.left_source, o.left_start, o.left_cut_a)?;
    let right_a = sequence(o.right_source, o.right_cut_a, o.right_end)?;
    let left_b = sequence(o.left_source, o.left_start, o.left_cut_b)?;
    let right_b = sequence(o.right_source, o.right_cut_b, o.right_end)?;
    let (a, _) = if o.full {
        genome_wide::profile_event_interior(
            &panel,
            &[left_a.as_slice(), right_a.as_slice()].concat(),
            150,
            genome_wide::MAX_FEATURES,
        )?
    } else {
        genome_wide::profile_event_seam(&panel, &left_a, &right_a, 150, genome_wide::MAX_FEATURES)?
    };
    let (b, _) = if o.full {
        genome_wide::profile_event_interior(
            &panel,
            &[left_b.as_slice(), right_b.as_slice()].concat(),
            150,
            genome_wide::MAX_FEATURES,
        )?
    } else {
        genome_wide::profile_event_seam(&panel, &left_b, &right_b, 150, genome_wide::MAX_FEATURES)?
    };
    let mut keys = a.keys().chain(b.keys()).cloned().collect::<Vec<_>>();
    keys.sort();
    keys.dedup();
    let mut rows = Vec::new();
    for feature in keys {
        let qa = a.get(&feature).copied().unwrap_or(0);
        let qb = b.get(&feature).copied().unwrap_or(0);
        let observed = sample.counts.count(&feature)?;
        let la = model.loss(qa, observed)?;
        let lb = model.loss(qb, observed)?;
        rows.push(Row {
            feature,
            predicted_a: qa,
            predicted_b: qb,
            observed,
            signal_a: qa as f64 * o.depth / 150.0,
            signal_b: qb as f64 * o.depth / 150.0,
            loss_a: la,
            loss_b: lb,
            a_minus_b: la - lb,
        });
    }
    rows.sort_by(|a, b| b.a_minus_b.abs().total_cmp(&a.a_minus_b.abs()));
    let sum_a: f64 = rows.iter().map(|r| r.loss_a).sum();
    let sum_b: f64 = rows.iter().map(|r| r.loss_b).sum();
    let categories = rows.iter().fold(
        BTreeMap::<String, (usize, u64, u64, f64)>::new(),
        |mut m, r| {
            let key = match (r.predicted_a > 0, r.predicted_b > 0, r.observed > 0) {
                (true, false, true) => "a_only_observed",
                (true, false, false) => "a_only_unobserved",
                (false, true, true) => "b_only_observed",
                (false, true, false) => "b_only_unobserved",
                (true, true, true) => "both_observed",
                (true, true, false) => "both_unobserved",
                _ => "neither",
            }
            .to_string();
            let e = m.entry(key).or_default();
            e.0 += 1;
            e.1 += r.observed;
            e.2 += r.predicted_a.abs_diff(r.predicted_b);
            e.3 += r.a_minus_b;
            m
        },
    );
    println!(
        "{}",
        serde_json::to_string_pretty(
            &serde_json::json!({"model":if o.full {"assessment-full-feature-decomposition-v1"} else {"assessment-seam-feature-decomposition-v1"},"a_loss":sum_a,"b_loss":sum_b,"a_minus_b":sum_a-sum_b,"a_features":a.len(),"b_features":b.len(),"categories":categories,"rows":rows})
        )?
    );
    Ok(())
}
