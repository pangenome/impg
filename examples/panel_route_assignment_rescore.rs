//! External-run full-subwalk rescore for two complete public Assignments.
#[path = "panel_route_diploid_search/mod.rs"]
mod search;

use clap::Parser;
use impg::{
    genome_inference::{panel_routes as routes, read_json, sample, PanelIdentity},
    graph::reverse_complement,
    syng::{SyncmerParams, SyngIndex},
};
use search::{genome_wide as genome, partition::ScoreModel};
use std::{fs, io, path::PathBuf};

#[derive(Parser)]
struct Options {
    #[arg(long)]
    panel: String,
    #[arg(long)]
    routes: PathBuf,
    #[arg(long)]
    sample: PathBuf,
    #[arg(long)]
    assignment_a: PathBuf,
    #[arg(long)]
    assignment_b: PathBuf,
    #[arg(long)]
    out_dir: PathBuf,
    #[arg(long, default_value_t = 10.0)]
    depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
}

fn main() -> io::Result<()> {
    let options = Options::parse();
    fs::create_dir_all(&options.out_dir)?;
    let identity = PanelIdentity::read(&options.panel)?;
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    let graph: routes::Graph = read_json(&options.routes.join("graph.json"))?;
    if graph.panel != identity {
        return Err(io::Error::other("panel identity mismatch"));
    }
    let assignments: [routes::Assignment; 2] = [
        read_json(&options.assignment_a)?,
        read_json(&options.assignment_b)?,
    ];
    let checksum = graph.digest()?;
    for assignment in &assignments {
        if assignment.graph_checksum != checksum
            || assignment.family >= graph.families.len()
            || assignment.routes.len() != graph.families[assignment.family].paths.len()
        {
            return Err(io::Error::other("assignment identity/shape mismatch"));
        }
    }
    let sources = routes::Sources::open(
        &graph.source_paths,
        graph
            .lanes
            .iter()
            .map(|lane| (lane.name.clone(), lane.length))
            .collect(),
    )?;
    let sample = sample::SampleIndex::load(&options.sample, &identity)?;
    let histogram = *sample
        .stats
        .read_lengths
        .get(&150)
        .ok_or_else(|| io::Error::other("sample lacks L150 histogram"))?;
    let model = ScoreModel {
        read_length: 150,
        histogram,
        denominator: (150 * histogram) as f64,
        depth: options.depth,
        background: options.background,
    };
    let components = assignments[0].routes.len();
    if assignments[1].routes.len() != components {
        return Err(io::Error::other("assignment component mismatch"));
    }
    let mut assembled = Vec::<[Vec<u8>; 2]>::with_capacity(components);
    for component in 0..components {
        let mut copies = [Vec::new(), Vec::new()];
        for copy in 0..2 {
            for segment in &assignments[copy].routes[component].segments {
                let mut sequence = sources.fetch(segment.source, segment.start, segment.end)?;
                if segment.reverse {
                    sequence = reverse_complement(&sequence);
                }
                copies[copy].extend(sequence);
            }
        }
        assembled.push(copies);
    }
    let views = assembled
        .iter()
        .map(|copies| [copies[0].as_slice(), copies[1].as_slice()])
        .collect::<Vec<_>>();
    let runs = options.out_dir.join("external-runs");
    if runs.exists() {
        fs::remove_dir_all(&runs)?;
    }
    let (loss, cost, initial_runs) =
        genome::external_rescore_components(&panel, &views, &sample.counts, &model, &runs)?;
    let result = serde_json::json!({
        "model": "experimental-assignment-pair-full-subwalk-rescore-v1",
        "joint_full_subwalk_loss": loss,
        "profile_work": cost.mem_queries,
        "integrated_windows": cost.integrated_windows,
        "external_initial_runs": initial_runs,
        "graph_checksum": checksum,
    });
    fs::write(
        options.out_dir.join("result.json"),
        serde_json::to_vec_pretty(&result).map_err(io::Error::other)?,
    )?;
    println!(
        "{}",
        serde_json::to_string_pretty(&result).map_err(io::Error::other)?
    );
    Ok(())
}
