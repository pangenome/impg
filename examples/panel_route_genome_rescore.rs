//! Aggregate chromosome-reset candidates and perform one external-run joint rescore.
#[path = "panel_route_diploid_search/mod.rs"]
mod search;

use clap::Parser;
use impg::{
    genome_inference::{panel_routes as routes, read_json, sample, PanelIdentity},
    graph::reverse_complement,
    syng::{SyncmerParams, SyngIndex},
};
use search::{genome_wide as genome, partition::ScoreModel};
use serde_json::Value;
use std::{collections::BTreeMap, fs, io, path::PathBuf};

#[derive(Parser)]
struct Options {
    #[arg(long)]
    panel: String,
    #[arg(long)]
    routes: PathBuf,
    #[arg(long)]
    sample: PathBuf,
    #[arg(long)]
    components: PathBuf,
    #[arg(long)]
    out_dir: PathBuf,
    #[arg(long, default_value_t = 216)]
    family: usize,
    #[arg(long, default_value_t = 10.0)]
    depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

fn main() -> io::Result<()> {
    let options = Options::parse();
    fs::create_dir_all(&options.out_dir)?;
    let identity = PanelIdentity::read(&options.panel)?;
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    let graph: routes::Graph = read_json(&options.routes.join("graph.json"))?;
    if graph.panel != identity || options.family >= graph.families.len() {
        return Err(invalid("incompatible graph/family"));
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
        .ok_or_else(|| invalid("sample lacks L150 histogram"))?;
    let model = ScoreModel {
        read_length: 150,
        histogram,
        denominator: (150 * histogram) as f64,
        depth: options.depth,
        background: options.background,
    };

    let mut by_component = BTreeMap::<String, ([routes::Route; 2], bool, String)>::new();
    for entry in fs::read_dir(&options.components)? {
        let entry = entry?;
        if !entry.file_type()?.is_dir() || !entry.path().join("run.log").is_file() {
            continue;
        }
        let text = fs::read_to_string(entry.path().join("run.log"))?;
        let start = text
            .find('{')
            .ok_or_else(|| invalid("component JSON missing"))?;
        let value: Value = serde_json::from_str(&text[start..]).map_err(io::Error::other)?;
        let Some(final_result) = value.get("final_result").filter(|value| !value.is_null()) else {
            continue;
        };
        let component = value["component"]
            .as_str()
            .ok_or_else(|| invalid("component name missing"))?
            .to_string();
        let public_routes: [routes::Route; 2] =
            serde_json::from_value(final_result["public_physical_routes"].clone())
                .map_err(io::Error::other)?;
        by_component.insert(
            component,
            (
                public_routes,
                value["complete"].as_bool().unwrap_or(false),
                value["stop_reason"].as_str().unwrap_or("unknown").into(),
            ),
        );
    }

    let family = &graph.families[options.family];
    let mut assignments = [Vec::new(), Vec::new()];
    let mut assembled = Vec::<[Vec<u8>; 2]>::new();
    let mut statuses = Vec::new();
    for &lane in &family.paths {
        let name = graph.lanes[lane].name.clone();
        let (component_routes, complete, stop_reason) = by_component
            .remove(&name)
            .ok_or_else(|| invalid(format!("missing component result for {name}")))?;
        let mut copies = [Vec::new(), Vec::new()];
        for copy in 0..2 {
            assignments[copy].push(component_routes[copy].clone());
            for segment in &component_routes[copy].segments {
                let mut sequence = sources.fetch(segment.source, segment.start, segment.end)?;
                if segment.reverse {
                    sequence = reverse_complement(&sequence);
                }
                copies[copy].extend(sequence);
            }
        }
        assembled.push(copies);
        statuses.push(serde_json::json!({
            "component": name,
            "complete": complete,
            "stop_reason": stop_reason,
        }));
    }
    let views = assembled
        .iter()
        .map(|copies| [copies[0].as_slice(), copies[1].as_slice()])
        .collect::<Vec<_>>();
    let run_directory = options.out_dir.join("external-runs");
    if run_directory.exists() {
        fs::remove_dir_all(&run_directory)?;
    }
    let (loss, cost, initial_runs) = genome::external_rescore_components(
        &panel,
        &views,
        &sample.counts,
        &model,
        &run_directory,
    )?;
    let checksum = graph.digest()?;
    for copy in 0..2 {
        let assignment = routes::Assignment {
            version: 1,
            model: routes::MODEL.into(),
            graph_checksum: checksum.clone(),
            family: options.family,
            routes: assignments[copy].clone(),
        };
        fs::write(
            options.out_dir.join(format!("assignment-{copy}.json")),
            serde_json::to_vec_pretty(&assignment).map_err(io::Error::other)?,
        )?;
    }
    let result = serde_json::json!({
        "model": "experimental-whole-axis-external-rescore-v1",
        "components": statuses,
        "all_components_exact": statuses.iter().all(|value| value["complete"] == true),
        "joint_full_subwalk_loss": loss,
        "profile_work": cost.mem_queries,
        "integrated_windows": cost.integrated_windows,
        "external_initial_runs": initial_runs,
        "graph_checksum": checksum,
        "family": options.family,
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
