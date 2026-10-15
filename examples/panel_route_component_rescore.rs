//! Assessment-side single-component two-copy rescore.
//!
//! Scores one component's two-copy route pair with the authoritative external
//! full-subwalk evaluator (the same `external_rescore` machinery the chain
//! runner uses for finalist rescoring), and reports the production port-word
//! seam legality of each route. Assessment tooling only: no production
//! behavior depends on this example.
#[path = "panel_route_diploid_search/mod.rs"]
mod search;

use clap::Parser;
use impg::{
    genome_inference::{panel_routes as routes, read_json, sample, PanelIdentity},
    graph::reverse_complement,
    syng::{SyncmerParams, SyngIndex},
};
use search::{genome_wide as genome, partition::ScoreModel};
use serde_json::json;
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
    route_a: PathBuf,
    #[arg(long)]
    route_b: PathBuf,
    #[arg(long)]
    out_dir: PathBuf,
    #[arg(long)]
    label: String,
    #[arg(long, default_value_t = 10.0)]
    depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
    /// Report port-word legality only; skip the scoring run.
    #[arg(long)]
    legality_only: bool,
}

#[derive(serde::Deserialize)]
struct RouteFile {
    segments: Vec<routes::Segment>,
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
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
    let checksum = graph.digest()?;
    let sources = routes::Sources::open(
        &graph.source_paths,
        graph
            .lanes
            .iter()
            .map(|lane| (lane.name.clone(), lane.length))
            .collect(),
    )?;
    let mut ports = routes::Ports::open(&options.routes, &graph)?;
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

    let route_files = [
        read_json::<RouteFile>(&options.route_a)?,
        read_json::<RouteFile>(&options.route_b)?,
    ];
    let mut legality = Vec::new();
    let mut assembled = [Vec::new(), Vec::new()];
    for (copy, file) in route_files.iter().enumerate() {
        let segments = &file.segments;
        if segments.is_empty() {
            return Err(io::Error::other("empty route"));
        }
        // Production native endpoint pairing check.
        let first = &segments[0];
        let last = segments.last().expect("nonempty");
        let endpoints_ok = first.start == 0
            && !first.reverse
            && !last.reverse
            && last.end == graph.lanes[last.source].length;
        let mut seam_reports = Vec::new();
        for pair in segments.windows(2) {
            let left = &pair[0];
            let right = &pair[1];
            let left_port = ports.at_cut(&graph, left.source, left.exit(), left.reverse)?;
            let right_port = ports.at_cut(&graph, right.source, right.entry(), right.reverse)?;
            let report = match (left_port, right_port) {
                (Some(l), Some(r)) => {
                    let equal = l.word == r.word;
                    json!({
                        "left": {"source": left.source, "exit": left.exit(), "reverse": left.reverse, "anchor": l.anchor, "word": hex(&l.word)},
                        "right": {"source": right.source, "entry": right.entry(), "reverse": right.reverse, "anchor": r.anchor, "word": hex(&r.word)},
                        "legal": equal,
                    })
                }
                (l, r) => json!({
                    "left_present": l.is_some(),
                    "right_present": r.is_some(),
                    "legal": false,
                }),
            };
            seam_reports.push(report);
        }
        let all_legal = seam_reports
            .iter()
            .all(|report| report["legal"].as_bool() == Some(true));
        legality.push(json!({
            "copy": copy,
            "segments": segments.len(),
            "native_endpoints_ok": endpoints_ok,
            "seams": seam_reports,
            "all_port_word_seams_legal": all_legal,
        }));
        for segment in segments {
            let mut sequence = sources.fetch(segment.source, segment.start, segment.end)?;
            if segment.reverse {
                sequence = reverse_complement(&sequence);
            }
            assembled[copy].extend(sequence);
        }
    }
    let all_legal = legality
        .iter()
        .all(|report| report["all_port_word_seams_legal"].as_bool() == Some(true))
        && legality
            .iter()
            .all(|report| report["native_endpoints_ok"].as_bool() == Some(true));

    let mut result = json!({
        "model": "assessment-component-pair-full-subwalk-rescore-v1",
        "label": options.label,
        "graph_checksum": checksum,
        "legality": legality,
        "all_legal": all_legal,
        "sequence_lengths": [assembled[0].len(), assembled[1].len()],
    });

    if !options.legality_only {
        let runs = options.out_dir.join("external-runs");
        if runs.exists() {
            fs::remove_dir_all(&runs)?;
        }
        let views = [assembled[0].as_slice(), assembled[1].as_slice()];
        let (loss, cost, initial_runs) =
            genome::external_rescore(&panel, views, &sample.counts, &model, &runs)?;
        result["joint_full_subwalk_loss"] = json!(loss);
        result["profile_work"] = json!(cost.mem_queries);
        result["integrated_windows"] = json!(cost.integrated_windows);
        result["external_initial_runs"] = json!(initial_runs);
    }

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
