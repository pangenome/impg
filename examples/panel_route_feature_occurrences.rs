//! Assessment-side panel source multiplicity for diagnostic MEM-subwalk features.
use clap::Parser;
use impg::{
    sample_mem_bwt::{encode_walk, reverse_complement},
    syng::{SyncmerParams, SyngIndex},
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io,
    path::PathBuf,
};

#[derive(Parser)]
struct Options {
    #[arg(long)]
    panel: String,
    #[arg(long)]
    diagnostic: PathBuf,
}

#[derive(Deserialize)]
struct Diagnostic {
    rows: Vec<Row>,
}
#[derive(Deserialize)]
struct Row {
    feature: Vec<u64>,
    predicted_a: u64,
    predicted_b: u64,
    observed: u64,
}
#[derive(Serialize)]
struct ResultRow {
    feature: Vec<u64>,
    observed: u64,
    predicted_b: u64,
    panel_occurrences: u64,
    panel_sources: Vec<String>,
}

fn main() -> io::Result<()> {
    let o = Options::parse();
    let panel = SyngIndex::load(&o.panel, SyncmerParams::default())?;
    let diagnostic: Diagnostic =
        serde_json::from_reader(File::open(o.diagnostic)?).map_err(io::Error::other)?;
    let targets = diagnostic
        .rows
        .into_iter()
        .filter(|r| r.predicted_a == 0 && r.predicted_b > 0 && r.observed > 0)
        .collect::<Vec<_>>();
    let mut by_first = BTreeMap::<u64, Vec<(usize, Vec<u64>)>>::new();
    for (i, target) in targets.iter().enumerate() {
        for pattern in [target.feature.clone(), reverse_complement(&target.feature)] {
            by_first.entry(pattern[0]).or_default().push((i, pattern));
        }
    }
    let mut occurrences = vec![0u64; targets.len()];
    let mut sources = vec![BTreeSet::<String>::new(); targets.len()];
    for (path, start) in panel.name_map.path_starts.iter().enumerate() {
        let Some(start) = start else { continue };
        let tape = encode_walk(&panel.walk_forward_path(start))?;
        let mut found = BTreeSet::new();
        for position in (0..tape.len()).step_by(2) {
            let Some(patterns) = by_first.get(&tape[position]) else {
                continue;
            };
            for (target, pattern) in patterns {
                if position + pattern.len() <= tape.len()
                    && tape[position..position + pattern.len()] == pattern[..]
                {
                    occurrences[*target] = occurrences[*target].saturating_add(1);
                    found.insert(*target);
                }
            }
        }
        for target in found {
            sources[target].insert(panel.name_map.path_to_name[path].clone());
        }
    }
    let rows = targets
        .into_iter()
        .enumerate()
        .map(|(i, r)| ResultRow {
            feature: r.feature,
            observed: r.observed,
            predicted_b: r.predicted_b,
            panel_occurrences: occurrences[i],
            panel_sources: sources[i].iter().cloned().collect(),
        })
        .collect::<Vec<_>>();
    println!(
        "{}",
        serde_json::to_string_pretty(
            &serde_json::json!({"model":"assessment-panel-feature-multiplicity-v1","rows":rows})
        )
        .map_err(io::Error::other)?
    );
    Ok(())
}
