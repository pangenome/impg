//! THE FRAME-OFFSET PROBE (the global-panel call mode's convention
//! oracle): does the AGC slice fetched at the record header's interval
//! align with the global panel's own walk over that interval? The
//! canonical anchor walk of the fetched content must place its anchors
//! at exactly the panel walk's steps (position + node) — any constant
//! offset names the convention, fail-closed.
use clap::Parser;
use impg::agc_index::AgcIndex;
use impg::genome_inference::mem_records;
use impg::syng::{SyncmerParams, SyngIndex};
use std::collections::BTreeMap;

#[derive(Parser)]
struct Options {
    #[arg(long)]
    panel: String,
    #[arg(long)]
    agc: std::path::PathBuf,
    /// The record header "path:s-e" (the territory-local frame).
    #[arg(long)]
    record: String,
}

fn main() -> std::io::Result<()> {
    let options = Options::parse();
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())?;
    let agc = AgcIndex::build_from_files(&[options.agc.display().to_string()])?;
    let (name, range) = options
        .record
        .rsplit_once(':')
        .ok_or_else(|| std::io::Error::other("no interval"))?;
    let (s, e): (u64, u64) = range
        .split_once('-')
        .and_then(|(a, b)| Some((a.parse().ok()?, b.parse().ok()?)))
        .ok_or_else(|| std::io::Error::other("bad interval"))?;
    let path_idx = panel
        .name_map
        .name_to_path
        .get(name)
        .copied()
        .ok_or_else(|| std::io::Error::other("path absent"))? as usize;
    let length = 4000u64;
    // fetch the content under BOTH candidate conventions: 1-based
    // half-open (delta=1 -> 0-based [s, s+len)) and 0-based half-open
    // (delta=0 -> 0-based [s, s+len) means fetch(s, s+len)).
    let variants: [(&str, i32, i32); 2] = [
        ("delta=1 (fetch 1-based half-open)", (s + 1) as i32, (s + 1 + length) as i32),
        ("delta=0 (fetch 0-based half-open)", s as i32, (s + length) as i32),
    ];
    for (label, lo, hi) in variants {
        let content = agc
            .fetch_sequence(name, lo, hi)
            .map_err(|e| std::io::Error::other(format!("{label}: {e}")))?;
        if content.len() != length as usize {
            eprintln!("{label}: fetched {} bases (expected {length})", content.len());
            continue;
        }
        let walk = mem_records::canonical_anchor_walk(&panel, &content)?;
        let steps: BTreeMap<u64, i32> = panel
            .walk_path_range(path_idx, s, s + length)?
            .into_iter()
            .map(|(node, bp)| (bp, node))
            .collect();
        let mut exact = 0usize;
        let mut off_by_one = 0usize;
        let mut absent = 0usize;
        for &(node, pos) in &walk {
            let global = s + pos;
            if steps.get(&global) == Some(&node) {
                exact += 1;
            } else if steps.get(&(global + 1)) == Some(&node)
                || steps.get(&(global.wrapping_sub(1))) == Some(&node)
            {
                off_by_one += 1;
            } else {
                absent += 1;
            }
        }
        eprintln!(
            "{label}: anchors {} | exact {} | off-by-one {} | absent {}",
            walk.len(),
            exact,
            off_by_one,
            absent
        );
    }
    Ok(())
}
