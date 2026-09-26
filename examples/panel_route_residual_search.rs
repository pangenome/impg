//! Standalone Gate B1: no production CLI, warm starts, truth inputs or emission.
#[path = "panel_route_residual_search/mod.rs"]
mod search;
use clap::Parser;
fn main() -> std::io::Result<()> {
    let result = search::run(search::Options::parse(), None)?;
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}
