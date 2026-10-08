//! Bounded automatic fixed-two-whole-copy experiment, not production inference.
#[path = "panel_route_diploid_search/mod.rs"]
mod search;
use clap::Parser;
fn main() -> std::io::Result<()> {
    let result = search::run(search::Options::parse(), None)?;
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}
