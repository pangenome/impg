//! No production CLI integration. Only explicitly bounded synthetic fixture input.
#[path = "panel_route_convex_diagnostic/mod.rs"]
mod diagnostic;
use clap::Parser;
use std::path::PathBuf;
#[derive(Parser)]
#[command(
    about = "Gate-A finite synthetic Poisson diagnostic; no sequence emission or genome-wide certificate"
)]
struct Args {
    /// Fresh retained finite513 or coupled64x64 test fixture, never biological data.
    #[arg(long)]
    fixture: PathBuf,
    /// Must not exist; failures retained here.
    #[arg(long)]
    out: PathBuf,
}
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let summary = diagnostic::run_fixture(&args.fixture, &args.out)?;
    println!("{}", serde_json::to_string_pretty(&summary)?);
    Ok(())
}
