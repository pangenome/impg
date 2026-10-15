//! Restricted finite diagnostic only; no production integration or truth inputs.
#[path = "panel_route_diploid_diagnostic/mod.rs"]
mod diagnostic;
use clap::Parser;
use std::path::PathBuf;
#[derive(Parser)]
#[command(about = "Fixed-two-copy finite diagnostic, not automatic diploid inference")]
struct Args {
    #[arg(long)]
    panel: PathBuf,
    #[arg(long)]
    routes: PathBuf,
    #[arg(long)]
    sample: PathBuf,
    /// JSON array of supplied complete public haploid assignments (max8).
    #[arg(long)]
    domain: PathBuf,
    /// Fresh output; partial evidence survives failures.
    #[arg(long)]
    out: PathBuf,
}
fn main() -> std::io::Result<()> {
    let a = Args::parse();
    let result = diagnostic::run(&a.panel, &a.routes, &a.sample, &a.domain, &a.out)?;
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}
