//! One-shot experimental completion proposals; never mutates a production snapshot.
#[path = "panel_route_completion_probe/mod.rs"]
mod probe;
use clap::Parser;
fn main() -> std::io::Result<()> {
    probe::run(probe::Options::parse())
}
