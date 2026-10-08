//! THE SCORING LAYER (slice B of the local-realignment evidence model;
//! owner direction 2026-11-06, slice A's anchor-projection layer stands
//! as the placement machinery of record). Assessment-side only.
//! THE THIN WRAPPER: the implementation lives in the library
//! (src/genome_inference/realign_score.rs) so the product CLI command
//! path and this instrument example share ONE implementation; this
//! wrapper preserves the committed example CLI surface byte-for-byte.

use clap::Parser;
use impg::genome_inference::realign_score::{self as instrument, Options};
use std::io;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    options: Options,
}

fn main() -> io::Result<()> {
    instrument::run(Cli::parse().options)
}
