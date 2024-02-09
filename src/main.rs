// main.rs

use clap::Parser;
use std::fs::File;
use std::io::{self, BufReader};
use impg::{Impg, parse_query};

mod paf;
mod impg;

/// Command line tool for querying PAF file overlaps
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Path to the PAF file
    #[clap(short, long, value_parser)]
    paf_file: String,

    /// Query in the format seq_name:start-end
    #[clap(short, long, value_parser)]
    query: String,

    /// Enable transitive overlap queries
    #[clap(long, action)]
    transitive: bool,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    let (target_name, query_range) = parse_query(&args.query).expect("Invalid query format");
    let (query_start, query_end) = query_range;

    let file = File::open(&args.paf_file)?;
    let reader = BufReader::new(file);

    let records = paf::parse_paf(reader).expect("Failed to parse PAF records");
    let impg = Impg::from_paf_records(&records).expect("Failed to create index");

    let overlaps = if args.transitive {
        impg.query_transitive(&target_name, query_start, query_end)
    } else {
        impg.query(&target_name, query_start, query_end)
    };

    println!("Overlapping regions:");
    for overlap in overlaps {
        println!("{:?}", overlap);
        //println!("{}", overlap);
    }

    Ok(())
}
