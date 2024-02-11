use clap::Parser;
use std::fs::File;
use std::io::{self, BufReader};
use impg::Impg;

mod paf;
mod impg;

/// Command-line tool for querying overlaps in PAF files.
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Path to the PAF file.
    #[clap(short, long, value_parser)]
    paf_file: String,

    /// Query in the format `seq_name:start-end`.
    #[clap(short, long, value_parser)]
    query: String,

    /// Enable transitive overlap queries.
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
    let mut seq_index = impg::SequenceIndex::new();
    for record in &records {
        seq_index.get_or_insert_id(&record.query_name);
        seq_index.get_or_insert_id(&record.target_name);
    }

    let target_id = seq_index.get_id(&target_name).expect("Target name not found in index");
    let impg = Impg::from_paf_records(&records, &seq_index).expect("Failed to create index");

    let overlaps = if args.transitive {
        impg.query_transitive(target_id, query_start, query_end)
    } else {
        impg.query(target_id, query_start, query_end)
    };

    println!("Overlapping regions:");
    for overlap in overlaps {
        println!("{:?}", overlap);
    }

    Ok(())
}

/// Parses the query string into a target sequence name and a range.
fn parse_query(query: &str) -> Result<(String, (i32, i32)), &'static str> {
    let parts: Vec<&str> = query.split(':').collect();
    if parts.len() != 2 {
        return Err("Query format should be `seq_name:start-end`");
    }
    let range_parts: Vec<&str> = parts[1].split('-').collect();
    if range_parts.len() != 2 {
        return Err("Query range format should be `start-end`");
    }

    let start = range_parts[0].parse::<i32>().map_err(|_| "Invalid start value")?;
    let end = range_parts[1].parse::<i32>().map_err(|_| "Invalid end value")?;

    Ok((parts[0].to_string(), (start, end)))
}
