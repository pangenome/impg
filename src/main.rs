use clap::Parser;
use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use flate2::read::GzDecoder;
use impg::impg::{Impg, SerializableImpg};
use impg::paf;
//use std::collections::HashMap;
//use bincode;

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
    #[clap(short='x', long, action)]
    transitive: bool,

    /// Path to the impg index file.
    #[clap(short, long, value_parser)]
    index_file: Option<String>,

    /// Force reindexing even if the index file already exists.
    #[clap(short='f', long, action)]
    force_reindex: bool,
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    let (target_name, target_range) = parse_query(&args.query).expect("Invalid query format");
    let (target_start, target_end) = target_range;

    let impg = if let Some(index_file) = &args.index_file {
        if std::path::Path::new(index_file).exists() && !args.force_reindex {
            // Load the index from the file
            let file = File::open(index_file)?;
            let reader = BufReader::new(file);
            let serializable: SerializableImpg = bincode::deserialize_from(reader).expect("Failed to deserialize index");
            Impg::from_serializable(serializable)
        } else {
            // Generate the index from the PAF file
            let file = File::open(&args.paf_file)?;
            let reader: Box<dyn io::Read> = if args.paf_file.ends_with(".gz") {
                Box::new(GzDecoder::new(file))
            } else {
                Box::new(file)
            };
            let reader = BufReader::new(reader);
            let records = paf::parse_paf(reader).expect("Failed to parse PAF records");
            let impg = Impg::from_paf_records(&records).expect("Failed to create index");

            // Serialize the index to the file
            if let Some(index_file) = &args.index_file {
                let serializable = impg.to_serializable();
                let file = File::create(index_file)?;
                let writer = BufWriter::new(file);
                bincode::serialize_into(writer, &serializable).expect("Failed to serialize index");
            }

            impg
        }
    } else {
        // Generate the index from the PAF file
        let file = File::open(&args.paf_file)?;
        let reader: Box<dyn io::Read> = if args.paf_file.ends_with(".gz") {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let reader = BufReader::new(reader);
        let records = paf::parse_paf(reader).expect("Failed to parse PAF records");
        Impg::from_paf_records(&records).expect("Failed to create index")
    };

    let target_id = impg.seq_index.get_id(&target_name).expect("Target name not found in index");
    let overlaps = if args.transitive {
        impg.query_transitive(target_id, target_start, target_end)
    } else {
        impg.query(target_id, target_start, target_end)
    };

    // write output in BED format relative to the queries that match the target range
    println!("{}\t{}\t{}", target_name, target_start, target_end);
    for overlap in overlaps {
        println!("{}\t{}\t{}",
                 impg.seq_index.get_name(overlap.metadata).unwrap(),
                 overlap.first, overlap.last);
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

    // assert that start and end are at least one apart
    if start >= end {
        return Err("Start value must be less than end value");
    }

    Ok((parts[0].to_string(), (start, end)))
}
