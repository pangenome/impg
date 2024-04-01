use clap::Parser;
use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use std::num::NonZeroUsize;
use noodles::bgzf;
use impg::impg::{Impg, SerializableImpg, QueryInterval};
use coitrees::IntervalTree;
use impg::paf;
use rayon::ThreadPoolBuilder;

/// Command-line tool for querying overlaps in PAF files.
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Path to the PAF file. If specified without an index, the tool will look for or generate an associated index file.
    #[clap(short='p', long, value_parser)]
    paf_file: Option<String>,

    /// Path to the index file. Use this to specify a custom index file or to force the use of an index.
    #[clap(short='i', long, value_parser)]
    index_file: Option<String>,

    /// Force the regeneration of the index, even if it already exists.
    #[clap(short='I', long, action)]
    force_reindex: bool,

    /// Query in the format `seq_name:start-end`.
    #[clap(short='q', long, value_parser)]
    query: Option<String>,

    /// Enable transitive overlap queries.
    #[clap(short='x', long, action)]
    transitive: bool,

    /// Print stats about the index.
    #[clap(short='s', long, action)]
    stats: bool,

    /// Number of threads for parallel processing.
    #[clap(short='t', long, value_parser, default_value_t = NonZeroUsize::new(1).unwrap())]
    num_threads: NonZeroUsize,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    // Configure the global thread pool to use the specified number of threads
    ThreadPoolBuilder::new().num_threads(args.num_threads.into()).build_global().unwrap();

    let impg = match args {
        Args { paf_file: Some(paf), index_file: None, force_reindex: false, .. } => load_or_generate_index(&paf, None, args.num_threads)?,
        Args { paf_file: Some(paf), index_file: None, force_reindex: true, .. } => generate_index(&paf, None, args.num_threads)?,
        Args { paf_file: Some(paf), index_file: Some(index), force_reindex: false, .. } => load_or_generate_index(&paf, Some(&index), args.num_threads)?,
        Args { paf_file: Some(paf), index_file: Some(index), force_reindex: true, .. } => generate_index(&paf, Some(&index), args.num_threads)?,
        Args { paf_file: None, index_file: Some(index), .. } => load_index(&index)?,
        _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Either a PAF file or an index file must be provided")),
    };

    if args.stats {
        print_stats(&impg);
    }

    if let Some(query) = args.query {
        let (target_name, target_range) = parse_query(&query)?;
        let results = perform_query(&impg, &target_name, target_range, args.transitive);
        output_results(&impg, results);
    }

    Ok(())
}


fn load_or_generate_index(paf_file: &str, index_file: Option<&str>, num_threads: NonZeroUsize) -> io::Result<Impg> {
    let index_file = index_file.map(|s| s.to_string());
    let index_file = index_file.unwrap_or_else(|| format!("{}.impg", paf_file));
    let index_file = index_file.as_str();
    if std::path::Path::new(index_file).exists() {
        load_index(index_file)
    } else {
        generate_index(paf_file, Some(index_file), num_threads)
    }
}

fn generate_index(paf_file: &str, index_file: Option<&str>, num_threads: NonZeroUsize) -> io::Result<Impg> {
    let file = File::open(paf_file)?;
    let reader: Box<dyn io::Read> = if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
        Box::new(bgzf::MultithreadedReader::with_worker_count(num_threads, file))
    } else {
        Box::new(file)
    };
    let reader = BufReader::new(reader);
    let records = paf::parse_paf(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to parse PAF records: {:?}", e)))?;
    let impg = Impg::from_paf_records(&records, paf_file).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to create index: {:?}", e)))?;

    if let Some(index_file) = index_file {
        let serializable = impg.to_serializable();
        let file = File::create(index_file)?;
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &serializable).map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Failed to serialize index: {:?}", e)))?;
    }

    Ok(impg)
}

fn load_index(index_file: &str) -> io::Result<Impg> {
    let file = File::open(index_file)?;
    let reader = BufReader::new(file);
    let serializable: SerializableImpg = bincode::deserialize_from(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to deserialize index: {:?}", e)))?;
    Ok(Impg::from_serializable(serializable))
}

fn parse_query(query: &str) -> io::Result<(String, (i32, i32))> {
    let parts: Vec<&str> = query.split(':').collect();
    if parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Query format should be `seq_name:start-end`"));
    }
    let range_parts: Vec<&str> = parts[1].split('-').collect();
    if range_parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Query range format should be `start-end`"));
    }

    let start = range_parts[0].parse::<i32>().map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid start value"))?;
    let end = range_parts[1].parse::<i32>().map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid end value"))?;

    if start >= end {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Start value must be less than end value"));
    }

    Ok((parts[0].to_string(), (start, end)))
}

fn perform_query(impg: &Impg, target_name: &str, target_range: (i32, i32), transitive: bool) -> Vec<QueryInterval> {
    let (target_start, target_end) = target_range;
    let target_id = impg.seq_index.get_id(target_name).expect("Target name not found in index");
    if transitive {
        impg.query_transitive(target_id, target_start, target_end)
    } else {
        impg.query(target_id, target_start, target_end)
    }
}

fn output_results(impg: &Impg, results: Vec<QueryInterval>) {
    for overlap in results {
        println!("{}\t{}\t{}",
                 impg.seq_index.get_name(overlap.metadata).unwrap(),
                 overlap.first, overlap.last);
    }
}

fn print_stats(impg: &Impg) {
    println!("Number of sequences: {}", impg.seq_index.len());
    println!("Number of overlaps: {}", impg.trees.values().map(|tree| tree.len()).sum::<usize>());
}
