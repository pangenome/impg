use std::env;
use std::fs;
use std::io;

use impg::resolution::path_sequences;
use impg::sequence_index::{SequenceIndex, UnifiedSequenceIndex};

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let Some(gfa_path) = args.next() else {
        usage();
    };
    let sequence_files: Vec<String> = args.collect();
    if sequence_files.is_empty() {
        usage();
    }

    let gfa = fs::read_to_string(&gfa_path)?;
    let paths = path_sequences(&gfa)?;
    let sequence_index = UnifiedSequenceIndex::from_files(&sequence_files)?;

    let mut unparsable_paths = 0usize;
    let mut forward_matches = 0usize;
    let mut reverse_complement_matches = 0usize;
    let mut spelling_mismatches = 0usize;
    let mut first_mismatch: Option<String> = None;

    for (name, observed) in &paths {
        let Some((source_name, start, end)) = parse_path_interval(name) else {
            unparsable_paths += 1;
            if first_mismatch.is_none() {
                first_mismatch = Some(format!("path name is not source interval-like: {name}"));
            }
            continue;
        };
        let expected = sequence_index.fetch_sequence(source_name, start, end)?;
        let observed_bytes = observed.as_bytes();
        if expected.eq_ignore_ascii_case(observed_bytes) {
            forward_matches += 1;
            continue;
        }
        let expected_rc = impg::graph::reverse_complement(&expected);
        if expected_rc.eq_ignore_ascii_case(observed_bytes) {
            reverse_complement_matches += 1;
            continue;
        }

        spelling_mismatches += 1;
        if first_mismatch.is_none() {
            first_mismatch = Some(describe_mismatch(name, &expected, observed_bytes));
        }
    }

    println!("paths\t{}", paths.len());
    println!("forward_matches\t{forward_matches}");
    println!("reverse_complement_matches\t{reverse_complement_matches}");
    println!("unparsable_paths\t{unparsable_paths}");
    println!("spelling_mismatches\t{spelling_mismatches}");
    if let Some(first) = first_mismatch {
        println!("first_mismatch\t{first}");
    }

    if unparsable_paths == 0 && spelling_mismatches == 0 {
        Ok(())
    } else {
        Err(io::Error::other("GFA path sources differ"))
    }
}

fn usage() -> ! {
    eprintln!("usage: validate_gfa_path_sources <graph.gfa> <sequence.fa|sequence.agc>...");
    std::process::exit(2);
}

fn parse_path_interval(name: &str) -> Option<(&str, i32, i32)> {
    let (source_name, interval) = name.rsplit_once(':')?;
    let (start, end) = interval.split_once('-')?;
    let start = start.parse::<i32>().ok()?;
    let end = end.parse::<i32>().ok()?;
    Some((source_name, start.min(end), start.max(end)))
}

fn describe_mismatch(name: &str, expected: &[u8], observed: &[u8]) -> String {
    let first_diff = expected
        .iter()
        .zip(observed.iter())
        .position(|(a, b)| !a.eq_ignore_ascii_case(b))
        .unwrap_or_else(|| expected.len().min(observed.len()));
    format!(
        "{} at bp {} expected_len={} observed_len={} expected={} observed={}",
        name,
        first_diff,
        expected.len(),
        observed.len(),
        preview(expected, first_diff),
        preview(observed, first_diff)
    )
}

fn preview(seq: &[u8], center: usize) -> String {
    let start = center.saturating_sub(40);
    let end = seq.len().min(center.saturating_add(40));
    String::from_utf8_lossy(&seq[start..end]).into_owned()
}
