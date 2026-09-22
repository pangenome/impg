use std::collections::HashMap;
use std::env;
use std::fs;
use std::io;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let Some(expected_path) = args.next() else {
        eprintln!("usage: compare_gfa_paths <expected.gfa> <observed.gfa>");
        std::process::exit(2);
    };
    let Some(observed_path) = args.next() else {
        eprintln!("usage: compare_gfa_paths <expected.gfa> <observed.gfa>");
        std::process::exit(2);
    };
    if args.next().is_some() {
        eprintln!("usage: compare_gfa_paths <expected.gfa> <observed.gfa>");
        std::process::exit(2);
    }

    let expected_gfa = fs::read_to_string(&expected_path)?;
    let observed_gfa = fs::read_to_string(&observed_path)?;
    let expected = path_map(&expected_gfa)?;
    let observed = path_map(&observed_gfa)?;

    let mut missing = 0usize;
    let mut extra = 0usize;
    let mut spelling_mismatches = 0usize;
    for (name, expected_seq) in &expected {
        match observed.get(name) {
            Some(observed_seq) if observed_seq == expected_seq => {}
            Some(_) => spelling_mismatches += 1,
            None => missing += 1,
        }
    }
    for name in observed.keys() {
        if !expected.contains_key(name) {
            extra += 1;
        }
    }

    println!("expected_paths\t{}", expected.len());
    println!("observed_paths\t{}", observed.len());
    println!("missing_paths\t{missing}");
    println!("extra_paths\t{extra}");
    println!("spelling_mismatches\t{spelling_mismatches}");

    if missing == 0 && extra == 0 && spelling_mismatches == 0 {
        Ok(())
    } else {
        Err(io::Error::other("GFA path spellings differ"))
    }
}

fn path_map(gfa: &str) -> io::Result<HashMap<String, String>> {
    impg::resolution::path_sequences(gfa).map(|paths| paths.into_iter().collect())
}
