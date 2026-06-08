use std::collections::{BTreeMap, BTreeSet};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

const IMPG_BIN: &str = env!("CARGO_BIN_EXE_impg");
const GFAFFIX_BIN: &str = env!("CARGO_BIN_EXE_gfaffix");

fn impg_binary() -> PathBuf {
    let impg = PathBuf::from(IMPG_BIN);
    let gfaffix = Path::new(GFAFFIX_BIN);
    assert!(
        gfaffix.exists(),
        "Cargo did not build gfaffix at {GFAFFIX_BIN}"
    );

    let sibling = impg.with_file_name("gfaffix");
    if !sibling.exists() && sibling != gfaffix {
        std::fs::copy(gfaffix, &sibling).unwrap_or_else(|err| {
            panic!(
                "failed to prepare gfaffix sibling {} for {}: {err}",
                sibling.display(),
                impg.display()
            )
        });
    }
    impg
}

fn path_map(gfa: &str) -> BTreeMap<String, String> {
    impg::resolution::path_sequences(gfa)
        .expect("parse graph-output GFA paths")
        .into_iter()
        .collect()
}

fn varied_dna_sequence(len: usize, seed: u32) -> String {
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state.wrapping_mul(1664525).wrapping_add(1013904223);
            b"ACGT"[((state >> 24) & 3) as usize] as char
        })
        .collect()
}

fn assert_path_sequence_matches_source_slice(
    path_name: &str,
    path_sequence: &str,
    source_sequences: &BTreeMap<String, String>,
) {
    let (source_name, range) = path_name
        .rsplit_once(':')
        .unwrap_or_else(|| panic!("path name should include source coordinates: {path_name}"));
    let (start, end) = range
        .split_once('-')
        .unwrap_or_else(|| panic!("path name should include start-end range: {path_name}"));
    let start = start
        .parse::<usize>()
        .unwrap_or_else(|err| panic!("path start should parse in {path_name}: {err}"));
    let end = end
        .parse::<usize>()
        .unwrap_or_else(|err| panic!("path end should parse in {path_name}: {err}"));
    let source = source_sequences
        .get(source_name)
        .unwrap_or_else(|| panic!("unexpected source path name {source_name} in {path_name}"));
    assert!(
        start <= end && end <= source.len(),
        "path coordinates should be within source sequence bounds: {path_name}"
    );
    assert_eq!(
        path_sequence,
        &source[start..end],
        "path sequence should match the source slice named by {path_name}"
    );
}

#[test]
fn query_poa_crush_preserves_source_path_spellings() {
    let dir = tempfile::tempdir().expect("create temp dir");
    let fasta_path = dir.path().join("input.fa");
    let paf_path = dir.path().join("input.paf");
    let index_path = dir.path().join("input.impg");
    let out_prefix = dir.path().join("query.poa.crush");

    let seq_a = "ACGT".repeat(16);
    let mut seq_b = seq_a.clone().into_bytes();
    seq_b[32] = b'T';
    let seq_b = String::from_utf8(seq_b).unwrap();

    {
        let mut fasta = std::fs::File::create(&fasta_path).expect("create FASTA");
        writeln!(fasta, ">HG001#1#chr6").unwrap();
        writeln!(fasta, "{seq_a}").unwrap();
        writeln!(fasta, ">HG002#1#chr6").unwrap();
        writeln!(fasta, "{seq_b}").unwrap();
    }
    std::fs::write(
        &paf_path,
        "HG001#1#chr6\t64\t0\t64\t+\tHG002#1#chr6\t64\t0\t64\t63\t64\t60\tcg:Z:32=1X31=\n",
    )
    .expect("write PAF");

    let bin = impg_binary();
    let index = Command::new(&bin)
        .args([
            "index",
            "-a",
            paf_path.to_str().unwrap(),
            "-i",
            index_path.to_str().unwrap(),
        ])
        .output()
        .expect("run impg index");
    assert!(
        index.status.success(),
        "impg index failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&index.stdout),
        String::from_utf8_lossy(&index.stderr)
    );

    let query = Command::new(&bin)
        .args([
            "query",
            "-d",
            "0",
            "-i",
            index_path.to_str().unwrap(),
            "-a",
            paf_path.to_str().unwrap(),
            "-r",
            "HG001#1#chr6:0-64",
            "--min-transitive-len",
            "0",
            "-o",
            "gfa:poa:crush,method=poa,max-rounds=1",
            "--sequence-files",
            fasta_path.to_str().unwrap(),
            "-O",
            out_prefix.to_str().unwrap(),
            "-t",
            "1",
        ])
        .output()
        .expect("run impg query -o gfa:poa:crush");
    assert!(
        query.status.success(),
        "impg query -o gfa:poa:crush failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&query.stdout),
        String::from_utf8_lossy(&query.stderr)
    );

    let gfa_path = PathBuf::from(format!("{}.gfa", out_prefix.to_string_lossy()));
    let gfa = std::fs::read_to_string(&gfa_path).expect("read query GFA");
    let paths = path_map(&gfa);
    assert_eq!(paths["HG001#1#chr6:0-64"], seq_a);
    assert_eq!(paths["HG002#1#chr6:0-64"], seq_b);
    assert_eq!(
        paths.keys().cloned().collect::<Vec<_>>(),
        vec![
            "HG001#1#chr6:0-64".to_string(),
            "HG002#1#chr6:0-64".to_string()
        ],
        "graph-output crush should preserve source-coordinate path names without synthetic local IDs"
    );
}

#[test]
fn query_syng_local_seed_driver_preserves_collected_path_spellings() {
    let dir = tempfile::tempdir().expect("create temp dir");
    let fasta_path = dir.path().join("input.fa");
    let syng_prefix = dir.path().join("input.syng");
    let out_prefix = dir.path().join("query.syng-local.seed");

    let seq_a = varied_dna_sequence(512, 17);
    let mut seq_b = seq_a.clone().into_bytes();
    seq_b[128] = if seq_b[128] == b'A' { b'C' } else { b'A' };
    seq_b[320] = if seq_b[320] == b'G' { b'T' } else { b'G' };
    let seq_b = String::from_utf8(seq_b).unwrap();
    let source_sequences = BTreeMap::from([
        ("HG001#1#chr6".to_string(), seq_a.clone()),
        ("HG002#1#chr6".to_string(), seq_b.clone()),
    ]);

    {
        let mut fasta = std::fs::File::create(&fasta_path).expect("create FASTA");
        writeln!(fasta, ">HG001#1#chr6").unwrap();
        writeln!(fasta, "{seq_a}").unwrap();
        writeln!(fasta, ">HG002#1#chr6").unwrap();
        writeln!(fasta, "{seq_b}").unwrap();
    }

    let bin = impg_binary();
    let syng = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            syng_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("run impg syng");
    assert!(
        syng.status.success(),
        "impg syng failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&syng.stdout),
        String::from_utf8_lossy(&syng.stderr)
    );

    let query = Command::new(&bin)
        .args([
            "query",
            "-d",
            "0",
            "-a",
            syng_prefix.to_str().unwrap(),
            "--sequence-files",
            fasta_path.to_str().unwrap(),
            "-r",
            "HG001#1#chr6:0-512",
            "--min-transitive-len",
            "0",
            "-o",
            "gfa:fastga:syng-local:nosort",
            "--min-match-len",
            "1",
            "--no-filter",
            "-O",
            out_prefix.to_str().unwrap(),
            "-t",
            "1",
            "-v",
            "2",
        ])
        .env("RUST_LOG", "info")
        .output()
        .expect("run impg query -o gfa:fastga:syng-local");
    assert!(
        query.status.success(),
        "impg query -o gfa:fastga:syng-local failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&query.stdout),
        String::from_utf8_lossy(&query.stderr)
    );
    let stderr = String::from_utf8_lossy(&query.stderr);
    assert!(
        stderr.contains("[local seed] starting route=whole-region-sweepga-seqwish"),
        "query should route plain syng-local through explicit local seed induction:\n{stderr}"
    );
    assert!(
        stderr.contains("impg-local-seed route=whole-region-sweepga-seqwish aligner=fastga"),
        "local seed run should record a reproducible command/config:\n{stderr}"
    );

    let gfa_path = PathBuf::from(format!("{}.gfa", out_prefix.to_string_lossy()));
    let gfa = std::fs::read_to_string(&gfa_path).expect("read query GFA");
    let paths = path_map(&gfa);
    for (path_name, sequence) in &paths {
        assert!(
            !path_name.starts_with("local_") && !path_name.starts_with("__impg"),
            "syng-local seed graph should not introduce synthetic local IDs: {path_name}"
        );
        assert_path_sequence_matches_source_slice(path_name, sequence, &source_sequences);
    }
    let sources = paths
        .keys()
        .map(|name| {
            name.rsplit_once(':')
                .map(|(source, _)| source)
                .unwrap_or(name)
        })
        .collect::<BTreeSet<_>>();
    assert_eq!(
        sources,
        BTreeSet::from(["HG001#1#chr6", "HG002#1#chr6"]),
        "SYNG-collected fixture should contribute both PanSN haplotype source names"
    );
    assert_eq!(
        paths.len(),
        2,
        "syng-local seed graph should preserve one named local path per collected haplotype"
    );
}
