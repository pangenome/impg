use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::Path;
use std::process::Command;

const IMPG_BIN: &str = env!("CARGO_BIN_EXE_impg");
const GFAFFIX_BIN: &str = env!("CARGO_BIN_EXE_gfaffix");

fn run_graph_poa(records: &[(&str, &str)]) -> String {
    assert!(
        Path::new(GFAFFIX_BIN).exists(),
        "Cargo did not build gfaffix at {GFAFFIX_BIN}"
    );

    let dir = tempfile::tempdir().expect("create test temp dir");
    let fasta_path = dir.path().join("input.fa");
    let gfa_path = dir.path().join("output.gfa");

    {
        let mut fasta = std::fs::File::create(&fasta_path).expect("create FASTA input");
        for (name, seq) in records {
            writeln!(fasta, ">{name}").expect("write FASTA header");
            writeln!(fasta, "{seq}").expect("write FASTA sequence");
        }
    }

    let output = Command::new(IMPG_BIN)
        .args([
            "graph",
            "--sequence-files",
            fasta_path.to_str().expect("FASTA path is UTF-8"),
            "--gfa-engine",
            "poa",
            "-g",
            gfa_path.to_str().expect("GFA path is UTF-8"),
        ])
        .output()
        .expect("run impg graph --gfa-engine poa");

    assert!(
        output.status.success(),
        "impg graph --gfa-engine poa failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    std::fs::read_to_string(&gfa_path).expect("read POA GFA output")
}

fn run_graph_pggb_empty_paf(records: &[(&str, &str)]) -> String {
    assert!(
        Path::new(GFAFFIX_BIN).exists(),
        "Cargo did not build gfaffix at {GFAFFIX_BIN}"
    );

    let dir = tempfile::tempdir().expect("create test temp dir");
    let fasta_path = dir.path().join("input.fa");
    let paf_path = dir.path().join("empty.paf");
    let gfa_path = dir.path().join("output.gfa");

    {
        let mut fasta = std::fs::File::create(&fasta_path).expect("create FASTA input");
        for (name, seq) in records {
            writeln!(fasta, ">{name}").expect("write FASTA header");
            writeln!(fasta, "{seq}").expect("write FASTA sequence");
        }
    }
    std::fs::write(&paf_path, "").expect("write empty PAF");

    let output = Command::new(IMPG_BIN)
        .args([
            "graph",
            "--sequence-files",
            fasta_path.to_str().expect("FASTA path is UTF-8"),
            "--paf-file",
            paf_path.to_str().expect("PAF path is UTF-8"),
            "--gfa-engine",
            "pggb",
            "--threads",
            "1",
            "--min-match-len",
            "1",
            "--target-poa-length",
            "64",
            "--max-node-length",
            "64",
            "--poa-padding-fraction",
            "0",
            "-g",
            gfa_path.to_str().expect("GFA path is UTF-8"),
        ])
        .output()
        .expect("run impg graph --gfa-engine pggb");

    assert!(
        output.status.success(),
        "impg graph --gfa-engine pggb failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    std::fs::read_to_string(&gfa_path).expect("read PGGB GFA output")
}

fn assert_gfa_paths_match_records(gfa: &str, records: &[(&str, &str)]) {
    let observed = impg::resolution::path_sequences(gfa)
        .expect("parse output GFA path spellings")
        .into_iter()
        .collect::<BTreeMap<_, _>>();
    let expected = records
        .iter()
        .map(|(name, seq)| ((*name).to_string(), (*seq).to_string()))
        .collect::<BTreeMap<_, _>>();

    assert_eq!(observed, expected);
}

#[test]
fn graph_poa_cli_preserves_pansn_coordinate_headers_and_spellings() {
    let records = [
        ("HG001#1#chr6:100-220(+)", "ACGTAAAACCCCGGGGTTTT"),
        ("HG002#2#chr6:101-221(-)", "ACGTAAAAGCCCGGGGTTTT"),
        ("HG003#1#chr6:102-222(+)", "ACGTAAAACCCCGGAGTTTT"),
    ];

    let gfa = run_graph_poa(&records);
    assert_gfa_paths_match_records(&gfa, &records);
}

#[test]
fn graph_poa_cli_preserves_homopolymer_loop_like_headers_and_spellings() {
    let records = [
        (
            "HG001#1#chr6:7068-7183(+)",
            "TTGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCGTACGTAGG",
        ),
        (
            "HG002#1#chr6:7068-7184(+)",
            "TTGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTACGTAGG",
        ),
        (
            "HG003#2#chr6:7069-7182(-)",
            "TTGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCGTACGTAGG",
        ),
    ];

    let gfa = run_graph_poa(&records);
    assert_gfa_paths_match_records(&gfa, &records);
}

#[test]
fn graph_pggb_cli_preserves_direct_fasta_coordinate_headers_after_smoothing() {
    let records = [
        ("ref#0#snp_bubble_3path:1-27", "ACGTACGTACGTCGATTACAGATTACA"),
        (
            "hapA#0#snp_bubble_3path:1-27",
            "ACGTACGTACGTTGATTACAGATTACA",
        ),
        (
            "hapB#0#snp_bubble_3path:1-27",
            "ACGTACGTACGTGGATTACAGATTACA",
        ),
    ];

    let gfa = run_graph_pggb_empty_paf(&records);
    assert_gfa_paths_match_records(&gfa, &records);
}
