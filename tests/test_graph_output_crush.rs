use std::collections::BTreeMap;
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
