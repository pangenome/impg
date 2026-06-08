use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::time::{Duration, Instant};

const IMPG_BIN: &str = env!("CARGO_BIN_EXE_impg");
const GFAFFIX_BIN: &str = env!("CARGO_BIN_EXE_gfaffix");

const C4_A_NAME: &str = "C4SMOKE_A#0#chr6";
const C4_B_NAME: &str = "C4SMOKE_B#0#chr6";
const C4_A_SEQ: &str = "CCTCGGTCTCGGTGTTTGTGGACCATCACCTGGCACCCTCCTTCTCTTTGTGGCCTTCTACTACCATGGAGACCACCAGTGGCCAACTCCCTGCGAGTGGATGTCCAGGCTGGGCCTGCGAGGGCAAGGTGACCGGGGTCAGGAGAGTGGCACTTGTGCCGAGGGGGTTGAGACAGGGTGATTGCCAACAGGGCTGGATTTAGCTTGGGGGCAGTGAGGATACCGG";
const C4_B_SEQ: &str = "CCTCGGTCTCGGTGTTTGTGGACCATCACCTGGCACCCTCCTTCTCTTTGTGGCCTTCTACTACCCTGGAGACCACCAGTGGCCAACTCCCTGCGAGTGGATGTCCAGGCTGGGCCTGCGAGGGCAAGGTGACCGGGGTCAGGAGAGTGGCACTTGTGCCGAGGGGGTTGAGACAGGGTGATTGCCAACAGGGCTGGATTTAGCTTGGGGGCAGTGAGGATACCGG";
const C4_ONE_ROW_PAF: &str = "C4SMOKE_A#0#chr6\t226\t0\t226\t+\tC4SMOKE_B#0#chr6\t226\t0\t226\t225\t226\t24\tgi:f:0.995575\tbi:f:0.995575\tmd:f:0.992962\tcg:Z:65=1X160=\tch:Z:chain_41\tst:Z:scaffold\n";

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

fn run_with_timeout(command: &mut Command, timeout: Duration) -> (Output, Duration) {
    command.stdout(Stdio::piped()).stderr(Stdio::piped());
    let mut child = command.spawn().expect("spawn impg graph");
    let start = Instant::now();

    loop {
        if let Some(_status) = child.try_wait().expect("poll impg graph") {
            let elapsed = start.elapsed();
            let output = child.wait_with_output().expect("collect impg graph output");
            return (output, elapsed);
        }

        if start.elapsed() >= timeout {
            let _ = child.kill();
            let output = child
                .wait_with_output()
                .expect("collect timed-out impg graph output");
            panic!(
                "impg graph exceeded {:?}\nstdout:\n{}\nstderr:\n{}",
                timeout,
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
        }

        std::thread::sleep(Duration::from_millis(10));
    }
}

fn run_c4_seqwish_graph(engine: &str, paf_contents: &str) -> (String, Duration) {
    let dir = tempfile::tempdir().expect("create test temp dir");
    let fasta_path = dir.path().join("c4.fa");
    let paf_path = dir.path().join("input.paf");
    let gfa_path = dir.path().join("output.gfa");
    let temp_path = dir.path().join("tmp");
    std::fs::create_dir(&temp_path).expect("create graph temp dir");

    {
        let mut fasta = std::fs::File::create(&fasta_path).expect("create FASTA");
        writeln!(fasta, ">{C4_A_NAME}").unwrap();
        writeln!(fasta, "{C4_A_SEQ}").unwrap();
        writeln!(fasta, ">{C4_B_NAME}").unwrap();
        writeln!(fasta, "{C4_B_SEQ}").unwrap();
    }
    std::fs::write(&paf_path, paf_contents).expect("write PAF");

    let mut command = Command::new(impg_binary());
    command
        .arg("graph")
        .arg("--sequence-files")
        .arg(&fasta_path)
        .arg("--paf-file")
        .arg(&paf_path)
        .arg("--gfa-engine")
        .arg(engine)
        .arg("--no-filter")
        .arg("--temp-dir")
        .arg(&temp_path)
        .arg("-g")
        .arg(&gfa_path)
        .arg("-t")
        .arg("1");

    let timeout = Duration::from_secs(10);
    let (output, elapsed) = run_with_timeout(&mut command, timeout);
    assert!(
        output.status.success(),
        "impg graph failed for engine {engine}\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        elapsed < timeout,
        "impg graph for engine {engine} took {:?}, expected under {:?}",
        elapsed,
        timeout
    );

    (
        std::fs::read_to_string(&gfa_path).expect("read output GFA"),
        elapsed,
    )
}

fn assert_valid_c4_gfa(gfa: &str) {
    assert!(
        gfa.lines().any(|line| line == "H\tVN:Z:1.0"),
        "missing GFA header:\n{gfa}"
    );
    assert!(
        gfa.lines().any(|line| line.starts_with("S\t")),
        "missing GFA segment:\n{gfa}"
    );
    assert!(
        gfa.lines().any(|line| line.starts_with("P\t")),
        "missing GFA path:\n{gfa}"
    );

    let paths = impg::resolution::path_sequences(gfa)
        .expect("parse GFA path spellings")
        .into_iter()
        .collect::<BTreeMap<_, _>>();
    assert_eq!(paths.len(), 2, "unexpected path set in GFA:\n{gfa}");
    assert_eq!(paths[C4_A_NAME], C4_A_SEQ);
    assert_eq!(paths[C4_B_NAME], C4_B_SEQ);
}

fn assert_segment_ids_are_sorted(gfa: &str) {
    let segment_ids = gfa
        .lines()
        .filter_map(|line| {
            let mut fields = line.split('\t');
            match fields.next() {
                Some("S") => fields.next()?.parse::<u64>().ok(),
                _ => None,
            }
        })
        .collect::<Vec<_>>();
    let mut sorted = segment_ids.clone();
    sorted.sort_unstable();
    assert_eq!(segment_ids, sorted, "segment IDs are not sorted:\n{gfa}");
}

#[test]
fn graph_seqwish_cli_c4_tiny_one_row_paf_finishes_under_10s() {
    let (gfa, _elapsed) = run_c4_seqwish_graph("seqwish", C4_ONE_ROW_PAF);
    assert_valid_c4_gfa(&gfa);
}

#[test]
fn graph_seqwish_crush_cli_c4_tiny_empty_paf_finishes_under_10s_and_sorts_gfa() {
    let (gfa, _elapsed) =
        run_c4_seqwish_graph("seqwish:crush,method=poa,max-rounds=1,polish-rounds=0", "");
    assert_valid_c4_gfa(&gfa);
    assert_segment_ids_are_sorted(&gfa);
}
