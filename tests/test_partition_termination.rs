//! Regression tests for `impg partition` termination.
//!
//! `partition_alignments` re-selects the longest missing region every round.
//! If a selected region produces no query overlaps, nothing masks it, the
//! missing map never changes, and the same windows are re-queried forever —
//! a silent single-core infinite loop (see docs/bug-infer-partition-hang.md).
//!
//! Covered here:
//! 1. The stall guard force-masks a re-selected zero-progress region so the
//!    loop always terminates.
//! 2. syng transitive queries tolerate anchorless paths (sequences shorter
//!    than the syncmer length w+k) instead of erroring out; the walk returns
//!    empty, the guard force-masks, and partitioning finishes.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::time::{Duration, Instant};

fn get_impg_binary() -> PathBuf {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return PathBuf::from(path);
    }
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let candidates = [
        manifest_dir.join("target/release/impg"),
        manifest_dir.join("target/debug/impg"),
    ];
    for path in &candidates {
        if path.exists() {
            return path.clone();
        }
    }
    PathBuf::from("impg")
}

/// Run impg with a hard wall clock. A pre-fix infinite loop shows up as a
/// killed, non-zero-exit run instead of a hung test suite.
fn run_impg_bounded(work_dir: &Path, args: &[&str], wall: Duration) -> (bool, String) {
    let impg = get_impg_binary();
    let mut child: Child = Command::new(&impg)
        .current_dir(work_dir)
        .args(args)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn impg");
    let deadline = Instant::now() + wall;
    let completed = loop {
        match child.try_wait().expect("try_wait failed") {
            Some(status) => break status.success(),
            None => {
                if Instant::now() >= deadline {
                    child.kill().expect("failed to kill hung impg");
                    let _ = child.wait();
                    return (false, String::from("KILLED: exceeded wall clock"));
                }
                std::thread::sleep(Duration::from_millis(100));
            }
        }
    };
    let mut stderr = String::new();
    if let Some(mut pipe) = child.stderr.take() {
        let _ = pipe.read_to_string(&mut stderr);
    }
    (completed, stderr)
}

fn write_fasta(path: &Path, records: &[(&str, String)]) {
    use std::io::Write;
    let mut f = File::create(path).unwrap();
    for (name, seq) in records {
        writeln!(f, ">{name}").unwrap();
        f.write_all(seq.as_bytes()).unwrap();
        writeln!(f).unwrap();
    }
}

use std::fs::File;

fn random_acgt(len: usize, seed: u64) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut state = seed | 1;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            bases[((state >> 33) % 4) as usize] as char
        })
        .collect()
}

#[test]
fn test_partition_syng_terminates_on_anchorless_sequence() {
    let temp_dir = tempfile::TempDir::new().unwrap();
    let work_dir = temp_dir.path();

    let backbone = random_acgt(200, 13);
    let seq_b = {
        let mut s = backbone.clone();
        s.replace_range(180..200, &random_acgt(20, 5));
        s
    };
    write_fasta(
        &work_dir.join("tiny.fa"),
        &[
            ("sT", "A".repeat(50)), // shorter than w+k: no syncmers, no homologs
            ("sA", backbone),
            ("sB", seq_b),
        ],
    );

    let build = run_impg_bounded(
        work_dir,
        &["syng", "-f", "tiny.fa", "-o", "tidx", "-t", "2"],
        Duration::from_secs(60),
    );
    assert!(build.0, "syng build failed: {}", build.1);

    let (ok, stderr) = run_impg_bounded(
        work_dir,
        &[
            "partition",
            "-a",
            "tidx",
            "-w",
            "1000",
            "-d",
            "100",
            "-o",
            "bed",
            "-t",
            "2",
        ],
        Duration::from_secs(120),
    );
    assert!(ok, "partition did not terminate cleanly: {}", stderr);
    assert!(
        stderr.contains("force-masking"),
        "expected stall guard warning, stderr: {stderr}"
    );
}
