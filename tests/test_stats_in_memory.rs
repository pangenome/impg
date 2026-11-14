//! `impg stats` must work on an in-memory index (no `-i`).
//!
//! `stats` walks every interval tree twice: once for overlap counts, once for
//! the bridge-genome coverage. It evicts each tree after use to keep memory
//! flat. An in-memory index has no index file, so the second walk used to call
//! `File::open("")` and panic in a worker thread:
//!
//! ```text
//! thread '<unnamed>' panicked at src/impg.rs: Failed to open index file:
//! Os { code: 2, kind: NotFound, message: "No such file or directory" }
//! ```

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

fn get_impg_binary() -> PathBuf {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return PathBuf::from(path);
    }
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    for candidate in [
        manifest_dir.join("target/release/impg"),
        manifest_dir.join("target/debug/impg"),
    ] {
        if candidate.exists() {
            return candidate;
        }
    }
    PathBuf::from("impg")
}

fn write_paf(path: &std::path::Path) {
    // Two records over three sequences: q -> t and t -> r.
    // Bidirectional indexing turns them into four overlaps.
    let mut file = File::create(path).unwrap();
    writeln!(
        file,
        "q\t1000\t0\t500\t+\tt\t1000\t0\t500\t500\t500\t60\tcg:Z:500="
    )
    .unwrap();
    writeln!(
        file,
        "t\t1000\t100\t600\t+\tr\t1000\t100\t600\t500\t500\t60\tcg:Z:500="
    )
    .unwrap();
}

#[test]
fn test_stats_on_in_memory_index_does_not_panic() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();

    write_paf(&work_dir.join("a.paf"));

    // No -i: build the index in memory, which is the case that used to panic.
    let output = Command::new(get_impg_binary())
        .current_dir(work_dir)
        .args(["stats", "-a", "a.paf"])
        .output()
        .expect("failed to run impg stats");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    assert!(
        !stderr.contains("panicked"),
        "impg stats panicked:\n{stderr}\n\nstdout:\n{stdout}"
    );
    assert!(
        output.status.success(),
        "impg stats failed with {:?}:\n{stderr}",
        output.status
    );

    // The overlap count comes from the first walk and must be complete.
    assert!(
        stdout.contains("Number of overlaps: 4"),
        "expected 4 overlaps (2 records, bidirectional), got:\n{stdout}"
    );

    // The bridge coverage comes from the second walk. Before the fix it was
    // computed from zero trees, so it must now see all three sequences.
    assert!(
        stdout.contains("Bridge sequences (both): 3"),
        "expected 3 bridge sequences, got:\n{stdout}"
    );
}
