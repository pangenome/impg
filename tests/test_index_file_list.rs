//! An index records the ordered list of alignment files it was built from.
//!
//! The trees store a file index and a byte offset for every CIGAR. Pass the same
//! files in a different order and every file index still resolves, so the query
//! reads CIGARs from the wrong file and returns wrong intervals with no error.
//! The index header now records the list, and loading refuses a mismatch.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
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

fn run(binary: &Path, dir: &Path, args: &[&str]) -> Output {
    Command::new(binary)
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to run {:?}: {e}", args))
}

/// Two alignment files with different content, so their sizes differ and a
/// swap is detectable.
fn write_fixtures(dir: &Path) {
    let mut one = File::create(dir.join("one.paf")).unwrap();
    writeln!(
        one,
        "q1\t2000\t0\t1000\t+\tt1\t2000\t0\t1000\t1000\t1000\t60\tcg:Z:1000="
    )
    .unwrap();

    let mut two = File::create(dir.join("two.paf")).unwrap();
    writeln!(
        two,
        "q2\t3000\t0\t1000\t+\tt2\t3000\t500\t1500\t1000\t1000\t60\tcg:Z:1000="
    )
    .unwrap();
    // A second record, so the two files do not have the same size.
    writeln!(
        two,
        "q2\t3000\t1000\t2000\t+\tt3\t3000\t0\t1000\t1000\t1000\t60\tcg:Z:1000="
    )
    .unwrap();
}

#[test]
fn test_index_refuses_a_reordered_alignment_file_list() {
    let temp_dir = TempDir::new().unwrap();
    let work_dir = temp_dir.path();
    let binary = get_impg_binary();
    write_fixtures(work_dir);

    let build = run(
        &binary,
        work_dir,
        &[
            "index",
            "-a",
            "one.paf",
            "two.paf",
            "-i",
            "idx.impg",
        ],
    );
    assert!(
        build.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    // The same list in the same order still works.
    let same = run(
        &binary,
        work_dir,
        &[
            "query", "-a", "one.paf", "two.paf", "-i", "idx.impg", "-r", "t1:0-1000", "-d", "0",
            "-o", "bed",
        ],
    );
    assert!(
        same.status.success(),
        "the original list must still load: {}",
        String::from_utf8_lossy(&same.stderr)
    );

    // A reordered list must be refused, not silently read wrong.
    let reordered = run(
        &binary,
        work_dir,
        &[
            "query", "-a", "two.paf", "one.paf", "-i", "idx.impg", "-r", "t1:0-1000", "-d", "0",
            "-o", "bed",
        ],
    );
    assert!(
        !reordered.status.success(),
        "a reordered file list must be refused"
    );
    let stderr = String::from_utf8_lossy(&reordered.stderr);
    assert!(
        stderr.contains("different alignment file list"),
        "the error must say the list changed: {stderr}"
    );

    // A subset must be refused as well: the file indices no longer line up.
    let subset = run(
        &binary,
        work_dir,
        &[
            "query", "-a", "two.paf", "-i", "idx.impg", "-r", "t1:0-1000", "-d", "0", "-o", "bed",
        ],
    );
    assert!(!subset.status.success(), "a shorter file list must be refused");
    assert!(
        String::from_utf8_lossy(&subset.stderr).contains("was built from 2"),
        "the error must give the recorded count"
    );
}
