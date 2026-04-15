//! End-to-end integration tests for the syng CLI path.
//!
//! Covers gaps that unit tests missed:
//! - `impg syng --agc` builds a non-empty index (caught a silent failure in the
//!   yeast235 workflow where only `.syng.names` got written).
//! - `impg syng --agc` round-trips: built index can be loaded and queried.
//! - `impg partition --syng` runs end-to-end and produces non-empty BED.

use ahash::AHashSet;
use ragc_core::{StreamingQueueCompressor, StreamingQueueConfig};
use std::path::PathBuf;
use std::process::Command;

// All syng tests must be serialized because the C library has non-thread-safe
// global state (see src/syng.rs test comments). tests/ and src/ test modules
// run in separate processes, but within this binary we still need a lock.
static SYNG_LOCK: std::sync::LazyLock<std::sync::Mutex<()>> =
    std::sync::LazyLock::new(|| std::sync::Mutex::new(()));

fn lock_syng() -> std::sync::MutexGuard<'static, ()> {
    SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

/// Deterministic pseudo-random sequence (numeric encoding, 0-3 = A/C/G/T) —
/// what ragc-core's push() expects.
fn make_sequence_numeric(len: usize, seed: u8) -> Vec<u8> {
    let mut seq = Vec::with_capacity(len);
    let mut state: u32 = seed as u32;
    for _ in 0..len {
        state = state.wrapping_mul(1103515245).wrapping_add(12345);
        seq.push(((state >> 16) % 4) as u8);
    }
    seq
}

/// Same sequence but as ASCII for FASTA output.
fn numeric_to_ascii(numeric: &[u8]) -> Vec<u8> {
    numeric
        .iter()
        .map(|&b| match b {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        })
        .collect()
}

/// Build a small test AGC archive with multiple samples sharing a backbone.
fn create_test_agc(path: &str) {
    let config = StreamingQueueConfig {
        queue_capacity: 10 * 1024 * 1024,
        num_threads: 1,
        verbosity: 0,
        ..Default::default()
    };

    // Shared backbone (same seed) so samples share syncmer nodes.
    let backbone = make_sequence_numeric(800, 42);

    let samples: Vec<(&str, &str, Vec<u8>)> = vec![
        ("sampleA", "chr1", {
            let mut s = backbone.clone();
            s.extend(make_sequence_numeric(400, 1));
            s
        }),
        ("sampleB", "chr1", {
            let mut s = backbone.clone();
            s.extend(make_sequence_numeric(400, 2));
            s
        }),
        ("sampleC", "chr1", {
            // Completely different — should NOT share nodes
            make_sequence_numeric(1200, 99)
        }),
    ];

    let mut compressor =
        StreamingQueueCompressor::with_splitters(path, config, AHashSet::new())
            .expect("Failed to create AGC compressor");

    for (sample, contig, data) in samples {
        compressor
            .push(sample.to_string(), contig.to_string(), data)
            .expect("Failed to push contig");
    }
    compressor.finalize().expect("Failed to finalize AGC");
}

/// Locate the compiled `impg` binary for CLI tests.
fn impg_binary() -> Option<PathBuf> {
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    for candidate in [
        manifest_dir.join("target/debug/impg"),
        manifest_dir.join("target/release/impg"),
    ] {
        if candidate.exists() {
            return Some(candidate);
        }
    }
    None
}

#[test]
fn test_syng_agc_build_produces_non_empty_index() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_agc_build");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let agc_path = dir.join("test.agc");
    create_test_agc(agc_path.to_str().unwrap());

    let out_prefix = dir.join("idx");

    let output = Command::new(&bin)
        .args([
            "syng",
            "--agc",
            agc_path.to_str().unwrap(),
            "-o",
            out_prefix.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run impg syng --agc");

    assert!(
        output.status.success(),
        "impg syng --agc failed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let gbwt_path = format!("{}.1gbwt", out_prefix.to_str().unwrap());
    let khash_path = format!("{}.1khash", out_prefix.to_str().unwrap());
    let names_path = format!("{}.syng.names", out_prefix.to_str().unwrap());

    assert!(std::path::Path::new(&gbwt_path).exists(), ".1gbwt missing");
    assert!(std::path::Path::new(&khash_path).exists(), ".1khash missing");
    assert!(std::path::Path::new(&names_path).exists(), ".syng.names missing");

    // Key assertion: the GBWT must contain actual vertex data, not just a
    // schema header. The yeast235 bug slipped through because nothing
    // asserted on file size. An empty GBWT is ~1300 bytes; a populated one
    // for even this tiny test should be several KB.
    let gbwt_size = std::fs::metadata(&gbwt_path).unwrap().len();
    assert!(
        gbwt_size > 2000,
        ".1gbwt is only {} bytes — suggests no syncmer vertices were written. \
         This is the bug that slipped through on yeast235.agc.",
        gbwt_size
    );

    let khash_size = std::fs::metadata(&khash_path).unwrap().len();
    assert!(
        khash_size > 1000,
        ".1khash is only {} bytes — suggests no k-mer hash entries",
        khash_size
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_agc_roundtrip_query() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_agc_query");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let agc_path = dir.join("test.agc");
    create_test_agc(agc_path.to_str().unwrap());

    let out_prefix = dir.join("idx");
    let out_prefix_str = out_prefix.to_str().unwrap();

    let build_output = Command::new(&bin)
        .args(["syng", "--agc", agc_path.to_str().unwrap(), "-o", out_prefix_str])
        .output()
        .expect("Failed to run impg syng --agc");
    assert!(
        build_output.status.success(),
        "impg syng --agc failed. stderr: {}",
        String::from_utf8_lossy(&build_output.stderr)
    );

    // Load the index back and run a query on the shared backbone region.
    let index = impg::syng::SyngIndex::load(
        out_prefix_str,
        impg::syng::SyncmerParams::default(),
    )
    .expect("Failed to load built AGC-derived syng index");

    assert_eq!(
        index.name_map.path_to_name.len(),
        3,
        "Expected 3 sequences in name map"
    );

    // Find one of the contigs that contains the shared backbone.
    let query_name = index
        .name_map
        .path_to_name
        .iter()
        .find(|n| n.contains("sampleA"))
        .expect("sampleA contig should be in name map")
        .clone();

    let intervals = index
        .query_region(&query_name, 0, 500, 0)
        .expect("query_region failed");

    // Should find at least self-hit + sampleB (shared backbone)
    assert!(
        intervals.len() >= 2,
        "Expected at least 2 homologous intervals (self + sampleB), got {}: {:?}",
        intervals.len(),
        intervals.iter().map(|iv| &iv.genome).collect::<Vec<_>>()
    );

    let genomes: Vec<&str> = intervals.iter().map(|iv| iv.genome.as_str()).collect();
    assert!(
        genomes.iter().any(|g| g.contains("sampleA")),
        "Should find self-hit, got: {:?}",
        genomes
    );
    assert!(
        genomes.iter().any(|g| g.contains("sampleB")),
        "Should find sampleB (shared backbone), got: {:?}",
        genomes
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_fasta_build_produces_non_empty_index() {
    // Same non-empty-GBWT assertion but via the FASTA input path, to catch
    // regressions in either build path symmetrically.
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_fasta_nonempty");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        let backbone = numeric_to_ascii(&make_sequence_numeric(800, 42));
        let tail_a = numeric_to_ascii(&make_sequence_numeric(400, 1));
        let tail_b = numeric_to_ascii(&make_sequence_numeric(400, 2));
        let other = numeric_to_ascii(&make_sequence_numeric(1200, 99));

        writeln!(f, ">sampleA#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleC#chr1").unwrap();
        f.write_all(&other).unwrap();
        writeln!(f).unwrap();
    }

    let out_prefix = dir.join("idx");
    let out = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            out_prefix.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run impg syng -f");
    assert!(
        out.status.success(),
        "impg syng -f failed. stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let gbwt_size =
        std::fs::metadata(format!("{}.1gbwt", out_prefix.to_str().unwrap()))
            .unwrap()
            .len();
    assert!(
        gbwt_size > 2000,
        ".1gbwt is only {} bytes from FASTA build",
        gbwt_size
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_identical_sequences_build_and_query() {
    // Regression test for the vendored syng hash.c REMOVED-sentinel bug: two
    // byte-identical sequences (like AAA#0#chrIII and SGDref#0#chrIII in
    // yeast235) would both get stored as (start_node=1, start_count=0),
    // colliding in the GBWT startCount table. On query, the C layer would
    // die with "syngBWTpathStartOld startNode 1 count 0 >= startCount 0".
    //
    // This test builds two identical sequences via the CLI, loads the index,
    // and queries — it crashes at the C layer if the hash bug comes back.
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_ident_seqs_cli");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        // Two records, byte-identical sequence content
        let seq = numeric_to_ascii(&make_sequence_numeric(1500, 42));
        writeln!(f, ">sampleA#0#chrIII").unwrap();
        f.write_all(&seq).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chrIII").unwrap();
        f.write_all(&seq).unwrap();
        writeln!(f).unwrap();
    }

    let out_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args(["syng", "-f", fasta_path.to_str().unwrap(), "-o", out_prefix.to_str().unwrap()])
        .output()
        .expect("impg syng failed to run");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    // Load and query — this is where the bug used to crash at C level.
    let index = impg::syng::SyngIndex::load(
        out_prefix.to_str().unwrap(),
        impg::syng::SyncmerParams::default(),
    )
    .expect("load failed");

    // Verify start_counts are distinct (the bug used to save both as 0)
    let start_a = index.name_map.path_starts[0].as_ref().expect("seqA missing");
    let start_b = index.name_map.path_starts[1].as_ref().expect("seqB missing");
    assert_eq!(
        start_a.start_node, start_b.start_node,
        "identical sequences should share start_node"
    );
    assert_ne!(
        start_a.start_count, start_b.start_count,
        "two paths at same start_node must have distinct start_count indices \
         (vendored syng hash.c REMOVED-sentinel bug regression)"
    );

    // Query both sequences — used to crash with startCount 0 >= 0
    let intervals_a = index
        .query_region("sampleA#0#chrIII", 0, 1000, 0)
        .expect("query for sampleA should succeed");
    let intervals_b = index
        .query_region("sampleB#0#chrIII", 0, 1000, 0)
        .expect("query for sampleB should succeed");

    assert!(!intervals_a.is_empty());
    assert!(!intervals_b.is_empty());

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_partition_syng_end_to_end_bed() {
    // End-to-end: build syng index from FASTA, then run `impg partition --syng`
    // and verify non-empty BED output. This is the path I added but never
    // exercised via CLI until now.
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_partition_syng_bed");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    // Build a FASTA with enough content to produce multiple partition windows.
    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        // Two sequences sharing a 2kbp backbone so partition finds homology
        let backbone = numeric_to_ascii(&make_sequence_numeric(2000, 42));
        let tail1 = numeric_to_ascii(&make_sequence_numeric(1000, 1));
        let tail2 = numeric_to_ascii(&make_sequence_numeric(1000, 2));

        writeln!(f, ">sampleA#1#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail1).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#1#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail2).unwrap();
        writeln!(f).unwrap();
    }

    // Step 1: build syng index
    let syng_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            syng_prefix.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed. stderr: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    // Step 2: run partition --syng with BED output
    let out_folder = dir.join("parts");
    std::fs::create_dir_all(&out_folder).unwrap();
    let part = Command::new(&bin)
        .args([
            "partition",
            "--syng",
            syng_prefix.to_str().unwrap(),
            "-w",
            "1500",
            "-o",
            "bed",
            "--output-folder",
            out_folder.to_str().unwrap(),
            "--min-missing-size",
            "100",
            "--min-boundary-distance",
            "0",
            "-t",
            "1",
        ])
        .output()
        .expect("Failed to run impg partition --syng");
    assert!(
        part.status.success(),
        "impg partition --syng failed. stdout: {} stderr: {}",
        String::from_utf8_lossy(&part.stdout),
        String::from_utf8_lossy(&part.stderr)
    );

    // Single-file BED mode writes to `partitions.bed` in the output folder.
    let bed_files: Vec<_> = std::fs::read_dir(&out_folder)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("bed"))
        .collect();

    assert!(
        !bed_files.is_empty(),
        "No BED files produced in {:?}",
        out_folder
    );

    // At least one BED should have content (non-empty partition)
    let any_non_empty = bed_files
        .iter()
        .any(|e| std::fs::metadata(e.path()).unwrap().len() > 0);
    assert!(any_non_empty, "All partition BED files are empty");

    std::fs::remove_dir_all(&dir).ok();
}

/// Regression test: `impg query --syng -o gfa --gfa-engine seqwish:X` must
/// treat `X` as a sub-window size and split the query range into per-window
/// partitions, not silently collapse into a flat single-engine run.
///
/// Before this fix, `pggb:X`/`seqwish:X` in the syng+gfa path set
/// `partition_size = Some(X)` but was only used as a `skip_normalize` boolean
/// flag — the actual `X` was discarded and one flat engine call ran on the
/// entire query range. This broke for large regions because syng's
/// `query_region` returned whole-chromosome context spans.
///
/// We verify the new behavior by checking that the per-sub-window log lines
/// show up in stderr, one per sub-window.
#[test]
fn test_query_syng_gfa_subwindow_splitter() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_query_syng_subwindow");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    // FASTA with two haplotypes sharing a 15 kbp backbone. FastGA needs
    // several kbp of alignable sequence per partition to find hits, so we
    // size this to 15 kbp shared + 2 kbp unique tail = 17 kbp/hap total.
    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        let backbone = numeric_to_ascii(&make_sequence_numeric(15000, 42));
        let tail1 = numeric_to_ascii(&make_sequence_numeric(2000, 1));
        let tail2 = numeric_to_ascii(&make_sequence_numeric(2000, 2));
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail1).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail2).unwrap();
        writeln!(f).unwrap();
    }

    // Build syng index
    let syng_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args(["syng", "-f", fasta_path.to_str().unwrap(), "-o", syng_prefix.to_str().unwrap()])
        .output()
        .expect("Failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    // Query 15 kbp of sampleA's backbone with seqwish:5000 — three sub-windows
    // of 5 kbp each. Each sub-window pulls homologs from both haplotypes,
    // producing ~10 kbp of partition sequence for FastGA/seqwish to chew on.
    let out_prefix = dir.join("region");
    let out = Command::new(&bin)
        .args([
            "query",
            "--syng", syng_prefix.to_str().unwrap(),
            "--sequence-files", fasta_path.to_str().unwrap(),
            "-r", "sampleA#0#chr1:0-15000",
            "-o", "gfa",
            "--gfa-engine", "seqwish:5000",
            "-O", out_prefix.to_str().unwrap(),
            "-t", "1",
            "-v", "2",  // info-level logging to capture sub-window log lines
        ])
        .env("RUST_LOG", "info")
        .output()
        .expect("Failed to run impg query");

    let stderr = String::from_utf8_lossy(&out.stderr);

    // Primary assertion: sub-window wiring — three log lines must appear.
    // This is the regression check for the pggb:X / seqwish:X no-op bug.
    let subwindow_count = stderr
        .lines()
        .filter(|l| l.contains("[syng sub-window"))
        .count();
    assert_eq!(
        subwindow_count, 3,
        "expected 3 sub-window log lines for 15kbp query with seqwish:5000, got {}. stderr: {}",
        subwindow_count, stderr
    );

    // If the downstream pipeline succeeded, also verify the GFA is non-empty.
    // A failure at the FastGA/seqwish stage doesn't invalidate the
    // sub-windowing wiring test above.
    if out.status.success() {
        let gfa_path = format!("{}.gfa", out_prefix.to_str().unwrap());
        let gfa_size = std::fs::metadata(&gfa_path)
            .expect("GFA output missing")
            .len();
        assert!(
            gfa_size > 100,
            "GFA output is only {} bytes — expected a non-trivial graph",
            gfa_size
        );
        let gfa_content = std::fs::read_to_string(&gfa_path).unwrap();
        let s_lines = gfa_content.lines().filter(|l| l.starts_with("S\t")).count();
        assert!(s_lines > 0, "GFA has no S lines");
    } else {
        eprintln!(
            "Sub-windowing OK but downstream engine pipeline failed; \
             this test only asserts sub-window wiring. stderr: {}",
            stderr
        );
    }

    std::fs::remove_dir_all(&dir).ok();
}
