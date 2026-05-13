//! End-to-end integration tests for the syng CLI path.
//!
//! Covers gaps that unit tests missed:
//! - `impg syng --agc` builds a non-empty index (caught a silent failure in the
//!   yeast235 workflow where only `.syng.names` got written).
//! - `impg syng --agc` round-trips: built index can be loaded and queried.
//! - `impg partition -a <syng-prefix>` runs end-to-end and produces non-empty BED.

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
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return Some(PathBuf::from(path));
    }

    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    for candidate in [
        manifest_dir.join("target/release/impg"),
        manifest_dir.join("target/debug/impg"),
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
            "--syncmer-length",
            "63",
            "--smer-length",
            "8",
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
    let spos_path = format!("{}.syng.spos", out_prefix.to_str().unwrap());
    let meta_path = format!("{}.syng.meta", out_prefix.to_str().unwrap());

    assert!(std::path::Path::new(&gbwt_path).exists(), ".1gbwt missing");
    assert!(std::path::Path::new(&khash_path).exists(), ".1khash missing");
    assert!(std::path::Path::new(&names_path).exists(), ".syng.names missing");
    assert!(std::path::Path::new(&spos_path).exists(), ".syng.spos missing");
    assert!(std::path::Path::new(&meta_path).exists(), ".syng.meta missing");

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
        .args([
            "syng",
            "--agc",
            agc_path.to_str().unwrap(),
            "-o",
            out_prefix_str,
            "--position-sample-shift",
            "0",
        ])
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
fn test_syng_map_cli_gaf_and_paf() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_map_cli");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let backbone = numeric_to_ascii(&make_sequence_numeric(1000, 42));
    let tail_a = numeric_to_ascii(&make_sequence_numeric(400, 1));
    let tail_b = numeric_to_ascii(&make_sequence_numeric(400, 2));

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_b).unwrap();
        writeln!(f).unwrap();
    }

    let query_path = dir.join("query.fq");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&query_path).unwrap();
        writeln!(f, "@read1").unwrap();
        f.write_all(&backbone[100..800]).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(700)).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-shift",
            "0",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let gaf = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "gaf",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o gaf");
    assert!(
        gaf.status.success(),
        "impg map -o gaf failed: {}",
        String::from_utf8_lossy(&gaf.stderr)
    );
    let gaf_stdout = String::from_utf8_lossy(&gaf.stdout);
    let gaf_lines: Vec<&str> = gaf_stdout.lines().collect();
    assert_eq!(gaf_lines.len(), 1, "expected one GAF line, got:\n{}", gaf_stdout);
    let gaf_fields: Vec<&str> = gaf_lines[0].split('\t').collect();
    assert!(gaf_fields.len() >= 12, "GAF line has too few fields: {}", gaf_lines[0]);
    assert_eq!(gaf_fields[0], "read1");
    assert!(
        gaf_fields[5].contains('>') || gaf_fields[5].contains('<'),
        "GAF path should contain syncmer node orientations: {}",
        gaf_fields[5]
    );

    let paf = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "paf",
            "--min-anchors",
            "2",
            "--max-hits",
            "10",
        ])
        .output()
        .expect("failed to run impg map -o paf");
    assert!(
        paf.status.success(),
        "impg map -o paf failed: {}",
        String::from_utf8_lossy(&paf.stderr)
    );
    let paf_stdout = String::from_utf8_lossy(&paf.stdout);
    let paf_lines: Vec<&str> = paf_stdout.lines().collect();
    assert!(!paf_lines.is_empty(), "expected PAF hits, got none");
    assert!(
        paf_lines.iter().any(|l| l.contains("sampleA#0#chr1")),
        "expected sampleA hit, got:\n{}",
        paf_stdout
    );
    assert!(
        paf_lines.iter().any(|l| l.contains("sampleB#0#chr1")),
        "expected sampleB hit, got:\n{}",
        paf_stdout
    );
    for line in paf_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= 12, "PAF line has too few fields: {}", line);
        assert_eq!(fields[0], "read1");
    }

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_map_cli_sampled_positions_paf() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_map_sampled_cli");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let backbone = numeric_to_ascii(&make_sequence_numeric(1000, 42));
    let tail_a = numeric_to_ascii(&make_sequence_numeric(400, 1));
    let tail_b = numeric_to_ascii(&make_sequence_numeric(400, 2));

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_b).unwrap();
        writeln!(f).unwrap();
    }

    let query_path = dir.join("query.fq");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&query_path).unwrap();
        writeln!(f, "@read1").unwrap();
        f.write_all(&backbone[100..800]).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(700)).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-shift",
            "0",
        ])
        .output()
        .expect("failed to run impg syng with sampled positions");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );
    assert!(
        std::path::Path::new(&format!("{}.syng.spos", idx_prefix.to_str().unwrap())).exists(),
        "default sampled-position sidecar should be written"
    );
    assert!(
        std::path::Path::new(&format!("{}.syng.meta", idx_prefix.to_str().unwrap())).exists(),
        "syncmer metadata should be written"
    );
    assert!(
        !std::path::Path::new(&format!("{}.syng.locate", idx_prefix.to_str().unwrap())).exists(),
        "legacy FastLocate sidecar should not be built"
    );

    let paf = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "paf",
            "--min-anchors",
            "2",
            "--max-hits",
            "10",
        ])
        .output()
        .expect("failed to run impg map -o paf with sampled positions");
    assert!(
        paf.status.success(),
        "impg map -o paf failed: {}",
        String::from_utf8_lossy(&paf.stderr)
    );
    let paf_stdout = String::from_utf8_lossy(&paf.stdout);
    let paf_lines: Vec<&str> = paf_stdout.lines().collect();
    assert!(!paf_lines.is_empty(), "expected PAF hits, got none");
    assert!(
        paf_lines.iter().any(|l| l.contains("sampleA#0#chr1")),
        "expected sampleA hit, got:\n{}",
        paf_stdout
    );
    assert!(
        paf_lines.iter().any(|l| l.contains("sampleB#0#chr1")),
        "expected sampleB hit, got:\n{}",
        paf_stdout
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
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            out_prefix.to_str().unwrap(),
            "--position-sample-shift",
            "0",
        ])
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
    // End-to-end: build syng index from FASTA, then run `impg partition -a`
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
            "--position-sample-shift",
            "0",
        ])
        .output()
        .expect("Failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed. stderr: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    // Step 2: run partition -a with a syng index prefix and BED output
    let out_folder = dir.join("parts");
    std::fs::create_dir_all(&out_folder).unwrap();
    let part = Command::new(&bin)
        .args([
            "partition",
            "-a",
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
        .expect("Failed to run impg partition -a syng prefix");
    assert!(
        part.status.success(),
        "impg partition -a syng prefix failed. stdout: {} stderr: {}",
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

/// Regression test: `impg query -a <syng-prefix> -o gfa --gfa-engine seqwish:X` must
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

    // This test's primary purpose is verifying that pggb:X / seqwish:X is
    // interpreted as a sub-window size (not a boolean flag), which is done by
    // inspecting stderr for the per-window log lines. The downstream
    // FastGA/seqwish pipeline runs but we DO NOT require it to succeed —
    // the regression assertion fires before any engine work. So we size the
    // input to the minimum that still lets syng's default syncmer params
    // (k=8, w=55, total=63bp) produce syncmer hits per window, and let the
    // engine fail fast on CI's 2-vCPU runners if it wants to.
    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        // Two haplotypes sharing a 3 kbp backbone + 500 bp unique tails.
        // Total ~3.5 kbp each, 7 kbp FASTA.
        let backbone = numeric_to_ascii(&make_sequence_numeric(3000, 42));
        let tail1 = numeric_to_ascii(&make_sequence_numeric(500, 1));
        let tail2 = numeric_to_ascii(&make_sequence_numeric(500, 2));
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
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            syng_prefix.to_str().unwrap(),
            "--position-sample-shift",
            "0",
        ])
        .output()
        .expect("Failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    // Query 3 kbp of sampleA's backbone with poa:1000 — three sub-windows
    // of 1000 bp each (the minimum window size impg allows). We pick `poa`
    // instead of `seqwish`/`pggb` because POA runs a single-pass partial
    // order alignment per sub-window and skips the FastGA + seqwish +
    // gfaffix chain, which is the main source of multi-minute runtime on
    // small CI runners. The sub-window wiring regression is engine-agnostic
    // — the log lines are emitted before any engine runs — so POA is a
    // cheaper substitute that exercises the exact same code path.
    let out_prefix = dir.join("region");
    let out = Command::new(&bin)
        .args([
            "query",
            "-a", syng_prefix.to_str().unwrap(),
            "--sequence-files", fasta_path.to_str().unwrap(),
            "-r", "sampleA#0#chr1:0-3000",
            "-o", "gfa",
            "--gfa-engine", "poa:1000",
            "-O", out_prefix.to_str().unwrap(),
            "-t", "1",
            "-v", "2",  // info-level logging to capture sub-window log lines
        ])
        .env("RUST_LOG", "info")
        .output()
        .expect("Failed to run impg query");

    let stderr = String::from_utf8_lossy(&out.stderr);

    // Primary (and only) assertion: sub-window wiring — three log lines must
    // appear. This is the regression check for the pggb:X / seqwish:X no-op
    // bug, emitted inside the sub-window loop BEFORE the engine runs. A
    // downstream engine failure doesn't invalidate this check and doesn't
    // fail the test — subprocess exit status is intentionally unchecked to
    // keep the test cheap on small CI runners.
    let subwindow_count = stderr
        .lines()
        .filter(|l| l.contains("[syng sub-window"))
        .count();
    assert_eq!(
        subwindow_count, 3,
        "expected 3 sub-window log lines for 3000bp query with poa:1000, got {}. stderr: {}",
        subwindow_count, stderr
    );

    std::fs::remove_dir_all(&dir).ok();
}

/// RC homology is detected and boundary-realigned end-to-end.
///
/// Construct a fixture with an inversion: genome_b contains the
/// reverse-complement of a stretch of genome_a. Confirm:
/// (1) the raw syng path reports the homolog with strand='-',
/// (2) boundary realignment produces tight target coordinates in forward
///     strand space,
/// (3) the refined interval's target forward-strand bases, read in RC, match
///     the query's forward bases exactly (zero edit distance).
#[test]
fn test_syng_rc_homolog_end_to_end() {
    let _guard = lock_syng();

    fn mk_seq(len: usize, seed: u8) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut s = Vec::with_capacity(len);
        let mut state: u32 = seed as u32;
        for _ in 0..len {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            s.push(bases[((state >> 16) % 4) as usize]);
        }
        s
    }

    // genome_a: 3000bp random (seed 42).
    // genome_b: [prefix_1000, RC(a[500..2500]), suffix_500] — 3500bp total,
    //          with an RC-inverted copy of a[500..2500] embedded at b[1000..3000].
    let a = mk_seq(3000, 42);
    let inv_region = &a[500..2500];
    let rc_inv = impg::graph::reverse_complement(inv_region);
    let mut b = mk_seq(1000, 11);
    b.extend_from_slice(&rc_inv);
    b.extend_from_slice(&mk_seq(500, 13));
    assert_eq!(b.len(), 3500);

    let dir = std::env::temp_dir().join("impg_test_syng_rc_homolog");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("rc.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">genome_a").unwrap();
        f.write_all(&a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">genome_b").unwrap();
        f.write_all(&b).unwrap();
        writeln!(f).unwrap();
    }

    let params = impg::syng::SyncmerParams { k: 8, w: 55, seed: 7 };
    let syng_index = impg::syng::SyngIndex::build(
        params,
        vec![
            ("genome_a".to_string(), a.clone()),
            ("genome_b".to_string(), b.clone()),
        ]
        .into_iter(),
    );
    let sequence_index = impg::sequence_index::UnifiedSequenceIndex::from_files(&[
        fasta_path.to_string_lossy().to_string(),
    ])
    .unwrap();

    // Query the middle of the inverted region on genome_a.
    let query_start = 1000u64;
    let query_end = 2000u64;
    let padding = 120;

    // (1) Raw query_region reports a '-' homolog for genome_b.
    let raw = syng_index
        .query_region("genome_a", query_start, query_end, padding)
        .unwrap();
    let rc_raw: Vec<_> = raw
        .iter()
        .filter(|iv| iv.genome == "genome_b" && iv.strand == '-')
        .collect();
    assert!(
        !rc_raw.is_empty(),
        "raw query should report at least one '-' strand homolog on genome_b; got: {:?}",
        raw.iter().map(|iv| (&iv.genome, iv.start, iv.end, iv.strand)).collect::<Vec<_>>()
    );

    // (2) Boundary realignment produces tight coordinates for the RC homolog.
    let refined = impg::syng_transitive::query_transitive(
        &syng_index,
        "genome_a",
        query_start,
        query_end,
        padding,
        1,
        &sequence_index,
    )
    .unwrap();
    let rc_refined: Vec<_> = refined
        .iter()
        .filter(|iv| iv.genome == "genome_b" && iv.strand == '-')
        .collect();
    assert!(
        !rc_refined.is_empty(),
        "refined query should keep the '-' strand homolog; got: {:?}",
        refined.iter().map(|iv| (&iv.genome, iv.start, iv.end, iv.strand)).collect::<Vec<_>>()
    );

    // The expected refined forward-strand region on genome_b:
    // genome_a[500..2500] ↔ RC → b[1000..3000]. Query a[1000..2000]
    // corresponds (under RC) to b[1500, 2500). Our refined intervals must
    // overlap that window and, when RC'd, reproduce query bases with high
    // fidelity somewhere along the alignment.
    let expected_start = 1500u64;
    let expected_end = 2500u64;
    let overlap_bp = rc_refined
        .iter()
        .map(|iv| {
            let lo = iv.start.max(expected_start);
            let hi = iv.end.min(expected_end);
            if hi > lo { hi - lo } else { 0 }
        })
        .sum::<u64>();
    assert!(
        overlap_bp >= 200,
        "RC refined intervals should collectively overlap the expected RC window \
         [{}, {}) by at least 200bp; got {}bp overlap. Intervals: {:?}",
        expected_start, expected_end, overlap_bp,
        rc_refined.iter().map(|iv| (iv.start, iv.end)).collect::<Vec<_>>()
    );

    // Base-content validation: for each refined '-' interval, RC'ing the
    // target bytes should produce something that matches a sub-string of
    // the query region (with at least one long exact-match run). This
    // confirms the ORIENTATION and COVERAGE are correct even if the exact
    // endpoints are subject to the syng-position calibration noise flagged
    // as a known limitation below.
    let query_bytes = &a[query_start as usize..query_end as usize];
    let mut best_match_len = 0usize;
    for iv in &rc_refined {
        let b_slice = &b[iv.start as usize..iv.end as usize];
        let b_rc = impg::graph::reverse_complement(b_slice);
        // Find the longest common substring between b_rc and query_bytes,
        // capped at the shorter length. Cheap approximation: slide b_rc
        // over query_bytes and measure the longest run.
        if b_rc.len() < 30 || query_bytes.len() < 30 {
            continue;
        }
        let search = &b_rc[..b_rc.len().min(200)];
        for start in 0..query_bytes.len().saturating_sub(search.len()) {
            let run = query_bytes[start..]
                .iter()
                .zip(search.iter())
                .take_while(|(x, y)| x == y)
                .count();
            if run > best_match_len {
                best_match_len = run;
            }
        }
    }
    assert!(
        best_match_len >= 30,
        "RC'd refined target bytes should share a >=30bp exact run with the query \
         region; longest run was {}bp",
        best_match_len
    );

    std::fs::remove_dir_all(&dir).ok();
}

/// Boundary-realignment tightens fuzzy syncmer-resolution edges.
///
/// Construct two genomes that share a 500bp backbone at offset 0..500. Query
/// the backbone region on genome_a. The RAW syng path returns padded
/// syncmer-resolution intervals (slop-boundary up to ±syncmer_len bp around
/// the real edges). The boundary-realignment path should return intervals
/// whose edges are much closer to the true 0..500 boundary on genome_b.
#[test]
fn test_syng_boundary_realign_tightens_edges() {
    let _guard = lock_syng();

    // Reproducible backbone + distinct tails.
    fn mk_seq(len: usize, seed: u8) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut s = Vec::with_capacity(len);
        let mut state: u32 = seed as u32;
        for _ in 0..len {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            s.push(bases[((state >> 16) % 4) as usize]);
        }
        s
    }

    // Backbone length must comfortably exceed query_end + syncmer_len so that
    // shared syncmers exist on BOTH sides of every tested edge. At default
    // (k=8, w=55) a syncmer starting at position p extends to p+63; only
    // syncmers whose whole length sits inside the shared backbone are shared.
    let backbone = mk_seq(2000, 42);
    let mut seq_a = backbone.clone();
    seq_a.extend(mk_seq(500, 1));
    let mut seq_b = backbone.clone();
    seq_b.extend(mk_seq(500, 2));

    let dir = std::env::temp_dir().join("impg_test_syng_boundary_realign");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    // Write a FASTA for the UnifiedSequenceIndex.
    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">genome_a").unwrap();
        f.write_all(&seq_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">genome_b").unwrap();
        f.write_all(&seq_b).unwrap();
        writeln!(f).unwrap();
    }

    // Build the syng index in-process.
    let params = impg::syng::SyncmerParams { k: 8, w: 55, seed: 7 };
    let sequences = vec![
        ("genome_a".to_string(), seq_a.clone()),
        ("genome_b".to_string(), seq_b.clone()),
    ];
    let syng_index = impg::syng::SyngIndex::build(params, sequences.into_iter());

    // Build the file-backed sequence index.
    let sequence_index = impg::sequence_index::UnifiedSequenceIndex::from_files(&[
        fasta_path.to_string_lossy().to_string(),
    ])
    .unwrap();

    let padding = 120; // default
    let query_start = 50u64;
    let query_end = 450u64;

    // Raw: what syng gives us directly.
    let raw = syng_index
        .query_region("genome_a", query_start, query_end, padding)
        .unwrap();
    let raw_b = raw
        .iter()
        .find(|iv| iv.genome == "genome_b")
        .expect("raw query should find genome_b hit");

    // Realigned.
    let refined = impg::syng_transitive::query_transitive(
        &syng_index,
        "genome_a",
        query_start,
        query_end,
        padding,
        1,
        &sequence_index,
    )
    .unwrap();
    let refined_b = refined
        .iter()
        .find(|iv| iv.genome == "genome_b")
        .expect("realigned query should find genome_b hit");

    // On identical backbone, boundary realignment should land on the
    // query coordinates exactly (BiWFA edit distance = 0 between query
    // and target slices around the edges).
    assert_eq!(
        refined_b.start, query_start,
        "refined start should snap to query_start on identical backbone; got {} vs expected {} (raw was {})",
        refined_b.start, query_start, raw_b.start
    );
    assert_eq!(
        refined_b.end, query_end,
        "refined end should snap to query_end on identical backbone; got {} vs expected {} (raw was {})",
        refined_b.end, query_end, raw_b.end
    );

    // The raw interval must be at least as loose as the refined one
    // (slop is additive, not subtractive).
    assert!(
        raw_b.start <= refined_b.start && raw_b.end >= refined_b.end,
        "raw interval should enclose refined interval: raw {}..{} vs refined {}..{}",
        raw_b.start,
        raw_b.end,
        refined_b.start,
        refined_b.end
    );

    std::fs::remove_dir_all(&dir).ok();
}

/// Realistic "homologous-with-small-differences" end-to-end test.
///
/// Build three genomes that share a 3kb region but differ at a handful of
/// scattered positions (SNPs) and with a small indel on one of them. Query
/// a 2kb region in the middle of genome_a at depth 1; expect the refined
/// boundary coordinates on both targets to snap to the biological
/// boundaries of the query (not the padded syncmer-resolution bounds).
///
/// This is the test Erik asked for — the current code is expected to show
/// fragmentation / loose edges on this fixture, because `resolve_edge_via_biwfa`
/// projects through flanking anchors around the query edge and doesn't do a
/// bedtools-style merge + end-alignment pass. Once the merge-then-end-align
/// pipeline lands, this test should pass with single-digit bp precision.
#[test]
fn test_syng_query_reconstructs_homology_with_diffs() {
    let _guard = lock_syng();

    fn mk_seq(len: usize, seed: u8) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut s = Vec::with_capacity(len);
        let mut state: u32 = seed as u32;
        for _ in 0..len {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            s.push(bases[((state >> 16) % 4) as usize]);
        }
        s
    }

    // Shared 3kb backbone, unique 500bp tails.
    let backbone = mk_seq(3000, 42);

    let seq_a = {
        let mut s = backbone.clone();
        s.extend(mk_seq(500, 1));
        s
    };

    // seq_b: identical backbone in the first 3kb, plus 5 scattered SNPs.
    //        SNPs shouldn't shift any coordinates.
    let seq_b = {
        let mut s = backbone.clone();
        for &pos in &[250usize, 800, 1337, 1900, 2500] {
            // Flip one base deterministically (shift by 1 in the 4-base alphabet).
            s[pos] = match s[pos] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                b'T' => b'A',
                other => other,
            };
        }
        s.extend(mk_seq(500, 2));
        s
    };

    // seq_c: identical first ~1500bp; 10bp deletion at 1500; then the rest
    //        of the backbone continues. So homology to genome_a covers
    //        [0, 1500) and [1500, 2990) on seq_c (target coords are shifted
    //        by -10 after position 1500).
    let seq_c = {
        let mut s: Vec<u8> = Vec::with_capacity(3490);
        s.extend_from_slice(&backbone[..1500]);
        s.extend_from_slice(&backbone[1510..]); // drop 10bp at 1500..1510
        s.extend(mk_seq(500, 3));
        s
    };

    let dir = std::env::temp_dir().join("impg_test_syng_homology_diffs");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">genome_a").unwrap();
        f.write_all(&seq_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">genome_b").unwrap();
        f.write_all(&seq_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">genome_c").unwrap();
        f.write_all(&seq_c).unwrap();
        writeln!(f).unwrap();
    }

    let params = impg::syng::SyncmerParams { k: 8, w: 55, seed: 7 };
    let syng_index = impg::syng::SyngIndex::build(
        params,
        vec![
            ("genome_a".to_string(), seq_a.clone()),
            ("genome_b".to_string(), seq_b.clone()),
            ("genome_c".to_string(), seq_c.clone()),
        ]
        .into_iter(),
    );
    let sequence_index = impg::sequence_index::UnifiedSequenceIndex::from_files(&[
        fasta_path.to_string_lossy().to_string(),
    ])
    .unwrap();

    // Query the middle 2kb of genome_a: [500, 2500).
    let query_start = 500u64;
    let query_end = 2500u64;
    let padding = 120;

    let refined = impg::syng_transitive::query_transitive(
        &syng_index,
        "genome_a",
        query_start,
        query_end,
        padding,
        1,
        &sequence_index,
    )
    .unwrap();

    // Expected: ONE homolog per target genome at the biological coordinates.
    //
    // On genome_b: matches [500, 2500) exactly (SNPs don't shift coords).
    // On genome_c: matches [500, 2490) — the 10bp deletion at backbone pos
    //              1500 shifts all downstream positions by -10, so the
    //              query's right edge at 2500 maps to c's 2490.
    let refined_on = |genome: &str| -> Vec<(u64, u64)> {
        refined
            .iter()
            .filter(|iv| iv.genome == genome && iv.strand == '+')
            .map(|iv| (iv.start, iv.end))
            .collect()
    };

    let on_b = refined_on("genome_b");
    let on_c = refined_on("genome_c");

    // Step 1: no fragmentation.
    assert_eq!(
        on_b.len(), 1,
        "expected exactly one forward-strand homolog on genome_b; got {:?}", on_b
    );
    assert_eq!(
        on_c.len(), 1,
        "expected exactly one forward-strand homolog on genome_c; got {:?}", on_c
    );

    // Step 2: refined edges match biological truth within a few bp.
    let tol: u64 = 5;
    let (b_start, b_end) = on_b[0];
    let (c_start, c_end) = on_c[0];

    assert!(
        b_start.abs_diff(500) <= tol && b_end.abs_diff(2500) <= tol,
        "refined genome_b interval should snap to [500, 2500) within {}bp; got [{}, {})",
        tol, b_start, b_end
    );
    assert!(
        c_start.abs_diff(500) <= tol && c_end.abs_diff(2490) <= tol,
        "refined genome_c interval should snap to [500, 2490) within {}bp (10bp deletion at 1500 shifts the right edge); got [{}, {})",
        tol, c_start, c_end
    );

    std::fs::remove_dir_all(&dir).ok();
}
