//! Regression test for the `sv-classify` subcommand (see src/commands/sv_classify.rs).
//!
//! Fixture: `tests/test_data/yeast-toy-35-sim-td.paf` — S288C#0 aligned against
//! TOY288#0 (as target=TOY288#0, query=S288C#0) on the toy yeast genome, carrying
//! known DEL, INS, strand-based INV, MUM&Co-style TCON/TDUP, and a genuine
//! cross-chromosome TRA (chrX <-> chrXI) at fixed loci. This locks in the whole
//! classification pipeline (gap-event merging, strand-based INV detection,
//! tandem-overlap detection, cross-chromosome TRA) end to end through the CLI,
//! not just the underlying functions.

use std::path::PathBuf;
use std::process::Command;

fn impg_binary() -> PathBuf {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return PathBuf::from(path);
    }
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    for candidate in [
        manifest_dir.join("target/release/impg"),
        manifest_dir.join("target/debug/impg"),
    ] {
        if candidate.exists() {
            return candidate;
        }
    }
    panic!("could not locate the impg binary (set CARGO_BIN_EXE_impg or build it first)");
}

const EXPECTED_OUTPUT: &str = "\
chrom\tstart\tend\tsv_type\tsize\tsupport\tquery_chrom\tquery_start\tquery_end
TOY288#0#chrI\t190000\t195000\tINV\t5000\t1\tS288C#0#chrI\t190000\t195000
TOY288#0#chrI\t196900\t197102\tDEL\t202\t1\tS288C#0#chrI\t196900\t196900
TOY288#0#chrII\t699797\t699797\tINS\t251\t1\tS288C#0#chrII\t699797\t700048
TOY288#0#chrIII\t250049\t250300\tDEL\t251\t1\tS288C#0#chrIII\t250049\t250049
TOY288#0#chrIX\t393762\t393762\tTCON\t4104\t1\tS288C#0#chrIX\t389658\t393762
TOY288#0#chrIX\t397866\t397866\tTCON\t4104\t1\tS288C#0#chrIX\t389658\t393762
TOY288#0#chrIX\t401970\t401970\tTCON\t4104\t1\tS288C#0#chrIX\t389658\t393762
TOY288#0#chrV\t459797\t459797\tINS\t1001\t1\tS288C#0#chrV\t459797\t460798
TOY288#0#chrVI\t162939\t262700\tINV\t99761\t1\tS288C#0#chrVI\t162856\t262696
TOY288#0#chrX\t701749\t753611\tTRA\t51862\t1\tS288C#0#chrXI\t615000\t666862
TOY288#0#chrXI\t614915\t664862\tTRA\t49947\t1\tS288C#0#chrX\t701664\t751611
";

#[test]
fn toy_sv_are_stable() {
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let paf_path = manifest_dir.join("tests/test_data/yeast-toy-35-sim-td.paf");
    assert!(paf_path.exists(), "missing fixture: {}", paf_path.display());

    // A stale sidecar index from a previous run (or another test) would make
    // this test validate the cached index instead of a fresh one; force a
    // rebuild so the test exercises the real query -> classify pipeline.
    let index_path = manifest_dir.join("tests/test_data/yeast-toy-35-sim-td.paf.impg");
    std::fs::remove_file(&index_path).ok();

    let run = Command::new(impg_binary())
        .args([
            "sv-classify",
            "--target-name",
            "TOY288#0",
            "-a",
            paf_path.to_str().unwrap(),
            "--query-name",
            "S288C#0",
        ])
        .output()
        .expect("failed to run impg sv-classify");

    std::fs::remove_file(&index_path).ok();

    assert!(
        run.status.success(),
        "impg sv-classify failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let stdout = String::from_utf8(run.stdout).expect("stdout was not valid UTF-8");
    assert_eq!(stdout, EXPECTED_OUTPUT);
}
