//! Crush perf benches.  Real input, real algorithms, real timing.
//!
//! Default input: `data/c4_crush_eval_20260523T140141Z/C4A.blunt.gfa`
//! (the C4 reference workload from `docs/crush-design.md` §5.2).
//!
//! Override with: `CRUSH_BENCH_GFA=/path/to/blunt.gfa cargo bench --bench crush`.
//!
//! Benches measure the three optimizable hot paths called out in
//! `docs/crush-design.md` §2.1 ((2) per-round POVU input prep, (3) per-round
//! render/parse, (4) path-sequence validation).  End-to-end crush timing
//! (including alignment) is reported separately in `docs/crush-perf-report.md`
//! and is out of scope here.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::fs;
use std::path::PathBuf;
use std::time::Duration;

const DEFAULT_BENCH_GFA: &str = "data/c4_crush_eval_20260523T140141Z/C4A.blunt.gfa";

fn bench_gfa_path() -> PathBuf {
    if let Ok(env) = std::env::var("CRUSH_BENCH_GFA") {
        return PathBuf::from(env);
    }
    // Fall back to absolute path if relative is not present (the workspace's
    // data/ dir lives in the parent worktree and is .gitignore'd).
    let rel = PathBuf::from(DEFAULT_BENCH_GFA);
    if rel.exists() {
        return rel;
    }
    let abs = PathBuf::from("/home/erikg/impg").join(DEFAULT_BENCH_GFA);
    abs
}

fn load_input() -> Option<(PathBuf, String)> {
    let path = bench_gfa_path();
    match fs::read_to_string(&path) {
        Ok(text) => Some((path, text)),
        Err(err) => {
            eprintln!(
                "crush bench: skipping (input GFA not found at {}: {}). \
                 Set CRUSH_BENCH_GFA=/path/to/blunt.gfa to override.",
                path.display(),
                err
            );
            None
        }
    }
}

fn bench_crush(c: &mut Criterion) {
    let Some((path, gfa_text)) = load_input() else {
        return;
    };
    let bytes = gfa_text.len() as u64;
    eprintln!(
        "crush bench input: {} ({:.1} MiB)",
        path.display(),
        bytes as f64 / (1024.0 * 1024.0)
    );

    // Parse once up front so the segment / path summary is printed.
    let (n_seg, n_path) = impg::resolution::bench_parse_gfa(&gfa_text)
        .expect("bench input must be a valid blunt GFA");
    eprintln!("crush bench shape: {} segments, {} paths", n_seg, n_path);

    let mut group = c.benchmark_group("crush");
    group
        .sample_size(10)
        .warm_up_time(Duration::from_secs(2))
        .measurement_time(Duration::from_secs(20))
        .throughput(Throughput::Bytes(bytes));

    group.bench_function(BenchmarkId::new("parse_gfa", "C4A.blunt"), |b| {
        b.iter(|| {
            impg::resolution::bench_parse_gfa(&gfa_text).expect("parse must succeed in bench")
        })
    });

    group.bench_function(
        BenchmarkId::new("parse_render_roundtrip", "C4A.blunt"),
        |b| {
            b.iter(|| {
                impg::resolution::bench_parse_render_gfa(&gfa_text)
                    .expect("parse+render must succeed in bench")
            })
        },
    );

    group.bench_function(
        BenchmarkId::new("validate_streaming_self", "C4A.blunt"),
        |b| {
            b.iter(|| {
                impg::resolution::bench_validate_streaming_self(&gfa_text)
                    .expect("validate must succeed in bench")
            })
        },
    );

    group.finish();
}

criterion_group!(benches, bench_crush);
criterion_main!(benches);
