use std::{env, error::Error, fs, io, path::PathBuf, process::Command, time::Instant};

use impg::smooth::{smooth_gfa, SmoothBlockSource, SmoothConfig};

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args = env::args().collect::<Vec<_>>();
    if args.len() < 3 || args.len() > 6 {
        eprintln!(
            "usage: {} <input.gfa> <output-prefix> [iterations=3] [target-bp=10000] [threads=32]",
            args.first()
                .map(String::as_str)
                .unwrap_or("neighbor_merge_existing_gfa")
        );
        std::process::exit(2);
    }

    let input = &args[1];
    let output_prefix = PathBuf::from(&args[2]);
    let iterations = parse_positive_arg(args.get(3), 3, "iterations")?;
    let target_bp = parse_positive_arg(args.get(4), 10_000, "target-bp")?;
    let threads = parse_positive_arg(args.get(5), 32, "threads")?;

    let t0 = Instant::now();
    let gfa = fs::read_to_string(input)?;
    let n_haps = count_pansn_haplotypes(&gfa).max(1);
    log::info!(
        "[neighbor-merge-existing-gfa] input={} n_haps={} iterations={} target_bp={} threads={}",
        input,
        n_haps,
        iterations,
        target_bp,
        threads
    );

    let config = SmoothConfig {
        num_threads: threads,
        target_poa_lengths: vec![target_bp; iterations],
        max_node_length: 100,
        poa_padding_fraction: 0.001,
        pre_sorted: false,
        block_source: SmoothBlockSource::NeighborMergePoasta,
        ..SmoothConfig::new(n_haps)
    };

    let smoothed = smooth_gfa(&gfa, &config)?;
    let raw_path = output_prefix.with_extension("raw-smooth.gfa");
    fs::write(&raw_path, smoothed.as_bytes())?;
    log::info!(
        "[neighbor-merge-existing-gfa] wrote {} in {:.3}s",
        raw_path.display(),
        t0.elapsed().as_secs_f64()
    );

    let gfaffix_start = Instant::now();
    let fixed = run_gfaffix(&smoothed, threads)?;
    let gfaffix_path = output_prefix.with_extension("gfaffix.gfa");
    fs::write(&gfaffix_path, fixed.as_bytes())?;
    log::info!(
        "[neighbor-merge-existing-gfa] wrote {} after gfaffix in {:.3}s",
        gfaffix_path.display(),
        gfaffix_start.elapsed().as_secs_f64()
    );

    let sort_start = Instant::now();
    let sorted = impg::graph::sort_gfa_pipeline(&fixed, "Ygs", threads)?;
    let sorted_path = output_prefix.with_extension("Ygs.gfa");
    fs::write(&sorted_path, sorted.as_bytes())?;
    log::info!(
        "[neighbor-merge-existing-gfa] wrote {} after Ygs sort in {:.3}s (total {:.3}s)",
        sorted_path.display(),
        sort_start.elapsed().as_secs_f64(),
        t0.elapsed().as_secs_f64()
    );

    Ok(())
}

fn run_gfaffix(gfa: &str, threads: usize) -> Result<String, Box<dyn Error>> {
    let exe = env::current_exe()?;
    let mut candidates = vec![exe.with_file_name("gfaffix")];
    if let Some(release_dir) = exe.parent().and_then(|dir| dir.parent()) {
        candidates.push(release_dir.join("gfaffix"));
    }
    let gfaffix = candidates
        .into_iter()
        .find(|path| path.exists())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                "gfaffix binary not found next to the current executable or its parent directory",
            )
        })?;

    let input = tempfile::Builder::new().suffix(".gfa").tempfile()?;
    fs::write(input.path(), gfa.as_bytes())?;
    let output = tempfile::Builder::new().suffix(".gfa").tempfile()?;

    let status = Command::new(gfaffix)
        .arg(input.path())
        .arg("-o")
        .arg(output.path())
        .arg("-p")
        .arg(threads.to_string())
        .status()?;
    if !status.success() {
        return Err(format!("gfaffix exited with status {:?}", status.code()).into());
    }

    Ok(fs::read_to_string(output.path())?)
}

fn parse_positive_arg(
    value: Option<&String>,
    default: usize,
    name: &str,
) -> Result<usize, Box<dyn Error>> {
    let Some(value) = value else {
        return Ok(default);
    };
    let parsed = value.parse::<usize>()?;
    if parsed == 0 {
        return Err(format!("{name} must be > 0").into());
    }
    Ok(parsed)
}

fn count_pansn_haplotypes(gfa: &str) -> usize {
    sweepga::pansn::count_pansn_keys(
        gfa.lines().filter_map(|line| {
            let mut fields = line.split('\t');
            if fields.next() == Some("P") {
                fields.next()
            } else {
                None
            }
        }),
        sweepga::pansn::PanSnLevel::Haplotype,
    )
}
