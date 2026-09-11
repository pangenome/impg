//! Experimental inference CLI; deliberately independent of the legacy infer stitcher.
use crate::genome_inference::{self as genome, calling, catalog, genotype, sample, threading};
use crate::sample_mem_bwt::invalid;
use crate::syng::{SyncmerParams, SyngIndex};
use clap::{Args, Subcommand};
use std::io;
use std::path::PathBuf;

#[derive(Debug, Args)]
pub struct Common {
    /// Exact syng build prefix (not a sidecar path); freshly built metadata v2
    #[arg(long)]
    panel: String,
    /// New directory; even existing empty directories are refused
    #[arg(long)]
    out_dir: PathBuf,
}
#[derive(Debug, Args)]
pub struct CatalogSource {
    /// Directory of completed BED3 ownership groups, one group per file
    #[arg(long, required_unless_present = "groups", conflicts_with = "groups")]
    bed_dir: Option<PathBuf>,
    /// Versioned JSON ownership groups with optional explicit scaffold/strand
    #[arg(long)]
    groups: Option<PathBuf>,
}
#[derive(Debug, Args)]
pub struct CallOptions {
    /// Acknowledge bounds/anchors are NOT validated homology, alleles or copy structure
    #[arg(long)]
    allow_unvalidated_catalog: bool,
    /// Only 1 supported; interpretation depends on the selected command/model
    #[arg(long, default_value_t = 1)]
    ploidy: usize,
    /// Evaluation-only JSON; opened only after frozen calls are written
    #[arg(long)]
    truth: Option<PathBuf>,
}
#[derive(Debug, Subcommand)]
pub enum Command {
    /// Build reusable weighted MEM-substring FM counts (no read IDs or locations)
    BuildSample {
        #[command(flatten)]
        common: Common,
        #[arg(long, required = true, num_args = 1..)]
        reads: Vec<PathBuf>,
    },
    /// Preserve source occurrences and audit global two-anchor feature ownership
    BuildCatalog {
        #[command(flatten)]
        common: Common,
        #[command(flatten)]
        source: CatalogSource,
    },
    /// Provisional presence-compatibility of source occurrences, NOT genotyping
    Call {
        #[command(flatten)]
        common: Common,
        #[arg(long)]
        sample: PathBuf,
        #[arg(long)]
        catalog: PathBuf,
        #[command(flatten)]
        options: CallOptions,
    },
    /// Experimental quantitative haploid source bundles, optionally threaded on an explicit axis
    Genotype {
        #[command(flatten)]
        common: Common,
        #[arg(long)]
        sample: PathBuf,
        #[arg(long)]
        catalog: PathBuf,
        #[command(flatten)]
        options: CallOptions,
        /// Independently declared haploid sequencing depth, never estimated from truth
        #[arg(long)]
        haploid_depth: f64,
        /// Positive expected background MEM-substring count per factor
        #[arg(long, default_value_t = 0.1)]
        background: f64,
        /// Working mean Poisson deviance cutoff; not a calibrated significance level
        #[arg(long, default_value_t = 10.0)]
        max_mean_deviance: f64,
        /// Version-1 reference axis; omitted groups still have independent partition calls
        #[arg(long)]
        axis: Option<PathBuf>,
        /// Cost of source-path/orientation/nonmonotonic changes; reset at chromosomes
        #[arg(long, default_value_t = 10.0)]
        switch_penalty: f64,
        /// Debug only: emit reproducible per-factor counts/exposures and bundle multiplicities (potentially huge)
        #[arg(long)]
        include_feature_details: bool,
    },
    /// One-command bootstrap: reads -> sample index -> catalog -> diagnostics -> evaluation
    Run {
        #[command(flatten)]
        common: Common,
        #[arg(long, required = true, num_args = 1..)]
        reads: Vec<PathBuf>,
        #[command(flatten)]
        source: CatalogSource,
        #[command(flatten)]
        options: CallOptions,
    },
}
fn input(source: CatalogSource) -> io::Result<catalog::CatalogInput> {
    match (source.bed_dir, source.groups) {
        (Some(path), None) => catalog::import_beds(&path),
        (None, Some(path)) => genome::read_json(&path),
        _ => Err(invalid("exactly one of --bed-dir or --groups required")),
    }
}
fn finish_calls(
    out: &std::path::Path,
    catalog: &catalog::Catalog,
    sample: &sample::SampleIndex,
    options: CallOptions,
) -> io::Result<serde_json::Value> {
    let calls = calling::call(
        catalog,
        sample,
        options.allow_unvalidated_catalog,
        options.ploidy,
    )?;
    genome::write_json(&out.join("calls.json"), &calls)?;
    if let Some(path) = options.truth {
        let truth = genome::read_json(&path)?;
        let evaluation = calling::evaluate(catalog, &calls, truth)?;
        genome::write_json(&out.join("evaluation.json"), &evaluation)?;
    }
    Ok(
        serde_json::json!({"groups": calls.calls.len(), "global_factors_used": calls.global_factors_used,
        "sample": sample.stats, "catalog_accepted": false, "scope": "bootstrap-source-occurrence-diagnostics-only"}),
    )
}
pub fn run(command: Command) -> io::Result<()> {
    match command {
        Command::BuildSample { common, reads } => genome::with_output(&common.out_dir, || {
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let panel = SyngIndex::load(&common.panel, SyncmerParams::default())?;
            let sample = sample::build(&panel, identity, &reads)?;
            let path = common.out_dir.join("sample.membwt");
            sample.save(&path)?;
            Ok(
                serde_json::json!({"sample": sample.stats, "sample_index_bytes": std::fs::metadata(path)?.len()}),
            )
        }),
        Command::BuildCatalog { common, source } => genome::with_output(&common.out_dir, || {
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let panel = SyngIndex::load(&common.panel, SyncmerParams::default())?;
            let catalog = catalog::build(&panel, identity, input(source)?)?;
            let path = common.out_dir.join("catalog.json");
            genome::save_catalog(&path, &catalog)?;
            Ok(
                serde_json::json!({"groups": catalog.groups.len(), "occurrences": catalog.occurrences.len(), "features": catalog.features.len(), "catalog_bytes": std::fs::metadata(path)?.len(), "catalog_accepted": false}),
            )
        }),
        Command::Call {
            common,
            sample,
            catalog,
            options,
        } => genome::with_output(&common.out_dir, || {
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let catalog = genome::load_catalog(&catalog)?;
            if catalog.panel != identity {
                return Err(invalid("catalog/current panel fingerprint mismatch"));
            }
            let sample = sample::SampleIndex::load(&sample, &identity)?;
            finish_calls(&common.out_dir, &catalog, &sample, options)
        }),
        Command::Genotype {
            common,
            sample,
            catalog,
            options,
            haploid_depth,
            background,
            max_mean_deviance,
            axis,
            switch_penalty,
            include_feature_details,
        } => genome::with_output_model(&common.out_dir, genotype::MODEL, || {
            let parameters = genotype::Parameters {
                haploid_depth,
                background,
                max_mean_deviance,
            };
            parameters.validate()?;
            if !switch_penalty.is_finite() || switch_penalty < 0.0 {
                return Err(invalid("switch penalty must be finite and nonnegative"));
            }
            if options.ploidy != 1 || !options.allow_unvalidated_catalog {
                return Err(invalid(
                    "quantitative mode requires ploidy 1 and --allow-unvalidated-catalog",
                ));
            }
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let (catalog, catalog_checksum) = genome::load_catalog_with_checksum(&catalog)?;
            if catalog.panel != identity {
                return Err(invalid("catalog/current panel fingerprint mismatch"));
            }
            let (sample, sample_checksum) =
                sample::SampleIndex::load_with_checksum(&sample, &identity)?;
            let started = std::time::Instant::now();
            let mut calls = genotype::call(
                &catalog,
                &sample,
                parameters,
                options.allow_unvalidated_catalog,
                options.ploidy,
            )?;
            calls.sample_payload_checksum = Some(sample_checksum);
            calls.catalog_payload_checksum = Some(catalog_checksum);
            if !include_feature_details {
                calls.omit_feature_details();
            }
            genome::write_json(&common.out_dir.join("calls.json"), &calls)?;
            let call_seconds = started.elapsed().as_secs_f64();
            let started = std::time::Instant::now();
            let threads = axis
                .map(|path| {
                    let axis = genome::read_json(&path)?;
                    let threads = threading::thread(&catalog, &calls, axis, switch_penalty)?;
                    genome::write_json(&common.out_dir.join("threads.json"), &threads)?;
                    Ok::<_, io::Error>(threads)
                })
                .transpose()?;
            let thread_seconds = started.elapsed().as_secs_f64();
            // Both artifacts are frozen before even opening truth.
            if let Some(path) = options.truth {
                let truth = genome::read_json(&path)?;
                let evaluation = genotype::evaluate(&catalog, &calls, threads.as_ref(), truth)?;
                genome::write_json(&common.out_dir.join("evaluation.json"), &evaluation)?;
            }
            Ok(
                serde_json::json!({"model": genotype::MODEL, "groups": calls.calls.len(),
                "global_factors_used": calls.global_factors_used, "excluded_features": calls.excluded_features,
                "call_and_save_seconds": call_seconds, "thread_and_save_seconds": thread_seconds,
                "resolved_axis_intervals": threads.as_ref().map(|t| t.resolved_intervals),
                "unresolved_axis_intervals": threads.as_ref().map(|t| t.unresolved_intervals),
                "calls_bytes": std::fs::metadata(common.out_dir.join("calls.json"))?.len(),
                "sample": sample.stats, "catalog_accepted": false}),
            )
        }),
        Command::Run {
            common,
            reads,
            source,
            options,
        } => genome::with_output(&common.out_dir, || {
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let panel = SyngIndex::load(&common.panel, SyncmerParams::default())?;
            let started = std::time::Instant::now();
            let sample = sample::build(&panel, identity.clone(), &reads)?;
            let sample_path = common.out_dir.join("sample.membwt");
            sample.save(&sample_path)?;
            let sample_seconds = started.elapsed().as_secs_f64();
            let started = std::time::Instant::now();
            let catalog = catalog::build(&panel, identity, input(source)?)?;
            let catalog_path = common.out_dir.join("catalog.json");
            genome::save_catalog(&catalog_path, &catalog)?;
            let catalog_seconds = started.elapsed().as_secs_f64();
            let started = std::time::Instant::now();
            let mut stats = finish_calls(&common.out_dir, &catalog, &sample, options)?;
            stats["call_and_evaluation_seconds"] = started.elapsed().as_secs_f64().into();
            stats["sample_build_and_save_seconds"] = sample_seconds.into();
            stats["catalog_build_and_save_seconds"] = catalog_seconds.into();
            stats["sample_index_bytes"] = std::fs::metadata(sample_path)?.len().into();
            stats["catalog_bytes"] = std::fs::metadata(catalog_path)?.len().into();
            Ok(stats)
        }),
    }
}
