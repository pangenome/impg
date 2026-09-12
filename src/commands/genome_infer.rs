//! Experimental inference CLI; deliberately independent of the legacy infer stitcher.
use crate::genome_inference::{
    self as genome, calling, catalog, genotype, joint, observations, panel_routes, sample,
    threading,
};
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
#[derive(Debug, Args)]
pub struct JointEvidence {
    #[arg(long)]
    compiled: PathBuf,
    #[arg(long)]
    sample: PathBuf,
    #[arg(long)]
    haploid_depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
}
#[derive(Debug, Args)]
pub struct RouteEvidence {
    #[arg(long)]
    routes: PathBuf,
    #[arg(long)]
    sample: PathBuf,
    #[arg(long)]
    haploid_depth: f64,
    #[arg(long, default_value_t = 0.1)]
    background: f64,
    /// Memory resource failure, never partial feature selection
    #[arg(long, default_value_t = 50000000)]
    max_feature_terms: usize,
    #[arg(long, default_value_t = 1000000)]
    cache_terms: usize,
}
#[derive(Debug, Subcommand)]
pub enum Command {
    /// Automatic all-source lanes, all identity endpoint families, DNA hubs and indexed native start profiles
    BuildPanelRoutes {
        #[command(flatten)]
        common: Common,
        #[arg(long)]
        catalog: PathBuf,
        #[arg(long, required=true, num_args=1..)]
        sources: Vec<String>,
        #[arg(long, required=true, value_delimiter=',', num_args=1..)]
        read_lengths: Vec<u64>,
        #[arg(long, default_value_t = 65536)]
        core_bp: u64,
        #[arg(long, default_value_t = 2000000)]
        max_profile_terms: usize,
    },
    /// Exact lazy full-genome route evaluation; fixed-universe background-relative objective only
    EvaluatePanelRoutes {
        #[command(flatten)]
        common: Common,
        #[command(flatten)]
        evidence: RouteEvidence,
        #[arg(
            long,
            required_unless_present = "native_family",
            conflicts_with = "native_family"
        )]
        assignment: Option<PathBuf>,
        /// Automatically select this identity's complete native path inventory
        #[arg(long)]
        native_family: Option<String>,
    },
    /// All-native initializations plus cross-source routes; identity mixing reported separately (no emission)
    SearchPanelRoutes {
        #[command(flatten)]
        common: Common,
        #[command(flatten)]
        evidence: RouteEvidence,
        #[arg(long, default_value_t = 100000)]
        max_work: u64,
        #[arg(long, default_value_t = 10000)]
        max_evaluations: u64,
        #[arg(long, default_value_t = 1024)]
        max_frontier: usize,
        #[arg(long, default_value_t = 1000)]
        max_optima: usize,
        #[arg(long, default_value_t = 1e-9)]
        tie_epsilon: f64,
    },
    /// Fresh exact native replay of explicit mixed-source linear molecule alternatives
    CompileJointWalks {
        #[command(flatten)]
        common: Common,
        #[arg(long)]
        layout: PathBuf,
        #[arg(long, required=true, num_args=1..)]
        sources: Vec<String>,
        #[arg(long, required=true, value_delimiter=',', num_args=1..)]
        read_lengths: Vec<u64>,
        /// Optional full legacy catalog: definitions only, never conditional profiles
        #[arg(long)]
        registry_catalog: Option<PathBuf>,
    },
    /// Evaluate one complete physical-slot assignment; not unary genotypes or sequence
    EvaluateJointWalks {
        #[command(flatten)]
        common: Common,
        #[command(flatten)]
        evidence: JointEvidence,
        #[arg(long)]
        assignment: PathBuf,
    },
    /// Bounded exhaustive coupled search in a finite explicit layout (not genome-wide completion)
    SolveJointWalks {
        #[command(flatten)]
        common: Common,
        #[command(flatten)]
        evidence: JointEvidence,
        /// Includes infeasible complete assignments; no alternatives silently pruned
        #[arg(long, default_value_t = 100000)]
        max_assignments: u64,
        #[arg(long, default_value_t = 1000)]
        max_optima: usize,
        #[arg(long, default_value_t = 1e-9)]
        tie_epsilon: f64,
    },
    /// Build reusable weighted MEM-substring FM counts (no read IDs or locations)
    BuildSample {
        #[command(flatten)]
        common: Common,
        #[arg(long, required = true, num_args = 1..)]
        reads: Vec<PathBuf>,
    },
    /// Compile reusable partition-major physical incidences and exact source-read profiles
    BuildPartitionObservations {
        #[command(flatten)]
        common: Common,
        #[arg(long)]
        catalog: PathBuf,
        #[arg(long, required=true, num_args=1..)]
        sources: Vec<String>,
        #[arg(long, required=true, value_delimiter=',', num_args=1..)]
        read_lengths: Vec<u64>,
        /// Select all originally owned features, never restrict the panel/source scan
        #[arg(long)]
        feature_group: Vec<String>,
    },
    /// Occurrence-weighted partition calls; partial feature scopes cannot be threaded
    GenotypePartitions {
        #[command(flatten)]
        common: Common,
        #[arg(long)]
        observations: PathBuf,
        #[arg(long)]
        sample: PathBuf,
        #[arg(long)]
        haploid_depth: f64,
        #[arg(long, default_value_t = 0.1)]
        background: f64,
        #[arg(long, default_value_t = 10.0)]
        max_mean_deviance: f64,
        #[arg(long)]
        allow_unvalidated_catalog: bool,
        #[arg(long)]
        axis: Option<PathBuf>,
        #[arg(long, default_value_t = 10.0)]
        switch_penalty: f64,
    },
    /// Rethread frozen full-scope partition calls without reloading legacy catalog/sample
    ThreadPartitions {
        #[command(flatten)]
        common: Common,
        #[arg(long)]
        observations: PathBuf,
        #[arg(long)]
        calls: PathBuf,
        #[arg(long)]
        axis: PathBuf,
        #[arg(long, default_value_t = 10.0)]
        switch_penalty: f64,
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
    /// Spell frozen source-supported blocks (no truth input; not a chromosome assembly)
    Reconstruct {
        #[arg(long)]
        calls: PathBuf,
        #[arg(long)]
        threads: PathBuf,
        /// Exact panel .names sidecar; verified against calls, no graph/catalog reload
        #[arg(long)]
        panel_names: PathBuf,
        #[arg(long, required=true, num_args=1..)]
        sources: Vec<String>,
        /// Split at positive source/reference gaps, or copy same-source gaps as imputed bp
        #[arg(long, default_value="split", value_parser=["split", "copy-source"])]
        gap_policy: String,
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// Independently replay a base-level PAF against whole reconstruction and truth FASTAs
    EvaluateSequence {
        #[arg(long)]
        query: PathBuf,
        #[arg(long)]
        truth: PathBuf,
        #[arg(long)]
        paf: PathBuf,
        #[arg(long)]
        aligner_metadata: Option<PathBuf>,
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// Run inspected wfmash b55cf75 with explicit non-masking settings, then evaluate
    AlignSequence {
        #[arg(long)]
        query: PathBuf,
        #[arg(long)]
        truth: PathBuf,
        /// Exact executable path (fingerprinted); use evaluate-sequence for other tools
        #[arg(long)]
        wfmash: PathBuf,
        #[arg(long, default_value_t = 4)]
        threads: usize,
        #[arg(long)]
        out_dir: PathBuf,
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
fn joint_run(
    common: Common,
    evidence: JointEvidence,
    operation: impl FnOnce(&joint::Problem<'_>) -> io::Result<serde_json::Value>,
) -> io::Result<()> {
    genome::with_output_model(&common.out_dir, joint::MODEL, || {
        let identity = genome::PanelIdentity::read(&common.panel)?;
        let (compiled, compiled_checksum) = joint::Compiled::load(&evidence.compiled, &identity)?;
        let (sample, sample_checksum) =
            sample::SampleIndex::load_with_checksum(&evidence.sample, &identity)?;
        let problem = joint::Problem::new(
            &compiled,
            &sample,
            evidence.haploid_depth,
            evidence.background,
        )?;
        genome::write_json(&common.out_dir.join("factors.json"), &problem.factors)?;
        let result = operation(&problem)?;
        let artifact = serde_json::json!({"version":joint::VERSION,"model":joint::MODEL,
            "compiled_checksum":compiled_checksum,"sample_payload_checksum":sample_checksum,
            "compiler_identity":joint::compiler_identity(),"haploid_depth":evidence.haploid_depth,
            "background":evidence.background,"count_policy":genome::COUNT_POLICY,
            "feature_universe":"registry-union-candidate-generated-union-sample-positive",
            "factor_omissions":"none; constants and unsupported sample residuals explicit",
            "result":result});
        genome::write_json(&common.out_dir.join("result.json"), &artifact)?;
        Ok(
            serde_json::json!({"factors":problem.factors.len(),"scope":"finite-explicit-layout-only",
            "sequence_emission_authorized":false}),
        )
    })
}
fn route_run(
    common: Common,
    evidence: RouteEvidence,
    operation: impl FnOnce(
        &mut panel_routes::Evaluator<'_>,
        &std::path::Path,
    ) -> io::Result<serde_json::Value>,
) -> io::Result<()> {
    genome::with_output_model(&common.out_dir, panel_routes::MODEL, || {
        let identity = genome::PanelIdentity::read(&common.panel)?;
        let graph = panel_routes::Graph::load(&evidence.routes, &identity)?;
        let panel = SyngIndex::load(&common.panel, SyncmerParams::default())?;
        let (sample, sample_checksum) =
            sample::SampleIndex::load_with_checksum(&evidence.sample, &identity)?;
        let mut evaluator = panel_routes::Evaluator::new(
            &evidence.routes,
            &graph,
            &panel,
            &sample,
            evidence.haploid_depth,
            evidence.background,
            evidence.max_feature_terms,
            evidence.cache_terms,
        )?;
        let result = operation(&mut evaluator, &common.out_dir)?;
        let value = serde_json::json!({"version":panel_routes::VERSION,"model":panel_routes::MODEL,
            "graph_checksum":graph.digest()?,"sample_payload_checksum":sample_checksum,
            "compiler_identity":panel_routes::compiler_identity(),"count_policy":genome::COUNT_POLICY,
            "haploid_depth":evidence.haploid_depth,"background":evidence.background,
            "objective_kind":"fixed-universe-background-relative-NLL-only","result":result});
        genome::write_json(&common.out_dir.join("result.json"), &value)?;
        Ok(
            serde_json::json!({"result":genome::reconstruction::fingerprint(&common.out_dir.join("result.json"))?,"sequence_emission_authorized":false}),
        )
    })
}
pub fn run(command: Command) -> io::Result<()> {
    match command {
        Command::BuildPanelRoutes {
            common,
            catalog,
            sources,
            read_lengths,
            core_bp,
            max_profile_terms,
        } => genome::with_output_model(&common.out_dir, panel_routes::MODEL, || {
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let panel = SyngIndex::load(&common.panel, SyncmerParams::default())?;
            panel_routes::build(
                &panel,
                identity,
                &catalog,
                &sources,
                read_lengths,
                &common.out_dir,
                core_bp,
                max_profile_terms,
            )
        }),
        Command::EvaluatePanelRoutes {
            common,
            evidence,
            assignment,
            native_family,
        } => route_run(common, evidence, |e, _out| {
            let assignment = if let Some(path) = assignment {
                genome::read_json(&path)?
            } else {
                let family = e
                    .graph
                    .families
                    .iter()
                    .position(|f| Some(&f.identity) == native_family.as_ref())
                    .ok_or_else(|| invalid("unknown native topology family"))?;
                e.graph.native_assignment(family)?
            };
            serde_json::to_value(e.evaluate(&assignment)?).map_err(io::Error::other)
        }),
        Command::SearchPanelRoutes {
            common,
            evidence,
            max_work,
            max_evaluations,
            max_frontier,
            max_optima,
            tie_epsilon,
        } => route_run(common, evidence, |e, out| {
            let result = panel_routes::search(
                e,
                panel_routes::SearchBudget {
                    max_work,
                    max_evaluations,
                    max_frontier,
                    max_optima,
                    tie_epsilon,
                },
                out,
            )?;
            if let Some(incumbent) = &result.incumbent {
                genome::write_json(
                    &out.join("incumbent-assignment.json"),
                    &incumbent.assignment,
                )?;
                genome::write_json(
                    &out.join("incumbent-evaluation.json"),
                    &e.evaluate(&incumbent.assignment)?,
                )?;
            }
            serde_json::to_value(result).map_err(io::Error::other)
        }),
        Command::CompileJointWalks {
            common,
            layout,
            sources,
            read_lengths,
            registry_catalog,
        } => genome::with_output_model(&common.out_dir, joint::MODEL, || {
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let panel = SyngIndex::load(&common.panel, SyncmerParams::default())?;
            let compiled = joint::compile(
                &panel,
                identity,
                genome::read_json(&layout)?,
                &sources,
                read_lengths,
                registry_catalog.as_deref(),
            )?;
            compiled.save(&common.out_dir.join("joint-profiles.json"))?;
            Ok(
                serde_json::json!({"slots":compiled.layout.slots.len(),"features":compiled.definitions.len(),
                    "scope":"finite-explicit-layout-only","context_complete":true}),
            )
        }),
        Command::EvaluateJointWalks {
            common,
            evidence,
            assignment,
        } => joint_run(common, evidence, |problem| {
            serde_json::to_value(joint::evaluate(problem, &genome::read_json(&assignment)?)?)
                .map_err(io::Error::other)
        }),
        Command::SolveJointWalks {
            common,
            evidence,
            max_assignments,
            max_optima,
            tie_epsilon,
        } => joint_run(common, evidence, |problem| {
            serde_json::to_value(joint::solve(
                problem,
                max_assignments,
                max_optima,
                tie_epsilon,
            )?)
            .map_err(io::Error::other)
        }),
        Command::BuildPartitionObservations {
            common,
            catalog,
            sources,
            read_lengths,
            feature_group,
        } => genome::with_output_model(&common.out_dir, observations::MODEL, || {
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let panel = SyngIndex::load(&common.panel, SyncmerParams::default())?;
            observations::build(
                &panel,
                identity,
                &catalog,
                &sources,
                read_lengths,
                &feature_group,
                &common.out_dir,
            )
        }),
        Command::GenotypePartitions {
            common,
            observations: directory,
            sample: sample_path,
            haploid_depth,
            background,
            max_mean_deviance,
            allow_unvalidated_catalog,
            axis,
            switch_penalty,
        } => genome::with_output_model(&common.out_dir, observations::MODEL, || {
            if !allow_unvalidated_catalog {
                return Err(invalid("explicit --allow-unvalidated-catalog required"));
            }
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let (metadata, checksum) = observations::load(&directory, &identity)?;
            if axis.is_some() && !metadata.full_feature_scope {
                return Err(invalid(
                    "partial feature scope cannot be threaded/reconstructed",
                ));
            }
            let (sample, sample_checksum) =
                sample::SampleIndex::load_with_checksum(&sample_path, &identity)?;
            let mut calls = observations::call(
                &metadata,
                checksum,
                &directory,
                &sample,
                sample_checksum,
                genotype::Parameters {
                    haploid_depth,
                    background,
                    max_mean_deviance,
                },
                &common.out_dir,
            )?;
            let axis: Option<threading::Axis> =
                axis.map(|path| genome::read_json(&path)).transpose()?;
            if let Some(axis) = &axis {
                observations::prepare_orientations(&metadata, &directory, &mut calls, axis)?;
            }
            genome::write_json(&common.out_dir.join("calls.json"), &calls)?;
            if let Some(axis) = axis {
                let threads = observations::thread(&metadata, &calls, axis, switch_penalty)?;
                genome::write_json(&common.out_dir.join("threads.json"), &threads)?;
            }
            Ok(
                serde_json::json!({"model": observations::MODEL, "full_feature_scope": metadata.full_feature_scope,
                    "global_factors_used": calls.global_factors_used, "excluded_features":calls.excluded_features,
                    "exclusions":calls.excluded_features_by_reason, "observations":calls.observations.as_ref().map(|p| &p.metadata_checksum)}),
            )
        }),
        Command::ThreadPartitions {
            common,
            observations: directory,
            calls,
            axis,
            switch_penalty,
        } => genome::with_output_model(&common.out_dir, observations::THREAD_MODEL, || {
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let (metadata, checksum) = observations::load(&directory, &identity)?;
            let calls_path = calls;
            let original = std::fs::read(&calls_path)?;
            let mut calls: genotype::Genotypes =
                serde_json::from_slice(&original).map_err(io::Error::other)?;
            if calls.model != observations::MODEL
                || calls
                    .observations
                    .as_ref()
                    .is_none_or(|p| p.metadata_checksum != checksum)
            {
                return Err(invalid("calls/observation linkage mismatch"));
            }
            let axis: threading::Axis = genome::read_json(&axis)?;
            observations::prepare_orientations(&metadata, &directory, &mut calls, &axis)?;
            observations::write_derived_calls(
                &common.out_dir.join("calls.json"),
                &calls_path,
                &original,
                &mut calls,
            )?;
            let threads = observations::thread(&metadata, &calls, axis, switch_penalty)?;
            genome::write_json(&common.out_dir.join("threads.json"), &threads)?;
            Ok(
                serde_json::json!({"resolved_intervals":threads.resolved_intervals,"unresolved_intervals":threads.unresolved_intervals}),
            )
        }),
        Command::Reconstruct {
            calls,
            threads,
            panel_names,
            sources,
            gap_policy,
            out_dir,
        } => genome::with_output_model(&out_dir, genome::reconstruction::MODEL, || {
            genome::reconstruction::run(
                &calls,
                &threads,
                &panel_names,
                &sources,
                &out_dir,
                gap_policy == "copy-source",
            )
        }),
        Command::EvaluateSequence {
            query,
            truth,
            paf,
            aligner_metadata,
            out_dir,
        } => genome::with_output_model(&out_dir, genome::sequence_evaluation::MODEL, || {
            genome::sequence_evaluation::run(
                &query,
                &truth,
                &paf,
                &out_dir,
                aligner_metadata.as_deref(),
            )
        }),
        Command::AlignSequence {
            query,
            truth,
            wfmash,
            threads,
            out_dir,
        } => genome::with_output_model(&out_dir, genome::sequence_evaluation::MODEL, || {
            genome::sequence_evaluation::align(&query, &truth, &wfmash, &out_dir, threads)
        }),
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
