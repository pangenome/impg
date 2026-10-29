//! Experimental inference CLI; deliberately independent of the legacy infer stitcher.
mod panel_route_search_policy;
use crate::genome_inference::{
    self as genome, calling, catalog, diplotype_chain, diplotype_dosage, genotype, joint,
    observations, panel_routes, realign_score, reanchor, sample, threading,
};
use crate::sample_mem_bwt::invalid;
use crate::syng::{SyncmerParams, SyngIndex};
use clap::{Args, Subcommand};
use std::io;
use std::io::BufRead;
use std::path::{Path, PathBuf};

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
    /// THE PANEL-SIDE PER-LOCUS PRECOMPUTE STEP (SS10.6.23 slice 1, the
    /// Frontier-2 increment's own STEP: the once-per-panel price, paid at
    /// panel-prep or standalone, NEVER per sample): precompute the row
    /// context (the deduplicated member rows' stored + canonical-scheme
    /// walks, the sequences left in the sources and keyed by their FNV-64
    /// content keys) and the per-window fold groupings into the panel
    /// asset's new section — <partition-graphs>/panel-precompute-loci/ by
    /// the default — keyed on CONTENT (the panel identity, the maps digest,
    /// the routes graph digest, the component, the frame + record origin,
    /// the stored-frame convention, k). The call auto-discovers the asset
    /// beside its --partition-graphs and serves it when every key matches
    /// (IMPG_PANEL_PRECOMPUTE=off restores the live path verbatim).
    PrecomputePanelLoci {
        /// The panel syng prefix (the panel the asset keys on)
        #[arg(long)]
        panel: String,
        /// The route graph directory (graph.json + sources)
        #[arg(long)]
        routes: PathBuf,
        /// The locality map directory (window<N>.gfa.map.json — the asset is
        /// written beside it, into <dir>/panel-precompute-loci)
        #[arg(long)]
        partition_graphs: PathBuf,
        /// The component lane — the axis path whose member rows define the
        /// window grid
        #[arg(long)]
        component: String,
        /// THE GLOBAL-PANEL FRAME (the committed global-panel call mode's
        /// frame statement: the maps translate to the global frame at load,
        /// exactly as the call's --global-frame does)
        #[arg(long, default_value_t = false)]
        global_frame: bool,
        /// THE RECORD-ORIGIN CONVENTION (the global frame's content-origin
        /// parameter, the call's own --record-origin)
        #[arg(long, default_value_t = 0, requires = "global_frame")]
        record_origin: u64,
        /// The asset output directory (default:
        /// <partition-graphs>/panel-precompute-loci)
        #[arg(long)]
        out_dir: Option<PathBuf>,
    },
    /// THE PRODUCT OBJECTIVE (the owner's product swap): per-locus diplotype calls
    /// with the model-internal quality — the realignment likelihood +
    /// territory-image extents + the alignment-induced locality domains as the
    /// selection path. Default output: calls.jsonl, one truth-free record per
    /// locus (the two called paths, the class identity, the QUAL, the E/gap
    /// fields). The called-vs-truth sequence QV lives ONLY behind --truth-qv
    /// (the test mode, the truth haplotypes as input; never in the default
    /// output).
    CallDiplotypes {
        /// The panel syng prefix
        #[arg(long)]
        panel: String,
        /// The route graph directory (graph.json + sources)
        #[arg(long)]
        routes: PathBuf,
        /// The locality map directory (window<N>.gfa.map.json)
        #[arg(long)]
        partition_graphs: PathBuf,
        /// The component lane (e.g. S288C#0#chrI) — the axis path whose
        /// member rows define the window-to-partition mapping
        #[arg(long)]
        component: String,
        /// The multi-matching census receipt (the sample's read records and
        /// placements against the panel; truth-free)
        #[arg(long)]
        census: PathBuf,
        /// The sample FASTQ (read identity + the quality-byte assertion)
        #[arg(long)]
        reads: PathBuf,
        /// The record-derivation cache (a pure function of the FASTQ and the panel)
        #[arg(long)]
        derive_cache: PathBuf,
        /// The scored loci (window ids of the component, comma-separated,
        /// distinct and sorted, e.g. 2,4,7) — or state --region instead
        /// (the two are mutually exclusive)
        #[arg(long, required_unless_present = "region", conflicts_with = "region")]
        loci: Option<String>,
        /// THE REGION-AS-INTERFACE FILTER (the owner's ruling: "detect
        /// any locus, any position, any allele"): a region chr:start-end
        /// (0-based half-open, the component's own reference frame)
        /// resolved against the component's axis rows to the covered
        /// windows — the scored loci, in place of --loci
        #[arg(long, required_unless_present = "loci", conflicts_with = "loci")]
        region: Option<String>,
        /// The selection convention (default: the alignment-induced locality
        /// domain — the current best measured selection)
        #[arg(long, value_enum, default_value_t = realign_score::Selection::Locality)]
        selection: realign_score::Selection,
        /// Resident-set guard in GiB (the 64 GiB discipline; 0 = no guard)
        #[arg(long, default_value_t = 64.0)]
        rss_budget_gib: f64,
        /// THE DIPLOTYPE SEQUENCE EMISSION (the write-up codification,
        /// the product): every called locus's reconstructed diploid
        /// sequence - the chosen pair's walks from the local graph
        /// spelled into per-locus hap1/hap2 records (the tied-class
        /// calls emit the coalesced called set per the owner's
        /// ambiguity ruling), each record carrying the model-internal
        /// QUAL. Written to calls.jsonl.sequences.jsonl beside the
        /// product calls. Truth-free product output.
        #[arg(long)]
        emit_sequences: bool,
        /// THE MASKING OPERATOR (the held-out experiment's stranger
        /// filter, folded onto the call path): candidate rows whose
        /// member path names match any of these fnmatch-style globs
        /// (e.g. 'HG002#*') leave the candidate universe - loudly
        /// counted per locus in the product record. THE AXIS ROW (the
        /// --component lane) is never excluded: it is the coordinate
        /// frame. Default (empty): no row is excluded, the committed
        /// candidate universe verbatim - the default product path is
        /// byte-identical to the committed receipts.
        #[arg(long, num_args = 1.., value_delimiter = ',')]
        exclude_paths: Vec<String>,
        /// THE FLANK-OVERLAP CANDIDATE FILTER (the owner's design, the
        /// presence semantics): N in bp (e.g. 10000) or a percent of the
        /// window length (e.g. 10% -> ceil(pct x length)). A member row is
        /// kept only if it carries aligned material over BOTH FLANK
        /// intervals of the window's axis extent - the first and last
        /// N bp, in the focal reference's own frame; a window shorter
        /// than N degenerates to simple presence; the axis row never
        /// excluded; the excluded rows loudly counted per locus.
        /// Default: OFF (the committed candidate universe verbatim -
        /// the default product path is byte-identical to the committed
        /// receipts).
        #[arg(long)]
        flank_overlap: Option<realign_score::FlankOverlap>,
        /// THE GLOBAL-PANEL FRAME (the campaign's lever 1, the
        /// global-panel call mode:
        /// docs/syng-gmem-bwt/global-panel-call-mode.md). When set,
        /// --panel names the ONE global panel and --partition-graphs
        /// the territory-local substrate (member rows "path:start-end",
        /// slice-relative coordinates): the maps are translated to
        /// the global frame at load, the machinery runs against the
        /// global grid, and the product emission is translated back
        /// to the partition's own frame (the committed
        /// reference-anchored output shape). Default (unset): the
        /// per-territory frame verbatim — the committed path
        /// byte-for-byte.
        #[arg(long, default_value_t = false)]
        global_frame: bool,
        /// THE RECORD-ORIGIN CONVENTION (the global-panel call mode's
        /// content-origin parameter): the committed substrate's
        /// locus.fa content need not sit at its header's own interval
        /// — the L1-L4 substrates' content sits ONE BASE LATER (d=1,
        /// the uniform delta=1 AGC-identity receipts), the pilot's at
        /// the header interval exactly (d=0). The global frame, which
        /// deletes locus.fa by design, must be TOLD the convention:
        /// the maps and the component lane translate at the content
        /// origin start+d, and the emission's back-translation
        /// cancels it exactly (the product output is
        /// origin-independent). Default 0. Requires --global-frame
        /// (the committed local frame is closed — it derives the
        /// convention itself, per record, against the AGC).
        #[arg(long, default_value_t = 0, requires = "global_frame")]
        record_origin: u64,
        /// THE FAR-BIN PRUNE (the owner's own routing principle, the
        /// §10.6.17-named census lever, enforced at the call's census
        /// view: "we used to decide where to put each read, so it's
        /// basically being routed by its best MEMs. If they go into
        /// different bins unless they're basically adjacent, we
        /// don't want to be thinking about them."). Census records
        /// whose placements scatter into NON-ADJACENT windows (the
        /// measured far-bin repeat class — L1's one 43,686-mult
        /// single-anchor pattern at windows [12,16], 31.1% of the
        /// restricted census mass) leave the census view entirely;
        /// boundary-straddlers keep BOTH adjacent bins (pure interval
        /// overlap, the routing's own rule). The adjacency test is
        /// geometry only (the touched-window diameter vs the census's
        /// own pair-bins constant 1) — no tuned thresholds.
        /// DEFAULT: ON at the global frame (the owner's ruling on the
        /// §10.6.18 gate receipts: zero truth losses at L1-L4, the
        /// pilot's FOUR truth WINS — rank-1 4/22 -> 8/22, median QV
        /// 36.3 -> 45.0 — L1 loci 74.8 -> 34.7s). The opt-out form
        /// stays for the record: `--far-bin-prune=false` (the
        /// committed unpruned census view, byte-identical to the
        /// pre-ruling receipts). The LOCAL frame is byte-frozen and
        /// refuses the prune loudly — state `--far-bin-prune=false`
        /// for the committed local form.
        #[arg(
            long,
            value_parser = clap::value_parser!(bool),
            default_value_t = true,
            default_missing_value = "true",
            num_args = 0..=1
        )]
        far_bin_prune: bool,
        /// THE TEST MODE (the owner's ruling): the truth haplotypes as panel
        /// path names, exactly two (e.g. S288C#0#chrMT SK1#0#chrMT) —
        /// computes the called-vs-truth sequence QV per locus (the committed
        /// convention) into calls.jsonl.truth-qv.jsonl; never mixed into the
        /// default product output
        #[arg(long, num_args = 2)]
        truth_qv: Option<Vec<String>>,
        /// The pinned biWFA alignment helper for the test mode (gap-affine
        /// End2End 0/4/6/2, Medium memory, no heuristic — the same machinery
        /// the validated QV receipts used)
        #[arg(long)]
        qv_biwfa: Option<PathBuf>,
        /// New output directory (refused if it exists): calls.jsonl + the
        /// instrument receipts beside it
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// THE CHAIN/PHASING LAYER (the owner's product ruling: the diploid
    /// product is TWO whole-chromosome molecules per component): the
    /// thin orientation DP over the CLI's CLOSED per-locus calls
    /// (calls.jsonl) — the junction-spanning reads' votes (the edge
    /// layer, the census's crossing records) plus the called folds'
    /// shared-segment continuation (the adjacency layer) orient each
    /// boundary; the accumulated chain emits the two molecules (the
    /// per-locus orientation, the spelled sequences, the per-boundary
    /// link confidence). The chain NEVER writes the calls — the
    /// per-locus product is bit-identical pre/post chain (the thin-
    /// layer property, proven by the emitted calls fingerprint).
    /// The switch evaluation against the truth haplotypes lives ONLY
    /// behind --truth-qv-file (the test mode, the product run's own
    /// truth-QV artifact as input; never in the default output).
    ChainMolecules {
        /// The CLI product's closed per-locus calls (calls.jsonl)
        #[arg(long)]
        calls: PathBuf,
        /// The panel syng prefix
        #[arg(long)]
        panel: String,
        /// The route graph directory (graph.json + sources)
        #[arg(long)]
        routes: PathBuf,
        /// The locality map directory (window<N>.gfa.map.json)
        #[arg(long)]
        partition_graphs: PathBuf,
        /// The multi-matching census receipt (the same truth-free
        /// input the product run consumed)
        #[arg(long)]
        census: PathBuf,
        /// THE TEST MODE (the owner's ruling, assessment-side only):
        /// the product run's calls.jsonl.truth-qv.jsonl — consumed
        /// read-only into molecules.jsonl.truth-qv.jsonl (the switch
        /// table and the per-locus accuracy, separate columns, never
        /// conflated)
        #[arg(long)]
        truth_qv_file: Option<PathBuf>,
        /// Resident-set guard in GiB (the 64 GiB discipline; 0 = no guard)
        #[arg(long, default_value_t = 64.0)]
        rss_budget_gib: f64,
        /// New output directory (refused if it exists): molecules.jsonl
        /// + the test-mode artifact beside it
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// THE DOSAGE SURFACE (the copy-number layer of the two-molecule
    /// emission): the per-panel-segment copy count implied by the
    /// chain's two molecules — derived from molecules.jsonl and
    /// nothing else product-side, with REPEAT-VISIT MULTIPLICITY
    /// (a molecule traversing a segment twice carries its true
    /// multiplicity), the territory/axis extents per segment, and
    /// the DEPTH-CONSISTENCY QC (the observed per-segment census
    /// mass vs the expected under the emitted dosage — the
    /// likelihood's own per-covered-copy expectation, the ratio
    /// stated per segment, never binned, no thresholds). The
    /// truth-referenced dosage comparison lives ONLY behind
    /// --truth-qv-file (the test mode, the product run's own
    /// truth-QV artifact + the receipt's fold identities as input;
    /// never in the default output).
    EmitDosage {
        /// The chain layer's two-molecule product (molecules.jsonl)
        #[arg(long)]
        molecules: PathBuf,
        /// The panel syng prefix
        #[arg(long)]
        panel: String,
        /// The locality map directory (window<N>.gfa.map.json; the
        /// axis window extents)
        #[arg(long)]
        partition_graphs: PathBuf,
        /// The multi-matching census receipt (the same truth-free
        /// input the product and chain runs consumed; the observed
        /// mass)
        #[arg(long)]
        census: PathBuf,
        /// THE TEST MODE (the owner's ruling, assessment-side only):
        /// the product run's calls.jsonl.truth-qv.jsonl — the truth
        /// pair per locus, consumed read-only into
        /// dosage.jsonl.truth-qv.jsonl (the agreement table and the
        /// per-locus accuracy, separate columns, never conflated)
        #[arg(long, requires = "folds")]
        truth_qv_file: Option<PathBuf>,
        /// The product run's instrument receipt (the per-locus fold
        /// identities — the truth folds' member rows; required by
        /// the test mode)
        #[arg(long, requires = "truth_qv_file")]
        folds: Option<PathBuf>,
        /// Resident-set guard in GiB (the 64 GiB discipline; 0 = no guard)
        #[arg(long, default_value_t = 64.0)]
        rss_budget_gib: f64,
        /// New output directory (refused if it exists): dosage.jsonl
        /// + the test-mode artifact beside it
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// THE K-REANCHOR TRANSFORM (the owner's dynamic local-graph design,
    /// stated verbatim: "dynamically, take the sequences and build a local
    /// syng graph using different kmer sizes and project the mappings we
    /// have from the original set onto those"): the committed multi-matching
    /// census records re-anchored onto a target panel's node grid by pure
    /// COORDINATE LOOKUP — no census re-run, no re-alignment. Per merged
    /// read-matched interval [lo, hi): the target panel's walk steps with
    /// bp >= lo and bp + k_target <= hi (a mapping at sequence offset X
    /// lands on whichever k node contains X); the `partitions` window-touch
    /// field re-anchors by the routing's own territory-touch rule onto the
    /// target tiling's beds. The occurrence structure (path, start,
    /// orientation, multiplicity, t_r, share, anchors, total_occurrences,
    /// pair_bins) is preserved verbatim; the empty-interval read-edge-anchor
    /// class keeps its window attribution verbatim (the named zero-coverage
    /// boundary of the pure coordinate lookup). The walks are loaded in
    /// memory from the two syng indexes (no dump round trip); the census is
    /// streamed; the transform is seconds-class for millions of occurrences.
    ReanchorCensus {
        /// The committed multi-matching census receipt (jsonl or jsonl.gz)
        #[arg(long)]
        census: PathBuf,
        /// The source panel syng prefix (the grid the census records were
        /// derived against; supplies the coordinate recovery)
        #[arg(long)]
        source_panel: String,
        /// The target panel syng prefix (the new grid)
        #[arg(long)]
        panel: String,
        /// The target tiling's partition bed directory (partitionN.bed)
        #[arg(long)]
        bed_dir: PathBuf,
        /// The target locality map directory (windowN.gfa.map.json)
        #[arg(long)]
        partition_graphs: PathBuf,
        /// The component lane (the axis path whose member rows define the
        /// window-to-partition mapping)
        #[arg(long)]
        component: String,
        /// The re-anchored census output (jsonl or jsonl.gz)
        #[arg(long)]
        out: PathBuf,
        /// Optional stats receipt path (the transform's measured counters)
        #[arg(long)]
        stats: Option<PathBuf>,
    },
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
    /// The call-path route inventory by lookup: graph.json's lanes from the panel's own name map, the census routing cross-check; no native ports/profiles (the calling path never reads them)
    DeriveRoutes {
        #[command(flatten)]
        common: Common,
        /// The re-anchored census receipt (jsonl or jsonl.gz) — every
        /// occurrence's path id is cross-checked against the target
        /// panel's name map (the routing bind)
        #[arg(long)]
        census: PathBuf,
        #[arg(long, required=true, num_args=1..)]
        sources: Vec<String>,
        #[arg(long, required=true, value_delimiter=',', num_args=1..)]
        read_lengths: Vec<u64>,
        #[arg(long, default_value_t = 65536)]
        core_bp: u64,
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
    /// Fair sample-ranked lazy coupled search; explicit immutable-parent resume (no emission)
    SearchPanelRoutesGuided {
        #[command(flatten)]
        common: Common,
        #[command(flatten)]
        evidence: RouteEvidence,
        #[command(flatten)]
        options: panel_route_search_policy::Options,
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
    /// THE THREE-STEP UI - PARTITION (the owner's UI doctrine: map -> partition -> call):
    /// the PANEL-side step, built once per panel territory and REUSED by every sample.
    /// Emits the partition's row into the ONE partition artifact - the four-column BED
    /// (chrom, start, end, PARTITION_NAME) - and builds the panel-side precomputes keyed
    /// by the partition name under <out-dir>/<name>/: the local syng (the committed
    /// recipe: k=63, s=8, seed 7, pos64), the window beds (regenerated verbatim from the
    /// committed locality maps - the region grid), and the call-path routes inventory
    /// (derive-routes, the lookup lever). The member subset per partition is DERIVED from
    /// the reference range (the --sources fasta IS the range's member subset); the
    /// locality maps are the durable panel-side precompute and are REGISTERED, never
    /// rebuilt here.
    Partition {
        /// The partition name (the BED4 fourth column; keys every panel-side precompute)
        #[arg(long)]
        name: String,
        /// The reference range chromosome (the BED4 first column)
        #[arg(long)]
        chrom: String,
        /// The reference range start, 0-based half-open (the BED4 second column)
        #[arg(long)]
        start: u64,
        /// The reference range end (the BED4 third column)
        #[arg(long)]
        end: u64,
        /// The member-subset fasta for the range (the panel sequences overlapping it;
        /// the substrate of the local syng)
        #[arg(long)]
        sources: PathBuf,
        /// The component lane - the reference axis path (e.g. HG002#1#chr1:start-end)
        #[arg(long)]
        component: String,
        /// The reference axis JSON (the committed panel-side asset)
        #[arg(long)]
        axis: PathBuf,
        /// The locality map directory (windowN.gfa.map.json - the committed window grid;
        /// the window beds are regenerated from it verbatim)
        #[arg(long)]
        partition_graphs: PathBuf,
        /// The census receipt for the routes cross-check (any sample's census against
        /// this panel; the routing bind)
        #[arg(long)]
        census: PathBuf,
        /// The read lengths for the routes inventory
        #[arg(long, value_delimiter = ',', num_args = 1.., default_values_t = vec![148u64])]
        read_lengths: Vec<u64>,
        /// THE GLOBAL-PANEL CALL MODE (the campaign's lever 1:
        /// docs/syng-gmem-bwt/global-panel-call-mode.md): the ONE
        /// global panel syng prefix. When set (with --global-routes),
        /// the per-territory local syng build and the per-territory
        /// derive-routes are SKIPPED (the cut and the call bind to
        /// the global panel; prep's panel-side pole collapses to the
        /// window-bed regeneration) and the registry records the
        /// global pair with "frame": "global". Default (unset): the
        /// committed per-territory precomputes verbatim.
        #[arg(long, requires = "global_routes")]
        global_panel: Option<String>,
        /// The whole-panel routes inventory (derive-routes --panel
        /// <global syng> --sources <panel AGC>, built once,
        /// panel-side) — the call reads only its lane inventory.
        #[arg(long, requires = "global_panel")]
        global_routes: Option<PathBuf>,
        /// THE GLOBAL PANEL'S AGC (the fleet call path's own input,
        /// RECORDED in the registry so the three-step CALL step needs
        /// no manual flag for it: the fleet children build their AGC
        /// indexes from the panel's own AGC file). Requires
        /// --global-panel.
        #[arg(long, requires = "global_panel")]
        global_agc: Option<PathBuf>,
        /// THE RECORD-ORIGIN CONVENTION (the global-panel call mode's
        /// content-origin parameter, recorded in the registry so the
        /// CALL step states the same convention): the committed
        /// substrate's locus.fa content need not sit at its header's
        /// own interval — the L1-L4 substrates' content sits ONE BASE
        /// LATER (d=1, the uniform delta=1 AGC-identity receipts),
        /// the pilot's at the header interval exactly (d=0). The
        /// global frame, which deletes locus.fa by design, must be
        /// TOLD the convention; the committed local frame derives it
        /// per record against the AGC and needs no statement.
        /// Default 0. Requires --global-panel.
        #[arg(long, default_value_t = 0, requires = "global_panel")]
        record_origin: u64,
        /// The partitions root output directory (partitions.bed + <name>/)
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// THE THREE-STEP UI - PANEL-PREP (the owner's 2-3 step doctrine:
    /// "make it like 2-3 steps to do whole genome analysis" — STEP 1,
    /// once per panel, reusable forever): the ONE panel-prep step that
    /// does BOTH panel-side pieces in sequence — (1) the whole-panel
    /// routes inventory (derive-routes over the ONE global panel; the
    /// ~25-min parallel derive, built at the FIRST invocation and
    /// reused by every later one — an existing inventory is skipped,
    /// loudly), and (2) the PARTITION step's global-panel form (the
    /// BED4 row + the window beds + the registry that records the
    /// global pair, the content origin, and the panel AGC so the
    /// CALL step needs no manual flags). One invocation per partition
    /// row (the routes inventory builds once, at the first).
    PanelPrep {
        /// The partition name (the BED4 fourth column; keys every panel-side precompute)
        #[arg(long)]
        name: String,
        /// The reference range chromosome (the BED4 first column)
        #[arg(long)]
        chrom: String,
        /// The reference range start, 0-based half-open (the BED4 second column)
        #[arg(long)]
        start: u64,
        /// The reference range end (the BED4 third column)
        #[arg(long)]
        end: u64,
        /// The member-subset fasta for the range (the panel sequences overlapping it)
        #[arg(long)]
        sources: PathBuf,
        /// The component lane - the reference axis path (e.g. HG002#1#chr1:start-end)
        #[arg(long)]
        component: String,
        /// The reference axis JSON (the committed panel-side asset)
        #[arg(long)]
        axis: PathBuf,
        /// The locality map directory (windowN.gfa.map.json - the committed window grid;
        /// the window beds are regenerated from it verbatim)
        #[arg(long)]
        partition_graphs: PathBuf,
        /// The census receipt for the routes cross-check (any sample's census against
        /// this panel; the routing bind)
        #[arg(long)]
        census: PathBuf,
        /// The read lengths for the routes inventory
        #[arg(long, value_delimiter = ',', num_args = 1.., default_values_t = vec![148u64])]
        read_lengths: Vec<u64>,
        /// THE GLOBAL PANEL syng prefix (the ONE panel; the routes
        /// inventory and the registry both name it)
        #[arg(long)]
        global_panel: String,
        /// THE GLOBAL PANEL'S AGC (derive-routes' --sources input;
        /// recorded in the registry for the CALL step's fleet children)
        #[arg(long)]
        global_agc: PathBuf,
        /// The whole-panel routes inventory output directory (built at
        /// the first invocation, reused forever after)
        #[arg(long)]
        global_routes: PathBuf,
        /// THE RECORD-ORIGIN CONVENTION (the content-origin parameter,
        /// recorded in the registry so the CALL step states the same
        /// convention; the committed L1-L4 substrates sit at d=1, the
        /// pilot at d=0). Default 0.
        #[arg(long, default_value_t = 0)]
        record_origin: u64,
        /// The partitions root output directory (partitions.bed + <name>/)
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// THE BATCHED PANEL-PREP (one global panel load, N territories) —
    /// the `cut-batch-seq` pattern at the panel-prep step.
    /// WHERE THE LOAD COSTS ANYTHING: in STEADY STATE it does not — an
    /// already-prepped territory takes the per-locus reuse path in both
    /// forms (~0.02s, nothing to batch). The multiplication is the
    /// ONE-TIME BUILD of a panel's territory backlog: every
    /// substrate-done-but-unprepped territory pays `SyngIndex::load` of
    /// the same panel in its own process. Measured at
    /// RAYON_NUM_THREADS=8 on three genuinely un-prepped territories
    /// (T00643, T00645, T00648): the single form 271.25 / 271.77 /
    /// 258.35s (801.37s total), this form 372.81s — one 189.4s load
    /// plus 58.0 / 56.2 / 56.8s of per-territory work. The load is 71%
    /// of a single-territory wall and is paid once per batch instead of
    /// once per territory; per-territory cost becomes
    /// load/B + ~57s for a batch of B. This form loads the panel ONCE
    /// and then runs the committed single-territory panel-prep row for
    /// every manifest line, so each territory's outputs are
    /// byte-identical by construction (same function, same inputs):
    /// verified on those three territories — partition.json identical
    /// but for its out-dir spelling, all 74 window beds and all 12
    /// precompute-asset files identical by sha256. Like the single form
    /// it REUSES a whole-panel routes inventory and a current-format
    /// per-locus precompute asset rather than rebuilding either.
    PanelPrepBatch {
        /// The JSONL manifest: one territory per line (name, chrom,
        /// start, end, sources, component, axis, partition_graphs,
        /// census, read_lengths, record_origin)
        #[arg(long)]
        manifest: PathBuf,
        /// THE GLOBAL PANEL syng prefix (the ONE panel: loaded once for
        /// the whole manifest instead of once per row)
        #[arg(long)]
        global_panel: String,
        /// THE GLOBAL PANEL'S AGC (derive-routes' --sources input)
        #[arg(long)]
        global_agc: PathBuf,
        /// The whole-panel routes inventory (built by the first row that
        /// needs it, reused by every later one)
        #[arg(long)]
        global_routes: PathBuf,
        /// The partitions root output directory (partitions.bed + <name>/):
        /// the rows share it — one BED4 universe, one immutability check
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// THE THREE-STEP UI - MAP (the owner's UI doctrine: map -> partition -> call): the
    /// SAMPLE-side step, run ONCE per sample: the projection. Cuts every (selected)
    /// partition's territory evidence out of the sample's projection - the per-partition
    /// census, record-derivation cache and reads FASTA the CALL step consumes. The
    /// projection itself (the shard directory) is the sample's map at the whole-sample
    /// level and is built once per sample (sample_projection project + dedup); pass its
    /// directory to reuse a committed projection. --local-derive selects the corrected
    /// cut recipe (route globally, attribute locally) for territories that require it.
    Map {
        /// The GLOBAL panel syng prefix (the projection's panel)
        #[arg(long)]
        panel: String,
        /// The global panel AGC (the projection's --sequence-files input)
        #[arg(long)]
        agc: PathBuf,
        /// The sample's projection shard directory (the sample's map, built once:
        /// sample_projection project + dedup)
        #[arg(long)]
        shards_dir: PathBuf,
        /// The partitions root (the PARTITION step's out-dir: partitions.bed + <name>/)
        #[arg(long)]
        partitions: PathBuf,
        /// The sample FASTQ (recorded in the map manifest; the CALL step's reads input)
        #[arg(long)]
        reads: PathBuf,
        /// The region filter: a comma-list of partition names, or chr:start-end
        /// (the BED rows overlapping the range)
        #[arg(long)]
        region: Option<String>,
        /// The corrected cut recipe (route globally, attribute locally - the
        /// pilot-class territory); default: the committed original recipe
        #[arg(long, default_value_t = false)]
        local_derive: bool,
        /// New output directory (the map root: <name>/ per partition + manifest.json)
        #[arg(long)]
        out_dir: PathBuf,
    },
    /// THE THREE-STEP UI - CALL (the owner's UI doctrine: map -> partition -> call): the
    /// per-sample call over the partitions: consumes the MAP step's cuts + the PARTITION
    /// step's panel-side precomputes and calls every (selected) partition - the per-locus
    /// diplotype calls with the REFERENCE-ANCHORED outputs (the per-locus hap1/hap2
    /// diplotype sequences - the emission - the model-internal QUAL, the dosage fields).
    /// --region selects on demand: partition names, or any chr:start-end range resolved
    /// against the partitions' own reference frames ("detect any locus, any position,
    /// any allele"); --exclude-paths is the masking operator; --truth-qv the test mode.
    /// Named call-region because the literal `call` name is held by the legacy
    /// provisional presence-compatibility command (the committed CLI paths are immutable).
    /// That collision reaches INSIDE the `infer` subtree too - `impg infer call`
    /// already resolves to the legacy provisional command (verified: adding
    /// `name = "call"` here made clap route BOTH `call` and `call-region` to the
    /// legacy variant, i.e. it would silently run the wrong command). So this
    /// stays `call-region` until the legacy provisional command is renamed, which
    /// is a schema decision, not a rename.
    CallRegion {
        /// The MAP step's output directory (the per-partition cuts + the manifest)
        #[arg(long)]
        map: PathBuf,
        /// The partitions root (the PARTITION step's out-dir: partitions.bed + <name>/)
        #[arg(long)]
        partitions: PathBuf,
        /// The region filter: a comma-list of partition names, or chr:start-end (the
        /// covered windows of every overlapping partition)
        #[arg(long)]
        region: Option<String>,
        /// THE FAR-BIN PRUNE (the product default, the §10.6.18 ruling
        /// executed at the call-diplotypes engine): applies at the
        /// global-frame call path (the registry's frame:global
        /// partitions and the --global-panel form); the committed
        /// local-frame path is byte-frozen and never prunes. DEFAULT:
        /// ON at the global frame; the opt-out form stays for the
        /// record: `--far-bin-prune=false`.
        #[arg(
            long,
            value_parser = clap::value_parser!(bool),
            default_value_t = true,
            default_missing_value = "true",
            num_args = 0..=1
        )]
        far_bin_prune: bool,
        /// THE MASKING OPERATOR (the held-out experiment's stranger filter): candidate
        /// rows whose member path names match any fnmatch-style glob leave the candidate
        /// universe, loudly counted (the axis row never excluded). Default: none.
        #[arg(long, num_args = 1.., value_delimiter = ',')]
        exclude_paths: Vec<String>,
        /// THE FLANK-OVERLAP CANDIDATE FILTER (the owner's design, the presence
        /// semantics): N in bp (e.g. 10000) or a percent of the window length
        /// (e.g. 10% -> ceil(pct x length)). A member row is kept only if it
        /// carries aligned material over BOTH FLANK intervals of the window's
        /// axis extent - the first and last N bp, in the focal reference's own
        /// frame; a window shorter than N degenerates to simple presence; the
        /// axis row never excluded; the excluded rows loudly counted per
        /// partition. Default: OFF (the committed candidate universe verbatim).
        #[arg(long)]
        flank_overlap: Option<realign_score::FlankOverlap>,
        /// THE TEST MODE (assessment-side only): the truth haplotypes as panel path
        /// names, exactly two - the called-vs-truth sequence QV into
        /// calls.jsonl.truth-qv.jsonl per partition
        #[arg(long, num_args = 2)]
        truth_qv: Option<Vec<String>>,
        /// The pinned biWFA alignment helper for the --truth-qv test mode
        #[arg(long)]
        qv_biwfa: Option<PathBuf>,
        /// The selection convention (default: the alignment-induced locality domain)
        #[arg(long, value_enum, default_value_t = realign_score::Selection::Locality)]
        selection: realign_score::Selection,
        /// THE GLOBAL-PANEL CALL MODE (the campaign's lever 1:
        /// docs/syng-gmem-bwt/global-panel-call-mode.md): the ONE
        /// global panel syng prefix. When set (with --global-routes),
        /// every selected partition is called against the global
        /// panel + the ONE whole-panel routes inventory instead of
        /// the registry's per-territory syng/routes — the maps are
        /// translated to the global frame at load and the product
        /// emission back to the partition's own frame. Default
        /// (unset): the committed per-territory path byte-for-byte.
        #[arg(long, requires = "global_routes")]
        global_panel: Option<String>,
        /// The whole-panel routes inventory (derive-routes --panel
        /// <global syng> --sources <panel AGC>, built once,
        /// panel-side) — the call reads only its lane inventory.
        #[arg(long, requires = "global_panel")]
        global_routes: Option<PathBuf>,
        /// Resident-set guard in GiB (the 64 GiB discipline; 0 = no guard)
        #[arg(long, default_value_t = 64.0)]
        rss_budget_gib: f64,
        /// New output directory (the calls root: <name>/ per partition)
        #[arg(long)]
        out_dir: PathBuf,
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
        &str,
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
        let result = operation(&mut evaluator, &common.out_dir, &sample_checksum)?;
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
        Command::PrecomputePanelLoci {
            panel,
            routes,
            partition_graphs,
            component,
            global_frame,
            record_origin,
            out_dir,
        } => {
            // (SS10.6.23 slice 1 + SS10.6.25 slice 2 — the once-per-panel
            // precompute step. An existing asset of the CURRENT format is
            // REUSED loudly, never rebuilt per invocation (the
            // routes-inventory precedent: a panel-side asset is built once;
            // remove it to rebuild). A format-/1 asset (slice 1's, without
            // the bind-scaffolding section) is UPGRADED loudly — rebuilt
            // in place so the bind scaffolding is keyed and served.
            let out_dir = out_dir
                .unwrap_or_else(|| partition_graphs.join("panel-precompute-loci"));
            if out_dir.join("manifest.json").is_file() {
                let manifest_format = std::fs::read_to_string(out_dir.join("manifest.json"))
                    .ok()
                    .and_then(|text| {
                        serde_json::from_str::<genome::realign_score::PanelPrecomputeManifest>(&text)
                            .ok()
                    })
                    .map(|manifest| manifest.format_text());
                if manifest_format.as_deref() == Some("impg-panel-precompute-loci/2") {
                    println!(
                        "=== panel-precompute: the per-locus precompute asset {} exists — REUSED \
                         (built once per panel; remove it to rebuild)",
                        out_dir.display()
                    );
                    return Ok(());
                }
                println!(
                    "=== panel-precompute: the existing asset {} is format {:?} (slice 1's) — \
                     REBUILDING with the bind-scaffolding section (the once-per-panel price)",
                    out_dir.display(),
                    manifest_format.unwrap_or_else(|| "unreadable".to_string())
                );
            } else {
                println!(
                    "=== panel-precompute: building the per-locus precompute asset into {} \
                     (the row context + the fold layouts + the bind scaffolding; the \
                     once-per-panel price)",
                    out_dir.display()
                );
            }
            genome::realign_score::precompute_panel_loci(
                &panel,
                &routes,
                &partition_graphs,
                &component,
                global_frame,
                record_origin,
                &out_dir,
                None,
            )
        }
        Command::CallDiplotypes {
            panel,
            routes,
            partition_graphs,
            component,
            census,
            reads,
            derive_cache,
            loci,
            region,
            selection,
            rss_budget_gib,
            truth_qv,
            qv_biwfa,
            emit_sequences,
            exclude_paths,
            flank_overlap,
            global_frame,
            record_origin,
            far_bin_prune,
            out_dir,
        } => {
            // (The output manifest's model label — the class-scan
            // order states itself here too (the owner's Ruling 3):
            // a sparse-skip run's manifest can never be mistaken for
            // an exact-mode one. THE FAR-BIN PRUNE likewise labels
            // its own runs. Env-unset keeps the committed string
            // verbatim.)
            let model_label: std::borrow::Cow<str> =
                match std::env::var("IMPG_CLASS_SCAN_ORDER") {
                    Ok(value) if value == "sparse-skip" => {
                        std::borrow::Cow::Owned(format!(
                            "{}-classscan-sparse",
                            realign_score::MODEL,
                        ))
                    }
                    _ => std::borrow::Cow::Borrowed(realign_score::MODEL),
                };
            let model_label: std::borrow::Cow<str> = if far_bin_prune {
                std::borrow::Cow::Owned(format!("{model_label}-farbin-prune"))
            } else {
                model_label
            };
            genome::with_output_model(
            &out_dir,
            &model_label,
            || {
                let options = realign_score::Options {
                    panel,
                    routes,
                    partition_graphs,
                    component,
                    census,
                    reads,
                    derive_cache,
                    loci: loci.unwrap_or_default(),
                    // The instrument receipts land beside the product
                    // calls under the CLI's own output directory (the
                    // translation-gate and checker artifacts).
                    out: out_dir.join("instrument-receipt.jsonl"),
                    exactness_out: out_dir.join("instrument-receipt.exactness.jsonl"),
                    exactness_sample: 100,
                    // The committed runner's artifact set verbatim: the
                    // skeleton sidecar is part of the receipt-of-record
                    // layout the checker consumes.
                    skeleton_out: Some(out_dir.join("instrument-receipt.skeleton.jsonl")),
                    anatomy_out: None,
                    frame_audit_only: false,
                    rss_budget_gib,
                    selection,
                    product_calls: Some(out_dir.join("calls.jsonl")),
                    truth_qv: truth_qv.unwrap_or_default(),
                    product_sequences: emit_sequences,
                    exclude_paths,
                    flank_overlap,
                    region,
                    global_frame,
                    record_origin,
                    far_bin_prune,
                    qv_biwfa: qv_biwfa.unwrap_or_else(|| {
                        "genome/instrumented/qv-biwfa/target/release/qv-biwfa".into()
                    }),
                };
                let truth_qv_enabled = !options.truth_qv.is_empty();
                realign_score::run(options)?;
                let calls = out_dir.join("calls.jsonl");
                let truth_qv_path = out_dir.join("calls.jsonl.truth-qv.jsonl");
                let sequences_path = out_dir.join("calls.jsonl.sequences.jsonl");
                Ok(serde_json::json!({
                    "calls": genome::reconstruction::fingerprint(&calls)?,
                    "calls_bytes": std::fs::metadata(&calls)?.len(),
                    "instrument_receipt": genome::reconstruction::fingerprint(
                        &out_dir.join("instrument-receipt.jsonl"),
                    )?,
                    "sequences_emission": if emit_sequences {
                        genome::reconstruction::fingerprint(&sequences_path)?
                    } else {
                        serde_json::Value::Null
                    },
                    "truth_qv_test_mode": if truth_qv_enabled {
                        genome::reconstruction::fingerprint(&truth_qv_path)?
                    } else {
                        serde_json::Value::Null
                    },
                    "selection": format!("{selection:?}"),
                    "scope": "per-locus-diplotype-calls-with-model-internal-qual",
                    "truth_in_default_output": false,
                }))
            },
            )
        },
        Command::ChainMolecules {
            calls,
            panel,
            routes,
            partition_graphs,
            census,
            truth_qv_file,
            rss_budget_gib,
            out_dir,
        } => genome::with_output_model(&out_dir, diplotype_chain::MODEL, || {
            let options = diplotype_chain::Options {
                calls: calls.clone(),
                panel,
                routes,
                partition_graphs,
                census,
                truth_qv: truth_qv_file.clone(),
                rss_budget_gib,
                out_dir: out_dir.clone(),
            };
            diplotype_chain::run(options)?;
            let truth_out = out_dir.join("molecules.jsonl.truth-qv.jsonl");
            Ok(serde_json::json!({
                "molecules": genome::reconstruction::fingerprint(
                    &out_dir.join("molecules.jsonl"),
                )?,
                "calls_consumed": genome::reconstruction::fingerprint(&calls)?,
                "calls_bytes": std::fs::metadata(&calls)?.len(),
                "truth_qv_test_mode": if truth_qv_file.is_some() {
                    genome::reconstruction::fingerprint(&truth_out)?
                } else {
                    serde_json::Value::Null
                },
                "scope": "two-molecules-per-component-chain-over-closed-calls",
                "truth_in_default_output": false,
            }))
        }),
        Command::EmitDosage {
            molecules,
            panel,
            partition_graphs,
            census,
            truth_qv_file,
            folds,
            rss_budget_gib,
            out_dir,
        } => genome::with_output_model(&out_dir, diplotype_dosage::MODEL, || {
            let options = diplotype_dosage::Options {
                molecules: molecules.clone(),
                panel,
                partition_graphs,
                census,
                truth_qv_file: truth_qv_file.clone(),
                folds: folds.clone(),
                rss_budget_gib,
                out_dir: out_dir.clone(),
            };
            diplotype_dosage::run(options)?;
            let truth_out = out_dir.join("dosage.jsonl.truth-qv.jsonl");
            Ok(serde_json::json!({
                "dosage": genome::reconstruction::fingerprint(
                    &out_dir.join("dosage.jsonl"),
                )?,
                "molecules_consumed": genome::reconstruction::fingerprint(&molecules)?,
                "truth_qv_test_mode": if truth_qv_file.is_some() {
                    genome::reconstruction::fingerprint(&truth_out)?
                } else {
                    serde_json::Value::Null
                },
                "scope": "per-panel-segment-copy-counts-derived-from-the-two-molecules",
                "truth_in_default_output": false,
            }))
        }),
        Command::ReanchorCensus {
            census,
            source_panel,
            panel,
            bed_dir,
            partition_graphs,
            component,
            out,
            stats,
        } => reanchor::run(&reanchor::Options {
            census,
            source_panel,
            panel,
            bed_dir,
            partition_graphs,
            component,
            out,
            stats,
        }),
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
        Command::DeriveRoutes {
            common,
            census,
            sources,
            read_lengths,
            core_bp,
        } => genome::with_output_model(&common.out_dir, panel_routes::MODEL, || {
            let identity = genome::PanelIdentity::read(&common.panel)?;
            let panel = SyngIndex::load(&common.panel, SyncmerParams::default())?;
            let mut receipt = panel_routes::lookup(
                &panel,
                identity,
                &census,
                &sources,
                read_lengths,
                &common.out_dir,
                core_bp,
            )?;
            // THE LEGACY-PANEL OPT-IN RECEIPT (the route graph's own
            // manifest): the panel's metadata version and the opt-in
            // state are part of the artifact's provenance — a v1 panel
            // is admissible ONLY under IMPG_PANEL_LEGACY_V1_OK=1 for
            // committed-evidence continuity (the graph's identity
            // checksums bind it to the exact sidecars either way).
            if let Some(value) = receipt.as_object_mut() {
                value.insert(
                    "panel_metadata_version".into(),
                    serde_json::json!(genome::panel_metadata_version(&common.panel)?),
                );
                value.insert(
                    "panel_legacy_opt_in".into(),
                    serde_json::json!(
                        std::env::var("IMPG_PANEL_LEGACY_V1_OK").ok().as_deref() == Some("1")
                    ),
                );
            }
            Ok(receipt)
        }),
        Command::EvaluatePanelRoutes {
            common,
            evidence,
            assignment,
            native_family,
        } => route_run(common, evidence, |e, _out, _sample| {
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
        } => route_run(common, evidence, |e, out, _sample| {
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
        Command::SearchPanelRoutesGuided {
            common,
            evidence,
            options,
        } => {
            let cache_terms = evidence.cache_terms;
            route_run(common, evidence, |e, out, sample_checksum| {
                panel_route_search_policy::run(e, sample_checksum.into(), cache_terms, options, out)
            })
        }
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

        // ------------------------------------------------------------------
        // THE THREE-STEP UI (the owner's doctrine: map -> partition ->
        // call). Thin wrappers over the existing machinery - every
        // step composes the committed commands, so every committed
        // gate applies verbatim underneath.
        // ------------------------------------------------------------------
        Command::Partition {
            name,
            chrom,
            start,
            end,
            sources,
            component,
            axis,
            partition_graphs,
            census,
            read_lengths,
            global_panel,
            global_routes,
            global_agc,
            record_origin,
            out_dir,
        } => {
            ensure(start < end, format!("the partition range is empty: {chrom}:{start}-{end}"))?;
            std::fs::create_dir_all(&out_dir)?;
            // THE ONE PARTITION ARTIFACT: the four-column BED
            // (chrom, start, end, PARTITION_NAME). Each row a named
            // reference range; a name is admitted once, and same-
            // chromosome rows never overlap (the territory grid).
            let bed = out_dir.join("partitions.bed");
            let mut bed_rows = read_partitions_bed(&bed)?;
            ensure(
                !bed_rows.iter().any(|row| row.3 == name),
                format!("the partition name {name:?} is already in {}", bed.display()),
            )?;
            for row in &bed_rows {
                ensure(
                    row.0 != chrom || row.2 <= start || end <= row.1,
                    format!(
                        "the partition {chrom}:{start}-{end} overlaps the existing \
                         partition {} ({}:{}-{})",
                        row.3,
                        row.0,
                        row.1,
                        row.2
                    ),
                )?;
            }
            let dir = out_dir.join(&name);
            ensure(
                !dir.exists(),
                format!(
                    "the partition directory {} already exists (the committed-receipt \
                     immutability discipline); remove it first if it is a discarded run",
                    dir.display()
                ),
            )?;
            std::fs::create_dir_all(&dir)?;
            // THE GLOBAL-PANEL CALL MODE (lever 1): under
            // --global-panel/--global-routes the per-territory local
            // syng and per-territory routes inventory are NOT BUILT —
            // the cut and the call bind to the ONE global panel and
            // the ONE whole-panel routes inventory (the registry
            // records them with the frame marker); the window beds
            // are regenerated verbatim either way (the locality maps
            // are the durable panel-side asset).
            let (syng, syng_seconds, routes, routes_seconds, frame) = match (
                &global_panel,
                &global_routes,
            ) {
                (Some(panel), Some(routes_dir)) => (
                    panel.clone(),
                    0.0,
                    routes_dir.display().to_string(),
                    0.0,
                    "global",
                ),
                _ => {
                    // The panel-side precompute 1: the local syng (the
                    // committed recipe verbatim: k=63, s=8, seed 7, pos64).
                    let syng = dir.join("locus-pos64.syng");
                    let started = std::time::Instant::now();
                    run_self(&[
                        "syng",
                        "-f",
                        &sources.display().to_string(),
                        "--syncmer-length",
                        "63",
                        "--smer-length",
                        "8",
                        "--syncmer-seed",
                        "7",
                        "--position-sample-rate",
                        "64",
                        "-o",
                        &syng.display().to_string(),
                    ])?;
                    let syng_seconds = started.elapsed().as_secs_f64();
                    // The panel-side precompute 3 (local frame only):
                    // the call-path routes inventory (the lookup lever;
                    // the census cross-check).
                    let routes = dir.join("routes");
                    let read_lengths = read_lengths
                        .iter()
                        .map(|l| l.to_string())
                        .collect::<Vec<_>>()
                        .join(",");
                    let started = std::time::Instant::now();
                    run_self(&[
                        "genome-infer",
                        "derive-routes",
                        "--panel",
                        &syng.display().to_string(),
                        "--sources",
                        &sources.display().to_string(),
                        "--read-lengths",
                        &read_lengths,
                        "--census",
                        &census.display().to_string(),
                        "--out-dir",
                        &routes.display().to_string(),
                    ])?;
                    let routes_seconds = started.elapsed().as_secs_f64();
                    (
                        syng.display().to_string(),
                        syng_seconds,
                        routes.display().to_string(),
                        routes_seconds,
                        "local",
                    )
                }
            };
            // The panel-side precompute 2: the window beds - the
            // region grid WITHIN the territory, regenerated verbatim
            // from the committed locality maps (the maps are the
            // durable panel-side asset; the member rows per window in
            // the map's own order).
            //
            // THE BED FILENAME IS THE AXIS INTERVAL'S OWN GROUP NAME
            // (the sparse-community fix, found by the whole-genome
            // run at T00017/T00032/T00000: the export-profile axis
            // names its windows by the ORIGINAL w10k/d1k community
            // ids - partition0, partition21, partition22, partition23
            // - while the locality maps renumber their windows 0..k;
            // naming the regenerated bed by the map index leaves the
            // cut's group resolution (axis group -> bed file) with a
            // NotFound at every sparse territory. The pairing is BY
            // ORDER with the strongest cross-check: the i-th map
            // (axis order) must carry the component's own row at the
            // i-th interval's extent, fail-closed loudly - verified
            // across BOTH substrate classes (L1's dense 25; T00017's
            // sparse 4; T00032/T00000 mixed).
            let beds = dir.join("beds");
            std::fs::create_dir_all(&beds)?;
            let mut windows = 0u64;
            let started = std::time::Instant::now();
            let mut maps: Vec<(u32, Vec<(String, u64, u64)>)> = Vec::new();
            for entry in std::fs::read_dir(&partition_graphs)? {
                let path = entry?.path();
                let file = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
                if !file.ends_with(".gfa.map.json") {
                    continue;
                }
                let map: serde_json::Value = genome::read_json(&path)?;
                let partition = map["partition"]
                    .as_u64()
                    .ok_or_else(|| invalid(format!("a locality map carries no partition id: {file}")))?
                    as u32;
                let mut rows = Vec::new();
                for member in map["members"].as_array().unwrap_or(&Vec::new()) {
                    let (path_name, mstart, mend) = (
                        member["path_name"].as_str(),
                        member["start"].as_u64(),
                        member["end"].as_u64(),
                    );
                    if let (Some(path_name), Some(mstart), Some(mend)) = (path_name, mstart, mend)
                    {
                        rows.push((path_name.to_string(), mstart, mend));
                    }
                }
                ensure(
                    !rows.is_empty(),
                    format!("a locality map carries no member rows: {file}"),
                )?;
                windows += 1;
                maps.push((partition, rows));
            }
            ensure(
                !maps.is_empty(),
                format!("no locality maps in {}", partition_graphs.display()),
            )?;
            maps.sort_by_key(|(partition, _)| *partition);
            // the axis intervals: the authoritative window group names,
            // in the call's own window order (extent-sorted)
            let axis_json: serde_json::Value = genome::read_json(&axis)?;
            let mut intervals: Vec<(u64, u64, String)> = axis_json["intervals"]
                .as_array()
                .ok_or_else(|| invalid("the axis carries no intervals array"))?
                .iter()
                .map(|interval| {
                    Ok((
                        interval["start"].as_u64().unwrap_or(0),
                        interval["end"].as_u64().unwrap_or(0),
                        interval["group"]
                            .as_str()
                            .ok_or_else(|| invalid("an axis interval carries no group"))?
                            .to_string(),
                    ))
                })
                .collect::<io::Result<Vec<_>>>()?;
            intervals.sort_by_key(|(start, end, _)| (*start, *end));
            ensure(
                intervals.len() == maps.len(),
                format!(
                    "the axis carries {} intervals but the locality maps carry {} windows \
                     (the sparse-community pairing is by order; a count mismatch is a \
                     substrate error, loudly named)",
                    intervals.len(),
                    maps.len()
                ),
            )?;
            for ((interval_start, interval_end, group), (_, rows)) in intervals.iter().zip(&maps) {
                ensure(
                    rows.iter().any(|(path_name, mstart, mend)| {
                        path_name == &component && *mstart == *interval_start && *mend == *interval_end
                    }),
                    format!(
                        "the locality map for axis group {group} carries no member row for \
                         the component {component} at its own extent \
                         {interval_start}-{interval_end} (the pairing is by order with the \
                         strongest cross-check; a mismatch is a substrate error, loudly named)"
                    ),
                )?;
                let mut text = String::new();
                for (path_name, mstart, mend) in rows {
                    text.push_str(&format!("{path_name}\t{mstart}\t{mend}\n"));
                }
                std::fs::write(beds.join(format!("{group}.bed")), text)?;
            }
            let beds_seconds = started.elapsed().as_secs_f64();
            // THE REGISTRY (keyed by the partition name): everything
            // the MAP and CALL steps consume. Under the global-panel
            // mode the syng/routes entries name the ONE global pair
            // and the frame marker states it; the committed local
            // registry is unchanged otherwise.
            let mut manifest = serde_json::json!({
                "version": 1,
                "name": name,
                "chrom": chrom,
                "start": start,
                "end": end,
                "sources": sources,
                "component": component,
                "axis": axis,
                "partition_graphs": partition_graphs,
                "beds": beds,
                "syng": syng,
                "routes": routes,
                "frame": frame,
                "record_origin": record_origin,
                "windows": windows,
            });
            // THE FLEET CALL PATH'S OWN REGISTRY ENTRY: the global
            // panel's AGC is recorded ONLY in the global frame (the
            // committed local registry's bytes are unchanged).
            if let Some(agc) = &global_agc {
                manifest["agc"] = serde_json::json!(agc);
            }
            genome::write_json(&dir.join("partition.json"), &manifest)?;
            let row = format!("{chrom}\t{start}\t{end}\t{name}\n");
            let mut text = std::fs::read_to_string(&bed).unwrap_or_default();
            text.push_str(&row);
            std::fs::write(&bed, text)?;
            bed_rows.push((chrom.clone(), start, end, name.clone()));
            let routes_manifest = std::path::Path::new(&routes).join("manifest.json");
            let receipt = serde_json::json!({
                "partition": name,
                "bed": bed,
                "bed4_rows": bed_rows.len(),
                "windows": windows,
                "frame": frame,
                "syng_seconds": syng_seconds,
                "beds_seconds": beds_seconds,
                "routes_seconds": routes_seconds,
                "routes_receipt": if routes_manifest.is_file() {
                    genome::reconstruction::fingerprint(&routes_manifest)?
                } else {
                    serde_json::Value::Null
                },
                "manifest": manifest,
                "scope": "panel-side-precomputes-keyed-by-partition-name",
            });
            println!("{}", serde_json::to_string_pretty(&receipt)?);
            Ok(())
        }
        Command::PanelPrep {
            name,
            chrom,
            start,
            end,
            sources,
            component,
            axis,
            partition_graphs,
            census,
            read_lengths,
            global_panel,
            global_agc,
            global_routes,
            record_origin,
            out_dir,
        } => run_panel_prep_row(
            // The single form loads the panel inside the precompute
            // piece's own subprocess, exactly as it always did.
            None,
            &global_panel,
            &global_agc,
            &global_routes,
            &out_dir,
            PanelPrepRow {
                name,
                chrom,
                start,
                end,
                sources,
                component,
                axis,
                partition_graphs,
                census,
                read_lengths,
                record_origin,
            },
        ),
        Command::PanelPrepBatch {
            manifest,
            global_panel,
            global_agc,
            global_routes,
            out_dir,
        } => {
            let mut rows: Vec<PanelPrepRow> = Vec::new();
            for line in std::io::BufReader::new(std::fs::File::open(&manifest)?).lines() {
                let line = line?;
                if line.trim().is_empty() {
                    continue;
                }
                rows.push(serde_json::from_str(&line).map_err(|e| {
                    invalid(format!(
                        "panel-prep-batch manifest row parse in {}: {e}",
                        manifest.display()
                    ))
                })?);
            }
            ensure(!rows.is_empty(), "empty panel-prep-batch manifest")?;
            let started = std::time::Instant::now();
            let panel = SyngIndex::load(&global_panel, SyncmerParams::default())?;
            eprintln!(
                "[panel-prep-batch] panel loaded once in {:.1}s; {} territories",
                started.elapsed().as_secs_f64(),
                rows.len()
            );
            for row in rows {
                let one = std::time::Instant::now();
                let name = row.name.clone();
                run_panel_prep_row(
                    Some(&panel),
                    &global_panel,
                    &global_agc,
                    &global_routes,
                    &out_dir,
                    row,
                )?;
                eprintln!(
                    "[panel-prep-batch] {} done in {:.1}s",
                    name,
                    one.elapsed().as_secs_f64()
                );
            }
            Ok(())
        }
        Command::Map {
            panel,
            agc,
            shards_dir,
            partitions,
            reads,
            region,
            local_derive,
            out_dir,
        } => {
            let bed_rows = read_partitions_bed(&partitions.join("partitions.bed"))?;
            ensure(
                !bed_rows.is_empty(),
                format!("no partitions in {}", partitions.display()),
            )?;
            let selected = select_partition_rows(&bed_rows, region.as_deref())?;
            ensure(
                !out_dir.exists(),
                format!(
                    "the map directory {} already exists (the committed-receipt \
                     immutability discipline)",
                    out_dir.display()
                ),
            )?;
            std::fs::create_dir_all(&out_dir)?;
            let mut cut_receipts = Vec::new();
            // The batch cut's manifest rows: the committed single cut,
            // run per row after ONE panel load instead of N loads per
            // territory (the ~190s global panel load is 51% of the map
            // leg's wall).
            let mut cut_specs: Vec<serde_json::Value> = Vec::new();
            for index in &selected {
                let (_, _, _, name) = &bed_rows[*index];
                let manifest: serde_json::Value =
                    genome::read_json(&partitions.join(name).join("partition.json"))?;
                let cut_dir = out_dir.join(name);
                ensure(
                    !cut_dir.exists(),
                    format!("the cut directory {} already exists", cut_dir.display()),
                )?;
                // THE GLOBAL-FRAME REGISTRY'S OWN MAP CUT (the
                // panel-prep step's form): the registry carries no
                // per-territory local syng (prep ZERO deleted it by
                // design), but the fleet call step's RESTRICTED
                // census profile demands the territory's kept-read
                // set — the LOCAL cut's own evidence mass (the
                // §10.6.17 basis) — so the map step builds the local
                // syng HERE, into the map's own partition directory
                // (the PARTITION step's local-form recipe verbatim,
                // ~3s), and runs the committed LOCAL-frame cut against
                // it. The committed local-frame registry path is
                // unchanged (the registry's own syng, byte-for-byte).
                let frame = manifest["frame"].as_str().unwrap_or("local");
                let cut_panel: std::borrow::Cow<str> = if frame == "global" {
                    std::fs::create_dir_all(&cut_dir)?;
                    let syng = cut_dir.join("locus-pos64.syng");
                    let sources = manifest["sources"].as_str()
                        .ok_or_else(|| invalid("the partition manifest carries no sources"))?;
                    run_self(&[
                        "syng",
                        "-f",
                        sources,
                        "--syncmer-length",
                        "63",
                        "--smer-length",
                        "8",
                        "--syncmer-seed",
                        "7",
                        "--position-sample-rate",
                        "64",
                        "-o",
                        &syng.display().to_string(),
                    ])?;
                    std::borrow::Cow::Owned(syng.display().to_string())
                } else {
                    std::borrow::Cow::Borrowed(manifest["syng"].as_str().ok_or_else(
                        || invalid("the partition manifest carries no syng"))?)
                };
                cut_specs.push(serde_json::json!({
                    "name": name,
                    "panel": cut_panel.as_ref(),
                    "axis": manifest["axis"].as_str().ok_or_else(
                        || invalid("the partition manifest carries no axis"))?,
                    "bed_directory": manifest["beds"].as_str().ok_or_else(
                        || invalid("the partition manifest carries no beds"))?,
                    "component": manifest["component"].as_str().ok_or_else(
                        || invalid("the partition manifest carries no component"))?,
                    "sources": manifest["sources"].as_str().ok_or_else(
                        || invalid("the partition manifest carries no sources"))?,
                    "out_dir": cut_dir,
                }));
            }
            // ONE panel load for the whole batch (the committed single
            // cut, run per manifest row): the frame and record-origin
            // stay the helper's defaults, exactly as the per-territory
            // exec above left them.
            let batch_manifest = out_dir.join("cut-batch-seq.manifest.jsonl");
            let mut batch_rows = String::new();
            for row in &cut_specs {
                batch_rows.push_str(&serde_json::to_string(row)?);
                batch_rows.push('\n');
            }
            std::fs::write(&batch_manifest, batch_rows)?;
            let mut command = std::process::Command::new(sample_projection_bin()?);
            command
                .arg("cut-batch-seq")
                .arg("--global-panel").arg(&panel)
                .arg("--agc").arg(&agc)
                .arg("--shards-dir").arg(&shards_dir)
                .arg("--manifest").arg(&batch_manifest);
            if local_derive {
                command.env("IMPG_CUT_LOCAL_DERIVE", "1");
            }
            let started = std::time::Instant::now();
            let status = command
                .status()
                .map_err(|e| invalid(format!("the sample_projection helper failed to launch: {e}")))?;
            ensure(
                status.success(),
                format!(
                    "the map's batch cut of {} partitions failed (status {status:?})",
                    cut_specs.len()
                ),
            )?;
            let batch_seconds = started.elapsed().as_secs_f64();
            println!(
                "=== map: batch cut of {} partitions ({batch_seconds:.1}s, one panel load)",
                cut_specs.len()
            );
            for index in &selected {
                let (chrom, start, end, name) = &bed_rows[*index];
                let cut_dir = out_dir.join(name);
                let record = serde_json::json!({
                    "partition": name,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "seconds": batch_seconds,
                    "wall_scope": "batch",
                    "census": genome::reconstruction::fingerprint(
                        &cut_dir.join("cut-census.jsonl"))?,
                    "cache": {
                        "path": cut_dir.join("cut-cache.bin"),
                        "bytes": std::fs::metadata(cut_dir.join("cut-cache.bin"))?.len(),
                    },
                });
                println!("=== map: partition {name} cut (batch wall {batch_seconds:.1}s)");
                cut_receipts.push(record);
            }
            let manifest = serde_json::json!({
                "version": 1,
                "panel": panel,
                "shards_dir": shards_dir,
                "reads": reads,
                "partitions": selected
                    .iter()
                    .map(|&i| bed_rows[i].3.clone())
                    .collect::<Vec<_>>(),
                "cuts": cut_receipts,
                "local_derive": local_derive,
            });
            genome::write_json(&out_dir.join("manifest.json"), &manifest)?;
            println!("{}", serde_json::to_string_pretty(&manifest)?);
            Ok(())
        }
        Command::CallRegion {
            map,
            partitions,
            region,
            far_bin_prune,
            exclude_paths,
            flank_overlap,
            truth_qv,
            qv_biwfa,
            selection,
            global_panel,
            global_routes,
            rss_budget_gib,
            out_dir,
        } => {
            let map_manifest: serde_json::Value = genome::read_json(&map.join("manifest.json"))?;
            let reads: PathBuf = map_manifest["reads"]
                .as_str()
                .ok_or_else(|| invalid("the map manifest carries no reads"))?
                .into();
            let bed_rows = read_partitions_bed(&partitions.join("partitions.bed"))?;
            let selected = select_partition_rows(&bed_rows, region.as_deref())?;
            ensure(
                !out_dir.exists(),
                format!(
                    "the calls directory {} already exists (the committed-receipt \
                     immutability discipline)",
                    out_dir.display()
                ),
            )?;
            std::fs::create_dir_all(&out_dir)?;
            let mut calls_receipts = Vec::new();
            // ------------------------------------------------------------------
            // THE FLEET FORM (STEP 3's product default, the owner's 2-3
            // step doctrine): every GLOBAL-frame partition (the
            // PANEL-PREP registry) is called by the fork-after-load
            // fleet over the ONE shared panel — the parent loads the
            // panel once and stays the load owner; width-32 children
            // (the batched cut amortization) each build their own AGC
            // index, batch-cut their chunk's territories (the
            // RESTRICTED census profile, the basis AUTO-DERIVED from
            // the MAP step's own cut-reads set) and run their chunk's
            // calls concurrently at the T16 rung with the FAR-BIN
            // PRUNE default and the record origin AUTO-DERIVED from
            // the partition registry — no manual flags, no example
            // binary, no magic constants. The committed serial
            // per-partition path stays byte-for-byte for the LOCAL
            // frame. (The fleet forks BEFORE any parent-side rayon
            // touch — the prefork hygiene the fork driver's own audit
            // receipts — so the serial local-frame calls run AFTER
            // the children are reaped.)
            // ------------------------------------------------------------------
            let shards_dir: PathBuf = map_manifest["shards_dir"]
                .as_str()
                .ok_or_else(|| invalid("the map manifest carries no shards_dir"))?
                .into();
            let mut fleet_indexes: Vec<usize> = Vec::new();
            let mut local_indexes: Vec<usize> = Vec::new();
            for &index in &selected {
                let (_, _, _, name) = &bed_rows[index];
                let manifest: serde_json::Value =
                    genome::read_json(&partitions.join(name).join("partition.json"))?;
                if manifest["frame"].as_str() == Some("global") {
                    fleet_indexes.push(index);
                } else {
                    local_indexes.push(index);
                }
            }
            if !fleet_indexes.is_empty() {
                ensure(
                    truth_qv.is_none()
                        && exclude_paths.is_empty()
                        && flank_overlap.is_none()
                        && matches!(selection, realign_score::Selection::Locality),
                    "the fleet call path runs the product defaults (the test-mode and \
                     masking flags are the committed serial per-partition path's own; \
                     rerun without them, or call the partitions by the serial form)",
                )?;
                let mut specs: Vec<genome::sample_projection::CutBatchGlobalSpec> = Vec::new();
                let mut fleet_names: Vec<String> = Vec::new();
                let mut fleet_panel: Option<String> = None;
                let mut fleet_routes: Option<PathBuf> = None;
                let mut fleet_agc: Option<PathBuf> = None;
                for &index in &fleet_indexes {
                    let (_, _, _, name) = &bed_rows[index];
                    let manifest: serde_json::Value =
                        genome::read_json(&partitions.join(name).join("partition.json"))?;
                    let record_origin: u64 =
                        manifest["record_origin"].as_u64().unwrap_or(0);
                    let component = manifest["component"].as_str()
                        .ok_or_else(|| invalid("the partition manifest carries no component"))?
                        .to_string();
                    let graphs: PathBuf = manifest["partition_graphs"].as_str()
                        .ok_or_else(|| invalid("the partition manifest carries no graphs"))?
                        .into();
                    let axis: PathBuf = manifest["axis"].as_str()
                        .ok_or_else(|| invalid("the partition manifest carries no axis"))?
                        .into();
                    let beds: PathBuf = manifest["beds"].as_str()
                        .ok_or_else(|| invalid("the partition manifest carries no beds"))?
                        .into();
                    let panel = manifest["syng"].as_str()
                        .ok_or_else(|| invalid("the partition manifest carries no syng"))?
                        .to_string();
                    let routes: PathBuf = manifest["routes"].as_str()
                        .ok_or_else(|| invalid("the partition manifest carries no routes"))?
                        .into();
                    let agc: PathBuf = manifest["agc"].as_str()
                        .ok_or_else(|| invalid(
                            "the partition manifest carries no agc (the PANEL-PREP step \
                             records the panel AGC for the fleet call path; re-run \
                             panel-prep for this partition)",
                        ))?
                        .into();
                    match (&fleet_panel, &fleet_routes, &fleet_agc) {
                        (None, None, None) => {
                            fleet_panel = Some(panel.clone());
                            fleet_routes = Some(routes.clone());
                            fleet_agc = Some(agc.clone());
                        }
                        (Some(p), Some(r), Some(a)) => {
                            ensure(
                                p == &panel && r == &routes && a == &agc,
                                format!(
                                    "the fleet form loads ONE panel (partition {name} \
                                     names a different global pair than the earlier \
                                     partitions)",
                                ),
                            )?;
                        }
                        _ => unreachable!("the global pair is set atomically"),
                    }
                    // THE RESTRICTED CENSUS PROFILE'S BASIS, AUTO-DERIVED:
                    // the MAP step's own kept-read set for the partition
                    // (the §10.6.17 default demands the evidence-set
                    // statement; the map's cut IS the product's own).
                    let readset = map.join(name).join("cut-reads.fasta");
                    ensure(
                        readset.is_file(),
                        format!(
                            "the map carries no cut-reads.fasta for partition {name} \
                             (the restricted census profile demands the evidence-set \
                             statement; re-run the MAP step for this partition)",
                        ),
                    )?;
                    // THE SCORED LOCI: every window, or the --region
                    // range resolved against the partition's frame.
                    let axis_rows = realign_score::component_axis_rows(&graphs, &component)?;
                    let loci: Vec<u32> = match region_filter_range(region.as_deref())? {
                        Some(range) => {
                            let resolved = realign_score::resolve_region_loci(
                                &component,
                                &axis_rows,
                                range,
                            )?;
                            eprintln!(
                                "[call] partition {name}: region {range:?} resolved to \
                                 loci {resolved:?} (frame {component})"
                            );
                            resolved
                        }
                        None => (0..axis_rows.len() as u32).collect(),
                    };
                    if loci.is_empty() {
                        println!(
                            "=== call: partition {name} skipped (the region covers no \
                             window)"
                        );
                        continue;
                    }
                    let loci = loci
                        .iter()
                        .map(|l| l.to_string())
                        .collect::<Vec<_>>()
                        .join(",");
                    println!(
                        "=== call: partition {name} ({loci} windows, the FLEET form: \
                         fork-after-load, width-32 cut, T16 calls)"
                    );
                    let partition_out = out_dir.join(name);
                    // The fleet row's cut out_dir is a SIBLING of the
                    // call out_dir (never its parent): the call's own
                    // with_output_model CREATE_DIR-refuses an existing
                    // directory, and the cut phase materializes its
                    // out_dir (and the parent chain) first — <name>-cut
                    // beside <name>/ keeps both disciplines intact.
                    let cut_out = out_dir.join(format!("{name}-cut"));
                    specs.push(genome::sample_projection::CutBatchGlobalSpec {
                        name: name.clone(),
                        axis,
                        bed_directory: beds,
                        component,
                        out_dir: cut_out,
                        record_origin,
                        readset: Some(readset),
                        call: Some(genome::sample_projection::CutBatchGlobalCall {
                            reads: reads.clone(),
                            partition_graphs: graphs,
                            loci,
                            out_dir: partition_out.clone(),
                        }),
                    });
                    fleet_names.push(name.clone());
                }
                if !specs.is_empty() {
                    let fleet_panel = fleet_panel
                        .expect("a non-empty fleet row set carries the global pair");
                    let fleet_routes =
                        fleet_routes.expect("a non-empty fleet row set carries the routes");
                    let fleet_agc =
                        fleet_agc.expect("a non-empty fleet row set carries the agc");
                    // The explicit --global-panel/--global-routes flags,
                    // when passed, must name the registry's own pair
                    // (loud on mismatch; the registry is the record).
                    if let Some(flag_panel) = &global_panel {
                        ensure(
                            flag_panel == &fleet_panel,
                            format!(
                                "--global-panel {flag_panel} does not name the registry's \
                                 own panel {fleet_panel}"
                            ),
                        )?;
                    }
                    if let Some(flag_routes) = &global_routes {
                        ensure(
                            flag_routes == &fleet_routes,
                            format!(
                                "--global-routes {} does not name the registry's own \
                                 routes {}",
                                flag_routes.display(),
                                fleet_routes.display()
                            ),
                        )?;
                    }
                    genome::sample_projection::fork_cut_call_batch_global_specs(
                        &fleet_panel,
                        &fleet_agc,
                        &shards_dir,
                        specs,
                        &fleet_routes,
                        0.0,
                        32,
                        Some(16),
                        far_bin_prune,
                    )?;
                    for name in &fleet_names {
                        let partition_out = out_dir.join(name);
                        calls_receipts.push(serde_json::json!({
                            "partition": name,
                            "form": "fleet",
                            "calls": genome::reconstruction::fingerprint(
                                &partition_out.join("calls.jsonl"))?,
                            "sequences": genome::reconstruction::fingerprint(
                                &partition_out.join("calls.jsonl.sequences.jsonl"))?,
                            "cut_census": genome::reconstruction::fingerprint(
                                &out_dir
                                    .join(format!("{name}-cut"))
                                    .join("cut-census.jsonl"))?,
                            "truth_qv_test_mode": serde_json::Value::Null,
                        }));
                    }
                }
            }
            for index in &local_indexes {
                let (_, _, _, name) = &bed_rows[*index];
                let manifest: serde_json::Value =
                    genome::read_json(&partitions.join(name).join("partition.json"))?;
                // THE RECORD-ORIGIN CONVENTION (the registry's own
                // statement): the partition step recorded the
                // substrate's content-origin convention; the call
                // states the same one (the committed local frame
                // recorded 0 and never consumes it — the global
                // frame's maps translate at start+d).
                let record_origin: u64 = manifest["record_origin"].as_u64().unwrap_or(0);
                let component = manifest["component"].as_str()
                    .ok_or_else(|| invalid("the partition manifest carries no component"))?
                    .to_string();
                let graphs: PathBuf = manifest["partition_graphs"].as_str()
                    .ok_or_else(|| invalid("the partition manifest carries no graphs"))?
                    .into();
                // THE SCORED LOCI: every window of the partition, or
                // the --region range resolved against the partition's
                // own reference frame ("any locus, any position").
                let axis_rows = realign_score::component_axis_rows(&graphs, &component)?;
                let loci: Vec<u32> = match region_filter_range(region.as_deref())? {
                    Some(range) => {
                        let resolved =
                            realign_score::resolve_region_loci(&component, &axis_rows, range)?;
                        eprintln!(
                            "[call] partition {name}: region {range:?} resolved to loci \
                             {resolved:?} (frame {component})"
                        );
                        resolved
                    }
                    None => (0..axis_rows.len() as u32).collect(),
                };
                if loci.is_empty() {
                    println!(
                        "=== call: partition {name} skipped (the region covers no window)"
                    );
                    continue;
                }
                let loci = loci
                    .iter()
                    .map(|l| l.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                println!(
                    "=== call: partition {name} ({loci} windows, frame {component})"
                );
                let partition_out = out_dir.join(name);
                // THE GLOBAL-PANEL CALL MODE (lever 1): the ONE global
                // panel + the ONE whole-panel routes inventory replace
                // the registry's per-territory syng/routes; the
                // engine's --global-frame switch carries the frame
                // translation (maps on load, product rows on emission).
                let (panel, routes, global_frame) = match (&global_panel, &global_routes) {
                    (Some(panel), Some(routes)) => {
                        println!(
                            "=== call: partition {name} against the GLOBAL panel {}",
                            panel
                        );
                        (panel.clone(), routes.clone(), true)
                    }
                    _ => (
                        manifest["syng"].as_str()
                            .ok_or_else(|| invalid("the partition manifest carries no syng"))?
                            .into(),
                        manifest["routes"].as_str()
                            .ok_or_else(|| invalid("the partition manifest carries no routes"))?
                            .into(),
                        false,
                    ),
                };
                run(Command::CallDiplotypes {
                    panel,
                    routes,
                    partition_graphs: graphs,
                    component,
                    census: map.join(name).join("cut-census.jsonl"),
                    reads: reads.clone(),
                    derive_cache: map.join(name).join("cut-cache.bin"),
                    loci: Some(loci),
                    region: None,
                    selection: selection.clone(),
                    rss_budget_gib,
                    truth_qv: truth_qv.clone(),
                    qv_biwfa: qv_biwfa.clone(),
                    // THE CALL STEP'S REFERENCE-ANCHORED PRODUCT: the
                    // per-locus hap1/hap2 diplotype sequences are
                    // part of the call (the owner's doctrine), so
                    // the emission is always on here.
                    emit_sequences: true,
                    exclude_paths: exclude_paths.clone(),
                    flank_overlap,
                    global_frame,
                    record_origin,
                    // THE PRODUCT DEFAULT (the §10.6.18 ruling
                    // executed): the global-frame call path prunes
                    // the far-bin class by default (the opt-out
                    // --far-bin-prune=false stays); the committed
                    // LOCAL-frame path is byte-frozen and never
                    // prunes (the frame gate in run_injected holds
                    // both forms).
                    far_bin_prune: global_frame && far_bin_prune,
                    out_dir: partition_out.clone(),
                })?;
                calls_receipts.push(serde_json::json!({
                    "partition": name,
                    "calls": genome::reconstruction::fingerprint(
                        &partition_out.join("calls.jsonl"))?,
                    "sequences": genome::reconstruction::fingerprint(
                        &partition_out.join("calls.jsonl.sequences.jsonl"))?,
                    "truth_qv_test_mode": if truth_qv.is_some() {
                        genome::reconstruction::fingerprint(
                            &partition_out.join("calls.jsonl.truth-qv.jsonl"))?
                    } else {
                        serde_json::Value::Null
                    },
                }));
            }
            let receipt = serde_json::json!({
                "version": 1,
                "map": map,
                "partitions": partitions,
                "region": region,
                "exclude_paths": exclude_paths,
                "flank_overlap": flank_overlap.map(|f| f.setting()),
                "far_bin_prune": far_bin_prune,
                "global_panel": global_panel,
                "global_routes": global_routes,
                "partitions_called": calls_receipts.len(),
                "calls": calls_receipts,
                "scope": "per-partition-diplotype-calls-reference-anchored",
            });
            genome::write_json(&out_dir.join("manifest.json"), &receipt)?;
            println!("{}", serde_json::to_string_pretty(&receipt)?);
            Ok(())
        }
    }
}

fn ensure(ok: bool, message: impl Into<String>) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(message))
    }
}

/// THE BED4 READER (the one partition artifact: chrom, start, end,
/// PARTITION_NAME - each row a named reference range).
fn read_partitions_bed(path: &Path) -> io::Result<Vec<(String, u64, u64, String)>> {
    let text = match std::fs::read_to_string(path) {
        Ok(text) => text,
        Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(Vec::new()),
        Err(e) => return Err(e),
    };
    let mut rows = Vec::new();
    for (line_number, line) in text.lines().enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        ensure(
            fields.len() == 4 && !fields.iter().any(|f| f.is_empty()),
            format!(
                "invalid BED4 row {} in {}: expected chrom\\tstart\\tend\\tPARTITION_NAME",
                line_number + 1,
                path.display()
            ),
        )?;
        let (start, end) = (fields[1].parse::<u64>(), fields[2].parse::<u64>());
        let (start, end) = match (start, end) {
            (Ok(start), Ok(end)) => (start, end),
            _ => {
                return Err(invalid(format!(
                    "invalid BED4 coordinates in row {}: {}",
                    line_number + 1,
                    line
                )))
            }
        };
        rows.push((fields[0].to_string(), start, end, fields[3].to_string()));
    }
    Ok(rows)
}

/// THE REGION FILTER (the owner's ruling: the interface is
/// coordinates): a comma-list of partition names selects those BED
/// rows; a chr:start-end range selects the BED rows overlapping it;
/// absent selects every row.
fn select_partition_rows(
    bed_rows: &[(String, u64, u64, String)],
    region: Option<&str>,
) -> io::Result<Vec<usize>> {
    let Some(region) = region else {
        return Ok((0..bed_rows.len()).collect());
    };
    if let Some(range) = region_filter_range(Some(region))? {
        let bad_region = || invalid(&format!("the region must be chr:start-end: {region:?}"));
        let (chrom, bounds) = range.split_once(':').ok_or_else(bad_region)?;
        let (lo, hi) = bounds
            .split_once('-')
            .and_then(|(s, e)| Some((s.parse::<u64>().ok()?, e.parse::<u64>().ok()?)))
            .ok_or_else(bad_region)?;
        let selected: Vec<usize> = bed_rows
            .iter()
            .enumerate()
            .filter(|(_, row)| row.0 == chrom && row.1 < hi && row.2 > lo)
            .map(|(index, _)| index)
            .collect();
        ensure(
            !selected.is_empty(),
            format!("the region {region:?} covers no partition of the BED"),
        )?;
        return Ok(selected);
    }
    let names: Vec<&str> = region.split(',').map(|s| s.trim()).collect();
    // NAME the absent partitions. A literal "{region:?}" here prints the
    // braces (ensure takes impl Into<String>, so it compiles), and the
    // failure it hides is usually a RACE: the PANEL-PREP step appends each
    // territory's row to partitions.bed as it completes, so a map leg
    // started while a panel-prep backlog is running can see a partition
    // whose directory exists but whose BED row is not yet written.
    let missing: Vec<&&str> = names
        .iter()
        .filter(|name| !bed_rows.iter().any(|row| row.3 == **name))
        .collect();
    ensure(
        missing.is_empty(),
        format!(
            "the region names {missing:?} name partitions absent from the BED \
             ({} rows present; a concurrent PANEL-PREP writes these rows)",
            bed_rows.len()
        ),
    )?;
    let selected: Vec<usize> = bed_rows
        .iter()
        .enumerate()
        .filter(|(_, row)| names.contains(&row.3.as_str()))
        .map(|(index, _)| index)
        .collect();
    ensure(
            !selected.is_empty(),
            format!("the region {region:?} selects no partition"),
        )?;
    Ok(selected)
}

/// The region's coordinate-range form, when it is one: a spec
/// carrying `chrom:` is a chr:start-end range, not a name list.
fn region_filter_range(region: Option<&str>) -> io::Result<Option<&str>> {
    match region {
        None => Ok(None),
        Some(region) => Ok(region.contains(':').then_some(region)),
    }
}

/// Run this binary with the given subcommand args (the thin-wrapper
/// composition: the committed commands, as they are, in-process
/// semantics identical to typing them).
/// ONE PANEL-PREP TERRITORY'S ROW: the per-territory half of the
/// panel-prep step (name, range, substrate, component, axis, locality
/// maps, census). The other half — the global panel, its AGC and the
/// whole-panel routes inventory — is shared by every row of a batch and
/// is the caller's. The single form builds one row from its flags;
/// `panel-prep-batch` reads one row per manifest line, and BOTH call
/// `run_panel_prep_row`, so the two forms cannot drift apart.
#[derive(serde::Deserialize)]
struct PanelPrepRow {
    name: String,
    chrom: String,
    start: u64,
    end: u64,
    sources: PathBuf,
    component: String,
    axis: PathBuf,
    partition_graphs: PathBuf,
    census: PathBuf,
    /// The routes inventory's read lengths. Default 148, the committed
    /// default of the single form's flag.
    #[serde(default = "panel_prep_default_read_lengths")]
    read_lengths: Vec<u64>,
    /// The content-origin convention the registry records. Default 0,
    /// the single form's own default.
    #[serde(default)]
    record_origin: u64,
}

fn panel_prep_default_read_lengths() -> Vec<u64> {
    vec![148]
}

/// ONE TERRITORY'S PANEL-PREP, the committed three pieces in sequence.
/// THE COST CENTRE (measured: 271.6s for T00000 at RAYON_NUM_THREADS=8):
/// piece (3) shells out to `genome-infer precompute-panel-loci`, whose
/// `SyngIndex::load` of the ONE global panel costs 195.9s while the
/// derivation behind it is 21.4s. That load is a pure function of the
/// panel, so `panel` lets the batched form hand in the ONE load it
/// already paid; `None` is the single form, whose precompute piece keeps
/// its own subprocess and its own load, unchanged.
/// THE REUSES this must not defeat (both stated, both loud): the
/// whole-panel routes inventory is skipped when `global_routes`
/// already carries a manifest.json, and the per-locus precompute asset
/// is skipped when its manifest names the current format. A batch whose
/// rows are already prepped therefore does no panel work at all beyond
/// the one load.
fn run_panel_prep_row(
    panel: Option<&SyngIndex>,
    global_panel: &str,
    global_agc: &Path,
    global_routes: &Path,
    out_dir: &Path,
    row: PanelPrepRow,
) -> io::Result<()> {
    let PanelPrepRow {
        name,
        chrom,
        start,
        end,
        sources,
        component,
        axis,
        partition_graphs,
        census,
        read_lengths,
        record_origin,
    } = row;
    // STEP 1's OWN TWO PIECES, IN SEQUENCE (the owner's 2-3
    // step doctrine: ONE panel-prep step, once per panel,
    // reusable forever):
    // (1) THE WHOLE-PANEL ROUTES INVENTORY (derive-routes over
    // the ONE global panel): built at the FIRST invocation,
    // REUSED by every later one — an existing inventory is
    // skipped, loudly (the routes inventory is a panel-side
    // asset, never rebuilt per partition row).
    if global_routes.join("manifest.json").is_file() {
        println!(
            "=== panel-prep: the whole-panel routes inventory {} exists — REUSED \
             (built once per panel)",
            global_routes.display()
        );
    } else {
        println!(
            "=== panel-prep: building the whole-panel routes inventory into {} \
             (derive-routes over the ONE global panel; the once-per-panel step)",
            global_routes.display()
        );
        let read_lengths = read_lengths
            .iter()
            .map(|l| l.to_string())
            .collect::<Vec<_>>()
            .join(",");
        run_self(&[
            "genome-infer",
            "derive-routes",
            "--panel",
            global_panel,
            "--census",
            &census.display().to_string(),
            "--sources",
            &global_agc.display().to_string(),
            "--read-lengths",
            &read_lengths,
            "--out-dir",
            &global_routes.display().to_string(),
        ])?;
    }
    // (2) THE PARTITION GLOBAL FORM (the BED4 row + the window
    // beds + the registry): the committed PARTITION step's own
    // machinery, verbatim, with the global pair + the panel
    // AGC recorded so the CALL step needs no manual flags.
    run(Command::Partition {
        name,
        chrom,
        start,
        end,
        sources,
        component: component.clone(),
        axis,
        partition_graphs: partition_graphs.clone(),
        census,
        read_lengths,
        global_panel: Some(global_panel.to_string()),
        global_routes: Some(global_routes.to_path_buf()),
        global_agc: Some(global_agc.to_path_buf()),
        record_origin,
        out_dir: out_dir.to_path_buf(),
    })?;
    // (3) THE PANEL-SIDE PER-LOCUS PRECOMPUTE (SS10.6.23 slice 1
    // + SS10.6.25 slice 2 — the Frontier-2 increment's own
    // panel-prep pieces): the row context + the per-window fold
    // groupings + the bind scaffolding, keyed on content into
    // <partition-graphs>/panel-precompute-loci — the
    // once-per-panel price the CALL step's serve amortizes to
    // zero per sample. An asset of the CURRENT format is REUSED
    // loudly (the routes-inventory precedent); a slice-1 /1
    // asset is UPGRADED loudly (rebuilt with the
    // bind-scaffolding section); the call falls back to the
    // live path loudly on any content-key mismatch.
    let precompute_dir = partition_graphs.join("panel-precompute-loci");
    let precompute_format = std::fs::read_to_string(precompute_dir.join("manifest.json"))
        .ok()
        .and_then(|text| {
            serde_json::from_str::<genome::realign_score::PanelPrecomputeManifest>(&text).ok()
        })
        .map(|manifest| manifest.format_text());
    if precompute_format.as_deref() == Some("impg-panel-precompute-loci/2") {
        println!(
            "=== panel-prep: the per-locus precompute asset {} exists — REUSED \
             (built once per panel)",
            precompute_dir.display()
        );
    } else {
        if let Some(format) = precompute_format {
            println!(
                "=== panel-prep: the existing asset {} is format {:?} (slice 1's) — \
                 REBUILDING with the bind-scaffolding section",
                precompute_dir.display(),
                format
            );
        } else {
            println!(
                "=== panel-prep: building the per-locus precompute asset into {} \
                 (the row context + the fold layouts + the bind scaffolding; the \
                 once-per-panel price)",
                precompute_dir.display()
            );
        }
        match panel {
            // The batched form's ONE load, handed to the SAME function
            // the subprocess would have run — the same panel, the same
            // arguments, so the same asset.
            Some(panel) => genome::realign_score::precompute_panel_loci(
                global_panel,
                global_routes,
                &partition_graphs,
                &component,
                true,
                record_origin,
                &precompute_dir,
                Some(panel),
            )?,
            None => run_self(&[
                "genome-infer",
                "precompute-panel-loci",
                "--panel",
                global_panel,
                "--routes",
                &global_routes.display().to_string(),
                "--partition-graphs",
                &partition_graphs.display().to_string(),
                "--component",
                &component,
                "--global-frame",
                "--record-origin",
                &record_origin.to_string(),
            ])?,
        }
    }
    Ok(())
}

fn run_self(args: &[&str]) -> io::Result<()> {
    let exe = std::env::current_exe()?;
    let status = std::process::Command::new(&exe)
        .args(args)
        .status()
        .map_err(|e| invalid(format!("failed to launch {}: {e}", exe.display())))?;
    ensure(
        status.success(),
        format!("the command {:?} failed (status {status:?})", args),
    )
}

/// The sample_projection helper binary (the projection machinery:
/// project + dedup + cut), shipped beside the impg binary.
fn sample_projection_bin() -> io::Result<PathBuf> {
    let exe = std::env::current_exe()?;
    let dir = exe
        .parent()
        .ok_or_else(|| invalid("no parent of the running binary"))?;
    let candidate = dir.join("examples").join("sample_projection");
    ensure(
        candidate.is_file(),
        format!(
            "the sample_projection helper binary is absent at {} \
             (build it: cargo build --release --bin impg --example sample_projection)",
            candidate.display()
        ),
    )?;
    Ok(candidate)
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Subcommand as _;

    /// Parse a `genome-infer` subcommand invocation (the CLI surface
    /// itself, the interface-change tests' own form: the flag defaults
    /// are product behavior and get their own proof).
    fn parse_subcommand(args: &[&str]) -> Command {
        let cmd =
            Command::augment_subcommands(clap::Command::new("genome-infer"));
        let matches = cmd
            .try_get_matches_from(
                std::iter::once("genome-infer").chain(args.iter().copied()),
            )
            .expect("the invocation parses");
        <Command as clap::FromArgMatches>::from_arg_matches(&matches)
            .expect("the subcommand is known")
    }

    /// The call-diplotypes base form (every required input, the
    /// committed local-frame shape).
    fn call_diplotypes_args<'a>(extra: &'a [&'a str]) -> Vec<&'a str> {
        let mut args: Vec<&str> = [
            "call-diplotypes",
            "--panel", "p",
            "--routes", "r",
            "--partition-graphs", "g",
            "--component", "c",
            "--census", "n",
            "--reads", "q",
            "--derive-cache", "d",
            "--loci", "1",
            "--out-dir", "o",
        ]
        .to_vec();
        args.extend_from_slice(extra);
        args
    }

    #[test]
    fn far_bin_prune_defaults_on_and_opts_out_by_explicit_false() {
        // THE PRODUCT DEFAULT FLIP (the owner's ruling on the §10.6.18
        // gate receipts): the working form IS the default — the bare
        // invocation prunes; the opt-out form `--far-bin-prune=false`
        // stays for the record (the committed unpruned receipts);
        // the legacy bare `--far-bin-prune` keeps meaning ON.
        let Command::CallDiplotypes { far_bin_prune, .. } =
            parse_subcommand(&call_diplotypes_args(&[]))
        else {
            panic!("expected call-diplotypes");
        };
        assert!(far_bin_prune, "the default is ON");
        let Command::CallDiplotypes { far_bin_prune, .. } = parse_subcommand(
            &call_diplotypes_args(&["--far-bin-prune=false"]),
        )
        else {
            panic!("expected call-diplotypes");
        };
        assert!(!far_bin_prune, "the explicit opt-out is OFF");
        let Command::CallDiplotypes { far_bin_prune, .. } = parse_subcommand(
            &call_diplotypes_args(&["--far-bin-prune"]),
        )
        else {
            panic!("expected call-diplotypes");
        };
        assert!(far_bin_prune, "the bare flag stays ON");
    }

    #[test]
    fn call_region_far_bin_prune_defaults_on() {
        // The three-step CALL wrapper carries the same product
        // default (the frame gate in the engine keeps the committed
        // local-frame path byte-frozen; the effective value at the
        // local frame is forced OFF by the wrapper).
        let args = [
            "call-region",
            "--map", "m",
            "--partitions", "p",
            "--out-dir", "o",
        ];
        let Command::CallRegion { far_bin_prune, .. } = parse_subcommand(&args)
        else {
            panic!("expected call-region");
        };
        assert!(far_bin_prune, "the default is ON");
    }

    #[test]
    fn panel_prep_batch_row_defaults_are_the_single_forms_flag_defaults() {
        // A manifest row that omits them must state the SAME panel-prep
        // parameters the single form's flags default to -- read_lengths
        // [148] (--read-lengths) and record_origin 0 (--record-origin) --
        // or the batched form would silently precompute against a
        // different parameter pair than the territory's single-form run.
        let row: PanelPrepRow = serde_json::from_str(
            r#"{"name":"T1","chrom":"chr1","start":0,"end":25000,
                "sources":"locus.fa","component":"HG002#1#chr1:1-2",
                "axis":"axis.json","partition_graphs":"maps",
                "census":"census.jsonl"}"#,
        )
        .expect("a minimal row parses");
        assert_eq!(row.name, "T1");
        assert_eq!(row.read_lengths, vec![148]);
        assert_eq!(row.record_origin, 0);

        let explicit: PanelPrepRow = serde_json::from_str(
            r#"{"name":"T1","chrom":"chr1","start":0,"end":25000,
                "sources":"locus.fa","component":"HG002#1#chr1:1-2",
                "axis":"axis.json","partition_graphs":"maps",
                "census":"census.jsonl",
                "read_lengths":[148,250],"record_origin":1}"#,
        )
        .expect("an explicit row parses");
        assert_eq!(explicit.read_lengths, vec![148, 250]);
        assert_eq!(explicit.record_origin, 1);
    }
}
