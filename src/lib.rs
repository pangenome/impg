// lib.rs
pub mod agc_index;
pub mod alignment_record;
pub mod commands;
pub mod faidx;
pub mod forest_map;
pub mod graph;
pub mod impg;
pub mod impg_index;
pub mod multi_impg;
pub mod onealn;
pub mod paf;
pub mod realize;
pub mod tpa_parser;
pub mod seqidx;
pub mod sequence_index;
pub mod subset_filter;

/// GFA engine selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum GfaEngine {
    /// Recursive partitioning + POA + lacing
    Recursive,
    /// Seqwish graph induction via transitive closure
    Seqwish,
    /// Flat single-pass partial order alignment
    Poa,
}

/// Resolved engine configuration passed to subcommand functions.
pub struct EngineOpts {
    pub engine: GfaEngine,
    pub recursive_config: Option<realize::RealizeConfig>,
    pub num_threads: usize,
}
