//! Automatic all-source routing; fixed implicit-universe background-relative counts.
mod build;
mod evaluate;
mod graph;
mod import;
mod profiles;
mod search;
mod storage;
use super::{checksum, read_json, write_json, PanelIdentity, COUNT_POLICY};
use crate::{
    sample_mem_bwt::{canonical, encode_walk, invalid},
    syng::SyngIndex,
};
pub use build::build;
pub use evaluate::{Assignment, Evaluation, Evaluator, Route, Segment};
pub use graph::{sources_fetch_probe_stats, Family, Graph, Lane, Ports, Sources};
pub use search::{search, SearchBudget, SearchResult};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};
pub use storage::Port;
pub const MODEL: &str = "panel-route-relative-poisson-v1";
pub const VERSION: u32 = 1;
pub const RULE: &str = "all-native-identity-end-pairings;all-donors;raw-view-coordinate-union;equal-oriented-full-DNA;cut=floor(k/2);positive-traversals;canonical-source-span-capacity-v1";
pub const UNIVERSE: &str = "all-two-node-features-realizable-under-frozen-route-and-read-length-rules UNION full-registry UNION complete-sample-positive;uninstantiated-zero-signal-zero-count=symbolic-zero";
pub(crate) fn require(ok: bool, message: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(message))
    }
}
fn add(n: &mut u64, x: u64) -> io::Result<()> {
    *n = n
        .checked_add(x)
        .ok_or_else(|| invalid("route integer overflow"))?;
    Ok(())
}
fn digest(value: &impl Serialize) -> io::Result<String> {
    Ok(format!(
        "{:016x}",
        checksum(&serde_json::to_vec(value).map_err(io::Error::other)?)
    ))
}
pub fn compiler_identity() -> String {
    let mut hash = checksum(super::joint::compiler_identity().as_bytes());
    for s in [
        include_str!("mod.rs"),
        include_str!("build.rs"),
        include_str!("graph.rs"),
        include_str!("import.rs"),
        include_str!("profiles.rs"),
        include_str!("evaluate.rs"),
        include_str!("search.rs"),
        include_str!("storage.rs"),
    ] {
        super::hash_update(&mut hash, s.as_bytes());
    }
    format!("panel-route-native-start-v1-{hash:016x}")
}
#[cfg(test)]
mod tests;
