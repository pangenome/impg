//! Compatible task payload shapes only; no production scheduler or score import.
use super::*;
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Context {
    pub(super) family: usize,
    pub(super) slot: usize,
    pub(super) completed: Vec<Route>,
    pub(super) segments: Vec<Segment>,
    pub(super) source: usize,
    pub(super) cut: u64,
    pub(super) reverse: bool,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) enum AfterCheck {
    Close,
    Hub(Port),
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) enum Op {
    Start,
    SourceBound {
        lo: u64,
        hi: u64,
    },
    SourceScan {
        base: u64,
        permutation: SourcePermutation,
    },
    Check {
        piece: Segment,
        index: usize,
        next: AfterCheck,
    },
    Close {
        piece: Segment,
    },
    HubBound {
        port: Port,
        lower: Option<u64>,
        lo: u64,
        hi: u64,
    },
    HubScan {
        port: Port,
        base: u64,
        permutation: Permutation,
    },
    Child {
        piece: Segment,
        donor: Port,
    },
    Probe {
        routes: Vec<Route>,
        next_slot: usize,
    },
    ProbeCheck {
        assignment: Assignment,
        i: usize,
        j: usize,
    },
    Evaluate {
        assignment: Assignment,
    },
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Task {
    pub(super) id: u64,
    pub(super) ready: u64,
    pub(super) depth: usize,
    pub(super) prev: Option<u64>,
    pub(super) next: Option<u64>,
    pub(super) context: Context,
    pub(super) op: Op,
}
