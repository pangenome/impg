//! THE ONE WHOLE-SAMPLE PROJECTION — the thin example shell.
//!
//! The implementation lives in the library
//! (`impg::genome_inference::sample_projection`) since the product
//! CLI's three-step CALL step (the fleet form) consumes the same
//! machinery; this shell is the committed example entry point
//! verbatim (project + dedup + cut + the batch/fleet/trace receipt
//! modes — the same binary the committed receipts name, now a thin
//! delegate: the semantics are the library's own, unchanged).

fn main() -> std::io::Result<()> {
    impg::genome_inference::sample_projection::main()
}
