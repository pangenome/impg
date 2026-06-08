# Quality Pass: Localized SmoothXG-Like Polishing Workgroup

Task: `quality-pass-for`
Date: 2026-06-08

## Final Task Graph

The workgroup is now structured as:

1. `quality-pass-for`
2. `audit-smoothxg-pggb` and `design-localized-smoothxg` in parallel
3. `implement-explicit-local`
4. `implement-graph-badness`
5. `implement-localized-resolver`
6. `integrate-iterative-localized`
7. `validate-localized-polishing`
8. `synthesize-localized-polishing`

The parallel stage is limited to independent audit/design work. Implementation tasks are serialized because they are likely to touch overlapping graph-output, replacement, path-preservation, and reporting code.

## Dependency Changes

- Added `implement-graph-badness` after `implement-explicit-local`.
- Added `implement-localized-resolver` after `implement-graph-badness`.
- Removed redundant direct audit/design dependencies from `implement-graph-badness` and `implement-localized-resolver`; those now flow transitively through `implement-explicit-local`.
- Removed redundant direct `implement-explicit-local` and `implement-graph-badness` dependencies from `integrate-iterative-localized`; integration now waits on `implement-localized-resolver`, which carries the implementation chain transitively.

## Spec Tightening

Downstream task descriptions were tightened to keep the workgroup aimed at localized smoothing over SYNG-collected local sequence sets:

- Audit and design tasks now explicitly map PGGB/SmoothXG mechanics onto `SYNG sequence collection -> region-level seed induction -> localized dirty-region smoothing`.
- Code tasks now require `cargo build`, `cargo test`, targeted tests, `git diff --check`, and a generated-artifact scan.
- Validation and synthesis tasks now treat graph-quality metrics as diagnostic evidence, not hidden acceptance or rejection gates.
- Exact path corruption is called out as the hard invariant across implementation, integration, and validation.
- Every downstream task now includes no-big-artifact constraints: no generated GFA/GFA.gz/GFA.zst/PNG/PDF, no bulky per-fixture logs/artifacts, and no staged blob over 1 MiB.
- C4 remains a required diagnostic validation case, but tasks now explicitly avoid broad blind C4-only algorithm churn and require non-C4 SVR coverage when local data is available.

No source/code files were changed for this quality pass.
