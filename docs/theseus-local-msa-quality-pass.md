# Theseus Local MSA Batch Quality Pass

Task: `quality-pass-theseus`

## Batch Order

The Theseus local MSA integration batch is intentionally serialized:

1. `sync-branch-with`
2. `research-theseus-lib`
3. `implement-theseus-backed`
4. `validate-and-document`

This ordering enforces a clean branch prep from latest `origin/main` before
research, design, implementation, testing, or documentation work begins. The
batch must not build on the current `eg/c4-crush-resolution-controls` lineage
unless a human explicitly justifies that dependency in the task log. The
serialized order also prevents implementation and validation agents from editing
the same source, test, dependency, or documentation files in parallel.

## Scope Decisions

- `sync-branch-with` is limited to fetching `origin/main`, preserving the
  current branch and unrelated local/untracked files, and creating or switching
  to a clean Theseus feature branch such as `eg/theseus-local-msa`. It must not
  start Theseus research or implementation.
- `research-theseus-lib` is limited to upstream/repository research and the
  artifact `docs/theseus-local-msa-integration-notes.md`.
- `implement-theseus-backed` owns code, dependency, parser/config, adapter, and
  focused test changes needed for the selectable Theseus-backed path.
- `validate-and-document` runs after implementation and owns final validation,
  documentation, and small test hardening only.

## Acceptance Criteria Added

The refined downstream tasks require explicit criteria for integrating
`https://github.com/albertjimenezbl/theseus-lib`:

- Inspect the upstream Theseus repository at execution time rather than relying
  on assumptions.
- Decide and document the API, binary wrapper, vendoring, or dependency strategy.
- Pin the selected source, revision, or version where possible.
- Review license compatibility and attribution obligations.
- Make the build reproducible from a clean checkout, or document the exact
  external installation path if embedding is not practical.
- Add a selectable `theseus` engine or method without silently changing existing
  `poa`, `seqwish`, `pggb`, or local resolver behavior.
- Convert repository sequence/window/traversal inputs into Theseus inputs and
  convert Theseus outputs into valid local graph/GFA output.
- Preserve source path names, path spellings, sequence lengths, and
  coordinate-like headers.
- Add failing-first tests, parser/config coverage, structural graph validation,
  path/name preservation coverage, and at least one end-to-end Theseus invocation.

## Dependency Check

At this pass, the WG dependencies form the required pipeline:

- `sync-branch-with` is after `quality-pass-theseus` and must prepare the clean
  Theseus feature branch from latest `origin/main`.
- `research-theseus-lib` is after `sync-branch-with`.
- `implement-theseus-backed` is after `research-theseus-lib`.
- `validate-and-document` is after `implement-theseus-backed`.

No two Theseus integration tasks with overlapping likely file scope are left as
parallel open work.
