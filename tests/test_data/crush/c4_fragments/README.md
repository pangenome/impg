# C4 fragment regression fixtures

These fixtures are small C4-derived traversal fragments used by
`tests/test_crush_integration.rs`.

Source artifact:

```text
/home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z/debug
```

Each `*.fa`/`*.paf` pair was extracted from a real top-flubble SweepGA debug
block by keeping a small subset of FASTA records and only PAF lines whose query
and target names both remain in that subset. The PAF names and lengths are
therefore the actual sequence IDs handed to seqwish in the C4 run.

Fixtures:

- `easy_shared_flank.fa` / `easy_shared_flank.paf`:
  six 226 bp traversals from `graph_build_0000`; two compact allelic classes
  with obvious shared sequence. Protects the simple shared-node induction case.
- `easy_shared_flank.gfa`:
  the same six traversals wrapped in a minimal shared-flank bubble. This is a
  committed lacing fixture for the top-flubble resolver.
- `bounded_multi_bubble.fa` / `bounded_multi_bubble.paf`:
  twelve 193-204 bp traversals from `graph_build_0001`; three compact allelic
  classes. Protects bounded local complexity without requiring the full C4 run.
- `short_floor.fa` / `short_floor.paf`:
  nine 159-300 bp traversals from `graph_build_0002`. Protects the aligner and
  seqwish length-floor path where the configured 311 bp floor must not erase
  homologous short alignments.
- `duplicated_repeat.fa` / `duplicated_repeat.paf`:
  twelve 222 bp traversals from `graph_build_0004` with repeated/duplicated
  sequence classes. Protects filtering and graph induction on duplicated C4
  sequence.

The `unfolded_minrun` case in the test table uses the existing committed
fixtures:

```text
tests/test_data/crush/top_flubble_seqwish_minrun.fa
tests/test_data/crush/top_flubble_seqwish_minrun.paf
```

That fragment is the previous top-flubble failure class: nonempty, name- and
length-consistent C4 PAF reached seqwish, but the effective min-match threshold
was above every exact CIGAR run, so the replacement graph unfolded into one
path-specific segment per traversal.
