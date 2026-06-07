# Nested bubbles where top-level boundary over-merges

- Fixture ID: `nested_top_level_wrong`
- Class: `nested_bubbles_top_level_wrong`
- Tier: `ci`
- Assertion: `nested_parent_boundary_wrong`
- Expected shape: two independent events inside an overbroad parent
- Known failure mode: Top-level selection over-merges independent events and creates bad loops, long links, or excess path depth.
