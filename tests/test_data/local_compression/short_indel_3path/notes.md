# Short insertion/deletion bubble

- Fixture ID: `short_indel_3path`
- Class: `short_indel`
- Tier: `ci`
- Assertion: `short_indel_bubble`
- Expected shape: single insertion/deletion bubble
- Known failure mode: Resolver clips inserted bases, creates a dangling tip, or collapses deletion into the wrong path.
