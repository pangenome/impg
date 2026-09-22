# Dispersed repeat glue should be broken or ignored

- Fixture ID: `dispersed_repeat_glue_break`
- Class: `dispersed_repeat_glue_break_or_ignore`
- Tier: `local`
- Assertion: `dispersed_repeat_no_long_glue`
- Expected shape: two dispersed repeats without long glue
- Known failure mode: seqwish/SYNG transitive closure glues distant repeat copies into long links.
