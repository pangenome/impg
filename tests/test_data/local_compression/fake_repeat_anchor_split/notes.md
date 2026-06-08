# Bubble split by fake internal repeat anchor

- Fixture ID: `fake_repeat_anchor_split`
- Class: `bubble_split_by_fake_repeat_anchor`
- Tier: `ci`
- Assertion: `fake_repeat_anchor_not_boundary`
- Expected shape: bubble with internal fake repeat anchor
- Known failure mode: Candidate discovery treats the repeat as a hard anchor and under-compresses the true event.
