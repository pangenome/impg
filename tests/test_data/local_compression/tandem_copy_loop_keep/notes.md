# Tandem copy-number loop should be kept

- Fixture ID: `tandem_copy_loop_keep`
- Class: `tandem_copy_number_loop_cyclic`
- Tier: `ci`
- Assertion: `tandem_copy_loop_required`
- Expected shape: copy-number loop over a tandem motif
- Known failure mode: Normalization or crush linearizes the copy-number loop and loses cyclic topology.
