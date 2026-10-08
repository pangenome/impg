# Exact supported v1 source archive

Five byte-identical source files from the owner-verified a520 baseline (the sixth
identity entry is the version literal, not a sixth file). These are data inputs
for source authentication; no historical search module is compiled. `genome_infer.rs`
is the **entire** historical CLI. The archive was SHA256-verified against the owner
manifest before copying; this local manifest documents those bytes. Do not format
or edit this directory. `transition::supported_identity` recomputes the historical
bindings from these real inputs and their version declaration. The live v2 identity
also binds each archived source independently. The external declaration required
for conversion pins checkpoints/ledgers, not replacement source identities.
