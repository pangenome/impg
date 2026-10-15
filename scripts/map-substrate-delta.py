#!/usr/bin/env python3
"""THE HONEST SUBSTRATE DELTA STATEMENT (the three-step map gate):
the fresh map step's cut stats named against the committed cut's own
stats. The committed window-bed substrate of the L1-L4 territories
was disk-cleaned (the owner's disk ruling); the partition step
regenerates the window beds verbatim from the committed locality
maps, so the fresh census is NOT assumed identical to the committed
census - the delta is measured and stated, field by field."""
import json
import sys

old, new = (json.load(open(a)) for a in sys.argv[1:3])
keys = ["reads_scanned", "reads_kept", "touching_patterns",
        "routed_patterns", "census_mult", "territory_rows",
        "territory_entries"]
print("the fresh map cut vs the committed cut (the committed window-bed")
print("substrate was disk-cleaned; the partition step regenerates the")
print("beds verbatim from the committed locality maps):")
for k in keys:
    o, n = old.get(k), new.get(k)
    print(f"  {k}: committed {o} -> fresh {n}"
          + (" (DELTA)" if o != n else " (identical)"))
