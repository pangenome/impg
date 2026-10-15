#!/usr/bin/env python3
"""THE HONEST CALLS-DELTA STATEMENT (the three-step call gate): the
three-step path's calls named against the committed receipt of
record, record by record - the fresh-cut substrate's measured effect
on the product, never assumed identical, never silently accepted."""
import json
import sys

old = [json.loads(l) for l in open(sys.argv[1])]
new = [json.loads(l) for l in open(sys.argv[2])]
print(f"the three-step calls vs the committed receipt: {len(new)} records")
same = sum(1 for o, n in zip(old, new) if o == n)
print(f"  byte-equal records: {same}/{len(new)}")
for o, n in zip(old, new):
    if o != n:
        print(f"  locus {n.get('locus')}: DELTA (the fresh-cut substrate)")
if len(old) != len(new):
    print(f"  record-count delta: committed {len(old)} vs fresh {len(new)}")
