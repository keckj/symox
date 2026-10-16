"""Group fuzz findings by minimal case (ignoring which check found it)."""
import json
import sys
from collections import defaultdict

path = sys.argv[1]
classes = sys.argv[2].split(",") if len(sys.argv) > 2 else None
groups = defaultdict(list)
with open(path) as fh:
    for line in fh:
        r = json.loads(line)
        if classes and r["cls"] not in classes:
            continue
        groups[(r["cls"], r["min_src"])].append(r)
order = {"died": 0, "timeout": 1, "crash": 2, "wrong": 3, "slow": 4, "weak": 5, "form": 6}
for (cls, src), rs in sorted(groups.items(), key=lambda kv: (order.get(kv[0][0], 9), len(kv[0][1]))):
    checks = ",".join(sorted({r["check"] for r in rs}))
    print(f"[{cls}] {src}   <{checks}> x{len(rs)}")
    print(f"      {rs[0]['min_msg'][:300]}")
