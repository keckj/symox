"""Median per (mode, case) of raw.tsv, before vs after."""
import collections, statistics, sys
rows = collections.defaultdict(list); tags = {}
for line in open(sys.argv[1] if len(sys.argv) > 1 else 'raw.tsv'):
    label, mode, case, seconds, tag = line.rstrip('\n').split('\t')
    rows[(mode, case, label)].append(float(seconds)); tags[(mode, case, label)] = tag
cases = []
for mode, case, _ in rows:
    if (mode, case) not in cases: cases.append((mode, case))
print('| mode | case | before (ms) | after (ms) | before result | after result |')
print('|---|---|---|---|---|---|')
for mode, case in cases:
    b = statistics.median(rows[(mode, case, 'before')]) * 1e3
    a = statistics.median(rows[(mode, case, 'after')]) * 1e3
    print(f"| {mode} | {case} | {b:.3f} | {a:.3f} | {tags[(mode, case, 'before')]} | {tags[(mode, case, 'after')]} |")
