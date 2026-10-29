"""Median over processes of the first-call and repeated-call times in raw.jsonl."""
import collections
import json
import pathlib
import statistics

data = collections.defaultdict(lambda: collections.defaultdict(list))
errors = {}
for line in (pathlib.Path(__file__).parent / 'raw.jsonl').read_text().splitlines():
    r = json.loads(line)
    side = r['label'].split('-')[0]
    data[r['case']][side].append((r['first_us'], r.get('median_us')))
    if r.get('error'):
        errors[(r['case'], side)] = r['error']
print('| case | before first | before repeated | after first | after repeated | errors |')
print('|---|---|---|---|---|---|')
for case, sides in data.items():
    def med(side, i):
        xs = [t[i] for t in sides[side] if t[i] is not None]
        return f'{statistics.median(xs):.0f}' if xs else '-'
    note = ', '.join(f'{s}: {errors[(case, s)]}' for s in ('before', 'after') if (case, s) in errors)
    print(f"| {case} | {med('before', 0)} | {med('before', 1)} | {med('after', 0)} | {med('after', 1)} | {note} |")
