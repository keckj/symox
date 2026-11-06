"""Compare corpus outputs: python compare.py raw/sympy.jsonl raw/X.jsonl [raw/Y.jsonl ...]"""
import json, sys
def load(p): return {r['case']: r for r in map(json.loads, open(p))}
ref = load(sys.argv[1]); others = [(p, load(p)) for p in sys.argv[2:]]
fmt = lambda r: ('ERR' if r['error'] else ('-' if r['seconds'] is None else f"{r['seconds']*1000:.1f}"))
print('| case | SymPy ms | ' + ' | '.join(f'{p.split("/")[-1].split(".")[0]} ms / same' for p, _ in others) + ' |')
print('|---|---:|' + '---:|' * len(others))
for case, r in ref.items():
    cells = []
    for _, o in others:
        q = o.get(case)
        if q is None: cells.append('missing'); continue
        same = (q['result'] == r['result']) if r['error'] is None else (q['error'] or '').split(':')[0] == r['error'].split(':')[0]
        cells.append(f"{fmt(q)} / {'yes' if same else 'no'}")
    print(f"| {case} | {fmt(r)} | " + ' | '.join(cells) + ' |')
