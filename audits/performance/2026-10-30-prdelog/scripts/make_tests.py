"""Generate tests/test_native_prdelog.py cases from calls.json and the
SymPy 1.14 replay (replay_sympy.json): every corpus call, plus, per kind,
the first 12 random calls with a None answer and the first 12 without.

usage: make_tests.py CALLS.json SYMPY.json OUT.py
"""
import json
import sys

calls = json.load(open(sys.argv[1]))
ref = json.load(open(sys.argv[2]))['results']
n = calls['corpus']
picked = []
quota = {}
for i, (rec, res) in enumerate(zip(calls['records'], ref)):
    if i < n:
        picked.append([rec[0], rec[1], res[1]])
        continue
    none = res[1] is None or res[1][:2] == ['tuple', None]
    key = (rec[0], none)
    if quota.get(key, 0) < 12:
        quota[key] = quota.get(key, 0) + 1
        picked.append([rec[0], rec[1], res[1]])
with open(sys.argv[3], 'w') as fh:
    fh.write('CASES = [\n')
    for c in picked:
        fh.write('    %s,\n' % json.dumps(c))
    fh.write(']\n')
print(len(picked), 'cases')
