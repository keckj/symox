"""pick_cases.py <sympy.jsonl> <symox.jsonl> [per_pass]: prints the cases
pinned by tests/test_native_arr2mat.py. A case is kept when symox agrees
with SymPy 1.14 on the native path, the result differs from the input (or
is an exception) and the input string rebuilds the same expression with the
plain (non-canonicalizing) constructors in SymPy."""
import json
import re
import sys
import sympy as S
from sympy.tensor.array.expressions import array_expressions as ae

k = S.Symbol('k')
x, y = S.symbols('x y')
NS = dict(ArrayContraction=ae.ArrayContraction, ArrayDiagonal=ae.ArrayDiagonal,
          ArrayTensorProduct=ae.ArrayTensorProduct, ArrayAdd=ae.ArrayAdd, OneArray=ae.OneArray,
          ArraySymbol=ae.ArraySymbol, k=k, x=x, y=y, I=S.Identity(k), I3=S.Identity(3))
for s in 'MNP':
    NS[s] = S.MatrixSymbol(s, k, k)
NS.update(a=S.MatrixSymbol('a', k, 1), b=S.MatrixSymbol('b', k, 1), c=S.MatrixSymbol('c', 1, k),
          X=S.MatrixSymbol('X', 3, 3), Y=S.MatrixSymbol('Y', 3, 3), u=S.MatrixSymbol('u', 3, 1),
          A=ae.ArraySymbol('A', (k, k)))
ref = [json.loads(l) for l in open(sys.argv[1])]
got = [json.loads(l) for l in open(sys.argv[2])]
per = int(sys.argv[3]) if len(sys.argv) > 3 else 6
out = {}
for r, g in zip(ref, got):
    inp = r['input']
    if inp != g['input'] or 'PermuteDims' in inp or 'EXC' in inp or 'Trace' in inp or '.T' in inp:
        continue
    try:
        if str(eval(inp.replace('Identity(3)', 'I3'), dict(NS))) != inp:
            continue
    except Exception:
        continue
    for p, v in r.items():
        if p in ('case', 'input', 'rtd', 'remove_id', 'multi_id') or g.get(p) != v:
            continue
        if v == inp or (p == 'editor'):
            continue
        lst = out.setdefault(p, [])
        nexc = sum(1 for c in lst if c[1].startswith('EXC'))
        if v.startswith('EXC') and nexc >= 2:
            continue
        if len(lst) < per and all(c[0] != inp for c in lst):
            lst.append((inp, v))
for p in ('rtd', 'remove_id', 'multi_id', 'editor'):
    for r, g in zip(ref, got):
        inp = r['input']
        if p not in r or inp != g['input'] or 'PermuteDims' in inp or 'Trace' in inp or '.T' in inp or 'EXC' in inp:
            continue
        v = r[p]
        if g.get(p) != v or not (v.endswith('[]]') is False or p == 'editor'):
            continue
        try:
            if str(eval(inp.replace('Identity(3)', 'I3'), dict(NS))) != inp:
                continue
        except Exception:
            continue
        lst = out.setdefault(p, [])
        if v.startswith('EXC') and sum(1 for c in lst if c[1].startswith('EXC')) >= 2:
            continue
        if len(lst) < per and all(c[0] != inp for c in lst):
            lst.append((inp, v))
print('CASES = {')
for p, lst in out.items():
    print('    %r: [' % p)
    for inp, v in lst:
        print('        (%r,\n         %r),' % (inp, v))
    print('    ],')
print('}')
