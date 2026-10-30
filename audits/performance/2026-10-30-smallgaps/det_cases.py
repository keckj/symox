"""Float determinant cases: prints one line per case with srepr of det.

Usage: det_cases.py <module> <out.json>   (module: symox or sympy)
"""
import json
import random
import sys

mod = __import__(sys.argv[1])
Matrix, Float, Rational, Symbol, srepr = mod.Matrix, mod.Float, mod.Rational, mod.Symbol, mod.srepr
x, y = Symbol('x'), Symbol('y')


def cases():
    rnd = random.Random(20261030)
    out = []
    # Markov example (rows sum to 0.99999999999999989)
    P = [[0.5, 0.3, 0.2, 0.0], [0.1, 0.6, 0.2, 0.1], [0.25, 0.25, 0.25, 0.25], [0.0, 0.0, 0.3, 0.7]]
    out.append(('markov4', [[1 - P[i][j] if i == j else -P[i][j] for j in range(4)] for i in range(4)]))
    for n in (2, 3, 4, 5, 6, 7):
        for k in range(12):
            kind = k % 4
            rows = []
            for i in range(n):
                row = []
                for j in range(n):
                    if kind == 0:
                        v = Float(rnd.uniform(-5, 5))
                    elif kind == 1:
                        v = rnd.choice([Float(rnd.uniform(-2, 2)), rnd.randint(-3, 3), Rational(rnd.randint(-5, 5), rnd.randint(1, 4))])
                    elif kind == 2:
                        v = rnd.choice([Float(rnd.randint(-3, 3)), 0, 0, Float(rnd.uniform(-1, 1), 30)])
                    else:
                        v = rnd.choice([Float(rnd.uniform(-2, 2)), x, 2 * y, 0, 1])
                    row.append(v)
                rows.append(row)
            out.append(('n%d_k%d' % (n, k), rows))
    # singular-ish Float
    out.append(('sing4', [[1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0], [0.1, 0.2, 0.3, 0.5], [1.5, 0.0, 2.0, 1.0]]))
    out.append(('third4', [[1/3., 1/7., 0.1, 1.0], [0.2, 1/9., 1/11., 0.0], [0.3, 0.7, 1/13., 1.0], [0.0, 1.0, 0.0, 0.3]]))
    out.append(('prec50', [[Float('0.1', 50), Float('0.3', 50), 1, 0], [Float('0.7', 50), 2, 0, 1], [0, 1, Float('0.2', 15), 3], [1, 0, 1, Float('1.1', 50)]]))
    out.append(('block', [[1.5, 0, 0, 0], [2.0, 0.5, 0, 0], [0, 0, 0.25, 3.0], [0, 0, 1.0, 0.125]]))
    return out


import os
SKIP = set(filter(None, os.environ.get('DET_SKIP', '').split(',')))
res = {}
for name, rows in cases():
    M = Matrix(rows)
    r = {}
    print(name, file=sys.stderr, flush=True)
    for meth in (None, 'bareiss', 'berkowitz', 'lu', 'bird', 'laplace', 'domain-ge'):
        print(' ', meth, file=sys.stderr, flush=True)
        if '%s:%s' % (name, meth) in SKIP:
            r[str(meth)] = ['SKIPPED', '']
            continue
        try:
            d = M.det() if meth is None else M.det(method=meth)
            r[str(meth)] = [srepr(d), type(d).__name__]
        except Exception as e:
            r[str(meth)] = ['EXC ' + type(e).__name__ + ': ' + str(e)[:80], '']
    res[name] = r
json.dump(res, open(sys.argv[2], 'w'), indent=0, sort_keys=True)
