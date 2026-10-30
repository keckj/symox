"""array_derive, convert_indexed_to_array and convert_array_to_indexed timings.

usage: python bench.py <sympy|symox> <reps>
Prints one JSON line per case {case, first_ms, median_ms, result}. The
modules are imported before timing; first_ms is the first call, median_ms
the median of `reps` further calls. A failing call prints its exception.
"""
import importlib
import json
import statistics
import sys
import time

mod, reps = sys.argv[1], int(sys.argv[2])
m = importlib.import_module(mod)
ae = importlib.import_module(mod + '.tensor.array.expressions')
try:
    array_derive = importlib.import_module(mod + '.tensor.array.expressions.arrayexpr_derivatives').array_derive
except ImportError:  # the base tree only exports array_derive from the package
    array_derive = ae.array_derive
c = importlib.import_module(mod + '.tensor.array.expressions.from_indexed_to_array').convert_indexed_to_array
a2i = importlib.import_module(mod + '.tensor.array.expressions.from_array_to_indexed').convert_array_to_indexed
tp, con, dia, pd = ae.ArrayTensorProduct, ae.ArrayContraction, ae.ArrayDiagonal, ae.PermuteDims
i, j, k, l, d = m.symbols('i j k l d')
X, Y = m.MatrixSymbol('X', k, k), m.MatrixSymbol('Y', k, k)
M, N, P = (m.MatrixSymbol(s, d, d) for s in 'MNP')
A3, B = ae.ArraySymbol('A', (k, k, k)), ae.ArraySymbol('B', (k, k))
E, F = ae.ArraySymbol('E', (3, 2)), ae.ArraySymbol('F', (2, 3))
Sum = m.Sum

CASES = [
    ('array_derive(X, X)', lambda: array_derive(X, X)),
    ('array_derive(X.T, X)', lambda: array_derive(X.T, X)),
    ('array_derive(X*Y*X, X)', lambda: array_derive(X*Y*X, X)),
    ('array_derive(X**-1, X)', lambda: array_derive(X**-1, X)),
    ('array_derive(tp(A3, B), B)', lambda: array_derive(tp(A3, B), B)),
    ('array_derive(con(tp(X, Y, X), (1, 2), (3, 4)), X)', lambda: array_derive(con(tp(X, Y, X), (1, 2), (3, 4)), X)),
    ('c(Sum(M[i,j]*N[j,k], j))', lambda: c(Sum(M[i, j]*N[j, k], (j, 0, d - 1)))),
    ('c(Sum(M[i,j]*N[j,k]*P[k,l], j, k))',
     lambda: c(Sum(M[i, j]*N[j, k]*P[k, l], (j, 0, d - 1), (k, 0, d - 1)))),
    ('c(M[i,j]*N[k,l] + N[i,j]*M[k,l], [i,j,k,l])', lambda: c(M[i, j]*N[k, l] + N[i, j]*M[k, l], [i, j, k, l])),
    ('c(E[i,j] + F[j,i])', lambda: c(E[i, j] + F[j, i])),
    ('a2i(con(tp(M, N), (1, 2)), [i, j])', lambda: a2i(con(tp(M, N), (1, 2)), [i, j])),
    ('a2i(pd(tp(M, N), [1, 0, 3, 2]), [i, j, k, l])', lambda: a2i(pd(tp(M, N), [1, 0, 3, 2]), [i, j, k, l])),
    ('a2i(dia(tp(M, N), (0, 2)), [i, j, k])', lambda: a2i(dia(tp(M, N), (0, 2)), [i, j, k])),
]

for name, f in CASES:
    try:
        t = time.perf_counter(); r = f(); first = (time.perf_counter() - t)*1e3
        times = []
        for _ in range(reps):
            t = time.perf_counter(); f(); times.append((time.perf_counter() - t)*1e3)
        out = dict(case=name, first_ms=first, median_ms=statistics.median(times), result=str(r))
    except Exception as e:
        out = dict(case=name, first_ms=None, median_ms=None, result='%s: %s' % (type(e).__name__, e))
    print(json.dumps(out), flush=True)
