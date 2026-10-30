"""Rank-four matrix derivatives and controls.

usage: python bench.py <sympy|symox> <case> <reps>
Prints one JSON line {case, first_ms, median_ms, result}. The module is
imported (and symox._matrix_derivatives loaded) before timing, so the first
call compares the algorithms, not lazy imports.
"""
import json
import statistics
import sys
import time

mod, case, reps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = __import__(mod)
if mod == 'symox':
    import symox._matrix_derivatives  # noqa: F401
    import symox.tensor.array.expressions  # noqa: F401
else:
    import sympy.tensor.array.expressions.from_array_to_matrix  # noqa: F401
    import sympy.tensor.array.expressions.arrayexpr_derivatives  # noqa: F401
MatrixSymbol, Trace, symbols = S.MatrixSymbol, S.Trace, S.symbols
A, B, X = (MatrixSymbol(s, 3, 3) for s in 'ABX')
x = MatrixSymbol('x', 3, 1)
n = symbols('n')
M, N = MatrixSymbol('M', n, n), MatrixSymbol('N', n, n)
CASES = [
    ('(A*X).diff(X)', lambda: (A*X).diff(X)),
    ('(X*A).diff(X)', lambda: (X*A).diff(X)),
    ('(A*X*B).diff(X)', lambda: (A*X*B).diff(X)),
    ('X.T.diff(X)', lambda: X.T.diff(X)),
    ('(X*X).diff(X)', lambda: (X*X).diff(X)),
    ('(X*A*X).diff(X)', lambda: (X*A*X).diff(X)),
    ('(Trace(X)*A).diff(X)', lambda: (Trace(X)*A).diff(X)),
    ('(A*X**-1*B).diff(X)', lambda: (A*X**-1*B).diff(X)),
    ('(N.T*M*N).diff(N)', lambda: (N.T*M*N).diff(N)),
    ('(X**3).diff(X)', lambda: (X**3).diff(X)),
    ('(x.T*x).diff(x) (rank 0, control)', lambda: (x.T*x).diff(x)),
    ('Trace(A*X*B).diff(X) (rank 2, control)', lambda: Trace(A*X*B).diff(X)),
]
name, f = CASES[case]
t0 = time.perf_counter(); r = f(); first = time.perf_counter() - t0
times = []
for _ in range(reps):
    t0 = time.perf_counter(); f(); times.append(time.perf_counter() - t0)
print(json.dumps({'case': case, 'name': name, 'first_ms': first*1e3,
                  'median_ms': statistics.median(times)*1e3, 'result': str(r)}))
