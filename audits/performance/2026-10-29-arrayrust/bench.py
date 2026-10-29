"""Matrix-derivative and conversion timings.

usage: bench.py <sympy|symox> <native|python> <case index|all> [reps]
  native: symox as built; python: the same binary with the native
  matrix-line kernel hidden from symox._matrix_derivatives (the Python
  algorithm of the parent tree).
Prints one JSON line per case: first call (fresh process when a single
case is given) and the median of repeated calls, in milliseconds.
"""
import json, statistics, sys, time

lib, mode, which = sys.argv[1], sys.argv[2], sys.argv[3]
reps = int(sys.argv[4]) if len(sys.argv) > 4 else 50
if lib == 'symox':
    import symox as L
    from symox.tensor.array.expressions import convert_matrix_to_array, convert_array_to_matrix
    # Imported before timing in both modes so the first call compares the
    # algorithms, not the lazy import of this module.
    import symox._matrix_derivatives as md
    if mode == 'python':
        real = md._core

        class Proxy:
            def __getattr__(self, name):
                if name == '_matrix_lines_derivative':
                    return lambda e, x: None
                return getattr(real, name)
        md._core = Proxy()
else:
    import sympy as L
    from sympy.tensor.array.expressions import convert_matrix_to_array, convert_array_to_matrix
k = L.Symbol('k')
X = L.MatrixSymbol('X', k, k); A = L.MatrixSymbol('A', k, k); B = L.MatrixSymbol('B', k, k)
x = L.MatrixSymbol('x', k, 1); y = L.MatrixSymbol('y', k, 1)
T = L.Trace
CASES = [
    ('(x.T*x).diff(x)', lambda: (x.T*x).diff(x)),
    ('(x.T*A*x).diff(x)', lambda: (x.T*A*x).diff(x)),
    ('Trace(A*X*B).diff(X)', lambda: T(A*X*B).diff(X)),
    ('Trace(X**3).diff(X)', lambda: T(X**3).diff(X)),
    ('Trace(X.T*A*X).diff(X)', lambda: T(X.T*A*X).diff(X)),
    ('(x.T*X**-1*y).diff(X)', lambda: (x.T*X**-1*y).diff(X)),
    ('(Trace(X)*Trace(A*X)).diff(X)', lambda: (T(X)*T(A*X)).diff(X)),
    ('Trace(A*X*B*X).diff(X)', lambda: T(A*X*B*X).diff(X)),
    ('(A*X).diff(X) [rank 4]', lambda: (A*X).diff(X)),
    ('convert_array_to_matrix(convert_matrix_to_array(A*X*B))',
     lambda: convert_array_to_matrix(convert_matrix_to_array(A*X*B))),
]
indices = range(len(CASES)) if which == 'all' else [int(which)]
for i in indices:
    name, f = CASES[i]
    t0 = time.perf_counter(); r = f(); first = (time.perf_counter() - t0)*1e3
    times = []
    for _ in range(reps):
        t0 = time.perf_counter(); f(); times.append((time.perf_counter() - t0)*1e3)
    print(json.dumps({'case': i, 'name': name, 'result': str(r), 'first_ms': first,
                      'median_ms': statistics.median(times)}))
