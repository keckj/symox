"""usage: bench_canon.py <sympy|symox> <reps>: median ms per call of the
array-expression canonical forms and of convert_matrix_to_array."""
import importlib, json, statistics, sys, time
mod, reps = sys.argv[1], int(sys.argv[2])
m = importlib.import_module(mod)
try:
    ae = importlib.import_module(mod + '.tensor.array.expressions.array_expressions')
except ImportError:  # the base tree has only the package namespace
    ae = importlib.import_module(mod + '.tensor.array.expressions')
cm = importlib.import_module(mod + '.tensor.array.expressions.from_matrix_to_array').convert_matrix_to_array
tp, con, dia, pd = (getattr(ae, n, None) for n in ('_array_tensor_product', '_array_contraction', '_array_diagonal', '_permute_dims'))
k = m.symbols('k')
M, N, P, Q = [m.MatrixSymbol(s, k, k) for s in 'MNPQ']
a = m.MatrixSymbol('a', k, 1)
cases = [
    ('con(con(tp(M,N,P),(1,2)),(0,1))', lambda: con(con(tp(M, N, P), (1, 2)), (0, 1))),
    ('pd(con(tp(M,N,P),(1,2)),[3,2,1,0])', lambda: pd(con(tp(M, N, P), (1, 2)), [3, 2, 1, 0])),
    ('tp(Q,con(tp(M,N),(1,2)))', lambda: tp(Q, con(tp(M, N), (1, 2)))),
    ('dia(pd(tp(M,N),[1,2,3,0]),(0,1))', lambda: dia(pd(tp(M, N), [1, 2, 3, 0]), (0, 1))),
    ('cm(M*N*P)', lambda: cm(M*N*P)),
    ('cm((M*N).T)', lambda: cm((M*N).T)),
    ('cm(Trace(M*N))', lambda: cm(m.Trace(M*N))),
    ('cm(a.T*M*a)', lambda: cm(a.T*M*a)),
    ('cm(hadamard(M,N,P))', lambda: cm(m.hadamard_product(M, N, P))),
]
for i, (name, f) in enumerate(cases):
    try:
        t0 = time.perf_counter(); r = f(); first = (time.perf_counter() - t0) * 1e3
    except Exception as e:
        print(json.dumps({'case': i, 'name': name, 'error': type(e).__name__})); continue
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter(); f(); ts.append((time.perf_counter() - t0) * 1e3)
    print(json.dumps({'case': i, 'name': name, 'first_ms': first, 'median_ms': statistics.median(ts), 'result': str(r)}))
