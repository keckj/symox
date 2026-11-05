"""`bench_a2m.py <sympy|symox> <reps>`: first call and median of `reps`
calls (ms) per case, one JSON line per case with the printed result."""
import importlib
import json
import statistics
import sys
import time

mode, reps = sys.argv[1], int(sys.argv[2])
S = importlib.import_module(mode)
ae = importlib.import_module(mode + '.tensor.array.expressions.array_expressions')
a2m = importlib.import_module(mode + '.tensor.array.expressions.from_array_to_matrix')
m2a = importlib.import_module(mode + '.tensor.array.expressions.from_matrix_to_array')
cm, ca = m2a.convert_matrix_to_array, a2m.convert_array_to_matrix
tp, con, dia, pd = ae._array_tensor_product, ae._array_contraction, ae._array_diagonal, ae._permute_dims
k = S.Symbol('k')
M, N, P = [S.MatrixSymbol(s, k, k) for s in 'MNP']
a, b = S.MatrixSymbol('a', k, 1), S.MatrixSymbol('b', k, 1)
I = S.Identity(k)


def editor_roundtrip(e):
    ed = ae._EditArrayContraction(e)
    return ed.to_array_contraction()


CASES = [
    ('ca(cm(M*N))', cm(M * N), ca),
    ('ca(cm((M*N).T))', cm((M * N).T), ca),
    ('ca(cm(M*N*P))', cm(M * N * P), ca),
    ('ca(cm(a.T*M*b))', cm(a.T * M * b), ca),
    ('ca(cm(Trace(M*N)))', cm(S.Trace(M * N)), ca),
    ('ca(cm(hadamard(M,N)))', cm(S.HadamardProduct(M, N)), ca),
    ('ca(cm(M*N + N.T*P))', cm(M * N + N.T * P), ca),
    ('ca(tp(a,b))', tp(a, b), ca),
    ('ca(dia(tp(M,a),(1,2)))', dia(tp(M, a), (1, 2)), ca),
    ('ca(con(tp(M,I,N),(1,2),(3,4)))', con(tp(M, I, N), (1, 2), (3, 4)), ca),
    ('ca(con(tp(M,N,P),(0,2),(1,5)))', con(tp(M, N, P), (0, 2), (1, 5)), ca),
    ('editor(con(tp(M,N,P),(1,2),(3,4)))', con(tp(M, N, P), (1, 2), (3, 4)), editor_roundtrip),
]

for name, arg, fn in CASES:
    t0 = time.perf_counter()
    try:
        res = fn(arg)
    except Exception as e:
        print(json.dumps({'mode': mode, 'case': name, 'first_ms': None, 'median_ms': None, 'result': 'EXC:' + type(e).__name__}))
        continue
    first = (time.perf_counter() - t0) * 1e3
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn(arg)
        ts.append((time.perf_counter() - t0) * 1e3)
    print(json.dumps({'mode': mode, 'case': name, 'first_ms': first, 'median_ms': statistics.median(ts), 'result': str(res)}))
