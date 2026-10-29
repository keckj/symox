"""Before/after timings for the tracearr wave (symbolic-order matrix/array
derivatives, held Derivative of explicit matrices/arrays, constructor
piecewise_fold). Usage: PYTHONPATH=<tree>/python python harness.py <label> [module]
(module defaults to symox; pass sympy, run from /tmp, for the reference).
Each case runs in a fresh process: first call, then the median of 5 batches
of repeated calls. Output: one JSON line per case."""
import json, statistics, subprocess, sys, time

SETUP = """
from symox import *
x, y, k = symbols('x y k'); n = Symbol('n', integer=True)
X = MatrixSymbol('X', 2, 2); A = MatrixSymbol('A', 2, 2); B = MatrixSymbol('B', 2, 2)
"""
CASES = {
    'matsym_order_n': 'X.diff((x, n))',
    'immutable_matrix_order_n': 'ImmutableMatrix([x]).diff((x, n))',
    'array_order_n': 'Array([x, x**2]).diff((x, n))',
    'derivative_matrix_doit': 'Derivative(Matrix([x]), x).doit()',
    'derivative_array_order_n': 'Derivative(Array([x, x**2]), (x, n))',
    'sum_piecewise_fold': 'Sum(x*Piecewise((n, x > 0), (1, True)), (n, 0, k))',
    'integral_piecewise_fold': 'Integral(y*Piecewise((x, x > 0), (1, True)), x)',
    'sum_plain': 'Sum(n**2, (n, 0, k))',
    'integral_plain': 'Integral(x**2, (x, 0, 1))',
    'matsym_diff_x': '(x*X).diff(x)',
    'matrix_diff_x': 'Matrix([[x**2, x]]).diff(x)',
    'nested_trace': 'Trace(A*Trace(B*X)*X).diff(X)',
}
CHILD = """
import time, statistics, json, sys
t0 = time.perf_counter()
exec(compile(sys.argv[1], 's', 'exec'))
setup = time.perf_counter() - t0
code = compile(sys.argv[2], 'c', 'eval')
t0 = time.perf_counter()
try:
    r = eval(code); out = str(r)[:120]
except Exception as e:
    out = 'EXC %s: %s' % (type(e).__name__, str(e)[:80])
first = time.perf_counter() - t0
reps = int(sys.argv[3]); batches = []
for _ in range(5):
    t0 = time.perf_counter()
    for _ in range(reps):
        try: eval(code)
        except Exception: pass
    batches.append((time.perf_counter() - t0) / reps)
print(json.dumps({'first_s': first, 'repeat_s': statistics.median(batches), 'result': out}))
"""
label = sys.argv[1]
module = sys.argv[2] if len(sys.argv) > 2 else 'symox'
SETUP = SETUP.replace('from symox import', 'from %s import' % module)
for name, expr in CASES.items():
    runs = []
    for _ in range(3):
        p = subprocess.run([sys.executable, '-c', CHILD, SETUP, expr, '50'], capture_output=True, text=True)
        runs.append(json.loads(p.stdout.strip().splitlines()[-1]))
    print(json.dumps({'label': label, 'case': name, 'expr': expr,
                      'first_s': statistics.median(r['first_s'] for r in runs),
                      'repeat_s': statistics.median(r['repeat_s'] for r in runs),
                      'result': runs[0]['result']}), flush=True)
