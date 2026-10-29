"""lambdify before/after timings: python bench.py <label> <out.json>.

Run with PYTHONPATH pointing at the symox tree to measure. Every case runs
in a fresh process for the first call, then repeated calls (median of runs).
"""
import json, subprocess, sys, time, statistics

CASES = {
    'build_numpy_sin_poly': ("f=None", "f=S.lambdify(x, S.sin(x)+x**2, 'numpy')"),
    'call_numpy_sin_float': ("f=S.lambdify(x, S.sin(x)+x**2, 'numpy')", "f(0.5)"),
    'call_numpy_sin_npfloat': ("f=S.lambdify(x, S.sin(x)+x**2, 'numpy'); v=np.float64(0.5)", "f(v)"),
    'call_numpy_sin_int': ("f=S.lambdify(x, S.sin(x)+x**2, 'numpy')", "f(2)"),
    'call_numpy_sin_array1000': ("f=S.lambdify(x, S.sin(x)+x**2, 'numpy'); v=np.linspace(0,1,1000)", "f(v)"),
    'call_math_sin_float': ("f=S.lambdify(x, S.sin(x)+x**2, 'math')", "f(0.5)"),
    'call_default_gamma_float': ("f=S.lambdify(x, S.gamma(x))", "f(0.5)"),
    'call_default_rel_float': ("f=S.lambdify(x, x > 0)", "f(0.5)"),
    'call_default_piecewise_float': ("f=S.lambdify(x, S.Piecewise((x, x > 0), (-x, True)))", "f(0.5)"),
    'call_numpy_poly2_float': ("f=S.lambdify([x, y], x**2*y + 3*x - y, 'numpy')", "f(0.5, 1.5)"),
    'getsource_numpy_sin': ("import inspect; f=S.lambdify(x, S.sin(x), 'numpy')", "inspect.getsource(f)"),
}

PRE = "import symox as S, numpy as np; x, y = S.symbols('x y')\n"


def child(name):
    setup, stmt = CASES[name]
    ns = {}
    exec(PRE, ns)
    t0 = time.perf_counter()
    exec(setup, ns)
    try:
        exec(stmt, ns)
        ok = True
    except Exception as e:
        ok = type(e).__name__
    first = time.perf_counter() - t0
    code = compile(stmt, '<s>', 'exec')
    reps = []
    if ok is True:
        for _ in range(7):
            n = 2000 if 'build' not in name else 200
            t = time.perf_counter()
            for _ in range(n):
                exec(code, ns)
            reps.append((time.perf_counter() - t) / n)
    print(json.dumps({'first': first, 'repeat': statistics.median(reps) if reps else None, 'ok': ok}))


if __name__ == '__main__':
    if sys.argv[1] == '--child':
        child(sys.argv[2]); sys.exit()
    label, out = sys.argv[1], sys.argv[2]
    res = {}
    for name in CASES:
        runs = [json.loads(subprocess.check_output([sys.executable, __file__, '--child', name]))
                for _ in range(5)]
        res[name] = {'first_median': statistics.median(r['first'] for r in runs),
                     'repeat_median': statistics.median(r['repeat'] for r in runs) if runs[0]['ok'] is True else None,
                     'ok': runs[0]['ok']}
        print(label, name, res[name], flush=True)
    json.dump(res, open(out, 'w'), indent=1)
