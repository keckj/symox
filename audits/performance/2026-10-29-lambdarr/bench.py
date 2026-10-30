"""lambdify array timings: python bench.py <label> <out.json> [case-filter].

Run with PYTHONPATH pointing at the symox tree to measure (BENCH_MODULE=sympy
measures the reference SymPy instead; run it from /tmp). Each case runs in
PROCS fresh processes: "first" is lambdify + the first call (symox/numpy
imports excluded), "repeat" the median per-call time of 7 batches.
"""
import json, os, subprocess, sys, time, statistics

EXPRS = {
    'sin': ("x", "S.sin(x)", "'numpy'"),
    'poly': ("x", "x**2 + 3*x + 1", "'numpy'"),
    'expcos': ("x", "S.exp(-x**2)*S.cos(x)", "'numpy'"),
    'hypot': ("(x, y)", "S.sqrt(x**2 + y**2)", "'numpy'"),
    'piecewise': ("x", "S.Piecewise((x**2, x < 1), (S.sin(x), True))", "'numpy'"),
    'gamma': ("x", "S.gamma(x)", "None"),
    'poly3': ("(x, y, z)", "x**2*y + 3*x*z - y*z**3 + 2", "'numpy'"),
}
SIZES = [10, 1000, 100000]
DTYPES = ['f64', 'c128']
SCALARS = {
    'gamma_float': ("x", "S.gamma(x)", "None", "0.5"),
    'poly2_float': ("[x, y]", "x**2*y + 3*x - y", "'numpy'", "0.5, 1.5"),
}

CASES = {}
for name, (args, expr, mod) in EXPRS.items():
    nargs = 1 if args == 'x' else args.count(',') + 1
    for dt in DTYPES:
        for n in SIZES:
            data = "np.linspace(0.1, 2, %d)" % n
            if dt == 'c128':
                data = "(%s + 0.25j)" % data
            vals = ", ".join("V%d" % i for i in range(nargs))
            setup_data = "; ".join("V%d = %s*%s" % (i, data, 1 + i / 4) for i in range(nargs))
            CASES['%s_%s_%d' % (name, dt, n)] = (setup_data, "f=S.lambdify(%s, %s, %s)" % (args, expr, mod), "f(%s)" % vals)
for name, (args, expr, mod, vals) in SCALARS.items():
    CASES[name] = ("", "f=S.lambdify(%s, %s, %s)" % (args, expr, mod), "f(%s)" % vals)

PRE = "import %s as S, numpy as np; x, y, z = S.symbols('x y z')\n" % os.environ.get('BENCH_MODULE', 'symox')
PROCS = int(os.environ.get('BENCH_PROCS', '3'))


def child(name):
    data, build, stmt = CASES[name]
    ns = {}
    exec(PRE, ns)
    exec(data, ns)
    t0 = time.perf_counter()
    try:
        exec(build, ns)
        exec(stmt, ns)
        ok = True
    except Exception as e:
        ok = type(e).__name__
    first = time.perf_counter() - t0
    code = compile(stmt, '<s>', 'exec')
    reps = []
    if ok is True:
        t = time.perf_counter(); exec(code, ns); one = time.perf_counter() - t
        n = max(1, min(20000, int(0.05 / max(one, 1e-7))))
        for _ in range(7):
            t = time.perf_counter()
            for _ in range(n):
                exec(code, ns)
            reps.append((time.perf_counter() - t) / n)
    print(json.dumps({'first': first, 'repeat': statistics.median(reps) if reps else None, 'ok': ok}))


if __name__ == '__main__':
    if sys.argv[1] == '--child':
        child(sys.argv[2]); sys.exit()
    label, out = sys.argv[1], sys.argv[2]
    pattern = sys.argv[3] if len(sys.argv) > 3 else ''
    res = {}
    for name in CASES:
        if pattern not in name:
            continue
        runs = [json.loads(subprocess.check_output([sys.executable, __file__, '--child', name]))
                for _ in range(PROCS)]
        ok = runs[0]['ok']
        res[name] = {'first_median': statistics.median(r['first'] for r in runs),
                     'repeat_median': statistics.median(r['repeat'] for r in runs) if ok is True else None,
                     'ok': ok}
        print(label, name, res[name], flush=True)
    json.dump(res, open(out, 'w'), indent=1)
