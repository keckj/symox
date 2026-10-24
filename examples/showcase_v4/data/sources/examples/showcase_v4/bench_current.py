"""Reproducible inputs for the refreshed showcase (no rendering dependency).

Workers run in fresh processes. Setup and a smaller warm-up precede timing;
serialization and correctness checks follow it. See refresh.py for the driver.
"""
import bisect
import importlib
import json
import random
import signal
import sys
import time

INTEGRANDS = [
    'x**5*exp(x)', 'sin(x)**6', 'x**3*log(x)**2', 'x**2*cos(x)*sin(x)', 'exp(2*x)*sin(3*x)',
    '(x**3+2*x)/(x**2+1)**2', 'x*atan(x)', 'x/(x**3-1)', 'sin(x)**3*cos(x)**4',
    'x**4*exp(-x)*cos(x)', '1/(x**6-1)', 'log(x**2+1)', 'sinh(x)**4', 'tan(x)**3',
    'x**7*exp(-x**2)', 'cos(x)**8', 'x*sin(x)**2', 'exp(x)*cos(x)**2', '1/(x**3+1)',
    'x**2*atan(x)', 'log(x)**3', 'x**3*sinh(x)', 'asin(x)', '1/(x**2+2*x+5)',
    '(2*x+3)/(x**2+x+1)', 'x**2*exp(3*x)*sin(x)', 'sqrt(1-x**2)', '1/sqrt(x**2+1)',
    'sin(3*x)*cos(5*x)', 'sin(x)**2*cos(x)**2', 'tan(x)**4', 'x*exp(x)*sin(x)', 'exp(-x)*sinh(2*x)',
    'cosh(x)**3', 'x**2*log(x+1)', 'log(x)/x**2', 'x**4/(x**2+1)', '1/(x**4-1)',
    '(x+1)/(x*(x**2+1))', '1/(x*(x+1)**3)', 'x**5/(x**3-1)', 'exp(sqrt(x))',
    'x**3*cos(x**2)', 'x*exp(-x)*cos(2*x)', '1/(x**2*(x**2+1))', 'x**9*exp(x)', 'x**6*sin(x)',
    'cos(x)**5*sin(x)**5', 'x**5*log(x)**2', 'x*log(x**2+1)', 'x**2*sinh(x)*cosh(x)',
    'atan(sqrt(x))', '1/((x-1)**2*(x+2))', 'x**3*exp(-2*x)*sin(x)', 'sqrt(x**2+1)',
    '1/(x**5+x)', 'cos(x)/(1+sin(x)**2)', 'exp(3*x)*cos(2*x)**2', 'x**4*log(x)**3',
    'tan(x)/cos(x)**4', 'x*sinh(x)**2', 'x**8/(x**2+1)', '(3*x**2+1)/(x**3+x+1)',
]


def timed(fn):
    def expired(signum, frame):
        raise TimeoutError('timed operation exceeded 120 seconds')
    signal.signal(signal.SIGALRM, expired)
    start = time.perf_counter()
    signal.setitimer(signal.ITIMER_REAL, 120)
    try:
        result = fn()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    return result, time.perf_counter() - start


def worker(lib, case):
    sp = importlib.import_module(lib)
    x, y = sp.symbols('x y')
    if case == 'integrals':
        ns = vars(sp).copy()
        ns['x'] = x
        fs = [eval(src, ns) for src in INTEGRANDS]
        sp.integrate(x, x)
        items = []
        for src, f in zip(INTEGRANDS, fs):
            result, dt = timed(lambda: sp.integrate(f, x))
            items.append(dict(time=dt, input=src, result=str(result)))
        return dict(items=items)
    if case.startswith('det'):
        n = int(case[3:])
        matrix = sp.Matrix(n, n, lambda i, j: sp.Symbol(f'a{i}{j}'))
        sp.Matrix([[x, 1], [1, y]]).det()
        result, dt = timed(lambda: sp.expand(matrix.det()))
        return dict(time=dt, terms=len(sp.Add.make_args(result)), result=str(result))
    if case == 'taylor_exp':
        sp.series(sp.exp(sp.sin(x)), x, 0, 4)
        result, dt = timed(lambda: sp.series(sp.exp(sp.sin(x)), x, 0, 25))
        return dict(time=dt, result=str(result))
    if case == 'solve':
        equations = [x*x + y*y - 4, y - x**4 + 1]
        sp.solve([x+y-1, x-y], [x, y])
        result, dt = timed(lambda: sp.solve(equations, [x, y]))
        roots = [[str(a), str(b)] for a, b in result]
        numeric = [[str(a.evalf(60)), str(b.evalf(60))] for a, b in result]
        real = [[float(sp.re(a.evalf())), float(sp.re(b.evalf()))] for a,b in result
                if abs(complex(a.evalf()).imag) < 1e-12 and abs(complex(b.evalf()).imag) < 1e-12]
        return dict(time=dt, count=len(result), real=real, roots=roots, numeric=numeric)
    if case.startswith('pendulum'):
        sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]))
        from pendulum_race import equations_of_motion, probe_point
        n = int(case[8:])
        equations_of_motion(sp, 1)
        start = time.perf_counter()
        args, accels = equations_of_motion(sp, n)
        f = sp.lambdify(args, accels)
        dt = time.perf_counter() - start
        probes = [probe_point(n), [0.15*(i+1) for i in range(n)] + [0]*n + [10]]
        values = [[float(v) for v in f(*point)] for point in probes]
        return dict(time=dt, probes=values)
    if case == 'throughput':
        # The identical deterministic stream is used by both libraries. Every
        # input is distinct; constructing it is included in the elapsed time.
        rng = random.Random(20261022)
        specs = [(rng.randint(1,5),rng.randint(1,4),rng.randint(1,4)) for _ in range(256)]
        sp.expand(sp.diff((x+1)**3, x, 2))
        stamps, results, sources = [], [], []
        start = time.perf_counter()
        i = 0
        while time.perf_counter()-start < 6:
            a,b,k = specs[i % len(specs)]
            e = (x + a)**k * sp.sin(b*x + sp.Rational(i+1, 1009))
            r = sp.expand(sp.diff(e, x, 2))
            elapsed = time.perf_counter()-start
            if elapsed > 6:
                break
            stamps.append(elapsed)
            if i < 40:
                results.append(str(r))
                sources.append(str(e))
            i += 1
        return dict(count=len(stamps), elapsed=time.perf_counter()-start,
                    samples=[bisect.bisect_right(stamps, j/60) for j in range(361)],
                    results=results, srcs=sources, seed=20261022)
    raise ValueError(case)


if __name__ == '__main__':
    try:
        result = worker(*sys.argv[1:])
    except TimeoutError:
        result = dict(timeout=120)
    print(json.dumps(result))
