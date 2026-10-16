"""Benchmarks: symox vs SymPy (and SymEngine when installed).

    python benches/bench_vs_sympy.py
"""

from __future__ import annotations

import time

import numpy as np

import symox
import sympy

try:
    import symengine
except ImportError:  # pragma: no cover
    symengine = None


def timeit(fn, min_time=0.3, max_runs=50):
    fn()  # warm-up
    runs, total, best = 0, 0.0, float("inf")
    while total < min_time and runs < max_runs:
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        best = min(best, dt)
        total += dt
        runs += 1
    return best


def bench_expand(m):
    x, y, z, w = m.symbols("x y z w")
    e = (x + y + z + w + 1) ** 15
    return lambda: m.expand(e)


def bench_expand2(m):
    x, y, z, w = m.symbols("x y z w")
    e = (x + y + z + w) ** 15
    f = e * (e + w)
    return lambda: m.expand(m.expand(e) * (m.expand(e) + w))


def bench_build_sum(m):
    x = m.Symbol("x")

    def run():
        s = 0
        for i in range(2000):
            s = s + i * x**i
        return s

    return run


def bench_diff(m):
    x, y = m.symbols("x y")
    e = m.sin(x * y) ** 10 * m.exp(x**2) * m.log(x + y) / (x**2 + y**2 + 1)

    def run():
        r = e
        for _ in range(4):
            r = m.diff(r, x)
        return r

    return run


def bench_subs(m):
    x, y, z = m.symbols("x y z")
    p = m.expand((x + y + z + 1) ** 12)
    return lambda: p.subs({x: y + 1, z: 2})


def bench_sympify(m):
    s = " + ".join(f"{i}*x**{i}*sin(y)/(z + {i})" for i in range(300))
    return lambda: m.sympify(s)


def bench_lambdify(m):
    x, y = m.symbols("x y")
    e = m.sin(x) ** 2 + m.cos(y) * m.exp(-x * y) + m.sqrt(x**2 + y**2 + 1)
    xs = np.linspace(0, 1, 1_000_000)
    ys = np.linspace(1, 2, 1_000_000)
    if m is sympy:
        f = sympy.lambdify((x, y), e, "numpy")
    elif m is symox:
        f = symox.lambdify((x, y), e)
    else:
        f = symengine.Lambdify([x, y], [e])
        return lambda: f(np.stack([xs, ys], axis=1))
    return lambda: f(xs, ys)


BENCHES = [
    ("expand((x+y+z+w+1)**15)", bench_expand),
    ("expand(e*(e+w)), e=(x+y+z+w)**15", bench_expand2),
    ("build sum_{i<2000} i*x**i", bench_build_sum),
    ("4th derivative of a product", bench_diff),
    ("subs on 455-term polynomial", bench_subs),
    ("sympify 300-term string", bench_sympify),
    ("lambdify + eval on 1e6 points", bench_lambdify),
]


def main():
    mods = [("sympy", sympy), ("symox", symox)]
    if symengine is not None:
        mods.insert(1, ("symengine", symengine))
    header = f"{'benchmark':38s}" + "".join(f"{n:>12s}" for n, _ in mods) + f"{'speedup':>10s}"
    print(header)
    print("-" * len(header))
    for name, b in BENCHES:
        times = {}
        for mn, m in mods:
            try:
                fn = b(m)
                slow = mn == "sympy"
                times[mn] = timeit(fn, min_time=0.1 if slow else 0.3, max_runs=3 if slow else 50)
            except Exception as e:  # noqa: BLE001
                times[mn] = None
                print(f"  ({mn} failed: {type(e).__name__}: {e})")
        row = f"{name:38s}"
        for mn, _ in mods:
            t = times[mn]
            row += f"{'n/a':>12s}" if t is None else f"{t * 1e3:10.2f}ms"
        if times.get("sympy") and times.get("symox"):
            row += f"{times['sympy'] / times['symox']:9.0f}x"
        print(row)


if __name__ == "__main__":
    main()
