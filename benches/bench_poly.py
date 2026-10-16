"""Benchmarks of the polynomial subsystem: symox vs SymPy.

    python benches/bench_poly.py [--quick]

Each benchmark builds its input outside of the timed region (expansions
are done once) and times the polynomial operation itself.
"""

from __future__ import annotations

import sys
import time

import symox
import sympy


def timeit(fn, min_time=0.3, max_runs=30, max_total=20.0):
    t0 = time.perf_counter()
    fn()  # warm-up (also the first measurement)
    first = time.perf_counter() - t0
    if first > max_total / 2:
        return first
    runs, total, best = 1, first, first
    while total < min_time and runs < max_runs:
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        best = min(best, dt)
        total += dt
        runs += 1
    return best


def bench_factor_x50(m):
    x = m.Symbol("x")
    e = x**50 - 1
    return lambda: m.factor(e)


def bench_factor_wilkinson(m):
    x = m.Symbol("x")
    p = 1
    for i in range(1, 21):
        p = p * (x - i)
    e = m.expand(p)
    return lambda: m.factor(e)


def bench_factor_dense(m):
    x = m.Symbol("x")
    e = m.expand((x**5 + 3 * x**2 - 7 * x + 11) * (2 * x**7 - x**3 + 5) * (x**6 + x + 1) * (3 * x**4 - 2))
    return lambda: m.factor(e)


def bench_factor_swinnerton_dyer(m):
    # x^8 - 40x^6 + 352x^4 - 960x^2 + 576: irreducible, splits into
    # linear/quadratic factors modulo every prime
    x = m.Symbol("x")
    e = x**8 - 40 * x**6 + 352 * x**4 - 960 * x**2 + 576
    return lambda: m.factor(e)


def bench_factor_xyz10(m):
    x, y, z = m.symbols("x y z")
    e = m.expand((x + y + z) ** 10)
    return lambda: m.factor(e)


def bench_factor_multi(m):
    x, y, z = m.symbols("x y z")
    e = m.expand((x**2 * y + 3 * z - 1) * (x * y * z + y**2 - 2) * (x - y + 2 * z))
    return lambda: m.factor(e)


def bench_factor_bivariate(m):
    x, y = m.symbols("x y")
    e = m.expand((x**4 + x * y**3 + 2 * y - 1) * (x**3 * y**2 - 3 * x + y**4 + 5) * (x**2 - y**2 + x * y + 7))
    return lambda: m.factor(e)


def bench_gcd_multi(m):
    x, y, z = m.symbols("x y z")
    g = m.expand((x + y + z + 1) ** 6)
    a = m.expand(g * (x - y + 2 * z + 3) ** 4)
    b = m.expand(g * (x * y - z + 2) ** 3)
    return lambda: m.gcd(a, b)


def bench_gcd_uni(m):
    x = m.Symbol("x")
    g = m.expand((x**10 + 3 * x**7 - x + 5) ** 3)
    a = m.expand(g * (x**20 - 7 * x**3 + 1))
    b = m.expand(g * (x**15 + x**8 + 2 * x - 9))
    return lambda: m.gcd(a, b)


def bench_cancel(m):
    x, y = m.symbols("x y")
    g = m.expand((x**3 + x * y + y**2 + 1) ** 4)
    num = m.expand(g * (x - y + 1) ** 3)
    den = m.expand(g * (x + 2 * y - 3) ** 2)
    e = num / den
    return lambda: m.cancel(e)


def bench_cancel_sum(m):
    x = m.Symbol("x")
    e = sum(1 / (x + i) for i in range(1, 16))
    return lambda: m.cancel(e)


def bench_apart(m):
    x = m.Symbol("x")
    d = 1
    for i in range(1, 9):
        d = d * (x + i) ** (1 + i % 2)
    e = (x**3 + 1) / m.expand(d)
    return lambda: m.apart(e)


def bench_resultant(m):
    x, y = m.symbols("x y")
    f = x**6 + y * x**3 + y**2 * x + 1
    g = x**5 - y**3 * x**2 + 2 * x - y
    return lambda: m.resultant(f, g, x)


def bench_roots(m):
    x = m.Symbol("x")
    e = m.expand((x - 1) * (x + 2) * (x - 3) * (x**2 - 2) * (x**2 + x + 1))
    return lambda: m.roots(e)


def bench_nroots(m):
    x = m.Symbol("x")
    e = x**30 - 3 * x**7 + x**2 - 1
    return lambda: m.nroots(e)


def bench_sqf(m):
    x, y = m.symbols("x y")
    e = m.expand((x + y) ** 4 * (x**2 - y) ** 3 * (x * y + 1) ** 2)
    return lambda: m.sqf_list(e)


BENCHES = [
    ("factor(x**50 - 1)", bench_factor_x50),
    ("factor(Wilkinson 20)", bench_factor_wilkinson),
    ("factor(dense deg 22, 4 factors)", bench_factor_dense),
    ("factor(Swinnerton-Dyer S3)", bench_factor_swinnerton_dyer),
    ("factor(expand((x+y+z)**10))", bench_factor_xyz10),
    ("factor(trivariate, 3 factors)", bench_factor_multi),
    ("factor(bivariate, 3 factors)", bench_factor_bivariate),
    ("gcd(trivariate, deg 10)", bench_gcd_multi),
    ("gcd(univariate, deg 50)", bench_gcd_uni),
    ("cancel(bivariate rational)", bench_cancel),
    ("cancel(sum 1/(x+i), 15 terms)", bench_cancel_sum),
    ("apart(deg 12 denominator)", bench_apart),
    ("resultant(bivariate)", bench_resultant),
    ("roots(deg 8, mixed)", bench_roots),
    ("nroots(deg 30)", bench_nroots),
    ("sqf_list(bivariate)", bench_sqf),
]


def main():
    quick = "--quick" in sys.argv
    print(f"{'benchmark':36s} {'sympy':>12s} {'symox':>12s} {'speedup':>9s}")
    for name, mk in BENCHES:
        t_x = timeit(mk(symox))
        if quick:
            t_s = float("nan")
        else:
            t_s = timeit(mk(sympy), min_time=0.2, max_runs=5)
        print(f"{name:36s} {t_s*1e3:10.3f}ms {t_x*1e3:10.3f}ms {t_s/t_x:8.1f}x", flush=True)


if __name__ == "__main__":
    main()
