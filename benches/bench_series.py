"""Benchmarks of series expansions and limits: symox vs SymPy.

    python benches/bench_series.py [sympy_time_limit_seconds]

SymPy is run once per case with a cold cache (its cache would otherwise
turn repeated calls into a dictionary lookup) and is interrupted after the
time limit (default 120 s). symox timings are the best of several runs.
"""

from __future__ import annotations

import signal
import sys
import time

import sympy
from sympy.core.cache import clear_cache

import symox

SERIES = [
    "series(exp(sin(x))*cos(x), x, 0, 30)",
    "series(1/cos(x), x, 0, 40)",
    "series(tan(x), x, 0, 30)",
    "series(exp(sin(x)), x, 0, 10)",
    "series(log(1+x)/x, x, 0, 20)",
    "series(1/(1-x-x**2), x, 0, 30)",
    "series(sin(x)/x**3, x, 0, 20)",
    "series(atan(x), x, oo, 10)",
    "series(sqrt(1+sin(x)), x, 0, 12)",
    "series(log(x), x, 1, 10)",
]

LIMITS = [
    "limit(sin(x)/x, x, 0)",
    "limit(x*log(x), x, 0)",
    "limit(x**10*exp(-x), x, oo)",
    "limit((1+1/x)**x, x, oo)",
    "limit((x**2+1)/(2*x**2-3), x, oo)",
    "limit((1-cos(x))/x**2, x, 0)",
    "limit(x**x, x, 0)",
    "limit(sqrt(x**2+x)-x, x, oo)",
    "limit((3**x+5**x)**(1/x), x, oo)",
    "limit(exp(x)*(exp(1/x - exp(-x)) - exp(1/x)), x, oo)",
    "limit(log(log(x))/log(x), x, oo)",
    "limit(1/x - 1/sin(x), x, 0)",
]


def timeit(fn, min_time=0.3, max_runs=200):
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


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def time_sympy(fn, limit):
    """One SymPy run (cold cache), or None after ``limit`` seconds."""
    clear_cache()
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(limit)
    try:
        t0 = time.perf_counter()
        fn()
        return time.perf_counter() - t0
    except _Timeout:
        return None
    finally:
        signal.alarm(0)


def namespace(m):
    ns = {k: getattr(m, k) for k in dir(m) if not k.startswith("_")}
    ns["x"] = m.Symbol("x")
    return ns


def fmt(ts, tx, limit):
    if ts is None:
        return f"{'>' + str(limit) + 's':>11s} {tx*1e3:9.3f}ms {'>' + str(int(limit / tx)) + 'x':>9s}"
    return f"{ts*1e3:9.1f}ms {tx*1e3:9.3f}ms {ts/tx:8.0f}x"


def main(limit=120):
    ns_s, ns_x = namespace(sympy), namespace(symox)
    print(f"{'case':58s} {'sympy':>11s} {'symox':>11s} {'speedup':>9s}")
    for group, cases in (("series", SERIES), ("limits", LIMITS)):
        print(f"--- {group}")
        for src in cases:
            code = compile(src, "<bench>", "eval")
            ts = time_sympy(lambda: eval(code, ns_s), limit)
            tx = timeit(lambda: eval(code, ns_x))
            print(f"{src:58s} {fmt(ts, tx, limit)}", flush=True)
    codes = [compile(s, "<bench>", "eval") for s in LIMITS]
    ts = time_sympy(lambda: [eval(c, ns_s) for c in codes], limit)
    tx = timeit(lambda: [eval(c, ns_x) for c in codes])
    print(f"--- batch of {len(LIMITS)} limits")
    print(f"{'all limits above':58s} {fmt(ts, tx, limit)}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 120)
