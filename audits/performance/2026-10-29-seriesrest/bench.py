"""Fresh-process and repeated timings of the seriesrest calls.

Usage: PYTHONPATH=<python dir> python bench.py fresh|repeat <index>
Prints the wall time (seconds) of the call, and its result.
"""
import sys
import time

from symox import (Mod, Piecewise, Rational, S, Symbol, ceiling, floor,
                   limit, oo, symbols)

x, y = symbols("x y")
b = Symbol("b", positive=True)
n = Symbol("n", integer=True)

CALLS = [
    ("floor(x+1/2)-floor(x) at oo", lambda: limit(floor(x + S.Half) - floor(x), x, oo)),
    ("ceiling(x+1/2)-ceiling(x) at oo", lambda: limit(ceiling(x + S.Half) - ceiling(x), x, oo)),
    ("floor(x)-x at oo", lambda: limit(floor(x) - x, x, oo)),
    ("Mod(x, 3) at oo", lambda: limit(Mod(x, 3), x, oo)),
    ("floor(x)/x at oo", lambda: limit(floor(x)/x, x, oo)),
    ("x*floor(3/x)/2 at 0+", lambda: limit(x*floor(3/x)/2, x, 0, "+")),
    ("Mod(x, b) at n*b+", lambda: limit(Mod(x, b), x, n*b, "+")),
    ("(-1)**x*2**(-x) at oo", lambda: limit((-1)**x*2**(-x), x, oo)),
    ("(-1/4)**x at oo", lambda: limit((-Rational(1, 4))**x, x, oo)),
    ("Piecewise series at oo", lambda: Piecewise((x, y < 1), (x**5, True)).series(x, oo, 3)),
]

if __name__ == "__main__":
    mode, i = sys.argv[1], int(sys.argv[2])
    name, f = CALLS[i]
    reps = 1 if mode == "fresh" else 7
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        try:
            r = f()
        except Exception as e:  # the base raises for the Piecewise series
            r = "%s: %s" % (type(e).__name__, e)
        ts.append(time.perf_counter() - t0)
    ts.sort()
    print("%s\t%.6f\t%s" % (name, ts[len(ts)//2], r))
