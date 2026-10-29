"""Before/after timing of public solve calls touched by solvesuite.

usage: python bench.py <label>   (PYTHONPATH selects the symox tree)
For each case: first call in a fresh process (subprocess) and the median of
repeated calls (cache cleared between calls when available).
"""
import json, statistics, subprocess, sys, time

CASES = [
    "solve(Eq(x**2 - 1/(x**2 - 4), 4 - 1/(x**2 - 4)), x)",
    "solve(Eq(x + 1/x, 1/x), x)",
    "solve(Eq(x**2, 4), x)",
    "solve(Eq(1/x + x, 2), x)",
    "solve(x**2 - 4, x)",
    "solve(-a*x + 2*x*log(x), x)",
    "solve(1 - oo*x, x)",
    "solve(x - 3, x)",
    "solve([x*y - x - y, x*y - x - y], [x, y])",
    "solve([x**2 + y**2 - 1, x - y], [x, y])",
    "solve((7**x)**x - 2, x)",
    "solve(x**x - 4, x)",
    "solve(exp(x) - 2, x)",
    "solve(sqrt(A**2 + B**2) - 3, A)",
]
SETUP = ("from symox import *\nfrom symox.abc import a, x, y\n"
         "A, B = symbols('a b', imaginary=True)\n")

def one(case, reps):
    g = {}
    exec(SETUP, g)
    try:
        from symox.core.cache import clear_cache
    except Exception:
        clear_cache = lambda: None
    times = []
    for _ in range(reps):
        clear_cache()
        t = time.perf_counter()
        try:
            r = repr(eval(case, g))
        except Exception as e:
            r = "EXC " + type(e).__name__
        times.append(time.perf_counter() - t)
    return r, times

if __name__ == "__main__":
    if sys.argv[1] == "--one":
        r, t = one(sys.argv[2], 1)
        print(json.dumps([r, t[0]]))
        sys.exit()
    out = {}
    for case in CASES:
        firsts = []
        for _ in range(5):
            p = subprocess.run([sys.executable, __file__, "--one", case],
                               capture_output=True, text=True, timeout=300)
            r, t = json.loads(p.stdout.strip().splitlines()[-1])
            firsts.append(t)
        r2, reps = one(case, 15)
        out[case] = {"result": r2, "first_median_s": statistics.median(firsts),
                     "repeat_median_s": statistics.median(reps)}
        print(f"{statistics.median(firsts)*1e3:9.3f} ms {statistics.median(reps)*1e3:9.3f} ms  {case} -> {r2[:70]}")
    json.dump(out, open(sys.argv[1] + ".json", "w"), indent=1)
