"""limheld timings: python bench.py OUT.json (run with PYTHONPATH at the tree).

For each call: first call in a fresh process (cache cleared) and the
median of repeated calls; the printed value is recorded as well.
"""
import json, subprocess, sys, time

CASES = [
    "limit(y/x/log(x), x, 0)",
    "limit(-y*z*exp(x)/x, x, oo)",
    "limit(xr**n - xr**(n - 5), xr, oo)",
    "limit(xp**n - xp**(n + 2), xp, oo)",
    "limit(Piecewise((x**2/2, x <= 0.5), (x/2 - 0.125, True)), x, 0)",
    "limit(Piecewise(((x - 9)/5, x < -1), ((x - 9)/5, x > 4), (sqrt(Abs(x - 3)), True)), x, -1)",
    "limit(2*sqrt(x)*Piecewise(((4*x - 2)/Abs(sqrt(4 - 4*(2*x - 1)**2)), 4*x - 2 >= 0), ((2 - 4*x)/Abs(sqrt(4 - 4*(2*x - 1)**2)), True)), x, 0)",
    "(1 + Order(x)).limit(x, 0)",
    "((x + Order(x**2))/x).limit(x, 0)",
    "limit(sin(x)/x, x, 0)",
    "limit(exp(x)*(exp(1/x - exp(-x)) - exp(1/x)), x, oo)",
]
SETUP = """
from symox import *
x, y, z = symbols('x y z'); n = Symbol('n', integer=True, positive=True)
xr = Symbol('x', real=True); xp = Symbol('x', positive=True)
"""

def child(case, reps):
    ns = {}
    exec(SETUP, ns)
    t0 = time.perf_counter(); v = eval(case, ns); first = time.perf_counter() - t0
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter(); eval(case, ns); ts.append(time.perf_counter() - t0)
    ts.sort()
    return {"value": str(v), "first": first, "median": ts[len(ts)//2]}

if __name__ == "__main__":
    if sys.argv[1] == "--child":
        print(json.dumps(child(CASES[int(sys.argv[2])], 7)))
        sys.exit()
    out = {}
    for i, c in enumerate(CASES):
        runs = []
        for _ in range(3):
            r = subprocess.run([sys.executable, __file__, "--child", str(i)], capture_output=True, text=True, timeout=300)
            runs.append(json.loads(r.stdout.strip().splitlines()[-1]))
        firsts = sorted(r["first"] for r in runs); meds = sorted(r["median"] for r in runs)
        out[c] = {"value": runs[0]["value"], "first": firsts[1], "median": meds[1]}
        print(c, out[c], flush=True)
    json.dump(out, open(sys.argv[1], "w"), indent=1)
